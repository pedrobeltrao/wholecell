"""The localisation combiner: evidence.tsv -> per-gene compartment weights. One plain function per step (STEPS, a dict
name -> function); the preset's `combiner.steps` lists them in order. Every step reads and writes the per-gene state:

  Gene.src[source]   {compartment: score} of one source (after within-source combination)
  Gene.cls[source]   its evidence class
  Gene.w             the combined weights (set by `combine`, revised by the later steps)
  Gene.label         sources used, joined with '+' (inventory loc_source)

Steps (unification plan §3a):
  drop_unreliable     rows with reliability 0 are dropped (HPA 'Uncertain')
  restrict_membrane   transmembrane genes keep only membrane keys; a source with none left keeps all or is dropped
                      (per class: if_none keep | drop)
  split               multi-localisation within a source: score (score-proportional; equal when scores tie) | rank
                      (1, 1/2, 1/4, ... by score) | equal
  combine             classes -> one distribution: blend (fixed class weights, renormalised over the classes present)
                      | precedence (the first class present in `classes` wins); no evidence -> fallback key
  new_from_imaging    a key supported only by imaging classes is kept only if its source score >= min_score, and never
                      when it would make the protein nucleus+cytosol only
  round               shares rounded to `digits`, shares < min_share dropped (the legacy tables round to 3)
  nc_model            nucleus / cytoplasm split of the soluble pool: fraction (fn from a fractionation class, background
                      b: true = max(obs - b, 0) / (1 - b)) | ratio_exp (r = exp(k (N - C)) from imaging scores,
                      fn = r Vn / (r Vn + Vc)); no evidence -> unchanged
  membrane_mechanism  a membrane / shell key that only imaging supports needs a mechanism feature or two imaging sources

Common-preset steps (preset `common`, unification plan §3a / step B5; criteria docs/unification_criteria/B5_localisation.md):
  annotation_evidence the key list chosen by the cell's key-selection steps (Gene.order; yeast ordered-annotation steps)
                      becomes ONE annotation source with every key scored 1 (class of the source it came from); the
                      consumed raw annotation sources are removed. No key-selection steps (Jurkat): no-op
  resolve_imaging     imaging pseudo-keys resolved per protein (yeast 'periphery' -> plasma membrane for membrane / PM-
                      annotated proteins, cytosol otherwise); imaging ignored for proteins with skip_if_keys (cell wall)
  combine (classes)   precedence classes: an exclusive class (same-cell imaging) alone; else the first annotation class
                      present and the imaging class blended 0.5 : 0.5 at the imaging's resolution (groups: cell data;
                      within a group the annotation's shares decide)
  new_from_imaging    (groups option) a group named only by imaging kept if its best score >= min_score, never a pool
                      key, flagged groups only with the flag
  nc_model (ratio)    r from the cell's imaging calibration (category table | exp(k (N - C))), else 1 when both pool
                      sides are present; fn = r Vn / (r Vn + Vc) with the build's geometry measures; pool sides named by
                      imaging (score >= add_min) or with share >= new_min may be added
Revision 2026-11-18 (docs/unification_criteria/B5_localisation_revision_2026-11-18.md; the frozen B5 gate unchanged):
  imaging_background  HPA nuclear background b: observed nuclear share = b + (1 - b) true of the HPA split distribution
                      -> true = clip((obs - b) / (1 - b), 0, 1) (after split, before combine; HPA only, b from cell data)
  new_from_imaging    background_sources = all: keys of a background-corrected imaging source (HPA) count as supported
                      (blended as in v7); CYCLoPs keeps the min_score / nc-key guards
  nc_model (ratio)    chromatin-bound genes (cell data regexes): whole pool in the nucleus; prior: the fractionation
                      (Itzhak 2016) ratio blended with the imaging ratio in log space (weight w; with no imaging the
                      equal-concentration default takes its place); bound: nuclear share in [1 - m, m] unless >= 2
                      independent sources are exclusive on that side (nc_target, shared with the common v6 R1)
Curated weight overrides (preset common, cell data; 2026-11-19, yeast B5 audit):
  weight_overrides    override rows (field `comp_weights`, ops set / drop / set_keep_other) replace the combined weights of
                      a gene after every model step (before round): curated localisation the evidence model gets wrong
                      (e.g. yeast Sir2/3/4: CYCLoPs scores their nuclear silencing foci as endosome / cytoplasm)

Ordered-annotation steps (yeast legacy_x5; Gene.order = the ranked key list, Gene.unranked, Gene.flags):
  ordered_annotation  the first primary source present gives the ordered list (keys resolved per gene, generic key refined
                      from another source), else the fallback source (unranked; membrane proteins: its membrane terms);
                      membrane proteins with only soluble keys take the membrane source's keys
  annotation_overrides  override rows (field `annotation`, op set) replace the ordered list (yeast LOC_FIX)
  mito_fractionation  organelle fractionation: a measured sub-compartment replaces every organelle key; a proteome
                      table removes the organelle from proteins outside its high-confidence set that are depleted from
                      the pure fraction (and vetoes the measured source), or gives the generic key its sub-class; sets
                      the flag (mito_hc) of proteins in the high-confidence proteome
  rank_weights        rank split by list order (base^i; all equal when unranked), secondary-location down-weighting
                      rules, the pooled keys (cytosol + nucleoplasm) re-split by volume (equal concentration, r = 1)
  imaging_groups      imaging LOC-scores set the share of each compartment group (summed scores), the annotation's
                      rank weights decide within a group; new groups only from new.min_score (never `never`, flagged
                      groups only with the flag, e.g. mitochondrion only when mito_hc); the nucleus / cytosol split of
                      the nc group from r = exp(k (N - C)) of the raw scores; not used for skip_if_keys, no shared
                      group, or a membrane protein without a non-nc group
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd

from .evidence import CLASSES


class Gene:
    __slots__ = ("gene", "src", "cls", "rows", "w", "label", "tm", "features", "order", "unranked", "flags")

    def __init__(self, gene, tm=False, features=None):
        self.gene, self.tm, self.features = gene, bool(tm), features or {}
        self.src, self.cls, self.rows = {}, {}, {}
        self.w, self.label = None, []
        self.order, self.unranked, self.flags = [], False, {}   # ordered annotation (ordered-annotation steps)


def build_state(E: pd.DataFrame, F: pd.DataFrame, genes, within: str = "max") -> dict:
    """evidence rows (file order) -> {gene: Gene}; within-source scores combined by `max` (legacy) or `sum`"""
    S = {g: Gene(g, F.at[g, "transmembrane"] if g in F.index else False,
                 {c: F.at[g, c] for c in F.columns} if g in F.index else {}) for g in genes}
    for r in E.itertuples(index=False):
        st = S.get(r.gene)
        if st is None:
            continue
        q = st.src.setdefault(r.source, {})
        st.cls[r.source] = r.evidence_class
        st.rows.setdefault(r.source, []).append(r)
        if within == "max":
            q[r.compartment] = max(q.get(r.compartment, 0), r.score)
        else:
            q[r.compartment] = q.get(r.compartment, 0) + r.score
    return S


# ------------------------------------------------------------------------------------------------------------- steps
def drop_unreliable(S, o, ctx=None):
    """sources whose rows carry reliability 0 are dropped for that gene (HPA 'Uncertain': reliability is per gene)"""
    for st in S.values():
        for s in list(st.src):
            if any(r.reliability == 0 for r in st.rows[s]):
                del st.src[s]


def restrict_membrane(S, o, ctx=None):
    keys = set(o["keys"])
    pol = o["if_none"]
    for st in S.values():
        if not st.tm:
            continue
        for s in list(st.src):
            q = st.src[s]
            r = {k: v for k, v in q.items() if k in keys}
            if r:
                st.src[s] = r
            elif pol.get(st.cls[s], "keep") == "drop":
                del st.src[s]


def _split(q: dict, rule: str) -> dict:
    if rule == "score":
        t = sum(q.values())
        return {k: v / t for k, v in q.items()}
    if rule == "equal":
        return {k: 1.0 / len(q) for k in q}
    if rule == "rank":                       # yeast X5: 1, 1/2, 1/4 by score (ties keep evidence order), normalised
        order = sorted(q, key=lambda k: -q[k])
        raw = {k: 0.5 ** i for i, k in enumerate(order)}
        t = sum(raw.values())
        return {k: raw[k] / t for k in q}
    raise ValueError(f"split rule '{rule}' (allowed: score, rank, equal)")


def split(S, o, ctx=None):
    """fractionation sources are not distributions (nc_model input only) and are left as scores"""
    for st in S.values():
        for s, q in st.src.items():
            if q and st.cls[s] != "fractionation":
                st.src[s] = _split(q, o)


def _class_dist(st, ss, sw):
    """one distribution of several sources of a class (average with source weights, default 1)"""
    if len(ss) == 1:
        return dict(st.src[ss[0]])
    tw = sum(sw.get(s, 1.0) for s in ss)
    keys = list(dict.fromkeys(k for s in ss for k in st.src[s]))
    return {k: sum(sw.get(s, 1.0) * st.src[s].get(k, 0.0) for s in ss) / tw for k in keys}


def combine_classes(S, o, ctx=None):
    """precedence classes (preset common; plan §3a step 1): o = {exclusive: [classes used alone when present],
    annotation: [classes, first present wins], imaging: [classes, first present wins], weights: {annotation, imaging},
    groups: {key: group} (imaging resolution; keys not listed are their own group), fallback: {compartment, label}}.
    Annotation and imaging present: group share = wa x annotation group share + wm x imaging group share (weights
    renormalised); within a group the annotation's shares, else the imaging's."""
    sw = o.get("source_weights") or {}
    fb = o["fallback"]
    grp = o.get("groups") or {}
    G = lambda k: grp.get(k, k)
    ex, ann, img = list(o.get("exclusive") or []), list(o["annotation"]), list(o["imaging"])
    wts = o["weights"]
    for st in S.values():
        by_cls = {}
        for s, q in st.src.items():
            if q:
                by_cls.setdefault(st.cls[s], []).append(s)
        pick = lambda classes: next((c for c in classes if c in by_cls), None)
        ce = pick(ex)
        if ce is not None:
            st.w, st.label = _class_dist(st, by_cls[ce], sw), list(by_cls[ce])
            continue
        ca, ci = pick(ann), pick(img)
        if ca is None and ci is None:
            st.w, st.label = {fb["compartment"]: 1.0}, [fb["label"]]
            continue
        if ca is None or ci is None:
            c = ca or ci
            st.w, st.label = _class_dist(st, by_cls[c], sw), list(by_cls[c])
            continue
        A, M = _class_dist(st, by_cls[ca], sw), _class_dist(st, by_cls[ci], sw)
        wa, wm = wts["annotation"], wts["imaging"]
        z = wa + wm
        wa, wm = wa / z, wm / z
        gA, gM = {}, {}
        for k, v in A.items():
            gA[G(k)] = gA.get(G(k), 0.0) + v
        for k, v in M.items():
            gM[G(k)] = gM.get(G(k), 0.0) + v
        w = {}
        for g in dict.fromkeys(list(gA) + list(gM)):
            share = wa * gA.get(g, 0.0) + wm * gM.get(g, 0.0)
            src = A if gA.get(g, 0.0) > 0 else M
            inner = {k: v for k, v in src.items() if G(k) == g}
            zi = sum(inner.values())
            for k, v in inner.items():
                w[k] = w.get(k, 0.0) + share * v / zi
        st.w = w
        st.label = list(by_cls[ca]) + list(by_cls[ci])


def combine(S, o, ctx=None):
    """classes in o['classes'] order; within a class several sources are averaged with weight o['source_weights'] (1);
    mode `classes`: combine_classes (preset common)"""
    if o["mode"] == "classes":
        return combine_classes(S, o, ctx)
    classes, mode = o["classes"], o["mode"]
    sw = o.get("source_weights") or {}
    fb = o["fallback"]
    for st in S.values():
        by_cls = {}
        for s, q in st.src.items():
            if q and st.cls[s] in classes:
                by_cls.setdefault(st.cls[s], []).append(s)
        present = [c for c in classes if c in by_cls]
        if not present:
            st.w, st.label = {fb["compartment"]: 1.0}, [fb["label"]]
            continue
        cq = {}
        for c in present:
            ss = by_cls[c]
            if len(ss) == 1:
                cq[c] = st.src[ss[0]]
            else:
                tw = sum(sw.get(s, 1.0) for s in ss)
                keys = list(dict.fromkeys(k for s in ss for k in st.src[s]))
                cq[c] = {k: sum(sw.get(s, 1.0) * st.src[s].get(k, 0.0) for s in ss) / tw for k in keys}
        if mode == "precedence":
            present = present[:1]
        st.label = [s for c in present for s in by_cls[c]]
        if len(present) == 1:
            st.w = dict(cq[present[0]])
            continue
        wt = {c: o["weights"][c] for c in present}
        tw = sum(wt.values())
        wt = {c: v / tw for c, v in wt.items()}
        if o.get("key_order", "first_seen") == "set_union":       # legacy: set(hpa) | set(uniprot) (PYTHONHASHSEED=0)
            keys = set()
            for c in present:
                keys = keys | set(cq[c])
        else:
            keys = list(dict.fromkeys(k for c in present for k in cq[c]))
        st.w = {k: sum(wt[c] * cq[c].get(k, 0) for c in present) for k in keys}


def new_from_imaging(S, o, ctx=None):
    """keys that only imaging classes support: kept if the best imaging score >= min_score and the key is not one of
    nc_keys (nucleus / cytosol come from annotation or the nc_model, never from imaging alone); then renormalised.
    Options of preset common: groups ({key: group}: support and scores are counted per group, the imaging resolution),
    needs_flag ({group: flag}: a new group only for genes with that flag, e.g. yeast mitochondrion / mito_hc),
    no_support_classes (evidence classes that never count as support: common = fractionation, a tie-breaker only),
    background_sources (revision 2026-11-18; 'all': the keys of an imaging source corrected by the step
    imaging_background count as supported, i.e. that source is blended as in v7 (main and additional terms, nucleus and
    cytosol included); the score threshold and the nc-key exclusion were the X5 guards against CYCLoPs LOC-score noise
    and the HPA nuclear floor, which imaging_background now removes explicitly; other imaging sources, e.g. yeast
    CYCLoPs, keep both guards; null: off)"""
    img = set(o["imaging_classes"])
    nc_keys = set(o.get("nc_keys", ["nucleus", "cytosol"]))
    grp = o.get("groups") or {}
    G = lambda k: grp.get(k, k)
    flags = o.get("needs_flag") or {}
    nosup = set(o.get("no_support_classes") or [])   # classes that do not count as support (common: fractionation)
    allbg = o.get("background_sources") == "all"
    for st in S.values():
        if not st.w:
            continue
        other = {G(k) for s, q in st.src.items() if st.cls[s] not in img and st.cls[s] not in nosup for k in q}
        if allbg:                           # keys of background-corrected imaging sources count as supported
            other |= {G(k) for s, q in st.src.items() if st.cls[s] in img and s in st.flags.get("bg_sources", ()) for k in q}
        if not other:
            continue                        # imaging is the only evidence: nothing is "new"
        drop = []
        for k in st.w:
            if G(k) in other:
                continue
            best = max((r.score for s in st.src if st.cls[s] in img for r in st.rows[s] if G(r.compartment) == G(k)),
                       default=0.0)
            fl = flags.get(G(k))
            if best < o["min_score"] or k in nc_keys or (fl is not None and not st.flags.get(fl, False)):
                drop.append(k)
        keep = {k: v for k, v in st.w.items() if k not in drop}
        if drop and keep:
            t = sum(keep.values())
            st.w = {k: v / t for k, v in keep.items()}


def round_shares(S, o, ctx=None):
    d, m = o["digits"], o["min_share"]
    for st in S.values():
        st.w = {k: round(v, d) for k, v in st.w.items() if v >= m}


def _resplit(w: dict, fn: float, o: dict, new_min: float):
    """soluble pool (o.pool keys) re-split: nucleus share fn, cytoplasm 1 - fn kept in proportion over the cytoplasm keys
    present (default key when none); a protein without the nucleus key gets one only if fn >= new_min. None = unchanged."""
    nuc = o["nucleus"]
    w = dict(w)
    a = w.get(nuc, 0.0)
    pool = sum(w.get(k, 0.0) for k in o["pool"])
    if pool == 0 or (a == 0 and fn < new_min):
        return None
    cyt = {k: w[k] for k in o["cytoplasm"] if w.get(k, 0) > 0} or {o["default_cytoplasm"]: 1.0}
    cs = sum(cyt.values())
    for k in o["pool"]:
        w.pop(k, None)
    for k, v in cyt.items():
        w[k] = pool * (1 - fn) * v / cs
    w[nuc] = pool * fn
    return w


def nc_model(S, o, ctx=None):
    meth = o["method"]
    if meth == "ratio":
        return nc_ratio(S, o, ctx)
    n = 0
    for st in S.values():
        if (o.get("skip_transmembrane", True) and st.tm) or not st.w:
            continue
        if meth == "fraction":
            ss = [s for s in st.rows if st.cls[s] == o["class"]]
            if not ss:
                continue
            r = st.rows[ss[0]][0]
            b = o["background"]
            fn = max(r.score - b, 0.0) / (1 - b)
            label = ss[0]
        elif meth == "ratio_exp":
            sc = {}
            src = None
            for s in st.rows:
                if st.cls[s] in o["classes"]:
                    for r in st.rows[s]:
                        if r.compartment in (o["nucleus"], o["cytosol_evidence"]):
                            sc[r.compartment] = r.score
                            src = src or s
            if not sc:
                continue
            ratio = math.exp(o["k"] * (sc.get(o["nucleus"], 0.0) - sc.get(o["cytosol_evidence"], 0.0)))
            vn, vc = o["volumes"]["nucleus"], o["volumes"]["cytoplasm"]
            fn = ratio * vn / (ratio * vn + vc)
            label = src
        else:
            raise ValueError(f"nc_model.method '{meth}' (allowed: fraction, ratio_exp)")
        nw = _resplit(st.w, fn, o, o.get("new_min", 0.0))
        if nw is None:
            continue
        st.w = nw
        st.label = st.label + [label]
        n += 1
    if ctx is not None:
        ctx.say(f"nc_model {meth}: nucleus / cytoplasm re-split for {n} soluble proteins")


def membrane_mechanism(S, o, ctx=None):
    keys, img = set(o["keys"]), set(o["imaging_classes"])
    feat = o.get("mechanism_feature", "mechanism")
    for st in S.values():
        if not st.w or st.tm or st.features.get(feat):
            continue
        drop = []
        for k in st.w:
            if k not in keys:
                continue
            sup = [s for s, q in st.src.items() if k in q]
            if sup and all(st.cls[s] in img for s in sup) and len(sup) < 2:
                drop.append(k)
        keep = {k: v for k, v in st.w.items() if k not in drop}
        if drop and keep:
            t = sum(keep.values())
            st.w = {k: v / t for k, v in keep.items()}


# ------------------------------------------------------------------- ordered-annotation steps (yeast legacy_x5)
def _keys(st, source):
    return list(st.src.get(source) or {})


def _resolve_keys(st, keys, res):
    """per-gene key resolution: {key: {to, genes, feature}}: `key` -> `to` for genes whose feature (gene_name) is listed,
    dropped otherwise; deduplicated
    # origin: yeast pipeline/inventory.py:_spindle_pole (SPB_CORE) @ e3b7448"""
    out = []
    for c in keys:
        r = res.get(c)
        if r is not None:
            c = r["to"] if st.features.get(r.get("feature", "gene_name")) in set(r["genes"]) else None
        if c:
            out.append(c)
    return list(dict.fromkeys(out))


def _refine_keys(keys, extra, generic, subs):
    if generic not in keys:
        return keys
    sub = [c for c in dict.fromkeys(extra) if c in subs]
    if not sub:
        return keys
    out = []
    for c in keys:
        out += sub if c == generic else [c]
    return list(dict.fromkeys(out))


def ordered_annotation(S, o, ctx=None):
    """Gene.order from the ordered annotation sources (rows in source order = rank):
      primary            sources in precedence order; the first with keys gives the list (resolve applied)
      resolve            {key: {to, genes, feature}} per-gene key resolution (yeast spindle pole -> SPB core genes only)
      refine             {generic, sub_keys, source}: when the list has the generic key and `source` names exactly one
                         sub-key, the generic key becomes that sub-key
      fallback           {source, unranked, membrane_source, unless_keys}: no primary keys -> the fallback source (resolve
                         applied; unranked = its order carries no rank); a transmembrane protein then takes the keys of
                         membrane_source when it has any (not when the list has one of unless_keys)
      tm_soluble         {keys, source}: a transmembrane protein whose keys are all in `keys` takes the keys of `source`
    # origin: yeast pipeline/inventory.py:build (locs / src / unranked loop), GO_UNRANKED @ e3b7448"""
    res = o.get("resolve") or {}
    ref = o.get("refine") or {}
    fb = o.get("fallback") or {}
    tms = o.get("tm_soluble") or {}
    for st in S.values():
        a, lab = [], []
        for s in o["primary"]:
            a = _resolve_keys(st, _keys(st, s), res)
            if a:
                lab = [s]
                break
        if a and ref and ref["generic"] in a:
            gc = _keys(st, ref["source"])
            if sum(c in set(ref["sub_keys"]) for c in gc) == 1:
                a = _refine_keys(a, gc, ref["generic"], set(ref["sub_keys"]))
        unranked = False
        if not a and fb:
            a = _resolve_keys(st, _keys(st, fb["source"]), res)
            if a:
                lab = [fb["source"]]
                unranked = bool(fb.get("unranked"))
                if st.tm and fb.get("membrane_source") and not set(fb.get("unless_keys") or []) & set(a):
                    a = _keys(st, fb["membrane_source"]) or a
        if tms and st.tm and all(c in set(tms["keys"]) for c in a):
            gm = _keys(st, tms["source"])
            if gm:
                a, lab = gm, lab + [tms["source"]]
        st.order, st.unranked, st.label = a, unranked, lab or ["none"]


def annotation_overrides(S, o, ctx=None):
    """override rows (o['rows']: the cell's override table; rule_set / rule_ids select) with field `annotation`, op set
    (value = the ordered key list) replace Gene.order, through overrides.apply (yeast LOC_FIX)"""
    from . import overrides as OV
    T = OV.check(pd.DataFrame(o["rows"]).astype(str), "data.overrides")
    T = OV.select(T, o["rule_set"], o["rule_ids"])
    T = T[T.field == "annotation"]

    def put(g, r, new):
        S[g].order, S[g].unranked = list(new), False
        S[g].label = [f"override:{r.rule_id}"]

    n = OV.apply(T, lambda kt, k: [k] if k in S else [], lambda g, f: S[g].order, put)
    if ctx is not None:
        ctx.say(f"annotation_overrides: {n} rows applied ({', '.join(o['rule_ids'])})")


def mito_fractionation(S, o, ctx=None):
    """organelle fractionation on Gene.order (keys = the organelle's keys, generic first):
      measured           source whose (single) key replaces every organelle key of the list (Vogtle 2017)
      proteome           source with one row per gene: compartment = sub-class key (generic when none), score 1 = in the
                         high-confidence proteome, reliability 0 = depleted from the pure organelle (Morgenstern 2017)
      veto               a protein outside the high-confidence proteome and depleted is not in the organelle: the
                         measured source is not used, the organelle keys are removed; an empty list then takes
                         not_mito_fallback ({sources: [{source, transmembrane_only}], default}) without organelle keys
      tm_map             sub-class key -> key for transmembrane proteins (soluble matrix / IMS class -> inner membrane)
      flag               Gene.flags[flag] = in the high-confidence proteome and not removed
    # origin: yeast pipeline/inventory.py:apply_vogtle, apply_morgenstern, build (VOGTLE_NEEDS_MITO) @ e3b7448"""
    keys = list(o["keys"])
    org = set(keys)
    generic = o["generic"]
    tm_map = o.get("tm_map") or {}
    nmf = o.get("not_mito_fallback") or {}
    n = {"measured": 0, "removed": 0, "sub": 0}
    for st in S.values():
        prow = (st.rows.get(o["proteome"]) or [None])[0]
        hc = prow is not None and prow.score == 1.0
        depleted = prow is not None and prow.reliability == 0.0
        vetoed = bool(o.get("veto")) and prow is not None and not hc and depleted
        mk = _keys(st, o["measured"])
        not_org = False
        if mk and not vetoed:
            if any(c in org for c in st.order):
                new = list(dict.fromkeys(mk[0] if c in org else c for c in st.order))
                if new != st.order:
                    st.order = new
                    st.label = st.label + [o["measured"]]
                    n["measured"] += 1
        if prow is not None and (not mk or vetoed) and any(c in org for c in st.order):
            if not hc and depleted:
                rest = [c for c in st.order if c not in org]
                if not rest:
                    for f in nmf.get("sources") or []:
                        if f.get("transmembrane_only") and not st.tm:
                            continue
                        rest = [c for c in _keys(st, f["source"]) if c not in org]
                        if rest:
                            break
                    rest = rest or [nmf["default"]]
                st.order = rest
                st.label = st.label + [o["proteome"] + "_not_organelle"]
                not_org = True
                n["removed"] += 1
            elif generic in st.order and prow.compartment != generic:
                tgt = tm_map.get(prow.compartment, prow.compartment) if st.tm else prow.compartment
                st.order = list(dict.fromkeys(tgt if c == generic else c for c in st.order))
                st.label = st.label + [o["proteome"] + "_sub"]
                n["sub"] += 1
        st.flags[o["flag"]] = hc and not not_org
    if ctx is not None:
        ctx.say(f"mito_fractionation: measured sub-compartment {n['measured']}, removed {n['removed']}, sub-class {n['sub']}")


def _rank(order, o, flags=None, equal=False):
    """rank weights of an ordered key list (see rank_weights)
    # origin: yeast pipeline/inventory.py:weights (MITO_SECONDARY, equal-concentration pool) @ e3b7448"""
    if not order:
        return dict(o["empty"])
    flags = flags or {}
    base = o.get("base", 0.5)
    raw = {c: (1.0 if equal else base ** i) for i, c in enumerate(order)}
    for r in o.get("secondary") or []:
        if r.get("flag") and not flags.get(r["flag"]):
            continue
        if order[0] in set(r["primary"]):
            for c in order[1:]:
                if (c in set(r["keys"])) if "keys" in r else (c not in set(r["keys_not"])):
                    raw[c] *= r["factor"]
    pool = o.get("pool")
    if pool:
        vol = pool["volumes"]
        cn = [c for c in order if c in set(pool["keys"])]
        if cn:
            pooled = sum(raw[c] for c in cn)
            tot = sum(vol[c] for c in cn)
            for c in cn:
                raw[c] = pooled * vol[c] / tot
    z = sum(raw.values())
    return {c: v / z for c, v in raw.items()}


def rank_weights(S, o, ctx=None):
    """Gene.w from Gene.order: base^i by list position (1, 1/2, 1/4 ...; every key 1 when Gene.unranked), then the
    `secondary` rules ({primary, keys | keys_not, factor, flag}: when the first key is in `primary`, later keys in `keys`
    (or not in `keys_not`) are multiplied by factor; a rule with `flag` only for genes with that flag), then the `pool`
    keys ({keys, volumes}) pooled and re-split by volume (equal concentration across the pore, r = 1), normalised;
    no keys -> `empty`"""
    for st in S.values():
        st.w = _rank(st.order, o, st.flags, st.unranked)


def imaging_groups(S, o, ctx=None):
    """imaging LOC-scores x annotation, group-wise (options: see the module docstring and the legacy_x5 preset)
    # origin: yeast pipeline/inventory.py:cyclops_weights, nucleocyto_split, build (CYCLoPs loop) @ e3b7448"""
    groups = o.get("groups") or {}
    G = lambda c: groups.get(c, c)
    per = o.get("periphery") or {}
    skip = set(o.get("skip_if_keys") or [])
    new = o["new"]
    nc = o["nc"]
    rk = o["rank"]
    pm_keys = set(per.get("annotated") or [per.get("membrane")])
    n = 0
    for st in S.values():
        rows = st.rows.get(o["source"]) or []
        if not rows:
            continue
        comps = st.order
        if skip & set(comps):
            continue
        pm_ok = st.tm or bool(pm_keys & set(comps))
        cls_c = {}
        for r in rows:
            c = r.compartment
            if per and c == per["key"]:
                c = per["membrane"] if pm_ok else per["other"]
            cls_c[c] = cls_c.get(c, 0.0) + r.score
        gs = {}
        for c, v in cls_c.items():
            gs[G(c)] = gs.get(G(c), 0.0) + v
        ann_g = {}
        for c, f in (_rank(comps, rk) if comps else {}).items():
            ann_g.setdefault(G(c), {})[c] = f
        if comps and not set(ann_g) & set(gs):
            continue

        def new_ok(g, v):
            fl = (new.get("needs_flag") or {}).get(g)
            return v >= new["min_score"] and g not in set(new.get("never") or []) and (fl is None or st.flags.get(fl, False))

        gs = {g: v for g, v in gs.items() if g in ann_g or new_ok(g, v)}
        if not gs or (st.tm and not any(g != nc["group"] for g in gs)):
            continue
        z = sum(gs.values())
        w = {}
        for g, v in gs.items():
            inner = ann_g[g] if g in ann_g else _rank([c for c in cls_c if G(c) == g], rk)
            if g == nc["group"]:
                inner = _nc_split(rows, inner, nc)
            zi = sum(inner.values())
            for c, f in inner.items():
                w[c] = w.get(c, 0.0) + v / z * f / zi
        st.w = w
        st.label = st.label + [o["source"]]
        n += 1
    if ctx is not None:
        ctx.say(f"imaging_groups: {o['source']} sets the group shares of {n} proteins")


def _nc_split(rows, inner, nc):
    """nucleus / cytosol re-split of a group's weights: r = exp(k (N - C)) of the raw scores (note = term), shares
    r Vn / (r Vn + Vc); numpy exp as the legacy code
    # origin: yeast pipeline/inventory.py:nucleocyto_split (NUC_K) @ e3b7448"""
    nu, cy, k = nc["nucleus"], nc["cytosol"], nc["k"]
    if not rows or cy not in inner or nu not in inner or k <= 0:
        return inner
    sc = {r.note: r.score for r in rows}
    d = sc.get(nc["nucleus_term"], 0.0) - sc.get(nc["cytosol_term"], 0.0)
    tot = inner[cy] + inner[nu]
    en = float(np.exp(k * d))
    cn = en * nc["volumes"][nu]
    cc = nc["volumes"][cy]
    out = dict(inner)
    out[nu] = tot * cn / (cn + cc)
    out[cy] = tot * cc / (cn + cc)
    return out


# ---------------------------------------------------------------------------- common-preset steps (B5, preset common)
def weight_overrides(S, o, ctx=None):
    """override rows (o['rows']: the cell's override table; rule_set / rule_ids select) with field `comp_weights`
    (set / drop / set_keep_other) replace Gene.w through overrides.apply; shares renormalised; label override:<rule_id>
    (preset common: after nc_model, before round, so the curated value is final)"""
    from . import overrides as OV
    T = OV.check(pd.DataFrame(o["rows"]).astype(str), "data.overrides")
    T = OV.select(T, o["rule_set"], o["rule_ids"])
    T = T[T.field == "comp_weights"]

    def put(g, r, new):
        z = sum(v for v in new.values() if v > 0)
        S[g].w = {k: v / z for k, v in new.items() if v > 0}
        S[g].label = list(S[g].label) + [f"override:{r.rule_id}"]

    n = OV.apply(T, lambda kt, k: [k] if k in S else [], lambda g, f: dict(S[g].w), put)
    if ctx is not None:
        ctx.say(f"weight_overrides: {n} rows applied ({', '.join(o['rule_ids'])})")


def annotation_evidence(S, o, ctx=None):
    """Gene.order (the key list of the key-selection steps) -> one annotation source o['source'] with every key scored 1
    and the evidence class of the source the list came from (o['classes']: {source: class}; overrides and fractionation
    labels keep the class of the first label); the raw sources in o['consume'] are removed. Genes without a list keep
    their evidence. Rows of the new source are synthetic (score 1, note 'selected')."""
    import collections
    Row = collections.namedtuple("Row", ["gene", "uniprot", "source", "evidence_class", "compartment", "score",
                                         "reliability", "n_obs", "note"])
    name, cls_of, consume = o["source"], o["classes"], set(o["consume"])
    default_cls = o.get("default_class", "annotation_curated")
    for st in S.values():
        for s in list(st.src):
            if s in consume:
                del st.src[s]
        if not st.order:
            continue
        c = next((cls_of[x] for x in st.label if x in cls_of), default_cls)
        st.src[name] = {k: 1.0 for k in st.order}
        st.cls[name] = c
        st.rows[name] = [Row(st.gene, "", name, c, k, 1.0, float("nan"), 1, "selected") for k in st.order]


def resolve_imaging(S, o, ctx=None):
    """imaging pseudo-keys -> keys per protein (o['pseudo']: {key: {membrane, other, annotated}}: `membrane` for
    transmembrane proteins and proteins whose other evidence names a key in `annotated`, else `other`); imaging sources
    are dropped for proteins whose other evidence names a key in o['skip_if_keys'] (yeast: cell wall)
    # origin: yeast pipeline/inventory.py:build (CYCLoPs 'periphery', cell-wall skip) @ e3b7448"""
    img = set(o["imaging_classes"])
    pseudo = o.get("pseudo") or {}
    skip = set(o.get("skip_if_keys") or [])
    for st in S.values():
        ann = {k for s, q in st.src.items() if st.cls[s] not in img for k in q}
        for s in [s for s in st.src if st.cls[s] in img]:
            if skip & ann:
                del st.src[s]
                continue
            q, rows = {}, []
            for r in st.rows[s]:
                c = r.compartment
                p = pseudo.get(c)
                if p is not None:
                    c = p["membrane"] if (st.tm or set(p.get("annotated") or []) & ann) else p["other"]
                rows.append(r._replace(compartment=c))
            for k, v in st.src[s].items():
                p = pseudo.get(k)
                if p is not None:
                    k = p["membrane"] if (st.tm or set(p.get("annotated") or []) & ann) else p["other"]
                q[k] = max(q.get(k, 0.0), v)
            st.src[s], st.rows[s] = q, rows


def _volumes(o, ctx):
    """Vn, Vc: o['volumes'] = {nucleus: <region>, cytoplasm: <region>, measures: 'geometry'} read from the build's
    geometry/measures.json (region_volumes_mc.volume_um3), or numbers {nucleus: x, cytoplasm: y}"""
    v = o["volumes"]
    if isinstance(v.get("nucleus"), (int, float)):
        return float(v["nucleus"]), float(v["cytoplasm"]), "numbers"
    import json
    import os
    bdir = getattr(ctx, "bdir", None) or v.get("build")
    p = os.path.join(bdir, "geometry", "measures.json") if bdir else None
    if not p or not os.path.exists(p):
        raise FileNotFoundError(f"nc_model ratio: geometry measures needed for the volumes ({p}); run the geometry stage first")
    vol = json.load(open(p))["region_volumes_mc"]["volume_um3"]
    return float(vol[v["nucleus"]]), float(vol[v["cytoplasm"]]), p


def _hpa_category(st, cal):
    """v6 R1 imaging category from one imaging source's rows (main = score >= main_min): N, Nc, NC, Cn, C or None;
    None also when a main key is an organelle key (cal['skip_main_keys'])"""
    rows = [r for s in st.src if s == cal["source"] for r in st.rows[s]]
    if not rows:
        return None
    nk, ck = set(cal["nucleus_keys"]), set(cal["cytoplasm_keys"])
    main = {r.compartment for r in rows if r.score >= cal["main_min"]}
    add = {r.compartment for r in rows if r.score < cal["main_min"]} - main
    if main & set(cal.get("skip_main_keys") or []):
        return None
    nm, cm, na, ca = bool(main & nk), bool(main & ck), bool(add & nk), bool(add & ck)
    if nm and not cm and not ca:
        return "N"
    if nm and not cm:
        return "Nc"
    if nm and cm:
        return "NC"
    if cm and na:
        return "Cn"
    if cm:
        return "C"
    return None


def _imaging_ratio(st, cal, vn, vc):
    """(r, source label, {side: best imaging score}, info) from the cell's calibration, or (None, None, {}, {}); info:
    cat (HPA category), fn_img (the imaging-only nuclear share), tie_break (the fractionation tie-breaker was used)"""
    m = cal["method"]
    if st.gene in cal.get("_exclude", ()):
        return None, None, {}, {}
    info = {}
    if m == "category":
        cat = _hpa_category(st, cal)
        if cat is None:
            return None, None, {}, {}
        r = float(cal["ratios"][cat])
        info["cat"] = cat
        tb = cal.get("tie_break") or {}
        if tb and cat in set(tb.get("categories") or []):
            fr = [x for s in st.rows if st.cls[s] == tb["class"] for x in st.rows[s]]
            if fr:
                b = tb["background"]
                fn = max(fr[0].score - b, 0.0) / (1 - b)
                lo, hi = tb["clip"]
                r = float(min(max((fn / vn) / max((1 - fn) / vc, 1e-9), lo), hi)) if 0 < fn < 1 else (hi if fn >= 1 else lo)
                info["tie_break"] = True
        rows = [x for x in st.rows.get(cal["source"], [])]
    elif m == "exp":
        rows = st.rows.get(cal["source"], []) if cal["source"] in st.src else []
        sc = {r.note: r.score for r in rows}
        if cal["nucleus_term"] not in sc and cal["cytosol_term"] not in sc:
            return None, None, {}, {}
        r = float(np.exp(cal["k"] * (sc.get(cal["nucleus_term"], 0.0) - sc.get(cal["cytosol_term"], 0.0))))
    else:
        raise ValueError(f"nc_model calibration method '{m}' (allowed: category, exp, none)")
    info["fn_img"] = r * vn / (r * vn + vc)
    named = {}
    dropped = st.flags.get("bg_drop_n") or set()      # nuclear terms removed by the background term name no side
    for x in rows:
        side = "n" if x.compartment in set(cal["nucleus_keys"]) else "c" if x.compartment in set(cal["cytoplasm_keys"]) else None
        if side == "n" and cal["source"] in dropped:
            continue
        if side:
            named[side] = max(named.get(side, 0.0), x.score)
    return r, cal["source"], named, info


# ----------------------------------------------- N/C model revision 2026-11-18 (docs/unification_criteria/B5_localisation_revision_2026-11-18.md)
def chromatin_match(st, ch):
    """True when the gene is in the cell's chromatin-bound pool: ch = {regex: [patterns], feature: <feature> | None (the
    gene key)}; patterns are full-match regexes on the gene symbol (Jurkat: the histone families consumed by the
    nucleosome pre-placement + linker H1 + histone variants; yeast: the histone classes + HHO1)"""
    if not ch or not ch.get("regex"):
        return False
    import re
    if "_re" not in ch:
        ch["_re"] = re.compile("|".join(f"(?:{p})" for p in ch["regex"]))
    name = st.features.get(ch["feature"]) if ch.get("feature") else st.gene
    return isinstance(name, str) and ch["_re"].fullmatch(name) is not None


def prior_fraction(st, pr):
    """the fractionation prior's nuclear share fn_f = clip((obs - background) / (1 - background), 0, 1) of the first row of
    source pr['source'] (Itzhak 2016 N/(N+C)), or None"""
    if not pr:
        return None
    rows = st.rows.get(pr["source"]) or []
    if not rows:
        return None
    b = float(pr.get("background", 0.0))
    return min(max((rows[0].score - b) / (1 - b), 0.0), 1.0)


def ratio_of_share(fn, vn, vc, clip):
    """nuclear share of the pool -> N/C concentration ratio r = (fn / Vn) / ((1 - fn) / Vc), clipped to clip = [lo, hi]"""
    lo, hi = clip
    if fn <= 0:
        return lo
    if fn >= 1:
        return hi
    return float(min(max((fn / vn) / ((1 - fn) / vc), lo), hi))


def nc_target(r_img, fn_f, ex, o, vn, vc, tie_break=False, method="category"):
    """the revised N/C target (one function for the v0 N/C model and the re-expressed v6 R1 of preset common):
      prior (o['prior'], weight w): log r = (log r_img + w log r_f) / (1 + w) with r_f = the fractionation ratio clipped
        to prior.clip; without imaging the equal-concentration default (r = 1) takes the imaging term's place:
        log r = w log r_f / (1 + w); not applied when the imaging ratio already used the fractionation tie-breaker (NC)
      bound (o['bound'] = {max_share m, min_sources k, calibrations}): for a categorical imaging calibration (method in
        bound.calibrations; yeast's continuous CYCLoPs exp(k (N - C)) fit is not a category) fn <= m unless >= k
        independent sources are exclusive on the nucleus side (ex['n']), fn >= 1 - m unless >= k are exclusive on the
        cytoplasm side (ex['c'])
    -> (r, fn, info) or (None, None, info) when neither imaging nor the prior gives a ratio"""
    pr, bd = o.get("prior") or {}, o.get("bound") or {}
    info = {}
    r_f = ratio_of_share(fn_f, vn, vc, pr["clip"]) if (pr and fn_f is not None) else None
    if r_f is not None and tie_break:
        r_f = None                                   # NC: the fractionation already decided r (tie-breaker)
    w = float(pr.get("weight", 1.0)) if pr else 0.0
    if r_img is None and r_f is None:
        return None, None, info
    if r_f is None:
        r = r_img
    elif r_img is None:
        r = math.exp(w * math.log(r_f) / (1 + w))
        info["prior_only"] = True
    else:
        r = math.exp((math.log(r_img) + w * math.log(r_f)) / (1 + w))
    if r_f is not None:
        info["r_prior"] = r_f
    fn = r * vn / (r * vn + vc)
    if bd and method in (bd.get("calibrations") or [method]):
        m, k = float(bd["max_share"]), int(bd.get("min_sources", 2))
        fb = fn
        if ex.get("n", 0) < k:
            fn = min(fn, m)
        if ex.get("c", 0) < k:
            fn = max(fn, 1 - m)
        if fn != fb:
            info["bounded"] = (round(fb, 4), round(fn, 4))
            r = (fn / vn) / ((1 - fn) / vc) if 0 < fn < 1 else r
        info["bound_n"], info["bound_c"] = ex.get("n", 0) < k, ex.get("c", 0) < k
    return r, fn, info


def _exclusive_sources(st, nuc, cyt, img_info, fn_f, cal_src, m):
    """independent sources exclusive on one pool side: {'n': count, 'c': count} over annotation sources (pool keys all on
    one side), the imaging source (HPA category N / C, or an imaging-only nuclear share >= m / <= 1 - m) and the
    fractionation prior (fn_f >= m / <= 1 - m); the imaging and the prior are counted separately in 'img' / 'other'"""
    ex = {"n": 0, "c": 0}
    other = {"n": 0, "c": 0}
    for s, q in st.src.items():
        if s == cal_src or st.cls[s] in ("fractionation",) or st.cls[s].startswith("imaging"):
            continue
        ks = {k for k, v in q.items() if v > 0 and (k == nuc or k in cyt)}
        if ks and ks <= {nuc}:
            other["n"] += 1
        elif ks and not (ks & {nuc}):
            other["c"] += 1
    if fn_f is not None:
        if fn_f >= m:
            other["n"] += 1
        elif fn_f <= 1 - m:
            other["c"] += 1
    ex["n"], ex["c"] = other["n"], other["c"]
    cat, fi = img_info.get("cat"), img_info.get("fn_img")
    if cat == "N" or (cat is None and fi is not None and fi >= m):
        ex["n"] += 1
    if cat == "C" or (cat is None and fi is not None and fi <= 1 - m):
        ex["c"] += 1
    return ex, other


def nc_ratio(S, o, ctx=None):
    """N/C model of preset common (plan §3a step 4; revision 2026-11-18): see the module docstring. o = {pool: {nucleus,
    cytoplasm: [...], default_cytoplasm}, calibration: {method: category | exp | none, ...} (cell data), volumes, add_min,
    new_min, skip_transmembrane, chromatin, prior, bound}. Order per soluble protein with a pool share:
      1. chromatin (cell data): a chromatin-bound gene with a nucleus key gets the whole pool in the nucleus (label
         'nc:chromatin'); without a nucleus key it is unchanged
      2. r_img from the imaging calibration; prior + bound by nc_target; r = 1 when neither and both sides are present
      3. a missing side is added when imaging names it (>= add_min), its share would be >= new_min, or the bound applies
         against the present side (the bound guarantees the other side >= 1 - m)
    Gene.flags['nc'] records the per-gene inputs and result (written to localisation/nc_model.tsv)."""
    pool = o["pool"]
    nuc, cyt = pool["nucleus"], list(pool["cytoplasm"])
    cal = dict(o.get("calibration") or {"method": "none"})
    vn, vc, vsrc = _volumes(o, ctx)
    if cal.get("exclude_genes"):
        cal["_exclude"] = set(cal["exclude_genes"])
    cal.setdefault("nucleus_keys", [nuc])
    cal.setdefault("cytoplasm_keys", cyt)
    ch = dict(o.get("chromatin") or {})
    bd = o.get("bound") or {}
    m = float(bd.get("max_share", 1.0)) if bd else 1.0
    n = {"imaging": 0, "equal": 0, "added": 0, "prior": 0, "bounded": 0, "chromatin": 0}
    for st in S.values():
        if (o.get("skip_transmembrane", True) and st.tm) or not st.w:
            continue
        w = st.w
        P = w.get(nuc, 0.0) + sum(w.get(k, 0.0) for k in cyt)
        if P <= 0:
            continue
        has_n, has_c = w.get(nuc, 0.0) > 0, any(w.get(k, 0.0) > 0 for k in cyt)
        rec = dict(pool=round(P, 6))
        if chromatin_match(st, ch):
            rec["chromatin"] = True
            st.flags["nc"] = rec
            if has_n:
                st.w = {k: v for k, v in w.items() if k != nuc and k not in cyt}
                st.w[nuc] = P
                st.label = st.label + ["nc:chromatin"]
                n["chromatin"] += 1
                rec.update(r=None, fn=1.0)
            continue
        r_img, src, named, info = (None, None, {}, {}) if cal["method"] == "none" else _imaging_ratio(st, cal, vn, vc)
        fn_f = prior_fraction(st, o.get("prior"))
        ex, other = _exclusive_sources(st, nuc, cyt, info, fn_f, src, m)
        r, fn, ti = nc_target(r_img, fn_f, ex, o, vn, vc, tie_break=bool(info.get("tie_break")), method=cal["method"])
        rec.update(cat=info.get("cat"), r_img=r_img, fn_prior=fn_f, ex_n=ex["n"], ex_c=ex["c"], ex_n_other=other["n"],
                   ex_c_other=other["c"])
        if r is None:
            if not (has_n and has_c):
                st.flags["nc"] = rec
                continue
            r, lab = 1.0, "nc:equal"
            fn = r * vn / (r * vn + vc)
        else:
            parts = ([src] if r_img is not None else []) + (["prior"] if "r_prior" in ti else [])
            lab = "nc:" + "+".join(parts)
        # bounded mixture: when imaging set the ratio and the opposite side lacks k exclusive sources, the bound gives this
        # side >= 1 - m, so it may be added (no single imaging category leaves a protein at 100 % of one pool side)
        allow_n = has_n or named.get("n", 0.0) >= o["add_min"] or fn >= o["new_min"] or (bool(ti.get("bound_c")) and r_img is not None)
        allow_c = has_c or named.get("c", 0.0) >= o["add_min"] or 1 - fn >= o["new_min"] or (bool(ti.get("bound_n")) and r_img is not None)
        if not allow_n or not allow_c:
            st.flags["nc"] = rec
            continue                        # a side neither present, named nor bounded: unchanged (one side only)
        cw = {k: w[k] for k in cyt if w.get(k, 0.0) > 0} or {pool["default_cytoplasm"]: 1.0}
        cs = sum(cw.values())
        nw = {k: v for k, v in w.items() if k != nuc and k not in cyt}
        if fn > 0:
            nw[nuc] = P * fn
        if fn < 1:
            for k, v in cw.items():
                nw[k] = P * (1 - fn) * v / cs
        if (not has_n and fn > 0) or (not has_c and fn < 1):
            n["added"] += 1
        st.w = nw
        st.label = st.label + [lab]
        n["imaging" if r_img is not None else "prior" if lab != "nc:equal" else "equal"] += 1
        n["bounded"] += "bounded" in ti
        rec.update(r=r, fn=fn, label=lab, bounded=ti.get("bounded"))
        st.flags["nc"] = rec
    if ctx is not None and hasattr(ctx, "say"):
        ctx.say(f"nc_model ratio ({cal['method']}; Vn {vn:.4g} / Vc {vc:.4g} um3 from {vsrc}): imaging ratio {n['imaging']}, "
                f"prior only {n['prior']}, equal concentration {n['equal']}, pool side added {n['added']}, bounded "
                f"{n['bounded']}, chromatin {n['chromatin']}")


def imaging_background(S, o, ctx=None):
    """imaging nuclear background term b (plan §3c R8c; revision 2026-11-18): observed nuclear share = b + (1 - b) x true
    of an imaging source's split distribution -> true = clip((obs - b) / (1 - b), 0, 1); the nuclear keys
    (o['nucleus_keys']) are rescaled to `true`, the other keys of the source to 1 - true; obs 0 or 1 unchanged. Applied to
    o['sources'] only (Jurkat: HPA antibody imaging; never to live-GFP imaging such as OpenCell / CYCLoPs). Marks the
    source background-corrected (Gene.flags['bg_sources']) and records sources whose nuclear terms were removed
    (Gene.flags['bg_drop_n']: they name no nucleus side in the N/C model). o = {sources, b, nucleus_keys, provenance}"""
    b = float(o["b"])
    if not 0 <= b < 1:
        raise ValueError(f"imaging_background.b = {b}: must be in [0, 1)")
    nk, srcs = set(o["nucleus_keys"]), set(o["sources"])
    n = 0
    for st in S.values():
        for s in [s for s in st.src if s in srcs]:
            st.flags.setdefault("bg_sources", set()).add(s)
            q = st.src[s]
            obs = sum(v for k, v in q.items() if k in nk)
            tot = sum(q.values())
            if b == 0 or obs <= 0 or obs >= tot or tot <= 0:
                continue
            obs /= tot
            tru = min(max((obs - b) / (1 - b), 0.0), 1.0)
            q2 = {}
            for k, v in q.items():
                x = v / tot * (tru / obs if k in nk else (1 - tru) / (1 - obs))
                if x > 0:
                    q2[k] = x
            st.src[s] = q2
            n += 1
            if tru == 0:
                st.flags.setdefault("bg_drop_n", set()).add(s)
    if ctx is not None and hasattr(ctx, "say"):
        ctx.say(f"imaging_background: b = {b} on {sorted(srcs)} ({o.get('provenance', '')}): {n} gene distributions corrected")


def nc_table(S) -> pd.DataFrame:
    """per-gene record of the N/C model (Gene.flags['nc']): the v0 table localisation/nc_model.tsv (read by v6 R1 of preset
    common for the prior / exclusivity / chromatin inputs)"""
    cols = ["gene", "pool", "chromatin", "cat", "r_img", "fn_prior", "ex_n", "ex_c", "ex_n_other", "ex_c_other", "r", "fn",
            "label", "bounded"]
    rows = [dict(gene=g, **st.flags["nc"]) for g, st in S.items() if st.flags.get("nc")]
    T = pd.DataFrame(rows, columns=cols)
    T["chromatin"] = T["chromatin"].fillna(False).astype(bool)
    return T


STEPS = {"drop_unreliable": drop_unreliable, "restrict_membrane": restrict_membrane, "split": split, "combine": combine,
         "new_from_imaging": new_from_imaging, "round": round_shares, "nc_model": nc_model,
         "membrane_mechanism": membrane_mechanism, "ordered_annotation": ordered_annotation,
         "annotation_overrides": annotation_overrides, "mito_fractionation": mito_fractionation,
         "rank_weights": rank_weights, "imaging_groups": imaging_groups, "annotation_evidence": annotation_evidence,
         "resolve_imaging": resolve_imaging, "imaging_background": imaging_background,
         "weight_overrides": weight_overrides}
OPTIONAL = ("ordered_annotation", "annotation_overrides", "mito_fractionation", "annotation_evidence", "resolve_imaging",
            "imaging_background", "weight_overrides")


def run(E: pd.DataFrame, F: pd.DataFrame, genes, o: dict, ctx=None) -> dict:
    """the preset's combiner options -> {gene: Gene}; o['steps'] names STEPS in order, options o[<step>]"""
    for c in set(E.evidence_class):
        if c not in CLASSES:
            raise ValueError(f"evidence_class '{c}' not in {list(CLASSES)}")
    S = build_state(E, F, genes, o.get("within_source", "max"))
    for s in o["steps"]:
        if o.get(s) is None and s in OPTIONAL:      # a cell-data step the cell does not declare (preset common)
            continue
        STEPS[s](S, o.get(s), ctx)
    return S


def compartments_string(w: dict) -> str:
    """'key:w;key:w', largest first (stable for ties), as in the legacy inventory_v0.tsv"""
    return ";".join(f"{k}:{v}" for k, v in sorted(w.items(), key=lambda kv: -kv[1]))
