# origin: jurkat/validate/loc_compare.py (legacy workspace, 2026-11-15)
# Changed: steps model / compare / holdout ported (numerics unchanged); the reference-building steps (fetch, celllines,
# jurkat images) and the rule-development steps (rules, report) are NOT ported: their outputs are declared frozen
# inputs (val_*) or are inventory-rule tooling. Every list and map moved out of the code: region names and kinds from the
# build's region table, localisation-key map (SOL/TMR) from pack.region_map, region->compartment grouping, pre-pass
# library regions, consumed keys, histone families, hold-out exclusions from config validate.loc_compare; source
# vocabularies from the vocabulary file it names (config/validate/loc_vocab_human.yaml). Version switches (LOC_VERSION) are
# gone: the code is the v7 code path. Inputs are the build's own pack table, inventory_v0 and tiles; reference tables are
# frozen inputs with checksums. Added: an 'hpa_consensus_all' set in the hold-out summary (all genes with HPA data).
"""Annotation-level localisation comparison for every model gene.

  model    per-gene intended (pack table comp_weights x copies, complexes expanded, regions via the localisation-key
           map) and placed (tile job tables + pre-pass) copies per (source key, region) -> loc_model.parquet
  compare  mapping table, per-source scores, cause classes -> loc_compare.tsv, loc_mapping.tsv
  holdout  scores against reference sets (cell-line HPA annotations and its frozen hold-out half, OpenCell, HPA md5 test
           half, HPA consensus) -> loc_holdout.json, loc_holdout_pergene.tsv
"""
from __future__ import annotations

import fnmatch
import hashlib
import json
import os
import re
import time

import numpy as np
import pandas as pd
import yaml

from . import common, regions as R_


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


class Ctx:
    """configuration of one loc_compare run (all maps from config / vocabulary / region table)"""

    def __init__(self, cfg, bdir, inp: common.Inputs):
        self.cfg, self.bdir, self.inp = cfg, bdir, inp
        lc = (cfg.get("validate", {}) or {}).get("loc_compare", {}) or {}
        self.lc = lc
        voc = yaml.safe_load(open(lc["vocab"]))
        self.voc = voc
        self.COMPS = list(voc["comps"]); self.MEMCOMP = set(voc["memcomp"])
        self.SPECIAL = {k: set(v) for k, v in (voc.get("special_support") or {}).items()}
        self.SPERM = set(voc.get("sperm_terms", []))
        self.HPA_MAP = {k: tuple(v) for k, v in voc["hpa_map"].items()}
        self.HPA_EXTRA = {k: tuple(v) for k, v in (voc.get("hpa_extra_support") or {}).items()}
        self.OC_MAP = dict(voc["opencell_map"])
        self.UP_RULES = [tuple(x) for x in voc["uniprot_rules"]]
        self.UP_RULES_C = [(re.compile(a), b, c) for a, b, c in self.UP_RULES]
        self.GO_SKIP = re.compile(voc["go_skip"])
        self.REL_W = dict(voc["reliability_weight"]); self.W_UNIPROT = float(voc["w_uniprot"]); self.W_GO = float(voc["w_go"])
        self.DISAGREE = float(voc["disagree"])
        regs = R_.load_regions(cfg, inp)
        self.RN = R_.names(regs)
        self.SOLUBLE = {r["name"] for r in regs if r["kind"] in ("volume", "virtual")}
        rmap = R_.region_map(cfg, inp)
        self.SOL = {k: v["soluble"] for k, v in rmap.items() if v.get("soluble")}
        self.TMR = {k: v["tm"] for k, v in rmap.items() if v.get("tm")}
        self.SOL_DEFAULT = lc.get("soluble_default", "cytosol"); self.TM_DEFAULT = lc.get("tm_default", "ER")
        self.REG2COMP = dict(lc["reg2comp"])
        self.SRC2COMP = dict(lc["src2comp"])
        self.PM_REGS = set(lc.get("pm_regions", []))
        self.PRE_ING = dict(lc.get("preplaced_ingredient", {}) or {})
        self.PRE_REG = dict(lc.get("preplaced_region", {}) or {})
        self.PRE_DEFAULT = lc.get("preplaced_region_default", "cytosol")
        self.CONSUMED_KEY = dict(lc.get("consumed_key", {}) or {})
        self.LUMEN = lc.get("lumen_rule", {"regions": [], "field": "signal_frac", "min": 0.5, "else": "cytosol"})
        self.TMOVR = lc.get("tm_override", {}) or {}
        self.FAMILIES = list(lc.get("families", []) or [])
        self.FAM_RE = lc.get("family_regex", r"^{f}(C|-)\d")
        self.HOLDOUT_EXCLUDE = set(lc.get("holdout_exclude", []) or [])
        self.CELL_LINE = lc.get("cell_line")
        self.FRACT_TAG = lc.get("fractionation_source_tag", "")
        self.REG_DESC = dict(lc.get("mapping_descriptions", {}) or {})
        refs = lc.get("references", {}) or {}
        fi = cfg.get("frozen_inputs", {}) or {}
        self.REF = {k: (fi[v]["path"] if v in fi else v) for k, v in refs.items()}
        self.out = common.outdir(bdir)

    def pre_region(self, lib):
        for pat, reg in self.PRE_REG.items():
            if fnmatch.fnmatchcase(lib, pat):
                return reg
        return self.PRE_DEFAULT


# ----------------------------------------------------------------- model
def model_genes(C):
    I = pd.read_csv(C.inp["pack_table"], sep="\t")
    toks = {g for s in I.loc[I.copies > 0, "genes"].astype(str) for g in s.split(";")}
    V0 = pd.read_csv(C.inp["inventory_v0"], sep="\t")
    fam = {f: list(V0.gene[(V0.copies > 0) & V0.gene.astype(str).str.match(C.FAM_RE.replace("{f}", f))]) for f in C.FAMILIES}
    genes = set()
    for t in toks:
        genes |= set(fam[t]) if t in fam else {t}
    return sorted(genes), fam


def model(C, overwrite=False, nproc=8):
    inp = C.inp
    I = pd.read_csv(inp["pack_table"], sep="\t")
    V0 = pd.read_csv(inp["inventory_v0"], sep="\t").set_index("gene")
    L = json.load(open(inp["libs"])); libs = L["libs"]; consumed = L.get("consumed", {})
    f_vol, f_pm = L.get("f_vol", 1.0), L.get("f_pm", 1.0)
    S = common.tile_sums(inp, common.cache_path(C.bdir, "tile_sums.npz"), nproc=nproc)
    req, pla, nreg = S["req"], S["pla"], S["nreg"]
    PL = json.load(open(inp["preplaced_libs"])) if inp.has("preplaced_libs") else []
    lc = common.lib_cache(inp, list(libs) + list(I.ing_id) + PL, common.cache_path(C.bdir, "lib_all.npz"), nproc=nproc)
    memb = {i: (lc.present[lc.idx[i]] and lc.membrane[lc.idx[i]] != "") for i in I.ing_id}
    lidx = {l: i for i, l in enumerate(libs)}
    pc = np.bincount(np.load(inp["preplaced_types"], mmap_mode="r"), minlength=len(PL)) if PL else []
    pre = {}
    for l, c in zip(PL, pc):
        ing = C.PRE_ING.get(l, l); r_ = C.pre_region(l); pre.setdefault(ing, {}); pre[ing][r_] = pre[ing].get(r_, 0) + float(c)
    genes, fam = model_genes(C)
    famw = {f: (V0.loc[m, "copies"] / V0.loc[m, "copies"].sum()) for f, m in fam.items()}
    lumen_regs = set(C.LUMEN.get("regions", []))
    rows = []
    for r in I.itertuples():
        if r.copies <= 0: continue
        cw = json.loads(r.comp_weights) if isinstance(r.comp_weights, str) else {}
        if not cw: cw = {"unassigned": 1.0}
        cons = float(consumed.get(r.ing_id, 0)); rest = max(r.copies - cons, 0.0)
        ck = C.CONSUMED_KEY.get(r.ing_id)
        if ck and cons > 0 and rest > 0 and ck in cw:
            # the pre-pass consumption comes out of this key only (jwc consumed_key)
            C_ = float(r.copies); cw = {k: (C_ * v / rest if k != ck else max(C_ * v - cons, 0.0) / rest) for k, v in cw.items()}
        tm = bool(r.transmembrane) and memb.get(r.ing_id, False)
        ent = {}
        li = lidx.get(r.ing_id)
        for comp, wt in cw.items():
            if tm:
                reg = C.TMR.get(comp, C.TM_DEFAULT)
                ov = C.TMOVR.get(reg)
                if ov and getattr(r, ov["field"]) >= ov["min"]: reg = ov["region"]
            else:
                reg = C.SOL.get(comp, C.SOL_DEFAULT)
                if reg in lumen_regs and getattr(r, C.LUMEN["field"]) < C.LUMEN["min"]: reg = C.LUMEN["else"]
            ri = C.RN.index(reg)
            f = (pla[li, ri] / req[li, ri]) if (li is not None and req[li, ri] > 0) else 0.0
            fs = f_pm if reg in C.PM_REGS else f_vol
            x = ent.setdefault((comp, reg), [0.0, 0.0]); x[0] += rest * wt * fs; x[1] += rest * wt * f * fs
        for reg, c in pre.get(r.ing_id, {}).items():
            x = ent.setdefault(("prepass", reg), [0.0, 0.0]); x[0] += c; x[1] += c
        gs = str(r.genes).split(";"); st = [float(s) for s in str(r.stoich).split(";")]
        for g, s in zip(gs, st):
            members = famw[g].items() if g in famw else [(g, 1.0)]
            for gg, fw in members:
                for (comp, reg), (a, b) in ent.items():
                    rows.append((gg, r.ing_id, r.kind, s * fw, comp, reg, a * s * fw, b * s * fw, bool(r.transmembrane), tm,
                                 float(r.mass_da), float(r.signal_frac)))
    M = pd.DataFrame(rows, columns=["gene", "ing_id", "kind", "stoich", "src_comp", "region", "intended", "placed", "ing_tm",
                                    "ing_memlib", "ing_mass_da", "signal_frac"])
    common.write_df(os.path.join(C.out, "loc_model.parquet"), M, overwrite)
    log(f"model table: {M.gene.nunique()} genes, {M.ing_id.nunique()} ingredients, intended {M.intended.sum():,.0f} "
        f"gene copies, placed-equivalent {M.placed.sum():,.0f}")
    return M


# ----------------------------------------------------------------- scoring
def up_comp(C, head):
    h = head.strip().lower()
    for rx, c, fl in C.UP_RULES_C:
        if rx.search(h): return c, fl
    return None, ""


def parse_uniprot(C, text):
    if not isinstance(text, str): return []
    t = re.sub(r"\{[^}]*\}", "", text).replace("SUBCELLULAR LOCATION:", " ")
    t = re.sub(r"Note=.*?(?=SUBCELLULAR|$)", "", t, flags=re.S)
    out = []
    for st in re.split(r"\.\s+|\.$", t):
        st = re.sub(r"^\s*\[[^\]]*\]:\s*", "", st.strip())
        if not st: continue
        parts = [x.strip() for x in st.split(";")]
        head = parts[0]; rest = " ".join(parts[1:]).lower()
        anchor = "gpi" if "gpi-anchor" in rest or "gpi-like" in rest else "lipid" if "lipid-anchor" in rest else \
            "peripheral" if "peripheral membrane protein" in rest else "tm" if ("pass membrane" in rest) else ""
        c, fl = up_comp(C, head)
        out.append((head, c, fl, anchor))
    return out


def dist(terms):
    q = {}
    for c, w in terms:
        if c: q[c] = max(q.get(c, 0.0), w)
    return q


def score(C, p, q, main_needed=True):
    """unsup = model mass in compartments the source does not annotate; miss = share of the source's main compartments
    with < 5% of the model's copies; score = 0.5 unsup + 0.5 miss; tvd = total variation to the normalised source.
    One definition for both cells (A8 phase 2): wholecell/validate/locmetrics.score."""
    from .locmetrics import score as _score
    return _score(p, q, C.COMPS, C.SPECIAL, main_needed)


def mapping_table(C, overwrite=False):
    rows = []
    reg_of = C.REG_DESC
    H = pd.read_csv(C.REF["hpa_consensus"], sep="\t")
    hterms = sorted({t for col in ["Main location", "Additional location"] for s_ in H[col].dropna() for t in s_.split(";")})
    for t in hterms:
        c, fl = C.HPA_MAP.get(t.lower(), (None, "sperm_excluded" if t.lower() in C.SPERM else "unmapped"))
        rows.append(dict(source="HPA", term=t, comp=c or "", model_regions=reg_of.get(c, ""), flag=fl,
                         weight="main 1 / additional 0.3 x reliability (E 1, S 0.8, A 0.5, Uncertain 0)"))
    rows.append(dict(source="HPA", term="Predicted to be secreted (extracellular column)", comp="", model_regions="", flag="ignored", weight=""))
    for t, c in C.OC_MAP.items():
        rows.append(dict(source="OpenCell", term=t, comp=c, model_regions=reg_of.get(c, ""),
                         flag="no_nucleolus" if "nucleolus" in t else "", weight="grade 3: 1, grade 2: 0.6, grade 1: 0.3 (held out)"))
    for rx, c, fl in C.UP_RULES:
        rows.append(dict(source="UniProt CC / GO-CC (regex on location head)", term=rx, comp=c, model_regions=reg_of.get(c, ""), flag=fl,
                         weight=f"UniProt: every statement main, x{C.W_UNIPROT}; GO-CC: support only, x{C.W_GO}"))
    rows.append(dict(source="GO-CC", term=C.GO_SKIP.pattern, comp="", model_regions="", flag="skipped (uninformative)", weight=""))
    for k, c in C.SRC2COMP.items():
        rows.append(dict(source="model inventory comp_weights", term=k, comp=c or "pre-pass region",
                         model_regions=f"soluble: {C.SOL.get(k, '-')}; TM: {C.TMR.get(k, '-')}", flag="", weight="intended copies"))
    for r_, c in C.REG2COMP.items():
        rows.append(dict(source="model region (pack)", term=r_, comp=c, model_regions=r_, flag="", weight="placed copies"))
    T = pd.DataFrame(rows)
    common.write_df(os.path.join(C.out, "loc_mapping.tsv"), T, overwrite, index=False)
    return T


def gene_table(C, M=None):
    if M is None:
        M = pd.read_parquet(os.path.join(C.out, "loc_model.parquet"))
    V0 = pd.read_csv(C.inp["inventory_v0"], sep="\t").drop_duplicates("gene").set_index("gene")
    COMPS = C.COMPS
    M["comp"] = M.region.map(C.REG2COMP); M["soluble"] = M.region.isin(C.SOLUBLE)
    M["src_c"] = M.src_comp.map(C.SRC2COMP)
    M["lost_mem"] = M.src_c.isin(C.MEMCOMP) & M.soluble & (M.src_c != M.comp)
    gi = M.groupby(["gene", "comp"]).intended.sum().unstack(fill_value=0).reindex(columns=COMPS, fill_value=0)
    gp = M.groupby(["gene", "comp"]).placed.sum().unstack(fill_value=0).reindex(columns=COMPS, fill_value=0)
    gr = M.groupby(["gene", "region"]).intended.sum().unstack(fill_value=0)
    G = pd.DataFrame(index=gi.index)
    G["copies"] = gi.sum(1); G["copies_placed"] = gp.reindex(G.index).sum(1)
    G["placed_frac"] = G.copies_placed / G.copies
    G["mw_kda"] = V0.mass_da.reindex(G.index) / 1000
    G["mass_pg"] = G.copies * G.mw_kda * 1000 * 1.66054e-12
    G["uniprot"] = V0.uniprot.reindex(G.index); G["loc_source"] = V0.loc_source.reindex(G.index)
    G["inv_compartments"] = V0.compartments.reindex(G.index)
    for c in COMPS:
        G[f"m_{c}"] = (gi[c] / G.copies).round(4)
        G[f"pl_{c}"] = (gp.reindex(G.index)[c] / G.copies_placed.replace(0, np.nan)).round(4)
    G["model_regions"] = [";".join(f"{r_}:{v / tot:.2f}" for r_, v in row.sort_values(ascending=False).items() if v / tot >= 0.01)
                          for (g, row), tot in zip(gr.iterrows(), gr.sum(1))]
    G["soluble_frac"] = M[M.soluble].groupby("gene").intended.sum().reindex(G.index).fillna(0) / G.copies
    G["tm_ing_soluble_frac"] = M[M.soluble & M.ing_tm].groupby("gene").intended.sum().reindex(G.index).fillna(0) / G.copies
    G["lost_mem_frac"] = M[M.lost_mem].groupby("gene").intended.sum().reindex(G.index).fillna(0) / G.copies
    sol = M[M.soluble & ~M.ing_tm]
    G["unit_kda"] = (sol.intended * sol.ing_mass_da).groupby(sol.gene).sum() / sol.groupby("gene").intended.sum() / 1000
    G["n_ingredients"] = M.groupby("gene").ing_id.nunique()
    G["in_complex"] = M[M.kind != "protein"].groupby("gene").intended.sum().reindex(G.index).fillna(0) / G.copies
    return G, M


def model_p(C, G, prefix="m_"):
    return {g: {c: float(v) for c, v in zip(C.COMPS, row) if v > 0} for g, row in zip(G.index, G[[f"{prefix}{c}" for c in C.COMPS]].values)}


def hpa_table(C):
    return pd.read_csv(C.REF["hpa_consensus"], sep="\t").drop_duplicates("Gene name").set_index("Gene name")


def _hpa_terms(C, txt_main, txt_add):
    terms, flags = [], set()
    for txt, w in ((txt_main, 1.0), (txt_add, 0.3)):
        for t in str(txt).split(";") if isinstance(txt, str) else []:
            c, fl = C.HPA_MAP.get(t.lower(), (None, ""))
            if c: terms.append((c, w))
            ex = C.HPA_EXTRA.get(t.lower())
            if ex: terms.append((ex[0], ex[1] * w))
            if fl: flags.add(fl)
    return terms, flags


def sources(C, G):
    hpa = hpa_table(C) if C.REF.get("hpa_consensus") else pd.DataFrame()
    up = pd.read_csv(C.REF["uniprot"], sep="\t", usecols=["Entry", "Gene Names (primary)", "Subcellular location [CC]",
                     "Gene Ontology (cellular component)", "Transmembrane", "Intramembrane", "Lipidation", "Signal peptide", "Mass"]).set_index("Entry") \
        if C.REF.get("uniprot") else pd.DataFrame()
    oc = pd.read_csv(C.REF["opencell"]).drop_duplicates("target_name").set_index("target_name") if C.REF.get("opencell") else pd.DataFrame()
    S = {}
    for g in G.index:
        d = {}
        if g in hpa.index:
            h = hpa.loc[g]; rel = h["Reliability"]
            terms, flags = _hpa_terms(C, h["Main location"], h["Additional location"])
            d["hpa"] = dict(q=dist(terms), rel=rel, w=C.REL_W.get(rel, 0.0), main=h["Main location"], add=h["Additional location"],
                            flags=flags)
        u = G.at[g, "uniprot"]
        if isinstance(u, str) and u in up.index:
            r = up.loc[u]; st = parse_uniprot(C, r["Subcellular location [CC]"])
            q = dist([(c, 1.0) for _, c, _, _ in st])
            lip = r["Lipidation"] if isinstance(r["Lipidation"], str) else ""
            lipid = sorted({k for k, rx in (("myristoyl", "myristoyl"), ("palmitoyl", "palmitoyl"), ("prenyl", "farnesyl|geranylgeranyl"),
                                            ("GPI", "GPI-anchor amidated")) if re.search(rx, lip)})
            anchors = sorted({f"{a}@{c}" for _, c, _, a in st if a in ("lipid", "gpi", "peripheral") and c})
            go = []; gocc = str(r["Gene Ontology (cellular component)"]).lower()
            actin = bool(re.search(r"actin|cell cortex|cortical", gocc + " " + str(r["Subcellular location [CC]"]).lower()))
            if isinstance(r["Gene Ontology (cellular component)"], str):
                for t in r["Gene Ontology (cellular component)"].split(";"):
                    name = re.sub(r"\s*\[GO:\d+\]", "", t).strip().lower()
                    if not name or C.GO_SKIP.search(name): continue
                    c, _ = up_comp(C, name); go.append((c, 1.0))
            d["uniprot"] = dict(q=q, w=C.W_UNIPROT if q else 0.0, text=";".join(sorted({h_ for h_, c, _, _ in st if c})),
                                flags={fl for _, _, fl, _ in st if fl}, lipid=lipid, anchors=anchors, actin=actin,
                                tm=(isinstance(r["Transmembrane"], str) and "TRANSMEM" in r["Transmembrane"]) or
                                   (isinstance(r["Intramembrane"], str) and "INTRAMEM" in r["Intramembrane"]),
                                sp=isinstance(r["Signal peptide"], str) and len(r["Signal peptide"]) > 0)
            d["go"] = dict(q=dist(go), w=C.W_GO if go else 0.0)
        if g in oc.index:
            o = oc.loc[g]; terms = []
            for col, w in (("annotations_grade_3", 1.0), ("annotations_grade_2", 0.6), ("annotations_grade_1", 0.3)):
                if isinstance(o[col], str):
                    for t in o[col].split(";"): terms.append((C.OC_MAP.get(t.strip()), w))
            q = dist(terms)
            if q and max(q.values()) < 1:
                mx = max(q.values()); q = {c: (1.0 if w == mx else w) for c, w in q.items()}
            d["opencell"] = dict(q=q, w=1.0, text=" | ".join(f"g{k}:{o[f'annotations_grade_{k}']}" for k in (3, 2, 1)
                                                             if isinstance(o[f"annotations_grade_{k}"], str)))
        inv = G.at[g, "inv_compartments"]
        if isinstance(inv, str) and inv:
            tm = bool(d.get("uniprot", {}).get("tm")); sp = bool(d.get("uniprot", {}).get("sp")); pi = {}
            lumen_regs = set(C.LUMEN.get("regions", []))
            for kv in inv.split(";"):
                k, v = kv.rsplit(":", 1)
                reg = C.TMR.get(k, C.TM_DEFAULT) if tm else C.SOL.get(k, C.SOL_DEFAULT)
                if not tm and reg in lumen_regs and not sp: reg = C.LUMEN["else"]
                c = C.REG2COMP[reg]; pi[c] = pi.get(c, 0) + float(v)
            t_ = sum(pi.values()); d["_pinv"] = {c: v / t_ for c, v in pi.items()} if t_ else None
        S[g] = d
    return S


def cellline_refs(C):
    f = C.REF.get("hpa_celllines")
    if not f or not os.path.exists(f) or not C.CELL_LINE: return {}
    T = pd.read_csv(f, sep="\t", keep_default_na=False)
    J = T[(T.cell_line.str.upper() == C.CELL_LINE.upper())]
    out = {}
    for g, sub in J.groupby("gene"):
        locs = [t for x in sub.locations for t in x.split(";") if t]
        terms = [(C.HPA_MAP.get(t.lower(), (None, ""))[0], 1.0) for t in locs]
        terms += [(C.HPA_EXTRA[t.lower()][0], C.HPA_EXTRA[t.lower()][1]) for t in locs if t.lower() in C.HPA_EXTRA]
        out[g] = dict(q=dist([(c, w) for c, w in terms if c]), text=";".join(sorted(set(locs))),
                      rel=";".join(sorted(set(sub.ab_reliability))), n_ab=sub.antibody.nunique())
    return out


def combine(C, p, S_g):
    num = den = 0.0; per = {}
    for src in ("hpa", "uniprot", "go"):
        if src not in S_g or S_g[src]["w"] <= 0 or not S_g[src]["q"]: continue
        sc = score(C, p, S_g[src]["q"], main_needed=(src != "go")); per[src] = sc
        x = sc["unsup"] if src == "go" else sc["score"]
        num += S_g[src]["w"] * x; den += S_g[src]["w"]
    return (num / den if den else np.nan), per


def classify(C, g, G, S_g, per, comb):
    r = G.loc[g]; up = S_g.get("uniprot", {}); hp = S_g.get("hpa", {})
    p = {c: r[f"m_{c}"] for c in C.COMPS}
    contrib = {k: S_g[k]["w"] * (v["unsup"] if k == "go" else v["score"]) for k, v in per.items()}
    drv = max(contrib, key=contrib.get); q = S_g[drv]["q"]; sc = per[drv]
    main = [c for c, w in q.items() if w >= 0.999 and c in C.COMPS]
    miss_mem = [c for c in main if c in C.MEMCOMP and p.get(c, 0) < 0.05]
    anchored = bool(up.get("lipid") or up.get("anchors"))
    lipidated = bool(up.get("lipid")) or any(a.startswith(("lipid@", "gpi@")) for a in up.get("anchors", []))
    tgts = C.MEMCOMP if lipidated else ({a.split("@")[1] for a in up.get("anchors", [])} & C.MEMCOMP or C.MEMCOMP)
    if (r.tm_ing_soluble_frac >= 0.2) or (up.get("tm") and r.soluble_frac >= 0.5):
        return "c_tm_in_soluble", drv
    if anchored and r.soluble_frac >= 0.5 and set(miss_mem) & tgts:
        return "a_peripheral_or_lipid_anchored", drv
    if (p.get("nucleus", 0) < 0.05 and p.get("cytosol", 0) >= 0.8 and (r.unit_kda if r.unit_kda == r.unit_kda else r.mw_kda) < 60
            and q.get("nucleus", 0) > 0):
        return "b_small_soluble_no_nuclear_pool", drv
    pinv = S_g.get("_pinv")
    if pinv and r.in_complex >= 0.3 and sc["score"] - score(C, pinv, q, main_needed=(drv != "go"))["score"] >= 0.2:
        return "h_assembly_localisation", drv
    if C.FRACT_TAG and C.FRACT_TAG in str(r.loc_source) and set(main) <= {"nucleus", "cytosol"} and p.get("nucleus", 0) + p.get("cytosol", 0) >= 0.9:
        return "i_fractionation_nc_share", drv
    if "hpa" in per and "uniprot" in per:
        hq, uq = hp["q"], up["q"]; cs = set(hq) | set(uq)
        h_n = {c: hq.get(c, 0) / sum(hq.values()) for c in cs}; u_n = {c: uq.get(c, 0) / sum(uq.values()) for c in cs}
        if 0.5 * sum(abs(h_n[c] - u_n[c]) for c in cs) >= 0.5 and min(per["hpa"]["score"], per["uniprot"]["score"]) < C.DISAGREE:
            return "d_source_conflict", drv
    if miss_mem and r.soluble_frac >= 0.8:
        return "g_membrane_main_no_anchor_packed_soluble", drv
    if sc["unsup"] < 0.2 or len([c for c in q if c in C.COMPS]) >= 3:
        return "e_multilocalisation_shares", drv
    return "f_single_protein", drv


def hpa_split(g):
    """deterministic HPA fit/test half (md5 of the gene symbol)"""
    return "test" if int(hashlib.md5(g.encode()).hexdigest(), 16) % 2 else "fit"


def compare(C, overwrite=False):
    mapping_table(C, overwrite)
    G, M = gene_table(C); S = sources(C, G); J = cellline_refs(C)
    P = model_p(C, G); PP = model_p(C, G, "pl_")
    cl = (C.CELL_LINE or "cellline").lower()
    rows = []
    for g in G.index:
        d = S[g]; p = P[g]
        comb, per = combine(C, p, d)
        combp, _ = combine(C, PP[g], d) if PP[g] else (np.nan, {})
        r = dict(gene=g, combined_score=comb, combined_score_placed=combp, hpa_split=hpa_split(g))
        for k in ("hpa", "uniprot", "go", "opencell"):
            if k in d and d[k]["q"]:
                sc = per.get(k) or score(C, p, d[k]["q"], main_needed=(k != "go"))
                for m_ in ("unsup", "miss", "score", "tvd", "tophit"): r[f"{k}_{m_}"] = sc[m_]
        if g in J and J[g]["q"]:
            sc = score(C, p, J[g]["q"])
            for m_ in ("unsup", "miss", "score", "tvd", "tophit"): r[f"{cl}_{m_}"] = sc[m_]
        hp = d.get("hpa", {}); up = d.get("uniprot", {})
        r["actin_cortex_annot"] = bool(up.get("actin")) or "Actin filaments" in str(hp.get("main")) + str(hp.get("add"))
        r["tm_in_soluble_flag"] = bool(G.at[g, "tm_ing_soluble_frac"] >= 0.2 or (up.get("tm") and G.at[g, "soluble_frac"] >= 0.5))
        r.update(hpa_main=hp.get("main"), hpa_add=hp.get("add"), hpa_reliability=hp.get("rel"),
                 hpa_flags=";".join(sorted(hp.get("flags", set()) | up.get("flags", set()))),
                 uniprot_loc=up.get("text"), uniprot_lipid=";".join(up.get("lipid", [])),
                 uniprot_anchor=";".join(up.get("anchors", [])), uniprot_tm=up.get("tm"), uniprot_sp=up.get("sp"),
                 opencell=d.get("opencell", {}).get("text"))
        r[f"hpa_{cl}"] = J.get(g, {}).get("text"); r[f"hpa_{cl}_rel"] = J.get(g, {}).get("rel")
        if comb == comb and comb >= C.DISAGREE: r["cause"], r["driver"] = classify(C, g, G, d, per, comb)
        else: r["cause"], r["driver"] = ("agree" if comb == comb else "no_fit_source"), ""
        rows.append(r)
    R = pd.DataFrame(rows).set_index("gene")
    T = G.join(R)
    T["disagree_mass_pg"] = (T.combined_score * T.mass_pg)
    T = T.sort_values("disagree_mass_pg", ascending=False)
    common.write_df(os.path.join(C.out, "loc_compare.tsv"), T.round(4), overwrite)
    log(f"loc_compare.tsv: {len(T)} genes, {T.mass_pg.sum():.2f} pg")
    return T


# ----------------------------------------------------------------- hold-out
def _read_list(f):
    return {l.strip() for l in open(f) if l.strip() and not l.startswith("#")} if f and os.path.exists(f) else set()


def _agg(df, sc_col, mass):
    m = df[sc_col].notna(); w = mass[m]
    return dict(unweighted=round(float(df[sc_col][m].mean()), 4), mass_weighted=round(float((df[sc_col][m] * w).sum() / w.sum()), 4)) if m.any() else None


def holdout(C, overwrite=False):
    T = pd.read_csv(os.path.join(C.out, "loc_compare.tsv"), sep="\t", index_col=0)
    S = sources(C, T); J = cellline_refs(C)
    cl = (C.CELL_LINE or "cellline").lower(); CL = C.CELL_LINE or ""
    Jrel = {}
    if C.REF.get("hpa_celllines") and CL:
        Ct = pd.read_csv(C.REF["hpa_celllines"], sep="\t", keep_default_na=False)
        Jrel = Ct[Ct.cell_line == CL].groupby("gene").ab_reliability.agg(lambda x: set(x))
    hold, fit = _read_list(C.REF.get("holdout_genes")), _read_list(C.REF.get("fit_genes"))
    half = lambda g: "holdout" if g in hold else "fit" if g in fit else None
    P = model_p(C, T); PP = model_p(C, T, "pl_")
    rows = []
    for g in T.index:
        if g in C.HOLDOUT_EXCLUDE: continue
        jh = half(g)
        hq = S[g].get("hpa", {})
        for name, q in ((f"hpa_{cl}", J.get(g, {}).get("q")), ("opencell", S[g].get("opencell", {}).get("q")),
                        (f"hpa_{cl}_{jh}_half" if jh else None, J.get(g, {}).get("q")),
                        ("hpa_test_half", hq.get("q") if (hpa_split(g) == "test" and hq.get("w", 0) > 0) else None),
                        ("hpa_consensus_all", hq.get("q") if hq.get("w", 0) > 0 else None)):
            if not q or not name: continue
            if name.startswith(f"hpa_{cl}") and Jrel.get(g, set()) <= {"uncertain"}: continue
            for mode, pp in (("intended", P[g]), ("placed", PP[g])):
                if not pp: continue
                sc = score(C, pp, q)
                rows.append(dict(gene=g, set=name, mode=mode, mass_pg=T.at[g, "mass_pg"], agree=1 - sc["score"],
                                 agree_tvd=1 - sc["tvd"], tophit=sc["tophit"], supported_mass=1 - sc["unsup"], main_present=1 - sc["miss"]))
    R = pd.DataFrame(rows)
    common.write_df(os.path.join(C.out, "loc_holdout_pergene.tsv"), R, overwrite, index=False)
    out = dict(model=f"build {C.bdir} (pack table {C.inp['pack_table']}; tiles {C.inp['tiles_dir']})", date=time.strftime("%Y-%m-%d"),
               definitions=dict(
                   agree="1 - (0.5 x model share in compartments the reference does not annotate + 0.5 x share of the reference's "
                         "main compartments holding < 5% of the model's copies)",
                   agree_tvd="1 - total-variation distance between model and normalised reference distribution",
                   tophit="model's largest compartment is annotated by the reference",
                   supported_mass="model share in compartments the reference annotates",
                   main_present="share of reference main compartments with >= 5% of the model's copies",
                   compartments=C.COMPS, mass="per-gene protein mass in the model (pg)"),
               sets={})
    desc = {f"hpa_{cl}": f"HPA ICC/IF annotations made on {CL} cells (antibodies with only 'uncertain' reliability excluded)",
            "opencell": "OpenCell (HEK293T endogenous tags), grades 3/2/1 = 1/0.6/0.3; best grade = main",
            "hpa_test_half": "HPA consensus (main 1 / additional 0.3, reliability-weighted) for the md5 test half of genes",
            "hpa_consensus_all": "HPA consensus for every gene with HPA data and reliability weight > 0 (not a hold-out)",
            f"hpa_{cl}_holdout_half": f"CLEAN hold-out: HPA {CL} annotations for the frozen hold-out half (val list)",
            f"hpa_{cl}_fit_half": f"HPA {CL} annotations for the frozen fit half (not a hold-out)"}
    clean = set(C.lc.get("clean_sets", ["opencell", f"hpa_{cl}_holdout_half"]))
    for (name, mode), sub in R.groupby(["set", "mode"]):
        sub = sub.set_index("gene"); m = sub.mass_pg
        d = out["sets"].setdefault(name, dict(description=desc.get(name, name), n_genes=int(len(sub)), mass_pg=round(float(m.sum()), 3),
                                              clean=name in clean))
        d[mode] = {k: _agg(sub, k, m) for k in ("agree", "agree_tvd", "tophit", "supported_mass", "main_present")}
    out["leakage_notes"] = list(C.lc.get("leakage_notes", []) or [])
    if C.HOLDOUT_EXCLUDE:
        out["excluded_genes"] = sorted(C.HOLDOUT_EXCLUDE)
    out["references"] = {k: dict(path=v, md5=common.md5(v) if v and os.path.exists(v) else None) for k, v in C.REF.items()}
    common.write_json(os.path.join(C.out, "loc_holdout.json"), out, overwrite)
    log(json.dumps({k: {m: v[m]["agree"] for m in ("intended", "placed") if m in v} | {"n": v["n_genes"], "pg": v["mass_pg"]}
                    for k, v in out["sets"].items()}, indent=1))
    return out
