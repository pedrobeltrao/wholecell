"""Surface proteomics in the inventory stage: adapter.surface `surface_proteomics` (canonical surface.tsv) and the
factor table for an activated-cell variant of the copy-number step. OFF by default: nothing runs unless a cell config
sets `inventory.adapters.surface: surface_proteomics` (and, for the variant, adds the factor table to
`inventory.copies.per_protein_scaling`). Spec: docs/reference/inventory_schema.md §4c. copies.py is not modified: the
copy-number slot `future_sources: surface_proteomics` stays declared-not-enabled; the two supported uses go through the
existing hooks (canonical surface.tsv; copies `per_protein_scaling: factor_table`).

Input: a long protein-level table of a surface-capture experiment (e.g. glyco-capture), one row per (protein, run[,
comparison]) with a measured log2 intensity; rows repeated per statistical comparison are collapsed (the measured value of
a (protein, run) must be identical in every repeat, else ValueError); rows whose condition is not in `conditions` are
dropped (corrupt / shifted rows) and counted. Per run the values can be re-centred on the median of the proteins measured
in every run (`recentre: complete_median`), which removes a global yield drift between conditions (shift = median over those proteins of the run's
deviation from the protein's mean).

Why intensities are not used as copy numbers by default: in surface-capture data the protein intensity depends on the
number and detectability of the captured peptides per protein, so it is not comparable ACROSS proteins; per-protein
changes over conditions are the quantitative part (Jurkat glyco time course 2026, data/derived/surface_jurkat_2026/).

Modes (`inventory.params.surface_proteomics.mode`):
  presence    (default) evidence only: per protein and baseline condition, mean log2 intensity (relative units),
              detection n/N, topology (UniProt topological domains: exoplasmic side, n_term_side) -> surface.tsv;
              pm_fraction left empty.
  pm_share    as presence, plus pm_fraction from a precomputed per-protein table (`pm_share_table`: columns
              uniprot|gene, pm_fraction [, pm_fraction_lo, pm_fraction_hi, model]) produced by an offline fitted model;
              refuses to run unless the table declares `validated: true` in a '# validated: true' header line or
              `allow_unvalidated: true` is set (see DECISION 2026-11-19 in PLANNING_NOTES for the Jurkat 2026 fit).
  timepoint   rows for `timepoint.reference` and `timepoint.condition` with the per-protein log2 change in `evidence`;
              `factor_table(...)` turns the robust changes into a copies factor table (gene, factor) for an activated
              variant (copies: per_protein_scaling: [{kind: factor_table, rule_id: surface_<cond>, path: <table>}]).

Params (defaults in DEFAULTS):
  path            frozen input name | path of the long table (confidential data: declare as frozen input with
                  redistributable: no; never commit the data)
  source_id       label in surface.tsv
  columns         {protein, run, condition, replicate, value, features, gene} -> column names of the long table
  conditions      list of valid condition labels (others dropped) | null (keep all)
  recentre        complete_median | none
  baseline        condition(s) describing the cell of the build (Jurkat resting: [NAct, Act_0])
  min_detected    replicates (per condition) needed to call a protein detected (default 3)
  annotation      UniProt TSV (Entry, Gene Names (primary), Topological domain, Transmembrane, Lipidation, Signal peptide)
                  for gene names and topology; null = genes from columns.gene, topology empty
  timepoint       {condition, reference, min_abs_log2: 0.585, min_pairs_agree: 3, min_detected: 3, pm_only_genes: null,
                   stats_table: null (optional per-protein q-values: columns uniprot, contrast, q), max_q: 0.05}
"""
from __future__ import annotations

import copy
import re

import numpy as np
import pandas as pd

from ..registry import register
from . import schema

DEFAULTS = {
    "path": None,
    "source_id": "surface_proteomics",
    "columns": {"protein": "entry", "run": "run", "condition": "condition", "replicate": "bio_replicate",
                "value": "normalised_protein_intensity_log2", "features": "n_features_per_sample", "gene": None},
    "conditions": None,
    "recentre": "complete_median",
    "baseline": None,
    "min_detected": 3,
    "mode": "presence",
    "pm_share_table": None,
    "allow_unvalidated": False,
    "annotation": None,
    "timepoint": {"condition": None, "reference": None, "min_abs_log2": 0.585, "min_pairs_agree": 3, "min_detected": 3,
                  "pm_only_genes": None, "stats_table": None, "max_q": 0.05},
}
MODES = ("presence", "pm_share", "timepoint")
UNIT = "log2 intensity (surface capture; relative within a protein, not comparable across proteins)"


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def options(params: dict | None) -> dict:
    o = _merge(DEFAULTS, params or {})
    if o["mode"] not in MODES:
        raise ValueError(f"surface_proteomics: unknown mode {o['mode']!r} (one of {MODES})")
    if o["recentre"] not in ("complete_median", "none"):
        raise ValueError(f"surface_proteomics: unknown recentre {o['recentre']!r}")
    return o


# ------------------------------------------------------------------------------------------------ reading
def read_long(path: str, o: dict):
    """-> (M: protein x run log2, F: protein x run features | None, samples: run, condition, replicate; info dict)"""
    c = o["columns"]
    d = pd.read_csv(path, low_memory=False)
    n0 = len(d)
    if o["conditions"]:
        d = d[d[c["condition"]].isin(o["conditions"])]
    info = {"rows": n0, "rows_dropped_condition": int(n0 - len(d))}
    d = d.copy()
    d[c["value"]] = pd.to_numeric(d[c["value"]], errors="coerce")
    meas = d.dropna(subset=[c["value"]])
    nun = meas.groupby([c["protein"], c["run"]])[c["value"]].nunique()
    if (nun > 1).any():
        bad = nun[nun > 1].index[:3].tolist()
        raise ValueError(f"surface_proteomics: measured value differs between repeated rows for {bad}")
    M = meas.groupby([c["protein"], c["run"]])[c["value"]].first().unstack()
    F = None
    if c.get("features") and c["features"] in d:
        meas = meas.assign(_f=pd.to_numeric(meas[c["features"]], errors="coerce"))
        F = meas.groupby([c["protein"], c["run"]])._f.first().unstack().reindex(index=M.index, columns=M.columns)
    S = d[[c["run"], c["condition"], c["replicate"]]].drop_duplicates()
    S.columns = ["run", "condition", "replicate"]
    if S.run.duplicated().any():
        raise ValueError("surface_proteomics: a run maps to several conditions / replicates")
    S = S.set_index("run").reindex(M.columns)
    genes = None
    if c.get("gene") and c["gene"] in d:
        genes = d.groupby(c["protein"])[c["gene"]].first()
    info.update(proteins=int(M.shape[0]), runs=int(M.shape[1]), pairs_measured=int(M.notna().sum().sum()))
    return M, F, S, genes, info


def recentre(M: pd.DataFrame, how: str):
    """shift each run by the median, over the proteins measured in every run, of (value - that protein's mean over runs),
    centred on the median run (-> M, shift)"""
    if how == "none":
        return M, pd.Series(0.0, index=M.columns)
    full = M.dropna()
    if len(full) < 20:
        raise ValueError(f"surface_proteomics: only {len(full)} proteins measured in every run; cannot re-centre")
    med = full.sub(full.mean(axis=1), axis=0).median()
    shift = med - med.median()
    return M - shift, shift


# ------------------------------------------------------------------------------------------------ topology


def infer_sides(L, tm, topo, sp_end=0, gpi=False, nglyc=()):
    """per-residue membrane side: 'e' exoplasmic, 'c' cytoplasmic, 'm' TM helix, '' unknown. Non-TM segments between TM
    helices alternate sides; the phase comes from (in order) annotated TOPO_DOM (Extracellular / Lumenal = e,
    Cytoplasmic = c; majority vote over segments, parity-propagated), annotated N-glycosylation sites (exoplasmic), a
    signal peptide (first segment exoplasmic); GPI-anchored or secreted (signal peptide, no TM): mature chain exoplasmic."""
    side = [""] * L
    tm = sorted((a, b) for a, b in tm if 1 <= a <= b <= L)
    for a, b in tm:
        for i in range(a - 1, b):
            side[i] = "m"
    if not tm:
        if gpi or sp_end:
            for i in range(sp_end, L):
                side[i] = "e"
        for a, b, note in topo:     # soluble with annotation (rare)
            s = "e" if note.startswith(("Extracellular", "Lumenal")) else ("c" if note.startswith("Cytoplasmic") else "")
            for i in range(max(a - 1, 0), min(b, L)):
                side[i] = s or side[i]
        return side
    bounds, prev = [], 0
    for a, b in tm:
        bounds.append((prev, a - 1)); prev = b
    bounds.append((prev, L))                         # segment k = residues [start, end) 0-based
    vote = 0                                         # + : segment 0 exoplasmic
    for a, b, note in topo:
        s = 1 if note.startswith(("Extracellular", "Lumenal")) else (-1 if note.startswith("Cytoplasmic") else 0)
        if not s: continue
        mid = (a + b) / 2 - 1
        for k, (x, y) in enumerate(bounds):
            if x <= mid < y:
                vote += s * (b - a + 1) * (1 if k % 2 == 0 else -1)
    src = "topo"
    if vote == 0:
        for p in nglyc:
            for k, (x, y) in enumerate(bounds):
                if x <= p - 1 < y:
                    vote += 1 if k % 2 == 0 else -1
        src = "nglyc"
    if vote == 0 and sp_end:
        vote = 1; src = "signal"
    if vote == 0:
        return side
    for k, (x, y) in enumerate(bounds):
        exo = (k % 2 == 0) == (vote > 0)
        for i in range(max(x, sp_end if k == 0 else x), y):
            side[i] = "e" if exo else "c"
    return side


def _ranges(field, key):
    if not isinstance(field, str):
        return []
    return [(int(a), int(b), (n or "").split(";")[0]) for a, b, n in
            re.findall(key + r" <?(\d+)\.\.>?(\d+);(?: /note=\"([^\"]*)\")?", field)]


def _topology_rows(u: pd.DataFrame) -> pd.DataFrame:
    col = {k.lower(): k for k in u.columns}
    get = lambda r, name: r.get(col.get(name.lower()), np.nan)  # noqa: E731
    out = []
    for r in u.to_dict("records"):
        g = get(r, "Gene Names (primary)"); L = get(r, "Length"); seq = get(r, "Sequence")
        doms = _ranges(get(r, "Topological domain"), "TOPO_DOM")
        tms = [(a, b) for a, b, _ in _ranges(get(r, "Transmembrane"), "TRANSMEM")]
        lip = get(r, "Lipidation"); gpi = isinstance(lip, str) and "GPI-anchor" in lip
        sp = _ranges(get(r, "Signal peptide"), "SIGNAL"); sp_end = sp[0][1] if sp else 0
        gly = get(r, "Glycosylation")
        ng = [int(m) for m in re.findall(r"CARBOHYD (\d+); /note=\"N-linked", gly)] if isinstance(gly, str) else []
        if not isinstance(L, (int, float)) or pd.isna(L):
            L = len(seq) if isinstance(seq, str) else max([b for _, b in tms] + [b for _, b, _ in doms] + [sp_end, 1])
        sd = infer_sides(int(L), tms, doms, sp_end, gpi, ng)
        first = next((x for x in sd if x in ("e", "c")), "")
        nts = {"e": "extracellular", "c": "cytoplasmic"}.get(first, "")
        n_tm = len(tms)
        if gpi:
            cls = "GPI"
        elif n_tm == 0:
            cls = "secreted_or_peripheral" if sp_end else "no_TM"
        elif n_tm == 1:
            cls = {"extracellular": "single_pass_typeI", "cytoplasmic": "single_pass_typeII"}.get(nts, "single_pass")
        else:
            cls = f"multi_pass_{n_tm}"
        exo = []
        i = 0
        while i < len(sd):
            if sd[i] == "e":
                j = i
                while j < len(sd) and sd[j] == "e":
                    j += 1
                exo.append(f"{i + 1}-{j}"); i = j
            else:
                i += 1
        out.append(dict(uniprot=r[col["entry"]], gene=(str(g).split("; ")[0] if isinstance(g, str) else np.nan),
                        topology=cls, n_term_side=nts, exo_segments=";".join(exo),
                        exo_frac=(sum(x == "e" for x in sd) / len(sd)) if sd else np.nan))
    return pd.DataFrame(out).drop_duplicates("uniprot").set_index("uniprot")


def topology_table(path: str) -> pd.DataFrame:
    """UniProt TSV -> per Entry: gene, topology class, n_term_side, exoplasmic segments and fraction (infer_sides: TM
    helices alternate sides, phase from TOPO_DOM, else N-glycosylation sites, else signal peptide; GPI-anchored or
    signal peptide without TM: the mature chain). Partial annotations are completed (UniProt TRAC / TRBC1 / TRBC2 only
    annotate the cytoplasmic tail; their constant domains are exoplasmic and glyco-captured)."""
    return _topology_rows(pd.read_csv(path, sep="\t"))


# ------------------------------------------------------------------------------------------------ core table
def _cond_stats(M: pd.DataFrame, S: pd.DataFrame, cond: str):
    r = S.index[S.condition == cond]
    if not len(r):
        raise ValueError(f"surface_proteomics: condition {cond!r} not in the data ({sorted(S.condition.unique())})")
    return M[r].mean(axis=1), M[r].notna().sum(axis=1), len(r)


def paired_change(M: pd.DataFrame, S: pd.DataFrame, cond: str, ref: str) -> pd.DataFrame:
    """per protein: log2 change cond - ref (difference of means), paired replicate differences agreeing in sign."""
    a, na, Na = _cond_stats(M, S, cond)
    b, nb, Nb = _cond_stats(M, S, ref)
    ra = S[S.condition == cond].reset_index().set_index("replicate").run
    rb = S[S.condition == ref].reset_index().set_index("replicate").run
    reps = sorted(set(ra.index) & set(rb.index))
    D = pd.DataFrame({k: M[ra[k]] - M[rb[k]] for k in reps})
    lfc = a - b
    agree = (np.sign(D).eq(np.sign(lfc), axis=0) & D.notna()).sum(axis=1)
    return pd.DataFrame({"log2_change": lfc, "n_cond": na, "n_ref": nb, "N_cond": Na, "N_ref": Nb,
                         "n_pairs": D.notna().sum(axis=1), "n_pairs_agree": agree})


def surface_table(ctx_resolve, params: dict) -> tuple[pd.DataFrame, dict]:
    """the canonical surface.tsv rows (schema.SURFACE) + info; ctx_resolve maps file references to paths"""
    o = options(params)
    if not o["path"]:
        raise ValueError("surface_proteomics: params.path is not set")
    M, F, S, genes, info = read_long(ctx_resolve(o["path"]), o)
    M, shift = recentre(M, o["recentre"])
    info["recentre_shift_by_condition"] = shift.groupby(S.condition).mean().round(4).to_dict()
    topo = topology_table(ctx_resolve(o["annotation"])) if o["annotation"] else None
    gene = pd.Series(np.nan, index=M.index, dtype=object)
    if topo is not None:
        gene = gene.fillna(topo.gene.reindex(M.index))
    if genes is not None:
        gene = gene.fillna(genes.reindex(M.index))
    feats = F.max(axis=1) if F is not None else pd.Series(np.nan, index=M.index)
    if o["mode"] == "timepoint":
        tp = o["timepoint"]
        if not (tp["condition"] and tp["reference"]):
            raise ValueError("surface_proteomics: mode timepoint needs timepoint.condition and timepoint.reference")
        conds = [tp["reference"], tp["condition"]]
        ch = paired_change(M, S, tp["condition"], tp["reference"])
    else:
        conds = o["baseline"] or sorted(S.condition.unique())
        conds = [conds] if isinstance(conds, str) else list(conds)
        ch = None
    pmf = None
    if o["mode"] == "pm_share":
        pmf = read_pm_share(ctx_resolve(o["pm_share_table"]) if o["pm_share_table"] else None, o["allow_unvalidated"])
    rows = []
    for cnd in conds:
        mu, n, N = _cond_stats(M, S, cnd)
        for p in M.index:
            if n[p] == 0:
                continue
            det = n[p] >= min(o["min_detected"], N)
            ev = [f"detected {n[p]}/{N}" + ("" if det else " (below min_detected)")]
            if not np.isnan(feats[p]):
                ev.append(f"max features {int(feats[p])}")
            if ch is not None:
                c = ch.loc[p]
                if not np.isnan(c.log2_change):
                    ev.append(f"log2 {o['timepoint']['condition']}/{o['timepoint']['reference']} {c.log2_change:+.3f} "
                              f"(pairs agree {int(c.n_pairs_agree)}/{int(c.n_pairs)})")
            frac = np.nan
            if pmf is not None:
                key = p if p in pmf.index else gene.get(p)
                frac = pmf.pm_fraction.get(key, np.nan)
            t = topo.loc[p] if (topo is not None and p in topo.index) else None
            rows.append(dict(uniprot=p, gene=gene.get(p), source_id=o["source_id"], condition=cnd, replicate="mean",
                             pm_abundance=float(mu[p]), unit=UNIT, pm_fraction=frac,
                             topology=(t.topology if t is not None else ""),
                             n_term_side=(t.n_term_side if t is not None else ""),
                             evidence="; ".join(ev), notes=("detected" if det else "weak") + f"; mode {o['mode']}"))
    out = pd.DataFrame(rows, columns=schema.SURFACE)
    info.update(mode=o["mode"], conditions=conds, rows_out=int(len(out)),
                detected={c: int(out[(out.condition == c) & out.notes.str.startswith("detected")].shape[0]) for c in conds})
    return out, info


def read_pm_share(path, allow_unvalidated=False) -> pd.DataFrame:
    if not path:
        raise ValueError("surface_proteomics: mode pm_share needs pm_share_table")
    with open(path) as fh:
        head = [ln for ln in (fh.readline() for _ in range(20)) if ln.startswith("#")]
    validated = any(re.match(r"#\s*validated:\s*true", ln, re.I) for ln in head)
    if not validated and not allow_unvalidated:
        raise ValueError(f"surface_proteomics: {path} is not marked '# validated: true'; the PM-share model must pass its "
                         "cross-validation before it sets fractions (set allow_unvalidated: true only for diagnostics)")
    t = pd.read_csv(path, sep="\t", comment="#")
    key = "uniprot" if "uniprot" in t else "gene"
    if "pm_fraction" not in t:
        raise ValueError(f"surface_proteomics: {path} has no pm_fraction column")
    if ((t.pm_fraction < 0) | (t.pm_fraction > 1)).any():
        raise ValueError("surface_proteomics: pm_fraction outside [0, 1]")
    return t.set_index(key)


# ------------------------------------------------------------------------------------------------ activated variant
def factor_table(ctx_resolve, params: dict) -> pd.DataFrame:
    """robust per-protein surface changes of timepoint.condition vs timepoint.reference -> (gene, uniprot, factor,
    log2_change, evidence) for copies per_protein_scaling: factor_table. A protein qualifies when it is detected
    (>= min_detected) in both conditions, |log2 change| >= min_abs_log2, >= min_pairs_agree paired replicate differences
    agree in sign, (stats_table given) q <= max_q, and (pm_only_genes given) its gene is in that list (e.g. genes with
    >= 75% of their copies at the PM, so a surface change stands for a whole-cell change)."""
    o = options(_merge(params, {"mode": "timepoint"}))
    tp = o["timepoint"]
    M, F, S, genes, info = read_long(ctx_resolve(o["path"]), o)
    M, _ = recentre(M, o["recentre"])
    ch = paired_change(M, S, tp["condition"], tp["reference"])
    gene = pd.Series(np.nan, index=M.index, dtype=object)
    if o["annotation"]:
        gene = gene.fillna(topology_table(ctx_resolve(o["annotation"])).gene.reindex(M.index))
    if genes is not None:
        gene = gene.fillna(genes.reindex(M.index))
    ok = ((ch.n_cond >= tp["min_detected"]) & (ch.n_ref >= tp["min_detected"]) &
          (ch.log2_change.abs() >= tp["min_abs_log2"]) & (ch.n_pairs_agree >= tp["min_pairs_agree"]))
    if tp.get("stats_table"):
        st = pd.read_csv(ctx_resolve(tp["stats_table"]), sep="\t")
        st = st[st.contrast == f"{tp['condition']}_vs_{tp['reference']}"].set_index("uniprot").q
        ok &= st.reindex(ch.index).le(tp["max_q"]).fillna(False)
    if tp.get("pm_only_genes") is not None:
        ok &= gene.isin(set(tp["pm_only_genes"]))
    x = ch[ok]
    return pd.DataFrame({"gene": gene[x.index].values, "uniprot": x.index, "factor": 2.0 ** x.log2_change.values,
                         "log2_change": x.log2_change.values,
                         "evidence": [f"surface {tp['condition']} vs {tp['reference']}: pairs agree {a}/{n}"
                                      for a, n in zip(x.n_pairs_agree, x.n_pairs)]}).dropna(subset=["gene"]).reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ adapter
@register("adapter.surface", "surface_proteomics")
def surface_proteomics(ctx) -> pd.DataFrame:
    out, info = surface_table(ctx.resolve, ctx.params("surface_proteomics"))
    ctx.report.append(f"surface_proteomics: {info}")
    return out
