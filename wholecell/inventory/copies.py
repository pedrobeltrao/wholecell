"""Generic copy-number step `copies` (cell-agnostic, config-driven): abundance sources -> copies per cell.

    python -m wholecell inventory --cell C --build B --step copies        -> inventory/copies/*

One step for every cell (human Jurkat, yeast, others); everything cell-specific is in `inventory.copies` of the cell config
(spec: docs/reference/inventory_schema.md §4b; README: wholecell/inventory/COPIES.md). Pipeline, in this order:

  1 annotation    protein universe candidates (key, uniprot, mass_da [, gene, orf]); reader by name
  2 sources       typed abundance sources (relative_ibaq | ppb | ppm | lfq | absolute_copies), each with its reader, replicate
                  aggregation and protein-group splitting rule; values in the source's own unit (`column`), n_rep
  3 universe      annotation keys, optionally restricted to keys with a source value or an RNA row; sorted or annotation order
  4 combine       which source sets each protein: `priority` (first source with a value) or `geometric_mean`; secondary
                  sources are put on the primary's scale (`harmonise: median_ratio`) -> value column + `source` label + set_by
  5 imputation    none | rna_loglog (log10 value ~ log10 TPM fitted on detected proteins, capped at a percentile of the
                  detected values, `low` calls x low_factor, off / no RNA -> 0); otherwise unquantified proteins -> 0
  6 absent        zero lists with reasons (gene lists, regex patterns), applied before the mass scaling
  7 per-protein   optional factors on the relative values: size_slope (value x ratio^slope; yeast G1 scaling, Lanz 2024)
                  or factor_table
  8 total mass    molar fraction x N, N such that sum(copies x mass) = total protein (pg): mode fixed (value, or
                  inventory.total_protein_pg), anchors (implied totals of structural anchors combined by a rule), or
                  source_sum (an absolute source keeps its mass)
  9 adjustments   optional, after the scaling: family_target (scale a protein set so that its family statistic equals a
                  structure count, e.g. yeast ribosomes = cryo-ET density x cytosol volume; the other proteins keep theirs)
 10 rounding      optional (Jurkat: integers; yeast: floats)
 11 anchor check  diagnostic: per anchor (histones vs nucleosomes, nucleoporins vs pores, r-proteins vs ribosomes ...)
                  the family ratios copies / (stoich x count) and the implied total protein = total / median ratio
 12 ledger        copy_ledger.tsv: one row per protein and rule that set or changed its value (rule ids), global rows '*'

Outputs (inventory/copies/): copies.tsv (one row per protein: annotation, per-source values, imputation inputs, `source`
label, value, note, factors, copies, set_by, rule_ids), copy_ledger.tsv, copy_sources.tsv (declared sources incl. future
slots and their status), anchors.tsv (per anchor x count x family), summary.json, report.txt.
The proteome adapter `copies` (registered here) hands copies.tsv to the v0 step (columns `inventory.copies.v0_columns`).

Origins: the Jurkat recipe is the arithmetic of plugins/jurkat/inventory/inventory_v0.py (legacy jurkat/inventory/
inventory_v0.py); the yeast recipe (PaxDb ppm, NOT_IN_THIS_CELL, G1 size scaling, ribosomes from the cryo-ET density,
anchors) follows yeast_cell_viewer pipeline/inventory.py, ingredients.py and validation/anchors.py at commit 34d123f
(G. Magnusson; read-only reference, no licence yet: code derived from it is marked "origin: yeast_cell_viewer").
"""
from __future__ import annotations

import json
import re

import numpy as np
import pandas as pd

from .. import registry
from ..registry import register
from . import schema

OUTPUTS = ["copies/copies.tsv", "copies/copy_ledger.tsv", "copies/copy_sources.tsv", "copies/anchors.tsv",
           "copies/summary.json", "copies/report.txt"]
LEDGER = ["key", "rule_id", "field", "old", "new", "factor", "evidence"]
ANCHORS = ["anchor_id", "structure", "count_label", "count", "family", "members", "stoich", "copies", "need", "ratio"]
SOURCE_TYPES = {
    "relative_ibaq": "molar (iBAQ ~ moles): value used as is",
    "ppb": "molar parts per billion: value used as is",
    "ppm": "molar parts per million: value used as is",
    "lfq": "LFQ / summed intensity ~ mass: value / mass_da (molar)",
    "absolute_copies": "copies per cell: value used as is (molar); keeps its mass with total_protein.mode source_sum",
}
FUTURE_TYPES = {"secms_abundance": "SEC-MS summed intensity (cross-check of the relative abundances)",
                "surface_proteomics": "surface proteomics (plasma-membrane copy numbers)"}


def _cfg(ctx):
    c = ctx.inv.get("copies")
    if not c:
        raise ValueError("inventory.copies is not configured")
    return c


def _read_tsv(path, **kw):
    return pd.read_csv(path, sep="\t", comment=kw.pop("comment", "#"), **kw)


def _fmt(x):
    if isinstance(x, (float, np.floating)):
        return "" if np.isnan(x) else repr(float(x))
    return "" if x is None else str(x)


# ------------------------------------------------------------------------------------------------ readers (registered)
def reader(name):
    return registry.get("copies.reader", name)


def uniprot_reviewed_genes(ctx, spec):
    """UniProt TSV export keyed by primary gene name (reviewed entries; first of several names; first entry per gene).
    Columns: Entry, Reviewed, Gene Names (primary), Mass. Origin: plugins/jurkat/inventory/inventory_v0.py::_uniprot."""
    up = pd.read_csv(ctx.resolve(spec["path"]), sep="\t")
    up = up[up.Reviewed == "reviewed"].rename(columns={"Gene Names (primary)": "gene"})
    up = up.dropna(subset=["gene"]).drop_duplicates("gene")
    up["gene"] = up.gene.str.split("; ").str[0]
    up = up.drop_duplicates("gene").set_index("gene")
    out = pd.DataFrame(index=up.index)
    out["uniprot"] = up.Entry
    out["mass_da"] = up.Mass.astype(float)
    out["gene_name"] = up.index
    return out


register("copies.reader", "uniprot_reviewed_genes")(uniprot_reviewed_genes)


def uniprot_ordered_locus(ctx, spec):
    """UniProt TSV export keyed by gene name with the ordered-locus name (ORF) as fallback, one row per ORF (first token
    of 'Gene Names (ordered locus)'); column `orf` carried for ORF-keyed sources. Positional columns as in the yeast
    export (Entry, Reviewed, ordered locus, primary, Length, Mass, ...) unless spec.columns names them.
    origin: yeast_cell_viewer pipeline/inventory.py build() (34d123f), the UniProt part."""
    u = pd.read_csv(ctx.resolve(spec["path"]), sep="\t")
    cols = spec.get("columns") or {"uniprot": 0, "orf": 2, "gene": 3, "mass_da": 5}
    pick = {k: (u.columns[v] if isinstance(v, int) else v) for k, v in cols.items()}
    out = pd.DataFrame({k: u[c] for k, c in pick.items()})
    out["orf"] = out.orf.astype(str).str.split(r"[ ;]").str[0]
    out = out.drop_duplicates("orf")
    out["gene_name"] = out.gene                      # UniProt primary name (NaN for unnamed ORFs)
    out["gene"] = out.gene.fillna(out.orf)
    if out.gene.duplicated().any():
        raise ValueError(f"uniprot_ordered_locus: duplicated keys {out.gene[out.gene.duplicated()].tolist()[:5]}")
    out["mass_da"] = out.mass_da.astype(float)
    return out.set_index("gene")[["uniprot", "mass_da", "orf", "gene_name"]]


register("copies.reader", "uniprot_ordered_locus")(uniprot_ordered_locus)


def annotation_table(ctx, spec):
    """plain TSV: columns key (default gene), uniprot, mass_da [, orf]"""
    T = _read_tsv(ctx.resolve(spec["path"]))
    k = spec.get("key_column", "gene")
    T = T.drop_duplicates(k).set_index(k)
    T.index.name = "gene"
    T["mass_da"] = T.mass_da.astype(float)
    if "gene_name" not in T.columns:
        T["gene_name"] = T.index
    return T


register("copies.reader", "annotation_table")(annotation_table)


def _split_group(genes, rule):
    """protein-group splitting: which of a group's genes share its value (equally)"""
    rule = rule or {"rule": "first"}
    r = rule.get("rule", "first")
    if r == "family_prefix":          # Jurkat: genes sharing the first gene's n-letter prefix (histone / tubulin families)
        n = int(rule.get("n", 3))
        return [x for x in genes if x[:n] == genes[0][:n]]
    if r == "first":
        return genes[:1]
    if r == "equal_all":
        return list(genes)
    raise ValueError(f"unknown protein-group split rule {r}")


def _aggregate(P, how):
    """replicate aggregation of a gene x replicate frame -> (value, n_rep)"""
    nrep = P.notna().sum(axis=1)
    if how in ("median_detected", "median"):
        v = P.median(axis=1, skipna=True)
    elif how in ("mean_detected", "mean"):
        v = P.mean(axis=1, skipna=True)
    elif how == "geomean_detected":
        v = np.exp(np.log(P.where(P > 0)).mean(axis=1, skipna=True))
    elif how == "single":
        v = P.iloc[:, 0]
    else:
        raise ValueError(f"unknown replicate aggregation {how}")
    return v, nrep


def maxquant_protein_groups(ctx, spec, ann):
    """MaxQuant proteinGroups-like table: one row per protein group (id column of ';'-separated accessions, isoform
    suffixes dropped) x replicate columns (names containing spec.sample_match, or spec.sample_columns). Accessions map to
    keys via the annotation's uniprot column; a group's per-replicate values are split equally among the genes chosen by
    spec.group_split; values summed per gene (min_count 1), then aggregated over replicates (spec.aggregate).
    Arithmetic of plugins/jurkat/inventory/inventory_v0.py (Geiger 2012 iBAQ ppb)."""
    g = pd.read_csv(ctx.resolve(spec["path"]), sep="\t")
    cols = spec.get("sample_columns") or [c for c in g.columns if spec["sample_match"] in c]
    acc2gene = dict(zip(ann.uniprot, ann.index))
    rows = []
    for ids, vals in zip(g[spec.get("id_column", "Majority.protein.IDs")], g[cols].apply(pd.to_numeric, errors="coerce").values):
        if np.all(np.isnan(vals)) or np.nansum(vals) <= 0:
            continue
        gs = list(dict.fromkeys(acc2gene[a.split("-")[0]] for a in str(ids).split(";") if a.split("-")[0] in acc2gene))
        if not gs:
            continue
        fam = _split_group(gs, spec.get("group_split"))
        for x in fam:
            rows.append((x, *(vals / len(fam))))
    P = pd.DataFrame(rows, columns=["gene"] + cols).groupby("gene").sum(min_count=1)
    v, nrep = _aggregate(P, spec.get("aggregate", "median_detected"))
    v = v[v > 0]
    ctx.say(f"source {spec['id']}: protein groups mapped: {len(rows)} genes: {len(v)}")
    return v, nrep


register("copies.reader", "maxquant_protein_groups")(maxquant_protein_groups)


def paxdb(ctx, spec, ann):
    """PaxDb abundance file (columns gene/internal id, string id '<taxon>.<id>', abundance; '#' header lines). The id after
    the taxon prefix maps to the annotation column spec.map_on (default orf); duplicated ids are summed.
    origin: yeast_cell_viewer pipeline/inventory.py build() (34d123f), the PaxDb part."""
    pax = pd.read_csv(ctx.resolve(spec["path"]), sep="\t", comment="#", names=["pax_gene", "string_id", "value"])
    pre = spec.get("strip_prefix")
    pax["id"] = pax.string_id.str.replace(pre, "", regex=False) if pre else pax.string_id
    pax = pax.groupby("id", as_index=False).value.sum()
    on = spec.get("map_on", "orf")
    key_of = dict(zip(ann[on], ann.index)) if on != "key" else {k: k for k in ann.index}
    pax["key"] = pax.id.map(key_of)
    unm = pax.key.isna()
    ctx.say(f"source {spec['id']}: {len(pax)} ids, {int((~unm).sum())} matched to the annotation "
            f"({pax.value[~unm].sum():.6g} of {pax.value.sum():.6g} {spec.get('unit', '')})")
    m = pax[~unm]
    v = m.groupby("key").value.sum()
    return v, pd.Series(1, index=v.index)


register("copies.reader", "paxdb")(paxdb)


def abundance_table(ctx, spec, ann):
    """plain TSV: id column (spec.id_column, mapped through spec.map_on of the annotation, default the key) and one or
    more value columns (spec.value_columns, or spec.value_regex): replicates aggregated by spec.aggregate"""
    T = _read_tsv(ctx.resolve(spec["path"]))
    idc = spec.get("id_column", "gene")
    vc = spec.get("value_columns") or [c for c in T.columns if re.search(spec.get("value_regex", r"^value$"), c)]
    on = spec.get("map_on", "key")
    key_of = {k: k for k in ann.index} if on == "key" else dict(zip(ann[on], ann.index))
    T["key"] = T[idc].map(key_of)
    T = T.dropna(subset=["key"])
    P = T.groupby("key")[vc].sum(min_count=1)
    v, nrep = _aggregate(P.apply(pd.to_numeric, errors="coerce"), spec.get("aggregate", "median_detected"))
    v = v[v > 0]
    ctx.say(f"source {spec['id']}: {len(v)} proteins with a value ({len(vc)} replicate columns)")
    return v, nrep


register("copies.reader", "table")(abundance_table)


# -------------------------------------------------------------------------------------------------------- steps
class _Ledger:
    def __init__(self):
        self.rows = []

    def add(self, keys, rule_id, field, old, new, factor, evidence):
        keys = list(keys)
        n = len(keys)
        bc = lambda x: list(x) if isinstance(x, (list, np.ndarray, pd.Series)) else [x] * n
        for k, o, nw, f, e in zip(keys, bc(old), bc(new), bc(factor), bc(evidence)):
            self.rows.append((k, rule_id, field, _fmt(o), _fmt(nw), _fmt(f), e))

    def frame(self):
        return pd.DataFrame(self.rows, columns=LEDGER)


def _match(T, spec, col=None):
    """boolean mask of rows named by spec: genes (list or {gene: reason}) and/or regex patterns, on column spec.match_on
    (default the key); returns (mask, reason per row)"""
    on = spec.get("match_on", "key")
    names = T.index.to_series() if on == "key" else T[on].fillna("").astype(str)
    genes = spec.get("genes") or []
    reasons = genes if isinstance(genes, dict) else {g: spec.get("reason", "") for g in genes}
    m = names.isin(list(reasons)).values
    pats = spec.get("patterns") or []
    if pats:
        m = m | names.str.match("|".join(pats)).values
    if spec.get("regex"):
        m = m | names.str.match(spec["regex"]).values
    return m, names.map(reasons)


def _load_list(ctx, spec):
    """absent-list entry with genes / patterns inline or from a YAML file (spec.file + genes_key / patterns_key)"""
    s = dict(spec)
    if s.get("file"):
        import yaml
        d = yaml.safe_load(open(ctx.resolve(s["file"])))
        if s.get("genes_key"):
            s["genes"] = d[s["genes_key"]]
        if s.get("patterns_key"):
            s["patterns"] = d[s["patterns_key"]]
        if s.get("categories_key"):           # {category: [genes]} -> one entry per category
            return [dict(s, id=f"{s['id']}.{cat}", reason=cat, genes=list(gs), file=None, categories_key=None)
                    for cat, gs in d[s["categories_key"]].items()]
    if s.get("categories"):
        return [dict(s, id=f"{s['id']}.{cat}", reason=cat, genes=list(gs), categories=None) for cat, gs in s["categories"].items()]
    return [s]


def _slopes(ctx, spec):
    """per-protein size-scaling slopes from a table (xlsx or TSV): replicate slopes give a standard error; slopes with
    >= min_replicates and SE <= max_se are kept (clipped), others take the default (median of the slope column).
    origin: yeast_cell_viewer pipeline/inventory.py g1_slopes() (34d123f; Lanz et al. 2024 Table S3)."""
    t = spec["table"]
    p = ctx.resolve(t["path"])
    if p.endswith((".xlsx", ".xls")):
        import warnings
        warnings.filterwarnings("ignore")
        d = pd.read_excel(p, sheet_name=t.get("sheet", 0))
    else:
        d = _read_tsv(p)
    e = d[[c for c in d.columns if str(c).startswith(t["replicate_prefix"])]] if t.get("replicate_prefix") else None
    sl_col = t["slope_column"]
    default = float(d[sl_col].median()) if t.get("default", "median") == "median" else float(t["default"])
    ok = pd.Series(True, index=d.index)
    if e is not None:
        se = e.std(axis=1, ddof=1) / np.sqrt(e.notna().sum(axis=1))
        ok = (e.notna().sum(axis=1) >= int(t.get("min_replicates", 2))) & (se <= float(t.get("max_se", np.inf)))
    sl = d[sl_col].clip(*t["clip"]) if t.get("clip") else d[sl_col]
    ids = d.loc[ok, t["id_column"]].astype(str)
    if t.get("upper", True):
        ids = ids.str.upper()
    return dict(zip(ids, sl[ok])), default, int(ok.sum()), len(d)


def _families(ctx, T, a):
    """anchor families: {family: (stoich, [keys])}"""
    on = a.get("match_on", "key")
    names = T.index.to_series() if on == "key" else T[on].fillna("").astype(str)
    by_name = {}
    for k, nm in zip(T.index, names):
        by_name.setdefault(nm, []).append(k)
    filt = re.compile(a["filter_regex"]) if a.get("filter_regex") else None
    excl = set(a.get("exclude") or [])
    keep = lambda g: g in by_name and g not in excl and (filt is None or filt.match(g))
    fams = {}
    if a.get("families"):
        for f, s in a["families"].items():
            gs = [g for g in (s.get("genes") or []) if keep(g)]
            if s.get("regex"):
                gs += [g for g in by_name if re.match(s["regex"], g) and keep(g) and g not in gs]
            fams[f] = (float(s.get("stoich", 1)), gs)
    if a.get("stoich_map"):
        sm = a["stoich_map"]
        if "file" in sm:                      # {file: <yaml>, key: <mapping key>}
            import yaml
            sm = yaml.safe_load(open(ctx.resolve(sm["file"])))[sm["key"]]
        for g, st in sm.items():
            if keep(g):
                fams[g] = (float(st), [g])
    if a.get("stoich_from_species"):
        s = a["stoich_from_species"]
        P = pd.read_csv(ctx.resolve(s["table"]), sep="\t").set_index(s.get("id_column", "ing_id"))
        r = P.loc[s["species_id"]]
        for g, st in zip(str(r.genes).split(";"), str(r.stoich).split(";")):
            if keep(g):
                fams[g] = (float(st), [g])
    if a.get("members"):
        m = a["members"]
        sub = a.get("family_sub") or [r"$^", ""]
        ov = a.get("stoich_overrides") or {}
        for g in by_name:
            if (re.match(m.get("regex", r"$^"), g) or g in (m.get("genes") or [])) and keep(g):
                f = re.sub(sub[0], sub[1], g)
                fams.setdefault(f, (float(ov.get(f, a.get("stoich", 1))), []))[1].append(g)
    return {f: (st, [k for g in gs for k in by_name[g]]) for f, (st, gs) in fams.items() if gs}


def anchor_table(ctx, T, copies, total_pg, anchors):
    """family ratios copies / (stoich x count) per anchor and count; summary: median ratio and implied total protein"""
    rows, summ = [], {}
    for a in anchors:
        fams = _families(ctx, T, a)
        counts = a.get("counts") or []
        for c in counts:
            rat = {}
            for f, (st, keys) in fams.items():
                cp = float(copies.reindex(keys).fillna(0).sum())
                need = st * float(c["value"])
                rat[f] = cp / need
                rows.append(dict(anchor_id=a["id"], structure=a.get("structure", ""), count_label=c["label"],
                                 count=float(c["value"]), family=f, members=";".join(map(str, keys)), stoich=st, copies=cp,
                                 need=need, ratio=rat[f]))
            v = np.array(list(rat.values()))
            stat = a.get("statistic", "median")
            r = float(np.median(v)) if stat == "median" else float(np.exp(np.mean(np.log(v[v > 0]))))
            summ.setdefault(a["id"], {})[c["label"]] = dict(
                structure=a.get("structure", ""), count=float(c["value"]), count_source=c.get("source", ""),
                n_families=int(len(v)), ratio=r, p10_p90=[float(np.percentile(v, 10)), float(np.percentile(v, 90))] if len(v) else None,
                implied_total_pg=(total_pg / r) if r > 0 else None, primary=bool(c.get("primary", c is counts[0])))
    return pd.DataFrame(rows, columns=ANCHORS), summ


def _combine_totals(vals, rule):
    v = np.array([x for x in vals if x is not None and np.isfinite(x)], dtype=float)
    if not len(v):
        raise ValueError("no anchor gave an implied total")
    if rule == "geometric_mean":
        return float(np.exp(np.mean(np.log(v))))
    if rule == "median":
        return float(np.median(v))
    if rule == "mean":
        return float(np.mean(v))
    if rule in ("min", "max"):
        return float(getattr(np, rule)(v))
    raise ValueError(f"unknown anchor combination rule {rule}")


def _anchor_summary(summ, total_pg, rules=("geometric_mean", "median", "min", "max")):
    prim = {a: next((x for x in d.values() if x["primary"]), next(iter(d.values()))) for a, d in summ.items()}
    imp = {a: x["implied_total_pg"] for a, x in prim.items()}
    vals = [v for v in imp.values() if v]
    out = dict(total_protein_pg_used=total_pg, implied_total_pg_primary=imp)
    if vals:
        out["spread_max_over_min"] = float(max(vals) / min(vals))
        out["combined"] = {r: _combine_totals(vals, r) for r in rules}
    return out


def step(ctx):
    C = _cfg(ctx)
    L = _Ledger()
    AV = float(C.get("avogadro", schema.AVOGADRO))
    summary = dict(method_id=C.get("method_id", ""), avogadro=AV)
    # 1 annotation
    an = C["annotation"]
    ann = reader(an["reader"])(ctx, an)
    # 2 sources (enabled, by priority) + future slots
    srcs, src_rows = [], []
    for s in C.get("sources") or []:
        if s.get("type") not in SOURCE_TYPES:
            raise ValueError(f"copies source {s.get('id')}: type must be one of {sorted(SOURCE_TYPES)}")
        on = s.get("enabled", True)
        src_rows.append(dict(source_id=s["id"], type=s["type"], unit=s.get("unit", ""), reader=s.get("reader", ""),
                             path=str(s.get("path", "")), role=s.get("role", "abundance"), status="enabled" if on else "disabled",
                             column=s.get("column", s["id"]), label=s.get("label", s["id"]), notes=s.get("notes", "")))
        if on:
            v, nrep = reader(s["reader"])(ctx, s, ann)
            srcs.append((s, v, nrep))
    for s in C.get("future_sources") or []:
        if s.get("enabled"):
            raise NotImplementedError(f"copies future source {s['id']} ({s.get('type')}) is declared but not implemented yet")
        src_rows.append(dict(source_id=s["id"], type=s.get("type", ""), unit=s.get("unit", ""), reader=s.get("reader", ""),
                             path=str(s.get("path", "")), role=s.get("role", ""), status="declared, not enabled",
                             column="", label="", notes=s.get("notes", FUTURE_TYPES.get(s.get("type"), ""))))
    if not srcs:
        raise ValueError("inventory.copies: no enabled abundance source")
    imp = C.get("imputation") or {"method": "none"}
    rna = None
    if imp.get("method") == "rna_loglog":
        rna = pd.read_csv(ctx.need(imp.get("rna", "rna/rna.tsv")), sep="\t", index_col=imp.get("rna_key", "gene"))
    # 3 universe
    U = C.get("universe") or {}
    keys = list(ann.index)
    req = U.get("require") or []
    if req:
        have = set()
        if "sources" in req:
            for _, v, _ in srcs:
                have |= set(v.index)
        if "rna" in req and rna is not None:
            have |= set(rna.index)
        keys = [k for k in keys if k in have]
    keys = sorted(set(keys)) if U.get("order", "sorted") == "sorted" else list(dict.fromkeys(keys))
    T = pd.DataFrame(index=pd.Index(keys, name=an.get("key_name", "gene")))
    T["uniprot"] = ann.uniprot.reindex(T.index)
    T["mass_da"] = ann.mass_da.reindex(T.index).astype(float)
    for c in an.get("carry") or []:
        T[c] = ann[c].reindex(T.index)
    for s, v, nrep in srcs:
        col = s.get("column", s["id"])
        T[col] = v.reindex(T.index)
        if s.get("n_rep_column"):
            T[s["n_rep_column"]] = nrep.reindex(T.index).fillna(0).astype(int)
    if rna is not None:
        T[imp.get("tpm_out", "tpm")] = rna[imp.get("tpm_column", "tpm_median")].reindex(T.index)
        T[imp.get("call_out", "rna_call")] = rna[imp.get("call_column", "call")].reindex(T.index).fillna(imp.get("no_rna", "no_rna"))
    # 4 combine
    cmb = C.get("combine") or {}
    order = cmb.get("order") or [s["id"] for s, _, _ in srcs]
    byid = {s["id"]: (s, v, n) for s, v, n in srcs}
    vcol = C.get("value_column", "value")
    prim = byid[order[0]][0]
    vals, harm = {}, {}
    for sid in order:
        s = byid[sid][0]
        x = T[s.get("column", sid)]
        if s["type"] == "lfq":
            x = x / T.mass_da
        vals[sid] = x
    for sid in order[1:]:
        if (byid[sid][0].get("harmonise") or cmb.get("harmonise", "median_ratio")) == "median_ratio":
            both = (vals[order[0]] > 0) & (vals[sid] > 0)
            f = float(np.median(vals[order[0]][both] / vals[sid][both])) if both.any() else 1.0
        else:
            f = 1.0
        harm[sid] = f
        if f != 1.0:
            vals[sid] = vals[sid] * f
    label = pd.Series(None, index=T.index, dtype=object)
    setby = pd.Series("", index=T.index, dtype=object)
    if cmb.get("mode", "priority") == "priority":
        value = pd.Series(np.nan, index=T.index)
        for sid in order:
            s = byid[sid][0]
            m = value.isna() & vals[sid].notna() & (vals[sid] > 0)
            value[m] = vals[sid][m]
            label[m] = s.get("label", sid); setby[m] = sid
            L.add(T.index[m], f"source.{sid}", vcol, "", value[m].values, harm.get(sid, 1.0),
                  f"{s['type']} ({s.get('unit', '')}); {s.get('aggregate', '')}; priority {order.index(sid) + 1}")
    elif cmb["mode"] == "geometric_mean":
        V = pd.DataFrame({sid: vals[sid].where(vals[sid] > 0) for sid in order})
        value = np.exp(np.log(V).mean(axis=1, skipna=True))
        n = V.notna().sum(axis=1)
        m = n > 0
        label[m] = V[m].notna().apply(lambda r: "+".join(V.columns[r.values]), axis=1)
        setby[m] = "combine.geometric_mean"
        L.add(T.index[m], "combine.geometric_mean", vcol, "", value[m].values, "", "geometric mean of harmonised sources")
    else:
        raise ValueError(f"unknown combine mode {cmb['mode']}")
    summary["sources"] = {sid: dict(n_values=int(vals[sid].notna().sum()), n_set=int((setby == sid).sum()),
                                    harmonise_factor=harm.get(sid, 1.0), type=byid[sid][0]["type"]) for sid in order}
    # 5 imputation (or zero for unquantified)
    und = value.isna()
    if imp.get("method") == "rna_loglog":
        fitcol = byid[imp.get("fit_source", order[0])][0].get("column", imp.get("fit_source", order[0]))
        tp, call = T[imp.get("tpm_out", "tpm")], T[imp.get("call_out", "rna_call")]
        d = T[T[fitcol].notna() & (tp > 0)]
        b, a = np.polyfit(np.log10(d[imp.get("tpm_out", "tpm")]), np.log10(d[fitcol]), 1)
        floor = np.percentile(T[fitcol].dropna(), imp.get("floor_percentile", 5))
        fit = 10 ** (a + b * np.log10(tp.clip(lower=imp.get("tpm_clip", 1e-3))))
        lab = {"expressed": "imputed_rna", "low": "imputed_rna_low", "zero": "zero_off"} | (imp.get("labels") or {})
        label[und & (call == "expressed")] = lab["expressed"]
        label[und & (call == "low")] = lab["low"]
        zc = imp.get("zero_calls", ["off", "no_rna"])
        if not all(isinstance(z, str) for z in zc):
            raise ValueError(f"imputation.zero_calls must be strings (YAML reads a bare off/no as false): {zc}")
        label[und & call.isin(zc)] = lab["zero"]
        m = label == lab["expressed"]
        value[m] = np.minimum(fit[m], floor)
        L.add(T.index[m], "impute.rna_loglog", vcol, "", value[m].values, "",
              [f"RNA call expressed, TPM {t:.4g}: min(fit, floor {floor:.6g})" for t in tp[m]])
        m = label == lab["low"]
        value[m] = np.minimum(fit[m], floor) * imp.get("low_factor", 0.25)
        L.add(T.index[m], "impute.rna_loglog_low", vcol, "", value[m].values, imp.get("low_factor", 0.25),
              [f"RNA call low, TPM {t:.4g}: min(fit, floor) x {imp.get('low_factor', 0.25)}" for t in tp[m]])
        m = label == lab["zero"]
        L.add(T.index[m], "impute.rna_off", vcol, "", 0.0, "", [f"undetected, RNA call {c}" for c in call[m]])
        rest = und & label.isna()
        label[rest] = imp.get("missing_label", "not_quantified")
        L.add(T.index[rest], "missing.zero", vcol, "", 0.0, "", [f"undetected, RNA call {c} (no rule)" for c in call[rest]])
        summary["imputation"] = dict(method="rna_loglog", fit_source=fitcol, intercept=float(a), slope=float(b),
                                     floor=float(floor), floor_percentile=imp.get("floor_percentile", 5),
                                     n_fit=int(len(d)), n_imputed=int((label == lab["expressed"]).sum()),
                                     n_imputed_low=int((label == lab["low"]).sum()), n_zero=int((label == lab["zero"]).sum()))
        ctx.say(f"impute fit: log10 {vcol} = {a:.2f} + {b:.2f} log10 TPM; floor {floor:.0f}")
    else:
        label[und] = imp.get("missing_label", "not_quantified")
        L.add(T.index[und], "missing.zero", vcol, "", 0.0, "", "no source value, no imputation")
        summary["imputation"] = dict(method="none", n_zero=int(und.sum()))
    value = value.fillna(0.0)
    rule_of = {"imputed_rna": "impute.rna_loglog", "imputed_rna_low": "impute.rna_loglog_low", "zero_off": "impute.rna_off"}
    if imp.get("method") == "rna_loglog":
        rule_of = {lab[k]: v for k, v in (("expressed", "impute.rna_loglog"), ("low", "impute.rna_loglog_low"),
                                          ("zero", "impute.rna_off"))}
    setby[und] = label[und].map(lambda x: rule_of.get(x, "missing.zero"))
    note = pd.Series("", index=T.index, dtype=object)
    # 6 absent lists
    total_before = float(value.sum())
    summary["absent"] = {}
    absent_any = np.zeros(len(T), bool)
    for spec0 in C.get("absent") or []:
        for s in _load_list(ctx, spec0):
            m, why = _match(T, s)
            share = float(value[m].sum() / total_before) if total_before else 0.0
            L.add(T.index[m], f"absent.{s['id']}", vcol, value[m].values, 0.0, "",
                  [f"{s.get('reason', '')}: {w}" if isinstance(w, str) and w else s.get("reason", "") for w in why[m]])
            value[m] = 0.0
            label[m] = s.get("label", "absent")
            setby[m] = f"absent.{s['id']}"
            if s.get("note", True):
                note[m] = why[m].fillna(s.get("reason", "")).values if s.get("note") in (True, None) else s["note"]
            absent_any |= m
            summary["absent"][s["id"]] = dict(proteins=int(m.sum()), share_of_value=share, reason=s.get("reason", ""))
    if imp.get("method") == "rna_loglog" and imp.get("flag_detected_off", "detected despite RNA off"):
        m = (label == prim.get("label", prim["id"])) & (T[imp.get("call_out", "rna_call")] == "off")
        note[m] = imp.get("flag_detected_off", "detected despite RNA off")
    # 7 per-protein scaling (relative values)
    pfac = pd.Series(1.0, index=T.index)
    summary["per_protein_scaling"] = {}
    for s in C.get("per_protein_scaling") or []:
        on = s.get("match_on", "key")
        ids = T.index.to_series() if on == "key" else T[on].astype(str)
        if s["kind"] == "size_slope":
            sl, default, n_ok, n_tab = _slopes(ctx, s)
            slope = ids.map(sl).fillna(default)
            r = s["ratio"]
            ratio = float(r["numerator"]) / float(r["denominator"]) if isinstance(r, dict) else float(r)
            f = ratio ** slope
            T[s.get("column", "size_slope")] = slope
            summary["per_protein_scaling"][s["rule_id"]] = dict(kind="size_slope", ratio=ratio, default_slope=default,
                                                                n_table=n_tab, n_kept=n_ok, n_matched=int(ids.isin(sl).sum()))
        elif s["kind"] == "factor_table":
            ft = _read_tsv(ctx.resolve(s["path"]))
            f = ids.map(dict(zip(ft[s.get("id_column", "gene")], ft[s.get("factor_column", "factor")]))).fillna(1.0)
            summary["per_protein_scaling"][s["rule_id"]] = dict(kind="factor_table", n_matched=int((f != 1).sum()))
        else:
            raise ValueError(f"unknown per-protein scaling {s['kind']}")
        ch = (f != 1.0) & (value > 0)
        L.add(T.index[ch], f"scale.{s['rule_id']}", vcol, value[ch].values, (value * f)[ch].values, f[ch].values,
              s.get("evidence", s["kind"]))
        value = value * f
        pfac = pfac * f
    # 8 total protein mass
    frac = value / value.sum()
    mean_mass = (frac * T.mass_da.fillna(T.mass_da.median())).sum()
    n_of = lambda pg: pg * 1e-12 / (mean_mass / AV)
    tp = C.get("total_protein") or {"mode": "fixed"}
    anchors = C.get("anchors") or []
    if tp.get("mode", "fixed") == "fixed":
        total = float(tp["value_pg"]) if tp.get("value_pg") is not None else float(ctx.inv["total_protein_pg"])
        how = tp.get("source") or ctx.inv.get("total_protein_source", "")
    elif tp["mode"] == "anchors":
        prov = float(tp.get("provisional_pg", 1.0))
        use = [a for a in anchors if a["id"] in tp["anchors"]]
        _, s0 = anchor_table(ctx, T, frac * n_of(prov), prov, use)
        imp_t = {a: next((x for x in d.values() if x["primary"]), next(iter(d.values())))["implied_total_pg"] for a, d in s0.items()}
        total = _combine_totals(list(imp_t.values()), tp.get("combine", "geometric_mean"))
        if tp.get("round_to"):
            total = round(total / float(tp["round_to"])) * float(tp["round_to"])
        how = f"anchors {sorted(imp_t)} -> {tp.get('combine', 'geometric_mean')} = {total:.6g} pg ({imp_t})"
    elif tp["mode"] == "source_sum":
        sid = tp.get("source_id", order[0])
        col = byid[sid][0].get("column", sid)
        total = float((T[col].fillna(0) * T.mass_da).sum() / AV * 1e12)
        how = f"sum(copies x mass) of absolute source {sid}"
    else:
        raise ValueError(f"unknown total_protein.mode {tp['mode']}")
    N = n_of(total)
    cf = frac * N
    L.add(["*"], f"total.{tp.get('mode', 'fixed')}", "copies", "", total, N, f"total protein {total:.6g} pg ({how}); "
          f"mean mass {mean_mass:.6g} Da; {N:.6g} molecules")
    summary["total_protein"] = dict(mode=tp.get("mode", "fixed"), total_protein_pg=total, source=how, molecules=float(N),
                                    mean_mass_da=float(mean_mass))
    ctx.say(f"total {total} pg -> {N:.3g} molecules; mean mass {mean_mass / 1e3:.1f} kDa ({how})")
    # 9 adjustments (structure counts)
    afac = pd.Series(1.0, index=T.index)
    summary["adjustments"] = {}
    for s in C.get("adjustments") or []:
        if s["kind"] != "family_target":
            raise ValueError(f"unknown adjustment {s['kind']}")
        on = s.get("match_on", "key")
        names = T.index.to_series() if on == "key" else T[on].fillna("").astype(str)
        mem = s["members"]
        m = names.str.match(mem.get("regex", r"$^")) | names.isin(mem.get("genes") or [])
        sub = s.get("family_sub") or [r"$^", ""]
        fam = cf[m].groupby(names[m].str.replace(sub[0], sub[1], regex=True)).sum()
        tg = s["target"]
        target = float(tg["count"]) if "count" in tg else float(tg["density_per_um3"]) * float(tg["volume_um3"])
        before = float(np.median(fam[fam > 0]))
        f = target / before
        old = cf[m].values
        cf = cf.copy()
        cf[m] *= f
        afac[m] *= f
        L.add(T.index[m], f"adjust.{s['rule_id']}", "copies", old, cf[m].values, f,
              f"{s.get('structure', '')}: median family {before:.6g} -> target {target:.6g} ({tg.get('source', '')})")
        pg = float((cf * T.mass_da.fillna(T.mass_da.median())).sum() / AV * 1e12)
        summary["adjustments"][s["rule_id"]] = dict(factor=f, target=target, before=before, n_proteins=int(m.sum()),
                                                    n_families=int(len(fam)), protein_pg_after=pg)
        ctx.say(f"adjustment {s['rule_id']}: {s.get('structure', '')} {target:,.0f} (median family {before:,.0f}; x {f:.3f}); "
                f"protein total {pg:.3f} pg")
    # 10 rounding
    copies = cf.round().astype(np.int64) if C.get("round_copies", True) else cf
    # 11 anchor check (diagnostic)
    ac = C.get("anchor_check") or {}
    A, asum = anchor_table(ctx, T, copies.astype(float) if ac.get("use", "final") == "final" else frac * N, total, anchors)
    summary["anchors"] = asum
    summary["anchor_summary"] = _anchor_summary(asum, total) if asum else {}
    for a, d in asum.items():
        for lab_, x in d.items():
            ctx.say(f"anchor {a} [{lab_}: {x['count']:,.0f}]: median ratio {x['ratio']:.3f} over {x['n_families']} families "
                    f"-> implied total {x['implied_total_pg'] if x['implied_total_pg'] is None else round(x['implied_total_pg'], 3)} pg")
    # 12 table + ledger
    T["source"] = label
    T[vcol] = value
    T["note"] = note
    if C.get("per_protein_scaling"):
        T["scale_factor"] = pfac
    if C.get("adjustments"):
        T["adjust_factor"] = afac
    T["copies"] = copies
    T["set_by"] = setby
    Lf = L.frame()
    rid = Lf[Lf.key != "*"].groupby("key").rule_id.agg(";".join)
    T["rule_ids"] = rid.reindex(T.index).fillna("").values
    summary["labels"] = {k: int(v) for k, v in label.value_counts().items()}
    summary["n_proteins"] = int(len(T))
    summary["copies_total"] = float(copies.sum())
    summary["protein_pg"] = float((copies * T.mass_da.fillna(T.mass_da.median())).sum() / AV * 1e12)
    summary["future_sources"] = [r["source_id"] for r in src_rows if r["status"].startswith("declared")]
    ctx.write_tsv(T, "copies/copies.tsv", index=True)
    ctx.write_tsv(Lf, "copies/copy_ledger.tsv")
    ctx.write_tsv(pd.DataFrame(src_rows, columns=["source_id", "type", "unit", "reader", "path", "role", "status", "column",
                                                  "label", "notes"]), "copies/copy_sources.tsv")
    ctx.write_tsv(A, "copies/anchors.tsv")
    from .. import build
    p = ctx.out("copies/summary.json")
    build.atomic_save(p, lambda t: json.dump(summary, open(t, "w"), indent=1, default=float))
    ctx.write_report("copies/report.txt")
    return T, summary


def frozen_names(cfg):
    """frozen inputs named anywhere in inventory.copies"""
    names = set((cfg.get("frozen_inputs") or {}).keys())
    out = set()

    def scan(o):
        if isinstance(o, str):
            if o in names:
                out.add(o)
        elif isinstance(o, dict):
            for v in o.values():
                scan(v)
        elif isinstance(o, list):
            for v in o:
                scan(v)
    scan((cfg.get("inventory") or {}).get("copies") or {})
    return sorted(out)


def needs_rna(cfg):
    c = (cfg.get("inventory") or {}).get("copies") or {}
    return (c.get("imputation") or {}).get("method") == "rna_loglog"


# ------------------------------------------------------------------------------------------- proteome adapter `copies`
def copies_adapter(ctx, total_protein_pg=None):
    """adapter.proteome: the copies step's table (inventory/copies/copies.tsv), columns inventory.copies.v0_columns
    (default: every column except the ledger helpers); floats read back exactly (round_trip)"""
    C = _cfg(ctx)
    T = pd.read_csv(ctx.need("copies/copies.tsv"), sep="\t", index_col=0, keep_default_na=False, na_values=[""],
                    float_precision="round_trip")
    cols = C.get("v0_columns") or [c for c in T.columns if c not in ("set_by", "rule_ids")]
    out = T[cols].copy()
    if "note" in out.columns:
        out["note"] = out.note.fillna("")
    out.index.name = "gene"
    return out


copies_adapter.PREV = ["copies/copies.tsv"]
register("adapter.proteome", "copies")(copies_adapter)
