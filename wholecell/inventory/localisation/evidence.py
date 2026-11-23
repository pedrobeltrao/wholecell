"""Localisation evidence: the long evidence table and the source adapters that write it.

evidence.tsv (inventory/localisation/evidence.tsv; one row per gene x source x compartment term, in source order):
  gene, uniprot      the protein (gene = the inventory key; uniprot = accession when known)
  source             short source label (config `sources[].source`, e.g. hpa); also the label in loc_source
  evidence_class     one of CLASSES (precedence / blend classes of the combiner)
  compartment        a localisation key of the cell (mapped by the adapter's term map)
  score              0-1 (annotation term weight, imaging score, or the observed nuclear fraction for fractionation)
  reliability        0-1 or empty (HPA Enhanced 1 / Supported 0.8 / Approved 0.5 / Uncertain 0)
  n_obs              observations behind the row (fractionation profiles; 1 for a term)
  note               the source term / class that gave the row

gene_features.tsv (inventory/localisation/gene_features.tsv): gene, uniprot, transmembrane (UniProt Transmembrane
feature present) - the per-gene features the combiner rules read.

Adapters (SOURCES, a dict name -> function f(genes, spec, ctx) -> list of rows): hpa_subcellular (HPA
subcellular_location.tsv), uniprot_cc (UniProt 'Subcellular location [CC]' text, regex rules), fractionation_nc (Itzhak
2016-style N/C partition table), localisation_table (a TSV gene, compartments 'key:w;key:w', as used by the testcell);
keyed by ordered locus (yeast legacy_x5, rows in source order = rank): uniprot_cc_statements, go_cc_terms,
class_columns, proteome_classes, imaging_locscore (genes.reader uniprot_ordered_locus -> locus_features).
Term maps are cell data (TSVs with a reason column); adapters contain no compartment lists.
"""
from __future__ import annotations

import io
import os
import re

import pandas as pd

COLUMNS = ["gene", "uniprot", "source", "evidence_class", "compartment", "score", "reliability", "n_obs", "note"]
FEATURES = ["gene", "uniprot", "transmembrane"]
CLASSES = ("imaging_cell", "imaging_other", "annotation_curated", "annotation_ht", "fractionation", "proteomics_surface")
HPA_RELIABILITY = {"Enhanced": 1.0, "Supported": 0.8, "Approved": 0.5, "Uncertain": 0.0}


def read_curated(path: str) -> pd.DataFrame:
    """curated TSV (repo-relative or absolute; '#' lines are provenance comments), all columns as strings"""
    from ...config import REPO
    p = path if os.path.isabs(path) else os.path.join(REPO, path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"curated localisation table {path} not found (looked at {p})")
    txt = "".join(line for line in open(p) if not line.startswith("#"))
    return pd.read_csv(io.StringIO(txt), sep="\t", dtype=str, keep_default_na=False)


def uniprot_genes(path: str) -> pd.DataFrame:
    """UniProt dump -> reviewed entries indexed by primary gene name (first of 'Gene Names (primary)'; first entry per
    gene). The gene -> accession map and features of the Jurkat inventory (legacy inventory_v0._uniprot)."""
    up = pd.read_csv(path, sep="\t")
    up = up[up.Reviewed == "reviewed"].rename(columns={"Gene Names (primary)": "gene"})
    up = up.dropna(subset=["gene"]).drop_duplicates("gene")
    up["gene"] = up.gene.str.split("; ").str[0]
    return up.drop_duplicates("gene").set_index("gene")


def gene_features(genes, up: pd.DataFrame) -> pd.DataFrame:
    F = pd.DataFrame(index=pd.Index(genes, name="gene"))
    F["uniprot"] = up.Entry.reindex(F.index)
    F["transmembrane"] = up.Transmembrane.reindex(F.index).notna().values
    return F


def _row(g, acc, spec, comp, score, rel, n, note):
    return (g, acc, spec["source"], spec["evidence_class"], comp, score, rel, n, note)


# ----------------------------------------------------------------------------------------------------------- adapters
def hpa_subcellular(genes, spec, ctx):
    """HPA subcellular_location.tsv: Main location terms score main_score, Additional location terms additional_score;
    term -> key from data term_map (TSV term, compartment, reason); reliability from the HPA Reliability column."""
    tm = read_curated(spec["term_map"])
    term2comp = dict(zip(tm.term, tm.compartment))
    rel_map = dict(HPA_RELIABILITY, **(spec.get("reliability") or {}))
    hpa = pd.read_csv(ctx.resolve(spec["input"]), sep="\t").drop_duplicates("Gene name").set_index("Gene name")
    cols = (("Main location", float(spec.get("main_score", 1.0))), ("Additional location", float(spec.get("additional_score", 0.3))))
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        if g not in hpa.index:
            continue
        rl = hpa.at[g, "Reliability"]
        rel = rel_map.get(rl, float("nan"))
        for col, sc in cols:
            for t in str(hpa.at[g, col]).split(";"):
                c = term2comp.get(t.strip())
                if c:
                    rows.append(_row(g, acc.get(g), spec, c, sc, rel, 1, f"{col.split()[0].lower()}:{t.strip()} ({rl})"))
    return rows


def _uniprot_cc_text(s: pd.Series) -> pd.Series:
    """evidence tags {...} and the 'SUBCELLULAR LOCATION:' heads removed, notes (Note=...) cut"""
    return (s.str.replace(r"\{[^}]*\}", "", regex=True).str.replace("SUBCELLULAR LOCATION:", "", regex=False)
            .str.split("Note=").str[0])


def uniprot_cc(genes, spec, ctx):
    """UniProt 'Subcellular location [CC]' text of the gene's reviewed entry: every rule (TSV compartment, regex, reason;
    searched in the lower-cased text, in file order) that matches gives a row with score 1; keys in drop_if_other are
    dropped when another key matched (e.g. 'secreted' next to a cellular location)."""
    rules = read_curated(spec["rules"])
    comp_rx = list(zip(rules.compartment, rules.regex))
    drop_if_other = set(spec.get("drop_if_other") or [])
    txt = _uniprot_cc_text(ctx.uniprot_table["Subcellular location [CC]"].reindex(genes).fillna(""))
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        t = txt[g].lower()
        hit = [c for c, pat in comp_rx if re.search(pat, t)]
        if hit and set(hit) - drop_if_other:
            hit = [c for c in hit if c not in drop_if_other]
        for c in hit:
            rows.append(_row(g, acc.get(g), spec, c, 1.0, float("nan"), 1, "cc"))
    return rows


def fractionation_nc(genes, spec, ctx):
    """N/C partition table (Itzhak et al. 2016 layout: gene, frac_nuclear, frac_cytosol, n_good_mnc_profiles, mnc_class,
    primary_for_gene): one row per gene, compartment nucleus, score = observed nuclear fraction N / (N + C); kept when
    n_good_mnc_profiles >= min_profiles and mnc_class in classes. The background correction is the combiner's nc_model."""
    Z = pd.read_csv(ctx.resolve(spec["input"]), sep="\t")
    Z = Z[Z.primary_for_gene & (Z.n_good_mnc_profiles >= spec.get("min_profiles", 3))
          & Z.mnc_class.isin(spec.get("classes", ["N", "C", "N/C"]))].drop_duplicates("gene").set_index("gene")
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        if g not in Z.index:
            continue
        nc, cy = Z.at[g, "frac_nuclear"], Z.at[g, "frac_cytosol"]
        rows.append(_row(g, acc.get(g), spec, "nucleus", nc / (nc + cy), float("nan"), int(Z.at[g, "n_good_mnc_profiles"]),
                         f"mnc_class {Z.at[g, 'mnc_class']}"))
    return rows


def localisation_table(genes, spec, ctx):
    """TSV gene, compartments ('key:w;key:w'): one row per key, score = w (generic adapter, e.g. the testcell table)"""
    T = pd.read_csv(ctx.resolve(spec["input"]), sep="\t", comment="#").drop_duplicates("gene").set_index("gene")
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        if g not in T.index or not isinstance(T.at[g, "compartments"], str):
            continue
        for kv in T.at[g, "compartments"].split(";"):
            k, v = kv.rsplit(":", 1)
            rows.append(_row(g, acc.get(g), spec, k, float(v), float("nan"), 1, "table"))
    return rows


# ------------------------------------------------- ordered-locus adapters (yeast legacy_x5; usable by any UniProt dump)
# Keyed by the first ordered-locus name of the UniProt dump (yeast: the ORF), one entry per key (first row), as the
# yeast inventory. Every row keeps the order of the source (statement / term order = the rank of the X5 split).
_LOCUS_CACHE: dict = {}


def locus_table(path: str, key_column: str = "Gene Names (ordered locus)") -> pd.DataFrame:
    """UniProt dump -> one row per ordered-locus key (first token, split at blank or ';'; first row per key), file order
    # origin: yeast pipeline/inventory.py:build (u.orf ... drop_duplicates('orf')) @ e3b7448"""
    k = (path, key_column)
    if k not in _LOCUS_CACHE:
        u = pd.read_csv(path, sep="\t")
        u["_key"] = u[key_column].astype(str).str.split(r"[ ;]").str[0]
        _LOCUS_CACHE[k] = (u, u.drop_duplicates("_key").set_index("_key"))
    return _LOCUS_CACHE[k][1]


def _locus_raw(path: str, key_column: str = "Gene Names (ordered locus)") -> pd.DataFrame:
    locus_table(path, key_column)
    return _LOCUS_CACHE[(path, key_column)][0]


def locus_features(path: str, key_column: str = "Gene Names (ordered locus)", tm_column: str = "Transmembrane",
                   cc_column: str = "Subcellular location [CC]", barrel_text: str = "Multi-pass membrane protein",
                   genes=None) -> pd.DataFrame:
    """gene_features keyed by ordered locus: uniprot, gene_name, n_tm (TRANSMEM ranges; a multi-pass protein without
    ranges = an outer-membrane beta-barrel counts 1), transmembrane = n_tm > 0
    # origin: yeast pipeline/inventory.py:build (n_tm, barrel) @ e3b7448"""
    u = locus_table(path, key_column)
    F = pd.DataFrame(index=pd.Index(list(u.index) if genes is None else list(genes), name="gene"))
    F["uniprot"] = u.Entry.reindex(F.index).values
    F["gene_name"] = u["Gene Names (primary)"].reindex(F.index).values
    ntm = u[tm_column].fillna("").str.count("TRANSMEM")
    barrel = u[cc_column].fillna("").str.contains(barrel_text, regex=False) & (ntm == 0)
    ntm = ntm.where(~barrel, 1)
    F["n_tm"] = ntm.reindex(F.index).fillna(0).astype(int).values
    F["transmembrane"] = F.n_tm > 0
    return F


def _first_match(rules, t):
    for c, pat in rules:
        if re.search(pat, t):
            return c
    return None


def _refine(keys, extra, generic, subs):
    """generic key -> the sub-keys named in extra (order of extra), deduplicated
    # origin: yeast pipeline/inventory.py:refine_mito @ e3b7448"""
    if generic not in keys:
        return keys
    sub = [c for c in dict.fromkeys(extra) if c in subs]
    if not sub:
        return keys
    out = []
    for c in keys:
        out += sub if c == generic else [c]
    return list(dict.fromkeys(out))


def uniprot_cc_statements(genes, spec, ctx):
    """UniProt 'Subcellular location [CC]' (before 'Note='), statement by statement in UniProt order: each statement ->
    the first rule (TSV compartment, regex, reason) matching its lower-cased text; rows in statement order, score 1.
      drop_only_pmids    statements whose PubMed evidence is only these papers (and no ECO:0000250) are dropped
                         (yeast: the Sickmann 2003 / Reinders 2006 organelle proteomes)
      precedence         {pmid, dedicated_eco, ht_min_entries}: when any statement cites `pmid` (Huh 2003 GFP), only the
                         statements citing it or with `dedicated_eco` evidence from a paper cited in < ht_min_entries
                         entries of the dump are kept (note 'precedence:<pmid>')
      refine             {generic, sub_keys}: the generic key is replaced by the sub-keys of ALL statements
    # origin: yeast pipeline/inventory.py:parse_uniprot_cc, _pmid_entries @ e3b7448"""
    import collections
    rules = read_curated(spec["rules"])
    comp_rx = list(zip(rules.compartment, rules.regex))
    path = ctx.resolve(spec["input"])
    kc = spec.get("key_column", "Gene Names (ordered locus)")
    col = spec.get("column", "Subcellular location [CC]")
    u = locus_table(path, kc)
    drop_only = {str(p) for p in spec.get("drop_only_pmids") or []}
    keep_eco = spec.get("keep_eco", "ECO:0000250")
    pr = spec.get("precedence") or {}
    ref = spec.get("refine") or {}
    cnt = collections.Counter()
    if pr:
        for s in _locus_raw(path, kc)[col].dropna():
            cnt.update(set(re.findall(r"PubMed:(\d+)", s.split("Note=")[0])))
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        s = u[col].get(g) if g in u.index else None
        if not isinstance(s, str):
            continue
        s = s.split("Note=")[0].replace("SUBCELLULAR LOCATION:", "")
        out, kept, has = [], [], False
        for m in re.finditer(r"([^.;{}]+?)\s*(\{[^}]*\})?\s*[.;]", s):
            t, ev = m.group(1).strip().lower(), (m.group(2) or "")
            if not t:
                continue
            pm = set(re.findall(r"PubMed:(\d+)", ev))
            if pm and pm <= drop_only and keep_eco not in ev:
                continue
            c = _first_match(comp_rx, t)
            if c is None:
                continue
            out.append(c)
            if pr:
                has |= str(pr["pmid"]) in pm
                if str(pr["pmid"]) in pm or (pr["dedicated_eco"] in ev and any(cnt[p] < pr["ht_min_entries"] for p in pm)):
                    kept.append(c)
        keys = list(dict.fromkeys(kept if has else out))
        if ref:
            keys = _refine(keys, out, ref["generic"], set(ref["sub_keys"]))
        note = f"precedence:{pr['pmid']}" if has else "cc"
        for c in keys:
            rows.append(_row(g, acc.get(g), spec, c, 1.0, float("nan"), 1, note))
    return rows


def go_cc_terms(genes, spec, ctx):
    """GO cellular component terms ('term [GO:id]; ...', lower-cased) in the dump's order: each term -> the first rule
    (extra_rules first, then the rules TSV); rows in term order, score 1 (GO terms come sorted alphabetically: the order
    carries no rank, see the combiner's ordered_annotation.fallback_unranked).
      skip               {regex, unless}: terms matching regex and not `unless` are skipped (complexes)
      exclude_keys       matched terms of these keys give no row (the term still stops at its first rule)
      restrict_if        {key, keep}: when `key` is present, only the `keep` keys stay
    # origin: yeast pipeline/inventory.py:parse_go_cc, tm_go_fallback @ e3b7448"""
    rules = read_curated(spec["rules"])
    comp_rx = [(r["compartment"], r["regex"]) for r in spec.get("extra_rules") or []] + list(zip(rules.compartment, rules.regex))
    u = locus_table(ctx.resolve(spec["input"]), spec.get("key_column", "Gene Names (ordered locus)"))
    col = spec.get("column", "Gene Ontology (cellular component)")
    skip = spec.get("skip") or {}
    excl = set(spec.get("exclude_keys") or [])
    rif = spec.get("restrict_if") or {}
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        s = u[col].get(g) if g in u.index else None
        if not isinstance(s, str):
            continue
        out = []
        for t in s.lower().split(";"):
            t = re.sub(r"\[go:\d+\]", "", t).strip()
            if skip and re.search(skip["regex"], t) and not re.search(skip["unless"], t):
                continue
            c = _first_match(comp_rx, t)
            if c is not None and c not in excl:
                out.append(c)
        out = list(dict.fromkeys(out))
        if rif and rif["key"] in out:
            out = [c for c in out if c in rif["keep"]]
        for c in out:
            rows.append(_row(g, acc.get(g), spec, c, 1.0, float("nan"), 1, "go"))
    return rows


def _read_any(path, **kw):
    """Excel (.xls/.xlsx; pandas.read_excel) or a lossless TSV conversion of it (.tsv: the table below the header row)"""
    if path.endswith((".tsv", ".txt")):
        return pd.read_csv(path, sep="\t", float_precision="round_trip")
    kw = {k: v for k, v in kw.items() if v is not None}
    return pd.read_excel(path, **kw)


def class_columns(genes, spec, ctx):
    """a table with one column per class (a non-empty text cell = assigned; Vogtle 2017 sub-mitochondrial classes):
    genes with exactly one assigned class (exactly_one) -> one row, compartment = classes[column], score 1. Later rows of
    the same id replace earlier ones.
    # origin: yeast pipeline/inventory.py:vogtle_submito @ e3b7448"""
    v = _read_any(ctx.resolve(spec["input"]), header=spec.get("header", 0), sheet_name=spec.get("sheet"))
    cls = spec["classes"]
    out = {}
    for _, r in v.iterrows():
        subs = [(k, c) for k, c in cls.items() if isinstance(r[k], str)]
        if len(subs) == 1 or (subs and not spec.get("exactly_one", True)):
            out[r[spec.get("id_column", "ORF")]] = subs
    acc = ctx.uniprot_of
    return [_row(g, acc.get(g), spec, c, 1.0, float("nan"), 1, f"class:{k}") for g in genes if g in out for k, c in out[g]]


def proteome_classes(genes, spec, ctx):
    """an organelle proteome table (Morgenstern 2017: pure vs crude mitochondria): one row per gene in the table,
    compartment = sub_classes[sub-class column] (the organelle key when none / unmapped), score 1 = in the high-confidence
    proteome (hc_column == 1) else 0, reliability 0 = depleted from the pure organelle (class_column == depleted_class)
    else 1; note '<class>; <sub-class>'. The combiner's mito_fractionation reads these.
    # origin: yeast pipeline/inventory.py:morgenstern @ e3b7448"""
    m = _read_any(ctx.resolve(spec["input"]), header=spec.get("header", 0), sheet_name=spec.get("sheet"))
    sub_map = spec["sub_classes"]
    out = {}
    for _, r in m.iterrows():
        cls = r[spec["class_column"]] if isinstance(r[spec["class_column"]], str) else None
        sub = r[spec["sub_column"]] if isinstance(r[spec["sub_column"]], str) else None
        out[r[spec["id_column"]]] = (r[spec["hc_column"]] == 1, cls, sub)
    acc = ctx.uniprot_of
    rows = []
    for g in genes:
        if g not in out:
            continue
        hc, cls, sub = out[g]
        rows.append(_row(g, acc.get(g), spec, sub_map.get(sub, spec["organelle_key"]), 1.0 if hc else 0.0,
                         0.0 if cls == spec["depleted_class"] else 1.0, 1, f"{cls}; {sub}"))
    return rows


def _prenylated(spec, ctx):
    """ids with a prenylated cysteine (S-farnesyl / S-geranylgeranyl) within the last `last_n` residues (UniProt
    Lipidation features; lengths from the main dump)
    # origin: yeast pipeline/inventory.py:prenylated @ e3b7448"""
    u = pd.read_csv(ctx.resolve(spec["features"]), sep="\t", dtype=str)
    ln = pd.read_csv(ctx.resolve(spec["lengths"]), sep="\t", usecols=[0, 4])
    ln = dict(zip(ln.iloc[:, 0], ln.iloc[:, 1]))
    out = set()
    for a, orf, lip in zip(u["Entry"], u["Gene Names (ordered locus)"], u["Lipidation"]):
        if isinstance(orf, str) and isinstance(lip, str) and a in ln:
            pos = [int(p) for p in re.findall(r"LIPID (\d+); /note=\"S-(?:farnesyl|geranylgeranyl)", lip)]
            if any(p >= ln[a] - spec.get("last_n", 4) for p in pos):
                out.add(orf.split()[0].upper())
    return out


def imaging_locscore(genes, spec, ctx):
    """imaging LOC-scores (CYCLoPs, Koh et al. 2015: fraction of cells showing the protein in a class), mean over the
    screens (inputs) that imaged the protein; classes with mean >= min_score give one row each, in column order:
    compartment = term_map[class] (TSV term, compartment, reason; may be a pseudo-key such as 'periphery' that the
    combiner resolves), score = mean LOC-score, note = the class. Exclusions (source not used for the gene):
    exclude.gene_regex on the primary gene name (GFP-tagged r-proteins), exclude.prenylated (C-terminal GFP blocks
    prenylation: {features, lengths, last_n}).
    # origin: yeast pipeline/inventory.py:cyclops_scores, build (RPROT, PRENYL_NO_CYCLOPS) @ e3b7448"""
    tm = read_curated(spec["term_map"])
    term2comp = dict(zip(tm.term, tm.compartment))
    idc = spec.get("id_column", "ORF")
    d = pd.concat([_read_any(ctx.resolve(f), header=spec.get("header", 0)).rename(columns={idc: "_id"})
                   for f in spec["inputs"]])
    cls = [c for c in d.columns if c in term2comp]
    m = d.groupby("_id")[cls].mean()
    mn = float(spec.get("min_score", 0.05))
    ex = spec.get("exclude") or {}
    skip = set()
    if ex.get("gene_regex") or ex.get("prenylated"):
        u = locus_table(ctx.resolve(ex.get("genes_input", spec.get("genes_input"))))
        if ex.get("gene_regex"):
            rx = re.compile(ex["gene_regex"])
            gn = u["Gene Names (primary)"]
            skip |= {g for g in genes if g in u.index and rx.match(gn[g] if isinstance(gn[g], str) else "")}
        if ex.get("prenylated"):
            skip |= _prenylated(ex["prenylated"], ctx) & set(genes)
    acc = ctx.uniprot_of
    rows = []
    idx = set(m.index)
    for g in genes:
        if g not in idx or g in skip:
            continue
        for c, v in m.loc[g].items():
            if v >= mn:
                rows.append(_row(g, acc.get(g), spec, term2comp[c], float(v), float("nan"), 1, c))
    return rows


SOURCES = {"hpa_subcellular": hpa_subcellular, "uniprot_cc": uniprot_cc, "fractionation_nc": fractionation_nc,
           "localisation_table": localisation_table, "uniprot_cc_statements": uniprot_cc_statements,
           "go_cc_terms": go_cc_terms, "class_columns": class_columns, "proteome_classes": proteome_classes,
           "imaging_locscore": imaging_locscore}


def to_frame(rows) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=COLUMNS)


def read(path: str) -> pd.DataFrame:
    """evidence.tsv -> DataFrame with exact float round trip (scores feed exact arithmetic)"""
    return pd.read_csv(path, sep="\t", keep_default_na=False, na_values={"score": [""], "reliability": [""], "n_obs": [""]},
                       dtype={"gene": str, "uniprot": str, "source": str, "evidence_class": str, "compartment": str,
                              "note": str}, float_precision="round_trip")
