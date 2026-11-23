"""SEC-MS input adapters -> one long peptide table or one long protein table.

Long peptide table (`PEP`): protein, peptide, condition, replicate, fraction, intensity
Long protein table (`PROT`): protein, condition, replicate, fraction, intensity [, sem]

Formats (dataset key `format`):
  openswath       OpenSWATH/TRIC feature_alignment (tab separated): filename, FullPeptideName, ProteinName, Intensity,
                  m_score, decoy; plus a fraction map (injection file -> fraction [, condition, replicate])
  peptide_long    any TSV with columns named in `columns` {protein, peptide, fraction, intensity [, condition,
                  replicate, qvalue]}
  protein_matrix  wide protein x fraction TSV; value columns matched by `value_regex` with named groups `fraction`
                  [, `condition`, `replicate`]; optional `sem_regex` (same groups) for a standard error per value
  protein_long    TSV with columns {protein, fraction, intensity [, condition, replicate, sem]}
Protein-level dataset options (protein_matrix, protein_long; dataset properties, not parameters):
  duplicate_rows  (protein_long) sum (default: rows of one protein x fraction are summed) | unique (exports at precursor
                  level that repeat the protein-group quantity on every precursor row, e.g. Spectronaut PG.IBAQ:
                  the repeated value is taken once; differing values for one protein x fraction raise)
  exclude_ids_regex  protein ids to drop (e.g. '^NA$|^CON_' for unmapped groups and contaminants)
  id_groups       {strip: '[]', sep: ';', mode: expand | first | drop}: protein-group ids such as '[P1;P2]' (paralog
                  groups); expand = the group's profile for every member (shape-based outputs are unchanged; summed
                  abundances then count the group once per member), first = first accession, drop = remove the group
Missing values (NA / blank) are 'not detected' = 0 intensity, as in every reader.
No dataset-specific constant lives here: column names, regexes, skipped header lines come from the dataset config.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

PEP = ["protein", "peptide", "condition", "replicate", "fraction", "intensity"]
PROT = ["protein", "condition", "replicate", "fraction", "intensity"]


def read_fraction_map(spec: dict, resolve) -> pd.DataFrame:
    """fraction map -> DataFrame file, fraction, condition, replicate.
    spec: {path, skiprows (int, default 0), file_col, fraction_col (names or 0-based positions), condition, replicate}.
    Rows whose fraction is not an integer are dropped (sample sheets often mix library / result rows)."""
    path = resolve(spec["path"])
    T = pd.read_csv(path, sep=spec.get("sep", "\t"), skiprows=int(spec.get("skiprows", 0))).dropna()
    fc, rc = spec.get("file_col", 0), spec.get("fraction_col", 1)
    fcol = T.columns[fc] if isinstance(fc, int) else fc
    rcol = T.columns[rc] if isinstance(rc, int) else rc
    M = pd.DataFrame({"file": T[fcol].astype(str), "fraction": T[rcol]})
    M = M[M.fraction.astype(str).str.isdigit()].copy()
    M["fraction"] = M.fraction.astype(int)
    for k in ("condition", "replicate"):
        col = spec.get(f"{k}_col")
        M[k] = T.loc[M.index, col].values if col else spec.get(k, "default" if k == "condition" else 1)
    return M.reset_index(drop=True)


def read_openswath(ds: dict, q: dict, resolve) -> pd.DataFrame:
    fm = read_fraction_map(ds["fraction_map"], resolve)
    cols = ["decoy", "filename", "FullPeptideName", "Intensity", "m_score", "ProteinName"]
    fa = pd.read_csv(resolve(ds["path"]), sep="\t", usecols=cols)
    fa = fa[(fa.decoy == 0) & (fa.m_score < float(q["qvalue_max"]))]
    if q.get("exclude_protein_regex"):
        fa = fa[~fa.ProteinName.str.contains(q["exclude_protein_regex"])]
    fa = fa.assign(file=fa.filename.str.split("/").str[-1])
    fa = fa.merge(fm, on="file")
    return pd.DataFrame({"protein": fa.ProteinName.values, "peptide": fa.FullPeptideName.values,
                         "condition": fa.condition.values, "replicate": fa.replicate.values,
                         "fraction": fa.fraction.values, "intensity": fa.Intensity.values})


def read_peptide_long(ds: dict, q: dict, resolve) -> pd.DataFrame:
    c = {"protein": "protein", "peptide": "peptide", "fraction": "fraction", "intensity": "intensity"} | (ds.get("columns") or {})
    T = pd.read_csv(resolve(ds["path"]), sep="\t", comment="#")
    if c.get("qvalue") and c["qvalue"] in T.columns:
        T = T[T[c["qvalue"]] < float(q["qvalue_max"])]
    if q.get("exclude_protein_regex"):
        T = T[~T[c["protein"]].astype(str).str.contains(q["exclude_protein_regex"])]
    out = pd.DataFrame({k: T[c[k]].values for k in ("protein", "peptide", "fraction", "intensity")})
    out["condition"] = T[c["condition"]].values if c.get("condition") in T.columns else ds.get("condition", "default")
    out["replicate"] = T[c["replicate"]].values if c.get("replicate") in T.columns else 1
    out["fraction"] = out.fraction.astype(int)
    return out[PEP]


def clean_ids(out: pd.DataFrame, ds: dict) -> pd.DataFrame:
    """dataset options exclude_ids_regex and id_groups on a protein-level table (column `protein`)"""
    out = out[out.protein.notna()].copy()
    out["protein"] = out.protein.astype(str)
    if ds.get("exclude_ids_regex"):
        out = out[~out.protein.str.contains(ds["exclude_ids_regex"])]
    g = ds.get("id_groups")
    if g:
        strip, sep, mode = g.get("strip", "[]"), g.get("sep", ";"), g.get("mode", "expand")
        if mode not in ("expand", "first", "drop"):
            raise ValueError(f"id_groups.mode must be expand | first | drop (got {mode!r})")
        ids = out.protein.str.strip(strip)
        multi = ids.str.contains(sep, regex=False)
        if mode == "drop":
            out, ids = out[~multi], ids[~multi]
            out = out.assign(protein=ids.values)
        elif mode == "first":
            out = out.assign(protein=ids.str.split(sep, regex=False).str[0].values)
        else:
            out = out.assign(protein=ids.str.split(sep, regex=False)).explode("protein")
        # a member that is also quantified on its own keeps its own profile
        if mode != "drop" and multi.any():
            own = set(out.protein[~multi.reindex(out.index).values])
            out = out[~(multi.reindex(out.index).values & out.protein.isin(own).values)]
    return out.reset_index(drop=True)


def read_protein_matrix(ds: dict, resolve) -> pd.DataFrame:
    T = pd.read_csv(resolve(ds["path"]), sep="\t", comment="#")
    idc = ds.get("id_col", "protein")
    rx = re.compile(ds["value_regex"])
    sem_of = {}
    if ds.get("sem_regex"):
        srx = re.compile(ds["sem_regex"])
        for c in T.columns:
            if (sm := srx.match(c)):
                sem_of[tuple(sorted(sm.groupdict().items()))] = c
    parts = []
    for col in T.columns:
        m = rx.match(col)
        if not m:
            continue
        g = m.groupdict()
        cond = g.get("condition") or ds.get("condition", "default")
        rep = g.get("replicate") or ds.get("replicate", 1)
        d = pd.DataFrame({"protein": T[idc].values, "condition": cond, "replicate": rep, "fraction": int(g["fraction"]),
                          "intensity": pd.to_numeric(T[col], errors="coerce").fillna(0).values})
        sc = sem_of.get(tuple(sorted(g.items())))
        if sem_of:
            d["sem"] = pd.to_numeric(T[sc], errors="coerce").fillna(0).values if sc else 0.0
        parts.append(d)
    if not parts:
        raise ValueError(f"protein_matrix: no column matches value_regex {ds['value_regex']!r}")
    out = pd.concat(parts, ignore_index=True)
    keep = ds.get("conditions")
    if keep:
        out = out[out.condition.isin(keep)]
    # a protein with no signal in one run (condition x replicate) of a multi-run matrix is not detected there: drop its
    # rows for that run (an all-zero profile would otherwise be scored with its apex at fraction 1)
    tot = out.groupby(["protein", "condition", "replicate"], sort=False).intensity.transform("sum")
    out = out[tot > 0].reset_index(drop=True)
    return clean_ids(out, ds)


def read_protein_long(ds: dict, resolve) -> pd.DataFrame:
    c = {"protein": "protein", "fraction": "fraction", "intensity": "intensity"} | (ds.get("columns") or {})
    T = pd.read_csv(resolve(ds["path"]), sep="\t", comment="#")
    out = pd.DataFrame({k: T[c[k]].values for k in ("protein", "fraction", "intensity")})
    out["condition"] = T[c["condition"]].values if c.get("condition") in T.columns else ds.get("condition", "default")
    out["replicate"] = T[c["replicate"]].values if c.get("replicate") in T.columns else 1
    if c.get("sem") in T.columns:
        out["sem"] = T[c["sem"]].values
    out["fraction"] = out.fraction.astype(int)
    out = clean_ids(out, ds)
    dup = ds.get("duplicate_rows", "sum")
    if dup == "unique":
        out = out[out.intensity.notna()]
        key = ["protein", "condition", "replicate", "fraction"]
        nu = out.groupby(key, dropna=False).intensity.nunique()
        if (nu > 1).any():
            bad = nu[nu > 1].index[:3].tolist()
            raise ValueError(f"duplicate_rows: unique, but {int((nu > 1).sum())} protein x fraction cells hold different "
                             f"values (e.g. {bad}): not a repeated protein-group quantity")
        out = out.drop_duplicates(key).reset_index(drop=True)
    elif dup != "sum":
        raise ValueError(f"duplicate_rows must be sum | unique (got {dup!r})")
    out["intensity"] = pd.to_numeric(out.intensity, errors="coerce").fillna(0)
    return out


def read_annotation(spec: dict, resolve) -> pd.DataFrame:
    """protein annotation -> DataFrame indexed by accession: mono_kda, gene, tm (bool).
    spec: {path, id, mass_da | mass_kda, gene, tm (column whose non-empty text means transmembrane, or a bool column)}"""
    T = pd.read_csv(resolve(spec["path"]), sep="\t", comment="#")
    A = pd.DataFrame(index=pd.Index(T[spec.get("id", "Entry")].values, name="uniprot"))
    if spec.get("mass_kda"):
        A["mass_da"] = T[spec["mass_kda"]].values.astype(float) * 1e3
        A["mono_kda"] = T[spec["mass_kda"]].values.astype(float)
    else:
        A["mass_da"] = T[spec.get("mass_da", "Mass")].values.astype(float)
        A["mono_kda"] = T[spec.get("mass_da", "Mass")].values / 1e3
    A["gene"] = T[spec["gene"]].values if spec.get("gene") in T.columns else None
    tm = T[spec["tm"]] if spec.get("tm") in T.columns else pd.Series(False, index=T.index)
    A["tm"] = tm.map(lambda v: bool(v) if isinstance(v, (bool, np.bool_)) else isinstance(v, str)).values
    return A[~A.index.duplicated()]


def load(ds: dict, q: dict, resolve) -> tuple[str, pd.DataFrame]:
    """-> ('peptide', PEP table) or ('protein', PROT table)"""
    f = ds["format"]
    if f == "openswath":
        return "peptide", read_openswath(ds, q, resolve)
    if f == "peptide_long":
        return "peptide", read_peptide_long(ds, q, resolve)
    if f == "protein_matrix":
        return "protein", read_protein_matrix(ds, resolve)
    if f == "protein_long":
        return "protein", read_protein_long(ds, resolve)
    raise ValueError(f"unknown SEC-MS format {f!r} (openswath | peptide_long | protein_matrix | protein_long)")
