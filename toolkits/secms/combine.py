"""Combine the SEC-MS runs of one cell (several conditions) into ONE assembly-state estimate per protein.

Rule (Pedro 2026-11-17; DECISION in PLANNING_NOTES "SEC-MS lab 2026"):
  per protein x threshold k, over the conditions that score the protein (each with its own calibration):
    assembled_frac     = mean of the condition values (equal weight per condition; a protein seen in one condition
                         takes that condition's value)
    assembled_frac_sd  = SD over the conditions (ddof 1; NaN when n_conditions = 1): condition (dis)agreement as the
                         uncertainty of the combined value
    assembled_frac_lo/_hi = mean of the conditions' calibration bounds, widened to include every condition's value
    cond_range         = max - min of the condition values; n_conditions; conditions (labels)
    apex_app_kda       = geometric mean; apex_fraction = from the condition with the largest summed intensity
  pairs: r = mean over conditions where both partners are profiled (r_sd = SD); complexes / subunits: numeric columns
  averaged over the conditions that measured them (measured = any). Rows carry source_id = <cell label>, condition =
  'combined', replicate = 'combined', so that SecView selects them with assembly_options.secms.
Optional `min_fractions`: a condition whose profile of the protein has fewer non-zero raw fractions is ignored for that
protein (unless it is the only one); default 0 = off.

Combination spec (optional; 2026-11-17, final SEC-MS estimates with literature data; `--spec <yaml>` on the CLI):
  label: <cell>_final
  datasets:                                  # by source_id; datasets of the run not listed are left out
    - {source_id: A, weight: 1.0, role: primary}
    - {source_id: B, weight: 1.0, role: fill}   # used only for proteins (pairs, complexes) that no primary dataset scores
    - {source_id: C, weight: 1.0, replicates: split_weight}   # (2026-11-17) the weight is split over the dataset's
                                             # replicate runs that score the key (default: each replicate run gets it)
  discordant_range: 0.3                      # confidence 'discordant' when the primary values span more than this
  per protein x k: assembled_frac = weighted mean over the primary datasets scoring the protein (weights renormalised over
  those present), assembled_frac_sd = reliability-weighted SD over them (= SD ddof 1 for equal weights), range and
  lo/hi as above; proteins scored only by fill datasets take the (weighted) fill value; confidence = high (>= 2 primary,
  range <= discordant_range) | discordant (>= 2 primary, larger range) | single (one primary dataset) | low (fill only)
  | unmeasured (2026-11-23: no dataset can resolve k x monomer; assembled_frac NaN, protein.py measurable range).
  Without a spec every condition is primary with weight 1 (= the equal-weight rule above).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

KEY = ["uniprot", "threshold_k"]


def _prep(T: pd.DataFrame, keys: list, spec: dict | None):
    """single-run rows of the spec's datasets, `_cond`, `_w` (weight), `_fill` (role fill); fill rows dropped for keys that
    a primary dataset covers. Returns (rows, fill_only Series by keys)."""
    T = T[T.replicate.astype(str) != "mean"].copy()
    T["_cond"] = T.source_id.astype(str) + "|" + T.condition.astype(str)
    if spec and spec.get("datasets"):
        D = {str(d["source_id"]): d for d in spec["datasets"]}
        T = T[T.source_id.astype(str).isin(D)].copy()
        T["_w"] = T.source_id.astype(str).map(lambda x: float(D[x].get("weight", 1.0)))
        T["_fill"] = T.source_id.astype(str).map(lambda x: D[x].get("role", "primary") == "fill")
        split = [x for x, d in D.items() if d.get("replicates") == "split_weight"]
        if split:                       # the dataset's weight is shared by its replicate runs scoring the key
            m = T.source_id.astype(str).isin(split)
            n = T[m].groupby([T[m][k] for k in keys] + [T[m].source_id, T[m].condition])["_w"].transform("size")
            T.loc[m, "_w"] = T.loc[m, "_w"] / n
    else:
        T["_w"], T["_fill"] = 1.0, False
    has_primary = (~T["_fill"]).groupby([T[k] for k in keys]).transform("any")
    T = T[~T["_fill"] | ~has_primary].copy()
    fill_only = T.groupby(keys, sort=True)["_fill"].all()
    return T, fill_only


def _wmean(T: pd.DataFrame, keys: list, cols: list) -> pd.DataFrame:
    """weighted mean of cols per key over rows with a value (weights renormalised over those present)"""
    out = {}
    for c in cols:
        x = T[c].astype(float)
        ok = x.notna()
        w = T["_w"].where(ok, 0.0)
        g = pd.DataFrame({"wx": (w * x.fillna(0)), "w": w}).groupby([T[k] for k in keys], sort=True).sum()
        out[c] = (g.wx / g.w).where(g.w > 0)
    return pd.DataFrame(out)


def _wsd(T: pd.DataFrame, keys: list, col: str, mean: pd.Series) -> pd.Series:
    """reliability-weighted SD over datasets (equal weights: SD ddof 1); NaN for a single dataset"""
    m = mean.reindex(pd.MultiIndex.from_frame(T[keys]) if len(keys) > 1 else T[keys[0]]).values
    d = pd.DataFrame({"a": T["_w"] * (T[col].astype(float) - m) ** 2, "v1": T["_w"], "v2": T["_w"] ** 2})
    g = d.groupby([T[k] for k in keys], sort=True).sum()
    den = g.v1 - g.v2 / g.v1
    return np.sqrt(g.a / den).where(den > 0)


def combine_protein(P: pd.DataFrame, label: str, nnz: pd.DataFrame | None = None, intensity: pd.DataFrame | None = None,
                    min_fractions: int = 0, spec: dict | None = None) -> pd.DataFrame:
    """P: secms_protein rows of single runs (one row per uniprot x condition x k). nnz / intensity: optional DataFrames
    uniprot x condition (non-zero raw fractions, summed raw intensity). spec: optional weights / roles (module doc).
    Unmeasured rows (assembled_frac NaN: k x monomer outside a run's calibrated range, protein.py) do not score the
    protein: the other datasets (fill datasets included) decide; a key no dataset measures gets one NaN row."""
    P = P[P.replicate.astype(str) != "mean"]
    um = P[P.assembled_frac.isna()]
    P = P[P.assembled_frac.notna()]
    P, fill_only = _prep(P, KEY, spec)
    if nnz is not None and min_fractions:
        n = P.apply(lambda r: nnz.at[r.uniprot, r._cond] if (r.uniprot in nnz.index and r._cond in nnz.columns) else 0, axis=1)
        nc = P.groupby(KEY)._cond.transform("nunique")
        P = P[(n >= min_fractions) | (nc == 1)]
    g = P.groupby(KEY, sort=True)
    M = _wmean(P, KEY, ["assembled_frac", "assembled_frac_lo", "assembled_frac_hi", "void_frac", "assembled_frac_peaks"]
               + [c for c in ("assembled_frac_mid",) if c in P.columns])
    M["assembled_frac_sd"] = _wsd(P, KEY, "assembled_frac", M.assembled_frac)
    M["cond_range"] = g.assembled_frac.max() - g.assembled_frac.min()
    M["assembled_frac_lo"] = np.minimum(M.assembled_frac_lo, g.assembled_frac.min())
    M["assembled_frac_hi"] = np.maximum(M.assembled_frac_hi, g.assembled_frac.max())
    M["n_conditions"] = g._cond.nunique()
    M["conditions"] = g._cond.agg(lambda s: ";".join(sorted(set(s))))
    M["weight_sum"] = g["_w"].sum()
    M["fill_only"] = fill_only.reindex(M.index).astype(bool)
    rng_max = float((spec or {}).get("discordant_range", 0.3))
    M["confidence"] = np.where(M.fill_only, "low", np.where(M.n_conditions == 1, "single",
                               np.where(M.cond_range > rng_max, "discordant", "high")))
    lw = P["_w"] * np.log(P.apex_app_kda.astype(float))
    gl = pd.DataFrame({"a": lw, "w": P["_w"].where(P.apex_app_kda.notna(), 0)}).groupby([P[k] for k in KEY], sort=True).sum()
    M["apex_app_kda"] = np.exp(gl.a / gl.w)
    for c in ("gene", "mono_kda", "tm", "species", "n_peptides"):
        if c in P.columns:
            M[c] = g[c].first()
    if intensity is not None:
        P["_int"] = [intensity.at[u, c] if (u in intensity.index and c in intensity.columns) else 0 for u, c in zip(P.uniprot, P._cond)]
        best = P.sort_values("_int", ascending=False).drop_duplicates(KEY).set_index(KEY)
    else:
        best = P.drop_duplicates(KEY).set_index(KEY)
    M["apex_fraction"] = best.apex_fraction.reindex(M.index)
    M["n_fractions"] = best.n_fractions.reindex(M.index)
    M["n_peaks"] = best.n_peaks.reindex(M.index)
    M = M.reset_index()
    M["free_frac"] = 1 - M.assembled_frac
    M["source_id"] = label
    M["condition"] = "combined"
    M["replicate"] = "combined"
    M["n_rep"] = M.n_conditions
    M["calibration_id"] = "per_condition"
    M["notes"] = np.where(M.fill_only, "fill dataset only (low confidence)",
                          np.where(M.n_conditions > 1, "weighted mean over datasets" if spec else "mean over conditions",
                                   "single condition"))
    if "measurable" in P.columns:
        M["measurable"] = True
    if len(um):                                    # keys that no dataset of the combination measures: NaN rows
        if spec and spec.get("datasets"):
            um = um[um.source_id.astype(str).isin(spec_ids(spec))]
        U = um.drop_duplicates(KEY).set_index(KEY)
        U = U[~U.index.isin(pd.MultiIndex.from_frame(M[KEY]))]
        if len(U):
            N = pd.DataFrame(index=U.index).reset_index()
            for c in ("gene", "mono_kda", "tm", "species", "n_peptides", "apex_fraction", "apex_app_kda", "n_fractions",
                      "n_peaks"):
                if c in U.columns:
                    N[c] = U[c].values
            N["n_conditions"] = 0
            N["confidence"] = "unmeasured"
            N["measurable"] = False
            N["source_id"], N["condition"], N["replicate"], N["n_rep"] = label, "combined", "combined", 0
            N["calibration_id"] = "per_condition"
            N["notes"] = "unmeasured: k x monomer outside the calibrated range in every dataset"
            M = pd.concat([M, N.reindex(columns=M.columns)], ignore_index=True)
    return M


def combine_pairs(T: pd.DataFrame, label: str, spec: dict | None = None) -> pd.DataFrame:
    if not len(T):
        return T
    cols = list(T.columns)
    K = ["uniprot_a", "uniprot_b"]
    T, _ = _prep(T, K, spec)
    g = T.groupby(K, sort=True)
    M = _wmean(T, K, ["r", "apex_delta"])
    M["r_sd"] = _wsd(T, K, "r", M.r)
    M["n_rep"] = g.size()
    M["pair_source"] = g.pair_source.first()
    M = M.reset_index()
    M["source_id"], M["condition"], M["replicate"] = label, "combined", "combined"
    return M[cols]


def combine_complex(C: pd.DataFrame, S: pd.DataFrame, label: str, spec: dict | None = None):
    ccols, scols = list(C.columns), list(S.columns)
    num = ["coelution_peak", "coelution_height", "mean_corr", "assembled_frac", "app_kda", "n_subunits_detected"]
    if spec and spec.get("datasets"):                           # measured rows: primary datasets first, fill otherwise
        Cm0 = C[C.measured.astype(str) == "True"]
        Cp, _ = _prep(Cm0, ["complex_id"], spec) if len(Cm0) else (Cm0.assign(_w=1.0), None)
        keep = set(spec_ids(spec))
        C = C[(C.replicate.astype(str) != "mean") & C.source_id.astype(str).isin(keep)]
        S = S[(S.replicate.astype(str) != "mean") & S.source_id.astype(str).isin(keep)]
    else:
        C = C[C.replicate.astype(str) != "mean"]
        S = S[S.replicate.astype(str) != "mean"]
        Cp = (C[C.measured.astype(str) == "True"] if "measured" in C.columns else C).assign(_w=1.0)
    first = C.groupby("complex_id", sort=False).first()
    M = first.copy()
    M[num] = M[num].astype(float)
    if len(Cp):
        W = _wmean(Cp, ["complex_id"], num)
        M.loc[W.index, num] = W[num].values
    M["measured"] = C.groupby("complex_id").measured.apply(lambda s: (s.astype(str) == "True").any())
    M["free_frac"] = 1 - M.assembled_frac
    M = M.reset_index()
    M["source_id"], M["condition"], M["replicate"] = label, "combined", "combined"
    Sd = S[S.detected.astype(str) == "True"]
    Sd = _prep(Sd, ["complex_id", "uniprot"], spec)[0] if len(Sd) else Sd.assign(_w=1.0)
    Sf = S.groupby(["complex_id", "uniprot"], sort=False).first()
    Sf[["in_complex_frac", "apex_app_kda"]] = Sf[["in_complex_frac", "apex_app_kda"]].astype(float)
    if len(Sd):
        W = _wmean(Sd, ["complex_id", "uniprot"], ["in_complex_frac", "apex_app_kda"])
        Sf.loc[W.index, ["in_complex_frac", "apex_app_kda"]] = W.values
    Sf["detected"] = S.groupby(["complex_id", "uniprot"]).detected.apply(lambda s: (s.astype(str) == "True").any())
    Sf["free_frac"] = 1 - Sf.in_complex_frac
    Sf = Sf.reset_index()
    Sf["source_id"], Sf["condition"], Sf["replicate"] = label, "combined", "combined"
    return M[ccols], Sf[scols]


def spec_ids(spec: dict | None) -> list:
    return [str(d["source_id"]) for d in ((spec or {}).get("datasets") or [])]


def combine_tables(res: dict, label: str, min_fractions: int = 0, spec: dict | None = None) -> dict:
    """res: toolkit tables of one cell ({name: DataFrame}, as pipeline.run or read back from a run dir) ->
    the same tables with the combined rows appended"""
    out = {k: v for k, v in res.items()}
    pf = res.get("profiles")
    nnz = inten = None
    if pf is not None and len(pf):
        raw = pf[pf.kind == "raw"]
        fc = [c for c in raw.columns if c.startswith("f") and c[1:].isdigit()]
        key = raw.source_id.astype(str) + "|" + raw.condition.astype(str)
        nnz = pd.DataFrame({"uniprot": raw.uniprot, "c": key, "v": (raw[fc] > 0).sum(1)}).pivot_table(
            index="uniprot", columns="c", values="v", aggfunc="max")
        inten = pd.DataFrame({"uniprot": raw.uniprot, "c": key, "v": raw[fc].sum(1)}).pivot_table(
            index="uniprot", columns="c", values="v", aggfunc="max")
    P = res["secms_protein"]
    M = combine_protein(P, label, nnz, inten, min_fractions, spec)
    out["secms_protein"] = pd.concat([P, M.reindex(columns=list(P.columns) + [c for c in M.columns if c not in P.columns])],
                                     ignore_index=True)
    if len(res.get("secms_pairs", [])):
        out["secms_pairs"] = pd.concat([res["secms_pairs"], combine_pairs(res["secms_pairs"], label, spec)], ignore_index=True)
    if len(res.get("secms_complex", [])):
        Cm, Sm = combine_complex(res["secms_complex"], res["secms_complex_subunit"], label, spec)
        out["secms_complex"] = pd.concat([res["secms_complex"], Cm], ignore_index=True)
        out["secms_complex_subunit"] = pd.concat([res["secms_complex_subunit"], Sm], ignore_index=True)
    if len(res.get("secms_calibrants", [])):                    # diagnostic table: the first condition's, relabelled
        C = res["secms_calibrants"]
        s0 = spec_ids(spec)[0] if spec_ids(spec) else C.source_id.iloc[0]
        C0 = C[C.source_id == s0]
        c0 = C0[C0.condition == C0.condition.iloc[0]].copy()
        c0["source_id"], c0["condition"], c0["replicate"] = label, "combined", "combined"
        out["secms_calibrants"] = pd.concat([C, c0], ignore_index=True)
    return out
