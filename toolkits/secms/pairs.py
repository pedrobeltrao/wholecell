"""Pairwise co-elution (heterodimer evidence): Pearson r of the smoothed profiles of two proteins in one run, for
(a) candidate pairs from a table (e.g. predicted heterodimers; orientation kept as given, r computed exactly as the
legacy np.corrcoef(prof[a], prof[b])[0, 1]) and (b) optionally every detected pair with r >= all_pairs_min_r.
apex_delta = |apex fraction difference|. Replicates: one row per run, plus 'mean' rows (r mean, r_sd, n_rep)."""
from __future__ import annotations

import numpy as np
import pandas as pd


def read_candidates(spec: dict, resolve) -> list[tuple[str, str]]:
    T = pd.read_csv(resolve(spec["path"]), sep=spec.get("sep", "\t"), comment=spec.get("comment"))
    if spec.get("require_true"):
        T = T[T[spec["require_true"]].astype(str).str.lower() == "true"]
    a, b = spec.get("a", "uniprot_a"), spec.get("b", "uniprot_b")
    return list(dict.fromkeys(zip(T[a].astype(str), T[b].astype(str))))


def pair_rows(prof, cand, meta, all_min_r=None, source="candidates"):
    det = prof.index
    pos = {a: i for i, a in enumerate(det)}
    V = prof.values
    apex = V.argmax(axis=1)
    rows = []
    for a, b in cand:
        if a not in pos or b not in pos:
            continue
        with np.errstate(invalid="ignore", divide="ignore"):
            r = float(np.corrcoef(V[pos[a]], V[pos[b]])[0, 1])
        rows.append(dict(uniprot_a=a, uniprot_b=b, **meta, r=r, r_sd=np.nan, n_rep=1,
                         apex_delta=int(abs(apex[pos[a]] - apex[pos[b]])), pair_source=source))
    if all_min_r is not None:
        seen = {(a, b) for a, b in cand} | {(b, a) for a, b in cand}
        ok = V.sum(axis=1) > 0
        with np.errstate(invalid="ignore", divide="ignore"):
            C = np.corrcoef(V[ok])
        ids = np.array(det)[ok]
        ap = apex[ok]
        ii, jj = np.nonzero(np.triu(C >= float(all_min_r), 1))
        for i, j in zip(ii, jj):
            if (ids[i], ids[j]) in seen:
                continue
            rows.append(dict(uniprot_a=ids[i], uniprot_b=ids[j], **meta, r=float(C[i, j]), r_sd=np.nan, n_rep=1,
                             apex_delta=int(abs(ap[i] - ap[j])), pair_source="all_pairs"))
    return rows


def combine_replicates(T: pd.DataFrame) -> pd.DataFrame:
    if T.empty or T.replicate.nunique() < 2:
        return T.iloc[0:0]
    g = T.groupby(["uniprot_a", "uniprot_b", "condition"], sort=True)
    M = g.r.mean().to_frame()
    M["r_sd"] = g.r.std(ddof=1)
    M["n_rep"] = g.size()
    M["apex_delta"] = g.apex_delta.median()
    first = g.first()
    for c in ("source_id", "pair_source"):
        M[c] = first[c]
    M = M.reset_index()
    M["replicate"] = "mean"
    return M.reindex(columns=list(T.columns) + [c for c in M.columns if c not in T.columns])
