"""Protein quantification (top-N peptides) and profile smoothing.

profiles(kind, table, nf, quant, profile) -> {(condition, replicate): {"raw": DataFrame protein x fraction 1..nf,
                                             "smooth": same after the running mean, "sem": DataFrame or None}},
                                            n_peptides (Series per protein)
Top-N (CCprofiler proteinQuantification): per protein, peptides ranked by summed intensity (over all runs of the
dataset with rank_scope=all_runs, the default, so every replicate uses the same peptides; per_run ranks each run
separately); the N best are summed per fraction. Ranks use pandas' 'average' tie rule (legacy-identical).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _topn_one(pep, keep_pairs, nf):
    sel = pep.merge(keep_pairs)
    prof = sel.groupby(["protein", "fraction"]).intensity.sum().unstack(fill_value=0)
    return prof.reindex(columns=range(1, nf + 1), fill_value=0)


def smooth(prof: pd.DataFrame, window: int, mode: str = "reflect") -> pd.DataFrame:
    if not window or int(window) <= 1:
        return prof.astype(float)
    from scipy.ndimage import uniform_filter1d
    return pd.DataFrame(uniform_filter1d(prof.values.astype(float), int(window), axis=1, mode=mode),
                        index=prof.index, columns=prof.columns)


def profiles(kind: str, T: pd.DataFrame, nf: int, quant: dict, profile: dict):
    out, npep = {}, None
    runs = sorted({(c, r) for c, r in zip(T.condition, T.replicate)}, key=lambda x: (str(x[0]), str(x[1])))
    if kind == "peptide":
        n = int(quant.get("top_n", 3))
        # per run: summed intensity per protein x peptide x fraction (legacy: several precursors / charge states summed)
        per_run = {}
        for c, r in runs:
            sub = T[(T.condition == c) & (T.replicate == r)]
            per_run[(c, r)] = sub.groupby(["protein", "peptide", "fraction"]).intensity.sum().reset_index()
        scope = quant.get("rank_scope", "all_runs")

        def ranked(pep):
            tot = pep.groupby(["protein", "peptide"]).intensity.sum().rename("t").reset_index()
            tot["r"] = tot.groupby("protein").t.rank(ascending=False)
            return tot

        allpep = pd.concat(list(per_run.values()), ignore_index=True) if len(per_run) > 1 else next(iter(per_run.values()))
        tot_all = ranked(allpep)
        npep = tot_all.groupby("protein").size()
        ok = set(npep.index[npep >= int(quant.get("min_peptides", 1))])
        for key, pep in per_run.items():
            tot = tot_all if scope == "all_runs" else ranked(pep)
            keep = tot[(tot.r <= n) & tot.protein.isin(ok)][["protein", "peptide"]]
            raw = _topn_one(pep, keep, nf)
            out[key] = {"raw": raw, "sem": None}
    else:
        npep = pd.Series(dtype=float)
        for c, r in runs:
            sub = T[(T.condition == c) & (T.replicate == r)]
            raw = sub.pivot_table(index="protein", columns="fraction", values="intensity", aggfunc="sum", fill_value=0)
            raw = raw.reindex(columns=range(1, nf + 1), fill_value=0)
            sem = None
            if "sem" in sub.columns:
                sem = sub.pivot_table(index="protein", columns="fraction", values="sem", aggfunc="sum", fill_value=0)
                sem = sem.reindex(index=raw.index, columns=range(1, nf + 1), fill_value=0)
            out[(c, r)] = {"raw": raw, "sem": sem}
    for v in out.values():
        v["smooth"] = smooth(v["raw"], profile.get("smooth_window", 1), profile.get("smooth_mode", "reflect"))
        if v["sem"] is not None:
            v["sem"] = v["sem"].reindex(v["raw"].index)
    return out, npep
