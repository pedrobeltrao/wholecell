"""Fraction -> apparent MW calibration with diagnostics.

Model: log10(MW kDa) = icpt + slope x fraction, robust Theil-Sen fit (CCprofiler calibrateMW fits a line on log MW too).
Point sources (config `calibration.sources`, combined in one fit):
  standards  column standards {fraction, kda} (the usual globular protein ladder of the column / of the run)
  anchors    curated proteins of well-established native mass (organism file: gene [, uniprot], native_kda); point =
             apex fraction of the protein's (smoothed) profile in this dataset
Optional outlier rejection (|residual| > outlier_mad x MAD/0.6745 of the first fit, then refit). Diagnostics: n, Spearman
rho, residual SD (log10 and fold), leave-one-out residuals (prediction error of each point by a fit without it), the
apparent MW at the first / last fraction, the void fraction's apparent MW. The complex co-elution calibrants
(complexes.py) are a cross-check only: in mild lysis they co-elute partially (Jurkat 2026-11-08: rho -0.34).
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _fit(x, y):
    from scipy.stats import theilslopes
    slope, icpt, _, _ = theilslopes(y, x)
    return float(icpt), float(slope)


def anchor_points(prof: pd.DataFrame, anchors: pd.DataFrame, annot: pd.DataFrame) -> pd.DataFrame:
    """apex fraction of each anchor protein in `prof`. Gene -> accession over the profiled proteins with an annotation
    (the last accession wins for a gene with several, as in the legacy code); an explicit `uniprot` column wins.
    Native mass: `native_kda`, or `oligomer_n` x the annotation's monomer mass (homomers)."""
    g2a = {annot.gene.get(a): a for a in prof.index if a in annot.index}
    rows = []
    for r in anchors.itertuples(index=False):
        acc = getattr(r, "uniprot", None)
        acc = acc if isinstance(acc, str) and acc in prof.index else g2a.get(r.gene)
        if acc is None:
            continue
        kda = getattr(r, "native_kda", np.nan)
        n = getattr(r, "oligomer_n", np.nan)
        if not (isinstance(kda, (int, float)) and kda == kda) and n == n and acc in annot.index:
            kda = float(n) * float(annot.mono_kda[acc])            # homomer: n x monomer mass of this annotation
        if not (kda == kda):
            continue
        rows.append(dict(kind="anchor", id=acc, gene=r.gene, fraction=int(np.argmax(prof.loc[acc].values)) + 1,
                         kda=float(kda)))
    return pd.DataFrame(rows, columns=["kind", "id", "gene", "fraction", "kda"])


def standard_points(std: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(dict(kind="standard", id=std["name"] if "name" in std.columns else std.index.astype(str),
                             gene="", fraction=std.fraction.astype(float), kda=std.kda.astype(float)))


def fit(points: pd.DataFrame, nf: int, cfg: dict) -> tuple[dict, pd.DataFrame]:
    if len(points) < int(cfg.get("min_points", 3)):
        raise ValueError(f"calibration: only {len(points)} points (min_points {cfg.get('min_points', 3)})")
    from scipy.stats import spearmanr
    x = points.fraction.values.astype(float)
    y = np.log10(points.kda.values.astype(float))
    icpt, slope = _fit(x, y)
    res = y - (icpt + slope * x)
    kept = np.ones(len(x), bool)
    mad = cfg.get("outlier_mad")
    if mad:
        kept = np.abs(res) < float(mad) * np.median(np.abs(res)) / 0.6745
        icpt, slope = _fit(x[kept], y[kept])
        res = y - (icpt + slope * x)
    loo = np.full(len(x), np.nan)
    idx = np.nonzero(kept)[0]
    if len(idx) > 3:
        for i in idx:
            m = kept.copy(); m[i] = False
            a, b = _fit(x[m], y[m])
            loo[i] = y[i] - (a + b * x[i])
    P = points.copy()
    P["log10_fit"] = icpt + slope * x
    P["residual"] = res
    P["kept"] = kept
    P["loo_residual"] = loo
    rk = res[kept]
    sd = float(rk.std())
    vf = cfg.get("void_fraction")
    info = dict(method="theil_sen", sources=list(cfg.get("sources", ["anchors"])), icpt=icpt, slope=slope,
                n_points=int(len(x)), n_kept=int(kept.sum()),
                spearman_rho=float(spearmanr(x[kept], y[kept])[0]) if kept.sum() > 2 else None,
                resid_sd_log10=sd, resid_fold=float(10 ** sd),
                loo_rmse_log10=float(np.sqrt(np.nanmean(loo ** 2))) if np.isfinite(loo).any() else None,
                kda_first_fraction=float(10 ** (icpt + slope * 1)), kda_last_fraction=float(10 ** (icpt + slope * nf)),
                void_fraction=vf, kda_void=float(10 ** (icpt + slope * vf)) if vf else None,
                fraction_range_of_points=[float(x.min()), float(x.max())],
                kda_range_of_points=[float(points.kda.min()), float(points.kda.max())])
    if info["loo_rmse_log10"] is not None:
        info["loo_fold"] = float(10 ** info["loo_rmse_log10"])
    return info, P


def app_kda(info: dict, nf: int, shift_log10: float = 0.0) -> np.ndarray:
    """apparent MW (kDa) at fraction centres 1..nf (legacy expression, bit-identical)"""
    a = 10 ** (info["icpt"] + info["slope"] * np.arange(1, nf + 1))
    return a if not shift_log10 else a * 10 ** shift_log10
