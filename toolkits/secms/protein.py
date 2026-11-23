"""Per-protein assembled fractions, peak detection and uncertainty.

assembled_frac(k) = intensity at apparent MW >= k x monomer MW / total intensity (fraction centres; legacy expression),
for every k in `assembled.thresholds_k`. void_frac = intensity at fractions <= the void fraction (aggregates, membranes,
very large assemblies; counted as assembled, reported separately). Peak-level estimator: peaks of the smoothed profile
(scipy find_peaks: height >= min_rel_height x max, prominence >= min_rel_prominence x max, distance >= min_distance),
the profile split at the minima between consecutive peaks; assembled_frac_peaks = area of peaks whose apex is at
>= primary k x monomer. Uncertainty: calibration (assembled_frac_lo / _hi with the apparent MW scale shifted by -/+ one
residual SD of the calibration), replicates (assembled_frac_sd over replicate profiles, written on the 'mean' rows), and,
for matrices with a standard error per value, a seeded parametric bootstrap (assembled_frac_sd on single-run rows).
Measurable range (config `assembled.range_margin_sd`, 2026-11-23; null = off, the legacy behaviour): a_k is defined only
when the column can separate the monomer from a k-mer, i.e. the monomer mass is at least 10^(margin x calibration
residual SD) x the apparent mass of the last fraction (the monomer peak elutes inside the collected fractions) and k x
monomer is at most 10^-(margin x SD) x the apparent mass of the first non-void fraction. Otherwise the row is
`measurable` = False and every a_k column is NaN ("unmeasured"): without the rule a protein with k x monomer below the
last fraction's apparent mass scores a_k = 1 by construction (yeast lab NaCl: 2 x monomer < 36.9 kDa).
Midpoint share (2026-11-23): assembled_frac_mid = intensity at apparent MW >= sqrt(k (k - 1)) x monomer, the geometric
midpoint between a (k-1)-mer and a k-mer (k = 2: 1.41 x monomer), i.e. the share closer on the log-mass scale to >= k
copies than to k - 1. a_k at k x monomer puts the boundary on the k-mer itself, so a pure dimer reads ~0.5 at k = 2 on an
unbiased calibration; the assembly rules can take the oligomer share of a protein of known homo-oligomeric state from the
midpoint (wholecell/inventory/secms.py, share source secms_midpoint). The bound fraction stays a_k.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def _afrac(p, app, m, k):
    t = p.sum()
    return p[app >= k * m].sum() / t


def measurable(app, m, k, sd_log10, margin_sd, vf=None) -> bool:
    """True when k x monomer can be resolved by this calibrated column (module doc); margin_sd None = always True"""
    if margin_sd is None:
        return True
    w = 10 ** (float(margin_sd) * float(sd_log10 or 0.0))
    top = app[int(vf)] if vf and int(vf) < len(app) else app[0]
    return bool(m >= app[-1] * w and k * m <= top / w)


def protein_rows(prof, annot, app, info, cfg, meta, npep=None, sem=None, rng=None):
    ac = cfg.get("assembled", {}) or {}
    ks = list(ac.get("thresholds_k", [2, 3, 4, 6]))
    k0 = ac.get("primary_k", ks[0])
    vf = (cfg.get("calibration", {}) or {}).get("void_fraction")
    pk = cfg.get("peaks", {}) or {}
    unc = cfg.get("uncertainty", {}) or {}
    sd = info.get("resid_sd_log10", 0.0) if unc.get("calibration_sd", True) else 0.0
    margin = ac.get("range_margin_sd")
    app_lo, app_hi = app * 10 ** (-sd), app * 10 ** sd
    nb = int((unc.get("bootstrap", {}) or {}).get("n", 0)) if sem is not None else 0
    rows, peaks = [], []
    for a in prof.index:
        if a not in annot.index:
            continue
        m = annot.mono_kda[a]
        p = prof.loc[a].values
        t = p.sum()
        if t <= 0:
            continue
        ia = int(np.argmax(p))
        base = dict(uniprot=a, gene=annot.gene.get(a), **meta, apex_fraction=ia + 1, apex_app_kda=app[ia], mono_kda=m)
        vfrac = float(p[:int(vf)].sum() / t) if vf else np.nan
        pl = _peaks(p, app, m, k0, vf, pk)
        for j, q in enumerate(pl):
            peaks.append(dict(uniprot=a, **meta, peak=j + 1, **q))
        apk = float(sum(q["area_frac"] for q in pl if q["cls"] != "monomer")) if pl else np.nan
        for k in ks:
            ok = measurable(app, m, k, info.get("resid_sd_log10", 0.0), margin, vf)
            f = _afrac(p, app, m, k) if ok else np.nan
            r = dict(base, assembled_frac=f, free_frac=1 - f, threshold_k=k,
                     assembled_frac_lo=_afrac(p, app_lo, m, k) if ok else np.nan,
                     assembled_frac_hi=_afrac(p, app_hi, m, k) if ok else np.nan,
                     assembled_frac_sd=np.nan, n_rep=1, n_peptides=int(npep.get(a, 0)) if npep is not None and len(npep) else np.nan,
                     void_frac=vfrac, n_peaks=len(pl), assembled_frac_peaks=apk if (k == k0 and ok) else np.nan,
                     tm=bool(annot.tm.get(a, False)), measurable=ok,
                     assembled_frac_mid=_afrac(p, app, m, np.sqrt(k * (k - 1))) if ok else np.nan)
            if ok and nb and sem is not None and a in sem.index:
                s = sem.loc[a].values
                draws = [_afrac(np.clip(p + s * rng.standard_normal(len(p)), 0, None) + 1e-300, app, m, k) for _ in range(nb)]
                r["assembled_frac_sd"] = float(np.std(draws))
            rows.append(r)
    return rows, peaks


def _peaks(p, app, m, k0, vf, pk):
    from scipy.signal import find_peaks
    mx = p.max()
    if mx <= 0:
        return []
    idx, _ = find_peaks(np.r_[0, p, 0], height=float(pk.get("min_rel_height", 0.05)) * mx,
                        prominence=float(pk.get("min_rel_prominence", 0.1)) * mx, distance=int(pk.get("min_distance", 2)))
    idx = idx - 1
    if not len(idx):
        idx = np.array([int(np.argmax(p))])
    cuts = [0] + [int(idx[i] + np.argmin(p[idx[i]:idx[i + 1] + 1])) for i in range(len(idx) - 1)] + [len(p)]
    t = p.sum()
    out = []
    for i, ap in enumerate(idx):
        lo, hi = cuts[i], cuts[i + 1]
        cls = "void" if vf and ap + 1 <= int(vf) else ("assembled" if app[ap] >= k0 * m else "monomer")
        out.append(dict(apex_fraction=int(ap) + 1, left=lo + 1, right=hi, area_frac=float(p[lo:hi].sum() / t),
                        app_kda=float(app[ap]), ratio_to_mono=float(app[ap] / m), cls=cls))
    return out


def combine_replicates(P: pd.DataFrame, keys=("uniprot", "condition", "threshold_k")) -> pd.DataFrame:
    """mean rows (replicate='mean') over the replicate rows of each protein x condition x threshold"""
    if P.replicate.nunique() < 2:
        return P.iloc[0:0]
    num = ["assembled_frac", "assembled_frac_lo", "assembled_frac_hi", "void_frac", "assembled_frac_peaks", "apex_app_kda"]
    num += [c for c in ("assembled_frac_mid",) if c in P.columns]
    g = P.groupby(list(keys), sort=True)
    M = g[num].mean()
    M["assembled_frac_sd"] = g.assembled_frac.std(ddof=1)
    M["n_rep"] = g.size()
    if "measurable" in P.columns:                     # measurable in any replicate (the mean is over those that are)
        M["measurable"] = g.measurable.any()
    first = g.first()
    for c in P.columns:
        if c not in M.columns and c not in keys:
            M[c] = first[c]
    M = M.reset_index()
    M["replicate"] = "mean"
    M["free_frac"] = 1 - M.assembled_frac
    return M[P.columns]
