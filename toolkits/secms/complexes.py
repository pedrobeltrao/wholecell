"""Complex definitions, per-complex co-elution, per-subunit free pools, Complex Portal calibrants.

Complex sources (membership; stoichiometry where known):
  complex_portal  Complex Portal ComplexTab export (species file, e.g. 9606.tsv / 559292.tsv); stoichiometry from the
                  'Expanded participant list' ("P12345(2)|..."; 0 = unknown); isoform suffixes dropped
  humap           hu.MAP 3.0 complexes CSV/TSV (membership only: an id column + a column of space-separated UniProt
                  accessions; optional confidence column / minimum); stoichiometry unknown (1 per member)

Co-elution of a complex (detected subunits only): subunit profiles normalised to their maximum; the co-elution signal is
their minimum per fraction (all subunits present); coelution_peak = its argmax, coelution_height = its maximum;
mean_corr = mean pairwise Pearson r. Free pool of a subunit = intensity eluting below the complex window,
window = app. MW >= max(window_lo x expected complex MW, k_floor x subunit monomer MW), capped at the void fraction's
apparent MW so that complexes larger than the column range count their void intensity as in-complex.

Complex Portal calibrants (column MW calibration, species-specific): heteromeric complexes with full stoichiometry,
every subunit detected and annotated, no TM subunit, expected MW within [mw_min, mw_max]; co-elution height >= h_min
and mean r >= r_min; one per name key (variants merged; highest co-elution kept); point = (coelution_peak,
expected MW). Same selection as the Jurkat 2026-11-08 script (sec_calibrants_peak.tsv).
"""
from __future__ import annotations

import itertools
import re

import numpy as np
import pandas as pd

_CP_RX = re.compile(r"([A-Z][A-Z0-9]{5,9})(?:-\w+)?\((\d+)\)")


def read_complex_portal(path: str) -> list[dict]:
    cp = pd.read_csv(path, sep="\t")
    out = []
    for _, r in cp.iterrows():
        parts = _CP_RX.findall(str(r["Expanded participant list"]))
        if not parts:
            continue
        sub = {}
        for a, s in parts:
            sub[a] = sub.get(a, 0) + int(s)
        out.append(dict(complex_id=r["#Complex ac"], name=r["Recommended name"], source="complex_portal", subunits=sub,
                        stoich_known=not any(int(s) == 0 for _, s in parts)))
    return out


def read_humap(path: str, id_col: str = "HuMAP3_ID", members_col: str = "Uniprot_ACCs", sep: str | None = None,
               confidence_col: str | None = None, confidence_min=None) -> list[dict]:
    T = pd.read_csv(path, sep=sep or ("," if path.endswith(".csv") else "\t"))
    if confidence_col and confidence_min is not None:
        T = T[pd.to_numeric(T[confidence_col], errors="coerce") <= float(confidence_min)]   # hu.MAP: 1 = best
    out = []
    for r in T.itertuples(index=False):
        accs = [a for a in re.split(r"[\s;,]+", str(getattr(r, members_col))) if a]
        if len(accs) >= 2:
            out.append(dict(complex_id=str(getattr(r, id_col)), name=str(getattr(r, id_col)), source="humap",
                            subunits={a: 1 for a in dict.fromkeys(accs)}, stoich_known=False))
    return out


def read_sources(specs: list[dict], resolve) -> list[dict]:
    out = []
    for s in specs or []:
        if s["format"] == "complex_portal":
            out += read_complex_portal(resolve(s["path"]))
        elif s["format"] == "humap":
            out += read_humap(resolve(s["path"]), s.get("id_col", "HuMAP3_ID"), s.get("members_col", "Uniprot_ACCs"),
                              s.get("sep"), s.get("confidence_col"), s.get("confidence_min"))
        else:
            raise ValueError(f"unknown complex source format {s['format']!r} (complex_portal | humap)")
    return out


def coelution(prof: pd.DataFrame, accs: list[str]):
    """legacy-identical: P normalised to max, co = min, peak/height, mean pairwise r"""
    P = prof.loc[accs].values
    with np.errstate(invalid="ignore", divide="ignore"):
        P = P / P.max(axis=1, keepdims=True)
        co = P.min(axis=0)
        corr = np.nanmean([np.corrcoef(P[i], P[j])[0, 1] for i, j in itertools.combinations(range(len(P)), 2)]) \
            if len(P) > 1 else np.nan
    return int(np.argmax(co)) + 1, float(co.max()), float(corr)


def complex_tables(prof, complexes, annot, app, nf, cfg, meta):
    """-> (complex rows DataFrame, subunit rows DataFrame). meta: source_id, condition, replicate"""
    det = set(prof.index)
    cc = cfg.get("complexes", {}) or {}
    win_lo, kfl = float(cc.get("window_lo", 0.5)), float(cc.get("k_floor", 2.0))
    vf = (cfg.get("calibration", {}) or {}).get("void_fraction")
    cap = app[int(vf) - 1] if vf else np.inf
    min_det, min_frac = int(cc.get("min_subunits_detected", 2)), float(cc.get("min_detected_frac", 0.5))
    crow, srow = [], []
    for c in complexes:
        sub = c["subunits"]
        if len(sub) < 2 and not (len(sub) == 1 and list(sub.values())[0] >= 2):
            continue
        d = [a for a in sub if a in det and a in annot.index]
        if len(d) < 1:
            continue
        stoich = {a: (s if s > 0 else 1) for a, s in sub.items()}
        known = all(a in annot.index for a in sub)
        exp = float(sum(annot.mono_kda[a] * stoich[a] for a in sub if a in annot.index)) if known else np.nan
        peak = height = corr = np.nan
        if len(d) >= 2:
            peak, height, corr = coelution(prof, d)
        lo_c = min(win_lo * exp, cap) if np.isfinite(exp) else np.nan
        fr = []
        for a in sub:
            row = dict(complex_id=c["complex_id"], complex_source=c["source"], uniprot=a,
                       gene=annot.gene.get(a) if a in annot.index else None, stoich=sub[a], **meta,
                       detected=a in d, mono_kda=annot.mono_kda.get(a, np.nan) if a in annot.index else np.nan)
            if a in d and np.isfinite(lo_c):
                p = prof.loc[a].values; t = p.sum()
                lo = max(lo_c, min(kfl * annot.mono_kda[a], cap))
                f = float(p[app >= lo].sum() / t) if t > 0 else np.nan
                row.update(in_complex_frac=f, free_frac=1 - f if f == f else np.nan, window_lo_kda=lo,
                           apex_fraction=int(np.argmax(p)) + 1, apex_app_kda=float(app[np.argmax(p)]))
                fr.append(f)
            srow.append(row)
        n_ok = len(d) >= min_det and len(d) / len(sub) >= min_frac
        fa = float(np.nanmedian(fr)) if fr and np.isfinite(fr).any() else np.nan
        crow.append(dict(complex_id=c["complex_id"], species_id="", **meta, n_subunits_detected=len(d),
                         n_subunits=len(sub), coelution_peak=peak, coelution_height=height, mean_corr=corr,
                         assembled_frac=fa, free_frac=1 - fa if fa == fa else np.nan,
                         app_kda=float(app[int(peak) - 1]) if peak == peak else np.nan, expected_kda=exp,
                         notes="" if n_ok else "coverage below min_subunits_detected/min_detected_frac",
                         complex_source=c["source"], name=c["name"], stoich_known=c["stoich_known"],
                         measured=bool(n_ok and fa == fa), window_lo_kda=lo_c,
                         tm_any=any(bool(annot.tm.get(a, False)) for a in sub if a in annot.index)))
    return pd.DataFrame(crow), pd.DataFrame(srow)


def calibrants(prof, complexes, annot, cfg):
    """Complex Portal calibrant points (legacy-identical selection). -> DataFrame cpx, name, n, mw_kda, peak, co_height,
    mcorr, key, kept (kept = inside 2.5 x MAD of a first Theil-Sen fit)."""
    from scipy.stats import theilslopes
    c = (cfg.get("calibration", {}) or {}).get("complex_calibrants", {}) or {}
    det = set(prof.index)
    rows = []
    for x in complexes:
        if x["source"] != "complex_portal" or not x["stoich_known"]:
            continue
        sub = x["subunits"]
        if len(sub) < 2 or not all(a in det for a in sub) or not all(a in annot.index for a in sub):
            continue
        if any(bool(annot.tm[a]) for a in sub):
            continue
        mw = sum(annot.mass_da[a] * s for a, s in sub.items()) / 1e3
        if not (float(c.get("mw_min", 40)) <= mw <= float(c.get("mw_max", 5000))):
            continue
        pk, h, corr = coelution(prof, list(sub))
        rows.append(dict(cpx=x["complex_id"], name=x["name"], n=len(sub), mw_kda=mw, peak=pk, co_height=h, mcorr=corr))
    cal = pd.DataFrame(rows, columns=["cpx", "name", "n", "mw_kda", "peak", "co_height", "mcorr"])
    cal = cal[(cal.co_height >= float(c.get("height_min", 0.5))) & (cal.mcorr >= float(c.get("corr_min", 0.6)))]
    cal["key"] = cal.name.str.replace(c.get("name_key_regex", r",.*variant.*|\(.*|variant.*"), "", regex=True).str.strip()
    cal = cal.sort_values("co_height", ascending=False).drop_duplicates("key")
    if len(cal) >= 3:
        x, y = cal.peak.values, np.log10(cal.mw_kda.values)
        slope, icpt, _, _ = theilslopes(y, x)
        res = y - (icpt + slope * x)
        cal["kept"] = np.abs(res) < 2.5 * np.median(np.abs(res)) / 0.6745
    else:
        cal["kept"] = False
    return cal
