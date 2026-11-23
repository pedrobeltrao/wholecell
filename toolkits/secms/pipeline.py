"""SEC-MS toolkit pipeline: one run config -> the standard SEC-MS tables (docs/reference/inventory_schema.md §4a).

Run config (YAML or dict; the inventory adapter `secms_toolkit` passes inventory.params.secms_toolkit):
  species: human | yeast                       species settings (species/<name>.yaml): Complex Portal taxon, anchors, ...
  preset: default | legacy_v7                  presets/<name>.yaml (default = the same parameters for every cell)
  overrides: {}                                deep-merged on the preset; must be justified by a DECISION (recorded)
  annotation: {path, id, mass_da, gene, tm}    UniProt TSV of the species (column names default from the species file)
  complexes: [{format: complex_portal, path}, {format: humap, path, ...}]
  pair_candidates: {path, a, b, require_true}  optional (heterodimer evidence for these pairs)
  datasets:                                    one entry per SEC-MS dataset (all processed with the same parameters)
    - source_id, format (openswath | peptide_long | protein_matrix | protein_long), path, n_fractions,
      fraction_map (openswath), value_regex / sem_regex (protein_matrix), condition, void_fraction (column property),
      standards: {path} (optional column standards: fraction, kda[, name]), conditions (subset), notes,
      calibration_from: <source_id of an earlier dataset>   (dataset property: same column and fraction numbering;
      that dataset's calibration is used, this dataset's own fits are reported as checks, its QC does not stop the run)
Steps per dataset and run (condition x replicate): load -> top-N quantification -> smoothing -> calibration (primary,
QC gate, fallback, cross-checks) -> per-protein assembled fractions (k x monomer; NaN = unmeasured outside the calibrated
range, protein.py) + peaks + uncertainty -> per-complex
co-elution and per-subunit free pools -> pairwise co-elution of candidate pairs; then replicate means.
`run(cfg, resolve)` returns {table name: DataFrame} + diagnostics; `write(result, out_dir)` writes them.
"""
from __future__ import annotations

import copy
import json
import os
import zlib

import numpy as np
import pandas as pd
import yaml

from wholecell.inventory import schema
from . import calib, complexes as cx, io, pairs, protein, quant

HERE = os.path.dirname(os.path.abspath(__file__))
SPECIES = ("human", "yeast")
TABLES = {"secms_protein": schema.SECMS_PROTEIN + schema.SECMS_PROTEIN_EXTRA,
          "secms_complex": schema.SECMS_COMPLEX + schema.SECMS_COMPLEX_EXTRA,
          "secms_complex_subunit": schema.SECMS_COMPLEX_SUBUNIT,
          "secms_pairs": schema.SECMS_PAIRS,
          "secms_peaks": schema.SECMS_PEAKS,
          "secms_calibration": schema.SECMS_CALIBRATION,
          "secms_calibrants": schema.SECMS_CALIBRANTS}


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def load_preset(name: str) -> dict:
    p = yaml.safe_load(open(os.path.join(HERE, "presets", f"{name}.yaml")))
    parent = p.pop("extends", None)
    return _merge(load_preset(parent.replace(".yaml", "")), p) if parent else p


def load_species(name: str) -> dict:
    if name not in SPECIES:
        raise ValueError(f"SEC-MS toolkit species must be one of {SPECIES} (got {name!r}); another species needs a "
                         f"species file, a Complex Portal export and (optionally) curated anchors, see the toolkit README")
    return yaml.safe_load(open(os.path.join(HERE, "species", f"{name}.yaml")))


def params_of(run_cfg: dict) -> dict:
    return _merge(load_preset(run_cfg.get("preset", "default")), run_cfg.get("overrides") or {})


def anchor_file(sp, P=None):
    """the species' anchor table of the preset's anchor-set version (calibration.anchor_set; none = species `anchors`)"""
    ver = ((P or {}).get("calibration") or {}).get("anchor_set")
    if isinstance(ver, dict):                                     # per-species versions {species: version}
        ver = ver.get(sp.get("species"))
    if ver is None:
        return sp.get("anchors")
    sets = sp.get("anchor_sets") or {}
    if ver not in sets:
        raise ValueError(f"anchor_set {ver!r} not defined for species {sp.get('species')} (have {sorted(sets)})")
    return sets[ver]


def _anchor_set_of(sp, P=None):
    ver = ((P or {}).get("calibration") or {}).get("anchor_set")
    return ver.get(sp.get("species")) if isinstance(ver, dict) else ver


def _anchors(sp, P=None):
    f = anchor_file(sp, P)
    if not f:
        return pd.DataFrame(columns=["gene", "native_kda"])
    return pd.read_csv(os.path.join(HERE, "data", f), sep="\t", comment="#")


def _points(method, prof, P, ds, sp, annot, cplx, resolve, cal_cache):
    if method == "anchors":
        return calib.anchor_points(prof, _anchors(sp, P), annot)
    if method == "complex_portal":
        cal = cal_cache["cal"]
        pts = pd.DataFrame(dict(kind="complex_calibrant", id=cal.cpx.values, gene="", fraction=cal.peak.values.astype(float),
                                kda=cal.mw_kda.values))
        vf = ds.get("void_fraction")
        return pts[pts.fraction > float(vf)] if vf else pts
    if method == "standards":
        if not ds.get("standards"):
            return pd.DataFrame(columns=["kind", "id", "gene", "fraction", "kda"])
        return calib.standard_points(pd.read_csv(resolve(ds["standards"]["path"]), sep="\t", comment="#"))
    raise ValueError(f"unknown calibration method {method!r} (complex_portal | anchors | standards)")


def _qc(info, qc):
    why = []
    if info["n_kept"] < int(qc.get("min_points") or 0):
        why.append(f"n_kept {info['n_kept']} < {qc.get('min_points')}")
    if qc.get("max_spearman_rho") is not None and (info["spearman_rho"] is None or info["spearman_rho"] > qc["max_spearman_rho"]):
        why.append(f"rho {info['spearman_rho']} > {qc['max_spearman_rho']}")
    if qc.get("max_resid_fold") is not None and info["resid_fold"] > qc["max_resid_fold"]:
        why.append(f"residual x{info['resid_fold']:.2f} > x{qc['max_resid_fold']}")
    return why


def calibrate(prof, P, ds, sp, annot, cplx, resolve, nf, cal_cache, shared=None):
    C = dict(P["calibration"]); C["void_fraction"] = ds.get("void_fraction")
    qc = C.get("qc", {}) or {}
    fits, pts_all, log = {}, [], []
    methods = [C["primary"]] + [m for m in [C.get("fallback")] + list(C.get("checks") or []) if m and m != C["primary"]]
    for m in dict.fromkeys(methods):
        pts = _points(m, prof, P, ds, sp, annot, cplx, resolve, cal_cache)
        if len(pts) < 3:
            log.append(f"{m}: {len(pts)} points, not fitted")
            continue
        info, pp = calib.fit(pts, nf, dict(C, min_points=3, sources=[m]))
        info["qc_fail"] = _qc(info, qc)
        fits[m] = info
        pts_all.append(pp.assign(method=m))
    use = C["primary"]
    if shared is not None:                                        # calibration_from: the column's calibration of
        src = ds["calibration_from"]                              # another dataset (same column, same fractions)
        fits[f"from:{src}"] = dict(shared, qc_fail=[])
        pts_all.append(shared["_points"].assign(method=f"from:{src}"))
        log.append(f"calibration_from {src} (own {use} fit: {'; '.join(fits[use]['qc_fail']) or 'QC ok'})"
                   if use in fits else f"calibration_from {src} (own {use}: not fitted)")
        use = f"from:{src}"
    elif use not in fits or fits[use]["qc_fail"]:
        why = fits[use]["qc_fail"] if use in fits else ["not fitted"]
        fb = C.get("fallback")
        if fb and fb in fits and not fits[fb]["qc_fail"]:
            log.append(f"primary {use} failed QC ({'; '.join(why)}): fallback {fb}")
            use = fb
        elif use in fits and not qc.get("max_spearman_rho") and not qc.get("max_resid_fold"):
            pass                                                  # no QC gate configured (legacy preset)
        else:
            raise RuntimeError(f"SEC calibration failed for {ds['source_id']}: primary {C['primary']} ({'; '.join(why)}), "
                               f"fallback {fb} {'QC: ' + '; '.join(fits[fb]['qc_fail']) if fb in fits else 'not available'}")
    info = {k: v for k, v in dict(fits[use], used=use, log=log).items() if k != "_points"}
    app = calib.app_kda(info, nf)
    # cross-check: every other method's points under the calibration used
    checks = {}
    for m, f in fits.items():
        pp = [x for x in pts_all if x.method.iloc[0] == m][0]
        pred = np.log10(app[np.clip(np.round(pp.fraction.values).astype(int), 1, nf) - 1])
        r = np.log10(pp.kda.values) - pred
        checks[m] = dict(n=int(len(pp)), own_fit={k: f[k] for k in ("icpt", "slope", "spearman_rho", "resid_fold", "loo_fold", "n_kept") if k in f},
                         qc_fail=f["qc_fail"], median_abs_fold_under_used=float(10 ** np.median(np.abs(r))))
    info["checks"] = checks
    pts = pd.concat(pts_all, ignore_index=True) if pts_all else pd.DataFrame()
    if len(pts):
        pts["role"] = np.where(pts.method == use, "used", "check")
    return info, app, pts


def run_dataset(ds, run_cfg, P, sp, annot, cplx, resolve, cand, seed_base=0, cal_store=None):
    q = P["quant"]
    kind, T = io.load(ds, q, resolve)
    nf = int(ds["n_fractions"])
    prof_runs, npep = quant.profiles(kind, T, nf, q, P["profile"])
    out = {k: [] for k in TABLES}
    out["profiles"] = []
    diag = dict(source_id=ds["source_id"], input_kind=kind, n_fractions=nf, runs={})
    for (cond, rep), pr in prof_runs.items():
        prof = pr["smooth"]
        meta = dict(source_id=ds["source_id"], condition=ds.get("condition", cond) if cond == "default" else cond,
                    replicate=rep)
        cal_cache = {"cal": cx.calibrants(prof, cplx, annot, P)}
        shared = None
        if ds.get("calibration_from"):
            shared = (cal_store or {}).get(ds["calibration_from"])
            if shared is None:
                raise ValueError(f"{ds['source_id']}: calibration_from {ds['calibration_from']!r} is not an earlier dataset "
                                 f"of this run config with a single run")
        info, app, pts = calibrate(prof, P, ds, sp, annot, cplx, resolve, nf, cal_cache, shared)
        if cal_store is not None and len(prof_runs) == 1:
            u = pts[pts.role == "used"] if len(pts) else pts
            cal_store[ds["source_id"]] = dict({k: v for k, v in info.items() if k not in ("checks", "log", "used")},
                                              _points=u.drop(columns=["method", "role"], errors="ignore"))
        cal_id = f"{info['used']}:{ds['source_id']}:{meta['condition']}:{rep}"
        if len(pts):
            out["secms_calibration"].append(pts.assign(**meta, calibration_id=cal_id))
        rng = np.random.default_rng([int((P["uncertainty"].get("bootstrap") or {}).get("seed", 0)), seed_base,
                                     zlib.crc32(str((cond, rep)).encode())]) if pr["sem"] is not None else None
        Pd = dict(P, calibration=dict(P["calibration"], void_fraction=ds.get("void_fraction")))
        rows, peaks = protein.protein_rows(prof, annot, app, info, Pd, meta, npep, pr["sem"], rng)
        Pt = pd.DataFrame(rows)
        if len(Pt):
            Pt["n_fractions"] = nf
            Pt["calibration_id"] = cal_id
            Pt["species"] = sp["species"]
            Pt["notes"] = np.where(Pt.tm, "TM protein", "")
            um = ~Pt.measurable.astype(bool)
            Pt.loc[um, "notes"] = (Pt.notes[um] + "; ").str.lstrip("; ") + "unmeasured: k x monomer outside the calibrated range"
        out["secms_protein"].append(Pt)
        out["secms_peaks"].append(pd.DataFrame(peaks))
        Ct, St = cx.complex_tables(prof, cplx, annot, app, nf, Pd, meta)
        cal = cal_cache["cal"]
        if len(Ct):
            rank = {c: i for i, c in enumerate(cal.cpx)}
            Ct["calibrant"] = Ct.complex_id.isin(rank)
            Ct["calibrant_rank"] = Ct.complex_id.map(rank)
            Ct["calibrant_kept"] = Ct.complex_id.map(dict(zip(cal.cpx, cal.kept)))
            Ct["calibration_id"] = cal_id
            Ct["species"] = sp["species"]
        out["secms_complex"].append(Ct)
        out["secms_complex_subunit"].append(St)
        out["secms_calibrants"].append(cal.assign(**meta))
        out["secms_pairs"].append(pd.DataFrame(pairs.pair_rows(prof, cand, meta, (P.get("pairs") or {}).get("all_pairs_min_r"))))
        for knd, M in (("smooth", prof), ("raw", pr["raw"])):
            W = M.copy(); W.columns = [f"f{c}" for c in W.columns]
            W.insert(0, "kind", knd); W.insert(0, "replicate", rep); W.insert(0, "condition", meta["condition"])
            W.insert(0, "source_id", ds["source_id"]); W.index.name = "uniprot"
            out["profiles"].append(W.reset_index())
        k0 = P["assembled"].get("primary_k", 2)
        a0 = Pt[Pt.threshold_k == k0] if len(Pt) else Pt
        diag["runs"][f"{meta['condition']}|{rep}"] = dict(
            calibration=info, n_proteins_profiled=int(len(prof)), n_proteins_scored=int(a0.uniprot.nunique()) if len(a0) else 0,
            median_assembled_frac_k=float(a0.assembled_frac.median()) if len(a0) else None,
            n_unmeasured_k=int((~a0.measurable.astype(bool)).sum()) if len(a0) else 0,
            n_complexes_measured=int(Ct.measured.sum()) if len(Ct) else 0, n_calibrants=int(len(cal)),
            n_pairs=int(sum(len(x) for x in out["secms_pairs"][-1:])))
    res = {}
    for k, parts in out.items():
        parts = [x for x in parts if x is not None and len(x)]
        res[k] = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=TABLES.get(k, []))
    # replicate means
    if len(res["secms_protein"]) and res["secms_protein"].replicate.nunique() > 1:
        res["secms_protein"] = pd.concat([res["secms_protein"], protein.combine_replicates(res["secms_protein"])], ignore_index=True)
    for k, keys, num in (("secms_complex", ["complex_id", "condition"], ["coelution_height", "mean_corr", "assembled_frac",
                                                                          "app_kda"]),
                         ("secms_complex_subunit", ["complex_id", "uniprot", "condition"], ["in_complex_frac", "apex_app_kda"])):
        T = res[k]
        if len(T) and T.replicate.nunique() > 1:
            g = T.groupby(keys, sort=False)
            M = g.first()
            M[num] = g[num].mean()
            if "measured" in T.columns:
                M["measured"] = g.measured.all()
            M = M.reset_index()
            M["replicate"] = "mean"
            if "free_frac" in M.columns:
                M["free_frac"] = 1 - M[num[2] if k == "secms_complex" else "in_complex_frac"]
            res[k] = pd.concat([T, M[T.columns]], ignore_index=True)
    if len(res["secms_pairs"]) and res["secms_pairs"].replicate.nunique() > 1:
        res["secms_pairs"] = pd.concat([res["secms_pairs"], pairs.combine_replicates(res["secms_pairs"])], ignore_index=True)
    return res, diag


def run(run_cfg: dict, resolve=lambda x: x) -> dict:
    sp = load_species(run_cfg.get("species"))
    P = params_of(run_cfg)
    ann = dict(sp.get("annotation_defaults", {})) | dict(run_cfg["annotation"])
    annot = io.read_annotation(ann, resolve)
    cplx = []
    for s in run_cfg.get("complexes") or []:
        if s["format"] == "humap" and not sp.get("humap_allowed"):
            raise ValueError(f"hu.MAP is human-only (species {sp['species']})")
        if s["format"] == "complex_portal":
            tx = pd.read_csv(resolve(s["path"]), sep="\t", usecols=["Taxonomy identifier"])["Taxonomy identifier"]
            if len(tx) and not (tx.astype(str) == str(sp["complex_portal_taxon"])).all():
                raise ValueError(f"Complex Portal file {s['path']} is not taxon {sp['complex_portal_taxon']} ({sp['species']})")
        cplx += cx.read_sources([s], resolve)
    cand = pairs.read_candidates(run_cfg["pair_candidates"], resolve) if run_cfg.get("pair_candidates") else []
    allres, diags, cal_store = {}, [], {}
    for i, ds in enumerate(run_cfg["datasets"]):
        res, d = run_dataset(ds, run_cfg, P, sp, annot, cplx, resolve, cand, seed_base=i, cal_store=cal_store)
        diags.append(d)
        for k, v in res.items():
            allres.setdefault(k, []).append(v)
    out = {k: pd.concat(v, ignore_index=True) for k, v in allres.items()}
    for k, cols in TABLES.items():
        out[k] = out[k].reindex(columns=cols + [c for c in out[k].columns if c not in cols]) if k in out else pd.DataFrame(columns=cols)
    pf = out["profiles"]
    out["profiles"] = pf[schema.SECMS_PROFILE_ID + [c for c in pf.columns if c not in schema.SECMS_PROFILE_ID]]
    out["_diagnostics"] = dict(species=sp["species"], preset=run_cfg.get("preset", "default"),
                               anchor_set=_anchor_set_of(sp, P), anchor_file=anchor_file(sp, P),
                               overrides=run_cfg.get("overrides") or {}, params=P, datasets=diags,
                               n_complex_definitions=len(cplx), n_pair_candidates=len(cand))
    return out


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def diagnostics_json(d) -> str:
    return json.dumps(d, indent=1, default=_json_default)


def write(result: dict, out_dir: str, writer=None):
    """writer(df, rel) defaults to a plain TSV writer into out_dir (the inventory adapter passes ctx.write_tsv)"""
    os.makedirs(out_dir, exist_ok=True)
    if writer is None:
        def writer(df, rel):
            df.to_csv(os.path.join(out_dir, rel), sep="\t", index=False)
    for k, df in result.items():
        if k.startswith("_"):
            continue
        writer(df, f"{k}.tsv")
