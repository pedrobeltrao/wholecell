"""Score the yeast X5 reproduction under the criteria addendum of 2026-11-17 (seed-spread rule).

    python plugins/yeast/repro/seedspread_x5.py --s0 <build> --s1 <build> --s2 <build> --x5 <X5 run dir> --out <dir>
                                                 [--w4 <W4 work> <W4 work> <W4 work>] [--nproc 4]
                                                 [--extra <V2 file>=<path relative to each port build>]
                                                 [--accepted <json file of accepted deviations>]

Criteria: docs/reference/history/yeast_x5_reproduction_criteria.md (frozen) +
docs/reference/history/yeast_x5_criteria_addendum_2026-11-17.md (this script implements the addendum).
Runs in the engine environment (conda_envs/yeast_x5). Reads X5, W4 and the port builds READ-ONLY; writes only <out>/.
  - S3, L5, P2: re-scored on s0 vs X5 with the addendum rules (timing tokens, sig ignored, partner CSR as sets);
  - P3, V1, V2, V3 (must) and P4, P7, V4, X2 (report): seed-spread rule over s0, s1, s2;
  - every other row: its level from <s0>/validate/compare_x5_final/verdict.json.
Extraction helpers (leaves, num, VOLATILE, npz_dict, arrays_equal, _strip, s3_lines, tolerance constants) come from
compare_x5.py.

POST-HOC CODE FIX 2026-11-17 (Pedro approved 2026-11-17; docs/reference/history/yeast_x5_scorer_fix_2026-11-17.md),
made after the seed-spread results had been seen, implementing the frozen addendum text section 4 ("Lists (per-bin
vectors, clash.top_type_pairs[*] entries) and strings are report-only (not gating)"):
  - a scalar whose path contains a list index ([k]) is REPORT-ONLY: still scored and listed (report_only = true), never
    gating (row_level ignores it). Non-list scalars missing in a port seed stay NA, as the addendum says.
  - --extra FILE=REL: for the port seeds, V2 reads FILE (an X5-relative validation path) from <port build>/REL instead
    of <port build>/engine/work/FILE (used for the vacuole_surface.json re-run of the port defect).
  - --accepted JSON: deviations accepted by Pedro (scalar names + decision); the strict verdict is unchanged; the
    verdict with those scalars treated as passed is reported separately (verdict_after_accepted_deviations).
Margins, tolerances and every other rule are unchanged.
"""
import argparse
import glob
import json
import math
import os
import re
import sys
from contextlib import contextmanager
from multiprocessing import Pool

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_x5 as cx  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
from wholecell.config import load_site  # noqa: E402
W4_ROOT = os.path.join(str(load_site().get("yeast_root") or "<yeast_root>"), "runs")   # {yeast_root}/runs (read-only)
W4_DEFAULT = [f"{W4_ROOT}/W4_seed{i}/work" for i in range(3)]
ADDENDUM = "docs/reference/history/yeast_x5_criteria_addendum_2026-11-17.md"
SEED_ROWS = {"P3": "must", "V1": "must", "V2": "must", "V3": "must", "P4": "report", "P7": "report", "V4": "report", "X2": "report"}
SHARE_RULE = 0.99            # P4 per-tile counts and P7 per-type totals: >= 99 % of scalars pass (addendum section 4)
TIMING = re.compile(r"(?<![\w.])\d+(?:\.\d+)?s\b")
LIST_INDEX = re.compile(r"\[\d+\]")      # a path through a list element: report-only (addendum section 4; fix 2026-11-17)
SCORER_FIX = "docs/reference/history/yeast_x5_scorer_fix_2026-11-17.md"
EXTRA = {}                                  # --extra: V2 file (X5-relative) -> path relative to each port build

# diffusion benchmark test RMSE (log10), addendum table (Pedro, 2026-11-17); used in addition to the W4 files
RMSE_TABLE = {"x5": {"crowding_abs": 0.207, "crowding_ratio": 0.091, "binding_abs": 0.241, "binding_ratio": 0.163},
              "w4": [{"crowding_abs": 0.207, "crowding_ratio": 0.091, "binding_abs": 0.241, "binding_ratio": 0.163},
                     {"crowding_abs": 0.212, "crowding_ratio": 0.095, "binding_abs": 0.248, "binding_ratio": 0.170},
                     {"crowding_abs": 0.212, "crowding_ratio": 0.090, "binding_abs": 0.254, "binding_ratio": 0.187}],
              "port_s0": {"crowding_abs": 0.218, "crowding_ratio": 0.093, "binding_abs": 0.263, "binding_ratio": 0.196}}
RMSE_KEYS = {"crowding_abs": ("+ crowding(R_h) from BD", "rmse_log10_abs"), "crowding_ratio": ("+ crowding(R_h) from BD", "rmse_log10_ratio"),
             "binding_abs": ("+ crowding + binding", "rmse_log10_abs"), "binding_ratio": ("+ crowding + binding", "rmse_log10_ratio")}


# ---------------------------------------------------------------- tolerance widths (absolute, at scale n)
def tol_count(rel):
    """count tolerance max(3 sqrt n, rel n) (P3, P5-like), as a function of n and the values"""
    return lambda n, vals: max(3.0 * math.sqrt(n), rel * n)


def tol_rel(rel):
    """relative tolerance rel * n"""
    return lambda n, vals: rel * n


def tol_abs(c):
    """absolute tolerance c"""
    return lambda n, vals: c


def tol_v2(n, vals):
    """V2 scalar tolerance max(1 % relative, 0.01 absolute)"""
    return max(cx.V2_REL * n, cx.V2_ABS)


def tol_v1(path):
    """V1 leaf-class tolerance (same classes as compare_x5._v_tol_v1), as a width at scale n"""
    p = path.lower()

    def f(n, vals):
        if all(isinstance(v, int) for v in vals):
            return max(3.0, 0.01 * n)
        if "clash" in p and "rate" in p or "outside" in p:
            return max(cx.V1_RATE_REL * n, cx.V1_RATE_ABS)
        if "median_angle" in p:
            return cx.V1_ANGLE
        if "density" in p:
            return cx.V1_DENS
        if all(0 <= v <= 1 for v in vals):
            return cx.V1_FRAC
        return max(cx.V2_REL * n, cx.V2_ABS)
    return f


# ---------------------------------------------------------------- the seed-spread test for one scalar
def finite(v):
    """True for a finite int/float that is not a bool"""
    return cx.num(v) and math.isfinite(v)


def score_scalar(name, port, x5, w4, tolfn):
    """seed-spread test of one scalar: port = [s0, s1, s2], x5 = value or None, w4 = list of values/None"""
    rec = dict(name=str(name), port=port, x5=x5, w4=w4, report_only=bool(LIST_INDEX.search(str(name))))
    if not all(finite(v) for v in port):
        return dict(rec, status="NA", why="metric missing in a port seed", passed=False)
    w4v = [v for v in (w4 or []) if finite(v)]
    x5v = x5 if finite(x5) else None
    if x5v is None and not w4v:
        return dict(rec, status="NA", why="no X5 or W4 value", passed=False)
    lo_p, hi_p = min(port), max(port)
    n = max([abs(lo_p), abs(hi_p)] + ([abs(x5v)] if x5v is not None else []))
    tol = tolfn(n, port + ([x5v] if x5v is not None else []))
    margin = max(0.5 * (hi_p - lo_p), tol)
    lo, hi = lo_p - margin, hi_p + margin
    a = all(lo <= v <= hi for v in port)
    x5_ok = x5v is not None and lo <= x5v <= hi
    w4_ok = bool(w4v) and min(w4v) <= hi and max(w4v) >= lo
    matched = "+".join(k for k, ok in (("X5", x5_ok), ("W4", w4_ok)) if ok) or "none"
    ok = a and (x5_ok or w4_ok)
    return dict(rec, status="pass" if ok else "fail", passed=ok, port_min=lo_p, port_max=hi_p, tol=tol, margin=margin, lo=lo, hi=hi,
                w4_min=min(w4v) if w4v else None, w4_max=max(w4v) if w4v else None, matched=matched,
                ref="X5+W4" if w4v else "X5 only (no W4 counterpart)", why="" if ok else ("port seed outside" if not a else "X5 and W4 outside"))


def score_dicts(prefix, port_d, x5_d, w4_d, tolfn_for, keys=None):
    """score every key of the port/X5 dicts (key union of port and X5 unless keys given); w4_d entries may be None"""
    keys = keys if keys is not None else sorted(set().union(*[set(d) for d in port_d + [x5_d]]), key=str)
    out = []
    for k in keys:
        w4 = [None if d is None else d.get(k) for d in w4_d]
        out.append(score_scalar(f"{prefix}{k}", [d.get(k) for d in port_d], x5_d.get(k), w4, tolfn_for(k)))
    return out


def row_level(scalars, share_groups=None):
    """row pass rule: every scalar passes, except groups in share_groups (prefix -> min share of passing scalars)"""
    share_groups = share_groups or {}
    n_ro = sum(1 for s in scalars if s.get("report_only"))
    scalars = [s for s in scalars if not s.get("report_only")]          # list-indexed scalars never gate (section 4)
    if not scalars:
        return "NA", f"no gating scalars ({n_ro} report-only)"
    if all(s["status"] == "NA" for s in scalars):
        return "NA", "no scalar evaluable"
    notes = []
    ok = True
    for pre, share in share_groups.items():
        grp = [s for s in scalars if s["name"].startswith(pre)]
        frac = sum(s["passed"] for s in grp) / max(len(grp), 1)
        notes.append(f"{pre}* {frac:.4f} pass (need >= {share})")
        ok &= bool(grp) and frac >= share
    rest = [s for s in scalars if not any(s["name"].startswith(p) for p in share_groups)]
    nbad = sum(not s["passed"] for s in rest)
    notes.append(f"other scalars: {len(rest) - nbad}/{len(rest)} pass")
    notes.append(f"{n_ro} list-indexed scalars report-only")
    ok &= nbad == 0
    return ("SS" if ok else "fail"), "; ".join(notes)


# ---------------------------------------------------------------- extractors (one work dir -> {name: value})
def lib_names(work):
    """library type names in the order of yc3d/lib_ids.npy"""
    return [str(x) for x in np.load(f"{work}/yc3d/lib_ids.npy", allow_pickle=True)]


def ex_p3(work):
    """P3 scalars: requested copies per (type name, region) from tile_counts.parquet; thinning.json numbers"""
    names = lib_names(work)
    T = pd.read_parquet(f"{work}/yc3d/tile_counts.parquet", columns=["lib", "region", "n"])
    g = T.groupby(["lib", "region"]).n.sum()
    pairs = {f"{names[int(l)]}|{int(r)}": int(n) for (l, r), n in g.items()}
    thin = {k: v for k, v in cx.leaves(json.load(open(f"{work}/yc3d/thinning.json"))).items() if finite(v)}
    return dict(pairs=pairs, thin=thin, names=set(names))


def numeric_leaves(path):
    """non-volatile, non-bool numeric leaves of a json file (lists, strings and flags are not scored)"""
    L = cx.leaves(json.load(open(path)))
    return {k: v for k, v in L.items() if finite(v) and not cx.VOLATILE.search(k.split(".")[-1])}


def ex_v1(work):
    """V1 scalars: numeric leaves of results/wholecell3d_validation.json"""
    return numeric_leaves(f"{work}/results/wholecell3d_validation.json")


def ex_v2(work, files):
    """V2 scalars: numeric leaves per validation json (None for a missing file)"""
    out = {}
    for rel in files:
        p = os.path.join(work, rel)
        try:
            out[rel] = numeric_leaves(p) if os.path.exists(p) else None
        except Exception:  # noqa: BLE001
            out[rel] = None
    return out


def ex_v3(work):
    """V3 scalars: contact pair count, encounter rows/means, sweep phi, diffusion-benchmark test RMSE (4)"""
    out = {}
    with np.load(f"{work}/contacts/contacts.npz") as z:
        out["contacts.pairs"] = sum(int(z[k].sum()) for k in z.files if re.fullmatch(r"C\d+", k))
    e = pd.read_parquet(f"{work}/contacts/encounters_kinases.parquet")
    out["encounters.rows"] = len(e)
    for c in ("p_near_2", "p_near_10", "c_local_uM"):
        if c in e.columns:
            out[f"encounters.mean_{c}"] = float(e[c].mean())
    for j in ("cyto120_mf", "gem200_mf", "nuc100_mf"):
        out[f"sweeps.{j}.phi"] = json.load(open(f"{work}/sweeps/{j}.json")).get("phi")
    d = json.load(open(f"{work}/validation/diffusion_benchmark.json"))
    for k, (model, key) in RMSE_KEYS.items():
        out[f"diffusion.test_{k}"] = d.get(model, {}).get("test", {}).get(key)
    return out


def ex_v4(work):
    """V4 scalars: numbers of the condensate log-phase / stress lines (by position) and microscopy spot counts"""
    out = {}
    lines = [cx._strip(l) for l in open(f"{work}/condensates/cond_test.log") if re.match(r"(log phase|stress) ", l)]
    for i, l in enumerate(lines):
        for j, v in enumerate(cx._nums(l)):
            out[f"cond.line{i}.n{j}"] = v
    for k, v in json.load(open(f"{work}/microscopy/spot_counts.json")).items():
        if finite(v):
            out[f"spots.{k}"] = v
    return out


def ex_p4(work):
    """P4 scalars: placed instances per tile (pack_summary.parquet) and their total"""
    S = pd.read_parquet(f"{work}/yc3d/pack_summary.parquet", columns=["ti", "placed"])
    out = {f"tile.{int(t)}": int(p) for t, p in zip(S.ti, S.placed)}
    out["total_placed"] = int(S.placed.sum())
    return out


def _tile_jobs(path):
    """per-library placed counts of one tiles_raw tile (job table)"""
    with np.load(path) as z:
        return pd.Series(z["job_placed"].astype(np.int64)).groupby(z["job_lib"].astype(np.int64)).sum()


def ex_p7(work, nproc):
    """P7 scalars: number of tiles_raw files, total placed (job tables), placed per library type name"""
    fs = sorted(glob.glob(f"{work}/yc3d/tiles_raw/*.npz"))
    with Pool(nproc) as pool:
        r = pool.map(_tile_jobs, fs, chunksize=32)
    g = pd.concat(r).groupby(level=0).sum() if r else pd.Series(dtype=np.int64)
    names = lib_names(work)
    out = {f"type.{names[int(l)]}": int(n) for l, n in g.items()}
    out["n_tile_files"] = len(fs)
    out["total_placed"] = int(g.sum())
    return out


def ex_x2(work):
    """X2 scalars: number of viewer tile files and of viewer library files"""
    d = f"{work}/yc3d/web/data"
    return {"n_tile_files": len(glob.glob(f"{d}/tiles/*")), "n_library_files": len(glob.glob(f"{d}/library/*"))}


def md5(path):
    """md5 of a file"""
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


# ---------------------------------------------------------------- extraction over all builds (memoised per path)
CACHE = {}


def get(kind, work, *args):
    """memoised extractor call; returns None (= metric missing) on any error"""
    key = (kind, work)
    if key not in CACHE:
        fn = {"p3": ex_p3, "v1": ex_v1, "v2": ex_v2, "v3": ex_v3, "v4": ex_v4, "p4": ex_p4, "p7": ex_p7, "x2": ex_x2}[kind]
        try:
            CACHE[key] = fn(work, *args)
        except Exception as e:  # noqa: BLE001
            print(f"   extraction {kind} failed for {work}: {type(e).__name__}: {e}", flush=True)
            CACHE[key] = None
    return CACHE[key]


def collect(kind, works, *args):
    """extract one metric family for the port seeds, X5 and the W4 seeds: -> (ports, x5, w4s)"""
    vals = [get(kind, w, *args) for w in works["port"]]
    return vals, get(kind, works["x5"], *args), [get(kind, w, *args) for w in works["w4"]]


def dict_or_empty(d):
    """None -> {} so that a missing build gives NA scalars instead of a crash"""
    return d if d is not None else {}


# ---------------------------------------------------------------- seed-spread rows
def row_p3(works):
    """P3: requested per (type, region) with max(3 sqrt n, 0.5 % n); thinning factors +-0.005; all scalars pass"""
    ports, x5, w4 = collect("p3", works)
    pp = [dict_or_empty(p).get("pairs", {}) for p in ports]
    xp = dict_or_empty(x5).get("pairs", {})
    keys = sorted(set().union(*[set(d) for d in pp + [xp]]))
    w4p = []
    for w in w4:
        if w is None:
            w4p.append(None)
            continue
        # a pair whose type is not in the W4 library has no W4 counterpart (None); a missing pair of a known type is 0
        w4p.append({k: w["pairs"].get(k, 0) if k.split("|")[0] in w["names"] else None for k in keys})
    fill = [{k: d.get(k, 0) for k in keys} if p is not None else {} for d, p in zip(pp, ports)]
    xf = {k: xp.get(k, 0) for k in keys} if x5 is not None else {}
    sc = score_dicts("pair.", fill, xf, w4p, lambda k: tol_count(cx.P3_REL), keys)
    th = [dict_or_empty(p).get("thin", {}) for p in ports]
    sc += score_dicts("thin.", th, dict_or_empty(x5).get("thin", {}), [None if w is None else w["thin"] for w in w4],
                      lambda k: tol_abs(cx.P3_THIN))
    return sc, {}


def row_v1(works):
    """V1: numeric leaves of wholecell3d_validation.json with the V1 leaf classes; all scalars pass"""
    ports, x5, w4 = collect("v1", works)
    return score_dicts("", [dict_or_empty(p) for p in ports], dict_or_empty(x5), w4, tol_v1), {}


def v2_files(x5work):
    """the V2 file list: X5 validation/*.json + validation/latest/*.json minus MANIFEST and diffusion_benchmark"""
    fs = sorted(glob.glob(f"{x5work}/validation/*.json")) + sorted(glob.glob(f"{x5work}/validation/latest/*.json"))
    skip = {"MANIFEST.json"} | cx.GPU_FILES
    return [os.path.relpath(f, x5work) for f in fs if os.path.basename(f) not in skip]


def row_v2(works):
    """V2: per file, numeric leaves within max(1 %, 0.01); a file passes if all its scalars pass; row = all files"""
    files = v2_files(works["x5"])
    ports, x5, w4 = collect("v2", works, files)
    extra_used = {}
    if EXTRA:
        ports = [dict(p) if p is not None else None for p in ports]
        for i, w in enumerate(works["port"]):
            b = os.path.dirname(os.path.dirname(w))                      # <build>/engine/work -> <build>
            for f, rel in EXTRA.items():
                if f not in files or ports[i] is None:
                    raise SystemExit(f"--extra {f}: not a V2 file of X5 or no V2 extraction for {w}")
                q = os.path.join(b, rel)
                ports[i][f] = numeric_leaves(q) if os.path.exists(q) else None
                extra_used.setdefault(f, []).append(dict(seed=i, path=q, exists=os.path.exists(q), md5=md5(q) if os.path.exists(q) else None))
    sc = []
    for f in files:
        pd_ = [dict_or_empty(dict_or_empty(p).get(f)) for p in ports]
        w4f = [None if w is None else w.get(f) for w in w4]
        xf = dict_or_empty(dict_or_empty(x5).get(f))
        if not any(pd_) and not xf:
            sc.append(score_scalar(f"{f}:<file>", [None] * 3, None, None, tol_v2))
            continue
        sc += score_dicts(f"{f}:", pd_, xf, w4f, lambda k: tol_v2)
    per_file = {}
    for s in sc:
        f = s["name"].split(":")[0]
        per_file.setdefault(f, [0, 0])
        per_file[f][0] += 1
        per_file[f][1] += int(not s["passed"])
    gating = {}
    for s in sc:
        if not s.get("report_only"):
            f = s["name"].split(":")[0]
            gating.setdefault(f, 0)
            gating[f] += int(not s["passed"])
    return sc, {"V2_files_failing": sorted(f for f, b in gating.items() if b), "V2_files": per_file,
                "V2_files_failing_incl_report_only": sorted(f for f, (n, b) in per_file.items() if b), "V2_extra_inputs": extra_used}


def row_v3(works):
    """V3: contacts 1 %, encounters 5 %, sweep phi 2 x SD / nuc100 0.01, diffusion test RMSE 10 %; all scalars pass"""
    ports, x5, w4s = collect("v3", works)
    tolmap = {"contacts.": tol_rel(cx.V3_CONTACT_REL), "encounters.": tol_rel(cx.V3_ENC_REL),
              "sweeps.cyto120_mf": tol_abs(2 * cx.V3_PHI_SD["cyto120_mf"]), "sweeps.gem200_mf": tol_abs(2 * cx.V3_PHI_SD["gem200_mf"]),
              "sweeps.nuc100_mf": tol_abs(cx.V3_NUC_PHI), "diffusion.": tol_rel(cx.V3_DIFF_REL)}
    tolfor = lambda k: next(v for p, v in tolmap.items() if k.startswith(p))  # noqa: E731
    pd_ = [dict_or_empty(p) for p in ports]
    xd = dict(dict_or_empty(x5))
    keys = sorted(set().union(*[set(d) for d in pd_ + [xd]]))
    sc = []
    for k in keys:
        w4 = [None if w is None else w.get(k) for w in w4s]
        if k.startswith("diffusion.test_"):
            t = k[len("diffusion.test_"):]
            w4 = w4 + [r[t] for r in RMSE_TABLE["w4"]]          # addendum table, in addition to the W4 files
            if not finite(xd.get(k)):
                xd[k] = RMSE_TABLE["x5"][t]
        sc.append(score_scalar(k, [d.get(k) for d in pd_], xd.get(k), w4, tolfor(k)))
    return sc, {"V3_rmse_table": RMSE_TABLE}


def row_v4(works):
    """V4: condensate line numbers max(1 %, 0.05); spot counts max(1, 1 %); all scalars pass"""
    ports, x5, w4 = collect("v4", works)
    tolfor = lambda k: (lambda n, v: max(1.0, cx.V4_REL * n)) if k.startswith("spots.") else (lambda n, v: max(cx.V4_REL * n, cx.V4_ONSET))  # noqa: E731
    return score_dicts("", [dict_or_empty(p) for p in ports], dict_or_empty(x5), w4, tolfor), {}


def row_p4(works):
    """P4: per-tile placed counts max(3 sqrt n, 1 % n) for >= 99 % of tiles; total placed within 0.5 % (required)"""
    ports, x5, w4 = collect("p4", works)
    tolfor = lambda k: tol_rel(cx.P6_INST_REL) if k == "total_placed" else tol_count(cx.P5_REL)  # noqa: E731
    return score_dicts("", [dict_or_empty(p) for p in ports], dict_or_empty(x5), w4, tolfor), {"share_groups": {"tile.": SHARE_RULE}}


def row_p7(works, nproc):
    """P7: tiles_raw file count (exact), total placed 0.5 %, per-type placed max(3 sqrt n, 1 % n) for >= 99 % of types"""
    ports, x5, w4 = collect("p7", works, nproc)
    tolfor = lambda k: tol_abs(0) if k == "n_tile_files" else (tol_rel(cx.P6_INST_REL) if k == "total_placed" else tol_count(cx.P5_REL))  # noqa: E731
    # a type not placed in a W4 seed has no W4 counterpart (None); a type absent from a port/X5 build counts 0
    pd_ = [dict_or_empty(p) for p in ports]
    xd = dict_or_empty(x5)
    keys = sorted(set().union(*[set(d) for d in pd_ + [xd]]))
    fill = lambda d, ok: {k: d.get(k, 0) for k in keys} if ok else {}  # noqa: E731
    return score_dicts("", [fill(d, p is not None) for d, p in zip(pd_, ports)], fill(xd, x5 is not None), w4, tolfor, keys), \
        {"share_groups": {"type.": SHARE_RULE}}


def row_x2(works):
    """X2: viewer tile-file count and library-file count, exact; both pass"""
    ports, x5, w4 = collect("x2", works)
    return score_dicts("", [dict_or_empty(p) for p in ports], dict_or_empty(x5), w4, lambda k: tol_abs(0)), {}


# ---------------------------------------------------------------- S3, L5, P2 (addendum rules, s0 vs X5)
@contextmanager
def no_timing_tokens():
    """temporarily make compare_x5._nums drop timing tokens (a number followed by 's')"""
    orig = cx._nums
    cx._nums = lambda s: orig(TIMING.sub(" ", s))
    try:
        yield
    finally:
        cx._nums = orig


def rescore_s3(R, P):
    """S3 with timing tokens excluded (compare_x5.s3_lines, numbers without e.g. '3s')"""
    with no_timing_tokens():
        cx.s3_lines(R, P)
    r = cx.ROWS.pop()
    r["note"] = (r["note"] + "; " if r["note"] else "") + "timing tokens excluded (addendum 2026-11-17)"
    return r, {"S3_lines": cx.DETAILS.get("S3_lines")}


def npz_tol(A, B, tol, rel, skip=()):
    """compare two npz dicts ignoring keys in skip: 'E', 'T' or 'fail (...)' with float |d| <= tol (abs) or rel"""
    ks = [k for k in A if k not in skip]
    if sorted(ks) != sorted(k for k in B if k not in skip):
        return "fail (keys)", math.inf
    if all(cx.arrays_equal(A[k], B[k]) for k in ks):
        return "E", 0.0
    md = 0.0
    for k in ks:
        if cx.arrays_equal(A[k], B[k]):
            continue
        if A[k].shape != B[k].shape or A[k].dtype.kind not in "fc":
            return f"fail ({k} non-float differs)", math.inf
        d = np.abs(A[k].astype(float) - B[k].astype(float))
        md = max(md, float(np.max(d / np.maximum(np.abs(A[k]), 1e-6) if rel else d)))
    return ("T" if md <= tol else f"fail ({md:.3g})"), md


def rescore_l5(R, P):
    """L5: lib_ids identical; libpack_meta / pack_trees_15 centres and radii within 0.05 nm; 'sig' ignored"""
    ia = np.load(f"{R}/yc3d/lib_ids.npy", allow_pickle=True).astype(str)
    ib = np.load(f"{P}/yc3d/lib_ids.npy", allow_pickle=True).astype(str)
    ids_ok = cx.arrays_equal(ia, ib)
    res = {f: npz_tol(cx.npz_dict(f"{R}/yc3d/{f}", True), cx.npz_dict(f"{P}/yc3d/{f}", True), cx.L2_TREE_NM, False, skip=("sig",))[0]
           for f in ("libpack_meta.npz", "pack_trees_15.npz")}
    lev = level_of(ids_ok, res.values())
    return dict(id="L5", metric="lib_ids.npy identical; libpack_meta / pack_trees_15 (sig ignored)", x5=len(ia), port=len(ib), delta=cx._s(res),
                tolerance=f"ids identical; radii/centres <= {cx.L2_TREE_NM} nm", level=lev, cls="must",
                note="signature strings ignored (addendum 2026-11-17)"), {"L5": res}


def level_of(ok, levels):
    """combine part levels: E if all E, T if all E/T, else fail"""
    levels = list(levels)
    if ok and all(v == "E" for v in levels):
        return "E"
    return "T" if ok and all(v in ("E", "T") for v in levels) else "fail"


def csr_sets(Z):
    """partner CSR rows sorted by (row, idx, kind, strength): tie order inside a row is not content"""
    off = Z["off"]
    rid = np.repeat(np.arange(len(off) - 1), np.diff(off))
    o = np.lexsort((Z["strength"], Z["kind"], Z["idx"], rid))
    return Z["idx"][o], Z["kind"][o], Z["strength"][o]


def rescore_p2(R, P):
    """P2: network.npz floats rel <= 1e-5; partners.npz compared as per-row sets of (idx, kind) with strength rel <= 1e-5"""
    res = {}
    A, B = cx.npz_dict(f"{R}/organisation/network.npz"), cx.npz_dict(f"{P}/organisation/network.npz")
    res["network.npz"] = npz_tol(A, B, cx.P2_REL, True)[0]
    A, B = cx.npz_dict(f"{R}/organisation/partners.npz"), cx.npz_dict(f"{P}/organisation/partners.npz")
    info = {}
    if sorted(A) != sorted(B) or not cx.arrays_equal(A["off"], B["off"]) or not cx.arrays_equal(A["ids"].astype(str), B["ids"].astype(str)):
        res["partners.npz"] = "fail (ids/off differ)"
    elif all(cx.arrays_equal(A[k], B[k]) for k in A):
        res["partners.npz"] = "E"
    else:
        (ia, ka, sa), (ib, kb, sb) = csr_sets(A), csr_sets(B)
        sets_ok = np.array_equal(ia, ib) and np.array_equal(ka, kb)
        md = float(np.max(np.abs(sa - sb) / np.maximum(np.abs(sa), 1e-6))) if sets_ok else math.inf
        info = dict(sets_identical=bool(sets_ok), max_rel_strength=md,
                    positions_differing_raw=int((A["idx"] != B["idx"]).sum()))
        res["partners.npz"] = "T" if sets_ok and md <= cx.P2_REL else f"fail ({md:.3g})"
    lev = level_of(True, res.values())
    return dict(id="P2", metric="organisation network.npz / partners.npz (partner CSR as per-row sets)", x5="", port="", delta=cx._s(res),
                tolerance=f"partner sets identical; floats rel <= {cx.P2_REL}", level=lev, cls="must",
                note=f"tie order ignored (addendum 2026-11-17); {info}"), {"P2": res, "P2_partners": info}


# ---------------------------------------------------------------- assembly
def seed_row(rid, sc, extra, metric):
    """turn the scalar list of one seed-spread row into a result row + details"""
    lev, rule_note = row_level(sc, extra.get("share_groups"))
    g = [s for s in sc if not s.get("report_only")]
    ro = [s for s in sc if s.get("report_only")]
    bad = [s for s in g if not s["passed"]]
    ro_bad = [s for s in ro if not s["passed"]]
    ok_sc = [s for s in g if s["passed"]]
    no_w4 = sum(1 for s in g if s.get("ref", "").startswith("X5 only"))
    show = bad[:5] if bad else g[:12]
    r = dict(id=rid, metric=metric, x5="", port=f"{len(g)} gating scalars (+{len(ro)} report-only)",
             delta=f"{len(bad)} not passed ({sum(s['status'] == 'NA' for s in bad)} NA); report-only not passed {len(ro_bad)}",
             tolerance="seed spread: margin = max(0.5 x port range, criteria tolerance)", level=lev, cls=SEED_ROWS[rid],
             note=f"{rule_note}; {no_w4} scalars tested against X5 only (no W4 counterpart)")
    summ = dict(n_scalars=len(g), n_pass=len(ok_sc), n_fail=sum(s["status"] == "fail" for s in g), n_na=sum(s["status"] == "NA" for s in g),
                n_x5_only=no_w4, rule=rule_note, shown_scalars=[slim(s) for s in show],
                matched={m: sum(1 for s in g if s.get("matched") == m) for m in ("X5", "W4", "X5+W4", "none")},
                report_only=dict(n_scalars=len(ro), n_pass=sum(s["passed"] for s in ro), n_fail=sum(s["status"] == "fail" for s in ro),
                                 n_na=sum(s["status"] == "NA" for s in ro),
                                 list_paths=sorted({LIST_INDEX.sub("[*]", s["name"]) for s in ro})))
    det = dict(summary=summ, not_passed=[slim(s) for s in bad[:2000]], report_only_not_passed=[slim(s) for s in ro_bad[:2000]],
               **{k: v for k, v in extra.items() if k != "share_groups"})
    return r, summ, det


def slim(s):
    """compact json-able form of a scalar record"""
    keys = ("name", "status", "report_only", "port", "port_min", "port_max", "margin", "tol", "lo", "hi", "x5", "w4_min", "w4_max", "matched", "ref", "why")
    return {k: s.get(k) for k in keys if k in s}


def verdict_of(rows):
    """section 7 verdict of a list of result rows"""
    passed = lambda r: r["level"] in ("E", "T", "SS")  # noqa: E731
    must = [r for r in rows if r["cls"] == "must"]
    rep = [r for r in rows if r["cls"] == "report"]
    return "NOT REPRODUCED" if not all(passed(r) for r in must) else ("REPRODUCED WITH DEVIATIONS" if not all(passed(r) for r in rep) else "REPRODUCED")


def after_accepted(rows, scalars, share, accepted):
    """rows re-levelled with the accepted scalars treated as passed (scoring itself unchanged); -> (verdict, block)"""
    names = {n for a in accepted for n in a["scalars"]}
    found, rows2, changed, remaining = set(), [], {}, {}
    for r in rows:
        sc = scalars.get(r["id"])
        if sc is None:
            rows2.append(r)
            continue
        sc2 = []
        for s in sc:
            if s["name"] in names:
                found.add(s["name"])
                sc2.append(dict(s, passed=True, status="accepted"))
            else:
                sc2.append(s)
        lev, note = row_level(sc2, share.get(r["id"]))
        if lev != r["level"]:
            changed[r["id"]] = f"{r['level']} -> {lev}"
        rows2.append(dict(r, level=lev))
        if any(s["name"] in names for s in sc):
            remaining[r["id"]] = [slim(s) for s in sc2 if not s["passed"] and not s.get("report_only")]
    passed = lambda r: r["level"] in ("E", "T", "SS")  # noqa: E731
    v = verdict_of(rows2)
    return v, dict(verdict=v, rows_changed=changed, remaining_not_passed_gating=remaining, missing_accepted_scalars=sorted(names - found),
                   failed_must=[r["id"] for r in rows2 if r["cls"] == "must" and not passed(r)],
                   failed_report=[r["id"] for r in rows2 if r["cls"] == "report" and not passed(r)],
                   counts={c: {l: sum(1 for r in rows2 if r["cls"] == c and r["level"] == l) for l in ("E", "T", "SS", "fail", "NA")}
                           for c in ("must", "report")})


def write_outputs(out, rows, details, scalars, args, share=None, accepted=None):
    """verdict.json, criteria_results.tsv, details.json, scalars.tsv"""
    passed = lambda r: r["level"] in ("E", "T", "SS")  # noqa: E731
    must = [r for r in rows if r["cls"] == "must"]
    rep = [r for r in rows if r["cls"] == "report"]
    v = verdict_of(rows)
    res = dict(verdict=v, verdict_kind="strict (frozen criteria + addendum, scorer fix 2026-11-17)",
               criteria="docs/reference/history/yeast_x5_reproduction_criteria.md", addendum=ADDENDUM, scorer_fix=SCORER_FIX,
               builds=dict(s0=args.s0, s1=args.s1, s2=args.s2, x5=args.x5, w4=args.w4), extra_inputs=EXTRA,
               counts={c: {l: sum(1 for r in rows if r["cls"] == c and r["level"] == l) for l in ("E", "T", "SS", "fail", "NA")} for c in ("must", "report")},
               failed_must=[r["id"] for r in must if not passed(r)], failed_report=[r["id"] for r in rep if not passed(r)])
    if accepted:
        v2, blk = after_accepted(rows, scalars, share or {}, accepted)
        res["accepted_deviations"] = accepted
        res["verdict_after_accepted_deviations"] = v2
        res["after_accepted_deviations"] = dict(blk, note="the accepted scalars are treated as passed; nothing else changes; the strict "
                                                          "verdict above stays the verdict of the frozen criteria + addendum")
    res["rows"] = [dict(r, passed=passed(r), seed_spread=details.get(r["id"], {}).get("summary")) for r in rows]
    os.makedirs(out, exist_ok=True)
    cols = ["id", "metric", "x5", "port", "delta", "tolerance", "level", "cls", "source", "note"]
    pd.DataFrame(rows).reindex(columns=cols).to_csv(os.path.join(out, "criteria_results.tsv"), sep="\t", index=False)
    json.dump(res, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    json.dump(details, open(os.path.join(out, "details.json"), "w"), indent=1, default=str)
    pd.DataFrame([dict(row=rid, **slim(s)) for rid, sc in scalars.items() for s in sc]).to_csv(os.path.join(out, "scalars.tsv.gz"), sep="\t", index=False)
    print("VERDICT:", v, json.dumps(res["counts"]), "failed must:", res["failed_must"], "failed report:", res["failed_report"])
    if accepted:
        print("VERDICT AFTER ACCEPTED DEVIATIONS:", res["verdict_after_accepted_deviations"], json.dumps(res["after_accepted_deviations"]["counts"]),
              "failed must:", res["after_accepted_deviations"]["failed_must"], "failed report:", res["after_accepted_deviations"]["failed_report"])
    return res


def main():
    """parse arguments, rescore S3/L5/P2, score the seed-spread rows, copy the rest from s0, write the verdict"""
    ap = argparse.ArgumentParser()
    for s in ("s0", "s1", "s2"):
        ap.add_argument(f"--{s}", required=True, help="port build dir")
    ap.add_argument("--x5", required=True, help="X5 run dir (contains work/)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--w4", nargs=3, default=W4_DEFAULT, help="three W4 work dirs")
    ap.add_argument("--nproc", type=int, default=4)
    ap.add_argument("--extra", action="append", default=[], metavar="FILE=REL",
                    help="V2: read FILE (X5-relative, e.g. validation/vacuole_surface.json) from <port build>/REL for every port seed")
    ap.add_argument("--accepted", default=None, help="json list of accepted deviations [{scalars: [...], decision, note}] (verdict after them "
                                                    "reported separately; strict verdict unchanged)")
    a = ap.parse_args()
    for x in a.extra:
        f, _, rel = x.partition("=")
        if not rel:
            raise SystemExit(f"--extra {x}: expected FILE=REL")
        EXTRA[f] = rel
    accepted = json.load(open(a.accepted)) if a.accepted else None
    R = os.path.join(a.x5, "work")
    pw = [os.path.join(b, "engine", "work") for b in (a.s0, a.s1, a.s2)]
    works = dict(port=pw, x5=R, w4=a.w4)
    base = json.load(open(os.path.join(a.s0, "validate", "compare_x5_final", "verdict.json")))["rows"]
    new, details, scalars, share = {}, {}, {}, {}
    for rid, fn in (("S3", rescore_s3), ("L5", rescore_l5), ("P2", rescore_p2)):
        try:
            r, d = fn(R, pw[0])
        except Exception as e:  # noqa: BLE001
            r, d = dict(id=rid, metric=rid, level="NA", cls="must", note=f"error {type(e).__name__}: {e}"), {}
        new[rid] = dict(r, source="addendum rule, s0 vs X5")
        details[rid] = d
        print(f"{rid:4s} {new[rid]['level']:4s} {new[rid].get('delta', '')}", flush=True)
    seeders = {"P3": (row_p3, "requested per (type, region); thinning"), "V1": (row_v1, "stage-4 validation scalars"),
               "V2": (row_v2, "validate_all CPU comparison scalars"), "V3": (row_v3, "GPU steps: contacts, encounters, sweeps, diffusion test RMSE"),
               "P4": (row_p4, "placed per tile; total placed"), "P7": (lambda w: row_p7(w, a.nproc), "tiles_raw: file count, total, per-type placed"),
               "V4": (row_v4, "condensate lines; spot counts"), "X2": (row_x2, "viewer tile / library file counts")}
    for rid, (fn, metric) in seeders.items():
        try:
            sc, extra = fn(works)
            share[rid] = extra.get("share_groups")
            r, summ, det = seed_row(rid, sc, extra, metric)
        except Exception as e:  # noqa: BLE001
            import traceback
            sc, r, det = [], dict(id=rid, metric=metric, level="NA", cls=SEED_ROWS[rid], note=f"error {type(e).__name__}: {e}"), {"error": traceback.format_exc()}
        new[rid] = dict(r, source="seed spread (addendum)")
        details[rid] = det
        scalars[rid] = sc
        print(f"{rid:4s} {r['level']:4s} {r.get('delta', '')} | {r.get('note', '')}", flush=True)
    rows = [new.get(r["id"], dict(r, source="s0 compare_x5_final")) for r in base]
    write_outputs(a.out, rows, details, scalars, a, share, accepted)


if __name__ == "__main__":
    main()
