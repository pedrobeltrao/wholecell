"""Compare a wholecell yeast build with yeast build X5 (criteria: docs/reference/history/yeast_x5_reproduction_criteria.md, frozen at 941ee95).

    python plugins/yeast/repro/compare_x5.py --bdir <build> --ref <yeast_root>/runs/X5_seed0 [--nproc 16]
    python plugins/yeast/repro/compare_x5.py --self-test --ref <X5>       # X5 against itself: every row must be E

Runs in the engine environment. Reads X5 READ-ONLY. Writes <build>/validate/compare_x5/{criteria_results.tsv,
verdict.json, details/*.json}. Level E = identical content; T = within the tolerances below (constants = the criteria
tables); fail otherwise. Interpretation notes (fixed at commit, before the run):
  - npz / json / parquet are compared by content (zip timestamps, gzip headers, build ids, absolute paths, timings and
    log timestamps are not content);
  - V1/V2 "every scalar": numeric values in (nested) dicts; list-valued entries (per-bin / per-protein vectors) are
    summarised (share of elements outside the tolerance) in the details, not gating;
  - V1 leaf classes: integers = counts (max(3, 1 %)); keys with clash/outside = rates (max(25 % rel, 1e-5));
    median_angle_deg (0.01 deg); keys under 'density' (2 mg/ml); other values in [0, 1] = fractions (0.005);
    other floats 1 % relative or 0.01 absolute (the V2 default).
"""
import argparse
import glob
import gzip
import hashlib
import json
import math
import os
import re
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- tolerances (= criteria tables)
G1_TOL_NM = 0.01
G2_TOL_NM = 0.01
G4_TOT_REL = 1e-3
G4_TILE_FRAC_VOL = 0.005 * 200.0 ** 3
G4_TILE_FRAC_AREA = 0.005 * 200.0 ** 2
G4_TILE_SHARE = 0.99
G5_AREA_REL = 1e-3
I1_LOG2 = 0.01
I1_W = 1e-6
I2_COPIES_REL = 1e-9
I2_OTHER_REL = 1e-6
I3_REL = 1e-12
I4_REL = 1e-3
A1_REL = 1e-9
L2_BEADS_NM = 1e-4
L2_TREE_NM = 0.05
L2_SHARE = 0.999
S2_ABS, S2_REL = 3, 0.01
S3_REL = 0.01
S4_CENTROID_NM, S4_RMS_REL = 20.0, 0.01
P1_REL = 1e-6
P2_REL = 1e-5
P3_REL = 0.005
P3_THIN = 0.005
P5_REL = 0.01
P5_SMALL_ABS, P5_SMALL_SHARE = 3, 0.99
P6_FRAC, P6_INST_REL = 0.003, 0.005
V1_RATE_REL, V1_RATE_ABS = 0.25, 1e-5
V1_ANGLE = 0.01
V1_DENS = 2.0
V1_FRAC = 0.005
V2_REL, V2_ABS = 0.01, 0.01
V3_CONTACT_REL = 0.01
V3_ENC_REL = 0.05
V3_PHI_SD = {"cyto120_mf": 0.013, "gem200_mf": 0.010}
V3_NUC_PHI = 0.01
V3_DIFF_REL = 0.10
V4_REL, V4_ONSET = 0.01, 0.05

ROWS = []
DETAILS = {}


def row(rid, metric, x5, port, delta, tol, level, cls, note=""):
    ROWS.append(dict(id=rid, metric=metric, x5=_s(x5), port=_s(port), delta=_s(delta), tolerance=tol, level=level, cls=cls, note=note))
    print(f"{rid:5s} {level:5s} {cls:6s} {metric} | x5={_s(x5)} port={_s(port)} d={_s(delta)} {note}", flush=True)


def _s(v):
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    s = str(v)
    return s if len(s) < 200 else s[:197] + "..."


def guard(rid, metric, cls):
    """decorator: a crash in one comparison is a 'fail' row with the exception, never a silent skip"""
    def deco(fn):
        def run(*a, **k):
            try:
                return fn(*a, **k)
            except Exception as e:  # noqa: BLE001
                import traceback
                row(rid, metric, None, None, None, "-", "fail", cls, f"comparison error: {type(e).__name__}: {e}")
                DETAILS[f"{rid}_error"] = traceback.format_exc()
        return run
    return deco


def npz_dict(p, pickle=False):
    with np.load(p, allow_pickle=pickle) as z:
        return {k: z[k] for k in z.files}


def arrays_equal(a, b):
    if a.dtype.kind in "fc" and b.dtype.kind in "fc":
        return a.shape == b.shape and np.array_equal(a, b, equal_nan=True)
    return a.shape == b.shape and np.array_equal(a, b)


def npz_equal(pa, pb, pickle=False):
    A, Bz = npz_dict(pa, pickle), npz_dict(pb, pickle)
    return sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A)


def leaves(o, pre=""):
    """flatten json -> {path: scalar}, lists of scalars kept as lists"""
    out = {}
    if isinstance(o, dict):
        for k, v in o.items():
            out.update(leaves(v, f"{pre}.{k}" if pre else str(k)))
    elif isinstance(o, list) and o and all(isinstance(x, (dict, list)) for x in o):
        for i, v in enumerate(o):
            out.update(leaves(v, f"{pre}[{i}]"))
    else:
        out[pre] = o
    return out


def num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


VOLATILE = re.compile(r"(build_id|build|checked|modified|started|t0|host|time|seconds|wall_s|elapsed|path|file|dir|commit|date|version)$", re.I)


# ---------------------------------------------------------------- geometry
def geometry(R, P):
    @guard("G1", "geometry.json (G, pores, primitives)", "must")
    def g1():
        A = json.load(open(f"{R}/yc3d/geometry.json")); Bj = json.load(open(f"{P}/yc3d/geometry.json"))
        if A == Bj:
            return row("G1", "geometry.json content", "identical", "identical", 0.0, "E", "E", "must")
        same = sorted(A) == sorted(Bj) and len(A["PT"]) == len(Bj["PT"]) and len(A["pores"]) == len(Bj["pores"]) and A["PT"] == Bj["PT"] and A["PC"] == Bj["PC"]
        md = 0.0; diff = {}
        for k in A:
            if A[k] == Bj[k]:
                continue
            try:
                d = float(np.max(np.abs(np.asarray(A[k], float) - np.asarray(Bj[k], float)))); diff[k] = d; md = max(md, d)
            except Exception:
                diff[k] = "non-numeric difference"; md = math.inf
        DETAILS["G1_diff_keys"] = diff
        row("G1", "geometry.json max |d| (nm) over differing numeric keys", f"{len(A['PT'])} prims, {len(A['pores'])} pores",
            f"{len(Bj['PT'])} prims, {len(Bj['pores'])} pores", md, f"<= {G1_TOL_NM} nm, same counts", "T" if same and md <= G1_TOL_NM else "fail", "must",
            f"differing keys: {sorted(diff)}")
    g1()

    @guard("G2", "geometry_fields.npz", "must")
    def g2():
        A = npz_dict(f"{R}/yc3d/geometry_fields.npz"); Bz = npz_dict(f"{P}/yc3d/geometry_fields.npz")
        if sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A):
            return row("G2", "real-shape SDF bricks (cell, vacuole, nucleus)", "identical", "identical", 0.0, "E", "E", "must")
        ok = sorted(A) == sorted(Bz); q = 0.0
        for k in A:
            if k.startswith("x") and ok and A[k].shape == Bz[k].shape:
                q = max(q, float(np.quantile(np.abs(A[k].astype(float) - Bz[k].astype(float)), 0.999)))
            elif not arrays_equal(A[k], Bz[k]):
                ok = False
        row("G2", "SDF bricks 99.9th pct |d| (nm)", "", "", q, f"<= {G2_TOL_NM} nm, same bricks", "T" if ok and q <= G2_TOL_NM else "fail", "must")
    g2()

    @guard("G4", "tile measures", "must")
    def g4():
        res = []
        for f in ("tile_measure", "tile_labvol"):
            A = pd.read_parquet(f"{R}/yc3d/{f}.parquet"); Bp = pd.read_parquet(f"{P}/yc3d/{f}.parquet")
            if A.equals(Bp):
                res.append((f, "E", 0.0, 1.0)); continue
            if list(A.index) != list(Bp.index) or list(A.columns) != list(Bp.columns):
                res.append((f, "fail", math.inf, 0.0)); continue
            tot = (A.sum() - Bp.sum()).abs() / A.sum().abs().clip(lower=1e-9)
            share = []
            for c in A.columns:
                thr = G4_TILE_FRAC_AREA if c.startswith("mem_") else G4_TILE_FRAC_VOL
                share.append(float((np.abs(A[c].values - Bp[c].values) <= thr).mean()))
            ok = float(tot.max()) <= G4_TOT_REL and min(share) >= G4_TILE_SHARE
            res.append((f, "T" if ok else "fail", float(tot.max()), min(share)))
        lev = "E" if all(r[1] == "E" for r in res) else ("T" if all(r[1] in "ET" for r in res) else "fail")
        n = len(pd.read_parquet(f"{P}/yc3d/tile_measure.parquet"))
        row("G4", "tile_measure + tile_labvol (max rel. total diff; min share of tiles within 0.5 %)", "5912 tiles", f"{n} tiles",
            "; ".join(f"{f}: {r:.2e} / {s:.4f}" for f, _, r, s in res), f"totals <= {G4_TOT_REL}, >= {G4_TILE_SHARE} tiles", lev, "must")
    g4()

    @guard("G5", "per-tile membrane meshes", "report")
    def g5():
        fa = sorted(os.path.basename(p) for p in glob.glob(f"{R}/yc3d/geom/*.npz"))
        fb = sorted(os.path.basename(p) for p in glob.glob(f"{P}/yc3d/geom/*.npz"))
        if fa != fb:
            return row("G5", "mesh tile set", len(fa), len(fb), None, "same set", "fail", "report")
        with Pool(NPROC) as pool:
            r = pool.starmap(_mesh_cmp, [(f"{R}/yc3d/geom/{f}", f"{P}/yc3d/geom/{f}") for f in fa], chunksize=32)
        same = sum(x[0] for x in r); aa = np.sum([x[1] for x in r], 0); ab = np.sum([x[2] for x in r], 0)
        rel = float(np.max(np.abs(aa - ab) / np.maximum(aa, 1e-9)))
        lev = "E" if same == len(fa) else ("T" if rel <= G5_AREA_REL else "fail")
        row("G5", "mesh tiles identical; max rel. membrane area diff", f"{len(fa)} tiles", f"{same} identical", rel, f"area <= {G5_AREA_REL}", lev, "report")
    g5()


def _mesh_cmp(a, b):
    A = npz_dict(a); Bz = npz_dict(b)
    return (sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A)), A["area"], Bz["area"]


def g3(bdir):
    @guard("G3", "cellgeom_v1 package validation", "must")
    def run():
        pj = json.load(open(os.path.join(bdir, "geometry", "package.json")))
        v = json.load(open(os.path.join(pj["package_dir"], "validation.json")))
        failed = [k for k, c in v["checks"].items() if not c["pass"]]
        row("G3", "cellgeom_v1 package: checksums, schema, native round trip, prims/fields/pores vs native", "-", pj["package_id"], None,
            "all checks pass", "E" if v["pass"] else "fail", "must", f"failed: {failed}" if failed else f"waived: {v.get('waived')}")
    run()


# ---------------------------------------------------------------- inventory, assembly
def inventory(R, P, bdir):
    @guard("I1", "ingredients.csv", "must")
    def i1():
        A = pd.read_csv(f"{R}/results/inventory/ingredients.csv"); Bc = pd.read_csv(f"{P}/results/inventory/ingredients.csv")
        if A.equals(Bc):
            return row("I1", "ingredients.csv", len(A), len(Bc), 0, "E", "E", "must")
        sa, sb = set(A.ing_id), set(Bc.ing_id)
        if sa != sb:
            return row("I1", "ingredient set", len(sa), len(sb), len(sa ^ sb), "same set", "fail", "must")
        A = A.set_index("ing_id"); Bc = Bc.set_index("ing_id").loc[A.index]
        l2 = np.abs(np.log2(np.maximum(Bc.copies_per_cell, 1e-12) / np.maximum(A.copies_per_cell, 1e-12)))
        l2 = float(l2[(A.copies_per_cell > 0) | (Bc.copies_per_cell > 0)].max()) if len(A) else 0.0
        wd = 0.0
        for x, y in zip(A.comp_weights.astype(str), Bc.comp_weights.astype(str)):
            if x != y:
                try:
                    dx, dy = json.loads(x.replace("'", '"')), json.loads(y.replace("'", '"'))
                    wd = max(wd, max(abs(dx.get(k, 0) - dy.get(k, 0)) for k in set(dx) | set(dy)))
                except Exception:
                    wd = math.inf
        row("I1", "ingredients: max |log2 copies ratio|; max comp-weight diff", len(A), len(Bc), f"{l2:.3g}; {wd:.3g}",
            f"<= {I1_LOG2}; <= {I1_W}", "T" if l2 <= I1_LOG2 and wd <= I1_W else "fail", "must")
    i1()

    @guard("I2", "protein_features.parquet", "must")
    def i2():
        A = pd.read_parquet(f"{R}/features/protein_features.parquet"); Bp = pd.read_parquet(f"{P}/features/protein_features.parquet")
        if A.equals(Bp):
            return row("I2", "protein_features.parquet", len(A), len(Bp), 0, "E", "E", "must")
        if list(A.columns) != list(Bp.columns) or len(A) != len(Bp) or set(A.orf) != set(Bp.orf):
            return row("I2", "protein_features layout", A.shape, Bp.shape, None, "same layout", "fail", "must")
        # rows keyed by ORF (row order is not content: ties of the copy-number sort may swap with ULP-level differences)
        A = A.set_index("orf"); Bp = Bp.set_index("orf").loc[A.index]
        worst = {}
        for c in A.columns:
            if A[c].equals(Bp[c]):
                continue
            if A[c].dtype.kind in "fi" and Bp[c].dtype.kind in "fi":
                d = float(np.nanmax(np.abs(A[c] - Bp[c]) / np.maximum(np.abs(A[c]), 1e-12)))
                worst[c] = d
            else:
                worst[c] = math.inf
        ok = all(v <= (I2_COPIES_REL if c == "copies_per_cell" else I2_OTHER_REL) for c, v in worst.items())
        DETAILS["I2_worst"] = worst
        row("I2", "protein_features: worst rel. diff per differing column", len(A), len(Bp), max(worst.values()) if worst else 0,
            f"copies <= {I2_COPIES_REL}, other <= {I2_OTHER_REL}", "T" if ok else "fail", "must", f"columns: {sorted(worst)[:8]}")
    i2()

    @guard("I3", "core copies step vs engine", "must")
    def i3():
        if bdir is None:
            return
        T = pd.read_csv(os.path.join(bdir, "inventory", "copies", "copies.tsv"), sep="\t")
        key = "orf" if "orf" in T.columns else T.columns[0]
        T = T.set_index(key)
        Bp = pd.read_parquet(f"{P}/features/protein_features.parquet").set_index("orf")
        A = pd.read_parquet(f"{R}/features/protein_features.parquet").set_index("orf")
        same = set(T.index) == set(Bp.index) and len(T) == len(Bp)
        c = T.copies.reindex(Bp.index)
        rel = float(np.nanmax(np.abs(c - Bp.copies_per_cell) / np.maximum(Bp.copies_per_cell, 1e-12)))
        rel5 = float(np.nanmax(np.abs(c - A.copies_per_cell.reindex(Bp.index)) / np.maximum(A.copies_per_cell.reindex(Bp.index), 1e-12)))
        row("I3", "core copies step vs port engine copies_per_cell (max rel.)", f"vs X5: {rel5:.2e}", f"{len(T)} proteins", rel,
            f"<= {I3_REL}, same 6,067 proteins", "E" if same and rel == 0 else ("T" if same and rel <= I3_REL else "fail"), "must")
    i3()

    @guard("I4", "s1 log numbers", "must")
    def i4():
        pats = [r"not in this cell .*", r"ribosomes from the cryo-ET density.*", r"curated complexes .*", r"pre-ribosomes: .*"]
        _log_lines("I4", f"{R}/s1_inventory.log", f"{P}/s1_inventory.log", pats, I4_REL, "must", "s1 log: NOT_IN_THIS_CELL, ribosomes, complexes, pre-ribosomes")
    i4()

    @guard("A1", "assembly tables", "must")
    def a1():
        res = []
        for f in ("oligo_map.json", "f_actin.json"):
            A = json.load(open(f"{R}/results/inventory/{f}")); Bj = json.load(open(f"{P}/results/inventory/{f}"))
            res.append((f, "E" if A == Bj else _json_tol(A, Bj, A1_REL)))
        for f in ("complexes_tried.csv", "complexes_added.csv", "heterodimers.csv", "oligomers.csv"):
            A = pd.read_csv(f"{R}/results/inventory/{f}"); Bc = pd.read_csv(f"{P}/results/inventory/{f}")
            res.append((f, "E" if A.equals(Bc) else _df_tol(A, Bc, A1_REL)))
        lev = "E" if all(r[1] == "E" for r in res) else ("T" if all(r[1] in "ET" for r in res) else "fail")
        row("A1", "assembly-state tables (oligo_map, f_actin, complexes tried/added, heterodimers, oligomers)", "6 files", "6 files",
            "; ".join(f"{f}:{l}" for f, l in res), f"same rows, rel <= {A1_REL}", lev, "must")
    a1()


def _json_tol(A, Bj, rel):
    la, lb = leaves(A), leaves(Bj)
    if sorted(la) != sorted(lb):
        return "fail"
    for k in la:
        x, y = la[k], lb[k]
        if x == y:
            continue
        if num(x) and num(y) and abs(x - y) <= rel * max(abs(x), 1e-12):
            continue
        if isinstance(x, list) and isinstance(y, list) and len(x) == len(y) and all(
                (a == b) or (num(a) and num(b) and abs(a - b) <= rel * max(abs(a), 1e-12)) for a, b in zip(x, y)):
            continue
        return "fail"
    return "T"


def _df_tol(A, Bc, rel):
    if list(A.columns) != list(Bc.columns) or len(A) != len(Bc):
        return "fail"
    k = A.columns[0]
    if set(A[k].astype(str)) != set(Bc[k].astype(str)):
        return "fail"
    # rows keyed by the first column (row order is not content); duplicated keys: sort by all columns as strings
    if A[k].is_unique:
        A = A.set_index(k); Bc = Bc.set_index(k).loc[A.index]
    else:
        A = A.sort_values(list(A.columns), key=lambda c: c.astype(str)).reset_index(drop=True)
        Bc = Bc.sort_values(list(Bc.columns), key=lambda c: c.astype(str)).reset_index(drop=True)
    for c in A.columns:
        if A[c].equals(Bc[c]):
            continue
        if A[c].dtype.kind in "fi" and Bc[c].dtype.kind in "fi":
            if float(np.nanmax(np.abs(A[c] - Bc[c]) / np.maximum(np.abs(A[c]), 1e-12))) > rel:
                return "fail"
        else:
            return "fail"
    return "T"


def _strip(line):
    return re.sub(r"^\d\d:\d\d:\d\d\s+", "", line.rstrip("\n"))


NUMRE = re.compile(r"-?\d+(?:,\d{3})*(?:\.\d+)?(?:e[-+]?\d+)?")


def _nums(s):
    return [float(x.replace(",", "")) for x in NUMRE.findall(s)]


def _log_lines(rid, fa, fb, pats, rel, cls, metric, absmin=0.0):
    la = [_strip(l) for l in open(fa)]; lb = [_strip(l) for l in open(fb)]
    worst = 0.0; lev = "E"; missing = []; per = {}
    for p in pats:
        ma = [l for l in la if re.match(p, l)]; mb = [l for l in lb if re.match(p, l)]
        if len(ma) != len(mb) or not ma:
            missing.append(p); lev = "fail"; continue
        for x, y in zip(ma, mb):
            if x == y:
                continue
            nx, ny = _nums(x), _nums(y)
            if len(nx) != len(ny):
                lev = "fail"; per[x[:60]] = "different number of values"; continue
            r = max((abs(a - b) / max(abs(a), 1e-12) for a, b in zip(nx, ny) if abs(a - b) > absmin), default=0.0)
            worst = max(worst, r); per[x[:60]] = r
            if lev == "E":
                lev = "T"
            if r > rel:
                lev = "fail"
    DETAILS[f"{rid}_lines"] = per
    row(rid, metric, f"{len(pats)} line patterns", "", worst, f"each number within {rel:g} rel.", lev, cls,
        f"missing/unmatched: {missing}" if missing else "")


# ---------------------------------------------------------------- library
def _lib_cmp(a, b, name):
    try:
        if name.endswith(".npz"):
            A = npz_dict(a, True); Bz = npz_dict(b, True)
            if sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A):
                return name, "E", ""
            if sorted(A) != sorted(Bz) or any(A[k].shape != Bz[k].shape for k in A):
                return name, "fail", "keys/shapes"
            for k in A:
                if arrays_equal(A[k], Bz[k]):
                    continue
                if A[k].dtype.kind not in "fc":
                    return name, "fail", f"{k} non-float differs"
                d = float(np.max(np.abs(A[k].astype(float) - Bz[k].astype(float))))
                tol = L2_BEADS_NM if k == "beads" else L2_TREE_NM
                if d > tol:
                    return name, "fail", f"{k} {d:.3g}"
            return name, "T", ""
        if name.endswith(".npy"):
            return name, ("E" if arrays_equal(np.load(a, allow_pickle=True), np.load(b, allow_pickle=True)) else "fail"), ""
        if name.endswith(".json"):
            A = json.load(open(a)); Bj = json.load(open(b))
            return name, ("E" if A == Bj else _json_tol(A, Bj, 1e-6)), ""
        return name, ("E" if open(a, "rb").read() == open(b, "rb").read() else "fail"), ""
    except Exception as e:  # noqa: BLE001
        return name, "fail", f"{type(e).__name__}: {e}"


def library(Rlib, Plib, start_manifest, R, P):
    start = {}
    for line in open(start_manifest):
        md5, size, mtime, rel = line.rstrip("\n").split("  ", 3)
        start[rel[2:] if rel.startswith("./") else rel] = md5

    @guard("L1", "library file set", "must")
    def l1():
        fa, fb = set(os.listdir(Rlib)), set(os.listdir(Plib))
        DETAILS["L1_only_x5"] = sorted(fa - fb); DETAILS["L1_only_port"] = sorted(fb - fa)
        row("L1", "library file set", len(fa), len(fb), len(fa ^ fb), "identical set", "E" if fa == fb else "fail", "must",
            f"only X5 {len(fa - fb)}, only port {len(fb - fa)}")
    l1()

    @guard("L2", "rebuilt library entries", "must")
    def l2():
        names = sorted((set(os.listdir(Rlib)) & set(os.listdir(Plib))) - set(start))
        with Pool(NPROC) as pool:
            r = pool.starmap(_lib_cmp, [(f"{Rlib}/{n}", f"{Plib}/{n}", n) for n in names], chunksize=64)
        nE = sum(x[1] == "E" for x in r); nT = sum(x[1] == "T" for x in r); nF = [x for x in r if x[1] == "fail"]
        DETAILS["L2_fail"] = nF[:200]; DETAILS["L2_T"] = [x[0] for x in r if x[1] == "T"][:500]
        shape_fail = any(x[2] == "keys/shapes" for x in nF)
        lev = "E" if nE == len(r) else ("T" if not shape_fail and (len(r) - len(nF)) / max(len(r), 1) >= L2_SHARE else "fail")
        row("L2", "rebuilt entries: identical / within tolerance / failing", len(r), f"E {nE}, T {nT}, fail {len(nF)}", len(nF),
            f"beads <= {L2_BEADS_NM} nm, trees <= {L2_TREE_NM} nm for >= {L2_SHARE:.1%}", lev, "must",
            f"first failures: {[x[0] + ' ' + x[2] for x in nF[:5]]}")
    l2()

    @guard("L3", "library_af.json", "must")
    def l3():
        A = json.load(open(f"{R}/results/library_af.json")); Bj = json.load(open(f"{P}/results/library_af.json"))
        la, lb = leaves(A), leaves(Bj)
        la = {k: v for k, v in la.items() if not VOLATILE.search(k)}; lb = {k: v for k, v in lb.items() if not VOLATILE.search(k)}
        row("L3", "library_af.json (AF models written / failed / missing)", _s({k: v for k, v in la.items() if num(v)}), _s({k: v for k, v in lb.items() if num(v)}),
            None, "identical", "E" if la == lb else "fail", "must")
    l3()

    @guard("L4", "start entries unchanged", "report")
    def l4():
        bad = []
        for n, m in start.items():
            p = f"{Plib}/{n}"
            if not os.path.exists(p) or hashlib.md5(open(p, "rb").read()).hexdigest() != m:
                bad.append(n)
        DETAILS["L4_changed"] = bad
        row("L4", "342 start entries unchanged at the end (md5)", len(start), len(start) - len(bad), len(bad), "identical", "E" if not bad else "fail", "report",
            f"changed: {bad[:8]}")
    l4()

    @guard("L5", "lib_ids, libpack_meta, pack trees", "must")
    def l5():
        ia = np.load(f"{R}/yc3d/lib_ids.npy", allow_pickle=True); ib = np.load(f"{P}/yc3d/lib_ids.npy", allow_pickle=True)
        ids_ok = arrays_equal(ia.astype(str), ib.astype(str))
        res = {}
        for f in ("libpack_meta.npz", "pack_trees_15.npz"):
            A = npz_dict(f"{R}/yc3d/{f}", True); Bz = npz_dict(f"{P}/yc3d/{f}", True)
            if sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A):
                res[f] = "E"; continue
            ok = sorted(A) == sorted(Bz) and all(A[k].shape == Bz[k].shape for k in A)
            md = 0.0
            for k in A:
                if not ok or arrays_equal(A[k], Bz[k]):
                    continue
                if A[k].dtype.kind in "fc":
                    md = max(md, float(np.max(np.abs(A[k] - Bz[k]))))
                else:
                    ok = False
            res[f] = "T" if ok and md <= L2_TREE_NM else f"fail ({md:.3g})"
        lev = "E" if ids_ok and all(v == "E" for v in res.values()) else ("T" if ids_ok and all(v in ("E", "T") for v in res.values()) else "fail")
        row("L5", "lib_ids.npy identical; libpack_meta / pack_trees_15", len(ia), len(ib), _s(res), f"ids identical; radii/centres <= {L2_TREE_NM} nm", lev, "must")
    l5()


# ---------------------------------------------------------------- pre-pass
PRE_GROUPS = {"npc": ["npc"], "chromatin_rabl": ["nuc", "seg", "chain"], "rdna_nucleolus": ["rpol", "r90", "rnuc", "rdna", "prerrna"],
              "mt_spb_kinetochores": ["mt", "spb", "kt", "kto"], "npc_fg_receptors": ["fg", "kap"],
              "pol2_origins_sites": ["pol2", "pol2dna", "pol2rna", "mcm", "orc", "tf", "lm"], "eisosomes_mcc": ["eis", "mcc"],
              "actin_endocytic": ["act", "tpm", "coat", "rvs", "pat"], "mito_atp_mtdna": ["atp", "mtdna"], "vacuole_cvt_nvj": ["cvt", "nvj"],
              "polysomes_mrna": ["poly", "mrna", "rnp"]}


def prepass(R, P):
    A = npz_dict(f"{R}/yc3d/prepass.npz"); Bz = npz_dict(f"{P}/yc3d/prepass.npz")
    row("S1", "prepass.npz array key set", len(A), len(Bz), len(set(A) ^ set(Bz)), "identical", "E" if sorted(A) == sorted(Bz) else "fail", "must")
    # S2 counts per array
    worst = 0.0; lev = "E"; bad = []
    for k in sorted(set(A) & set(Bz)):
        na, nb = (len(A[k]) if A[k].ndim else 1), (len(Bz[k]) if Bz[k].ndim else 1)
        if na != nb:
            lev = "T" if lev == "E" else lev
            d = abs(na - nb)
            worst = max(worst, d / max(na, 1))
            if d > max(S2_ABS, S2_REL * na):
                bad.append((k, na, nb)); lev = "fail"
    DETAILS["S2_count_diffs"] = {k: [len(A[k]) if A[k].ndim else 1, len(Bz[k]) if Bz[k].ndim else 1] for k in sorted(set(A) & set(Bz))
                                 if (len(A[k]) if A[k].ndim else 1) != (len(Bz[k]) if Bz[k].ndim else 1)}
    row("S2", "instance/bead counts per pre-pass array (max rel. diff)", "", "", worst, f"each within max({S2_ABS}, {S2_REL:.0%})", lev, "must",
        f"outside: {bad[:6]}")
    # S4 positions per structure group
    same = {}; worst_c = 0.0; worst_r = 0.0; lev4 = "E"
    for g, keys in PRE_GROUPS.items():
        ks = [k for k in A if any(k == p or k.startswith(p + "_") for p in keys)]
        eq = all(k in Bz and arrays_equal(A[k], Bz[k]) for k in ks)
        same[g] = eq
        if eq:
            continue
        lev4 = "T" if lev4 == "E" else lev4
        for k in ks:
            if (k.endswith("_pos") or k.endswith("_pts")) and k in Bz and len(A[k]) and len(Bz[k]):
                ca, cb = A[k].astype(float).mean(0), Bz[k].astype(float).mean(0)
                ra, rb = np.sqrt((A[k].astype(float) ** 2).sum(1).mean()), np.sqrt((Bz[k].astype(float) ** 2).sum(1).mean())
                dc = float(np.linalg.norm(ca - cb)); dr = float(abs(ra - rb) / ra)
                worst_c = max(worst_c, dc); worst_r = max(worst_r, dr)
                if dc > S4_CENTROID_NM or dr > S4_RMS_REL:
                    lev4 = "fail"
    DETAILS["S4_identical_groups"] = same
    row("S4", "pre-pass positions per structure group (identical groups; worst centroid nm / RMS rel.)", f"{len(same)} groups",
        f"{sum(same.values())} identical", f"{worst_c:.3g} / {worst_r:.3g}", f"centroid <= {S4_CENTROID_NM} nm, RMS <= {S4_RMS_REL:.0%}", lev4, "report",
        f"differing: {[g for g, v in same.items() if not v]}")


def s3_lines(R, P):
    pats = [r"microtubules: .*", r"kinetochores: .*", r"outer kinetochores: .*", r"Rabl layout: .*", r"nucleosome map: .*", r"chromatin: \d.*",
            r"pore receptors: .*", r"FG mesh: .*", r"rDNA: .*", r"rDNA chromatin: .*", r"engaged Pol II: .*", r"origins: .*",
            r"chromatin sites: .*", r"tRNA genes, telomeres, silent loci: .*", r"eisosomes: .*", r"MCC proteins .*", r"actin cables: .*",
            r"tropomyosin: .*", r"endocytic coat: .*", r"amphiphysin: .*", r"actin patches: \d.*", r"ATP synthase dimers: .*",
            r"mtDNA: .*", r"Cvt complex: .*", r"nucleus-vacuole junction: .*", r"polysomes: .*", r"RAC on polysome ribosomes: .*",
            r"mRNA: .*", r"mRNP ends: .*"]
    _log_lines("S3", f"{R}/s3_pack.log", f"{P}/s3_pack.log", pats, S3_REL, "must", "s3 log structure numbers (29 builder lines)")


# ---------------------------------------------------------------- packing
def _tile_cmp(a, b):
    A = npz_dict(a); Bz = npz_dict(b)
    same = sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A)
    def jobs(Z):
        return pd.DataFrame(dict(lib=Z["job_lib"].astype(np.int64), region=Z["job_region"].astype(np.int64), placed=Z["job_placed"].astype(np.int64)))
    return same, jobs(A).groupby(["lib", "region"]).placed.sum(), jobs(Bz).groupby(["lib", "region"]).placed.sum(), len(A["types"]), len(Bz["types"])


def packing(R, P):
    @guard("P1", "expected_totals", "must")
    def p1():
        A = pd.read_parquet(f"{R}/yc3d/expected_totals.parquet"); Bp = pd.read_parquet(f"{P}/yc3d/expected_totals.parquet")
        if A.equals(Bp):
            return row("P1", "expected_totals.parquet", len(A), len(Bp), 0, "E", "E", "must")
        k = ["lib_id", "ing_id", "region"]
        M = A.merge(Bp, on=k, how="outer", suffixes=("_x", "_p"))
        ok = not M.isna().any().any()
        rel = float(np.nanmax(np.abs(M.expected_x - M.expected_p) / np.maximum(M.expected_x.abs(), 1e-12))) if ok else math.inf
        row("P1", "expected copies per (type, region): max rel. diff", len(A), len(Bp), rel, f"<= {P1_REL}", "T" if ok and rel <= P1_REL else "fail", "must")
    p1()

    @guard("P2", "organisation layer", "must")
    def p2():
        res = {}
        for f in ("network.npz", "partners.npz"):
            A = npz_dict(f"{R}/organisation/{f}"); Bz = npz_dict(f"{P}/organisation/{f}")
            if sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A):
                res[f] = "E"; continue
            ok = sorted(A) == sorted(Bz); md = 0.0
            for k in A:
                if not ok or arrays_equal(A[k], Bz[k]):
                    continue
                if f == "partners.npz" and k in ("off", "idx", "kind", "ids"):
                    ok = False; continue
                if A[k].shape != Bz[k].shape or A[k].dtype.kind not in "fc":
                    ok = False; continue
                md = max(md, float(np.max(np.abs(A[k] - Bz[k]) / np.maximum(np.abs(A[k]), 1e-6))))
            res[f] = "T" if ok and md <= P2_REL else f"fail ({md:.3g})"
        lev = "E" if all(v == "E" for v in res.values()) else ("T" if all(v in ("E", "T") for v in res.values()) else "fail")
        row("P2", "organisation network.npz / partners.npz", "", "", _s(res), f"partner CSR identical; floats rel <= {P2_REL}", lev, "must")
    p2()

    @guard("P3", "tile_counts + thinning", "must")
    def p3():
        A = pd.read_parquet(f"{R}/yc3d/tile_counts.parquet"); Bp = pd.read_parquet(f"{P}/yc3d/tile_counts.parquet")
        ta = json.load(open(f"{R}/yc3d/thinning.json")); tb = json.load(open(f"{P}/yc3d/thinning.json"))
        lt = leaves(ta); lb = leaves(tb)
        thin_d = max((abs(lt[k] - lb.get(k, math.inf)) for k in lt if num(lt[k])), default=0.0)
        if A.equals(Bp) and ta == tb:
            return row("P3", "tile_counts.parquet + thinning.json", len(A), len(Bp), 0, "E", "E", "must")
        ga = A.groupby(["lib", "region"]).n.sum(); gb = Bp.groupby(["lib", "region"]).n.sum()
        J = pd.concat([ga.rename("a"), gb.rename("b")], axis=1).fillna(0)
        tol = np.maximum(3 * np.sqrt(J.a), P3_REL * J.a)
        nbad = int((np.abs(J.a - J.b) > tol).sum())
        DETAILS["P3_outside"] = J[np.abs(J.a - J.b) > tol].head(50).reset_index().values.tolist()
        row("P3", "requested per (type, region) outside tolerance; thinning max diff", f"{len(ga)} pairs, {int(J.a.sum())} copies",
            f"{len(gb)} pairs, {int(J.b.sum())} copies", f"{nbad} pairs; thin {thin_d:.3g}", f"max(3 sqrt n, {P3_REL:.1%} n); thinning +-{P3_THIN}",
            "T" if nbad == 0 and thin_d <= P3_THIN else "fail", "must")
    p3()

    @guard("P5", "placed per (type, region)", "must")
    def p5():
        fa = sorted(os.path.basename(p) for p in glob.glob(f"{R}/yc3d/tiles_raw/*.npz"))
        fb = set(os.path.basename(p) for p in glob.glob(f"{P}/yc3d/tiles_raw/*.npz"))
        common = [f for f in fa if f in fb]
        with Pool(NPROC) as pool:
            r = pool.starmap(_tile_cmp, [(f"{R}/yc3d/tiles_raw/{f}", f"{P}/yc3d/tiles_raw/{f}") for f in common], chunksize=16)
        same = sum(x[0] for x in r)
        na = np.array([x[3] for x in r]); nb = np.array([x[4] for x in r])
        ga = pd.concat([x[1] for x in r]).groupby(level=[0, 1]).sum(); gb = pd.concat([x[2] for x in r]).groupby(level=[0, 1]).sum()
        J = pd.concat([ga.rename("a"), gb.rename("b")], axis=1).fillna(0)
        big = J[J.a >= 10]; small = J[J.a < 10]
        bad_big = int((np.abs(big.a - big.b) > np.maximum(3 * np.sqrt(big.a), P5_REL * big.a)).sum())
        share_small = float((np.abs(small.a - small.b) <= P5_SMALL_ABS).mean()) if len(small) else 1.0
        DETAILS["P5_outside_big"] = big[np.abs(big.a - big.b) > np.maximum(3 * np.sqrt(big.a), P5_REL * big.a)].head(100).reset_index().values.tolist()
        lev7 = "E" if same == len(fa) and len(fb) == len(fa) else "fail"
        row("P7", "packed tiles: arrays identical (types, pos, quat, job table)", f"{len(fa)} tiles", f"{same} identical", len(fa) - same, "identical", lev7, "report")
        n_tiles_eq = int((na == nb).sum())
        row("P4", "instances per tile identical", f"{len(fa)} tiles", f"{n_tiles_eq} equal", int(len(fa) - n_tiles_eq), "identical",
            "E" if n_tiles_eq == len(fa) and len(fb) == len(fa) else "fail", "report", f"sum |d| {int(np.abs(na - nb).sum())}")
        lev = "E" if ga.equals(gb) else ("T" if bad_big == 0 and share_small >= P5_SMALL_SHARE else "fail")
        row("P5", "placed per (type, region): pairs n>=10 outside; share of n<10 pairs within +-3", f"{len(big)} / {len(small)} pairs; {int(J.a.sum())} placed",
            f"{int(J.b.sum())} placed", f"{bad_big}; {share_small:.4f}", f"max(3 sqrt n, {P5_REL:.0%} n); >= {P5_SMALL_SHARE:.0%}", lev, "must")
    p5()

    @guard("P6", "placed totals", "must")
    def p6():
        def tot(f):
            for l in open(f):
                m = re.search(r"placed ([\d.]+)M / ([\d.]+)M = ([\d.]+)", l)
                if m and "total" in l:
                    return float(m.group(3)), l.strip()
            return None, None
        fa, la = tot(f"{R}/s3_pack.log"); fb, lb = tot(f"{P}/s3_pack.log")
        ia = json.load(open(f"{R}/results/wholecell3d_validation.json"))["instances_total"]
        ib = json.load(open(f"{P}/results/wholecell3d_validation.json"))["instances_total"]
        lev = "E" if (la and lb and _strip(la).split("total ")[1].split(";")[1] == _strip(lb).split("total ")[1].split(";")[1] and ia == ib) else (
            "T" if abs(fa - fb) <= P6_FRAC and abs(ia - ib) <= P6_INST_REL * ia else "fail")
        row("P6", "placed fraction (s3 log); instances_total", f"{fa}; {ia}", f"{fb}; {ib}", f"{fb - fa:+.4f}; {(ib - ia) / ia:+.4%}",
            f"+-{P6_FRAC}; +-{P6_INST_REL:.1%}", lev, "must")
    p6()


# ---------------------------------------------------------------- validation
def _v_tol_v1(path, x, y):
    if isinstance(x, int) and isinstance(y, int):
        return abs(x - y) <= max(3, 0.01 * abs(x))
    p = path.lower()
    if "clash" in p and "rate" in p or "outside" in p:
        return abs(x - y) <= max(V1_RATE_REL * abs(x), V1_RATE_ABS)
    if "median_angle" in p:
        return abs(x - y) <= V1_ANGLE
    if "density" in p:
        return abs(x - y) <= V1_DENS
    if 0 <= x <= 1 and 0 <= y <= 1:
        return abs(x - y) <= V1_FRAC
    return abs(x - y) <= max(V2_REL * abs(x), V2_ABS)


def _v_tol_v2(path, x, y):
    return abs(x - y) <= max(V2_REL * abs(x), V2_ABS)


def json_cmp(a, b, tolfn):
    """-> (level, n_scalars, bad scalars, list summary)"""
    A = json.load(open(a)); Bj = json.load(open(b))
    la = {k: v for k, v in leaves(A).items() if not VOLATILE.search(k.split(".")[-1])}
    lb = {k: v for k, v in leaves(Bj).items() if not VOLATILE.search(k.split(".")[-1])}
    if la == lb:
        return "E", len(la), [], {}
    bad = []; lists = {}
    for k in sorted(set(la) | set(lb)):
        x, y = la.get(k), lb.get(k)
        if x == y:
            continue
        if num(x) and num(y):
            if not (math.isfinite(x) and math.isfinite(y)) or not tolfn(k, x, y):
                bad.append((k, x, y))
        elif isinstance(x, list) and isinstance(y, list):
            if len(x) != len(y):
                lists[k] = f"length {len(x)} vs {len(y)}"
            else:
                out = sum(1 for p, q in zip(x, y) if p != q and not (num(p) and num(q) and _v_tol_v2(k, p, q)))
                lists[k] = f"{out}/{len(x)} elements outside"
        elif isinstance(x, str) and isinstance(y, str):
            continue                                  # labels / free text: not scalars
        else:
            bad.append((k, x, y))
    return ("T" if not bad else "fail"), len(la), bad, lists


GPU_FILES = {"diffusion_benchmark.json"}


def validation(R, P):
    @guard("V1", "wholecell3d_validation.json", "must")
    def v1():
        lev, n, bad, lists = json_cmp(f"{R}/results/wholecell3d_validation.json", f"{P}/results/wholecell3d_validation.json", _v_tol_v1)
        DETAILS["V1_bad"] = bad; DETAILS["V1_lists"] = lists
        row("V1", "stage-4 validation: numeric scalars outside tolerance", n, n - len(bad), len(bad), "per leaf class (header)", lev, "must",
            f"{bad[:6]}")
    v1()
    files = sorted(glob.glob(f"{R}/validation/*.json")) + sorted(glob.glob(f"{R}/validation/latest/*.json"))
    res = {}
    for f in files:
        rel = os.path.relpath(f, R)
        if os.path.basename(f) in ("MANIFEST.json",) or os.path.basename(f) in GPU_FILES:
            continue
        g = os.path.join(P, rel)
        if not os.path.exists(g):
            res[rel] = ("fail", 0, [("missing in port", None, None)], {}); continue
        try:
            res[rel] = json_cmp(f, g, _v_tol_v2)
        except Exception as e:  # noqa: BLE001
            res[rel] = ("fail", 0, [(f"{type(e).__name__}: {e}", None, None)], {})
    DETAILS["V2"] = {k: dict(level=v[0], n=v[1], bad=v[2][:30], lists=v[3]) for k, v in res.items()}
    for k, v in res.items():
        print(f"   V2 {v[0]:4s} {k}  bad={v[2][:3]}", flush=True)
    lev = "E" if all(v[0] == "E" for v in res.values()) else ("T" if all(v[0] in "ET" for v in res.values()) else "fail")
    row("V2", "validate_all CPU comparisons (files E / T / fail)", f"{len(res)} files",
        f"E {sum(v[0] == 'E' for v in res.values())}, T {sum(v[0] == 'T' for v in res.values())}, fail {sum(v[0] == 'fail' for v in res.values())}",
        None, f"scalars within max({V2_REL:.0%}, {V2_ABS})", lev, "must", f"failing: {[k for k, v in res.items() if v[0] == 'fail']}")
    # anchors as its own visible sub-row (named in the criteria)
    if "validation/anchors.json" in res:
        v = res["validation/anchors.json"]
        row("V2a", "anchors.json (sub-row of V2)", "", "", len(v[2]), "as V2", v[0], "info", f"{v[2][:4]}")

    @guard("V3", "GPU steps", "must")
    def v3():
        parts = {}
        A = npz_dict(f"{R}/contacts/contacts.npz"); Bz = npz_dict(f"{P}/contacts/contacts.npz")
        if sorted(A) == sorted(Bz) and all(arrays_equal(A[k], Bz[k]) for k in A):
            parts["contacts"] = ("E", 0.0)
        else:
            ca = sum(int(A[k].sum()) for k in A if re.fullmatch(r"C\d+", k)); cb = sum(int(Bz[k].sum()) for k in Bz if re.fullmatch(r"C\d+", k))
            d = abs(ca - cb) / max(ca, 1); parts["contacts"] = ("T" if d <= V3_CONTACT_REL else "fail", d)
        ea = pd.read_parquet(f"{R}/contacts/encounters_kinases.parquet"); eb = pd.read_parquet(f"{P}/contacts/encounters_kinases.parquet")
        if ea.equals(eb):
            parts["encounters"] = ("E", 0.0)
        else:
            cols = [c for c in ("p_near_2", "p_near_10", "c_local_uM") if c in ea.columns]
            d = max(abs(ea[c].mean() - eb[c].mean()) / max(abs(ea[c].mean()), 1e-12) for c in cols) if len(ea) and len(eb) else math.inf
            d = max(d, abs(len(ea) - len(eb)) / max(len(ea), 1))
            parts["encounters"] = ("T" if d <= V3_ENC_REL else "fail", d)
        for j, sd in V3_PHI_SD.items():
            a = json.load(open(f"{R}/sweeps/{j}.json")); b = json.load(open(f"{P}/sweeps/{j}.json"))
            d = abs(a["phi"] - b["phi"]); parts[j] = ("E" if a == b else ("T" if d <= 2 * sd else "fail"), d)
        a = json.load(open(f"{R}/sweeps/nuc100_mf.json")); b = json.load(open(f"{P}/sweeps/nuc100_mf.json"))
        d = abs(a.get("phi", 0) - b.get("phi", 0)); parts["nuc100_mf"] = ("E" if a == b else ("T" if d <= V3_NUC_PHI else "fail"), d)
        a = json.load(open(f"{R}/validation/diffusion_benchmark.json")); b = json.load(open(f"{P}/validation/diffusion_benchmark.json"))
        if a == b:
            parts["diffusion_benchmark"] = ("E", 0.0)
        else:
            la, lb = leaves(a), leaves(b)
            ks = [k for k in la if num(la[k]) and num(lb.get(k)) and re.search(r"(D_over_D0|ratio|D/D0|pred)", k)]
            d = max((abs(la[k] - lb[k]) / max(abs(la[k]), 1e-12) for k in ks), default=None)
            if d is None:
                lev2, _, bad, _ = json_cmp(f"{R}/validation/diffusion_benchmark.json", f"{P}/validation/diffusion_benchmark.json",
                                           lambda p, x, y: abs(x - y) <= max(V3_DIFF_REL * abs(x), V2_ABS))
                parts["diffusion_benchmark"] = (lev2, len(bad))
            else:
                parts["diffusion_benchmark"] = ("T" if d <= V3_DIFF_REL else "fail", d)
        DETAILS["V3"] = parts
        lev = "E" if all(v[0] == "E" for v in parts.values()) else ("T" if all(v[0] in "ET" for v in parts.values()) else "fail")
        row("V3", "GPU steps: contacts, kinase encounters, sweeps (phi), diffusion benchmark", "", "",
            "; ".join(f"{k}:{v[0]}({_s(v[1])})" for k, v in parts.items()), "contacts 1 %, encounters 5 %, phi 2 SD, D/D0 10 %", lev, "must")
    v3()

    @guard("V4", "condensates, microscopy", "report")
    def v4():
        pats = [r"log phase .*", r"stress .*"]
        la = [l for l in open(f"{R}/condensates/cond_test.log") if re.match("|".join(pats), l)]
        lb = [l for l in open(f"{P}/condensates/cond_test.log") if re.match("|".join(pats), l)]
        a = json.load(open(f"{R}/microscopy/spot_counts.json")); b = json.load(open(f"{P}/microscopy/spot_counts.json"))
        lev = "E" if la == lb and a == b else "T"
        for x, y in zip(la, lb):
            nx, ny = _nums(x), _nums(y)
            if len(nx) != len(ny) or any(abs(p - q) > max(V4_REL * abs(p), V4_ONSET) for p, q in zip(nx, ny)):
                lev = "fail"
        if len(la) != len(lb) or not la:
            lev = "fail"
        if a != b and any(abs(a[k] - b.get(k, -99)) > max(1, V4_REL * a[k]) for k in a):
            lev = "fail" if lev != "fail" else lev
        row("V4", "condensate negative control lines; microscopy spot counts", _s(a), _s(b), None, f"counts {V4_REL:.0%}, onset +-{V4_ONSET}", lev, "report")
    v4()

    @guard("V5", "build stamp check", "must")
    def v5():
        M = json.load(open(f"{P}/validation/MANIFEST.json"))
        stale = [k for k, v in M["files"].items() if not v["of_this_build"]]
        row("V5", "build_stamp check: result files of this build / older", len(json.load(open(f"{R}/validation/MANIFEST.json"))["files"]),
            len(M["files"]), len(stale), "0 stale", "E" if not stale else "fail", "must", f"stale: {stale[:5]}")
    v5()


# ---------------------------------------------------------------- shared-pipeline outputs
def export_checks(bdir, R, P):
    @guard("X1", "YCT2 export", "must")
    def x1():
        web = os.path.join(bdir, "export", "web")
        ec = json.load(open(os.path.join(web, "encode_check.json")))
        idx = json.load(open(os.path.join(web, "tiles", "index.json")))
        regs = json.load(open(os.path.join(web, "regions.json")))["regions"]
        S = json.load(open(os.path.join(bdir, "pack", "summary.json")))
        sys.path.insert(0, REPO)
        from wholecell.export import yct2
        # quaternion sign + region ids on a sample of tiles
        rng = np.random.default_rng(0); rows_ = idx["tiles"]; sample = [rows_[i] for i in rng.choice(len(rows_), min(200, len(rows_)), replace=False)]
        wneg = 0; maxreg = 0; nreg = len(regs)
        for t in sample:
            with open(os.path.join(web, "tiles", t[3]), "rb") as fh:
                fh.seek(t[4]); raw = gzip.decompress(fh.read(t[5]))
            D = yct2.decode_tile(raw)
            wneg += int((D["quat"][:, 0] < -1e-6).sum()); maxreg = max(maxreg, int(D["region"].max()) if len(D["region"]) else 0)
        ok = (ec["count_equal"] and ec["pos_ok"] and ec["rot_ok"] and ec["instances_source"] == S["instances"] and wneg == 0
              and maxreg < nreg and idx["grid"]["origin_nm"] == [0.0, 0.0, 0.0] and idx["tile_nm"] == 200.0)
        DETAILS["X1"] = dict(encode_check=ec, n_regions=nreg, max_region_id_seen=maxreg, w_negative=wneg, grid=idx["grid"])
        row("X1", "YCT2: counts, pos/rot error, w >= 0, regions (u16 palette), grid from config/geometry", f"source {S['instances']}",
            f"decoded {ec['instances_decoded']}; regions {nreg}; max id {maxreg}", f"pos {ec['max_pos_err_nm']:.4g} nm, rot {ec['max_rot_err_rad']:.3g}",
            "all checks", "E" if ok else "fail", "must", f"grid {idx['grid']}")
    x1()

    @guard("X2", "yeast stage-5 viewer encoding", "report")
    def x2():
        wa, wb = f"{R}/yc3d/web/data", f"{P}/yc3d/web/data"
        ta = sorted(os.path.basename(p) for p in glob.glob(f"{wa}/tiles/*.bin.gz"))
        tb = set(os.path.basename(p) for p in glob.glob(f"{wb}/tiles/*.bin.gz"))
        with Pool(NPROC) as pool:
            r = pool.starmap(_gz_eq, [(f"{wa}/tiles/{f}", f"{wb}/tiles/{f}") for f in ta if f in tb], chunksize=32)
        lib = {}
        for f in sorted(glob.glob(f"{wa}/library/*")):
            g = os.path.join(wb, "library", os.path.basename(f))
            lib[os.path.basename(f)] = os.path.exists(g) and _gz_eq(f, g)
        n_eq = sum(r)
        DETAILS["X2_library"] = lib
        row("X2", "viewer tiles (decompressed) identical; library files identical", f"{len(ta)} tiles, {len(lib)} lib files",
            f"{n_eq} tiles, {sum(lib.values())} lib files", None, "identical", "E" if n_eq == len(ta) and len(tb) == len(ta) and all(lib.values()) else "fail", "report")
    x2()


def _gz_eq(a, b):
    def rd(p):
        d = open(p, "rb").read()
        return gzip.decompress(d) if d[:2] == b"\x1f\x8b" else d
    try:
        return rd(a) == rd(b)
    except Exception:
        return False


# ---------------------------------------------------------------- main
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
NPROC = 8


def verdict(out):
    must = [r for r in ROWS if r["cls"] == "must"]; rep = [r for r in ROWS if r["cls"] == "report"]
    v = "NOT REPRODUCED" if any(r["level"] == "fail" for r in must) else (
        "REPRODUCED WITH DEVIATIONS" if any(r["level"] == "fail" for r in rep) else "REPRODUCED")
    lev = "E" if all(r["level"] == "E" for r in must) else "T"
    res = dict(verdict=v, must_level=lev if v != "NOT REPRODUCED" else "fail", criteria="docs/reference/history/yeast_x5_reproduction_criteria.md",
               counts={c: {l: sum(1 for r in ROWS if r["cls"] == c and r["level"] == l) for l in ("E", "T", "fail")} for c in ("must", "report")},
               failed_must=[r["id"] for r in must if r["level"] == "fail"], failed_report=[r["id"] for r in rep if r["level"] == "fail"],
               first_not_exact=[r["id"] for r in ROWS if r["level"] != "E" and r["cls"] in ("must", "report")][:1], rows=ROWS)
    os.makedirs(out, exist_ok=True)
    pd.DataFrame(ROWS, columns=["id", "metric", "x5", "port", "delta", "tolerance", "level", "cls", "note"]).to_csv(
        os.path.join(out, "criteria_results.tsv"), sep="\t", index=False)
    json.dump(res, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    json.dump(DETAILS, open(os.path.join(out, "details.json"), "w"), indent=1, default=str)
    print("VERDICT:", v, json.dumps(res["counts"]))
    return res


def main():
    global NPROC
    ap = argparse.ArgumentParser()
    ap.add_argument("--bdir", default=None)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--self-test", action="store_true", help="compare X5 with itself (engine-only rows)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--only", default=None, help="comma list of sections (geometry,inventory,library,prepass,packing,validation,export)")
    a = ap.parse_args()
    NPROC = a.nproc
    R = os.path.join(a.ref, "work")
    if a.self_test:
        P, Plib, bdir = R, os.path.join(a.ref, "lib"), None
    else:
        bdir = a.bdir; P = os.path.join(bdir, "engine", "work"); Plib = os.path.join(bdir, "engine", "lib")
    out = a.out or os.path.join(bdir, "validate", "compare_x5")
    only = set(a.only.split(",")) if a.only else None
    sec = lambda s: only is None or s in only
    frozen = None
    try:
        cfg = json.load(open(os.path.join(bdir, "config.resolved.json"))) if bdir else None
        frozen = cfg["pipeline_engine"]["frozen_dir"] if cfg else None
    except Exception:
        pass
    try:                     # core yeast builds (A9: no pipeline_engine block): the folder of the start-entry manifest
        frozen = frozen or os.path.dirname(cfg["frozen_inputs"]["yeast_x5_lib_start"]["manifest"])
    except Exception:
        pass
    frozen = frozen or os.path.join(REPO, "..", "data", "derived", "wholecell_frozen", "yeast_x5")
    if sec("geometry"):
        geometry(R, P)
        if bdir:
            g3(bdir)
    if sec("inventory"):
        inventory(R, P, bdir)
    if sec("library"):
        library(os.path.join(a.ref, "lib"), Plib, os.path.join(frozen, "manifest_lib_start.tsv"), R, P)
    if sec("prepass"):
        guard("S", "prepass", "must")(prepass)(R, P)
        guard("S3", "s3 log", "must")(s3_lines)(R, P)
    if sec("packing"):
        packing(R, P)
    if sec("validation"):
        validation(R, P)
    if sec("export") and bdir:
        export_checks(bdir, R, P)
    verdict(out)


if __name__ == "__main__":
    main()
