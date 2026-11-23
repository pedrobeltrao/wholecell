"""A2 check: the core tiling + request steps on the yeast config vs the yeast port build (criteria G4, G5, P3 of
docs/reference/history/yeast_x5_reproduction_criteria.md; level E = identical arrays).

    python -m wholecell pack workers=3 --cell <dev yaml> --build B --step tiles     # core tiling (pack.tiling)
    python -m wholecell pack workers=3 --cell <dev yaml> --build B --step counts    # pack.requests legacy_x5
    python tests/check_a2_yeast_tiles.py --core B --ref <port build> [--out result.json] [--workers 3]

The dev yaml extends config/cells/yeast.yaml with the request inputs pointed at the port build (declared inputs:
expected_totals = criterion P1, lib_ids, pack_trees_15, thinning.json, organisation loadings).
  G4  tiles.parquet (named by tile, columns = pack.regions) vs yc3d/tile_measure.parquet; tile_labvol.parquet vs
      yc3d/tile_labvol.parquet (label columns matched by engine code)
  G5  tile_geom/<tile>.npz vs yc3d/geom/<tile>.npz: same tile set, same keys, identical arrays (dtype included)
  P3  counts.npz vs yc3d/tile_counts.parquet: identical rows in order; else identical as a multiset; per (type, region)
      totals; thinning factors used vs yc3d/thinning.json
"""
import argparse
import glob
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


def _name(ijk):
    return f"t_{int(ijk[0]):+03d}_{int(ijk[1]):+03d}_{int(ijk[2]):+03d}"


def g4(core, ref):
    meta = json.load(open(f"{core}/pack/tiles_meta.json"))
    M = pd.read_parquet(f"{core}/pack/tiles.parquet")
    off = np.round(np.asarray(meta["lo0"]) / meta["TILE"]).astype(np.int64)
    names = [_name(t) for t in M[["ti", "tj", "tk"]].values + off]
    A = M[meta["columns"]].set_axis(names, axis=0)
    R = pd.read_parquet(f"{ref}/engine/work/yc3d/tile_measure.parquet")
    out = dict(n_core=len(A), n_ref=len(R), same_tiles=list(A.index) == list(R.index),
               same_columns=list(A.columns) == list(R.columns))
    out["measure_identical"] = bool(out["same_tiles"] and out["same_columns"] and A.equals(R))
    if out["same_tiles"] and out["same_columns"] and not out["measure_identical"]:
        d = np.abs(A.values - R.values)
        out["measure_max_abs_diff"] = {c: float(x) for c, x in zip(A.columns, d.max(0)) if x > 0}
        out["measure_values_differing"] = int((d > 0).sum())
    LV = pd.read_parquet(f"{core}/pack/tile_labvol.parquet")
    RL = pd.read_parquet(f"{ref}/engine/work/yc3d/tile_labvol.parquet")
    out["labvol_columns"] = dict(zip(RL.columns, LV.columns))          # yeast cell3d name -> package name (same code)
    out["labvol_identical"] = bool(LV.shape == RL.shape and np.array_equal(LV.values, RL.values)
                                   and LV.values.dtype == RL.values.dtype)
    if LV.shape == RL.shape and not out["labvol_identical"]:
        out["labvol_values_differing"] = int((LV.values != RL.values).sum())
        out["labvol_max_abs_diff"] = float(np.abs(LV.values - RL.values).max())
    out["pass"] = out["measure_identical"] and out["labvol_identical"]
    return out


def _cmp(a, b):
    A, B = np.load(a), np.load(b)
    if sorted(A.files) != sorted(B.files):
        return False, None, None, 0.0
    same = all(A[k].dtype == B[k].dtype and np.array_equal(A[k], B[k]) for k in A.files)
    dv = 0.0
    if not same and A["verts"].shape == B["verts"].shape:
        dv = float(np.abs(A["verts"].astype(np.float64) - B["verts"].astype(np.float64)).max()) if len(A["verts"]) else 0.0
    return same, A["area"], B["area"], dv


def g5(core, ref, workers):
    fa = sorted(os.path.basename(p) for p in glob.glob(f"{core}/pack/tile_geom/*.npz"))
    fb = sorted(os.path.basename(p) for p in glob.glob(f"{ref}/engine/work/yc3d/geom/*.npz"))
    out = dict(n_core=len(fa), n_ref=len(fb), same_set=fa == fb)
    if not out["same_set"]:
        out["pass"] = False
        return out
    with Pool(workers) as pool:
        r = pool.starmap(_cmp, [(f"{core}/pack/tile_geom/{f}", f"{ref}/engine/work/yc3d/geom/{f}") for f in fa], chunksize=32)
    same = [x[0] for x in r]
    aa = np.sum([x[1] for x in r if x[1] is not None], 0); ab = np.sum([x[2] for x in r if x[2] is not None], 0)
    out.update(identical_tiles=int(sum(same)), differing=[f for f, s in zip(fa, same) if not s][:20],
               max_abs_vert_diff_nm=float(max(x[3] for x in r)),
               area_rel_diff=[float(x) for x in np.abs(aa - ab) / np.maximum(ab, 1e-9)])
    out["pass"] = out["identical_tiles"] == len(fa)
    return out


def p3(core, ref):
    z = np.load(f"{core}/pack/counts.npz")
    C = pd.DataFrame({k: z[k] for k in ("tile", "lib", "region", "n")})
    R = pd.read_parquet(f"{ref}/engine/work/yc3d/tile_counts.parquet")
    out = dict(rows_core=len(C), rows_ref=len(R), copies_core=int(C.n.sum()), copies_ref=int(R.n.sum()))
    out["identical_in_order"] = bool(C.equals(R.reset_index(drop=True)))
    key = ["tile", "lib", "region"]
    Cs = C.sort_values(key, kind="stable").reset_index(drop=True); Rs = R.sort_values(key, kind="stable").reset_index(drop=True)
    out["identical_as_multiset"] = bool(len(Cs) == len(Rs) and all(np.array_equal(Cs[k].values, Rs[k].values) for k in C.columns))
    ga = R.groupby(["lib", "region"]).n.sum(); gb = C.groupby(["lib", "region"]).n.sum()
    J = pd.concat([ga.rename("ref"), gb.rename("core")], axis=1).fillna(0)
    out["type_region_pairs"] = len(J)
    out["type_region_totals_differing"] = int((J.ref != J.core).sum())
    if not out["identical_as_multiset"]:
        M = Cs.merge(Rs, on=key, how="outer", suffixes=("_core", "_ref")).fillna(0)
        out["cells_differing"] = int((M.n_core != M.n_ref).sum())
        out["sum_abs_diff_copies"] = int(np.abs(M.n_core - M.n_ref).sum())
    L = json.load(open(f"{core}/pack/libs.json"))
    tref = json.load(open(f"{ref}/engine/work/yc3d/thinning.json"))
    out["thinning_used"] = L.get("thinning"); out["thinning_ref"] = tref
    out["thinning_identical"] = L.get("thinning") == tref
    out["lib_order_identical"] = L["libs"] == [str(x) for x in np.load(f"{ref}/engine/work/yc3d/lib_ids.npy", allow_pickle=True)]
    out["pass"] = out["identical_in_order"] and out["thinning_identical"] and out["lib_order_identical"]
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--core", required=True, help="build dir of the core run (pack/ inside)")
    ap.add_argument("--ref", required=True, help="yeast port build (engine/work/yc3d inside)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--skip", default="", help="comma list of G4,G5,P3 to skip")
    a = ap.parse_args(argv)
    skip = set(a.skip.split(",")) - {""}
    res = {}
    if "G4" not in skip:
        res["G4"] = g4(a.core, a.ref); print("G4", json.dumps(res["G4"]), flush=True)
    if "G5" not in skip:
        res["G5"] = g5(a.core, a.ref, a.workers); print("G5", json.dumps(res["G5"]), flush=True)
    if "P3" not in skip and os.path.exists(f"{a.core}/pack/counts.npz"):
        res["P3"] = p3(a.core, a.ref); print("P3", json.dumps(res["P3"]), flush=True)
    res["pass"] = all(v["pass"] for v in res.values() if isinstance(v, dict))
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)
    print("A2 check", "PASS (level E)" if res["pass"] else "NOT IDENTICAL")
    return 0 if res["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
