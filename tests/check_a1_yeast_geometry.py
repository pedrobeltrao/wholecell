"""A1 (yeast half): cellgeom_eval_v1 + config/geometry/rules/yeast_cell3d_v1.yaml vs the vendored yeast cell3d on a
yeast cellgeom_v1 package (level E: labels, the 10 membrane SDFs and the packer clearance identical).

The vendored evaluator runs in the yeast engine env, the core evaluator in vcell (and optionally in the yeast env too,
to separate numba-version effects); points and results go through npz files in --out:

    python tests/check_a1_yeast_geometry.py points   --package P --out D [--n-uniform 2000000 --n-near 500000 --n-pore 50000]
    $YEAST_PY tests/check_a1_yeast_geometry.py vendored --package P --out D           # -> D/vendored.npz (reference)
    python tests/check_a1_yeast_geometry.py vendored --package P --out D --tag vcell  # -> D/vendored_vcell.npz (optional)
    python tests/check_a1_yeast_geometry.py core     --package P --out D [--tag vcell] # -> D/core_<tag>.npz
    python tests/check_a1_yeast_geometry.py compare  --out D --pairs core_vcell:vendored,core_yeastenv:vendored
                                                                                      # -> D/result.json

points: seeded (numpy default_rng(--seed)): uniform in the enclosing cube [-R, R]^3 (R = enclosing radius), points
within 50 nm of any membrane (|SDF| < 50 for some membrane, candidates uniform in the cube), points near pores
(centre + N(0, 40 nm)), and box tests (--n-boxes boxes of 600 nm with the vendored Cell.select prim / pore subsets).
"""
import argparse
import json
import os
import sys
import tempfile
import time

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
sys.dont_write_bytecode = True
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

SKIPS = (-1, 0, 1, 2, 4, 9)        # clearance: no skip, PM, ER, NE, MITO_OM, CELL_OUTER


def _core(pkg):
    from wholecell.geometry import cellgeom
    return cellgeom.load(pkg)


def cmd_points(a):
    cell = _core(a.package)
    R = float(cell.manifest["domain"]["enclosing_radius_nm"])
    rng = np.random.default_rng(a.seed)
    U = rng.uniform(-R, R, (a.n_uniform, 3))
    near = []
    n = 0
    while n < a.n_near:
        C = rng.uniform(-R, R, (2_000_000, 3))
        M = cell.mems64(C)
        C = C[(np.abs(M) < 50.0).any(1)]
        near.append(C); n += len(C)
    N = np.concatenate(near)[:a.n_near]
    k = rng.integers(0, len(cell.pores), a.n_pore)
    Q = cell.pores[k, :3] + rng.normal(0.0, 40.0, (a.n_pore, 3))
    boxes = []
    for b in range(a.n_boxes):
        c = rng.uniform(-0.8 * R, 0.8 * R, 3)
        boxes.append(c + rng.uniform(-300.0, 300.0, (a.n_box_points, 3)))
    np.savez(os.path.join(a.out, "points.npz"), uniform=U, near=N, pore=Q, boxes=np.stack(boxes))
    print(f"points: uniform {len(U):,} near {len(N):,} pore {len(Q):,} boxes {a.n_boxes} x {a.n_box_points}")


def cmd_vendored(a):
    import numba as nb
    from plugins.yeast import geom_engine
    cell = geom_engine.load(a.package)
    import cell3d as c3

    @nb.njit(parallel=True)
    def mems(P, G, pores, pidx, PT, PQ, PC):
        g = np.zeros((P.shape[0], c3.NMEM))
        for i in nb.prange(P.shape[0]):
            out = np.zeros(c3.NMEM)
            c3.mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
            for m in range(c3.NMEM):
                g[i, m] = out[m]
        return g

    @nb.njit(parallel=True)
    def clear(P, G, pores, pidx, PT, PQ, PC, skip):
        g = np.zeros(P.shape[0])
        for i in nb.prange(P.shape[0]):
            out = np.zeros(c3.NMEM)
            c3.mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
            g[i] = c3.min_mem_dist(out, skip)
        return g

    Z = np.load(os.path.join(a.out, "points.npz"))
    allp = np.arange(len(cell.PT), dtype=np.int64)
    res = {}
    t = time.time()
    for key in ("uniform", "near", "pore"):
        P = Z[key]
        res[f"{key}_lab"] = c3.label_points(P, cell.G, cell.pores, allp, cell.PT, cell.PQ, cell.PC)
        res[f"{key}_mem"] = mems(P, cell.G, cell.pores, allp, cell.PT, cell.PQ, cell.PC)
        for s in SKIPS:
            res[f"{key}_clear{s}"] = clear(P, cell.G, cell.pores, allp, cell.PT, cell.PQ, cell.PC, s)
    pidxs, pores, blab, bmem = [], [], [], []
    for B in Z["boxes"]:
        pidx, pr = cell.select(B.min(0), B.max(0))
        pidxs.append(pidx); pores.append(pr)
        blab.append(c3.label_points(B, cell.G, pr, pidx, cell.PT, cell.PQ, cell.PC))
        bmem.append(mems(B, cell.G, pr, pidx, cell.PT, cell.PQ, cell.PC))
    res["box_lab"] = np.stack(blab); res["box_mem"] = np.stack(bmem)
    res["box_pidx"] = np.array(pidxs, dtype=object); res["box_pores"] = np.array(pores, dtype=object)
    np.savez(os.path.join(a.out, "vendored.npz" if a.tag == "ref" else f"vendored_{a.tag}.npz"), **res)
    print(f"vendored cell3d done in {time.time() - t:.0f} s (numba {nb.__version__})")


def cmd_core(a):
    import numba as nb
    cell = _core(a.package)
    ev = cell.engine
    Z = np.load(os.path.join(a.out, "points.npz"))
    V = np.load(os.path.join(a.out, "vendored.npz"), allow_pickle=True)
    res = {}
    t = time.time()
    for key in ("uniform", "near", "pore"):
        P = Z[key]
        res[f"{key}_lab"] = cell.labels(P)
        res[f"{key}_mem"] = cell.mems64(P)
        res[f"{key}_mem32"] = cell.mems(P)
        for s in SKIPS:
            res[f"{key}_clear{s}"] = cell.clearance(P, s)
    blab, bmem = [], []
    for B, pidx, pr in zip(Z["boxes"], V["box_pidx"], V["box_pores"]):
        pidx = np.asarray(pidx, np.int64); pr = np.ascontiguousarray(pr, np.float64).reshape(-1, 6)
        blab.append(ev.label_points(B, cell.G, pr, pidx, cell.PT, cell.PQ, cell.PC))
        bmem.append(ev.mem_points64(B, cell.G, pr, pidx, cell.PT, cell.PQ, cell.PC))
    res["box_lab"] = np.stack(blab); res["box_mem"] = np.stack(bmem)
    np.savez(os.path.join(a.out, f"core_{a.tag}.npz"), **res)
    print(f"core cellgeom_eval_v1 [{a.tag}] done in {time.time() - t:.0f} s (numba {nb.__version__})")


def compare(V, C):
    """level-E comparison of a result C with a reference V (labels, membrane SDFs, clearance identical)"""
    out = {}
    for key in ("uniform", "near", "pore", "box"):
        lv, lc = V[f"{key}_lab"], C[f"{key}_lab"]
        mv, mc = V[f"{key}_mem"].reshape(-1, 10), C[f"{key}_mem"].reshape(-1, 10)
        d = np.abs(mv - mc)
        r = dict(n=int(lv.size), label_mismatches=int((lv != lc).sum()),
                 mem_identical=bool(np.array_equal(mv, mc)),
                 max_abs_dsdf=[float(x) for x in d.max(0)],
                 n_mem_values_differing=[int(x) for x in (mv != mc).sum(0)])
        if f"{key}_mem32" in C:
            r["mem_points_float32_identical"] = bool(np.array_equal(mv.astype(np.float32), C[f"{key}_mem32"]))
        if key != "box":
            r["clearance_identical"] = {str(s): bool(np.array_equal(V[f"{key}_clear{s}"], C[f"{key}_clear{s}"]))
                                        for s in SKIPS}
            r["label_counts"] = {int(k): int(v) for k, v in zip(*np.unique(lv, return_counts=True))}
        out[key] = r
    ok = all(r["label_mismatches"] == 0 and r["mem_identical"] and all(r.get("clearance_identical", {1: True}).values())
             and r.get("mem_points_float32_identical", True) for r in out.values())
    return ok, out


def cmd_compare(a):
    res = {"membranes": ["PM", "ER", "NE", "VAC", "MITO_OM", "MITO_IM", "GOLGI", "PEROX", "LD", "CELL_OUTER"]}
    allok = True
    for tag in a.pairs.split(","):
        cf, vf = tag.split(":")
        V = np.load(os.path.join(a.out, f"{vf}.npz"), allow_pickle=True)
        ok, r = compare(V, np.load(os.path.join(a.out, f"{cf}.npz")))
        res[tag] = dict(pass_=ok, **r)
        allok &= ok
        print(f"[{tag}] level E: {ok}; " + "; ".join(
            f"{k}: {v['label_mismatches']} label mismatches, max|dSDF| {max(v['max_abs_dsdf']):.3g}" for k, v in r.items()))
    res["pass"] = allok
    json.dump(res, open(os.path.join(a.out, "result.json"), "w"), indent=1)
    return 0 if allok else 1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["points", "vendored", "core", "compare"])
    ap.add_argument("--package")
    ap.add_argument("--out", required=True)
    ap.add_argument("--tag", default=None, help="core: result tag (default vcell); vendored: default ref")
    ap.add_argument("--pairs", default="core_vcell:vendored")
    ap.add_argument("--seed", type=int, default=20261117)
    ap.add_argument("--n-uniform", type=int, default=2_000_000)
    ap.add_argument("--n-near", type=int, default=500_000)
    ap.add_argument("--n-pore", type=int, default=50_000)
    ap.add_argument("--n-boxes", type=int, default=20)
    ap.add_argument("--n-box-points", type=int, default=5_000)
    a = ap.parse_args(argv)
    a.tag = a.tag or ("ref" if a.cmd == "vendored" else "vcell")
    os.makedirs(a.out, exist_ok=True)
    return {"points": cmd_points, "vendored": cmd_vendored, "core": cmd_core, "compare": cmd_compare}[a.cmd](a) or 0


if __name__ == "__main__":
    sys.exit(main())
