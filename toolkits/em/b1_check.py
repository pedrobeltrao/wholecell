"""Criterion B1 (docs/reference/history/v7_reproduction_criteria.md): an evaluator's G vs legacy cell3d on jgeom_cell1_v7.npz.

run(Gn, pores_n, data_root, ..., eng=<evaluator module>) -> dict; eng = sdf_sections_v1 (default) or cellgeom_eval_v1
(the table-driven evaluator; Gn must then carry its rule tables, e.g. from cellgeom.load). Points (seeded, wholecell.seeds stream 'geometry.b1' of master 20261115):
n_uniform uniform in the legacy cell-field box, n_near within 50 nm of a membrane (rejection on the legacy membrane
SDFs; the ER field clamps at 48 nm, clamped values are not 'near'), n_pore within 120 nm of pore centres. Compared:
engine labels (int8) and the 8 membrane SDFs in float64. Legacy code is imported read-only: no bytecode written and
the numba cache must be redirected (env NUMBA_CACHE_DIR, set before numba is imported; tests/check_b1_geometry.py does).
"""
import os
import sys
import time

import numba as nb
import numpy as np

from wholecell import seeds
from wholecell.geometry.engines import sdf_sections_v1 as e

LEGACY_NPZ = "{data_root}/data/derived/jgeom_cell1_v7.npz"


def legacy_module(data_root):
    if not os.environ.get("NUMBA_CACHE_DIR"):
        raise RuntimeError("set NUMBA_CACHE_DIR before importing legacy cell3d (its numba cache must not be written "
                           "into the legacy tree)")
    sys.dont_write_bytecode = True
    sys.path.insert(0, os.path.join(data_root, "jurkat", "yc", "pipeline"))
    try:
        import cell3d as lc3
    finally:
        sys.path.pop(0)
    return lc3


def run(Gn, pores_n, data_root, legacy=LEGACY_NPZ, n_uniform=2_000_000, n_near=500_000, n_pore=50_000,
        threads=8, source="", log=print, eng=e):
    nb.set_num_threads(threads)
    lc3 = legacy_module(data_root)
    lpath = legacy.format(data_root=data_root)
    z = np.load(lpath)
    Gl = z["G"]; pores_l = z["pores"]
    log(f"B1: legacy G {Gl.shape[0]:,} new G {Gn.shape[0]:,} ({source})")
    if not np.array_equal(pores_l, pores_n):
        return {"criterion": "B1", "pass": False, "error": "pores differ"}
    E = np.zeros(0, np.int64); PQ = np.zeros((0, 8))

    @nb.njit(parallel=True)
    def legacy_eval(P, G, pores):
        n = P.shape[0]; lab = np.zeros(n, np.int8); M = np.zeros((n, 8))
        for i in nb.prange(n):
            out = np.zeros(8)
            mat = lc3.mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, E, E, PQ, E, out)
            lab[i] = lc3.label_from(P[i, 0], P[i, 1], P[i, 2], G, pores, out, mat)
            for m in range(8):
                M[i, m] = out[m]
        return lab, M

    new_mem, new_label = eng.mem_sdfs, eng.label_from

    @nb.njit(parallel=True)
    def new_eval(P, G, pores):
        n = P.shape[0]; lab = np.zeros(n, np.int8); M = np.zeros((n, 8))
        for i in nb.prange(n):
            out = np.zeros(8)
            mat = new_mem(P[i, 0], P[i, 1], P[i, 2], G, pores, E, E, PQ, E, out)
            lab[i] = new_label(P[i, 0], P[i, 1], P[i, 2], G, pores, out, mat)
            for m in range(8):
                M[i, m] = out[m]
        return lab, M

    rng = seeds.rng({"seeds": {"master": 20261115}}, "geometry.b1")
    o = 9  # legacy cell field header
    lo = Gl[o:o + 3]; hi = lo + Gl[o + 3] * Gl[o + 7] * Gl[o + 4:o + 7]
    U = lo + (hi - lo) * rng.random((n_uniform, 3))
    oe = 9 + 12 * lc3.F_ER; er_clamp = min(Gl[oe + 8], Gl[oe + 9])
    near = []; nn = 0; drawn = 0
    while nn < n_near:
        C = lo + (hi - lo) * rng.random((2_000_000, 3)); drawn += len(C)
        _, M = legacy_eval(C, Gl, pores_l)
        A = np.abs(M[:, :7]); A[:, lc3.M_ER] = np.where(A[:, lc3.M_ER] < er_clamp, A[:, lc3.M_ER], 1e9)
        k = C[A.min(1) < 50.0]
        near.append(k); nn += len(k)
    Nr = np.concatenate(near)[:n_near]
    pi = rng.integers(0, len(pores_l), n_pore)
    d = rng.normal(size=(n_pore, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    Pp = pores_l[pi, :3] + d * (120.0 * rng.random(n_pore) ** (1 / 3))[:, None]
    res = {"criterion": "B1", "legacy": legacy, "new_source": source, "engine": eng.ENGINE_NAME,
           "seed_stream": "geometry.b1/master 20261115", "box_lo": lo.tolist(), "box_hi": hi.tolist(),
           "near_candidates_drawn": drawn, "sets": {}}
    allok = True
    for name, P in (("uniform", U), ("near_membrane_50nm", Nr), ("near_pore_120nm", Pp)):
        t1 = time.time()
        l0, m0 = legacy_eval(P, Gl, pores_l)
        l1, m1 = new_eval(P, Gn, pores_n)
        dm = np.abs(m0 - m1)
        nlab = int((l0 != l1).sum())
        r = {"n": int(len(P)), "label_mismatch": nlab, "label_identical_frac": 1 - nlab / len(P),
             "mem_sdf_max_absdiff_nm": float(dm.max()), "mem_sdf_exact_equal": bool(np.array_equal(m0, m1)),
             "label_counts": {lc3.LABEL_NAMES[k]: int((l0 == k).sum()) for k in range(len(lc3.LABEL_NAMES))},
             "seconds": round(time.time() - t1, 1)}
        r["pass"] = nlab == 0 and r["mem_sdf_max_absdiff_nm"] < 1e-9
        allok &= r["pass"]
        res["sets"][name] = r
        log(name + " " + str({k: v for k, v in r.items() if k != "label_counts"}))
    res["pass"] = bool(allok)
    res["n_total"] = int(sum(r["n"] for r in res["sets"].values()))
    log(f"B1 {'PASS' if allok else 'FAIL'} ({res['n_total']:,} points)")
    return res
