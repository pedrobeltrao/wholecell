# new in wholecell (B4 J6 fix, 2026-11-18)
"""B4 J6 (ATP synthase dimer rows on the crista rims) after the dense rim sampling: measurements for
docs/unification_criteria/B4_preplaced_addendum_J6_2026-11-18.md, on the module files of a Jurkat serial pre-placed run.

  dimers, monomers, share of the CPX-6151 pool; density per um of rim with two rim-length measures: the rim lines the
  builder traced (its log line "traced rim length") and the B6b band estimate of the whole cell (4,040 um,
  B6b_crista_rims_addendum_2026-11-17.md R6); spacing of neighbouring dimers along the rows;
  overlaps: sphere trees of the ATP monomers (library CPX-6151) between different dimers and within a dimer, and
  against every other pre-placed layer of the run near the dimers (instance layers: library sphere trees; bead layers
  mrna / mtdna: their beads), overlap = d < r_i + r_j - tol for tol 0 / 0.3 (pack TOL) / 0.5 nm (J5);
  membrane distances on a seeded sample of dimers: |IM midsurface distance| of the transmembrane beads (|z| < 2 nm in
  the library frame), the per-dimer maximum, the share of every-third soluble beads within 0.6 nm of a bilayer and
  outside the matrix / IMS (the site test of the builder, re-evaluated independently in Python).

  python scripts/b4_j6_atp_check.py --serial <run>/preplaced/serial --build <bdir: library, inventory> \
         --log <stage_serial.log> --out <json> [--sample 4000]
NUMBA_NUM_THREADS / OMP_NUM_THREADS set the cores (3 for checks).
"""
import argparse
import json
import os
import re
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
# float32 export rounding of the pose (B4_preplaced_addendum_J6_2026-11-18.md §3: the one TM exception is 4.5000005 nm,
# the one soluble-bead exception 3e-5 nm inside 0.6 nm; positions up to ~10 um have a float32 step <= 1e-3 nm)
F32_EPS_NM = 1e-3
RIM_BAND_UM = 4040.0          # B6b_crista_rims_addendum_2026-11-17.md R6 (band method, em_cell1_v8r)
LIT = (50.0, 83.0)            # dimers per um of rim: Strauss et al. 2008 EMBO J; Davies et al. 2011 PNAS


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _rot(q):
    from scipy.spatial.transform import Rotation
    q = np.asarray(q, np.float64)
    return Rotation.from_quat(np.c_[q[:, 1:], q[:, :1]]).as_matrix()


def _spheres(R, pos, stc, str_):
    C = (np.einsum("nij,sj->nsi", R, stc) + pos[:, None, :]).reshape(-1, 3)
    return C, np.tile(str_, len(R)), np.repeat(np.arange(len(R)), len(str_))


def _owner_pairs(CA, RA, OA, CB, RB, OB, tols, rmax, chunk=400_000, tree=None):
    """instance pairs (owner of an A sphere, owner of a B sphere) with a sphere pair d < ra + rb - tol, per tol
    (unique rows; A x A: pass the same arrays, pairs with oa < ob kept). Chunked over A (memory)."""
    from scipy.spatial import cKDTree
    tb = tree if tree is not None else cKDTree(CB)
    same = CA is CB
    acc = {t: [] for t in tols}
    for c0 in range(0, len(CA), chunk):
        sl = slice(c0, min(c0 + chunk, len(CA)))
        lst = tb.query_ball_point(CA[sl], RA[sl] + rmax)
        ia = np.repeat(np.arange(sl.start, sl.stop), [len(x) for x in lst])
        ib = np.fromiter((j for x in lst for j in x), np.int64, len(ia))
        del lst
        if not len(ia):
            continue
        oa, ob = OA[ia], OB[ib]
        k = oa < ob if same else np.ones(len(ia), bool)
        ia, ib, oa, ob = ia[k], ib[k], oa[k], ob[k]
        d = np.linalg.norm(CA[ia] - CB[ib], axis=1); sr = RA[ia] + RB[ib]
        for t in tols:
            m = d < sr - t
            if m.any():
                acc[t].append(np.unique(np.c_[oa[m], ob[m]], axis=0))
    return {t: (np.unique(np.concatenate(v), axis=0) if v else np.zeros((0, 2), np.int64)) for t, v in acc.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--serial", required=True); ap.add_argument("--build", required=True); ap.add_argument("--log")
    ap.add_argument("--out", required=True); ap.add_argument("--sample", type=int, default=4000)
    ap.add_argument("--skip-overlaps", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    import pandas as pd
    from scipy.spatial import cKDTree
    from wholecell import config, registry
    from wholecell.preplaced.context import context
    sys.path[:0] = [os.path.join(REPO, "scripts", "dev"), os.path.join(REPO, "scripts")]
    from b4_preplaced_check import _lib_trees
    cfg = config.load("jurkat"); registry.load_plugins(cfg)
    LIB = _lib_trees(a.build)
    mods = os.path.join(a.serial, "modules")
    z = np.load(os.path.join(mods, "atp_synthase.npz"))
    mono = str(z["lib_ids"][0]); pos = z["pos"].astype(np.float64); Rm = _rot(z["quat"]); dim = z["dimer"]
    I = pd.read_csv(os.path.join(a.build, "inventory", "pack_table.tsv"), sep="\t")
    pool = float(I.loc[I.ing_id == mono, "copies"].iloc[0])
    nd = int(dim.max()) + 1 if len(dim) else 0
    res = dict(serial=a.serial, monomer=mono, monomers=int(len(pos)), dimers=nd, pool=pool,
               monomer_share_of_pool=round(len(pos) / pool, 4), dimers_le_pool_half=bool(0 < nd <= pool / 2))
    traced = None
    if a.log and os.path.exists(a.log):
        for line in open(a.log):
            m = re.search(r"rows\): (\d+) of (\d+) rim sites tried \((\d+) rim points on (\d+) rim lines, ([\d.]+) um traced", line)
            if m:
                traced = float(m.group(5))
                res["builder_log"] = dict(dimers=int(m.group(1)), sites_tried=int(m.group(2)), rim_points=int(m.group(3)),
                                          rim_lines=int(m.group(4)), traced_um=traced)
    res["density"] = dict(per_um_traced=round(nd / traced, 2) if traced else None, rim_traced_um=traced,
                          per_um_band=round(nd / RIM_BAND_UM, 2), rim_band_um=RIM_BAND_UM, literature_per_um=list(LIT))
    res["density"]["in_literature_traced"] = bool(traced and LIT[0] <= nd / traced <= LIT[1])
    res["density"]["in_literature_band"] = bool(LIT[0] <= nd / RIM_BAND_UM <= LIT[1])
    # dimer centres (mean of the two monomers) and spacing to the nearest dimer
    DC = np.zeros((nd, 3)); np.add.at(DC, dim, pos); DC /= np.maximum(np.bincount(dim, minlength=nd), 1)[:, None]
    if nd > 1:
        dd, _ = cKDTree(DC).query(DC, k=2)
        res["dimer_nn_nm"] = dict(median=round(float(np.median(dd[:, 1])), 2), p05=round(float(np.quantile(dd[:, 1], 0.05)), 2),
                                  p95=round(float(np.quantile(dd[:, 1], 0.95)), 2), min=round(float(dd[:, 1].min()), 2))
    log("dimers", nd, res["density"], res.get("dimer_nn_nm"))
    L = np.load(LIB[mono]); stc = L["st_c"].astype(np.float64); str_ = L["st_r"].astype(np.float64)
    if a.skip_overlaps and os.path.exists(a.out):          # keep the overlap results of an earlier run of this file
        old = json.load(open(a.out)); res.update({k: old[k] for k in ("overlap_atp_atp", "overlap_atp_other") if k in old})
    else:
        _overlaps(res, a, cfg, LIB, mods, pos, Rm, dim, stc, str_)

    # membrane distances (seeded sample of dimers)
    ctx = context(cfg, a.build, dict(libs=[mono]))
    K, cell = ctx.K, ctx.geom.cell
    Lb = L["beads"].astype(np.float64); tm = np.abs(Lb[:, 2]) < 2.0
    rng = np.random.default_rng(20261118)
    pick = np.sort(rng.choice(nd, min(a.sample, nd), replace=False)) if nd else np.zeros(0, int)
    keep_lab = [K.lab("mito_matrix"), K.lab("mito_IMS")]
    tm_dev, per_dim_max, sol_in, sol_bad = [], {}, {}, {}
    beads = _make_beads(K)
    args = ctx.cell.args()
    TM2 = np.r_[tm, tm]
    sol2 = ~TM2 & (np.arange(2 * len(Lb)) % 3 == 0)    # the builder's dimer bead array: copy 0 then copy 1, every third
    for d_ in pick:                                    # soluble residue of the concatenation
        ii = np.flatnonzero(dim == d_)
        X = np.concatenate([Lb @ Rm[i].T + pos[i] for i in ii])
        im, lab, mm = beads(X, *args)
        d = np.abs(im[TM2]); tm_dev.append(d); per_dim_max[d_] = float(d.max())
        sol_in[d_] = [int((mm[sol2] < 0.6).sum()), int(sol2.sum()), int((mm[sol2] < 0.6 - F32_EPS_NM).sum())]
        sol_bad[d_] = int((~np.isin(lab[sol2], keep_lab)).sum())
    tm_all = np.concatenate(tm_dev) if tm_dev else np.zeros(0)
    pdm = np.array(list(per_dim_max.values()))
    fr = np.array([v[0] / max(v[1], 1) for v in sol_in.values()])
    fr32 = np.array([v[2] / max(v[1], 1) for v in sol_in.values()])
    res["membrane"] = dict(
        sample_dimers=int(len(pick)), tm_beads=int(len(tm_all)),
        tm_abs_im_nm=dict(median=round(float(np.median(tm_all)), 3), p90=round(float(np.quantile(tm_all, 0.9)), 3),
                          p99=round(float(np.quantile(tm_all, 0.99)), 3), max=round(float(tm_all.max()), 3)) if len(tm_all) else None,
        per_dimer_max_tm_nm=dict(median=round(float(np.median(pdm)), 3), max=round(float(pdm.max()), 3),
                                 share_le={str(t): round(float((pdm <= t).mean()), 4) for t in (2.5, 3.0, 3.5, 4.5)},
                                 share_le_4_5_f32=round(float((pdm <= 4.5 + F32_EPS_NM).mean()), 6),
                                 n_over_4_5=int((pdm > 4.5).sum()), max_exact=float(pdm.max())) if len(pdm) else None,
        soluble_in_membrane_share=dict(median=round(float(np.median(fr)), 4), max=round(float(fr.max()), 4),
                                       dimers_over_2pct=int((fr > 0.02).sum()), dimers_over_2pct_f32=int((fr32 > 0.02).sum()),
                                       f32_eps_nm=F32_EPS_NM) if len(fr) else None,
        soluble_outside_matrix_ims=int(sum(sol_bad.values())))
    log("membrane", res["membrane"])
    res["seconds"] = round(time.time() - t0)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1, default=str)
    log("done", res["seconds"], "s")


def _overlaps(res, a, cfg, LIB, mods, pos, Rm, dim, stc, str_):
    from scipy.spatial import cKDTree
    from b4_preplaced_check import _layers_new
    CA, RA, OA = _spheres(Rm, pos, stc, str_)
    tols = (0.0, 0.3, 0.5)
    pr = _owner_pairs(CA, RA, OA, CA, RA, OA, tols, float(str_.max()))
    ov = {}
    for t, inst in pr.items():
        same = dim[inst[:, 0]] == dim[inst[:, 1]] if len(inst) else np.zeros(0, bool)
        ov[str(t)] = dict(monomer_pairs_between_dimers=int((~same).sum()), monomer_pairs_within_dimer=int(same.sum()))
    res["overlap_atp_atp"] = ov
    log("ATP-ATP overlaps", ov)
    # overlaps against the other pre-placed layers near the dimers
    lay = _layers_new(a.serial, cfg)
    bmax = float((np.linalg.norm(stc, axis=1) + str_).max())
    tree_d = cKDTree(pos)
    other = {}
    for name, (ids, types, P, Q) in lay.items():
        if name == "atp_synthase" or not len(P):
            continue
        P = np.asarray(P, np.float64)
        brs = []
        for s in ids:
            s = str(s)
            if s in LIB:
                q = np.load(LIB[s]); brs.append(float((np.linalg.norm(q["st_c"], axis=1) + q["st_r"]).max()))
            else:
                brs.append(np.nan)
        brs = np.array(brs)
        rr = float(np.nanmax(brs)) + bmax
        box = np.flatnonzero(((P >= pos.min(0) - rr) & (P <= pos.max(0) + rr)).all(1)) if len(pos) else np.zeros(0, np.int64)
        idx = box[tree_d.query_ball_point(P[box], rr, return_length=True) > 0] if len(box) else box
        row = dict(instances=int(len(P)), near_atp=int(len(idx)))
        if len(idx):
            Cs, Rs, Os = [], [], []
            for ti in np.unique(np.asarray(types)[idx]):
                sel = idx[np.asarray(types)[idx] == ti]
                q = np.load(LIB[str(ids[ti])])
                c, r, o = _spheres(_rot(np.asarray(Q)[sel]), P[sel], q["st_c"].astype(np.float64), q["st_r"].astype(np.float64))
                Cs.append(c); Rs.append(r); Os.append(sel[o])
            CB, RB, OB = np.concatenate(Cs), np.concatenate(Rs), np.concatenate(Os)
            prb = _owner_pairs(CA, RA, OA, CB, RB, OB, tols, float(RB.max()))
            for t, pairs in prb.items():
                row[f"pairs_tol{t}"] = int(len(pairs))
        other[name] = row
    for name in ("mtdna", "mrna"):
        f = os.path.join(mods, name + ".npz")
        if not os.path.exists(f):
            continue
        q = np.load(f); B = q["pts"].astype(np.float64); r = float(q["bead_radius"])
        box = ((B >= pos.min(0) - bmax - r) & (B <= pos.max(0) + bmax + r)).all(1) if len(pos) else np.zeros(len(B), bool)
        B = B[box]
        tbb = cKDTree(B) if len(B) else None
        near = tbb.query_ball_point(pos, bmax + r, return_length=True) if tbb is not None else np.zeros(len(pos))
        row = dict(beads=int(box.size), beads_in_atp_box=int(len(B)), atp_monomers_near=int((near > 0).sum()))
        if (near > 0).any():
            sub = B[np.unique(np.concatenate(tbb.query_ball_point(pos[near > 0], bmax + r)))]
            pb = _owner_pairs(CA, RA, OA, sub, np.full(len(sub), r), np.arange(len(sub)), tols, r)
            for t, pairs in pb.items():
                row[f"atp_monomers_overlapping_tol{t}"] = int(len(np.unique(pairs[:, 0]))) if len(pairs) else 0
                row[f"bead_pairs_tol{t}"] = int(len(pairs))
        other[name] = row
    res["overlap_atp_other"] = other
    log("ATP-other overlaps", other)


def _make_beads(K):
    """per bead: IM midsurface distance, label, clearance to the nearest bilayer (K.mmd, as the builder's site test)"""
    import numba as nb
    c3 = K.c3; NMEM = int(K.nmem); M_IM = K.mem("MITO_IM"); mmd = K.mmd

    @nb.njit(cache=False)
    def beads(X, G, pores, pidx, PT, PQ, PC):
        out = np.zeros(NMEM); im = np.empty(X.shape[0]); lab = np.empty(X.shape[0], np.int64); mm = np.empty(X.shape[0])
        for k in range(X.shape[0]):
            mtx = c3.mem_sdfs(X[k, 0], X[k, 1], X[k, 2], G, pores, pidx, PT, PQ, PC, out)
            im[k] = out[M_IM]
            lab[k] = c3.label_from(X[k, 0], X[k, 1], X[k, 2], G, pores, out, mtx)
            mm[k] = mmd(out, -1, G)
        return im, lab, mm

    return beads


if __name__ == "__main__":
    main()
