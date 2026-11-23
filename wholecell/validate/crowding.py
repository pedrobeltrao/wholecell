# origin: jurkat/validate/crowding.py (legacy workspace, 2026-11-15)
# Changed: pack files, library folder and region list from the build / config (validate.crowding: regions, n_boxes,
# box, vox, probes, seeds); numerics unchanged. Seeds: the legacy seeds (box choice 1, calibration 0, walk seed = tile id)
# are used when validate.crowding.legacy_seeds is true, otherwise named wholecell.seeds streams. Output to
# <build>/validate/crowding.{json,tsv} (no live files). Fallback box choice (no fully filled tile) uses a stable sort
# (legacy quicksort tie order is CPU-dependent).
"""Crowding and tracer diffusion predicted from the packed cell.

For sample boxes (BOX nm cubes at tile centres that lie fully inside one region):
  1. occupancy at residue level on a VOX nm grid: every placed instance (random-packed tiles, pre-pass polymers, mRNA
     beads) is expanded to its residue/nucleotide beads (library beads rotated by the instance quaternion, w-x-y-z);
     a voxel is occupied if within R_BEAD of a bead. R_BEAD is calibrated so that the bead union of a library sample
     has the protein's volume (1.21 nm3/kDa), i.e. the occupied fraction phi is a macromolecular volume fraction.
  2. probe-accessible volume fraction f(a): voxels farther than a from any occupied voxel.
  3. obstruction-limited long-time diffusion D_obs/D0 for a probe of radius a: lattice random walk in the accessible space.
  4. hydrodynamic factor: Tokuyama-Oppenheim short-time self-diffusion for hard spheres at the box's phi.
  Predicted D/D0 = (D_obs/D0) x (D_hyd/D0); an upper bound for an inert tracer.
"""
from __future__ import annotations

import json
import os
import time

import numba as nb
import numpy as np
import pandas as pd
from scipy import ndimage

from . import common

_beads = {}
_LIB = [None]


def lib_beads(lid):
    if lid not in _beads:
        z = np.load(os.path.join(_LIB[0], f"{lid}.npz")); _beads[lid] = z["beads"].astype(np.float32)
    return _beads[lid]


def qmat(q):
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)], [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def calibrate_r(names, mass_kda, vox, n=40, seed=0):
    """bead radius r such that the bead union volume matches 1.21 nm3/kDa (median over a library sample)"""
    rng = np.random.default_rng(seed); idx = rng.choice(len(names), min(n, len(names)), replace=False)
    rs = np.arange(0.0, 1.5, 0.05); ratios = []
    for i in idx:
        b = lib_beads(names[i]); lo = b.min(0) - 2; dims = np.ceil((b.max(0) + 2 - lo) / vox).astype(int) + 1
        if np.prod(dims) > 4e7: continue
        g = np.ones(dims, bool); ij = np.floor((b - lo) / vox).astype(int); g[ij[:, 0], ij[:, 1], ij[:, 2]] = False
        d = ndimage.distance_transform_edt(g) * vox
        ratios.append([(d <= r).sum() * vox ** 3 / (mass_kda[i] * 1.21) for r in rs])
    med = np.median(np.array(ratios), 0)
    return float(np.interp(1.0, med, rs)), dict(zip(np.round(rs, 3).tolist(), np.round(med, 3).tolist()))


def gather(box_lo, A, key, L, PL, P_):
    """world-space beads of every instance near the box"""
    BOX, MARGIN, TILE = P_["box"], P_["margin"], P_["tile"]
    lo, hi = box_lo - MARGIN, box_lo + BOX + MARGIN
    pts = []
    t = np.floor((box_lo + BOX / 2 - A["lo0"]) / TILE).astype(int)
    for d in np.stack(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"), -1).reshape(-1, 3):
        j = key.get(tuple(t + d))
        if j is None: continue
        srcs = []
        f = os.path.join(A["tiles_dir"], f"{j}.npz")
        if os.path.exists(f):
            z = np.load(f); srcs.append((z["types"], z["pos"], z["quat"], L))
        if A["poff"] is not None:
            p0, p1 = A["poff"][j], A["poff"][j + 1]
            if p1 > p0:
                srcs.append((A["ptypes"][p0:p1], A["ppos"][p0:p1], A["pquat"][p0:p1], PL))
        for ty, pos, qu, names in srcs:
            pos = np.asarray(pos); k = np.where(np.all((pos > lo - 30) & (pos < hi + 30), axis=1))[0]
            for i in k:
                b = lib_beads(names[ty[i]]) @ qmat(np.asarray(qu[i], float)).T + pos[i]
                b = b[np.all((b > lo) & (b < hi), axis=1)]
                if len(b): pts.append(b)
        if A["moff"] is not None:
            m0, m1 = A["moff"][j], A["moff"][j + 1]
            if m1 > m0:
                mp = np.asarray(A["mpts"][m0:m1], np.float32); mp = mp[np.all((mp > lo) & (mp < hi), axis=1)]
                if len(mp): pts.append(("mrna", mp))
    prot = [p for p in pts if not isinstance(p, tuple)]; mr = [p[1] for p in pts if isinstance(p, tuple)]
    return (np.concatenate(prot) if prot else np.zeros((0, 3), np.float32)), (np.concatenate(mr) if mr else np.zeros((0, 3), np.float32))


@nb.njit(cache=True)
def walk(acc, starts, nsteps, rec_every, seed):
    np.random.seed(seed)
    nx, ny, nz = acc.shape; nw = starts.shape[0]; nrec = nsteps // rec_every
    pos = starts.copy(); disp = np.zeros((nw, 3), np.int64); msd = np.zeros(nrec)
    D = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]])
    for s in range(nsteps):
        for w in range(nw):
            k = np.random.randint(6)
            x = pos[w, 0] + D[k, 0]; y = pos[w, 1] + D[k, 1]; z = pos[w, 2] + D[k, 2]
            if x < 0 or y < 0 or z < 0 or x >= nx or y >= ny or z >= nz: continue
            if acc[x, y, z]:
                pos[w, 0] = x; pos[w, 1] = y; pos[w, 2] = z
                disp[w, 0] += D[k, 0]; disp[w, 1] += D[k, 1]; disp[w, 2] += D[k, 2]
        if (s + 1) % rec_every == 0:
            m = 0.0
            for w in range(nw): m += disp[w, 0] ** 2 + disp[w, 1] ** 2 + disp[w, 2] ** 2
            msd[(s + 1) // rec_every - 1] = m / nw
    return msd


def tokuyama(phi):
    b = np.sqrt(9 * phi / 8); c = 11 * phi / 16
    H = 2 * b * b / (1 - b) - c / (1 + 2 * c) - b * c * (2 + c) / ((1 + c) * (1 - b + c))
    return 1.0 / (1.0 + H)


def run(cfg, bdir, inp: common.Inputs, out_json, out_tsv, overwrite=False, n_boxes=None):
    from .. import seeds
    cc = (cfg.get("validate", {}) or {}).get("crowding", {}) or {}
    _LIB[0] = inp.lib_dir()
    meta = json.load(open(inp["tiles_meta"]))
    P_ = dict(tile=float(meta.get("TILE", meta.get("tile", 500.0))), box=float(cc.get("box", 150.0)), vox=float(cc.get("vox", 0.5)),
              margin=float(cc.get("margin", 25.0)))
    TILE, BOX, VOX = P_["tile"], P_["box"], P_["vox"]
    PROBES = [float(x) for x in cc.get("probes", [0.0, 1.0, 1.6, 2.3, 3.0])]
    nbox = int(n_boxes or cc.get("n_boxes", 4))
    legacy = bool(cc.get("legacy_seeds", False))
    M = pd.read_parquet(inp["tiles_parquet"])
    REGIONS = [r for r in cc.get("regions", []) if r in M.columns]
    if not REGIONS:
        vols = [c for c in M.columns if c not in ("tile", "ti", "tj", "tk") and not c.startswith("mem_")]
        REGIONS = vols[:1]
    L = json.load(open(inp["libs"])); names = L["libs"]
    mass = common.lib_masses(inp, bdir)[0] / 1e3          # libs.json mass_da, else library bead mass (yeast; A8 phase 2)
    has_pre = inp.has("preplaced_off")
    has_mr = inp.has("mrna_off")
    A = dict(lo0=np.array(meta["lo0"]), tiles_dir=inp["tiles_dir"],
             poff=np.load(inp["preplaced_off"]) if has_pre else None,
             ptypes=np.load(inp["preplaced_types"], mmap_mode="r") if has_pre else None,
             ppos=np.load(inp["preplaced_pos"], mmap_mode="r") if has_pre else None,
             pquat=np.load(inp["preplaced_quat"], mmap_mode="r") if has_pre else None,
             moff=np.load(inp["mrna_off"]) if has_mr else None, mpts=np.load(inp["mrna_pts"], mmap_mode="r") if has_mr else None)
    PL = json.load(open(inp["preplaced_libs"])) if has_pre else []
    mrad = json.load(open(inp["mrna_meta"]))["bead_radius"] if inp.has("mrna_meta") else 1.0
    key = {tuple(t): i for i, t in enumerate(M[["ti", "tj", "tk"]].values)}
    cal_seed = 0 if legacy else seeds.seed_int(cfg, "validate.crowding.calibrate")
    r_bead, curve = calibrate_r(names, mass, VOX, seed=cal_seed)
    rng = np.random.default_rng(1 if legacy else seeds.seed_int(cfg, "validate.crowding.boxes"))
    rows = []
    pad = 10
    for reg in REGIONS:
        full = np.where(M[reg].values >= 0.98 * TILE ** 3)[0]
        if len(full) == 0:
            # legacy: np.argsort(-x) (quicksort, SIMD-dispatched): the order of tied tiles depends on the CPU, so the legacy
            # box sample is not reproducible across machines. Port: stable sort (deterministic).
            full = np.argsort(-M[reg].values, kind="stable")[:50]
        for ti in rng.choice(full, min(nbox, len(full)), replace=False):
            t0 = time.time()
            box_lo = M.loc[ti, ["ti", "tj", "tk"]].values * TILE + A["lo0"] + (TILE - BOX) / 2
            P, Rm = gather(box_lo, A, key, L["libs"], PL, P_)
            n = int(round((BOX + 2 * pad) / VOX)); glo = box_lo - pad
            g = np.ones((n, n, n), bool)
            for pts in (P, Rm):
                ij = np.floor((pts - glo) / VOX).astype(int); ok = np.all((ij >= 0) & (ij < n), axis=1); ij = ij[ok]
                g[ij[:, 0], ij[:, 1], ij[:, 2]] = False
            d_b = ndimage.distance_transform_edt(g) * VOX
            occ = d_b <= r_bead
            if len(Rm):
                gm = np.ones((n, n, n), bool); ij = np.floor((Rm - glo) / VOX).astype(int); ok = np.all((ij >= 0) & (ij < n), axis=1)
                gm[ij[ok, 0], ij[ok, 1], ij[ok, 2]] = False; occ |= ndimage.distance_transform_edt(gm) * VOX <= mrad
            c = slice(int(pad / VOX), int((pad + BOX) / VOX)); phi = float(occ[c, c, c].mean())
            d_o = ndimage.distance_transform_edt(~occ) * VOX
            row = dict(region=reg, tile=int(ti), box_lo=box_lo.round(1).tolist(), phi=round(phi, 4), n_beads=int(len(P)),
                       D_hyd=round(float(tokuyama(phi)), 3))
            for a in PROBES:
                acc = d_o > a; f = float(acc[c, c, c].mean()); row[f"f_acc_{a}"] = round(f, 4)
                cand = np.argwhere(acc[c, c, c]) + int(pad / VOX)
                if len(cand) < 50: row[f"Dobs_{a}"] = 0.0; continue
                st = cand[rng.choice(len(cand), 1500, replace=len(cand) < 1500)].astype(np.int64)
                nsteps, rec = int(cc.get("walk_steps", 16000)), int(cc.get("walk_record", 400))
                wseed = int(ti) if legacy else seeds.seed_int(cfg, "validate.crowding.walk", int(ti))
                msd = walk(acc, st, nsteps, rec, wseed)
                tt = np.arange(1, len(msd) + 1) * rec; h = len(msd) // 2
                slope = np.polyfit(tt[h:], msd[h:], 1)[0]
                row[f"Dobs_{a}"] = round(float(slope), 3); row[f"D_{a}"] = round(float(slope * row["D_hyd"]), 3)
            row["s"] = round(time.time() - t0, 1); rows.append(row); print(row, flush=True)
    T = pd.DataFrame(rows)
    common.write_df(out_tsv, T, overwrite, index=False)
    S = T.groupby("region")[[c for c in T.columns if c.startswith(("phi", "f_acc", "Dobs", "D_"))]].mean().round(3)
    out = dict(r_bead=r_bead, calib=curve, summary=S.reset_index().to_dict("records"),
               settings=dict(P_, probes=PROBES, n_boxes=nbox, regions=REGIONS, legacy_seeds=legacy, lib_dir=_LIB[0]))
    common.write_json(out_json, out, overwrite)
    return out
