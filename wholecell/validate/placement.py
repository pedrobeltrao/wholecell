"""Placement integrity on a seeded sample of tiles (new in wholecell; ideas from the yeast validation/placement_checks.py:
per-species placed share, shell-band test).

Sample: validate.placement.n_tiles tile ids (default 300) drawn from the named seed stream 'validate.placement.tiles'
(or validate.placement.seed), so every build of one config is checked on the same tiles.

  (a) clash rate      instance pairs whose sphere trees overlap by more than `clash_nm` (0.5 nm) within a tile, counting
                      pairs with at least one packed instance (packed-packed + packed-pre-pass), per 1e6 packed instances.
                      Packed-mRNA bead pairs and pre-pass-pre-pass pairs are reported separately (not in the rate).
  (b) containment     packed instances whose centre is not allowed for their region, per 1e6 packed instances:
                      volume / virtual regions: centre label in the allowed labels (virtual distance rules not re-tested);
                      membranes: |SDF of the membrane| at the centre <= membrane_tol_nm (packer projects to < 0.1 nm);
                      shells: allowed label and side*SDF in [d0, d1 (+ bounding radius)] (+- shell_tol_nm), and within
                      pore_lat of a pore axis (|axial| <= pore_axial_nm) where the shell has a pore limit.
  (c) orientation     membrane-region instances whose library +z axis (rotated by the instance quaternion, w-x-y-z) points
                      against the membrane SDF gradient (cytoplasmic / N side), per 1e6 membrane instances.
  species share       per region: placed / requested over all tiles, number of species (>= min_req requested) below 0.5
                      and 0.8 of requested (from the whole-cell job tables).
Output: <build>/validate/placement.json (+ placement_tiles.tsv per sampled tile)
"""
from __future__ import annotations

import json
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numba as nb
import numpy as np
import pandas as pd

from . import common, regions as R_
from .geom import label_code_map, open_geometry

_G = {}


@nb.njit(cache=True)
def _qmat(q):
    a, b, c, d = q[0], q[1], q[2], q[3]
    M = np.empty((3, 3))
    M[0, 0] = a*a+b*b-c*c-d*d; M[0, 1] = 2*(b*c-a*d); M[0, 2] = 2*(b*d+a*c)
    M[1, 0] = 2*(b*c+a*d); M[1, 1] = a*a-b*b+c*c-d*d; M[1, 2] = 2*(c*d-a*b)
    M[2, 0] = 2*(b*d-a*c); M[2, 1] = 2*(c*d+a*b); M[2, 2] = a*a-b*b-c*c+d*d
    return M


@nb.njit(cache=True)
def _spheres(types, pos, quat, st_off, STC, STR, inst0):
    n = 0
    for i in range(types.shape[0]):
        n += st_off[types[i] + 1] - st_off[types[i]]
    C = np.empty((n, 3)); R = np.empty(n); I = np.empty(n, np.int64)
    k = 0
    for i in range(types.shape[0]):
        M = _qmat(quat[i].astype(np.float64))
        a0 = st_off[types[i]]; a1 = st_off[types[i] + 1]
        for s in range(a0, a1):
            for a in range(3):
                C[k, a] = M[a, 0] * STC[s, 0] + M[a, 1] * STC[s, 1] + M[a, 2] * STC[s, 2] + pos[i, a]
            R[k] = STR[s]; I[k] = inst0 + i
            k += 1
    return C, R, I


@nb.njit(cache=True)
def _clash(C, R, I, cat, h, tol, cap):
    """sphere pairs of different instances overlapping by > tol; skips pairs where both instances have cat == 2 (mRNA)
    or both cat == 1 and skip11. Returns (inst_a, inst_b, depth) arrays (truncated at cap) and the total count."""
    n = C.shape[0]
    out_a = np.empty(cap, np.int64); out_b = np.empty(cap, np.int64); out_d = np.empty(cap)
    if n == 0:
        return out_a[:0], out_b[:0], out_d[:0], 0
    lo = np.empty(3); dims = np.empty(3, np.int64)
    for a in range(3):
        lo[a] = C[:, a].min() - 1e-6
        dims[a] = int((C[:, a].max() - lo[a]) / h) + 1
    cid = np.empty(n, np.int64)
    for s in range(n):
        i = int((C[s, 0] - lo[0]) / h); j = int((C[s, 1] - lo[1]) / h); k = int((C[s, 2] - lo[2]) / h)
        cid[s] = (i * dims[1] + j) * dims[2] + k
    order = np.argsort(cid)
    sc = cid[order]
    ncell = dims[0] * dims[1] * dims[2]
    start = np.searchsorted(sc, np.arange(ncell + 1))
    rmax = R.max()
    m = 0; tot = 0
    for s in range(n):
        x = C[s, 0]; y = C[s, 1]; z = C[s, 2]; r = R[s]
        reach = r + rmax - tol
        if reach <= 0:
            continue
        i0 = max(int((x - reach - lo[0]) / h), 0); i1 = min(int((x + reach - lo[0]) / h), dims[0] - 1)
        j0 = max(int((y - reach - lo[1]) / h), 0); j1 = min(int((y + reach - lo[1]) / h), dims[1] - 1)
        k0 = max(int((z - reach - lo[2]) / h), 0); k1 = min(int((z + reach - lo[2]) / h), dims[2] - 1)
        cs = cat[I[s]]
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    c = (i * dims[1] + j) * dims[2] + k
                    for e in range(start[c], start[c + 1]):
                        q = order[e]
                        if q <= s or I[q] == I[s]:
                            continue
                        cq = cat[I[q]]
                        if cs == 2 and cq == 2:
                            continue
                        rr = r + R[q] - tol
                        if rr <= 0:
                            continue
                        dx = x - C[q, 0]; dy = y - C[q, 1]; dz = z - C[q, 2]
                        d2 = dx * dx + dy * dy + dz * dz
                        if d2 < rr * rr:
                            tot += 1
                            if m < cap:
                                out_a[m] = I[s]; out_b[m] = I[q]; out_d[m] = rr - np.sqrt(d2) + tol
                                m += 1
    return out_a[:m], out_b[:m], out_d[:m], tot


def _tile_clash(t):
    g = _G
    f = os.path.join(g["tiles_dir"], f"{t}.npz")
    if not os.path.exists(f):
        return dict(tile=int(t), n_packed=0, n_prepass=0, n_mrna=0, pp=0, pq=0, qq=0, pm=0, max_depth=0.0, worst="")
    z = np.load(f)
    ty, pos, qu = z["types"], z["pos"].astype(np.float64), z["quat"]
    n0 = len(ty)
    C1, R1, I1 = _spheres(ty.astype(np.int64), pos, qu, g["st_off"], g["stc"], g["str"], 0)
    parts_C, parts_R, parts_I = [C1], [R1], [I1]
    npre = 0; nmr = 0
    if g["poff"] is not None:
        p0, p1 = int(g["poff"][t]), int(g["poff"][t + 1])
        if p1 > p0:
            pt = np.asarray(g["ptypes"][p0:p1]).astype(np.int64)
            C2, R2, I2 = _spheres(pt, np.asarray(g["ppos"][p0:p1], np.float64), np.asarray(g["pquat"][p0:p1]),
                                  g["pst_off"], g["pstc"], g["pstr"], n0)
            parts_C.append(C2); parts_R.append(R2); parts_I.append(I2); npre = p1 - p0
    if g["moff"] is not None and g["with_mrna"]:
        m0, m1 = int(g["moff"][t]), int(g["moff"][t + 1])
        if m1 > m0:
            mp = np.asarray(g["mpts"][m0:m1], np.float64); nmr = m1 - m0
            parts_C.append(mp); parts_R.append(np.full(nmr, g["mrad"])); parts_I.append(np.arange(n0 + npre, n0 + npre + nmr))
    C = np.concatenate(parts_C); R = np.concatenate(parts_R); I = np.concatenate(parts_I)
    cat = np.concatenate([np.zeros(n0, np.int8), np.ones(npre, np.int8), np.full(nmr, 2, np.int8)])
    a, b, d, tot = _clash(C, R, I, cat, g["h"], g["tol"], 2_000_000)
    res = dict(tile=int(t), n_packed=int(n0), n_prepass=int(npre), n_mrna=int(nmr), sphere_pairs=int(tot), truncated=bool(tot > len(a)))
    if len(a):
        key = np.minimum(a, b) * (len(cat) + 1) + np.maximum(a, b)
        u, first = np.unique(key, return_index=True)
        dep = np.zeros(len(u)); np.maximum.at(dep, np.searchsorted(u, key), d)
        ia, ib = a[first], b[first]
        ca, cb = cat[ia], cat[ib]
        lo_, hi_ = np.minimum(ca, cb), np.maximum(ca, cb)
        res.update(pp=int(((lo_ == 0) & (hi_ == 0)).sum()), pq=int(((lo_ == 0) & (hi_ == 1)).sum()),
                   qq=int(((lo_ == 1) & (hi_ == 1)).sum()), pm=int(((lo_ == 0) & (hi_ == 2)).sum()),
                   qm=int(((lo_ == 1) & (hi_ == 2)).sum()))
        packed = lo_ == 0
        if packed.any():
            k = np.argmax(np.where(packed, dep, -1))
            nm = lambda i: (g["libs"][ty[i]] if i < n0 else (g["plibs"][int(g["ptypes"][int(g["poff"][t]) + i - n0])] if i < n0 + npre else "mRNA"))
            res.update(max_depth=float(dep[k]), worst=f"{nm(ia[k])}|{nm(ib[k])}")
        else:
            res.update(max_depth=0.0, worst="")
    else:
        res.update(pp=0, pq=0, qq=0, pm=0, qm=0, max_depth=0.0, worst="")
    return res


def _pack_trees(inp):
    """collision shapes of the packer when they are not the library sphere trees (Inputs key pack_trees: npz with ids,
    off, C, R; yeast wc3d.pack_trees: 15 residues per sphere). C7(a) is defined on the shapes the packer avoided
    overlaps with; the library trees of the yeast library are 1.9-3.8x fatter (wc3d.py) and overlap by construction."""
    p = inp.get("pack_trees")
    if not p or not os.path.exists(p):
        return None
    z = np.load(p, allow_pickle=False)
    off = z["off"]; C = z["C"]; R = z["R"]
    return {str(i): (C[off[k]:off[k + 1]], R[off[k]:off[k + 1]]) for k, i in enumerate(z["ids"])}


def _override(trees, names, st_off, stc, strr, br):
    off = [0]; Cs = []; Rs = []; brn = np.array(br, np.float64)
    for j, n in enumerate(names):
        if n in trees:
            c, r = trees[n]
            brn[j] = float((np.linalg.norm(c, axis=1) + r).max()) if len(r) else 0.0
        else:
            c, r = stc[st_off[j]:st_off[j + 1]], strr[st_off[j]:st_off[j + 1]]
        Cs.append(np.asarray(c, np.float64)); Rs.append(np.asarray(r, np.float64)); off.append(off[-1] + len(r))
    return (np.array(off, np.int64), np.concatenate(Cs) if Cs else np.zeros((0, 3)), np.concatenate(Rs) if Rs else np.zeros(0), brn)


def _grad(geo, P, m, h=0.5):
    g = np.zeros((len(P), 3))
    for a in range(3):
        d = np.zeros(3); d[a] = h
        g[:, a] = (geo.mems(P + d)[:, m].astype(np.float64) - geo.mems(P - d)[:, m].astype(np.float64)) / (2 * h)
    return g


def orientation(geo, pos, qu, m):
    """membrane instances of membrane m: (sign error bool, tilt angle in degrees). Sign error: the library +z axis
    (third column of R(q), w-x-y-z) points against the SDF gradient of the membrane (+ side). Tilt (A8 phase 2, C9):
    angle between +z and the unit SDF gradient at the centre, arccos(ez . n), 0 = upright, 180 = flipped (the yeast
    stage-4 angle uses the normal at the projected membrane point; centres sit on the membrane, < 0.1 nm)."""
    gr = _grad(geo, pos, m)
    a_, b_, c_, d_ = qu[:, 0], qu[:, 1], qu[:, 2], qu[:, 3]
    ez = np.stack([2 * (b_ * d_ + a_ * c_), 2 * (c_ * d_ - a_ * b_), a_ * a_ - b_ * b_ - c_ * c_ + d_ * d_], 1)
    dot = (ez * gr).sum(1)
    n = np.linalg.norm(gr, axis=1) * np.linalg.norm(ez, axis=1)
    ang = np.degrees(np.arccos(np.clip(dot / np.where(n > 0, n, 1.0), -1, 1)))
    return dot < 0, ang


def tilt_summary(tilts):
    a = np.concatenate(tilts) if tilts else np.zeros(0)
    if not len(a):
        return dict(tilt_median_deg=None, tilt_frac_lt_1deg=None, tilt_n=0)
    return dict(tilt_median_deg=float(np.median(a)), tilt_frac_lt_1deg=float((a < 1).mean()), tilt_n=int(len(a)),
                tilt_p99_deg=float(np.percentile(a, 99)))


def orientation_only(cfg, bdir, inp, tiles, out_json=None):
    """signs and tilt of the membrane instances of the given tiles (the placement sample), without the clash scan:
    for builds whose placement.json predates the tilt metric (A8 phase 2). Same code path as run()."""
    regs = R_.load_regions(cfg, inp)
    geo = open_geometry(cfg, inp, bdir)
    mcode = {n: i for i, n in enumerate(geo.mem_names)}
    nreg = len(regs)
    n_by = np.zeros(nreg, np.int64); bad = np.zeros(nreg, np.int64); tilts = []
    for t in tiles:
        f = os.path.join(inp["tiles_dir"], f"{int(t)}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f)
        pos = z["pos"].astype(np.float64); qu = z["quat"].astype(np.float64)
        reg = np.repeat(z["job_region"].astype(np.int64), z["job_placed"])
        for r in np.unique(reg):
            if regs[r]["kind"] != "membrane":
                continue
            k = np.where(reg == r)[0]
            ob, ang = orientation(geo, pos[k], qu[k], mcode[regs[r]["membrane"]])
            n_by[r] += len(k); bad[r] += int(ob.sum()); tilts.append(ang)
    mem_i = R_.by_kind(regs, "membrane")
    nm = int(n_by[mem_i].sum())
    out = dict(n_tiles_sampled=int(len(tiles)), n_membrane_instances=nm,
               orientation_error_rate_per_1e6=(float(bad.sum()) / nm * 1e6) if nm else None,
               orientation_by_region={regs[r]["name"]: dict(n=int(n_by[r]), bad=int(bad[r])) for r in mem_i if n_by[r]},
               **tilt_summary(tilts), geometry=geo.source)
    if out_json:
        common.write_json(out_json, out, True)
    return out


def _pore_ok(P, pores, lim, zwin):
    ok = np.zeros(len(P), bool)
    for k in range(len(pores)):
        v = P - pores[k, :3]; z = v @ pores[k, 3:6]
        r = v - z[:, None] * pores[k, 3:6]
        ok |= (np.abs(z) <= zwin) & ((r * r).sum(1) < lim * lim)
    return ok


def run(cfg, bdir, inp: common.Inputs, out_json, out_tsv, overwrite=False, nproc=8, n_tiles=None, tiles=None):
    pc = (cfg.get("validate", {}) or {}).get("placement", {}) or {}
    t0 = time.time()
    n_all = inp.n_tiles()
    if tiles is None:
        n = int(n_tiles or pc.get("n_tiles", 300))
        tiles = common.sample_tiles(cfg, n_all, n, seed=pc.get("seed"))
    tiles = np.asarray(tiles, np.int64)
    regs = R_.load_regions(cfg, inp)
    L = json.load(open(inp["libs"]))["libs"]
    PL = json.load(open(inp["preplaced_libs"])) if inp.has("preplaced_libs") else []
    lc = common.lib_cache(inp, list(L) + PL, common.cache_path(bdir, "lib_all.npz"), nproc=nproc)
    st_off, stc, strr, br, memflag = lc.sub(L)
    trees = _pack_trees(inp)
    if trees:                      # the packer's own collision shapes (yeast: pack_trees_15), not the library trees
        st_off, stc, strr, br = _override(trees, L, st_off, stc, strr, br)
    has_pre = inp.has("preplaced_off") and len(PL) > 0
    _G.clear()
    _G.update(tiles_dir=inp["tiles_dir"], st_off=st_off, stc=stc, str=strr, libs=L, plibs=PL,
              h=float(pc.get("grid_nm", 6.0)), tol=float(pc.get("clash_nm", 0.5)), with_mrna=bool(pc.get("mrna", True)),
              poff=np.load(inp["preplaced_off"]) if has_pre else None,
              ptypes=np.load(inp["preplaced_types"], mmap_mode="r") if has_pre else None,
              ppos=np.load(inp["preplaced_pos"], mmap_mode="r") if has_pre else None,
              pquat=np.load(inp["preplaced_quat"], mmap_mode="r") if has_pre else None,
              moff=np.load(inp["mrna_off"]) if inp.has("mrna_off") else None,
              mpts=np.load(inp["mrna_pts"], mmap_mode="r") if inp.has("mrna_pts") else None,
              mrad=float(json.load(open(inp["mrna_meta"]))["bead_radius"]) if inp.has("mrna_meta") else 1.0)
    if has_pre:
        p_off, p_c, p_r, p_br, _ = lc.sub(PL)
        if trees:
            p_off, p_c, p_r, _ = _override(trees, PL, p_off, p_c, p_r, p_br)
        _G.update(pst_off=p_off, pstc=p_c, pstr=p_r)
    import multiprocessing as mp
    with ProcessPoolExecutor(nproc, mp_context=mp.get_context("fork")) as ex:
        clash = list(ex.map(_tile_clash, tiles.tolist(), chunksize=2))
    t_clash = time.time() - t0
    # containment and orientation (geometry, main process)
    geo = open_geometry(cfg, inp, bdir)
    lcode = label_code_map(geo)
    mcode = {n: i for i, n in enumerate(geo.mem_names)}
    mem_tol = float(pc.get("membrane_tol_nm", 0.5)); sh_tol = float(pc.get("shell_tol_nm", 0.5))
    zwin = float(pc.get("pore_axial_nm", 150.0))
    nreg = len(regs)
    bad_by_reg = np.zeros(nreg, np.int64); n_by_reg = np.zeros(nreg, np.int64); orient_bad = np.zeros(nreg, np.int64)
    unknown = set()
    per_tile = []
    tilts = []
    for c in clash:
        t = c["tile"]
        f = os.path.join(inp["tiles_dir"], f"{t}.npz")
        nb_t = 0; no_t = 0
        if os.path.exists(f):
            z = np.load(f)
            ty = z["types"]; pos = z["pos"].astype(np.float64); qu = z["quat"].astype(np.float64)
            reg = np.repeat(z["job_region"].astype(np.int64), z["job_placed"])
            if len(ty):
                lab = geo.labels(pos)
                need_m = np.array([regs[r]["kind"] in ("membrane", "shell") for r in range(nreg)])[reg]
                mem = geo.mems(pos[need_m]) if need_m.any() else None
                midx = -np.ones(len(ty), np.int64); midx[need_m] = np.arange(need_m.sum())
                for r in np.unique(reg):
                    R = regs[r]; k = np.where(reg == r)[0]; n_by_reg[r] += len(k)
                    if R["kind"] in ("volume", "virtual"):
                        allowed = [lcode[x] for x in R["labels"] if x in lcode]
                        if not allowed: unknown.add(R["name"]); continue
                        bad = ~np.isin(lab[k], allowed)
                    elif R["kind"] == "membrane":
                        m = mcode[R["membrane"]]; d = mem[midx[k], m]
                        bad = np.abs(d) > mem_tol
                        ob, ang = orientation(geo, pos[k], qu[k], m)
                        orient_bad[r] += int(ob.sum()); no_t += int(ob.sum()); tilts.append(ang)
                    elif R["kind"] == "shell":
                        m = mcode[R["membrane"]]; side = float(R.get("side", 1.0))
                        d1 = R_.resolve_d1(R, geo.mem_half[m])
                        u = mem[midx[k], m].astype(np.float64)
                        u = np.abs(u) if side == 0.0 else side * u          # side 0: both sides (|SDF|, yeast bands)
                        hi = d1 + (br[ty[k]] if R["add_radius"] else 0.0)
                        allowed = [lcode[x] for x in R["labels"] if x in lcode]
                        bad = (u < float(R.get("d0", 0.0)) - sh_tol) | (u > hi + sh_tol)
                        if allowed:
                            bad |= ~np.isin(lab[k], allowed)
                        if R["pore_lat"] > 0 and len(geo.pores):
                            bad |= ~_pore_ok(pos[k], geo.pores, R["pore_lat"] + sh_tol, zwin)
                    else:
                        unknown.add(R["name"]); continue
                    bad_by_reg[r] += int(bad.sum()); nb_t += int(bad.sum())
        per_tile.append(dict(c, containment_bad=nb_t, orientation_bad=no_t))
    T = pd.DataFrame(per_tile)
    common.write_df(out_tsv, T, overwrite, index=False)
    npk = int(T.n_packed.sum())
    mem_i = R_.by_kind(regs, "membrane")
    n_mem = int(n_by_reg[mem_i].sum()) if mem_i else 0
    per_m = lambda x, n: (float(x) / n * 1e6) if n else None
    # whole-cell per-species placed share (job tables)
    S = common.tile_sums(inp, common.cache_path(bdir, "tile_sums.npz"), nproc=nproc)
    min_req = float(pc.get("min_req", 100))
    share = {}
    for r in range(nreg):
        rq = S["req"][:, r]; pl = S["pla"][:, r]
        if rq.sum() <= 0: continue
        k = rq >= min_req; f = np.where(rq > 0, pl / np.maximum(rq, 1), np.nan)
        share[regs[r]["name"]] = dict(requested=int(rq.sum()), placed=int(pl.sum()), frac=round(float(pl.sum() / rq.sum()), 4),
                                      n_species=int(k.sum()), n_below_0_5=int((f[k] < 0.5).sum()), n_below_0_8=int((f[k] < 0.8).sum()))
    out = dict(
        n_tiles_sampled=int(len(tiles)), tiles_md5=__import__("hashlib").md5(tiles.tobytes()).hexdigest(),
        n_packed=npk, n_prepass=int(T.n_prepass.sum()), n_mrna_beads=int(T.n_mrna.sum()), n_membrane_instances=n_mem,
        clash_rate_per_1e6=per_m(T.pp.sum() + T.pq.sum(), npk),
        clash_pairs=dict(packed_packed=int(T.pp.sum()), packed_prepass=int(T.pq.sum()), prepass_prepass=int(T.qq.sum()),
                         packed_mrna=int(T.pm.sum()), prepass_mrna=int(T.get("qm", pd.Series([0])).sum())),
        clash_rate_mrna_per_1e6=per_m(T.pm.sum(), npk),
        max_clash_depth_nm=float(T.max_depth.max()) if len(T) else 0.0,
        worst_pair=str(T.loc[T.max_depth.idxmax(), "worst"]) if len(T) and T.max_depth.max() > 0 else "",
        containment_rate_per_1e6=per_m(bad_by_reg.sum(), npk),
        containment_by_region={regs[r]["name"]: dict(n=int(n_by_reg[r]), bad=int(bad_by_reg[r])) for r in range(nreg) if n_by_reg[r]},
        orientation_error_rate_per_1e6=per_m(orient_bad.sum(), n_mem),
        orientation_by_region={regs[r]["name"]: dict(n=int(n_by_reg[r]), bad=int(orient_bad[r])) for r in mem_i if n_by_reg[r]},
        **tilt_summary(tilts),
        regions_not_checked=sorted(unknown), species_share=share,
        settings=dict(collision_shapes=inp.get("pack_trees") or "library sphere trees (st_c, st_r)", clash_nm=_G["tol"], grid_nm=_G["h"], membrane_tol_nm=mem_tol, shell_tol_nm=sh_tol, pore_axial_nm=zwin,
                      mrna=_G["with_mrna"], geometry=geo.source, tiles=tiles.tolist()),
        seconds=dict(clash=round(t_clash, 1), total=round(time.time() - t0, 1)))
    common.write_json(out_json, out, overwrite)
    return out
