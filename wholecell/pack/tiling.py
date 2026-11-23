# origin: tile measures from wholecell/pack/stage.py (_measure / step_tiles, jurkat jwc.py); tile list, voxel region masks,
# label volumes and per-tile membrane meshes from yeast yeast_cell_dev_bundle/pipeline/tiles3d.py @ e3b7448
# (tile_list, _voxel_regions, region_masks, tile_geometry) and the prim / pore box selection of yeast
# pipeline/cell3d.py:Cell.select @ e3b7448. Changed: region ids, label masks, bands and the wall depth come from config
# `pack.regions` (no REG_MASK / BANDS tables); the evaluator is the core cellgeom_eval_v1; tile size, voxel size,
# selection margins and the mesh membranes are `pack:` config values (A2, 2026-11-17).
"""Tiling (pack step `tiles`): tile list, per-tile region measures, label volumes and membrane meshes.

Config (`pack:`; defaults = the Jurkat v7 tiling):
  tile_nm               tile edge (Jurkat 500, yeast 200)
  meas_h_nm             membrane-area lattice spacing (marching cubes) and, unless meas_vox_nm is set, the voxel size
                        of the volume measures (Jurkat 8)
  meas_vox_nm           voxel size of the volume measures (yeast 4)
  tile_keep_half_diag_f select centre_depth: keep tiles whose centre lies within this many tile edges outside the cell
  tiling:
    select              centre_depth (v7: centre near the cell, tiles without content dropped) | enclosing_sphere
                        (yeast: tiles meeting the enclosing sphere minus tiles whose 27 sample points all lie more than
                        min(cell-field clamp - 1, sample_reach_nm) outside the cell field; empty tiles kept)
    sample_reach_nm     enclosing_sphere: distance within which every tile point has a sample (88: 27 samples 100 nm
                        apart reach every point within 86.6 nm)
    box_margin_nm       evaluate each tile with only the prims whose bounding sphere meets the tile grown by this
                        margin (null: all prims and pores; v7)
    pore_reach_nm       ... and the pores within this distance of that grown box
    labvol              true: also write tile_labvol.parquet (volume per engine label per tile)
    meshes              true: also write tile_geom/<tile name>.npz (lo, measure, labvol, verts, faces, face_mem, area)
    mesh_membranes      membranes meshed into tile_geom (default: the membrane regions)
Region measures (columns of tiles.parquet, Regions.measures): volume-like regions = voxel count x voxel volume of the
voxel centres whose label is one of the region's labels (+ shell band / sdf_lt restrictions, + virtual providers);
membranes = marching-cubes area of the membrane SDF on the meas_h_nm lattice (open face clipped); preplaced ids = 0.
"""
from __future__ import annotations

import os

import numpy as np

_S = {}          # state shared with forked workers (set by stage.step_tiles)


def tile_name(ijk) -> str:
    """t_+00_-03_+11: tile with lower corner ijk * tile_nm (model frame)"""
    return f"t_{int(ijk[0]):+03d}_{int(ijk[1]):+03d}_{int(ijk[2]):+03d}"


def grid(cell, TILE):
    """the tile lattice covering the enclosing sphere: lower corner lo0 (nm) and the number of tiles per axis"""
    R = cell.G[0]
    return -np.full(3, np.ceil(R / TILE) * TILE), int(2 * np.ceil(R / TILE))


def select_centre_depth(cell, c3, TILE, keep_f, outer):
    """v7: grid indices of the tiles whose centre lies within keep_f tile edges outside the outer membrane"""
    lo0, nt = grid(cell, TILE)
    I, J, K = np.meshgrid(*[np.arange(nt)] * 3, indexing="ij"); T = np.stack([I.ravel(), J.ravel(), K.ravel()], 1)
    C = T * TILE + lo0 + TILE / 2
    dep = cell.mems(C.astype(np.float64))[:, outer]
    return T[dep > -TILE * keep_f], lo0, nt


def select_enclosing_sphere(cell, c3, TILE, reach):
    """origin: yeast pipeline/tiles3d.py:tile_list @ e3b7448. Grid indices (same lattice and order as select_centre_depth)
    of the tiles that meet the enclosing sphere, minus tiles certainly outside the cell field"""
    lo0, nt = grid(cell, TILE)
    R = cell.G[0]
    n = nt // 2
    out = []
    has_field = c3.n_fields(cell.G) > c3.F_CELL
    if has_field:
        clamp = cell.G[c3.G_N + 1 + 12 * c3.F_CELL + 9]       # the cell field is exact to `clamp` outside the cell
        S = np.array([[i, j, k] for i in (0, .5, 1) for j in (0, .5, 1) for k in (0, .5, 1)])
    for ix in range(-n, n):
        for iy in range(-n, n):
            for iz in range(-n, n):
                lo = np.array([ix, iy, iz]) * TILE; hi = lo + TILE
                q = np.clip(0.0, lo, hi)
                if np.linalg.norm(q) < R:
                    if has_field:
                        g = lo + TILE * S
                        if max(c3.outer_depth(p[0], p[1], p[2], cell.G) for p in g) <= -min(clamp - 1.0, reach):
                            continue
                    out.append((ix + n, iy + n, iz + n))
    return np.array(out, np.int64).reshape(-1, 3), lo0, nt


def box_subset(cell, lo, hi, margin, pore_reach):
    """origin: yeast pipeline/cell3d.py:Cell.select @ e3b7448. Prims whose bounding sphere meets the box [lo, hi] grown
    by margin, and the pores whose centre lies within pore_reach of that grown box -> (pidx, pores)"""
    lo = np.asarray(lo) - margin; hi = np.asarray(hi) + margin
    bs = cell.prim_bsphere
    bc, br = bs[:, :3], bs[:, 3]
    q = np.clip(bc, lo, hi)
    pidx = np.where(((q - bc) ** 2).sum(1) <= br ** 2)[0].astype(np.int64)
    pc = cell.pores[:, :3]
    q = np.clip(pc, lo, hi)
    pk = np.where(((q - pc) ** 2).sum(1) <= pore_reach ** 2)[0]
    return pidx, np.ascontiguousarray(cell.pores[pk])


class Evaluator:
    """labels / membrane SDFs at points, with all prims (v7) or the prims and pores near one tile (box_margin_nm)"""

    def __init__(self, cell, c3, lo, TILE, margin, pore_reach):
        self.cell, self.c3 = cell, c3
        if margin is None:
            self.args = None
        else:
            pidx, pores = box_subset(cell, lo, lo + TILE, float(margin), float(pore_reach))
            self.args = (cell.G, pores, pidx, cell.PT, cell.PQ, cell.PC)

    def labels(self, P):
        return self.cell.labels(P) if self.args is None else self.c3.label_points(np.asarray(P, np.float64), *self.args)

    def mems(self, P):          # float32 (marching-cubes lattice)
        return self.cell.mems(P) if self.args is None else self.c3.mem_points(np.asarray(P, np.float64), *self.args)

    def mems64(self, P):        # float64 (band / depth tests)
        return self.cell.mems64(P) if self.args is None else self.c3.mem_points64(np.asarray(P, np.float64), *self.args)


def _restriction(r, reg):
    """(membrane code, side, d0, d1) of a volume-measured shell, (membrane code, 'lt', None, nm) of sdf_lt, or None"""
    return reg.restriction(r["name"])


def region_voxels(reg, ev, Pm, lab, vset=None):
    """{volume-like region name: bool mask over the voxel centres Pm (labels lab)}: the region's labels, restricted by
    its shell band / sdf_lt (float64 membrane SDFs of the evaluator ev); virtual regions from vset. The same rule as the
    tile measures (measure_tiles), used for the packing seeds (A7)."""
    out = {}
    restr = {n: reg.restriction(n) for n, _ in reg.VOL_REG if reg.spec[n]["kind"] != "virtual"}
    need = sorted({x[0] for x in restr.values() if x is not None})
    S64 = None
    if need:
        cand = np.zeros(len(lab), bool)
        for n, labs in reg.VOL_REG:
            if restr.get(n) is not None:
                cand |= np.isin(lab, labs)
        S64 = np.zeros((len(lab), reg_nmems(ev)))
        if cand.any():
            S64[cand] = ev.mems64(Pm[cand])
    for n, labs in reg.VOL_REG:
        if reg.spec[n]["kind"] == "virtual":
            continue
        m = np.isin(lab, labs)
        x = restr[n]
        if x is not None:
            s = S64[:, x[0]]
            if x[1] == "lt":
                m &= s < x[3]
            else:
                u = np.abs(s) if x[1] == 0.0 else x[1] * s
                k = u < x[3]
                if x[2] > 0.0:
                    k &= u >= x[2]
                m &= k
        out[n] = m
    if vset:
        mk = vset.masks(Pm, lab)
        for n in reg.VIRTUAL:
            if n in mk:
                out[n] = mk[n] if mk[n] is not None else np.zeros(len(lab), bool)
    return out


def reg_nmems(ev):
    return ev.cell.n_mems


def measure_tiles(idx):
    """worker: measures of the tiles idx (rows of _S['T']) -> [[i] + column values]; writes labvol rows / meshes"""
    from skimage.measure import marching_cubes
    from .stage import face_keep
    p, reg, cell, c3 = _S["p"], _S["reg"], _S["cell"], _S["c3"]
    tc = _S["tcfg"]
    TILE, MEAS_H, VOX = p.TILE, p.MEAS_H, p.MEAS_VOX
    cols = reg.measures()
    vcodes = {ri: tuple(reg.label_codes[l] for l in r["labels"]) for ri, _, kind, r in cols if kind == "volume"}
    restr = {ri: _restriction(r, reg) for ri, _, kind, r in cols if kind == "volume"}
    need64 = sorted({x[0] for x in restr.values() if x is not None})
    area_mems = [reg.mem_codes[r["membrane"]] for _, _, kind, r in cols if kind == "area"]
    mesh_mems = [reg.mem_codes[m] for m in tc.get("mesh_membranes", [])] if tc.get("meshes") else []
    if tc.get("meshes") and not tc.get("mesh_membranes"):
        mesh_mems = list(area_mems)
    lat_mems = sorted(set(area_mems) | set(mesh_mems))
    nlab = _S["nlab"]
    out = []
    for i in idx:
        lo = _S["T"][i] * TILE + _S["lo0"]
        ev = Evaluator(cell, c3, lo, TILE, tc.get("box_margin_nm"), tc.get("pore_reach_nm"))
        n = int(TILE / VOX); g = (np.arange(n) + 0.5) * VOX
        X, Y, Z = np.meshgrid(lo[0] + g, lo[1] + g, lo[2] + g, indexing="ij")
        Pm = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1); lab = ev.labels(Pm)
        # membrane SDFs (float64) at the voxels that a restricted region may use
        S64 = None
        if need64:
            cand = np.zeros(len(lab), bool)
            for ri, x in restr.items():
                if x is not None:
                    cand |= np.isin(lab, vcodes[ri])
            S64 = np.zeros((len(lab), cell.n_mems))
            if cand.any():
                S64[cand] = ev.mems64(Pm[cand])
        vals = {}
        for ri, col, kind, r in cols:
            if kind != "volume":
                continue
            x = restr[ri]
            if x is None:
                vals[col] = float(np.isin(lab, vcodes[ri]).sum() * VOX ** 3)
                continue
            s = S64[:, x[0]]
            if x[1] == "lt":
                keep = s < x[3]
            else:
                u = np.abs(s) if x[1] == 0.0 else x[1] * s
                keep = u < x[3]
                if x[2] > 0.0:
                    keep &= u >= x[2]
            vals[col] = float((np.isin(lab, vcodes[ri]) & keep).sum() * VOX ** 3)
        vs = _S["vset"]
        if vs:
            mk = vs.masks(Pm, lab)
            for ri, col, kind, r in cols:
                if r["kind"] == "virtual" and r["name"] in mk:
                    vals[col] = float(mk[r["name"]].sum() * VOX ** 3) if mk[r["name"]] is not None else 0.0
        # membranes: marching cubes on the meas_h_nm lattice
        m = int(TILE / MEAS_H) + 1; gm = np.arange(m) * MEAS_H
        X, Y, Z = np.meshgrid(lo[0] + gm, lo[1] + gm, lo[2] + gm, indexing="ij")
        S = ev.mems(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(m, m, m, -1)
        area = np.zeros(cell.n_mems)
        verts, faces, fmem = [], [], []
        nv = 0
        for mi in lat_mems:
            s = S[..., mi]
            if not (s.min() < 0 < s.max()):
                continue
            v, f, _, _ = marching_cubes(s, 0.0, spacing=(MEAS_H,) * 3)
            tt = face_keep((v + lo)[f].reshape(-1, 9), cell.G, c3)
            area[mi] = float(0.5 * np.linalg.norm(np.cross(tt[:, 3:6] - tt[:, 0:3], tt[:, 6:9] - tt[:, 0:3]), axis=1).sum())
            if mi in mesh_mems and len(f):
                verts.append((v + lo).astype(np.float32)); faces.append((f + nv).astype(np.int32))
                fmem.append(np.full(len(f), mi, np.int8)); nv += len(v)
        for ri, col, kind, r in cols:
            if kind == "area":
                vals[col] = float(area[reg.mem_codes[r["membrane"]]])
            elif kind == "zero":
                vals[col] = 0.0
        row = [vals[col] for _, col, _, _ in cols]
        labvol = np.bincount(lab.ravel(), minlength=nlab) * VOX ** 3 if (tc.get("labvol") or tc.get("meshes")) else None
        if tc.get("meshes"):
            from .stage import _atomic
            nm = tile_name(_S["T"][i] + np.round(_S["lo0"] / TILE).astype(np.int64))
            d = dict(lo=lo, measure=np.array(row, np.float64), labvol=labvol,
                     verts=np.concatenate(verts) if verts else np.zeros((0, 3), np.float32),
                     faces=np.concatenate(faces) if faces else np.zeros((0, 3), np.int32),
                     face_mem=np.concatenate(fmem) if fmem else np.zeros(0, np.int8), area=area)
            _atomic(os.path.join(p.out, "tile_geom", nm + ".npz"), lambda tmp: np.savez_compressed(tmp, **d))
        out.append([i] + row + ([labvol] if labvol is not None else []))
    return out
