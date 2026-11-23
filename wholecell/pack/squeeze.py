# origin: yeast squeeze.py @ e3b7448 (kernels _grid, _rotq, _wall, relax, accept; driver from wc3d.py:squeeze_tile)
# Changed: 'import cell3d as c3' -> core engine; HALF / min_mem_dist -> per-package values (c3.half(G), membrane
# clearance |d_m| - mem_half(G, m) as in pack3d.check_volume); accept() takes the core shell arrays instead of yeast
# bands and only volume-region molecules are movable (shells are fixed obstacles here); the driver is a function of
# explicit arrays (squeeze_tile) with a numpy Generator from a named seed stream; constants are arguments.
"""Squeeze: overlap relaxation after random sequential addition (pack.squeeze, unification plan P5).

RSA stops near its jamming limit (~35-40 % sphere-tree volume fraction). After a tile has been packed, copies that did
not fit are inserted anyway at random seed voxels of their region, and the movable molecules of the tile are relaxed as
rigid bodies: overlapping spheres of different molecules push each other apart (translation and a small rotation),
spheres that cross a membrane or leave their region are pushed back along the membrane distance gradient, molecule
centres stay in the tile. Membrane proteins, shells, pre-placed structures and neighbour tiles are fixed obstacles.
The result is then made exact again (accept): molecules are kept one by one (fixed ones first, then movable ones in
order) only if they overlap nothing kept and lie in their region, as in RSA. A round that does not raise the count is
undone.
"""
import numpy as np
import numba as nb

from wholecell.geometry import engine as c3
from . import pack3d as p3

H_MIN = 1.0             # nm: smallest neighbour-grid cell (the cell is max movable radius + max radius: exact search)
STEP_MAX = 0.6          # nm per iteration
ROT_MAX = 0.06          # rad per iteration
WALL_EVERY = 4          # wall (membrane / region) forces every this many iterations


# origin: yeast squeeze.py:_grid @ e3b7448
@nb.njit(cache=True)
def _grid(lo, dims, h, P, R):
    head = -np.ones(dims[0] * dims[1] * dims[2], np.int64); nxt = np.empty(P.shape[0], np.int64)
    for s in range(P.shape[0]):
        i = min(max(int((P[s, 0] - lo[0]) / h), 0), dims[0] - 1)
        j = min(max(int((P[s, 1] - lo[1]) / h), 0), dims[1] - 1)
        k = min(max(int((P[s, 2] - lo[2]) / h), 0), dims[2] - 1)
        c = (i * dims[1] + j) * dims[2] + k
        nxt[s] = head[c]; head[c] = s
    return head, nxt


# origin: yeast squeeze.py:_rotq @ e3b7448
@nb.njit(cache=True, fastmath=True)
def _rotq(q, wx, wy, wz):
    """q <- exp(w/2) * q for a small rotation vector w"""
    a = np.sqrt(wx * wx + wy * wy + wz * wz)
    if a < 1e-12:
        return
    s = np.sin(0.5 * a) / a; c = np.cos(0.5 * a)
    b0, b1, b2, b3 = c, wx * s, wy * s, wz * s
    a0, a1, a2, a3 = q[0], q[1], q[2], q[3]
    q[0] = b0 * a0 - b1 * a1 - b2 * a2 - b3 * a3
    q[1] = b0 * a1 + b1 * a0 + b2 * a3 - b3 * a2
    q[2] = b0 * a2 - b1 * a3 + b2 * a0 + b3 * a1
    q[3] = b0 * a3 + b1 * a2 - b2 * a1 + b3 * a0
    n = np.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2 + q[3] ** 2)
    for t in range(4):
        q[t] /= n


# origin: yeast squeeze.py:_wall @ e3b7448
# Changed: clearance of membrane k = |d_k| - mem_half(G, k) (package values), as pack3d.check_volume tests it.
@nb.njit(cache=True, fastmath=True)
def _wall(x, y, z, r, mask, G, pores, pidx, PT, PQ, PC, out, g):
    """push (g) for a sphere that crosses a membrane or sits outside its allowed labels; returns the violation (nm)"""
    m = c3.mem_sdfs(x, y, z, G, pores, pidx, PT, PQ, PC, out)
    l = c3.label_from(x, y, z, G, pores, out, m)
    dmin = 1e9; mi = -1
    for k in range(c3.n_mems(G)):                 # A7: every membrane of the package (yeast: 10, as the yeast _wall)
        dk = abs(out[k]) - c3.mem_half(G, k)
        if dk < dmin:
            dmin = dk; mi = k
    ok_label = (mask >> l) & 1
    viol = r - dmin
    g[0] = 0.0; g[1] = 0.0; g[2] = 0.0
    if ok_label and viol <= 0:
        return 0.0
    e = 0.5
    for a in range(3):
        xp = x + (e if a == 0 else 0.0); yp = y + (e if a == 1 else 0.0); zp = z + (e if a == 2 else 0.0)
        c3.mem_sdfs(xp, yp, zp, G, pores, pidx, PT, PQ, PC, out); dp = out[mi]
        xm = x - (e if a == 0 else 0.0); ym = y - (e if a == 1 else 0.0); zm = z - (e if a == 2 else 0.0)
        c3.mem_sdfs(xm, ym, zm, G, pores, pidx, PT, PQ, PC, out); dm = out[mi]
        g[a] = (dp - dm) / (2 * e)
    gn = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2) + 1e-12
    c3.mem_sdfs(x, y, z, G, pores, pidx, PT, PQ, PC, out); d = out[mi]
    # on the allowed side: move away from the membrane; on the wrong side: move across it
    sgn = 1.0 if d > 0 else -1.0
    if not ok_label:
        sgn = -sgn
        viol = abs(d) + r + c3.half(G)
    for a in range(3):
        g[a] = g[a] / gn * sgn * viol
    return viol


# origin: yeast squeeze.py:relax @ e3b7448
# Changed: near-membrane test with the package clearance (c3.min_mem_dist_g, c3.half(G)); step_max, rot_max,
# wall_every are arguments (module constants above); the neighbour-grid cell is the largest possible contact distance
# (largest movable sphere radius + largest radius, at least h_min) instead of a fixed 8 nm: exact for any tree and
# 5-10x fewer pair tests with the ~1.5 nm library spheres.
@nb.njit(cache=True, fastmath=True)
def relax(lo, size, halo, G, pores, pidx, PT, PQ, PC, st_off, STC, STR, regmask,
          typ, pos, quat, reg, movable, OC, OR, n_iter, h_min, step_max, rot_max, wall_every):
    """typ/pos/quat/reg: all tile instances (movable ones are relaxed in place); OC/OR: fixed spheres. Deterministic
    (no random draws). Returns the summed overlap + wall violation (nm) of the last iteration."""
    n = typ.shape[0]
    nsm = 0; nsf = 0
    for i in range(n):
        k = st_off[typ[i] + 1] - st_off[typ[i]]
        if movable[i]:
            nsm += k
        else:
            nsf += k
    NS = nsm + nsf + OR.shape[0]
    P = np.empty((NS, 3)); R = np.empty(NS); own = np.empty(NS, np.int64)
    s0 = np.empty(n, np.int64); s1 = np.empty(n, np.int64)
    s = 0
    for pas in range(2):
        for i in range(n):
            if (pas == 0) != movable[i]:
                continue
            a0 = st_off[typ[i]]; a1 = st_off[typ[i] + 1]
            s0[i] = s; s1[i] = s + a1 - a0
            for k in range(a1 - a0):
                R[s + k] = STR[a0 + k]; own[s + k] = i
            s += a1 - a0
    for k in range(OR.shape[0]):
        P[s + k, 0] = OC[k, 0]; P[s + k, 1] = OC[k, 1]; P[s + k, 2] = OC[k, 2]; R[s + k] = OR[k]; own[s + k] = -1
    M = np.empty((3, 3)); qq = np.empty(4)
    for i in range(n):
        _place(i, typ, pos, quat, st_off, STC, s0, s1, P, M, qq)
    rm = 0.0; ra = 0.0
    for k in range(NS):
        ra = max(ra, R[k])
        if k < nsm:
            rm = max(rm, R[k])
    h = max(h_min, rm + ra)
    glo = lo - halo
    dims = np.empty(3, np.int64)
    for a in range(3):
        dims[a] = int(np.ceil((size + 2 * halo) / h)) + 1
    F = np.zeros((n, 3)); T = np.zeros((n, 3)); I = np.ones(n)
    for i in range(n):
        if movable[i]:
            acc = 1.0
            for k in range(s0[i], s1[i]):
                acc += (P[k, 0] - pos[i, 0]) ** 2 + (P[k, 1] - pos[i, 1]) ** 2 + (P[k, 2] - pos[i, 2]) ** 2
            I[i] = acc / max(s1[i] - s0[i], 1)
    out = np.zeros(c3.NMEM_MAX); g = np.zeros(3)
    near_mem = np.zeros(n, np.bool_)
    total = 0.0
    for it in range(n_iter):
        if it % wall_every == 0:          # wall forces only for molecules whose bounding sphere may reach a membrane
            for i in range(n):
                if movable[i]:
                    c3.mem_sdfs(pos[i, 0], pos[i, 1], pos[i, 2], G, pores, pidx, PT, PQ, PC, out)
                    br = 0.0
                    for k in range(s0[i], s1[i]):
                        br = max(br, np.sqrt((P[k, 0] - pos[i, 0]) ** 2 + (P[k, 1] - pos[i, 1]) ** 2 + (P[k, 2] - pos[i, 2]) ** 2) + R[k])
                    near_mem[i] = c3.min_mem_dist_g(out, -1, G) < br + c3.half(G) + 2.0
        head, nxt = _grid(glo, dims, h, P, R)
        F[:, :] = 0.0; T[:, :] = 0.0; total = 0.0
        for i in range(n):
            if not movable[i]:
                continue
            for sp in range(s0[i], s1[i]):
                x, y, z, r = P[sp, 0], P[sp, 1], P[sp, 2], R[sp]
                ci = min(max(int((x - glo[0]) / h), 0), dims[0] - 1)
                cj = min(max(int((y - glo[1]) / h), 0), dims[1] - 1)
                ck = min(max(int((z - glo[2]) / h), 0), dims[2] - 1)
                for di in range(max(ci - 1, 0), min(ci + 2, dims[0])):
                    for dj in range(max(cj - 1, 0), min(cj + 2, dims[1])):
                        for dk in range(max(ck - 1, 0), min(ck + 2, dims[2])):
                            e = head[(di * dims[1] + dj) * dims[2] + dk]
                            while e >= 0:
                                if own[e] != i:
                                    dx = x - P[e, 0]; dy = y - P[e, 1]; dz = z - P[e, 2]
                                    d2 = dx * dx + dy * dy + dz * dz; rr = r + R[e]
                                    if d2 < rr * rr:
                                        d = np.sqrt(d2) + 1e-9; ov = rr - d
                                        w = 0.5 if (own[e] >= 0 and movable[own[e]]) else 1.0
                                        fx = w * ov * dx / d; fy = w * ov * dy / d; fz = w * ov * dz / d
                                        F[i, 0] += fx; F[i, 1] += fy; F[i, 2] += fz
                                        rx = x - pos[i, 0]; ry = y - pos[i, 1]; rz = z - pos[i, 2]
                                        T[i, 0] += ry * fz - rz * fy; T[i, 1] += rz * fx - rx * fz; T[i, 2] += rx * fy - ry * fx
                                        total += ov
                                e = nxt[e]
                if it % wall_every == 0 and near_mem[i]:
                    v = _wall(x, y, z, r, regmask[reg[i]], G, pores, pidx, PT, PQ, PC, out, g)
                    if v > 0:
                        F[i, 0] += g[0]; F[i, 1] += g[1]; F[i, 2] += g[2]; total += v
        if total < 1e-3 * max(nsm, 1):
            break
        for i in range(n):
            if not movable[i]:
                continue
            fn = np.sqrt(F[i, 0] ** 2 + F[i, 1] ** 2 + F[i, 2] ** 2)
            sc = 1.0 / max(1.0, s1[i] - s0[i]) ** 0.5
            stp = min(step_max, fn * sc) / (fn + 1e-12)
            for a in range(3):
                pos[i, a] += F[i, a] * stp
                pos[i, a] = min(max(pos[i, a], lo[a] + 1e-3), lo[a] + size - 1e-3)
            tn = np.sqrt(T[i, 0] ** 2 + T[i, 1] ** 2 + T[i, 2] ** 2)
            ang = min(rot_max, 0.5 * tn / (I[i] * max(1.0, s1[i] - s0[i])))
            if tn > 1e-9:
                for t in range(4):
                    qq[t] = quat[i, t]
                _rotq(qq, T[i, 0] / tn * ang, T[i, 1] / tn * ang, T[i, 2] / tn * ang)
                for t in range(4):
                    quat[i, t] = qq[t]
            _place(i, typ, pos, quat, st_off, STC, s0, s1, P, M, qq)
    return total


# origin: yeast squeeze.py:relax.place_spheres @ e3b7448 (closure made a function: explicit arguments)
@nb.njit(cache=True, fastmath=True)
def _place(i, typ, pos, quat, st_off, STC, s0, s1, P, M, qq):
    a0 = st_off[typ[i]]
    for t in range(4):
        qq[t] = quat[i, t]
    p3.quat_mat(qq, M)
    for k in range(s1[i] - s0[i]):
        for a in range(3):
            P[s0[i] + k, a] = M[a, 0] * STC[a0 + k, 0] + M[a, 1] * STC[a0 + k, 1] + M[a, 2] * STC[a0 + k, 2] + pos[i, a]


# origin: yeast squeeze.py:accept @ e3b7448
# Changed: the band test (band_m, band_w) applies to volume-measured shells only (Regions.band_arrays; empty = none).
@nb.njit(cache=True, fastmath=True)
def accept(lo, size, halo, G, pores, pidx, PT, PQ, PC, st_off, STC, STR, BR, regmask,
           typ, pos, quat, reg, movable, OC, OR, cap_s, cap_e, band_m, band_w):
    """keep[i]: fixed instances always; movable ones in order if they overlap nothing kept and lie in their region
    (and, for a volume-measured shell, within band_w of its membrane band_m: the relaxation must not push them out)"""
    glo = lo - halo; ghi = lo + size + halo
    dims, head, nxt, sid, cen, rad, cnt = p3.grid_new(glo, ghi, p3.H_GRID, cap_s, cap_e)
    p3.g_insert(glo, p3.H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    n = typ.shape[0]; keep = np.zeros(n, np.bool_)
    M = np.empty((3, 3)); q = np.empty(4); out = np.zeros(c3.NMEM_MAX)
    for pas in range(2):
        for i in range(n):
            if (pas == 0) == bool(movable[i]):
                continue
            a0 = st_off[typ[i]]; a1 = st_off[typ[i] + 1]
            for t in range(4):
                q[t] = quat[i, t]
            p3.quat_mat(q, M); Pi = np.empty((a1 - a0, 3)); p3.transform(STC[a0:a1], M, pos[i], Pi)
            if movable[i]:
                if p3.g_overlap(glo, p3.H_GRID, dims, head, nxt, sid, cen, rad, Pi, STR[a0:a1], p3.TOL):
                    continue
                if not p3.check_volume(Pi, STR[a0:a1], pos[i, 0], pos[i, 1], pos[i, 2], BR[typ[i]], regmask[reg[i]],
                                       G, pores, pidx, PT, PQ, PC, out):
                    continue
                if band_m.shape[0] > reg[i] and band_m[reg[i]] >= 0:
                    c3.mem_sdfs(pos[i, 0], pos[i, 1], pos[i, 2], G, pores, pidx, PT, PQ, PC, out)
                    if abs(out[band_m[reg[i]]]) > band_w[reg[i]]:
                        continue
            if not p3.g_insert(glo, p3.H_GRID, dims, head, nxt, sid, cen, rad, cnt, Pi, STR[a0:a1]):
                break
            keep[i] = True
    return keep


# origin: yeast wc3d.py:squeeze_tile @ e3b7448
# Changed: explicit arguments (cell arrays, region tables, squeeze settings); `sq_jobs` marks the jobs whose missing
# copies are inserted and whose molecules move (targeted squeeze: only underfilled regions); exclusion zones re-tested
# with the core rule bitmasks.
def squeeze_tile(rng, lo, size, halo, cell, st_off, STC, STR, BR, regmask, jl, jr, jn, placed, types, pos, quat,
                 vox_off, VOXP, vox_h, OC, OR, sq_jobs, n_rounds, n_iter, insert_frac, EXC=None, EXC_TYPE=None,
                 args=None, band=None):
    """insert the copies RSA could not place (at most insert_frac x movable molecules, >= 50, per round) into the jobs
    marked in sq_jobs, relax the molecules of those jobs, keep the overlap-free ones; a round that does not raise the
    count is undone. Instances are grouped by job (as pack_tile returns them). Returns types, pos, quat, placed.
    args: evaluator arguments (G, pores, pidx, PT, PQ, PC) of the tile (default: the whole cell); band: (band_m,
    band_w) per region (Regions.band_arrays; default none)."""
    typ = types.astype(np.int64); P = pos.astype(np.float64); Q = quat.astype(np.float64)
    J = np.repeat(np.arange(len(jl)), placed)
    vol_job = np.array([bool(sq_jobs[j]) and vox_off[jr[j] + 1] > vox_off[jr[j]] for j in range(len(jl))], bool)
    if args is None:
        args = (cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC)
    bm, bw = band if band is not None else (np.zeros(0, np.int64), np.zeros(0))
    for rnd in range(n_rounds):
        have = np.bincount(J, minlength=len(jl)); deficit = np.where(vol_job, jn - have, 0).clip(min=0)
        if deficit.sum() == 0:
            break
        n_mov = int(vol_job[J].sum())
        take = deficit * min(1.0, max(50, insert_frac * n_mov) / deficit.sum())
        take = np.floor(take + rng.random(len(take))).astype(np.int64)
        nt, nj, npos = [], [], []
        for j in np.where(take > 0)[0]:
            r = jr[j]; v = rng.integers(vox_off[r], vox_off[r + 1], take[j])
            npos.append(VOXP[v, :3] + (rng.random((take[j], 3)) - 0.5) * vox_h); nt.append(np.full(take[j], jl[j])); nj.append(np.full(take[j], j))
        if not nt:
            break
        nq = rng.normal(size=(sum(map(len, nt)), 4)); nq /= np.linalg.norm(nq, axis=1, keepdims=True)
        typ2 = np.concatenate([typ, np.concatenate(nt)]); P2 = np.concatenate([P, np.concatenate(npos)])
        Q2 = np.concatenate([Q, nq]); J2 = np.concatenate([J, np.concatenate(nj)])
        reg = jr[J2].astype(np.int64); movable = vol_job[J2]
        relax(lo, size, halo, *args, st_off, STC, STR, regmask, typ2, P2, Q2, reg, movable, OC, OR, n_iter,
              H_MIN, STEP_MAX, ROT_MAX, WALL_EVERY)
        cs = len(OR) + int(np.diff(st_off)[typ2].sum()) + 1000
        keep = accept(lo, size, halo, *args, st_off, STC, STR, BR, regmask, typ2, P2, Q2, reg, movable, OC, OR, cs, cs * 14,
                      bm, bw)
        if EXC is not None and len(EXC):                  # exclusion zones also hold after the squeeze
            bits = EXC_TYPE[typ2]
            fl = (bits != 0) & keep & movable
            for i in np.where(fl)[0]:
                d2 = ((EXC[:, :3] - P2[i]) ** 2).sum(1)
                hit = (d2 < EXC[:, 3] ** 2) & (((bits[i] >> EXC[:, 4].astype(np.int64)) & 1) == 1)
                if hit.any():
                    keep[i] = False
        if keep.sum() <= len(typ):
            break                                         # no gain: keep the previous configuration
        typ, P, Q, J = typ2[keep], P2[keep], Q2[keep], J2[keep]
    o = np.argsort(J, kind="stable")
    return (typ[o].astype(types.dtype), P[o].astype(pos.dtype), Q[o].astype(quat.dtype),
            np.bincount(J, minlength=len(jl)).astype(placed.dtype))
