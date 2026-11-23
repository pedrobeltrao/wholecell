# origin: jurkat/yc/pipeline/pack3d.py (legacy workspace, 2026-11-15)
# Changed: 'import cell3d as c3' -> 'from wholecell.geometry import engine as c3'; engine constants that became per-package
# values in G are read with the engine accessors inside numba code (c3.HALF -> c3.half(G), c3.PERI_W -> c3.peri_w(G),
# c3.PORE_R -> c3.pore_r(G), c3.min_mem_dist(out, s) -> c3.min_mem_dist_g(out, s, G)).
# A7 (2026-11-16): pack_tile takes four more arguments, per-region try budgets (REG_TRIES, REG_GIVEUP) and exclusion
# zones (EXC, EXC_TYPE), ported from the yeast packer; empty arrays = the v7 kernel (no extra random draws).
# A7 phase 2 (2026-11-17): a partner-attraction proposal for a shell job keeps its point and is tested against the
# shell band (origin: yeast pack3d.py:pack_tile band_m / band_w); only reached with a partner CSR (v7 passes none).
# A7 yeast side (2026-11-17): BAND_M / BAND_W: the same band test for volume-measured shells (yeast bands); empty = v7.
# perox-fix (2026-11-17): MCODE, the package's membrane codes of the membrane-protein tests (outer wall, its exempt
# membrane, pore membrane) instead of the sdf_sections_v1 constants c3.M_OUTER / c3.M_PM / c3.M_NE (yeast: code 7 =
# PEROX was taken as the outer wall, so no peroxisomal membrane protein could be placed); empty = those constants (v7).
# freeze fix (2026-11-19): partner-proposed membrane candidates must pass normal_ok (|grad| >= 0.5 and normal within 90
# degrees of the gradient, measured as the C9 validator does); only reached with a partner CSR (v7 unchanged).
"""3D random-sequential-addition packing of the whole cell in 200 nm tiles (numba core).

Collision geometry: per-ingredient sphere trees (st_c, st_r) on a local uniform grid (tile + halo).
Region tests: exact analytic SDFs from cell3d (label + distance to every membrane).
Membrane proteins: sampled on the tile's marching-cubes mesh by area, projected onto the exact surface,
library +z aligned to the SDF gradient (cytoplasmic / N-side), random spin.
Shell jobs (v5): soluble ingredients restricted to a band of signed distance from one membrane (lipid-anchored proteins,
cortical layer, lamina): candidates are sampled on that membrane's mesh, projected, offset along the normal into the band,
checked as volume jobs (label mask + membrane clearance) plus the band test on the centre; half of the tries align the
library's smallest principal axis with the membrane normal (as in thin lumens)."""
import numpy as np
import numba as nb
from wholecell.geometry import engine as c3

H_GRID = 6.0
TOL = 0.3


# ------------------------------------------------------------------ grid
@nb.njit(cache=True)
def grid_new(lo, hi, h, cap_s, cap_e):
    dims = np.empty(3, np.int64)
    for a in range(3):
        dims[a] = int(np.ceil((hi[a] - lo[a]) / h)) + 1
    head = -np.ones(dims[0] * dims[1] * dims[2], np.int64)
    nxt = np.empty(cap_e, np.int64); sid = np.empty(cap_e, np.int64)
    cen = np.empty((cap_s, 3), np.float64); rad = np.empty(cap_s, np.float64)
    cnt = np.zeros(2, np.int64)   # n_sph, n_ent
    return dims, head, nxt, sid, cen, rad, cnt


@nb.njit(cache=True, fastmath=True)
def _rng3(x, r, lo, h, d):
    i0 = int(np.floor((x - r - lo) / h)); i1 = int(np.floor((x + r - lo) / h))
    return max(i0, 0), min(i1, d - 1)


@nb.njit(cache=True, fastmath=True)
def g_overlap(lo, h, dims, head, nxt, sid, cen, rad, P, R, tol):
    for s in range(P.shape[0]):
        x, y, z, r = P[s, 0], P[s, 1], P[s, 2], R[s]
        i0, i1 = _rng3(x, r, lo[0], h, dims[0]); j0, j1 = _rng3(y, r, lo[1], h, dims[1]); k0, k1 = _rng3(z, r, lo[2], h, dims[2])
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    e = head[(i * dims[1] + j) * dims[2] + k]
                    while e >= 0:
                        q = sid[e]
                        dx = x - cen[q, 0]; dy = y - cen[q, 1]; dz = z - cen[q, 2]
                        rr = r + rad[q] - tol
                        if rr > 0 and dx * dx + dy * dy + dz * dz < rr * rr:
                            return True
                        e = nxt[e]
    return False


@nb.njit(cache=True, fastmath=True)
def g_insert(lo, h, dims, head, nxt, sid, cen, rad, cnt, P, R):
    for s in range(P.shape[0]):
        x, y, z, r = P[s, 0], P[s, 1], P[s, 2], R[s]
        i0, i1 = _rng3(x, r, lo[0], h, dims[0]); j0, j1 = _rng3(y, r, lo[1], h, dims[1]); k0, k1 = _rng3(z, r, lo[2], h, dims[2])
        if i1 < i0 or j1 < j0 or k1 < k0:
            continue
        need = (i1 - i0 + 1) * (j1 - j0 + 1) * (k1 - k0 + 1)
        if cnt[0] + 1 > cen.shape[0] or cnt[1] + need > nxt.shape[0]:
            return False
        q = cnt[0]; cnt[0] += 1
        cen[q, 0] = x; cen[q, 1] = y; cen[q, 2] = z; rad[q] = r
        for i in range(i0, i1 + 1):
            for j in range(j0, j1 + 1):
                for k in range(k0, k1 + 1):
                    c = (i * dims[1] + j) * dims[2] + k
                    e = cnt[1]; cnt[1] += 1
                    sid[e] = q; nxt[e] = head[c]; head[c] = e
    return True


# ------------------------------------------------------------------ rotations
@nb.njit(cache=True, fastmath=True)
def quat_mat(q, M):
    a, b, c, d = q[0], q[1], q[2], q[3]
    M[0, 0] = a*a+b*b-c*c-d*d; M[0, 1] = 2*(b*c-a*d); M[0, 2] = 2*(b*d+a*c)
    M[1, 0] = 2*(b*c+a*d); M[1, 1] = a*a-b*b+c*c-d*d; M[1, 2] = 2*(c*d-a*b)
    M[2, 0] = 2*(b*d-a*c); M[2, 1] = 2*(c*d+a*b); M[2, 2] = a*a-b*b-c*c+d*d


@nb.njit(cache=True)
def rand_quat(q):
    s = 0.0
    for i in range(4):
        q[i] = np.random.normal(); s += q[i] * q[i]
    s = np.sqrt(s)
    for i in range(4):
        q[i] /= s
    if q[0] < 0:
        for i in range(4):
            q[i] = -q[i]


@nb.njit(cache=True)
def quat_z_to(n, spin, q):
    """quaternion = (rotation taking +z to unit n) * (rotation by spin about z)."""
    # q1: z -> n
    cx, cy, cz = -n[1], n[0], 0.0          # z x n
    w = 1.0 + n[2]
    if w < 1e-8:
        a1, b1, c1, d1 = 0.0, 1.0, 0.0, 0.0
    else:
        s = np.sqrt(w * w + cx * cx + cy * cy)
        a1, b1, c1, d1 = w / s, cx / s, cy / s, cz / s
    a2, b2, c2, d2 = np.cos(spin / 2), 0.0, 0.0, np.sin(spin / 2)
    q[0] = a1*a2 - b1*b2 - c1*c2 - d1*d2
    q[1] = a1*b2 + b1*a2 + c1*d2 - d1*c2
    q[2] = a1*c2 - b1*d2 + c1*a2 + d1*b2
    q[3] = a1*d2 + b1*c2 - c1*b2 + d1*a2
    if q[0] < 0:
        for i in range(4):
            q[i] = -q[i]


@nb.njit(cache=True, fastmath=True)
def transform(stc, M, t, P):
    for s in range(stc.shape[0]):
        for a in range(3):
            P[s, a] = M[a, 0] * stc[s, 0] + M[a, 1] * stc[s, 1] + M[a, 2] * stc[s, 2] + t[a]


# ------------------------------------------------------------------ region tests
@nb.njit(cache=True, fastmath=True)
def check_volume(P, R, cx, cy, cz, br, mask, G, pores, pidx, PT, PQ, PC, out):
    m = c3.mem_sdfs(cx, cy, cz, G, pores, pidx, PT, PQ, PC, out)
    l = c3.label_from(cx, cy, cz, G, pores, out, m)
    if not (mask >> l) & 1:
        return False
    if c3.min_mem_dist_g(out, -1, G) >= br + c3.half(G):
        return True
    for s in range(P.shape[0]):
        m = c3.mem_sdfs(P[s, 0], P[s, 1], P[s, 2], G, pores, pidx, PT, PQ, PC, out)
        l = c3.label_from(P[s, 0], P[s, 1], P[s, 2], G, pores, out, m)
        if not (mask >> l) & 1:
            return False
        if c3.min_mem_dist_g(out, -1, G) < R[s] + c3.half(G):
            return False
    return True


@nb.njit(cache=True, fastmath=True)
def check_membrane(P, R, ZL, own, G, pores, pidx, PT, PQ, PC, out, outer, exempt):
    """outer: code of the outer wall membrane (-1: none); exempt: the membrane whose proteins span that wall (-1: none),
    tested against the open face instead (Jurkat: CELL_OUTER and PM; yeast: CELL_OUTER and none)"""
    for s in range(P.shape[0]):
        c3.mem_sdfs(P[s, 0], P[s, 1], P[s, 2], G, pores, pidx, PT, PQ, PC, out)
        if c3.min_mem_dist_g(out, own, G) < R[s] + c3.half(G):
            return False
        if own != exempt and outer >= 0 and out[outer] < R[s]:   # outside the cell (PM proteins span it: exo domains are extracellular)
            return False
        if own == exempt and c3.face_dist(P[s, 2], G) < R[s]:   # open face: no PM, nothing beyond the imaged volume
            return False
        z = ZL[s]
        if abs(z) >= R[s] + c3.half(G):
            d = out[own]
            if d * z <= 0 or abs(d) < R[s] + c3.half(G) - 1.0:
                return False
    return True


@nb.njit(cache=True, fastmath=True)
def mem_value(px, py, pz, m, G, pores, pidx, PT, PQ, PC, out):
    c3.mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out)
    return out[m]


@nb.njit(cache=True, fastmath=True)
def project_to_membrane(p, m, n, G, pores, pidx, PT, PQ, PC, out):
    """Newton projection onto the zero set of membrane m; n receives the unit gradient (cytoplasmic side)."""
    e = 0.4
    for it in range(4):
        d = mem_value(p[0], p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out)
        gx = (mem_value(p[0] + e, p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out) - mem_value(p[0] - e, p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out)) / (2 * e)
        gy = (mem_value(p[0], p[1] + e, p[2], m, G, pores, pidx, PT, PQ, PC, out) - mem_value(p[0], p[1] - e, p[2], m, G, pores, pidx, PT, PQ, PC, out)) / (2 * e)
        gz = (mem_value(p[0], p[1], p[2] + e, m, G, pores, pidx, PT, PQ, PC, out) - mem_value(p[0], p[1], p[2] - e, m, G, pores, pidx, PT, PQ, PC, out)) / (2 * e)
        gn = np.sqrt(gx * gx + gy * gy + gz * gz)
        if gn < 1e-6:
            return False
        n[0] = gx / gn; n[1] = gy / gn; n[2] = gz / gn
        p[0] -= d * n[0]; p[1] -= d * n[1]; p[2] -= d * n[2]
        if abs(d) < 0.02:
            break
    return abs(mem_value(p[0], p[1], p[2], m, G, pores, pidx, PT, PQ, PC, out)) < 0.1


# freeze fix (2026-11-19): partner proposals of the organisation layer start a membrane candidate anywhere around a
# placed partner, so Newton projection can end where the membrane SDF is not a distance (gradient far below 1: the
# medial ridges and seams of the min / max composed sections) and its normal is ill-defined. The validator (C9,
# wholecell/validate/placement.py:orientation) judges the sign of the library +z against the SDF gradient at the stored
# (float32) centre by central differences with h = 0.5 nm; the packer's normal is the gradient at the last Newton
# iterate (h = 0.4), which can differ there (Jurkat combined build: 3 flips at validator |grad| 0.52-0.73; subset jb2:
# 9 of 10 flips at 0.02-0.47). The test below evaluates the gradient as the validator does (same stencil, at the
# float32-rounded centre) and rejects the candidate when |grad| < GRAD_MIN or when the placement normal is more than
# 90 degrees from that gradient (a sign flip by the validator's measure). Applied to partner proposals only (on-mesh
# proposals start within a voxel of the surface; v7, without a partner CSR, stays bitwise unchanged); no random draws.
GRAD_MIN = 0.5
GRAD_H = 0.5


@nb.njit(cache=True, fastmath=True)
def normal_ok(p, m, n, G, pores, pidx, PT, PQ, PC, out):
    """|grad SDF_m| >= GRAD_MIN at the float32-rounded centre p (central differences, h = GRAD_H, as the validator) and
    n . grad > 0 (unit normal n within 90 degrees of the gradient)"""
    x = float(np.float32(p[0])); y = float(np.float32(p[1])); z = float(np.float32(p[2])); h = GRAD_H
    gx = (mem_value(x + h, y, z, m, G, pores, pidx, PT, PQ, PC, out) - mem_value(x - h, y, z, m, G, pores, pidx, PT, PQ, PC, out)) / (2 * h)
    gy = (mem_value(x, y + h, z, m, G, pores, pidx, PT, PQ, PC, out) - mem_value(x, y - h, z, m, G, pores, pidx, PT, PQ, PC, out)) / (2 * h)
    gz = (mem_value(x, y, z + h, m, G, pores, pidx, PT, PQ, PC, out) - mem_value(x, y, z - h, m, G, pores, pidx, PT, PQ, PC, out)) / (2 * h)
    if gx * gx + gy * gy + gz * gz < GRAD_MIN * GRAD_MIN:
        return False
    return n[0] * gx + n[1] * gy + n[2] * gz > 0.0


@nb.njit(cache=True, fastmath=True)
def near_pore(p, G, pores, lim):
    return c3.in_pore(p[0], p[1], p[2], pores, lim) >= 0


# ------------------------------------------------------------------ tile packer
KRES = 32     # anchors kept per library type (uniform reservoir sample of everything placed in the tile + halo)


@nb.njit(cache=True)
def _res_add(res_pos, res_n, ty, p):
    n = res_n[ty]
    slot = n if n < KRES else int(np.random.random() * (n + 1))
    if slot < KRES:
        for k in range(3):
            res_pos[ty, slot, k] = p[k]
    res_n[ty] = n + 1


@nb.njit(cache=True)
def _near(res_pos, res_n, pt, d, t):
    """t = a random kept anchor of type pt + d along a random direction"""
    slot = int(np.random.random() * min(res_n[pt], KRES))
    dx = np.random.normal(); dy = np.random.normal(); dz = np.random.normal()
    dn = np.sqrt(dx * dx + dy * dy + dz * dz) + 1e-12
    t[0] = res_pos[pt, slot, 0] + d * dx / dn; t[1] = res_pos[pt, slot, 1] + d * dy / dn; t[2] = res_pos[pt, slot, 2] + d * dz / dn


LUMEN_GIVEUP_F = 3


# origin: yeast pack3d.py:_in_exc @ e3b7448
# Changed: zones carry a rule id (column 4) and a type sees only the rules in its bitmask (several exclusion rules).
@nb.njit(cache=True)
def _in_exc(t, EXC, bits):
    """True if point t lies inside an exclusion sphere (EXC rows: x, y, z, r, rule) of a rule set in bitmask `bits`"""
    for i in range(EXC.shape[0]):
        if (bits >> int(EXC[i, 4])) & 1:
            d = (t[0] - EXC[i, 0]) ** 2 + (t[1] - EXC[i, 1]) ** 2 + (t[2] - EXC[i, 2]) ** 2
            if d < EXC[i, 3] * EXC[i, 3]:
                return True
    return False


@nb.njit(cache=True, fastmath=True)
def sample_tri(t0, t1, CUMA, TRI, t):
    """random point on triangles [t0, t1) of the tile mesh, by area"""
    u = np.random.random() * CUMA[t1 - 1] if t0 == 0 else CUMA[t0 - 1] + np.random.random() * (CUMA[t1 - 1] - CUMA[t0 - 1])
    lo_i = t0; hi_i = t1 - 1
    while lo_i < hi_i:
        mid = (lo_i + hi_i) // 2
        if CUMA[mid] < u:
            lo_i = mid + 1
        else:
            hi_i = mid
    a = np.random.random(); b = np.random.random()
    if a + b > 1:
        a = 1 - a; b = 1 - b
    for k in range(3):
        t[k] = TRI[lo_i, k] + a * (TRI[lo_i, 3 + k] - TRI[lo_i, k]) + b * (TRI[lo_i, 6 + k] - TRI[lo_i, k])


@nb.njit(cache=True, fastmath=True)
def lumen_quat(t, m, qf, G, pores, pidx, PT, PQ, PC, out, nrm, q):
    """q = (z -> local normal of membrane m at t, random spin) * qf; False if the normal is undefined"""
    e = 1.0; n2 = 0.0
    for k in range(3):
        t[k] += e; a = mem_value(t[0], t[1], t[2], m, G, pores, pidx, PT, PQ, PC, out)
        t[k] -= 2 * e; b = mem_value(t[0], t[1], t[2], m, G, pores, pidx, PT, PQ, PC, out)
        t[k] += e; nrm[k] = a - b; n2 += nrm[k] * nrm[k]
    if n2 < 1e-8:
        return False
    n2 = np.sqrt(n2)
    for k in range(3):
        nrm[k] /= n2
    qa = np.empty(4)
    quat_z_to(nrm, np.random.random() * 2 * np.pi, qa)
    q[0] = qa[0]*qf[0] - qa[1]*qf[1] - qa[2]*qf[2] - qa[3]*qf[3]
    q[1] = qa[0]*qf[1] + qa[1]*qf[0] + qa[2]*qf[3] - qa[3]*qf[2]
    q[2] = qa[0]*qf[2] - qa[1]*qf[3] + qa[2]*qf[0] + qa[3]*qf[1]
    q[3] = qa[0]*qf[3] + qa[1]*qf[2] - qa[2]*qf[1] + qa[3]*qf[0]
    if q[0] < 0:
        for i in range(4):
            q[i] = -q[i]
    return True


@nb.njit(cache=True)
def pack_tile(seed, lo, size, halo,
              G, pores, pidx, PT, PQ, PC,
              st_off, STC, STR, BR, regmask,
              jl, jr, jn, n_mem_regions0,
              vox_off, VOXP, vox_h,
              tri_off, TRI, CUMA,
              OC, OR,
              n_mrna, mrna_mean, cyto_region,
              max_tries, give_up, cap_s, cap_e,
              AT, AP, ptn_off, ptn_idx, ptn_str, ptn_kind, prm, QFLAT, REGMEM,
              n_shell0, SH_MEM, SH_MASK, SH_SIGN, SH_D0, SH_D1, SH_RAD,
              REG_TRIES, REG_GIVEUP, EXC, EXC_TYPE, BAND_M, BAND_W, MCODE):
    """... AT, AP: types and positions of instances already placed around the tile (anchors); ptn_off/ptn_idx/ptn_str/ptn_kind:
    partner CSR over library types and prm (A_contact, A_loose, K_half, P_max, attract tries, loose gap min, max),
    see organisation.py. Empty ptn_off disables attraction.
    QFLAT (n_types, 4): rotation taking each library's smallest principal axis to +z; REGMEM (n volume regions): membrane
    whose normal a thin lumen follows (-1 none). In such lumens half of the tries align the smallest axis with the local
    membrane normal (random spin), and a job gives up after LUMEN_GIVEUP_F x give_up failures in a row. Empty QFLAT
    disables this.
    Regions >= n_shell0 are shells (SH_RAD bit 0: add the bounding radius to the outer bound; SH_RAD >> 1: pore-axis limit, nm): shell s = reg - n_shell0 places a soluble ingredient with its centre at signed distance
    SH_SIGN[s] * u, u in [SH_D0[s], SH_D1[s] (+ BR if SH_RAD[s])], from membrane SH_MEM[s] (SDF sign convention of cell3d:
    + = cytoplasmic / nucleoplasmic / matrix side), all beads in labels SH_MASK[s] (bitmask) with membrane clearance.
    A7 (origin: yeast pack3d.py:pack_tile @ e3b7448, reg_tries / reg_give_up / EXC / exc_type): REG_TRIES, REG_GIVEUP
    (per region, 0 = as above): tries per copy and failures in a row, replacing the global and thin-lumen budgets;
    EXC (m, 5: centre, radius, rule id), EXC_TYPE (per library type, bitmask of rules): soluble copies (volume and
    shell jobs) of a type may not have their centre in a sphere of one of its rules. Empty arrays disable both.
    A7 yeast side (origin: yeast pack3d.py:pack_tile band_m / band_w): BAND_M, BAND_W (per region; empty = none): a
    partner proposal for a volume job of a region with BAND_M[reg] >= 0 (a shell measured as a volume, packed from its
    voxels) must lie within BAND_W[reg] of membrane BAND_M[reg] (|SDF|).
    MCODE (3 engine membrane codes of the package, -1 = none; empty = the sdf_sections_v1 constants M_OUTER, M_PM, M_NE):
    outer wall (no membrane-protein bead beyond it), the membrane exempt from that wall (bounded by the open face
    instead), the pore membrane (its proteins keep clear of the pores). Regions.membrane_codes builds it."""
    np.random.seed(seed)
    n_types = BR.shape[0]
    res_pos = np.empty((n_types, KRES, 3)); res_n = np.zeros(n_types, np.int64)
    for ia in range(AT.shape[0]):
        _res_add(res_pos, res_n, AT[ia], AP[ia])
    use_att = ptn_off.shape[0] == n_types + 1
    wk = np.empty(256); att_stats = np.zeros(2, np.int64)
    att_tries = int(prm[4]) if use_att else 0
    glo = lo - halo; ghi = lo + size + halo
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    ok = g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for j in range(jn.shape[0]):
        ntot += jn[j]
    otype = np.empty(ntot, np.int32); opos = np.empty((ntot, 3), np.float32); oquat = np.empty((ntot, 4), np.float32)
    placed = np.zeros(jn.shape[0], np.int64)
    out = np.zeros(c3.NMEM_MAX)                       # A7: packages with more membranes (yeast 10)
    m_outer, m_exempt, m_pore = c3.M_OUTER, c3.M_PM, c3.M_NE
    if MCODE.shape[0] >= 3:                           # package membrane codes (Regions.membrane_codes)
        m_outer, m_exempt, m_pore = MCODE[0], MCODE[1], MCODE[2]
    M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3); nrm = np.empty(3)
    pk_list = np.zeros(max(pores.shape[0], 1), np.int64)
    no = 0
    overflow = not ok
    for j in range(jn.shape[0]):
        if overflow:
            break
        li = jl[j]; reg = jr[j]
        a0 = st_off[li]; a1 = st_off[li + 1]
        stc = STC[a0:a1]; str_ = STR[a0:a1]
        P = np.empty((a1 - a0, 3)); ZL = np.empty(a1 - a0)
        is_shell = reg >= n_shell0
        is_mem = reg >= n_mem_regions0 and not is_shell
        mi = reg - n_mem_regions0
        sm = -1; sa = 0.0; sb = 0.0; ssg = 1.0; si = 0; plim = 0.0
        v0 = 0; v1 = 0; t0 = 0; t1 = 0; npk = 0
        if is_shell:
            si = reg - n_shell0; sm = SH_MEM[si]; ssg = SH_SIGN[si]
            sa = SH_D0[si]; sb = SH_D1[si] + (BR[li] if (SH_RAD[si] & 1) else 0.0)
            plim = float(SH_RAD[si] >> 1)                    # v6: > 0 = centre within plim nm of a pore axis (NPC shell)
            t0 = tri_off[sm]; t1 = tri_off[sm + 1]
            if t1 <= t0:
                continue
            # v7: narrow pore-limited shells (central channel, plim <= 60 nm) sample their centre directly in the pore
            # cylinder (radius plim, |z| <= sb + PERI_W/2 along the pore axis) of a pore near the tile; membrane-triangle
            # sampling almost never lands within a few tens of nm of an axis
            npk = 0
            if plim > 0.0 and plim <= 60.0:
                for k in range(pores.shape[0]):
                    inside = True
                    for d in range(3):
                        if pores[k, d] < lo[d] - plim - sb or pores[k, d] > lo[d] + size + plim + sb:
                            inside = False
                    if inside:
                        pk_list[npk] = k; npk += 1
                if npk == 0:
                    continue
        elif is_mem:
            t0 = tri_off[mi]; t1 = tri_off[mi + 1]
            if t1 <= t0:
                continue
            for s in range(a1 - a0):
                ZL[s] = stc[s, 2]
        else:
            v0 = vox_off[reg]; v1 = vox_off[reg + 1]
            if v1 <= v0:
                continue
        thin = (not is_mem) and (not is_shell) and QFLAT.shape[0] == n_types and REGMEM[reg] >= 0
        gu = give_up * LUMEN_GIVEUP_F if (thin or is_shell) else give_up
        mt_ = max_tries
        if thin and prm.shape[0] > 8:                        # v7: larger attempt budget in thin lumens (prm[7] tries x,
            mt_ = int(max_tries * prm[7]); gu = int(gu * prm[8])   # prm[8] give-up x): RSA there is attempt-limited
        if REG_TRIES.shape[0] > reg and REG_TRIES[reg] > 0:   # A7 per-region budgets (origin: yeast pack3d.py:pack_tile)
            mt_ = REG_TRIES[reg]
        if REG_GIVEUP.shape[0] > reg and REG_GIVEUP[reg] > 0:
            gu = REG_GIVEUP[reg]
        exb = EXC_TYPE[li] if EXC.shape[0] > 0 else 0          # A7 exclusion rules of this type (0: none)
        # v7 seed pruning in thin lumens: VOXP column 3 = membrane clearance of the voxel (seeds of the region sorted by it,
        # descending) and QFLAT column 4 = the library's minimal support half-width (min over directions of the extent from
        # its origin). A centre closer to a membrane than that cannot be valid in any orientation, so seeds are drawn only
        # from voxels with clearance >= QFLAT[li, 4] - VOX_JIT (jitter margin). Without the extra columns: as v6.
        vcut = v1 if not (is_mem or is_shell) else 0
        if thin and VOXP.shape[1] >= 4 and QFLAT.shape[1] >= 5:
            need = QFLAT[li, 4] - 0.87 * vox_h
            lo_i = v0; hi_i = v1
            while lo_i < hi_i:                                # first index with clearance < need (descending order)
                mid = (lo_i + hi_i) // 2
                if VOXP[mid, 3] >= need:
                    lo_i = mid + 1
                else:
                    hi_i = mid
            vcut = lo_i if lo_i > v0 else v1
        fails_in_row = 0
        for inst in range(jn[j]):
            success = False
            att_k = -1
            if use_att and ptn_off[li + 1] > ptn_off[li]:
                tot = 0.0
                for k in range(ptn_off[li], ptn_off[li + 1]):
                    nn = float(res_n[ptn_idx[k]])
                    a_ = ptn_str[k] * prm[ptn_kind[k]] * nn / (prm[2] + nn)
                    wk[k - ptn_off[li]] = a_; tot += a_
                if tot > 0 and np.random.random() < min(prm[3], tot):
                    x = np.random.random() * tot
                    for k in range(ptn_off[li], ptn_off[li + 1]):
                        x -= wk[k - ptn_off[li]]
                        if x <= 0:
                            att_k = k
                            break
                    if att_k < 0:
                        att_k = ptn_off[li + 1] - 1
                    att_stats[0] += 1
            for tr in range(mt_):
                near = att_k >= 0 and tr < att_tries
                if near:
                    pt = ptn_idx[att_k]
                    if ptn_kind[att_k] == 0:
                        dd = (BR[li] + BR[pt]) * (0.45 + 0.5 * np.random.random())
                    else:
                        dd = BR[li] + BR[pt] + prm[5] + (prm[6] - prm[5]) * np.random.random()
                    _near(res_pos, res_n, pt, dd, t)
                if is_shell and npk > 0:
                    kp = pk_list[int(np.random.random() * npk)]
                    zz = (2.0 * np.random.random() - 1.0) * (sb + 0.5 * c3.peri_w(G))
                    rr = plim * np.sqrt(np.random.random()); ph = 2.0 * np.pi * np.random.random()
                    ax_ = pores[kp, 3:6]
                    e1x = -ax_[1]; e1y = ax_[0]; e1z = 0.0
                    if abs(ax_[2]) > 0.9:
                        e1x = 0.0; e1y = -ax_[2]; e1z = ax_[1]
                    nn_ = np.sqrt(e1x * e1x + e1y * e1y + e1z * e1z)
                    e1x /= nn_; e1y /= nn_; e1z /= nn_
                    e2x = ax_[1] * e1z - ax_[2] * e1y; e2y = ax_[2] * e1x - ax_[0] * e1z; e2z = ax_[0] * e1y - ax_[1] * e1x
                    cph = np.cos(ph) * rr; sph = np.sin(ph) * rr
                    t[0] = pores[kp, 0] + zz * ax_[0] + cph * e1x + sph * e2x
                    t[1] = pores[kp, 1] + zz * ax_[1] + cph * e1y + sph * e2y
                    t[2] = pores[kp, 2] + zz * ax_[2] + cph * e1z + sph * e2z
                    if t[0] < lo[0] or t[1] < lo[1] or t[2] < lo[2] or t[0] >= lo[0] + size or t[1] >= lo[1] + size or t[2] >= lo[2] + size:
                        continue
                    u_ = ssg * mem_value(t[0], t[1], t[2], sm, G, pores, pidx, PT, PQ, PC, out)
                    if u_ < sa or u_ > sb:
                        continue
                    if exb != 0 and _in_exc(t, EXC, exb):
                        continue
                    rand_quat(q)
                    quat_mat(q, M)
                    transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_volume(P, str_, t[0], t[1], t[2], BR[li], SH_MASK[si], G, pores, pidx, PT, PQ, PC, out):
                        continue
                elif is_shell:
                    if not near:                             # A7 P6: a partner proposal keeps its point (band test below)
                        sample_tri(t0, t1, CUMA, TRI, t)
                        if not project_to_membrane(t, sm, nrm, G, pores, pidx, PT, PQ, PC, out):
                            continue
                        off = ssg * (sa + np.random.random() * (sb - sa))
                        for k in range(3):
                            t[k] += off * nrm[k]
                    if t[0] < lo[0] or t[1] < lo[1] or t[2] < lo[2] or t[0] >= lo[0] + size or t[1] >= lo[1] + size or t[2] >= lo[2] + size:
                        continue
                    u_ = ssg * mem_value(t[0], t[1], t[2], sm, G, pores, pidx, PT, PQ, PC, out)
                    if u_ < sa or u_ > sb:
                        continue
                    if plim > 0.0 and not near_pore(t, G, pores, plim):
                        continue
                    if exb != 0 and _in_exc(t, EXC, exb):
                        continue
                    if not (QFLAT.shape[0] == n_types and np.random.random() < 0.5 and lumen_quat(t, sm, QFLAT[li], G, pores, pidx, PT, PQ, PC, out, nrm, q)):
                        rand_quat(q)
                    quat_mat(q, M)
                    transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_volume(P, str_, t[0], t[1], t[2], BR[li], SH_MASK[si], G, pores, pidx, PT, PQ, PC, out):
                        continue
                elif is_mem and near:
                    pass
                elif is_mem:
                    u = np.random.random() * CUMA[t1 - 1] if t0 == 0 else CUMA[t0 - 1] + np.random.random() * (CUMA[t1 - 1] - CUMA[t0 - 1])
                    # binary search in [t0, t1)
                    lo_i = t0; hi_i = t1 - 1
                    while lo_i < hi_i:
                        mid = (lo_i + hi_i) // 2
                        if CUMA[mid] < u:
                            lo_i = mid + 1
                        else:
                            hi_i = mid
                    a = np.random.random(); b = np.random.random()
                    if a + b > 1:
                        a = 1 - a; b = 1 - b
                    for k in range(3):
                        t[k] = TRI[lo_i, k] + a * (TRI[lo_i, 3 + k] - TRI[lo_i, k]) + b * (TRI[lo_i, 6 + k] - TRI[lo_i, k])
                if is_shell:
                    pass
                elif is_mem:
                    if not project_to_membrane(t, mi, nrm, G, pores, pidx, PT, PQ, PC, out):
                        continue
                    if t[0] < lo[0] or t[1] < lo[1] or t[2] < lo[2] or t[0] >= lo[0] + size or t[1] >= lo[1] + size or t[2] >= lo[2] + size:
                        continue
                    if mi == m_pore and near_pore(t, G, pores, c3.pore_r(G) + 0.5 * c3.peri_w(G) + 25.0):
                        continue
                    if near and not normal_ok(t, mi, nrm, G, pores, pidx, PT, PQ, PC, out):   # freeze fix: see normal_ok
                        continue
                    quat_z_to(nrm, np.random.random() * 2 * np.pi, q)
                    quat_mat(q, M)
                    transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_membrane(P, str_, ZL, mi, G, pores, pidx, PT, PQ, PC, out, m_outer, m_exempt):
                        continue
                else:
                    if near:
                        if t[0] < lo[0] or t[1] < lo[1] or t[2] < lo[2] or t[0] >= lo[0] + size or t[1] >= lo[1] + size or t[2] >= lo[2] + size:
                            continue
                        if BAND_M.shape[0] > reg and BAND_M[reg] >= 0:
                            if abs(mem_value(t[0], t[1], t[2], BAND_M[reg], G, pores, pidx, PT, PQ, PC, out)) > BAND_W[reg]:
                                continue
                    else:
                        vi = v0 + int(np.random.random() * (vcut - v0))
                        for k in range(3):
                            t[k] = VOXP[vi, k] + (np.random.random() - 0.5) * vox_h
                    if exb != 0 and _in_exc(t, EXC, exb):
                        continue
                    if thin and np.random.random() < 0.5 and lumen_quat(t, REGMEM[reg], QFLAT[li], G, pores, pidx, PT, PQ, PC, out, nrm, q):
                        pass
                    else:
                        rand_quat(q)
                    quat_mat(q, M)
                    transform(stc, M, t, P)
                    if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL):
                        continue
                    if not check_volume(P, str_, t[0], t[1], t[2], BR[li], regmask[reg], G, pores, pidx, PT, PQ, PC, out):
                        continue
                if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                    overflow = True
                    break
                otype[no] = li
                for k in range(3):
                    opos[no, k] = t[k]
                for k in range(4):
                    oquat[no, k] = q[k]
                no += 1; placed[j] += 1
                if use_att:
                    _res_add(res_pos, res_n, li, t)
                    if near:
                        att_stats[1] += 1
                success = True
                break
            if overflow:
                break
            if success:
                fails_in_row = 0
            else:
                fails_in_row += 1
                if fails_in_row >= gu:
                    break
    # ---------------- mRNA threads (1 bead / nt, 0.6 nm rise) through free cytosol
    mr_pts = np.empty((0, 3), np.float32)
    mr_off = np.zeros(1, np.int64)
    if n_mrna > 0 and not overflow:
        v0 = vox_off[cyto_region]; v1 = vox_off[cyto_region + 1]
        maxb = 0
        lens = np.empty(n_mrna, np.int64)
        for c in range(n_mrna):
            lens[c] = max(300, int(np.random.exponential(mrna_mean)))
            maxb += lens[c]
        mr_pts = np.empty((maxb, 3), np.float32)
        mr_off = np.zeros(n_mrna + 1, np.int64)
        nb_ = 0
        bead = np.empty((1, 3)); br_ = np.full(1, 0.5)
        prev = np.empty((1, 3))
        mlo = lo - 10.0; mhi = lo + size + 10.0
        for c in range(n_mrna):
            if v1 <= v0:
                mr_off[c + 1] = nb_
                continue
            start_ok = False
            for at in range(200):
                vi = v0 + int(np.random.random() * (v1 - v0))
                for k in range(3):
                    bead[0, k] = VOXP[vi, k] + (np.random.random() - 0.5) * vox_h
                if not g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, bead, br_, 0.0):
                    m = c3.mem_sdfs(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, pidx, PT, PQ, PC, out)
                    if c3.label_from(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, out, m) == c3.CYTO and c3.min_mem_dist_g(out, -1, G) > 0.5 + c3.half(G):
                        start_ok = True
                        break
            if not start_ok:
                mr_off[c + 1] = nb_
                continue
            for k in range(3):
                mr_pts[nb_, k] = bead[0, k]
            nb_ += 1; nthis = 1
            dx = np.random.normal(); dy = np.random.normal(); dz = np.random.normal()
            dn = np.sqrt(dx * dx + dy * dy + dz * dz); dx /= dn; dy /= dn; dz /= dn
            stuck = 0
            pending = False
            while nthis < lens[c] and stuck < 200:
                ex = dx + np.random.normal() * 0.9; ey = dy + np.random.normal() * 0.9; ez = dz + np.random.normal() * 0.9
                en = np.sqrt(ex * ex + ey * ey + ez * ez); ex /= en; ey /= en; ez /= en
                bead[0, 0] = mr_pts[nb_ - 1, 0] + 0.6 * ex
                bead[0, 1] = mr_pts[nb_ - 1, 1] + 0.6 * ey
                bead[0, 2] = mr_pts[nb_ - 1, 2] + 0.6 * ez
                if bead[0, 0] < mlo[0] or bead[0, 1] < mlo[1] or bead[0, 2] < mlo[2] or bead[0, 0] > mhi[0] or bead[0, 1] > mhi[1] or bead[0, 2] > mhi[2]:
                    stuck += 1; continue
                if g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, bead, br_, 0.0):
                    stuck += 1; continue
                m = c3.mem_sdfs(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, pidx, PT, PQ, PC, out)
                if c3.label_from(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, out, m) != c3.CYTO or c3.min_mem_dist_g(out, -1, G) < 0.5 + c3.half(G):
                    stuck += 1; continue
                # insert the previous bead (lagged so consecutive beads do not self-collide)
                if pending:
                    if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, prev, br_):
                        overflow = True; break
                for k in range(3):
                    prev[0, k] = mr_pts[nb_ - 1, k]
                pending = True
                for k in range(3):
                    mr_pts[nb_, k] = bead[0, k]
                nb_ += 1; nthis += 1
                dx, dy, dz = ex, ey, ez; stuck = 0
            if overflow:
                break
            # insert the last two beads
            for k in range(3):
                prev[0, k] = mr_pts[nb_ - 1, k]
            g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, prev, br_)
            if nthis > 1:
                for k in range(3):
                    prev[0, k] = mr_pts[nb_ - 2, k]
                g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, prev, br_)
            mr_off[c + 1] = nb_
        mr_pts = mr_pts[:nb_]
        mr_off = mr_off[:n_mrna + 1]
    return otype[:no], opos[:no], oquat[:no], placed, mr_pts, mr_off, overflow, cnt[0], cnt[1], att_stats


# ------------------------------------------------------------------ helpers for obstacles
@nb.njit(cache=True)
def emit_spheres(types, pos, quat, st_off, STC, STR, BR, lo, hi):
    """Sphere-tree spheres of instances whose bounding sphere intersects box [lo, hi]."""
    n = 0
    keep = np.zeros(types.shape[0], np.bool_)
    for i in range(types.shape[0]):
        r = BR[types[i]]
        ok = True
        for a in range(3):
            if pos[i, a] + r < lo[a] or pos[i, a] - r > hi[a]:
                ok = False
        if ok:
            keep[i] = True
            n += st_off[types[i] + 1] - st_off[types[i]]
    C = np.empty((n, 3)); R = np.empty(n)
    M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    k = 0
    for i in range(types.shape[0]):
        if not keep[i]:
            continue
        for a in range(4):
            q[a] = quat[i, a]
        quat_mat(q, M)
        for a in range(3):
            t[a] = pos[i, a]
        a0 = st_off[types[i]]; a1 = st_off[types[i] + 1]
        transform(STC[a0:a1], M, t, C[k:k + a1 - a0])
        R[k:k + a1 - a0] = STR[a0:a1]
        k += a1 - a0
    return C, R


# ------------------------------------------------------------------ chromatin pre-pass (global, nucleus only)
@nb.njit(cache=True)
def chromatin_walk(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask,
                   OC, OR, chain_len, starts, r_start, step, cap_s, cap_e):
    """Persistent self-avoiding nucleosome walks, one per chromosome, first fragment starting within r_start of
    starts[c] (the SPB or the chromosome's kinetochore-MT tip: Rabl-like). Returns positions, quaternions, chain id
    and a 'break' flag (1 = fragment start)."""
    np.random.seed(seed)
    dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
    g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
    ntot = 0
    for c in range(chain_len.shape[0]):
        ntot += chain_len[c]
    pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32)
    chain = np.empty(ntot, np.int32); brk = np.zeros(ntot, np.int8)
    out = np.zeros(c3.NMEM_MAX); M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
    P = np.empty((stc.shape[0], 3))
    n = 0
    for c in range(chain_len.shape[0]):
        placed = 0; have_prev = False
        dx, dy, dz = 0.0, 0.0, 0.0
        px, py, pz = 0.0, 0.0, 0.0
        fails = 0; restarts = 0
        while placed < chain_len[c] and restarts < 20000:
            if not have_prev:
                rr = r_start if placed == 0 else 900.0
                cxx = starts[c, 0] if placed == 0 else G[3]
                cyy = starts[c, 1] if placed == 0 else G[4]
                czz = starts[c, 2] if placed == 0 else G[5]
                t[0] = cxx + (np.random.random() * 2 - 1) * rr
                t[1] = cyy + (np.random.random() * 2 - 1) * rr
                t[2] = czz + (np.random.random() * 2 - 1) * rr
            else:
                ex = np.random.normal(); ey = np.random.normal(); ez = np.random.normal()
                en = np.sqrt(ex * ex + ey * ey + ez * ez)
                ex = 0.5 * dx + ex / en; ey = 0.5 * dy + ey / en; ez = 0.5 * dz + ez / en
                en = np.sqrt(ex * ex + ey * ey + ez * ez)
                ex /= en; ey /= en; ez /= en
                t[0] = px + step * ex; t[1] = py + step * ey; t[2] = pz + step * ez
            rand_quat(q); quat_mat(q, M); transform(stc, M, t, P)
            good = not g_overlap(glo, H_GRID, dims, head, nxt, sid, cen, rad, P, str_, TOL)
            if good:
                good = check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out)
            if not good:
                fails += 1
                if fails % 40 == 0 or not have_prev:
                    have_prev = False
                    restarts += 1
                continue
            if not g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                break
            if have_prev:
                dx = (t[0] - px) / step; dy = (t[1] - py) / step; dz = (t[2] - pz) / step
            else:
                brk[n] = 1
            px, py, pz = t[0], t[1], t[2]
            have_prev = True
            for k in range(3):
                pos[n, k] = t[k]
            for k in range(4):
                qu[n, k] = q[k]
            chain[n] = c
            n += 1; placed += 1; fails = 0
    return pos[:n], qu[:n], chain[:n], brk[:n]
