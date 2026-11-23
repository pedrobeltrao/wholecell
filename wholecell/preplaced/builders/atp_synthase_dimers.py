# origin: yeast pipeline/atp_dimers.py:_site_ok, place @ e3b7448 (rim placement: copied, geometry through ctx.K)
# origin: yeast scripts/s3_pack.py:main ("ATP synthase dimer rows on cristae rims" block) @ e3b7448
# Changed (unification B4, 2026-11-17): one core module, no ctx.mod. Two rim sources (params.rims): `prims` = the
# yeast crista prims (C_CRISTA pills of the yeast geometry; the yeast dimer library entry ATP_dimer made by
# atp_dimers.model, a library builder still run by the yeast context) - yeast output bitwise; `lamellae` = the
# lamellar cristae of a core geometry (section `cristae_lamellae_v1`, Jurkat), new. Yeast constants are params with
# the yeast values as defaults (config/cells/yeast.yaml unchanged).
# Changed (B4 J6 fix, 2026-11-18): lamellae rim_sampling `rows` (dense rim lines, tilt about the rim tangent, faster site
# kernel with the same test); the default `random` and the yeast prims path are unchanged.
"""ATP synthase dimer rows on the cristae rims (serial pre-placed builder `core.atp_synthase_dimers`).

Dimers sit on the rim of a crista, one every SPACING nm along it, ridge along the rim tangent: library x -> crista
normal, y -> rim tangent, z -> outward (the dimer's two membrane normals straddle the rim; Davies et al. 2012 PNAS
109:13602). A site is kept when the transmembrane residues lie within TM_TOL of the inner-membrane midsurface, every
third soluble residue is in the matrix or the IMS with at most MAX_IN_MEMBRANE of them inside a bilayer, and the dimer's
sphere tree does not overlap a dimer placed before. Pool: the inner-membrane ATP synthase monomers, at most half of
them in dimers (two monomers per dimer are deducted).

rims: prims (yeast)  the crista pills of the yeast geometry (cell.PC == C_CRISTA, cell.PQ: centre, normal, radius):
                     sites on the rim circle. The dimer is one library type (ATP_dimer). Adds pre: atp_pos, atp_q, li_atp.
rims: lamellae       the lamellar cristae of a core geometry: every lamella is a stack of slabs (normal n, spacing S,
                     lumen half width CR_W) clipped CR_GAP inside the boundary inner membrane (b, the matrix SDF of the
                     boundary IM); its rim line is the slab mid-plane (u = 0) at b = -CR_GAP - rho. Rim points are found
                     by projecting seeded random matrix points onto a slab mid-plane and Newton steps on b within that
                     plane; outward = the in-plane gradient of b. rho is tried from params.rho_nm (first site test
                     that passes wins). Sites closer than SPACING to a kept site are skipped. The dimer is built from two
                     copies of the monomer library entry (params.monomer, Jurkat CPX-6151 = 8H9T) with the yeast dimer
                     geometry (params.dimer_angle_deg between the membrane normals, membrane centres on a circle of
                     params.rim_radius_nm; C2 about the bisector z; in-plane monomer angle psi: the clash-free psi
                     (no bead pair < clash_nm) with the most transmembrane bead pairs within contact_nm, then the
                     smallest monomer-monomer gap, then the smallest |psi|; 2-degree scan). NOTE (2026-11-17): the
                     lamellar cristae of em_cell1_v8 have flat caps (lumen 2 CR_W = 20 nm, square corners), not the
                     curved rims of the yeast prims; with the yeast dimer (84 deg, 8.8 nm) no site keeps every TM
                     residue within 2.5 nm (best site median 4.6 nm); tm_tol_nm / rim_radius_nm are params. B6b (2026-11-17):
                     em_cell1_v8r rounds the rims (CR_RIM = 10 nm, rim axis at b = -CR_GAP - 10, so rho ~ 10); the
                     flat TM footprint still exceeds 2.5 nm on any smooth rim (ideal cylinder: >= 2.5-2.7 nm).
                     params.rim_sampling (lamellae only):
                       random (default; B4)  rim points from params.samples random points of the id-grid box (~1 rim
                                             point per um of rim on em_cell1_v8r), tried in random order, skipped within
                                             SPACING of a kept site.
                       rows (J6 fix 2026-11-18)  dense rim lines: params.seeds_per_cell seeds per enabled cristae id cell
                                             are projected onto their rim lines, each rim line is traced every
                                             params.trace_step_nm (lamella_rim_rows) and walked in seeded line order: a
                                             dimer is tried at each trace point, after a placed dimer the walk moves one
                                             row spacing (spacing_nm) on. Poses per point: params.tilt_deg (rotations
                                             about the rim tangent; the rims of oblique lamellae are sheared) x rho_nm,
                                             in that order; the first pose clear of the placed dimers that passes the
                                             site test is kept. The site test is site_ok evaluated by first_ok (same
                                             beads and decision; adaptive bead order, early exits, poses in parallel
                                             chunks: the result does not depend on the thread count).
                     params.obstacles (lamellae): ctx.obs sets the dimers must not overlap (Jurkat: npc).
                     Exported as monomer instances (two per dimer): instance group
                     "atp_synthase" (dimer index), consumption {monomer: 2 n}; obs["atp"] = the dimers' spheres.
"""
import json

import numpy as np

from wholecell.preplaced import kern
from wholecell.registry import register

LIB_ID = "ATP_dimer"
MONOMER = "CPX_3281_Mitochondrial_proton_transportin"
SPACING = 13.0
TM_TOL = 2.5
MAX_IN_MEMBRANE = 0.02
H = 10.0


def _make_site_ok(K):
    import numba as nb
    c3 = K.c3
    NMEM = int(K.nmem); M_IM = K.mem("MITO_IM"); MATRIX = K.lab("mito_matrix"); IMS = K.lab("mito_IMS")
    mmd = K.mmd

    @nb.njit(cache=False)
    def site_ok(X, tm, G, pores, pidx, PT, PQ, PC, tm_tol, max_in):
        """one rim site (yeast atp_dimers._site_ok): transmembrane residues within tm_tol of the inner membrane; every
        third soluble residue in the matrix or the intermembrane space, at most max_in of them inside the bilayer"""
        out = np.zeros(NMEM); n_sol = 0; n_in = 0
        for k in range(X.shape[0]):
            mtx = c3.mem_sdfs(X[k, 0], X[k, 1], X[k, 2], G, pores, pidx, PT, PQ, PC, out)
            if tm[k]:
                if abs(out[M_IM]) > tm_tol:
                    return False
            elif k % 3 == 0:
                l = c3.label_from(X[k, 0], X[k, 1], X[k, 2], G, pores, out, mtx)
                if l != MATRIX and l != IMS:
                    return False
                n_sol += 1
                if mmd(out, -1, G) < 0.6:
                    n_in += 1
        return n_in <= max_in * max(n_sol, 1)

    return site_ok


def _make_site_rows(K):
    import numba as nb
    c3 = K.c3
    NMEM = int(K.nmem); M_IM = K.mem("MITO_IM"); MATRIX = K.lab("mito_matrix"); IMS = K.lab("mito_IMS")
    mmd = K.mmd

    @nb.njit(cache=False)
    def pose_ok(B, tm, sol, order, R, C, in_cap, G, pores, pidx, PT, PQ, PC, tm_tol, max_in):
        """the site test (= site_ok) of the dimer frame B rotated by R at C, beads transformed one at a time in `order`,
        which adapts (a failing bead moves halfway to the front; the decision does not depend on the order)"""
        out = np.zeros(NMEM)
        n_sol = 0; n_in = 0
        for jj in range(order.shape[0]):
            k = order[jj]
            if not (tm[k] or sol[k]):
                continue
            x0 = R[0, 0] * B[k, 0] + R[0, 1] * B[k, 1] + R[0, 2] * B[k, 2] + C[0]
            x1 = R[1, 0] * B[k, 0] + R[1, 1] * B[k, 1] + R[1, 2] * B[k, 2] + C[1]
            x2 = R[2, 0] * B[k, 0] + R[2, 1] * B[k, 1] + R[2, 2] * B[k, 2] + C[2]
            mtx = c3.mem_sdfs(x0, x1, x2, G, pores, pidx, PT, PQ, PC, out)
            if tm[k]:
                bad = abs(out[M_IM]) > tm_tol
            else:
                l = c3.label_from(x0, x1, x2, G, pores, out, mtx)
                bad = l != MATRIX and l != IMS
            if bad:
                h = jj // 2
                order[jj] = order[h]; order[h] = k
                return False
            if not tm[k]:
                n_sol += 1
                if mmd(out, -1, G) < 0.6:
                    n_in += 1
                    if n_in > in_cap:          # the final count can only be larger: same decision, earlier
                        h = jj // 2
                        order[jj] = order[h]; order[h] = k
                        return False
        return n_in <= max_in * max(n_sol, 1)

    @nb.njit(cache=False, parallel=True)
    def first_ok(B, tm, sol, orders, Rs, C, G, pores, pidx, PT, PQ, PC, tm_tol, max_in):
        """first pose j (dimer frame B rotated by Rs[j], at C[j]) that passes the site test; -1 if none. Poses are
        tested in chunks of orders.shape[0] in parallel (slot i keeps its own bead order); the first passing index is
        returned, so the result does not depend on the chunk size (thread count)"""
        n_test = 0
        for k in range(B.shape[0]):
            if sol[k] and not tm[k]:
                n_test += 1
        in_cap = max_in * max(n_test, 1)
        P = orders.shape[0]; nv = Rs.shape[0]
        flags = np.zeros(P, np.bool_)
        for c0 in range(0, nv, P):
            m = min(P, nv - c0)
            for i in nb.prange(m):
                flags[i] = pose_ok(B, tm, sol, orders[i], Rs[c0 + i], C[c0 + i], in_cap, G, pores, pidx, PT, PQ, PC, tm_tol, max_in)
            for i in range(m):
                if flags[i]:
                    return c0 + i
        return -1

    return first_ok


def _make_clear(K):
    import numba as nb
    p3 = K.p3

    @nb.njit(cache=False)
    def clear(stc, str_, Rs, C, glo, h, dims, head, nxt, sid, cen, rad, tol):
        """per pose (sphere tree stc / str_ rotated by Rs[j] at C[j]): True if it overlaps no sphere of the grid"""
        P = np.empty((stc.shape[0], 3)); ok = np.zeros(Rs.shape[0], np.bool_)
        for j in range(Rs.shape[0]):
            p3.transform(stc, Rs[j], C[j], P)
            ok[j] = not p3.g_overlap(glo, h, dims, head, nxt, sid, cen, rad, P, str_, tol)
        return ok

    return clear


def place_prims(ctx, n_max, seed, params, log=print):
    """yeast: rows on the rims of the crista prims (= atp_dimers.place; seed = the resolved integer)"""
    from scipy.spatial.transform import Rotation
    K, cell, lp = ctx.K, ctx.cell, ctx.lp
    p3, c3 = K.p3, K.c3
    lib_id = params.get("lib_id", LIB_ID)
    spacing, tm_tol, max_in = float(params.get("spacing_nm", SPACING)), float(params.get("tm_tol_nm", TM_TOL)), float(params.get("max_in_membrane", MAX_IN_MEMBRANE))
    site_ok = kern.factory("atp_synthase.site_ok", K, _make_site_ok)
    li = int(np.where(lp["ids"] == lib_id)[0][0])
    s0, s1 = lp["st_off"][li], lp["st_off"][li + 1]
    stc, str_ = lp["STC"][s0:s1], lp["STR"][s0:s1]
    z = np.load(ctx.lib_path(lib_id)); B = z["beads"].astype(np.float64); tm = z["tm"].astype(bool)
    G = cell.G; glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0)
    cap = int(n_max * len(str_) * 1.1) + 1000
    grid = p3.grid_new(glo, ghi, H, cap, cap * 14)
    rng = np.random.default_rng(seed)
    out = np.zeros(K.nmem); M = np.empty((3, 3)); P = np.empty((len(stc), 3))
    cr = [i for i, c in enumerate(cell.PC) if c == c3.C_CRISTA]
    rng.shuffle(cr)
    pos, quat, n_try = [], [], 0
    bmax = float(np.linalg.norm(B, axis=1).max())
    for i in cr:
        if len(pos) >= n_max:
            break
        q = cell.PQ[i]; c = q[:3]; ax = q[3:6] / np.linalg.norm(q[3:6]); R = q[6]
        rr = R + bmax + 5.0; a_loc = cell.args(c - rr, c + rr, margin=80.0)
        e1 = np.cross(ax, rng.normal(size=3)); e1 /= np.linalg.norm(e1); e2 = np.cross(ax, e1)
        th0 = rng.uniform(0, 2 * np.pi)
        for th in th0 + np.arange(0, 2 * np.pi - spacing / R, spacing / R):
            if len(pos) >= n_max:
                break
            n_try += 1
            u = np.cos(th) * e1 + np.sin(th) * e2; t = np.cross(ax, u)
            Rm = np.c_[ax, t, u]                      # library x -> crista normal, y -> rim tangent, z -> outward
            if np.linalg.det(Rm) < 0:
                Rm[:, 0] *= -1
            centre = c + R * u
            X = B @ Rm.T + centre
            if not site_ok(X, tm, *a_loc, tm_tol, max_in):
                continue
            qq = Rotation.from_matrix(Rm).as_quat(); qq = np.array([qq[3], qq[0], qq[1], qq[2]])
            p3.quat_mat(qq, M); p3.transform(stc, M, centre, P)
            if p3.g_overlap(glo, H, *grid[:6], P, str_, p3.TOL):
                continue
            p3.g_insert(glo, H, *grid, P, str_)
            pos.append(centre); quat.append(qq)
    log(f"ATP synthase dimers: {len(pos)} on cristae rims ({len(pos) / max(n_try, 1):.0%} of {n_try} rim sites tried; "
        f"{2 * len(pos)} monomers of {2 * n_max})")
    return dict(atp_pos=np.array(pos, np.float32).reshape(-1, 3), atp_q=np.array(quat, np.float32).reshape(-1, 4), li_atp=li)


# ------------------------------------------------------------------ lamellar cristae (core geometry)
def _rot_y(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])


def _rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def dimer_model(L, tm, params, log=print):
    """two monomers (library frame: membrane mid-plane z = 0, F1 at +z) in the yeast dimer geometry. Returns the
    monomer poses (R_i, t_i) in the dimer frame (origin on the rim circle axis, x across the ridge, y along it,
    z = bisector of the membrane normals) and the choice made."""
    from scipy.spatial import cKDTree
    alpha = np.radians(float(params.get("dimer_angle_deg", 84.0))) / 2.0
    clash, contact = float(params.get("clash_nm", 0.5)), float(params.get("contact_nm", 1.5))
    ctm = L[tm].mean(0); ctm[2] = 0.0
    rt0 = float(params.get("rim_radius_nm", 8.8))
    best = None
    for k_r in range(41):                       # rim radius grows by 0.25 nm only if no psi is clash-free
        rt = rt0 + 0.25 * k_r
        for psi in np.radians(np.r_[0:180:2.0, -180:0:2.0]):
            poses = []
            for extra in (0.0, np.pi):          # monomer 1 = C2 rotation (about the bisector z) of monomer 0
                R = _rot_z(extra) @ _rot_y(alpha) @ _rot_z(psi)
                n = R @ np.array([0, 0, 1.0])
                t = rt * n - R @ ctm
                poses.append((R, t))
            X0 = L @ poses[0][0].T + poses[0][1]; X1 = L @ poses[1][0].T + poses[1][1]
            d, _ = cKDTree(X0).query(X1)
            if d.min() < clash:
                continue
            T0 = X0[tm]; T1 = X1[tm]
            nc = int(sum(len(x) for x in cKDTree(T0).query_ball_point(T1, contact)))
            key = (nc, -round(float(d.min()), 2), -abs(np.degrees(psi)))
            if best is None or key > best[0]:
                best = (key, poses, dict(rim_radius_nm=rt, psi_deg=float(np.degrees(psi)), tm_contacts=nc, min_gap_nm=float(d.min())))
        if best is not None:
            break
    if best is None:
        raise RuntimeError("ATP synthase dimer: no clash-free monomer arrangement")
    info = dict(best[2], dimer_angle_deg=float(np.degrees(2 * alpha)))
    log(f"ATP synthase dimer model (from the monomer, yeast geometry): {info}")
    return best[1], info


def _make_rims(K):
    import numba as nb
    c3 = K.c3

    @nb.njit(cache=False)
    def bnd(px, py, pz, G, pores, pidx, PT, PQ, PC, out, a, off):
        """matrix SDF of the boundary inner membrane (negative inside; cristae not included)"""
        c3.mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out)
        return out[a] + off

    return bnd


def _im_row(K, G):
    """(membrane index a, offset) of the K_INNER_CRISTAE row of the core membrane table"""
    c3 = K.c3
    o = int(G[c3.G_TABLES])
    for r in range(int(G[o + 1])):
        q = o + c3.TH + c3.MROW * r
        if int(G[q]) == c3.K_INNER_CRISTAE:
            return int(G[q + 1]), float(G[q + 3])
    raise KeyError("no inner-cristae membrane row in this geometry")


def lamella_rim_sites(ctx, rng, n_samples, rho, log=print):
    """rim points (centre on the rim axis at depth rho behind the cap), lamella normal, outward direction"""
    K, cell = ctx.K, ctx.cell
    c3 = K.c3
    G = cell.G
    o = int(G[c3.G_CRISTAE])
    nl = int(G[o]); lo = G[o + 1:o + 4].copy(); h = float(G[o + 4]); dims = np.array([int(G[o + 5]), int(G[o + 6]), int(G[o + 7])])
    gap = float(G[c3.S_CR_GAP])
    a_im, off = _im_row(K, G)
    bnd = kern.factory("atp_synthase.bnd", K, _make_rims)
    args = cell.args(); out = np.zeros(K.nmem)
    X = lo + rng.random((n_samples, 3)) * dims * h
    lab = cell.label(X)
    X = X[lab == K.lab("mito_matrix")]
    sites = []
    def b_at(p):
        return bnd(p[0], p[1], p[2], *args, out, a_im, off)
    for p in X:
        ix = np.floor((p - lo) / h).astype(int)
        k = int(G[o + 8 + 12 * nl + (ix[0] * dims[1] + ix[1]) * dims[2] + ix[2]])
        if k < 0:
            continue
        q = o + 8 + 12 * k
        if G[q + 11] == 0.0:
            continue
        n = G[q:q + 3]; S = G[q + 4]
        t = n @ p - G[q + 3]
        u = t - S * np.floor(t / S + 0.5)
        x = p - u * n                         # nearest slab mid-plane
        target = -gap - rho
        ok = False
        for _ in range(12):                   # Newton on b within the plane
            bx = b_at(x)
            g = np.array([b_at(x + e) - b_at(x - e) for e in np.eye(3)]) / 2.0
            g -= (g @ n) * n; gn = np.linalg.norm(g)
            if gn < 1e-3:
                break
            if abs(bx - target) < 0.2:
                ok = True; break
            x = x - (bx - target) / (gn * gn) * g
        if not ok:
            continue
        w = g / gn                             # outward (toward the boundary IM)
        sites.append(np.r_[x, n, w, k])
    return np.array(sites).reshape(-1, 10)


def _make_trace(K):
    import numba as nb
    c3 = K.c3

    @nb.njit(cache=False)
    def bval(x0, x1, x2, G, pores, pidx, PT, PQ, PC, out, a, off):
        c3.mem_sdfs(x0, x1, x2, G, pores, pidx, PT, PQ, PC, out)
        return out[a] + off

    @nb.njit(cache=False)
    def project(x, n, target, G, pores, pidx, PT, PQ, PC, out, a, off, w):
        """x (in a slab mid-plane) -> the rim line b = target by Newton steps on b within the plane (as
        lamella_rim_sites); w = outward in-plane unit gradient of b. x, w updated in place."""
        g = np.empty(3)
        for _ in range(12):
            bx = bval(x[0], x[1], x[2], G, pores, pidx, PT, PQ, PC, out, a, off)
            g[0] = (bval(x[0] + 1.0, x[1], x[2], G, pores, pidx, PT, PQ, PC, out, a, off)
                    - bval(x[0] - 1.0, x[1], x[2], G, pores, pidx, PT, PQ, PC, out, a, off)) / 2.0
            g[1] = (bval(x[0], x[1] + 1.0, x[2], G, pores, pidx, PT, PQ, PC, out, a, off)
                    - bval(x[0], x[1] - 1.0, x[2], G, pores, pidx, PT, PQ, PC, out, a, off)) / 2.0
            g[2] = (bval(x[0], x[1], x[2] + 1.0, G, pores, pidx, PT, PQ, PC, out, a, off)
                    - bval(x[0], x[1], x[2] - 1.0, G, pores, pidx, PT, PQ, PC, out, a, off)) / 2.0
            gn_ = g[0] * n[0] + g[1] * n[1] + g[2] * n[2]
            for i in range(3):
                g[i] -= gn_ * n[i]
            gn = np.sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2])
            if gn < 1e-3:
                return False
            if abs(bx - target) < 0.2:
                for i in range(3):
                    w[i] = g[i] / gn
                return True
            for i in range(3):
                x[i] -= (bx - target) / (gn * gn) * g[i]
        return False

    @nb.njit(cache=False)
    def lam_id(x, G, o, nl, lo, h, d0, d1, d2):
        i0 = int(np.floor((x[0] - lo[0]) / h)); i1 = int(np.floor((x[1] - lo[1]) / h)); i2 = int(np.floor((x[2] - lo[2]) / h))
        if i0 < 0 or i1 < 0 or i2 < 0 or i0 >= d0 or i1 >= d1 or i2 >= d2:
            return -1
        k = int(G[o + 8 + 12 * nl + (i0 * d1 + i1) * d2 + i2])
        if k >= 0 and G[o + 8 + 12 * k + 11] == 0.0:
            return -1
        return k

    @nb.njit(cache=False)
    def trace(x0, w0, n, k, sgn, step, max_n, target, G, pores, pidx, PT, PQ, PC, out, a, off, o, nl, lo, h, d0, d1, d2, P, W):
        """walk the rim line of lamella k from x0 (outward w0) in the direction sgn * (n x w), one point every `step` nm
        (each projected back onto the rim line); stops where the rim line ends (projection fails, leaves lamella k,
        jumps) or closes on x0. Returns (points written to P / W, closed)."""
        x = x0.copy(); w = w0.copy(); xn = np.empty(3); wn = np.empty(3)
        for i in range(max_n):
            t0 = sgn * (n[1] * w[2] - n[2] * w[1]); t1 = sgn * (n[2] * w[0] - n[0] * w[2]); t2 = sgn * (n[0] * w[1] - n[1] * w[0])
            xn[0] = x[0] + step * t0; xn[1] = x[1] + step * t1; xn[2] = x[2] + step * t2
            if not project(xn, n, target, G, pores, pidx, PT, PQ, PC, out, a, off, wn):
                return i, False
            if lam_id(xn, G, o, nl, lo, h, d0, d1, d2) != k:
                return i, False
            dx = xn[0] - x[0]; dy = xn[1] - x[1]; dz = xn[2] - x[2]
            if dx * t0 + dy * t1 + dz * t2 < 0.3 * step or dx * dx + dy * dy + dz * dz > 4.0 * step * step:
                return i, False
            if i >= 3:
                ex = xn[0] - x0[0]; ey = xn[1] - x0[1]; ez = xn[2] - x0[2]
                if ex * ex + ey * ey + ez * ez < (0.75 * step) ** 2:
                    return i, True
            for j in range(3):
                P[i, j] = xn[j]; W[i, j] = wn[j]; x[j] = xn[j]; w[j] = wn[j]
        return max_n, False

    return project, lam_id, trace


def lamella_cells(G, c3):
    """enabled id cells of the cristae section (index triples) and the section header"""
    o = int(G[c3.G_CRISTAE]); nl = int(G[o]); lo = G[o + 1:o + 4].copy(); h = float(G[o + 4])
    dims = np.array([int(G[o + 5]), int(G[o + 6]), int(G[o + 7])])
    par = G[o + 8:o + 8 + 12 * nl].reshape(nl, 12)
    ids = G[o + 8 + 12 * nl:o + 8 + 12 * nl + dims.prod()].reshape(dims).astype(np.int64)
    en = (ids >= 0) & (par[np.maximum(ids, 0), 11] > 0)
    return np.argwhere(en), ids, dict(o=o, nl=nl, lo=lo, h=h, dims=dims)


def lamella_rim_rows(ctx, rng, params, cells=None, log=print):
    """dense rim lines (B4 J6 fix, 2026-11-18): seeds_per_cell seeds per enabled cristae id cell, drawn uniformly in
    those cells (`cells`: a subset of them, for checks), each seed projected onto its rim line (as lamella_rim_sites, rho = 0), then the whole rim
    line through it is traced in both directions every `trace_step_nm` nm (rim lines already traced are skipped: a
    seed within 2 steps of a traced point). Returns the polylines (points, outward, lamella normal, closed) and the
    traced rim length (nm)."""
    K, cell = ctx.K, ctx.cell
    c3 = K.c3
    G = cell.G
    allc, ids, hd = lamella_cells(G, c3)
    cells = allc if cells is None else cells
    gap = float(G[c3.S_CR_GAP])
    a_im, off = _im_row(K, G)
    project, lam_id, trace = kern.factory("atp_synthase.trace", K, _make_trace)
    args = cell.args(); out = np.zeros(K.nmem)
    step = float(params.get("trace_step_nm", 2.0))
    n_seeds = int(np.ceil(float(params.get("seeds_per_cell", 1.0)) * len(cells)))
    max_n = int(params.get("trace_max_points", 200000))
    o, nl, lo, h, dims = hd["o"], hd["nl"], hd["lo"], hd["h"], hd["dims"]
    X = lo + (cells[rng.integers(0, len(cells), n_seeds)] + rng.random((n_seeds, 3))) * h
    X = X[cell.label(X) == K.lab("mito_matrix")]
    target = -gap
    Pb = np.empty((max_n, 3)); Wb = np.empty((max_n, 3))
    lines, length = [], 0.0
    vox = set(); vh = 2.0 * step
    def seen(p):
        c = np.floor(p / vh).astype(np.int64)
        for d0 in (-1, 0, 1):
            for d1 in (-1, 0, 1):
                for d2 in (-1, 0, 1):
                    if (c[0] + d0, c[1] + d1, c[2] + d2) in vox:
                        return True
        return False
    for p in X:
        k = lam_id(p, G, o, nl, lo, h, dims[0], dims[1], dims[2])
        if k < 0:
            continue
        q = o + 8 + 12 * k
        n = G[q:q + 3].copy(); S = G[q + 4]
        t = n @ p - G[q + 3]
        x = p - (t - S * np.floor(t / S + 0.5)) * n
        w = np.empty(3)
        if not project(x, n, target, *args, out, a_im, off, w) or lam_id(x, G, o, nl, lo, h, dims[0], dims[1], dims[2]) != k:
            continue
        if seen(x):
            continue
        nf, closed = trace(x, w, n, k, 1.0, step, max_n, target, *args, out, a_im, off, o, nl, lo, h, dims[0], dims[1], dims[2], Pb, Wb)
        F, FW = Pb[:nf].copy(), Wb[:nf].copy()
        if closed:
            B, BW = np.zeros((0, 3)), np.zeros((0, 3))
        else:
            nb_, _ = trace(x, w, n, k, -1.0, step, max_n, target, *args, out, a_im, off, o, nl, lo, h, dims[0], dims[1], dims[2], Pb, Wb)
            B, BW = Pb[:nb_][::-1].copy(), Wb[:nb_][::-1].copy()
        L = np.r_[B, x[None], F]; LW = np.r_[BW, w[None], FW]
        keep = np.array([not seen(y) for y in L])           # parts of rim lines traced from an earlier seed
        for y in L[keep]:
            vox.add(tuple(np.floor(y / vh).astype(np.int64)))
        # split into runs of unseen points
        idx = np.flatnonzero(keep)
        if not len(idx):
            continue
        cuts = np.flatnonzero(np.diff(idx) > 1) + 1
        runs = np.split(idx, cuts)
        for r in runs:
            seg = L[r]
            ln = float(np.linalg.norm(np.diff(seg, axis=0), axis=1).sum()) if len(seg) > 1 else 0.0
            if closed and len(runs) == 1 and len(r) == len(L):
                ln += float(np.linalg.norm(seg[-1] - seg[0]))
            length += ln
            lines.append(dict(P=seg, W=LW[r], n=n, closed=bool(closed and len(runs) == 1 and len(r) == len(L)), length=ln))
    log(f"ATP synthase rim lines: {len(lines)} lines traced from {len(X)} matrix seeds ({n_seeds} drawn in "
        f"{len(cells)} cristae id cells); traced rim length {length / 1e3:.1f} um (step {step} nm)")
    return lines, length


def place_lamellae(ctx, n_max, seed, params, log=print, cells=None):
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    K, cell = ctx.K, ctx.cell
    p3 = K.p3
    rng = np.random.default_rng(seed)
    mono = params.get("monomer", "CPX-6151")
    z = np.load(ctx.lib_path(mono)); L = z["beads"].astype(np.float64)
    tm = np.abs(L[:, 2]) < float(params.get("tm_slab_nm", 2.0))      # TM residues: the membrane slab of the library frame
    poses, dinfo = dimer_model(L, tm, params, log=log)
    B = np.concatenate([L @ R.T + t for R, t in poses]); TM = np.r_[tm, tm]
    stc0, str0 = z["st_c"].astype(np.float64), z["st_r"].astype(np.float64)
    stc = np.concatenate([stc0 @ R.T + t for R, t in poses]); str_ = np.r_[str0, str0]
    spacing, tm_tol, max_in = float(params.get("spacing_nm", SPACING)), float(params.get("tm_tol_nm", TM_TOL)), float(params.get("max_in_membrane", MAX_IN_MEMBRANE))
    site_ok = kern.factory("atp_synthase.site_ok", K, _make_site_ok)
    rhos = [float(r) for r in params.get("rho_nm", [6.5, 4.5, 8.5, 2.5, 10.5])]
    sampling = params.get("rim_sampling", "random")
    tilts = [0.0] if sampling != "rows" else [float(a) for a in params.get("tilt_deg", [0.0])]
    if sampling == "rows":
        # same test, faster: only the tested beads, transmembrane beads farthest from the ridge axis (y) first
        SOL = (~TM) & (np.arange(len(B)) % 3 == 0)
        zl = np.abs(np.r_[L[:, 2], L[:, 2]])                  # distance to the monomer's membrane plane (library frame)
        ltm = np.r_[np.hypot(L[:, 0] - L[tm, 0].mean(), L[:, 1] - L[tm, 1].mean())] if tm.any() else np.zeros(len(L))
        ltm = np.r_[ltm, ltm]
        it, isol = np.flatnonzero(TM), np.flatnonzero(SOL & ~TM)
        ORD = np.r_[it[np.argsort(-ltm[it], kind="stable")], isol[np.argsort(zl[isol], kind="stable")]]
        import numba
        ORDS = np.tile(ORD, (max(1, numba.get_num_threads()), 1))
        RT = np.array([np.eye(3) if a == 0.0 else Rotation.from_rotvec(np.array([0.0, 1.0, 0.0]) * np.radians(a)).as_matrix()
                       for a in tilts])
        RHO = np.tile(np.array(rhos), len(tilts)); TLT = np.repeat(np.array(tilts), len(rhos))
        clear = kern.factory("atp_synthase.clear", K, _make_clear)
        first_ok = kern.factory("atp_synthase.site_rows", K, _make_site_rows)
    args = cell.args()
    M = np.empty((3, 3)); P = np.empty((len(stc), 3))
    pos, quat, kept_pts, n_try, rho_used, tilt_used = [], [], [], 0, [], []

    def try_site_rows(x, n, w):
        """rows mode: poses (tilt, rho) in order; the first pose clear of the placed dimers that passes the site test is
        kept (overlap test first: same pose as testing the site first)"""
        t_ = np.cross(n, w)
        Rm0 = np.c_[n, t_, w]
        if np.linalg.det(Rm0) < 0:
            Rm0[:, 0] *= -1
        Rs = np.repeat(np.einsum("ij,tjk->tik", Rm0, RT), len(rhos), axis=0)   # tilt t: rotation about the rim tangent
        Cs = x - RHO[:, None] * Rs[:, :, 2]               # dimer origin (rim circle axis) rho behind the cap point
        idx = np.flatnonzero(clear(stc, str_, Rs, Cs, glo, H, *grid[:6], p3.TOL))
        if not len(idx):
            return False
        j = first_ok(B, TM, SOL, ORDS, Rs[idx], Cs[idx], *args, tm_tol, max_in)
        if j < 0:
            return False
        j = idx[j]
        p3.transform(stc, Rs[j], Cs[j], P)
        p3.g_insert(glo, H, *grid, P, str_)
        pos.append(Cs[j]); quat.append(Rs[j]); kept_pts.append(x); rho_used.append(RHO[j]); tilt_used.append(TLT[j])
        return True

    def try_site(x, n, w):
        """random mode: dimer at rim point x (lamella normal n, outward w), the first rho of rhos whose site passes"""
        t_ = np.cross(n, w)
        Rm = np.c_[n, t_, w]
        if np.linalg.det(Rm) < 0:
            Rm[:, 0] *= -1
        for rho in rhos:
            centre = x - rho * w                      # dimer origin (rim circle axis) rho behind the cap point
            X = B @ Rm.T + centre
            if not site_ok(X, TM, *args, tm_tol, max_in):
                continue
            qq = Rotation.from_matrix(Rm).as_quat(); qq = np.array([qq[3], qq[0], qq[1], qq[2]])
            p3.quat_mat(qq, M); p3.transform(stc, M, centre, P)
            if p3.g_overlap(glo, H, *grid[:6], P, str_, p3.TOL):
                continue
            p3.g_insert(glo, H, *grid, P, str_)
            pos.append(centre); quat.append(Rm); kept_pts.append(x); rho_used.append(rho)
            return True
        return False

    if sampling == "rows":
        # dense rows (B4 J6 fix 2026-11-18): every traced rim line is walked from one end (closed lines: from their
        # seed point), a dimer is tried at each trace point and, once one is placed, the walk jumps one row spacing
        # (spacing_nm) along the rim before trying again; lines in seeded random order.
        lines, rim_len = lamella_rim_rows(ctx, rng, params, cells=cells, log=log)
        order = rng.permutation(len(lines))
        allP = np.concatenate([ln["P"] for ln in lines]) if lines else np.zeros((1, 3))
        glo = allP.min(0) - 100.0; ghi = allP.max(0) + 100.0
        reach = float(np.linalg.norm(B, axis=1).max() + max(rhos) + float(params.get("trace_step_nm", 2.0)) + 1.0)
        OC, OR = _obstacles(ctx, params.get("obstacles") or [], glo, ghi, near=allP, reach=reach, log=log)
        cap = int(min(n_max, rim_len / spacing + 100) * len(str_) * 1.1) + 1000 + len(OR)
        grid = p3.grid_new(glo, ghi, H, cap, cap * 14)
        if len(OR):
            p3.g_insert(glo, H, *grid, OC, OR)
        for li in order:
            if len(pos) >= n_max:
                break
            ln = lines[li]; Pl, Wl, n = ln["P"], ln["W"], ln["n"]
            arc = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(Pl, axis=0), axis=1))]
            first, nxt = None, 0.0
            for i in range(len(Pl)):
                if len(pos) >= n_max:
                    break
                if arc[i] < nxt:
                    continue
                if ln["closed"] and first is not None and arc[i] > first + ln["length"] - spacing:
                    break
                n_try += 1
                if try_site_rows(Pl[i], n, Wl[i]):
                    first = arc[i] if first is None else first
                    nxt = arc[i] + spacing
        npts = int(sum(len(ln["P"]) for ln in lines))
        log(f"ATP synthase dimers (lamellar cristae, rows): {len(pos)} of {n_try} rim sites tried ({npts} rim points on "
            f"{len(lines)} rim lines, {rim_len / 1e3:.1f} um traced, {len(pos) / max(rim_len / 1e3, 1e-9):.1f} per um); "
            f"{2 * len(pos)} monomers of {2 * n_max}; rho used "
            f"{dict(zip(*np.unique(rho_used, return_counts=True))) if rho_used else {}}; tilt used "
            f"{dict(zip(*np.unique(tilt_used, return_counts=True))) if tilt_used else {}}")
        return _finish(pos, quat, poses, stc, str_, dict(dinfo, rim_sampling="rows", sites_tried=n_try, rim_points=npts,
                                                         rim_lines=len(lines), rim_length_traced_um=round(rim_len / 1e3, 3),
                                                         dimers=len(pos)))
    if sampling != "random":
        raise ValueError(f"core.atp_synthase_dimers: rim_sampling {sampling!r} (random | rows)")
    n_samples = int(params.get("samples", 400000))
    S = lamella_rim_sites(ctx, rng, n_samples, 0.0, log=log)      # rim line points (rho = 0: on the cap)
    order = rng.permutation(len(S)); S = S[order]
    glo = S[:, :3].min(0) - 100.0 if len(S) else np.zeros(3); ghi = S[:, :3].max(0) + 100.0 if len(S) else np.ones(3)
    OC, OR = _obstacles(ctx, params.get("obstacles") or [], glo, ghi, log=log)
    cap = int(min(n_max, len(S)) * len(str_) * 1.1) + 1000 + len(OR)
    grid = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if len(OR):
        p3.g_insert(glo, H, *grid, OC, OR)
    tree, tree_n = None, 0
    for s in S:
        if len(pos) >= n_max:
            break
        x, n, w = s[:3], s[3:6], s[6:9]
        if kept_pts:
            if tree is None or len(kept_pts) - tree_n > 2000:
                tree = cKDTree(np.array(kept_pts)); tree_n = len(kept_pts)
            if tree.query_ball_point(x, spacing * 0.999):
                continue
            if len(kept_pts) > tree_n and (np.linalg.norm(np.asarray(kept_pts[tree_n:]) - x, axis=1) < spacing * 0.999).any():
                continue
        n_try += 1
        try_site(x, n, w)
    log(f"ATP synthase dimers (lamellar cristae): {len(pos)} of {n_try} rim sites tried ({len(S)} rim points from "
        f"{n_samples} samples); {2 * len(pos)} monomers of {2 * n_max}; rho used {dict(zip(*np.unique(rho_used, return_counts=True))) if rho_used else {}}")
    return _finish(pos, quat, poses, stc, str_, dict(dinfo, sites_tried=n_try, rim_points=int(len(S)), dimers=len(pos)))


def _obstacles(ctx, names, glo, ghi, near=None, reach=0.0, log=print):
    """obstacle spheres (ctx.obs sets named in params.obstacles, e.g. npc: the pores of the nuclear envelope next to a
    mitochondrion) inside the dimer grid box and, with `near` (rim points), within reach + their radius of a rim point;
    the dimers must not overlap them (lamellae; J6 fix 2026-11-18)"""
    from scipy.spatial import cKDTree
    C, R = [np.zeros((0, 3))], [np.zeros(0)]
    tree = cKDTree(near) if near is not None and len(near) else None
    for nm in names:
        if nm not in ctx.obs:
            raise KeyError(f"core.atp_synthase_dimers: obstacle set {nm!r} not built before this builder")
        c, r = ctx.obs[nm]
        c, r = np.asarray(c, np.float64), np.asarray(r, np.float64)
        k = np.flatnonzero(((c >= glo - r[:, None]) & (c <= ghi + r[:, None])).all(1))
        if tree is not None and len(k):
            d, _ = tree.query(c[k], distance_upper_bound=reach + float(r[k].max()))
            k = k[d <= reach + r[k]]
        C.append(c[k]); R.append(r[k])
        log(f"ATP synthase dimers: obstacles {nm}: {len(k)} of {len(r)} spheres within reach of the rims")
    return np.ascontiguousarray(np.concatenate(C)), np.ascontiguousarray(np.concatenate(R))


def _finish(pos, quat, poses, stc, str_, info):
    from scipy.spatial.transform import Rotation
    # monomer instances
    mp, mq, did = [], [], []
    for d, (c, Rm) in enumerate(zip(pos, quat)):
        for R, t in poses:
            Rw = Rm @ R
            q = Rotation.from_matrix(Rw).as_quat(); mq.append([q[3], q[0], q[1], q[2]])
            mp.append(c + Rm @ t); did.append(d)
    dq = [Rotation.from_matrix(Rm).as_quat() for Rm in quat]
    dq = np.array([[q[3], q[0], q[1], q[2]] for q in dq], np.float32).reshape(-1, 4)
    return dict(pos=np.array(mp, np.float32).reshape(-1, 3), quat=np.array(mq, np.float32).reshape(-1, 4),
                dimer=np.array(did, np.int32), dimer_pos=np.array(pos, np.float32).reshape(-1, 3), dimer_q=dq,
                stc=stc, str_=str_, info=info)


@register("preplaced", "core.atp_synthase_dimers")
def build(ctx, params):
    exp3 = ctx.exp3
    rims = params.get("rims", "prims")
    if rims == "prims":                      # yeast
        mono = params.get("monomer", MONOMER)
        im = (exp3.ing_id == mono) & (exp3.region == params.get("region", "mem_MITO_IM"))
        dim = place_prims(ctx, int(exp3.loc[im, "expected"].sum() // 2), ctx.seed(params.get("seed", 7)), params, log=ctx.log)
        ctx.pre.update(dim)
        ctx.consume([mono], 2 * len(dim["atp_pos"]))
        return
    if rims != "lamellae":
        raise ValueError(f"core.atp_synthase_dimers: rims {rims!r} (prims | lamellae)")
    from wholecell.preplaced.context import add_group
    mono = params.get("monomer", "CPX-6151")
    rows = exp3.ing_id == mono
    if params.get("region"):
        rows &= exp3.region == params["region"]
    n_pool = float(exp3.loc[rows, "expected"].sum())
    n_max = int(n_pool * float(params.get("max_dimer_share", 1.0)) // 2)
    r = place_lamellae(ctx, n_max, ctx.seed(params.get("seed", "preplaced.atp_synthase")), params, log=ctx.log)
    li = int(np.where(ctx.lp["ids"] == mono)[0][0])
    n = len(r["pos"])
    add_group(ctx, params.get("group", "atp_synthase"), np.full(n, li, np.int32), r["pos"], r["quat"], dimer=r["dimer"],
              consumed={mono: n})
    ctx.pre.update(atp_pos=r["dimer_pos"], atp_q=r["dimer_q"])
    # obstacles: the dimers' spheres (dimer tree in the dimer frame)
    from scipy.spatial.transform import Rotation
    if len(r["dimer_pos"]):
        q = r["dimer_q"].astype(np.float64)
        Rm = Rotation.from_quat(np.c_[q[:, 1:], q[:, :1]]).as_matrix()
        C = (np.einsum("nij,sj->nsi", Rm, r["stc"]) + r["dimer_pos"][:, None, :]).reshape(-1, 3)
        ctx.obs["atp"] = (C, np.tile(r["str_"], len(Rm)))
    else:
        ctx.obs["atp"] = (np.zeros((0, 3)), np.zeros(0))
    ctx.consume([mono], n)
    ctx.state["atp_synthase"] = dict(r["info"], pool=n_pool, monomers=n)
