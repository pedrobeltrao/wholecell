# origin: yeast pipeline/microtubules.py:place @ e3b7448 (procedural track source: copied, geometry through ctx.K)
# origin: yeast scripts/s3_pack.py:main (closure microtubules(), first block; tubulin and gamma-TuRC pool lines after
#         the chromatin-bound proteins) @ e3b7448 (builders core.microtubules_mtoc / core.microtubule_pools)
# origin: jurkat/em/cytoskel/mt_stage3.py:build_lattice / transport_frames / lattice / quats (= jcytoskeleton.lattice,
#         legacy workspace, 2026-11-10) (traced track source)
# Changed (unification B4, 2026-11-17): one core module; numerics no longer called through ctx.mod. The yeast library
# entries and lattice parameters (yeast microtubules.model, a library builder) are still made by the yeast context and
# read here from the library (ctx.lib_path); the yeast constants (ingredient ids, ring genes, region codes, lengths)
# are params with the yeast values as defaults, so config/cells/yeast.yaml is unchanged. The yeast seed int is the one
# the vendored place() computed (ctx.seed(params.seed) = config.seed(k)).
"""Microtubules (serial pre-placed builders `core.microtubules`, `core.microtubules_mtoc`, `core.microtubule_pools`).

One MT builder, two track sources (params.tracks):

  procedural  nucleation from an MTOC (yeast: the spindle pole body in the nuclear envelope): kinetochore, nuclear and
              cytoplasmic MTs grown layer by layer from gamma-TuRC rings on the SPB plaques until their target length,
              a membrane or an obstacle (yeast pipeline/microtubules.py). Adds pre: mt_*; obs["nucleus"], obs["mt_all"],
              obs["mt_cyto"]; state["mtd"], state["chrom_starts"] (builder core.microtubules_mtoc, the yeast name).
  traced      tracks are polylines of a geometry layer (Jurkat: `mt_traces`, the EM-traced centerlines, stage-3 variant
              b): each polyline is oriented by its stored polarity flag (minus end = the end nearer the NE), resampled
              at the dimer repeat, given parallel-transported frames (one random perpendicular per MT, stream `seed`)
              and filled with the 13-protofilament lattice of the library's lattice json (13_3, seam). Dimers whose
              sphere tree leaves the allowed labels / touches a membrane (check_volume) or overlaps an obstacle set
              (params.obstacles: ctx.obs names) are dropped one by one (as v7 dropped dimers at the centrioles and
              NPCs); MT-MT contacts are not tested (as v7). Records the instance group "microtubules" (mt_id) with
              consumption {params.ingredient: dimers}, obs["mt_all"].

core.microtubule_pools: deducts the placed dimers and ring subunits (yeast order: after the chromatin-bound proteins).
"""
import json
import time

import numpy as np

from wholecell.preplaced import kern
from wholecell.registry import register

# yeast defaults (pipeline/microtubules.py constants)
LIB_MT, LIB_RING_NUC, LIB_RING_CYTO = "MT_tubulin", "gTuRC_nuc", "gTuRC_cyto"
NPF = 13
N_KMT, KMT_LEN = 16, (250.0, 60.0)
N_NUC, NUC_LEN = 4, (400.0, 100.0)
N_CYTO, CYTO_LEN = 3, (1200.0, 300.0)
MIN_LEN = 60.0
TRIES = 400
POLYMER_CAP = 0.7
PLAQUE_R = 90.0
PLAQUE_GAP = 25.0
CONE, CONE_CYTO = np.radians(30.0), np.radians(70.0)
BEND = 8.2 / 300.0
LOOK, SAFE = 40.0, 60.0
H = 20.0
TUBULIN_ING = "CPX_1424_Tubulin_alpha_beta_heterodimeric"
TUSC_ING = "CPX_1198_Gamma_tubulin_small_complex"
RING_GENES = {"TUB4": 13, "SPC97": 7, "SPC98": 6}
R_CYTO, R_NUCP = 0, 1          # yeast tiles3d region codes of mt_reg


def _lattice_json(ctx, lib_id):
    return json.load(open(ctx.lib_path(lib_id)[:-4] + "_lattice.json"))


def _rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def _quat(Rm):
    from scipy.spatial.transform import Rotation
    q = Rotation.from_matrix(Rm).as_quat()
    return np.array([q[3], q[0], q[1], q[2]])


def _frame(w, rng):
    e1 = np.cross(w, rng.normal(size=3)); e1 /= np.linalg.norm(e1)
    return np.c_[e1, np.cross(w, e1), w]


# ------------------------------------------------------------------ procedural (yeast SPB)
def place_procedural(ctx, n_dimers, obstacles=None, seed=0, params=None, log=print):
    """MTs + rings from the SPB. n_dimers: alpha/beta-tubulin dimer pool; seed: the integer seed (already resolved).
    Returns instances (mt_pos/q/type/reg), kMT tips. = yeast microtubules.place with the geometry through ctx.K."""
    P_ = params or {}
    K, cell, lp = ctx.K, ctx.cell, ctx.lp
    p3, c3 = K.p3, K.c3
    lib_mt, lib_rn, lib_rc = P_.get("lib_mt", LIB_MT), P_.get("lib_ring_nuc", LIB_RING_NUC), P_.get("lib_ring_cyto", LIB_RING_CYTO)
    prm = _lattice_json(ctx, lib_mt)
    li = {k: int(np.where(lp["ids"] == k)[0][0]) for k in (lib_mt, lib_rn, lib_rc)}
    phi, dz, D, r0 = np.radians(prm["phi_deg"]), prm["dz"], prm["dimer"], prm["r0"]
    ring_c = np.array(prm["ring_c"])
    rng = np.random.default_rng(seed)
    G = cell.G; cn = G[3:6]
    nrm = G[7:10] - cn; nrm /= np.linalg.norm(nrm)
    spb_dir = -nrm                                   # SPB opposite the nucleolus (as in the chromatin pre-pass)
    _, n_out = c3.nuc_surface(G, spb_dir)            # NE normal at the SPB
    r_nucp, r_cyto = int(P_.get("region_nucleus", R_NUCP)), int(P_.get("region_cytosol", R_CYTO))
    plaques = {"nuc": (c3.nuc_surface(G, spb_dir, -(c3.PERI_W + PLAQUE_GAP))[0], -n_out, 1 << K.lab("nucleoplasm"), r_nucp, lib_rn),
               "cyto": (c3.nuc_surface(G, spb_dir, PLAQUE_GAP)[0], n_out, 1 << K.lab("cytosol"), r_cyto, lib_rc)}
    plan = [("nuc", "kMT", rng.normal(*KMT_LEN)) for _ in range(N_KMT)] + \
           [("nuc", "nMT", rng.normal(*NUC_LEN)) for _ in range(N_NUC)] + [("cyto", "cMT", rng.normal(*CYTO_LEN)) for _ in range(N_CYTO)]
    need = sum(max(L, MIN_LEN) / D * NPF for _, _, L in plan)
    scale = min(1.0, POLYMER_CAP * n_dimers / need)
    glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0)
    cap = int(need * 1.2) * 40 + 100000 + (len(obstacles[1]) if obstacles is not None else 0)
    grid = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        OC, OR = obstacles
        k = np.all((OC > glo) & (OC < ghi), axis=1)
        p3.g_insert(glo, H, *grid, OC[k], OR[k])
    sub = {}
    for k, v in li.items():
        s0, s1 = lp["st_off"][v], lp["st_off"][v + 1]
        sub[k] = (lp["STC"][s0:s1], lp["STR"][s0:s1], float(lp["BR"][v]))
    out = np.zeros(K.nmem); M = np.empty((3, 3))

    def body(lid, centre, Rm, mask):
        stc, str_, br = sub[lid]
        P = np.empty((len(stc), 3)); p3.quat_mat(_quat(Rm), M); p3.transform(stc, M, centre, P)
        if p3.g_overlap(glo, H, *grid[:6], P, str_, p3.TOL):
            return None
        if not p3.check_volume(P, str_, centre[0], centre[1], centre[2], br, mask, *a, out):
            return None
        return P

    def clearance(p):
        c3.mem_sdfs(p[0], p[1], p[2], *a, out)
        return K.mmd(out, -1, G)

    def steer(p, w):
        """bend toward more membrane clearance (at most BEND per layer) when a membrane is within SAFE ahead"""
        q = p + w * LOOK
        if clearance(q) > SAFE:
            return w
        g = np.array([clearance(q + e * 4.0) - clearance(q - e * 4.0) for e in np.eye(3)])
        g -= (g @ w) * w; n = np.linalg.norm(g)
        if n < 1e-6:
            return w
        w2 = np.cos(BEND) * w + np.sin(BEND) * g / n
        return w2 / np.linalg.norm(w2)

    def grow(O, w, ring, mask, nlay):
        W = _frame(w, rng)
        Pr = body(ring, O + W @ ring_c, W, mask)
        if Pr is None:
            return None, O, W
        got = [(ring, O + W @ ring_c, W, Pr)]; Oj = O.copy()
        for j in range(nlay):
            layer = []
            for k in range(NPF):
                Rk = W @ _rz(k * phi)
                cen = Oj + Rk @ np.array([r0, 0.0, 0.0]) + W[:, 2] * (k * dz)
                P = body(lib_mt, cen, Rk, mask)
                if P is None:
                    return got, Oj, W
                layer.append((lib_mt, cen, Rk, P))
            got += layer
            Oj = Oj + W[:, 2] * D
            w2 = steer(Oj, W[:, 2])                  # parallel-transport the frame to the new direction
            ax = np.cross(W[:, 2], w2); sn = np.linalg.norm(ax)
            if sn > 1e-9:
                from scipy.spatial.transform import Rotation
                W = Rotation.from_rotvec(ax / sn * np.arcsin(min(sn, 1.0))).as_matrix() @ W
        return got, Oj, W

    pos, rot, typ, reg, mt_id = [], [], [], [], []
    tips, tdir, kinds, lengths = [], [], [], []
    for side, kind, L in plan:
        L = max(MIN_LEN, L * scale); nlay = int(round(L / D))
        centre0, normal, mask, region, ring = plaques[side]
        Wn = _frame(normal, rng)
        best = None
        for _ in range(TRIES):                       # random minus-end site; MTs splay out from the plaque (bouquet)
            r = PLAQUE_R * np.sqrt(rng.random()); th = rng.uniform(0, 2 * np.pi)
            radial = np.cos(th) * Wn[:, 0] + np.sin(th) * Wn[:, 1]
            O = centre0 + r * radial + normal * rng.uniform(0, 15.0)
            cone = CONE if side == "nuc" else CONE_CYTO
            w = normal + np.tan(cone) * (r / PLAQUE_R) * radial + 0.15 * rng.normal(size=3); w /= np.linalg.norm(w)
            got, Oe, We = grow(O, w, ring, mask, nlay)
            if got is None:
                continue
            if best is None or len(got) > len(best[0]):
                best = (got, Oe, We)
            if len(got) - 1 >= nlay * NPF:
                break
        if best is None or (len(best[0]) - 1) / NPF * D < MIN_LEN:
            continue
        got, Oe, We = best
        for lid, cen, Rm, P in got:
            p3.g_insert(glo, H, *grid, P, sub[lid][1])
            pos.append(cen); rot.append(_quat(Rm)); typ.append(li[lid]); reg.append(region); mt_id.append(len(lengths))
        n_l = (len(got) - 1) // NPF
        tips.append(Oe + We[:, 2] * 12.0); tdir.append(We[:, 2].copy()); kinds.append(kind); lengths.append(n_l * D)
    # chromosomes without a kinetochore MT start next to the inner plaque
    while sum(k == "kMT" for k in kinds) < N_KMT:
        tips.append(plaques["nuc"][0] + plaques["nuc"][1] * 40.0); tdir.append(np.asarray(plaques["nuc"][1], float)); kinds.append("kMT"); lengths.append(np.nan)
    lengths = np.array(lengths); kinds = np.array(kinds)
    kmt_tips = np.array([t for t, k in zip(tips, kinds) if k == "kMT"], np.float64).reshape(-1, 3)
    kmt_dirs = np.array([t for t, k in zip(tdir, kinds) if k == "kMT"], np.float64).reshape(-1, 3)
    keep = ~np.isnan(lengths); lengths, kinds = lengths[keep], kinds[keep]
    n_tub = int(sum(t == li[lib_mt] for t in typ))
    log(f"microtubules: {len(lengths)} MTs ({', '.join(f'{k} {int((kinds == k).sum())}' for k in ('kMT', 'nMT', 'cMT'))}), "
        f"{n_tub} tubulin dimers ({n_tub / max(n_dimers, 1):.0%} of pool, length scale {scale:.2f}); "
        f"mean length kMT {lengths[kinds == 'kMT'].mean() if (kinds == 'kMT').any() else 0:.0f} nm, "
        f"cMT {lengths[kinds == 'cMT'].mean() if (kinds == 'cMT').any() else 0:.0f} nm; lattice {prm}")
    return dict(mt_pos=np.array(pos, np.float32), mt_q=np.array(rot, np.float32), mt_type=np.array(typ, np.int32),
                mt_reg=np.array(reg, np.int16), mt_id=np.array(mt_id, np.int32), mt_len=lengths.astype(np.float32),
                mt_kind=kinds.astype("U3"), kmt_tips=kmt_tips, kmt_dirs=kmt_dirs,
                n_rings_nuc=int(sum(t == li[lib_rn] for t in typ)), n_rings_cyto=int(sum(t == li[lib_rc] for t in typ)),
                n_tub=n_tub)


@register("preplaced", "core.microtubules_mtoc")
def build(ctx, params):
    """procedural MTs from the MTOC (yeast stage-3 block; the name of the yeast config)"""
    exp3 = ctx.exp3
    tub_ing = params.get("tubulin_ing", TUBULIN_ING)
    tub = (exp3.ing_id == tub_ing) & (exp3.region == params.get("pool_region", "cytosol"))
    mtd = {}
    mtd.update(place_procedural(ctx, float(exp3.loc[tub, "expected"].sum()), obstacles=ctx.obs["npc"],
                                seed=ctx.seed(params.get("seed", 31)), params=params, log=ctx.log))
    nucl = mtd["mt_reg"] == int(params.get("region_nucleus", R_NUCP))
    ctx.obs["nucleus"] = ctx.emit(mtd["mt_type"][nucl], mtd["mt_pos"][nucl], mtd["mt_q"][nucl])
    ctx.obs["mt_all"] = ctx.emit(mtd["mt_type"], mtd["mt_pos"], mtd["mt_q"])
    cyto = mtd["mt_reg"] == int(params.get("region_cytosol", R_CYTO))
    ctx.obs["mt_cyto"] = ctx.emit(mtd["mt_type"][cyto], mtd["mt_pos"][cyto], mtd["mt_q"][cyto])
    ctx.state["mtd"] = mtd
    ctx.state["chrom_starts"] = mtd["kmt_tips"]
    spindle_to_pre(ctx)


def spindle_to_pre(ctx):
    """the spindle arrays of state['mtd'] (prefixes mt, kt, kto, spb) into pre, as stage 3 copied them"""
    ctx.pre.update({k: v for k, v in ctx.state["mtd"].items() if k.split("_")[0] in ("mt", "kt", "kto", "spb")})


@register("preplaced", "core.microtubule_pools")
def pools(ctx, params):
    exp3 = ctx.exp3; mtd = ctx.state["mtd"]
    tub = (exp3.ing_id == params.get("tubulin_ing", TUBULIN_ING)) & (exp3.region == params.get("pool_region", "cytosol"))
    exp3.loc[tub, "expected"] = (exp3.loc[tub, "expected"] - mtd["n_tub"]).clip(lower=0)
    # ring subunits: gamma-TuSC complexes first, then the free subunit pools (SPC110 from its CMD1 heterodimer)
    tusc = int(params.get("tusc_per_ring", 6)); tusc_mult = params.get("tusc_subunits", {"TUB4": 2})
    n_all = mtd["n_rings_nuc"] + mtd["n_rings_cyto"]
    m = ctx.consume([params.get("tusc_ing", TUSC_ING)], tusc * n_all)
    for g, k in (params.get("ring_genes") or RING_GENES).items():
        ctx.consume([g], k * n_all - tusc_mult.get(g, 1) * m)
    att = params.get("ring_attachment", {"ings": ["SPC110", "CMD1_SPC110_het"], "per_nuclear_ring": 14})
    ctx.consume(list(att["ings"]), att["per_nuclear_ring"] * mtd["n_rings_nuc"])


# ------------------------------------------------------------------ traced (EM polylines)
def _resample(P, step):
    """arc-length resample to ~step spacing, endpoints kept (mt_stage3.resample)"""
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]
    if d[-1] < 1e-6:
        return P[:1]
    n = max(int(round(d[-1] / step)), 1); t = np.linspace(0, d[-1], n + 1)
    return np.stack([np.interp(t, d, P[:, k]) for k in range(3)], 1)


def _transport_frames(A, rng):
    t = np.gradient(A, axis=0) if len(A) > 1 else np.array([[0, 0, 1.0]]); t /= np.linalg.norm(t, axis=1, keepdims=True) + 1e-12
    x = np.cross(t[0], rng.normal(size=3)); x /= np.linalg.norm(x)
    F = np.zeros((len(A), 3, 3))
    for i in range(len(A)):
        x = x - (x @ t[i]) * t[i]; x /= np.linalg.norm(x)
        F[i] = np.c_[x, np.cross(t[i], x), t[i]]
    return F


def _lattice(axis_pts, frames, prm):
    from scipy.spatial.transform import Rotation
    phi, dz, r0 = np.radians(prm["phi_deg"]), prm["dz"], prm["r0"]
    n = len(axis_pts); ang = np.arange(NPF) * phi
    c, s = np.cos(ang), np.sin(ang)
    Rz = np.zeros((NPF, 3, 3)); Rz[:, 0, 0] = c; Rz[:, 0, 1] = -s; Rz[:, 1, 0] = s; Rz[:, 1, 1] = c; Rz[:, 2, 2] = 1
    Rk = np.einsum("nij,kjl->nkil", frames, Rz).reshape(-1, 3, 3)
    cen = (axis_pts[:, None, :] + np.einsum("nkij,j->nki", Rk.reshape(n, NPF, 3, 3), np.array([r0, 0, 0]))
           + frames[:, None, :, 2] * (np.arange(NPF) * dz)[None, :, None]).reshape(-1, 3)
    q = Rotation.from_matrix(Rk).as_quat(); q = np.c_[q[:, 3], q[:, 0], q[:, 1], q[:, 2]]; q[q[:, 0] < 0] *= -1
    return cen, q


def _make_traced(K):
    import numba as nb
    p3 = K.p3
    TOL = float(p3.TOL)
    K_NMEM = int(K.nmem)

    @nb.njit(cache=False)
    def screen(pos, quat, stc, str_, br, mask, G, pores, pidx, PT, PQ, PC, glo, h, grids_n, dims, heads, nxts, sids, cens, rads):
        """per dimer: 0 kept, 1 geometry (label / membrane), 2 + j overlap with obstacle set j"""
        n = pos.shape[0]; code = np.zeros(n, np.int8)
        out = np.zeros(K_NMEM); M = np.empty((3, 3)); P = np.empty((stc.shape[0], 3)); q = np.empty(4); t = np.empty(3)
        for i in range(n):
            for k in range(4):
                q[k] = quat[i, k]
            for k in range(3):
                t[k] = pos[i, k]
            p3.quat_mat(q, M); p3.transform(stc, M, t, P)
            if not p3.check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out):
                code[i] = 1
                continue
            for j in range(grids_n):
                if p3.g_overlap(glo, h, dims[j], heads[j], nxts[j], sids[j], cens[j], rads[j], P, str_, TOL):
                    code[i] = 2 + j
                    break
        return code

    return screen


def place_traced(ctx, params, log=print):
    from numba.typed import List
    K, cell = ctx.K, ctx.cell
    p3 = K.p3
    lib_mt = params.get("lib", "MT_TUBA1B-TUBB4B")
    prm = _lattice_json(ctx, lib_mt)
    layer = ctx.geom.layer(params.get("layer", "mt_traces"))
    z = np.load(layer, allow_pickle=True)
    pts, off = z["points_model_nm"].astype(np.float64), z["offsets"]
    flip = z["flipped_for_lattice"] if "flipped_for_lattice" in z.files else np.zeros(len(off) - 1, bool)
    rng = np.random.default_rng(ctx.seed(params.get("seed", "preplaced.microtubules")))
    t0 = time.time()
    C, Q, mid, Ltot = [], [], [], 0.0
    for i in range(len(off) - 1):
        P = pts[off[i]:off[i + 1]]
        Ltot += float(np.linalg.norm(np.diff(P, axis=0), axis=1).sum()) if len(P) > 1 else 0.0
        if flip[i]:
            P = P[::-1]
        A = _resample(P, prm["dimer"])
        c_, q_ = _lattice(A, _transport_frames(A, rng), prm)
        C.append(c_); Q.append(q_); mid.append(np.full(len(c_), i, np.int32))
    C = np.vstack(C); Q = np.vstack(Q); mid = np.concatenate(mid)
    n0 = len(C)
    # collision screen: geometry and obstacle sets
    li = int(np.where(ctx.lp["ids"] == lib_mt)[0][0])
    s0, s1 = ctx.lp["st_off"][li], ctx.lp["st_off"][li + 1]
    stc, str_ = np.ascontiguousarray(ctx.lp["STC"][s0:s1]), np.ascontiguousarray(ctx.lp["STR"][s0:s1]); br = float(ctx.lp["BR"][li])
    glo = C.min(0) - 100.0; ghi = C.max(0) + 100.0; h = float(params.get("grid_nm", 30.0))
    names = [o for o in (params.get("obstacles") or []) if o in ctx.obs]
    dims, heads, nxts, sids, cens, rads = List(), List(), List(), List(), List(), List()
    for o in names:
        OC, OR = ctx.obs[o]
        OC = np.asarray(OC, np.float64); OR = np.asarray(OR, np.float64)
        k = np.all((OC > glo - OR[:, None]) & (OC < ghi + OR[:, None]), axis=1)
        OC, OR = np.ascontiguousarray(OC[k]), np.ascontiguousarray(OR[k])
        cap = len(OR) + 10
        g = p3.grid_new(glo, ghi, h, cap, cap * 27 + 10)
        p3.g_insert(glo, h, *g, OC, OR)
        d_, hd, nx, sd, ce, ra, _ = g
        dims.append(d_); heads.append(hd); nxts.append(nx); sids.append(sd); cens.append(ce); rads.append(ra)
    if not names:          # numba needs typed lists with one element type
        g = p3.grid_new(glo, glo + 1.0, h, 1, 1)
        dims.append(g[0]); heads.append(g[1]); nxts.append(g[2]); sids.append(g[3]); cens.append(g[4]); rads.append(g[5])
    mask = 0
    for lab in params.get("labels", ["cytosol"]):
        mask |= 1 << K.lab(lab)
    screen = kern.factory("microtubules.traced", K, _make_traced)
    code = screen(C, Q, stc, str_, br, mask, *cell.args(), glo, h, len(names), dims, heads, nxts, sids, cens, rads)
    keep = code == 0
    drops = {"geometry": int((code == 1).sum())}
    for j, o in enumerate(names):
        drops[f"obstacle:{o}"] = int((code == 2 + j).sum())
    C, Q, mid = C[keep], Q[keep], mid[keep]
    log(f"microtubules (traced): {len(off) - 1} tracks of {layer.split('/')[-1]}, {Ltot / 1e3:.1f} um; lattice "
        f"{n0:,} dimers, {len(C):,} kept; dropped {drops} ({time.time() - t0:.0f} s)")
    return dict(pos=C.astype(np.float32), quat=Q.astype(np.float32), mt_id=mid, n_lattice=n0, drops=drops,
                n_tracks=int(len(off) - 1), length_um=round(Ltot / 1e3, 1), li=li)


@register("preplaced", "core.microtubules")
def build_tracks(ctx, params):
    """one MT builder: params.tracks = procedural (MTOC; = core.microtubules_mtoc) | traced (EM polylines)"""
    tracks = params.get("tracks", "procedural")
    if tracks == "procedural":
        return build(ctx, params)
    if tracks != "traced":
        raise ValueError(f"core.microtubules: tracks {tracks!r} (procedural | traced)")
    from wholecell.preplaced.context import add_group
    r = place_traced(ctx, params, log=ctx.log)
    n = len(r["pos"])
    types = np.full(n, r["li"], np.int32)
    ing = params.get("ingredient", "TUBA1B-TUBB4B")
    add_group(ctx, params.get("group", "microtubules"), types, r["pos"], r["quat"], mt_id=r["mt_id"], kind=np.ones(n, np.int8),
              consumed={ing: n})
    ctx.pre.update(mt_pos=r["pos"], mt_q=r["quat"], mt_id=r["mt_id"])
    if params.get("emit_obstacles", True):
        if params.get("obstacle_shape", "tree") == "bounding":     # one bounding sphere per dimer (memory: ~70 M tree spheres)
            ctx.obs["mt_all"] = (np.asarray(r["pos"], np.float64), np.full(n, float(ctx.lp["BR"][r["li"]])))
        else:
            ctx.obs["mt_all"] = ctx.emit(types, r["pos"], r["quat"])
    ctx.state["microtubules"] = dict(dimers=n, lattice_dimers=r["n_lattice"], drops=r["drops"], tracks=r["n_tracks"],
                                     length_um=r["length_um"])
