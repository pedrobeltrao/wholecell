# origin: yeast pipeline/mtdna.py:place @ e3b7448 (copied, geometry through ctx.K)
# origin: yeast pipeline/pack3d.py:arm_walks @ e3b7448 (copied into this module, geometry through ctx.K)
# origin: yeast scripts/s3_pack.py:main ("mtDNA nucleoids in the matrix" block) @ e3b7448
# Changed (unification B4, 2026-11-17): one core module, no ctx.mod; the yeast constants are params with the yeast
# values as defaults (genomes per fL 0.45 x cell volume 34 fL = inventory.CELL_VOL; config/cells/yeast.yaml unchanged);
# nucleoid sites: params.sites = geometry (yeast: the geometry's mito_nucleoids / mito_segments) | matrix (core
# geometry: seeded random points in mito_matrix, >= min_gap_nm apart; walks per nucleoid in a local box).
"""mtDNA nucleoids in the mitochondrial matrix (serial pre-placed builder `core.mtdna_nucleoids`).

Each nucleoid is a cluster of LOOP_KB DNA loops: self-avoiding persistent walks (10 bp per 3.4 nm step) that start near
the nucleoid centre and are drawn back to it (arm_walks), drawn at 2 bp per bead. Genomes: one per nucleoid plus the
rest spread at random over the nucleoids. Avoids the ATP synthase dimers (yeast: pre atp_*; core: obs["atp"]).
Adds pre: mtdna_pts, mtdna_off, mtdna_centres, mtdna_genomes.
"""
import numpy as np

from wholecell.preplaced import kern
from wholecell.registry import register

N_NUCLEOIDS, GENOME_KB = 15, 85.8
GENOMES_PER_FL = 0.45
CELL_VOL_FL = 34.0            # yeast inventory.CELL_VOL (um^3)
LOOP_KB = 5.0
STEP, BP_PER_STEP = 3.4, 10
DNA_R = 1.0
DNA_STEP = 0.68


def _make_walks(K):
    import numba as nb
    p3, c3 = K.p3, K.c3
    H_GRID = float(p3.H_GRID); TOL = float(p3.TOL); NMEM = int(K.nmem)
    grid_new, g_insert, g_overlap = p3.grid_new, p3.g_insert, p3.g_overlap
    rand_quat, quat_mat, transform, check_volume = p3.rand_quat, p3.quat_mat, p3.transform, p3.check_volume

    @nb.njit(cache=False)
    def arm_walks(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask,
                  OC, OR, seg_len, seg_start, seg_target, r_start, step, cap_s, cap_e):
        """= yeast pack3d.arm_walks (persistent self-avoiding walks from seg_start towards seg_target)"""
        np.random.seed(seed)
        dims, head, nxt, sid, cen, rad, cnt = grid_new(glo, ghi, H_GRID, cap_s, cap_e)
        g_insert(glo, H_GRID, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
        ntot = 0
        for c in range(seg_len.shape[0]):
            ntot += seg_len[c]
        pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32)
        chain = np.empty(ntot, np.int32); brk = np.zeros(ntot, np.int8)
        out = np.zeros(NMEM); M = np.empty((3, 3)); q = np.empty(4); t = np.empty(3)
        P = np.empty((stc.shape[0], 3))
        n = 0
        for c in range(seg_len.shape[0]):
            placed = 0; have_prev = False
            dx, dy, dz = 0.0, 0.0, 0.0
            px, py, pz = seg_start[c, 0], seg_start[c, 1], seg_start[c, 2]
            fails = 0; restarts = 0
            while placed < seg_len[c] and restarts < 20000:
                if not have_prev:
                    rr = r_start if placed == 0 else 30.0 + 2.0 * min(restarts, 200)
                    t[0] = px + (np.random.random() * 2 - 1) * rr
                    t[1] = py + (np.random.random() * 2 - 1) * rr
                    t[2] = pz + (np.random.random() * 2 - 1) * rr
                else:
                    gx = seg_target[c, 0] - px; gy = seg_target[c, 1] - py; gz = seg_target[c, 2] - pz
                    d = np.sqrt(gx * gx + gy * gy + gz * gz) + 1e-9
                    rem = seg_len[c] - placed
                    w = d / (2.0 * step * np.sqrt(rem + 1.0)) - 0.25
                    w = min(max(w, 0.0), 1.5)
                    ex = np.random.normal(); ey = np.random.normal(); ez = np.random.normal()
                    en = np.sqrt(ex * ex + ey * ey + ez * ez)
                    ex = 0.5 * dx + ex / en + w * gx / d; ey = 0.5 * dy + ey / en + w * gy / d; ez = 0.5 * dz + ez / en + w * gz / d
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

    return arm_walks


def _beads(pos, chain, brk):
    """2-bp beads along each unbroken piece of each loop"""
    pts, off = [], [0]
    cut = np.r_[np.where((np.diff(chain) != 0) | (brk[1:] == 1))[0] + 1, len(pos)]
    a0 = 0
    for b0 in cut:
        P = pos[a0:b0].astype(np.float64); a0 = b0
        if len(P) < 2:
            continue
        seg = np.linalg.norm(np.diff(P, axis=0), axis=1); s = np.r_[0, np.cumsum(seg)]
        t = np.arange(0, s[-1], DNA_STEP)
        Q = np.c_[[np.interp(t, s, P[:, j]) for j in range(3)]].T
        pts.append(Q); off.append(off[-1] + len(Q))
    return pts, off


def place(ctx, obstacles=None, seed=0, params=None, log=print):
    """= yeast mtdna.place with the geometry through ctx.K (sites: the geometry's)"""
    P_ = params or {}
    K, cell = ctx.K, ctx.cell
    MATRIX = K.lab("mito_matrix")
    arm_walks = kern.factory("mtdna.arm_walks", K, _make_walks)
    rng = np.random.default_rng(seed)
    segs = np.array(cell.d["mito_segments"], np.float64)
    seg_len = np.linalg.norm(segs[:, 1] - segs[:, 0], axis=1)
    sites = cell.d.get("mito_nucleoids")
    if sites:
        centres = [np.array(s_[0]) for s_ in sites if cell.label(np.array(s_[0]))[0] == MATRIX]
    else:
        centres = []
        for _ in range(20000):
            if len(centres) >= N_NUCLEOIDS:
                break
            k = rng.choice(len(segs), p=seg_len / seg_len.sum())
            p = segs[k, 0] + (segs[k, 1] - segs[k, 0]) * rng.uniform(0.1, 0.9) + rng.normal(0, 25.0, 3)
            if cell.label(p)[0] != MATRIX or any(np.linalg.norm(p - c) < 500.0 for c in centres):
                continue
            centres.append(p)
    centres = np.array(centres)
    GENOMES = max(len(centres), int(round(float(P_.get("genomes_per_fl", GENOMES_PER_FL)) * float(P_.get("cell_vol_fl", CELL_VOL_FL)))))
    genomes = np.bincount(rng.integers(0, len(centres), GENOMES - len(centres)), minlength=len(centres)) + 1
    loops = [(i, int(round(g * GENOME_KB / LOOP_KB))) for i, g in enumerate(genomes)]
    steps = int(LOOP_KB * 1000 / BP_PER_STEP)
    starts = np.array([centres[i] for i, n in loops for _ in range(n)])
    L = np.full(len(starts), steps, np.int64)
    G = cell.G; glo = centres.min(0) - 600; ghi = centres.max(0) + 600
    a = cell.args(glo, ghi, margin=100.0)
    OC, OR = obstacles if obstacles is not None else (np.zeros((0, 3)), np.zeros(0))
    k = np.all((OC > glo) & (OC < ghi), axis=1); OC, OR = OC[k].astype(np.float64), OR[k].astype(np.float64)
    nsph = int(L.sum()) + len(OR)
    pos, _, chain, brk = arm_walks(seed, glo, ghi, *a, np.zeros((1, 3)), np.array([DNA_R]), DNA_R, 1 << MATRIX,
                                   OC, OR, L, starts, starts, 40.0, STEP, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    pts, off = _beads(pos, chain, brk)
    pts = np.concatenate(pts).astype(np.float32)
    rad = [np.percentile(np.linalg.norm(pos[np.isin(chain, np.where(np.repeat(np.arange(len(loops)), [n for _, n in loops]) == i)[0])] - centres[i], axis=1), 90)
           for i in range(len(centres))]
    log(f"mtDNA: {len(centres)} nucleoids, {GENOMES} genomes ({genomes.min()}-{genomes.max()} per nucleoid), "
        f"{len(pos) * BP_PER_STEP / 1e6:.2f} of {L.sum() * BP_PER_STEP / 1e6:.2f} Mb placed in {len(off) - 1} pieces; "
        f"90% of each nucleoid's DNA within {np.median(rad):.0f} nm (median) of its centre")
    return dict(mtdna_pts=pts, mtdna_off=np.array(off, np.int64), mtdna_centres=centres.astype(np.float32),
                mtdna_genomes=genomes.astype(np.int32))


def place_matrix(ctx, obstacles, seed, params, log=print):
    """core geometry: nucleoid centres = seeded random mito_matrix points >= min_gap_nm apart; genomes = params.genomes
    (1 per nucleoid + the rest at random); the walks of each nucleoid in a local box (centre +- box_nm)"""
    from scipy.spatial import cKDTree
    K, cell = ctx.K, ctx.cell
    MATRIX = K.lab("mito_matrix")
    arm_walks = kern.factory("mtdna.arm_walks", K, _make_walks)
    rng = np.random.default_rng(seed)
    genome_kb = float(params.get("genome_kb", 16.569))
    n_gen = int(params["genomes"]); per = float(params.get("genomes_per_nucleoid", 1.4))
    n_target = max(1, int(round(n_gen / per)))
    gap = float(params.get("min_gap_nm", 500.0)); box = float(params.get("box_nm", 600.0))
    # candidate points: uniform in the bounding box of the matrix (label test), in batches
    lo, hi = params.get("lo"), params.get("hi")
    if lo is None:
        R = float(cell.G[0]); lo, hi = np.full(3, -R), np.full(3, R)
    lo, hi = np.asarray(lo, float), np.asarray(hi, float)
    centres = []
    tries = 0
    while len(centres) < n_target and tries < int(params.get("max_samples", 20_000_000)):
        X = lo + rng.random((200000, 3)) * (hi - lo); tries += len(X)
        X = X[cell.label(X) == MATRIX]
        for p in X:
            if len(centres) >= n_target:
                break
            if centres and (np.linalg.norm(np.asarray(centres) - p, axis=1) < gap).any():
                continue
            centres.append(p)
    centres = np.array(centres).reshape(-1, 3)
    if not len(centres):
        raise RuntimeError("mtDNA: no matrix point found")
    n_gen = max(n_gen, len(centres))
    genomes = np.bincount(rng.integers(0, len(centres), n_gen - len(centres)), minlength=len(centres)) + 1
    steps = int(LOOP_KB * 1000 / BP_PER_STEP)
    OC, OR = obstacles if obstacles is not None else (np.zeros((0, 3)), np.zeros(0))
    OC, OR = np.asarray(OC, np.float64), np.asarray(OR, np.float64)
    ot = cKDTree(OC) if len(OC) else None
    seeds = rng.integers(0, 2 ** 31 - 1, len(centres))
    args = cell.args()
    pts_all, off_all, placed_steps, total_steps, rad = [], [0], 0, 0, []
    for i, c in enumerate(centres):
        n_loops = int(round(genomes[i] * genome_kb / LOOP_KB))
        n_loops = max(n_loops, 1)
        starts = np.repeat(c[None], n_loops, 0)
        L = np.full(n_loops, steps, np.int64)
        glo, ghi = c - box, c + box
        if ot is not None:
            k = ot.query_ball_point(c, box * 1.75)
            oc, orr = np.ascontiguousarray(OC[k]), np.ascontiguousarray(OR[k])
        else:
            oc, orr = np.zeros((0, 3)), np.zeros(0)
        nsph = int(L.sum()) + len(orr)
        pos, _, chain, brk = arm_walks(int(seeds[i]), glo, ghi, *args, np.zeros((1, 3)), np.array([DNA_R]), DNA_R, 1 << MATRIX,
                                       oc, orr, L, starts, starts, 40.0, STEP, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
        placed_steps += len(pos); total_steps += int(L.sum())
        pts, off = _beads(pos, chain, brk)
        for p in pts:
            pts_all.append(p); off_all.append(off_all[-1] + len(p))
        if len(pos):
            rad.append(np.percentile(np.linalg.norm(pos - c, axis=1), 90))
    pts = np.concatenate(pts_all).astype(np.float32) if pts_all else np.zeros((0, 3), np.float32)
    log(f"mtDNA (matrix sites): {len(centres)} nucleoids (target {n_target}, >= {gap:.0f} nm apart), {n_gen} genomes of "
        f"{genome_kb} kb ({genomes.min()}-{genomes.max()} per nucleoid), {placed_steps * BP_PER_STEP / 1e6:.2f} of "
        f"{total_steps * BP_PER_STEP / 1e6:.2f} Mb placed in {len(off_all) - 1} pieces; 90% of each nucleoid's DNA within "
        f"{np.median(rad) if rad else 0:.0f} nm (median) of its centre")
    return dict(mtdna_pts=pts, mtdna_off=np.array(off_all, np.int64), mtdna_centres=centres.astype(np.float32),
                mtdna_genomes=genomes.astype(np.int32)), dict(nucleoids=int(len(centres)), target=n_target, genomes=n_gen,
                                                               mb_placed=placed_steps * BP_PER_STEP / 1e6,
                                                               mb_target=total_steps * BP_PER_STEP / 1e6)


@register("preplaced", "core.mtdna_nucleoids")
def build(ctx, params):
    sites = params.get("sites", "geometry")
    if sites == "geometry":                 # yeast
        pre = ctx.pre
        AC, AR = ctx.emit(np.full(len(pre["atp_pos"]), pre["li_atp"], np.int32), pre["atp_pos"], pre["atp_q"])
        pre.update(place(ctx, obstacles=(AC, AR), seed=ctx.seed(params.get("seed", 81)), params=params, log=ctx.log))
        return
    if sites != "matrix":
        raise ValueError(f"core.mtdna_nucleoids: sites {sites!r} (geometry | matrix)")
    obs = [ctx.obs[o] for o in (params.get("obstacles") or ["atp"]) if o in ctx.obs]
    O = (np.concatenate([o[0] for o in obs]), np.concatenate([o[1] for o in obs])) if obs else None
    r, info = place_matrix(ctx, O, ctx.seed(params.get("seed", "preplaced.mtdna")), params, log=ctx.log)
    ctx.pre.update(r)
    ctx.state["mtdna"] = info
