# origin: yeast pipeline/polysome.py (walk, _resample, _hermite, place, _qmul, _mat2quat), pipeline/mrna_genes.py
#         (transcripts, _transcripts_inferred) and scripts/s3_pack.py:main ("polysomes: a share of the cytosolic
#         ribosomes" block) @ e3b7448
# Changed (B4, 2026-11-17): one core module; the geometry tests of the walk come from the cell's kernels (ctx.K; the
# yeast kernels are its X5 evaluator, so the yeast output is bitwise unchanged); the mRNA channel, disome pose,
# ribosome library entry, start radius, constants of the mRNA model and the cell-specific data (gene table, tethers,
# exclusion zones) are params / registered providers; NEW option: EM point priors (every retained EM ribosome is placed
# at its measured position, +- a seeded nudge, and starts or continues a chain; Jurkat).
# Changed (walk fix, 2026-11-18): NEW option params.walk (candidate search for the grown ribosomes, priors first);
# default = the legacy walk above (yeast bitwise unchanged).
"""Polysomes (serial pre-placed builder `core.polysomes`).

A share of the cytosolic ribosomes threaded on mRNAs. Each translated mRNA is one chain with its own ribosome count
and linker lengths (free nucleotides between neighbours ~ exponential, linker = end-to-end distance of a single-
stranded RNA of that contour); each next ribosome is placed so that its mRNA entry lies the linker beyond the previous
exit (mRNA channel of the 80S library entry), in a cone around the exit direction, clear of everything before it and
inside the allowed labels; neighbours with <= collide_nt free nucleotides are placed as a collided disome. Chains can
be tethered to a membrane (first ribosome near a seed point, all within tmax of the membrane).

Transcripts (params.transcripts):
  provider: <registered preplaced_data name>   per-gene table (gene, orf, copies, cds, utr5, utr3, polya, te) of the
            cell (yeast: Pelechano/SGD/TIF-seq/PAL-seq tables); the mRNA set and ribosome loads follow the measured
            density and occupancy (mode inferred: as many mRNAs as carry the polysome budget) or a count (mode count)
  table: per-mRNA rows (gene, length) of a frozen transcript set (Jurkat: the v7 mRNA layer); the ribosome load per
            translated mRNA is Poisson(length x density) with the density solved so that the ribosomes in chains of >= 2
            equal the polysome budget (polysome_fraction x the RIBO80S pool)
Point priors (params.point_priors; Jurkat): EM ribosome centres (parquet; filters; voxel size; frame centre from the
  geometry), duplicates closer than dedup_nm merged, grouped by single linkage at link_nm (the v7 polysome-nucleus
  rule). Groups take the largest chains (rank matched; a group larger than every chain is split), EM singletons seed the
  remaining chains, then the rest are placed as monosomes (chains of 1, not on an mRNA). Every prior is placed at its
  position (first try) or within nudge_nm of it (up to prior_tries tries) with a random orientation, else dropped; ER-
  class groups are tethered to the ER for their grown ribosomes.

Linker length (params.linker_nt_mean; B4 addendum 2026-11-17): the mean number of free nucleotides between
  neighbouring ribosomes that the walk draws (exponential) per link. Absent / null = per chain from the transcript
  (cds / ribosomes - footprint_nt; yeast: the literature value of its gene table, unchanged). A number = the same mean
  for every chain (a per-cell value measured from EM ribosome data; Jurkat: config/cells/jurkat.yaml, provenance there
  and in docs/unification_criteria/B4_preplaced_addendum_2026-11-17.md). Only the 3D walk uses it: the ribosome loads
  and the mRNA lengths are unchanged.

Walk search (params.walk; walk fix 2026-11-18): {candidates, rounds, pick, dir_cone, orient_cone, priors_first}.
  candidates 0 (default) = the legacy walk (one pose per try, entry direction exactly against the link, the link's
  free nt redrawn after K_ORIENT failed orientations, the chain ends after max_fails; yeast: bitwise unchanged).
  candidates K > 0: every grown ribosome after the first (not a collided disome, not a prior) is chosen among rounds
  of K candidate poses on the shell the link can reach (entry = previous exit + L D; L the end-to-end draw of the
  round's free nt; D within dir_cone of the exit direction; entry direction within orient_cone of -D; random roll);
  each round draws the free nt anew; of the candidates that fit (no overlap, allowed labels, tether) of the first round
  with any fit, `pick` places: nearest (centre nearest the previous ribosome: contact spacing), compact (nearest the
  centroid of the chain so far) or first; the chain ends only when rounds + 1 rounds give no fit. priors_first (with
  point priors): the priors of all chains are placed before any chain grows. Jurkat: config/cells/jurkat.yaml,
  docs/unification_criteria/B4_preplaced_addendum_polysomes_2026-11-18.md.

Adds pre: poly_*, li_ribo; state: transcripts (tx, tx_chains, tx_free), poly, n_mrna, channel. Core context: instance
group "polysomes" (chain, origin 0 = EM prior / 1 = grown, er) with its RIBO80S consumption.
"""
import os

import numpy as np
import numba as nb

from wholecell.registry import register
from wholecell.preplaced import kern
from ._trees import sphere_tree

# defaults = the yeast values (pipeline/polysome.py, mrna_genes.py; their sources are in those files)
POLY_FRACTION = 0.835
MEAN_SIZE = 5.0
MRNA_SHARE = 0.85
MAX_SIZE = 20
LINK_MIN, LINK_MAX = 3.0, 12.0
CONE = 0.7
SS_LP = 1.5
SS_NM_PER_NT = 0.6
FOOTPRINT_NT = 28
RANDOM_SPACING = True
MAX_FAILS_FLEX = 1200
K_ORIENT = 40
CONE_FLEX = 1.2
H_WALK = 20.0
NT_STEP = 0.5
COLLIDE_NT = 5.0
N_MRNA = 3.0e4               # wc3d.N_MRNA (generic mRNAs only)
# walk fix 2026-11-18 (params.walk): candidate search for the grown ribosomes; candidates 0 = the legacy walk (yeast,
# bitwise unchanged). dir_cone None = the legacy cone; orient_cone: entry direction within this angle (rad) of -D;
# priors_first: with point priors, every chain's priors are placed before any chain grows (grown ribosomes of earlier
# chains no longer take the positions of later EM priors)
WALK_SEARCH = dict(candidates=0, rounds=0, pick="nearest", dir_cone=None, orient_cone=np.pi, priors_first=False)


@nb.njit(cache=True)
def _qmul(a, b, out):
    out[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3]
    out[1] = a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2]
    out[2] = a[0] * b[2] - a[1] * b[3] + a[2] * b[0] + a[3] * b[1]
    out[3] = a[0] * b[3] + a[1] * b[2] - a[2] * b[1] + a[3] * b[0]


def _make(K):
    p3, c3, NMEM = K.p3, K.c3, K.nmem

    @nb.njit
    def walk(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, stc, str_, br, mask, sizes, ex, tex, en, q_en, OC, OR, cap_s, cap_e,
             lrms, cone, max_fails, seeds, tmem, tmax, q_rel, t_rel, collide_nt, R0, npri, poff, PP, nudge, pri_tries,
             wk, wrounds, wpick, wcone, wocone, pfirst):
        """Grow one chain per entry of sizes (see polysome.py walk). New (B4): npri[c] > 0: the first npri[c] ribosomes of
        chain c are priors at PP[poff[c] + k] (first try exactly there, then up to pri_tries tries within nudge nm,
        random orientation; a prior that does not fit is dropped and the chain grows instead); the tether test applies
        to grown ribosomes only. npri all 0 = the yeast walk (same random draws). New (walk fix 2026-11-18): wk > 0 =
        candidate search for every grown (not collided, not prior) ribosome after the first: rounds of wk candidate poses
        on the shell reachable through the link (entry = exit + L D, L the end-to-end draw of the round's free nt, D within
        wcone of the exit direction, entry direction within wocone of -D, random roll); among the candidates that fit,
        the one nearest the previous centre (wpick 0), the first (wpick 1) or the one nearest the centroid of the chain
        so far (wpick 2, compact) is placed; the chain ends only when wrounds
        rounds give no fit. wk = 0 = the legacy walk (no extra random draws). pfirst: the priors of all chains are
        placed first (phase 0, chain order), then every chain grows from its last prior (phase 1); the caller re-sorts
        the output by chain. Returns pos, quat, chain index, collided flag,
        prior flag, dropped priors."""
        np.random.seed(seed)
        dims, head, nxt, sid, cen, rad, cnt = p3.grid_new(glo, ghi, H_WALK, cap_s, cap_e)
        if OC.shape[0]:
            p3.g_insert(glo, H_WALK, dims, head, nxt, sid, cen, rad, cnt, OC, OR)     # e.g. actin cables
        ntot = sizes.sum()
        pos = np.empty((ntot, 3), np.float32); qu = np.empty((ntot, 4), np.float32); chain = np.empty(ntot, np.int32)
        coll = np.zeros(ntot, np.bool_); pri = np.zeros(ntot, np.bool_); qprev = np.empty(4); tprev = np.empty(3); Mp = np.empty((3, 3))
        out = np.zeros(NMEM); M = np.empty((3, 3)); q = np.empty(4); q2 = np.empty(4); t = np.empty(3)
        P = np.empty((stc.shape[0], 3)); D = np.empty(3); negD = np.empty(3)
        exw = np.empty(3); texw = np.empty(3)
        E = np.empty(3); tb = np.empty(3); qb = np.empty(4)
        n = 0; ndrop = 0
        nch = sizes.shape[0]
        nph = 2 if pfirst else 1
        c_pl = np.zeros(nch, np.int64); c_t = np.zeros((nch, 3)); c_q = np.zeros((nch, 4))
        c_sum = np.zeros((nch, 3) if wpick == 2 else (1, 3))     # pick compact: sum of the chain's centres
        for it in range(nph * nch):
            ph = it // nch; c = it - ph * nch       # pfirst: phase 0 = the priors of every chain, phase 1 = growth
            placed = 0; fails = 0; tm = tmem[c]; cur_L = -1.0; tries_L = 0; x_pend = -1.0; coll_tried = False; pk = 0
            if pfirst and ph == 1:
                placed = c_pl[c]; pk = npri[c]
                if placed > 0:                       # continue from the chain's last prior
                    for k in range(4):
                        q[k] = c_q[c, k]; qprev[k] = q[k]
                    p3.quat_mat(q, M)
                    for k in range(3):
                        tprev[k] = c_t[c, k]
                        exw[k] = c_t[c, k] + M[k, 0] * ex[0] + M[k, 1] * ex[1] + M[k, 2] * ex[2]
                        texw[k] = M[k, 0] * tex[0] + M[k, 1] * tex[1] + M[k, 2] * tex[2]
            while placed < sizes[c]:
                if pfirst and ph == 0 and pk >= npri[c]:
                    break
                coll_now = False
                prior = pk < npri[c]
                if prior:
                    kk = poff[c] + pk
                    if fails == 0:
                        for k in range(3):
                            t[k] = PP[kk, k]
                    else:
                        while True:
                            for k in range(3):
                                D[k] = np.random.random() * 2 - 1
                            if D[0] * D[0] + D[1] * D[1] + D[2] * D[2] <= 1.0:
                                break
                        for k in range(3):
                            t[k] = PP[kk, k] + nudge * D[k]
                    p3.rand_quat(q)
                elif placed == 0 and tm >= 0:
                    if fails > 300:
                        tm = -1
                    for k in range(3):
                        t[k] = seeds[c, k] + (np.random.random() * 2 - 1) * 15.0
                    p3.rand_quat(q)
                elif placed == 0:
                    while True:
                        for k in range(3):
                            t[k] = (np.random.random() * 2 - 1) * R0
                        if t[0] * t[0] + t[1] * t[1] + t[2] * t[2] < R0 * R0:
                            break
                    p3.rand_quat(q)
                elif collide_nt > 0.0 and lrms[c] > 0.0 and RANDOM_SPACING and cur_L < 0.0 and not coll_tried:
                    # draw this link's free nucleotides; a collision places the trailing ribosome in the disome pose
                    nt = lrms[c] ** 2 / (2.0 * SS_LP * SS_NM_PER_NT)
                    x_pend = -np.log(1.0 - np.random.random()) * nt
                    coll_tried = True
                    if x_pend <= collide_nt:
                        coll_now = True
                        _qmul(qprev, q_rel, q); p3.quat_mat(qprev, Mp)
                        for k in range(3):
                            t[k] = tprev[k] + Mp[k, 0] * t_rel[0] + Mp[k, 1] * t_rel[1] + Mp[k, 2] * t_rel[2]
                search = wk > 0 and placed > 0 and not coll_now and not prior
                if search:
                    # candidate search on the reachable shell of this link (free nt fixed per link, x_pend if drawn)
                    nt = lrms[c] ** 2 / (2.0 * SS_LP * SS_NM_PER_NT) if lrms[c] > 0.0 else 0.0
                    best = 1e30; found = False
                    for rr in range(wrounds + 1):
                        # every round draws the link's free nt anew (round 0: the pending draw of the collision test),
                        # as the legacy walk does after K_ORIENT failed orientations
                        xd = x_pend if x_pend >= 0.0 else (-np.log(1.0 - np.random.random()) * nt if nt > 0.0 else 0.0)
                        x_pend = -1.0
                        rms = np.sqrt(2.0 * SS_LP * SS_NM_PER_NT * max(xd, 5.0)) if nt > 0.0 else 0.0
                        for jj in range(wk):
                            while True:
                                for k in range(3):
                                    D[k] = np.random.normal()
                                dn = np.sqrt(D[0] ** 2 + D[1] ** 2 + D[2] ** 2)
                                s = 0.0
                                for k in range(3):
                                    D[k] = D[k] / dn
                                    s += D[k] * texw[k]
                                if s > np.cos(wcone):
                                    break
                            while True:     # entry direction within wocone of -D
                                for k in range(3):
                                    E[k] = np.random.normal()
                                dn = np.sqrt(E[0] ** 2 + E[1] ** 2 + E[2] ** 2)
                                s = 0.0
                                for k in range(3):
                                    E[k] = E[k] / dn
                                    s -= E[k] * D[k]
                                if s > np.cos(wocone):
                                    break
                            if rms > 0.0:
                                sg = rms / np.sqrt(3.0)
                                a1 = np.random.normal() * sg; a2 = np.random.normal() * sg; a3 = np.random.normal() * sg
                                L = min(max(np.sqrt(a1 * a1 + a2 * a2 + a3 * a3), 2.0), 4.0 * rms)
                            else:
                                L = LINK_MIN + np.random.random() * (LINK_MAX - LINK_MIN)
                            p3.quat_z_to(E, np.random.random() * 2 * np.pi, q2)
                            _qmul(q2, q_en, q)
                            p3.quat_mat(q, M)
                            for k in range(3):
                                t[k] = exw[k] + L * D[k] - (M[k, 0] * en[0] + M[k, 1] * en[1] + M[k, 2] * en[2])
                            sc = 0.0
                            if wpick == 2:          # compact: nearest the centroid of the chain's ribosomes so far
                                for k in range(3):
                                    sc += (t[k] - c_sum[c, k] / placed) ** 2
                            else:
                                for k in range(3):
                                    sc += (t[k] - tprev[k]) ** 2
                            if found and (wpick == 1 or sc >= best):     # pick 0 / 2: only a better candidate is tested
                                continue
                            p3.transform(stc, M, t, P)
                            ok = not p3.g_overlap(glo, H_WALK, dims, head, nxt, sid, cen, rad, P, str_, p3.TOL)
                            if ok:
                                ok = p3.check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out)
                            if ok and tm >= 0:
                                c3.mem_sdfs(t[0], t[1], t[2], G, pores, pidx, PT, PQ, PC, out)
                                ok = abs(out[tm]) <= tmax
                            if ok:
                                found = True; best = sc
                                for k in range(3):
                                    tb[k] = t[k]
                                for k in range(4):
                                    qb[k] = q[k]
                        if found:
                            break
                    if not found:
                        break            # no candidate fits in any round: the polysome ends shorter
                    for k in range(3):
                        t[k] = tb[k]
                    for k in range(4):
                        q[k] = qb[k]
                if placed > 0 and not coll_now and not prior and not search:
                    while True:   # direction within CONE of the previous exit direction
                        for k in range(3):
                            D[k] = np.random.normal()
                        dn = np.sqrt(D[0] ** 2 + D[1] ** 2 + D[2] ** 2)
                        s = 0.0
                        for k in range(3):
                            D[k] = D[k] / dn
                            s += D[k] * texw[k]
                        if s > np.cos(cone):
                            break
                    if cur_L > 0.0 and tries_L < K_ORIENT:
                        L = cur_L                    # the linker is set by the mRNA: try other orientations first
                    elif lrms[c] > 0.0:
                        rms = lrms[c]
                        if RANDOM_SPACING:
                            nt = lrms[c] ** 2 / (2.0 * SS_LP * SS_NM_PER_NT)
                            xd = x_pend if x_pend >= 0.0 else -np.log(1.0 - np.random.random()) * nt
                            x_pend = -1.0
                            rms = np.sqrt(2.0 * SS_LP * SS_NM_PER_NT * max(xd, 5.0))
                        sg = rms / np.sqrt(3.0)
                        a1 = np.random.normal() * sg; a2 = np.random.normal() * sg; a3 = np.random.normal() * sg
                        L = min(max(np.sqrt(a1 * a1 + a2 * a2 + a3 * a3), 2.0), 4.0 * rms)
                    else:
                        L = LINK_MIN + np.random.random() * (LINK_MAX - LINK_MIN)
                    if L != cur_L:
                        cur_L = L; tries_L = 0
                    tries_L += 1
                    for k in range(3):
                        negD[k] = -D[k]
                    p3.quat_z_to(negD, np.random.random() * 2 * np.pi, q2)   # outward entry direction -> -D
                    _qmul(q2, q_en, q)
                    p3.quat_mat(q, M)
                    for k in range(3):
                        t[k] = exw[k] + L * D[k] - (M[k, 0] * en[0] + M[k, 1] * en[1] + M[k, 2] * en[2])
                p3.quat_mat(q, M); p3.transform(stc, M, t, P)
                good = search or not p3.g_overlap(glo, H_WALK, dims, head, nxt, sid, cen, rad, P, str_, p3.TOL)
                if good and not search:
                    good = p3.check_volume(P, str_, t[0], t[1], t[2], br, mask, G, pores, pidx, PT, PQ, PC, out)
                if good and tm >= 0 and not prior and not search:
                    c3.mem_sdfs(t[0], t[1], t[2], G, pores, pidx, PT, PQ, PC, out)
                    good = abs(out[tm]) <= tmax
                if not good:
                    fails += 1
                    if prior and fails > pri_tries:
                        pk += 1; fails = 0; ndrop += 1        # the prior does not fit: dropped
                        continue
                    if placed > 0 and fails > max_fails:
                        break            # polysome ends shorter
                    continue
                if not p3.g_insert(glo, H_WALK, dims, head, nxt, sid, cen, rad, cnt, P, str_):
                    return pos[:n], qu[:n], chain[:n], coll[:n], pri[:n], ndrop
                for k in range(3):
                    exw[k] = t[k] + M[k, 0] * ex[0] + M[k, 1] * ex[1] + M[k, 2] * ex[2]
                    texw[k] = M[k, 0] * tex[0] + M[k, 1] * tex[1] + M[k, 2] * tex[2]
                    pos[n, k] = t[k]
                for k in range(4):
                    qu[n, k] = q[k]; qprev[k] = q[k]
                for k in range(3):
                    tprev[k] = t[k]
                chain[n] = c; coll[n] = coll_now; pri[n] = prior
                if wpick == 2:
                    for k in range(3):
                        c_sum[c, k] += t[k]
                if prior:
                    pk += 1
                cur_L = -1.0; x_pend = -1.0; coll_tried = False
                n += 1; placed += 1; fails = 0
            if pfirst and ph == 0:
                c_pl[c] = placed
                for k in range(3):
                    c_t[c, k] = tprev[k]
                for k in range(4):
                    c_q[c, k] = qprev[k]
        return pos[:n], qu[:n], chain[:n], coll[:n], pri[:n], ndrop

    return walk


def _resample(P, step=NT_STEP):
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1)
    s = np.r_[0, np.cumsum(seg)]
    u = np.linspace(0, s[-1], max(2, int(np.ceil(s[-1] / step)) + 1))
    return np.c_[[np.interp(u, s, P[:, k]) for k in range(3)]].T


def _hermite(a, ta, b, tb):
    L = np.linalg.norm(b - a)
    s = np.linspace(0, 1, max(3, int(L * 2)))[:, None]
    h00 = 2 * s ** 3 - 3 * s ** 2 + 1; h10 = s ** 3 - 2 * s ** 2 + s; h01 = -2 * s ** 3 + 3 * s ** 2; h11 = s ** 3 - s ** 2
    return h00 * a + h10 * ta * L + h01 * b + h11 * tb * L


def _mat2quat(R):
    from scipy.spatial.transform import Rotation
    x, y, z, w = Rotation.from_matrix(R).as_quat()
    return np.array([w, x, y, z])          # pack3d convention (w, x, y, z)


def rot_to_z(v):
    # origin: yeast pipeline/build_structs.py:rot_to_z @ e3b7448 (unchanged)
    v = v / np.linalg.norm(v)
    z = np.array([0, 0, 1.0])
    c = np.dot(v, z)
    if c > 0.999999:
        return np.eye(3)
    if c < -0.999999:
        return np.diag([1, -1, -1.0])
    ax = np.cross(v, z); s = np.linalg.norm(ax); ax /= s
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def place(K, cell, lp, ribo, beads, ch, disome, n_ribo, seed, obstacles=None, log=print, n_mrna=None, sizes=None, linker_nt=None,
          tether=None, mask=None, R0=None, priors=None, collide_nt=COLLIDE_NT, nudge=10.0, prior_tries=60, walk_per=550,
          mrna_step=NT_STEP, walk_search=None):
    """origin: polysome.place (yeast; see its docstring). B4: the ribosome library id, its beads, the channel ch and the
    disome pose (q, t) are arguments; seed is the resolved integer; priors = (npri per chain, positions in chain order)."""
    from scipy.spatial.transform import Rotation
    li = int(np.where(lp["ids"] == ribo)[0][0])
    stc, str_ = sphere_tree(beads.astype(np.float64), per=walk_per)  # yeast: 550 residues per sphere, ~30 coarse spheres
    ch = np.asarray(ch, np.float64)
    ex, en = ch[0], ch[-1]
    tex = ch[0] - ch[3]; tex /= np.linalg.norm(tex)
    ten = ch[-1] - ch[-4]; ten /= np.linalg.norm(ten)
    q_en = _mat2quat(rot_to_z(ten))
    rng = np.random.default_rng(seed)
    if sizes is not None:
        sizes = np.asarray(sizes, np.int64); mean = float(sizes.mean())
    else:
        mean = MEAN_SIZE if not n_mrna else max(MEAN_SIZE, n_ribo / (MRNA_SHARE * n_mrna))
        if mean > MEAN_SIZE:
            log(f"polysomes: mean size raised to {mean:.1f} ribosomes so that {n_ribo:.0f} ribosomes fit on {MRNA_SHARE:.0%} of {n_mrna} mRNAs")
        sizes = []
        while sum(sizes) < n_ribo:
            sizes.append(int(min(MAX_SIZE, 1 + rng.poisson(mean - 1))))
        sizes = np.array(sizes, np.int64)
    G = cell.G
    glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0)
    OC, OR = obstacles if obstacles is not None else (np.zeros((0, 3)), np.zeros(0))
    lrms = np.zeros(len(sizes)) if linker_nt is None else np.sqrt(2.0 * SS_LP * np.maximum(np.asarray(linker_nt, float), 5.0) * SS_NM_PER_NT)
    ns = len(str_) * len(sizes) * int(mean + 2) + len(OR)
    if tether is None:
        tether = (np.full(len(sizes), -1, np.int64), np.zeros((len(sizes), 3)), 0.0)
    if collide_nt > 0 and linker_nt is not None:
        q_rel, t_rel = np.asarray(disome[0], np.float64), np.asarray(disome[1], np.float64)
    else:
        q_rel, t_rel = np.array([1.0, 0, 0, 0]), np.zeros(3)
    if priors is None:
        npri = np.zeros(len(sizes), np.int64); PP = np.zeros((1, 3))
    else:
        npri = np.asarray(priors[0], np.int64); PP = np.asarray(priors[1], np.float64).reshape(-1, 3)
    poff = np.r_[0, np.cumsum(npri)[:-1]].astype(np.int64)
    ws = dict(WALK_SEARCH); ws.update(walk_search or {})
    cone = CONE_FLEX if linker_nt is not None else CONE
    wcone = cone if ws["dir_cone"] is None else float(ws["dir_cone"])
    wpick = {"nearest": 0, "first": 1, "compact": 2}[ws["pick"]]
    if ws["candidates"]:
        log(f"polysomes: walk candidate search {ws['candidates']} poses x {ws['rounds'] + 1} rounds per grown ribosome, "
            f"pick {ws['pick']}, direction cone {wcone:.2f} rad, entry orientation cone {float(ws['orient_cone']):.2f} rad")
    walk = kern.factory("polysomes.walk", K, _make)
    pos, qu, chn, coll, pri, ndrop = walk(seed, glo, ghi, *a, stc, str_, float(np.max(np.linalg.norm(stc, axis=1) + str_)),
                                         K.mask("cytosol") if mask is None else mask,
                                         sizes, ex, tex, en, q_en, OC.astype(np.float64), OR.astype(np.float64), ns + 1000, ns * 14 + 10000,
                                         lrms, cone, MAX_FAILS_FLEX if linker_nt is not None else 60,
                                         np.asarray(tether[1], np.float64), np.asarray(tether[0], np.int64), float(tether[2]),
                                         np.asarray(q_rel, np.float64), np.asarray(t_rel, np.float64),
                                         collide_nt if linker_nt is not None else 0.0, float(G[2] - 30.0 if R0 is None else R0),
                                         npri, poff, PP, float(nudge), int(prior_tries),
                                         int(ws["candidates"]), int(ws["rounds"]), wpick, wcone, float(ws["orient_cone"]),
                                         bool(ws["priors_first"]) and priors is not None)
    if ws["priors_first"] and priors is not None:     # phase order -> chain order (priors before grown within a chain)
        o = np.argsort(chn, kind="stable")
        pos, qu, chn, coll, pri = pos[o], qu[o], chn[o], coll[o], pri[o]
    # mRNA polylines
    Rm = Rotation.from_quat(np.c_[qu[:, 1:], qu[:, :1]].astype(np.float64)).as_matrix()
    pieces, links = [], []
    for i in range(len(pos)):
        P = ch @ Rm[i].T + pos[i]
        if i + 1 < len(pos) and chn[i + 1] == chn[i]:
            Pn = ch @ Rm[i + 1].T + pos[i + 1]
            lk = _hermite(P[0], Rm[i] @ tex, Pn[-1], -(Rm[i + 1] @ ten))
            links.append(lk[1:-1]); P = np.vstack([P[::-1], lk[1:]])
        else:
            P = P[::-1]
        pieces.append(_resample(P, mrna_step))
    off = np.r_[0, np.cumsum([len(p) for p in pieces])].astype(np.int64)
    n_chain = len(np.unique(chn))
    log(f"polysomes: {n_chain} chains, {len(pos)}/{sizes.sum()} ribosomes (mean {len(pos) / max(n_chain, 1):.1f} per mRNA), "
        f"{off[-1]} mRNA beads; {int(coll.sum())} ribosomes placed as collided behind their neighbour")
    if priors is not None:
        log(f"polysomes: point priors {int(pri.sum())} placed at (or within {nudge:.0f} nm of) their positions, {ndrop} dropped "
            f"(no fit within {prior_tries} tries)")
    out = dict(poly_pos=pos, poly_q=qu, poly_chain=chn, poly_coll=coll, li_ribo=li, poly_mrna=np.concatenate(pieces).astype(np.float32),
               poly_mrna_off=off, poly_link=np.concatenate(links).astype(np.float32) if links else np.zeros((0, 3), np.float32),
               poly_chains=n_chain)
    return out, pri, int(ndrop)


# ------------------------------------------------------------------ transcripts (origin: mrna_genes.py)
def transcripts(rng, g, n_ribo, n_total=None, log=print, P=None):
    """one row per mRNA molecule from a per-gene table g (origin: mrna_genes.transcripts; constants in P)"""
    if P["mrna_mode"] == "inferred" and n_total is None:
        return _transcripts_inferred(rng, g, n_ribo, log, P)
    if n_total is None:
        n_total = P["n_total_default"]
    scale = n_total / g.copies.sum()
    lam = g.copies.values * scale
    n = np.floor(lam).astype(np.int64) + (rng.random(len(lam)) < lam - np.floor(lam))      # stochastic rounding
    m = g.loc[np.repeat(np.arange(len(g)), n)].reset_index(drop=True)
    w = (m.cds * m.te).values
    f0 = 0.0
    if P["occupancy_mean"]:
        gw = (g.cds * g.te).values[n > 0]; sp = w.sum() / max(n_ribo, 1.0)
        for _ in range(30):
            occ1 = (1.0 - np.exp(-gw / sp)).mean()
            f0 = float(np.clip(1.0 - P["occupancy_mean"] / occ1, 0.0, 0.9))
            sp = (1.0 - f0) * w.sum() / max(n_ribo, 1.0)
    if P["spacing_nt"]:
        spacing = P["spacing_nt"] * float(g.te.values[n > 0].mean())
        occ1 = (1.0 - np.exp(-(g.cds * g.te).values[n > 0] / spacing)).mean()
        f0 = float(np.clip(1.0 - P["occupancy_mean"] / occ1, 0.0, 0.9)) if P["occupancy_mean"] else 0.0
    act = rng.random(len(m)) >= f0
    if not P["spacing_nt"]:
        spacing = (w * act).sum() / max(n_ribo, 1.0)
    m["ribosomes"] = np.where(act, np.minimum(P["max_ribo"], rng.poisson(w / spacing)), 0)
    m["length"] = (m.utr5 + m.cds + m.utr3 + m.polya).round().astype(np.int64)
    log(f"mRNA genes: {len(m)} mRNAs of {int((n > 0).sum())} genes (data total {g.copies.sum():.0f}, scale {scale:.2f}); "
        f"untranslated share {f0:.3f}; ribosome spacing {spacing:.0f} nt of CDS; {int((m.ribosomes > 0).sum())} translated, mean {m.ribosomes[m.ribosomes > 0].mean():.1f} "
        f"ribosomes, {int(m.ribosomes.sum())} on mRNAs ({m.ribosomes.sum() / max(n_ribo, 1.0):.2f} x the budget of {n_ribo:.0f}); "
        f"median length {int(m.length.median())} nt")
    return m


def _transcripts_inferred(rng, g, n_ribo, log, P):
    p = (g.copies / g.copies.sum()).values; cds = g.cds.values.astype(float); te = g.te.values.astype(float)
    spacing = float((p * cds * te).sum() / (P["density_cw"] / 100.0 * (p * cds).sum()))
    lam = np.minimum(cds * te / spacing, P["max_ribo"])
    f0 = float(np.clip(1.0 - P["occupancy_cw"] / (p * (1.0 - np.exp(-lam))).sum(), 0.0, 0.9))
    per_molecule = (1.0 - f0) * float((p * lam).sum())
    n_total = n_ribo / per_molecule
    copies = g.copies.values * n_total / g.copies.sum()
    n = np.floor(copies).astype(np.int64) + (rng.random(len(copies)) < copies - np.floor(copies))
    m = g.loc[np.repeat(np.arange(len(g)), n)].reset_index(drop=True)
    act = rng.random(len(m)) >= f0
    m["ribosomes"] = np.where(act, np.minimum(P["max_ribo"], rng.poisson((m.cds * m.te).values / spacing)), 0)
    m["length"] = (m.utr5 + m.cds + m.utr3 + m.polya).round().astype(np.int64)
    tr = m.ribosomes > 0
    log(f"mRNA genes (inferred total): {len(m)} mRNAs of {int((n > 0).sum())} genes ({len(m) * P['pop_mean_vol'] / P['cell_vol']:.0f} "
        f"per average cell; Miura 2008 36,139, Zenklusen 2008 about 60,000); untranslated share {f0:.3f}, on ribosomes {tr.mean():.3f}; "
        f"spacing {spacing:.0f} nt of CDS at TE 1; {int(tr.sum())} translated, mean {m.ribosomes[tr].mean():.2f} ribosomes; "
        f"{int(m.ribosomes.sum())} on mRNAs of a budget of {n_ribo:.0f}; median length {int(m.length.median())} nt")
    return m


def transcripts_table(rng, T, n_poly, occupancy, max_ribo, log=print):
    """B4 (Jurkat): one row per mRNA of a frozen transcript set T (gene, length nt); a share `occupancy` of them is
    translated; ribosomes per translated mRNA ~ Poisson(length / spacing) (capped), the spacing solved (fixed point)
    so that the expected ribosomes in chains of >= 2 equal n_poly (the polysome budget)."""
    m = T[["gene", "length"]].copy().reset_index(drop=True)
    act = rng.random(len(m)) < occupancy
    L = m.length.values.astype(float)
    sp = (L[act]).sum() / max(n_poly, 1.0)
    for _ in range(60):                       # E[k; k >= 2] = lam - lam exp(-lam) for k ~ Poisson(lam)
        lam = np.minimum(L[act] / sp, max_ribo)
        e2 = float((lam - lam * np.exp(-lam)).sum())
        sp *= e2 / n_poly
    lam = np.minimum(L / sp, max_ribo)
    m["ribosomes"] = np.where(act, np.minimum(max_ribo, rng.poisson(lam)), 0)
    m["cds"] = m.length; m["utr5"] = 0.0
    k = m.ribosomes.values
    log(f"mRNA table: {len(m)} mRNAs, translated share {act.mean():.3f}; spacing {sp:.0f} nt per ribosome; "
        f"{int((k >= 2).sum())} chains of >= 2 holding {int(k[k >= 2].sum())} ribosomes (budget {n_poly:.0f}); "
        f"{int((k == 1).sum())} single ribosomes; median length {int(m.length.median())} nt")
    return m


# ------------------------------------------------------------------ point priors (B4)
def point_priors(ctx, pp, log):
    """EM ribosome centres in the model frame, de-duplicated, with single-linkage groups and an ER flag"""
    import pandas as pd
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    from wholecell.preplaced.stage import _resolve
    path = _resolve(ctx.cfg, ctx.bdir, pp["table"])
    E = pd.read_parquet(path)
    for c, v in (pp.get("filter_eq") or {}).items():
        E = E[E[c] == v]
    for c, v in (pp.get("filter_not_in") or {}).items():
        E = E[~E[c].isin(v)]
    centre = np.asarray(ctx.geom.centre_nm, np.float64)
    P = E[["x", "y", "z"]].values * np.asarray(pp.get("voxel_nm", [4.0, 4.0, 3.44])) - centre
    er = E[pp.get("class_col", "cls")].isin(pp.get("er_classes", [])).values
    n0 = len(P)
    t = cKDTree(P); keep = np.ones(len(P), bool)
    for i, j in sorted(t.query_pairs(float(pp.get("dedup_nm", 20.0)))):
        if keep[i] and keep[j]:
            keep[j] = False
    P, er = P[keep], er[keep]
    pr = cKDTree(P).query_pairs(float(pp.get("link_nm", 34.0)), output_type="ndarray")
    A = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(P), len(P)))
    ng, lab = connected_components(A, directed=False)
    log(f"point priors {os.path.basename(path)}: {n0:,} centres, {n0 - len(P):,} duplicates (< {pp.get('dedup_nm', 20.0)} nm) merged "
        f"-> {len(P):,}; {ng:,} groups at {pp.get('link_nm', 34.0)} nm, {int((np.bincount(lab) == 1).sum()):,} singletons; "
        f"ER class {int(er.sum()):,}")
    return P, lab, er


def _order_group(P):
    """visit order of a group's centres: nearest-neighbour path from an end (farthest point from the centroid)"""
    if len(P) <= 2:
        return np.arange(len(P))
    from scipy.spatial.distance import cdist
    D = cdist(P, P); left = set(range(len(P)))
    cur = int(np.argmax(np.linalg.norm(P - P.mean(0), axis=1))); order = [cur]; left.remove(cur)
    while left:
        L = list(left); cur = L[int(np.argmin(D[cur, L]))]; order.append(cur); left.remove(cur)
    return np.array(order)


def assign_priors(P, lab, er, sizes, rng, max_size):
    """chains (translated mRNAs, ribosome counts `sizes`) seeded by EM groups: groups (largest first, split at max_size)
    take the largest chains (rank matched, a chain grows to at least its group size), singletons seed the remaining
    chains in random order, leftover singletons become monosome chains. Returns (order of chains, sizes, npri,
    prior positions, er flag per chain, n monosome chains appended)"""
    groups = [np.where(lab == g)[0] for g in range(lab.max() + 1)]
    multi = [g for g in groups if len(g) >= 2]; single = np.array([g[0] for g in groups if len(g) == 1], np.int64)
    parts = []
    for g in sorted(multi, key=len, reverse=True):
        g = g[_order_group(P[g])]
        for s in range(0, len(g), max_size):
            parts.append(g[s:s + max_size])
    parts.sort(key=len, reverse=True)
    ci = np.argsort(-np.asarray(sizes), kind="stable")          # chains by size, largest first
    sizes = np.asarray(sizes, np.int64).copy()
    order, npri, pts, erf = [], [], [], []
    used = 0
    for k, g in enumerate(parts):
        if k >= len(ci):
            break
        c = ci[k]; order.append(c); sizes[c] = max(sizes[c], len(g)); npri.append(len(g)); pts.append(P[g])
        erf.append(bool(er[g].mean() >= 0.5)); used += 1
    rest = ci[used:]
    perm = rng.permutation(len(single))
    ns = min(len(rest), len(single))
    for k in range(ns):
        order.append(rest[k]); npri.append(1); pts.append(P[single[perm[k]]][None]); erf.append(bool(er[single[perm[k]]]))
    for c in rest[ns:]:
        order.append(c); npri.append(0); erf.append(False)
    mono = single[perm[ns:]]
    return (np.array(order, np.int64), sizes, np.array(npri, np.int64), (np.concatenate(pts) if pts else np.zeros((0, 3))),
            np.array(erf, bool), mono)


# ------------------------------------------------------------------ the builder
def _pool(ctx, params, ribo):
    exp3 = ctx.exp3
    r = exp3.lib_id == ribo
    if params.get("pool_region", "cytosol") is not None:
        r &= exp3.region == params.get("pool_region", "cytosol")
    return r


def _constants(params, ctx):
    d = dict(mrna_mode="inferred", density_cw=0.468, occupancy_cw=0.773, occupancy_mean=0.71, spacing_nt=156.0, max_ribo=40,
             n_total_default=None, pop_mean_vol=1.0, cell_vol=1.0)
    d.update(params.get("mrna_model") or {})
    return d


@register("preplaced", "core.polysomes")
def build(ctx, params):
    from wholecell import registry
    cell, lp, exp3, log, K = ctx.cell, ctx.lp, ctx.exp3, ctx.log, ctx.K
    ribo = params.get("ribosome", "ribosome_80S")
    beads = np.load(ctx.lib_path(ribo))["beads"]
    ch = np.load(_data_path(ctx, params.get("channel", f"lib:{ribo}_mrna_channel.npy"))).astype(np.float64)
    dz = np.load(_data_path(ctx, params.get("disome", f"lib:{ribo}_disome_6I7O.npy")))
    disome = (dz[:4], dz[4:7])
    obs_names = params.get("obstacles", ["other", "cortex", "fg"])
    OCs = [ctx.obs[k][0] for k in obs_names if k in ctx.obs]; ORs = [ctx.obs[k][1] for k in obs_names if k in ctx.obs]
    rib = _pool(ctx, params, ribo)
    pool = float(exp3.loc[rib, "expected"].sum())
    frac = float(params.get("polysome_fraction", POLY_FRACTION))
    tx = tx_chains = tx_free = None
    T = params.get("transcripts") or {}
    priors = mono = None; erf = None
    if T.get("provider"):      # gene-specific mRNAs: each translated mRNA is one polysome chain of its own size
        prov = registry.get("preplaced_data", T["provider"])(ctx, T)
        ctx.state["gene_table"] = prov["genes"]          # per-gene poly(A) etc. for the mRNP builder
        P = _constants(params, ctx); P.update(prov.get("constants") or {})
        tx = transcripts(np.random.default_rng(ctx.seed(params.get("transcripts_seed", 23))), prov["genes"], frac * pool,
                         log=log, P=P)
        tx_chains, tx_free = tx[tx.ribosomes > 0].reset_index(drop=True), tx[tx.ribosomes == 0].reset_index(drop=True)
    elif T.get("table"):
        import pandas as pd
        TT = registry.get("preplaced_data", T["table"])(ctx, T)
        tx = transcripts_table(np.random.default_rng(ctx.seed(params.get("transcripts_seed", 23))), TT, frac * pool,
                               float(T.get("occupancy", 0.8)), int(T.get("max_ribo", 40)), log=log)
        tx_chains, tx_free = tx[tx.ribosomes > 0].reset_index(drop=True), tx[tx.ribosomes == 0].reset_index(drop=True)
    exz = params.get("exclusion_zones")       # ribosome exclusion zones (yeast: endocytic sites), a registered provider
    EXZ = registry.get("preplaced_data", exz)(ctx, {}) if exz else np.zeros((0, 4))
    obstacles = (np.concatenate(OCs + [EXZ[:, :3]]) if OCs or len(EXZ) else np.zeros((0, 3)),
                 np.concatenate(ORs + [EXZ[:, 3]]) if ORs or len(EXZ) else np.zeros(0))
    n_ribo = frac * pool
    n_mrna = int(params.get("n_mrna", N_MRNA))
    sizes = tx_chains.ribosomes.values if tx is not None else None
    linker_nt = None
    lnt = params.get("linker_nt_mean")       # B4 addendum: per-cell linker (None = from the transcripts, yeast)
    if tx is not None:
        fp = float(params.get("footprint_nt", FOOTPRINT_NT))
        linker_nt = (tx_chains.cds / tx_chains.ribosomes - fp).values
        if lnt is not None:
            linker_nt = np.full(len(tx_chains), float(lnt))
    tether = None
    if tx is not None and params.get("tether"):
        tether = registry.get("preplaced_data", params["tether"])(ctx, dict(chains=tx_chains, rng=np.random.default_rng(
            ctx.seed(params.get("tether_seed", 29)))))
    pp = params.get("point_priors")
    if pp:
        if tx is None:
            raise ValueError("core.polysomes: point priors need transcripts (chains)")
        PEM, lab, er = point_priors(ctx, pp, log)
        rngp = np.random.default_rng(ctx.seed(pp.get("seed", "preplaced.polysomes.priors")))
        order, sizes2, npri, PP, erf, mono = assign_priors(PEM, lab, er, tx_chains.ribosomes.values, rngp, int(pp.get("max_group", 40)))
        tx_chains = tx_chains.iloc[order].reset_index(drop=True)
        tx_chains["ribosomes"] = sizes2[order]
        # monosomes at the remaining EM singletons: chains of one ribosome, not on an mRNA (appended last)
        sizes = np.r_[tx_chains.ribosomes.values, np.ones(len(mono), np.int64)]
        npri = np.r_[npri, np.ones(len(mono), np.int64)]
        PP = np.concatenate([PP, PEM[mono]]) if len(mono) else PP
        linker_nt = np.r_[(tx_chains.cds / tx_chains.ribosomes - float(params.get("footprint_nt", FOOTPRINT_NT))).values
                          if lnt is None else np.full(len(tx_chains), float(lnt)), np.full(len(mono), 5.0)]
        if lnt is not None:
            log(f"polysomes: linker {float(lnt):.1f} free nt per link (mean; params.linker_nt_mean) for every chain")
        erf = np.r_[erf, er[mono]]
        tm = np.where(erf, K.mem(pp.get("er_membrane", "ER")), -1).astype(np.int64)
        tether = (tm, np.zeros((len(sizes), 3)), float(pp.get("er_tmax_nm", 45.0)))
        priors = (npri, PP)
        log(f"point priors: {int((npri[:len(order)] > 1).sum()):,} chains seeded by EM groups, "
            f"{int((npri[:len(order)] == 1).sum()):,} by EM singletons, {int((npri[:len(order)] == 0).sum()):,} unseeded; "
            f"{len(mono):,} EM singletons as monosomes; {int(erf.sum()):,} ER-class chains tethered (grown ribosomes)")
    R0 = None
    if params.get("start_radius_slot") is not None:
        R0 = float(cell.G[int(params["start_radius_slot"])] - 30.0)
    poly, pri, ndrop = place(K, cell, lp, ribo, beads, ch, disome, n_ribo, ctx.seed(params.get("seed", 17)), obstacles=obstacles,
                             log=log, n_mrna=n_mrna, sizes=sizes, linker_nt=linker_nt, tether=tether,
                             mask=K.mask(*params.get("labels", ["cytosol"])), R0=R0, priors=priors,
                             collide_nt=float(params.get("collide_nt", COLLIDE_NT)), nudge=float((pp or {}).get("nudge_nm", 10.0)),
                             prior_tries=int((pp or {}).get("prior_tries", 60)), walk_per=int(params.get("walk_tree_per", 550)),
                             mrna_step=float(params.get("mrna_step_nm", NT_STEP)), walk_search=params.get("walk"))
    ctx.pre.update(poly)
    ctx.state.update(tx=tx, tx_chains=tx_chains, tx_free=tx_free, poly=poly, n_mrna=n_mrna, channel=ch)
    if mono is not None and len(mono):
        # EM monosomes (chains appended after the mRNA chains, so their ribosomes are the trailing ones) carry no mRNA:
        # the mRNA builder sees the mRNA chains only; the monosomes are an obstacle set of their own
        nk = int(np.searchsorted(poly["poly_chain"], len(tx_chains)))
        off = poly["poly_mrna_off"]
        view = dict(poly, poly_pos=poly["poly_pos"][:nk], poly_q=poly["poly_q"][:nk], poly_chain=poly["poly_chain"][:nk],
                    poly_coll=poly["poly_coll"][:nk], poly_mrna=poly["poly_mrna"][:off[nk]], poly_mrna_off=off[:nk + 1],
                    poly_chains=int(len(np.unique(poly["poly_chain"][:nk]))))
        ctx.state["poly"] = view
        # coarse walk tree (550 residues per sphere), not the library tree (3,509 spheres per 80S: memory)
        from scipy.spatial.transform import Rotation
        rc, rr = sphere_tree(beads.astype(np.float64), per=int(params.get("walk_tree_per", 550)))
        Rm = Rotation.from_quat(np.c_[poly["poly_q"][nk:, 1:], poly["poly_q"][nk:, :1]].astype(np.float64)).as_matrix()
        ctx.obs["monosomes"] = ((np.einsum("nij,sj->nsi", Rm, rc) + poly["poly_pos"][nk:, None, :]).reshape(-1, 3),
                                np.tile(rr, len(Rm)))
    if hasattr(ctx, "groups"):              # core context: instance group + consumption, polysome statistics
        from wholecell.preplaced.context import add_group
        chn = poly["poly_chain"]
        cnt = np.bincount(chn, minlength=len(sizes))
        n_in = int(cnt[chn][cnt[chn] >= 2].size)
        erc = erf[chn] if erf is not None else np.zeros(len(chn), bool)
        add_group(ctx, "polysomes", np.full(len(chn), poly["li_ribo"], np.int32), poly["poly_pos"], poly["poly_q"],
                  chain=chn.astype(np.int32), origin=np.where(pri, 0, 1).astype(np.int8), er=erc.astype(np.int8),
                  collided=poly["poly_coll"], consumed={ribo: int(len(chn))})
        st = dict(pool=pool, placed=int(len(chn)), in_chains_ge2=n_in, polysome_fraction_of_pool=n_in / max(pool, 1.0),
                  priors_placed=int(pri.sum()), priors_dropped=ndrop, er_chain_ribosomes=int(erc.sum()),
                  chains=int(len(np.unique(chn))), packed_free=float(pool - len(chn)))
        ctx.state["polysome_stats"] = st
        log(f"polysomes: {n_in:,} ribosomes in chains of >= 2 = {st['polysome_fraction_of_pool']:.3f} of the pool {pool:,.0f}; "
            f"{len(chn):,} pre-placed, {pool - len(chn):,.0f} left to the packer")


def _data_path(ctx, v):
    """'lib:<file>' = a file in the library folder of the context (yeast layout: the channel / disome caches next to
    the entries); 'repo:<rel>' = a repository file; else a config path (build: / frozen names resolved)"""
    if v.startswith("lib:"):
        return os.path.join(ctx.lib_dir, v[4:])
    if v.startswith("repo:"):
        from wholecell.config import REPO
        return os.path.join(REPO, v[5:])
    from wholecell.preplaced.stage import _resolve
    return _resolve(ctx.cfg, ctx.bdir, v)
