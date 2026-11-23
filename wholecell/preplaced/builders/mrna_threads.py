# origin: yeast scripts/s3_pack.py:main ("full-length mRNAs, threaded before packing" block) @ e3b7448
# origin: yeast pipeline/mrna.py:thread_all, place @ e3b7448
# origin: yeast pipeline/rac.py:rac_spheres @ e3b7448
# Changed (unification B4, 2026-11-17): the numerics live here (no ctx.mod); geometry primitives come from ctx.K
# (the cell's evaluator: yeast vendored cell3d/pack3d, core cellgeom_eval_v1/pack3d), c3.min_mem_dist + c3.HALF ->
# K.mmd(out, skip, G) + K.half(G) (same values for yeast); seeds are the integers of ctx.seed (place() no longer
# re-seeds through the yeast config); library files through ctx.lib_path; cell data (ribosome id, RAC id, bead
# representation, obstacle sets, allowed labels) are params with the yeast values as defaults. New: nt_per_bead
# (yeast 1: one bead per nt, 0.6 nm, radius 0.5; Jurkat 5: 3.0 nm, radius 1.0 = the v7 mRNA representation); with 1
# every number is the yeast one.
"""Full-length mRNA threads (serial pre-placed builder `core.mrna_threads`).

Polysome mRNAs through the ribosome mRNA channels plus their 5'/3' tails, and the untranslated mRNAs, as bead
threads avoiding every obstacle before them (cytoskeleton, organelle structures, NPCs, ribosome sphere trees, RAC at
the tunnel exits, polysome linkers, FG mesh). Adds pre: mrna_pts, mrna_off, mrna_kind, mrna_target, mrna_gene,
poly_chain_gene; state: the obstacle set (mrna_obstacles) the mRNP builder reuses, mrna_bead_r, mrna_nt_per_bead.

Inputs from ctx.state (polysome builder): poly (poly_pos, poly_q, poly_chain, poly_mrna, poly_mrna_off, poly_link,
poly_chains), n_mrna, tx / tx_chains / tx_free (transcript tables: gene, length, cds, utr5, ribosomes), channel (mRNA
channel path in the ribosome frame, 5'->3'; default: the library file <ribosome>_mrna_channel.npy).

params:
  seed                 stream (yeast 29)
  ribosome             library id of the ribosome (yeast ribosome_80S)
  rac                  true: RAC spheres at the tunnel exits of ribosomes typed `rac_id` in pre poly_type (yeast)
  rac_id               ribosome_80S_RAC
  obstacles            ordered obstacle sets: names of ctx.obs, or "@npc" (NPCs emitted from pre npc_pos/npc_q),
                       "@ribosomes" (polysome ribosome trees + RAC), "@links" (polysome linker beads, radius 0.5)
                       default (yeast order): [other, cortex, "@npc", "@ribosomes", "@links", fg]
  labels               allowed labels of the beads (default [cytosol])
  nt_per_bead          nucleotides per bead (default 1)
  step_nm, bead_r      bead spacing / radius (default 0.6 x nt_per_bead, 0.5 if nt_per_bead == 1 else 1.0)
  mean_nt, log_sigma, min_nt, persist, tail5_fraction   generic lengths (no transcript table)
  grid_h               collision grid cell of the threads (yeast 10 nm; core layout default 20 nm)
"""
import os

import numpy as np
import numba as nb
import pandas as pd

from wholecell.preplaced import kern
from wholecell.registry import register

# origin: yeast pipeline/mrna.py constants @ e3b7448
MEAN_NT = 1500.0
LOG_SIGMA = 0.5
MIN_NT = 300
STEP = 0.6            # nm per nucleotide bead
BEAD_R = 0.5
PERSIST = 1.0
TAIL5_FRACTION = (0.1, 0.4)
H = 10.0


def sphere_tree(xyz, per=30, seed=0):
    """structlib.sphere_tree (core port in _trees.py)"""
    from ._trees import sphere_tree as _st
    return _st(xyz, per=per, seed=seed)


def rac_spheres(ctx, pos, qu, typ, lp, rac_id="ribosome_80S_RAC", ribo="ribosome_80S"):
    """sphere tree of the RAC part of every decorated ribosome (types `typ`, ribosome frame), as obstacles for the
    mRNA threads that are placed after the decoration"""
    # origin: yeast pipeline/rac.py:rac_spheres @ e3b7448 (library paths through ctx.lib_path)
    from scipy.spatial.transform import Rotation
    ids = list(lp["ids"])
    if rac_id not in ids:
        return np.zeros((0, 3)), np.zeros(0)
    k = np.asarray(typ) == ids.index(rac_id)
    z = np.load(ctx.lib_path(rac_id)); n0 = len(np.load(ctx.lib_path(ribo))["beads"])
    rc, rr = sphere_tree(z["beads"][n0:].astype(np.float64), per=30)
    q = np.asarray(qu, np.float64)[k]
    Rm = Rotation.from_quat(np.c_[q[:, 1:], q[:, :1]]).as_matrix()
    C = (np.einsum("nij,sj->nsi", Rm, rc) + np.asarray(pos, np.float64)[k][:, None, :]).reshape(-1, 3)
    return C, np.tile(rr, int(k.sum()))


def _make(K):
    p3, c3, mmd, half = K.p3, K.c3, K.mmd, K.half
    NMEM = K.nmem

    @nb.njit(cache=False)
    def thread_all(seed, glo, ghi, G, pores, pidx, PT, PQ, PC, OC, OR, starts, dirs, lens, R0, cap_s, cap_e,
                   step, bead_r, persist, mask, maxd, nstraight, H):
        """One self-avoiding thread per row of lens. starts[i] with |dirs[i]| > 0: thread anchored there (polysome
        tails, FG domains; it first steps straight out along dirs[i] until clear, <= nstraight beads); dirs[i] == 0:
        random start in an allowed compartment. Beads must lie in compartments of the label bitmask `mask`, clear of
        membranes, and (maxd > 0) within maxd of the thread's start. H: collision grid cell (nm; yeast 10). Returns
        points and offsets."""
        # origin: yeast pipeline/mrna.py:thread_all @ e3b7448 (c3.min_mem_dist(out, -1) + c3.HALF -> mmd + half; the
        # grid cell H is an argument)
        np.random.seed(seed)
        dims, head, nxt, sid, cen, rad, cnt = p3.grid_new(glo, ghi, H, cap_s, cap_e)
        p3.g_insert(glo, H, dims, head, nxt, sid, cen, rad, cnt, OC, OR)
        n = lens.shape[0]
        pts = np.empty((lens.sum() + nstraight * n, 3), np.float32); off = np.zeros(n + 1, np.int64)
        out = np.zeros(NMEM); bead = np.empty((1, 3)); prev = np.empty((1, 3)); br = np.full(1, bead_r)
        nb_ = 0; nst = 0
        for c in range(n):
            anchored = dirs[c, 0] * dirs[c, 0] + dirs[c, 1] * dirs[c, 1] + dirs[c, 2] * dirs[c, 2] > 0.0
            ok = False
            if anchored:
                nst = 0
                for j in range(1, nstraight + 1):
                    for k in range(3):
                        bead[0, k] = starts[c, k] + step * j * dirs[c, k]
                    if not p3.g_overlap(glo, H, dims, head, nxt, sid, cen, rad, bead, br, 0.0):
                        ok = True; nst = j
                        break
                if ok:
                    for j in range(1, nst):
                        for k in range(3):
                            pts[nb_, k] = starts[c, k] + step * j * dirs[c, k]
                            prev[0, k] = pts[nb_, k]
                        if j <= nst - 2:
                            p3.g_insert(glo, H, dims, head, nxt, sid, cen, rad, cnt, prev, br)
                        nb_ += 1
                dx, dy, dz = dirs[c, 0], dirs[c, 1], dirs[c, 2]
            else:
                for at in range(500):
                    while True:
                        for k in range(3):
                            bead[0, k] = (np.random.random() * 2 - 1) * R0
                        if bead[0, 0] ** 2 + bead[0, 1] ** 2 + bead[0, 2] ** 2 < R0 * R0:
                            break
                    if p3.g_overlap(glo, H, dims, head, nxt, sid, cen, rad, bead, br, 0.0):
                        continue
                    m = c3.mem_sdfs(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, pidx, PT, PQ, PC, out)
                    if (mask >> c3.label_from(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, out, m)) & 1 and mmd(out, -1, G) > bead_r + half(G):
                        ok = True
                        break
                dx = np.random.normal(); dy = np.random.normal(); dz = np.random.normal()
                dn = np.sqrt(dx * dx + dy * dy + dz * dz); dx /= dn; dy /= dn; dz /= dn
            if not ok:
                off[c + 1] = nb_
                continue
            for k in range(3):
                pts[nb_, k] = bead[0, k]
            nb_ += 1; nthis = 1 + (nst - 1 if anchored else 0); stuck = 0; pending = False
            while nthis < lens[c] and stuck < 300:
                ex = persist * dx + np.random.normal(); ey = persist * dy + np.random.normal(); ez = persist * dz + np.random.normal()
                en = np.sqrt(ex * ex + ey * ey + ez * ez); ex /= en; ey /= en; ez /= en
                bead[0, 0] = pts[nb_ - 1, 0] + step * ex
                bead[0, 1] = pts[nb_ - 1, 1] + step * ey
                bead[0, 2] = pts[nb_ - 1, 2] + step * ez
                if maxd > 0.0:
                    ddx = bead[0, 0] - pts[off[c], 0]; ddy = bead[0, 1] - pts[off[c], 1]; ddz = bead[0, 2] - pts[off[c], 2]
                    if ddx * ddx + ddy * ddy + ddz * ddz > maxd * maxd:
                        stuck += 1; continue
                if p3.g_overlap(glo, H, dims, head, nxt, sid, cen, rad, bead, br, 0.0):
                    stuck += 1; continue
                m = c3.mem_sdfs(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, pidx, PT, PQ, PC, out)
                if not ((mask >> c3.label_from(bead[0, 0], bead[0, 1], bead[0, 2], G, pores, out, m)) & 1) or mmd(out, -1, G) < bead_r + half(G):
                    stuck += 1; continue
                if pending:
                    if not p3.g_insert(glo, H, dims, head, nxt, sid, cen, rad, cnt, prev, br):
                        off[c + 1] = nb_
                        return pts[:nb_], off
                for k in range(3):
                    prev[0, k] = pts[nb_ - 1, k]
                pending = True
                for k in range(3):
                    pts[nb_, k] = bead[0, k]
                nb_ += 1; nthis += 1
                dx, dy, dz = ex, ey, ez; stuck = 0
            for back in range(min(2, nthis)):
                for k in range(3):
                    prev[0, k] = pts[nb_ - 1 - back, k]
                p3.g_insert(glo, H, dims, head, nxt, sid, cen, rad, cnt, prev, br)
            off[c + 1] = nb_
        return pts[:nb_], off

    return thread_all


def thread_all(K):
    """the compiled thread kernel of kernel set K (shared with the FG mesh builder)"""
    return kern.factory("mrna_threads.thread_all", K, _make)


def domain(K, cell, grid_h=None):
    """whole-cell box (glo, ghi), radius R0 of the random thread starts and the collision grid cell of the threads:
    yeast X5 layout as the vendored code (G[0] enclosing radius, G[2] cell radius, 10 nm cells); core layout: G[0]
    (G_R, enclosing radius) for both (starts are rejected by label), grid cell grid_h (default 20 nm: the Jurkat cell
    is ~3x the yeast cell across, a 10 nm grid would need ~30 GB)"""
    G = cell.G
    glo = np.full(3, -G[0] - 10.0); ghi = -glo
    if K.name == "yeast_x5":
        return glo, ghi, G[2] - 30.0, float(H if grid_h is None else grid_h)
    return glo, ghi, float(G[0]), float(20.0 if grid_h is None else grid_h)


def lengths(rng, n, mean_nt=MEAN_NT, log_sigma=LOG_SIGMA, min_nt=MIN_NT):
    # origin: yeast pipeline/mrna.py:lengths @ e3b7448
    mu = np.log(mean_nt) - log_sigma ** 2 / 2
    return np.maximum(min_nt, rng.lognormal(mu, log_sigma, n)).astype(np.int64)


def place(ctx, poly, ch, n_free, obstacles, seed, log=print, tx_chains=None, tx_free=None, nt_per_bead=1, step=None,
          bead_r=None, persist=PERSIST, mask=None, tail5=TAIL5_FRACTION, mean_nt=MEAN_NT, log_sigma=LOG_SIGMA, min_nt=MIN_NT,
          grid_h=None):
    """poly: polysome dict; ch: mRNA channel path in the ribosome frame (5'->3'); n_free: free mRNAs; obstacles:
    (centres, radii) of everything pre-placed in the cytosol; seed: integer seed. tx_chains / tx_free: transcript
    tables (gene-specific lengths and 5'/3' tails; tx_chains in the order of the polysome sizes); chains whose walk
    placed no ribosome become free mRNAs. None: lognormal lengths, random 5' share. Lengths are in nucleotides;
    threads have one bead per nt_per_bead nucleotides."""
    # origin: yeast pipeline/mrna.py:place @ e3b7448 (seed given as an integer; bead representation params)
    from scipy.spatial.transform import Rotation
    K, cell = ctx.K, ctx.cell
    npb = int(nt_per_bead)
    step = STEP * npb if step is None else float(step)
    bead_r = (BEAD_R if npb == 1 else 1.0) if bead_r is None else float(bead_r)
    mask = (1 << K.lab("cytosol")) if mask is None else mask
    rng = np.random.default_rng(seed)
    tex = ch[0] - ch[3]; tex /= np.linalg.norm(tex)
    ten = ch[-1] - ch[-4]; ten /= np.linalg.norm(ten)
    chn = poly["poly_chain"]; Rm = Rotation.from_quat(np.c_[poly["poly_q"][:, 1:], poly["poly_q"][:, :1]].astype(np.float64)).as_matrix()
    off = poly["poly_mrna_off"]
    starts, dirs, lens, kind = [], [], [], []
    chains = np.unique(chn)
    genes = []
    if tx_chains is not None:
        lost = np.setdiff1d(np.arange(len(tx_chains)), chains)
        free = pd.concat([tx_free, tx_chains.iloc[lost]], ignore_index=True) if tx_free is not None else tx_chains.iloc[lost]
        n_free = len(free)
        L = np.r_[tx_chains.length.values[chains], free.length.values].astype(np.int64)
        genes = list(tx_chains.gene.values[chains]) + list(free.gene.values)
    else:
        L = lengths(rng, len(chains) + n_free, mean_nt, log_sigma, min_nt)
    drawn = 0
    for ci, c in enumerate(chains):
        idx = np.where(chn == c)[0]; first, last = idx[0], idx[-1]
        P = poly["poly_mrna"][off[first]:off[last + 1]]
        have = int(np.linalg.norm(np.diff(P, axis=0), axis=1).sum() / step) * npb      # nt through the ribosomes
        drawn += have
        rest = max(0, L[ci] - have)
        if tx_chains is not None:
            r = tx_chains.iloc[c]
            covered = min(r.cds, len(idx) * r.cds / max(r.ribosomes, 1))
            f5 = float(np.clip((r.utr5 + 0.5 * max(0.0, r.cds - covered)) / max(rest, 1), 0.0, 1.0))
        else:
            f5 = rng.uniform(*tail5)
        starts.append(Rm[first] @ ch[-1] + poly["poly_pos"][first]); dirs.append(Rm[first] @ ten)
        lens.append(int(rest * (1 - f5)) if npb == 1 else int(rest * (1 - f5) / npb)); kind.append(1)
        starts.append(Rm[last] @ ch[0] + poly["poly_pos"][last]); dirs.append(Rm[last] @ tex)
        lens.append(int(rest * f5) if npb == 1 else int(rest * f5 / npb)); kind.append(2)
    for i in range(n_free):
        starts.append(np.zeros(3)); dirs.append(np.zeros(3))
        lens.append(int(L[len(chains) + i]) if npb == 1 else int(L[len(chains) + i] / npb)); kind.append(0)
    starts = np.array(starts, np.float64).reshape(-1, 3); dirs = np.array(dirs, np.float64).reshape(-1, 3)
    lens = np.maximum(np.array(lens, np.int64), 1)
    OC, OR = obstacles
    glo, ghi, R0, h = domain(K, cell, grid_h)
    a = cell.args(glo, ghi, margin=50.0)
    cap_s = int(lens.sum()) + len(OR) + 1000
    pts, toff = thread_all(K)(seed, glo, ghi, *a, OC.astype(np.float64), OR.astype(np.float64), starts, dirs, lens, R0,
                              cap_s, int(cap_s * 3 + 8 * len(OR)), step, bead_r, persist, mask, 0.0, 20, h)
    got = np.diff(toff)
    kind = np.array(kind, np.int8)
    log(f"mRNA: {len(chains)} polysome mRNAs (+{got[kind > 0].sum() * npb / 1e6:.1f}M nt of tails to {drawn / 1e6:.1f}M nt through "
        f"ribosomes) and {n_free} free mRNAs ({got[kind == 0].sum() * npb / 1e6:.1f}M nt); threads reached "
        f"{got.sum() / lens.sum():.1%} of their target length; total {(got.sum() * npb + drawn) / 1e6:.1f}M nt")
    out = dict(mrna_pts=pts, mrna_off=toff, mrna_kind=kind, mrna_target=lens.astype(np.int32))
    if genes:
        g_thread = [g for g in genes[:len(chains)] for _ in (0, 1)] + genes[len(chains):]
        out["mrna_gene"] = np.array(g_thread, dtype="U12"); out["poly_chain_gene"] = np.array(genes[:len(chains)], dtype="U12")
    return out


def channel(ctx, ribo):
    """mRNA channel path of the ribosome (library frame, 5'->3'): ctx.state['channel'] or the library-side cache"""
    ch = ctx.state.get("channel")
    if ch is None:
        p = os.path.join(os.path.dirname(ctx.lib_path(ribo)), f"{ribo}_mrna_channel.npy")
        if not os.path.exists(p):
            raise FileNotFoundError(f"mrna_threads: no ctx.state['channel'] and no channel file {p}")
        ch = np.load(p)
    return np.asarray(ch).astype(np.float64)


DEFAULT_OBSTACLES = ["other", "cortex", "@npc", "@ribosomes", "@links", "fg"]


@register("preplaced", "core.mrna_threads")
def build(ctx, params):
    from scipy.spatial.transform import Rotation
    pre, lp, st = ctx.pre, ctx.lp, ctx.state
    ribo = params.get("ribosome", "ribosome_80S")
    poly = st["poly"]
    ch = channel(ctx, ribo)
    parts = []
    for name in params.get("obstacles", DEFAULT_OBSTACLES):
        if name == "@npc":
            if "npc_pos" in pre:
                npc_t = np.full(len(pre["npc_pos"]), pre["li_npc"], np.int32)
                parts.append(ctx.emit(npc_t, pre["npc_pos"], pre["npc_q"]))
        elif name == "@ribosomes":
            zr = np.load(ctx.lib_path(ribo))["beads"].astype(np.float64)
            rc, rr = sphere_tree(zr, per=int(params.get("ribosome_tree_per", 100)))   # Jurkat: coarser (memory)
            Rm = Rotation.from_quat(np.c_[poly["poly_q"][:, 1:], poly["poly_q"][:, :1]].astype(np.float64)).as_matrix()
            RC = (np.einsum("nij,sj->nsi", Rm, rc) + poly["poly_pos"][:, None, :]).reshape(-1, 3)
            rr_all = np.tile(rr, len(Rm))
            if params.get("rac", True) and "poly_type" in pre:          # RAC at the tunnel exits (yeast)
                RAC_C, RAC_R = rac_spheres(ctx, poly["poly_pos"], poly["poly_q"], pre["poly_type"], lp,
                                           params.get("rac_id", "ribosome_80S_RAC"), ribo)
                RC = np.concatenate([RC, RAC_C]); rr_all = np.concatenate([rr_all, RAC_R])
            parts.append((RC, rr_all))
        elif name == "@links":
            parts.append((poly["poly_link"], np.full(len(poly["poly_link"]), 0.5)))
        elif name in ctx.obs:
            parts.append(ctx.obs[name])
        elif not params.get("missing_ok", False):
            raise KeyError(f"core.mrna_threads: obstacle set {name!r} not in ctx.obs ({sorted(ctx.obs)})")
    obstacles = (np.concatenate([p[0] for p in parts]), np.concatenate([p[1] for p in parts]))
    npb = int(params.get("nt_per_bead", 1))
    mask = ctx.K.mask(*params.get("labels", ["cytosol"]))
    tx = st.get("tx")
    mr = place(ctx, poly, ch, st["n_mrna"] - poly["poly_chains"], obstacles, ctx.seed(params.get("seed", 29)), log=ctx.log,
               tx_chains=st["tx_chains"] if tx is not None else None, tx_free=st["tx_free"] if tx is not None else None,
               nt_per_bead=npb, step=params.get("step_nm"), bead_r=params.get("bead_r"),
               persist=params.get("persist", PERSIST), mask=mask, tail5=tuple(params.get("tail5_fraction", TAIL5_FRACTION)),
               mean_nt=params.get("mean_nt", MEAN_NT), log_sigma=params.get("log_sigma", LOG_SIGMA), min_nt=params.get("min_nt", MIN_NT),
               grid_h=params.get("grid_h"))
    pre.update(mr)
    st["mrna_obstacles"] = obstacles
    st["mrna_nt_per_bead"] = npb
    st["mrna_bead_r"] = params.get("bead_r", BEAD_R if npb == 1 else 1.0)
