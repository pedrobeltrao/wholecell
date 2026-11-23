# origin: yeast scripts/s3_pack.py:main ("FG-repeat domains of the nucleoporins" block) @ e3b7448
# origin: yeast pipeline/fg_mesh.py:place @ e3b7448
# Changed (unification B4, 2026-11-17): the numerics live here (no ctx.mod for the placement): the FG chains are
# threaded by the common thread kernel (mrna_threads.thread_all on ctx.K); seed = the integer of ctx.seed; the FG
# chain constants are params with the yeast values as defaults. The anchors (FG chain starts on the NPC scaffold) are a
# library-side cache, <library dir>/<npc lib>_fg_anchors.npz (yeast: NPC_fg_anchors.npz, made by the vendored
# fg_mesh.anchors() from 7N9F); when the cache is missing the yeast context's library builder makes it.
"""FG-repeat domains threaded into the pores (serial pre-placed builder `core.fg_mesh`).

Bead chains of the nucleoporin FG domains from their anchors on the NPC scaffold, avoiding every pre-placed instance
(NPC, chromatin, MTs, receptors). Adds pre: fg_pts, fg_off; obstacles: obs["fg"] (bead centres, bead radius).

params: seed (yeast 37), npc_lib (NPC), ca_step (0.38 nm), bead_r (0.3), persist (0.6), max_reach (60 nm),
labels (allowed: [nuclear_pore_channel, cytosol, nucleoplasm]), grid_h (thread grid cell; yeast 10 nm, core 20 nm).
"""
import os

import numpy as np

from wholecell.registry import register

# origin: yeast pipeline/fg_mesh.py constants @ e3b7448
CA_STEP = 0.38
BEAD_R = 0.3
PERSIST = 0.6
MAX_REACH = 60.0


def anchors(ctx, npc_lib="NPC"):
    """per FG chain of the NPC library entry: anchor point, inward direction, residues to add (library frame)"""
    p = os.path.join(os.path.dirname(ctx.lib_path(npc_lib)), f"{npc_lib}_fg_anchors.npz")
    if os.path.exists(p):
        z = np.load(p)
        return z["pos"], z["dir"], z["n"], z["gene"]
    if getattr(ctx, "_mod", None) is not None:          # yeast: the library-side cache builder (7N9F)
        return ctx.mod("fg_mesh").anchors()
    raise FileNotFoundError(f"core.fg_mesh: FG anchor cache {p} missing (library-side model of the NPC entry)")


def place(ctx, pre, obstacles, seed, log=print, npc_lib="NPC", ca_step=CA_STEP, bead_r=BEAD_R, persist=PERSIST,
          max_reach=MAX_REACH, mask=None, grid_h=None):
    # origin: yeast pipeline/fg_mesh.py:place @ e3b7448 (seed given as an integer)
    from scipy.spatial.transform import Rotation
    from .mrna_threads import thread_all
    K, cell = ctx.K, ctx.cell
    P, D, N, Gn = anchors(ctx, npc_lib)
    starts, dirs, lens = [], [], []
    for p, q in zip(pre["npc_pos"], pre["npc_q"]):
        R = Rotation.from_quat([q[1], q[2], q[3], q[0]]).as_matrix()
        starts.append(P @ R.T + p); dirs.append(D @ R.T); lens.append(N)
    starts = np.concatenate(starts).astype(np.float64); dirs = np.concatenate(dirs).astype(np.float64)
    lens = np.concatenate(lens).astype(np.int64)
    OC, OR = obstacles
    from .mrna_threads import domain
    G = cell.G
    glo, ghi, R0, h = domain(K, cell, grid_h)
    a = cell.args(glo, ghi, margin=50.0)
    cap_s = int(lens.sum()) + len(OR) + 1000
    if mask is None:
        mask = K.mask("nuclear_pore_channel", "cytosol", "nucleoplasm")
    pts, off = thread_all(K)(seed, glo, ghi, *a, OC.astype(np.float64), OR.astype(np.float64), starts, dirs, lens,
                             R0, cap_s, int(cap_s * 3 + 8 * len(OR)), ca_step, bead_r, persist, mask, max_reach, 60, h)
    got = np.diff(off)
    per = {g: (got[np.tile(Gn, len(pre['npc_pos'])) == g].sum(), lens[np.tile(Gn, len(pre['npc_pos'])) == g].sum()) for g in np.unique(Gn)}
    log(f"FG mesh: {len(lens)} FG domains on {len(pre['npc_pos'])} pores, {got.sum() / 1e6:.2f}M of {lens.sum() / 1e6:.2f}M residues "
        f"({got.sum() / lens.sum():.0%}); " + ", ".join(f"{g} {a_ / b_:.0%}" for g, (a_, b_) in per.items()))
    return dict(fg_pts=pts, fg_off=off)


@register("preplaced", "core.fg_mesh")
def build(ctx, params):
    gt0, gp0, gq0 = ctx.pre_instances()
    FC, FR = ctx.emit(gt0, gp0, gq0)
    labels = params.get("labels")
    mask = ctx.K.mask(*labels) if labels else None
    bead_r = float(params.get("bead_r", BEAD_R))
    ctx.pre.update(place(ctx, ctx.pre, (FC, FR), ctx.seed(params.get("seed", 37)), log=ctx.log,
                         npc_lib=params.get("npc_lib", "NPC"), ca_step=float(params.get("ca_step", CA_STEP)), bead_r=bead_r,
                         persist=float(params.get("persist", PERSIST)), max_reach=float(params.get("max_reach", MAX_REACH)),
                         mask=mask, grid_h=params.get("grid_h")))
    ctx.obs["fg"] = (ctx.pre["fg_pts"].astype(np.float64), np.full(len(ctx.pre["fg_pts"]), bead_r))
