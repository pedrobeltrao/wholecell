# origin: yeast pipeline/wc3d.py:prepass (chromatin part, from "cn = cell.G[3:6]" to the returned dict) @ e3b7448
# Changed: own registered serial builder (was the second half of wc3d.prepass, after the before_chromatin callback);
# the NPCs come from pre (core.npc_pores), the extra obstacles / start points from obs["nucleus"] /
# state["chrom_starts"] (the spindle builders) as the callback returned them. rabl_segments, fiber_steps, the fibre
# walks and the nucleosome map are the vendored functions (unchanged); stream number is a param.
"""Rabl chromatin (serial pre-placed builder `yeast.rabl_chromatin`).

Chromosome arms as nucleosome fibres from the kinetochores (centromeres at the SPB) to telomere foci at the nuclear
envelope, XII-R through the nucleolus; one nucleosome per measured dyad (nucleosome map), walk order = genome order.
Stream `seed` (5: the same stream number as the NPC spins, as in wc3d.prepass): rabl layout rng seed+100, fibre walks
seed. Adds pre: li_nuc, nuc_pos, nuc_q, nuc_chain, nuc_brk, chain_target, seg_start, seg_target, seg_chrom, seg_kind,
nuc_bp.
"""
import time

import numpy as np

from wholecell.registry import register


@register("preplaced", "yeast.rabl_chromatin")
def build(ctx, params):
    wc3d = ctx.mod("wc3d"); c3 = ctx.mod("cell3d"); p3 = ctx.mod("pack3d")
    cell, lp, log = ctx.cell, ctx.lp, ctx.log
    seed = ctx.seed(params.get("seed", 5))
    li_npc = int(ctx.pre["li_npc"]); li_nuc = int(np.where(lp["ids"] == "nucleosome")[0][0])
    npc_pos, npc_q = ctx.pre["npc_pos"], ctx.pre["npc_q"]
    types = np.full(len(npc_pos), li_npc, np.int32)
    cn = cell.G[3:6]; rn = float(cell.G[11:14].max())
    glo = cn - rn - 80; ghi = cn + rn + 80
    OC, OR = p3.emit_spheres(types, npc_pos.astype(np.float32), npc_q.astype(np.float32), lp["st_off"], lp["STC"], lp["STR"], lp["BR"], glo, ghi)
    n_nuc = ctx.state["n_nuc"]
    lens = np.array(wc3d.CHROM_KB, float); lens = np.round(lens / lens.sum() * n_nuc).astype(np.int64)
    dirn = cell.G[7:10] - cn; dirn /= np.linalg.norm(dirn)
    spb = c3.nuc_surface(cell.G, -dirn, -(c3.PERI_W + 40))[0]
    starts, r_start = np.tile(spb, (len(lens), 1)), 150.0
    extra = {}
    if "nucleus" in ctx.obs:                  # obstacles placed before the chromatin (spindle)
        C, R = ctx.obs["nucleus"]; k = np.all((C > glo) & (C < ghi), axis=1)
        OC = np.concatenate([OC, C[k]]); OR = np.concatenate([OR, R[k]])
    st = ctx.state.get("chrom_starts")
    if st is not None and len(st) == len(lens):     # chromosome walks start at the kinetochore-MT plus ends
        starts, r_start = np.asarray(st, np.float64), 20.0
    a = cell.args(glo, ghi, margin=200.0)
    s0, s1 = lp["st_off"][li_nuc], lp["st_off"][li_nuc + 1]
    mask = 1 << c3.NUCP
    nsph = (s1 - s0) * lens.sum() + len(OR)
    t0 = time.time()
    sl, ss, st_, sc, sk = wc3d.rabl_segments(cell, lens, starts, np.random.default_rng(seed + 100), log=log)
    dy = None
    if wc3d.FIBER and wc3d.NUC_MAP:
        nucleosome_map, genome_map = ctx.mod("nucleosome_map"), ctx.mod("genome_map")
        dy = nucleosome_map.segment_dyads(genome_map.segments())
        sl = np.array([len(p) for p in dy], np.int64)
        log(f"nucleosome map: {sl.sum()} measured dyads on the arms (inventory {int(lens.sum())})")
    if wc3d.FIBER:
        sp = None if dy is None else np.concatenate([np.r_[wc3d.NRL_BP, np.abs(np.diff(p))] for p in dy])
        delta, lk = wc3d.fiber_steps(sl, sp)
        pos, qu, chain, brk = p3.fiber_walks(seed, glo, ghi, *a, lp["STC"][s0:s1], lp["STR"][s0:s1], lp["BR"][li_nuc], mask,
                                             OC, OR, sl, ss, st_, r_start, delta, wc3d.FIBER_LP, wc3d.FIBER_R,
                                             lk, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    else:
        pos, qu, chain, brk = p3.arm_walks(seed, glo, ghi, *a, lp["STC"][s0:s1], lp["STR"][s0:s1], lp["BR"][li_nuc], mask,
                                           OC, OR, sl, ss, st_, r_start, 12.5, int(nsph * 1.1) + 1000, int(nsph * 10) + 10000)
    log(f"chromatin: {len(pos)}/{sl.sum()} nucleosomes, {int(brk.sum())} fragments, {time.time()-t0:.0f}s")
    if dy is not None:          # walks place their nucleosomes in order: the first n of each segment's dyads
        extra["nuc_bp"] = np.concatenate([dy[c][:int((chain == c).sum())] for c in range(len(dy))]).astype(np.int64)
        assert len(extra["nuc_bp"]) == len(pos) and np.all(np.diff(chain) >= 0)
    ctx.pre.update(li_nuc=li_nuc, nuc_pos=pos, nuc_q=qu, nuc_chain=chain, nuc_brk=brk, chain_target=sl,
                   seg_start=ss, seg_target=st_, seg_chrom=sc, seg_kind=sk.astype("U8"))
    ctx.pre.update(extra)
