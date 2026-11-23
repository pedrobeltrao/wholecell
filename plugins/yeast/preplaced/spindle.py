# origin: yeast scripts/s3_pack.py:main (closure microtubules(): kinetochores, outer kinetochore, SPB plaques; the
#         struct_counts / n_ccan / n_mind deductions after wc3d.prepass) @ e3b7448
# Changed: three registered serial builders sharing state["mtd"] and obs["nucleus"] as the closure shared C, R, mtd;
# stream numbers and the ring Spc110 share are params. Numerics: pipeline/kinetochores.py, spb.py (unchanged).
"""Yeast spindle structures (serial pre-placed builders `yeast.kinetochores`, `yeast.spb`, `yeast.spindle_pools`).

yeast.kinetochores: CCAN + MIND on the 16 kinetochore-MT plus ends (the chromosome arms then start beyond them:
  state["chrom_starts"]), then the outer kinetochore (DASH rings, Ndc80 complexes) from their pools.
yeast.spb: the SPB plaques (Spc110 ... Cmd1 layers) scaled to the Spc110 pool minus the gamma-TuRC attachment share.
yeast.spindle_pools: deducts the placed kinetochore/SPB subunits (after the chromosome arms, as stage 3 did).
"""
import numpy as np

from wholecell.registry import register
from wholecell.preplaced.builders.microtubules import spindle_to_pre


def _pool(exp3):
    return lambda genes: {g: float(exp3.loc[exp3.lib_id == g, "expected"].sum()) for g in genes}


@register("preplaced", "yeast.kinetochores")
def kinetochores(ctx, params):
    kin = ctx.mod("kinetochores")
    exp3, mtd = ctx.exp3, ctx.state["mtd"]
    npc_spheres = ctx.obs["npc"]; C, R = ctx.obs["nucleus"]
    pcc = (exp3.lib_id == kin.CCAN_ID) & (exp3.region == "nucleoplasm"); pmi = (exp3.lib_id == kin.MIND_ID) & (exp3.region == "nucleoplasm")
    kt = kin.place(ctx.cell, ctx.lp, mtd["kmt_tips"], mtd["kmt_dirs"], (np.concatenate([C, npc_spheres[0]]), np.concatenate([R, npc_spheres[1]])),
                   float(exp3.loc[pcc, "expected"].sum()), float(exp3.loc[pmi, "expected"].sum()), seed=params.get("seed", 61), log=ctx.log)
    starts = mtd["kmt_tips"]
    if kt:
        mtd.update({k: v for k, v in kt.items() if k.startswith("kt_")}); mtd.update(n_ccan=kt["n_ccan"], n_mind=kt["n_mind"])
        C, R = np.concatenate([C, kt["spheres"][0]]), np.concatenate([R, kt["spheres"][1]])
        starts = kt["starts"]
    ko, ko_n, ko_sph = kin.outer(ctx.cell, ctx.lp, mtd["kmt_tips"], mtd["kmt_dirs"], (np.concatenate([C, npc_spheres[0]]),
                                 np.concatenate([R, npc_spheres[1]])), _pool(exp3)(kin.DASH + kin.NDC80C),
                                 seed=params.get("outer_seed", 67), log=ctx.log)
    C, R = np.concatenate([C, ko_sph[0]]), np.concatenate([R, ko_sph[1]])
    mtd.update(ko)
    ctx.obs["nucleus"] = (C, R)
    ctx.state["chrom_starts"] = starts
    ctx.state["struct_counts"] = dict(ko_n)
    spindle_to_pre(ctx)


@register("preplaced", "yeast.spb")
def spb_plaques(ctx, params):
    spb = ctx.mod("spb")
    cell, exp3, mtd = ctx.cell, ctx.exp3, ctx.state["mtd"]
    npc_spheres = ctx.obs["npc"]; C, R = ctx.obs["nucleus"]; CA, RA = ctx.obs["mt_all"]
    per_ring = params.get("ring_spc110_per_nuclear_ring", 14)
    nrm = cell.G[7:10] - cell.G[3:6]
    ps_ = _pool(exp3)(spb.LAYERS); ps_["SPC110"] = max(0.0, ps_["SPC110"] - per_ring * mtd["n_rings_nuc"])     # minus the ring share
    sp, sp_n, sp_sph = spb.place(cell, ctx.lp, -nrm / np.linalg.norm(nrm), ps_,
                                 (np.concatenate([CA, C, npc_spheres[0]]), np.concatenate([RA, R, npc_spheres[1]])),
                                 seed=params.get("seed", 71), log=ctx.log, ring_spc110=per_ring * mtd["n_rings_nuc"])
    ctx.obs["nucleus"] = (np.concatenate([C, sp_sph[0]]), np.concatenate([R, sp_sph[1]]))
    mtd.update(sp)
    ctx.state.setdefault("struct_counts", {}).update(sp_n)
    spindle_to_pre(ctx)


@register("preplaced", "yeast.spindle_pools")
def spindle_pools(ctx, params):
    kin = ctx.mod("kinetochores")
    exp3, mtd = ctx.exp3, ctx.state["mtd"]
    for g_, n_ in ctx.state.get("struct_counts", {}).items():
        ctx.consume(list(exp3.loc[exp3.lib_id == g_, "ing_id"].unique()), n_)
    pcc = (exp3.lib_id == kin.CCAN_ID) & (exp3.region == "nucleoplasm"); pmi = (exp3.lib_id == kin.MIND_ID) & (exp3.region == "nucleoplasm")
    for m_, k_ in ((pcc, "n_ccan"), (pmi, "n_mind")):
        if k_ in mtd:
            exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - mtd[k_]).clip(lower=0)
