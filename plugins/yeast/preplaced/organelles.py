# origin: yeast scripts/s3_pack.py:main (blocks "Cvt complex next to the vacuole", "nucleus-vacuole junction",
#         "RAC rides on ribosomes") @ e3b7448
# Changed: registered serial builders; stream numbers are params. Numerics: pipeline/cvt.py, nvj.py, rac.py
# (unchanged).
"""Yeast organelle-contact and ribosome structures (serial pre-placed builders).

yeast.cvt   Cvt complex next to the vacuole: precursor Ape1 dodecamers and Ams1 (from the vacuolar pools), Atg19 coat.
            pre: cvt_*
yeast.nvj   nucleus-vacuole junction: Vac8-Nvj1 on the vacuolar membrane facing the nucleus (all Nvj1, one Vac8
            each). pre: nvj_*, li_nvj
yeast.rac   RAC on polysome ribosomes (pre poly_type) and on free ribosomes in proportion; deducts the polysome
            ribosomes from the cytosolic pool and adds the free ribosome-RAC rows (same regions)
"""
import numpy as np
import pandas as pd

from wholecell.registry import register


@register("preplaced", "yeast.cvt")
def cvt_complex(ctx, params):
    cvt = ctx.mod("cvt")
    exp3 = ctx.exp3
    OC, OR = ctx.obs["other"]; MC, MR = ctx.obs["cortex"]; FGC, FGR = ctx.obs["fg"]
    tot = lambda i: float(exp3.loc[exp3.ing_id == i, "expected"].sum())
    cv, used = cvt.place(ctx.cell, ctx.lp, tot("APE1"), tot("AMS1"), int(tot("ATG19")),
                         obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])), seed=params.get("seed", 61), log=ctx.log)
    if cv:
        ctx.pre.update(cv)
        ctx.consume(["APE1"], used[cvt.APE1]); ctx.consume(["AMS1"], used[cvt.AMS1]); ctx.consume(["ATG19"], used[cvt.ATG19])
        KC, KR = ctx.emit(cv["cvt_type"], cv["cvt_pos"], cv["cvt_q"])
        ctx.obs["other"] = (np.concatenate([OC, KC]), np.concatenate([OR, KR]))


@register("preplaced", "yeast.nvj")
def nvj_junction(ctx, params):
    nvj = ctx.mod("nvj")
    exp3 = ctx.exp3
    OC, OR = ctx.obs["other"]; MC, MR = ctx.obs["cortex"]; FGC, FGR = ctx.obs["fg"]
    tot = lambda i: float(exp3.loc[exp3.ing_id == i, "expected"].sum())
    nv, n_nvj = nvj.place(ctx.cell, ctx.lp, tot("NVJ1"), obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])),
                          seed=params.get("seed", 63), log=ctx.log)
    if n_nvj:
        ctx.pre.update(nv); ctx.consume(["NVJ1"], n_nvj); ctx.consume(["VAC8"], n_nvj)
        JC, JR = ctx.emit(np.full(n_nvj, nv["li_nvj"], np.int32), nv["nvj_pos"], nv["nvj_q"])
        ctx.obs["other"] = (np.concatenate([OC, JC]), np.concatenate([OR, JR]))


@register("preplaced", "yeast.rac")
def rac_on_ribosomes(ctx, params):
    rac = ctx.mod("rac"); ps = ctx.mod("polysome")
    exp3, pre, poly = ctx.exp3, ctx.pre, ctx.state["poly"]
    OC, OR = ctx.obs["other"]; MC, MR = ctx.obs["cortex"]; FGC, FGR = ctx.obs["fg"]
    rib = (exp3.lib_id == ps.RIBO) & (exp3.region == "cytosol")
    # RAC rides on ribosomes (all RAC complexes, spread over polysome and free ribosomes in proportion)
    n_rib = float(exp3.loc[rib, "expected"].sum())
    racr = exp3.ing_id == rac.RAC_ING; n_rac = float(exp3.loc[racr, "expected"].sum()); f_rac = min(1.0, n_rac / max(n_rib, 1.0))
    pre["poly_type"] = rac.decorate(ctx.cell, ctx.lp, poly, f_rac * len(poly["poly_pos"]),
                                    obstacles=(np.concatenate([OC, MC, FGC]), np.concatenate([OR, MR, FGR])),
                                    seed=params.get("seed", 71), log=ctx.log)
    on_poly = int((pre["poly_type"] != poly["li_ribo"]).sum())
    exp3.loc[rib, "expected"] = exp3.loc[rib, "expected"] - len(poly["poly_pos"])
    n_free = float(exp3.loc[rib, "expected"].sum())
    free_rac = min(n_rac - on_poly, f_rac * n_free)
    new = exp3.loc[rib].assign(lib_id=rac.LIB_ID, ing_id=rac.LIB_ID, expected=lambda d: d.expected * free_rac / max(n_free, 1.0))
    exp3.loc[rib, "expected"] *= 1.0 - free_rac / max(n_free, 1.0)
    ctx.consume([rac.RAC_ING], on_poly + free_rac)
    ctx.exp3 = pd.concat([exp3, new], ignore_index=True)        # free ribosomes with RAC, in the same regions
    ctx.log(f"  RAC: {n_rac:.0f} complexes for {n_rib:.0f} cytosolic ribosomes; {on_poly} on polysome ribosomes, {free_rac:.0f} on free ones")
