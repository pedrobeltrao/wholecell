# origin: yeast scripts/s3_pack.py:main (blocks "eisosome coats", "MCC marker proteins", "F-actin reserved in s1:
#         cables", "tropomyosin", "endocytic coat", "amphiphysin", "Arp2/3-branched patches") @ e3b7448
# Changed: registered serial builders (one per block), the growing obstacle sets of stage 3 (MC/MR, CC/CR, OC/OR) as
# ctx.obs["cortex"], ["actin"], ["other"]; stream numbers are params. Numerics: pipeline/eisosomes.py, mcc.py,
# cytoskeleton.py, amphiphysin.py, actin_patches.py, inventory.py (unchanged).
"""Yeast cell-cortex structures (serial pre-placed builders).

yeast.eisosomes       Pil1/Lsp1 coats in the PM furrows (every Pil1/Lsp1 form reduced by the share in coats).
                      pre: eis_*; obs["cortex"] = cytoplasmic MTs + coats
yeast.mcc             Sur7/Nce102 on the furrow walls (MCC_FRACTION of their PM pools). pre: mcc_*
yeast.actin_cables    F-actin cables (copies reserved in s1, f_actin.json) + Tpm1 along both long-pitch strands.
                      pre: act_*, tpm_*, li_actin, li_tpm; obs["actin"]
yeast.endocytic_coat  Sla2/Sla1/End3 at the invagination tips (Picco 2015). pre: coat_*
yeast.amphiphysin     Rvs161-Rvs167 on the deeper endocytic tubes. pre: rvs_*, li_rvs
yeast.actin_patches   Arp2/3-branched actin patches at the endocytic sites (capping protein); F-actin that does not
                      fit returns to the ACT1 monomer pool. pre: pat_*; obs["other"] = actin + patches
"""
import json

import numpy as np

from wholecell.registry import register


@register("preplaced", "yeast.eisosomes")
def eisosomes(ctx, params):
    eis = ctx.mod("eisosomes")
    exp3, oligo_map, log = ctx.exp3, ctx.oligo_map, ctx.log
    MC, MR = ctx.obs["mt_cyto"]
    # eisosome coats in the PM furrows, from the Pil1/Lsp1 pools (monomers, homodimers, heterodimer); the coat takes
    # the pools' Pil1 : Lsp1 ratio, so every form is reduced by the same fraction
    genes = ctx.ing.set_index("ing_id").genes.fillna("").str.split(";")
    units = {}
    for g in ("PIL1", "LSP1"):
        r = exp3.ing_id.map(lambda i: g in genes.get(i, [])).astype(bool)
        mult = exp3.ing_id.map(lambda i: float(oligo_map[i][1]) if i in oligo_map else 1.0)   # homo-oligomers: k per instance
        units[g] = float((exp3.loc[r, "expected"] * mult[r]).sum())
    eisd, used = eis.place(ctx.cell, ctx.lp, units["PIL1"] / 2, units["LSP1"] / 2, obstacles=(MC, MR), seed=params.get("seed", 53), log=log)
    ctx.pre.update(eisd)
    f_eis = 2.0 * (used["Pil1_dimer"] + used["Lsp1_dimer"]) / max(units["PIL1"] + units["LSP1"], 1.0)
    r = exp3.ing_id.map(lambda i: bool({"PIL1", "LSP1"} & set(genes.get(i, []))))
    exp3.loc[r, "expected"] *= 1.0 - f_eis
    log(f"  Pil1/Lsp1: {units['PIL1']:.0f} / {units['LSP1']:.0f} monomers, {f_eis:.1%} in eisosome coats")
    EC, ER_ = ctx.emit(eisd["eis_type"], eisd["eis_pos"], eisd["eis_q"])
    ctx.obs["cortex"] = (np.concatenate([MC, EC]), np.concatenate([MR, ER_]))          # obstacles for everything placed after


@register("preplaced", "yeast.mcc")
def mcc_proteins(ctx, params):
    mcc = ctx.mod("mcc")
    exp3 = ctx.exp3
    MC, MR = ctx.obs["cortex"]
    mr_ = exp3.ing_id.isin(mcc.GENES) & (exp3.region == "mem_PM")
    want = {lib: mcc.MCC_FRACTION * n for lib, n in exp3.loc[mr_].groupby("lib_id").expected.sum().items()}
    md_, used_m = mcc.place(ctx.cell, ctx.lp, want, obstacles=(MC, MR), seed=params.get("seed", 73), log=ctx.log)
    if len(md_["mcc_pos"]):
        ctx.pre.update(md_)
        for lib, k in used_m.items():
            r_ = mr_ & (exp3.lib_id == lib); exp3.loc[r_, "expected"] = (exp3.loc[r_, "expected"] - k).clip(lower=0)
        QC, QR = ctx.emit(md_["mcc_type"], md_["mcc_pos"], md_["mcc_q"])
        ctx.obs["cortex"] = (np.concatenate([MC, QC]), np.concatenate([MR, QR]))


@register("preplaced", "yeast.actin_cables")
def actin_cables(ctx, params):
    cs = ctx.mod("cytoskeleton")
    MC, MR = ctx.obs["cortex"]
    fa = json.load(open(f"{ctx.inv_dir}/f_actin.json"))
    cab = cs.place(ctx.cell, ctx.lp, int(fa["cables"]), seed=params.get("seed", 23), obstacles=(MC, MR), log=ctx.log)
    ctx.pre.update(cab)
    # tropomyosin (TPM1 coiled coils) along both long-pitch strands of the cable filaments
    tpm = cs.decorate(ctx.cell, ctx.lp, cab, obstacles=(MC, MR), log=ctx.log)
    ctx.pre.update(tpm)
    ctx.consume([cs.TPM_ING], 2 * len(tpm["tpm_pos"]))
    CC, CR = ctx.emit(np.r_[np.full(len(cab["act_pos"]), cab["li_actin"]), np.full(len(tpm["tpm_pos"]), tpm["li_tpm"])].astype(np.int32),
                      np.r_[cab["act_pos"], tpm["tpm_pos"]].astype(np.float32), np.r_[cab["act_q"], tpm["tpm_q"]].astype(np.float32))
    ctx.obs["actin"] = (CC, CR)
    ctx.state.update(f_actin=fa, n_cable_actin=len(cab["act_pos"]))


@register("preplaced", "yeast.endocytic_coat")
def endocytic_coat(ctx, params):
    amph = ctx.mod("amphiphysin")
    exp3 = ctx.exp3
    MC, MR = ctx.obs["cortex"]; CC, CR = ctx.obs["actin"]
    pools = {g: float(exp3.loc[exp3.ing_id == g, "expected"].sum()) for g in amph.COAT}
    co, used_c = amph.coat(ctx.cell, ctx.lp, pools, obstacles=(np.concatenate([MC, CC]), np.concatenate([MR, CR])),
                           seed=params.get("seed", 69), log=ctx.log)
    if len(co["coat_pos"]):
        ctx.pre.update(co)
        for g, k in used_c.items():
            ctx.consume([g], k)
        KC2, KR2 = ctx.emit(co["coat_type"], co["coat_pos"], co["coat_q"])
        ctx.obs["actin"] = (np.concatenate([CC, KC2]), np.concatenate([CR, KR2]))


RVS_FREE = ("RVS161", "RVS167")


def rvs_pool(exp3, het_id):
    """Rvs161-Rvs167 heterodimers available to the amphiphysin builder: the assembled heterodimer ingredient (het_id)
    plus pairs of the free monomers. The stage-1 heterodimer step writes het_id only when the two proteins share a
    compartment: with the legacy localisation (Rvs167 74% PM) it never did and the pool was the free monomers; with the
    B5 `common` preset both are cytosol and every Rvs167 is in het_id, so a monomer-only pool was 0 (fix 2026-11-19)."""
    het = float(exp3.loc[exp3.ing_id == het_id, "expected"].sum())
    free = min(float(exp3.loc[exp3.ing_id == g, "expected"].sum()) for g in RVS_FREE)
    return het + free


def consume_rvs(ctx, het_id, n):
    """take n dimers: from the heterodimer ingredient first, the rest as one Rvs161 + one Rvs167 each"""
    took = ctx.consume([het_id], n) if (ctx.exp3.ing_id == het_id).any() else 0.0
    rest = n - took
    if rest > 0.5:
        for g in RVS_FREE:
            ctx.consume([g], rest)


@register("preplaced", "yeast.amphiphysin")
def amphiphysin(ctx, params):
    amph = ctx.mod("amphiphysin")
    exp3 = ctx.exp3
    MC, MR = ctx.obs["cortex"]; CC, CR = ctx.obs["actin"]
    rv_pool = rvs_pool(exp3, amph.LIB_ID)
    rv, n_rv = amph.place(ctx.cell, ctx.lp, rv_pool, obstacles=(np.concatenate([MC, CC]), np.concatenate([MR, CR])),
                          seed=params.get("seed", 67), log=ctx.log)
    if n_rv:
        ctx.pre.update(rv); consume_rvs(ctx, amph.LIB_ID, n_rv)
        VC, VR = ctx.emit(np.full(n_rv, rv["li_rvs"], np.int32), rv["rvs_pos"], rv["rvs_q"])
        ctx.obs["actin"] = (np.concatenate([CC, VC]), np.concatenate([CR, VR]))


@register("preplaced", "yeast.actin_patches")
def actin_patches(ctx, params):
    patches = ctx.mod("actin_patches"); inv = ctx.mod("inventory")
    cell, exp3, log = ctx.cell, ctx.exp3, ctx.log
    MC, MR = ctx.obs["cortex"]; CC, CR = ctx.obs["actin"]
    fa = ctx.state["f_actin"]
    patches.model(dict(zip(*[inv.build()[0][c] for c in ("gene", "uniprot")])))
    avail = lambda i: float(exp3.loc[(exp3.ing_id == i) & (exp3.region == "cytosol"), "expected"].sum())
    pat, used = patches.place(cell, ctx.lp, fa["patches"], avail(patches.ARP23), avail(patches.CP_ING),
                              obstacles=(np.concatenate([MC, CC]), np.concatenate([MR, CR])), seed=params.get("seed", 41), log=log,
                              sites=[(p, u) for p, u, _ in cell.endocytic_sites()] or None)
    ctx.pre.update(pat)
    ctx.consume([patches.ARP23], used[patches.ARP23]); ctx.consume([patches.CP_ING], used[patches.CP])
    back = fa["cables"] + fa["patches"] - ctx.state["n_cable_actin"] - used[patches.ACTIN]
    act = (exp3.ing_id == "ACT1") & (exp3.region == "cytosol")
    if act.any():
        exp3.loc[act, "expected"] += back
    else:
        exp3.loc[len(exp3)] = ["ACT1", "ACT1", "cytosol", back]
    log(f"  F-actin: {ctx.state['n_cable_actin']} in cables, {used[patches.ACTIN]} in patches, {back:.0f} returned to the monomer pool")
    PC_, PR_ = ctx.emit(pat["pat_type"], pat["pat_pos"], pat["pat_q"])
    ctx.obs["other"] = (np.concatenate([CC, PC_]), np.concatenate([CR, PR_]))
