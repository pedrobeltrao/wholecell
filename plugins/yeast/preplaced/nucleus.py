# origin: yeast scripts/s3_pack.py:main (blocks "rDNA in the nucleolus", "transcribing Pol II", "licensed origins",
#         "chromatin proteins on their ChIP-exo sites", "tRNA genes ... silent chromatin") and pipeline/wc3d.py:
#         rdna_chromatin call @ e3b7448
# Changed: registered serial builders (one per block); stream numbers are params. Numerics: pipeline/rdna.py,
# wc3d.rdna_chromatin, pol2.py, origins.py, chromatin_sites.py, trna_telomeres.py (unchanged).
"""Yeast nuclear structures on the chromatin (serial pre-placed builders).

yeast.rdna_transcription  active rDNA genes in the nucleolus: Pol I, nascent pre-rRNA (90S on uncleaved transcripts),
                          DNA beads; engaged Pol I and tethered 90S leave their nucleolar pools. pre: rdna_*, prerrna_*,
                          rpol_*, r90_*, li_rpol, li_r90
yeast.rdna_chromatin      nucleosomal rDNA (inactive repeats, spacers) as walks in the nucleolus. pre: rnuc_*
yeast.pol2                engaged Pol II on its genes (Pelechano 2010 counts), nucleosomes evicted around each, DNA
                          through the cleft and nascent RNA. pre: pol2*, li_pol2; nuc_* thinned
yeast.origins             MCM double hexamer + ORC on confirmed origins (flanking nucleosomes evicted). pre: mcm_*,
                          orc_*, li_mcm, li_orc; nuc_* thinned
yeast.chromatin_sites     chromatin proteins on their ChIP-exo peaks (Rossi 2021), then tRNA genes / telomeres /
                          silent loci (Pol III, TFIIIC, Rap1, Sir2/3/4) on the sites the peaks leave out. pre: tf_*, lm_*
All need the measured nucleosome map (pre nuc_bp) except the rDNA builders.
"""
import numpy as np

from wholecell.registry import register


def _thin_nucleosomes(pre, keep):
    nb = pre["nuc_brk"].copy(); nb[np.r_[False, ~keep[:-1]] & keep] = 1
    for key in ("nuc_pos", "nuc_q", "nuc_chain", "nuc_bp"):
        if key in pre:
            pre[key] = pre[key][keep]
    pre["nuc_brk"] = nb[keep]


@register("preplaced", "yeast.rdna_transcription")
def rdna_transcription(ctx, params):
    if not ctx.state["n_rnuc"]:
        return
    rdna = ctx.mod("rdna")
    cell, exp3 = ctx.cell, ctx.exp3
    FGC, FGR = ctx.obs["fg"]
    gt0, gp0, gq0 = ctx.pre_instances()
    cn, rn = cell.G[3:6], float(cell.G[11:14].max())
    RC1, RR1 = ctx.emit(gt0, gp0, gq0, cn - rn - 100, cn + rn + 100)
    pol = (exp3.lib_id == "PolI") & (exp3.region == "nucleolus")
    r90 = (exp3.lib_id == "preribo_90S") & (exp3.region == "nucleolus")
    rd = rdna.place(cell, ctx.lp, float(exp3.loc[pol, "expected"].sum()), float(exp3.loc[r90, "expected"].sum()),
                    (np.concatenate([RC1, FGC]), np.concatenate([RR1, FGR])), seed=params.get("seed", 41), log=ctx.log)
    ctx.pre.update(rd)
    exp3.loc[pol, "expected"] = (exp3.loc[pol, "expected"] - len(rd["rpol_pos"])).clip(lower=0)
    exp3.loc[r90, "expected"] = (exp3.loc[r90, "expected"] - len(rd["r90_pos"])).clip(lower=0)


@register("preplaced", "yeast.rdna_chromatin")
def rdna_chromatin(ctx, params):
    n_rnuc = ctx.state["n_rnuc"]
    if not n_rnuc:
        return
    rdna = ctx.mod("rdna"); wc3d = ctx.mod("wc3d")
    cell, pre = ctx.cell, ctx.pre
    FGC, FGR = ctx.obs["fg"]
    cn, rn = cell.G[3:6], float(cell.G[11:14].max())
    gt1, gp1, gq1 = ctx.pre_instances()
    RC2, RR2 = ctx.emit(gt1, gp1, gq1, cn - rn - 100, cn + rn + 100)
    pre.update(wc3d.rdna_chromatin(cell, ctx.lp, n_rnuc, (np.concatenate([RC2, FGC, pre["rdna_pts"], pre["prerrna_pts"]]),
                                                         np.concatenate([RR2, FGR, np.full(len(pre["rdna_pts"]), 1.0),
                                                                         np.full(len(pre["prerrna_pts"]), rdna.RNA_R)])),
                                   n_chains=params.get("n_chains", 8), seed=params.get("seed", 9), log=ctx.log))


@register("preplaced", "yeast.pol2")
def pol2_engaged(ctx, params):
    pol2 = ctx.mod("pol2")
    pre, exp3 = ctx.pre, ctx.exp3
    FGC, FGR = ctx.obs["fg"]
    gtx, gpx, gqx = ctx.pre_instances({k: v for k, v in pre.items() if not k.startswith("nuc_")})
    XC, XR = ctx.emit(gtx, gpx, gqx)
    pr2 = (exp3.ing_id == "PolII") & (exp3.region == "nucleoplasm")
    r2 = pol2.place(ctx.cell, ctx.lp, pre, float(exp3.loc[pr2, "expected"].sum()), (np.concatenate([XC, FGC]), np.concatenate([XR, FGR])),
                    seed=params.get("seed", 91), log=ctx.log)
    _thin_nucleosomes(pre, r2.pop("keep")); pre.update(r2)
    exp3.loc[pr2, "expected"] = (exp3.loc[pr2, "expected"] - len(r2["pol2_pos"])).clip(lower=0)


@register("preplaced", "yeast.origins")
def licensed_origins(ctx, params):
    origins = ctx.mod("origins")
    pre, exp3 = ctx.pre, ctx.exp3
    FGC, FGR = ctx.obs["fg"]
    pmc = (exp3.lib_id == origins.MCM_ID) & (exp3.region == "nucleoplasm")
    por = (exp3.lib_id == origins.ORC_ID) & (exp3.region == "nucleoplasm")
    gtx, gpx, gqx = ctx.pre_instances({k: v for k, v in pre.items() if not k.startswith("nuc_")})
    XC, XR = ctx.emit(gtx, gpx, gqx)
    ro = origins.place(ctx.cell, ctx.lp, pre, float(exp3.loc[pmc, "expected"].sum()), float(exp3.loc[por, "expected"].sum()),
                       (np.concatenate([XC, FGC, pre["pol2dna_pts"], pre["pol2rna_pts"]]),
                        np.concatenate([XR, FGR, np.full(len(pre["pol2dna_pts"]) + len(pre["pol2rna_pts"]), 1.0)])),
                       seed=params.get("seed", 77), log=ctx.log)
    _thin_nucleosomes(pre, ro.pop("keep")); pre.update(ro)
    for m_, k_ in ((pmc, "mcm_pos"), (por, "orc_pos")):
        if k_ in ro:
            exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - len(ro[k_])).clip(lower=0)


@register("preplaced", "yeast.chromatin_sites")
def chromatin_sites(ctx, params):
    chs = ctx.mod("chromatin_sites"); ttl = ctx.mod("trna_telomeres")
    cell, lp, pre, exp3, log = ctx.cell, ctx.lp, ctx.pre, ctx.exp3, ctx.log
    FGC, FGR = ctx.obs["fg"]
    nucp = exp3.region == "nucleoplasm"
    gtx, gpx, gqx = ctx.pre_instances()
    XC, XR = ctx.emit(gtx, gpx, gqx)
    beads = np.concatenate([pre[k] for k in ("pol2dna_pts", "pol2rna_pts", "rdna_pts", "prerrna_pts") if k in pre])
    rt, tf_counts = chs.place(cell, lp, pre, exp3[nucp].groupby("lib_id").expected.sum().to_dict(),
                              (np.concatenate([XC, FGC, beads]), np.concatenate([XR, FGR, np.full(len(beads), 1.0)])),
                              seed=params.get("seed", 53), log=log)
    pre.update(rt)
    for g_, n_ in tf_counts.items():
        m_ = nucp & (exp3.lib_id == g_)
        exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - n_ * exp3.loc[m_, "expected"] / max(exp3.loc[m_, "expected"].sum(), 1e-9)).clip(lower=0)
    # tRNA genes (Pol III, TFIIIC), telomere ends (Rap1), silent chromatin (Sir2/3/4): sites the ChIP-exo peaks leave out
    if params.get("loci", True) and "tf_pos" in pre and "nuc_bp" in pre:
        nuc_any = exp3.region.isin(["nucleoplasm", "nucleolus"])
        TC, TR = ctx.emit(pre["tf_type"], pre["tf_pos"], pre["tf_q"])
        rl, lm_counts = chs.place(cell, lp, pre, exp3[nuc_any].groupby("lib_id").expected.sum().to_dict(),
                                  (np.concatenate([XC, TC, FGC, beads]), np.concatenate([XR, TR, FGR, np.full(len(beads), 1.0)])),
                                  sites=ttl.sites(pre), prefix="lm", what="tRNA genes, telomeres, silent loci",
                                  seed=params.get("loci_seed", 59), log=log)
        pre.update(rl)
        for g_, n_ in lm_counts.items():
            m_ = nuc_any & (exp3.lib_id == g_)
            exp3.loc[m_, "expected"] = (exp3.loc[m_, "expected"] - n_ * exp3.loc[m_, "expected"] / max(exp3.loc[m_, "expected"].sum(), 1e-9)).clip(lower=0)
