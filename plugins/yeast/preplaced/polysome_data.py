# origin: yeast scripts/s3_pack.py:main (polysome block: mrna_genes.transcripts inputs, wc3d.exclusion_zones and its
#         log line, mrna_localization.tether) @ e3b7448
# Changed (B4, 2026-11-17): the yeast data the common polysome builder (wholecell/preplaced/builders/polysomes.py)
# reads, as registered `preplaced_data` providers; values and log lines unchanged.
"""Yeast data of the common polysome builder: the per-gene mRNA table and the mRNA-model constants
(pipeline/mrna_genes.py: Pelechano 2010 copies, SGD CDS, TIF-seq UTRs, PAL-seq poly(A), TE; Arava 2003 density and
occupancy), the ribosome exclusion zones of the endocytic sites (wc3d.exclusion_zones: Picco 2018) and the
membrane tethers of the polysomes (mrna_localization.tether: Williams 2014 mitochondrial and Jan 2014 ER genes on
the cell3d primitives). These are yeast data and yeast-geometry seeding rules, read from the vendored modules."""
from wholecell.registry import register


@register("preplaced_data", "yeast.mrna_gene_table")
def gene_table(ctx, params):
    mg = ctx.mod("mrna_genes"); inv = ctx.mod("inventory")
    n_total = mg.MRNA_PER_CELL * inv.CELL_VOL / inv.POP_MEAN_VOL if mg.MRNA_PER_CELL else None
    return dict(genes=mg.gene_table(),
                constants=dict(mrna_mode=mg.MRNA_MODE, density_cw=mg.DENSITY_CW, occupancy_cw=mg.OCCUPANCY_CW,
                               occupancy_mean=mg.OCCUPANCY_MEAN, spacing_nt=mg.SPACING_NT, max_ribo=mg.MAX_RIBO,
                               n_total_default=n_total, pop_mean_vol=inv.POP_MEAN_VOL, cell_vol=inv.CELL_VOL))


@register("preplaced_data", "yeast.ribosome_exclusion")
def exclusion(ctx, params):
    wc3d = ctx.mod("wc3d")
    EXZ = wc3d.exclusion_zones(ctx.cell)
    ctx.log(f"ribosome exclusion zones: {len(EXZ)} spheres of {wc3d.EXCL_R:.0f} nm around {len(ctx.cell.endocytic_sites())} endocytic sites")
    return EXZ


@register("preplaced_data", "yeast.mrna_tether")
def tether(ctx, params):
    return ctx.mod("mrna_localization").tether(ctx.cell, params["chains"], params["rng"], log=ctx.log)
