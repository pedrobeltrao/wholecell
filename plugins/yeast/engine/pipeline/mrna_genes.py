# origin: yeast yeast_cell_dev_bundle/pipeline/mrna_genes.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Gene-specific mRNAs (IMPROVEMENTS M1-M2): every mRNA molecule is a transcript of a named gene, with the gene's
length (5' UTR + CDS + 3' UTR + poly(A)) and a ribosome load proportional to its CDS.

Inputs
  copies per gene   mRNA molecules per cell (RA) of Pelechano et al. 2010 PLoS One 5:e15442, Table S1 (compiled from
                    Nagalakshmi et al. 2008 and Miura et al. 2008); 5,085 genes, 25,900 mRNAs per cell
  CDS               summed CDS segments of each SGD ORF (R64), so introns are left out (until 2026-11-12 the ORF span,
                    which counted the introns of 208 genes, mostly ribosomal-protein genes, as coding sequence)
  UTRs, poly(A)     per gene: TIF-seq UTRs in YPD (Pelechano et al. 2013 Nature 497:127) and PAL-seq tails (Subtelny
                    et al. 2014 Nature 508:66); genes without data take UTR5_NT, UTR3_NT, POLYA_NT. TIF-seq UTRs
                    are genomic distances, so a 5' UTR intron (SGD five_prime_UTR_intron, 21 genes here) is
                    subtracted where the UTR is longer than the intron (the transcript starts upstream of it)
  ribosome load     ribosomes on an mRNA ~ Poisson(CDS x TE / spacing), TE the per-gene translational efficiency from
                    ribosome profiling (Subtelny et al. 2014). Since 2026-11-13 (round U3) the spacing is the measured
                    mean density, one ribosome per 156 nt of ORF averaged over genes (Arava et al. 2003), and the number
                    of ribosomes in polysomes is an outcome; before, one spacing was chosen so that the mRNAs carried an
                    assumed 70% of the cytosolic ribosomes (polysome.POLY_FRACTION), which needed 122 nt once a quarter
                    of the mRNAs were untranslated. SPACING_NT = None restores that.
  untranslated      (2026-11-13, round U3) on top of that, every mRNA is untranslated with one probability for all genes,
                    chosen so that the mean ribosome occupancy over genes is 71%, the average of polysome profiling
                    over 5,701 genes (Arava et al. 2003 PNAS 100:3889). Without it only 3% of mRNAs were free, since
                    a Poisson draw with a mean of 7.6 is rarely zero. Only that one mean is used: the per-gene
                    occupancies and the 14% of mRNAs with a single ribosome stay a test (validation/mrna_occupancy.py).
Output of transcripts(): DataFrame, one row per mRNA molecule: orf, gene, cds, utr5, utr3, polya, length, ribosomes.
"""
import os
import numpy as np
import pandas as pd

from config import DATA

PELECHANO = f"{DATA}/pelechano2010/files/pone.0015442.s002.xls"
SGD = f"{DATA}/sgd/SGD_features.tab"
UTR_TABLE = f"{DATA}/mrna/utr_per_gene.tsv"          # TIF-seq, YPD (Pelechano 2013); data/mrna/README.md
POLYA_TABLE = f"{DATA}/mrna/polya_per_gene.tsv"      # PAL-seq, total RNA (Subtelny 2014)
TE_TABLE = f"{DATA}/mrna/te_per_gene.tsv"            # ribosome profiling TE (Subtelny 2014)
UTR5_NT, UTR3_NT, POLYA_NT = 50, 120, 33     # where no per-gene value exists (near the data medians)
SPACING_NT = 156.0       # nt of ORF per ribosome, mean over genes (Arava 2003: 0.64 ribosomes per 100 nt); None: budget
# mRNAs per cell (2026-11-13, round U3): Miura et al. 2008 (BMC Genomics 9:574; S288C, YPD, competitive PCR) count
# 36,139 per average cell; most mRNAs scale with cell size in G1 (Chen et al. 2020 Mol Cell 78:359; Swaffer et al. 2021),
# so the G1 cell gets that times its volume share (inventory.CELL_VOL / POP_MEAN_VOL). Pelechano 2010 gives the genes'
# shares. Before: Pelechano's own total (25,900) times the protein share 0.70 = 18,100, which with the measured
# occupancy and density put only 36% of the ribosomes in polysomes (82-85% by polysome profiles: Arava 2003 in YPD,
# Metzl-Raz 2017 in SC). Other totals: about 60,000 (Zenklusen 2008, by rescaling an older 15,000), 12,200 (von der Haar 2008).
MRNA_PER_CELL = 36139.0
# Inferred total (2026-11-15, MRNA_MODE "inferred", the default): ribosomes on mRNAs = polysome share x cytosolic
# ribosomes (polysome.POLY_FRACTION, measured 0.82-0.85), and the mRNAs are as many as it takes to carry them at the
# measured load. Arava et al. 2003 Table S3 (5,643 genes, with copies per cell): averaged over translated molecules
# 0.468 ribosomes per 100 nt of ORF (one per 214 nt; the gene average is 0.517, the text's 0.64 is not reproduced by the
# table) and 5.35 ribosomes per translated mRNA (Riba et al. 2019, independently: 5.5); averaged over molecules 77.3%
# of mRNAs are on ribosomes (gene average 71.7%). With these the G1 cell needs about 50,000 mRNAs, about 61,000 per
# average cell: near Zenklusen et al. 2008 (about 60,000) and 1.7 x Miura et al. 2008 (36,139), which become tests.
# MRNA_MODE "count" uses MRNA_PER_CELL x the volume share (round U, 2026-11-13: 29,200 mRNAs, 57% of ribosomes on
# mRNAs, 7.2 per translated mRNA).
MRNA_MODE = "inferred"
DENSITY_CW = 0.468       # ribosomes per 100 nt of ORF, averaged over translated molecules (Arava 2003 Table S3)
OCCUPANCY_CW = 0.773     # share of mRNA molecules on ribosomes, averaged over molecules (Arava 2003 Table S3)
OCCUPANCY_MEAN = 0.71    # mean over genes of the share of a gene's mRNAs on ribosomes (Arava 2003); None: off
MAX_RIBO = 40             # cap on ribosomes per mRNA (the walk stops earlier when it runs out of space)
GENE_MRNA = True          # False: the earlier generic mRNAs (lognormal lengths, Poisson polysome sizes)


def gene_table():
    import warnings
    warnings.filterwarnings("ignore")
    p = pd.read_excel(PELECHANO, engine="xlrd")
    p = pd.DataFrame(dict(orf=p.iloc[:, 0].astype(str).str.upper(), gene=p.iloc[:, 1].astype(str),
                          copies=pd.to_numeric(p.iloc[:, 2], errors="coerce"))).dropna(subset=["copies"])
    sgd = pd.read_csv(SGD, sep="\t", header=None, dtype=str)

    def seg_nt(kind):          # summed length of a feature's segments per parent ORF
        x = sgd[sgd[1] == kind]
        n = (pd.to_numeric(x[10], errors="coerce") - pd.to_numeric(x[9], errors="coerce")).abs() + 1
        return n.groupby(x[6].str.upper().values).sum()
    f = sgd[sgd[1] == "ORF"]
    a, b = pd.to_numeric(f[9], errors="coerce"), pd.to_numeric(f[10], errors="coerce")
    span = pd.Series((abs(b - a) + 1).values, index=f[3].str.upper().values)
    cds = pd.DataFrame(dict(orf=span.index, cds=seg_nt("CDS").reindex(span.index).fillna(span).values))
    g = p.merge(cds, on="orf", how="inner")
    # PORT PATCH 5 (2026-11-17): pandas >= 3 keeps NaN through astype(str); test isna() too (same result as pandas 2)
    g["gene"] = np.where(g.gene.isin(["nan", ""]) | g.gene.isna(), g.orf, g.gene)
    g["utr5"], g["utr3"], g["polya"], g["te"] = float(UTR5_NT), float(UTR3_NT), float(POLYA_NT), 1.0
    for tab, cols in ((UTR_TABLE, ("utr5", "utr3")), (POLYA_TABLE, ("polya",)), (TE_TABLE, ("te",))):
        if tab is not None and os.path.exists(tab):
            t = pd.read_csv(tab, sep="\t")
            t["orf"] = t.orf.str.upper()
            g = g.drop(columns=list(cols)).merge(t[["orf", *cols]], on="orf", how="left")
            for c in cols:
                g[c] = g[c].fillna({"utr5": UTR5_NT, "utr3": UTR3_NT, "polya": POLYA_NT, "te": 1.0}[c])
    ui = g.orf.map(seg_nt("five_prime_UTR_intron")).fillna(0.0)
    g["utr5"] = np.where(g.utr5 > ui, g.utr5 - ui, g.utr5)
    return g


def transcripts(rng, n_ribo, n_total=None, log=print):
    """one row per mRNA molecule. n_ribo: ribosomes to put on mRNAs (the polysome budget; used only with SPACING_NT
    None, and for the log line); n_total: scale the
    per-gene copies to this many mRNAs (None: MRNA_PER_CELL times the cell's volume share)"""
    import inventory
    g = gene_table()
    if MRNA_MODE == "inferred" and n_total is None:
        return _transcripts_inferred(rng, g, n_ribo, log)
    if n_total is None:
        n_total = MRNA_PER_CELL * inventory.CELL_VOL / inventory.POP_MEAN_VOL if MRNA_PER_CELL else inventory.CELL_SHARE * g.copies.sum()
    scale = n_total / g.copies.sum()
    lam = g.copies.values * scale
    n = np.floor(lam).astype(np.int64) + (rng.random(len(lam)) < lam - np.floor(lam))      # stochastic rounding
    m = g.loc[np.repeat(np.arange(len(g)), n)].reset_index(drop=True)
    # one global spacing (nt of CDS per ribosome) so that the expected ribosomes on the translated mRNAs equal n_ribo
    w = (m.cds * m.te).values
    f0 = 0.0
    if OCCUPANCY_MEAN:
        # per-gene occupancy (1 - f0)(1 - exp(-lambda_g)); f0 and the spacing are solved together (fixed point)
        gw = (g.cds * g.te).values[n > 0]; sp = w.sum() / max(n_ribo, 1.0)
        for _ in range(30):
            occ1 = (1.0 - np.exp(-gw / sp)).mean()
            f0 = float(np.clip(1.0 - OCCUPANCY_MEAN / occ1, 0.0, 0.9))
            sp = (1.0 - f0) * w.sum() / max(n_ribo, 1.0)
    if SPACING_NT:          # measured density: the gene-average ribosomes per nt of ORF is 1 / SPACING_NT
        spacing = SPACING_NT * float(g.te.values[n > 0].mean())
        occ1 = (1.0 - np.exp(-(g.cds * g.te).values[n > 0] / spacing)).mean()
        f0 = float(np.clip(1.0 - OCCUPANCY_MEAN / occ1, 0.0, 0.9)) if OCCUPANCY_MEAN else 0.0
    act = rng.random(len(m)) >= f0
    if not SPACING_NT:
        spacing = (w * act).sum() / max(n_ribo, 1.0)
    m["ribosomes"] = np.where(act, np.minimum(MAX_RIBO, rng.poisson(w / spacing)), 0)
    m["length"] = (m.utr5 + m.cds + m.utr3 + m.polya).round().astype(np.int64)
    log(f"mRNA genes: {len(m)} mRNAs of {int((n > 0).sum())} genes (data total {g.copies.sum():.0f}, scale {scale:.2f}); "
        f"untranslated share {f0:.3f}; ribosome spacing {spacing:.0f} nt of CDS; {int((m.ribosomes > 0).sum())} translated, mean {m.ribosomes[m.ribosomes > 0].mean():.1f} "
        f"ribosomes, {int(m.ribosomes.sum())} on mRNAs ({m.ribosomes.sum() / max(n_ribo, 1.0):.2f} x the budget of {n_ribo:.0f}); "
        f"median length {int(m.length.median())} nt")
    return m


def _transcripts_inferred(rng, g, n_ribo, log=print):
    """MRNA_MODE "inferred": spacing from the molecule-averaged density, untranslated share from the molecule-averaged
    occupancy, and as many mRNAs as carry n_ribo ribosomes (expected values over the genes' shares of molecules)."""
    p = (g.copies / g.copies.sum()).values; cds = g.cds.values.astype(float); te = g.te.values.astype(float)
    spacing = float((p * cds * te).sum() / (DENSITY_CW / 100.0 * (p * cds).sum()))
    lam = np.minimum(cds * te / spacing, MAX_RIBO)
    f0 = float(np.clip(1.0 - OCCUPANCY_CW / (p * (1.0 - np.exp(-lam))).sum(), 0.0, 0.9))
    per_molecule = (1.0 - f0) * float((p * lam).sum())
    n_total = n_ribo / per_molecule
    copies = g.copies.values * n_total / g.copies.sum()
    n = np.floor(copies).astype(np.int64) + (rng.random(len(copies)) < copies - np.floor(copies))
    m = g.loc[np.repeat(np.arange(len(g)), n)].reset_index(drop=True)
    act = rng.random(len(m)) >= f0
    m["ribosomes"] = np.where(act, np.minimum(MAX_RIBO, rng.poisson((m.cds * m.te).values / spacing)), 0)
    m["length"] = (m.utr5 + m.cds + m.utr3 + m.polya).round().astype(np.int64)
    import inventory
    tr = m.ribosomes > 0
    log(f"mRNA genes (inferred total): {len(m)} mRNAs of {int((n > 0).sum())} genes ({len(m) * inventory.POP_MEAN_VOL / inventory.CELL_VOL:.0f} "
        f"per average cell; Miura 2008 36,139, Zenklusen 2008 about 60,000); untranslated share {f0:.3f}, on ribosomes {tr.mean():.3f}; "
        f"spacing {spacing:.0f} nt of CDS at TE 1; {int(tr.sum())} translated, mean {m.ribosomes[tr].mean():.2f} ribosomes; "
        f"{int(m.ribosomes.sum())} on mRNAs of a budget of {n_ribo:.0f}; median length {int(m.length.median())} nt")
    return m
