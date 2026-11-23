# origin: yeast yeast_cell_dev_bundle/validation/mrna_occupancy.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Ribosome loading of the packed cell's mRNAs against polysome profiling and ribosome profiling (round U3, 2026-11-13).

The build uses two global means of Arava et al. 2003 (PNAS 100:3889): 71% of a gene's mRNAs on ribosomes on average,
and one ribosome per 156 nt of ORF (mrna_genes.OCCUPANCY_MEAN, SPACING_NT). What it does not use, and is tested here:
  per-gene occupancy        Spearman, model against Arava Table S3 (genes with >= 5 mRNAs in the model)
  per-gene density          ribosomes per 100 nt of ORF on translated mRNAs, Spearman against Arava
  single-ribosome mRNAs     share of all mRNA molecules with exactly one ribosome (Arava: 14%, from profiles and
                            abundances; the model's monosomes are chains of one)
  ribosomes per mRNA        per gene, Spearman against Riba et al. 2019 (PNAS 116:15023, Dataset S1: real_density x
                            orf_length, an independent ribosome-profiling estimate), and the copy-weighted means
  mRNA total                since 2026-11-15 inferred from the polysome share and the load (mrna_genes.MRNA_MODE): per
                            average cell against Miura et al. 2008 (36,139) and Zenklusen et al. 2008 (about 60,000)
  polysome share            the share of cytosolic 80S on mRNAs (an input since 2026-11-15; what placement achieves), against 82-85% by polysome
                            profiles (literature/ribosomes/README_polysome_share_mrna_count_2026-11-13.md)
Output: work/validation/mrna_occupancy.json
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import DATA, WORK, W3D  # noqa: E402
import mrna_genes as mg  # noqa: E402

ARAVA = f"{DATA}/validation/mrna/arava2003/arava2003_suppTable3_per_gene_occupancy_density.xls"
RIBA = f"{DATA}/validation/mrna/riba2019_density/pnas.1817299116.sd01.tsv"
MIN_MRNA = 5


def arava():
    warnings.filterwarnings("ignore")
    x = pd.read_excel(ARAVA, header=None)
    h = int(np.where(x.iloc[:, 0].astype(str).str.strip() == "YORF")[0][0])
    d = x.iloc[h + 1:, [0, 2, 6, 7]].copy(); d.columns = ["orf", "length", "density", "occupancy"]
    d["orf"] = d.orf.astype(str).str.strip().str.upper()
    for c in ("length", "density", "occupancy"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.dropna()


def model():
    z = np.load(f"{W3D}/prepass.npz")
    chn, cg, kind, mgene = z["poly_chain"], z["poly_chain_gene"], z["mrna_kind"], z["mrna_gene"]
    per = np.bincount(np.unique(chn, return_inverse=True)[1])        # ribosomes per chain, chain order
    free_gene = mgene[kind == 0]
    g = mg.gene_table()[["orf", "gene", "cds"]]
    name2orf = dict(zip(g.gene, g.orf)); name2orf.update({o: o for o in g.orf})
    rows = pd.DataFrame(dict(gene=np.r_[cg, free_gene], ribosomes=np.r_[per, np.zeros(len(free_gene), int)]))
    rows["orf"] = rows.gene.map(name2orf)
    return rows.dropna(subset=["orf"]).merge(g[["orf", "cds"]], on="orf", how="left"), len(chn)


def main():
    rows, n_on = model()
    by = rows.groupby("orf").agg(n=("ribosomes", "size"), occ=("ribosomes", lambda r: (r > 0).mean()),
                                 rib=("ribosomes", "mean"), cds=("cds", "first"))
    tr = rows[rows.ribosomes > 0].groupby("orf").agg(rib_tr=("ribosomes", "mean"))
    by = by.join(tr); by["dens"] = 100 * by.rib_tr / by.cds
    a = arava().set_index("orf"); j = by[by.n >= MIN_MRNA].join(a, how="inner")
    res = dict(mrnas=int(len(rows)), free_share=float((rows.ribosomes == 0).mean()),
               single_ribosome_share=float((rows.ribosomes == 1).mean()), single_ribosome_share_arava=0.14,
               mean_occupancy_over_genes=float(by.occ.mean()), mean_occupancy_arava=0.71,
               occupancy_spearman_vs_arava=dict(rho=float(spearmanr(j.occ, j.occupancy)[0]), n=int(len(j))),
               density_spearman_vs_arava=dict(rho=float(spearmanr(j.dens, j.density)[0]), n=int(j.dens.notna().sum())),
               density_mean_per_100nt=dict(model=float(by.dens.mean()), arava=float(a.density.mean())))
    if os.path.exists(RIBA):
        r = pd.read_csv(RIBA, sep="\t"); r["orf"] = r.ensembl_gene_id.str.upper()
        r["rpm"] = r.real_density * r.orf_length
        r = r.set_index("orf"); k = by[by.n >= MIN_MRNA].join(r[["rpm", "real_tpm"]], how="inner").dropna()
        res["ribosomes_per_mrna_vs_riba2019"] = dict(rho=float(spearmanr(k.rib, k.rpm)[0]), n=int(len(k)),
                                                     copy_weighted_mean_model=float(np.average(k.rib, weights=k.n)),
                                                     copy_weighted_mean_riba=float(np.average(r.rpm.dropna(), weights=r.real_tpm[r.rpm.notna()])))
    import inventory
    res["mrnas_per_average_cell"] = dict(model=float(len(rows) * inventory.POP_MEAN_VOL / inventory.CELL_VOL),
                                         miura2008=36139, zenklusen2008_about=60000, mode=mg.MRNA_MODE)
    st = f"{WORK}/validation/mrna_stats.json"
    if os.path.exists(st):
        res["fraction_cytosolic_80S_on_mrna"] = json.load(open(st)).get("fraction_cytosolic_80S_on_mrna")
        res["fraction_in_polysomes_measured"] = "0.82-0.85 (A254 profiles: Arava 2003 YPD 85 +- 5%, Metzl-Raz 2017 SC 1 - 0.18)"
    os.makedirs(f"{WORK}/validation", exist_ok=True)
    json.dump(res, open(f"{WORK}/validation/mrna_occupancy.json", "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
