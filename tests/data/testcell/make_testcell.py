"""Generator of the synthetic testcell inventory sources (run once; outputs are committed):
  proteome.tsv       per gene: gene, uniprot, mass_da, copies              (adapter.proteome   proteome_table)
  localisation.tsv   per gene: gene, compartments, transmembrane, n_tm, signal_peptide, mito_om (localisation_table)
  complexes.tsv      complex_id, subunits, copies, struct_src              (adapter.assembly   complex_table)
  library_ids.txt    species ids the packing agent copies from the legacy library ({data_root}/jurkat/yc/pipeline/lib)

Species ids, subunit composition, gene masses / accessions are taken from the legacy Jurkat tables (ingredients_v7.tsv,
inventory_v0.tsv; read-only); copies are synthetic: soluble species fill their compartment volume (analytic testcell,
toolkits/analytic/specs/testcell.yaml: cytosol ~5.1, nucleoplasm ~0.53, matrix ~0.035, ER lumen ~0.010 um3) to the
target density by mass fraction; membrane species get a few hundred to a few thousand copies.

    PYTHONHASHSEED=0 python tests/data/testcell/make_testcell.py
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
from wholecell.config import load_site  # noqa: E402

AV = 6.02214076e23
VOL = {"cytosol": 5.1, "nucleus": 0.53, "mitochondrion": 0.035, "er": 0.010}          # um3 (approximate)
DENS = {"cytosol": 85.0, "nucleus": 90.0, "mitochondrion": 90.0, "er": 90.0}          # mg/ml targets
SOLUBLE = {   # compartment -> [(species id, mass fraction, compartments string)]
    "cytosol": [("CPX-5993", .16, None), ("CPX-6030", .14, None), ("GAPDH_x4", .12, None), ("ALDOA_x4", .07, None),
                ("ENO1_x2", .08, None), ("HSP90AB1", .08, None), ("PGK1", .04, "cytosol:0.9;er:0.1"),   # no signal peptide:
                ("EEF2", .07, None), ("YWHAZ_afx2", .04, None), ("CPX-3288", .06, None),             # generic.signal_lumen
                ("ACTB", .03, None), ("PFN1", .01, None), ("PPIA", .01, None),                      # moves its er share
                ("LCP1", .04, "cortex:0.75;cytosol:0.25")],
    "nucleus": [("CPX-2940", .20, None), ("HNRNPA1", .20, None), ("HMGB1", .15, None), ("NPM1", .15, None),
                ("PARP1", .15, None), ("HNRNPC", .15, None)],
    "mitochondrion": [("HSPD1_x7", .30, None), ("CS_x2", .20, None), ("GOT2", .20, None), ("MDH2", .15, None), ("ACO2", .15, None)],
    "er": [("HSPA5", .40, None), ("CALR", .25, None), ("PDIA3", .20, None), ("HSP90B1_afx2", .15, None)],
}
MEMBRANE = [  # (species id, copies, compartments)
    ("PTPRC", 3000, "plasma_membrane:1.0"), ("CD47", 2500, "plasma_membrane:1.0"), ("SLC2A1", 2000, "plasma_membrane:1.0"),
    ("ATP1A1", 2500, "plasma_membrane:1.0"), ("CD5", 800, "plasma_membrane:1.0"),
    ("CANX", 1500, "er:1.0"), ("RPN1", 800, "er:1.0"), ("SSR1", 600, "er:1.0"),
    ("CPX-6151", 300, "mitochondrion:1.0"), ("SLC25A5", 2000, "mitochondrion:1.0"), ("SLC25A3", 1000, "mitochondrion:1.0"),
    ("VDAC1", 2500, "mitochondrion:1.0"), ("TOMM40", 500, "mitochondrion:1.0"),
    ("LBR", 1000, "nuclear_envelope:1.0"), ("EMD", 500, "nuclear_envelope:1.0"),
    ("LCK", 5000, "pm_anchor_cyto:1.0"),                                    # lipid-anchored peripheral shell species
]
OM = {"VDAC1", "TOMM40"}


def main():
    root = load_site()["data_root"]
    I = pd.read_csv(f"{root}/data/derived/ingredients_v7.tsv", sep="\t").set_index("ing_id")
    V0 = pd.read_csv(f"{root}/data/derived/inventory_v0.tsv", sep="\t").set_index("gene")
    up = pd.read_csv(f"{root}/data/reference/uniprot_human_full.tsv.gz", sep="\t", usecols=["Entry", "Signal peptide"]).set_index("Entry")
    lib = f"{root}/jurkat/yc/pipeline/lib"
    copies, comp, cpx = {}, {}, []

    def add_species(sid, n, cstr):
        r = I.loc[sid]
        gs, st = r.genes.split(";"), [int(float(x)) for x in str(r.stoich).split(";")]
        if len(gs) > 1 or r.kind != "protein":
            cpx.append(dict(complex_id=sid, subunits=";".join(f"{g}:{s}" for g, s in zip(gs, st)), copies=int(n),
                            struct_src=r.struct_src))
        for g, s in zip(gs, st):
            copies[g] = copies.get(g, 0) + int(n) * s
            comp.setdefault(g, cstr)

    for c, lst in SOLUBLE.items():
        tot = sum(f for _, f, _ in lst)
        pg = VOL[c] * DENS[c] * 1e-3
        for sid, f, cstr in lst:
            n = int(round(pg * f / tot * 1e-12 / (I.at[sid, "mass_da"] / AV)))
            add_species(sid, n, cstr or (c + ":1.0"))
    for sid, n, cstr in MEMBRANE:
        add_species(sid, n, cstr)
    genes = sorted(copies)
    P = pd.DataFrame(dict(gene=genes, uniprot=[V0.at[g, "uniprot"] for g in genes], mass_da=[float(V0.at[g, "mass_da"]) for g in genes],
                          copies=[copies[g] for g in genes]))
    tm = [bool(V0.at[g, "transmembrane"]) for g in genes]
    L = pd.DataFrame(dict(gene=genes, compartments=[comp[g] for g in genes], transmembrane=tm, n_tm=[int(t) for t in tm],
                          signal_peptide=[isinstance(up["Signal peptide"].get(V0.at[g, "uniprot"]), str) for g in genes],
                          mito_om=[g in OM for g in genes]))
    C = pd.DataFrame(cpx)
    hdr = ("# testcell synthetic {what} (generated by tests/data/testcell/make_testcell.py; species ids / subunits / masses from the\n"
           "# legacy Jurkat tables, copies synthetic: ~85-90 mg/ml in the analytic testcell volumes, membranes 300-5000 copies)\n")
    for df, name, what in ((P, "proteome.tsv", "proteome source"), (L, "localisation.tsv", "localisation source"),
                           (C, "complexes.tsv", "complex source")):
        with open(os.path.join(HERE, name), "w") as f:
            f.write(hdr.format(what=what)); df.to_csv(f, sep="\t", index=False)
    free_genes = {g for g in genes if g in {s for s, _, _ in sum(SOLUBLE.values(), [])} | {s for s, _, _ in MEMBRANE}}
    ids = sorted({c["complex_id"] for c in cpx} | free_genes)
    missing = [i for i in ids if not os.path.exists(f"{lib}/{i}.npz")]
    assert not missing, missing
    open(os.path.join(HERE, "library_ids.txt"), "w").write(
        "# species ids of the testcell inventory = file names <id>.npz in the legacy library {data_root}/jurkat/yc/pipeline/lib\n"
        + "\n".join(ids) + "\n")
    mass = (P.copies * P.mass_da).sum() / AV * 1e12
    print(f"{len(genes)} genes, {len(ids)} species, {P.copies.sum():,} gene copies, total {mass:.5f} pg")


if __name__ == "__main__":
    main()
