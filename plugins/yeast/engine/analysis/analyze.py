# origin: yeast yeast_cell_dev_bundle/analysis/analyze.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Accuracy gap analysis of the yeast cell model vs. available structures/annotations."""
import sys, json, re, collections
import numpy as np, pandas as pd

import os  # PORT PATCH 5 (release 2026-11-23): author-machine paths -> env (see PORT_PATCHES.md)
B = os.environ.get("YC_BUNDLE", "yeast_cell_dev_bundle")
V = os.environ.get("YC_VIEWER_TYPES", "data/library/types.json")
DIM = os.environ.get("YC_AFDB_DIMERS", "afdb_yeast_dimers")
sys.path.insert(0, f"{B}/pipeline")
import inventory
inventory.D = f"{B}/pipeline/data"

df, g = inventory.build()
types = json.load(open(V))["types"]
T = pd.DataFrame(types)
T["mc"] = T.mass * T.copies
TOTM = T.mc.sum()
af = T[T.src == "AF"].copy()
acc2gene = dict(zip(df.uniprot, df.gene))
cop = dict(zip(df.gene, df.copies_per_cell))
ntm = dict(zip(df.gene, df.n_tm))
mass = dict(zip(df.gene, df.mass_da))
olig = dict(zip(df.gene, df.oligomer_n))
af_by_acc = {r.uniprot: r for r in af.itertuples()}
pdb_genes = set(g_ for t in T[T.src == "PDB"].itertuples() for g_ in re.split(r"[;\s]+", t.genes) if g_)

def pct(x): return f"{100 * x:.1f}%"
print("== overall")
print("types", len(T), "AF", len(af), "copies", T.copies.sum(), "AF mass share", pct(af.mc.sum() / TOTM))

print("\n== 1. oligomeric state (AF monomers)")
af["uni_olig"] = af.gene.map(olig).fillna(1)
o = af[af.uni_olig > 1]
print(f"AF types UniProt says homo-oligomer: {len(o)} types, {int(o.copies.sum())} copies, {pct(o.mc.sum() / TOTM)} of cell protein mass")
print(o.groupby("uni_olig").agg(n=("id", "size"), copies=("copies", "sum")).to_string())
hom = pd.read_csv(f"{DIM}/homodimers/candidates.csv")
hacc = set(hom.uniprotAccession_1)
h_af = af[af.uniprot.isin(hacc)]
print(f"homodimer models: {len(hom)}; matching placed AF monomer types: {len(h_af)}, copies {int(h_af.copies.sum())}, mass {pct(h_af.mc.sum() / TOTM)}")
print("  of these, UniProt says monomer/unknown:", int((h_af.uni_olig == 1).sum()), " dimer:", int((h_af.uni_olig == 2).sum()), " higher:", int((h_af.uni_olig > 2).sum()))
print("  ipTM range", hom.ipTM.min(), "-", hom.ipTM.max(), " top by copies:")
print(h_af.sort_values("copies", ascending=False)[["gene", "copies", "uni_olig", "plddt_mean"]].head(15).to_string(index=False))
o_nomodel = o[~o.uniprot.isin(hacc)]
print(f"UniProt oligomers WITHOUT a homodimer model: {len(o_nomodel)} types, mass {pct(o_nomodel.mc.sum() / TOTM)}")

print("\n== 2. heterodimers")
het = pd.read_csv(f"{DIM}/heterodimers/candidates.csv")
rows = []
for r in het.itertuples():
    a, b = af_by_acc.get(r.uniprotAccession_1), af_by_acc.get(r.uniprotAccession_2)
    inA = acc2gene.get(r.uniprotAccession_1) in pdb_genes or acc2gene.get(r.uniprotAccession_2) in pdb_genes
    rows.append(dict(g1=r.gene_1, g2=r.gene_2, ipTM=r.ipTM, pDockQ2=r.max_pDockQ2, both_af=a is not None and b is not None, in_pdb_cpx=inA,
                     c1=a.copies if a is not None else 0, c2=b.copies if b is not None else 0,
                     m=(a.mass + b.mass) if (a is not None and b is not None) else 0))
H = pd.DataFrame(rows)
H["pair_copies"] = np.minimum(H.c1, H.c2)
H["mc"] = H.pair_copies * H.m
print(f"heterodimer models {len(H)}; both partners placed as AF monomers: {int(H.both_af.sum())}; touching an existing PDB complex: {int(H.in_pdb_cpx.sum())}")
print(f"  upper bound mass that could be paired: {pct(H.mc.sum() / TOTM)} (sum over pairs, overlaps not resolved)")
print(f"  ipTM>=0.8: {int((H.ipTM >= 0.8).sum())}, 0.6-0.8: {int(((H.ipTM >= 0.6) & (H.ipTM < 0.8)).sum())}")
deg = collections.Counter(list(H.g1) + list(H.g2))
print("  proteins in >1 pair (hubs / larger complexes):", sum(1 for v in deg.values() if v > 1), " top:", deg.most_common(8))
print(H.sort_values("mc", ascending=False).head(15)[["g1", "g2", "ipTM", "pDockQ2", "pair_copies"]].to_string(index=False))

print("\n== 3. Complex Portal coverage")
cp = pd.read_csv(f"{B}/pipeline/data/complexportal_559292.tsv", sep="\t")
crow = []
for r in cp.itertuples(index=False):
    sub = {}
    for a, s in re.findall(r"([A-Z0-9]{6,10})(?:-\w+)?\((\d+)\)", r[4]):
        gg = acc2gene.get(a)
        if gg: sub[gg] = max(1, int(s))
    if len(sub) < 2: continue
    c = float(np.median([cop.get(x, 0) / s for x, s in sub.items()]))
    m = sum(mass.get(x, 0) * s for x, s in sub.items())
    modeled = len(set(sub) & pdb_genes) / len(sub)
    crow.append(dict(cpx=r[0], name=r[1][:45], nsub=len(sub), copies=c, mass_kda=m / 1e3, mc=c * m,
                     tm=any(ntm.get(x, 0) > 0 for x in sub), modeled=modeled,
                     pdbs=len(re.findall(r"wwpdb:(\w{4})", str(r[8])))))
C = pd.DataFrame(crow).sort_values("mc", ascending=False)
print(f"Complex Portal multi-subunit complexes: {len(C)}; fully modeled as assemblies: {int((C.modeled == 1).sum())}; partly: {int(((C.modeled > 0) & (C.modeled < 1)).sum())}")
um = C[C.modeled == 0]
print(f"unmodeled: {len(um)} ({int(um.tm.sum())} membrane); their subunit mass ~{pct(um.mc.sum() / TOTM)} of cell protein (overlapping complexes double-count)")
print(um.head(45)[["cpx", "name", "nsub", "copies", "mass_kda", "tm", "pdbs"]].to_string(index=False))
C.to_csv(f"{sys.argv[1]}/complex_coverage.csv", index=False)

print("\n== 4. membrane proteins")
mem_ids = {str(i) for i in range(14, 22)}
names = ["PM", "ER", "NE", "VAC", "MITO_OM", "MITO_IM", "GOLGI", "PEROX"]
mt = collections.Counter()
for t in types:
    for k, v in t["regions"].items():
        if k in mem_ids: mt[names[int(k) - 14]] += v
print("placed membrane-embedded copies:", dict(mt))
tmdf = df[df.n_tm > 0]
print(f"TM proteins in inventory: {len(tmdf)} genes, {int(tmdf.copies_per_cell.sum())} copies")
locs = collections.Counter()
for r in tmdf.itertuples():
    for c_, w in json.loads(r.comp_weights).items(): locs[c_] += w * r.copies_per_cell
print("TM copies by annotated compartment:", {k: int(v) for k, v in locs.most_common()})
print("TM types by oligomer annotation:", tmdf.oligomer_n.value_counts().to_dict())
# TM genes in Complex Portal membrane complexes that are unmodeled
memc = um[um.tm]
print(f"unmodeled membrane complexes: {len(memc)}; top by mass:")
print(memc.head(20)[["name", "nsub", "copies", "mass_kda", "pdbs"]].to_string(index=False))

print("\n== 5. pLDDT trimming (mass not drawn)")
af["lost"] = (1 - af.plddt_frac_kept) * af.mc
print(f"AF residue mass trimmed: {pct(af.lost.sum() / TOTM)} of all protein mass; types with <50% kept: {int((af.plddt_frac_kept < 0.5).sum())}")
reg_names = ["cytosol", "nucleoplasm", "nucleolus", "ER_lumen", "vacuole", "mito_matrix", "mito_IMS", "golgi_endosome", "peroxisome", "cell_wall"]
lost = collections.Counter(); tot = collections.Counter()
for t in types:
    if t["src"] != "AF": continue
    n = sum(t["regions"].values()) or 1
    for k, v in t["regions"].items():
        nm = reg_names[int(k)] if int(k) < 10 else "other"
        lost[nm] += (1 - t["plddt_frac_kept"]) * t["mass"] * v; tot[nm] += t["mass"] * v
print({k: pct(lost[k] / tot[k]) for k in tot if tot[k] > 0})
print(af.sort_values("lost", ascending=False).head(10)[["gene", "copies", "plddt_frac_kept"]].to_string(index=False))

print("\n== 6. cytoskeleton & polymer genes (how they are placed now)")
for gname in ["ACT1", "TUB1", "TUB2", "TUB3", "ABP1", "CRN1", "ARC40", "ARP2", "ARP3", "CAP1", "TPM1", "TPM2", "COF1", "PFY1", "SAC6", "MYO1", "MYO2", "BNI1", "PIL1", "LSP1", "SUR7", "CDC3", "CDC10", "CDC11", "CDC12", "SHS1", "SPC42", "SPC110"]:
    hit = T[T.genes.str.contains(rf"(?:^|;){gname}(?:;|$)", regex=True)]
    print(f"  {gname:7s} copies(inv)={cop.get(gname, 0):9.0f}  as: {', '.join(f'{h.id}({h.src})' for h in hit.itertuples())}")

print("\n== 7. nucleic acids")
print(T[T.nnt > 0][["id", "pdb", "nnt", "copies"]].to_string(index=False))
