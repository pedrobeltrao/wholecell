# origin: jurkat/yc/pipeline/jingredients.py (legacy workspace, 2026-11-15)
# changes: inventory_builder "jurkat.ingredients_v0"; reads the build's inventory/assembly/* instead of data/derived;
# absolute paths removed (frozen inputs); RIBO80 / NPC pdb ids, ribosome subunit complex ids, RIBO_W and the histone
# families come from plugins/jurkat/data/ingredients_v0.yaml + params histone_families (TSV, since 2026-11-16 A4). Returns the table (the stage writes it).
# Numerics identical to legacy.
"""Jurkat: inventory + assembly tables -> packing ingredients (pack-table columns).
  - complexes kept as one species if a local PDB assembly covers >= 50% of subunit mass, else returned to the pools;
    60S + 40S merged into 80S (6QZP) at min(60S, 40S), remainder as free subunits; NPC uses 7R5K.
  - AFDB dimers only when the model file is local (frozen fetch status ok/cached).
  - comp_weights of an assembly = mass-weighted mean of its subunits; transmembrane if any subunit is TM.
  - free monomers = inventory copies - copies consumed by all assemblies (clipped at 0).
  - signal_frac / mito_om_frac = subunit-mass fraction with a UniProt signal peptide / 'Mitochondrion outer membrane'."""
import json
import os

import numpy as np
import pandas as pd
import yaml

from wholecell.registry import register

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def parse_w(s):
    if not isinstance(s, str) or not s:
        return {"cytosol": 1.0}
    out = {}
    for t in s.split(";"):
        k, v = t.split(":"); out[k] = float(v)
    return out


def ingredients_v0(ctx):
    p = ctx.params("jurkat.ingredients_v0")
    cur = yaml.safe_load(open(os.path.join(DATA, "ingredients_v0.yaml")))
    from wholecell.inventory.assembly.steps import read_curated
    H = read_curated(p["histone_families"])
    fam = dict(zip(H.family, H.gene_regex))
    RIBO80, NPC_PDB = cur["ribosome_80s_pdb"], cur["npc_pdb"]
    C60, C40 = cur["ribosome_60s_complex"], cur["ribosome_40s_complex"]
    NPC_ID = cur["npc_complex"]
    inv = pd.read_csv(ctx.need("assembly/inventory_v0_complexed.tsv"), sep="\t", index_col="gene")
    inv["w"] = inv.compartments.map(parse_w)
    up = pd.read_csv(ctx.frozen(p["uniprot_full"]), sep="\t",
                     usecols=["Entry", "Signal peptide", "Subcellular location [CC]"]).set_index("Entry")
    sig = inv.uniprot.map(lambda a: isinstance(up["Signal peptide"].get(a), str)).fillna(False).astype(bool)
    om = inv.uniprot.map(lambda a: "Mitochondrion outer membrane" in str(up["Subcellular location [CC]"].get(a, ""))).fillna(False).astype(bool)
    tm = inv.transmembrane.fillna(False).astype(bool)
    mass = inv.mass_da.fillna(0.0)
    used = pd.Series(0.0, index=inv.index)
    rows = []

    def asm_weights(sub):
        tot = sum(mass[g] * s for g, s in sub.items()) or 1.0
        w = {}
        for g, s in sub.items():
            for k, v in inv.w[g].items():
                w[k] = w.get(k, 0.0) + v * mass[g] * s / tot
        return {k: round(v, 3) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 0.01}

    def mfrac(sub, flag):
        tot = sum(mass[g] * s for g, s in sub.items()) or 1.0
        return round(sum(mass[g] * s for g, s in sub.items() if flag[g]) / tot, 3)

    def add(ing_id, kind, sub, copies, src):
        if copies <= 0:
            return
        for g, s in sub.items():
            used[g] += copies * s
        rows.append(dict(ing_id=ing_id, kind=kind, genes=";".join(sub), stoich=";".join(str(s) for s in sub.values()),
                         copies=int(copies), mass_da=float(sum(mass[g] * s for g, s in sub.items())),
                         comp_weights=json.dumps(asm_weights(sub)), transmembrane=bool(any(tm[g] for g in sub)),
                         n_tm=int(sum(tm[g] for g in sub)), struct_src=src, uniprot="",
                         signal_frac=mfrac(sub, sig), mito_om_frac=mfrac(sub, om)))

    C = pd.read_csv(ctx.need("assembly/complexes_v0.tsv"), sep="\t")
    C = C[C.copies_built > 0].set_index("cpx")
    nuc = cur["nucleosome"]
    for cpx, r in C.iterrows():
        if cpx == nuc["complex"]:
            n = int(r.copies_built)
            for f, pat in fam.items():
                gs = inv.index[inv.index.str.match(pat)]
                tot = inv.copies[gs].sum()
                for g in gs:
                    used[g] += 2 * n * inv.copies[g] / tot
            rows.append(dict(ing_id=nuc["ing_id"], kind="complex", genes="H2A;H2B;H3;H4", stoich="2;2;2;2", copies=n,
                             mass_da=float(nuc["mass_da"]), comp_weights=json.dumps(nuc["comp_weights"]), transmembrane=False, n_tm=0,
                             struct_src=nuc["struct_src"], uniprot="", signal_frac=0.0, mito_om_frac=0.0))
            continue
        if cpx in (C60, C40):
            continue
        sub = {x.split(":")[0]: int(x.split(":")[1]) for x in r.subunits.split(";")}
        sub = {g: s for g, s in sub.items() if g in inv.index and inv.copies[g] > 0}
        if cpx == NPC_ID:
            add(cur["npc_ing_id"], "complex", sub, int(r.copies_built), f"pdb:{NPC_PDB}")
            continue
        if isinstance(r.best_pdb, str) and r.pdb_mass_cov >= p.get("pdb_mass_cov_min", 0.5):
            add(cpx, "complex", sub, int(r.copies_built), f"pdb:{r.best_pdb}")
    s60 = {x.split(":")[0]: int(x.split(":")[1]) for x in C.loc[C60, "subunits"].split(";")}
    s40 = {x.split(":")[0]: int(x.split(":")[1]) for x in C.loc[C40, "subunits"].split(";")}
    n60, n40 = int(C.loc[C60, "copies_built"]), int(C.loc[C40, "copies_built"])
    n80 = min(n60, n40)
    sub80 = {g: s for g, s in {**s60, **s40}.items() if g in inv.index}
    add(cur["ribosome_80s_ing_id"], "complex", sub80, n80, f"pdb:{RIBO80}")
    add(C60, "complex", {g: s for g, s in s60.items() if g in inv.index}, n60 - n80, f"pdb:{C.loc[C60, 'best_pdb']}")
    add(C40, "complex", {g: s for g, s in s40.items() if g in inv.index}, n40 - n80, f"pdb:{C.loc[C40, 'best_pdb']}")
    H = pd.read_csv(ctx.need("assembly/homo_oligomers_v0.tsv"), sep="\t")
    for r in H.itertuples():
        add(f"{r.gene}_x{r.state}", "oligomer", {r.gene: int(r.state)}, int(r.n_oligomers), f"pdb:{r.example_pdb}")
    st = pd.read_csv(ctx.frozen(p["afdb_fetch_status"]), sep="\t").set_index("model").status
    okm = set(st.index[st.isin(["ok", "cached"])])
    for r in pd.read_csv(ctx.need("assembly/afdb_homodimers_v0.tsv"), sep="\t").itertuples():
        if r.model in okm:
            add(f"{r.gene}_afx2", "afdb_homodimer", {r.gene: 2}, int(r.n_dimers), f"afdb_complex:{r.model}")
    for r in pd.read_csv(ctx.need("assembly/afdb_heterodimers_v0.tsv"), sep="\t").itertuples():
        if r.model in okm and r.gene1 != r.gene2:
            add(f"{r.gene1}-{r.gene2}", "afdb_heterodimer", {r.gene1: 1, r.gene2: 1}, int(r.n), f"afdb_complex:{r.model}")
    free = (inv.copies - used).clip(lower=0).round().astype(np.int64)
    over = (used - inv.copies)[used > inv.copies * 1.02 + 10]
    for g in inv.index[free > 0]:
        rows.append(dict(ing_id=g, kind="protein", genes=g, stoich="1", copies=int(free[g]), mass_da=float(mass[g]),
                         comp_weights=json.dumps(inv.w[g]), transmembrane=bool(tm[g]), n_tm=int(tm[g]),
                         struct_src=f"afdb_monomer:{inv.uniprot[g]}", uniprot=inv.uniprot[g],
                         signal_frac=float(sig[g]), mito_om_frac=float(om[g])))
    I = pd.DataFrame(rows)
    for k, w in cur["ribosome_weights"].items():
        I.loc[I.ing_id == k, "comp_weights"] = json.dumps(w)
    I["mass_pg"] = I.copies * I.mass_da / 6.022e23 * 1e12
    ctx.say(f"ingredients: {len(I)} species, {I.copies.sum():.3g} copies, {I.mass_pg.sum():.1f} pg; genes overdrawn (>2%): {len(over)}")
    return I


ingredients_v0.PREV = ["assembly/inventory_v0_complexed.tsv", "assembly/complexes_v0.tsv", "assembly/homo_oligomers_v0.tsv",
                       "assembly/afdb_homodimers_v0.tsv", "assembly/afdb_heterodimers_v0.tsv"]
register("inventory_builder", "jurkat.ingredients_v0")(ingredients_v0)
