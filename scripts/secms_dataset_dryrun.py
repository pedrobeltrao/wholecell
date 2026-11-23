"""Dry run of a different SEC-MS dataset feeding the cell's assembly rules (no repack, no rule sets).

    PYTHONHASHSEED=0 python scripts/secms_dataset_dryrun.py --cell jurkat --reference <build dir> \
        --secms-dir <dir with secms_*.tsv> --secms-select source_id=X,condition=Y,replicate=Z --out <dry-run dir>

Copies the reference build's inventory/{rna,v0} and the SEC-MS tables of --secms-dir (as inventory/secms), re-runs the
steps `assembly` and `ingredients` with assembly_steps.overrides.secms = the selection (everything else as in the cell
config: Jurkat = preset legacy_v7, complex_copies stoichiometry_cap), and compares with the reference build. Only the
SEC-fed rules can change: homo-oligomer share (a_k), AFDB homodimer share (a2), heterodimer SEC evidence (pair r).
Writes <out>/dryrun.json, <out>/dryrun_proteins.tsv (per protein: copies in each assembled category, ref vs new,
mass moved) and <out>/dryrun_species.tsv.
"""
import argparse
import json
import os
import shutil
import sys

os.environ.setdefault("PYTHONHASHSEED", "0")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.config import load  # noqa: E402
from wholecell.inventory import stage  # noqa: E402

AV = 6.02214076e23
CATS = {"homo_oligomers_v0.tsv": "homo_oligomers", "afdb_homodimers_v0.tsv": "afdb_homodimers",
        "afdb_heterodimers_v0.tsv": "afdb_heterodimers"}


def _rd(base, rel, **k):
    return pd.read_csv(os.path.join(base, rel), sep="\t", **k)


def per_protein_copies(base):
    """copies of each protein (uniprot) placed in each assembled category"""
    inv = _rd(base, "assembly/inventory_v0_complexed.tsv")
    out = pd.DataFrame({"gene": inv.gene.values, "mass_da": inv.mass_da.values, "copies": inv.copies.values,
                        "complex_portal": inv.copies_in_complexes.fillna(0).values}, index=inv.uniprot.values)
    out = out[~out.index.duplicated()]
    g2u = dict(zip(inv.gene, inv.uniprot))
    H = _rd(base, "assembly/homo_oligomers_v0.tsv")
    D = _rd(base, "assembly/afdb_homodimers_v0.tsv")
    X = _rd(base, "assembly/afdb_heterodimers_v0.tsv")
    het = pd.concat([pd.Series(X.n.values, index=X.gene1.map(g2u).values), pd.Series(X.n.values, index=X.gene2.map(g2u).values)])
    for cat, ser in (("homo_oligomers", H.groupby("uniprot").copies_used.sum()),
                     ("afdb_homodimers", D.groupby("uniprot").n_dimers.sum() * 2),
                     ("afdb_heterodimers", het.groupby(level=0).sum())):
        out[cat] = ser.reindex(out.index).fillna(0).values
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--secms-dir", required=True)
    ap.add_argument("--secms-select", required=True, help="source_id=..,condition=..,replicate=..")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=2)
    a = ap.parse_args()
    sel = dict(kv.split("=", 1) for kv in a.secms_select.split(","))
    cfg0 = load(a.cell)
    name = cfg0["inventory"]["adapters"]["assembly"]
    if name != "assembly_steps":
        sys.exit(f"assembly adapter {name}: this script handles assembly_steps (overrides.secms)")
    cfg = load(a.cell, overrides={"inventory": {"params": {name: {"overrides": {"secms": sel}, "workers": a.workers}}}})
    registry.load_plugins(cfg)
    ref = os.path.join(a.reference, "inventory")
    out = os.path.join(a.out, "inventory")
    if os.path.exists(out):
        sys.exit(f"{out} exists (write-once)")
    for sub in ("rna", "v0"):
        if os.path.isdir(os.path.join(ref, sub)):
            shutil.copytree(os.path.join(ref, sub), os.path.join(out, sub))
    os.makedirs(os.path.join(out, "secms"))
    for f in os.listdir(a.secms_dir):
        if f.endswith(".tsv") or f.endswith(".json"):
            shutil.copy(os.path.join(a.secms_dir, f), os.path.join(out, "secms", f))
    stage.run(cfg, a.out, step="assembly")
    stage.run(cfg, a.out, step="ingredients")

    res = {"cell": a.cell, "reference": a.reference, "secms_dir": a.secms_dir, "secms_select": sel}
    tot = float(_rd(ref, "assembly/inventory_v0_complexed.tsv").mass_pg.sum())
    cx = {}
    for lab, base in (("ref", ref), ("new", out)):
        d = {"complex_portal": float(_rd(base, "assembly/complexes_v0.tsv").mass_built_pg.sum())}
        for f, cat in CATS.items():
            d[cat] = float(_rd(base, f"assembly/{f}").mass_pg.sum())
        d["complexed_total"] = sum(d.values())
        d["complexed_share"] = d["complexed_total"] / tot
        cx[lab] = d
    res["mass_pg"] = cx
    res["total_protein_pg"] = tot
    Pr, Pn = per_protein_copies(ref), per_protein_copies(out)
    cats = ["complex_portal"] + list(CATS.values())
    P = Pr[["gene", "mass_da", "copies"]].copy()
    for c in cats:
        P[f"{c}_ref"] = Pr[c]
        P[f"{c}_new"] = Pn[c].reindex(P.index).fillna(0)
    P["assembled_ref"] = Pr[cats].sum(1)
    P["assembled_new"] = Pn[cats].sum(1).reindex(P.index).fillna(0)
    P["delta_copies"] = P.assembled_new - P.assembled_ref
    P["delta_mass_pg"] = P.delta_copies * P.mass_da.fillna(0) / AV * 1e12
    P["share_ref"] = P.assembled_ref / P.copies.replace(0, np.nan)
    P["share_new"] = P.assembled_new / P.copies.replace(0, np.nan)
    ch = P[P.delta_copies.abs() > 0.5]
    ch.reindex(ch.delta_mass_pg.abs().sort_values(ascending=False).index).to_csv(
        os.path.join(a.out, "dryrun_proteins.tsv"), sep="\t", index_label="uniprot")
    res["proteins"] = dict(n_changed=int(len(ch)), n_more_assembled=int((ch.delta_copies > 0).sum()),
                           n_less_assembled=int((ch.delta_copies < 0).sum()),
                           mass_to_assembled_pg=float(ch.delta_mass_pg.clip(lower=0).sum()),
                           mass_to_free_pg=float(-ch.delta_mass_pg.clip(upper=0).sum()))
    res["proteins"]["mass_moved_pg"] = res["proteins"]["mass_to_assembled_pg"] + res["proteins"]["mass_to_free_pg"]
    res["proteins"]["mass_moved_pct_total"] = 100 * res["proteins"]["mass_moved_pg"] / tot
    Sr, Sn = _rd(ref, "ingredients/ingredients.tsv"), _rd(out, "ingredients/ingredients.tsv")
    s = Sr[["ing_id", "copies", "mass_pg"]].merge(Sn[["ing_id", "copies", "mass_pg"]], on="ing_id", how="outer",
                                                   suffixes=("_ref", "_new")).fillna(0)
    s["d"] = s.mass_pg_new - s.mass_pg_ref
    s[s.copies_ref != s.copies_new].sort_values("d").to_csv(os.path.join(a.out, "dryrun_species.tsv"), sep="\t", index=False)
    res["species"] = dict(n_ref=int(len(Sr)), n_new=int(len(Sn)), n_added=int((s.copies_ref == 0).sum()),
                          n_gone=int((s.copies_new == 0).sum()), n_copies_changed=int((s.copies_ref != s.copies_new).sum()),
                          mass_moved_pg=float(s.d.abs().sum() / 2))
    top = ch.reindex(ch.delta_mass_pg.abs().sort_values(ascending=False).index).head(15)
    res["top_proteins"] = [dict(uniprot=u, gene=r.gene, delta_mass_pg=float(r.delta_mass_pg),
                                share_ref=float(r.share_ref), share_new=float(r.share_new)) for u, r in top.iterrows()]
    json.dump(res, open(os.path.join(a.out, "dryrun.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "top_proteins"}, indent=1, default=str))


if __name__ == "__main__":
    main()
