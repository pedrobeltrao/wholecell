# origin: yeast yeast_cell_dev_bundle/validation/membrane_proteomics.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Membrane localization against organelle-membrane proteomics (2026-11-16): held-out data.

Lists of proteins found in purified or enriched membrane fractions (pipeline/data/validation/membrane_proteomics/,
README.md there): plasma membrane (Villers 2017 enriched; Delom 2006, Navarre 2002 present), vacuolar membrane
(Wiederhold 2009 enriched), mitochondrial outer membrane (Zahedi 2006), Golgi (Inadome 2005), tubular ER (Wang 2017),
peroxisome (Schafer 2001, Yi 2002, Marelli 2004), lipid droplets (Athenstaedt 1999, Currie 2014). None comes from the
GFP collection or CYCLoPs, the model's localization inputs (the TSVs carry SGD-2009 and Huh-2003 columns; they are not
read here). Membrane fractions carry contaminants (the README estimates them per list), so recall is the robust number;
AUROC asks whether the model's weight for that compartment ranks the listed proteins above all other expressed ones.
Model side: the compartment weights of each protein (features table of the build: comp_weights, shares of its copies).
Output: work/validation/membrane_proteomics.json
"""
import os, sys, json
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import DATA, WORK  # noqa: E402

D = f"{DATA}/validation/membrane_proteomics"
SETS = [  # file, model compartment(s), kind of list
    ("villers_2017/villers2017_plasma_membrane.tsv", ("plasma_membrane",), "enriched"),
    ("wiederhold_2010_review/delom2006_plasma_membrane.tsv", ("plasma_membrane",), "present"),
    ("wiederhold_2010_review/navarre2002_plasma_membrane.tsv", ("plasma_membrane",), "present"),
    ("wiederhold_2010_review/wiederhold2009_vacuolar_membrane.tsv", ("vacuole",), "enriched"),
    ("wiederhold_2010_review/zahedi2006_mito_outer_membrane.tsv", ("mito_OM",), "present"),
    ("wiederhold_2010_review/inadome2005_golgi.tsv", ("golgi_endosome",), "present"),
    ("wang_2017/wang2017_tubular_er.tsv", ("ER", "nuclear_envelope"), "enriched"),
    ("wiederhold_2010_review/schafer2001_peroxisome_membrane.tsv", ("peroxisome",), "present"),
    ("wiederhold_2010_review/yi2002_peroxisome.tsv", ("peroxisome",), "present"),
    ("wiederhold_2010_review/marelli2004_peroxisome_icat.tsv", ("peroxisome",), "enriched"),
    ("wiederhold_2010_review/athenstaedt1999_lipid_droplet.tsv", ("lipid_droplet",), "present"),
    ("currie_2014/currie2014_lipid_droplet.tsv", ("lipid_droplet",), "enriched"),
]
MITO = ("mito_OM", "mito_IM", "mito_IMS", "mito_matrix")


def auroc(pos, neg):
    s = np.r_[pos, neg]; r = pd.Series(s).rank().values
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))) if len(pos) and len(neg) else float("nan")


def main(inventory=False, out=True):
    if inventory:     # the inventory's weights, rebuilt from the current code (to compare rule variants without a build)
        import inventory as I
        f = I.build()[0]
    else:
        f = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    f = f[f.copies_per_cell > 0].copy()
    W = [json.loads(w) if isinstance(w, str) else (w or {}) for w in f.comp_weights]
    f["top"] = [max(w, key=w.get) if w else "cytosol" for w in W]
    res = {}
    for fn, comps, kind in SETS:
        d = pd.read_csv(f"{D}/{fn}", sep="\t"); orfs = set(d.orf.astype(str).str.upper())
        score = np.array([sum(w.get(c, 0.0) for c in comps) for w in W])
        is_pos = f.orf.astype(str).str.upper().isin(orfs).values
        top_ok = f.top.isin(comps).values
        name = os.path.basename(fn).replace(".tsv", "")
        miss = f[is_pos & ~top_ok].sort_values("copies_per_cell", ascending=False)
        res[name] = dict(compartment=list(comps), list_kind=kind, listed=len(orfs), in_model=int(is_pos.sum()),
                         recall_top=float(top_ok[is_pos].mean()) if is_pos.any() else None,
                         recall_weight_0p2=float((score[is_pos] >= 0.2).mean()) if is_pos.any() else None,
                         auroc=auroc(score[is_pos], score[~is_pos]),
                         model_top_of_listed=f.top[is_pos].value_counts().head(6).to_dict(),
                         largest_misses=[(g, t, int(c)) for g, t, c in zip(miss.gene.head(12), miss.top.head(12), miss.copies_per_cell.head(12))])
        r = res[name]
        print(f"{name:42s} {kind:8s} n {r['in_model']:4d}  recall top {r['recall_top']:.2f}  w>=0.2 {r['recall_weight_0p2']:.2f}  AUROC {r['auroc']:.2f}")
    if out:
        os.makedirs(f"{WORK}/validation", exist_ok=True)
        json.dump(res, open(f"{WORK}/validation/membrane_proteomics.json", "w"), indent=1)
    return res


if __name__ == "__main__":
    main(inventory="--inventory" in sys.argv)
