#!/usr/bin/env python
"""Yeast species table for the core library stage (A3): the species the yeast s0 library builds, from the vendored
engine's inventory and ingredients (inventory.build(), ingredients.build_ingredients(), exactly as s0_library.py), in
ingredient order. In a build the plugin step yeast_library_species (bridge.library_species) makes it from the build's
protein table (core copies + localisation) with species_table(); this CLI is the A3 check's standalone version.

    <yeast env python> plugins/yeast/library_species.py --out <build>/inventory/library_species.tsv \
        --data-dir <yeast pipeline/data> --work <scratch dir>

Runs the READ-ONLY vendored engine (plugins/yeast/engine/pipeline) with YC_ROOT = <work>/yc_root (links to the data
folder and the engine's pipeline files) and YC_WORK = <work>/yc_work, so nothing is written into the yeast tree or the
engine. Columns: ing_id, kind, gene (genes of a complex, ';'-joined), uniprot, pdb, membrane.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PIPE = os.path.join(HERE, "engine", "pipeline")


def species_table(ing):
    """ingredient table (engine ingredients.build_ingredients) -> species table (ing_id, kind, gene, uniprot, pdb, membrane)"""
    cols = ["ing_id", "kind", "genes" if "genes" in ing.columns else "gene", "uniprot", "pdb", "membrane"]
    out = ing.reindex(columns=[c for c in cols if c in ing.columns]).rename(columns={"genes": "gene"})
    if "gene" not in out.columns:
        out["gene"] = ""
    return out[["ing_id", "kind", "gene", "uniprot", "pdb", "membrane"]].fillna("").astype(str)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--data-dir", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--pdb-dir", default=None)
    a = ap.parse_args()
    root = os.path.join(a.work, "yc_root")
    os.makedirs(root, exist_ok=True)
    links = {"data": a.data_dir}
    for f in os.listdir(PIPE):                     # non-code engine files read via config.ROOT (cell_meshes.npz, ...)
        if not f.endswith(".py") and os.path.isfile(os.path.join(PIPE, f)):
            links[f] = os.path.join(PIPE, f)
    for k, v in links.items():
        p = os.path.join(root, k)
        if not os.path.lexists(p):
            os.symlink(v, p)
    os.environ.update(YC_ROOT=root, YC_WORK=os.path.join(a.work, "yc_work"), YC_LIB=os.path.join(a.work, "yc_lib_unused"),
                      PYTHONHASHSEED=os.environ.get("PYTHONHASHSEED", "0"))
    if a.pdb_dir:
        os.environ["YC_PDB"] = a.pdb_dir
    sys.path.insert(0, PIPE)
    import inventory as inv
    import ingredients as ingm
    df, _geo = inv.build()
    ing, _ = ingm.build_ingredients(df)
    out = species_table(ing)
    pma1 = df.set_index("gene").loc["PMA1", "uniprot"]
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    out.to_csv(a.out, sep="\t", index=False)
    print(f"{len(out)} species ({(out.kind == 'protein').sum()} protein, {(out.kind == 'complex').sum()} complex) -> {a.out}; "
          f"PMA1 uniprot {pma1}", flush=True)


if __name__ == "__main__":
    main()
