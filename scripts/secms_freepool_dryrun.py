"""Dry run of assembly_options.complex_copies = measured_free_pool, or of any other assembly option with --keep-cap
(no repack, no rule sets).

    PYTHONHASHSEED=0 python scripts/secms_freepool_dryrun.py --cell jurkat --reference <build dir> --out <dry-run dir>
        [--set measured_free_pool.exclude_tm=false ...]

Copies the reference build's rna/, v0/ and secms/ outputs (the option changes nothing upstream of the assembly step),
re-runs the steps `assembly` and `ingredients` with the option ON (config override, the cell config is not edited) and
compares with the reference build: complexes whose copies change, protein mass moved between complexes and free pools
(per protein), the knock-on effect on homo-oligomers / AFDB dimers (they consume the free pools), species-level copies
and mass. Writes <out>/dryrun_freepool.json and <out>/dryrun_freepool_complexes.tsv, _proteins.tsv.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--adapter", default=None, help="assembly adapter params key (default: inventory.adapters.assembly)")
    ap.add_argument("--set", nargs="*", default=[], help="extra assembly_options, dotted key=value (YAML values)")
    ap.add_argument("--keep-cap", action="store_true", help="do not switch complex_copies (test another option alone)")
    a = ap.parse_args()
    import yaml
    cfg0 = load(a.cell)
    name = a.adapter or cfg0["inventory"]["adapters"]["assembly"]
    opts = {} if a.keep_cap else {"complex_copies": "measured_free_pool"}
    for kv in a.set:
        k, v = kv.split("=", 1)
        cur = opts
        ks = k.split(".")
        for x in ks[:-1]:
            cur = cur.setdefault(x, {})
        cur[ks[-1]] = yaml.safe_load(v)
    cfg = load(a.cell, overrides={"inventory": {"params": {name: {"assembly_options": opts}}}})
    registry.load_plugins(cfg)
    ref = os.path.join(a.reference, "inventory")
    out = os.path.join(a.out, "inventory")
    if os.path.exists(out):
        sys.exit(f"{out} exists (write-once)")
    for sub in ("rna", "v0", "secms"):
        if os.path.isdir(os.path.join(ref, sub)):
            shutil.copytree(os.path.join(ref, sub), os.path.join(out, sub))
    stage.run(cfg, a.out, step="assembly")
    stage.run(cfg, a.out, step="ingredients")

    rd = lambda base, rel, **k: pd.read_csv(os.path.join(base, rel), sep="\t", **k)
    res = {"cell": a.cell, "reference": a.reference, "assembly_options": opts}
    # complexes
    Cr, Cn = rd(ref, "assembly/complexes_v0.tsv"), rd(out, "assembly/complexes_v0.tsv")
    D = rd(out, "assembly/complex_copies_rule.tsv")
    m = Cr[["cpx", "name", "copies_built", "mass_kda", "mass_built_pg"]].merge(
        Cn[["cpx", "copies_built", "mass_built_pg"]], on="cpx", suffixes=("_cap", "_measured")).merge(
        D[["cpx", "rule", "n_measured", "measured_note"]], on="cpx", how="left")
    m["delta_copies"] = m.copies_built_measured - m.copies_built_cap
    m["delta_mass_pg"] = m.mass_built_pg_measured - m.mass_built_pg_cap
    ch = m[m.delta_copies != 0]
    res["complexes"] = dict(n_parsed=int(len(m)), n_built_cap=int((m.copies_built_cap > 0).sum()),
                            n_built_measured=int((m.copies_built_measured > 0).sum()),
                            n_rule_measured=int((m.rule == "measured_free_pool").sum()), n_changed=int(len(ch)),
                            n_up=int((ch.delta_copies > 0).sum()), n_down=int((ch.delta_copies < 0).sum()),
                            mass_in_complexes_pg_cap=float(m.mass_built_pg_cap.sum()),
                            mass_in_complexes_pg_measured=float(m.mass_built_pg_measured.sum()),
                            mass_up_pg=float(ch.delta_mass_pg.clip(lower=0).sum()),
                            mass_down_pg=float(-ch.delta_mass_pg.clip(upper=0).sum()))
    m.sort_values("delta_mass_pg").to_csv(os.path.join(a.out, "dryrun_freepool_complexes.tsv"), sep="\t", index=False)
    # proteins: copies in Complex Portal complexes
    Ir = rd(ref, "assembly/inventory_v0_complexed.tsv", index_col="gene")
    In = rd(out, "assembly/inventory_v0_complexed.tsv", index_col="gene")
    dc = (In.copies_in_complexes.fillna(0) - Ir.copies_in_complexes.fillna(0))
    mpg = dc * Ir.mass_da.fillna(0) / AV * 1e12
    tot = float(Ir.mass_pg.sum())
    P = pd.DataFrame(dict(uniprot=Ir.uniprot, copies=Ir.copies, frac_in_cpx_cap=Ir.frac_in_complexes,
                          frac_in_cpx_measured=In.frac_in_complexes, delta_copies=dc, delta_mass_pg=mpg))
    P[P.delta_copies != 0].sort_values("delta_mass_pg").to_csv(os.path.join(a.out, "dryrun_freepool_proteins.tsv"), sep="\t")
    res["proteins_complex_portal"] = dict(n_changed=int((dc != 0).sum()), n_more_complexed=int((dc > 0).sum()),
                                          n_less_complexed=int((dc < 0).sum()),
                                          mass_to_complexes_pg=float(mpg.clip(lower=0).sum()),
                                          mass_to_free_pg=float(-mpg.clip(upper=0).sum()), total_protein_pg=tot)
    # knock-on: homo-oligomers / AFDB dimers
    for f, col in (("homo_oligomers_v0.tsv", "mass_pg"), ("afdb_homodimers_v0.tsv", "mass_pg"), ("afdb_heterodimers_v0.tsv", "mass_pg")):
        r, n = rd(ref, f"assembly/{f}"), rd(out, f"assembly/{f}")
        res[f.replace(".tsv", "")] = dict(rows_cap=int(len(r)), rows_measured=int(len(n)), mass_pg_cap=float(r[col].sum()),
                                          mass_pg_measured=float(n[col].sum()))
    cx = lambda base: float(rd(base, "assembly/complexes_v0.tsv").mass_built_pg.sum() +
                            sum(rd(base, f"assembly/{f}").mass_pg.sum() for f in
                                ("homo_oligomers_v0.tsv", "afdb_homodimers_v0.tsv", "afdb_heterodimers_v0.tsv")))
    res["complexed_fraction_of_protein_mass"] = dict(cap=cx(ref) / tot, measured=cx(out) / tot)
    # species level (ingredients step)
    Sr, Sn = rd(ref, "ingredients/ingredients.tsv"), rd(out, "ingredients/ingredients.tsv")
    s = Sr[["ing_id", "copies", "mass_pg"]].merge(Sn[["ing_id", "copies", "mass_pg"]], on="ing_id", how="outer",
                                                   suffixes=("_cap", "_measured")).fillna(0)
    s["d"] = s.mass_pg_measured - s.mass_pg_cap
    res["species"] = dict(n_cap=int(len(Sr)), n_measured=int(len(Sn)), n_new=int((s.copies_cap == 0).sum()),
                          n_gone=int((s.copies_measured == 0).sum()), n_copies_changed=int((s.copies_cap != s.copies_measured).sum()),
                          mass_moved_pg=float(s.d.abs().sum() / 2), total_pg_cap=float(Sr.mass_pg.sum()),
                          total_pg_measured=float(Sn.mass_pg.sum()))
    res["top_complex_changes"] = m.reindex(m.delta_mass_pg.abs().sort_values(ascending=False).index).head(15)[
        ["cpx", "name", "copies_built_cap", "copies_built_measured", "delta_mass_pg", "measured_note"]].to_dict("records")
    json.dump(res, open(os.path.join(a.out, "dryrun_freepool.json"), "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "top_complex_changes"}, indent=1, default=str))


if __name__ == "__main__":
    main()
