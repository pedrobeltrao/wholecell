"""Yeast copy numbers through the shared copy-number step vs yeast build W (docs/reference/history/copies_regression_criteria.md Y1-Y4).

    PYTHONHASHSEED=0 python scripts/copies_yeast_check.py --out <new dir> [--set k=v ...]

Runs only the `copies` step with config/cells/yeast.yaml (inputs read in place from {yeast_root}, read-only; the yeast
code is NOT imported or run) into <out>/inventory/copies/, then compares with build W's per-protein copies
(protein_features.parquet), its s1 log numbers and anchors.json. Writes <out>/copies_yeast_verdict.json.
--set copies.total_protein.mode=anchors etc. runs a variant (diagnostic; the verdict then only reports differences).
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.config import load  # noqa: E402
from wholecell.inventory import copies  # noqa: E402
from wholecell.inventory.context import StepContext  # noqa: E402

# build W s1 log (runs/P2_seed0/build.log, 15:14:35): NOT_IN_THIS_CELL and the ribosome scaling
W_ABSENT = {"absent in G1": (8, 4.509649784355788e-05), "alpha-specific": (8, 7.153518537785694e-06),
            "genotype": (4, 0.0008410272791149524), "glucose-repressed": (11, 0.0017065072888808497),
            "meiosis and sporulation": (42, 7.910830095971548e-05)}
W_RIBO = dict(target=262722, factor=1.185, protein_pg=3.61)


def _set(d, dotted, v):
    ks = dotted.split(".")
    for k in ks[:-1]:
        d = d.setdefault(k, {})
    try:
        v = json.loads(v)
    except ValueError:
        pass
    d[ks[-1]] = v


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", action="append", default=[])
    a = ap.parse_args()
    ov = {"inventory": {}}
    for s in a.set:
        k, v = s.split("=", 1)
        _set(ov["inventory"], k, v)
    cfg = load("yeast", overrides=ov if a.set else None)
    registry.load_plugins(cfg)
    os.makedirs(os.path.join(a.out, "inventory"), exist_ok=True)
    ctx = StepContext(cfg, a.out, "copies")
    T, S = copies.step(ctx)
    W = pd.read_parquet(cfg.frozen("val_yeast_buildW_protein_features"))
    W = W.set_index("orf")
    M = T.reset_index().set_index("orf")
    res = dict(variant=a.set)
    res["Y1_same_proteins"] = set(W.index) == set(M.index) and len(W) == len(M)
    m = M.copies.reindex(W.index).astype(float)
    w = W.copies_per_cell.astype(float)
    rel = (m - w).abs() / w.where(w > 0)
    res["Y2_max_rel_diff"] = float(rel.max())
    res["Y2_zero_match"] = bool(((w == 0) == (m == 0)).all())
    res["Y2_pass"] = bool(res["Y2_max_rel_diff"] <= 1e-9 and res["Y2_zero_match"])
    worst = rel.sort_values(ascending=False).head(5)
    res["Y2_worst"] = {k: [float(m[k]), float(w[k])] for k in worst.index}
    ab = {k.split(".", 1)[1]: (v["proteins"], v["share_of_value"]) for k, v in S["absent"].items()}
    res["Y3_absent"] = {k: dict(ours=ab.get(k), W=list(v)) for k, v in W_ABSENT.items()}
    ab_ok = all(ab.get(k) and ab[k][0] == v[0] and abs(ab[k][1] - v[1]) <= 1e-9 * max(v[1], 1e-12) + 1e-15
                for k, v in W_ABSENT.items())
    adj = S["adjustments"].get("ribosomes_cryoet_density", {})
    res["Y3_ribosomes"] = dict(target=adj.get("target"), factor=adj.get("factor"), protein_pg=adj.get("protein_pg_after"), W=W_RIBO)
    rib_ok = bool(adj) and round(adj["target"]) == W_RIBO["target"] and round(adj["factor"], 3) == W_RIBO["factor"] \
        and round(adj["protein_pg_after"], 2) == W_RIBO["protein_pg"]
    res["Y3_pass"] = bool(ab_ok and rib_ok)
    WA = json.load(open(cfg.frozen("val_yeast_buildW_anchors")))
    an = S["anchors"]
    y4 = {}
    if "nucleoporins" in an:
        ours = an["nucleoporins"]["model_pores_W"]["ratio"]
        y4["nucleoporins_median"] = [ours, WA["nucleoporins_per_pore_need"]["median"]]
    if "core_histones" in an:
        A = pd.read_csv(os.path.join(a.out, "inventory", "copies", "anchors.tsv"), sep="\t")
        h = A[(A.anchor_id == "core_histones") & (A.count_label == "build_W_nucleosomes")].set_index("family").ratio
        y4.update({f"histone_{k}": [float(h[k]), v] for k, v in WA["core_histones_per_nucleosome_need"].items()})
    res["Y4"] = y4
    res["Y4_pass"] = bool(y4) and all(abs(x - y) <= 1e-6 * abs(y) for x, y in y4.values())
    res["anchor_summary"] = S.get("anchor_summary")
    res["total_protein"] = S["total_protein"]
    res["verdict"] = "REPRODUCED" if (res["Y1_same_proteins"] and res["Y2_pass"] and res["Y3_pass"] and res["Y4_pass"]) else "DIFFERS"
    json.dump(res, open(os.path.join(a.out, "copies_yeast_verdict.json"), "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
