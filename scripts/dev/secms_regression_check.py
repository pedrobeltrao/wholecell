"""Criteria S1-S11 of docs/reference/history/secms_regression_criteria.md on a Jurkat inventory build (toolkit preset legacy_v7).
Config: jurkat_v7_port (since Step 3, 2026-11-19, jurkat.yaml runs the final combined SEC-MS with preset default; the
v7 SEC-MS run lives in the v7 port), so the build must be a jurkat_v7_port inventory build.

    PYTHONHASHSEED=0 python scripts/secms_regression_check.py --build <build dir> [--v7repro <port build dir>]

Writes <build>/secms_regression_verdict.json and prints it. Reference tables: {data_root}/data/derived (legacy v7).
S8 uses the orphan sec_assembly_thresholds.tsv: a_k must be equal; its `rule` / `mass_pg` columns are checked against
an 80 pg inventory only when --inv80 <inventory_v0_complexed.tsv of an 80 pg run> is given (the orphan was written
from the stale 80 pg default inventory; 2026-11-16 finding), else reported as differences against this build.
"""
import argparse
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell.config import load  # noqa: E402
from wholecell.inventory import secms  # noqa: E402

MD5 = {"S1": ("assembly/sec_assembly.tsv", "12785dddc48fdb917f19a0eeeb190180"),
       "S2": ("assembly/sec_calibrants_peak.tsv", "8d2bc00fbc9d75db9dcfdeed10f0e628"),
       "S3": ("assembly/homo_oligomers_v0.tsv", "47c8a0a62da9cf22776e23689fe72ac1"),
       "S4": ("assembly/afdb_homodimers_v0.tsv", "205a99f82885b97c0decf9a22929c068"),
       "S5": ("assembly/afdb_heterodimers_v0.tsv", "2beee72eb8410c71990d7d9b04ce0161"),
       "S6": ("assembly/complexes_v0.tsv", "0c95b67246b26e85407f5e6989d2317e"),
       "S9": ("pack_table.tsv", "d1a0099caad792d0ba42cbc52b52f0a5")}


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def make_inv80(build):
    import shutil
    from wholecell import registry
    from wholecell.inventory import stage
    from wholecell.inventory.context import StepContext
    b = build.rstrip("/") + "_s8inv80"
    cfg = load("jurkat_v7_port", overrides={"inventory": {"total_protein_pg": 80}})
    registry.load_plugins(cfg)
    for sub in ("rna", "secms"):
        shutil.copytree(os.path.join(build, "inventory", sub), os.path.join(b, "inventory", sub))
    if "copies" in stage.steps(cfg):                 # copies step (since 2026-11-16) feeds v0
        stage.run(cfg, b, step="copies")
    stage.run(cfg, b, step="v0")
    from wholecell.inventory import assembly as AS
    from wholecell.inventory.assembly import steps as S
    ctx = StepContext(cfg, b, "assembly")
    p = ctx.params("assembly_steps")
    o = AS.options(p)
    st = S.State(ctx, p, o, secms.SecView(ctx, o["secms"]))
    S.curated_machines(st)
    S.complex_portal(st)
    S.write_complex_tables(st)
    return os.path.join(b, "inventory", "assembly", "inventory_v0_complexed.tsv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", required=True)
    ap.add_argument("--v7repro", default=None)
    ap.add_argument("--inv80", default=None)
    ap.add_argument("--make-inv80", action="store_true",
                    help="recompute v0 + complexes at total_protein_pg 80 in <build>_s8inv80/ (for S8 rule/mass_pg)")
    a = ap.parse_args()
    if a.make_inv80:
        a.inv80 = make_inv80(a.build)
    cfg = load("jurkat_v7_port")
    data = cfg["_site"]["data_root"]
    inv = os.path.join(a.build, "inventory")
    v7r = a.v7repro or os.path.join(cfg["_site"]["results_root"], "jurkat", "20261115-ec0b403-s20261115-v7repro2")
    out = {}
    for k, (rel, ref) in MD5.items():
        m = md5(os.path.join(inv, rel))
        out[k] = dict(table=rel, md5=m, reference=ref, pass_=m == ref)
    # S7: content identical to legacy up to order-only compartments cells (= the v7repro2 port table)
    rel = "assembly/inventory_v0_complexed.tsv"
    m, r = md5(os.path.join(inv, rel)), md5(os.path.join(v7r, "inventory", rel))
    cmp = json.load(open(os.path.join(inv, "compare_legacy.json")))["tables"].get(rel, {})
    out["S7"] = dict(table=rel, md5=m, v7repro2_md5=r, pass_=m == r,
                     legacy_content_order_only=cmp.get("content", {}).get("order_only_cells"))
    # S8: orphan thresholds
    P = pd.read_csv(os.path.join(inv, "secms", "secms_protein.tsv"), sep="\t", low_memory=False)
    W = P.pivot(index="uniprot", columns="threshold_k", values="assembled_frac")
    T = pd.read_csv(os.path.join(data, "data", "derived", "sec_assembly_thresholds.tsv"), sep="\t").set_index("uniprot")
    present = bool(T.index.isin(W.index).all())
    dk = {f"a{k}": float((W[k].reindex(T.index) - T[f"a{k}"]).abs().max()) for k in (2, 3, 4, 6)}
    I = pd.read_csv(a.inv80 or os.path.join(inv, "assembly", "inventory_v0_complexed.tsv"), sep="\t")
    I = I[I.copies > 0].drop_duplicates("uniprot").set_index("uniprot")
    rule = I.frac_in_complexes.reindex(T.index).fillna(0)
    nd = int((rule - T.rule).abs().gt(0).sum())
    out["S8"] = dict(all_present=present, max_abs_diff=dk, a_k_pass=present and max(dk.values()) <= 1e-12,
                     rule_rows_differing=nd, rule_max_abs_diff=float((rule - T.rule).abs().max()),
                     rule_inventory=a.inv80 or "this build (100 pg)", pass_=present and max(dk.values()) <= 1e-12 and nd == 0)
    # S10: compare_legacy status = v7repro2
    def status(x):
        return {f"{s}:{k}": (v.get("identical"), json.dumps(v.get("content", {}).get("real_cells"), sort_keys=True))
                for s in ("tables", "ingredients") for k, v in x[s].items()}
    A = status(cmp_all := json.load(open(os.path.join(inv, "compare_legacy.json"))))
    B = status(json.load(open(os.path.join(v7r, "inventory", "compare_legacy.json"))))
    out["S10"] = dict(entries=len(A), differing=[k for k in A if A.get(k) != B.get(k)],
                      summary_equal=cmp_all.get("summary") == json.load(open(os.path.join(v7r, "inventory", "compare_legacy.json"))).get("summary"))
    out["S10"]["pass_"] = not out["S10"]["differing"] and out["S10"]["summary_equal"]
    # S11
    used = json.load(open(os.path.join(inv, "assembly", "assembly_options.json")))
    from wholecell.inventory import assembly as AS
    out["S11"] = dict(config=AS.options(cfg["inventory"]["params"]["assembly_steps"])["complex_portal"]["copies"],
                      build=used["complex_copies"])
    out["S11"]["pass_"] = out["S11"]["config"] == out["S11"]["build"] == "stoichiometry_cap"
    failed = [k for k, v in out.items() if not v["pass_"]]
    out["verdict"] = "REPRODUCED" if not failed else f"NOT REPRODUCED: {failed}"
    out["build"] = a.build
    json.dump(out, open(os.path.join(a.build, "secms_regression_verdict.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
