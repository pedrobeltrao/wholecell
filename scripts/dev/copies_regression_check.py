"""Criteria C1-C5, C7 of docs/reference/history/copies_regression_criteria.md on a Jurkat inventory build made with the copies step
(C6 is tests/test_copies.py).

    PYTHONHASHSEED=0 python scripts/copies_regression_check.py --build <build dir> [--v7repro <port build dir>]

Writes <build>/copies_regression_verdict.json and prints it.
"""
import argparse
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pandas as pd  # noqa: E402

from wholecell.config import load  # noqa: E402
from wholecell.inventory import compare  # noqa: E402

V0_PORT_MD5 = "b40f1ba260f8a2460d62160bbbe68342"
PACK_MD5 = "d1a0099caad792d0ba42cbc52b52f0a5"


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def status(cj):
    """per entry: (identical, content identical, any-row-order identical, B3 pass)"""
    out = {}
    for sec in ("tables", "ingredients"):
        for k, v in cj.get(sec, {}).items():
            c = v.get("content") or {}
            out[k] = (v.get("identical"), c.get("identical_content"), c.get("identical_content_any_row_order"),
                      (v.get("B3") or {}).get("pass") if isinstance(v.get("B3"), dict) else v.get("B3"))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", required=True)
    ap.add_argument("--v7repro", default=None)
    a = ap.parse_args()
    cfg = load("jurkat")
    dr = cfg["_site"]["data_root"]
    rs = cfg["_site"]["results_root"]
    ref = a.v7repro or os.path.join(rs, "jurkat", "20261115-ec0b403-s20261115-v7repro2")
    inv = os.path.join(a.build, "inventory")
    res = {}
    v0 = os.path.join(inv, "v0", "inventory_v0.tsv")
    res["C1"] = dict(md5=md5(v0), expected=V0_PORT_MD5, pass_=md5(v0) == V0_PORT_MD5)
    leg = os.path.join(dr, "data", "derived", "inventory_v0.tsv")
    A = pd.read_csv(v0, sep="\t", keep_default_na=False, dtype=str)
    B = pd.read_csv(leg, sep="\t", keep_default_na=False, dtype=str)
    cd = compare._content_diff(v0, leg, rtol=1e-14)
    res["C2"] = dict(same_genes_order=list(A.gene) == list(B.gene), copies_identical=bool((A.copies == B.copies).all()),
                     content_identical_rtol_1e14=cd.get("identical_content"), order_only_cells=cd.get("order_only_cells"),
                     float_cells=cd.get("float_cells"), max_rel_diff=cd.get("max_rel_diff"), real_cells=cd.get("real_cells"))
    res["C2"]["pass_"] = bool(res["C2"]["same_genes_order"] and res["C2"]["copies_identical"] and cd.get("identical_content"))
    C = pd.read_csv(os.path.join(inv, "copies", "copies.tsv"), sep="\t", keep_default_na=False, dtype=str)
    res["C3"] = dict(copies_equal=bool((C.copies.values == A.copies.values).all()), source_equal=bool((C.source.values == A.source.values).all()))
    res["C3"]["pass_"] = res["C3"]["copies_equal"] and res["C3"]["source_equal"]
    pt = md5(os.path.join(inv, "pack_table.tsv"))
    res["C4"] = dict(md5=pt, expected=PACK_MD5, pass_=pt == PACK_MD5)
    s_new = status(json.load(open(os.path.join(inv, "compare_legacy.json"))))
    s_ref = status(json.load(open(os.path.join(ref, "inventory", "compare_legacy.json"))))
    diff = {k: dict(new=s_new.get(k), v7repro2=s_ref.get(k)) for k in set(s_new) | set(s_ref) if s_new.get(k) != s_ref.get(k)}
    cj = json.load(open(os.path.join(inv, "compare_legacy.json")))
    res["C5"] = dict(entries=len(s_new), status_differences=diff, B3_all_pass=cj["summary"].get("B3_all_pass"),
                     pass_=not diff and bool(cj["summary"].get("B3_all_pass")))
    L = pd.read_csv(os.path.join(inv, "copies", "copy_ledger.tsv"), sep="\t", keep_default_na=False, dtype=str)
    rules_by_key = L[L.key != "*"].groupby("key").rule_id.agg(lambda x: set(x))
    ok_led = all(r.set_by in rules_by_key.get(r.gene, set()) or f"source.{r.set_by}" in rules_by_key.get(r.gene, set())
                 for r in C.itertuples())
    res["C7"] = dict(all_set_by=bool((C.set_by != "").all()), ledger_rows=len(L), rule_ids=sorted(set(L.rule_id)),
                     every_rule_change_logged=ok_led, set_by_counts=C.set_by.value_counts().to_dict())
    res["C7"]["pass_"] = bool(res["C7"]["all_set_by"] and ok_led)
    res["verdict"] = "REPRODUCED" if all(res[k]["pass_"] for k in ("C1", "C2", "C3", "C4", "C5", "C7")) else (
        "REPRODUCED (copies)" if all(res[k]["pass_"] for k in ("C2", "C3", "C4", "C5")) else "FAILED")
    json.dump(res, open(os.path.join(a.build, "copies_regression_verdict.json"), "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))


if __name__ == "__main__":
    main()
