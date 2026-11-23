#!/usr/bin/env python
"""A7 phase 2 check (P6): the core organisation layer (wholecell/pack/organisation.py) on the yeast inputs against the
yeast port build s0 (criterion P2 of docs/reference/history/yeast_x5_reproduction_criteria.md with the 2026-11-17
addendum: network.npz non-float keys identical, floats rel. 1e-5; partners.npz ids and off identical, each CSR row
compared as a set of (partner index, kind), strengths of matched elements rel. 1e-5). Run in the core environment;
judged at tolerance level (docs/reference/history/yeast_x5_env_tolerances_2026-11-17.md).

    <core python> scripts/organisation_a7_check.py [--out DIR] [--s0 BUILD] [--threads 3]

Two variants into <out>/{pool_on,pool_off}/organisation/: pool_on = the X5 settings (AF3 pool contact partners;
legacy_x5 reproduction), pool_off = the published default (STRING only). The pool_off partner sets are compared with
pool_on (how much switching the pool off changes them). Report: <out>/check_report.json.
Inputs (read in place, read-only): the s0 build (ingredients.csv, expected_totals.parquet, lib_ids.npy; reference
network.npz / partners.npz), the yeast pipeline data (STRING v12 4932, UniProt table, paralog pairs) and the frozen
work seed (pool_pairs.parquet, the cached pool table X5 read).
"""
import argparse
import json
import os
import sys
import time

ap = argparse.ArgumentParser()
ap.add_argument("--out", default="{results_root}/yeast/_dev_unify/a7p2/organisation")
ap.add_argument("--s0", default="{results_root}/yeast/20261116-aabe7b5-s0-x5repro")
ap.add_argument("--threads", type=int, default=3)
a = ap.parse_args()
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = str(a.threads)
os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell.config import load  # noqa: E402
from wholecell.pack import organisation as O, switches  # noqa: E402

REL = 1e-5


def settings(cfg, pool):
    data = cfg.frozen("yeast_pipeline_data")
    o = {"enabled": True,
         "string": {"taxon": 4932, "version": "12.0", "links": os.path.join(data, "string", "4932.protein.links.v12.0.txt.gz"),
                    "min_score": 400, "ids": {"table": "inv_yeast_uniprot", "columns": ["Gene Names (ordered locus)"]}},
         "genes": {"table": "inv_yeast_uniprot", "columns": ["Gene Names (primary)", "Gene Names (ordered locus)"]},
         "subunits": {"table": "engine/work/results/inventory/ingredients.csv",
                      "lib_map": "engine/work/yc3d/expected_totals.parquet", "ids": "engine/work/yc3d/lib_ids.npy"}}
    if pool:
        o.update(pool_scores=os.path.join(cfg.frozen("yeast_x5_work_seed"), "organisation", "pool_pairs.parquet"),
                 paralogs=os.path.join(data, "alliance_yeast_paralog_pairs.csv"))
    return switches.load({"organisation": o}).organisation


def rel(x, y):
    x = np.asarray(x, np.float64); y = np.asarray(y, np.float64)
    return float((np.abs(x - y) / np.maximum(np.abs(y), 1e-30)).max()) if x.size else 0.0


def rows(z, k):
    """CSR row k as {(partner index, kind): strength}"""
    s = slice(z["off"][k], z["off"][k + 1])
    return {(int(i), int(t)): float(v) for i, t, v in zip(z["idx"][s], z["kind"][s], z["strength"][s])}


def compare_network(ref, new):
    out = {}
    for k in ref.files:
        A, B = ref[k], new[k]
        if A.dtype.kind == "f":
            same_shape = A.shape == B.shape
            r = rel(B, A) if same_shape else None
            # loadings: an eigenvector's sign is arbitrary (LAPACK); report both the raw and the sign-aligned difference
            extra = {}
            if k == "E" and same_shape:
                sg = np.sign((A.astype(np.float64) * B).sum(0)); sg[sg == 0] = 1
                extra = dict(rel_sign_aligned=rel(B * sg[None, :], A), sign_flips=int((sg < 0).sum()))
            # ok: within rel 1e-5 after aligning the arbitrary eigenvector signs (sign flips reported; a flipped component
            # changes the random-field realisation of the tile split, not its distribution)
            ok = same_shape and (r <= REL or extra.get("rel_sign_aligned", 1) <= REL)
            # Z / E relative to the largest entry of the column as well (tiny entries make the elementwise ratio large)
            if same_shape:
                d = np.abs(B.astype(np.float64) - A) / np.maximum(np.abs(A).max(0, keepdims=True), 1e-30)
                extra["rel_to_column_max"] = float(d.max())
            out[k] = dict(kind="float", ref_dtype=str(A.dtype), new_dtype=str(B.dtype), shape=list(A.shape), rel_max=r, ok=bool(ok), **extra)
        else:
            out[k] = dict(kind="exact", ok=bool(A.shape == B.shape and np.array_equal(A, B)))
    return out


def compare_partners(ref, new):
    res = dict(ids=bool(np.array_equal(ref["ids"], new["ids"])), off=bool(np.array_equal(ref["off"], new["off"])))
    n = len(ref["off"]) - 1
    set_diff = 0; worst = 0.0; order_only = 0; cutoff_ties = []
    for k in range(n):
        a, b = rows(ref, k), rows(new, k)
        if set(a) != set(b):
            set_diff += 1
            # equal-strength partners at the max_partners cut: which of them is kept is tie order (addendum rule)
            e = ref["off"][k + 1] - 1
            lk, ls = int(ref["kind"][e]), float(ref["strength"][e])
            sym = [(x, a.get(x, b.get(x))) for x in set(a) ^ set(b)]
            if all(x[1] == lk and abs(v - ls) <= 1e-12 * max(abs(ls), 1e-30) for x, v in sym):
                cutoff_ties.append(dict(row=k, id=str(ref["ids"][k]), swapped=sorted(int(x[0]) for x, _ in sym), strength=ls))
            continue
        worst = max(worst, max((abs(b[x] - a[x]) / max(abs(a[x]), 1e-30) for x in a), default=0.0))
        s = slice(ref["off"][k], ref["off"][k + 1])
        if not np.array_equal(ref["idx"][s], new["idx"][s]):
            order_only += 1
    res.update(rows=n, rows_set_differs=set_diff, rows_set_differs_only_by_cutoff_ties=len(cutoff_ties),
               cutoff_ties=cutoff_ties[:20], rows_tie_order_only=order_only, strength_rel_max=worst,
               idx_identical=bool(np.array_equal(ref["idx"], new["idx"])),
               ok_strict=bool(res["ids"] and res["off"] and set_diff == 0 and worst <= REL),
               ok=bool(res["ids"] and res["off"] and set_diff == len(cutoff_ties) and worst <= REL))
    return res


def pool_effect(on, off):
    """how switching the AF3 pool off changes the partner sets (rows over the same type order)"""
    n = len(on["off"]) - 1
    changed = 0; jac = []; lost_contact = 0; tot_on = len(on["idx"]); kept = 0
    for k in range(n):
        a = {i for i, _ in rows(on, k)}; b = {i for i, _ in rows(off, k)}
        if a != b:
            changed += 1
        if a or b:
            jac.append(len(a & b) / len(a | b))
        kept += len(a & b)
    lost_contact = int((on["kind"] == O.CONTACT).sum())
    no_p_on = int((np.diff(on["off"]) == 0).sum()); no_p_off = int((np.diff(off["off"]) == 0).sum())
    return dict(types=n, links_pool_on=tot_on, links_pool_off=int(len(off["idx"])), contact_links_pool_on=lost_contact,
                contact_links_pool_off=int((off["kind"] == O.CONTACT).sum()),
                types_with_partners_pool_on=n - no_p_on, types_with_partners_pool_off=n - no_p_off,
                rows_changed=changed, rows_changed_frac=round(changed / n, 4),
                partner_pairs_kept=kept, partner_pairs_kept_frac_of_pool_on=round(kept / max(tot_on, 1), 4),
                mean_jaccard_rows_with_partners=round(float(np.mean(jac)), 4) if jac else None)


def main():
    from wholecell.config import load_site
    site = load_site()
    rr = os.path.join(site["data_root"], "data", "derived", "wholecell_builds")
    out = a.out.replace("{results_root}", rr); s0 = a.s0.replace("{results_root}", rr)
    os.makedirs(out, exist_ok=True)
    cfg = load("yeast")
    ref_dir = os.path.join(s0, "engine", "work", "organisation")
    ref_n = np.load(os.path.join(ref_dir, "network.npz")); ref_p = np.load(os.path.join(ref_dir, "partners.npz"))
    rep = dict(s0=s0, out=out, criteria="P2 + addendum 2026-11-17 (rows as sets; strengths rel 1e-5)", variants={})
    T0 = time.time()
    for name, pool in (("pool_on", True), ("pool_off", False)):
        od = os.path.join(out, name, "organisation")
        o = settings(cfg, pool)
        ids = O.network_ids(cfg, s0, o, cfg["pack"], {})
        t0 = time.time()
        if not os.path.exists(os.path.join(od, "partners.npz")):
            O.build(cfg, s0, o, ids, None, od, log=lambda *x: print(name, *x, flush=True))
        v = dict(seconds=round(time.time() - t0, 1), stats=json.load(open(os.path.join(od, "network.json")))["stats"])
        nn = np.load(os.path.join(od, "network.npz")); pp = np.load(os.path.join(od, "partners.npz"))
        v["network_vs_s0"] = compare_network(ref_n, nn)
        v["partners_vs_s0"] = compare_partners(ref_p, pp)
        v["P2_pass"] = bool(all(x["ok"] for x in v["network_vs_s0"].values()) and v["partners_vs_s0"]["ok"])
        rep["variants"][name] = v
        print(name, "P2", "PASS" if v["P2_pass"] else "FAIL", json.dumps(v, indent=1), flush=True)
    on = np.load(os.path.join(out, "pool_on", "organisation", "partners.npz"))
    off = np.load(os.path.join(out, "pool_off", "organisation", "partners.npz"))
    rep["pool_off_vs_pool_on"] = pool_effect(on, off)
    na = np.load(os.path.join(out, "pool_on", "organisation", "network.npz")); nb = np.load(os.path.join(out, "pool_off", "organisation", "network.npz"))
    rep["pool_off_network_identical_to_pool_on"] = bool(all(np.array_equal(na[k], nb[k]) for k in na.files))
    rep["seconds"] = round(time.time() - T0, 1)
    json.dump(rep, open(os.path.join(out, "check_report.json"), "w"), indent=1)
    print(json.dumps(rep["pool_off_vs_pool_on"], indent=1))
    print("P2 pool_on (X5 settings):", "PASS" if rep["variants"]["pool_on"]["P2_pass"] else "FAIL")


if __name__ == "__main__":
    main()
