"""Code check (criteria section 4): outputs of the ported validation run on the archived reference pack vs the stored
legacy validation files (validate.compare_reference.stored_validation). Separate from the verdict.

Expected: loc_holdout / loc_compare identical within 0.001; crowding identical; densities different only through the
documented legacy defect (the v7 run used a fixed mRNA file and put all RNA in the cytosol; the port uses the build's own
RNA sources, incl. nascent RNA in the nucleoplasm); ribosome_nn_v7.json is a byte copy of v6 (not comparable).
Output: <build>/validate/code_check_vs_stored.tsv (file, item, stored, ported, delta, tolerance, status, note)
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from . import common

TOL_LOC = 0.001


def _row(R, f, item, a, b, tol, note="", exact=False, status=None):
    try:
        d = float(b) - float(a)
    except (TypeError, ValueError):
        d = ""
    if status is None:
        if a is None or b is None:
            status = "MISSING"
        elif exact:
            status = "identical" if a == b else "DIFFERS"
        else:
            status = "identical" if a == b else ("within" if isinstance(d, float) and abs(d) <= tol else "DIFFERS")
    R.append(dict(file=f, item=item, stored=a, ported=b, delta=d, tolerance=tol if not exact else "identical", status=status, note=note))


def run(cfg, bdir, **kw):
    cr = (cfg.get("validate", {}) or {}).get("compare_reference", {}) or {}
    st = cr.get("stored_validation", {}) or {}
    v = os.path.join(bdir, "validate")
    R = []
    # loc_holdout
    if st.get("loc_holdout") and os.path.exists(os.path.join(v, "loc_holdout.json")):
        A = json.load(open(st["loc_holdout"])); B = json.load(open(os.path.join(v, "loc_holdout.json")))
        f = os.path.basename(st["loc_holdout"])
        for name, sa in A["sets"].items():
            sb = B["sets"].get(name)
            if sb is None:
                _row(R, f, f"{name}", "present", None, TOL_LOC); continue
            _row(R, f, f"{name}.n_genes", sa["n_genes"], sb["n_genes"], 0, exact=True)
            _row(R, f, f"{name}.mass_pg", sa["mass_pg"], sb["mass_pg"], TOL_LOC)
            for mode in ("intended", "placed"):
                for m, va in (sa.get(mode) or {}).items():
                    vb = (sb.get(mode) or {}).get(m) or {}
                    for agg in ("unweighted", "mass_weighted"):
                        _row(R, f, f"{name}.{mode}.{m}.{agg}", (va or {}).get(agg), vb.get(agg), TOL_LOC)
    # loc_compare.tsv
    if st.get("loc_compare") and os.path.exists(os.path.join(v, "loc_compare.tsv")):
        A = pd.read_csv(st["loc_compare"], sep="\t", index_col=0, low_memory=False)
        B = pd.read_csv(os.path.join(v, "loc_compare.tsv"), sep="\t", index_col=0, low_memory=False)
        f = os.path.basename(st["loc_compare"])
        _row(R, f, "genes", len(A), len(B), 0, note=f"only stored {len(set(A.index) - set(B.index))}, only ported "
             f"{len(set(B.index) - set(A.index))}", exact=True)
        g = A.index.intersection(B.index)
        miss_cols = [c for c in A.columns if c not in B.columns]
        if miss_cols:
            _row(R, f, "columns missing in port", "", ";".join(miss_cols), 0, status="DIFFERS")
        num = [c for c in A.columns if c in B.columns and pd.api.types.is_numeric_dtype(A[c]) and pd.api.types.is_numeric_dtype(B[c])]
        txt = [c for c in A.columns if c in B.columns and c not in num]
        nbad_genes = np.zeros(len(g), bool)
        for c in num:
            a = A.loc[g, c].astype(float).values; b = B.loc[g, c].astype(float).values
            both_nan = np.isnan(a) & np.isnan(b)
            d = np.where(both_nan, 0.0, np.abs(a - b)); d = np.where(np.isnan(d), np.inf, d)
            nbad_genes |= d > TOL_LOC
            md = float(d.max()) if len(d) else 0.0
            if md > 0:
                _row(R, f, f"col {c}: max |delta| (genes > {TOL_LOC})", 0.0, md, TOL_LOC, note=f"{int((d > TOL_LOC).sum())} genes")
        _row(R, f, f"numeric columns ({len(num)}): genes with any |delta| > {TOL_LOC}", 0, int(nbad_genes.sum()), 0, exact=True)
        nt = 0
        for c in txt:
            a = A.loc[g, c].astype(object).where(A.loc[g, c].notna(), "").astype(str).values
            b = B.loc[g, c].astype(object).where(B.loc[g, c].notna(), "").astype(str).values
            k = int((a != b).sum()); nt += k
            if k:
                _row(R, f, f"col {c}: genes differing (text)", 0, k, 0, exact=True)
        _row(R, f, f"text columns ({len(txt)}): differing cells", 0, nt, 0, exact=True)
    # crowding
    if st.get("crowding") and os.path.exists(os.path.join(v, "crowding.json")):
        A = json.load(open(st["crowding"])); B = json.load(open(os.path.join(v, "crowding.json")))
        f = os.path.basename(st["crowding"])
        _row(R, f, "r_bead", A["r_bead"], B["r_bead"], 0, exact=True)
        sa = {r["region"]: r for r in A["summary"]}; sb = {r["region"]: r for r in B["summary"]}
        n_diff = 0
        ta, tb = st["crowding"].replace(".json", ".tsv"), os.path.join(v, "crowding.tsv")
        expl = ""
        if os.path.exists(ta) and os.path.exists(tb):
            la = pd.read_csv(ta, sep="\t")[["region", "tile"]].values.tolist(); lb = pd.read_csv(tb, sep="\t")[["region", "tile"]].values.tolist()
            same_boxes = la == lb
            _row(R, f, "sampled boxes (region, tile)", len(la), sum(x == y for x, y in zip(la, lb)), 0, exact=True,
                 status="identical" if same_boxes else "DIFFERS (explained)",
                 note="" if same_boxes else "legacy fallback box choice used np.argsort quicksort: tied tiles ordered differently "
                 "on another CPU; port uses a stable sort -> other boxes, hence other values")
            if not same_boxes:
                expl = "different sampled boxes (see row above)"
        for reg, ra in sa.items():
            for k, x in ra.items():
                if k == "region": continue
                y = sb.get(reg, {}).get(k)
                if x != y:
                    n_diff += 1; _row(R, f, f"{reg}.{k}", x, y, 0, exact=True, note=expl, status="DIFFERS (explained)" if expl else None)
        _row(R, f, "summary values differing", 0, n_diff, 0, exact=True, note=expl,
             status="DIFFERS (explained)" if (n_diff and expl) else None)
    # densities
    if st.get("densities") and os.path.exists(os.path.join(v, "densities.json")):
        A = json.load(open(st["densities"])); B = json.load(open(os.path.join(v, "densities.json")))
        f = os.path.basename(st["densities"])
        ta = {r["region"]: r for r in A["rows"]}; tb = {r["region"]: r for r in B["rows"]}
        note_rna = "legacy defect: fixed mRNA file, all RNA in cytosol; port: build's own RNA sources (mRNA -> cytosol, nascent -> nucleoplasm)"
        for reg, ra in ta.items():
            for k, x in ra.items():
                if k == "region": continue
                y = tb.get(reg, {}).get(k)
                expl = reg in (B.get("mrna_pg_by_region") or {}) and k in ("prepass_pg", "mg_per_ml")
                _row(R, f, f"{reg}.{k}", x, y, 0, note=note_rna if (x != y and expl) else "", exact=True,
                     status=None if x == y or not expl else "DIFFERS (explained)")
        for k in ("nuc_cyto_ratio", "mrna_pg"):
            x, y = A.get(k), B.get(k)
            _row(R, f, k, x, y, 0, exact=True, note=note_rna if x != y else "", status=None if x == y else "DIFFERS (explained)")
        for grp in ("mem_pg", "shell_pg"):
            nd = sum(A[grp][r] != B.get(grp, {}).get(r) for r in A[grp])
            _row(R, f, f"{grp}: regions differing", 0, nd, 0, exact=True)
    # ribosome_nn
    if st.get("ribosome_nn") and os.path.exists(os.path.join(v, "ribosome_nn.json")):
        A = json.load(open(st["ribosome_nn"])); B = json.load(open(os.path.join(v, "ribosome_nn.json")))
        f = os.path.basename(st["ribosome_nn"])
        ra = {r["set"]: r for r in A["rows"]}; rb = {r["set"]: r for r in B["rows"]}
        for s in ra:
            if s in rb:
                _row(R, f, f"{s}: n / nn_median", f"{ra[s]['n']} / {ra[s]['nn_median']}", f"{rb[s]['n']} / {rb[s]['nn_median']}", "",
                     status="not comparable", note="stored file is a byte copy of ribosome_nn_v6.json (legacy defect: VAL_POLY not set)")
    T = pd.DataFrame(R, columns=["file", "item", "stored", "ported", "delta", "tolerance", "status", "note"])
    p = os.path.join(v, "code_check_vs_stored.tsv")
    common.write_df(p, T, overwrite=True, index=False)
    print(T.status.value_counts().to_string())
    return T
