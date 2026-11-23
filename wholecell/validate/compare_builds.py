"""Build-vs-reference comparison implementing docs/reference/history/v7_reproduction_criteria.md (rows B1-B12, V1-V8).

    python -m wholecell validate --cell jurkat --build <port build> --step compare_v7 [port_pack_dir=...] [ref_validate_dir=...]

Reference = config validate.compare_reference (legacy file map `files`, the ported validation run on the reference pack
`validate_dir`, expected checksums/totals); port = the build directory (pack/, inventory/, preplaced/, geometry/, validate/; a build made before
the A0c stage rename is read in its old names through wholecell.build.ref_path / common.Inputs).
Outputs: <port build>/validate/compare_v7/criteria_results.tsv (id, metric, v7, port, delta, tolerance, level, class, note)
and verdict.json. Levels: E (identical), T (within tolerance), FAIL, NA (not evaluated: input missing or not ported).
Verdict (criteria section 1): REPRODUCED if every must row is E or T; REPRODUCED WITH DEVIATIONS if every must row passes
but a report row fails; NOT REPRODUCED if a must row fails. A must row that is NA is listed in verdict.json
('must_not_evaluated') and prevents REPRODUCED (counted as not passed).
Rows with sub-items (per region, per set) carry ids like 'V1:cytosol'; every sub-row counts.

The thresholds below are the constants of the criteria document (reviewed at commit). Do not change them without a
dated addendum in the document.
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd

from . import common
from ..build import ref_path

# ------------------------------------------------------------------ thresholds (docs/reference/history/v7_reproduction_criteria.md)
B1_LABEL_IDENTICAL_FRAC = 0.9999     # labels >= 99.99 % identical
B1_SDF_ABS_NM = 0.01                 # membrane SDF |d| < 0.01 nm (99.9th pct)
B1_N_UNIFORM, B1_N_NEAR = 2_000_000, 500_000
B2_REL = 0.001                       # region volume / membrane area within 0.1 %
B3_LOG2 = 0.01                       # per-species copies |log2 ratio| <= 0.01
B3_CW = 1e-6                         # comp_weights max |delta| <= 1e-6
B5_REL = 0.001                       # lamina count +- 0.1 %
B6_TOTAL_REL = 0.001                 # total requested instances +- 0.1 %
B6_PAIR_SIGMA, B6_PAIR_REL = 3.0, 0.005    # per (species, region) totals within max(3 sqrt n, 0.5 % n)
B9_SIGMA, B9_REL = 3.0, 0.01         # within max(3 sqrt n_v7, 1 % n_v7) for n_v7 >= 10
B9_NMIN = 10
B9_SMALL_FRAC, B9_SMALL_ABS = 0.99, 3   # >= 99 % of pairs with n_v7 < 10 within +- 3
B10_FRAC = 0.005                     # frac within +- 0.005 per region
B11_INST_REL = 0.005                 # instances +- 0.5 %
B11_MASS_FRAC = 0.003                # mass_frac +- 0.003
B12_SAMPLE = 200                     # tiles compared array by array (plus full count equality, B8)
V1_MGML, V1_RATIO = 1.0, 0.01        # +- 1 mg/ml per region; nuclear/cytoplasmic ratio +- 0.01
V2_PHI = 0.005                       # excluded-volume fraction +- 0.005
V3_TOL = 0.005                       # agree / agree_tvd / tophit, unweighted and mass-weighted, +- 0.005
V3_GENE_TVD, V3_GENE_FRAC = 0.02, 0.99    # >= 99 % of genes within 0.02 TVD (placed compartment shares)
V3_SETS = ("opencell", "hpa_consensus_all", "hpa_jurkat")
V4_TOL = 0.01                        # hold-out half +- 0.01
V4_SET = "hpa_jurkat_holdout_half"
V6_REL, V6_ABS = 0.10, 1.0           # max(10 % relative, 1 per 1e6 instances absolute)
V6_MIN_TILES = 300
V7_NN = 0.3                          # median +- 0.3 nm
SCORE_METRICS = ("agree", "agree_tvd", "tophit")
CLASS = dict(B1="must", B2="must", B3="must", B4="must", B5="report", B6="must", B7="must", B8="report", B9="must",
             B10="must", B11="must", B12="report", V1="must", V2="report", V3="must", V4="must", V5="must", V6="must",
             V7="report", V8="report")


class Rows(list):
    def add(self, id_, metric, v7, port, tol, level, note="", delta=None):
        if delta is None:
            try:
                delta = float(port) - float(v7)
            except (TypeError, ValueError):
                delta = ""
        self.append(dict(id=id_, metric=metric, v7=v7, port=port, delta=delta, tolerance=tol, level=level,
                         **{"class": CLASS[id_.split(":")[0]]}, note=note))


def _lvl(identical, within):
    return "E" if identical else ("T" if within else "FAIL")


def _load_json(p):
    return json.load(open(p)) if p and os.path.exists(p) else None


def _md5_or_none(p):
    return common.md5(p) if p and os.path.exists(p) else None


# ------------------------------------------------------------------ B rows
def b1(rows, port_geom_val):
    V = _load_json(port_geom_val)
    if V is None:
        return rows.add("B1", "geometry engine vs legacy cell3d (labels, membrane SDFs)", "", "", "", "NA",
                        f"geometry validation output missing ({port_geom_val})")
    b = (V.get("checks", {}) or {}).get("b1_legacy_equality") or V.get("b1_legacy_equality")
    if not b:
        return rows.add("B1", "geometry engine vs legacy cell3d", "", "", "", "NA", "no b1_legacy_equality in geometry validation")
    sets = b.get("sets", {})
    n_u = sets.get("uniform", {}).get("n", 0); n_n = sets.get("near_membrane_50nm", {}).get("n", 0)
    fr = min(s.get("label_identical_frac", 0.0) for s in sets.values())
    md = max(s.get("mem_sdf_max_absdiff_nm", np.inf) for s in sets.values())
    exact = all(s.get("label_mismatch", 1) == 0 and s.get("mem_sdf_exact_equal", False) for s in sets.values())
    enough = n_u >= B1_N_UNIFORM and n_n >= B1_N_NEAR
    lvl = _lvl(exact and enough, enough and fr >= B1_LABEL_IDENTICAL_FRAC and md < B1_SDF_ABS_NM)
    rows.add("B1", "geometry: min label-identical fraction / max |membrane SDF diff| (nm)", "1.0 / 0.0", f"{fr} / {md}",
             f"labels >= {B1_LABEL_IDENTICAL_FRAC}; |dSDF| < {B1_SDF_ABS_NM} nm", lvl,
             f"from geometry validation ({os.path.basename(port_geom_val)}); sets {', '.join(f'{k} n={v.get('n')}' for k, v in sets.items())}",
             delta=md)


def b2(rows, ref, port):
    pr, pp = ref["tiles_parquet"], port["tiles_parquet"]
    if not os.path.exists(pp):
        return rows.add("B2", "tiles.parquet (region volumes, membrane areas)", "", "", "", "NA", f"missing {pp}")
    if common.md5(pr) == common.md5(pp):
        return rows.add("B2", "tiles.parquet md5", common.md5(pr), common.md5(pp), "identical", "E", "", delta=0)
    A = pd.read_parquet(pr); Bp = pd.read_parquet(pp)
    cols = [c for c in A.columns if c not in ("tile", "ti", "tj", "tk")]
    worst = 0.0; bad = []
    for c in cols:
        a = float(A[c].sum()); b = float(Bp[c].sum()) if c in Bp.columns else float("nan")
        rel = abs(b - a) / a if a > 0 else (0.0 if b == 0 else float("inf"))
        if not rel <= B2_REL: bad.append(c)
        worst = max(worst, rel) if rel == rel else float("inf")
    rows.add("B2", "region volume / membrane area, max relative difference", 0.0, worst, f"<= {B2_REL}",
             "T" if not bad else "FAIL", f"md5 differs; out of tolerance: {bad}" if bad else "md5 differs; all measures within tolerance",
             delta=worst)


def _species_table_check(a_path, b_path):
    A = pd.read_csv(a_path, sep="\t"); B = pd.read_csv(b_path, sep="\t")
    sa, sb = set(A.ing_id), set(B.ing_id)
    A = A.set_index("ing_id"); B = B.set_index("ing_id")
    common_ = sorted(sa & sb)
    ca = A.loc[common_, "copies"].astype(float); cb = B.loc[common_, "copies"].astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        lr = np.abs(np.log2(np.where((ca > 0) & (cb > 0), cb / ca, np.where(ca == cb, 1.0, np.inf))))
    max_lr = float(np.nanmax(lr)) if len(lr) else 0.0
    max_cw = 0.0
    for s in common_:
        x = json.loads(A.at[s, "comp_weights"]) if isinstance(A.at[s, "comp_weights"], str) else {}
        y = json.loads(B.at[s, "comp_weights"]) if isinstance(B.at[s, "comp_weights"], str) else {}
        for k in set(x) | set(y):
            max_cw = max(max_cw, abs(float(x.get(k, 0)) - float(y.get(k, 0))))
    return dict(same_set=sa == sb, n_only_ref=len(sa - sb), n_only_port=len(sb - sa), max_log2=max_lr, max_cw=max_cw)


def b3(rows, ref, port, exp_md5=None):
    pr, pp = ref["pack_table"], port["pack_table"]
    if not os.path.exists(pp):
        return rows.add("B3", "inventory pack table vs reference ingredients", exp_md5 or "", "", "", "NA", f"missing {pp}")
    mr, mp = common.md5(pr), common.md5(pp)
    if mr == mp and (exp_md5 is None or mp == exp_md5):
        return rows.add("B3", "pack table md5", mr, mp, "identical", "E", "", delta=0)
    c = _species_table_check(pr, pp)
    ok = c["same_set"] and c["max_log2"] <= B3_LOG2 and c["max_cw"] <= B3_CW
    rows.add("B3", "pack table: max |log2 copies ratio| / max |d comp_weights| / species set", "0 / 0 / identical",
             f"{c['max_log2']:.3g} / {c['max_cw']:.3g} / {'identical' if c['same_set'] else 'differs'}",
             f"<= {B3_LOG2} / <= {B3_CW} / identical", "T" if ok else "FAIL",
             f"md5 {mp} != {mr}; only in ref {c['n_only_ref']}, only in port {c['n_only_port']}", delta=c["max_log2"])


def b4(rows, bdir, cr):
    pf = cr.get("port_files", {}) or {}
    ic = os.path.join(bdir, pf.get("inventory_compare", "inventory/compare_legacy.json"))
    if os.path.exists(ic) and ic.endswith(".json"):          # inventory stage comparison (wholecell.inventory.compare)
        J = json.load(open(ic))
        tabs = dict(J.get("tables", {}) or {}); tabs.update(J.get("ingredients", {}) or {})
        same = {k for k, v in tabs.items() if v.get("identical") or (v.get("content", {}) or {}).get("identical_content")
                or (v.get("content", {}) or {}).get("identical_content_any_row_order")}
        rest = sorted(set(tabs) - same)
        b3ok = bool((J.get("summary", {}) or {}).get("B3_all_pass"))
        lvl = "E" if not rest else ("T" if b3ok else "FAIL")
        return rows.add("B4", "intermediate inventory tables: identical content (any row order) / all", len(tabs), len(same),
                        "identical content; else B3 tolerances (change tables: documented differences only)", lvl,
                        (f"content differs: {rest} (B3 tolerances pass; check the documented explanation in "
                         f"{os.path.relpath(ic, bdir)} notes)" if rest else f"from {os.path.relpath(ic, bdir)}"), delta=len(tabs) - len(same))
    if os.path.exists(ic):
        T = pd.read_csv(ic, sep="\t")
        lc = next((c for c in ("level", "result", "status") if c in T.columns), None)
        if lc is not None:
            lv = T[lc].astype(str).str.upper()
            n_e = int(lv.isin(["E", "IDENTICAL", "TRUE"]).sum()); n_t = int(lv.isin(["T"]).sum()); n = len(T)
            lvl = "E" if n_e == n else ("T" if n_e + n_t == n else "FAIL")
            return rows.add("B4", "intermediate inventory tables (inventory stage comparison)", n, f"E {n_e}, T {n_t}",
                            "identical; else B3 tolerances per table", lvl, f"from {os.path.relpath(ic, bdir)}", delta="")
    pairs = cr.get("inventory_pairs", {}) or {}
    if not pairs:
        return rows.add("B4", "intermediate inventory tables", "", "", "", "NA",
                        "no inventory comparison table and no validate.compare_reference.inventory_pairs")
    n_e = n_t = 0; notes = []
    for rel, legacy in pairs.items():
        p = os.path.join(bdir, rel)
        if not os.path.exists(p):
            notes.append(f"missing {rel}"); continue
        if common.md5(p) == common.md5(legacy):
            n_e += 1; continue
        try:
            c = _species_table_check(legacy, p)
            if c["same_set"] and c["max_log2"] <= B3_LOG2 and c["max_cw"] <= B3_CW:
                n_t += 1; continue
        except Exception as e:  # noqa: BLE001 - tables without ing_id/copies columns
            notes.append(f"{rel}: {type(e).__name__}")
        notes.append(f"{rel} differs")
    n = len(pairs)
    rows.add("B4", "intermediate inventory tables (md5 / B3 tolerances)", n, f"E {n_e}, T {n_t}", "identical; else B3 tolerances",
             "E" if n_e == n else ("T" if n_e + n_t == n else "FAIL"), "; ".join(notes), delta="")


def b5(rows, ref, bdir, cr):
    pl = ref_path(bdir, (cr.get("port_files", {}) or {}).get("lamina", "preplaced/lamina_cell1.npz"))
    rl = ref.get("lamina")
    if not rl or not os.path.exists(pl):
        return rows.add("B5", "lamina rebuilt vs frozen lamina_cell1.npz", "", "", "", "NA", f"missing {pl}")
    A = np.load(rl); B = np.load(pl)
    same = set(A.files) == set(B.files) and all(np.array_equal(A[k], B[k]) for k in A.files)
    k = "types" if "types" in A.files else A.files[0]
    na, nb = len(A[k]), len(B[k]) if k in B.files else 0
    rows.add("B5", "lamina instance count", na, nb, f"+- {B5_REL * 100:.1f} %", _lvl(same, abs(nb - na) <= B5_REL * na),
             "arrays identical" if same else "arrays differ (positions not required)")


def _counts_table(cpath, libs_path, rnames):
    z = np.load(cpath); L = json.load(open(libs_path))["libs"]
    key = z["lib"].astype(np.int64) * 64 + z["region"].astype(np.int64)
    u, inv = np.unique(key, return_inverse=True); n = np.bincount(inv, z["n"].astype(np.float64))
    return {(L[k // 64], rnames[k % 64] if rnames else int(k % 64)): v for k, v in zip(u.tolist(), n.tolist())}, float(z["n"].sum())


def b6(rows, ref, port, rnames_ref, rnames_port):
    pr, pp = ref["counts"], port["counts"]
    if not os.path.exists(pp):
        return rows.add("B6", "counts.npz", "", "", "", "NA", f"missing {pp}")
    if common.md5(pr) == common.md5(pp):
        return rows.add("B6", "counts.npz md5", common.md5(pr), common.md5(pp), "identical", "E", "", delta=0)
    za, zb = np.load(pr), np.load(pp)
    if set(za.files) == set(zb.files) and all(np.array_equal(za[k], zb[k]) for k in za.files):
        return rows.add("B6", "counts.npz arrays", "identical", "identical", "identical", "E", "file bytes differ (compression)", delta=0)
    A, ta = _counts_table(pr, ref["libs"], rnames_ref); B, tb = _counts_table(pp, port["libs"], rnames_port)
    bad = 0
    for k in set(A) | set(B):
        n = A.get(k, 0.0); d = abs(B.get(k, 0.0) - n)
        if d > max(B6_PAIR_SIGMA * math.sqrt(n), B6_PAIR_REL * n): bad += 1
    ok = abs(tb - ta) <= B6_TOTAL_REL * ta and bad == 0
    rows.add("B6", "counts: total requested instances (pairs out of tolerance)", ta, tb,
             f"total +- {B6_TOTAL_REL * 100:.1f} %; pairs within max({B6_PAIR_SIGMA:g} sqrt n, {B6_PAIR_REL * 100:.1f} % n)",
             "T" if ok else "FAIL", f"{bad} (species, region) pairs out of tolerance of {len(set(A) | set(B))}")


def b7(rows, ref, port):
    keys = ["preplaced_types", "preplaced_pos", "preplaced_quat", "preplaced_off", "preplaced_libs", "mrna_off", "mrna_pts"]
    miss = [k for k in keys if not os.path.exists(port[k])]
    if miss:
        return rows.add("B7", "pre-pass index + mRNA bead index", "", "", "", "NA", f"missing {miss}")
    diff = [k for k in keys if common.md5(ref[k]) != common.md5(port[k])]
    if not diff:
        return rows.add("B7", "pre-pass index + mRNA index (7 files, md5)", "identical", "identical", "identical", "E", "", delta=0)
    pa, pb = np.diff(np.load(ref["preplaced_off"])), np.diff(np.load(port["preplaced_off"]))
    ma, mb = np.diff(np.load(ref["mrna_off"])), np.diff(np.load(port["mrna_off"]))
    ok = np.array_equal(pa, pb) and np.array_equal(ma, mb)
    rows.add("B7", "pre-pass / mRNA instance counts per tile", int(pa.sum()), int(pb.sum()), "per-tile counts identical",
             "T" if ok else "FAIL", f"files differing: {diff}")


def _pairs_placed(S, libs, rnames):
    pla = S["pla"]; nreg = S["nreg"]
    out = {}
    nz = np.nonzero(pla)
    for li, ri in zip(*nz):
        out[(libs[li], rnames[ri] if ri < len(rnames) else int(ri))] = float(pla[li, ri])
    return out


def b8_b9_b12(rows, cfg, ref, port, Sref, Sport, libs_ref, libs_port, rn_ref, rn_port, subset_note=""):
    # B8: instance count per tile
    ta = dict(zip(Sref["tile"].tolist(), Sref["ninst"].tolist())); tb = dict(zip(Sport["tile"].tolist(), Sport["ninst"].tolist()))
    tiles = sorted(set(ta) | set(tb))
    ndiff = sum(ta.get(t, -1) != tb.get(t, -1) for t in tiles)
    rows.add("B8", "packed instances per tile (tiles differing)", len(ta), len(tb), "identical for all tiles",
             "E" if ndiff == 0 else "FAIL", f"{ndiff} of {len(tiles)} tiles differ{subset_note}", delta=ndiff)
    # B9: placed copies per (species, region)
    A = _pairs_placed(Sref, libs_ref, rn_ref); B = _pairs_placed(Sport, libs_port, rn_port)
    keys = set(A) | set(B)
    same = all(A.get(k, 0.0) == B.get(k, 0.0) for k in keys)
    big = [k for k in keys if A.get(k, 0.0) >= B9_NMIN]; small = [k for k in keys if A.get(k, 0.0) < B9_NMIN]
    bad_big = [k for k in big if abs(B.get(k, 0.0) - A[k]) > max(B9_SIGMA * math.sqrt(A[k]), B9_REL * A[k])]
    ok_small = sum(abs(B.get(k, 0.0) - A.get(k, 0.0)) <= B9_SMALL_ABS for k in small)
    frac_small = ok_small / len(small) if small else 1.0
    lvl = _lvl(same, not bad_big and frac_small >= B9_SMALL_FRAC)
    rows.add("B9", "placed copies per (species, region): pairs n_v7>=10 out of tolerance / small pairs within +-3",
             f"{len(big)} / {len(small)}", f"{len(bad_big)} / {frac_small:.4f}",
             f"0 out of max({B9_SIGMA:g} sqrt n, {B9_REL * 100:g} % n) / >= {B9_SMALL_FRAC}", lvl,
             (f"{len(keys)} pairs; worst: {sorted(bad_big, key=lambda k: -abs(B.get(k, 0) - A[k]))[:3]}" if bad_big else f"{len(keys)} pairs")
             + subset_note, delta=len(bad_big))
    # B12: positions / orientations, sampled tiles
    common_t = [t for t in tiles if t in ta and t in tb]
    rng = np.random.default_rng(__import__("wholecell.seeds", fromlist=["seed_int"]).seed_int(cfg, "validate.compare.b12"))
    samp = sorted(rng.choice(common_t, min(B12_SAMPLE, len(common_t)), replace=False).tolist()) if common_t else []
    neq = 0
    for t in samp:
        za = np.load(os.path.join(ref["tiles_dir"], f"{t}.npz")); zb = np.load(os.path.join(port["tiles_dir"], f"{t}.npz"))
        if not all(np.array_equal(za[k], zb[k]) for k in ("types", "pos", "quat")):
            neq += 1
    rows.add("B12", "positions/orientations identical (sampled tiles differing)", len(samp), neq,
             f"identical arrays ({len(samp)} sampled tiles + B8 counts)", "E" if (neq == 0 and ndiff == 0) else "FAIL",
             "not required for T" + subset_note, delta=neq)


def b10_b11(rows, ref, port, exp_tot):
    A = _load_json(ref["summary"]); B = _load_json(port["summary"])
    if B is None:
        rows.add("B10", "placed mass per region (summary.json)", "", "", "", "NA", f"missing {port['summary']}")
        return rows.add("B11", "cell totals", "", "", "", "NA", f"missing {port['summary']}")
    ra, rb = A["regions"], B.get("regions", {})
    for r in ra:
        a, b = ra[r], rb.get(r)
        if b is None:
            rows.add(f"B10:{r}", "region frac / placed_pg", a["frac"], "", f"+- {B10_FRAC}", "FAIL", "region missing in port"); continue
        same = a["frac"] == b["frac"] and a["placed_pg"] == b["placed_pg"]
        rows.add(f"B10:{r}", "region placed mass fraction (placed_pg)", f"{a['frac']} ({a['placed_pg']})", f"{b['frac']} ({b['placed_pg']})",
                 f"frac +- {B10_FRAC}", _lvl(same, abs(b["frac"] - a["frac"]) <= B10_FRAC), "", delta=b["frac"] - a["frac"])
    for r in set(rb) - set(ra):
        rows.add(f"B10:{r}", "region placed mass fraction", "", rb[r]["frac"], f"+- {B10_FRAC}", "FAIL", "region not in reference")
    ia, ib = A["instances"], B.get("instances")
    exp_ok = all(A.get(k) == v for k, v in (exp_tot or {}).items())
    same = ia == ib and A["placed_pg"] == B.get("placed_pg") and A["mass_frac"] == B.get("mass_frac")
    within = ib is not None and abs(ib - ia) <= B11_INST_REL * ia and abs(B.get("mass_frac", -1) - A["mass_frac"]) <= B11_MASS_FRAC
    rows.add("B11", "cell totals: instances / placed_pg / mass_frac", f"{ia} / {A['placed_pg']} / {A['mass_frac']}",
             f"{ib} / {B.get('placed_pg')} / {B.get('mass_frac')}", f"instances +- {B11_INST_REL * 100:.1f} %; mass_frac +- {B11_MASS_FRAC}",
             _lvl(same, within), "" if exp_ok else "reference summary differs from the expected totals in config", delta=(ib or 0) - ia)


# ------------------------------------------------------------------ V rows
def v1(rows, ra, rb):
    A = _load_json(os.path.join(ra, "densities.json")); B = _load_json(os.path.join(rb, "densities.json"))
    if A is None or B is None:
        return rows.add("V1", "densities", "", "", "", "NA", "densities.json missing (ref or port)")
    ta = {r["region"]: r for r in A["rows"]}; tb = {r["region"]: r for r in B["rows"]}
    for r in ta:
        a = ta[r]["mg_per_ml"]; b = tb.get(r, {}).get("mg_per_ml")
        rows.add(f"V1:{r}", "density mg/ml (packed + pre-pass)", a, b, f"+- {V1_MGML}",
                 "FAIL" if b is None else _lvl(a == b, abs(b - a) <= V1_MGML), "" if b is not None else "region missing in port")
    a, b = A.get("nuc_cyto_ratio"), B.get("nuc_cyto_ratio")
    rows.add("V1:nuc_cyto_ratio", "nuclear / cytoplasmic density ratio", a, b, f"+- {V1_RATIO}",
             "NA" if a is None or b is None else _lvl(a == b, abs(b - a) <= V1_RATIO))


def v2(rows, ra, rb):
    A = _load_json(os.path.join(ra, "crowding.json")); B = _load_json(os.path.join(rb, "crowding.json"))
    if A is None or B is None:
        return rows.add("V2", "crowding", "", "", "", "NA", "crowding.json missing (ref or port)")
    sa = {r["region"]: r for r in A["summary"]}; sb = {r["region"]: r for r in B["summary"]}
    for r in sa:
        a = sa[r]["phi"]; b = sb.get(r, {}).get("phi")
        rows.add(f"V2:{r}", "excluded-volume fraction phi", a, b, f"+- {V2_PHI}", "FAIL" if b is None else _lvl(a == b, abs(b - a) <= V2_PHI))


def _set_scores(S, name):
    s = S.get("sets", {}).get(name)
    if not s:
        return None
    out = {}
    for mode in ("intended", "placed"):
        for m in SCORE_METRICS:
            v = (s.get(mode, {}) or {}).get(m) or {}
            for agg in ("unweighted", "mass_weighted"):
                out[(mode, m, agg)] = v.get(agg)
    return out


def v3_v4(rows, ra, rb):
    A = _load_json(os.path.join(ra, "loc_holdout.json")); B = _load_json(os.path.join(rb, "loc_holdout.json"))
    if A is None or B is None:
        rows.add("V3", "loc_compare agreement", "", "", "", "NA", "loc_holdout.json missing (ref or port)")
        return rows.add("V4", "loc_compare hold-out half", "", "", "", "NA", "loc_holdout.json missing (ref or port)")
    for sets, rid, tol in ((V3_SETS, "V3", V3_TOL), ((V4_SET,), "V4", V4_TOL)):
        for name in sets:
            a, b = _set_scores(A, name), _set_scores(B, name)
            if a is None or b is None:
                rows.add(f"{rid}:{name}", "agreement scores", "", "", f"+- {tol}", "NA", "set missing"); continue
            d = [abs(b[k] - a[k]) for k in a if a[k] is not None and b.get(k) is not None]
            md = max(d) if d else float("nan")
            same = all(a[k] == b.get(k) for k in a)
            ka = ("placed", "agree", "mass_weighted")
            rows.add(f"{rid}:{name}", "agree/agree_tvd/tophit x intended/placed x unweighted/mass-weighted: max |delta| "
                     "(shown: placed agree, mass-weighted)", a[ka], b[ka], f"each +- {tol}", _lvl(same, md <= tol),
                     f"max |delta| over 12 scores = {md:.4g}", delta=md)
    # per-gene placed shares
    pa, pb = os.path.join(ra, "loc_compare.tsv"), os.path.join(rb, "loc_compare.tsv")
    if not (os.path.exists(pa) and os.path.exists(pb)):
        return rows.add("V3:genes", "per-gene placed compartment shares", "", "", "", "NA", "loc_compare.tsv missing")
    TA = pd.read_csv(pa, sep="\t", index_col=0); TB = pd.read_csv(pb, sep="\t", index_col=0)
    cols = [c for c in TA.columns if c.startswith("pl_")]
    genes = sorted(set(TA.index) | set(TB.index))
    XA = TA.reindex(genes)[cols].fillna(0.0).values; XB = TB.reindex(genes)[cols].reindex(columns=cols).fillna(0.0).values
    present = np.isin(genes, TA.index) & np.isin(genes, TB.index)
    tvd = 0.5 * np.abs(XA - XB).sum(1)
    ok = (tvd <= V3_GENE_TVD) & present
    fr = float(ok.mean()) if len(ok) else 1.0
    rows.add("V3:genes", "genes with placed-share TVD <= 0.02 (fraction)", 1.0, round(fr, 5), f">= {V3_GENE_FRAC}",
             _lvl(bool(np.all(tvd == 0) and present.all()), fr >= V3_GENE_FRAC),
             f"{len(genes)} genes, {int((~present).sum())} present on one side only, max TVD {float(tvd.max()) if len(tvd) else 0:.4g}",
             delta=fr - 1.0)


def v5(rows, ra, rb, exp_md5):
    A = _load_json(os.path.join(ra, "chromatin_frozen.json")); B = _load_json(os.path.join(rb, "chromatin_frozen.json"))
    if B is None:
        return rows.add("V5", "chromatin metrics (frozen)", "", "", "", "NA", "chromatin_frozen.json missing in port")
    ok_sum = B["prepass"]["md5_actual"] == (exp_md5 or B["prepass"]["md5_expected"])
    same = A is not None and A["prepass"]["md5_actual"] == B["prepass"]["md5_actual"] and \
        {k: v["md5"] for k, v in A["files"].items()} == {k: v["md5"] for k, v in B["files"].items()} and A["files"] == B["files"]
    rows.add("V5", "chromatin pre-pass checksum; frozen metrics copied verbatim", A["prepass"]["md5_actual"] if A else "",
             B["prepass"]["md5_actual"], "checksum identical; values verbatim", "E" if (ok_sum and same) else "FAIL",
             "; ".join(f"{k}: {v.get('v8')}" for k, v in (B.get("headline") or {}).items())[:300], delta="")


def v6(rows, ra, rb):
    A = _load_json(os.path.join(ra, "placement.json")); B = _load_json(os.path.join(rb, "placement.json"))
    if A is None or B is None:
        return rows.add("V6", "placement integrity", "", "", "", "NA", "placement.json missing (ref or port)")
    note = "" if A.get("tiles_md5") == B.get("tiles_md5") else "DIFFERENT TILE SAMPLES"
    if A.get("n_tiles_sampled", 0) < V6_MIN_TILES or B.get("n_tiles_sampled", 0) < V6_MIN_TILES:
        note += f" fewer than {V6_MIN_TILES} tiles sampled"
    for key, lab in (("clash_rate_per_1e6", "a clash rate per 1e6 instances"),
                     ("containment_rate_per_1e6", "b containment violations per 1e6 instances"),
                     ("orientation_error_rate_per_1e6", "c membrane orientation sign errors per 1e6 membrane instances")):
        a, b = A.get(key), B.get(key)
        if a is None or b is None:
            rows.add(f"V6:{key.split('_')[0]}", lab, a, b, "", "NA", "rate undefined"); continue
        tol = max(V6_REL * abs(a), V6_ABS)
        lvl = _lvl(a == b and not note, abs(b - a) <= tol and "DIFFERENT" not in note)
        rows.add(f"V6:{key.split('_')[0]}", lab, round(a, 4), round(b, 4), f"max(10 % rel, 1 per 1e6) = {tol:.3g}", lvl, note.strip())


def v7(rows, ra, rb):
    A = _load_json(os.path.join(ra, "ribosome_nn.json")); B = _load_json(os.path.join(rb, "ribosome_nn.json"))
    if A is None or B is None:
        return rows.add("V7", "ribosome NN median", "", "", "", "NA", "ribosome_nn.json missing (ref or port)")

    def pick(J):
        rr = {r["set"]: r for r in J["rows"]}
        k = next((s for s in rr if s.startswith("model pre-pass + packed") and "thinned" not in s), "model (all)")
        return k, rr[k]["nn_median"]
    ka, a = pick(A); kb, b = pick(B)
    rows.add("V7", f"ribosome NN median (nm), row '{ka}'", a, b, f"+- {V7_NN} nm", _lvl(a == b, abs(b - a) <= V7_NN),
             "" if ka == kb else f"port row '{kb}'")


# ------------------------------------------------------------------ driver
def verdict(rows):
    must = [r for r in rows if r["class"] == "must"]; rep = [r for r in rows if r["class"] == "report"]
    must_fail = [r["id"] for r in must if r["level"] == "FAIL"]
    must_na = [r["id"] for r in must if r["level"] == "NA"]
    rep_fail = [r["id"] for r in rep if r["level"] == "FAIL"]
    if must_fail or must_na:
        v = "NOT REPRODUCED"
    elif rep_fail:
        v = "REPRODUCED WITH DEVIATIONS"
    else:
        v = "REPRODUCED"
    lv = pd.Series([r["level"] for r in rows]).value_counts().to_dict()
    return dict(verdict=v, must_failed=must_fail, must_not_evaluated=must_na, report_failed=rep_fail, levels=lv,
                n_rows=len(rows), all_exact=all(r["level"] == "E" for r in rows if r["id"] != "V8"))


def run(cfg, bdir, **kw):
    cr = (cfg.get("validate", {}) or {}).get("compare_reference", {}) or {}
    if not cr:
        raise RuntimeError("no validate.compare_reference in the config")
    nproc = common.nprocs(cfg, kw)
    ref_name = cr.get("name", "ref")
    ref_vdir = kw.get("ref_validate_dir") or cr["validate_dir"]
    ref_bdir = os.path.dirname(ref_vdir)
    ref = common.Inputs(cfg, ref_bdir, ref=ref_name)
    pkw = {k[5:]: v for k, v in kw.items() if k.startswith("port_")}   # port_pack_dir=, port_tiles_dir=, ...
    port = common.Inputs(cfg, bdir, **pkw)
    port_vdir = os.path.join(bdir, "validate")
    out = common.outdir(bdir, f"compare_{ref_name}")
    from . import regions as R_
    rn_ref = R_.names(R_.load_regions(cfg, ref))
    try:
        rn_port = R_.names(R_.load_regions(cfg, port))
    except RuntimeError:
        rn_port = rn_ref
    rows = Rows()
    gv = kw.get("geometry_validation") or os.path.join(bdir, (cr.get("port_files", {}) or {}).get("geometry_validation", "geometry/validation.json"))
    if not os.path.exists(gv):
        pk = (cfg.get("geometry", {}) or {}).get("package") or os.path.join(bdir, "geometry", "package")
        if pk and os.path.exists(os.path.join(pk, "validation.json")):
            gv = os.path.join(pk, "validation.json")
    b1(rows, gv)
    b2(rows, ref, port)
    b3(rows, ref, port, (cr.get("expected_md5", {}) or {}).get("pack_table"))
    b4(rows, bdir, cr)
    b5(rows, ref, bdir, cr)
    b6(rows, ref, port, rn_ref, rn_port)
    b7(rows, ref, port)
    # B8/B9/B12 from the job tables of every tile (cached per side)
    pfiles = port.tile_files()
    rfiles = ref.tile_files()
    subset_note = ""
    if pfiles and len(pfiles) < len(rfiles):
        ids = {os.path.basename(f) for f in pfiles}
        rfiles = [f for f in rfiles if os.path.basename(f) in ids]
        subset_note = f" (port has {len(pfiles)} of {ref.n_tiles()} tiles: compared on that subset)"
        Sref = common.tile_sums(ref, None, nproc=nproc, files=rfiles)
    else:
        Sref = common.tile_sums(ref, os.path.join(ref_vdir, "cache", "tile_sums.npz"), nproc=nproc)
    if pfiles:
        Sport = common.tile_sums(port, None if subset_note else os.path.join(port_vdir, "cache", "tile_sums.npz"), nproc=nproc, files=pfiles)
        libs_ref = json.load(open(ref["libs"]))["libs"]; libs_port = json.load(open(port["libs"]))["libs"]
        b8_b9_b12(rows, cfg, ref, port, Sref, Sport, libs_ref, libs_port, rn_ref, rn_port, subset_note)
    else:
        for i, m in (("B8", "instances per tile"), ("B9", "placed copies per (species, region)"), ("B12", "positions")):
            rows.add(i, m, "", "", "", "NA", f"no port tiles in {port['tiles_dir']}")
    b10_b11(rows, ref, port, cr.get("expected_totals"))
    v1(rows, ref_vdir, port_vdir)
    v2(rows, ref_vdir, port_vdir)
    v3_v4(rows, ref_vdir, port_vdir)
    v5(rows, ref_vdir, port_vdir, (cfg.get("validate", {}).get("chromatin", {}) or {}).get("expected_md5"))
    v6(rows, ref_vdir, port_vdir)
    v7(rows, ref_vdir, port_vdir)
    rows.add("V8", "Tier-2 imaging and fluorescence panels", "", "", "", "NA", "not ported in phase 1 (stated as not evaluated)")
    T = pd.DataFrame(rows, columns=["id", "metric", "v7", "port", "delta", "tolerance", "level", "class", "note"])
    common.write_df(os.path.join(out, "criteria_results.tsv"), T, overwrite=True, index=False)
    V = verdict(rows)
    V.update(reference=dict(name=ref_name, validate_dir=ref_vdir, files=ref.paths), port=dict(build=bdir, files=port.paths),
             criteria_doc=cr.get("criteria_doc"), criteria_doc_md5=_md5_or_none(cr.get("criteria_doc")), subset=subset_note.strip())
    common.write_json(os.path.join(out, "verdict.json"), V, overwrite=True)
    print(T.to_string(index=False, max_colwidth=60))
    print(json.dumps({k: V[k] for k in ("verdict", "must_failed", "must_not_evaluated", "report_failed", "levels")}, indent=1))
    return V
