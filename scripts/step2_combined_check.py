# new in wholecell (Step-2 combined build, 2026-11-18)
"""Step-2 combined build (B6 + B6b, B3 [Jurkat], B4, B1, B2, B5 in ONE build per cell): build-level verdict.

    python scripts/step2_combined_check.py --cell jurkat|yeast --build <bdir> [--complete 0|1] [--out <dir>]
                                           [--nproc 3] [--quick] [--reuse-from <dir>] [--accepted <json> ...]

Writes <out>/verdict.json (default <bdir>/validate/compare_step2/) with one section per step plus the combined suite
no_regression, and <out>/regressions.tsv. Nothing is written into reference builds or packages.

Criteria (build-level items that are measurable on a full build; accepted deviations included):
  B6   docs/unification_criteria/B6_geometry.md Part B (B1, B3, B4; B2 superseded) + B6_accepted_deviations_2026-11-17.json
  B6b  docs/unification_criteria/B6b_crista_rims.md (build-level EXPECTED, last paragraph) + addendum 2026-11-17
  B3   docs/unification_criteria/B3_library.md Part B (B1, B3, B4; L1-L7 via the library) + addendum 2026-11-17 (L5a)
  B4   docs/unification_criteria/B4_preplaced.md (J1-J8 / Y1-Y6, Part B) + addenda 2026-11-17 (J5b per class),
       2026-11-18 J6 (4.5 nm TM tolerance, rim density) and 2026-11-18 polysomes (J3b / J4 accepted model limitation)
  B1   docs/unification_criteria/B1_packer.md Part B (K1-K9) + addendum 2026-11-18 (Part B applies at the combined build).
       DECISION pending/applied: COMMON_DEFAULTS = {} (v7 values), yeast keeps preset legacy_x5 (no C3 pin)
  B2   docs/unification_criteria/B2_organisation.md §4 (O0, O2, O3, O5; O1 is vs the B1 build -> NA here; O4 superseded)
  B5   docs/unification_criteria/B5_localisation.md §4 (L1, L3, L4, L5; L2 / L6 superseded; L7 = completeness) +
       addendum 2026-11-18 (C1 gate on independent references only)

Reference build of the combined no_regression: Jurkat = the accepted B3 build (geometry em_cell1_v8, library build mode,
frozen v7 pre-placed layers, B2/B5 off); yeast = the port s0 (+ PATCH 4; suite in yeast/_dev_unify/a8p2/). Each step's
own suite item is superseded by ONE combined no_regression in which a regression counts as EXPECTED only when a step's
criteria list it (the step is named), else UNEXPECTED.

Item status: PASS | FAIL | NA (cannot be judged here: reason given; never counted as passed) | REPORT (not gated) |
ACCEPTED_DEVIATION (fails the frozen bound; accepted by a dated decision cited in `reason`; not counted as failing).
Step status: ERROR (the section crashed) | FAIL (any FAIL) | PASS (>= 1 PASS, no FAIL) | NA (nothing judged).
Overall: INCOMPLETE (build not complete / suite missing) | ACCEPTED (every step PASS or NA-free of FAIL and
no_regression PASS) | NOT ACCEPTED. NA items are listed separately in `na_items`.
CPU: --nproc processes (default 3), NUMBA / OMP threads 3.

Fixes 2026-11-19 (diagnosis of the first Jurkat combined verdict; re-score with --out .../compare_step2_v2 --reuse-from
.../compare_step2): B1.K8 compared pack_done.json `switches` (Switches.as_dict) with the whole effective `pack:` section
(always unequal) -> compare with switches.load(effective(...)), Jurkat: v7 values except organisation; B4.J6 TM / soluble
items ignored the float32 export rounding that the J6 addendum §3 accepts -> eps 1e-3 nm (b4_j6_atp_check.py fields);
B1.K6 was REPORT and its mass-consistency rule (B1 EXPECTED "C11 densities: per K6") was not applied in the combined
no_regression -> gated + C11 rows classified "B1 K6" (the placed deficit stays judged by K2); the B5 revision's
documented L2 deviation (OpenCell intended tophit -0.014) is EXPECTED.

--accepted <json> (2026-11-19; e.g. docs/unification_criteria/step2_combined_accepted_deviations_2026-11-19.json): dated
build-level accepted deviations. An entry applies only to the builds it lists and only while the value is not worse than
its `limit`. kind "item": a FAIL step item (step, item, key = failing sub-metric / region) becomes ACCEPTED_DEVIATION when
EVERY failing key of it is covered; kind "regression": an UNEXPECTED no_regression row (check, metric) moves to
`regressions_accepted` (class from the entry: KNOWN / ACCEPTED / EXPECTED). no_regression with only accepted rows left is
ACCEPTED_DEVIATION. Verdict: ACCEPTED if every remaining failure is an accepted deviation.
"""
import argparse
import copy
import csv
import glob
import json
import os
import re
import sys
import time
import traceback

for _k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_k, "3")
os.environ.setdefault("PYTHONHASHSEED", "0")
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path[:0] = [os.path.join(REPO, "scripts", "dev"), os.path.join(REPO, "scripts")]
from wholecell import config as wcfg  # noqa: E402

SITE = wcfg.load_site()
RES, GEOM = SITE["results_root"], SITE["geom_root"]
REF = {"jurkat": f"{RES}/jurkat/20261117-bc121a6-s20261115-b3", "yeast": f"{RES}/yeast/20261116-aabe7b5-s0-x5repro"}
REF_SUITE = {"jurkat": f"{REF['jurkat']}/validate/suite.json",
             "yeast": f"{RES}/yeast/_dev_unify/a8p2/20261116-aabe7b5-s0-x5repro/suite.json"}
YEAST_SEED_SUITES = sorted(glob.glob(f"{RES}/yeast/_dev_unify/b6/seedspread_suite/*/suite.json"))
V7REPRO2 = f"{RES}/jurkat/20261115-ec0b403-s20261115-v7repro2"          # J5b reference tiles (B4 verdict)
PKG = {"jurkat": f"{GEOM}/jurkat/em_cell1_v8r", "yeast": f"{GEOM}/yeast/analytic_x5_s0_f5f3a5c723b2_memhalf1"}
PKG_REF_JURKAT = f"{GEOM}/jurkat/em_cell1_v8"                          # package of the Jurkat reference build
JURKAT_V8R_SHA = "6c86b1688722f4320520abbb46a66fb7e9e97c95576e0186ebf3909b9c562fa1"
B6_ACCEPTED = os.path.join(REPO, "docs/unification_criteria/B6_accepted_deviations_2026-11-17.json")
AV = 6.02214076e23
# B4 Y1 in the combined build (2026-11-19): library order of the PATCH-4 reference pre-pass (= s0 / port lib_ids), and the
# yeast serial builders that draw counts or placements from the inventory's localisation shares / pools (B5 changes them:
# chromatin_sites = TF / silencing nucleoplasm pools; actin_patches = patch proteins (B5 builder_anchored); polysomes =
# ER / mitochondrial tethering of the translated mRNAs, then rac, mrna, mrnp on the polysome chains)
YREF_LIB_IDS = f"{REF['yeast']}/engine/work/yc3d/lib_ids.npy"
Y1_LOCALISATION_BUILDERS = ("chromatin_sites", "actin_patches", "polysomes", "rac", "mrna", "mrnp")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def jload(p):
    return json.load(open(p)) if p and os.path.exists(p) else None


def _f(x):
    try:
        x = float(x)
        return None if x != x else x
    except (TypeError, ValueError):
        return None


def item(id_, status, value=None, rule="", reason="", source=""):
    assert status in ("PASS", "FAIL", "NA", "REPORT", "ACCEPTED_DEVIATION"), status
    return dict(id=id_, status=status, value=value, rule=rule, reason=reason, source=source)


def judged(id_, ok, value, rule, source="", reason=""):
    if ok is None:
        return item(id_, "NA", value, rule, reason or "value missing", source)
    return item(id_, "PASS" if ok else "FAIL", value, rule, reason, source)


def suite_rows(p):
    J = jload(p)
    return {(r["check"], r["metric"]): r for r in J["metrics"]} if J else {}


def sval(S, check, metric):
    r = S.get((check, metric))
    return _f(r.get("value")) if r else None


# ------------------------------------------------------------------------------------------------ accepted deviations
def load_accepted(paths, bdir):
    """entries of --accepted files that apply to this build (by build directory name)"""
    out = []
    for p in paths or []:
        J = json.load(open(p))
        if os.path.basename(bdir) not in J.get("builds", []):
            log("accepted", p, "does not list this build: ignored")
            continue
        for e in J["accepted"]:
            out.append(dict(e, _file=os.path.relpath(os.path.abspath(p), REPO), _source=J.get("decision_source", "")))
    return out


def _within(e, v):
    v = _f(v)
    return v is not None and (v >= e["limit"] if e["direction"] == "up" else v <= e["limit"])


def _acc_note(e):
    return f"{e['class']}: {e['cause']} [{e['decision']}; {e['_file']}]" + (f"; fix deferred to {e['fix_deferred_to']}" if e.get("fix_deferred_to") else "")


def _failing_keys(it):
    """{key: value} of the failing parts of a FAIL item (sub-metric lists with `ok`, B1.K2 region dict); None = not keyed"""
    v = it.get("value")
    if isinstance(v, list) and v and all(isinstance(x, dict) and "ok" in x for x in v):
        return {x["metric"]: x["value"] for x in v if x["ok"] is False}
    if isinstance(v, dict) and "below_0_98_and_below_reference" in v:
        return {k: x["build"] for k, x in v["below_0_98_and_below_reference"].items()}
    return None


def apply_accepted_items(steps, acc):
    for name, s in steps.items():
        for it in s.get("items", []):
            if it["status"] != "FAIL":
                continue
            ent = {e["key"]: e for e in acc if e["kind"] == "item" and e["step"] == name and e["item"] == it["id"]}
            keys = _failing_keys(it)
            if not ent or not keys:
                continue
            cov = {k: ent[k] for k, v in keys.items() if k in ent and _within(ent[k], v)}
            if set(cov) == set(keys):
                it["status"] = "ACCEPTED_DEVIATION"
                it["reason"] = "; ".join(f"{k} = {_f(keys[k]):.4g}: {_acc_note(e)}" for k, e in cov.items())
            else:
                it["reason"] = (it.get("reason") or "") + f" [not covered by --accepted: {sorted(set(keys) - set(cov))}]"


def apply_accepted_regressions(nr, acc):
    ent = {(e["check"], e["metric"]): e for e in acc if e["kind"] == "regression"}
    keep, moved = [], []
    for r in nr.get("regressions_unexpected", []):
        e = ent.get((r["check"], r["metric"]))
        if e and _within(e, r["value"]):
            moved.append(dict(r, step=e.get("step") or r.get("step"), accepted_class=e["class"], explained=_acc_note(e)))
        else:
            keep.append(r)
    nr["regressions_unexpected"], nr["regressions_accepted"] = keep, moved
    if nr.get("status") == "FAIL" and not keep:
        nr["status"] = "ACCEPTED_DEVIATION"


# ------------------------------------------------------------------------------------------------ combined suite
def combined_no_regression(cell, bdir, ctx):
    """suite rows of the build vs the reference's suite.json (wholecell.validate.suite.no_regression), each regression
    classified EXPECTED (step named) or UNEXPECTED"""
    from wholecell.validate.suite import no_regression
    S = jload(os.path.join(bdir, "validate", "suite.json"))
    R = jload(REF_SUITE[cell])
    if R:
        ctx["suite_ref"] = {(r["check"], r["metric"]): r for r in R["metrics"]}
    if S is None or R is None:
        return dict(status="NA", reason=f"suite.json missing (build: {S is not None}, reference {REF_SUITE[cell]}: {R is not None})")
    table, _ = no_regression([copy.deepcopy(r) for r in S["metrics"]], R["metrics"], REF_SUITE[cell])
    ctx["suite_new"] = {(r["check"], r["metric"]): r for r in S["metrics"]}
    ctx["suite_ref"] = {(r["check"], r["metric"]): r for r in R["metrics"]}
    band = {}
    if cell == "yeast":
        for sp in YEAST_SEED_SUITES:
            for r in json.load(open(sp))["metrics"]:
                x = _f(r.get("value"))
                if x is not None:
                    band.setdefault((r["check"], r["metric"]), []).append(x)
    acc = {(a["check"], a["metric"]): a for a in json.load(open(B6_ACCEPTED))["accepted"]} if cell == "jurkat" else {}
    regs = [r for r in table if r["status"] == "regress"]
    expected, unexpected = [], []
    for r in regs:
        why = classify(cell, r, band, acc, ctx)
        rec = dict(check=r["check"], metric=r["metric"], value=r["value"], reference=r["reference_build_value"],
                   tolerance=r["tolerance"], direction=r["direction"], step=why[0] if why else None,
                   explained=why[1] if why else None)
        (expected if why else unexpected).append(rec)
    counts = {}
    for r in table:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    return dict(status="PASS" if not unexpected else "FAIL", reference_suite=REF_SUITE[cell], reference_build=REF[cell],
                counts=counts, regressions_unexpected=unexpected, regressions_expected=expected,
                rule="no regression (plan §4 tolerances) outside the EXPECTED lists of the steps' criteria")


GATED_B5 = {"jurkat": ("opencell", "hpa_jurkat_holdout_half"), "yeast": ("SWAT_native", "SWAT_NOP1", "hyperLOPIT")}
# B5 L2 rows of GATED sets that the B5 revision reported outside ±0.01 before B5 was accepted (B5_localisation_revision_
# 2026-11-18.md §2 "L2 ... Outside the band: OpenCell tophit −0.014", intended level of the dry run dryrun_7b9880f:
# 0.8858 vs 0.8997; B5_localisation_addendum_2026-11-18.md: "B5 accepted (Pedro, 2026-11-18)" on that revision).
# Only the documented metric; mass-weighted / placed variants are not listed there and stay UNEXPECTED (decision).
B5_DOCUMENTED_L2 = {"jurkat": {"opencell.intended.tophit.unweighted": (
    -0.014, "documented before acceptance: B5 revision 2026-11-18 §2 L2 'OpenCell tophit -0.014' (dry run 0.8858 vs 0.8997); "
            "B5 accepted by Pedro 2026-11-18 (addendum)")}}
B4_JURKAT_KEYS = ("perinuclear", "density.cytosol", "density.matrix", "prepass_mrna", "prepass_prepass", "nn_median",
                  "c12", "nup", "npc")      # b4_preplaced_check.EXPECTED["jurkat"]: judged by J1-J7 instead


def classify(cell, r, band, acc, ctx):
    """(step, explanation) when a criteria file lists this regression as expected, else None"""
    c, m = r["check"], str(r["metric"])
    key = (c, m); full = f"{c}.{m}".lower()
    val, ref, tol = _f(r["value"]), _f(r["reference_build_value"]), _f(r["tolerance"]) or 0.0
    if key in acc:
        return "B6", f"accepted deviation: {acc[key]['decision']}"
    # B5: C1 rows of non-gated sets (addendum 2026-11-18), C2-C4 (L6 scope C5-C17; C3 judged by L3)
    if c == "C1":
        s = m.split(".")[0]
        if s not in GATED_B5[cell]:
            return "B5", "C1 set reported, not gated (B5 addendum 2026-11-18 / CYCLoPs = model input)"
        doc = B5_DOCUMENTED_L2.get(cell, {}).get(m)
        if doc is not None and val is not None and ref is not None and val - ref >= doc[0] - 0.002:
            return "B5", doc[1]
    if c in ("C2", "C3", "C4"):
        return "B5", "outside the B5 L6 no-regression scope (C5-C17); localisation judged by L1 / L3"
    # B1: C10 placed fractions judged by K1 / K2
    if c == "C10" and m in ("placed_mass_fraction", "placed_copies_fraction", "species_below_50pct_placed"):
        return "B1", "C10 placed fractions judged by K1 / K2 (B1 EXPECTED)"
    # B2: yeast C15 judged by O2
    if c == "C15":
        return "B2", "C15 judged by O2 (B2 EXPECTED)"
    if cell == "jurkat":
        if c == "C11" and ("mito_ims" in full or "ims" in m.lower().split(".")[-1]) and val is not None and ref:
            if abs(val - ref) <= 0.025 * abs(ref) + tol:
                return "B6b", "IMS density change <= 2.5 % (B6b R6 volume moved)"
        if any(k in full for k in B4_JURKAT_KEYS):
            if "matrix" in full and val is not None and ref and abs(val - ref) <= 0.025 * abs(ref) + 2.0:
                return "B4/B6b", "matrix density: new ATP synthase / mtDNA content and the B6b volume change (B4 EXPECTED ±2 mg/ml or mass-consistent)"
            if "matrix" not in full:
                return "B4", "B4 EXPECTED (judged by J1-J7)"
        if c == "C11" and m.startswith("density."):
            d = ctx.get("b5_density_intended")
            if d is None:
                d = ctx["b5_density_intended"] = _intended_density(cell, ctx)
            di = d.get(m.split(".", 1)[1]) if isinstance(d, dict) else None
            if di is not None and val is not None and ref is not None and abs((val - ref) - di) <= 2.0:
                return "B5", f"density change explained by the requested-mass change {di:+.2f} mg/ml (B5 L6)"
            # B1 EXPECTED "C11 densities: per K6": mass-consistent with the placed-mass change -> the density row itself is
            # explained; the placed-mass deficit behind it is NOT exempted: it is judged by B1.K2 (and listed there)
            dp = _placed_density_change(cell, ctx).get(m.split(".", 1)[1])
            if dp is not None and val is not None and ref is not None and abs((val - ref) - dp) <= 2.0:
                return "B1 K6", (f"mass-consistent with the placed packed-mass change {dp:+.2f} mg/ml (B1 K6); the placed deficit "
                                 f"is judged by B1.K2 (intended change {di:+.2f})" if di is not None else
                                 f"mass-consistent with the placed packed-mass change {dp:+.2f} mg/ml (B1 K6); deficit judged by B1.K2")
    else:
        b_ = band.get(key)
        if b_ and val is not None:
            b_ = b_ + ([ref] if ref is not None else [])
            if min(b_) - tol <= val <= max(b_) + tol:
                return "B6", f"within the port seed spread s0/s1/s2 [{min(b_):.6g}, {max(b_):.6g}] +- {tol} (B6 yeast EXPECTED)"
        if c == "C11" and m.startswith("density."):
            d = ctx.get("b5_density_intended")
            if d is None:
                d = ctx["b5_density_intended"] = _intended_density(cell, ctx)
            di = d.get(m.split(".", 1)[1]) if isinstance(d, dict) else None
            if di is not None and val is not None and ref is not None and abs((val - ref) - di) <= 2.0:
                return "B5", f"density change explained by the requested-mass change {di:+.2f} mg/ml (B5 L6)"
    return None


def _placed_density_change(cell, ctx):
    """B1 K6: {region: (placed packed mass, build - reference) / region volume} in mg/ml from the K1/K2 tile sums; {} if
    not computable (yeast: the reference is the port layout)"""
    if "placed_density_change" in ctx:
        return ctx["placed_density_change"]
    res = {}
    try:
        if cell == "jurkat":
            for tag, b in (("build", ctx["bdir"]), ("reference", REF[cell])):
                ctx.setdefault("c10", {})[tag] = ctx.get("c10", {}).get(tag) or _c10_whole_cell(cell, b, ctx, tag)
            vol = {r["region"]: r["volume_um3"] for r in _densities_json(ctx["bdir"]).get("rows") or []}
            nb, rf = ctx["c10"]["build"]["regions"], ctx["c10"]["reference"]["regions"]
            res = {k: (nb[k]["placed_pg"] - (rf.get(k) or {}).get("placed_pg", 0.0)) / vol[k] * 1e3 for k in nb if vol.get(k)}
    except Exception as e:
        ctx.setdefault("notes", []).append(f"B1 K6 placed-mass density change not computable: {type(e).__name__}: {e}")
    ctx["placed_density_change"] = res
    return res


def _intended_density(cell, ctx):
    """B5 L6: requested-mass density change per region (build - reference), b5_localisation_check helper; {} on failure"""
    try:
        import b5_localisation_check as B5
        return B5.requested_density_change(cell, ctx["bdir"], REF[cell], os.path.join(ctx["out"], "b5_density"), ctx["nproc"])
    except Exception as e:
        ctx.setdefault("notes", []).append(f"B5 requested-mass density change not computable: {type(e).__name__}: {e}")
        return {}


# ------------------------------------------------------------------------------------------------ B6 / B6b
def sec_b6(cell, bdir, ctx):
    items = []
    pj = jload(os.path.join(bdir, "geometry", "package.json")) or {}
    gv = jload(os.path.join(bdir, "geometry", "validation.json")) or {}
    pre = jload(f"{RES}/{cell}/_dev_unify/b6/geometry_check.json") or {}
    want = os.path.basename(PKG[cell])
    pid = pj.get("package_id")
    items.append(judged("B6.B1", bool(pid and want in str(pid) or want in str(pj.get("package_dir", ""))) and bool(gv.get("pass")),
                        dict(package_id=pid, package_dir=pj.get("package_dir"), validation_pass=gv.get("pass")),
                        f"geometry/package.json names {want}; geometry validation pass", "geometry/package.json, validation.json"))
    items.append(judged("B6.pre_build", pre.get("pass"), pre.get("pass"), "pre-build B6 checks passed",
                        f"{RES}/{cell}/_dev_unify/b6/geometry_check.json"))
    meas = (jload(os.path.join(bdir, "geometry", "measures.json")) or {}).get("region_volumes_mc")
    if cell == "jurkat":
        man = jload(os.path.join(PKG[cell], "manifest.json"))
        import hashlib
        sha = hashlib.sha256(open(os.path.join(PKG[cell], "manifest.json"), "rb").read()).hexdigest() if man else None
        items.append(judged("B6b.package", sha == JURKAT_V8R_SHA, sha, "em_cell1_v8r manifest sha256 = 6c86b168... (B6b R2)",
                            f"{PKG[cell]}/manifest.json"))
        refm = (jload(os.path.join(REF[cell], "geometry", "measures.json")) or {}).get("region_volumes_mc")
        pv8 = jload(os.path.join(PKG_REF_JURKAT, "validation.json"))["checks"]["volumes_mc"]["volume_um3"]
        pv8r = jload(os.path.join(PKG[cell], "validation.json"))["checks"]["volumes_mc"]["volume_um3"]
        if meas and refm:
            v, r, ci = meas["volume_um3"], refm["volume_um3"], meas["ci95_um3"]
            rho = float(meas.get("n", 2_000_000)) / float(meas.get("box_um3") or 1.0)      # MC points per um3
            ptol = lambda d, e: 3.0 * np.sqrt((abs(d) + abs(e)) / rho) + 2.0 / rho           # paired-delta noise (3 sd)
            dev = {}
            for k in v:
                exp = pv8r.get(k, 0.0) - pv8.get(k, 0.0)
                d = v[k] - r.get(k, 0.0)
                dev[k] = dict(delta=round(d, 4), expected=round(exp, 4), tol=round(ptol(d, exp), 4), ok=bool(abs(d - exp) <= ptol(d, exp)))
            items.append(judged("B6.B3/B6b.R4/B4.J1b volumes", all(d["ok"] for d in dev.values()), dev,
                                "every MC label volume change vs the reference build (em_cell1_v8; same seed stream, paired points) = "
                                "the package change v8 -> v8r of the validator MC (B4 J1b pores, B6b rims) within 3 sd of the "
                                "paired-difference noise", "geometry/measures.json (build, reference), validation.json of both packages"))
            per = v.get("perinuclear_space")
            items.append(judged("B6.J4 perinuclear", per is not None and 7.7 - 0.43 - 2 * ci.get("perinuclear_space", 0) <= per <= 8.5 - 0.43 + 2 * ci.get("perinuclear_space", 0),
                                per, "B6 J4 window 7.7-8.5 um3 shifted by the B4 EXPECTED pore openings -0.43 um3 (1,696 x 2.54e-4), +-2 CI",
                                "geometry/measures.json"))
            dm = v["mito_matrix"] - r["mito_matrix"]
            items.append(judged("B6b.R6 matrix", 0.10 - ptol(dm, 0) <= dm <= 0.25 + ptol(dm, 0), round(dm, 4),
                                "mito_matrix up by 0.10-0.25 um3 vs the reference (paired points; IMS -> matrix at the rims)", "geometry/measures.json"))
        else:
            items.append(item("B6.B3 volumes", "NA", None, "MC volumes vs reference", "geometry/measures.json missing (build or reference)"))
        lg = os.path.join(bdir, "logs", "pack_index.log")
        txt = open(lg).read() if os.path.exists(lg) else ""
        mm = re.search(r"index drop 'ne_lumen_b6': ([\d,]+)", txt)
        n = int(mm.group(1).replace(",", "")) if mm else None
        exp = (pre.get("C3") or {}).get("expected_index_drop_instances")
        items.append(judged("B6.B4 index drop", None if n is None else bool(exp is not None and n <= exp), dict(dropped=n, pre_build_expected=exp),
                            "pack index log `index drop 'ne_lumen_b6'` <= the pre-build expected count (<= 1 % per module)",
                            "logs/pack_index.log", "" if n is not None else "no `index drop 'ne_lumen_b6'` line in logs/pack_index.log"))
        c1 = (pre.get("C1") or {}).get("pass_")
        items.append(judged("B6.C1 chromatin inside INM", c1, (pre.get("C1") or {}).get("share", c1), ">= 99.9 % (pre-build; rep5 frozen, INM unchanged)",
                            "b6 geometry_check.json C1"))
    else:
        ref_v = (pre.get("measures_reference") or {}).get("volume_um3")
        items.append(judged("B6.B3 volumes identical", None if not (meas and ref_v) else meas["volume_um3"] == ref_v,
                            None, "MC volumes identical to the old package (Y5)", "geometry/measures.json vs pre-build measures_reference"))
    S = ctx.get("suite_new") or {}
    c16 = sval(S, "C16", "depends_on_geometry_mismatch")
    items.append(judged("B6.C4/C16 depends_on_geometry", None if c16 is None else c16 == 0, c16, "C16 depends_on_geometry_mismatch = 0", "suite.json"))
    if cell == "jurkat":
        c13 = sval(S, "C13", "prepass_md5_ok")
        items.append(judged("B6.B4 C13 prepass_md5_ok", None if c13 is None else c13 == 1, c13, "C13 prepass_md5_ok = 1 (frozen chromatin)", "suite.json"))
    return dict(items=items, superseded="B6 B2 (suite no_regression) -> combined no_regression")


# ------------------------------------------------------------------------------------------------ B3 (Jurkat)
def sec_b3(cell, bdir, ctx):
    import pandas as pd
    items = []
    meta = jload(os.path.join(bdir, "library", "library_meta.json")) or {}
    steps = meta.get("steps", [])
    items.append(judged("B3.B1 library mode", None if not meta else (meta.get("mode") == "build" and "tm_route" in steps and "adaptive_trees" in steps),
                        dict(mode=meta.get("mode"), steps=steps), "library stage mode build, tm_route + adaptive_trees steps",
                        "library/library_meta.json"))
    mp, rp = os.path.join(bdir, "library", "library_manifest.tsv"), os.path.join(REF[cell], "library", "library_manifest.tsv")
    if os.path.exists(mp):
        M = pd.read_csv(mp, sep="\t", keep_default_na=False); Rm = pd.read_csv(rp, sep="\t", keep_default_na=False)
        items.append(judged("B3.B4 provenance", bool((M.origin != "unknown").all()) and bool((meta.get("results") or {}).get("inputs")),
                            int((M.origin == "unknown").sum()), "every library file has an origin; inputs recorded in library_meta.json",
                            "library/library_manifest.tsv"))
        a, b = dict(zip(Rm.species, Rm.md5)), dict(zip(M.species, M.md5))
        missing = sorted(set(a) - set(b)); diff = sorted(k for k in a if k in b and a[k] != b[k]); extra = sorted(set(b) - set(a))
        same = not missing and not diff
        refv = jload(os.path.join(REF[cell], "validate", "compare_B3", "verdict.json")) or {}
        lib_items = {k: v for k, v in (refv.get("items") or {}).items() if k.startswith("L")}
        if same:
            items.append(judged("B3.L1-L7", bool(lib_items) and all(lib_items.values()), dict(reference_items=lib_items, extra_species=extra[:50],
                                n_extra=len(extra)), "every species of the accepted B3 library present with the same md5 -> L1-L7 of the "
                                "accepted B3 verdict carry over (extra species listed)", f"{rp}; {REF[cell]}/validate/compare_B3/verdict.json"))
        elif ctx["quick"]:
            items.append(item("B3.L1-L7", "NA", dict(missing=len(missing), differing=len(diff)), "L1-L7 vs v7 library",
                              "library differs from the accepted B3 library and --quick skips `library_b3_check.py compare`"))
        else:
            import subprocess
            o = os.path.join(ctx["out"], "b3_compare")
            subprocess.run([sys.executable, os.path.join(REPO, "scripts", "dev", "library_b3_check.py"), "compare", "--build", bdir, "--out", o,
                            "--nproc", str(ctx["nproc"])], cwd=REPO, check=False)
            L = jload(os.path.join(o, "library_compare.json")) or {}
            got = {k: (L.get(k) or {}).get("pass") for k in ("L1", "L2", "L3", "L4", "L5a", "L6", "L7")}
            items.append(judged("B3.L1-L7", None if not L else all(bool(v) for v in got.values()), dict(items=got, missing=len(missing), differing=len(diff)),
                                "L1-L7 (L5a addendum) of the build's library vs v7", f"{o}/library_compare.json"))
    else:
        items.append(item("B3.L1-L7", "NA", None, "library vs accepted B3 library", "library/library_manifest.tsv missing"))
    S, R = ctx.get("suite_new") or {}, ctx.get("suite_ref") or {}
    from wholecell.validate.suite import TOL
    ok = []
    for m in ("sign_errors_per_1e6", "median_tilt_deg"):
        v, r = sval(S, "C9", m), sval(R, "C9", m)
        rr = R.get(("C9", m))
        tol = TOL[rr["tol_kind"]][0](r) if rr and r is not None else None
        ok.append(dict(metric=m, value=v, reference=r, tolerance=tol, ok=None if v is None or r is None or tol is None else v <= r + tol))
    items.append(judged("B3.B3 C9", None if any(x["ok"] is None for x in ok) else all(x["ok"] for x in ok), ok,
                        "C9 sign errors and median tilt <= reference + §4 tolerance", "suite.json"))
    return dict(items=items, superseded="B3 B2 (suite no_regression) -> combined no_regression")


# ------------------------------------------------------------------------------------------------ B4
def _find_log(bdir, pat):
    for f in sorted(glob.glob(os.path.join(bdir, "logs", "*.log")) + glob.glob(os.path.join(bdir, "preplaced", "**", "*.log"), recursive=True)):
        try:
            if re.search(pat, open(f, errors="replace").read()):
                return f
        except OSError:
            pass
    return None


def sec_b4(cell, bdir, ctx):
    import b4_preplaced_check as B4
    items = []
    serial = os.path.join(bdir, "preplaced", "serial")
    complete = os.path.exists(os.path.join(serial, "serial_manifest.json"))
    items.append(judged("B4.B1 serial stage", bool(ctx["complete"]) and complete, complete,
                        "build complete; preplaced/serial/serial_manifest.json lists the B4 builders", "preplaced/serial"))
    S = ctx.get("suite_new") or {}
    if cell == "jurkat":
        if not complete:
            items.append(item("B4.J1-J7", "NA", None, "B4 J items", "preplaced/serial missing"))
            return dict(items=items)
        out = os.path.join(ctx["out"], "b4"); os.makedirs(out, exist_ok=True)
        pkp = os.path.join(out, "packed_ribo80s.npy")              # cached for re-runs of this check (reads every tile)
        if not os.path.exists(pkp):
            np.save(pkp, B4._packed_ribosomes(bdir))
        a = argparse.Namespace(cell="jurkat", serial=serial, build=bdir, ref_build=V7REPRO2, out=os.path.join(out, "measure.json"),
                               packed=os.path.join(out, "packed_ribo80s.npy"), no_overlaps=bool(ctx["quick"]),
                               lamina=B4._selected_npz(bdir, "lamina"))
        mj = os.path.join(out, "measure.json")
        m = jload(mj) if ctx.get("reuse") and os.path.exists(mj) else None      # --reuse-from: same build, same measure code
        if not m or (not ctx["quick"] and "J5" not in m):
            m = B4.measure(a)
        items.append(judged("B4.J1 pores", bool(m["J1"]["pass"] and m["J1"]["inventory_ok"]), m["J1"],
                            "pores = config = package = placed (2,696), anchor ratio 0.5-2, min separation >= 150 nm, inventory NPC = config", "measure.json"))
        items.append(judged("B4.J2 MT dimers", m["J2"]["pass_"], m["J2"]["dimers"], "2,325,171 +- 1 %", "measure.json"))
        items.append(judged("B4.J3a polysome fraction", m["J3"]["pass_a"], m["J3"]["a"], "0.68-0.72", "measure.json"))
        lim = "accepted model limitation: B4_preplaced_addendum_polysomes_2026-11-18.md §4 (main session decision 2026-11-18)"
        b = m["J3"].get("b_with_packed")
        items.append(item("B4.J3b linked share", "PASS" if m["J3"]["pass_b"] else ("ACCEPTED_DEVIATION" if b is not None else "NA"), b,
                          "|linked - 0.697| <= 0.05", "" if m["J3"]["pass_b"] else (lim if b is not None else "no packed ribosomes"), "measure.json"))
        nn = ((m["J4"].get("with_packed") or {}).get("nn_median_thinned_to_EM"))
        items.append(item("B4.J4 ribosome NN", "PASS" if m["J4"].get("pass") else ("ACCEPTED_DEVIATION" if nn is not None else "NA"), nn,
                          "27.24-29.24 nm (thinned to EM)", "" if m["J4"].get("pass") else (lim if nn is not None else "no packed ribosomes"), "measure.json"))
        pl = jload(os.path.join(bdir, "validate", "placement.json"))
        if pl:
            rate = pl["clash_pairs"]["prepass_prepass"] / max(pl.get("n_prepass", 1), 1) * 1e6
            ref = B4.V7["prepass_prepass_per_1e6"]
            items.append(judged("B4.J5a prepass-prepass", rate <= ref + max(0.1 * ref, 1.0), rate, f"<= {ref:.1f} + 10 %", "validate/placement.json"))
        else:
            items.append(item("B4.J5a prepass-prepass", "NA", None, "<= v7 + 10 %", "validate/placement.json missing"))
        if "J5" in m:
            ad = m["J5"].get("addendum") or {}
            items.append(judged("B4.J5b per class (addendum 2026-11-17)", ad.get("pass_b"), ad.get("failing"),
                                "every B4 pair class <= reference + max(10 %, 1) per 1e6 B4 instances (absent class = new)", "measure.json"))
        else:
            items.append(item("B4.J5b per class (addendum 2026-11-17)", "NA", None, "per class rule", "--quick: overlaps skipped"))
        j6 = m["J6"]; mon = j6.get("atp_monomers") or 0; pool = j6.get("pool_cpx6151") or 0
        md = j6.get("mtdna") or {}
        items.append(judged("B4.J6 counts", bool(mon > 0 and pool and mon / 2 <= pool / 2 and md), dict(monomers=mon, pool=pool, mtdna=bool(md)),
                            "dimers > 0 and <= pool / 2; mtDNA layer present", "measure.json"))
        lgp = _find_log(bdir, r"rim sites tried")
        jo = os.path.join(out, "j6_atp.json")
        import subprocess
        cmd = [sys.executable, os.path.join(REPO, "scripts", "dev", "b4_j6_atp_check.py"), "--serial", serial, "--build", bdir, "--out", jo,
               "--sample", "1000" if ctx["quick"] else "4000"] + (["--log", lgp] if lgp else []) + (["--skip-overlaps"] if ctx["quick"] else [])
        if not os.path.exists(jo):                                  # cached for re-runs of this check
            subprocess.run(cmd, cwd=REPO, check=False)
        J = jload(jo)
        if J:
            mb = J.get("membrane") or {}
            # addendum J6 §3: "J6 geometry tests at the addendum tolerance: 100 % (up to float32 export rounding)" -> the
            # exported float32 poses are judged with F32_EPS_NM (b4_j6_atp_check.py); strict values kept in `value`
            pdmd = mb.get("per_dimer_max_tm_nm") or {}
            pdm, pdm32 = pdmd.get("share_le", {}).get("4.5"), pdmd.get("share_le_4_5_f32")
            sold = mb.get("soluble_in_membrane_share") or {}
            sol, sol32 = sold.get("dimers_over_2pct"), sold.get("dimers_over_2pct_f32")
            f32 = "float32 export rounding (eps 1e-3 nm) allowed, addendum J6 §3"
            items.append(judged("B4.J6 TM tolerance (addendum J6: 4.5 nm)", None if pdm32 is None else pdm32 >= 1.0 - 1e-9,
                                dict(share_le_4_5=pdm, share_le_4_5_f32=pdm32, n_over_4_5=pdmd.get("n_over_4_5"), max_nm=pdmd.get("max_exact")),
                                f"every sampled dimer's TM beads within 4.5 nm of the IM midsurface; {f32}", jo,
                                "" if pdm32 is not None else "j6_atp.json predates the float32-tolerance fields: re-run b4_j6_atp_check.py"))
            items.append(judged("B4.J6 soluble in membrane", None if sol32 is None else sol32 == 0, dict(dimers_over_2pct=sol, dimers_over_2pct_f32=sol32,
                                max_share=sold.get("max")), f"<= 2 % of tested soluble beads in a bilayer, per dimer; {f32}", jo,
                                "" if sol32 is not None else "j6_atp.json predates the float32-tolerance fields: re-run b4_j6_atp_check.py"))
            dens = J.get("density") or {}
            items.append(judged("B4.J6 rim density (addendum J6)", bool(dens.get("in_literature_band") or dens.get("in_literature_traced")),
                                dens, "dimers per um of rim (traced or band) in 50-83 (literature)", jo))
            ov = J.get("overlap_atp_atp")
            items.append(item("B4.J6 overlaps", "REPORT", dict(atp_atp=ov, atp_other=J.get("overlap_atp_other")) if ov else None,
                              "reported (addendum J6)", "" if ov else "--quick: overlaps skipped", jo))
        else:
            items.append(item("B4.J6 TM tolerance (addendum J6: 4.5 nm)", "NA", None, "4.5 nm", f"b4_j6_atp_check.py produced no {jo}"))
        items.append(item("B4.J7 consumption", "REPORT", m.get("J7"), "module json consumed = placed", "not automated beyond the report", "measure.json"))
        items.append(judged("B4.J8 depends_on_geometry", None if sval(S, "C16", "depends_on_geometry_mismatch") is None else sval(S, "C16", "depends_on_geometry_mismatch") == 0,
                            sval(S, "C16", "depends_on_geometry_mismatch"), "0 mismatches", "suite.json C16"))
    else:
        if not complete:
            items.append(item("B4.Y1-Y3", "NA", None, "", "preplaced/serial missing"))
        else:
            ref = np.load(os.path.join(B4.YREF_DIR, "prepass.npz")); new = np.load(os.path.join(serial, "preplaced.npz"))
            # library-index arrays (li_*, *_type) compared by library NAME: the combined build's library order differs
            # from s0 whenever the species list does (B5: 88 heterodimer / oligomer species differ), so raw indices
            # differ even where the same structures are placed (fix 2026-11-19; the reference = the s0 library order)
            rl = np.array([str(x) for x in np.load(YREF_LIB_IDS)]) if os.path.exists(YREF_LIB_IDS) else None
            nl = np.array([str(x) for x in np.load(os.path.join(serial, "lib_ids.npy"))])

            def _same(k):
                if k not in new.files:
                    return False
                a, b = ref[k], new[k]
                if rl is not None and (k.startswith("li_") or k.endswith("_type")) and a.dtype.kind in "iu" and b.dtype.kind in "iu":
                    if a.size and int(a.max()) >= len(rl) or b.size and int(b.max()) >= len(nl):
                        return False
                    a, b = rl[a], nl[b]
                return a.shape == b.shape and np.array_equal(a, b)
            diff = [k for k in ref.files if k not in ("mrna_gene", "poly_chain_gene") and not _same(k)]
            man = jload(os.path.join(serial, "serial_manifest.json")) or {}
            owner = {}
            for bld in man.get("builders", []):
                for k in bld.get("keys", []):
                    owner.setdefault(k, bld["name"])
            cfgr = jload(os.path.join(bdir, "config.resolved.json")) or {}
            b5_on = (((cfgr.get("inventory") or {}).get("params") or {}).get("localisation") or {}).get("preset") == "common"
            loc_dep = [k for k in diff if owner.get(k) in Y1_LOCALISATION_BUILDERS] if b5_on else []
            hard = [k for k in diff if k not in loc_dep]
            items.append(judged("B4.Y1 arrays", not hard, dict(differing=hard, by_name=rl is not None),
                                "every pre-pass array identical to the PATCH-4 reference by library name (PATCH-5 gene labels "
                                "excepted); with B5 on, the arrays of the localisation-fed builders are reported below", B4.YREF_DIR))
            if b5_on:
                items.append(item("B4.Y1 localisation-fed builders (B5)", "REPORT",
                                  {b: [k for k in loc_dep if owner.get(k) == b] for b in sorted({owner[k] for k in loc_dep})},
                                  "reported: builders whose inputs are the inventory's localisation shares / pools (B5 changes them)",
                                  "EXPECTED from B5 (accepted; yeast curated / builder-anchored overrides, "
                                  "B5_localisation_record_yeast_overrides_2026-11-19.md): " + ", ".join(Y1_LOCALISATION_BUILDERS),
                                  "preplaced/serial/serial_manifest.json"))
            items.append(judged("B4.Y2 pores", int(len(new["npc_pos"])) == B4.YREF["pores"], int(len(new["npc_pos"])), "93", "preplaced.npz"))
            items.append(judged("B4.Y3 MTs", bool(np.array_equal(ref["mt_pos"], new["mt_pos"])), None, "mt_pos identical", "preplaced.npz"))
        ic = sval(S, "C12", "polysome_in_chains_fraction")
        items.append(item("B4.Y4 polysome fraction", "REPORT", ic, "reported vs 0.82-0.85 (reference 0.784)", "", "suite.json"))
        nn = sval(S, "C12", "nn_median_model")
        items.append(judged("B4.Y5 ribosome NN", None if nn is None else abs(nn - B4.YREF["nn"]) <= 0.5, nn, f"{B4.YREF['nn']:.2f} +- 0.5 nm", "suite.json"))
        pp = sval(S, "C7", "a.prepass_prepass_per_1e6")
        items.append(judged("B4.Y6 prepass-prepass", None if pp is None else pp <= B4.YREF["prepass_prepass_per_1e6"] * 1.1 + 1, pp,
                            "<= 642,511 + 10 %", "suite.json"))
    return dict(items=items, superseded="B4 B2 (suite no_regression) -> combined no_regression")


# ------------------------------------------------------------------------------------------------ B1
def _c10_whole_cell(cell, bdir, ctx, tag):
    """A4 placed_of_unthinned over all tiles + per-region placed mass fraction (core layout: pack/tiles, counts.npz)"""
    from wholecell.validate import common
    cfg = wcfg.load(cell)
    common.CACHE_DIR = os.path.join(ctx["out"], f"cache_{tag}")
    os.makedirs(common.CACHE_DIR, exist_ok=True)
    try:
        inp = common.Inputs(cfg, bdir)
        own = os.path.join(bdir, "validate", "cache", "tile_sums.npz")
        T = common.tile_sums(inp, own if os.path.exists(own) else os.path.join(common.CACHE_DIR, "tile_sums.npz"), nproc=ctx["nproc"])
        L = json.load(open(inp["libs"]))
        m, msrc = common.lib_masses(inp, bdir, nproc=ctx["nproc"])
    finally:
        common.CACHE_DIR = None
    rn = L.get("rnames") or []
    nl = len(m)
    req = (T["req"][:nl] * m[:, None]).sum(0); pla = (T["pla"][:nl] * m[:, None]).sum(0)
    cz = os.path.join(bdir, "pack", "counts.npz")
    unth = None
    if os.path.exists(cz):
        c = np.load(cz)
        unth = np.zeros(T["req"].shape[1]); np.add.at(unth, c["region"].astype(np.int64), c["n"] * m[c["lib"].astype(np.int64)])
    reg = {rn[i]: dict(placed_mass_frac=float(pla[i] / req[i]) if req[i] else None, placed_pg=float(pla[i] / AV * 1e12),
                       placed_of_unthinned=float(pla[i] / unth[i]) if unth is not None and unth[i] else None) for i in range(len(rn))}
    return dict(placed_mass_fraction=float(pla.sum() / req.sum()), placed_of_unthinned=float(pla.sum() / unth.sum()) if unth is not None else None,
                regions=reg, mass_source=msrc)


GATED_K2 = {"jurkat": ("cytosol", "nucleoplasm", "mito_matrix", "er_lumen", "nucleolus", "mito_ims", "golgi_lumen",
                       "endolysosome_lumen", "perinuclear_space"),
            "yeast": None}


def _donor(bdir):
    """donor build of a downstream-reuse build (scripts/build_downstream_reuse.py: upstream rules ran there), else None"""
    return (jload(os.path.join(bdir, "DOWNSTREAM_REUSE.json")) or {}).get("donor")


def _rule_time(bdir, rule):
    import b2_organisation_check as B2
    t = B2._rule_time(bdir, rule)
    if t is None and _donor(bdir):                 # upstream rule of a downstream-reuse build: its donor's snakemake.log
        t = B2._rule_time(_donor(bdir), rule)
    return t


def _host(bdir):
    t = open(os.path.join(bdir, "launch.log")).read() if os.path.exists(os.path.join(bdir, "launch.log")) else ""
    m = re.search(r"^host:\s*(\S+)", t, re.M)
    if not m and _donor(bdir):
        return _host(_donor(bdir))
    return m.group(1) if m else None


def _densities_json(bdir):
    """region volumes: validate/densities.json, else validate/core/densities.json (yeast: the suite runs the core
    densities step into validate/core; fix 2026-11-19, K2 / K6 found no volumes for yeast)"""
    for f in (os.path.join(bdir, "validate", "densities.json"), os.path.join(bdir, "validate", "core", "densities.json")):
        if os.path.exists(f):
            return jload(f) or {}
    return {}


def _pack_code(bdir):
    """code checkout that ran the pack of bdir: DOWNSTREAM_REUSE.json donor's _code (scripts/build_downstream_reuse.py),
    else <bdir>/_code, else this repository"""
    dr = jload(os.path.join(bdir, "DOWNSTREAM_REUSE.json")) or {}
    for c in ([os.path.join(dr["donor"], "_code")] if dr.get("donor") else []) + [os.path.join(bdir, "_code"), REPO]:
        if os.path.isdir(os.path.join(c, "wholecell")):
            return c
    return REPO


def sec_b1(cell, bdir, ctx):
    items = []
    S, R = ctx.get("suite_new") or {}, ctx.get("suite_ref") or {}
    # K1 / K2
    try:
        new = (ctx.get("c10") or {}).get("build") or _c10_whole_cell(cell, bdir, ctx, "build")
        json.dump(new, open(os.path.join(ctx["out"], "c10_build.json"), "w"), indent=1)
    except Exception as e:
        new = None; ctx.setdefault("notes", []).append(f"B1 K1/K2 (build): {type(e).__name__}: {e}")
    ref = None
    if cell == "jurkat":
        try:
            ref = (ctx.get("c10") or {}).get("reference") or _c10_whole_cell(cell, REF[cell], ctx, "reference")
        except Exception as e:
            ctx.setdefault("notes", []).append(f"B1 K1 (reference): {type(e).__name__}: {e}")
    if new and ref and new["placed_of_unthinned"] is not None and ref["placed_of_unthinned"] is not None:
        items.append(judged("B1.K1 placed of unthinned", new["placed_of_unthinned"] >= ref["placed_of_unthinned"] - 0.003,
                            dict(build=new["placed_of_unthinned"], reference=ref["placed_of_unthinned"]), ">= reference - 0.003", "pack/tiles, counts.npz"))
    else:
        v, r = sval(S, "C10", "placed_mass_fraction"), sval(R, "C10", "placed_mass_fraction")
        items.append(item("B1.K1 placed of unthinned", "NA", dict(build_unthinned=(new or {}).get("placed_of_unthinned"), c10_build=v, c10_reference=r),
                          ">= reference - 0.003",
                          "reference placed_of_unthinned not computable (yeast reference = port layout); suite C10 placed_mass_fraction shown "
                          "for information (same thinning preset legacy_x5 in both)" if cell == "yeast" else "tile sums not computable"))
    if new:
        vols = _densities_json(bdir).get("rows") or []
        big = {r["region"] for r in vols if (r.get("volume_um3") or 0) >= 1.0}
        reg = {k: v for k, v in new["regions"].items() if k in big and v["placed_mass_frac"] is not None}
        low = {k: v["placed_mass_frac"] for k, v in reg.items() if v["placed_mass_frac"] < 0.98}
        th = jload(os.path.join(bdir, "pack", "thinning.json"))
        refreg = (ref or {}).get("regions") or {}
        worse = {k: dict(build=x, reference=(refreg.get(k) or {}).get("placed_mass_frac")) for k, x in low.items()
                 if (refreg.get(k) or {}).get("placed_mass_frac") is None or x < refreg[k]["placed_mass_frac"] - 0.003}
        pre = {k: dict(build=x, reference=refreg[k]["placed_mass_frac"]) for k, x in low.items() if k not in worse}
        items.append(judged("B1.K2 regions >= 1 um3", None if not reg else not worse, dict(below_0_98_and_below_reference=worse, regions=len(reg),
                            thinning=bool(th)), "placed mass fraction >= 0.98 per volume region >= 1 um3; a region below 0.98 that is "
                            "not below the reference - 0.003 is listed in `B1.K2 causes` (pre-existing), not failed here",
                            "pack/tiles; validate/densities.json", "" if reg else "no region volumes (densities.json) or tile sums"))
        if pre:
            items.append(item("B1.K2 causes", "NA", pre, "regions < 0.98 need a listed cause (geometry limit / thinning)",
                              "below 0.98 also in the reference build: the cause list is a review item, not automated"))
    else:
        items.append(item("B1.K2 regions >= 1 um3", "NA", None, ">= 0.98", "tile sums not computable"))
    si = sval(S, "C10", "seam_index")
    items.append(judged("B1.K3 seam", None if si is None else abs(si - 1) <= 0.1, si, "|seam_index - 1| <= 0.1", "suite.json C10"))
    pp = sval(S, "C7", "a.packed_packed_per_1e6")
    k4 = [None if pp is None else bool(pp <= 1)]
    v4 = dict(packed_packed=pp)
    if cell == "yeast":
        b0, b1, rb = sval(S, "C7", "b.bead_clash_rate_without_actin_structures"), sval(S, "C7", "b.bead_clash_rate"), sval(R, "C7", "b.bead_clash_rate")
        v4.update(bead_wo_actin=b0, bead=b1, bead_ref=rb)
        k4 += [None if b0 is None else bool(b0 <= 1e-5), None if b1 is None or rb is None else bool(b1 <= rb + max(0.1 * rb, 1e-6))]
    # a missing component (yeast C7(b) comes from the yeast validation scripts) makes K4 NA, not FAIL (fix 2026-11-19)
    k4_ok = False if any(x is False for x in k4) else (None if pp is None or any(x is None for x in k4) else True)
    items.append(judged("B1.K4 clashes", k4_ok, v4,
                        "C7(a) packed-packed <= 1 per 1e6; yeast C7(b) without actin <= 1e-5 and bead rate <= ref + max(10 %, 1e-6)", "suite.json C7",
                        "" if k4_ok is not None else "value missing"))
    k5 = []
    for mtr in ("a.containment_rate_per_1e6",) + (("b.beads_outside_region_fraction",) if cell == "yeast" else ()):
        v, r = sval(S, "C8", mtr), sval(R, "C8", mtr)
        tol = None if r is None else (max(0.1 * r, 1.0) if mtr.startswith("a.") else max(0.1 * r, 1e-6))
        k5.append(dict(metric=mtr, value=v, reference=r, ok=None if v is None or r is None else v <= r + tol))
    items.append(judged("B1.K5 containment", None if any(x["ok"] is None for x in k5) else all(x["ok"] for x in k5), k5,
                        "C8 <= reference + max(10 %, 1 per 1e6)", "suite.json C8"))
    dens = []
    for (c, mtr), r in R.items():
        if c != "C11" or not mtr.startswith("density."):
            continue
        v, rv = sval(S, c, mtr), _f(r.get("value"))
        if v is None or rv is None:
            continue
        dens.append(dict(metric=mtr, value=v, reference=rv, delta=round(v - rv, 3), within_3=abs(v - rv) <= 3.0))
    dpl = _placed_density_change(cell, ctx)
    for d in dens:
        x = dpl.get(d["metric"].split(".", 1)[1]) if dpl else None
        d["d_rho_placed"] = None if x is None else round(x, 3)
        d["mass_consistent"] = None if x is None else abs(d["delta"] - x) <= 2.0
    bad = [d for d in dens if not d["within_3"]]
    items.append(judged("B1.K6 densities", None if not dens or any(d["mass_consistent"] is None for d in bad) else
                        all(d["mass_consistent"] for d in bad), dict(outside_3_mg_ml=bad, n=len(dens)),
                        "|d rho| <= 3 mg/ml; or, when the region's placed mass changed, mass-consistent |d rho - d m_placed / V| <= 2 mg/ml "
                        "(B1 K6; d m_placed = packed mass placed per region, all tiles); a placed-mass deficit is judged by K2",
                        "suite.json C11; pack/tiles", "" if dens else "no C11 densities"))
    rref = "pack_pack" if cell == "jurkat" else "engine_pack"      # yeast port reference: one engine_pack rule
    tw, tr = _rule_time(bdir, "pack_pack"), _rule_time(REF[cell], rref)
    hn, hr = _host(bdir), _host(REF[cell])
    if tw is None or tr is None:
        items.append(item("B1.K7 wall time", "NA", dict(build_s=tw, reference_s=tr), "pack_pack <= 2 x reference",
                          f"rule time not in snakemake.log (build pack_pack: {tw}, reference {rref}: {tr})"))
    elif hn is None or hn != hr:
        items.append(item("B1.K7 wall time", "REPORT", dict(build_s=tw, reference_s=tr, ratio=round(tw / tr, 3), host=hn, ref_host=hr),
                          "pack_pack <= 2 x reference (same node only)", "different or unknown node: reported, not gated (B1 Part B)"))
    else:
        items.append(judged("B1.K7 wall time", tw <= 2 * tr, dict(build_s=tw, reference_s=tr, ratio=round(tw / tr, 3), host=hn, reference_rule=rref),
                            "pack_pack <= 2 x reference (same node)", "snakemake.log",
                            "" if cell == "jurkat" else "yeast reference rule engine_pack also contains the port pre-pass (lenient comparison)"))
    pdj = jload(os.path.join(bdir, "pack", "pack_done.json")) or {}
    try:
        cfgr = jload(os.path.join(bdir, "config.resolved.json")) or {}
        # the switches of the code that PACKED this build (its _code worktree; a downstream-reuse build: the donor's), not
        # of this script's checkout: later pack defaults (e.g. a new organisation key) are not settings of this build
        code = _pack_code(bdir)
        cmd = ("import json,sys; from wholecell.pack import switches as S; e=S.effective(json.load(sys.stdin)); "
               "print(json.dumps(dict(eff=e, want=S.load(e).as_dict(), cd=getattr(S,'COMMON_DEFAULTS',None), "
               "defaults=S.DEFAULTS), default=str))")
        import subprocess
        r = subprocess.run([sys.executable, "-c", cmd], input=json.dumps(cfgr.get("pack") or {}), capture_output=True, text=True,
                           cwd=code, env=dict(os.environ, PYTHONPATH=code, PYTHONHASHSEED="0"), check=True)
        SW = json.loads(r.stdout.strip().splitlines()[-1])
        eff, cd = SW["eff"], SW["cd"]
        switches = argparse.Namespace(DEFAULTS=SW["defaults"])
        if "switches" in pdj:
            # pack_done.json records Switches.as_dict() (the algorithm keys of switches.DEFAULTS), not the whole `pack:`
            # section: compare like with like (the 2026-11-19 verdict compared the record with effective() -> always FAIL)
            want = SW["want"]
            same = json.dumps(pdj["switches"], sort_keys=True, default=str) == json.dumps(want, sort_keys=True, default=str)
            diff = sorted(k for k in set(want) | set(pdj["switches"]) if json.dumps(want.get(k), sort_keys=True, default=str)
                          != json.dumps(pdj["switches"].get(k), sort_keys=True, default=str))
            v7 = None
            if cell == "jurkat":        # B1 outcome (addendum2 2026-11-18): COMMON_DEFAULTS = {} -> every switch v7 except B2 organisation
                v7 = sorted(k for k in switches.DEFAULTS if k != "organisation" and json.dumps(pdj["switches"].get(k), sort_keys=True, default=str)
                            != json.dumps(switches.DEFAULTS[k], sort_keys=True, default=str))
            items.append(judged("B1.K8 settings applied", same and not v7 and cd == {} if cell == "jurkat" else same,
                                dict(common_defaults=cd, preset=eff.get("_preset"), differing_keys=diff, non_v7_keys_jurkat=v7,
                                     organisation_enabled=(pdj["switches"].get("organisation") or {}).get("enabled")),
                                "pack_done.json switches = switches.load(effective(pack)) (COMMON_DEFAULTS merged on the cell); Jurkat: "
                                "COMMON_DEFAULTS = {} and every switch at its v7 value except organisation (B2 on) (B1 addendum2 2026-11-18)",
                                f"pack/pack_done.json, config.resolved.json; switches code {code}"))
        else:
            items.append(item("B1.K8 settings applied", "NA", dict(common_defaults=cd, preset=(cfgr.get("pack") or {}).get("preset")),
                              "pack_done.json switches = effective settings", "pack/pack_done.json has no `switches` record"))
    except Exception as e:
        items.append(item("B1.K8 settings applied", "NA", None, "", f"{type(e).__name__}: {e}"))
    return dict(items=items, superseded="B1 K9 (suite no_regression) -> combined no_regression",
                decision="COMMON_DEFAULTS = {} (v7 values); yeast keeps preset legacy_x5 (C3 yeast pin not applied)")


# ------------------------------------------------------------------------------------------------ B2
def sec_b2(cell, bdir, ctx):
    import b2_organisation_check as B2
    items = []
    J = jload(os.path.join(bdir, "pack", "organisation", "network.json")) or {}
    st, sg = J.get("stats", {}), J.get("string", {})
    taxon = 9606 if cell == "jurkat" else 4932
    md5 = ((J.get("inputs") or {}).get("string_links") or {}).get("md5")
    chk = dict(network_json=bool(J), pool_pairs_0=st.get("pool_pairs") == 0, contact_links_0=st.get("contact_links") == 0,
               taxon=sg.get("taxon") in (taxon, str(taxon)), version=str(sg.get("version")) == "12.0",
               min_score=int(sg.get("min_score") or 0) == 400, links_md5=(md5 == B2.STRING_MD5[taxon]) if taxon in B2.STRING_MD5 else bool(md5))
    items.append(judged("B2.O0 settings", all(chk.values()), chk, "network.json: pool_pairs 0, contact_links 0, STRING taxon / 12.0 / >= 400, "
                        "links md5 (config.resolved diff vs a B1 build: n/a, no B1 build)", "pack/organisation/network.json"))
    items.append(item("B2.O1 compartment totals", "NA", None, "vs the B1 build of the same cell",
                      "no B1 build exists and the reference build differs by B4 / B5 (requests change): totals are confounded"))
    if cell == "yeast":
        S, R = ctx.get("suite_new") or {}, ctx.get("suite_ref") or {}
        rows = []
        for (c, mtr), r in R.items():
            if c == "C15" and mtr.endswith("median_auroc.snapshot_assembly"):
                v, rv = sval(S, c, mtr), _f(r.get("value"))
                rows.append(dict(metric=mtr, value=v, reference=rv, ok=None if v is None or rv is None else v >= rv - 0.01))
        items.append(judged("B2.O2 C15", None if not rows or any(x["ok"] is None for x in rows) else all(x["ok"] for x in rows), rows,
                            "each C15 median_auroc.snapshot_assembly >= reference - 0.01", "suite.json C15",
                            "" if rows else "no C15 rows in the reference suite"))
    pf = B2.partner_fraction(bdir)
    if cell == "jurkat":
        f = pf.get("fraction_next_to_partner")
        items.append(judged("B2.O3 partner placement", None if f is None else f > 0, pf, "fraction placed next to a STRING partner > 0", "pack/summary.json"))
    else:
        items.append(item("B2.O3 partner placement", "REPORT", pf, "reported (pool off: contact partners gone)", "", "pack/summary.json"))
    S = ctx.get("suite_new") or {}
    c16 = [sval(S, "C16", "manifest_stale_files"), sval(S, "C16", "depends_on_geometry_mismatch")]
    items.append(judged("B2.O5 complete", bool(ctx["complete"]) and all(x in (0, None) for x in c16) and any(x is not None for x in c16),
                        dict(complete=ctx["complete"], manifest_stale_files=c16[0], depends_on_geometry_mismatch=c16[1]),
                        "build complete; C16 0 stale, 0 depends_on_geometry mismatches", "suite.json C16"))
    try:
        items.append(item("B2.wall time", "REPORT", {r: _rule_time(bdir, r) for r in ("pack_organisation", "pack_counts", "pack_pack")},
                          "reported, not gated", "", "snakemake.log"))
    except Exception:
        pass
    return dict(items=items, superseded="B2 O4 (suite no_regression) -> combined no_regression")


# ------------------------------------------------------------------------------------------------ B5
def sec_b5(cell, bdir, ctx):
    import b5_localisation_check as B5
    items = []
    S, R = ctx.get("suite_new") or {}, ctx.get("suite_ref") or {}
    cfgr = jload(os.path.join(bdir, "config.resolved.json")) or {}
    preset = (((cfgr.get("inventory") or {}).get("params") or {}).get("localisation") or {}).get("preset")
    items.append(judged("B5.switch", None if preset is None else preset == "common", preset, "inventory.params.localisation.preset = common",
                        "config.resolved.json"))
    sm = {k[1]: r for k, r in S.items()}; rm = {k[1]: r for k, r in R.items()}
    for lvl in ("intended", "placed"):
        for s in B5.C1_SETS[cell]:
            k = f"{s}.{lvl}.{B5.C1_HEAD[cell]}.unweighted"
            v, r = _f((sm.get(k) or {}).get("value")), _f((rm.get(k) or {}).get("value"))
            if s in GATED_B5[cell]:
                items.append(judged(f"B5.L1 {k}", None if v is None or r is None else v >= r - 0.005, dict(value=v, reference=r),
                                    ">= reference - 0.005 (gated: independent reference, addendum 2026-11-18)", "suite.json C1"))
            else:
                items.append(item(f"B5.L1 {k}", "REPORT", dict(value=v, reference=r), "reported, not gated (addendum 2026-11-18)", "", "suite.json C1"))
    if cell == "yeast":
        for mtr in ("image_nc.test.spearman", "image_nc.test.auroc_nuclear"):
            v, r = _f((sm.get(mtr) or {}).get("value")), _f((rm.get(mtr) or {}).get("value"))
            items.append(judged(f"B5.L3 {mtr}", None if v is None or r is None else v >= r - 0.01, dict(value=v, reference=r), ">= reference - 0.01", "suite.json C3"))
        for s in ("SWAT_native", "hyperLOPIT"):
            k = f"{s}.intended.tophit.unweighted"; v = _f((sm.get(k) or {}).get("value"))
            items.append(judged(f"B5.L4 {k}", None if v is None else v >= B5.YEAST_BENCH[s], v, f">= {B5.YEAST_BENCH[s]}", "suite.json C1"))
        v = _f((sm.get("image_nc.test.spearman") or {}).get("value"))
        items.append(judged("B5.L4 image_nc.test.spearman", None if v is None else v >= B5.YEAST_BENCH["image_nc"], v, ">= 0.80", "suite.json C3"))
    else:
        cfg = wcfg.load(cell)
        res = {}
        for tag, b in (("build", bdir), ("reference", REF[cell])):
            meas = jload(os.path.join(b, "geometry", "measures.json")); lm = os.path.join(b, "validate", "loc_model.parquet")
            res[tag] = B5.jurkat_c3i(cfg, lm, meas)[0] if meas and os.path.exists(lm) else {}
        json.dump(res, open(os.path.join(ctx["out"], "b5_c3i.json"), "w"), indent=1)
        for lvl in ("intended", "placed"):
            k = f"c3i.{lvl}.holdout.spearman"; v, r = res["build"].get(k), res["reference"].get(k)
            items.append(judged(f"B5.L3 {k}", None if v is None or r is None else v >= r - 0.01, dict(value=v, reference=r),
                                "C3i hold-out Spearman >= reference - 0.01", "validate/loc_model.parquet"))
    # L5 mass moved (build's own evidence: legacy preset -> common ablation)
    try:
        inv = os.path.join(bdir, "inventory")
        cfg_c, cfg_l = B5.cfg_for(cell, "common"), B5.cfg_for(cell, B5.LEGACY[cell])
        mass, total = B5.gene_mass(inv, cell)
        stages, genes = B5.run_ablation(cell, cfg_l, cfg_c, bdir, inv)
        rows = B5.mass_rows(stages, mass, total) + (B5.rule_rows(inv, "common") if cell == "jurkat" else [])
        B5.flag(rows, total)
        import pandas as pd
        pd.DataFrame(rows).to_csv(os.path.join(ctx["out"], "b5_mass_moved.tsv"), sep="\t", index=False)
        review = [r["rule"] for r in rows if r["review"]]
        items.append(item("B5.L5 mass moved", "REPORT", dict(max_pg=max(r["mass_moved_pg"] for r in rows), total_pg=total, review=review),
                          f"<= {B5.REVIEW_PG} pg and <= {B5.REVIEW_FRAC:.0%} of total protein per rule, else review (B5 accepted with the "
                          "revision's review list, Pedro 2026-11-18)", "", "b5_mass_moved.tsv"))
    except Exception as e:
        items.append(item("B5.L5 mass moved", "NA", None, "mass moved per rule", f"{type(e).__name__}: {e}"))
    return dict(items=items, superseded="B5 L2 / L6 (suite no_regression) -> combined no_regression; L7 = completeness")


# ------------------------------------------------------------------------------------------------ main
SECTIONS = [("B6", sec_b6), ("B3", sec_b3), ("B4", sec_b4), ("B1", sec_b1), ("B2", sec_b2), ("B5", sec_b5)]


def step_status(items):
    st = [i["status"] for i in items]
    if "FAIL" in st:
        return "FAIL"
    if "PASS" in st or "ACCEPTED_DEVIATION" in st:
        return "PASS"
    return "NA"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True, choices=["jurkat", "yeast"])
    ap.add_argument("--build", required=True)
    ap.add_argument("--complete", type=int, default=1)
    ap.add_argument("--out")
    ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--quick", action="store_true", help="skip the slow overlap measures (B4 J5b, J6 overlaps) and the B3 library compare")
    ap.add_argument("--reuse-from", help="earlier out dir of this check on the SAME build: reuse its caches (tile sums, B4 packed ribosomes, "
                    "B4 measure.json) instead of recomputing; the verdict is always recomputed")
    ap.add_argument("--accepted", action="append", default=[], help="dated accepted-deviations json (repeatable): matching "
                    "FAIL items / unexpected regressions become ACCEPTED_DEVIATION (see module docstring)")
    a = ap.parse_args()
    bdir = os.path.abspath(a.build.rstrip("/"))
    if os.path.realpath(bdir) == os.path.realpath(REF[a.cell]) and not a.out:
        raise SystemExit("refusing to write into the reference build: give --out")
    out = a.out or os.path.join(bdir, "validate", "compare_step2")
    os.makedirs(out, exist_ok=True)
    if a.reuse_from:
        import shutil
        src = os.path.abspath(a.reuse_from)
        assert os.path.realpath(src) != os.path.realpath(out), "--reuse-from must differ from --out"
        for d in ("cache_build", "cache_reference", os.path.join("b5_density", "cache_build"), os.path.join("b5_density", "cache_reference")):
            if os.path.isdir(os.path.join(src, d)) and not os.path.exists(os.path.join(out, d)):
                shutil.copytree(os.path.join(src, d), os.path.join(out, d))
        for f in (os.path.join("b4", "packed_ribo80s.npy"), os.path.join("b4", "measure.json")):
            if os.path.exists(os.path.join(src, f)) and not os.path.exists(os.path.join(out, f)):
                os.makedirs(os.path.dirname(os.path.join(out, f)), exist_ok=True); shutil.copy2(os.path.join(src, f), os.path.join(out, f))
    from wholecell.validate import common
    common.MAX_PROCS = a.nproc
    t0 = time.time()
    ctx = dict(bdir=bdir, out=out, nproc=a.nproc, quick=a.quick, complete=bool(a.complete), reuse=bool(a.reuse_from))
    V = dict(step="Step 2 combined", cell=a.cell, build=bdir, reference_build=REF[a.cell], reference_suite=REF_SUITE[a.cell],
             criteria="docs/unification_criteria/ (B6, B6b, B3, B4, B1, B2, B5 + addenda)", complete=bool(a.complete), quick=a.quick,
             date=time.strftime("%Y-%m-%d %H:%M"))
    try:
        V["no_regression"] = combined_no_regression(a.cell, bdir, ctx)
    except Exception as e:
        V["no_regression"] = dict(status="ERROR", error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc()[-2000:])
    acc = load_accepted(a.accepted, bdir)
    V["accepted_files"] = [os.path.relpath(os.path.abspath(p), REPO) for p in a.accepted]
    if acc:
        apply_accepted_regressions(V["no_regression"], acc)
    log("no_regression", V["no_regression"]["status"], {k: len(V["no_regression"].get(k, [])) for k in ("regressions_unexpected", "regressions_expected")})
    V["steps"] = {}
    for name, fn in SECTIONS:
        if name == "B3" and a.cell != "jurkat":
            V["steps"][name] = dict(status="NA", items=[], reason="B3 is the Jurkat library build mode (yeast: unchanged library path)")
            continue
        t1 = time.time()
        try:
            s = fn(a.cell, bdir, ctx)
            if acc:
                apply_accepted_items({name: s}, acc)
            s["status"] = step_status(s["items"])
        except Exception as e:
            s = dict(status="ERROR", items=[], error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc()[-3000:])
        s["seconds"] = round(time.time() - t1, 1)
        V["steps"][name] = s
        log(name, s["status"], [(i["id"], i["status"]) for i in s["items"]], f"{s['seconds']} s")
        json.dump(V, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)     # partial results survive a crash
    V["na_items"] = [f"{n}: {i['id']}: {i['reason']}" for n, s in V["steps"].items() for i in s["items"] if i["status"] == "NA"]
    V["failing"] = [f"{n}: {i['id']}" for n, s in V["steps"].items() for i in s["items"] if i["status"] == "FAIL"] + \
                   [f"{n}: ERROR" for n, s in V["steps"].items() if s["status"] == "ERROR"] + \
                   ([f"no_regression: {V['no_regression']['status']} ({len(V['no_regression'].get('regressions_unexpected', []))} unexpected; "
                     f"{V['no_regression'].get('reason') or V['no_regression'].get('error') or ''})"]
                    if V["no_regression"]["status"] not in ("PASS", "ACCEPTED_DEVIATION") else [])
    V["accepted_deviations"] = [f"{n}: {i['id']}: {i['reason']}" for n, s in V["steps"].items() for i in s["items"] if i["status"] == "ACCEPTED_DEVIATION"] + \
        [f"no_regression: {r['check']} {r['metric']} = {r['value']} (reference {r['reference']}): {r['explained']}"
         for r in V["no_regression"].get("regressions_accepted", [])]
    V["verdict_rule"] = "ACCEPTED if no FAIL item, no ERROR and no_regression PASS or ACCEPTED_DEVIATION (every remaining failure an accepted deviation)"
    V["notes"] = ctx.get("notes", [])
    if not a.complete or not os.path.exists(os.path.join(bdir, "validate", "suite.json")):
        V["verdict"] = "INCOMPLETE"
    else:
        V["verdict"] = "ACCEPTED" if not V["failing"] else "NOT ACCEPTED"
    V["steps_status"] = {n: s["status"] for n, s in V["steps"].items()}
    V["seconds"] = round(time.time() - t0, 1)
    json.dump(V, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    nr = V["no_regression"]
    with open(os.path.join(out, "regressions.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["class", "check", "metric", "value", "reference", "tolerance", "step", "explained"])
        for cls in ("regressions_unexpected", "regressions_accepted", "regressions_expected"):
            for r in nr.get(cls, []):
                w.writerow([cls.split("_")[1] + (f":{r['accepted_class']}" if r.get("accepted_class") else ""), r["check"], r["metric"], r["value"], r["reference"], r["tolerance"], r["step"], r["explained"]])
    log("VERDICT", V["verdict"], V["steps_status"], "failing", V["failing"], "->", out)


if __name__ == "__main__":
    main()
