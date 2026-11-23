# new in wholecell (Step 3 freeze preparation, 2026-11-19)
"""Freeze check for the paper builds (Jurkat v8, yeast next): docs/unification_criteria/freeze_v8.md (DRAFT).

Preflight (before a freeze build; run from a checkout of COMMIT):
    python scripts/freeze_check.py --preflight --cell jurkat|yeast --commit <sha> [--main-repo <repo>]
        [--criteria-md5 <md5>] [--waive gate=<decision> ...] [--satisfied gate=<evidence> ...] [--out <dir>]
  -> <out>/preflight.json; exit 0 only if every check passes (COMMIT on main and clean, required commits contained,
     criteria md5, config facts of §2b, pending gates of §2a satisfied or waived).

Verdict (on a finished freeze build):
    python scripts/freeze_check.py --cell jurkat|yeast --build <bdir> [--reference <bdir>] [--complete 0|1]
        [--commit <sha>] [--out <dir>] [--nproc 3] [--quick] [--accepted <json> ...] [--waive ...] [--satisfied ...]
        [--replicates 1]
  -> <out>/verdict.json (default <bdir>/validate/compare_freeze_v8/), <out>/regressions.tsv. Gates G0-G9 of the
     criteria §6. Nothing is written into the reference build.

Reuse: scripts/step2_combined_check.py (S2): the Step-2 step sections (G2), the accepted-deviation mechanism
(--accepted, same json format; an entry applies only to the builds it lists), C10 whole-cell tile sums, the B1 K6 /
B5 L6 density helpers; wholecell.validate.suite.no_regression (plan §4 tolerances) for G8.

Gates (criteria §6): G0 integrity | G1 pending inputs (secms_jurkat, surface, pm_topology, copy_fixes / secms_yeast:
satisfied in the config, or --satisfied gate=<evidence> (addendum), or --waive gate=<dated DECISION>) | G2 Step-2 items
(ACCEPTED / EXPECTED deviations of the Step-2 record carried with their limits; KNOWN not carried) | G3 orientation
(C9 = 0) | G4 nucleolus cap | G5 yeast B5 overrides | G6 SEC-MS | G7 surface / PM topology / copy fixes / ledger
rule_id | G8 no regression vs the freeze reference | G9 replicates (report).
Verdict: ACCEPTED (no FAIL / ERROR, no un-waived NA gate) | INCOMPLETE (build/suite incomplete or reference pending) |
NOT ACCEPTED. CPU: --nproc processes (default 3), NUMBA / OMP threads 3.
"""
import argparse
import contextlib
import copy
import csv
import hashlib
import json
import os
import re
import subprocess
import sys
import time
import traceback

for _k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS"):
    os.environ.setdefault(_k, "3")
os.environ.setdefault("PYTHONHASHSEED", "0")

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path[:0] = [os.path.join(REPO, "scripts", "dev"), os.path.join(REPO, "scripts")]
import step2_combined_check as S2  # noqa: E402

item, judged, jload, _f, sval, log = S2.item, S2.judged, S2.jload, S2._f, S2.sval, S2.log
RES = S2.RES
AV = 6.02214076e23

CRITERIA = "docs/unification_criteria/freeze_v8.md"
STEP2_RECORD = os.path.join(REPO, "docs/unification_criteria/step2_combined_accepted_deviations_2026-11-19.json")
# freeze reference = the accepted combined Step-2 build of the cell (criteria §1); yeast: pending its Step-2 verdict
FREEZE_REF = {"jurkat": f"{RES}/jurkat/20261118-40d7967-s20261115-step2", "yeast": None}
# commits every freeze COMMIT must contain (criteria §2a, merged inputs)
REQUIRED = {"jurkat": ["40d7967", "257ec59", "ab9927c", "3227713"],
            "yeast": ["40d7967", "257ec59", "ab9927c", "3227713", "b16bd09", "b691e9a", "6978b43", "a1fb253", "9df5a94"]}
GATES = {"jurkat": ("secms_jurkat", "surface", "pm_topology", "copy_fixes"), "yeast": ("secms_yeast",)}
# frozen inputs that may differ from the reference build (criteria §2c)
ALLOWED_FROZEN = {"jurkat": (r"^inv_sec_", r"^inv_surface_", r"^inv_copy_fixes"), "yeast": (r"^inv_sec_", r"^inv_secms_")}
GEOM_SHA = {"jurkat": "6c86b1688722f4320520abbb46a66fb7e9e97c95576e0186ebf3909b9c562fa1",
            "yeast": "230b96d91559b6aca5557bae68cfafeed227c24be3aaa76879831acd29dd8f5b"}
# E2 / E3 / E4 / E8 bounds (criteria §4)
NUCLEOLUS_MIN = {"jurkat": 0.86, "yeast": None}          # yeast: >= reference - 0.003
E3_YEAST = [("C1", "SWAT_native.intended.tophit.unweighted", 0.6682 - 0.005),
            ("C1", "hyperLOPIT.intended.tophit.unweighted", 0.8483 - 0.005),
            ("C3", "image_nc.test.spearman", 0.8345 - 0.01)]
# E4 (revision 2026-11-19 evening): wiring caee4f2 + HSPD1 apex_check_exclude (step3-secms 906b9b9), dry run
# jurkat/_dev_unify/step3_secms_hspd1/dryrun.json (wiring alone: 0.5352 / 2.72 pg)
E4_JURKAT = dict(ref_share=0.5253, new_share=0.5371, share_tol=0.001, moved_pg=2.53, moved_tol=0.10, cp_tol=1e-6, cp_tol_fixes=0.05)
# E8 (revision 2026-11-19 evening): PDGFRB capped to its RNA-based estimate (copies = round(estimate in copies));
# SOS1 / SOS2 NOT raised: v0 copies = reference (the post-scaling adjustment leaves every other protein unchanged)
E8 = dict(pdgfrb_max=100, sos_abs=0)
HSPD1 = "P10809"                                            # E4: apex_check_exclude (heptamer kept)
CARRIED_CLASSES = ("ACCEPTED", "EXPECTED")                 # Step-2 record classes carried to the freeze (KNOWN are fixed)


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def git(*a, repo=REPO):
    r = subprocess.run(["git", "-C", repo, *a], capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def kv(lst):
    out = {}
    for s in lst or []:
        k, _, v = s.partition("=")
        out[k.strip()] = v.strip() or "(no reference given)"
    return out


def dig(d, *path):
    for p in path:
        if not isinstance(d, dict):
            return None
        d = d.get(p)
    return d


# ------------------------------------------------------------------------------------------------ config facts / gates
def config_facts(cfg, cell):
    """criteria §2b: list of (name, ok, value)"""
    inv, pk = cfg.get("inventory") or {}, cfg.get("pack") or {}
    par = inv.get("params") or {}
    org = pk.get("organisation") or {}
    loc = par.get("localisation") or {}
    out = [("geometry manifest_sha256", dig(cfg, "geometry", "manifest_sha256") == GEOM_SHA[cell], dig(cfg, "geometry", "manifest_sha256")),
           ("organisation enabled", bool(org.get("enabled")) if isinstance(org, dict) else bool(org), org.get("enabled") if isinstance(org, dict) else org),
           ("field_cap_regions = [nucleolus]", list(org.get("field_cap_regions") or []) == ["nucleolus"], org.get("field_cap_regions")),
           ("localisation preset common", loc.get("preset") == "common", loc.get("preset")),
           ("build.allow_dirty false", not dig(cfg, "build", "allow_dirty"), dig(cfg, "build", "allow_dirty"))]
    if cell == "jurkat":
        mods = [m for m in (dig(cfg, "preplaced", "serial", "builders") or []) if isinstance(m, dict) and m.get("name") == "npc"]
        npc = (mods[0].get("params") or {}).get("n_pores") if mods else None
        out += [("pre-placed npc n_pores = 2696 (B4 J1)", npc == 2696, npc),
                ("assembly preset legacy_v7", dig(par, "assembly_steps", "preset") == "legacy_v7", dig(par, "assembly_steps", "preset")),
                ("no pack preset pin (B1 COMMON_DEFAULTS = {})", not pk.get("preset"), pk.get("preset")),
                ("library mode build (B3)", dig(cfg, "library", "mode") == "build", dig(cfg, "library", "mode")),
                ("export YCT2 tiles enabled (DECISION 2026-11-19)", dig(cfg, "export", "tiles", "enabled") is True and
                 (dig(cfg, "export", "enabled") is not False), dig(cfg, "export", "tiles", "enabled"))]
    else:
        out += [("pack preset [legacy_x5, b1_s4y] (B1 C3 pin)", list(pk.get("preset") or []) == ["legacy_x5", "b1_s4y"], pk.get("preset")),
                ("weight_overrides (B5 yeast overrides)", bool(dig(loc, "common", "weight_overrides")), bool(dig(loc, "common", "weight_overrides")))]
    return out


def gate_state(cfg, cell):
    """criteria §2a pending inputs: {gate: (satisfied_in_config, evidence)}"""
    inv = cfg.get("inventory") or {}
    par = inv.get("params") or {}
    sec = dig(par, "assembly_steps", "overrides", "secms", "source_id")
    st = {}
    if cell == "jurkat":
        ex = [str(r.get("uniprot")) for r in (dig(par, "assembly_steps", "data", "apex_check_exclude") or []) if isinstance(r, dict)]
        st["secms_jurkat"] = (sec == "jurkat_final" and HSPD1 in ex,
                              f"assembly_steps.overrides.secms.source_id = {sec}; apex_check_exclude = {ex} (HSPD1 {HSPD1} required)")
        smode = dig(par, "surface_proteomics", "mode") or "presence"
        sad = dig(inv, "adapters", "surface")
        st["surface"] = (bool(sad) and smode == "presence",
                         f"inventory.adapters.surface = {sad}; surface_proteomics.mode = {smode} (presence only: DECISION 2026-11-19)")
        fx = [x for x in (dig(cfg, "library", "steps") or []) if isinstance(x, dict) and x.get("step") == "pm_topology_fix"]
        st["pm_topology"] = (dig(cfg, "library", "mode") == "build" and len(fx) == 1 and bool(fx[0].get("table")),
                             f"library.steps pm_topology_fix = {fx} (freeze v8 E7, wholecell/library/topology.py)")
        adj = [s.get("rule_id") for s in (dig(inv, "copies", "adjustments") or []) if isinstance(s, dict)]
        st["copy_fixes"] = ("v8.copy_fixes" in adj, f"inventory.copies.adjustments rule ids = {adj}")
    else:
        st["secms_yeast"] = (sec == "yeast_final" and bool(dig(inv, "adapters", "secms")),
                             f"adapters.secms = {dig(inv, 'adapters', 'secms')}; assembly_steps.overrides.secms.source_id = {sec} "
                             "(other key in the yeast SEC-MS commit: pass --satisfied secms_yeast=<commit>)")
    return st


def g1_items(cfg, cell, waive, satisfied):
    items, open_ = [], []
    for g, (ok, ev) in gate_state(cfg, cell).items():
        if ok:
            items.append(item(f"G1 {g}", "PASS", ev, "pending input merged (config)"))
        elif g in satisfied:
            items.append(item(f"G1 {g}", "PASS", ev, "pending input merged", f"--satisfied: {satisfied[g]}"))
        elif g in waive:
            items.append(item(f"G1 {g}", "REPORT", ev, "pending input", f"WAIVED: {waive[g]}"))
        else:
            items.append(item(f"G1 {g}", "NA", ev, "pending input merged or waived by a dated DECISION", "open (blocks the freeze)"))
            open_.append(g)
    return items, open_


def gate_on(g, cfg, cell, satisfied):
    st = gate_state(cfg, cell).get(g)
    return bool(st and st[0]) or g in satisfied


# ------------------------------------------------------------------------------------------------ preflight
def preflight(a):
    from wholecell import config as wcfg
    res = dict(cell=a.cell, commit=a.commit, date=time.strftime("%Y-%m-%d %H:%M"), checks=[])
    chk = res["checks"]
    rc, head = git("rev-parse", "HEAD")
    rc2, full = git("rev-parse", a.commit)
    chk.append(dict(id="checkout = COMMIT", ok=rc == 0 and rc2 == 0 and head == full, value=dict(head=head, commit=full)))
    rc, st = git("status", "--porcelain", "--untracked-files=no")
    chk.append(dict(id="checkout clean", ok=rc == 0 and st == "", value=st[:500]))
    main_repo = a.main_repo or REPO
    rc, _ = git("merge-base", "--is-ancestor", full or a.commit, "main", repo=main_repo)
    chk.append(dict(id="COMMIT on main", ok=rc == 0, value=main_repo))
    for c in REQUIRED[a.cell]:
        rc, _ = git("merge-base", "--is-ancestor", c, full or a.commit)
        chk.append(dict(id=f"contains {c}", ok=rc == 0, value=None))
    crit = os.path.join(REPO, CRITERIA)
    m = md5(crit) if os.path.exists(crit) else None
    chk.append(dict(id="criteria md5", ok=(m is not None) and (a.criteria_md5 is None or m == a.criteria_md5),
                    value=dict(file=CRITERIA, md5=m, expected=a.criteria_md5)))
    cfg = wcfg.load(a.cell)
    for name, ok, v in config_facts(cfg, a.cell):
        chk.append(dict(id=f"config: {name}", ok=bool(ok), value=v))
    waive, sat = kv(a.waive), kv(a.satisfied)
    gi, open_ = g1_items(cfg, a.cell, waive, sat)
    for i in gi:
        chk.append(dict(id=i["id"], ok=i["status"] in ("PASS", "REPORT"), value=i["value"], note=i["reason"]))
    if a.replicates != 1:
        chk.append(dict(id="replicates", ok=False, value=a.replicates,
                        note="1 seed per cell (DECISION 2026-11-19 evening, criteria §0 D1 / §3): no replicate builds"))
    res["waived"], res["satisfied"], res["open_gates"] = waive, sat, open_
    res["ok"] = all(c["ok"] for c in chk)
    out = a.out or os.getcwd()
    os.makedirs(out, exist_ok=True)
    json.dump(res, open(os.path.join(out, "preflight.json"), "w"), indent=1, default=str)
    for c in chk:
        log("PREFLIGHT", "ok " if c["ok"] else "FAIL", c["id"], c.get("note", "") or "")
    log("PREFLIGHT", "PASS" if res["ok"] else "FAIL", "->", os.path.join(out, "preflight.json"))
    return 0 if res["ok"] else 1


# ------------------------------------------------------------------------------------------------ helpers
def frozen_rows(bdir):
    p = os.path.join(bdir, "frozen_inputs.tsv")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return {r["name"]: r for r in csv.DictReader(f, delimiter="\t")}


def c10_of(cell, bdir, ctx, tag):
    """C10 whole-cell tile sums of a build: its own Step-2 / freeze check output if present, else computed"""
    for d in ("compare_freeze_v8", "compare_step2_v3", "compare_step2_v2", "compare_step2"):
        p = os.path.join(bdir, "validate", d, "c10_build.json")
        if os.path.exists(p) and os.path.realpath(os.path.dirname(p)) != os.path.realpath(ctx["out"]):
            return jload(p), p
    return S2._c10_whole_cell(cell, bdir, ctx, tag), "computed"


@contextlib.contextmanager
def s2_reference(cell, ref):
    """S2 helpers (density changes) read S2.REF[cell]: point them at the freeze reference for the duration"""
    old = S2.REF[cell]
    S2.REF[cell] = ref
    try:
        yield
    finally:
        S2.REF[cell] = old


def complexed(inv):
    """SEC-MS dry-run metric (scripts/secms_dataset_dryrun.py): assembled mass by category and complexed share"""
    import pandas as pd
    rd = lambda f: pd.read_csv(os.path.join(inv, "assembly", f), sep="\t")
    tot = float(rd("inventory_v0_complexed.tsv").mass_pg.sum())
    d = {"complex_portal": float(rd("complexes_v0.tsv").mass_built_pg.sum())}
    for f, c in (("homo_oligomers_v0.tsv", "homo_oligomers"), ("afdb_homodimers_v0.tsv", "afdb_homodimers"),
                 ("afdb_heterodimers_v0.tsv", "afdb_heterodimers")):
        d[c] = float(rd(f).mass_pg.sum())
    d["complexed_total"] = sum(d.values())
    d["total_pg"] = tot
    d["complexed_share"] = d["complexed_total"] / tot
    return d


def secms_moved(ref_inv, new_inv):
    import secms_dataset_dryrun as D
    Pr, Pn = D.per_protein_copies(ref_inv), D.per_protein_copies(new_inv)
    cats = ["complex_portal", "homo_oligomers", "afdb_homodimers", "afdb_heterodimers"]
    ar = Pr[cats].sum(1)
    an = Pn[cats].sum(1).reindex(Pr.index).fillna(0)
    dm = (an - ar) * Pr.mass_da.fillna(0) / AV * 1e12
    ch = dm[(an - ar).abs() > 0.5]
    return dict(n_changed=int(len(ch)), moved_pg=float(ch.abs().sum()), to_assembled_pg=float(ch.clip(lower=0).sum()),
                to_free_pg=float(-ch.clip(upper=0).sum()))


def pdgfrb_estimate(bdir):
    """(RNA-based estimate in copies, copies) of PDGFRB in a build: the copy step's undetected-`low`-call value
    min(10^(a + b log10 TPM), floor) x low_factor (ppb; fit from inventory/copies/summary.json, low_factor 0.25 = the
    cell config) at the build's scale: ppb / sum(ppb) x molecules (summary total_protein)"""
    import math
    import pandas as pd
    S = jload(os.path.join(bdir, "inventory", "copies", "summary.json")) or {}
    im = S.get("imputation") or {}
    p = os.path.join(bdir, "inventory", "copies", "copies.tsv")
    if not im.get("slope") or not os.path.exists(p):
        return None, None
    T = pd.read_csv(p, sep="\t", usecols=["gene", "tpm", "ppb", "copies"]).set_index("gene")
    if "PDGFRB" not in T.index:
        return None, None
    tpm = float(T.tpm["PDGFRB"])
    est = min(10 ** (im["intercept"] + im["slope"] * math.log10(max(tpm, 1e-3))), im["floor"]) * 0.25
    n = float((S.get("total_protein") or {}).get("molecules") or 0)
    return est / float(T.ppb.fillna(0).sum()) * n, float(T.copies["PDGFRB"])


def copies_of(bdir, genes):
    import pandas as pd
    p = os.path.join(bdir, "inventory", "copies", "copies.tsv")
    if not os.path.exists(p):
        return {}
    T = pd.read_csv(p, sep="\t", usecols=["gene", "copies"]).set_index("gene").copies
    return {g: (float(T[g]) if g in T.index else None) for g in genes}


# ------------------------------------------------------------------------------------------------ gates
def g0(cell, bdir, ref, cfg, a, ctx):
    items = []
    items.append(judged("G0 build complete", bool(a.complete), a.complete, "snakemake 100 % (freeze_build.sh)"))
    items.append(judged("G0 suite.json", os.path.exists(os.path.join(bdir, "validate", "suite.json")), None, "validate/suite.json present"))
    bs = jload(os.path.join(bdir, "BUILD_STAMP.json")) or {}
    same = a.commit is None or bool(bs.get("commit")) and (bs["commit"].startswith(a.commit) or a.commit.startswith(bs["commit"]))
    items.append(judged("G0 clean commit", bs.get("dirty") is False and same,
                        dict(commit=bs.get("commit"), dirty=bs.get("dirty"), expected=a.commit), "BUILD_STAMP dirty false; commit = COMMIT"))
    rc, _ = git("merge-base", "--is-ancestor", bs.get("commit") or "none", "main", repo=a.main_repo or REPO)
    items.append(judged("G0 commit on main", rc == 0, bs.get("commit"), "the build commit is an ancestor of main"))
    facts = config_facts(cfg, cell)
    items.append(judged("G0 config facts", all(ok for _, ok, _ in facts), [dict(metric=n, value=v, ok=bool(ok)) for n, ok, v in facts],
                        "criteria §2b", "config.resolved.json"))
    fb = frozen_rows(bdir)
    fr = frozen_rows(ref) if ref else None
    if fb is None:
        items.append(item("G0 frozen inputs", "FAIL", None, "frozen_inputs.tsv present, all ok", "missing"))
    else:
        bad_ok = sorted(k for k, r in fb.items() if str(r.get("ok")) != "True")
        if fr is None:
            items.append(item("G0 frozen inputs", "NA" if not bad_ok else "FAIL", dict(not_ok=bad_ok), "= reference ± §2c allowances",
                              "reference pending: only the per-file checks judged"))
        else:
            sig = lambda r: (r.get("algo"), r.get("expected"))
            diff = sorted(k for k in set(fb) | set(fr) if k not in fb or k not in fr or sig(fb[k]) != sig(fr[k]))
            allowed = [k for k in diff if any(re.search(p, k) for p in ALLOWED_FROZEN[cell])]
            other = [k for k in diff if k not in allowed]
            items.append(judged("G0 frozen inputs", not bad_ok and not other, dict(not_ok=bad_ok, allowed_changes=allowed, other_changes=other),
                                "every file ok; differences vs the reference only in the §2c allowances", "frozen_inputs.tsv"))
    if cell == "jurkat":
        ec = jload(os.path.join(bdir, "export", "web", "encode_check.json"))
        keys = ("count_equal", "pos_ok", "rot_ok")
        items.append(judged("G0 YCT2 tiles", None if ec is None else all(ec.get(k) is True for k in keys) and ec.get("blocks") == "all",
                            None if ec is None else {k: ec.get(k) for k in keys + ("n_tiles", "instances_source", "instances_decoded", "blocks")},
                            "export/web/encode_check.json: whole cell (blocks all), count_equal, pos_ok, rot_ok", "export/web"))
    S = ctx.get("suite_new") or {}
    c16 = [("manifest_stale_files", 0), ("depends_on_geometry_mismatch", 0)] + ([("yeast.stale_result_files", 0)] if cell == "yeast" else [])
    vals = [dict(metric=m, value=sval(S, "C16", m), ok=(sval(S, "C16", m) or 0) <= lim) for m, lim in c16]
    items.append(judged("G0 C16 build integrity", None if not S else all(v["ok"] for v in vals), vals, "stale 0, geometry mismatch 0 / None", "suite.json C16"))
    return items


def g2(cell, bdir, a, ctx, acc_extra):
    """Step-2 step sections on the freeze build (their own references) + carried Step-2 deviations (items; the carried
    regression entries are vs the B3 reference of Step 2 and do not apply to G8, whose reference is the Step-2 build)"""
    rec = json.load(open(STEP2_RECORD)) if cell == "jurkat" else {"accepted": []}
    carried = [dict(e, _file=os.path.relpath(STEP2_RECORD, REPO), _source=rec.get("decision_source", ""))
               for e in rec["accepted"] if e.get("class") in CARRIED_CLASSES]
    steps = {}
    try:
        S2.combined_no_regression(cell, bdir, ctx)       # fills ctx suite_new / suite_ref (Step-2 reference); not gated here
    except Exception as e:
        ctx.setdefault("notes", []).append(f"S2 no_regression (context only): {type(e).__name__}: {e}")
    for name, fn in S2.SECTIONS:
        if name == "B3" and cell != "jurkat":
            steps[name] = dict(status="NA", items=[], reason="B3 is Jurkat only")
            continue
        t1 = time.time()
        try:
            s = fn(cell, bdir, ctx)
            S2.apply_accepted_items({name: s}, carried + [e for e in acc_extra if e["kind"] == "item"])
            s["status"] = S2.step_status(s["items"])
        except Exception as e:
            s = dict(status="ERROR", items=[], error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc()[-3000:])
        s["seconds"] = round(time.time() - t1, 1)
        steps[name] = s
        log("G2", name, s["status"], [(i["id"], i["status"]) for i in s["items"]], f"{s['seconds']} s")
    return steps


def g3_g4(cell, bdir, ref, ctx):
    S, R = ctx["suite_new"], ctx.get("suite_freeze_ref") or {}
    items = []
    se = sval(S, "C9", "sign_errors_per_1e6")
    items.append(judged("G3 E1 C9 sign errors = 0", None if se is None else se == 0, dict(build=se, reference=sval(R, "C9", "sign_errors_per_1e6")),
                        "C9 sign_errors_per_1e6 = 0 (orientation fix 3227713)", "suite.json C9"))
    mt, mr = sval(S, "C9", "median_tilt_deg"), sval(R, "C9", "median_tilt_deg")
    items.append(judged("G3 E1 C9 median tilt", None if mt is None or mr is None else mt <= mr + 0.5, dict(build=mt, reference=mr),
                        "<= reference + 0.5 deg", "suite.json C9"))
    cb, cr = ctx.get("c10_build"), ctx.get("c10_ref")
    if not cb:
        items.append(item("G4 E2 nucleolus", "NA", None, "nucleolus placed mass fraction", "tile sums not computable"))
        return items
    nb = (cb["regions"].get("nucleolus") or {}).get("placed_mass_frac")
    nr = ((cr or {}).get("regions", {}).get("nucleolus") or {}).get("placed_mass_frac")
    lim = NUCLEOLUS_MIN[cell] if NUCLEOLUS_MIN[cell] is not None else (nr - 0.003 if nr is not None else None)
    items.append(judged("G4 E2 nucleolus placed fraction", None if nb is None or lim is None else nb >= lim,
                        dict(build=nb, reference=nr, bound=lim), ">= 0.86 (Jurkat) / >= reference - 0.003 (yeast)", "C10 tile sums"))
    if cr:
        ob, orf = cb.get("placed_mass_fraction"), cr.get("placed_mass_fraction")
        items.append(judged("G4 E1/E2 overall placed fraction", None if ob is None or orf is None else ob >= orf - 0.003,
                            dict(build=ob, reference=orf), ">= reference - 0.003", "C10 tile sums"))
        worse = {k: dict(build=v["placed_mass_frac"], reference=(cr["regions"].get(k) or {}).get("placed_mass_frac"))
                 for k, v in cb["regions"].items() if k != "nucleolus" and v.get("placed_mass_frac") is not None
                 and (cr["regions"].get(k) or {}).get("placed_mass_frac") is not None
                 and v["placed_mass_frac"] < cr["regions"][k]["placed_mass_frac"] - 0.003}
        items.append(judged("G4 E2 other regions", not worse, [dict(metric=k, value=v["build"], reference=v["reference"], ok=False)
                                                              for k, v in worse.items()], "every other region >= reference - 0.003", "C10 tile sums"))
    else:
        items.append(item("G4 E1/E2 overall placed fraction", "NA", dict(build=cb.get("placed_mass_fraction")), ">= reference - 0.003", "reference pending"))
    return items


def g5(cell, ctx):
    if cell != "yeast":
        return []
    S = ctx["suite_new"]
    vals = [dict(metric=f"{c}.{m}", value=sval(S, c, m), bound=b, ok=None if sval(S, c, m) is None else sval(S, c, m) >= b) for c, m, b in E3_YEAST]
    return [judged("G5 E3 yeast B5 overrides (§3d)", None if any(v["ok"] is None for v in vals) else all(v["ok"] for v in vals), vals,
                   "SWAT native >= 0.6632, hyperLOPIT >= 0.8433, image N/C test >= 0.8245", "suite.json C1 / C3")]


def g6(cell, bdir, ref, cfg, satisfied, waive):
    items = []
    gname = "secms_jurkat" if cell == "jurkat" else "secms_yeast"
    on = gate_on(gname, cfg, cell, satisfied)
    if cell == "yeast":
        items.append(item("G6 E5 yeast SEC-MS", "NA" if on else "REPORT", None, "bounds TBD (yeast SEC-MS addendum)",
                          "merged, bounds not yet written" if on else f"gate not merged ({'waived: ' + waive[gname] if gname in waive else 'see G1'})"))
        return items
    try:
        nb = complexed(os.path.join(bdir, "inventory"))
        nr = complexed(os.path.join(ref, "inventory")) if ref else None
        mv = secms_moved(os.path.join(ref, "inventory"), os.path.join(bdir, "inventory")) if ref else None
    except Exception as e:
        return [item("G6 E4 Jurkat SEC-MS", "NA", None, "complexed share / mass moved", f"{type(e).__name__}: {e}")]
    val = dict(build=nb, reference=nr, moved=mv)
    if not on:
        items.append(item("G6 E4 Jurkat SEC-MS", "REPORT", val, "E4 bounds apply only with the SEC-MS input merged",
                          f"gate secms_jurkat not merged ({'waived: ' + waive[gname] if gname in waive else 'see G1'})"))
        return items
    fixes = gate_on("copy_fixes", cfg, cell, satisfied)
    E = E4_JURKAT
    ok_ref = nr is not None and abs(nr["complexed_share"] - E["ref_share"]) <= E["share_tol"]
    ok_new = abs(nb["complexed_share"] - E["new_share"]) <= E["share_tol"]
    ok_mv = mv is not None and abs(mv["moved_pg"] - E["moved_pg"]) <= E["moved_tol"]
    ok_cp = nr is not None and abs(nb["complex_portal"] - nr["complex_portal"]) <= (E["cp_tol_fixes"] if fixes else E["cp_tol"])
    items.append(judged("G6 E4 complexed share", ok_ref and ok_new, dict(build=nb["complexed_share"], reference=(nr or {}).get("complexed_share")),
                        f"reference {E['ref_share']} and build {E['new_share']} (± {E['share_tol']})", "inventory/assembly"))
    items.append(judged("G6 E4 mass moved", ok_mv, mv, f"{E['moved_pg']} ± {E['moved_tol']} pg free <-> assembled", "inventory/assembly"))
    items.append(judged("G6 E4 Complex Portal unchanged", ok_cp, dict(build=nb["complex_portal"], reference=(nr or {}).get("complex_portal")),
                        "|delta| <= 1e-6 pg (0.05 pg with the E8 copy fixes)", "assembly/complexes_v0.tsv"))
    return items


def g7(cell, bdir, ref, cfg, satisfied, waive):
    items = []
    # ledger rule_id (E9, backlog 5c): no numeric rule ids in version_changes.tsv
    p = os.path.join(bdir, "inventory", "version_changes.tsv")
    if os.path.exists(p):
        n = bad = 0
        with open(p) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                if r.get("rule_ids"):
                    n += 1
                    bad += any(re.fullmatch(r"-?[0-9.eE+-]+", x) for x in r["rule_ids"].split(","))
        items.append(judged("G7 E9 ledger rule_id", bad == 0, dict(rows_with_rule_ids=n, numeric=bad), "0 numeric rule ids", "inventory/version_changes.tsv"))
    else:
        items.append(item("G7 E9 ledger rule_id", "NA", None, "0 numeric rule ids", "version_changes.tsv missing"))
    if cell != "jurkat":
        return items
    cp = copies_of(bdir, ["PDGFRB", "SOS1", "SOS2"])
    if gate_on("copy_fixes", cfg, cell, satisfied):
        est, n = pdgfrb_estimate(bdir)
        ok_p = est is not None and n is not None and n <= E8["pdgfrb_max"] and n == round(est)
        items.append(judged("G7 E8 PDGFRB at its RNA-based estimate", ok_p, dict(copies=n, rna_estimate_copies=est),
                            f"copies = round(min(fit, floor) x low_factor / sum(ppb) x molecules) of this build; <= {E8['pdgfrb_max']}",
                            "inventory/copies/{copies.tsv, summary.json}"))
        rc = copies_of(ref, ["SOS1", "SOS2"]) if ref else {}
        ok_s = bool(rc) and all(cp.get(g) is not None and rc.get(g) is not None and abs(cp[g] - rc[g]) <= E8["sos_abs"] for g in ("SOS1", "SOS2"))
        items.append(judged("G7 E8 SOS not raised", ok_s if rc else None, dict(build={g: cp.get(g) for g in ("SOS1", "SOS2")}, reference=rc),
                            "SOS1 / SOS2 copies = reference (Jurkat measured value kept; simulations use 2,080 and the primary-T value explicitly)",
                            "inventory/copies/copies.tsv"))
    else:
        items.append(item("G7 E8 copy fixes", "REPORT", cp, "decision needed (criteria §5b)",
                          "waived: " + waive["copy_fixes"] if "copy_fixes" in waive else "not merged (see G1)"))
    # E6 surface proteomics, presence evidence only (revision 2026-11-19 evening): surface.tsv filled, pm_fraction empty;
    # copies / localisation untouched (judged by G7 E8 SOS, G6 Complex Portal and G8; the adapter writes surface.tsv only)
    if gate_on("surface", cfg, cell, satisfied):
        import pandas as pd
        p = os.path.join(bdir, "inventory", "surface.tsv")
        T = pd.read_csv(p, sep="\t") if os.path.exists(p) else None
        sid = dig(cfg, "inventory", "params", "surface_proteomics", "source_id")
        val = None if T is None else dict(rows=len(T), proteins=int(T.uniprot.nunique()), source_ids=sorted(T.source_id.dropna().unique().tolist()),
                                          pm_fraction_set=int(T.pm_fraction.notna().sum()))
        ok = T is not None and len(T) > 0 and val["source_ids"] == [sid] and val["pm_fraction_set"] == 0 and val["proteins"] <= 596
        items.append(judged("G7 E6 surface presence evidence", ok, val, "surface.tsv rows of the configured source only; "
                            "<= 596 proteins; pm_fraction empty (presence mode: no PM copy model)", "inventory/surface.tsv"))
    else:
        items.append(item("G7 E6 surface", "REPORT", None, "presence evidence (adapter mode presence)",
                          "waived: " + waive["surface"] if "surface" in waive else "not merged (see G1)"))
    on = gate_on("pm_topology", cfg, cell, satisfied)
    if on and cell == "jurkat":
        # E7 (fix written 2026-11-19, wholecell/library/topology.py): library diff vs the reference library + S0 orientation rule
        from wholecell.library import topology as TP
        fx = [x for x in (dig(cfg, "library", "steps") or []) if isinstance(x, dict) and x.get("step") == "pm_topology_fix"]
        lib, rlib = os.path.join(bdir, "library", "lib"), (os.path.join(ref, "library", "lib") if ref else None)
        tab = os.path.join(bdir, "_code", fx[0]["table"]) if fx else None
        if not (fx and os.path.isdir(lib) and rlib and os.path.isdir(rlib) and os.path.exists(tab)):
            items.append(judged("G7 E7 pm_topology", None, dict(lib=lib, reference_lib=rlib, table=tab),
                                "criteria §4 E7", "library / reference library / fix table missing"))
        else:
            res, rows = TP.check(lib, rlib, os.path.join(bdir, "inventory", "pack_table.tsv"), tab)
            import pandas as pd
            out_dir = os.path.join(bdir, "validate", "compare_freeze_v8")
            os.makedirs(out_dir, exist_ok=True)
            pd.DataFrame(rows).to_csv(os.path.join(out_dir, "e7_pm_units.tsv"), sep="\t", index=False)
            json.dump(res, open(os.path.join(out_dir, "e7_pm_topology.json"), "w"), indent=1, default=str)
            items.append(judged("G7 E7 only the fix units changed (exact rotations of the reference)", res["only_fix_units_changed"],
                                dict(changed=res["changed_units"], errors=res["rotation_errors"]),
                                "library npz md5 = reference except the rotated units; each = reference rotated 180 deg about x",
                                "library/lib vs reference library/lib"))
            items.append(judged("G7 E7 glyco units not flagged", not res["glyco_flagged"] and not res["glyco_missing"],
                                dict(units=res["glyco_units"], copies=res["glyco_copies"], flagged=res["glyco_flagged"], missing=res["glyco_missing"]),
                                "S0 orientation rule: 0 of the 13 glyco-captured units flagged", "validate/compare_freeze_v8/e7_pm_topology.json"))
            items.append(judged("G7 E7 TRBC1 by TM alternation", not res["alternation_flagged"], res["alternation_units"],
                                "tm_alternation units not flagged when TOPO_DOM gaps are filled by alternation", "e7_pm_units.tsv"))
            items.append(judged("G7 E7 residual flags unchanged, none new", res["residual_unchanged"] and not res["new_flags"],
                                dict(residual=len(res["residual_flags"]), new=res["new_flags"], reference_flagged=res["ref_flagged"]),
                                "other reference flags unchanged (no glyco evidence; left as known residuals); no new flag", "e7_pm_units.tsv"))
            items.append(item("G7 E7 known residuals (listed)", "REPORT", dict(units=res["residual_flags"], copies=res["residual_copies"]),
                              "flagged PM units without glyco evidence, left unchanged (DECISION 2026-11-19)", "", "e7_pm_units.tsv"))
    else:
        items.append(item("G7 E7 pm_topology", "REPORT", None,
                          "criteria §4 E7: the 13 glyco-captured units re-oriented (topology check: 0 of them flipped); others unchanged",
                          "waived: " + waive["pm_topology"] if "pm_topology" in waive else "not merged (see G1)"))
    return items


def g8(cell, bdir, ref, cfg, satisfied, ctx, acc_extra):
    from wholecell.validate.suite import no_regression
    S = jload(os.path.join(bdir, "validate", "suite.json"))
    rsuite = os.path.join(ref, "validate", "suite.json") if ref else None
    R = jload(rsuite) if rsuite else None
    if S is None or R is None:
        return dict(status="NA", reason=f"suite.json missing (build {S is not None}; reference {rsuite or 'pending'})")
    table, _ = no_regression([copy.deepcopy(r) for r in S["metrics"]], R["metrics"], rsuite)
    species_change = any(gate_on(g, cfg, cell, satisfied) for g in ("secms_jurkat", "surface", "copy_fixes", "pm_topology", "secms_yeast"))
    band = {}
    if cell == "yeast":
        for sp in S2.YEAST_SEED_SUITES:
            for r in json.load(open(sp))["metrics"]:
                x = _f(r.get("value"))
                if x is not None:
                    band.setdefault((r["check"], r["metric"]), []).append(x)
    c2 = dict(ctx, out=os.path.join(ctx["out"], "g8"), c10={"build": ctx.get("c10_build"), "reference": ctx.get("c10_ref")})
    for k in ("placed_density_change", "b5_density_intended"):
        c2.pop(k, None)
    os.makedirs(c2["out"], exist_ok=True)

    def why(r):
        c, m = r["check"], str(r["metric"])
        val, rv, tol = _f(r["value"]), _f(r["reference_build_value"]), _f(r["tolerance"]) or 0.0
        if c == "C10" and m in ("placed_copies_fraction", "species_below_50pct_placed") and species_change:
            return "E4/E6/E8", "species set changed by a merged freeze input (SEC-MS / surface / copy fixes); placed mass judged by G4"
        if c == "C11" and m.startswith("density."):
            reg = m.split(".", 1)[1]
            if cell == "jurkat":
                with s2_reference(cell, ref):
                    dp = S2._placed_density_change(cell, c2).get(reg)
                if dp is not None and val is not None and rv is not None and abs((val - rv) - dp) <= 2.0:
                    return "E2" if reg == "nucleolus" else "K6 mass-consistent", f"mass-consistent with the placed-mass change {dp:+.2f} mg/ml (B1 K6 rule)"
            if species_change:
                with s2_reference(cell, ref):
                    d = c2.get("b5_density_intended")
                    if d is None:
                        d = c2["b5_density_intended"] = S2._intended_density(cell, c2)
                di = d.get(reg) if isinstance(d, dict) else None
                if di is not None and val is not None and rv is not None and abs((val - rv) - di) <= 2.0:
                    return "E4/E6/E8", f"explained by the requested-mass change {di:+.2f} mg/ml"
        if cell == "yeast":
            b_ = band.get((c, m))
            if b_ and val is not None:
                b_ = b_ + ([rv] if rv is not None else [])
                if min(b_) - tol <= val <= max(b_) + tol:
                    return "seed", f"within the port seed spread [{min(b_):.6g}, {max(b_):.6g}] ± {tol}"
        return None

    regs = [r for r in table if r["status"] == "regress"]
    exp, unexp = [], []
    for r in regs:
        w = why(r)
        rec = dict(check=r["check"], metric=r["metric"], value=r["value"], reference=r["reference_build_value"], tolerance=r["tolerance"],
                   direction=r["direction"], step=w[0] if w else None, explained=w[1] if w else None)
        (exp if w else unexp).append(rec)
    counts = {}
    for r in table:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    nr = dict(status="PASS" if not unexp else "FAIL", reference_suite=rsuite, counts=counts, regressions_unexpected=unexp,
              regressions_expected=exp, rule="plan §4 tolerances vs the freeze reference; EXPECTED only via criteria §4")
    S2.apply_accepted_regressions(nr, [e for e in acc_extra if e["kind"] == "regression"])
    with open(os.path.join(ctx["out"], "no_regression_freeze.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["check", "metric", "value", "reference", "delta", "direction", "tolerance", "status"])
        for r in table:
            w.writerow([r["check"], r["metric"], r["value"], r["reference_build_value"], r["delta"], r["direction"], r["tolerance"], r["status"]])
    return nr


# ------------------------------------------------------------------------------------------------ main
def verdict(a):
    from wholecell.validate import common
    common.MAX_PROCS = a.nproc
    cell = a.cell
    bdir = os.path.abspath(a.build.rstrip("/"))
    ref = os.path.abspath(a.reference) if a.reference else FREEZE_REF[cell]
    for r_ in [ref, S2.REF[cell]]:
        if r_ and os.path.realpath(bdir) == os.path.realpath(r_):
            raise SystemExit("refusing: --build is a reference build")
    out = a.out or os.path.join(bdir, "validate", "compare_freeze_v8")
    os.makedirs(out, exist_ok=True)
    cfg = jload(os.path.join(bdir, "config.resolved.json"))
    if cfg is None:
        raise SystemExit(f"{bdir}/config.resolved.json missing")
    waive, sat = kv(a.waive), kv(a.satisfied)
    acc_extra = S2.load_accepted(a.accepted, bdir)
    t0 = time.time()
    ctx = dict(bdir=bdir, out=out, nproc=a.nproc, quick=a.quick, complete=bool(a.complete), reuse=False)
    V = dict(step="Freeze (paper build)", cell=cell, build=bdir, reference_build=ref or "PENDING (accepted yeast combined Step-2 build)",
             criteria=CRITERIA, criteria_md5=md5(os.path.join(REPO, CRITERIA)), complete=bool(a.complete), quick=a.quick,
             replicates=f"{a.replicates} ({'decided: 1 seed per cell, D1' if a.replicates == 1 else 'NOT the decided value'})", waived=waive, satisfied=sat,
             accepted_files=[os.path.relpath(os.path.abspath(p), REPO) for p in a.accepted], date=time.strftime("%Y-%m-%d %H:%M"))
    gates = {}

    def run(name, fn):
        t1 = time.time()
        try:
            its = fn()
            g = dict(items=its, status=S2.step_status(its))
        except Exception as e:
            g = dict(status="ERROR", items=[], error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc()[-3000:])
        g["seconds"] = round(time.time() - t1, 1)
        gates[name] = g
        log(name, g["status"], [(i["id"], i["status"]) for i in g["items"]], f"{g['seconds']} s")
        json.dump(dict(V, gates=gates), open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)

    # G2 first: fills ctx suite_new and the C10 tile sums of the build (S2.sec_b1 -> <out>/c10_build.json)
    steps = {}
    if not a.skip_steps:
        steps = g2(cell, bdir, a, ctx, acc_extra)
    else:
        ctx["suite_new"] = S2.suite_rows(os.path.join(bdir, "validate", "suite.json"))
    V["steps"] = steps
    ctx["suite_new"] = ctx.get("suite_new") or S2.suite_rows(os.path.join(bdir, "validate", "suite.json"))
    ctx["suite_freeze_ref"] = S2.suite_rows(os.path.join(ref, "validate", "suite.json")) if ref else {}
    try:
        ctx["c10_build"] = jload(os.path.join(out, "c10_build.json")) or c10_of(cell, bdir, ctx, "build")[0]
        if ref:
            ctx["c10_ref"], V["c10_reference_source"] = c10_of(cell, ref, ctx, "freeze_reference")
    except Exception as e:
        ctx.setdefault("notes", []).append(f"C10 tile sums: {type(e).__name__}: {e}")
    run("G0 integrity", lambda: g0(cell, bdir, ref, cfg, a, ctx))
    gi, open_ = g1_items(cfg, cell, waive, sat)
    run("G1 pending inputs", lambda: gi)
    run("G3/G4 orientation + nucleolus", lambda: g3_g4(cell, bdir, ref, ctx))
    run("G5 yeast B5 overrides", lambda: g5(cell, ctx) or [item("G5", "REPORT", None, "", "yeast only")])
    run("G6 SEC-MS", lambda: g6(cell, bdir, ref, cfg, sat, waive))
    run("G7 surface / PM topology / copy fixes / ledger", lambda: g7(cell, bdir, ref, cfg, sat, waive))
    try:
        V["no_regression"] = g8(cell, bdir, ref, cfg, sat, ctx, acc_extra)
    except Exception as e:
        V["no_regression"] = dict(status="ERROR", error=f"{type(e).__name__}: {e}", traceback=traceback.format_exc()[-2000:])
    gates["G9 replicates"] = dict(status="NA", items=[item("G9 replicates", "REPORT", a.replicates, "criteria §3", "1 seed per cell (D1, 2026-11-19 evening)")])
    V["gates"] = gates
    secs = {**{f"G2 {k}": s for k, s in steps.items()}, **gates}
    V["na_items"] = [f"{n}: {i['id']}: {i['reason']}" for n, s in secs.items() for i in s["items"] if i["status"] == "NA"]
    V["blocking_na"] = [x for x in V["na_items"] if x.startswith(("G0", "G1", "G3", "G4", "G5", "G6", "G7"))] + \
        (["G2: skipped (--skip-steps)"] if a.skip_steps else [])
    V["failing"] = [f"{n}: {i['id']}" for n, s in secs.items() for i in s["items"] if i["status"] == "FAIL"] + \
                   [f"{n}: ERROR" for n, s in secs.items() if s["status"] == "ERROR"] + \
                   ([f"no_regression: {V['no_regression']['status']} ({len(V['no_regression'].get('regressions_unexpected', []))} unexpected; "
                     f"{V['no_regression'].get('reason') or V['no_regression'].get('error') or ''})"]
                    if V["no_regression"]["status"] not in ("PASS", "ACCEPTED_DEVIATION") else [])
    V["accepted_deviations"] = [f"{n}: {i['id']}: {i['reason']}" for n, s in secs.items() for i in s["items"] if i["status"] == "ACCEPTED_DEVIATION"] + \
        [f"no_regression: {r['check']} {r['metric']} = {r['value']} (reference {r['reference']}): {r['explained']}"
         for r in V["no_regression"].get("regressions_accepted", [])]
    V["open_gates"] = open_
    V["notes"] = ctx.get("notes", [])
    V["verdict_rule"] = ("ACCEPTED if no FAIL / ERROR, no blocking NA (open gate or unmeasurable gated item) and no_regression PASS / "
                         "ACCEPTED_DEVIATION; INCOMPLETE if the build / suite is incomplete or the reference is pending")
    if not a.complete or not os.path.exists(os.path.join(bdir, "validate", "suite.json")) or not ref:
        V["verdict"] = "INCOMPLETE"
    else:
        V["verdict"] = "ACCEPTED" if not V["failing"] and not V["blocking_na"] else "NOT ACCEPTED"
    V["seconds"] = round(time.time() - t0, 1)
    json.dump(V, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    nr = V["no_regression"]
    with open(os.path.join(out, "regressions.tsv"), "w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["class", "check", "metric", "value", "reference", "tolerance", "change", "explained"])
        for cls in ("regressions_unexpected", "regressions_accepted", "regressions_expected"):
            for r in nr.get(cls, []):
                w.writerow([cls.split("_")[1], r["check"], r["metric"], r["value"], r["reference"], r["tolerance"], r.get("step"), r.get("explained")])
    log("VERDICT", V["verdict"], "failing", V["failing"], "blocking NA", V["blocking_na"], "->", out)
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", required=True, choices=["jurkat", "yeast"])
    ap.add_argument("--preflight", action="store_true")
    ap.add_argument("--commit")
    ap.add_argument("--main-repo", help="repository whose `main` must contain COMMIT (default: this checkout)")
    ap.add_argument("--criteria-md5")
    ap.add_argument("--build")
    ap.add_argument("--reference", help="freeze reference build (default: FREEZE_REF; yeast pending)")
    ap.add_argument("--complete", type=int, default=1)
    ap.add_argument("--out")
    ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--quick", action="store_true", help="S2 --quick: skip the slow overlap measures in G2")
    ap.add_argument("--skip-steps", action="store_true", help="skip G2 (development only; the verdict is then not a freeze verdict)")
    ap.add_argument("--accepted", action="append", default=[], help="dated accepted-deviations json for THIS build (S2 format)")
    ap.add_argument("--waive", action="append", default=[], help="gate=<dated DECISION>: pending input not included")
    ap.add_argument("--satisfied", action="append", default=[], help="gate=<commit / addendum>: pending input included (no config marker)")
    ap.add_argument("--replicates", type=int, default=1)
    a = ap.parse_args()
    if a.preflight:
        if not a.commit:
            ap.error("--preflight needs --commit")
        sys.exit(preflight(a))
    if not a.build:
        ap.error("--build required")
    if a.skip_steps:
        log("WARNING: --skip-steps: G2 not run (development only)")
    sys.exit(verdict(a))


if __name__ == "__main__":
    main()
