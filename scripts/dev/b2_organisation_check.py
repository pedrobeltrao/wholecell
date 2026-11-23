# new in wholecell (B2, 2026-11-18)
"""B2 organisation layer on (STRING; AF3 pool off): measurements and verdict
(criteria docs/unification_criteria/B2_organisation.md).

  network   statistics of pack/organisation of a build (or check dir): network.json stats, partner CSR (types with
            partners, partners per type, strength), share of packed mass in types with a profile / with partners
            python scripts/b2_organisation_check.py network --build <bdir> [--out <json>]
  totals    O1: placed (and requested) copies per (species, region) from pack/tiles/*.npz of two builds
            python scripts/b2_organisation_check.py totals --build <bdir> --ref-build <bdir> --out <json>
  verdict   compare_B2/verdict.json of a B2 build vs the B1 build of the same cell: O0 settings (config.resolved.json,
            network.json), O1 totals, O2 C15 (yeast), O3 partner placement, O4 suite no_regression vs the reference's
            validate/suite.json (suite re-run into compare_B2/suite), O5 completeness
            python scripts/b2_organisation_check.py verdict --cell jurkat|yeast --build <bdir> --ref-build <b1 bdir> \
                   [--complete 1] [--ref-status <b1 full_build.status>]
"""
import argparse
import glob
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
AVO = 6.02214076e23
C15_TOL = 0.01                       # §4: localisation / AUROC tolerance
C15_PLAN = {"pca": 0.85, "prox_label": 0.81, "xlms": 0.99}
STRING_MD5 = {9606: "5b7c91d02926e723ad30e81ab9d009a4"}     # frozen input inv_string_links


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def jload(p):
    return json.load(open(p)) if p and os.path.exists(p) else None


def jsave(p, obj):
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p + ".tmp", "w") as f:
        json.dump(obj, f, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    os.replace(p + ".tmp", p)


def _cell_id(cell):
    return "yeast" if cell.startswith("yeast") else "jurkat" if cell.startswith("jurkat") else cell


# ------------------------------------------------------------------ masses (as scripts/b1_benchmark.py A4)
def lib_masses(bdir, cell, L, cache):
    if L.get("mass_da"):
        return np.asarray(L["mass_da"], float), "libs.json mass_da"
    sys.path[:0] = [os.path.join(REPO, "scripts", "dev"), os.path.join(REPO, "scripts")]
    from b1_benchmark import _extra_inputs        # yeast layout inputs (engine/lib, packing trees, geometry)
    from wholecell import config as C
    from wholecell.validate import common
    cfg = C.load(cell)
    common.CACHE_DIR = cache
    m, src = common.lib_masses(common.Inputs(cfg, bdir, **_extra_inputs("yeast", bdir, cfg)), bdir, nproc=3)
    return np.asarray(m, float), src


# ------------------------------------------------------------------ O1 totals
def tile_totals(bdir):
    """{(library id, region name): [placed, requested]} summed over pack/tiles/*.npz, and libs.json"""
    L = jload(os.path.join(bdir, "pack", "libs.json"))
    libs, rn = [str(x) for x in L["libs"]], [str(x) for x in L["rnames"]]
    nl, nr = len(libs), len(rn)
    pla = np.zeros(nl * nr, np.int64); req = np.zeros(nl * nr, np.int64)
    files = glob.glob(os.path.join(bdir, "pack", "tiles", "*.npz"))
    for f in files:
        z = np.load(f)
        k = z["job_lib"].astype(np.int64) * nr + z["job_region"].astype(np.int64)
        np.add.at(pla, k, z["job_placed"].astype(np.int64)); np.add.at(req, k, z["job_n"].astype(np.int64))
    nz = np.nonzero((pla > 0) | (req > 0))[0]
    return {(libs[i // nr], rn[i % nr]): (int(pla[i]), int(req[i])) for i in nz}, L, len(files)


def totals(bdir, ref, cell, cache):
    T, L, nt = tile_totals(bdir)
    R, LR, ntr = tile_totals(ref)
    mass, msrc = lib_masses(ref, cell, LR, cache)
    mlib = {str(l): float(m) for l, m in zip(LR["libs"], mass)}
    mnew, _ = (np.asarray(L["mass_da"], float), "") if L.get("mass_da") else (None, "")
    if mnew is not None:
        for l, m in zip(L["libs"], mnew):
            mlib.setdefault(str(l), float(m))
    keys = sorted(set(T) | set(R))
    rows, out_mass, placed_mass = [], 0.0, 0.0
    n_in = 0
    req_diff = 0
    for k in keys:
        n, nq = R.get(k, (0, 0)); m, mq = T.get(k, (0, 0))
        tol = max(3.0 * np.sqrt(n), 0.01 * n) if n > 0 else 3.0
        d = m - n
        ok = bool(abs(d) <= tol)
        n_in += ok
        req_diff += (mq != nq)
        w = mlib.get(k[0], 0.0)
        placed_mass += n * w
        if not ok:
            out_mass += abs(d) * w
            rows.append(dict(lib=k[0], region=k[1], ref_placed=n, placed=m, delta=d, tol=round(tol, 2),
                             ratio=round(abs(d) / tol, 3), ref_requested=nq, requested=mq, mass_da=w))
    rows.sort(key=lambda r: -r["ratio"])
    frac_in = n_in / max(len(keys), 1)
    out_frac = out_mass / max(placed_mass, 1.0)
    # per region placed / requested mass (report)
    reg = {}
    for k in keys:
        w = mlib.get(k[0], 0.0); a = reg.setdefault(k[1], [0.0, 0.0, 0.0, 0.0])
        a[0] += R.get(k, (0, 0))[0] * w; a[1] += T.get(k, (0, 0))[0] * w; a[2] += R.get(k, (0, 0))[1] * w; a[3] += T.get(k, (0, 0))[1] * w
    regions = {r: dict(ref_placed_pg=a[0] / AVO * 1e12, placed_pg=a[1] / AVO * 1e12, ref_requested_pg=a[2] / AVO * 1e12,
                       requested_pg=a[3] / AVO * 1e12, rel_placed=(a[1] - a[0]) / a[0] if a[0] else None)
               for r, a in sorted(reg.items())}
    return dict(tiles=nt, ref_tiles=ntr, pairs=len(keys), pairs_within=n_in, fraction_within=frac_in,
                out_of_tol_pairs=len(rows), out_of_tol_mass_fraction=out_frac, mass_source=msrc,
                requested_pairs_differing=int(req_diff),
                pass_=bool(frac_in >= 0.995 and out_frac <= 0.001), worst=rows[:100], regions=regions)


# ------------------------------------------------------------------ network statistics (O3, report)
def network(bdir):
    od = os.path.join(bdir, "pack", "organisation")
    J = jload(os.path.join(od, "network.json"))
    if J is None:
        return None
    pr = np.load(os.path.join(od, "partners.npz"), allow_pickle=True)
    nz = np.load(os.path.join(od, "network.npz"), allow_pickle=True)
    ids = [str(x) for x in pr["ids"]]
    npart = np.diff(pr["off"])
    out = dict(stats=J["stats"], string=J.get("string"), inputs=J.get("inputs"),
               partners_per_type=dict(mean=float(npart.mean()), median=float(np.median(npart)),
                                      mean_with_partners=float(npart[npart > 0].mean()) if (npart > 0).any() else 0.0,
                                      at_cap=int((npart >= 16).sum())),
               strength=dict(mean=float(pr["strength"].mean()) if len(pr["strength"]) else None,
                             median=float(np.median(pr["strength"])) if len(pr["strength"]) else None))
    # share of the requested packed copies / mass in types with a profile / with partners (pack/counts.npz if any,
    # else the pack table copies)
    has_prof = np.zeros(len(ids), bool); has_prof[nz["tix"]] = True
    has_part = npart > 0
    L = jload(os.path.join(bdir, "pack", "libs.json"))
    cz = os.path.join(bdir, "pack", "counts.npz")
    if L and os.path.exists(cz):
        c = np.load(cz)
        pos = {u: i for i, u in enumerate(ids)}
        libs = [str(x) for x in L["libs"]]
        to = np.array([pos.get(u, -1) for u in libs])
        n = np.bincount(c["lib"].astype(np.int64), weights=c["n"].astype(float), minlength=len(libs))
        m = n * np.asarray(L.get("mass_da") or np.ones(len(libs)), float)
        sel = to >= 0
        out["requested_share"] = dict(
            copies_profile=float(n[sel][has_prof[to[sel]]].sum() / n.sum()), copies_partners=float(n[sel][has_part[to[sel]]].sum() / n.sum()),
            mass_profile=float(m[sel][has_prof[to[sel]]].sum() / m.sum()), mass_partners=float(m[sel][has_part[to[sel]]].sum() / m.sum()))
    else:
        import pandas as pd
        t = pd.read_csv(os.path.join(bdir, "inventory", "pack_table.tsv"), sep="\t")
        t = t.set_index(t.ing_id.astype(str))
        cop = np.array([float(t.copies.get(u, 0)) for u in ids]); ms = np.array([float(t.mass_da.get(u, 0)) for u in ids]) * cop
        out["pack_table_share"] = dict(copies_profile=float(cop[has_prof].sum() / cop.sum()), copies_partners=float(cop[has_part].sum() / cop.sum()),
                                       mass_profile=float(ms[has_prof].sum() / ms.sum()), mass_partners=float(ms[has_part].sum() / ms.sum()))
    return out


def partner_fraction(bdir):
    S = jload(os.path.join(bdir, "pack", "summary.json")) or {}
    o = S.get("organisation") or {}
    n = S.get("instances") or 0
    return dict(instances=n, placed_next_to_partner=o.get("placed_next_to_partner"),
                copies_with_partner_proposals=o.get("copies_with_partner_proposals"),
                fraction_next_to_partner=(o["placed_next_to_partner"] / n) if (n and o.get("placed_next_to_partner") is not None) else None,
                fraction_with_proposals=(o["copies_with_partner_proposals"] / n) if (n and o.get("copies_with_partner_proposals") is not None) else None)


# ------------------------------------------------------------------ O0 settings
def _flat(d, p=""):
    out = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flat(v, f"{p}.{k}" if p else str(k)))
    else:
        out[p] = json.dumps(d, sort_keys=True, default=str)
    return out


IGNORE = re.compile(r"^(_raw\.|_source|build_id|build_dir|cell$|_cell|_build)")


def settings(bdir, ref, cell):
    a = _flat(jload(os.path.join(ref, "config.resolved.json")) or {})
    b = _flat(jload(os.path.join(bdir, "config.resolved.json")) or {})
    rb, rr = os.path.realpath(bdir), os.path.realpath(ref)
    diff = sorted(k for k in set(a) | set(b) if not IGNORE.match(k)
                  and (a.get(k) or "").replace(rr, "<B>").replace(os.path.basename(rr), "<ID>")
                  != (b.get(k) or "").replace(rb, "<B>").replace(os.path.basename(rb), "<ID>"))
    outside = [k for k in diff if not k.startswith("pack.organisation.")]
    J = jload(os.path.join(bdir, "pack", "organisation", "network.json")) or {}
    st, sg = J.get("stats", {}), J.get("string", {})
    taxon = 9606 if _cell_id(cell) == "jurkat" else 4932
    md5 = ((J.get("inputs") or {}).get("string_links") or {}).get("md5")
    checks = dict(config_only_organisation=not outside, network_json=bool(J), pool_pairs_0=st.get("pool_pairs") == 0,
                  contact_links_0=st.get("contact_links") == 0, taxon=sg.get("taxon") in (taxon, str(taxon)),
                  version=str(sg.get("version")) == "12.0", min_score=int(sg.get("min_score") or 0) == 400,
                  links_md5=(md5 == STRING_MD5[taxon]) if taxon in STRING_MD5 else bool(md5))
    return dict(pass_=all(checks.values()), checks=checks, differing_keys=diff, outside_organisation=outside, links_md5=md5)


# ------------------------------------------------------------------ O2 / O4 suite
def _suite(cell, bdir, ref_suite, out):
    nr = os.path.join(out, "suite", "no_regression.tsv")
    if not os.path.exists(nr):
        subprocess.run([sys.executable, "-m", "wholecell", "validate", "--cell", cell, "--build", bdir, "--step", "suite",
                        "--out", os.path.join(out, "suite"), f"ref_suite={ref_suite}"], cwd=REPO, check=False)
    return nr


def _rows(suite_json):
    J = jload(suite_json) or {}
    return {(r.get("check"), r.get("metric")): r.get("value") for r in (J.get("metrics") or J.get("rows") or [])}


def c15(new_suite, ref_suite):
    A, B = _rows(ref_suite), _rows(new_suite)
    sets = sorted({m.split(".")[0] for (c, m) in A if c == "C15" and m.endswith(".median_auroc.snapshot_assembly")})
    res = {}
    for s in sets:
        k = ("C15", f"{s}.median_auroc.snapshot_assembly")
        ks = ("C15", f"{s}.median_auroc.snapshot")
        a, b = A.get(k), B.get(k)
        ok = a is not None and b is not None and float(b) >= float(a) - C15_TOL
        res[s] = dict(ref=a, new=b, delta=(float(b) - float(a)) if (a is not None and b is not None) else None, pass_=bool(ok),
                      strict_ge=(b is not None and a is not None and float(b) >= float(a)), plan=C15_PLAN.get(s),
                      snapshot_ref=A.get(ks), snapshot_new=B.get(ks))
    return dict(sets=res, pass_=bool(res) and all(v["pass_"] for v in res.values()), tol=C15_TOL)


EXPECTED = ()      # criteria §5: no metric is exempted (the standard §4 tolerances apply)


def _rule_time(bdir, rule):
    import datetime
    f = os.path.join(bdir, "snakemake.log")
    if not os.path.exists(f):
        return None
    ts = start = jid = None; want = False
    for line in open(f):
        m = re.match(r"\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]", line)
        if m:
            ts = datetime.datetime.strptime(re.sub(" +", " ", m.group(1)), "%a %b %d %H:%M:%S %Y"); continue
        if re.match(rf"\s*(local)?rule {rule}:", line):
            start = ts; want = True; continue
        if want and "jobid:" in line:
            jid = line.split("jobid:")[1].strip(); want = False; continue
        if start is not None and jid is not None and re.match(rf"Finished job {jid}\.", line.strip()):
            return (ts - start).total_seconds()
    return None


def verdict(args):
    import pandas as pd
    bdir, ref = args.build, args.ref_build
    out = os.path.join(bdir, "validate", "compare_B2")
    os.makedirs(out, exist_ok=True)
    cid = _cell_id(args.cell)
    items, extra = {}, {}
    ref_suite = os.path.join(ref, "validate", "suite.json")
    extra["reference"] = dict(build=ref, suite=ref_suite, status=args.ref_status,
                              b1_verdict=(jload(os.path.join(ref, "validate", "compare_B1", "verdict.json")) or {}).get("verdict"))
    s = settings(bdir, ref, args.cell); extra["O0"] = s; items["O0"] = s["pass_"]
    log("O0", s["pass_"], s["checks"], "outside:", s["outside_organisation"][:10])
    t = totals(bdir, ref, args.cell, os.path.join(out, "cache")); extra["O1"] = t; items["O1"] = t["pass_"]
    log("O1", t["pass_"], f"{t['pairs_within']}/{t['pairs']} within ({t['fraction_within']:.4f}); out-of-tol mass {t['out_of_tol_mass_fraction']:.2e}")
    nr = _suite(args.cell, bdir, ref_suite, out) if os.path.exists(ref_suite) else None
    new_suite = os.path.join(out, "suite", "suite.json")
    if cid == "yeast":
        c = c15(new_suite, ref_suite); extra["O2"] = c; items["O2"] = c["pass_"]
        log("O2", c["pass_"], {k: (v["ref"], v["new"]) for k, v in c["sets"].items()})
    else:
        extra["O2"] = dict(status="n/a", why="no Jurkat model-neighbour code (C15 suite proximity_na); BioGRID 5.0.263 human "
                           "(BioPlex Huttlin 2015/2017/2021, Go 2021 BioID) and IntAct are on disk for a future adapter")
    pf, pfr = partner_fraction(bdir), partner_fraction(ref)
    extra["O3"] = dict(b2=pf, reference=pfr, network=network(bdir), network_reference=network(ref))
    items["O3"] = bool(pf["fraction_next_to_partner"]) if cid == "jurkat" else True
    log("O3 partner fraction", pf["fraction_next_to_partner"], "ref", pfr["fraction_next_to_partner"])
    regress = []
    if nr and os.path.exists(nr):
        N = pd.read_csv(nr, sep="\t", keep_default_na=False)
        regress = N[N.status == "regress"].to_dict("records")
        regress = [r for r in regress if not any(k in (str(r.get("check", "")) + "." + str(r.get("metric", ""))).lower() for k in EXPECTED)]
    items["O4"] = bool(nr and os.path.exists(nr)) and not regress
    items["O5"] = bool(args.complete)
    extra["wall_s"] = {r: dict(b2=_rule_time(bdir, r), reference=_rule_time(ref, r))
                       for r in ("pack_organisation", "pack_counts", "pack_pack")}
    fail = [k for k, v in items.items() if not v]
    Vd = dict(step="B2", cell=args.cell, criteria="docs/unification_criteria/B2_organisation.md", build=bdir,
              no_regression=nr, items=items, regressions=regress[:50], n_regressions=len(regress), **extra,
              note="C15 gains from STRING-driven placement are partly circular (STRING imports BioGRID / IntAct records)",
              verdict="ACCEPTED" if not fail else "NOT ACCEPTED", failing=fail)
    jsave(os.path.join(out, "verdict.json"), Vd)
    log("verdict", Vd["verdict"], fail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["network", "totals", "verdict"])
    ap.add_argument("--cell", default="jurkat"); ap.add_argument("--build"); ap.add_argument("--ref-build")
    ap.add_argument("--out"); ap.add_argument("--complete", type=int, default=1); ap.add_argument("--ref-status")
    a = ap.parse_args()
    if a.what == "network":
        r = network(a.build)
        print(json.dumps(r, indent=1, default=str))
        if a.out:
            jsave(a.out, r)
    elif a.what == "totals":
        r = totals(a.build, a.ref_build, a.cell, os.path.join(os.path.dirname(a.out or "."), "cache"))
        if a.out:
            jsave(a.out, r)
        print(json.dumps({k: v for k, v in r.items() if k != "worst"}, indent=1, default=str))
    else:
        verdict(a)


if __name__ == "__main__":
    main()
