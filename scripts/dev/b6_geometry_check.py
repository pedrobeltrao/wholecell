"""B6 geometry fixes: checks of docs/unification_criteria/B6_geometry.md (frozen criteria, md5 in PLANNING_NOTES).

    PYTHONHASHSEED=0 python scripts/b6_geometry_check.py jurkat [--out DIR]   pre-build: J1-J8, C1-C4 (packages v7 vs v8,
                                                                               frozen pre-placed layers)
    PYTHONHASHSEED=0 python scripts/b6_geometry_check.py yeast [--out DIR]    pre-build: Y1-Y5 (old vs new package)
    PYTHONHASHSEED=0 python scripts/b6_geometry_check.py verdict --cell C --build BDIR [--complete 0|1] [--seed-suites S1 S2 ...]
        build level: B1-B4 + the pre-build results + suite no_regression with the EXPECTED list
        -> BDIR/validate/compare_B6/{verdict.json, regressions.tsv}
Pre-build outputs: {results_root}/<cell>/_dev_unify/b6/geometry_check.json (default --out). Nothing is written into
packages or reference builds. Threads: NUMBA_NUM_THREADS (default here 3).
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import subprocess
import sys
import time

os.environ.setdefault("NUMBA_NUM_THREADS", "3")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell import config as wcfg  # noqa: E402
from wholecell.geometry import cellgeom  # noqa: E402

SITE = wcfg.load_site()
GEOM, RES = SITE["geom_root"], SITE["results_root"]
V7, V8 = f"{GEOM}/jurkat/em_cell1_v7", f"{GEOM}/jurkat/em_cell1_v8"
Y_OLD, Y_NEW = f"{GEOM}/yeast/analytic_x5_s0_f5f3a5c723b2", f"{GEOM}/yeast/analytic_x5_s0_f5f3a5c723b2_memhalf1"
REF = {"jurkat": f"{RES}/jurkat/20261115-ec0b403-s20261115-v7repro2", "yeast": f"{RES}/yeast/20261116-aabe7b5-s0-x5repro"}
NUC_INSIDE = ("nucleoplasm", "nucleolus", "nuclear_speck")
CRITERIA = "docs/unification_criteria/B6_geometry.md"


def jdump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + f".tmp{os.getpid()}"
    json.dump(obj, open(tmp, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    os.replace(tmp, path)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def names(cell, lab):
    return np.array([cell.code_to_name.get(int(c), f"?{c}") for c in range(64)], dtype=object)[lab]


def transitions(cell_a, cell_b, P, chunk=5_000_000):
    """labels of points P under two packages -> (names_a, names_b) (chunked)"""
    la, lb = [], []
    for i in range(0, len(P), chunk):
        Q = np.ascontiguousarray(P[i:i + chunk], np.float64)
        la.append(cell_a.labels(Q)); lb.append(cell_b.labels(Q))
    la, lb = np.concatenate(la) if la else np.zeros(0, np.int8), np.concatenate(lb) if lb else np.zeros(0, np.int8)
    return names(cell_a, la), names(cell_b, lb)


def ttable(na, nb):
    out = {}
    for a, b in zip(*np.unique(np.c_[na, nb].astype(str), axis=0, return_counts=True)):
        out[f"{a[0]} -> {a[1]}"] = int(b)
    return out


# ------------------------------------------------------------------ Jurkat
def mesh(cell):
    import trimesh
    m = trimesh.load(os.path.join(cell.package_dir, cell.manifest["layers"]["ne_mesh"]["files"][0]), process=False)
    V = np.asarray(m.vertices, np.float64) - np.asarray(cell.manifest["frame"]["centre_label_nm"], np.float64)
    return V, np.asarray(m.faces, np.int64), m


def mesh_dist(V, F, P, k=32):
    from scipy.spatial import cKDTree
    from toolkits.em import cellshape as cs
    _, cand = cKDTree(V[F].mean(1)).query(P, k=k, workers=3)
    return cs._min_tri_dist(np.ascontiguousarray(P, np.float64), cand.astype(np.int64), V, F)


def check_jurkat(out):
    t0 = time.time()
    a, b = cellgeom.load(V7, verify="quick"), cellgeom.load(V8, verify="quick")
    ma, mb = a.manifest, b.manifest
    va, vb = json.load(open(f"{V7}/validation.json")), json.load(open(f"{V8}/validation.json"))
    R = {"criteria": CRITERIA, "packages": {"v7": ma["package_id"], "v8": mb["package_id"],
                                            "v8_manifest_sha256": cellgeom.sha256_file(f"{V8}/manifest.json")}}
    # J1 only fields/ne_mid.npz differs; params equal; v7 untouched
    diff = sorted(k for k in set(ma["files"]) | set(mb["files"]) if ma["files"].get(k, {}).get("sha256") != mb["files"].get(k, {}).get("sha256"))
    v7_sha = cellgeom.sha256_file(f"{V7}/manifest.json")
    R["J1"] = dict(files_differing=diff, params_equal=ma["engine"]["params"] == mb["engine"]["params"], v7_manifest_sha256=v7_sha,
                   pass_=diff == ["fields/ne_mid.npz"] and ma["engine"]["params"] == mb["engine"]["params"]
                   and v7_sha == "eb16248b7505f17f135c8c822bcba6e8dc02cff8ad67b7723eeef6e5419762c7")
    # J2 validator
    C = vb["checks"]
    e = C["eikonal"]["fields"]["ne_mid"]; fa = C["field_volumes_areas"]["fields"]["ne_mid"]
    R["J2"] = dict(pass_validator=vb["pass"], failed=vb["failed_checks"], waived=vb["waived_failures"], pores=C["pores"],
                   eikonal_ne_mid=e, eikonal_ne_mid_v7=va["checks"]["eikonal"]["fields"]["ne_mid"], field_ne_mid=fa,
                   label_rule=C["label_rule"]["agreement"], layers_pass=C["layers"]["pass"],
                   pass_=bool(vb["pass"] and C["pores"]["pass"] and "pores" not in mb.get("validation_waivers", {})
                              and 0.95 <= e["median_grad"] <= 1.05 and e["p95"] <= 1.10 and fa["pass"]
                              and vb["waived_failures"] == ["field_volumes_areas"] and C["label_rule"]["agreement"] == 1.0
                              and C["layers"]["pass"]))
    # J3 field = true distance (200,000 points within 100 nm of the mesh)
    V, F, m = mesh(b)
    rng = np.random.default_rng(20261117)
    area = np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1)
    tri = rng.choice(len(F), 200_000, p=area / area.sum())
    r1, r2 = rng.random(200_000), rng.random(200_000); s = np.sqrt(r1)
    X = (1 - s)[:, None] * V[F[tri, 0]] + (s * (1 - r2))[:, None] * V[F[tri, 1]] + (s * r2)[:, None] * V[F[tri, 2]]
    nrm = np.cross(V[F[tri, 1]] - V[F[tri, 0]], V[F[tri, 2]] - V[F[tri, 0]]); nrm /= np.linalg.norm(nrm, axis=1, keepdims=True)
    if np.mean(np.einsum("ij,ij->i", nrm, X - V.mean(0))) < 0:
        nrm = -nrm
    off = rng.uniform(-100, 100, 200_000)
    P = X + off[:, None] * nrm
    d = mesh_dist(V, F, P)
    fb = b.engine.field_points(P, b.G, mb["fields"]["ne_mid"]["slot"])
    fa_ = a.engine.field_points(P, a.G, ma["fields"]["ne_mid"]["slot"])
    big = np.abs(off) > 20
    pz = b.engine.field_points(b.pores[:, :3], b.G, mb["fields"]["ne_mid"]["slot"])
    R["J3"] = dict(n=200_000, abs_err_p50=float(np.median(np.abs(np.abs(fb) - d))), abs_err_p99=float(np.percentile(np.abs(np.abs(fb) - d), 99)),
                   abs_err_max=float(np.abs(np.abs(fb) - d).max()), sign_agree_offset_gt20=float((np.sign(fb[big]) == np.sign(off[big])).mean()),
                   v7_abs_err_p99=float(np.percentile(np.abs(np.abs(fa_) - d), 99)), pore_centre_abs_median=float(np.median(np.abs(pz))),
                   pore_centre_abs_max=float(np.abs(pz).max()))
    R["J3"]["pass_"] = R["J3"]["abs_err_p99"] <= 1.0 and R["J3"]["pore_centre_abs_median"] <= 0.5
    log("J3", {k: v for k, v in R["J3"].items()})
    # J4 / J5 region volumes (validator MC: 2,000,000 points, seed 0, stream geometry.mc_volumes)
    A, B = va["checks"]["volumes_mc"], vb["checks"]["volumes_mc"]
    vol = {k: dict(v7=A["volume_um3"][k], v8=B["volume_um3"][k], delta=round(B["volume_um3"][k] - A["volume_um3"][k], 4),
                   ci95=B["ci95_um3"][k]) for k in A["volume_um3"]}
    R["J4"] = dict(perinuclear_space=vol["perinuclear_space"], expected=8.08, accept=[7.7, 8.5],
                   pass_=7.7 <= vol["perinuclear_space"]["v8"] <= 8.5)
    ok = {}
    ok["cytosol"] = -7.0 <= vol["cytosol"]["delta"] <= -4.0
    ok["nucleoplasm"] = abs(vol["nucleoplasm"]["delta"]) <= 1.0
    ok["nuclear_pore_channel"] = 0.30 <= vol["nuclear_pore_channel"]["v8"] <= 0.60
    for k in ("er_lumen", "mito_ims", "mito_matrix", "golgi_lumen", "golgi_intercisternal", "endolysosome_lumen"):
        x = vol[k]; ok[k] = (x["delta"] <= 2 * x["ci95"]) and (x["delta"] >= -(0.02 * x["v7"] + 2 * x["ci95"]))
    for k in ("nucleolus", "nuclear_speck", "exterior", "outside_imaged_volume"):
        ok[k] = abs(vol[k]["delta"]) <= 2 * vol[k]["ci95"]
    R["J5"] = dict(volumes=vol, ok=ok, pass_=all(ok.values()))
    # J6 change confined to the NE (2,000,000 uniform points)
    lo, hi = cellgeom._box(b)
    U = lo + (hi - lo) * np.random.default_rng(6).random((2_000_000, 3))
    la, lb = transitions(a, b, U)
    ch = la != lb
    dch = mesh_dist(V, F, U[ch]) if ch.any() else np.zeros(0)
    R["J6"] = dict(n=2_000_000, changed=int(ch.sum()), changed_max_dist_to_mesh_nm=float(dch.max()) if len(dch) else 0.0,
                   transitions=ttable(la[ch], lb[ch]), pass_=bool(len(dch) == 0 or dch.max() <= 150.0))
    log("J6", R["J6"]["changed"], R["J6"]["changed_max_dist_to_mesh_nm"])
    # J7 areas, J8 pores
    fva, fvb = va["checks"]["field_volumes_areas"]["fields"], vb["checks"]["field_volumes_areas"]["fields"]
    R["J7"] = dict(ne_mid_area_shell_um2=fvb["ne_mid"]["area_um2_shell"], ne_mid_area_mesh_um2=round(m.area / 1e6, 2),
                   other_fields_identical={k: ma["files"][f"fields/{k}.npz"] == mb["files"][f"fields/{k}.npz"] for k in fva if k != "ne_mid"},
                   other_fields_mc_v8={k: fvb[k] for k in fvb if k != "ne_mid"})
    R["J7"]["pass_"] = abs(R["J7"]["ne_mid_area_shell_um2"] / R["J7"]["ne_mid_area_mesh_um2"] - 1) <= 0.05 and all(R["J7"]["other_fields_identical"].values())
    pk = ("PORE_R", "PERI_W", "PORE_ZWIN", "PORE_LW")
    R["J8"] = dict(pores_sha_equal=ma["files"]["pores.npz"] == mb["files"]["pores.npz"], counts=mb["pores"]["counts"],
                   params={k: mb["engine"]["params"][k] for k in pk})
    R["J8"]["pass_"] = R["J8"]["pores_sha_equal"] and all(ma["engine"]["params"][k] == mb["engine"]["params"][k] for k in pk)
    # C1-C3 frozen pre-placed layers
    cfg = wcfg.load("jurkat")
    fi = cfg["frozen_inputs"]
    mods = {}
    for mod in cfg["preplaced"]["modules"]:
        key = (mod.get("frozen") or {}).get("npz") if mod.get("builder") else mod["npz"]
        mods[mod["name"]] = (mod.get("kind", "instances"), fi[key]["path"], fi[key].get("md5"))
    lay = {}
    for name, (kind, path, md5) in mods.items():
        z = np.load(path)
        P = z["pos"] if kind == "instances" else z["pts"]
        la, lb = transitions(a, b, P)
        ent = dict(kind=kind, file=path, md5_declared=md5, n=int(len(P)), changed=int((la != lb).sum()),
                   transitions=ttable(la[la != lb], lb[la != lb]),
                   in_perinuclear_v7=int((la == "perinuclear_space").sum()), in_perinuclear_v8=int((lb == "perinuclear_space").sum()))
        ent["perinuclear_v8_frac"] = ent["in_perinuclear_v8"] / max(ent["n"], 1)
        if kind == "instances":
            T = z["types"]; ids = [str(x) for x in z["lib_ids"]]
            ent["perinuclear_v8_by_lib"] = {ids[t]: int(c) for t, c in zip(*np.unique(T[lb == "perinuclear_space"], return_counts=True))}
        if name == "chromatin":
            nuc = T == ids.index("NUCLEOSOME")
            ins_b, ins_a = np.isin(lb, NUC_INSIDE), np.isin(la, NUC_INSIDE)
            ent["C1"] = dict(n_nucleosomes=int(nuc.sum()), inside_frac_v8=float(ins_b[nuc].mean()), inside_frac_v7=float(ins_a[nuc].mean()),
                             outside_v8=int((~ins_b[nuc]).sum()), all_types_inside_frac_v8=float(ins_b.mean()),
                             outside_labels_v8={k: int(v) for k, v in zip(*np.unique(lb[nuc & ~ins_b].astype(str), return_counts=True))},
                             threshold=0.999)
            ent["C1"]["pass_"] = ent["C1"]["inside_frac_v8"] >= 0.999
            ent["C2_trim_needed"] = not ent["C1"]["pass_"]
        lay[name] = ent
        log("layer", name, ent["n"], "changed", ent["changed"], "perinuclear v8", ent["in_perinuclear_v8"])
    R["C1"] = lay["chromatin"]["C1"]
    R["C2"] = dict(trim_applied=lay["chromatin"]["C2_trim_needed"], note="trim only if C1 fails (pre-stated rule)")
    R["C3"] = dict(layers=lay, index_drop="ne_lumen_b6 (inside_labels perinuclear_space, all pre-placed libraries)",
                   expected_index_drop_instances=int(sum(v["in_perinuclear_v8"] for v in lay.values() if v["kind"] == "instances")),
                   pass_=all(v["perinuclear_v8_frac"] <= 0.01 for v in lay.values() if v["kind"] == "instances"))
    # C4 depends_on_geometry (config check without a build: the package of each config)
    from wholecell.preplaced.stage import check_depends_on_geometry
    c4 = {}
    for cname in ("jurkat", "jurkat_v7_port"):
        c = wcfg.load(cname)
        rows, bad = check_depends_on_geometry(c, "/nonexistent")
        c4[cname] = dict(package=rows[0]["package_id"] if rows else None, mismatches=bad, n_modules=len(rows))
    R["C4"] = dict(configs=c4, pass_=all(not v["mismatches"] for v in c4.values()) and c4["jurkat"]["package"] == "jurkat/em_cell1_v8"
                   and c4["jurkat_v7_port"]["package"] == "jurkat/em_cell1_v7")
    R["seconds"] = round(time.time() - t0, 1)
    R["pass"] = all(R[k]["pass_"] for k in ("J1", "J2", "J3", "J4", "J5", "J6", "J7", "J8", "C1", "C3", "C4"))
    jdump(R, os.path.join(out, "geometry_check.json"))
    log("jurkat pre-build checks:", {k: R[k]["pass_"] for k in R if isinstance(R[k], dict) and "pass_" in R[k]}, "->", out)
    return R


# ------------------------------------------------------------------ yeast
def check_yeast(out):
    t0 = time.time()
    from wholecell import seeds
    mo, mn = json.load(open(f"{Y_OLD}/manifest.json")), json.load(open(f"{Y_NEW}/manifest.json"))
    vn = json.load(open(f"{Y_NEW}/validation.json"))
    sys.path.insert(0, REPO)
    from toolkits.analytic.yeast_mem_half import half_table
    tab = half_table("yeast")
    names_ = [m["name"] for m in sorted(mn["membranes"], key=lambda m: m["engine_code"])]
    R = {"criteria": CRITERIA, "packages": {"old": mo["package_id"], "new": mn["package_id"],
                                            "new_manifest_sha256": cellgeom.sha256_file(f"{Y_NEW}/manifest.json")}}
    R["Y1"] = dict(MEM_HALF=dict(zip(names_, mn["engine"]["params"]["MEM_HALF"])), table={k: v[0] for k, v in tab.items()},
                   HALF=mn["engine"]["params"]["HALF"])
    R["Y1"]["pass_"] = all(R["Y1"]["MEM_HALF"][k] == tab[k][0] for k in names_) and R["Y1"]["HALF"] == 2.0 and \
        all(m["half_thickness_nm"] == mn["engine"]["params"]["MEM_HALF"][m["engine_code"]] for m in mn["membranes"])
    old_sha = cellgeom.sha256_file(f"{Y_OLD}/manifest.json")
    R["Y2"] = dict(content_identical=mn["files"] == mo["files"], parent_manifest_sha256=old_sha,
                   pass_=mn["files"] == mo["files"] and old_sha == "8ee3c1bc813f0051de8d08955f38e7b6866aef98eb04317315f1c159b3675302")
    # Y3: old package compiles to the plain rule; G identical to the A9 rule table (main 8f6be0d); new: clearance check
    old_rules = os.path.join(out, "yeast_cell3d_v1.a9.yaml")
    os.makedirs(out, exist_ok=True)
    open(old_rules, "w").write(subprocess.run(["git", "-C", REPO, "show", "8f6be0d:config/geometry/rules/yeast_cell3d_v1.yaml"],
                                              capture_output=True, text=True, check=True).stdout)
    c_old_a9 = cellgeom.load(Y_OLD, verify="quick", rules=old_rules)
    c_old = cellgeom.load(Y_OLD, verify="quick")
    c_new = cellgeom.load(Y_NEW, verify="quick")
    from wholecell.geometry.engines import cellgeom_eval_v1 as ev
    mode = lambda c: int(c.G[int(c.G[ev.G_TABLES]) + 4])
    R["Y3"] = dict(old_G_identical_to_A9=bool(np.array_equal(c_old.G, c_old_a9.G)), old_mode=["thickness", "plain"][mode(c_old)],
                   new_mode=["thickness", "plain"][mode(c_new)], vs_parent=vn["checks"]["vs_parent"])
    R["Y3"]["pass_"] = R["Y3"]["old_G_identical_to_A9"] and R["Y3"]["old_mode"] == "plain" and R["Y3"]["new_mode"] == "thickness" \
        and vn["checks"]["vs_parent"]["clearance_max_abs_delta_nm"] <= 0.11 + 1e-9 and vn["checks"]["vs_parent"]["pass"]
    R["Y4"] = dict(validation_pass=vn["pass"], checks={k: v["pass"] for k, v in vn["checks"].items()}, pass_=vn["pass"])
    # Y5: labels (validation: 1.1 M points) + MC volumes (2,000,000, the measures stream of yeast.yaml)
    seed = seeds.seed_int(wcfg.load("yeast"), "geometry.measures")
    vo = cellgeom.mc_volumes(c_old, n=2_000_000, seed=seed, stream="geometry.measures")
    vw = cellgeom.mc_volumes(c_new, n=2_000_000, seed=seed, stream="geometry.measures")
    R["Y5"] = dict(label_mismatch=vn["checks"]["vs_parent"]["label_mismatch"], mems_identical=vn["checks"]["vs_parent"]["mems_identical"],
                   mc_volumes_identical=vo["volume_um3"] == vw["volume_um3"], volumes_um3=vw["volume_um3"], seed=seed,
                   max_boundary_shift_nm=vn["checks"]["vs_parent"]["max_abs_mem_half_change_nm"])
    R["Y5"]["pass_"] = R["Y5"]["label_mismatch"] == 0 and R["Y5"]["mems_identical"] and R["Y5"]["mc_volumes_identical"] \
        and R["Y5"]["max_boundary_shift_nm"] <= 2.0
    R["measures_reference"] = vo                     # for the build verdict (B3: identical volumes)
    R["seconds"] = round(time.time() - t0, 1)
    R["pass"] = all(R[k]["pass_"] for k in ("Y1", "Y2", "Y3", "Y4", "Y5"))
    jdump(R, os.path.join(out, "geometry_check.json"))
    log("yeast pre-build checks:", {k: R[k]["pass_"] for k in R if isinstance(R[k], dict) and "pass_" in R[k]}, "->", out)
    return R


# ------------------------------------------------------------------ build verdict
EXPECTED_JURKAT = {   # (check, metric) -> rule on (value, reference) explained in the criteria, Part B
    ("C11", "density.cytosol"): "up <= +3.0 mg/ml and mass-consistent (density x cytosol volume within 1 %)",
    ("C11", "density.perinuclear"): "report (no species maps to the NE lumen)",
    ("C11", "nuc_cyto_ratio"): "report",
    ("C11", "phi.cytosol"): "report",
    ("C11", "density.ER_lumen"): "mass-consistent with a <= 2 % volume loss",
    ("C11", "density.matrix"): "mass-consistent with a <= 2 % volume loss",
    ("C11", "density.golgi_lumen"): "mass-consistent with a <= 2 % volume loss",
    ("C11", "density.endo_lumen"): "mass-consistent with a <= 2 % volume loss",
    ("C16", "depends_on_geometry_mismatch"): "new row (n/a in the reference): must be 0",
}


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def verdict(cell, bdir, complete, seed_suites, outname="compare_B6", accepted=None):
    out = os.path.join(bdir, "validate", outname)
    os.makedirs(out, exist_ok=True)
    acc = {(a["check"], a["metric"]): a for a in (json.load(open(accepted))["accepted"] if accepted else [])}
    pre = json.load(open(f"{RES}/{cell}/_dev_unify/b6/geometry_check.json"))
    V = {"cell": cell, "build": bdir, "criteria": CRITERIA, "reference_build": REF[cell], "pre_build_pass": pre["pass"]}
    exp_pkg = pre["packages"]["v8" if cell == "jurkat" else "new"]
    # B1
    pj, gv = (json.load(open(os.path.join(bdir, "geometry", f))) if os.path.exists(os.path.join(bdir, "geometry", f)) else {}
              for f in ("package.json", "validation.json"))
    V["B1"] = dict(build_complete=complete, package_id=pj.get("package_id"), geometry_validation_pass=gv.get("pass"),
                   pass_=bool(complete and pj.get("package_id") == exp_pkg and gv.get("pass")))
    # B3 measures vs reference (Jurkat: J4/J5 bounds on the build's own MC; yeast: identical to the old package)
    mp = os.path.join(bdir, "geometry", "measures.json")
    meas = json.load(open(mp))["region_volumes_mc"] if os.path.exists(mp) else None
    if cell == "jurkat":
        refm = json.load(open(os.path.join(REF[cell], "geometry", "measures.json")))["region_volumes_mc"]
        if meas:
            v, r, ci = meas["volume_um3"], refm["volume_um3"], meas["ci95_um3"]
            d = {k: v[k] - r[k] for k in v}
            ok = dict(perinuclear=7.7 <= v["perinuclear_space"] <= 8.5, cytosol=-7.0 <= d["cytosol"] <= -4.0,
                      nucleoplasm=abs(d["nucleoplasm"]) <= 1.0, pore=0.30 <= v["nuclear_pore_channel"] <= 0.60)
            for k in ("er_lumen", "mito_ims", "mito_matrix", "golgi_lumen", "golgi_intercisternal", "endolysosome_lumen"):
                ok[k] = d[k] <= 2 * ci[k] and d[k] >= -(0.02 * r[k] + 2 * ci[k])
            for k in ("nucleolus", "nuclear_speck", "exterior", "outside_imaged_volume"):
                ok[k] = abs(d[k]) <= 2 * ci[k]
            V["B3"] = dict(volumes=v, reference=r, delta={k: round(x, 4) for k, x in d.items()}, ok=ok, pass_=all(ok.values()))
        else:
            V["B3"] = dict(pass_=False, why=f"missing {mp}")
    else:
        same = meas is not None and meas["volume_um3"] == pre["measures_reference"]["volume_um3"]
        V["B3"] = dict(identical_to_old_package=same, pass_=bool(same))
    # B4 index drop + chromatin
    if cell == "jurkat":
        lg = os.path.join(bdir, "logs", "pack_index.log")
        m = re.search(r"index drop 'ne_lumen_b6': ([\d,]+)", open(lg).read()) if os.path.exists(lg) else None
        n = int(m.group(1).replace(",", "")) if m else None
        V["B4"] = dict(index_drop_count=n, expected=pre["C3"]["expected_index_drop_instances"], chromatin_C1=pre["C1"],
                       # frozen B4 text: "count <= 1 % per module". The pre-build check C3 showed every module <= 1 % at the expected
                       # drop; the build may drop fewer (instances already removed by earlier rules), never more (fix 2026-11-17: was ==)
                       pass_=bool(n is not None and n <= pre["C3"]["expected_index_drop_instances"] and pre["C3"]["pass_"] and pre["C1"]["pass_"]))
    else:
        V["B4"] = dict(pass_=True, note="yeast: no frozen pre-placed layers (serial builders in the build)")
    # B2 suite no_regression
    nr = os.path.join(bdir, "validate", "no_regression.tsv")
    rows = list(csv.DictReader(open(nr), delimiter="\t")) if os.path.exists(nr) else []
    suite = json.load(open(os.path.join(bdir, "validate", "suite.json")))["metrics"] if os.path.exists(os.path.join(bdir, "validate", "suite.json")) else []
    srow = {(r["check"], r["metric"]): r for r in suite}
    band = {}
    for sp in seed_suites or []:                          # yeast: seed spread of the existing port builds
        for r in json.load(open(sp))["metrics"]:
            x = _f(r.get("value"))
            if x is not None:
                band.setdefault((r["check"], r["metric"]), []).append(x)
    regs, expl = [], []
    for r in rows:
        if r.get("status") != "regress":
            continue
        key = (r["check"], r["metric"]); val, ref, tol = _f(r["value"]), _f(r["reference_build_value"]), _f(r.get("tolerance")) or 0.0
        why = None
        if cell == "jurkat" and key in EXPECTED_JURKAT:
            rule = EXPECTED_JURKAT[key]
            if key == ("C11", "density.cytosol") and meas and val is not None and ref is not None:
                vr = refm["volume_um3"]["cytosol"] / meas["volume_um3"]["cytosol"]
                if val - ref <= 3.0 and abs(val - ref * vr) <= 0.01 * ref * vr:
                    why = f"expected: {rule} (reference x V_ref/V_new = {ref * vr:.2f})"
            elif rule.startswith("report"):
                why = f"expected: {rule}"
            elif rule.startswith("mass-consistent") and val is not None and ref is not None and val <= ref * 1.02 + tol:
                why = f"expected: {rule}"
        if cell == "yeast" and key in band and val is not None:
            b_ = band[key] + ([ref] if ref is not None else [])
            if min(b_) - tol <= val <= max(b_) + tol:
                why = f"within the port seed spread s0/s1/s2 [{min(b_):.6g}, {max(b_):.6g}] +- {tol}"
        if why is None and key in acc:                                  # dated DECISION (accepted deviation), see --accepted file
            why = f"accepted: {acc[key]['decision']} ({acc[key]['reason']})"
        (expl if why else regs).append(dict(check=key[0], metric=key[1], value=val, reference=ref, tolerance=tol, explained=why))
    c16 = srow.get(("C16", "depends_on_geometry_mismatch"), {})
    V["B2"] = dict(no_regression=nr if rows else None, regressions_unexplained=regs, regressions_expected=expl,
                   c16_depends_on_geometry_mismatch=c16.get("value"), c13_prepass_md5_ok=srow.get(("C13", "prepass_md5_ok"), {}).get("value"),
                   pass_=bool(rows) and not regs and _f(c16.get("value")) == 0.0)
    V["accepted"] = bool(pre["pass"] and all(V[k]["pass_"] for k in ("B1", "B2", "B3", "B4")))
    V["verdict"] = "ACCEPTED" if V["accepted"] else "NOT ACCEPTED"
    V["failed"] = [k for k in ("B1", "B2", "B3", "B4") if not V[k]["pass_"]] + ([] if pre["pass"] else ["pre_build"])
    jdump(V, os.path.join(out, "verdict.json"))
    with open(os.path.join(out, "regressions.tsv"), "w") as f:
        f.write("check\tmetric\tvalue\treference\ttolerance\texplained\n")
        for r in regs + expl:
            f.write("\t".join(str(r[k]) for k in ("check", "metric", "value", "reference", "tolerance", "explained")) + "\n")
    log(cell, V["verdict"], "failed", V["failed"], "->", out)
    return V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["jurkat", "yeast", "verdict"])
    ap.add_argument("--out"); ap.add_argument("--cell"); ap.add_argument("--build")
    ap.add_argument("--complete", type=int, default=1); ap.add_argument("--seed-suites", nargs="*", default=[])
    ap.add_argument("--outname", default="compare_B6"); ap.add_argument("--accepted", default=None, help="json of dated accepted deviations")
    a = ap.parse_args()
    if a.what == "verdict":
        verdict(a.cell, a.build, bool(a.complete), a.seed_suites, a.outname, a.accepted)
    else:
        out = a.out or f"{RES}/{a.what}/_dev_unify/b6"
        (check_jurkat if a.what == "jurkat" else check_yeast)(out)


if __name__ == "__main__":
    main()
