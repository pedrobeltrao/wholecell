"""Yeast package step (B6, 2026-11-17): per-membrane half thickness. Package N -> N+1, the parent is never edited.

Copies every content file of the parent package (native/, fields/, prims.npz, pores.npz) byte for byte and writes a new
manifest whose engine.params.MEM_HALF = membranes[].half_thickness_nm come from config/geometry/membrane_half_thickness.tsv
(rows cell = yeast; D_P-P / 2 of Pogozheva et al. 2022 SI Table S20). Labels and SDFs do not depend on the half
thickness; the packer clearance does (rule table clearance: thickness_corrected). Criteria Y1-Y5 of
docs/unification_criteria/B6_geometry.md; the validation here = parent checks re-run on the copy + schema + label and
clearance comparison with the parent under the core evaluator.

    python -m toolkits.analytic.yeast_mem_half [parent={geom_root}/yeast/analytic_x5_s0_f5f3a5c723b2] [geometry_id]
"""
import copy
import csv
import datetime as _dt
import hashlib
import json
import os
import shutil
import stat
import sys

import numpy as np

from wholecell.build import git_info
from wholecell.config import REPO, load_site
from wholecell.geometry import cellgeom

TABLE = os.path.join(REPO, "config", "geometry", "membrane_half_thickness.tsv")
SOURCE = "Pogozheva et al. 2022 J Chem Inf Model 62:1036, SI Table S20 (D_P-P / 2); config/geometry/membrane_half_thickness.tsv"
CODE = ("toolkits/analytic/yeast_mem_half.py", "config/geometry/membrane_half_thickness.tsv")


def half_table(cell="yeast"):
    rows = [r for r in csv.DictReader((l for l in open(TABLE) if not l.startswith("#")), delimiter="\t") if r["cell"] == cell]
    return {r["membrane"]: (float(r["half_thickness_nm"]), r["pogozheva_system"] or None, r["reason"]) for r in rows}


def _sha(rel):
    return hashlib.sha256(open(os.path.join(REPO, rel), "rb").read()).hexdigest()


def compare(parent, new, n=2_000_000, seed=0):
    """labels / membrane SDFs / clearance of the new package vs its parent (core evaluator, seeded points: uniform in
    the enclosing box + points within 30 nm of a membrane)"""
    a, b = cellgeom.load(parent, verify="quick"), cellgeom.load(new, verify="quick")
    rng = np.random.default_rng(seed)
    R = a.manifest["domain"]["enclosing_radius_nm"]
    U = rng.uniform(-R, R, (n, 3))
    M = a.mems64(U)
    near = U[(np.abs(M[:, :-1]) < 30.0).any(1)]
    P = np.vstack([U[: n // 2], near[: n // 2]])
    la, lb = a.labels(P), b.labels(P)
    ma, mb = a.mems64(P), b.mems64(P)
    ca, cb = a.clearance(P, -1), b.clearance(P, -1)
    mh_a = np.asarray(a.manifest["engine"]["params"]["MEM_HALF"]); mh_b = np.asarray(b.manifest["engine"]["params"]["MEM_HALF"])
    half = float(b.manifest["engine"]["params"]["HALF"])
    expect = np.min(np.abs(mb) - mh_b[None, :], axis=1) + half     # min_m(|d_m| - MEM_HALF[m]) + HALF, all membranes
    return dict(n=int(len(P)), n_near=int(min(len(near), n // 2)), label_mismatch=int((la != lb).sum()),
                mems_identical=bool(np.array_equal(ma, mb)), clearance_max_abs_delta_nm=float(np.abs(cb - ca).max()),
                clearance_vs_formula_max_err_nm=float(np.abs(cb - expect).max()),
                clearance_old_equals_plain=bool(np.array_equal(ca, np.min(np.abs(ma), axis=1))),
                max_abs_mem_half_change_nm=float(np.abs(mh_b - mh_a).max()))


def main(parent="{geom_root}/yeast/analytic_x5_s0_f5f3a5c723b2", geometry_id=None):
    parent = parent.format(**load_site())
    man0 = json.load(open(os.path.join(parent, "manifest.json")))
    cellgeom.verify_checksums(parent, man0)
    gid = geometry_id or f"{man0['geometry_id']}_memhalf1"
    dest = os.path.join(os.path.dirname(parent), gid)
    if os.path.exists(dest):
        raise FileExistsError(f"{dest} exists (packages are immutable)")
    tmp = f"{dest}.tmp-{os.getpid()}"
    shutil.rmtree(tmp, ignore_errors=True)
    for rel in man0["files"]:
        os.makedirs(os.path.dirname(os.path.join(tmp, rel)), exist_ok=True)
        shutil.copyfile(os.path.join(parent, rel), os.path.join(tmp, rel))
    tab = half_table("yeast")
    man = copy.deepcopy(man0)
    names = [m["name"] for m in sorted(man["membranes"], key=lambda m: m["engine_code"])]
    missing = [n for n in names if n not in tab]
    if missing:
        raise KeyError(f"{TABLE}: no yeast row for {missing}")
    mh = [tab[n][0] for n in names]
    man["engine"]["params"]["MEM_HALF"] = mh
    man["engine"].setdefault("params_source", {})["MEM_HALF"] = SOURCE
    for m in man["membranes"]:
        h, sysname, why = tab[m["name"]]
        m["half_thickness_nm"] = h
        m["half_thickness_source"] = f"Pogozheva 2022 Table S20 {sysname}" if sysname else why
    gi = git_info()
    man.update(geometry_id=gid, package_id=f"{man['cell_id']}/{gid}", created=_dt.datetime.now().isoformat(timespec="seconds"),
               description=man0["description"] + "; per-membrane half thickness (B6)")
    man["generator"] = dict(toolkit="toolkits/analytic/yeast_mem_half.py", repo_commit=gi["commit"], dirty=gi["dirty"],
                            parent=dict(package_id=man0["package_id"], manifest_sha256=cellgeom.sha256_file(os.path.join(parent, "manifest.json")),
                                        generator=man0["generator"]),
                            method="parent content files copied byte for byte; manifest MEM_HALF / half_thickness_nm from the table",
                            code_sha256={f: _sha(f) for f in CODE}, criteria="docs/unification_criteria/B6_geometry.md")
    man["files"] = {rel: dict(sha256=cellgeom.sha256_file(os.path.join(tmp, rel)), bytes=os.path.getsize(os.path.join(tmp, rel)))
                    for rel in sorted(man0["files"])}
    json.dump(man, open(os.path.join(tmp, "manifest.json"), "w"), indent=1, default=float)
    # validation: parent's checks carried (content files identical), schema, comparison with the parent
    v0 = json.load(open(os.path.join(parent, "validation.json")))
    chk = {}
    chk["content_identical_to_parent"] = dict(pass_=man["files"] == man0["files"], n=len(man["files"]))
    errs = [f"{m['name']}: half_thickness_nm != MEM_HALF" for m in man["membranes"] if m["half_thickness_nm"] != mh[m["engine_code"]]]
    chk["schema_half_thickness"] = dict(pass_=not errs, errors=errs, MEM_HALF=dict(zip(names, mh)))
    for k, c in v0["checks"].items():                 # parent checks hold for identical content files
        chk[f"parent.{k}"] = dict(pass_=bool(c.get("pass", True)), carried_from=man0["package_id"])
    cmp_ = compare(parent, tmp)
    cmp_["pass_"] = (cmp_["label_mismatch"] == 0 and cmp_["mems_identical"]
                     and cmp_["clearance_max_abs_delta_nm"] <= cmp_["max_abs_mem_half_change_nm"] + 1e-9
                     and cmp_["clearance_vs_formula_max_err_nm"] <= 1e-9 and cmp_["clearance_old_equals_plain"])
    chk["vs_parent"] = cmp_
    checks = {k: {("pass" if kk == "pass_" else kk): vv for kk, vv in v.items()} for k, v in chk.items()}
    res = dict(spec_version=cellgeom.SPEC_VERSION, package_id=man["package_id"], engine=man["engine"]["name"], checks=checks,
               waived=list(man.get("validation_waivers", {})), manifest_sha256=cellgeom.sha256_file(os.path.join(tmp, "manifest.json")),
               **{"pass": all(c["pass"] for c in checks.values())})
    json.dump(res, open(os.path.join(tmp, "validation.json"), "w"), indent=1, default=float)
    if not res["pass"]:
        raise SystemExit(f"validation failed; left in {tmp}: {[k for k, c in checks.items() if not c['pass']]}")
    os.rename(tmp, dest)
    for root, _, fs in os.walk(dest):
        for f in fs:
            p = os.path.join(root, f)
            os.chmod(p, os.stat(p).st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
    print(json.dumps(dict(package=dest, manifest_sha256=res["manifest_sha256"], MEM_HALF=dict(zip(names, mh)),
                          vs_parent={k: v for k, v in cmp_.items() if k != "pass_"}), indent=1))
    return dest, res


if __name__ == "__main__":
    main(*sys.argv[1:3])
