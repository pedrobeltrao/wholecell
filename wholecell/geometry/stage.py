"""Geometry stage: python -m wholecell geometry --cell C --build B [revalidate=1] [mc_points=N]

Loads the cellgeom_v1 package named in config `geometry.package` (checksums verified) with the evaluator
`geometry.evaluator` (cellgeom_eval_v1) and rule table `geometry.rules`, checks the manifest sha256 against
`geometry.manifest_sha256` and the package layout against `geometry.engine`, reuses the package's validation.json when
it is fresh (same manifest sha256, passed) or re-runs the validator, and writes (write-once) into <bdir>/geometry/:
  package.json     package_dir / path, package_id, manifest_sha256, engine (layout), spec_version, evaluator,
                   rules (path) + rules_sha256
  validation.json  the validator result (+ reused flag, source)
  measures.json    region volumes by seeded Monte Carlo (stream geometry.measures), membrane areas and per-field
                   volumes from the package summary, open-face factors
"""
from __future__ import annotations

import json
import os

from wholecell import build, seeds
from wholecell.geometry import cellgeom


def _write_json(path, obj):
    build.write_once(path)
    build.atomic_save(path, lambda p: json.dump(obj, open(p, "w"), indent=1, default=float))


def measures(cell, cfg, n=2_000_000):
    seed = seeds.seed_int(cfg, "geometry.measures")
    vol = cellgeom.mc_volumes(cell, n=n, seed=seed, stream="geometry.measures")
    summ = cell.summary or {}
    fields = {k: v for k, v in summ.items() if isinstance(v, dict) and ("area_um2" in v or "vol_sdf_um3" in v
                                                                         or "vol_inside_um3" in v)}
    mem_area = {}
    for m in cell.manifest["membranes"]:
        f = summ.get(m.get("field"), {}) if isinstance(summ.get(m.get("field")), dict) else {}
        if "area_um2" in f and m["name"] != "CELL_OUTER":
            note = ""
            if m["name"] == "NE":
                note = "midsurface area of the field (NE: one midsurface; ONM/INM each ~ this)"
            elif m["name"] == "MITO_IM":
                note = "field (OM) area: boundary IM only; cristae (section) not included" if "cristae" in cell.manifest["sections"] else "field (OM) area"
            elif m["name"] == "GOLGI" and "golgi_cis" in cell.manifest["sections"]:
                note = "summary blob area; cisternal membranes (section golgi_cis) not included"
            mem_area[m["name"]] = {"area_um2": f["area_um2"], "from_field": m["field"], "note": note}
    of = summ.get("open_face", {}) if isinstance(summ.get("open_face"), dict) else {}
    return {"package_id": cell.package_id, "manifest_sha256": cell.manifest_sha256,
            "region_volumes_mc": vol, "membrane_areas_summary": mem_area, "field_summary": fields,
            "open_face": {"f_vol": of.get("f_vol", 1.0), "f_pm": of.get("f_pm", 1.0)},
            "seed_stream": "geometry.measures", "seed": seed}


def run(cfg, bdir, step=None, **kw):
    g = cfg.get("geometry") or {}
    pkg = g.get("package")
    if not pkg:
        raise KeyError("config geometry.package is not set")
    out = os.path.join(bdir, "geometry")
    os.makedirs(out, exist_ok=True)
    cell = cellgeom.load_cfg(cfg, verify=True)
    man = cell.manifest
    exp = g.get("manifest_sha256")
    if exp and exp != cell.manifest_sha256:
        raise RuntimeError(f"{pkg}: manifest sha256 {cell.manifest_sha256} != config geometry.manifest_sha256 {exp}")
    if g.get("engine") and g["engine"] != man["engine"]["name"]:
        raise RuntimeError(f"config geometry.engine {g['engine']} != package engine {man['engine']['name']}")
    pj = {"package_dir": cell.package_dir, "path": cell.package_dir, "package_id": cell.package_id,
          "manifest_sha256": cell.manifest_sha256, "engine": man["engine"]["name"], "spec_version": man["spec_version"],
          "expected_sha256_ok": (exp == cell.manifest_sha256) if exp else None,
          "evaluator": cell.evaluator, "rules": cell.rules_path, "rules_sha256": cell.rules_sha256}
    vpath = os.path.join(cell.package_dir, "validation.json")
    val = json.load(open(vpath)) if os.path.exists(vpath) else None
    fresh = (val is not None and val.get("manifest_sha256") == cell.manifest_sha256 and val.get("pass")
             and not int(kw.get("revalidate", 0)))
    if fresh:
        val = dict(val, reused=True, source=vpath)
    else:
        val = cellgeom.validate(cell, seed=seeds.seed_int(cfg, "geometry.validate"))
        val.update(reused=False, source="wholecell.geometry.cellgeom.validate", manifest_sha256=cell.manifest_sha256)
    meas = measures(cell, cfg, n=int(kw.get("mc_points", g.get("mc_points", 2_000_000))))
    _write_json(os.path.join(out, "package.json"), pj)
    _write_json(os.path.join(out, "validation.json"), val)
    _write_json(os.path.join(out, "measures.json"), meas)
    if not val["pass"]:
        raise RuntimeError(f"geometry validation failed: {val.get('failed_checks')}")
    return pj
