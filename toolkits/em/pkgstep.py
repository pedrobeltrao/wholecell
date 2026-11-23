"""Package-to-package steps for the EM toolkit: read package N, edit its parts, write package N+1 (never in place).

    cell, parts = pkgstep.open_pkg("{geom_root}/jurkat/em_cell1_s2")    # parts = split G (R, fields, face, sections)
    parts["speckle"] = sec                                              # edit (legacy algorithm on the parts)
    pkgstep.publish(cell, parts, "em_cell1_s3", step="jspeckles", summary=summ, section_meta={...})

publish() carries over everything the step did not change: manifest metadata (frame, labels, membranes, label rule,
engine params unless given), section metadata, pores, meshes/ and layers/ files (layers re-declared with depends_on =
the new package id and carried_from = the parent). publish_new() writes a first package from a G vector (jgeom_build).
"""
import copy
import json
import os

import numpy as np

from wholecell.build import git_info
from wholecell.config import load_site
from wholecell.geometry import cellgeom
from wholecell.geometry.engines import sdf_sections_v1 as e

SECTION_TYPE = {"cristae": "cristae_lamellae_v1", "speckle": "speckle_spheres_v1", "golgi_cis": "golgi_cisternae_v1"}
_CARRY = ("description", "organism", "cell_line", "source", "frame", "domain", "labels", "membranes", "label_rule",
          "validation_waivers")


def _resolve(p):
    return p.format(**load_site())


def open_pkg(src):
    cell = cellgeom.load(_resolve(src))
    return cell, e.split_legacy_G(cell.G, g_n=e.G_N)


def _lossless32(a):
    """float32 when that is lossless (the parent's storage dtype; load() reads float64), so a carried section keeps
    its file bytes (B6 fix: re-saving as float64 changed the golgi_cis sha256 without changing a value)"""
    a = np.asarray(a)
    a32 = a.astype(np.float32)
    return a32 if np.array_equal(a32.astype(np.float64), a.astype(np.float64)) else a


def _write_parts(w, parts, field_names, sec_meta):
    for k, fld in enumerate(parts["fields"]):
        name = field_names[k]
        if int(np.prod(fld[0][4:7])) == 0:
            w.add_empty_field(name, k)
        else:
            w.add_field(name, k, *fld)
    for key, typ in SECTION_TYPE.items():
        if parts.get(key) is None:
            continue
        meta = dict(sec_meta.get(key, {})); meta.pop("type", None); meta.pop("file", None)
        if key == "golgi_cis":
            W, S, (h10, idx, dat) = parts[key]
            arrays = {"W": np.float64(W), "S": np.float64(S), "header": np.r_[h10, 0.0, float(len(idx))],
                      "index": idx.astype(np.int32), "data": _lossless32(dat)}
        else:
            arrays = {"section": parts[key]}
        w.add_section(key, typ, arrays, **meta)


def publish(src_cell, parts, geometry_id, step, summary=None, params=None, pores=None, pore_source=None,
            section_meta=None, generator_extra=None, validate_kw=None, waivers=None):
    man0 = src_cell.manifest
    dest = os.path.join(load_site()["geom_root"], man0["cell_id"], geometry_id)
    w = cellgeom.PackageWriter(dest)
    names = {v: k for k, v in man0["engine"]["field_slots"].items()}
    sec_meta = {k: v for k, v in man0["sections"].items()}
    sec_meta.update(section_meta or {})
    _write_parts(w, parts, [names[k] for k in range(len(parts["fields"]))], sec_meta)
    for k, f in man0["fields"].items():                      # carry field metadata (description, band, provenance)
        if k in w.fields:
            for mk, mv in f.items():
                w.fields[k].setdefault(mk, mv)
    if pores is None and man0.get("pores"):
        z = np.load(os.path.join(src_cell.package_dir, man0["pores"]["file"]))
        pores, pore_source = z["pores"], z["source"]
    if pores is not None:
        pm = {k: v for k, v in (man0.get("pores") or {}).items() if k not in ("file", "n", "columns", "source_codes", "counts")}
        w.set_pores(pores, pore_source, **pm)
    for rel in man0["files"]:
        if rel.startswith(("layers/", "meshes/")):
            w.add_file(rel, os.path.join(src_cell.package_dir, rel))
    for name, L in man0.get("layers", {}).items():
        L = copy.deepcopy(L)
        w.add_layer(name, L.pop("kind"), L.pop("schema_version"), files=L.pop("files"), refs=L.pop("refs"),
                    provenance=L.pop("provenance"), depends_on=None,
                    carried_from=man0["package_id"], **{k: v for k, v in L.items() if k not in ("depends_on", "carried_from")})
    man = {k: copy.deepcopy(man0[k]) for k in _CARRY if k in man0}
    man.update(cell_id=man0["cell_id"], geometry_id=geometry_id)
    man["engine"] = copy.deepcopy(man0["engine"])
    if params:
        man["engine"]["params"].update(params)
    man["domain"]["enclosing_radius_nm"] = parts["R"]
    man["summary"] = summary if summary is not None else copy.deepcopy(man0.get("summary", {}))
    if waivers is not None:
        man["validation_waivers"] = waivers
    gi = git_info()
    man["generator"] = {"toolkit": f"toolkits/em/{step}.py", "repo_commit": gi["commit"], "dirty": gi["dirty"],
                        "parent": {"package_id": man0["package_id"], "manifest_sha256": src_cell.manifest_sha256},
                        **(generator_extra or {})}
    return w.finalize(man, validate_kw=validate_kw)


def publish_new(geometry_id, G, pores, pore_source, centre_label_nm, summary, cell_id="jurkat", step="jgeom_build",
                field_meta=None, generator_extra=None):
    """first package of a chain from a jgeom_build G vector (sdf_sections_v1 layout); metadata as convert_legacy"""
    from toolkits.em import convert_legacy as cl
    parts = e.split_legacy_G(G, g_n=e.G_N)
    dest = os.path.join(load_site()["geom_root"], cell_id, geometry_id)
    w = cellgeom.PackageWriter(dest)
    _write_parts(w, parts, e.FIELD_NAMES[:len(parts["fields"])], {})
    for k, m in (field_meta or cl.FIELD_META).items():
        if k in w.fields:
            w.fields[k].update(m)
    w.set_pores(pores, pore_source, min_sep_nm=150.0, frame="model")
    params = dict(e.DEFAULT_PARAMS)
    mems, _ = cl.membranes_table(params["MEM_HALF"], "Pogozheva et al. 2022")
    c = [float(x) for x in centre_label_nm]
    clip = []
    if parts["face"] is not None:
        of = summary.get("open_face", {})
        clip = [{"name": "open_face", "normal": [0.0, 0.0, 1.0], "point": [0.0, 0.0, parts["face"]],
                 "kept_side": "+normal", "label": "outside_imaged_volume", "cap_method": of.get("method"),
                 "f_vol": of.get("f_vol"), "f_pm": of.get("f_pm")}]
    gi = git_info()
    man = {"cell_id": cell_id, "geometry_id": geometry_id,
           "generator": {"toolkit": f"toolkits/em/{step}.py", "repo_commit": gi["commit"], "dirty": gi["dirty"],
                         **(generator_extra or {})},
           "engine": {"name": e.ENGINE_NAME, "layout_id": e.LAYOUT_ID, "params": params,
                      "field_slots": {n: k for k, n in enumerate(e.FIELD_NAMES)}},
           "frame": {"units": "nm", "handedness": "right", "centre_label_nm": c,
                     "source_affines": {"seg_s0_index_to_model": cl.affine(cl.VOX_S0, -np.asarray(c)),
                                        "pred_s0_index_to_model": cl.affine(cl.VOX_S0, -np.asarray(c), flip_y=cl.Y_S0 - 1),
                                        "raw_s0_index_to_model": cl.affine(cl.VOX_S0, -np.asarray(c) - cl.RAW_OFFSET_S0 * cl.VOX_S0),
                                        "label_frame_nm_to_model": cl.affine([1, 1, 1], -np.asarray(c))}},
           "domain": {"enclosing_radius_nm": parts["R"], "clip_planes": clip},
           "labels": cl.labels_table(), "membranes": mems, "label_rule": cl.label_rule(), "summary": summary}
    return w.finalize(man)


def summary_copy(cell):
    return json.loads(json.dumps(cell.summary))
