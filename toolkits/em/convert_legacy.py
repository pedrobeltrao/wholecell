"""Wrap the published Jurkat v7 geometry into a cellgeom_v1 package (no re-segmentation, no recomputation).

    PYTHONHASHSEED=0 NUMBA_CACHE_DIR=<tmp> python -m toolkits.em.convert_legacy [--dest {geom_root}/jurkat/em_cell1_v7]
        [--no-b1] [--threads 8]

Inputs (read-only, {data_root}/data/derived): jgeom_cell1_v7.npz (G, pores, centre_nm, R, summary),
cell1_NE_hmff_midsurface.ply, cytoskel/stage3/mt_stage3_b_centerlines.npz, prepass_v7/cytoskeleton_cell1.json,
cytoskel/microvilli_summary.json; referenced by path+md5 only: ribo_seg_instances.parquet, np_seg_instances.parquet,
jgeom_cell1_v7_emgolgi.npz. The fields, sections and pores are split out of the legacy G vector unchanged
(sdf_sections_v1.split_legacy_G; field data stored float32 where lossless, which holds for all v7 fields), the legacy
summary is kept verbatim as manifest.summary, and criterion B1 is run on the package before it is published.
"""
import argparse
import json
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
sys.dont_write_bytecode = True

import numpy as np  # noqa: E402

from wholecell.build import git_info  # noqa: E402
from wholecell.config import load_site  # noqa: E402
from wholecell.geometry import cellgeom  # noqa: E402
from wholecell.geometry.engines import sdf_sections_v1 as e  # noqa: E402

GEOMETRY_ID = "em_cell1_v7"
VOX_S0 = np.array([4.0, 4.0, 3.44])          # jrc_jurkat-1 s0 voxel (nm), x y z
Y_S0 = 3000                                   # s0 extent along y (for the *_pred y-flip)
RAW_OFFSET_S0 = np.array([0, 9, -2])          # raw = label + offset (s0 voxels), jurkat/em/io.py
N_PORES_EM = 300                              # first 300 pores = np_seg detections (jgeom_build log), rest random (rng 7)

FIELD_META = {
    "cell": dict(description="plasma-membrane midsurface (refined cell-1 mask, extruded below the open face)",
                 band_nm=[160.0, 64.0], smoothing_sigma_nm=16.0, source_labels="cell1_labels_s2 ids 1-6"),
    "ne_mid": dict(description="perinuclear midsurface: exact distance to the HMFF NE mesh", band_nm=[120.0, 120.0],
                   source="cell1_NE_hmff_midsurface.ply"),
    "mito": dict(description="mitochondrial outer membrane (EM mito_seg)", band_nm=[64.0, 64.0], smoothing_sigma_nm=12.0),
    "er": dict(description="ER membrane (EM er_seg + er_pred, thresholded)", band_nm=[48.0, 48.0], smoothing_sigma_nm=10.0),
    "golgi": dict(description="Golgi envelope: juxta-centrosomal ribbon blob (jgolgi_ribbon.py, v7); the EM Golgi blob "
                  "is kept as reference layer golgi_em_blob", band_nm=[64.0, 64.0]),
    "endo": dict(description="endosome + lysosome + vesicle lumen (EM)", band_nm=[64.0, 64.0], smoothing_sigma_nm=12.0),
    "nucleolus": dict(description="nucleolus (EM nucleolus_seg inside the nucleus), no membrane", band_nm=[64.0, 64.0],
                      smoothing_sigma_nm=16.0),
}

# engine code -> core label name (sdf_sections_v1 label codes)
ENGINE_LABELS = {e.EXT: "exterior", e.CYTO: "cytosol", e.PERI: "perinuclear_space", e.NUCP: "nucleoplasm",
                 e.NUCL: "nucleolus", e.IMS: "mito_ims", e.MATRIX: "mito_matrix", e.ERL: "er_lumen",
                 e.GOLGI: "golgi_lumen", e.ENDO: "endolysosome_lumen", e.PORE: "nuclear_pore_channel",
                 e.OUTV: "outside_imaged_volume", e.SPK: "nuclear_speck", e.GIC: "golgi_intercisternal"}


def labels_table(engine_labels=ENGINE_LABELS):
    core = {L["name"]: L for L in cellgeom.core_labels_json()}
    out = []
    for code, name in sorted(engine_labels.items()):
        L = dict(core[name]); L["engine_code"] = int(code); L["engine_name"] = e.LABEL_NAMES[code]
        out.append(L)
    return out


def membranes_table(mem_half, half_source):
    """sdf_sections_v1 membranes. sign: engine membrane SDF negative on p_side; gradient points to n_side."""
    mh = list(mem_half)
    return [
        dict(name="PM", engine_code=e.M_PM, field="cell", derivation="-cell", offset_nm=0.0, half_thickness_nm=mh[0],
             n_side="cytosol", p_side="exterior", leaflets="bilayer", go="GO:0005886"),
        dict(name="NE", engine_code=e.M_NE, field="ne_mid",
             derivation="|ne_mid| - PERI_W/2; within pores: torus of radius PORE_R + PERI_W/2 around the pore axis",
             offset_nm=["+PERI_W/2 (ONM)", "-PERI_W/2 (INM)"], half_thickness_nm=mh[1],
             n_side=["cytosol (ONM)", "nucleoplasm (INM)"], p_side="perinuclear_space", leaflets="bilayer",
             go="GO:0005635", surfaces=[{"name": "ONM", "go": "GO:0005640"}, {"name": "INM", "go": "GO:0005637"},
                                        {"name": "pore_membrane", "go": "GO:0044611"}]),
        dict(name="ER", engine_code=e.M_ER, field="er", derivation="er", offset_nm=0.0, half_thickness_nm=mh[2],
             n_side="cytosol", p_side="er_lumen", leaflets="bilayer", go="GO:0005789"),
        dict(name="MITO_OM", engine_code=e.M_OM, field="mito", derivation="mito", offset_nm=0.0,
             half_thickness_nm=mh[3], n_side="cytosol", p_side="mito_ims", leaflets="bilayer", go="GO:0005741"),
        dict(name="MITO_IM", engine_code=e.M_IM, field="mito",
             derivation="-(max(mito + IMS_W, -crista_lumen)) (boundary IM + cristae section)", offset_nm="-IMS_W",
             half_thickness_nm=mh[4], n_side="mito_matrix", p_side="mito_ims", leaflets="bilayer", go="GO:0005743"),
        dict(name="GOLGI", engine_code=e.M_GOLGI, field="golgi",
             derivation="golgi blob, or max(blob, cisternal slab) with section golgi_cis", offset_nm=0.0,
             half_thickness_nm=mh[5], n_side="cytosol", p_side="golgi_lumen", leaflets="bilayer", go="GO:0000139"),
        dict(name="ENDO", engine_code=e.M_ENDO, field="endo", derivation="endo", offset_nm=0.0,
             half_thickness_nm=mh[6], n_side="cytosol", p_side="endolysosome_lumen", leaflets="bilayer",
             go="GO:0010008"),
        dict(name="CELL_OUTER", engine_code=e.M_OUTER, field="cell", derivation="min(-cell, z - z_face) (domain wall)",
             offset_nm=0.0, half_thickness_nm=mh[7], n_side="cytosol", p_side="exterior", leaflets="none (hard wall)",
             go=None),
    ], half_source


def label_rule():
    """declarative mirror of sdf_sections_v1.label_from (validated against the engine)"""
    return [
        {"if": [{"clip": "open_face"}], "label": "outside_imaged_volume"},
        {"if": [{"mem": "CELL_OUTER", "lt": 0}], "label": "exterior"},
        {"if": [{"mem": "NE", "lt": 0}], "label": "perinuclear_space"},
        {"if": [{"field": "ne_mid", "lt": "-0.5*PERI_W"}], "then": [
            {"if": [{"field": "ne_mid", "gt": "-PERI_W-PORE_LW"}, {"in_pore": "PORE_R+0.5*PERI_W", "zwin": "PORE_ZWIN"}],
             "label": "nuclear_pore_channel"},
            {"if": [{"field": "nucleolus", "lt": 0}], "label": "nucleolus"},
            {"if": [{"in_section": "speckle"}], "label": "nuclear_speck"},
            {"label": "nucleoplasm"}]},
        {"if": [{"field": "ne_mid", "lt": 0}], "label": "nuclear_pore_channel"},
        {"if": [{"mem": "MITO_OM", "lt": 0}], "then": [
            {"if": [{"var": "matrix", "lt": 0}], "label": "mito_matrix"},
            {"label": "mito_ims"}]},
        {"if": [{"mem": "ER", "lt": 0}], "label": "er_lumen"},
        {"if": [{"mem": "GOLGI", "lt": 0}], "label": "golgi_lumen"},
        {"if": [{"mem": "ENDO", "lt": 0}], "label": "endolysosome_lumen"},
        {"if": [{"has_section": "golgi_cis"}, {"field": "golgi", "lt": 0}], "label": "golgi_intercisternal"},
        {"label": "cytosol"},
    ]


def affine(scale, translate, flip_y=None):
    A = np.eye(4); A[:3, :3] = np.diag(scale); A[:3, 3] = translate
    if flip_y is not None:
        A[1, 1] = -scale[1]; A[1, 3] = scale[1] * flip_y + translate[1]
    return A.tolist()


def ref(path, role):
    return {"path": path, "md5": cellgeom.md5_file(path), "bytes": os.path.getsize(path), "role": role}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--dest", default="{geom_root}/jurkat/" + GEOMETRY_ID)
    ap.add_argument("--legacy", default="{data_root}/data/derived/jgeom_cell1_v7.npz")
    ap.add_argument("--no-b1", action="store_true")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args(argv)
    import numba
    numba.set_num_threads(a.threads)
    site = load_site()
    D = os.path.join(site["data_root"], "data", "derived")
    dest = a.dest.format(**site); lpath = a.legacy.format(**site)
    z = np.load(lpath)
    G = z["G"]; pores = z["pores"]; centre = z["centre_nm"]; summ = json.loads(str(z["summary"]))
    parts = e.split_legacy_G(G)
    del G
    w = cellgeom.PackageWriter(dest)
    for k, (name, fld) in enumerate(zip(e.FIELD_NAMES, parts["fields"])):
        w.add_field(name, k, *fld, **FIELD_META[name])
    sp = summ.get("cristae", {})
    w.add_section("cristae", "cristae_lamellae_v1", {"section": parts["cristae"]},
                  params={k: sp.get(k) for k in ("n_mito", "with_cristae", "S_range")},
                  description="procedural lamellar cristae (jcristae.py, seed 17): n, id-grid origin xyz, spacing, "
                              "dims xyz, n x 12 lamella params, dense id grid (relative offsets)",
                  engine_params=["CR_W", "CR_GAP", "CR_JW", "IMS_W"])
    ss = summ.get("speckles", {})
    w.add_section("speckle", "speckle_spheres_v1", {"section": parts["speckle"]},
                  params={k: ss.get(k) for k in ("n", "r_nm", "seed", "d_edge_nm", "sep_nm", "source")},
                  description="nuclear speckles: n, n x (cx, cy, cz, r) nm (jspeckles.py, seed 5, v8 chromatin)",
                  label="nuclear_speck")
    W, S, (h10, idx, dat) = parts["golgi_cis"]
    w.add_section("golgi_cis", "golgi_cisternae_v1",
                  {"W": np.float64(W), "S": np.float64(S), "header": np.r_[h10, 0.0, float(len(idx))],
                   "index": idx.astype(np.int32), "data": dat.astype(np.float32)},
                  params={"W_nm": W, "S_nm": S}, engine_params=["CIS_CLR"],
                  description="cisternal Golgi (jgolgi.py seed 11 on the jgolgi_ribbon blob): lumen = layers parallel to "
                              "the blob surface from a deep SDF (inside clamp 320 nm)", label="golgi_intercisternal")
    src = np.zeros(len(pores), np.int8); src[:N_PORES_EM] = 1
    w.set_pores(pores, src, min_sep_nm=150.0, frame="model",
                provenance="jgeom_build.py: 300 np_seg centroids projected on the HMFF NE mesh (triangle centroid + "
                           "outward normal), 700 random area-weighted triangle centroids (rng 7), >= 150 nm apart")
    mesh_rel = w.add_file("meshes/ne_mid.ply", os.path.join(D, "cell1_NE_hmff_midsurface.ply"))
    w.add_layer("ne_mesh", "mesh_midsurface", 1, files=[mesh_rel],
                provenance={"source": "data/derived/cell1_NE_hmff_midsurface.ply",
                            "method": "HMFF fit (FreeDTS, scratch_checks/hmff_pilot_nucleus.py + hmff_export.py), not "
                                      "bitwise reproducible -> deposited"},
                frame="label_frame_nm", to_model="subtract frame.centre_label_nm", membrane="NE (midsurface)")
    mt_rel = w.add_file("layers/mt_traces/mt_stage3_b_centerlines.npz",
                        os.path.join(D, "cytoskel/stage3/mt_stage3_b_centerlines.npz"))
    w.add_layer("mt_traces", "cytoskeleton_tracks", 1, files=[mt_rel],
                provenance={"source": "data/derived/cytoskel/stage3/mt_stage3_b_centerlines.npz",
                            "method": "U-Net MT tracing + stage-3 variant b (jurkat/em/cytoskel/mt_stage3.py)",
                            "built_on": "jgeom_cell1.npz (same model frame / centre_nm as v7)"},
                frame="model", type="microtubule", keys="points_model_nm, offsets, flipped_for_lattice, ...")
    cj = json.load(open(os.path.join(D, "prepass_v7/cytoskeleton_cell1.json")))
    cen = cj["log"]["centrioles"]
    cen_rel = w.add_file("layers/centrosome/centrosome.json", text=json.dumps(
        {"centre_nm": cen["new_centre_nm"], "mtoc_candidate_nm": cen["mtoc_candidate_nm"],
         "old_centre_nm": cen["old_centre_nm"], "extent_radius_nm": cen["extent_radius_nm"],
         "frame": "model", "key_used_by_jvirtual": "log.centrioles.new_centre_nm"}, indent=1))
    w.add_layer("centrosome", "centrosome", 1, files=[cen_rel],
                refs=[ref(os.path.join(D, "prepass_v7/cytoskeleton_cell1.json"), "source json (v7 cytoskeleton pre-pass)")],
                provenance={"method": "centriole pair placed by jcytoskel_v7.py near the EM MTOC candidate"})
    w.add_layer("npc_positions", "npc", 1, files=["pores.npz"],
                provenance={"em_detections": "np_seg instances (scratch_checks/count_instances_per_cell.py)"},
                refs=[ref(os.path.join(D, "np_seg_instances.parquet"), "EM np_seg instance centroids (s0 voxels)")],
                counts={"em": N_PORES_EM, "random": int(len(pores) - N_PORES_EM)})
    w.add_layer("ribosome_priors", "point_priors", 1,
                refs=[ref(os.path.join(D, "ribo_seg_instances.parquet"),
                          "EM ribo_seg instance centroids, s0 voxels label frame, columns cell/cell2")],
                provenance={"method": "scratch_checks/count_instances_per_cell.py + jurkat/em/assign_instances.py"})
    mv_rel = w.add_file("layers/microvilli/microvilli_summary.json", os.path.join(D, "cytoskel/microvilli_summary.json"))
    w.add_layer("microvilli", "microvilli", 1, files=[mv_rel],
                provenance={"method": "jurkat/em/cytoskel/microvilli.py + microvilli_summary.py (8 nm PM patches)"})
    w.add_layer("golgi_em_blob", "reference_geometry", 1,
                refs=[ref(os.path.join(D, "jgeom_cell1_v7_emgolgi.npz"), "pre-ribbon geometry with the EM Golgi blob")],
                provenance={"note": "v7 replaced the EM Golgi blob by the juxta-centrosomal ribbon (jgolgi_ribbon.py)"})

    params = dict(e.DEFAULT_PARAMS)
    mems, _ = membranes_table(params["MEM_HALF"], "Pogozheva et al. 2022 SI Table S20 D_P-P/2 (data/lipids_jurkat.json)")
    c = [float(x) for x in centre]
    gi = git_info()
    man = {
        "cell_id": "jurkat", "geometry_id": GEOMETRY_ID,
        "description": "Jurkat E6-1 cell 1 of OpenOrganelle jrc_jurkat-1: the published pack-v7 geometry, wrapped",
        "organism": {"name": "Homo sapiens", "taxon": 9606}, "cell_line": "Jurkat E6-1",
        "source": {"dataset": "OpenOrganelle jrc_jurkat-1", "licence": "CC BY 4.0", "doi": "10.25378/janelia.13114211",
                   "voxel_s0_nm": VOX_S0.tolist(), "cell": 1,
                   "legacy_file": {"path": "data/derived/jgeom_cell1_v7.npz", "md5": cellgeom.md5_file(lpath)}},
        "generator": {"toolkit": "toolkits/em/convert_legacy.py", "repo_commit": gi["commit"], "dirty": gi["dirty"],
                      "method": "wrap (no recomputation)"},
        "engine": {"name": e.ENGINE_NAME, "layout_id": e.LAYOUT_ID, "params": params,
                   "params_source": {"MEM_HALF": "Pogozheva et al. 2022 SI Table S20 (D_P-P / 2)",
                                     "others": "cell3d.py v7 constants"},
                   "field_slots": {n: k for k, n in enumerate(e.FIELD_NAMES)}},
        "frame": {"units": "nm", "handedness": "right", "axes": "x, y, z of the OpenOrganelle label (_seg) frame",
                  "origin": "bounding-box centre of the refined cell-1 mask (isotropic 16 nm label grid)",
                  "centre_label_nm": c,
                  "label_frame_convention": "label-frame nm = s2 index x (16, 16, 13.76) = s0 index x (4, 4, 3.44); "
                                            "OME physical voxel centre of s2 = label-frame nm + (6, 6, 5.16)",
                  "source_affines": {
                      "seg_s0_index_to_model": affine(VOX_S0, -np.asarray(c)),
                      "pred_s0_index_to_model": affine(VOX_S0, -np.asarray(c), flip_y=Y_S0 - 1),
                      "raw_s0_index_to_model": affine(VOX_S0, -np.asarray(c) - RAW_OFFSET_S0 * VOX_S0),
                      "label_frame_nm_to_model": affine([1, 1, 1], -np.asarray(c))},
                  "notes": "*_pred layers are y-flipped vs *_seg (y_seg = 2999 - y_pred at s0); raw EM = label + (0, 9, -2) s0 voxels"},
        "domain": {"enclosing_radius_nm": parts["R"],
                   "clip_planes": [{"name": "open_face", "normal": [0.0, 0.0, 1.0], "point": [0.0, 0.0, parts["face"]],
                                    "kept_side": "+normal", "label": "outside_imaged_volume",
                                    "reason": "cell 1 truncated by the EM volume z = 0 face; fields extruded below it",
                                    "cap_method": summ["open_face"]["method"],
                                    "f_vol": summ["open_face"]["f_vol"], "f_pm": summ["open_face"]["f_pm"]}]},
        "labels": labels_table(),
        "membranes": mems,
        "label_rule": label_rule(),
        "summary": summ,
        "validation_waivers": {
            "pores": "legacy v7 defect kept for exact reproduction (B1): the ne_mid field is a distance only on its "
                     "nuclear side; outside the zero crossing it rises with slope ~5-7 to the 120 nm clamp, and the "
                     "zero lies ~8-10 nm inside the HMFF mesh, so pore centres (mesh triangle centroids) sit at "
                     "ne_mid ~ +60 nm and the ONM only a few nm outside the zero (perinuclear volume ~2.7 um3 instead "
                     "of ~8.3 um3). Fix = rebuild ne_mid (new geometry version).",
            "field_volumes_areas": "summary is the legacy jgeom summary verbatim, measured by the legacy generator on "
                                   "other quantities (mask volumes; marching-cubes and mesh areas; golgi = EM blob "
                                   "before the v7 ribbon; ne_mid area from the mesh, field defect above); MC values "
                                   "are reported, not compared"},
        "summary_note": "legacy jgeom summary, verbatim (packer reads open_face.f_vol/f_pm and <field>.area_um2); "
                        "golgi.* is the EM-blob value used by v7, golgi_ribbon.* the ribbon",
    }

    def b1(cell):
        if a.no_b1:
            return {}
        from toolkits.em import b1_check
        return {"b1_legacy_equality": b1_check.run(cell.G, cell.pores, site["data_root"], legacy=a.legacy,
                                                   threads=a.threads, source="cellgeom.load(package, tmp dir)")}
    man, res = w.finalize(man, extra_checks=b1)
    print(json.dumps({"dest": dest, "pass": res["pass"], "manifest_sha256": res["manifest_sha256"],
                      "checks": {k: v.get("pass") for k, v in res["checks"].items()}}, indent=1))
    return 0 if res["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
