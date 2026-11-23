"""Geometry engine `yeast_cell3d_v1`: evaluator of the yeast cellgeom_v1 packages (toolkits/analytic/yeast.py).

The evaluator is the vendored yeast `pipeline/cell3d.py` (numba: sdf_prim, sdf_cer, field_sdf, nuc_sdf, mem_sdfs,
label_from, label_of). A package stores (a) the engine-native files `native/geometry.json` + `native/geometry_fields.npz`
that cell3d.Cell.load reads, and (b) the same content in the shared cellgeom_v1 layout: `fields/{cell,vacuole,nucleus}.npz`
(sparse-brick SDFs, spec §3), `prims.npz` (spec §6.2 primitive table: PT/PQ/PC/bounding spheres), `pores.npz`,
labels / membranes with engine codes. The toolkit validates that both describe the same cell (native round trip).

Runs in the engine environment (config pipeline_engine.python_env), where the vendored pipeline is importable:
    from plugins.yeast import geom_engine; cell = geom_engine.load(package_dir)    # -> cell3d.Cell
"""
import json
import os
import sys

ENGINE_NAME = "yeast_cell3d_v1"
LABELS = ["exterior", "wall_mannan", "wall_glucan", "cytosol", "cER_lumen", "perinuclear", "nucleoplasm", "nucleolus",
          "vacuole", "mito_IMS", "mito_matrix", "lipid_droplet", "golgi_endosome", "peroxisome", "pore_channel"]   # cell3d codes
MEMBRANES = ["PM", "ER", "NE", "VAC", "MITO_OM", "MITO_IM", "GOLGI", "PEROX", "LD", "CELL_OUTER"]
FIELD_SLOTS = {"cell": 0, "vacuole": 1, "nucleus": 2}                     # cell3d.F_CELL, F_VAC, F_NUC


def engine_code_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(here, "engine", "pipeline")


def load(pkg):
    """package dir -> vendored cell3d.Cell (verifies the manifest checksums of the native files first)"""
    import hashlib
    man = json.load(open(os.path.join(pkg, "manifest.json")))
    if man["engine"]["name"] != ENGINE_NAME:
        raise ValueError(f"{pkg}: engine {man['engine']['name']} != {ENGINE_NAME}")
    for rel in ("native/geometry.json", "native/geometry_fields.npz"):
        h = hashlib.sha256(open(os.path.join(pkg, rel), "rb").read()).hexdigest()
        if h != man["files"][rel]["sha256"]:
            raise RuntimeError(f"{pkg}/{rel}: sha256 mismatch")
    if engine_code_dir() not in sys.path:
        sys.path.insert(0, engine_code_dir())
    import cell3d as c3
    return c3.Cell.load(os.path.join(pkg, "native", "geometry.json"))
