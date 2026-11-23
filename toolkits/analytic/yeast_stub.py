"""Yeast-style analytic primitive generator: INTERFACE STUB (cellgeom_v1 §6.2; future engine `analytic_prims`).

Not implemented: building a yeast package raises NotImplementedError. This module fixes the primitive-table layout the
yeast port will write, so the format, the validator and a future `analytic_prims` engine can be developed against it.

Reference implementation (read-only, NO LICENCE, local only; never write, run or lock anything there):
  <yeast_root>/yeast_cell_dev_bundle/pipeline/cell3d.py   (prims + label_of)
  <yeast_root>/yeast_cell_dev_bundle/scripts/s2_geometry.py (generator)
  <yeast_root>/yeast_cell_dev_bundle/pipeline/cellshape.py (real-shape fields)
Audit: data/derived/repro_audit/yeast/formats.md §3a.

prims.npz (one row per primitive):
  type      int16 (n,)    PRIM_TYPES code
  params    float64 (n,8) PRIM_PARAMS layout per type (yeast PQ)
  comp      int16 (n,)    class -> vocabulary label/membrane via manifest.prims.classes (yeast PC)
  field_bind int16 (n,)   engine field slot the primitive modifies, -1 = standalone
  op        int8 (n,)     OPS code; smooth unions carry k in manifest.prims.smooth_k[comp] (yeast EIS_K, SCAR_K)
  bsphere   float64 (n,4) bounding sphere (yeast bc, br)
Every prim-carrying package also stores baked fields (§3) so sdf_sections_v1 consumers work (validator check:
prim vs baked |Δd| < 0.5 nm within 40 nm of a surface).
"""
import numpy as np

PRIM_TYPES = {"sphere": 0, "capsule": 1, "pill": 2, "cerpatch": 3, "slot": 4, "torus": 5,
              "ellipsoid": 6, "sheet": 7, "lamella": 8, "cisterna_stack": 9}     # 0-5 = yeast T_*; 6-9 reserved
PRIM_PARAMS = {
    "sphere": ["cx", "cy", "cz", "r", "-", "-", "-", "-"],
    "capsule": ["ax", "ay", "az", "bx", "by", "bz", "r", "-"],
    "pill": ["cx", "cy", "cz", "nx", "ny", "nz", "disc_r", "half_thickness"],
    "cerpatch": ["ax", "ay", "az", "cos_half_angle", "half_angle", "mid_r", "half_thickness", "-"],
    "slot": ["ax", "ay", "az", "bx", "by", "bz", "r", "height"],
    "torus": ["cx", "cy", "cz", "nx", "ny", "nz", "R", "r"],
}
OPS = {"union": 0, "subtract": 1, "smooth_union": 2}
YEAST_CLASSES = ["vacuole", "lipid_droplet", "golgi", "peroxisome", "mito_om", "mito_im", "crista", "cortical_er",
                 "pm_invagination", "bud_scar"]                                   # yeast C_* order (PC codes 0-9)


def empty_table(n=0):
    return {"type": np.zeros(n, np.int16), "params": np.zeros((n, 8)), "comp": np.zeros(n, np.int16),
            "field_bind": np.full(n, -1, np.int16), "op": np.zeros(n, np.int8), "bsphere": np.zeros((n, 4))}


def from_yeast_geometry(geometry_json):
    """convert a yeast geometry.json (G, PT, PQ, PC, bc, br, pores, consts) into (prims table, baked fields, manifest
    fragment). Planned mapping: PT -> type, PQ -> params, PC -> comp (YEAST_CLASSES), bc/br -> bsphere; yeast fields
    F_CELL/F_VAC/F_NUC -> cellgeom fields; labels -> core ids (wall_mannan / wall_glucan as extensions >= 1000 with
    parent cell_wall); per-membrane half thicknesses adopted from Pogozheva instead of the fixed HALF = 2."""
    raise NotImplementedError("yeast analytic packages are not built yet (stub; see module docstring)")


def generate(spec, dest=None):
    raise NotImplementedError("yeast analytic packages are not built yet (stub; see module docstring)")
