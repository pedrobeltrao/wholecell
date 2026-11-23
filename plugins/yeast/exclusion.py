# origin: yeast pipeline/wc3d.py:exclusion_zones @ e3b7448 (constants as parameters; sites from the vendored cell3d)
"""Exclusion-zone provider `yeast.endocytic_sites` (registry kind exclusion_zones; core pack.exclusion[].around.provider).

Ribosome exclusion zones around the endocytic sites (2026-11-15 in the yeast model): ribosomes are absent from a zone
181 +- 41 nm across around the actin network of an endocytic site, reaching about 100 nm beyond the invagination's tip
(Picco et al. 2018 Mol Biol Cell 29:1346; Buser & Drubin 2013 Microsc Microanal 19:381). Modelled as a capsule of
radius r_nm along the invagination axis, from the plasma membrane to beyond_nm - r_nm past the tip, as spheres every
step_nm. Sites (point 17.5 nm inside the PM on the axis, outward direction, tip depth) come from the vendored cell3d
`Cell.endocytic_sites` of the build's geometry package (native/geometry.json; the yeast geometry engine).

    params: {r_nm: 90.5, beyond_nm: 100.0, step_nm: 30.0}  (wc3d.EXCL_R, EXCL_BEYOND, EXCL_STEP)
Returns (m, 4) spheres x, y, z, r in model nm (the yeast frame = the package frame).
"""
import numpy as np

from wholecell.registry import register


def zones(sites, r_nm, beyond_nm, step_nm):
    """capsule spheres along each site axis (origin: wc3d.exclusion_zones)"""
    out = []
    for p0, u, depth in sites:
        p0 = np.asarray(p0, float); u = np.asarray(u, float)
        for s_ in np.arange(0.0, depth + beyond_nm - r_nm + 1e-6, step_nm):
            out.append([*(p0 - u * s_), r_nm])
    return np.array(out, np.float64).reshape(-1, 4)


@register("exclusion_zones", "yeast.endocytic_sites")
def endocytic_sites(cfg, bdir, geom, params):
    from . import geom_engine
    cell = geom_engine.load(geom.source)
    return zones(cell.endocytic_sites(), float(params.get("r_nm", 90.5)), float(params.get("beyond_nm", 100.0)),
                 float(params.get("step_nm", 30.0)))
