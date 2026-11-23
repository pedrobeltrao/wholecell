# new in wholecell (B6, 2026-11-17): no legacy origin; reuses jgeom_build.mesh_sdf / bricks unchanged
"""Recompute the perinuclear midsurface field `ne_mid` of a Jurkat package as the true signed distance to the HMFF NE
midsurface mesh on BOTH sides (package N -> N+1; B6 of the unification plan, criteria docs/unification_criteria/B6_geometry.md).

The v7 defect: `ne_mid` of em_cell1_v7 is a distance only on the nuclear side. Grid nodes outside the mesh hold the
120 nm clamp, so the zero lies ~8-10 nm inside the mesh and the ONM (ne_mid = +PERI_W/2) only a few nm outside it
(perinuclear volume 2.74 um3 instead of ~8.3). This step keeps everything else of the parent package (fields, sections,
pores, mesh, layers: identical files) and rewrites only fields/ne_mid.npz on the SAME grid (origin, bricks, h = 16 nm,
clamps +-120 nm), computed by jgeom_build.mesh_sdf (ray-parity occupancy of the closed mesh for the sign, exact
point-triangle distance within 200 nm, EDT beyond) on the mesh bounding box + MARGIN; nodes outside that box are
outside the mesh by more than the clamp.

usage: NUMBA_NUM_THREADS=3 python -m toolkits.em.jne_fix [src={geom_root}/jurkat/em_cell1_v7] [geometry_id=em_cell1_v8]
"""
import copy
import os
import sys
import time

import numpy as np

CODE = ("toolkits/em/jne_fix.py", "toolkits/em/jgeom_build.py", "toolkits/em/cellshape.py")
MARGIN = 320.0          # nm around the mesh bounding box (> clamp 120 + 2 brick widths of slack)


def _sha(rel):
    import hashlib
    from wholecell.config import REPO
    return hashlib.sha256(open(os.path.join(REPO, rel), "rb").read()).hexdigest()


def ne_field(cell, mesh_path, h10):
    """new ne_mid bricks (h10, index, data) on the grid of header h10, from the mesh in the label frame"""
    import trimesh
    from toolkits.em import jgeom_build as jb
    lo_f = np.asarray(h10[:3], np.float64); h = float(h10[3]); nb = np.asarray(h10[4:7], int); B = int(h10[7])
    cin, cout = float(h10[8]), float(h10[9])
    assert h == jb.H, f"grid h {h} != jgeom_build.H {jb.H}"
    m = trimesh.load(mesh_path, process=False)
    V = np.asarray(m.vertices, np.float64) - np.asarray(cell.manifest["frame"]["centre_label_nm"], np.float64)
    F = np.asarray(m.faces, np.int64)
    n = nb * B + 1                                                     # nodes of the full grid
    i0 = np.maximum(np.floor((V.min(0) - MARGIN - lo_f) / h).astype(int), 0)
    i1 = np.minimum(np.ceil((V.max(0) + MARGIN - lo_f) / h).astype(int) + 1, n)
    lo_s = lo_f + i0 * h
    t = time.time()
    sub = jb.mesh_sdf(V, F, lo_s, tuple(i1 - i0))                      # toolkit function, unchanged
    print(f"mesh_sdf on sub-grid {tuple(i1 - i0)} ({time.time() - t:.0f} s)", flush=True)
    full = np.full(tuple(n), cout, np.float32)                         # outside the box: outside the mesh by > clamp
    full[i0[0]:i1[0], i0[1]:i1[1], i0[2]:i1[2]] = sub
    vol_inside = float((full < 0).sum() * h ** 3 / 1e9)
    hh, idx, dat = jb.bricks(full, lo_f, (cin, cout))
    assert np.array_equal(hh[:10], np.asarray(h10[:10], np.float64)), "brick header changed"
    return (hh, idx, dat), vol_inside, float(m.area / 1e6), float(m.volume / 1e9)


def main(src="{geom_root}/jurkat/em_cell1_v7", geometry_id="em_cell1_v8"):
    from toolkits.em import pkgstep
    t0 = time.time()
    cell, parts = pkgstep.open_pkg(src)
    man = cell.manifest
    k = man["engine"]["field_slots"]["ne_mid"]
    mesh_rel = man["layers"]["ne_mesh"]["files"][0]
    new, vol_in, area, mvol = ne_field(cell, os.path.join(cell.package_dir, mesh_rel), parts["fields"][k][0])
    parts["fields"][k] = new
    summ = pkgstep.summary_copy(cell)
    summ["ne_mid"] = dict(vol_inside_um3=round(vol_in, 1), area_um2=round(area, 1), mesh_volume_um3=round(mvol, 2),
                          method="jne_fix: exact signed distance to the HMFF mesh, both sides (B6)")
    cell.manifest = man = copy.deepcopy(man)                           # carried into the new package by pkgstep.publish
    man["description"] = ("Jurkat E6-1 cell 1 of OpenOrganelle jrc_jurkat-1: the pack-v7 geometry with the perinuclear "
                          "midsurface field ne_mid recomputed as a true signed distance on both sides (B6)")
    man["fields"]["ne_mid"]["description"] = "perinuclear midsurface: exact signed distance to the HMFF NE mesh, both sides"
    man["fields"]["ne_mid"]["recomputed_by"] = "toolkits/em/jne_fix.py (B6; v7 field had the clamp on the cytoplasmic side)"
    waivers = {"field_volumes_areas": "summary of the unchanged fields er / endo / golgi is the legacy jgeom summary "
                                      "(mask volumes, marching-cubes areas; golgi = EM blob before the v7 ribbon); MC "
                                      "values are reported, not compared. ne_mid passes against its mesh."}
    man_new, res = pkgstep.publish(cell, parts, geometry_id, step="jne_fix", summary=summ, waivers=waivers,
                                   generator_extra={"method": "ne_mid recomputed (jgeom_build.mesh_sdf on the mesh box "
                                                              f"+ {MARGIN:.0f} nm); all other files carried",
                                                    "criteria": "docs/unification_criteria/B6_geometry.md",
                                                    "code_sha256": {f: _sha(f) for f in CODE}})
    print(f"{man_new['package_id']}: validation pass {res['pass']}, failed {res.get('failed_checks')}, "
          f"waived {res.get('waived_failures')} ({time.time() - t0:.0f} s)", flush=True)
    return man_new, res


if __name__ == "__main__":
    main(*sys.argv[1:3])
