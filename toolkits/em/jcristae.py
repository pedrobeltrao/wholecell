# origin: jurkat/yc/pipeline/jcristae.py (legacy workspace, 2026-11-15)
# restructured: reads package N, writes package N+1 (toolkits/em/pkgstep.py); paths from config; NOT RUN in wholecell
"""Procedural lamellar cristae for the mitochondria of the Jurkat cell geometry (cell3d cristae section of G).

The EM segmentation gives each mitochondrion as one solid object (mito_seg; no cristae label), so packs v1-v3 had only
the inner boundary membrane (IM = OM inset by IMS_W) and accounted for cristae with an area factor (CRISTAE_F = 3).
Here every mitochondrion (connected component of the same smoothed mask the mito field was built from) gets a stack of
lamellar cristae perpendicular to its long axis (PCA), as in most mammalian mitochondria:
  - crista lumen (intracristal space, labelled IMS): slab |u| < CR_W around planes spaced S apart along the axis
    (S drawn per mitochondrion from S_RANGE, random phase), ending CR_GAP inside the inner boundary membrane;
  - crista junctions: each lamella is joined to the boundary IMS by a ribbon of half-width CR_JW across the slab along the
    second principal axis through the centroid (two junctions per lamella);
  - matrix SDF = max(boundary, -lumen); IM SDF = -matrix (cell3d.mem_sdfs).
Mitochondria thinner than ~2 x CR_GAP get no lamellae (the gap rule removes them). Each component gets an id on a dense
ID_H grid (nearest component within ID_REACH) so the evaluator can find its parameters.
Values: intracristal lumen ~10-20 nm wide, IM midsurface separation across it ~20 nm (CR_W = 10); lamellae ~70-100 nm
apart; junction necks ~15-25 nm (Mannella 2006; Kühlbrandt 2015). The resulting IM/OM area ratio is reported and
replaces CRISTAE_F in the pack accounting (jwc.step_counts uses the measured IM area).
usage (wholecell): python -m toolkits.em.jcristae <src package dir> <new geometry_id> [seed=17]
(legacy: rewrote data/derived/jgeom_cell1.npz in place) """
import os, sys, json, numpy as np
from scipy import ndimage as ndi
from wholecell.geometry.engines import sdf_sections_v1 as c3
from toolkits.em import jgeom_build as jb
from toolkits.em import pkgstep
from toolkits.em._paths import inp
S_RANGE = (70.0, 100.0)
ID_H = 48.0
ID_REACH = 96.0
MIN_VOL = 0.004          # um3: smaller fragments get no cristae


def main(src, geometry_id, seed=17, cell=1):
    src_cell, parts = pkgstep.open_pkg(src)
    parts["cristae"] = None                                   # replace an earlier cristae section (legacy: from backup)
    G0 = c3.build_G(**parts, params=src_cell.manifest["engine"]["params"])
    z = {"pores": src_cell.pores}; centre = np.asarray(src_cell.manifest["frame"]["centre_label_nm"])
    # the mask of jgeom_build (isotropic 16 nm, 2-voxel pad, model frame), smoothed as for the mito field
    L = np.load(inp(f"cell{cell}_labels_s2.npz")); lab = L["labels"]; o2 = L["origin_s2_xyz"]
    Li = np.pad(ndi.zoom(lab, (1, 1, jb.VOXS2[2] / jb.H), order=0), 2)
    lo = o2 * jb.VOXS2 - centre - 2 * jb.H
    occ = ndi.gaussian_filter((Li == 3).astype(np.float32), jb.SPEC["mito"][1] / jb.H) > 0.5
    comp, n = ndi.label(occ)
    print(f"mitochondria: {n} components, {occ.sum() * jb.H ** 3 / 1e9:.1f} um3", flush=True)
    rng = np.random.default_rng(seed); par = np.zeros((n, 12))
    objs = ndi.find_objects(comp); nvox = ndi.sum(np.ones_like(comp), comp, np.arange(1, n + 1))
    for k in range(n):
        sl = objs[k]; idx = np.argwhere(comp[sl] == k + 1) + [s.start for s in sl]
        P = lo + idx * jb.H; c = P.mean(0)
        if len(P) > 3:
            w, U = np.linalg.eigh(np.cov((P - c).T)); a, v = U[:, 2], U[:, 1]
        else:
            a, v = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
        S = rng.uniform(*S_RANGE)
        par[k] = [*a, rng.uniform(0, S), S, *v, *c, float(nvox[k] * jb.H ** 3 / 1e9 >= MIN_VOL)]
    # id grid: nearest component within ID_REACH, on an ID_H grid over the mitochondria bounding box (+ reach)
    f = int(round(ID_H / jb.H))
    idx = np.argwhere(comp > 0); b0 = np.maximum(idx.min(0) - int(ID_REACH / jb.H) - f, 0); b1 = idx.max(0) + int(ID_REACH / jb.H) + f
    sub = comp[b0[0]:b1[0], b0[1]:b1[1], b0[2]:b1[2]]
    d, ind = ndi.distance_transform_edt(sub == 0, sampling=jb.H, return_indices=True)
    near = sub[tuple(ind)]; near[d > ID_REACH] = 0
    ids = near[f // 2::f, f // 2::f, f // 2::f].astype(np.float64) - 1.0          # -1 = none
    gorg = lo + (b0 + f // 2) * jb.H - 0.5 * ID_H                               # cell (i,j,k) covers gorg + [i, i+1) * ID_H
    sec = np.r_[n, gorg, ID_H, ids.shape, par.ravel(), ids.ravel()]
    parts["cristae"] = sec
    G = c3.build_G(**parts, params=src_cell.manifest["engine"]["params"])
    print(f"cristae section: {n} mitochondria ({int(par[:, 11].sum())} with cristae), id grid {ids.shape} at {ID_H} nm, "
          f"{sec.size / 1e6:.1f} M values", flush=True)
    # area check on a sample of mitochondria: IM vs OM area by marching cubes of the SDFs at 4 nm
    from skimage.measure import marching_cubes
    C = c3.Cell(G, z["pores"]); C0 = c3.Cell(G0, z["pores"]); rows = []
    big = np.argsort(-nvox)[:12]
    for k in big:
        sl = objs[k]; a0 = lo + np.array([s.start for s in sl]) * jb.H - 40; a1 = lo + np.array([s.stop for s in sl]) * jb.H + 40
        h = 4.0; ax = [np.arange(a0[i], a1[i], h) for i in range(3)]
        if np.prod([len(x) for x in ax]) > 3e7: continue
        X = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3); sh = tuple(len(x) for x in ax)
        r = {}
        for nm, cc in (("new", C), ("old", C0)):
            S_ = cc.mems(X).reshape(*sh, c3.NMEM)
            for m in (c3.M_OM, c3.M_IM):
                s = S_[..., m]
                if s.min() < 0 < s.max():
                    v, fc, _, _ = marching_cubes(s, 0.0, spacing=(h,) * 3)
                    r[(nm, m)] = float(np.linalg.norm(np.cross(v[fc[:, 1]] - v[fc[:, 0]], v[fc[:, 2]] - v[fc[:, 0]]), axis=1).sum() / 2)
            if nm == "new":
                labs = cc.labels(X); r["matrix_frac"] = float((labs == c3.MATRIX).sum() / max(((labs == c3.MATRIX) | (labs == c3.IMS)).sum(), 1))
        rows.append(dict(comp=int(k), vol_um3=round(float(nvox[k] * jb.H ** 3 / 1e9), 3), S=round(par[k, 4], 1),
                         IM_over_OM_new=round(r.get(("new", c3.M_IM), 0) / max(r.get(("new", c3.M_OM), 1), 1), 2),
                         IM_over_OM_old=round(r.get(("old", c3.M_IM), 0) / max(r.get(("old", c3.M_OM), 1), 1), 2),
                         matrix_frac_of_mito=round(r.get("matrix_frac", 0), 3)))
        print(rows[-1], flush=True)
    p = src_cell.manifest["engine"]["params"]
    summ = pkgstep.summary_copy(src_cell); summ["cristae"] = dict(n_mito=int(n), with_cristae=int(par[:, 11].sum()), S_range=S_RANGE,
                                                                  CR_W=p["CR_W"], CR_GAP=p["CR_GAP"], CR_JW=p["CR_JW"], sample=rows)
    man, res = pkgstep.publish(src_cell, parts, geometry_id, step="jcristae", summary=summ,
                               section_meta={"cristae": {"params": {"n_mito": int(n), "S_range": S_RANGE, "seed": seed},
                                                         "engine_params": ["CR_W", "CR_GAP", "CR_JW", "IMS_W"]}},
                               generator_extra={"seed": seed})
    print("written", geometry_id, "validation pass", res["pass"], flush=True)


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], *(int(x) for x in a[2:3]))
