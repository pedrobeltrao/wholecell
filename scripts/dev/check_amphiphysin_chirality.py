"""Chirality / placement check of the Rvs161-Rvs167 (amphiphysin) dimers of the yeast pre-pass (PORT PATCH 4, 2026-11-17).

The vendored builder (plugins/yeast/engine/pipeline/amphiphysin.py, place) oriented each dimer with R = W Rl^T, W =
[around, u, -radial] (left-handed: around = u x radial) and Rl = the library's principal axes (right-handed for this
entry), so det R = -1: an improper rotation, i.e. a mirror image. The residue beads it tested for clashes/containment
were R z + t (the mirrored dimer), and the stored quaternion was scipy 1.15 Rotation.from_matrix(R) (scipy 1.18 raises),
a proper rotation that is not the intended orientation. Patch 4: W = [around, -u, -radial], Rl right-handed by
construction -> det R = +1, same long axis (around the tube) and concave face (towards the membrane).

For every placed dimer of a pre-pass (s0 = before the patch; patched runs) this script reconstructs the placement frame
from the geometry (tube of the endocytic site, the dimer's centre) and reports:
  - det of the stored orientation (quaternion -> matrix; always +1) and of the builder's matrix (old: -1, new: +1),
  - that the stored quaternion is from_matrix(builder matrix) (old matrix: needs scipy < 1.18),
  - chirality vs the source structure (AFDB model AF-0000000210771937, the library entry = its pLDDT >= 50 residues,
    centred): (a) signed volume of the four half-chain centroids (Rvs161 N/C half, Rvs167 N/C half), (b) sign agreement
    of the signed volumes of 5000 fixed random bead quadruplets; internal distances (max |d_ij - d_ij(source)|),
  - placement: cos(long axis, around the tube), cos(concave axis, towards the tube axis = membrane), closest bead to
    the membrane (distance from the tube axis - tube radius - bilayer half thickness), sites, dimers per site.

    <yeast_x5 python> scripts/check_amphiphysin_chirality.py <out_dir> label=<prepass.npz or layers/amphiphysin.npz> ...
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

from _devsite import VC, RESULTS, YEAST as _YEAST_ROOT, PY as _SITE_PY, SMK as _SITE_SMK, GIT as _DEV_GIT  # noqa: E402,F401  (release: site paths)
S0 = f"{RESULTS}/yeast/20261116-aabe7b5-s0-x5repro"
LIB_F = f"{S0}/engine/lib/RVS161_RVS167_het.npz"
GEO_F = f"{S0}/engine/work/yc3d/geometry.json"
SRC_F = f"{_YEAST_ROOT}/inputs/afdb_yeast_dimers/heterodimers/AF-0000000210771937.pdb"
T_SLOT, HALF = 4, 2.0                                     # cell3d.T_SLOT, cell3d.HALF


def source():
    """CA coordinates (nm), chain and pLDDT of the source model; the library entry must be its pLDDT >= 50 part, centred"""
    xyz, ch, pl = [], [], []
    for line in open(SRC_F):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])]); ch.append(line[21]); pl.append(float(line[60:66]))
    xyz, ch, pl = np.array(xyz) / 10.0, np.array(ch), np.array(pl)
    keep = pl >= 50
    return xyz[keep], ch[keep]


def lib_frames(z):
    """the builder's library frame, original (Rl_old: eigenvectors as eigh returns them) and patched (mid = short x long)"""
    c0 = z.mean(0); w, V = np.linalg.eigh(np.cov((z - c0).T))
    al, am, ash = V[:, 2], V[:, 1], V[:, 0]
    s = (z - c0) @ al; t = (z - c0) @ ash
    ends, mid = np.abs(s) > 0.35 * np.abs(s).max(), np.abs(s) < 0.15 * np.abs(s).max()
    if t[ends].mean() < t[mid].mean():
        ash = -ash
    return c0, np.c_[al, am, ash], np.c_[al, np.cross(ash, al), ash]


def quat_mat(q):                                          # pack3d convention (w, x, y, z)
    w, x, y, zq = q / np.linalg.norm(q)
    return np.array([[1 - 2 * (y * y + zq * zq), 2 * (x * y - w * zq), 2 * (x * zq + w * y)],
                     [2 * (x * y + w * zq), 1 - 2 * (x * x + zq * zq), 2 * (y * zq - w * x)],
                     [2 * (x * zq - w * y), 2 * (y * zq + w * x), 1 - 2 * (x * x + y * y)]])


def chir_centroids(X, seg):
    c = [X[s].mean(0) for s in seg]
    return float(np.linalg.det(np.c_[c[1] - c[0], c[2] - c[0], c[3] - c[0]]))


def quad_signs(X, Q4):
    a, b, c, d = (X[Q4[:, k]] for k in range(4))
    return np.sign(np.einsum("ij,ij->i", b - a, np.cross(c - a, d - a)))


def analyse(label, path, z, seg, Q4, D0, sites, out):
    from scipy.spatial.transform import Rotation
    import scipy
    zz = np.load(path)
    pos, qs = zz["rvs_pos"].astype(np.float64), zz["rvs_q"].astype(np.float64)
    c0, Rl_old, Rl_new = lib_frames(z)
    s_src_c = chir_centroids(z, seg); s_src_q = quad_signs(z, Q4)
    tips, rt = sites[:, :3], sites[:, 6]; U = tips / np.linalg.norm(tips, axis=1)[:, None]
    rows = []
    for t, q in zip(pos, qs):
        # site: nearest tube axis (line from the tip outwards along u)
        rel = t[None] - tips; h = np.einsum("ij,ij->i", rel, U); perp = rel - h[:, None] * U
        dist = np.linalg.norm(perp, axis=1); dist[h < -5] = np.inf
        k = int(np.argmin(dist)); u = U[k]; radial = perp[k] / dist[k]; around = np.cross(u, radial)
        R_old = np.c_[around, u, -radial] @ Rl_old.T
        R_new = np.c_[around, -u, -radial] @ Rl_new.T
        Qm = quat_mat(q)
        r = dict(site=k, det_stored=float(np.linalg.det(Qm)), det_old=float(np.linalg.det(R_old)), det_new=float(np.linalg.det(R_new)))
        # which builder matrix produced the stored quaternion
        for nm, Rm in (("new", R_new), ("old", R_old)):
            try:
                qq = Rotation.from_matrix(Rm).as_quat(); qq = np.array([qq[3], qq[0], qq[1], qq[2]])
                r[f"q_from_{nm}"] = float(min(np.abs(qq - q).max(), np.abs(qq + q).max()))
            except ValueError:
                r[f"q_from_{nm}"] = None                   # scipy >= 1.18 rejects det -1
        r["angle_stored_vs_intended_deg"] = float(np.degrees(np.arccos(np.clip((np.trace(Qm.T @ R_new) - 1) / 2, -1, 1))))
        Bs = z @ Qm.T + t                                  # stored dimer as rendered/packed (proper)
        Bo = z @ R_old.T + t                               # what the old builder tested (det -1)
        for nm, B in (("stored", Bs), ("oldtest", Bo)):
            r[f"chir_centroids_{nm}"] = chir_centroids(B, seg)
            r[f"quad_agree_{nm}"] = float((quad_signs(B, Q4) == s_src_q).mean())
            d = np.linalg.norm(B[D0[0]] - B[D0[1]], axis=1); r[f"dist_err_{nm}"] = float(np.abs(d - D0[2]).max())
            rel_b = B - tips[k]; hb = rel_b @ u
            rb = np.linalg.norm(rel_b - hb[:, None] * u, axis=1)
            r[f"min_clear_{nm}"] = float((rb - rt[k] - HALF).min())
        r["cos_long_around"] = float((Qm @ Rl_new[:, 0]) @ around)
        r["cos_concave_membrane"] = float((Qm @ Rl_new[:, 2]) @ -radial)
        rows.append(r)
    import pandas as pd
    df = pd.DataFrame(rows); df.to_csv(f"{out}/instances_{label}.tsv", sep="\t", index=False)
    f = lambda c: [float(df[c].min()), float(df[c].median()), float(df[c].max())] if len(df) and df[c].notna().any() else None
    summ = dict(label=label, path=path, scipy=scipy.__version__, n=int(len(df)), sites_used=int(df.site.nunique()) if len(df) else 0,
                per_site=[int(x) for x in np.bincount(df.site, minlength=len(sites))] if len(df) else [],
                source_chir_centroids=s_src_c)
    for c in df.columns:
        if c != "site":
            summ[c] = f(c)
    summ["n_builder_matrix_improper"] = int((df.det_old < 0).sum()) if label.startswith("s0") else int((df.det_new < 0).sum())
    summ["n_stored_mirror"] = int((df.quad_agree_stored < 0.5).sum())
    summ["n_tested_beads_mirror"] = int((df.quad_agree_oldtest < 0.5).sum()) if label.startswith("s0") else 0
    return summ


def main():
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    z = np.load(LIB_F)["beads"].astype(np.float64)
    xs, ch = source()
    assert len(xs) == len(z), (len(xs), len(z))
    lib_vs_src = float(np.abs((xs - xs.mean(0)) - z).max())           # library = source, centred (no rotation)
    seg = []
    for c in sorted(set(ch)):
        idx = np.where(ch == c)[0]; seg += [idx[: len(idx) // 2], idx[len(idx) // 2:]]
    rng = np.random.default_rng(4)
    Q4 = rng.integers(0, len(z), (5000, 4)); Q4 = Q4[np.array([len(set(r)) == 4 for r in Q4])]
    P = rng.integers(0, len(z), (5000, 2)); D0 = (P[:, 0], P[:, 1], np.linalg.norm(z[P[:, 0]] - z[P[:, 1]], axis=1))
    geo = json.load(open(GEO_F)); PT = np.array(geo["PT"]); PQ = np.array(geo["PQ"], np.float64).reshape(-1, 8)
    sites = np.array([PQ[i] for i in np.where(PT == T_SLOT)[0] if np.linalg.norm(PQ[i, 3:6] - PQ[i, :3]) <= 1.0])
    c0, Rl_old, Rl_new = lib_frames(z)
    zm = (z - c0) @ np.diag([1, 1, -1.0])                             # a mirror image of the source, for the sign check
    rep = dict(source=SRC_F, library=LIB_F, lib_minus_centred_source_max_abs_nm=lib_vs_src,
               det_Rl_original=float(np.linalg.det(Rl_old)), det_Rl_patched=float(np.linalg.det(Rl_new)),
               Rl_patched_minus_original_max_abs=float(np.abs(Rl_new - Rl_old).max()),
               chir_centroids_source=chir_centroids(z, seg), chir_centroids_mirror_control=chir_centroids(zm, seg),
               quad_agree_mirror_control=float((quad_signs(zm, Q4) == quad_signs(z, Q4)).mean()),
               n_endocytic_sites=int(len(sites)), runs=[])
    for a in sys.argv[2:]:
        label, path = a.split("=", 1)
        rep["runs"].append(analyse(label, path, z, seg, Q4, D0, sites, out))
    json.dump(rep, open(f"{out}/chirality_check.json", "w"), indent=1)
    L = [f"source {os.path.basename(SRC_F)}: library = centred source (max |diff| {lib_vs_src:.2e} nm); det Rl original "
         f"{rep['det_Rl_original']:+.3f}, patched {rep['det_Rl_patched']:+.3f}",
         f"signed centroid volume: source {rep['chir_centroids_source']:+.2f} nm^3, mirror control {rep['chir_centroids_mirror_control']:+.2f}; "
         f"quadruplet sign agreement of the mirror control {rep['quad_agree_mirror_control']:.3f}", ""]
    cols = ["n", "sites_used", "det_stored", "det_old", "det_new", "q_from_old", "q_from_new", "chir_centroids_stored",
            "quad_agree_stored", "dist_err_stored", "chir_centroids_oldtest", "quad_agree_oldtest", "cos_long_around",
            "cos_concave_membrane", "angle_stored_vs_intended_deg", "min_clear_stored", "min_clear_oldtest",
            "n_builder_matrix_improper", "n_tested_beads_mirror", "n_stored_mirror"]
    for r in rep["runs"]:
        L.append(f"== {r['label']} ({r['path']}; scipy {r['scipy']})")
        for c in cols:
            v = r.get(c)
            L.append(f"  {c:30s} " + (("min/median/max " + " / ".join(f"{x:.4g}" for x in v)) if isinstance(v, list) else str(v)))
        L.append(f"  dimers per site: {r['per_site']}")
    open(f"{out}/chirality_check.txt", "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
