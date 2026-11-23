# origin: jurkat/yc/pipeline/jspeckles.py (legacy workspace, 2026-11-15)
# restructured: reads package N, writes package N+1 (toolkits/em/pkgstep.py); paths from config; NOT RUN in wholecell
"""Nuclear speckles for the Jurkat cell (pack v6, PLANNING_NOTES "v6 implementation", rule R6b).

There is no speckle label in the EM data, so speckles are procedural: N_SPK spheres of radius R_SPK placed in the
interchromatin space of the nucleoplasm (label NUCP), at the most chromatin-poor sites (nucleosome pre-pass density
smoothed over the sphere volume), >= D_EDGE from anything that is not nucleoplasm (NE, pores, nucleolus),
>= D_SEP apart (centre to centre), every sphere fully inside NUCP (100 surface probes).
The section n, n x (cx, cy, cz, r) is appended to G and G[cell3d.G_SPECKLE] points at it; cell3d labels SPK inside.
Literature: 20-50 speckles per mammalian nucleus, 0.5-2 um, in interchromatin regions (Spector & Lamond 2011, CSH Perspect
Biol 3:a000646). We use the small end (0.5 um) as the volume the packer can fill at a moderate density.

usage (wholecell): NUMBA_NUM_THREADS=8 python -m toolkits.em.jspeckles <src package> <new geometry_id> <chromatin npz>
       [n=30] [r_nm=250] [seed=5]   (v7: chromatin = chromatin/v8/prepass/chromatin_cell1.npz)
legacy usage: NUMBA_NUM_THREADS=8 python jspeckles.py [n=30] [r_nm=250] [seed=5]
v7: JGEOM / JCHROM env select the geometry file (rewritten in place) and the chromatin pre-pass; a Golgi cisternae
section (jgolgi.py) after the speckle section is preserved. Backs up nothing itself (run after copying jgeom_cell1.npz to jgeom_cell1_v5.npz); idempotent: an existing speckle
section is replaced (G is rebuilt from the part before it)."""
import os, sys, json, time
import numpy as np
from wholecell.geometry.engines import sdf_sections_v1 as c3
from toolkits.em import pkgstep
from toolkits.em._paths import inp

A = sys.argv[1:]
N_SPK = int(A[3]) if len(A) > 3 else 30
R_SPK = float(A[4]) if len(A) > 4 else 250.0
SEED = int(A[5]) if len(A) > 5 else 5
D_EDGE, D_SEP, H = 350.0, 1000.0, 100.0


def main(src, geometry_id, chrom=None):
    t0 = time.time()
    CHROM = chrom or inp("prepass/chromatin_cell1.npz")
    src_cell, parts = pkgstep.open_pkg(src); summ = pkgstep.summary_copy(src_cell)
    params = src_cell.manifest["engine"]["params"]
    parts["speckle"] = None                                    # replace an earlier speckle section (Golgi section kept)
    G = c3.build_G(**parts, params=params)
    z = {"pores": src_cell.pores}
    cell = c3.Cell(G, z["pores"])
    ch = np.load(CHROM); ids = list(ch["lib_ids"])
    nuc = ch["pos"][ch["types"] == ids.index("NUCLEOSOME")].astype(np.float64)
    lo = nuc.min(0) - 500; hi = nuc.max(0) + 500; n = np.ceil((hi - lo) / H).astype(int)
    print(f"{len(nuc):,} nucleosomes; grid {n.tolist()} at {H:.0f} nm", flush=True)
    hist = np.histogramdd(nuc, bins=n, range=list(zip(lo, lo + n * H)))[0].astype(np.float32)
    from scipy.ndimage import uniform_filter
    k = int(round(2 * R_SPK / H)) | 1
    dens = uniform_filter(hist, size=k, mode="constant")       # nucleosomes per H^3 averaged over the sphere's cube
    g = [lo[i] + (np.arange(n[i]) + 0.5) * H for i in range(3)]
    X, Y, Zg = np.meshgrid(*g, indexing="ij"); P = np.stack([X.ravel(), Y.ravel(), Zg.ravel()], 1)
    lab = cell.labels(P)
    cand = np.nonzero(lab == c3.NUCP)[0]
    # clearance from everything that is not nucleoplasm (NE / perinuclear / pores, nucleolus, cytoplasm): exact EDT on the
    # 100 nm label grid (the stored SDFs are clamped at 64-160 nm, too short for this test)
    from scipy.ndimage import distance_transform_edt
    clr = distance_transform_edt((lab == c3.NUCP).reshape(tuple(n))).ravel() * H
    cand = cand[clr[cand] >= D_EDGE + R_SPK]
    d = dens.ravel()[cand]
    rng = np.random.default_rng(SEED)
    order = cand[np.lexsort((rng.random(len(cand)), d))]       # most chromatin-poor first, random tie-break
    print(f"{len(cand):,} candidate centres; nucleosome density (per {H:.0f} nm voxel, sphere-averaged) "
          f"min {d.min():.2f}, median {np.median(d):.2f}, all-nucleus median {np.median(dens.ravel()[lab == c3.NUCP]):.2f}", flush=True)
    u = rng.normal(size=(100, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    sel = []
    for i in order:
        c = P[i]
        if any(np.linalg.norm(c - P[j]) < D_SEP for j in sel): continue
        if not np.all(cell.labels(c + 0.98 * R_SPK * u) == c3.NUCP): continue
        sel.append(i)
        if len(sel) == N_SPK: break
    C = P[sel]
    sec = np.r_[float(len(sel)), np.c_[C, np.full(len(sel), R_SPK)].ravel()]
    parts["speckle"] = sec
    G2 = c3.build_G(**parts, params=params)
    cell2 = c3.Cell(G2, z["pores"])
    # check: label volume by Monte Carlo inside the spheres and chromatin load
    u3 = rng.normal(size=(len(C), 400, 3)); q = (C[:, None, :] + R_SPK * np.cbrt(rng.random((len(C), 400, 1))) * u3 /
         np.linalg.norm(u3, axis=2, keepdims=True)).reshape(-1, 3)
    frac_spk = float((cell2.labels(q) == c3.SPK).mean())
    vol = len(sel) * 4 / 3 * np.pi * (R_SPK / 1000) ** 3
    inside = np.zeros(len(nuc), bool)
    from scipy.spatial import cKDTree
    tr = cKDTree(C); dd, _ = tr.query(nuc, k=1); inside = dd < R_SPK
    nuc_dens_spk = inside.sum() / vol
    nuc_dens_all = len(nuc) / (summ.get("nucleus_um3") or 254.5)
    S = dict(n=len(sel), r_nm=R_SPK, volume_um3=round(vol, 3), label_check_frac=round(frac_spk, 4), seed=SEED,
             d_edge_nm=D_EDGE, sep_nm=D_SEP, centres_nm=np.round(C, 1).tolist(),
             nucleosomes_inside=int(inside.sum()), nucleosome_density_inside_per_um3=round(float(nuc_dens_spk), 1),
             nucleosome_density_nucleus_per_um3=round(float(nuc_dens_all), 1),
             source="procedural (no EM label): most chromatin-poor NUCP sites; jspeckles.py")
    summ["speckles"] = S
    man, res = pkgstep.publish(src_cell, parts, geometry_id, step="jspeckles", summary=summ,
                               section_meta={"speckle": {"params": {k: S[k] for k in ("n", "r_nm", "seed", "d_edge_nm", "sep_nm", "source")},
                                                         "label": "nuclear_speck", "chromatin": CHROM}},
                               generator_extra={"seed": SEED, "chromatin": CHROM})
    print(json.dumps({k: v for k, v in S.items() if k != "centres_nm"}, indent=1), f"\n{time.time() - t0:.0f}s",
          geometry_id, "validation pass", res["pass"], flush=True)


if __name__ == "__main__":
    main(A[0], A[1], A[2] if len(A) > 2 else None)
