# new in wholecell (B4, 2026-11-17): no legacy origin; the random fill is jgeom_build.main's pore fill (same rule)
"""Set the nuclear pore count of a Jurkat package (package N -> N+1; B4 of the unification plan, criteria
docs/unification_criteria/B4_preplaced.md J1/J1b).

The pore count follows the literature density checked against the nucleoporin anchor (DECISION 2026-11-16, Pedro):
N = round(rho x A_NE), rho = 11 NPC/um2 (Otsuka et al. 2016 eLife 5:e19071, HeLa interphase by EM; Maul & Deaven 1977),
A_NE = area of the NE midsurface mesh of the package. The pores of the parent are kept unchanged and first (EM np_seg
detections, then random); N - n_parent more are drawn as in jgeom_build: area-weighted random triangles of the NE
mesh, centre = triangle centroid, axis = triangle normal (outward), >= MIN_SEP from every pore, source code 0
(random). Everything else of the parent (fields, sections, mesh, layers) is carried unchanged (pkgstep.publish).

usage: NUMBA_NUM_THREADS=3 python -m toolkits.em.jpores [src={geom_root}/jurkat/em_cell1_v8] [geometry_id=em_cell1_v8p]
       [n_pores=2696] [seed=20261117]
"""
import copy
import os
import sys
import time

import numpy as np

RHO_PER_UM2 = 11.0          # Otsuka et al. 2016 eLife 5:e19071 (HeLa, steady-state interphase density by EM)
MIN_SEP = 150.0             # nm, as jgeom_build (MIN_SEP)
CODE = ("toolkits/em/jpores.py", "toolkits/em/pkgstep.py")


def fill(V, F, pores0, n_target, seed, min_sep=MIN_SEP, max_tries=2_000_000):
    """pores0 (k, 6) kept; random fill to n_target (jgeom_build rule); returns (pores (n, 6), tries)"""
    from scipy.spatial import cKDTree
    tri_c = V[F].mean(1)
    nrm = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    area_w = np.linalg.norm(nrm, axis=1)
    nrm = nrm / area_w[:, None]
    if np.mean(np.einsum("ij,ij->i", nrm, tri_c - V.mean(0))) < 0:          # outward (jgeom_build)
        nrm = -nrm
    area_w = area_w / area_w.sum()
    rng = np.random.default_rng(seed)
    sel = [np.asarray(p, np.float64) for p in pores0]
    pts = np.array([p[:3] for p in sel])
    grid = {}                                        # 150 nm cells: neighbour test without rebuilding trees
    h = min_sep

    def key(c):
        return tuple(np.floor(c / h).astype(int))

    for p in pts:
        grid.setdefault(key(p), []).append(p)
    tries = 0
    while len(sel) < n_target and tries < max_tries:
        k = rng.choice(len(F), p=area_w); tries += 1
        c = tri_c[k]; kc = key(c); ok = True
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for q in grid.get((kc[0] + dx, kc[1] + dy, kc[2] + dz), ()):
                        if np.linalg.norm(q - c) < min_sep:
                            ok = False; break
                    if not ok: break
                if not ok: break
            if not ok: break
        if not ok:
            continue
        sel.append(np.r_[c, nrm[k]]); grid.setdefault(kc, []).append(c)
    P = np.array(sel)
    d, _ = cKDTree(P[:, :3]).query(P[:, :3], k=2)
    return P, tries, float(d[:, 1].min())


def main(src="{geom_root}/jurkat/em_cell1_v8", geometry_id="em_cell1_v8p", n_pores=None, seed=20261117):
    import trimesh
    from toolkits.em import pkgstep
    t0 = time.time()
    cell, parts = pkgstep.open_pkg(src)
    man = cell.manifest
    mesh_rel = man["layers"]["ne_mesh"]["files"][0]
    m = trimesh.load(os.path.join(cell.package_dir, mesh_rel), process=False)
    V = np.asarray(m.vertices, np.float64) - np.asarray(man["frame"]["centre_label_nm"], np.float64)
    F = np.asarray(m.faces, np.int64)
    area = float(m.area / 1e6)
    n = int(n_pores) if n_pores else int(round(RHO_PER_UM2 * area))
    z = np.load(os.path.join(cell.package_dir, man["pores"]["file"]))
    p0, s0 = z["pores"], z["source"]
    P, tries, dmin = fill(V, F, p0, n, int(seed))
    assert np.array_equal(P[:len(p0)], p0), "parent pores changed"
    src_codes = np.r_[s0, np.zeros(len(P) - len(p0), s0.dtype)]
    print(f"pores: {len(p0)} kept ({int((s0 == 1).sum())} EM), {len(P) - len(p0)} added -> {len(P)} "
          f"(target {n} = {RHO_PER_UM2} /um2 x {area:.1f} um2; {tries} draws; min separation {dmin:.1f} nm)", flush=True)
    if len(P) != n:
        raise RuntimeError(f"pore fill reached {len(P)} of {n}")
    cell.manifest = man = copy.deepcopy(man)
    man["description"] = (man["description"].rstrip(".") + f"; nuclear pores set to {n} (B4: 11 NPC/um2 x NE "
                          "midsurface area, checked against the nucleoporin anchor)")
    L = man["layers"].get("npc_positions")
    if L is not None:
        L["counts"] = {"em": int((src_codes == 1).sum()), "random": int((src_codes == 0).sum())}
        L.setdefault("provenance", {})["b4_fill"] = (f"toolkits/em/jpores.py: + {len(P) - len(p0)} random pores (seed "
                                                     f"{seed}, >= {MIN_SEP:.0f} nm apart), density {RHO_PER_UM2}/um2")
    summ = pkgstep.summary_copy(cell)
    summ["pores"] = dict(n=int(len(P)), em=int((src_codes == 1).sum()), random=int((src_codes == 0).sum()),
                         density_per_um2=round(len(P) / area, 3), ne_area_um2=round(area, 1),
                         rule="11 NPC/um2 (Otsuka et al. 2016 eLife 5:e19071) x NE midsurface area; Nup anchor check in "
                              "docs/unification_criteria/B4_preplaced.md J1")
    man_new, res = pkgstep.publish(cell, parts, geometry_id, step="jpores", summary=summ, pores=P, pore_source=src_codes,
                                   waivers=man.get("validation_waivers"),
                                   generator_extra={"method": f"pores set to {n}: parent pores kept, random fill (jgeom_build "
                                                              f"rule, seed {seed}, min separation {MIN_SEP:.0f} nm)",
                                                    "criteria": "docs/unification_criteria/B4_preplaced.md",
                                                    "pore_seed": int(seed), "n_pores": int(n),
                                                    "code_sha256": {f: _sha(f) for f in CODE}})
    print(f"{man_new['package_id']}: validation pass {res['pass']}, failed {res.get('failed_checks')}, "
          f"waived {res.get('waived_failures')} ({time.time() - t0:.0f} s)", flush=True)
    return man_new, res


def _sha(rel):
    import hashlib
    from wholecell.config import REPO
    return hashlib.sha256(open(os.path.join(REPO, rel), "rb").read()).hexdigest()


if __name__ == "__main__":
    a = sys.argv[1:]
    main(*(a[:2]), *(int(x) for x in a[2:4]))
