# origin: jurkat/yc/pipeline/jlamina.py (legacy workspace, 2026-11-15)
# Changed: registered preplaced builder "lamina_lattice" (generic; mammalian-lamin numbers are config params): species
# list, patch size, pore exclusion, lateral factor, priority (env LAMINA_PRIORITY / LAMINA_LAT_F before) and the other
# literals are params with the legacy values as defaults; seed from wholecell.seeds stream "prepass.lamina" (Jurkat:
# legacy 20261112); library npz from the build's library manifest; geometry (G, pores, centre) from the build geometry,
# MEM_HALF[NE] and PERI_W from G; NE mesh path from inputs; outputs written to the module folder of the build, json
# paths absolute. Numerics unchanged.
"""Nuclear lamina lattice pre-pass (two crossing layers of lamin rods under the INM; see legacy docstring, rule L2).

Each library is used in its principal frame (smallest axis along the local INM normal); the INM is divided into
patches (Voronoi cells of seeds patch_nm apart on the NE midsurface mesh, random in-plane orientation per patch);
layer 1 rods along e1 at depth D1 below the INM midsurface, layer 2 along e2 = n x e1 at D2; rods head-to-tail every
PERIOD nm, filaments LAT nm apart (staggered rows); sites projected onto the depth isosurface of the NE field; dropped
when centre or ends are not nucleoplasm, near a pore axis, too close to the INM, or overlapping a rod of a neighbouring
patch; species assigned at random by copies, thinned to the target. priority=chromatin drops rods clashing with
chromatin pre-pass instances; priority=lamina keeps them and lists the chromatin instances to drop at indexing.
Outputs (module folder): lamina_cell1.npz (lib_ids, types, pos, quat, layer), lamina_cell1.json (consumed,
consumed_key, ...), chromatin_drop_lamina.npy (+ .json: target module, n, why).
"""
import json
import os
import time

import numpy as np
from scipy.spatial import cKDTree

from wholecell.registry import register

DEFAULTS = dict(genes=["LMNB1", "LMNB2", "LMNA"], inventory_key="lamina", patch_nm=500.0, pore_x_nm=80.0,
                pore_z_nm=150.0, lat_f=1.0, priority="lamina", period_pad_nm=0.5, depth_pad_nm=0.5,
                border_nm=35.0, end_inset_nm=1.5, chrom_inm_nm=45.0, capsule_inset_nm=1.2, clash_tol_nm=0.5,
                overlap_f=0.95, newton_it=3, newton_h=0.5, chunk=2_000_000, out_prefix="lamina_cell1",
                drop_name="chromatin_drop_lamina")


def frames(names, LIB):
    out = {}
    for n in names:
        z = np.load(LIB[n]); c = z["st_c"].astype(np.float64); r = z["st_r"].astype(np.float64)
        w = r ** 3; m = (c * w[:, None]).sum(0) / w.sum(); X = c - m
        _, U = np.linalg.eigh((X * w[:, None]).T @ X); U = U[:, [2, 1, 0]]
        if np.linalg.det(U) < 0: U[:, 2] *= -1
        Y = X @ U
        bc = 0.5 * ((Y + r[:, None]).max(0) + (Y - r[:, None]).min(0))   # bounding-box centre in the principal frame
        Y = Y - bc; m = m + U @ bc                                       # rods are framed on their bounding box
        ext = np.abs(Y).max(0) + r.max()
        out[n] = dict(m=m, U=U, Y=Y, r=r, half_len=float(ext[0]), half_w=float(ext[1]), half_thick=float(ext[2]))
    return out


def project(cell, M_NE, P, depth, it=3, h=0.5):
    """move points onto the NE-distance isosurface out[M_NE] = depth on the nucleoplasm side; returns points and unit
    normals pointing into the nucleus (direction of increasing distance)"""
    P = P.copy()
    for _ in range(it):
        S = cell.mems(P)[:, M_NE].astype(np.float64)
        g = np.stack([(cell.mems(P + h * e)[:, M_NE] - cell.mems(P - h * e)[:, M_NE]) / (2 * h) for e in np.eye(3)], 1)
        gn = np.linalg.norm(g, axis=1, keepdims=True); g /= np.maximum(gn, 1e-9)
        P -= (S - depth)[:, None] * g
    S = cell.mems(P)[:, M_NE].astype(np.float64)
    good = (np.abs(S - depth) < 1.0) & (gn[:, 0] > 0.5)          # converged, and inside the stored NE band
    return P, g, good


def build(geom, LIB, chrom, ing, ne_mesh, out_dir, seed, params, chrom_name="chromatin"):
    """geom: wholecell.pack.geomload.Geom; LIB: {species: npz}; chrom: chromatin pre-pass npz; ing: inventory table
    (pre-lamina rule step); ne_mesh: NE midsurface mesh (source-volume frame); returns the json dict"""
    import pandas as pd
    import trimesh
    from wholecell.pack.regions import label_map
    pp = dict(DEFAULTS); pp.update(params or {})
    GENES = list(pp["genes"]); PATCH = float(pp["patch_nm"]); PORE_X = float(pp["pore_x_nm"]); LAT_F = float(pp["lat_f"])
    PRIORITY = pp["priority"]; KEY = pp["inventory_key"]
    c3 = geom.engine; cell = geom.cell; gpar = geom.params
    M_NE = c3.M_NE; NUCP = label_map(geom)["nucleoplasm"]
    MH_NE = float(gpar["MEM_HALF"][M_NE]); PERI_W = float(gpar["PERI_W"])
    t0 = time.time(); rng = np.random.default_rng(seed)
    centre = geom.centre_nm
    I = pd.read_csv(ing, sep="\t").set_index("ing_id")
    tgt = {g: int(round(I.at[g, "copies"] * json.loads(I.at[g, "comp_weights"]).get(KEY, 0.0))) for g in GENES}
    ntot = sum(tgt.values())
    Fm = frames(GENES, LIB)
    hlm = max(Fm[g]["half_len"] for g in GENES); hwm = max(Fm[g]["half_w"] for g in GENES); htm = max(Fm[g]["half_thick"] for g in GENES)
    PERIOD = 2 * hlm + pp["period_pad_nm"]; LAT_MIN = LAT_F * 2 * hwm               # box-framed neighbours
    D1 = MH_NE + htm + pp["depth_pad_nm"]; D2 = D1 + 2 * htm + pp["depth_pad_nm"]
    m = trimesh.load(ne_mesh, process=False)
    Vt = np.asarray(m.vertices, np.float64) - centre; area = float(m.area)            # NE midsurface ~ INM area
    Vn = np.asarray(m.vertex_normals, np.float64)
    if np.mean(np.einsum("ij,ij->i", Vn, Vt - Vt.mean(0))) < 0: Vn = -Vn           # outward (to the cytoplasm)
    lat = max(LAT_MIN, 2 * area / ntot / PERIOD)                                     # two layers
    print(f"lamina copies {tgt} (total {ntot:,}); NE midsurface {area / 1e6:.1f} um2; lateral spacing {lat:.2f} nm "
          f"(needed {2 * area / ntot / PERIOD:.2f})", flush=True)
    # patch seeds: greedy subsample of mesh vertices PATCH apart
    order = rng.permutation(len(Vt)); seeds = []; snorm = []
    tr_v = None
    for i in order:
        p = Vt[i]
        if seeds and tr_v is not None and tr_v.query(p)[0] < PATCH: continue
        if seeds and np.min(np.linalg.norm(np.array(seeds[-200:]) - p, axis=1)) < PATCH: continue
        seeds.append(p); snorm.append(Vn[i])
        if len(seeds) % 200 == 0: tr_v = cKDTree(np.array(seeds))
    seeds = np.array(seeds); stree = cKDTree(seeds)
    print(f"{len(seeds)} patches ({time.time() - t0:.0f}s)", flush=True)
    sites, dirs, layer, patch = [], [], [], []
    half = 1.1 * PATCH
    for k, p in enumerate(seeds):
        n = -snorm[k]                                                              # into the nucleus
        a = rng.normal(size=3); e1 = a - (a @ n) * n; e1 /= np.linalg.norm(e1); e2 = np.cross(n, e1)
        for L, (d, u) in enumerate(((D1, e1), (D2, e2))):
            v = np.cross(n, u)
            us = np.arange(-half, half, PERIOD); vs = np.arange(-half, half, lat)
            U_, V_ = np.meshgrid(us, vs, indexing="ij"); U_ = U_ + (np.arange(len(vs))[None, :] % 2) * PERIOD / 2
            Q = p + (0.5 * PERI_W + d) * n + U_.ravel()[:, None] * u + V_.ravel()[:, None] * v
            sites.append(Q); dirs.append(np.repeat(u[None], len(Q), 0)); layer.append(np.full(len(Q), L)); patch.append(np.full(len(Q), k))
    S = np.concatenate(sites); Dr = np.concatenate(dirs); Ly = np.concatenate(layer); Pt = np.concatenate(patch)
    print(f"{len(S):,} raw sites ({time.time() - t0:.0f}s)", flush=True)
    # project to depth, Voronoi assignment
    P = np.empty_like(S); N = np.empty_like(S); G_ok = np.zeros(len(S), bool)
    ch_ = int(pp["chunk"])
    for L, d in ((0, D1), (1, D2)):
        s = Ly == L
        for a in range(0, int(s.sum()), ch_):
            idx = np.nonzero(s)[0][a:a + ch_]
            P[idx], N[idx], G_ok[idx] = project(cell, M_NE, S[idx], d, int(pp["newton_it"]), float(pp["newton_h"]))
    dd, kk = stree.query(P, k=2)
    keep = (kk[:, 0] == Pt) & G_ok
    P, N, Dr, Ly, Pt, dd = P[keep], N[keep], Dr[keep], Ly[keep], Pt[keep], dd[keep]
    Dr = Dr - np.einsum("ij,ij->i", Dr, N)[:, None] * N; Dr /= np.linalg.norm(Dr, axis=1, keepdims=True)
    border = (dd[:, 1] - dd[:, 0]) < 2 * pp["border_nm"]
    print(f"{len(P):,} sites after Voronoi assignment ({int(border.sum()):,} at patch borders) ({time.time() - t0:.0f}s)", flush=True)
    # species (random by copies), thinning to the target
    sp = rng.choice(len(GENES), size=len(P), p=np.array([tgt[g] for g in GENES]) / ntot)
    hl = np.array([Fm[g]["half_len"] for g in GENES])[sp]; ht = np.array([Fm[g]["half_thick"] for g in GENES])[sp]
    # geometry checks: centre and both rod ends nucleoplasm, depth, pores
    ok = np.ones(len(P), bool); reason = {}
    ei = pp["end_inset_nm"]
    for nm_, Q in (("centre", P), ("end+", P + (hl - ei)[:, None] * Dr), ("end-", P - (hl - ei)[:, None] * Dr)):
        lab = cell.labels(Q); Sq = cell.mems(Q)
        bad = (lab != NUCP) | (Sq[:, M_NE] - np.float64(MH_NE) < ht - pp["clash_tol_nm"])   # float64 as legacy (np scalar)
        reason[nm_] = int((bad & ok).sum()); ok &= ~bad
    pr = cell.pores; ptree = cKDTree(pr[:, :3])
    near = ptree.query_ball_point(P, PORE_X + 60.0)
    bad = np.zeros(len(P), bool)
    for i, js in enumerate(near):
        if not js or not ok[i]: continue
        v = P[i] - pr[js, :3]; zz = np.einsum("ij,ij->i", v, pr[js, 3:]); rho = np.sqrt(np.maximum((v ** 2).sum(1) - zz ** 2, 0))
        bad[i] = bool(((rho < PORE_X) & (np.abs(zz) < pp["pore_z_nm"])).any())
    reason["pore"] = int(bad.sum()); ok &= ~bad
    # chromatin: instances near the INM; capsule (rod axis +-(hl - inset), radius ht) vs bounding sphere
    ch = np.load(chrom); cid = list(ch["lib_ids"]); cbr = []
    for l in cid:
        zl = np.load(LIB[l]); cbr.append(float((np.linalg.norm(zl["st_c"], axis=1) + zl["st_r"]).max()))
    cbr = np.array(cbr)
    cp = ch["pos"].astype(np.float64); sel = np.zeros(len(cp), bool)
    lo, hi = P.min(0) - 50, P.max(0) + 50
    inb = np.all((cp > lo) & (cp < hi), axis=1)
    idx = np.nonzero(inb)[0]
    for a in range(0, len(idx), 5_000_000):
        j = idx[a:a + 5_000_000]; sel[j] = cell.mems(cp[j])[:, M_NE] < pp["chrom_inm_nm"]
    cP = cp[sel]; cR = cbr[ch["types"][sel]]
    print(f"chromatin instances within {pp['chrom_inm_nm']} nm of the INM: {int(sel.sum()):,} of {len(cp):,} ({time.time() - t0:.0f}s)", flush=True)
    cidx = np.nonzero(sel)[0]
    ctree = cKDTree(cP); bad = np.zeros(len(P), bool); hits = [[] for _ in range(len(P))]
    rq = float(hl.max() + cR.max() + ht.max())
    ci = pp["capsule_inset_nm"]
    for i, js in enumerate(ctree.query_ball_point(P, rq)):
        if not js or not ok[i]: continue
        v = cP[js] - P[i]; t = np.clip(v @ Dr[i], -(hl[i] - ci), hl[i] - ci)
        d = np.linalg.norm(v - t[:, None] * Dr[i], axis=1)
        hb = d < cR[js] + ht[i] - pp["clash_tol_nm"]
        if hb.any(): bad[i] = True; hits[i] = list(cidx[np.asarray(js)[hb]])
    if PRIORITY == "chromatin":
        reason["chromatin"] = int(bad.sum()); ok &= ~bad; chrom_drop = np.zeros(0, np.int64)
    else:                                   # lamina first: clashing chromatin instances are dropped at pack indexing
        reason["chromatin_clash_sites_kept"] = int(bad.sum()); chrom_drop = None
    hits = [hits[i] for i in np.nonzero(ok)[0]]
    P, N, Dr, Ly, Pt, sp, border = P[ok], N[ok], Dr[ok], Ly[ok], Pt[ok], sp[ok], border[ok]
    W1 = Dr; W3 = N; W2 = np.cross(W3, W1)
    Rw = np.stack([W1, W2, W3], 2)                          # columns = world images of the principal axes
    # overlap resolution, bead level: patches in order, each rod accepted when it overlaps neither an accepted rod of an
    # earlier patch nor an earlier-accepted rod of its own patch
    of = pp["overlap_f"]
    Yb = {s_: (Fm[g]["Y"], Fm[g]["r"]) for s_, g in enumerate(GENES)}
    rmx = max(float(Fm[g]["r"].max()) for g in GENES)
    acc = np.zeros(len(P), bool); accB = {}; order_p = np.argsort(Pt, kind="stable"); bounds = np.searchsorted(Pt[order_p], np.arange(len(seeds) + 1))
    for kp in range(len(seeds)):
        ii = order_p[bounds[kp]:bounds[kp + 1]]
        if not len(ii): continue
        rng.shuffle(ii)
        B = np.concatenate([P[i] + Yb[sp[i]][0] @ Rw[i].T for i in ii]); Rb = np.concatenate([Yb[sp[i]][1] for i in ii])
        Ob = np.concatenate([np.full(len(Yb[sp[i]][1]), j) for j, i in enumerate(ii)])
        bad = np.zeros(len(ii), bool)
        nb_ = [q for q in stree.query_ball_point(seeds[kp], 3.0 * PATCH) if q in accB]
        if nb_:
            AB = np.concatenate([accB[q][0] for q in nb_]); AR = np.concatenate([accB[q][1] for q in nb_])
            tA = cKDTree(AB); d0, _ = tA.query(B, k=1)
            for j in np.nonzero(d0 < (Rb + rmx) * of)[0]:
                if bad[Ob[j]]: continue
                js = tA.query_ball_point(B[j], Rb[j] + rmx)
                if js and (np.linalg.norm(AB[js] - B[j], axis=1) < (AR[js] + Rb[j]) * of).any(): bad[Ob[j]] = True
        prs = cKDTree(B).query_pairs(2 * rmx, output_type="ndarray")
        if len(prs):
            x_, y_ = prs[:, 0], prs[:, 1]
            hit = (Ob[x_] != Ob[y_]) & (np.linalg.norm(B[x_] - B[y_], axis=1) < (Rb[x_] + Rb[y_]) * of)
            conf = {}
            for u_, v_ in zip(Ob[x_][hit], Ob[y_][hit]):
                conf.setdefault(max(u_, v_), set()).add(min(u_, v_))
            ok_ = ~bad
            for j in range(len(ii)):                          # greedy in the shuffled order
                if ok_[j] and any(ok_[q] for q in conf.get(j, ())): ok_[j] = False
        else:
            ok_ = ~bad
        acc[ii[ok_]] = True
        keepb = ok_[Ob]; accB[kp] = (B[keepb], Rb[keepb])
    reason["rod_overlap"] = int((~acc).sum())
    hits = [hits[i] for i in np.nonzero(acc)[0]]
    P, N, Dr, Ly, Pt, sp, Rw = P[acc], N[acc], Dr[acc], Ly[acc], Pt[acc], sp[acc], Rw[acc]
    # thin to the target per species
    keep = np.zeros(len(P), bool)
    for s_, g in enumerate(GENES):
        ii = np.nonzero(sp == s_)[0]; rng.shuffle(ii); keep[ii[:tgt[g]]] = True
    hits = [hits[i] for i in np.nonzero(keep)[0]]
    P, N, Dr, Ly, Pt, sp, Rw = P[keep], N[keep], Dr[keep], Ly[keep], Pt[keep], sp[keep], Rw[keep]
    Ulib = np.stack([Fm[g]["U"] for g in GENES])[sp]; mlib = np.stack([Fm[g]["m"] for g in GENES])[sp]
    R = Rw @ np.transpose(Ulib, (0, 2, 1))                    # world <- library
    pos = P - np.einsum("nij,nj->ni", R, mlib)
    k = np.ones(len(P), bool)
    if chrom_drop is None:
        chrom_drop = np.unique(np.array([x for i in np.nonzero(k)[0] for x in hits[i]], np.int64))
    pos, R, sp, Ly = pos[k], R[k], sp[k], Ly[k]
    from scipy.spatial.transform import Rotation
    q = Rotation.from_matrix(R).as_quat(); quat = np.c_[q[:, 3:], q[:, :3]].astype(np.float32)
    os.makedirs(out_dir, exist_ok=True)
    pre, dn = pp["out_prefix"], pp["drop_name"]
    np.savez_compressed(f"{out_dir}/{pre}.npz", lib_ids=np.array(GENES), types=sp.astype(np.int32), pos=pos.astype(np.float32),
                        quat=quat, layer=Ly.astype(np.int8))
    placed = {g: int((sp == s_).sum()) for s_, g in enumerate(GENES)}
    np.save(f"{out_dir}/{dn}.npy", chrom_drop)
    json.dump(dict(target=chrom_name, target_npz=os.path.abspath(chrom), n=int(len(ch["types"])),
                   drop=f"{dn}.npy",                     # relative to this json's folder
                   why="overlap with the lamina lattice (lamina_lattice, lamina priority)"),
              open(f"{out_dir}/{dn}.json", "w"), indent=1)
    cl = list(ch["lib_ids"]); ct = ch["types"][chrom_drop] if len(chrom_drop) else np.zeros(0, int)
    reason["chromatin_instances_dropped"] = int(len(chrom_drop))
    reason["nucleosomes_dropped"] = int((ct == cl.index("NUCLEOSOME")).sum()) if "NUCLEOSOME" in cl else 0
    J = dict(version="wholecell lamina_lattice", consumed=placed, consumed_key=KEY, target=tgt,
             lib_to_ingredient={g: g for g in GENES}, geometry=str(geom.source),
             chromatin=os.path.abspath(chrom), chromatin_n=int(len(ch["types"])), priority=PRIORITY,
             chromatin_drop=f"{dn}.npy (indices into the chromatin file above; pack index drops them only when the instance count matches)",
             ingredients=os.path.abspath(ing), lateral_spacing_nm=round(lat, 2),
             params=dict(LAT_F=LAT_F, PERIOD=PERIOD, LAT_MIN=LAT_MIN, PATCH=PATCH, D1=D1, D2=D2, PORE_X=PORE_X, seed=int(seed)),
             dropped=reason, patches=int(len(seeds)), placed_frac=round(sum(placed.values()) / max(ntot, 1), 4),
             rods={g: dict(half_len=round(Fm[g]["half_len"], 2), half_thick=round(Fm[g]["half_thick"], 2)) for g in GENES},
             density_per_um2=round(sum(placed.values()) / (area / 1e6), 1), runtime_s=round(time.time() - t0))
    json.dump(J, open(f"{out_dir}/{pre}.json", "w"), indent=1)
    print(json.dumps(J, indent=1), flush=True)
    return J


register("preplaced", "lamina_lattice")(build)
