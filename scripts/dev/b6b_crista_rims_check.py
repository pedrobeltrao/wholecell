# new in wholecell (B6b, 2026-11-17)
"""B6b rounded crista rims: measurements for docs/unification_criteria/B6b_crista_rims.md (Jurkat, pre-build).

  geom   R2-R8: package em_cell1_v8r vs its parent em_cell1_v8p (files, manifest, validator, MC volumes, paired label /
         membrane changes, volume moved, IM area).
         python scripts/b6b_crista_rims_check.py geom --out <json> [--n 20000000]
  atp    R9 / R10: core.atp_synthase_dimers (rims lamellae) on a package; cell-wide at the given sampling, and saturated
         on a sub-region (mitochondria --comps) with the rim length of those lamellae (band method).
         python scripts/b6b_crista_rims_check.py atp --bdir <dir with inventory/ library/> --out <json> \
                --package em_cell1_v8r [--params JSON] [--comps 3] [--sat-per-site 5]
NUMBA_NUM_THREADS / OMP_NUM_THREADS set the cores (3 for checks).
"""
import argparse
import copy
import hashlib
import json
import os
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _load(pid):
    from wholecell import config
    from wholecell.config import load_site
    from wholecell.geometry import cellgeom
    cfg = config.load("jurkat")
    return cellgeom.load_cfg(cfg, f"{load_site()['geom_root']}/jurkat/{pid}")


class Lamellae:
    """the cristae section of a G: id grid, per-lamella params, enabled cells"""
    def __init__(self, cell):
        e = cell.engine; G = cell.G
        o = int(G[e.G_CRISTAE]); n = int(G[o])
        self.lo = G[o + 1:o + 4].copy(); self.h = float(G[o + 4]); self.dims = G[o + 5:o + 8].astype(int)
        self.par = G[o + 8:o + 8 + 12 * n].reshape(n, 12).copy()
        self.ids = G[o + 8 + 12 * n:o + 8 + 12 * n + self.dims.prod()].reshape(self.dims).astype(int)
        self.W, self.GAP, self.JW, self.IMSW = (float(G[s]) for s in (e.S_CR_W, e.S_CR_GAP, e.S_CR_JW, e.S_IMS_W))
        en = (self.ids >= 0) & (self.par[np.maximum(self.ids, 0), 11] > 0)
        self.cells = np.argwhere(en)

    def sample(self, rng, m, cells=None):
        c = self.cells if cells is None else cells
        return self.lo + (c[rng.integers(0, len(c), m)] + rng.random((m, 3))) * self.h

    def coords(self, X, om):
        """lamella id, u (distance to the mid-plane), e = b + CR_GAP, w = |junction coordinate| - CR_JW"""
        ii = np.floor((X - self.lo) / self.h).astype(int); k = self.ids[ii[:, 0], ii[:, 1], ii[:, 2]]
        q = self.par[k]; S = q[:, 4]
        t = (q[:, :3] * X).sum(1) - q[:, 3]; u = t - S * np.floor(t / S + 0.5)
        w = np.abs((q[:, 5:8] * (X - q[:, 8:11])).sum(1)) - self.JW
        return k, u, om + self.IMSW + self.GAP, w


def geom(a):
    from wholecell.config import load_site
    gr = load_site()["geom_root"]
    P0, P1 = f"{gr}/jurkat/{a.parent}", f"{gr}/jurkat/{a.package}"
    m0, m1 = json.load(open(f"{P0}/manifest.json")), json.load(open(f"{P1}/manifest.json"))
    res = {"parent": a.parent, "package": a.package}
    # R2
    f0 = sorted(m0["files"]); f1 = sorted(m1["files"])
    same = {f: (_sha(f"{P0}/{f}") == _sha(f"{P1}/{f}")) for f in f0 if f in f1}
    diff_keys = sorted(k for k in set(m0) | set(m1) if m0.get(k) != m1.get(k))
    p0, p1 = m0["engine"]["params"], m1["engine"]["params"]
    res["R2"] = dict(content_files=len(same), identical=int(sum(same.values())), different=[f for f, s in same.items() if not s],
                     only_parent=sorted(set(f0) - set(f1)), only_new=sorted(set(f1) - set(f0)),
                     manifest_keys_changed=diff_keys,
                     engine_params_changed={k: (p0.get(k), p1.get(k)) for k in set(p0) | set(p1) if p0.get(k) != p1.get(k)},
                     engine_other_changed=[k for k in m0["engine"] if k != "params" and m0["engine"][k] != m1["engine"].get(k)],
                     parent_manifest_sha256=_sha(f"{P0}/manifest.json"))
    res["R2"]["pass"] = bool(not res["R2"]["different"] and not res["R2"]["only_parent"] and not res["R2"]["only_new"]
                             and set(res["R2"]["engine_params_changed"]) == {"CR_RIM"} and not res["R2"]["engine_other_changed"]
                             and res["R2"]["parent_manifest_sha256"].startswith("4bfc54cb"))
    # R3 / R4
    v0, v1 = json.load(open(f"{P0}/validation.json")), json.load(open(f"{P1}/validation.json"))
    res["R3"] = dict(pass_=v1["pass"], failed=v1["failed_checks"], waived=v1["waived_failures"],
                     parent_failed=v0["failed_checks"], label_rule=v1["checks"]["label_rule"]["agreement"])
    res["R3"]["pass"] = bool(v1["pass"] and v1["failed_checks"] == v0["failed_checks"] and res["R3"]["label_rule"] == 1.0)
    V0, V1 = v0["checks"]["volumes_mc"]["volume_um3"], v1["checks"]["volumes_mc"]["volume_um3"]
    dv = {k: round(V1[k] - V0[k], 4) for k in V0}
    res["R4"] = dict(parent=V0, new=V1, delta=dv, ci95=v1["checks"]["volumes_mc"]["ci95_um3"])
    others = {k: d for k, d in dv.items() if k not in ("mito_matrix", "mito_ims") and d != 0}
    res["R4"]["pass"] = bool(not others and dv["mito_matrix"] >= 0 and abs(dv["mito_matrix"] + dv["mito_ims"]) < 1e-3)
    # R5-R8 paired points
    c0, c1 = _load(a.parent), _load(a.package)
    e = c0.engine
    lam = Lamellae(c0)
    vol_cells = len(lam.cells) * lam.h ** 3
    rng = np.random.default_rng(1)
    tot = dict(n=0, changed=0, ims_to_matrix=0, other_change=0, outside_zone=0, mem_other_diff=0, im_diff_outside=0,
               shell0=0, shell1=0, band=0)
    zone_max = 0.0
    for c in range(0, a.n, 2_000_000):
        m = min(2_000_000, a.n - c)
        X = lam.sample(rng, m)
        L0, L1 = c0.labels(X), c1.labels(X)
        M0, M1 = c0.mems64(X), c1.mems64(X)
        k, u, eb, w = lam.coords(X, M0[:, e.M_OM])
        zone = (np.abs(u) <= lam.W + 0.5) & (eb >= -lam.W - 0.5) & (eb <= 0.5)
        ch = L0 != L1
        tot["n"] += m; tot["changed"] += int(ch.sum())
        tot["ims_to_matrix"] += int((ch & (L0 == e.IMS) & (L1 == e.MATRIX)).sum())
        tot["other_change"] += int((ch & ~((L0 == e.IMS) & (L1 == e.MATRIX))).sum())
        tot["outside_zone"] += int((ch & ~zone).sum())
        oth = [j for j in range(M0.shape[1]) if j != e.M_IM]
        tot["mem_other_diff"] += int((M0[:, oth] != M1[:, oth]).any(1).sum())
        tot["im_diff_outside"] += int(((M0[:, e.M_IM] != M1[:, e.M_IM]) & ~zone).sum())
        tot["shell0"] += int((np.abs(M0[:, e.M_IM]) < 1.0).sum()); tot["shell1"] += int((np.abs(M1[:, e.M_IM]) < 1.0).sum())
        tot["band"] += int(((np.abs(u) < lam.W) & (eb > -lam.W) & (eb < 0) & (w > 0)).sum())
        if ch.any():
            # distance of changed points to the rim line (u = 0, e = -W): <= W*sqrt(2)
            zone_max = max(zone_max, float(np.sqrt(u[ch] ** 2 + (eb[ch] + lam.W) ** 2).max()))
        log(f"paired points {tot['n']:,}: changed {tot['changed']:,}")
    f = vol_cells / tot["n"]
    L_rim = tot["band"] * f / (2 * lam.W ** 2)
    res["R5"] = dict(tot, max_dist_changed_to_rim_axis_nm=round(zone_max, 3), rim_length_um=round(L_rim / 1e3, 1))
    # uniform points over the cell box
    rng2 = np.random.default_rng(2); R = c0.G[0]
    X = rng2.uniform(-R, R, (2_000_000, 3))
    L0, L1 = c0.labels(X), c1.labels(X)
    M0, M1 = c0.mems64(X), c1.mems64(X)
    chb = L0 != L1
    oth = [j for j in range(M0.shape[1]) if j != e.M_IM]
    res["R5"]["box"] = dict(n=2_000_000, changed=int(chb.sum()),
                            ims_to_matrix=int((chb & (L0 == e.IMS) & (L1 == e.MATRIX)).sum()),
                            mem_other_diff=int((M0[:, oth] != M1[:, oth]).any(1).sum()))
    if chb.any():
        k, u, eb, w = lam.coords(X[chb], M0[chb, e.M_OM])
        res["R5"]["box"]["outside_zone"] = int((~((np.abs(u) <= lam.W + 0.5) & (eb >= -lam.W - 0.5) & (eb <= 0.5))).sum())
    else:
        res["R5"]["box"]["outside_zone"] = 0
    res["R5"]["pass"] = bool(tot["other_change"] == 0 and tot["outside_zone"] == 0 and res["R5"]["box"]["outside_zone"] == 0
                             and res["R5"]["box"]["changed"] == res["R5"]["box"]["ims_to_matrix"])
    dV = tot["ims_to_matrix"] * f / 1e9
    exp = L_rim * (2 - np.pi / 2) * lam.W ** 2 / 1e9
    res["R6"] = dict(delta_um3=round(dV, 4), expected_um3=round(float(exp), 4), window=[0.10, 0.25],
                     pass_=bool(0.10 <= dV <= 0.25))
    res["R6"]["pass"] = res["R6"].pop("pass_")
    res["R7"] = dict(other_membranes_diff_points=tot["mem_other_diff"] + res["R5"]["box"]["mem_other_diff"],
                     im_diff_outside_zone=tot["im_diff_outside"])
    res["R7"]["pass"] = bool(res["R7"]["other_membranes_diff_points"] == 0 and tot["im_diff_outside"] == 0)
    A0, A1 = tot["shell0"] * f / 2.0 / 1e6, tot["shell1"] * f / 2.0 / 1e6
    res["R8"] = dict(IM_area_parent_um2=round(A0, 1), IM_area_new_um2=round(A1, 1), delta_um2=round(A1 - A0, 1),
                     expected_um2=round(float(-L_rim * (4 - np.pi) * lam.W / 1e6), 1), window=[-50, -20])
    res["R8"]["pass"] = bool(-50 <= A1 - A0 <= -20)
    return res


def atp(a):
    from wholecell import config, registry
    from wholecell.config import load_site
    from wholecell.preplaced import builders  # noqa: F401
    from wholecell.preplaced.builders import atp_synthase_dimers as AD
    from wholecell.preplaced.context import context
    cfg = config.load("jurkat"); registry.load_plugins(cfg)
    cfg = copy.deepcopy(cfg); cfg["geometry"]["package"] = f"{load_site()['geom_root']}/jurkat/{a.package}"
    cfg["geometry"].pop("manifest_sha256", None)
    base = {}
    for m in (cfg.get("preplaced") or {}).get("serial", {}).get("builders", []):
        if m.get("name") == "atp_synthase":
            base = dict(m["params"])
    if not base:
        base = dict(rims="lamellae", monomer="CPX-6151", seed="preplaced.atp_synthase")
    base.update(json.loads(a.params or "{}"))
    res = {"package": a.package, "params": base}
    ctx = context(cfg, a.bdir, dict(libs=[base.get("monomer", "CPX-6151")]))
    if not a.skip_cell:
        t = time.time()
        registry.get("preplaced", "core.atp_synthase_dimers")(ctx, dict(base))
        st = ctx.state["atp_synthase"]
        res["R9"] = dict(st, seconds=round(time.time() - t))
        res["R9"]["pass"] = bool(st["dimers"] > 0 and st["monomers"] <= st["pool"] / 2)
        log("R9", res["R9"])
    # R10 saturated sub-region
    lam = Lamellae(ctx.geom.cell)
    e = ctx.geom.cell.engine
    k_of = lam.ids[tuple(lam.cells.T)]
    comps, cnt = np.unique(k_of, return_counts=True)
    order = comps[np.argsort(cnt)]
    pick = order[len(order) // 2: len(order) // 2 + a.comps]           # mid-size mitochondria
    cells = lam.cells[np.isin(k_of, pick)]
    rng = np.random.default_rng(3)
    # rim length of the picked lamellae (band method)
    nb_ = 4_000_000
    X = lam.sample(rng, nb_, cells)
    om = ctx.geom.cell.mems64(X)[:, e.M_OM]
    k, u, eb, w = lam.coords(X, om)
    band = ((np.abs(u) < lam.W) & (eb > -lam.W) & (eb < 0) & (w > 0)).sum()
    L_sub = band * len(cells) * lam.h ** 3 / nb_ / (2 * lam.W ** 2)
    sp = float(base.get("spacing_nm", AD.SPACING))
    # rim points: same algorithm as AD.lamella_rim_sites, matrix samples drawn in the picked id cells only
    frac = float((ctx.geom.cell.labels(X[:200000]) == ctx.K.lab("mito_matrix")).mean())
    n_samp = int(a.sat_per_site * L_sub / sp / max(frac, 1e-3) * 1.3)
    orig = AD.lamella_rim_sites

    def sub_sites(ctx_, rng_, n_samples, rho, log=print):
        G = ctx_.cell.G; c3 = ctx_.K.c3
        o = int(G[c3.G_CRISTAE]); nl = int(G[o]); h = lam.h; dims = lam.dims; lo = lam.lo
        gap = float(G[c3.S_CR_GAP]); a_im, off = AD._im_row(ctx_.K, G)
        from wholecell.preplaced import kern
        bnd = kern.factory("atp_synthase.bnd", ctx_.K, AD._make_rims)
        args = ctx_.cell.args(); out = np.zeros(ctx_.K.nmem)
        Xs = lam.sample(rng_, n_samples, cells)
        Xs = Xs[ctx_.cell.label(Xs) == ctx_.K.lab("mito_matrix")]
        sites = []
        b_at = lambda p: bnd(p[0], p[1], p[2], *args, out, a_im, off)  # noqa: E731
        for p in Xs:
            ix = np.floor((p - lo) / h).astype(int)
            kk = int(G[o + 8 + 12 * nl + (ix[0] * dims[1] + ix[1]) * dims[2] + ix[2]])
            if kk < 0 or G[o + 8 + 12 * kk + 11] == 0.0:
                continue
            q = o + 8 + 12 * kk; n = G[q:q + 3]; S = G[q + 4]
            t = n @ p - G[q + 3]; uu = t - S * np.floor(t / S + 0.5); x = p - uu * n
            target = -gap - rho; ok = False
            for _ in range(12):
                bx = b_at(x)
                g = np.array([b_at(x + ee) - b_at(x - ee) for ee in np.eye(3)]) / 2.0
                g -= (g @ n) * n; gn = np.linalg.norm(g)
                if gn < 1e-3:
                    break
                if abs(bx - target) < 0.2:
                    ok = True; break
                x = x - (bx - target) / (gn * gn) * g
            if ok:
                sites.append(np.r_[x, n, g / gn, kk])
        return np.array(sites).reshape(-1, 10)

    AD.lamella_rim_sites = sub_sites
    try:
        t = time.time()
        r = AD.place_lamellae(ctx, 10 ** 9, ctx.seed(base.get("seed", "preplaced.atp_synthase")),
                              dict(base, samples=n_samp), log=log)
    finally:
        AD.lamella_rim_sites = orig
    nd = r["info"]["dimers"]
    dens = nd / (L_sub / 1e3)
    res["R10"] = dict(mitochondria=[int(x) for x in pick], id_cells=int(len(cells)), rim_length_um=round(L_sub / 1e3, 2),
                      samples=n_samp, rim_points=r["info"]["rim_points"], rim_points_per_spacing=round(r["info"]["rim_points"] * sp / L_sub, 2),
                      sites_tried=r["info"]["sites_tried"], dimers=nd, dimers_per_um=round(dens, 1),
                      full_row_per_um=round(1000 / sp, 1), literature_per_um=[50, 83], window=[38, 83],
                      seconds=round(time.time() - t))
    res["R10"]["pass"] = bool(38 <= dens <= 83)
    if "R9" in res and a.rim_length_um:
        res["R10"]["cell_wide_per_um"] = round(res["R9"]["dimers"] / a.rim_length_um, 2)
        res["R10"]["cell_wide_coverage_of_saturated"] = round(res["R9"]["dimers"] / a.rim_length_um / dens, 4)
    log("R10", res["R10"])
    return res


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("geom"); g.add_argument("--parent", default="em_cell1_v8p"); g.add_argument("--package", default="em_cell1_v8r")
    g.add_argument("--n", type=int, default=20_000_000); g.add_argument("--out", required=True)
    t = sp.add_parser("atp"); t.add_argument("--package", default="em_cell1_v8r"); t.add_argument("--bdir", required=True)
    t.add_argument("--params", default=None); t.add_argument("--comps", type=int, default=3)
    t.add_argument("--sat-per-site", type=float, default=5.0); t.add_argument("--rim-length-um", type=float, default=None)
    t.add_argument("--skip-cell", action="store_true"); t.add_argument("--out", required=True)
    a = ap.parse_args()
    t0 = time.time()
    res = geom(a) if a.cmd == "geom" else atp(a)
    res["seconds"] = round(time.time() - t0)
    json.dump(res, open(a.out, "w"), indent=1, default=str)
    log(json.dumps({k: v.get("pass") for k, v in res.items() if isinstance(v, dict) and "pass" in v}))


if __name__ == "__main__":
    main()
