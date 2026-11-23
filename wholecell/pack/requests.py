"""Tile requests, method `legacy_x5` (pack.requests): the yeast X5 split of requested copies over tiles.

The v7 method (`multinomial`, default) lives in stage.step_counts. `legacy_x5` reproduces yeast stage 3: each
(library type, region) total is rounded stochastically and split over the tiles by a multinomial proportional to the
tile's region measure, times exp(log weight) of the organisation field for types with a network profile and a crowding
correction that keeps each tile's requested macromolecular volume proportional to its region measure. Region ids =
positions in pack.regions; tile ids = rows of tiles.parquet.

Declared inputs (paths relative to the build directory, or absolute), until the core makes them (A5/A6/A7):
  totals         parquet (lib_id, ing_id, region, expected): requested copies per type and region after the pre-pass
                 consumption and shell cap (yeast yc3d/expected_totals.parquet; X5 criterion P1)
  lib_order      npy of library type ids: the lib index of the output (yeast yc3d/lib_ids.npy)
  shape_volumes  npz (ids, off, R): packing sphere trees; summed sphere volume per type (crowding correction)
  thinning       json {region: factor} (pilot-calibrated thinning of the packer; optional)
  organisation   {loadings: npz (ids, tix, E): field loadings per type with a profile, beta, field_l_nm, n_features}
                 (optional). Without it and with pack.organisation on, the loadings of the core organisation step
                 (pack/organisation/network.npz, float64) and the pack.organisation beta / field_l_nm / n_features
                 are used (A7 yeast side).
With pack.shell_cap on, the shell cap (origin: yeast wc3d.shell_capacity) is applied to the totals before the split
(a no-op on totals that are already capped, as the yeast expected_totals are). With pack.thinning on and resplit, the
pack step calls legacy_x5 again with the pilot factors (yeast wc3d.distribute(thin(...)): same seed, thinned totals).
Seed streams: pack.counts (the split), pack.organisation.field (the random fields), pack.counts.field_cap (the split
of the regions listed in pack.organisation.field_cap_regions; freeze fix 2026-11-19, field_cap).
"""
from __future__ import annotations

import json
import os

import numpy as np

KEYS = ("method", "totals", "lib_order", "shape_volumes", "thinning", "organisation")


def _path(bdir, v):
    return v if os.path.isabs(v) else os.path.join(bdir, v)


# origin: yeast pipeline/organisation.py:field_values @ e3b7448 (seed and constants as arguments)
def field_values(points, n_fields, seed, field_l, n_features):
    """unit-variance smooth random fields (correlation length field_l) at points (nm); n x n_fields"""
    rng = np.random.default_rng(seed)
    K = rng.normal(0, 1.0 / field_l, (n_fields, n_features, 3)); ph = rng.uniform(0, 2 * np.pi, (n_fields, n_features))
    return np.sqrt(2.0 / n_features) * np.cos(np.einsum("pc,mfc->pmf", points, K) + ph[None]).sum(2)


# origin: yeast pipeline/organisation.py:tile_logmod @ e3b7448
def tile_logmod(lib_ids, tix, E, centres, beta, seed, field_l, n_features):
    """per library type with a network profile, the log weight over tiles (centres, nm): beta * E_i . F"""
    F = field_values(centres, E.shape[1], seed, field_l, n_features)
    return {lib_ids[i]: beta * (F @ E[k]) for k, i in enumerate(tix)}


# origin: yeast pipeline/wc3d.py:crowding_factors @ e3b7448
def crowding_factors(e, M, reg_idx, logmod, vol, iters=30):
    """Per region, tile factors c_t for the modulated types such that the expected macromolecular volume per tile
    stays proportional to the region's measure in the tile (composition varies between tiles, crowding does not)"""
    out = {}
    for reg, g in e[e.lib_id.isin(logmod)].groupby("region"):
        ri = reg_idx[reg]; W = M[:, ri]; nz = np.nonzero(W)[0]
        if len(nz) == 0:
            continue
        A = np.stack([np.exp(logmod[l][nz] - logmod[l][nz].max()) for l in g.lib_id])         # rows x tiles
        V = (g.expected.values * np.array([vol[l] for l in g.lib_id]))[:, None]
        c = np.ones(len(nz))
        for _ in range(iters):
            Wc = A * (W[nz] * c)[None]
            occ = (V * Wc / Wc.sum(1, keepdims=True)).sum(0)
            ratio = (W[nz] / W[nz].sum()) / (occ / occ.sum())
            c *= ratio
        cf = np.zeros(len(W)); cf[nz] = c / c.mean(); out[ri] = cf
    return out


# freeze fix (2026-11-19, main-session decision): field cap in crowded volume regions (pack.organisation.field_cap_regions)
def field_cap(e, M, reg_idx, logmod, cf, mass, regions, iters=500, tol=1e-6):
    """Per listed region, tile factors (replacing cf[region]) for the modulated types such that no tile is asked for more
    macromolecular MASS than the unmodulated split gives it: U_t <= R_t, R_t = sum_l x_l m_l W_t / sum W (x_l: expected
    copies of type l in the region, m_l its mass, W_t the tile's region measure, U_t the expected mass under the field).

    Why mass: the crowding correction (crowding_factors) evens the summed packing-sphere volume per tile, but sphere
    volume per Da differs between types (overlapping spheres of large complexes: Jurkat nucleolus 12.9 vs 19.7 Da / nm^3),
    so the field still moves mass between tiles (Jurkat nucleolus, B2: expected mass per tile 0.71-1.41 x the
    unmodulated split); tiles asked for more fill only part of it (0.60-0.72) and the region loses mass (0.864 -> 0.807).
    The excess of the capped tiles is redistributed to the tiles below their unmodulated request (not dropped: the cell
    totals and the type's composition field stay; dropping would lower the requested mass the K2 fill is judged against).
    With the totals conserved, U_t <= R_t on every tile forces U_t = R_t, so the fixed point is reached by a multiplicative
    per-tile scaling of the modulated types' weights (iterative proportional fitting from cf, the minimal change in the
    Kullback-Leibler sense); unmodulated types keep the measure split. Types keep their field contrast between tiles of
    equal measure (composition varies, mass per tile does not).
    e: rows (lib_id, region, expected); M: tiles x regions measures; mass: {lib_id: Da}; regions: names.
    Returns ({region index: factor per tile}, {region name: dict(max_ratio_before, max_ratio_after, iterations)})."""
    out = dict(cf); rep = {}
    for reg in regions:
        if reg not in reg_idx:
            continue
        ri = reg_idx[reg]; W = M[:, ri]; nz = np.nonzero(W)[0]
        g = e[e.region == reg]
        if len(nz) == 0 or not len(g):
            continue
        Wn = W[nz]; ref = Wn / Wn.sum()
        x = g.expected.values.astype(np.float64) * np.array([float(mass[l]) for l in g.lib_id])
        mod = np.array([l in logmod for l in g.lib_id], bool)
        R = ref * x.sum()
        base = ref * x[~mod].sum()                                # unmodulated types: the measure split
        if not mod.any():
            continue
        c0 = cf[ri][nz] if ri in cf else np.ones(len(nz))
        A = np.stack([np.exp(logmod[l][nz] - logmod[l][nz].max()) for l in g.lib_id[mod]]) * (Wn * c0)[None]
        xm = x[mod][:, None]
        a = np.ones(len(nz)); before = None; it = 0
        for it in range(1, iters + 1):
            P = A * a[None]
            U = base + (xm * P / P.sum(1, keepdims=True)).sum(0)
            ratio = U / R
            if before is None:
                before = float(ratio.max())
            if ratio.max() <= 1.0 + tol:
                break
            a /= ratio
        P = A * a[None]
        U = base + (xm * P / P.sum(1, keepdims=True)).sum(0)
        f = np.zeros(len(W)); f[nz] = c0 * a; f[nz] /= f[nz].mean()
        out[ri] = f
        rep[reg] = dict(max_ratio_before=round(before, 4), max_ratio_after=round(float((U / R).max()), 6), iterations=it,
                        tiles=int(len(nz)), types_modulated=int(mod.sum()))
    return out, rep


# origin: yeast pipeline/wc3d.py:thin @ e3b7448
def thin(exp, factors):
    """the packing request: expected copies x the region's factor"""
    exp = exp.copy()
    for name, f in factors.items():
        exp.loc[exp.region == name, "expected"] *= f
    return exp


# origin: yeast pipeline/wc3d.py:distribute @ e3b7448 (region names, measure matrix and seed as arguments)
def distribute(exp, M, rnames, lib_index, seed, logmod=None, vol=None, cap=None, report=None, cap_seed=None):
    """Multinomial split of each (lib, region) total over tiles proportional to the tile's region measure (M: tiles x
    regions in region order), times exp(logmod[lib]) per tile if given and the crowding correction for those types.
    cap: (mass {lib_id: Da}, region names) for the field cap (field_cap; freeze fix 2026-11-19), its report goes to
    report["field_cap"] if a dict is given; the tile split of a capped region is drawn from Generator(cap_seed)
    (default: seed + 1) after the uncapped draw from `seed` (discarded), so every other region keeps its draws bitwise.
    Returns long table tile, lib, region, n (cell totals conserved after stochastic rounding)."""
    import pandas as pd
    rng = np.random.default_rng(seed)
    reg_idx = {n: i for i, n in enumerate(rnames)}
    e = exp[exp.region.isin(reg_idx)].groupby(["lib_id", "region"], as_index=False).expected.sum()
    cf = crowding_factors(e, M, reg_idx, logmod, vol) if logmod is not None else {}
    capped = {}
    if logmod is not None and cap is not None and cap[1]:
        cf2, rep = field_cap(e, M, reg_idx, logmod, cf, cap[0], cap[1])
        capped = {reg_idx[r]: cf2[reg_idx[r]] for r in rep}
        if report is not None:
            report["field_cap"] = rep
    rng_cap = np.random.default_rng(seed + 1 if cap_seed is None else cap_seed) if capped else None
    T, Li, Ri, Nn = [], [], [], []
    for r in e.itertuples():
        ri = reg_idx[r.region]
        w = M[:, ri]
        if w.sum() <= 0:
            continue
        N = int(np.floor(r.expected) + (rng.random() < r.expected - np.floor(r.expected)))
        if N == 0:
            continue
        nz = np.nonzero(w)[0]
        w0 = w
        if logmod is not None and r.lib_id in logmod:
            lw = logmod[r.lib_id][nz]
            w = np.zeros_like(w); w[nz] = M[nz, ri] * np.exp(lw - lw.max()) * (cf[ri][nz] if ri in cf else 1.0)
        cnt = rng.multinomial(N, w[nz] / w[nz].sum())
        if ri in capped:          # field cap: own stream; the uncapped draw above keeps stream `seed` of the other regions
            if logmod is not None and r.lib_id in logmod:
                w = np.zeros_like(w0); w[nz] = M[nz, ri] * np.exp(lw - lw.max()) * capped[ri][nz]
            cnt = rng_cap.multinomial(N, w[nz] / w[nz].sum())
        k = cnt > 0
        T.append(nz[k]); Li.append(np.full(k.sum(), lib_index[r.lib_id])); Ri.append(np.full(k.sum(), ri)); Nn.append(cnt[k])
    d = pd.DataFrame(dict(tile=np.concatenate(T).astype(np.int32), lib=np.concatenate(Li).astype(np.int32),
                          region=np.concatenate(Ri).astype(np.int16), n=np.concatenate(Nn).astype(np.int32)))
    return d.sort_values(["tile"]).reset_index(drop=True)


def shape_volumes(path, lib_ids):
    """{library type: summed sphere volume of its packing tree}"""
    z = np.load(path, allow_pickle=True)
    pos = {str(i): k for k, i in enumerate(z["ids"])}
    off, R = z["off"], z["R"]
    miss = [l for l in lib_ids if l not in pos]
    if miss:
        raise ValueError(f"pack.requests.shape_volumes {path}: no tree for {len(miss)} types, e.g. {miss[:5]}")
    return {l: float((4 / 3 * np.pi * R[off[pos[l]]:off[pos[l] + 1]] ** 3).sum()) for l in lib_ids}


def _shell_cap_totals(exp, reg, Mx, vol, log=print):
    """pack.shell_cap on a totals table (lib_id, ing_id, region, expected) with config region names; origin: yeast
    wc3d.shell_capacity @ e3b7448 (each capped shell's excess goes to the bulk region holding its first label, every
    species of the shell by the same factor). Returns (table, report)."""
    import pandas as pd
    from .algorithms import shell_cap
    rid = {n: i for i, n in enumerate(reg.RNAMES)}
    e = exp[exp.region.isin(rid)]
    libs = sorted(set(e.lib_id)); li = {l: k for k, l in enumerate(libs)}
    tot = {}
    for r in e.itertuples():
        k = (li[r.lib_id], rid[r.region]); tot[k] = tot.get(k, 0.0) + float(r.expected)
    vol_lib = np.array([vol[l] for l in libs])
    vol_regions = {rid[n]: labs for n, labs in reg.VOL_REG if reg.spec[n]["kind"] == "volume"}
    region_volume = {i: float(Mx[:, i].sum()) for i, n in enumerate(reg.RNAMES) if reg.spec[n]["kind"] in ("volume", "virtual")
                     or (reg.spec[n]["kind"] == "shell" and reg.spec[n].get("measure", "area") == "volume")}
    shells = []
    for n in reg.RNAMES:
        r = reg.spec[n]
        if r["kind"] != "shell":
            continue
        x = reg.restriction(n) if r.get("measure", "area") == "volume" else None
        shells.append(dict(ri=rid[n], area_nm2=0.0, d0=float(r.get("d0", 0.0)), d1=float(r["d1"]), add_radius=False,
                           pore_limit=int(r.get("pore_limit_nm", 0)), labels=tuple(reg.label_codes[l] for l in r["labels"]),
                           volume_nm3=float(Mx[:, rid[n]].sum()) if x is not None else None))
    shells = [s_ for s_ in shells if s_["volume_nm3"] is not None]          # area shells: core counts step (v7 method)
    tot2, rep = shell_cap(tot, vol_lib, np.zeros(len(libs)), vol_regions, region_volume, shells, log=log)
    if not rep:
        return exp, {}
    ing = e.groupby(["lib_id", "region"]).ing_id.first().to_dict()
    rows = []
    for (l, ri), x in tot2.items():
        rows.append((libs[l], ing.get((libs[l], reg.RNAMES[ri]), libs[l]), reg.RNAMES[ri], x))
    out = pd.concat([exp[~exp.region.isin(rid)], pd.DataFrame(rows, columns=["lib_id", "ing_id", "region", "expected"])],
                    ignore_index=True)
    return out, {reg.RNAMES[ri]: dict(factor=round(f, 4), moved=round(mv, 1), to=reg.RNAMES[t]) for ri, (f, mv, t) in rep.items()}


def legacy_x5(cfg, bdir, p, reg, M, lo0, sw=None, factors=None, log=print):
    """counts table (tile, lib, region, n) + libs.json content for pack.requests method legacy_x5. sw: packer switches
    (organisation loadings from the core step, shell cap); factors: {region name: thinning factor} (pilot; instead of
    the requests.thinning file)"""
    import pandas as pd
    from wholecell import seeds
    rq = p.pc["requests"]
    bad = [k for k in rq if k not in KEYS]
    if bad:
        raise ValueError(f"pack.requests: unknown key(s) {bad} (known: {', '.join(KEYS)})")
    for k in ("totals", "lib_order"):
        if not rq.get(k):
            raise ValueError(f"pack.requests method legacy_x5 needs `{k}` (see wholecell/pack/requests.py)")
    exp = pd.read_parquet(_path(bdir, rq["totals"]))
    lib_ids = [str(x) for x in np.load(_path(bdir, rq["lib_order"]), allow_pickle=True)]
    rnames = reg.RNAMES
    Mx = np.stack([M[reg.tile_column(ri)].values.astype(np.float64) for ri in range(len(rnames))], 1)
    logmod = vol = None
    og = rq.get("organisation")
    if not og and sw is not None and sw.organisation["enabled"]:          # the core organisation step (A7)
        o = sw.organisation
        og = dict(loadings=os.path.join(p.out, "organisation", "network.npz"), beta=o["beta"], field_l_nm=o["field_l_nm"],
                  n_features=o["n_features"])
        if not os.path.exists(og["loadings"]):
            raise FileNotFoundError(f"pack.organisation is on but {og['loadings']} is missing: run `pack --step organisation`")
    if og:
        z = np.load(_path(bdir, og["loadings"]), allow_pickle=True)
        if [str(x) for x in z["ids"]] != lib_ids:
            raise ValueError("pack.requests.organisation.loadings: ids differ from lib_order")
        ijk = M[["ti", "tj", "tk"]].values + np.round(np.asarray(lo0) / p.TILE).astype(np.int64)
        centres = (ijk.astype(np.float64) + 0.5) * p.TILE
        logmod = tile_logmod(lib_ids, z["tix"], np.asarray(z["E"], np.float64), centres, float(og.get("beta", 0.8)),
                             seeds.seed_int(cfg, "pack.organisation.field"), float(og.get("field_l_nm", 500.0)),
                             int(og.get("n_features", 128)))
        if not rq.get("shape_volumes"):
            raise ValueError("pack.requests.organisation needs shape_volumes (crowding correction)")
    if og or (sw is not None and sw.shell_cap["enabled"]):
        if not rq.get("shape_volumes"):
            raise ValueError("pack.requests legacy_x5 with the organisation layer or the shell cap needs shape_volumes")
        vol = shape_volumes(_path(bdir, rq["shape_volumes"]), lib_ids)
    cap = {}
    if sw is not None and sw.shell_cap["enabled"]:
        exp, cap = _shell_cap_totals(exp, reg, Mx, vol, log=log)
        log(f"shell cap (legacy_x5 totals): {len(cap)} shells capped" + (f": {cap}" if cap else " (totals already capped)"))
    if factors is None:
        factors = json.load(open(_path(bdir, rq["thinning"]))) if rq.get("thinning") else {}
    if factors:
        exp = thin(exp, factors)
    seed = seeds.seed_int(cfg, "pack.counts")
    fcap = None; rep = {}
    cap_regions = list(sw.organisation["field_cap_regions"]) if (og and sw is not None) else []
    if cap_regions:                                       # freeze fix 2026-11-19 (field_cap): masses from the pack table
        bad = [r for r in cap_regions if r not in {n for n, _ in reg.VOL_REG}]
        if bad:
            raise ValueError(f"config pack.organisation.field_cap_regions: {bad} are not volume regions")
        pt = pd.read_csv(p.pack_table(), sep="\t", usecols=["ing_id", "mass_da"])
        md = dict(zip(pt.ing_id.astype(str), pt.mass_da.astype(np.float64)))
        sub = exp[exp.region.isin(cap_regions)]
        miss = sorted(set(sub.ing_id.astype(str)) - set(md))
        if miss:
            raise ValueError(f"pack.organisation.field_cap_regions: {len(miss)} ingredients without mass_da in {p.pack_table()}, e.g. {miss[:5]}")
        w = sub.assign(m=sub.ing_id.astype(str).map(md) * sub.expected)
        gm = w.groupby("lib_id")[["m", "expected"]].sum()
        fcap = ((gm.m / gm.expected.where(gm.expected > 0, 1.0)).to_dict(), cap_regions)   # copy-weighted mean mass per type
    d = distribute(exp, Mx, rnames, {k: i for i, k in enumerate(lib_ids)}, seed, logmod=logmod if og else None, vol=vol,
                   cap=fcap, report=rep, cap_seed=seeds.seed_int(cfg, "pack.counts.field_cap") if fcap else None)
    if rep.get("field_cap"):
        log(f"organisation field cap (mass per tile <= the unmodulated split): {rep['field_cap']}")
    info = dict(libs=lib_ids, rnames=list(rnames), requests=dict(rq), thinning=factors, counts_seed=seed,
                organisation_seed=seeds.seed_int(cfg, "pack.organisation.field") if og else None,
                organisation_loadings=og.get("loadings") if og else None, shell_cap=cap,
                **({"field_cap": rep["field_cap"]} if rep.get("field_cap") else {}))
    return d, info
