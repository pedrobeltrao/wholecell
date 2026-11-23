#!/usr/bin/env python
"""A7 yeast side: the yeast pack through the core packer (pack.preset legacy_x5) against the yeast port build s0, X5
criteria P3 / P5 / P6 at level T (docs/reference/history/yeast_x5_reproduction_criteria.md) with the 2026-11-17
addendum (seed-spread intervals when the port builds s1 / s2 are given), plus clash / containment / exclusion checks.
Run in the core environment.

    <core python> scripts/pack_a7_yeast_check.py --build <core build> --s0 <s0 build> [--s1 B --s2 B]
                  [--tiles <npy>|all] [--sample 300] [--nproc 3] --out <dir>

  --tiles npy   subset mode: only these tile ids (rows of tiles.parquet); per-tile and per-(type, region) comparisons
                over the subset (requests differ tile by tile: the core split draws differently, A2)
  --tiles all   whole cell (default): P3 (requested per (type, region), thinning factors), P5 (placed per (type,
                region)), P6 (placed / requested, instances incl. pre-placed), P4-style per-tile placed counts (report)
  --ref-totals  reference "s0 + PORT PATCH 4" (env tolerance addendum 2026-11-17): pairs whose pool the patch changed
                are compared with the s0 (s1, s2) value x E_ref / E_s0 (expected totals after the pre-pass and shell cap)
Clash check (core definition, wholecell/validate/placement.py:_clash): sphere pairs of different instances overlapping
by > 0.5 nm, with the packing trees (the kernel's own collision shapes); pairs with at least one packed instance of the
tile, against packed instances of the tile and its 26 neighbours, the pre-placed instances and the RNA / DNA beads of
those tiles. Containment: packed centres in an allowed label of their region (volume-like regions), bands within
d1 + 4 nm of their membrane, membrane proteins |SDF| <= 0.5 nm; ribosome centres outside the endocytic exclusion zones.
Clash / containment on the subset, or on --sample seeded tiles of the whole cell; for both builds.
Writes <out>/verdict.json (+ per_tile.tsv, pairs.tsv.gz).
"""
import argparse
import json
import os
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell.build import ref_path  # noqa: E402  (reads builds made before the A0c rename in their old names)

ap = argparse.ArgumentParser()
ap.add_argument("--build", required=True)
ap.add_argument("--cell", default=None, help="config of the core build (default: <build>/config.yaml, else yeast)")
ap.add_argument("--s0", required=True)
ap.add_argument("--s1"); ap.add_argument("--s2")
ap.add_argument("--tiles", default="all")
ap.add_argument("--ref-totals", default=None, help="expected totals of the reference 's0 + PORT PATCH 4' (parquet lib_id, "
                "region, expected; e.g. <amphiphysin_fix>/ref_yeastenv_patched/preplaced/serial/"
                "expected_totals_after_shell_capacity.parquet): s0 / s1 / s2 per-(type, region) values are scaled by "
                "E_ref / E_s0 where the pre-pass patch changed the pool")
ap.add_argument("--sample", type=int, default=300)
ap.add_argument("--nproc", type=int, default=3)
ap.add_argument("--out", required=True)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
T0 = time.time()


def tol_count(n, rel):                                     # frozen count tolerance max(3 sqrt n, rel n)
    n = np.abs(np.asarray(n, np.float64))
    return np.maximum(3 * np.sqrt(n), rel * n)


def seed_interval(core, v0, vs, tol):
    """addendum 2026-11-17 s4: I = [min v - m, max v + m], m = max(0.5 (max v - min v), tol); pass if core in I"""
    V = np.stack([v0] + vs) if vs else v0[None]
    lo, hi = V.min(0), V.max(0)
    m = np.maximum(0.5 * (hi - lo), tol)
    return (core >= lo - m) & (core <= hi + m), lo - m, hi + m


def raw_dir(b):
    return os.path.join(b, "engine", "work", "yc3d")


# ------------------------------------------------------------------ inputs
M = pd.read_parquet(os.path.join(a.build, "pack", "tiles.parquet"))
nt = len(M)
LJ = json.load(open(os.path.join(a.build, "pack", "libs.json")))
L = LJ["libs"]
ids0 = [str(x) for x in np.load(os.path.join(raw_dir(a.s0), "lib_ids.npy"))]
assert L == ids0, "core libs.json != s0 lib_ids"
seeds_b = [b for b in (a.s1, a.s2) if b]
for b in seeds_b:
    assert [str(x) for x in np.load(os.path.join(raw_dir(b), "lib_ids.npy"))] == ids0, f"{b}: other library order"
names = list(json.load(open(os.path.join(a.build, "pack", "tiles_meta.json")))["columns"])
s0_view = pd.read_parquet(os.path.join(a.s0, "pack", "tiles.parquet"))
s0_view = s0_view[~s0_view.prepass_only].reset_index(drop=True)
off = np.round(-np.asarray(json.load(open(os.path.join(a.build, "pack", "tiles_meta.json")))["lo0"]) / 200.0).astype(int)
assert len(s0_view) == nt and np.array_equal(s0_view[["ti", "tj", "tk"]].values + off, M[["ti", "tj", "tk"]].values)
s0_names = s0_view.name.tolist()
if a.tiles == "all":
    tiles = np.arange(nt); mode = "whole_cell"
else:
    tiles = np.load(a.tiles).astype(np.int64); mode = "subset"


def tile_core(t):
    f = os.path.join(a.build, "pack", "tiles", f"{t}.npz")
    return np.load(f) if os.path.exists(f) else None


def tile_ref(b, t):
    f = os.path.join(raw_dir(b), "tiles_raw", s0_names[t] + ".npz")
    return np.load(f) if os.path.exists(f) else None


def job_tables(loader, tl):
    """per (lib, region): requested, placed over tiles tl; per tile requested, placed, instances"""
    req, plc, rows = {}, {}, []
    for t in tl:
        z = loader(t)
        if z is None:
            rows.append((t, 0, 0, 0)); continue
        jl, jr, jn, jp = (z[k].astype(np.int64) for k in ("job_lib", "job_region", "job_n", "job_placed"))
        k = jl * 64 + jr
        for kk, n_, p_ in zip(k, jn, jp):
            req[kk] = req.get(kk, 0) + int(n_); plc[kk] = plc.get(kk, 0) + int(p_)
        rows.append((t, int(jn.sum()), int(jp.sum()), len(z["types"])))
    return req, plc, pd.DataFrame(rows, columns=["tile", "requested", "placed", "instances"]).set_index("tile")


def series(d):
    return pd.Series(d, dtype=np.float64)


out = dict(mode=mode, build=a.build, s0=a.s0, seeds=seeds_b, n_tiles=int(len(tiles)), criteria={})
rq_c, pl_c, tc = job_tables(tile_core, tiles)
rq_0, pl_0, t0 = job_tables(lambda t: tile_ref(a.s0, t), tiles)
ref_seeds = [job_tables(lambda t, b=b: tile_ref(b, t), tiles) for b in seeds_b] if mode == "whole_cell" else []
print(f"job tables read ({time.time() - T0:.0f}s)", flush=True)
# ------------------------------------------------------------------ per tile (P4-style; report)
pt = t0.add_suffix("_s0").join(tc.add_suffix("_core"))
tl_ = tol_count(pt.placed_s0, 0.01)
pt["placed_ok"] = (pt.placed_core - pt.placed_s0).abs() <= tl_
pt["req_ok"] = (pt.requested_core - pt.requested_s0).abs() <= tol_count(pt.requested_s0, 0.01)
lg = os.path.join(a.build, "pack", "pack_log.tsv")
if os.path.exists(lg):
    plog = pd.read_csv(lg, sep="\t", header=None, names=["phase", "tile", "req", "placed", "overflow", "sec"])
    pt["sec_core"] = plog.groupby("tile").sec.sum().reindex(pt.index).fillna(0.0)
    pt["overflow_core"] = plog.groupby("tile").overflow.max().reindex(pt.index).fillna(0).astype(int)
ps = pd.read_parquet(os.path.join(raw_dir(a.s0), "pack_summary.parquet")).set_index("ti")
pt["sec_s0"] = ps.sec.reindex(pt.index).values
pt.to_csv(os.path.join(a.out, "per_tile.tsv"), sep="\t")
out["criteria"]["P4_per_tile_placed"] = dict(
    cls="report", tiles=int(len(pt)), within_tol=int(pt.placed_ok.sum()), share=round(float(pt.placed_ok.mean()), 4),
    requested_within_tol=int(pt.req_ok.sum()), passes=bool(pt.placed_ok.mean() >= 0.99),
    tol="max(3 sqrt n, 1 % n) of s0 placed", worst=pt.assign(d=(pt.placed_core - pt.placed_s0)).sort_values("d").iloc[[0, -1]][
        ["requested_s0", "requested_core", "placed_s0", "placed_core"]].reset_index().to_dict("records"))
if "sec_core" in pt:
    out["timing"] = dict(core_cpu_s=round(float(pt.sec_core.sum()), 1), s0_cpu_s=round(float(pt.sec_s0.sum()), 1),
                         ratio=round(float(pt.sec_core.sum() / max(pt.sec_s0.sum(), 1e-9)), 3),
                         per_tile_ratio_median=round(float((pt.sec_core / pt.sec_s0.clip(lower=1e-3)).median()), 3),
                         note="per-tile CPU seconds summed over the 6 rounds incl. the squeeze (core: pack_log.tsv; s0: "
                              "pack_summary.parquet sec); s0 ran 24 workers on eu-g4, core this machine")


# ------------------------------------------------------------------ P5 / P3 per (type, region)
def pair_frame(dc, d0, dseeds):
    k = sorted(set(dc) | set(d0) | {x for d in dseeds for x in d})
    f = pd.DataFrame(dict(core=[dc.get(x, 0) for x in k], s0=[d0.get(x, 0) for x in k]), index=k)
    for i, d in enumerate(dseeds):
        f[f"s{i + 1}"] = [d.get(x, 0) for x in k]
    f["lib"] = [L[x // 64] for x in k]; f["region"] = [names[x % 64] for x in k]
    return f


ADJ = {}
if a.ref_totals and mode == "whole_cell":
    e0 = pd.read_parquet(os.path.join(raw_dir(a.s0), "expected_totals.parquet")).groupby(["lib_id", "region"]).expected.sum()
    er = pd.read_parquet(a.ref_totals).groupby(["lib_id", "region"]).expected.sum()
    u = e0.index.union(er.index); e0 = e0.reindex(u, fill_value=0.0); er = er.reindex(u, fill_value=0.0)
    li = {l: i for i, l in enumerate(L)}; ri = {n: i for i, n in enumerate(names)}
    for (l, r), x0, xr in zip(u, e0.values, er.values):
        if l in li and r in ri and abs(xr - x0) > 1e-9 * max(abs(x0), 1.0):
            ADJ[li[l] * 64 + ri[r]] = xr / x0 if x0 > 0 else np.nan


def adjust(f):
    """reference 's0 + patch 4': scale the s* columns of pairs whose pool the patch changed"""
    if not ADJ:
        return f
    f = f.copy(); k = [x for x in f.index if x in ADJ]
    sc = [c for c in f.columns if c.startswith("s") and c[1:].isdigit()]
    f[sc] = f[sc].astype(np.float64)        # pandas 3 refuses float values in int64 columns (scaled counts are floats)
    for x in k:
        r = ADJ[x]
        if np.isnan(r):
            f.loc[x, sc] = np.nan
        else:
            f.loc[x, sc] = f.loc[x, sc].astype(np.float64) * r
    f["adjusted"] = [x in ADJ for x in f.index]
    return f[~f[sc].isna().any(axis=1)]


def pair_test(f, rel, small_rule):
    f = adjust(f)
    sc = [c for c in f.columns if c.startswith("s") and c[1:].isdigit() and c != "s0"]
    n = f[["s0", "core"] + sc].abs().max(1)
    tol = tol_count(n, rel)
    ok, lo, hi = seed_interval(f.core.values, f.s0.values, [f[c].values for c in sc], tol)
    f = f.assign(tol=tol, lo=lo, hi=hi, ok=ok)
    if small_rule:                                         # P5: n_s0 >= 10 within tol; >= 99 % of n_s0 < 10 within +- 3
        big = f.s0 >= 10
        small = f[~big]
        ok_small = (small.core - small.s0).abs() <= np.maximum(3, small.hi - small.s0)
        res = dict(pairs=int(len(f)), big=int(big.sum()), big_fail=int((~f.ok[big]).sum()),
                   small=int(len(small)), small_within_3=round(float(ok_small.mean()) if len(small) else 1.0, 4))
        res["passes"] = res["big_fail"] == 0 and res["small_within_3"] >= 0.99
    else:
        res = dict(pairs=int(len(f)), fail=int((~f.ok).sum()), differ=int((f.core != f.s0).sum()))
        res["passes"] = res["fail"] == 0
    if "adjusted" in f:
        res["adjusted_pairs_patch4"] = int(f.adjusted.sum()); res["adjusted_fail"] = int((~f.ok[f.adjusted]).sum())
    bad = f[~f.ok].copy(); bad["d"] = bad.core - bad.s0
    res["worst"] = bad.reindex(bad.d.abs().sort_values(ascending=False).index).head(8)[
        ["lib", "region", "s0", "core", "tol"] + sc].to_dict("records")
    return f, res


fP5, r5 = pair_test(pair_frame(pl_c, pl_0, [s[1] for s in ref_seeds]), 0.01, True)
r5.update(cls="must" if mode == "whole_cell" else "subset", rule="X5 P5 (max(3 sqrt n, 1 % n); +-3 for n < 10)"
          + (" + addendum seed spread s0-s2" if ref_seeds else ""))
out["criteria"]["P5_placed_type_region"] = r5
fP3, r3 = pair_test(pair_frame(rq_c, rq_0, [s[0] for s in ref_seeds]), 0.005, False)
r3.update(cls="must" if mode == "whole_cell" else "subset", rule="X5 P3 requested per (type, region), max(3 sqrt n, 0.5 % n)"
          + (" + addendum seed spread" if ref_seeds else ""))
out["criteria"]["P3_requested_type_region"] = r3
pd.concat([fP5.assign(kind="placed"), fP3.assign(kind="requested")]).to_csv(os.path.join(a.out, "pairs.tsv.gz"), sep="\t")
# thinning factors (P3)
thc = os.path.join(a.build, "pack", "thinning.json")
th0 = json.load(open(os.path.join(raw_dir(a.s0), "thinning.json")))
if mode == "whole_cell":
    fc = json.load(open(thc)).get("factors_exact", {}) if os.path.exists(thc) else LJ.get("thinning", {})
    ths = [json.load(open(os.path.join(raw_dir(b), "thinning.json"))) for b in seeds_b]
    rows = {}
    for r in sorted(set(fc) | set(th0)):
        v0 = np.array([th0.get(r, 1.0)]); vs = [np.array([t.get(r, 1.0)]) for t in ths]
        ok, lo, hi = seed_interval(np.array([fc.get(r, 1.0)]), v0, vs, 0.005)
        rows[r] = dict(core=fc.get(r, 1.0), s0=th0.get(r, 1.0), seeds=[t.get(r, 1.0) for t in ths], interval=[lo[0], hi[0]], ok=bool(ok[0]))
    out["criteria"]["P3_thinning"] = dict(cls="must", factors=rows, passes=all(v["ok"] for v in rows.values()), tol=0.005)
# ------------------------------------------------------------------ P6 totals
Pc, Rc = sum(pl_c.values()), sum(rq_c.values())
P0, R0 = sum(pl_0.values()), sum(rq_0.values())
npre_c = int(np.load(ref_path(a.build, "pack/preplaced_off.npy"))[-1]) if mode == "whole_cell" else 0
npre_0 = json.load(open(os.path.join(a.s0, "pack", "summary.json")))["prepass"] if mode == "whole_cell" else 0
fr_c, fr_0 = Pc / max(Rc, 1), P0 / max(R0, 1)
p6 = dict(cls="must" if mode == "whole_cell" else "subset", placed_core=int(Pc), requested_core=int(Rc), frac_core=round(fr_c, 5),
          placed_s0=int(P0), requested_s0=int(R0), frac_s0=round(fr_0, 5), instances_core=int(Pc + npre_c),
          instances_s0=int(P0 + npre_0))
fr_s = [sum(s[1].values()) / max(sum(s[0].values()), 1) for s in ref_seeds]
ins_s = [sum(s[1].values()) + json.load(open(os.path.join(b, "pack", "summary.json")))["prepass"] for s, b in zip(ref_seeds, seeds_b)]
okf, lof, hif = seed_interval(np.array([fr_c]), np.array([fr_0]), [np.array([x]) for x in fr_s], 0.003)
oki, loi, hii = seed_interval(np.array([float(p6["instances_core"])]), np.array([float(p6["instances_s0"])]),
                              [np.array([float(x)]) for x in ins_s], 0.005 * p6["instances_s0"])
p6.update(frac_interval=[round(float(lof[0]), 5), round(float(hif[0]), 5)], instances_interval=[float(loi[0]), float(hii[0])],
          frac_seeds=[round(x, 5) for x in fr_s], instances_seeds=[int(x) for x in ins_s],
          passes=bool(okf[0] and oki[0]), rule="placed fraction +- 0.003, instances +- 0.5 % (+ addendum seed spread)")
out["criteria"]["P6_totals"] = p6
print(json.dumps({k: v.get("passes") for k, v in out["criteria"].items()}), f"({time.time() - T0:.0f}s)", flush=True)

# ------------------------------------------------------------------ clash / containment / exclusion (core definitions)
from wholecell.config import load  # noqa: E402
from wholecell.pack import geomload, switches  # noqa: E402
from wholecell.pack.regions import Regions  # noqa: E402
from wholecell.validate.placement import _clash, _spheres  # noqa: E402

cfg = load(a.cell or "yeast")
pc = switches.effective(cfg["pack"], "legacy_x5")
g = geomload.load(cfg, a.build)
reg = Regions(pc, g)
tr = np.load(os.path.join(raw_dir(a.s0), "pack_trees_15.npz"), allow_pickle=True)
pos_ = {str(i): k for k, i in enumerate(tr["ids"])}
so = [0]; C_, R_ = [], []
for l in L:
    k = pos_[l]; c = tr["C"][tr["off"][k]:tr["off"][k + 1]]; r = tr["R"][tr["off"][k]:tr["off"][k + 1]]
    C_.append(c); R_.append(r); so.append(so[-1] + len(r))
so = np.array(so, np.int64); STC = np.concatenate(C_).astype(np.float64); STR = np.concatenate(R_).astype(np.float64)


def index(b):
    d = os.path.join(b, "pack")
    PL = json.load(open(ref_path(d, "preplaced_libs.json")))
    lidx = np.array([L.index(x) for x in PL], np.int64)
    mo = os.path.join(d, "mrna_off.npy")
    return dict(pt=np.load(ref_path(d, "preplaced_types.npy"), mmap_mode="r"), pp=np.load(ref_path(d, "preplaced_pos.npy"), mmap_mode="r"),
                pq=np.load(ref_path(d, "preplaced_quat.npy"), mmap_mode="r"), po=np.load(ref_path(d, "preplaced_off.npy")), lidx=lidx,
                mp=np.load(os.path.join(d, "mrna_pts.npy"), mmap_mode="r") if os.path.exists(mo) else None,
                mo=np.load(mo) if os.path.exists(mo) else None,
                mr=json.load(open(os.path.join(d, "mrna_meta.json")))["bead_radius"] if os.path.exists(mo) else 0.5)


IDX_CORE = index(a.build)
IDX_S0 = IDX_CORE if mode == "subset" else index(a.s0)       # whole cell: s0 keeps its own (unpatched) pre-pass
Tij = M[["ti", "tj", "tk"]].values
KEY = {tuple(t): i for i, t in enumerate(Tij)}
lo0 = np.asarray(json.load(open(os.path.join(a.build, "pack", "tiles_meta.json")))["lo0"])
HALO = 40.0


def neigh(t):
    out_ = []
    for d in np.stack(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"), -1).reshape(-1, 3):
        j = KEY.get(tuple(Tij[t] + d))
        if j is not None and j != t:
            out_.append(j)
    return out_


def clash_tile(t, loader, IDX):
    z = loader(t)
    if z is None or not len(z["types"]):
        return dict(n=0, pp=0, po=0, pb=0, depth=0.0)
    lo = Tij[t] * 200.0 + lo0; glo, ghi = lo - HALO, lo + 200.0 + HALO
    ty = z["types"].astype(np.int64); n0 = len(ty)
    parts = [_spheres(ty, z["pos"].astype(np.float64), z["quat"], so, STC, STR, 0)]
    nxt = n0; cats = [np.zeros(n0, np.int8)]
    for j in [t] + neigh(t):
        if j != t:
            zj = loader(j)
            if zj is not None and len(zj["types"]):
                k = np.all((zj["pos"] >= glo - 40) & (zj["pos"] <= ghi + 40), axis=1)
                if k.any():
                    parts.append(_spheres(zj["types"][k].astype(np.int64), zj["pos"][k].astype(np.float64), zj["quat"][k], so, STC, STR, nxt))
                    cats.append(np.ones(int(k.sum()), np.int8)); nxt += int(k.sum())
        p0, p1 = int(IDX["po"][j]), int(IDX["po"][j + 1])
        if p1 > p0:
            pp = np.asarray(IDX["pp"][p0:p1], np.float64); k = np.all((pp >= glo - 70) & (pp <= ghi + 70), axis=1)
            if k.any():
                parts.append(_spheres(IDX["lidx"][np.asarray(IDX["pt"][p0:p1])[k].astype(np.int64)], pp[k],
                                      np.asarray(IDX["pq"][p0:p1])[k], so, STC, STR, nxt))
                cats.append(np.ones(int(k.sum()), np.int8)); nxt += int(k.sum())
        if IDX["mo"] is not None:
            m0, m1 = int(IDX["mo"][j]), int(IDX["mo"][j + 1])
            if m1 > m0:
                mp = np.asarray(IDX["mp"][m0:m1], np.float64); k = np.all((mp >= glo - 1) & (mp <= ghi + 1), axis=1)
                if k.any():
                    parts.append((mp[k], np.full(int(k.sum()), IDX["mr"]), np.arange(nxt, nxt + int(k.sum()))))
                    cats.append(np.full(int(k.sum()), 2, np.int8)); nxt += int(k.sum())
    C = np.concatenate([p[0] for p in parts]); R = np.concatenate([p[1] for p in parts]); I = np.concatenate([p[2] for p in parts])
    cat = np.concatenate(cats)
    # only pairs with a packed instance of this tile: other-other pairs are dropped below
    ia, ib, d, tot = _clash(C, R, I, np.where(cat == 0, 0, 2).astype(np.int8), 6.0, 0.5, 5_000_000)
    own = (ia < n0) | (ib < n0)
    ia, ib, d = ia[own], ib[own], d[own]
    if len(ia):
        key = np.minimum(ia, ib) * (nxt + 1) + np.maximum(ia, ib); u, first = np.unique(key, return_index=True)
        other = np.where(ia[first] < n0, ib[first], ia[first]); co = cat[other]
        return dict(n=n0, pp=int((co == 0).sum()), po=int((co == 1).sum()), pb=int((co == 2).sum()), depth=float(d.max()))
    return dict(n=n0, pp=0, po=0, pb=0, depth=0.0)


def contain_tile(t, loader, EXZ, ribo):
    z = loader(t)
    if z is None or not len(z["types"]):
        return dict(n=0, out=0, exc=0)
    rg = np.repeat(z["job_region"].astype(np.int64), z["job_placed"].astype(np.int64))
    P = z["pos"].astype(np.float64); ty = z["types"].astype(np.int64)
    assert len(rg) == len(P)
    lab = g.cell.labels(P); S = g.cell.mems64(P) if hasattr(g.cell, "mems64") else g.cell.mems(P).astype(np.float64)
    bad = np.zeros(len(P), bool)
    for ri in np.unique(rg):
        k = rg == ri; nm = reg.RNAMES[ri]; r = reg.spec[nm]
        if r["kind"] == "membrane":
            bad[k] = np.abs(S[k, reg.mem_codes[r["membrane"]]]) > 0.5
            continue
        labs = [reg.label_codes[l] for l in r["labels"]]
        bad[k] = ~np.isin(lab[k], labs)
        if r["kind"] == "shell":
            x = reg.restriction(nm)
            bad[k] |= np.abs(S[k, x[0]]) > x[3] + float(pc.get("volume_shell_tol_nm", 4.0))
    exc = 0
    if len(EXZ):
        kr = ribo[ty]
        if kr.any():
            d2 = ((P[kr][:, None, :] - EXZ[None, :, :3]) ** 2).sum(-1)
            exc = int((d2 < EXZ[None, :, 3] ** 2).any(1).sum())
    return dict(n=len(P), out=int(bad.sum()), exc=exc)


if mode == "subset":
    ctiles = tiles
else:                                                    # seeded sample of tiles with packed content
    rng = np.random.default_rng(20261117)
    have = np.array([t for t in range(nt) if pt.placed_core.get(t, 0) > 0])
    ctiles = np.sort(rng.choice(have, size=min(a.sample, len(have)), replace=False))
from plugins.yeast.exclusion import zones  # noqa: E402
from plugins.yeast import geom_engine  # noqa: E402
EXZ = zones(geom_engine.load(g.source).endocytic_sites(), 90.5, 100.0, 30.0)
ribo = np.array([l.startswith(("ribosome_80S", "preribo_40S", "preribo_60S_cytoplasm")) for l in L])
os.environ.setdefault("NUMBA_NUM_THREADS", str(a.nproc))
rows = []
for t in ctiles:
    t = int(t)
    cc = clash_tile(t, tile_core, IDX_CORE); c0 = clash_tile(t, lambda j: tile_ref(a.s0, j), IDX_S0)
    kc = contain_tile(t, tile_core, EXZ, ribo); k0 = contain_tile(t, lambda j: tile_ref(a.s0, j), EXZ, ribo)
    rows.append(dict(tile=t, **{f"{k}_core": v for k, v in cc.items()}, **{f"{k}_s0": v for k, v in c0.items()},
                     **{f"c{k}_core": v for k, v in kc.items()}, **{f"c{k}_s0": v for k, v in k0.items()}))
CT = pd.DataFrame(rows).set_index("tile")
CT.to_csv(os.path.join(a.out, "clash_containment.tsv"), sep="\t")
per_m = lambda x, n: round(1e6 * float(x) / max(float(n), 1.0), 3)
cl = {}
for b in ("core", "s0"):
    n = CT[f"n_{b}"].sum()
    cl[b] = dict(packed=int(n), packed_packed=int(CT[f"pp_{b}"].sum()), packed_other=int(CT[f"po_{b}"].sum()),
                 packed_bead=int(CT[f"pb_{b}"].sum()), clash_rate_per_1e6=per_m(CT[f"pp_{b}"].sum() + CT[f"po_{b}"].sum(), n),
                 max_depth_nm=round(float(CT[f"depth_{b}"].max()), 3), outside=int(CT[f"cout_{b}"].sum()),
                 outside_per_1e6=per_m(CT[f"cout_{b}"].sum(), n), ribosomes_in_exclusion_zones=int(CT[f"cexc_{b}"].sum()))
out["clash_containment"] = dict(tiles=int(len(CT)), definition="packing trees (15 res/sphere), overlap > 0.5 nm; pairs with a "
                                "packed instance of the tile vs packed (tile + neighbours), pre-placed and beads", **cl,
                                passes=bool(cl["core"]["packed_packed"] + cl["core"]["packed_other"] == 0 and cl["core"]["outside"] == 0
                                            and cl["core"]["ribosomes_in_exclusion_zones"] == 0))
must = [v for v in out["criteria"].values() if v.get("cls") == "must"]
if mode == "whole_cell":
    out["verdict"] = "PASS" if all(v["passes"] for v in must) and out["clash_containment"]["passes"] else "FAIL"
else:
    out["verdict"] = "PASS (subset)" if all(v["passes"] for k, v in out["criteria"].items() if k != "P4_per_tile_placed") \
        and out["clash_containment"]["passes"] else "SEE DETAILS (subset)"
out["seconds"] = round(time.time() - T0, 1)
json.dump(out, open(os.path.join(a.out, "verdict.json"), "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
print(json.dumps(out, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))[:6000])
