"""Adapter: yeast engine pack outputs -> the core pack layout (<build>/pack/), so the core exporter (wholecell.export,
tile_format_v2 / YCT2) and core validation readers work on a yeast build unchanged.

Runs in the engine environment (vendored yeast modules importable):
    python plugins/yeast/pack_view.py --bdir <build> --config-json <build>/config.resolved.json

Writes (nothing is copied; tiles are links to the engine's tiles_raw):
  pack/tiles.parquet      tile (row), ti, tj, tk (200 nm grid, origin 0) + the 27 yeast region measures; extra rows
                          for tiles that hold only pre-placed instances (prepass_only = True)
  pack/tiles_meta.json    TILE 200, lo0 [0, 0, 0]
  pack/tiles/<row>.npz -> engine work/yc3d/tiles_raw/t_..npz (types index lib_ids, pos model nm, quat w,x,y,z, job table)
  pack/libs.json          libs = lib_ids order (library type names)
  pack/preplaced_{types,pos,quat,off}.npy, preplaced_libs.json, preplaced_groups.json
                          pre-placed instances of prepass.npz (wc3d.pre_instances order) bucketed by tile; one
                          preplaced_libs entry per (structure group, library), group = config export.preplaced_group_regions
  pack/summary.json       instance counts (packed, pre-placed), source files
Quaternions stay w-first (tile_format_v2 §2); the exporter canonicalises the sign (w >= 0).
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bdir", required=True)
    ap.add_argument("--config-json", required=True)
    a = ap.parse_args()
    cfg = json.load(open(a.config_json))
    code = cfg["pipeline_engine"]["code"]
    sys.path.insert(0, os.path.join(code, "pipeline"))
    import config as C
    import wc3d
    y = C.W3D
    out = os.path.join(a.bdir, "pack"); os.makedirs(os.path.join(out, "tiles"), exist_ok=True)
    TILE = float(cfg["pack"].get("tile_nm", 200.0))
    meas = pd.read_parquet(f"{y}/tile_measure.parquet")
    names = list(meas.index)
    ijk = np.array([[int(x) for x in re.findall(r"[+-]\d+", n)] for n in names], np.int64)
    key = {tuple(t): i for i, t in enumerate(ijk.tolist())}
    lib_ids = [str(x) for x in np.load(f"{y}/lib_ids.npy", allow_pickle=True)]
    # pre-placed instances, bucketed by the tile of their centre
    pre = dict(np.load(f"{y}/prepass.npz", allow_pickle=False))
    gmap = cfg.get("export", {}).get("preplaced_group_regions", {}) or {}
    parts = [(k, li) for k, li in wc3d.PRE_GROUPS if f"{k}_pos" in pre]
    gt, gp, gq = wc3d.pre_instances(pre)
    grp = np.concatenate([np.full(len(pre[f"{k}_pos"]), gi, np.int64) for gi, (k, _) in enumerate(parts)])
    gnames = [gmap.get(k, k) for k, _ in parts]
    pair = np.c_[grp, gt.astype(np.int64)]
    uniq, inv = np.unique(pair, axis=0, return_inverse=True)
    inv = inv.reshape(-1)
    PL = [lib_ids[int(t)] for _, t in uniq]
    PG = [gnames[int(g)] for g, _ in uniq]
    tk = np.floor(gp.astype(np.float64) / TILE).astype(np.int64)
    extra = []
    tile_of = np.empty(len(gp), np.int64)
    for n, t in enumerate(map(tuple, tk.tolist())):
        i = key.get(t)
        if i is None:
            i = len(names) + len(extra); key[t] = i; extra.append(t)
        tile_of[n] = i
    nt = len(names) + len(extra)
    o = np.argsort(tile_of, kind="stable")
    off = np.zeros(nt + 1, np.int64); np.add.at(off, tile_of + 1, 1); off = np.cumsum(off)
    np.save(os.path.join(out, "preplaced_types.npy"), inv[o].astype(np.int32))
    np.save(os.path.join(out, "preplaced_pos.npy"), gp[o].astype(np.float32))
    np.save(os.path.join(out, "preplaced_quat.npy"), gq[o].astype(np.float32))
    np.save(os.path.join(out, "preplaced_off.npy"), off)
    json.dump(PL, open(os.path.join(out, "preplaced_libs.json"), "w"))
    json.dump(PG, open(os.path.join(out, "preplaced_groups.json"), "w"))
    # tile table + links
    allijk = np.r_[ijk, np.array(extra, np.int64).reshape(-1, 3)]
    T = pd.DataFrame(dict(tile=np.arange(nt), ti=allijk[:, 0], tj=allijk[:, 1], tk=allijk[:, 2],
                          name=names + ["t_%+03d_%+03d_%+03d" % t for t in extra], prepass_only=[False] * len(names) + [True] * len(extra)))
    M = meas.reset_index(drop=True)
    for c in M.columns:
        T[c] = np.r_[M[c].values.astype(np.float64), np.zeros(len(extra))]
    T.to_parquet(os.path.join(out, "tiles.parquet"))
    n_packed = 0
    for i, n in enumerate(names):
        src = os.path.join(C.RAW, f"{n}.npz")
        dst = os.path.join(out, "tiles", f"{i}.npz")
        if os.path.exists(src):
            if not os.path.lexists(dst):
                os.symlink(src, dst)
            with np.load(src) as z:
                n_packed += len(z["types"])
    json.dump(dict(TILE=TILE, lo0=[0.0, 0.0, 0.0], source="yeast tiles3d (200 nm, tile key floor(p/200))", n_tiles=len(names),
                   n_prepass_only_tiles=len(extra)), open(os.path.join(out, "tiles_meta.json"), "w"), indent=1)
    json.dump(dict(libs=lib_ids, mass_da=[], source=f"{y}/lib_ids.npy"), open(os.path.join(out, "libs.json"), "w"))
    S = dict(instances=int(n_packed + len(gp)), packed=int(n_packed), prepass=int(len(gp)), tiles=len(names),
             prepass_entries=len(PL), prepass_groups=sorted(set(PG)), engine_work=C.WORK)
    json.dump(S, open(os.path.join(out, "summary.json"), "w"), indent=1)
    print(json.dumps(S))


if __name__ == "__main__":
    main()
