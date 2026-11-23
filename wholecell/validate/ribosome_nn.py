# origin: jurkat/validate/ribosome_nn.py (legacy workspace, 2026-11-15)
# Changed: VAL_GEOM/VAL_POLY defect fixed: the frame centre comes from the build's geometry and the polysome pre-pass
# from the build's own inputs (validate.inputs.polysomes / preplaced manifest), never from default legacy files. The EM
# reference table, its class filter and voxel size are config (validate.ribosome_nn) and a declared val_* frozen input.
# Added: a row with the packed 80S instances of the tiles (library validate.ribosome_nn.packed_lib) merged with the
# pre-pass (criteria V7: polysome pre-pass + packed 80S). Thinning seed: legacy 0 or a wholecell.seeds stream.
"""Ribosome spatial statistics: model vs EM.

EM: ribosome instances of the cell (classes per config), model frame nm = voxel index x voxel size - centre_nm.
Model: the polysome pre-pass (+ packed 80S). EM detection recall is ~0.5, so the model is thinned at random to the EM
count before comparing nearest-neighbour (NN) distances. Output: <build>/validate/ribosome_nn.json
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from . import common

R_BINS = [20, 25, 30, 35, 40, 50, 60, 80, 100]


def stats(P, name):
    t = cKDTree(P); d, _ = t.query(P, k=2); nn = d[:, 1]
    cnt = {f"n_within_{r}": round(float(np.mean([len(x) - 1 for x in t.query_ball_point(P[:200000], r)])), 3) for r in (30, 50, 100)}
    h = {f"frac_nn_lt_{r}": round(float((nn < r).mean()), 4) for r in R_BINS}
    return dict(set=name, n=int(len(P)), nn_median=round(float(np.median(nn)), 2), nn_p10=round(float(np.percentile(nn, 10)), 2),
                nn_p90=round(float(np.percentile(nn, 90)), 2), **h, **cnt)


def packed_positions(inp, lib_name):
    L = json.load(open(inp["libs"]))["libs"]
    if lib_name not in L:
        return np.zeros((0, 3))
    li = L.index(lib_name); out = []
    for f in inp.tile_files():
        z = np.load(f)
        if not (z["job_lib"] == li).any():
            continue
        t = z["types"]; k = t == li
        if k.any():
            out.append(z["pos"][k].astype(np.float64))
    return np.concatenate(out) if out else np.zeros((0, 3))


def run(cfg, bdir, inp: common.Inputs, geo_centre, out_path, overwrite=False):
    from .. import seeds
    rc = (cfg.get("validate", {}) or {}).get("ribosome_nn", {}) or {}
    em_path = inp.lookup_input(rc.get("em_input", "val_ribo_em")) or (cfg.get("frozen_inputs", {}).get(rc.get("em_input", "val_ribo_em"), {}) or {}).get("path")
    E = pd.read_parquet(em_path)
    for col, val in (rc.get("em_filter_eq") or {}).items():
        E = E[E[col] == val]
    for col, vals in (rc.get("em_filter_not_in") or {}).items():
        E = E[~E[col].isin(vals)]
    PE = E[["x", "y", "z"]].values * np.array(rc.get("em_voxel_nm", [4.0, 4.0, 3.44])) - geo_centre
    poly = inp.lookup_input("polysomes")
    z = np.load(poly); PM = z["pos"].astype(np.float64); org = z["origin"]
    legacy = bool(rc.get("legacy_seeds", False))
    rng = np.random.default_rng(0 if legacy else seeds.seed_int(cfg, "validate.ribosome_nn"))
    thin = PM[rng.random(len(PM)) < len(PE) / len(PM)]
    rows = [stats(PE, rc.get("em_label", "EM")), stats(PM, "model (all)"), stats(thin, "model thinned to EM count"),
            stats(PM[org == org.max()], f"model origin={org.max()}"), stats(PM[org == org.min()], f"model origin={org.min()}")]
    jit = thin + rng.normal(0, 120, thin.shape)
    rows.append(stats(jit, "thinned model, 120 nm jitter (no-clustering reference)"))
    extra = {}
    if rc.get("packed_lib"):
        PK = packed_positions(inp, rc["packed_lib"])
        if len(PK):
            ALL = np.concatenate([PM, PK])
            rng2 = np.random.default_rng(seeds.seed_int(cfg, "validate.ribosome_nn.packed"))
            thin2 = ALL[rng2.random(len(ALL)) < len(PE) / len(ALL)]
            rows.append(stats(ALL, f"model pre-pass + packed {rc['packed_lib']}"))
            rows.append(stats(thin2, f"model pre-pass + packed {rc['packed_lib']}, thinned to EM count"))
            extra["n_packed_80S"] = int(len(PK))
    out = dict(rows=rows, origin_codes=np.unique(org).tolist(), note="EM recall ~0.5: compare the EM row with 'thinned' row",
               inputs=dict(polysomes=poly, polysomes_md5=common.md5(poly), em=em_path, centre_nm=list(map(float, geo_centre))), **extra)
    common.write_json(out_path, out, overwrite)
    return out
