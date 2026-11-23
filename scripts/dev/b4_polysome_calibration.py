# new in wholecell (B4 addendum, 2026-11-17)
"""Per-cell polysome linker (core.polysomes params.linker_nt_mean) measured from EM ribosome data on a calibration
subset (docs/unification_criteria/B4_preplaced_addendum_2026-11-17.md §1).

  run         the Jurkat serial builders npc, microtubules, polysomes (config/cells/jurkat.yaml; the obstacles the
              polysome walk sees) with params.linker_nt_mean = <value> ("transcripts" = the B4 default: per chain from
              the transcript table) -> <out>.npz (ribosome centres, chain, EM-prior flag) + <out>.json (builder stats)
                python scripts/b4_polysome_calibration.py run --bdir <dir: library, inventory> --linker 40 --out <prefix>
              [--walk '{"candidates": 32, "rounds": 8}']  (params.walk of the walk fix 2026-11-18; '' = legacy walk)
  complement  uniform random cytosol points (seeded), the proxy for the free 80S the packer adds (pool - pre-placed)
                python scripts/b4_polysome_calibration.py complement --bdir <dir> --n 700000 --out <npy>
  score       EM-comparable statistics (J4 NN median of the model thinned to the EM count; J3b 34-nm single-linkage
              share of the model thinned to the de-duplicated EM count; same thinning streams as b4_preplaced_check)
              on the CALIBRATION SUBSET only, against the EM ribosomes in the same subset; picks the value with the
              lowest |NN - EM| / 1 nm + |linked - EM| / 0.05 among those keeping J3a (polysome fraction) in 0.68-0.72
                python scripts/b4_polysome_calibration.py score --runs a.npz,b.npz --complement c.npy --out <json>

Calibration subset (fixed, declared before any scan): model frame (nm, EM coordinates x voxel minus the geometry
centre); 1-um blocks i = floor(x / 1000) (same for j, k); a point is in the subset iff (i + 2 j + 3 k) mod 5 == 0
(one block in five, spread over the whole cell). Statistics are medians / shares over the points inside the subset
(neighbours searched among all points, so no edge effect). The evaluation (frozen J3b / J4) is the whole cell;
its complement (the 4 blocks in 5 never used for the choice) is reported too.
"""
import argparse
import copy
import json
import os
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
from scripts.dev.b4_preplaced_check import EM_N, EM_DEDUP_N, EM_LINKED, EM_NN  # noqa: E402

BLOCK_NM, MOD, KEEP = 1000.0, 5, 0
J3A = (0.68, 0.72)                      # frozen J3a window (EM polysome fraction 0.70)
TOL_NN, TOL_LINK = 1.0, 0.05            # the frozen J4 (+- 1 nm) and J3b (+- 0.05) tolerances: units of the score


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def in_subset(P):
    b = np.floor(np.asarray(P, np.float64) / BLOCK_NM).astype(np.int64)
    return (b[:, 0] + 2 * b[:, 1] + 3 * b[:, 2]) % MOD == KEEP


def _cfg():
    from wholecell import config as C, registry
    cfg = C.load("jurkat")
    registry.load_plugins(cfg)
    return cfg


def run(a):
    from wholecell import registry
    from wholecell.preplaced import serial, builders, context  # noqa: F401
    cfg = _cfg()
    sc = copy.deepcopy(cfg["preplaced"]["serial"])
    bl = [b for b in sc["builders"] if b["name"] in ("npc", "microtubules", "polysomes")]
    for b in bl:
        if b["name"] == "polysomes":
            if a.linker == "transcripts":
                b["params"].pop("linker_nt_mean", None)
            else:
                b["params"]["linker_nt_mean"] = float(a.linker)
            if a.walk is not None:      # walk fix 2026-11-18: params.walk (candidate search) override, JSON
                b["params"]["walk"] = json.loads(a.walk) if a.walk else None
    t0 = time.time()
    ctx = registry.get("preplaced_context", sc["context"])(cfg, a.bdir, dict(sc.get("context_params") or {}))
    serial.run_builders(ctx, bl, out_dir=None)
    pos = ctx.pre["poly_pos"]; chn = ctx.pre["poly_chain"]
    st = dict(ctx.state.get("polysome_stats") or {}, linker=a.linker, walk=bl[-1]["params"].get("walk"), seconds=round(time.time() - t0))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez(a.out + ".npz", pos=pos, chain=chn, coll=ctx.pre["poly_coll"])
    json.dump(st, open(a.out + ".json", "w"), indent=1, default=float)
    log("STATS", json.dumps(st, default=float))


def complement(a):
    from wholecell.pack import geomload
    from wholecell.preplaced.context import CoreCell
    cfg = _cfg()
    cell = CoreCell(geomload.load(cfg, a.bdir))
    K = __import__("wholecell.preplaced.kern", fromlist=["core"]).core()
    cyto = K.lab("cytosol")
    rng = np.random.default_rng(a.seed)
    R = float(cell.G[0]) + 10.0
    out = []; n = 0
    while n < a.n:
        X = rng.uniform(-R, R, (2_000_000, 3))
        X = X[cell.label(X) == cyto]
        out.append(X); n += len(X)
    P = np.concatenate(out)[:a.n]
    np.save(a.out, P.astype(np.float32))
    log(f"complement: {len(P):,} uniform cytosol points (seed {a.seed}) -> {a.out}")


def em_points():
    """EM ribosome centres (raw, as J4 EM row) and de-duplicated (< 20 nm, as J3 / the point priors) in the model frame"""
    import pandas as pd
    from scipy.spatial import cKDTree
    from wholecell.pack import geomload
    from wholecell.preplaced.stage import _resolve
    cfg = _cfg()
    pp = [b for b in cfg["preplaced"]["serial"]["builders"] if b["name"] == "polysomes"][0]["params"]["point_priors"]
    E = pd.read_parquet(_resolve(cfg, None, pp["table"]))
    for c, v in (pp.get("filter_eq") or {}).items():
        E = E[E[c] == v]
    for c, v in (pp.get("filter_not_in") or {}).items():
        E = E[~E[c].isin(v)]
    centre = np.asarray(geomload.load(cfg, None).centre_nm, np.float64)
    P = E[["x", "y", "z"]].values * np.asarray(pp.get("voxel_nm")) - centre
    keep = np.ones(len(P), bool)
    for i, j in sorted(cKDTree(P).query_pairs(float(pp.get("dedup_nm", 20.0)))):
        if keep[i] and keep[j]:
            keep[j] = False
    return P, P[keep]


def _nn(P):
    from scipy.spatial import cKDTree
    d, _ = cKDTree(P).query(P, k=2)
    return d[:, 1]


def _linked(P, r=34.0):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    pr = cKDTree(P).query_pairs(r, output_type="ndarray")
    A = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(P), len(P)))
    _, lab = connected_components(A, directed=False)
    return np.bincount(lab)[lab] >= 2


def _within_group_nn(P, r=34.0):
    """NN distance of the points that have a neighbour within r (centre spacing inside polysome groups)"""
    d = _nn(P)
    return d, d <= r


def stats_sets(P_all, Pd, sub_all, sub_d):
    """median NN over the subset points (raw set), linked share over the subset points (dedup set), and the median
    within-group NN (dedup set, pairs <= 34 nm) over the subset points"""
    nn = _nn(P_all); lk = _linked(Pd); dg, ing = _within_group_nn(Pd)
    out = {}
    for nm, s1, s2 in (("subset", sub_all, sub_d), ("complement", ~sub_all, ~sub_d), ("whole", np.ones_like(sub_all), np.ones_like(sub_d))):
        out[nm] = dict(n=int(s1.sum()), nn_median=round(float(np.median(nn[s1])), 3), linked=round(float(lk[s2].mean()), 4),
                       group_nn_median=round(float(np.median(dg[s2 & ing])), 3))
    return out


def model_stats(PM, seed=0):
    """thinning exactly as b4_preplaced_check.ribosome_stats (streams seed, seed + 1)"""
    rng = np.random.default_rng(seed)
    thin = PM[rng.random(len(PM)) < EM_N / len(PM)]
    rng2 = np.random.default_rng(seed + 1)
    thin2 = PM[rng2.random(len(PM)) < EM_DEDUP_N / len(PM)]
    return stats_sets(thin, thin2, in_subset(thin), in_subset(thin2))


def score(a):
    E, Ed = em_points()
    em = stats_sets(E, Ed, in_subset(E), in_subset(Ed))
    log("EM", json.dumps(em))
    C = np.load(a.complement).astype(np.float64) if a.complement else None
    rows = []
    for f in a.runs.split(","):
        z = np.load(f); st = json.load(open(f[:-4] + ".json"))
        P = z["pos"].astype(np.float64)
        ch = z["chain"]
        same = ch[1:] == ch[:-1]
        cons = np.linalg.norm(P[1:] - P[:-1], axis=1)[same]
        r = dict(run=f, linker=st["linker"], stats=st, consecutive_centre_median=round(float(np.median(cons)), 3),
                 consecutive_centre_q10_q90=[round(float(x), 2) for x in np.quantile(cons, [0.1, 0.9])],
                 preplaced_only=model_stats(P))
        if C is not None:
            n_c = int(round(st["pool"] - len(P)))
            r["with_complement"] = model_stats(np.concatenate([P, C[:n_c]]))
            r["n_complement"] = n_c
        key = "with_complement" if C is not None else "preplaced_only"
        s = r[key]["subset"]
        r["score_subset"] = round(abs(s["nn_median"] - em["subset"]["nn_median"]) / TOL_NN
                                  + abs(s["linked"] - em["subset"]["linked"]) / TOL_LINK, 4)
        rows.append(r)
        log(r["linker"], "cons", r["consecutive_centre_median"], key, json.dumps(r[key]["subset"]), "score", r["score_subset"])
    for r in rows:      # the linker must keep the EM-set polysome fraction (J3a window, frozen): chains must not end short
        r["eligible_j3a"] = bool(J3A[0] <= r["stats"]["polysome_fraction_of_pool"] <= J3A[1])
    elig = [r for r in rows if r["eligible_j3a"]] or rows
    best = min(elig, key=lambda r: r["score_subset"])
    out = dict(method="argmin over the scanned values of |NN - NN_EM| / 1 nm + |linked - linked_EM| / 0.05 on the calibration "
                      "subset (model = pre-placed + random cytosolic complement for the packed free 80S, thinned as J3b / J4), "
                      "among the values that keep the polysome fraction in the J3a window 0.68-0.72",
               unconstrained_argmin=min(rows, key=lambda r: r["score_subset"])["linker"],
               subset=dict(block_nm=BLOCK_NM, rule="(i + 2 j + 3 k) mod 5 == 0, i = floor(x / 1000 nm), model frame"),
               em=em, em_frozen=dict(nn=EM_NN, linked=EM_LINKED), rows=rows, chosen=best["linker"])
    json.dump(out, open(a.out, "w"), indent=1, default=float)
    log("chosen", best["linker"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["run", "complement", "score"])
    ap.add_argument("--bdir"); ap.add_argument("--linker", default="transcripts"); ap.add_argument("--out")
    ap.add_argument("--n", type=int, default=700000); ap.add_argument("--seed", type=int, default=20261117)
    ap.add_argument("--runs"); ap.add_argument("--complement")
    ap.add_argument("--walk", default=None, help="JSON for params.walk of the polysome builder ('' = none; default: config)")
    a = ap.parse_args()
    {"run": run, "complement": complement, "score": score}[a.what](a)


if __name__ == "__main__":
    main()
