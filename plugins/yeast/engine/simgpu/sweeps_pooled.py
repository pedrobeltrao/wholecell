# origin: yeast yeast_cell_dev_bundle/simgpu/sweeps_pooled.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Diffusion sweeps of this build, pooled over several boxes (2026-11-14), for validation/diffusion.py.

Until build S the crowding term of the diffusion model came from sweeps of one box per compartment of an earlier build,
copied into every run. One 120 nm box is too small a sample: re-run on builds S and U2 it differed by more than the
change being tested (its volume fraction was 0.133 against 0.178 averaged over 8 boxes). Now every build runs the
cytoplasm and GEM sweeps in N_BOX boxes centred on tiles that are >99.5% cytosol (seeded draw), and pools them into
work/sweeps/{cyto120_mf,gem200_mf}.json; the nucleoplasm sweep (nuc100_mf) is run once at its usual centre.
  python simgpu/sweeps_pooled.py [--gpus 4]
"""
import os, sys, json, argparse, subprocess
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
import config as C  # noqa: E402
import webenc as we  # noqa: E402

N_BOX, JOBS = 8, ("gem200_mf", "cyto120_mf")


def centres():
    lv = pd.read_parquet(f"{C.W3D}/tile_labvol.parquet"); tiles = we.tile_index_list()
    v = lv.values; cyt = list(lv.columns).index("cytosol")
    ok = np.where(v[:, cyt] > 0.995 * v.sum(1))[0]
    rng = np.random.default_rng(C.seed(11)); picks = []
    for i in rng.permutation(ok):
        c = np.array(tiles[i], float) * we.TILE + 0.5 * we.TILE
        if all(np.linalg.norm(c - q) > 400.0 for q in picks):
            picks.append(c)
        if len(picks) == N_BOX:
            break
    return picks


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--gpus", type=int, default=4); a = ap.parse_args()
    out = f"{C.WORK}/sweeps"; tmp = f"{C.WORK}/sweeps_multi"; os.makedirs(out, exist_ok=True); os.makedirs(tmp, exist_ok=True)
    cs = centres()
    jobs = [(j, k) for k in range(len(cs)) for j in JOBS] + [("nuc100_mf", -1)]
    py, sw = sys.executable, os.path.join(HERE, "sweep.py")
    running = []
    for n, (j, k) in enumerate(jobs):
        env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(n % a.gpus))
        if k >= 0:
            env["YC_SWEEP_CENTER"] = ",".join(f"{x:.1f}" for x in cs[k])
        fn = f"{tmp}/{j}_{k}.json" if k >= 0 else f"{out}/{j}.json"
        running.append(subprocess.Popen([py, sw, j, fn], env=env, stdout=open(fn.replace(".json", ".log"), "w"), stderr=subprocess.STDOUT))
        if len(running) == a.gpus:
            running.pop(0).wait()
    for p in running:
        p.wait()
    for j in JOBS:
        P = [json.load(open(f"{tmp}/{j}_{k}.json")) for k in range(len(cs))]
        res = dict(job=j, n=int(sum(p["n"] for p in P)), phi=float(np.mean([p["phi"] for p in P])), phi_sd=float(np.std([p["phi"] for p in P])),
                   wall_s=float(sum(p["wall_s"] for p in P)), boxes=len(P), centres=[list(map(float, c)) for c in cs])
        for key in ("rh", "D", "D0", "tracer", "ids"):
            res[key] = [x for p in P for x in p[key]]
        json.dump(res, open(f"{out}/{j}.json", "w"))
        print(f"{j}: {len(P)} boxes, volume fraction {res['phi']:.3f} +- {res['phi_sd']:.3f}")


if __name__ == "__main__":
    main()
