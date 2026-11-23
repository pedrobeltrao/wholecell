# new in wholecell (2026-11-20: Jurkat freeze build died at inventory_copies after 2 h; check the inventory stage alone)
"""Run the whole inventory stage of a cell's CURRENT config end to end (inventory only, minutes) in a scratch build:

    python scripts/inventory_e2e.py --cell jurkat --out <scratch dir> [--upstream-from <build>]

frozen (md5 of every frozen input) + every step of wholecell.inventory.stage.steps(cfg) in plan order, each as
`python -m wholecell inventory --step <s>` (the build's own command line, NUMBA env as workflow WC()). geometry/ (measures) is symlinked and preplaced/ (Jurkat rule v7,
lamina consumption) copied from --upstream-from (a build of the same geometry package). Exit 0 only if every step succeeds and
the last step's outputs exist; <out>/inventory_e2e.json records steps, seconds, return codes.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--upstream-from", default=None, help="build whose geometry/ (symlinked, read only) and preplaced/ (copied) are reused")
    ap.add_argument("--skip", nargs="*", default=["compare"])
    a = ap.parse_args(argv)
    from wholecell import registry
    from wholecell.config import load
    cfg = load(a.cell)
    registry.load_plugins(cfg)
    from wholecell.inventory import copies, generic, secms  # noqa: F401
    from wholecell.inventory import stage as inv
    plan = inv.steps(cfg)
    if os.path.exists(a.out):
        raise SystemExit(f"{a.out} exists")
    os.makedirs(a.out)
    if a.upstream_from:
        os.symlink(os.path.join(os.path.abspath(a.upstream_from), "geometry"), os.path.join(a.out, "geometry"))
        shutil.copytree(os.path.join(a.upstream_from, "preplaced"), os.path.join(a.out, "preplaced"), symlinks=True)
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=REPO, NUMBA_NUM_THREADS="8", NUMBA_THREADING_LAYER="workqueue",
               OMP_NUM_THREADS=os.environ.get("OMP_NUM_THREADS", "2"))
    rec = []
    py = sys.executable
    for stage, step in [("frozen", None)] + [("inventory", s) for s in plan if s not in a.skip]:
        cmd = [py, "-m", "wholecell", stage, "--cell", a.cell, "--build", a.out] + (["--step", step] if step else [])
        t = time.time()
        os.makedirs(os.path.join(a.out, "logs"), exist_ok=True)
        lf = os.path.join(a.out, "logs", f"{stage}_{(step or '').replace(':', '_')}.log")
        with open(lf, "w") as fh:
            rc = subprocess.run(cmd, cwd=REPO, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
        rec.append(dict(stage=stage, step=step, rc=rc, seconds=round(time.time() - t, 1), log=lf))
        print(stage, step, "rc", rc, f"{time.time() - t:.0f} s", flush=True)
        if rc:
            print(open(lf).read()[-3000:], flush=True)
            break
    last = [s for s in plan if s not in a.skip][-1]
    missing = [x for x in plan[last]["outputs"] if not os.path.exists(os.path.join(a.out, x))]
    ok = all(r["rc"] == 0 for r in rec) and len(rec) == 1 + len([s for s in plan if s not in a.skip]) and not missing
    json.dump(dict(cell=a.cell, ok=ok, steps=rec, missing_final_outputs=missing), open(os.path.join(a.out, "inventory_e2e.json"), "w"), indent=1)
    print("inventory e2e", "OK" if ok else "FAILED", flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
