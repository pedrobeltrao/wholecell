#!/usr/bin/env python
"""B1 reduced benchmark (docs/unification_criteria/B1_packer_addendum_2026-11-18.md; frozen B1_packer.md A4/A5 made
concrete for the frozen Step-2 subsets): runs R-J, C1-C3 (Jurkat) and C1, R-Y, C2, C3 (yeast) on the Step-2 subset
tiles with the Step-2 neighbour rule, through the Step-2 harness (scripts/step2_subset_test.py of --harness: its
`_pack` (subset pack with the reference context) and `metrics` (subset metrics, harness code) subcommands and its DAG
helpers), then the addendum selection rule -> <out>/selection.json (+ results.json, results.tsv, status.txt).

    <core python> scripts/b1_reduced.py --commit SHA --harness DIR --out DIR [--workers 24]

Run folders: {results_root}/<cell>/_dev_unify/step2/b1r_<commit7>_<cell>_<id>/ (new; a folder of this driver is
resumed). Upstream (whole-cell stages up to counts + pre-placed index) is made ONCE per cell in the base run (Jurkat
R-J, yeast R-Y) as the harness does it (snapshot DAG jobs, library reused from the frame build when the library
section is identical, tiles when the geometry is identical), or linked from a donor Step-2 run whose code and
non-pack config are identical (recorded); other runs link the base run's upstream files and recompute the counts
when the shell cap differs (the only algorithm key the Jurkat / yeast counts depend on here). Phase A: Jurkat
upstream (12 cores) + Jurkat packs (12 workers) in parallel with yeast upstream (8) + yeast packs (12 workers).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time

import numpy as np

from _devsite import VC, RESULTS, YEAST as _YEAST_ROOT, PY as _SITE_PY, SMK as _SITE_SMK, GIT as _DEV_GIT  # noqa: E402,F401  (release: site paths)
R = f"{RESULTS}"
PY = _SITE_PY
GIT = _DEV_GIT
FRAME = {"jurkat": f"{R}/jurkat/20261117-bc121a6-s20261115-b3", "yeast": f"{R}/yeast/20261118-9977542-s0-a9-core"}
GATED = {"jurkat": ["cytosol", "nucleoplasm", "matrix", "ER_lumen", "golgi_lumen", "endo_lumen", "nucleolus", "speckle"],
         "yeast": ["cytosol", "nucleoplasm", "vacuole", "cell_wall"]}
SEAM_GATED = {"jurkat": "nucleoplasm", "yeast": "cytosol"}
# (cell, id, config, role); order of the runs (addendum §2)
RUNS = {"jurkat": [("s1", "config/cells/dev/b1r_jurkat_s1.yaml", "ref"), ("c1", "config/cells/dev/b1r_jurkat_c1.yaml", "cand"),
                   ("c2", "config/cells/dev/b1r_jurkat_c2.yaml", "cand"), ("c3", "config/cells/dev/b1r_jurkat_c3.yaml", "cand")],
        "yeast": [("c1", "config/cells/dev/b1r_yeast_c1.yaml", "cand"), ("r", "config/cells/yeast.yaml", "ref"),
                  ("c2", "config/cells/dev/b1r_yeast_c2.yaml", "cand"), ("c3", "config/cells/dev/b1r_yeast_c3.yaml", "cand")]}
BASE = {"jurkat": "s1", "yeast": "r"}
S1_OF = {"jurkat": "s1", "yeast": "c1"}            # condition 3 reference (frozen S1)
REF_OF = {"jurkat": "s1", "yeast": "r"}            # conditions 4, 5 reference arm
CANDS = ["c1", "c2", "c3"]
COUNT_FILES = ("counts.npz", "libs.json")
BUDGET_S = 5 * 3600; LAST_START_S = 4.5 * 3600
H = None; OUT = None; T_START = time.time(); LOCK = threading.Lock()


def log(*a):
    msg = time.strftime("%F %T ") + " ".join(str(x) for x in a)
    with LOCK:
        print(msg, flush=True)
        with open(os.path.join(OUT, "status.txt"), "a") as f:
            f.write(msg + "\n")


def jload(p, d=None):
    try:
        return json.load(open(p))
    except Exception:  # noqa: BLE001
        return d


def jsave(p, o):
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(o, f, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    os.replace(tmp, p)


def load_harness(hdir):
    spec = importlib.util.spec_from_file_location("step2_subset_test", os.path.join(hdir, "scripts", "step2_subset_test.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def rdir(commit, cell, rid):
    return os.path.join(R, cell, "_dev_unify", "step2", f"b1r_{commit[:7]}_{cell}_{rid}")


def snapshot(run, commit):
    code = os.path.join(run, "_code")
    if not os.path.exists(os.path.join(code, "SNAPSHOT_SHA")):
        os.makedirs(code, exist_ok=True)
        p1 = subprocess.Popen(["git", "-C", GIT, "archive", commit], stdout=subprocess.PIPE)
        subprocess.run(["tar", "-x", "-C", code], stdin=p1.stdout, check=True); p1.wait()
        full = subprocess.run(["git", "-C", GIT, "rev-parse", commit], capture_output=True, text=True).stdout.strip()
        open(os.path.join(code, "SNAPSHOT_SHA"), "w").write(full + "\n")
    cfgp = os.path.join(run, "cell.yaml")
    return code, cfgp


def new_run(commit, cell, rid, cfg_rel):
    run = rdir(commit, cell, rid)
    info = jload(os.path.join(run, "run.json"))
    if os.path.exists(run) and not (info and info.get("driver") == "b1_reduced"):
        raise SystemExit(f"{run} exists and is not a b1_reduced run folder")
    os.makedirs(os.path.join(run, "logs"), exist_ok=True)
    code, cfgp = snapshot(run, commit)
    if not os.path.exists(cfgp):
        with open(cfgp, "w") as f:
            f.write(f"# B1 reduced benchmark run {cell} {rid}: {cfg_rel} at {commit} (snapshot {code})\n"
                    f"extends: {os.path.join(code, cfg_rel)}\nbuild:\n  allow_dirty: true\n")
    if not info:
        info = dict(driver="b1_reduced", name=os.path.basename(run), cell=cell, run_id=rid, cell_id=cell, config=cfg_rel,
                    commit=commit, host=os.uname().nodename, harness=H.__file__, subset=H.SUBSETS[cell],
                    subset_md5=H.SUBSET_MD5[cell], context=H.CONTEXT[cell], started=time.strftime("%FT%T"), stages={})
        jsave(os.path.join(run, "run.json"), info)
    return run, cfgp, info


def eff_pack(code, cfgp):
    r = subprocess.run([PY, "-c", "import json,sys; from wholecell.config import load; from wholecell.pack import switches; "
                        "e=switches.effective(load(sys.argv[1])['pack']); print(json.dumps(dict(shell_cap=bool((e.get('shell_cap') or {}).get('enabled') "
                        "if isinstance(e.get('shell_cap'), dict) else e.get('shell_cap')), preset=e.get('_preset'), "
                        "org=bool((e.get('organisation') or {}).get('enabled')))))", cfgp],
                       cwd=code, capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=code, PYTHONHASHSEED="0"))
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


# ------------------------------------------------------------------ upstream (as the harness's cmd_run, steps 1-2)
def find_donor(cell, run, cfgp, commit):
    """a complete Step-2 run of another agent whose code (outside docs/scripts/tests/config/subsets/b1r overlays) and
    non-pack config sections are identical to this run's (e.g. the harness's B4 run of main): its upstream is linked"""
    base = os.path.join(R, cell, "_dev_unify", "step2")
    for d in sorted(os.listdir(base)) if os.path.isdir(base) else []:
        dd = os.path.join(base, d)
        if dd == run or d.startswith("b1r_"):
            continue
        ri = jload(os.path.join(dd, "run.json")) or {}
        st = ri.get("stages") or {}
        if not st or any(v.get("status") not in ("ok", "reused") for k, v in st.items() if k not in ("subset_pack", "metrics", "compare")):
            continue
        if st.get("upstream", {}).get("status") == "reused":
            continue                                      # only a run that made its own upstream
        if not all(os.path.exists(os.path.join(dd, "pack", f)) for f in ("counts.npz", "preplaced_off.npy", "tiles.parquet")):
            continue
        dc = (open(os.path.join(dd, "_code", "SNAPSHOT_SHA")).read().strip() if os.path.exists(os.path.join(dd, "_code", "SNAPSHOT_SHA")) else "")
        if not dc:
            continue
        diff = subprocess.run(["git", "-C", GIT, "diff", "--name-only", dc, commit], capture_output=True, text=True)
        if diff.returncode != 0:
            continue
        files = [f for f in diff.stdout.split() if f]
        if any(not (f.startswith(("docs/", "scripts/", "tests/", "config/subsets/")) or re.match(r"config/cells/dev/b1r_", f)) for f in files):
            continue
        ok = True
        for sec in ("geometry", "inventory", "library", "preplaced", "seeds"):
            same, _, _ = H.check_same_section(run, cfgp, dd, sec)
            ok = ok and same
        if ok:
            return dd, dc
    return None, None


def upstream(cell, run, cfgp, info, commit, workers):
    if info["stages"].get("upstream_done"):
        return
    env = dict(os.environ, PYTHONHASHSEED="0", NUMBA_CACHE_DIR=os.path.join(run, "numba_cache"))
    donor, dc = find_donor(cell, run, cfgp, commit)
    if donor:
        H.link_upstream(run, donor)
        for f in os.listdir(os.path.join(donor, "pack")):
            if f in ("subset_tiles.npy", "thinning.json", "counts_thinned.npz"):
                p = os.path.join(run, "pack", f)
                if os.path.islink(p):
                    os.unlink(p)
        same_pack, _, _ = H.check_same_section(run, cfgp, donor, "pack")
        info["stages"]["upstream"] = dict(status="reused", donor=donor, donor_commit=dc, pack_section_identical=same_pack,
                                          note="code identical outside docs/scripts/tests/subsets/overlays; geometry, inventory, "
                                               "library, preplaced, seeds sections identical")
        log(cell, "upstream linked from donor", donor, dc, "pack section identical:", same_pack)
        if not same_pack:                                 # counts of this run's own pack config
            for f in COUNT_FILES:
                p = os.path.join(run, "pack", f)
                if os.path.islink(p):
                    os.unlink(p)
            dt = counts(run, cfgp, os.path.join(run, "_code"))
            info["stages"]["pack_counts"] = dict(status="ok", seconds=round(dt, 1))
    else:
        targets = [os.path.join(run, "pack", "counts.npz"), os.path.join(run, "pack", "preplaced_off.npy")]
        jobs = H.dry_run_jobs(run, cfgp, workers, targets)
        info["dag_jobs"] = [j["rule"] for j in jobs]
        fr = FRAME[cell]
        same_lib, _, _ = H.check_same_section(run, cfgp, fr, "library")
        for j in jobs:
            name = j["rule"]
            if info["stages"].get(name, {}).get("status") in ("ok", "reused"):
                continue
            if same_lib and name == "library":         # only the core library; plugin steps library_* (yeast species table,
                                                      # engine/lib copy) write into this run's own folders and are run
                if not os.path.lexists(os.path.join(run, "library")):
                    os.symlink(os.path.join(fr, "library"), os.path.join(run, "library"))
                info["stages"][name] = dict(status="reused", donor=fr); continue
            if name == "pack_tiles":
                pj = lambda b: (jload(os.path.join(b, "geometry", "package.json")) or {}).get("manifest_sha256")
                if pj(run) is not None and pj(run) == pj(fr):
                    for f in ("tiles.parquet", "tiles_meta.json", "tile_labvol.parquet", "tile_geom"):
                        s = os.path.join(fr, "pack", f)
                        if os.path.exists(s) and not os.path.lexists(os.path.join(run, "pack", f)):
                            os.makedirs(os.path.join(run, "pack"), exist_ok=True)
                            os.symlink(s, os.path.join(run, "pack", f))
                    info["stages"][name] = dict(status="reused", donor=fr); continue
            t0 = time.time()
            rc = subprocess.run(["bash", "-c", j["cmd"]], env=env).returncode
            info["stages"][name] = dict(status="ok" if rc == 0 else "failed", rc=rc, seconds=round(time.time() - t0, 1))
            jsave(os.path.join(run, "run.json"), info)
            log(cell, f"  upstream {name}: rc {rc}, {time.time() - t0:.0f} s")
            if rc != 0:
                raise RuntimeError(f"{cell} upstream stage {name} failed (see {run}/logs)")
    info["upstream_pack_files"] = sorted(os.listdir(os.path.join(run, "pack")))
    info["stages"]["upstream_done"] = dict(status="ok", at=time.strftime("%FT%T"))
    jsave(os.path.join(run, "run.json"), info)


def counts(run, cfgp, code):
    t0 = time.time()
    with open(os.path.join(run, "logs", "pack_counts.log"), "w") as lf:
        rc = subprocess.run([PY, "-m", "wholecell", "pack", "--cell", cfgp, "--build", run, "--step", "counts"], cwd=code,
                            env=dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=code, NUMBA_NUM_THREADS="8", OMP_NUM_THREADS="1",
                                     NUMBA_CACHE_DIR=os.path.join(run, "numba_cache")), stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise RuntimeError(f"counts failed in {run}")
    return time.time() - t0


def prepare_cell(cell, commit, workers):
    """upstream of the base run; every other run of the cell linked to it, counts per shell-cap key"""
    rid_b = BASE[cell]
    cfg_b = dict((r, c) for r, c, _ in RUNS[cell])[rid_b]
    base, cfgp_b, info_b = new_run(commit, cell, rid_b, cfg_b)
    upstream(cell, base, cfgp_b, info_b, commit, workers)
    key_b = eff_pack(os.path.join(base, "_code"), cfgp_b)["shell_cap"]
    done_keys = {key_b: base}
    for rid, cfg_rel, _ in RUNS[cell]:
        if rid == rid_b:
            continue
        run, cfgp, info = new_run(commit, cell, rid, cfg_rel)
        if info["stages"].get("prepared"):
            done_keys.setdefault(eff_pack(os.path.join(run, "_code"), cfgp)["shell_cap"], run)
            continue
        for e in H.UPSTREAM_TOP:
            s = os.path.join(base, e)
            if os.path.lexists(s) and not os.path.lexists(os.path.join(run, e)):
                os.symlink(os.path.realpath(s), os.path.join(run, e))
        os.makedirs(os.path.join(run, "pack"), exist_ok=True)
        for e in info_b["upstream_pack_files"]:
            if e in COUNT_FILES or e in H.PACK_OUT or ".tmp" in e:
                continue
            d = os.path.join(run, "pack", e)
            if not os.path.lexists(d):
                os.symlink(os.path.realpath(os.path.join(base, "pack", e)), d)
        key = eff_pack(os.path.join(run, "_code"), cfgp)["shell_cap"]
        if key in done_keys:
            for f in COUNT_FILES:
                shutil.copy2(os.path.join(done_keys[key], "pack", f), os.path.join(run, "pack", f))
            info["stages"]["pack_counts"] = dict(status="copied", donor=done_keys[key], shell_cap=key)
        else:
            dt = counts(run, cfgp, os.path.join(run, "_code"))
            info["stages"]["pack_counts"] = dict(status="ok", seconds=round(dt, 1), shell_cap=key)
            done_keys[key] = run
            log(cell, rid, f"counts (shell cap {key}) {dt:.0f} s")
        info["stages"]["upstream"] = dict(status="linked", donor=base)
        info["stages"]["prepared"] = dict(status="ok")
        jsave(os.path.join(run, "run.json"), info)
    log(cell, "prepared: base", base, "counts keys", {str(k): os.path.basename(v) for k, v in done_keys.items()})


# ------------------------------------------------------------------ subset pack + metrics (harness)
def parse_pack_log(p):
    t = open(p).read() if os.path.exists(p) else ""
    pilot = sum(float(x) for x in re.findall(r"thinning pilot: .*?, (\d+)s", t))
    phases = sum(float(x) for x in re.findall(r"^phase \d+ done: .*?, (\d+)s$", t, re.M))
    tot = re.findall(r"subset pack done in (\d+) s", t)
    fac = re.findall(r"thinning: .* factors (\{.*\})", t)
    return dict(pilot_s=pilot, phases_s=phases, harness_pack_s=float(tot[-1]) if tot else None, thinning_factors_log=fac[-1] if fac else None)


def pack_and_metrics(cell, rid, commit, workers, timeout):
    run = rdir(commit, cell, rid)
    info = jload(os.path.join(run, "run.json"))
    st = info["stages"]
    if st.get("metrics", {}).get("status") == "ok" or st.get("subset_pack", {}).get("status") in ("timeout", "failed", "not_run"):
        return info
    if time.time() - T_START > LAST_START_S:
        st["subset_pack"] = dict(status="not_run", reason="budget (4 h 30 min elapsed)"); jsave(os.path.join(run, "run.json"), info)
        log(cell, rid, "NOT RUN (budget)"); return info
    code = os.path.join(run, "_code")
    if st.get("subset_pack", {}).get("status") != "ok":
        for e in H.PACK_OUT | {"thinning.json", "counts_thinned.npz", "subset_tiles.npy"}:   # an interrupted pack: start again
            p = os.path.join(run, "pack", e)
            if os.path.islink(p):
                raise RuntimeError(f"{p} is a link (would write into another folder)")
            if os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.exists(p):
                os.remove(p)
        lp = os.path.join(run, "logs", "subset_pack.log")
        t0 = time.time()
        with open(lp, "w") as lf:
            pr = subprocess.Popen([PY, H.__file__, "_pack", "--run", run, "--cell", cell, "--workers", str(workers)], cwd=code,
                                  env=dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=code, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                                           MKL_NUM_THREADS="1", NUMBA_CACHE_DIR=os.path.join(run, "numba_cache"),
                                           # the values the pack step sets before numba loads (the harness imports numba
                                           # earlier: without them a JIT compile in a forked worker clashes)
                                           NUMBA_NUM_THREADS="1", NUMBA_THREADING_LAYER="workqueue"),
                                  stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)
            try:
                rc = pr.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                os.killpg(pr.pid, signal.SIGKILL); pr.wait(); rc = "timeout"
        dt = time.time() - t0
        pl = parse_pack_log(lp)
        wall = pl["phases_s"] if rc == 0 else None        # rounds + squeeze (sum of the phase times; no input loading, no pilot)
        st["subset_pack"] = dict(status="ok" if rc == 0 else ("timeout" if rc == "timeout" else "failed"), rc=rc, seconds=round(dt, 1),
                                 timeout_s=timeout, workers=workers, wall_cond_s=wall, **pl)
        jsave(os.path.join(run, "run.json"), info)
        log(cell, rid, f"subset pack {st['subset_pack']['status']}: {dt:.0f} s total, pilot {pl['pilot_s']:.0f} s, "
                       f"wall (rounds + squeeze) {wall}")
        if rc != 0:
            return info
    t0 = time.time()
    with open(os.path.join(run, "logs", "metrics.log"), "w") as lf:
        rc = subprocess.run([PY, H.__file__, "metrics", "--run", run], cwd=os.path.dirname(os.path.dirname(H.__file__)),
                            env=dict(os.environ, PYTHONHASHSEED="0", OMP_NUM_THREADS="1", NUMBA_NUM_THREADS="3",
                                     NUMBA_CACHE_DIR=os.path.join(run, "numba_cache")), stdout=lf, stderr=subprocess.STDOUT).returncode
    info = jload(os.path.join(run, "run.json")); st = info["stages"]
    st["metrics"] = dict(status="ok" if rc == 0 else "failed", rc=rc, seconds=round(time.time() - t0, 1))
    jsave(os.path.join(run, "run.json"), info)
    if cell == "jurkat" and rid != "s1" and rc == 0:
        packtree_c7(run)
    log(cell, rid, f"metrics rc {rc}, {time.time() - t0:.0f} s")
    return info


def packtree_c7(run):
    """Jurkat runs with packing shapes: C7(a) also on the packing trees (report; the gate is on the library trees)"""
    code = f"""
import json, glob, os, sys, numpy as np
sys.path.insert(0, {os.path.dirname(os.path.dirname(H.__file__))!r})
from wholecell.config import load
from wholecell.validate import common, placement
run = {run!r}
cfg = load(os.path.join(run, 'cell.yaml'))
ids = np.load(os.path.join(run, 'subset_tiles.npy'))
pt = sorted(glob.glob(os.path.join(run, 'library', 'pack_shapes', 'pack_trees_*.npz')))[-1]
common.CACHE_DIR = os.path.join(run, 'metrics', 'cache', 'packtrees'); os.makedirs(common.CACHE_DIR, exist_ok=True)
inp = common.Inputs(cfg, run, tiles_dir=os.path.join(run, 'metrics', 'views', 'run'), pack_trees=pt)
P = placement.run(cfg, run, inp, os.path.join(run, 'metrics', 'placement_packtrees.json'), os.path.join(run, 'metrics', 'placement_packtrees_tiles.tsv'),
                  overwrite=True, nproc=8, tiles=ids)
n = P['n_packed']
json.dump(dict(pack_trees=pt, packed_packed_per_1e6=P['clash_pairs']['packed_packed'] / n * 1e6 if n else None,
               clash_rate_per_1e6=P['clash_rate_per_1e6'], containment_rate_per_1e6=P['containment_rate_per_1e6']),
          open(os.path.join(run, 'c7_packtrees.json'), 'w'), indent=1)
"""
    r = subprocess.run([PY, "-c", code], capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED="0", NUMBA_NUM_THREADS="3"))
    if r.returncode != 0:
        log(os.path.basename(run), "pack-tree C7 (report) failed:", r.stderr[-500:])


# ------------------------------------------------------------------ selection (addendum §3)
def summarise(cell, rid, commit):
    run = rdir(commit, cell, rid)
    info = jload(os.path.join(run, "run.json")) or {}
    res = jload(os.path.join(run, "result.json"))
    sp = (info.get("stages") or {}).get("subset_pack", {})
    s = dict(run=run, status=sp.get("status", "missing") if (info.get("stages") or {}).get("metrics", {}).get("status") == "ok" or sp.get("status") != "ok"
             else "metrics_failed", wall_s=sp.get("wall_cond_s"), pilot_s=sp.get("pilot_s"), pack_total_s=sp.get("seconds"),
             tile_cpu_s=((res or {}).get("timing") or {}).get("tile_cpu_s"))
    if s["status"] != "ok" or not res:
        if s["status"] == "ok":
            s["status"] = "no_result"
        return s
    c10 = res["c10"]; reg = c10["regions"]
    s.update(placed_mass_fraction=c10["placed_mass_fraction"], placed_of_unthinned=c10["placed_of_unthinned"],
             regions={r: dict(placed_mass_frac=reg[r]["placed_mass_frac"], placed_of_unthinned=reg[r]["placed_of_unthinned"],
                              requested=reg[r]["requested"]) for r in GATED[cell] if r in reg},
             seam=res.get("seam"), c7a=res.get("c7a"), c8=res.get("c8"),
             thinning=(res.get("stage") or {}).get("thinning_factors"), instances=c10.get("instances"))
    pt = jload(os.path.join(run, "c7_packtrees.json"))
    if pt:
        s["c7a_packtrees_report"] = pt
    return s


def conditions(cell, s, s1, ref):
    c = {}
    bad, rep = {}, {}
    for r, v in s["regions"].items():
        n = v["requested"]; f = v["placed_mass_frac"]
        if f is None:
            continue
        if n < 100:
            rep[r] = dict(value=f, n=n, note="< 100 requested copies on the subset: reported, not gated"); continue
        lim = 0.98 - 3 * np.sqrt(0.98 * 0.02 / n)
        if f < lim:
            bad[r] = dict(value=round(f, 5), limit=round(lim, 5), n=n)
    worst = min(((v["placed_mass_frac"], r) for r, v in s["regions"].items() if v["placed_mass_frac"] is not None), default=(None, None))
    c["1_regions"] = dict(passed=not bad, failing=bad, reported_only=rep, worst=dict(region=worst[1], value=worst[0]))
    sm = (s.get("seam") or {}).get(SEAM_GATED[cell]) or {}
    if sm.get("index") is None or (sm.get("n_interior") or 0) < 1000:
        c["2_seam"] = dict(passed=True, gated=False, value=sm.get("index"), n_interior=sm.get("n_interior"), note="not measurable: deferred to K3")
    else:
        c["2_seam"] = dict(passed=abs(sm["index"] - 1) <= 0.1, gated=True, value=sm["index"], n_tiles=sm.get("n_tiles"), n_interior=sm.get("n_interior"))
    if s1 and s1.get("status") == "ok":
        lim = s1["placed_of_unthinned"] - 0.003
        c["3_unthinned"] = dict(passed=s["placed_of_unthinned"] >= lim, value=s["placed_of_unthinned"], limit=lim)
    else:
        c["3_unthinned"] = dict(passed=False, value=s["placed_of_unthinned"], limit=None, note="S1 reference missing")
    pp = s["c7a"]["packed_packed_per_1e6"]
    ok4 = pp is not None and pp <= 1.0; lims = {}
    if ref and ref.get("status") == "ok":
        for k, v, rv in (("clash_rate_per_1e6", s["c7a"]["clash_rate_per_1e6"], ref["c7a"]["clash_rate_per_1e6"]),
                         ("containment_rate_per_1e6", s["c8"]["containment_rate_per_1e6"], ref["c8"]["containment_rate_per_1e6"])):
            if rv is None:
                continue
            lim = rv + max(0.1 * abs(rv), 1.0); lims[k] = dict(value=v, limit=lim); ok4 = ok4 and v is not None and v <= lim
    else:
        ok4 = False; lims["note"] = "reference arm missing"
    c["4_clash_containment"] = dict(passed=ok4, packed_packed_per_1e6=pp, **lims)
    rw = (ref or {}).get("wall_s")
    c["5_wall"] = dict(passed=bool(s["wall_s"] is not None and rw and s["wall_s"] <= 2 * rw), value=s["wall_s"], limit=2 * rw if rw else None,
                       note="subset pack wall = rounds + squeeze (sum of phase times; input loading and the whole-cell pilot "
                            "reported separately); binding check K7 at the combined build")
    return c


def select(commit):
    sys.path.insert(0, os.path.join(rdir(commit, "jurkat", "s1"), "_code"))
    from wholecell.pack import switches
    S = {cell: {rid: summarise(cell, rid, commit) for rid, _, _ in RUNS[cell]} for cell in RUNS}
    table = {}
    for cell in S:
        s1 = S[cell][S1_OF[cell]]; ref = S[cell][REF_OF[cell]]
        for rid in CANDS:
            s = S[cell][rid]
            if s["status"] != "ok":
                table.setdefault(cell, {})[rid] = dict(status=s["status"], eligible=False, n_failed=5, cost=None); continue
            c = conditions(cell, s, s1, ref)
            rw = ref.get("wall_s") if ref.get("status") == "ok" else None
            table.setdefault(cell, {})[rid] = dict(status="ok", eligible=all(v["passed"] for v in c.values()),
                                                   n_failed=sum(not v["passed"] for v in c.values()), conditions=c,
                                                   cost=(s["wall_s"] / rw) if (s["wall_s"] and rw) else None)
    cand = []
    for rid in CANDS:
        tj, ty = table["jurkat"][rid], table["yeast"][rid]
        cost = tj["cost"] + ty["cost"] if (tj["cost"] is not None and ty["cost"] is not None) else None
        cand.append(dict(candidate=rid, eligible_both=bool(tj["eligible"] and ty["eligible"]), cost=cost, n_failed=tj["n_failed"] + ty["n_failed"],
                         eligible=dict(jurkat=tj["eligible"], yeast=ty["eligible"])))
    el = [c for c in cand if c["eligible_both"] and c["cost"] is not None]
    if el:
        chosen = min(el, key=lambda c: c["cost"])["candidate"]; verdict = "SELECTED"
    else:
        chosen = min(cand, key=lambda c: (c["n_failed"], c["cost"] if c["cost"] is not None else 1e9))["candidate"]; verdict = "NO ELIGIBLE SETTING"
    # P8 shapes (addendum §3): Jurkat C1 vs R-J
    p, q = S["jurkat"]["c1"], S["jurkat"]["s1"]
    why = []
    if p["status"] == "ok" and q["status"] == "ok":
        if not (p["c7a"]["packed_packed_per_1e6"] is not None and p["c7a"]["packed_packed_per_1e6"] <= 1.0):
            why.append(f"C7(a) packed-packed on library trees {p['c7a']['packed_packed_per_1e6']} > 1")
        if p["placed_of_unthinned"] < q["placed_of_unthinned"] - 0.003:
            why.append(f"placed_of_unthinned {p['placed_of_unthinned']:.4f} < R-J {q['placed_of_unthinned']:.4f} - 0.003")
        for r, v in q["regions"].items():
            a_, b_ = (p["regions"].get(r) or {}).get("placed_mass_frac"), v["placed_mass_frac"]
            n = v["requested"]
            if a_ is not None and b_ is not None and n >= 100 and a_ < b_ - 0.003 - 3 * np.sqrt(0.98 * 0.02 / n):
                why.append(f"{r} {a_:.4f} < R-J {b_:.4f} - 0.003 - binomial")
        if (p["wall_s"] or 1e18) > (q["wall_s"] or 0):
            why.append(f"wall {p['wall_s']} > R-J {q['wall_s']}")
        shapes = dict(decision="library" if why else "pack_shapes", reason="; ".join(why) or "all C1-vs-R-J conditions met")
    else:
        shapes = dict(decision="undecided", reason=f"C1 {p['status']} / R-J {q['status']}")
    preset = {"c1": ["b1_s1"], "c2": ["b1_s4"], "c3": ["b1_s4y"]}[chosen]
    cd = {}
    for n_ in preset:
        cd = switches._merge(cd, switches.B1_BENCHMARK[n_])
    cd.pop("shapes", None)
    cd["budgets"] = {"per_region": {}}
    if shapes["decision"] == "pack_shapes":
        cd["shapes"] = {"source": "pack_shapes"}
    shapes_note = None
    if shapes["decision"] != "pack_shapes":
        shapes_note = ("the candidates were run with packing shapes (Jurkat) but the shapes check did not pass: the shapes value "
                       "of the common default is for Pedro to decide (dict given without shapes = library trees)")
    sel = dict(verdict=verdict, chosen=chosen, candidates=cand, table=table, runs=S, shapes=shapes, shapes_note=shapes_note,
               common_defaults=cd if verdict == "SELECTED" else None, proposal=cd,
               apply="nothing applied automatically: Pedro decides; then paste into wholecell/pack/switches.py COMMON_DEFAULTS",
               criteria=["docs/unification_criteria/B1_packer.md (md5 f544eb73724b2accdd9c2fd0bec88b89)",
                         "docs/unification_criteria/B1_packer_addendum_2026-11-18.md"],
               commit=commit, harness=H.__file__, elapsed_s=round(time.time() - T_START), written=time.strftime("%FT%T"))
    jsave(os.path.join(OUT, "results.json"), S)
    rows = ["cell\trun\tstatus\twall_s\tpilot_s\tplaced_mass_fraction\tplaced_of_unthinned\tworst_region\tworst\tseam\tc7a_pp_per_1e6\tc7a_clash_per_1e6\tc8_per_1e6\teligible\tn_failed\tcost"]
    for cell in S:
        for rid, s in S[cell].items():
            t = table[cell].get(rid, {})
            if s["status"] == "ok":
                w = min(((v["placed_mass_frac"], r) for r, v in s["regions"].items() if v["placed_mass_frac"] is not None), default=(None, None))
                sm = (s.get("seam") or {}).get(SEAM_GATED[cell]) or {}
                rows.append("\t".join(str(x) for x in (cell, rid, "ok", s["wall_s"], s["pilot_s"], round(s["placed_mass_fraction"], 5),
                                                         round(s["placed_of_unthinned"], 5), w[1], w[0] and round(w[0], 5), sm.get("index"),
                                                         s["c7a"]["packed_packed_per_1e6"], s["c7a"]["clash_rate_per_1e6"],
                                                         s["c8"]["containment_rate_per_1e6"], t.get("eligible", "ref"), t.get("n_failed", ""), t.get("cost", ""))))
            else:
                rows.append(f"{cell}\t{rid}\t{s['status']}" + "\t" * 13)
    open(os.path.join(OUT, "results.tsv"), "w").write("\n".join(rows) + "\n")
    jsave(os.path.join(OUT, "selection.json"), sel)
    log(f"selection: {verdict}: {chosen}; shapes {shapes['decision']} ({shapes['reason']})")
    return sel


# ------------------------------------------------------------------ main
def main():
    global H, OUT
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--commit", required=True); ap.add_argument("--harness", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=24); ap.add_argument("--select-only", action="store_true")
    ap.add_argument("--retry-failed", action="store_true", help="runs whose subset pack FAILED (not timed out) are run again")
    a = ap.parse_args()
    OUT = os.path.abspath(a.out); os.makedirs(OUT, exist_ok=True)
    H = load_harness(a.harness)
    commit = subprocess.run(["git", "-C", GIT, "rev-parse", a.commit], capture_output=True, text=True).stdout.strip()
    if a.select_only:
        select(commit); return
    log(f"START b1 reduced (commit {commit}, harness {a.harness}, host {os.uname().nodename})")
    if a.retry_failed:                                    # failed packs, and timeouts set from a failed reference arm
        for cell in RUNS:
            rf = jload(os.path.join(rdir(commit, cell, REF_OF[cell]), "run.json")) or {}
            ref_failed = (rf.get("stages", {}).get("subset_pack") or {}).get("status") == "failed"
            for rid, _, _ in RUNS[cell]:
                f = os.path.join(rdir(commit, cell, rid), "run.json"); info = jload(f)
                stt = ((info or {}).get("stages", {}).get("subset_pack") or {}).get("status")
                if stt == "failed" or (stt == "timeout" and ref_failed):
                    info.setdefault("retries", []).append(info["stages"].pop("subset_pack")); info["stages"].pop("metrics", None)
                    jsave(f, info); log(cell, rid, f"subset pack {stt}: retried")
    errs = {}

    def jurkat_all():                                     # Jurkat: upstream (12 cores), then its packs (12 workers)
        try:
            prepare_cell("jurkat", commit, 12)
            tref = None
            for rid, _, role in RUNS["jurkat"]:
                to = 3600.0 if role == "ref" else max(1800.0, 2.2 * (tref or 3600.0 / 2.2))
                info = pack_and_metrics("jurkat", rid, commit, 12, to)
                if role == "ref":
                    tref = (info["stages"].get("subset_pack") or {}).get("seconds")
        except Exception as e:  # noqa: BLE001
            errs["jurkat"] = repr(e); log("jurkat FAILED", repr(e))

    th = threading.Thread(target=jurkat_all); th.start()
    # yeast: upstream (8 cores), then its packs (12 workers), in parallel with Jurkat (<= 24 cores in all; every subset
    # colour has <= 12 tiles, so 12 workers pack a colour in one wave)
    try:
        prepare_cell("yeast", commit, 8)
        tref = None
        for rid, _, role in RUNS["yeast"]:
            if role == "ref" or tref is None:
                to = 5400.0 if role == "ref" else 3600.0
            else:
                to = max(1800.0, 2.2 * tref)
            info = pack_and_metrics("yeast", rid, commit, 12, to)
            if role == "ref":
                tref = (info["stages"].get("subset_pack") or {}).get("seconds") or 5400.0 / 2.2
    except Exception as e:  # noqa: BLE001
        errs["yeast"] = repr(e); log("yeast FAILED", repr(e))
    th.join()
    sel = select(commit)
    sel["errors"] = errs; jsave(os.path.join(OUT, "selection.json"), sel)
    log(f"DONE in {time.time() - T_START:.0f} s")


if __name__ == "__main__":
    main()
