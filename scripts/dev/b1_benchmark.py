#!/usr/bin/env python
"""B1 packer benchmark (criteria docs/unification_criteria/B1_packer.md, Part A): the benchmark settings of
wholecell/pack/switches.py B1_BENCHMARK on a 500-tile Jurkat subset and on the whole yeast cell, from copies of the
inputs of the latest complete build of each cell; metrics; the frozen selection rule -> selection.json.

    <core python> scripts/b1_benchmark.py all --out DIR [--workers 24] [--cells jurkat,yeast] [--settings ...]
    <core python> scripts/b1_benchmark.py {prepare | run | metrics | select} --out DIR [--cell C --setting S]

Run from a code snapshot (the queue script _queue/b1_benchmark.sh). Writes only under --out; input builds are read
(their files linked into fresh run folders, A1). Resumable: a run whose results.json exists is not repeated.
Outputs: inputs.json, <cell>/<setting>/{cell.yaml, pack/, run.json, results.json, logs}, results.tsv, results.json,
selection.json, status.txt.
"""
from __future__ import annotations

import argparse
import copy
import glob
import hashlib
import json
import os
import re
import subprocess
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from _devsite import VC, RESULTS, YEAST as _YEAST_ROOT, PY as _SITE_PY, SMK as _SITE_SMK, GIT as _DEV_GIT  # noqa: E402,F401  (release: site paths)
R = f"{RESULTS}"
SUBSET_MD5 = "0877fd0647b42b82db7d14bf5a64998e"          # criteria A2 (B6, B3, v7repro2 builds)
JURKAT_GATED = ["cytosol", "nucleoplasm", "matrix", "ER_lumen", "golgi_lumen", "endo_lumen", "nucleolus", "speckle"]
ORDER = {"jurkat": ["b1_s1", "b1_s2", "b1_s3", "b1_s4", "b1_s4y", "b1_s4p"],
         "yeast": ["b1_s1", "b1_s2", "b1_s3", "b1_s4", "b1_s4y", "b1_s5"]}
COMMON_CANDIDATES = ["b1_s1", "b1_s2", "b1_s3", "b1_s4", "b1_s4y"]
PACK_OUTPUTS = {"counts.npz", "libs.json", "counts_thinned.npz", "thinning.json", "tiles", "rounds", "pilot",
                "pack_log.tsv", "pack_done.json", "summary.json"}
SKIP_TOP = {"pack", "validate", "logs", "export", "snakemake.log", "launch.log", "_code", "numba_cache", "compare",
            ".snakemake", "launch.nohup"}
PY = sys.executable


def log(out, *a):
    msg = time.strftime("%F %T ") + " ".join(str(x) for x in a)
    print(msg, flush=True)
    with open(os.path.join(out, "status.txt"), "a") as f:
        f.write(msg + "\n")


def jload(p, default=None):
    try:
        return json.load(open(p))
    except Exception:
        return default


def jsave(p, obj):
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    os.replace(tmp, p)


# ------------------------------------------------------------------ inputs (criteria A1)
def _status_build(status_file):
    """the build dir named in a full_build.status ('RUNNING build X' or 'build X (snakemake')"""
    if not os.path.exists(status_file):
        return None
    hits = re.findall(r"(?:RUNNING build|\bbuild) (/\S+)", open(status_file).read())
    return hits[-1] if hits else None


def _complete(b):
    return bool(b) and os.path.isdir(b) and os.path.exists(os.path.join(b, "pack", "pack_done.json")) and \
        os.path.exists(os.path.join(b, "snakemake.log")) and "(100%) done" in open(os.path.join(b, "snakemake.log")).read()


def _pack_wall_snakemake(b):
    """wall time (s) of the pack_pack rule of a snakemake build (rule start -> 'Finished job <id>')"""
    import datetime
    f = os.path.join(b, "snakemake.log")
    if not os.path.exists(f):
        return None
    ts = None; start = None; jid = None; want = False
    for line in open(f):
        m = re.match(r"\[(\w{3} \w{3} +\d+ \d\d:\d\d:\d\d \d{4})\]", line)
        if m:
            ts = datetime.datetime.strptime(re.sub(" +", " ", m.group(1)), "%a %b %d %H:%M:%S %Y")
            continue
        if re.match(r"\s*(local)?rule pack_pack:", line):
            start = ts; want = True; continue
        if want and "jobid:" in line:
            jid = line.split("jobid:")[1].strip(); want = False; continue
        if start is not None and jid is not None and re.match(rf"Finished job {jid}\.", line.strip()):
            return (ts - start).total_seconds()
    return None


def resolve_inputs():
    """{cell: {build, kind, base_config, complete, pack_wall_s, host}} (criteria A1 order)"""
    out = {}
    for cell, chain in (("jurkat", [("b4", f"{R}/jurkat/_dev_unify/b4/full_build.status"),
                                    ("b3", f"{R}/jurkat/_dev_unify/b3/full_build.status"),
                                    ("b6", f"{R}/jurkat/20261117-04dec6f-s20261115-b6")]),
                        ("yeast", [("b4", f"{R}/yeast/_dev_unify/b4/full_build.status"),
                                   ("b6", f"{R}/yeast/_dev_unify/b6/full_build.status"),
                                   ("a9", f"{R}/yeast/_dev_unify/a9/full_build.status")])):
        for tag, ref in chain:
            b = ref if os.path.isdir(ref) else _status_build(ref)
            if _complete(b):
                host = ""
                lf = os.path.join(b, "launch.log")
                if os.path.exists(lf):
                    host = (re.findall(r"^host: (\S+)", open(lf).read(), re.M) or [""])[0]
                out[cell] = dict(build=b, kind=tag, base_config=os.path.join(b, "_code", "config", "cells", f"{cell}.yaml"),
                                 complete=True, pack_wall_s=_pack_wall_snakemake(b), host=host)
                bc = open(lf).read() if os.path.exists(lf) else ""
                m = re.search(r"^cell: (\S+)", bc, re.M)
                if m:                                     # config the build was made with (e.g. yeast_x5_core)
                    out[cell]["base_config"] = os.path.join(b, "_code", "config", "cells", f"{m.group(1)}.yaml")
                break
        else:
            if cell == "yeast":                            # A7 core run (not a snakemake build): its inputs + snapshot code
                a7 = f"{R}/yeast/_dev_unify/a7y/full_7a70243"
                st = open(os.path.join(a7, "logs", "status.txt")).read()
                m = re.search(r"pack_pack done (\d+) s", st)
                out[cell] = dict(build=a7, kind="a7_run", base_config=os.path.join(a7, "dev_yeast_a7y_full.yaml"),
                                 complete=os.path.exists(os.path.join(a7, "pack", "pack_done.json")),
                                 pack_wall_s=float(m.group(1)) if m else None, host="eu-g4-002")
            else:
                raise SystemExit(f"no complete input build for {cell}")
    return out


# ------------------------------------------------------------------ Jurkat subset (criteria A2)
def jurkat_subset(tiles_parquet, n=500, w=5):
    import pandas as pd
    t = pd.read_parquet(tiles_parquet)
    tot = t[JURKAT_GATED].sum()
    best = None
    for a in range(int(t.ti.min()), int(t.ti.max()) + 1):
        for b in range(int(t.tj.min()), int(t.tj.max()) + 1):
            s = t[(t.ti >= a) & (t.ti < a + w) & (t.tj >= b) & (t.tj < b + w)]
            if len(s) < n:
                continue
            o = np.lexsort((s.tile.values, np.abs(s.tk.values - s.tk.mean())))
            s = s.iloc[o[:n]]
            cov = float((s[JURKAT_GATED].sum() / tot).min())
            if best is None or cov > best[0]:
                best = (cov, a, b, s)
    cov, a, b, s = best
    ids = np.sort(s.tile.values.astype(np.int64))
    return ids, dict(a=a, b=b, tk=[int(s.tk.min()), int(s.tk.max())], n=len(ids), coverage_min=round(cov, 4),
                     coverage={r: round(float(s[r].sum() / tot[r]), 4) for r in JURKAT_GATED},
                     md5=hashlib.md5(ids.tobytes()).hexdigest(), md5_expected=SUBSET_MD5,
                     md5_as_frozen=hashlib.md5(ids.tobytes()).hexdigest() == SUBSET_MD5)


def prepare(out):
    os.makedirs(out, exist_ok=True)
    if os.path.exists(os.path.join(out, "inputs.json")):
        return jload(os.path.join(out, "inputs.json"))
    inp = resolve_inputs()
    os.makedirs(os.path.join(out, "jurkat"), exist_ok=True)
    ids, info = jurkat_subset(os.path.join(inp["jurkat"]["build"], "pack", "tiles.parquet"))
    np.save(os.path.join(out, "jurkat", "subset_tiles.npy"), ids)
    inp["jurkat"]["subset"] = info
    inp["snapshot"] = open(os.path.join(REPO, "SNAPSHOT_SHA")).read().strip() if os.path.exists(os.path.join(REPO, "SNAPSHOT_SHA")) else None
    inp["host"] = os.uname().nodename
    import shutil
    for cell, v in inp.items():                       # library caches of the input build's validation (same library)
        if not (isinstance(v, dict) and "build" in v):
            continue
        cd = os.path.join(out, cell, "cache"); os.makedirs(cd, exist_ok=True)
        for name in ("lib_all.npz", "lib_packed.npz"):
            hits = sorted(glob.glob(os.path.join(v["build"], "validate", "**", name), recursive=True), key=len)
            if hits and not os.path.exists(os.path.join(cd, name)):
                shutil.copy2(hits[0], os.path.join(cd, name)); v.setdefault("seeded_cache", []).append(hits[0])
    jsave(os.path.join(out, "inputs.json"), inp)
    log(out, "inputs:", json.dumps({c: (v["kind"], v["build"]) for c, v in inp.items() if isinstance(v, dict) and "build" in v}),
        "| jurkat subset md5", info["md5"], "(as frozen)" if info["md5_as_frozen"] else "(DIFFERS from the frozen md5: the rule decides)")
    return inp


# ------------------------------------------------------------------ one run (criteria A1, A3)
def _rename_legacy(pc):
    """pre-A0c key names inside presets (config.legacy_keys renames only top-level pack keys)"""
    for p in (pc.get("presets") or {}).values():
        if isinstance(p, dict) and "prepass_index" in p and "preplaced_index" not in p:
            p["preplaced_index"] = p.pop("prepass_index")


def run_config(base_config, setting):
    from wholecell.config import load_raw
    from wholecell.pack import switches
    raw = load_raw(base_config)
    pc = raw.setdefault("pack", {})
    _rename_legacy(pc)
    pinned = [n for n in switches.preset_names(pc.get("preset")) if n not in (switches.DEFAULT_PRESET, "legacy_v7")]
    pc.setdefault("presets", {})
    if setting != "R":
        pc["presets"][setting] = copy.deepcopy(switches.B1_BENCHMARK[setting])
        pc["preset"] = pinned + [setting]
    else:
        pc["preset"] = pinned or "none"
    raw.setdefault("build", {})["allow_dirty"] = True
    return raw


def link_inputs(src, rdir):
    os.makedirs(os.path.join(rdir, "pack"), exist_ok=True)
    for e in sorted(os.listdir(src)):
        if e in SKIP_TOP or e.startswith("compare") or e.endswith((".nohup", ".status")):
            continue
        d = os.path.join(rdir, e)
        if not os.path.lexists(d):
            os.symlink(os.path.join(src, e), d)
    for e in sorted(os.listdir(os.path.join(src, "pack"))):
        if e in PACK_OUTPUTS or ".tmp" in e:
            continue
        d = os.path.join(rdir, "pack", "preplaced_" + e[len("prepass_"):] if e.startswith("prepass_") else e)   # pre-A0c names
        if not os.path.lexists(d):
            os.symlink(os.path.join(src, "pack", e), d)


def _env(out):
    env = dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=REPO, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMBA_CACHE_DIR=os.path.join(out, "numba_cache"))
    return env


def _step(out, rdir, cfgp, step, extra=(), timeout=None, numba_threads="3"):
    env = _env(out)
    if step != "pack":
        env["NUMBA_NUM_THREADS"] = numba_threads
    cmd = [PY, "-m", "wholecell", "pack", "--cell", cfgp, "--build", rdir, "--step", step, *extra]
    t0 = time.time()
    with open(os.path.join(rdir, "logs", f"{step}.log"), "w") as lf:
        try:
            rc = subprocess.run(cmd, cwd=REPO, env=env, stdout=lf, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            rc = "timeout"
    return rc, time.time() - t0


def timeout_for(out, cell, inp):
    if cell == "jurkat":
        s1 = jload(os.path.join(out, "jurkat", "b1_s1", "run.json"), {})
        return max(1800.0, 2.2 * s1["pack_wall_s"]) if s1.get("pack_wall_s") else None
    w = inp["yeast"].get("pack_wall_s")
    return 2.2 * w if w else 36000.0


def run(out, cell, setting, workers, inp):
    rdir = os.path.join(out, cell, setting)
    rj = os.path.join(rdir, "run.json")
    if jload(rj, {}).get("status") in ("ok", "timeout", "failed", "skipped"):
        return jload(rj)
    src = inp[cell]["build"]
    if setting == "b1_s4p" and not glob.glob(os.path.join(src, "library", "pack_shapes", "pack_trees_*.npz")):
        info = dict(status="skipped", reason="no library/pack_shapes in the input build (B3 or later needed)")
        os.makedirs(rdir, exist_ok=True); jsave(rj, info); log(out, cell, setting, "SKIPPED:", info["reason"])
        return info
    if os.path.exists(os.path.join(rdir, "pack")):           # an interrupted run of this benchmark: start it again
        import shutil
        shutil.rmtree(os.path.join(rdir, "pack"))
    os.makedirs(os.path.join(rdir, "logs"), exist_ok=True)
    link_inputs(src, rdir)
    import yaml
    cfgp = os.path.join(rdir, "cell.yaml")
    with open(cfgp, "w") as f:
        f.write(f"# B1 benchmark run {setting} on {src} (config {inp[cell]['base_config']}, extends resolved)\n")
        yaml.safe_dump(run_config(inp[cell]["base_config"], setting), f, sort_keys=False, width=160)
    info = dict(setting=setting, cell=cell, input=src, workers=workers, host=os.uname().nodename, started=time.strftime("%FT%T"))
    # the counts step depends on the algorithm keys only through the shell cap (and the shapes): reuse the counts of an
    # earlier run of this benchmark with the same key (criteria A4: the counts time is reported, not part of the wall time)
    from wholecell.pack import switches
    key = lambda s: (bool(switches.B1_BENCHMARK[s]["shell_cap"]["enabled"]), (switches.B1_BENCHMARK[s].get("shapes") or {}).get("source"))
    donor = next((s for s in ORDER[cell] if s != setting and jload(os.path.join(out, cell, s, "run.json"), {}).get("status") == "ok"
                  and key(s) == key(setting) and os.path.exists(os.path.join(out, cell, s, "pack", "counts.npz"))), None)
    if donor:
        import shutil
        for f in ("counts.npz", "libs.json"):
            shutil.copy2(os.path.join(out, cell, donor, "pack", f), os.path.join(rdir, "pack", f))
        info.update(counts_from=donor, counts_s=jload(os.path.join(out, cell, donor, "run.json"))["counts_s"])
        rc = 0
    else:
        rc, dt = _step(out, rdir, cfgp, "counts")
        info["counts_s"] = round(dt, 1)
    if rc != 0:
        info.update(status="failed", failed_step="counts", rc=rc); jsave(rj, info); log(out, cell, setting, "FAILED counts"); return info
    extra = [f"workers={workers}"]
    if cell == "jurkat":
        extra.append(f"tiles={os.path.join(out, 'jurkat', 'subset_tiles.npy')}")
    to = timeout_for(out, cell, inp)
    info["timeout_s"] = to
    rc, dt = _step(out, rdir, cfgp, "pack", extra, timeout=to)
    info["pack_wall_s"] = round(dt, 1)
    if rc != 0:
        info.update(status="timeout" if rc == "timeout" else "failed", failed_step="pack", rc=rc)
        jsave(rj, info); log(out, cell, setting, info["status"].upper(), f"pack after {dt:.0f} s"); return info
    rc, dt2 = _step(out, rdir, cfgp, "summary")
    info["summary_s"] = round(dt2, 1)
    lt = os.path.join(rdir, "pack", "pack_log.tsv")
    if os.path.exists(lt):
        a = np.loadtxt(lt, ndmin=2)
        info["tile_cpu_s"] = float(a[:, 5].sum()) if len(a) else 0.0
    info.update(status="ok", finished=time.strftime("%FT%T"))
    jsave(rj, info)
    log(out, cell, setting, f"packed: wall {info['pack_wall_s']:.0f} s (counts {info['counts_s']:.0f} s)")
    return info


# ------------------------------------------------------------------ metrics (criteria A4)
def _masses(cfg, bdir, L, out, cell, extra):
    from wholecell.validate import common
    if L.get("mass_da"):
        return np.asarray(L["mass_da"], float), "libs.json mass_da"
    common.CACHE_DIR = os.path.join(out, cell, "cache")
    m, src = common.lib_masses(common.Inputs(cfg, bdir, **extra), bdir, nproc=3)
    return m, src


def _extra_inputs(cell, bdir, cfg):
    """core-step inputs of the yeast layout (as plugins/yeast/validate_adapter.core_inputs)"""
    if cell != "yeast":
        return {}
    ex = {}
    if os.path.isdir(os.path.join(bdir, "engine", "lib")):
        ex["lib_dir"] = os.path.join(bdir, "engine", "lib")
    from wholecell.pack import switches
    sh = switches.effective(cfg.get("pack") or {}).get("shapes") or {}
    f = sh.get("file") if isinstance(sh, dict) else None
    if f:
        ex["pack_trees"] = f if os.path.isabs(f) else os.path.join(bdir, f)
    pj = os.path.join(bdir, "geometry", "package.json")
    if os.path.exists(pj):
        ex["geom"] = json.load(open(pj))["package_dir"]
    return ex


def _seam(pack_dir, tiles_sel, T, lo0, tp):
    from wholecell.validate.suite import seam_counts, seam_index
    nb = ni = nt = 0
    for _, r in tp[tp.tile.isin(tiles_sel)].iterrows():
        f = os.path.join(pack_dir, "tiles", f"{int(r['tile'])}.npz")
        if not os.path.exists(f):
            continue
        ijk = np.array([int(r.ti), int(r.tj), int(r.tk)])
        b, i = seam_counts(np.load(f)["pos"].astype(np.float64) - (lo0 + ijk * T), T)
        nb += b; ni += i; nt += 1
    return dict(index=seam_index(nb, ni, T, nt) if ni else None, n_tiles=nt, n_band=int(nb), n_interior=int(ni))


def metrics(out, cell, setting, inp):
    import pandas as pd
    from wholecell.config import load
    from wholecell.validate import common, placement
    from wholecell.pack import switches
    rdir = os.path.join(out, cell, setting)
    res_f = os.path.join(rdir, "results.json")
    if os.path.exists(res_f):
        return jload(res_f)
    os.makedirs(rdir, exist_ok=True)
    if setting == "R":
        bdir = inp[cell]["build"]
        cfgp = os.path.join(rdir, "cell.yaml")
        if not os.path.exists(cfgp):
            import yaml
            with open(cfgp, "w") as f:
                yaml.safe_dump(run_config(inp[cell]["base_config"], "R"), f, sort_keys=False, width=160)
        ri = dict(status="ok", pack_wall_s=inp[cell].get("pack_wall_s"), input=bdir, note="input build's own pack, not re-packed")
        jsave(os.path.join(rdir, "run.json"), ri)
    else:
        bdir = rdir
        cfgp = os.path.join(rdir, "cell.yaml")
        ri = jload(os.path.join(rdir, "run.json"), {})
        if ri.get("status") != "ok":
            return dict(setting=setting, cell=cell, status=ri.get("status", "missing"), run=ri)
    cfg = load(cfgp)
    common.CACHE_DIR = os.path.join(out, cell, "cache")
    os.makedirs(common.CACHE_DIR, exist_ok=True)
    pd_ = os.path.join(bdir, "pack")
    regs = [r["name"] for r in cfg["pack"]["regions"]]
    kinds = {r["name"]: r["kind"] for r in cfg["pack"]["regions"]}
    tp = pd.read_parquet(os.path.join(pd_, "tiles.parquet")); meta = json.load(open(os.path.join(pd_, "tiles_meta.json")))
    T = float(meta["TILE"]); lo0 = np.asarray(meta["lo0"], float)
    sel = np.load(os.path.join(out, "jurkat", "subset_tiles.npy")) if cell == "jurkat" else tp.tile.values.astype(np.int64)
    L = json.load(open(os.path.join(pd_, "libs.json")))
    extra = _extra_inputs(cell, bdir, cfg)
    mass, msrc = _masses(cfg, bdir, L, out, cell, extra)
    nreg = len(regs)
    req = np.zeros(nreg); pla = np.zeros(nreg); n_inst = 0; packed = []
    for t in sel:
        f = os.path.join(pd_, "tiles", f"{int(t)}.npz")
        if not os.path.exists(f):
            continue
        z = np.load(f); packed.append(int(t)); n_inst += len(z["types"])
        np.add.at(req, z["job_region"].astype(np.int64), z["job_n"] * mass[z["job_lib"]])
        np.add.at(pla, z["job_region"].astype(np.int64), z["job_placed"] * mass[z["job_lib"]])
    cz = np.load(os.path.join(pd_, "counts.npz"))
    k = np.isin(cz["tile"], np.asarray(packed))
    unth = np.zeros(nreg)
    np.add.at(unth, cz["region"][k].astype(np.int64), cz["n"][k] * mass[cz["lib"][k]])
    vol_um3 = {r: float(tp[r].sum()) / 1e9 for r in regs if r in tp.columns}
    if cell == "jurkat":
        gated = [r for r in JURKAT_GATED if r in regs]
    else:
        gated = [r for r in regs if kinds[r] == "volume" and vol_um3.get(r, 0) >= 1.0 and unth[regs.index(r)] > 0]
    fr = {r: dict(placed_frac=float(pla[i] / req[i]) if req[i] > 0 else None,
                  placed_of_unthinned=float(pla[i] / unth[i]) if unth[i] > 0 else None,
                  request_mass_pg=float(req[i] / 6.022e23 * 1e12), unthinned_mass_pg=float(unth[i] / 6.022e23 * 1e12),
                  volume_um3=vol_um3.get(r), gated=r in gated, kind=kinds[r]) for i, r in enumerate(regs) if unth[i] > 0 or req[i] > 0}
    worst = min(((fr[r]["placed_frac"], r) for r in gated if fr.get(r) and fr[r]["placed_frac"] is not None), default=(None, None))
    res = dict(setting=setting, cell=cell, status="ok", input=bdir, tiles_packed=len(packed), instances=n_inst, mass_source=msrc,
               placed_mass_fraction=float(pla.sum() / req.sum()), placed_of_unthinned=float(pla.sum() / unth.sum()),
               worst_region=worst[1], worst_region_fraction=worst[0], regions=fr, gated=gated,
               thinning=(jload(os.path.join(pd_, "thinning.json"), {}) or {}).get("factors"),
               run=ri, preset=switches.effective(cfg["pack"]).get("_preset"))
    # seam (criteria A4)
    if cell == "jurkat":
        S = set(map(tuple, tp[tp.tile.isin(sel)][["ti", "tj", "tk"]].values.tolist()))
        sub = tp[tp.tile.isin(sel)]
        inner = sub[[all((i + a, j + b, kk + c) in S for a, b, c in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)))
                     for i, j, kk in zip(sub.ti, sub.tj, sub.tk)]]
        q = inner.tile[(inner["nucleoplasm"] / T ** 3) >= 0.95].values
        res["seam"] = dict(region="nucleoplasm", min_frac=0.95, **_seam(pd_, q, T, lo0, tp))
        q2 = inner.tile[(inner["cytosol"] / T ** 3) >= 0.8].values
        res["seam_cytosol_report"] = dict(region="cytosol", min_frac=0.8, **_seam(pd_, q2, T, lo0, tp))
    else:
        sc = ((cfg.get("validate") or {}).get("suite") or {}).get("seam") or {}
        reg_ = sc.get("region", "cytosol"); mf = float(sc.get("min_frac", 0.95)); n = int(sc.get("n_tiles", 300))
        cand = tp[(tp[reg_] / T ** 3) >= mf]
        rng = np.random.default_rng(int(sc.get("seed", 0)))
        if len(cand) > n:
            cand = cand.iloc[np.sort(rng.choice(len(cand), n, replace=False))]
        res["seam"] = dict(region=reg_, min_frac=mf, n_qualifying=int(((tp[reg_] / T ** 3) >= mf).sum()),
                           **_seam(pd_, cand.tile.values, T, lo0, tp))
    # C7(a) / C8(a) (core placement step, criteria A4)
    if cell == "jurkat":
        ptiles = np.sort(np.random.default_rng(0).choice(np.asarray(packed), min(300, len(packed)), replace=False))
    else:
        ptiles = None
    variants = [("", {})]
    if setting == "b1_s4p":
        pt = sorted(glob.glob(os.path.join(bdir, "library", "pack_shapes", "pack_trees_*.npz")))
        if pt:
            variants.append(("_pack_trees", {"pack_trees": pt[-1]}))
    for suf, ex2 in variants:
        inp_ = common.Inputs(cfg, bdir, **{**extra, **ex2})
        pj = os.path.join(rdir, f"placement{suf}.json")
        P = jload(pj)
        if P is None:
            t0 = time.time()
            placement.run(cfg, bdir, inp_, pj, os.path.join(rdir, f"placement{suf}_tiles.tsv"), overwrite=True, nproc=int(os.environ.get("B1_METRIC_PROCS", "8")), tiles=ptiles)
            P = jload(pj); P["_seconds_b1"] = time.time() - t0; jsave(pj, P)
        npk = P["n_packed"]
        res[f"c7{suf}"] = dict(clash_rate_per_1e6=P["clash_rate_per_1e6"],
                               packed_packed_per_1e6=P["clash_pairs"]["packed_packed"] / npk * 1e6 if npk else None,
                               packed_prepass_per_1e6=P["clash_pairs"]["packed_prepass"] / npk * 1e6 if npk else None,
                               max_clash_depth_nm=P["max_clash_depth_nm"], n_tiles=P["n_tiles_sampled"], n_packed=npk,
                               collision_shapes=str((P.get("settings") or {}).get("collision_shapes", "library sphere trees")))
        res[f"c8{suf}"] = dict(containment_rate_per_1e6=P["containment_rate_per_1e6"])
    jsave(res_f, res)
    log(out, cell, setting, f"metrics: placed {res['placed_mass_fraction']:.4f} (of unthinned {res['placed_of_unthinned']:.4f}), "
        f"worst {res['worst_region']} {res['worst_region_fraction']}, seam {res['seam'].get('index')}, "
        f"C7a pp {res['c7']['packed_packed_per_1e6']}, C8a {res['c8']['containment_rate_per_1e6']}")
    return res


# ------------------------------------------------------------------ selection (criteria A5)
def _tol(v):
    return max(0.1 * abs(v), 1.0)


def conditions(res, s1, ref, ref_wall):
    """criteria A5 conditions 1-5 for one run -> {name: (pass, value, limit)}"""
    c = {}
    bad = {r: res["regions"][r]["placed_frac"] for r in res["gated"]
           if res["regions"].get(r) and res["regions"][r]["placed_frac"] is not None and res["regions"][r]["placed_frac"] < 0.98}
    c["1_regions_ge_0.98"] = (not bad, bad or "all gated >= 0.98", 0.98)
    si = res["seam"].get("index")
    c["2_seam"] = (si is not None and abs(si - 1) <= 0.1, si, "|index - 1| <= 0.1")
    lim = s1["placed_of_unthinned"] - 0.003 if s1 else None
    c["3_placed_of_unthinned"] = (lim is None or res["placed_of_unthinned"] >= lim, res["placed_of_unthinned"], lim)
    pp = res["c7"]["packed_packed_per_1e6"]
    ok4 = pp is not None and pp <= 1.0
    lims = {}
    if ref and ref.get("status") == "ok":
        for k, v in (("clash_rate_per_1e6", res["c7"]["clash_rate_per_1e6"]), ("containment_rate_per_1e6", res["c8"]["containment_rate_per_1e6"])):
            rv = ref["c7" if k.startswith("clash") else "c8"][k]
            if rv is None:
                continue
            lims[k] = (v, rv + _tol(rv)); ok4 = ok4 and v is not None and v <= rv + _tol(rv)
    c["4_clash_containment"] = (ok4, dict(packed_packed_per_1e6=pp, **{k: v[0] for k, v in lims.items()}),
                                dict(packed_packed_per_1e6=1.0, **{k: v[1] for k, v in lims.items()}))
    w = res["run"].get("pack_wall_s")
    c["5_wall"] = (w is not None and ref_wall is not None and w <= 2 * ref_wall, w, 2 * ref_wall if ref_wall else None)
    return c


def select(out, inp):
    from wholecell.pack import switches
    allres = {}
    for cell in ("jurkat", "yeast"):
        for s in ORDER[cell] + ["R"]:
            f = os.path.join(out, cell, s, "results.json")
            r = jload(f) or dict(status=(jload(os.path.join(out, cell, s, "run.json"), {}) or {}).get("status", "missing"))
            allres.setdefault(cell, {})[s] = r
    table = {}
    for cell, rs in allres.items():
        s1 = rs.get("b1_s1") if rs.get("b1_s1", {}).get("status") == "ok" else None
        ref = rs.get("R")
        ref_wall = (s1 or {}).get("run", {}).get("pack_wall_s") if cell == "jurkat" else inp["yeast"].get("pack_wall_s")
        s1w = (s1 or {}).get("run", {}).get("pack_wall_s")
        for s, r in rs.items():
            if s == "R":
                continue
            if r.get("status") != "ok":
                table.setdefault(cell, {})[s] = dict(eligible=False, status=r.get("status"), n_failed=5, cost=None)
                continue
            c = conditions(r, s1, ref, ref_wall)
            w = r["run"].get("pack_wall_s")
            table.setdefault(cell, {})[s] = dict(eligible=all(v[0] for v in c.values()), status="ok",
                                                 n_failed=sum(not v[0] for v in c.values()),
                                                 conditions={k: dict(passed=v[0], value=v[1], limit=v[2]) for k, v in c.items()},
                                                 wall_s=w, cost_cell=(w / s1w) if (w and s1w) else None)
    cand = []
    for s in COMMON_CANDIDATES:
        tj, ty = table.get("jurkat", {}).get(s, {}), table.get("yeast", {}).get(s, {})
        cost = (tj.get("cost_cell") or 0) + (ty.get("cost_cell") or 0) if tj.get("cost_cell") and ty.get("cost_cell") else None
        cand.append(dict(setting=s, eligible_both=bool(tj.get("eligible") and ty.get("eligible")), cost=cost,
                         n_failed=(tj.get("n_failed", 5) or 0) + (ty.get("n_failed", 5) or 0)))
    elig = [c for c in cand if c["eligible_both"] and c["cost"] is not None]
    if elig:
        chosen = min(elig, key=lambda c: c["cost"])["setting"]; verdict = "SELECTED"
    else:
        ok = [c for c in cand if c["cost"] is not None] or cand
        chosen = min(ok, key=lambda c: (c["n_failed"], c["cost"] if c["cost"] is not None else 1e9))["setting"]
        verdict = "NO ELIGIBLE SETTING"
    yeast_only_s5 = bool(table.get("yeast", {}).get("b1_s5", {}).get("eligible")) and \
        not any(table.get("yeast", {}).get(s, {}).get("eligible") for s in COMMON_CANDIDATES)
    # packing shapes (criteria A5, P8): Jurkat S4p vs S4
    shapes = dict(decision="library", reason="")
    p, q = allres.get("jurkat", {}).get("b1_s4p", {}), allres.get("jurkat", {}).get("b1_s4", {})
    if p.get("status") == "ok" and q.get("status") == "ok":
        why = []
        if not (p["c7"]["packed_packed_per_1e6"] is not None and p["c7"]["packed_packed_per_1e6"] <= 1.0):
            why.append(f"C7(a) packed-packed on library trees {p['c7']['packed_packed_per_1e6']} > 1")
        if p["placed_of_unthinned"] < q["placed_of_unthinned"] - 0.003:
            why.append("placed_of_unthinned below S4 - 0.003")
        for r in q["gated"]:
            a_, b_ = p["regions"].get(r, {}).get("placed_frac"), q["regions"].get(r, {}).get("placed_frac")
            if a_ is not None and b_ is not None and a_ < b_ - 0.003:
                why.append(f"{r} {a_:.4f} < S4 {b_:.4f} - 0.003")
        if (p["run"].get("pack_wall_s") or 1e18) > (q["run"].get("pack_wall_s") or 0):
            why.append("wall time above S4")
        shapes = dict(decision="library" if why else "pack_shapes", reason="; ".join(why) or "all S4p conditions met")
    else:
        shapes["reason"] = f"S4p {p.get('status')} / S4 {q.get('status')}: not decidable, library kept"
    cd = copy.deepcopy(switches.B1_BENCHMARK[chosen])
    cd.pop("shapes", None)
    cd["budgets"] = {"per_region": {}}          # a None in a preset removes the key; as a default: no per-region budgets
    if shapes["decision"] == "pack_shapes":
        cd["shapes"] = {"source": "pack_shapes"}
    sel = dict(verdict=verdict, chosen=chosen, candidates=cand, table=table, shapes=shapes, yeast_only_full_squeeze=yeast_only_s5,
               common_defaults=cd if verdict == "SELECTED" else None, proposal=cd if verdict != "SELECTED" else None,
               apply=("paste common_defaults into wholecell/pack/switches.py COMMON_DEFAULTS (criteria Part C); yeast: remove the "
                      "legacy_x5 pin" if verdict == "SELECTED" else "nothing applied: Pedro decides (proposal = the setting with "
                      "the fewest failed conditions)"),
               inputs={c: dict(build=v["build"], kind=v["kind"], pack_wall_s=v.get("pack_wall_s")) for c, v in inp.items()
                       if isinstance(v, dict) and "build" in v},
               criteria="docs/unification_criteria/B1_packer.md", written=time.strftime("%FT%T"))
    jsave(os.path.join(out, "selection.json"), sel)
    # flat table
    rows = ["cell\tsetting\tstatus\tpack_wall_s\tplaced_mass_fraction\tplaced_of_unthinned\tworst_region\tworst_fraction\tseam\t"
            "c7a_packed_packed_per_1e6\tc7a_clash_rate_per_1e6\tc8a_per_1e6\teligible\tn_failed"]
    for cell, rs in allres.items():
        for s, r in rs.items():
            t = table.get(cell, {}).get(s, {})
            if r.get("status") == "ok" and "placed_mass_fraction" in r:
                rows.append("\t".join(str(x) for x in (cell, s, "ok", r["run"].get("pack_wall_s"), round(r["placed_mass_fraction"], 5),
                                                         round(r["placed_of_unthinned"], 5), r["worst_region"], r["worst_region_fraction"],
                                                         r["seam"].get("index"), r["c7"]["packed_packed_per_1e6"], r["c7"]["clash_rate_per_1e6"],
                                                         r["c8"]["containment_rate_per_1e6"], t.get("eligible", ""), t.get("n_failed", ""))))
            else:
                rows.append(f"{cell}\t{s}\t{r.get('status')}" + "\t" * 12)
    with open(os.path.join(out, "results.tsv"), "w") as f:
        f.write("\n".join(rows) + "\n")
    jsave(os.path.join(out, "results.json"), allres)
    log(out, f"selection: {verdict}: {chosen}; shapes {shapes['decision']} ({shapes['reason']})")
    return sel


# ------------------------------------------------------------------ CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["all", "prepare", "run", "metrics", "select"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--cells", default="jurkat,yeast")
    ap.add_argument("--cell"); ap.add_argument("--setting")
    a = ap.parse_args()
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)
    inp = prepare(out)
    if a.cmd == "prepare":
        return
    if a.cmd == "run":
        run(out, a.cell, a.setting, a.workers, inp); return
    if a.cmd == "metrics":
        metrics(out, a.cell, a.setting, inp); return
    if a.cmd == "select":
        select(out, inp); return
    for cell in a.cells.split(","):
        for s in ORDER[cell]:
            try:
                r = run(out, cell, s, a.workers, inp)
                if r.get("status") == "ok":
                    _metrics_sub(out, cell, s)
            except Exception as e:                       # recorded; the benchmark goes on (criteria A3)
                log(out, cell, s, "ERROR", repr(e))
    for cell in a.cells.split(","):
        try:
            _metrics_sub(out, cell, "R")
        except Exception as e:
            log(out, cell, "R", "ERROR", repr(e))
    select(out, inp)


def _metrics_sub(out, cell, s):
    """metrics in a subprocess (memory of the placement step is returned after each run)"""
    with open(os.path.join(out, cell, s, "logs_metrics.log") if os.path.isdir(os.path.join(out, cell, s)) else os.devnull, "a") as lf:
        rc = subprocess.run([PY, os.path.abspath(__file__), "metrics", "--out", out, "--cell", cell, "--setting", s], cwd=REPO,
                            env=dict(_env(out), NUMBA_NUM_THREADS="3"), stdout=lf, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        log(out, cell, s, f"metrics FAILED rc={rc}")


if __name__ == "__main__":
    main()
