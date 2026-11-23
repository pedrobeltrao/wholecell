"""Yeast pipeline engine `yeast_x5` and packer `yeast_rsa_v1` (registered in plugins/yeast/__init__.py).

The yeast stage scripts (vendored, plugins/yeast/engine = git archive e3b7448 + PORT_PATCHES.md) run unchanged in
their own environment (config pipeline_engine.python_env). This module reproduces what the yeast
scripts/build_isolated.sh + run_all.sh did for X5, step by step, inside a wholecell build directory:

  <build>/engine/work/      YC_WORK (seeded from the frozen work_seed: RCSB cache, pool pairs, sweeps, image tables)
  <build>/engine/lib/       YC_LIB  (starts as the 342 frozen X5 start entries; every other entry is rebuilt)
  <build>/engine/pdb/       YC_PDB  (link farm to the yeast bundle work/pdb, read-only; new downloads land here)
  <build>/engine/numba_cache/  NUMBA_CACHE_DIR (fresh per build; run_all.sh cleared the in-tree caches)
  <code>/.venv -> python_env, <code>/pipeline/data -> the yeast bundle pipeline/data (read in place, read-only)

Environment of every step = build_isolated.sh: YC_WORK, YC_LIB, YC_PDB, YC_INPUTS, YC_SEED, YC_NPROC,
OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1, PYTHONHASHSEED=0, YC_GPU_LOCK (per build), plus YC_NGPU and
CUDA_VISIBLE_DEVICES for the GPU validation steps, YC_GEOM_PACKAGE for stage 2.
Nothing is written, run or locked in the yeast tree (its files are read-only for this user in any case).
"""
from __future__ import annotations

import glob
import json
import os
import re
import shutil
import subprocess
import sys
import time

STAGE_SCRIPTS = {"s0_library", "s1_inventory", "s2_geometry", "s3_pack", "s4_validate", "s5_encode"}


def _ec(cfg):
    return cfg["pipeline_engine"]


def paths(cfg, bdir):
    e = _ec(cfg)
    code = e["code"]
    if not os.path.isabs(code):
        code = os.path.join(cfg["_site"]["repo"], code)
    eb = os.path.join(bdir, "engine")
    return dict(code=code, engine=eb, work=os.path.join(eb, "work"), lib=os.path.join(eb, "lib"),
                pdb=os.path.join(eb, "pdb"), numba=os.path.join(eb, "numba_cache"), done=os.path.join(eb, "_done"),
                py=os.path.join(code, ".venv", "bin", "python"), frozen=e["frozen_dir"])


def nproc(cfg):
    return int(os.environ.get("WHOLECELL_ENGINE_NPROC") or _ec(cfg).get("nproc") or 8)


def env(cfg, bdir, extra=None):
    e = _ec(cfg); P = paths(cfg, bdir)
    out = dict(os.environ)
    for k in list(out):                       # never inherit another build's settings
        if k.startswith("YC_") or k.startswith("NUMBA_"):
            out.pop(k)
    out.update(YC_WORK=P["work"], YC_LIB=P["lib"], YC_PDB=P["pdb"], YC_INPUTS=e["inputs_dir"], YC_SEED=str(int(e.get("seed", 0))),
               YC_NPROC=str(nproc(cfg)), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0",
               YC_GPU_LOCK=os.path.join(P["engine"], "gpu.lock"), YC_NGPU=str(int(e.get("gpus", 1))),
               CUDA_VISIBLE_DEVICES=",".join(str(i) for i in range(int(e.get("gpus", 1)))),
               NUMBA_CACHE_DIR=P["numba"], PYTHONPATH="", PYTHONDONTWRITEBYTECODE="")
    out.pop("PYTHONPATH")
    if extra:
        out.update(extra)
    return out


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sh(cmd, cwd, environ, tee=None):
    """run a shell command (bash, pipefail); stream to stdout and optionally tee into a file (as run_all.sh)"""
    full = f"set -o pipefail; {cmd}" + (f" 2>&1 | tee {tee}" if tee else "")
    log("$", cmd, f"(cwd {cwd})")
    r = subprocess.run(["bash", "-c", full], cwd=cwd, env=environ)
    if r.returncode != 0:
        raise SystemExit(f"engine step failed (exit {r.returncode}): {cmd}")


# ------------------------------------------------------------------ manifests of frozen inputs
def read_manifest(path):
    rows = []
    for line in open(path):
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("  ", 3)                          # make_manifests.sh: md5  size  mtime  ./relpath
        if parts[-1] == "" or not parts[0] or len(parts[0]) != 32:
            # file not readable by this user when the manifest was made: md5sum failed, printf shifted the fields
            size, mtime, rel = [p for p in parts if p != ""][:3]
            md5 = None
        else:
            md5, size, mtime, rel = parts
        rows.append((md5, int(size), int(float(mtime)), rel[2:] if rel.startswith("./") else rel))
    return rows


def _md5(p):
    import hashlib
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def verify_manifest(root, manifest, full_md5=False):
    """(n files, [(relpath, problem)]) of `root` against a make_manifests.sh manifest. full_md5: md5 of every file;
    else size + mtime, and a file whose mtime differs is accepted when its md5 matches the manifest (a downloaded or
    copied tree has new mtimes); problems: missing | size | md5 | mtime (mtime differs and the manifest has no md5)"""
    bad = []
    rows = read_manifest(manifest)
    for md5, size, mtime, rel in rows:
        p = os.path.join(root, rel)
        if not os.path.exists(p):
            bad.append((rel, "missing")); continue
        st = os.stat(p)
        if st.st_size != size:
            bad.append((rel, "size")); continue
        if (full_md5 or int(st.st_mtime) != mtime) and md5 is not None:
            if _md5(p) != md5:
                bad.append((rel, "md5"))
        elif int(st.st_mtime) != mtime:
            bad.append((rel, "mtime"))
    return len(rows), bad


# ------------------------------------------------------------------ steps
def step_setup(cfg, bdir, kw):
    e = _ec(cfg); P = paths(cfg, bdir)
    for d in ("engine", "work", "lib", "pdb", "numba", "done"):
        os.makedirs(P[d], exist_ok=True)
    fz = P["frozen"]
    # 1. frozen inputs: start library + work seed (md5 of every file), read-in-place trees (size + mtime vs manifest)
    checks = {}
    for name, root, full in (("lib_start", os.path.join(fz, "lib_start"), True), ("work_seed", os.path.join(fz, "work_seed"), True),
                             ("pipeline_data", e["data_dir"], False), ("pdb", e["pdb_dir"], False),
                             ("afdb_dimers", os.path.join(e["inputs_dir"], "afdb_yeast_dimers"), False)):
        man = os.path.join(fz, f"manifest_{name}.tsv")
        n, bad = verify_manifest(root, man, full_md5=full)
        # files of the read-in-place yeast data changed after X5 started are listed (not read by the build) but not fatal
        changed_ok = [b for b in bad if name == "pipeline_data" and b[0].startswith("validation/literature/endocytosis/")]
        hard = [b for b in bad if b not in changed_ok]
        checks[name] = dict(root=root, manifest=man, files=n, mismatches=bad, fatal=hard, level="md5" if full else "size+mtime (md5 when mtime differs)")
        log(f"frozen {name}: {n} files, {len(bad)} mismatches ({len(hard)} fatal)")
        if hard:
            raise SystemExit(f"frozen input {name} does not match its manifest: {hard[:10]}")
    # 2. code links (untracked, inside the build's own worktree)
    venv = os.path.join(P["code"], ".venv")
    if not os.path.lexists(venv):
        os.symlink(e["python_env"], venv)
    data = os.path.join(P["code"], "pipeline", "data")
    if not os.path.lexists(data):
        os.symlink(e["data_dir"], data)
    # 3. library start set (copy, as build_isolated.sh copied the bundle library) and the PDB link farm
    for f in sorted(os.listdir(os.path.join(fz, "lib_start"))):
        dst = os.path.join(P["lib"], f)
        if not os.path.exists(dst):
            shutil.copy2(os.path.join(fz, "lib_start", f), dst)
    for f in sorted(os.listdir(e["pdb_dir"])):
        dst = os.path.join(P["pdb"], f)
        if not os.path.lexists(dst):
            os.symlink(os.path.join(e["pdb_dir"], f), dst)
    # 4. work seed (build_isolated.sh: rcsb cache, pool pairs, sweeps, image/proximity tables, pair summary)
    ws = os.path.join(fz, "work_seed")
    for root, dirs, fs in os.walk(ws):
        dirs.sort()
        rel = os.path.relpath(root, ws)
        os.makedirs(os.path.join(P["work"], rel), exist_ok=True)
        for f in sorted(fs):
            dst = os.path.join(P["work"], rel, f)
            if os.path.exists(dst):
                continue
            if rel == "interactions":                                   # linked, as build_isolated.sh did
                os.symlink(os.path.join(root, f), dst)
            else:
                shutil.copy2(os.path.join(root, f), dst)
    # 5. yeast build stamp (run_all.sh: build_stamp.py start <first>)
    sh(f"{P['py']} build_stamp.py start {int(e.get('first_stage', 0))}", os.path.join(P["code"], "scripts"), env(cfg, bdir))
    pip = subprocess.run([P["py"], "-m", "pip", "list", "--format=freeze"], capture_output=True, text=True, env=env(cfg, bdir)).stdout
    json.dump(dict(checks=checks, nproc=nproc(cfg), python=P["py"], python_env=e["python_env"],
                   python_version=subprocess.run([P["py"], "-c", "import sys; print(sys.version)"], capture_output=True, text=True).stdout.strip(),
                   packages=pip.splitlines(), host=os.uname().nodename),
              open(os.path.join(P["engine"], "setup.json"), "w"), indent=1)


def step_stage(cfg, bdir, stage, extra_env=None):
    """one yeast stage script, exactly as run_all.sh: cd scripts; $PY sN.py 2>&1 | tee $WORK/sN.log"""
    P = paths(cfg, bdir)
    E = env(cfg, bdir, extra_env)
    sh(f"{P['py']} {stage}.py", os.path.join(P["code"], "scripts"), E, tee=os.path.join(P["work"], f"{stage}.log"))
    if stage == "s1_inventory":       # run_all.sh: refresh the per-protein feature table after the inventory
        sh(f"{P['py']} features.py > {os.path.join(P['work'], 'features.log')} 2>&1", os.path.join(P["code"], "pipeline"), E)


def step_geometry_package(cfg, bdir, kw):
    P = paths(cfg, bdir); g = cfg.get("geometry") or {}
    out = os.path.join(bdir, "geometry", "package.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tk = os.path.join(cfg["_site"]["repo"], "toolkits", "analytic", "yeast.py")
    sh(f"{P['py']} {tk} --engine-code {P['code']} --dest-root {g['package_root']} --out {out} --log {os.path.join(P['work'], 'geometry_package.log')}",
       os.path.join(P["code"], "scripts"), env(cfg, bdir))


def geom_package_dir(bdir):
    return json.load(open(os.path.join(bdir, "geometry", "package.json")))["package_dir"]


VAC_MARK = "# port (plugins/yeast/runner.py): vacuole lines of the geometry generator, copied from"
VAC_LINE = re.compile(r"^(?:\d\d:\d\d:\d\d)?\s*vacuole")      # generator lines that start with "vacuole" (after the time)


def append_generator_vacuole_lines(work, bdir):
    """X5 wrote the vacuole lines ("vacuole: N compartment(s), V um3, surface S um2", ...) into s2_geometry.log, which
    validation/vacuole_surface.py parses. In the port the geometry is made by the toolkit step (geometry_package), so
    those lines are in work/geometry_package.log (fallback <build>/logs/engine_geometry_package.log): append them to
    work/s2_geometry.log after a marker line (port defect fixed 2026-11-17). Returns the number of lines appended."""
    dst = os.path.join(work, "s2_geometry.log")
    src = next((f for f in (os.path.join(work, "geometry_package.log"), os.path.join(bdir, "logs", "engine_geometry_package.log"))
                if os.path.exists(f)), None)
    if src is None:
        log("no geometry generator log found: s2_geometry.log has no vacuole lines (vacuole_surface.py will fail)")
        return 0
    lines = [l if l.endswith("\n") else l + "\n" for l in open(src, errors="replace")
             if VAC_LINE.match(l)]
    if os.path.exists(dst) and VAC_MARK in open(dst, errors="replace").read():
        return 0                                         # already appended (step re-run without a fresh stage log)
    with open(dst, "a") as f:
        f.write(f"{VAC_MARK} {src}\n")
        f.writelines(lines)
    log(f"appended {len(lines)} generator vacuole lines from {src} to {dst}")
    return len(lines)


def step_s2(cfg, bdir, kw):
    step_stage(cfg, bdir, "s2_geometry", {"YC_GEOM_PACKAGE": geom_package_dir(bdir)})
    append_generator_vacuole_lines(paths(cfg, bdir)["work"], bdir)


def step_validate_all(cfg, bdir, kw):
    P = paths(cfg, bdir)
    sh(f"bash {os.path.join(P['code'], 'scripts', 'validate_all.sh')}", P["code"], env(cfg, bdir))


def step_core_copies(cfg, bdir, kw):
    """the shared core copy-number step on this config (no yeast code), for criterion I3"""
    from wholecell.inventory import copies
    from wholecell.inventory.context import StepContext
    os.makedirs(os.path.join(bdir, "inventory"), exist_ok=True)
    ctx = StepContext(cfg, bdir, "copies")
    T, S = copies.step(ctx)
    log("core copies step:", len(T), "proteins")


def step_pack_view(cfg, bdir, kw):
    P = paths(cfg, bdir)
    sh(f"{P['py']} {os.path.join(cfg['_site']['repo'], 'plugins', 'yeast', 'pack_view.py')} --bdir {bdir} --config-json {os.path.join(bdir, 'config.resolved.json')}",
       os.path.join(P["code"], "scripts"), env(cfg, bdir))


def step_compare_x5(cfg, bdir, kw):
    P = paths(cfg, bdir)
    sh(f"{P['py']} {os.path.join(cfg['_site']['repo'], 'plugins', 'yeast', 'repro', 'compare_x5.py')} --bdir {bdir} --ref {_ec(cfg)['reference']} "
       f"--nproc {min(nproc(cfg), 16)}", os.path.join(P["code"], "scripts"), env(cfg, bdir))


STEPS = {"setup": step_setup, "s0_library": lambda c, b, k: step_stage(c, b, "s0_library"),
         "s1_inventory": lambda c, b, k: step_stage(c, b, "s1_inventory"), "geometry_package": step_geometry_package,
         "s2_geometry": step_s2, "s4_validate": lambda c, b, k: step_stage(c, b, "s4_validate"),
         "s5_encode": lambda c, b, k: step_stage(c, b, "s5_encode"), "validate_all": step_validate_all,
         "core_copies": step_core_copies, "pack_view": step_pack_view, "compare_x5": step_compare_x5}


class YeastX5Engine:
    """registry object for kind pipeline_engine"""

    name = "yeast_x5"

    @staticmethod
    def steps(cfg):
        out = []
        for s in _ec(cfg)["steps"]:
            out.append(dict(id=s["id"], stage=s.get("stage", "engine"), step=s.get("step", s["id"]), after=list(s.get("after", [])),
                            gpu=s["id"] in ("validate_all",), heavy=s["id"] in ("pack", "s0_library", "s1_inventory", "s2_geometry", "s4_validate", "s5_encode")))
        return out

    @staticmethod
    def run(cfg, bdir, step, **kw):
        if step not in STEPS:
            raise SystemExit(f"yeast_x5: unknown engine step '{step}' (known: {sorted(STEPS)})")
        t0 = time.time()
        STEPS[step](cfg, bdir, kw)
        log(f"engine step {step} done in {time.time() - t0:.0f} s")


def pack_yeast_rsa(cfg, bdir, step=None, **kw):
    """pack engine yeast_rsa_v1 = yeast stage 3 (serial pre-pass builders, organisation layer, pilot thinning, 6 size-class
    rounds x 8 checkerboard phases, squeeze); outputs in <build>/engine/work/yc3d (tiles_raw, prepass.npz, ...)"""
    if step not in (None, "pack", "all"):
        raise SystemExit(f"pack engine yeast_rsa_v1 has one step 'pack' (got {step})")
    t0 = time.time()
    step_stage(cfg, bdir, "s3_pack")
    log(f"yeast_rsa_v1 pack done in {time.time() - t0:.0f} s")
