"""Build discipline: isolated build directories, build stamp, frozen-input verification, write-once outputs.

Layout of one build (spec: docs/reference/cell_config.md, section "build"):
  <results_root>/<cell_id>/<build_id>/
     BUILD_STAMP.json         commit, dirty flag, config sha256, master seed (+ legacy seeds), host, env, started/ended
     config.resolved.json     the resolved config used
     frozen_inputs.tsv        name, path, expected checksum, actual checksum, ok
     <stage>/...              one folder per stage; a stage writes only inside its own folder

build_id = <YYYYMMDD>-<commit7>[-dirty]-s<master_seed>[-<tag>]
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import socket
import subprocess
import sys

from .config import REPO


def git_info(repo: str = REPO) -> dict:
    def g(*a):
        try:
            return subprocess.run(["git", "-C", repo, *a], capture_output=True, text=True, timeout=30).stdout.strip()
        except Exception:
            return ""
    commit = g("rev-parse", "HEAD")
    dirty = bool(g("status", "--porcelain", "--untracked-files=no"))
    return {"commit": commit, "commit7": commit[:7], "dirty": dirty}


def make_build_id(cfg, tag: str | None = None, date: str | None = None) -> str:
    gi = git_info()
    date = date or _dt.date.today().strftime("%Y%m%d")
    s = int(cfg.get("seeds", {}).get("master", 0))
    bid = f"{date}-{gi['commit7'] or 'nogit'}{'-dirty' if gi['dirty'] else ''}-s{s}"
    return bid + (f"-{tag}" if tag else "")


def build_dir(cfg, build_id: str) -> str:
    root = cfg["_site"].get("results_root") or os.path.join(REPO, "results")
    return os.path.join(root, cfg["cell_id"], build_id)


def file_checksum(path: str, algo: str = "md5", chunk: int = 1 << 24) -> str:
    h = hashlib.new(algo)
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def visible_names(folder: str) -> list[str]:
    """sorted entry names of a folder without hidden ones (".*": fetch markers, editor and OS files are not data)"""
    return sorted(n for n in os.listdir(folder) if not n.startswith("."))


def listing_md5(folder: str) -> str:
    """md5 of the sorted visible file-name listing of a folder (the pdb_listing_md5 of the assembly inputs)"""
    return hashlib.md5("\n".join(visible_names(folder)).encode()).hexdigest()


def _listing_checks(cfg):
    """(frozen input name, expected listing md5) of every config block that declares pdb_assemblies + pdb_listing_md5"""
    out, stack = [], [dict(cfg)]
    while stack:
        d = stack.pop()
        if isinstance(d, dict):
            if isinstance(d.get("pdb_assemblies"), str) and d.get("pdb_listing_md5"):
                out.append((d["pdb_assemblies"], str(d["pdb_listing_md5"])))
            stack.extend(v for k, v in d.items() if k != "frozen_inputs")
        elif isinstance(d, list):
            stack.extend(d)
    return sorted(set(out))


def _plugin_checks(cfg):
    """rows of the plugin steps that register a frozen_check (same name as the step's `run`): file-by-file checks of
    input trees that a plugin step otherwise makes only inside a build (yeast: the pipeline data read set)"""
    from . import registry
    registry.load_plugins(cfg)
    rows = []
    for st in cfg.get("plugin_steps") or []:
        if ("frozen_check", st.get("run")) in registry._REG:
            rows += registry.get("frozen_check", st["run"])(cfg, st.get("params") or {})
    return rows


def verify_frozen(cfg, out_tsv: str | None = None, strict: bool = True) -> list[dict]:
    """check every declared frozen input against its checksum (md5 or sha256 key). Directories: checksum of a
    manifest file named in `manifest` (if given) else existence only. Folders with a declared listing md5
    (pdb_listing_md5 of the assembly inputs) are also checked by their visible file-name listing, the check
    inventory_assembly makes, so a changed folder fails here and not hours into a build; plugin steps that register a
    frozen_check add their file-by-file checks (yeast: the data tree the engine reads)."""
    rows = []
    for name, spec in sorted((cfg.get("frozen_inputs") or {}).items()):
        p = spec["path"]
        algo = "sha256" if "sha256" in spec else "md5"
        exp = spec.get(algo)
        if os.path.isdir(p):
            act = file_checksum(spec["manifest"], algo) if spec.get("manifest") else "dir"
        elif os.path.exists(p):
            act = file_checksum(p, algo) if exp else "unchecked"
        else:
            act = "MISSING"
        ok = act != "MISSING" and (exp is None or exp == act or act == "dir")
        rows.append(dict(name=name, path=p, algo=algo, expected=exp or "", actual=act, ok=ok))
    rows += _plugin_checks(cfg)
    fi = cfg.get("frozen_inputs") or {}
    for name, exp in _listing_checks(cfg):
        p = (fi.get(name) or {}).get("path", name)
        act = listing_md5(p) if os.path.isdir(p) else "MISSING"
        rows.append(dict(name=f"{name}:listing", path=p, algo="md5", expected=exp, actual=act, ok=act == exp))
    if out_tsv:
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        with open(out_tsv, "w") as f:
            f.write("name\tpath\talgo\texpected\tactual\tok\n")
            for r in rows:
                f.write("\t".join(str(r[k]) for k in ("name", "path", "algo", "expected", "actual", "ok")) + "\n")
    bad = [r for r in rows if not r["ok"]]
    if bad and strict:
        raise RuntimeError("frozen inputs failed verification: " + ", ".join(f"{r['name']} ({r['actual']})" for r in bad))
    return rows


def write_stamp(cfg, bdir: str, extra: dict | None = None) -> dict:
    """write BUILD_STAMP.json + config.resolved.json (refuses a dirty tree unless cfg build.allow_dirty)"""
    gi = git_info()
    if gi["dirty"] and not cfg.get("build", {}).get("allow_dirty", False):
        raise RuntimeError("working tree is dirty; commit first or set build.allow_dirty: true (recorded in the stamp)")
    os.makedirs(bdir, exist_ok=True)
    stamp = dict(
        cell_id=cfg["cell_id"], build_dir=bdir, **gi,
        config_source=cfg.get("_source"), config_sha256=cfg.config_hash(),
        seeds=cfg.get("seeds", {}), pythonhashseed=os.environ.get("PYTHONHASHSEED"),
        host=socket.gethostname(), platform=platform.platform(), python=sys.version.split()[0],
        executable=sys.executable, started=_dt.datetime.now().isoformat(timespec="seconds"))
    if extra:
        stamp.update(extra)
    p = os.path.join(bdir, "BUILD_STAMP.json")
    if os.path.exists(p):
        old = json.load(open(p))
        if old.get("config_sha256") != stamp["config_sha256"] or old.get("commit") != stamp["commit"]:
            raise RuntimeError(f"{p} exists for a different commit/config; use a new build id (no in-place rebuilds)")
        return old
    json.dump(stamp, open(p, "w"), indent=1)
    json.dump({k: v for k, v in cfg.items() if k != "_site"} | {"_site": cfg["_site"]},
              open(os.path.join(bdir, "config.resolved.json"), "w"), indent=1, default=str)
    return stamp


def check_env():
    """PYTHONHASHSEED must be 0 for every build process (set by the Snakefile / launcher)"""
    if os.environ.get("PYTHONHASHSEED") != "0":
        raise RuntimeError("PYTHONHASHSEED must be set to 0 for build steps")


def write_once(path: str):
    """guard for stage outputs: raise if the file exists (stages never rewrite outputs in place)"""
    if os.path.exists(path):
        raise FileExistsError(f"{path} exists; stage outputs are write-once (new build id for a rebuild)")
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    return path


def atomic_save(path: str, saver, *args, **kw):
    """write to a hidden temp file in the same folder, then os.replace (a killed job never leaves partial files)"""
    d, b = os.path.split(path)
    tmp = os.path.join(d, f".{b}.tmp{os.getpid()}")
    saver(tmp, *args, **kw)
    if b.endswith(".npz") and not tmp.endswith(".npz") and os.path.exists(tmp + ".npz"):
        tmp = tmp + ".npz"
    os.replace(tmp, path)


# ---------------------------------------------------------------------------------------------- legacy names (A0c)
# Stage names follow the paper's Methods (unification plan s8). The pre-placed stage was called "prepass" until A0c
# (2026-11-17): builds made before that have <build>/prepass/ and pack/prepass_*.npy. New builds use the names on the
# left. LEGACY_NAMES is the ONE place that knows the old names; it is used only to READ builds made before the rename
# (reference builds, "vs previous" CI comparisons): ref_path() for files, scripts/ci_compare_builds.py for paths and
# text. Each entry maps one path component (folder or file name) of a new build to its name in an old build.
LEGACY_NAMES = {
    "preplaced": "prepass",                                       # stage folder <build>/preplaced/
    "preplaced_manifest.json": "prepass_manifest.json",           # preplaced/: what the pack index reads
    "preplaced.npz": "prepass.npz",                               # preplaced/serial/: all arrays of the serial run
    "expected_after_preplaced.parquet": "expected_after_prepass.parquet",   # preplaced/serial/
    "preplaced_types.npy": "prepass_types.npy",                   # pack/: pre-placed index
    "preplaced_pos.npy": "prepass_pos.npy",
    "preplaced_quat.npy": "prepass_quat.npy",
    "preplaced_off.npy": "prepass_off.npy",
    "preplaced_libs.json": "prepass_libs.json",
    "preplaced_groups.json": "prepass_groups.json",
    "lib_preplaced.npz": "lib_prepass.npz",                       # validate/cache/
}


def legacy_rel(rel: str) -> str:
    """a build-relative path in the names of builds made before A0c (each component through LEGACY_NAMES)"""
    parts = rel.replace(os.sep, "/").split("/")
    return os.path.join(*[LEGACY_NAMES.get(p, p) for p in parts]) if rel else rel


def ref_path(base: str, rel: str) -> str:
    """<base>/<rel> in the current names; if that does not exist but the pre-A0c name does (a reference build made
    before the rename), the old path. Use it wherever a reader may be pointed at an old build."""
    p = os.path.join(base, rel)
    if not os.path.exists(p):
        q = os.path.join(base, legacy_rel(rel))
        if os.path.exists(q):
            return q
    return p
