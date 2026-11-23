"""Yeast plugin steps of the core stage list (config plugin_steps; registry kind `plugin_step`; unification A9).

The yeast cell builds through the core stages. Some core steps still call vendored yeast modules (plugins/yeast/engine,
never edited: assembly plugin steps and hooks, the pre-placed context and builders, the yeast validation scripts).
Those modules read and write a yeast-style file layout (YC_ROOT, YC_WORK, YC_LIB, YC_PDB). Inside a core build that
layout is <build>/engine/ (config: inventory.params.assembly_steps.engine, anchor `yeast_engine`). The steps here fill it
from core outputs (links, copies, format views: never the other way round) or run vendored yeast scripts there:

  yeast.engine_setup       engine/{root,work,pdb}: data links, work seed copy, PDB link farm, manifests re-checked,
                           the yeast build stamp (work/build.json)
  yeast.protein_table      core copies + core localisation in the engine's stage-1 protein table layout (the columns
                           the core does not make yet come from the vendored inventory.build(); copies and
                           comp_weights must agree with it, else the step fails)
  yeast.library_species    species of the base library (vendored ingredients.build_ingredients, as s0_library.py)
  yeast.engine_lib         engine/lib = writable copy of the sealed library/lib
  yeast.engine_view        engine/work/yc3d/{geometry.json, geometry_fields.npz} (package native files),
                           {tile_measure, tile_labvol}.parquet (core tiles, yeast names), engine/work/results/inventory
                           (core assembly tables), work/s1_inventory.log (the core inventory rule logs)
  yeast.engine_pack_view   engine/work/yc3d after the core pack: tiles_raw/<tile>.npz (engine tile layout: core tile
                           arrays + empty in-tile mRNA arrays; a link when the core tile has them), geom -> pack/tile_geom,
                           prepass.npz, lib_ids.npy,
                           expected_totals.parquet (wc3d.shell_capacity of the pre-placed totals, as stage 3 saved it),
                           tile_counts.parquet, pack_summary.parquet, thinning.json, mrna_counts.npy;
                           work/s3_pack.log (the core pre-placed + pack rule logs)
  yeast.engine_script      one vendored script (features.py, s4_validate.py, s5_encode.py) in the core environment
  yeast.validate_all       scripts/validate_all.sh (GPU steps under a per-build lock)
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time

from wholecell.registry import register

from .assembly import engine_env, resolve_layout


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def layout(bdir, params):
    E = resolve_layout(params.get("engine") or {}, bdir)
    for k in ("root", "work", "lib", "pdb", "code"):
        if not E.get(k):
            raise SystemExit(f"plugin step params.engine.{k} is required (config anchor yeast_engine)")
    return E


def environ(E, extra=None):
    """environment of a vendored yeast script (as plugins/yeast/runner.env, core interpreter)"""
    # NUMBA_NUM_THREADS / NUMBA_THREADING_LAYER (exported by the core since 2026-11-20, wholecell/numba_env.py) are not
    # passed on: the vendored scripts keep their own numba defaults, as before
    out = {k: v for k, v in os.environ.items() if not k.startswith(("YC_", "NUMBA_CACHE", "NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER"))}
    out.update(YC_ROOT=E["root"], YC_WORK=E["work"], YC_LIB=E["lib"], YC_PDB=E["pdb"], YC_SEED=str(int(E.get("seed", 0))),
               YC_NPROC=str(int(E.get("nproc", 8))), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", PYTHONHASHSEED="0",
               NUMBA_CACHE_DIR=os.path.join(E["work"], "numba_cache"), PYTHONDONTWRITEBYTECODE="1")
    if E.get("inputs"):
        out["YC_INPUTS"] = E["inputs"]
    out.pop("PYTHONPATH", None)
    out.update(extra or {})
    return out


def sh(cmd, cwd, env, logfile=None):
    full = f"set -o pipefail; {cmd}" + (f" 2>&1 | tee {logfile}" if logfile else "")
    log("$", cmd, f"(cwd {cwd})")
    r = subprocess.run(["bash", "-c", full], cwd=cwd, env=env)
    if r.returncode != 0:
        raise SystemExit(f"yeast step failed (exit {r.returncode}): {cmd}")


def _link(src, dst):
    if not os.path.lexists(dst):
        try:
            os.symlink(src, dst)
        except FileExistsError:
            pass


def _cat_logs(bdir, names, dst):
    """engine-style stage log (e.g. work/s3_pack.log, read by plugins/yeast/repro/compare_x5.py) = the core rule logs of
    this build, concatenated in order (missing logs skipped)"""
    with open(dst, "w") as out:
        for n in names:
            f = os.path.join(bdir, "logs", f"{n}.log")
            if os.path.exists(f):
                out.write(open(f, errors="replace").read())


def _jdump(obj, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    json.dump(obj, open(tmp, "w"), indent=1, default=str)
    os.replace(tmp, path)


# ------------------------------------------------------------------ inventory
def engine_setup(cfg, bdir, p):
    """origin: plugins/yeast/runner.py step_setup (X5 port), without the engine venv"""
    from .runner import verify_manifest
    E = layout(bdir, p)
    fi = cfg["frozen_inputs"]
    for k in ("root", "work", "pdb"):
        os.makedirs(E[k], exist_ok=True)
    checks = {}
    for name in p.get("check") or []:
        spec = fi[name]
        full = name == p.get("work_seed")
        n, bad = verify_manifest(spec["path"], spec["manifest"], full_md5=full)
        # files of the read-in-place yeast data changed after X5 started (not read by the build): listed, not fatal
        hard = [b for b in bad if not (name == p.get("data") and b[0].startswith("validation/literature/endocytosis/"))]
        checks[name] = dict(files=n, mismatches=bad[:50], n_mismatches=len(bad), fatal=hard[:50], level="md5" if full else "size+mtime")
        log(f"frozen {name}: {n} files, {len(bad)} mismatches ({len(hard)} fatal)")
        if hard:
            raise SystemExit(f"frozen input {name} does not match its manifest: {hard[:10]}")
    # root: the yeast pipeline data + the engine's non-code files (cell_meshes.npz, vacuole_template.npz, ...)
    _link(fi[p["data"]]["path"], os.path.join(E["root"], "data"))
    pipe = os.path.join(E["code"], "pipeline")
    for f in sorted(os.listdir(pipe)):
        if not f.endswith(".py") and os.path.isfile(os.path.join(pipe, f)):
            _link(os.path.join(pipe, f), os.path.join(E["root"], f))
    # PDB link farm (read in place; new downloads would land in the build)
    src = fi[p["pdb"]]["path"]
    for f in sorted(os.listdir(src)):
        _link(os.path.join(src, f), os.path.join(E["pdb"], f))
    # work seed (copied; interactions linked, as build_isolated.sh of the yeast tree did)
    ws = fi[p["work_seed"]]["path"]
    for root, dirs, fs in os.walk(ws):
        dirs.sort()
        rel = os.path.relpath(root, ws)
        os.makedirs(os.path.join(E["work"], rel), exist_ok=True)
        for f in sorted(fs):
            dst = os.path.join(E["work"], rel, f)
            if os.path.lexists(dst):
                continue
            if rel == "interactions":
                os.symlink(os.path.join(root, f), dst)
            else:
                shutil.copy2(os.path.join(root, f), dst)
    # yeast build stamp (run_all.sh: build_stamp.py start <first stage>): validation MANIFEST checks against it
    sh(f"{sys.executable} build_stamp.py start 0", os.path.join(E["code"], "scripts"), environ(E))
    _jdump(dict(layout=E, checks=checks, python=sys.executable, python_version=sys.version.split()[0],
                host=os.uname().nodename), os.path.join(bdir, "engine", "setup.json"))


def protein_table(cfg, bdir, p):
    """the stage-1 protein table of the assembly steps: rows, order and annotation columns of the vendored
    inventory.build(); copies_per_cell / ppm from the core copies step, comp_weights from the core localisation"""
    import numpy as np
    import pandas as pd
    E = layout(bdir, p)
    engine_env(E)
    import inventory as inv
    df, _geo = inv.build()
    C = pd.read_csv(os.path.join(bdir, "inventory", "copies", "copies.tsv"), sep="\t")
    V = pd.read_csv(os.path.join(bdir, "inventory", "v0", "inventory_v0.tsv"), sep="\t", keep_default_na=False)
    key = "orf" if "orf" in C.columns else C.columns[0]
    cc = [c for c in ("copies_per_cell", "copies") if c in C.columns][0]
    core_copies = C.set_index(key)[cc].astype(float)
    core_ppm = C.set_index(key)["ppm"].astype(float) if "ppm" in C.columns else None

    def parse(s):
        return {k: float(v) for k, v in (kv.rsplit(":", 1) for kv in s.split(";"))} if s else {}
    vkey = "orf" if "orf" in V.columns else V.columns[0]          # localisation evidence key (config genes.key)
    core_w = {g: parse(s) for g, s in zip(V[vkey], V["compartments"])}
    tol = float(p.get("tolerance_rel", 1e-9))
    rep = dict(rows=len(df), tolerance_rel=tol)
    a = df.orf.map(core_copies)
    if a.isna().any():
        raise SystemExit(f"core copies miss {int(a.isna().sum())} ORFs of the protein table (e.g. {list(df.orf[a.isna()][:5])})")
    d = (a - df.copies_per_cell).abs() / df.copies_per_cell.abs().clip(lower=1e-12)
    rep["copies_max_rel"] = float(d.max())
    rep["copies_identical"] = int((a.values == df.copies_per_cell.values).sum())
    if core_ppm is not None:
        b = df.orf.map(core_ppm)
        rep["ppm_max_rel"] = float(((b - df.ppm).abs() / df.ppm.abs().clip(lower=1e-12)).max())
    out_w, wmax, nkey = [], 0.0, 0
    for orf, s in zip(df.orf, df.comp_weights):
        e = json.loads(s)
        m = core_w.get(orf)
        if m is None:
            raise SystemExit(f"core localisation has no row for ORF {orf}")
        if set(e) != set(m):
            nkey += 1
            out_w.append(json.dumps({k: round(v, 4) for k, v in m.items()}))
            continue
        wmax = max(wmax, max(abs(e[k] - m[k]) for k in e))
        out_w.append(json.dumps({k: round(m[k], 4) for k in e}))          # the engine's key order, core values
    rep.update(comp_weights_max_abs=wmax, comp_weights_key_set_differences=nkey)
    # the localisation must equal the engine's only under the X5 preset; with another preset (B5: common) the core
    # comp_weights are the point, the differences are reported (copies must still agree)
    preset = ((cfg.get("inventory") or {}).get("params") or {}).get("localisation", {}).get("preset")
    rep["localisation_preset"] = preset
    rep["localisation_must_match_engine"] = strict = preset in (None, "legacy_x5")
    ok = rep["copies_max_rel"] <= tol and (not strict or (wmax <= 1e-6 and nkey == 0))
    rep["pass"] = bool(ok)
    out = df.copy()
    out["copies_per_cell"] = a.values
    if core_ppm is not None:
        out["ppm"] = df.orf.map(core_ppm).values
    out["comp_weights"] = out_w
    out = out.sort_values("copies_per_cell", ascending=False, kind="stable").reset_index(drop=True)
    rep["order_unchanged"] = bool(list(out.orf) == list(df.orf))
    od = os.path.join(bdir, "inventory", "yeast")
    os.makedirs(od, exist_ok=True)
    _jdump(rep, os.path.join(od, "protein_table_check.json"))
    log("protein table:", json.dumps(rep))
    if not ok:
        raise SystemExit("core copies / localisation differ from the engine's stage-1 table beyond tolerance "
                         "(inventory/yeast/protein_table_check.json)")
    out.attrs = {k: (json.loads(json.dumps(v, default=float))) for k, v in df.attrs.items()}
    tmp = os.path.join(od, f"protein_table.parquet.tmp{os.getpid()}")
    out.to_parquet(tmp)
    os.replace(tmp, os.path.join(od, "protein_table.parquet"))
    log(f"{len(out)} proteins, {float(np.sum(out.copies_per_cell)):.0f} copies -> inventory/yeast/protein_table.parquet")


# ------------------------------------------------------------------ library
def library_species(cfg, bdir, p):
    """origin: plugins/yeast/library_species.py (A3) on the protein table of this build"""
    import pandas as pd
    E = layout(bdir, p)
    engine_env(E)
    import ingredients as ingm
    df = pd.read_parquet(os.path.join(bdir, "inventory", "yeast", "protein_table.parquet"))
    ing, _ = ingm.build_ingredients(df)
    from .library_species import species_table
    out = species_table(ing)
    dst = os.path.join(bdir, "inventory", "yeast", "library_species.tsv")
    out.to_csv(dst + ".tmp", sep="\t", index=False)
    os.replace(dst + ".tmp", dst)
    log(f"{len(out)} species ({(out.kind == 'protein').sum()} protein, {(out.kind == 'complex').sum()} complex) -> {dst}")


def engine_lib(cfg, bdir, p):
    E = layout(bdir, p)
    src = os.path.join(bdir, "library", "lib")
    if os.path.exists(E["lib"]) and os.listdir(E["lib"]):
        raise SystemExit(f"{E['lib']} exists and is not empty (a build's engine library is made once)")
    os.makedirs(E["lib"], exist_ok=True)
    n = 0
    for f in sorted(os.listdir(src)):
        s = os.path.join(src, f)
        if os.path.isfile(s):
            shutil.copy2(s, os.path.join(E["lib"], f))
            os.chmod(os.path.join(E["lib"], f), 0o644)            # sealed files are read-only; engine stages rewrite some
            n += 1
    af = os.path.join(bdir, "library", "library_af.json")            # the engine's results/library_af.json (s0 wrote it)
    if os.path.exists(af):
        os.makedirs(os.path.join(E["work"], "results"), exist_ok=True)
        shutil.copy2(af, os.path.join(E["work"], "results", "library_af.json"))
    _jdump(dict(source=src, files=n, note="writable copy of the sealed library; assembly hooks and pre-placed models add "
                                          "entries"), os.path.join(bdir, "engine", "lib.json"))
    log(f"engine/lib: {n} files copied from library/lib")


# ------------------------------------------------------------------ pre-placed (stage-2 view)
def _tile_names(bdir):
    import numpy as np
    import pandas as pd
    meta = json.load(open(os.path.join(bdir, "pack", "tiles_meta.json")))
    M = pd.read_parquet(os.path.join(bdir, "pack", "tiles.parquet"))
    off = np.round(np.asarray(meta["lo0"]) / meta["TILE"]).astype(np.int64)
    ijk = M[["ti", "tj", "tk"]].values + off
    return M, meta, [f"t_{int(a):+03d}_{int(b):+03d}_{int(c):+03d}" for a, b, c in ijk]


def engine_view(cfg, bdir, p):
    import pandas as pd
    E = layout(bdir, p)
    engine_env(E)
    import cell3d as c3
    y = os.path.join(E["work"], "yc3d")
    os.makedirs(y, exist_ok=True)
    pkg = json.load(open(os.path.join(bdir, "geometry", "package.json")))["package_dir"]
    for f in ("geometry.json", "geometry_fields.npz"):                # the package's engine-native files (read-only)
        _link(os.path.join(pkg, "native", f), os.path.join(y, f))
    M, meta, names = _tile_names(bdir)
    A = M[meta["columns"]].set_axis(names, axis=0)
    A.to_parquet(os.path.join(y, "tile_measure.parquet"))
    LV = pd.read_parquet(os.path.join(bdir, "pack", "tile_labvol.parquet"))
    if LV.shape[1] != len(c3.LABEL_NAMES):
        raise SystemExit(f"tile_labvol: {LV.shape[1]} label columns, the engine has {len(c3.LABEL_NAMES)}")
    LV = LV.set_axis(list(c3.LABEL_NAMES), axis=1).set_axis(names, axis=0)   # same label order (package = cell3d codes)
    LV.to_parquet(os.path.join(y, "tile_labvol.parquet"))
    # stage-1 tables (results/inventory) = the core assembly outputs (copies: engine scripts must not write into them)
    src = os.path.join(bdir, "inventory", "assembly")
    dst = os.path.join(E["work"], "results", "inventory")
    os.makedirs(dst, exist_ok=True)
    n = 0
    for f in sorted(os.listdir(src)):
        if f.endswith((".csv", ".json")) and f != "assembly_steps.json":
            shutil.copy2(os.path.join(src, f), os.path.join(dst, f)); n += 1
    _cat_logs(bdir, ["inventory_copies", "inventory_v0", "inventory_yeast_protein_table", "inventory_assembly"],
              os.path.join(E["work"], "s1_inventory.log"))
    log(f"engine view: package {pkg}; {len(A)} tiles; {n} assembly tables -> {dst}")


# ------------------------------------------------------------------ validation (pack view)
# keys of an engine raw tile (wc3d.pack_one np.savez); mrna_pts / mrna_off = in-tile mRNA threads, empty since the
# threads are pre-placed (s3_pack: "no in-tile threads any more"; the pre-placed threads are in prepass.npz)
ENGINE_TILE_KEYS = ("types", "pos", "quat", "job_lib", "job_region", "job_n", "job_placed", "mrna_pts", "mrna_off", "overflow")


def _engine_tile(z, dst):
    """core pack tile -> engine raw tile file (same instance arrays; empty in-tile mRNA unless the core tile has them)"""
    import numpy as np
    d = {k: z[k] for k in z.files}
    d.setdefault("mrna_pts", np.zeros((0, 3), np.float32))
    d.setdefault("mrna_off", np.zeros(1, np.int64))
    tmp = dst + ".tmp.npz"
    if os.path.lexists(dst):
        os.remove(dst)
    np.savez(tmp, **d)
    os.replace(tmp, dst)


def engine_pack_view(cfg, bdir, p):
    import numpy as np
    import pandas as pd
    E = layout(bdir, p)
    engine_env(E)
    import wc3d
    y = os.path.join(E["work"], "yc3d")
    raw = os.path.join(y, "tiles_raw")
    os.makedirs(raw, exist_ok=True)
    M, meta, names = _tile_names(bdir)
    ser = os.path.join(bdir, "preplaced", "serial")
    lib_ids = [str(x) for x in np.load(os.path.join(ser, "lib_ids.npy"))]
    libs = json.load(open(os.path.join(bdir, "pack", "libs.json")))["libs"]
    if [str(x) for x in libs] != lib_ids:
        raise SystemExit("pack/libs.json order differs from preplaced/serial/lib_ids.npy")
    rows = []
    for i, nm in enumerate(names):
        f = os.path.join(bdir, "pack", "tiles", f"{i}.npz")
        z = np.load(f)
        dst = os.path.join(raw, f"{nm}.npz")
        if all(k in z.files for k in ENGINE_TILE_KEYS):
            _link(os.path.relpath(f, raw), dst)
        else:                                  # engine tile layout (wc3d.pack_one): the in-tile mRNA arrays the vendored
            _engine_tile(z, dst)               # readers index (val3d.load_instances, webenc._stage2); core tiles lack them
        z = np.load(dst)
        att = z["att_stats"].astype(np.int64) if "att_stats" in z.files else np.zeros(2, np.int64)
        rows.append((i, int(z["job_n"].sum()), int(z["job_placed"].sum()), bool(z["overflow"]), 0.0,
                     int(len(z["mrna_off"]) - 1) if "mrna_off" in z.files else 0, int(att[0]), int(att[1])))
    pd.DataFrame(rows, columns=["ti", "requested", "placed", "overflow", "sec", "n_mrna", "attract_tried", "attract_placed"]
                 ).to_parquet(os.path.join(y, "pack_summary.parquet"))
    _link(os.path.relpath(os.path.join(ser, "preplaced.npz"), y), os.path.join(y, "prepass.npz"))
    # per-tile membrane meshes (config GDIR = yc3d/geom; s5_encode / webenc read them): the core pack writes the same
    # files (tiles3d layout, engine tile names) to pack/tile_geom
    tg = os.path.join(bdir, "pack", "tile_geom")
    if os.path.isdir(tg):
        _link(os.path.relpath(tg, y), os.path.join(y, "geom"))
    np.save(os.path.join(y, "lib_ids.npy"), np.array(lib_ids))
    np.save(os.path.join(y, "mrna_counts.npy"), np.zeros(len(names), np.int64))
    # requests: the thinned counts when the pilot thinned, else the counts (yeast tile_counts layout)
    cf = os.path.join(bdir, "pack", "counts_thinned.npz")
    cz = np.load(cf if os.path.exists(cf) else os.path.join(bdir, "pack", "counts.npz"))
    pd.DataFrame({k: cz[k] for k in ("tile", "lib", "region", "n")}).to_parquet(os.path.join(y, "tile_counts.parquet"))
    tj = os.path.join(bdir, "pack", "thinning.json")
    if os.path.exists(tj):
        _jdump(json.load(open(tj)).get("factors_exact") or {}, os.path.join(y, "thinning.json"))
    # expected totals as stage 3 saved them: the pre-placed totals through wc3d.shell_capacity
    exp3 = pd.read_parquet(os.path.join(ser, "expected_after_preplaced.parquet"))
    if wc3d.SHELL_CAP:
        measure = pd.read_parquet(os.path.join(y, "tile_measure.parquet"))
        lp = wc3d.build_libpack(lib_ids)
        exp3 = wc3d.shell_capacity(exp3, measure, lp, log=log)
    exp3.to_parquet(os.path.join(y, "expected_totals.parquet"))
    _cat_logs(bdir, ["preplaced_serial", "pack_tiles", "pack_organisation", "pack_counts", "pack_index", "pack_pack",
                     "pack_summary"], os.path.join(E["work"], "s3_pack.log"))
    log(f"engine pack view: {len(names)} tiles, {sum(r[2] for r in rows):,} placed of {sum(r[1] for r in rows):,}")


# ------------------------------------------------------------------ vendored scripts
def _python(p):
    """interpreter of a vendored script: params.python (e.g. the yeast engine env "{data_root}/conda_envs/yeast_x5/bin/python"
    for s5_encode / validate_all: cupy, fast_simplification, pyBigWig, playwright, markdown, openmm are not in the core
    env), else the core interpreter"""
    py = p.get("python") or sys.executable
    if not os.path.exists(py):
        raise SystemExit(f"plugin step params.python {py} does not exist")
    return py


def _code_links(E, py=None):
    """validate_all.sh runs $B/.venv/bin/python with B = the engine code dir and reads pipeline/data: links inside the
    build's own code worktree (as plugins/yeast/runner.step_setup did); the interpreter = params.python (else core)"""
    venv = os.path.join(E["code"], ".venv")
    prefix = os.path.dirname(os.path.dirname(py)) if py else sys.prefix
    if os.path.islink(venv) and os.path.realpath(venv) != os.path.realpath(prefix):
        os.remove(venv)
    _link(prefix, venv)
    _link(os.path.join(E["root"], "data"), os.path.join(E["code"], "pipeline", "data"))


def engine_script(cfg, bdir, p):
    E = layout(bdir, p)
    env = environ(E, {"YC_NPROC": str(int(p.get("nproc") or E.get("nproc", 8)))})
    cwd = os.path.join(E["code"], p["cwd"])
    sh(f"{_python(p)} {p['script']}", cwd, env, logfile=os.path.join(E["work"], p.get("log") or p["script"] + ".log"))


def validate_all(cfg, bdir, p):
    E = layout(bdir, p)
    _code_links(E, _python(p))
    n = int(p.get("gpus", 1))
    env = environ(E, {"YC_NPROC": str(int(p.get("nproc") or E.get("nproc", 8))), "YC_NGPU": str(n),
                      "CUDA_VISIBLE_DEVICES": ",".join(str(i) for i in range(n)),
                      "YC_GPU_LOCK": os.path.join(bdir, "engine", "gpu.lock")})
    sh(f"bash {os.path.join(E['code'], 'scripts', 'validate_all.sh')}", E["code"], env)


for _n, _f in (("yeast.engine_setup", engine_setup), ("yeast.protein_table", protein_table),
               ("yeast.library_species", library_species), ("yeast.engine_lib", engine_lib),
               ("yeast.engine_view", engine_view), ("yeast.engine_pack_view", engine_pack_view),
               ("yeast.engine_script", engine_script), ("yeast.validate_all", validate_all)):
    register("plugin_step", _n)(_f)
