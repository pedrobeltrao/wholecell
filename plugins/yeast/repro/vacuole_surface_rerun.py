"""Re-run validation/vacuole_surface.py on a finished X5-port build without touching its engine/work (2026-11-17).

    conda_envs/yeast_x5/bin/python plugins/yeast/repro/vacuole_surface_rerun.py --bdir <build>

Port defect: vacuole_surface.py (vendored, never edited) parses the generator's vacuole line from work/s2_geometry.log,
but in the port that line is in <build>/logs/engine_geometry_package.log (the toolkit step makes the geometry), so the
script failed in validate_all (fixed for future builds in plugins/yeast/runner.append_generator_vacuole_lines).
This script runs ONLY vacuole_surface.py, with the build's own _code engine and the YC_* environment of
plugins/yeast/runner.env, against a SHADOW work dir <build>/validate/vacuole_surface_rerun/work:
  s2_geometry.log                    = the build's engine/work/s2_geometry.log + marker line + the vacuole lines of
                                       <build>/logs/engine_geometry_package.log (a real file)
  yc3d/expected_totals.parquet       -> symlink to the build's file
  results/inventory/ingredients.csv  -> symlink to the build's file
  validation/vacuole_surface.json    the output
The figure that the script writes next to its code (<code>/docs/figures) is redirected to <rerun>/docs/figures.
Provenance (inputs with md5, command, environment, python) -> <rerun>/provenance.json. Refuses to run if <rerun> exists.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

INPUTS = ["s2_geometry.log", "yc3d/expected_totals.parquet", "results/inventory/ingredients.csv"]
MARK = "# vacuole_surface_rerun (plugins/yeast/repro/vacuole_surface_rerun.py): vacuole lines of the geometry generator, copied from"
VAC_LINE = re.compile(r"^(?:\d\d:\d\d:\d\d)?\s*vacuole")
WRAP = """import importlib.util, os, sys
spec = importlib.util.spec_from_file_location("vacuole_surface", sys.argv[1])
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
m.HERE = sys.argv[2]          # figure -> <rerun>/docs/figures (the script writes os.path.join(HERE, '..', 'docs', 'figures'))
m.main()
"""


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bdir", required=True)
    ap.add_argument("--nproc", type=int, default=3)
    a = ap.parse_args()
    b = os.path.abspath(a.bdir)
    cfg = json.load(open(os.path.join(b, "config.resolved.json")))
    e = cfg["pipeline_engine"]
    code = e["code"]
    work = os.path.join(b, "engine", "work")
    gen = os.path.join(b, "logs", "engine_geometry_package.log")
    out = os.path.join(b, "validate", "vacuole_surface_rerun")
    if os.path.exists(out):
        raise SystemExit(f"{out} exists (never overwritten)")
    sw = os.path.join(out, "work")
    for d in ("yc3d", "results/inventory", "validation"):
        os.makedirs(os.path.join(sw, d))
    vac = [l if l.endswith("\n") else l + "\n" for l in open(gen, errors="replace") if VAC_LINE.match(l)]
    with open(os.path.join(sw, "s2_geometry.log"), "w") as f:
        f.write(open(os.path.join(work, "s2_geometry.log"), errors="replace").read())
        f.write(f"{MARK} {gen}\n")
        f.writelines(vac)
    for rel in INPUTS[1:]:
        os.symlink(os.path.join(work, rel), os.path.join(sw, rel))
    script = os.path.join(code, "validation", "vacuole_surface.py")
    wrap = os.path.join(out, "run_vacuole_surface.py")
    open(wrap, "w").write(WRAP)
    py = os.path.join(code, ".venv", "bin", "python")
    env = {k: v for k, v in os.environ.items() if not k.startswith(("YC_", "NUMBA_"))}
    env.pop("PYTHONPATH", None)
    yc = dict(YC_WORK=sw, YC_LIB=os.path.join(b, "engine", "lib"), YC_PDB=os.path.join(b, "engine", "pdb"), YC_INPUTS=e["inputs_dir"],
              YC_SEED=str(int(e.get("seed", 0))), YC_NPROC=str(a.nproc), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
              PYTHONHASHSEED="0", NUMBA_CACHE_DIR=os.path.join(out, "numba_cache"), MPLCONFIGDIR=os.path.join(out, "mplconfig"),
              PYTHONDONTWRITEBYTECODE="1")
    env.update(yc)
    cmd = [py, wrap, script, os.path.join(out, "code_here")]
    cwd = os.path.join(code, "validation")
    t0 = time.time()
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True)
    open(os.path.join(out, "run.log"), "w").write(r.stdout + r.stderr)
    res = os.path.join(sw, "validation", "vacuole_surface.json")
    prov = dict(
        purpose="re-run of validation/vacuole_surface.py only (port defect: generator vacuole line not in s2_geometry.log); "
                "the build's engine/work is not written",
        date=time.strftime("%Y-%m-%d %H:%M:%S"), build=b, script=script, script_md5=md5(script), python=py,
        python_version=subprocess.run([py, "-c", "import sys; print(sys.version)"], capture_output=True, text=True).stdout.strip(),
        command=" ".join(cmd), cwd=cwd, env=yc, returncode=r.returncode, wall_s=round(time.time() - t0, 1),
        shadow_work=sw,
        inputs={rel: dict(source=os.path.join(work, rel), md5=md5(os.path.join(work, rel)),
                          shadow="file (build log + appended vacuole lines)" if rel == "s2_geometry.log" else "symlink")
                for rel in INPUTS},
        generator_log=dict(path=gen, md5=md5(gen), vacuole_lines=[l.rstrip("\n") for l in vac]),
        shadow_s2_geometry_log_md5=md5(os.path.join(sw, "s2_geometry.log")),
        other_inputs_read=dict(chen2025_s8=os.path.join(code, "pipeline", "data", "validation", "literature", "chen2025_sxt", "mbc-36-ar132-s008.xlsx"),
                               meshes=os.path.join(e["inputs_dir"], "yeast_organelle_meshes.zip")),
        output=res, output_md5=md5(res) if os.path.exists(res) else None,
        figure=os.path.join(out, "docs", "figures", "vacuole_surface.png"))
    json.dump(prov, open(os.path.join(out, "provenance.json"), "w"), indent=1)
    print(json.dumps(dict(build=os.path.basename(b), returncode=r.returncode, output=res if os.path.exists(res) else None)))
    if r.returncode:
        print(r.stderr[-3000:])
        raise SystemExit(1)


if __name__ == "__main__":
    main()
