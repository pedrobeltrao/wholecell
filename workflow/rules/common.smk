# shared setup: the plan of this cell (computed in the build environment by `python -m wholecell.plan`), build dir,
# command template. The Snakemake env is orchestration only; no wholecell code is imported here.
import json
import os
import subprocess

import yaml

REPO = workflow.basedir
if "cell" not in config:
    raise ValueError("snakemake --config cell=<cell> is required")
CELL = config["cell"]
_site_file = os.environ.get("WHOLECELL_SITE") or os.path.join(REPO, "config", "site.yaml")
_SITE = yaml.safe_load(open(_site_file))
PY = os.environ.get("WHOLECELL_PYTHON") or _SITE["python"]
_cmd = [PY, "-m", "wholecell.plan", "--cell", CELL] + (["--build-tag", str(config["build_tag"])] if config.get("build_tag") else [])
_res = subprocess.run(_cmd, cwd=REPO, capture_output=True, text=True,
                      env=dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=REPO))
_line = [l for l in _res.stdout.splitlines() if l.startswith("WHOLECELL_PLAN ")]
if _res.returncode != 0 or not _line:
    raise RuntimeError(f"wholecell.plan failed:\n{_res.stdout}\n{_res.stderr}")
PLAN = json.loads(_line[-1][len("WHOLECELL_PLAN "):])
BUILD_ID = config.get("build_id") or PLAN["build_id"]
BDIR = config.get("build_dir") or os.path.join(PLAN["results_root"], PLAN["cell_id"], BUILD_ID)
MAXW = int(config.get("max_workers") or PLAN["max_workers"])
FINAL = []          # targets of rule all (appended by the stage rule files)


def B(*a):
    return os.path.join(BDIR, *a)


def LOG(name):
    return B("logs", f"{name}.log")


def WC(stage, step=None, extra=""):
    """shell command of one stage step (the build environment is fixed here: PYTHONHASHSEED=0, repo on PYTHONPATH)"""
    # numba threads fixed in the environment before Python starts (wholecell/numba_env.py; 2026-11-20)
    nt = 1 if stage == "pack" and step in ("tiles", "pack") else 8
    s = (f"cd {REPO} && PYTHONHASHSEED=0 PYTHONPATH={REPO} NUMBA_NUM_THREADS={nt} NUMBA_THREADING_LAYER=workqueue "
         f"{PY} -m wholecell {stage} --cell {CELL} --build {BDIR}")
    if step:
        s += f" --step {step}"
    return s + (f" {extra}" if extra else "")


def FROZEN(name):
    return PLAN["frozen"][name]


def SAFE(name):
    return "".join(c if c.isalnum() else "_" for c in name)


STAMP = B("BUILD_STAMP.json")
FROZEN_OK = B("frozen_inputs.tsv")
GEOM = B("geometry", "package.json")

# plugin steps (config plugin_steps, wholecell/plugin_steps.py): their sentinels, and BEFORE(rule) = the sentinels of
# the plugin steps that list `rule` in `before` (added to that rule's inputs; empty for cells without plugin steps)
PLUG = PLAN.get("plugin_steps") or []
PLUG_DONE = lambda s: B("_smk", f"{s['stage']}_{s['id']}.done")


def BEFORE(rule_name):
    return [PLUG_DONE(s) for s in PLUG if rule_name in s["before"]]
