"""Workflow plan: python -m wholecell.plan --cell C [--build-tag T] -> JSON on stdout.

The Snakemake environment is orchestration only (no numpy / pandas / numba), so the Snakefile asks the build
environment for everything it needs to lay out the DAG: resolved frozen-input paths, the inventory step plan
(wholecell.inventory.stage.steps), pre-placed modules / serial builders, plugin steps, validation / export steps, worker
count and the default build id. `pipeline_engine` (reproduction configs only, e.g. config/cells/yeast_x5_port.yaml):
the DAG is the engine's step list instead.
"""
from __future__ import annotations

import argparse
import json
import os
import sys


def plan(cell: str, build_tag: str | None = None) -> dict:
    from . import build, registry
    from .config import load
    cfg = load(cell)
    registry.load_plugins(cfg)
    out = dict(cell=cell, cell_id=cfg["cell_id"], python=cfg["_site"].get("python") or sys.executable,
               build_id=build.make_build_id(cfg, build_tag), results_root=cfg["_site"].get("results_root"),
               max_workers=int((cfg.get("pack") or {}).get("max_workers") or cfg["_site"].get("max_workers") or 8),
               frozen={k: v["path"] for k, v in (cfg.get("frozen_inputs") or {}).items()},
               geometry_package=(cfg.get("geometry") or {}).get("package", ""))
    pe = cfg.get("pipeline_engine") or {}
    out["pipeline_engine"] = pe.get("name")
    from .pack.switches import effective
    pce = effective(cfg.get("pack") or {})                # pack.preset merged (A7)
    out["pack_engine"] = pce.get("engine") or "v7"
    og = pce.get("organisation")
    out["pack_organisation"] = bool(og if isinstance(og, bool) else (og or {}).get("enabled", False))
    if pe.get("name"):
        # whole-pipeline engine (plugin): the DAG comes from its step list (workflow/rules/engine.smk)
        from .engine import stage as eng
        out["engine_steps"] = eng.steps(cfg)
        return out
    from .inventory import stage as inv
    from .inventory import copies, generic, secms  # noqa: F401  (core adapters/builders/rules register themselves, as in inventory run())
    out["inventory"] = inv.steps(cfg)
    lc = cfg.get("library") or {}
    out["library_nproc"] = int(lc.get("nproc", 8)) if lc.get("mode") == "build" else 1   # build mode: processes of the builders
    out["preplaced_modules"] = (cfg.get("preplaced") or {}).get("modules") or []
    ser = (cfg.get("preplaced") or {}).get("serial")
    out["preplaced_serial"] = {"out": ser.get("out", "preplaced/serial")} if ser else None
    from .plugin_steps import specs
    out["plugin_steps"] = specs(cfg)                    # cell-specific steps of the stage list (workflow/rules/plugin_steps.smk)
    try:
        from .validate import stage as val
        out["validate_steps"] = list(val.STEPS_ALL)
    except ImportError:
        out["validate_steps"] = []
    out["compare_reference"] = bool((cfg.get("validate") or {}).get("compare_reference"))
    ex = cfg.get("export") or {}
    out["export_steps"] = (ex.get("steps") or ["all"]) if (ex.get("enabled", True) and
                                                           os.path.exists(os.path.join(os.path.dirname(__file__), "export", "stage.py"))) else []
    out["geometry_stage"] = os.path.exists(os.path.join(os.path.dirname(__file__), "geometry", "stage.py"))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", required=True)
    ap.add_argument("--build-tag", default=None)
    a = ap.parse_args(argv)
    os.environ.setdefault("PYTHONHASHSEED", "0")
    p = plan(a.cell, a.build_tag)
    sys.stdout.write("WHOLECELL_PLAN " + json.dumps(p) + "\n")


if __name__ == "__main__":
    main()
