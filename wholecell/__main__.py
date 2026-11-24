"""Stage dispatcher: python -m wholecell <stage> --cell <cell|config.yaml> --build <build_dir> [--step S] [key=value ...]

--cell may be left out for an existing build (the build's config.resolved.json is used), e.g. the standalone web-tile
re-export: python -m wholecell export --build <build_dir> --step web --out <dir>

Each stage package exposes wholecell/<stage>/stage.py:run(cfg, bdir, step=None, **kw). Plugins named in the cell config
are imported first, so registered engines/adapters/rules are available. PYTHONHASHSEED must be 0.
Special stages: `stamp` (create the build dir, BUILD_STAMP.json), `frozen` (verify frozen inputs).
`--step plugin:<id>` runs a cell-specific step of that stage (config plugin_steps, wholecell/plugin_steps.py).
Deprecated stage names (STAGE_ALIASES) still run, with one line naming the new one.
"""
import argparse
import importlib
import os
import sys

from . import build, registry
from .config import load

STAGE_ALIASES = {"prepass": "preplaced"}       # A0c (2026-11-17): stage names follow the paper's Methods


def main(argv=None):
    ap = argparse.ArgumentParser(prog="wholecell")
    ap.add_argument("stage")
    ap.add_argument("--cell", default=None, help="cell id or config.yaml (default: the build's config.resolved.json)")
    ap.add_argument("--build", required=True, help="build directory (results_root/<cell>/<build_id>)")
    ap.add_argument("--step", default=None)
    ap.add_argument("--out", default=None, help="output folder for steps that support it (validate --step suite, export --step web)")
    ap.add_argument("kv", nargs="*", help="extra key=value arguments passed to the stage")
    a = ap.parse_args(argv)
    if a.stage in STAGE_ALIASES:
        print(f"stage '{a.stage}' is now '{STAGE_ALIASES[a.stage]}'", file=sys.stderr)
        a.stage = STAGE_ALIASES[a.stage]
    from . import numba_env                              # before any plugin imports numba (2026-11-20)
    numba_env.ensure(numba_env.wanted(a.stage, a.step), force=a.stage == "pack" and a.step in ("tiles", "pack"))
    build.check_env()
    if a.cell:
        cfg = load(a.cell)
    else:                                                # an existing build: its own resolved config
        rp = os.path.join(a.build, "config.resolved.json")
        if not os.path.exists(rp):
            ap.error(f"--cell not given and {rp} does not exist")
        from .config import CellConfig
        import json
        cfg = CellConfig(json.load(open(rp)))
    registry.load_plugins(cfg)
    kw = dict(x.split("=", 1) for x in a.kv)
    if a.out:
        kw["out"] = a.out
    if a.stage == "stamp":
        build.write_stamp(cfg, a.build)
        return
    if a.stage == "frozen":
        build.verify_frozen(cfg, os.path.join(a.build, "frozen_inputs.tsv"))
        return
    if a.step and a.step.startswith("plugin:"):          # a cell-specific step of this stage (config plugin_steps)
        from .plugin_steps import run_step
        run_step(cfg, a.build, a.stage, a.step.split(":", 1)[1], **kw)
        return
    mod = importlib.import_module(f"wholecell.{a.stage}.stage")
    mod.run(cfg, a.build, step=a.step, **kw)


if __name__ == "__main__":
    sys.exit(main())
