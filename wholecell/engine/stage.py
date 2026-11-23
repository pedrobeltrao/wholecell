"""External pipeline engines: python -m wholecell engine --cell C --build B --step S

A cell config may name `pipeline_engine: {name, steps, ...}`: a registered engine (registry kind `pipeline_engine`,
provided by a plugin) that runs a cell's own stage code inside a wholecell build directory, step by step, under the
shared workflow (workflow/rules/engine.smk builds the DAG from engine.steps(cfg)). Used for cells whose legacy
pipeline is ported as a whole and reproduced first (yeast X5: plugins/yeast, engine yeast_x5); stages move into the
shared core one at a time afterwards. Core stages (inventory copies, pack via pack.engine, export) are mixed into the
same DAG by the engine's step list.
"""
from .. import registry


def engine(cfg):
    pe = cfg.get("pipeline_engine") or {}
    if not pe.get("name"):
        raise SystemExit("cell config has no pipeline_engine")
    return registry.get("pipeline_engine", pe["name"])


def steps(cfg):
    return engine(cfg).steps(cfg)


def run(cfg, bdir, step=None, **kw):
    if not step:
        raise SystemExit("engine: --step is required")
    engine(cfg).run(cfg, bdir, step, **kw)
