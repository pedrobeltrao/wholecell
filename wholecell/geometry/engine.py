"""Active geometry evaluator shim: `from wholecell.geometry import engine as c3` (spec: docs/reference/cellgeom_v1.md §5).

Re-exports the namespace of the ONE core evaluator, cellgeom_eval_v1 (table-driven: rule tables in G; numba functions,
label / membrane codes, slot constants). Cell configs name it in `geometry.evaluator` (checked by cellgeom.load_cfg).
The legacy evaluators stay importable for the reproduction tests only: get_engine("sdf_sections_v1") or a dotted path.
"""
import importlib

DEFAULT_ENGINE = "cellgeom_eval_v1"


def get_engine(name=None):
    name = name or DEFAULT_ENGINE
    mod = name if "." in name else f"wholecell.geometry.engines.{name}"
    return importlib.import_module(mod)


_mod = get_engine()
ACTIVE_ENGINE = _mod.ENGINE_NAME
globals().update({k: v for k, v in vars(_mod).items() if not k.startswith("__")})
