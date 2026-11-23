"""Plugin registry: everything cell-specific is registered by name and selected from the cell config.

Kinds (open set; a new kind needs no change here):
  geometry_engine, adapter.proteome, adapter.rna, adapter.localisation, adapter.assembly, adapter.surface,
  inventory_rule, preplaced (+ preplaced_context), virtual_region, validator, exporter

    from wholecell.registry import register
    @register("inventory_rule", "jurkat.v4")
    def apply(ctx): ...

Plugins are Python packages under plugins/ (or anywhere importable) listed in the cell config `plugins:`; importing
them runs their @register decorators. Core components register themselves the same way when their module is imported.
"""
from __future__ import annotations

import importlib
import os
import sys

_REG: dict[tuple[str, str], object] = {}


def register(kind: str, name: str):
    def deco(obj):
        key = (kind, name)
        if key in _REG and _REG[key] is not obj:
            raise ValueError(f"duplicate registration {kind}:{name}")
        _REG[key] = obj
        return obj
    return deco


def get(kind: str, name: str):
    if (kind, name) not in _REG:
        known = sorted(n for k, n in _REG if k == kind)
        raise KeyError(f"no {kind} named '{name}' registered (known: {known}); is its plugin listed in the cell config?")
    return _REG[(kind, name)]


def names(kind: str) -> list[str]:
    return sorted(n for k, n in _REG if k == kind)


def load_plugins(cfg) -> list[str]:
    """import the plugin packages named in cfg['plugins'] (module paths; 'jurkat' -> plugins.jurkat)"""
    from .config import REPO
    if REPO not in sys.path:
        sys.path.insert(0, REPO)
    done = []
    for p in cfg.get("plugins", []) or []:
        mod = p if "." in p else f"plugins.{p}"
        importlib.import_module(mod)
        done.append(mod)
    return done
