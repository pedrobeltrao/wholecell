"""Cell configuration loading (spec: docs/reference/cell_config.md).

A cell is described by config/cells/<cell>.yaml. Nothing in the package knows a list of cells, compartments or species:
everything cell-specific comes from that file, the geometry package it points to, and the plugins it names.

Path placeholders resolved in every string value:
  {data_root}     site data root (config/site.yaml or env WHOLECELL_DATA_ROOT)
  {results_root}  where builds are written (config/site.yaml or env WHOLECELL_RESULTS_ROOT)
  {geom_root}     where immutable cellgeom_v1 packages live
  {yeast_root}    optional: the yeast model input tree (env WHOLECELL_YEAST_ROOT; yeast configs only)
  {afdb_root}     optional: folder with AlphaFold DB proteome tars (env WHOLECELL_AFDB_ROOT; Jurkat library build)
  {repo}          the repository root (the checkout running the build)
A placeholder whose site key is not set stays unresolved (the frozen-input check then reports the file MISSING).
No absolute path may appear in package code; site-specific absolute paths live only in config/site.yaml.
"""
from __future__ import annotations

import copy
import hashlib
import json
import os
import re

import yaml

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_PH = re.compile(r"\{(data_root|results_root|geom_root|repo|yeast_root|afdb_root)\}")


def load_site(path: str | None = None) -> dict:
    """site settings (absolute paths of this machine). Env vars override the file."""
    path = path or os.environ.get("WHOLECELL_SITE") or os.path.join(REPO, "config", "site.yaml")
    site = yaml.safe_load(open(path)) if os.path.exists(path) else {}
    for k in ("data_root", "results_root", "geom_root", "yeast_root", "afdb_root"):
        env = os.environ.get("WHOLECELL_" + k.upper())
        if env:
            site[k] = env
    site.setdefault("repo", REPO)
    return site


def _resolve(obj, site):
    if isinstance(obj, str):
        return _PH.sub(lambda m: str(site.get(m.group(1), m.group(0))), obj)
    if isinstance(obj, list):
        return [_resolve(x, site) for x in obj]
    if isinstance(obj, dict):
        return {k: _resolve(v, site) for k, v in obj.items()}
    return obj


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in over.items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def cell_config_path(cell: str) -> str:
    """config/cells/<cell>.yaml, or a path given directly"""
    return cell if cell.endswith((".yaml", ".yml")) else os.path.join(REPO, "config", "cells", f"{cell}.yaml")


# Config keys renamed with the stage "prepass" -> "preplaced" (A0c, 2026-11-17; unification plan s8). A config written
# before the rename still loads: the old keys are renamed on load, with one message naming them. Old key (dotted path)
# -> new key; `pack.regions[].kind: prepass` -> `preplaced` likewise. Seed stream names ("prepass.<module>") are NOT
# renamed: they define the seeds.
LEGACY_CONFIG_KEYS = {
    "prepass": "preplaced",
    "pack.prepass_index": "preplaced_index",
    "export.prepass_groups": "preplaced_groups",
    "export.prepass_group_regions": "preplaced_group_regions",
    "validate.densities.prepass_volume": "preplaced_volume",
    "validate.densities.prepass_group_volume": "preplaced_group_volume",
    "validate.loc_compare.prepass_ingredient": "preplaced_ingredient",
    "validate.loc_compare.prepass_region": "preplaced_region",
    "validate.loc_compare.prepass_region_default": "preplaced_region_default",
    "validate.chromatin.prepass_input": "preplaced_input",
    **{f"validate.compare_reference.files.prepass_{k}": f"preplaced_{k}"
       for k in ("types", "pos", "quat", "off", "libs", "groups")},
}


def legacy_keys(raw: dict, source: str = "config") -> dict:
    """rename pre-A0c keys (LEGACY_CONFIG_KEYS) in place; prints one line naming them"""
    done = []
    for old, new in LEGACY_CONFIG_KEYS.items():
        *parents, key = old.split(".")
        cur = raw
        for p in parents:
            cur = cur.get(p) if isinstance(cur, dict) else None
        if not isinstance(cur, dict) or key not in cur:
            continue
        if new in cur:
            raise ValueError(f"{source}: both '{old}' (old name) and '{new}' are set; keep only '{new}'")
        cur[new] = cur.pop(key)
        done.append(f"{old} -> {new}")
    for r in ((raw.get("pack") or {}).get("regions") or []) if isinstance(raw.get("pack"), dict) else []:
        if isinstance(r, dict) and r.get("kind") == "prepass":
            r["kind"] = "preplaced"
            done.append(f"pack.regions {r.get('name')}: kind prepass -> preplaced")
    if done:
        import sys
        print(f"{source}: stage 'prepass' is now 'preplaced'; old config keys renamed on load (update the file): "
              + "; ".join(done), file=sys.stderr)
    return raw


def load_raw(cell: str) -> dict:
    """raw YAML with `extends:` resolved (a config may extend another config file, relative to its own folder);
    pre-A0c key names are renamed (legacy_keys)"""
    path = cell_config_path(cell)
    raw = legacy_keys(yaml.safe_load(open(path)), path)
    if raw.get("extends"):
        parent = load_raw(os.path.join(os.path.dirname(path), raw.pop("extends")))
        raw = _merge(parent, raw)
    return raw


class CellConfig(dict):
    """resolved cell config (a dict) with helpers"""

    @property
    def cell_id(self) -> str:
        return self["cell_id"]

    def get_path(self, dotted: str, default=None):
        cur = self
        for part in dotted.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return default
            cur = cur[part]
        return cur

    def frozen(self, name: str) -> str:
        """path of a declared frozen input (checksum verified by the check_frozen stage, see build.verify_frozen)"""
        fi = self.get("frozen_inputs", {})
        if name not in fi:
            raise KeyError(f"frozen input '{name}' is not declared in config frozen_inputs")
        return fi[name]["path"]

    def config_hash(self) -> str:
        """sha256 of the canonical JSON of the UNRESOLVED config (placeholders kept), so the hash does not depend on the
        machine's site paths; volatile keys (_site, _source, _raw) excluded"""
        c = self.get("_raw") or {k: v for k, v in self.items() if not k.startswith("_")}
        return hashlib.sha256(json.dumps(c, sort_keys=True, default=str).encode()).hexdigest()


def load(cell: str, site_path: str | None = None, overrides: dict | None = None) -> CellConfig:
    site = load_site(site_path)
    raw = load_raw(cell)
    if overrides:
        raw = _merge(raw, overrides)
    cfg = CellConfig(_resolve(raw, site))
    cfg["_site"] = site
    cfg["_raw"] = raw
    cfg["_source"] = os.path.abspath(cell_config_path(cell))
    for req in ("cell_id", "schema_version"):
        if req not in cfg:
            raise ValueError(f"cell config {cell}: missing required key '{req}'")
    return cfg
