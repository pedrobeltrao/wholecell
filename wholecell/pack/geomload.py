"""Geometry access for the pack and preplaced stages.

Resolution order of the geometry used by a build:
  1. <bdir>/geometry/package.json written by the geometry stage ({"package_dir": ...} or {"path": ...})
  2. config `geometry.package` (a cellgeom_v1 package directory)
  3. config `pack.dev_legacy_geometry` (development only: a legacy cell3d npz converted with the engine helpers
     split_legacy_G + build_G; used by the pack-port tests before the package existed)
Returns a Geom object with .cell (engine Cell), .engine (module), .params (engine params from G), .summary (legacy
jgeom-summary shape), .centre_nm (model-frame origin in the source volume, nm), .source (path), and optional name maps
(.label_aliases, .membrane_aliases) from the package manifest.
"""
from __future__ import annotations

import json
import os

import numpy as np


class Geom:
    def __init__(self, cell, engine, summary, centre_nm, source, label_aliases=None, membrane_aliases=None, manifest=None):
        self.cell = cell
        self.engine = engine
        self.summary = summary or {}
        self.centre_nm = None if centre_nm is None else np.asarray(centre_nm, np.float64)
        self.source = source
        self.label_aliases = label_aliases or {}
        self.membrane_aliases = membrane_aliases or {}
        self.manifest = manifest or {}

    @property
    def params(self):
        return self.engine.params_from_G(self.cell.G)

    def layer(self, name: str, i: int = 0) -> str:
        """path of a package layer file (cellgeom_v1 layers: ne_mesh, mt_traces, centrosome, ...)"""
        if not hasattr(self.cell, "layer_path"):
            raise KeyError(f"geometry {self.source} has no layers (legacy npz); give the input as a path")
        return self.cell.layer_path(name, i)

    def resolve_input(self, v):
        """'layer:<name>' -> package layer file; anything else is returned unchanged (a path)"""
        return self.layer(v.split(":", 1)[1]) if isinstance(v, str) and v.startswith("layer:") else v

    def flag(self, slot_name: str) -> bool:
        """engine header slot (e.g. 'G_FACE') > 0"""
        return bool(self.cell.G[getattr(self.engine, slot_name)] > 0)


def package_dir(cfg, bdir=None):
    if bdir:
        pj = os.path.join(bdir, "geometry", "package.json")
        if os.path.exists(pj):
            j = json.load(open(pj))
            return j.get("package_dir") or j.get("path")
    return (cfg.get("geometry") or {}).get("package")


def _aliases(manifest, key):
    """optional name -> engine code maps from a package manifest (several tolerated shapes)"""
    m = manifest.get(key) or {}
    out = {}
    if isinstance(m, dict):
        for k, v in m.items():
            if isinstance(v, int):
                out[k] = v
            elif isinstance(v, dict) and "engine_code" in v:
                out[k] = int(v["engine_code"])
    elif isinstance(m, list):
        for e in m:
            if isinstance(e, dict) and "name" in e and ("engine_code" in e or "code" in e):
                c = e.get("engine_code", e.get("code"))
                if c is not None:                 # vocabulary entries the package does not use have no code
                    out[e["name"]] = int(c)
    return out


def load(cfg, bdir=None) -> Geom:
    from wholecell.geometry import engine as c3
    pdir = package_dir(cfg, bdir)
    if pdir and os.path.isdir(pdir):
        from wholecell.geometry import cellgeom
        cell = cellgeom.load_cfg(cfg, pdir)
        man = getattr(cell, "manifest", {}) or {}
        summ = getattr(cell, "summary", None) or man.get("summary") or {}
        centre = getattr(cell, "centre_nm", None)
        if centre is None:
            fr = man.get("frame") or {}
            centre = fr.get("centre_label_nm", fr.get("centre_nm", summ.get("centre_nm")))
        return Geom(cell, c3, summ, centre, pdir, _aliases(man, "labels"), _aliases(man, "membranes"), man)
    leg = (cfg.get("pack") or {}).get("dev_legacy_geometry")
    if leg:
        z = np.load(leg)
        G = c3.build_G(**c3.split_legacy_G(z["G"]))
        cell = c3.Cell(G, z["pores"])
        return Geom(cell, c3, json.loads(str(z["summary"])), z["centre_nm"] if "centre_nm" in z.files else None, leg)
    raise FileNotFoundError("no geometry: <bdir>/geometry/package.json, config geometry.package or pack.dev_legacy_geometry")
