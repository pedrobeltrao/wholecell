"""Geometry access for validation: labels(P), mems(P), label/membrane names, centre, pores.

Sources (first that applies):
  geom=<file.npz>        legacy geometry npz (G, pores, centre_nm), evaluated with the legacy cell3d of
                         validate.compare_reference.legacy_pipeline_dir (only to validate the archived v7 pack)
  geom=<package dir>     a cellgeom_v1 package
  default                the build's geometry package: <build>/geometry/package (link written by the geometry stage),
                         else config geometry.package / geometry.package_dir
Label and membrane names are those of the engine (LABEL_NAMES / MEM_NAMES); region tables refer to these names or to
package vocabulary names mapped to engine codes by the package manifest.
"""
from __future__ import annotations

import json
import os

import numpy as np


class PackageGeometry:
    def __init__(self, pkg_dir, cfg=None):
        from ..geometry import cellgeom
        self.cell = cellgeom.load_cfg(cfg or {}, pkg_dir)
        from ..geometry import engine as c3
        self.c3 = c3
        man = getattr(self.cell, "manifest", {}) or {}
        self.label_names = list(getattr(c3, "LABEL_NAMES", []))
        self.mem_names = list(getattr(c3, "MEM_NAMES", ["PM", "NE", "ER", "MITO_OM", "MITO_IM", "GOLGI", "ENDO", "CELL_OUTER"]))
        # A8 phase 2: the package's own engine vocabulary (manifest labels[].engine_code/engine_name, membranes[].engine_code)
        # names the codes its evaluator returns. For sdf_sections_v1 packages (Jurkat, testcell) these are the engine's
        # LABEL_NAMES / MEM_NAMES (unchanged); a package of another layout (yeast_cell3d_v1: 15 labels, 10 membranes in
        # another order) was read with the human names before, which mislabelled every yeast centre.
        labs = [L for L in man.get("labels", []) if L.get("engine_code") is not None]
        if labs:
            n = 1 + max(int(L["engine_code"]) for L in labs)
            names = [""] * n
            for L in labs:
                names[int(L["engine_code"])] = L.get("engine_name") or L["name"]
            if names[:len(self.label_names)] != self.label_names[:n] or n != len(self.label_names):
                self.label_names = names
        mems = man.get("membranes") or []
        own_mems = False
        if mems and all(m.get("engine_code") is not None for m in mems):
            mn = [m["name"] for m in sorted(mems, key=lambda m: int(m["engine_code"]))]
            if mn != self.mem_names:
                self.mem_names = mn; own_mems = True
        G = self.cell.G
        if own_mems:
            self.mem_half = [float(m.get("half_thickness_nm") or 2.0) for m in sorted(mems, key=lambda m: int(m["engine_code"]))]
        elif hasattr(c3, "mem_half"):
            self.mem_half = [float(c3.mem_half(G, m)) for m in range(len(self.mem_names))]
        else:
            self.mem_half = [float(x) for x in getattr(c3, "MEM_HALF", np.full(len(self.mem_names), 2.0))]
        self.vocab = {}
        c2n = getattr(self.cell, "code_to_name", None) or {}
        for code, name in c2n.items():                   # package vocabulary names -> engine codes
            self.vocab[name] = int(code)
        for code, lid in (getattr(self.cell, "code_to_id", None) or {}).items():
            self.vocab[str(lid)] = int(code)
        fr = man.get("frame", {}) or {}
        c = fr.get("centre_label_nm") or fr.get("centre_nm") or man.get("centre_nm")
        self.centre_nm = np.asarray(c, float) if c is not None else None
        self.pores = np.asarray(getattr(self.cell, "pores", np.zeros((0, 6))))
        self.source = f"cellgeom_v1 package {pkg_dir}"

    def labels(self, P):
        return np.asarray(self.cell.labels(np.asarray(P, np.float64)))

    def mems(self, P):
        return np.asarray(self.cell.mems(np.asarray(P, np.float64)))


def open_geometry(cfg, inp, bdir):
    g = inp.get("geom")
    cache = os.path.join(bdir, "validate", "cache", "numba")
    if g and g.endswith(".npz"):
        from .legacy_v7 import LegacyGeometry
        pdir = (cfg.get("validate", {}).get("compare_reference", {}) or {}).get("legacy_pipeline_dir")
        if not pdir:
            raise RuntimeError("legacy geometry npz needs validate.compare_reference.legacy_pipeline_dir")
        return LegacyGeometry(pdir, g, cache)
    pj = os.path.join(bdir, "geometry", "package.json")          # yeast port: package dir recorded, no link
    pj = json.load(open(pj)).get("package_dir") if os.path.exists(pj) else None
    for cand in (g, os.path.join(bdir, "geometry", "package"), (cfg.get("geometry", {}) or {}).get("package"),
                 (cfg.get("geometry", {}) or {}).get("package_dir"), pj):
        if cand and os.path.exists(cand):
            return PackageGeometry(cand, cfg)
    raise RuntimeError("no geometry: pass geom=<package dir | legacy npz> or build the geometry stage first")


def label_code_map(geo):
    """name -> engine label code (engine names, plus package vocabulary names if the manifest maps them)"""
    m = {n: i for i, n in enumerate(geo.label_names)}
    voc = getattr(geo, "vocab", {}) or {}
    if isinstance(voc, dict):
        for k, v in voc.items():
            if isinstance(v, dict) and "engine_code" in v:
                m.setdefault(k, int(v["engine_code"]))
            elif isinstance(v, int):
                m.setdefault(k, v)
    return m
