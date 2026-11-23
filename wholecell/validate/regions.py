"""Region table for validation (normalised view of config `pack.regions`, see docs/ARCHITECTURE.md section 6).

Normalised region record (list order = job_region index):
  name, kind ('volume' | 'membrane' | 'shell' | 'virtual'),
  labels   [label names]   volume / virtual / shell: allowed labels of the instance centre
  membrane name            membrane / shell: the membrane (engine membrane name)
  side     +1 / -1         shell: side of the membrane SDF (+ = cytoplasmic / nucleoplasmic / matrix side)
  d0, d1   nm              shell: band of signed distance (side * sdf) of the centre; d1 a number (+ the membrane half
                           thickness if d1_mem_half is true, as in pack.regions), 'anchor' (2 + mem half) or 'mem_half+X'
  add_radius bool          shell: the library bounding radius is added to d1
  pore_lat nm              shell: centre within this distance of a pore axis (0 = no limit; pack.regions: pore_limit_nm)
Sources, in order: regions_json= (a JSON list of such records, e.g. the dump of the legacy v7 table made by
legacy_v7.dump_regions), config pack.regions (list of dicts or {name: dict}), else the libs.json `rnames` with unknown kinds
(only steps that need no kinds can run then).
"""
from __future__ import annotations

import json
import re


def _norm_one(name, r):
    r = dict(r or {})
    r.setdefault("name", name)
    kind = r.get("kind") or r.get("type") or "volume"
    r["kind"] = kind
    if isinstance(r.get("labels"), str):
        r["labels"] = [r["labels"]]
    r.setdefault("labels", [])
    r.setdefault("side", 1.0)
    r.setdefault("d0", 0.0)
    r["add_radius"] = bool(r.get("add_radius", r.get("rad", False)))
    r["pore_lat"] = float(r.get("pore_lat", r.get("pore_limit_nm", r.get("pore_limit", 0.0))) or 0.0)
    return r


def load_regions(cfg, inp=None):
    p = inp.get("regions_json") if inp is not None else None
    if p:
        R = json.load(open(p))
        R = R["regions"] if isinstance(R, dict) and "regions" in R else R
        return [_norm_one(r.get("name"), r) for r in R]
    reg = (cfg.get("pack", {}) or {}).get("regions")
    if reg:
        if isinstance(reg, dict):
            return [_norm_one(k, v) for k, v in reg.items()]
        return [_norm_one(r.get("name"), r) for r in reg]
    if inp is not None and inp.has("libs"):
        L = json.load(open(inp["libs"]))
        if "rnames" in L:
            return [_norm_one(n, {"kind": "unknown"}) for n in L["rnames"]]
    raise RuntimeError("no region table: config pack.regions missing and no regions_json= given")


def region_map(cfg, inp=None):
    """localisation key -> {soluble: region, tm: region} (config pack.region_map, or the regions_json dump)"""
    p = inp.get("regions_json") if inp is not None else None
    if p:
        R = json.load(open(p))
        if isinstance(R, dict) and R.get("region_map"):
            return R["region_map"]
    return (cfg.get("pack", {}) or {}).get("region_map") or {}


def resolve_d1(r, mem_half):
    """numeric d1 of a shell (mem_half = that membrane's half thickness)"""
    d1 = r.get("d1", 0.0)
    if isinstance(d1, (int, float)):
        return float(d1) + (mem_half if r.get("d1_mem_half") else 0.0)   # pack.regions form (wholecell.pack.regions)
    s = str(d1).strip()
    if s == "anchor":
        return 2.0 + mem_half
    m = re.fullmatch(r"mem_half\s*\+\s*([0-9.]+)", s)
    if m:
        return mem_half + float(m.group(1))
    return float(s)


def names(regions):
    return [r["name"] for r in regions]


def by_kind(regions, *kinds):
    return [i for i, r in enumerate(regions) if r["kind"] in kinds]
