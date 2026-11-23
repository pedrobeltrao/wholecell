"""Virtual volume regions: generic interface (regions of kind `virtual` in pack.regions).

A virtual region is a volume region whose voxels are selected by a rule instead of a label (e.g. a band around
microtubule axes). Providers are registered plugins of kind `virtual_region`, declared in config:

    pack:
      virtual_providers:
        - {id: cyto_v7, plugin: jurkat.cytoskeleton_v7, params: {...}, inputs: {name: path | 'layer:<package layer>'}}
      regions:
        - {name: mt_lattice, kind: virtual, labels: [cytosol], provider: cyto_v7, select: mt_lattice}

A provider object exposes
    load_context(params: dict, inputs: dict, geom) -> ctx          (called once, before the worker fork)
    masks(P, lab, ctx, selects: list[str], base_codes) -> {select: bool mask over P}
where `base_codes` are the engine label codes of the regions' base labels. The core calls masks only for points sets
that contain at least one base-label voxel (the legacy behaviour: a tile without cytosol has empty virtual regions);
the tile measure of a virtual region is its selected volume.
"""
from __future__ import annotations

import numpy as np

from wholecell import registry


class VirtualSet:
    def __init__(self, pcfg: dict, regions, geom):
        self.groups = []                               # (provider, ctx, [(region, select)], base_codes)
        if not regions.VIRTUAL:
            return
        provs = {p["id"]: p for p in (pcfg.get("virtual_providers") or [])}
        by = {}
        for name in regions.VIRTUAL:
            sp = regions.virtual_spec[name]
            by.setdefault(sp["provider"], []).append((name, sp.get("select", name)))
        vol = dict(regions.VOL_REG)
        for pid, items in by.items():
            p = provs[pid]
            obj = registry.get("virtual_region", p["plugin"])
            inputs = {k: geom.resolve_input(v) for k, v in (p.get("inputs") or {}).items()}
            ctx = obj.load_context(p.get("params") or {}, inputs, geom)
            base = np.array(sorted({c for n, _ in items for c in vol[n]}), np.int64)
            self.groups.append((obj, ctx, items, base))

    def __bool__(self):
        return bool(self.groups)

    def masks(self, P, lab):
        """{region name: mask} for every virtual region (None when the point set has no base-label voxel)"""
        out = {}
        for obj, ctx, items, base in self.groups:
            if np.isin(lab, base).any():
                mk = obj.masks(P, lab, ctx, [s for _, s in items], base)
                for name, s in items:
                    out[name] = mk[s]
            else:
                for name, _ in items:
                    out[name] = None
        return out
