"""Adapters to the archived legacy Jurkat pipeline, used ONLY to run the ported validation on the archived v7 pack.

The legacy code tree is read-only: it is imported with bytecode writing disabled and the numba cache redirected to the
validation cache folder, so nothing is written next to it. The pipeline folder comes from config
(validate.compare_reference.legacy_pipeline_dir); no path is hard-coded here.

  dump_regions(pipeline_dir, out_json)  the v7 region table (jpack_test with JWC_REGIONS=v7) + localisation-key map
                                        (SOL/TMR) + label / membrane names, in the normalised form of regions.py
  LegacyGeometry(pipeline_dir, npz)     cell3d evaluator on a legacy geometry npz (jgeom_cell1_v7.npz)
"""
from __future__ import annotations

import importlib
import json
import os
import sys


def _import_legacy(pipeline_dir, modname, cache_dir=None, env=None):
    sys.dont_write_bytecode = True
    cd = os.environ.get("NUMBA_CACHE_DIR") or cache_dir
    if not cd:
        raise RuntimeError("legacy import needs a numba cache folder outside the legacy tree")
    os.environ["NUMBA_CACHE_DIR"] = cd
    if "numba" in sys.modules:             # numba already imported: its cache locator reads config.CACHE_DIR at compile time
        sys.modules["numba"].config.CACHE_DIR = cd
    for k, v in (env or {}).items():
        os.environ[k] = v
    if pipeline_dir not in sys.path:
        sys.path.insert(0, pipeline_dir)
    return importlib.import_module(modname)


def dump_regions(pipeline_dir, out_json, cache_dir=None):
    jp = _import_legacy(pipeline_dir, "jpack_test", cache_dir, {"JWC_REGIONS": "v7"})
    c3 = sys.modules["cell3d"]
    LN = list(c3.LABEL_NAMES); MN = list(c3.MEM_NAMES)
    regs = []
    for name, labs in jp.VOL_REG:
        regs.append(dict(name=name, kind="virtual" if name in jp.VIRTUAL else "volume", labels=[LN[l] for l in labs]))
    for name, m in jp.MEM_REG:
        regs.append(dict(name=name, kind="membrane", membrane=MN[m]))
    for name, m, labs, sg, d0, d1, rd in jp.SHELL_REG:
        regs.append(dict(name=name, kind="shell", membrane=MN[m], labels=[LN[l] for l in labs], side=float(sg), d0=float(d0),
                         d1=float(d1), add_radius=bool(rd & 1), pore_lat=float(rd >> 1)))
    assert [r["name"] for r in regs] == list(jp.RNAMES)
    rmap = {k: {"soluble": jp.SOL.get(k), "tm": jp.TMR.get(k)} for k in sorted(set(jp.SOL) | set(jp.TMR))}
    out = dict(source="legacy jpack_test (JWC_REGIONS=v7), dumped by wholecell.validate.legacy_v7", regions=regs,
               region_map=rmap, label_names=LN, mem_names=MN, mem_half=[float(x) for x in c3.MEM_HALF],
               lumen_mem={k: MN[v] for k, v in jp.LUMEN_MEM.items()})
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    json.dump(out, open(out_json, "w"), indent=1)
    return out


class LegacyGeometry:
    """labels(P) / mems(P) on a legacy geometry npz with the legacy cell3d evaluator"""

    def __init__(self, pipeline_dir, npz, cache_dir=None):
        import numpy as np
        self.c3 = _import_legacy(pipeline_dir, "cell3d", cache_dir)
        z = np.load(npz)
        self.cell = self.c3.Cell(z["G"], z["pores"])
        self.centre_nm = np.asarray(z["centre_nm"], float) if "centre_nm" in z.files else None
        self.pores = self.cell.pores
        self.label_names = list(self.c3.LABEL_NAMES)
        self.mem_names = list(self.c3.MEM_NAMES)
        self.mem_half = [float(x) for x in self.c3.MEM_HALF]
        self.source = f"legacy cell3d on {os.path.basename(npz)}"

    def labels(self, P):
        return self.cell.labels(P)

    def mems(self, P):
        return self.cell.mems(P)
