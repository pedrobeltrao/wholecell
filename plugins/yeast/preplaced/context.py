# origin: yeast scripts/s3_pack.py:main (stage-3 set-up before the pre-pass, lines "cell = ..." to "n_nuc -= n_rnuc") @ e3b7448
# Changed: wrapped as the registered preplaced_context "yeast.s3"; the inventory tables are read as scripts/_common.py
# load_inventory did; nothing else.
"""Yeast context of the serial pre-placed stage (wholecell/preplaced/serial.py): what the vendored stage 3 sets up
before its pre-pass, unchanged and in the same order:

  geometry (yc3d/geometry.json, cell3d.Cell), canonical tiles and tile measures, the s1 inventory (ingredients.csv,
  oligo_map.json), expected copies per (library type, region) (wc3d.expected_totals), the builder library entries
  (filament/tropomyosin/MT/ATP-dimer/eisosome/RAC/NVJ/amphiphysin models: loaded from the library, built only if
  missing), the library pack (wc3d.build_libpack over the same id set) and the nucleosome split arms / rDNA.

The numerical modules are the vendored engine's (plugins/yeast/engine/pipeline, never edited), pointed at the build
layout by assembly.engine_env (YC_* variables, before the first import). context_params.engine = {root, work, lib,
pdb, code, inputs, seed} ("build:<rel>" = inside the build); YC_WORK must hold yc3d/{geometry.json, geometry_fields.npz, tile_measure.parquet} (stage 2)
and results/inventory/{ingredients.csv, oligo_map.json, f_actin.json} (stage 1).
"""
from __future__ import annotations

import json
import os

from wholecell.registry import register


def _engine(E, bdir):
    from plugins.yeast.assembly import engine_env
    E = {k: (os.path.join(bdir, v.split(":", 1)[1]) if isinstance(v, str) and v.startswith("build:") else v) for k, v in E.items()}
    engine_env(E)
    import importlib
    return importlib.import_module


@register("preplaced_context", "yeast.s3")
def context(cfg, bdir, params):
    import numpy as np
    import pandas as pd
    from wholecell.preplaced.serial import Context, log
    imp = _engine(params["engine"], bdir)
    C = imp("config"); c3 = imp("cell3d"); t3 = imp("tiles3d"); p3 = imp("pack3d"); wc3d = imp("wc3d")
    cs = imp("cytoskeleton"); mt = imp("microtubules"); atp = imp("atp_dimers"); eis = imp("eisosomes"); rac = imp("rac")
    nvj = imp("nvj"); amph = imp("amphiphysin"); patches = imp("actin_patches"); rdna = imp("rdna")

    cell = c3.Cell.load(f"{C.W3D}/geometry.json")
    names = pd.read_parquet(f"{C.W3D}/tile_measure.parquet").index          # _common.canonical_tiles
    import re
    tiles = [tuple(int(x) for x in re.findall(r"[+-]\d+", n)) for n in names]
    measure = pd.read_parquet(f"{C.W3D}/tile_measure.parquet")
    assert list(measure.index) == [t3.tile_name(t) for t in tiles]
    inv_dir = f"{C.RESULTS}/inventory"                                       # _common.load_inventory
    ing = pd.read_csv(f"{inv_dir}/ingredients.csv")
    oligo_map = {k: tuple(v) for k, v in json.load(open(f"{inv_dir}/oligo_map.json")).items()}
    exp3 = wc3d.expected_totals(ing, oligo_map)
    log("expected per region:", exp3.groupby("region").expected.sum().round(0).to_dict())
    dropped = exp3[~exp3.region.isin(t3.REGION_NAMES)].groupby("region").expected.sum()
    log("not packed:", dropped.round(0).to_dict())
    cs.filament_model(); cs.tropomyosin_model(log=log); mt.model(); atp.model(log=log); eis.model(log=log); rac.build(log=log); nvj.model(log=log); amph.model(log=log)
    lib_ids = sorted(set(exp3[exp3.region.isin(t3.REGION_NAMES)].lib_id) | {"NPC", "nucleosome", cs.LIB_ID, cs.TPM_ID, atp.LIB_ID, *mt.LIB_IDS, *patches.LIB_IDS, *eis.LIB_IDS, rac.LIB_ID, nvj.LIB_ID, amph.LIB_ID}
                     | ({"ACT1"} if os.path.exists(f"{C.LIB}/ACT1.npz") else set()))
    lp = wc3d.build_libpack(lib_ids)
    n_nuc = int(round(ing.set_index("ing_id").loc["nucleosome", "copies_per_cell"]))
    n_rnuc = rdna.nucleosomes_in_rdna() if "preribo_90S" in set(lp["ids"]) else 0
    n_nuc -= n_rnuc

    BIG = np.full(3, 1e4)

    def emit(types, pos, q, lo=None, hi=None):
        return p3.emit_spheres(types, pos, q, lp["st_off"], lp["STC"], lp["STR"], lp["BR"],
                               np.full(3, -1e4) if lo is None else lo, BIG.copy() if hi is None else hi)

    from wholecell.preplaced import kern
    ctx = Context(cell=cell, lp=lp, exp3=exp3, log=log, mod=imp, emit=emit, pre_instances=wc3d.pre_instances,
                  ing=ing, oligo_map=oligo_map, tiles=tiles, measure=measure, inv_dir=inv_dir, lib_dir=C.LIB,
                  w3d=C.W3D, seed=C.seed, K=kern.yeast(p3, c3), lib_path=lambda i: f"{C.LIB}/{i}.npz")
    ctx.state.update(n_nuc=n_nuc, n_rnuc=n_rnuc)
    return ctx
