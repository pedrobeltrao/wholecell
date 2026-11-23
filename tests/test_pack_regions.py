"""Tests of the region schema and the core tiling on yeast (unification plan A2), about a minute on 1 CPU:

    PYTHONHASHSEED=0 python tests/test_pack_regions.py        (or pytest tests/test_pack_regions.py)

  testcell   pack.regions is canonical (core packer); measure columns = volumes then mem_<membrane>
  schema     unknown kind / label / membrane / shell side give errors naming the region
  yeast      (needs the yeast port build s0 package) config/cells/yeast.yaml: 27 regions in the tiles3d id order, the
             kernel order / id maps / band arrays of the core packer (A7); the core tiling of 4 tiles == yc3d/geom/<tile>.npz of the port build s0
             (measure, labvol; when that build exists). Full check: tests/check_a2_yeast_tiles.py.
"""
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

YEAST_PKG = "yeast/analytic_x5_s0_f5f3a5c723b2"
YEAST_S0 = "yeast/20261116-aabe7b5-s0-x5repro"
# yeast tiles3d.REGION_NAMES @ e3b7448 (ids 0-26)
YEAST_NAMES = ["cytosol", "nucleoplasm", "nucleolus", "ER_lumen", "vacuole", "mito_matrix", "mito_IMS", "golgi_endosome",
               "peroxisome", "cell_wall", "PM_peripheral", "NE_peripheral", "ER_peripheral", "LD_peripheral", "mem_PM",
               "mem_ER", "mem_NE", "mem_VAC", "mem_MITO_OM", "mem_MITO_IM", "mem_GOLGI", "mem_PEROX",
               "chromatin_nucleosome", "NE_pore_NPC", "VAC_peripheral", "GOLGI_peripheral", "PEROX_peripheral"]


class _G:
    """minimal geometry stand-in for Regions (names -> codes, MEM_HALF)"""
    def __init__(self, cell):
        from wholecell.pack.geomload import _aliases
        from wholecell.geometry import engine
        self.engine = engine
        self.label_aliases = _aliases(cell.manifest, "labels")
        self.membrane_aliases = _aliases(cell.manifest, "membranes")
        self.params = engine.params_from_G(cell.G)


def _testcell():
    from wholecell.config import load
    from wholecell.geometry import cellgeom
    cfg = load("testcell")
    return cfg, cellgeom.load_cfg(cfg, verify="quick")


def test_testcell_canonical():
    from wholecell.pack.regions import Regions
    cfg, cell = _testcell()
    reg = Regions(cfg["pack"], _G(cell))
    assert not reg.not_packable and reg.RNAMES == [r["name"] for r in cfg["pack"]["regions"]] and reg.identity
    assert (reg.k_of == np.arange(len(reg.RNAMES))).all() and len(reg.band_arrays()[0]) == 0
    reg.require_core_packer()
    cols = [c for _, c, _, _ in reg.measures()]
    assert cols == [n for n, _ in reg.VOL_REG] + ["mem_" + n for n, _ in reg.MEM_REG], cols
    # membrane-protein test codes = the sdf_sections_v1 constants the v7 kernel had built in (outer wall, PM exempt, NE)
    c3 = cell.engine
    assert reg.membrane_codes(cell.G, c3).tolist() == [c3.M_OUTER, c3.M_PM, c3.M_NE]


def test_schema_errors():
    import copy
    from wholecell.pack.regions import Regions
    cfg, cell = _testcell()
    g = _G(cell)
    for bad, msg in [({"name": "x", "kind": "bogus"}, "unknown kind"), ({"name": "x", "kind": "volume", "labels": ["nope"]}, "label"),
                     ({"name": "x", "kind": "shell", "membrane": "PM", "labels": ["cytosol"], "d1": 5, "side": 2}, "side"),
                     ({"name": "x", "kind": "shell", "membrane": "XX", "labels": ["cytosol"], "d1": 5}, "membrane"),
                     ({"name": "x", "kind": "shell", "membrane": "PM", "labels": ["cytosol"]}, "missing d1")]:
        pc = copy.deepcopy(cfg["pack"]); pc["regions"] = pc["regions"] + [bad]
        try:
            Regions(pc, g)
        except ValueError as e:
            assert msg in str(e) and "x" in str(e), str(e)
        else:
            raise AssertionError(f"no error for {bad}")


def test_yeast_regions_and_tiles():
    from wholecell.config import load, load_site
    from wholecell.geometry import cellgeom
    from wholecell.pack.regions import Regions
    site = load_site()
    pkg = os.path.join(site["geom_root"], YEAST_PKG)
    if not os.path.isdir(pkg):
        print(f"SKIP: yeast package {pkg} missing"); return
    cfg = load("yeast")
    cell = cellgeom.load(pkg, verify="quick")
    reg = Regions(cfg["pack"], _G(cell))
    assert reg.RNAMES == YEAST_NAMES and [c for _, c, _, _ in reg.measures()] == YEAST_NAMES
    reg.require_core_packer()                            # A7: any order (kernel order + id maps)
    assert not reg.identity and len(reg.KNAMES) == 25 and reg.KNAMES[:17] == YEAST_NAMES[:14] + YEAST_NAMES[24:]
    assert [reg.RNAMES[c] for c in reg.c_of] == reg.KNAMES and all(reg.c_of[reg.k_of[i]] == i for i in range(27) if reg.k_of[i] >= 0)
    assert (reg.k_of[[22, 23]] == -1).all() and reg.NV == 17 and not reg.SHELL_REG
    bm, bw = reg.band_arrays()
    assert bm[reg.KNAMES.index("PM_peripheral")] == reg.mem_codes["PM"] and bw[reg.KNAMES.index("PM_peripheral")] == 22.0 + 4.0
    assert (reg.pack_class[:17] == 2).all() and (reg.pack_class[17:] == 0).all()
    assert reg.restriction("cell_wall")[1] == "lt" and reg.restriction("NE_peripheral")[1:] == (0.0, 0.0, 20.0)
    # membrane-protein tests (perox-fix): yeast codes, not the sdf_sections_v1 constants (7 = PEROX in yeast); the yeast
    # wall CELL_OUTER exempts no membrane (yeast pack3d.check_membrane), pores on NE
    mc = reg.membrane_codes(cell.G, _G(cell).engine).tolist()
    assert mc == [reg.mem_codes["CELL_OUTER"], -1, reg.mem_codes["NE"]] and mc[0] != reg.mem_codes["PEROX"], mc
    s0 = os.path.join(site["results_root"], YEAST_S0, "engine", "work", "yc3d", "geom")
    if not os.path.isdir(s0):
        print(f"SKIP: yeast port build s0 tiles {s0} missing"); return
    from wholecell.pack import stage, tiling
    from wholecell.pack.geomload import Geom
    g = Geom(cell, cell.engine, {}, None, pkg, _G(cell).label_aliases, _G(cell).membrane_aliases, cell.manifest)
    p = stage.P(cfg, tempfile.mkdtemp(), {"workers": 1})
    T, lo0, nt = tiling.select_enclosing_sphere(cell, cell.engine, p.TILE, 88.0)
    assert len(T) == 5912
    stage._S.update(p=p, reg=Regions(cfg["pack"], g), cell=cell, c3=cell.engine, T=T, lo0=lo0, vset=None)
    tc = dict(cfg["pack"]["tiling"]); tc["meshes"] = False
    tiling._S.update(stage._S, tcfg=tc, nlab=1 + max(cell.code_to_name))
    idx = [0, 2000, 2956, 4100]
    for r in tiling.measure_tiles(idx):
        nm = tiling.tile_name(T[r[0]] + np.round(lo0 / p.TILE).astype(np.int64))
        z = np.load(os.path.join(s0, nm + ".npz"))
        assert np.array_equal(np.array(r[1:28]), z["measure"]) and np.array_equal(r[-1], z["labvol"]), nm
    print(f"  yeast: 4 tiles identical to the port build s0 ({', '.join(str(i) for i in idx)})")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("pack region tests OK")
