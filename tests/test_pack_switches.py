"""Tests of the packer switches (wholecell/pack/switches.py; unification plan A7 phase 1), a few minutes on 1 CPU:

    PYTHONHASHSEED=0 python tests/test_pack_switches.py        (or pytest tests/test_pack_switches.py)

  switches       defaults are the v7 values; config errors name the key and suggest the nearest valid one
  presets (B1)   COMMON_DEFAULTS applied to unpinned cells; pins legacy_v7 / legacy_x5 / none; preset lists; benchmark
                 settings; the reproduction presets keep their effective pack: section (md5)
  algorithms     size-class rounds, pilot tile choice, binomial thinning, shell cap, budgets, squeeze job choice,
                 exclusion tables (plain Python; small hand-made inputs)
  kernel         pack3d._in_exc numba == pure-Python reference; pack_tile on a synthetic 120 nm cytosol tile of the
                 testcell: empty switch arrays == explicit zero arrays (bit-identical), an exclusion zone covering the
                 tile keeps its species out, a per-region budget of 1 try places no more than the default
  squeeze        relax / accept numba == their pure-Python source (py_func) on a few molecules; squeeze_tile on an
                 over-asked tile places at least as many copies as RSA and its result has no overlap (accept keeps all)
"""
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell.pack import algorithms as A  # noqa: E402
from wholecell.pack import switches  # noqa: E402

LIBDIR = os.path.join(REPO, "tests", "data", "testcell", "lib")


# ------------------------------------------------------------------ switches
def test_switch_defaults():
    sw = switches.load({}, ["cytosol"])
    assert sw.legacy
    assert sw.rounds["n"] == 1 and not sw.thinning["enabled"] and not sw.shell_cap["enabled"]
    assert sw.budgets["per_region"] == {} and sw.squeeze["scope"] == "none" and sw.exclusion == []
    on = switches.load({"rounds": {"n": 6}, "thinning": {"enabled": True}, "shell_cap": True,
                        "squeeze": {"scope": "underfilled", "threshold": 0.97}}, ["cytosol"])
    assert not on.legacy and on.rounds["alternate"] is True and on.thinning["pilot_tiles"] == 16
    for bad, words in (({"rounds": {"nn": 2}}, "did you mean 'n'"), ({"squeeze": {"scope": "underfiled"}}, "did you mean 'underfilled'"),
                       ({"budgets": {"per_region": {"cytosl": {"max_tries": 5}}}}, "did you mean 'cytosol'"),
                       ({"exclusion": [{"libs": ["A"], "r_nm": 5}]}, "around")):
        try:
            switches.load(bad, ["cytosol"])
        except ValueError as e:
            assert words in str(e), (bad, str(e))
        else:
            raise AssertionError(f"no error for {bad}")


# ------------------------------------------------------------------ plain-Python algorithms
def test_rounds():
    so = np.array([0, 1, 2, 3, 4]); STR = np.array([4.0, 3.0, 2.0, 1.0]); BR = STR.copy()
    pack_class = np.array([2, 0])                         # region 0 volume, region 1 membrane
    lib = np.array([0, 1, 2, 3, 0]); region = np.array([0, 0, 0, 0, 1]); n = np.array([1, 1, 1, 64, 5])
    assert len(A.round_edges(lib, region, n, so, STR, BR, pack_class, 1)) == 0
    e = A.round_edges(lib, region, n, so, STR, BR, pack_class, 2)
    # sphere volumes r^3 x n: 64, 27, 8, 64 -> half of 163 is reached at lib 1 (r 3)
    assert e.tolist() == [3.0], e
    r = A.job_rounds(BR[[0, 1, 2, 3, 0]], np.array([False, False, False, False, True]), e)
    assert r.tolist() == [0, 0, 1, 1, 0], r
    assert A.phase_order(range(8), 1, True) == list(range(7, -1, -1)) and A.phase_order(range(8), 1, False) == list(range(8))


def test_thinning_parts():
    rng = np.random.default_rng(1)
    Tij = np.array([[i, 0, 0] for i in range(20)])
    ch = A.pilot_tiles([0, 1], {0: np.ones(20), 1: np.r_[np.zeros(10), np.ones(10)]}, Tij, 4, rng)
    allc = [t for v in ch.values() for t in v]
    assert len(allc) == len(set(allc)) and all(10 <= t for t in ch[1])
    assert all(abs(a - b) > 1 for i, a in enumerate(allc) for b in allc[:i])
    assert A.thinning_factor(0.99, 0.97, 0.98) == 1.0 and abs(A.thinning_factor(0.8, 0.97, 0.98) - 0.776) < 1e-12
    t = np.arange(1000); l = np.zeros(1000, int); r = np.r_[np.zeros(500, int), np.ones(500, int)]; n = np.full(1000, 100)
    t2, l2, r2, n2 = A.thin_counts(t, l, r, n, {}, np.random.default_rng(0))
    assert np.array_equal(n2, n)
    t2, l2, r2, n2 = A.thin_counts(t, l, r, n, {1: 0.5}, np.random.default_rng(0))
    assert n2[r2 == 0].sum() == 50000 and abs(n2[r2 == 1].sum() / 50000 - 0.5) < 0.01
    vol_lib = np.array([10.0]); crowd, phi = A.crowded_regions(l, r, n, vol_lib, np.array([1e6, 1e5]), [0, 1], 0, 1.1)
    assert list(crowd) == [1] and abs(phi[1] / phi[0] - 10) < 1e-9


def test_shell_cap():
    tot = {(0, 0): 1000.0, (0, 2): 1000.0, (1, 2): 500.0}
    vol_lib = np.array([100.0, 100.0]); br = np.array([3.0, 3.0])
    shells = [dict(ri=2, area_nm2=1e4, d0=0.0, d1=10.0, add_radius=False, pore_limit=0, labels=(1,))]
    out, rep = A.shell_cap(tot, vol_lib, br, {0: (1,)}, {0: 1e6}, shells, log=lambda *a: None)
    assert abs(sum(out.values()) - sum(tot.values())) < 1e-6           # copies conserved
    f, moved, target = rep[2]
    # capped at the bulk region's requested fraction before the spill (as in the yeast packer): 1000 x 100 / 1e6
    assert target == 0 and abs((out[(0, 2)] + out[(1, 2)]) * 100 / 1e5 - 1000 * 100 / 1e6) < 1e-12
    assert abs(out[(0, 2)] / out[(1, 2)] - 2.0) < 1e-12                 # every species of the shell by the same factor
    out2, rep2 = A.shell_cap(tot, vol_lib, br, {0: (1,)}, {0: 1e6}, [dict(shells[0], pore_limit=200)], log=lambda *a: None)
    assert out2 == tot and not rep2


def test_budgets_squeeze_jobs_exclusion():
    tr, gu = A.budget_arrays({}, ["a", "b"])
    assert len(tr) == 0 and len(gu) == 0
    tr, gu = A.budget_arrays({"b": {"max_tries": 7}}, ["a", "b"])
    assert tr.tolist() == [0, 7] and gu.tolist() == [0, 0]
    jl = np.array([0, 1, 0]); jr = np.array([0, 0, 1]); jn = np.array([100, 100, 10]); pl = np.array([100, 98, 5])
    mass = np.array([1.0, 1.0]); lv = np.array([True, True])
    assert A.squeeze_jobs("underfilled", 0.97, jl, jr, jn, pl, mass, lv).tolist() == [False, False, True]
    assert A.squeeze_jobs("all", 0.97, jl, jr, jn, pl, mass, lv).tolist() == [True, True, True]
    assert not A.squeeze_jobs("none", 0.97, jl, jr, jn, pl, mass, lv).any()
    et, around, rad = A.exclusion_tables([{"libs": ["R80"], "around": {"lib_prefixes": ["ACT"]}, "r_nm": 50},
                                          {"lib_prefixes": ["R"], "around": {"libs": ["X"]}, "r_nm": 10}],
                                         ["R80", "P1"], ["ACTB", "X"])
    assert et.tolist() == [3, 0] and around.tolist() == [[True, False], [False, True]]
    E = A.exclusion_spheres(np.array([0, 1, 0]), np.array([[0, 0, 0], [5, 5, 5], [900, 0, 0]], float), around, rad,
                            np.full(3, -100.0), np.full(3, 100.0))
    assert E.shape == (2, 5) and sorted(E[:, 4].tolist()) == [0.0, 1.0]


def test_preset_and_providers():
    from wholecell.pack import switches
    pc = {"max_tries": 200, "requests": {"method": "legacy_x5", "thinning": "t.json", "lib_order": "l.npy"}, "preset": "x",
          "presets": {"x": {"max_tries": 2000, "requests": {"thinning": None}, "rounds": {"n": 6}}}}
    e = switches.effective(pc)
    assert e["max_tries"] == 2000 and e["requests"] == {"method": "legacy_x5", "lib_order": "l.npy"} and e["rounds"]["n"] == 6
    assert "presets" not in e and e["_preset"] == "x" and pc["max_tries"] == 200            # input unchanged
    assert switches.effective(pc, "none")["max_tries"] == 200
    try:
        switches.effective(pc, "y")
    except ValueError as err:
        assert "pack.presets" in str(err)
    else:
        raise AssertionError("unknown preset accepted")
    sw = switches.load({"exclusion": [{"libs": ["R80"], "around": {"provider": "p", "params": {"r_nm": 5}}}],
                        "thinning": {"pilot_squeeze": True}})
    assert sw.thinning["pilot_squeeze"] and not sw.thinning["resplit"]
    et, around, rad = A.exclusion_tables(sw.exclusion, ["R80", "P1"], ["ACTB"])
    assert et.tolist() == [1, 0] and not around.any()
    Z = np.array([[0, 0, 0, 10.0, 0], [500, 0, 0, 10.0, 0]])
    assert len(A.exclusion_tile_spheres(Z, np.full(3, -50.0), np.full(3, 50.0))) == 1
    # shell cap of a volume-measured shell (yeast band): measured volume instead of area x width
    tot = {(0, 0): 1000.0, (0, 2): 1000.0}
    out, rep = A.shell_cap(tot, np.array([100.0]), np.array([3.0]), {0: (1,)}, {0: 1e6},
                           [dict(ri=2, area_nm2=0.0, d0=0.0, d1=20.0, add_radius=False, pore_limit=0, labels=(1,), volume_nm3=1e5)],
                           log=lambda *a: None)
    assert abs(out[(0, 2)] * 100 / 1e5 - 0.1) < 1e-12 and rep[2][2] == 0


# Effective pack: sections of the reproduction presets (md5 of the canonical JSON of switches.effective on the raw
# config, `_preset` removed), recorded at B1 (2026-11-17) = the values before the B1 mechanism and the yeast split.
# A change here changes a reproduction preset: only on purpose (and then update the number with the reason).
REPRO_PACK_MD5 = {("jurkat_v7_port", None): "cead23259b5d46793144ae7e856256ee",
                  ("jurkat", "legacy_v7"): "5946890e58cc630024d4d2e939689677",   # v9 2026-11-23: organisation off (was 533b53b6...); freeze fix 2026-11-19: + organisation.field_cap_regions (was b2cdc84e...; B2: + organisation, was 0296332c...)
                  ("testcell", "legacy_v7"): "f89bbca208dcbeb9e5172c6d5b96fef5",
                  ("yeast", None): "1fa4390edf952330d7531ae341aa090c",   # u2 2026-11-23: organisation off (was 200dec5d...); freeze fix 2026-11-19: + organisation.field_cap_regions (was 7edc6905...); B1 2026-11-18: C3 pin (was ac4fdfb1... B2: AF3 pool off; 0e58bf80... = yeast_x5_core)
                  ("yeast_x5_core", None): "fdb7204589bf29880f22b422cab8b894",   # freeze fix 2026-11-19: field_cap_regions [] spelled out (= the default; was 0e58bf80...)
                  ("yeast_x5_port", None): "d53e76b814e6a6421a4eee5110cd73b8"}


def _eff_md5(cell, preset=None):
    import hashlib
    import json
    from wholecell.config import load_raw
    e = switches.effective(load_raw(cell).get("pack") or {}, preset)
    e.pop("_preset", None)
    if "organisation" in e and e["organisation"] is None:   # `organisation: null` = absent (switches.load: defaults, off;
        e.pop("organisation")                                # jurkat_v7_port turns off the B2 layer of jurkat.yaml so)
    return hashlib.md5(json.dumps(e, sort_keys=True, default=str).encode()).hexdigest(), e


def test_common_defaults_mechanism():
    """B1: COMMON_DEFAULTS as the built-in preset common_defaults for unpinned cells; pins; preset lists; benchmark
    settings (criteria docs/unification_criteria/B1_packer.md Part A3 / C)"""
    from wholecell.config import load_raw
    alg = set(switches.ALGORITHM_KEYS)
    assert set(switches.COMMON_DEFAULTS) <= alg, f"COMMON_DEFAULTS may only set {sorted(alg)}"
    switches.load(switches.effective({}), ["cytosol"])                                   # the defaults validate
    # unpinned -> common_defaults; legacy_v7 / none -> nothing merged; lists merge in order; config presets win
    pc = {"max_tries": 200, "regions": [], "presets": {"a": {"max_tries": 7, "rounds": {"n": 2}}, "b": {"max_tries": 9}}}
    e = switches.effective(pc)
    assert e["_preset"] == "common_defaults"
    exp = switches._merge({k: v for k, v in pc.items() if k != "presets"}, switches.COMMON_DEFAULTS)
    assert {k: v for k, v in e.items() if k != "_preset"} == exp
    for pin in ("legacy_v7", "none"):
        e = switches.effective(dict(pc, preset=pin))
        assert e["max_tries"] == 200 and "presets" not in e and "rounds" not in e
    assert switches.effective(pc, ["a", "b"])["max_tries"] == 9 and switches.effective(pc, "b,a")["max_tries"] == 7
    assert switches.effective(dict(pc, preset=["a", "b"]))["_preset"] == "a,b"
    assert switches.effective(dict(pc, presets={"legacy_v7": {"give_up": 3}}), "legacy_v7")["give_up"] == 3
    try:
        switches.effective(pc, ["a", "legacy_x6"])
    except ValueError as err:
        assert "legacy_x6" in str(err)
    else:
        raise AssertionError("unknown preset in a list accepted")
    # benchmark settings: algorithm keys only, every one of them set (shapes only in S4p), and they load
    for name, s in switches.B1_BENCHMARK.items():
        assert set(s) <= alg and set(s) >= alg - {"shapes"}, name
        sw = switches.load(switches.effective({"thin_lumens": {}}, ["legacy_v7", name]), ["cytosol"])
        assert sw.budgets["per_region"] == {}
    B = switches.B1_BENCHMARK
    assert switches.load(switches.effective({}, "b1_s1")).legacy
    assert B["b1_s2"]["rounds"]["n"] == 6 and not B["b1_s2"]["thinning"]["enabled"]
    assert B["b1_s3"]["thinning"]["enabled"] and B["b1_s3"]["shell_cap"]["enabled"] and B["b1_s3"]["squeeze"]["scope"] == "none"
    assert B["b1_s4"]["squeeze"]["scope"] == "underfilled" and B["b1_s4"]["thinning"]["pilot_squeeze"]
    assert B["b1_s5"]["squeeze"]["scope"] == "all" and B["b1_s4y"]["max_tries"] == 2000 and B["b1_s4p"]["shapes"]["source"] == "pack_shapes"
    # pins of the cell configs at B1
    pins = {c: (load_raw(c).get("pack") or {}).get("preset") for c in
            ("jurkat", "testcell", "jurkat_v7_port", "yeast", "yeast_x5_core", "yeast_x5_port")}
    assert pins == {"jurkat": None, "testcell": None, "jurkat_v7_port": "legacy_v7", "yeast": ["legacy_x5", "b1_s4y"],
                    "yeast_x5_core": "legacy_x5", "yeast_x5_port": "none"}, pins
    yp = load_raw("yeast")["pack"]["presets"]["legacy_x5"]
    assert set(yp) <= alg, f"yeast legacy_x5 holds non-algorithm keys {sorted(set(yp) - alg)} (cell keys go in pack:)"
    # reproduction presets unchanged; unpinned cells unchanged while COMMON_DEFAULTS is empty
    for (cell, pre), md5 in REPRO_PACK_MD5.items():
        assert _eff_md5(cell, pre)[0] == md5, f"effective pack: of {cell} (preset {pre}) changed"
    if not switches.COMMON_DEFAULTS:
        for cell in ("jurkat", "testcell"):
            assert _eff_md5(cell)[1] == _eff_md5(cell, "legacy_v7")[1]
    # yeast benchmark settings ride on legacy_x5: cell keys kept, algorithm keys replaced
    _, e = _eff_md5("yeast", ["legacy_x5", "b1_s1"])
    sw = switches.load(e)
    assert sw.legacy is False and sw.rounds["n"] == 1 and sw.squeeze["scope"] == "none" and sw.budgets["per_region"] == {}
    assert e["max_tries"] == 200 and not sw.organisation["enabled"] and sw.shapes["source"] == "pack_shapes" and sw.exclusion
    assert sw.thinning["resplit"] and not sw.thinning["enabled"]
    _, e = _eff_md5("yeast_x5_core")  # pinned legacy_x5: X5 values
    sw = switches.load(e)
    assert sw.squeeze["scope"] == "all" and sw.thinning["enabled"] and len(sw.budgets["per_region"]) == 11
    _, e = _eff_md5("yeast")          # pinned C3 since B1 (2026-11-18): targeted squeeze, 2000 / 1200 everywhere
    sw = switches.load(e)
    assert sw.squeeze["scope"] == "underfilled" and sw.thinning["enabled"] and sw.rounds["n"] == 6 and e["max_tries"] == 2000


# ------------------------------------------------------------------ numba kernels on a synthetic tile
_T = {}


def _tile():
    """kernel inputs for a 120 nm cube of testcell cytosol (no membranes in it), two library species"""
    if _T:
        return _T
    from wholecell.config import load
    from wholecell.pack import geomload
    from wholecell.pack.regions import Regions
    cfg = load("testcell")
    g = geomload.load(cfg); cell = g.cell; reg = Regions(cfg["pack"], g)
    size, vox = 120.0, 4.0
    cyto = reg.RNAMES.index("cytosol")
    rng = np.random.default_rng(3)
    R = float(cell.G[0])
    for _ in range(2000):                                  # a cube wholly in cytosol, 40 nm from any membrane
        lo = rng.uniform(-0.6 * R, 0.6 * R, 3)
        P = np.stack(np.meshgrid(*[lo[a] + (np.arange(int(size / vox)) + 0.5) * vox for a in range(3)], indexing="ij"), -1).reshape(-1, 3)
        Pe = np.stack(np.meshgrid(*[lo[a] - 40 + np.arange(11) * 20.0 for a in range(3)], indexing="ij"), -1).reshape(-1, 3)
        if np.all(np.isin(cell.labels(P), reg.VOL_REG[cyto][1])) and np.all(np.isin(cell.labels(Pe), reg.VOL_REG[cyto][1])) \
                and np.abs(cell.mems(Pe)).min() > 40:
            break
    else:
        raise RuntimeError("no cytosol cube found")
    libs = ["ACTB", "ALDOA_x4"]
    so = [0]; C, Rr, B = [], [], []
    for li in libs:
        z = np.load(os.path.join(LIBDIR, li + ".npz")); c = z["st_c"].astype(np.float64); r = z["st_r"].astype(np.float64)
        C.append(c); Rr.append(r); so.append(so[-1] + len(r)); B.append(float((np.linalg.norm(c, axis=1) + r).max()))
    _T.update(cfg=cfg, cell=cell, reg=reg, lo=lo, size=size, vox=vox, VOXP=P, cyto=cyto, so=np.array(so, np.int64),
              STC=np.concatenate(C), STR=np.concatenate(Rr), BR=np.array(B))
    return _T


def _pack(jn, REG_TRIES=None, REG_GIVEUP=None, EXC=None, EXC_TYPE=None, seed=11, BAND_M=None, BAND_W=None):
    from wholecell.pack import pack3d as p3
    T = _tile(); cell, reg = T["cell"], T["reg"]
    e = np.zeros(0, np.int64)
    nv = reg.NV
    vox_off = np.array([0] * (T["cyto"] + 1) + [len(T["VOXP"])] * (nv - T["cyto"]) + [len(T["VOXP"])] * len(reg.MEM_REG), np.int64)
    jl = np.arange(len(jn), dtype=np.int64); jr = np.full(len(jn), T["cyto"], np.int64)
    nsph = int((np.diff(T["so"])[jl] * np.asarray(jn)).sum()); cap_s = nsph + 1000
    return p3.pack_tile(seed, T["lo"], T["size"], 40.0, cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC,
                        T["so"], T["STC"], T["STR"], T["BR"], reg.regmask(), jl, jr, np.asarray(jn, np.int64), nv, vox_off,
                        T["VOXP"], T["vox"], np.zeros(len(reg.MEM_REG) + 1, np.int64), np.zeros((0, 9)), np.zeros(1),
                        np.zeros((0, 3)), np.zeros(0), 0, 1500.0, 0, 200, 25, cap_s, cap_s * 14,
                        e, np.zeros((0, 3)), e, e, np.zeros(0), e, np.r_[np.zeros(7), 10.0, 16.0], np.zeros((0, 4)),
                        reg.regmem(), *reg.shell_arrays(),
                        e if REG_TRIES is None else REG_TRIES, e if REG_GIVEUP is None else REG_GIVEUP,
                        np.zeros((0, 5)) if EXC is None else EXC, e if EXC_TYPE is None else EXC_TYPE,
                        e if BAND_M is None else BAND_M, np.zeros(0) if BAND_W is None else BAND_W,
                        reg.membrane_codes(cell.G, p3.c3))


def test_in_exc_reference():
    from wholecell.pack.pack3d import _in_exc
    rng = np.random.default_rng(0)
    EXC = np.c_[rng.uniform(-50, 50, (20, 3)), rng.uniform(5, 20, 20), rng.integers(0, 3, 20)]
    for _ in range(2000):
        t = rng.uniform(-60, 60, 3); bits = int(rng.integers(0, 8))
        ref = bool(np.any((((EXC[:, :3] - t) ** 2).sum(1) < EXC[:, 3] ** 2) & (((bits >> EXC[:, 4].astype(int)) & 1) == 1)))
        assert _in_exc(t, EXC, bits) == ref == _in_exc.py_func(t, EXC, bits)


def test_kernel_switch_arrays():
    T = _tile(); nreg = len(T["reg"].RNAMES)
    a = _pack([60, 30])
    b = _pack([60, 30], REG_TRIES=np.zeros(nreg, np.int64), REG_GIVEUP=np.zeros(nreg, np.int64),
              EXC=np.zeros((0, 5)), EXC_TYPE=np.zeros(2, np.int64), BAND_M=np.full(nreg, -1, np.int64), BAND_W=np.zeros(nreg))
    for x, y in zip(a[:4], b[:4]):
        assert x.dtype == y.dtype and np.array_equal(x, y)
    assert a[3].sum() > 0
    c = T["lo"] + T["size"] / 2
    ex = _pack([60, 30], EXC=np.array([[c[0], c[1], c[2], 1e4, 0.0]]), EXC_TYPE=np.array([1, 0], np.int64))
    assert ex[3][0] == 0 and ex[3][1] > 0, ex[3]
    tr = np.zeros(nreg, np.int64); tr[T["cyto"]] = 1
    bud = _pack([4000, 0], REG_TRIES=tr, REG_GIVEUP=tr)
    full = _pack([4000, 0])
    assert bud[3][0] <= full[3][0] and bud[3][0] < 4000, (bud[3], full[3])


def test_squeeze():
    from wholecell.pack import squeeze as S
    T = _tile(); cell = T["cell"]; reg = T["reg"]
    jn = np.array([3000, 400])
    types, pos, quat, placed = _pack(jn)[:4]
    assert placed.sum() < jn.sum()                          # over-asked: RSA leaves a deficit
    args = (cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC)
    # relax / accept: numba vs their pure-Python source on a few molecules
    k = np.r_[0:6, placed[0]:placed[0] + 4]
    typ = types[k].astype(np.int64); P0 = T["lo"] + T["size"] / 2 + np.random.default_rng(1).normal(0, 4.0, (len(k), 3))   # overlapping
    Q0 = quat[k].astype(np.float64); regs = np.full(len(k), T["cyto"], np.int64); mov = np.ones(len(k), bool)
    out = []
    for fn in (S.relax, S.relax.py_func):
        P, Q = P0.copy(), Q0.copy()
        fn(T["lo"], T["size"], 40.0, *args, T["so"], T["STC"], T["STR"], reg.regmask(), typ, P, Q, regs, mov,
           np.zeros((0, 3)), np.zeros(0), 20, S.H_MIN, S.STEP_MAX, S.ROT_MAX, S.WALL_EVERY)
        out.append((P, Q))
    assert np.allclose(out[0][0], out[1][0], atol=1e-6) and np.allclose(out[0][1], out[1][1], atol=1e-8)
    assert not np.allclose(out[0][0], P0)
    cs = 20000
    ka = S.accept(T["lo"], T["size"], 40.0, *args, T["so"], T["STC"], T["STR"], T["BR"], reg.regmask(), typ, P0, Q0, regs, mov,
                  np.zeros((0, 3)), np.zeros(0), cs, cs * 14, np.zeros(0, np.int64), np.zeros(0))
    kb = S.accept.py_func(T["lo"], T["size"], 40.0, *args, T["so"], T["STC"], T["STR"], T["BR"], reg.regmask(), typ, P0, Q0, regs, mov,
                          np.zeros((0, 3)), np.zeros(0), cs, cs * 14, np.zeros(0, np.int64), np.zeros(0))
    assert np.array_equal(ka, kb)
    # squeeze_tile: more copies, still exact
    jl = np.arange(2, dtype=np.int64); jr = np.full(2, T["cyto"], np.int64)
    vox_off = np.array([0] * (T["cyto"] + 1) + [len(T["VOXP"])] * (reg.NV - T["cyto"]), np.int64)
    t2, p2, q2, pl2 = S.squeeze_tile(np.random.default_rng(5), T["lo"], T["size"], 40.0, cell, T["so"], T["STC"], T["STR"],
                                     T["BR"], reg.regmask(), jl, jr, jn, placed, types, pos, quat, vox_off, T["VOXP"], T["vox"],
                                     np.zeros((0, 3)), np.zeros(0), np.ones(2, bool), 2, 30, 0.15)
    assert pl2.sum() >= placed.sum() and len(t2) == pl2.sum() and np.all(pl2 <= jn)
    assert np.all(np.diff(np.repeat(np.arange(2), pl2)) >= 0) and np.array_equal(t2, np.repeat(jl, pl2).astype(t2.dtype))
    cs = int(np.diff(T["so"])[t2].sum()) + 1000
    keep = S.accept(T["lo"], T["size"], 40.0, *args, T["so"], T["STC"], T["STR"], T["BR"], reg.regmask(), t2.astype(np.int64),
                    p2.astype(np.float64), q2.astype(np.float64), np.full(len(t2), T["cyto"], np.int64), np.ones(len(t2), bool),
                    np.zeros((0, 3)), np.zeros(0), cs, cs * 14, np.zeros(0, np.int64), np.zeros(0))
    assert keep.all(), f"{(~keep).sum()} squeezed molecules overlap or leave their region"
    print(f"  squeeze: RSA placed {placed.sum()}, after squeeze {pl2.sum()} of {jn.sum()}")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("pack switch tests OK")
