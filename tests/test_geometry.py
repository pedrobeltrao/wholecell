"""Fast geometry tests (engine smoke tests, rule table, python reference vs numba, testcell package round trip,
determinism, geometry stage).

    PYTHONHASHSEED=0 python tests/test_geometry.py        (or pytest tests/test_geometry.py)

The B1 equality test against legacy cell3d (2.55 M points) is tests/check_b1_geometry.py (--engine cellgeom_eval_v1).
"""
import json
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import yaml  # noqa: E402

from wholecell.config import load as load_config  # noqa: E402
from wholecell.geometry import cellgeom  # noqa: E402
from wholecell.geometry import eval_ref, rules  # noqa: E402
from wholecell.geometry.engines import cellgeom_eval_v1 as ev  # noqa: E402
from wholecell.geometry.engines import sdf_sections_v1 as e  # noqa: E402

SPEC = os.path.join(REPO, "toolkits", "analytic", "specs", "testcell.yaml")
SMALL_VAL = dict(n_rule=40_000, n_mc=100_000, n_area=400_000)


def _sphere_field(r, centre=(0, 0, 0), h=16.0, pad=200.0, band=(64.0, 64.0)):
    from toolkits.analytic.generate import bricks
    lo = np.asarray(centre, float) - r - pad
    n = int(np.ceil(2 * (r + pad) / h)) + 1
    X = np.stack(np.meshgrid(*[lo[a] + np.arange(n) * h for a in range(3)], indexing="ij"), -1).reshape(-1, 3)
    d = (np.linalg.norm(X - np.asarray(centre), axis=1) - r).reshape(n, n, n)
    return bricks(d, lo, h, 8, band)


def test_engine_smoke():
    cell = _sphere_field(1000.0)
    nuc = _sphere_field(400.0, band=(120.0, 120.0))
    fields = [cell, nuc] + [e.empty_field() for _ in range(5)]
    params = dict(e.DEFAULT_PARAMS, PERI_W=30.0, HALF=2.5)
    G = e.build_G(1200.0, fields, params)
    assert G[e.G_LAYOUT] == e.LAYOUT_ID and e.peri_w(G) == 30.0 and e.half(G) == 2.5
    assert e.params_from_G(G)["PERI_W"] == 30.0
    pores = np.array([[400.0, 0, 0, 1.0, 0, 0]])
    c = e.Cell(G, pores)
    P = np.array([[0, 0, 0], [700, 0, 0], [1500, 0, 0], [0, 400, 0], [0, 0, 600], [390, 0, 0], [0, 0, -400]], float)
    lab = c.labels(P)
    assert list(lab[:5]) == [e.NUCP, e.CYTO, e.EXT, e.PERI, e.CYTO], lab
    assert lab[5] == e.PORE and lab[6] == e.PERI
    M = c.mems64(P)
    assert M[1, e.M_PM] == 64.0 and M[1, e.M_ER] == 64.0               # inside clamp of cell; empty field = outside clamp
    assert abs(c.mems64(np.array([[980.0, 0, 0]]))[0, e.M_PM] - 20.0) < 0.5
    out = np.zeros(e.NMEM); e.mem_sdfs(700.0, 0.0, 0.0, G, pores, c.pidx, c.PT, c.PQ, c.PC, out)
    assert abs(e.min_mem_dist_g(out, -1, G) - (min(abs(out) - G[e.S_MEM_HALF:e.S_MEM_HALF + 8]) + 2.5)) < 1e-12
    # split / build round trip and a legacy-layout G
    G2 = e.build_G(**e.split_legacy_G(G, g_n=e.G_N), params=params)
    assert np.array_equal(G, G2)
    leg = np.r_[1200.0, np.zeros(7), 7.0]
    off = len(leg) + 12 * 7
    hdr, secs = [], []
    for h10, idx, dat in fields:
        hdr.append(np.r_[h10, off, off + len(idx)]); secs += [idx, dat]; off += len(idx) + len(dat)
    Gl = np.concatenate([leg] + hdr + secs)
    assert np.array_equal(e.build_G(**e.split_legacy_G(Gl), params=params), G)
    try:
        e.Cell(Gl, pores)
        raise AssertionError("legacy G accepted")
    except ValueError:
        pass


def test_engine_shim():
    from wholecell.geometry import engine as c3
    assert c3.ACTIVE_ENGINE == "cellgeom_eval_v1" and c3.NMEM == 8 and c3.M_OUTER == 7 and hasattr(c3, "min_mem_dist_g")


def _smoke_geometry():
    cell = _sphere_field(1000.0)
    nuc = _sphere_field(400.0, band=(120.0, 120.0))
    mito = _sphere_field(150.0, centre=(650.0, 0, 0), band=(64.0, 64.0))
    fields = [cell, nuc, mito] + [e.empty_field() for _ in range(4)]
    params = dict(e.DEFAULT_PARAMS, PERI_W=30.0, HALF=2.5)
    pores = np.array([[400.0, 0, 0, 1.0, 0, 0]])
    return fields, params, pores


def test_eval_engine_smoke():
    """cellgeom_eval_v1 with the default rule table == sdf_sections_v1 (labels and membrane SDFs, exact)"""
    fields, params, pores = _smoke_geometry()
    G = ev.build_G(1200.0, fields, params, face=-900.0)
    Gl = e.build_G(1200.0, fields, params, face=-900.0)
    assert G[ev.G_LAYOUT] == ev.LAYOUT_ID and G[ev.G_TABLES] > 0 and ev.peri_w(G) == 30.0
    assert np.array_equal(G[:int(G[ev.G_TABLES])][ev.G_N:], Gl[e.G_N:])     # same layout, tables appended
    c, cl = ev.Cell(G, pores), e.Cell(Gl, pores)
    P = np.array([[0, 0, 0], [700, 0, 0], [1500, 0, 0], [0, 400, 0], [0, 0, 600], [390, 0, 0], [0, 0, -400],
                  [650, 0, 0], [650, 0, 146], [0, 0, -950]], float)
    lab = c.labels(P)
    assert list(lab) == [ev.NUCP, ev.MATRIX, ev.EXT, ev.PERI, ev.CYTO, ev.PORE, ev.PERI, ev.MATRIX, ev.IMS, ev.OUTV], lab
    rng = np.random.default_rng(1)
    Q = rng.uniform(-1200, 1200, (50_000, 3))
    assert np.array_equal(c.labels(Q), cl.labels(Q)) and np.array_equal(c.mems64(Q), cl.mems64(Q))
    try:
        e.Cell(G, pores)
        raise AssertionError("cellgeom_eval_v1 G accepted by the legacy Cell")
    except ValueError:
        pass
    G2 = ev.build_G(**ev.split_legacy_G(G, g_n=ev.G_N), params=params)
    assert np.array_equal(G, G2)


def test_rule_table():
    """the rule YAML compiles; requires-dropping; malformed tables are refused"""
    rt = rules.load_rules(rules.default_rules_path("sdf_sections_v1"))
    assert [m["name"] for m in rt["membranes"]] == ev.MEM_NAMES
    allp = {"field:" + f for f in ev.FIELD_NAMES} | {"section:speckle", "section:golgi_cis", "section:cristae",
                                                     "clip:open_face"}
    M, L = rules.compile_tables(rt, e.DEFAULT_PARAMS, allp)
    assert M.shape == (ev.NMEM, ev.MROW) and L.shape == (15, ev.RROW)
    assert list(M[:, 0]) == [ev.K_FIELD, ev.K_SHELL_PORES, ev.K_FIELD, ev.K_FIELD, ev.K_INNER_CRISTAE, ev.K_GOLGI_CIS,
                             ev.K_FIELD, ev.K_WALL]
    assert M[0, 2] == -1.0 and list(M[1, 3:6]) == [34.0, 45.0, 150.0] and M[4, 3] == 8.0
    r = L[3]                                    # pore channel just inside the INM: 3 conditions
    assert r[0] == ev.PORE and r[1] == 3 and list(r[2:6]) == [ev.T_FIELD, ev.F_NE, ev.OP_LT, -17.0]
    assert list(r[6:10]) == [ev.T_FIELD, ev.F_NE, ev.OP_GT, -42.0] and list(r[10:14]) == [ev.T_PORE, ev.M_NE, ev.OP_LT, -1 + 1]
    assert L[-1, 0] == ev.CYTO and L[-1, 1] == 0
    # testcell-like package: no golgi / endo fields, no speckle / Golgi-cisternae sections, no open face
    small = {"field:cell", "field:ne_mid", "field:mito", "field:er", "field:nucleolus", "section:cristae"}
    rows = rules.active_rows(rt, small)
    assert [x["label"] for x in rows] == ["exterior", "perinuclear_space", "nuclear_pore_channel", "nucleolus",
                                         "nucleoplasm", "nuclear_pore_channel", "mito_matrix", "mito_ims", "er_lumen",
                                         "cytosol"]
    bad = [dict(rt, membranes=rt["membranes"][::-1]),
           dict(rt, labels=rt["labels"][:-1]),
           dict(rt, labels=[{"label": "cytosol", "if": ["nonsense < 0"]}, {"label": "cytosol"}]),
           dict(rt, labels=[{"label": "no_such_label"}]),
           dict(rt, membranes=[dict(m, of="ENDO") if m["name"] == "MITO_IM" else m for m in rt["membranes"]])]
    for b in bad:
        try:
            rules.compile_tables(b, e.DEFAULT_PARAMS, allp)
            raise AssertionError("malformed rule table accepted")
        except (ValueError, KeyError):
            pass


def _near_points(cell, n, rng, lo, hi):
    """half uniform in the box, half within 50 nm of a membrane (clamped values excluded)"""
    C = lo + (hi - lo) * rng.random((n * 40, 3))
    A = np.abs(cell.mems64(C)[:, :ev.NMEM - 1])
    near = C[(A < 50.0).any(1)][:n // 2]
    return np.vstack([lo + (hi - lo) * rng.random((n - len(near), 3)), near])


def _crosscheck(cell, n=10_000, seed=0):
    """python reference (rule YAML interpreted directly) vs numba kernel; numba vs legacy sdf_sections_v1"""
    rng = np.random.default_rng(seed)
    lo, hi = cellgeom._box(cell)
    P = _near_points(cell, n, rng, lo, hi)
    if cell.pores.shape[0]:
        k = rng.integers(0, len(cell.pores), n // 10)
        P = np.vstack([P, cell.pores[k, :3] + rng.normal(0, 40.0, (len(k), 3))])
    lab = cell.labels(P); M = cell.mems64(P)
    ref = eval_ref.RefCell(cell, cell.rules)
    lr = np.array([ref.label_of(p) for p in P]); Mr = np.array([ref.mem_sdfs(p)[0] for p in P])
    assert np.array_equal(lab, lr), int((lab != lr).sum())
    assert np.abs(M - Mr).max() < 1e-9, np.abs(M - Mr).max()
    leg = cellgeom.load(cell.package_dir, verify=False, evaluator="sdf_sections_v1", check_active_engine=False)
    assert np.array_equal(lab, leg.labels(P)) and np.array_equal(M, leg.mems64(P))
    return len(P), np.unique(lab)


def test_reference_crosscheck(tmp=None):
    """~10^4 points per package: python reference == numba (labels; SDF < 1e-9 nm) and == sdf_sections_v1 (exact)"""
    from wholecell.config import load_site
    site = load_site()
    for cid, rel in (("testcell", "testcell/analytic_v1"), ("jurkat", "jurkat/em_cell1_v7")):
        pkg = os.path.join(site["geom_root"], rel)
        if not os.path.isdir(pkg):
            print(f"SKIP: {pkg} missing"); continue
        cell = cellgeom.load_cfg(load_config(cid), verify="quick")
        assert cell.evaluator == "cellgeom_eval_v1" and cell.rules_path.endswith("sdf_sections_v1.yaml")
        n, codes = _crosscheck(cell)
        print(f"  {cid}: {n} points, labels {sorted(cell.code_to_name[int(c)] for c in codes)}")


YEAST_PKG = "yeast/analytic_x5_s0_f5f3a5c723b2"     # the yeast PORT build s0 package (20261116-aabe7b5-s0-x5repro)


def _yeast_vendored(pkg, no_fields=False):
    from plugins.yeast import geom_engine
    vc = geom_engine.load(pkg)
    import cell3d as c3
    if no_fields:
        d = json.load(open(os.path.join(pkg, "native", "geometry.json"))); d.pop("fields", None)
        vc = c3.Cell(d)
    return vc, c3


def _yeast_vendored_eval(vc, c3, P):
    allp = np.arange(len(vc.PT), dtype=np.int64)
    lab = c3.label_points(P, vc.G, vc.pores, allp, vc.PT, vc.PQ, vc.PC)
    M = np.zeros((len(P), c3.NMEM)); C = np.zeros(len(P)); out = np.zeros(c3.NMEM)
    for i, p in enumerate(P):
        c3.mem_sdfs(p[0], p[1], p[2], vc.G, vc.pores, allp, vc.PT, vc.PQ, vc.PC, out)
        M[i] = out; C[i] = c3.min_mem_dist(out, -1)
    return lab, M, C


def _yeast_points(cell, n, rng):
    R = float(cell.manifest["domain"]["enclosing_radius_nm"])
    C = rng.uniform(-R, R, (n * 20, 3))
    near = C[(np.abs(cell.mems64(C)) < 50.0).any(1)][:n // 2]
    k = rng.integers(0, len(cell.pores), n // 10)
    return np.vstack([rng.uniform(-R, R, (n - len(near), 3)), near, cell.pores[k, :3] + rng.normal(0, 40.0, (len(k), 3))])


def test_yeast_package(tmp=None):
    """cellgeom_eval_v1 + config/geometry/rules/yeast_cell3d_v1.yaml on the yeast s0 package vs the vendored cell3d
    (~10^4 points): labels identical; membrane SDFs and packer clearance identical under the vendored env's numba
    (0.61.2; full check tests/check_a1_yeast_geometry.py), within 1e-9 nm otherwise (numba 0.68 compiles cell3d's
    sdf_cer with a different fastmath reassociation, ulp level); python reference == numba (labels; SDF < 1e-9).
    Also the analytic fallbacks (package without fields) against cell3d without fields."""
    import numba
    from wholecell.config import load_site
    pkg = os.path.join(load_site()["geom_root"], YEAST_PKG)
    if not os.path.isdir(pkg):
        print(f"SKIP: yeast package {pkg} missing (port build s0)"); return
    cell = cellgeom.load(pkg, verify="quick")
    assert cell.evaluator == "cellgeom_eval_v1" and cell.rules_path.endswith("yeast_cell3d_v1.yaml")
    assert cell.n_mems == 10 and len(cell.PT) == cell.manifest["prims"]["n"]
    rng = np.random.default_rng(7)
    P = _yeast_points(cell, 10_000, rng)
    lab, M, Cl = cell.labels(P), cell.mems64(P), cell.clearance(P, -1)
    vc, c3 = _yeast_vendored(pkg)
    lv, Mv, Cv = _yeast_vendored_eval(vc, c3, P)
    exact = numba.__version__ == "0.61.2"
    tol = 0.0 if exact else 1e-9
    assert np.array_equal(lab, lv), int((lab != lv).sum())
    assert np.abs(M - Mv).max() <= tol and np.abs(Cl - Cv).max() <= tol, (np.abs(M - Mv).max(0), np.abs(Cl - Cv).max())
    assert len(np.unique(lab)) == 15, np.unique(lab)
    # python reference (rule YAML interpreted directly)
    ref = eval_ref.RefCell(cell, cell.rules)
    Q = P[rng.choice(len(P), 1500, replace=False)]
    lr = np.array([ref.label_of(p) for p in Q]); Mr = np.array([ref.mem_sdfs(p)[0] for p in Q])
    assert np.array_equal(cell.labels(Q), lr), int((cell.labels(Q) != lr).sum())
    assert np.abs(cell.mems64(Q) - Mr).max() < 1e-9, np.abs(cell.mems64(Q) - Mr).max()
    # analytic fallbacks: the same package without its fields (sphere cell, ellipsoid nucleus, no vacuole field)
    tmp = tmp or tempfile.mkdtemp(prefix="wc_yeast_nofields_")
    man = json.load(open(os.path.join(pkg, "manifest.json")))
    for n, f in man["fields"].items():
        man["fields"][n] = dict(slot=f["slot"], empty=True)
        man["files"].pop(f.get("file"), None)
    json.dump(man, open(os.path.join(tmp, "manifest.json"), "w"))
    for rel in ("native", "prims.npz", "pores.npz"):
        os.symlink(os.path.join(pkg, rel), os.path.join(tmp, rel))
    c0 = cellgeom.load(tmp, verify=False)
    assert len(c0.G[int(c0.G[ev.G_TABLES]):]) > 0 and int(c0.G[int(c0.G[ev.G_TABLES]) + 3]) == 2   # 2 fallbacks
    v0, _ = _yeast_vendored(pkg, no_fields=True)
    P0 = _yeast_points(c0, 4_000, rng)
    l0, M0, C0 = _yeast_vendored_eval(v0, c3, P0)
    assert np.array_equal(c0.labels(P0), l0), int((c0.labels(P0) != l0).sum())
    assert np.abs(c0.mems64(P0) - M0).max() <= tol and np.abs(c0.clearance(P0) - C0).max() <= tol
    ref0 = eval_ref.RefCell(c0, c0.rules)
    Q0 = P0[:800]
    assert np.array_equal(c0.labels(Q0), np.array([ref0.label_of(p) for p in Q0]))
    assert np.abs(c0.mems64(Q0) - np.array([ref0.mem_sdfs(p)[0] for p in Q0])).max() < 1e-9
    print(f"  yeast: {len(P)} points, 15 labels, max|dSDF| vs cell3d {np.abs(M - Mv).max():.2g} (numba "
          f"{numba.__version__}); no-field fallback {len(P0)} points ok")


def test_testcell_roundtrip_and_determinism(tmp=None):
    from toolkits.analytic.generate import generate
    spec = yaml.safe_load(open(SPEC))
    tmp = tmp or tempfile.mkdtemp(prefix="wc_geom_test_")
    dest = os.path.join(tmp, "testcell_rt")
    man, res = generate(spec, dest, validate_kw=SMALL_VAL)
    assert res["pass"], res["failed_checks"]
    cell = cellgeom.load(dest)
    assert cell.package_id == "testcell/analytic_v1" and cell.pores.shape == (30, 6)
    rng = np.random.default_rng(0)
    P = rng.uniform(-1300, 1300, (20_000, 3))
    names = np.array([cell.code_to_name[int(c)] for c in cell.labels(P)], dtype=object)
    assert (names == cellgeom.eval_label_rule(cell, P)).all()
    assert cell.summary["open_face"]["f_vol"] == 1.0 and cell.summary["cell"]["vol_sdf_um3"] > 5.0
    # determinism: same data as the published package (same spec + seed). Integer arrays must be identical; float
    # arrays may differ in the last bits between CPU generations (pores: 1.4e-14 nm on eu-g4 vs the node that wrote
    # the package, 2026-11-16), so they are compared to 1e-9 nm. Identical sha256 is the usual case, not required.
    cfg = load_config("testcell")
    pub = cfg["geometry"]["package"]
    if os.path.isdir(pub):
        for rel in man["files"]:
            if rel.endswith(".npz"):
                a, b = np.load(os.path.join(dest, rel)), np.load(os.path.join(pub, rel))
                assert sorted(a.files) == sorted(b.files), rel
                for k in a.files:
                    x, y = a[k], b[k]
                    ok = x.shape == y.shape and (np.allclose(x, y, rtol=0, atol=1e-9) if x.dtype.kind == "f"
                                                 else np.array_equal(x, y))
                    assert ok, f"{rel}:{k} differs from the published package {pub}"
    # immutability: a second write to the same folder is refused
    try:
        cellgeom.PackageWriter(dest)
        raise AssertionError("existing package folder accepted")
    except FileExistsError:
        pass


def test_published_testcell_and_stage(tmp=None):
    cfg = load_config("testcell")
    g = cfg["geometry"]
    if not os.path.isdir(g["package"]):
        print("SKIP: published testcell package missing"); return
    from wholecell.geometry import stage
    cell = cellgeom.load(g["package"])
    assert cell.manifest_sha256 == g["manifest_sha256"]
    tmp = tmp or tempfile.mkdtemp(prefix="wc_geom_stage_")
    pj = stage.run(cfg, tmp, mc_points=100_000)
    assert pj["manifest_sha256"] == g["manifest_sha256"]
    v = json.load(open(os.path.join(tmp, "geometry", "validation.json")))
    m = json.load(open(os.path.join(tmp, "geometry", "measures.json")))
    assert v["pass"] and v["reused"]
    assert abs(m["region_volumes_mc"]["volume_um3"]["cytosol"] - 5.08) < 0.2
    assert "PM" in m["membrane_areas_summary"]


if __name__ == "__main__":
    os.environ.setdefault("PYTHONHASHSEED", "0")
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name, flush=True)
