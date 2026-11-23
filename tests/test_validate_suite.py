"""Tests of the common validation suite (wholecell/validate/suite.py, plan s4, A8 phase 1), under a minute on 1-3 CPUs:

    PYTHONHASHSEED=0 python tests/test_validate_suite.py        (or pytest tests/test_validate_suite.py)

  seam        synthetic tiles: uniform centres give index 1 (within sampling noise); an over-filled 10 nm face band gives
              the expected ratio; seam_index() on a fake build with tiles.parquet + tile npz files
  regression  compare() for every direction; no_regression() against a reference with every status (ok, regress,
              improve, n/a without reference / without value); plan s4 tolerances
  registry    every check C1-C17 is in CHECKS and gives at least one row (none silently skipped)
  testcell    (needs the latest testcell CI build) the suite runs into a temp folder, writes nothing into the build,
              reports n/a (with a reason) where data are absent, and C10 placed fraction equals placement.json; the
              orientation re-run (signs + tilt) reproduces placement.json's sign error rate
  phase 2     locmetrics.score = the legacy Jurkat loc_compare.score (random inputs) and = yeast localization.py:score
              (agree, top-1, mass-weighted) on random yeast-style sets; C3 AUROC = sklearn roc_auc_score; C6 copies_vs
"""
import glob
import json
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell.validate import suite as S  # noqa: E402


# ------------------------------------------------------------------ seam index
def _uniform(n, T, rng):
    return rng.uniform(0, T, (n, 3))


def test_seam_uniform():
    rng = np.random.default_rng(1)
    T = 200.0
    b, i = S.seam_counts(_uniform(400000, T, rng), T)
    x = S.seam_index(b, i, T)
    assert abs(x - 1) < 0.02, x


def test_seam_overfilled_band():
    rng = np.random.default_rng(2)
    T = 500.0
    p = _uniform(300000, T, rng)
    vb, vi = S.seam_volumes(T)
    # add extra centres uniformly inside the band so that the band density is 1.5x the bulk density
    dens = len(p) / T ** 3
    extra = []
    need = int(0.5 * dens * vb)
    while sum(len(e) for e in extra) < need:
        q = rng.uniform(0, T, (200000, 3))
        d = np.minimum(q, T - q).min(1)
        extra.append(q[d < S.SEAM_BAND_NM])
    p = np.vstack([p, np.vstack(extra)[:need]])
    b, i = S.seam_counts(p, T)
    x = S.seam_index(b, i, T)
    assert abs(x - 1.5) < 0.03, x
    # depleted band: half of the band centres removed
    q = _uniform(300000, T, rng)
    dq = np.minimum(q, T - q).min(1)
    q = q[(dq >= S.SEAM_BAND_NM) | (rng.uniform(size=len(q)) < 0.5)]
    b, i = S.seam_counts(q, T)
    assert abs(S.seam_index(b, i, T) - 0.5) < 0.03
    assert S.seam_counts(np.zeros((0, 3)), T) == (0, 0)


def test_seam_build():
    """seam() on a fake build: 4 cytosol tiles (one below the volume threshold), global coordinates with lo0"""
    rng = np.random.default_rng(3)
    T = 100.0
    lo0 = np.array([-150.0, 0.0, 50.0])
    with tempfile.TemporaryDirectory() as d:
        os.makedirs(os.path.join(d, "pack", "tiles"))
        json.dump(dict(TILE=T, lo0=lo0.tolist()), open(os.path.join(d, "pack", "tiles_meta.json"), "w"))
        rows = []
        for t, ijk, frac in ((0, (0, 0, 0), 1.0), (1, (1, 0, 0), 0.99), (2, (0, 1, 1), 0.96), (3, (2, 2, 2), 0.5)):
            pos = lo0 + np.array(ijk) * T + rng.uniform(0, T, (60000, 3))
            if t == 3:                     # not a seam tile: an empty interior would dominate if it were used
                pos = pos[np.minimum(pos - lo0 - np.array(ijk) * T, T - (pos - lo0 - np.array(ijk) * T)).min(1) < 10]
            np.savez(os.path.join(d, "pack", "tiles", f"{t}.npz"), pos=pos.astype(np.float32))
            rows.append(dict(tile=t, ti=ijk[0], tj=ijk[1], tk=ijk[2], cytosol=frac * T ** 3))
        pd.DataFrame(rows).to_parquet(os.path.join(d, "pack", "tiles.parquet"))
        ctx = type("C", (), dict(sc={}, pack=lambda self, n: os.path.join(d, "pack", n)))()
        r = S.seam(ctx)
        assert r["n_tiles"] == 3 and r["n_qualifying"] == 3, r
        assert abs(r["index"] - 1) < 0.03, r
        r = S.seam(ctx, min_frac=0.99)
        assert r["n_tiles"] == 2, r


# ------------------------------------------------------------------ no-regression
def test_compare_directions():
    c = S.compare
    assert c(0.80, 0.82, "up", 0.01) == "regress"
    assert c(0.815, 0.82, "up", 0.01) == "ok"
    assert c(0.84, 0.82, "up", 0.01) == "improve"
    assert c(3.0, 1.0, "down", 1.0) == "regress"
    assert c(1.9, 1.0, "down", 1.0) == "ok"
    assert c(0.0, 1.5, "down", 1.0) == "improve"
    assert c(1.10, 1.02, "target:1", 0.02) == "regress"
    assert c(0.97, 1.02, "target:1", 0.02) == "ok"           # |d| 0.03 vs 0.02: within tolerance
    assert c(1.00, 1.05, "target:1", 0.02) == "improve"
    assert c(160.0, 163.0, "stable", 2.0) == "regress"
    assert c(164.0, 163.0, "stable", 2.0) == "ok"
    assert c(5.0, 1.0, "report", None) == "ok"
    assert c(None, 1.0, "up", 0.01) == "n/a"
    assert c(1.0, None, "up", 0.01) == "n/a"


def test_tolerances_plan_s4():
    T = {k: v[0] for k, v in S.TOL.items()}
    assert T["loc"](0.8) == 0.01
    assert T["rate1e6"](0.36) == 1.0 and T["rate1e6"](50.0) == 5.0
    assert T["ratefrac"](2.5e-5) == 2.5e-6 and T["ratefrac"](1e-6) == 1e-6
    assert T["density"](160) == 2.0 and T["nn"](28) == 0.5 and T["placed"](0.96) == 0.003


def test_no_regression_table():
    mk = lambda m, v, d, t: S.row("C1", m, v, direction=d, tol=t)
    ref = [mk("a", 0.83, "up", "loc"), mk("b", 0.0, "down", "rate1e6"), mk("c", 162.8, "stable", "density"),
           mk("d", 1.04, "target:1", "seam"), mk("e", 0.96, "up", "placed")]
    cur = [mk("a", 0.81, "up", "loc"), mk("b", 0.5, "down", "rate1e6"), mk("c", 166.0, "stable", "density"),
           mk("d", 1.00, "target:1", "seam"), mk("e", 0.958, "up", "placed"), mk("f", 1.0, "up", "loc"),
           S.na("C1", "g", "n/a (no data): test")]
    table, rows = S.no_regression(cur, ref, "ref")
    st = {t["metric"]: t["status"] for t in table}
    assert st == dict(a="regress", b="ok", c="regress", d="improve", e="ok", f="n/a", g="n/a"), st
    assert next(t for t in table if t["metric"] == "f")["reason"].startswith("n/a (metric has no value in the reference")
    assert next(r for r in rows if r["metric"] == "g")["reason"] == "n/a (no data): test"
    table, _ = S.no_regression([mk("a", 0.81, "up", "loc")], None, "none")
    assert table[0]["status"] == "n/a" and "no reference build" in table[0]["reason"]
    for t in table:
        assert t["doc"].startswith("docs/stages/validation.md#c")


def test_legacy_match_rounding():
    assert S._match(0.9643817, 0.964) and not S._match(0.9653, 0.964)
    assert S._match(0.1, 0.1) and not S._match(92, 186)


# ------------------------------------------------------------------ phase 2: one metric definition
def _legacy_jurkat_score(COMPS, SPECIAL, p, q, main_needed=True):
    """verbatim copy of the legacy Jurkat loc_compare.score (before A8 phase 2)"""
    sup = {c for c in q if c in COMPS}
    for sp, comps in SPECIAL.items():
        if sp in q: sup |= comps
    unsup = sum(v for c, v in p.items() if c not in sup)
    main = [c for c, w in q.items() if w >= 0.999 and c in COMPS]
    miss = (sum(p.get(c, 0.0) < 0.05 for c in main) / len(main)) if (main and main_needed) else 0.0
    qq = {c: w for c, w in q.items() if c in COMPS}; s = sum(qq.values())
    tvd = 0.5 * sum(abs(p.get(c, 0.0) - qq.get(c, 0.0) / s) for c in COMPS) if s else np.nan
    top = max(p, key=p.get) if p else None
    return dict(unsup=unsup, miss=miss, score=0.5 * unsup + 0.5 * miss, tvd=tvd, tophit=float(top in sup) if top else np.nan)


def _rand_p(rng, comps):
    k = rng.integers(1, 4)
    c = rng.choice(comps, k, replace=False)
    w = rng.dirichlet(np.ones(k)) if rng.uniform() < 0.8 else np.r_[1.0, np.zeros(k - 1)]
    return {cc: float(x) for cc, x in zip(c, w)}


def test_locmetrics_equals_legacy_jurkat():
    from wholecell.validate import locmetrics as LM
    rng = np.random.default_rng(5)
    comps = ["cytosol", "nucleoplasm", "ER", "PM", "mito", "golgi"]
    special = {"vesicles": {"golgi", "ER"}}
    for _ in range(3000):
        p = _rand_p(rng, comps)
        q = {c: float(rng.choice([1.0, 0.3, 0.6])) for c in rng.choice(comps + ["vesicles"], rng.integers(1, 4), replace=False)}
        a = _legacy_jurkat_score(comps, special, p, q); b = LM.score(p, q, comps, special)
        for k in a:
            assert (a[k] == b[k]) or (np.isnan(a[k]) and np.isnan(b[k])), (k, a, b)


def test_locmetrics_equals_yeast_formula():
    """yeast localization.py:score (share mode) agree / top1 / mass-weighted, re-implemented from its text"""
    from wholecell.validate import locmetrics as LM
    rng = np.random.default_rng(7)
    CL = ["cytosol", "nucleus", "nucleolus", "ER", "mitochondrion", "vacuole"]
    model, ref = {}, {}
    for g in range(400):
        sh = rng.dirichlet(np.ones(len(CL))) * (rng.uniform(size=len(CL)) < 0.5)
        if sh.sum() == 0:
            sh[0] = 1.0
        sh = sh / sh.sum()
        model[f"g{g}"] = ({c: float(x) for c, x in zip(CL, sh)}, float(rng.uniform(1e3, 1e8)))
        ref[f"g{g}"] = set(rng.choice(CL, rng.integers(1, 3), replace=False))
    agree, hit, w = [], [], []
    for g, r in ref.items():
        p, m = model[g]
        sh = np.array([p[c] for c in CL]); inref = np.array([c in r for c in CL])
        agree.append(1 - 0.5 * sh[~inref].sum() - 0.5 * (sh[inref] < 0.05).mean())
        hit.append(CL[int(np.argmax(sh))] in r); w.append(m)
    agree, w = np.array(agree), np.array(w)
    out, _ = LM.score_set(model, ref, CL)
    assert abs(out["agree"]["unweighted"] - agree.mean()) < 1e-12
    assert abs(out["agree"]["mass_weighted"] - (agree * w).sum() / w.sum()) < 1e-12
    assert abs(out["tophit"]["unweighted"] - np.mean(hit)) < 1e-12


def test_c3_c6_metrics():
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(3)
    model = rng.normal(0, 2, 500); meas = model * 0.5 + rng.normal(0, 1, 500)
    r = S.image_nc_metrics(meas, model)
    nucl, cyto = model > 1, model < -1
    y = np.r_[np.ones(nucl.sum()), np.zeros(cyto.sum())]
    assert abs(r["auroc_nuclear"] - roc_auc_score(y, np.r_[meas[nucl], meas[cyto]])) < 1e-12
    assert S.image_nc_metrics(meas, model - 10)["auroc_nuclear"] is None        # no model-nuclear genes: n/a
    a = pd.Series([12.0, 100, 1000, 0, 5]); b = pd.Series([20.0, 150, 900, 7, 0])
    c = S.copies_vs(a, b)
    assert c["n"] == 3 and abs(c["spearman"] - 1.0) < 1e-12 and c["within_2x"] == 1.0


# ------------------------------------------------------------------ registry
def test_registry_complete():
    assert list(S.CHECKS) == [f"C{i}" for i in range(1, 18)]


def _latest_testcell():
    from wholecell.config import load_site
    root = load_site()["results_root"]
    c = sorted(glob.glob(os.path.join(root, "testcell", "*-ci-w8")))
    c = [x for x in c if os.path.exists(os.path.join(x, "validate", "MANIFEST.tsv"))]
    return c[-1] if c else None


def test_testcell_build():
    b = _latest_testcell()
    if not b:
        print("  (no testcell CI build: skipped)")
        return
    from wholecell.config import load
    from wholecell import registry
    cfg = load("testcell"); registry.load_plugins(cfg)
    before = {p: os.path.getmtime(p) for p in glob.glob(os.path.join(b, "**"), recursive=True) if os.path.isfile(p)}
    with tempfile.TemporaryDirectory() as out:
        res = S.run(cfg, b, out=out, nproc=2)
        for f in ("suite.json", "no_regression.tsv", "legacy_check.tsv"):
            assert os.path.exists(os.path.join(out, f)), f
        after = {p: os.path.getmtime(p) for p in glob.glob(os.path.join(b, "**"), recursive=True) if os.path.isfile(p)}
        assert after == before, "the suite wrote into the build"
        M = res["metrics"]
        assert {r["check"] for r in M} == set(S.CHECKS)
        for r in M:
            assert r["status"] in S.STATUS, r
            if r["value"] is None:
                assert r["status"] == "n/a" and r["reason"].startswith("n/a ("), r
        by = {(r["check"], r["metric"]): r for r in M}
        for c in ("C1", "C12", "C13"):                       # testcell has no localisation / EM / chromatin references
            assert all(r["value"] is None for r in M if r["check"] == c), c
        P = json.load(open(os.path.join(b, "validate", "placement.json")))
        ss = P["species_share"]
        assert abs(by[("C10", "placed_copies_fraction")]["value"] - sum(v["placed"] for v in ss.values()) / sum(v["requested"] for v in ss.values())) < 1e-12
        assert by[("C7", "a.clash_rate_per_1e6")]["value"] == P["clash_rate_per_1e6"]
        # signs + tilt: from placement.json when the build already has tilt (builds after A8 phase 2), else re-run
        O = P if "tilt_median_deg" in P else json.load(open(os.path.join(out, "core", "orientation.json")))
        assert O["orientation_error_rate_per_1e6"] == P["orientation_error_rate_per_1e6"], O
        assert by[("C9", "median_tilt_deg")]["value"] == O["tilt_median_deg"]
        si = by[("C10", "seam_index")]
        assert si["value"] is None or 0.5 < si["value"] < 1.5, si


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("validate suite tests OK")
