# new in wholecell (freeze v8 E7, 2026-11-19)
"""PM topology fix (wholecell/library/topology.py): rotation, alternation rule, S0 flag rule, library step apply."""
import os
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from wholecell.library import topology as t  # noqa: E402
from wholecell.library import build  # noqa: E402


def _unit(path, z_ecto=+5.0, n=40):
    rng = np.random.default_rng(0)
    b = np.concatenate([rng.normal([0, 0, z_ecto], 1, (n, 3)), rng.normal([0, 0, 0], 0.5, (10, 3)),
                        rng.normal([0, 0, -z_ecto / 2], 0.5, (5, 3))]).astype(np.float32)
    np.savez_compressed(path, beads=b, kind=np.zeros(len(b), np.int8), st_c=b[::5].copy(), st_r=np.ones(len(b[::5]), np.float32),
                        st_c30=b[::9].copy(), st_r30=np.ones(len(b[::9]), np.float32), src=np.str_("pdb:test"), membrane=np.str_("TM"))
    return b


def test_rotation_proper_and_exact():
    with tempfile.TemporaryDirectory() as d:
        b = _unit(os.path.join(d, "u.npz"))
        z = np.load(os.path.join(d, "u.npz"))
        r = t.rotate(z)
        assert np.array_equal(r["beads"][:, 0], b[:, 0]) and np.array_equal(r["beads"][:, 1], -b[:, 1])
        assert np.array_equal(r["beads"][:, 2], -b[:, 2]) and r["beads"].dtype == np.float32
        assert np.array_equal(r["st_c"], b[::5] * np.array([1, -1, -1], np.float32))
        assert np.array_equal(r["st_c30"][:, 2], -b[::9, 2]) and np.array_equal(r["st_r"], z["st_r"])
        assert np.linalg.det(np.diag(t.ROT)) == 1.0                     # proper rotation: chirality kept
        rr = t.rotate(type("Z", (), dict(files=list(r), __getitem__=lambda s, k: r[k]))())
        assert np.array_equal(rr["beads"], b)                           # twice = identity


def test_alternation_trbc1():
    topo = 'TOPO_DOM 171..176; /note="Cytoplasmic"; /evidence="ECO:0000255"'
    tm = 'TRANSMEM 150..170; /note="Helical"; /evidence="ECO:0000255"'
    assert t.alternation_counts(topo, tm, 176) == (149, 6)
    # two TMs: N-term cyto annotated -> loop extracellular, C-term cyto
    assert t.alternation_counts('TOPO_DOM 1..10; /note="Cytoplasmic"', 'TRANSMEM 11..30; /note="Helical"; TRANSMEM 51..70; /note="Helical"', 80) == (20, 20)
    assert t.alternation_counts("", 'TRANSMEM 11..30; /note="Helical"', 80) is None


def test_balance_flag_rule():
    assert t.balance_flag(0.95, 0.08) and t.balance_flag(0.0, 1.0)
    assert not t.balance_flag(0.95, 0.92) and not t.balance_flag(0.6, 0.1) and not t.balance_flag(np.nan, 0.1)


def test_apply_only_listed_units():
    import pandas as pd
    with tempfile.TemporaryDirectory() as d:
        _unit(os.path.join(d, "A.npz")); _unit(os.path.join(d, "B.npz"), z_ecto=-5.0)
        mb = build._md5(os.path.join(d, "B.npz"))
        T = pd.DataFrame([dict(lib="A", rule="glyco_ectodomain_exoplasmic", genes="X", uniprot="", n_PM="10", evidence="test")])
        lo0 = t.bead_balance(np.load(os.path.join(d, "A.npz"))["beads"])[2]
        rec = t.apply(d, T, log=lambda *a: None)
        assert rec[0]["rotated"] and lo0 < 0.2 and rec[0]["obs_lo_after"] > 0.8      # ectodomain now on -z (exoplasmic)
        z = np.load(os.path.join(d, "A.npz"))
        assert int(z["topology_fixed"]) == 1 and str(z["topology_rule"]) == "glyco_ectodomain_exoplasmic"
        assert build._md5(os.path.join(d, "B.npz")) == mb                               # unlisted unit untouched
        try:
            t.apply(d, T, log=lambda *a: None)
            raise AssertionError("second apply must fail")
        except RuntimeError:
            pass
        bad = T.assign(lib="C")
        try:
            t.apply(d, bad, log=lambda *a: None)
            raise AssertionError("missing unit must fail")
        except FileNotFoundError:
            pass


def test_table_and_config():
    import yaml
    T = t.read_table(os.path.join(REPO, "plugins", "jurkat", "data", "library", "pm_topology_fix.tsv"))
    g = T[T.rule == "glyco_ectodomain_exoplasmic"]
    assert len(g) == 13 and g.n_PM.astype(int).sum() == 168347
    assert list(T[T.rule == "tm_alternation"].lib) == ["TRBC1_x2"]
    cfg = yaml.safe_load(open(os.path.join(REPO, "config", "cells", "jurkat.yaml")))
    steps = [s["step"] for s in cfg["library"]["steps"]]
    assert steps[-1] == "pm_topology_fix" and steps.index("adaptive_trees") == len(steps) - 2
    assert "pm_topology_fix" in build.STEPS


if __name__ == "__main__":
    for k, fn in list(globals().items()):
        if k.startswith("test_"):
            fn()
            print("ok", k)
