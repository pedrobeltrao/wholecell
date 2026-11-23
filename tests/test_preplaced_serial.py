"""Tests of the serial pre-placed runner (wholecell/preplaced/serial.py; unification A6), seconds on 1 CPU:

    python tests/test_preplaced_serial.py        (or pytest tests/test_preplaced_serial.py)

  order        builders run in the configured order on one context; each sees what the ones before it left (pre,
               obstacles, pools)
  layers       per-builder layer = the arrays it added or replaced (identity), written as layers/<name>.npz; the
               union of the layers in order reassembles ctx.pre
  consume      Context.consume = the stage-3 consume (proportional over regions, capped at the pool, log line)
  registry     every builder of config/cells/yeast.yaml preplaced.serial is registered (core + yeast plugin), and
               the yeast order is the stage-3 order (fixed list)
  legacy       the A0c rename "prepass" -> "preplaced": old config keys load (renamed, one message), old build files
               are found through wholecell.build.ref_path, the CLI accepts `prepass`, and scripts/ci_compare_builds.py
               matches an old and a new build file by file
The exactness against the yeast port build is checked by scripts/preplaced_a6_check.py (full builders, ~1 h).
"""
import os
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.preplaced import serial  # noqa: E402

YEAST_ORDER = ["core.npc_pores", "core.microtubules_mtoc", "yeast.kinetochores", "yeast.spb", "yeast.rabl_chromatin",
               "yeast.spindle_pools", "core.npc_receptors", "core.fg_mesh", "yeast.rdna_transcription", "yeast.rdna_chromatin",
               "yeast.pol2", "yeast.origins", "yeast.chromatin_sites", "core.microtubule_pools", "yeast.eisosomes", "yeast.mcc",
               "yeast.actin_cables", "yeast.endocytic_coat", "yeast.amphiphysin", "yeast.actin_patches", "core.atp_synthase_dimers",
               "core.mtdna_nucleoids", "yeast.cvt", "yeast.nvj", "core.polysomes", "yeast.rac", "core.mrna_threads", "core.mrnp_ends"]


@registry.register("preplaced", "test.a")
def _a(ctx, params):
    ctx.pre["a_pos"] = np.arange(params["n"] * 3, dtype=np.float32).reshape(-1, 3)
    ctx.obs["a"] = (ctx.pre["a_pos"].astype(np.float64), np.ones(params["n"]))
    ctx.consume(["X"], params["n"])


@registry.register("preplaced", "test.b")
def _b(ctx, params):
    C, R = ctx.obs["a"]                                  # sees the builder before it
    ctx.pre["b_pos"] = C[:1] + 1.0
    ctx.pre["a_pos"] = ctx.pre["a_pos"][1:]              # replaces an array of a: part of b's layer
    ctx.state["seen"] = float(ctx.exp3.expected.sum())


def _ctx(lines):
    exp3 = pd.DataFrame(dict(lib_id=["X", "X", "Y"], ing_id=["X", "X", "Y"], region=["cytosol", "nucleoplasm", "cytosol"],
                             expected=[30.0, 10.0, 5.0]))
    return serial.Context(exp3=exp3, log=lambda *a: lines.append(" ".join(map(str, a))))


def test_order_layers_consume():
    lines = []
    ctx = _ctx(lines)
    with tempfile.TemporaryDirectory() as d:
        rows = serial.run_builders(ctx, [{"name": "a", "builder": "test.a", "params": {"n": 4}}, {"name": "b", "builder": "test.b"}], d)
        assert [r["name"] for r in rows] == ["a", "b"]
        assert rows[0]["keys"] == ["a_pos"] and sorted(rows[1]["keys"]) == ["a_pos", "b_pos"]
        assert rows[1]["counts"]["a_pos"] == 3
        U = {}
        for r in rows:
            with np.load(os.path.join(d, "layers", f"{r['name']}.npz")) as z:
                U.update({k: z[k] for k in z.files})
        assert sorted(U) == sorted(ctx.pre) and all(np.array_equal(U[k], ctx.pre[k]) for k in U)
    # consume: 4 of 40 X taken proportionally (30 -> 27, 10 -> 9), Y untouched; b saw the reduced pool
    assert np.allclose(ctx.exp3.expected.values, [27.0, 9.0, 5.0]) and ctx.state["seen"] == 41.0
    assert lines == ["  X: 4 of 40 used (4 needed)"]
    # capped at the pool
    lines.clear(); assert ctx.consume(["Y"], 9) == 5.0 and ctx.exp3.expected.values[2] == 0.0


def test_amphiphysin_pool_heterodimer():
    """yeast.amphiphysin draws Rvs161-Rvs167 from the heterodimer ingredient when stage 1 assembled it (B5 `common`
    localisation: both cytosol, every Rvs167 in RVS161_RVS167_het) and from free monomer pairs otherwise (legacy);
    a monomer-only pool placed 0 of 586 dimers under B5 (2026-11-19)"""
    from plugins.yeast.preplaced import cortex
    het = "RVS161_RVS167_het"
    for rows, pool, after in (
            # B5: heterodimer + leftover Rvs161 only -> pool = het; consumed from het, Rvs161 untouched
            ([(het, "cytosol", 9240.0), ("RVS161", "cytosol", 2765.0)], 9240.0, [8947.0, 2765.0]),
            # legacy: free monomers only -> pool = min; one of each consumed per dimer
            ([("RVS161", "cytosol", 7386.0), ("RVS167", "PM_peripheral", 6893.0), ("RVS167", "cytosol", 2348.0)], 7386.0,
             [7093.0, 6893.0 * (1 - 293 / 9241), 2348.0 * (1 - 293 / 9241)]),
            # both: het first, the rest from monomer pairs
            ([(het, "cytosol", 100.0), ("RVS161", "cytosol", 500.0), ("RVS167", "cytosol", 400.0)], 500.0, [0.0, 307.0, 207.0])):
        exp3 = pd.DataFrame(dict(lib_id=[r[0] for r in rows], ing_id=[r[0] for r in rows], region=[r[1] for r in rows],
                                 expected=[r[2] for r in rows]))
        ctx = serial.Context(exp3=exp3, log=lambda *a: None)
        assert np.isclose(cortex.rvs_pool(exp3, het), pool), (rows, cortex.rvs_pool(exp3, het))
        cortex.consume_rvs(ctx, het, 293)
        assert np.allclose(ctx.exp3.expected.values, after, atol=1e-6), (rows, ctx.exp3.expected.values)


def test_yeast_registry_and_order():
    import yaml
    from wholecell.preplaced import builders  # noqa: F401
    import plugins.yeast  # noqa: F401
    cfg = yaml.safe_load(open(os.path.join(REPO, "config", "cells", "yeast.yaml")))
    bl = cfg["preplaced"]["serial"]["builders"]
    assert [b["builder"] for b in bl] == YEAST_ORDER
    assert len({b["name"] for b in bl}) == len(bl)
    for b in bl:
        registry.get("preplaced", b["builder"])
    registry.get("preplaced_context", cfg["preplaced"]["serial"]["context"])


def test_legacy_names():
    import contextlib
    import io
    import json
    from wholecell import __main__ as cli
    from wholecell.build import legacy_rel, ref_path
    from wholecell.config import legacy_keys
    sys.path.insert(0, os.path.join(REPO, "scripts"))
    import ci_compare_builds as cc
    raw = {"prepass": {"modules": []}, "pack": {"prepass_index": {"npz": "x"}, "regions": [{"name": "n", "kind": "prepass"}]},
           "seeds": {"legacy": {"prepass.lamina": 1}}}
    err = io.StringIO()
    with contextlib.redirect_stderr(err):
        legacy_keys(raw, "old.yaml")
    assert "preplaced" in raw and "prepass" not in raw and raw["pack"]["preplaced_index"] == {"npz": "x"}
    assert raw["pack"]["regions"][0]["kind"] == "preplaced" and raw["seeds"]["legacy"] == {"prepass.lamina": 1}
    assert err.getvalue().count("\n") == 1 and "stage 'prepass' is now 'preplaced'" in err.getvalue()
    assert legacy_rel("pack/preplaced_pos.npy") == os.path.join("pack", "prepass_pos.npy")
    assert legacy_rel("preplaced/preplaced_manifest.json") == os.path.join("prepass", "prepass_manifest.json")
    assert cli.STAGE_ALIASES["prepass"] == "preplaced"
    with tempfile.TemporaryDirectory() as d:
        old, new = os.path.join(d, "old"), os.path.join(d, "new")
        for b, stage in ((old, "prepass"), (new, "preplaced")):
            os.makedirs(os.path.join(b, stage)); os.makedirs(os.path.join(b, "pack"))
            json.dump({"instances": []}, open(os.path.join(b, stage, f"{stage}_manifest.json"), "w"))
            np.save(os.path.join(b, "pack", f"{stage}_off.npy"), np.arange(3))
            json.dump({stage: {"modules": []}, "f": f"{b}/pack/{stage}_libs.json"}, open(os.path.join(b, "config.resolved.json"), "w"))
        assert ref_path(old, "pack/preplaced_off.npy") == os.path.join(old, "pack", "prepass_off.npy")
        assert ref_path(new, "pack/preplaced_off.npy") == os.path.join(new, "pack", "preplaced_off.npy")
        r = cc.compare(old, new)
        assert r["identical"] and len(r["aliased"]) == 2, r
        np.save(os.path.join(new, "pack", "preplaced_off.npy"), np.arange(4))
        assert not cc.compare(old, new)["identical"]


if __name__ == "__main__":
    for n, f in sorted(globals().items()):
        if n.startswith("test_") and callable(f):
            f(); print("ok", n)
