"""Tests of the B4 pieces of the common pre-placed builders (seconds on 1 CPU):

    python tests/test_preplaced_b4.py        (or pytest tests/test_preplaced_b4.py)

  kernels        core kernel set: label / membrane codes by common name, label masks
  priors         EM point priors -> chains: groups take the largest chains (rank matched, split at max_group),
                 singletons seed the rest, leftovers become monosomes; every prior is used once
  table          transcripts from a per-mRNA table: ribosomes in chains of >= 2 match the polysome budget (2 %)
  subset         frozen_subset exports: instance rows and bead threads selected by a key, consumption per ingredient
  registry       the common builders and the core context are registered
  linker         addendum 2026-11-17: params.linker_nt_mean is a per-cell config value (Jurkat: a number with provenance;
                 yeast: absent = per chain from its transcripts, unchanged)
  walk_option    walk fix 2026-11-18: params.walk of core.polysomes (candidate search for grown ribosomes) defaults to
                 the legacy walk (yeast absent); Jurkat sets it
  frozen_drops   npc-overlap fix 2026-11-17: instances fixed by the geometry (NPCs at pores) drop the overlapping frozen
                 obstacle instances (sphere trees, selected kinds only): out of the obstacle set and the frozen_subset export
  j5_classes     addendum 2026-11-17 J5b rule: per pair class per 1e6 instances of the class's B4 modules, absent class
                 = rate 0 in the reference, classes of a new module reported (pass None), tolerance max(10 %, 1 per 1e6)
  atp_rows       J6 addendum 2026-11-18: Jurkat ATP synthase dimers use the dense rim rows (rim_sampling rows, tm_tol 4.5,
                 literature row spacing); yeast keeps the crista prims (no rows params); rows mode is off by default
The exactness of the yeast builders against the port build: scripts/preplaced_a6_check.py; Jurkat numbers:
scripts/b4_preplaced_check.py measure.
"""
import os
import sys
import tempfile
import types

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


def test_kernels():
    from wholecell.preplaced import kern
    K = kern.core()
    assert K.lab("cytosol") == 1 and K.lab("mito_matrix") == 6
    assert K.lab("nuclear_pore_channel") == K.lab("pore_channel")
    assert K.mask("cytosol", "nucleoplasm") == (1 << 1) | (1 << 3)
    assert K.mem("ER") == 2 and K.nmem >= 8


def test_priors():
    from wholecell.preplaced.builders.polysomes import assign_priors
    rng = np.random.default_rng(0)
    P = np.r_[np.zeros((5, 3)), np.ones((3, 3)) * 100, np.ones((1, 3)) * 300, np.ones((1, 3)) * 400, np.ones((1, 3)) * 500]
    P = P + rng.normal(0, 1, P.shape)
    lab = np.array([0, 0, 0, 0, 0, 1, 1, 1, 2, 3, 4])
    er = np.zeros(11, bool); er[5:8] = True
    sizes = np.array([2, 9, 4, 3])
    order, s2, npri, PP, erf, mono = assign_priors(P, lab, er, sizes, np.random.default_rng(1), max_size=4)
    # groups of 5 (split 4 + 1) and 3: parts 4, 3, 1 take the chains of 9, 4, 3 ribosomes; one singleton seeds the
    # chain of 2; two singletons are left as monosomes
    assert list(order[:2]) == [1, 2] and list(npri[:2]) == [4, 3] and erf[1] and not erf[0]
    assert npri.sum() + len(mono) == 11 and len(PP) == npri.sum()
    assert (s2[order] >= npri).all()


def test_table():
    from wholecell.preplaced.builders.polysomes import transcripts_table
    T = pd.DataFrame(dict(gene=[f"g{i}" for i in range(5000)], length=np.random.default_rng(2).integers(500, 4000, 5000)))
    m = transcripts_table(np.random.default_rng(3), T, 20000.0, 0.8, 40, log=lambda *a: None)
    k = m.ribosomes.values
    assert abs(k[k >= 2].sum() / 20000.0 - 1) < 0.02
    assert {"gene", "length", "cds", "utr5", "ribosomes"} <= set(m.columns)


def test_subset():
    from wholecell.preplaced import context as C
    d = tempfile.mkdtemp()
    f = os.path.join(d, "inst.npz")
    np.savez(f, lib_ids=np.array(["MT_X", "F_A"]), types=np.array([0, 0, 1, 1, 1], np.int32), pos=np.zeros((5, 3), np.float32),
             quat=np.zeros((5, 4), np.float32), kind=np.array([1, 0, 2, 2, 2], np.int8))
    b = os.path.join(d, "beads.npz")
    np.savez(b, pts=np.arange(30, dtype=np.float32).reshape(10, 3), off=np.array([0, 3, 7, 10]), kind=np.array([0, 2, 1]),
             inrib=np.arange(10) % 2, bead_radius=np.float32(1.0))
    ctx = types.SimpleNamespace(cfg={}, bdir=d)
    import wholecell.preplaced.stage as S
    old = S._resolve
    S._resolve = lambda cfg, bdir, v: v
    try:
        j = C._frozen_subset(ctx, dict(npz=f, key="kind", values=[0, 2], lib_to_ingredient={"MT_X": "TUB", "F_A": "ACT"}), os.path.join(d, "o.npz"))
        z = np.load(os.path.join(d, "o.npz"))
        assert j["consumed"] == {"TUB": 1, "ACT": 3} and len(z["types"]) == 4
        j = C._frozen_subset(ctx, dict(npz=b, key="kind", values=[2]), os.path.join(d, "ob.npz"))
        z = np.load(os.path.join(d, "ob.npz"))
        assert j["n_threads"] == 1 and list(z["off"]) == [0, 4] and np.array_equal(z["pts"][:, 0], [9, 12, 15, 18])
        assert list(z["inrib"]) == [1, 0, 1, 0]
    finally:
        S._resolve = old


def test_registry():
    from wholecell import registry
    from wholecell.preplaced import builders, context  # noqa: F401
    for n in ("core.npc_pores", "core.microtubules", "core.microtubules_mtoc", "core.microtubule_pools", "core.polysomes",
              "core.mrna_threads", "core.mrnp_ends", "core.npc_receptors", "core.fg_mesh", "core.atp_synthase_dimers",
              "core.mtdna_nucleoids"):
        registry.get("preplaced", n)
    registry.get("preplaced_context", "core.cellgeom")


def test_linker():
    from wholecell import config as C
    def poly(cell):
        return [b for b in C.load(cell)["preplaced"]["serial"]["builders"] if b["builder"] == "core.polysomes"][0]["params"]
    assert "linker_nt_mean" not in poly("yeast")
    j = poly("jurkat")
    assert isinstance(j.get("linker_nt_mean"), (int, float)) and j["linker_nt_mean"] >= 5.0


def test_walk_option():
    """walk fix 2026-11-18: params.walk (candidate search) is a generic option of core.polysomes whose default is the
    legacy walk (candidates 0: yeast bitwise unchanged, scripts/preplaced_a6_check.py); Jurkat sets it in its config"""
    from wholecell import config as C
    from wholecell.preplaced.builders.polysomes import WALK_SEARCH
    assert WALK_SEARCH["candidates"] == 0 and WALK_SEARCH["priors_first"] is False
    def poly(cell):
        return [b for b in C.load(cell)["preplaced"]["serial"]["builders"] if b["builder"] == "core.polysomes"][0]["params"]
    assert not (poly("yeast").get("walk") or {}).get("candidates")
    w = poly("jurkat").get("walk") or {}
    assert w.get("candidates", 0) > 0 and w.get("rounds", 0) >= 1
    assert w.get("pick", "nearest") in ("nearest", "first", "compact") and set(w) <= set(WALK_SEARCH)


def test_frozen_drops():
    from wholecell.preplaced import context as C, kern
    from wholecell.preplaced.serial import Context
    d = tempfile.mkdtemp()
    LIB = {}
    for n, r in (("NPC", 10.0), ("F_A", 2.0)):
        LIB[n] = os.path.join(d, f"{n}.npz"); np.savez(LIB[n], st_c=np.zeros((1, 3)), st_r=np.array([r]))
    f = os.path.join(d, "cyto.npz")                         # kind 2 at x = 11 (overlaps the NPC sphere), 13 (clear); kind 0 at 11.5
    np.savez(f, lib_ids=np.array(["F_A"]), types=np.zeros(3, np.int32), pos=np.array([[11, 0, 0], [13, 0, 0], [11.5, 0, 0]], np.float32),
             quat=np.tile(np.array([1, 0, 0, 0], np.float32), (3, 1)), kind=np.array([2, 2, 0], np.int8))
    import wholecell.preplaced.stage as S
    old = S._resolve
    S._resolve = lambda cfg, bdir, v: v
    try:
        ctx = Context(cell=None, lp=C.libpack(["NPC"], LIB), exp3=None, log=lambda *a: None, K=kern.core(), cfg={}, bdir=d, LIB=LIB)
        ctx.frozen_specs = [dict(name="cortical_actin", npz=f, kind_key="kind", kinds=[2], obs="actin", shape="bounding")]
        ctx.frozen_paths, ctx.frozen_drops = [f], {}
        C._set_frozen_obs(ctx)
        assert len(ctx.obs["actin"][1]) == 2
        out = C.drop_frozen_overlaps(ctx, [0], np.zeros((1, 3)), np.array([[1.0, 0, 0, 0]]), {"actin"}, why="test")
        assert out == {"cortical_actin": 1} and list(ctx.frozen_drops[f]) == [0]
        assert np.allclose(ctx.obs["actin"][0], [[13, 0, 0]]) and ctx.frozen_layers[0]["dropped"] == 1
        j = C._frozen_subset(ctx, dict(npz=f, key="kind", values=[0, 2]), os.path.join(d, "o.npz"))
        z = np.load(os.path.join(d, "o.npz"))
        assert j["n"] == 2 and j["dropped"] == 1 and j["consumed"] == {"F_A": 2} and np.allclose(z["pos"][:, 0], [13, 11.5])
    finally:
        S._resolve = old


def test_j5_classes():
    sys.path.insert(0, os.path.join(REPO, "scripts", "dev"))
    from b4_preplaced_check import j5_addendum
    b4 = ("polysomes", "microtubules", "npc", "atp_synthase")
    def side(inst, pairs):
        rates = {}
        for c, v in pairs.items():
            ms = sorted(set(m for m in c.split("|") if m in b4))
            den = sum(inst[m] for m in ms)
            rates[c] = dict(pairs=v, den_modules=ms, den=den, per_1e6=v / den * 1e6)
        return dict(instances=inst, pairs=pairs, class_rates_addendum=rates, prepass_prepass_per_1e6=0.0, prepass_prepass_pairs=sum(pairs.values()))
    old = side(dict(microtubules=1000, polysomes=2000, npc=10, chromatin=5000, cytoskeleton_rest=800),
               {"microtubules|microtubules": 100, "polysomes|polysomes": 50, "chromatin|npc": 1})
    new = side(dict(microtubules=1000, polysomes=500, npc=30, chromatin=5000, cytoskeleton_rest=800, atp_synthase=10),
               {"microtubules|microtubules": 105, "polysomes|polysomes": 20, "chromatin|npc": 3, "cytoskeleton_rest|npc": 1,
                "atp_synthase|polysomes": 1})
    r = j5_addendum(new, old, b4)["classes"]
    assert r["microtubules|microtubules"]["pass_"] is True              # +5 %: within 10 %
    assert r["polysomes|polysomes"]["pass_"] is False                   # 40,000 vs 25,000 per 1e6 polysome instances
    assert r["chromatin|npc"]["pass_"] is True                          # per NPC 100,000 = ref (absolute pairs x3, NPCs x3)
    assert r["cytoskeleton_rest|npc"]["pass_"] is False and r["cytoskeleton_rest|npc"]["per_1e6_ref"] == 0.0
    assert r["atp_synthase|polysomes"]["pass_"] is None and r["atp_synthase|polysomes"]["new_module"]


def test_atp_rows():
    import inspect
    from wholecell import config as C
    from wholecell.preplaced.builders import atp_synthase_dimers as AD
    def atp(cell):
        return [b for b in C.load(cell)["preplaced"]["serial"]["builders"] if b["builder"] == "core.atp_synthase_dimers"][0]["params"]
    y = atp("yeast")
    assert y.get("rims", "prims") == "prims" and not {"rim_sampling", "tilt_deg", "seeds_per_cell"} & set(y)
    j = atp("jurkat")
    assert j["rims"] == "lamellae" and j["rim_sampling"] == "rows" and j["tm_tol_nm"] == 4.5
    assert 12.0 <= j["spacing_nm"] <= 20.0 and j["tilt_deg"][0] == 0                    # literature row spacing 12 to > 20 nm
    assert 'params.get("rim_sampling", "random")' in inspect.getsource(AD.place_lamellae)  # default: the B4 random sampling


if __name__ == "__main__":
    for t in (test_kernels, test_priors, test_table, test_subset, test_registry, test_linker, test_walk_option, test_frozen_drops,
              test_j5_classes, test_atp_rows):
        t(); print("ok", t.__name__)
