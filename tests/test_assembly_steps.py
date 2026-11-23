"""Core assembly step list + presets (wholecell/inventory/assembly): preset loading, config errors, a testcell run.

    PYTHONHASHSEED=0 python tests/test_assembly_steps.py        (or pytest tests/test_assembly_steps.py)
"""
import json
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.config import load as load_config  # noqa: E402
from wholecell.inventory import assembly as AS  # noqa: E402
from wholecell.inventory import secms, stage  # noqa: E402


def _raises(fn, *needles):
    try:
        fn()
    except AS.AssemblyConfigError as e:
        msg = str(e)
        for n in needles:
            assert n in msg, (n, msg)
        return msg
    raise AssertionError(f"no AssemblyConfigError (expected {needles})")


def test_presets():
    assert "legacy_v7" in AS.presets() and "legacy_x5" in AS.presets()
    o = AS.options({"preset": "legacy_v7"})
    assert o["steps"] == ["curated_machines", "complex_portal", "secms_fractions", "pdb_homo_oligomers",
                          "afdb_homodimers", "afdb_heterodimers"]
    assert set(o["steps"]) == set(AS.DEFAULT_ORDER)
    # the SEC-MS rule options of legacy_v7 are the secms.py defaults (the v7 rules)
    from wholecell.inventory.assembly import steps as S
    assert S.legacy_options(o) == secms.ASSEMBLY_DEFAULTS
    _raises(lambda: AS.options({"preset": "legacy_v8"}), "inventory.params.assembly_steps.preset", "did you mean 'legacy_v7'")
    _raises(lambda: AS.options({}), "preset is required")


def test_config_errors():
    v7 = {"preset": "legacy_v7"}
    _raises(lambda: AS.options(dict(v7, overrides={"complex_portl": {}})), "overrides.complex_portl", "did you mean 'complex_portal'")
    _raises(lambda: AS.options(dict(v7, overrides={"complex_portal": {"copies": "measured_freepool"}})),
            "complex_portal.copies", "stoichiometry_cap", "did you mean 'measured_free_pool'")
    _raises(lambda: AS.options(dict(v7, overrides={"steps": ["complex_portal", "afdb_homodimer"]})),
            "unknown step 'afdb_homodimer'", "did you mean 'afdb_homodimers'")
    _raises(lambda: AS.options(dict(v7, overrides={"steps": ["secms_fractions", "complex_portal"]})), "come first")
    o = AS.options(v7)
    full = {k: "x" for s in AS.NEEDS.values() for k in s[0]}
    data = {"histone_families": "x", "nucleosome": {}}
    typo = {("sift" if k == "sifts" else k): v for k, v in full.items()}
    _raises(lambda: AS.check_params(dict(v7, inputs=typo, data=data), o, True),
            "inputs.sifts is required by step complex_portal", "did you mean 'sift'")
    _raises(lambda: AS.check_params(dict(v7, inputs=full, data=data), o, False), "secms_fractions needs inventory.adapters.secms")
    _raises(lambda: AS.check_params(dict(v7, inputs=full, data={"histone_familes": "x"}), o, True),
            "data.histone_families is required", "did you mean 'histone_familes'")


def test_jurkat_uses_legacy_v7():
    cfg = load_config("jurkat")
    registry.load_plugins(cfg)
    assert cfg["inventory"]["adapters"]["assembly"] == "assembly_steps"
    p = cfg["inventory"]["params"]["assembly_steps"]
    o = AS.options(p)
    assert o["preset"] == "legacy_v7" and o["complex_portal"]["copies"] == "stoichiometry_cap"
    AS.check_params(p, o, True)
    outs = stage.steps(cfg)["assembly"]["outputs"]
    for t in ("complexes_v0", "inventory_v0_complexed", "sec_assembly", "sec_calibrants_peak", "homo_oligomers_v0",
              "afdb_homodimers_v0", "afdb_heterodimers_v0"):
        assert f"inventory/assembly/{t}.tsv" in outs, t


def test_apex_check_exclude():
    """data.apex_check_exclude (Step 3, 2026-11-19): Jurkat lists HSPD1 with its share read from Fossati 2021; the v7 port
    lists nothing; rows need uniprot + reason; a share_source_id is a SEC-MS run with a row for the protein"""
    from wholecell.inventory.assembly import steps as S
    cfg = load_config("jurkat")
    p = cfg["inventory"]["params"]["assembly_steps"]
    o = AS.options(p)
    AS.check_params(p, o, True)                                  # a declared optional data entry of pdb_homo_oligomers
    ex = p["data"]["apex_check_exclude"]
    assert [(r["uniprot"], r["gene"], r["share_source_id"]) for r in ex] == [("P10809", "HSPD1", "fossati2021_jurkat")]
    assert all(r["reason"].strip() for r in ex)
    assert "fossati2021_jurkat" in {d["source_id"] for d in cfg["inventory"]["params"]["secms_toolkit"]["datasets"]}
    port = load_config("jurkat_v7_port")["inventory"]["params"]["assembly_steps"]
    assert port["data"]["apex_check_exclude"] == []              # v7 reference unchanged

    class View:
        def __init__(self, src, rows):
            self.source_id, self.condition, self.replicate = src, "c", "1"
            self.thresholds = pd.DataFrame(rows).T

        def has(self, a): return a in self.thresholds.index

        def a(self, a, k): return float(self.thresholds.at[a, k])

    views = {"fos": View("fos", {"P1": {2: 0.9, 3: 0.9, 4: 0.84}})}
    orig = secms.SecView
    secms.SecView = lambda ctx, sel: views[sel["source_id"]]
    try:
        class St:
            ctx, view = None, View("final", {"P1": {2: 0.9, 3: 0.6, 4: 0.47}})
            data = {"apex_check_exclude": [{"uniprot": "P1", "gene": "G1", "share_source_id": "fos", "reason": "artefact"},
                                           {"uniprot": "P2", "gene": "G2", "reason": "artefact, rules' run share"}]}
        e = S.apex_check_exclude(St())
        assert e["P1"]["view"] is views["fos"] and e["P2"]["view"] is None
        lo = S.legacy_options(o)
        assert secms.homo_share(e["P1"]["view"], "P1", 7, lo) == (0.84, 4)          # heptamer: a_4 of the named run
        assert secms.homo_share(St.view, "P1", 7, lo) == (0.47, 4)
        St.data = {"apex_check_exclude": [{"uniprot": "P1", "share_source_id": "fos"}]}
        try:
            S.apex_check_exclude(St())
            raise AssertionError("row without reason accepted")
        except ValueError as err:
            assert "reason" in str(err)
        St.data = {"apex_check_exclude": [{"uniprot": "P9", "share_source_id": "fos", "reason": "x"}]}
        try:
            S.apex_check_exclude(St())
            raise AssertionError("protein without a row in the share source accepted")
        except ValueError as err:
            assert "no SEC-MS row" in str(err)
    finally:
        secms.SecView = orig


def test_yeast_legacy_x5():
    """legacy_x5: pools ingredients, yeast plugin steps registered, every input / data entry declared, curated TSVs
    readable with a reason column, config errors name the plugin step and the hook"""
    cfg = load_config("yeast")
    registry.load_plugins(cfg)
    p = cfg["inventory"]["params"]["assembly_steps"]
    assert p["inputs"]["protein_table"] == "build:inventory/yeast/protein_table.parquet"   # A9: plugin step yeast_protein_table
    q = dict(p, inputs={k: v for k, v in p["inputs"].items() if k != "protein_table"})
    _raises(lambda: AS.check_params(q, AS.options(q), False), "inputs.protein_table is required by step curated_machines")
    o = AS.options(p)
    assert o.get("pools") == "ingredients" and o["steps"][0] == "curated_machines"
    assert o["steps"][1] == "yeast_barrel_monomers" and "yeast_preribosomes" in o["steps"] and o["steps"][-2:] == ["afdb_heterodimers", "pdb_homo_oligomers"]
    AS.check_params(p, o, False)
    from wholecell.inventory.assembly.steps import read_curated
    for k, v in p["data"].items():
        assert "reason" in read_curated(v).columns, k
    outs = AS.outputs({"inventory": {"params": {"assembly_steps": p}}})
    for f in ("ingredients.csv", "oligo_map.json", "complexes_tried.csv", "complexes_added.csv", "heterodimers.csv",
              "oligomers.csv", "f_actin.json"):
        assert f"assembly/{f}" in outs, f
    x5 = {"preset": "legacy_x5"}
    _raises(lambda: AS.options(dict(x5, overrides={"steps": ["curated_machines", "yeast_preribosome"]})),
            "unknown step 'yeast_preribosome'", "did you mean 'yeast_preribosomes'")
    _raises(lambda: AS.options(dict(x5, overrides={"steps": ["complex_portal"]})), "starts with curated_machines")
    _raises(lambda: AS.options(dict(x5, overrides={"afdb_heterodimers": {"builder": "yeast_structure"}})),
            "not a registered assembly_hook", "did you mean 'yeast_structures'")
    _raises(lambda: AS.options(dict(x5, overrides={"complex_portal": {"overlap": "skip_if_use"}})),
            "complex_portal.overlap", "did you mean 'skip_if_used'")
    _raises(lambda: AS.check_params(dict(p, inputs={k: v for k, v in p["inputs"].items() if k != "rcsb_cache"}), o, False),
            "inputs.rcsb_cache is required by step complex_portal")


def test_yeast_x5_secms():
    """x5_secms (Step 3) = legacy_x5 + the legacy_v7 SEC-fed rules with the same cut-offs; the share rule splits an
    ingredient into '<id>__homo' (oligomer, parent's comp_weights) + the monomer remainder"""
    import numpy as np
    from wholecell.inventory.assembly import ingredient_steps as IS
    cfg = load_config("yeast")
    registry.load_plugins(cfg)
    o = AS.options(cfg["inventory"]["params"]["assembly_steps"])
    x5 = AS.options({"preset": "legacy_x5"})
    assert o["preset"] == "x5_secms" and o["secms"]["source_id"] == "yeast_final"
    v7 = AS.options({"preset": "legacy_v7"})
    sh = {k: v for k, v in o["pdb_homo_oligomers"]["share"].items() if k not in ("apex_k", "apex_check_side")}
    assert sh == v7["pdb_homo_oligomers"]["share"] and o["pdb_homo_oligomers"]["af_share"] == v7["afdb_homodimers"]["share"]
    assert o["pdb_homo_oligomers"]["share"]["apex_check_side"] == "below"                   # addendum A1
    fb = dict(o["pdb_homo_oligomers"]["afdb_fallback"])
    assert fb.pop("ipsae_only_min") == v7["afdb_homodimers"]["min_ipsae"]                    # addendum A2
    assert fb == x5["pdb_homo_oligomers"]["afdb_fallback"]
    ev = o["afdb_heterodimers"]["evidence"]
    assert ev["sec_r_min"] == v7["afdb_heterodimers"]["evidence"]["sec_r_min"] and ev["string_min"] == v7["afdb_heterodimers"]["evidence"]["string_min"]
    strip = lambda d: {k: ({kk: vv for kk, vv in v.items() if kk not in ("copies", "share", "af_share", "evidence", "afdb_fallback")}
                           if isinstance(v, dict) else v) for k, v in d.items() if k not in ("secms", "preset")}
    assert strip(o) == strip(x5)                                                      # everything else = legacy_x5
    assert AS.options(load_config("yeast_x5_core")["inventory"]["params"]["assembly_steps"]) == x5

    class View:
        thresholds = pd.DataFrame({2: [0.8, 0.6, 0.9], 3: [0.7, 0.5, 0.8], 4: [0.4, 0.1, 0.7]}, index=["P1", "P2", "P4"])

        def has(self, a): return a in self.thresholds.index

        def a(self, a, k): return float(self.thresholds.at[a, k]) if self.has(a) else np.nan

        def legacy_sec_assembly(self, k):
            return pd.DataFrame(dict(uniprot=["P1", "P2", "P4"], apex_app_kda=[200.0, 100.0, 1000.0], mono_kda=[50.0, 50.0, 50.0]))

    class St:
        view = View()
        ing = pd.DataFrame(dict(ing_id=["A", "B", "C", "D"], kind="protein", uniprot=["P1", "P2", "P3", "P4"],
                                label=["A", "B", "C", "D"], copies_per_cell=[1000.0] * 4, comp_weights=['{"cytosol": 1.0}'] * 4))
    for side, kept in (("below", True), ("both", False)):     # D: dimer, apex 1000 / (2 x 50) = 10 > 2.5, a_4 0.7 > 0.2
        st = St()
        st.ing = St.ing.copy()
        oo = dict(o["pdb_homo_oligomers"], share=dict(o["pdb_homo_oligomers"]["share"], apex_check_side=side))
        g = IS._homo_share_rule(st, oo)
        od, rd = {"D": ("D_oligo", 2)}, dict(decision="PDB")
        g(next(st.ing[st.ing.ing_id == "D"].itertuples()), od, rd, from_af=False)
        assert ("D" in od) == kept and (rd.get("share") == 0.9 if kept else rd["share"] == 0.0)
    st = St()
    st.ing = St.ing[St.ing.ing_id != "D"].reset_index(drop=True)
    f = IS._homo_share_rule(st, o["pdb_homo_oligomers"])
    oligo = {"A": ("A_oligo", 4), "B": ("B_oligo", 2), "C": ("C_oligo", 2)}
    R = {r.ing_id: r for r in st.ing.itertuples()}
    rows = {}
    for i, af in (("A", False), ("B", True), ("C", True)):
        rows[i] = dict(decision="PDB" if not af else "AF dimer")
        f(R[i], oligo, rows[i], from_af=af)
    f.apply(oligo)
    c = st.ing.set_index("ing_id").copies_per_cell
    assert rows["A"]["share"] == 0.7 and rows["A"]["sec_k"] == 3          # tetramer: a_3; apex 200 / (4 x 50) = 1 ok
    assert rows["B"]["share"] == 0.6 and rows["C"]["share"] == 0.3        # AF dimer: a_2; no SEC row: default 0.3
    assert np.isclose(c["A__homo"], 700) and np.isclose(c["A"], 300) and np.isclose(c["C__homo"], 300)
    assert oligo["A__homo"] == ("A_oligo", 4) and "A" not in oligo
    assert (st.ing.comp_weights == '{"cytosol": 1.0}').all()


def test_ingredient_helpers():
    """ported yeast helpers on small cases (exact values of the yeast code)"""
    from wholecell.inventory.assembly import ingredient_steps as IS
    pool = {"A": 100.0, "B": 50.0, "COB": 1.0}
    assert IS._complex_copies(pool, {"A": 1, "B": 1, "COB": 1}, {"COB"}) == 75.0     # mtDNA subunit left out
    assert IS._complex_copies(pool, {"A": 1, "B": 1, "COB": 1}, set()) == 2.0        # 2 x limiting
    assert IS._k3("ABCDE", "ABCDF") == 2 / 3
    sc = dict(min_shared=0.3, facing={"ER"}, inner="mito_IM", inner_side=["mito_matrix", "mito_IM"])
    assert IS._compatible({"cytosol": 1.0}, {"cytosol": 0.5, "nucleoplasm": 0.5}, False, False, **sc) == {"cytosol": 1.0}
    assert IS._compatible({"ER": 1.0}, {"cytosol": 1.0}, True, False, **sc) == {"ER": 1.0}
    assert IS._compatible({"ER": 1.0}, {"nucleoplasm": 1.0}, True, False, **sc) is None


def test_testcell_step_list():
    """the step list on the testcell (no Complex Portal / PDB inputs): steps [secms_fractions] only"""
    ov = {"inventory": {"adapters": {"assembly": "assembly_steps"},
                        "params": {"assembly_steps": {"preset": "legacy_v7", "overrides": {"steps": ["secms_fractions"]}}}}}
    cfg = load_config("testcell", overrides=ov)
    registry.load_plugins(cfg)
    plan = stage.steps(cfg)
    assert plan["assembly"]["outputs"] == [f"inventory/assembly/{x}" for x in AS.ALWAYS + AS.OUTS["secms_fractions"]]
    with tempfile.TemporaryDirectory() as d:
        for s in ("v0", "secms", "assembly"):
            stage.run(cfg, d, step=s)
        a = os.path.join(d, "inventory", "assembly")
        I = pd.read_csv(os.path.join(a, "inventory_v0_complexed.tsv"), sep="\t")
        assert len(I) > 0 and (I.copies_in_complexes == 0).all()
        assert len(pd.read_csv(os.path.join(a, "complexes_v0.tsv"), sep="\t")) == 0
        S = pd.read_csv(os.path.join(a, "sec_assembly.tsv"), sep="\t")
        assert len(S) > 0 and S.assembled_frac.between(0, 1).all()
        j = json.load(open(os.path.join(a, "assembly_steps.json")))
        assert j["preset"] == "legacy_v7" and j["steps"] == ["secms_fractions"]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "ok")
