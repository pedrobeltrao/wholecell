"""SEC-MS toolkit + assembly adapter tests (seconds; synthetic testcell SEC-MS with known truth).

    PYTHONHASHSEED=0 python tests/test_secms.py        (or pytest tests/test_secms.py)
"""
import json
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from toolkits.secms import pipeline  # noqa: E402
from wholecell import registry  # noqa: E402
from wholecell.config import load as load_config  # noqa: E402
from wholecell.inventory import schema, secms, stage  # noqa: E402

TC = os.path.join(REPO, "tests", "data", "testcell")
_RES = {}


def _toolkit():
    if "r" not in _RES:
        cfg = load_config("testcell")
        p = secms._frozen_paths(dict(cfg["inventory"]["params"]["secms_toolkit"]), type("C", (), {"cfg": cfg, "frozen": cfg.frozen})())
        _RES["r"] = pipeline.run(p)
    return _RES["r"]


def test_calibration_recovers_column_law():
    d = _toolkit()["_diagnostics"]["datasets"][0]
    for run, r in d["runs"].items():
        c = r["calibration"]
        assert c["used"] == "anchors", run                       # species anchors are the calibration (v7 approach)
        assert abs(c["slope"] + 0.055) < 0.003 and abs(c["icpt"] - 3.6) < 0.06, (c["icpt"], c["slope"])
        assert c["spearman_rho"] < -0.95 and c["resid_fold"] < 1.2
        assert c["checks"]["standards"]["own_fit"]["resid_fold"] < 1.05       # standards check fitted and reported
        assert "complex_portal" in c["checks"]                                # diagnostic only


def test_assembled_fraction_vs_truth():
    P = _toolkit()["secms_protein"]
    T = pd.read_csv(os.path.join(TC, "secms_truth.tsv"), sep="\t", comment="#")
    m = P[(P.threshold_k == 2) & (P.replicate.astype(str) == "mean")].merge(T, on="uniprot")
    assert len(m) == len(T)
    well = m[(m.min_ratio >= 2.5) | m.min_ratio.isna()]                     # states well above 2x monomer, or monomers
    err = (well.assembled_frac - well.assembled_true).abs()
    assert err.max() < 0.1, well.loc[err.idxmax()].to_dict()
    assert (m.assembled_frac_lo <= m.assembled_frac + 1e-12).all() and (m.assembled_frac_hi >= m.assembled_frac - 1e-12).all()
    assert m.assembled_frac_sd.notna().all() and (m.n_rep == 2).all()       # replicate SD on the mean rows
    assert set(P.columns[:len(schema.SECMS_PROTEIN)]) == set(schema.SECMS_PROTEIN)
    assert not P.uniprot.str.contains("/|DECOY").any()                     # shared / decoy peptides removed
    g = P[(P.threshold_k == 2) & (P.gene == "PGK1")]
    assert (g.apex_fraction > 10).all()                                     # q-value 0.05 void rows were filtered


def test_complexes_and_pairs():
    r = _toolkit()
    C = r["secms_complex"]
    c = C[(C.complex_id == "CPX-6030") & (C.replicate.astype(str) == "1")].iloc[0]
    assert bool(c.measured) and abs(c.assembled_frac - 0.75) < 0.1 and c.n_subunits_detected == 8
    assert (C.complex_source == "humap").any()                              # hu.MAP membership parsed (human)
    S = r["secms_complex_subunit"]
    assert list(S.columns[:len(schema.SECMS_COMPLEX_SUBUNIT)]) == schema.SECMS_COMPLEX_SUBUNIT
    X = r["secms_pairs"].set_index(["uniprot_a", "uniprot_b", "replicate"])
    acc = pd.read_csv(os.path.join(TC, "secms_annotation.tsv"), sep="\t", comment="#").set_index("gene").uniprot
    assert X.loc[(acc["TCP1"], acc["CCT2"], "mean"), "r"] > 0.95
    assert X.loc[(acc["GAPDH"], acc["PGK1"], "mean"), "r"] < 0.5
    assert X.loc[(acc["TCP1"], acc["CCT2"], "mean"), "n_rep"] == 2


def test_species_switch_and_presets():
    for bad in ("mouse", None):
        try:
            pipeline.load_species(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"species {bad} accepted")
    y = pipeline.load_species("yeast")
    assert y["taxon"] == 559292 and not y["humap_allowed"]
    A = pd.read_csv(os.path.join(REPO, "toolkits", "secms", "data", y["anchors"]), sep="\t", comment="#")
    assert len(A) >= 20 and A.oligomer_n.notna().sum() + A.native_kda.notna().sum() == len(A)
    d, l = pipeline.load_preset("default"), pipeline.load_preset("legacy_v7")
    assert d["calibration"]["primary"] == l["calibration"]["primary"] == "anchors"
    assert d["calibration"]["fallback"] is None
    assert l["quant"]["min_peptides"] == 1 and d["quant"]["min_peptides"] == 2


def _run_testcell(d, overrides=None):
    cfg = load_config("testcell", overrides=overrides)
    registry.load_plugins(cfg)
    stage.run(cfg, d)
    return cfg, os.path.join(d, "inventory")


def test_measured_free_pool_option_testcell():
    with tempfile.TemporaryDirectory() as d0, tempfile.TemporaryDirectory() as d1:
        _, a = _run_testcell(d0)
        _, b = _run_testcell(d1, {"inventory": {"params": {"complex_table": {"assembly_options": {"complex_copies": "measured_free_pool"}}}}})
        Ca = pd.read_csv(os.path.join(a, "assembly", "complexes.tsv"), sep="\t").set_index("complex_id")
        Cb = pd.read_csv(os.path.join(b, "assembly", "complexes.tsv"), sep="\t").set_index("complex_id")
        v0 = pd.read_csv(os.path.join(b, "v0", "inventory_v0.tsv"), sep="\t", index_col="gene")
        S = pd.read_csv(os.path.join(b, "secms", "secms_complex_subunit.tsv"), sep="\t")
        S = S[(S.complex_id == "CPX-6030") & (S.replicate.astype(str) == "mean")] if (S.replicate.astype(str) == "mean").any() \
            else S[(S.complex_id == "CPX-6030") & (S.replicate.astype(str) == "1")]
        exp = np.median([v0.copies[g] * f / 2 for g, f in zip(S.gene, S.in_complex_frac)])
        assert Cb.copies_built["CPX-6030"] != Ca.copies_built["CPX-6030"]
        assert abs(Cb.copies_built["CPX-6030"] - exp) <= 1, (Cb.copies_built["CPX-6030"], exp)
        assert Cb.copies_built["CPX-6151"] == Ca.copies_built["CPX-6151"]   # TM complex keeps the table copies
        assert Cb.copies_built["GAPDH_x4"] == Ca.copies_built["GAPDH_x4"]   # not a measured complex: table copies
        for name, cols in secms.canonical_tables().items():
            T = pd.read_csv(os.path.join(a, name), sep="\t")
            assert list(T.columns) == cols and len(T) > 0, name


def test_jurkat_secms_plan_and_freepool_off():
    """S11: the new option is OFF for Jurkat; the secms step feeds the assembly step. Step 3 (2026-11-19): jurkat.yaml
    runs the final combined SEC-MS (preset default, 3 datasets, label jurkat_final read by the rules); the v7 port keeps
    the v7 run (preset legacy_v7, Fossati only, no combination)"""
    from wholecell.inventory import assembly as AS
    for cell, preset, label, frozen in (("jurkat", "default", "jurkat_final",
                                         {"inv_sec_fossati_features", "inv_sec_lab2026_jurkat_nt", "inv_sec_lab2026_jurkat_tcract24h"}),
                                        ("jurkat_v7_port", "legacy_v7", "inv_sec_fossati_features", {"inv_sec_fossati_features"})):
        cfg = load_config(cell)
        registry.load_plugins(cfg)
        p = cfg["inventory"]["params"]["assembly_steps"]
        o = AS.options(p)
        assert o["complex_portal"]["copies"] == "stoichiometry_cap", cell
        assert o["secms"]["source_id"] == label, cell
        tk = cfg["inventory"]["params"]["secms_toolkit"]
        assert tk["preset"] == preset, cell
        assert (tk.get("combine") or {}).get("label") == (label if cell == "jurkat" else None), cell
        plan = stage.steps(cfg)
        ks = list(plan)
        assert ks.index("secms") < ks.index("assembly")
        assert "inventory/secms/secms_protein.tsv" in plan["assembly"]["inputs"]
        assert set(plan["secms"]["frozen"]) >= frozen and not ({"inv_sec_lab2026_jurkat_nt"} & set(plan["secms"]["frozen"]) - frozen), cell


def test_yeast_secms_tables_only():
    """Step 3: yeast.yaml runs the SEC-MS step (final combined estimate, label yeast_final); the legacy_x5 assembly has
    no SEC-MS reader; the A9 reference yeast_x5_core has no secms step"""
    cfg = load_config("yeast")
    registry.load_plugins(cfg)
    assert cfg["inventory"]["adapters"]["secms"] == "secms_toolkit"
    tk = cfg["inventory"]["params"]["secms_toolkit"]
    assert tk["species"] == "yeast" and tk["preset"] == "default" and tk["combine"]["label"] == "yeast_final"
    assert {d["source_id"] for d in tk["combine"]["datasets"]} == {d["source_id"] for d in tk["datasets"]}
    from wholecell.inventory import assembly as AS
    o = AS.options(cfg["inventory"]["params"]["assembly_steps"])
    assert o["pools"] == "ingredients" and "secms_fractions" not in o["steps"]
    assert "secms" in stage.steps(cfg)
    ref = load_config("yeast_x5_core")
    registry.load_plugins(ref)
    assert not ref["inventory"]["adapters"].get("secms") and "secms" not in stage.steps(ref)


def test_inventory_step_combine():
    """inventory step `secms` with secms_toolkit.combine: combined rows (label, condition = replicate = 'combined')
    appended to every table, selectable by SecView"""
    import tempfile
    cfg = load_config("testcell", overrides={"inventory": {"params": {"secms_toolkit": {"combine": {
        "label": "tc_final", "datasets": [{"source_id": "testcell_sec_synthetic", "weight": 1.0, "role": "primary"}]}}}}})
    registry.load_plugins(cfg)
    with tempfile.TemporaryDirectory() as d:
        stage.run(cfg, d, step="secms")
        P = pd.read_csv(os.path.join(d, "inventory", "secms", "secms_protein.tsv"), sep="\t", dtype={"replicate": str})
        C = P[P.source_id == "tc_final"]
        assert len(C) and set(C.condition) == {"combined"} and set(C.replicate) == {"combined"}
        assert set(C.uniprot) == set(P[P.source_id == "testcell_sec_synthetic"].uniprot)
        diag = json.load(open(os.path.join(d, "inventory", "secms", "diagnostics.json")))
        assert diag["combine"]["label"] == "tc_final"
        from wholecell.inventory.context import StepContext
        v = secms.SecView(StepContext(cfg, d, "assembly"), {"source_id": "tc_final"})
        assert v.condition == "combined" and v.replicate == "combined" and len(v.thresholds)
    bad = load_config("testcell", overrides={"inventory": {"params": {"secms_toolkit": {"combine": {
        "label": "x", "datasets": [{"source_id": "nope"}]}}}}})
    with tempfile.TemporaryDirectory() as d:
        try:
            stage.run(bad, d, step="secms")
            raise AssertionError("unknown combine dataset accepted")
        except ValueError as e:
            assert "not in secms_toolkit.datasets" in str(e)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "ok")
