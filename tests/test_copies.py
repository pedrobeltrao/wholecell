"""Copy-number step tests (seconds; synthetic tables only, no Jurkat or yeast inputs needed).

    PYTHONHASHSEED=0 python tests/test_copies.py        (or pytest tests/test_copies.py)

  test_jurkat_recipe_equals_legacy_adapter   criterion C6: the Jurkat recipe of config/cells/jurkat.yaml (inventory.copies,
                                             paths swapped for synthetic MaxQuant / UniProt / RNA tables) gives the same
                                             table (md5 of the TSV) as the ported legacy adapter jurkat.geiger_ibaq_v0
  test_testcell_chain_through_copies         the testcell inventory with the proteome routed through the copies step
                                             (absolute source, total = source sum) gives the same pack_table.tsv
  test_options                               typed sources (lfq, harmonised secondary source), priority provenance, absent
                                             categories, size-slope scaling, family-target adjustment, anchors (check and
                                             total protein from anchors), ledger rule ids, future-source slots
  test_configs                               Jurkat config = legacy settings (anchors diagnostic only); yeast config = PaxDb,
                                             NOT_IN_THIS_CELL, G1 scaling, ribosome adjustment, 3.5 pg fixed
"""
import hashlib
import json
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.config import load as load_config, load_raw  # noqa: E402
from wholecell.inventory import copies, stage  # noqa: E402
from wholecell.inventory.context import StepContext  # noqa: E402


def _md5_frame(df):
    return hashlib.md5(df.to_csv(sep="\t").encode()).hexdigest()


# ------------------------------------------------------------------------------------------- synthetic Jurkat-like data
def _jurkat_like(d, seed=7):
    rng = np.random.default_rng(seed)
    genes = [f"G{i:03d}" for i in range(60)] + ["H2AC1", "H2AC2", "H4C1", "TUBB", "PTEN", "IFNA1", "IFNB1", "MTAP"]
    acc = {g: f"P{i:05d}" for i, g in enumerate(genes)}
    up = pd.DataFrame({"Entry": [acc[g] for g in genes], "Reviewed": "reviewed",
                       "Gene Names (primary)": [g if g != "G005" else "G005; G005B" for g in genes],
                       "Mass": rng.integers(8000, 250000, len(genes)).astype(float)})
    up.loc[len(up)] = ["Q99999", "unreviewed", "G999", 30000.0]
    up.to_csv(os.path.join(d, "uniprot.tsv"), sep="\t", index=False)
    rows = []
    for g in genes[:40] + ["H4C1", "PTEN", "IFNA1"]:
        ids = acc[g] + (";" + acc[g] + "-2" if g == "G003" else "")
        rows.append([ids] + list(np.where(rng.random(3) < 0.2, np.nan, 10 ** rng.normal(4, 1, 3))) + [1.0])
    rows.append([";".join([acc["H2AC1"], acc["H2AC2"], acc["TUBB"]]), 5e5, 6e5, np.nan, 1.0])   # family split -> H2AC1/2
    rows.append([acc["G041"], np.nan, np.nan, np.nan, 3.0])                                          # not in Jurkat columns
    rows.append(["X00001", 100.0, 100.0, 100.0, 1.0])                                                # unmapped accession
    mq = pd.DataFrame(rows, columns=["Majority.protein.IDs", "iBAQ.Jurkat_1", "iBAQ.Jurkat_2", "iBAQ.Jurkat_3", "iBAQ.HeLa_1"])
    mq.to_csv(os.path.join(d, "mq.tsv"), sep="\t", index=False)
    rg = genes[:30] + genes[40:58] + ["TUBB", "MTAP", "IFNB1"]
    calls = (["expressed"] * 30 + ["expressed"] * 6 + ["low"] * 6 + ["off"] * 6 + ["expressed", "low", "off"])
    rna = pd.DataFrame({"gene": rg, "tpm_median": np.round(10 ** rng.normal(0.5, 1, len(rg)), 6), "call": calls})
    rna.loc[rna.call == "off", "tpm_median"] = 0.0
    os.makedirs(os.path.join(d, "inventory", "rna"), exist_ok=True)
    rna.to_csv(os.path.join(d, "inventory", "rna", "rna.tsv"), sep="\t", index=False)
    return dict(uniprot=os.path.join(d, "uniprot.tsv"), mq=os.path.join(d, "mq.tsv"))


def test_jurkat_recipe_equals_legacy_adapter():
    with tempfile.TemporaryDirectory() as d:
        f = _jurkat_like(d)
        C = json.loads(json.dumps(load_raw("jurkat")["inventory"]["copies"]))
        C["annotation"]["path"] = "inv_t_uniprot"
        C["sources"][0]["path"] = "inv_t_mq"
        C["anchors"] = []
        C["adjustments"] = []                  # the legacy v0 recipe (= jurkat_v7_port); v8.copy_fixes: test_jurkat_copy_fixes
        ov = {"plugins": ["jurkat"],
              "frozen_inputs": {"inv_t_uniprot": {"path": f["uniprot"]}, "inv_t_mq": {"path": f["mq"]}},
              "inventory": {"total_protein_pg": 0.25, "copies": C,
                            "params": {"jurkat.geiger_ibaq_v0": {"proteome": "inv_t_mq", "uniprot": "inv_t_uniprot",
                                                                 "sample_match": "Jurkat", "source_label": "geiger",
                                                                 "floor_percentile": 5, "low_factor": 0.25}}}}
        cfg = load_config("testcell", overrides=ov)
        registry.load_plugins(cfg)
        ctx = StepContext(cfg, d, "copies")
        T, S = copies.step(ctx)
        new = copies.copies_adapter(StepContext(cfg, d, "v0"), 0.25)
        legacy = registry.get("adapter.proteome", "jurkat.geiger_ibaq_v0")(StepContext(cfg, d, "v0"), 0.25)
        assert list(new.columns) == list(legacy.columns)
        assert _md5_frame(new) == _md5_frame(legacy), new.compare(legacy)
        assert set(T.source) >= {"geiger", "imputed_rna", "imputed_rna_low", "zero_off", "zero_override"}
        assert T.at["PTEN", "note"].startswith("absent protein") and T.at["IFNA1", "note"] == ""   # legacy: no note for pattern hits
        assert T.at["H2AC1", "ppb_geiger"] == T.at["H2AC2", "ppb_geiger"] and np.isnan(T.at["TUBB", "ppb_geiger"])
        L = pd.read_csv(os.path.join(d, "inventory", "copies", "copy_ledger.tsv"), sep="\t", keep_default_na=False)
        assert set(L[L.key == "PTEN"].rule_id) == {"source.inv_geiger2012_pg", "absent.jurkat_zero_genes"}
        assert (T.set_by != "").all() and T.loc["PTEN", "rule_ids"] == "source.inv_geiger2012_pg;absent.jurkat_zero_genes"
        assert abs(S["protein_pg"] - 0.25) < 1e-3


def test_testcell_chain_through_copies():
    with tempfile.TemporaryDirectory() as d0, tempfile.TemporaryDirectory() as d1:
        cfg0 = load_config("testcell")
        registry.load_plugins(cfg0)
        stage.run(cfg0, d0)
        ov = {"inventory": {"adapters": {"proteome": "copies"}, "copies": {
            "method_id": "testcell_absolute",
            "annotation": {"reader": "annotation_table", "path": "inv_testcell_proteome"},
            "universe": {"require": ["sources"], "order": "sorted"},
            "sources": [{"id": "inv_testcell_proteome", "type": "absolute_copies", "unit": "copies per cell", "reader": "table",
                         "path": "inv_testcell_proteome", "id_column": "gene", "value_columns": ["copies"], "aggregate": "single",
                         "column": "copies_source", "label": "inv_testcell_proteome"}],
            "imputation": {"method": "none"}, "total_protein": {"mode": "source_sum"}, "round_copies": True,
            "v0_columns": ["uniprot", "mass_da", "source", "copies"]}}}
        cfg1 = load_config("testcell", overrides=ov)
        plan = stage.steps(cfg1)
        assert list(plan)[:2] == ["copies", "v0"] and "inventory/copies/copies.tsv" in plan["v0"]["inputs"]
        stage.run(cfg1, d1)
        a = open(os.path.join(d0, "inventory", "pack_table.tsv"), "rb").read()
        b = open(os.path.join(d1, "inventory", "pack_table.tsv"), "rb").read()
        assert a == b
        for o in plan["copies"]["outputs"]:
            assert os.path.exists(os.path.join(d1, o)), o


# --------------------------------------------------------------------------------------------- synthetic option tests
def _yeastlike(d, seed=3):
    rng = np.random.default_rng(seed)
    orfs = [f"Y{i:04d}W" for i in range(80)]
    names = [None if i % 9 == 0 else f"GENE{i}" for i in range(80)]
    names[1:9] = ["HTA1", "HTB1", "HHT1", "HHF1", "RPL1A", "RPL1B", "RPS2", "ASC1"]
    names[10:13] = ["NUP1", "NUP2", "SEC13"]
    names[20:22] = ["SPO11", "GAL1"]
    up = pd.DataFrame({"Entry": [f"Q{i:05d}" for i in range(80)], "Reviewed": "reviewed",
                       "Gene Names (ordered locus)": orfs, "Gene Names (primary)": names,
                       "Length": 300, "Mass": rng.integers(10000, 120000, 80).astype(float)})
    up.to_csv(os.path.join(d, "up.tsv"), sep="\t", index=False)
    ppm = 10 ** rng.normal(1.5, 0.8, 80)
    ppm[1:5] = 2000.0                       # histones: equal ppm
    ppm[[0, 30]] = np.nan                   # not in PaxDb
    with open(os.path.join(d, "paxdb.txt"), "w") as fh:
        fh.write("#name: synthetic\n")
        for i in range(80):
            if not np.isnan(ppm[i]):
                fh.write(f"{names[i] or orfs[i]}\t9999.{orfs[i]}\t{ppm[i]}\n")
        fh.write(f"X\t9999.{orfs[5]}\t1.0\n")       # duplicated id: summed
    lfq = pd.DataFrame({"gene": [names[i] or orfs[i] for i in (0, 30, 31)], "lfq_1": [5e6, 8e6, 1e6], "lfq_2": [6e6, np.nan, 2e6]})
    lfq.to_csv(os.path.join(d, "lfq.tsv"), sep="\t", index=False)
    pd.DataFrame({"orf": orfs[:40], "slope": np.linspace(-1, 0.5, 40)}).to_csv(os.path.join(d, "slopes.tsv"), sep="\t", index=False)
    return d


def _opt_cfg(d, **over):
    C = {"method_id": "synthetic_yeastlike", "avogadro": 6.02214e23,
         "annotation": {"reader": "uniprot_ordered_locus", "path": os.path.join(d, "up.tsv"), "carry": ["orf", "gene_name"]},
         "universe": {"order": "annotation"}, "value_column": "ppm",
         "sources": [{"id": "pax", "type": "ppm", "reader": "paxdb", "path": os.path.join(d, "paxdb.txt"), "strip_prefix": "9999.",
                      "map_on": "orf", "column": "ppm_pax", "label": "paxdb"},
                     {"id": "lfq", "type": "lfq", "reader": "table", "path": os.path.join(d, "lfq.tsv"), "id_column": "gene",
                      "value_columns": ["lfq_1", "lfq_2"], "aggregate": "mean_detected", "column": "lfq", "label": "lfq",
                      "n_rep_column": "n_rep_lfq"}],
         "combine": {"mode": "priority", "harmonise": "median_ratio"},
         "imputation": {"method": "none", "missing_label": "not_quantified"},
         "absent": [{"id": "nitc", "categories": {"meiosis": ["SPO11"], "glucose-repressed": ["GAL1"]}, "match_on": "gene_name",
                     "label": "not_in_this_cell"}],
         "per_protein_scaling": [{"rule_id": "size", "kind": "size_slope", "match_on": "orf", "ratio": {"numerator": 34, "denominator": 42},
                                  "table": {"path": os.path.join(d, "slopes.tsv"), "id_column": "orf", "slope_column": "slope",
                                            "default": "median"}}],
         "total_protein": {"mode": "fixed", "value_pg": 2.0},
         "adjustments": [{"rule_id": "ribo", "kind": "family_target", "match_on": "gene_name",
                          "members": {"regex": "^RP[LS]\\d+[AB]?$", "genes": ["ASC1"]}, "family_sub": ["[AB]$", ""],
                          "target": {"density_per_um3": 1000.0, "volume_um3": 5.0}}],
         "round_copies": False,
         "anchors": [{"id": "hist", "structure": "nucleosome", "families": {"H2A": {"stoich": 2, "genes": ["HTA1"]},
                                                                             "H2B": {"stoich": 2, "genes": ["HTB1"]},
                                                                             "H3": {"stoich": 2, "genes": ["HHT1"]},
                                                                             "H4": {"stoich": 2, "genes": ["HHF1"]}},
                      "counts": [{"label": "n", "value": 1000, "primary": True}, {"label": "n2", "value": 2000}]},
                     {"id": "nup", "stoich_map": {"NUP1": 8, "NUP2": 16, "SEC13": 16}, "filter_regex": "^NUP\\d+$",
                      "counts": [{"label": "pores", "value": 50}]}],
         "future_sources": [{"id": "sec", "type": "secms_abundance", "enabled": False}]}
    C.update(over)
    return load_config("testcell", overrides={"inventory": {"copies": C}})


def test_options():
    with tempfile.TemporaryDirectory() as d:
        _yeastlike(d)
        T, S = copies.step(StepContext(_opt_cfg(d), d, "copies"))
        T = T.reset_index().set_index("orf")
        g = T.set_index("gene")
        assert len(T) == 80 and list(T.index) == [f"Y{i:04d}W" for i in range(80)]          # annotation order
        assert g.at["GENE31", "source"] == "paxdb" and g.at["Y0000W", "source"] == "lfq" and g.at["GENE30", "source"] == "lfq"
        # lfq is mass-proportional: molar = mean intensity / mass, harmonised by the median ratio on shared proteins
        h = S["sources"]["lfq"]["harmonise_factor"]
        assert abs(h - g.at["GENE31", "ppm_pax"] / (1.5e6 / g.at["GENE31", "mass_da"])) < 1e-9 * h
        assert g.at["GENE30", "n_rep_lfq"] == 1 and g.at["GENE31", "n_rep_lfq"] == 2
        assert g.at["SPO11", "copies"] == 0 and g.at["GAL1", "copies"] == 0 and g.at["SPO11", "source"] == "not_in_this_cell"
        assert set(S["absent"]) == {"nitc.meiosis", "nitc.glucose-repressed"}
        # size slope: factor (34/42)^slope for ORFs in the table, median slope otherwise
        sl = pd.read_csv(os.path.join(d, "slopes.tsv"), sep="\t").set_index("orf").slope
        assert abs(T.at["Y0003W", "scale_factor"] - (34 / 42) ** sl["Y0003W"]) < 1e-12
        assert abs(T.at["Y0050W", "scale_factor"] - (34 / 42) ** sl.median()) < 1e-12
        # family target: median family copies (RPL1 = RPL1A + RPL1B, RPS2, ASC1) = 5000
        fam = [g.loc[["RPL1A", "RPL1B"], "copies"].sum(), g.at["RPS2", "copies"], g.at["ASC1", "copies"]]
        assert abs(np.median(fam) - 5000.0) < 1e-6 and S["adjustments"]["ribo"]["n_families"] == 3
        # histones: equal ppm and stoich -> their copies follow ppm x size factor; anchor ratio and implied total
        A = pd.read_csv(os.path.join(d, "inventory", "copies", "anchors.tsv"), sep="\t")
        r = A[(A.anchor_id == "hist") & (A.count_label == "n")].set_index("family").ratio
        assert abs(S["anchors"]["hist"]["n"]["ratio"] - float(np.median(r))) < 1e-12
        assert abs(S["anchors"]["hist"]["n"]["implied_total_pg"] - 2.0 / float(np.median(r))) < 1e-9
        assert set(A[A.anchor_id == "nup"].family) == {"NUP1", "NUP2"}                       # SEC13 filtered out
        L = pd.read_csv(os.path.join(d, "inventory", "copies", "copy_ledger.tsv"), sep="\t", keep_default_na=False)
        assert {"source.pax", "source.lfq", "absent.nitc.meiosis", "scale.size", "adjust.ribo", "total.fixed"} <= set(L.rule_id)
        assert (L[L.rule_id == "adjust.ribo"].key.isin(g.loc[["RPL1A", "RPL1B", "RPS2", "ASC1"]].index)).all()
        src = pd.read_csv(os.path.join(d, "inventory", "copies", "copy_sources.tsv"), sep="\t")
        assert src.set_index("source_id").at["sec", "status"] == "declared, not enabled"
    # total protein from anchors: with one anchor the scaled ratio becomes 1 (before adjustments; no adjustment here)
    with tempfile.TemporaryDirectory() as d:
        _yeastlike(d)
        cfg = _opt_cfg(d, adjustments=[], total_protein={"mode": "anchors", "anchors": ["hist"], "combine": "geometric_mean"})
        T, S = copies.step(StepContext(cfg, d, "copies"))
        assert abs(S["anchors"]["hist"]["n"]["ratio"] - 1.0) < 1e-9
        assert S["total_protein"]["mode"] == "anchors"
    # a future source cannot be switched on before it is implemented
    with tempfile.TemporaryDirectory() as d:
        _yeastlike(d)
        try:
            copies.step(StepContext(_opt_cfg(d, future_sources=[{"id": "sec", "type": "secms_abundance", "enabled": True}]), d, "copies"))
            raise AssertionError("enabled future source accepted")
        except NotImplementedError:
            pass


def test_configs():
    J = load_config("jurkat")["inventory"]
    C = J["copies"]
    assert J["adapters"]["proteome"] == "copies"
    assert C["total_protein"]["mode"] == "fixed" and J["total_protein_pg"] == 100 and C["round_copies"] is True
    assert C["sources"][0]["group_split"] == {"rule": "family_prefix", "n": 3} and C["imputation"]["method"] == "rna_loglog"
    assert not C.get("per_protein_scaling") and [r["rule_id"] for r in C["adjustments"]] == ["v8.copy_fixes"]
    assert not load_config("jurkat_v7_port")["inventory"]["copies"]["adjustments"]      # v7 recipe unchanged
    assert {a["id"] for a in C["anchors"]} == {"core_histones", "nucleoporins", "ribosomal_proteins"}
    assert all(not s.get("enabled") for s in C["future_sources"])
    Y = load_config("yeast")["inventory"]
    C = Y["copies"]
    assert Y["total_protein_pg"] == 3.5 and C["total_protein"]["mode"] == "fixed" and C["round_copies"] is False
    assert C["sources"][0]["reader"] == "paxdb" and C["per_protein_scaling"][0]["kind"] == "size_slope"
    assert C["adjustments"][0]["target"]["density_per_um3"] == 14250.0
    import yaml
    nitc = yaml.safe_load(open(os.path.join(REPO, "plugins", "yeast", "data", "inventory_copies.yaml")))["not_in_this_cell"]
    assert len(nitc) == 5 and sum(map(len, nitc.values())) == 4 + 8 + 43 + 11 + 8


def test_jurkat_copy_fixes():
    """v8.copy_fixes (DECISION 2026-11-19 evening): one post-scaling adjustment (family_target on PDGFRB alone, so no other
    protein changes); target = the RNA-based estimate (undetected `low` call value min(fit, floor) x low_factor, Step-2
    fit) in copies at 100 pg; SOS1 / SOS2 not touched; on a synthetic cell only the member's copies change"""
    C = load_config("jurkat")["inventory"]["copies"]
    (rule,) = C["adjustments"]
    assert rule["kind"] == "family_target" and rule["match_on"] == "key" and rule["members"] == {"genes": ["PDGFRB"]}
    a, b, floor, low, tpm = 2.7320674797595044, 0.7679479131312847, 97.35650449372179, 0.25, 0.19534589452783302
    est = min(10 ** (a + b * np.log10(tpm)), floor) * low                      # ppb
    assert C["imputation"]["low_factor"] == low and abs(est - 24.339126123430447) < 1e-9
    n_mol, sum_ppb = 1595375332.6162028, 952062291.1797417                     # Step-2 copies summary / copies.tsv
    assert abs(rule["target"]["count"] - est / sum_ppb * n_mol) < 1e-6 and round(rule["target"]["count"]) == 41
    assert load_config("jurkat_v7_port")["inventory"]["copies"]["adjustments"] == []
    with tempfile.TemporaryDirectory() as d:
        _yeastlike(d)
        T0 = copies.step(StepContext(_opt_cfg(d), d, "copies"))[0]
        adj = [dict(a) for a in _opt_cfg(d)["inventory"]["copies"]["adjustments"]]
    cand = T0[(T0.copies > 0) & ~T0.gene_name.fillna("").str.match(r"^(RP[LS]|ASC1)")]
    g = cand.index[0]
    with tempfile.TemporaryDirectory() as d:
        _yeastlike(d)
        T1 = copies.step(StepContext(_opt_cfg(d, adjustments=adj + [dict(rule, members={"genes": [g]}, target={"count": 3.0})]),
                                     d, "copies"))[0]
    assert abs(T1.copies[g] - 3.0) < 1e-9
    o = T0.index.drop(g)
    assert (T1.copies[o] == T0.copies[o]).all()                                 # nobody else changes


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
