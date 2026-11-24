"""Inventory stage tests: the testcell chain end to end (seconds) and the schema of every canonical table.

    PYTHONHASHSEED=0 python tests/test_inventory.py        (or pytest tests/test_inventory.py)
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
from wholecell.inventory import schema, stage  # noqa: E402


def _run_testcell(bdir):
    cfg = load_config("testcell")
    registry.load_plugins(cfg)
    stage.run(cfg, bdir)
    return cfg, os.path.join(bdir, "inventory")


def test_testcell_chain():
    with tempfile.TemporaryDirectory() as d:
        cfg, inv = _run_testcell(d)
        plan = stage.steps(cfg)
        assert list(plan) == ["v0", "secms", "assembly", "ingredients", "rules:generic.signal_lumen", "pack_table", "canonical"]
        for s, spec in plan.items():
            for o in spec["outputs"]:
                assert os.path.exists(os.path.join(d, o)), (s, o)
        P = pd.read_csv(os.path.join(inv, "pack_table.tsv"), sep="\t")
        assert list(P.columns) == schema.PACK_TABLE
        ids = {l.strip() for l in open(os.path.join(REPO, "tests", "data", "testcell", "library_ids.txt")) if l.strip() and not l.startswith("#")}
        assert set(P.ing_id) == ids, set(P.ing_id) ^ ids                 # every species has a library entry
        assert P.copies.sum() < 5_000_000
        assert abs(P.mass_pg.sum() - cfg["inventory"]["total_protein_pg"]) / cfg["inventory"]["total_protein_pg"] < 0.02
        for w in P.comp_weights:
            assert abs(sum(json.loads(w).values()) - 1) < 2e-3
        pgk = json.loads(P.set_index("ing_id").at["PGK1", "comp_weights"])
        assert "er" not in pgk and pgk["cytosol"] == 1.0                  # generic.signal_lumen applied
        hspa5 = json.loads(P.set_index("ing_id").at["HSPA5", "comp_weights"])
        assert hspa5 == {"er": 1.0}                                        # signal peptide: lumen kept
        ch = pd.read_csv(os.path.join(inv, "rules", "generic.signal_lumen", "changes.tsv"), sep="\t")
        assert list(ch.columns[:len(schema.CHANGES)]) == schema.CHANGES and set(ch.species_id) == {"PGK1"}
        assert (ch.rule_id == "generic.signal_lumen").all()
        # canonical tables: present, columns as in the schema
        for name, cols in (("species.tsv", schema.SPECIES), ("species_loc.tsv", schema.SPECIES_LOC),
                           ("field_provenance.tsv", schema.FIELD_PROVENANCE), ("sources.tsv", schema.SOURCES),
                           ("proteome_sources.tsv", schema.PROTEOME_SOURCES), ("localisation_sources.tsv", schema.LOCALISATION_SOURCES),
                           ("gene_abundance.tsv", schema.GENE_ABUNDANCE), ("secms_protein.tsv", schema.SECMS_PROTEIN),
                           ("secms_complex.tsv", schema.SECMS_COMPLEX), ("surface.tsv", schema.SURFACE),
                           ("changes_all.tsv", schema.CHANGES), ("version_changes.tsv", schema.VERSION_CHANGES),
                           ("secms_bound_fraction.tsv", schema.SECMS_BOUND_FRACTION)):
            T = pd.read_csv(os.path.join(inv, name), sep="\t")
            assert list(T.columns[:len(cols)]) == cols, name
        S = pd.read_csv(os.path.join(inv, "species.tsv"), sep="\t").set_index("species_id")
        assert S.at["CPX-5993", "species_key"] == "cpx:CPX-5993" and S.at["GAPDH_x4", "kind"] == "complex"
        assert S.at["PGK1", "species_key"].startswith("prot:")
        L = pd.read_csv(os.path.join(inv, "species_loc.tsv"), sep="\t")
        assert L[(L.species_id == "PGK1") & (L.loc_key == "cytosol")].rule_id.item() == "generic.signal_lumen"
        assert abs(L.groupby("species_id").weight.sum() - 1).max() < 2e-3
        VC = pd.read_csv(os.path.join(inv, "version_changes.tsv"), sep="\t")
        assert "generic.signal_lumen" in set(VC.rule_ids.dropna())
        for name in ("secms_protein.tsv", "secms_complex.tsv", "secms_complex_subunit.tsv", "secms_pairs.tsv",
                     "secms_calibration.tsv"):
            assert len(pd.read_csv(os.path.join(inv, name), sep="\t")) > 0         # synthetic SEC-MS (toolkits/secms)
        assert len(pd.read_csv(os.path.join(inv, "surface.tsv"), sep="\t")) == 0   # empty but present


def test_total_protein_required():
    with tempfile.TemporaryDirectory() as d:
        cfg = load_config("testcell")
        del cfg["inventory"]["total_protein_pg"]
        try:
            stage.run(cfg, d, step="v0")
        except ValueError as e:
            assert "total_protein_pg" in str(e)
        else:
            raise AssertionError("v0 ran without inventory.total_protein_pg")


def test_jurkat_plan():
    """the Jurkat step plan is buildable from the config (no data are read)"""
    cfg = load_config("jurkat")
    registry.load_plugins(cfg)
    plan = stage.steps(cfg)
    assert [s for s in plan if s.startswith("rules:")] == ["rules:jurkat.v4", "rules:jurkat.v5", "rules:jurkat.v6", "rules:jurkat.v7",
                                                           "rules:jurkat.v7_lamina"]
    v7 = plan["rules:jurkat.v7"]
    assert "inv_lamina_consumption_v7" not in v7["frozen"] and "inventory/rules/jurkat.v6/ingredients.tsv" in v7["inputs"]
    # v9: the build's own lamina (no legacy lattice); since faf17e1 rule L3 (+ L2) is the rule set jurkat.v7_lamina, after
    # the library (no DAG cycle library -> lamina -> rules -> pack table)
    lam = plan["rules:jurkat.v7_lamina"]
    assert "preplaced/lamina/consumption.json" in lam["inputs"] and "inventory/rules/jurkat.v7/ingredients.tsv" in lam["inputs"]
    assert "preplaced/lamina/consumption.json" not in v7["inputs"]
    assert "inv_v4pack_loc_compare" in plan["rules:jurkat.v5"]["frozen"]
    assert {"inv_v5pack_loc_compare", "inv_tier2_v5"} <= set(plan["rules:jurkat.v6"]["frozen"])
    for s in plan.values():
        for n in s["frozen"]:
            assert n in cfg["frozen_inputs"], n


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "ok")
