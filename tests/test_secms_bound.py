"""Per-protein SEC-MS bound fraction (wholecell/inventory/secms_bound.py; canonical secms_bound_fraction.tsv; spec
docs/reference/inventory_schema.md §4c): synthetic tables with known answers + the testcell chain + the simbox columns.

    PYTHONHASHSEED=0 python tests/test_secms_bound.py        (or pytest tests/test_secms_bound.py)
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

from wholecell import registry  # noqa: E402
from wholecell.config import load as load_config  # noqa: E402
from wholecell.inventory import schema, secms_bound, stage  # noqa: E402
from wholecell.inventory.context import StepContext  # noqa: E402


def _w(d, rel, df):
    p = os.path.join(d, "inventory", rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    df.to_csv(p, sep="\t", index=False)


def _synthetic(d, oligo_map=None, extra_source=False):
    """A: monomer only (a_2 0.8); B: 60 % in a heterodimer A-less complex; C: homo-tetramer (k 3 by state) at 50 %;
    D: not measured; E: family-held (no species names it); F: zero copies; 'G1; G2': yeast-style merged label, monomer
    row + complex naming G1 only; H: yeast-style library oligomer (oligo_map) 100 %"""
    P = pd.DataFrame([
        ("A", "protein", "A", "1", 1000.0), ("B", "protein", "B", "1", 400.0), ("BX", "complex", "B;X", "1;1", 600.0),
        ("C", "protein", "C", "1", 500.0), ("C_x4", "oligomer", "C", "4", 125.0), ("D", "protein", "D", "1", 100.0),
        ("FAM", "complex", "FAMILY", "2", 50.0), ("G1; G2", "protein", "G1; G2", "1", 300.0),
        ("GX", "complex", "G1;X", "1;1", 100.0), ("H_oligo", "protein", "H", "1", 800.0), ("X", "protein", "X", "1", 10.0),
    ], columns=["ing_id", "kind", "genes", "stoich", "copies"])
    for c in schema.PACK_TABLE:
        if c not in P.columns:
            P[c] = 0
    _w(d, "pack_table.tsv", P[schema.PACK_TABLE])
    V = pd.DataFrame([("A", "PA", 1000.0), ("B", "PB", 1000.0), ("C", "PC", 1000.0), ("D", "PD", 100.0),
                      ("E", "PE", 100.0), ("F", "PF", 0.0), ("G1; G2", "PG", 400.0), ("H", "PH", 800.0),
                      ("X", "PX", 710.0)], columns=["gene", "uniprot", "copies"]).assign(mass_da=50000.0)
    _w(d, "v0/inventory_v0.tsv", V)
    rows = []
    a = {"PA": {2: 0.8, 3: 0.7, 4: 0.6}, "PB": {2: 0.9, 3: 0.5, 4: 0.4}, "PC": {2: 0.95, 3: 0.5, 4: 0.3},
         "PE": {2: 0.9, 3: 0.9, 4: 0.9}, "PF": {2: 0.5, 3: 0.5, 4: 0.5}, "PG": {2: 0.6, 3: 0.6, 4: 0.6},
         "PH": {2: 0.7, 3: 0.7, 4: 0.7}, "PX": {2: np.nan, 3: np.nan, 4: np.nan}}
    srcs = [("run1", "c1", "1")] + ([("run2", "c2", "1")] if extra_source else [])
    for src, cond, rep in srcs:
        for u, ks in a.items():
            for k, v in ks.items():
                rows.append(dict(uniprot=u, gene=u[1:], source_id=src, condition=cond, replicate=rep, threshold_k=k,
                                 assembled_frac=v, assembled_frac_sd=0.05, assembled_frac_lo=max(v - 0.1, 0) if v == v else v,
                                 assembled_frac_hi=min(v + 0.1, 1) if v == v else v, n_rep=2, apex_app_kda=500.0, mono_kda=50.0))
    _w(d, "secms/secms_protein.tsv", pd.DataFrame(rows))
    if oligo_map is not None:
        p = os.path.join(d, "inventory", "assembly", "oligo_map.json")
        os.makedirs(os.path.dirname(p), exist_ok=True)
        json.dump(oligo_map, open(p, "w"))


def _compute(d, overrides=None):
    cfg = load_config("testcell", overrides=overrides)
    ctx = StepContext(cfg, d, "canonical")
    P = pd.read_csv(os.path.join(d, "inventory", "pack_table.tsv"), sep="\t")
    V = pd.read_csv(os.path.join(d, "inventory", "v0", "inventory_v0.tsv"), sep="\t", dtype={"gene": str}).set_index("gene")
    return secms_bound.compute(ctx, P, V)


def test_synthetic_known_answers():
    with tempfile.TemporaryDirectory() as d:
        _synthetic(d, oligo_map={"H_oligo": ["H_lib", 2]})
        F, S = _compute(d)
        assert list(F.columns) == schema.SECMS_BOUND_FRACTION
        F = F.set_index("gene")
        # A: monomer only, k 2 (rules default), gap = a_2
        assert F.at["A", "model_state"] == "monomer_only" and F.at["A", "structure_represented"] == False  # noqa: E712
        assert F.at["A", "threshold_k"] == 2 and abs(F.at["A", "unrepresented_bound_fraction"] - 0.8) < 1e-12
        assert abs(F.at["A", "unrepresented_share_of_free"] - 0.8) < 1e-12
        assert F.at["A", "secms_bound_fraction_sd"] == 0.05 and abs(F.at["A", "secms_bound_fraction_lo"] - 0.7) < 1e-12
        assert F.at["A", "datasets"] == "run1|c1|1" and F.at["A", "secms_source_id"] == "run1"
        # B: 60 % modelled in a complex; a_2 0.9 -> gap 0.3; share of the 400 free copies 0.3 / 0.4
        assert F.at["B", "model_state"] == "in_assembly" and F.at["B", "model_species"] == "BX"
        assert abs(F.at["B", "modelled_assembled_fraction"] - 0.6) < 1e-12
        assert abs(F.at["B", "unrepresented_bound_fraction"] - 0.3) < 1e-12
        assert abs(F.at["B", "unrepresented_share_of_free"] - 0.75) < 1e-12
        # C: homo-tetramer -> k 3 (k_by_state n <= 5); a_3 0.5 = modelled 0.5 -> no gap
        assert F.at["C", "threshold_k"] == 3 and F.at["C", "model_homo_state"] == 4
        assert abs(F.at["C", "secms_bound_fraction"] - 0.5) < 1e-12 and F.at["C", "unrepresented_bound_fraction"] == 0
        # D: not measured -> NA with reason; model columns still filled
        assert np.isnan(F.at["D", "secms_bound_fraction"]) and "not quantified" in F.at["D", "na_reason"]
        assert F.at["D", "modelled_assembled_fraction"] == 0
        # E: no species names it, copies held elsewhere -> represented, modelled 1, no gap
        assert F.at["E", "model_state"] == "not_in_pack" and F.at["E", "structure_represented"] == True  # noqa: E712
        assert F.at["E", "modelled_assembled_fraction"] == 1 and F.at["E", "unrepresented_bound_fraction"] == 0
        # F: zero copies
        assert F.at["F", "model_state"] == "absent" and np.isnan(F.at["F", "modelled_assembled_fraction"])
        assert np.isnan(F.at["F", "unrepresented_bound_fraction"]) and F.at["F", "na_reason"]
        # 'G1; G2': merged label matched to its monomer row; complex naming G1 counts; modelled 1 - 300/400
        assert F.at["G1; G2", "model_state"] == "in_assembly" and abs(F.at["G1; G2", "modelled_assembled_fraction"] - 0.25) < 1e-12
        assert abs(F.at["G1; G2", "unrepresented_bound_fraction"] - 0.35) < 1e-12
        # H: library oligomer of oligo_map (kind protein) is an assembly, state 2 -> k 2
        assert F.at["H", "model_state"] == "in_assembly" and F.at["H", "model_homo_state"] == 2
        assert F.at["H", "modelled_assembled_fraction"] == 1 and F.at["H", "unrepresented_bound_fraction"] == 0
        # X: a_k NaN
        assert "NaN" in F.at["X", "na_reason"] and np.isnan(F.at["X", "unrepresented_bound_fraction"])
        # summary: unrepresented mass = sum(gap x mass) over measured proteins
        m = F.copies * 50000.0 / schema.AVOGADRO * 1e12
        assert abs(S["unrepresented_bound_pg"] - float((F.unrepresented_bound_fraction.fillna(0) * m).sum())) < 1e-15
        assert S["top"][0]["gene"] == "A"


def test_synthetic_k_override_and_selection():
    with tempfile.TemporaryDirectory() as d:
        _synthetic(d, extra_source=True)
        F, S = _compute(d)                         # two runs, no selection: NA everywhere with a reason (no failure)
        assert F.secms_bound_fraction.isna().all() and F.na_reason.str.contains("not selected").all()
        F, S = _compute(d, {"inventory": {"canonical": {"secms_bound_fraction": {"source_id": "run2", "k": 4}}}})
        F = F.set_index("gene")
        assert (F.threshold_k == 4).all() and F.at["A", "secms_source_id"] == "run2"
        assert abs(F.at["A", "secms_bound_fraction"] - 0.6) < 1e-12
        assert S["selection"]["source_id"] == "run2"
    with tempfile.TemporaryDirectory() as d:
        _synthetic(d)
        F, _ = _compute(d, {"inventory": {"adapters": {"secms": None}}})
        assert F.secms_bound_fraction.isna().all() and F.na_reason.str.startswith("no SEC-MS step").all()
        assert F.modelled_assembled_fraction.notna().sum() == len(F) - 1   # model side filled without SEC-MS


def test_testcell_chain_and_simbox_columns():
    from wholecell.export.stage import secms_bound_columns
    with tempfile.TemporaryDirectory() as d:
        cfg = load_config("testcell")
        registry.load_plugins(cfg)
        stage.run(cfg, d)
        inv = os.path.join(d, "inventory")
        F = pd.read_csv(os.path.join(inv, "secms_bound_fraction.tsv"), sep="\t", dtype={"gene": str})
        V = pd.read_csv(os.path.join(inv, "v0", "inventory_v0.tsv"), sep="\t", dtype={"gene": str})
        assert list(F.columns) == schema.SECMS_BOUND_FRACTION and list(F.gene) == list(V.gene.drop_duplicates())
        F = F.set_index("gene")
        # measured value = the secms_protein row the assembly rules read (testcell: the only run, replicate mean) at k
        T = pd.read_csv(os.path.join(inv, "secms", "secms_protein.tsv"), sep="\t", dtype={"replicate": str})
        T = T[T.replicate == "mean"].set_index(["uniprot", "threshold_k"]).assembled_frac
        for g in ("GAPDH", "PGK1", "ENO1", "CCT2"):
            assert abs(F.at[g, "secms_bound_fraction"] - T[(F.at[g, "uniprot"], int(F.at[g, "threshold_k"]))]) < 1e-12, g
        assert F.at["GAPDH", "threshold_k"] == 3 and F.at["ENO1", "threshold_k"] == 2   # k_by_state of the rules
        assert F.at["GAPDH", "modelled_assembled_fraction"] == 1 and F.at["PGK1", "modelled_assembled_fraction"] == 0
        x = F[F.secms_bound_fraction.notna() & (F.copies > 0)]
        gap = (x.secms_bound_fraction - x.modelled_assembled_fraction).clip(lower=0)
        assert np.allclose(gap, x.unrepresented_bound_fraction)
        assert (x.structure_represented == (x.model_state != "monomer_only")).all()
        S = json.load(open(os.path.join(inv, "secms_bound_fraction_summary.json")))
        assert S["selection"]["replicate"] == "mean" and 0 < S["unrepresented_bound_frac"] < 1
        # simbox columns: monomer species get their protein's values, assemblies NaN
        sb, meta = secms_bound_columns(["PGK1", "GAPDH_x4", "HSP90AB1", "prepass_thing"], inv)
        assert meta["status"] == "ok" and meta["threshold_k"]
        assert abs(sb["HSP90AB1"][1] - F.at["HSP90AB1", "unrepresented_share_of_free"]) < 1e-12
        assert np.isnan(sb["GAPDH_x4"][0]) and np.isnan(sb["prepass_thing"][0])
        sb, meta = secms_bound_columns(["PGK1"], d)        # no table -> NaN + status
        assert np.isnan(sb["PGK1"][0]) and meta["status"].startswith("unavailable")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "ok")
