"""Surface-proteomics adapter tests (seconds; synthetic tables only, no lab data needed).

    PYTHONHASHSEED=0 python tests/test_surface.py        (or pytest tests/test_surface.py)

  test_read_dedup_and_corrupt_rows   rows repeated per comparison collapse to one value per (protein, run); a shifted row
                                     with an unknown condition is dropped; conflicting repeats raise
  test_recentre_removes_drift        a global per-condition drift is removed (complete-median re-centring)
  test_presence_mode                 surface.tsv columns, detection evidence, topology from the UniProt table, no fractions
  test_topology_partial_annotation   TM sides alternate; phase from TOPO_DOM / N-glycans / signal peptide (TRBC-like
                                     tail-only annotations completed)
  test_pm_share_requires_validation  pm_share refuses an unvalidated table unless allow_unvalidated; validated table sets
                                     pm_fraction
  test_timepoint_factor_table        only robust changes (size, replicate agreement, detection, PM gene filter) become
                                     factors 2^log2_change
  test_registered_and_presence_only            adapter registered; Jurkat config keeps it OFF (adapters.surface null, slot disabled)
                                     and the config fragment parses and names the adapter
"""
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.inventory import schema, surface  # noqa: E402

COND = ["NAct", "Act_0", "Act_24h"]
RNG = np.random.default_rng(0)
TMP = tempfile.mkdtemp(prefix="surface_test_")


def _synthetic(drift=0.0, n=60):
    """long table: n proteins x 3 conditions x 4 replicates, repeated for 2 comparisons; protein P0 up 4x at 24h, P1 down
    4x, P2 up 4x but replicate-inconsistent, P3 absent at NAct; + one corrupt shifted row"""
    base = RNG.normal(23, 1.5, n)
    rows = []
    for ci, c in enumerate(COND):
        for rep in range(1, 5):
            run = f"run_{c}_{rep}"
            for i in range(n):
                v = base[i] + RNG.normal(0, 0.1) + (drift if c == "Act_24h" else 0.0)
                if c == "Act_24h" and i == 0: v += 2
                if c == "Act_24h" and i == 1: v -= 2
                if c == "Act_24h" and i == 2: v += (2 if rep <= 2 else -2.2)
                if c == "NAct" and i == 3: v = np.nan
                for comp in ("A_vs_NAct", "A_vs_Act_0"):
                    rows.append(dict(entry=f"P{i}", condition=c, run=run, bio_replicate=rep, comparison=comp,
                                     normalised_protein_intensity_log2=v, n_features_per_sample=3 + i % 5))
    d = pd.DataFrame(rows)
    d = pd.concat([d, pd.DataFrame([dict(entry="440", condition="21.0145", run="HIDE1_HUMAN", bio_replicate="Protein",
                                         comparison="no", normalised_protein_intensity_log2=np.nan)])])
    p = os.path.join(TMP, f"long_{drift}.csv")
    d.to_csv(p, index=False)
    return p


def _annotation():
    rows = []
    for i in range(60):
        if i % 3 == 0:     # type I
            rows.append(dict(Entry=f"P{i}", **{"Gene Names (primary)": f"G{i}", "Topological domain":
                         'TOPO_DOM 20..200; /note="Extracellular"; TOPO_DOM 222..300; /note="Cytoplasmic"',
                         "Transmembrane": 'TRANSMEM 201..221; /note="Helical"', "Lipidation": "", "Signal peptide": "SIGNAL 1..19"}))
        elif i % 3 == 1:   # GPI
            rows.append(dict(Entry=f"P{i}", **{"Gene Names (primary)": f"G{i}", "Topological domain": "", "Transmembrane": "",
                         "Lipidation": 'LIPID 100; /note="GPI-anchor amidated serine"', "Signal peptide": "SIGNAL 1..22"}))
        else:              # type II
            rows.append(dict(Entry=f"P{i}", **{"Gene Names (primary)": f"G{i}", "Topological domain":
                         'TOPO_DOM 1..40; /note="Cytoplasmic"; TOPO_DOM 62..400; /note="Extracellular"',
                         "Transmembrane": 'TRANSMEM 41..61; /note="Helical; Signal-anchor for type II membrane protein"',
                         "Lipidation": "", "Signal peptide": ""}))
    p = os.path.join(TMP, "uniprot.tsv")
    pd.DataFrame(rows).to_csv(p, sep="\t", index=False)
    return p


def _params(path, **kw):
    p = {"path": path, "conditions": COND, "annotation": _annotation(), "baseline": ["NAct", "Act_0"]}
    p.update(kw)
    return p


ident = lambda x: x  # noqa: E731


def test_read_dedup_and_corrupt_rows():
    path = _synthetic()
    o = surface.options(_params(path))
    M, F, S, genes, info = surface.read_long(path, o)
    assert M.shape == (60, 12) and info["rows_dropped_condition"] == 1
    assert info["pairs_measured"] == 60 * 12 - 4                       # P3 missing in the 4 NAct runs
    assert set(S.condition) == set(COND) and F.loc["P4"].max() == 3 + 4 % 5
    d = pd.read_csv(path)
    d.loc[(d.entry == "P5") & (d.comparison == "A_vs_Act_0") & (d.run == "run_NAct_1"), "normalised_protein_intensity_log2"] += 1
    bad = os.path.join(TMP, "bad.csv"); d.to_csv(bad, index=False)
    try:
        surface.read_long(bad, o)
    except ValueError as e:
        assert "differs" in str(e)
    else:
        raise AssertionError("conflicting repeats accepted")


def test_recentre_removes_drift():
    path = _synthetic(drift=0.8)
    o = surface.options(_params(path))
    M, _, S, _, _ = surface.read_long(path, o)
    Mc, shift = surface.recentre(M, "complete_median")
    s = shift.groupby(S.condition).mean()
    assert abs((s["Act_24h"] - s["NAct"]) - 0.8) < 0.1
    ch = surface.paired_change(Mc, S, "Act_24h", "NAct")
    assert abs(ch.log2_change.drop(["P0", "P1", "P2", "P3"]).median()) < 0.05


def test_presence_mode():
    out, info = surface.surface_table(ident, _params(_synthetic()))
    assert list(out.columns) == schema.SURFACE
    assert set(out.condition) == {"NAct", "Act_0"} and out.pm_fraction.isna().all()
    assert len(out[out.condition == "NAct"]) == 59                       # P3 not measured at NAct -> no row
    r = out[(out.uniprot == "P0") & (out.condition == "Act_0")].iloc[0]
    assert r.topology == "single_pass_typeI" and r.n_term_side == "extracellular" and r.gene == "G0"
    assert out[out.uniprot == "P1"].topology.iloc[0] == "GPI"
    r2 = out[out.uniprot == "P2"].iloc[0]
    assert r2.topology == "single_pass_typeII" and r2.n_term_side == "cytoplasmic"
    assert r.evidence.startswith("detected 4/4") and info["detected"]["NAct"] == 59


def test_topology_partial_annotation():
    # TRBC2-like: one TM, only the C-terminal cytoplasmic tail annotated -> N-terminal constant domain exoplasmic
    sd = surface.infer_sides(178, [(145, 167)], [(168, 178, "Cytoplasmic")])
    assert set(sd[:144]) == {"e"} and set(sd[167:]) == {"c"}
    # no TOPO_DOM: phase from an annotated N-glycan, else from the signal peptide; multi-pass alternation
    sd = surface.infer_sides(300, [(50, 70), (100, 120)], [], nglyc=[30])
    assert sd[29] == "e" and sd[85] == "c" and sd[200] == "e"
    sd = surface.infer_sides(300, [(50, 70)], [], sp_end=20)
    assert sd[10] == "" and sd[30] == "e" and sd[100] == "c"
    t = surface._topology_rows(pd.DataFrame([{"Entry": "X", "Gene Names (primary)": "TRBX", "Length": 178,
                                              "Topological domain": 'TOPO_DOM 168..178; /note="Cytoplasmic"',
                                              "Transmembrane": 'TRANSMEM 145..167; /note="Helical"'}]))
    assert t.loc["X", "topology"] == "single_pass_typeI" and t.loc["X", "n_term_side"] == "extracellular"
    assert t.loc["X", "exo_segments"] == "1-144"


def test_pm_share_requires_validation():
    path = _synthetic()
    t = os.path.join(TMP, "pm.tsv")
    pd.DataFrame({"uniprot": ["P0", "P1"], "pm_fraction": [0.8, 0.3]}).to_csv(t, sep="\t", index=False)
    try:
        surface.surface_table(ident, _params(path, mode="pm_share", pm_share_table=t))
    except ValueError as e:
        assert "validated" in str(e)
    else:
        raise AssertionError("unvalidated PM-share table accepted")
    out, _ = surface.surface_table(ident, _params(path, mode="pm_share", pm_share_table=t, allow_unvalidated=True))
    assert out.set_index(["uniprot", "condition"]).pm_fraction[("P0", "NAct")] == 0.8
    tv = os.path.join(TMP, "pm_valid.tsv")
    with open(tv, "w") as fh:
        fh.write("# validated: true\n")
        pd.DataFrame({"uniprot": ["P0"], "pm_fraction": [0.5]}).to_csv(fh, sep="\t", index=False)
    out, _ = surface.surface_table(ident, _params(path, mode="pm_share", pm_share_table=tv))
    assert out[out.uniprot == "P0"].pm_fraction.iloc[0] == 0.5 and out[out.uniprot == "P5"].pm_fraction.isna().all()


def test_timepoint_factor_table():
    path = _synthetic(drift=0.8)
    p = _params(path, timepoint={"condition": "Act_24h", "reference": "NAct"})
    ft = surface.factor_table(ident, p)
    assert set(ft.uniprot) == {"P0", "P1"}, ft                         # P2 inconsistent, P3 not detected at NAct
    f = ft.set_index("uniprot").factor
    assert abs(np.log2(f["P0"]) - 2) < 0.2 and abs(np.log2(f["P1"]) + 2) < 0.2
    p["timepoint"]["pm_only_genes"] = ["G1"]
    assert surface.factor_table(ident, p).gene.tolist() == ["G1"]
    out, _ = surface.surface_table(ident, _params(path, mode="timepoint",
                                                  timepoint={"condition": "Act_24h", "reference": "NAct"}))
    assert set(out.condition) == {"NAct", "Act_24h"}
    import re
    ev = out[(out.uniprot == "P0")].evidence.iloc[0]
    assert abs(float(re.search(r"log2 Act_24h/NAct ([+-][\d.]+)", ev).group(1)) - 2) < 0.2, ev


def test_registered_and_presence_only():
    import yaml
    assert registry.get("adapter.surface", "surface_proteomics") is surface.surface_proteomics
    cfg = yaml.safe_load(open(os.path.join(REPO, "config", "cells", "jurkat.yaml")))
    inv = cfg["inventory"]
    # freeze v8 D7 (2026-11-19): enabled, presence evidence only (never pm_share)
    assert inv["adapters"]["surface"] == "surface_proteomics"
    assert inv["params"]["surface_proteomics"]["mode"] == "presence"
    assert "per_protein_scaling" not in str(inv["copies"].get("per_protein_scaling", "")) or "surface" not in str(inv["copies"].get("per_protein_scaling"))
    slot = [s for s in inv["copies"]["future_sources"] if s["type"] == "surface_proteomics"]
    assert slot and not slot[0]["enabled"]          # copies slot never enabled: presence acts via adapter.surface only
    frag = yaml.safe_load(open(os.path.join(REPO, "config", "fragments", "jurkat_surface_proteomics.yaml")))
    assert frag["inventory"]["adapters"]["surface"] == "surface_proteomics"
    sp = frag["inventory"]["params"]["surface_proteomics"]
    assert sp["mode"] == "presence" and surface.options(sp)["mode"] == "presence"
    fi = frag["frozen_inputs"]["inv_surface_jurkat_glyco2026"]
    assert fi["redistributable"] in (True, "yes") and "CC BY-NC" in fi["licence"]   # lab data public since 2026-11-23


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
