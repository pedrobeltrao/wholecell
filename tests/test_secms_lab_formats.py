"""SEC-MS toolkit: protein-level formats of the lab 2026 datasets, shared calibration, condition combination.
Synthetic only (the testcell SEC-MS); the lab data are confidential and never used in tests.

    PYTHONHASHSEED=0 python tests/test_secms_lab_formats.py        (or pytest)
"""
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from toolkits.secms import combine, io, pipeline  # noqa: E402
from wholecell.config import load as load_config  # noqa: E402
from wholecell.inventory import secms  # noqa: E402

_R = {}


def _base():
    """testcell toolkit config + its peptide-level run (reference) + protein-level raw profiles of replicate 1"""
    if "b" not in _R:
        cfg = load_config("testcell")
        p = secms._frozen_paths(dict(cfg["inventory"]["params"]["secms_toolkit"]),
                                type("C", (), {"cfg": cfg, "frozen": cfg.frozen})())
        p.pop("pair_candidates", None)
        p["complexes"] = [c for c in p["complexes"] if c["format"] == "complex_portal"]
        r = pipeline.run(p)
        pf = r["profiles"]
        raw = pf[(pf.kind == "raw") & (pf.replicate.astype(str) == "1")].set_index("uniprot")
        fc = [c for c in raw.columns if c.startswith("f") and c[1:].isdigit()]
        _R["b"] = (p, r, raw[fc].fillna(0))
    return _R["b"]


def _ds(p, **kw):
    d = {k: v for k, v in p["datasets"][0].items() if k not in ("columns", "format", "path")}
    d.update(kw)
    return d


def _prot(res, sid, k=2):
    P = res["secms_protein"]
    return P[(P.source_id == sid) & (P.threshold_k == k)].set_index("uniprot").assembled_frac


def test_protein_long_repeated_pg_quantity():
    """precursor-level export: the protein-group value repeated on n rows is taken once; unmapped 'NA' ids dropped"""
    p, _, raw = _base()
    rng = np.random.default_rng(0)
    rows = []
    for acc, v in raw.iterrows():
        for i, x in enumerate(v.values):
            if x > 0:
                rows += [(acc, i + 1, x)] * int(rng.integers(1, 6))
    rows += [("NA", 3, 1e9)] * 4 + [("NA", 4, np.nan)] * 2
    T = pd.DataFrame(rows, columns=["input_id.y", "r_fraction", "pg_ibaq"])
    T.insert(0, "r_condition", "C1")
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "long.tsv")
        T.to_csv(f, sep="\t", index=False, na_rep="NA")
        cols = {"protein": "input_id.y", "fraction": "r_fraction", "intensity": "pg_ibaq", "condition": "r_condition"}
        ds = _ds(p, source_id="long", format="protein_long", path=f, columns=cols, duplicate_rows="unique",
                 exclude_ids_regex="^NA$")
        kind, X = io.load(ds, {}, lambda x: x)
        assert kind == "protein" and "NA" not in set(X.protein)
        piv = X.pivot_table(index="protein", columns="fraction", values="intensity", aggfunc="sum")
        assert np.allclose(piv.reindex(index=raw.index, columns=range(1, raw.shape[1] + 1)).fillna(0).values, raw.values)
        res = pipeline.run(dict(p, datasets=[ds]))
        ds_sum = dict(ds, duplicate_rows="sum")                  # summing the repeats would distort the profiles
        bad = pipeline.run(dict(p, datasets=[dict(ds_sum, source_id="sum")]))
        a, b = _prot(res, "long"), _prot(bad, "sum")
        assert (a - b).abs().max() > 0.01
        T2 = T.copy()
        T2.loc[T2.index[1], "pg_ibaq"] = T2.pg_ibaq.iloc[1] + 1  # two different values in one cell -> error
        T2.to_csv(f, sep="\t", index=False, na_rep="NA")
        try:
            io.load(ds, {}, lambda x: x)
        except ValueError as e:
            assert "different" in str(e)
        else:
            raise AssertionError("conflicting duplicate values accepted")


def test_protein_matrix_blanks_and_groups():
    """wide ID x 1..n matrix, blank = missing (0), '[P1;P2]' groups expanded unless a member has its own row"""
    p, _, raw = _base()
    accs = list(raw.index)
    M = raw.replace(0, np.nan).copy()
    M.columns = [str(i + 1) for i in range(M.shape[1])]
    g = M.loc[[accs[0]]].copy(); g.index = [f"[{accs[0]};Q99999]"]   # group whose first member is also a row
    h = M.loc[[accs[1]]].copy(); h.index = ["[Q88888;Q77777]"]       # group of two otherwise absent accessions
    M = pd.concat([M, g, h]); M.index.name = "ID"
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "wide.tsv")
        M.to_csv(f, sep="\t", na_rep="")
        ds = _ds(p, source_id="wide", format="protein_matrix", path=f, id_col="ID", value_regex=r"^(?P<fraction>\d+)$",
                 id_groups={"strip": "[]", "sep": ";", "mode": "expand"})
        _, X = io.load(ds, {}, lambda x: x)
        assert {"Q99999", "Q88888", "Q77777"} <= set(X.protein) and not X.protein.str.contains(r"\[|;").any()
        assert np.isclose(X[X.protein == accs[0]].intensity.sum(), raw.loc[accs[0]].sum())   # own row kept, not doubled
        assert (X.intensity >= 0).all() and X.intensity.notna().all()
        _, Xf = io.load(dict(ds, id_groups={"mode": "first"}), {}, lambda x: x)
        assert "Q88888" in set(Xf.protein) and "Q77777" not in set(Xf.protein)
        res = pipeline.run(dict(p, datasets=[ds]))
        assert len(_prot(res, "wide")) > 0
        # two replicate runs in one matrix: a protein absent from run b is dropped for b only (no all-zero profile)
        R = pd.concat([M.add_prefix("a."), M.add_prefix("b.")], axis=1)
        R.loc[accs[2], [c for c in R.columns if c.startswith("b.")]] = np.nan
        f2 = os.path.join(d, "wide_reps.tsv"); R.to_csv(f2, sep="\t", na_rep="")
        _, Y = io.load(dict(ds, path=f2, value_regex=r"^(?P<replicate>\w)\.(?P<fraction>\d+)$"), {}, lambda x: x)
        assert set(Y[Y.protein == accs[2]].replicate) == {"a"} and set(Y[Y.protein == accs[3]].replicate) == {"a", "b"}


def test_calibration_from_and_combine():
    """dataset B reuses dataset A's calibration (same column); conditions combine to one estimate per protein"""
    p, _, raw = _base()
    with tempfile.TemporaryDirectory() as d:
        dss = []
        for sid, scale, drop in (("A", 1.0, 0), ("B", 0.5, 5)):
            M = raw.iloc[drop:].replace(0, np.nan) * scale
            M.columns = [str(i + 1) for i in range(M.shape[1])]
            M.index.name = "ID"
            f = os.path.join(d, f"{sid}.tsv")
            M.to_csv(f, sep="\t", na_rep="")
            dss.append(_ds(p, source_id=sid, condition=sid, format="protein_matrix", path=f, id_col="ID",
                           value_regex=r"^(?P<fraction>\d+)$"))
        dss[1]["calibration_from"] = "A"
        res = pipeline.run(dict(p, datasets=dss))
        cal = {ds["source_id"]: list(ds["runs"].values())[0]["calibration"] for ds in res["_diagnostics"]["datasets"]}
        assert cal["B"]["used"] == "from:A" and cal["B"]["icpt"] == cal["A"]["icpt"] and cal["B"]["slope"] == cal["A"]["slope"]
        assert "anchors" in cal["B"]["checks"]                   # B's own anchor fit is still reported
        a, b = _prot(res, "A"), _prot(res, "B")
        sh = a.index.intersection(b.index)
        assert np.allclose(a[sh], b[sh], equal_nan=True)        # same profiles (scaled) + same calibration
        sh = sh[a[sh].notna().values]                           # measurable (k x monomer inside the calibrated range)
        out = combine.combine_tables({k: v for k, v in res.items() if not k.startswith("_")}, "cell_combined")
        C = out["secms_protein"]
        C = C[(C.source_id == "cell_combined") & (C.threshold_k == 2)].set_index("uniprot")
        assert set(C.index) == set(a.index) | set(b.index)     # unmeasured keys keep a NaN row
        um = C[C.assembled_frac.isna()]
        assert len(um) and (um.confidence == "unmeasured").all() and (um.n_conditions == 0).all()
        assert set(um.index) == set(a.index[a.isna()]) | set(b.index[b.isna()])
        assert (C.loc[sh, "n_conditions"] == 2).all() and np.allclose(C.loc[sh, "assembled_frac"], a[sh])
        only = a.index.difference(b.index)
        only = only[a[only].notna().values]
        assert len(only) >= 1 and (C.loc[only, "n_conditions"] == 1).all() and C.loc[only, "assembled_frac_sd"].isna().all()
        Cm = C[C.assembled_frac.notna()]
        assert (Cm.assembled_frac_lo <= Cm.assembled_frac + 1e-12).all() and (Cm.assembled_frac_hi >= Cm.assembled_frac - 1e-12).all()
        # a disagreeing pair: the mean and the SD over conditions
        P = res["secms_protein"].copy()
        i = P.index[(P.source_id == "B") & (P.uniprot == sh[0]) & (P.threshold_k == 2)]
        P.loc[i, "assembled_frac"] = P.loc[i, "assembled_frac"] * 0 + 0.0
        m = combine.combine_protein(P, "x").set_index(["uniprot", "threshold_k"]).loc[(sh[0], 2)]
        assert np.isclose(m.assembled_frac, a[sh[0]] / 2) and np.isclose(m.assembled_frac_sd, a[sh[0]] / np.sqrt(2))
        # SecView-compatible selection keys
        assert set(out["secms_pairs"].condition) >= {"combined"} if len(out["secms_pairs"]) else True


def test_lab_dataset_configs_parse():
    """the lab 2026 run configs (paths only; the data themselves are not read here)"""
    d = os.path.join(REPO, "toolkits", "secms", "datasets")
    for f, sp, n in (("lab2026_jurkat.yaml", "human", 2), ("lab2026_yeast.yaml", "yeast", 2),
                     ("final2026_jurkat.yaml", "human", 3), ("final2026_yeast.yaml", "yeast", 3),
                     ("update2026_yeast.yaml", "yeast", 4)):
        c = yaml.safe_load(open(os.path.join(d, f)))
        assert c["species"] == sp and c["preset"] == "default" and not c.get("overrides") and len(c["datasets"]) == n
        pipeline.load_species(c["species"])
        ids = [x["source_id"] for x in c["datasets"]]
        for i, x in enumerate(c["datasets"]):
            assert x["format"] in ("protein_long", "protein_matrix", "openswath") and x["n_fractions"] in (27, 54, 60, 70)
            if x.get("calibration_from"):
                assert x["calibration_from"] in ids[:i]


def test_anchor_sets():
    """anchor-set versions: default = human v2 / yeast v3 (v1 minus the 4 consistent outliers per species; yeast v3 also
    minus PDC1, EFT1, review 2026-11-19), legacy_v7 = v1 (the v7 file); v2 files unchanged"""
    for sp, out, ver in (("human", {"PFN1", "HSP90AB1", "YWHAZ", "HSPA8"}, "v2"),
                         ("yeast", {"COF1", "HSP82", "SSA1", "BMH1", "PDC1", "EFT1"}, "v3")):
        S = pipeline.load_species(sp)
        a1 = pipeline._anchors(S, pipeline.params_of({"preset": "legacy_v7"}))
        a2 = pipeline._anchors(S, pipeline.params_of({"preset": "default"}))
        assert pipeline.anchor_file(S, pipeline.params_of({"preset": "default"})) == S["anchor_sets"][ver]
        v2 = pipeline._anchors(S, {"calibration": {"anchor_set": "v2"}})
        assert len(v2) == len(a1) - 4                                             # v2 file unchanged
        assert pipeline.anchor_file(S, pipeline.params_of({"preset": "legacy_v7"})) == S["anchors"]   # v1 file unchanged
        assert set(a1.gene) - set(a2.gene) == out and set(a2.gene) < set(a1.gene)
        m = a1.set_index("gene").drop(sorted(out)).reset_index()
        assert m.fillna(-1).equals(a2.fillna(-1))                                 # other rows identical
    try:
        pipeline.anchor_file(pipeline.load_species("human"), {"calibration": {"anchor_set": "v9"}})
        raise AssertionError("unknown anchor set accepted")
    except ValueError:
        pass


def test_combine_weights_and_fill():
    """spec: weighted mean over primary datasets (renormalised over those present), reliability-weighted SD, fill
    datasets only where no primary dataset scores the protein (confidence low)"""
    rows = []
    for sid, vals in (("A", {"p1": 0.2, "p2": 0.4}), ("B", {"p1": 0.8, "p3": 0.6}), ("F", {"p1": 0.0, "p3": 0.1, "p4": 0.9})):
        for u, v in vals.items():
            rows.append(dict(uniprot=u, threshold_k=2, source_id=sid, condition=sid, replicate="1", assembled_frac=v,
                             assembled_frac_lo=v - 0.05, assembled_frac_hi=v + 0.05, void_frac=0.0, assembled_frac_peaks=v,
                             apex_app_kda=100.0, apex_fraction=10, n_fractions=60, n_peaks=1))
    P = pd.DataFrame(rows)
    spec = {"datasets": [{"source_id": "A", "weight": 1.0}, {"source_id": "B", "weight": 3.0},
                         {"source_id": "F", "weight": 1.0, "role": "fill"}], "discordant_range": 0.3}
    M = combine.combine_protein(P, "x", spec=spec).set_index("uniprot")
    assert np.isclose(M.loc["p1", "assembled_frac"], (0.2 + 3 * 0.8) / 4) and M.loc["p1", "n_conditions"] == 2
    v1, v2 = 4.0, 10.0
    sd = np.sqrt((1 * (0.2 - 0.65) ** 2 + 3 * (0.8 - 0.65) ** 2) / (v1 - v2 / v1))
    assert np.isclose(M.loc["p1", "assembled_frac_sd"], sd) and M.loc["p1", "confidence"] == "discordant"
    assert np.isclose(M.loc["p2", "assembled_frac"], 0.4) and M.loc["p2", "confidence"] == "single"
    assert np.isclose(M.loc["p3", "assembled_frac"], 0.6) and "F|F" not in M.loc["p3", "conditions"]
    assert np.isclose(M.loc["p4", "assembled_frac"], 0.9) and M.loc["p4", "confidence"] == "low" and M.loc["p4", "fill_only"]
    E = combine.combine_protein(P[P.source_id != "F"], "x").set_index("uniprot")   # no spec = equal-weight mean
    assert np.isclose(E.loc["p1", "assembled_frac"], 0.5) and np.isclose(E.loc["p1", "assembled_frac_sd"], np.std([0.2, 0.8], ddof=1))



def test_combine_replicates_split_weight():
    """spec replicates: split_weight -> a dataset's weight is shared by its replicate runs scoring the protein (a dataset
    with 4 replicates counts as much as a single-run dataset of the same weight); without it each run gets the weight"""
    rows = []
    def r(sid, rep, u, v):
        return dict(uniprot=u, threshold_k=2, source_id=sid, condition=sid, replicate=rep, assembled_frac=v,
                    assembled_frac_lo=v, assembled_frac_hi=v, void_frac=0.0, assembled_frac_peaks=v, apex_app_kda=100.0,
                    apex_fraction=10, n_fractions=27, n_peaks=1)
    rows += [r("A", "1", "p1", 0.0), r("A", "1", "p2", 0.0)]
    rows += [r("S", rep, "p1", 1.0) for rep in ("a", "b", "c", "d")] + [r("S", rep, "p2", 1.0) for rep in ("a", "b")]
    P = pd.DataFrame(rows)
    split = {"datasets": [{"source_id": "A"}, {"source_id": "S", "replicates": "split_weight"}]}
    M = combine.combine_protein(P, "x", spec=split).set_index("uniprot")
    assert np.allclose(M.assembled_frac, 0.5) and (M.n_conditions == 2).all() and np.allclose(M.weight_sum, 2.0)
    N = combine.combine_protein(P, "x", spec={"datasets": [{"source_id": "A"}, {"source_id": "S"}]}).set_index("uniprot")
    assert np.isclose(N.loc["p1", "assembled_frac"], 0.8) and np.isclose(N.loc["p2", "assembled_frac"], 2 / 3)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "ok")
