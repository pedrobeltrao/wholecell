"""Core localisation framework (wholecell/inventory/localisation): combiner on synthetic evidence (precedence, blend, split
rules, membrane restriction, new-from-imaging, N/C models; the ordered-annotation steps of legacy_x5: ordered list,
annotation overrides, organelle fractionation, rank weights, imaging groups), override tables (checks, ops, grouping),
presets and config errors, the Jurkat and yeast preset wiring, and a testcell run through the inventory v0 step. The
full yeast check against the port build s0 is scripts/dev/localisation_x5_check.py (needs the yeast data).

    PYTHONHASHSEED=0 python tests/test_localisation.py        (or pytest tests/test_localisation.py)
"""
import json
import math
import os
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.config import load as load_config  # noqa: E402
from wholecell.inventory import localisation as L  # noqa: E402
from wholecell.inventory import stage  # noqa: E402
from wholecell.inventory.localisation import combiner as C  # noqa: E402
from wholecell.inventory.localisation import evidence as EV  # noqa: E402
from wholecell.inventory.localisation import overrides as OV  # noqa: E402

NAN = float("nan")


def _raises(fn, exc, *needles):
    try:
        fn()
    except exc as e:
        for n in needles:
            assert n in str(e), (n, str(e))
        return str(e)
    raise AssertionError(f"no {exc.__name__} (expected {needles})")


def _ev(rows):
    return EV.to_frame([(g, "", s, c, k, sc, rel, 1, "") for g, s, c, k, sc, rel in rows])


def _feat(genes, tm=()):
    return pd.DataFrame({"uniprot": "", "transmembrane": [g in tm for g in genes]}, index=pd.Index(genes, name="gene"))


def _opts(**kw):
    o = {"steps": ["drop_unreliable", "split", "combine"], "within_source": "max", "drop_unreliable": {}, "split": "score",
         "combine": {"mode": "blend", "classes": ["imaging_other", "annotation_curated"],
                     "weights": {"imaging_other": 0.5, "annotation_curated": 0.5}, "key_order": "first_seen",
                     "fallback": {"compartment": "unassigned", "label": "none"}}}
    o.update(kw)
    return o


E1 = [("A", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("A", "hpa", "imaging_other", "cytosol", 0.3, 1.0),
      ("A", "up", "annotation_curated", "cytosol", 1.0, NAN),
      ("B", "up", "annotation_curated", "er", 1.0, NAN), ("B", "up", "annotation_curated", "golgi", 1.0, NAN),
      ("C", "hpa", "imaging_other", "mitochondrion", 1.0, 0.0), ("C", "up", "annotation_curated", "cytosol", 1.0, NAN)]


def test_blend_split_precedence():
    S = C.run(_ev(E1), _feat(["A", "B", "C", "D"]), ["A", "B", "C", "D"], _opts())
    a = S["A"].w                                  # HPA nucleus 1 / 0.3 -> 1/1.3, 0.3/1.3; UniProt cytosol 1; 50:50
    assert abs(a["nucleus"] - 0.5 / 1.3) < 1e-12 and abs(a["cytosol"] - (0.5 * 0.3 / 1.3 + 0.5)) < 1e-12
    assert S["A"].label == ["hpa", "up"]
    assert S["B"].w == {"er": 0.5, "golgi": 0.5} and S["B"].label == ["up"]          # equal scores: equal split
    assert S["C"].w == {"cytosol": 1.0} and S["C"].label == ["up"]                     # reliability 0 dropped
    assert S["D"].w == {"unassigned": 1.0} and S["D"].label == ["none"]               # no evidence: fallback
    o = _opts(combine=dict(_opts()["combine"], mode="precedence"))
    S = C.run(_ev(E1), _feat(["A"]), ["A"], o)
    assert abs(S["A"].w["nucleus"] - 1 / 1.3) < 1e-12 and S["A"].label == ["hpa"]     # imaging_other first in classes
    o = _opts(combine=dict(_opts()["combine"], mode="precedence", classes=["annotation_curated", "imaging_other"]))
    assert C.run(_ev(E1), _feat(["A"]), ["A"], o)["A"].w == {"cytosol": 1.0}
    # split rules on one source with scores 1, 0.6, 0.3
    q = {"x": 1.0, "y": 0.6, "z": 0.3}
    assert C._split(q, "rank") == {"x": 4 / 7, "y": 2 / 7, "z": 1 / 7}
    assert C._split(q, "equal") == {"x": 1 / 3, "y": 1 / 3, "z": 1 / 3}
    assert abs(C._split(q, "score")["y"] - 0.6 / 1.9) < 1e-15
    assert C._split({"x": 0.5, "y": 0.5}, "score") == {"x": 0.5, "y": 0.5}           # ties -> equal


def test_membrane_and_new_from_imaging():
    rows = [("T", "hpa", "imaging_other", "cytosol", 1.0, 1.0), ("T", "up", "annotation_curated", "er", 1.0, NAN),
            ("T", "up", "annotation_curated", "cytosol", 1.0, NAN),
            ("U", "hpa", "imaging_other", "plasma_membrane", 1.0, 1.0), ("U", "up", "annotation_curated", "cytosol", 1.0, NAN)]
    o = _opts(steps=["restrict_membrane", "split", "combine"],
              restrict_membrane={"keys": ["er", "plasma_membrane"], "if_none": {"imaging_other": "drop", "annotation_curated": "keep"}})
    S = C.run(_ev(rows), _feat(["T", "U"], tm={"T"}), ["T", "U"], o)
    assert S["T"].w == {"er": 1.0} and S["T"].label == ["up"]        # TM: HPA cytosol only -> HPA dropped; UniProt er kept
    assert S["U"].w == {"plasma_membrane": 0.5, "cytosol": 0.5}      # not TM: unrestricted
    rows = [("N", "img", "imaging_other", "golgi", 0.4, NAN), ("N", "img", "imaging_other", "er", 0.6, NAN),
            ("N", "img", "imaging_other", "nucleus", 0.9, NAN), ("N", "up", "annotation_curated", "cytosol", 1.0, NAN)]
    o = _opts(steps=["split", "combine", "new_from_imaging"],
              new_from_imaging={"imaging_classes": ["imaging_other"], "min_score": 0.5, "nc_keys": ["nucleus", "cytosol"]})
    w = C.run(_ev(rows), _feat(["N"]), ["N"], o)["N"].w
    assert set(w) == {"cytosol", "er"} and abs(sum(w.values()) - 1) < 1e-12   # golgi < 0.5 and nucleus (N/C key) not new
    # membrane_mechanism: a membrane key only one imaging source supports needs a mechanism feature
    rows = [("M", "img", "imaging_other", "plasma_membrane", 1.0, NAN), ("M", "up", "annotation_curated", "cytosol", 1.0, NAN)]
    o = _opts(steps=["split", "combine", "membrane_mechanism"],
              membrane_mechanism={"keys": ["plasma_membrane"], "imaging_classes": ["imaging_other"]})
    assert C.run(_ev(rows), _feat(["M"]), ["M"], o)["M"].w == {"cytosol": 1.0}
    F = _feat(["M"]); F["mechanism"] = True
    assert C.run(_ev(rows), F, ["M"], o)["M"].w == {"plasma_membrane": 0.5, "cytosol": 0.5}


def test_nc_models():
    nc = {"pool": ["nucleus", "cytosol", "cytoskeleton"], "nucleus": "nucleus", "cytoplasm": ["cytosol", "cytoskeleton"],
          "default_cytoplasm": "cytosol", "skip_transmembrane": True}
    rows = [("F", "up", "annotation_curated", "cytosol", 1.0, NAN), ("F", "up", "annotation_curated", "cytoskeleton", 1.0, NAN),
            ("F", "up", "annotation_curated", "er", 1.0, NAN), ("F", "itz", "fractionation", "nucleus", 0.62, NAN),
            ("G", "up", "annotation_curated", "cytosol", 1.0, NAN), ("G", "itz", "fractionation", "nucleus", 0.3, NAN)]
    o = _opts(steps=["split", "combine", "nc_model"],
              nc_model=dict(nc, method="fraction", **{"class": "fractionation", "background": 0.05, "new_min": 0.5}))
    S = C.run(_ev(rows), _feat(["F", "G"]), ["F", "G"], o)
    fn = (0.62 - 0.05) / 0.95                       # background-corrected nuclear share of the soluble pool 2/3
    w = S["F"].w
    assert abs(w["nucleus"] - 2 / 3 * fn) < 1e-12 and abs(w["cytosol"] - 1 / 3 * (1 - fn)) < 1e-12 and w["er"] == 1 / 3
    assert S["F"].label == ["up", "itz"]
    assert S["G"].w == {"cytosol": 1.0}             # no nucleus key and fn < new_min: unchanged
    rows = [("R", "up", "annotation_curated", "cytosol", 1.0, NAN), ("R", "up", "annotation_curated", "nucleus", 1.0, NAN),
            ("R", "cy", "imaging_other", "nucleus", 0.7, NAN), ("R", "cy", "imaging_other", "cytosol", 0.2, NAN)]
    o = _opts(steps=["split", "combine", "nc_model"], combine=dict(_opts()["combine"], classes=["annotation_curated"]),
              nc_model=dict(nc, method="ratio_exp", classes=["imaging_other"], k=5.0, cytosol_evidence="cytosol",
                            volumes={"nucleus": 3.0, "cytoplasm": 30.0}, new_min=0.0))
    w = C.run(_ev(rows), _feat(["R"]), ["R"], o)["R"].w
    r = math.exp(5.0 * (0.7 - 0.2))
    assert abs(w["nucleus"] - r * 3 / (r * 3 + 30)) < 1e-12 and abs(sum(w.values()) - 1) < 1e-12


def _species():
    return pd.DataFrame({"ing_id": ["P1", "CPX-1", "P2"], "genes": ["G1", "G1;G2", "G2"],
                         "comp_weights": ['{"cytosol": 0.6, "mitochondrion": 0.3, "nucleus": 0.1}', '{"er": 1.0}',
                                          '{"mitochondrion": 1.0}'],
                         "transmembrane": [True, False, False], "n_tm": [1, 0, 0]})


def _ov(rows):
    return OV.check(pd.DataFrame(rows, columns=OV.OVERRIDE_COLUMNS), "test")


def test_overrides():
    I = _species()
    T = _ov([["s", "x.clic", "gene", "G1", "transmembrane", "set", "false", "soluble", ""],
             ["s", "x.clic", "gene", "G1", "comp_weights", "drop", '{"keys": ["mitochondrion"], "fallback": {"cytosol": 1.0}}', "no MTS", ""],
             ["s", "x.mito", "gene", "G2", "comp_weights", "drop", '{"keys": ["mitochondrion"], "fallback": {"cytosol": 1.0}}', "no MTS", ""],
             ["s", "x.pm", "species", "CPX-1", "comp_weights", "set_keep_other",
              '{"shares": {"plasma_membrane": 1.0}, "keep": 0.2, "exclude": ["cytosol"]}', "surface", "ref"]])
    assert [len(g) for g in OV.groups(T)] == [2, 1, 1]
    log = []

    def find(kt, key):
        if kt == "species":
            return [i for i in I.index if I.at[i, "ing_id"] == key][:1]
        return [i for i in I.index if set(I.at[i, "genes"].split(";")) == {key}]

    def get(i, f):
        return json.loads(I.at[i, f]) if f == "comp_weights" else I.at[i, f]

    def put(i, r, new):
        log.append((I.at[i, "ing_id"], r.rule_id, r.field))
        I.at[i, r.field] = json.dumps(new) if r.field == "comp_weights" else new

    assert OV.apply(T, find, get, put) == 4
    assert log == [("P1", "x.clic", "transmembrane"), ("P1", "x.clic", "comp_weights"), ("P2", "x.mito", "comp_weights"),
                   ("CPX-1", "x.pm", "comp_weights")]
    assert json.loads(I.at[0, "comp_weights"]) == {"cytosol": 0.6, "nucleus": 0.1} and not I.at[0, "transmembrane"]
    assert json.loads(I.at[2, "comp_weights"]) == {"cytosol": 1.0}                   # nothing left: fallback
    w = json.loads(I.at[1, "comp_weights"])
    assert abs(w["plasma_membrane"] - 0.8) < 1e-15 and abs(w["er"] - 0.2) < 1e-15
    bad = lambda **kw: [dict(zip(OV.OVERRIDE_COLUMNS, ["s", "x", "gene", "G", "comp_weights", "set", '{"er": 1}', "why", ""]), **kw)]
    chk = lambda rows: OV.check(pd.DataFrame(rows), "t.tsv")
    _raises(lambda: chk(bad(key_type="genes")), OV.OverrideError, "key_type 'genes'", "did you mean 'gene'")
    _raises(lambda: chk(bad(field="comp_weight")), OV.OverrideError, "field 'comp_weight'", "did you mean 'comp_weights'")
    _raises(lambda: chk(bad(op="drp")), OV.OverrideError, "op 'drp'", "did you mean 'drop'")
    _raises(lambda: chk(bad(value="{er: 1}")), OV.OverrideError, "value is not JSON")
    _raises(lambda: chk(bad(reason=" ")), OV.OverrideError, "reason is empty")
    _raises(lambda: OV.check(pd.DataFrame({"rule_id": []}), "t.tsv"), OV.OverrideError, "missing columns")


def test_presets_and_config_errors():
    assert {"legacy_v7", "legacy_x5"} <= set(L.presets())
    o = L.options({"preset": "legacy_v7"})
    assert o["combiner"]["steps"] == ["drop_unreliable", "restrict_membrane", "split", "combine", "round", "nc_model", "round"]
    kinds = {s["id"]: s["kind"] for steps in o["rules"].values() for s in steps}
    assert kinds["v5.P_nuclear_pool"] == kinds["v6.R1_nc"] == "frozen_overrides"
    assert all(s.get("reason") for steps in o["rules"].values() for s in steps if s["kind"] == "legacy")
    LE = L.LocalisationConfigError
    _raises(lambda: L.options({"preset": "legacy_v8"}), LE, "did you mean 'legacy_v7'")
    _raises(lambda: L.options({}), LE, "preset is required")
    _raises(lambda: L.options({"preset": "legacy_v7", "overrides": {"combiner": {"splt": "rank"}}}), LE,
            "overrides.combiner.splt", "did you mean 'split'")
    _raises(lambda: L.options({"preset": "legacy_v7", "overrides": {"combiner": {"split": "rnak"}}}), LE,
            "combiner.split = 'rnak'", "did you mean 'rank'")
    _raises(lambda: L.options({"preset": "legacy_v7", "overrides": {"combiner": {"steps": ["split", "combin"]}}}), LE,
            "unknown step 'combin'", "did you mean 'combine'")
    _raises(lambda: L.options({"preset": "legacy_v7", "overrides": {"combiner": {"combine": {"classes": ["imaging"]}}}}), LE,
            "'imaging' is not an evidence class")


def test_jurkat_preset_and_overrides():
    """Jurkat: core adapter + preset legacy_v7; the override TSV parses, every row's rule id is an `overrides` step, and
    every `overrides` step has rows; the rule-set code's step lists equal the preset's"""
    cfg = load_config("jurkat_v7_port")          # B5: jurkat.yaml runs preset common; the v7 port keeps legacy_v7
    registry.load_plugins(cfg)
    assert cfg["inventory"]["adapters"]["localisation"] == "localisation"
    p = cfg["inventory"]["params"]["localisation"]
    o = L.options(p)
    assert o["preset"] == "legacy_v7" and p["label"] == "jurkat.hpa_uniprot_itzhak_v0"
    T = OV.check(EV.read_curated(p["data"]["overrides"]), "jurkat")
    ov_steps = {(rs, s["id"]) for rs, steps in o["rules"].items() for s in steps if s["kind"] == "overrides"}
    assert set(zip(T.rule_set, T.rule_id)) == ov_steps
    from plugins.jurkat.inventory import rules_v4, rules_v5, rules_v6, rules_v7
    for rs, mod in (("jurkat.v4", rules_v4), ("jurkat.v5", rules_v5), ("jurkat.v6", rules_v6), ("jurkat.v7", rules_v7)):
        got = [(s["id"], s["kind"]) for s in o["rules"][rs]]
        assert [g[0] for g in got] == [c[0] for c in mod.STEPS], rs
        assert all(g[1] in c[1] if isinstance(c[1], tuple) else g[1] == c[1] for g, c in zip(got, mod.STEPS)), rs
    v0 = stage.steps(cfg)["v0"]
    assert "inventory/localisation/evidence.tsv" in v0["outputs"]
    assert set(v0["frozen"]) >= {"inv_uniprot_full", "inv_hpa_subcellular", "inv_itzhak2016"}
    for k in ("term_map", "rules"):
        for sp in p["sources"]:
            if k in sp:
                assert len(EV.read_curated(sp[k])) > 0 and "reason" in EV.read_curated(sp[k]).columns


def test_testcell_v0():
    """the testcell's synthetic table as evidence (one annotation class, score split) through the inventory v0 step gives
    the same weights as the generic localisation_table adapter"""
    loc = {"preset": "legacy_v7", "label": "testcell_framework",
           "overrides": {"combiner": {"steps": ["split", "combine"],
                                      "combine": {"classes": ["annotation_curated"], "key_order": "first_seen",
                                                  "fallback": {"compartment": "unassigned", "label": "none"}}}},
           "genes": {"table": "inv_testcell_localisation"},
           "sources": [{"source": "table", "adapter": "localisation_table", "evidence_class": "annotation_curated",
                        "input": "inv_testcell_localisation"}]}
    ov = {"inventory": {"adapters": {"localisation": "localisation"}, "params": {"localisation": loc}}}
    cfg = load_config("testcell", overrides=ov)
    registry.load_plugins(cfg)
    ref = load_config("testcell")
    with tempfile.TemporaryDirectory() as d, tempfile.TemporaryDirectory() as d0:
        stage.run(cfg, d, step="v0")
        stage.run(ref, d0, step="v0")
        A = pd.read_csv(os.path.join(d, "inventory", "v0", "inventory_v0.tsv"), sep="\t", index_col="gene")
        B = pd.read_csv(os.path.join(d0, "inventory", "v0", "inventory_v0.tsv"), sep="\t", index_col="gene")
        E = EV.read(os.path.join(d, "inventory", "localisation", "evidence.tsv"))
    assert len(E) > 0 and list(E.columns) == EV.COLUMNS
    w = lambda s: {k: float(v) for k, v in (x.rsplit(":", 1) for x in s.split(";"))}
    assert list(A.index) == list(B.index)
    for g in A.index:
        a, b = w(A.at[g, "compartments"]), w(B.at[g, "compartments"])
        tb = sum(b.values())
        assert set(a) == set(b) and all(abs(a[k] - b[k] / tb) < 1e-12 for k in a), (g, a, b)
    assert (A.transmembrane == B.transmembrane).all()
    assert set(A.loc_source) <= {"table", "none"}


MITO = ["mitochondrion", "mito_IM", "mito_OM", "mito_IMS", "mito_matrix"]
RANK = {"base": 0.5, "empty": {"cytosol": 1.0},
        "secondary": [{"primary": ["cytosol", "nucleoplasm"], "keys": MITO, "factor": 0.1},
                      {"primary": MITO, "keys_not": MITO, "factor": 0.1, "flag": "mito_hc"}],
        "pool": {"keys": ["cytosol", "nucleoplasm"], "volumes": {"cytosol": 3.0, "nucleoplasm": 1.0}}}


def _ev_notes(rows):
    """rows (gene, source, class, compartment, score, reliability, note), file order"""
    return EV.to_frame([(g, "", s, c, k, sc, rel, 1, n) for g, s, c, k, sc, rel, n in rows])


def _state(rows, genes, tm=(), names=None):
    F = _feat(genes, tm)
    F["gene_name"] = [(names or {}).get(g, g) for g in genes]
    return C.build_state(_ev_notes(rows), F, genes)


def _close(a, b, tol=1e-15):
    assert set(a) == set(b) and all(abs(a[k] - b[k]) <= tol for k in a), (a, b)


def test_ordered_annotation_and_overrides():
    up, go, gm = "annotation_curated", "annotation_ht", "annotation_ht"
    R = [("A", "up", up, "cytosol", 1.0, NAN, ""), ("A", "up", up, "spindle_pole", 1.0, NAN, ""),
         ("A", "up", up, "mitochondrion", 1.0, NAN, ""), ("A", "go", go, "mito_IM", 1.0, NAN, ""),
         ("A", "go", go, "cytosol", 1.0, NAN, ""),
         ("B", "up", up, "spindle_pole", 1.0, NAN, ""), ("B", "go", go, "vacuole", 1.0, NAN, ""),
         ("B", "go", go, "ER", 1.0, NAN, ""), ("B", "gom", gm, "plasma_membrane", 1.0, NAN, ""),
         ("C", "up", up, "cytosol", 1.0, NAN, ""), ("C", "gom", gm, "ER", 1.0, NAN, ""),
         ("E", "up", up, "mitochondrion", 1.0, NAN, ""), ("E", "go", go, "mito_IM", 1.0, NAN, ""),
         ("E", "go", go, "mito_matrix", 1.0, NAN, ""),
         ("G", "go", go, "golgi_endosome", 1.0, NAN, ""), ("G", "go", go, "ER", 1.0, NAN, "")]
    S = _state(R, list("ABCDEG"), tm={"B", "C"}, names={"A": "SPC42", "B": "CMD1"})
    o = {"primary": ["up"], "resolve": {"spindle_pole": {"to": "nuclear_envelope", "genes": ["SPC42"], "feature": "gene_name"}},
         "refine": {"generic": "mitochondrion", "sub_keys": MITO[1:], "source": "go"},
         "fallback": {"source": "go", "unranked": True, "membrane_source": "gom", "unless_keys": ["cell_wall"]},
         "tm_soluble": {"keys": ["cytosol", "nucleoplasm", "nucleolus"], "source": "gom"}}
    C.ordered_annotation(S, o)
    assert S["A"].order == ["cytosol", "nuclear_envelope", "mito_IM"] and not S["A"].unranked   # SPB gene; one GO sub-key
    assert S["B"].order == ["plasma_membrane"] and S["B"].unranked and S["B"].label == ["go"]   # spindle dropped -> GO, TM
    assert S["C"].order == ["ER"] and S["C"].label == ["up", "gom"]                            # TM with soluble keys only
    assert S["D"].order == [] and S["D"].label == ["none"]
    assert S["E"].order == ["mitochondrion"]                                                   # two GO sub-keys: no refine
    assert S["G"].order == ["golgi_endosome", "ER"] and S["G"].unranked                        # GO fallback (unranked)
    rows = [dict(rule_set="y", rule_id="y.FIX", key_type="gene", key="D", field="annotation", op="set",
                 value='["plasma_membrane"]', reason="curated", citation="test")]
    C.annotation_overrides(S, {"rows": rows, "rule_set": "y", "rule_ids": ["y.FIX"]})
    assert S["D"].order == ["plasma_membrane"] and S["D"].label == ["override:y.FIX"]
    _raises(lambda: _ov([("y", "y.FIX", "gene", "D", "annotation", "set", '{"plasma_membrane": 1}', "r", "")]),
            OV.OverrideError, "annotation set needs a non-empty JSON list")


def test_mito_fractionation():
    up, fr = "annotation_curated", "fractionation"
    R = [("M1", "up", up, "mitochondrion", 1.0, NAN, ""), ("M1", "up", up, "cytosol", 1.0, NAN, ""),
         ("M1", "vog", fr, "mito_matrix", 1.0, NAN, ""), ("M1", "mg", fr, "mitochondrion", 1.0, 1.0, ""),
         ("M2", "up", up, "mitochondrion", 1.0, NAN, ""), ("M2", "vog", fr, "mito_IM", 1.0, NAN, ""),
         ("M2", "mg", fr, "mitochondrion", 0.0, 0.0, ""), ("M2", "go", "annotation_ht", "golgi_endosome", 1.0, NAN, ""),
         ("M2", "go", "annotation_ht", "mitochondrion", 1.0, NAN, ""),
         ("M3", "up", up, "mitochondrion", 1.0, NAN, ""), ("M3", "up", up, "ER", 1.0, NAN, ""),
         ("M3", "mg", fr, "mito_matrix", 1.0, 1.0, ""),
         ("M4", "up", up, "cytosol", 1.0, NAN, ""), ("M4", "mg", fr, "mitochondrion", 1.0, 1.0, ""),
         ("M5", "up", up, "mitochondrion", 1.0, NAN, ""), ("M5", "mg", fr, "mitochondrion", 0.0, 0.0, "")]
    genes = ["M1", "M2", "M3", "M4", "M5"]
    S = _state(R, genes, tm={"M3"})
    C.ordered_annotation(S, {"primary": ["up"]})
    C.mito_fractionation(S, {"keys": MITO, "generic": "mitochondrion", "measured": "vog", "proteome": "mg", "veto": True,
                             "tm_map": {"mito_matrix": "mito_IM", "mito_IMS": "mito_IM"}, "flag": "mito_hc",
                             "not_mito_fallback": {"sources": [{"source": "gom", "transmembrane_only": True},
                                                               {"source": "go"}], "default": "cytosol"}})
    assert S["M1"].order == ["mito_matrix", "cytosol"] and S["M1"].flags["mito_hc"]     # measured sub-compartment
    assert S["M2"].order == ["golgi_endosome"] and not S["M2"].flags["mito_hc"]         # veto: not HC + depleted
    assert S["M3"].order == ["mito_IM", "ER"] and S["M3"].flags["mito_hc"]              # sub-class, TM -> inner membrane
    assert S["M4"].order == ["cytosol"] and S["M4"].flags["mito_hc"]
    assert S["M5"].order == ["cytosol"]                                                 # nothing left -> default


def test_rank_weights():
    _close(C._rank(["cytosol", "mito_IM", "ER"], RANK), {"cytosol": 1 / 1.3, "mito_IM": 0.05 / 1.3, "ER": 0.25 / 1.3})
    _close(C._rank(["mito_IM", "vacuole"], RANK, {"mito_hc": True}), {"mito_IM": 1 / 1.05, "vacuole": 0.05 / 1.05})
    _close(C._rank(["mito_IM", "vacuole"], RANK), {"mito_IM": 2 / 3, "vacuole": 1 / 3})          # flag off: plain rank
    _close(C._rank(["nucleoplasm", "cytosol"], RANK), {"nucleoplasm": 0.25, "cytosol": 0.75})    # pooled, split by volume
    _close(C._rank(["vacuole", "ER"], RANK, equal=True), {"vacuole": 0.5, "ER": 0.5})
    assert C._rank([], RANK) == {"cytosol": 1.0}
    S = _state([], ["X"])
    S["X"].order = ["ER", "vacuole", "peroxisome"]
    C.rank_weights(S, RANK)
    _close(S["X"].w, {"ER": 4 / 7, "vacuole": 2 / 7, "peroxisome": 1 / 7})


def test_imaging_groups():
    im = "imaging_other"
    R = [("I1", "cy", im, "cytosol", 0.6, NAN, "Cytoplasm"), ("I1", "cy", im, "nucleoplasm", 0.2, NAN, "Nucleus"),
         ("I1", "cy", im, "periphery", 0.2, NAN, "Cell Periphery"), ("I1", "cy", im, "ER", 0.6, NAN, "ER"),
         ("I1", "cy", im, "mitochondrion", 0.7, NAN, "Mitochondria"),
         ("I2", "cy", im, "nucleoplasm", 0.5, NAN, "Nucleus"), ("I2", "cy", im, "cytosol", 0.3, NAN, "Cytoplasm"),
         ("I3", "cy", im, "cytosol", 0.9, NAN, "Cytoplasm"),
         ("I4", "cy", im, "cytosol", 0.9, NAN, "Cytoplasm"),
         ("I5", "cy", im, "cytosol", 0.9, NAN, "Cytoplasm"),
         ("I6", "cy", im, "mitochondrion", 0.6, NAN, "Mitochondria"), ("I6", "cy", im, "cytosol", 0.4, NAN, "Cytoplasm"),
         ("I7", "cy", im, "periphery", 0.8, NAN, "Cell Periphery")]
    genes = ["I1", "I2", "I3", "I4", "I5", "I6", "I7"]
    S = _state(R, genes, tm={"I4", "I5"})
    order = {"I1": ["cytosol"], "I2": ["nucleoplasm", "cytosol"], "I3": ["cell_wall"], "I4": ["ER"], "I5": ["cytosol"],
             "I6": ["cytosol"], "I7": ["plasma_membrane"]}
    for g in genes:
        S[g].order = order[g]
    S["I6"].flags["mito_hc"] = True
    C.rank_weights(S, RANK)
    before = {g: dict(S[g].w) for g in genes}
    o = {"source": "cy", "groups": {"cytosol": "nc", "nucleoplasm": "nc", "nucleolus": "nc", "nuclear_envelope": "ER"},
         "periphery": {"key": "periphery", "membrane": "plasma_membrane", "other": "cytosol"},
         "skip_if_keys": ["cell_wall"], "new": {"min_score": 0.5, "never": ["nc"], "needs_flag": {"mitochondrion": "mito_hc"}},
         "nc": {"group": "nc", "nucleus": "nucleoplasm", "cytosol": "cytosol", "nucleus_term": "Nucleus",
                "cytosol_term": "Cytoplasm", "k": 1.0, "volumes": {"cytosol": 3.0, "nucleoplasm": 1.0}}, "rank": RANK}
    C.imaging_groups(S, o)
    # group shares from summed scores (periphery -> cytosol for a soluble protein); new ER group (0.6 >= 0.5); the
    # mitochondrion group needs the mito_hc flag; within the nc group the annotation decides (cytosol only)
    _close(S["I1"].w, {"cytosol": 1.0 / 1.6, "ER": 0.6 / 1.6})
    en = math.exp(0.2)                                       # r = exp(k (N - C)) = exp(0.5 - 0.3), Vn 1, Vc 3
    _close(S["I2"].w, {"nucleoplasm": en / (en + 3), "cytosol": 3 / (en + 3)}, 1e-15)
    for g in ("I3", "I4", "I5"):                             # cell wall / no shared group / TM with the nc group only
        assert S[g].w == before[g], g
    _close(S["I6"].w, {"cytosol": 0.4, "mitochondrion": 0.6})
    assert S["I7"].w == {"plasma_membrane": 1.0}             # periphery = PM for a protein annotated there


def test_yeast_preset_and_data():
    """legacy_x5 validates; the yeast localisation block names known adapters / classes and its data files parse"""
    import yaml
    o = L.options({"preset": "legacy_x5"})
    assert o["combiner"]["steps"] == ["ordered_annotation", "annotation_overrides", "mito_fractionation", "rank_weights",
                                      "imaging_groups", "round"]
    assert o["combiner"]["imaging_groups"]["nc"]["k"] == 5.0 and o["combiner"]["imaging_groups"]["new"]["min_score"] == 0.5
    frag = yaml.safe_load(open(os.path.join(REPO, "config/cells/yeast.yaml")))      # the localisation block (A9: merged)
    p = frag["inventory"]["params"]["localisation"]
    assert p["preset"] == "common"                   # B5 switch; the A9 reference yeast_x5_core keeps legacy_x5
    assert load_config("yeast_x5_core")["inventory"]["params"]["localisation"]["preset"] == "legacy_x5"
    for sp in p["sources"]:
        assert sp["adapter"] in EV.SOURCES and sp["evidence_class"] in EV.CLASSES, sp["source"]
        for k in ("rules", "term_map"):
            if k in sp:
                assert len(EV.read_curated(sp[k])) > 0
    cyc = [sp for sp in p["sources"] if sp["source"] == "cyclops"][0]
    assert cyc["min_score"] == 0.05
    cb = L._resolve_data(o["combiner"], p)
    assert "SPC42" in cb["ordered_annotation"]["resolve"]["spindle_pole"]["genes"]
    T = OV.check(pd.DataFrame(cb["annotation_overrides"]["rows"]), "yeast overrides")
    X = OV.select(T, "yeast.x5", ["x5.LOC_FIX"])       # legacy_x5 reads only these rows
    assert list(X.rule_id) == ["x5.LOC_FIX"] and json.loads(X.value[0]) == ["plasma_membrane"]
    # preset common: curated comp_weights rows (2026-11-19 yeast B5 audit, approved): Sir2/3/4 stay nuclear
    wo = p["common"]["weight_overrides"]
    assert wo["rule_set"] == "yeast.common" and "b5.nuclear_silencing" in wo["rule_ids"]
    W = OV.select(T, wo["rule_set"], ["b5.nuclear_silencing"])
    assert dict(zip(W.key, W.value.map(json.loads))) == {"YDL042C": {"nucleolus": 1.0}, "YLR442C": {"nucleoplasm": 1.0},
                                                          "YDR227W": {"nucleoplasm": 1.0}}
    assert (W.field == "comp_weights").all() and W.citation.str.len().gt(0).all()
    assert set(frag["frozen_inputs"]) >= {"loc_yeast_vogtle2017", "loc_yeast_morgenstern2017", "loc_yeast_cyclops_wt1"}


# ------------------------------------------------------------------------------------------------ preset common (B5)
def test_common_weight_overrides():
    """step weight_overrides (2026-11-19): comp_weights rows of the selected rule ids replace the combined weights
    (renormalised, label override:<id>); other rule sets / fields untouched; the step sits after nc_model in common"""
    o = L.options({"preset": "common"})
    st = o["combiner"]["steps"]
    assert st.index("weight_overrides") == st.index("nc_model") + 1 and st[-1] == "round"
    assert "weight_overrides" in C.OPTIONAL
    E = _ev_rows([("G1", "a", "annotation_curated", "nucleoplasm", 1.0, None, ""),
                  ("G2", "a", "annotation_curated", "cytosol", 1.0, None, "")])
    F = pd.DataFrame({"uniprot": [None, None], "transmembrane": [False, False]}, index=pd.Index(["G1", "G2"], name="gene"))
    S = C.build_state(E, F, ["G1", "G2"])
    S["G1"].w, S["G1"].label = {"cytosol": 0.88, "nucleoplasm": 0.12}, ["a"]
    S["G2"].w, S["G2"].label = {"cytosol": 1.0}, ["a"]
    rows = [dict(rule_set="y.c", rule_id="y.keep", key_type="gene", key="G1", field="comp_weights", op="set",
                 value='{"nucleoplasm": 2.0}', reason="r", citation="c"),
            dict(rule_set="y.c", rule_id="y.other", key_type="gene", key="G2", field="comp_weights", op="set",
                 value='{"ER": 1.0}', reason="r", citation="c"),
            dict(rule_set="y.x5", rule_id="y.keep", key_type="gene", key="G2", field="annotation", op="set",
                 value='["ER"]', reason="r", citation="c")]
    C.weight_overrides(S, {"rows": rows, "rule_set": "y.c", "rule_ids": ["y.keep"]})
    assert S["G1"].w == {"nucleoplasm": 1.0} and S["G1"].label == ["a", "override:y.keep"]
    assert S["G2"].w == {"cytosol": 1.0}


def _ev_rows(rows):
    """rows: (gene, source, class, key, score, reliability, note)"""
    return EV.to_frame([(g, "", s, c, k, sc, rel, 1, n) for g, s, c, k, sc, rel, n in rows])


def test_common_combine_classes_and_new_from_imaging():
    o = {"steps": ["split", "combine", "new_from_imaging"], "within_source": "max", "split": "score",
         "combine": {"mode": "classes", "exclusive": ["imaging_cell"], "annotation": ["annotation_curated", "annotation_ht"],
                     "imaging": ["imaging_other"], "weights": {"annotation": 0.5, "imaging": 0.5},
                     "groups": {"mito_IM": "mito", "mito_matrix": "mito", "mitochondrion": "mito"},
                     "fallback": {"compartment": "cytosol", "label": "none"}},
         "new_from_imaging": {"imaging_classes": ["imaging_cell", "imaging_other"], "min_score": 0.5, "nc_keys": ["cytosol"],
                              "groups": {"mito_IM": "mito", "mito_matrix": "mito", "mitochondrion": "mito"},
                              "needs_flag": {"mito": "hc"}, "no_support_classes": ["fractionation"]}}
    rows = [("A", "up", "annotation_curated", "mito_matrix", 1.0, NAN, ""), ("A", "up", "annotation_curated", "cytosol", 1.0, NAN, ""),
            ("A", "cy", "imaging_other", "mitochondrion", 0.9, NAN, ""), ("A", "cy", "imaging_other", "cytosol", 0.1, NAN, ""),
            ("B", "go", "annotation_ht", "cytosol", 1.0, NAN, ""), ("B", "up", "annotation_curated", "ER", 1.0, NAN, ""),
            ("C", "go", "annotation_ht", "ER", 1.0, NAN, ""), ("C", "cy", "imaging_other", "mitochondrion", 0.8, NAN, ""),
            ("C", "cy", "imaging_other", "ER", 0.2, NAN, ""),
            ("D", "img", "imaging_cell", "nucleus", 1.0, NAN, ""), ("D", "up", "annotation_curated", "cytosol", 1.0, NAN, "")]
    S = C.run(_ev_rows(rows), _feat(list("ABCDE")), list("ABCDE"), o)
    # A: groups mito 0.5*0.5 + 0.5*0.9 = 0.7 (annotation's key mito_matrix), cytosol 0.5*0.5 + 0.5*0.1 = 0.3
    _close(S["A"].w, {"mito_matrix": 0.7, "cytosol": 0.3}, 1e-12)
    assert S["B"].w == {"ER": 1.0} and S["B"].label == ["up"]          # curated before high-throughput annotation
    # C: imaging-only group mito (0.8 >= 0.5) dropped without the flag; ER kept
    assert S["C"].w == {"ER": 1.0}
    S["C"].flags["hc"] = True
    assert S["D"].w == {"nucleus": 1.0} and S["D"].label == ["img"]    # same-cell imaging is exclusive
    assert S["E"].w == {"cytosol": 1.0} and S["E"].label == ["none"]
    # with the flag the imaging-only group stays
    from wholecell.inventory.localisation.combiner import build_state, STEPS
    St = build_state(_ev_rows(rows), _feat(["C"]), ["C"])
    St["C"].flags["hc"] = True
    for s in o["steps"]:
        STEPS[s](St, o.get(s))
    _close(St["C"].w, {"mitochondrion": 0.4, "ER": 0.6}, 1e-12)


class _B:
    def __init__(self, d):
        self.bdir = d

    def say(self, *a):
        pass


def test_common_nc_ratio():
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "geometry"))
    json.dump({"region_volumes_mc": {"volume_um3": {"nucleoplasm": 2.0, "cytosol": 8.0}}}, open(os.path.join(d, "geometry", "measures.json"), "w"))
    cal = {"method": "category", "source": "hpa", "main_min": 1.0, "skip_main_keys": ["er"],
           "ratios": {"N": 10.0, "Nc": 3.0, "NC": 1.0, "Cn": 0.3, "C": 0.02},
           "tie_break": {"categories": ["NC"], "class": "fractionation", "background": 0.05, "clip": [1 / 3, 3.0]},
           "exclude_genes": ["X"]}
    nc = {"method": "ratio", "pool": {"nucleus": "nucleus", "cytoplasm": ["cytosol", "cytoskeleton"], "default_cytoplasm": "cytosol"},
          "calibration": cal, "volumes": {"nucleus": "nucleoplasm", "cytoplasm": "cytosol", "measures": "geometry"},
          "add_min": 0.3, "new_min": 0.5, "skip_transmembrane": True}
    o = _opts(steps=["split", "combine", "nc_model"], nc_model=nc)
    o["combine"]["classes"] = ["annotation_curated"]
    rows = [("N", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("N", "up", "annotation_curated", "cytosol", 1.0, NAN),
            ("Q", "hpa", "imaging_other", "cytosol", 1.0, 1.0), ("Q", "hpa", "imaging_other", "nucleus", 0.3, 1.0),
            ("Q", "up", "annotation_curated", "cytosol", 1.0, NAN), ("Q", "up", "annotation_curated", "er", 1.0, NAN),
            ("M", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("M", "hpa", "imaging_other", "cytosol", 1.0, 1.0),
            ("M", "itz", "fractionation", "nucleus", 0.62, NAN), ("M", "up", "annotation_curated", "nucleus", 1.0, NAN),
            ("M", "up", "annotation_curated", "cytosol", 1.0, NAN),
            ("E", "up", "annotation_curated", "nucleus", 1.0, NAN), ("E", "up", "annotation_curated", "cytosol", 1.0, NAN),
            ("X", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("X", "up", "annotation_curated", "cytosol", 1.0, NAN),
            ("O", "hpa", "imaging_other", "er", 1.0, 1.0), ("O", "hpa", "imaging_other", "nucleus", 1.0, 1.0),
            ("O", "up", "annotation_curated", "cytosol", 1.0, NAN), ("O", "up", "annotation_curated", "nucleus", 1.0, NAN)]
    S = C.run(_ev(rows), _feat(list("NQMEXO")), list("NQMEXO"), o, _B(d))
    fn = lambda r: r * 2 / (r * 2 + 8)
    _close(S["N"].w, {"nucleus": fn(10), "cytosol": 1 - fn(10)}, 1e-12)       # category N; nucleus side named by imaging
    assert S["N"].label[-1] == "nc:hpa"
    _close(S["Q"].w, {"nucleus": 0.5 * fn(0.3), "cytosol": 0.5 * (1 - fn(0.3)), "er": 0.5}, 1e-12)   # Cn (additional 0.3)
    f = (0.62 - 0.05) / 0.95
    r = min(max((f / 2) / ((1 - f) / 8), 1 / 3), 3.0)                       # NC: Itzhak tie-breaker
    _close(S["M"].w, {"nucleus": fn(r), "cytosol": 1 - fn(r)}, 1e-12)
    _close(S["E"].w, {"nucleus": fn(1), "cytosol": 1 - fn(1)}, 1e-12)        # no imaging: r = 1
    assert S["E"].label[-1] == "nc:equal"
    assert S["X"].w == {"cytosol": 1.0}                                      # hold-out gene: no imaging calibration
    _close(S["O"].w, {"nucleus": fn(1), "cytosol": 1 - fn(1)}, 1e-12)      # organelle main: no category, r = 1
    # exp calibration (yeast): r = exp(k (N - C)) of the terms
    cal2 = {"method": "exp", "source": "cy", "k": 5.0, "nucleus_term": "Nucleus", "cytosol_term": "Cytoplasm"}
    o2 = _opts(steps=["split", "combine", "nc_model"], nc_model=dict(nc, calibration=cal2,
                                                                     pool={"nucleus": "nucleoplasm", "cytoplasm": ["cytosol"], "default_cytoplasm": "cytosol"}))
    o2["combine"]["classes"] = ["annotation_curated"]
    rows2 = [("Y", "cy", "imaging_other", "nucleoplasm", 0.7, NAN, "Nucleus"), ("Y", "cy", "imaging_other", "cytosol", 0.2, NAN, "Cytoplasm"),
             ("Y", "up", "annotation_curated", "cytosol", 1.0, NAN, ""), ("Y", "up", "annotation_curated", "nucleoplasm", 1.0, NAN, "")]
    w = C.run(_ev_rows(rows2), _feat(["Y"]), ["Y"], o2, _B(d))["Y"].w
    rr = math.exp(5 * 0.5)
    _close(w, {"nucleoplasm": fn(rr), "cytosol": 1 - fn(rr)}, 1e-12)


def test_common_annotation_evidence_and_resolve_imaging():
    rows = [("G", "uniprot", "annotation_curated", "ER", 1.0, NAN, ""), ("G", "go", "annotation_ht", "cytosol", 1.0, NAN, ""),
            ("G", "cy", "imaging_other", "periphery", 0.6, NAN, "Cell Periphery"),
            ("H", "go", "annotation_ht", "cell_wall", 1.0, NAN, ""), ("H", "cy", "imaging_other", "cytosol", 0.9, NAN, "Cytoplasm")]
    S = C.build_state(_ev_rows(rows), _feat(["G", "H"], tm={"G"}), ["G", "H"])
    S["G"].order, S["G"].label = ["ER"], ["uniprot"]
    S["H"].order, S["H"].label = ["cell_wall"], ["go"]
    C.annotation_evidence(S, {"source": "annotation", "classes": {"uniprot": "annotation_curated", "go": "annotation_ht"},
                              "consume": ["uniprot", "go"]})
    assert S["G"].src["annotation"] == {"ER": 1.0} and S["G"].cls["annotation"] == "annotation_curated"
    assert S["H"].cls["annotation"] == "annotation_ht" and "go" not in S["H"].src
    C.resolve_imaging(S, {"imaging_classes": ["imaging_other"], "skip_if_keys": ["cell_wall"],
                          "pseudo": {"periphery": {"membrane": "plasma_membrane", "other": "cytosol", "annotated": ["plasma_membrane"]}}})
    assert S["G"].src["cy"] == {"plasma_membrane": 0.6} and S["G"].rows["cy"][0].compartment == "plasma_membrane"   # TM
    assert "cy" not in S["H"].src                                            # cell-wall protein: no imaging


def test_common_preset_cells():
    """preset common validates; both cells' common data resolve; rule kinds; the switch script edits only what it names"""
    o = L.options({"preset": "common"})
    kinds = {s["id"]: s["kind"] for steps in o["rules"].values() for s in steps}
    assert kinds["v5.P_nuclear_pool"] == kinds["v6.R1_nc"] == "rule"
    assert o["combiner"]["split"] == "score" and o["combiner"]["combine"]["mode"] == "classes"
    for cell in ("jurkat", "yeast"):
        cfg = load_config(cell, overrides={"inventory": {"params": {"localisation": {"preset": "common"}}}})
        p = cfg["inventory"]["params"]["localisation"]
        cb = L._resolve_data(L.options(p)["combiner"], p)
        assert cb["nc_model"]["calibration"]["method"] == {"jurkat": "category", "yeast": "exp"}[cell]
        assert cb["combine"]["fallback"]["compartment"] == {"jurkat": "unassigned", "yeast": "cytosol"}[cell]
        assert len(cb["restrict_membrane"]["keys"]) >= 8
        assert cb["nc_model"]["bound"] == {"max_share": 0.9, "min_sources": 2, "calibrations": ["category"]}   # revision 2026-11-18
        assert len(cb["nc_model"]["chromatin"]["regex"]) >= 2
        if cell == "jurkat":
            assert cb["ordered_annotation"] is None and cb["combine"]["groups"] is None
            assert cb["imaging_background"]["sources"] == ["hpa"] and 0 <= cb["imaging_background"]["b"] < 1
            assert cb["nc_model"]["prior"]["source"] == "itzhak" and cb["new_from_imaging"]["background_sources"] == "all"
        else:
            assert cb["imaging_background"] is None and cb["nc_model"]["prior"] is None
            assert cb["nc_model"]["chromatin"]["feature"] == "gene_name"
            assert "SPC42" in cb["ordered_annotation"]["resolve"]["spindle_pole"]["genes"]
            assert cb["combine"]["groups"]["mito_IM"] == "mitochondrion"
        # B5 switch applied to the cell configs (step 2 combined build); the reference configs keep the legacy presets
        assert load_config(cell)["inventory"]["params"]["localisation"]["preset"] == "common"
        ref = {"jurkat": "jurkat_v7_port", "yeast": "yeast_x5_core"}[cell]
        assert load_config(ref)["inventory"]["params"]["localisation"]["preset"] == {"jurkat": "legacy_v7", "yeast": "legacy_x5"}[cell]
    import shutil
    import subprocess
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "config", "cells"))
    for c in ("jurkat", "yeast"):                     # the switch script on the configs before the switch (dev repo commit a209ad6)
        # copies of `git show a209ad6:config/cells/<c>.yaml` (release repo has a fresh history; AFDB tar path -> {afdb_root})
        txt = open(os.path.join(REPO, "tests", "data", "localisation_switch", f"{c}_a209ad6.yaml")).read()
        open(os.path.join(d, "config", "cells", f"{c}.yaml"), "w").write(txt)
    pre = {}
    import yaml
    for c in ("jurkat", "yeast"):
        pre[c] = yaml.safe_load(open(os.path.join(d, "config", "cells", f"{c}.yaml")))
    subprocess.run([sys.executable, os.path.join(REPO, "scripts", "dev", "localisation_switch.py"), "--repo", d, "jurkat", "yeast"],
                   check=True, capture_output=True)
    for c in ("jurkat", "yeast"):
        a = pre[c]
        b = yaml.safe_load(open(os.path.join(d, "config", "cells", f"{c}.yaml")))
        la, lb = a["inventory"]["params"]["localisation"], b["inventory"]["params"]["localisation"]
        assert lb["preset"] == "common" and lb["label"] == f"{c}.localisation_common"
        assert {k: v for k, v in la.items() if k not in ("preset", "label")} == {k: v for k, v in lb.items() if k not in ("preset", "label")}
        a["inventory"]["params"]["localisation"] = b["inventory"]["params"]["localisation"] = None
        a["inventory"].pop("localisation_sources", None); b["inventory"].pop("localisation_sources", None)
        assert a == b, c                                    # nothing else changed
    b = yaml.safe_load(open(os.path.join(d, "config", "cells", "jurkat.yaml")))
    tm = [s["term_map"] for s in b["inventory"]["localisation_sources"]]
    assert tm[0] == "plugins/jurkat/data/localisation/hpa_term_map.tsv" and tm[1] == "plugins/jurkat/data/localisation/uniprot_cc_rules.tsv"
    for t in tm[:2]:
        assert os.path.exists(os.path.join(REPO, t))


def test_common_revision_2026_11_18():
    """revision 2026-11-18 (docs/unification_criteria/B5_localisation_revision_2026-11-18.md): HPA background b, pool keys
    from background-corrected imaging, Itzhak prior, bounded mixtures, chromatin-bound pools, nc_model.tsv record"""
    d = tempfile.mkdtemp()
    os.makedirs(os.path.join(d, "geometry"))
    json.dump({"region_volumes_mc": {"volume_um3": {"nucleoplasm": 2.0, "cytosol": 3.0}}}, open(os.path.join(d, "geometry", "measures.json"), "w"))
    fn = lambda r: r * 2 / (r * 2 + 3)
    # imaging_background: obs = 1/1.3 nuclear (main nucleus 1, additional cytosol 0.3) -> true = (obs - b) / (1 - b)
    b = 0.25
    rows = [("A", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("A", "hpa", "imaging_other", "cytosol", 0.3, 1.0),
            ("B", "hpa", "imaging_other", "mitochondrion", 1.0, 1.0), ("B", "hpa", "imaging_other", "nucleus", 0.3, 1.0),
            ("C", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("D", "cy", "imaging_other", "nucleus", 0.5, 1.0),
            ("D", "cy", "imaging_other", "cytosol", 0.5, 1.0)]
    S = C.build_state(_ev(rows), _feat(list("ABCD")), list("ABCD"))
    C.split(S, "score")
    C.imaging_background(S, {"sources": ["hpa"], "b": b, "nucleus_keys": ["nucleus"]})
    obs = 1 / 1.3
    t = (obs - b) / (1 - b)
    _close(S["A"].src["hpa"], {"nucleus": t, "cytosol": 1 - t}, 1e-12)
    assert S["B"].src["hpa"] == {"mitochondrion": 1.0} and S["B"].flags["bg_drop_n"] == {"hpa"}   # 0.23 < b: removed
    assert S["C"].src["hpa"] == {"nucleus": 1.0}                             # obs 1 -> true 1
    assert S["D"].src["cy"] == {"nucleus": 0.5, "cytosol": 0.5} and "bg_sources" not in S["D"].flags   # not an HPA source
    _raises(lambda: C.imaging_background(S, {"sources": ["hpa"], "b": 1.0, "nucleus_keys": ["nucleus"]}), ValueError, "[0, 1)")
    # new_from_imaging background_sources = all: a background-corrected source is blended as in v7; without the step the
    # X5 guards apply (nc keys never new from imaging, min_score 0.5)
    nfi = {"imaging_classes": ["imaging_other"], "min_score": 0.5, "nc_keys": ["nucleus", "cytosol"],
           "no_support_classes": ["fractionation"], "background_sources": "all"}
    o = {"steps": ["split", "imaging_background", "combine", "new_from_imaging"], "within_source": "max", "split": "score",
         "imaging_background": {"sources": ["hpa"], "b": 0.25, "nucleus_keys": ["nucleus"]},
         "combine": {"mode": "classes", "annotation": ["annotation_curated"], "imaging": ["imaging_other"],
                     "weights": {"annotation": 0.5, "imaging": 0.5}, "fallback": {"compartment": "unassigned", "label": "none"}},
         "new_from_imaging": nfi}
    rw = [("W", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("W", "up", "annotation_curated", "endolysosomal", 1.0, NAN),
          ("V", "hpa", "imaging_other", "plasma_membrane", 1.0, 1.0), ("V", "hpa", "imaging_other", "golgi", 0.3, 1.0),
          ("V", "up", "annotation_curated", "plasma_membrane", 1.0, NAN)]
    S = C.run(_ev(rw), _feat(["W", "V"]), ["W", "V"], o)
    assert S["W"].w == {"endolysosomal": 0.5, "nucleus": 0.5}               # HPA-only nucleus kept (background-corrected)
    _close(S["V"].w, {"plasma_membrane": 0.5 + 0.5 / 1.3, "golgi": 0.5 * 0.3 / 1.3}, 1e-12)   # additional term kept (v7)
    o0 = dict(o, steps=["split", "combine", "new_from_imaging"])
    assert C.run(_ev(rw), _feat(["W"]), ["W"], o0)["W"].w == {"endolysosomal": 1.0}
    # nc_target: prior blend in log space, prior alone shrunk to r = 1, bound unless 2 exclusive sources
    pr = {"source": "itz", "background": 0.05, "weight": 1.0, "clip": [0.01, 100.0]}
    oo = {"prior": pr, "bound": {"max_share": 0.9, "min_sources": 2}}
    rf = C.ratio_of_share(0.2, 2, 3, pr["clip"])
    _close({"r": rf}, {"r": (0.2 / 2) / (0.8 / 3)}, 1e-12)
    r, f, ti = C.nc_target(10.0, 0.2, {"n": 1, "c": 0}, oo, 2, 3)
    _close({"r": r}, {"r": math.sqrt(10 * rf)}, 1e-12)
    r, f, ti = C.nc_target(None, 0.2, {"n": 0, "c": 0}, oo, 2, 3)
    _close({"r": r}, {"r": math.sqrt(rf)}, 1e-12)
    assert ti["prior_only"]
    r, f, ti = C.nc_target(0.02, None, {"n": 0, "c": 1}, oo, 2, 3)
    assert abs(f - 0.1) < 1e-12 and ti["bounded"][1] == 0.1                  # single C category: >= 1 - m nuclear
    r, f, ti = C.nc_target(0.02, None, {"n": 0, "c": 2}, oo, 2, 3)
    assert abs(f - fn(0.02)) < 1e-12 and "bounded" not in ti                 # two exclusive sources: unbounded
    r, f, ti = C.nc_target(3.0, 0.5, {"n": 0, "c": 0}, oo, 2, 3, tie_break=True)
    assert r == 3.0 and "r_prior" not in ti                                  # NC tie-break already used the fractionation
    assert C.nc_target(None, None, {}, oo, 2, 3)[0] is None
    oo2 = {"prior": pr, "bound": dict(oo["bound"], calibrations=["category"])}
    r, f, ti = C.nc_target(0.02, None, {"n": 0, "c": 1}, oo2, 2, 3, method="exp")
    assert abs(f - fn(0.02)) < 1e-12 and "bounded" not in ti                 # continuous calibration (yeast): no bound
    # the N/C model end to end: prior, bound (single imaging source adds the missing side), chromatin
    cal = {"method": "category", "source": "hpa", "main_min": 1.0, "skip_main_keys": ["er"],
           "ratios": {"N": 10.0, "Nc": 3.0, "NC": 1.0, "Cn": 0.3, "C": 0.02}}
    nc = {"method": "ratio", "pool": {"nucleus": "nucleus", "cytoplasm": ["cytosol"], "default_cytoplasm": "cytosol"},
          "calibration": cal, "volumes": {"nucleus": "nucleoplasm", "cytoplasm": "cytosol", "measures": "geometry"},
          "add_min": 0.3, "new_min": 0.5, "skip_transmembrane": True, "prior": pr, "bound": oo["bound"],
          "chromatin": {"regex": ["^H2BC\\d+$", "^H1-\\d+$"]}}
    o = _opts(steps=["split", "combine", "nc_model"], nc_model=nc)
    o["combine"]["classes"] = ["imaging_other", "annotation_curated"]
    o["combine"]["mode"] = "precedence"
    rows = [("TKT", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("TKT", "itz", "fractionation", "nucleus", 0.07, NAN),
            ("SOLO", "hpa", "imaging_other", "nucleus", 1.0, 1.0),
            ("CY", "hpa", "imaging_other", "cytosol", 1.0, 1.0),
            ("H2BC12", "hpa", "imaging_other", "nucleus", 1.0, 1.0), ("H2BC12", "hpa", "imaging_other", "cytosol", 0.3, 1.0),
            ("H1-9", "up", "annotation_curated", "cytosol", 1.0, NAN)]
    g = ["TKT", "SOLO", "CY", "H2BC12", "H1-9"]
    S = C.run(_ev(rows), _feat(g), g, o, _B(d))
    ff = (0.07 - 0.05) / 0.95
    rt = math.sqrt(10 * C.ratio_of_share(ff, 2, 3, pr["clip"]))
    _close(S["TKT"].w, {"nucleus": fn(rt), "cytosol": 1 - fn(rt)}, 1e-12)  # imaging N x Itzhak C: a mixture, not 100 %
    assert S["TKT"].label[-1] == "nc:hpa+prior"
    _close(S["SOLO"].w, {"nucleus": min(fn(10), 0.9), "cytosol": 1 - min(fn(10), 0.9)}, 1e-12)   # bound: cytosol added
    _close(S["CY"].w, {"nucleus": 0.1, "cytosol": 0.9}, 1e-12)              # single C category: 0.1 nuclear
    _close(S["H2BC12"].w, {"nucleus": 1.0}, 1e-12)
    assert S["H2BC12"].label[-1] == "nc:chromatin"
    assert S["H1-9"].w == {"cytosol": 1.0}                                   # chromatin gene without a nucleus key: unchanged
    T = C.nc_table(S)
    assert set(T.gene) == set(g) and bool(T.set_index("gene").at["H2BC12", "chromatin"])
    assert T.set_index("gene").at["TKT", "ex_c_other"] == 1                  # Itzhak C is exclusive on the cytoplasm side


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "ok")
