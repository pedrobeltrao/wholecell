"""No silent copy losses in the pack counts (testcell; about a minute on 2 CPUs, no pytest needed):

    PYTHONHASHSEED=0 python tests/test_pack_ledger.py        (or: pytest tests/test_pack_ledger.py)

The testcell frozen library and a synthetic pack table (every test library species, 400 copies) plus one species
without a library entry (`NOLIB`, an assembly of two library species):
  - no ledger reason: the counts step fails (unexplained loss above pack.ledger.max_unexplained_mass_frac)
  - library/decomposed.tsv (library step `fallback`, decompose): NOLIB's copies go to its subunits (copies x
    stoichiometry), pack/copy_ledger.tsv records the move, nothing is lost
  - pack.ledger.explain_no_structure (reference presets): NOLIB is dropped and recorded as explained
  - an index drop with return_to_pool: the dropped copies (pack/index_drops.tsv) are added back before the counts
  - an index drop without a reason is a config error
"""
import json
import os
import shutil
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("NUMBA_NUM_THREADS", "2")
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")

OVR = {"pack": {"meas_h_nm": 16.0}}          # coarse tile measures (as test_pack_smoke)
NOLIB_COPIES = 400


def _cfg(extra=None):
    from wholecell.config import load
    o = json.loads(json.dumps(OVR))
    for k, v in (extra or {}).items():
        o["pack"][k] = v
    return load("testcell", overrides=o)


def _pack_table(cfg, path):
    """every testcell library species (400 copies) + NOLIB (complex of the first two soluble species, no library entry)"""
    import pandas as pd
    fi = cfg["frozen_inputs"][cfg["library"]["frozen"]]
    rows = []
    for f in sorted(os.listdir(fi["path"])):
        if not f.endswith(".npz"):
            continue
        z = np.load(os.path.join(fi["path"], f)); tm = "membrane" in z.files
        rows.append(dict(ing_id=f[:-4], kind="protein", genes=f[:-4], stoich="1", copies=400, mass_da=50000.0, mass_pg=0.0,
                         comp_weights=json.dumps({"plasma_membrane": 1.0} if tm else {"cytosol": 1.0}), transmembrane=tm,
                         n_tm=1 if tm else 0, struct_src="test", uniprot="", signal_frac=0.0, mito_om_frac=0.0))
    sol = [r["ing_id"] for r in rows if not r["transmembrane"]][:2]
    rows.append(dict(ing_id="NOLIB", kind="complex", genes=";".join(sol), stoich="2;1", copies=NOLIB_COPIES, mass_da=150000.0,
                     mass_pg=0.0, comp_weights=json.dumps({"cytosol": 1.0}), transmembrane=False, n_tm=0, struct_src="pdb:none",
                     uniprot="", signal_frac=0.0, mito_om_frac=0.0))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return sol


def _bdir(cfg, tmp, name, tiles_from=None):
    from wholecell.library import stage as lib
    b = os.path.join(tmp, name)
    lib.run(cfg, b)
    os.makedirs(os.path.join(b, "preplaced"), exist_ok=True)
    json.dump({"instances": [], "beads": [], "drops": []}, open(os.path.join(b, "preplaced", "preplaced_manifest.json"), "w"))
    if tiles_from:
        os.makedirs(os.path.join(b, "pack"), exist_ok=True)
        for f in ("tiles.parquet", "tiles_meta.json"):
            shutil.copy(os.path.join(tiles_from, "pack", f), os.path.join(b, "pack", f))
    return b


def _requested(b):
    L = json.load(open(os.path.join(b, "pack", "libs.json")))
    cz = np.load(os.path.join(b, "pack", "counts.npz"))
    per = np.bincount(cz["lib"], weights=cz["n"], minlength=len(L["libs"]))
    return dict(zip(L["libs"], per.astype(np.int64))), L


def test_pack_ledger():
    import pandas as pd
    from wholecell.pack import stage
    cfg = _cfg()
    tmp = tempfile.mkdtemp(prefix="wc_pack_ledger_")
    try:
        pt = os.path.join(tmp, "pack_table.tsv")
        sol = _pack_table(cfg, pt)
        b0 = _bdir(cfg, tmp, "fail")
        stage.STEPS["tiles"](cfg, b0, {"workers": "2"})
        # 1. a species without a structure and no reason: the counts step fails
        try:
            stage.STEPS["counts"](cfg, b0, {"pack_table": pt})
            raise AssertionError("an unexplained copy loss above the limit must fail the counts step")
        except RuntimeError as e:
            assert "max_unexplained_mass_frac" in str(e) and "NOLIB" in str(e)
        L = pd.read_csv(os.path.join(b0, "pack", "copy_ledger.tsv"), sep="\t")
        assert L[["event", "ing_id", "copies"]].values.tolist() == [["no_structure", "NOLIB", NOLIB_COPIES]] and not L.explained[0]
        # 2. decomposed by the library: copies moved to the subunits, nothing lost
        b1 = _bdir(cfg, tmp, "decomposed", tiles_from=b0)
        pd.DataFrame([dict(parent="NOLIB", subunit=sol[0], gene=sol[0], uniprot="", n=2, mass_da=50000.0, transmembrane=False, built="species"),
                      dict(parent="NOLIB", subunit=sol[1], gene=sol[1], uniprot="", n=1, mass_da=50000.0, transmembrane=False, built="species")]
                     ).to_csv(os.path.join(b1, "library", "decomposed.tsv"), sep="\t", index=False)
        stage.STEPS["counts"](cfg, b1, {"pack_table": pt})
        req, info = _requested(b1)
        assert "NOLIB" not in req and info["ledger"]["unexplained_mass_frac"] == 0.0
        L = pd.read_csv(os.path.join(b1, "pack", "copy_ledger.tsv"), sep="\t")
        assert sorted(zip(L.event, L.ing_id, L.copies)) == sorted([("decomposed_in", sol[0], 2 * NOLIB_COPIES),
                                                                   ("decomposed_in", sol[1], NOLIB_COPIES),
                                                                   ("decomposed_out", "NOLIB", NOLIB_COPIES)])
        # 3. a reference preset that declares the old drop: recorded as explained, no failure
        cfg3 = _cfg({"ledger": {"explain_no_structure": "test reference: dropped as before"}})
        b3 = _bdir(cfg3, tmp, "explained", tiles_from=b0)
        stage.STEPS["counts"](cfg3, b3, {"pack_table": pt})
        req3, info3 = _requested(b3)
        assert "NOLIB" not in req3 and info3["ledger"]["lost_explained_pg"] > 0 and info3["ledger"]["lost_unexplained_pg"] == 0
        # random rounding / tile split: the decomposition adds exactly 3 x 400 requested copies (expected, then rounded)
        assert abs((sum(req.values()) - sum(req3.values())) - 3 * NOLIB_COPIES) <= len(req)
        # 4. an index drop returned to the pool: its copies come back before the counts
        drop = {"name": "t_drop", "libs": [sol[0]], "mode": "inside_labels", "labels": ["cytosol"], "reason": "test",
                "return_to_pool": {sol[0]: sol[0]}}
        cfg4 = _cfg({"ledger": {"explain_no_structure": "test"}, "index_drops": [drop]})
        b4 = _bdir(cfg4, tmp, "returned", tiles_from=b0)
        pd.DataFrame([dict(name="t_drop", lib=sol[0], n=150, returned_to=sol[0], reason="test")]).to_csv(
            os.path.join(b4, "pack", "index_drops.tsv"), sep="\t", index=False)
        stage.STEPS["counts"](cfg4, b4, {"pack_table": pt})
        L = pd.read_csv(os.path.join(b4, "pack", "copy_ledger.tsv"), sep="\t")
        assert ("index_drop_returned", sol[0], 150) in set(zip(L.event, L.ing_id, L.copies))
        req4, _ = _requested(b4)
        assert abs(req4[sol[0]] - req3[sol[0]] - 150) <= 2
        # 5. an index drop without a reason is a config error
        cfg5 = _cfg({"index_drops": [{"name": "no_reason", "libs": [sol[0]], "mode": "inside_labels", "labels": ["cytosol"]}]})
        b5 = _bdir(cfg5, tmp, "no_reason", tiles_from=b0)
        try:
            stage.STEPS["index"](cfg5, b5, {})
            raise AssertionError("an index drop without a reason must be refused")
        except ValueError as e:
            assert "reason" in str(e)
        print("pack ledger OK")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_pack_ledger()
