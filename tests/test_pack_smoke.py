"""Pack-stage smoke test on the synthetic test cell (a few minutes on 2 CPUs; no pytest needed):

    PYTHONHASHSEED=0 python tests/test_pack_smoke.py        (or: pytest tests/test_pack_smoke.py)

library (frozen test library) -> empty pre-pass manifest -> tiles (coarse measures) -> counts (synthetic pack table:
every test library species, 400 copies, cytosol or plasma membrane) -> index -> pack of a few colour-0 tiles, twice:
with 1 and with 2 workers. Checks: region table consistency, instances placed, results identical between the two runs
(determinism and worker-count independence), placed instances inside their tile + halo.
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

OVR = {"pack": {"meas_h_nm": 16.0}}          # coarse tile measures: the test is about the pipeline, not the numbers


def _cfg():
    from wholecell.config import load
    return load("testcell", overrides=OVR)


def _pack_table(cfg, path):
    import pandas as pd
    lc = cfg["library"]
    fi = cfg["frozen_inputs"][lc["frozen"]] if lc.get("frozen") else {"path": lc["frozen_dir"]}
    rows = []
    for f in sorted(os.listdir(fi["path"])):
        if not f.endswith(".npz"):
            continue
        z = np.load(os.path.join(fi["path"], f)); tm = "membrane" in z.files
        rows.append(dict(ing_id=f[:-4], kind="protein", genes=f[:-4], stoich="1", copies=400, mass_da=50000.0, mass_pg=0.0,
                         comp_weights=json.dumps({"plasma_membrane": 1.0} if tm else {"cytosol": 1.0}), transmembrane=tm,
                         n_tm=1 if tm else 0, struct_src="test", uniprot="", signal_frac=0.0, mito_om_frac=0.0))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return len(rows)


def _prep(cfg, bdir):
    from wholecell.library import stage as lib
    lib.run(cfg, bdir)
    os.makedirs(os.path.join(bdir, "preplaced"), exist_ok=True)
    json.dump({"instances": [], "beads": [], "drops": []}, open(os.path.join(bdir, "preplaced", "preplaced_manifest.json"), "w"))


def test_pack_smoke():
    os.environ["NUMBA_NUM_THREADS"] = "1"
    os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
    cfg = _cfg()
    from wholecell.pack import stage, geomload
    from wholecell.pack.regions import Regions
    g = geomload.load(cfg)
    reg = Regions(cfg["pack"], g)
    assert reg.RNAMES == [r["name"] for r in cfg["pack"]["regions"]]
    assert reg.NSH0 == reg.NV + len(reg.MEM_REG)
    tmp = tempfile.mkdtemp(prefix="wc_pack_smoke_")
    try:
        b1, b2 = os.path.join(tmp, "w1"), os.path.join(tmp, "w2")
        pt = os.path.join(tmp, "pack_table.tsv"); n_sp = _pack_table(cfg, pt)
        assert n_sp > 0
        for b in (b1, b2):
            _prep(cfg, b)
        stage.STEPS["tiles"](cfg, b1, {"workers": "2"})
        stage.STEPS["counts"](cfg, b1, {"pack_table": pt})
        stage.STEPS["index"](cfg, b1, {})
        os.makedirs(os.path.join(b2, "pack"), exist_ok=True)
        for f in ("tiles.parquet", "tiles_meta.json", "counts.npz", "libs.json", "preplaced_types.npy", "preplaced_pos.npy",
                  "preplaced_quat.npy", "preplaced_off.npy", "preplaced_libs.json"):
            shutil.copy(os.path.join(b1, "pack", f), os.path.join(b2, "pack", f))
        import pandas as pd
        M = pd.read_parquet(os.path.join(b1, "pack", "tiles.parquet"))
        cz = np.load(os.path.join(b1, "pack", "counts.npz"))
        assert cz["n"].sum() > 0
        Tij = M[["ti", "tj", "tk"]].values
        col0 = np.nonzero(((Tij[:, 0] % 2) * 4 + (Tij[:, 1] % 2) * 2 + (Tij[:, 2] % 2)) == 0)[0]
        njob = np.array([int(cz["n"][cz["tile"] == t].sum()) for t in col0])
        busy = [int(t) for t in col0[np.argsort(-njob, kind="stable")][:3] if njob[list(col0).index(t)] > 0]
        assert busy, "no colour-0 tile with jobs"
        sel = ",".join(str(t) for t in busy)
        stage.STEPS["pack"](cfg, b1, {"workers": "1", "tiles": sel, "phases": "0"})
        stage.STEPS["pack"](cfg, b2, {"workers": "2", "tiles": sel, "phases": "0"})
        meta = json.load(open(os.path.join(b1, "pack", "tiles_meta.json")))
        tot = 0
        for t in busy:
            z1 = np.load(os.path.join(b1, "pack", "tiles", f"{t}.npz")); z2 = np.load(os.path.join(b2, "pack", "tiles", f"{t}.npz"))
            for k in z1.files:
                assert np.array_equal(z1[k], z2[k]), f"tile {t} array {k} differs between 1 and 2 workers"
            lo = Tij[t] * meta["TILE"] + np.array(meta["lo0"])
            if len(z1["pos"]):
                assert np.all(z1["pos"] >= lo - 40.0 - 1e-3) and np.all(z1["pos"] <= lo + meta["TILE"] + 40.0 + 1e-3)
            assert int(z1["job_placed"].sum()) == len(z1["types"])
            tot += len(z1["types"])
        assert tot > 0, "nothing placed"
        print(f"pack smoke OK: {len(M)} tiles, {int(cz['n'].sum())} requested, {tot} placed in tiles {busy}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    test_pack_smoke()
