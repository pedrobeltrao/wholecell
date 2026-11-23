"""Library stage tests (seconds; synthetic inputs + the testcell frozen library; no yeast or Jurkat inputs needed).

    PYTHONHASHSEED=0 python tests/test_library.py        (or pytest tests/test_library.py)

  test_sphere_tree_and_rotation   k = ceil(n / per) spheres, single-sphere case, rot_to_z maps v to +z
  test_build_mode_synthetic       library.mode build end to end: frozen start entry (md5-checked), an AlphaFold model from
                                  a synthetic AFDB tar (pLDDT filter, library_af.json with a missing protein), sealed
                                  read-only lib/, manifest origins, pack_shapes (pack_trees_15 + pack_meta); a rerun is
                                  skipped; a lib/ without manifest is refused (never edited in place)
  test_jurkat_steps_registered    the Jurkat builder steps (B3) are library steps
  test_start_subset               library.start as a subset of a frozen library (TSV file, md5): only those files copied,
                                  a wrong md5 is refused
  test_adaptive_trees_synthetic   adaptive_trees (legacy lib_refine_st through the core step): an extended chain gets a
                                  finer tree (<= 2.5x its protein volume), st_c30 keeps the per-30 tree; a compact one is kept
  test_jurkat_config              jurkat.yaml: mode build, steps known, start files in the v7 manifest with its md5,
                                  species table rows; jurkat_v7_port.yaml keeps the frozen v7 library
  test_testcell_frozen_unchanged  the testcell frozen import gives the committed manifest's md5s
  test_placeholder_geometry       fallback placeholder: TM helices across the slab, soluble parts outside, sides from the
                                  topology; soluble ball of the protein volume; deterministic
  test_fallback_and_ledger        step fallback (placeholders, an assembly split into subunits) and the species ledger;
                                  a species without a structure above library.ledger.max_missing_mass_frac fails the stage
"""
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("OMP_NUM_THREADS", "1")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell.library import build as LB, builders as BL, stage as LS  # noqa: E402


def _md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def test_sphere_tree_and_rotation():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(95, 3))
    c, r = BL.sphere_tree(x, per=30)
    assert len(c) == 4 and len(r) == 4 and (r > BL.BEAD_R).all()
    c1, r1 = BL.sphere_tree(x[:20], per=30)
    assert c1.shape == (1, 3) and np.allclose(c1[0], x[:20].mean(0))
    v = np.array([0.3, -2.0, 0.5])
    assert np.allclose(BL.rot_to_z(v) @ (v / np.linalg.norm(v)), [0, 0, 1])
    assert BL.parse_ranges('TRANSMEM 5..25; /note="Helical"', "TRANSMEM") == [(5, 25, "Helical")]


def _af_cif(n=60, seed=3):
    import gemmi
    rng = np.random.default_rng(seed)
    st = gemmi.Structure(); st.name = "AF"
    try:
        m = gemmi.Model(1)                                         # gemmi >= 0.7
    except TypeError:
        m = gemmi.Model("1")
    ch = gemmi.Chain("A")
    xyz = np.cumsum(rng.normal(scale=2.2, size=(n, 3)), axis=0)
    for i in range(n):
        res = gemmi.Residue(); res.name = "ALA"; res.seqid = gemmi.SeqId(i + 1, " ")
        a = gemmi.Atom(); a.name = "CA"; a.element = gemmi.Element("C"); a.pos = gemmi.Position(*xyz[i])
        a.b_iso = 30.0 if i < 12 else 90.0                       # 12 low-confidence residues
        res.add_atom(a); ch.add_residue(res)
    m.add_chain(ch); st.add_model(m); st.setup_entities()
    return st.make_mmcif_document().as_string()


def _setup(tmp):
    start = os.path.join(tmp, "start"); os.makedirs(start)
    np.savez_compressed(os.path.join(start, "TEMPLATE.npz"), beads=np.zeros((3, 3), np.float32), kind=np.zeros(3, np.int8),
                        st_c=np.zeros((1, 3), np.float32), st_r=np.ones(1, np.float32))
    man = os.path.join(tmp, "start_manifest.tsv")
    open(man, "w").write(f"file\tmd5\tbytes\nTEMPLATE.npz\t{_md5(os.path.join(start, 'TEMPLATE.npz'))}\t0\n")
    tar = os.path.join(tmp, "afdb.tar")
    with tarfile.open(tar, "w") as tf:
        data = gzip.compress(_af_cif().encode())
        ti = tarfile.TarInfo("AF-P00001-F1-model_v6.cif.gz"); ti.size = len(data)
        tf.addfile(ti, io.BytesIO(data))
    sp = os.path.join(tmp, "species.tsv")
    open(sp, "w").write("ing_id\tkind\tgene\tuniprot\tpdb\tmembrane\nPROT1\tprotein\tPROT1\tP00001\t\t\n"
                        "PROT2\tprotein\tPROT2\tP00002\t\t\n")
    topo = os.path.join(tmp, "topo.tsv")
    open(topo, "w").write("Entry\tTopological domain\tTransmembrane\tIntramembrane\nP00001\t\t\t\n")
    cfg = {"_site": {"repo": REPO},
           "frozen_inputs": {"start_x": {"path": start, "manifest": man}},
           "library": {"mode": "build", "start": "start_x", "species": sp, "afdb_tar": tar, "pdb_dir": tmp,
                       "topology": {"path": topo}, "nproc": 1, "ledger": {"max_missing_mass_frac": 0.5},
                       "steps": [{"step": "af_models", "plddt_min": 50, "per": 30}], "pack_shapes": {"per": 15}}}
    return cfg


def test_build_mode_synthetic():
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _setup(tmp)
        bdir = os.path.join(tmp, "build")
        LS.run(cfg, bdir, step="import")
        lib = os.path.join(bdir, "library", "lib")
        assert sorted(os.listdir(lib)) == ["PROT1.npz", "TEMPLATE.npz"]
        z = np.load(os.path.join(lib, "PROT1.npz"))
        assert z["beads"].shape == (48, 3) and abs(float(z["plddt_frac_kept"]) - 48 / 60) < 1e-6
        assert np.allclose(z["beads"].mean(0), 0, atol=1e-5) and len(z["st_r"]) == 2
        assert not os.access(os.path.join(lib, "PROT1.npz"), os.W_OK) or os.geteuid() == 0
        af = json.load(open(os.path.join(bdir, "library", "library_af.json")))
        assert af == {"written": 1, "failed": {}, "missing": ["PROT2"]}
        rows = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(bdir, "library", "library_manifest.tsv"))]
        assert rows[0] == ["species", "path", "md5", "origin"]
        assert {r[0]: r[3] for r in rows[1:]} == {"PROT1": "af_models", "TEMPLATE": "start"}
        led = {r[0]: r[4:6] for r in (ln.rstrip("\n").split("\t") for ln in open(os.path.join(bdir, "library", "species_ledger.tsv")))}
        assert led["PROT1"] == ["built", "afdb_monomer"] and led["PROT2"] == ["missing", ""]
        pt = np.load(os.path.join(bdir, "library", "pack_shapes", "pack_trees_15.npz"))
        assert list(pt["ids"]) == ["PROT1", "TEMPLATE"] and pt["off"][-1] == len(pt["R"]) and pt["off"][1] == 4
        pm = np.load(os.path.join(bdir, "library", "pack_shapes", "pack_meta.npz"))
        assert list(pm["src"]) == ["AF", "AF"] and pm["BR"][0] > 0
        LS.run(cfg, bdir)                                         # sealed: skipped
        os.remove(os.path.join(bdir, "library", "library_manifest.tsv"))
        try:
            LS.run(cfg, bdir)
            raise AssertionError("an existing lib/ without manifest must be refused")
        except RuntimeError as e:
            assert "never edited in place" in str(e)


def test_jurkat_steps_registered():
    for name in ("jurkat_assemblies", "jurkat_monomers", "reorient_tm", "tm_route", "adaptive_trees"):
        assert name in LB.STEPS and callable(LB.STEPS[name])


def _ctx(tmp, lc, cfg_extra=None):
    cfg = dict(library=lc, frozen_inputs={}, _site={"repo": REPO})
    cfg.update(cfg_extra or {})
    X = LB.context(cfg, os.path.join(tmp, "b"))
    os.makedirs(X.partial)
    return X


def test_start_subset():
    with tempfile.TemporaryDirectory() as tmp:
        src = os.path.join(tmp, "frozen")
        os.makedirs(src)
        for n in ("A.npz", "B.npz", "C.json"):
            open(os.path.join(src, n), "w").write(n)
        LS.make_manifest(src, os.path.join(tmp, "man.tsv"))
        open(os.path.join(tmp, "sub.tsv"), "w").write(f"file\tmd5\treason\nB.npz\t{_md5(os.path.join(src, 'B.npz'))}\tt\n"
                                                    f"C.json\t{_md5(os.path.join(src, 'C.json'))}\tt\n")
        fi = {"F": {"path": src, "manifest": os.path.join(tmp, "man.tsv")}}
        X = _ctx(tmp, {"start": {"frozen": "F", "files": os.path.join(tmp, "sub.tsv")}}, {"frozen_inputs": fi})
        LB.step_start(X, {})
        assert sorted(os.listdir(X.partial)) == ["B.npz", "C.json"] and X.origin == {"B.npz": "start", "C.json": "start"}
        open(os.path.join(tmp, "bad.tsv"), "w").write("file\tmd5\treason\nA.npz\t0\tt\n")
        X2 = _ctx(os.path.join(tmp, "x2"), {"start": {"frozen": "F", "files": os.path.join(tmp, "bad.tsv")}}, {"frozen_inputs": fi})
        try:
            LB.step_start(X2, {})
            raise AssertionError("a start file with a wrong md5 must be refused")
        except RuntimeError as e:
            assert "md5" in str(e)


def test_adaptive_trees_synthetic():
    with tempfile.TemporaryDirectory() as tmp:
        sp = os.path.join(tmp, "species.tsv")
        open(sp, "w").write("ing_id\tkind\tgenes\taccs\tstoich\tstruct_src\tpdb\tuniprot\ttransmembrane\tsignal_frac\tmass_da\n"
                            "EXT\tprotein\tEXT\tP1\t\tafdb_monomer:P1\t\tP1\tFalse\t0\t22000\n"
                            "CMP\tprotein\tCMP\tP2\t\tafdb_monomer:P2\t\tP2\tFalse\t0\t22000\n")
        X = _ctx(tmp, {"species": sp})
        t = np.arange(200) * 0.38                                     # 200 residues: extended chain vs compact ball
        ext = np.c_[t, 0.3 * np.sin(t), 0.3 * np.cos(t)] - [t.mean(), 0, 0]
        rng = np.random.default_rng(0)
        cmp_ = rng.normal(size=(200, 3)); cmp_ *= (1.9 * rng.random(200) ** (1 / 3) / np.linalg.norm(cmp_, axis=1))[:, None]
        for name, b in (("EXT", ext), ("CMP", cmp_)):
            c, r = BL.sphere_tree(b, per=30)
            np.savez_compressed(os.path.join(X.partial, name + ".npz"), beads=b.astype(np.float32), kind=np.zeros(200, np.int8),
                                st_c=c.astype(np.float32), st_r=r.astype(np.float32))
        X["nproc"] = 1
        LB.STEPS["adaptive_trees"](X, {"max_ratio": 2.5, "pers": [15, 8, 5]})
        e = np.load(os.path.join(X.partial, "EXT.npz")); k = np.load(os.path.join(X.partial, "CMP.npz"))
        assert "st_c30" in e.files and len(e["st_r"]) > len(e["st_r30"]) == 7
        assert "st_c30" not in k.files and len(k["st_r"]) == 7
        assert X.origin == {"EXT.npz": "adaptive_trees"} and X.results["adaptive_trees"]["action"].get("kept") == 1


def test_jurkat_config():
    import pandas as pd
    from wholecell.config import load
    cfg = load("jurkat")
    lc = cfg["library"]
    assert lc["mode"] == "build" and all(s["step"] in LB.STEPS for s in lc["steps"])
    st = pd.read_csv(os.path.join(REPO, lc["start"]["files"]), sep="\t")
    man = dict(LB.read_any_manifest(cfg["frozen_inputs"]["lib_structures_v7"]["manifest"]))
    assert len(st) == 4 and all(man[f] == m for f, m in zip(st.file, st.md5))
    assert lc["species"]["from"] == "pack_table"                  # rebuild v9: species = this build's pack table
    sp = pd.read_csv(os.path.join(REPO, "plugins/jurkat/data/library/species_v7.tsv.gz"), sep="\t", keep_default_na=False)
    assert len(sp) == 16279 and sp.ing_id.is_unique and not set(st.file.str[:-4]) & set(sp.ing_id)   # v7 table (history)
    steps = [s["step"] for s in lc["steps"]]
    assert steps.index("fallback") < steps.index("adaptive_trees")
    assert load("jurkat_v7_port")["library"]["mode"] == "frozen"


def test_testcell_frozen_unchanged():
    from wholecell.config import load
    cfg = load("testcell")
    assert cfg["library"]["mode"] == "frozen"
    with tempfile.TemporaryDirectory() as tmp:
        LS.run(cfg, tmp)
        got = [ln.rstrip("\n").split("\t") for ln in open(os.path.join(tmp, "library", "library_manifest.tsv"))][1:]
        ref = LS.read_manifest(os.path.join(REPO, "tests", "data", "testcell", "library_manifest.tsv.gz"))
        assert [(g[0] + ".npz", g[2]) for g in got] == [(r["file"], r["md5"]) for r in ref]



def test_placeholder_geometry():
    """placeholder of a 2-pass TM protein: helices across the slab, soluble stretches outside |z| >= 2.3 nm, N-terminus on
    the cytoplasmic side (+z) by the topology, ball of the protein volume for a soluble protein"""
    tm = [(41, 61, "Helical"), (81, 101, "Helical")]
    topo = [(1, 40, "Cytoplasmic"), (62, 80, "Extracellular"), (102, 150, "Cytoplasmic")]
    xyz, memb = BL.placeholder_beads(150, tm, topo)
    assert memb and xyz.shape == (150, 3)
    seg = np.zeros(150, bool); seg[40:61] = seg[80:101] = True
    assert np.abs(xyz[seg, 2]).max() < 1.6 and np.abs(xyz[~seg, 2]).min() >= 2.3 - 1e-9
    assert xyz[:40, 2].min() > 0 and xyz[61:80, 2].max() < 0 and xyz[101:, 2].min() > 0
    xyz2, _ = BL.placeholder_beads(150, tm, [(1, 40, "Lumenal"), (62, 80, "Cytoplasmic"), (102, 150, "Lumenal")])
    assert xyz2[:40, 2].max() < 0                                   # N-terminus lumenal -> -z
    b, memb = BL.placeholder_beads(200)
    R = (3 * 200 * BL.RES_VOL_NM3 / (4 * np.pi)) ** (1 / 3)
    assert not memb and np.linalg.norm(b, axis=1).max() <= 2 * R and np.allclose(b.mean(0), 0, atol=1e-9)
    assert np.allclose(b, BL.placeholder_beads(200)[0])            # deterministic


def test_fallback_and_ledger():
    """step fallback: proteins without a model get a placeholder (TM one oriented, membrane flag), an assembly without a
    structure is split into its subunits (decomposed.tsv), the species ledger records every status; a missing species
    above library.ledger.max_missing_mass_frac fails the stage"""
    import pandas as pd
    with tempfile.TemporaryDirectory() as tmp:
        cfg = _setup(tmp)
        sp = os.path.join(tmp, "species_j.tsv")
        open(sp, "w").write("ing_id\tkind\tgenes\taccs\tstoich\tstruct_src\tpdb\tuniprot\ttransmembrane\tsignal_frac\tmass_da\tcopies\n"
                            "PROT1\tprotein\tPROT1\tP00001\t1\tafdb_monomer:P00001\t\tP00001\tFalse\t0\t6600\t1000\n"
                            "SOL\tprotein\tSOL\tP00002\t1\tafdb_monomer:P00002\t\tP00002\tFalse\t0\t13200\t500\n"
                            "TMP\tprotein\tTMP\tP00003\t1\tafdb_monomer:P00003\t\tP00003\tTrue\t0\t16500\t200\n"
                            "CPX\tcomplex\tPROT1;NEW\tP00001;P00004\t2;1\tpdb:none\tnone\t\tFalse\t0\t30000\t100\n")
        topo = os.path.join(tmp, "topo_j.tsv")
        open(topo, "w").write("Entry\tTopological domain\tTransmembrane\tIntramembrane\n"
                              'P00003\tTOPO_DOM 1..60; /note="Cytoplasmic"; TOPO_DOM 82..150; /note="Lumenal"\t'
                              'TRANSMEM 61..81; /note="Helical"\t\n')
        lens = os.path.join(tmp, "lens.tsv")
        open(lens, "w").write("Entry\tLength\tMass\nP00002\t120\t13200\nP00003\t150\t16500\nP00004\t90\t9900\n")
        cfg["library"].update(species=sp, ledger={"max_missing_mass_frac": 1e-4}, steps=[
            {"step": "af_models", "plddt_min": 50, "per": 30},
            {"step": "fallback", "topology": {"path": topo, "columns": ["acc", "topo", "tm", "intramem"]},
             "lengths": {"path": lens, "columns": {"acc": "Entry", "length": "Length", "mass": "Mass"}}, "decompose": True}])
        bdir = os.path.join(tmp, "build")
        LS.run(cfg, bdir, step="import")
        lib = os.path.join(bdir, "library", "lib")
        assert sorted(os.listdir(lib)) == ["NEW.npz", "PROT1.npz", "SOL.npz", "TEMPLATE.npz", "TMP.npz"]
        s = np.load(os.path.join(lib, "SOL.npz")); t = np.load(os.path.join(lib, "TMP.npz"))
        assert len(s["beads"]) == 120 and "membrane" not in s.files and str(s["src"]) == "fallback:placeholder"
        assert len(t["beads"]) == 150 and str(t["membrane"]) == "TM" and t["beads"][:60, 2].min() > 0
        L = pd.read_csv(os.path.join(bdir, "library", "species_ledger.tsv"), sep="\t", keep_default_na=False).set_index("ing_id")
        assert L.status.to_dict() == {"PROT1": "built", "SOL": "fallback", "TMP": "fallback", "CPX": "decomposed"}
        assert L.source.to_dict() == {"PROT1": "afdb_monomer", "SOL": "placeholder", "TMP": "placeholder", "CPX": ""}
        D = pd.read_csv(os.path.join(bdir, "library", "decomposed.tsv"), sep="\t")
        assert list(zip(D.subunit, D.n)) == [("PROT1", 2), ("NEW", 1)] and D.mass_da.tolist() == [6600.0, 9900.0]
        # without the fallback step: 3 of 4 species (by mass) have no structure -> the stage fails before sealing
        cfg["library"]["steps"] = [{"step": "af_models", "plddt_min": 50, "per": 30}]
        try:
            LS.run(cfg, os.path.join(tmp, "build2"), step="import")
            raise AssertionError("species without a structure above the ledger limit must fail the library stage")
        except RuntimeError as e:
            assert "max_missing_mass_frac" in str(e)

if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
