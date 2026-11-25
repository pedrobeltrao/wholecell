"""Library chain rule, v7 structure pick and the ATP synthase guard (seconds; synthetic mmCIF / SIFTS, offline):

    python tests/test_library_chains.py        (or pytest tests/test_library_chains.py)

  test_decide_rules            the three classes: members kept (all copies; gene match; declarations), foreign removed
                               (another human protein of a complete species, antibody / Fab by description, fusion
                               residues), ambiguous kept and flagged (unmapped / peptide, undeclared nucleic acid,
                               non-human chain, paralogue of a member, species with an empty accession entry); a monomer
                               keeps exactly one chain
  test_assembly_filter_check   jlib_build.build_assembly on a synthetic PDB assembly with the rule on: the foreign chains
                               (a partner protein, a Fab) out of the npz, the flagged ones in it, its chain record;
                               chains.check passes it and FAILS a model that holds a class-2 chain, a structure model
                               without a record, and one with more beads than its chains; a complex without a member
                               chain is not built (no_allowed_chains)
  test_afdb_chain_identity     AFDB complex model: chain accessions from the model's struct_ref; a chain of another protein
                               is removed
  test_core_builders_filter    builders.build_complex / build_pdb_protein with a chain filter (core steps)
  test_struct_src_from         library.species.struct_src_from: the v7 pick for the species the table covers, the
                               inventory struct_src for the others; both recorded
  test_atp_guard               ATP synthase dimers: placed < min_placed_frac x min(requested, rim capacity) fails; v8-like
                               numbers pass; null switches it off; jurkat.yaml sets it
  test_compare_geom            compare report: ectodomain height / Rg of a synthetic membrane model, flags
"""
import gzip
import json
import os
import shutil
import sys
import tempfile

os.environ.setdefault("PYTHONHASHSEED", "0")
os.environ.setdefault("OMP_NUM_THREADS", "1")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell.library import chains as CF  # noqa: E402

LEGACY = os.path.join(REPO, "wholecell", "library", "legacy")


# ------------------------------------------------------------------ synthetic structures
def _chain(gemmi, name, n, x0, nucleic=False):
    ch = gemmi.Chain(name)
    for i in range(n):
        r = gemmi.Residue()
        r.name = "DA" if nucleic else "ALA"
        r.seqid = gemmi.SeqId(i + 1, " ")
        r.label_seq = i + 1
        r.entity_type = gemmi.EntityType.Polymer
        a = gemmi.Atom()
        a.name = "P" if nucleic else "CA"
        a.element = gemmi.Element("P" if nucleic else "C")
        a.pos = gemmi.Position(x0 + 3.8 * (i % 10), 4.0 * (i // 10), 0.0)
        r.add_atom(a)
        ch.add_residue(r)
    return ch


def write_assembly(path, chains, desc=None):
    """chains: [(name, n_residues, nucleic)] -> mmCIF (gz) with entities, label_seq and the entity descriptions desc
    {chain: text}; chain i at x = 60 i A"""
    import gemmi
    st = gemmi.Structure()
    st.name = "TEST"
    m = gemmi.Model("1")
    for i, (nm, n, na) in enumerate(chains):
        m.add_chain(_chain(gemmi, nm, n, 60.0 * i, na))
    st.add_model(m)
    st.setup_entities()
    doc = st.make_mmcif_document()
    blk = doc.sole_block()
    ents = [(e.name, e.subchains) for e in st.entities]
    ch_of = {sc.subchain_id(): ch.name for ch in st[0] for sc in ch.subchains()}
    blk.set_mmcif_category("_entity.", {"id": [e for e, _ in ents], "type": ["polymer"] * len(ents),
                                        "pdbx_description": [(desc or {}).get(ch_of.get(sc[0], ""), "?") for _, sc in ents]})
    txt = doc.as_string()
    with gzip.open(path, "wt") as fh:
        fh.write(txt)
    return path


def sifts_table(path, rows):
    with gzip.open(path, "wt") as fh:
        fh.write("# test\nPDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END\n")
        for p, c, a, b, e in rows:
            fh.write(f"{p}\t{c}\t{a}\t{b}\t{e}\t{b}\t{e}\t{b}\t{e}\n")
    return path


def _info(spec):
    """chain_info from {chain: [(acc, n)] or ('na', n)}"""
    acc, chn, kind = [], [], []
    for c, parts in spec.items():
        for a, n in parts:
            acc += ["" if a in ("na", "-") else a] * n; chn += [c] * n; kind += [1 if a == "na" else 0] * n
    return np.array(acc), np.array(chn), np.array(kind, np.int8)


def test_decide_rules():
    pol = CF.Policy()
    acc, chn, kind = _info({"A": [("P1", 30)], "A2": [("P1", 30)], "B": [("P2", 20)], "C": [("na", 10)], "D": [("-", 5)],
                            "E": [("P1", 10), ("P00720", 6)], "F": [("Q9", 12)]})
    spec = CF.species_spec(dict(ing_id="CPX-1", kind="complex", genes="G1;G2;G3", accs="P1;P2;P3"))
    xyz = np.zeros((len(acc), 3)); pos = np.arange(len(acc))
    run = lambda sp, p, **k: CF.apply(xyz, kind, acc, pos, chn, sp, p, "pdb:test", **k)   # noqa: E731
    *_, c2, rec = run(spec, pol)
    assert sorted(set(c2)) == ["A", "A2", "B", "C", "D", "E"] and len(c2) == 30 + 30 + 20 + 10 + 5 + 10, (set(c2), len(c2))
    cls = {d["chain"]: (d["cls"], d["why"]) for d in rec["kept"]}
    assert cls["A2"] == ("member", "accession") and cls["C"] == ("flagged", "nucleic_acid_not_declared")
    assert cls["D"] == ("flagged", "unmapped_peptide"), cls
    rem = {(d["chain"], d["why"]) for d in rec["removed"]}
    assert rem == {("E", "fusion_residues"), ("F", "other_protein")}, rem
    assert rec["missing"] == ["P3"] and rec["n_struct"] == 105
    # gene match: F (Q9) maps to gene G3 of the species -> member
    pg = CF.Policy(gene_match=True, acc_genes={"Q9": frozenset({"G3", "ALIAS"})})
    *_, c3, rec3 = run(spec, pg)
    assert "F" in set(c3) and rec3["gene_matched"][0]["genes"] == "G3" and rec3["missing"] == ["P3"]
    # ambiguous non-members are kept and flagged: non-human chain, paralogue of a member, incomplete accession list
    *_, c4, r4 = run(spec, pol, tax={"F": {"8355"}})
    assert "F" in set(c4) and [d["why"] for d in r4["kept"] if d["chain"] == "F"] == ["non_human_possible_ortholog:8355"]
    pp = CF.Policy(paralogs={"G2": {"G2L"}}, acc_primary={"Q9": "G2L", "P2": "G2"})
    *_, c5, r5 = run(spec, pp)
    assert [d["why"] for d in r5["kept"] if d["chain"] == "F"] == ["paralogue_of_member:G2L"]
    inc = CF.species_spec(dict(ing_id="CPX-1", kind="complex", genes="G1;G2;G3", accs="P1;P2;"))
    *_, c6, r6 = run(inc, pol)
    assert [d["why"] for d in r6["kept"] if d["chain"] == "F"] == ["species_accessions_incomplete:G3"] and r6["unverifiable"] == ["G3"]
    # an antibody fragment (description) is foreign even without a SIFTS mapping
    *_, c7, r7 = run(spec, pol, desc_fn=lambda: {"D": "Fab heavy chain", "C": "DNA (147-MER)"})
    assert "D" not in set(c7) and ("D", "antibody_or_binder_by_description") in {(d["chain"], d["why"]) for d in r7["removed"]}
    assert not pol.foreign_by_desc("HLA-A2.1-RESTRICTED INFLUENZA A MATRIX EPITOPE") and not pol.foreign_by_desc("Myosin light chain")
    assert pol.foreign_by_desc("T-cell receptor beta chain") and pol.foreign_by_desc("Nanobody Nb35")
    # identified by name (UniProt reviewed name / gene symbol): a non-human ortholog of a member is a member, a
    # non-human chain of another complex is foreign, an unknown non-human one stays flagged
    pn = CF.Policy(name_accs={"subunit two": {"P2"}, "large ribosomal subunit protein ul2": {"Q7"}},
                   acc_primary={"P1": "G1", "P2": "G2", "Q7": "RPL8"})
    assert CF.classify_protein(spec, pn, ["X1"], "Subunit two", {"9986"})[:2] == ("member", "name_match:P2")
    assert CF.classify_protein(spec, pn, ["X2"], "uL2", {"9986"})[:2] == ("foreign", "other_protein_by_name:Q7")
    assert CF.classify_protein(spec, pn, [], "uL2", None)[0] == "foreign"
    assert CF.classify_protein(spec, pn, ["X3"], "something", {"9986"})[0] == "flagged"
    assert CF.protein_names("Large ribosomal subunit protein uL2 (60S ribosomal protein L8) [Cleaved into: x]") == [
        "Large ribosomal subunit protein uL2", "60S ribosomal protein L8"]
    # declarations: DNA of the species (regex on the description), an unmapped chain, an extra accession -> members
    pd_ = CF.Policy(nucleic={"CPX-1": "DNA"}, allow_chains={"CPX-1": {"D"}}, allow_accs={"CPX-1": {"Q9"}})
    *_, c8, r8 = run(spec, pd_, desc_fn=lambda: {"C": "DNA (147-MER)"})
    assert all(d["cls"] == "member" for d in r8["kept"]) and set(c8) == {"A", "A2", "B", "C", "D", "E", "F"}
    # a chain whose residues SIFTS does not map (other numbering) takes its chain-level SIFTS accession
    *_, c10, r10 = run(spec, pol, chain_accs={"D": {"Q8"}, "A": {"P1"}})
    assert "D" not in set(c10) and ("D", "other_protein") in {(d["chain"], d["why"]) for d in r10["removed"]}
    *_, c11, r11 = run(spec, pol, chain_accs={"D": {"P3"}})
    assert [d["cls"] for d in r11["kept"] if d["chain"] == "D"] == ["member"] and r11["missing"] == []
    # a monomer model keeps exactly one chain of the protein (the largest)
    acc, chn, kind = _info({"A": [("P1", 30)], "B": [("P1", 40)], "C": [("P2", 10)]})
    sp = CF.species_spec(dict(ing_id="P1G", kind="protein", genes="G1", accs="P1", uniprot="P1"))
    *_, c9, rec9 = CF.apply(np.zeros((80, 3)), kind, acc, np.arange(80), chn, sp, pol, "pdb:test")
    assert set(c9) == {"B"} and {d["why"] for d in rec9["removed"]} == {"extra_copy_monomer", "other_protein"}


def _legacy(tmp, pdb_dir, sifts):
    if LEGACY not in sys.path:
        sys.path.insert(0, LEGACY)
    import build_structs as bs, jlib_build as jb  # noqa: E401
    bs.LIB = jb.LIB = tmp
    jb.PDBDIR = pdb_dir
    jb.pdb_file = lambda pid: os.path.join(pdb_dir, f"{pid}-assembly1.cif.gz")
    import pandas as pd
    S = pd.read_csv(sifts, sep="\t", comment="#", usecols=["PDB", "CHAIN", "SP_PRIMARY", "RES_BEG", "RES_END", "SP_BEG", "SP_END"])
    jb.sifts = lambda: S
    return bs, jb


def test_assembly_filter_check():
    tmp = tempfile.mkdtemp()
    try:
        pdb_dir = os.path.join(tmp, "pdb"); os.makedirs(pdb_dir)
        lib = os.path.join(tmp, "lib"); os.makedirs(lib)
        write_assembly(os.path.join(pdb_dir, "9xyz-assembly1.cif.gz"),
                       [("A", 30, False), ("B", 25, False), ("H", 20, False), ("L", 18, False), ("D", 12, True), ("P", 9, False)],
                       desc={"A": "subunit alpha", "B": "subunit beta", "H": "partner kinase", "L": "Fab light chain",
                             "D": "DNA (12-MER)", "P": "bound peptide"})
        sifts = sifts_table(os.path.join(tmp, "sifts.tsv.gz"), [("9xyz", "A", "P11111", 1, 30), ("9xyz", "B", "P22222", 1, 25),
                                                                ("9xyz", "H", "P99999", 1, 20)])   # L (Fab), P unmapped
        bs, jb = _legacy(lib, pdb_dir, sifts)
        jb.CHAIN_POLICY = CF.Policy()
        row = dict(ing_id="CPX-T", struct_src="pdb:9xyz", transmembrane=False, accs=["P11111", "P22222"], kind="complex",
                   genes="GA;GB", uniprot="")
        ing, status, note = jb.build_assembly(row)
        assert status == "ok", (status, note)
        z = np.load(os.path.join(lib, "CPX-T.npz"))
        assert len(z["beads"]) == 30 + 25 + 12 + 9 and int(z["n_struct_beads"]) == 76 and int(z["beads_added"]) == 0
        rec = json.loads(str(z["chain_set"]))
        assert [(d["chain"], d["cls"]) for d in rec["kept"]] == [("A", "member"), ("B", "member"), ("D", "flagged"), ("P", "flagged")]
        assert {d["chain"]: d["why"] for d in rec["removed"]} == {"H": "other_protein", "L": "antibody_or_binder_by_description"}
        rows = [dict(ing_id="CPX-T", kind="complex", genes="GA;GB", accs="P11111;P22222", uniprot="")]
        ok, err = CF.check(lib, rows, CF.Policy(), sifts)
        assert len(ok) == 1 and ok[0]["flagged"] == 2 and not err, err
        # off: every chain (legacy); the check fails such a model (no record)
        jb.CHAIN_POLICY = None
        jb.build_assembly(dict(row, ing_id="CPX-U"))
        assert len(np.load(os.path.join(lib, "CPX-U.npz"))["beads"]) == 30 + 25 + 20 + 18 + 12 + 9
        _, err = CF.check(lib, rows + [dict(rows[0], ing_id="CPX-U")], CF.Policy(), sifts)
        assert [e[0] for e in err] == ["CPX-U"] and "no chain record" in err[0][2]
        # a record keeping the partner chain (class 2), and a model with more beads than its chains: both fail
        d = dict(np.load(os.path.join(lib, "CPX-T.npz")))
        r2 = json.loads(str(d["chain_set"]))
        r2["kept"].append(dict(chain="H", acc="P99999", n=20, cls="flagged", why="unmapped", type="protein", desc=""))
        r2["n_struct"] = 96
        np.savez(os.path.join(lib, "CPX-V.npz"), **dict(d, chain_set=np.array(json.dumps(r2)), n_struct_beads=np.int32(96)))
        np.savez(os.path.join(lib, "CPX-W.npz"), **dict(d, beads=np.zeros((90, 3), np.float32)))
        _, err = CF.check(lib, [dict(rows[0], ing_id=i) for i in ("CPX-T", "CPX-V", "CPX-W")], CF.Policy(), sifts)
        E = {e[0]: e[2] for e in err}
        assert set(E) == {"CPX-V", "CPX-W"} and "H: ['P99999'] other_protein (class 2)" in E["CPX-V"] and "90 beads > 76" in E["CPX-W"], E
        # a declared nucleic acid is a member; a complex without a member chain is not built
        jb.CHAIN_POLICY = CF.Policy(nucleic={"CPX-N": "DNA"})
        assert jb.build_assembly(dict(row, ing_id="CPX-N"))[1] == "ok"
        rn = json.loads(str(np.load(os.path.join(lib, "CPX-N.npz"))["chain_set"]))
        assert [d["cls"] for d in rn["kept"] if d["chain"] == "D"] == ["member"]
        _, err = CF.check(lib, [dict(rows[0], ing_id="CPX-N")], jb.CHAIN_POLICY, sifts)
        assert not err
        st = jb.build_assembly(dict(row, ing_id="CPX-Z", accs=["Q00000"]))
        assert st[1] == "no_allowed_chains" and not os.path.exists(os.path.join(lib, "CPX-Z.npz")), st
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_afdb_chain_identity():
    import gemmi
    tmp = tempfile.mkdtemp()
    try:
        p = os.path.join(tmp, "AF-0000000000000001.cif.gz")
        write_assembly(p, [("A", 20, False), ("B", 15, False)])
        txt = gzip.open(p, "rt").read().rstrip("\n") + "\n"
        txt += ("loop_\n_struct_ref.id\n_struct_ref.db_code\n_struct_ref.db_name\n_struct_ref.entity_id\n"
                "_struct_ref.pdbx_db_accession\n1 X_HUMAN UNP 1 P11111\n2 Y_HUMAN UNP 2 P33333\n"
                "loop_\n_struct_ref_seq.align_id\n_struct_ref_seq.ref_id\n_struct_ref_seq.pdbx_PDB_id_code\n"
                "_struct_ref_seq.pdbx_strand_id\n1 1 AF A\n2 2 AF B\n")
        with gzip.open(p, "wt") as fh:
            fh.write(txt)
        fa = CF.afdb_chain_accs(p)
        assert fa == {"A": "P11111", "B": "P33333"}, fa
        st = gemmi.read_structure(p)
        chn = np.array([ch.name for ch in st[0] for _ in ch])
        acc_legacy = np.array(["P11111"] * 20 + ["P22222"] * 15)          # chain order of the species' accession list
        spec = CF.species_spec(dict(ing_id="A-B", kind="afdb_heterodimer", genes="GA;GB", accs="P11111;P22222"))
        x, k, a, q, c, rec = CF.apply(np.zeros((35, 3)), np.zeros(35, np.int8), acc_legacy, np.arange(35), chn, spec,
                                      CF.Policy(), "afdb_complex:AF-0000000000000001",
                                      acc_decide=np.array([fa[i] for i in chn]))
        assert set(c) == {"A"} and rec["missing"] == ["P22222"] and rec["removed"][0]["acc"] == "P33333"
        assert rec["removed"][0]["why"] == "other_protein"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_core_builders_filter():
    from wholecell.library import builders as BL
    tmp = tempfile.mkdtemp()
    try:
        write_assembly(os.path.join(tmp, "1abc-assembly1.cif.gz"), [("A", 30, False), ("B", 30, False), ("C", 10, False)])
        cacc = {"A": "P1", "B": "P1", "C": "P7"}
        spec = CF.species_spec(dict(ing_id="G_oligo", kind="oligomer", genes="G", accs="P1"))

        def flt(chn, kind):
            acc = np.array([cacc[c] for c in chn])
            dec = CF.decide(spec, CF.chain_info(acc, chn, kind), CF.Policy())
            m = CF.mask(acc, chn, kind, dec)
            return m, CF.npz_extra(CF.record(dec, int(m.sum()), "pdb:1abc", len(m)))
        n, _, _ = BL.build_pdb_protein(tmp, "G_oligo", "1abc", tmp, chain_filter=flt)
        assert n == 60
        n, _, _ = BL.build_complex(tmp, "G_cx", "1abc", tmp, chain_filter=flt)
        z = np.load(os.path.join(tmp, "G_cx.npz"))
        assert n == 60 and str(z["src"]) == "pdb:1abc" and int(z["n_struct_beads"]) == 60
        n, _, _ = BL.build_pdb_protein(tmp, "G_all", "1abc", tmp)                   # off: every chain (yeast s0)
        assert n == 70
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_struct_src_from():
    import pandas as pd
    from wholecell.library import build as LB
    tmp = tempfile.mkdtemp()
    try:
        os.makedirs(os.path.join(tmp, "inventory"))
        pd.DataFrame([dict(ing_id="CPX-1", kind="complex", genes="A;B", stoich="1;1", struct_src="pdb:6eqa", uniprot="",
                           transmembrane=True, signal_frac=1.0, mass_da=5e4, copies=10),
                      dict(ing_id="CPX-2", kind="complex", genes="C;D", stoich="1;1", struct_src="pdb:1aaa", uniprot="",
                           transmembrane=False, signal_frac=0.0, mass_da=5e4, copies=5),
                      dict(ing_id="P", kind="protein", genes="P", stoich="1", struct_src="afdb_monomer:Q1", uniprot="Q1",
                           transmembrane=False, signal_frac=0.0, mass_da=1e4, copies=3)]).to_csv(
            os.path.join(tmp, "inventory", "pack_table.tsv"), sep="\t", index=False)
        v7 = os.path.join(tmp, "species_v7.tsv")
        pd.DataFrame([dict(ing_id="CPX-1", struct_src="pdb:2x4n"), dict(ing_id="P", struct_src="afdb_monomer:Q1")]).to_csv(
            v7, sep="\t", index=False)
        cfg = {"library": {"species": {"from": "pack_table", "struct_src_from": {"table": v7, "column": "struct_src"}}},
               "_site": {"repo": REPO}}
        X = LB.context(cfg, tmp)
        os.makedirs(X.out)
        LB.species_from_pack_table(X)
        S = pd.read_csv(X.species_path, sep="\t", keep_default_na=False).set_index("ing_id")
        assert S.struct_src.to_dict() == {"CPX-1": "pdb:2x4n", "CPX-2": "pdb:1aaa", "P": "afdb_monomer:Q1"}
        assert S.pdb["CPX-1"] == "2x4n" and S.struct_src_inventory["CPX-1"] == "pdb:6eqa"
        assert S.struct_src_from.to_dict() == {"CPX-1": "struct_src_from", "CPX-2": "inventory", "P": "struct_src_from"}
        assert X.results["species"]["struct_src_differs_from_inventory"] == 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def test_atp_guard():
    from wholecell import config as C
    from wholecell.preplaced.builders import atp_synthase_dimers as AD
    info = dict(rim_length_traced_um=3902.3, psi_deg=-80.0)
    p = dict(spacing_nm=12.0, min_placed_frac=0.5)
    r = AD.placed_guard(205152, 434344, info, p, log=lambda *a: None)              # v8: 63% of the rim capacity
    assert r["rim_capacity_dimers"] == 325191 and r["possible_dimers"] == 325191 and r["placed_frac_of_possible"] > 0.6
    try:
        AD.placed_guard(599, 434344, dict(info, psi_deg=-30.0), p, log=lambda *a: None)   # v9
        raise AssertionError("the guard did not fail")
    except RuntimeError as e:
        assert "only 599 dimers placed" in str(e) and "psi_deg" in str(e)
    AD.placed_guard(599, 434344, info, dict(p, min_placed_frac=None), log=lambda *a: None)    # off
    AD.placed_guard(10, 15, {}, p, log=lambda *a: None)                                   # random sampling: of requested
    j = [b for b in C.load("jurkat")["preplaced"]["serial"]["builders"] if b["builder"] == "core.atp_synthase_dimers"][0]
    assert j["params"]["min_placed_frac"] == 0.5


def test_compare_geom():
    from wholecell.library import compare as CP
    tmp = tempfile.mkdtemp()
    try:
        tmb = np.c_[np.zeros(9), np.zeros(9), np.linspace(-2, 2, 9)]
        ecto = np.c_[np.zeros(10), np.zeros(10), -2.0 - np.arange(1, 11) * 0.8]   # 8 nm rod below the outer leaflet
        cyto = np.c_[np.zeros(5), np.zeros(5), 3.0 + np.arange(5)]
        b = np.vstack([tmb, ecto, cyto]).astype(np.float32)
        np.savez(os.path.join(tmp, "U.npz"), beads=b, membrane=np.array("TM"), src=np.array("pdb:1aaa"))
        np.savez(os.path.join(tmp, "U0.npz"), beads=np.vstack([tmb, ecto[:5], cyto]).astype(np.float32), membrane=np.array("TM"),
                 src=np.array("pdb:2bbb"))
        g = CP.geom(os.path.join(tmp, "U.npz"))
        assert abs(g["ecto_h"] - 8.0) < 1e-5 and abs(g["ecto_depth"] - 8.0) < 1e-5 and g["beads"] == 24
        R = CP.compare([dict(ing_id="U", kind="complex")], {"U": os.path.join(tmp, "U.npz")}, {"U": os.path.join(tmp, "U0.npz")},
                       pm={"U"})
        assert set(R["flags"][0].split(";")) == {"source", "beads", "pm_height", "rg"}, R["flags"][0]
        assert abs(R.d_ecto_h[0] - 4.0) < 1e-5
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    for t in (test_decide_rules, test_assembly_filter_check, test_afdb_chain_identity, test_core_builders_filter,
              test_struct_src_from, test_atp_guard, test_compare_geom):
        t(); print("ok", t.__name__, flush=True)
