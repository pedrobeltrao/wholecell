"""Tests of the organisation layer (P6) and the packing shapes (P8) of the core packer (unification plan A7 phase 2),
about a minute on 1 CPU:

    PYTHONHASHSEED=0 python tests/test_pack_organisation.py        (or pytest tests/test_pack_organisation.py)

  switches      organisation and shapes default to v7 (off / library); config errors name the key
  network       PPR profiles: dense inverse == Neumann series == a pure-Python power iteration (rel 1e-10); nodiag
                rows sum to 1 with a zero seed entry; type profiles are standardised (Z Z^T = Pearson r)
  inputs        token map (first row wins), STRING edges (score cut, a < b, max per pair), subunits (uniprot field,
                else genes), pool edges with the paralog filter
  partners      contact before loose, hub down-weighting, shared-subunit pairs skipped, r_loose cut, max_partners;
                remap to another type order
  kernel        partner attraction on a synthetic cytosol tile of the testcell: copies are proposed next to a placed
                partner and end up closer to it than without the CSR
  shapes        pack.shapes trees replace the library trees; a stale tree (sig != manifest md5) and a missing tree of a
                species with beads are errors
  b2 overlays   config/cells/{jurkat,yeast}_b2.yaml differ from their cell configs only in pack.organisation (B2)
  field cap     (freeze fix 2026-11-19) in a listed region the field may not raise a tile's expected mass above the
                unmodulated split: after requests.field_cap every tile is at it (totals conserved), other regions and
                unmodulated types unchanged, the type contrast kept; the multinomial split (distribute) follows, the
                other regions' draws unchanged (own seed stream); config errors name the key
  normal test   (freeze fix 2026-11-19) pack3d.normal_ok: on a membrane point the gradient normal passes, its reverse
                fails (a sign flip by the validator's measure), a flat-gradient point fails for either normal
"""
import json
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "tests"))

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell.pack import organisation as O  # noqa: E402
from wholecell.pack import switches  # noqa: E402

LIBDIR = os.path.join(REPO, "tests", "data", "testcell", "lib")


# ------------------------------------------------------------------ switches
def test_switch_defaults():
    sw = switches.load({}, ["cytosol"])
    assert sw.legacy and sw.organisation["enabled"] is False and sw.shapes == {"source": "library", "file": None}
    assert not switches.load({"shapes": "pack_shapes"}).legacy
    on = switches.load({"organisation": {"enabled": True, "string": {"links": "x.txt.gz", "ids": {"aliases": "a.txt.gz"}}}})
    assert not on.legacy and on.organisation["string"]["min_score"] == 400 and on.organisation["string"]["ids"]["source"] == "UniProt_AC"
    for bad, words in (({"shapes": "pack_shape"}, "did you mean 'pack_shapes'"),
                       ({"organisation": {"enabled": True}}, "STRING links file"),
                       ({"organisation": {"enabled": True, "string": {"links": "x"}}}, "aliases"),
                       ({"organisation": {"dampin": 0.5}}, "did you mean 'damping'"),
                       ({"organisation": {"string": {"min_scor": 1}}}, "did you mean 'min_score'")):
        try:
            switches.load(bad)
        except ValueError as e:
            assert words in str(e), (bad, str(e))
        else:
            raise AssertionError(f"no error for {bad}")


# ------------------------------------------------------------------ network propagation
def _graph(n=40, m=120, seed=0):
    rng = np.random.default_rng(seed)
    a = rng.integers(0, n, m); b = rng.integers(0, n, m); k = a != b
    names = [f"P{i:03d}" for i in range(n)]
    d = pd.DataFrame(dict(a=[names[min(x, y)] for x, y in zip(a[k], b[k])], b=[names[max(x, y)] for x, y in zip(a[k], b[k])],
                          score=rng.integers(400, 1000, int(k.sum()))))
    return d.groupby(["a", "b"], as_index=False).score.max()


def _ppr_python(edges, d, nodes):
    """pure-Python reference: power iteration x <- (1 - d) e_s + d W x per seed, lists only"""
    ix = {u: i for i, u in enumerate(nodes)}; N = len(nodes)
    nb = [dict() for _ in range(N)]
    for a, b, s in zip(edges.a, edges.b, edges.score):
        nb[ix[a]][ix[b]] = s / 1000.0; nb[ix[b]][ix[a]] = s / 1000.0
    deg = [sum(v.values()) for v in nb]
    out = []
    for s in range(N):
        x = [0.0] * N; x[s] = 1.0 - d
        for _ in range(80):
            y = [0.0] * N; y[s] = 1.0 - d
            for j in range(N):                           # W[i, j] = A[i, j] / deg[j] (column-stochastic)
                if x[j]:
                    for i, w in nb[j].items():
                        y[i] += d * w / deg[j] * x[j]
            x = y
        out.append(x)
    return np.array(out)


def test_ppr():
    e = _graph()
    nodes, Xd, rd = O.ppr_profiles(e, 0.5, "raw")
    ref = _ppr_python(e, 0.5, nodes)
    assert np.allclose(Xd, ref, rtol=1e-10, atol=1e-14), np.abs(Xd - ref).max()
    nodes, Xd, rd = O.ppr_profiles(e, 0.5, "nodiag")
    need = nodes[::3]
    nodes2, Xn, rn = O.ppr_profiles(e, 0.5, "nodiag", need=need, dense_max=0)
    assert nodes2 == nodes and len(Xn) == len(need)
    for u in need:
        i = nodes.index(u)
        assert np.allclose(Xn[rn[i]], Xd[rd[i]], rtol=1e-10, atol=1e-15)
    assert np.allclose(Xd.sum(1), 1.0) and np.all(np.diag(Xd) == 0)
    subunits = {"T0": {nodes[0], nodes[1]}, "T1": {nodes[2]}, "T2": {"NOT_IN_NETWORK"}}
    tix, Z = O.type_profiles(["T0", "T1", "T2"], subunits, nodes, Xd, rd)
    assert tix.tolist() == [0, 1]
    r = np.corrcoef(np.stack([Xd[[0, 1]].mean(0), Xd[2]]))[0, 1]
    assert abs((Z @ Z.T)[0, 1] - r) < 1e-12
    E, lam = O.field_loadings(Z @ Z.T, 1)
    assert E.shape == (2, 1) and lam[0] >= 0


# ------------------------------------------------------------------ inputs
def test_inputs():
    up = pd.DataFrame({"Entry": ["U1", "U2", "U3"], "Gene Names (primary)": ["AAA1", "BBB1", "AAA1"],
                       "Gene Names (ordered locus)": ["YA01; YA01B", "YB02", None]})
    m = O.token_map(up, ["Gene Names (primary)", "Gene Names (ordered locus)"])
    assert m == {"AAA1": "U1", "YA01": "U1", "YA01B": "U1", "BBB1": "U2", "YB02": "U2"}
    tmp = tempfile.mkdtemp(prefix="wc_org_")
    f = os.path.join(tmp, "4932.protein.links.v12.0.txt")
    with open(f, "w") as fh:
        fh.write("protein1 protein2 combined_score\n4932.YA01 4932.YB02 500\n4932.YB02 4932.YA01 500\n"
                 "4932.YA01B 4932.YB02 700\n4932.YA01 4932.YZZZ 900\n4932.YB02 4932.YA01 300\n")
    loc = O.token_map(up, ["Gene Names (ordered locus)"])
    ed = O.string_edges(f, loc, 400)
    assert ed.values.tolist() == [["U1", "U2", 700]], ed
    tab = pd.DataFrame(dict(ing_id=["i1", "i2", "i3"], uniprot=["U1;U2", None, None], genes=[None, "AAA1;XXX", None]))
    su = O.lib_subunits([("i1", "L1"), ("i2", "L1"), ("i2", "L2"), ("i3", "L3")], tab, m)
    assert su == {"L1": {"U1", "U2"}, "L2": {"U1"}}
    pairs = pd.DataFrame(dict(a=["A", "A2", "A", "C"], b=["B", "B", "C", "D"], pae=[5.0, 3.0, 2.0, 9.0]))
    par = pd.DataFrame(dict(gene1_uniprot=["A"], gene2_uniprot=["A2"]))
    pe = O.pool_edges(pairs, par, 8.0)
    assert sorted(zip(pe.a, pe.b)) == [("A", "C"), ("A2", "B")]       # A-B dropped (paralog A2 scores better); C-D > 8


# ------------------------------------------------------------------ partners
def test_partners():
    ids = ["A", "B", "C", "D", "AB"]
    sub = {"A": {"pa"}, "B": {"pb"}, "C": {"pc"}, "D": {"pd"}, "AB": {"pa", "pb"}}
    pool = pd.DataFrame(dict(a=["pa"], b=["pb"], pae=[2.0]))
    tix = np.array([0, 1, 2, 3]); R = np.array([[1, .2, .9, .5], [.2, 1, .1, .3], [.9, .1, 1, .45], [.5, .3, .45, 1.0]])
    pr = O.partners(ids, sub, pool, tix, R, 8.0, 0.4, 16)
    rows = {ids[i]: [(ids[pr["idx"][k]], int(pr["kind"][k]), round(float(pr["strength"][k]), 6))
                     for k in range(pr["off"][i], pr["off"][i + 1])] for i in range(len(ids))}
    # A: contact B (PAE 2 -> 1, network r 0.2 -> (1 + 0.2) / 2), then loose C (0.9) before D (0.5); AB shares subunits
    assert rows["A"] == [("B", 0, 0.6), ("C", 1, round(0.2 + 0.8 * 0.5 / 0.6, 6)), ("D", 1, round(0.2 + 0.8 * 0.1 / 0.6, 6))], rows["A"]
    assert [x[0] for x in rows["B"]] == ["A"] and rows["AB"] == []
    assert [x[0] for x in rows["C"]] == ["A", "D"]
    pr1 = O.partners(ids, sub, pool, tix, R, 8.0, 0.4, 1)
    assert np.diff(pr1["off"]).max() == 1
    pr0 = O.partners(ids, sub, None, tix, R, 8.0, 0.4, 16)            # pool off: loose partners only
    assert (pr0["kind"] == O.LOOSE).all() and "B" not in [ids[j] for j in pr0["idx"][pr0["off"][0]:pr0["off"][1]]]
    same = O.remap_partners(pr, ids, ids)
    assert all(np.array_equal(same[k], pr[k]) for k in same)
    sub_ids = ["C", "A", "X"]
    rm = O.remap_partners(pr, ids, sub_ids)
    got = {sub_ids[i]: [sub_ids[j] for j in rm["idx"][rm["off"][i]:rm["off"][i + 1]]] for i in range(3)}
    assert got == {"C": ["A"], "A": ["C"], "X": []}, got
    assert O.params(switches.DEFAULTS["organisation"]).tolist() == [0.6, 0.25, 8.0, 0.6, 40, 2.0, 15.0]


# ------------------------------------------------------------------ kernel: partner attraction
def test_kernel_attraction():
    from test_pack_switches import _tile
    from wholecell.pack import pack3d as p3
    T = _tile(); cell, reg = T["cell"], T["reg"]
    e = np.zeros(0, np.int64)
    nv = reg.NV
    vox_off = np.array([0] * (T["cyto"] + 1) + [len(T["VOXP"])] * (nv - T["cyto"]) + [len(T["VOXP"])] * len(reg.MEM_REG), np.int64)
    jn = np.array([60, 200], np.int64); jl = np.arange(2, dtype=np.int64); jr = np.full(2, T["cyto"], np.int64)
    cap = int((np.diff(T["so"])[jl] * jn).sum()) + 1000

    def run(ptn):
        return p3.pack_tile(7, T["lo"], T["size"], 40.0, cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC,
                            T["so"], T["STC"], T["STR"], T["BR"], reg.regmask(), jl, jr, jn, nv, vox_off, T["VOXP"], T["vox"],
                            np.zeros(len(reg.MEM_REG) + 1, np.int64), np.zeros((0, 9)), np.zeros(1), np.zeros((0, 3)), np.zeros(0),
                            0, 1500.0, 0, 200, 25, cap, cap * 14, *ptn, np.zeros((0, 4)), reg.regmem(), *reg.shell_arrays(),
                            e, e, np.zeros((0, 5)), e, e, np.zeros(0), e)
    v7 = run((e, np.zeros((0, 3)), e, e, np.zeros(0), e, np.r_[np.zeros(7), 10.0, 16.0]))
    prm = np.r_[O.params(dict(switches.DEFAULTS["organisation"], p_max=1.0, a_contact=1.0, k_half=0.5)), 10.0, 16.0]
    att = run((e, np.zeros((0, 3)), np.array([0, 0, 1], np.int64), np.array([0], np.int64), np.array([1.0]),
               np.array([0], np.int64), prm))                # type 1 (ALDOA tetramer) attracted to type 0 (ACTB), contact

    def dist(r):
        t, p = r[0], r[1].astype(np.float64)
        a, b = p[t == 0], p[t == 1]
        return float(np.mean(np.min(np.linalg.norm(b[:, None] - a[None], axis=2), axis=1)))
    assert att[9][0] > 0 and att[9][1] > 0, att[9]
    assert v7[9].tolist() == [0, 0]
    assert dist(att) < 0.8 * dist(v7), (dist(att), dist(v7))
    print(f"  attraction: {att[9][1]} of {att[3][1]} copies placed next to a partner; mean nearest-partner distance "
          f"{dist(att):.1f} nm vs {dist(v7):.1f} nm without")


# ------------------------------------------------------------------ packing shapes
def test_shapes():
    import shutil
    from wholecell.library import pack_shapes as PS
    from wholecell.pack import stage
    tmp = tempfile.mkdtemp(prefix="wc_shapes_")
    try:
        libs = ["ACTB", "ALDOA_x4"]
        man = os.path.join(tmp, "library"); os.makedirs(man)
        rows = [(i, os.path.join(LIBDIR, i + ".npz")) for i in libs]
        md5 = {i: PS._md5(p) for i, p in rows}
        pd.DataFrame(dict(species=libs, path=[p for _, p in rows], md5=[md5[i] for i in libs])).to_csv(
            os.path.join(man, "library_manifest.tsv"), sep="\t", index=False)
        PS.run(dict(rows), os.path.join(man, "pack_shapes"), per=15, nproc=1, log=lambda *a: None)

        class _P:
            bdir = tmp

            def library(self):
                return dict(rows)
        sw = switches.load({"shapes": "pack_shapes"})
        sh = stage._shapes(_P(), sw, libs)
        z = np.load(os.path.join(man, "pack_shapes", "pack_trees_15.npz"))
        for k, i in enumerate(z["ids"]):
            assert np.array_equal(sh[str(i)][1], z["R"][z["off"][k]:z["off"][k + 1]])
        so, STC, STR, BR = stage._libpack(libs, dict(rows), sh)
        so0, _, STR0, BR0 = stage._libpack(libs, dict(rows))
        assert len(STR) != len(STR0) and np.all(BR > 0)       # finer trees (15 residues per sphere)
        assert stage._shapes(_P(), switches.load({}), libs) is None
        pd.DataFrame(dict(species=libs, path=[p for _, p in rows], md5=["x", md5["ALDOA_x4"]])).to_csv(
            os.path.join(man, "library_manifest.tsv"), sep="\t", index=False)
        try:
            stage._shapes(_P(), sw, libs)
        except ValueError as e:
            assert "another library file" in str(e)
        else:
            raise AssertionError("stale tree not detected")
        rows.append(("ACO2", os.path.join(LIBDIR, "ACO2.npz")))
        try:
            stage._shapes(_P(), sw, libs + ["ACO2"])
        except ValueError as e:
            assert "have no tree" in str(e)
        else:
            raise AssertionError("missing tree not detected")
        # pack.shapes.library (yeast: engine/lib): trees made from another directory (engine wc3d.pack_trees signature
        # size:mtime), including a species that is not in the manifest (an assembly product); checked against that dir
        rows.pop()
        eng = os.path.join(tmp, "engine_lib"); os.makedirs(eng)
        for i, p_ in rows + [("ACTB_oligo", os.path.join(LIBDIR, "ACTB.npz"))]:
            shutil.copy2(p_, os.path.join(eng, i + ".npz"))
        ids2 = libs + ["ACTB_oligo"]
        tr = PS.pack_trees({i: os.path.join(eng, i + ".npz") for i in ids2}, 15, nproc=1)
        off = np.r_[0, np.cumsum([len(tr[i][1]) for i in ids2])]

        def sm(i):
            st = os.stat(os.path.join(eng, i + ".npz")); return f"{st.st_size}:{st.st_mtime:.0f}"
        tf = os.path.join(tmp, "engine_trees.npz")
        np.savez(tf, ids=np.array(ids2), sig=np.array([sm(i) for i in ids2]), off=off,
                 C=np.concatenate([tr[i][0] for i in ids2]), R=np.concatenate([tr[i][1] for i in ids2]))
        sw2 = switches.load({"shapes": {"source": "pack_shapes", "file": tf, "library": eng}})
        # without pack.shapes.library a size:mtime signature is checked against the manifest file itself (copy2 kept
        # size and mtime: same file); the assembly product has no manifest entry, so no path for the packer
        stage._shapes(_P(), switches.load({"shapes": {"source": "pack_shapes", "file": tf}}), libs)
        assert "ACTB_oligo" not in stage._shapes_library(_P(), switches.load({"shapes": {"source": "pack_shapes", "file": tf}}))[0]
        sh2 = stage._shapes(_P(), sw2, ids2)
        assert sorted(sh2) == sorted(ids2) and stage._shapes_library(_P(), sw2)[0]["ACTB_oligo"].startswith(eng)
        os.utime(os.path.join(eng, "ACTB.npz"), (1.0e9, 1.0e9))     # the engine rewrote a library entry after the trees
        try:
            stage._shapes(_P(), sw2, ids2)
        except ValueError as e:
            assert "1 trees were made from another library file (e.g. ['ACTB']" in str(e), e
        else:
            raise AssertionError("stale engine tree not detected")
        print("  shapes: md5 and size:mtime signatures, pack.shapes.library overlay, stale trees refused")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ------------------------------------------------------------------ B2 in the cell configs (criteria docs/unification_criteria/B2_organisation.md)
def test_b2_cell_configs():
    """B2 settings kept in jurkat / yeast (STRING v12 >= 400, AF3 pool and paralogs off) but the layer is OFF since the
    v9/u2 rebuild (2026-11-23); the reference configs keep theirs (jurkat_v7_port off, yeast_x5_core on with the pool)"""
    from wholecell.config import load_raw
    for c in ("jurkat", "yeast"):
        p = switches.effective(load_raw(c)["pack"])
        o = switches.load(p, [r["name"] for r in p["regions"]]).organisation
        assert not o["enabled"] and o["pool_scores"] is None and o["paralogs"] is None and o["string"]["min_score"] == 400, c
        assert str(o["string"]["version"]) == "12.0" and o["string"]["taxon"] == (9606 if c == "jurkat" else 4932), c
    j = load_raw("jurkat")
    assert j["pack"]["organisation"]["string"]["links"] in j["frozen_inputs"]
    assert j["pack"]["organisation"]["string"]["ids"]["table"] in j["frozen_inputs"]
    assert not (load_raw("jurkat_v7_port")["pack"].get("organisation") or {}).get("enabled")      # v7 port: off
    yo = load_raw("yeast_x5_core")["pack"]["organisation"]                                       # X5 setting: on, pool on
    assert yo["enabled"] and yo["pool_scores"] and yo["paralogs"]
    print("  B2: organisation off in jurkat / yeast (settings kept); reference configs keep their settings")


# ------------------------------------------------------------------ freeze fix: field cap in crowded regions
def _cap_case(seed=0, nt=60, ntype=8):
    from wholecell.pack.requests import crowding_factors
    rng = np.random.default_rng(seed)
    M = np.zeros((nt, 2)); M[:, 0] = rng.uniform(0.2, 1.0, nt); M[: nt // 2, 1] = 1.0      # region 0: crowded, 1: other
    libs = list(range(ntype))
    logmod = {l: rng.normal(0, 0.8, nt) for l in libs[:-1]}                                 # last type: unmodulated
    vol = {l: float(rng.uniform(1e3, 1e5)) for l in libs}
    mass = {l: vol[l] * float(rng.uniform(5, 20)) for l in libs}                            # Da per sphere volume varies
    e = pd.DataFrame([(l, r, float(rng.uniform(1e3, 1e4))) for l in libs for r in ("crowded", "other")],
                     columns=["lib_id", "region", "expected"])
    ri = {"crowded": 0, "other": 1}
    return M, libs, logmod, vol, mass, e, ri, crowding_factors(e, M, ri, logmod, vol)


def _tile_mass(e, M, logmod, cf, mass, reg, ri):
    W = M[:, ri[reg]]; U = np.zeros(len(W)); per = {}
    for r in e[e.region == reg].itertuples():
        w = W * np.exp(logmod[r.lib_id] - logmod[r.lib_id][W > 0].max()) * cf[ri[reg]] if r.lib_id in logmod else W.copy()
        per[r.lib_id] = r.expected * w / w.sum(); U += per[r.lib_id] * mass[r.lib_id]
    return U, W / W.sum() * U.sum(), per


def test_field_cap():
    from wholecell.pack.requests import field_cap, distribute
    M, libs, logmod, vol, mass, e, ri, cf = _cap_case()
    U0, R, per0 = _tile_mass(e, M, logmod, cf, mass, "crowded", ri)
    assert U0.max() / R[U0.argmax()] > 1.05                     # the volume crowding correction leaves mass moved
    cf2, rep = field_cap(e, M, ri, logmod, cf, mass, ["crowded"])
    U1, R1, per1 = _tile_mass(e, M, logmod, cf2, mass, "crowded", ri)
    assert np.allclose(R1, R) and np.all(U1 <= R * (1 + 1e-5)) and np.allclose(U1, R, rtol=1e-5), (U1 / R).max()
    assert rep["crowded"]["max_ratio_before"] > 1.05 and rep["crowded"]["max_ratio_after"] <= 1 + 1e-5, rep
    assert np.array_equal(cf2[1], cf[1])                        # other regions untouched
    l_un = libs[-1]
    assert np.allclose(per1[l_un], per0[l_un])                  # unmodulated types keep the measure split
    for l in libs[:-1]:                                         # totals conserved, field contrast kept (same sign)
        assert np.isclose(per1[l].sum(), per0[l].sum())
        nz = M[:, 0] > 0
        assert np.corrcoef(np.log(per1[l][nz] / M[nz, 0]), logmod[l][nz])[0, 1] > 0.3
    same, _ = field_cap(e, M, ri, logmod, cf, mass, [])
    assert all(np.array_equal(same[k], cf[k]) for k in cf)
    # the multinomial split with the cap: per-tile mass of a large draw at the unmodulated split
    big = e.assign(expected=e.expected * 200, region=e.region.map({"crowded": "crowded", "other": "other"}))
    rep2 = {}
    d = distribute(big, M, ["crowded", "other"], {l: l for l in libs}, 5, logmod=logmod, vol=vol, cap=(mass, ["crowded"]), report=rep2)
    d0 = distribute(big, M, ["crowded", "other"], {l: l for l in libs}, 5, logmod=logmod, vol=vol)
    def mt(dd):
        x = dd[dd.region == 0]; return np.bincount(x.tile, weights=x.n * x.lib.map(mass), minlength=len(M))
    Rb = M[:, 0] / M[:, 0].sum() * mt(d).sum()
    assert np.abs(mt(d) / Rb - 1).max() < 0.03 < np.abs(mt(d0) / Rb - 1).max(), (np.abs(mt(d) / Rb - 1).max(), np.abs(mt(d0) / Rb - 1).max())
    tot1, tot0 = d.groupby(["lib", "region"]).n.sum(), d0.groupby(["lib", "region"]).n.sum()
    assert "crowded" in rep2["field_cap"] and (tot1 - tot0).abs().max() == 0     # same totals (rounding stream)
    o1 = d[d.region == 1].sort_values(["tile", "lib"]).reset_index(drop=True)
    o0 = d0[d0.region == 1].sort_values(["tile", "lib"]).reset_index(drop=True)
    assert o1.equals(o0)                                        # other regions: the uncapped draws (own cap stream)
    # config
    assert switches.load({}, ["cytosol"]).organisation["field_cap_regions"] == []
    assert switches.load({"organisation": {"field_cap_regions": ["nucleolus"]}}, ["cytosol", "nucleolus"]).organisation["field_cap_regions"] == ["nucleolus"]
    for bad, words in (({"organisation": {"field_cap_regions": ["nucleolos"]}}, "did you mean 'nucleolus'"),
                       ({"organisation": {"field_cap_regions": "nucleolus"}}, "expected a list")):
        try:
            switches.load(bad, ["cytosol", "nucleolus"])
        except ValueError as ex:
            assert words in str(ex) and "field_cap_regions" in str(ex), str(ex)
        else:
            raise AssertionError(f"no error for {bad}")
    print(f"  field cap: max tile mass / unmodulated {rep['crowded']['max_ratio_before']:.3f} -> "
          f"{rep['crowded']['max_ratio_after']:.6f} ({rep['crowded']['iterations']} iterations); multinomial split "
          f"max |dev| {np.abs(mt(d) / Rb - 1).max():.3f} vs {np.abs(mt(d0) / Rb - 1).max():.3f} without")


# ------------------------------------------------------------------ freeze fix: membrane normal test of partner proposals
def test_normal_ok():
    from wholecell.config import load
    from wholecell.pack import geomload
    from wholecell.pack import pack3d as p3
    from wholecell.geometry import engine as c3
    g = geomload.load(load("testcell")); cell = g.cell
    out = np.zeros(c3.NMEM_MAX); n = np.zeros(3)
    rng = np.random.default_rng(1)
    R = float(cell.G[0]); m = 0                                  # PM
    ok = flips = 0
    S = rng.uniform(-1.1 * R, 1.1 * R, (20000, 3))
    S = S[np.abs(cell.mems(S)[:, m]) < 20.0][:60]                # starts within 20 nm of the membrane
    for p in S:
        p = p.astype(np.float64)
        if not p3.project_to_membrane(p, m, n, cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC, out):
            continue
        a = p3.normal_ok(p, m, n, cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC, out)
        b = p3.normal_ok(p, m, -n, cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC, out)
        ok += a; flips += b
        assert a != b or not a, "a normal and its reverse both pass"
    assert ok >= 10 and flips == 0, (ok, flips)
    # flat gradient: points where the validator's |grad| of membrane m is < GRAD_MIN fail for any normal
    P = rng.uniform(-1.1 * R, 1.1 * R, (20000, 3)); h = p3.GRAD_H
    G = np.stack([(cell.mems(P + h * np.eye(3)[k])[:, m] - cell.mems(P - h * np.eye(3)[k])[:, m]) / (2 * h) for k in range(3)], 1)
    flat = P[np.linalg.norm(G, axis=1) < 0.4]
    if not len(flat):
        flat = np.zeros((1, 3))                                  # centre of the (convex) cell: the PM SDF has a ridge there
    for q in flat[:20]:
        for s in (1.0, -1.0):
            nn = s * np.array([0.0, 0.0, 1.0])
            assert not p3.normal_ok(q.astype(np.float64), m, nn, cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC, out)
    print(f"  normal test: {ok} membrane points pass with their normal, 0 with the reverse; {min(len(flat), 20)} flat-gradient points rejected")


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("pack organisation / shapes tests OK")
