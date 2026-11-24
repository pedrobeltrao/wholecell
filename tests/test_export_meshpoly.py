"""Tests of the YCT2 MESH / GHOST / POLY sections (wholecell/export/meshpoly.py, yct2.py), synthetic data, a few seconds:

    PYTHONHASHSEED=0 python tests/test_export_meshpoly.py        (or pytest tests/test_export_meshpoly.py)

  mesh_roundtrip       marching-cubes mesh of a sphere membrane in one tile: encode -> decode within the position
                       quantisation, normals to positive SDF, area ~ 4 pi R^2, u32 triangle indices when nv > 65535
  poly_roundtrip       polylines (mRNA / DNA / extra with kinds): quantise -> pack -> unpack exact to q/2
  tile_sections        a tile with ghosts, three POLY sections and a MESH decodes back; PALETTE prefix and INST bytes are
                       identical to the same tile without the extra sections (no regression of the instance content)
  determinism          the same inputs give the same raw and gzip bytes twice
  sources_and_store    thread resampling (3 nm beads -> 0.6 nm steps), linker DNA between two nucleosomes, PolyStore
                       sorted per tile: the per-tile sections decode back to the source points
"""
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell.export import meshpoly as MP, quat as Q, yct2  # noqa: E402

TILE, QO, QR = 200.0, 80.0, 360.0


class SphereCell:
    """two membranes: code 0 = sphere (SDF positive outside), code 1 = plane z = c (never crossed in the tile)"""

    def __init__(self, c, R):
        self.c, self.R = np.asarray(c, float), float(R)

    def mems(self, P):
        P = np.asarray(P, float)
        return np.stack([np.linalg.norm(P - self.c, axis=1) - self.R, P[:, 2] - 1e4], 1).astype(np.float32)


def _instances(rng, n, lo):
    sp = rng.integers(0, 5, n); rg = rng.integers(0, 3, n)
    pos = lo + rng.random((n, 3)) * TILE
    q = rng.normal(size=(n, 4)); q /= np.linalg.norm(q, axis=1, keepdims=True); q[q[:, 0] < 0] *= -1
    return sp, rg, np.zeros(n, np.uint32), pos, q


def _mesh(lo):
    cell = SphereCell(lo + 100.0, 70.0)
    return cell, MP.TileMesher(cell, [(0, "S"), (1, "P")], TILE)(lo)


def test_mesh_roundtrip():
    lo = np.array([200.0, -400.0, 600.0])
    cell, (V, F, FM) = _mesh(lo)
    assert len(F) > 1000 and set(FM.tolist()) == {0}
    b, nv, nt, big = MP.encode_mesh(V, F, FM, lo, QO, QR)
    assert not big and len(b) % 4 == 0
    V2, F2, FM2 = MP.decode_mesh(b, nt, big, lo, QO, QR)
    assert np.array_equal(F2, F) and np.array_equal(FM2, FM) and len(V2) == nv
    assert np.abs(V2 - V).max() <= QR / 65535 / 2 + 1e-9
    tri = V2[F2]; n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    c = tri.mean(1)
    assert (np.einsum("ij,ij->i", n, c - cell.c) > 0).mean() > 0.999          # normals to positive SDF (outside)
    area = 0.5 * np.linalg.norm(n, axis=1).sum()
    assert abs(area / (4 * np.pi * 70.0 ** 2) - 1) < 0.03, area
    sd = np.abs(cell.mems(V2)[:, 0])
    assert np.median(sd) < 0.5 and sd.max() < 2.0
    # u32 indices when the vertex count exceeds 65535
    Vb = np.repeat(V, 1 + 70000 // len(V), axis=0)[:70000]
    b2, nv2, nt2, big2 = MP.encode_mesh(Vb, F, FM, lo, QO, QR)
    assert big2 and nv2 == 70000
    V3, F3, _ = MP.decode_mesh(b2, nt2, True, lo, QO, QR)
    assert len(V3) == 70000 and np.array_equal(F3, F)


def _chains(rng, lo, n, step, kinds=None):
    out = []
    for _ in range(n):
        m = int(rng.integers(2, 64))
        d = rng.normal(size=(m - 1, 3)); d = d / np.linalg.norm(d, axis=1, keepdims=True) * step
        p0 = lo + 20 + rng.random(3) * 160
        out.append(np.vstack([p0, p0 + np.cumsum(d, 0)]))
    return out


def test_poly_roundtrip():
    rng = np.random.default_rng(1); lo = np.array([0.0, 200.0, -200.0])
    for sec in (4, 5, 6):
        q = MP.SECTIONS[sec]["poly_q"]
        ch = _chains(rng, lo, 30, MP.SECTIONS[sec]["step"])
        lens, st, d = MP.quantise_chains(ch, lo, q)
        kinds = rng.integers(0, 6, len(ch)) if sec == 6 else None
        b = MP.pack_polylines(lens, st, d, kinds)
        assert len(b) % 4 == 0
        res = MP.unpack_polylines(b, len(ch), lo, q, extra=sec == 6)
        got = res[0] if sec == 6 else res
        if sec == 6:
            assert np.array_equal(res[1], kinds)
        for a, g in zip(ch, got):
            assert a.shape == g.shape and np.abs(a - g).max() <= q / 2 + 1e-9


def _tile(extra=True, seed=2):
    rng = np.random.default_rng(seed); lo = np.array([400.0, 0.0, -200.0]); ijk = [2, 0, -1]
    sp, rg, fl, pos, q = _instances(rng, 500, lo)
    kw = {}
    if extra:
        g = _instances(rng, 40, lo - 15)
        g = (g[0] + 3, g[1], g[2] | yct2.FLAG_GHOST | yct2.FLAG_NEAR_MEMBRANE, g[3], g[4])    # species 3..7: new palette entries
        polys = {}
        for sec in (4, 5, 6):
            ch = _chains(rng, lo, 12, MP.SECTIONS[sec]["step"])
            lens, st, d = MP.quantise_chains(ch, lo, MP.SECTIONS[sec]["poly_q"])
            polys[sec] = (MP.pack_polylines(lens, st, d, np.full(len(ch), 5) if sec == 6 else None), len(ch))
        _, (V, F, FM) = _mesh(lo)
        mb, nv, nt, big = MP.encode_mesh(V, F, FM, lo, QO, QR)
        kw = dict(ghost=g, polys=polys, mesh=(mb, big, nt))
    raw, n_i, n_p = yct2.encode_tile(ijk, lo, sp, rg, fl, pos, q, TILE, QO, QR, 123, 456, **kw)
    return raw, lo, kw


def _secs(raw):
    h = np.frombuffer(raw[:128], "<u4"); ns = int(h[12])
    tab = np.frombuffer(raw[128:128 + 16 * ns], "<u4").reshape(ns, 4)
    return h, {int(k): raw[o:o + l] for k, o, l, c in tab}


def test_tile_sections():
    raw0, lo, _ = _tile(extra=False)
    raw1, lo, kw = _tile(extra=True)
    h0, s0 = _secs(raw0); h1, s1 = _secs(raw1)
    assert s1[yct2.K_INST] == s0[yct2.K_INST]                                   # instance words unchanged
    assert s1[yct2.K_PALETTE][:len(s0[yct2.K_PALETTE])] == s0[yct2.K_PALETTE]    # own runs first, ghost entries appended
    assert int(h1[11]) > int(h0[11]) and int(h1[10]) == 40
    assert all(len(b) % 4 == 0 for b in s1.values()) and all(int(o) % 4 == 0 for o in np.frombuffer(raw1[128:128 + 16 * int(h1[12])], "<u4")[1::4])
    D = yct2.decode_tile(raw1)
    gs, gr, gf, gp, gq = kw["ghost"]
    o = np.lexsort((gr, gs))
    G = D["ghost"]
    assert np.array_equal(G["species"], gs[o]) and np.array_equal(G["region"], gr[o]) and np.array_equal(G["flags"], gf[o])
    assert np.abs(lo + G["pos_rel"] - gp[o]).max() <= QR / 131070 + 1e-9
    assert Q.angle_between(G["quat"], gq[o]).max() < 1e-4
    for sec in (4, 5, 6):
        assert D["sections"][sec] == kw["polys"][sec][0] and D["counts"][sec] == kw["polys"][sec][1]
    nt = D["counts"][yct2.K_MESH]
    V, F, FM = MP.decode_mesh(D["sections"][yct2.K_MESH], nt, D["mesh_u32"], lo, QO, QR)
    assert len(F) == kw["mesh"][2] and len(V) > 0


def test_determinism():
    a, _, _ = _tile(); b, _, _ = _tile()
    assert a == b and yct2.gz(a) == yct2.gz(b)
    lo = np.array([0.0, 0.0, 0.0])
    assert MP.encode_mesh(*_mesh(lo)[1], lo, QO, QR)[0] == MP.encode_mesh(*_mesh(lo)[1], lo, QO, QR)[0]


def test_sources_and_store():
    rng = np.random.default_rng(3)
    # two threads of 3 nm beads (second one with a 50 nm gap: split there)
    t1 = np.cumsum(np.r_[[[10.0, 10.0, 10.0]], np.tile([[3.0, 0.0, 0.0]], (150, 1))], 0)
    t2 = np.cumsum(np.r_[[[50.0, 300.0, 50.0]], np.tile([[0.0, 3.0, 0.0]], (40, 1)), [[0.0, 50.0, 0.0]],
                         np.tile([[0.0, 0.0, 3.0]], (30, 1))], 0)
    pts = np.vstack([t1, t2]).astype(np.float32); off = np.array([0, len(t1), len(t1) + len(t2)])
    src = MP.ThreadSource("mrna", 4, pts, off)
    P, ln = zip(*src.chunks()); P = np.concatenate(P); ln = np.concatenate(ln)
    assert ln.max() <= MP.PIECE and ln.min() >= 2
    first = np.r_[0, np.cumsum(ln)[:-1]]
    st = np.linalg.norm(np.diff(P, axis=0), axis=1)
    inner = np.ones(len(P) - 1, bool); inner[first[1:] - 1] = False
    assert st[inner].max() <= 0.6 + 1e-6 and st[inner].min() > 0.59         # resampled to 0.6 nm, never across a gap
    assert abs(len(P) * 0.6 - (450 + 120 + 90)) < 10
    # linker DNA between two nucleosomes 20 nm apart (toy library: two antiparallel strands along a helix)
    th = np.linspace(0, 3.5 * np.pi, 40)
    I = np.c_[4.2 * np.cos(th), 4.2 * np.sin(th), np.linspace(-2, 2, 40)]; J = (I + [0, 0, 0.9])[::-1]
    beads = np.vstack([np.zeros((5, 3)), I, J]); kind = np.r_[np.zeros(5, int), np.ones(80, int)]
    ends = MP.nucleosome_ends(beads, kind)
    pos = np.array([[100.0, 100.0, 100.0], [120.0, 100.0, 100.0], [300.0, 300.0, 300.0]])
    qq = np.tile([[1.0, 0, 0, 0]], (3, 1))
    lk = MP.LinkerSource("linker", pos, qq, np.array([0, 0, 0]), np.array([0, 0, 0]), ends)
    LP, Lln = zip(*lk.chunks()); LP = np.concatenate(LP); Lln = np.concatenate(Lln)
    assert len(Lln) == 1                                                         # the 3rd nucleosome is > 40 nm away
    a = pos[0] + ends[1]; b = pos[1] + ends[0]
    assert np.allclose(LP[0], a) and np.allclose(LP[-1], b)
    assert np.linalg.norm(np.diff(LP, axis=0), axis=1).max() < 1.27
    # store: two tiles' worth, sorted per tile, decodes back
    with tempfile.TemporaryDirectory() as d:
        store = MP.PolyStore(d, [0.0, 0.0, 0.0], TILE, [0, 0, 0], [3, 3, 3])
        store.build([src, lk], log=lambda *a: None)
        store.open()
        tiles = store.tile_keys()
        assert tiles == {tuple(int(x) for x in np.floor(p[len(p) // 2] / TILE)) for p in np.split(P, np.cumsum(ln)[:-1])} | {(0, 0, 0)}
        got = []
        for t in sorted(tiles):
            secs, npts = store.tile_sections(t)
            lo = np.array(t) * TILE
            for sec, (bts, n) in secs.items():
                got += MP.unpack_polylines(bts, n, lo, MP.SECTIONS[sec]["poly_q"])
        allp = np.concatenate(got)
        assert len(allp) == len(P) + len(LP)
        ref = np.concatenate([P, LP])
        k1 = np.lexsort(np.round(allp, 1).T); k2 = np.lexsort(np.round(ref, 1).T)
        assert np.abs(allp[k1] - ref[k2]).max() <= 0.006


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("export mesh/poly tests OK")
