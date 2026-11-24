"""YCT2 MESH, GHOST and POLY sections (docs/reference/tile_format_v2.md section 4) and their build-side sources.

MESH  per viewer tile, marching cubes (MESH_H = 8 nm) of every lipid membrane of the build's geometry package (membrane
      SDFs of cellgeom_eval_v1, `cell.mems`; membranes with leaflets "bilayer"/"monolayer"), on the tile's own grid
      (lo .. lo + tile_nm, so neighbouring tiles share their edge vertices). Triangles are wound so that the face normal
      (b - a) x (c - a) points to positive SDF (the cytoplasmic / N side; the viewer reads it so). Triangles whose centroid
      lies below an open face of the imaged volume (G_FACE / G_ZFACE) are dropped. Layout: u16 vert[nv*3] (instance
      position quantisation), pad4; u16 tri[nt*3] (u32 when nv > 65535, header flag bit 0), pad4; u16 mem[nt] (membrane
      engine code), pad4. Section count = nt; nv = floor((len - tri bytes - pad4(2 nt)) / 6), and index.json n_vert.
POLY  u32 len[n]; i32 start[n*3] (poly_q units from the tile origin lo); i8 delta[(sum len - n)*3]; pad4. POLY_EXTRA is
      prefixed by u16 kind[n] (pad4). Section count = n polylines. Pieces of <= PIECE points belong to the tile of their
      middle point. Threads whose steps do not fit int8 deltas are resampled along their arc length.
GHOST instances of the 26 neighbours in the lipid-mask regions (membrane regions + config) whose bounding sphere + 3 nm
      reaches the tile: INST encoding, palette entries appended after the tile's own entries, flags ghost | near_membrane.

origin: tile meshes from jurkat/viewer/jenc.py (MESH_H, marching cubes per tile, clip_mesh) and jenc_full.py tile_mesh
(9^3 pre-test); polylines (pack_polylines, linker_chains, piece cutting of polymer_lists) from the yeast engine
plugins/yeast/engine/pipeline/webenc.py and scripts/s5_encode.py (yeast_cell_dev_bundle @ e3b7448), re-implemented
vectorised.
"""
from __future__ import annotations

import json
import os

import numpy as np

MESH_H = 8.0
PIECE = 64                      # points per polyline piece (yeast MRNA_PIECE)
GAP_NM = 8.0                    # a step longer than this splits a thread (no bead bridges a gap)
LINKER_STEP = 0.68              # nm between linker DNA points (2 bp per bead)
LINKER_MAX = 40.0               # nm: longer nucleosome gaps get no linker (yeast webenc.linker_chains)
SECTIONS = {4: dict(name="mRNA", poly_q=0.005, step=0.6, free_kind=0),
            5: dict(name="DNA", poly_q=0.01, step=0.68, free_kind=1),
            6: dict(name="extra", poly_q=0.005, step=0.6, free_kind=None)}
EXTRA_KINDS = {0: "nascent RNA", 4: "pre-rRNA", 5: "FG repeats (disordered protein)"}


def pad4(b: bytes) -> bytes:
    return b + b"\0" * ((-len(b)) % 4)


# ------------------------------------------------------------------ MESH
def lipid_membranes(manifest, skip=()):
    """[(engine_code, name)] of the package membranes that carry lipids (leaflets bilayer / monolayer)"""
    out = []
    for m in sorted(manifest["membranes"], key=lambda x: x["engine_code"]):
        lf = str(m.get("leaflets", "bilayer")).lower()
        if m["name"] in skip or not (lf.startswith("bilayer") or lf.startswith("mono")):
            continue
        out.append((int(m["engine_code"]), m["name"]))
    return out


def clip_mesh(V, F, FM, zf):
    """drop triangles whose centroid lies below the open face (z < zf) and the vertices no triangle uses any more"""
    if zf is None or not len(F):
        return V, F, FM
    keep = V[F].mean(1)[:, 2] >= zf
    F, FM = F[keep], FM[keep]
    used = np.zeros(len(V), bool); used[F.ravel()] = True
    remap = np.cumsum(used) - 1
    return V[used], remap[F].astype(F.dtype), FM


class TileMesher:
    """membrane meshes of one viewer tile from the geometry package (cell.mems = membrane SDFs, columns = engine codes)"""

    def __init__(self, cell, membranes, tile_nm, h=MESH_H, zf=None):
        self.cell, self.mem, self.tile_nm, self.h, self.zf = cell, list(membranes), float(tile_nm), float(h), zf
        self.codes = np.array([c for c, _ in self.mem], np.int64)
        n = int(round(self.tile_nm / self.h))
        if abs(n * self.h - self.tile_nm) > 1e-9:
            raise ValueError(f"mesh_h_nm {h} must divide tile_nm {tile_nm}")
        self.g = np.arange(n + 1) * self.h
        self.g9 = np.linspace(0.0, self.tile_nm, 9)

    def _grid(self, lo, g):
        X = np.stack(np.meshgrid(lo[0] + g, lo[1] + g, lo[2] + g, indexing="ij"), -1).reshape(-1, 3)
        return np.asarray(self.cell.mems(X), np.float64)[:, self.codes]

    def __call__(self, lo):
        """-> V (nv,3) float64 model nm, F (nt,3) int64, FM (nt,) uint16 engine code"""
        from skimage.measure import marching_cubes
        lo = np.asarray(lo, np.float64)
        empty = (np.zeros((0, 3)), np.zeros((0, 3), np.int64), np.zeros(0, np.uint16))
        if not len(self.codes):
            return empty
        m9 = self._grid(lo, self.g9)
        near = ((m9.min(0) < 0) & (m9.max(0) > 0)) | (np.abs(m9).min(0) < 30.0)
        if not near.any():
            return empty
        m_ = len(self.g)
        S = self._grid(lo, self.g).reshape(m_, m_, m_, len(self.codes))
        Vs, Fs, FM, nv = [], [], [], 0
        for j in np.flatnonzero(near):
            s = np.clip(S[..., j], -1e3, 1e3)
            if not (s.min() < 0 < s.max()):
                continue
            v, f, _, _ = marching_cubes(s, 0.0, spacing=(self.h,) * 3)
            v = v + lo; f = f.astype(np.int64)
            f = self._orient(v, f, int(self.codes[j]))
            Vs.append(v); Fs.append(f + nv); FM.append(np.full(len(f), self.codes[j], np.uint16)); nv += len(v)
        if not Vs:
            return empty
        V, F, FMm = np.concatenate(Vs), np.concatenate(Fs), np.concatenate(FM)
        return clip_mesh(V, F, FMm, self.zf)

    def _orient(self, v, f, code):
        """wind triangles so that (b - a) x (c - a) points to positive SDF: test the SDF 1 nm along the normals of up to 32
        faces (deterministic choice); flip all if most point to negative SDF"""
        if not len(f):
            return f
        n = np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]])
        ln = np.linalg.norm(n, axis=1)
        ok = np.flatnonzero(ln > 1e-6)
        if not len(ok):
            return f
        pick = ok[np.linspace(0, len(ok) - 1, min(32, len(ok))).astype(np.int64)]
        c = v[f[pick]].mean(1)
        u = n[pick] / ln[pick, None]
        sp = np.asarray(self.cell.mems(c + u), np.float64)[:, code]
        sm = np.asarray(self.cell.mems(c - u), np.float64)[:, code]
        if np.sum(sp > sm) * 2 < len(pick):
            f = f[:, [0, 2, 1]]
        return f


def encode_mesh(V, F, FM, lo, q_org, q_range):
    """-> (section bytes, nv, nt, big)"""
    nv, nt = len(V), len(F)
    if not nt:
        return b"", 0, 0, False
    u = np.clip(np.rint((np.asarray(V, np.float64) - (np.asarray(lo, np.float64) - q_org)) / q_range * 65535.0), 0, 65535)
    big = nv > 65535
    b = pad4(u.astype("<u2").tobytes()) + pad4(np.asarray(F).astype("<u4" if big else "<u2").tobytes()) \
        + pad4(np.asarray(FM).astype("<u2").tobytes())
    return b, nv, nt, big


def decode_mesh(sec: bytes, nt: int, big: bool, lo, q_org, q_range):
    """inverse of encode_mesh -> V (model nm), F, FM"""
    tb = nt * 12 if big else len(pad4(b"\0" * (6 * nt)))
    mb = len(pad4(b"\0" * (2 * nt)))
    nv = (len(sec) - tb - mb) // 6
    o = 0
    u = np.frombuffer(sec, "<u2", nv * 3, o).reshape(-1, 3).astype(np.float64); o += len(pad4(b"\0" * (6 * nv)))
    F = np.frombuffer(sec, "<u4" if big else "<u2", nt * 3, o).reshape(-1, 3).astype(np.int64); o += tb
    FM = np.frombuffer(sec, "<u2", nt, o).copy()
    V = np.asarray(lo, np.float64) - q_org + u * (q_range / 65535.0)
    return V, F, FM


# ------------------------------------------------------------------ POLY
def pack_polylines(lens, starts, deltas, kinds=None) -> bytes:
    """pieces already quantised: lens (n,), starts (n,3) int, deltas (sum lens - n, 3) int8; kinds (n,) for POLY_EXTRA"""
    if not len(lens):
        return b""
    b = np.asarray(lens, "<u4").tobytes() + np.asarray(starts, "<i4").tobytes() + np.asarray(deltas, "i1").tobytes()
    if kinds is not None:
        b = pad4(np.asarray(kinds, "<u2").tobytes()) + b
    return pad4(b)


def quantise_chains(chains, lo, q):
    """list of (m,3) point arrays -> (lens, starts, deltas) relative to lo in q units (exact delta coding)"""
    lens, starts, deltas = [], [], []
    for c in chains:
        P = np.rint((np.asarray(c, np.float64) - lo) / q).astype(np.int64)
        d = np.diff(P, axis=0)
        if len(d) and np.abs(d).max() > 127:
            raise ValueError(f"polyline step {np.abs(d).max()} > 127 units of {q} nm")
        lens.append(len(P)); starts.append(P[0]); deltas.append(d)
    return (np.array(lens, np.uint32), np.array(starts, np.int64).reshape(-1, 3),
            np.concatenate(deltas).astype(np.int8) if deltas else np.zeros((0, 3), np.int8))


def unpack_polylines(sec: bytes, n: int, lo, q, extra=False):
    """-> list of (m,3) float arrays (and kinds for POLY_EXTRA)"""
    o = 0; kinds = None
    if extra:
        kinds = np.frombuffer(sec, "<u2", n, 0).copy(); o = len(pad4(b"\0" * (2 * n)))
    lens = np.frombuffer(sec, "<u4", n, o).astype(np.int64); o += 4 * n
    st = np.frombuffer(sec, "<i4", 3 * n, o).reshape(-1, 3).astype(np.int64); o += 12 * n
    nd = int(lens.sum() - n)
    d = np.frombuffer(sec, "i1", 3 * nd, o).reshape(-1, 3).astype(np.int64)
    out = []; k = 0
    for i in range(n):
        m = int(lens[i])
        P = np.vstack([st[i][None], st[i][None] + np.cumsum(d[k:k + m - 1], 0)]) if m > 1 else st[i][None]
        k += m - 1
        out.append(np.asarray(lo, np.float64) + P * q)
    return (out, kinds) if extra else out


# ------------------------------------------------------------------ polymer sources of a build
def _quat_mat(q):
    """(n,4) w,x,y,z -> (n,3,3)"""
    q = np.asarray(q, np.float64)
    q = q / np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q[:, 0], q[:, 1], q[:, 2], q[:, 3]
    return np.stack([np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)], -1),
                     np.stack([2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)], -1),
                     np.stack([2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)], -1)], 1)


def nucleosome_ends(beads, kind):
    """entry / exit points and tangents of the nucleosomal DNA (library frame; yeast webenc.linker_chains): the DNA beads
    (kind 1) are two antiparallel strands I, J"""
    d = np.asarray(beads, np.float64)[np.asarray(kind) == 1]; n = len(d) // 2; I, J = d[:n], d[n:]
    e_in, e_out = 0.5 * (I[0] + J[-1]), 0.5 * (I[-1] + J[0])
    t_in = -(0.5 * (I[3] + J[-4]) - e_in); t_out = e_out - 0.5 * (I[-4] + J[3])
    return e_in, e_out, t_in / np.linalg.norm(t_in), t_out / np.linalg.norm(t_out)


class ThreadSource:
    """polylines given as points + offsets (RNA / DNA threads of the pre-placed stage)"""

    def __init__(self, name, section, pts, off, extra_kind=0, keep=None):
        self.name, self.section, self.extra_kind = name, int(section), int(extra_kind)
        self.pts, self.off, self.keep = pts, np.asarray(off, np.int64), keep
        self.step = SECTIONS[self.section]["step"]
        self.qmax = 120 * SECTIONS[self.section]["poly_q"]          # longest step that int8 deltas hold (with margin)

    def chunks(self, max_pts=4_000_000):
        """yield (points (N,3) float64, piece lengths), deterministic, thread by thread (resampled where needed)"""
        off = self.off; nt = len(off) - 1; i = 0
        while i < nt:
            j = i; tot = 0
            while j < nt and (tot == 0 or tot + off[j + 1] - off[j] <= max_pts):
                tot += off[j + 1] - off[j]; j += 1
            P = np.asarray(self.pts[off[i]:off[j]], np.float64)
            keep = None if self.keep is None else np.asarray(self.keep[off[i]:off[j]])
            out = []
            for t in range(i, j):
                a, b = off[t] - off[i], off[t + 1] - off[i]
                out += self._thread(P[a:b], None if keep is None else keep[a:b])
            if out:
                yield np.concatenate(out), np.array([len(x) for x in out], np.int64)
            i = j

    def _thread(self, P, keep):
        if len(P) < 2:
            return []
        brk = np.zeros(len(P), bool); brk[0] = True
        st = np.linalg.norm(np.diff(P, axis=0), axis=1)
        brk[1:] |= st > GAP_NM
        if keep is not None:
            brk[1:] |= ~keep[:-1]
            brk &= keep
            if not keep.any():
                return []
        idx = np.flatnonzero(keep) if keep is not None else np.arange(len(P))
        seg_start = np.flatnonzero(brk[idx])
        out = []
        for a, b in zip(seg_start, list(seg_start[1:]) + [len(idx)]):
            R = P[idx[a:b]]
            if len(R) < 2:
                continue
            s = np.linalg.norm(np.diff(R, axis=0), axis=1)
            if s.max() > self.qmax:                      # resample along the arc length (jenc_full mrna_slab_pieces)
                c = np.r_[0.0, np.cumsum(s)]
                t = np.arange(0.0, c[-1], self.step)
                R = np.stack([np.interp(t, c, R[:, k]) for k in range(3)], 1)
            for p0 in range(0, len(R), PIECE):
                q = R[p0:p0 + PIECE]
                if len(q) >= 2:
                    out.append(q)
        return out


class LinkerSource:
    """linker DNA between consecutive nucleosomes of a chain (cubic Hermite from the exit of one to the entry of the next;
    schematic geometry without clash test, as yeast webenc.linker_chains), one piece per linker"""

    def __init__(self, name, pos, quat, chain, brk, ends, keep=None):
        self.name, self.section, self.extra_kind = name, 5, 0
        self.pos, self.quat, self.chain, self.brk, self.ends, self.keep = pos, quat, chain, brk, ends, keep

    def chunks(self, max_pairs=200_000):
        e_in, e_out, t_in, t_out = self.ends
        n = len(self.pos)
        for i0 in range(0, max(n - 1, 0), max_pairs):
            i1 = min(n, i0 + max_pairs + 1)
            p = np.asarray(self.pos[i0:i1], np.float64); R = _quat_mat(self.quat[i0:i1])
            ch = np.asarray(self.chain[i0:i1]); bk = np.asarray(self.brk[i0:i1]).astype(bool)
            ok = (ch[1:] == ch[:-1]) & ~bk[1:]
            if self.keep is not None:
                kp = np.asarray(self.keep[i0:i1]); ok &= kp[1:] & kp[:-1]
            a = np.einsum("nij,j->ni", R[:-1], e_out) + p[:-1]; ta = np.einsum("nij,j->ni", R[:-1], t_out)
            b = np.einsum("nij,j->ni", R[1:], e_in) + p[1:]; tb = np.einsum("nij,j->ni", R[1:], t_in)
            L = np.linalg.norm(b - a, axis=1)
            ok &= L <= LINKER_MAX
            k = np.flatnonzero(ok)
            m = np.maximum(2, np.ceil(L[k] * 1.15 / LINKER_STEP).astype(np.int64))
            pid = np.repeat(np.arange(len(k)), m)
            first = np.r_[0, np.cumsum(m)[:-1]]
            s = ((np.arange(len(pid)) - first[pid]) / (m[pid] - 1))[:, None]
            h00 = 2 * s ** 3 - 3 * s ** 2 + 1; h10 = s ** 3 - 2 * s ** 2 + s; h01 = -2 * s ** 3 + 3 * s ** 2; h11 = s ** 3 - s ** 2
            kk = k[pid]; Lk = L[kk][:, None]
            P = h00 * a[kk] + h10 * (ta[kk] * Lk) + h01 * b[kk] + h11 * (tb[kk] * Lk)
            if len(m):
                yield P, m


def build_sources(cfg, bdir, lib_dir, log=print):
    """polymer sources of a build (pre-placed stage outputs), as [ThreadSource | LinkerSource] in a fixed order.

    preplaced/serial/preplaced.npz (both cells): polysome mRNA (poly_mrna), mRNA threads (mrna_pts: free mRNA + polysome
    5'/3' parts), FG repeats, pre-rRNA, nascent Pol II RNA, rDNA, DNA through Pol II, mtDNA, nucleosome chains (nuc_*,
    rnuc_*). preplaced_manifest.json: further RNA bead sets packed as obstacles (Jurkat mrna_nuclear, mrna_nascent) and
    instance sets with nucleosome chains (Jurkat chromatin, with the drop lists of later pre-pass modules).
    Config export.polylines: {enabled, linker_dna (default true), skip: [source names]}."""
    pc = ((cfg.get("export", {}) or {}).get("polylines", {}) or {})
    skip = set(pc.get("skip", []) or [])
    srcs = []
    f = os.path.join(bdir, "preplaced", "serial", "preplaced.npz")
    Z = np.load(f, mmap_mode="r") if os.path.exists(f) else None
    keys = set(Z.files) if Z is not None else set()
    used_names = set()
    table = [("polysome_mrna", "poly_mrna", "poly_mrna_off", 4, 0), ("mrna", "mrna_pts", "mrna_off", 4, 0),
             ("fg_repeats", "fg_pts", "fg_off", 6, 5), ("prerrna", "prerrna_pts", "prerrna_off", 6, 4),
             ("pol2rna", "pol2rna_pts", "pol2rna_off", 6, 0), ("rdna", "rdna_pts", "rdna_off", 5, 0),
             ("pol2dna", "pol2dna_pts", "pol2dna_off", 5, 0), ("mtdna", "mtdna_pts", "mtdna_off", 5, 0)]
    for name, kp, ko, sec, xk in table:
        if kp in keys and ko in keys and name not in skip:
            srcs.append(ThreadSource(name, sec, Z[kp], Z[ko], xk)); used_names.add(name)
    man_f = os.path.join(bdir, "preplaced", "preplaced_manifest.json")
    man = json.load(open(man_f)) if os.path.exists(man_f) else {}
    for b in man.get("beads", []) or []:
        nm = b["name"]
        if nm in used_names or nm in skip or nm in ("mrna", "mtdna"):          # already read from preplaced.npz
            continue
        z = np.load(b["npz"], mmap_mode="r")
        if "pts" not in z.files or "off" not in z.files:
            continue
        keep = None
        if b.get("drop_npy"):
            keep = np.ones(len(z["pts"]), bool); keep[np.load(b["drop_npy"])] = False
        sec, xk = (6, 0) if "nascent" in nm else ((5, 0) if "dna" in nm else (4, 0))
        srcs.append(ThreadSource(nm, sec, z["pts"], z["off"], xk, keep)); used_names.add(nm)
    if pc.get("linker_dna", True) and "linker_dna" not in skip:
        ends = None
        def _ends(lib_name):
            p = os.path.join(lib_dir, f"{lib_name}.npz")
            z = np.load(p)
            return nucleosome_ends(z["beads"], z["kind"])
        lib_ids = None
        lp = os.path.join(bdir, "preplaced", "serial", "lib_ids.npy")
        if os.path.exists(lp):
            lib_ids = [str(x) for x in np.load(lp)]
        for g in ("nuc", "rnuc"):
            if all(f"{g}_{k}" in keys for k in ("pos", "q", "chain", "brk")):
                li = int(Z["li_nuc"]) if "li_nuc" in keys else None
                lname = lib_ids[li] if (lib_ids is not None and li is not None and li < len(lib_ids)) else "nucleosome"
                ends = _ends(lname)
                srcs.append(LinkerSource(f"linker_{g}", Z[f"{g}_pos"], Z[f"{g}_q"], Z[f"{g}_chain"], Z[f"{g}_brk"], ends))
        drops = {}
        for d in man.get("drops", []) or []:
            drops.setdefault(d["target"], []).append(d)
        for inst in man.get("instances", []) or []:
            z = np.load(inst["npz"], mmap_mode="r")
            if not {"lib_ids", "types", "pos", "quat", "chain", "brk"} <= set(z.files):
                continue
            ids = [str(x) for x in z["lib_ids"]]
            nucl = [i for i, x in enumerate(ids) if x.upper().startswith("NUCLEOSOME")]
            if not nucl:
                continue
            t = np.asarray(z["types"])
            keep_all = np.ones(len(t), bool)
            for d in drops.get(inst["name"], []):
                if int(d["n"]) == len(t):
                    keep_all[np.load(d["npy"])] = False
            sel = np.flatnonzero(np.isin(t, nucl))
            if "nuc_index" in z.files:                   # chain order (types interleave nucleosomes and other chromatin)
                ni = np.asarray(z["nuc_index"])[sel]; ch = np.asarray(z["chain"])[sel]
                sel = sel[np.lexsort((ni, ch))]
            srcs.append(LinkerSource(f"linker_{inst['name']}", np.asarray(z["pos"])[sel], np.asarray(z["quat"])[sel],
                                     np.asarray(z["chain"])[sel], np.asarray(z["brk"])[sel], _ends(ids[nucl[0]]),
                                     keep_all[sel]))
            log(f"   linker DNA from {inst['name']}: {len(sel):,} nucleosomes ({int((~keep_all[sel]).sum()):,} dropped)")
    log("polymer sources: " + ", ".join(f"{s.name} (sec {s.section})" for s in srcs))
    return srcs


class PolyStore:
    """all polyline pieces of a build, quantised per viewer tile and sorted by tile, in memory-mapped work files:
    per section: key (n,) int64 tile key, lens (n,) u32, starts (n,3) i32, kinds (n,) u16, delta offsets (n+1,) int64,
    deltas (D,3) int8. Two passes over the sources (count, then fill), so memory stays bounded."""

    def __init__(self, work, origin, tile_nm, kmin, kdim):
        self.work, self.origin, self.tile_nm = work, np.asarray(origin, np.float64), float(tile_nm)
        self.kmin, self.kdim = np.asarray(kmin, np.int64), np.asarray(kdim, np.int64)
        self.stats = {}

    def key_of(self, ijk):
        r = ijk - self.kmin
        inside = np.all((r >= 0) & (r < self.kdim), axis=1)
        k = (r[:, 0] * self.kdim[1] + r[:, 1]) * self.kdim[2] + r[:, 2]
        return np.where(inside, k, -1)

    def _pieces(self, P, lens):
        first = np.r_[0, np.cumsum(lens)[:-1]].astype(np.int64)
        mid = P[first + lens // 2]
        ijk = np.floor((mid - self.origin) / self.tile_nm).astype(np.int64)
        return ijk, self.key_of(ijk), first

    def build(self, sources, log=print):
        os.makedirs(self.work, exist_ok=True)
        meta = {}
        for sec in (4, 5, 6):
            S = [s for s in sources if s.section == sec]
            keys, lens, srcn = [], [], []
            for si, s in enumerate(S):                                  # pass 1: tile key and length of every piece
                n0 = sum(len(k) for k in keys)
                for P, ln in s.chunks():
                    _, k, _ = self._pieces(P, ln)
                    keys.append(k); lens.append(ln)
                n1 = sum(len(k) for k in keys)
                srcn.append((s.name, n1 - n0))
            key = np.concatenate(keys) if keys else np.zeros(0, np.int64)
            ln = np.concatenate(lens) if lens else np.zeros(0, np.int64)
            inside = key >= 0
            order = np.argsort(np.where(inside, key, np.iinfo(np.int64).max), kind="stable")
            n_in = int(inside.sum())
            dest = np.full(len(key), -1, np.int64); dest[order[:n_in]] = np.arange(n_in)
            ln_s = ln[order[:n_in]]
            doff = np.r_[0, np.cumsum(ln_s - 1)].astype(np.int64)
            q = SECTIONS[sec]["poly_q"]
            pre = os.path.join(self.work, f"poly{sec}_")
            D = np.lib.format.open_memmap(pre + "delta.npy", "w+", np.int8, (max(int(doff[-1]), 1), 3))
            ST = np.zeros((n_in, 3), np.int32); KD = np.zeros(n_in, np.uint16)
            gi = 0; maxd = 0
            for si, s in enumerate(S):                                  # pass 2: quantise into the sorted slots
                for P, ln_c in s.chunks():
                    ijk, k, first = self._pieces(P, ln_c)
                    d_ = dest[gi:gi + len(ln_c)]; gi += len(ln_c)
                    pid = np.repeat(np.arange(len(ln_c)), ln_c)
                    Q_ = np.rint((P - (self.origin + ijk[pid] * self.tile_nm)) / q).astype(np.int64)
                    within = np.arange(len(P)) - first[pid]
                    nf = (within > 0) & (d_[pid] >= 0)
                    dd = Q_[1:] - Q_[:-1]
                    dd = dd[nf[1:]]
                    if len(dd):
                        mx = int(np.abs(dd).max())
                        if mx > 127:
                            raise ValueError(f"{s.name}: polyline step {mx} > 127 units of {q} nm")
                        maxd = max(maxd, mx)
                    rows = doff[d_[pid[nf]]] + within[nf] - 1
                    D[rows] = dd.astype(np.int8)
                    ok = d_ >= 0
                    ST[d_[ok]] = Q_[first[ok]]; KD[d_[ok]] = s.extra_kind
            D.flush(); del D
            np.save(pre + "key.npy", key[order[:n_in]]); np.save(pre + "lens.npy", ln_s.astype(np.uint32))
            np.save(pre + "starts.npy", ST); np.save(pre + "kinds.npy", KD); np.save(pre + "doff.npy", doff)
            meta[sec] = dict(pieces=n_in, points=int(ln_s.sum()), outside=int((~inside).sum()), max_delta=maxd,
                             sources=[dict(name=a, pieces=int(b)) for a, b in srcn])
            log(f"POLY {sec} ({SECTIONS[sec]['name']}): {n_in:,} pieces, {int(ln_s.sum()):,} points, "
                f"{int((~inside).sum()):,} pieces outside the tile grid; " + ", ".join(f"{a} {b:,}" for a, b in srcn))
        self.stats = meta
        json.dump(meta, open(os.path.join(self.work, "poly_meta.json"), "w"), indent=1)
        return meta

    def open(self):
        self.sec = {}
        for sec in (4, 5, 6):
            pre = os.path.join(self.work, f"poly{sec}_")
            if not os.path.exists(pre + "key.npy"):
                continue
            self.sec[sec] = dict(key=np.load(pre + "key.npy", mmap_mode="r"), lens=np.load(pre + "lens.npy", mmap_mode="r"),
                                 starts=np.load(pre + "starts.npy", mmap_mode="r"), kinds=np.load(pre + "kinds.npy", mmap_mode="r"),
                                 doff=np.load(pre + "doff.npy", mmap_mode="r"), delta=np.load(pre + "delta.npy", mmap_mode="r"))
        return self

    def tile_sections(self, ijk):
        """-> {section: (bytes, n polylines)}, n points"""
        k = int(self.key_of(np.asarray(ijk, np.int64)[None])[0])
        out = {}; npts = 0
        if k < 0:
            return out, 0
        for sec, A in self.sec.items():
            a = int(np.searchsorted(A["key"], k, "left")); b = int(np.searchsorted(A["key"], k, "right"))
            if b <= a:
                continue
            lens = np.asarray(A["lens"][a:b]); d = np.asarray(A["delta"][int(A["doff"][a]):int(A["doff"][b])])
            out[sec] = (pack_polylines(lens, np.asarray(A["starts"][a:b]), d, np.asarray(A["kinds"][a:b]) if sec == 6 else None),
                        b - a)
            npts += int(lens.sum())
        return out, npts

    def tile_keys(self):
        """set of tile ijk (tuples) holding at least one piece"""
        out = set()
        for A in self.sec.values():
            u = np.unique(np.asarray(A["key"]))
            r = np.stack([u // (self.kdim[1] * self.kdim[2]), (u // self.kdim[2]) % self.kdim[1], u % self.kdim[2]], 1) + self.kmin
            out |= {tuple(int(x) for x in t) for t in r}
        return out
