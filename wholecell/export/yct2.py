"""YCT2 stream-tile encoder/decoder (docs/reference/tile_format_v2.md sections 3-5).

Scope: LOD 0 instances (packed + pre-pass, flag bit 2 for pre-pass), PALETTE + INST sections, sharded block files
with index.json byte ranges, deterministic gzip (mtime=0); since 2026-11-24 also GHOST, POLY (mRNA, DNA, extra) and MESH
sections (payloads built in meshpoly.py; PALETTE + INST bytes unchanged: ghost-only palette entries are appended after
the tile's own runs). Not yet written (declared in index.json not_written): COLUMN sections and the LOD-1 overview.
Encoder references: jurkat/viewer/jenc.py, jenc_full.py (human YCT1, sharded) and the yeast webenc.py (read-only; not
copied: the word layouts are re-implemented from the spec).
"""
from __future__ import annotations

import gzip
import hashlib
import struct
import zlib

import numpy as np

from . import quat as Q

MAGIC = struct.unpack("<I", b"YCT2")[0]
VERSION = (2 << 16) | 0
HDR_WORDS = 32
K_PALETTE, K_INST, K_GHOST, K_POLY_MRNA, K_POLY_LINK, K_POLY_EXTRA, K_MESH, K_COLUMN = 1, 2, 3, 4, 5, 6, 7, 8
FLAG_NEAR_MEMBRANE, FLAG_FIXED, FLAG_PREPASS, FLAG_GHOST = 1, 2, 4, 8
CORE_FLAG_BITS = [dict(bit=0, name="near_membrane"), dict(bit=1, name="fixed"), dict(bit=2, name="prepass"),
                  dict(bit=3, name="ghost")]
TILE_COLUMNS = ["ix", "iy", "iz", "shard", "off", "gz_bytes", "n_inst", "n_ghost", "n_poly_beads", "n_vert", "n_tri",
                "label_mask", "crc32"]


def tag4(s: str | bytes) -> int:
    b = s.encode() if isinstance(s, str) else s
    return struct.unpack("<I", hashlib.sha256(b).digest()[:4])[0]


def gz(blob: bytes) -> bytes:
    return gzip.compress(blob, compresslevel=6, mtime=0)


def encode_tile(ijk, lo, species, region, flags, pos, quat, tile_nm, q_org, q_range, lib_tag, build_tag, ghost=None,
                polys=None, mesh=None):
    """one tile -> (raw bytes, n_inst, n_palette). Instances are sorted by (species, region) (stable).
    ghost = (species, region, flags, pos, quat) of neighbour instances (written as GHOST, palette entries appended);
    polys = {section kind: (bytes, n polylines)} (meshpoly.pack_polylines); mesh = (bytes, big, n_tri)
    (meshpoly.encode_mesh)."""
    species = np.asarray(species, np.uint32); region = np.asarray(region, np.uint16)
    flags = np.asarray(flags, np.uint32); pos = np.asarray(pos, np.float64); quat = np.asarray(quat, np.float64)
    o = np.lexsort((region, species))
    species, region, flags, pos, quat = species[o], region[o], flags[o], pos[o], quat[o]
    key = species.astype(np.uint64) << np.uint64(16) | region.astype(np.uint64)
    chg = np.r_[True, key[1:] != key[:-1]] if len(key) else np.zeros(0, bool)
    pal_idx = np.cumsum(chg) - 1
    pal_s, pal_r = species[chg], region[chg]
    inst = _inst_words(pal_idx, flags, pos, quat, lo, q_org, q_range)
    secs = [(K_PALETTE, None, 0), (K_INST, inst, len(species))]
    n_ghost = 0
    if ghost is not None and len(ghost[0]):
        gs, gr, gf, gp, gq = (np.asarray(x) for x in ghost)
        go = np.lexsort((gr, gs))
        gs, gr, gf, gp, gq = gs[go].astype(np.uint32), gr[go].astype(np.uint16), gf[go].astype(np.uint32), gp[go], gq[go]
        pk = {int(k): i for i, k in enumerate(key[chg].tolist())}
        gkey = gs.astype(np.uint64) << np.uint64(16) | gr.astype(np.uint64)
        new = [int(k) for k in np.unique(gkey) if int(k) not in pk]
        for k in new:
            pk[k] = len(pk)
        if new:
            nk = np.array(new, np.uint64)
            pal_s = np.r_[pal_s, (nk >> np.uint64(16)).astype(np.uint32)]
            pal_r = np.r_[pal_r, (nk & np.uint64(0xFFFF)).astype(np.uint16)]
        gidx = np.array([pk[int(k)] for k in gkey.tolist()], np.int64)
        secs.append((K_GHOST, _inst_words(gidx, gf, gp, gq, lo, q_org, q_range), len(gs)))
        n_ghost = len(gs)
    if len(pal_s) > 65535:
        raise ValueError("tile palette > 65535 entries: split the tile")
    pal = np.stack([pal_s.astype(np.uint32), pal_r.astype(np.uint32)], 1).astype("<u4").tobytes()
    secs[0] = (K_PALETTE, pal, len(pal_s))
    for k in (K_POLY_MRNA, K_POLY_LINK, K_POLY_EXTRA):
        b = (polys or {}).get(k)
        if b and b[0]:
            secs.append((k, b[0], int(b[1])))
    big = False
    if mesh is not None and mesh[0]:
        secs.append((K_MESH, mesh[0], mesh[2])); big = bool(mesh[1])
    off = HDR_WORDS * 4 + len(secs) * 16
    table = []; body = b""
    for kind, data, cnt in secs:
        table += [kind, off + len(body), len(data), cnt]
        body += data
    hdr = np.zeros(HDR_WORDS, "<u4")
    hdr[0] = MAGIC; hdr[1] = VERSION
    hdr[2:5] = np.array(ijk, np.int32).view(np.uint32)
    hdr[5] = 0
    hdr[6:9] = np.array([tile_nm, q_org, q_range], np.float32).view(np.uint32)
    hdr[9] = len(species); hdr[10] = n_ghost; hdr[11] = len(pal_s); hdr[12] = len(secs); hdr[13] = 1 if big else 0
    hdr[14] = lib_tag; hdr[15] = build_tag
    raw = hdr.tobytes() + np.array(table, "<u4").tobytes() + body
    return raw, len(species), len(pal_s)


def _inst_words(pal_idx, flags, pos, quat, lo, q_org, q_range):
    u = np.clip(np.rint((np.asarray(pos, np.float64) - (np.asarray(lo) - q_org)) / q_range * 65535.0), 0, 65535).astype(np.uint32)
    qidx, q3 = Q.encode_smallest3(np.asarray(quat, np.float64))
    q3u = q3.view(np.uint16).astype(np.uint32)
    w0 = np.asarray(pal_idx).astype(np.uint32) | ((np.asarray(flags, np.uint32) & 0x3FFF) << 16) | (qidx << 30)
    w1 = u[:, 0] | (u[:, 1] << 16)
    w2 = u[:, 2] | (q3u[:, 0] << 16)
    w3 = q3u[:, 1] | (q3u[:, 2] << 16)
    return np.stack([w0, w1, w2, w3], 1).astype("<u4").tobytes()


def _decode_inst(words: bytes, pal, q_org, q_range):
    w = np.frombuffer(words, "<u4").reshape(-1, 4)
    pidx = w[:, 0] & 0xFFFF; flags = (w[:, 0] >> 16) & 0x3FFF; qidx = w[:, 0] >> 30
    u = np.stack([w[:, 1] & 0xFFFF, w[:, 1] >> 16, w[:, 2] & 0xFFFF], 1).astype(np.float64)
    q3 = np.stack([w[:, 2] >> 16, w[:, 3] & 0xFFFF, w[:, 3] >> 16], 1).astype(np.uint16).view(np.int16)
    return dict(species=pal[pidx, 0], region=(pal[pidx, 1] & 0xFFFF).astype(np.uint16), flags=flags,
                pos_rel=u * (q_range / 65535.0) - q_org, quat=Q.decode_smallest3(qidx, q3))


def decode_tile(raw: bytes):
    """-> dict(ijk, species, region, flags, pos_rel, quat, header, sections {kind: bytes}, counts {kind: count},
    ghost (same fields as the instances, or None)); positions relative to the tile origin lo_rel (grid frame: add the
    grid origin). MESH / POLY payloads: meshpoly.decode_mesh / unpack_polylines on sections[kind]."""
    h = np.frombuffer(raw[:HDR_WORDS * 4], "<u4")
    if h[0] != MAGIC:
        raise ValueError("not a YCT2 tile")
    ijk = h[2:5].view(np.int32).tolist()
    tile_nm, q_org, q_range = h[6:9].view(np.float32).astype(np.float64).tolist()
    ns = int(h[12])
    tab = np.frombuffer(raw[HDR_WORDS * 4: HDR_WORDS * 4 + ns * 16], "<u4").reshape(ns, 4)
    sec = {int(k): raw[o:o + l] for k, o, l, c in tab}
    cnt = {int(k): int(c) for k, o, l, c in tab}
    pal = np.frombuffer(sec[K_PALETTE], "<u4").reshape(-1, 2)
    D = _decode_inst(sec[K_INST], pal, q_org, q_range)
    lo = np.array(ijk, np.float64) * tile_nm   # relative to the grid origin; caller adds origin
    G = _decode_inst(sec[K_GHOST], pal, q_org, q_range) if K_GHOST in sec else None
    return dict(ijk=ijk, lo_rel=lo, header=h, tile_nm=tile_nm, q_org=q_org, q_range=q_range, sections=sec, counts=cnt,
                ghost=G, mesh_u32=bool(h[13] & 1), **D)


def crc32(raw: bytes) -> int:
    return zlib.crc32(raw) & 0xFFFFFFFF
