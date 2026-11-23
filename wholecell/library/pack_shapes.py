"""Packing shapes (unification plan §2b P8): sphere trees made from the library beads for packing only, separate from
the library's own trees (about 30 residues per sphere, 40 for curated complexes; those stay for the viewer simulation
and contact statistics). Yeast policy (pipeline/wc3d.py PACK_PER = 15, 2026-11-13): the library trees are fat for
packing (union 1.9-3.8x the bead union), 15 residues per sphere placed 93.1% vs 87.9% of the matrix request.

Writes <library>/pack_shapes/:
  pack_trees_<per>.npz   ids (sorted), sig (md5 of the library file the tree was made from), off, C (float64), R
  pack_meta.npz          ids (in the order asked), BR (bounding radius from the pack trees), src, pdb (structure
                         source per type: npz `src` key, the curated PDB id table, npz `pdb` key, else AF), memb, nbeads
The adaptive refinement of the Jurkat library (lib_refine_st: refine until the tree volume is <= 2.5x the bead
volume) is the library step `adaptive_trees` (B3); it changes the library's own trees, not these.
max_beads (optional): species with more beads keep the library's own tree here (KMeans at 15 residues per sphere is
too costly for the NPC, 624,240 beads; the same bound as adaptive_trees' skip_large).
"""
from __future__ import annotations

import hashlib
import os

import numpy as np

from .builders import sphere_tree


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _tree(i, path, per):
    c, r = sphere_tree(np.load(path)["beads"].astype(np.float64), per=per)
    return i, np.asarray(c, np.float64), np.asarray(r, np.float64)


# origin: yeast pipeline/wc3d.py:pack_trees @ e3b7448 (signature = md5 of the library file instead of size:mtime;
# no cache across calls: one folder per build)
def pack_trees(entries, per, nproc=8, keep=()):
    """entries: {id: library npz path}. Returns {id: (centres, radii)} at `per` residues per sphere (ids in `keep`:
    the library's own tree)."""
    from joblib import Parallel, delayed
    ids = sorted(i for i in entries if i not in keep)
    out = Parallel(n_jobs=min(nproc, 32))(delayed(_tree)(i, entries[i], per) for i in ids)
    res = {i: (c, r) for i, c, r in out}
    for i in keep:
        z = np.load(entries[i])
        res[i] = (np.asarray(z["st_c"], np.float64), np.asarray(z["st_r"], np.float64))
    return res


def write_pack_trees(trees, entries, per, out_dir):
    ids = sorted(trees)
    off = np.r_[0, np.cumsum([len(trees[i][1]) for i in ids])]
    fn = os.path.join(out_dir, f"pack_trees_{per}.npz")
    np.savez(fn, ids=np.array(ids), sig=np.array([_md5(entries[i]) for i in ids]), off=off,
             C=np.concatenate([trees[i][0] for i in ids]), R=np.concatenate([trees[i][1] for i in ids]))
    return fn


# origin: yeast scripts/s4_validate.py:src_of @ e3b7448 (curated id -> PDB table passed in)
def source_of(z, lid, curated_pdb):
    """(source, pdb id) of a library type: npz `src` ('AF-multimer', 'pred', ...), a curated PDB entry, an npz `pdb`
    key ('PDB'), else 'AF' (AlphaFold DB monomer)"""
    if "src" in z.files:
        return str(z["src"]), str(z["pdb"]) if "pdb" in z.files else ""
    if lid in curated_pdb:
        return "PDB", curated_pdb[lid]
    if "pdb" in z.files:
        return "PDB", str(z["pdb"])
    return "AF", ""


# origin: yeast pipeline/wc3d.py:build_libpack @ e3b7448 (+ s4_validate.py libpack_meta.npz: ids, BR, src, pdb)
def pack_meta(ids, entries, trees, curated_pdb):
    BR, MEM, NB, SRC, PDB = [], [], [], [], []
    for i in ids:
        z = np.load(entries[i])
        c, r = trees[i]
        BR.append(float((np.linalg.norm(c, axis=1) + r).max()))
        MEM.append("membrane" in z.files); NB.append(len(z["beads"]))
        s, p = source_of(z, i, curated_pdb)
        SRC.append(s); PDB.append(p)
    return dict(ids=np.array([str(i) for i in ids]), BR=np.array(BR), src=np.array(SRC), pdb=np.array(PDB),
                memb=np.array(MEM), nbeads=np.array(NB))


def run(entries, out_dir, per=15, ids=None, curated_pdb=None, nproc=8, log=print, max_beads=None):
    """pack shapes for `ids` (default: every entry with beads), entries {id: npz path}"""
    os.makedirs(out_dir, exist_ok=True)
    if ids is None:
        ids = []
        for i in sorted(entries):
            with np.load(entries[i]) as z:
                if "beads" in z.files and len(z["beads"]):
                    ids.append(i)
    missing = [i for i in ids if i not in entries]
    if missing:
        raise KeyError(f"pack_shapes: {len(missing)} ids without a library entry (e.g. {missing[:5]})")
    sub = {i: entries[i] for i in ids}
    keep = []
    if max_beads:
        for i in ids:
            with np.load(sub[i]) as z:
                if len(z["beads"]) > int(max_beads):
                    keep.append(i)
        if keep:
            log(f"pack shapes: {len(keep)} species with > {int(max_beads)} beads keep their library tree ({keep[:5]})")
    trees = pack_trees(sub, per, nproc, keep=keep)
    fn = write_pack_trees(trees, sub, per, out_dir)
    meta = pack_meta(ids, sub, trees, curated_pdb or {})
    np.savez(os.path.join(out_dir, "pack_meta.npz"), **meta)
    n_sph = sum(len(t[1]) for t in trees.values())
    log(f"pack shapes: {len(ids)} types, {n_sph:,} spheres at {per} residues per sphere -> {fn}")
    return fn
