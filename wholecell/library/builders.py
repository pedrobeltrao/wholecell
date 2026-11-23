"""Structure-library builders (library mode `build`): one residue bead per residue (CA of amino acids, P / C4' of
nucleotides, nm) plus a sphere tree for overlap tests, written as <library>/<id>.npz with arrays
beads (n,3 float32, centred), kind (n int8: 0 protein, 1 nucleic acid), st_c (k,3 float32), st_r (k float32) and
optional extras (membrane label, plddt_frac_kept, pdb).

Generic parts of the yeast library build (yeast pipeline/structlib.py, build_structs.py, membrane_orient.py) ported
with their numerics unchanged; organism data (which complexes, which PDB entries, how a membrane complex is oriented)
come from config / TSV tables (plugins/<cell>/data/library/). Every function takes explicit paths (no module-level
configuration, no directories created at import).
"""
from __future__ import annotations

import gzip
import os
import re
import tarfile

import numpy as np

BEAD_R = 0.42  # nm, residue bead radius used for rendering / sphere-tree padding
# origin: yeast pipeline/membrane_orient.py:P_SIDE,N_SIDE @ e3b7448
P_SIDE = re.compile(r"lumen|extracellular|vacuol|intermembrane|perinuclear|periplasm|peroxisom", re.I)
N_SIDE = re.compile(r"cytoplasm|mitochondrial matrix|nuclear|nucleoplasm", re.I)


# ------------------------------------------------------------------ beads and sphere trees
# origin: yeast pipeline/structlib.py:beads_from_structure @ e3b7448
def beads_from_structure(st, chains=None, plddt_min=None):
    import gemmi
    st.setup_entities()
    model = st[0]
    xyz, kind, chain_ids = [], [], []
    for ch in model:
        if chains is not None and ch.name not in chains:
            continue
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            atom = None
            if info is not None and info.is_amino_acid():
                atom = res.find_atom("CA", "*"); k = 0
            elif info is not None and info.is_nucleic_acid():
                atom = res.find_atom("P", "*") or res.find_atom("C4'", "*"); k = 1
            if atom is None:
                continue
            if plddt_min is not None and atom.b_iso < plddt_min:
                continue
            xyz.append([atom.pos.x, atom.pos.y, atom.pos.z]); kind.append(k); chain_ids.append(ch.name)
    return np.asarray(xyz, float) / 10.0, np.asarray(kind, np.int8), np.asarray(chain_ids)


# origin: yeast pipeline/structlib.py:sphere_tree @ e3b7448
def sphere_tree(xyz, per=30, pad=BEAD_R, seed=0):
    """k = ceil(n / per) spheres: KMeans clusters of the beads, each sphere at the cluster mean with the 95th percentile
    bead distance + pad as radius (one sphere: all beads)"""
    from sklearn.cluster import KMeans
    n = len(xyz)
    if n == 0:
        return np.zeros((0, 3)), np.zeros(0)
    k = max(1, int(np.ceil(n / per)))
    if k == 1:
        c = xyz.mean(0, keepdims=True)
        r = np.array([np.percentile(np.linalg.norm(xyz - c, axis=1), 95) + pad])
        return c, r
    km = KMeans(n_clusters=k, n_init=1, random_state=seed, max_iter=50).fit(xyz)
    cs, rs = [], []
    for i in range(k):
        p = xyz[km.labels_ == i]
        if len(p) == 0:
            continue
        c = p.mean(0)
        cs.append(c); rs.append(np.percentile(np.linalg.norm(p - c, axis=1), 95) + pad)
    return np.array(cs), np.array(rs)


# origin: yeast pipeline/build_structs.py:rot_to_z @ e3b7448
def rot_to_z(v):
    v = v / np.linalg.norm(v)
    z = np.array([0, 0, 1.0])
    c = np.dot(v, z)
    if c > 0.999999:
        return np.eye(3)
    if c < -0.999999:
        return np.diag([1, -1, -1.0])
    ax = np.cross(v, z); s = np.linalg.norm(ax); ax /= s
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


# origin: yeast pipeline/build_structs.py:parse_ranges @ e3b7448
def parse_ranges(s, key):
    """UniProt feature string -> [(start, end, note)] for features `key` (TRANSMEM, TOPO_DOM, ...)"""
    if not isinstance(s, str):
        return []
    out = []
    for m in re.finditer(key + r" (\d+)\.\.(\d+);\s*/note=\"([^\"]*)\"", s):
        out.append((int(m.group(1)), int(m.group(2)), m.group(3)))
    return out


# origin: yeast pipeline/build_structs.py:save @ e3b7448 (library folder an argument)
def save(lib, ing_id, beads, kind, extra=None, per=30):
    beads = np.asarray(beads, float)
    c, r = sphere_tree(beads, per=per)
    d = dict(beads=beads.astype(np.float32), kind=kind.astype(np.int8), st_c=c.astype(np.float32), st_r=r.astype(np.float32))
    if extra:
        d.update(extra)
    np.savez_compressed(os.path.join(lib, f"{ing_id}.npz"), **d)
    return len(beads), len(r), float(np.linalg.norm(beads, axis=1).max())


# origin: yeast pipeline/build_structs.py:residue_beads_with_numbers @ e3b7448
def residue_beads_with_numbers(st):
    """beads + (chain, seqnum, subchain) for TM mapping"""
    import gemmi
    st.setup_entities()
    out = []
    for ch in st[0]:
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None:
                continue
            if info.is_amino_acid():
                a = res.find_atom("CA", "*"); k = 0
            elif info.is_nucleic_acid():
                a = res.find_atom("P", "*") or res.find_atom("C4'", "*"); k = 1
            else:
                continue
            if a is None:
                continue
            out.append((a.pos.x / 10, a.pos.y / 10, a.pos.z / 10, k, ch.name, res.seqid.num, res.subchain))
    return out


def assembly_path(pdb_dir, pdb):
    """RCSB biological assembly 1 (<pdb_dir>/<ID>-assembly1.cif.gz). Offline: a missing file is an error that names
    the file and how to obtain it."""
    p = os.path.join(pdb_dir, f"{pdb}-assembly1.cif.gz")
    if not os.path.exists(p):
        raise FileNotFoundError(f"PDB assembly {p} missing: download https://files.rcsb.org/download/{pdb}-assembly1.cif.gz "
                                f"into the declared pdb_dir (the library stage never downloads)")
    return p


# ------------------------------------------------------------------ membrane orientation of curated complexes
# The yeast build chose the method by complex id (NPC, Pma1_hexamer, V-ATPase); here the method is a column of the
# organism's complex table (orient = ring_axis | tm_slab | long_axis_cring | none).
# origin: yeast pipeline/build_structs.py:process_complex (NPC branch) @ e3b7448
def orient_ring_axis(xyz, rows, extra):
    w, V = np.linalg.eigh(np.cov(xyz.T)); ax = V[:, 0]  # smallest variance = ring axis
    R = rot_to_z(ax); xyz = xyz @ R.T
    xyz[:, 2] -= np.median(xyz[:, 2])
    return xyz


# origin: yeast pipeline/build_structs.py:process_complex (Pma1_hexamer branch) @ e3b7448
def orient_tm_slab(xyz, rows, tm_ranges):
    seq = np.array([r[5] for r in rows])
    tm = np.zeros(len(rows), bool)
    for a, b, _ in tm_ranges:
        tm |= (seq >= a) & (seq <= b)
    p = xyz[tm]; w, V = np.linalg.eigh(np.cov(p.T)); ax = V[:, 0]  # TM layer = slab -> smallest variance
    R = rot_to_z(ax); xyz = xyz @ R.T
    mid = xyz[tm, 2].mean(); xyz[:, 2] -= mid
    if xyz[~tm, 2].mean() < 0:  # cytoplasmic domains -> +z
        xyz[:, 1:] *= -1
    return xyz


# origin: yeast pipeline/build_structs.py:process_complex (V-ATPase branch) @ e3b7448
def orient_long_axis_cring(xyz, rows):
    w, V = np.linalg.eigh(np.cov(xyz.T)); ax = V[:, 2]
    R = rot_to_z(ax); xyz = xyz @ R.T
    # membrane ring (c-ring): chains of ~160 aa (140-220) occurring in the ring; its median height = midplane;
    # fallback: the end with fewer residues
    z = xyz[:, 2]; lo, hi = np.percentile(z, 5), np.percentile(z, 95)
    n_lo = (z < lo + 6).sum(); n_hi = (z > hi - 6).sum()
    chains = np.array([r[4] for r in rows])
    import pandas as pd
    lens = pd.Series(chains).value_counts()
    cchains = set(lens[(lens > 140) & (lens < 220)].index)
    cmask = np.isin(chains, list(cchains))
    mid = np.median(z[cmask]) if cmask.sum() > 300 else (lo + 2 if n_lo < n_hi else hi - 2)
    xyz[:, 2] -= mid
    if np.median(xyz[~cmask, 2]) < 0:  # soluble head is cytoplasmic (+z)
        xyz[:, 1:] *= -1
    return xyz


# origin: yeast pipeline/build_structs.py:process_complex @ e3b7448 (orientation method and membrane from the table)
def build_complex(lib, ing_id, pdb, pdb_dir, orient="none", membrane="", tm_ranges=None, per=40):
    import gemmi
    st = gemmi.read_structure(assembly_path(pdb_dir, pdb))
    st.remove_hydrogens()
    rows = residue_beads_with_numbers(st)
    xyz = np.array([r[:3] for r in rows]); kind = np.array([r[3] for r in rows])
    xyz -= xyz.mean(0)
    extra = {}
    if orient == "ring_axis":
        xyz = orient_ring_axis(xyz, rows, extra)
    elif orient == "tm_slab":
        xyz = orient_tm_slab(xyz, rows, tm_ranges or [])
    elif orient == "long_axis_cring":
        xyz = orient_long_axis_cring(xyz, rows)
    elif orient not in ("none", "", None):
        raise ValueError(f"complex {ing_id}: unknown orient '{orient}' (ring_axis | tm_slab | long_axis_cring | none)")
    if membrane:
        extra["membrane"] = membrane
    return save(lib, ing_id, xyz, kind, extra, per=per)


# ------------------------------------------------------------------ AlphaFold DB monomers
# origin: yeast pipeline/build_structs.py:process_af @ e3b7448
def build_af(lib, ing_id, st, tm_ranges=None, topo=None, plddt_min=50, min_kept=10, per=30):
    """CA beads of residues with pLDDT >= plddt_min (TM segments always kept; < min_kept kept: keep all); TM proteins
    rotated so that the summed TM-segment direction is +z, centred on the TM residues, N side (cytoplasm, matrix,
    nucleoplasm) to +z from the topology annotation"""
    st.remove_hydrogens()
    rows = []
    for ch in st[0]:
        for res in ch:
            a = res.find_atom("CA", "*")
            if a is not None:
                rows.append((a.pos.x / 10, a.pos.y / 10, a.pos.z / 10, res.seqid.num, a.b_iso))
    arr = np.array(rows)
    xyz, seq, pl = arr[:, :3], arr[:, 3].astype(int), arr[:, 4]
    keep = pl >= plddt_min
    if tm_ranges:
        tm = np.zeros(len(seq), bool)
        for a, b, _ in tm_ranges:
            tm |= (seq >= a) & (seq <= b)
        keep |= tm
    if keep.sum() < min_kept:  # mostly disordered: keep everything (rendered as compact blob)
        keep[:] = True
    xyz, seq = xyz[keep], seq[keep]
    extra = {"plddt_frac_kept": np.float32(keep.mean())}
    if tm_ranges:
        tm = np.zeros(len(seq), bool); vecs = []
        for a, b, _ in tm_ranges:
            m = (seq >= a) & (seq <= b); tm |= m
            if m.sum() >= 5:
                p = xyz[m]; vecs.append(p[-1] - p[0])
        vecs = np.array(vecs)
        ref = vecs[0]
        n = sum(v if np.dot(v, ref) >= 0 else -v for v in vecs)
        R = rot_to_z(n); xyz = (xyz - xyz[tm].mean(0)) @ R.T
        nsd = np.zeros(len(seq), bool); psd = np.zeros(len(seq), bool)
        for a, b, note in (topo or []):
            m = (seq >= a) & (seq <= b)
            if P_SIDE.search(note):
                psd |= m
            elif N_SIDE.search(note):
                nsd |= m
        if nsd.any() or psd.any():
            zn = xyz[nsd, 2].mean() if nsd.any() else 0.0
            zp = xyz[psd, 2].mean() if psd.any() else 0.0
            if zn - zp < 0:
                xyz[:, 1:] *= -1
        extra["membrane"] = "TM"
    else:
        xyz = xyz - xyz.mean(0)
    return save(lib, ing_id, xyz, np.zeros(len(xyz), np.int8), extra, per=per)


# origin: yeast pipeline/build_structs.py:_af_worker @ e3b7448
def _af_worker(lib, job, plddt_min, min_kept, per):
    import gemmi
    acc, ing_id, name, data, tm, tp = job
    if name.endswith(".cif.gz"):
        st = gemmi.make_structure_from_block(gemmi.cif.read_string(data)[0])
    else:
        st = gemmi.read_pdb_string(data)
    try:
        return ing_id, build_af(lib, ing_id, st, tm, tp, plddt_min=plddt_min, min_kept=min_kept, per=per)
    except Exception as e:  # noqa: BLE001 - recorded in library_af.json (failed)
        return ing_id, ("ERR", str(e))


# origin: yeast pipeline/build_structs.py:process_all_af @ e3b7448 (tar read in chunks: same jobs, bounded memory)
def build_all_af(lib, want, topo, tar_path, nproc=8, plddt_min=50, min_kept=10, per=30, chunk=600, log=print):
    """want: {uniprot: ing_id} (later rows of the species table win, as the yeast dict). topo: DataFrame indexed by
    accession with columns tm, topo (UniProt TRANSMEM / TOPO_DOM strings). Returns ({ing_id: result}, missing ids)."""
    from joblib import Parallel, delayed
    res, seen, jobs = {}, set(), []

    def flush():
        if jobs:
            out = Parallel(n_jobs=nproc, batch_size=8)(delayed(_af_worker)(lib, j, plddt_min, min_kept, per) for j in jobs)
            res.update(dict(out))
            jobs.clear()
            log(f"  AF models: {len(res)} done")
    with tarfile.open(tar_path) as tf:
        for m in tf:
            name = os.path.basename(m.name)
            if not name.endswith(".cif.gz") or "-F1-" not in name:
                continue
            acc = name.split("-")[1]
            if acc in want:
                data = gzip.decompress(tf.extractfile(m).read()).decode()
                tm = parse_ranges(topo.loc[acc, "tm"], "TRANSMEM") if acc in topo.index else []
                tp = parse_ranges(topo.loc[acc, "topo"], "TOPO_DOM") if acc in topo.index else []
                jobs.append((acc, want[acc], name, data, tm, tp)); seen.add(want[acc])
                if len(jobs) >= chunk:
                    flush()
    flush()
    return res, set(want.values()) - seen


# ------------------------------------------------------------------ curated PDB entries
# origin: yeast pipeline/build_structs.py:process_homo_pdb @ e3b7448 (offline: the assembly must be in pdb_dir)
def build_pdb_protein(lib, ing_id, pdb, pdb_dir, per=30):
    """protein residues of assembly 1 of a PDB entry (homo-oligomers of the curated table), centred"""
    import gemmi
    st = gemmi.read_structure(assembly_path(pdb_dir, pdb))
    xyz, kind, _ = beads_from_structure(st)
    xyz = xyz[kind == 0]
    return save(lib, ing_id, xyz - xyz.mean(0), np.zeros(len(xyz), np.int8), {"pdb": pdb}, per=per)


# origin: yeast scripts/s0_library.py:main (trna branch) @ e3b7448
def build_pdb_rna(lib, ing_id, pdb, pdb_dir, per=15):
    """all residue beads (protein and nucleic acid) of assembly 1 without ligands and waters, centred (tRNA-Phe 1EHZ:
    reference of the tRNA species)"""
    import gemmi
    st = gemmi.read_structure(assembly_path(pdb_dir, pdb)); st.remove_ligands_and_waters()
    xyz, knd, _ = beads_from_structure(st); xyz = xyz - xyz.mean(0)
    return save(lib, ing_id, xyz, knd, per=per)


# ------------------------------------------------------------------ fallback structures (library step `fallback`)
# Species without any model from the steps before (no AlphaFold DB entry: selenoproteins, proteins longer than the AFDB
# limit, recent UniProt entries) get one by a declared rule instead of being dropped. Two tiers, both reproducible
# offline: (1) the best-covered chain of an experimental structure of the same protein in the RCSB cache (SIFTS),
# (2) a coarse placeholder sized from the sequence length, with the UniProt TRANSMEM segments as ideal helices
# across the membrane. Every fallback entry carries `src` = "fallback:..." and `fallback` = its tier.
RES_VOL_NM3 = 1.21 * 0.110       # nm3 per residue: 1.21 nm3/kDa (protein volume of adaptive_trees) x 110 Da per residue
HELIX_RISE_NM = 0.15             # alpha helix: rise per residue
HELIX_RADIUS_NM = 0.23           # CA radius of an alpha helix
HELIX_TURN_DEG = 100.0           # rotation per residue
HELIX_SPACING_NM = 1.0           # axis distance of packed TM helices
SOLUBLE_Z_NM = 2.3               # soluble parts start outside |z| >= 2.3 nm (the slab of reorient_tm)


def _ball(m, centre, seed):
    """m points uniform in a ball of the protein volume of m residues (deterministic: Generator(seed))"""
    rng = np.random.default_rng(seed)
    R = (3.0 * m * RES_VOL_NM3 / (4.0 * np.pi)) ** (1.0 / 3.0)
    out = np.zeros((0, 3))
    while len(out) < m:
        p = rng.uniform(-1.0, 1.0, size=(2 * m + 8, 3))
        out = np.concatenate([out, p[(p ** 2).sum(1) <= 1.0]])
    return out[:m] * R + np.asarray(centre, float), R


def _hex_positions(k):
    """k helix axis positions (x, y) on a hexagonal lattice, filled shell by shell around the origin, centred"""
    pts = [(0.0, 0.0)]
    ring = 1
    while len(pts) < k:
        for side in range(6):
            a0, a1 = np.deg2rad(60 * side), np.deg2rad(60 * (side + 1))
            p0 = ring * np.array([np.cos(a0), np.sin(a0)]); p1 = ring * np.array([np.cos(a1), np.sin(a1)])
            for j in range(ring):
                pts.append(tuple(p0 + (p1 - p0) * j / ring))
        ring += 1
    P = np.array(pts[:k]) * HELIX_SPACING_NM
    return P - P.mean(0)


def placeholder_beads(n_res, tm_ranges=None, topo=None):
    """coarse placeholder of an n_res-residue protein (one bead per residue, nm): soluble = a ball of the protein volume;
    with TRANSMEM segments = ideal helices across the slab (axes on a hexagonal lattice, 1 nm apart, midplane z = 0) and
    the soluble stretches as balls outside |z| >= 2.3 nm on their side, cytoplasmic side +z (library convention). The
    side of the first stretch comes from the UniProt TOPO_DOM annotation (N side: cytoplasm, matrix, nucleoplasm; P side:
    lumen, extracellular, ...); without one the N-terminus is taken as cytoplasmic. Returns (beads, membrane flag)."""
    n_res = int(n_res)
    tms = sorted((int(a), int(b)) for a, b, _ in (tm_ranges or []) if 1 <= int(a) <= int(b) <= n_res)
    if not tms:
        xyz, _ = _ball(n_res, (0.0, 0.0, 0.0), seed=0)
        return xyz - xyz.mean(0), False
    side0 = +1.0                                          # side of residue 1: +z cytoplasmic (N side)
    for a, b, note in sorted(topo or []):
        if int(a) < tms[0][0]:
            if P_SIDE.search(note):
                side0 = -1.0
            elif N_SIDE.search(note):
                side0 = +1.0
            break
    axes = _hex_positions(len(tms))
    xyz = np.zeros((n_res, 3))
    side = side0
    # stretches: [1, tm1), tm1, (tm1, tm2), tm2, ..., (tmk, n]
    pos = 1
    for k, (a, b) in enumerate(tms):
        if a > pos:                                       # soluble stretch before this helix, on the current side
            m = a - pos
            c = (axes[k][0], axes[k][1], 0.0)
            ball, R = _ball(m, (0.0, 0.0, 0.0), seed=k + 1)
            xyz[pos - 1:a - 1] = ball + np.array([c[0], c[1], side * (SOLUBLE_Z_NM + R)])
        L = b - a + 1
        j = np.arange(L)
        z = side * (HELIX_RISE_NM * ((L - 1) / 2.0 - j))  # from the current side to the other
        ang = np.deg2rad(HELIX_TURN_DEG) * j
        xyz[a - 1:b, 0] = axes[k][0] + HELIX_RADIUS_NM * np.cos(ang)
        xyz[a - 1:b, 1] = axes[k][1] + HELIX_RADIUS_NM * np.sin(ang)
        xyz[a - 1:b, 2] = z
        side = -side
        pos = b + 1
    if pos <= n_res:                                      # C-terminal stretch
        m = n_res - pos + 1
        ball, R = _ball(m, (0.0, 0.0, 0.0), seed=len(tms) + 1)
        xyz[pos - 1:] = ball + np.array([axes[-1][0], axes[-1][1], side * (SOLUBLE_Z_NM + R)])
    xyz[:, :2] -= xyz[:, :2].mean(0)
    return xyz, True


def build_placeholder(lib, ing_id, n_res, tm_ranges=None, topo=None, per=30):
    xyz, memb = placeholder_beads(n_res, tm_ranges, topo)
    extra = {"src": "fallback:placeholder", "fallback": "placeholder_seq_length", "n_res": np.int32(n_res)}
    if memb:
        extra["membrane"] = "TM"
    return save(lib, ing_id, xyz, np.zeros(len(xyz), np.int8), extra, per=per)


def pdb_chain_candidates(sifts, acc, have_pdb, length, min_cov):
    """[(pdb, chain, residues mapped)] of the cached entries whose chain maps >= min_cov x length UniProt residues of
    acc (SIFTS ranges), best first (ties: PDB id, chain)"""
    s = sifts[sifts.SP_PRIMARY == acc]
    if not len(s):
        return []
    n = (s.SP_END.astype(int) - s.SP_BEG.astype(int) + 1).groupby([s.PDB, s.CHAIN]).sum()
    out = [(p, c, int(k)) for (p, c), k in n.items() if p in have_pdb and k >= min_cov * length]
    return sorted(out, key=lambda x: (-x[2], x[0], x[1]))


def build_pdb_chain(lib, ing_id, pdb, chain, pdb_dir, min_beads, per=30):
    """CA beads of one chain of assembly 1 (an experimental structure of a protein without an AlphaFold model), centred;
    None if fewer than min_beads residues are observed"""
    import gemmi
    st = gemmi.read_structure(assembly_path(pdb_dir, pdb)); st.remove_hydrogens()
    names = [ch.name for ch in st[0] if ch.name.split("-")[0] == chain]
    if not names:
        return None
    xyz, kind, _ = beads_from_structure(st, chains=names[:1])
    xyz = xyz[kind == 0]
    if len(xyz) < min_beads:
        return None
    extra = {"src": f"fallback:pdb_chain:{pdb}:{chain}", "fallback": "pdb_chain", "pdb": pdb}
    return save(lib, ing_id, xyz - xyz.mean(0), np.zeros(len(xyz), np.int8), extra, per=per)
