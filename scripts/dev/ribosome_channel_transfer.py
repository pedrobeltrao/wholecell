# new in wholecell (B4, 2026-11-17). Superposition helpers follow yeast pipeline/polysome.py:channel_path and
# pipeline/rac.py:_trimmed @ e3b7448 (Kabsch with iterative trimming of flexible parts).
"""mRNA channel and collided-disome pose of the human 80S library entry (Jurkat RIBO80S, PDB 6QZP), transferred from
the yeast 80S (4V88 library frame; channel from 7ZW0, disome from 6I7O, both computed by the yeast pipeline) by a
superposition of the rRNA phosphates of the two ribosomes, paired by sequence alignment chain by chain (25S/28S,
18S, 5.8S, 5S; identical aligned nucleotides only), trimmed iteratively (expansion segments and flexible parts).

Writes (derived data of the common polysome builder, config `preplaced.serial` polysomes params):
  plugins/jurkat/data/preplaced/RIBO80S_mrna_channel.npy     (n, 3) nm, human library frame, 5' -> 3'
  plugins/jurkat/data/preplaced/RIBO80S_disome_6I7O.npy      [qw, qx, qy, qz, tx, ty, tz, push] (trailing in leading frame)
  plugins/jurkat/data/preplaced/RIBO80S_channel_transfer.json  provenance and fit statistics

usage: python scripts/ribosome_channel_transfer.py <yeast engine lib dir> <yeast pdb dir> <human RIBO80S npz> <6QZP assembly cif>
"""
import json
import os
import sys

import gemmi
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "plugins", "jurkat", "data", "preplaced")


def beads(st):
    """per-residue beads as structlib.beads_from_structure (CA / P or C4'), with chain and residue index"""
    st.setup_entities()
    xyz, kind, key = [], [], []
    for ch in st[0]:
        for i, res in enumerate(ch):
            info = gemmi.find_tabulated_residue(res.name)
            atom = None
            if info is not None and info.is_amino_acid():
                atom = res.find_atom("CA", "*"); k = 0
            elif info is not None and info.is_nucleic_acid():
                atom = res.find_atom("P", "*") or res.find_atom("C4'", "*"); k = 1
            if atom is None:
                continue
            xyz.append([atom.pos.x, atom.pos.y, atom.pos.z]); kind.append(k); key.append((ch.name, i))
    return np.asarray(xyz, float) / 10.0, np.asarray(kind, np.int8), key


def rna_chains(st):
    out = {}
    for ch in st[0]:
        seq, P = [], []
        for res in ch:
            info = gemmi.find_tabulated_residue(res.name)
            if info is None or not info.is_nucleic_acid():
                continue
            a = res.find_atom("P", "*")
            if a is None:
                continue
            seq.append(info.one_letter_code.upper()); P.append([a.pos.x, a.pos.y, a.pos.z])
        if len(seq) > 100:
            out[ch.name] = ("".join(seq), np.asarray(P) / 10.0)
    return out


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); U, S, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, qc - R @ pc


def trimmed(P, Q, rounds=6):
    R, t = kabsch(P, Q); keep = np.ones(len(P), bool)
    for _ in range(rounds):
        res = np.linalg.norm(P @ R.T + t - Q, axis=1)
        keep = res < max(0.3, np.percentile(res, 75)); R, t = kabsch(P[keep], Q[keep])
    res = np.linalg.norm(P @ R.T + t - Q, axis=1)
    return R, t, float(np.median(res[keep])), float(keep.mean())


def aligned_pairs(s1, s2):
    al = gemmi.align_string_sequences(list(s1), list(s2), [])
    i = j = 0; pairs = []
    for c in al.cigar_str().replace("M", "M ").replace("I", "I ").replace("D", "D ").split():
        n, op = int(c[:-1]), c[-1]
        if op == "M":
            for _ in range(n):
                if s1[i] == s2[j]:
                    pairs.append((i, j))
                i += 1; j += 1
        elif op == "I":
            i += n
        else:
            j += n
    return pairs


def quat_of(R):
    from scipy.spatial.transform import Rotation
    x, y, z, w = Rotation.from_matrix(R).as_quat()
    return np.array([w, x, y, z])


def mat_of(q):
    from scipy.spatial.transform import Rotation
    return Rotation.from_quat([q[1], q[2], q[3], q[0]]).as_matrix()


def main(ylib, ypdb, hnpz, hcif):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
    from wholecell.preplaced.builders._trees import sphere_tree
    st_y = gemmi.read_structure(os.path.join(ypdb, "4V88-assembly1.cif.gz")); st_y.remove_ligands_and_waters()
    st_h = gemmi.read_structure(hcif); st_h.remove_ligands_and_waters()
    xy, _, _ = beads(st_y)
    off_y = xy.mean(0)                                        # yeast library frame = 4V88 beads - their mean
    xh, _, _ = beads(st_h)
    zh = np.load(hnpz); Lh = zh["beads"].astype(np.float64)
    assert len(xh) == len(Lh), f"6QZP beads {len(xh)} != library beads {len(Lh)}"
    Rh, th = kabsch(xh, Lh)                                   # 6QZP -> human library frame
    dev_h = float(np.abs(xh @ Rh.T + th - Lh).max())
    # rRNA pairs yeast -> human
    cy, ch_ = rna_chains(st_y), rna_chains(st_h)
    P, Q, rows = [], [], []
    for a in sorted(cy, key=lambda k: -len(cy[k][0])):
        best = max(ch_, key=lambda b: len(aligned_pairs(cy[a][0], ch_[b][0])) if 0.5 < len(ch_[b][0]) / len(cy[a][0]) < 2.0 else -1)
        pr = aligned_pairs(cy[a][0], ch_[best][0])
        if len(pr) < 50:
            continue
        I, J = np.array(pr).T
        P.append(cy[a][1][I]); Q.append(ch_[best][1][J]); rows.append(dict(yeast=a, human=best, n_y=len(cy[a][0]), n_h=len(ch_[best][0]), pairs=len(pr)))
    P, Q = np.vstack(P), np.vstack(Q)
    Rs, ts, med, kept = trimmed(P, Q)                          # 4V88 -> 6QZP (nm)
    # composite map: yeast library frame -> human library frame
    Ra = Rh @ Rs; ta = Rh @ (Rs @ off_y + ts) + th
    ch_y = np.load(os.path.join(ylib, "ribosome_80S_mrna_channel.npy")).astype(np.float64)
    ch_h = ch_y @ Ra.T + ta
    d_y = np.min(np.linalg.norm(ch_y[:, None] - (xy - off_y)[None], axis=2), axis=1)
    d_h = np.min(np.linalg.norm(ch_h[:, None] - Lh[None], axis=2), axis=1)
    # collided disome: yeast pose (pushed) -> unpushed -> human frame -> pushed to clear the human trees
    a = np.load(os.path.join(ylib, "ribosome_80S_disome_6I7O.npy"))
    M = mat_of(a[:4]); t = a[4:7]; push_y = float(a[7]); u = t / np.linalg.norm(t)
    t0 = t - push_y * u
    Mh = Ra @ M @ Ra.T; th0 = Ra @ t0 + ta - Mh @ ta
    stc, str_ = sphere_tree(Lh, per=550)
    trees = [(stc, str_), (zh["st_c"].astype(np.float64), zh["st_r"].astype(np.float64))]
    uh = th0 / np.linalg.norm(th0); push = 0.0
    TOL = 0.3
    while any((np.linalg.norm(C[:, None] - (C @ Mh.T + th0 + push * uh)[None], axis=2) - Rr[:, None] - Rr[None]).min() < -TOL
              for C, Rr in trees):
        push += 0.25
    th1 = th0 + push * uh
    exit_entry = float(np.linalg.norm(Mh @ ch_h[-1] + th1 - ch_h[0]))
    os.makedirs(OUT, exist_ok=True)
    np.save(os.path.join(OUT, "RIBO80S_mrna_channel.npy"), ch_h.astype(np.float32))
    np.save(os.path.join(OUT, "RIBO80S_disome_6I7O.npy"), np.r_[quat_of(Mh), th1, push])
    info = dict(method="rRNA phosphate superposition yeast 4V88 -> human 6QZP (sequence-aligned identical nucleotides, "
                       "trimmed Kabsch), then 6QZP -> library frame; yeast channel (7ZW0) and disome (6I7O) carried over",
                chains=rows, rrna_pairs=int(len(P)), kept_fraction=round(kept, 3), median_dev_nm=round(med, 3),
                lib_frame_max_dev_nm=round(dev_h, 4), channel_points=int(len(ch_h)),
                channel_min_dist_to_beads_nm=dict(yeast_median=round(float(np.median(d_y)), 3), human_median=round(float(np.median(d_h)), 3),
                                                  human_min=round(float(d_h.min()), 3)),
                disome=dict(push_yeast_nm=push_y, push_human_nm=push, centre_distance_nm=round(float(np.linalg.norm(th1)), 2),
                            trailing_entry_to_leading_exit_nm=round(exit_entry, 2)),
                inputs=dict(yeast_lib=ylib, yeast_pdb=ypdb, human_npz=hnpz, human_cif=hcif))
    json.dump(info, open(os.path.join(OUT, "RIBO80S_channel_transfer.json"), "w"), indent=1)
    print(json.dumps(info, indent=1))


if __name__ == "__main__":
    main(*sys.argv[1:5])
