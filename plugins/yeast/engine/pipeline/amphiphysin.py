# origin: yeast yeast_cell_dev_bundle/pipeline/amphiphysin.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Amphiphysin (Rvs161-Rvs167) on the endocytic invaginations (open item of L8).

The BAR-domain heterodimer Rvs161/Rvs167 assembles on the tubular part of the invagination once it is ~50 nm deep and
peaks at ~100 Rvs167 at scission (~100 nm deep), enough to cover the tube from its base to the coated tip (Picco et
al. 2015 eLife 4:e04535, quantitative live imaging). Here each endocytic tube deeper than RVS_START carries
RVS_PEAK x (depth - RVS_START) / (100 - RVS_START) heterodimers (linear assembly, as observed) on the tube between the
neck and the tip region (the last TIP_COAT nm belong to the coat), BAR long axis around the tube, concave face on the
membrane. Since 2026-11-15 on a lattice: rings of as many dimers as the circumference holds (the dimer's own length),
stacked along the axis at the dimer's width and each turned by half a dimer, filled in random order (a few jittered
tries per slot); before, random positions with 40 tries per dimer jammed at about 40% of what was asked. Structure: the AlphaFold heterodimer RVS161_RVS167_het. Copies come from the free
Rvs161 and Rvs167 pools.
"""
import numpy as np
import pack3d as p3
import cell3d as c3
from config import LIB
from config import seed as _seed

LIB_ID = "RVS161_RVS167_het"
AF_MODEL, ACC_A, ACC_B = "AF-0000000210771937", "P25343", "P39743"     # AFDB predicted heterodimer Rvs161-Rvs167
RVS_START, RVS_PEAK, TIP_COAT = 50.0, 100, 35.0
LATTICE_GAP, SLOT_TRIES = 0.6, 6         # nm added to the dimer's residue extent between lattice slots; jittered tries per slot
RADIAL_EXTRA = 0.0                        # nm further out than membrane contact (0.5-1.5 tested 2026-11-15: no gain)
BEAD_GAP = 0.6                           # nm: closest allowed approach of two dimers' residue beads
BEAD_R = np.full(1, 0.4)                 # nm, residue radius for the containment test (broadcast below)
H = 20.0


def model(log=print):
    """the library entry of the heterodimer. The heterodimer step of stage 1 writes it only when it accepts the pair,
    which it has not done since the two proteins lost their shared compartment (2026-11-07); the entry on disk was a
    leftover of an earlier build, and a build from an empty library stopped here (2026-11-13)."""
    import os
    if not os.path.exists(f"{LIB}/{LIB_ID}.npz"):
        import heterodimers as hd
        n, _ = hd.build(LIB_ID, AF_MODEL, ACC_A, ACC_B, False)
        log(f"amphiphysin: {LIB_ID} built from {AF_MODEL} ({n} residues)")


def place(cell, lp, n_max, obstacles=None, seed=67, log=print):
    from scipy.spatial.transform import Rotation
    seed = _seed(seed); rng = np.random.default_rng(seed)
    ids = list(lp["ids"])
    if LIB_ID not in ids:
        log("amphiphysin: no heterodimer model"); return {}, 0
    li = ids.index(LIB_ID)
    stc, str_ = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]], lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]]
    z = np.load(f"{LIB}/{LIB_ID}.npz")["beads"].astype(np.float64)
    c0 = z.mean(0); w, V = np.linalg.eigh(np.cov((z - c0).T))
    ax_long, ax_mid, ax_short = V[:, 2], V[:, 1], V[:, 0]
    # concave side: the ends of the crescent lie towards it along the short axis, relative to the middle
    s = (z - c0) @ ax_long; t = (z - c0) @ ax_short
    ends, mid = np.abs(s) > 0.35 * np.abs(s).max(), np.abs(s) < 0.15 * np.abs(s).max()
    if t[ends].mean() < t[mid].mean():
        ax_short = -ax_short
    # PORT PATCH 4 (2026-11-17): right-handed library frame (the sign of an eigenvector is arbitrary; mid = short x long
    # makes det(Rl) = +1 whatever eigh returns). The mid axis carries no placement meaning (only long = around the tube
    # and short = concave side towards the membrane do), so its sign is free.
    ax_mid = np.cross(ax_short, ax_long)
    Rl = np.c_[ax_long, ax_mid, ax_short]          # library axes -> (around, along, outward-concave)
    # lattice spacing from the residues (the packing sphere trees, one sphere per 30 residues, make the dimer 14 x 9 nm
    # instead of 10 x 5.5): neighbouring dimers are checked against each other residue by residue, everything else
    # against the sphere trees
    zl = (z - c0) @ Rl
    length = float(np.percentile(zl[:, 0], 97) - np.percentile(zl[:, 0], 3))
    width = float(np.percentile(zl[:, 1], 97) - np.percentile(zl[:, 1], 3))
    thick = float(((stc - c0) @ ax_short + str_).max() - ((stc - c0) @ ax_short - str_).min())
    sites = []
    for i in np.where(cell.PT == c3.T_SLOT)[0]:
        q = cell.PQ[i]
        if np.linalg.norm(q[3:6] - q[:3]) <= 1.0:
            sites.append(q)
    G = cell.G; glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0); out = np.zeros(c3.NMEM)
    cap = 100000 + (len(obstacles[1]) if obstacles is not None else 0)
    G1 = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        k = np.all((obstacles[0] > glo) & (obstacles[0] < ghi), axis=1)
        p3.g_insert(glo, H, *G1, obstacles[0][k], obstacles[1][k])
    M = np.empty((3, 3)); pos, rot = [], []
    want_total = 0; slots_total = 0; placed_P = []
    for q in sites:
        tip, r_tube = q[:3], q[6]
        u = tip / np.linalg.norm(tip)                               # outward radial = tube axis towards the mouth
        mouth = u * c3.ray_depth(G, u, G[0] - G[2])
        depth = float(np.linalg.norm(mouth) - np.linalg.norm(tip) + r_tube)
        n = int(round(RVS_PEAK * min(1.0, (depth - RVS_START) / (100.0 - RVS_START)))) if depth > RVS_START else 0
        want_total += n
        h0, h1 = TIP_COAT - r_tube, depth - r_tube - 10.0         # along the axis from the tip-sphere centre
        if h1 <= h0 or n == 0:
            continue
        e1 = np.cross(u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        done = 0
        # lattice slots: rings around the tube at the dimer's length, every dimer width along the axis
        r_c = r_tube + c3.HALF + 0.5 + RADIAL_EXTRA + 0.5 * thick
        per_ring = max(1, int(np.floor(2 * np.pi * r_c / (length + LATTICE_GAP))))
        hs = np.arange(h0 + 0.5 * width, h1 - 0.5 * width + 1e-6, width + LATTICE_GAP)
        slots = [(h, 2 * np.pi * (k + 0.5 * (ri % 2)) / per_ring) for ri, h in enumerate(hs) for k in range(per_ring)]
        order = rng.permutation(len(slots))
        tries = [(0.0, 0.0)] + [(rng.normal(0, 1.0), rng.normal(0, 0.15)) for _ in range(SLOT_TRIES - 1)]
        from scipy.spatial import cKDTree
        site_beads = []; tree = None
        for si in order:
            if done >= n or len(pos) >= n_max:
                break
            for dh, dph in tries:
                h = slots[si][0] + dh; ph = slots[si][1] + dph
                radial = np.cos(ph) * e1 + np.sin(ph) * e2; around = np.cross(u, radial)
                # library (around, along, concave->membrane) = world (around, -u, -radial)
                # PORT PATCH 4 (2026-11-17): was np.c_[around, u, -radial], a left-handed frame (det -1 since around =
                # u x radial), so R = W Rl^T was an improper rotation: the residue beads checked below were the mirror image
                # of the dimer, and Rotation.from_matrix (scipy 1.15 silently; 1.18 raises) turned it into an unrelated
                # proper rotation for the stored orientation. -u makes W right-handed (det +1); long axis around the tube and
                # concave face on the membrane unchanged.
                W = np.c_[around, -u, -radial]
                R = W @ Rl.T
                ctr = tip + u * h + radial * r_c
                t = ctr - R @ c0
                qq = Rotation.from_matrix(R).as_quat(); qq = np.array([qq[3], qq[0], qq[1], qq[2]])
                P = np.empty((len(stc), 3)); p3.quat_mat(qq, M); p3.transform(stc, M, t, P)
                if p3.g_overlap(glo, H, *G1[:6], P, str_, p3.TOL):          # other molecules (sphere trees)
                    continue
                Bz = (z - c0) @ R.T + ctr                                      # this dimer's residues
                if tree is not None and np.isfinite(tree.query(Bz, distance_upper_bound=BEAD_GAP)[0]).any():
                    continue                                                    # another Rvs dimer, residue by residue
                # containment on the residues: on a tube this narrow the sphere trees nearest the membrane cut into it
                if not p3.check_volume(np.ascontiguousarray(Bz), np.full(len(Bz), 0.4), ctr[0], ctr[1], ctr[2], float(np.linalg.norm(z - c0, axis=1).max()) + 0.5, 1 << c3.CYTO, *a, out):
                    continue
                placed_P.append((P, str_)); pos.append(t); rot.append(qq); done += 1
                site_beads.append(Bz); tree = cKDTree(np.concatenate(site_beads))
                break
        slots_total += len(slots)
        for P_, r_ in placed_P:                    # this site's dimers become obstacles for the next sites
            p3.g_insert(glo, H, *G1, P_, r_)
        placed_P = []
    log(f"amphiphysin: {len(pos)} of {want_total} Rvs161-Rvs167 on {sum(1 for q in sites)} endocytic tubes; lattice {slots_total} slots "
        f"(dimer {length:.1f} x {width:.1f} nm) "
        f"(from {RVS_START:.0f} nm depth, {RVS_PEAK} at 100 nm; Picco et al. 2015)")
    return dict(rvs_pos=np.array(pos, np.float32).reshape(-1, 3), rvs_q=np.array(rot, np.float32).reshape(-1, 4), li_rvs=li), len(pos)


COAT = {"SLA2": (40, False), "SLA1": (90, True), "END3": (60, True)}   # per site; True: leaves during invagination


def coat(cell, lp, pools, obstacles=None, seed=69, log=print):
    """endocytic coat proteins on the tip of each invagination (the last TIP_COAT nm; flat sites: the membrane patch
    under the actin). Numbers per site from Picco et al. 2015: Sla2 ~40 throughout, Sla1 ~90 and End3 ~60 at the start,
    disassembling during invagination (here: linearly to none at 100 nm depth, assumption). Orientations random (no
    structure of the coat); each molecule as close to the membrane as it fits. pools: gene -> free copies."""
    from scipy.spatial.transform import Rotation
    seed = _seed(seed); rng = np.random.default_rng(seed)
    ids = list(lp["ids"])
    sub = {g: (ids.index(g), lp["STC"][lp["st_off"][ids.index(g)]:lp["st_off"][ids.index(g) + 1]],
               lp["STR"][lp["st_off"][ids.index(g)]:lp["st_off"][ids.index(g) + 1]], float(lp["BR"][ids.index(g)]))
           for g in COAT if g in ids}
    G = cell.G; glo = np.full(3, -G[0] - 10.0); ghi = -glo
    a = cell.args(glo, ghi, margin=50.0); out = np.zeros(c3.NMEM)
    cap = 100000 + (len(obstacles[1]) if obstacles is not None else 0)
    G1 = p3.grid_new(glo, ghi, H, cap, cap * 14)
    if obstacles is not None and len(obstacles[1]):
        k = np.all((obstacles[0] > glo) & (obstacles[0] < ghi), axis=1)
        p3.g_insert(glo, H, *G1, obstacles[0][k], obstacles[1][k])
    M = np.empty((3, 3)); pos, rot, typ = [], [], []
    left = {g: float(pools.get(g, 0)) for g in sub}; want = {g: 0 for g in sub}
    for i in np.where(cell.PT == c3.T_SLOT)[0]:
        q = cell.PQ[i]
        if np.linalg.norm(q[3:6] - q[:3]) > 1.0:
            continue
        tip, r = q[:3], q[6]; u = tip / np.linalg.norm(tip)
        depth = float(np.linalg.norm(u * c3.ray_depth(G, u, G[0] - G[2])) - np.linalg.norm(tip) + r)
        e1 = np.cross(u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        hc = TIP_COAT - r                                           # cylinder part of the coated region
        a_hemi, a_cyl = 2 * np.pi * r * r, 2 * np.pi * r * max(hc, 0.0)
        for g, (n0, leaves) in COAT.items():
            if g not in sub:
                continue
            n = int(round(n0 * (max(0.0, 1.0 - depth / 100.0) if leaves else 1.0)))
            want[g] += n
            li, stc, str_, br = sub[g]
            done = 0
            for _ in range(n * 30):
                if done >= n or left[g] < 1:
                    break
                ph = rng.uniform(0, 2 * np.pi); radial_ip = np.cos(ph) * e1 + np.sin(ph) * e2
                if rng.random() < a_hemi / (a_hemi + a_cyl):         # on the tip hemisphere (cosine-uniform)
                    ct = rng.random(); d = -u * ct + radial_ip * np.sqrt(1 - ct * ct); base = tip
                else:
                    d = radial_ip; base = tip + u * rng.uniform(0, hc)
                qq = Rotation.random(random_state=rng.integers(1 << 31)).as_quat(); qq = np.array([qq[3], qq[0], qq[1], qq[2]])
                p3.quat_mat(qq, M)
                for off in np.arange(c3.HALF + 1.0, c3.HALF + 1.0 + 2 * br, 1.5):
                    t = base + d * (r + off)
                    P = np.empty((len(stc), 3)); p3.transform(stc, M, t, P)
                    if p3.g_overlap(glo, H, *G1[:6], P, str_, p3.TOL):
                        continue
                    if p3.check_volume(P, str_, t[0], t[1], t[2], br, 1 << c3.CYTO, *a, out):
                        p3.g_insert(glo, H, *G1, P, str_); pos.append(t); rot.append(qq); typ.append(li); done += 1; left[g] -= 1
                        break
    used = {g: int(sum(1 for t in typ if t == sub[g][0])) for g in sub}
    log("endocytic coat: " + ", ".join(f"{g} {used[g]} of {want[g]}" for g in sub) + " at the invagination tips (Picco et al. 2015)")
    return dict(coat_pos=np.array(pos, np.float32).reshape(-1, 3), coat_q=np.array(rot, np.float32).reshape(-1, 4),
                coat_type=np.array(typ, np.int32)), used
