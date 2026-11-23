# origin: yeast scripts/s3_pack.py:main ("transport receptors in the pore channels" block) @ e3b7448
# origin: yeast pipeline/pore_receptors.py:place @ e3b7448; pipeline/pack3d.py:LocalGrid @ e3b7448
# Changed (unification B4, 2026-11-17): the numerics live here (no ctx.mod); geometry primitives from ctx.K; seed =
# the integer of ctx.seed; receptor list, receptors per pore, channel size and tries are params with the yeast values
# as defaults; the search box is the nucleus of the yeast X5 layout (G[3:6] +- max G[11:14] + 100 nm) for the yeast
# kernels and the pores' bounding box + 300 nm otherwise (only the obstacle crop and the evaluator's prim selection
# depend on it).
"""Transport receptors in the NPC channels (serial pre-placed builder `core.npc_receptors`).

Karyopherins, NTF2, Ran, mRNA export receptors in the central channel of every pore, from their free pools, avoiding
every instance pre-placed so far. Adds pre: kap_*; deducts the placed copies.

params: seed (yeast 73), receptors (library ids; yeast list below), n_per_pore (20), ch_r (22 nm), ch_z (25 nm),
tries (300), labels (allowed: [nuclear_pore_channel, nucleoplasm, cytosol]).
"""
import numpy as np

from wholecell.registry import register

# origin: yeast pipeline/pore_receptors.py constants @ e3b7448
RECEPTORS = ["KAP95", "SRP1", "KAP123", "PSE1", "CPX_26769_RAN_RANBP1_CRM1_nuclear_export_c",
             "CPX_1142_MEX67_MTR2_mRNA_nuclear_export_f", "NTF2_oligo"]
N_PER_PORE = 20
CH_R, CH_Z = 22.0, 25.0
TRIES = 300


class LocalGrid:
    """overlap test for pre-pass placements in a box (spatial hash of the kernels' pack3d): obstacles inside the box,
    then every placed molecule"""
    # origin: yeast pipeline/pack3d.py:LocalGrid @ e3b7448 (grid functions from K.p3)

    def __init__(self, K, glo, ghi, OC, OR, n_new_spheres=200000):
        self.p3 = p3 = K.p3
        self.lo = np.asarray(glo, np.float64); hi = np.asarray(ghi, np.float64)
        k = np.all((OC > self.lo) & (OC < hi), axis=1)
        OC, OR = np.ascontiguousarray(OC[k], np.float64), np.ascontiguousarray(OR[k], np.float64)
        cap_s = len(OR) + n_new_spheres
        self.g = p3.grid_new(self.lo, hi, p3.H_GRID, cap_s, cap_s * 14)
        if len(OR):
            p3.g_insert(self.lo, p3.H_GRID, *self.g, OC, OR)

    def free(self, P, R):
        dims, head, nxt, sid, cen, rad, cnt = self.g
        return not self.p3.g_overlap(self.lo, self.p3.H_GRID, dims, head, nxt, sid, cen, rad, np.ascontiguousarray(P, np.float64),
                                     np.ascontiguousarray(R, np.float64), self.p3.TOL)

    def add(self, P, R):
        return self.p3.g_insert(self.lo, self.p3.H_GRID, *self.g, np.ascontiguousarray(P, np.float64), np.ascontiguousarray(R, np.float64))


def _box(K, cell):
    G = cell.G
    if K.name == "yeast_x5":           # nucleus centre and largest radius of the X5 layout header
        cn = G[3:6]; rn = float(G[11:14].max())
        return cn - rn - 100, cn + rn + 100
    P = np.asarray(cell.pores, np.float64)[:, :3]
    return P.min(0) - 300.0, P.max(0) + 300.0


def place(ctx, pools, obstacles, seed, log=print, receptors=RECEPTORS, n_per_pore=N_PER_PORE, ch_r=CH_R, ch_z=CH_Z,
          tries=TRIES, mask=None):
    """pools: {lib_id: free copies}. Returns pre-pass group 'kap' (kap_type, kap_pos, kap_q), per-type counts and
    the placed spheres."""
    # origin: yeast pipeline/pore_receptors.py:place @ e3b7448 (seed given as an integer)
    from scipy.spatial.transform import Rotation
    K, cell, lp = ctx.K, ctx.cell, ctx.lp
    ids = list(lp["ids"])
    rng = np.random.default_rng(seed)
    types_ok = [g for g in receptors if g in ids and pools.get(g, 0) >= 1]
    if not types_ok or not len(cell.pores):
        return {}, {}, (np.zeros((0, 3)), np.zeros(0))
    w = np.array([pools[g] for g in types_ok], float); w /= w.sum()
    glo, ghi = _box(K, cell)
    a = cell.args(glo, ghi, margin=200.0)
    OC, OR = obstacles
    grid = LocalGrid(K, glo, ghi, OC.astype(np.float64), OR.astype(np.float64))
    out = np.zeros(K.nmem)
    if mask is None:
        mask = K.mask("nuclear_pore_channel", "nucleoplasm", "cytosol")
    types, pos, quat, counts, sph_c, sph_r = [], [], [], {}, [], []
    for pr in cell.pores:
        c0, u = pr[:3].astype(float), pr[3:6].astype(float); u /= np.linalg.norm(u)
        e1 = np.cross(u, [0.3, 0.5, 0.8]); e1 /= np.linalg.norm(e1); e2 = np.cross(u, e1)
        for g in rng.choice(types_ok, n_per_pore, p=w):
            if counts.get(g, 0) >= pools[g]:
                continue
            li = ids.index(g)
            C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
            R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
            for _ in range(tries):
                r = ch_r * np.sqrt(rng.random()); ph = rng.uniform(0, 2 * np.pi); z = rng.uniform(-ch_z, ch_z)
                ctr = c0 + z * u + r * (np.cos(ph) * e1 + np.sin(ph) * e2)
                Rm = Rotation.random(random_state=rng).as_matrix()
                t = ctr - Rm @ C.mean(0); S = C @ Rm.T + t
                if grid.free(S, R) and K.p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), mask, *a, out):
                    grid.add(S, R); sph_c.append(S); sph_r.append(R)
                    types.append(li); pos.append(t); quat.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1))
                    counts[g] = counts.get(g, 0) + 1
                    break
    log(f"pore receptors: {len(types)} of {n_per_pore * len(cell.pores)} in the channels of {len(cell.pores)} pores ("
        + ", ".join(f"{g.split('_')[0] if not g.startswith('CPX') else g.split('_')[2]} {n}" for g, n in counts.items()) + ")")
    if not types:
        return {}, counts, (np.zeros((0, 3)), np.zeros(0))
    return (dict(kap_type=np.array(types, np.int32), kap_pos=np.array(pos, np.float32), kap_q=np.array(quat, np.float32)),
            counts, (np.concatenate(sph_c), np.concatenate(sph_r)))


@register("preplaced", "core.npc_receptors")
def build(ctx, params):
    exp3 = ctx.exp3
    gt0, gp0, gq0 = ctx.pre_instances()
    KC, KR = ctx.emit(gt0, gp0, gq0)
    labels = params.get("labels")
    rk, k_counts, _ = place(ctx, exp3.groupby("lib_id").expected.sum().to_dict(), (KC, KR), ctx.seed(params.get("seed", 73)),
                            log=ctx.log, receptors=list(params.get("receptors", RECEPTORS)),
                            n_per_pore=int(params.get("n_per_pore", N_PER_PORE)), ch_r=float(params.get("ch_r", CH_R)),
                            ch_z=float(params.get("ch_z", CH_Z)), tries=int(params.get("tries", TRIES)),
                            mask=ctx.K.mask(*labels) if labels else None)
    ctx.pre.update(rk)
    for g_, n_ in k_counts.items():
        ctx.consume(list(exp3.loc[exp3.lib_id == g_, "ing_id"].unique()), n_)
