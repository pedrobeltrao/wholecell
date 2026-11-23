# origin: yeast scripts/s3_pack.py:main ("proteins on the mRNA ends" block) @ e3b7448
# origin: yeast pipeline/mrnp.py:ends, place @ e3b7448
# Changed (unification B4, 2026-11-17): the numerics live here (no ctx.mod); geometry primitives from ctx.K; seed =
# the integer of ctx.seed; the poly(A)-binding / cap-binding library ids, footprint and caps are params with the yeast
# values as defaults; the mRNA bead radius comes from the mRNA builder (ctx.state mrna_bead_r); per-gene poly(A)
# lengths come from ctx.state ("polya" {gene: nt}, or "gene_table" with gene, polya columns) or the context's
# gene_table() data function, not from a gene module.
"""Proteins on the mRNA ends (serial pre-placed builder `core.mrnp_ends`).

Poly(A)-binding protein on the 3' tails (per-gene poly(A) length when gene-specific mRNAs are used) and eIF4E on the
5' caps, from their free cytosolic pools, avoiding the mRNA-thread obstacles. Adds pre: rnp_*; deducts the placed
copies.

params: seed (yeast 37), pab (PAB1), cap (CDC33), pool_region (cytosol; null = all rows of the pool), footprint_nt (27),
pab_max (3), touch_nm (1.5), tries (60), own_nt (40), labels ([cytosol]).
"""
import numpy as np

from wholecell.registry import register

# origin: yeast pipeline/mrnp.py constants @ e3b7448
PAB_FOOTPRINT_NT, PAB_MAX = 27, 3
TOUCH_NM, TRIES = 1.5, 60
OWN_NT = 40


def ends(pre, own_nt=OWN_NT):
    """(3' end points, 3' genes, 5' end points) of all threads"""
    # origin: yeast pipeline/mrnp.py:ends @ e3b7448
    P, off, kind = pre["mrna_pts"], pre["mrna_off"], pre["mrna_kind"]
    g = pre.get("mrna_gene")
    e3, g3, e5, own3, own5 = [], [], [], [], []
    for i, k in enumerate(kind):
        a, b = off[i], off[i + 1]
        if b - a < 2:
            continue
        if k in (0, 1):
            e3.append(b - 1); g3.append(g[i] if g is not None else ""); own3.append((max(a, b - own_nt), b))
        if k == 0:
            e5.append(a); own5.append((a, min(b, a + own_nt)))
        if k == 2:
            e5.append(b - 1); own5.append((max(a, b - own_nt), b))
    return np.array(e3, np.int64), np.array(g3), np.array(e5, np.int64), own3, own5


def place(ctx, n_pab, n_cap, obstacles, seed, polya=None, log=print, pab_id="PAB1", cap_id="CDC33", bead_r=0.5,
          footprint_nt=PAB_FOOTPRINT_NT, pab_max=PAB_MAX, touch_nm=TOUCH_NM, tries=TRIES, own_nt=OWN_NT, mask=None):
    """obstacles: (centres, radii) of everything placed so far except the mRNA beads (added here).
    polya: {gene: poly(A) nt}. Returns pre-pass group 'rnp' (rnp_type, rnp_pos, rnp_q)."""
    # origin: yeast pipeline/mrnp.py:place @ e3b7448 (seed given as an integer)
    from scipy.spatial import cKDTree
    from scipy.spatial.transform import Rotation
    K, cell, lp, pre = ctx.K, ctx.cell, ctx.lp, ctx.pre
    p3 = K.p3
    ids = list(lp["ids"])
    if "mrna_pts" not in pre or pab_id not in ids:
        return {}
    rng = np.random.default_rng(seed)
    P = pre["mrna_pts"].astype(np.float64)
    e3, g3, e5, own3, own5 = ends(pre, own_nt)
    OC, OR = obstacles
    AC = np.concatenate([OC, P]); AR = np.concatenate([OR, np.full(len(P), bead_r)])
    is_bead = np.r_[np.zeros(len(OC), bool), np.ones(len(P), bool)]
    bead_idx = np.r_[np.full(len(OC), -1), np.arange(len(P))]
    tree = cKDTree(AC); rmax = AR.max()
    G = cell.G; glo = np.full(3, -G[0] - 10.0); a = cell.args(glo, -glo, margin=50.0)
    out = np.zeros(K.nmem)
    if mask is None:
        mask = 1 << K.lab("cytosol")
    placed_c, placed_r = [], []
    st = dict(tree=None, n=0, R=None)

    def free(S, R, own):
        for nb, s, r in zip(tree.query_ball_point(S, R.max() + rmax), S, R):
            if not nb:
                continue
            nb = np.array(nb)
            nb = nb[~(is_bead[nb] & (bead_idx[nb] >= own[0]) & (bead_idx[nb] < own[1]))]
            if len(nb) and (np.linalg.norm(AC[nb] - s, axis=1) < AR[nb] + r - p3.TOL).any():
                return False
        if len(placed_c) - st["n"] > 200:
            st.update(tree=cKDTree(np.concatenate(placed_c)), n=len(placed_c), R=np.concatenate(placed_r))
        if st["tree"] is not None:
            for nb, s, r in zip(st["tree"].query_ball_point(S, R.max() + st["R"].max()), S, R):
                if nb and (np.linalg.norm(st["tree"].data[nb] - s, axis=1) < st["R"][nb] + r - p3.TOL).any():
                    return False
        for C_, R_ in zip(placed_c[st["n"]:], placed_r[st["n"]:]):
            if np.linalg.norm(C_.mean(0) - S.mean(0)) < 60 and (np.linalg.norm(S[:, None] - C_[None], axis=2) < R[:, None] + R_[None] - p3.TOL).any():
                return False
        return True

    types, pos, quat = [], [], []

    def put(li, E, own):
        C = lp["STC"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        R = lp["STR"][lp["st_off"][li]:lp["st_off"][li + 1]].astype(np.float64)
        for _ in range(tries):
            Rm = Rotation.random(random_state=rng).as_matrix()
            S0 = C @ Rm.T
            k = rng.integers(len(S0))
            u = rng.normal(size=3); u /= np.linalg.norm(u)
            t = E - S0[k] + u * (R[k] + bead_r + rng.uniform(0, touch_nm))
            S = S0 + t
            if not free(S, R, own):
                continue
            if not p3.check_volume(S, R, t[0], t[1], t[2], float(lp["BR"][li]), mask, *a, out):
                continue
            placed_c.append(S); placed_r.append(R)
            types.append(li); pos.append(t); quat.append(np.roll(Rotation.from_matrix(Rm).as_quat(), 1))
            return True
        return False

    li_pab = ids.index(pab_id); li_cap = ids.index(cap_id) if cap_id in ids else -1
    n3 = np.ones(len(e3), np.int64)
    if polya is not None:
        n3 = np.clip(np.array([polya.get(g, footprint_nt) for g in g3]) // footprint_nt, 1, pab_max).astype(np.int64)
    got_p = 0
    for i in rng.permutation(len(e3)):
        for _ in range(n3[i]):
            if got_p >= n_pab:
                break
            if put(li_pab, P[e3[i]], own3[i]):
                got_p += 1
    got_c = 0
    if li_cap >= 0:
        for i in rng.permutation(len(e5)):
            if got_c >= n_cap:
                break
            if put(li_cap, P[e5[i]], own5[i]):
                got_c += 1
    log(f"mRNP ends: {got_p} Pab1 on {len(e3)} 3' ends ({int(n3.sum())} wanted), {got_c} eIF4E (Cdc33) on {len(e5)} 5' ends")
    if not types:
        return {}
    return dict(rnp_type=np.array(types, np.int32), rnp_pos=np.array(pos, np.float32), rnp_q=np.array(quat, np.float32),
                n_pab=got_p, n_cap=got_c)


def _polya(ctx):
    """per-gene poly(A) lengths {gene: nt} when gene-specific mRNAs were threaded, else None"""
    st = ctx.state
    if st.get("tx") is None:
        return None
    if st.get("polya") is not None:
        return st["polya"]
    gt = st.get("gene_table")
    if gt is None and callable(getattr(ctx, "gene_table", None)):
        gt = ctx.gene_table()
    if gt is None:
        raise KeyError("core.mrnp_ends: gene-specific mRNAs but no poly(A) table (ctx.state polya / gene_table, or "
                       "ctx.gene_table())")
    return gt.drop_duplicates("gene").set_index("gene").polya.to_dict()


@register("preplaced", "core.mrnp_ends")
def build(ctx, params):
    exp3 = ctx.exp3
    pab_id, cap_id = params.get("pab", "PAB1"), params.get("cap", "CDC33")
    reg = params.get("pool_region", "cytosol")
    inreg = (exp3.region == reg) if reg else np.ones(len(exp3), bool)
    pab = (exp3.lib_id == pab_id) & inreg; cap = (exp3.lib_id == cap_id) & inreg
    labels = params.get("labels")
    rnp = place(ctx, float(exp3.loc[pab, "expected"].sum()), float(exp3.loc[cap, "expected"].sum()),
                ctx.state["mrna_obstacles"], ctx.seed(params.get("seed", 37)), polya=_polya(ctx), log=ctx.log,
                pab_id=pab_id, cap_id=cap_id, bead_r=float(ctx.state.get("mrna_bead_r", 0.5)),
                footprint_nt=int(params.get("footprint_nt", PAB_FOOTPRINT_NT)), pab_max=int(params.get("pab_max", PAB_MAX)),
                touch_nm=float(params.get("touch_nm", TOUCH_NM)), tries=int(params.get("tries", TRIES)),
                own_nt=int(params.get("own_nt", OWN_NT)), mask=ctx.K.mask(*labels) if labels else None)
    if rnp:
        ctx.pre.update({k: v for k, v in rnp.items() if k.startswith("rnp_")})
        exp3.loc[pab, "expected"] = (exp3.loc[pab, "expected"] - rnp["n_pab"]).clip(lower=0)
        exp3.loc[cap, "expected"] = (exp3.loc[cap, "expected"] - rnp["n_cap"]).clip(lower=0)
