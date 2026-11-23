"""Plain-Python parts of the packer switches (wholecell/pack/switches.py): size-class rounds (P1), pilot thinning (P2),
shell cap (P3), per-region budgets (P4), squeeze job selection (P5) and exclusion rules (P7). The numba parts are in
pack3d.py (budgets, exclusion test) and squeeze.py. Each function is ported from the yeast X5 packer (origin comments)
and works on explicit arrays, so that it can be tested on its own (tests/test_pack_switches.py).
"""
from __future__ import annotations

import numpy as np


def sphere_volumes(so, STR):
    """per library: summed sphere volume of its collision tree (nm^3; the spheres of one molecule overlap)"""
    return np.array([float((4.0 / 3.0 * np.pi * STR[a:b] ** 3).sum()) for a, b in zip(so[:-1], so[1:])])


# ------------------------------------------------------------------ P1 size-class rounds
# origin: yeast wc3d.py:round_edges @ e3b7448
# Changed: explicit arrays; membrane jobs are those of pack_class 0 (core region kinds).
def round_edges(lib, region, n, so, STR, BR, pack_class, n_rounds):
    """bounding-radius thresholds (descending, n_rounds - 1 values) that split the requested non-membrane copies
    (count entries lib, region, n) into n_rounds size classes of about equal summed sphere volume, largest first"""
    if n_rounds <= 1:
        return np.zeros(0)
    k = pack_class[region] != 0
    lib, n = lib[k], n[k]
    if not len(lib):
        return np.zeros(n_rounds - 1)
    br = BR[lib]
    v3 = np.array([float((STR[a:b] ** 3).sum()) for a, b in zip(so[:-1], so[1:])])
    v = v3[lib] * n
    o = np.argsort(-br, kind="stable"); cum = np.cumsum(v[o]) / v.sum()
    return np.array([br[o][min(np.searchsorted(cum, j / n_rounds), len(o) - 1)] for j in range(1, n_rounds)])


# origin: yeast wc3d.py:job_rounds @ e3b7448
def job_rounds(br, is_mem, edges):
    """packing round of each job: membrane jobs and the largest species in round 0, smaller size classes later"""
    if edges is None or len(edges) == 0:
        return np.zeros(len(br), np.int64)
    return np.where(is_mem, 0, np.searchsorted(-np.asarray(edges), -np.asarray(br), side="left")).astype(np.int64)


def phase_order(phases, rnd, alternate):
    """checkerboard phase order of round rnd (reversed in every second round when alternate)"""
    return list(phases) if (not alternate or rnd % 2 == 0) else list(phases)[::-1]


# ------------------------------------------------------------------ P2 pilot-calibrated thinning
# origin: yeast wc3d.py:pilot_capacity @ e3b7448 (region test and tile choice)
def crowded_regions(lib, region, n, vol_lib, region_volume, vol_regions, ref, min_rel):
    """{region index: requested sphere volume fraction} of the volume regions asking for more than min_rel x the
    reference region's requested fraction (vol_regions: indices of volume regions; region_volume: nm^3 per index)"""
    req = np.bincount(region, weights=vol_lib[lib] * n, minlength=len(region_volume))
    phi = {int(r): req[r] / region_volume[r] for r in vol_regions if region_volume[r] > 0}
    p_ref = phi.get(ref, 0.0)
    return {r: f for r, f in phi.items() if r != ref and req[r] > 0 and f > min_rel * p_ref}, phi


# origin: yeast wc3d.py:pilot_capacity @ e3b7448 (tile choice: in proportion to the region's volume, not adjacent)
def pilot_tiles(regions, w_by_region, Tij, n_tiles, rng):
    """{region: [tile ids]}: per region up to n_tiles tiles drawn in proportion to the region's volume in them, no two
    chosen tiles (over all regions) adjacent (Chebyshev distance > 1)"""
    chosen, out = [], {}
    for r in regions:
        w = np.asarray(w_by_region[r], np.float64); mine = []
        nz = int((w > 0).sum())
        cand = rng.choice(len(w), size=min(len(w), 40 * n_tiles), replace=False, p=w / w.sum()) if nz > 40 * n_tiles \
            else rng.permutation(np.nonzero(w)[0])
        for ti in cand:
            if len(mine) == n_tiles:
                break
            if all(np.abs(Tij[ti] - Tij[tj]).max() > 1 for tj in chosen):
                chosen.append(int(ti)); mine.append(int(ti))
        out[r] = mine
    return out


def thinning_factor(held, safety, keep_above):
    """request factor of a region whose pilot tiles held `held` of the requested sphere volume"""
    return 1.0 if held >= keep_above else safety * held


# origin: yeast wc3d.py:thin @ e3b7448
# Changed: the yeast packer re-draws the whole tile split from thinned totals; here each (tile, lib, region) count is
# thinned binomially (same distribution: a thinned multinomial is a multinomial of the binomially thinned total).
def thin_counts(tile, lib, region, n, factors, rng):
    """counts with every entry of region r kept binomially with probability factors[r]; zero entries dropped"""
    n2 = n.astype(np.int64).copy()
    for r, f in sorted(factors.items()):
        k = region == r
        n2[k] = rng.binomial(n2[k], f)
    keep = n2 > 0
    return tile[keep], lib[keep], region[keep], n2[keep]


# ------------------------------------------------------------------ P3 shell cap
# origin: yeast wc3d.py:shell_capacity @ e3b7448
# Changed: core shells (membrane area x width) or measured volumes (volume-measured shells = yeast bands); the cap is the
# requested density of the bulk volume region on the shell's side (yeast: always the cytosol); pore-limited shells
# are skipped.
def shell_cap(tot, vol_lib, br_lib, vol_regions, region_volume, shells, log=print):
    """tot: {(lib, region): expected copies} (edited copy returned). vol_regions: {region index: label codes} of label
    volume regions; region_volume: {region index: nm^3}; shells: list of dict(ri, area_nm2, d0, d1, add_radius,
    pore_limit, labels, optional volume_nm3: measured volume of a volume-measured shell, used instead of area x width
    (yeast bands)). Returns (tot, report {shell region: (factor, copies moved, target region)})."""
    tot = dict(tot)
    req = {}
    for (li, ri), x in tot.items():
        req[ri] = req.get(ri, 0.0) + x * vol_lib[li]
    phi = {ri: req.get(ri, 0.0) / v for ri, v in region_volume.items() if v > 0}
    rep = {}
    for s in shells:
        vs = s.get("volume_nm3")
        if s["pore_limit"] > 0 or (s["area_nm2"] <= 0 if vs is None else vs <= 0):
            continue
        keys = [k for k in tot if k[1] == s["ri"] and tot[k] > 0]
        if not keys:
            continue
        target = next((ri for ri, labs in vol_regions.items() if s["labels"][0] in labs), None)
        if target is None or target not in phi:
            continue
        w = s["d1"] - s["d0"]
        if s["add_radius"]:
            ww = np.array([tot[k] * vol_lib[k[0]] for k in keys])
            w += float((ww * np.array([br_lib[k[0]] for k in keys])).sum() / ww.sum())
        p_s = sum(tot[k] * vol_lib[k[0]] for k in keys) / (s["area_nm2"] * w if vs is None else vs)
        cap = phi[target]
        if p_s <= cap:
            continue
        f = cap / p_s; moved = 0.0
        for k in keys:
            x = tot[k]; tot[k] = x * f; moved += x * (1 - f)
            tk = (k[0], target); tot[tk] = tot.get(tk, 0.0) + x * (1 - f)
        rep[s["ri"]] = (f, moved, target)
        log(f"  shell cap {s['ri']}: asked for sphere volume fraction {p_s:.3f}, capped at {cap:.3f} (region {target}): "
            f"x {f:.3f} for every species, {moved:,.0f} copies moved")
    return tot, rep


# ------------------------------------------------------------------ P4 per-region budgets
def budget_arrays(per_region, rnames):
    """(REG_TRIES, REG_GIVEUP) int64 arrays over regions (0 = global); empty when no region has a budget"""
    if not per_region:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    tr = np.zeros(len(rnames), np.int64); gu = np.zeros(len(rnames), np.int64)
    for name, b in per_region.items():
        i = rnames.index(name)
        tr[i] = int(b.get("max_tries", 0)); gu[i] = int(b.get("give_up", 0))
    return tr, gu


# ------------------------------------------------------------------ P5 squeeze job selection
def squeeze_jobs(scope, threshold, jl, jr, jn, placed, mass, is_label_volume):
    """bool per job: jobs whose region (label volume regions only) is squeezed in this tile. underfilled: regions placed
    below `threshold` of the requested mass in the tile; all: every label volume job of a tile with any deficit (origin:
    yeast wc3d.squeeze_tile: all soluble molecules move; copies are inserted only where a job has a deficit)"""
    sel = np.zeros(len(jl), bool)
    if scope == "none" or not len(jl):
        return sel
    lv = is_label_volume[jr]
    if scope == "all":                                   # yeast: every label-volume molecule moves (inserts: deficits)
        return lv.copy() if (lv & (placed < jn)).any() else sel
    for r in np.unique(jr[lv]):
        k = jr == r
        need = float((jn[k] * mass[jl[k]]).sum()); got = float((placed[k] * mass[jl[k]]).sum())
        if need <= 0 or got >= need:
            continue
        if scope == "all" or got / need < threshold:
            sel |= k
    return sel


# ------------------------------------------------------------------ P7 exclusion rules
def _match(ids, libs, prefixes):
    ex = set(libs or []); pf = tuple(prefixes or [])
    return np.array([(str(i) in ex) or (bool(pf) and str(i).startswith(pf)) for i in ids], bool)


def exclusion_tables(rules, packed_libs, preplaced_libs):
    """EXC_TYPE (bitmask of rules per packed library), around (bool (n_rules, n preplaced libs)), radii (n_rules,)"""
    et = np.zeros(len(packed_libs), np.int64)
    around = np.zeros((len(rules), len(preplaced_libs)), bool); rad = np.zeros(len(rules))
    for i, e in enumerate(rules):
        et[_match(packed_libs, e.get("libs"), e.get("lib_prefixes"))] |= np.int64(1) << i
        if e["around"].get("provider"):                  # zone spheres from a provider (exclusion_provider_spheres)
            continue
        around[i] = _match(preplaced_libs, e["around"].get("libs"), e["around"].get("lib_prefixes"))
        rad[i] = float(e["r_nm"])
    return et, around, rad


def exclusion_spheres(ptypes, ppos, around, rad, glo, ghi):
    """(m, 5) spheres (x, y, z, r, rule) around pre-placed instances (types into preplaced libs) that reach box [glo, ghi]"""
    out = []
    for i in range(len(rad)):
        k = around[i][ptypes]
        if not k.any():
            continue
        P = np.asarray(ppos[k], np.float64)
        k2 = np.all((P > glo - rad[i]) & (P < ghi + rad[i]), axis=1)
        if k2.any():
            out.append(np.c_[P[k2], np.full(int(k2.sum()), rad[i]), np.full(int(k2.sum()), float(i))])
    return np.ascontiguousarray(np.concatenate(out)) if out else np.zeros((0, 5))


def exclusion_tile_spheres(Z, glo, ghi):
    """rows (x, y, z, r, rule) of provider zone spheres Z that reach the box [glo, ghi] (origin: yeast
    wc3d._tile_exclusion: centre within the box grown by the sphere radius)"""
    if not len(Z):
        return np.zeros((0, 5))
    k = np.all((Z[:, :3] > glo - Z[:, 3:4]) & (Z[:, :3] < ghi + Z[:, 3:4]), axis=1)
    return np.ascontiguousarray(Z[k])
