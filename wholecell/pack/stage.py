# origin: jurkat/yc/pipeline/jwc.py (legacy workspace, 2026-11-15)
# Changed: steps are functions of (cfg, bdir) writing to <bdir>/pack/; every JWC_* env switch and module constant is a
# `pack:` config value; regions / region map / area factors / thin lumens come from config (wholecell.pack.regions);
# virtual regions through registered providers (wholecell.pack.virtual); the face-actin and speckle chromatin drops
# are generic config-driven index drops (`pack.index_drops`); inputs come from the inventory pack table, the library
# manifest and the preplaced manifest of the build; seeds from wholecell.seeds (legacy v7 seeds via config);
# step outputs are written atomically and existing outputs are kept (tiles are resumable). Numerics unchanged.
# A7 (2026-11-16): yeast packer algorithms as `pack:` switches with v7 defaults (wholecell/pack/switches.py): size-class
# rounds, pilot thinning, shell cap, per-region budgets, squeeze, exclusion zones (origin comments at each function).
# A7 phase 2 (2026-11-17): organisation layer (P6, step `organisation`, wholecell/pack/organisation.py) and packing
# shapes (P8, pack.shapes) as switches with v7 defaults.
# A7 yeast side (2026-11-17): region tables in any order (kernel order + id maps, Regions.k_of / c_of); shells measured
# as volumes packed from their voxels; the pack evaluator uses the tiling box subset (pack.tiling.box_margin_nm);
# pre-pass index from a grouped npz (pack.preplaced_index, yeast prepass.npz); exclusion zones from providers; pilot
# squeeze and legacy_x5 re-split for thinning; presets (pack.preset / pack.presets, switches.effective).
"""Pack stage: python -m wholecell pack --cell C --build B --step {tiles,organisation,counts,index,pack,summary} [key=value]

  tiles    tile list (tile_nm cubes near the cell) + per-tile region volumes (meas_h_nm label voxels; virtual regions:
           selected volume) and membrane areas (meas_h_nm marching cubes per membrane SDF)    -> tiles.parquet
  organisation  (pack.organisation on) STRING network profiles, field loadings and partner CSR
           -> organisation/{network.npz, partners.npz, network.json}
  counts   expected copies per (tile, library, region): pack-table copies minus pre-pass consumption, mapped to regions
           (pack.region_map + mapping rules, area factors, open-face scaling) and split over tiles by a multinomial
           proportional to the tile measure (seed stream pack.counts)                 -> counts.npz, libs.json
  index    pre-pass instances and RNA beads bucketed per tile (index drops applied) -> preplaced_*.npy, mrna_*.npy
  pack     checkerboard phases (2x2x2 colours); each tile packed by pack3d.pack_tile (seed stream pack.tile, key = tile
           index) with obstacles = pre-pass + beads + already packed neighbour tiles  -> tiles/<i>.npz, pack_log.tsv
  summary  placed / expected mass per region                                          -> summary.json
Extra key=value arguments: workers=N (default: config pack.max_workers, else site max_workers); preset=<name>
(pack.presets entry, overrides pack.preset); for `pack`:
tiles=i,j,k | tiles=<file.npy> (tile subset), phases=0,1 (colour subset), box=x0,y0,z0,x1,y1,z1 (model nm);
pack_table=<tsv> (override of <bdir>/inventory/pack_table.tsv, development only).
Worker count never changes results: tiles of one colour are independent, colours are packed in order.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

_S = {}          # state shared with forked workers


# ------------------------------------------------------------------ config / inputs
class P:
    """pack parameters (config `pack:` with legacy v7 defaults)"""

    def __init__(self, cfg, bdir, kw):
        from .switches import effective
        pc = effective(cfg.get("pack") or {}, kw.get("preset"))
        self.cfg, self.pc, self.bdir, self.kw = cfg, pc, bdir, kw
        self.out = os.path.join(bdir, "pack")
        self.TILE = float(pc.get("tile_nm", 500.0))
        self.VOX = float(pc.get("vox_nm", 4.0))
        self.MESH_H = float(pc.get("mesh_h_nm", 8.0))
        self.MEAS_H = float(pc.get("meas_h_nm", 8.0))
        self.MEAS_VOX = float(pc.get("meas_vox_nm", self.MEAS_H))     # voxel size of the tile volume measures
        self.HALO = float(pc.get("halo_nm", 40.0))
        self.KEEP_F = float(pc.get("tile_keep_half_diag_f", 0.87))
        self.MAX_TRIES = int(pc.get("max_tries", 200))
        self.GIVE_UP = int(pc.get("give_up", 25))
        self.THIN_TRY_F = float(pc.get("thin_try_f", 10.0))
        self.THIN_GIVEUP_F = float(pc.get("thin_giveup_f", 16.0))
        self.SEED_PRUNE = bool(pc.get("seed_prune", False))
        self.LUMEN_ALIGN = bool(pc.get("lumen_align", True))
        self.SHELL_AREA_F = bool(pc.get("shell_area_factor", True))
        self.HMIN_NDIR = int(pc.get("seed_prune_ndir", 200)); self.HMIN_SEED = int(pc.get("seed_prune_seed", 5))
        self.MRNA_MEAN = float(pc.get("mrna_mean", 1500.0))
        site = cfg.get("_site") or {}
        self.workers = int(kw.get("workers") or pc.get("max_workers") or site.get("max_workers") or 8)

    def path(self, *a):
        return os.path.join(self.out, *a)

    def pack_table(self):
        return self.kw.get("pack_table") or self.pc.get("dev_pack_table") or os.path.join(self.bdir, "inventory", "pack_table.tsv")

    def library(self):
        """{species: npz path} from <bdir>/library/library_manifest.tsv (npz entries only)"""
        import pandas as pd
        f = os.path.join(self.bdir, "library", "library_manifest.tsv")
        L = pd.read_csv(f, sep="\t")
        L = L[L.path.str.endswith(".npz")]
        return dict(zip(L.species, L.path))

    def preplaced(self):
        f = os.path.join(self.bdir, "preplaced", "preplaced_manifest.json")
        return json.load(open(f))


def _setup(cfg, bdir, kw):
    from wholecell import registry
    registry.load_plugins(cfg)
    from . import geomload
    from .regions import Regions
    p = P(cfg, bdir, kw)
    g = geomload.load(cfg, bdir)
    reg = Regions(p.pc, g)
    return p, g, reg


def _virtual(p, reg, g):
    from .virtual import VirtualSet
    return VirtualSet(p.pc, reg, g)


def face_keep(tri, G, c3):
    """triangles (n, 9) whose centroid lies in the modelled domain (above the open face)"""
    if G[c3.G_FACE] <= 0 or len(tri) == 0:
        return tri
    return tri[(tri[:, 2] + tri[:, 5] + tri[:, 8]) / 3.0 >= G[c3.G_ZFACE]]


def _atomic(path, saver, *a, **k):
    from wholecell.build import atomic_save
    os.makedirs(os.path.dirname(path), exist_ok=True)
    atomic_save(path, saver, *a, **k)


def _json_save(path, obj, **k):
    def sv(tmp):
        with open(tmp, "w") as f:
            json.dump(obj, f, **k)
    _atomic(path, sv)


def _npy_save(path, arr):
    def sv(tmp):
        with open(tmp, "wb") as f:
            np.save(f, arr)
    _atomic(path, sv)


def _done(*paths):
    if all(os.path.exists(x) for x in paths):
        print("outputs exist, step skipped:", ", ".join(os.path.basename(x) for x in paths), flush=True)
        return True
    return False


# ------------------------------------------------------------------ tiles + measures (wholecell/pack/tiling.py)
TILING_KEYS = ("select", "sample_reach_nm", "box_margin_nm", "pore_reach_nm", "labvol", "meshes", "mesh_membranes")


def step_tiles(cfg, bdir, kw):
    import pandas as pd
    from . import tiling
    p, g, reg = _setup(cfg, bdir, kw)
    if _done(p.path("tiles.parquet"), p.path("tiles_meta.json")):
        return
    tc = dict(p.pc.get("tiling") or {})
    bad = [k for k in tc if k not in TILING_KEYS]
    if bad:
        raise ValueError(f"pack.tiling: unknown key(s) {bad} (known: {', '.join(TILING_KEYS)})")
    sel = tc.get("select", "centre_depth")
    cell, c3, TILE = g.cell, g.engine, p.TILE
    if tc.get("box_margin_nm") is not None and not hasattr(cell, "prim_bsphere"):
        raise ValueError("pack.tiling.box_margin_nm needs prim bounding spheres (prims.npz 'bsphere') in the package")
    if sel == "centre_depth":
        outer = reg.mem_codes.get("CELL_OUTER", c3.M_OUTER)
        T, lo0, nt = tiling.select_centre_depth(cell, c3, TILE, p.KEEP_F, outer)
    elif sel == "enclosing_sphere":
        T, lo0, nt = tiling.select_enclosing_sphere(cell, c3, TILE, float(tc.get("sample_reach_nm", 88.0)))
    else:
        raise ValueError(f"pack.tiling.select: {sel!r} (centre_depth | enclosing_sphere)")
    print(f"candidate tiles {nt ** 3:,} -> {len(T):,} near the cell ({sel})", flush=True)
    c2n = getattr(cell, "code_to_name", None) or dict(enumerate(c3.LABEL_NAMES))     # engine label code -> name
    nlab = 1 + max(c2n)
    _S.update(p=p, reg=reg, cell=cell, c3=c3, T=T, lo0=lo0, vset=_virtual(p, reg, g))
    tiling._S.update(_S, tcfg=tc, nlab=nlab)
    import multiprocessing as mp
    chunks = np.array_split(np.arange(len(T)), p.workers * 8)
    t0 = time.time()
    with mp.get_context("fork").Pool(p.workers) as pool:
        res = [r for part in pool.imap_unordered(tiling.measure_tiles, chunks) for r in part]
    res.sort(key=lambda r: r[0])
    meas = reg.measures()
    cols = [c for _, c, _, _ in meas]
    vcols = [c for _, c, k, _ in meas if k == "volume"]; acols = [c for _, c, k, _ in meas if k == "area"]
    M = pd.DataFrame([r[:1 + len(cols)] for r in res], columns=["tile"] + cols)
    M[["ti", "tj", "tk"]] = T[M.tile.values]
    keep = np.ones(len(M), bool)
    if sel == "centre_depth":                            # v7: tiles without cell content are dropped
        keep = ((M[vcols].sum(axis=1) > 0) | (M[acols].sum(axis=1) > 0)).values
    M = M[keep].reset_index(drop=True)
    M["tile"] = np.arange(len(M))
    _atomic(p.path("tiles.parquet"), lambda tmp: M.to_parquet(tmp))
    if tc.get("labvol"):
        names = [c2n.get(k, str(k)) for k in range(nlab)]
        LV = pd.DataFrame(np.array([r[-1] for r in res])[keep], columns=names)
        _atomic(p.path("tile_labvol.parquet"), lambda tmp: LV.to_parquet(tmp))
    _json_save(p.path("tiles_meta.json"), dict(lo0=lo0.tolist(), TILE=TILE, **({} if sel == "centre_depth" else
                                                                             dict(select=sel, columns=cols))))
    print(f"{len(M):,} tiles in {time.time() - t0:.0f}s; volumes um3:", (M[vcols].sum() / 1e9).round(3).to_dict(),
          "| areas um2:", (M[acols].sum() / 1e6).round(2).to_dict(), flush=True)


# ------------------------------------------------------------------ organisation layer (P6)
def step_organisation(cfg, bdir, kw):
    from . import switches, organisation
    p = P(cfg, bdir, kw)
    sw = switches.load(p.pc)
    o = sw.organisation
    if not o["enabled"]:
        print("pack.organisation is off: nothing to do", flush=True)
        return
    out = p.path("organisation")
    if _done(os.path.join(out, "network.npz"), os.path.join(out, "partners.npz"), os.path.join(out, "network.json")):
        return
    LIB = p.library()
    ids = organisation.network_ids(cfg, bdir, o, p.pc, dict(pack_table=p.pack_table(), library=LIB))
    t0 = time.time()
    organisation.build(cfg, bdir, o, ids, p.pack_table(), out, log=lambda *a: print(*a, flush=True))
    print(f"organisation layer in {time.time() - t0:.0f}s -> {out}", flush=True)


def _partners(p, sw, L, pl):
    """pack.organisation (P6), small scale: partner CSR over the packed libraries L, kernel parameters, and the map of
    pre-pass library indices into L (pre-placed instances of packed types are anchors too; -1 = not an anchor)"""
    from .organisation import remap_partners, params
    f = p.path("organisation", "partners.npz")
    if not os.path.exists(f):
        raise FileNotFoundError(f"pack.organisation is on but {f} is missing: run `pack --step organisation`")
    z = np.load(f, allow_pickle=True)
    pr = remap_partners({k: z[k] for k in ("off", "idx", "strength", "kind")}, z["ids"], L)
    pos = {u: i for i, u in enumerate(L)}
    p2l = np.array([pos.get(u, -1) for u in pl], np.int64)
    print(f"organisation partners: {int((np.diff(pr['off']) > 0).sum())} of {len(L)} packed types with partners, "
          f"{len(pr['idx']):,} links ({int((pr['kind'] == 0).sum()):,} contact); {int((p2l >= 0).sum())} of {len(pl)} pre-pass "
          f"species are anchors", flush=True)
    return pr, params(sw.organisation), p2l


# ------------------------------------------------------------------ expected counts
def _area_factors(p, g, carea):
    fac = {}
    for a in p.pc.get("area_factors") or []:
        mem, F = a["membrane"], float(a["factor"])
        fac[mem] = F
        ms = a.get("measured")
        if ms and g.flag(ms["when_flag"]) and g.engine.n_fields(g.cell.G):
            if ms["mode"] == "ref_membrane":
                fac[mem] = max(1.0, F * carea[ms["ref"]] / max(carea[mem], 1.0))
            elif ms["mode"] == "summary_area_um2":
                v = g.summary
                for k in ms["key"]:
                    v = v[k]
                fac[mem] = max(1.0, F * v * 1e6 / max(carea[mem], 1.0))
            else:
                raise ValueError(f"area factor mode {ms['mode']}")
    return fac


# ------------------------------------------------------------------ copy ledger (no silent losses)
# Every copy of the pack table that is not requested from the packer is a row of pack/copy_ledger.tsv (counts step)
# or pack/index_drops.tsv (index step) with a reason. Rows marked unexplained (species without a structure, unless a
# reference preset declares pack.ledger.explain_no_structure) fail the counts step above
# pack.ledger.max_unexplained_mass_frac of the pack-table mass. Copies moved, not lost (an assembly without a structure
# split into its subunits, index-dropped pre-pass instances returned to the pool) are recorded as well.
LEDGER_COLS = ("event", "ing_id", "copies", "mass_pg", "explained", "reason")
AVOGADRO = 6.022e23


class CopyLedger:
    def __init__(self, p):
        self.p, self.rows = p, []
        self.cfg = p.pc.get("ledger") or {}

    def add(self, event, ing_id, copies, mass_da, explained, reason):
        self.rows.append((event, str(ing_id), float(copies), float(copies) * float(mass_da) / AVOGADRO * 1e12,
                          bool(explained), reason))

    def write(self, total_mass_pg, log=print):
        """pack/copy_ledger.tsv + the check; returns the summary dict (libs.json key `ledger`)"""
        import pandas as pd
        L = pd.DataFrame(self.rows, columns=list(LEDGER_COLS))
        _atomic(self.p.path("copy_ledger.tsv"), lambda tmp: L.to_csv(tmp, sep="\t", index=False))
        lost = L[L.event.isin(LOSS_EVENTS)]
        un = lost[~lost.explained]
        lim = float(self.cfg.get("max_unexplained_mass_frac", 1e-4))
        frac = float(un.mass_pg.sum() / max(total_mass_pg, 1e-30))
        out = dict(rows=len(L), by_event={k: dict(rows=int(len(g)), copies=round(float(g.copies.sum()), 1),
                                                   mass_pg=round(float(g.mass_pg.sum()), 6)) for k, g in L.groupby("event")},
                   lost_explained_pg=round(float(lost[lost.explained].mass_pg.sum()), 6),
                   lost_unexplained_pg=round(float(un.mass_pg.sum()), 6), pack_table_pg=round(total_mass_pg, 4),
                   unexplained_mass_frac=frac, max_unexplained_mass_frac=lim)
        log(f"copy ledger: {out['by_event']}; unexplained {frac:.2e} of the pack-table mass (limit {lim:g})")
        if frac > lim:
            raise RuntimeError(f"pack counts: unexplained copy losses {frac:.2e} of the pack-table mass > pack.ledger."
                               f"max_unexplained_mass_frac {lim:g}: {un.sort_values('mass_pg').ing_id.tolist()[-20:]} "
                               f"(see {self.p.path('copy_ledger.tsv')})")
        return out


LOSS_EVENTS = ("no_structure", "not_requested")


def _decompose(p, I, LIB, led, log=print):
    """assemblies without a structure -> their subunit monomers (library/decomposed.tsv, library step `fallback` with
    decompose): copies x stoichiometry added to the monomer rows (comp_weights mixed copy-weighted; a monomer that is not
    a pack-table species gets a row with the assembly's compartments)"""
    import pandas as pd
    f = os.path.join(p.bdir, "library", "decomposed.tsv")
    if not os.path.exists(f):
        return I
    D = pd.read_csv(f, sep="\t", keep_default_na=False)
    if not len(D):
        return I
    I = I.copy()
    for parent, grp in D.groupby("parent", sort=False):
        ii = I.index[I.ing_id == parent]
        if not len(ii) or int(I.at[ii[0], "copies"]) <= 0 or parent in LIB:
            continue
        i = ii[0]; C = int(I.at[i, "copies"]); w = json.loads(I.at[i, "comp_weights"])
        for r in grp.itertuples():
            add = C * int(r.n)
            jj = I.index[I.ing_id == r.subunit]
            if len(jj):
                j = jj[0]; c0 = int(I.at[j, "copies"]); w0 = json.loads(I.at[j, "comp_weights"])
                keys = list(dict.fromkeys(list(w0) + list(w)))
                I.at[j, "comp_weights"] = json.dumps({k: (c0 * w0.get(k, 0.0) + add * w.get(k, 0.0)) / (c0 + add) for k in keys})
                I.at[j, "copies"] = c0 + add
            else:
                tm = str(r.transmembrane) in ("True", "true", "1")
                row = I.loc[i].to_dict()
                row.update(ing_id=r.subunit, kind="protein", genes=r.gene, stoich="1", copies=add, mass_da=float(r.mass_da),
                           transmembrane=tm, n_tm=int(tm), struct_src=f"decomposed:{parent}", uniprot=r.uniprot)
                if "mass_pg" in row:
                    row["mass_pg"] = add * float(r.mass_da) / AVOGADRO * 1e12
                I = pd.concat([I, pd.DataFrame([row])], ignore_index=True)
            led.add("decomposed_in", r.subunit, add, r.mass_da, True, f"subunit of {parent} (no structure of the assembly)")
        led.add("decomposed_out", parent, C, I.at[i, "mass_da"], True, "assembly without a structure: copies moved to its subunits")
        I.at[i, "copies"] = 0
        log(f"  {parent}: {C:,} copies without a structure returned to {len(grp)} subunit monomers", flush=True)
    return I


def _index_returns(p, I, led, log=print):
    """pre-pass instances dropped at indexing that pack.index_drops sends back to the pool (return_to_pool: {library id:
    pack-table species}): their copies are added again (pack/index_drops.tsv of the index step)"""
    import pandas as pd
    ret = {d["name"]: d["return_to_pool"] for d in p.pc.get("index_drops") or [] if d.get("return_to_pool")}
    if not ret:
        return I
    f = p.path("index_drops.tsv")
    if not os.path.exists(f):
        raise FileNotFoundError(f"pack.index_drops {sorted(ret)} return copies to the pool: run `pack --step index` before counts ({f})")
    D = pd.read_csv(f, sep="\t")
    I = I.copy()
    for r in D[D.name.isin(ret)].itertuples():
        tgt = ret[r.name].get(r.lib)
        if tgt is None:
            continue
        ii = I.index[I.ing_id == tgt]
        if not len(ii):
            raise ValueError(f"pack.index_drops {r.name}: return_to_pool target {tgt} is not a pack-table species")
        I.at[ii[0], "copies"] = int(I.at[ii[0], "copies"]) + int(r.n)
        led.add("index_drop_returned", tgt, r.n, I.at[ii[0], "mass_da"], True, f"{r.name}: {r.reason} (returned to the pool)")
        log(f"  {tgt}: {int(r.n):,} copies of index drop '{r.name}' ({r.lib}) returned to the pool", flush=True)
    return I


def _ledger_x5(cfg, bdir, p, reg):
    """copy ledger of pack.requests method legacy_x5: per-type totals whose region is not a pack region are not requested
    (the yeast engine labels them, e.g. no_structure, dropped_TM_in_soluble_compartment); recorded with that label as
    the reason (explained: a stated engine rule), except no_structure (unexplained unless pack.ledger.explain_no_structure)"""
    import pandas as pd
    from .requests import _path
    rq = p.pc["requests"]
    exp = pd.read_parquet(_path(bdir, rq["totals"]))
    pt = pd.read_csv(p.pack_table(), sep="\t", usecols=["ing_id", "copies", "mass_da"])
    mass = dict(zip(pt.ing_id.astype(str), pt.mass_da.astype(float)))
    led = CopyLedger(p)
    why = (p.pc.get("ledger") or {}).get("explain_no_structure")
    off = exp[~exp.region.isin(reg.RNAMES)].groupby(["region", "ing_id"], as_index=False).expected.sum()
    for r in off.itertuples():
        unexpl = r.region == "no_structure" and not why
        led.add("not_requested", r.ing_id, r.expected, mass.get(str(r.ing_id), 0.0), not unexpl,
                (why if r.region == "no_structure" and why else f"engine region label '{r.region}'"))
    return led.write(float((pt.copies * pt.mass_da).sum() / AVOGADRO * 1e12))


def step_counts(cfg, bdir, kw):
    import pandas as pd
    from wholecell import seeds
    p, g, reg = _setup(cfg, bdir, kw)
    if _done(p.path("counts.npz"), p.path("libs.json")):
        return
    method = (p.pc.get("requests") or {}).get("method", "multinomial")
    if method == "legacy_x5":                             # yeast X5 request split (wholecell/pack/requests.py)
        from . import requests, switches
        d, info = requests.legacy_x5(cfg, bdir, p, reg, pd.read_parquet(p.path("tiles.parquet")),
                                     json.load(open(p.path("tiles_meta.json")))["lo0"], sw=switches.load(p.pc, reg.RNAMES),
                                     log=lambda *a: print(*a, flush=True))
        info["ledger"] = _ledger_x5(cfg, bdir, p, reg)
        _atomic(p.path("counts.npz"), lambda tmp: np.savez_compressed(tmp, **{k: d[k].values for k in ("tile", "lib", "region", "n")}))
        _json_save(p.path("libs.json"), info)
        print(f"counts (legacy_x5): {len(d):,} (tile, lib, region) entries, {int(d.n.sum()):,} instances", flush=True)
        return
    if method != "multinomial":
        raise ValueError(f"pack.requests.method {method!r} (multinomial | legacy_x5)")
    reg.require_core_packer("the v7 counts step")
    RNAMES, VOL_REG, MEM_REG = reg.RNAMES, reg.VOL_REG, reg.MEM_REG
    M = pd.read_parquet(p.path("tiles.parquet"))
    ING = p.pack_table()
    I = pd.read_csv(ING, sep="\t")
    man = p.preplaced()
    pre = man.get("instances", [])
    consumed = {}
    for m in pre:
        for k, v in json.load(open(m["json"])).get("consumed", {}).items(): consumed[k] = consumed.get(k, 0) + v
    ckey = {}                                             # pre-pass consumption taken out of one key (lamina lattice)
    for m in pre:
        js = json.load(open(m["json"]))
        if js.get("consumed_key"):
            for k, v in js["consumed"].items(): ckey[k] = (js["consumed_key"], v)
    for k, (key, c_) in ckey.items():
        ii = I.index[I.ing_id == k]
        if not len(ii): continue
        i = ii[0]; w = json.loads(I.at[i, "comp_weights"]); C_ = float(I.at[i, "copies"]); left = C_ - c_
        if left <= 0 or key not in w: continue
        nw = {kk: C_ * vv / left for kk, vv in w.items() if kk != key}; nw[key] = max(C_ * w[key] - c_, 0.0) / left
        I.at[i, "comp_weights"] = json.dumps(nw)
        print(f"  {k}: {c_:,} pre-pass copies taken from key '{key}' ({C_ * w[key] - c_:,.0f} left there)", flush=True)
    total_pg = float((I.copies * I.mass_da).sum() / AVOGADRO * 1e12)
    I["copies"] = (I.copies - I.ing_id.map(consumed).fillna(0)).clip(lower=0).round().astype(np.int64)
    LIB = p.library()
    led = CopyLedger(p)
    I = _index_returns(p, I, led)
    I = _decompose(p, I, LIB, led)
    if p.pc.get("libs_from"):                              # keep an earlier library order (tiles store indices into libs.json)
        order = json.load(open(p.pc["libs_from"]))["libs"]
        I = I.set_index("ing_id").loc[order].reset_index()
    else:
        why = (p.pc.get("ledger") or {}).get("explain_no_structure")
        for r in I[~I.ing_id.isin(set(LIB)) & (I.copies > 0)].itertuples():   # no library entry: not packed, recorded
            led.add("no_structure", r.ing_id, r.copies, r.mass_da, bool(why), why or "no library entry (see library/species_ledger.tsv)")
        I = I[I.ing_id.isin(set(LIB)) & (I.copies > 0)].reset_index(drop=True)
    ledger = led.write(total_pg)
    memb = np.array([("membrane" in np.load(LIB[i]).files) for i in I.ing_id])
    cvol = {r: float(M[r].sum()) for r, _ in VOL_REG}
    carea = {m: float(M["mem_" + m].sum()) for m, _ in MEM_REG}
    fac = _area_factors(p, g, carea)
    print("membrane area factors (inventory / modelled):", {k: round(v, 2) for k, v in fac.items()}, flush=True)
    ofc = p.pc.get("open_face") or {}
    of = g.summary.get(ofc.get("summary_key", "open_face")) if (ofc and g.flag(ofc.get("when_flag", "G_FACE"))) else None
    f_vol, f_pm = (of["f_vol"], of["f_pm"]) if of else (1.0, 1.0)
    pm_mems = set(ofc.get("membranes", []))
    pm_regs = {m for m, _ in MEM_REG if m in pm_mems} | {s for s, mn in reg.SHELL_MEMNAME.items() if mn in pm_mems}
    print(f"open-face scaling of cell copy numbers: f_vol {f_vol}, f_pm {f_pm} (PM regions {sorted(pm_regs)})", flush=True)
    seed = seeds.seed_int(cfg, "pack.counts")
    rng = np.random.default_rng(seed); rows = []
    tot = {}
    for li, r in I.iterrows():
        for comp, wt in json.loads(r.comp_weights).items():
            if memb[li] and r.transmembrane:
                rg = reg.map_tm(comp, r)
                share = 1.0 / fac.get(rg, 1.0)              # membrane not fully built: only the modelled share
            else:
                rg = reg.map_soluble(comp, r)
                share = 1.0 / fac.get(reg.SHELL_MEMNAME.get(rg), 1.0) if p.SHELL_AREA_F else 1.0
            share *= f_pm if rg in pm_regs else f_vol
            key = (li, RNAMES.index(rg))
            tot[key] = tot.get(key, 0.0) + r.copies * wt * share
    from . import switches
    sw = switches.load(p.pc, RNAMES)
    shell_cap_rep = None
    if sw.shell_cap["enabled"]:                           # P3 (pack.shell_cap): before the stochastic rounding
        tot, shell_cap_rep = _shell_cap(tot, I, LIB, reg, M, _shapes(p, sw, list(I.ing_id)))
    field = _organisation_field(cfg, p, sw, reg, M, I, LIB, tot) if sw.organisation["enabled"] else None
    pcache = {}
    # freeze fix 2026-11-19: the tile split of a field-capped region is drawn from its own stream (pack.counts.field_cap);
    # stream pack.counts still makes the uncapped draw (discarded), so every other region keeps its draws bitwise
    capped = field[3] if field is not None else {}
    rng_cap = seeds.rng(cfg, "pack.counts.field_cap") if capped else None
    for (li, ri), x in tot.items():
        n = int(x) + (rng.random() < x - int(x))
        if n <= 0: continue
        if ri not in pcache:
            w = M[reg.tile_column(ri)].values.astype(np.float64); s = w.sum()
            pcache[ri] = (np.nonzero(w)[0], w[w > 0] / s) if s > 0 else (None, None)
        nz, pr = pcache[ri]
        if nz is None: continue
        if field is not None and li in field[0]:          # pack.organisation: x exp(log weight) x crowding correction
            lw = field[0][li][nz]
            wf = M[reg.tile_column(ri)].values.astype(np.float64)[nz] * np.exp(lw - lw.max()) * (field[1][ri][nz] if ri in field[1] else 1.0)
            pr = wf / wf.sum()
        k = rng.multinomial(n, pr)
        if ri in capped:                                  # field cap (unmodulated types: the same measure split)
            if li in field[0]:
                wc = M[reg.tile_column(ri)].values.astype(np.float64)[nz] * np.exp(lw - lw.max()) * capped[ri][nz]
                pr = wc / wc.sum()
            k = rng_cap.multinomial(n, pr)
        sel = k > 0
        rows.append(np.stack([nz[sel], np.full(sel.sum(), li), np.full(sel.sum(), ri), k[sel]], 1))
    C = np.concatenate(rows).astype(np.int64)
    C = C[np.lexsort((C[:, 1], C[:, 0]))]
    _atomic(p.path("counts.npz"), lambda tmp: np.savez_compressed(tmp, tile=C[:, 0].astype(np.int32), lib=C[:, 1].astype(np.int32),
                                                                   region=C[:, 2].astype(np.int16), n=C[:, 3].astype(np.int32)))
    _json_save(p.path("libs.json"), dict(libs=list(I.ing_id), mass_da=I.mass_da.tolist(), consumed=consumed, f_vol=f_vol, f_pm=f_pm,
                                         area_factors=fac, rnames=list(RNAMES), ingredients=ING, geometry=g.source,
                                         prepass=[m["npz"] for m in pre], counts_seed=seed, ledger=ledger,   # output key kept at A0c
                                         **({"shell_cap": shell_cap_rep} if shell_cap_rep is not None else {}),
                                         **({"organisation": dict(network=p.path("organisation", "network.npz"),
                                                                  field_seed=seeds.seed_int(cfg, "pack.organisation.field"),
                                                                  types_modulated=len(field[0]),
                                                                  **({"field_cap": field[2]} if field[2] is not None else {}))}
                                            if field is not None else {})))
    print(f"counts: {len(C):,} (tile, lib, region) entries, {C[:, 3].sum():,} instances; per region:",
          {RNAMES[ri]: int(C[C[:, 2] == ri, 3].sum()) for ri in range(len(RNAMES))}, flush=True)


def _organisation_field(cfg, p, sw, reg, M, I, LIB, tot):
    """pack.organisation (P6), large scale: ({lib index: log weight per tile}, {region index: crowding factor per tile})
    from the field loadings of pack/organisation/network.npz (origin: yeast wc3d.distribute, organisation.tile_logmod)"""
    import pandas as pd
    from wholecell import seeds
    from .requests import tile_logmod, crowding_factors
    o = sw.organisation
    f = p.path("organisation", "network.npz")
    if not os.path.exists(f):
        raise FileNotFoundError(f"pack.organisation is on but {f} is missing: run `pack --step organisation` before counts")
    z = np.load(f, allow_pickle=True)
    net = [str(x) for x in z["ids"]]; pos = {u: li for li, u in enumerate(I.ing_id)}
    keep = [k for k, i in enumerate(z["tix"]) if net[i] in pos]
    tix = np.array([pos[net[z["tix"][k]]] for k in keep], np.int64)
    E = np.asarray(z["E"], np.float64)[keep]
    centres = np.asarray(json.load(open(p.path("tiles_meta.json")))["lo0"]) + (M[["ti", "tj", "tk"]].values + 0.5) * p.TILE
    logmod = tile_logmod(list(range(len(I))), tix, E, centres, float(o["beta"]), seeds.seed_int(cfg, "pack.organisation.field"),
                         float(o["field_l_nm"]), int(o["n_features"]))
    shapes = _shapes(p, sw, list(I.ing_id))
    vol = {}
    for li, i in enumerate(I.ing_id):
        if li in logmod:
            vol[li] = float((4 / 3 * np.pi * _tree(i, LIB, shapes)[1] ** 3).sum())
    e = pd.DataFrame([(li, reg.RNAMES[ri], x) for (li, ri), x in tot.items()], columns=["lib_id", "region", "expected"])
    e = e.groupby(["lib_id", "region"], as_index=False).expected.sum()
    Mx = np.stack([M[reg.tile_column(ri)].values.astype(np.float64) for ri in range(len(reg.RNAMES))], 1)
    cf = crowding_factors(e, Mx, {n: i for i, n in enumerate(reg.RNAMES)}, logmod, vol)
    print(f"organisation field: {len(logmod)} of {len(I)} types modulated over {len(M)} tiles (beta {o['beta']}, "
          f"{E.shape[1]} components, {o['field_l_nm']} nm); crowding correction in {len(cf)} regions", flush=True)
    cap_rep = None
    if o["field_cap_regions"]:                            # freeze fix 2026-11-19 (requests.field_cap)
        from .requests import field_cap
        vol_names = {r for r, _ in reg.VOL_REG}
        bad = [r for r in o["field_cap_regions"] if r not in vol_names]
        if bad:
            raise ValueError(f"config pack.organisation.field_cap_regions: {bad} are not volume regions")
        cf2, cap_rep = field_cap(e, Mx, {n: i for i, n in enumerate(reg.RNAMES)}, logmod, cf,
                                 dict(enumerate(I.mass_da.astype(np.float64))), o["field_cap_regions"])
        print(f"organisation field cap (mass per tile <= the unmodulated split): {cap_rep}", flush=True)
        return logmod, cf, cap_rep, {ri: cf2[ri] for ri in cf2 if reg.RNAMES[ri] in cap_rep}
    return logmod, cf, cap_rep, {}


def _shell_cap(tot, I, LIB, reg, M, shapes=None):
    """pack.shell_cap (P3): shells ask for at most the requested density of the bulk region on their side"""
    from .algorithms import shell_cap
    vol_lib, br_lib = np.zeros(len(I)), np.zeros(len(I))
    for li, i in enumerate(I.ing_id):
        c, r = _tree(i, LIB, shapes)
        vol_lib[li] = float((4.0 / 3.0 * np.pi * r ** 3).sum()); br_lib[li] = float((np.linalg.norm(c, axis=1) + r).max())
    rid = {n: i for i, n in enumerate(reg.RNAMES)}            # config region ids (counts keys)
    vol_regions = {rid[name]: labs for name, labs in reg.VOL_REG if reg.spec[name]["kind"] == "volume"}
    region_volume = {rid[name]: float(M[name].sum()) for name, _ in reg.VOL_REG}
    shells = [dict(ri=rid[name], area_nm2=float(M["mem_" + reg.SHELL_MEMNAME[name]].sum()), d0=d0, d1=d1,
                   add_radius=bool(rad & 1), pore_limit=rad >> 1, labels=labs)
              for name, m, labs, side, d0, d1, rad in reg.SHELL_REG]
    for name, labs in reg.VOL_REG:                            # volume-measured shells: measured volume
        x = reg.restriction(name) if reg.spec[name]["kind"] == "shell" else None
        if x is not None:
            shells.append(dict(ri=rid[name], area_nm2=0.0, d0=x[2], d1=x[3], add_radius=False, pore_limit=0, labels=labs,
                               volume_nm3=float(M[name].sum())))
    tot, rep = shell_cap(tot, vol_lib, br_lib, vol_regions, region_volume, shells)
    print(f"shell cap: {len(rep)} shells capped", flush=True)
    return tot, {reg.RNAMES[ri]: dict(factor=round(f, 4), moved=round(mv, 1), to=reg.RNAMES[t]) for ri, (f, mv, t) in rep.items()}


# ------------------------------------------------------------------ pre-pass index
def _lib_match(libs, d):
    ex = set(d.get("libs", [])); pf = tuple(d.get("lib_prefixes", []))
    return np.array([(l in ex) or (bool(pf) and l.startswith(pf)) for l in libs], bool)


def _index_drop(d, libs, T_, P_, g):
    """boolean drop mask over the indexed pre-pass instances for one `pack.index_drops` entry"""
    c3, cell = g.engine, g.cell
    drop = np.zeros(len(T_), bool)
    if d.get("when_flag") and not g.flag(d["when_flag"]):
        return drop
    sel = _lib_match(libs, d)[T_]
    mode = d["mode"]
    if mode == "face_band":                 # within z_above_face_nm of the open face and deeper than min_sdf_nm below a membrane
        near = sel & (P_[:, 2] < cell.G[c3.G_ZFACE] + float(d["z_above_face_nm"]))
        if near.any():
            from .regions import membrane_map
            dep = cell.mems(P_[near].astype(np.float64))[:, membrane_map(g)[d["membrane"]]]
            drop[np.nonzero(near)[0][dep > float(d["min_sdf_nm"])]] = True
    elif mode == "g_spheres":               # inside the spheres of a G section (n, then n x (x, y, z, r)), e.g. speckles
        from scipy.spatial import cKDTree
        o_ = int(cell.G[getattr(c3, d["section_slot"])]); n_ = int(cell.G[o_]); S_ = cell.G[o_ + 1:o_ + 1 + 4 * n_].reshape(n_, 4)
        d_, k_ = cKDTree(S_[:, :3]).query(P_[sel].astype(np.float64), k=1)
        drop[np.nonzero(sel)[0][d_ < S_[k_, 3]]] = True
    elif mode == "inside_labels":           # centre label in the given label names
        from .regions import label_map
        lm = label_map(g); codes = [lm[x] for x in d["labels"]]
        ii = np.nonzero(sel)[0]
        if len(ii):
            drop[ii[np.isin(cell.labels(P_[ii].astype(np.float64)), codes)]] = True
    else:
        raise ValueError(f"index drop mode {mode}")
    return drop


INDEX_DROP_KEYS = ("name", "libs", "lib_prefixes", "when_flag", "mode", "z_above_face_nm", "membrane", "min_sdf_nm",
                   "section_slot", "labels", "reason", "return_to_pool")


def _index_drops(p, libs, T_, P_, Q_, g):
    """apply pack.index_drops in order; every entry states its `reason` (and optionally return_to_pool {library id:
    pack-table species}, read by the counts step); pack/index_drops.tsv lists the dropped instances per (drop, library)"""
    import pandas as pd
    rows = []
    for d in p.pc.get("index_drops") or []:
        bad = [k for k in d if k not in INDEX_DROP_KEYS]
        if bad:
            raise ValueError(f"config pack.index_drops {d.get('name')}: unknown key(s) {bad} (known: {', '.join(INDEX_DROP_KEYS)})")
        if not d.get("reason"):
            raise ValueError(f"config pack.index_drops {d.get('name')}: `reason` missing (every dropped copy is recorded with a reason)")
        drop = _index_drop(d, libs, T_, P_, g)
        print(f"index drop '{d['name']}': {int(drop.sum()):,} pre-pass instances dropped", flush=True)
        if drop.any():
            ids, n = np.unique(T_[drop], return_counts=True)
            ret = d.get("return_to_pool") or {}
            for t, k in zip(ids, n):
                rows.append((d["name"], libs[int(t)], int(k), ret.get(libs[int(t)], ""), d["reason"]))
            T_, P_, Q_ = T_[~drop], P_[~drop], Q_[~drop]
    D = pd.DataFrame(rows, columns=["name", "lib", "n", "returned_to", "reason"])
    D = D.groupby(["name", "lib", "returned_to", "reason"], as_index=False, sort=False).n.sum()[["name", "lib", "n", "returned_to", "reason"]]
    _atomic(p.path("index_drops.tsv"), lambda tmp: D.to_csv(tmp, sep="\t", index=False))
    return T_, P_, Q_


def _index_outside(p, libs, T_, keep):
    """pre-pass instances whose centre lies in no tile (outside the packed volume): appended to pack/index_drops.tsv"""
    import pandas as pd
    if keep.all():
        return
    ids, n = np.unique(T_[~keep], return_counts=True)
    D = pd.read_csv(p.path("index_drops.tsv"), sep="\t", keep_default_na=False)
    D = pd.concat([D, pd.DataFrame([("outside_tiles", libs[int(t)], int(k), "", "centre outside every tile (outside the packed volume)")
                                    for t, k in zip(ids, n)], columns=D.columns)], ignore_index=True)
    _atomic(p.path("index_drops.tsv"), lambda tmp: D.to_csv(tmp, sep="\t", index=False))


def _tile_ids(tij, M):
    dims = M[["ti", "tj", "tk"]].values.max(0) + 2
    lut = -np.ones(dims, np.int64); lut[M.ti, M.tj, M.tk] = np.arange(len(M))
    ok = np.all((tij >= 0) & (tij < dims), axis=1)
    out = -np.ones(len(tij), np.int64); out[ok] = lut[tij[ok, 0], tij[ok, 1], tij[ok, 2]]
    return out


PREPLACED_INDEX_KEYS = ("npz", "lib_ids", "groups", "beads", "bead_radius_nm", "group_names")


def _index_grouped(p, cfg, bdir, M, lo0):
    """pack.preplaced_index (A7 yeast side): pre-placed instances from ONE grouped npz (yeast s3 prepass.npz layout: per
    group k arrays k_pos, k_q and k_type, or one library index in the scalar named by the group's second entry), library
    indices into `lib_ids` (npy); RNA / DNA bead arrays named in `beads` (radius bead_radius_nm). Index entries = one
    per (group, library) pair (origin: plugins/yeast/pack_view.py, wc3d.pre_instances order); preplaced_groups.json gives
    each entry's group name (group_names map, default config export.preplaced_group_regions, else the group key).
    Returns (libs, groups, T_, P_, Q_, bead points or None, bead radius)."""
    from .requests import _path
    ci = p.pc["preplaced_index"]
    bad = [k for k in ci if k not in PREPLACED_INDEX_KEYS]
    if bad:
        raise ValueError(f"pack.preplaced_index: unknown key(s) {bad} (known: {', '.join(PREPLACED_INDEX_KEYS)})")
    pre = np.load(_path(bdir, ci["npz"]), allow_pickle=False)
    lib_ids = [str(x) for x in np.load(_path(bdir, ci["lib_ids"]), allow_pickle=True)]
    gmap = ci.get("group_names") or (cfg.get("export") or {}).get("preplaced_group_regions") or {}
    parts = [(k, li) for k, li in ci["groups"] if f"{k}_pos" in pre.files]
    if not parts:
        raise ValueError(f"pack.preplaced_index: none of the groups is in {ci['npz']}")
    gt = np.concatenate([pre[f"{k}_type"].astype(np.int64) if f"{k}_type" in pre.files else
                         np.full(len(pre[f"{k}_pos"]), int(pre[li]), np.int64) for k, li in parts])
    grp = np.concatenate([np.full(len(pre[f"{k}_pos"]), gi, np.int64) for gi, (k, _) in enumerate(parts)])
    P_ = np.concatenate([pre[f"{k}_pos"] for k, _ in parts]).astype(np.float32)
    Q_ = np.concatenate([pre[f"{k}_q"] for k, _ in parts]).astype(np.float32)
    uniq, inv = np.unique(np.c_[grp, gt], axis=0, return_inverse=True)
    libs = [lib_ids[int(t)] for _, t in uniq]
    groups = [gmap.get(parts[int(g_)][0], parts[int(g_)][0]) for g_, _ in uniq]
    print(f"pre-pass groups ({ci['npz']}): " + ", ".join(f"{k} {len(pre[f'{k}_pos']):,}" for k, _ in parts), flush=True)
    bk = [k for k in (ci.get("beads") or []) if k in pre.files and len(pre[k])]
    pts = np.concatenate([pre[k].astype(np.float32) for k in bk]) if bk else None
    if bk:
        print("RNA / DNA beads: " + ", ".join(f"{k} {len(pre[k]):,}" for k in bk), flush=True)
    return libs, groups, inv.reshape(-1).astype(np.int32), P_, Q_, pts, float(ci.get("bead_radius_nm", 0.5))


def step_index(cfg, bdir, kw):
    import pandas as pd
    p, g, reg = _setup(cfg, bdir, kw)
    outs = [p.path(f"preplaced_{n}.npy") for n in ("types", "pos", "quat", "off")] + [p.path("preplaced_libs.json"),
                                                                                       p.path("index_drops.tsv")]
    if p.pc.get("preplaced_index"):                       # grouped npz source (yeast prepass.npz)
        outs += [p.path("preplaced_groups.json")] + ([p.path("mrna_pts.npy"), p.path("mrna_off.npy"), p.path("mrna_meta.json")]
                                                    if p.pc["preplaced_index"].get("beads") else [])
        if _done(*outs):
            return
        meta = json.load(open(p.path("tiles_meta.json"))); lo0 = np.array(meta["lo0"]); TILE = p.TILE
        M = pd.read_parquet(p.path("tiles.parquet"))
        libs, groups, T_, P_, Q_, pts, brad = _index_grouped(p, cfg, bdir, M, lo0)
        T_, P_, Q_ = _index_drops(p, libs, T_, P_, Q_, g)
        tid = _tile_ids(np.floor((P_ - lo0) / TILE).astype(np.int64), M)
        keep = tid >= 0
        _index_outside(p, libs, T_, keep)
        o = np.argsort(tid[keep], kind="stable")
        for name, arr in (("types", T_[keep][o]), ("pos", P_[keep][o]), ("quat", Q_[keep][o])):
            _npy_save(p.path(f"preplaced_{name}.npy"), arr)
        _npy_save(p.path("preplaced_off.npy"), np.searchsorted(tid[keep][o], np.arange(len(M) + 1)))
        _json_save(p.path("preplaced_libs.json"), libs); _json_save(p.path("preplaced_groups.json"), groups)
        print(f"pre-pass instances {len(T_):,} ({len(libs)} (group, library) entries), indexed {int(keep.sum()):,} into "
              f"{len(M):,} tiles ({int((~keep).sum()):,} outside every tile)", flush=True)
        if pts is not None:
            tid = _tile_ids(np.floor((pts - lo0) / TILE).astype(np.int64), M)
            k = tid >= 0; o = np.argsort(tid[k], kind="stable")
            _npy_save(p.path("mrna_pts.npy"), pts[k][o])
            _npy_save(p.path("mrna_off.npy"), np.searchsorted(tid[k][o], np.arange(len(M) + 1)))
            _json_save(p.path("mrna_meta.json"), dict(bead_radius=brad))
            print(f"beads indexed: {int(k.sum()):,} of {len(pts):,}", flush=True)
        return
    man = p.preplaced()
    if man.get("beads"):
        outs += [p.path("mrna_pts.npy"), p.path("mrna_off.npy"), p.path("mrna_meta.json")]
    if _done(*outs):
        return
    meta = json.load(open(p.path("tiles_meta.json"))); lo0 = np.array(meta["lo0"]); TILE = p.TILE
    M = pd.read_parquet(p.path("tiles.parquet"))
    key = {tuple(t): i for i, t in enumerate(M[["ti", "tj", "tk"]].values)}
    libs, T_, P_, Q_ = [], [], [], []
    drops = man.get("drops", [])
    for m in man.get("instances", []):
        z = np.load(m["npz"])
        ids = list(z["lib_ids"]); base = len(libs); libs += ids
        keep_f = np.ones(len(z["types"]), bool)
        for dj in drops:                                 # instances displaced by a later pre-pass (lamina over chromatin)
            if dj["target"] != m["name"]:
                continue
            if int(dj["n"]) == len(keep_f):
                keep_f[np.load(dj["npy"])] = False
                print(f"   {m['name']}: {int((~keep_f).sum()):,} instances dropped ({dj.get('why', '')})", flush=True)
            else:
                print(f"   WARNING {m['name']}: drop list {dj['npy']} was made for {dj['n']} instances, file has {len(keep_f)}; "
                      "not applied (re-run its builder)", flush=True)
        T_.append(z["types"][keep_f].astype(np.int32) + base); P_.append(z["pos"][keep_f].astype(np.float32))
        Q_.append(z["quat"][keep_f].astype(np.float32))
        print("  ", m["name"], m["npz"], len(z["types"]), flush=True)
    if T_:
        T_ = np.concatenate(T_); P_ = np.concatenate(P_); Q_ = np.concatenate(Q_)
    else:                                                 # no pre-pass modules
        T_ = np.zeros(0, np.int32); P_ = np.zeros((0, 3), np.float32); Q_ = np.zeros((0, 4), np.float32)
    T_, P_, Q_ = _index_drops(p, libs, T_, P_, Q_, g)
    tij = np.floor((P_ - lo0) / TILE).astype(np.int64)
    tid = np.array([key.get(tuple(t), -1) for t in tij]) if len(tij) < 5_000_000 else _tile_ids(tij, M)
    keep = tid >= 0
    _index_outside(p, libs, T_, keep)
    o = np.argsort(tid[keep], kind="stable")
    for name, arr in (("types", T_[keep][o]), ("pos", P_[keep][o]), ("quat", Q_[keep][o])):
        _npy_save(p.path(f"preplaced_{name}.npy"), arr)
    off = np.searchsorted(tid[keep][o], np.arange(len(M) + 1))
    _npy_save(p.path("preplaced_off.npy"), off); _json_save(p.path("preplaced_libs.json"), libs)
    print(f"pre-pass instances {len(T_):,}, indexed {int(keep.sum()):,} into {len(M):,} tiles", flush=True)
    beads = man.get("beads", [])
    if beads:                                             # RNA beads (obstacles; beads inside ribosome channels skipped)
        b0 = beads[0]
        z = np.load(b0["npz"]); pts = z["pts"]; keepb = z["inrib"] == 0 if "inrib" in z.files else np.ones(len(pts), bool)
        if b0.get("drop_npy"):
            dr = np.load(b0["drop_npy"]); keepb[dr] = False; print(f"RNA beads dropped (clash list): {len(dr):,}", flush=True)
        pts = pts[keepb]
        for bx in beads[1:]:                              # further bead sets (same format and bead radius)
            zx = np.load(bx["npz"]); px = zx["pts"]; kx = zx["inrib"] == 0 if "inrib" in zx.files else np.ones(len(px), bool)
            assert abs(float(zx["bead_radius"]) - float(z["bead_radius"])) < 1e-6, "bead radii differ"
            if bx.get("drop_npy"):
                kx[np.load(bx["drop_npy"])] = False
            pts = np.concatenate([pts, px[kx]]); print(f"extra RNA beads {bx['name']}: {int(kx.sum()):,}", flush=True)
        tid = _tile_ids(np.floor((pts - lo0) / TILE).astype(np.int64), M)
        k = tid >= 0; o = np.argsort(tid[k], kind="stable")
        _npy_save(p.path("mrna_pts.npy"), pts[k][o].astype(np.float32))
        _npy_save(p.path("mrna_off.npy"), np.searchsorted(tid[k][o], np.arange(len(M) + 1)))
        _json_save(p.path("mrna_meta.json"), dict(bead_radius=float(z["bead_radius"])))
        print(f"RNA beads indexed: {int(k.sum()):,}", flush=True)


# ------------------------------------------------------------------ packing
def _tree(i, LIB, shapes=None):
    """collision sphere tree (centres, radii) of species i: the library's own tree, or its packing shape (pack.shapes)"""
    if shapes is not None and i in shapes:
        return shapes[i]
    z = np.load(LIB[i])
    return z["st_c"].astype(np.float64), z["st_r"].astype(np.float64)


def _libpack(ids, LIB, shapes=None):
    if not len(ids):
        return np.zeros(1, np.int64), np.zeros((0, 3)), np.zeros(0), np.zeros(0)
    so = [0]; C, R, B = [], [], []
    for i in ids:
        c, r = _tree(i, LIB, shapes)
        C.append(c); R.append(r); so.append(so[-1] + len(r)); B.append(float((np.linalg.norm(c, axis=1) + r).max()))
    return np.array(so, np.int64), np.concatenate(C), np.concatenate(R), np.array(B)


def _sig_ok(sig, path, ref_md5=None):
    """does a tree signature match the library file `path`? md5 of the file (core library pack_shapes; `ref_md5`: the
    manifest's md5 of that file) or 'size:mtime' (the yeast engine's wc3d.pack_trees cache, made from engine/lib)"""
    sig = str(sig)
    if ":" in sig:
        st = os.stat(path)
        return sig == f"{st.st_size}:{st.st_mtime:.0f}"
    if ref_md5 is None:
        from wholecell.library.pack_shapes import _md5
        ref_md5 = _md5(path)
    return sig == ref_md5


def _shapes_library(p, sw):
    """({species: npz path} the packer reads, pack.shapes.library dir or None): the manifest, overridden / extended by
    the npz files of pack.shapes.library (yeast: engine/lib) when pack.shapes is on and that key is set"""
    LIB = dict(p.library())
    d = sw.shapes.get("library") if sw.shapes["source"] != "library" else None
    if d:
        d = d if os.path.isabs(d) else os.path.join(p.bdir, d)
        if not os.path.isdir(d):
            raise FileNotFoundError(f"pack.shapes.library: {d} is not a directory")
        LIB.update({n[:-4]: os.path.join(d, n) for n in os.listdir(d) if n.endswith(".npz")})
    return LIB, d


def _shapes(p, sw, ids):
    """pack.shapes (P8): {species: (centres, radii)} from the packing-shape trees for `ids`, or None (library trees).
    A species without a tree must have no beads (nothing to make a tree from: its library tree is used); a tree made
    from another version of the library file than the one the packer sees (sig != manifest md5, or != the file in
    pack.shapes.library) is an error. pack.shapes.library: directory whose npz files override / extend the manifest
    (yeast: engine/lib, the library after assembly and pre-placed models; its trees carry the engine's size:mtime)."""
    if sw.shapes["source"] == "library":
        return None
    import glob
    import pandas as pd
    f = sw.shapes.get("file")
    if f:
        f = f if os.path.isabs(f) else os.path.join(p.bdir, f)
    else:
        cand = sorted(glob.glob(os.path.join(p.bdir, "library", "pack_shapes", "pack_trees_*.npz")))
        if len(cand) != 1:
            raise FileNotFoundError(f"pack.shapes: expected one library/pack_shapes/pack_trees_*.npz in {p.bdir}, found "
                                    f"{len(cand)} (run `library --step pack_shapes` or set pack.shapes.file)")
        f = cand[0]
    z = np.load(f, allow_pickle=True)
    pos = {str(i): k for k, i in enumerate(z["ids"])}
    off, C, R, sig = z["off"], z["C"], z["R"], z["sig"]
    man = pd.read_csv(os.path.join(p.bdir, "library", "library_manifest.tsv"), sep="\t")
    md5 = dict(zip(man.species, man.md5)) if "md5" in man.columns else {}
    LIB, d = _shapes_library(p, sw)
    if d:                                  # the files the trees were made from (and the packer reads): the manifest md5
        md5 = {}                           # does not describe them
    out, nolib, stale, missing = {}, [], [], []
    for i in dict.fromkeys(ids):
        k = pos.get(i)
        if k is None:
            if i not in LIB:
                missing.append(i); continue
            with np.load(LIB[i]) as zz:
                (missing if ("beads" in zz.files and len(zz["beads"])) else nolib).append(i)
            continue
        if i in LIB:
            if (d or i in md5) and not _sig_ok(sig[k], LIB[i], md5.get(i)):
                stale.append(i)
        elif d:                            # a tree of a species that is not in the library the packer sees
            stale.append(i)
        out[i] = (np.asarray(C[off[k]:off[k + 1]], np.float64), np.asarray(R[off[k]:off[k + 1]], np.float64))
    if missing or stale:
        raise ValueError(f"pack.shapes {f}: {len(missing)} species with beads have no tree (e.g. {missing[:3]}), "
                         f"{len(stale)} trees were made from another library file (e.g. {stale[:3]}"
                         + (f"; checked against {d}" if d else "; checked against the library manifest md5")
                         + "); rerun library --step pack_shapes (or the step that wrote the trees file)")
    print(f"packing shapes {os.path.basename(f)}: {len(out)} species, {sum(len(v[1]) for v in out.values()):,} spheres"
          + (f" (signatures checked against {d})" if d else "")
          + (f"; {len(nolib)} species without beads keep their library tree" if nolib else ""), flush=True)
    return out


def _hmin(so, STC, STR, ndir=200, seed=5):
    """per library: minimal support half-width min_u max_b (u . c_b + r_b) over ndir directions (pack3d seed pruning)"""
    u = np.random.default_rng(seed).normal(size=(ndir, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    out = np.zeros(len(so) - 1)
    for i in range(len(so) - 1):
        c = STC[so[i]:so[i + 1]]; r = STR[so[i]:so[i + 1]]
        out[i] = float((c @ u.T + r[:, None]).max(0).min()) if len(r) else 0.0
    return out


def _neigh(t):
    out = []
    for d in np.stack(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"), -1).reshape(-1, 3):
        if d.any():
            j = _S["key"].get(tuple(t + d))
            if j is not None: out.append(j)
    return out


def _tile_npz(ti):
    """packed state of tile ti visible to its neighbours: the finished tile, else its rounds file (pack.rounds n > 1)"""
    p = _S["p"]
    f = p.path("tiles", f"{ti}.npz")
    if os.path.exists(f):
        return f
    if _S["n_rounds"] > 1:
        f = p.path("rounds", f"{ti}.npz")
        if os.path.exists(f):
            return f
    return None


def _pack_tile(ti):
    from skimage.measure import marching_cubes
    from . import pack3d as p3
    p, reg, cell, c3 = _S["p"], _S["reg"], _S["cell"], _S["c3"]
    TILE, VOX, MESH_H, HALO = p.TILE, p.VOX, p.MESH_H, p.HALO
    VOL_REG, MEM_REG, NV, NSH0 = reg.VOL_REG, reg.MEM_REG, reg.NV, reg.NSH0
    t = _S["Tij"][ti]; lo = t * TILE + _S["lo0"]
    glo, ghi = lo - HALO, lo + TILE + HALO
    pilot = _S.get("mode") == "pilot"
    rnd, n_rounds = _S["rnd"], _S["n_rounds"]
    last = pilot or rnd == n_rounds - 1
    out_f = p.path("pilot", f"{ti}.npz") if pilot else p.path("tiles", f"{ti}.npz")
    if os.path.exists(out_f): return ti, -1, -1, False, 0.0
    rnd_f = p.path("rounds", f"{ti}.npz")
    own = None
    if not pilot and n_rounds > 1:
        if os.path.exists(rnd_f):
            own = dict(np.load(rnd_f))
            if int(own["round"]) >= rnd:                 # this round is already done (resumed run)
                return ti, -1, -1, False, 0.0
            if not reg.identity:                         # files hold config region ids; the kernel its own
                own["job_region"] = reg.k_of[own["job_region"].astype(np.int64)]
    t0 = time.time()
    # jobs: membranes first, then shells (+ virtual volumes), then volumes; largest bounding radius first
    a, b = _S["coff"][ti], _S["coff"][ti + 1]
    jl, jr, jn = _S["clib"][a:b].astype(np.int64), _S["creg"][a:b].astype(np.int64), _S["cn"][a:b].astype(np.int64)
    if len(jl):
        cls = reg.pack_class[jr]
        o = np.lexsort((-_S["BR"][jl], cls)); jl, jr, jn = jl[o], jr[o], jn[o]
    if n_rounds > 1 and not pilot:                       # pack.rounds: this round's size class only
        from .algorithms import job_rounds
        k = job_rounds(_S["BR"][jl], reg.pack_class[jr] == 0, _S["edges"]) == rnd
        jl, jr, jn = jl[k], jr[k], jn[k]
        if not len(jl) and not last:
            return ti, 0, 0, False, 0.0                   # nothing to add in this round
    # seeds per volume region (vox_nm voxels); evaluator: the whole cell, or the prims / pores near the tile
    # (pack.tiling.box_margin_nm, as the tile measures; yeast cell.args(lo, hi, margin))
    from .tiling import Evaluator, region_voxels
    ev = Evaluator(cell, c3, lo, TILE, _S["box"][0], _S["box"][1])
    kargs = ev.args if ev.args is not None else (cell.G, cell.pores, cell.pidx, cell.PT, cell.PQ, cell.PC)
    n = int(TILE / VOX); g = (np.arange(n) + 0.5) * VOX
    X, Y, Z = np.meshgrid(lo[0] + g, lo[1] + g, lo[2] + g, indexing="ij")
    P = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1); lab = ev.labels(P)
    vox_off, VP = [0], []
    vmask = region_voxels(reg, ev, P, lab, _S["vset"])
    thin_r = reg.regmem()
    gp = cell.params
    for ri_, (name, labs) in enumerate(VOL_REG):
        m = vmask.get(name)
        if m is None: m = np.zeros(len(P), bool)
        Pr = P[m]
        if p.SEED_PRUNE:                                 # membrane clearance per seed voxel, thin lumens sorted by it
            c_ = np.full(len(Pr), 1e3)
            if thin_r[ri_] >= 0 and len(Pr):
                Sm = cell.mems(Pr); c_ = np.min(np.abs(Sm) - np.asarray(gp["MEM_HALF"])[None, :], axis=1) + gp["HALF"]
                o_ = np.argsort(-c_, kind="stable"); Pr, c_ = Pr[o_], c_[o_]
            Pr = np.c_[Pr, c_]
        VP.append(Pr); vox_off.append(vox_off[-1] + int(m.sum()))
    vox_off += [vox_off[-1]] * len(MEM_REG); VOXP = np.concatenate(VP); del P, X, Y, Z, lab
    # membrane triangles
    m_ = int(TILE / MESH_H) + 1; gm = np.arange(m_) * MESH_H
    X, Y, Z = np.meshgrid(lo[0] + gm, lo[1] + gm, lo[2] + gm, indexing="ij")
    S = ev.mems(np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1)).reshape(m_, m_, m_, -1)
    TRI, tri_off = [], [0]
    for name, mi in MEM_REG:
        s = S[..., mi]
        if s.min() < 0 < s.max():
            v, fc, _, _ = marching_cubes(s, 0.0, spacing=(MESH_H,) * 3); tt = face_keep((v + lo)[fc].reshape(-1, 9), cell.G, c3)
        else:
            tt = np.zeros((0, 9))
        TRI.append(tt); tri_off.append(tri_off[-1] + len(tt))
    TRI = np.concatenate(TRI); del S
    ar = 0.5 * np.linalg.norm(np.cross(TRI[:, 3:6] - TRI[:, 0:3], TRI[:, 6:9] - TRI[:, 0:3]), axis=1)
    CUMA = np.cumsum(ar) if len(ar) else np.zeros(1)
    # obstacles: pre-pass instances in this tile + neighbours, RNA beads, and packed neighbour tiles (pilot: none)
    OCs, ORs = [], []
    org = _S["org"]
    ATs, APs = [], []                                    # pack.organisation: placed instances around the tile (anchors)
    nb_ = _neigh(t)
    for j in [ti] + nb_:
        p0, p1 = _S["poff"][j], _S["poff"][j + 1]
        if p1 > p0 and org is not None:
            pt_ = org[2][_S["ptypes"][p0:p1]]; pp_ = _S["ppos"][p0:p1]
            k_ = (pt_ >= 0) & np.all((pp_ >= glo) & (pp_ <= ghi), axis=1)
            ATs.append(pt_[k_]); APs.append(pp_[k_])
        if p1 > p0:
            C0, R0 = p3.emit_spheres(_S["ptypes"][p0:p1].astype(np.int64), _S["ppos"][p0:p1].astype(np.float64), _S["pquat"][p0:p1].astype(np.float64),
                                     _S["pso"], _S["pSTC"], _S["pSTR"], _S["pBR"], glo, ghi)
            OCs.append(C0); ORs.append(R0)
        if "moff" in _S:
            m0, m1 = _S["moff"][j], _S["moff"][j + 1]
            if m1 > m0:
                mp_ = np.asarray(_S["mpts"][m0:m1], np.float64)
                kk = np.all((mp_ >= glo - 1) & (mp_ <= ghi + 1), axis=1)
                OCs.append(mp_[kk]); ORs.append(np.full(int(kk.sum()), _S["mrad"]))
        f = _tile_npz(j) if (j != ti and not pilot) else None
        if f is not None:
            z = np.load(f)
            C1, R1 = p3.emit_spheres(z["types"].astype(np.int64), z["pos"].astype(np.float64), z["quat"].astype(np.float64),
                                     _S["so"], _S["STC"], _S["STR"], _S["BR"], glo, ghi)
            OCs.append(C1); ORs.append(R1)
            if org is not None:
                k_ = np.all((z["pos"] >= glo) & (z["pos"] <= ghi), axis=1); ATs.append(z["types"][k_]); APs.append(z["pos"][k_])
    OCx = np.concatenate(OCs) if OCs else np.zeros((0, 3)); ORx = np.concatenate(ORs) if ORs else np.zeros(0)
    if own is not None and len(own["types"]):          # earlier rounds of this tile
        C1, R1 = p3.emit_spheres(own["types"].astype(np.int64), own["pos"].astype(np.float64), own["quat"].astype(np.float64),
                                 _S["so"], _S["STC"], _S["STR"], _S["BR"], glo, ghi)
        OC, OR = np.concatenate([OCx, C1]), np.concatenate([ORx, R1])
        if org is not None:
            ATs.append(own["types"]); APs.append(own["pos"])
    else:
        OC, OR = OCx, ORx
    nsph = int((np.diff(_S["so"])[jl] * jn).sum()) if len(jl) else 0
    cap_s = nsph + len(OR) + 1000; cap_e = cap_s * 14
    e = np.zeros(0, np.int64)
    EXC = _exclusion(ti, nb_, glo, ghi)
    if pilot:
        seed = _S["seed_fn"]("pack.thinning.pilot", int(ti))
    elif rnd == 0:
        seed = _S["tile_seed"](int(ti))
    else:
        seed = _S["seed_fn"]("pack.round", int(ti), int(rnd))
    if org is not None:                                  # partner attraction (kernel arguments AT ... prm[0:7])
        AT = np.concatenate(ATs).astype(np.int64) if ATs else e
        AP = np.concatenate(APs).astype(np.float64) if APs else np.zeros((0, 3))
        pr = org[0]; ptn = (AT, AP, pr["off"], pr["idx"], pr["strength"], pr["kind"], np.r_[org[1], p.THIN_TRY_F, p.THIN_GIVEUP_F])
    else:
        ptn = (e, np.zeros((0, 3)), e, e, np.zeros(0), e, np.r_[np.zeros(7), p.THIN_TRY_F, p.THIN_GIVEUP_F])
    res = p3.pack_tile(seed, lo, TILE, HALO, *kargs,
                       _S["so"], _S["STC"], _S["STR"], _S["BR"], reg.regmask(), jl, jr, jn, NV, np.array(vox_off, np.int64), VOXP, VOX,
                       np.array(tri_off, np.int64), TRI, CUMA, OC, OR, 0, p.MRNA_MEAN, 0, p.MAX_TRIES, p.GIVE_UP, cap_s, cap_e,
                       *ptn, _S["QF"],
                       reg.regmem(), *reg.shell_arrays(), _S["REG_TRIES"], _S["REG_GIVEUP"], EXC,
                       _S["EXC_TYPE"] if len(EXC) else e, *_S["band"], _S["mcode"])
    types, pos, quat, placed, _, _, overflow, ns, ne, att = res
    n_req, n_new = int(jn.sum()), int(placed.sum())
    extra = {}
    if org is not None:                                  # attraction proposals / copies placed next to a partner (all rounds)
        extra["att_stats"] = att.astype(np.int64) + (own["att_stats"] if own is not None and "att_stats" in own else 0)
    if own is not None:                                  # instances stay grouped by job: earlier rounds' jobs, then this round's
        types = np.concatenate([own["types"], types]); pos = np.concatenate([own["pos"], pos]); quat = np.concatenate([own["quat"], quat])
        jl = np.concatenate([own["job_lib"].astype(np.int64), jl]); jr = np.concatenate([own["job_region"].astype(np.int64), jr])
        jn = np.concatenate([own["job_n"].astype(np.int64), jn]); placed = np.concatenate([own["job_placed"].astype(np.int64), placed])
        overflow = bool(overflow) or bool(own["overflow"])
    sw = _S["sw"]
    if last and (not pilot or sw.thinning["pilot_squeeze"]) and sw.squeeze["scope"] != "none" and len(jl) and not overflow:
        from .algorithms import squeeze_jobs
        sel = squeeze_jobs(sw.squeeze["scope"], float(sw.squeeze["threshold"]), jl, jr, jn, placed, _S["mass"], reg.pack_class == 2)
        if sel.any():
            from .squeeze import squeeze_tile
            before = int(placed.sum())
            types, pos, quat, placed = squeeze_tile(_S["seed_rng"]("pack.thinning.squeeze" if pilot else "pack.squeeze", int(ti)),
                                                    lo, TILE, HALO, cell,
                                                    _S["so"], _S["STC"], _S["STR"], _S["BR"], reg.regmask(), jl, jr, jn, placed,
                                                    types, pos, quat, np.array(vox_off, np.int64), VOXP, VOX, OCx, ORx, sel,
                                                    int(sw.squeeze["rounds"]), int(sw.squeeze["iterations"]),
                                                    float(sw.squeeze["insert_frac"]), EXC, _S["EXC_TYPE"], args=kargs,
                                                    band=_S["band"])
            n_new += int(placed.sum()) - before
    if not reg.identity:                                 # tile files hold config region ids
        jr = reg.c_of[jr]
    if last:
        _atomic(out_f, lambda tmp: np.savez(tmp, types=types.astype(np.int32), pos=pos, quat=quat, job_lib=jl.astype(np.int32),
                                            job_region=jr.astype(np.int16), job_n=jn.astype(np.int32),
                                            job_placed=placed.astype(np.int32), overflow=overflow, **extra))
        if own is not None:
            os.remove(rnd_f)
    else:
        _atomic(rnd_f, lambda tmp: np.savez(tmp, types=types.astype(np.int32), pos=pos, quat=quat, job_lib=jl.astype(np.int32),
                                            job_region=jr.astype(np.int16), job_n=jn.astype(np.int32),
                                            job_placed=placed.astype(np.int32), overflow=overflow, round=rnd, **extra))
    return ti, n_req, n_new, bool(overflow), time.time() - t0


def _exclusion(ti, nb_, glo, ghi):
    """exclusion spheres (pack.exclusion) reaching this tile + halo, from the pre-pass index of the tile and neighbours"""
    if not _S["exc_rules"]:
        return np.zeros((0, 5))
    from .algorithms import exclusion_spheres, exclusion_tile_spheres
    out = []
    if len(_S["exc_zones"]):                             # provider zones (geometry features)
        E = exclusion_tile_spheres(_S["exc_zones"], glo, ghi)
        if len(E):
            out.append(E)
    for j in [ti] + nb_:
        p0, p1 = _S["poff"][j], _S["poff"][j + 1]
        if p1 > p0:
            E = exclusion_spheres(_S["ptypes"][p0:p1], _S["ppos"][p0:p1], _S["exc_around"], _S["exc_rad"], glo, ghi)
            if len(E):
                out.append(E)
    return np.ascontiguousarray(np.concatenate(out)) if out else np.zeros((0, 5))


def _tile_selection(p, M, meta):
    kw = p.kw
    sel = np.arange(len(M))
    if kw.get("box"):
        b = np.array([float(x) for x in kw["box"].split(",")])
        c = M[["ti", "tj", "tk"]].values * p.TILE + np.array(meta["lo0"]) + p.TILE / 2
        sel = np.where(np.all((c >= b[:3]) & (c <= b[3:]), axis=1))[0]; print("sub-box tiles:", len(sel), flush=True)
    if kw.get("tiles"):
        t_ = kw["tiles"]
        sel = np.load(t_) if t_.endswith(".npy") else np.array([int(x) for x in t_.split(",")]); print("tile list:", len(sel), flush=True)
    return sel


def _run_tiles(p, tiles, log, label):
    """pack `tiles` (one colour: independent) with p.workers forked workers; returns (instances, placed, overflow)"""
    import multiprocessing as mp
    t0 = time.time(); done = 0; inst = plc = 0; ovf = 0
    with mp.get_context("fork").Pool(p.workers) as pool:
        for ti, ni, np_, ov, dt in pool.imap_unordered(_pack_tile, tiles, chunksize=1):
            done += 1
            if ni >= 0:
                inst += ni; plc += np_; ovf += ov
                if log is not None:
                    log.write(f"{label}\t{ti}\t{ni}\t{np_}\t{int(ov)}\t{dt:.1f}\n")
            if done % 200 == 0:
                if log is not None:
                    log.flush()
                print(f"  phase {label}: {done}/{len(tiles)} tiles, placed {plc:,}/{inst:,}, {time.time() - t0:.0f}s", flush=True)
    if log is not None:
        log.flush()
    return inst, plc, ovf, time.time() - t0


def _thinning(p, cfg, sw, reg, M):
    """pack.thinning (P2): pilot-pack crowded volume regions, thin their request; writes thinning.json and
    counts_thinned.npz (counts.npz is kept as the unthinned request) and returns the thinned count arrays"""
    import shutil
    from wholecell import seeds
    from .algorithms import crowded_regions, pilot_tiles, thinning_factor, thin_counts, sphere_volumes
    th = sw.thinning
    fj, fc = p.path("thinning.json"), p.path("counts_thinned.npz")
    if os.path.exists(fj) and os.path.exists(fc):
        z = np.load(fc)
        print(f"thinning: {fj} exists, using {fc}", flush=True)
        return z["tile"], z["lib"], _kreg(reg, z["region"]), z["n"]
    clib, creg, cn, ctile = _S["clib"], _S["creg"], _S["cn"], _S["ctile"]          # kernel region ids
    vol_lib = sphere_volumes(_S["so"], _S["STR"])
    K = reg.KNAMES
    rvol = np.array([float(M[r].sum()) if r in M.columns else 0.0 for r in K])
    ref = K.index(th["reference_region"])
    crowd, phi = crowded_regions(clib.astype(np.int64), creg.astype(np.int64), cn.astype(np.int64), vol_lib, rvol,
                                 range(reg.NV), ref, float(th["min_rel"]))
    info = dict(reference_region=th["reference_region"], requested_phi={K[r]: round(f, 4) for r, f in phi.items()},
                tested={}, factors={}, settings=th)
    factors = {}
    if crowd:
        rng = seeds.rng(cfg, "pack.thinning")
        chosen = pilot_tiles(sorted(crowd), {r: M[K[r]].values for r in crowd}, _S["Tij"], int(th["pilot_tiles"]), rng)
        tl = sorted({t for v in chosen.values() for t in v})
        shutil.rmtree(p.path("pilot"), ignore_errors=True); os.makedirs(p.path("pilot"))
        _S.update(mode="pilot", rnd=0, coff=np.searchsorted(ctile, np.arange(len(M) + 1)))
        try:
            inst, plc, _, dt = _run_tiles(p, tl, None, "pilot")
            print(f"thinning pilot: {len(tl)} tiles, placed {plc:,}/{inst:,}, {dt:.0f}s", flush=True)
            for r in sorted(crowd):
                rq = pl = 0.0
                for ti in chosen[r]:
                    z = np.load(p.path("pilot", f"{ti}.npz")); k = z["job_region"] == reg.c_of[r]     # files: config ids
                    rq += float((z["job_n"][k] * vol_lib[z["job_lib"][k]]).sum()); pl += float((z["job_placed"][k] * vol_lib[z["job_lib"][k]]).sum())
                held = pl / max(rq, 1e-9)
                f = thinning_factor(held, float(th["safety"]), float(th["keep_above"]))
                info["tested"][K[r]] = dict(phi=round(crowd[r], 4), tiles=chosen[r], held=round(held, 4), factor=round(f, 4))
                print(f"  pilot {K[r]}: requested sphere volume fraction {crowd[r]:.3f}; {len(chosen[r])} tiles held "
                      f"{held:.3f} of the requested volume -> every species x {f:.3f}", flush=True)
                if f < 1.0:
                    factors[r] = f
        finally:
            _S.update(mode=None)
            shutil.rmtree(p.path("pilot"), ignore_errors=True)
    method = (p.pc.get("requests") or {}).get("method", "multinomial")
    if factors and th["resplit"] and method == "legacy_x5":     # yeast: re-split from the thinned totals (same seed)
        from . import requests
        d, _ = requests.legacy_x5(cfg, p.bdir, p, reg, M, json.load(open(p.path("tiles_meta.json")))["lo0"], sw=sw,
                                  factors={K[r]: float(f) for r, f in factors.items()}, log=lambda *a: print(*a, flush=True))
        tile, lib, region, n = (d.tile.values.astype(np.int64), d.lib.values.astype(np.int64),
                                _kreg(reg, d.region.values), d.n.values.astype(np.int64))
        info["resplit"] = True
    else:
        tile, lib, region, n = thin_counts(ctile, clib, creg, cn, factors, seeds.rng(cfg, "pack.thinning.draw"))
    info["factors"] = {K[r]: round(f, 4) for r, f in factors.items()}
    info["factors_exact"] = {K[r]: float(f) for r, f in factors.items()}
    info["instances"] = dict(before=int(cn.sum()), after=int(n.sum()))
    _atomic(fc, lambda tmp: np.savez_compressed(tmp, tile=tile.astype(np.int32), lib=lib.astype(np.int32),
                                                region=_creg(reg, region).astype(np.int16), n=n.astype(np.int32)))
    _json_save(fj, info, indent=1)
    print(f"thinning: {int(cn.sum()):,} -> {int(n.sum()):,} instances requested; factors {info['factors']}", flush=True)
    return tile, lib, region, n


def _kreg(reg, r):
    """config region ids -> kernel ids (identity tables: unchanged)"""
    return r if reg.identity else reg.k_of[np.asarray(r).astype(np.int64)]


def _creg(reg, r):
    """kernel region ids -> config ids"""
    return r if reg.identity else reg.c_of[np.asarray(r).astype(np.int64)]


def _exclusion_zones(cfg, bdir, g, rules):
    """(m, 5) zone spheres (x, y, z, r, rule) of the provider rules of pack.exclusion (registry kind exclusion_zones)"""
    from wholecell import registry
    out = []
    for i, e in enumerate(rules):
        a = e["around"]
        if not a.get("provider"):
            continue
        Z = np.asarray(registry.get("exclusion_zones", a["provider"])(cfg, bdir, g, dict(a.get("params") or {})), np.float64)
        Z = Z.reshape(-1, 4)
        print(f"exclusion rule {e.get('name', i)}: {len(Z)} zone spheres from provider {a['provider']}", flush=True)
        out.append(np.c_[Z, np.full(len(Z), float(i))])
    return np.ascontiguousarray(np.concatenate(out)) if out else np.zeros((0, 5))


def step_pack(cfg, bdir, kw):
    import pandas as pd
    from wholecell import seeds
    from . import switches
    from .algorithms import round_edges, phase_order, budget_arrays, exclusion_tables
    p, g, reg = _setup(cfg, bdir, kw)
    reg.require_core_packer("the core pack step")
    sw = switches.load(p.pc, reg.RNAMES)
    if not reg.identity:
        print(f"region table in config order ({len(reg.RNAMES)} ids); kernel order: {reg.KNAMES}", flush=True)
    M = pd.read_parquet(p.path("tiles.parquet")); meta = json.load(open(p.path("tiles_meta.json")))
    cell = g.cell
    LIB, _ = _shapes_library(p, sw)                       # manifest (+ pack.shapes.library files)
    LJ = json.load(open(p.path("libs.json")))
    L = LJ["libs"]
    pl = json.load(open(p.path("preplaced_libs.json")))
    shapes = _shapes(p, sw, list(L) + list(pl))           # pack.shapes (P8): None = library trees (v7)
    so, STC, STR, BR = _libpack(L, LIB, shapes)
    pso, pSTC, pSTR, pBR = _libpack(pl, LIB, shapes)
    cz = np.load(p.path("counts.npz"))
    Tij = M[["ti", "tj", "tk"]].values
    from .regions import qflat
    if not p.LUMEN_ALIGN:
        QF = np.zeros((0, 4))
    else:
        QF = qflat(so, STC, STR)
        if p.SEED_PRUNE:
            QF = np.c_[QF, _hmin(so, STC, STR, p.HMIN_NDIR, p.HMIN_SEED)]
    REG_TRIES, REG_GIVEUP = budget_arrays(sw.budgets["per_region"], reg.KNAMES)
    exc_type, exc_around, exc_rad = exclusion_tables(sw.exclusion, L, pl) if sw.exclusion else (np.zeros(0, np.int64), None, None)
    if sw.exclusion:
        for i, e_ in enumerate(sw.exclusion):
            if not e_["around"].get("provider") and float(e_["r_nm"]) > p.TILE:
                raise ValueError(f"config pack.exclusion[{i}].r_nm = {e_['r_nm']} > tile_nm {p.TILE}: zones are read from the tile and its neighbours only")
            if not e_["around"].get("provider"):
                print(f"exclusion rule {e_.get('name', i)}: {int(((exc_type >> i) & 1).sum())} packed species kept "
                      f"{float(e_['r_nm'])} nm from {int(exc_around[i].sum())} pre-placed species", flush=True)
            else:
                print(f"exclusion rule {e_.get('name', i)}: {int(((exc_type >> i) & 1).sum())} packed species", flush=True)
    _S.update(p=p, reg=reg, cell=cell, c3=g.engine, Tij=Tij, lo0=np.array(meta["lo0"]), key={tuple(t): i for i, t in enumerate(Tij)},
              vset=_virtual(p, reg, g), so=so, STC=STC, STR=STR, BR=BR, pso=pso, pSTC=pSTC, pSTR=pSTR, pBR=pBR, QF=QF,
              clib=cz["lib"], creg=_kreg(reg, cz["region"]), cn=cz["n"], ctile=cz["tile"],
              tile_seed=lambda k: seeds.seed_int(cfg, "pack.tile", k),
              seed_fn=lambda name, *k: seeds.seed_int(cfg, name, *k), seed_rng=lambda name, *k: seeds.rng(cfg, name, *k),
              sw=sw, mode=None, rnd=0, n_rounds=int(sw.rounds["n"]), edges=np.zeros(0),
              mass=np.array(LJ["mass_da"], np.float64) if LJ.get("mass_da") else np.ones(len(L)),     # legacy_x5: by count
              box=((p.pc.get("tiling") or {}).get("box_margin_nm"), (p.pc.get("tiling") or {}).get("pore_reach_nm", 0.0)),
              band=reg.band_arrays(), mcode=reg.membrane_codes(cell.G, g.engine), exc_zones=_exclusion_zones(cfg, bdir, g, sw.exclusion) if sw.exclusion else np.zeros((0, 5)),
              REG_TRIES=REG_TRIES, REG_GIVEUP=REG_GIVEUP, EXC_TYPE=exc_type, exc_rules=bool(sw.exclusion),
              exc_around=exc_around, exc_rad=exc_rad, org=_partners(p, sw, L, pl) if sw.organisation["enabled"] else None,
              # loaded fully into RAM before the fork (shared copy-on-write)
              ptypes=np.load(p.path("preplaced_types.npy")), ppos=np.load(p.path("preplaced_pos.npy")),
              pquat=np.load(p.path("preplaced_quat.npy")), poff=np.load(p.path("preplaced_off.npy")))
    if os.path.exists(p.path("mrna_off.npy")):
        _S.update(mpts=np.load(p.path("mrna_pts.npy")), moff=np.load(p.path("mrna_off.npy")),
                  mrad=json.load(open(p.path("mrna_meta.json")))["bead_radius"])
    mn_ = {int(c): n for n, c in reg.mem_codes.items()}
    print("membrane-protein tests (membrane codes of the package):", {k: f"{int(c)} {mn_.get(int(c), 'none')}" for k, c in
          zip(("outer wall", "exempt from it", "pores"), _S["mcode"])}, flush=True)
    if not sw.legacy:
        print("packer switches (non-default):", {k: v for k, v in sw.as_dict().items() if v != switches.DEFAULTS[k]}, flush=True)
    if sw.thinning["enabled"]:
        tile_, lib_, reg_, n_ = _thinning(p, cfg, sw, reg, M)
        _S.update(clib=lib_, creg=reg_, cn=n_, ctile=tile_)
    _S["coff"] = np.searchsorted(_S["ctile"], np.arange(len(M) + 1))
    n_rounds = _S["n_rounds"]
    if n_rounds > 1:
        _S["edges"] = round_edges(_S["clib"].astype(np.int64), _S["creg"].astype(np.int64), _S["cn"].astype(np.int64),
                                  so, STR, BR, reg.pack_class, n_rounds)
        os.makedirs(p.path("rounds"), exist_ok=True)
        print(f"{n_rounds} packing rounds by size class; bounding-radius thresholds (nm) {np.round(_S['edges'], 1).tolist()}", flush=True)
    os.makedirs(p.path("tiles"), exist_ok=True)
    sel = _tile_selection(p, M, meta)
    phases = [int(x) for x in str(kw["phases"]).split(",")] if kw.get("phases") not in (None, "") else list(range(8))
    colour = (Tij[sel, 0] % 2) * 4 + (Tij[sel, 1] % 2) * 2 + (Tij[sel, 2] % 2)
    log = open(p.path("pack_log.tsv"), "a")
    T0 = time.time()
    for rnd in range(n_rounds):
        _S["rnd"] = rnd
        r_inst = r_plc = 0
        for ph in phase_order(phases, rnd, bool(sw.rounds["alternate"])):
            tiles = sel[colour == ph]
            inst, plc, ovf, dt = _run_tiles(p, tiles, log, ph)
            r_inst += inst; r_plc += plc
            print(f"phase {ph} done: {len(tiles)} tiles, placed {plc:,}/{inst:,} ({plc / max(inst, 1):.1%}), overflow {ovf}, "
                  f"{dt:.0f}s", flush=True)
        if n_rounds > 1:
            print(f"round {rnd} done: placed {r_plc:,}/{r_inst:,}, {time.time() - T0:.0f}s since the start", flush=True)
    log.close()
    if kw.get("tiles") is None and kw.get("box") is None and kw.get("phases") in (None, ""):
        missing = [i for i in range(len(M)) if not os.path.exists(p.path("tiles", f"{i}.npz"))]
        if missing:
            raise RuntimeError(f"{len(missing)} tiles missing after the pack (first {missing[:5]})")
        done = dict(tiles=len(M), finished=time.strftime("%Y-%m-%dT%H:%M:%S"))
        if not sw.legacy:
            done["switches"] = sw.as_dict()
        _json_save(p.path("pack_done.json"), done)


def step_summary(cfg, bdir, kw):
    import pandas as pd
    p, g, reg = _setup(cfg, bdir, kw)
    RNAMES = reg.RNAMES
    M = pd.read_parquet(p.path("tiles.parquet")); L = json.load(open(p.path("libs.json")))
    by_count = not L.get("mass_da")                      # legacy_x5 libs.json: no masses -> instance counts
    mass = np.ones(len(L["libs"])) if by_count else np.array(L["mass_da"])
    exp = {}; plc = {}; ninst = 0; nt = 0; att = None
    for ti in range(len(M)):
        f = p.path("tiles", f"{ti}.npz")
        if not os.path.exists(f): continue
        z = np.load(f); nt += 1; ninst += len(z["types"])
        if "att_stats" in z.files:                       # pack.organisation: partner proposals, copies placed next to a partner
            att = z["att_stats"].astype(np.int64) + (0 if att is None else att)
        for l, r, n, pp in zip(z["job_lib"], z["job_region"], z["job_n"], z["job_placed"]):
            exp[r] = exp.get(r, 0.0) + n * mass[l]; plc[r] = plc.get(r, 0.0) + pp * mass[l]
    tot_e, tot_p = sum(exp.values()), sum(plc.values())
    if by_count:
        out = {RNAMES[r]: dict(requested=int(exp[r]), placed=int(plc[r]), frac=round(plc[r] / max(exp[r], 1), 4)) for r in sorted(exp)}
        S = dict(tiles_packed=nt, instances=ninst, requested=int(tot_e), placed=int(tot_p),
                 placed_frac=round(tot_p / max(tot_e, 1), 4), regions=out, note="no library masses in libs.json: by count")
    else:
        out = {RNAMES[r]: dict(expected_pg=round(exp[r] / 6.022e23 * 1e12, 2), placed_pg=round(plc[r] / 6.022e23 * 1e12, 2),
                               frac=round(plc[r] / max(exp[r], 1), 3)) for r in sorted(exp)}
        S = dict(tiles_packed=nt, instances=ninst, placed_pg=round(tot_p / 6.022e23 * 1e12, 2), expected_pg=round(tot_e / 6.022e23 * 1e12, 2),
                 mass_frac=round(tot_p / max(tot_e, 1), 3), regions=out)
    if att is not None:
        S["organisation"] = dict(copies_with_partner_proposals=int(att[0]), placed_next_to_partner=int(att[1]))
    if L.get("ledger"):                                  # copies not requested (pack/copy_ledger.tsv, counts step)
        S["copy_ledger"] = L["ledger"]
    if os.path.exists(p.path("index_drops.tsv")):        # pre-pass instances dropped at indexing (pack/index_drops.tsv)
        D = pd.read_csv(p.path("index_drops.tsv"), sep="\t", keep_default_na=False)
        S["index_drops"] = {k: dict(instances=int(g.n.sum()), returned_to_pool=int(g.n[g.returned_to != ""].sum()),
                                    by_lib={str(a): int(b) for a, b in zip(g.lib, g.n)}, reason=str(g.reason.iloc[0]))
                            for k, g in D.groupby("name", sort=False)}
    if os.path.exists(p.path("thinning.json")):          # pack.thinning: the request above is the thinned one
        cz = np.load(p.path("counts.npz")); e0 = float((cz["n"] * mass[cz["lib"]]).sum())
        S["thinning"] = json.load(open(p.path("thinning.json")))["factors"]
        if by_count:
            S["requested_unthinned"] = int(e0); S["placed_frac_of_unthinned"] = round(tot_p / max(e0, 1), 4)
        else:
            S["expected_unthinned_pg"] = round(e0 / 6.022e23 * 1e12, 2)
            S["mass_frac_of_unthinned"] = round(tot_p / max(e0, 1), 3)
    _json_save(p.path("summary.json"), S, indent=1)
    print(json.dumps(S, indent=1), flush=True)


STEPS = {"tiles": step_tiles, "organisation": step_organisation, "counts": step_counts, "index": step_index, "pack": step_pack, "summary": step_summary}


def run(cfg, bdir, step=None, **kw):
    from .switches import effective
    eng = effective(cfg.get("pack") or {}, kw.get("preset")).get("engine")
    if eng and eng not in ("v7", "jurkat_v7", "core_v7", "core") and step in (None, "pack", "all"):
        # (the core tiles / counts steps run for every cell: A2)
        # alternative packer selected by config (registry kind pack_engine, e.g. yeast_rsa_v1 from plugins/yeast)
        from .. import registry
        return registry.get("pack_engine", eng)(cfg, bdir, step=step, **kw)
    if step not in STEPS:
        raise SystemExit(f"pack: --step must be one of {list(STEPS)}")
    numba_threading(step, effective(cfg.get("pack") or {}, kw.get("preset")).get("numba_threads", 8))
    STEPS[step](cfg, bdir, kw)


def numba_threading(step, numba_threads=8):
    """numba threading of a pack step: tiles / pack workers are forked and single-threaded (fork-safe workqueue layer),
    counts / index use pack.numba_threads (an exported NUMBA_NUM_THREADS wins).
    numba may already be imported here (`python -m wholecell` loads the cell's plugins first, and they import numba), with
    its config read while NUMBA_NUM_THREADS was unset (= the CPU count, 32). Setting only the environment then leaves
    config 32 / environment 1: threads launched without a recompile keep that state, and the first compile in a forked
    worker reloads the config and raises 'Cannot set NUMBA_NUM_THREADS to a different value once the threads have been
    launched (currently have 1, trying to set 32)' (yeast C3 thinning pilot, build 20261119-40d7967-s0-step2). So the
    numba config is reloaded here, before any thread is launched."""
    os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
    if step in ("tiles", "pack"):
        os.environ["NUMBA_NUM_THREADS"] = "1"
    else:
        os.environ.setdefault("NUMBA_NUM_THREADS", str(int(numba_threads)))
    if "numba" not in sys.modules:
        return
    from numba.core import config as nconfig
    from numba.np.ufunc import parallel
    want = int(os.environ["NUMBA_NUM_THREADS"])
    if parallel._is_initialized:
        if int(nconfig.NUMBA_NUM_THREADS) != want:
            raise RuntimeError(f"pack --step {step}: numba threads were launched with NUMBA_NUM_THREADS="
                               f"{nconfig.NUMBA_NUM_THREADS} before the pack step set {want}; export NUMBA_NUM_THREADS={want}")
        return
    nconfig.reload_config()
