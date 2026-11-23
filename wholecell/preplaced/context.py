# new in wholecell (unification B4, 2026-11-17)
"""Core context of the serial pre-placed stage (`core.cellgeom`; wholecell/preplaced/serial.py) for any cell built on
a cellgeom_v1 package with the core evaluator (Jurkat; every new cell). The yeast X5 layout has its own context
(plugins/yeast/preplaced/context.py, `yeast.s3`); both hand the builders the same fields:

  ctx.cell     geometry: .G, .pores (n, 6: centre on the NE midsurface, outward normal), .npc_pos, .args(lo, hi, margin)
               -> (G, pores, pidx, PT, PQ, PC), .label(P) -> label codes, .geom (the core Geom: layers, sections)
  ctx.K        the evaluator's kernels (wholecell/preplaced/kern.py)
  ctx.lp       library pack of the types the builders place: ids, st_off, STC, STR, BR (the library's sphere trees,
               as the pack uses them: pack.shapes source library)
  ctx.lib_path(id)   library npz of a species
  ctx.exp3     pools: one row per pack-table ingredient (lib_id = ing_id, region "all", expected = copies)
  ctx.seed(s)  integer seed of stream s (a string: wholecell.seeds stream of this build)
  ctx.obs      named obstacle sphere sets; context_params.obstacles adds frozen pre-placed layers (e.g. the
               centrioles and cortical actin of the frozen v7 cytoskeleton, selected by kind)
  ctx.groups   instance groups placed so far (name, types, pos, quat): ctx.pre_instances() and the module exports

context_params:
  libs: [<species>, ...]                  types the builders may place (library pack order = type index)
  obstacles: [{name, npz, kinds: [..], kind_key: kind, obs: <obs name>, shape: tree|bounding}]   frozen instance layers
                                          as obstacles (bounding: one sphere of the library bounding radius per instance)

A builder that places instances at positions fixed by the geometry (core.npc_pores: the package pores) cannot move
them away from a frozen layer; it drops the overlapping frozen instances instead (drop_frozen_overlaps: logged, left
out of the obstacle sets and of the frozen_subset exports of that file, recorded in serial_manifest frozen_obstacles).
"""
from __future__ import annotations

import json
import os

import numpy as np

from wholecell.registry import register


class CoreCell:
    """the builders' view of a core geometry (cellgeom_v1 package, core evaluator)"""

    def __init__(self, geom):
        c = geom.cell
        self.geom = geom
        self.G, self.pores = c.G, np.asarray(c.pores, np.float64)
        self.pidx, self.PT, self.PQ, self.PC = c.pidx, c.PT, c.PQ, c.PC
        self.npc_pos = self.pores[:, :3].copy()
        self._cell = c

    def args(self, lo=None, hi=None, margin=0.0):
        return self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC

    def label(self, P):
        return np.asarray(self._cell.labels(np.atleast_2d(np.asarray(P, np.float64))), np.int64)

    def mems(self, P):
        return self._cell.mems(np.atleast_2d(np.asarray(P, np.float64)))


def library_manifest(bdir):
    import pandas as pd
    L = pd.read_csv(os.path.join(bdir, "library", "library_manifest.tsv"), sep="\t")
    L = L[L.path.str.endswith(".npz")]
    return dict(zip(L.species, L.path))


def libpack(ids, LIB):
    """library pack (as yeast wc3d.build_libpack with the library's own sphere trees)"""
    st_off = [0]; C, R, BR = [], [], []
    for i in ids:
        z = np.load(LIB[i])
        c = z["st_c"].astype(np.float64); r = z["st_r"].astype(np.float64)
        C.append(c); R.append(r); st_off.append(st_off[-1] + len(r))
        BR.append(float((np.linalg.norm(c, axis=1) + r).max()))
    return dict(ids=np.array(list(ids)), st_off=np.array(st_off, np.int64), STC=np.concatenate(C), STR=np.concatenate(R),
                BR=np.array(BR))


def _frozen_obstacles(cfg, bdir, K, LIB, spec, log, drop=None):
    """obstacle spheres of one frozen instance layer (rows of spec.npz selected by kind); drop: row indices of that file
    left out (instances dropped by a builder, drop_frozen_overlaps)"""
    from wholecell.preplaced.stage import _resolve
    p = _resolve(cfg, bdir, spec["npz"])
    z = np.load(p, allow_pickle=True)
    keep = np.isin(z[spec.get("kind_key", "kind")], spec["kinds"]) if spec.get("kinds") is not None else np.ones(len(z["types"]), bool)
    nd = 0
    if drop is not None and len(drop):
        drop = np.asarray(drop, np.int64); nd = int(keep[drop].sum()); keep[drop] = False
    ids = [str(x) for x in z["lib_ids"]]
    t = z["types"][keep]
    lp = libpack(ids, LIB)
    pos, q = z["pos"][keep], z["quat"][keep]
    if spec.get("shape", "tree") == "bounding":            # one bounding sphere per instance (memory; e.g. 14 M actin subunits)
        C, R = pos.astype(np.float64), lp["BR"][t].astype(np.float64)
    else:
        C, R = K.p3.emit_spheres(t.astype(np.int32), pos, q, lp["st_off"], lp["STC"], lp["STR"], lp["BR"], np.full(3, -1e7), np.full(3, 1e7))
    log(f"obstacles {spec['name']}: {int(keep.sum()):,} instances of {os.path.basename(p)} (kinds {spec.get('kinds')}), {len(R):,} spheres"
        + (f" ({nd:,} dropped)" if nd else ""))
    return (np.asarray(C, np.float64), np.asarray(R, np.float64)), dict(name=spec["name"], npz=p, n=int(keep.sum()),
                                                                       lib_ids=[ids[k] for k in np.unique(t)])


def _set_frozen_obs(ctx, names=None):
    """(re)build the obstacle sets fed by the frozen layers (context_params.obstacles, in config order), leaving out the
    dropped rows (ctx.frozen_drops); names: only these obstacle sets"""
    built = {}
    ctx.frozen_layers = [] if names is None else ctx.frozen_layers
    for k, spec in enumerate(ctx.frozen_specs):
        name = spec.get("obs", spec["name"])
        if names is not None and name not in names:
            continue
        p = ctx.frozen_paths[k]
        (C, R), info = _frozen_obstacles(ctx.cfg, ctx.bdir, ctx.K, ctx.LIB, spec, ctx.log, ctx.frozen_drops.get(p))
        if name in built:
            C0, R0 = built[name]; C, R = np.concatenate([C0, C]), np.concatenate([R0, R])
        built[name] = (C, R)
        if names is None:
            ctx.frozen_layers.append(info)
        else:
            ctx.frozen_layers[k].update(n=info["n"])
    ctx.obs.update(built)


def drop_frozen_overlaps(ctx, types, pos, quat, obs_names, tol=0.0, why="", chunk=64):
    """Instances of the frozen obstacle layers (context_params.obstacles) feeding the obstacle sets obs_names that
    overlap the given instances (types in ctx.lp; library sphere trees of both, sphere overlap > tol nm) are dropped:
    the rows go to ctx.frozen_drops[<frozen npz>] (left out of these obstacle sets, rebuilt, and of every frozen_subset
    export of the same file: its consumption shrinks, the copies return to the pack). For instances whose position is
    fixed by the geometry (NPCs at the package pores) that a frozen layer of another geometry version did not avoid.
    Returns {layer name: rows dropped}."""
    from scipy.spatial import cKDTree
    from wholecell.validate.placement import _spheres
    lp = ctx.lp
    types = np.asarray(types, np.int64); pos = np.asarray(pos, np.float64); quat = np.asarray(quat, np.float64)
    out = {}
    for k, spec in enumerate(ctx.frozen_specs):
        if spec.get("obs", spec["name"]) not in obs_names:
            continue
        p = ctx.frozen_paths[k]
        z = np.load(p, allow_pickle=True)
        keep = np.isin(z[spec.get("kind_key", "kind")], spec["kinds"]) if spec.get("kinds") is not None else np.ones(len(z["types"]), bool)
        prev = ctx.frozen_drops.get(p, np.zeros(0, np.int64))
        keep[prev] = False
        rows = np.nonzero(keep)[0]
        flp = libpack([str(x) for x in z["lib_ids"]], ctx.LIB)
        ft = z["types"][rows].astype(np.int64); fP = z["pos"][rows].astype(np.float64)
        tree = cKDTree(fP)
        cand = tree.query_ball_point(pos, lp["BR"][types] + float(flp["BR"].max()) + tol)
        hit = set()
        for i in np.nonzero([len(c) > 0 for c in cand])[0]:
            c = np.asarray(cand[i], np.int64)
            c = c[np.linalg.norm(fP[c] - pos[i], axis=1) < lp["BR"][types[i]] + flp["BR"][ft[c]]]   # bounding spheres meet
            if not len(c):
                continue
            Ca, Ra, _ = _spheres(types[i:i + 1], pos[i:i + 1], quat[i:i + 1], lp["st_off"], lp["STC"], lp["STR"], 0)
            Cb, Rb, Ib = _spheres(ft[c], fP[c], z["quat"][rows[c]].astype(np.float64), flp["st_off"], flp["STC"], flp["STR"], 0)
            tb = cKDTree(Cb)
            for a, js in enumerate(tb.query_ball_point(Ca, Ra + Rb.max())):
                if js:
                    js = np.asarray(js)
                    ov = np.linalg.norm(Cb[js] - Ca[a], axis=1) < Ra[a] + Rb[js] - tol
                    hit.update(int(x) for x in c[Ib[js[ov]]])
        drop = rows[np.array(sorted(hit), np.int64)] if hit else np.zeros(0, np.int64)
        ctx.frozen_drops[p] = np.union1d(prev, drop).astype(np.int64)
        ctx.frozen_layers[k].setdefault("dropped", 0)
        ctx.frozen_layers[k]["dropped"] += int(len(drop))
        ctx.frozen_layers[k].setdefault("dropped_why", why)
        out[spec["name"]] = int(len(drop))
        ctx.log(f"  frozen {spec['name']}: {len(drop):,} instances dropped ({why or 'overlap'}; sphere overlap > {tol} nm)")
    if any(out.values()):
        _set_frozen_obs(ctx, {s.get("obs", s["name"]) for s in ctx.frozen_specs if s["name"] in out})
    return out


@register("preplaced_context", "core.cellgeom")
def context(cfg, bdir, params):
    import pandas as pd
    from wholecell import seeds
    from wholecell.pack import geomload
    from wholecell.preplaced import kern
    from wholecell.preplaced.serial import Context, log
    geom = geomload.load(cfg, bdir)
    cell = CoreCell(geom)
    K = kern.core()
    LIB = library_manifest(bdir)
    ids = list(params.get("libs") or [])
    missing = [i for i in ids if i not in LIB]
    if missing:
        raise KeyError(f"preplaced context: library species not in {bdir}/library/library_manifest.tsv: {missing}")
    lp = libpack(ids, LIB)
    I = pd.read_csv(os.path.join(bdir, "inventory", "pack_table.tsv"), sep="\t")
    exp3 = pd.DataFrame(dict(lib_id=I.ing_id.astype(str), ing_id=I.ing_id.astype(str), region="all",
                             expected=I.copies.astype(float)))
    log(f"core context: geometry {geom.manifest.get('package_id')}, {len(cell.pores)} pores, {len(ids)} library types, "
        f"{len(exp3)} pack-table rows")

    def emit(types, pos, q, lo=None, hi=None):
        return K.p3.emit_spheres(np.asarray(types, np.int32), np.asarray(pos), np.asarray(q), lp["st_off"], lp["STC"],
                                 lp["STR"], lp["BR"], np.full(3, -1e7) if lo is None else lo, np.full(3, 1e7) if hi is None else hi)

    def seed(s):
        return seeds.seed_int(cfg, s) if isinstance(s, str) else int(s)

    ctx = Context(cell=cell, lp=lp, exp3=exp3, log=log, emit=emit, K=K, seed=seed, cfg=cfg, bdir=bdir, geom=geom,
                  lib_path=lambda i: LIB[i], pre_instances=None)
    ctx.groups = []
    ctx._pre_instances = lambda pre=None: _pre_instances(ctx)
    from wholecell.preplaced.stage import _resolve
    ctx.LIB = LIB
    ctx.frozen_specs = list(params.get("obstacles") or [])
    ctx.frozen_paths = [_resolve(cfg, bdir, spec["npz"]) for spec in ctx.frozen_specs]
    ctx.frozen_drops = {}                   # frozen npz -> row indices dropped by builders (drop_frozen_overlaps)
    _set_frozen_obs(ctx)
    return ctx


def add_group(ctx, name, types, pos, quat, **extra):
    """record an instance group placed by a builder (pre_instances, module exports)"""
    ctx.groups = [g for g in getattr(ctx, "groups", []) if g["name"] != name]
    ctx.groups.append(dict(name=name, types=np.asarray(types, np.int32), pos=np.asarray(pos, np.float32),
                           quat=np.asarray(quat, np.float32), **extra))


def _pre_instances(ctx):
    G = getattr(ctx, "groups", [])
    if not G:
        return np.zeros(0, np.int32), np.zeros((0, 3), np.float32), np.zeros((0, 4), np.float32)
    return (np.concatenate([g["types"] for g in G]), np.concatenate([g["pos"] for g in G]),
            np.concatenate([g["quat"] for g in G]))


def _frozen_subset(ctx, fs, f):
    """rows of a frozen pre-placed file selected by a per-instance (or per-thread) key: instance layers (lib_ids, types,
    pos, quat, ...) or bead sets (pts, off, per-thread and per-bead arrays)"""
    from wholecell.preplaced.stage import _resolve
    src = _resolve(ctx.cfg, ctx.bdir, fs["npz"])
    z = np.load(src, allow_pickle=False)
    key, vals = fs["key"], fs["values"]
    if "pts" in z.files:                                    # bead set: threads selected, offsets rebuilt
        off = z["off"]; nthr = len(off) - 1
        sel = np.isin(z[key], vals)
        lens = np.diff(off)[sel]
        bead_sel = np.repeat(sel, np.diff(off))
        out = dict(pts=z["pts"][bead_sel], off=np.r_[0, np.cumsum(lens)].astype(np.int64))
        for k in z.files:
            if k in ("pts", "off"):
                continue
            a = z[k]
            if a.ndim and len(a) == nthr:
                out[k] = a[sel]
            elif a.ndim and len(a) == len(bead_sel):
                out[k] = a[bead_sel]
            else:
                out[k] = a
        np.savez(f, **out)
        return dict(n_threads=int(sel.sum()), n_beads=int(bead_sel.sum()), consumed={}, frozen=src, select={key: vals})
    sel = np.isin(z[key], vals)
    dropped = np.asarray(getattr(ctx, "frozen_drops", {}).get(src, np.zeros(0, np.int64)), np.int64)
    dropped = dropped[sel[dropped]] if len(dropped) else dropped
    sel[dropped] = False                                    # rows dropped by a builder (drop_frozen_overlaps)
    ids = [str(x) for x in z["lib_ids"]]
    out = {}
    for k in z.files:                                       # read each array once (lazy npz)
        a = z[k]
        out[k] = a[sel] if (a.ndim and len(a) == len(sel)) else a
    np.savez(f, **out)
    l2i = fs.get("lib_to_ingredient") or {}
    cnt = np.bincount(out["types"], minlength=len(ids))
    consumed = {}
    for i, n_ in zip(ids, cnt):
        if n_:
            consumed[l2i.get(i, i)] = consumed.get(l2i.get(i, i), 0) + int(n_)
    j = dict(n=int(sel.sum()), consumed=consumed, frozen=src, select={key: vals})
    if len(dropped):
        np.save(f[:-4] + "_dropped_rows.npy", dropped)
        j.update(dropped=int(len(dropped)), dropped_rows=os.path.basename(f[:-4]) + "_dropped_rows.npy (row indices of the frozen file)",
                 dropped_why=[x.get("dropped_why") for x in getattr(ctx, "frozen_layers", []) if x.get("npz") == src and x.get("dropped")])
    return j


def export_modules(ctx, sc, out_dir):
    """serial `exports`: one pre-placed module file pair per entry, in the module contract of the pack index
    (npz: lib_ids, types, pos, quat [+ extras]; json: consumed copies per ingredient) or a bead set (npz: pts, off,
    bead_radius [+ extras]). Entry: {module, group: <instance group name> | beads: {pts, off, radius, extra keys}}"""
    os.makedirs(os.path.join(out_dir, "modules"), exist_ok=True)
    rows = []
    groups = {g["name"]: g for g in getattr(ctx, "groups", [])}
    consumed = getattr(ctx, "consumed", {})
    for e in sc.get("exports") or []:
        f = os.path.join(out_dir, "modules", f"{e['module']}.npz")
        if e.get("group"):
            g = groups[e["group"]]
            ids_all = [str(x) for x in ctx.lp["ids"]]
            used = np.unique(g["types"])
            remap = {int(t): k for k, t in enumerate(used)}
            types = np.array([remap[int(t)] for t in g["types"]], np.int32)
            extra = {k: v for k, v in g.items() if k not in ("name", "types", "pos", "quat", "consumed")}
            np.savez(f, lib_ids=np.array([ids_all[t] for t in used]), types=types, pos=g["pos"], quat=g["quat"], **extra)
            cons = g.get("consumed") or {}
            j = dict(module=e["module"], n=int(len(types)), consumed={k: int(round(v)) for k, v in cons.items()},
                     builder_group=e["group"], source="preplaced serial (core builders)")
        elif e.get("frozen_subset"):
            j = _frozen_subset(ctx, e["frozen_subset"], f)
            j.update(module=e["module"], source="frozen subset (selected at build time; file not copied)")
        else:
            b = e["beads"]
            np.savez(f, pts=ctx.pre[b["pts"]], off=ctx.pre[b["offsets"]], bead_radius=np.float32(b["radius"]),
                     **{k: ctx.pre[v] for k, v in (b.get("extra") or {}).items()})
            j = dict(module=e["module"], n_beads=int(len(ctx.pre[b["pts"]])), consumed={}, source="preplaced serial (core builders)")
        json.dump(j, open(f[:-4] + ".json", "w"), indent=1)
        rows.append(dict(module=e["module"], npz=f, **{k: j[k] for k in ("consumed",)}))
    return rows
