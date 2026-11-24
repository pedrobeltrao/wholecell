"""Export stage: python -m wholecell export --cell C --build B --step S [key=value ...]

Steps (outputs in <build>/export/ unless out_dir= is given):
  regions       region table of the build (config pack.regions order, then pre-pass groups) + display data from config
                export.regions {name: {colour, group, display_as}} -> regions.json. Replaces the legacy hard-coded
                RNAMES_V5/6/7 and jenc.region_map.
  cell_entry    this build's entry for the viewer registry cells.json (docs/reference/viewer_multi_cell.md section 2):
                id, name, organism, cell_type, state, data_base, build {build_id, stamp, stamp_sha256, aliases}, files,
                stats, licence, status, default_bookmark, model_doc -> cell.json
  tiles_v2      YCT2 stream tiles (docs/reference/tile_format_v2.md): web/tiles/s_BX_BY_BZ[.k].bin shards + index.json(.gz),
                web/regions.json, web/library/types.json, web/encode_meta.json, web/encode_check.json (round trip, counts),
                web/version.json, web/MANIFEST.sha256. Optional target for large cells; blocks=<n> or blocks=bx,by,bz;...
                encodes a subset (test), selftest=1 re-encodes into a temp dir and compares the manifests (determinism).
                Tiles carry PALETTE + INST, GHOST, POLY (mRNA, DNA, extra) and MESH sections (meshpoly.py; switches
                mesh=0 / polylines=0 / ghosts=0 or config export.tiles.{mesh, polylines, ghosts, mesh_h_nm, mesh_skip,
                lipid_mask_regions, workers}); blocks are encoded in parallel (workers=, default 8), one task per
                (bx, by) column; output is independent of the worker count.
  web           = tiles_v2. Standalone re-export of an existing build: python -m wholecell export --build B --step web
                --out DIR (--cell defaults to the build's config.resolved.json); DIR is the web folder (tiles/, library/,
                regions.json, ...); the build folder is only read.
  sim_manifest  one small simulation box (docs/reference/simulation_workflows.md section 3, schema wholecell.simbox/1;
                stub writer): boxes/<name>/{box.json, molecules.parquet, bodies.npz, walls.npz, MANIFEST.sha256}.
                Monomer molecules carry the SEC-MS bound-fraction attribute of their protein (inventory
                secms_bound_fraction.tsv; columns secms_bound_fraction, secms_unrepresented_share_of_free; box.json
                `secms_bound_fraction` = provenance + intended use); metadata only, the box is otherwise unchanged.
                Box: config export.boxes[0] or box_tile=<pack tile id> edge_nm=<nm> (default: the fullest tile of the
                first volume region, 200 nm). No pickles; quaternions (w,x,y,z).
  all           regions cell_entry sim_manifest (+ tiles_v2 when config export.tiles.enabled is true)
Common arguments: as the validation stage (pack_dir=, tiles_dir=, ref=, geom=, lib_dir=, regions_json=), out_dir=,
overwrite=1.
"""
from __future__ import annotations

import os

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), ".numba_cache"))   # never into a legacy __pycache__

import fnmatch  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import shutil  # noqa: E402
import tempfile  # noqa: E402
import time  # noqa: E402
import zlib  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from ..validate import common, regions as R_  # noqa: E402
from . import meshpoly as MP, quat as Q, yct2  # noqa: E402

KIND_GROUP = {"volume": "soluble", "virtual": "soluble", "membrane": "membrane", "shell": "shell", "preplaced": "prepass"}
# Viewer format (docs/reference/tile_format_v2.md) names: pre-placed regions, groups, flags and the `src` of box molecules
# stay "prepass" (the stage was renamed at A0c, the published format was not). Config region kind -> viewer kind:
VIEWER_KIND = {"preplaced": "prepass"}
DEFAULT_COLOURS = ["#d9d2c3", "#9fb7d9", "#e3a6a6", "#a6d3a0", "#d6b3e6", "#f2d08c", "#8cc7c7", "#c9a27e", "#b0b0b0",
                   "#e6c3d6", "#a3c4bc", "#d4a5a5", "#9c9ede", "#cedb9c", "#e7ba52", "#ad494a"]


def _out(bdir, kw, *sub):
    d = kw.get("out_dir") or os.path.join(bdir, "export")
    d = os.path.join(d, *sub)
    os.makedirs(d, exist_ok=True)
    return d


def _build_id(bdir):
    return os.path.basename(os.path.normpath(bdir))


def _preplaced_groups(cfg, PL):
    """pre-pass library -> group name (config export.preplaced_groups {group: [fnmatch patterns]}; default one group)"""
    groups = (cfg.get("export", {}) or {}).get("preplaced_groups") or {"prepass": ["*"]}
    out = {}
    for lib in PL:
        for g, pats in groups.items():
            if any(fnmatch.fnmatchcase(lib, p) for p in pats):
                out[lib] = g
                break
        else:
            out[lib] = "prepass"
    return out


def region_table(cfg, inp):
    ec = cfg.get("export", {}) or {}
    disp = ec.get("regions", {}) or {}
    regs = R_.load_regions(cfg, inp)
    rows = []
    for i, r in enumerate(regs):
        d = disp.get(r["name"], {}) or {}
        refs = {}
        if r.get("labels"): refs["labels"] = list(r["labels"])
        if r.get("membrane"): refs["membrane"] = r["membrane"]
        rows.append(dict(id=i, name=r["name"], kind=VIEWER_KIND.get(r["kind"], r["kind"]), geometry_refs=refs, display_as=d.get("display_as", r["name"]),
                         colour=d.get("colour", DEFAULT_COLOURS[i % len(DEFAULT_COLOURS)]), group=d.get("group", KIND_GROUP.get(r["kind"], r["kind"]))))
    PL = json.load(open(inp["preplaced_libs"])) if inp.has("preplaced_libs") else []
    if inp.has("preplaced_groups"):
        # per pre-pass entry group names written by the pack stage / an engine adapter (one library can belong to several
        # structures, e.g. yeast F-actin in cables and in endocytic patches); a group may reuse a pack region name
        pgl = list(json.load(open(inp["preplaced_groups"])))
        if len(pgl) != len(PL):
            raise ValueError("preplaced_groups.json must align with preplaced_libs.json")
    else:
        pg = _preplaced_groups(cfg, PL)
        pgl = [pg[l] for l in PL]
    have = {r["name"] for r in rows}
    for g in dict.fromkeys(pgl):
        libs_g = sorted({l for l, gg in zip(PL, pgl) if gg == g})
        if g in have:                                  # pre-pass group mapped onto an existing (kind preplaced) region
            for r in rows:
                if r["name"] == g:
                    r["libraries"] = libs_g
            continue
        i = len(rows); d = disp.get(g, {}) or {}
        rows.append(dict(id=i, name=g, kind=VIEWER_KIND["preplaced"], geometry_refs={}, display_as=d.get("display_as", g),
                         colour=d.get("colour", DEFAULT_COLOURS[i % len(DEFAULT_COLOURS)]), group=d.get("group", "prepass"),
                         libraries=libs_g))
    return rows, pgl


def step_regions(cfg, bdir, inp, kw):
    rows, _ = region_table(cfg, inp)
    out = dict(schema="wholecell.regions/1", cell_id=cfg["cell_id"], build_id=_build_id(bdir), regions=rows)
    p = os.path.join(_out(bdir, kw), "regions.json")
    common.write_json(p, out, overwrite=common.flag(kw, "overwrite"))
    print(f"regions.json: {len(rows)} regions")
    return out


def step_cell_entry(cfg, bdir, inp, kw):
    ec = cfg.get("export", {}) or {}
    bid = _build_id(bdir)
    sp = os.path.join(bdir, "BUILD_STAMP.json")
    stats = {}
    if inp.has("summary"):
        S = json.load(open(inp["summary"]))
        stats["instances"] = S.get("instances")
    if inp.has("libs"):
        stats["species"] = len(json.load(open(inp["libs"]))["libs"])
    if inp.has("tiles_parquet"):
        M = pd.read_parquet(inp["tiles_parquet"])
        regs = R_.load_regions(cfg, inp)
        vol = [r["name"] for r in regs if r["kind"] == "volume" and r["name"] in M.columns]
        stats["volume_um3"] = round(float(M[vol].sum().sum()) / 1e9, 1)
    files = {"index": "tiles/index.json", "regions": "regions.json", "geometry": "geometry.json", "fields": "geometry_fields.bin",
             "meshes": "organelle_meshes.bin.gz", "bookmarks": "bookmarks.json", "model_doc": "MODEL.md",
             "features": "features.json", "library": "library/types.json"}
    files.update(ec.get("files", {}) or {})
    entry = dict(id=cfg["cell_id"], name=cfg.get("name", cfg["cell_id"]), organism=cfg.get("organism"),
                 cell_type=ec.get("cell_type", ""), state=ec.get("state", ""),
                 data_base=str(ec.get("data_base", "./")).replace("{build_id}", bid).replace("{cell_id}", cfg["cell_id"]),
                 build=dict(build_id=bid, stamp="BUILD_STAMP.json",
                            stamp_sha256=hashlib.sha256(open(sp, "rb").read()).hexdigest() if os.path.exists(sp) else None,
                            aliases=sorted((cfg.get("aliases") or {}).keys())),
                 files=files, thumbnail=ec.get("thumbnail", "thumb.webp"), stats=stats,
                 licence=ec.get("licence", {"data": "unspecified"}), status=ec.get("status", "preview"),
                 default_bookmark=ec.get("default_bookmark", "Whole cell"),
                 model_doc=os.path.basename(str(ec.get("model_doc", "MODEL.md"))),
                 sim_defaults=ec.get("sim_defaults", {}))
    out = dict(schema="wholecell.cells/1", entry=entry)
    p = os.path.join(_out(bdir, kw), "cell.json")
    common.write_json(p, out, overwrite=common.flag(kw, "overwrite"))
    print(f"cell.json: {entry['id']} build {bid}")
    return out


# ------------------------------------------------------------------ tiles_v2
class _TileCache:
    def __init__(self, inp, n=64):
        self.inp = inp; self.n = n; self.d = {}

    def get(self, t):
        if t in self.d:
            return self.d[t]
        f = os.path.join(self.inp["tiles_dir"], f"{t}.npz")
        if os.path.exists(f):
            z = np.load(f)
            reg = z["region"] if "region" in z.files else np.repeat(z["job_region"].astype(np.int64), z["job_placed"])
            v = (z["types"].astype(np.int64), z["pos"].astype(np.float64), z["quat"].astype(np.float64), reg.astype(np.int64))
        else:
            v = None
        if len(self.d) >= self.n:
            self.d.pop(next(iter(self.d)))
        self.d[t] = v
        return v


def _flag_opt(kw, tc, name, default=True):
    v = kw.get(name)
    if v is not None:
        return str(v).lower() not in ("0", "false", "no", "off")
    return bool(tc.get(name, default))


def _lib_br(cfg, bdir, inp, kw, types, web):
    """library cache (bounding radius, bead mass) of the export types; the build's validation cache is read, never
    rewritten from a standalone export (out=): a missing / incomplete cache is then written next to the output"""
    cp = os.path.join(bdir, "validate", "cache", "lib_all.npz")
    if kw.get("out"):
        ok = False
        if os.path.exists(cp):
            z = np.load(cp, allow_pickle=False); ok = set(types) <= set(str(x) for x in z["names"])
        if not ok:
            cp = os.path.join(_work_dir(web), "lib_all.npz")
    try:
        return common.lib_cache(inp, types, cp, nproc=min(8, common.nprocs(cfg, kw)))
    except Exception as e:  # noqa: BLE001 - library absent: types without radii
        print("library cache unavailable:", e)
        return None


def _work_dir(web):
    return os.path.normpath(web) + "_work"


_W = {}


def _encode(cfg, bdir, inp, kw, web):
    ec = cfg.get("export", {}) or {}
    tc = ec.get("tiles", {}) or {}
    tile_nm = float(tc.get("tile_nm", 200)); bt = int(tc.get("block_tiles", 5))
    q_org = float(tc.get("q_org_nm", 80)); q_range = float(tc.get("q_range_nm", 360)); shard_max = float(tc.get("shard_max_mb", 45)) * 1e6
    meta = json.load(open(inp["tiles_meta"]))
    PT = float(meta.get("TILE", meta.get("tile", 500.0))); lo0 = np.asarray(meta["lo0"], np.float64)
    origin = np.asarray(tc.get("origin_nm", lo0), np.float64)
    B = bt * tile_nm
    if abs(B / PT - round(B / PT)) > 1e-9 or np.any(np.abs((origin - lo0) / PT - np.round((origin - lo0) / PT)) > 1e-9):
        raise ValueError(f"block edge {B} nm / grid origin must align with the pack tile ({PT} nm)")
    M = pd.read_parquet(inp["tiles_parquet"], columns=["ti", "tj", "tk"])
    key = {tuple(t): i for i, t in enumerate(M[["ti", "tj", "tk"]].values.tolist())}
    L = json.load(open(inp["libs"])); libs = list(L["libs"])
    mass = list(L.get("mass_da") or [])     # yeast (pack.requests legacy_x5) libs.json has no mass_da: library bead mass below
    PL = json.load(open(inp["preplaced_libs"])) if inp.has("preplaced_libs") else []
    types = libs + [l for l in PL if l not in libs]
    tidx = {n: i for i, n in enumerate(types)}
    pmap = np.array([tidx[l] for l in PL], np.int64) if PL else np.zeros(0, np.int64)
    rows, pg = region_table(cfg, inp)
    rid = {r["name"]: r["id"] for r in rows}
    pre_reg = np.array([rid[g] for g in pg], np.int64) if PL else np.zeros(0, np.int64)     # pg: group per pre-pass entry
    has_pre = inp.has("preplaced_off") and len(PL) > 0
    lib_tag = yct2.tag4(json.dumps(types)); bid = _build_id(bdir); build_tag = yct2.tag4(bid)
    # blocks: from the pack tile grid
    tij = M[["ti", "tj", "tk"]].values.astype(np.int64)
    pt_lo = tij * PT + lo0
    blk_all = np.unique(np.floor((pt_lo - origin) / B + 1e-9).astype(np.int64), axis=0)
    sel = kw.get("blocks")
    if sel:
        if ";" in str(sel) or "," in str(sel):
            want = [tuple(int(x) for x in s.split(",")) for s in str(sel).split(";") if s]
            blocks = np.array(want, np.int64)
        else:
            # first n blocks by pack-tile occupancy (deterministic): blocks with the most pack tiles
            cnt = {}
            for b in np.floor((pt_lo - origin) / B + 1e-9).astype(np.int64):
                cnt[tuple(b)] = cnt.get(tuple(b), 0) + 1
            blocks = np.array(sorted(cnt, key=lambda b: (-cnt[b], b))[:int(sel)], np.int64)
    else:
        blocks = blk_all
    lc = _lib_br(cfg, bdir, inp, kw, types, web)
    # MESH / GHOST / POLY sections (meshpoly.py); switches: config export.tiles.{mesh, ghosts, polylines} or key=0/1
    do_mesh = _flag_opt(kw, tc, "mesh"); do_poly = _flag_opt(kw, tc, "polylines"); do_ghost = do_mesh and _flag_opt(kw, tc, "ghosts")
    mask_names = [r["name"] for r in rows if r["kind"] == "membrane"] + [n for n in (tc.get("lipid_mask_regions") or []) if n in rid]
    mask_ids = np.array(sorted({rid[n] for n in mask_names}), np.int64)
    br = np.array([float(lc.br[lc.idx[n]]) if (lc is not None and n in lc.idx) else 15.0 for n in types], np.float64)
    cell = mesher = None; mesh_info = {}
    if do_mesh:
        from ..geometry import cellgeom
        from ..geometry.engines import cellgeom_eval_v1 as ev
        pj = json.load(open(os.path.join(bdir, "geometry", "package.json")))
        cell = cellgeom.load(pj["package_dir"], verify="quick")
        mems = MP.lipid_membranes(cell.manifest, skip=set(tc.get("mesh_skip", []) or []))
        zf = float(cell.G[ev.G_ZFACE]) if cell.G[ev.G_FACE] > 0 else None
        mesh_h = float(tc.get("mesh_h_nm", MP.MESH_H))
        mesher = MP.TileMesher(cell, mems, tile_nm, mesh_h, zf)
        mesh_info = dict(mesh_h_nm=mesh_h, membranes=[dict(code=c, name=n) for c, n in mems], open_face_z_nm=zf,
                         package_id=pj.get("package_id"), source="geometry package membrane SDFs (cellgeom_eval_v1 mems)")
    kmin = blk_all.min(0) * bt; kdim = (blk_all.max(0) + 1) * bt - kmin
    poly = None; poly_tiles = set()
    if do_poly:
        # polyline work files (memory-mapped, scattered writes): on a local temp disk (tmp_dir= or $TMPDIR), never on the
        # build's / output's network file system; removed after the export. poly_work=DIR: keep and reuse (dev)
        work = kw.get("poly_work") or tempfile.mkdtemp(prefix="wholecell_poly_", dir=kw.get("tmp_dir") or None)
        _W["poly_tmp"] = None if kw.get("poly_work") else work
        poly = MP.PolyStore(work, origin, tile_nm, kmin, kdim)
        if kw.get("poly_work") and os.path.exists(os.path.join(work, "poly_meta.json")):
            poly.stats = {int(k): v for k, v in json.load(open(os.path.join(work, "poly_meta.json"))).items()}
            print(f"polylines reused from {work}", flush=True)
        else:
            poly.build(MP.build_sources(cfg, bdir, inp.lib_dir()))
        poly.open(); poly_tiles = poly.tile_keys()
    tdir = os.path.join(web, "tiles"); os.makedirs(tdir, exist_ok=True)
    poly_tmp = _W.get("poly_tmp")
    _W.clear()
    _W.update(poly_tmp=poly_tmp, inp=inp, key=key, PT=PT, lo0=lo0, origin=origin, B=B, bt=bt, tile_nm=tile_nm, q_org=q_org, q_range=q_range,
              shard_max=shard_max, has_pre=has_pre, pmap=pmap, pre_reg=pre_reg, lib_tag=lib_tag, build_tag=build_tag,
              tdir=tdir, mesher=mesher, poly=poly, poly_tiles=poly_tiles, do_ghost=do_ghost, mask_ids=mask_ids, br=br)
    blist = [tuple(x) for x in blocks.tolist()]
    groups = {}
    for b in blist:                                       # one task per (bx, by) column: consecutive blocks share pack tiles
        groups.setdefault(b[:2], []).append(b)
    tasks = [sorted(v, key=lambda b: b[2]) for _, v in sorted(groups.items())]
    nw = max(1, min(int(kw.get("workers", tc.get("workers", 8))), len(tasks)))
    res = {}
    if nw == 1:
        for t in tasks:
            for r in _encode_column(t):
                res[r["block"]] = r
    else:
        import multiprocessing as mp_
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(nw, mp_context=mp_.get_context("fork"), initializer=_worker_init) as ex:
            for out in ex.map(_encode_column, tasks):
                for r in out:
                    res[r["block"]] = r
    tiles_rows = []; counts = {}; files = []; src_n = 0
    check = dict(max_pos_err_nm=0.0, max_rot_err_rad=0.0, n_checked=0)
    sec_stats = dict(tiles_mesh=0, n_vert=0, n_tri=0, n_ghost=0, n_poly_tiles=0, n_poly_beads=0, mesh_area_um2={})
    names = {c: n for c, n in (mesher.mem if mesher else [])}
    for b in blist:                                       # assemble in block order (deterministic, independent of workers)
        r = res.get(b)
        if r is None:
            continue
        tiles_rows += r["rows"]; files += r["files"]; src_n += r["src_n"]
        for k, v in r["counts"].items():
            counts[k] = counts.get(k, 0) + v
        check["max_pos_err_nm"] = max(check["max_pos_err_nm"], r["check"]["max_pos_err_nm"])
        check["max_rot_err_rad"] = max(check["max_rot_err_rad"], r["check"]["max_rot_err_rad"])
        check["n_checked"] += r["check"]["n_checked"]
        for k in ("tiles_mesh", "n_vert", "n_tri", "n_ghost", "n_poly_tiles", "n_poly_beads"):
            sec_stats[k] += r["sec"][k]
        for c, a in r["sec"]["area"].items():
            nm_ = names.get(c, str(c))
            sec_stats["mesh_area_um2"][nm_] = sec_stats["mesh_area_um2"].get(nm_, 0.0) + a / 1e6
    sec_stats["mesh_area_um2"] = {k: round(v, 3) for k, v in sorted(sec_stats["mesh_area_um2"].items())}
    check["sections"] = sec_stats
    if poly is not None:
        check["polylines"] = poly.stats
    allt = np.array([r[:3] for r in tiles_rows]) if tiles_rows else np.zeros((0, 3), int)
    written = ["PALETTE", "INST"] + (["GHOST"] if do_ghost else []) + (["POLY_MRNA", "POLY_LINK", "POLY_EXTRA"] if do_poly else []) \
        + (["MESH"] if do_mesh else [])
    not_written = [x for x in ("GHOST", "POLY", "MESH") if not any(w.startswith(x) for w in written)] + ["COLUMN", "LOD1 overview"]
    pk = {}
    if do_poly:
        for sec_, d in MP.SECTIONS.items():
            pk[str(sec_)] = dict(name=d["name"], poly_q=d["poly_q"], free_kind=d["free_kind"], piece_points=MP.PIECE,
                                 sources=[x["name"] for x in poly.stats.get(sec_, {}).get("sources", [])])
        pk["6"]["kinds"] = {str(k): v for k, v in MP.EXTRA_KINDS.items()}
    index = dict(format="YCT2", version="2.0", build_id=bid, tile_nm=tile_nm, block_tiles=bt,
                 grid=dict(origin_nm=origin.tolist(), min_ijk=allt.min(0).tolist() if len(allt) else None,
                           max_ijk=allt.max(0).tolist() if len(allt) else None),
                 q_org=q_org, q_range=q_range, lod=[0], flag_bits=yct2.CORE_FLAG_BITS + list(ec.get("flag_bits", []) or []),
                 columns=[], polyline_kinds=pk, sections_written=written,
                 not_written=not_written, tile_columns=yct2.TILE_COLUMNS,
                 lib_tag=lib_tag, build_tag=build_tag)
    if do_mesh:
        index["mesh"] = mesh_info
    if do_ghost:
        index["lipid_mask_regions"] = [int(x) for x in mask_ids]
    index["tiles"] = tiles_rows
    ij = json.dumps(index, separators=(",", ":"))
    open(os.path.join(tdir, "index.json"), "w").write(ij)
    open(os.path.join(tdir, "index.json.gz"), "wb").write(yct2.gz(ij.encode()))
    rr = dict(schema="wholecell.regions/1", cell_id=cfg["cell_id"], build_id=bid, regions=rows)
    json.dump(rr, open(os.path.join(web, "regions.json"), "w"), indent=1)
    os.makedirs(os.path.join(web, "library"), exist_ok=True)
    # mass per type: libs.json mass_da (Jurkat, testcell); else, and for pre-pass types, the library bead mass (aa 110 Da,
    # nt 330 Da), as the yeast engine computed it (webenc.encode_library, val3d.lib_beads) and validate.common.lib_masses
    tj = [dict(id=i, name=n, mass_da=(mass[i] if i < len(mass) else (float(lc.bead_mass[lc.idx[n]]) if lc else None)),
               br_nm=(round(float(lc.br[lc.idx[n]]), 3) if lc else None), prepass=n in PL and n not in libs) for i, n in enumerate(types)]
    json.dump(dict(schema="wholecell.types/1", ids=types, types=tj), open(os.path.join(web, "library", "types.json"), "w"))
    json.dump(dict(v=bid), open(os.path.join(web, "version.json"), "w"))
    return index, counts, check, src_n, files


def _worker_init():
    try:
        import numba
        numba.set_num_threads(1)
    except Exception:  # noqa: BLE001
        pass


def _encode_column(blocks):
    cache = _TileCache(_W["inp"])
    return [_encode_block(b, cache) for b in blocks]


def _encode_block(b, cache):
    W = _W
    inp, key, PT, lo0, origin, B, bt, tile_nm = W["inp"], W["key"], W["PT"], W["lo0"], W["origin"], W["B"], W["bt"], W["tile_nm"]
    q_org, q_range, shard_max, tdir = W["q_org"], W["q_range"], W["shard_max"], W["tdir"]
    has_pre, pmap, pre_reg = W["has_pre"], W["pmap"], W["pre_reg"]
    mesher, poly, br, mask_ids = W["mesher"], W["poly"], W["br"], W["mask_ids"]
    if has_pre and "poff" not in W:
        W["poff"] = np.load(inp["preplaced_off"]); W["ptypes"] = np.load(inp["preplaced_types"], mmap_mode="r")
        W["ppos"] = np.load(inp["preplaced_pos"], mmap_mode="r"); W["pquat"] = np.load(inp["preplaced_quat"], mmap_mode="r")
    poff = W.get("poff"); ptypes = W.get("ptypes"); ppos = W.get("ppos"); pquat = W.get("pquat")
    ppb = int(round(B / PT))
    out = dict(block=b, rows=[], files=[], counts={}, src_n=0, check=dict(max_pos_err_nm=0.0, max_rot_err_rad=0.0, n_checked=0),
               sec=dict(tiles_mesh=0, n_vert=0, n_tri=0, n_ghost=0, n_poly_tiles=0, n_poly_beads=0, area={}))
    check = out["check"]; counts = out["counts"]; sec = out["sec"]
    blo = origin + np.array(b) * B
    p0 = np.floor((blo - lo0) / PT + 1e-9).astype(np.int64)
    parts = []
    for di in range(-1, ppb + 1):                       # block pack tiles + a one-tile margin (centres may stray)
        for dj in range(-1, ppb + 1):
            for dk in range(-1, ppb + 1):
                t = key.get((int(p0[0] + di), int(p0[1] + dj), int(p0[2] + dk)))
                if t is None: continue
                v = cache.get(t)
                if v is not None and len(v[0]):
                    parts.append((v[0], v[3], np.zeros(len(v[0]), np.uint32), v[1], v[2]))
                if has_pre and poff[t + 1] > poff[t]:
                    a, c = int(poff[t]), int(poff[t + 1]); pt = np.asarray(ptypes[a:c]).astype(np.int64)
                    parts.append((pmap[pt], pre_reg[pt], np.full(c - a, yct2.FLAG_PREPASS, np.uint32),
                                  np.asarray(ppos[a:c], np.float64), np.asarray(pquat[a:c], np.float64)))
    if parts:
        sp = np.concatenate([p[0] for p in parts]); rg = np.concatenate([p[1] for p in parts]); fl = np.concatenate([p[2] for p in parts])
        ps = np.concatenate([p[3] for p in parts]); qu = np.concatenate([p[4] for p in parts])
    else:
        sp = np.zeros(0, np.int64); rg = np.zeros(0, np.int64); fl = np.zeros(0, np.uint32)
        ps = np.zeros((0, 3)); qu = np.zeros((0, 4))
    ijk_all = np.floor((ps - origin) / tile_nm).astype(np.int64)
    # lipid-mask pool for ghosts: instances of the lipid-mask regions in the block + margin
    gm = np.isin(rg, mask_ids) if (W["do_ghost"] and len(rg)) else np.zeros(len(rg), bool)
    g_sp, g_rg, g_fl, g_ps, g_qu, g_ijk = sp[gm], rg[gm], fl[gm], ps[gm], qu[gm], ijk_all[gm]
    g_r = br[g_sp] + 3.0 if len(g_sp) else np.zeros(0)
    inb = np.all((ijk_all >= np.array(b) * bt) & (ijk_all < (np.array(b) + 1) * bt), axis=1)
    sp, rg, fl, ps, qu, ijk = sp[inb], rg[inb], fl[inb], ps[inb], qu[inb], ijk_all[inb]
    out["src_n"] = int(len(sp))
    o = np.lexsort((ijk[:, 2], ijk[:, 1], ijk[:, 0]))
    sp, rg, fl, ps, qu, ijk = sp[o], rg[o], fl[o], ps[o], qu[o], ijk[o]
    brk = np.r_[0, np.nonzero(np.any(np.diff(ijk, axis=0) != 0, axis=1))[0] + 1, len(sp)] if len(sp) else np.zeros(1, np.int64)
    inst_rng = {tuple(ijk[a].tolist()): (a, c) for a, c in zip(brk[:-1], brk[1:]) if c > a}
    # candidate tiles: instances, polyline pieces, and (mesh) every tile of the block
    b0 = np.array(b) * bt
    cand = set(inst_rng)
    cand |= {t for t in W["poly_tiles"] if all(b0[i] <= t[i] < b0[i] + bt for i in range(3))}
    if mesher is not None:
        cand |= {(int(b0[0] + i), int(b0[1] + j), int(b0[2] + k)) for i in range(bt) for j in range(bt) for k in range(bt)}
    shard_k = 0; buf = []; size = 0; members = []
    for t_ijk in sorted(cand):
        lo = origin + np.array(t_ijk) * tile_nm
        a, c = inst_rng.get(t_ijk, (0, 0))
        mesh = None; nv = nt = 0
        if mesher is not None:
            V, F, FM = mesher(lo)
            if len(F):
                mb, nv, nt, big = MP.encode_mesh(V, F, FM, lo, q_org, q_range)
                mesh = (mb, big, nt)
                tri = V[F]; ar = 0.5 * np.linalg.norm(np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)
                for code in np.unique(FM):
                    sec["area"][int(code)] = sec["area"].get(int(code), 0.0) + float(ar[FM == code].sum())
        polys, npb = poly.tile_sections(t_ijk) if poly is not None else ({}, 0)
        if c <= a and mesh is None and not polys:
            continue
        ghost = None
        if mesh is not None and len(g_sp):
            hi = lo + tile_nm
            ok = np.all((g_ps > lo - g_r[:, None]) & (g_ps < hi + g_r[:, None]), axis=1) & ~np.all(g_ijk == np.array(t_ijk), axis=1)
            ok &= np.all((g_ps > lo - q_org + 1) & (g_ps < lo - q_org + q_range - 1), axis=1)
            if ok.any():
                ghost = (g_sp[ok], g_rg[ok], g_fl[ok] | yct2.FLAG_GHOST | yct2.FLAG_NEAR_MEMBRANE, g_ps[ok], g_qu[ok])
        raw, n_i, n_p = yct2.encode_tile(list(t_ijk), lo, sp[a:c], rg[a:c], fl[a:c], ps[a:c], qu[a:c], tile_nm, q_org, q_range,
                                         W["lib_tag"], W["build_tag"], ghost=ghost, polys=polys, mesh=mesh)
        if c > a:                                         # round trip on every tile
            D = yct2.decode_tile(raw)
            pos_d = origin + D["lo_rel"] + D["pos_rel"]
            o2 = np.lexsort((rg[a:c], sp[a:c]))
            check["max_pos_err_nm"] = max(check["max_pos_err_nm"], float(np.abs(pos_d - ps[a:c][o2]).max()))
            check["max_rot_err_rad"] = max(check["max_rot_err_rad"], float(Q.angle_between(D["quat"], qu[a:c][o2]).max()))
            check["n_checked"] += int(c - a)
            for s_, r_ in zip(D["species"].tolist(), D["region"].tolist()):
                counts[(s_, r_)] = counts.get((s_, r_), 0) + 1
        n_g = 0 if ghost is None else len(ghost[0])
        sec["n_ghost"] += n_g; sec["n_vert"] += nv; sec["n_tri"] += nt; sec["tiles_mesh"] += int(nt > 0)
        sec["n_poly_beads"] += npb; sec["n_poly_tiles"] += int(npb > 0)
        g = yct2.gz(raw)
        if buf and size + len(g) > shard_max:
            members.append((shard_k, buf)); shard_k += 1; buf = []; size = 0
        buf.append((list(t_ijk), g, n_i, n_g, npb, nv, nt, yct2.crc32(raw))); size += len(g)
    if buf:
        members.append((shard_k, buf))
    single = len(members) == 1
    for k, mem in members:
        fn = f"s_{b[0]}_{b[1]}_{b[2]}.bin" if single else f"s_{b[0]}_{b[1]}_{b[2]}.{k}.bin"
        off = 0
        with open(os.path.join(tdir, fn), "wb") as fh:
            for t_ijk, g, n_i, n_g, npb, nv, nt, crc in mem:
                fh.write(g)
                out["rows"].append([*t_ijk, fn, off, len(g), n_i, n_g, npb, nv, nt, 0, crc])
                off += len(g)
        out["files"].append(fn)
    return out


def _manifest(web):
    rows = []
    for root, dirs, fs in os.walk(web):
        dirs.sort()
        for f in sorted(fs):
            if f == "MANIFEST.sha256": continue
            p = os.path.join(root, f)
            rows.append(f"{hashlib.sha256(open(p, 'rb').read()).hexdigest()}  {os.path.relpath(p, web)}")
    open(os.path.join(web, "MANIFEST.sha256"), "w").write("\n".join(rows) + "\n")
    data = [r for r in rows if not r.endswith(("encode_meta.json", "encode_check.json"))]   # run-specific files
    return hashlib.sha256("\n".join(data).encode()).hexdigest()


def step_tiles_v2(cfg, bdir, inp, kw):
    if kw.get("out"):                                     # standalone re-export of an existing build: out= is the web folder
        web = os.path.abspath(kw["out"]); os.makedirs(web, exist_ok=True)
        if os.path.realpath(web).startswith(os.path.realpath(bdir) + os.sep) and not common.flag(kw, "overwrite"):
            raise ValueError(f"out={web} lies inside the build folder (builds are write-once): choose another folder")
    else:
        web = _out(bdir, kw, "web")
    if os.listdir(web) and not common.flag(kw, "overwrite"):
        raise FileExistsError(f"{web} is not empty (export outputs are write-once; overwrite=1 to replace)")
    if os.listdir(web):                                   # overwrite: start from an empty folder (no stale shards)
        shutil.rmtree(web); os.makedirs(web)
    t0 = time.time()
    index, counts, check, src_n, files = _encode(cfg, bdir, inp, kw, web)
    if not common.flag(kw, "keep_work"):
        shutil.rmtree(_work_dir(web), ignore_errors=True)
        if _W.get("poly_tmp"):
            shutil.rmtree(_W["poly_tmp"], ignore_errors=True)
    n_dec = sum(counts.values())
    tc = (cfg.get("export", {}) or {}).get("tiles", {}) or {}
    q_range = float(tc.get("q_range_nm", 360))
    ok = dict(count_equal=(n_dec == src_n), pos_ok=check["max_pos_err_nm"] <= q_range / 131070 + 1e-4,
              rot_ok=check["max_rot_err_rad"] <= 1e-4)
    shard_mb = max((os.path.getsize(os.path.join(web, "tiles", f)) for f in files), default=0) / 1e6
    ec = dict(n_tiles=len(index["tiles"]), n_shards=len(files), instances_source=int(src_n), instances_decoded=int(n_dec),
              max_shard_mb=round(shard_mb, 3), **check, **ok, blocks=kw.get("blocks", "all"))
    json.dump(ec, open(os.path.join(web, "encode_check.json"), "w"), indent=1)
    meta = dict(build_id=_build_id(bdir), cell_id=cfg["cell_id"], config_sha256=cfg.config_hash() if hasattr(cfg, "config_hash") else None,
                encoder="wholecell.export.yct2", encoder_commit=__import__("wholecell.build", fromlist=["git_info"]).git_info()["commit7"],
                zlib_version=zlib.ZLIB_VERSION, gzip=dict(level=6, mtime=0), seconds=round(time.time() - t0, 1),
                inputs=dict(tiles_dir=inp["tiles_dir"], libs=inp["libs"]), sections=index.get("sections_written"),
                mesh=index.get("mesh"), polyline_kinds=index.get("polyline_kinds"),
                workers=int(kw.get("workers", ((cfg.get("export", {}) or {}).get("tiles", {}) or {}).get("workers", 8))))
    json.dump(meta, open(os.path.join(web, "encode_meta.json"), "w"), indent=1)
    digest = _manifest(web)
    res = dict(encode_check=ec, manifest_digest_excl_meta=digest)
    if common.flag(kw, "selftest"):
        tmp = tempfile.mkdtemp(prefix="yct2_selftest_", dir=kw.get("tmp_dir") or None)
        try:
            kw2 = dict(kw, out_dir=tmp, selftest="0", overwrite="1")
            kw2.pop("out", None)
            w2 = os.path.join(tmp, "web"); os.makedirs(w2, exist_ok=True)
            _encode(cfg, bdir, inp, kw2, w2)
            shutil.rmtree(_work_dir(w2), ignore_errors=True)
            if _W.get("poly_tmp"):
                shutil.rmtree(_W["poly_tmp"], ignore_errors=True)
            d2 = _manifest(w2)
            ec["determinism_reencode_identical"] = d2 == digest
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        json.dump(ec, open(os.path.join(web, "encode_check.json"), "w"), indent=1)
        _manifest(web)
    print(json.dumps(ec, indent=1))
    return res


# ------------------------------------------------------------------ simulation box (stub)
SECMS_BOUND_USE = ("metadata only, no simulation change: SEC-MS measured assembled fraction of the protein (a_k) and the "
                   "share of the modelled free (monomer) copies that SEC-MS measures as bound in assemblies the structural "
                   "model does not represent. Intended use (docs/methods_secms_bound_fraction.md): treat each monomer "
                   "molecule as a mixture of a free state and a bound state with probability secms_unrepresented_share_of_free "
                   "(e.g. larger effective hydrodynamic radius from apex_app_kda, or reduced mobility / transient tethering "
                   "for the bound share); NaN = not a monomer species, not measured, or no SEC-MS step.")


def secms_bound_columns(species, inv_dir):
    """per box species (library id): (secms_bound_fraction, secms_unrepresented_share_of_free) of monomer species, and the
    box.json metadata block; NaN for assemblies, pre-pass structures and unmeasured proteins"""
    tab = os.path.join(inv_dir, "secms_bound_fraction.tsv")
    nan = {s: (np.nan, np.nan) for s in species}
    if not os.path.exists(tab):
        return nan, dict(status="unavailable (inventory has no secms_bound_fraction.tsv)", intended_use=SECMS_BOUND_USE)
    F = pd.read_csv(tab, sep="\t", dtype={"gene": str}).drop_duplicates("gene").set_index("gene")
    P = pd.read_csv(os.path.join(inv_dir, "pack_table.tsv"), sep="\t", usecols=["ing_id", "kind", "genes"])
    omp = os.path.join(inv_dir, "assembly", "oligo_map.json")
    om = json.load(open(omp)) if os.path.exists(omp) else {}
    mono = P[(P.kind == "protein") & ~P.ing_id.isin(list(om))].set_index("ing_id").genes.astype(str).to_dict()
    out = {}
    for s in species:
        g = mono.get(s)
        if g is not None and g in F.index:
            out[s] = (float(F.at[g, "secms_bound_fraction"]), float(F.at[g, "unrepresented_share_of_free"]))
        else:
            out[s] = (np.nan, np.nan)
    m = F[F.secms_bound_fraction.notna()]
    meta = dict(status="ok", table="inventory/secms_bound_fraction.tsv",
                table_sha256=hashlib.sha256(open(tab, "rb").read()).hexdigest(),
                secms_source_id=";".join(sorted({str(x) for x in m.secms_source_id.dropna()})),
                threshold_k=sorted({int(x) for x in m.threshold_k.dropna()}),
                columns=dict(secms_bound_fraction="measured assembled fraction a_k of the protein (whole protein, 0-1)",
                             secms_unrepresented_share_of_free="max(0, measured - modelled) / (1 - modelled): share of "
                             "this protein's monomer copies measured bound in unmodelled assemblies (0-1)"),
                species_with_value=int(sum(1 for v in out.values() if np.isfinite(v[0]))), intended_use=SECMS_BOUND_USE)
    return out, meta


def step_sim_manifest(cfg, bdir, inp, kw):
    from ..validate.geom import open_geometry
    ec = cfg.get("export", {}) or {}
    sd = dict(temperature_K=310.15, viscosity_Pa_s=6.9e-4, ionic_strength_M=0.15); sd.update(ec.get("sim_defaults", {}) or {})
    box_cfg = (ec.get("boxes") or [{}])[0]
    M = pd.read_parquet(inp["tiles_parquet"])
    meta = json.load(open(inp["tiles_meta"])); PT = float(meta.get("TILE", 500.0)); lo0 = np.asarray(meta["lo0"], float)
    regs = R_.load_regions(cfg, inp); rn = R_.names(regs)
    edge = float(kw.get("edge_nm", box_cfg.get("edge_nm", 200.0)))
    if box_cfg.get("centre_nm") is not None:
        centre = np.asarray(box_cfg["centre_nm"], float); how = "config export.boxes[0]"
    else:
        if kw.get("box_tile") is not None:
            t = int(kw["box_tile"]); how = f"box_tile={t}"
        else:
            vol0 = next(r["name"] for r in regs if r["kind"] == "volume" and r["name"] in M.columns)
            t = int(np.argmax(M[vol0].values)); how = f"fullest tile of {vol0} (tile {t})"
        centre = M.loc[t, ["ti", "tj", "tk"]].values.astype(float) * PT + lo0 + PT / 2
    name = kw.get("box_name", box_cfg.get("name", "box0"))
    lo = centre - edge / 2; hi = centre + edge / 2
    out = _out(bdir, kw, "boxes", name)
    if os.listdir(out) and not common.flag(kw, "overwrite"):
        raise FileExistsError(f"{out} not empty (overwrite=1 to replace)")
    key = {tuple(x): i for i, x in enumerate(M[["ti", "tj", "tk"]].values.tolist())}
    t0 = np.floor((lo - lo0) / PT).astype(int); t1 = np.floor((hi - lo0) / PT).astype(int)
    L = json.load(open(inp["libs"]))["libs"]
    PL = json.load(open(inp["preplaced_libs"])) if inp.has("preplaced_libs") else []
    has_pre = inp.has("preplaced_off") and len(PL) > 0
    poff = np.load(inp["preplaced_off"]) if has_pre else None
    recs = []
    for i in range(t0[0], t1[0] + 1):
        for j in range(t0[1], t1[1] + 1):
            for k in range(t0[2], t1[2] + 1):
                t = key.get((i, j, k))
                if t is None: continue
                f = os.path.join(inp["tiles_dir"], f"{t}.npz")
                if os.path.exists(f):
                    z = np.load(f); reg = np.repeat(z["job_region"].astype(np.int64), z["job_placed"])
                    inb = np.all((z["pos"] >= lo) & (z["pos"] < hi), axis=1)
                    for ty, p, q, r in zip(z["types"][inb], z["pos"][inb], z["quat"][inb], reg[inb]):
                        recs.append((L[ty], rn[r], int(r), p, q, "packed"))
                if has_pre and poff[t + 1] > poff[t]:
                    a, c = int(poff[t]), int(poff[t + 1])
                    pp = np.load(inp["preplaced_pos"], mmap_mode="r")[a:c]; pq = np.load(inp["preplaced_quat"], mmap_mode="r")[a:c]
                    pt = np.load(inp["preplaced_types"], mmap_mode="r")[a:c]
                    inb = np.all((pp >= lo) & (pp < hi), axis=1)
                    for ty, p, q in zip(pt[inb], pp[inb], pq[inb]):
                        recs.append((PL[ty], "prepass", -1, p, q, f"prepass:{PL[ty]}"))
    species = sorted({r[0] for r in recs})
    lc = common.lib_cache(inp, species, os.path.join(bdir, "validate", "cache", f"lib_box_{name}.npz"), nproc=4) if species else None
    kT = 1.380649e-23 * sd["temperature_K"]; eta = sd["viscosity_Pa_s"]
    rh = {}
    zc = np.load(os.path.join(bdir, "validate", "cache", f"lib_box_{name}.npz"), allow_pickle=False) if species else None
    for s in species:
        i = lc.idx[s]; nres = max(int(zc["n_beads"][i]), 1)
        rh[s] = 0.475 * nres ** 0.29                 # nm; Wilkins 1999 folded proteins, R_h(A) = 4.75 N^0.29
    sidx = {s: i for i, s in enumerate(species)}
    mem_regs = {r["name"] for r in regs if r["kind"] == "membrane"}
    rows = []
    for m, (s, rname, ri, p, q, src) in enumerate(recs):
        q = Q.canonical(np.asarray(q, float)[None])[0]
        R = rh[s] * 1e-9
        rows.append(dict(mol_id=m, species=s, species_idx=sidx[s], region=rname, region_id=ri, x_nm=float(p[0]), y_nm=float(p[1]),
                         z_nm=float(p[2]), qw=q[0], qx=q[1], qy=q[2], qz=q[3], R_h_nm=round(rh[s], 4),
                         D_t_um2_s=kT / (6 * np.pi * eta * R) * 1e12, D_r_per_us=kT / (8 * np.pi * eta * R ** 3) * 1e-6,
                         fixed=src != "packed", tethered=rname in mem_regs, src=src, polymer_id=-1, polymer_pos=-1))
    T = pd.DataFrame(rows)
    sb, sb_meta = secms_bound_columns(species, os.path.dirname(inp["pack_table"]))
    if len(T):
        T["secms_bound_fraction"] = T.species.map(lambda s: sb[s][0]).astype(float)
        T["secms_unrepresented_share_of_free"] = T.species.map(lambda s: sb[s][1]).astype(float)
    T.to_parquet(os.path.join(out, "molecules.parquet"))
    if species:
        off, C, Rr, br, _ = lc.sub(species)
        np.savez(os.path.join(out, "bodies.npz"), species=np.array(species), off=off, centres=C.astype(np.float32), radii=Rr.astype(np.float32),
                 level=np.array("sphere_tree"), rh=np.array([rh[s] for s in species]), br=br, mass_da=np.array([float(zc["bead_mass"][lc.idx[s]]) for s in species]))
    h = float(kw.get("wall_h_nm", 4.0))
    walls = {}
    try:
        geo = open_geometry(cfg, inp, bdir)
        g = np.arange(lo[0], hi[0] + 1e-9, h); gy = np.arange(lo[1], hi[1] + 1e-9, h); gz_ = np.arange(lo[2], hi[2] + 1e-9, h)
        P = np.stack(np.meshgrid(g, gy, gz_, indexing="ij"), -1).reshape(-1, 3)
        S = geo.mems(P)
        for mi, mn in enumerate(geo.mem_names):
            walls[f"sdf_{mn}"] = S[:, mi].reshape(len(g), len(gy), len(gz_)).astype(np.float32)
        hl = 10.0
        gl = [np.arange(lo[a] + hl / 2, hi[a], hl) for a in range(3)]
        PL_ = np.stack(np.meshgrid(*gl, indexing="ij"), -1).reshape(-1, 3)
        walls["label_grid"] = geo.labels(PL_).reshape(*[len(x) for x in gl]).astype(np.uint8)
        walls["label_names"] = np.array(geo.label_names)
        walls["label_h_nm"] = np.array(hl)
        geo_src = geo.source
    except Exception as e:  # noqa: BLE001 - geometry unavailable: box without walls
        geo_src = f"unavailable ({type(e).__name__}: {e})"
    if walls:
        np.savez(os.path.join(out, "walls.npz"), origin_nm=lo, h_nm=np.array(h), **walls)
    box = dict(schema="wholecell.simbox/1", build_id=_build_id(bdir), cell_id=cfg["cell_id"], stub=True,
               box=dict(lo_nm=lo.tolist(), hi_nm=hi.tolist(), boundary="walls" if walls else "periodic", chosen_by=how),
               frame="model frame of the build geometry (nm)", units=dict(length="nm", time="ns", D_t="um2/s", D_r="1/us"),
               quaternion_order="w,x,y,z", **sd,
               diffusion_model=dict(name="Stokes-Einstein, R_h = 0.475 nm x N_res^0.29 (Wilkins 1999)", fixed_inputs=None),
               counts=dict(molecules=len(T), packed=int((T.src == "packed").sum()) if len(T) else 0,
                           prepass=int((T.src != "packed").sum()) if len(T) else 0, species=len(species)),
               geometry=geo_src, not_written=["beads.cif", "polymers.npz"], secms_bound_fraction=sb_meta,
               exporter_commit=__import__("wholecell.build", fromlist=["git_info"]).git_info()["commit7"])
    json.dump(box, open(os.path.join(out, "box.json"), "w"), indent=1)
    _manifest(out)
    print(json.dumps(box["counts"]), how)
    return box


def run(cfg, bdir, step=None, **kw):
    step = step or "all"
    inp = common.Inputs(cfg, bdir, **kw)
    if step == "all":
        for s in ("regions", "cell_entry", "sim_manifest"):
            run(cfg, bdir, s, **kw)
        if ((cfg.get("export", {}) or {}).get("tiles", {}) or {}).get("enabled"):
            run(cfg, bdir, "tiles_v2", **kw)
        return
    fn = dict(regions=step_regions, cell_entry=step_cell_entry, tiles_v2=step_tiles_v2, web=step_tiles_v2,
              sim_manifest=step_sim_manifest).get(step)
    if fn is None:
        raise ValueError(f"unknown export step '{step}' (regions, cell_entry, tiles_v2 (= web), sim_manifest, all)")
    return fn(cfg, bdir, inp, kw)
