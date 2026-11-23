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
from . import quat as Q, yct2  # noqa: E402

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
    poff = np.load(inp["preplaced_off"]) if has_pre else None
    ptypes = np.load(inp["preplaced_types"], mmap_mode="r") if has_pre else None
    ppos = np.load(inp["preplaced_pos"], mmap_mode="r") if has_pre else None
    pquat = np.load(inp["preplaced_quat"], mmap_mode="r") if has_pre else None
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
    cache = _TileCache(inp)
    ppb = int(round(B / PT))
    tiles_rows = []; counts = {}; files = []
    src_n = 0
    check = dict(max_pos_err_nm=0.0, max_rot_err_rad=0.0, n_checked=0)
    tdir = os.path.join(web, "tiles"); os.makedirs(tdir, exist_ok=True)
    for b in [tuple(x) for x in blocks.tolist()]:
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
                        parts.append((v[0], v[3], np.zeros(len(v[0]), np.uint32), v[1], v[2], 0 <= di < ppb and 0 <= dj < ppb and 0 <= dk < ppb))
                    if has_pre and poff[t + 1] > poff[t]:
                        a, c = int(poff[t]), int(poff[t + 1]); pt = np.asarray(ptypes[a:c]).astype(np.int64)
                        parts.append((pmap[pt], pre_reg[pt], np.full(c - a, yct2.FLAG_PREPASS, np.uint32),
                                      np.asarray(ppos[a:c], np.float64), np.asarray(pquat[a:c], np.float64), False))
        if not parts:
            continue
        sp = np.concatenate([p[0] for p in parts]); rg = np.concatenate([p[1] for p in parts]); fl = np.concatenate([p[2] for p in parts])
        ps = np.concatenate([p[3] for p in parts]); qu = np.concatenate([p[4] for p in parts])
        ijk = np.floor((ps - origin) / tile_nm).astype(np.int64)
        inb = np.all((ijk >= np.array(b) * bt) & (ijk < (np.array(b) + 1) * bt), axis=1)
        sp, rg, fl, ps, qu, ijk = sp[inb], rg[inb], fl[inb], ps[inb], qu[inb], ijk[inb]
        src_n += len(sp)
        o = np.lexsort((ijk[:, 2], ijk[:, 1], ijk[:, 0]))
        sp, rg, fl, ps, qu, ijk = sp[o], rg[o], fl[o], ps[o], qu[o], ijk[o]
        brk = np.r_[0, np.nonzero(np.any(np.diff(ijk, axis=0) != 0, axis=1))[0] + 1, len(sp)]
        shard_k = 0; buf = []; size = 0; members = []
        name = lambda k: f"s_{b[0]}_{b[1]}_{b[2]}.{k}.bin"
        for a, c in zip(brk[:-1], brk[1:]):
            if c <= a: continue
            t_ijk = ijk[a].tolist(); lo = origin + np.array(t_ijk) * tile_nm
            raw, n_i, n_p = yct2.encode_tile(t_ijk, lo, sp[a:c], rg[a:c], fl[a:c], ps[a:c], qu[a:c], tile_nm, q_org, q_range, lib_tag, build_tag)
            # round trip on every tile
            D = yct2.decode_tile(raw)
            pos_d = origin + D["lo_rel"] + D["pos_rel"]
            o2 = np.lexsort((rg[a:c], sp[a:c]))
            check["max_pos_err_nm"] = max(check["max_pos_err_nm"], float(np.abs(pos_d - ps[a:c][o2]).max()))
            check["max_rot_err_rad"] = max(check["max_rot_err_rad"], float(Q.angle_between(D["quat"], qu[a:c][o2]).max()))
            check["n_checked"] += int(c - a)
            for s_, r_ in zip(D["species"].tolist(), D["region"].tolist()):
                counts[(s_, r_)] = counts.get((s_, r_), 0) + 1
            g = yct2.gz(raw)
            if buf and size + len(g) > shard_max:
                members.append((shard_k, buf)); shard_k += 1; buf = []; size = 0
            buf.append((t_ijk, g, n_i, yct2.crc32(raw))); size += len(g)
        if buf:
            members.append((shard_k, buf))
        single = len(members) == 1
        for k, mem in members:
            fn = f"s_{b[0]}_{b[1]}_{b[2]}.bin" if single else name(k)
            off = 0
            with open(os.path.join(tdir, fn), "wb") as fh:
                for t_ijk, g, n_i, crc in mem:
                    fh.write(g)
                    tiles_rows.append([*t_ijk, fn, off, len(g), n_i, 0, 0, 0, 0, 0, crc])
                    off += len(g)
            files.append(fn)
    allt = np.array([r[:3] for r in tiles_rows]) if tiles_rows else np.zeros((0, 3), int)
    index = dict(format="YCT2", version="2.0", build_id=bid, tile_nm=tile_nm, block_tiles=bt,
                 grid=dict(origin_nm=origin.tolist(), min_ijk=allt.min(0).tolist() if len(allt) else None,
                           max_ijk=allt.max(0).tolist() if len(allt) else None),
                 q_org=q_org, q_range=q_range, lod=[0], flag_bits=yct2.CORE_FLAG_BITS + list(ec.get("flag_bits", []) or []),
                 columns=[], polyline_kinds={}, sections_written=["PALETTE", "INST"],
                 not_written=["GHOST", "POLY", "MESH", "COLUMN", "LOD1 overview"], tile_columns=yct2.TILE_COLUMNS,
                 lib_tag=lib_tag, build_tag=build_tag, tiles=tiles_rows)
    ij = json.dumps(index, separators=(",", ":"))
    open(os.path.join(tdir, "index.json"), "w").write(ij)
    open(os.path.join(tdir, "index.json.gz"), "wb").write(yct2.gz(ij.encode()))
    rr = dict(schema="wholecell.regions/1", cell_id=cfg["cell_id"], build_id=bid, regions=rows)
    json.dump(rr, open(os.path.join(web, "regions.json"), "w"), indent=1)
    os.makedirs(os.path.join(web, "library"), exist_ok=True)
    lc = None
    try:
        lc = common.lib_cache(inp, types, os.path.join(bdir, "validate", "cache", "lib_all.npz"), nproc=common.nprocs(cfg, kw))
    except Exception as e:  # noqa: BLE001 - library absent: types without radii
        print("library cache unavailable:", e)
    # mass per type: libs.json mass_da (Jurkat, testcell); else, and for pre-pass types, the library bead mass (aa 110 Da,
    # nt 330 Da), as the yeast engine computed it (webenc.encode_library, val3d.lib_beads) and validate.common.lib_masses
    tj = [dict(id=i, name=n, mass_da=(mass[i] if i < len(mass) else (float(lc.bead_mass[lc.idx[n]]) if lc else None)),
               br_nm=(round(float(lc.br[lc.idx[n]]), 3) if lc else None), prepass=n in PL and n not in libs) for i, n in enumerate(types)]
    json.dump(dict(schema="wholecell.types/1", ids=types, types=tj), open(os.path.join(web, "library", "types.json"), "w"))
    json.dump(dict(v=bid), open(os.path.join(web, "version.json"), "w"))
    return index, counts, check, src_n, files


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
    web = _out(bdir, kw, "web")
    if os.listdir(web) and not common.flag(kw, "overwrite"):
        raise FileExistsError(f"{web} is not empty (export outputs are write-once; overwrite=1 to replace)")
    if os.listdir(web):                                   # overwrite: start from an empty folder (no stale shards)
        shutil.rmtree(web); os.makedirs(web)
    t0 = time.time()
    index, counts, check, src_n, files = _encode(cfg, bdir, inp, kw, web)
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
                inputs=dict(tiles_dir=inp["tiles_dir"], libs=inp["libs"]))
    json.dump(meta, open(os.path.join(web, "encode_meta.json"), "w"), indent=1)
    digest = _manifest(web)
    res = dict(encode_check=ec, manifest_digest_excl_meta=digest)
    if common.flag(kw, "selftest"):
        tmp = tempfile.mkdtemp(prefix="yct2_selftest_", dir=kw.get("tmp_dir") or None)
        try:
            kw2 = dict(kw, out_dir=tmp, selftest="0", overwrite="1")
            w2 = os.path.join(tmp, "web"); os.makedirs(w2, exist_ok=True)
            _encode(cfg, bdir, inp, kw2, w2)
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
    fn = dict(regions=step_regions, cell_entry=step_cell_entry, tiles_v2=step_tiles_v2, sim_manifest=step_sim_manifest).get(step)
    if fn is None:
        raise ValueError(f"unknown export step '{step}' (regions, cell_entry, tiles_v2, sim_manifest, all)")
    return fn(cfg, bdir, inp, kw)
