"""Shared helpers of the validation stage: input resolution, library cache, whole-cell tile scan, small IO helpers.

Input resolution (`Inputs`): every validation step reads the pack outputs of ONE build. By default these are the files
of `<build>/pack/` in the pack-stage names (tiles/<i>.npz, tiles.parquet, tiles_meta.json, counts.npz, libs.json,
preplaced_{types,pos,quat,off}.npy, preplaced_libs.json, mrna_{pts,off}.npy, mrna_meta.json, summary.json), the inventory of
`<build>/inventory/`, the geometry package of the build and its structure library. Every one of them can be overridden:
  pack_dir=<dir>                 another pack folder with the same file names
  tiles_dir=<dir>                only the tile folder
  ref=<name>                     the file map of config `validate.compare_reference` (e.g. the archived legacy v7 pack,
                                 whose files carry _v7 suffixes); used to run the same code on a reference pack
  <key>=<path>                   any single key of PACK_FILES / INV_FILES, plus geom=, lib_dir=, regions_json=
No path is hard-coded here: defaults are build-relative, everything else comes from config or arguments.
Builds made before the stage rename (A0c: pack/prepass_*.npy, prepass/prepass_manifest.json) are read through
wholecell.build.ref_path (old file names when the new ones do not exist).
"""
from __future__ import annotations

import glob
import hashlib
import json
import os
from concurrent.futures import ProcessPoolExecutor

import numpy as np

from ..build import ref_path

PACK_FILES = {
    "tiles_dir": "tiles", "tiles_parquet": "tiles.parquet", "tiles_meta": "tiles_meta.json", "counts": "counts.npz",
    "libs": "libs.json", "preplaced_types": "preplaced_types.npy", "preplaced_pos": "preplaced_pos.npy",
    "preplaced_quat": "preplaced_quat.npy", "preplaced_off": "preplaced_off.npy", "preplaced_libs": "preplaced_libs.json",
    "preplaced_groups": "preplaced_groups.json",
    "mrna_pts": "mrna_pts.npy", "mrna_off": "mrna_off.npy", "mrna_meta": "mrna_meta.json", "summary": "summary.json",
}
INV_FILES = {"pack_table": "pack_table.tsv", "inventory_v0": "v0/inventory_v0.tsv"}   # inventory stage layout (stage.steps: v0/inventory_v0.tsv)
MAX_PROCS = 8          # shared node: never more than 8 worker processes
CACHE_DIR = None       # set by the common suite (A8 phase 2) when it runs core steps on a build it must not write into


def cache_path(bdir, name):
    """validation cache file: <build>/validate/cache/<name>, or CACHE_DIR/<name> when set (suite --out runs); a cache
    written before A0c under its old name (lib_prepass.npz) is reused"""
    return ref_path(CACHE_DIR or os.path.join(bdir, "validate", "cache"), name)


def lib_masses(inp, bdir, names=None, nproc=MAX_PROCS):
    """per-library mass (Da) of the packed library entries: libs.json mass_da when present (Jurkat, testcell), else the
    library bead mass (aa 110 Da, nt 330 Da: the densities rule for pre-pass entries; yeast pack/libs.json has no
    mass_da). Returns (mass array, source text)."""
    L = json.load(open(inp["libs"]))
    names = names or L["libs"]
    m = L.get("mass_da") or []
    if len(m) == len(names):
        return np.asarray(m, float), "libs.json mass_da"
    lc = lib_cache(inp, list(names), cache_path(bdir, "lib_packed.npz"), nproc=nproc)
    return np.array([lc.bead_mass[lc.idx[n]] for n in names], float), "library bead mass (aa 110 Da, nt 330 Da)"


def nprocs(cfg=None, kw=None) -> int:
    n = int((kw or {}).get("nproc") or (cfg or {}).get("validate", {}).get("nproc", MAX_PROCS))
    return max(1, min(n, MAX_PROCS))


class Inputs:
    """resolved input paths of one validation run (see module docstring)"""

    def __init__(self, cfg, bdir, **kw):
        self.cfg, self.bdir, self.kw = cfg, bdir, kw
        vcfg = cfg.get("validate", {}) or {}
        ref = kw.get("ref")
        refmap = {}
        if ref:
            cr = vcfg.get("compare_reference", {}) or {}
            if cr.get("name", ref) != ref and ref not in (cr.get("aliases") or []):
                raise KeyError(f"reference '{ref}' is not configured in validate.compare_reference")
            refmap = dict(cr.get("files", {}) or {})
        pack_dir = kw.get("pack_dir") or os.path.join(bdir, "pack")
        inv_dir = kw.get("inventory_dir") or os.path.join(bdir, "inventory")
        self.paths = {}
        for k, fn in PACK_FILES.items():
            self.paths[k] = kw.get(k) or refmap.get(k) or ref_path(pack_dir, fn)
        for k, fn in INV_FILES.items():
            self.paths[k] = kw.get(k) or refmap.get(k) or ref_path(inv_dir, fn)
        for k in ("geom", "lib_dir", "regions_json", "polysomes", "geometry_validation", "pack_trees"):
            v = kw.get(k) or refmap.get(k)
            if v:
                self.paths[k] = v
        for k, v in refmap.items():                     # extra reference-only keys (mrna sources etc.)
            self.paths.setdefault(k, v)
        self.ref = ref
        self.vcfg = vcfg

    def __getitem__(self, k):
        return self.paths[k]

    def get(self, k, default=None):
        return self.paths.get(k, default)

    def has(self, k):
        p = self.paths.get(k)
        return bool(p) and os.path.exists(p)

    def tile_files(self):
        d = self.paths["tiles_dir"]
        fs = glob.glob(os.path.join(d, "*.npz"))
        return sorted(fs, key=lambda f: int(os.path.basename(f)[:-4]))

    def n_tiles(self):
        import pandas as pd
        return len(pd.read_parquet(self.paths["tiles_parquet"], columns=["tile"]))

    def lib_dir(self):
        """structure library folder: lib_dir= > validate.lib_dir > library.dir > frozen input lib_structures
        ({build} in validate.lib_dir = this build's folder, e.g. yeast "{build}/engine/lib": the sealed library plus the
        entries the assembly hooks and pre-placed models add)"""
        p = self.paths.get("lib_dir") or self.vcfg.get("lib_dir") or (self.cfg.get("library", {}) or {}).get("dir")
        if p:
            p = str(p).replace("{build}", self.bdir)
        fi = self.cfg.get("frozen_inputs", {}) or {}
        fr = (self.cfg.get("library", {}) or {}).get("frozen")
        if not p and (self.cfg.get("library", {}) or {}).get("mode") == "build":    # built per build (yeast; Jurkat B3)
            p = os.path.join(self.bdir, "library", "lib")
        if not p and fr and fr in fi:
            p = fi[fr]["path"]
        if not p:
            for name in ("lib_structures_v7", "lib_structures", "lib_library", "lib_dir"):
                if name in fi:
                    p = fi[name]["path"]
                    break
        if not p:
            p = os.path.join(self.bdir, "library", "lib")
        return p

    def lookup_input(self, name):
        """path of a pack/preplaced input used by this build (polysomes, mRNA sources, ...):
        kw/reference map > <build>/preplaced/preplaced_manifest.json (entries of 'instances' / 'beads' / 'drops' by
        name, their 'npz') > validate.inputs.<name> or a frozen-input name. A build with a preplaced manifest that does
        not name the input falls back with a warning (a silent fallback once validated legacy v7 files)."""
        if self.paths.get(name):
            return self.paths[name]
        fi = self.cfg.get("frozen_inputs", {}) or {}
        man = ref_path(self.bdir, os.path.join("preplaced", "preplaced_manifest.json"))
        has_man = os.path.exists(man)
        if has_man:
            p = manifest_npz(json.load(open(man)), name)
            if p:
                return p
        if name in fi and name not in (self.vcfg.get("inputs", {}) or {}):    # a frozen input named directly (EM, tables)
            return fi[name]["path"]
        v = (self.vcfg.get("inputs", {}) or {}).get(name)
        if v and v in fi:
            v = fi[v]["path"]
        if v and has_man:
            print(f"validate: WARNING input '{name}' is not in {man}; using validate.inputs {v}", flush=True)
        return v


def manifest_npz(M, name):
    """npz of a named entry of a preplaced_manifest.json: the current format (lists 'instances', 'beads', 'drops' of
    {name, npz, ...}) or the older dict forms ({name: {path}}, {name: path}, {'modules': {...}})"""
    for sec in ("instances", "beads", "drops"):
        for e in M.get(sec, []) or []:
            if isinstance(e, dict) and e.get("name") == name and (e.get("npz") or e.get("path")):
                return e.get("npz") or e.get("path")
    ent = M.get(name)
    if ent is None and isinstance(M.get("modules"), dict):
        ent = M["modules"].get(name)
    if isinstance(ent, dict):
        return ent.get("npz") or ent.get("path")
    if isinstance(ent, str):
        return ent
    return None


# ----------------------------------------------------------------- IO helpers
def md5(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def outdir(bdir, *sub):
    d = os.path.join(bdir, "validate", *sub)
    os.makedirs(d, exist_ok=True)
    return d


def _jsonable(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, set):
        return sorted(o)
    return str(o)


def write_json(path, obj, overwrite=False):
    """atomic JSON write; refuses to overwrite unless overwrite=True (stage outputs are write-once)"""
    if os.path.exists(path) and not overwrite:
        raise FileExistsError(f"{path} exists (validation outputs are write-once; pass overwrite=1 to replace)")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=_jsonable)
    os.replace(tmp, path)
    return path


def write_df(path, df, overwrite=False, **kw):
    if os.path.exists(path) and not overwrite:
        raise FileExistsError(f"{path} exists (validation outputs are write-once; pass overwrite=1 to replace)")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = f"{path}.tmp{os.getpid()}"
    if path.endswith(".parquet"):
        df.to_parquet(tmp)
    else:
        df.to_csv(tmp, sep="\t", **kw)
    os.replace(tmp, path)
    return path


def flag(kw, name, default=False):
    v = kw.get(name, default)
    if isinstance(v, str):
        return v.lower() in ("1", "true", "yes", "y")
    return bool(v)


# ----------------------------------------------------------------- structure library cache
def _lib_one(args):
    lib_dir, name = args
    f = os.path.join(lib_dir, f"{name}.npz")
    if not os.path.exists(f):
        return name, None
    z = np.load(f)
    stc = z["st_c"].astype(np.float32) if "st_c" in z.files else np.zeros((0, 3), np.float32)
    strr = z["st_r"].astype(np.float32) if "st_r" in z.files else np.zeros(0, np.float32)
    kind = z["kind"] if "kind" in z.files else np.zeros(len(z["beads"]), np.int8)
    mass_beads = float(np.where(kind == 0, 110.0, 330.0).sum())        # legacy densities rule: aa 110 Da, nt 330 Da
    memb = str(z["membrane"]) if "membrane" in z.files else ""
    br = float((np.linalg.norm(stc, axis=1) + strr).max()) if len(strr) else 0.0
    return name, (stc, strr, mass_beads, memb, br, len(z["beads"]))


def lib_cache(inp: Inputs, names, cache_path, nproc=MAX_PROCS):
    """sphere trees, bead mass, membrane tag and bounding radius for the given library names, cached in one npz
    (the library is 16k small files on a slow file system: read once per validation build)"""
    names = list(dict.fromkeys(names))
    if os.path.exists(cache_path):
        z = np.load(cache_path, allow_pickle=False)
        have = list(z["names"])
        if set(names) <= set(have):
            return _LibCache(z)
    lib_dir = inp.lib_dir()
    with ProcessPoolExecutor(nproc) as ex:
        res = list(ex.map(_lib_one, [(lib_dir, n) for n in names], chunksize=64))
    off = [0]; C = []; R = []; mass = []; memb = []; br = []; nb = []; ok = []
    for name, r in res:
        if r is None:
            r = (np.zeros((0, 3), np.float32), np.zeros(0, np.float32), 0.0, "", 0.0, 0)
            ok.append(False)
        else:
            ok.append(True)
        C.append(r[0]); R.append(r[1]); off.append(off[-1] + len(r[1])); mass.append(r[2]); memb.append(r[3])
        br.append(r[4]); nb.append(r[5])
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    tmp = cache_path + f".tmp{os.getpid()}.npz"
    np.savez(tmp, names=np.array(names), st_off=np.array(off, np.int64), stc=np.concatenate(C) if C else np.zeros((0, 3)),
             str=np.concatenate(R) if R else np.zeros(0), bead_mass=np.array(mass), membrane=np.array(memb),
             br=np.array(br), n_beads=np.array(nb), present=np.array(ok), lib_dir=np.array(lib_dir))
    os.replace(tmp, cache_path)
    return _LibCache(np.load(cache_path, allow_pickle=False))


class _LibCache:
    def __init__(self, z):
        self.names = [str(x) for x in z["names"]]
        self.idx = {n: i for i, n in enumerate(self.names)}
        self.st_off = z["st_off"]; self.stc = z["stc"].astype(np.float64); self.str = z["str"].astype(np.float64)
        self.bead_mass = z["bead_mass"]; self.membrane = [str(x) for x in z["membrane"]]; self.br = z["br"]
        self.present = z["present"]

    def sub(self, names):
        """arrays re-indexed to `names` order: st_off, stc, str, br, membrane flags"""
        ii = [self.idx[n] for n in names]
        off = [0]; C = []; R = []
        for i in ii:
            a, b = self.st_off[i], self.st_off[i + 1]
            C.append(self.stc[a:b]); R.append(self.str[a:b]); off.append(off[-1] + b - a)
        return (np.array(off, np.int64), np.concatenate(C) if C else np.zeros((0, 3)), np.concatenate(R) if R else np.zeros(0),
                self.br[ii].astype(np.float64), np.array([self.membrane[i] != "" for i in ii]))


# ----------------------------------------------------------------- whole-cell tile scan (job tables only)
def _scan_chunk(args):
    files, nreg = args
    acc_req = {}; acc_pla = {}; ninst = []
    for f in files:
        z = np.load(f)
        jl, jr, jn, jp = z["job_lib"].astype(np.int64), z["job_region"].astype(np.int64), z["job_n"], z["job_placed"]
        k = jl * nreg + jr
        u, inv = np.unique(k, return_inverse=True)
        n = np.bincount(inv, jn); p = np.bincount(inv, jp)
        for a, b, c in zip(u.tolist(), n.tolist(), p.tolist()):
            acc_req[a] = acc_req.get(a, 0) + b; acc_pla[a] = acc_pla.get(a, 0) + c
        ninst.append((int(os.path.basename(f)[:-4]), int(jp.sum())))
    return acc_req, acc_pla, ninst


def tile_sums(inp: Inputs, cache_path=None, nproc=MAX_PROCS, nreg=None, files=None):
    """per (library, region) requested (job_n) and placed (job_placed) instances over all tiles, plus placed
    instances per tile. Reads only the job_* keys of each tile npz. Cached (npz) when cache_path is given."""
    if cache_path and os.path.exists(cache_path):
        z = np.load(cache_path)
        return dict(req=z["req"], pla=z["pla"], tile=z["tile"], ninst=z["ninst"], nreg=int(z["nreg"]))
    libs = json.load(open(inp["libs"]))
    nl = len(libs["libs"])
    nreg = nreg or 64
    files = files if files is not None else inp.tile_files()
    chunks = [(files[i::64], nreg) for i in range(min(64, len(files)))]
    req = np.zeros((nl, nreg)); pla = np.zeros((nl, nreg)); tiles = []
    with ProcessPoolExecutor(nproc) as ex:
        for acc_req, acc_pla, ninst in ex.map(_scan_chunk, chunks):
            for key, v in acc_req.items():
                req[key // nreg, key % nreg] += v
            for key, v in acc_pla.items():
                pla[key // nreg, key % nreg] += v
            tiles += ninst
    tiles.sort()
    out = dict(req=req, pla=pla, tile=np.array([t for t, _ in tiles], np.int64),
               ninst=np.array([n for _, n in tiles], np.int64), nreg=nreg)
    if cache_path:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        tmp = cache_path + f".tmp{os.getpid()}.npz"
        np.savez(tmp, **out)
        os.replace(tmp, cache_path)
    return out


def sample_tiles(cfg, n_tiles, n, stream="validate.placement.tiles", seed=None):
    """seeded tile sample (same ids for every build of one config): wholecell.seeds stream, or an explicit seed"""
    from .. import seeds
    rng = np.random.default_rng(seed) if seed is not None else seeds.rng(cfg, stream)
    return np.sort(rng.choice(n_tiles, min(n, n_tiles), replace=False))
