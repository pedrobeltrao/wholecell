"""cellgeom_v1 packages: write, load, validate (spec: docs/reference/cellgeom_v1.md).

    from wholecell.geometry import cellgeom
    w = cellgeom.PackageWriter(dest)                 # dest = {geom_root}/<cell>/<geometry_id> (must not exist)
    w.add_field("cell", h10, index, data); w.add_section(...); w.set_pores(pores, source); w.add_layer(...)
    man = w.finalize(manifest)                       # sha256 per file, manifest.json, validation.json, read-only files
    cell = cellgeom.load(dest)                       # verifies checksums, builds G for evaluator cellgeom_eval_v1
    cell = cellgeom.load_cfg(cfg)                    # the same with config geometry.package / .evaluator / .rules
    res = cellgeom.validate(cell)                    # dict of checks (also written by finalize)

The loaded object is a subclass of the engine's Cell: .G, .pores, .pidx/.PT/.PQ/.PC, .labels(P), .mems(P), plus
.manifest, .vocab (label id -> entry), .summary (legacy-shaped per-field area/volume dict + open_face f_vol/f_pm),
.package_dir, .package_id, .engine (evaluator module), .rules / .rules_path / .rules_sha256 (rule table, spec §4.3),
.code_to_name / .name_to_code (engine label codes), .mem_code(name).
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import shutil
import stat

import numpy as np

SPEC_VERSION = "cellgeom_v1"

# shared core label vocabulary (ids < 1000; stable). Cell-specific extensions use ids >= 1000 with a parent core id.
CORE_LABELS = [
    (0, "exterior", "GO:0005576", "outside the plasma membrane"),
    (1, "cytosol", "GO:0005829", "cytosol (incl. cytoskeleton space)"),
    (2, "nucleoplasm", "GO:0005654", ""),
    (3, "nucleolus", "GO:0005730", ""),
    (4, "perinuclear_space", "GO:0005641", "nuclear envelope lumen (between ONM and INM)"),
    (5, "er_lumen", "GO:0005788", ""),
    (6, "nuclear_pore_channel", "GO:0044613", "channel of the nuclear pore (NPC central channel space)"),
    (7, "mito_ims", "GO:0005758", "mitochondrial intermembrane space (incl. crista lumen when not separate)"),
    (8, "outside_imaged_volume", None, "beyond a crop / truncation face of the source volume (not modelled)"),
    (9, "mito_matrix", "GO:0005759", ""),
    (10, "golgi_lumen", "GO:0005796", ""),
    (11, "endolysosome_lumen", None, "endosome / lysosome lumen (EM endo+lyso classes merged)"),
    (12, "vacuole_lumen", "GO:0005775", ""),
    (13, "peroxisome_matrix", "GO:0005782", ""),
    (14, "lipid_droplet", "GO:0005811", ""),
    (15, "cell_wall", "GO:0005618", ""),
    (16, "nuclear_speck", "GO:0016607", ""),
    (17, "golgi_intercisternal", None, "cytosol between cisternae inside a Golgi stack envelope"),
    (18, "crista_lumen", "GO:0061617", "intracristal space when modelled separately from the IMS"),
    (19, "lysosome_lumen", "GO:0043202", ""),
]
CORE_BY_NAME = {n: i for i, n, _, _ in CORE_LABELS}
EXT_ID_MIN = 1000
MAX_ENGINE_CODE = 63          # packer label bitmasks are int64


def sha256_file(path, chunk=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def md5_file(path, chunk=1 << 24):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def savez_det(path, **arrays):
    """np.savez with fixed zip timestamps (byte-identical output for identical arrays -> stable sha256)"""
    import zipfile
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED, allowZip64=True) as zf:
        for k, a in arrays.items():
            zi = zipfile.ZipInfo(k + ".npy", date_time=(1980, 1, 1, 0, 0, 0))
            with zf.open(zi, "w", force_zip64=True) as f:
                np.lib.format.write_array(f, np.asanyarray(a), allow_pickle=False)


def core_labels_json():
    return [{"id": i, "name": n, "go": g, "description": d} for i, n, g, d in CORE_LABELS]


# ----------------------------------------------------------------- field files (sparse-brick SDF, one npz per field)
def write_field_npz(path, h10, index, data, dtype=None):
    """header float64[12]: origin xyz, h, nbx, nby, nbz, B, clamp_in, clamp_out, 0, len(index) (offsets of the index
    and data blocks relative to a concatenated [index, data] array); index int32 (-1 inside, -2 outside, k = brick);
    data (n_stored*(B+1)^3,) float32 when lossless else float64."""
    h10 = np.asarray(h10, np.float64); index = np.asarray(index); data = np.asarray(data, np.float64)
    if dtype is None:
        dtype = np.float32 if np.array_equal(data.astype(np.float32).astype(np.float64), data) else np.float64
    idx = index.astype(np.int32)
    assert np.array_equal(idx, index), "brick index not integer"
    savez_det(path, header=np.r_[h10, 0.0, float(len(idx))], index=idx, data=data.astype(dtype))
    return np.dtype(dtype).name


def read_field_npz(path):
    z = np.load(path)
    return z["header"][:10].astype(np.float64), z["index"].astype(np.float64), z["data"].astype(np.float64)


# ----------------------------------------------------------------- writer
class PackageWriter:
    """build a package in <dest>.tmp-<pid>; finalize() writes manifest + validation and renames to dest (immutable)"""

    def __init__(self, dest):
        self.dest = os.path.abspath(dest)
        if os.path.exists(self.dest):
            raise FileExistsError(f"{self.dest} exists: packages are immutable, use a new geometry_id")
        self.tmp = f"{self.dest}.tmp-{os.getpid()}"
        if os.path.exists(self.tmp):
            shutil.rmtree(self.tmp)
        os.makedirs(self.tmp)
        self.fields, self.sections, self.layers, self.files_meta = {}, {}, {}, {}
        self.pores_meta = None

    def path(self, rel):
        p = os.path.join(self.tmp, rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        return p

    def add_field(self, name, slot, h10, index, data, **meta):
        rel = f"fields/{name}.npz"
        dt = write_field_npz(self.path(rel), h10, index, data)
        n_st = int((np.asarray(index) >= 0).sum())
        self.fields[name] = dict(slot=int(slot), file=rel, dtype=dt, h_nm=float(h10[3]), brick=int(h10[7]),
                                 origin_nm=[float(x) for x in h10[:3]], bricks=[int(x) for x in h10[4:7]],
                                 clamp_in_nm=float(h10[8]), clamp_out_nm=float(h10[9]), bricks_stored=n_st,
                                 sign="negative inside", **meta)

    def add_empty_field(self, name, slot, **meta):
        self.fields[name] = dict(slot=int(slot), file=None, empty=True, sign="negative inside",
                                 note="absent: evaluates to +clamp everywhere (outside everywhere)", **meta)

    def add_section(self, name, type_id, arrays: dict, params: dict | None = None, **meta):
        rel = f"sections/{name}.npz"
        savez_det(self.path(rel), **arrays)
        self.sections[name] = dict(type=type_id, file=rel, params=params or {}, **meta)

    def set_pores(self, pores, source=None, source_names=("random", "em"), **meta):
        pores = np.asarray(pores, np.float64).reshape(-1, 6)
        src = np.zeros(len(pores), np.int8) if source is None else np.asarray(source, np.int8)
        savez_det(self.path("pores.npz"), pores=pores, source=src)
        self.pores_meta = dict(file="pores.npz", n=int(len(pores)), columns=["cx", "cy", "cz", "nx", "ny", "nz"],
                               source_codes=list(source_names),
                               counts={n: int((src == k).sum()) for k, n in enumerate(source_names)}, **meta)

    def add_file(self, rel, src=None, data=None, text=None):
        p = self.path(rel)
        if src is not None:
            shutil.copyfile(src, p)
        elif data is not None:
            savez_det(p, **data)
        else:
            open(p, "w").write(text)
        return rel

    def add_layer(self, name, kind, schema_version, files=(), refs=(), provenance=None, depends_on=None, **meta):
        """files: package-relative paths already added; refs: [{path, md5, size, role}] external files (not copied)"""
        self.layers[name] = dict(kind=kind, schema_version=int(schema_version), depends_on=depends_on,
                                 files=list(files), refs=list(refs), provenance=provenance or {}, **meta)

    def finalize(self, manifest: dict, validate_kw: dict | None = None, extra_checks=None, publish_failed=False):
        """extra_checks: callable(cell) -> {name: result dict with 'pass'} merged into validation.json (e.g. B1).
        A package that fails validation is left in the tmp folder (not published) unless publish_failed."""
        man = dict(manifest)
        man.setdefault("spec_version", SPEC_VERSION)
        pid = f"{man['cell_id']}/{man['geometry_id']}"
        man["package_id"] = pid
        man.setdefault("created", _dt.datetime.now().isoformat(timespec="seconds"))
        man["fields"] = self.fields
        man["sections"] = self.sections
        man["pores"] = self.pores_meta
        for L in self.layers.values():
            if L["depends_on"] is None:
                L["depends_on"] = pid
        man["layers"] = self.layers
        files = {}
        for root, _, fs in os.walk(self.tmp):
            for f in sorted(fs):
                p = os.path.join(root, f); rel = os.path.relpath(p, self.tmp)
                files[rel] = {"sha256": sha256_file(p), "bytes": os.path.getsize(p)}
        man["files"] = dict(sorted(files.items()))
        with open(os.path.join(self.tmp, "manifest.json"), "w") as f:
            json.dump(man, f, indent=1, default=float)
        cell = load(self.tmp, verify=True)
        res = validate(cell, **(validate_kw or {}))
        if extra_checks is not None:
            res["checks"].update(extra_checks(cell))
            finish(res, cell.manifest)
        res["manifest_sha256"] = sha256_file(os.path.join(self.tmp, "manifest.json"))
        with open(os.path.join(self.tmp, "validation.json"), "w") as f:
            json.dump(res, f, indent=1, default=float)
        del cell
        if not res["pass"] and not publish_failed:
            raise RuntimeError(f"validation failed; package left unpublished in {self.tmp} (see validation.json)")
        os.rename(self.tmp, self.dest)
        for root, _, fs in os.walk(self.dest):
            for f in fs:
                p = os.path.join(root, f)
                os.chmod(p, os.stat(p).st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
        return man, res


# ----------------------------------------------------------------- loader
def verify_checksums(pkg, man=None, quick=False):
    man = man or json.load(open(os.path.join(pkg, "manifest.json")))
    bad = []
    for rel, meta in man["files"].items():
        p = os.path.join(pkg, rel)
        if not os.path.exists(p):
            bad.append((rel, "missing")); continue
        if os.path.getsize(p) != meta["bytes"]:
            bad.append((rel, "size")); continue
        if not quick and sha256_file(p) != meta["sha256"]:
            bad.append((rel, "sha256"))
    if bad:
        raise RuntimeError(f"cellgeom package {pkg}: checksum failures {bad}")
    return len(man["files"])


def _engine_module(name):
    from .engine import get_engine
    return get_engine(name)


def make_cell_class(eng):
    class PackageCell(eng.Cell):
        pass
    return PackageCell


def load(pkg, verify=True, check_active_engine=True, evaluator=None, rules=None):
    """load a cellgeom_v1 package -> evaluator Cell with manifest/vocab/summary. verify: True (sha256), 'quick' (sizes)
    or False. evaluator: None = the active one (cellgeom_eval_v1, engine.py); it reads packages of the layouts it lists
    (LAYOUTS = manifest engine.name) with a rule table: `rules` (path), else config/geometry/rules/<layout>.yaml
    (spec §4.3). evaluator = the manifest engine name (e.g. sdf_sections_v1) loads with that legacy engine (no tables;
    reproduction tests only, needs check_active_engine=False)."""
    pkg = os.path.abspath(pkg)
    mpath = os.path.join(pkg, "manifest.json")
    man = json.load(open(mpath))
    if man.get("spec_version") != SPEC_VERSION:
        raise ValueError(f"{pkg}: spec_version {man.get('spec_version')} != {SPEC_VERSION}")
    if verify:
        verify_checksums(pkg, man, quick=(verify == "quick"))
    ename = man["engine"]["name"]
    from . import engine as active
    evaluator = evaluator or active.ACTIVE_ENGINE
    if check_active_engine and evaluator != active.ACTIVE_ENGINE:
        raise RuntimeError(f"evaluator {evaluator} != active evaluator {active.ACTIVE_ENGINE}")
    eng = _engine_module(evaluator)
    tables = rule_table = None
    if evaluator != ename:
        if ename not in getattr(eng, "LAYOUTS", ()):
            raise RuntimeError(f"evaluator {evaluator} cannot read packages of layout {ename}")
        from . import rules as _rules
        rule_table = _rules.load_rules(rules or _rules.default_rules_path(ename))
        if rule_table.get("layout") != ename:
            raise ValueError(f"rule table {rule_table['_path']} is for layout {rule_table.get('layout')}, package is {ename}")
        cf = (rule_table.get("constants") or {}).get("file")
        if cf and cf not in man["files"]:
            raise ValueError(f"rule table constants file {cf} is not a checksummed package file")
        tables = _rules.compile_for_manifest(rule_table, man, pkg)
    nslot = 1 + max(f["slot"] for f in man["fields"].values())
    fields = [None] * nslot
    for name, f in man["fields"].items():
        fields[f["slot"]] = eng.empty_field() if f.get("empty") else read_field_npz(os.path.join(pkg, f["file"]))
    for k in range(nslot):
        if fields[k] is None:
            fields[k] = eng.empty_field()
    kw = {}
    for name, s in man["sections"].items():
        key = eng.SECTION_TYPES.get(s["type"])
        if key is None:
            raise ValueError(f"engine {ename} cannot evaluate section type {s['type']} ({name})")
        z = np.load(os.path.join(pkg, s["file"]))
        if key == "golgi_cis":
            kw[key] = (float(z["W"]), float(z["S"]), (z["header"][:10].astype(np.float64),
                                                      z["index"].astype(np.float64), z["data"].astype(np.float64)))
        else:
            kw[key] = z["section"].astype(np.float64)
    dom = man["domain"]
    face = None
    for cp in dom.get("clip_planes", []):
        n = np.asarray(cp["normal"], float)
        if not np.allclose(n, [0, 0, 1]) or cp.get("kept_side") != "+normal":
            raise ValueError(f"engine {ename} supports only z clip planes with kept side +normal")
        face = float(cp["point"][2])
    if tables is not None:
        kw["tables"] = tables
    G = eng.build_G(dom["enclosing_radius_nm"], fields, man["engine"]["params"], face=face, **kw)
    pores = np.load(os.path.join(pkg, man["pores"]["file"]))["pores"] if man.get("pores") else None
    cell = make_cell_class(eng)(G, pores)
    if man.get("prims") and hasattr(cell, "set_prims"):
        z = np.load(os.path.join(pkg, man["prims"]["file"]))
        tmap = np.full(1 + max(man["prims"]["types"].values()), -1, np.int64)
        for tname, code in man["prims"]["types"].items():
            if tname not in eng.PRIM_TYPES:
                raise ValueError(f"evaluator {evaluator} has no prim type {tname}")
            tmap[code] = eng.PRIM_TYPES[tname]
        cell.set_prims(tmap[z["type"].astype(np.int64)], z["params"], z["comp"].astype(np.int64))
        if "bsphere" in z.files:                 # bounding spheres (centre, radius): tile-local prim selection
            cell.prim_bsphere = np.ascontiguousarray(z["bsphere"], np.float64).reshape(-1, 4)
    cell.rule_params = _rules.rule_params(rule_table, man, pkg) if rule_table else dict(man["engine"]["params"])
    cell.manifest = man
    cell.manifest_sha256 = sha256_file(mpath)
    cell.package_dir = pkg
    cell.package_id = man["package_id"]
    cell.engine = eng
    cell.evaluator = eng.ENGINE_NAME
    cell.rules = rule_table
    cell.rules_path = rule_table["_path"] if rule_table else None
    cell.rules_sha256 = rule_table["_sha256"] if rule_table else None
    cell.vocab = {L["id"]: L for L in man["labels"]}
    cell.code_to_name = {L["engine_code"]: L["name"] for L in man["labels"] if L.get("engine_code") is not None}
    cell.name_to_code = {v: k for k, v in cell.code_to_name.items()}
    cell.code_to_id = {L["engine_code"]: L["id"] for L in man["labels"] if L.get("engine_code") is not None}
    mcodes = {m["name"]: m["engine_code"] for m in man["membranes"]}
    cell.mem_code = mcodes.__getitem__
    cell.summary = man.get("summary", {})
    cell.layer_path = lambda name, i=0: os.path.join(pkg, man["layers"][name]["files"][i])
    return cell


def load_cfg(cfg, pkg=None, **kw):
    """load the package of a cell config: geometry.package (or pkg) with geometry.evaluator (default and only core
    value: cellgeom_eval_v1) and geometry.rules (rule table path; default config/geometry/rules/<layout>.yaml)"""
    g = cfg.get("geometry") or {}
    from . import engine as active
    ev = g.get("evaluator", active.ACTIVE_ENGINE)
    if ev != active.ACTIVE_ENGINE:
        raise RuntimeError(f"config geometry.evaluator {ev}: the core evaluator is {active.ACTIVE_ENGINE} (legacy "
                           f"evaluators are for the reproduction tests only)")
    return load(pkg or g["package"], evaluator=ev, rules=g.get("rules"), **kw)


# ----------------------------------------------------------------- declarative label rule
def _val(v, params):
    if isinstance(v, (int, float)):
        return float(v)
    return float(eval(v, {"__builtins__": {}}, {k: x for k, x in params.items() if not isinstance(x, list)}))


def eval_label_rule(cell, P, rule=None):
    """evaluate manifest['label_rule'] at points P -> label names array (object) ('' = no rule matched)"""
    man = cell.manifest; eng = cell.engine; P = np.asarray(P, np.float64)
    rule = rule or man["label_rule"]
    params = man["engine"]["params"]
    cache = {}

    def mems():
        if "m" not in cache:
            cache["m"] = cell.mems64(P)
        return cache["m"]

    def field(name):
        k = "f:" + name
        if k not in cache:
            f = man["fields"].get(name)
            cache[k] = (np.full(len(P), np.inf) if f is None or f.get("empty")
                        else eng.field_points(P, cell.G, f["slot"]))
        return cache[k]

    def cond(c):
        if "clip" in c:
            cp = next(x for x in man["domain"]["clip_planes"] if x["name"] == c["clip"])
            n = np.asarray(cp["normal"], float); q = np.asarray(cp["point"], float)
            s = (P - q) @ n if not np.allclose(n, [0, 0, 1]) else P[:, 2] - q[2]
            return s < 0
        if "mem" in c:
            x = mems()[:, cell.mem_code(c["mem"])]
        elif "field" in c:
            x = field(c["field"])
        elif c.get("var") == "matrix":
            x = -mems()[:, cell.mem_code(c.get("membrane", "MITO_IM"))]
        elif "in_pore" in c:
            if cell.pores.shape[0] == 0:
                return np.zeros(len(P), bool)
            return eng.pore_points(P, cell.pores, _val(c["in_pore"], params), _val(c.get("zwin", "PORE_ZWIN"), params)) >= 0
        elif "in_section" in c:
            if c["in_section"] not in man["sections"]:
                return np.zeros(len(P), bool)
            return eng.speckle_points(P, cell.G) >= 0
        elif "has_section" in c:
            return np.full(len(P), c["has_section"] in man["sections"])
        else:
            raise ValueError(f"unknown label-rule test {c}")
        if "lt" in c:
            return x < _val(c["lt"], params)
        return x > _val(c["gt"], params)

    out = np.full(len(P), "", dtype=object)

    def apply(rules, mask):
        rem = mask.copy()
        for r in rules:
            m = rem.copy()
            for c in r.get("if", []):
                if not m.any():
                    break
                m &= cond(c)
            if "label" in r:
                out[m] = r["label"]
            else:
                apply(r["then"], m)
            rem &= ~m
    apply(rule, np.ones(len(P), bool))
    return out


# ----------------------------------------------------------------- validator
def _box(cell):
    """box of the 'cell' field grid (falls back to the enclosing sphere's cube)"""
    f = cell.manifest["fields"].get("cell")
    if f and not f.get("empty"):
        lo = np.asarray(f["origin_nm"]); hi = lo + f["h_nm"] * f["brick"] * np.asarray(f["bricks"])
        return lo, hi
    R = cell.manifest["domain"]["enclosing_radius_nm"]
    return -np.full(3, R), np.full(3, R)


def mc_volumes(cell, n=2_000_000, seed=0, stream="geometry.mc_volumes"):
    """region volumes (um3) by Monte Carlo over the cell-field box, keyed by vocabulary name; + 95% CI half-width"""
    from wholecell import seeds
    rng = seeds.rng({"seeds": {"master": int(seed)}}, stream)
    lo, hi = _box(cell)
    P = lo + (hi - lo) * rng.random((n, 3))
    lab = cell.labels(P)
    V = float(np.prod(hi - lo)) / 1e9
    out, ci = {}, {}
    for code, name in sorted(cell.code_to_name.items()):
        p = float((lab == code).mean())
        out[name] = round(p * V, 4); ci[name] = round(1.96 * np.sqrt(p * (1 - p) / n) * V, 4)
    return dict(n=n, seed=int(seed), stream=stream, box_lo=lo.tolist(), box_hi=hi.tolist(), box_um3=V,
                volume_um3=out, ci95_um3=ci)


def validate(cell, n_rule=400_000, n_mc=2_000_000, n_area=4_000_000, seed=0):
    """reference validator (spec section 'Validation'). Returns a dict; 'pass' = all hard checks pass."""
    from wholecell import seeds
    if isinstance(cell, str):
        cell = load(cell)
    man = cell.manifest; eng = cell.engine
    rng = seeds.rng({"seeds": {"master": int(seed)}}, "geometry.validate")
    res = {"spec_version": SPEC_VERSION, "package_id": man["package_id"], "engine": man["engine"]["name"],
           "seed": int(seed), "created": _dt.datetime.now().isoformat(timespec="seconds"), "checks": {}}
    C = res["checks"]
    C["checksums"] = {"files": len(man["files"]), "pass": True}          # load() verified them
    # schema / vocabulary
    ids = [L["id"] for L in man["labels"]]; codes = [L["engine_code"] for L in man["labels"] if L.get("engine_code") is not None]
    errs = []
    if len(set(ids)) != len(ids): errs.append("duplicate label ids")
    if len(set(codes)) != len(codes): errs.append("duplicate engine codes")
    if any(c < 0 or c > MAX_ENGINE_CODE for c in codes): errs.append("engine code outside 0..63")
    for L in man["labels"]:
        if L["id"] < EXT_ID_MIN and CORE_BY_NAME.get(L["name"]) != L["id"]:
            errs.append(f"core label {L['name']} id {L['id']} not in the core table")
        if L["id"] >= EXT_ID_MIN and CORE_BY_NAME.get(L.get("parent")) is None:
            errs.append(f"extension label {L['name']} needs a core parent")
    mh = man["engine"]["params"].get("MEM_HALF")
    for m in man["membranes"]:
        if mh is not None and abs(float(m["half_thickness_nm"]) - float(mh[m["engine_code"]])) > 1e-12:
            errs.append(f"membrane {m['name']} half_thickness_nm != engine MEM_HALF")
    for r in _rule_labels(man["label_rule"]):
        if r not in {L["name"] for L in man["labels"]}: errs.append(f"label rule names unknown label {r}")
    C["schema"] = {"errors": errs, "pass": not errs}
    lo, hi = _box(cell)
    # label rule vs engine
    nu = n_rule // 2
    U = lo + (hi - lo) * rng.random((nu, 3))
    cand = lo + (hi - lo) * rng.random((n_rule * 4, 3))
    M = cell.mems64(cand)
    A = np.abs(M[:, :eng.NMEM - 1])
    near = cand[(A < 50.0).any(1)][:n_rule - nu]
    P = np.vstack([U, near])
    lab = cell.labels(P)
    names = np.array([cell.code_to_name.get(int(c), f"?{c}") for c in lab], dtype=object)
    rl = eval_label_rule(cell, P)
    agree = float((names == rl).mean())
    undefined = int(sum(1 for c in np.unique(lab) if int(c) not in cell.code_to_name))
    C["label_rule"] = {"n": int(len(P)), "agreement": agree, "mismatch": int((names != rl).sum()),
                       "undefined_engine_codes": undefined, "unmatched_rule": int((rl == "").sum()),
                       "pass": agree == 1.0 and undefined == 0}
    # eikonal sanity per stored field
    eik = {}
    for name, f in man["fields"].items():
        if f.get("empty"):
            continue
        flo = np.asarray(f["origin_nm"]); fhi = flo + f["h_nm"] * f["brick"] * np.asarray(f["bricks"])
        Q = flo + (fhi - flo) * rng.random((400_000, 3))
        d = eng.field_points(Q, cell.G, f["slot"])
        band = min(f["clamp_in_nm"], f["clamp_out_nm"]) - 2 * f["h_nm"]
        Q = Q[np.abs(d) < band][:50_000]
        if len(Q) < 100:
            eik[name] = {"n": int(len(Q)), "pass": True, "note": "too few band points"}; continue
        s = f["h_nm"] / 2
        g = np.zeros(len(Q))
        for ax in range(3):
            e = np.zeros(3); e[ax] = s
            g += ((eng.field_points(Q + e, cell.G, f["slot"]) - eng.field_points(Q - e, cell.G, f["slot"])) / (2 * s)) ** 2
        g = np.sqrt(g)
        med = float(np.median(g))
        eik[name] = {"n": int(len(Q)), "median_grad": round(med, 4), "p05": round(float(np.percentile(g, 5)), 4),
                     "p95": round(float(np.percentile(g, 95)), 4), "pass": 0.95 <= med <= 1.05}
    C["eikonal"] = {"fields": eik, "pass": all(v["pass"] for v in eik.values())}
    # nesting: organelle interiors inside the cell (and inside the kept domain)
    Q = lo + (hi - lo) * rng.random((n_mc, 3))
    cellf = eng.field_points(Q, cell.G, man["fields"]["cell"]["slot"]) if "cell" in man["fields"] else None
    nest = {}
    for name, f in man["fields"].items():
        if f.get("empty") or name == "cell" or cellf is None:
            continue
        d = eng.field_points(Q, cell.G, f["slot"])
        inside = d < 0
        frac = float((cellf[inside] < 0).mean()) if inside.any() else 1.0
        nest[name] = {"n_inside": int(inside.sum()), "frac_inside_cell": round(frac, 5), "pass": frac >= 0.99}
    C["nesting"] = {"fields": nest, "pass": all(v["pass"] for v in nest.values())}
    # volumes (MC labels) and field volumes / areas vs summary
    vol = mc_volumes(cell, n=n_mc, seed=seed)
    C["volumes_mc"] = vol
    zf = None
    for cp in man["domain"].get("clip_planes", []):
        zf = float(cp["point"][2])
    summ = man.get("summary", {})
    fa = {}
    for name, f in man["fields"].items():
        if f.get("empty"):
            continue
        flo = np.asarray(f["origin_nm"]); fhi = flo + f["h_nm"] * f["brick"] * np.asarray(f["bricks"])
        if zf is not None:
            flo = flo.copy(); flo[2] = max(flo[2], zf)
        Vb = float(np.prod(fhi - flo)) / 1e9
        Q = flo + (fhi - flo) * rng.random((n_area, 3))
        d = eng.field_points(Q, cell.G, f["slot"])
        p = float((d < 0).mean()); eps = 4.0
        nsh = int((np.abs(d) < eps).sum())
        a = nsh / n_area * Vb * 1e9 / (2 * eps) / 1e6
        ent = {"vol_um3": round(p * Vb, 4), "vol_ci95": round(1.96 * np.sqrt(p * (1 - p) / n_area) * Vb, 4),
               "area_um2_shell": round(a, 3), "area_rel_sd": round(1 / np.sqrt(max(nsh, 1)), 4)}
        s = summ.get(name, {}) if isinstance(summ.get(name), dict) else {}
        sv = s.get("vol_sdf_um3", s.get("vol_inside_um3")); sa = s.get("area_um2")
        ok = True
        if sv is not None:
            ent["summary_vol_um3"] = sv; tol = 4 * ent["vol_ci95"] / 1.96 + 0.01 * sv + 0.02
            ent["vol_ok"] = bool(abs(ent["vol_um3"] - sv) <= tol); ok &= ent["vol_ok"]
        if sa is not None:
            ent["summary_area_um2"] = sa; ent["area_rel_diff"] = round(a / sa - 1, 4) if sa else None
            ent["area_ok"] = bool(sa == 0 or abs(a / sa - 1) <= 0.05 + 4 * ent["area_rel_sd"]); ok &= ent["area_ok"]
        ent["pass"] = bool(ok)
        fa[name] = ent
    C["field_volumes_areas"] = {"fields": fa, "pass": all(v["pass"] for v in fa.values())}
    # clip planes: everything beyond is outside_imaged_volume
    cps = []
    for cp in man["domain"].get("clip_planes", []):
        Q = lo + (hi - lo) * rng.random((200_000, 3))
        below = Q[:, 2] < float(cp["point"][2])
        lb = cell.labels(Q[below])
        code = cell.name_to_code.get(cp.get("label", "outside_imaged_volume"))
        cps.append({"name": cp["name"], "n_beyond": int(below.sum()), "frac_labelled": float((lb == code).mean()) if len(lb) else 1.0})
    C["clip_planes"] = {"planes": cps, "pass": all(c["frac_labelled"] == 1.0 for c in cps)}
    # pores on the NE midsurface
    if cell.pores.shape[0] and "ne_mid" in man["fields"] and not man["fields"]["ne_mid"].get("empty"):
        from scipy.spatial import cKDTree
        k = man["fields"]["ne_mid"]["slot"]; Pc = cell.pores[:, :3]; Nn = cell.pores[:, 3:]
        d = eng.field_points(Pc, cell.G, k)
        s = 2.0; g = np.zeros_like(Pc)
        for ax in range(3):
            e = np.zeros(3); e[ax] = s
            g[:, ax] = (eng.field_points(Pc + e, cell.G, k) - eng.field_points(Pc - e, cell.G, k)) / (2 * s)
        g /= np.linalg.norm(g, axis=1, keepdims=True)
        ang = np.degrees(np.arccos(np.clip((g * Nn).sum(1), -1, 1)))
        dd, _ = cKDTree(Pc).query(Pc, k=2)
        minsep = man["pores"].get("min_sep_nm", 0.0)
        tol_d = man["fields"]["ne_mid"]["h_nm"] / 4
        C["pores"] = {"n": int(len(Pc)), "counts": man["pores"].get("counts"),
                      "max_abs_ne_mid_nm": round(float(np.abs(d).max()), 3), "median_abs_ne_mid_nm": round(float(np.median(np.abs(d))), 3),
                      "axis_vs_grad_deg_median": round(float(np.median(ang)), 2), "axis_vs_grad_deg_p95": round(float(np.percentile(ang, 95)), 2),
                      "min_spacing_nm": round(float(dd[:, 1].min()), 1), "declared_min_sep_nm": minsep,
                      "pass": bool(np.abs(d).max() < tol_d and np.median(ang) < 10 and dd[:, 1].min() >= minsep - 1e-6)}
    # layers
    lay = {}
    for name, L in man.get("layers", {}).items():
        ok = L.get("depends_on") == man["package_id"]
        missing = [f for f in L.get("files", []) if f not in man["files"]]
        refs = []
        for r in L.get("refs", []):
            ex = os.path.exists(r["path"])
            refs.append({"path": r["path"], "exists": ex,
                         "md5_ok": (md5_file(r["path"]) == r["md5"]) if ex and r.get("md5") else None})
        lay[name] = {"depends_on_ok": ok, "missing_files": missing, "refs": refs,
                     "pass": ok and not missing and all(x["exists"] and x["md5_ok"] is not False for x in refs)}
    C["layers"] = {"layers": lay, "pass": all(v["pass"] for v in lay.values())}
    finish(res, man)
    return res


def finish(res, man):
    """apply manifest validation_waivers ({check: reason}; documented legacy defects) and set res['pass']"""
    waiv = man.get("validation_waivers", {}) or {}
    C = res["checks"]
    for k, why in waiv.items():
        if k in C:
            C[k]["waived"] = why
    failed = [k for k, c in C.items() if isinstance(c, dict) and not c.get("pass", True)]
    res["failed_checks"] = failed
    res["waived_failures"] = [k for k in failed if k in waiv]
    res["pass"] = not [k for k in failed if k not in waiv]
    return res


def _rule_labels(rules):
    for r in rules:
        if "label" in r:
            yield r["label"]
        else:
            yield from _rule_labels(r["then"])
