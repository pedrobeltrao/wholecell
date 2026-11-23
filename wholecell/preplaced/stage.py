"""Pre-placed structures stage (`preplaced`; called "prepass" before A0c): fixed obstacles placed before tile packing (chromatin, polysomes, cytoskeleton, NPC, lamina, RNA beads).

    python -m wholecell preplaced --cell C --build B --step module:<name>   (build / compare one module with a builder)
    python -m wholecell preplaced --cell C --build B --step manifest        (write preplaced_manifest.json for the pack)
    python -m wholecell preplaced --cell C --build B --step serial          (ordered serial builders on one shared context:
                                                                           preplaced.serial; see serial.py; yeast, A6)

Config:
    preplaced:
      modules:                      # order = pre-pass library order in the pack index (legacy jwc config order)
        - {name: chromatin, kind: instances, mode: frozen, npz: <path>, json: <path>}
        - {name: lamina, kind: instances, mode: build|frozen, builder: lamina_lattice, seed_stream: prepass.lamina,
           params: {...}, inputs: {chromatin: "module:chromatin", table: "build:inventory/rules/v6/ingredients.tsv",
                                   ne_mesh: <path>},
           frozen: {npz: <path>, json: <path>, drop_json: <path>}}       # reference copy (B5 comparison / fallback)
        - {name: mrna, kind: beads, mode: frozen, npz: <path>, drop_npy: <path>}
        - {name: polysomes, kind: instances, mode: serial, npz: "build:preplaced/serial/modules/polysomes.npz", json: ...}
                                    # a module file written by the serial pass (preplaced.serial.exports; B4)
      drops:                        # frozen drop lists (instances of one module displaced by a later one)
        - {target: chromatin, json: <legacy drop json: n, drop, why>}
Paths are config strings (placeholders resolved; frozen files are declared with checksums in frozen_inputs and
verified by the `frozen` stage). `module:<name>` = the npz selected for that module; `build:<rel>` = a file of this
build. Builders are registered `preplaced` components: build(geom, LIB, chrom, ing, ne_mesh, out_dir, seed, params).

Module with a builder: `module:<name>` builds into <bdir>/preplaced/<name>/build/, compares with the frozen reference
(`compare_frozen.json`: identical arrays or not; criterion B5) and writes <bdir>/preplaced/<name>/selected.json plus
<bdir>/preplaced/<name>/consumption.json (byte copy of the selected module json: the explicit input of the inventory rule
that consumes the pre-pass copies, e.g. Jurkat v7 rule L3). mode: build -> the rebuilt files are used; mode: frozen ->
the frozen reference is used (the rebuild is still run and compared when `rebuild_check: true`).
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil


def _md5(p, chunk=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def _mods(cfg):
    return (cfg.get("preplaced") or {}).get("modules") or []


def _module(cfg, name):
    for m in _mods(cfg):
        if m["name"] == name:
            return m
    raise KeyError(f"preplaced module '{name}' not in config")


def _sel_path(bdir, name):
    return os.path.join(bdir, "preplaced", name, "selected.json")


def _resolve(cfg, bdir, v):
    if v is None:
        return None
    if v.startswith("module:"):
        return selected(cfg, bdir, v.split(":", 1)[1])["npz"]
    if v.startswith("build:"):
        return os.path.join(bdir, v.split(":", 1)[1])
    if v.startswith("frozen:"):
        return cfg.frozen(v.split(":", 1)[1])
    if v in (cfg.get("frozen_inputs") or {}):          # a bare frozen-input name
        return cfg.frozen(v)
    return v


def selected(cfg, bdir, name):
    """files used for a module: {npz, json, drops: [...]}"""
    m = _module(cfg, name)
    if m.get("builder"):
        p = _sel_path(bdir, name)
        if not os.path.exists(p):
            raise FileNotFoundError(f"{p}: run preplaced --step module:{name} first")
        return json.load(open(p))
    out = {"npz": _resolve(cfg, bdir, m["npz"]), "source": "serial" if m.get("mode") == "serial" else "frozen"}
    if m.get("json"):
        out["json"] = _resolve(cfg, bdir, m["json"])
    if m.get("drop_npy"):
        out["drop_npy"] = _resolve(cfg, bdir, m["drop_npy"])
    return out


def _compare_npz(a, b):
    import numpy as np
    za, zb = np.load(a), np.load(b)
    res = {"files_a": sorted(za.files), "files_b": sorted(zb.files), "arrays": {}}
    same = sorted(za.files) == sorted(zb.files)
    for k in sorted(set(za.files) & set(zb.files)):
        x, y = za[k], zb[k]
        eq = bool(x.shape == y.shape and x.dtype == y.dtype and np.array_equal(x, y))
        r = {"identical": eq, "shape_a": list(x.shape), "shape_b": list(y.shape)}
        if not eq and x.shape == y.shape and x.dtype.kind == "f":
            r["max_abs_diff"] = float(np.abs(x.astype(np.float64) - y.astype(np.float64)).max())
        res["arrays"][k] = r
        same &= eq
    res["identical"] = bool(same)
    return res


def _build_module(cfg, bdir, name):
    from wholecell import registry, seeds
    from wholecell.pack import geomload
    m = _module(cfg, name)
    mdir = os.path.join(bdir, "preplaced", name)
    sp = _sel_path(bdir, name)
    if os.path.exists(sp):
        print(f"preplaced {name}: selected.json exists, step skipped", flush=True)
        return
    bld = os.path.join(mdir, "build")
    fz = m.get("frozen") or {}
    do_build = m.get("mode", "build") == "build" or m.get("rebuild_check", False)
    J = None
    if do_build:
        pre = (m.get("params") or {}).get("out_prefix", "lamina_cell1")
        if not os.path.exists(os.path.join(bld, f"{pre}.json")):
            fn = registry.get("preplaced", m["builder"])
            geom = geomload.load(cfg, bdir)
            lib = _library(bdir)
            inp = {k: geom.resolve_input(_resolve(cfg, bdir, v)) for k, v in (m.get("inputs") or {}).items()}
            seed = seeds.seed_int(cfg, m.get("seed_stream", f"prepass.{name}"))   # stream name kept at A0c (= the seed)
            chrom_name = (m.get("inputs") or {}).get("chromatin", "module:chromatin").split(":", 1)[-1]
            tmpd = bld + f".tmp{os.getpid()}"
            if os.path.exists(tmpd):
                shutil.rmtree(tmpd)
            fn(geom, lib, inp["chromatin"], inp["table"], inp["ne_mesh"], tmpd, seed, m.get("params") or {}, chrom_name=chrom_name)
            if os.path.exists(bld):
                shutil.rmtree(bld)
            os.replace(tmpd, bld)
        J = json.load(open(os.path.join(bld, f"{pre}.json")))
        built = dict(npz=os.path.join(bld, f"{pre}.npz"), json=os.path.join(bld, f"{pre}.json"))
        dn = (m.get("params") or {}).get("drop_name", "chromatin_drop_lamina")
        dj = os.path.join(bld, f"{dn}.json")
        built["drops"] = [_drop_entry(json.load(open(dj)), dj)] if os.path.exists(dj) else []
        if fz.get("npz"):
            cmp = {"npz": _compare_npz(built["npz"], _resolve(cfg, bdir, fz["npz"]))}
            fj = json.load(open(_resolve(cfg, bdir, fz["json"]))) if fz.get("json") else {}
            cmp["consumed_identical"] = J.get("consumed") == fj.get("consumed")
            cmp["consumed_build"], cmp["consumed_frozen"] = J.get("consumed"), fj.get("consumed")
            cmp["dropped_build"], cmp["dropped_frozen"] = J.get("dropped"), fj.get("dropped")
            if fz.get("drop_json"):
                import numpy as np
                fdj = json.load(open(_resolve(cfg, bdir, fz["drop_json"])))
                a = np.load(built["drops"][0]["npy"]) if built["drops"] else np.zeros(0, np.int64)
                b = np.load(_drop_npy_path(cfg, fdj, _resolve(cfg, bdir, fz["drop_json"])))
                cmp["drop_list_identical"] = bool(a.shape == b.shape and np.array_equal(a, b))
            cmp["identical"] = bool(cmp["npz"]["identical"] and cmp["consumed_identical"] and cmp.get("drop_list_identical", True))
            json.dump(cmp, open(os.path.join(mdir, "compare_frozen.json"), "w"), indent=1)
            print(f"preplaced {name}: rebuilt vs frozen identical = {cmp['identical']}", flush=True)
    if m.get("mode", "build") == "build":
        sel = dict(built, source="build")
    else:
        sel = dict(npz=_resolve(cfg, bdir, fz["npz"]), json=_resolve(cfg, bdir, fz["json"]), source="frozen")
        if fz.get("drop_json"):
            p = _resolve(cfg, bdir, fz["drop_json"])
            sel["drops"] = [_drop_entry(json.load(open(p)), p, cfg)]
    os.makedirs(mdir, exist_ok=True)
    shutil.copyfile(sel["json"], os.path.join(mdir, "consumption.json"))
    tmp = sp + f".tmp{os.getpid()}"
    json.dump(sel, open(tmp, "w"), indent=1)
    os.replace(tmp, sp)
    print(f"preplaced {name}: using {sel['source']} files ({sel['npz']})", flush=True)


def _drop_npy_path(cfg, dj, jpath):
    """drop npy of a drop json: absolute, relative to the json folder, or (legacy) relative to data_root"""
    d = dj["drop"]
    if os.path.isabs(d):
        return d
    local = os.path.join(os.path.dirname(jpath), d)      # wholecell builders: relative to the json folder
    if os.path.exists(local):
        return local
    root = (cfg or {}).get("_site", {}).get("data_root", "") if cfg else ""
    cand = os.path.join(root, d)                         # legacy jsons: relative to the legacy workspace root
    return cand if os.path.exists(cand) else os.path.join(os.path.dirname(jpath), os.path.basename(d))


def _drop_entry(dj, jpath, cfg=None):
    tgt = dj["target"]
    if "/" in tgt or tgt.endswith(".npz"):           # legacy json: target is a path; module named by the drop config
        tgt = None
    return dict(target=tgt, n=int(dj["n"]), npy=_drop_npy_path(cfg, dj, jpath), why=dj.get("why", ""), json=jpath)


def _library(bdir):
    import pandas as pd
    L = pd.read_csv(os.path.join(bdir, "library", "library_manifest.tsv"), sep="\t")
    L = L[L.path.str.endswith(".npz")]
    return dict(zip(L.species, L.path))


def package_id(cfg, bdir):
    """geometry package id of this build (geometry/package.json, else the manifest of config geometry.package)"""
    pj = os.path.join(bdir, "geometry", "package.json")
    if os.path.exists(pj) and json.load(open(pj)).get("package_id"):
        return json.load(open(pj))["package_id"]
    pkg = (cfg.get("geometry") or {}).get("package")
    return json.load(open(os.path.join(pkg, "manifest.json")))["package_id"] if pkg else None


def check_depends_on_geometry(cfg, bdir):
    """`depends_on_geometry` of the frozen pre-placed modules (docs/reference/cell_config.md 3.4): a module whose files
    are used from a frozen input must declare the geometry package it was built on / checked against, and it must be
    this build's package. Modules built in this build depend on its package by construction. -> (rows, mismatches)"""
    pid = package_id(cfg, bdir)
    rows, bad = [], []
    for m in _mods(cfg):
        built = (bool(m.get("builder")) and m.get("mode", "frozen") == "build") or m.get("mode") == "serial"   # serial: B4
        decl = m.get("depends_on_geometry")
        dep = pid if built else decl
        ok = dep == pid
        rows.append(dict(name=m["name"], source="build" if built else "frozen", depends_on_geometry=dep, package_id=pid, ok=ok))
        if not ok:
            bad.append(f"{m['name']}: depends_on_geometry {decl!r} != package {pid!r}")
    return rows, bad


def _manifest(cfg, bdir):
    out = os.path.join(bdir, "preplaced", "preplaced_manifest.json")
    if os.path.exists(out):
        print("preplaced manifest exists, step skipped", flush=True)
        return
    deps, bad = check_depends_on_geometry(cfg, bdir)
    if bad:
        raise RuntimeError("preplaced: frozen modules not declared for this geometry (B6): " + "; ".join(bad))
    deps = {d["name"]: d["depends_on_geometry"] for d in deps}
    fi = cfg.get("frozen_inputs") or {}
    by_path = {v["path"]: (k, v.get("md5") or v.get("sha256")) for k, v in fi.items()}
    inst, beads, drops = [], [], []
    for m in _mods(cfg):
        s = selected(cfg, bdir, m["name"])
        e = dict(name=m["name"], npz=s["npz"], source=s.get("source", "frozen"), depends_on_geometry=deps[m["name"]])
        if s.get("json"):
            e["json"] = s["json"]
        if s.get("drop_npy"):
            e["drop_npy"] = s["drop_npy"]
        e["frozen_input"], e["checksum"] = by_path.get(s["npz"], (None, None))
        if e["checksum"] is None and s.get("source") in ("build", "serial"):
            e["checksum"] = _md5(s["npz"])
        (inst if m.get("kind", "instances") == "instances" else beads).append(e)
        for d in s.get("drops", []) or []:
            d = dict(d)
            d["target"] = d.get("target") or m.get("drops_target") or "chromatin"
            drops.append(d)
    for d in (cfg.get("preplaced") or {}).get("drops") or []:
        p = _resolve(cfg, bdir, d["json"])
        e = _drop_entry(json.load(open(p)), p, cfg); e["target"] = d["target"]
        drops.append(e)
    names = [m["name"] for m in inst]
    for d in drops:
        if d["target"] not in names:
            raise ValueError(f"preplaced drop list targets unknown module {d['target']}")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    tmp = out + f".tmp{os.getpid()}"
    json.dump(dict(instances=inst, beads=beads, drops=drops), open(tmp, "w"), indent=1)
    os.replace(tmp, out)
    print(f"preplaced manifest: {len(inst)} instance modules, {len(beads)} bead sets, {len(drops)} drop lists", flush=True)


def run(cfg, bdir, step=None, **kw):
    import os as _os
    _os.environ.setdefault("NUMBA_NUM_THREADS", str(int((cfg.get("preplaced") or {}).get("numba_threads", 8))))
    from wholecell import registry
    from . import builders  # noqa: F401  (registers core builders)
    registry.load_plugins(cfg)
    if step and step.startswith("module:"):
        _build_module(cfg, bdir, step.split(":", 1)[1])
    elif step == "manifest":
        _manifest(cfg, bdir)
    elif step == "serial":
        from . import serial
        serial.run(cfg, bdir)
    else:
        raise SystemExit("preplaced: --step module:<name> | serial | manifest")
