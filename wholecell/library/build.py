"""Library mode `build`: the yeast library policy in the core (unification plan §1, A3). See stage.py for the config.

One plain function per step (STEPS); each takes the context `X` (paths, config, species table, log) and writes into
the build's in-progress library folder <bdir>/library/lib.partial/. `seal` renames it to lib/ (made read-only) and
writes the manifest; nothing edits a sealed library. Port of yeast scripts/s0_library.py (generic part) @ e3b7448;
organism data (curated complexes and their orientation, curated PDB oligomers, reference RNA) are config / TSV.
Jurkat builder steps (B3: jurkat_assemblies, jurkat_monomers, reorient_tm, tm_route, adaptive_trees) live in jurkat.py.
"""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import stat
import time


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def read_any_manifest(path):
    """[(file, md5)] of a manifest: core TSV (header file, md5, bytes) or the yeast make_manifests.sh format
    (`md5  size  mtime  ./relpath`)"""
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if lines and lines[0].split("\t")[:2] == ["file", "md5"]:
        return [tuple(ln.split("\t")[:2]) for ln in lines[1:]]
    out = []
    for ln in lines:
        md5, _size, _mtime, rel = ln.split("  ", 3)
        out.append((rel[2:] if rel.startswith("./") else rel, md5))
    return out


def _path(cfg, bdir, v):
    """a path value: frozen input name -> its path; `{build}` -> the build directory; relative -> repo root"""
    if v is None:
        return None
    fi = cfg.get("frozen_inputs") or {}
    if isinstance(v, str) and v in fi:
        return fi[v]["path"]
    v = str(v).replace("{build}", bdir)
    if not os.path.isabs(v):
        v = os.path.join(cfg["_site"]["repo"] if cfg.get("_site", {}).get("repo") else os.getcwd(), v)
    return v


def _tsv(path):
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment=None)


class Ctx(dict):
    __getattr__ = dict.__getitem__


def context(cfg, bdir, nproc=None):
    lc = cfg["library"]
    out = os.path.join(bdir, "library")
    X = Ctx(cfg=cfg, lc=lc, bdir=bdir, out=out, lib=os.path.join(out, "lib"), partial=os.path.join(out, "lib.partial"),
            nproc=int(nproc or lc.get("nproc", 8)), origin={}, results={})
    X["pdb_dir"] = _path(cfg, bdir, lc.get("pdb_dir"))
    X["species_path"] = _path(cfg, bdir, lc.get("species"))
    return X


def species(X):
    if "species" not in X:
        if not X.species_path or not os.path.exists(X.species_path):
            raise FileNotFoundError(f"library.species: species table {X.species_path} missing (TSV ing_id, kind, uniprot, pdb; "
                                    f"written by the assembly stage, for yeast by plugins/yeast/library_species.py)")
        X["species"] = _tsv(X.species_path)
    return X["species"]


def topology(X):
    if "topology" not in X:
        import pandas as pd
        t = X.lc["topology"]
        topo = pd.read_csv(_path(X.cfg, X.bdir, t["path"]), sep="\t")
        topo.columns = list(t.get("columns", ["acc", "topo", "tm", "intramem"]))
        X["topology"] = topo.set_index("acc")
    return X["topology"]


def _mark(X, ids, step):
    for i in ids:
        X.origin[i] = step


# ------------------------------------------------------------------ steps
def step_start(X, p):
    """copy the declared frozen start entries (every file md5-checked against the manifest; copy2 keeps mtimes).
    library.start: <frozen input name> (the whole folder = its manifest) or {frozen: <name>, files: <TSV file, md5,
    reason>} (a subset of a frozen library, e.g. the v7 entries no scripted step writes; md5 checked against both)"""
    st = X.lc["start"]
    sub = None
    if isinstance(st, dict):
        name = st["frozen"]
        sub = {r["file"]: (r["md5"], r.get("reason", "")) for r in _tsv(_path(X.cfg, X.bdir, st["files"])).to_dict("records")}
    else:
        name = st
    fi = X.cfg["frozen_inputs"][name]
    src = fi["path"]
    man = read_any_manifest(fi["manifest"])
    if sub is None:
        have = sorted(f for f in os.listdir(src) if os.path.isfile(os.path.join(src, f)))
        listed = sorted(f for f, _ in man)
        if have != listed:
            raise RuntimeError(f"library start {src}: folder and manifest differ ({len(set(have) ^ set(listed))} files)")
    else:
        mm = dict(man)
        bad = [f for f, (m, _) in sub.items() if mm.get(f) != m]
        if bad:
            raise RuntimeError(f"library start {st['files']}: {len(bad)} files not in {fi['manifest']} with that md5 ({bad[:5]})")
        man = [(f, m) for f, (m, _) in sorted(sub.items())]
    bad = []
    for f, m in man:
        if _md5(os.path.join(src, f)) != m:
            bad.append(f)
            continue
        shutil.copy2(os.path.join(src, f), os.path.join(X.partial, f))
        X.origin[f] = "start"
    if bad:
        raise RuntimeError(f"library start {src}: {len(bad)} files differ from {fi['manifest']} (e.g. {bad[:5]})")
    X.results["start"] = dict(frozen_input=name, files=len(man), manifest=fi["manifest"], manifest_md5=_md5(fi["manifest"]))
    if sub is not None:
        X.results["start"].update(subset=st["files"], subset_md5=_md5(_path(X.cfg, X.bdir, st["files"])))
    log(f"start: {len(man)} frozen entries copied from {src}")


# origin: yeast scripts/s0_library.py:main (complexes branch) @ e3b7448
def step_complexes(X, p):
    """curated complexes (species rows kind == complex) from their RCSB assembly 1; membrane complexes oriented by the
    method of the organism's table (orient, membrane, tm_uniprot)"""
    from .builders import build_complex, parse_ranges
    sp = species(X)
    tab = {r["ing_id"]: r for r in _tsv(_path(X.cfg, X.bdir, p["table"])).to_dict("records")} if p.get("table") else {}
    per = int(p.get("per", 40))
    for r in sp[sp.kind == "complex"].itertuples():
        o = tab.get(r.ing_id, {})
        tm = []
        if o.get("tm_uniprot"):
            tm = parse_ranges(topology(X).loc[o["tm_uniprot"], "tm"], "TRANSMEM")
        n, ns, rmax = build_complex(X.partial, r.ing_id, r.pdb, X.pdb_dir, orient=o.get("orient", "none"),
                                    membrane=o.get("membrane", ""), tm_ranges=tm, per=per)
        X.origin[f"{r.ing_id}.npz"] = "complexes"
        log(f"complex {r.ing_id:16s} {r.pdb}: {n} beads, {ns} spheres, radius {rmax:.1f} nm")


# origin: yeast scripts/s0_library.py:main (af branch) @ e3b7448
def step_af_models(X, p):
    """one AlphaFold DB model per protein species (fragment F1); writes library_af.json (written, failed, missing)"""
    from .builders import build_all_af
    sp = species(X)
    want = {r.uniprot: r.ing_id for r in sp[(sp.kind == "protein") & (sp.uniprot != "")].itertuples()}
    tar = _path(X.cfg, X.bdir, X.lc["afdb_tar"])
    res, missing = build_all_af(X.partial, want, topology(X), tar, nproc=X.nproc, plddt_min=float(p.get("plddt_min", 50)),
                                min_kept=int(p.get("min_kept", 10)), per=int(p.get("per", 30)), log=log)
    err = {k: v[1] for k, v in res.items() if v[0] == "ERR"}
    _mark(X, [f"{k}.npz" for k in res if k not in err], "af_models")
    log(f"AlphaFold models: {len(res) - len(err)} written, {len(err)} failed, {len(missing)} proteins without a model in the tar")
    for k, v in list(err.items())[:20]:
        log("  failed", k, v)
    json.dump(dict(written=len(res) - len(err), failed=err, missing=sorted(missing)),
              open(os.path.join(X.out, "library_af.json"), "w"), indent=0)
    X.results["af_models"] = dict(tar=tar, written=len(res) - len(err), failed=len(err), missing=len(missing))


# origin: yeast scripts/s0_library.py:main (oligo branch, LEGACY_OLIGO_PDB) @ e3b7448
def step_pdb_oligomers(X, p):
    """curated PDB homo-oligomers (table: gene, pdb, reason), written as <gene>_oligo for genes that are species"""
    from .builders import build_pdb_protein
    have = set(species(X).ing_id)
    suffix = p.get("suffix", "_oligo")
    for r in _tsv(_path(X.cfg, X.bdir, p["table"])).itertuples():
        if r.gene not in have:
            continue
        n, ns, rmax = build_pdb_protein(X.partial, r.gene + suffix, r.pdb, X.pdb_dir, per=int(p.get("per", 30)))
        X.origin[f"{r.gene}{suffix}.npz"] = "pdb_oligomers"
        log(f"homo-oligomer {r.gene + suffix:12s} {r.pdb}: {n} beads, {ns} spheres")


# origin: yeast scripts/s0_library.py:main (trna branch) @ e3b7448
def step_pdb_rna(X, p):
    """reference RNA structures (entries: [{id, pdb, per}]; yeast: tRNA-Phe 1EHZ)"""
    from .builders import build_pdb_rna
    for e in p["entries"]:
        n, ns, rmax = build_pdb_rna(X.partial, e["id"], e["pdb"], X.pdb_dir, per=int(e.get("per", 15)))
        X.origin[f"{e['id']}.npz"] = "pdb_rna"
        log(f"{e['id']} ({e['pdb']}): {n} beads, {ns} spheres")


def _jurkat(name):
    """Jurkat builder steps (B3; wholecell/library/jurkat.py, legacy workers in wholecell/library/legacy/)"""
    def f(X, p):
        from . import jurkat
        return jurkat.STEPS[name](X, p)
    f.__name__ = f"step_{name}"
    f.__doc__ = f"Jurkat builder step '{name}' (wholecell/library/jurkat.py)"
    return f


step_jurkat_assemblies = _jurkat("jurkat_assemblies")   # option tm_graft: AF grafting of missing TM segments
step_jurkat_monomers = _jurkat("jurkat_monomers")
step_reorient_tm = _jurkat("reorient_tm")
step_tm_route = _jurkat("tm_route")
step_adaptive_trees = _jurkat("adaptive_trees")         # sphere trees refined until volume <= 2.5x the bead volume


def step_pm_topology_fix(X, p):
    """PM topology fix (freeze v8 E7): rotate the listed units so the ectodomain is exoplasmic (wholecell/library/topology.py)"""
    from . import topology
    return topology.step_pm_topology_fix(X, p)


def curated_pdb(X):
    """id -> PDB of curated entries (complexes of the species table, reference RNA): structure source 'PDB'"""
    out = {}
    if X.species_path and os.path.exists(X.species_path):
        sp = species(X)
        out.update({r.ing_id: r.pdb for r in sp[sp.kind == "complex"].itertuples() if r.pdb})
    for s in X.lc.get("steps", []):
        if s.get("step") == "pdb_rna":
            out.update({e["id"]: e["pdb"] for e in s["entries"]})
    return out


def step_seal(X, p=None):
    """lib.partial -> lib (files read-only), library_manifest.tsv (species, path, md5, origin), library_meta.json"""
    rows = []
    for f in sorted(os.listdir(X.partial)):
        fp = os.path.join(X.partial, f)
        rows.append((f, _md5(fp), X.origin.get(f, "unknown")))
    os.replace(X.partial, X.lib)
    for f, _, _ in rows:
        fp = os.path.join(X.lib, f)
        os.chmod(fp, os.stat(fp).st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
    by = {}
    for _, _, o in rows:
        by[o] = by.get(o, 0) + 1
    if X.get("_legacy"):                                 # inputs read by the Jurkat steps (md5 / size)
        from . import jurkat
        X.results["inputs"] = jurkat.inputs_record(X)
    meta = dict(mode="build", lib=X.lib, files=len(rows), species_npz=sum(1 for f, _, _ in rows if f.endswith(".npz")),
                by_origin=by, steps=[s.get("step") for s in X.lc.get("steps", [])], species=X.species_path,
                species_md5=_md5(X.species_path) if X.species_path and os.path.exists(X.species_path) else None,
                pdb_dir=X.pdb_dir, nproc=X.nproc, results=X.results)
    json.dump(meta, open(os.path.join(X.out, "library_meta.json"), "w"), indent=1)
    tmp = os.path.join(X.out, f"library_manifest.tsv.tmp{os.getpid()}")
    with open(tmp, "w") as fh:
        fh.write("species\tpath\tmd5\torigin\n" + "".join(f"{f.rsplit('.', 1)[0]}\t{os.path.join(X.lib, f)}\t{m}\t{o}\n"
                                                         for f, m, o in rows))
    os.replace(tmp, os.path.join(X.out, "library_manifest.tsv"))
    log(f"library sealed: {len(rows)} files in {X.lib} ({by})")


def entries(lib, extra_dirs=()):
    """{id: npz path} of a library folder (+ read-only extra folders searched after it)"""
    out = {}
    for d in list(extra_dirs)[::-1] + [lib]:
        for f in os.listdir(d):
            if f.endswith(".npz"):
                out[f[:-4]] = os.path.join(d, f)
    return out


def step_pack_shapes(X, p=None):
    """derived packing shapes (§2b P8): pack_shapes/pack_trees_<per>.npz + pack_meta.npz from the sealed library"""
    import numpy as np
    from . import pack_shapes as ps
    p = p or X.lc.get("pack_shapes") or {}
    if p.get("enabled", True) is False:
        return
    ids = None
    if p.get("ids"):
        f = _path(X.cfg, X.bdir, p["ids"])
        ids = [str(i) for i in np.load(f)] if f.endswith(".npy") else [ln.strip() for ln in open(f) if ln.strip()]
    ent = entries(X.lib, [_path(X.cfg, X.bdir, d) for d in p.get("extra_dirs", [])])
    ps.run(ent, os.path.join(X.out, "pack_shapes"), per=int(p.get("per", 15)), ids=ids, curated_pdb=curated_pdb(X),
           nproc=X.nproc, log=log, max_beads=p.get("max_beads"))


STEPS = dict(start=step_start, complexes=step_complexes, af_models=step_af_models, pdb_oligomers=step_pdb_oligomers,
             pdb_rna=step_pdb_rna, jurkat_assemblies=step_jurkat_assemblies, jurkat_monomers=step_jurkat_monomers,
             reorient_tm=step_reorient_tm, tm_route=step_tm_route, adaptive_trees=step_adaptive_trees,
             pm_topology_fix=step_pm_topology_fix,
             seal=step_seal, pack_shapes=step_pack_shapes)


def run(cfg, bdir, step=None, nproc=None, **kw):
    X = context(cfg, bdir, nproc)
    man = os.path.join(X.out, "library_manifest.tsv")
    if step in ("pack_shapes",):
        if not os.path.exists(man):
            raise RuntimeError("library pack_shapes: the library is not sealed yet (run the library stage first)")
        return step_pack_shapes(X)
    if os.path.exists(man):
        print("library manifest exists, step skipped", flush=True)
        return
    if os.path.exists(X.lib):
        raise RuntimeError(f"library {X.lib} exists without a manifest: a library folder is never edited in place; "
                           f"remove it (incomplete build) or start a new build")
    if os.path.exists(X.partial):
        raise RuntimeError(f"{X.partial} exists (interrupted run): remove it and rerun")
    os.makedirs(X.partial)
    t0 = time.time()
    seq = [("start", {})] + [(s["step"], s) for s in X.lc.get("steps", [])]
    for name, p in seq:
        if name not in STEPS:
            raise ValueError(f"library.steps: unknown step '{name}' (known: {', '.join(STEPS)})")
        t = time.time()
        STEPS[name](X, p)
        X.results.setdefault("timing_s", {})[name] = round(time.time() - t, 1)
    step_seal(X)
    step_pack_shapes(X)
    log(f"library build done in {time.time() - t0:.0f} s")
