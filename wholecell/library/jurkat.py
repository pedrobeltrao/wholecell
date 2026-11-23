"""Jurkat library builder steps (library mode `build`, unification B3): the five legacy scripts that wrote the v7
library in place (PLANNING_NOTES §7, §7b, v5 R3) as core steps, one function each, run once in this order on the
build's lib.partial/ (never on the frozen v7 library):

    jurkat_assemblies  jlib_build.assemblies with the final code of jlib_rebuild_tm (SIFTS SEQRES mapping; option
                       tm_graft: AlphaFold TM grafting, tm_graft.graft_missing_tm) -> every non-protein species
    jurkat_monomers    jlib_build.monomers (build_structs.process_af on the AFDB v6 F1 model of every protein species)
    reorient_tm        jlib_reorient_tm (1-2-pass TM monomers: soluble blocks out of the slab, incl. the matrix fix)
    tm_route           jlib_tm_route build (membrane libraries of the TM-routed species: orient / relaxed / AFDB
                       fragment / anchored construction)
    adaptive_trees     lib_refine_st with INCLUDE_MEMBRANE=1 over every species (per 30 -> 15, 8, 5 until the tree
                       volume is <= 2.5x the protein volume; the per-30 tree kept as st_c30 / st_r30)

The worker functions are the legacy code in wholecell/library/legacy/ (unchanged numerics); this module only sets
their module paths to the build (lib.partial/, declared inputs), makes them offline (a missing PDB assembly or
UniProt cache entry is an error, nothing is downloaded) and replaces their `main` drivers (job lists from the species
table, logs into <library>/work/). Running all steps once with the final code equals the v7 history: the legacy order
(refine before reorient / tm_route, then refine again for the rebuilt ids) leaves the same final trees, because every
rebuild dropped st_c30 and the refinement starts from the per-30 tree.

Config (library.steps entries; inputs are frozen-input names or paths):
    - {step: jurkat_assemblies, tm_graft: true, nproc: <n>}
    - {step: jurkat_monomers}
    - {step: reorient_tm}
    - {step: tm_route, ids: <TSV ing_id, group>}
    - {step: adaptive_trees, include_membrane: true, max_ratio: 2.5, pers: [15, 8, 5]}
    library.jurkat: {species: <species TSV>, sifts, uniprot_topology, uniprot_table, afdb_tar, afdb_tar_index,
                     afdb_fragment_index, uniprot_rest_cache, afdb_complex_models (folder; manifest md5-checked)}
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import shutil
import sys
import tarfile
import time

import numpy as np

LEGACY = os.path.join(os.path.dirname(os.path.abspath(__file__)), "legacy")
_MODS = ("config", "structlib", "build_structs", "membrane_orient", "tm_graft", "jlib_build", "jlib_reorient_tm",
         "jlib_tm_route", "lib_refine_st")


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def _limit_threads():
    """one BLAS / OpenMP thread per worker process (KMeans); results do not depend on it"""
    try:
        from threadpoolctl import threadpool_limits
        threadpool_limits(1)
    except Exception:  # noqa: BLE001
        pass


def _pool(n):
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing as mp
    return ProcessPoolExecutor(n, mp_context=mp.get_context("fork"), initializer=_limit_threads)


def _snap(d):
    return {f: (os.stat(os.path.join(d, f)).st_mtime_ns, os.path.getsize(os.path.join(d, f))) for f in os.listdir(d)}


def _written(X, before, step):
    after = _snap(X.partial)
    new = [f for f, v in after.items() if before.get(f) != v and not f.startswith(".")]
    for f in new:
        X.origin[f] = step
    return new


# ------------------------------------------------------------------ legacy modules bound to this build
def legacy(X):
    """import the legacy builder modules with their paths set to this build (once per process)"""
    if X.get("_legacy"):
        return X["_legacy"]
    from .build import _path
    jc = X.lc.get("jurkat") or {}
    P = lambda k: _path(X.cfg, X.bdir, jc[k]) if jc.get(k) else None   # noqa: E731
    work = os.path.join(X.out, "work")
    data = os.path.join(work, "data")                                  # legacy DATA folder: the two UniProt tables
    os.makedirs(data, exist_ok=True)
    for src, name in ((P("uniprot_topology"), "uniprot_topology.tsv"), (P("uniprot_table"), "uniprot_yeast.tsv")):
        if src and not os.path.exists(os.path.join(data, name)):
            shutil.copy2(src, os.path.join(data, name))
    os.environ.update(WHOLECELL_LIB_OUT=X.partial, WHOLECELL_LIB_DATA=data, WHOLECELL_AFDB_HUMAN_TAR=P("afdb_tar") or "",
                      WHOLECELL_PDB=X.pdb_dir or "", INCLUDE_MEMBRANE="1")
    for m in _MODS:                                                    # bare-name imports: no foreign module of the name
        if m in sys.modules and not getattr(sys.modules[m], "__file__", "").startswith(LEGACY):
            raise RuntimeError(f"module '{m}' already imported from {sys.modules[m].__file__}; legacy builders need theirs")
    if LEGACY not in sys.path:
        sys.path.insert(0, LEGACY)
    import config as lc_, structlib as sl, build_structs as bs, membrane_orient as mo, tm_graft as tg  # noqa: F401
    import jlib_build as jb, jlib_reorient_tm as jr, jlib_tm_route as jt, lib_refine_st as rs
    lc_.LIB = bs.LIB = jb.LIB = jr.LIB = rs.LIB = jt.LIB = jt.OUT = X.partial
    lc_.DATA = mo.DATA = data
    mo._TOPO = None; mo._BARREL = None
    jb.PDBDIR = X.pdb_dir
    jb.AFDIR = P("afdb_complex_models")
    jb.TAR = jr.TAR = jt.TAR = tg.TAR = P("afdb_tar")
    tg.IDX = P("afdb_tar_index"); tg._IDX = None
    jt.UPJ = P("uniprot_rest_cache"); jt.FRAG = P("afdb_fragment_index"); jt._UP = None
    jt.BAK = os.path.join(work, "tm_route_before")                     # pre-tm_route copies (legacy lib_bak_v4_tmroute)
    jt.D5 = jt.LOG = jt.IDS = jr.V = jb.V = rs.V = jt.V = "/nonexistent-legacy-root"   # nothing may read the legacy tree
    used_pdb = {}

    def pdb_file(pid):                                                 # offline: the assembly must be in the cache
        fn = f"{X.pdb_dir}/{pid}-assembly1.cif.gz"
        if not os.path.exists(fn):
            raise FileNotFoundError(f"PDB assembly {fn} not in the cache (offline build)")
        used_pdb[pid] = fn
        return fn
    jb.pdb_file = pdb_file
    S = None

    def sifts():
        nonlocal S
        if S is None:
            import pandas as pd
            S = pd.read_csv(P("sifts"), sep="\t", comment="#", usecols=["PDB", "CHAIN", "SP_PRIMARY", "RES_BEG", "RES_END",
                                                                            "SP_BEG", "SP_END"], low_memory=False)
        return S
    jb.sifts = sifts

    def fetch_uniprot(accs):                                           # offline: the REST cache is a frozen input
        have = json.load(open(jt.UPJ))
        miss = sorted(set(accs) - set(have) - {""})
        if miss:
            raise RuntimeError(f"tm_route: {len(miss)} accessions not in the UniProt cache {jt.UPJ} (e.g. {miss[:5]})")
    jt.fetch_uniprot = fetch_uniprot
    if jc.get("afdb_complex_models_manifest"):                         # AFDB complex models: md5-checked when read
        from .build import read_any_manifest
        X["_afc_md5"] = dict(read_any_manifest(_path(X.cfg, X.bdir, jc["afdb_complex_models_manifest"])))
    L = dict(sl=sl, bs=bs, mo=mo, tg=tg, jb=jb, jr=jr, jt=jt, rs=rs, used_pdb=used_pdb, work=work, P=P)
    X["_legacy"] = L
    return L


def table(X):
    """species table (ing_id, kind, genes, accs, stoich, struct_src, pdb, uniprot, transmembrane, signal_frac, mass_da)"""
    if "_jtable" not in X:
        import pandas as pd
        from .build import _path
        T = pd.read_csv(_path(X.cfg, X.bdir, X.lc["species"]), sep="\t", keep_default_na=False)
        T["transmembrane"] = T.transmembrane.astype(str).str.lower().isin(["true", "1"])
        for c in ("signal_frac", "mass_da"):
            T[c] = pd.to_numeric(T[c])
        X["_jtable"] = T
    return X["_jtable"]


def _accs(r):
    return r["accs"].split(";") if r["accs"] else [""] * len(r["genes"].split(";"))


def _check_afc(X, L, rows):
    want = X.get("_afc_md5")
    if want is None:
        return
    bad = []
    for r in rows:
        if r["struct_src"].startswith("afdb_complex:"):
            f = r["struct_src"].split(":", 1)[1] + ".cif.gz"
            if want.get(f) != _md5(os.path.join(L["jb"].AFDIR, f)):
                bad.append(f)
    if bad:
        raise RuntimeError(f"AFDB complex models differ from their manifest: {bad[:5]}")


# ------------------------------------------------------------------ steps
# origin: jurkat/yc/pipeline/jlib_build.py:assemblies + jlib_rebuild_tm.py (final code) @ legacy 2026-11-15
def step_jurkat_assemblies(X, p):
    """every non-protein species from its structure (struct_src pdb:<id> | afdb_complex:<model>); TM assemblies
    oriented by membrane_orient (option tm_graft: AF grafting of missing TM segments first); NPC ring axis"""
    L = legacy(X)
    jb, tg = L["jb"], L["tg"]
    T = table(X)
    A = T[T.kind != "protein"]
    rows = [dict(ing_id=r.ing_id, struct_src=r.struct_src, transmembrane=bool(r.transmembrane), accs=_accs(r._asdict()))
            for r in A.itertuples()]
    _check_afc(X, L, rows)
    graft = bool(p.get("tm_graft", True))
    orig = tg.graft_missing_tm
    if not graft:
        tg.graft_missing_tm = lambda xyz, kind, acc, pos, chn: (xyz, kind, acc, pos, "")
    for pid in sorted({r["struct_src"].split(":")[1] for r in rows if r["struct_src"].startswith("pdb:")}):
        jb.pdb_file(pid)
    jb.sifts()                                                          # loaded once, before the fork
    before = _snap(X.partial)
    n = int(p.get("nproc", X.nproc))
    order = sorted(range(len(rows)), key=lambda i: -_size(X, rows[i]))  # big assemblies first (NPC), same outputs
    try:
        with _pool(n) as ex:
            res = dict(zip(order, ex.map(jb.build_assembly, [rows[i] for i in order], chunksize=1)))
    finally:
        tg.graft_missing_tm = orig
    import pandas as pd
    R = pd.DataFrame([res[i] for i in range(len(rows))], columns=["ing_id", "status", "note"])
    R.to_csv(os.path.join(L["work"], "assemblies_log.tsv"), sep="\t", index=False)
    new = _written(X, before, "jurkat_assemblies")
    nt = R.note.fillna("")
    X.results["jurkat_assemblies"] = dict(rows=len(R), status=R.status.value_counts().to_dict(), written=len(new),
                                          oriented=int(nt.str.contains("oriented").sum()),
                                          orient_failed=int(nt.str.contains("orient_failed").sum()),
                                          grafted=int(nt.str.contains(r"graft: [1-9]").sum()), tm_graft=graft)
    log(f"assemblies: {X.results['jurkat_assemblies']}")


def _size(X, r):
    s = r["struct_src"]
    try:
        if s.startswith("pdb:"):
            return os.path.getsize(f"{X.pdb_dir}/{s[4:]}-assembly1.cif.gz")
    except OSError:
        pass
    return 0


# origin: jurkat/yc/pipeline/jlib_build.py:monomers + build_structs.process_all_af @ legacy 2026-11-15
# (the tar is read in chunks: same jobs, bounded memory; fork pool instead of joblib)
def step_jurkat_monomers(X, p):
    """one AFDB v6 F1 model per protein species (pLDDT < 50 trimmed, TRANSMEM kept; TM monomers oriented on the helix
    axes, cytoplasmic side +z from UniProt TOPO_DOM)"""
    import pandas as pd
    L = legacy(X)
    bs = L["bs"]
    T = table(X)
    P = T[T.kind == "protein"]
    topo = pd.read_csv(os.path.join(L["work"], "data", "uniprot_topology.tsv"), sep="\t").set_index("Entry")
    topo = topo.rename(columns={"Transmembrane": "tm", "Topological domain": "topo"})
    want = {r.uniprot: r.ing_id for r in P[P.uniprot != ""].itertuples()}
    before = _snap(X.partial)
    res, seen, jobs = {}, set(), []
    chunk = int(p.get("chunk", 600))

    def flush(ex):
        res.update(dict(ex.map(bs._af_worker, jobs, chunksize=8)))
        jobs.clear()
    with _pool(int(p.get("nproc", X.nproc))) as ex, tarfile.open(L["P"]("afdb_tar")) as tf:
        for m in tf:
            name = os.path.basename(m.name)
            if not name.endswith(".cif.gz") or "-F1-" not in name:
                continue
            acc = name.split("-")[1]
            if acc in want:
                data = gzip.decompress(tf.extractfile(m).read()).decode()
                tm = bs.parse_ranges(topo.loc[acc, "tm"], "TRANSMEM") if acc in topo.index else []
                tp = bs.parse_ranges(topo.loc[acc, "topo"], "TOPO_DOM") if acc in topo.index else []
                jobs.append((acc, want[acc], name, data, tm, tp)); seen.add(want[acc])
                if len(jobs) >= chunk:
                    flush(ex)
                    log(f"  monomers: {len(res)} done")
        flush(ex)
    missing = sorted(set(want.values()) - seen) + sorted(P.ing_id[P.uniprot == ""])
    rows = [(k, "ok" if not (isinstance(v, tuple) and v and v[0] == "ERR") else "error", str(v)[:150]) for k, v in res.items()]
    rows += [(m, "no_afdb_model", "") for m in missing]
    R = pd.DataFrame(rows, columns=["ing_id", "status", "note"])
    R.to_csv(os.path.join(L["work"], "monomers_log.tsv"), sep="\t", index=False)
    new = _written(X, before, "jurkat_monomers")
    X.results["jurkat_monomers"] = dict(proteins=len(P), status=R.status.value_counts().to_dict(), written=len(new))
    log(f"monomers: {X.results['jurkat_monomers']}")


# origin: jurkat/yc/pipeline/jlib_reorient_tm.py:main @ legacy 2026-11-15
def step_reorient_tm(X, p):
    """1-2-pass TM monomers (UniProt TRANSMEM) rebuilt with the soluble blocks (>= 30 residues) moved to |z| >= 2.3 nm"""
    import pandas as pd
    L = legacy(X)
    bs, jr = L["bs"], L["jr"]
    T = table(X)
    I = T[(T.kind == "protein") & T.transmembrane & (T.uniprot != "")]
    topo = pd.read_csv(os.path.join(L["work"], "data", "uniprot_topology.tsv"), sep="\t").set_index("Entry")
    want = {}
    for r in I.itertuples():
        if r.uniprot not in topo.index or not os.path.exists(f"{X.partial}/{r.ing_id}.npz"):
            continue
        tmr = bs.parse_ranges(topo.loc[r.uniprot, "Transmembrane"], "TRANSMEM")
        if 1 <= len(tmr) <= 2:
            want[r.uniprot] = (r.ing_id, tmr, bs.parse_ranges(topo.loc[r.uniprot, "Topological domain"], "TOPO_DOM"),
                               r.signal_frac >= 0.5)
    jobs = []
    with tarfile.open(L["P"]("afdb_tar")) as tf:
        for m in tf:
            name = os.path.basename(m.name)
            if name.endswith(".cif.gz") and "-F1-" in name and name.split("-")[1] in want:
                acc = name.split("-")[1]; ing, tmr, tp, sig = want[acc]
                jobs.append((acc, ing, gzip.decompress(tf.extractfile(m).read()).decode(), tmr, tp, sig))
    before = _snap(X.partial)
    with _pool(int(p.get("nproc", X.nproc))) as ex:
        res = list(ex.map(jr.rebuild, jobs, chunksize=4))
    R = pd.DataFrame(res, columns=["ing_id", "status", "in_slab_soluble_before", "in_slab_soluble_after", "beads_moved"])
    R.to_csv(os.path.join(L["work"], "reorient_tm_log.tsv"), sep="\t", index=False)
    new = _written(X, before, "reorient_tm")
    X.results["reorient_tm"] = dict(candidates=len(want), models=len(jobs), status=R.status.value_counts().to_dict(),
                                    in_slab_before=int(R.in_slab_soluble_before.sum()),
                                    in_slab_after=int(R.in_slab_soluble_after.sum()),
                                    moved=int((R.beads_moved > 0).sum()), written=len(new))
    log(f"reorient_tm: {X.results['reorient_tm']}")


# origin: jurkat/yc/pipeline/jlib_tm_route.py:main_build @ legacy 2026-11-15
def step_tm_route(X, p):
    """membrane libraries of the TM-routed species (ids TSV: ing_id, group): orient | orient_relaxed | af_fragment |
    anchored (see legacy/jlib_tm_route.py)"""
    import pandas as pd
    from .build import _path
    L = legacy(X)
    jt = L["jt"]
    S = table(X).set_index("ing_id")
    ids = pd.read_csv(_path(X.cfg, X.bdir, p["ids"]), sep="\t")
    jobs, allacc = [], set()
    for r in ids.itertuples():
        row = S.loc[r.ing_id].to_dict()
        row["accs"] = _accs(row)
        if row["kind"] == "protein" and row["uniprot"]:
            row["accs"] = [row["uniprot"]]
        allacc |= set(row["accs"]); jobs.append((r.ing_id, r.group, row))
    jt.fetch_uniprot(allacc)
    _check_afc(X, L, [j[2] for j in jobs])
    for row in (j[2] for j in jobs):
        if row["struct_src"].startswith("pdb:"):
            L["jb"].pdb_file(row["struct_src"].split(":")[1])
    L["jb"].sifts()
    before = _snap(X.partial)
    with _pool(int(p.get("nproc", X.nproc))) as ex:
        res = list(ex.map(jt.build, jobs, chunksize=1))
    R = pd.DataFrame(res)
    R.to_csv(os.path.join(L["work"], "tm_route_libs.tsv"), sep="\t", index=False)
    new = _written(X, before, "tm_route")
    err = R[R.note.astype(str).str.startswith("ERROR")]
    X.results["tm_route"] = dict(ids=len(R), method=R.method.value_counts().to_dict(), errors=len(err), written=len(new))
    log(f"tm_route: {X.results['tm_route']}")
    for r in err.itertuples():
        log("  ERROR", r.ing_id, r.note[:200])


# origin: jurkat/yc/pipeline/lib_refine_st.py:main (INCLUDE_MEMBRANE=1) @ legacy 2026-11-15
def step_adaptive_trees(X, p):
    """adaptive sphere trees: per 30 kept if its volume <= max_ratio x 1.21 nm3/kDa x mass, else per 15, 8, 5 (the
    first that meets the bound, else the finest); species with > 20,000 beads are skipped (compact assemblies)"""
    import pandas as pd
    L = legacy(X)
    rs = L["rs"]
    rs.MAX_RATIO = float(p.get("max_ratio", 2.5)); rs.PERS = tuple(int(x) for x in p.get("pers", (15, 8, 5)))
    if p.get("include_membrane", True):
        os.environ["INCLUDE_MEMBRANE"] = "1"
    else:
        os.environ.pop("INCLUDE_MEMBRANE", None)
    T = table(X)
    T = T[T.ing_id.map(lambda i: os.path.exists(f"{X.partial}/{i}.npz"))]
    jobs = list(zip(T.ing_id, T.mass_da / 1e3))
    before = _snap(X.partial)
    with _pool(int(p.get("nproc", X.nproc))) as ex:
        rows = list(ex.map(rs.refine, jobs, chunksize=16))
    R = pd.DataFrame(rows, columns=["ing_id", "kda", "ratio30", "ratio_new", "nsph30", "nsph_new", "action"])
    R.to_csv(os.path.join(L["work"], "adaptive_trees_log.tsv"), sep="\t", index=False)
    new = _written(X, before, "adaptive_trees")
    ok = R.ratio30.notna()
    X.results["adaptive_trees"] = dict(species=len(R), action=R.action.value_counts().to_dict(), written=len(new),
                                       ratio_median=[round(float(R.ratio30[ok].median()), 3), round(float(R.ratio_new[ok].median()), 3)],
                                       spheres=[int(R.nsph30.sum()), int(R.nsph_new.sum())])
    log(f"adaptive_trees: {X.results['adaptive_trees']}")


def inputs_record(X):
    """md5 / size of the inputs read by the Jurkat steps (library_meta.json results.inputs)"""
    L = X.get("_legacy")
    if not L:
        return None
    out = {}
    for k in ("afdb_tar", "afdb_tar_index", "afdb_fragment_index", "uniprot_rest_cache", "sifts", "uniprot_topology",
              "uniprot_table"):
        f = L["P"](k)
        if f and os.path.exists(f):
            out[k] = dict(path=f, bytes=os.path.getsize(f), md5=_md5(f) if os.path.getsize(f) < 2 << 30 else "frozen input (verified at build start)")
    out["species"] = dict(path=X.species_path, md5=_md5(X.species_path))
    pdb = L["used_pdb"]
    with open(os.path.join(L["work"], "pdb_assemblies_read.tsv"), "w") as fh:
        fh.write("pdb\tmd5\n" + "".join(f"{k}\t{_md5(v)}\n" for k, v in sorted(pdb.items())))
    out["pdb_assemblies"] = dict(n=len(pdb), list=os.path.join(L["work"], "pdb_assemblies_read.tsv"))
    return out


STEPS = dict(jurkat_assemblies=step_jurkat_assemblies, jurkat_monomers=step_jurkat_monomers, reorient_tm=step_reorient_tm,
             tm_route=step_tm_route, adaptive_trees=step_adaptive_trees)
