#!/usr/bin/env python
"""A3 check: the core library stage (library.mode build, config/cells/yeast.yaml `library:`) against the yeast
port build s0 (criteria L1-L5 of docs/reference/history/yeast_x5_reproduction_criteria.md with the 2026-11-17 addendum:
L5 ignores the `sig` strings of the pack trees).

    <python> scripts/library_a3_check.py --label yeastenv|vcell [--nproc 3] [--out DIR] [--s0 BUILD] [--smoke]

1. species table: <out>/species/library_species.tsv (plugins/yeast/library_species.py, yeast env; made once)
2. the core stage runs into <out>/run_<label>/library/ (a new build-like folder; inputs read in place, offline)
3. comparison with s0 (READ-ONLY): engine/lib, engine/work/results/library_af.json, engine/work/yc3d/{lib_ids.npy,
   libpack_meta.npz, pack_trees_15.npz}; written to <out>/run_<label>/check_report.json.
The s0 library holds stage-0 entries and entries written later by engine stages 1 and 3 (their mtime is after the s0
step ended, from the s0 log). L1/L2 compare the stage-0 part; L5 builds the pack trees of s0's 5,530 packed types from
the core stage's pack_shapes for entries of the canonical library and with the core pack_shapes code from the s0
files for entries the engine wrote later.
"""
import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import sys
import time

os.environ.setdefault("PYTHONHASHSEED", "0")
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
sys.path.insert(0, os.path.join(REPO, "plugins", "yeast", "repro"))

S0 = "{results_root}/yeast/20261116-aabe7b5-s0-x5repro"
OUT = "{results_root}/yeast/_dev_unify/a3"
# A9: the library block lives in config/cells/yeast.yaml (the A3 fragment was merged); species path redirected here
FRAG = {"library": {"species": "{build}/inventory/library_species.tsv"}}


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def s0_end(s0):
    """epoch of the end of the s0 library step (last log line time, on the date of the log file)"""
    lg = os.path.join(s0, "engine", "work", "s0_library.log")
    last = [ln for ln in open(lg) if re.match(r"\d\d:\d\d:\d\d ", ln)][-1][:8]
    day = datetime.datetime.fromtimestamp(os.path.getmtime(lg)).date()
    t = datetime.datetime.combine(day, datetime.time.fromisoformat(last))
    return t.timestamp() + 1.0, last


def main():
    import numpy as np
    import yaml
    from wholecell import config as C
    from wholecell.library import build as LB, pack_shapes as PS
    import compare_x5 as CX

    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--out", default=None)
    ap.add_argument("--s0", default=None)
    ap.add_argument("--smoke", action="store_true", help="tiny species subset, no comparison of counts")
    ap.add_argument("--compare-only", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    site = C.load_site()
    s0 = a.s0 or S0.replace("{results_root}", site["results_root"])
    out = a.out or OUT.replace("{results_root}", site["results_root"])
    run = os.path.join(out, f"run_{a.label}")
    cfg = C.load("yeast", overrides=FRAG)
    cfg["_site"]["repo"] = REPO
    sp_src = os.path.join(out, "species", "library_species.tsv")
    if not os.path.exists(sp_src):
        raise SystemExit(f"{sp_src} missing: run plugins/yeast/library_species.py in the yeast env first")
    os.makedirs(os.path.join(run, "inventory"), exist_ok=True)
    sp_dst = os.path.join(run, "inventory", "library_species.tsv")
    if a.smoke:
        import pandas as pd
        sp = pd.read_csv(sp_src, sep="\t", dtype=str, keep_default_na=False)
        keep = sp[sp.ing_id.isin(["PolII", "nucleosome", "Pma1_hexamer", "V-ATPase", "TDH3", "PMA1", "TOM40", "SEC61"])
                  | sp.index.isin(range(20, 40))]
        keep.to_csv(sp_dst, sep="\t", index=False)
    elif not os.path.exists(sp_dst):
        shutil.copy2(sp_src, sp_dst)
    import numpy, sklearn
    env = dict(python=sys.executable, numpy=numpy.__version__, sklearn=sklearn.__version__, nproc=a.nproc,
               omp=os.environ.get("OMP_NUM_THREADS"))
    if not a.compare_only:
        LB.run(cfg, run, nproc=a.nproc)
    t_build = time.time() - t0
    rep = dict(label=a.label, env=env, s0=s0, run=run, species=sp_dst, species_md5=md5(sp_dst), build_s=round(t_build))
    if a.smoke:
        json.dump(rep, open(os.path.join(run, "check_report.json"), "w"), indent=1)
        print(json.dumps(rep, indent=1))
        return

    lib = os.path.join(run, "library", "lib")
    plib = os.path.join(s0, "engine", "lib")
    pyc = os.path.join(s0, "engine", "work", "yc3d")
    tend, tlast = s0_end(s0)
    rep["s0_step_end"] = tlast
    man = {r.split("\t")[0]: r.rstrip("\n").split("\t") for r in list(open(os.path.join(run, "library", "library_manifest.tsv")))[1:]}
    origin = {v[0] + "." + v[1].rsplit(".", 1)[1]: v[3] for v in man.values()}
    canon = set(os.listdir(lib))
    port = set(os.listdir(plib))
    late = {f for f in port if os.path.getmtime(os.path.join(plib, f)) > tend}
    s0part = port - late
    # ---- L1 file set (stage-0 part of s0)
    start = dict((f, m) for f, m in LB.read_any_manifest(cfg["frozen_inputs"]["yeast_x5_lib_start"]["manifest"]))
    l1 = dict(canonical=len(canon), canonical_npz=sum(f.endswith(".npz") for f in canon), s0_lib=len(port),
              s0_stage0_files=len(s0part), s0_later_files=len(late),
              missing_in_canonical=sorted(s0part - canon)[:50], n_missing=len(s0part - canon),
              not_in_s0=sorted(canon - port)[:50], n_not_in_s0=len(canon - port),
              canonical_rewritten_later_by_engine=sorted(canon & late)[:50], n_rewritten_later=len(canon & late))
    lg = open(os.path.join(s0, "engine", "work", "s0_library.log")).read()
    m = re.search(r"library: (\d+) entries", lg)
    l1["s0_log_npz_entries"] = int(m.group(1)) if m else None
    l1["level"] = "E" if (not l1["n_missing"] and not l1["n_not_in_s0"] and l1["canonical_npz"] == l1["s0_log_npz_entries"]) else "fail"
    rep["L1"] = l1
    # ---- L2 rebuilt entries (stage-0 part, start entries excluded)
    names = sorted((canon & s0part) - set(start))
    from multiprocessing import Pool
    with Pool(a.nproc) as pool:
        r = pool.starmap(CX._lib_cmp, [(f"{lib}/{n}", f"{plib}/{n}", n) for n in names], chunksize=64)
    nE = sum(x[1] == "E" for x in r); nT = [x for x in r if x[1] == "T"]; nF = [x for x in r if x[1] == "fail"]
    # size of the differences of T entries
    worst = {}
    for n, _, _ in nT:
        A = CX.npz_dict(f"{lib}/{n}", True); B = CX.npz_dict(f"{plib}/{n}", True)
        for k in A:
            if A[k].dtype.kind in "fc" and not CX.arrays_equal(A[k], B[k]):
                worst[k] = max(worst.get(k, 0.0), float(np.max(np.abs(A[k].astype(float) - B[k].astype(float)))))
    shape_fail = any(x[2] == "keys/shapes" for x in nF)
    lev = "E" if nE == len(r) else ("T" if not shape_fail and (len(r) - len(nF)) / max(len(r), 1) >= CX.L2_SHARE else "fail")
    rep["L2"] = dict(compared=len(r), E=nE, T=len(nT), fail=len(nF), T_ids=[x[0] for x in nT][:100], fail_ids=nF[:50],
                     T_max_abs_diff=worst, level=lev,
                     not_compared_rewritten_later=len((canon & late) - set(start)))
    # ---- L3 library_af.json
    A = json.load(open(os.path.join(run, "library", "library_af.json")))
    B = json.load(open(os.path.join(s0, "engine", "work", "results", "library_af.json")))
    rep["L3"] = dict(core=dict(written=A["written"], failed=len(A["failed"]), missing=len(A["missing"])),
                     s0=dict(written=B["written"], failed=len(B["failed"]), missing=len(B["missing"])),
                     level="E" if A == B else "fail",
                     bytes_identical=open(os.path.join(run, "library", "library_af.json"), "rb").read()
                     == open(os.path.join(s0, "engine", "work", "results", "library_af.json"), "rb").read())
    # ---- L4 start entries
    bad = [f for f, m5 in start.items() if not os.path.exists(os.path.join(lib, f)) or md5(os.path.join(lib, f)) != m5]
    over = sorted(f for f in start if origin.get(f, "start") != "start")
    rep["L4"] = dict(start=len(start), unchanged=len(start) - len(bad), changed=bad[:20], rebuilt_by_core_step=over,
                     level="E" if not bad else "fail")
    # ---- L5 pack shapes for s0's packed types
    ids = [str(x) for x in np.load(os.path.join(pyc, "lib_ids.npy"))]
    z = np.load(os.path.join(run, "library", "pack_shapes", "pack_trees_15.npz"))
    off = z["off"]
    core_trees = {str(i): (z["C"][off[k]:off[k + 1]], z["R"][off[k]:off[k + 1]]) for k, i in enumerate(z["ids"])}
    src_core = [i for i in ids if f"{i}.npz" in canon and f"{i}.npz" not in late]
    src_port = [i for i in ids if i not in set(src_core)]
    ent = {i: os.path.join(lib, f"{i}.npz") for i in src_core}
    ent.update({i: os.path.join(plib, f"{i}.npz") for i in src_port})
    t5 = time.time()
    trees = {i: core_trees[i] for i in src_core}
    trees.update(PS.pack_trees({i: ent[i] for i in src_port}, 15, nproc=a.nproc))
    t5 = time.time() - t5
    sid = sorted(trees)
    offs = np.r_[0, np.cumsum([len(trees[i][1]) for i in sid])]
    mine = dict(ids=np.array(sid), off=offs, C=np.concatenate([trees[i][0] for i in sid]), R=np.concatenate([trees[i][1] for i in sid]))
    meta = PS.pack_meta(ids, ent, trees, LB.curated_pdb(LB.context(cfg, run, a.nproc)))
    res5 = {}
    for name, A5, keys in (("pack_trees_15.npz", mine, ["ids", "off", "C", "R"]), ("libpack_meta.npz", meta, ["ids", "BR", "src", "pdb"])):
        B5 = CX.npz_dict(os.path.join(pyc, name), True)
        ign = sorted(set(B5) - set(keys))
        d = {}
        for k in keys:
            if CX.arrays_equal(np.asarray(A5[k]), B5[k]):
                d[k] = "E"
            elif A5[k].shape == B5[k].shape and A5[k].dtype.kind in "fc":
                md = float(np.max(np.abs(A5[k] - B5[k])))
                d[k] = f"T {md:.3g}" if md <= CX.L2_TREE_NM else f"fail {md:.3g}"
            else:
                d[k] = "fail (shape/non-float)"
        res5[name] = dict(keys=d, ignored=ign)
    if any(v.startswith("fail") for x in res5.values() for v in x["keys"].values()):
        pt = mine; Bt = CX.npz_dict(os.path.join(pyc, "pack_trees_15.npz"), True)
        diff_ids = []
        if CX.arrays_equal(pt["off"], Bt["off"]):
            for k, i in enumerate(sid):
                a_, b_ = slice(pt["off"][k], pt["off"][k + 1]), slice(Bt["off"][k], Bt["off"][k + 1])
                if not (np.array_equal(pt["C"][a_], Bt["C"][b_]) and np.array_equal(pt["R"][a_], Bt["R"][b_])):
                    diff_ids.append(i)
        res5["differing_ids"] = diff_ids[:50]
    ids_ok = CX.arrays_equal(np.array(ids), np.load(os.path.join(pyc, "lib_ids.npy")).astype(str))
    vals = [v for x in res5.values() if isinstance(x, dict) for v in x["keys"].values()]
    rep["L5"] = dict(lib_ids=len(ids), lib_ids_identical=ids_ok, from_canonical=len(src_core), from_s0_later_engine=len(src_port),
                     trees_s=round(t5), files=res5,
                     note="lib_ids (the packed type list) is chosen by engine stage 3 and taken from s0; sig ignored (addendum)",
                     level="E" if ids_ok and all(v == "E" for v in vals) else
                     ("T" if ids_ok and all(v == "E" or v.startswith("T") for v in vals) else "fail"))
    # ---- builder log lines (beads / spheres / radius) vs s0 log
    mylog = os.path.join(out, f"run_{a.label}.log")
    if os.path.exists(mylog):
        pat = re.compile(r"(complex \S+|homo-oligomer \S+|AlphaFold models:.*|tRNA.*?:)\s*(.*)")
        def lines(t):
            return sorted(set(re.sub(r"^\d\d:\d\d:\d\d ", "", ln).strip() for ln in t.splitlines()
                              if re.search(r"^\d\d:\d\d:\d\d (complex |homo-oligomer |AlphaFold models:)", ln)))
        la, lb = lines(open(mylog).read()), lines(lg)
        rep["log_lines"] = dict(core=len(la), s0=len(lb), identical=la == lb, only_core=sorted(set(la) - set(lb))[:10],
                                only_s0=sorted(set(lb) - set(la))[:10])
    rep["verdict"] = {k: rep[k]["level"] for k in ("L1", "L2", "L3", "L4", "L5")}
    rep["wall_s"] = round(time.time() - t0)
    json.dump(rep, open(os.path.join(run, "check_report.json"), "w"), indent=1, default=str)
    print(json.dumps(rep["verdict"]), "->", os.path.join(run, "check_report.json"), flush=True)


if __name__ == "__main__":
    main()
