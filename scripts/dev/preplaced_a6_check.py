"""Check of the serial pre-placed stage (yeast settings) against the yeast port build s0 (unification A6, yeast half).

Runs `wholecell.preplaced.serial` with config/cells/yeast.yaml preplaced.serial (builders in the yeast stage-3 order,
context yeast.s3) on the inputs stage 3 of s0 read, and compares with what s0's stage 3 wrote:
  S1  prepass.npz array key set (ours: preplaced.npz)                identical
  S2  length of every array                                         identical
  S3  s3 log builder lines (29 patterns of compare_x5.s3_lines; timing tokens "<n>s" removed, addendum 2026-11-17),
      plus the whole pre-pass log sequence (every line from "microtubules:" to "mRNP ends:", consume lines included)
  S4  every array (dtype, shape, values) per structure group (compare_x5.PRE_GROUPS)    identical
  E   pools after the pre-pass: wc3d.shell_capacity(ctx.exp3) vs s0 yc3d/expected_totals.parquet (identical)
Level E = identical. Nothing is written into the s0 build or the yeast tree.

  1. prepare (once):  conda_envs/yeast_x5/bin/python scripts/preplaced_a6_check.py prepare
       OUT/engine/{root (data + cellshape files linked), work (results/, rcsb cache, yc3d/ stage-2 files: copies),
       lib (cp -a of the s0 library: unchanged since 18:33, before s3 started), pdb (link farm)}
  2. run + compare:   PYTHONHASHSEED=0 <python> scripts/preplaced_a6_check.py run <tag>    (tag: yeast_x5 | vcell)
       -> OUT/run-<tag>/{preplaced/serial/..., run.log, report.json, report.md}
  3. compare_ref <tag> <ref_dir> [n_s0]: per builder layer, the first n_s0 builders (default 18) against s0 (identical)
       and the rest against the serial pre-pass of another run (<ref_dir>/preplaced/serial; prepass/serial in runs
       before A0c; e.g. the yeast env with a port patch that changes output), at tolerance: per array max |difference| or count difference
       -> OUT/run-<tag>/{compare_ref.json, compare_ref.md}
  Env: A6_OUT = output folder (default _dev_unify/a6), A6_RUN_<tag> = run folder name (default run-<tag>).
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _devsite import VC, RESULTS, YEAST as _YEAST_ROOT, PY as _SITE_PY, SMK as _SITE_SMK, GIT as _DEV_GIT  # noqa: E402,F401  (release: site paths)
S0 = f"{RESULTS}/yeast/20261116-aabe7b5-s0-x5repro"
OUT = os.environ.get("A6_OUT", f"{RESULTS}/yeast/_dev_unify/a6")   # A6_OUT: other output folder
YEAST = _YEAST_ROOT
PIPE = f"{REPO}/plugins/yeast/engine/pipeline"
S3_PATS = [r"microtubules: .*", r"kinetochores: .*", r"outer kinetochores: .*", r"Rabl layout: .*", r"nucleosome map: .*", r"chromatin: \d.*",
           r"pore receptors: .*", r"FG mesh: .*", r"rDNA: .*", r"rDNA chromatin: .*", r"engaged Pol II: .*", r"origins: .*",
           r"chromatin sites: .*", r"tRNA genes, telomeres, silent loci: .*", r"eisosomes: .*", r"MCC proteins .*", r"actin cables: .*",
           r"tropomyosin: .*", r"endocytic coat: .*", r"amphiphysin: .*", r"actin patches: \d.*", r"ATP synthase dimers: .*",
           r"mtDNA: .*", r"Cvt complex: .*", r"nucleus-vacuole junction: .*", r"polysomes: .*", r"RAC on polysome ribosomes: .*",
           r"mRNA: .*", r"mRNP ends: .*"]          # = plugins/yeast/repro/compare_x5.py s3_lines
SER = "preplaced/serial"                                # serial output folder of our runs ("prepass/serial" before A0c)


def _rp(base, rel):
    """<base>/<rel>, or its pre-A0c name in an older run (wholecell.build.ref_path)"""
    sys.path.insert(0, REPO)
    from wholecell.build import ref_path
    return ref_path(base, rel)


TIMING = re.compile(r"(?<![\w.])\d+(?:\.\d+)?s\b")    # addendum 2026-11-17 section 3


def rundir(tag):
    return f"{OUT}/" + os.environ.get(f"A6_RUN_{tag}", f"run-{tag}")


def layout():
    return dict(root=f"{OUT}/engine/root", work=f"{OUT}/engine/work", lib=f"{OUT}/engine/lib", pdb=f"{OUT}/engine/pdb",
                code=f"{REPO}/plugins/yeast/engine", inputs=f"{YEAST}/inputs", seed=0, nproc=3)


def prepare():
    E = layout(); W0 = f"{S0}/engine/work"
    for k in ("root", "work", "pdb"):
        os.makedirs(E[k], exist_ok=True)
    links = {"data": os.readlink(f"{S0}/_code/plugins/yeast/engine/pipeline/data"),
             "cell_meshes.npz": f"{PIPE}/cell_meshes.npz", "vacuole_template.npz": f"{PIPE}/vacuole_template.npz"}
    for n, t in links.items():
        if not os.path.lexists(f"{E['root']}/{n}"):
            os.symlink(t, f"{E['root']}/{n}")
    if not os.path.exists(E["lib"]):
        subprocess.run(["cp", "-a", f"{S0}/engine/lib", E["lib"]], check=True)
    for d in ("results",):
        if not os.path.exists(f"{E['work']}/{d}"):
            subprocess.run(["cp", "-a", f"{W0}/{d}", f"{E['work']}/{d}"], check=True)
    if not os.path.exists(f"{E['work']}/rcsb_cache.json"):
        subprocess.run(["cp", "-a", f"{W0}/rcsb_cache.json", E["work"]], check=True)
    y = f"{E['work']}/yc3d"; os.makedirs(y, exist_ok=True)
    for f in ("geometry.json", "geometry_fields.npz", "tile_measure.parquet", "tile_labvol.parquet", "pack_trees_15.npz"):
        if not os.path.exists(f"{y}/{f}"):
            subprocess.run(["cp", "-a", f"{W0}/yc3d/{f}", y], check=True)
    pdb0 = f"{S0}/engine/pdb"
    for f in sorted(os.listdir(pdb0)):
        if not os.path.lexists(f"{E['pdb']}/{f}"):
            os.symlink(os.path.realpath(f"{pdb0}/{f}"), f"{E['pdb']}/{f}")
    print("prepared", OUT, flush=True)


def config():
    import yaml
    cfg = yaml.safe_load(open(f"{REPO}/config/cells/yeast.yaml"))
    sc = cfg["preplaced"]["serial"]
    sc["context_params"]["engine"] = layout()
    return cfg


class Cfg(dict):
    pass


def run(tag):
    if os.environ.get("PYTHONHASHSEED") != "0":
        raise SystemExit("run with PYTHONHASHSEED=0")
    rdir = rundir(tag)
    os.makedirs(rdir, exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = f"{rdir}/numba_cache"
    os.environ.setdefault("NUMBA_NUM_THREADS", "3")
    sys.path.insert(0, REPO)
    import numpy as np
    import numba
    from wholecell import registry
    from wholecell.preplaced import builders  # noqa: F401
    from wholecell.preplaced import serial
    cfg = Cfg(config())
    registry.load_plugins(cfg)
    t0 = time.time()
    ctx = serial.run(cfg, rdir, out=SER)
    if ctx is None:
        raise SystemExit(f"{rdir}/{SER} exists: remove it to rerun")
    # pools after the pre-pass, through the stage-3 shell capacity step (s0 saved expected_totals after it)
    wc3d = ctx.mod("wc3d")
    exp3 = wc3d.shell_capacity(ctx.exp3, ctx.measure, ctx.lp, log=ctx.log) if wc3d.SHELL_CAP else ctx.exp3
    exp3.to_parquet(f"{rdir}/{SER}/expected_totals_after_shell_capacity.parquet")
    json.dump(dict(seconds=round(time.time() - t0, 1), numba=numba.__version__, numpy=np.__version__, python=sys.version,
                   host=os.uname().nodename), open(f"{rdir}/env.json", "w"), indent=1)


def _strip(line):
    return re.sub(r"^\d\d:\d\d:\d\d\s+", "", line.rstrip("\n"))


def _preplaced_lines(path):
    L = [_strip(l) for l in open(path)]
    i0 = next(i for i, l in enumerate(L) if l.startswith("microtubules: "))
    i1 = next(i for i, l in enumerate(L) if l.startswith("mRNP ends: "))
    return [TIMING.sub("<t>", l) for l in L[i0:i1 + 1] if not l.startswith(("[preplaced]", "[prepass]"))]


def compare(tag):
    import numpy as np
    import pandas as pd
    sys.path.insert(0, f"{REPO}/plugins/yeast/repro")
    import compare_x5 as cx
    rdir = rundir(tag)
    A = cx.npz_dict(f"{S0}/engine/work/yc3d/prepass.npz"); B = cx.npz_dict(_rp(rdir, f"{SER}/preplaced.npz"))
    rep = {}
    rep["S1"] = dict(level="E" if sorted(A) == sorted(B) else "fail", n_s0=len(A), n_new=len(B),
                     only_s0=sorted(set(A) - set(B)), only_new=sorted(set(B) - set(A)))
    cnt = {k: [(len(A[k]) if A[k].ndim else 1), (len(B[k]) if B[k].ndim else 1)] for k in sorted(set(A) & set(B))}
    bad = {k: v for k, v in cnt.items() if v[0] != v[1]}
    rep["S2"] = dict(level="E" if not bad and rep["S1"]["level"] == "E" else "fail", differing=bad)
    arr = {}
    for k in sorted(set(A) & set(B)):
        x, y = A[k], B[k]
        same = x.dtype == y.dtype and cx.arrays_equal(x, y)
        if not same:
            r = dict(dtype=[str(x.dtype), str(y.dtype)], shape=[list(x.shape), list(y.shape)])
            if x.shape == y.shape and x.dtype.kind in "fc":
                d = np.abs(x.astype(np.float64) - y.astype(np.float64))
                r.update(max_abs=float(np.nanmax(d)), n_diff=int((d > 0).sum()))
            arr[k] = r
    groups = {}
    for g, keys in cx.PRE_GROUPS.items():
        ks = [k for k in A if any(k == p or k.startswith(p + "_") for p in keys)]
        groups[g] = all(k in B and k not in arr for k in ks)
    rep["S4"] = dict(level="E" if all(groups.values()) and not arr else "fail", groups=groups, arrays_differing=arr,
                     all_arrays_identical=not arr and rep["S1"]["level"] == "E")
    la = [_strip(l) for l in open(f"{S0}/engine/work/s3_pack.log")]
    lb = [_strip(l) for l in open(f"{rdir}/run.log")]
    pat = {}
    for p in S3_PATS:
        ma = [TIMING.sub("<t>", l) for l in la if re.match(p, l)]; mb = [TIMING.sub("<t>", l) for l in lb if re.match(p, l)]
        if ma != mb or not ma:
            pat[p] = dict(s0=ma, new=mb)
    sa, sb = _preplaced_lines(f"{S0}/engine/work/s3_pack.log"), _preplaced_lines(f"{rdir}/run.log")
    seq = dict(identical=sa == sb, n_s0=len(sa), n_new=len(sb))
    if sa != sb:
        import difflib
        seq["diff"] = list(difflib.unified_diff(sa, sb, "s0", "new", lineterm="", n=0))[:60]
    rep["S3"] = dict(level="E" if not pat else "fail", patterns=len(S3_PATS), differing=pat, full_sequence=seq)
    ea = pd.read_parquet(f"{S0}/engine/work/yc3d/expected_totals.parquet")
    eb = pd.read_parquet(_rp(rdir, f"{SER}/expected_totals_after_shell_capacity.parquet"))
    same = ea.reset_index(drop=True).equals(eb.reset_index(drop=True)) and list(ea.columns) == list(eb.columns)
    e = dict(level="E" if same else "fail", rows=[len(ea), len(eb)])
    if not same and len(ea) == len(eb):
        d = (ea.expected.values - eb.expected.values)
        e.update(max_abs=float(np.abs(d).max()), n_diff=int((d != 0).sum()))
    rep["E_pools"] = e
    rep["layers"] = _layers_union(rdir, B)
    rep["env"] = json.load(open(f"{rdir}/env.json"))
    json.dump(rep, open(f"{rdir}/report.json", "w"), indent=1, default=str)
    md = [f"# A6 serial pre-pass vs port s0 ({tag})", "", "| row | level | note |", "|---|---|---|"]
    md.append(f"| S1 | {rep['S1']['level']} | {len(A)} / {len(B)} arrays |")
    md.append(f"| S2 | {rep['S2']['level']} | {len(bad)} counts differ |")
    md.append(f"| S3 | {rep['S3']['level']} | {len(S3_PATS) - len(pat)}/{len(S3_PATS)} patterns identical (timing tokens removed); full pre-pass log sequence identical: {seq['identical']} ({seq['n_s0']}/{seq['n_new']} lines) |")
    md.append(f"| S4 | {rep['S4']['level']} | {sum(groups.values())}/{len(groups)} groups identical; arrays differing: {sorted(arr)[:12]} |")
    md.append(f"| E pools | {e['level']} | expected totals after the pre-pass (+ shell capacity) vs s0 expected_totals.parquet |")
    md.append(f"| layers | {'E' if rep['layers']['union_identical'] else 'fail'} | per-layer files reassemble preplaced.npz |")
    md.append(f"\nenv: {rep['env']}")
    open(f"{rdir}/report.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))
    return rep


def _layers_union(rdir, B):
    import numpy as np
    man = json.load(open(_rp(rdir, f"{SER}/serial_manifest.json")))
    U = {}
    for r in man["builders"]:
        p = _rp(rdir, f"{SER}/layers/{r['name']}.npz")
        if os.path.exists(p):
            with np.load(p) as z:
                for k in z.files:
                    U[k] = z[k]
    ok = sorted(U) == sorted(B) and all(U[k].dtype == B[k].dtype and np.array_equal(U[k], B[k]) for k in B)
    return dict(union_identical=bool(ok), order=[(r["name"], r["builder"], len(r["keys"]), r["seconds"]) for r in man["builders"]])


def compare_ref(tag, ref, n_s0=18):
    """layers 1..n_s0 vs s0 prepass.npz (exact), later layers vs ref/preplaced/serial/layers (tolerance report)"""
    import numpy as np
    sys.path.insert(0, f"{REPO}/plugins/yeast/repro")
    import compare_x5 as cx
    rdir = rundir(tag)
    A0 = cx.npz_dict(f"{S0}/engine/work/yc3d/prepass.npz")
    man = json.load(open(_rp(rdir, f"{SER}/serial_manifest.json")))
    mref = {r["name"]: r for r in json.load(open(_rp(ref, f"{SER}/serial_manifest.json")))["builders"]}
    last = {}                                   # key -> last builder writing it (s0 holds only the final arrays)
    for i, r in enumerate(man["builders"], 1):
        for k in r["keys"]:
            last[k] = i
    rows = []
    for i, r in enumerate(man["builders"], 1):
        p = _rp(rdir, f"{SER}/layers/{r['name']}.npz")
        B = cx.npz_dict(p) if os.path.exists(p) else {}
        q = _rp(ref, f"{SER}/layers/{r['name']}.npz")
        Aref = cx.npz_dict(q) if os.path.exists(q) else {}
        if i <= n_s0:
            # a key a later builder replaces (nuc_* after chromatin/pol2) is not in s0 in this state: vs ref's layer
            A = {k: (A0[k] if last[k] == i and k in A0 else Aref.get(k)) for k in B}
            A = {k: v for k, v in A.items() if v is not None}
            n_via = sum(1 for k in B if last[k] != i)
            against = "s0" + (f" ({n_via} superseded arrays vs ref)" if n_via else "")
        else:
            A, against = Aref, "ref"
        diffs = {}
        for k in sorted(set(A) | set(B)):
            if k not in A or k not in B:
                diffs[k] = dict(missing="ref" if k not in A else "new"); continue
            x, y = A[k], B[k]
            if x.dtype == y.dtype and cx.arrays_equal(x, y):
                continue
            d = dict(dtype=[str(x.dtype), str(y.dtype)], shape=[list(x.shape), list(y.shape)])
            if x.shape == y.shape and x.dtype.kind in "fciub":
                dd = np.abs(x.astype(np.float64) - y.astype(np.float64))
                d.update(max_abs=float(np.nanmax(dd)) if dd.size else 0.0, n_diff=int((dd > 0).sum()), n=int(dd.size))
            diffs[k] = d
        rows.append(dict(i=i, name=r["name"], against=against, arrays=len(B), identical=not diffs, diffs=diffs,
                         sec=r["seconds"], sec_ref=mref.get(r["name"], {}).get("seconds")))
    # log lines from the first differing builder on (amphiphysin onwards): ref vs new
    la, lb = _preplaced_lines(f"{ref}/run.log"), _preplaced_lines(f"{rdir}/run.log")
    import difflib
    rep = dict(tag=tag, ref=ref, n_s0=n_s0, layers=rows, log_identical=la == lb,
               log_diff=list(difflib.unified_diff(la, lb, "ref", "new", lineterm="", n=0))[:80])
    json.dump(rep, open(f"{rdir}/compare_ref.json", "w"), indent=1, default=str)
    md = [f"# serial pre-pass {tag}: builders 1-{n_s0} vs s0, later vs {ref}", "",
          "| # | builder | against | arrays | identical | differing arrays (max abs diff / n differing of n) |", "|---|---|---|---|---|---|"]
    for r in rows:
        dd = "; ".join(f"{k}: " + (f"missing in {v['missing']}" if "missing" in v else
                       (f"{v.get('max_abs', float('nan')):.3g} / {v.get('n_diff')} of {v.get('n')}" if "max_abs" in v else f"shape {v['shape']}"))
                       for k, v in r["diffs"].items())
        md.append(f"| {r['i']} | {r['name']} | {r['against']} | {r['arrays']} | {r['identical']} | {dd} |")
    md.append(f"\npre-pass log lines (timing removed) identical to ref: {rep['log_identical']}")
    if rep["log_diff"]:
        md += ["```"] + rep["log_diff"] + ["```"]
    open(f"{rdir}/compare_ref.md", "w").write("\n".join(md) + "\n")
    print("\n".join(md))
    return rep


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prepare":
        prepare()
    elif cmd == "run":
        run(sys.argv[2])
    elif cmd == "compare":
        compare(sys.argv[2])
    elif cmd == "compare_ref":
        compare_ref(sys.argv[2], sys.argv[3], *(int(a) for a in sys.argv[4:5]))
    else:
        raise SystemExit("prepare | run <tag> | compare <tag> | compare_ref <tag> <ref_dir> [n_s0]")
