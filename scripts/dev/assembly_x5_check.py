"""Check of the core assembly preset legacy_x5 against the yeast port build s0 (unification A4, yeast half).

Runs the core assembly step list (wholecell/inventory/assembly, preset legacy_x5, yeast plugin steps and structure
hooks) on the stage-1 inputs of the s0 build and compares every assembly table with the one s0 wrote (row-keyed, as
plugins/yeast/repro/compare_x5.py after 8f1d0e1): ingredients.csv (I1) and oligo_map.json, f_actin.json,
complexes_tried.csv, complexes_added.csv, heterodimers.csv, oligomers.csv (A1). Level E = identical files content.

  1. prepare (yeast env; once): engine scratch layout in OUT/engine (YC_ROOT with the bundle data, YC_WORK with the
     frozen RCSB cache, YC_LIB = the s0 library as it was when s1 started, YC_PDB = link farm) and the stage-1 input
     table (pipeline/inventory.build() of the vendored engine, the protein table s1 starts from) -> OUT/inputs/
        conda_envs/yeast_x5/bin/python scripts/assembly_x5_check.py prepare
  2. run + compare (core env):
        PYTHONHASHSEED=0 conda_envs/vcell/bin/python scripts/assembly_x5_check.py run vcell
     -> OUT/run-<tag>/inventory/assembly/*.{csv,json}, OUT/report-<tag>.{json,md} (tag: a label for the env;
     run it again with conda_envs/yeast_x5 and tag yeast_x5 to check the yeast env too)
Nothing is written into the s0 build or the yeast tree.
"""
from __future__ import annotations

import json
import os
import shutil
import sys
import time

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _devsite import VC, RESULTS, YEAST as _YEAST_ROOT, PY as _SITE_PY, SMK as _SITE_SMK, GIT as _DEV_GIT  # noqa: E402,F401  (release: site paths)
S0 = f"{RESULTS}/yeast/20261116-aabe7b5-s0-x5repro"
OUT = f"{RESULTS}/yeast/_dev_unify/a4"
YEAST = _YEAST_ROOT
FROZEN = f"{VC}/data/derived/wholecell_frozen/yeast_x5"
S1_START = "2026-11-16 18:28:00"         # s1_inventory of s0 started 18:28:05 (its log); library entries older = its start set
TABLES_JSON = ("oligo_map.json", "f_actin.json")
TABLES_CSV = ("ingredients.csv", "complexes_tried.csv", "complexes_added.csv", "heterodimers.csv", "oligomers.csv")


def engine_layout():
    return dict(root=f"{OUT}/engine/root", work=f"{OUT}/engine/work", lib=f"{OUT}/engine/lib", pdb=f"{OUT}/engine/pdb",
                code=f"{REPO}/plugins/yeast/engine", inputs=f"{YEAST}/inputs")


def prepare():
    E = engine_layout()
    for k in ("root", "work", "lib", "pdb"):
        os.makedirs(E[k], exist_ok=True)
    if not os.path.lexists(f"{E['root']}/data"):
        os.symlink(f"{YEAST}/yeast_cell_dev_bundle/pipeline/data", f"{E['root']}/data")
    shutil.copy2(f"{FROZEN}/work_seed/rcsb_cache.json", f"{E['work']}/rcsb_cache.json")
    t0 = time.mktime(time.strptime(S1_START, "%Y-%m-%d %H:%M:%S"))
    n = 0
    for f in sorted(os.listdir(f"{S0}/engine/lib")):
        p = f"{S0}/engine/lib/{f}"
        if os.path.getmtime(p) < t0 and not os.path.exists(f"{E['lib']}/{f}"):
            shutil.copy2(p, f"{E['lib']}/{f}"); n += 1
    pdb_dir = f"{YEAST}/yeast_cell_dev_bundle/work/pdb"
    for f in sorted(os.listdir(pdb_dir)):
        if not os.path.lexists(f"{E['pdb']}/{f}"):
            os.symlink(f"{pdb_dir}/{f}", f"{E['pdb']}/{f}")
    print(f"library start set: {n} entries copied (mtime < {S1_START})", flush=True)
    sys.path.insert(0, REPO)
    from plugins.yeast import assembly as YA
    YA.engine_env(E)
    import inventory as inv
    df, _ = inv.build()
    os.makedirs(f"{OUT}/inputs", exist_ok=True)
    df.to_parquet(f"{OUT}/inputs/protein_table.parquet")
    json.dump({k: v for k, v in df.attrs.items()}, open(f"{OUT}/inputs/protein_table.attrs.json", "w"), indent=1, default=str)
    import pandas as pd
    pf = pd.read_parquet(f"{S0}/engine/work/features/protein_features.parquet").set_index("orf")
    d = df.set_index("orf").copies_per_cell.reindex(pf.index)
    print(f"stage-1 input: {len(df)} proteins; copies_per_cell vs s0 protein_features: max |d| "
          f"{float((d - pf.copies_per_cell).abs().max()):.3g}", flush=True)


def config(run_dir):
    """the yeast cell's assembly_steps params (config/cells/yeast.yaml) with this check's inputs"""
    import yaml
    cfg = yaml.safe_load(open(f"{REPO}/config/cells/yeast.yaml"))
    p = cfg["inventory"]["params"]["assembly_steps"]
    p = json.loads(json.dumps(p))
    p["engine"] = engine_layout()
    return p


def compare(ref_dir, new_dir):
    import numpy as np
    import pandas as pd
    sys.path.insert(0, f"{REPO}/plugins/yeast/repro")
    res = {}
    for f in TABLES_JSON:
        A = json.load(open(f"{ref_dir}/{f}")); B = json.load(open(f"{new_dir}/{f}"))
        if A == B:
            res[f] = dict(level="E"); continue
        ka, kb = set(A), set(B)
        diff = sorted(k for k in ka & kb if A[k] != B[k])
        res[f] = dict(level="fail", only_ref=sorted(ka - kb)[:20], only_new=sorted(kb - ka)[:20], differ=diff[:20],
                      n_only_ref=len(ka - kb), n_only_new=len(kb - ka), n_differ=len(diff),
                      examples={k: [A[k], B[k]] for k in diff[:5]})
    for f in TABLES_CSV:
        A = pd.read_csv(f"{ref_dir}/{f}"); B = pd.read_csv(f"{new_dir}/{f}")
        if A.equals(B):
            res[f] = dict(level="E", rows=len(A)); continue
        r = dict(level="fail", rows_ref=len(A), rows_new=len(B))
        if list(A.columns) != list(B.columns):
            r["columns"] = [list(A.columns), list(B.columns)]
            res[f] = r; continue
        k = A.columns[0]
        if A[k].is_unique and B[k].is_unique:
            ka, kb = set(A[k].astype(str)), set(B[k].astype(str))
            r.update(n_only_ref=len(ka - kb), n_only_new=len(kb - ka), only_ref=sorted(ka - kb)[:20], only_new=sorted(kb - ka)[:20])
            Ai = A.set_index(A[k].astype(str)); Bi = B.set_index(B[k].astype(str))
            common = [x for x in Ai.index if x in kb]
            Ai, Bi = Ai.loc[common], Bi.loc[common]
            same_order = list(A[k].astype(str)) == list(B[k].astype(str))
            r["row_order_same"] = same_order
            cols = {}
            for c in A.columns:
                x, y = Ai[c], Bi[c]
                eq = (x == y) | (x.isna() & y.isna())
                if not eq.all():
                    bad = list(eq[~eq].index)
                    e = dict(n=len(bad), keys=bad[:10], examples=[[str(x[b]), str(y[b])] for b in bad[:3]])
                    if x.dtype.kind in "fi" and y.dtype.kind in "fi":
                        e["max_rel"] = float(np.nanmax(np.abs(x - y) / np.maximum(np.abs(x), 1e-12)))
                    cols[c] = e
            r["columns_differing"] = cols
            if not cols and not (ka ^ kb):
                r["level"] = "E-keyed" if not same_order else "E-values"
        else:
            Ak = A.astype(str).agg("|".join, axis=1); Bk = B.astype(str).agg("|".join, axis=1)
            sa, sb = set(Ak), set(Bk)
            r.update(n_rows_only_ref=len(sa - sb), n_rows_only_new=len(sb - sa), rows_only_ref=sorted(sa - sb)[:10],
                     rows_only_new=sorted(sb - sa)[:10])
            if sa == sb and len(A) == len(B):
                r["level"] = "E-keyed"
        res[f] = r
    return res


def run(tag="vcell"):
    os.environ.setdefault("PYTHONHASHSEED", "0")
    sys.path.insert(0, REPO)
    t0 = time.time()
    from wholecell.inventory.assembly import run_standalone
    run_dir = f"{OUT}/run-{tag}"
    if os.path.exists(f"{run_dir}/inventory/assembly"):
        shutil.rmtree(f"{run_dir}/inventory/assembly")
    p = config(run_dir)
    run_standalone(run_dir, p, protein_table=f"{OUT}/inputs/protein_table.parquet", plugins=["plugins.yeast"])
    wall = time.time() - t0
    res = compare(f"{S0}/engine/work/results/inventory", f"{run_dir}/inventory/assembly")
    lev = "E" if all(r["level"] == "E" for r in res.values()) else "fail"
    rep = dict(reference=f"{S0}/engine/work/results/inventory", run=f"{run_dir}/inventory/assembly", preset=p["preset"],
               level=lev, wall_s=round(wall, 1), tables=res, commit=os.popen(f"git -C {REPO} rev-parse --short HEAD").read().strip())
    json.dump(rep, open(f"{OUT}/report-{tag}.json", "w"), indent=1, default=str)
    lines = [f"# assembly legacy_x5 vs s0 ({time.strftime('%Y-%m-%d %H:%M')}; {tag}: Python {sys.version.split()[0]})", "",
             f"overall: **{lev}**; run {wall:.0f} s; worktree HEAD {rep['commit']}", "", "| table | level | notes |", "|---|---|---|"]
    for f, r in res.items():
        note = {k: v for k, v in r.items() if k != "level"}
        lines.append(f"| {f} | {r['level']} | {json.dumps(note, default=str)[:600]} |")
    open(f"{OUT}/report-{tag}.md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0 if lev == "E" else 1


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "run"
    if mode == "prepare":
        prepare()
    elif mode == "compare":
        r = compare(f"{S0}/engine/work/results/inventory", sys.argv[2])
        print(json.dumps(r, indent=1, default=str))
    else:
        sys.exit(run(sys.argv[2] if len(sys.argv) > 2 else "vcell"))
