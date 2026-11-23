#!/usr/bin/env python
"""Step-2 subset test harness (frozen definition: docs/unification_criteria/step2_subsets.md).

    <core python> scripts/step2_subset_test.py run --name N --config config/cells/<x>.yaml --commit SHA
                  [--reuse-library BUILD] [--reuse-tiles BUILD] [--reuse-upstream BUILD] [--workers 24]
                  [--ref RUN_DIR | --ref-build BUILD [--ref-port]] [--expect-bitwise] [--expect-same STAGES]
    <core python> scripts/step2_subset_test.py metrics --run RUN_DIR            (re-compute result.json)
    <core python> scripts/step2_subset_test.py compare --run RUN_DIR (--ref RUN_DIR | --ref-build BUILD [--ref-port])
    <core python> scripts/step2_subset_test.py select [--check]                  (re-derive the frozen tile lists)

`run` makes {results_root}/<cell_id>/_dev_unify/step2/<name>/ (refuses an existing one): a git-archive snapshot of
--commit in _code/, cell.yaml (extends the snapshot's --config, allow_dirty), the snapshot DAG's jobs up to the pack
counts and the pre-placed index executed one by one (Snakemake dry-run order), reused stages linked (declared and
checked), the subset pack with the reference-context neighbour rule (§3, `_pack`, run with the snapshot's code), the
subset metrics (§4, harness code) -> result.json, and with a reference the comparison (§5) -> compare.json.
Nothing is written outside the run folder (reference builds are read through symlink views inside it).
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from _devsite import VC, RESULTS, YEAST as _YEAST_ROOT, PY as _SITE_PY, SMK as _SITE_SMK, GIT as _DEV_GIT  # noqa: E402,F401  (release: site paths)
R = f"{RESULTS}"
PY = _SITE_PY
SMK = _SITE_SMK
GIT = _DEV_GIT                                 # main checkout (worktrees share its object store)
SUBSETS = {"jurkat": "config/subsets/jurkat_step2.tsv", "yeast": "config/subsets/yeast_step2.tsv"}
SUBSET_MD5 = {"jurkat": "dd7c4bf95419d7a4d5acc6358da8c359", "yeast": "cf4071d2756f76b84dbccb4f35ae083f"}   # ref_tile int64
QUICK_MD5 = {"jurkat": "0b4e174284050610db739649de115a19", "yeast": "51a1927c4cac25b275659e281277ceb6"}
FRAME = {"jurkat": f"{R}/jurkat/20261117-bc121a6-s20261115-b3", "yeast": f"{R}/yeast/20261118-9977542-s0-a9-core"}
CONTEXT = {"jurkat": dict(kind="core", build=f"{R}/jurkat/20261117-bc121a6-s20261115-b3"),
           "yeast": dict(kind="port", build=f"{R}/yeast/20261116-aabe7b5-s0-x5repro")}
SELECT_SEED = 20261118
SELECT_RULE = {"jurkat": dict(p75=True, busiest=6, membrane_rich=4, seam={"cytosol": 3, "nucleoplasm": 3}, components={}),
               "yeast": dict(p75=False, busiest=4, membrane_rich=3, seam={"cytosol": 3},
                             components={"peroxisome": ["peroxisome", "mem_PEROX", "PEROX_peripheral"],
                                         "lipid_droplet": ["LD_peripheral"]})}
SEAM_REGIONS = {"jurkat": ["cytosol", "nucleoplasm"], "yeast": ["cytosol"]}
RIBO_LIBS = ("RIBO80S", "ribosome_80S")
YEAST_EXTRA_LIBS = ("FAA4",)
AV = 6.022e23


def jload(p, default=None):
    try:
        return json.load(open(p))
    except Exception:  # noqa: BLE001
        return default


def jsave(p, obj):
    tmp = p + ".tmp"
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    os.replace(tmp, p)


def md5(p):
    if not p or not os.path.exists(p):
        return None
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def log(run, *a):
    msg = time.strftime("%F %T ") + " ".join(str(x) for x in a)
    print(msg, flush=True)
    with open(os.path.join(run, "status.txt"), "a") as f:
        f.write(msg + "\n")


# ------------------------------------------------------------------ subsets (§1)
def select(cell):
    """re-derive the frozen tile list of a cell from its frame build (rule §1)"""
    import pandas as pd
    from scipy import ndimage
    b = FRAME[cell]; ru = SELECT_RULE[cell]
    M = pd.read_parquet(f"{b}/pack/tiles.parquet")
    T = float(json.load(open(f"{b}/pack/tiles_meta.json"))["TILE"])
    cz = np.load(f"{b}/pack/counts.npz")
    req = np.bincount(cz["tile"], weights=cz["n"], minlength=len(M))
    cols = [c for c in M.columns if c not in ("tile", "ti", "tj", "tk") and (M[c] > 0).any()]
    mcols = [c for c in cols if c.startswith("mem_")]
    chosen, why = [], {}

    def add(t, r):
        t = int(t)
        if t not in why:
            chosen.append(t); why[t] = []
        why[t].append(r)

    def take(order, r, n=1):
        k = 0
        for t in order:
            if int(t) not in why:
                add(t, r); k += 1
                if k == n:
                    return
    for name, cs in ru["components"].items():
        s = M[(M[cs] > 0).any(axis=1)]
        ijk = s[["ti", "tj", "tk"]].values; lo = ijk.min(0)
        g = np.zeros(ijk.max(0) - lo + 1, np.int32); g[tuple((ijk - lo).T)] = 1
        lab, n = ndimage.label(g, structure=np.ones((3, 3, 3)))
        li = lab[tuple((ijk - lo).T)]
        sizes = np.bincount(li)[1:]
        k = 1 + int(np.lexsort((np.arange(n), sizes))[0])
        for t in sorted(s.tile.values[li == k]):
            add(t, f"whole_{name}")
    for c in cols:
        v = M[c].values
        take(np.lexsort((M.tile.values, -v)), f"top:{c}")
        if ru["p75"]:
            nz = np.nonzero(v > 0)[0]
            q = np.quantile(v[nz], 0.75)
            take(nz[np.lexsort((M.tile.values[nz], np.abs(v[nz] - q)))], f"p75:{c}")
    take(np.lexsort((M.tile.values, -req)), "busiest", ru["busiest"])
    take(np.lexsort((M.tile.values, -M[mcols].sum(axis=1).values)), "membrane_rich", ru["membrane_rich"])
    rng = np.random.default_rng(SELECT_SEED)
    for reg_, n in ru["seam"].items():
        q = M.tile.values[(M[reg_].values / T ** 3) >= 0.95]
        take(rng.permutation(np.sort(q)), f"seam:{reg_}", n)
    return np.array(sorted(chosen), np.int64)


QUICK = {"jurkat": "config/subsets/jurkat_quick.tsv", "yeast": "config/subsets/yeast_quick.tsv"}
QUICK_MAX = 18


def select_quick(cell):
    """quick development subset (docs/unification_criteria/step2_quick_subsets_2026-11-18.md): representative, not
    busy tiles. Candidates: tiles with requested instances <= the 75th percentile of the tiles with requests. Region
    measure columns in order of their cell total (largest first): the candidate whose measure is closest to the median
    of the candidates' non-zero measures; then one cytosol tile >= 95 % (median requested). Up to QUICK_MAX tiles, then
    tiles for any region of libs.json rnames still without requests (its median tile; ties: lowest id)."""
    import pandas as pd
    b = FRAME[cell]
    M = pd.read_parquet(f"{b}/pack/tiles.parquet")
    T = float(json.load(open(f"{b}/pack/tiles_meta.json"))["TILE"])
    cz = np.load(f"{b}/pack/counts.npz")
    req = np.bincount(cz["tile"], weights=cz["n"], minlength=len(M))
    q75 = np.percentile(req[req > 0], 75)
    cand = (req > 0) & (req <= q75)
    cols = [c for c in M.columns if c not in ("tile", "ti", "tj", "tk") and (M[c].values[cand] > 0).any()]
    cols = sorted(cols, key=lambda c: (-float(M[c].sum()), c))
    chosen, why = [], {}

    def median_tile(v, mask):
        nz = np.nonzero(mask & (v > 0))[0]
        if not len(nz):
            return None
        med = np.median(v[nz])
        for t in nz[np.lexsort((nz, np.abs(v[nz] - med)))]:
            if int(t) not in why:
                return int(t)
        return None
    for c in cols:
        if len(chosen) >= QUICK_MAX - 1:
            break
        t = median_tile(M[c].values, cand)
        if t is not None:
            chosen.append(t); why[t] = f"median:{c}"
    cyt = cand & (M["cytosol"].values / T ** 3 >= 0.95)
    t = median_tile(req, cyt)
    if t is not None:
        chosen.append(t); why[t] = "seam:cytosol"
    rn = json.load(open(f"{b}/pack/libs.json"))["rnames"]
    for _ in range(4):                                     # coverage of every requested region
        have = set(cz["region"][np.isin(cz["tile"], chosen)].tolist())
        miss = [r for r in sorted(set(cz["region"].tolist())) if r not in have]
        if not miss:
            break
        r = miss[0]
        rq = np.bincount(cz["tile"][cz["region"] == r], weights=cz["n"][cz["region"] == r], minlength=len(M))
        t = median_tile(rq, rq > 0)
        if t is None:
            break
        chosen.append(t); why[t] = f"cover:{rn[r]}"
    ids = np.array(sorted(chosen), np.int64)
    S = M.set_index("tile").loc[ids]
    return pd.DataFrame(dict(cell=cell, ref_tile=ids, ti=S.ti.values, tj=S.tj.values, tk=S.tk.values,
                             requested_ref=req[ids].astype(np.int64), reasons=[why[t] for t in ids]))


def load_subset(cell, kind="frozen"):
    import pandas as pd
    f, m = (SUBSETS[cell], SUBSET_MD5[cell]) if kind == "frozen" else (QUICK[cell], QUICK_MD5[cell])
    S = pd.read_csv(os.path.join(REPO, f), sep="\t")
    h = hashlib.md5(S.ref_tile.values.astype(np.int64).tobytes()).hexdigest()
    if h != m:
        raise SystemExit(f"{f}: ref_tile md5 {h} != {m}")
    return S


def frame_lo0(cell):
    return np.asarray(json.load(open(f"{FRAME[cell]}/pack/tiles_meta.json"))["lo0"], float)


def subset_ids(cell, bdir, kind="frozen"):
    """subset tile ids in the frame of a core-layout build / run (by grid index; refuses tiles missing there)"""
    import pandas as pd
    S = load_subset(cell, kind)
    meta = json.load(open(os.path.join(bdir, "pack", "tiles_meta.json")))
    T = float(meta["TILE"]); off = np.round((frame_lo0(cell) - np.asarray(meta["lo0"], float)) / T).astype(int)
    M = pd.read_parquet(os.path.join(bdir, "pack", "tiles.parquet"), columns=["tile", "ti", "tj", "tk"])
    key = {tuple(t): int(i) for i, t in zip(M.tile.values, M[["ti", "tj", "tk"]].values)}
    ids = [key.get(tuple(t + off)) for t in S[["ti", "tj", "tk"]].values]
    if any(i is None for i in ids):
        raise SystemExit(f"{bdir}: {sum(i is None for i in ids)} subset tiles are not tiles of this build")
    ids = np.array(ids, np.int64)
    return ids, bool(np.array_equal(ids, S.ref_tile.values))


# ------------------------------------------------------------------ run: snapshot, DAG jobs, reuse (§2)
def cell_id_of(cfg_path):
    sys.path.insert(0, REPO)
    from wholecell.config import load
    return load(cfg_path)["cell_id"]


def dry_run_jobs(run, cfgp, workers, targets):
    cmd = [SMK, "--snakefile", os.path.join(run, "_code", "Snakefile"), "--directory", run, "--cores", str(workers),
           "--config", f"cell={cfgp}", f"build_id={os.path.basename(run)}", f"build_dir={run}", "-n", "-p", *targets]
    r = subprocess.run(cmd, cwd=run, capture_output=True, text=True, env=dict(os.environ, PYTHONHASHSEED="0"), timeout=900)
    if r.returncode != 0:
        raise RuntimeError(f"snakemake dry run failed:\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}")
    jobs, cur = [], None
    for line in r.stdout.splitlines():
        m = re.match(r"^(?:local)?rule (\S+):$", line.strip())
        if m:
            cur = dict(rule=m.group(1), cmd=None); jobs.append(cur); continue
        if cur is not None and line.strip().startswith("output:"):
            cur["outputs"] = [x.strip() for x in line.strip()[len("output:"):].split(",") if x.strip()]
        if cur is not None and cur["cmd"] is None and line.strip().startswith("cd ") and " -m wholecell " in line:
            cur["cmd"] = line.strip()
    for j in jobs:
        if j["rule"] == "geometry" and j["cmd"] is None:        # run: block in geometry.smk
            j["cmd"] = (f"cd {run}/_code && PYTHONHASHSEED=0 PYTHONPATH={run}/_code {PY} -m wholecell geometry --cell {cfgp} "
                        f"--build {run} > {run}/logs/geometry.log 2>&1")
    bad = [j["rule"] for j in jobs if not j["cmd"]]
    if bad:
        raise RuntimeError(f"no shell command parsed for {bad}")
    return jobs


def _norm(x, prefixes):
    if isinstance(x, dict):
        return {k: _norm(v, prefixes) for k, v in x.items()}
    if isinstance(x, list):
        return [_norm(v, prefixes) for v in x]
    if isinstance(x, str):
        for p in prefixes:
            x = x.replace(p, "{code}")
    return x


def resolved_config(run, cfgp):
    r = subprocess.run([PY, "-c", "import json,sys; from wholecell.config import load; print(json.dumps(load(sys.argv[1]), default=str))", cfgp],
                       cwd=os.path.join(run, "_code"), capture_output=True, text=True,
                       env=dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=os.path.join(run, "_code")))
    if r.returncode != 0:
        raise RuntimeError(r.stderr[-2000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


def check_same_section(run, cfgp, donor, section):
    a = resolved_config(run, cfgp).get(section)
    b = (jload(os.path.join(donor, "config.resolved.json")) or {}).get(section)
    pre = [os.path.join(run, "_code"), os.path.join(donor, "_code")]
    return _norm(a, pre) == _norm(b, pre), a, b


UPSTREAM_CODE = ("wholecell", "plugins", "toolkits", "workflow", "Snakefile", "config/geometry", "config/site.yaml")


def effective_switches(run, pack_cfg):
    """pack switches after presets (snapshot code), for the upstream-equivalence check"""
    r = subprocess.run([PY, "-c", "import json,sys; from wholecell.pack import switches; "
                        "e = switches.effective(json.loads(sys.argv[1])); "
                        "x = {k: v for k, v in e.items() if k not in ('_preset', 'preset', 'presets') and k not in switches.DEFAULTS}; "
                        "print(json.dumps(dict(effective=x, switches=switches.load(e).as_dict()), default=str, sort_keys=True))",
                        json.dumps(pack_cfg or {}, default=str)], cwd=os.path.join(run, "_code"), capture_output=True, text=True,
                       env=dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=os.path.join(run, "_code")))
    return r.stdout.strip().splitlines()[-1] if r.returncode == 0 and r.stdout.strip() else f"ERROR {r.stderr[-300:]}"


UPSTREAM_TOP = ("geometry", "inventory", "library", "engine", "preplaced", "frozen_inputs.tsv", "BUILD_STAMP.json", "_smk")
PACK_OUT = {"tiles", "rounds", "pilot", "pack_log.tsv", "pack_done.json", "summary.json"}


NOT_INVENTORY = [":(exclude)wholecell/inventory", ":(exclude)plugins/*/inventory/*", ":(exclude)plugins/*/data/localisation/*"]
PART_CODE = {"geometry": ["wholecell/geometry", "toolkits", "config/geometry"],
             "library": list(UPSTREAM_CODE_BASE := ("wholecell", "plugins", "toolkits", "workflow", "Snakefile", "config/geometry")) + NOT_INVENTORY,
             "tiles": ["wholecell/geometry", "wholecell/pack", "toolkits", "config/geometry"]}
PART_SECTIONS = {"geometry": ["geometry"], "tiles": ["geometry"], "inventory": ["inventory", "geometry"],
                 "library": ["library"], "engine": ["inventory", "library", "geometry"],
                 "preplaced": ["preplaced", "inventory", "library", "geometry"], "_smk": ["inventory", "library", "preplaced", "geometry"]}


def donor_commit_of(d):
    c = (jload(os.path.join(d, "BUILD_STAMP.json")) or {}).get("commit") or ""
    if not c and os.path.exists(os.path.join(d, "_code", "SNAPSHOT_SHA")):
        c = open(os.path.join(d, "_code", "SNAPSHOT_SHA")).read().strip()
    return c


AUTO_PARTS = ["geometry", "tiles", "library", "inventory", "engine", "_smk", "preplaced", "organisation"]
PART_DEPS = {"engine": ["inventory", "library"], "_smk": ["inventory", "library"], "preplaced": ["inventory", "library", "geometry"],
             "organisation": ["inventory", "library"], "inventory": ["geometry"], "tiles": ["geometry"]}


def reuse_parts(run, cfgp, spec, info, auto=False):
    """--reuse-parts DONOR:part,part: link (copy for _smk) the named upstream parts of DONOR after checking that the code
    they depend on (git diff) and their config sections are identical; the DAG then skips jobs whose outputs exist.
    --reuse-auto DONOR (auto=True): every part of AUTO_PARTS that passes the check and whose prerequisites are reused
    (PART_DEPS) is linked; the others are rebuilt (cached upstream: only changed stages run)"""
    donor, parts = (spec, ",".join(AUTO_PARTS)) if auto else spec.rsplit(":", 1)
    parts = [x for x in parts.split(",") if x]
    full = open(os.path.join(run, "_code", "SNAPSHOT_SHA")).read().strip()
    dc = donor_commit_of(donor)
    ra, rb = resolved_config(run, cfgp), jload(os.path.join(donor, "config.resolved.json")) or {}
    pre = [os.path.join(run, "_code"), os.path.join(donor, "_code"), os.path.realpath(os.path.join(donor, "_code"))]
    chk, ok_parts = {}, []
    for part in parts:
        if part == "organisation" and not os.path.isdir(os.path.join(donor, "pack", "organisation")):
            continue
        if part not in ("tiles", "organisation", "_smk") and not os.path.exists(os.path.join(donor, part)):
            continue
        paths = PART_CODE.get(part, list(UPSTREAM_CODE))
        code_same = dc == full or subprocess.run(["git", "-C", GIT, "diff", "--quiet", full, dc, "--", *paths]).returncode == 0
        secs = [k for k in PART_SECTIONS.get(part, []) if _norm(ra.get(k), pre) != _norm(rb.get(k), pre)]
        keys = {"tiles": ("tiling", "tile_nm", "meas_h_nm", "meas_vox_nm", "regions"), "organisation": ("organisation",)}.get(part)
        if keys:
            ta = {k: (ra.get("pack") or {}).get(k) for k in keys}; tb = {k: (rb.get("pack") or {}).get(k) for k in keys}
            if _norm(ta, pre) != _norm(tb, pre):
                secs.append("pack." + "/".join(keys))
        deps_missing = [d for d in PART_DEPS.get(part, []) if d in AUTO_PARTS and d not in ok_parts] if auto else []
        chk[part] = dict(code_identical=code_same, sections_differ=secs, prerequisites_rebuilt=deps_missing)
        if code_same and not secs and not deps_missing:
            ok_parts.append(part)
        elif not auto:
            raise SystemExit(f"--reuse-parts {part} from {donor}: not equivalent {chk[part]}")
    for part in ok_parts:
        if part in ("tiles", "organisation"):
            os.makedirs(os.path.join(run, "pack"), exist_ok=True)
            fs = ("tiles.parquet", "tiles_meta.json", "tile_labvol.parquet", "tile_geom") if part == "tiles" else ("organisation",)
            for f in fs:
                s_ = os.path.join(donor, "pack", f)
                if os.path.exists(s_) and not os.path.lexists(os.path.join(run, "pack", f)):
                    os.symlink(os.path.realpath(s_), os.path.join(run, "pack", f))
        elif part == "_smk":
            if os.path.isdir(os.path.join(donor, "_smk")) and not os.path.exists(os.path.join(run, "_smk")):
                shutil.copytree(os.path.join(donor, "_smk"), os.path.join(run, "_smk"))
        else:
            s_ = os.path.join(donor, part)
            if os.path.exists(s_) and not os.path.lexists(os.path.join(run, part)):
                os.symlink(os.path.realpath(s_), os.path.join(run, part))
    info.setdefault("reuse_parts", {})[donor] = dict(donor_commit=dc, auto=auto, reused=ok_parts, checks=chk)
    log(run, f"reused parts {ok_parts} of {donor} (commit {dc[:7]})" + (f"; rebuilt {[p for p in chk if p not in ok_parts]}" if auto else ""))


def late_reuse_preplaced(run, donor, info):
    """--reuse-auto: the pre-placed stage of DONOR is reused when the rebuilt inputs it reads are content-identical
    (pack table, library manifest, geometry package) and the pre-placed code / config are unchanged"""
    ins = ["inventory/pack_table.tsv", "library/library_manifest.tsv"]
    same_in = all(content_md5(os.path.join(run, f)) == content_md5(os.path.join(donor, f)) for f in ins)
    g = lambda b: (jload(os.path.join(b, "geometry", "package.json")) or {}).get("manifest_sha256")
    full = open(os.path.join(run, "_code", "SNAPSHOT_SHA")).read().strip(); dc = donor_commit_of(donor)
    code_same = dc == full or subprocess.run(["git", "-C", GIT, "diff", "--quiet", full, dc, "--", "wholecell/preplaced", "plugins",
                                              *NOT_INVENTORY]).returncode == 0
    ra, rb = resolved_config(run, os.path.join(run, "cell.yaml")), jload(os.path.join(donor, "config.resolved.json")) or {}
    pre = [os.path.join(run, "_code"), os.path.join(donor, "_code"), os.path.realpath(os.path.join(donor, "_code"))]
    cfg_same = _norm(ra.get("preplaced"), pre) == _norm(rb.get("preplaced"), pre)
    ok = same_in and g(run) == g(donor) and code_same and cfg_same and os.path.isdir(os.path.join(donor, "preplaced"))
    info.setdefault("late_reuse", {})["preplaced"] = dict(donor=donor, inputs_identical=same_in, code_identical=code_same,
                                                          config_identical=cfg_same, reused=ok)
    if ok:
        os.symlink(os.path.realpath(os.path.join(donor, "preplaced")), os.path.join(run, "preplaced"))
        log(run, f"pre-placed stage reused from {donor} (inputs content-identical)")


def link_upstream(run, donor):
    for e in UPSTREAM_TOP:
        s = os.path.join(donor, e)
        if os.path.exists(s) and not os.path.lexists(os.path.join(run, e)):
            os.symlink(s, os.path.join(run, e))
    os.makedirs(os.path.join(run, "pack"), exist_ok=True)
    for e in sorted(os.listdir(os.path.join(donor, "pack"))):
        if e in PACK_OUT or ".tmp" in e:
            continue
        d = os.path.join(run, "pack", e)
        if not os.path.lexists(d):
            os.symlink(os.path.join(donor, "pack", e), d)


def cmd_run(a):
    t_all = time.time()
    cfg_rel = a.config
    cid = cell_id_of(os.path.join(REPO, cfg_rel))
    cell = "yeast" if cid.startswith("yeast") else "jurkat"
    run = os.path.join(R, cid, "_dev_unify", "step2", a.name)
    if os.path.exists(run) and not a.resume:
        raise SystemExit(f"{run} exists (no re-use of run folders; --resume continues an interrupted one)")
    os.makedirs(os.path.join(run, "logs"), exist_ok=True)
    info = jload(os.path.join(run, "run.json"), {}) if a.resume else {}
    info.update(name=a.name, cell=cell, cell_id=cid, config=cfg_rel, commit=a.commit, workers=a.workers, host=os.uname().nodename,
                harness_commit=(open(os.path.join(REPO, "SNAPSHOT_SHA")).read().strip() if os.path.exists(os.path.join(REPO, "SNAPSHOT_SHA"))
                                else subprocess.run(["git", "-C", REPO, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()),
                subset=SUBSETS[cell] if a.subset == "frozen" else QUICK[cell], subset_kind=a.subset,
                subset_md5=SUBSET_MD5[cell] if a.subset == "frozen" else QUICK_MD5[cell], context=CONTEXT[cell], started=info.get("started") or time.strftime("%FT%T"),
                reuse=dict(library=a.reuse_library, tiles=a.reuse_tiles, upstream=a.reuse_upstream, parts=a.reuse_parts),
                overlay=(open(a.overlay).read() if a.overlay else None))
    info.setdefault("stages", {})
    jsave(os.path.join(run, "run.json"), info)
    # 1. snapshot + config
    code = os.path.join(run, "_code")
    if not os.path.exists(os.path.join(code, "SNAPSHOT_SHA")):
        os.makedirs(code, exist_ok=True)
        p1 = subprocess.Popen(["git", "-C", GIT, "archive", a.commit], stdout=subprocess.PIPE)
        subprocess.run(["tar", "-x", "-C", code], stdin=p1.stdout, check=True); p1.wait()
        full = subprocess.run(["git", "-C", GIT, "rev-parse", a.commit], capture_output=True, text=True).stdout.strip()
        open(os.path.join(code, "SNAPSHOT_SHA"), "w").write(full + "\n")
        info["commit_full"] = full
    cfgp = os.path.join(run, "cell.yaml")
    if not os.path.exists(cfgp):
        with open(cfgp, "w") as f:
            f.write(f"# step-2 subset test {a.name}: {cfg_rel} at {a.commit} (snapshot {code})\n"
                    f"extends: {os.path.join(code, cfg_rel)}\nbuild:\n  allow_dirty: true\n")
            if a.overlay:                                 # extra top-level keys (e.g. the B5 localisation preset)
                f.write(f"# overlay {a.overlay}\n" + open(a.overlay).read())
    log(run, f"run {a.name}: {cfg_rel} @ {a.commit}, cell {cell}")
    env = dict(os.environ, PYTHONHASHSEED="0", NUMBA_CACHE_DIR=os.path.join(run, "numba_cache"))
    # 2. upstream
    if a.reuse_upstream:
        same, _, _ = check_same_section(run, cfgp, a.reuse_upstream, "pack")
        donor_commit = (jload(os.path.join(a.reuse_upstream, "BUILD_STAMP.json")) or {}).get("commit") or ""
        if not donor_commit and os.path.exists(os.path.join(a.reuse_upstream, "_code", "SNAPSHOT_SHA")):   # snapshot builds
            donor_commit = open(os.path.join(a.reuse_upstream, "_code", "SNAPSHOT_SHA")).read().strip()
        full = open(os.path.join(code, "SNAPSHOT_SHA")).read().strip()
        if full != donor_commit:
            # another commit: allowed only if the code and every upstream config section are identical (the donor may
            # differ by pack.preset / pack.presets with the same effective packer switches, e.g. a b1_s1 overlay)
            dc = subprocess.run(["git", "-C", GIT, "diff", "--quiet", full, donor_commit, "--", *UPSTREAM_CODE]).returncode
            ra, rb = resolved_config(run, cfgp), jload(os.path.join(a.reuse_upstream, "config.resolved.json")) or {}
            pre = [os.path.join(run, "_code"), os.path.join(a.reuse_upstream, "_code")]
            secs = sorted(k for k in (set(ra) | set(rb)) - {"pack", "build", "build_id", "cell_config", "config_source"}
                          if not k.startswith("_"))           # _raw / _source / _site: the overlay file itself, site paths
            diff = [k for k in secs if _norm(ra.get(k), pre) != _norm(rb.get(k), pre)]
            pa = {k: v for k, v in (ra.get("pack") or {}).items() if k not in ("preset", "presets")}
            pb = {k: v for k, v in (rb.get("pack") or {}).items() if k not in ("preset", "presets")}
            sw_same = effective_switches(run, ra.get("pack")) == effective_switches(run, rb.get("pack"))
            ok = dc == 0 and not diff and _norm(pa, pre) == _norm(pb, pre) and sw_same
            info["upstream_equivalence"] = dict(donor_commit=donor_commit, code_identical=dc == 0, config_sections_differ=diff,
                                                pack_without_presets_identical=_norm(pa, pre) == _norm(pb, pre),
                                                effective_switches_identical=sw_same, accepted=ok)
            if not ok:
                raise SystemExit(f"--reuse-upstream {a.reuse_upstream}: commit {donor_commit} != {full} and not equivalent: "
                                 f"{info['upstream_equivalence']}")
            same = True
        link_upstream(run, a.reuse_upstream)
        info["stages"]["upstream"] = dict(status="reused", donor=a.reuse_upstream, donor_commit=donor_commit, pack_config_identical=same)
        log(run, f"upstream reused from {a.reuse_upstream} (pack config identical: {same})")
    else:
        for spec in a.reuse_parts or []:
            reuse_parts(run, cfgp, spec, info)
        if a.reuse_auto:
            reuse_parts(run, cfgp, a.reuse_auto, info, auto=True)
        targets = [os.path.join(run, "pack", "counts.npz"), os.path.join(run, "pack", "preplaced_off.npy")]
        jobs = dry_run_jobs(run, cfgp, a.workers, targets)
        info["dag_jobs"] = [j["rule"] for j in jobs]
        skip = set()
        if a.reuse_library:
            same, x, y = check_same_section(run, cfgp, a.reuse_library, "library")
            if not same:
                raise SystemExit(f"--reuse-library: library sections differ ({x} vs {y})")
            if not os.path.lexists(os.path.join(run, "library")):
                os.symlink(os.path.join(a.reuse_library, "library"), os.path.join(run, "library"))
            lib_dir = os.path.join(run, "library") + os.sep        # only jobs that write nothing but library/ (yeast engine_lib runs)
            skip |= {j["rule"] for j in jobs if j["rule"] == "library" or (j["rule"].startswith("library_") and j.get("outputs") and
                     all(o.startswith(lib_dir) or "/_smk/" in o for o in j["outputs"]))}
            info["stages"]["library"] = dict(status="reused", donor=a.reuse_library, manifest_md5=md5(os.path.join(run, "library", "library_manifest.tsv")))
        for j in jobs:
            name = j["rule"]
            st = info["stages"].get(name, {})
            if st.get("status") in ("ok", "reused"):
                continue
            if name in skip:
                info["stages"][name] = dict(status="reused", donor=a.reuse_library); continue
            if name == "pack_tiles" and a.reuse_tiles:
                pj = lambda b: (jload(os.path.join(b, "geometry", "package.json")) or {}).get("manifest_sha256")
                if pj(run) is None or pj(run) != pj(a.reuse_tiles):
                    raise SystemExit(f"--reuse-tiles: geometry differs ({pj(run)} vs {pj(a.reuse_tiles)})")
                for f in ("tiles.parquet", "tiles_meta.json", "tile_labvol.parquet", "tile_geom"):
                    s = os.path.join(a.reuse_tiles, "pack", f)
                    if os.path.exists(s):
                        os.makedirs(os.path.join(run, "pack"), exist_ok=True)
                        os.symlink(s, os.path.join(run, "pack", f))
                info["stages"][name] = dict(status="reused", donor=a.reuse_tiles); continue
            if name.startswith("preplaced_") and a.reuse_auto and not os.path.lexists(os.path.join(run, "preplaced")):
                late_reuse_preplaced(run, a.reuse_auto, info)
            outs = j.get("outputs") or []
            if outs and all(os.path.exists(o) for o in outs):
                info["stages"][name] = dict(status="reused", note="outputs present (reused parts)"); continue
            outside = [o for o in outs if not os.path.realpath(o).startswith(os.path.realpath(run) + os.sep)]
            if outside:
                raise SystemExit(f"stage {name} would write outside the run folder (reused part): {outside[:3]}")
            t0 = time.time()
            rc = subprocess.run(["bash", "-c", j["cmd"]], env=env).returncode
            info["stages"][name] = dict(status="ok" if rc == 0 else "failed", rc=rc, seconds=round(time.time() - t0, 1))
            jsave(os.path.join(run, "run.json"), info)
            log(run, f"  {name}: rc {rc}, {time.time() - t0:.0f} s")
            if rc != 0:
                raise SystemExit(f"stage {name} failed (see {run}/logs)")
    jsave(os.path.join(run, "run.json"), info)
    # 3. subset pack (snapshot code)
    if info["stages"].get("subset_pack", {}).get("status") != "ok":
        t0 = time.time()
        rc = subprocess.run([PY, os.path.abspath(__file__), "_pack", "--run", run, "--cell", cell, "--workers", str(a.workers), "--subset", a.subset],
                            env=dict(env, PYTHONPATH=code), cwd=code,
                            stdout=open(os.path.join(run, "logs", "subset_pack.log"), "a"), stderr=subprocess.STDOUT).returncode
        info["stages"]["subset_pack"] = dict(status="ok" if rc == 0 else "failed", rc=rc, seconds=round(time.time() - t0, 1))
        jsave(os.path.join(run, "run.json"), info)
        log(run, f"  subset pack: rc {rc}, {time.time() - t0:.0f} s")
        if rc != 0:
            raise SystemExit(f"subset pack failed (see {run}/logs/subset_pack.log)")
    # 4. metrics + compare (harness code)
    t0 = time.time()
    res = run_metrics(run)
    info["stages"]["metrics"] = dict(status="ok", seconds=round(time.time() - t0, 1))
    info["finished"] = time.strftime("%FT%T"); info["wall_s"] = round(time.time() - t_all, 1)
    jsave(os.path.join(run, "run.json"), info)
    if a.ref or a.ref_build:
        t0 = time.time()
        cmp_ = do_compare(run, a.ref, a.ref_build, a.ref_port, a.expect_bitwise, a.expect_same)
        shutil.copy2(os.path.join(run, f"compare_{cmp_['reference_label']}.json"), os.path.join(run, "compare.json"))
        info["stages"]["compare"] = dict(status="ok", seconds=round(time.time() - t0, 1), verdict=cmp_["verdict"])
        jsave(os.path.join(run, "run.json"), info)
    log(run, f"DONE {a.name} in {time.time() - t_all:.0f} s; placed mass fraction {res['c10']['placed_mass_fraction']:.4f}")


# ------------------------------------------------------------------ subset pack with the reference context (§3)
_CTX = {}


def _ctx_open(j):
    """raw context tile of run tile j (dict types / pos / quat / job_*), or None"""
    c = _CTX
    f = c["files"].get(int(j))
    if f is None:
        return None
    z = np.load(f)
    return {k: z[k] for k in ("types", "pos", "quat", "job_lib", "job_region", "job_n", "job_placed") if k in z.files}


def _ctx_load(spec):
    """'ctx:<j>:<max size class>' -> the context tile's instances in the run's library indices, limited to the jobs
    of size class <= max (instances are grouped by job in tile files)"""
    import wholecell.pack.stage as st
    from wholecell.pack.algorithms import job_rounds
    _, j, mr = spec.split(":"); j = int(j); mr = int(mr)
    z = _ctx_open(j)
    lut = _CTX["lut"]
    ty = lut[z["types"].astype(np.int64)]
    keep = ty >= 0
    nr = st._S["n_rounds"]
    if nr > 1 and mr < nr - 1 and "job_placed" in z:
        jp = z["job_placed"].astype(np.int64)
        if int(jp.sum()) == len(ty):
            jl = lut[z["job_lib"].astype(np.int64)]
            jr = _CTX["rmap"][z["job_region"].astype(np.int64)]
            reg = st._S["reg"]
            jk = reg.k_of[np.maximum(jr, 0)] if not reg.identity else np.maximum(jr, 0)
            cls = job_rounds(st._S["BR"][np.maximum(jl, 0)], reg.pack_class[jk] == 0, st._S["edges"])
            cls = np.where((jl < 0) | (jr < 0), 99, cls)
            keep &= np.repeat(cls, jp) <= mr
    return dict(types=ty[keep].astype(np.int32), pos=z["pos"][keep], quat=z["quat"][keep])


def _wrapped_pack_tile(ti):
    _CTX["cur"] = int(ti)
    return _CTX["orig_pack"](ti)


def _wrapped_tile_npz(j):
    """neighbour j of the tile being packed: the run's own file (subset tiles), else the context tile when it would
    have been packed before (colour order of the current round), limited to the size classes done by then"""
    import wholecell.pack.stage as st
    from wholecell.pack.algorithms import phase_order
    f = _CTX["orig_npz"](j)
    if f is not None or int(j) in _CTX["sub"] or int(j) not in _CTX["files"]:
        return f
    rnd = st._S["rnd"]
    order = phase_order(list(range(8)), rnd, bool(st._S["sw"].rounds["alternate"]))
    col = _CTX["colour"]
    mr = rnd if order.index(int(col[j])) < order.index(int(col[_CTX["cur"]])) else rnd - 1
    return None if mr < 0 else f"ctx:{int(j)}:{mr}"


class _NP:
    """numpy proxy for wholecell.pack.stage: np.load('ctx:...') returns a context tile"""

    def __getattr__(self, k):
        return getattr(np, k)

    def load(self, f, *a, **k):
        if isinstance(f, str) and f.startswith("ctx:"):
            return _ctx_load(f)
        return np.load(f, *a, **k)


def cmd_pack(a):
    """runs with the snapshot's code on sys.path (PYTHONPATH=<run>/_code)"""
    # numba threading as `python -m wholecell pack --step pack` sets it (wholecell/pack/stage.py run), BEFORE numba is
    # imported (plugins / pandas helpers may import it): forked workers single-threaded
    os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
    os.environ["NUMBA_NUM_THREADS"] = "1"
    import pandas as pd
    run, cell = a.run, a.cell
    sys.path.insert(0, os.path.join(run, "_code"))
    from wholecell import registry
    from wholecell.build import check_env
    from wholecell.config import load
    import wholecell.pack.stage as st
    check_env()
    cfg = load(os.path.join(run, "cell.yaml"))
    registry.load_plugins(cfg)
    ids, same = subset_ids(cell, run, a.subset)
    np.save(os.path.join(run, "subset_tiles.npy"), ids)
    M = pd.read_parquet(os.path.join(run, "pack", "tiles.parquet"))
    meta = json.load(open(os.path.join(run, "pack", "tiles_meta.json")))
    T = float(meta["TILE"]); lo0 = np.asarray(meta["lo0"], float)
    Tij = M[["ti", "tj", "tk"]].values
    key = {tuple(t): i for i, t in enumerate(Tij)}
    L = json.load(open(os.path.join(run, "pack", "libs.json")))["libs"]
    rn_run = list(json.load(open(os.path.join(run, "pack", "libs.json"))).get("rnames") or
                  [r["name"] for r in (cfg.get("pack") or {}).get("regions", [])] or meta.get("columns"))
    ctx = CONTEXT[cell]; cb = ctx["build"]
    # context tile files by run tile id, library / region maps
    if ctx["kind"] == "core":
        CM = pd.read_parquet(os.path.join(cb, "pack", "tiles.parquet"), columns=["tile", "ti", "tj", "tk"])
        cmeta = json.load(open(os.path.join(cb, "pack", "tiles_meta.json")))
        coff = np.round((np.asarray(cmeta["lo0"], float) - lo0) / T).astype(int)
        cL = json.load(open(os.path.join(cb, "pack", "libs.json")))
        clibs, crn = cL["libs"], cL.get("rnames") or rn_run
        files = {}
        for i, t in zip(CM.tile.values, CM[["ti", "tj", "tk"]].values):
            j = key.get(tuple(t + coff))
            if j is not None:
                files[j] = os.path.join(cb, "pack", "tiles", f"{int(i)}.npz")
    else:                                                 # yeast port: engine/work/yc3d/tiles_raw/<name>.npz
        raw = os.path.join(cb, "engine", "work", "yc3d")
        CM = pd.read_parquet(os.path.join(cb, "pack", "tiles.parquet"))
        CM = CM[~CM.prepass_only].reset_index(drop=True)
        coff = np.round(-lo0 / T).astype(int)
        clibs = [str(x) for x in np.load(os.path.join(raw, "lib_ids.npy"))]
        crn = rn_run                                       # port job_region = the core region table (A7 check)
        files = {}
        for nm, t in zip(CM.name.values, CM[["ti", "tj", "tk"]].values):
            j = key.get(tuple(t + coff))
            if j is not None:
                files[j] = os.path.join(raw, "tiles_raw", f"{nm}.npz")
    pos_ = {u: i for i, u in enumerate(L)}
    lut = np.array([pos_.get(u, -1) for u in clibs], np.int64)
    rpos = {u: i for i, u in enumerate(rn_run)}
    rmap = np.array([rpos.get(u, -1) for u in crn], np.int64)
    sub = set(int(i) for i in ids)

    def neigh(t):
        out = []
        for d in np.stack(np.meshgrid([-1, 0, 1], [-1, 0, 1], [-1, 0, 1], indexing="ij"), -1).reshape(-1, 3):
            if d.any():
                j = key.get(tuple(Tij[t] + d))
                if j is not None:
                    out.append(int(j))
        return out
    ctx_tiles = sorted({j for t in ids for j in neigh(int(t)) if j not in sub})
    have = [j for j in ctx_tiles if j in files and os.path.exists(files[j])]
    n_inst = n_drop = 0; dropped = {}
    for j in (have if clibs != L else []):                 # identical library order: nothing to drop
        z = np.load(files[j])
        ty = z["types"].astype(np.int64); n_inst += len(ty)
        bad = lut[ty] < 0
        if bad.any():
            n_drop += int(bad.sum())
            for u, c in zip(*np.unique(ty[bad], return_counts=True)):
                dropped[clibs[u]] = dropped.get(clibs[u], 0) + int(c)
    cinfo = dict(kind=ctx["kind"], build=cb, tiles_subset=len(ids), subset_ids_equal_frame=same, context_tiles=len(ctx_tiles),
                 context_tiles_found=len(have), context_instances=n_inst, dropped_instances=n_drop,
                 dropped_species=dict(sorted(dropped.items(), key=lambda x: -x[1])[:30]), libs_identical=clibs == L,
                 regions_identical=list(crn) == list(rn_run))
    jsave(os.path.join(run, "context.json"), cinfo)
    print("context:", json.dumps({k: v for k, v in cinfo.items() if k != "dropped_species"}), flush=True)
    _CTX.update(files={j: files[j] for j in have}, lut=lut, rmap=rmap)
    # the neighbour rule: wrap the tile packer (current tile) and the neighbour file lookup (module-level functions:
    # the worker pool pickles them by name)
    _CTX.update(orig_pack=st._pack_tile, orig_npz=st._tile_npz, sub=sub,
                colour=(Tij[:, 0] % 2) * 4 + (Tij[:, 1] % 2) * 2 + (Tij[:, 2] % 2))
    st._pack_tile, st._tile_npz, st.np = _wrapped_pack_tile, _wrapped_tile_npz, _NP()
    tl = os.path.join(run, "pack", "subset_tiles.npy")
    np.save(tl, ids)
    t0 = time.time()
    st.run(cfg, run, step="pack", workers=str(a.workers), tiles=tl)
    print(f"subset pack done in {time.time() - t0:.0f} s", flush=True)


# ------------------------------------------------------------------ metrics (§4), harness code
def _views(out, label, src_tiles, ids):
    """folder of symlinks <id>.npz -> the source tile files of the subset (nothing else)"""
    d = os.path.join(out, "views", label, "tiles")
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)
    for i, f in zip(ids, src_tiles):
        if f and os.path.exists(f):
            os.symlink(f, os.path.join(d, f"{int(i)}.npz"))
    return d


def seed_cache(bdir, cache_dir):
    """copy the validation library caches of a build that reads the SAME library folder (Jurkat: the reused B3 library)"""
    lib = os.path.realpath(os.path.join(bdir, "library"))
    for b in [CONTEXT["jurkat"]["build"], FRAME["jurkat"]]:
        if os.path.realpath(os.path.join(b, "library")) == lib:
            for n in ("lib_all.npz", "lib_packed.npz", "lib_preplaced.npz"):
                hits = sorted(glob.glob(os.path.join(b, "validate", "**", n), recursive=True), key=len)
                if hits and not os.path.exists(os.path.join(cache_dir, n)):
                    shutil.copy2(hits[0], os.path.join(cache_dir, n))
            return


def _yeast_extra(cfg, bdir):
    from wholecell.pack import switches
    ex = {}
    if os.path.isdir(os.path.join(bdir, "engine", "lib")):
        ex["lib_dir"] = os.path.join(bdir, "engine", "lib")
    sh = switches.effective(cfg.get("pack") or {}).get("shapes") or {}
    f = sh.get("file") if isinstance(sh, dict) else None
    if f:
        ex["pack_trees"] = f if os.path.isabs(f) else os.path.join(bdir, f)
    pj = jload(os.path.join(bdir, "geometry", "package.json"))
    if pj:
        ex["geom"] = pj["package_dir"]
    return ex


def subset_metrics(cfgp, bdir, ids, out, label, cell, port=None, masses=None):
    """result dict of the subset metrics for a core-layout build / run (bdir) or, with port=<s0 build>, for the port
    tiles of the same grid (bdir = a core run of the same tiles, used for tiles.parquet / geometry / regions).
    masses: {library name: Da} used for both sides of a comparison (the test's; yeast masses are library bead masses
    of each build's own library otherwise); every mass used is written to masses_<label>.json"""
    sys.path.insert(0, REPO)
    import pandas as pd
    from wholecell.config import load
    from wholecell.validate import common, placement, regions as R_
    from wholecell.validate.densities import _match
    from wholecell.validate.suite import seam_counts, seam_index
    cfg = load(cfgp)
    os.makedirs(out, exist_ok=True)
    common.CACHE_DIR = os.path.join(out, "cache", label)
    os.makedirs(common.CACHE_DIR, exist_ok=True)
    seed_cache(bdir, common.CACHE_DIR)
    extra = _yeast_extra(cfg, bdir) if cell == "yeast" else {}
    pdir = os.path.join(bdir, "pack")
    if port is None:
        src = [os.path.join(pdir, "tiles", f"{int(i)}.npz") for i in ids]
        L = json.load(open(os.path.join(pdir, "libs.json")))
        ov = dict(tiles_dir=_views(out, label, src, ids))
    else:
        raw = os.path.join(port, "engine", "work", "yc3d")
        S0 = pd.read_parquet(os.path.join(port, "pack", "tiles.parquet"))
        S0 = S0[~S0.prepass_only].reset_index(drop=True)
        meta = json.load(open(os.path.join(pdir, "tiles_meta.json")))
        off = np.round(-np.asarray(meta["lo0"], float) / float(meta["TILE"])).astype(int)
        M_ = pd.read_parquet(os.path.join(pdir, "tiles.parquet"), columns=["tile", "ti", "tj", "tk"])
        nm = {tuple(t + off): n for n, t in zip(S0.name.values, S0[["ti", "tj", "tk"]].values)}
        src = [os.path.join(raw, "tiles_raw", nm[tuple(t)] + ".npz") for t in M_.set_index("tile").loc[ids][["ti", "tj", "tk"]].values]
        L = dict(libs=[str(x) for x in np.load(os.path.join(raw, "lib_ids.npy"))])
        vd = os.path.join(out, "views", label); os.makedirs(vd, exist_ok=True)
        jsave(os.path.join(vd, "libs.json"), L)
        ppk = os.path.join(port, "pack")
        ov = dict(tiles_dir=_views(out, label, src, ids), libs=os.path.join(vd, "libs.json"),
                  preplaced_types=os.path.join(ppk, "prepass_types.npy"), preplaced_pos=os.path.join(ppk, "prepass_pos.npy"),
                  preplaced_quat=os.path.join(ppk, "prepass_quat.npy"), preplaced_off=os.path.join(ppk, "prepass_off.npy"),
                  preplaced_libs=os.path.join(ppk, "prepass_libs.json"), preplaced_groups=os.path.join(ppk, "prepass_groups.json"),
                  mrna_pts="/nonexistent", mrna_off="/nonexistent", mrna_meta="/nonexistent", counts="/nonexistent",
                  lib_dir=os.path.join(port, "engine", "lib"), pack_trees=os.path.join(raw, "pack_trees_15.npz"))
        extra = {k: v for k, v in extra.items() if k == "geom"}
    inp = common.Inputs(cfg, bdir, **{**extra, **ov})
    regs = R_.load_regions(cfg, inp); RN = R_.names(regs)
    kinds = [r["kind"] for r in regs]
    nreg = len(RN)
    mass, msrc = common.lib_masses(inp, bdir, nproc=8)
    Lnames = L["libs"]
    if masses:
        mass = np.array([float(masses.get(n, m)) for n, m in zip(Lnames, mass)])
        msrc += " (overridden by name with the test's masses)"
    used = {n: float(m) for n, m in zip(Lnames, mass)}
    tp = pd.read_parquet(os.path.join(pdir, "tiles.parquet")).set_index("tile")
    meta = json.load(open(os.path.join(pdir, "tiles_meta.json"))); T = float(meta["TILE"]); lo0 = np.asarray(meta["lo0"], float)
    # S-C10, S-P3/P5
    reqc = np.zeros(nreg); plac = np.zeros(nreg); reqm = np.zeros(nreg); plam = np.zeros(nreg)
    pairs = {}; packed = []; ninst = 0; seam = {r: [0, 0, 0] for r in SEAM_REGIONS[cell]}
    for i, f in zip(ids, src):
        if not os.path.exists(f):
            continue
        z = np.load(f); packed.append(int(i)); ninst += len(z["types"])
        jl = z["job_lib"].astype(np.int64); jr = z["job_region"].astype(np.int64); jn = z["job_n"].astype(np.int64); jp = z["job_placed"].astype(np.int64)
        np.add.at(reqc, jr, jn); np.add.at(plac, jr, jp); np.add.at(reqm, jr, jn * mass[jl]); np.add.at(plam, jr, jp * mass[jl])
        for l_, r_, n_, p_ in zip(jl, jr, jn, jp):
            k = (Lnames[l_], RN[r_]); q = pairs.setdefault(k, [0, 0]); q[0] += int(n_); q[1] += int(p_)
        for r in seam:
            if r in tp.columns and tp.at[int(i), r] / T ** 3 >= 0.95:
                ijk = tp.loc[int(i), ["ti", "tj", "tk"]].values.astype(np.int64)
                b_, n_ = seam_counts(z["pos"].astype(np.float64) - (lo0 + ijk * T), T)
                seam[r][0] += b_; seam[r][1] += n_; seam[r][2] += 1
    unth = np.zeros(nreg)
    cpath = inp.get("counts")
    if cpath and os.path.exists(cpath):
        cz = np.load(cpath); k = np.isin(cz["tile"], np.asarray(packed))
        np.add.at(unth, cz["region"][k].astype(np.int64), cz["n"][k] * mass[cz["lib"][k]])
    vol = {RN[i]: float(tp.loc[packed, RN[i]].sum()) / 1e9 for i in range(nreg) if RN[i] in tp.columns}
    per_region = {RN[i]: dict(kind=kinds[i], requested=int(reqc[i]), placed=int(plac[i]),
                              placed_frac=float(plac[i] / reqc[i]) if reqc[i] else None,
                              placed_mass_frac=float(plam[i] / reqm[i]) if reqm[i] else None,
                              placed_of_unthinned=float(plam[i] / unth[i]) if unth[i] else None,
                              requested_pg=reqm[i] / AV * 1e12, volume_um3=vol.get(RN[i])) for i in range(nreg) if reqc[i] or unth[i]}
    c10 = dict(tiles_packed=len(packed), instances=ninst, mass_source=msrc, placed_copies_fraction=float(plac.sum() / max(reqc.sum(), 1)),
               placed_mass_fraction=float(plam.sum() / max(reqm.sum(), 1e-9)),
               placed_of_unthinned=float(plam.sum() / unth.sum()) if unth.sum() else None, regions=per_region)
    PT = pd.DataFrame([(l_, r_, q[0], q[1]) for (l_, r_), q in pairs.items()], columns=["lib", "region", "requested", "placed"])
    with gzip.open(os.path.join(out, f"pairs_{label}.tsv.gz"), "wt") as f:
        PT.to_csv(f, sep="\t", index=False)
    # S-C7a / S-C8 / S-C9 (core placement step on the subset tiles)
    pj = os.path.join(out, f"placement_{label}.json")
    P = placement.run(cfg, bdir, inp, pj, os.path.join(out, f"placement_{label}_tiles.tsv"), overwrite=True, nproc=8, tiles=np.asarray(packed, np.int64))
    npk = P["n_packed"]
    c7 = dict(n_packed=npk, clash_rate_per_1e6=P["clash_rate_per_1e6"],
              packed_packed_per_1e6=P["clash_pairs"]["packed_packed"] / npk * 1e6 if npk else None,
              packed_prepass_per_1e6=P["clash_pairs"]["packed_prepass"] / npk * 1e6 if npk else None,
              packed_mrna_per_1e6=P["clash_pairs"]["packed_mrna"] / npk * 1e6 if npk else None,
              max_clash_depth_nm=P["max_clash_depth_nm"], worst_pair=P.get("worst_pair"),
              collision_shapes=str((P.get("settings") or {}).get("collision_shapes")))
    c8 = dict(containment_rate_per_1e6=P["containment_rate_per_1e6"],
              by_region={k: v for k, v in (P.get("containment_by_region") or {}).items() if v.get("bad")})
    c9 = dict(orientation_error_rate_per_1e6=P.get("orientation_error_rate_per_1e6"), n_membrane_instances=P.get("n_membrane_instances"),
              **{k: P.get(k) for k in ("tilt_median_deg", "tilt_p99_deg", "tilt_frac_lt_1deg", "tilt_n") if k in P})
    seam_out = {r: dict(index=seam_index(v[0], v[1], T, v[2]) if v[1] else None, n_tiles=v[2], n_band=int(v[0]), n_interior=int(v[1]))
                for r, v in seam.items()}
    # S-C11 densities (packed + pre-placed in the subset tiles; RNA excluded)
    dc = (cfg.get("validate", {}) or {}).get("densities", {}) or {}
    vol_i = R_.by_kind(regs, "volume", "virtual"); sh_i = R_.by_kind(regs, "shell")
    default_region = dc.get("default_region", RN[vol_i[0]])

    def to_volume(x):
        x = x.copy()
        virt_to = dc.get("virtual_to", default_region)
        for i in R_.by_kind(regs, "virtual"):
            host = virt_to.get(RN[i], default_region) if isinstance(virt_to, dict) else virt_to
            x[RN.index(host)] += x[i]; x[i] = 0.0
        s2v = dc.get("shell_to_volume", {}) or {}; split = dc.get("shell_split", {}) or {}
        for k in sh_i:
            if RN[k] in s2v:
                x[RN.index(s2v[RN[k]])] += x[k]
            for r_, f_ in (split.get(RN[k], {}) or {}).items():
                x[RN.index(r_)] += f_ * x[k]
        return x
    pk_v = to_volume(plam); rq_v = to_volume(reqm)
    pre_v = np.zeros(nreg); pre_counts = {}
    if inp.has("preplaced_off") and inp.has("preplaced_libs"):
        PL = json.load(open(inp["preplaced_libs"])); po = np.load(inp["preplaced_off"]); pt = np.load(inp["preplaced_types"], mmap_mode="r")
        cnt = np.zeros(len(PL), np.int64)
        for i in packed:
            cnt += np.bincount(np.asarray(pt[po[i]:po[i + 1]]).astype(np.int64), minlength=len(PL))
        gv = dc.get("preplaced_group_volume") or {}
        grp = json.load(open(inp["preplaced_groups"])) if (gv and inp.has("preplaced_groups")) else None
        lc = common.lib_cache(inp, PL, common.cache_path(bdir, "lib_preplaced.npz"), nproc=4)
        for i, (l, c) in enumerate(zip(PL, cnt)):
            if not c:
                continue
            m = float((masses or {}).get("preplaced:" + l, lc.bead_mass[lc.idx[l]]))
            used["preplaced:" + l] = m
            k_ = f"{grp[i]}[{grp[:i].count(grp[i])}]" if grp is not None else None
            if grp is not None and (k_ in gv or grp[i] in gv):
                g_ = gv.get(k_, gv.get(grp[i])); spl = g_ if isinstance(g_, dict) else {g_: 1.0}
            else:
                spl = {_match(l, dc.get("preplaced_volume"), default_region): 1.0}
            for reg_, f_ in spl.items():
                if reg_ in RN:
                    pre_v[RN.index(reg_)] += f_ * c * m
            pre_counts[l] = pre_counts.get(l, 0) + int(c)
    c11 = {}
    for i in vol_i:
        v = vol.get(RN[i]) or 0.0
        if v <= 0 or kinds[i] != "volume":
            continue
        rho = lambda m_: m_ / AV / (v * 1e-12) * 1e3
        c11[RN[i]] = dict(volume_um3=round(v, 4), mg_per_ml=rho(pk_v[i] + pre_v[i]), packed_mg_per_ml=rho(pk_v[i]),
                          requested_mg_per_ml=rho(rq_v[i] + pre_v[i]))
    # S-C12 ribosome NN (packed + pre-placed 80S centres inside the subset tiles)
    from scipy.spatial import cKDTree
    rlib = (cfg.get("validate", {}) or {}).get("ribosome_nn", {}).get("packed_lib")
    rnames = tuple(x for x in ([rlib] if rlib else []) + list(RIBO_LIBS))
    rib = [i for i, u in enumerate(Lnames) if u.startswith(rnames)]
    Pr = []
    for f in src:
        if os.path.exists(f):
            z = np.load(f); k = np.isin(z["types"], rib)
            if k.any():
                Pr.append(z["pos"][k].astype(np.float64))
    n_pk = sum(len(x) for x in Pr)
    if inp.has("preplaced_off") and inp.has("preplaced_libs"):
        prib = [i for i, u in enumerate(PL) if u.startswith(rnames)]
        pp = np.load(inp["preplaced_pos"], mmap_mode="r")
        for i in packed:
            t_ = np.asarray(pt[po[i]:po[i + 1]]); k = np.isin(t_, prib)
            if k.any():
                Pr.append(np.asarray(pp[po[i]:po[i + 1]])[k].astype(np.float64))
    c12 = dict(n=0)
    if Pr:
        PR = np.concatenate(Pr)
        if len(PR) > 1:
            d, _ = cKDTree(PR).query(PR, k=2)
            c12 = dict(n=int(len(PR)), n_packed=int(n_pk), nn_median=float(np.median(d[:, 1])), nn_p10=float(np.percentile(d[:, 1], 10)))
    jsave(os.path.join(out, f"masses_{label}.json"), used)
    res = dict(label=label, bdir=bdir, port=port, config=cfgp, tiles=[int(x) for x in ids], c10=c10, c7a=c7, c8=c8, c9=c9,
               seam=seam_out, c11=c11, c12=c12, preplaced_in_subset=dict(sorted(pre_counts.items())),
               placement_seconds=P.get("seconds"))
    if cell == "yeast":
        res["yeast"] = dict(mem_PEROX=per_region.get("mem_PEROX"), LD_peripheral=per_region.get("LD_peripheral"),
                            FAA4=PT[PT.lib.isin(YEAST_EXTRA_LIBS)].to_dict("records"))
    return res


_BROOT = re.compile(re.escape(R) + r"/[^/]+/(?:_dev_unify/step2/)?[^/\"\s]+")


def content_md5(p):
    """md5 of a stage output's content: arrays of npy / npz (sorted keys), text with build folder paths replaced"""
    if not p or not os.path.exists(p):
        return None
    h = hashlib.md5()
    if p.endswith(".npy"):
        a = np.load(p, mmap_mode="r", allow_pickle=False)
        h.update(str(a.dtype).encode()); h.update(str(a.shape).encode()); h.update(np.ascontiguousarray(a).tobytes())
    elif p.endswith(".npz"):
        with np.load(p, allow_pickle=True) as z:
            for k in sorted(z.files):
                a = z[k]; h.update(k.encode()); h.update(str(a.dtype).encode())
                h.update(a.tobytes() if a.dtype != object else repr(a.tolist()).encode())
    else:
        h.update(_BROOT.sub("{build}", open(p, errors="replace").read()).encode())
    return h.hexdigest()


def stage_metrics(run):
    """whole-cell stage-level outputs (S-stage)"""
    out = {}
    f = lambda *p: os.path.join(run, *p)
    out["md5"] = {k: content_md5(f(*p.split("/"))) for k, p in {
        "pack_table": "inventory/pack_table.tsv", "library_manifest": "library/library_manifest.tsv",
        "preplaced_manifest": "preplaced/preplaced_manifest.json", "preplaced_serial": "preplaced/serial/preplaced.npz",
        "preplaced_serial_lib_ids": "preplaced/serial/lib_ids.npy", "counts": "pack/counts.npz", "counts_thinned": "pack/counts_thinned.npz",
        "libs": "pack/libs.json", "tiles": "pack/tiles.parquet", "preplaced_types": "pack/preplaced_types.npy",
        "preplaced_pos": "pack/preplaced_pos.npy", "preplaced_off": "pack/preplaced_off.npy", "mrna_pts": "pack/mrna_pts.npy",
        "organisation_partners": "pack/organisation/partners.npz"}.items()}
    gp = jload(f("geometry", "package.json")) or {}
    out["geometry"] = dict(package=gp.get("package_id"), manifest_sha256=gp.get("manifest_sha256"))
    man = jload(f("preplaced", "preplaced_manifest.json")) or {}
    out["preplaced_files_md5"] = {m.get("name", str(i)): content_md5(m.get("npz")) for i, m in enumerate(man.get("instances", []))}
    L = jload(f("pack", "libs.json")) or {}
    rn = L.get("rnames")
    if os.path.exists(f("pack", "counts.npz")):
        cz = np.load(f("pack", "counts.npz"))
        r = np.bincount(cz["region"].astype(np.int64), weights=cz["n"])
        out["requested_per_region"] = {(rn[i] if rn and i < len(rn) else str(i)): int(v) for i, v in enumerate(r) if v}
        out["requested_total"] = int(cz["n"].sum())
    if os.path.exists(f("pack", "preplaced_types.npy")):
        PL = jload(f("pack", "preplaced_libs.json")) or []
        c = np.bincount(np.load(f("pack", "preplaced_types.npy")).astype(np.int64), minlength=len(PL))
        agg = {}
        for u, n in zip(PL, c):
            agg[u] = agg.get(u, 0) + int(n)
        out["preplaced_indexed_per_lib"] = agg
        out["preplaced_indexed_total"] = int(c.sum())
    lg = f("logs", "pack_index.log")
    if os.path.exists(lg):
        out["index_log"] = [l.strip() for l in open(lg) if "drop" in l or "indexed" in l][:20]
    th = jload(f("pack", "thinning.json"))
    if th:
        out["thinning_factors"] = th.get("factors_exact") or th.get("factors")
    if os.path.exists(f("pack", "tiles.parquet")):
        import pandas as pd
        M = pd.read_parquet(f("pack", "tiles.parquet"))
        out["tiles"] = len(M)
        out["region_measure_totals"] = {c: float(M[c].sum()) for c in M.columns if c not in ("tile", "ti", "tj", "tk")}
    return out


def run_metrics(run):
    info = jload(os.path.join(run, "run.json"))
    cell = info["cell"]
    ids = np.load(os.path.join(run, "subset_tiles.npy"))
    res = subset_metrics(os.path.join(run, "cell.yaml"), run, ids, os.path.join(run, "metrics"), "run", cell)
    res["stage"] = stage_metrics(run)
    res["context"] = jload(os.path.join(run, "context.json"))
    lt = os.path.join(run, "pack", "pack_log.tsv")
    if os.path.exists(lt):
        a_ = np.loadtxt(lt, ndmin=2)
        res["timing"] = dict(tile_cpu_s=float(a_[:, 5].sum()) if len(a_) else 0.0, stages=info.get("stages"))
    jsave(os.path.join(run, "result.json"), res)
    return res


# ------------------------------------------------------------------ comparison (§5)
def _binom(f, n):
    return 3.0 * np.sqrt(max(f * (1 - f), 0.0) / max(n, 1))


def compare(t, r, mode, expect_bitwise=False, expect_same=(), bitwise=None, stage_r=None):
    items = {}
    # S-C10 per region and overall (one-sided)
    rows = {}
    for name, x in t["c10"]["regions"].items():
        y = r["c10"]["regions"].get(name)
        if not y or x["placed_mass_frac"] is None or y["placed_mass_frac"] is None or x["requested"] < 100:
            continue
        tol = max(0.003, _binom(y["placed_frac"] or 1.0, x["requested"]))
        d = x["placed_mass_frac"] - y["placed_mass_frac"]
        rows[name] = dict(test=round(x["placed_mass_frac"], 5), ref=round(y["placed_mass_frac"], 5), delta=round(d, 5), tol=round(tol, 5),
                          requested=x["requested"], ok=bool(d >= -tol))
    ov_t, ov_r = t["c10"]["placed_mass_fraction"], r["c10"]["placed_mass_fraction"]
    items["S-C10_regions"] = dict(passes=all(v["ok"] for v in rows.values()), failing=[k for k, v in rows.items() if not v["ok"]], rows=rows)
    items["S-C10_overall"] = dict(test=ov_t, ref=ov_r, delta=ov_t - ov_r, passes=bool(ov_t - ov_r >= -0.003))
    # S-P3 / S-P5 per (type, region)
    import pandas as pd
    a_ = pd.read_csv(t["_pairs"], sep="\t").set_index(["lib", "region"]); b_ = pd.read_csv(r["_pairs"], sep="\t").set_index(["lib", "region"])
    u = a_.index.union(b_.index); a_ = a_.reindex(u, fill_value=0); b_ = b_.reindex(u, fill_value=0)
    tolc = lambda n, rel: np.maximum(3 * np.sqrt(np.abs(n)), rel * np.abs(n))
    big = b_.placed >= 10
    d5 = (a_.placed - b_.placed).abs()
    p5_big_fail = int((d5[big] > tolc(np.maximum(a_.placed, b_.placed)[big], 0.01)).sum())
    p5_small = float((d5[~big] <= 3).mean()) if (~big).any() else 1.0
    p5 = dict(pairs=int(len(u)), big=int(big.sum()), big_fail=p5_big_fail, small=int((~big).sum()), small_within_3=round(p5_small, 4),
              passes=bool(p5_big_fail == 0 and p5_small >= 0.99))
    d3 = (a_.requested - b_.requested).abs()
    p3 = dict(pairs=int(len(u)), fail=int((d3 > tolc(np.maximum(a_.requested, b_.requested), 0.005)).sum()),
              differ=int((d3 > 0).sum()))
    p3["passes"] = p3["fail"] == 0
    w = (a_.placed - b_.placed); worst = w.reindex(w.abs().sort_values(ascending=False).index).head(10)
    p5["worst"] = [dict(lib=k[0], region=k[1], test=int(a_.placed[k]), ref=int(b_.placed[k]), req_test=int(a_.requested[k]),
                        req_ref=int(b_.requested[k])) for k in worst.index]
    k10 = (a_.requested >= 10) & (b_.requested >= 10)
    fa = a_.placed[k10] / a_.requested[k10]; fb = b_.placed[k10] / b_.requested[k10]
    tolf = np.maximum(0.01, 3 * np.sqrt((fb * (1 - fb)).clip(lower=0) / a_.requested[k10]))
    okf = (fa - fb).abs() <= tolf
    pf_min = 0.935 if mode == "port" else 0.99            # addendum 2026-11-18 A3 (core vs port, calibrated on full builds)
    pf = dict(pairs=int(k10.sum()), within=round(float(okf.mean()), 4) if k10.any() else 1.0, threshold=pf_min,
              passes=bool(okf.mean() >= pf_min) if k10.any() else True,
              worst=[dict(lib=k[0], region=k[1], test=round(float(fa[k]), 4), ref=round(float(fb[k]), 4), n=int(a_.requested[k]))
                     for k in (fa - fb).abs().sort_values(ascending=False).head(8).index])
    if mode == "port":
        p5["cls"] = p3["cls"] = "report (core vs port requests, A2)"; pf["cls"] = "judged"
    else:
        p5["cls"] = p3["cls"] = "judged"; pf["cls"] = "report"
    items["S-P5"] = p5; items["S-P3"] = p3; items["S-P5_frac"] = pf
    # S-C7a / S-C8 / S-C9
    c7t, c7r = t["c7a"], r["c7a"]
    rt = lambda x: max(0.1 * abs(x or 0.0), 1.0)
    items["S-C7a"] = dict(test=dict(pp=c7t["packed_packed_per_1e6"], rate=c7t["clash_rate_per_1e6"], depth=c7t["max_clash_depth_nm"]),
                          ref=dict(pp=c7r["packed_packed_per_1e6"], rate=c7r["clash_rate_per_1e6"], depth=c7r["max_clash_depth_nm"]),
                          passes=bool((c7t["packed_packed_per_1e6"] or 0) <= 1.0 and (c7t["clash_rate_per_1e6"] or 0) <= (c7r["clash_rate_per_1e6"] or 0) + rt(c7r["clash_rate_per_1e6"])))
    c8t, c8r = t["c8"]["containment_rate_per_1e6"] or 0.0, r["c8"]["containment_rate_per_1e6"] or 0.0
    items["S-C8"] = dict(test=c8t, ref=c8r, passes=bool(c8t <= c8r + rt(c8r)), test_bad=t["c8"]["by_region"])
    o_t, o_r = t["c9"].get("orientation_error_rate_per_1e6") or 0.0, r["c9"].get("orientation_error_rate_per_1e6") or 0.0
    mt, mr = t["c9"].get("tilt_median_deg"), r["c9"].get("tilt_median_deg")
    items["S-C9"] = dict(sign_test=o_t, sign_ref=o_r, tilt_test=mt, tilt_ref=mr,
                         passes=bool(o_t <= o_r + rt(o_r) and (mt is None or mr is None or abs(mt - mr) <= 0.5)))
    sm = {k: dict(test=v["index"], ref=(r["seam"].get(k) or {}).get("index"), n_tiles=v["n_tiles"]) for k, v in t["seam"].items()}
    items["S-seam"] = dict(rows=sm, passes=all(v["test"] is None or abs(v["test"] - 1) <= 0.1 for v in sm.values()))
    # S-C11
    rows = {}
    for k, x in t["c11"].items():
        y = r["c11"].get(k)
        if not y or x["volume_um3"] < 0.01:
            continue
        d = x["mg_per_ml"] - y["mg_per_ml"]; dq = x["requested_mg_per_ml"] - y["requested_mg_per_ml"]
        rows[k] = dict(test=round(x["mg_per_ml"], 2), ref=round(y["mg_per_ml"], 2), delta=round(d, 2), delta_requested=round(dq, 2),
                       ok=bool(abs(d) <= 2.0 or abs(d - dq) <= 2.0), mass_consistent_only=bool(abs(d) > 2.0 and abs(d - dq) <= 2.0))
    items["S-C11"] = dict(rows=rows, passes=all(v["ok"] for v in rows.values()), failing=[k for k, v in rows.items() if not v["ok"]])
    nt, nr_ = t["c12"].get("nn_median"), r["c12"].get("nn_median")
    items["S-C12"] = dict(test=nt, ref=nr_, n_test=t["c12"].get("n"), n_ref=r["c12"].get("n"),
                          passes=bool(nt is None or nr_ is None or abs(nt - nr_) <= 0.5))
    if bitwise is not None:
        items["bitwise"] = dict(identical_tiles=bitwise[0], tiles=bitwise[1], cls="judged" if expect_bitwise else "report",
                                passes=bool(bitwise[0] == bitwise[1]) if expect_bitwise else True)
    if stage_r is not None:
        mt_, mr_ = t["stage"]["md5"], stage_r["md5"]
        same = {k: (mt_.get(k) == mr_.get(k)) for k in mt_ if mt_.get(k) or mr_.get(k)}
        miss = [k for k in expect_same if not same.get(k, False)]
        rq = {}
        for k in set(t["stage"].get("requested_per_region", {})) | set(stage_r.get("requested_per_region", {})):
            x, y = t["stage"].get("requested_per_region", {}).get(k, 0), stage_r.get("requested_per_region", {}).get(k, 0)
            if x != y:
                rq[k] = dict(test=x, ref=y, rel=round((x - y) / max(y, 1), 5))
        pre = {}
        for k in set(t["stage"].get("preplaced_indexed_per_lib", {})) | set(stage_r.get("preplaced_indexed_per_lib", {})):
            x, y = t["stage"].get("preplaced_indexed_per_lib", {}).get(k, 0), stage_r.get("preplaced_indexed_per_lib", {}).get(k, 0)
            if x != y:
                pre[k] = dict(test=x, ref=y)
        items["S-stage"] = dict(md5_equal=same, expected_same=list(expect_same), not_same=miss, passes=not miss,
                                cls="judged" if expect_same else "report", requested_per_region_diff=rq, preplaced_per_lib_diff=pre,
                                thinning=dict(test=t["stage"].get("thinning_factors"), ref=stage_r.get("thinning_factors")))
    judged = {k: v for k, v in items.items() if v.get("cls", "judged") == "judged"}
    fail = [k for k, v in judged.items() if not v["passes"]]
    return dict(mode=mode, verdict="PASS" if not fail else "FAIL", failing=fail, items=items)


def _bitwise(run_tiles, ref_tiles):
    n = same = 0
    for a, b in zip(run_tiles, ref_tiles):
        if not (os.path.exists(a) and os.path.exists(b)):
            continue
        n += 1
        x, y = np.load(a), np.load(b)
        if all(np.array_equal(x[k], y[k]) for k in ("types", "pos", "quat")):
            same += 1
    return same, n


def do_compare(run, ref=None, ref_build=None, ref_port=False, expect_bitwise=False, expect_same=()):
    info = jload(os.path.join(run, "run.json")); cell = info["cell"]
    t = jload(os.path.join(run, "result.json")) or run_metrics(run)
    t["_pairs"] = os.path.join(run, "metrics", "pairs_run.tsv.gz")
    ids = np.load(os.path.join(run, "subset_tiles.npy"))
    out = os.path.join(run, "compare"); os.makedirs(out, exist_ok=True)
    run_tiles = [os.path.join(run, "pack", "tiles", f"{int(i)}.npz") for i in ids]
    stage_r = None
    tm = jload(os.path.join(run, "metrics", "masses_run.json"))
    if tm is None:                                        # result.json of an older harness: recompute (writes the masses)
        t = run_metrics(run); t["_pairs"] = os.path.join(run, "metrics", "pairs_run.tsv.gz")
        tm = jload(os.path.join(run, "metrics", "masses_run.json"))
    if ref:                                               # another harness run (its metrics recomputed with the test's masses)
        label = "run_" + os.path.basename(ref.rstrip("/"))
        rids = np.load(os.path.join(ref, "subset_tiles.npy"))
        if (jload(os.path.join(ref, "run.json")) or {}).get("subset_kind", "frozen") != info.get("subset_kind", "frozen"):
            raise SystemExit(f"--ref {ref}: other subset kind than the test")
        r = subset_metrics(os.path.join(ref, "cell.yaml"), ref, rids, out, label, cell, masses=tm)
        r["_pairs"] = os.path.join(out, f"pairs_{label}.tsv.gz"); r["stage"] = (jload(os.path.join(ref, "result.json")) or {}).get("stage")
        ref_tiles = [os.path.join(ref, "pack", "tiles", f"{int(i)}.npz") for i in rids]
        stage_r = r.get("stage"); mode = "run"
    elif ref_port:                                        # port s0 tiles of the same grid, read with this run's tables
        label = "port_" + os.path.basename(ref_build)
        r = subset_metrics(os.path.join(run, "cell.yaml"), run, ids, out, label, cell, port=ref_build, masses=tm)
        r["_pairs"] = os.path.join(out, f"pairs_{label}.tsv.gz")
        ref_tiles = [os.path.join(out, "views", label, "tiles", f"{int(i)}.npz") for i in ids]
        mode = "port"
    else:                                                 # a core full build (its own config, read through a view)
        label = "build_" + os.path.basename(ref_build)
        rids, _ = subset_ids(cell, ref_build, info.get("subset_kind", "frozen"))
        cfgr = glob.glob(os.path.join(ref_build, "_code", "config", "cells", f"{info['cell_id']}.yaml")) or \
            glob.glob(os.path.join(ref_build, "_code", "config", "cells", "*.yaml"))
        r = subset_metrics(cfgr[0], ref_build, rids, out, label, cell, masses=tm)
        r["_pairs"] = os.path.join(out, f"pairs_{label}.tsv.gz")
        ref_tiles = [os.path.join(ref_build, "pack", "tiles", f"{int(i)}.npz") for i in rids]
        r["stage"] = stage_metrics(ref_build); stage_r = r["stage"]
        mode = "build"
    jsave(os.path.join(out, f"result_{label}.json"), {k: v for k, v in r.items() if k != "_pairs"})
    bw = _bitwise(run_tiles, ref_tiles)
    c = compare(t, r, mode, expect_bitwise, [x for x in (expect_same or "").split(",") if x] if isinstance(expect_same, str) else expect_same,
                bitwise=bw, stage_r=stage_r)
    c.update(test=run, reference=ref or ref_build, reference_label=label)
    jsave(os.path.join(run, f"compare_{label}.json"), c)
    log(run, f"compare vs {label}: {c['verdict']} (failing: {c['failing']}); bitwise {bw[0]}/{bw[1]}")
    return c


# ------------------------------------------------------------------ CLI
def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    r = sp.add_parser("run")
    r.add_argument("--name", required=True); r.add_argument("--config", required=True); r.add_argument("--commit", required=True)
    r.add_argument("--workers", type=int, default=24)
    r.add_argument("--reuse-library"); r.add_argument("--reuse-tiles"); r.add_argument("--reuse-upstream")
    r.add_argument("--reuse-parts", action="append", help="DONOR:part,part (geometry, tiles, inventory, library, engine, preplaced, _smk)")
    r.add_argument("--reuse-auto", help="DONOR run / build: reuse every upstream part whose code and config are unchanged, rebuild the rest")
    r.add_argument("--subset", choices=["frozen", "quick"], default="frozen",
                   help="frozen (56 / 58 tiles; acceptance) or quick (19 tiles; development, step2_quick_subsets_2026-11-18.md)")
    r.add_argument("--overlay", help="YAML file appended to cell.yaml (top-level keys not in the base write)")
    r.add_argument("--ref"); r.add_argument("--ref-build"); r.add_argument("--ref-port", action="store_true")
    r.add_argument("--expect-bitwise", action="store_true"); r.add_argument("--expect-same", default="")
    r.add_argument("--resume", action="store_true")
    p = sp.add_parser("_pack")
    p.add_argument("--run", required=True); p.add_argument("--cell", required=True); p.add_argument("--workers", type=int, default=24)
    p.add_argument("--subset", default="frozen")
    m = sp.add_parser("metrics"); m.add_argument("--run", required=True)
    c = sp.add_parser("compare")
    c.add_argument("--run", required=True); c.add_argument("--ref"); c.add_argument("--ref-build"); c.add_argument("--ref-port", action="store_true")
    c.add_argument("--expect-bitwise", action="store_true"); c.add_argument("--expect-same", default="")
    s = sp.add_parser("select"); s.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.cmd == "run":
        cmd_run(a)
    elif a.cmd == "_pack":
        cmd_pack(a)
    elif a.cmd == "metrics":
        run_metrics(a.run)
    elif a.cmd == "compare":
        do_compare(a.run, a.ref, a.ref_build, a.ref_port, a.expect_bitwise, a.expect_same)
    elif a.cmd == "select":
        bad = 0
        for cell in ("jurkat", "yeast"):
            ids = select(cell); h = hashlib.md5(ids.tobytes()).hexdigest()
            print(cell, len(ids), h, "as frozen" if h == SUBSET_MD5[cell] else "DIFFERS")
            bad += h != SUBSET_MD5[cell]
            q = select_quick(cell); hq = hashlib.md5(q.ref_tile.values.astype(np.int64).tobytes()).hexdigest()
            print(cell, "quick", len(q), hq, "as recorded" if hq == QUICK_MD5[cell] else "DIFFERS")
            bad += hq != QUICK_MD5[cell]
        if a.check and bad:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
