# new in wholecell (downstream-only rerun on a finished upstream, 2026-11-19)
"""Run only the DOWNSTREAM rules (validate + export, incl. plugin steps of those stages) of a build into a NEW build
directory that reuses a donor build's finished upstream (stamp, frozen, geometry, inventory, library, preplaced, pack and
the upstream plugin steps), without re-packing:

    python scripts/build_downstream_reuse.py --donor <build> [--tag <tag>] [--cores 24] [--dry-run] [--launch]

  1. refuses a dirty tree; new build id = wholecell.plan (<date>-<commit7>-s<seed>-<tag>); refuses an existing dir
  2. <new>/_code = git worktree of HEAD (detached), as scripts/build_isolated.sh
  3. copies (cp -a: mtimes kept, never hardlinks: downstream steps must not be able to write into the donor) every
     donor file EXCEPT _code, .snakemake, validate/, export/, the outputs and logs of the downstream rules and their
     sentinels; absolute symlinks / text files naming the donor folder are re-pointed to the new folder (mtime kept)
  4. md5 manifest of the copied files (DOWNSTREAM_REUSE.md5) + DOWNSTREAM_REUSE.json (donor, donor commit, code commit,
     copied / excluded / re-pointed files, the rules snakemake will run)
  5. snakemake --dry-run with --rerun-triggers mtime: refuses unless ONLY downstream rules would run
  6. --launch: snakemake detached (setsid nohup), as build_isolated.sh, with --allowed-rules = the rules of the dry run
     (an upstream rule can never run, even if a frozen input's mtime changes); output <new>/snakemake.log

The donor is only read. BUILD_STAMP.json and config.resolved.json are the donor's (they describe the upstream that made
the reused outputs); the downstream code commit is in DOWNSTREAM_REUSE.json and in each validate MANIFEST.tsv.
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOWNSTREAM_STAGES = ("validate", "export")        # rule-name prefixes (core rules and plugin steps <stage>_<id>)
EXCLUDE_TOP = {"_code", ".snakemake", "validate", "export", "snakemake.log", "launch.log", "DOWNSTREAM_REUSE.json",
               "DOWNSTREAM_REUSE.md5"}
TEXT_EXT = (".json", ".tsv", ".csv", ".yaml", ".yml", ".txt")


def site():
    import yaml
    return yaml.safe_load(open(os.environ.get("WHOLECELL_SITE") or os.path.join(REPO, "config", "site.yaml")))


def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw).stdout


def downstream_outputs(cfg_cell, donor):
    """outputs (relative paths) of the downstream plugin steps + the engine files they write or link"""
    import yaml
    cfg = yaml.safe_load(open(cfg_cell))
    out = set()
    for s in cfg.get("plugin_steps") or []:
        if s.get("stage") in DOWNSTREAM_STAGES:
            out.update(s.get("outputs") or [])
            out.add(os.path.join("_smk", f"{s['stage']}_{s['id']}.done"))
            out.add(os.path.join("logs", f"{s['stage']}_{s['id']}.log"))
    # yeast engine layout written by the validate / export plugin steps (bridge.engine_pack_view, s4, s5, validate_all)
    eng = os.path.join("engine", "work")
    for f in ("yc3d/tiles_raw", "yc3d/geom", "yc3d/pack_summary.parquet", "yc3d/expected_totals.parquet", "yc3d/tile_counts.parquet",
              "yc3d/thinning.json", "yc3d/prepass.npz", "yc3d/lib_ids.npy", "yc3d/mrna_counts.npy", "yc3d/libpack_meta.npz",
              "yc3d/web", "yc3d/web_tmp", "s3_pack.log", "s4_validate.log", "s5_encode.log", "features.log",
              "features/protein_features.parquet", "validate_all.log", "results/wholecell3d_validation.json"):
        # (not removed: the work-seed folders validation/, sweeps/, interactions/, organisation/ and caches)
        out.add(os.path.join(eng, f))
    out.add(os.path.join("engine", "gpu.lock"))
    for f in os.listdir(os.path.join(donor, "logs")) if os.path.isdir(os.path.join(donor, "logs")) else []:
        if f.split("_")[0] in DOWNSTREAM_STAGES:
            out.add(os.path.join("logs", f))
    return sorted(out)


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", required=True)
    ap.add_argument("--cell")
    ap.add_argument("--tag", default="downstream")
    ap.add_argument("--cores", type=int, default=24)
    ap.add_argument("--dry-run", action="store_true", help="plan only: print the new folder and the excluded paths")
    ap.add_argument("--launch", action="store_true")
    ap.add_argument("--resume", action="store_true", help="continue an interrupted preparation of the same build id: entries "
                    "copied before the last one present are kept, the last one present is copied again")
    a = ap.parse_args()
    donor = os.path.realpath(a.donor.rstrip("/"))
    stamp = json.load(open(os.path.join(donor, "BUILD_STAMP.json")))
    cell = a.cell or stamp["cell_id"]
    S = site(); py, smk = S["python"], S["snakemake"]
    if sh(["git", "status", "--porcelain", "--untracked-files=no"], cwd=REPO).strip():
        raise SystemExit("refusing: working tree has uncommitted changes")
    commit = sh(["git", "rev-parse", "HEAD"], cwd=REPO).strip()
    plan = [l for l in sh([py, "-m", "wholecell.plan", "--cell", cell, "--build-tag", a.tag], cwd=REPO,
                          env=dict(os.environ, PYTHONHASHSEED="0", PYTHONPATH=REPO)).splitlines() if l.startswith("WHOLECELL_PLAN ")]
    P = json.loads(plan[-1][len("WHOLECELL_PLAN "):])
    bid = P["build_id"]; bdir = os.path.join(P["results_root"], P["cell_id"], bid)
    excl = downstream_outputs(os.path.join(REPO, "config", "cells", f"{cell}.yaml"), donor)
    print("new build:", bdir); print("excluded (downstream outputs):", len(excl))
    if a.dry_run:
        print("\n".join(excl)); return
    if os.path.exists(bdir) and not a.resume:
        raise SystemExit(f"refusing: {bdir} exists")
    t0 = time.time()
    if not a.resume:
        os.makedirs(bdir)
        sh(["git", "worktree", "add", "--detach", os.path.join(bdir, "_code"), commit], cwd=REPO)
    elif sh(["git", "-C", os.path.join(bdir, "_code"), "rev-parse", "HEAD"]).strip() != commit:
        raise SystemExit("--resume: the existing _code worktree is not at HEAD")
    # 3. copy (cp -a per top-level entry), then remove the excluded downstream outputs from the copy
    entries = [e for e in sorted(os.listdir(donor)) if e not in EXCLUDE_TOP]
    present = [e for e in entries if os.path.lexists(os.path.join(bdir, e))]
    keep = set(present[:-1]) if a.resume else set()          # the last entry present may be partial: copied again
    copied = []
    for e in entries:
        if e not in keep:
            if os.path.lexists(os.path.join(bdir, e)):
                sh(["rm", "-rf", os.path.join(bdir, e)])
            sh(["cp", "-a", os.path.join(donor, e), os.path.join(bdir, e)])
        copied.append(e)
    removed = []
    for r in excl:
        p = os.path.join(bdir, r)
        if os.path.lexists(p):
            sh(["rm", "-rf", p]); removed.append(r)
    # re-point absolute links and text files that name the donor folder (mtime kept)
    relinked, rewritten = [], []
    for root, dirs, files in os.walk(bdir):
        if root == bdir:
            dirs[:] = [d for d in dirs if d != "_code"]
        for n in files + [d for d in dirs if os.path.islink(os.path.join(root, d))]:
            p = os.path.join(root, n)
            if os.path.islink(p):
                t = os.readlink(p)
                if t.startswith(donor + "/"):
                    os.remove(p); os.symlink(bdir + t[len(donor):], p); relinked.append(os.path.relpath(p, bdir))
            elif n.endswith(TEXT_EXT) and os.path.getsize(p) < 200_000_000 and "/logs/" not in p + "/":
                s = open(p, "rb").read()
                if donor.encode() in s:
                    st = os.stat(p); before = hashlib.md5(s).hexdigest()
                    open(p, "wb").write(s.replace(donor.encode(), bdir.encode()))
                    os.utime(p, (st.st_atime, st.st_mtime)); rewritten.append(dict(file=os.path.relpath(p, bdir), donor_md5=before))
    # 4. md5 manifest of the copied regular files
    files = []
    for root, dirs, fs in os.walk(bdir):
        if root == bdir:
            dirs[:] = [d for d in dirs if d != "_code"]
        files += [os.path.join(root, f) for f in fs if not os.path.islink(os.path.join(root, f))]
    with ThreadPoolExecutor(8) as ex:
        sums = list(ex.map(md5, sorted(files)))
    with open(os.path.join(bdir, "DOWNSTREAM_REUSE.md5"), "w") as f:
        for p, m in zip(sorted(files), sums):
            f.write(f"{m}  {os.path.relpath(p, bdir)}\n")
    # 5. dry run: only downstream rules may run
    base = [smk, "--snakefile", os.path.join(bdir, "_code", "Snakefile"), "--directory", bdir, "--cores", str(a.cores),
            "--config", f"cell={cell}", f"build_id={bid}", f"build_dir={bdir}", "--rerun-triggers", "mtime", "--rerun-incomplete"]
    dr = subprocess.run(base + ["-n", "--quiet", "rules"], capture_output=True, text=True, cwd=bdir, env=dict(os.environ, PYTHONHASHSEED="0"))
    rules = sorted(set(re.findall(r"^\s*([A-Za-z0-9_]+)\s+\d+\s*$", dr.stdout, re.M)) - {"total", "job"})
    bad = [r for r in rules if r != "all" and r.split("_")[0] not in DOWNSTREAM_STAGES]
    info = dict(schema="wholecell.downstream_reuse/1", donor=donor, donor_commit=stamp.get("commit"), code_commit=commit,
                build_id=bid, build_dir=bdir, cell=cell, copied_top_level=copied, excluded_downstream=excl, removed_from_copy=removed,
                relinked=relinked, rewritten=rewritten, md5_manifest="DOWNSTREAM_REUSE.md5", n_files=len(files),
                rules_to_run=rules, dry_run_stdout=dr.stdout[-4000:], dry_run_stderr=dr.stderr[-2000:],
                note="BUILD_STAMP.json / config.resolved.json are the donor's (upstream provenance); downstream code = code_commit",
                prepared_s=round(time.time() - t0, 1), prepared=time.strftime("%Y-%m-%dT%H:%M:%S"))
    json.dump(info, open(os.path.join(bdir, "DOWNSTREAM_REUSE.json"), "w"), indent=1)
    print("rules to run:", rules)
    if dr.returncode != 0 or bad or not rules:
        raise SystemExit(f"refusing to launch: dry run rc {dr.returncode}, non-downstream rules {bad}\n{dr.stderr[-3000:]}")
    if a.launch:
        # only the downstream rules may run, also if a frozen input's mtime changes later (2026-11-19: a directory in
        # the read-only yeast tree was touched externally during the run -> mtime trigger re-ran the whole DAG)
        cmd = base + ["--allowed-rules"] + rules
        with open(os.path.join(bdir, "launch.log"), "w") as f:
            f.write(f"downstream reuse of {donor}\ncommit: {commit}\ncommand: {' '.join(cmd)}\nstarted: {time.strftime('%Y-%m-%dT%H:%M:%S')}\n")
        log = open(os.path.join(bdir, "snakemake.log"), "w")
        pr = subprocess.Popen(cmd, cwd=bdir, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                              start_new_session=True, env=dict(os.environ, PYTHONHASHSEED="0"))
        open(os.path.join(bdir, "launch.log"), "a").write(f"pid: {pr.pid}\n")
        print("launched pid", pr.pid)
    print(bdir)


if __name__ == "__main__":
    sys.exit(main())
