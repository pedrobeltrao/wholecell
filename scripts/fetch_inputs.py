"""Download the inputs of the paper cells from the Hugging Face inputs dataset into the site roots.

    export WHOLECELL_SITE=/path/to/site.yaml          # data_root, yeast_root, afdb_root, geom_root must be set
    python scripts/fetch_inputs.py --cell jurkat      # or yeast, or all; --group reference adds the reference presets' inputs
    python scripts/fetch_inputs.py --cell yeast --list   # what would be fetched, and what must come from elsewhere

Every dataset file is listed in the dataset's MANIFEST.tsv with its target root and path (the layout of
docs/data_inputs.md); files are checked against their md5 and tar shards are unpacked into their target folder (an
unpacked shard is recorded in <root>/.fetch_markers/<hf_path>.done, outside the data folders; markers that older
versions of this script left inside a data folder are moved there). Inputs that are not in the dataset (EXCLUDED.tsv:
AlphaFold DB proteome tars, third-party tables whose terms do not allow redistribution) are listed at the end with
their source; --source downloads those that have a direct URL.
Needs huggingface_hub (pip install huggingface_hub). Set HF_TOKEN only if the dataset is private.
"""
import argparse
import collections
import hashlib
import os
import shutil
import sys
import tarfile
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from wholecell.config import load_site  # noqa: E402

DATASET = "pedrobeltrao/wholecell-inputs"


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


ROUTES = os.path.join(REPO, "config", "source_routes.tsv")
UA = "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130.0 Safari/537.36"
BIG = 2 * 10**9                                 # files above this are checked by size, not md5, when already present


def read_routes(path=None):
    import io
    import pandas as pd
    body = "".join(ln for ln in open(path or ROUTES) if not ln.startswith("#"))     # header comments only ("#" may be in a URL)
    return pd.read_csv(io.StringIO(body), sep="\t", dtype=str, keep_default_na=False)


def route_of(routes, target_root, target_path):
    """the route of an EXCLUDED row: exact target path, else a git_dir folder prefix, else None"""
    for q in routes.itertuples():
        if q.target_root != target_root:
            continue
        if q.target_path == target_path or (q.method == "git_dir" and target_path.startswith(q.target_path)):
            return q
    return None


def _download(url, out):
    """url -> out (browser-like User-Agent; some hosts refuse scripted requests without one)"""
    os.makedirs(os.path.dirname(out), exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=120) as r, open(out, "wb") as f:
        shutil.copyfileobj(r, f, 1 << 22)


def _state(dst, md5, nbytes):
    if not os.path.isfile(dst):
        return "missing"
    if nbytes and os.path.getsize(dst) != int(nbytes):
        return "differs"
    if md5 and (os.path.getsize(dst) <= BIG or not nbytes) and md5_file(dst) != md5:
        return "differs"
    return "ok"


def _install(tmp, dst, md5):
    """move a download into place if its md5 matches (else keep it next to the target as <dst>.download)"""
    if md5 and md5_file(tmp) != md5:
        os.replace(tmp, dst + ".download")
        return f"md5 differs from the paper builds' ({md5}); kept as {dst}.download"
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    os.replace(tmp, dst)
    return "fetched"


def _git_dir(q, tmp_root, cache):
    """sparse checkout of folder q.member of repo q.url at q.commit (once per route); the checked-out folder"""
    import subprocess
    key = (q.url, q.commit, q.member)
    if key not in cache:
        d = os.path.join(tmp_root, f"git{len(cache)}")
        run = lambda *c: subprocess.run(c, check=True, capture_output=True)
        run("git", "clone", "-q", "--filter=blob:none", "--no-checkout", "--sparse", q.url, d)
        run("git", "-C", d, "sparse-checkout", "set", q.member)
        run("git", "-C", d, "checkout", "-q", q.commit)
        cache[key] = os.path.join(d, q.member)
    return cache[key]


def fetch_excluded(E, site, a, cells):
    """inputs that are not in the dataset (EXCLUDED.tsv): check what is in place; with --source fetch them by their
    route (config/source_routes.tsv; a plain https route of EXCLUDED.tsv is used as is); print what is left to do"""
    import tempfile
    import pandas as pd
    routes = read_routes()
    if a.with_optional and cells & {"yeast"}:
        opt = routes[routes.source.str.startswith("OPTIONAL")]
        E = pd.concat([E, pd.DataFrame({"input": "optional", "cell": "yeast", "target_root": opt.target_root,
                                        "target_path": opt.target_path, "bytes": "", "md5": "", "why": opt.source,
                                        "route": opt.url})], ignore_index=True)
    print(f"\nNot in the dataset ({len(E)} rows; EXCLUDED.tsv; routes: config/source_routes.tsv)")
    done, todo, manual, other = collections.Counter(), [], {}, []
    tmp_root = tempfile.mkdtemp(prefix="fetch_src_")
    git_cache = {}
    try:
        for r in E.itertuples():
            root = site.get(r.target_root)
            q = route_of(routes, r.target_root, r.target_path)
            method = q.method if q is not None else ("url" if r.route.startswith("https://") else None)
            if method is None:                      # builds dataset, not needed, made locally, user's own environment
                other.append(r)
                continue
            if not root:
                todo.append((r, f"site key {r.target_root} not set"))
                continue
            dst = os.path.join(root, r.target_path)
            st = _state(dst, r.md5, r.bytes)
            if st == "ok":
                done["in place"] += 1
                continue
            if method == "manual":
                manual[dst] = (q, st)
                continue
            if not a.source:
                todo.append((r, f"{st}; fetch with --source ({method})"))
                continue
            tmp = os.path.join(tmp_root, "dl")
            try:
                if method == "url":
                    _download(q.url if q is not None else r.route, tmp)
                elif method == "tar_member":
                    import tarfile as _tf
                    _download(q.url, tmp + ".tar.gz")
                    with _tf.open(tmp + ".tar.gz") as t, open(tmp, "wb") as f:
                        shutil.copyfileobj(t.extractfile(q.member), f)
                elif method == "git_dir":
                    src = os.path.join(_git_dir(q, tmp_root, git_cache), r.target_path[len(q.target_path):])
                    shutil.copyfile(src, tmp)
                msg = _install(tmp, dst, r.md5)
            except Exception as e:                  # network, host or repository error: report and go on
                msg = f"download failed: {e}"
            if msg == "fetched":
                done["fetched"] += 1
            else:
                todo.append((r, msg))
    finally:
        shutil.rmtree(tmp_root, ignore_errors=True)
    for k, n in sorted(done.items()):
        print(f"  {n} {k}")
    for r in other:
        print(f"  not fetched ({r.input}): {r.target_path} -- {r.route}")
    if todo:
        print(f"\n  {len(todo)} not in place:")
        for r, msg in todo:
            print(f"    {r.target_root}/{r.target_path}: {msg}")
    if manual:
        print(f"\n  {len(manual)} files need a browser (journal sites refuse scripted downloads); download each from the "
              f"page, save it under the path shown, and rerun this script to check its md5:")
        for dst, (q, st) in manual.items():
            print(f"    [{st}] {dst}\n        page: {q.url}\n        file: {q.member}\n        source: {q.source}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", choices=["jurkat", "yeast", "all"], default="all")
    ap.add_argument("--group", choices=["required", "reference", "all"], default="required")
    ap.add_argument("--repo-id", default=DATASET)
    ap.add_argument("--revision", default=None, help="dataset revision (tag or commit) to pin")
    ap.add_argument("--list", action="store_true", help="only list what would be fetched")
    ap.add_argument("--source", action="store_true",
                    help="also download the inputs that are not in the dataset, by their routes (config/source_routes.tsv)")
    ap.add_argument("--with-optional", action="store_true",
                    help="--source: also the optional large source files (YeastRGB raw image archive, 9.7 GB)")
    ap.add_argument("--excluded-only", action="store_true",
                    help="skip the dataset files: only check (and with --source fetch) the inputs that are not in it")
    a = ap.parse_args()
    import pandas as pd
    from huggingface_hub import hf_hub_download
    site = load_site()
    dl = lambda f: hf_hub_download(a.repo_id, f, repo_type="dataset", revision=a.revision, token=os.environ.get("HF_TOKEN"))
    F = pd.read_csv(dl("MANIFEST.tsv"), sep="\t", dtype=str, keep_default_na=False)
    E = pd.read_csv(dl("EXCLUDED.tsv"), sep="\t", dtype=str, keep_default_na=False)
    cells = {"jurkat", "yeast", "testcell"} if a.cell == "all" else {a.cell}
    groups = {"required", "optional"} | ({"reference"} if a.group in ("reference", "all") else set())
    F = F[F.cell.isin(cells) & F.group.isin(groups)]
    E = E[E.cell.isin(cells)]
    print(f"{len(F)} dataset files, {F.bytes.astype(int).sum() / 1e9:.2f} GB")
    missing_roots = sorted({r for r in F.target_root if not site.get(r)})
    if missing_roots:
        raise SystemExit(f"site keys not set: {missing_roots} (WHOLECELL_SITE)")
    for r in F.itertuples():
        if a.excluded_only:
            break
        dst = os.path.join(site[r.target_root], r.target_path)
        if a.list:
            print(r.kind, r.hf_path, "->", dst)
            continue
        if r.kind == "file":
            if os.path.exists(dst) and (not r.md5 or md5_file(dst) == r.md5):
                continue
            p = dl(r.hf_path)
            if r.md5 and md5_file(p) != r.md5:
                raise SystemExit(f"{r.hf_path}: md5 mismatch after download")
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(p, dst)
        else:                                   # tar shard: unpack into the target folder
            # the "unpacked" marker lives outside the data folders: a file inside one changes its listing (the PDB
            # folder listing md5 of inventory_assembly) and its manifest check
            done = os.path.join(site[r.target_root], ".fetch_markers", r.hf_path + ".done")
            legacy = os.path.join(dst, f".fetched.{os.path.basename(r.hf_path)}")
            if os.path.exists(legacy):          # marker of an earlier fetch_inputs.py inside the data folder: move it
                os.makedirs(os.path.dirname(done), exist_ok=True)
                os.replace(legacy, done)
                print("moved marker", legacy, "->", done)
            if os.path.exists(done):
                continue
            p = dl(r.hf_path)
            if r.md5 and md5_file(p) != r.md5:
                raise SystemExit(f"{r.hf_path}: md5 mismatch after download")
            os.makedirs(dst, exist_ok=True)
            with tarfile.open(p) as t:
                t.extractall(dst, filter="data")
            os.makedirs(os.path.dirname(done), exist_ok=True)
            open(done, "w").close()
        print("ok", r.hf_path)
    fetch_excluded(E, site, a, cells)
    print("\nThen check: PYTHONHASHSEED=0 python -m wholecell frozen --cell <cell> --build <scratch dir>")


if __name__ == "__main__":
    main()
