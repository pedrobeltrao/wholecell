"""Download the inputs of the paper cells from the Hugging Face inputs dataset into the site roots.

    export WHOLECELL_SITE=/path/to/site.yaml          # data_root, yeast_root, afdb_root, geom_root must be set
    python scripts/fetch_inputs.py --cell jurkat      # or yeast, or all; --group reference adds the reference presets' inputs
    python scripts/fetch_inputs.py --cell yeast --list   # what would be fetched, and what must come from elsewhere

Every dataset file is listed in the dataset's MANIFEST.tsv with its target root and path (the layout of
docs/data_inputs.md); files are checked against their md5 and tar shards are unpacked into their target folder. Inputs
that are not in the dataset (EXCLUDED.tsv: AlphaFold DB proteome tars, third-party tables whose terms do not allow
redistribution) are listed at the end with their source; --source downloads those that have a direct URL.
Needs huggingface_hub (pip install huggingface_hub). Set HF_TOKEN only if the dataset is private.
The repository id below is a PLACEHOLDER until the deposit exists.
"""
import argparse
import hashlib
import os
import shutil
import sys
import tarfile
import urllib.request

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
from wholecell.config import load_site  # noqa: E402

DATASET = "pbeltrao/wholecell-inputs"      # PLACEHOLDER


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cell", choices=["jurkat", "yeast", "all"], default="all")
    ap.add_argument("--group", choices=["required", "reference", "all"], default="required")
    ap.add_argument("--repo-id", default=DATASET)
    ap.add_argument("--revision", default=None, help="dataset revision (tag or commit) to pin")
    ap.add_argument("--list", action="store_true", help="only list what would be fetched")
    ap.add_argument("--source", action="store_true", help="also download excluded inputs that have a direct URL")
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
            done = os.path.join(dst, f".fetched.{os.path.basename(r.hf_path)}")
            if os.path.exists(done):
                continue
            p = dl(r.hf_path)
            if r.md5 and md5_file(p) != r.md5:
                raise SystemExit(f"{r.hf_path}: md5 mismatch after download")
            os.makedirs(dst, exist_ok=True)
            with tarfile.open(p) as t:
                t.extractall(dst, filter="data")
            open(done, "w").close()
        print("ok", r.hf_path)
    print(f"\nNot in the dataset ({len(E)} rows; EXCLUDED.tsv): get them from the source into the same path")
    for r in E.itertuples():
        dst = os.path.join(str(site.get(r.target_root) or f"<{r.target_root}>"), r.target_path)
        print(f"  {r.input}: {dst}\n      why: {r.why}\n      route: {r.route}")
        if a.source and r.route.startswith("https://") and not os.path.exists(dst):
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            urllib.request.urlretrieve(r.route, dst)
            if r.md5 and md5_file(dst) != r.md5:
                print(f"      WARNING: md5 of the download differs from the paper builds' ({r.md5})")
    print("\nThen check: PYTHONHASHSEED=0 python -m wholecell frozen --cell <cell> --build <scratch dir>")


if __name__ == "__main__":
    main()
