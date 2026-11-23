"""Build and upload the Hugging Face inputs dataset staged by scripts/deposit/stage_hf.py.

    HF_TOKEN=<write token> python upload_hf.py --staged <dir> --build <scratch dir> [--repo-id pbeltrao/wholecell-inputs]
    python upload_hf.py --staged <dir> --build <scratch dir> --dry-run      # build tars + final manifest, no upload

The token is read from the environment variable HF_TOKEN only (never from a file or the command line). Steps:
  1. tar shards of TARS.tsv are written to <build>/<hf path> (members checked against their md5 when known);
  2. the md5 of every dataset file is filled in, MANIFEST.tsv is written to <build>/ with README.md, fetch_inputs.py,
     TARS.tsv and EXCLUDED.tsv; the staged files are linked into <build>/ (no copy);
  3. unless --dry-run: the dataset repository is created (private) and <build>/ is uploaded with upload_large_folder.
Requires huggingface_hub (pip install huggingface_hub).
"""
import argparse
import hashlib
import os
import shutil
import tarfile

import pandas as pd


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--staged", required=True)
    ap.add_argument("--build", required=True)
    ap.add_argument("--repo-id", default="pbeltrao/wholecell-inputs")      # PLACEHOLDER repo id until the deposit exists
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--skip-member-md5", action="store_true", help="do not re-check tar members against their md5")
    a = ap.parse_args()
    F = pd.read_csv(os.path.join(a.staged, "MANIFEST.tsv"), sep="\t", dtype=str, keep_default_na=False)
    T = pd.read_csv(os.path.join(a.staged, "TARS.tsv"), sep="\t", dtype=str, keep_default_na=False)
    os.makedirs(a.build, exist_ok=True)
    for tar, g in T.groupby("tar"):
        dst = os.path.join(a.build, tar)
        if os.path.exists(dst):
            continue
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        with tarfile.open(dst + ".tmp", "w") as t:
            for r in g.itertuples():
                if r.md5 and not a.skip_member_md5 and md5_file(r.source) != r.md5:
                    raise SystemExit(f"{r.source}: md5 differs from its manifest; not uploading")
                t.add(r.source, arcname=r.member, recursive=False)
        os.replace(dst + ".tmp", dst)
        print("tar", tar, len(g), "members")
    for i, r in F.iterrows():
        p = os.path.join(a.build, r.hf_path)
        if r.kind == "file" and not os.path.lexists(p):
            os.makedirs(os.path.dirname(p), exist_ok=True)
            os.symlink(os.path.realpath(os.path.join(a.staged, "stage", r.hf_path)), p)
        got = md5_file(p)
        if r.md5 and got != r.md5:
            raise SystemExit(f"{r.hf_path}: md5 {got} != manifest {r.md5}; not uploading")
        F.at[i, "md5"] = got
        F.at[i, "bytes"] = str(os.path.getsize(p))
    F.to_csv(os.path.join(a.build, "MANIFEST.tsv"), sep="\t", index=False)
    for f in ("README.md", "TARS.tsv", "EXCLUDED.tsv", "SIZES.tsv", "fetch_inputs.py"):
        shutil.copyfile(os.path.join(a.staged, f), os.path.join(a.build, f))
    print(f"built {a.build}: {len(F)} files, {F.bytes.astype(int).sum() / 1e9:.2f} GB")
    if a.dry_run:
        return
    token = os.environ.get("HF_TOKEN")
    if not token:
        raise SystemExit("set HF_TOKEN in the environment (a write token for the dataset's owner)")
    from huggingface_hub import HfApi
    api = HfApi(token=token)
    api.create_repo(a.repo_id, repo_type="dataset", private=True, exist_ok=True)
    api.upload_large_folder(repo_id=a.repo_id, repo_type="dataset", folder_path=a.build)
    print(f"uploaded to https://huggingface.co/datasets/{a.repo_id} (private; make public after review)")


if __name__ == "__main__":
    main()
