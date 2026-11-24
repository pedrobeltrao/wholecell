"""Repository text hygiene, a second:

    python tests/test_repo_text.py        (or pytest tests/test_repo_text.py)

  no_carriage_returns   no tracked text file contains a carriage return (a stray \\r in envs/yeast_x5.pip.txt split
                        "tabulate==0.10.0" in two for uv pip install)
  env_pip_lines         every envs/*.pip.txt line is a comment, blank or a pinned requirement name==version
"""
import glob
import os
import re
import subprocess

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _tracked():
    out = subprocess.run(["git", "-C", REPO, "ls-files", "-z"], capture_output=True, check=True).stdout
    return [p for p in out.decode().split("\0") if p]


def test_no_carriage_returns():
    bad = []
    for rel in _tracked():
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            continue
        data = open(path, "rb").read()
        if b"\0" in data[:8192]:
            continue  # binary
        if b"\r" in data:
            bad.append(rel)
    assert not bad, f"carriage returns in {bad}"


def test_env_pip_lines():
    req = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*==[^\s=]+$")
    for path in sorted(glob.glob(os.path.join(REPO, "envs", "*.pip.txt"))):
        for i, line in enumerate(open(path, newline=""), 1):
            s = line.rstrip("\n")
            assert not s.strip() or s.startswith("#") or req.match(s), f"{os.path.basename(path)}:{i}: {s!r}"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("repo text tests OK")
