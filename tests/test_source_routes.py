"""config/source_routes.tsv and fetch_inputs.py --source (offline: file:// URLs and a local git repository, seconds):

    python tests/test_source_routes.py        (or pytest tests/test_source_routes.py)

  routes_table     the routes file parses; methods are known; url/tar_member/git_dir rows have what their method needs;
                   git_dir targets are folder prefixes; no two routes share a target
  fetch_excluded   every method against a fake EXCLUDED table: url, tar_member and git_dir downloads are installed when
                   their md5 matches; a download whose md5 differs is kept as <file>.download and reported; a manual row
                   is listed with its page and file; a second run finds everything in place and downloads nothing
"""
import argparse
import hashlib
import importlib.util
import io
import os
import subprocess
import sys
import tarfile
import tempfile

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)
_spec = importlib.util.spec_from_file_location("fetch_inputs", os.path.join(REPO, "scripts", "fetch_inputs.py"))
FX = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(FX)


def _md5(b):
    return hashlib.md5(b).hexdigest()


def test_routes_table():
    R = FX.read_routes()
    assert len(R) >= 19 and list(R.columns) == ["target_root", "target_path", "method", "url", "member", "commit", "source"]
    assert set(R.method) <= {"url", "tar_member", "git_dir", "manual"}, set(R.method)
    assert not R.duplicated(["target_root", "target_path"]).any()
    for q in R.itertuples():
        assert q.url.startswith("https://"), q
        if q.method in ("tar_member", "git_dir", "manual"):
            assert q.member, q
        if q.method == "git_dir":
            assert q.target_path.endswith("/") and len(q.commit) == 40, q


def test_fetch_excluded():
    with tempfile.TemporaryDirectory() as t:
        src, root = os.path.join(t, "src"), os.path.join(t, "root")
        os.makedirs(src)
        good, bad = b"table v1\n", b"table v2\n"
        open(os.path.join(src, "a.xls"), "wb").write(good)
        open(os.path.join(src, "b.xls"), "wb").write(bad)            # the source changed since the paper builds
        with tarfile.open(os.path.join(src, "t.tar.gz"), "w:gz") as tf:
            for n, b in (("x.fa", b">x\nACGU\n"), ("y.fa", b">y\n")):
                ti = tarfile.TarInfo(n); ti.size = len(b); tf.addfile(ti, io.BytesIO(b))
        g = os.path.join(t, "repo")
        os.makedirs(os.path.join(g, "peaks")); os.makedirs(os.path.join(g, "other"))
        beds = {f"p{i}.bed": f"chrI\t{i}\t{i + 9}\n".encode() for i in range(3)}
        for n, b in beds.items():
            open(os.path.join(g, "peaks", n), "wb").write(b)
        open(os.path.join(g, "other", "big.bin"), "wb").write(b"0" * 100)
        run = lambda *c: subprocess.run(c, cwd=g, check=True, capture_output=True)
        run("git", "init", "-q"); run("git", "add", ".")
        run("git", "-c", "user.name=t", "-c", "user.email=t@t", "commit", "-q", "-m", "x")
        commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=g, capture_output=True, text=True).stdout.strip()
        routes = os.path.join(t, "routes.tsv")
        rows = [("d/a.xls", "url", f"file://{src}/a.xls", "", ""), ("d/b.xls", "url", f"file://{src}/b.xls", "", ""),
                ("d/x.fa", "tar_member", f"file://{src}/t.tar.gz", "x.fa", ""),
                ("d/peaks/", "git_dir", f"file://{g}", "peaks", commit),
                ("d/m.xlsx", "manual", "https://example.org/article", "mmc9.xlsx", "")]
        open(routes, "w").write("# test\ntarget_root\ttarget_path\tmethod\turl\tmember\tcommit\tsource\n" +
                                "".join(f"yeast_root\t{p}\t{m}\t{u}\t{mem}\t{c}\tsrc\n" for p, m, u, mem, c in rows))
        E = pd.DataFrame([dict(input="i", cell="yeast", target_root="yeast_root", target_path=p, bytes=str(len(b)),
                               md5=_md5(b), why="w", route="fetch from source: x")
                          for p, b in [("d/a.xls", good), ("d/b.xls", good), ("d/x.fa", b">x\nACGU\n"),
                                       ("d/m.xlsx", b"manual")] + [(f"d/peaks/{n}", b) for n, b in beds.items()]])
        FX.ROUTES = routes
        a = argparse.Namespace(source=True, with_optional=False)
        out = io.StringIO(); sys.stdout, so = out, sys.stdout
        try:
            FX.fetch_excluded(E, {"yeast_root": root}, a, {"yeast"})
        finally:
            sys.stdout = so
        txt = out.getvalue()
        d = os.path.join(root, "d")
        assert open(os.path.join(d, "a.xls"), "rb").read() == good
        assert open(os.path.join(d, "x.fa"), "rb").read() == b">x\nACGU\n"
        assert {n: open(os.path.join(d, "peaks", n), "rb").read() for n in beds} == beds
        assert not os.path.exists(os.path.join(d, "b.xls")) and os.path.exists(os.path.join(d, "b.xls.download"))
        assert "5 fetched" in txt and "md5 differs" in txt and "https://example.org/article" in txt and "mmc9.xlsx" in txt, txt
        out = io.StringIO(); sys.stdout = out
        try:
            FX.fetch_excluded(E, {"yeast_root": root}, argparse.Namespace(source=False, with_optional=False), {"yeast"})
        finally:
            sys.stdout = so
        txt = out.getvalue()
        assert "5 in place" in txt and "fetched" not in txt, txt


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("source route tests OK")
