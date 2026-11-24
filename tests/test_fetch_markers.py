"""fetch_inputs.py markers and folder-listing checks (synthetic data, a few seconds):

    python tests/test_fetch_markers.py        (or pytest tests/test_fetch_markers.py)

  fetch_markers        fetch_inputs.py against a fake local dataset (huggingface_hub replaced): a tar shard is unpacked,
                       its marker is written to <root>/.fetch_markers/ and not into the data folder; a second run
                       downloads nothing; a marker an older version left inside the folder is moved out, without a
                       new download
  listing_hidden       listing_md5 / visible_names ignore hidden files: the PDB folder listing md5 is the same with and
                       without .fetched.* markers (the inventory_assembly failure of a fetched data root)
  frozen_listing       verify_frozen checks a declared pdb_listing_md5: OK on the right folder, a failing row when a
                       file is missing, unchanged by a hidden marker
"""
import hashlib
import io
import os
import subprocess
import sys
import tarfile
import tempfile
import types

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

from wholecell.build import listing_md5, verify_frozen, visible_names  # noqa: E402

NAMES = ["1abc-assembly1.cif.gz", "2xyz-assembly1.cif.gz", "3def-assembly1.cif.gz"]


def _pdb_folder(root):
    d = os.path.join(root, "data", "structures", "pdb", "assemblies")
    os.makedirs(d, exist_ok=True)
    for n in NAMES:
        open(os.path.join(d, n), "wb").write(n.encode())
    return d


def _want():
    return hashlib.md5("\n".join(sorted(NAMES)).encode()).hexdigest()


def test_listing_hidden():
    with tempfile.TemporaryDirectory() as t:
        d = _pdb_folder(t)
        before = listing_md5(d)
        for i in range(8):
            open(os.path.join(d, f".fetched.assemblies.part00{i}.tar"), "w").close()
        assert listing_md5(d) == before == _want()
        assert visible_names(d) == sorted(NAMES)


def test_frozen_listing():
    with tempfile.TemporaryDirectory() as t:
        d = _pdb_folder(t)
        cfg = {"frozen_inputs": {"inv_pdb_assemblies": {"path": d}},
               "inventory": {"assembly": {"inputs": {"pdb_assemblies": "inv_pdb_assemblies", "pdb_listing_md5": _want()}}}}
        rows = {r["name"]: r for r in verify_frozen(cfg, strict=False)}
        assert rows["inv_pdb_assemblies:listing"]["ok"]
        open(os.path.join(d, ".fetched.assemblies.part000.tar"), "w").close()
        assert all(r["ok"] for r in verify_frozen(cfg, strict=False))
        os.remove(os.path.join(d, NAMES[0]))
        rows = {r["name"]: r for r in verify_frozen(cfg, strict=False)}
        assert not rows["inv_pdb_assemblies:listing"]["ok"]
        try:
            verify_frozen(cfg)
        except RuntimeError as e:
            assert "inv_pdb_assemblies:listing" in str(e)
        else:
            raise AssertionError("strict verify_frozen did not fail on a changed listing")


def _fake_dataset(t):
    """local stand-in of the HF dataset: MANIFEST.tsv (one tar shard), EXCLUDED.tsv, and the shard"""
    ds = os.path.join(t, "dataset")
    hf_path = "cc0/data_root/data/structures/pdb/assemblies.part000.tar"
    os.makedirs(os.path.dirname(os.path.join(ds, hf_path)), exist_ok=True)
    with tarfile.open(os.path.join(ds, hf_path), "w") as tf:
        for n in NAMES:
            ti = tarfile.TarInfo(n); ti.size = len(n)
            tf.addfile(ti, io.BytesIO(n.encode()))
    open(os.path.join(ds, "MANIFEST.tsv"), "w").write(
        "hf_path\tkind\tbytes\tmd5\tlicence\ttarget_root\ttarget_path\tcell\tgroup\n"
        f"{hf_path}\ttar\t{os.path.getsize(os.path.join(ds, hf_path))}\t\tCC0\tdata_root\tdata/structures/pdb/assemblies\tjurkat\trequired\n")
    open(os.path.join(ds, "EXCLUDED.tsv"), "w").write("input\tcell\ttarget_root\ttarget_path\twhy\troute\tmd5\n")
    return ds, hf_path


DRIVER = r'''
import os, sys, types, runpy
ds, calls = sys.argv[1], sys.argv[2]
def hf_hub_download(repo_id, f, **kw):
    open(calls, "a").write(f + "\n")
    return os.path.join(ds, f)
sys.modules["huggingface_hub"] = types.SimpleNamespace(hf_hub_download=hf_hub_download)
sys.argv = ["fetch_inputs.py", "--cell", "jurkat"]
runpy.run_path(sys.argv[3] if len(sys.argv) > 3 else os.environ["FETCH"], run_name="__main__")
'''


def _run(t, ds):
    calls = os.path.join(t, "calls.txt")
    if os.path.exists(calls):
        os.remove(calls)
    env = dict(os.environ, WHOLECELL_SITE=os.path.join(t, "site.yaml"), FETCH=os.path.join(REPO, "scripts", "fetch_inputs.py"))
    r = subprocess.run([sys.executable, "-c", DRIVER, ds, calls], env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return [ln for ln in open(calls).read().split("\n") if ln.endswith(".tar")], r.stdout


def test_fetch_markers():
    with tempfile.TemporaryDirectory() as t:
        root = os.path.join(t, "root"); os.makedirs(root)
        open(os.path.join(t, "site.yaml"), "w").write(f"data_root: {root}\n")
        ds, hf_path = _fake_dataset(t)
        d = os.path.join(root, "data", "structures", "pdb", "assemblies")
        marker = os.path.join(root, ".fetch_markers", hf_path + ".done")
        tars, _ = _run(t, ds)
        assert tars == [hf_path], tars
        assert sorted(os.listdir(d)) == sorted(NAMES), os.listdir(d)       # no marker inside the data folder
        assert os.path.exists(marker) and listing_md5(d) == _want()
        tars, _ = _run(t, ds)
        assert tars == [], tars                                             # unpacked once
        os.replace(marker, os.path.join(d, ".fetched.assemblies.part000.tar"))   # layout of the old script
        tars, out = _run(t, ds)
        assert tars == [] and "moved marker" in out, (tars, out)
        assert sorted(os.listdir(d)) == sorted(NAMES) and os.path.exists(marker)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("fetch marker tests OK")
