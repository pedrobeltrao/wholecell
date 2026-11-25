"""validate Inputs.lookup_input reads the preplaced manifest the preplaced stage writes (lists of {name, npz}), so the
ribosome and density checks use the build's own modules, not the validate.inputs fallback (once: v7-era files).

    python tests/test_validate_inputs.py        (or pytest tests/test_validate_inputs.py)
"""
import contextlib
import io
import json
import os
import pathlib
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from wholecell.validate import common


def _tmp():
    return pathlib.Path(tempfile.mkdtemp(prefix="wc_vinp_"))


def _out(fn):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn()
    return r, buf.getvalue()


def _build(tmp, manifest):
    os.makedirs(tmp / "preplaced", exist_ok=True)
    json.dump(manifest, open(tmp / "preplaced" / "preplaced_manifest.json", "w"))
    return str(tmp)


CFG = {"validate": {"inputs": {"polysomes": "/legacy/polysomes_cell1.npz", "mrna_nascent": "pre_nascent"}},
       "frozen_inputs": {"pre_nascent": {"path": "/frozen/mrna_nascent.npz"},
                         "pre_mrna_v1": {"path": "/frozen/mrna_cell1.npz"}}}


def test_current_manifest_format():
    b = _build(_tmp(), {"instances": [{"name": "polysomes", "npz": "/b/modules/polysomes.npz", "source": "serial"}],
                          "beads": [{"name": "mrna", "npz": "/b/modules/mrna.npz"}], "drops": []})
    inp = common.Inputs(CFG, b)
    assert inp.lookup_input("polysomes") == "/b/modules/polysomes.npz"
    assert inp.lookup_input("mrna") == "/b/modules/mrna.npz"


def test_older_manifest_forms():
    assert common.manifest_npz({"polysomes": {"path": "/a.npz"}}, "polysomes") == "/a.npz"
    assert common.manifest_npz({"modules": {"polysomes": "/b.npz"}}, "polysomes") == "/b.npz"
    assert common.manifest_npz({"instances": []}, "polysomes") is None


def test_fallback_warns_and_resolves_frozen_names():
    b = _build(_tmp(), {"instances": [], "beads": [], "drops": []})
    inp = common.Inputs(CFG, b)
    p, out = _out(lambda: inp.lookup_input("mrna_nascent"))
    assert p == "/frozen/mrna_nascent.npz" and "WARNING input 'mrna_nascent'" in out     # validate.inputs -> frozen input
    assert _out(lambda: inp.lookup_input("pre_mrna_v1"))[0] == "/frozen/mrna_cell1.npz"  # a frozen-input name directly
    assert _out(lambda: inp.lookup_input("nothing"))[0] is None


def test_no_manifest_no_warning():
    inp = common.Inputs(CFG, str(_tmp()))
    p, out = _out(lambda: inp.lookup_input("polysomes"))
    assert p == "/legacy/polysomes_cell1.npz" and "WARNING" not in out


def test_jurkat_config_has_no_legacy_polysomes():
    from wholecell import config as C
    raw = C.load_raw("jurkat")
    assert "polysomes" not in ((raw.get("validate", {}) or {}).get("inputs", {}) or {})
    assert "mrna" not in ((raw.get("validate", {}) or {}).get("inputs", {}) or {})


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn(); print(name, "OK", flush=True)
    print("validate input tests OK")
