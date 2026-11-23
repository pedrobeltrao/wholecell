# new in wholecell (2026-11-20)
"""numba thread count fixed before numba is imported (wholecell/numba_env.py); regression for the freeze yeast build
validate_suite failure 'Cannot set NUMBA_NUM_THREADS ... currently have 8, trying to set 32'."""
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CODE = r"""
import os, sys
sys.path.insert(0, %r)
from wholecell import numba_env
numba_env.ensure(numba_env.wanted("validate", "suite"))
import numba as nb, numpy as np           # a plugin importing numba after ensure
import wholecell.validate.stage           # stage module setdefault: must be a no-op now
@nb.njit(parallel=True)
def f(a):
    s = 0.
    for i in nb.prange(len(a)): s += a[i]
    return s
f(np.ones(100))
os.environ["UNRELATED"] = "1"             # an env change makes the next compile re-read the config
@nb.njit
def g(x): return x + 1
g(1)
from numba.core import config
assert int(config.NUMBA_NUM_THREADS) == 8 == nb.get_num_threads(), (config.NUMBA_NUM_THREADS, nb.get_num_threads())
print("ok")
"""


def _run(code, env):
    e = {k: v for k, v in os.environ.items() if not k.startswith("NUMBA_NUM")}
    e.update(env)
    return subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=e, timeout=600)


def test_ensure_before_import():
    r = _run(CODE % REPO, {})
    assert r.returncode == 0 and "ok" in r.stdout, r.stderr[-2000:]


def test_wanted():
    sys.path.insert(0, REPO)
    from wholecell import numba_env
    assert numba_env.wanted("pack", "tiles") == 1 and numba_env.wanted("pack", "pack") == 1
    assert numba_env.wanted("pack", "counts") == 8 and numba_env.wanted("validate", "suite") == 8


def test_wc_exports_threads():
    s = open(os.path.join(REPO, "workflow", "rules", "common.smk")).read()
    assert "NUMBA_NUM_THREADS={nt}" in s
    m = open(os.path.join(REPO, "wholecell", "__main__.py")).read()
    assert m.index("numba_env.ensure") < m.index("registry.load_plugins")


if __name__ == "__main__":
    for k, fn in list(globals().items()):
        if k.startswith("test_"):
            fn()
            print("ok", k)
