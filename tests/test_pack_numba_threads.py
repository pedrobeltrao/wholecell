"""pack numba threading (wholecell/pack/stage.py numba_threading): numba imported before the pack step sets
NUMBA_NUM_THREADS=1 (plugins import it under `python -m wholecell`), threads launched without a recompile, then a compile
in a forked pool worker. Before the fix (environment set only) the worker raised 'Cannot set NUMBA_NUM_THREADS to a
different value once the threads have been launched' (yeast C3 thinning pilot, build 20261119-40d7967-s0-step2)."""
import os
import subprocess
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SCRIPT = r"""
import multiprocessing as mp, os, sys
import numba                                   # imported first, NUMBA_NUM_THREADS unset (config = CPU count)
MODE = sys.argv[1]
if MODE == "env_only":                         # the pack step's setup before the fix
    os.environ["NUMBA_THREADING_LAYER"] = "workqueue"; os.environ["NUMBA_NUM_THREADS"] = "1"
else:
    from wholecell.pack.stage import numba_threading
    numba_threading("pack")
from numba.np.ufunc import parallel
parallel._launch_threads()                     # e.g. a cached parallel kernel run in the parent: no config reload


def work(i):
    @numba.njit
    def f(x):
        return x + i
    return f(1)                                # a compile in the worker reloads the numba config


with mp.get_context("fork").Pool(2) as pool:
    print(sum(pool.map(work, range(4))))
"""


def _run(mode):
    env = {k: v for k, v in os.environ.items() if k not in ("NUMBA_NUM_THREADS", "NUMBA_THREADING_LAYER")}
    env["PYTHONPATH"] = REPO + os.pathsep + env.get("PYTHONPATH", "")
    env["PYTHONHASHSEED"] = "0"
    return subprocess.run([sys.executable, "-c", SCRIPT, mode], capture_output=True, text=True, env=env, timeout=600)


def test_numba_threading_forked_compile():
    if len(os.sched_getaffinity(0)) > 1:      # the mismatch needs a default thread count != 1
        r = _run("env_only")
        assert r.returncode != 0 and "Cannot set NUMBA_NUM_THREADS" in r.stderr, r.stderr[-2000:]
    r = _run("fixed")
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.strip() == "10"


if __name__ == "__main__":
    test_numba_threading_forked_compile()
    print("ok")
