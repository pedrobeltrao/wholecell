# new in wholecell (2026-11-20: freeze yeast build validate_suite failure 'Cannot set NUMBA_NUM_THREADS ... currently
# have 8, trying to set 32', the same class as the pack fix 786ab2e)
"""One place that fixes numba's thread count for a process, BEFORE numba is imported.

Failure class: numba reads NUMBA_NUM_THREADS when it is imported and again on every compile (config.reload_config,
processed only when os.environ changed). If the variable is set or changed AFTER numba was imported (a plugin imports
numba first, then a stage module calls os.environ.setdefault), the thread pool can be launched with one value and a
later compile re-reads another, raising 'Cannot set NUMBA_NUM_THREADS to a different value once the threads have been
launched'. Rule: the value is in the environment before Python starts (workflow/rules/common.smk WC() exports it per
rule) and `ensure()` is called by `python -m wholecell` before any plugin import, so it never changes afterwards.
"""
from __future__ import annotations

import os
import sys

DEFAULT = 8


def wanted(stage, step):
    """threads of one stage step: pack tiles / pack workers are forked and single-threaded, everything else 8"""
    return 1 if stage == "pack" and step in ("tiles", "pack") else DEFAULT


def ensure(n=DEFAULT, force=False):
    """set NUMBA_NUM_THREADS (setdefault unless force) + workqueue layer; if numba is already imported, reload its config
    while no thread was launched, else require the launched value to equal the environment (clear error)"""
    os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
    if force:
        os.environ["NUMBA_NUM_THREADS"] = str(int(n))
    else:
        os.environ.setdefault("NUMBA_NUM_THREADS", str(int(n)))
    want = int(os.environ["NUMBA_NUM_THREADS"])
    if "numba" not in sys.modules:
        return want
    from numba.core import config as nconfig
    from numba.np.ufunc import parallel
    if parallel._is_initialized:
        if int(nconfig.NUMBA_NUM_THREADS) != want:
            raise RuntimeError(f"numba threads already launched with NUMBA_NUM_THREADS={nconfig.NUMBA_NUM_THREADS}, environment "
                               f"now {want}: export NUMBA_NUM_THREADS before Python starts (wholecell/numba_env.py)")
        return want
    nconfig.reload_config()
    return want
