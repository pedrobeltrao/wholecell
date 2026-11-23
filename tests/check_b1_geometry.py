"""Criterion B1 (docs/reference/history/v7_reproduction_criteria.md): an evaluator vs legacy cell3d on jgeom_cell1_v7.npz.

    PYTHONHASHSEED=0 python tests/check_b1_geometry.py [--package {geom_root}/jurkat/em_cell1_v7] [--out b1.json]
        [--engine cellgeom_eval_v1|sdf_sections_v1] [--n-uniform 2000000] [--n-near 500000] [--n-pore 50000] [--threads 8]

--engine: cellgeom_eval_v1 (default; table-driven, rule table config/geometry/rules/sdf_sections_v1.yaml) or the legacy
sdf_sections_v1. Without --package the new G is assembled directly from the legacy G (split_legacy_G + build_G, default
Jurkat parameters and codes); with --package it is loaded through wholecell.geometry.cellgeom.load (checksums verified).
Core: toolkits/em/b1_check.py (also run by toolkits/em/convert_legacy.py into the package validation.json).
"""
import argparse
import json
import os
import sys
import tempfile

os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "wholecell_numba_cache"))
sys.dont_write_bytecode = True
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402

from wholecell.config import load_site  # noqa: E402
from wholecell.geometry.engine import get_engine  # noqa: E402
from toolkits.em import b1_check  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--package", default=None)
    ap.add_argument("--out", default=None)
    ap.add_argument("--n-uniform", type=int, default=2_000_000)
    ap.add_argument("--n-near", type=int, default=500_000)
    ap.add_argument("--n-pore", type=int, default=50_000)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--engine", default="cellgeom_eval_v1", choices=["cellgeom_eval_v1", "sdf_sections_v1"])
    a = ap.parse_args(argv)
    e = get_engine(a.engine)
    site = load_site()
    if a.package:
        from wholecell.geometry import cellgeom
        cell = cellgeom.load(a.package.format(**site), evaluator=a.engine, check_active_engine=False)
        Gn, pores, src = cell.G, cell.pores, f"cellgeom.load({a.package})"
    else:
        z = np.load(b1_check.LEGACY_NPZ.format(data_root=site["data_root"]))
        Gn = e.build_G(**e.split_legacy_G(z["G"])); pores = z["pores"]; src = "split_legacy_G + build_G"
    res = b1_check.run(Gn, pores, site["data_root"], n_uniform=a.n_uniform, n_near=a.n_near, n_pore=a.n_pore,
                       threads=a.threads, source=src, eng=e)
    res["rules"] = getattr(cell, "rules_path", None) if a.package else "default (build_G)"
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1)
    return 0 if res["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
