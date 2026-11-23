# new in the release staging (2026-11-23): site paths for the development / unification check scripts in scripts/dev.
"""These scripts were written against the development workspace (the legacy project tree under data_root, the
development builds under results_root, the yeast model tree under yeast_root). Their former absolute paths are read
from the site file instead (config/site.yaml or WHOLECELL_SITE; env overrides as in wholecell/config.py):

  VC       legacy workspace root = site data_root (env WHOLECELL_DEV_WORKSPACE overrides)
  RESULTS  site results_root (the development builds the checks compare against live there)
  YEAST    site yeast_root (READ-ONLY yeast model tree)
  PY, SMK  site python / snakemake
  GIT      the git checkout used for worktrees (env WHOLECELL_DEV_REPO, default: this checkout)

They reproduce development-time checks only and need that tree (see docs/release_notes.md)."""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if REPO not in sys.path:
    sys.path.insert(0, REPO)
from wholecell.config import load_site  # noqa: E402

SITE = load_site()
VC = os.environ.get("WHOLECELL_DEV_WORKSPACE") or str(SITE.get("data_root") or "<data_root>")
RESULTS = str(SITE.get("results_root") or os.path.join(VC, "data", "derived", "wholecell_builds"))
YEAST = str(SITE.get("yeast_root") or "<yeast_root>")
PY = str(SITE.get("python") or sys.executable)
SMK = str(SITE.get("snakemake") or "snakemake")
GIT = os.environ.get("WHOLECELL_DEV_REPO") or REPO
