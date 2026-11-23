# origin: jurkat/yc/pipeline/config.py (legacy workspace, 2026-11-15)
# ported, not run in phase 1 (library frozen). Changed: no directories are created at import; the library output folder
# is env WHOLECELL_LIB_OUT (default <this folder>/lib, never the frozen library); yeast-only entries removed.
"""Paths for the legacy library builders (wholecell/library/legacy). Env: WHOLECELL_LIB_OUT, WHOLECELL_LIB_DATA (the
legacy jurkat/yc/pipeline/data folder: uniprot_topology.tsv ...), WHOLECELL_PDB (PDB assembly cache),
WHOLECELL_AFDB_TAR, WHOLECELL_NPROC."""
import os

ROOT = os.path.dirname(os.path.abspath(__file__))
LIB = os.environ.get("WHOLECELL_LIB_OUT", os.path.join(ROOT, "lib"))
DATA = os.environ.get("WHOLECELL_LIB_DATA", os.path.join(ROOT, "data"))
PDB = os.environ.get("WHOLECELL_PDB", os.path.join(ROOT, "pdb"))
AFDB_TAR = os.environ.get("WHOLECELL_AFDB_TAR", "")
NPROC = int(os.environ.get("WHOLECELL_NPROC", 8))
