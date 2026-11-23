"""Paths for the EM toolkit (no absolute paths in code; everything from config/site.yaml).

inp(rel): a derived input: the toolkit work folder if a regenerated copy exists there, else the legacy (read-only)
          {data_root}/data/derived/<rel>.
out(rel): a NEW output file in the toolkit work folder {geom_root}/_work/<cell>_em/ (env WHOLECELL_EM_WORK overrides);
          refuses to overwrite (write-once). Geometry steps write packages instead (toolkits/em/pkgstep.py).
"""
import os

from wholecell.config import load_site

SITE = load_site()
DATA_ROOT = SITE["data_root"]
GEOM_ROOT = SITE["geom_root"]
LEGACY_DERIVED = os.path.join(DATA_ROOT, "data", "derived")
EM_ROOT = os.path.join(DATA_ROOT, "data", "em", "jrc_jurkat-1")
EM_LABELS = os.path.join(EM_ROOT, "jrc_jurkat-1.n5", "labels") + "/"
EXTERNAL = os.path.join(DATA_ROOT, "external")
WORK = os.environ.get("WHOLECELL_EM_WORK") or os.path.join(GEOM_ROOT, "_work", "jurkat_em")


def inp(rel):
    w = os.path.join(WORK, rel)
    return w if os.path.exists(w) else os.path.join(LEGACY_DERIVED, rel)


def out(rel, overwrite=False):
    p = os.path.join(WORK, rel)
    if os.path.exists(p) and not overwrite:
        raise FileExistsError(f"{p} exists (toolkit outputs are write-once; use a new WHOLECELL_EM_WORK)")
    os.makedirs(os.path.dirname(p), exist_ok=True)
    return p
