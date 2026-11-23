"""wholecell: shared, config-driven whole-cell structural model pipeline (yeast, Jurkat and further cells).

Stages (each a subpackage; outputs per build in <results_root>/<cell>/<build_id>/<stage>/):
  geometry   cellgeom_v1 packages (read, validate, engines for the packer)
  inventory  data adapters -> inventory v0 -> assembly state -> ordered rule sets -> pack table
  library    structure library (per-species bead models + sphere trees)
  preplaced  pre-placed structures (chromatin, cytoskeleton, polysomes, NPC, lamina, ...; "prepass" before A0c)
  pack       tile packing (RSA), regions from the cell config
  validate   build-level validation (densities, localisation, placement integrity, ...)
  export     viewer tiles (tile format v2), simulation boxes, cells.json entries
  provenance per-protein ledger across builds
"""
__version__ = "0.1.0"
