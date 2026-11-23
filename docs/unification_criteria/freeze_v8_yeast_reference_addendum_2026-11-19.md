# Addendum to freeze_v8.md (2026-11-19): yeast reference build; Jurkat PM-topology gate

`freeze_v8.md` (md5 9bbd657167cf754afe6e725d3cd95d24) is frozen and not edited. This addendum fills the "pending" items of §1.

1. **Yeast no-regression reference** = `yeast/20261119-87a2018-s0-step2ds`. This is the combined Step-2 build, a downstream rerun of
   `yeast/20261119-b16bd09-s0-step2`. Its Step-2 verdict (`validate/compare_step2_v2`) is NOT ACCEPTED only on B5.L1
   SWAT_native tophit and 4 mass-weighted C1 rows. That build predates the yeast B5 overrides (6978b43). Per Pedro's DECISION of
   2026-11-19, these rows are re-scored on the freeze build under §4 E3, so no accepted-deviation file is passed. The queue
   job passes the path as `REFERENCE`.
2. **Jurkat gate `pm_topology` (E7): fix written, not waived** (Pedro, 2026-11-19, late). An earlier version of
   this addendum said the gate was open; that is superseded. The library step `pm_topology_fix`
   (`wholecell/library/topology.py`) runs after `adaptive_trees` and reads the table
   `plugins/jurkat/data/library/pm_topology_fix.tsv`.
   - **13 glyco-captured units, 168,347 copies:** rule `glyco_ectodomain_exoplasmic`. Each unit is rotated 180° about x,
     which is a proper rotation, so chirality is kept. Beads and every sphere tree (`st_c*`) are rotated; radii and kinds
     are unchanged.
   - **TRBC1_x2:** rule `tm_alternation`. The side comes from UniProt TOPO_DOM, with gaps filled by alternating sides
     across the TM helices. The unit is rotated only if the bead balance contradicts that side.
   - **The other flagged units** have no glyco evidence. They are left unchanged and listed as known residuals.
   - **Copies:** not touched; no inventory table changes.

   The gate is satisfied by config (the `library.steps` marker). G7 E7 is judged by the freeze check:
   - only the fix units differ from the reference library, and each is an exact rotation of the reference unit;
   - S0 orientation rule: none of the 13 glyco units is flagged;
   - TRBC1 is not flagged under alternation;
   - the residual flags are unchanged, and there are no new flags.

   The check writes `validate/compare_freeze_v8/e7_pm_topology.json` and `e7_pm_units.tsv`. The test is
   `tests/test_pm_topology.py`. The dry run applied the step to a copy of the Step-2 library; its result is in
   PLANNING_NOTES (2026-11-19, late).
