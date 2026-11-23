# B5 record: yeast curated / builder-anchored localisation overrides (2026-11-19)

Approved 2026-11-19 (main session, Pedro informed). The frozen B5 criteria, the revision and the addendum are unchanged.

Under preset `common` the yeast inventory uses three override rule sets (step `weight_overrides`, after `nc_model`):
- `b5.nuclear_silencing` (SIR2 nucleolus; SIR3, SIR4 nucleoplasm): nuclear silencing proteins that the common
  combiner moved to the cytosol / golgi_endosome.
- `b5.audit_curated` (28 proteins): the common top compartment has no GO-experimental or UniProt-experimental support
  while the legacy (X5) top compartment has both; legacy weights restored.
- `b5.builder_anchored` (13 proteins): proteins whose placement is set by a pre-placed structure (endocytic coat,
  amphiphysin, actin patches, spindle, chromatin sites, nuclear export) and whose structural share fell by >= 0.2 under
  common (e.g. RVS161/RVS167 cortex 0.74, which must follow the structure, not the N/C model); legacy weights restored.

Why: B5's common combiner turns CYCLoPs "Cytoplasm" calls into ~99 % cytosol for nuclear proteins and moves
structure-anchored proteins out of their structures (audit `wholecell_builds/yeast/_dev_unify/b5_audit/`: 702 proteins,
29 flagged). The rows were chosen without the benchmark sets (SWAT, hyperLOPIT). Effect: 0.21 % of protein mass;
§3d SWAT native 0.6682, hyperLOPIT 0.8483, image N/C 0.8345 (all pass; improve or hold).

Known residuals (one curated source only; left as common): RRB1, MTR4, SMT3, TIF3, ARC18. CAP1/CAP2 kept common
(the legacy cytosol-only pool limited capping protein to 1,091; common gives 1,553).
These overrides enter the yeast freeze build together with the orientation and nucleolus-cap packer fixes (3227713).
