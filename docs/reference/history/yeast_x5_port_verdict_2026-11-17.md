# Yeast X5 port: final verdict (2026-11-17)

**Verdict: REPRODUCED (after accepted deviations).** Decided by the main session on 2026-11-17, in line with Pedro's
guidance to accept small non-content differences without extra runs; Pedro informed.

Basis (all files unchanged and kept):
- Frozen criteria `yeast_x5_reproduction_criteria.md`; addendum `yeast_x5_criteria_addendum_2026-11-17.md` (md5
  8ca1f90b…, frozen before seeds 1/2 were read); scorer fix `yeast_x5_scorer_fix_2026-11-17.md` (md5 adbfdd40…,
  post-hoc, implements the frozen §4 text "lists are report-only" — the literal reading is confirmed by the main
  session 2026-11-17; under it ~2,300 list-indexed scalars are report-only, all but 9 of which pass anyway);
  environment tolerances `yeast_x5_env_tolerances_2026-11-17.md`.
- Strict scored verdict (`yeast/20261116-aabe7b5-s0-x5repro/validate/compare_x5_seedspread_v2/verdict.json`):
  NOT REPRODUCED; must rows E 11 / T 11 / SS 3 / fail 1 (V2: 5 of 958 gating scalars).

Accepted deviations (the 5 V2 scalars):
1. ER periphery (2 scalars: `shells.ER_peripheral.median_nm`, `species.regions.ER_peripheral.species_below_0p8_of_region_share`;
   plus `regions[15].unplaced`, report-only under §4). Accepted by Pedro 2026-11-17 as practically equivalent. Diagnosis
   (`yeast/_dev_unify/er_periphery_diag/REPORT.md`): no port difference; the all-tile ER band median is 13.40–13.43 nm (port
   s0–s2) vs 13.43 nm (X5); X5's 400-tile sample sat 2.6 SD high; W4 is a different inventory, not a replicate.
2. PEROX halo (3 scalars `halos.membranes.PEROX.*`, NA in seed 2). Accepted by the main session 2026-11-17 as a non-content
   gap: the vendored `placement_checks.py` omits a membrane's halo when < 5e6 nm³ of cytosol next to it is sampled in its
   400 sampled tiles; the peroxisome sits at that cutoff (values present in s0/s1 and consistent with X5: 0.959 / 0.969 vs
   0.956). No re-run.

Also on record: PORT PATCH 4 (amphiphysin chirality bug, changes output) and PORT PATCH 5 (pandas 3 labels) were found
during the unification and are documented in `plugins/yeast/engine/PORT_PATCHES.md` with notes for the yeast team;
the reference for pre-placed builders 19–28 is port s0 + PATCH 4.
