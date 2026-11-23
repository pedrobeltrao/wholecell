# Yeast X5 seed-spread scoring: POST-HOC scorer fix and vacuole_surface re-run (2026-11-17)

Date: 2026-11-17. **Approved by Pedro 2026-11-17.** Separate dated file; the frozen criteria
`docs/reference/history/yeast_x5_reproduction_criteria.md` and the frozen addendum
`docs/reference/history/yeast_x5_criteria_addendum_2026-11-17.md` are **not edited** (addendum md5
`8ca1f90b64560d3e4e646a59546b339d`, unchanged; criteria md5 `dbed5ad461dee9187e7355be9469bd9c`, unchanged).

**This is a POST-HOC change, made after the seed-spread results had been seen** (verdict of
`<s0>/validate/compare_x5_seedspread/`, scorer `plugins/yeast/repro/seedspread_x5.py` as frozen in e2ed58f:
NOT REPRODUCED, must rows V1 and V2 failing). It changes (1) one input of V2 (a port defect, re-run) and (2) the
scorer code so that it implements the pre-registered text of addendum section 4. No margin, tolerance, file list, row
rule or verdict rule is changed. The first result stays on disk unchanged (`compare_x5_seedspread/`); the new one is
written next to it (`compare_x5_seedspread_v2/`).

## 1. Input change: validation/vacuole_surface.json (port defect, re-run)

`validation/vacuole_surface.py` (vendored, never edited) reads the model vacuole from
`{YC_WORK}/s2_geometry.log` ("vacuole: 1 compartment, V um3, surface S um2"). In X5 that line was written by stage 2.
In the port the geometry is made by the toolkit step (`toolkits/analytic/yeast.py`), which writes it to
`engine/work/geometry_package.log` / `<build>/logs/engine_geometry_package.log`; the script failed in every port seed
(`validate_all.log`: IndexError, "vacuole_surface failed"), and its 34 V2 scalars were NA in all three seeds.

- **Future builds**: `plugins/yeast/runner.py` `append_generator_vacuole_lines` (called by step `s2_geometry`) appends
  the generator's vacuole lines (lines starting with "vacuole" after the time stamp) to `engine/work/s2_geometry.log`
  after a marker line `# port (plugins/yeast/runner.py): vacuole lines of the geometry generator, copied from <log>`.
  Test: `tests/test_yeast_vacuole_log.py` (runs in conda_envs/yeast_x5; the vendored regex finds the numbers;
  idempotent; fallback to logs/engine_geometry_package.log). Not covered: the core (A9) yeast config
  `config/cells/yeast.yaml`, which reuses an existing geometry package and runs neither the generator nor stage 2,
  so its engine/work has no s2_geometry.log; vacuole_surface.py would still fail there (open item).
- **s0, s1, s2 (no rebuild)**: `plugins/yeast/repro/vacuole_surface_rerun.py --bdir <build>` ran ONLY
  vacuole_surface.py, with the build's own `_code` engine, the yeast_x5 python and the YC_* environment of
  `plugins/yeast/runner.env`, against a shadow work dir `<build>/validate/vacuole_surface_rerun/work`
  (s2_geometry.log = the build's own s2_geometry.log + marker + the 3 vacuole lines of the build's
  `logs/engine_geometry_package.log`; `yc3d/expected_totals.parquet` and `results/inventory/ingredients.csv` are
  symlinks to the build's files). The build's engine/work and _code were not written (the script's figure was
  redirected into the rerun folder). Provenance (inputs with md5, command, environment) in
  `<build>/validate/vacuole_surface_rerun/provenance.json`.

| Seed | output md5 | model: compartments, volume, area | area / sphere | vacuolar protein (fg) |
|---|---|---|---|---|
| s0 | 2349ac9a444f15c06b2c36e11dddbfac | 1, 5.84 um3, 16.2 um2 | 1.0330 | 46.2852 |
| s1 | 2349ac9a444f15c06b2c36e11dddbfac | 1, 5.84 um3, 16.2 um2 | 1.0330 | 46.2852 |
| s2 | 2349ac9a444f15c06b2c36e11dddbfac | 1, 5.84 um3, 16.2 um2 | 1.0330 | 46.2852 |
| X5 | 902d4583ef508d781cd761982dc63562 | 1, 5.84 um3, 16.2 um2 | 1.0330 | 46.2852 (last digits differ, ~1e-15 rel.) |

All 34 vacuole_surface scalars pass the seed-spread rule (matched X5+W4).

## 2. Scorer change: addendum section 4, lists are report-only

Frozen text (addendum section 4, last bullet): "Lists (per-bin vectors, `clash.top_type_pairs[*]` entries) and
strings are report-only (not gating), as in compare_x5." The e2ed58f scorer flattened lists of dicts into indexed paths
(`compare_x5.leaves`: `regions[15].unplaced`, `species.largest_species_below_half_of_region_share[10].requested`, ...)
and gated them, and treated an index present in only some seeds (a list of seed-dependent length) as NA.

Change (`plugins/yeast/repro/seedspread_x5.py`):
- A scalar whose path contains a list index (`[k]`) is **report-only**: scored as before and listed
  (`report_only = true`, verdict summaries `report_only`, details `report_only_not_passed`, scalars.tsv.gz), but never
  gating (`row_level` ignores it).
- Missing values: the text says "a missing port value makes the scalar NA" and "an NA scalar never passes". This is
  kept unchanged for every non-list scalar (e.g. `placement_checks.json:halos.membranes.PEROX.*`, absent in s2,
  stay NA and gating). No new rule for missing keys.
- `--extra FILE=REL`: V2 reads FILE (X5-relative) for each port seed from `<port build>/REL` (used:
  `validation/vacuole_surface.json=validate/vacuole_surface_rerun/work/validation/vacuole_surface.json`); md5 per seed
  recorded in details (`V2_extra_inputs`).
- `--accepted JSON`: Pedro's accepted deviations; the strict verdict is computed exactly as before; a separate field
  `verdict_after_accepted_deviations` treats the listed scalars as passed (nothing else changes).

Scope of the list rule (reading of "Lists ... are report-only" applied literally, as approved): it covers every
list-indexed path, not only the two examples. It made report-only 162 V1 scalars (`regions[*].*`, placed/expected per
region) and 2,158 V2 scalars (whole per-entry lists of `latest/localization_{dev,test}.json`, `singh2020_test.json`
`per_bait[*]`, `kinase_localization.json` `tests[*]`, `metabolism.json` `glycolysis[*]`, `literature_checks.json`
tables, `placement_checks.json` `largest_species...[*]`). Of these 2,320, 2,311 pass; the 9 not passed are
`regions[15].unplaced` (fail, V1) and 8 NA `largest_species_below_half_of_region_share[10..13].*`. Note for review:
the clause "as in compare_x5" could be read more narrowly (compare_x5 gated list-of-dict entries and only left
lists of scalars / lists of lists ungated, which the e2ed58f scorer already did); under that reading the fix would
change nothing. Pedro's approval (2026-11-17) names the list-index rule; it is implemented as approved, and no
de-gated scalar hides a failure other than the 9 listed.

## 3. Re-score (`<s0>/validate/compare_x5_seedspread_v2/`)

Command in `command.txt` there (yeast_x5 python, `--nproc 3`, X5 and default W4 dirs, the `--extra` above,
`--accepted accepted_deviations.json`); wall 4 min 32 s.

- **Strict verdict: NOT REPRODUCED** (must: E 11, T 11, SS 3, fail 1; report: E 2, T 1, SS 4). Failing must row: V2
  only (958 gating scalars, 5 not passed: `placement_checks.json:shells.ER_peripheral.median_nm` fail,
  `placement_checks.json:species.regions.ER_peripheral.species_below_0p8_of_region_share` fail, and 3 NA
  `placement_checks.json:halos.membranes.PEROX.{all_4_20nm_over_far, all_20_28nm_over_far, cytosol_region_4_20nm_over_far}`,
  absent in s2). V1 now SS (91/91 gating; `regions[15].unplaced` report-only).
- **Accepted deviations** (Pedro 2026-11-17: ER-periphery shift accepted as practically equivalent, diagnosis
  pending): `regions[15].unplaced` (V1, report-only anyway), `shells.ER_peripheral.median_nm`,
  `species.regions.ER_peripheral.species_below_0p8_of_region_share`. Scoring of these is unchanged.
- **Verdict after accepted deviations: NOT REPRODUCED** (V2 still fails on the 3 PEROX halo scalars, NA because s2
  has no PEROX entry in `halos.membranes`).
