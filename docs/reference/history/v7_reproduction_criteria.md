# Jurkat v7 reproduction: pre-declared criteria

Written 2026-11-15, **before** any wholecell Jurkat build was run. Committed to git before the run; its md5 is
recorded in `PLANNING_NOTES.md` (DECISION 2026-11-15, wholecell). Changes after the run are allowed only as a dated
addendum below the line at the end, never by editing the criteria.

## 0. What is compared

- **Reference (v7)**: the published pack v7. Tiles `data/derived/wc/tiles_v7/` (archived copy of `wc/tiles`),
  `wc/tiles_v7.parquet`, `wc/libs_v7.json`, `wc/summary_v7.json` (md5 6da1e0ee…), `wc/counts.npz` (md5 c8c1ade1…,
  the v7 counts), `ingredients_v7.tsv` (md5 d1a0099c…), `jgeom_cell1_v7.npz` (md5 0eea64ff…), validation files
  `validation/*_v7.*`. All read-only.
- **Port (P)**: the wholecell build of `config/cells/jurkat.yaml` with the legacy v7 seeds, into a new folder under
  `data/derived/wholecell_builds/jurkat/<build_id>/`.
- **Inputs declared frozen** (checksums in the config, verified before use; identical by construction): chromatin v8
  rep5 pre-pass (`chromatin/v8/prepass/chromatin_cell1.npz`, md5 cd4c3b4c…), nascent RNA, polysomes, cytoskeleton
  (traced MTs; `prepass_v7/cytoskeleton_cell1.npz` md5 1523616f…, `mt_stage3_b_centerlines.npz` md5 84b33cf2…), NPC,
  mRNA + drop list, the structure library (16,246 npz), the NE HMFF mesh, the pack-dependent rule tables (v4/v5
  `loc_compare`, `tier2_v5`, `tm_route_libs`) and the hold-out lists. Lamina: rebuilt by the ported `jlamina`
  (seed 20261112) **and** compared with the frozen `prepass_v7/lamina_cell1.npz` (md5 f6a75524…).
- **Validation code**: wherever the ported validation code is run, it is run on **both** the v7 tiles and the port
  tiles (same code, same geometry, same files), so differences measure the build, not the code. Known v7 defects
  (`ribosome_nn_v7.json` = v6 copy; v7 fluorescence and `densities_v7` used v6 inputs) are therefore not counted
  against the port; agreement of the ported validation code with the stored v7 numbers is reported separately (§4).

## 1. Expected outcome and the two-level verdict

With identical inputs, legacy seeds, the same numba (0.68.1) and the same phase order, the port is **expected to be
bit-identical** to v7 at every level (the legacy pack was never rebuilt, so this is itself a test).

- **Level E (exact)**: the artefact is identical (md5 / array equality).
- **Level T (tolerance)**: if not exact, the artefact is within the tolerances below.

Verdict: **REPRODUCED** if every *must* row passes at level E or T; **REPRODUCED WITH DEVIATIONS** if all *must*
rows pass but some *report* rows fail; **NOT REPRODUCED** if any *must* row fails. Each failure is explained.

## 2. Build artefacts

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| B1 | geometry: engine on the cellgeom_v1 package vs legacy `cell3d` on `jgeom_cell1_v7.npz`, 2,000,000 uniform points in the enclosing box + 500,000 points within 50 nm of a membrane (seeded) | labels identical, membrane SDFs identical | labels ≥ 99.99 % identical; membrane SDF \|Δ\| < 0.01 nm (99.9th pct) | must |
| B2 | geometry measures: region volumes, membrane areas (tiles step) | identical `tiles.parquet` | each region volume / membrane area within 0.1 % | must |
| B3 | inventory pack table vs `ingredients_v7.tsv` | md5 d1a0099c… | per species copies \|log2 ratio\| ≤ 0.01 and comp_weights max \|Δ\| ≤ 1e-6 for all species; species set identical | must |
| B4 | intermediate inventory tables (`inventory_v0`, v4/v5/v6 outputs vs `provenance/snapshots`, change tables) | identical content | same as B3 per table; change-table rows identical except the documented `rule_id` fix | must |
| B5 | lamina rebuilt vs frozen `lamina_cell1.npz` | identical arrays | count ± 0.1 %; positions not required | report |
| B6 | `counts.npz` (expected copies per tile, library, region) | arrays identical | total requested instances ± 0.1 %; per (species, region) totals within max(3√n, 0.5 % n) | must |
| B7 | pre-pass index (`prepass_types/pos/quat/off`, mRNA bead index) | arrays identical | instance counts per tile identical | must |
| B8 | packed tiles: instance count per tile | identical for all 7,485 tiles | — (see B9) | report |
| B9 | placed copies per (species, region), whole cell | identical | within max(3√n_v7, 1 % n_v7), n_v7 = v7 placed copies, for every pair with n_v7 ≥ 10; ≥ 99 % of pairs with n_v7 < 10 within ±3 | must |
| B10 | placed mass per region (summary.json `regions.*.placed_pg`, `frac`) | identical | `frac` within ± 0.005 (0.5 percentage points) per region | must |
| B11 | cell totals: `instances`, `placed_pg`, `mass_frac` | identical (880,040,974 / 87.72 / 0.964) | instances ± 0.5 %; mass_frac 0.964 ± 0.003 | must |
| B12 | positions/orientations of instances | identical arrays (tile by tile) | not required | report |

## 3. Validation scores (same ported code on v7 tiles and on port tiles)

| # | Metric | Tolerance (port vs v7, same code) | Class |
|---|---|---|---|
| V1 | densities per volume region (mg/ml, packed + pre-pass) | ± 1 mg/ml per region; nuclear/cytoplasmic ratio ± 0.01 | must |
| V2 | crowding (excluded-volume fraction per region) | ± 0.005 absolute per region | report |
| V3 | loc_compare, placed shares per gene: agreement vs OpenCell, HPA consensus and HPA-Jurkat (unweighted and mass-weighted `agree`, `agree_tvd`, `tophit`) | each ± 0.005; per-gene placed compartment shares: ≥ 99 % of genes within 0.02 TVD | must |
| V4 | loc_compare HPA-Jurkat hold-out half (frozen list, seed 20261110) | ± 0.01 | must |
| V5 | chromatin metrics (frozen input): rep5 file checksum identical; metrics copied from the frozen v8 validation outputs (EM voxel r at σ 250 nm 0.68 etc.) | checksum identical; values copied verbatim | must |
| V6 | placement integrity on a seeded sample of ≥ 300 tiles (same tile ids for both): (a) clash rate = fraction of instance pairs whose sphere trees overlap by > 0.5 nm, per 10⁶ instances; (b) containment rate = fraction of instances whose centre label is not allowed for their region; (c) membrane-instance orientation sign errors | each rate within max(10 % relative, 1 per 10⁶ instances absolute) | must |
| V7 | ribosome nearest-neighbour distance distribution (polysome pre-pass + packed 80S) | median ± 0.3 nm | report |
| V8 | Tier-2 imaging and fluorescence panels | not ported in phase 1 (stated as not evaluated) | report |

## 4. Ported validation code vs stored v7 numbers (code check, separate from the verdict)

Run on the v7 tiles. Expected identical to the stored v7 files except where the audit documents a defect:
`loc_holdout_v7.json` / `loc_compare_v7.tsv` (expected identical within ± 0.001), `crowding_v7.json` (identical),
`densities_v7.json` (differs only through the v6 mRNA file: cytosol/nucleoplasm packed+pre-pass mass; differences
explained), `ribosome_nn_v7.json` (= v6 copy; not comparable). Deviations are reported, not hidden.

## 5. Reporting

`validate/compare_v7/` in the port build holds `criteria_results.tsv` (row id, metric, v7, port, Δ, tolerance,
level E/T/fail, class) and `verdict.json`. The comparison script is `wholecell/validate/compare_builds.py`; it reads
this file's thresholds from its own constants, which must match the tables above (reviewed at commit).

---
(Addenda after the run go below this line, dated.)
