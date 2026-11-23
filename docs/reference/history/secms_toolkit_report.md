# SEC-MS toolkit: diagnostic report (2026-11-16)

Code: `toolkits/secms/` (README there), inventory step `secms` + `wholecell/inventory/secms.py` (SecView,
`assembly_options`), Jurkat assembly adapter `plugins/jurkat/inventory/assembly.py`. Spec: `docs/specs/inventory_schema.md`
§4a. Scope (Pedro 2026-11-16): human and yeast only (species switch); calibration = curated native-mass anchors (v7
approach); Complex Portal calibrants are a diagnostic only.

Runs (all CPU, no packing; `{results_root}/jurkat/_dev_secms/`):

| run | what |
|---|---|
| `regress-6da56bb` | full Jurkat inventory (rna -> compare) at the final commit, preset `legacy_v7`; `secms_regression_verdict.json` |
| `regress-6da56bb_s8inv80` | v0 + complexes at 80 pg (explains the orphan table, criterion S8) |
| `dryrun-freepool-6da56bb`, `dryrun-freepool-max1000-6da56bb` | `complex_copies: measured_free_pool` (with / without a 1 MDa size gate) |
| `option-{homo_const,homodimer_const,hetero_nosec}` | each existing SEC-dependent rule switched off alone |
| scratch toolkit runs | Jurkat with `default` vs `legacy_v7`; HeLa (Heusel 2020, PXD010288) with `default` as format test + calibration check |

## 1. Regression vs v7 (criteria frozen before the first run: `docs/secms_regression_criteria.md`, md5 31918e75, commit 78fdf83)

**Verdict: REPRODUCED** (11/11 criteria).

| id | result |
|---|---|
| S1-S6 | `sec_assembly`, `sec_calibrants_peak`, `homo_oligomers_v0`, `afdb_homodimers_v0`, `afdb_heterodimers_v0`, `complexes_v0`: md5 identical to legacy |
| S7 | `inventory_v0_complexed` identical to the v7repro2 port (only the known 2,218 order-only `compartments` cells vs legacy) |
| S8 | the orphan `sec_assembly_thresholds.tsv` now has a generator: toolkit a_k (k = 2/3/4/6) equal to the orphan for all 3,040 proteins (max abs diff 0). Its `rule` and `mass_pg` columns came from the stale **80 pg** inventory: recomputed at 80 pg they match exactly (rows, order, values). At 100 pg 46 `rule` values differ by <= 0.049; no consumer reads these columns |
| S9 | `pack_table.tsv` md5 d1a0099c (= ingredients_v7) |
| S10 | `compare_legacy.json`: same status as v7repro2 for every table / ingredient entry |
| S11 | `complex_copies: stoichiometry_cap` in the config and in the build (`assembly/assembly_options.json`) |

Tests: `tests/test_secms.py` (6 tests: calibration recovers the synthetic column law; assembled fraction vs truth;
complexes, hu.MAP, pairs; species switch and presets; measured free pool on testcell; Jurkat plan with the option off).
`tests/test_inventory.py` passes with the new `secms` step.

## 2. Calibration quality

| dataset | method | n | Spearman | residual (fold) | LOO (fold) | app. MW f20 / f40 / f60 |
|---|---|---|---|---|---|---|
| Jurkat Fossati 2021 (60 fr) | anchors (used) | 25 | -0.89 | x1.85 | x2.03 | 885 / 92 / 9.5 kDa (void f8 = 3.4 MDa) |
| | Complex Portal (check) | 52-60 | -0.28 to -0.35 | x2.2 | x2.3 | too flat (slope -0.011 to -0.013 vs -0.049) |
| HeLa Heusel 2020 interphase (65 fr) | anchors (used) | 23 | -0.83 | x2.47 (QC limit x2.5) | x3.1 | 1,452 / 181 / 23 kDa |
| HeLa mitosis | anchors (used) | 23 | -0.86 | x1.88 | x2.29 | |
| HeLa | Complex Portal (check) | 129-132 | -0.57 / -0.61 | | | |
| HeLa | standards (authors' calibration, from their `monomer_fraction`) | 8 | -1.00 | | | 997 / 102 / 10 kDa |
| testcell (synthetic) | anchors | 14 | -1.00 | x1.03 | x1.06 | recovers log10 MW = 3.6 - 0.055 f (icpt 3.58-3.59, slope -0.054) |

Findings:
- **Anchors vs globular standards (HeLa, same column):** anchors put each fraction **x1.7 higher in MW** (median over
  standard points; x1.5 at f20, x2.3 at f60), so native proteins elute later than globular standards of the same mass.
  The k = 2 assembled fraction is **+0.18 higher on average** under the anchor calibration (median 0.98 vs 0.94;
  Spearman 0.86 between the two). The choice of calibration moves assembled fractions much more than any other parameter.
- **Anchor outliers (Jurkat):** PFN1 (x0.18: elutes at ~80 kDa, presumably profilin-actin), HSPA8 (x0.18, client
  complexes), HSP90AB1 (x0.26), YWHAZ (x0.38) elute earlier than their nominal native mass; PPIA (x1.9) and EEF2
  (x1.8) elute later. On HeLa, NPM1 is in the void (fraction 3). These anchors are partly complexed in cells. Pruning
  them is a parameter change for every cell and needs a DECISION (not done; v7 set kept).
- Complex Portal peak-level calibration remains unusable as the primary method in Jurkat (rho ~-0.3). It is better on
  HeLa (-0.6) but still below the QC gate (-0.7).

## 3. Assembled fractions: toolkit defaults vs legacy (Jurkat)

`default` and `legacy_v7` use the same calibration and the same integration, so where both score a protein the
values are **identical** (2,274 proteins, max diff 0). The differences are elsewhere:

| | legacy_v7 | default |
|---|---|---|
| proteins scored | 3,110 | 2,274 (`min_peptides: 2` drops 836 single-peptide proteins) |
| mass-weighted assembled fraction (Jurkat inventory weights), k = 2 / 3 / 4 / 6 | 0.750 / 0.660 / 0.578 / 0.479 | 0.745 / 0.652 / 0.569 / 0.470 |
| calibration +/- 1 SD (k = 2, mass-weighted) | 0.60-0.87 | 0.59-0.87 |
| QC gate, uncertainty columns, peaks, void fraction | none | yes (void: 7.4% of mass-weighted intensity at fractions <= 8) |

- Peak-level vs threshold estimator: Spearman 0.78 (median 2 peaks per protein).
- **Exact dimers sit on the k = 2 boundary.** On synthetic data, states >= 2.5x monomer are recovered within 0.08
  (70 proteins). For exact dimers, about half of the dimer intensity falls below 2x monomer because of peak width, so
  measured 0.40-0.51 vs true 0.70-0.90. This holds for legacy and default alike, and the AFDB homodimer share (a2) is
  therefore biased low.
- Jurkat vs HeLa (default, both anchor-calibrated): Spearman 0.53 over 1,975 shared proteins.
- Uncertainty: the calibration SD dominates (lo-hi spread about 0.27 mass-weighted). Bootstrap from the HeLa SEM gives
  a median SD of 0.003. There is no replicate SD for Fossati (single run).

## 4. What each option moves (Jurkat dry runs, no repack; reference = v7 rules)

Reference: Complex Portal complexes 36.2 pg (incl. curated nucleosomes); complexed share of 101.1 pg protein = 52.5%
(complexes + homo-oligomers 6.13 + AFDB homodimers 6.21 + heterodimers 4.54 pg).

| option switched | proteins / species changed | mass moved | complexed share |
|---|---|---|---|
| homo-oligomer share: SEC a_k -> constant 0.5 | 691 species | 3.4 pg (homo-oligomers 6.13 -> 5.34 pg) | 52.5 -> 51.8% |
| AFDB homodimer share: SEC a2 -> constant 0.3 | 684 species | 2.7 pg (homodimers 6.21 -> 2.35 pg; heterodimers +0.85) | -> 49.6% |
| heterodimer SEC evidence off (STRING >= 700 only) | 129 species | 0.16 pg (1,402 -> 1,375 pairs) | -> 52.4% |
| **complex copies: measured free pool instead of the stoichiometry cap** | 356 complexes changed (110 up, 246 down; 535 measured); 733 proteins; 1,189 species | Complex Portal mass 36.2 -> 26.0 pg; 7.0 pg protein -> free pools (0.19 pg the other way); 12.6 pg species mass moved; downstream homo/dimers +0.6 pg | **52.5 -> 43.0%** |
| same, complexes > 1 MDa keep the cap (`max_expected_kda: 1000`) | 307 complexes; 566 proteins; 980 species | 36.2 -> 28.1 pg; 5.3 pg -> free | -> 45.0% |

Largest measured-free-pool changes: 40S subunit 2.37 M -> 0.22 M copies (-2.1 pg; the 100,000 g pre-SEC spin pellets
ribosomes, so SEC sees mostly the free pool), major spliceosome pre-B 294 k -> 0 (-1.3 pg), 26S proteasome 464 k -> 163
(-1.0 pg; dissociates into 20S/19S, which fall below the 0.5 x 26S window), CCT 954 k -> 567 k, LDH A1B3 halved,
DNA-PK halved, Arp2/3 781 k -> 18 k, COPI 295 k -> 41.

Assessment: as implemented, the option mostly measures what the Fossati sample preparation removed: ribosomes,
polysomes and chromatin are pelleted, and large machines dissociate under mild lysis. It is not ready for Jurkat
(kept OFF, as instructed). Before use it would need (a) a size or pelleting gate, (b) subcomplex-aware windows
(20S/19S count as assembled), and (c) a recovery correction (SEC input vs whole-cell copies). New datasets with
whole-lysate input and replicates would make it more credible.

## 5. Parameters a new dataset needs (dataset block of `params.secms_toolkit`)

Required: `source_id`, `format`, `path`, `n_fractions`, `condition`; for OpenSWATH a `fraction_map` (sample sheet
-> fraction, condition and replicate columns); for matrices a `value_regex` (+ `sem_regex`). Strongly recommended:
`void_fraction` (column void), `standards` (column standard fractions and kDa, for the cross-check), the replicate
structure, and the lysis / clarification conditions (they decide what the free pool means). Per species, not per
dataset: the annotation TSV (UniProt), the Complex Portal export, pair candidates. **Nothing in `presets/default.yaml`
should change per cell.**

## 6. Checklist for the coming Jurkat and yeast datasets

1. Freeze the raw files as `inv_*` frozen inputs (md5, origin, licence), including the sample sheet and column
   standards. Record the lysis buffer, clarification speed and column.
2. Add a dataset entry to `params.secms_toolkit.datasets` with `preset: default` (not `legacy_v7`). For Jurkat this is a
   new config (e.g. `jurkat_v8.yaml extends jurkat.yaml`), so the v7 reproduction config stays intact.
3. Run `python -m toolkits.secms run` stand-alone first, then read `diagnostics.json`: anchor calibration QC
   (rho <= -0.7, residual <= x2.5), LOO error, the standards check (expect about x1.7 vs anchors), the Complex Portal check, void.
4. **Yeast:** review `data/anchors_native_mass_yeast.tsv` (26 anchors, curated 2026-11-16, not yet tested on data)
   with G. Magnusson. Set the annotation's gene column (UniProt primary names) and use Complex Portal 559292 (yeast
   pipeline copy: `yeast_cell_dev_bundle/pipeline/data/complexportal_559292.tsv`, read-only; re-download for our tree).
   Yeast assembly state is currently stoichiometric with no SEC (median subunit pool with a 2x limiting cap, whole
   residual pool to homo-oligomers, heterodimers 0.5 or fully bound). Map these onto `assembly_options`:
   `homo_oligomer_share.source: constant` with value 1, and so on.
5. Compare with the Fossati / HeLa tables: assembled-fraction distribution, Spearman per shared protein, mass-weighted
   k = 2..6.
6. Decide (DECISION in PLANNING_NOTES) whether to prune the anchor outliers (PFN1, HSPA8, HSP90AB1, YWHAZ) for every
   cell, and whether a calibration from column standards should be the cross-check reported in papers.
7. Only after the above, consider `complex_copies: measured_free_pool` with the gates of §4 (size, subcomplex,
   recovery); dry-run it with `scripts/secms_freepool_dryrun.py`, never directly in a pack build.
8. Re-run `tests/test_secms.py`, and for Jurkat `scripts/secms_regression_check.py` on the v7 config (must stay REPRODUCED).

## 7. Notes

- HeLa: `monomer_fraction` in the Heusel table is the expected monomer fraction under the authors' standard
  calibration, not an assembled fraction (it gave the standards check above). No assembled-vs-monomer labels are local,
  so there is no per-protein ground truth for HeLa.
- No yeast SEC-MS data were downloaded (not trivially available offline). The yeast path is tested only through the
  species switch and the anchor file format.
- hu.MAP 3.0 is not local; the reader is tested on a synthetic file only.
