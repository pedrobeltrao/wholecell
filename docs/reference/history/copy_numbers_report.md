# Copy-number step: shared Jurkat + yeast method, regression, anchor check (2026-11-16)

Code: `wholecell/inventory/copies.py` (step `copies`, proteome adapter `copies`), README `wholecell/inventory/COPIES.md`,
spec `docs/specs/inventory_schema.md` §4b, criteria `docs/copies_regression_criteria.md` (committed 7f321d6, before any
run), tests `tests/test_copies.py`, checks `scripts/copies_regression_check.py`, `scripts/copies_yeast_check.py`.

## 1. Design

One generic step between `rna` and `v0`; a cell's method is the config block `inventory.copies`:

1. **Annotation and universe**: UniProt readers keyed by gene (Jurkat) or by gene with ORF fallback (yeast); universe
   restricted to proteins with a source value or an RNA row (Jurkat) or the whole annotation (yeast).
2. **Typed sources**: `relative_ibaq | ppb | ppm | lfq | absolute_copies`. Each source has its own reader (MaxQuant
   protein groups, PaxDb, plain table), replicate aggregation (median / mean / geomean of detected replicates) and
   protein-group split rule (`family_prefix` n = 3 = the Jurkat histone/tubulin split, `first`, `equal_all`). LFQ is
   divided by mass; secondary sources are put on the first source's scale by the median ratio.
3. **Combination and priority**: `priority` (first source with a value) or `geometric_mean`. `set_by` records which
   source or rule set each protein.
4. **Imputation**: `none` or `rna_loglog` (the Jurkat fit log10 ppb ~ log10 TPM, capped at the 5th percentile, `low`
   calls x 0.25, `off` / no RNA -> 0).
5. **Absent lists with reasons**: gene lists, regex patterns and categories, applied before the mass scaling. Jurkat
   ZERO + 9p21 patterns; yeast NOT_IN_THIS_CELL (5 categories, `plugins/yeast/data/inventory_copies.yaml`).
6. **Per-protein scaling** (optional): `size_slope` (yeast round U2 G1 scaling, (34/42)^slope, Lanz 2024) or
   `factor_table`.
7. **Absolute scale**: total protein `fixed` (value or `inventory.total_protein_pg`), `anchors` (implied totals of the
   chosen anchors combined by geometric mean / median / mean / min / max, optional rounding), or `source_sum`
   (absolute copies keep their mass).
8. **Structure-count adjustments** (optional, after scaling): `family_target`, e.g. yeast round U1: ribosomal proteins
   scaled so that their median family = 14,250 per µm³ x 18.44 µm³ cytosol; the other proteins keep their copies.
9. **Anchor check** (diagnostic, always run): for each anchor (histones vs nucleosomes, nucleoporins vs pores,
   r-proteins vs ribosomes, ...) the family ratios copies / (stoich x count), their median, and the implied total
   protein = total / median ratio, for every count given; the spread and combined values across anchors.
10. **Copy ledger**: `copy_ledger.tsv`, one row per protein and rule (`source.*`, `impute.*`, `absent.*`, `scale.*`,
    `adjust.*`, `missing.zero`; global `total.*`), plus `set_by` / `rule_ids` columns in `copies.tsv`.
11. **Future sources**: `future_sources` declares SEC-MS abundance (cross-check) and surface proteomics (PM copies)
    for both cells; listed in `copy_sources.tsv` as "declared, not enabled"; enabling raises until implemented.

## 2. Jurkat regression (criteria C1-C7)

Run `{results_root}/jurkat/_dev_copies/regress-f0e7e04` (inventory stage only, no pack, CPU, PYTHONHASHSEED=0; code = git
archive of the commit; `copies_regression_verdict.json`). A first run at cd970c7 gave the same tables.

| id | result |
|---|---|
| C1 `v0/inventory_v0.tsv` md5 = port v0 (b40f1ba2...) | pass |
| C2 vs legacy `inventory_v0.tsv` | pass: same genes and order, `copies` identical for all 19,894 genes; only the 2,218 known order-only `compartments` cells differ (no float or real differences at rtol 1e-14) |
| C3 `copies.tsv` copies and labels = v0 | pass |
| C4 `pack_table.tsv` md5 = d1a0099c... (= v7 `ingredients_v7.tsv`) | pass (byte-identical) |
| C5 `compare_legacy.json` statuses = v7repro2 | pass (19 entries, no status difference, B3 all pass) |
| C6 generic Jurkat recipe = ported legacy adapter on synthetic tables | pass (`tests/test_copies.py`) |
| C7 provenance | pass: every protein has `set_by`; ledger 19,917 rows |

**Verdict: REPRODUCED.** Jurkat now runs `adapters.proteome: copies`; `jurkat.geiger_ibaq_v0` stays registered for the
equivalence test. Expected, non-criterion change: canonical provenance strings name the adapter `copies`.

## 3. Yeast (config `config/cells/yeast.yaml`, criteria Y1-Y4)

The yeast method was expressed in the same schema (PaxDb ppm per ORF, NOT_IN_THIS_CELL, G1 size scaling, 3.5 pg fixed,
ribosomes from the cryo-ET density, anchors) and run read-only on the yeast inputs in place (`{yeast_root}`: UniProt
export, PaxDb file, Lanz 2024 Table S3; no yeast code imported, nothing written there).

| id | result |
|---|---|
| Y1 same 6,067 proteins as build W `protein_features.parquet` | pass |
| Y2 copies_per_cell, max relative difference | 4.3e-16 (float order only); zeros match; pass |
| Y3 NOT_IN_THIS_CELL counts and ppm shares; ribosome target 262,722, factor 1.185, total 3.61 pg | identical to the build W log; pass |
| Y4 anchors.json nucleoporin median 1.3989 and histone ratios H3 0.709 / H4 1.046 / H2A 0.861 / H2B 1.478 | identical; pass |

**Verdict: REPRODUCED** (`{results_root}/yeast/_dev_copies/check-f0e7e04`). Variant (diagnostic): total protein from
anchors (core histones + nucleoporins, geometric mean, rounded to 0.5 pg) on build-W inputs gives **3.0 pg**, not 3.5:
the 3.5 pg choice was made in round R (before the G1 scaling, at 5 pg: histones 3.9, Nups 3.2); with the round-U inputs
the anchors imply histones 3.67, nucleoporins 2.50, r-proteins 2.66 pg (vs 2 x 10^5 ribosomes; after the ribosome
adjustment). Not adopted: the yeast model is theirs to change; reported for the joint paper.

## 4. Jurkat structural anchors (diagnostic only; nothing in the model changed)

v7 inventory at 100 pg (v0 copies, before assembly). Implied total = 100 pg / median family ratio.

| anchor | count used | families | median ratio (p10-p90) | implied total protein (p10-p90 of families) |
|---|---|---|---|---|
| core histones (2 per nucleosome, all variants per type) | 28.19 M nucleosomes (genome-mapped, chromatin v8) | 4 | 0.93 (0.74-1.06) | **108 pg** (95-135) |
| | 31 M (DNA, 2n) | | 0.85 | 118 pg |
| nucleoporins (7R5K chains per pore, scaffold Nups) | 300 pores (EM np_seg, under-detected) | 22 | 7.46 (3.9-17) | **13 pg** (6-25) |
| | 1,000 pores (model geometry) | | 2.24 (1.2-5.1) | 45 pg (20-84) |
| cytosolic r-proteins (80S families, paralogs pooled) | 1.5 M ribosomes (EM recall-corrected) | 77 | 1.40 (0.58-2.33) | **71 pg** (43-173) |
| | 0.756 M (EM raw, lower bound) / 1.7 M (upper) | | 2.79 / 1.24 | 36 / 81 pg |

Spread of the primary anchors: 8.0 x (13-108 pg); without the nucleoporins 1.5 x (71-108 pg; geometric mean 88 pg).

Interpretation:
- **Histones agree with 100 pg** (108 pg against the 28.19 M placed nucleosomes, 118 pg against 31 M from DNA). The
  per-type spread (H3 0.67 to H2A 1.10) reflects how protein groups of identical histones were split among genes; the
  median is robust. This replaces the earlier "histone ruler 130-250 pg" (whole-cell DNA mass, 2n / 4n) with a direct
  count of the nucleosomes the model holds.
- **Ribosomes say ~70-80 pg** with the EM recall correction (1.4-1.7 M): the inventory's median r-protein family is
  2.1 M copies, 1.4 x the EM-corrected count. The known residual tension (PLANNING_NOTES 2026-11-08: inventory sets the
  count, 2.1 M vs 1.5 M) is the same finding seen as a mass. Small r-proteins (RPL36, RPL41, RPL37) are undetected by MS
  and imputed at the floor; they do not move the median.
- **Nucleoporins do not constrain the mass**: they constrain the pore count. At 100 pg the Nups make 2.2 pores' worth
  per model pore (7.5 per EM-detected pore), i.e. about 2,200 pores (~8.7 per µm² of the 258 µm² NE) if all Nups sat
  in pores. Nups have off-pore pools (NUP98, NUP62, NUP88, NUP54 and NUP210 have the highest ratios) and np_seg under-
  detects pores, so this anchor is not a mass ruler until the pore count is measured on raw EM. Yeast shows the same
  direction (1.4 Nups per pore need at 3.5 pg).
- **100 pg sits between the two usable anchors** (ribosomes 71-81 pg, histones 108-118 pg), as yeast's 3.5 pg sat
  between its histone and nucleoporin values in round R. No change proposed; if a single rule were adopted
  (geometric mean of histones and ribosomes) it would give ~88 pg, inside the 80-130 pg bounds.

## 5. Next steps

1. Measure the cell-1 NPC count on raw EM (the nucleoporin anchor and the 1,000 model pores both depend on it).
2. Implement the two declared sources: SEC-MS summed intensity as an abundance cross-check (both cells, from
   `inventory/secms/profiles.tsv`), surface proteomics as PM copy numbers (an adjustment kind) when the data arrive.
3. Yeast: adopt `inventory.copies` when the yeast port starts (the config already reproduces build W); discuss with the
   yeast session whether the scale should follow the anchors on the round-U inputs (3.0 pg) or stay fixed at 3.5 pg.
4. Optional for Jurkat v8 (Pedro's decision): ribosome-vs-histone anchor rule as the scale (~88 pg), or keep 100 pg.
