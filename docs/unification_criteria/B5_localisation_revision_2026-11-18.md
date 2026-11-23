# B5 localisation common defaults: revision 2026-11-18

**Status.** Design revision approved by Pedro on 2026-11-18, made before any B5 build.

**What stays fixed.** The frozen criteria `docs/unification_criteria/B5_localisation.md` (md5
`49a43723bfafb90ef5092b96f42bbf20`) are unchanged. That covers the gate (L1 C1 ≥ reference − 0.005 per reference set,
L2-L7) and the §5 report.

**What changes.** This file replaces §1 steps 5-7 and the §2 R1 target of the common defaults with the rules below.
Everything else in §1-§2 stands.

**Code.** Branch `b5-revise`, commit `7b9880f`:
- `wholecell/inventory/localisation/{combiner.py, presets/common.yaml, __init__.py}`
- `plugins/jurkat/inventory/rules_v6.py`
- the cell data in `config/cells/{jurkat,yeast}.yaml`
- `plugins/{jurkat,yeast}/data/localisation/chromatin_pool.tsv`

**Dry runs.**
- Jurkat: `wholecell_builds/jurkat/_dev_unify/b5/dryrun_7b9880f/`
- Yeast: `wholecell_builds/yeast/_dev_unify/b5/dryrun_7b9880f/`

## 0. Why

Pedro did not accept the previous dry run, `dryrun_e094024`:

- **N/C model mass.** The Jurkat N/C model moved 11.3 pg and `new_from_imaging` moved 2.6 pg.
- **C1 intended agree.** It failed on 3 of the 4 reference sets:

  | Reference set | Common − legacy |
  |---|---|
  | OpenCell | +0.0095 |
  | HPA consensus | −0.0061 |
  | HPA-JURKAT | −0.0156 |
  | HPA-JURKAT hold-out | −0.0093 |

- **TKT and PGAM1 became 100 % nuclear.** HPA nucleoplasm was the only imaging call. Itzhak 2016 fractionation, which
  says cytosolic, no longer counted.
- **Histones leaked to the cytosol.** H2B was about 1/3 cytosol, because an HPA additional "Cytosol" term raised
  category Nc (r = 3).
- **The nucleus was denser than the cytosol.** Intended nucleoplasm was 131 mg/ml against 115 mg/ml for the cytosol.

## 1. Revised defaults (preset `common`, one set for every cell; a cell differs only in data)

Step order (combiner, per protein):

`drop_unreliable` → key selection (cell data) → `resolve_imaging` → `restrict_membrane` → `split` →
**`imaging_background`** → `combine` → `new_from_imaging` → `nc_model` → `round`

### R-a. Imaging nuclear background term b (new step `imaging_background`)

**Model.** observed nuclear share = b + (1 − b)·true, so true = clip((obs − b)/(1 − b), 0, 1).

**What it acts on.**
- The input is an imaging source's split distribution, after score-proportional split and before `combine`.
- obs is the share of the source's nucleus keys.
- The nucleus keys are rescaled to `true` and the source's other keys to 1 − true.
- A distribution with obs = 0 or obs = 1 is unchanged.

**Which sources.** The sources listed in the cell data only. Jurkat lists `hpa`, the antibody immunofluorescence.
- It is never applied to live-GFP imaging (OpenCell, yeast CYCLoPs). There b = 0 by construction (estimate.json
  `b_OpenCell` = 0).
- Yeast declares no background term.

**b used: 0.1931.** Source: `data/derived/v8/hpa_background/estimate.json` (md5 `895e677fb5b4825759c3a71e35cbbbb2`;
`hpa_nuclear_background.py` @ 8e2452a, branch `hpa-background`), `decision.b_HPA`.
- The estimate falls in the **intermediate branch**, not the real-pool branch (b = 0).
- OpenCell live HEK293T non-nuclear genes have a median nuclear fraction of 0.139 in Jurkat geometry (CI 0.116-0.174).
  The HPA-Jurkat floor is 0.305.
- The fitted b is 0.1931 (CI 0.1154-0.2675).
- As the plan requires for this branch, the extremes b = 0 and b = 0.3051 are reported as a sensitivity (§3).

**Provenance.** The value sits in `inventory.params.localisation.common.imaging_background` of `config/cells/jurkat.yaml`
(b, provenance string). It is logged by the step.

**Stated approximation.** b was estimated on image-derived nuclear fractions. HPA enters the combiner as categorical
shares (main 1 / additional 0.3, split score-proportionally), and b is applied to those shares.

**Flags.** A source corrected by this step is marked background-corrected. When its nuclear terms are removed entirely
(true = 0), they no longer name a nucleus side in the N/C model.

### R-b. Background-corrected imaging is blended as in v7 (`new_from_imaging.background_sources: all`)

**What it does.** For an imaging source corrected by `imaging_background`, every key counts as supported. This includes
the HPA additional terms and nucleus or cytosol named by HPA alone.

**Why.** The two `new_from_imaging` guards were defences against imaging noise:
- min_score 0.5, the X5 CYCLOPS_NEW_MIN for continuous CYCLoPs LOC-scores;
- "never a nucleus / cytosol key from imaging alone", against the HPA nuclear floor.

For HPA that noise is now removed explicitly by R-a. In the e094024 dry run the guards were the largest single cause of
the HPA-set losses:
- an HPA "Nucleoplasm" main call on a UniProt-endolysosomal protein was dropped entirely (WASHC4, MON2, ANGEL1);
- HPA-JURKAT additional calls were dropped (IKZF1, CD2, MYO1G, GPR174).

Sources without a background term (yeast CYCLoPs) keep both guards unchanged.

**Effect in the dry run.**
- HPA-JURKAT agree +0.010 over keeping only HPA main pool keys (scratch variant, b = 0.27).
- HPA-JURKAT hold-out +0.007 over the same variant.

### R-c. N/C model: Itzhak 2016 fractionation as a prior / tie-breaker (`nc_model.prior`; Jurkat data)

**Inputs.**
- The fractionation share fn_f = clip((obs − 0.05)/0.95, 0, 1). 0.05 is the v7 nuclear base.
- Its concentration ratio r_f = (fn_f/Vn)/((1 − fn_f)/Vc), clipped to [0.01, 100].

**Combination, in log space:**
- With an imaging ratio: log r = (log r_img + w·log r_f)/(1 + w), with **w = 1**.
- Without imaging calibration: the equal-concentration default r = 1 takes the imaging term's place, so
  log r = w·log r_f/(1 + w). The fractionation alone never sets the split; this is the "prior" reading.
- Category NC: unchanged. The v6 R1 tie-breaker already replaces r by the clipped Itzhak ratio there, so the prior is
  not applied a second time.

**Why w = 1.** Imaging (HPA consensus over cell lines) and fractionation (HeLa) measure the same quantity with
different, comparable errors, and neither is a Jurkat measurement.

**Why the clip is [0.01, 100].** It is wider than the category range [0.02, 10] because fractionation is quantitative.
For example, HMGN2 (Itzhak 0.993 nuclear) stays 0.955 nuclear instead of being capped at r = 10.

**Evidence that it is the right direction.**
- OpenCell (live, no antibody, independent of the model inputs): on the 294 OpenCell genes that Itzhak calls cytosolic
  (fn < 0.15), agreement rises by +0.024 against legacy.
- HPA consensus: on the 1,504 such HPA genes it falls by 0.036 (b = 0.27 variant, scratch).
- The HPA sets are not independent of HPA as a model input, so they reward HPA dominance.

### R-d. Bounded mixtures (`nc_model.bound: {max_share: 0.9, min_sources: 2, calibrations: [category]}`)

**The bound.** When the split comes from a categorical imaging calibration (HPA category → r), the nuclear share of the
soluble pool is kept in [1 − m, m] with **m = 0.9**. The exception is when ≥ 2 independent sources are exclusive on that
side.

**Counting exclusive sources.**
- An annotation source counts if all its pool keys are on one side.
- The imaging source counts if its category is N or C.
- The fractionation counts if fn_f ≥ m or fn_f ≤ 1 − m.

**Effect on proteins with one pool side.** When the bound applies against that side, the other side is added. A single
imaging category therefore never leaves a protein at 100 % of the nucleus or of the cytoplasm. Examples: an HPA-only N
protein gets ≥ 10 % cytosol, and an HPA-only C protein gets 10 % nucleus.

**Why m = 0.9.**
- In this geometry the strongest calibrated category (N, r = 10) gives fn = 0.87. m = 0.9 sits just above it, so the
  bound does not alter calibrated N splits. It only clamps C (r = 0.02 → fn ≥ 0.1, i.e. r ≥ 0.17) and stops the
  "unchanged single side" 100 %.
- A single image cannot establish exclusivity:
  - nuclear proteins show a ~0.77 image nuclear fraction in HPA-Jurkat;
  - OpenCell cytosolic proteins show 0.14-0.21 in the nucleus (estimate.json).

**Calibrations.** The bound applies only to categorical calibrations. Yeast's CYCLoPs ratio exp(5(N − C)) is a
continuous fit of LOC-scores (NUC_K on SWAT), not a category. Bounding it lowered the yeast intended image N/C test
Spearman from 0.833 to 0.767, below the §3d 0.80 (scratch run). Yeast must keep §3d.

**Annotation-only and organelle proteins.** These are not bounded: the rule is about imaging categories. A protein whose
only evidence is a single HPA organelle term stays 100 % in that organelle. It is not "moved" by imaging; it is placed by
its only evidence.

### R-e. Chromatin-bound pools (`nc_model.chromatin`, cell data regexes)

**Rule.** A gene in the chromatin-bound pool that has a nucleus key gets its whole soluble pool (nucleus + cytoplasm
keys) in the nucleus, labelled `nc:chromatin`. No prior, bound or category applies to it. Without a nucleus key it is
unchanged. v6 R1 (common) skips these genes, so no rule spreads them back to the cytosol.

**Jurkat** (`plugins/jurkat/data/localisation/chromatin_pool.tsv`):
1. the histone families the inventory consumes into the pre-placed nucleosomes (= `histone_families.tsv`: H2A incl.
   H2A.Z/H2A.X/H2A.J/macroH2A, H2B, H3 incl. H3.3/H3.5, H4);
2. the linker histones `^H1-\d+$`;
3. the other histone-fold variants of GO CC nucleosome (GO:0000786) restricted to histone genes (CENPA, H3-4, H3-7,
   H3Y1/2, H2AB*, H2AL3, H2AP, H2BW*, H2BK1, H2BN1, H2BC*L, H4-16).

The GO term alone was not used. It also lists ACTB, RUVBL1/2, SPHK2 and the TIP60/SRCAP remodellers, which are not
chromatin-bound pools. HMGN/HMGB are not in the rule (dynamic binders). With the prior, HMGN2 is 0.955 nuclear
(legacy 0.993).

**Yeast** (`plugins/yeast/data/localisation/chromatin_pool.tsv`, on `gene_name`): the histone classes of the assembly
data (HTA1/2, HTB1/2, HHT1/2, HHF1/2, HTZ1, CSE4) plus the linker histone HHO1.

In the Jurkat run, 58 genes with a nucleus key take the chromatin rule. In yeast it changes 1 protein; the other yeast
histones were already 100 % nucleoplasm.

### R-f. v6 R1 (common) uses the same target

**The change.** The re-expressed R1 (criteria §2) computes its target with the same function as the v0 N/C model,
`combiner.nc_target`: HPA category or the hold-out non-JURKAT consensus, then the Itzhak prior, then the bound. Its
per-gene prior and exclusivity inputs come from the new v0 record `inventory/localisation/nc_model.tsv`. Chromatin genes
are skipped. The gates (2-fold, one-sided N/C, capacity) are unchanged.

**Why.** Without this, R1 would re-raise HPA-N proteins to r = 10 and undo R-c and R-d. TKT and PGAM1 are 0.87 nuclear in
the legacy final table because of legacy R1.

### Unchanged

- Score split, precedence classes with the 0.5 : 0.5 blend, add_min 0.3 / new_min 0.5, round 4 / 1e-4.
- v5.P, and all `legacy` rule steps.

### Dry-run report

The dry-run report adds ablation stages g (R-a + R-b), h (R-c), i (R-d) and j (R-e) between e and f, plus the report
items of §2 (`report.json` key `revision_2026_11_18`, `top_changes.tsv`).

## 2. Dry runs (inventory only, CPU 3 cores; commit 7b9880f; b = 0.1931)

Legacy = `legacy_v7` / `legacy_x5` from the same code and inputs (identical to the e094024 legacy runs). "Final" = the
intended species table after all rule sets (pack table).

### Watched proteins (Jurkat; nucleus share, rest cytosol unless stated)

| Protein | v0 legacy | v0 e094024 | v0 revised | final legacy | final e094024 | final revised |
|---|---|---|---|---|---|---|
| TKT | 0.022 | 1.0 | **0.280** | 0.869 | 1.0 | **0.280** |
| PGAM1 | 0 | 1.0 | **0.174** | 0.869 | 1.0 | **0.174** |
| H2B (H2BC12/17/9/4/18, all alike) | 0.885 (cyt 0.115) | 0.666 (cyt 0.334) | **1.0** | consumed into nucleosomes | – | – |
| H1-4 | 0.885 (cyt 0.115) | 0.666 (cyt 0.334) | **1.0** | 0.737 + speckle 0.148, cyt 0.115 | 0.558 + speckle 0.108, cyt 0.334 | **0.826 + speckle 0.174, cyt 0** |
| H1-0 | 0.812 (golgi/cytoskeleton 0.094 each) | 0.666 (cyt 0.334) | **0.884 (golgi 0.116, cyt 0)** | as v0 | 0.666 (cyt 0.334) | as v0 |
| H1-2, H2AZ1, H3-3A, MACROH2A1 | 1.0 | 1.0 | 1.0 | 0.764 + nucleolus 0.236 (H1-2) | 0.776 + nucleolus 0.224 (H1-2) | 0.770 + nucleolus 0.230 (H1-2) |
| HMGN2 (not protected) | 0.993 | 0.869 (cyt 0.131) | 0.955 (cyt 0.045) | 0.993 | 0.869 | 0.955 |
| HMGB1 (not protected) | 0.26 (cyt 0.41, 4 organelles 0.083) | 0.58 (cyt 0.087) | 0.449 (cyt 0.218) | 0.867 (cyt 0.072) | 0.580 (cyt 0.087) | 0.673 (cyt 0.129) |

### Intended soluble protein density

Pack table, nucleus + nucleolus + speckle over V nucleoplasm, against cytosol + cytoskeleton + unassigned over
V cytosol (mg/ml):

| | Nucleoplasm | Cytosol |
|---|---|---|
| Legacy | 122.7 | 116.3 |
| e094024 | 128.7 (131 at the P capacity check) | 112.1 (115) |
| **Revised** | **117.2** | **118.7** |

The nucleus is no longer denser than the cytosol.

### Mass moved per rule (Jurkat, criteria §5; ablation on the common run's evidence)

| Stage | Proteins | Top changes | pg |
|---|---|---|---|
| d new_from_imaging (frozen form) | 3,055 | 259 | 2.60 (review) |
| e nc_model (frozen form) | 4,908 | 1,215 | 11.30 (review) |
| g background b + HPA blended as v7 | 3,430 | 319 | 2.18 (review) |
| h Itzhak prior | 1,456 | 627 | 4.20 (review) |
| i bound 0.9 | 1,396 | 139 | 0.61 |
| j chromatin | 17 | 1 | 0.42 |
| f round | 6,804 | 1 | 0.03 |

Stages d, e, g and h are steps of one chain. Their net is what matters: **combiner total legacy → common
6.37 pg** (8,170 proteins, 735 top changes; e094024: 11.96 pg, 1,306).

Rule sets under common:
- v6.R1_nc 0.44 pg (e094024 0.07; legacy 6.06, review);
- v5.P 0.017 pg (legacy 1.08);
- every other rule ≤ 1.40 pg (v5.C_cortex), as before.

**Pack table total legacy → common: 5.86 pg** (8,414 species, 1,000 top changes; e094024 6.43 pg).

### Top 20 movers (v0, legacy → revised, pg moved)

ACTB 0.132, TKT 0.096, HMGB1 0.095, GAPDH 0.095, HSPA8 0.087, TRIM28 0.079, RAN 0.077, EEF2 0.073, RBBP7 0.071,
HSP90AB1 0.068, MCM3 0.065, H2BC18 0.063, H2BC9 0.063, HNRNPH1 0.062, MCM4 0.062, MCM6 0.062, DHX15 0.056, MCM5 0.053,
PCNA 0.052, NPM1 0.049.

Most are the legacy Itzhak re-split replaced by imaging × prior. Examples:
- TRIM28, MCM3-6, RBBP7 and PCNA go from about 0.4 nuclear to 1.0 (UniProt + HPA nuclear-only);
- GAPDH and EEF2 lose most of the legacy nucleus (C category, floor 0.1 of the pool).

### Proteins changing top compartment (v0, legacy → revised)

735 proteins, 6.39 pg (`top_changes.tsv`):

| Transition | Proteins |
|---|---|
| cytosol → nucleus | 230 |
| nucleus → cytosol | 124 |
| cytoskeleton → nucleus | 93 |
| cytosol → plasma membrane | 81 |

### C1 intended agree, unweighted (gate L1: ≥ legacy − 0.005)

| Set | Legacy | e094024 | Revised | Δ vs legacy | L1 |
|---|---|---|---|---|---|
| OpenCell | 0.8318 | 0.8413 | 0.8373 | +0.0055 | **pass** |
| HPA consensus | 0.8324 | 0.8263 | 0.8243 | −0.0081 | **fail** |
| HPA-JURKAT | 0.7487 | 0.7331 | 0.7414 | −0.0073 | **fail** |
| HPA-JURKAT hold-out | 0.7097 | 0.7004 | 0.7061 | −0.0036 | **pass** |

**L2, the other intended C1 metrics (± 0.01).** Outside the band:
- OpenCell tophit −0.014;
- HPA consensus agree_tvd −0.030 / −0.020 (unweighted / mass-weighted);
- HPA consensus tophit −0.028 / −0.020;
- HPA-JURKAT agree mass-weighted −0.043 (e094024 −0.045);
- HPA-JURKAT agree_tvd mass-weighted −0.016.

OpenCell agree mass-weighted is +0.024.

### Image N/C, C3i intended (Spearman with the HPA-Jurkat Tier-2 nc_hpa)

| | Legacy | e094024 | Revised |
|---|---|---|---|
| Hold-out (n 36; L3 gate) | 0.558 | 0.632 | 0.606 |
| All genes | 0.705 | 0.741 | 0.697 |

The hold-out gate (≥ legacy − 0.01) passes.

### Yeast (§3d benchmark L4; legacy_x5 → revised)

| Metric | Legacy | Revised | §3d threshold |
|---|---|---|---|
| SWAT_native intended tophit | 0.6671 | 0.6647 | ≥ 0.655 |
| hyperLOPIT intended tophit | 0.8462 | 0.8483 | ≥ 0.838 |
| Intended image N/C test Spearman | 0.811 | 0.833 | ≥ 0.80 |

All three pass, unchanged from e094024.

**Other yeast numbers.**
- SWAT_NOP1 0.7429 (legacy 0.7402); CYCLoPs 0.8153 (report only).
- Mass moved: combiner total 0.155 pg (4.4 % of the cell, review: stage c precedence classes 0.140 pg, as e094024). The
  revision stages move 0.0018 pg (chromatin, 1 protein; no background or prior for yeast; the bound does not apply to
  the CYCLoPs fit).
- Yeast histones stay 100 % nucleoplasm.

## 3. Sensitivity to b (estimate branch "intermediate": report both extremes)

Same code (7b9880f) with only b changed. Scratch runs: inventory + loc_compare, intended level.

| b | OpenCell | HPA consensus | HPA-JURKAT | HPA-JURKAT hold-out | C3i hold-out / all | Nucleoplasm / cytosol mg/ml |
|---|---|---|---|---|---|---|
| 0 (real-pool extreme) | 0.8378 | 0.8244 | 0.7429 | 0.7061 | **0.487** / 0.644 | 118.0 / 119.1 |
| **0.1931 (used)** | 0.8373 | 0.8243 | 0.7414 | 0.7061 | 0.606 / 0.697 | 117.2 / 118.7 |
| 0.3051 (full HPA floor) | 0.8372 | 0.8237 | 0.7409 | 0.7057 | 0.618 / 0.697 | 116.7 / 118.6 |
| legacy (gate − 0.005 / − 0.01) | 0.8318 | 0.8324 | 0.7487 | 0.7097 | 0.558 / 0.705 | 122.7 / 116.3 |

**What b does not change.** C1 moves by at most 0.002 across the range, and the L1 outcome is the same at every b.

**What b does change.** b decides C3i:
- with b = 0, the C3i hold-out (0.487) would fail L3 (≥ 0.548);
- at b ≥ 0.19 it passes.

This is the image-level evidence that the HPA nuclear floor is at least partly background.

## 4. Verdict against the frozen gate (dry-run level; placed metrics need the build)

**L1 (C1 intended).**
- Passes: OpenCell, HPA-JURKAT hold-out.
- Fails: HPA consensus and HPA-JURKAT.

**Other checks.**
- L3 C3i hold-out passes (it would fail at b = 0).
- L4 yeast passes.
- L5 needs review for combiner stages d/e/g/h. The chain net is 6.37 pg against the 2 pg per-rule limit.

**Why the two HPA sets fail.** It is by construction of R-c and R-d: legacy v6 R1 set HPA-N proteins to r = 10
(0.87 nuclear) regardless of fractionation, and the HPA sets reward that. The OpenCell improvement is kept but smaller:
+0.0055 against +0.0095 in e094024.

A B5 build with these defaults would be **NOT ACCEPTED** on L1 unless Pedro decides otherwise. The options:
- accept the HPA-set deviations as explained by the prior (an accepted-deviation record, like B6);
- lower the prior weight w;
- keep the frozen outcome.

## 5. Remaining concerns

1. **b is applied to categorical HPA shares.** It was estimated on image nuclear fractions. Per-gene imaging fractions
   (R8a) would be the natural input.
2. **R-b turns `new_from_imaging` off for HPA.** That step is effectively a no-op for Jurkat, a deviation from the
   frozen §1 step 6 for that cell. The step still runs for yeast.
3. **The bound is limited to categorical calibrations,** to keep yeast §3d. The CYCLoPs fit is unbounded.
4. **HMGN/HMGB and other chromatin binders are not protected.** The rule is histone-fold + linker histones only.
5. **Mass moved stays above 2 pg** in four ablation stages (review, L5).
