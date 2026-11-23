# B5 localisation common defaults: rules, acceptance criteria and reports (frozen before any B5 code, dry run or build)

Step B5 of the unification plan (`docs/reference/history/unification_plan.md` §3a-§3d, §5, §5b B5): "Localisation
common defaults (no pack-dependent rules; score-proportional split; N/C model with the current calibration)". Plan §5b
acceptance: "C1 ≥ reference − 0.005 per reference; C3 ≥ reference − 0.01; mass moved by any single rule ≤ 2 pg without
review". Order: B6 → B3 → B4 → B1 → B2 → **B5**; the B5 builds run after the B2 builds in the sequential queue
(`wholecell_builds/_queue/`). Written 2026-11-18 before any B5 code, dry run or build existed; the only numbers quoted
are those of existing builds (named). The md5 of this file goes to PLANNING_NOTES (coordinator). Later findings go to a
dated addendum, never into this file.

Binding decisions: A5 (presets `legacy_v7` / `legacy_x5` reproduce v7 / X5 exactly and stay for the reproductions);
DECISION 2026-11-17 "the `localisation_sources` term_map strings are cleaned up in B5"; DECISION 2026-11-16 (3) the HPA
30 % nuclear floor becomes a background term estimated from OpenCell live images (R8c) — **not part of B5** (not done
yet); plan §3c imaging-fitted fractions (R8a) — **not part of B5**.

## 1. The common defaults (preset `common`, one file for every cell)

`wholecell/inventory/localisation/presets/common.yaml`. A cell differs only in data: its sources (adapters + term maps),
its key vocabulary and lists, its key-selection data (§1 step 2) and its N/C calibration (§1 step 7), all under
`inventory.params.localisation` of the cell config. Steps, in this order (combiner, per protein):

1. **drop_unreliable** (unchanged): a source whose rows carry reliability 0 (HPA "Uncertain") is ignored for the gene.
2. **Key selection (cell data, no weights).** Steps that decide WHICH keys an annotation names, run only when the cell
   declares them: yeast keeps the X5 ordered-annotation selection (Huh 2003 precedence inside the adapter, spindle pole
   → nuclear envelope for the SPB core, generic mitochondrion refined from GO, GO CC fallback with GO membrane terms for
   membrane proteins, LOC_FIX override row, Vögtle 2017 measured sub-compartment / Morgenstern 2017 high-confidence
   proteome + veto). The selected key list becomes ONE annotation source (class of the source it came from: curated =
   UniProt, high-throughput = GO) with every key scored 1; its order carries no weight. Jurkat declares none.
3. **restrict_membrane** (common): a transmembrane protein keeps, per source, only the cell's membrane keys; an
   annotation source with no membrane key keeps all its keys, an imaging source with none is dropped (= the v7 rule).
   Yeast gets a membrane-key list (data; new).
4. **split = score**: within a source, shares proportional to the scores; equal when scores tie (all annotation terms
   score 1 → equal shares; HPA main 1 / additional 0.3; CYCLoPs LOC-scores). Replaces the yeast rank rule (1, ½, ¼) and
   the yeast MITO_SECONDARY 0.1 down-weighting (both dropped).
5. **combine = precedence classes**, one rule for all cells:
   - same-cell imaging (`imaging_cell`), when present, is used alone;
   - otherwise annotation = the first class present of `annotation_curated` > `annotation_ht`; imaging = `imaging_other`;
   - when both exist they are blended **0.5 : 0.5** at the resolution of the imaging (the cell's imaging groups: data;
     yeast X5 groups nucleocytoplasm / ER / mitochondrion, `periphery` resolved to the plasma membrane for membrane or
     PM-annotated proteins and to cytosol otherwise, no imaging for cell-wall proteins; Jurkat: every key its own
     group). Within a group the annotation's shares decide; a group named only by imaging uses the imaging shares;
   - `fractionation` never enters combine (N/C tie-breaker only, step 7; key selection for yeast, step 2);
   - no evidence: the cell's default key (Jurkat `unassigned`, yeast `cytosol`).
   For Jurkat this is the v7 50:50 HPA/UniProt blend. For yeast imaging no longer sets the group shares alone (X5) but
   carries half of them.
6. **new_from_imaging**: a group named only by imaging is kept only if its best imaging score ≥ **0.5** (`new_min`, the
   X5 CYCLOPS_NEW_MIN) and it is not a nucleus / cytoplasm pool key (these come from annotation or the N/C model, never
   from imaging alone); yeast keeps the X5 data condition "mitochondrion only for proteins of the high-confidence
   mitochondrial proteome". Then renormalised.
7. **nc_model** (common; replaces the Jurkat v0 Itzhak re-split and the yeast X5 NUC_K / volume pooling): for a soluble
   protein with a pool key (cell data: nucleus key; cytoplasm keys), the nucleus : cytoplasm concentration ratio is
   - r from the cell's imaging calibration ("the current calibration"):
     - Jurkat: HPA category (main / additional terms → N, Nc, NC, Cn, C exactly as v6 R1, organelle-main genes skipped)
       → r = 10 / 3 / 1 / 0.3 / 0.02 (`rules_v6.yaml` r_by_category, fit-half Tier-2 calibration); category NC uses the
       Itzhak 2016 tie-breaker (background 0.05, r clipped to [1/3, 3]); the frozen HPA-JURKAT hold-out genes
       (`inv_holdout_genes`) get no imaging calibration here (kept clean; v6 R1 handles them with its non-JURKAT consensus);
     - yeast: r = exp(k (N − C)) of the CYCLoPs Nucleus / Cytoplasm LOC-scores, k = 5 (SWAT dev-half fit);
   - else r = 1 (equal concentration) when the protein has both a nucleus and a cytoplasm key; unchanged when it has
     only one side;
   - shares fn = r·Vn / (r·Vn + Vc), Vn = `nucleoplasm`, Vc = `cytosol` of **this build's** `geometry/measures.json`
     (region_volumes_mc); the cytoplasm share stays in proportion over the cytoplasm keys present;
   - a missing pool side is added only when the imaging names it (score ≥ 0.3, i.e. any HPA term / CYCLoPs ≥ 0.3) or
     its share would be ≥ 0.5; transmembrane proteins are skipped.
8. **round**: 4 decimals, shares < 1e-4 dropped (the rule-set convention; was 3 / 0.01 Jurkat, 4 / 0 yeast).

Not in B5 (stated so nothing is read into the absence): plan §3a.5 membrane-mechanism step (no mechanism feature table
yet; v6 R2/R3/R5 keep their mechanism gates), R8a imaging-fitted fractions, the R8c background term b, surface
proteomics, reliability-weighted means within a class (each cell has one source per class).

## 2. No pack-dependent rules: v5 P and v6 R1 re-expressed (Jurkat; plan §3b)

Under preset `common` the two rules that gate on PLACED shares of earlier packs become rule steps of kind `rule` (no
frozen pack input; `frozen_overrides` stays for `legacy_v7` only):

- **v5.P nuclear pool**: gene universe = genes of the current table (was: v4-pack loc_compare); selection unchanged in
  spirit, on INTENDED shares of the species: single-gene, not transmembrane, mass < 60 kDa, signal_frac < 0.5, nuclear
  support (an HPA or UniProt row with key nucleus in this build's `localisation/evidence.tsv`), not HPA-JURKAT-imaged,
  not v4-curated, not actin/tubulin, **no imaging N/C calibration in step 7**, intended nucleus share < 0.05 and
  cytosol + unassigned ≥ 0.8; target nucleus = Vn / (Vn + Vc) of the cytosol + unassigned + nucleus pool (r = 1).
- **v6.R1 N/C from imaging categories**: gene universe = genes of the current table (was: v5-pack genes); the gene gate
  "placed nucleus + cytosol ≥ 0.8" becomes the species gate "intended nucleus + nucleolus + speckle + cytoplasmic keys
  ≥ 0.8 of the species' weights"; categories, hold-out consensus, Itzhak tie-breaker, the |log2| > 1 and direction gates
  unchanged.
- Both use Vn / Vc of this build's geometry measures (was the v4 constants 254.5 / 385.5 µm³).
- **Capacity**: a rule may not raise the intended nucleoplasm protein density above the intended cytosol density
  (mass of the nucleus / cytoplasm keys over Vn / Vc, before the rule); otherwise all its raises are scaled by one
  factor (reported). At pack time the B1 pilot thinning factor of `nucleoplasm` is reported in the verdict (capacity of
  the packer).

Still reading frozen v4/v5-pack tables after B5 (gene universes, annotation fields, Tier-2 calibration; out of B5 scope,
listed for review): v5.C_cortex (also gates on the v4 pack's placed PM share), v5.N_nucleolus, v6.R4, v6.R2/R3 (tier2_v5),
v6.R5, v6.R6a.

## 3. Mechanism (not effective on main until the B5 queue job)

`inventory.params.localisation.preset` selects the preset per cell; main keeps `legacy_v7` / `legacy_x5`; the cell data
of the common preset sit in the configs and are inert under the legacy presets. The B5 queue job checks out the merged
main commit, runs `scripts/localisation_switch.py` (preset → `common`; Jurkat `localisation_sources` term_map strings →
the TSV paths actually read), commits that change in the launch worktree ("switch commit", recorded in the verdict) and
builds it. After acceptance the same switch commit is merged into main.

## 4. Acceptance per cell (B5 build vs reference = the B2 build of the same cell)

Reference: the build named in `<cell>/_dev_unify/b2/full_build.status` ("RUNNING build …", complete); if absent, the
latest complete earlier B build (B1 → B4 → B3/B6 → A9) and the verdict says which. Metrics from `validate/suite.json`
(common suite) of both builds unless stated.

| # | Check | Pass |
|---|---|---|
| L1 | C1 headline per reference set, intended AND placed: Jurkat `agree.unweighted` for opencell, hpa_consensus_all, hpa_jurkat, hpa_jurkat_holdout_half; yeast `tophit.unweighted` for SWAT_native, SWAT_NOP1, hyperLOPIT (CYCLoPs = model input, report only) | ≥ reference − 0.005 each |
| L2 | every other C1 metric (agree, agree_tvd, tophit × unweighted / mass_weighted × intended / placed) | §4 no-regression (± 0.01) |
| L3 | C3 image N/C. Yeast: suite `image_nc.test.spearman` and `.auroc_nuclear` (placed copies, YeastRGB C-SWAT). Jurkat: the suite C3 reads the frozen v7 Tier-2 table and cannot respond, so B5 computes **C3i**: per Tier-2 gene (`tier2_v7.tsv`, measured `nc_hpa`), model log2 N/C = log2[(copies in nucleoplasm + nucleolus + speckle) / (Vnucleoplasm + Vnucleolus + Vspeckle)] − log2[copies in cytosol / Vcytosol] from `validate/loc_model.parquet`, intended and placed, Spearman on the hold-out split (gate) and all genes (report); identical code on both builds | ≥ reference − 0.01 |
| L4 | yeast benchmark (plan §3d): SWAT_native intended tophit, hyperLOPIT intended tophit, image N/C test Spearman | ≥ 0.655, ≥ 0.838, ≥ 0.80 |
| L5 | mass moved per rule (§5) | ≤ 2 pg AND ≤ 2 % of the cell's total protein; else **review** |
| L6 | no regression elsewhere: suite no_regression vs reference for C5-C17; C11 density changes are accepted when explained by the intended mass change of the region (|Δρ − Δρ_intended| ≤ 2 mg/ml); C10 placed mass fraction ≥ reference − 0.003 | as stated |
| L7 | build complete (snakemake 100 %), C16 build integrity, inventory compare not required | must |

The 2 % term in L5 is this file's reading of plan "2 pg" for a 3.5 pg yeast cell (Jurkat: 100 pg → 2 pg, unchanged;
yeast: 0.07 pg). Verdict: **ACCEPTED** (all pass, nothing to review), **NEEDS REVIEW** (all gates pass, ≥ 1 rule above
L5), **NOT ACCEPTED** (any of L1-L4, L6, L7 fails; reported per item). Written to
`<build>/validate/compare_B5/verdict.json`.

## 5. Mass-moved report (every B5 build and the dry run)

`compare_B5/mass_moved.tsv`: one row per rule: rule, level (combiner per protein / rule set per species), proteins or
species changed, proteins whose top compartment changes, mass moved (pg) = Σ mass × ½ Σ|w_after − w_before|, share of
total protein, review flag. Combiner attribution by a fixed ablation on this build's own evidence table, each step
switching one default on top of the previous: legacy preset → (a) key selection + annotation source (yeast) / membrane
restriction (yeast) → (b) split score → (c) combine precedence classes → (d) new_from_imaging → (e) nc_model →
(f) round; then v5.P and v6.R1 (from the rule sets' change tables). Protein mass = copies × mass_da of
`inventory/v0/inventory_v0.tsv`. Also: proteins whose top compartment differs between the reference and the B5 build.

Dry run (inventory only, no packing; CPU ≤ 3 cores): `wholecell_builds/{jurkat,yeast}/_dev_unify/b5/dryrun_<commit7>/`:
the inventory with `legacy_v7` / `legacy_x5` and with `common` from the same code and inputs, the §5 report, and the
intended-level metrics available without packing (Jurkat C1 intended for the four sets and C3i intended; yeast C1
intended for SWAT_native / SWAT_NOP1 / hyperLOPIT / CYCLoPs and the intended image N/C on the test split). Placed
metrics need the build.

## 6. Expected changes (written before the dry run)

- Jurkat v0: the Itzhak re-split of every soluble protein is replaced by r = 1 or the HPA category ratio (C2 history:
  the v2 Itzhak re-split lowered OpenCell agreement 0.786 → 0.752); HPA-only additional keys (0.3) and HPA-only nucleus /
  cytosol keys leave the blend (the N/C model re-adds pool keys named by HPA); expected C1 changes of order ± 0.01;
  nucleoplasm share of soluble protein may rise (r = 1 at Vn / (Vn + Vc) ≈ 0.40). v5.P / v6.R1 become small (the N/C
  model has set most targets at v0).
- Yeast: multi-localised proteins lose the rank order (equal shares), mitochondrial proteins with a cytosolic secondary
  location lose the 0.1 down-weighting, imaging carries half the group shares; SWAT / hyperLOPIT top-1 may move either
  way by ~0.01-0.03 (risk stated: L4 may fail). Image N/C Spearman should stay (same k, current volumes).
- No change in copies, assemblies, geometry or packer settings.
