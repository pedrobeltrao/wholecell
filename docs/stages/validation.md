# Stage 8: validation

**Purpose.** Measure the finished model and compare it with independent data, so a reader can judge each build
without trusting the pipeline.

**Today.** Folder and CLI stage `validate`. One Snakemake rule per step of `STEPS_ALL` in
`wholecell/validate/stage.py`: `densities`, `crowding`, `loc_compare`, `ribosome_nn`, `placement`, `chromatin`, then
`manifest`; `compare_v7` when `validate.compare_reference` is set (Jurkat); then `suite` (the common suite).

**Inputs.** Pack outputs, the pack table, the geometry, and reference data declared in `frozen_inputs` (`val_*`).

**Outputs** (`validate/`): `densities.json` (mg/ml per region), `crowding.json`/`.tsv`, `loc_compare.tsv`,
`loc_holdout.json` (Jurkat), `ribosome_nn.json`, `placement.json`/`placement_tiles.tsv` (clashes, containment),
`chromatin_frozen.json`, `MANIFEST.tsv`; `<step>.skipped.json` with a reason when a step has no configuration or
reference; `compare_v7/` (verdict against the v7 criteria).

**Config keys** (`validate:`): `nproc`, `steps` (`<step>: false` disables a step), `inputs`, and one block per step:
`densities` (`default_region`, `virtual_to`, `shell_to_volume`), `crowding` (`regions`, `n_boxes`, `box`, `probes`,
`legacy_seeds`), `loc_compare` (`vocab`, `cell_line`, `references`), `ribosome_nn` (`em_input`, ...), `placement`
(`n_tiles`, `clash_nm`, ...), `chromatin`, `compare_reference` (`criteria_doc`, `validate_dir`, `files`).

**Checks.** This stage is the checks. The Jurkat v7 reproduction criteria V1-V8 are in
[../reference/history/v7_reproduction_criteria.md](../reference/history/v7_reproduction_criteria.md); the build
compare is `wholecell/validate/compare_builds.py`.

**Code.** `wholecell/validate/*.py`; `workflow/rules/validate.smk`.

**Jurkat vs yeast today.** Yeast checks without a common definition yet run as plugin steps of this stage (A9):
`yeast_engine_pack_view` (the core pack in the engine file layout), `yeast_features`, `yeast_s4_validate` (stage 4),
`yeast_validate_all` (`validate_all.sh`, about 25 scripts, GPU); they list `before: [validate_suite]` and the suite
adapter `yeast_x5` reads their outputs. The core placement / crowding / densities steps run inside the suite for yeast
(`validate.steps` off as separate rules). `plugins/yeast/repro/compare_x5.py` compares a build with X5. Jurkat Tier-2 imaging checks are not
ported. The common suite (below) reads both.

## Common validation suite (A8, plan section 4)

`python -m wholecell validate --cell C --build B --step suite [--out DIR] [ref_suite=<suite.json>]` (workflow rule
`validate_suite`, after `manifest` and `compare_v7`). Code: `wholecell/validate/suite.py` (registry `CHECKS`: one
function per check id); reference adapters are registered as `validate_adapter` and named by
`validate.suite.adapter` (yeast: `plugins/yeast/validate_adapter.py`, `yeast_x5`).

**Outputs** (`validate/`, or `--out DIR`; with `--out` nothing is written into the build):

- `suite.json`: per metric `check`, `metric`, `value`, `unit`, `reference` (external reference value) and
  `reference_source`, `direction` (`up`, `down`, `target:<x>`, `stable` = two-sided until literature ranges are
  declared, `report`), `tolerance` and `tolerance_rule`, `status` (`ok`, `regress`, `improve`, `n/a`), `reason`
  (every `n/a` says why: `n/a (no data)`, `n/a (phase 2)`, `n/a (port defect)`, ...), `source`, `file`, `legacy` /
  `legacy_source` (the number the legacy code wrote for this build, where defined), `doc` (anchor below).
- `no_regression.tsv`: every metric against the declared reference build's `suite.json`
  (`validate.suite.reference: {build, suite}`; the reference build itself compares with itself). Rule of plan section
  4: a metric regresses if it moves in its bad direction by more than its tolerance. Tolerances: localisation, AUROC,
  Spearman ±0.01; rates max(10 % rel., 1 per 10⁶) (fractions: max(10 % rel., 10⁻⁶)); densities ±2 mg/ml; NN median
  ±0.5 nm; placed fraction ±0.003. Suite defaults where the plan gives none: φ ±0.005, N/C ratio ±0.01, anchor ratio
  ±0.05 in |ratio − 1|, seam index ±0.02 in |index − 1|, concentrations 10 % rel., tilt ±0.5°. A regression prints
  metric, value, reference-build value, tolerance and the doc anchor.
- `legacy_check.tsv`: suite value vs legacy value per metric (the A8 code check), with an explanation column.

**Sources.** `core` = read from this build's core validation step outputs (`validate/*.json`); `core (step run by the
suite)` = a core step the build's `validate/` lacks, run by the suite into `<out>/core/` (caches into `<out>/cache/`;
config `validate.suite.core_steps`, plus `orientation` by default; inputs the port build keeps elsewhere come from the
adapter's `core_inputs()`: yeast library `engine/lib`, packer collision shapes `pack_trees_15.npz`, geometry package);
`common` = computed by the suite with one definition for every cell (`locmetrics.py` for C1, `suite.py` for C3, C6,
C10); `suite (new)` = computed by the suite from the build's tiles; `adapter (yeast code)` = a number written by the
vendored yeast validation and read verbatim; `adapter (legacy Jurkat code)` = a legacy Jurkat table declared in config.
Jurkat legacy values come from the ported validation run on the archived v7 pack (`validate.compare_reference.validate_dir`,
itself code-checked against the stored v7 files). A legacy number of another definition is shown side by side in
`legacy` with the reason in `legacy_note` (column `explanation` of `legacy_check.tsv`); cell-specific findings are
attached to rows by `validate.suite.notes`.

**Seam index** (C10, new). Tiles whose volume is ≥ 0.95 in the seam region (`validate.suite.seam.region`, default
cytosol; volumes from `pack/tiles.parquet`), a seeded sample of ≤ 300 of them; number density of packed instance
centres in the 10 nm face band (centre < 10 nm from the nearest tile face) divided by the density in the tile interior
(centre ≥ 20 nm from every face). 1 = no seam, > 1 over-filled faces, < 1 depleted faces. Also reported: min and max
over the 8 tile parity classes. (Jurkat has no tile ≥ 0.98 cytosol because virtual regions overlap cytosol.) Yeast's
own `placement_checks` seam numbers (sphere-volume fraction 0–5 nm / 15–35 nm) are a different definition and are
reported separately.

| Check | Metric (common definition) | Source (phase 2) | Jurkat v7repro2 | Yeast port s0 |
|---|---|---|---|---|
| <a id="c1"></a>C1 localisation vs references | per gene: agree = 1 − (0.5 share outside the reference + 0.5 share of main reference classes with < 5 % of copies), agree_tvd, tophit; × intended/placed × unweighted/mass-weighted (`locmetrics.py`, used by Jurkat `loc_compare` too) | common; references per cell: Jurkat core `loc_compare` (OpenCell, HPA consensus, HPA JURKAT + held-out half); yeast adapter (model shares + SWAT N-GFP, SWAT NOP1, hyperLOPIT, test split; CYCLoPs report only, it is a model input) | OpenCell agree 0.8314, HPA 0.8326 (unchanged) | SWAT N-GFP tophit 0.6665, agree 0.6257; hyperLOPIT 0.8462 / 0.8187; NOP1 0.7399 / 0.6871: identical to `localization.py` (agree, top1, mass-weighted); agree_tvd new (uniform over reference classes); yeast macro F1 kept as `yeast_macro_f1` (extra) |
| <a id="c2"></a>C2 localisation across versions | agree, agree (mass), tophit per archived version v1–v7 | adapter (legacy Jurkat `imaging_across_versions.py` annot table) | 84 rows; 2 explained (inventory mass_da weights) | n/a (no archived yeast versions run with one code) |
| <a id="c3"></a>C3 image N/C | per gene log2 N/C of model vs image: Spearman; AUROC of the image ratio for model-nuclear (> 1) vs model-cytosolic (< −1) genes (≥ 5 each) | common on per-gene tables: yeast adapter (`image_nc.tsv`, YeastRGB C-SWAT); Jurkat legacy Tier-2 table `tier2_v7.tsv` (HPA JURKAT vs rendered model images, v7 pack) | held-out half ρ 0.476 (n 55), AUROC 0.90; all ρ 0.603 (= tier2_v7.json), AUROC 0.94 (new) | test ρ 0.807, AUROC 0.979 (= image_nc.json); dev 0.823 / 0.977 |
| <a id="c4"></a>C4 image compartment shares | 1 − TVD per imaged gene | n/a (not wired) | n/a | n/a |
| <a id="c5"></a>C5 copy-number anchors | per anchor: median over families of copies / (stoichiometry × structure count); implied total protein | core `inventory/copies/summary.json` | n/a (frozen v7 inventory has no copies step) | NPC 1.399, histones 0.954 (unchanged) |
| <a id="c6"></a>C6 copies vs an independent dataset | proteins with model > 0 and reference > 0: Spearman of copies, median ratio, share within 2-fold | common (core `inventory/copies/copies.tsv`) + adapter reference table | n/a (no second Jurkat absolute proteome on disk) | Das 2026 ρ 0.8083 (n 4,369); legacy 0.8080 on model ∩ Das ∩ Ho (n 4,345); Ho 2018 (not independent) 0.901 vs legacy 0.925 (same restriction) |
| <a id="c7"></a>C7 clashes | (a) overlap > 0.5 nm of the packer's collision shapes per 10⁶ packed (packed–packed + packed–pre-pass), 300 seeded tiles; pre-pass self-overlap reported; (b) bead overlap > 0.3 nm per sampled bead | (a) core `placement` (yeast: run by the suite; shapes = `pack_trees_15`, the yeast packer's; the library trees are 1.9–3.8× fatter and overlap by design: 1.8 × 10⁵ per 10⁶ on 6 test tiles); (b) adapter (stage 4) | (a) 0 | (a) 0 (new); pre-pass self-overlap 6.4 × 10⁵ per 10⁶ pre-pass (built chains, report); (b) 2.52 × 10⁻⁵ |
| <a id="c8"></a>C8 containment | (a) instance centres not allowed for their region per 10⁶ (+ by region kind); (b) beads outside region | (a) core `placement`; (b) adapter (stage 4) | (a) 0 | (a) 1,225 per 10⁶, all shells (3.6 % of shell instances): the yeast packer draws shell centres in 4 nm voxels with ±2 nm jitter, centres reach d1 + 3.5 nm (diagnosed; all in the cytosol), core tolerance d1 + 0.5 nm; (b) 5.78 × 10⁻⁵ |
| <a id="c9"></a>C9 membrane orientation | sign errors per 10⁶ membrane instances; tilt = angle(library +z, membrane normal = SDF gradient at the centre): median, share < 1°, 99th percentile | core `placement` (+ `orientation` re-run on the placement tiles when placement.json predates tilt; signs re-checked equal) | signs 0.357 (unchanged); median tilt 3.3 × 10⁻⁴°, 98.4 % < 1° (new) | signs 0 (new); median tilt 2.0 × 10⁻⁵° vs stage-4 1.6 × 10⁻³° (10,000 whole-cell instances, normal at the projected point), < 1°: 99.94 % vs 100 % |
| <a id="c10"></a>C10 placement integrity | placed / requested copies and mass (per library entry: `libs.json` mass_da, else library bead mass 110 Da/aa, 330 Da/nt); species < 50 % placed; seam index | common (`common.tile_sums` + tiles); yeast extras by adapter | copies 0.9826, mass 0.9644 (unchanged); seam 1.045 | copies 0.9924 (= stage 4); mass 0.978 (new, bead mass); seam 0.981 |
| <a id="c11"></a>C11 densities and crowding | mg/ml per volume region (packed + indexed pre-pass + RNA, region volumes of tiles.parquet); N/C = nucleus group / cytosol; excluded-volume φ (bead model, 4 boxes per region) | core `densities`, `crowding` (yeast: run by the suite; pre-pass by group, `validate.densities.preplaced_group_volume`) | cytosol 162.8, nucleoplasm 135.1 (unchanged); φ cytosol 0.074 | cytosol 179.7 vs stage-4 179.6; nucleolus 161.3 vs 166.8, nucleoplasm 111.7 vs 112.8, matrix 175.9 vs 177.6 (stage 4 adds DNA/RNA threads not in the pre-pass index), cell wall 5.8 vs 3.2 (region = outer 60 nm vs wall labels); N/C 0.662 vs 0.628 (legacy nucleoplasm only); φ cytosol 0.079, matrix 0.096 (new) |
| <a id="c12"></a>C12 ribosomes | NN median vs EM / cryo-ET; concentration | core `ribosome_nn` + suite; yeast adapter (`ribosome_cryoet.py`) | unchanged | unchanged (adapter) |
| <a id="c13"></a>C13 chromatin | Jurkat frozen v8 metrics; yeast Hi-C, WLC, nucleosomes | core `chromatin` (frozen); yeast adapter | unchanged | unchanged (adapter) |
| <a id="c14"></a>C14 diffusion benchmark | probe-level \|log10(D_model/D_meas)\| | n/a (BD v2): waits for the common BD engine (plan §6.1) | n/a | n/a |
| <a id="c15"></a>C15 proximity / interaction | median per-bait AUROC vs PPI sets | adapter (yeast `proximity_test.py`, gap 10 nm, contacts of the GPU contacts step) | n/a (no BioPlex / OpenCell IP-MS on disk; no Jurkat neighbour code) | PCA 0.815, proximity labelling 0.844, XL-MS 0.994 (snapshot + assembly) |
| <a id="c16"></a>C16 build integrity | build stamp; validation MANIFEST re-checked; yeast build_stamp MANIFEST.json; `depends_on_geometry_mismatch` (frozen pre-placed modules declared for another package; B6) | core MANIFEST; adapter (yeast) | unchanged | unchanged |
| <a id="c17"></a>C17 cell-specific extras | labelled, report only | adapter (yeast: tRNA–nucleolus, silent chromatin, DNase, kinase localisation, ER ribosome density, mRNA localisation, condensate negative control, model vacuole; `vacuole_surface` n/a port defect); Jurkat config `extras_na` | fluorescence panels, galleries: n/a (not ported) | 22 values |

**Geometry vocabulary fix (A8 phase 2).** The core validation geometry (`geom.PackageGeometry`) named the codes of every
package with the human engine names; the yeast package (yeast_cell3d_v1: 15 labels, 10 membranes in another order)
needs its own manifest vocabulary (`labels[].engine_code/engine_name`, `membranes[].engine_code`). Jurkat and testcell
packages give the same names as before.

**Still open after phase 2.** C2 for archived yeast builds; C3 AUROC/Spearman from a ported Jurkat Tier-2 (the table
is the legacy v7 output); C4; C5 for Jurkat (copies step); C6 for Jurkat (no second proteome); C7 (b) bead overlaps
for Jurkat; C14 and C15 neighbour statistics with the common BD engine (plan §6.1); the
simbox exporter and YCT2 for both cells (A8 items outside the suite).

**Methods text (draft).** Each build was validated against data not used to build it: protein density per
compartment, excluded-volume crowding, agreement of placed localisation with OpenCell and HPA (including a held-out
half of HPA Jurkat annotations frozen before model fitting), ribosome nearest-neighbour distances against cryo-EM
positions, and geometric integrity (overlaps and compartment containment on a seeded sample of tiles).
