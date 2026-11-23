# Freeze criteria: paper builds Jurkat v8 and yeast next (unification plan §5c) — DRAFT

Status: **DRAFT, 2026-11-19, revision 2026-11-19 evening** (Step 3 preparation; branch `freeze-prep`, not merged to
main; Pedro's decisions of the evening are folded in below and listed in §0). Written before any freeze build. No freeze build is launched or queued until (a) this draft is frozen (md5 recorded in PLANNING_NOTES by the
main session), (b) the surface-proteomics recommendation exists, and (c) the pending inputs of §2 are merged or waived
by a dated DECISION. Changes after freezing: dated addenda only (`freeze_v8_addendum_<date>.md`), as for every step.

Sources: `docs/reference/history/unification_plan.md` §4, §5b, §5c; PLANNING_NOTES 2026-11-16..19 (decisions cited
inline); `docs/unification_criteria/*` (B1–B6, B6b, addenda and records); `docs/unification_progress.md`.

Machinery (prepared with this draft; nothing writes into existing builds):
- `scripts/freeze_check.py` — preflight on a code commit, then the build-level verdict
  `<build>/validate/compare_freeze_v8/verdict.json` (gates G0–G9 below). Reuses `scripts/step2_combined_check.py`
  (Step-2 step items, accepted-deviation mechanism, C10 tile sums, mass-consistency of densities) and
  `wholecell.validate.suite.no_regression` (plan §4 tolerances).
- `wholecell_builds/{jurkat,yeast}/_dev_unify/freeze/freeze_build.sh` — queue-ready (refuses unless `QUEUED=1`):
  preflight → clean worktree of `COMMIT` → `scripts/build_isolated.sh <cell> paper 24` → suite (if the DAG did not
  run it) → `freeze_check.py` → release steps (§6) **only** on verdict ACCEPTED.

---

## 0. Revision 2026-11-19 evening — Pedro's decisions (PLANNING_NOTES DECISION 2026-11-19)

| # | Decision | Where it enters | Prepared |
|---|---|---|---|
| D1 | **Replicates: 1 seed per cell** (seed 0 = the paper build; no replicate-seed builds) | §3; G9 report `replicates: 1 (decided)`; preflight refuses `--replicates` ≠ 1 | config unchanged (nothing to add) |
| D2 | **Jurkat SEC-MS wiring merged** (final combined `jurkat_final`) **with HSPD1 excluded from the SEC apex check** (HSPD1 heptamer kept: share from Fossati 2021, = the Step-2 row) | §2a, §2b, §4 E4 (bounds recomputed with the exclusion), G1 `secms_jurkat` (requires the HSPD1 entry) | branch `step3-secms` (caee4f2 + 906b9b9); merge plan `step3_secms_merge_plan.md` on that branch |
| D3 | **Yeast anchor set v3** (v2 minus PDC1, EFT1) | §2a | branch `yeast-anchors-v3` b7f49dc; the yeast SEC-driven assembly (41d3ed5) waits for its fix (diagnose why the assembled share moves away from SEC) |
| D4 | **PDGFRB capped to its RNA-based estimate** (65,762 → 41 copies) | §2a, §4 E8, §5b; G1 `copy_fixes`, G7 E8 | this branch: rule `v8.copy_fixes` = one-gene post-scaling `adjustments` entry (config only; `copies.py` untouched); dry run §5b: nothing else changes |
| D5 | **SOS: keep the measured Jurkat value** (SOS1+SOS2 = 2,080 in the pack; v0 606 + 868) in the structural model, **not raised**. Simulations use both values explicitly: Jurkat measured (2,080) and primary T cells (central 6,500; literature 6,200–14,500, Rieckmann 2017 scaled to Jurkat) — the difference may be real Jurkat-vs-primary biology relevant to activation and is stated wherever SOS enters a simulation | §4 E8 (SOS unchanged is gated), §5b | nothing in the build; simulation configs carry both values (`data/reference/synapse_signalling/sos_copies.tsv`) |
| D6 | **Export the Jurkat viewer tiles (YCT2)** in the freeze build | §2b, §4 E10, G0 | this branch: `export.tiles.enabled: true` in jurkat.yaml (jurkat_v7_port keeps `false`) |
| D7 | **Surface proteomics only for (a) presence evidence and (b) the PM topology fix of 13 units (168,347 placed copies)**; surface intensities NOT used for PM copy numbers (the PM-share fit did not validate) | §2a, §4 E6 / E7, G1 `surface` (mode presence) / `pm_topology`, G7 | **pending its merge**: adapter on branch `surface-proteomics` (b78bb4b, e2ec4ab; not edited here); topology fix not written yet |

---

## 1. Definitions

| Item | Jurkat v8 | Yeast next |
|---|---|---|
| Content | A9 + B1–B6 (+B6b, B3) + logged fixes + SEC-MS + surface proteomics + imaging localisation (§3c of the plan; = B5 as accepted) | A9 + B1–B6 (+B6b) + yeast SEC-MS (+ G. Magnusson rounds only via the shared pipeline) |
| Config | `config/cells/jurkat.yaml` at `COMMIT` | `config/cells/yeast.yaml` at `COMMIT` |
| Build id | `<date>-<commit7>-s20261115-paper` (`build_isolated.sh jurkat paper 24`) | `<date>-<commit7>-s0-paper` |
| Paper seed ("seed 0") | master 20261115 with the pinned legacy streams of the config (`seeds.legacy`) | master 0 with the pinned X5 streams |
| Reference build (no-regression) | `jurkat/20261118-40d7967-s20261115-step2` = the **accepted** combined Step-2 build (verdict `validate/compare_step2_v3/verdict.json`, ACCEPTED with 7 recorded deviations: `step2_combined_accepted_deviations_2026-11-19.{md,json}`); suite `…/validate/suite.json` | **yeast reference = the accepted yeast combined Step-2 build (pending)**: the step2ds verdict (downstream rerun of `yeast/20261119-b16bd09-s0-step2`) does not exist yet; its path and any accepted deviations go into an addendum before the yeast freeze build |
| Tag | `paper-jurkat-v8` | `paper-yeast-<id>`; proposed `<id>` = `u1` (first unified yeast paper build); Pedro may rename before queueing (`YEAST_ID`) |

Freeze-only change = a change that is in the freeze build but not in its reference build. Every freeze-only change has
an EXPECTED entry in §4; a suite regression is EXPECTED only if a §4 entry names it and its bound holds.

## 2. Exact inputs

### 2a. Code commit

`COMMIT` = one commit **on main** (ancestor of main at launch; never a branch tip), clean (`BUILD_STAMP.dirty =
false`), containing every commit below. The preflight (`freeze_check.py --preflight`) checks ancestry and config facts.

| Input | Commit(s) | Status 2026-11-19 |
|---|---|---|
| Step-2 state of the accepted combined build | 40d7967 (Jurkat build commit) and its ancestors | merged |
| Step-2 verdict tooling + accepted-deviation records | 257ec59, ab9927c | merged |
| Freeze packer fixes: orientation (validator-consistent gradient, |grad| < 0.5 or normal > 90°) + nucleolus field cap | 3227713 | merged (accepted 2026-11-18/19) |
| Yeast amphiphysin pool fix | b16bd09 | merged |
| Yeast B5 overrides: Sir2/3/4 (b5.nuclear_silencing) + audit rows (b5.audit_curated 28, b5.builder_anchored 13) | b691e9a, 6978b43 (record a12bdac) | merged (approved 2026-11-19) |
| Yeast downstream integration fixes (export masses, engine layout, plugin env) | a1fb253, 9df5a94, f12684d, 845a2c2 | merged |
| **Jurkat + yeast SEC-MS wiring** (final combined SEC-MS: Jurkat `jurkat_final` read by the assembly rules; yeast tables) + **HSPD1 apex_check_exclude** (D2) | branch `step3-secms` caee4f2 (was main 24efb2c, reverted on main ef09448) + 906b9b9 | **agreed (D2), PENDING merge** (merge plan on the branch) |
| **Yeast anchor set v3** (D3) | branch `yeast-anchors-v3` b7f49dc | **agreed, PENDING merge** (with or before the yeast assembly commit) |
| **Yeast SEC-MS → assembly states** (same common rules / cut-offs as Jurkat) | branch `yeast-anchors-v3` 2fceeba, 41d3ed5 + the pending fix | **PENDING** (fix: why the assembled share moves away from SEC; then OK) |
| **Surface proteomics adapter, mode presence** (`adapter.surface: surface_proteomics`; evidence only, D7) | branch `surface-proteomics` b78bb4b, e2ec4ab | **agreed (D7), PENDING merge** + the config lines of §2b |
| **PM topology fix** of the 13 glyco-captured units (§5d, D7) | not written (library re-orientation; not `copies.py`) | **agreed (D7), TBD** (commit + check script) |
| **PDGFRB copy fix** `v8.copy_fixes` (§5b, D4) | this branch (jurkat.yaml `inventory.copies.adjustments`; config only) | **agreed (D4), prepared** |
| **SOS** (D5) | none (measured value kept) | decided: no change |
| **Jurkat YCT2 tiles** `export.tiles.enabled: true` (D6) | this branch | **agreed (D6), prepared** |
| Jurkat `validate.suite.reference` = the Step-2 build (this branch) | freeze-prep | prepared (§7) |

### 2b. Configuration facts checked by the preflight

Jurkat: geometry `em_cell1_v8r`, manifest sha256 `6c86b1688722f4320520abbb46a66fb7e9e97c95576e0186ebf3909b9c562fa1`;
`geometry.npc` / pre-placed `npc.n_pores` = 2,696 and assembly NPC cap 2,696 (B4 J1); `pack.organisation.enabled` =
true with `field_cap_regions: [nucleolus]`; localisation preset `common` (B5); assembly preset `legacy_v7` with
`complex_portal.copies` = stoichiometry_cap (measured_free_pool OFF); B1 `COMMON_DEFAULTS = {}` (no pack preset pin);
library mode build (B3); SEC-MS: `assembly_steps.overrides.secms.source_id = jurkat_final` and
`assembly_steps.data.apex_check_exclude` lists HSPD1 P10809 (**pending** gate `secms_jurkat`); surface:
`inventory.adapters.surface = surface_proteomics` with `surface_proteomics.mode = presence` (**pending** gate
`surface`); copy fixes: `inventory.copies.adjustments` has rule `v8.copy_fixes` (gate `copy_fixes`);
export: `export.tiles.enabled = true` (D6).

Yeast: geometry `analytic_x5_s0_f5f3a5c723b2_memhalf1`, manifest sha256
`230b96d91559b6aca5557bae68cfafeed227c24be3aaa76879831acd29dd8f5b`; `pack.preset = [legacy_x5, b1_s4y]` (B1 C3 pin);
`field_cap_regions: [nucleolus]`; localisation preset `common` with `weight_overrides`; SEC-MS adapter on and the
assembly reading it (**pending** gate `secms_yeast`).

### 2c. Frozen inputs

`<build>/frozen_inputs.tsv` must equal the reference build's list (name, algorithm, expected checksum; all `ok`)
except for these **allowed additions/changes** (any other difference fails G0):
- Jurkat SEC-MS: `inv_sec_lab2026_jurkat_nt`, `inv_sec_lab2026_jurkat_tcract24h` (+ any `inv_sec_*` listed by the
  merged SEC-MS wiring commit);
- surface proteomics: names `inv_surface_*` (dataset `data/proteomics/umpublished_datasets/surface_jurkat_2026/`,
  unpublished, local only);
- yeast SEC-MS: `inv_sec_*` / `inv_secms_*` names of the yeast SEC-MS commit;
- copy fixes: none — `v8.copy_fixes` is a config rule (no data file, no frozen input).
The geometry packages are not rebuilt; their manifests must match §2b (G0).

## 3. Seeds and replicates (decided parameter)

- The paper build is **seed 0** (§1). It alone carries the tag and is used by the simulations.
- `REPLICATES` = **1 seed per cell** (**DECISION Pedro 2026-11-19 evening, D1**; earlier default 2026-11-18: no
  replicate-seed packing runs only to estimate variance). No replicate builds are prepared or queued; the preflight
  refuses `--replicates` ≠ 1. The verdict records `replicates: 1 (decided)`.
- Seed noise of the reference: the Step-2 builds are single-seed; a regression that the yeast port seed spread
  (s0/s1/s2, `yeast/_dev_unify/b6/seedspread_suite/`) covers is EXPECTED for yeast, as in Step 2.

## 4. Expected changes vs the reference build (one entry per freeze-only change)

| # | Freeze-only change | Cell | Expected (gated unless "report") |
|---|---|---|---|
| E1 | Orientation fix (3227713) | both | C9 `sign_errors_per_1e6` = **0** (Jurkat reference 1.09 = 3 flips, KNOWN deviation #1 of the Step-2 record); C9 median tilt ≤ reference + 0.5°; placed mass fraction overall ≥ reference − 0.003; partner fraction report (Jurkat tiles 2857/6113/6381: 0.203 → 0.204) |
| E2 | Nucleolus field cap (3227713; `field_cap_regions: [nucleolus]`) | both | Jurkat nucleolus placed mass fraction **≥ 0.86** (reference build 0.807 = KNOWN deviation #2; B3 build 0.864; 43-tile check 0.866); yeast nucleolus ≥ reference − 0.003 (quick subset 0.926 → 0.949, report); C11 `density.nucleolus` change mass-consistent with the placed-mass change (± 2 mg/ml, B1 K6 rule); every other volume region placed fraction ≥ reference − 0.003 |
| E3 | Yeast B5 overrides (b691e9a, 6978b43; 44 rows, 0.21 % of protein mass) | yeast | C1 `SWAT_native.intended.tophit.unweighted` ≥ 0.6682 − 0.005; `hyperLOPIT.intended.tophit.unweighted` ≥ 0.8483 − 0.005; C3 `image_nc.test.spearman` ≥ 0.8345 − 0.01 (values of the approved dry run); B4.Y1 amphiphysin 293 of 586 dimers placed; chromatin sites ≥ 0.99 × legacy (report); SIR3 at silent loci (report, dry run 89); residuals RRB1, MTR4, SMT3, TIF3, ARC18 (report). The yeast reference (pending) predates these rows: C1/C3 changes in the good direction are improvements |
| E4 | Jurkat SEC-MS final combined (`jurkat_final`) **with HSPD1 excluded from the apex check** (branch step3-secms caee4f2 + 906b9b9, **pending merge**) | Jurkat | inventory complexed share **0.5253 → 0.5371 ± 0.001** (Complex Portal + homo-oligomer + AFDB homo/heterodimer mass / total); protein mass moved between free and assembled **2.53 ± 0.10 pg** (1,071 proteins; 1.86 to assembled, 0.67 to free); Complex Portal assembled mass unchanged (± 1e-6 pg; ± 0.05 pg with the E8 copy fix); HSPD1 heptamers = the reference row (594,934 × 7azp, share 0.8418; report); pack-table species: copies changed ~1,920, +97 / −200 species (report); top movers IMPDH2 +0.11, PRDX1 +0.09, KHDRBS1 +0.09, GAPDH −0.07 pg (report); `measured_free_pool` OFF. Numbers: dry run `jurkat/_dev_unify/step3_secms_hspd1/dryrun.json` (md5 db02d46f…; baseline = Step-2 inventory; the wiring alone gave 0.5352 / 2.72 pg with HSPD1 −0.19 pg) |
| E5 | Yeast SEC-MS → assembly (**pending**) | yeast | **TBD**: bounds written by the yeast SEC-MS addendum from its dry run (diagnostic today: model 58.3 % assembled vs SEC 83.6 % (k2) / 63.0 % (k4), per-protein r 0.07); Complex Portal stoichiometry cap unchanged; `measured_free_pool` OFF |
| E6 | Surface proteomics, **presence evidence only** (`adapter.surface: surface_proteomics`, mode presence; branch surface-proteomics, **pending merge**, D7) | Jurkat | `inventory/surface.tsv` filled (source `surface_glyco_jurkat2026`, ≤ 596 proteins, baseline NAct + Act_0, detection n/N, topology / n_term_side), `pm_fraction` empty; **nothing else changes**: copies, v0, localisation, assembly, ingredients and pack table identical to the build without the adapter (checked once on the merge commit by an inventory-only dry run on/off; in the freeze build by G7 E8 SOS, G6 Complex Portal and G8); new frozen input `inv_surface_jurkat_glyco2026` (md5 e092247b…, redistributable no; §2c) |
| E7 | PM topology fix of the **13 glyco-captured units** (§5d, D7; **TBD**, written after the surface merge) | Jurkat | library re-orientation (ectodomain exoplasmic) of exactly CPX-1802 ITGB1/ITGA4 55,700; TFRC_afx2 51,925; LNPEP_afx2 21,217; CPX-1803 ITGB1/ITGA6 13,443; CD7_x2 10,506; SLC4A7 7,727; CPX-1794 ITGB1/ITGA5 6,494; SLC44A1_x2 787; CPX-26677 IGF1R 168; CPX-9541 IRAG2/SUN2 167; PVR-CD226 158; NECTIN2-NECTIN1 31; CPX-937 NOTCH1 24 (= 168,347 placed copies, `data/derived/surface_jurkat_2026/topology_flipped_units.tsv` rows with `detected_surface`); topology check (S0 `unit_geom`) re-run on all PM units of the freeze library: **0 of the 13 flagged**, the other 29 v7 flags unchanged and listed (report; TRBC1_x2 13,875: side by TM alternation, outside the agreed 13 unless decided); copies unchanged (no `copies.py` change); C9 unaffected (sign = membrane normal); C11 PM / cortex densities no regression; placed PM mass ≥ reference − 0.003 |
| E8 | PDGFRB copy fix `v8.copy_fixes` (§5b, D4); **SOS not raised** (D5) | Jurkat | PDGFRB copies = round(its RNA-based estimate in copies) = **41** (65,762 before; estimate recomputed by G7 from the build's copy summary: min(fit, floor) × 0.25 / Σ ppb × molecules); **no other protein's copies change** (copies.tsv / v0 identical except PDGFRB; total protein 99.9865 pg); downstream only the PDGF family (CPX-2891 163 → 41, CPX-2886 2 → 0, CPX-2879 0 → 40, free PDGFC 0 → 81, PDGFA 0 → 1; Complex Portal mass −6e-5 pg); SOS1 / SOS2 copies = reference exactly (606 / 868; pack SOS1+SOS2 2,080, report) |
| E10 | Jurkat YCT2 viewer tiles (D6) | Jurkat | `export/web/` written by export step `tiles_v2` (via `all`): `encode_check.json` blocks `all`, `count_equal`, `pos_ok`, `rot_ok` true (G0); first whole-cell Jurkat encode (only a 1-block test exists: `_v7_reference/export/tiles_v2_block_test`, 1.13 M instances, 16 MB): expect ~0.9 G packed + pre-placed instances, roughly 12–15 GB of shards (report size and time); the MANIFEST of §8 covers it |
| E9 | Unchanged by design | both | geometry packages (sha §2b); NPC 2,696 pores = placed = inventory (B4 J1); ledger `rule_id` column (no numeric rule ids in `version_changes.tsv`); localisation inputs and B5 preset (Jurkat identical to the reference except species added/removed by E4/E6/E8); pre-placed builders (Jurkat B4 items J1–J8 as in the Step-2 verdict) |

Not in the freeze builds unless decided: TCR microvillus-tip enrichment (optional, Pedro). Jurkat export tiles: decided
(D6, E10).

## 5. v8 backlog check (2026-11-19, no builds; CPU ≤ 2)

### 5a. NPC density / pore count — DONE (B4 J1), in main's Jurkat config
Main's `jurkat.yaml`: geometry `em_cell1_v8r` (B6b on em_cell1_v8p = em_cell1_v8 + 1,696 pores), `npc.n_pores: 2696`,
`npc_copies: 2696`, `complex_copy_caps_b4.tsv` NPC = 2,696. Rule J1: 11 NPC/µm² (Otsuka 2016) × 245.1 µm² NE
midsurface = 2,696; Nup anchor = 2,240 pores, ratio 0.83 (∈ [0.5, 2]). Step-2 build verdict B4.J1 PASS: config =
package = placed = inventory NPC = 2,696, min separation 150.07 nm. Freeze: re-checked by G2 (B4.J1 must PASS).

### 5b. SOS copies and the PDGFRB 64k check — decided 2026-11-19 evening (D4, D5)
Current inventory (copy step of the Step-2 build, unchanged on main): **PDGFRB 65,762 copies** (v0; 65,432 free +
CPX-2891 / CPX-2886), iBAQ 39,244 ppb (Geiger 2012, 3 replicates) but **RNA TPM 0.195 (call `low`)**; Geiger detects
PDGFRB in all 11 cell lines (per-line medians ~37,000–168,000 ppb; incl. Jurkat, K562, HEK293, which do not express it) → a protein
group / peptide assignment artefact, not Jurkat biology (MS/RNA log-ratio z = +2.7). Its mass is small (0.0135 pg) but
it is a T-cell-implausible PM receptor (dasatinib off-target, application B3). Other MS-high/RNA-low genes are mostly
paralog peptide sharing (histone variants, HSPA1A/HSPA6, ACTA1, UBE2L5, EIF2S3B, LGALS9B/C; 76 genes with TPM < 1 and
> 10k copies) — not in scope.
**SOS:** SOS1 606 v0 copies (replicates 177 / 2,202 / 362 ppb: 12-fold spread) → 1,212 in CPX-395 (HRAS₂–SOS1;
stoichiometry cap = up to 2 × the limiting subunit, by design); SOS2 868 (838 free + 30 SOS2–HRAS). SOS1+SOS2 = 2,080
in the pack vs 6,200–14,500 in primary T cells scaled to Jurkat (Rieckmann 2017, reference panel), RNA log-log fit
~12,400; Geiger-scaled 1,697 (`data/reference/synapse_signalling/sos_copies.tsv`; S2 rerun central 6,500). Still low.
**Decision (Pedro 2026-11-19 evening):** PDGFRB capped to its RNA-based estimate (D4); SOS NOT raised — the measured
Jurkat value stays in the structural model (D5; simulations use 2,080 and the primary-T value explicitly). The earlier
SOS1/SOS2 × 4.41 proposal is withdrawn. RASGRP1 (734 copies, one replicate, RNA fit ~23k) noted, not changed.
**Implementation (this branch; config only, `copies.py` untouched, no new module needed):** one entry of the copy
step's existing `adjustments` hook (step 9, after the total-mass scaling) in jurkat.yaml `inventory.copies`: rule id
`v8.copy_fixes`, kind `family_target`, members `{genes: [PDGFRB]}` (a one-gene "family"), target count
40.785189997011194 = the RNA-based estimate in copies (24.339126 ppb = min(10^(2.73207 + 0.76795 log10 0.195346),
97.3565) × 0.25, the value the copy step gives an undetected RNA `low` call, / Σ ppb 952,062,291.18 × 1,595,375,332.6
molecules; derivation in the config comment). `jurkat_v7_port.yaml` sets `adjustments: []` (v7 recipe unchanged).
Why not the earlier `per_protein_scaling: factor_table` proposal: it acts before the total-mass scaling, so the 100 pg
total renormalises every other protein (+1.35e-4) — dry run `_dev_unify/freeze/copy_fixes_pdgfrb/`: PDGFRB at the
estimate but 5,913 other proteins ±rounding, 47 pack species added / 48 gone, 499 complexes changed. The adjustment
changes PDGFRB only; its 0.0135 pg are not redistributed (total 99.9865 pg, within the 80–130 pg bounds of the
100 pg DECISION). The target is a number in the config (valid for the frozen inputs); G7 E8 recomputes the estimate
from each build's own copy summary and fails if PDGFRB ≠ round(estimate). Tests:
`tests/test_copies.py::test_jurkat_copy_fixes` (target = estimate; one-gene adjustment leaves every other protein
unchanged on a synthetic cell), `test_configs`.
**Dry run** (inventory only, `jurkat/_dev_unify/freeze/copy_fixes_pdgfrb_adjust/`, dryrun.json md5 8d338a43…, vs the
Step-2 inventory `step3_secms/baseline_main_a12bdac`; freeze-prep differs from it in the inventory only by this rule):
PDGFRB 65,762 → **41** copies (0.01354 → 0.0000084 pg; ledger `adjust.v8.copy_fixes`, factor 6.2019e-4); the other
19,893 proteins: **0 copies changed**, ppb identical; SOS1 606 / SOS2 868 unchanged; imputation fit identical;
v0 differs in the PDGFRB row only; homo-oligomer / AFDB homo- and heterodimer tables identical; Complex Portal: 3
complexes change, all PDGF (CPX-2891 PDGFRB–PDGFC 163 → 41, CPX-2886 PDGFRB–PDGF-AB 2 → 0, CPX-2879 PDGFC₂ 0 → 40);
ingredients / pack table: only CPX-2891 (163 → 41), CPX-2886 and free PDGFRB (65,432) gone, free PDGFC (81) and
PDGFA (1) added; every other row identical.

### 5c. Ledger `rule_id` — DONE
`wholecell/inventory/canonical.py::_version_changes` reads the `rule_id` column (schema `inventory_schema.md` §
changes.tsv); Step-2 build `inventory/version_changes.tsv`: 3,145 rows with rule ids, 0 numeric (the legacy v6/v7 bug
put `mass_moved_pg` there). Freeze: G2 re-checks (0 numeric rule ids).

### 5d. PM topology — still wrong in the current library; what the surface work must provide
Re-run of the S0 topology check (`jurkat/applications/s0_geom.py::unit_geom`, UniProt TOPO_DOM, read-only) on the
Step-2 build library (`…-step2/library/lib`): **all 42 units are still flagged** (B3 rebuilt the trees but kept the
orientations), 207,337 placed copies in v7 terms. By copies: CPX-1802 ITGB1/ITGA4 55,700; TFRC_afx2 51,925; LNPEP_afx2
21,217; TMEM222 15,093; TRBC1_x2 13,875; CPX-1803 ITGB1/ITGA6 13,443; CD7_x2 10,506; SLC4A7 7,727; CPX-1794
ITGB1/ITGA5 6,494; CYBRD1_afx2 2,054; CPX-26464 SLC9A1/TESC 1,933; NDFIP1 1,722; SLC12A9 951; GDE1_afx2 923;
SLC44A1_x2 787; DSC1_x2 779; 26 more ≤ 168 each (IGF1R, CPX-9541, PVR-CD226, SLC24A1, GRIN3B, PAQR6, SLC44A5, PAQR8,
SCN complexes ×3, SLC4A8, ATP2C2, ABCA1, GPR135, CPX-3282, ABCA7, SLC12A4_x2, NECTIN2-NECTIN1, MADCAM1_x3, CPX-937,
MAL_afx2, TFR2_afx2, KCNK9_x2, GJA1_x12, CNTF-IL6R). Method: 14 "residue" (AF monomer residues), 28 "balance" (bead
balance vs UniProt ecto/cyto residue counts). Some flags are doubtful: TRBC1 / SLC44A1 have UniProt expected-ecto
fraction 0 although TRBC1's constant domain is extracellular (UniProt TOPO_DOM gap) → a blind flip would be wrong.
15 of the 51 genes are detected in the surface glyco-capture set (ITGB1, ITGA4, ITGA5, ITGA6, TFRC, LNPEP, CD7,
SLC4A7, SLC44A1, IGF1R, SUN2, PVR, CD226, NECTIN1, NOTCH1); TRBC1 is not.
**The surface-proteomics work must provide** (for the freeze): (1) per PM unit, the side evidence: glyco-capture
detection = extracellular N-glycan exposure (protein level; per-glycosite/peptide evidence if the raw data allow) to
arbitrate the 42 flags (fix list = flagged AND (surface-detected OR UniProt TOPO_DOM unambiguous); doubtful flags
listed); (2) the PM copy model (copies + intensities + properties + literature anchors) as `adapter.surface` output
`surface.tsv` with per-protein PM copies or PM share and a stated uncertainty, and which proteins it changes;
(3) the time-course use (resting = NAct for the paper build); (4) the topology check re-run on ALL PM units of the
freeze library (not only the 42 v7 units) with the list of units to re-orient. The fix itself (library re-orientation
of the listed units, a library/config change, not `copies.py`) is written after that recommendation (E7).

## 5e. Imaging localisation
Plan §5c item 4 (imaging localisation after the R8c nuclear-floor decision) is B5 as accepted (HPA background
b = 0.19 from OpenCell live images, addendum 2026-11-18); nothing further for the freeze. TCR microvillus-tip
enrichment: optional, not included unless Pedro decides.

## 6. Acceptance gates (verdict `compare_freeze_v8/verdict.json`)

| Gate | What | Rule |
|---|---|---|
| G0 integrity | build complete (snakemake 100 %), `BUILD_STAMP.dirty = false`, commit = `COMMIT` on main, config facts §2b (incl. tiles enabled), geometry manifests, frozen inputs = reference ± §2c allowances, suite C16 (`manifest_stale_files` 0, `depends_on_geometry_mismatch` 0 / None), suite.json present; Jurkat: YCT2 `export/web/encode_check.json` (E10) | all PASS |
| G1 pending inputs | gates `secms_jurkat` (source `jurkat_final` + HSPD1 in `apex_check_exclude`), `surface` (adapter on, mode presence), `pm_topology` (`--satisfied` with the fix commit), `copy_fixes` (rule `v8.copy_fixes`) (Jurkat) / `secms_yeast` (yeast) | each satisfied in the config **or** waived by a dated DECISION passed as `--waive <gate>=<decision ref>` |
| G2 step items | the Step-2 build-level items (B6/B6b, B3, B4, B1, B2, B5 of `step2_combined_check.py`, their own references) on the freeze build; Step-2 deviations of class ACCEPTED / EXPECTED (record 2026-11-19, #3–#7) carried forward with their limits; KNOWN deviations (#1 C9, #2 nucleolus) **not** carried (their fixes are in) | no FAIL left; any other FAIL needs a dated `--accepted` entry |
| G3 E1 orientation | §4 E1 | C9 sign errors = 0 |
| G4 E2 nucleolus | §4 E2 | Jurkat ≥ 0.86; yeast ≥ reference − 0.003; other regions ≥ reference − 0.003 |
| G5 E3 yeast B5 overrides | §4 E3 | bounds above |
| G6 E4/E5 SEC-MS | §4 E4 (Jurkat) / E5 (yeast, TBD) | bounds above; NA if the gate is waived |
| G7 E6/E7/E8 | surface presence evidence (E6), PM topology (E7), PDGFRB at its RNA estimate + SOS not raised (E8) | E6 / E8 bounds above; E7 judged by the check script of the fix commit (until then NA: blocks unless waived) |
| G8 no regression | suite vs the reference build's suite.json (plan §4 tolerances); each regression EXPECTED only via §4 (E1 C9 improves; E2 C10 placed-fraction rows and mass-consistent C11 nucleolus; E3 yeast C1/C3; E4 C11 rows mass-consistent with the requested-mass change; yeast seed-spread band), else UNEXPECTED | 0 UNEXPECTED, or each covered by a dated `--accepted` entry for this build |
| G9 replicates | `REPLICATES` | report (1 seed per cell, D1) |

Verdict: **ACCEPTED** if G0–G8 have no FAIL / ERROR and no un-waived NA gate; **INCOMPLETE** if the build or suite is
incomplete or the reference is pending (yeast until its Step-2 verdict exists); else **NOT ACCEPTED**.
Statuses as in Step 2: PASS | FAIL | NA (with reason; never counted as passed) | REPORT | ACCEPTED_DEVIATION.

## 7. Freeze config changes prepared on this branch

- `config/cells/jurkat.yaml` `validate.suite.reference` → the accepted Step-2 build and its `validate/suite.json`
  (the DAG's `validate/no_regression.tsv` is then vs the freeze reference). Yeast stays on the X5 port until its
  Step-2 verdict exists (then the same change, addendum).
- Revision 2026-11-19 evening: `inventory.copies.adjustments` = rule `v8.copy_fixes` (PDGFRB; D4) and
  `export.tiles.enabled: true` (D6) in jurkat.yaml; `jurkat_v7_port.yaml` pins `adjustments: []` and
  `export.tiles.enabled: false` (v7 port unchanged); `tests/test_copies.py` (copy-fix test; legacy-recipe test pins
  `adjustments: []`);
  `scripts/freeze_check.py` gates updated (E4 bounds with HSPD1, E6 presence, E8 PDGFRB/SOS, tiles in G0, replicates).
- Everything else that is accepted is already in main's configs (3227713 orientation + nucleolus cap; yeast B1 C3 pin;
  yeast overrides; B2/B5 folded in; B4 J1 pores). Inputs on other branches (§2a: step3-secms, yeast-anchors-v3,
  surface-proteomics) are not copied here; the surface config lines (§2b: `adapters.surface`, `params.surface_proteomics`
  from `config/fragments/jurkat_surface_proteomics.yaml`, mode presence; frozen input `inv_surface_jurkat_glyco2026`) are
  added to jurkat.yaml when that branch is merged.

## 8. Release steps (freeze_build.sh; only on ACCEPTED)

1. `export/MANIFEST.sha256`: sha256 of every file under `<build>/export/` (sub-manifests included, this file
   excluded), written once; its own sha256 into `<build>/FREEZE.json` and the status file, for PLANNING_NOTES (the
   main session records it).
2. Annotated tag on `COMMIT` in the main repository: `paper-jurkat-v8` / `paper-yeast-<id>`; message = build id,
   build dir, verdict path, MANIFEST sha256, criteria md5. Refuses if the tag exists (never moved).
3. `<build>/FREEZE.json` (cell, build id, commit, tag, criteria file + md5, verdict, waivers, accepted files,
   replicates, MANIFEST sha256, date).
4. Build folder read-only: `chmod -R a-w <build>` (the `_code` worktree included).
5. After the freeze any change needs a new build id; simulations read the pinned build (build id + MANIFEST sha256).
Not ACCEPTED → no tag, no MANIFEST, folder stays writable for diagnosis; the main session decides.

Decided 2026-11-19 evening (§0): replicates 1, Jurkat tiles, PDGFRB cap, SOS kept, HSPD1 apex exclusion, surface
presence + topology only. Still open (listed, not assumed): yeast tag id (`u1`); TCR microvillus enrichment; TRBC1_x2
topology (outside the agreed 13 units); the yeast SEC-driven assembly fix (E5).
