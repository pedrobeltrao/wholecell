# Step 3: SEC-MS drives the yeast assembly states (common SEC-fed rules, aligned with Jurkat)

Written 2026-11-19 (Pedro, Step 3 task) BEFORE any implementation or dry run. Frozen ONCE: what changes, the metrics
and the pass rules below do not change; later findings go to a dated addendum file, never into this one.
Branch `yeast-anchors-v3` (on `step3-secms` caee4f2; anchor set v3 commit b7f49dc). Per coordinator (2026-11-19):
prepared on the branch only, no merge into main in this round.

## 1. What changes (yeast only; Jurkat untouched)

The yeast assembly (`pools: ingredients`, preset legacy_x5) gets the SEC-fed rules of the Jurkat assembly (preset
legacy_v7, `wholecell/inventory/assembly/steps.py` + `wholecell/inventory/secms.py`) with the SAME cut-offs, fed by the
final combined yeast SEC-MS (`yeast_final` = lab NaCl + Swamy 2022 primary, YPD + Pang 2020 fill;
`toolkits/secms/datasets/update2026_yeast{,_combine}.yaml`), calibrated with yeast anchor set v3 (v2 if v3 is rejected):

| rule | Jurkat (legacy_v7) | yeast now (legacy_x5) | yeast new |
|---|---|---|---|
| R1 homo-oligomer share (PDB / audit state n) | oligomer copies = free x a_k, k = 2 (n <= 3), 3 (n <= 5), 4 (n > 5); no SEC row: 0.5; SEC apex check: rejected when apex / (n x monomer) outside [0.4, 2.5] and a_4 > 0.2 | all copies of the ingredient packed as the oligomer (`copies: none`) | Jurkat rule, same k map, default, apex check; the rest stays a monomer of the same ingredient |
| R2 AF homodimer share | dimer copies = free x a_2; no SEC row: 0.3 | AF dimer fallback: all copies as dimers | Jurkat rule (a_2, default 0.3); the ipTM / ipSAE gates of the yeast fallback unchanged |
| R3 heterodimer evidence | pair kept if SEC co-elution r >= 0.6 OR STRING >= 700 | no evidence gate (IntAct / Complex Portal only set the bound fraction) | Jurkat gate (SEC r of `secms_pairs`, candidate pairs = the yeast AFDB heterodimer table; STRING v12 4932 combined score); the yeast filters and the proportional allocation unchanged |
| complex copies | Complex Portal stoichiometry cap | same | unchanged; `measured_free_pool` OFF |
| curated machines, legacy oligomers, plugin steps | - | yeast curated | unchanged (not SEC-fed) |

a_k = combined `yeast_final` assembled fraction at apex >= k x monomer mass (`secms_protein.tsv`, condition = replicate
= combined). Implementation: a new assembly preset (legacy_x5 + the three rules, pools: ingredients), selected in
`config/cells/yeast.yaml`; the reference configs `yeast_x5_core` / `yeast_x5_port` and preset legacy_x5 unchanged.
R1/R2 split an ingredient: the oligomer part becomes a new ingredient row with the parent's comp_weights, the remainder
keeps the parent row (monomer).

## 2. Metrics (inventory-only dry run into `wholecell_builds/yeast/_dev_unify/step3_secms_assembly/`)

Reference = the same code and inputs with preset legacy_x5 (equivalently the assembly tables of build
`yeast/20261119-87a2018-s0-step2ds`; the dry run says which). Total protein mass M = sum over the protein table of
copies x mass.

- **Model assembled share per protein** m_p = 1 - (copies of p in monomer ingredients, i.e. kind protein and not in
  oligo_map) / (copies of p in the protein table), clipped to [0, 1]; SEC a_p = yeast_final a_2 (and a_4).
- **Agreement**: (i) per-protein Spearman(m_p, a_p) over proteins with a yeast_final row and copies > 0, at k = 2 and
  k = 4; (ii) mass-weighted mean absolute difference MAE_w = sum w_p |m_p - a_p| / sum w_p (w = copies x mass), k = 2
  and 4; (iii) aggregate gap: mass-weighted model share vs SEC share (state of 2026-11-17: model 58.3 % vs SEC 83.6 %
  k2 / 63.0 % k4, r 0.07 - recomputed for the reference with the definition above).
- **Mass moved** per rule (R1, R2, R3): mass that changes assembly state (oligomer/dimer/heterodimer <-> monomer), pg
  and % of M; top 20 per-protein changes.
- **comp_weights**: every ingredient row of the reference present in the new table has an identical comp_weights
  string; every new (split) row carries its parent's string; `inventory/v0` (per-protein localisation) identical.
- **Yeast §3d benchmark** (plan §3d, B5 L4): SWAT_native intended tophit, hyperLOPIT intended tophit, intended image
  N/C test Spearman (`scripts/b5_localisation_check.py` yeast_metrics on the dry run's per-protein compartments).

## 3. Pass rules

| # | check | pass |
|---|---|---|
| S1 | calibration QC of every yeast run (anchor set used) | all pass (n >= 8, rho <= -0.7, residual <= x2.5) |
| S2 | comp_weights / localisation | identical as defined in §2 (no exception) |
| S3 | §3d benchmark | equal to the reference values (and >= 0.655 / 0.838 / 0.80) |
| S4 | agreement with SEC | Spearman(k2) AND MAE_w(k2) both better than the reference; k4 and the aggregate gaps reported |
| S5 | mass moved per rule | any rule moving > 2 % of M (0.07 pg of ~3.6 pg) is **reviewed**: listed with its top proteins in the report (verdict NEEDS REVIEW, not a failure) |
| S6 | unchanged parts | Complex Portal / curated machine copies (complexes_added.csv) identical; legacy_x5, yeast_x5_core, yeast_x5_port option sets identical before/after (dry-run options check) |
| S7 | pack, quick subset (Step-2 harness, `--subset quick`, `--reuse-auto` from the latest yeast step2 run, <= 8 cores), only if it fits in ~30 min | placed fraction >= reference - 0.005, clashes not increased; skipped = reported as not run |

Verdict: **PASS** (S1-S4, S6 pass, S7 pass or not run, no S5 review), **NEEDS REVIEW** (gates pass, >= 1 rule over
2 %), **FAIL** (any of S1-S4, S6 fails, or S7 run and failed). Merge into main only on PASS / reviewed NEEDS REVIEW and
a separate go-ahead.
