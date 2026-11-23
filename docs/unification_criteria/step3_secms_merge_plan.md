# Step 3 SEC-MS: merge plan (branches `step3-secms`, then `yeast-anchors-v3`)

Status: **PREPARED 2026-11-19 evening, not merged.** Pedro agreed the content (PLANNING_NOTES DECISION 2026-11-19:
Jurkat SEC-MS wiring, anchor set v3, HSPD1 excluded from the SEC apex check); he OKs the merges with the next step.
Nothing below has been run on main. No builds, no queue jobs.

## 1. What is merged

### 1a. `step3-secms` (base a12bdac; first)

| Commit | Content |
|---|---|
| caee4f2 | Jurkat + yeast SEC-MS wiring: `secms` inventory step takes an optional `combine` spec; jurkat.yaml = `final2026_jurkat` (lab NT + TCRAct24h + Fossati 2021, preset default, anchor set v2, QC gate) + `final2026_jurkat_combine` (label `jurkat_final`), assembly rules read `jurkat_final` (preset legacy_v7, `complex_copies: stoichiometry_cap`, measured_free_pool OFF); jurkat_v7_port.yaml keeps the v7 run; yeast.yaml SEC-MS tables only; yeast_x5_core `adapters.secms: null`. Same content as main 24efb2c, which was reverted on main (ef09448) |
| this branch, "HSPD1 apex_check_exclude" (the commit that adds this file) | generic rule `assembly_steps.data.apex_check_exclude` (list of {uniprot, gene, reason, share_source_id}; `wholecell/inventory/assembly/steps.py::apex_check_exclude`, declared optional data of `pdb_homo_oligomers`): a listed protein skips the SEC apex check and keeps its PDB state; with `share_source_id` its oligomer share is a_k of that SEC-MS run. Jurkat lists HSPD1 (P10809) with `share_source_id: fossati2021_jurkat` (reason in jurkat.yaml); jurkat_v7_port lists nothing (`[]`). Test `tests/test_assembly_steps.py::test_apex_check_exclude`. Presets, `ASSEMBLY_DEFAULTS` and the yeast ingredient path are untouched |

Why HSPD1: the HSP60 heptamer (PDB 7azp, 6/8 entries; 7 x 61.1 = 427 kDa) elutes at 147 kDa (lab NT) and 192 kDa
(lab TCRAct24h) — dissociated on the column — but at 480 kDa in Fossati 2021. The combined run puts its apex at 238 kDa
and a_4 = 0.466 (lab 0.205 / 0.353, Fossati 0.842), which would free 1.86 M HSPD1 copies (−0.188 pg, the largest
single mover of Step 3). With the exclusion HSPD1 is exactly as in the Step-2 / v7 state (row of
`homo_oligomers_v0.tsv` identical).

### 1b. `yeast-anchors-v3` (on top of caee4f2; second, only after its fix is OK'd)

| Commit | Content |
|---|---|
| b7f49dc | **SEC-MS yeast calibration anchor set v3** (`toolkits/secms/data/anchors_native_mass_yeast_v3.tsv`, review `ANCHORS_yeast_v3_2026-11-19.md`): v2 minus PDC1, EFT1 (ADH1, CDC19, 20S, CCT flagged, kept). Agreed by Pedro 2026-11-19 |
| 2fceeba | criteria `step3_yeast_secms_assembly.md` (frozen before implementation, md5 e3153ceb…) |
| 41d3ed5 | preset `x5_secms` (= legacy_x5 + the Jurkat SEC-fed rules), yeast.yaml reads `yeast_final`; record `step3_yeast_secms_assembly_record_2026-11-19.md` (verdict NEEDS REVIEW: R2 AF-homodimer share moves 2.4 % > 2 %) |
| (pending) | the fix Pedro asked for (why the assembled share moves away from SEC, 62.1 → 58.5 % vs 82.6 % k2); merged only when OK'd. Owned by the yeast SEC-MS agent — not prepared here |

`git merge-tree` (2026-11-19 evening): `step3-secms` (with the HSPD1 commit) and `yeast-anchors-v3` (41d3ed5) both
merge into main 12a1e7a without conflicts, and `yeast-anchors-v3` merges cleanly on top of the merged `step3-secms`
(the HSPD1 test is placed away from the yeast branch's new test; `apex_check_exclude` is cell data, not a preset key,
so `test_yeast_x5_secms`'s `share == legacy_v7 share` assertion is unaffected). Re-run `git merge-tree` before merging
if either branch moves.

## 2. Order and commands (main repository; Pedro's OK first)

```
cd wholecell                                   # main checkout, clean
git merge --no-ff step3-secms                  # 1. Jurkat/yeast SEC-MS wiring + HSPD1 exclusion
#   -> run §3a checks; only then:
git merge --no-ff yeast-anchors-v3             # 2. anchor set v3 + yeast SEC-driven assembly (after its fix is OK'd)
#   -> run §3b checks
git merge --no-ff freeze-prep                  # 3. freeze criteria (its E4 bound assumes step 1; §2a status lines)
```
Commit messages end with the attribution lines of the session. If `yeast-anchors-v3` is not OK'd, only b7f49dc
(anchors v3, data + README) can be cherry-picked alone; the yeast assembly stays legacy_x5 (yeast.yaml SEC-MS tables
only, as after step 1).

## 3. Checks after each merge

### 3a. After `step3-secms`

Tests (core env `conda_envs/vcell/bin/python`, `PYTHONHASHSEED=0`; run as scripts, pytest is not installed):
`tests/test_assembly_steps.py` (7 incl. `test_apex_check_exclude`), `tests/test_secms.py` (8),
`tests/test_secms_lab_formats.py` (7), `tests/test_inventory.py` (3), `tests/test_copies.py`, `tests/test_localisation.py`,
`tests/test_validate_suite.py`; testcell CI (`scripts/ci_testcell.sh`, PASS; vs previous: identical expected — the
testcell has no `apex_check_exclude`). All of the first four passed on the branch (2026-11-19 evening).

Reference-config dry runs (inventory only, into `_dev_unify/step3_merge/`, never into a build dir):
1. **jurkat_v7_port** inventory → `scripts/secms_regression_check.py --build <dir>`: verdict REPRODUCED (S1–S11; v7
   assembly tables byte-identical; `apex_check_exclude: []`).
2. **yeast_x5_core**: `AS.options` = legacy_x5 unchanged (test), `adapters.secms: null`; assembly tables unchanged vs
   the A9 reference (`scripts/assembly_x5_check.py` only if the yeast engine code changed — it does not in this merge).
3. **yeast.yaml** inventory: SEC-MS tables written (`secms/`, label `yeast_final`), assembly still legacy_x5 (tables
   identical to the Step-2 yeast inventory).
4. **Jurkat inventory dry run** (jurkat.yaml; geometry copied from the Step-2 build; reference = Step-2 inventory
   `jurkat/_dev_unify/step3_secms/baseline_main_a12bdac`) must reproduce
   `jurkat/_dev_unify/step3_secms_hspd1/dryrun.json` (md5 db02d46fe84dcb3c4e8c2adcf85a3e8d):

| Quantity | Step-2 baseline | wiring caee4f2 alone | **with HSPD1 exclusion (to reproduce)** |
|---|---|---|---|
| complexed share (CP + homo-oligomer + AFDB homo/hetero mass / 101.086 pg) | 0.52529 | 0.53523 | **0.53709** |
| homo-oligomer mass (pg); rows | 6.1307; 812 | 6.6872; 751 | **6.8754; 751** |
| AFDB homodimers (pg); rows | 6.2090; 1,162 | 6.7024; 1,198 | **6.7024; 1,198** |
| AFDB heterodimers (pg); rows | 4.5405; 1,402 | 4.4947; 1,436 | **4.4947; 1,436** |
| Complex Portal mass (pg) | 36.21987 | 36.21987 | **36.21987** (0 complexes changed) |
| proteins changed vs baseline; mass moved (pg) | – | 1,072; 2.721 | **1,071; 2.533** (1.863 to assembled, 0.670 to free) |
| HSPD1 heptamers (7azp): share, n, copies | 0.8418, 594,934, 4,164,538 | 0.4664, 329,658, 2,307,606 | **0.8418, 594,934, 4,164,538** |
| pack_table species: n; copies changed vs baseline | 16,249 | 16,146; 1,922 | **16,146; 1,920** |

Versus the wiring alone, the exclusion changes exactly one protein (HSPD1: +1,856,932 copies in heptamers, +0.188 pg
assembled) and two pack-table rows (HSPD1_x7 329,658 → 594,934; HSPD1 free 2,639,686 → 782,754); copies.tsv,
inventory_v0, complexes_v0, sec_assembly, AFDB homo/heterodimer tables identical; homo_oligomers_v0 identical except
the HSPD1 row, which equals the Step-2 row. Top movers vs baseline are now IMPDH2 +0.113, PRDX1 +0.092, KHDRBS1 +0.086,
GAPDH −0.071, TOP2A −0.058 pg. Run: `step3_secms_hspd1/{dryrun_hspd1.sh, run_inv4.py, hspd1_cmp.py, new.log}` (284 s,
4 workers). Tolerance for the post-merge rerun: tables identical (same code + frozen inputs); else explain.

### 3b. After `yeast-anchors-v3`

Tests: as 3a + `tests/test_assembly_steps.py::test_yeast_x5_secms`, `tests/test_secms_lab_formats.py` (anchor sets
incl. v3). Dry runs: `scripts/step3_yeast_secms_dryrun.py` reproduces the record of the fixed branch (S1–S6;
currently `step3_yeast_secms_assembly_record_2026-11-19.md`, to be replaced by the fix's record), yeast_x5_core and
the legacy_x5 pseudo-build byte-identical, and 3a item 4 (Jurkat) unchanged (the yeast commit touches
`wholecell/inventory/secms.py::_frozen_paths` and the ingredient path only).

## 4. After the merges
- freeze criteria (`freeze-prep`): G1 `secms_jurkat` satisfied by config; E4 bound = the 3a numbers (share 0.5371 ±
  0.001, moved 2.53 ± 0.10 pg); `secms_yeast` satisfied once the yeast commit is in.
- Re-run the B5 yeast dry run on the merged state (unification_progress "next steps").
