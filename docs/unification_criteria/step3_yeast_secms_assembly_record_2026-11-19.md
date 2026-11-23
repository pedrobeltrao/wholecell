# Step 3 yeast SEC-MS assembly: dry-run record (2026-11-19)

Criteria: `step3_yeast_secms_assembly.md` (md5 e3153ceb96f7d79b09536bbdfdd256b2, commit 2fceeba). Dry run:
`scripts/step3_yeast_secms_dryrun.py`, out `wholecell_builds/yeast/_dev_unify/step3_secms_assembly/dryrun/`
(report.json, mass_moved.tsv, top20_changes.tsv, heterodimer_changes.tsv). Reference build
`yeast/20261119-87a2018-s0-step2ds`; the legacy_x5 pseudo-build reproduces its assembly tables, ingredients.tsv and
pack_table.tsv byte for byte. §3d reference sets: `_dev_unify/a8p2/20261116-aabe7b5-s0-x5repro/adapter/loc_inputs`.

| check | result |
|---|---|
| S1 calibration (anchor set v3) | pass: NaCl n 19 rho -0.852 x1.51; Pang 17 / -0.816 / x1.68; Swamy 20 / -0.966 / x1.37, 19 / -0.962 / x1.48, 20 / -0.957 / x1.37, 20 / -0.826 / x1.91; YPD = NaCl |
| S2 comp_weights | pass: 5,506 common ingredient rows identical (assembly, ingredients.tsv, pack_table); 529 `__homo` rows = parent's string; v0 identical |
| S3 §3d (intended) | pass, unchanged: SWAT_native 0.6647, hyperLOPIT 0.8483, image N/C test 0.8332 |
| S4 agreement (k2) | pass: Spearman 0.032 -> 0.092, MAE_w 0.357 -> 0.337 (k4: -0.003 -> 0.066, 0.405 -> 0.380). Aggregate: model 62.1 % -> 58.5 % vs SEC 82.6 % (k2) / 61.9 % (k4) |
| S5 mass moved (M = 3.605 pg) | R1 PDB homo share 0.041 pg (1.1 %, 110 proteins, 21 apex-check monomers incl. ACC1, DPS1); **R2 AF homodimer share 0.086 pg (2.4 %, 370 proteins: review)**; R3 heterodimer evidence 0.007 pg (0.2 %, 484 -> 443 pairs); net 0.133 pg (3.7 %) to free |
| S6 | pass: complexes_added identical; legacy_x5 / yeast_x5_core / yeast_x5_port option sets unchanged |
| S7 quick subset | not run (a 58-tile subset run took ~2.5 h; coordinator: skip if > ~30 min) |

Verdict: **NEEDS REVIEW** (R2 > 2 %: AF-dimer metabolic enzymes go from 100 % dimer to a_2, e.g. ILV5 0.50, FBA1 0.84,
PGI1 0.59, RNR1 0.52). Not merged (coordinator 2026-11-19: branch only).
