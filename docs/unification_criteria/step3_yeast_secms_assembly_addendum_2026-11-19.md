# Step 3 yeast SEC-MS assembly: addendum (2026-11-19, evening)

Addendum to `step3_yeast_secms_assembly.md` (md5 e3153ceb96f7d79b09536bbdfdd256b2, frozen, unchanged). Trigger: Pedro
(2026-11-19) - the first dry run (`dryrun/`, commit 41d3ed5) improved per-protein agreement but moved the aggregate model
assembled share away from SEC (62.1 % -> 58.5 % vs SEC 82.6 % k2 / 61.9 % k4). Pass rules S1-S7 unchanged.

## 1. Diagnosis (dry run `dryrun/`, k2-covered mass 3.582 pg, % points of that mass)

| cause | contribution to the -3.63 pt aggregate change |
|---|---|
| R2 AF homodimer share (legacy: 100 % dimer; now a_2) | -2.39 |
| R1 apex check: 21 PDB oligomers -> monomer (ACC1, DPS1, ...) | -0.73 |
| R1 PDB share (a_k by state) | -0.42 |
| R3 heterodimer evidence gate (484 -> 443 pairs) | -0.09 |

- (ii) asymmetry - the main cause. legacy_x5 puts 100 % of every oligomer-annotated protein into the oligomer, so a share
  a_k <= 1 can only lower it, and R3 can only remove pairs; no rule raises a protein. In Jurkat the same rules also raise
  (legacy_v7 baseline shares are constants < 1: +1.86 / -0.86 pg). The rules touch 16.6 % of the mass, where the model
  was ABOVE SEC (95.3 % vs SEC a_2 78.6 %). The aggregate gap is elsewhere: proteins the model keeps fully monomeric (34 %
  of mass) have SEC a_2 0.75 (median apex 31 x monomer: EF-1A/EF3/Hsp70/PAB1/CHC1 etc. in assemblies the model has no
  structure for), 25.1 of the 28.1 pt under-gap. Only ~6 pts of it carries any structural evidence (rejected AFDB
  heterodimer candidates 4.9, UniProt n > 1 without a structure 1.6, rejected AF homodimers 0.3).
- (iv) apex check - wrong direction. All 21 demotions had the apex ABOVE 2.5 x the state (none below 0.4 x); their SEC a_2
  is 0.53-1.00 (mass-weighted 0.91), and 15 of 21 apexes are above the top calibrant (1,236 kDa; Pang only to 268 kDa), i.e.
  extrapolated. An apex above the oligomer is evidence of a larger assembly, not of a monomer.
- (i) definitions. Model m_p = 1 - monomer-ingredient copies / protein copies: complexes (29.5 pts), heterodimers (5.4),
  homo-oligomers (27.1 -> 23.6) all count as assembled, copies not placed (surplus) also count (0.006 pg: negligible).
  SEC a_k = profile intensity at apparent mass >= k x monomer. Same definition on the model side (copies in an assembly of
  mass >= k x monomer; unplaced copies left out): legacy 58.2 % / new 55.5 % at k2 (SEC 82.6 %), 37.2 % / 36.9 % at k4 (SEC
  61.9 %); the gap is larger, not smaller, when computed the same way. The SEC aggregate itself is not tight: k2 band at
  -/+ one calibration SD = 67.2-95.1 %. For dimers the k = 2 cut sits on the state's own mass: a dimer whose apex is at
  ~2 x reads a_2 ~ 0.5 (ILV5 apex 1.93 x, a_2 0.50; RNR1 2.09 x, 0.51; ILV3 1.95 x, 0.44; the monomer PGK1 1.84 x reads
  0.40). This is the k_by_state map of legacy_v7 (n <= 3 -> k 2) and is kept (a cut below 2 x needs a new toolkit
  threshold, a shared default; not done).
- (iii) no sign / denominator error: a_k is a copy (intensity) fraction, model m_p a copy fraction, both mass-weighted by
  copies x mass. The share is applied to the residual pool while a_k is a whole-protein fraction (76 R1/R2 proteins with
  copies in complexes / heterodimers): a consistent denominator would move -0.1 pt; not changed.

## 2. Rule changes (x5_secms only; legacy_v7 / Jurkat and legacy_x5 unchanged)

- **A1** `pdb_homo_oligomers.share.apex_check_side: below`: the SEC apex check demotes a PDB oligomer to monomer only when
  the apex is BELOW apex_ratio[0] x (state x monomer); an apex above the window keeps the oligomer at its SEC share
  (decision suffix `SEC apex above state: kept`). legacy_v7 default `both`.
- **A2** `pdb_homo_oligomers.afdb_fallback.ipsae_only_min: 0.6`: the AF homodimer fallback also accepts max ipSAE > 0.6
  alone = the legacy_v7 `afdb_homodimers.min_ipsae` gate (the symmetric, Jurkat-semantics rule: structural evidence of a
  dimer gets its SEC share a_2). Replaces "ipTM / ipSAE gates of the yeast fallback unchanged" of the criteria §1.

Reported in addition (dry-run report.json `aggregate`, not pass rules): change by rule, rule-set model vs SEC, SEC k2
calibration band, same-definition shares (§1).

## 3. Dry run `dryrun_fix1/` (commit 4e6c921)

| | legacy_x5 | x5_secms (41d3ed5) | x5_secms + A1, A2 |
|---|---|---|---|
| aggregate model k2 / k4 (SEC 82.6 / 61.9 %) | 62.1 % | 58.5 % | 59.4 % |
| rule set (16.6 % of mass) model vs SEC a_2 | 95.3 vs 78.6 % | 73.3 vs 78.7 % | 79.1 vs 78.6 % |
| Spearman k2 / k4 | 0.032 / -0.003 | 0.092 / 0.066 | 0.096 / 0.070 |
| MAE_w k2 / k4 | 0.357 / 0.405 | 0.337 / 0.380 | 0.328 / 0.373 |
| change by rule (pts) | - | R2 -2.39, apex -0.73, R1 -0.42, R3 -0.09 | R2 -2.39, A2 +0.28, R1 -0.42, apex-kept -0.08, R3 -0.09 = -2.70 |
| mass moved R1 / R2 / R3 | - | 1.1 / 2.4 / 0.2 % | 0.5 / 2.7 / 0.2 % (R2 review) |

S1-S4, S6 pass; S2: 5,497 common rows identical comp_weights, 581 `__homo` rows = parent, v0 identical; S3 SWAT 0.6647,
hyperLOPIT 0.8483, image N/C 0.8332 (unchanged); legacy pseudo-build = reference build tables. Jurkat inventory-only run
(`jurkat/_dev_unify/step3_secms_yfix/new_4e6c921`): pack_table 4def9ebd = step3-secms state (`step3_secms/new_a12bdac_wiring`),
58 / 58 inventory tables byte-identical.

Conclusion: on the proteins the SEC-fed rules can act on, the model now equals SEC (79.1 vs 78.6 %). The aggregate stays
below legacy because the legacy 100 % assembly of those proteins was above SEC and the rest of the gap (~25 pts) sits in
proteins with no assembled structure in the model; no SEC-share rule can close it. Verdict: NEEDS REVIEW (R2 2.7 %).
