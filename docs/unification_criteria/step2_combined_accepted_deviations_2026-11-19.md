# Step-2 combined build (Jurkat): accepted deviations, 2026-11-19

Date: 2026-11-19. Decision source: **main session 2026-11-19, Pedro informed** (C9 also main session 2026-11-18).
Build: `jurkat/20261118-40d7967-s20261115-step2`. Previous verdict: `validate/compare_step2_v2/verdict.json`
(NOT ACCEPTED; `scripts/step2_combined_check.py` at main 31c1ac7, i.e. after the (a) verdict-script fixes). Reference
build of the combined no_regression: `jurkat/20261117-bc121a6-s20261115-b3` (accepted B3 build).

Machine-readable twin: `step2_combined_accepted_deviations_2026-11-19.json`, read by
`scripts/step2_combined_check.py --accepted <json>`. An entry applies only to the build it lists and only while the
value is not worse than its `limit`; a failing item is ACCEPTED_DEVIATION only if every failing key of it is covered.
Verdict rule: ACCEPTED if every remaining failure is an accepted deviation. Frozen criteria files are not edited.

Classes: **KNOWN** = known defect, fix accepted and deferred to the freeze builds; **ACCEPTED** = accepted deviation of
the combined build (no fix planned); **EXPECTED** = consequence of an accepted step (B5), recorded in
`B5_localisation_addendum_record_2026-11-19.md`.

| # | Check / metric | Value | Reference | Bound | Class | Cause | Fix |
|---|---|---|---|---|---|---|---|
| 1 | B3.B3 C9 and suite C9 `sign_errors_per_1e6` | 1.0915 | 0.0 | ≤ ref + 1.0 | KNOWN | 3 orientation sign flips (2 ER, 1 endosome) at positions with a weak SDF gradient | deferred to the freeze builds: accepted orientation fix (decisions 2026-11-18 + 2026-11-19) |
| 2 | B1.K2 nucleolus placed mass fraction | 0.8074 | 0.8636 | ≥ 0.98 or ≥ ref − 0.003 | KNOWN | B2 field split: field-driven per-tile requests exceed packer capacity in the crowded nucleolus | deferred to the freeze builds: cap field-driven per-tile requests in crowded volume regions at the reference per-tile request |
| 3 | B1.K2 speckle placed mass fraction | 0.9340 | 0.9612 | same | ACCEPTED | B5 mass moved into the speckle; packer capacity | — |
| 4 | B1.K2 cytosol placed mass fraction | 0.9705 | 0.9771 | same | ACCEPTED | B5 mass moved into the cytosol; packer capacity | — |
| 5 | B1.K2 golgi_lumen placed mass fraction | 0.9755 | 0.9802 | same | ACCEPTED | B5 mass / packer capacity (−0.0047, just beyond −0.003) | — |
| 6 | suite C1 `opencell.placed.tophit.unweighted` | 0.8804 | 0.8943 | ≥ ref − 0.01 | EXPECTED (B5) | tophit never gated; placed level of the B5 intended tophit shift already in the B5 dry run | — |
| 7 | suite C1 `opencell.intended.tophit.mass_weighted` | 0.9539 | 0.9641 | ≥ ref − 0.01 | EXPECTED (B5) | mass-weighted variant of the same shift (−0.0102) | — |

Context for 6–7: the gated OpenCell `agree` passes (intended 0.8375 vs 0.8314, placed 0.8369 vs 0.8308); the intended
unweighted tophit (−0.014) was already EXPECTED (B5 revision §2 L2).

Main-session wording for 3–5: "expect ACCEPTED with documented deviations". The freeze builds are judged afresh (this
record lists only this build).

Re-score: `validate/compare_step2_v3/` (this record + `--accepted`).
