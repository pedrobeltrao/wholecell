# B5 localisation addendum: record 2026-11-19 (combined Step-2 build)

Date: 2026-11-19. Decision source: **main session 2026-11-19, Pedro informed**. Separate dated record; the frozen
`B5_localisation_addendum_2026-11-18.md` (md5 4bfe8e52df0fba0ecdc7f94c938c3719) is not edited.

On the Jurkat combined Step-2 build (`jurkat/20261118-40d7967-s20261115-step2`, re-score
`validate/compare_step2_v2/`) two OpenCell C1 rows of the combined no_regression fell outside ±0.01 that the B5
revision had not listed. Both are classified **EXPECTED consequences of B5**:

| Suite row | Build | Reference (B3 build) | Change |
|---|---|---|---|
| C1 `opencell.placed.tophit.unweighted` | 0.8804 | 0.8943 | −0.0139 |
| C1 `opencell.intended.tophit.mass_weighted` | 0.9539 | 0.9641 | −0.0102 |

Why: tophit was never a B5 gate (addendum 2026-11-18 gates `agree` on independent references); the gated OpenCell
`agree` passes (intended 0.8375 vs 0.8314, placed 0.8369 vs 0.8308); the intended unweighted tophit shift (−0.014,
0.8858 vs 0.8997) was already in the B5 dry run (revision §2 L2) and these two rows are its placed / mass-weighted
variants.

Applied via `step2_combined_accepted_deviations_2026-11-19.json` (class EXPECTED, step B5).
