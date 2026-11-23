# Step-2 subsets: addendum 2026-11-18 (harness validation, one mass table, S-P5_frac threshold)

Addendum to `step2_subsets.md` (frozen, md5 `ff9f79baf2681eed42b790ef6e8943e3`; unchanged). Written after the
harness smoke run and the yeast pre-B4 run, BEFORE any result of tests a-d (Jurkat B4, yeast core / B6 / B4) existed.

## A1. Harness validation (§3 neighbour rule)

- Jurkat: commit bc121a6 + the B3 build's upstream, context = B3 build (`jurkat/_dev_unify/step2/j0_b3_repro_smoke`):
  the 56 subset tiles are **bitwise identical** to the B3 build's tiles (56/56). The rule reproduces a full build.
- Yeast: commit 9977542 + the upstream of the stopped A9 build (`yeast/_dev_unify/step2/yd_preb4_9977542`): the
  round-0 log of the subset run equals the A9 build's own round-0 log tile by tile where the A9 build got to (phase 0:
  12/12 tiles identical requested / placed; phases 1-2: 8/9 identical, one tile 488 vs 478 placed: the context is the
  port s0 there, not the A9 build). Rounds 1-5 cannot be checked (A9 stopped in round 0).

## A2. One mass table for both sides of a comparison

Yeast `libs.json` has no masses; the metrics use library bead masses, which come from each build's own library
(the port s0 library and a core library give different bead masses for some assembled species). From now on the
reference side of a comparison uses the TEST's masses by library name (`metrics/masses_run.json`; pre-placed entries
`preplaced:<lib>`). A measurement fix; no rule changes. Jurkat uses `libs.json` mass_da on both sides (unchanged).

## A3. S-P5_frac (core vs port, judged): threshold calibrated on full builds

The frozen rule (≥ 99 % of (type, region) pairs with ≥ 10 requested copies within max(0.01, 3·√(f(1−f)/n))) is not
reached by full builds: copies placed in one tile are correlated, so the binomial 3σ is too narrow. Measured on the
frozen yeast subset (pairs with ≥ 10 requested copies, mem_PEROX excluded where the pre-fix A7 core build is used):

| Pair of packs | share within |
|---|---|
| A7 full core build (7a70243) vs port s0 | 0.9454 |
| subset run 9977542 vs A7 full core build | 0.9562 |
| subset run 9977542 vs port s0 | 0.9472 |

New threshold for S-P5_frac: share within ≥ **0.935** (= full core build vs port − 0.01). Everything else in §5 is
unchanged.

## A4. Finding recorded with the calibration (not a rule)

Nucleolus, yeast: the core at 9977542 places every large pre-ribosome (preribo_60S_nucleolar, preribo_90S) of the
fully nucleolar tile 1679 in round 0 (532/532 requested; the A9 full build's own round-0 log shows the same 532/532),
the port s0 places 85/110 and the pre-fix core A7 94/130 of these particles there. Hence on the subset the nucleolus
placed MASS fraction is 0.934 (core) vs 0.868 (s0) while the copies fraction is 0.905 vs 0.918, and S-C11 nucleolus
differs by +8.7 mg/ml. This is a core-vs-port difference of the packer at/after the kernel membrane-code fix, not a
harness effect.
