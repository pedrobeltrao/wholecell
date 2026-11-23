# Freeze v8 — Jurkat acceptance addendum (2026-11-22)

Decision: Pedro Beltrao, 2026-11-22 (post hoc; made after seeing the verdict). Frozen criteria `freeze_v8.md`
(md5 9bbd657167cf754afe6e725d3cd95d24) are unchanged; this file records the decision separately.

Build: `data/derived/wholecell_builds/jurkat/20261120-0a3ba79-s20261115-paper` (commit 0a3ba79; complete, 20/20 steps;
resumed on eu-g4-019 after the 2026-11-20 node move). Automated verdict: NOT ACCEPTED
(`validate/compare_freeze_v8/verdict.json`, 2026-11-21 22:27). Accepted by decision, with the following:

1. **G2 B3 L1-L7 (L6 orientation vs v7)** — 13 flipped TM species (tilt 180°). They are exactly the 13 glyco-captured
   units of the PM topology fix (E7 / D7: CD7_x2, CPX-1794, CPX-1802, CPX-1803, CPX-26677, CPX-937, CPX-9541,
   LNPEP_afx2, NECTIN2-NECTIN1, PVR-CD226, SLC44A1_x2, SLC4A7, TFRC_afx2) = a listed cause; the check script did not
   read the E7 list. Bookkeeping. (L5 is scored under the 2026-11-17 addendum: L5a passes, 21 = v7.)
2. **G7 E7 glyco units not flagged** — 0 of the 13 flagged, rotations exact, only fix units changed; the item fails
   because 3 units (CD7_x2, NECTIN2-NECTIN1, SLC44A1_x2) have no placed PM copies to score. Bookkeeping.
3. **G4 E2 nucleolus placed fraction** — accepted deviation: 0.844 vs bound 0.86 (reference 0.807; improved by
   +0.037 but short of the bound set from the B3 / 43-tile checks).
4. **G4 E2 other regions / G2 B1.K2** — accepted deviations (placed mass fraction vs reference − 0.003): speckle
   0.920 vs 0.934, pericentriolar material 0.838 vs 0.880, Golgi cytosolic shell 0.817 vs 0.822, Golgi lumen 0.983 vs
   0.989, endosome lumen 0.926 vs 0.953; cytosol 0.972 vs 0.977. Small fill losses in minor regions; overall placed
   97.7 % (108.9 M of 111.4 M), no overflow. Not diagnosed further (decision).
5. **no_regression C1** — accepted deviation: HPA-JURKAT hold-out top-hit 0.596 vs 0.615 (tolerance 0.01).

Downstream: the v8 simulations (applications/diffusion: BD v2 verification PASS; applications/synapse: S2.0 x CRISPR)
were run on this build before release by decision (2026-11-21) and pin its fingerprint; the build content is unchanged
by the release.

Release: tag `paper-jurkat-v8` on 0a3ba79, `export/MANIFEST.sha256`, `FREEZE.json`, build folder read-only.
