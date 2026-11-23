# Freeze v8 — yeast acceptance addendum (2026-11-20)

Decision: Pedro Beltrao, 2026-11-20 (post hoc; made after seeing the verdict). Frozen criteria `freeze_v8.md`
(md5 9bbd657167cf754afe6e725d3cd95d24) are unchanged; this file records the decision separately.

Build: `data/derived/wholecell_builds/yeast/20261120-0a3ba79-s0-paperds` (commit 0a3ba79; downstream reuse of the packed
donor `yeast/20261119-bedcea9-s0-paper`). Automated verdict: NOT ACCEPTED (`validate/compare_freeze_v8/verdict.json`).
Accepted by decision, with the following:

1. **G0 frozen inputs** — `inv_yeast_complex_portal` changed vs the reference. Expected: SEC-MS wiring (step3-secms) uses
   the yeast Complex Portal set (taxon 559292). Added to the §2c allowances.
2. **G4 E2 other regions** — accepted deviations: LD_peripheral 0.831 vs 0.838, VAC_peripheral 0.988 vs 0.993
   (tolerance 0.003). Small fill losses on two minor membrane surfaces, consistent with the SEC-driven assembly changes.
3. **G5 E3 yeast B5 overrides (§3d)** — accepted deviation: SWAT native top-hit 0.66297 vs bound 0.6632 (miss 0.0002).
   hyperLOPIT (0.846 >= 0.8433) and image N/C pass.
4. **G6 E5 yeast SEC-MS** — not evaluable: bounds were never written into the criteria. Recorded as an omission; the SEC
   wiring was checked by its own dry runs and merge-plan tests (step3-secms, yeast-anchors-v3).

Release: tag `paper-yeast-u1` on 0a3ba79, `export/MANIFEST.sha256`, `FREEZE.json`, build folder read-only.
