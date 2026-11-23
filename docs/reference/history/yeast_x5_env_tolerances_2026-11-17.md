# Yeast X5 legacy checks in the core environment: accepted numerical tolerances (2026-11-17)

Date: 2026-11-17. **Decided by the main session, Pedro informed (reversible).** A separate addendum, like
`yeast_x5_criteria_addendum_2026-11-17.md`: the frozen criteria `yeast_x5_reproduction_criteria.md` are unchanged.

## Decision

One core environment for all cells (`envs/`, conda env `vcell`: Python 3.12, numpy 2.5.4, numba 0.68.1, scipy 1.18.1,
+ xlrd 2.0.2). The yeast X5 port builds ran in the vendored yeast environment (`conda_envs/yeast_x5`: Python 3.11,
numpy 2.2.6, numba 0.61.2, scipy 1.15). Core steps re-expressing yeast stages (unification Step 1, items A1-A9) are
checked against the yeast port build s0 **in the core environment, at tolerance level**. Exactness in the yeast
environment is reported where it was shown, but it is not required.

## Accepted environment effects (measured 2026-11-17, unification Step 1)

| Effect | Where seen | Size | Accepted tolerance |
|---|---|---|---|
| numba code generation (0.61.2 vs 0.68.1): float operations reordered | geometry evaluator, ER membrane SDF (A1); the vendored cell3d shifts by the same amount between the two numba versions | <= 3.8e-12 nm | membrane SDF \|Δ\| <= 1e-9 nm; labels identical except points within 1e-9 nm of a boundary |
| numpy multinomial / random draws (2.2.6 vs 2.5.4) from the same state | tile request split P3 (A2): 236 / 7,364 (type, region) totals | ±1 copy | per (type, region) total within max(3√n, 0.5 % n) (the frozen P3 tolerance); everything downstream of a changed draw is judged at the frozen level T / seed-spread rules |
| Python 3.12 float `sum()` (compensated) vs 3.11 | assembly tables (A4): 4 values | 1 ULP | relative 1e-12 (ported core code uses left-to-right sums, so the core tables are identical) |
| scipy (1.15 vs 1.18) | `Rotation.from_matrix` rejects improper matrices in 1.18 (vendored amphiphysin.py): exposed a real chirality bug, fixed as PORT PATCH 4 (changes output; see `plugins/yeast/engine/PORT_PATCHES.md`) | n/a | not a tolerance: the reference for the affected builders is "port s0 + PORT PATCH 4" |
| other scipy numerical differences | none measured yet | — | relative 1e-9 on floats, identical integers; anything larger is reported, not accepted silently |

## Consequence for the reference

- Builders and stages unaffected by PORT PATCH 4: reference = port build s0 (`yeast/20261116-aabe7b5-s0-x5repro`).
- Serial pre-placed builders from amphiphysin (builder 19 of 28) onward, and everything downstream: reference =
  port s0 + PORT PATCH 4, made in the yeast environment and compared with the core environment at the tolerances above.
