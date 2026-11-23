# B1 packer: second addendum 2026-11-18 (outcome: per-cell settings)

Date: 2026-11-18. **Decided by the main session, Pedro informed.** Post-hoc, after the reduced benchmark
(`B1_packer_addendum_2026-11-18.md`, md5 59e97f2e…; results `wholecell_builds/_benchmarks/b1r_e890bf4/`). The frozen
criteria `B1_packer.md` (md5 f544eb73…) and the first addendum are not edited.

## Finding

No single candidate setting passed in both cells on the frozen busy-tile Step-2 subsets. The selection rule's
"placed ≥ 0.98 by mass in every volume region ≥ 1 µm³" is not met even by the v7 settings on the Jurkat subset
(0.887 overall, ER lumen 0.83), whose tiles are deliberately the crowded ones; the same v7 settings reach 0.964 by mass
over the whole Jurkat cell. Packing shapes were rejected for Jurkat (26,177 packed–packed clashes per 1e6 against
library trees, lower fill in every region, slower); the squeeze arm timed out in the thinning pilot on 500 nm tiles.
Yeast C3 (legacy_x5 + rounds 6, pilot thinning, shell cap, targeted squeeze 0.97, packing shapes, 2000/1200 try
budgets) was the only eligible run: 0.971 overall, worst region nucleoplasm 0.979, seam 1.01, 0 clashes, 1,503 s =
0.47× the legacy_x5 wall.

## Decision (outcome of B1)

Per-cell budgets instead of one common setting:
- `wholecell/pack/switches.py` `COMMON_DEFAULTS = {}`: the v7 values; Jurkat unchanged.
- Yeast pins C3 (`config/cells/yeast.yaml` `pack.preset: [legacy_x5, b1_s4y]`): the same packer algorithm with
  different budgets.
- Reference configs keep their legacy settings (`yeast_x5_core.yaml` pins `legacy_x5`; `jurkat_v7_port.yaml` pins
  `legacy_v7`).
The build-level B1 checks (K1–K9) are judged in the combined Step-2 builds.
