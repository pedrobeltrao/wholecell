# B1 packer: addendum 2026-11-18 (reduced benchmark on the frozen Step-2 subsets)

Addendum to `docs/unification_criteria/B1_packer.md` (frozen, md5 f544eb73724b2accdd9c2fd0bec88b89; not edited).
Written 2026-11-18 BEFORE any run of the reduced benchmark (Pedro 2026-11-18: the 500-tile Jurkat bar + whole yeast
cell benchmark of S1-S5 (~17 h) is cancelled; shrink it to the serious candidates on the frozen Step-2 subsets).
Subsets, neighbour rule, seeds and subset metric definitions are those of
`docs/unification_criteria/step2_subsets.md` (frozen, md5 ff9f79baf2681eed42b790ef6e8943e3) and
`config/subsets/{jurkat,yeast}_step2.tsv` (Jurkat 56, yeast 58 tiles); nothing of those is changed here. Where this
addendum and the frozen B1 file differ, this addendum governs the reduced benchmark only; Part B (build-level
acceptance K1-K9) and Part C (mechanism) of the frozen file stand unchanged and are applied at the combined build.

## 1. Evidence used to choose the candidates (existing runs only)

| Run | Settings | Result |
|---|---|---|
| testcell A7 phase 1, defaults (`jurkat/_dev_unify/a7/testcell_defaults_ef905b0`) | v7 (single visit, 200 / 25, library trees) | mass 0.979; nucleoplasm **0.79**, ER_lumen 0.911; pack 326 s |
| testcell A7 phase 1, switches on (`.../testcell_switches_on_ef905b0`) | rounds 6, thinning (+ pilot squeeze), shell cap, targeted squeeze 0.97, library trees | of request 0.999, of unthinned 0.978 (nucleoplasm × 0.782); nucleoplasm 0.995, ER_lumen 0.916; pack **5,405 s (16.6 × v7)** |
| testcell A7 phase 2, 2 busiest tiles (`testcell/_dev_unify/a7p2`) | as above, library trees vs `pack_shapes` | pack 673 s → 237 s (nucleus-tile squeeze 636 → 203 s); before the squeeze the medium size-class round of the nucleus tiles places only 96.3 % / 93.2 %; with packing shapes the pilot thins nucleoplasm more (× 0.62 vs × 0.78): nucleus tile placed copies 93,702 → 77,715 |
| Jurkat whole cell v7 (plan §2b) | v7 | placed 0.964 |
| yeast history (`plugins/yeast/engine/model_validation_and_tuning.md` §16) | single pass, 10 × budgets (builds N-Q) vs build S (6 rounds, 15-residue shapes, shell caps, thinning, squeeze) | 96.5-96.7 % → 99.3 % |
| yeast A7 whole cell through the core (`yeast/_dev_unify/a7y/full_7a70243`) | legacy_x5 (2000 / 1200, rounds 6, thinning, shell cap, full squeeze, pack_trees_15) | 0.9918 (s0 0.9924); X5 squeeze = 89 % of pack time (plan) |

Reading. (a) Single visit with v7 budgets misses 0.98 in the crowded volume regions of both cells (testcell nucleoplasm
0.79; yeast 96.5 %). (b) Rounds + thinning without a squeeze still leave crowded tiles below 0.98 (96.3 % / 93.2 % in
the medium class of the nucleus tiles after thinning) → frozen S2 and S3 are not serious candidates. (c) The targeted
squeeze reaches the bar, but with the library trees it costs 16.6 × v7 on the testcell, far beyond the 2 × limit →
frozen S4 (library trees) is not a serious candidate for Jurkat; packing shapes cut the squeeze ~3 ×. (d) Full squeeze
(S5) is excluded as a common default by Pedro (targeted squeeze); it stays the yeast reference (legacy_x5). (e) Risk
recorded before the run: bulkier packing trees make the pilot thin more, so condition 3 (placed of the unthinned
request ≥ S1 − 0.003) may fail for the shape-based settings in Jurkat.

## 2. Reduced design

Candidates (every algorithm key set explicitly, as the frozen A3; all other `pack:` keys are the cell's; yeast keeps
its pinned `legacy_x5` with the candidate merged after it, preset list `[legacy_x5, <setting>]`):

| Id | Definition | = frozen | Jurkat | Yeast |
|---|---|---|---|---|
| C1 v7 + shapes | rounds 1; 200 / 25, thin lumens × 10 / × 16; no thinning, no shell cap; squeeze none; `shapes: pack_shapes` | S1 + shapes (yeast: = S1, its shapes are already pack_trees_15) | yes | yes |
| C2 rounds + thinning + targeted squeeze + shapes | rounds 6 equal_volume alternate; thinning on (pilot_tiles 16, min_rel 1.1, safety 0.97, keep_above 0.98, reference cytosol, pilot_squeeze on); shell cap on; squeeze underfilled, threshold 0.97, rounds 4, iterations 60, insert_frac 0.15; 200 / 25; `shapes: pack_shapes` | S4p (yeast: = S4) | yes | yes |
| C3 = C2 with the yeast budgets | as C2, budgets 2000 / 1200 in every region (thin-lumen factors 1) | S4y + shapes (yeast: = S4y) | yes | yes |

Reference arms (not candidates; needed by the frozen conditions 3-5):

| Id | Cell | Definition | Used for |
|---|---|---|---|
| R-J = S1 | Jurkat | v7 settings, library trees (the frozen S1) | condition 3 (placed of unthinned), condition 4 (R clash / containment), condition 5 (wall), P8 shapes check |
| R-Y = legacy_x5 | yeast | the yeast build's own settings (preset `legacy_x5` alone) | condition 4 (R) and condition 5 (wall); yeast condition 3 uses C1 (= frozen S1) |

Runs: 8 (Jurkat R-J, C1, C2, C3; yeast C1, R-Y, C2, C3), in this order, one at a time, ≤ 24 pack workers, the same
node, through the sequential build queue (`_queue/jobs.tsv`, job `b1_reduced`, script `_queue/b1_reduced.sh`).
Seeds: the configs' own (Jurkat master 20261115, yeast master 0), the same for every run; no replicate seeds.

Inputs. Jurkat: the accepted B3 build `jurkat/20261117-bc121a6-s20261115-b3` (library reused; context tiles per the
Step-2 neighbour rule) with main's `config/cells/jurkat.yaml`, i.e. main's B4 pre-placed configuration (whole-cell
pre-placed stage, or the upstream of an existing Step-2 Jurkat B4 run of the same commit and config if one exists: the
script records which). Yeast: the core yeast pipeline of main (`config/cells/yeast.yaml`) on its subset, context = the
port s0 as frozen in Step-2 §3. Harness: `scripts/step2_subset_test.py` (Step-2 §2) if committed in time; otherwise
`scripts/b1_benchmark.py` adapted to pack only the Step-2 subset tiles with the Step-2 neighbour rule. The report states
which.

Budget: ≤ 5 h wall in total on ≤ 24 cores. Timeouts (a run killed at its timeout is ineligible): R-J 60 min; Jurkat
candidates 2.2 × wall(R-J) (at least 30 min); R-Y 90 min; yeast candidates 2.2 × wall(R-Y) (at least 30 min). A run
that cannot start before 4 h 30 min of the budget have elapsed is recorded as `not run (budget)` and is ineligible.
A failed run is recorded with its log; the benchmark continues.

## 3. Metrics and selection rule on the subset (frozen A4 / A5, made concrete for a subset)

Metrics are the frozen A4 ones restricted to the subset tiles (Step-2 §4 definitions, one code version for all runs):

1. **Regions.** Gated regions as frozen (volume regions ≥ 1 µm³ in the cell with requested copies: Jurkat cytosol,
   nucleoplasm, matrix, ER_lumen, golgi_lumen, endo_lumen, nucleolus, speckle; yeast cytosol, nucleoplasm, vacuole,
   cell_wall). A region's placed_mass_fraction on the subset = Σ placed × mass / Σ request × mass **over the subset
   tiles of that region** (request = the packer's `job_n`, after thinning). Condition 1: ≥ 0.98 −
   3·√(0.98·0.02 / n), n = requested copies of the region on the subset (Step-2 binomial term; for n ≥ 2,000 the
   allowance is < 0.01). A gated region with < 100 requested copies on the subset is reported, not gated.
2. **Seam.** |index − 1| ≤ 0.1 where measurable: Step-2 S-seam tiles (≥ 95 % one region), with ≥ 1,000 interior
   instances summed over them; Jurkat gated on nucleoplasm (the frozen Jurkat choice), cytosol reported; yeast gated on
   cytosol. Not measurable → condition reported `n/a` and not gated (deferred to K3).
3. **Unthinned.** placed_of_unthinned (subset, by mass) ≥ that of the cell's S1 − 0.003 (Jurkat R-J; yeast C1).
4. **Clashes / containment.** S-C7a packed-packed ≤ 1 per 10⁶ packed; clash rate and S-C8 ≤ R + max(10 % rel.,
   1 per 10⁶) (R = R-J / R-Y). Jurkat clashes are measured on the library trees (frozen A4: the shapes the viewer and
   simulations use); the packing-tree value is reported.
5. **Wall time.** Subset pack wall (rounds + squeeze, the harness's subset-pack timer) ≤ 2 × that of the reference arm
   (Jurkat R-J, yeast R-Y). The whole-cell thinning pilot is a fixed cost that does not scale with the subset: its
   wall is reported separately, not in the ratio. The binding wall-time check stays K7 at the combined build (Step-2 §6:
   whole-cell pack wall time cannot be judged on a subset).

Selection (frozen A5, unchanged in substance): a candidate is eligible for a cell when conditions 1-5 hold; cost =
Σ over the cells of wall(s, cell) / wall(reference arm, cell); **the common default is the candidate eligible in both
cells with the least cost** among C1, C2, C3. If none: verdict `NO ELIGIBLE SETTING`, proposal = the candidate with the
fewest failed conditions summed over both cells (ties: least cost), for Pedro to decide; nothing is applied
automatically. The reference arms are not candidates.

Packing shapes (P8, adapted): the frozen check compares S4p against S4 (library trees); S4 is not run (§1 (c)). Instead
`shapes: pack_shapes` is part of the common default iff Jurkat C1 against R-J (the same settings with library trees)
has: C7(a) packed-packed ≤ 1 per 10⁶ on the library trees, placed_of_unthinned ≥ R-J − 0.003, every gated region ≥ its
R-J value − 0.003 (the binomial allowance of 1.), and subset pack wall ≤ wall(R-J). If the check fails while the chosen
candidate carries packing shapes, `selection.json` says so and the shapes value is for Pedro to decide.

`selection.json` (written by the script): per cell and run the conditions with values, eligibility, cost; the chosen
candidate; the shapes decision; and `common_defaults` = the exact dict for `COMMON_DEFAULTS`
(`wholecell/pack/switches.py`), to be pasted only after Pedro's decision.

## 4. What cannot be judged on a subset (deferred to the combined build, Part B K1-K9)

- Whole-cell placed mass and placed of unthinned vs the reference (K1) and per-region fractions over all tiles (K2);
  regions too small on the subset for condition 1 (< 100 requested copies).
- Seam index on the suite's whole-cell sample (K3); C7(b) yeast bead clashes (K4); whole-cell C8 (K5); C11 densities
  (K6); the binding pack wall time on the same node (K7) incl. the thinning pilot; K8, K9 no-regression.
- Thinning factors are whole cell (pilot) and are reported here, but their effect on whole-cell densities is K6.
- Frozen S2, S3, S5 and S4 with library trees are not run (§1); the frozen S4p-vs-S4 shapes comparison is replaced by
  C1-vs-R-J (§3).
- Seed noise: no replicate seeds (Pedro's rule); judged at the combined build with the §4 tolerances.
