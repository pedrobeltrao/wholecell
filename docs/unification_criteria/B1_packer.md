# B1 packer common defaults: benchmark protocol, selection rule and acceptance criteria (frozen before any B1 run)

Step B1 of the unification plan (`docs/reference/history/unification_plan.md` §2b P1-P9 and "Benchmark before freezing
defaults", §5, §5b B1): the yeast packer algorithms ported in A7 (size-class rounds, pilot thinning, shell cap,
per-region budgets, squeeze, packing shapes) get ONE set of cell-independent defaults, chosen by a benchmark. Order of
steps: B6 → B3 → B4 → **B1** → B2 → B5. Written 2026-11-17 before any benchmark run, B1 code or B1 build existed. The only
numbers used are those of existing builds (named with their paths). The md5 of this file goes to PLANNING_NOTES
(coordinator). Later findings go to a dated addendum file, never into this one.

Pedro's decisions that bind this step: the squeeze default is **targeted** (only underfilled (tile, region) pairs); **no
replicate-seed packing runs** just to estimate variance (seed noise is judged with the §4 tolerances and, for yeast, the
existing X5 s0/s1/s2 spread).

Plan §2b: "on a 500-tile Jurkat subset and on the whole yeast cell, measure placed fraction, worst-region fraction, seam
index, clash rate and wall time for: 1 v7 settings; 2 + rounds; 3 + rounds + thinning; 4 + targeted squeeze; 5 + full
squeeze (yeast only). Choose the cheapest setting that reaches placed ≥ 0.98 by mass in every volume region with ≥ 1 µm³
and seam index ≤ 1.1." Plan §5b B1: "placed mass fraction ≥ reference − 0.003 overall, and ≥ 0.98 in every volume region
≥ 1 µm³ (else explained); seam index ≤ 1.1; C7(a) ≤ 1 per 10⁶ packed–packed and C7(b) ≤ 1e-5; C8 ≤ reference; densities
per region ± 3 mg/ml unless placed mass changed (then mass-consistent); pack wall time ≤ 2× reference on the same node".
Made concrete below.

## Part A: benchmark protocol

### A1. Inputs (copies of the inputs of the latest complete build of each cell)

The benchmark runs in the sequential build queue (`wholecell_builds/_queue/jobs.tsv`, job `b1_benchmark`), after the
B3 / A9 / B6 / B4 builds queued before it, from a `git archive` snapshot of the B1 mechanism commit. It reads, never
writes, an input build; everything is written under `{results_root}/_benchmarks/b1_<commit7>/`.

| Cell | Input build (first that is complete: snakemake "100%" + `pack/pack_done.json`) |
|---|---|
| Jurkat | the B4 build (`jurkat/_dev_unify/b4/full_build.status`) → the B3 build (`.../b3/full_build.status`) → the B6 build `jurkat/20261117-04dec6f-s20261115-b6` |
| Yeast | the B4 build (`yeast/_dev_unify/b4/full_build.status`) → the B6 build → the A9 core build (`yeast/_dev_unify/a9/full_build.status`) → the A7 core run `yeast/_dev_unify/a7y/full_7a70243` (its inputs; packed with the snapshot code, which contains the membrane-code kernel fix 9977542) |

The input used is written to `inputs.json` (path, commit, complete flag). Per setting a fresh run folder holds a real
`pack/` folder; the input build's read-only inputs are symbolic links (`geometry`, `library`, `preplaced`, `inventory`,
`engine` for yeast, and the files of its `pack/` folder except the outputs of the counts / pack / summary steps:
`tiles.parquet`, `tiles_meta.json`, `tile_geom/`, `tile_labvol.parquet`, `organisation/`, the pre-placed index
`preplaced_*`, `mrna_*`). The config is the input build's own config (`<build>/_code/config/cells/<cell>.yaml`, extends
resolved) plus the setting (A3) as a preset; steps run: `counts`, `pack` (24 workers), `summary`, then the metrics (A4).
Jurkat packing shapes (setting S4p) need `library/pack_shapes/` in the input build (B3 and later); without it S4p is
skipped and reported.

### A2. The Jurkat 500-tile subset (deterministic; no seed needed)

Gated Jurkat volume regions (≥ 1 µm³ and copies requested in the reference builds): cytosol, nucleoplasm, matrix,
ER_lumen, golgi_lumen, endo_lumen, nucleolus, speckle (`perinuclear`, 7.9 µm³, has 0 requested copies: reported only).

Rule, on the input build's `pack/tiles.parquet`: among all bars of 5 × 5 tiles along the cell's long axis k (tiles with
ti ∈ [a, a+5), tj ∈ [b, b+5), any tk) holding ≥ 500 tiles, keep each bar's 500 tiles nearest its mean tk (|tk − mean|,
ties by tile id) and take the bar whose minimum over the gated regions of (region volume in the 500 tiles / region volume
in the cell) is largest (ties: smallest a, then b). The bar crosses PM, cortex, cytoplasm, mitochondria, ER, Golgi,
endosomes, NE, nucleoplasm, nucleolus and speckles; boundary tiles have unpacked neighbours (the same for every setting).

On the B6, B3 and v7repro2 builds this gives a = 24, b = 22, tk 12-31, **500 tiles, md5 of the sorted int64 tile ids
`0877fd0647b42b82db7d14bf5a64998e`**; region coverage 7.0 % (cytosol) to 30.7 % (nucleolus). The benchmark recomputes the
subset on its input build, writes `jurkat/subset_tiles.npy` and records whether the md5 is this one (a different md5 is
reported, not an error: the rule decides). v7 CPU time of these tiles in the B6 build: 31,019 CPU-s (≈ 22 min on 24
workers).

Yeast: the whole cell (5,912 tiles of 200 nm).

### A3. Settings

Only the algorithm keys change; every other `pack:` key (regions, tile size, halo, requests, pre-placed index, thin
lumens, lumen alignment, exclusion zones, organisation layer, packing shapes of yeast, seeds) is the input build's. Yeast
keeps its pinned preset `legacy_x5` and the setting is merged after it (preset list `[legacy_x5, b1_sN]`); Jurkat has no
pinned preset. Each setting sets every algorithm key explicitly (definitions: `wholecell/pack/switches.py`
`B1_BENCHMARK`, the same values as below).

| Setting | Budgets (tries per copy / failures in a row) | rounds | thinning (+ shell cap) | squeeze | shapes | Cells |
|---|---|---|---|---|---|---|
| S1 v7 settings | 200 / 25, thin lumens × 10 / × 16 (`thin_try_f`, `thin_giveup_f`), no per-region budgets | n 1 | off; shell cap off | none | cell's | both |
| S2 + rounds | as S1 | n 6, equal_volume, alternate | off | none | cell's | both |
| S3 + thinning | as S1 | as S2 | on: pilot_tiles 16, min_rel 1.1, safety 0.97, keep_above 0.98, reference cytosol, pilot_squeeze off; **shell cap on** | none | cell's | both |
| S4 + targeted squeeze | as S1 | as S2 | as S3, pilot_squeeze on | underfilled, threshold 0.97, rounds 4, iterations 60, insert_frac 0.15 | cell's | both |
| S5 + full squeeze | as S1 | as S2 | as S4 | all (rounds 4, iterations 60, insert_frac 0.15) | cell's | yeast only |
| S4y budget variant | 2000 / 1200 in every region (thin-lumen factors 1) | as S2 | as S4 | as S4 | cell's | both |
| S4p shapes variant | as S4 | as S2 | as S4 | as S4 | `pack_shapes` (library/pack_shapes, 15 residues per sphere) | Jurkat only |
| R reference | the input build's own pack, measured with the same code (not re-packed) | | | | | both |

The shell cap (P3) goes with thinning (both make "requested = placeable"). `thinning.resplit` stays the cell's (yeast
legacy_x5 requests: true). Yeast "cell's shapes" = its packing trees (`pack_trees_15`, legacy_x5); Jurkat's = the library
trees. Tile size is not benchmarked (stays 500 nm Jurkat, 200 nm yeast; the plan's 300 nm option is out of B1's scope).

Order of runs: Jurkat S1, S2, S3, S4, S4y, S4p; yeast S1, S2, S3, S4, S4y, S5; then R of both. One run at a time, 24
pack workers, the same node (`eu-g4-002`, 32 cores, shared). **Timeouts** (a run killed at its timeout is ineligible: it
already exceeds the 2 × wall-time limit): Jurkat 2.2 × the wall time of Jurkat S1 (at least 30 min); yeast 2.2 × the
input build's pack wall time (its `pack_pack` log; 10 h if unknown). A failed run is recorded with its log and the
benchmark continues.

### A4. Metrics (per run; `results.json`, `results.tsv`)

All placed fractions are by mass (Jurkat `libs.json` mass_da; yeast library bead mass, aa 110 Da / nt 330 Da, as suite
C10) and over the packed tiles of the run (Jurkat: the 500 subset tiles; yeast: all tiles).

- **placed_mass_fraction** = Σ placed × mass / Σ request × mass, request = the packer's request (`job_n` of the tiles:
  after thinning when thinning is on). Also **placed_of_unthinned** = Σ placed × mass / Σ `counts.npz` copies × mass of
  the same tiles (the request before thinning; = placed_mass_fraction without thinning).
- **region fractions**: the same two per volume region; **worst_region** = the minimum over the gated regions (Jurkat
  A2; yeast: volume regions with ≥ 1 µm³ in the cell by `tiles.parquet` and copies requested: cytosol, nucleoplasm,
  vacuole, cell_wall in the s0 build; nucleolus 0.5 µm³ and mito_matrix 0.8 µm³ are reported, not gated).
  Thinning factors and the request mass removed by thinning are reported.
- **seam index** (suite C10 definition and code, `wholecell/validate/suite.py` `seam_counts` / `seam_index`: packed
  instance centres per nm³ within 10 nm of a tile face / per nm³ ≥ 20 nm from every face). Yeast: the suite's own sample
  (`validate.suite.seam`: cytosol ≥ 0.95 of the tile, 300 tiles, seed 0). Jurkat subset: the subset tiles whose six
  face neighbours are all in the subset, with ≥ 0.95 of their volume in **nucleoplasm** (the subset holds only 2 interior
  tiles ≥ 0.95 cytosol; ~50 nucleoplasm), all of them; the same on cytosol ≥ 0.8 is reported.
- **C7(a)** (core `validate.placement`: sphere-tree overlap > 0.5 nm, 6 nm grid): `packed_packed_per_1e6` and
  `clash_rate_per_1e6` (packed-packed + packed-pre-pass). **C8(a)** `containment_rate_per_1e6`. Tiles: yeast the
  config's seeded sample (stream `validate.placement.tiles`, 300 tiles, as the suite); Jurkat 300 subset tiles drawn with
  `numpy.random.default_rng(0)` from the sorted subset ids. Collision shapes as the suite: Jurkat library trees (also for
  S4p: the gate is on the shapes the viewer and simulations use; S4p is also measured on its packing trees, reported),
  yeast `pack_trees_15`. C7(b) (yeast bead overlap, adapter code) is not measured in the benchmark (build level only).
- **wall time**: the pack step's elapsed time (`/usr/bin/time`, includes the thinning pilot and the squeeze); the counts
  step reported separately. CPU-s summed over tiles (`pack_log.tsv`) reported.

### A5. Selection rule (`selection.json`, written by the benchmark script)

A setting s is **eligible for a cell** when all hold:

1. every gated volume region has placed_mass_fraction ≥ 0.98 (of the packer's request; plan §2b);
2. |seam index − 1| ≤ 0.1 (plan "≤ 1.1"; a depleted face, < 0.9, is a seam too; the reference builds have 1.042
   Jurkat B6 and 0.981 yeast s0);
3. placed_of_unthinned ≥ that of S1 − 0.003 (thinning may not buy (1) by dropping mass: plan §5b overall rule);
4. C7(a) packed_packed_per_1e6 ≤ 1, and clash_rate_per_1e6 and C8(a) ≤ R's value + max(10 % rel., 1 per 10⁶);
5. wall time ≤ 2 × the reference wall time: Jurkat 2 × wall(S1) (S1 = the reference's settings, same tiles); yeast
   2 × the input build's pack wall time (legacy_x5, the reference settings, whole cell).

Cost of s = Σ over the cells of wall(s, cell) / wall(S1, cell). **The common default is the eligible-in-both-cells
setting of least cost** among S1, S2, S3, S4, S4y. S5 (full squeeze) is not a common candidate (Pedro: targeted
squeeze); if yeast is eligible only with S5, `selection.json` says so and the yeast value is a yeast pin to decide by
Pedro. If no setting is eligible in both cells: verdict `NO ELIGIBLE SETTING`, and the proposal is the setting with the
fewest failed conditions summed over both cells (ties: least cost), for Pedro to decide; nothing is applied automatically.

**Packing shapes (P8).** `shapes: pack_shapes` becomes part of the common default iff Jurkat S4p, against Jurkat S4, has:
C7(a) packed_packed ≤ 1 per 10⁶ on the library trees, placed_of_unthinned ≥ S4 − 0.003, every gated region ≥ its S4
value − 0.003, and wall time ≤ wall(S4). Otherwise the common default keeps `library` trees and yeast keeps its packing
trees as a cell value (as now).

`selection.json` holds: per cell and setting the conditions with values, eligibility, cost; the chosen setting; the
shapes decision; and **`common_defaults`**: the exact dict to paste into `COMMON_DEFAULTS` in
`wholecell/pack/switches.py` (the ONE place that sets the common defaults; Part C).

## Part B: build-level acceptance (one build per cell from the B1 commit, `scripts/build_isolated.sh <cell> b1 24`, in the queue)

**References.** The B4 build of the cell (`<cell>/_dev_unify/b4/full_build.status`, its `validate/suite.json`) when it is
complete; else the latest accepted build of the cell: Jurkat B3 → B6 (`20261117-04dec6f-s20261115-b6`, accepted with the
B6 addendum) → v7repro2; yeast B6 → A9 core → s0 (`yeast/_dev_unify/a8p2/20261116-aabe7b5-s0-x5repro/suite.json`). The
verdict names the reference used. Wall times are compared only between builds run on the same node with 24 cores
(`launch.log` host; otherwise the time check is reported, not gated).

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| K1 | Overall placed mass | placed_of_unthinned (A4 definition, all tiles) of the B1 build vs the same on the reference build | ≥ reference − 0.003 |
| K2 | Volume regions ≥ 1 µm³ | placed_mass_fraction per gated region (A4: of the packer's request), all tiles | ≥ 0.98 each, else listed with its cause (explained: a geometry limit, e.g. a thin lumen narrower than the species; or a thinning factor < 1 recorded in thinning.json) |
| K3 | Seam | suite C10 `seam_index` (cytosol ≥ 0.95, 300 tiles, seed 0; whole cell) | \|index − 1\| ≤ 0.1 |
| K4 | Clashes | suite C7 `a.packed_packed_per_1e6`; C7(b) `b.bead_clash_rate` where computed (yeast) | C7(a) ≤ 1; C7(b) ≤ 1e-5 on `b.bead_clash_rate_without_actin_structures` (the pre-placed actin structures are not packer output; s0: 9.86e-6) and `b.bead_clash_rate` ≤ reference + max(10 % rel., 1e-6) (s0: 2.52e-5). Jurkat C7(b): n/a (no bead-level metric) |
| K5 | Containment | suite C8 `a.containment_rate_per_1e6` (yeast also `b.beads_outside_region_fraction`) | ≤ reference + max(10 % rel., 1 per 10⁶) (rate tolerance of §4) |
| K6 | Densities | suite C11 `density.<region>` (mg/ml) per region | \|Δρ\| ≤ 3 mg/ml; or, when the region's placed mass changed, mass-consistent: \|Δρ − Δm_placed / V_region\| ≤ 2 mg/ml (§4 density tolerance) |
| K7 | Wall time | `pack_pack` rule wall time (snakemake rule start → finish, same node, 24 cores) | ≤ 2 × reference |
| K8 | Settings applied | `pack/pack_done.json` switches = `COMMON_DEFAULTS` merged on the cell (yeast: unpinned, or pinned values listed by selection.json) | yes |
| K9 | No regression elsewhere | suite `no_regression.tsv` vs the reference | no regression outside the EXPECTED list below |

**EXPECTED changes (both cells; caused by B1, so not regressions within these bounds).**

- Placed copies and mass per region: up (rounds, squeeze) or down by thinning in thinned regions (factors in
  `pack/thinning.json`); C10 placed fractions are judged by K1/K2, not the ± 0.003 per-metric rule.
- Shell cap: copies move from capped shells to their bulk volume region (`libs.json: shell_cap`); shell placed fractions
  and C8 shell rows change (reported).
- C11 densities: per K6. C11 crowding (φ, probe volumes): report, consistent with the density change.
- C1/C2/C3/C6 localisation and copy metrics: species in thinned regions lose copies by the thinning factor (compartment
  totals scale; rule-level shares unchanged); ± 0.01 standard otherwise.
- C12 ribosome NN (free 80S are packed in Jurkat since B4): ± 0.5 nm standard; up to + 1 nm explained if cytosol is
  thinned.
- Tile RNG trajectories diverge everywhere (rounds change the order of draws): pack metrics differ at seed-noise level.
- Yeast, if the selected default lacks the full squeeze: placed fractions in yeast regions fall from the legacy_x5 level
  (s0: 0.9924 by count of the thinned request); judged by K1/K2 (K1 vs the yeast reference). Yeast with S4y-type budgets:
  budgets unchanged from legacy_x5.
- Jurkat with packing shapes on (if selected): C7(a) is measured on the library trees (the suite's Jurkat setting): any
  overlap of the finer library trees counts.
- Wall time: per K7.

**Verdict.** `<build>/validate/compare_B1/verdict.json` per cell: `ACCEPTED` if K1-K9 pass (K2 failures only with a
listed cause); otherwise `NOT ACCEPTED` with the failing items. The benchmark results (`selection.json`) are attached by
path.

## Part C: the default-switching mechanism (code, before the benchmark; no effective change)

`wholecell/pack/switches.py`: `COMMON_DEFAULTS` (empty until the benchmark is read) is applied as the built-in preset
`common_defaults` to every cell that does not pin a preset (`pack.preset: legacy_v7 | legacy_x5 | none`, or a list of
presets merged in order). Pins at this commit: `jurkat_v7_port` → `legacy_v7`, `yeast` (and `yeast_x5_core`) →
`legacy_x5`, `yeast_x5_port` → `none`; `jurkat` and `testcell` unpinned (they take the common defaults once they are
set). With `COMMON_DEFAULTS = {}` every cell's effective `pack:` section is unchanged (tests: `tests/test_pack_switches.py`;
testcell CI check `tests/test_pack_smoke.py`). Applying the selection = pasting `selection.json: common_defaults` into
`COMMON_DEFAULTS`; for yeast also removing the `legacy_x5` pin (its non-algorithm keys already sit in its base `pack:`
section).
