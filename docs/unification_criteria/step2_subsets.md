# Step-2 subset tests: frozen subset definition, packing rule, metrics and pass rules

Written 2026-11-18 (Pedro: replace the per-step full builds of Step 2 by small, fast tests) before any subset test was
run. Frozen ONCE: the tile lists, the neighbour rule, the seeds and the metric / pass definitions below do not change
for the rest of Step 2. Later findings go to a dated addendum file, never into this one. Harness:
`scripts/step2_subset_test.py` (written after this file; it re-derives the tile lists from the rule in §1 and refuses
to run when the md5 differs from §1).

## 1. Subsets (data files)

| Cell | File | Tiles | md5 (file) | md5 (`ref_tile`, int64 bytes) | Frame |
|---|---|---|---|---|---|
| Jurkat (500 nm tiles) | `config/subsets/jurkat_step2.tsv` | 56 | `22ca415ad8e61a73d8937640ca4961c6` | `dd7c4bf95419d7a4d5acc6358da8c359` | tiles of the accepted B3 build `jurkat/20261117-bc121a6-s20261115-b3` (7,485 tiles) |
| Yeast (200 nm tiles) | `config/subsets/yeast_step2.tsv` | 58 | `54276fc4c849130cc66eff98ddf3d408` | `cf4071d2756f76b84dbccb4f35ae083f` | core tiles of `yeast/20261118-9977542-s0-a9-core` (5,912 tiles; the same tile set, row by row, as the port s0 with grid offset `-lo0 / 200`) |

Columns: `cell, ref_tile` (row of the frame's `tiles.parquet`), `ti, tj, tk` (tile grid index in that frame; THE
identity of a tile: tests map it to their own tile id by grid index), `requested_ref` (instances requested in the frame
build's `counts.npz`), `reasons`.

Rule (deterministic, frame build's `tiles.parquet` + `counts.npz`; "take" = the next tile of a ranking that is not
chosen yet, ties by lowest tile id):

1. Whole small regions (yeast): every tile of the smallest 26-connected component of the tiles with
   `peroxisome | mem_PEROX | PEROX_peripheral > 0` (15 tiles: one whole peroxisome) and of the tiles with
   `LD_peripheral > 0` (8 tiles: one whole lipid droplet; FAA4 / LD_peripheral check).
2. Every region measure column with a non-zero tile (volumes, shells measured as volumes, virtual regions, membrane
   areas; Jurkat 20, yeast 25 columns): take its richest tile; Jurkat also the tile closest to the 75th percentile of
   the non-zero measure (a typical, not only an extreme, tile).
3. Busiest tiles (requested instances): take 6 (Jurkat) / 4 (yeast). Membrane-rich tiles (sum of `mem_*` areas): take
   4 / 3.
4. Seam tiles: seeded permutation (`numpy.random.default_rng(20261118)`, sorted candidates) of the tiles ≥ 95 % one
   region: take 3 cytosol + 3 nucleoplasm (Jurkat), 3 cytosol (yeast).

Coverage: every region of `libs.json` `rnames` with requests in the frame build has requests in the subset (Jurkat 37
regions, yeast 27 incl. the virtual `chromatin_nucleosome` / `NE_pore_NPC`). Subset share of the cell's requested
instances: Jurkat 1.1 %, yeast 2.0 %. Fewest subset tiles per measure: Jurkat speckle 2, nucleolus 3, mem_PM 8; yeast
nucleolus 2, nucleoplasm 5, mem_NE 6. Region measure shares in the subset: Jurkat golgi 27 %, pcm 54 %, others 0.3-4 %;
yeast peroxisome 50 %, mem_PEROX 36 %, LD_peripheral 12 %, others 0.3-6 %.

## 2. What a test runs (harness)

Into a NEW folder `{results_root}/<cell_id>/_dev_unify/step2/<name>/` (never into an existing build): a `git archive`
snapshot of the commit under test (`_code/`), a config `cell.yaml` = `extends:` the snapshot's cell config (or overlay,
e.g. `jurkat_b2.yaml`) + `build.allow_dirty: true`; then the build DAG of that config (Snakemake dry run of the
snapshot's `Snakefile`, targets `pack/counts.npz` + the pre-placed index) executed job by job:

- geometry (links the package), inventory, library, pack tiles, pre-placed (whole cell: builders are not
  region-restricted), pack organisation / counts / index: all whole cell, as in a build;
- **reuse**, declared per run and recorded in `result.json`: (i) the library of the latest accepted build
  (`--reuse-library <build>`; Jurkat: the B3 build) unless the test changes the library: allowed only when the
  `library:` section of the resolved configs is identical; (ii) the whole upstream of an existing build of the SAME
  commit and config (`--reuse-upstream <build>`: symlinks, only the subset pack is run); (iii) the pack tiles
  (`tiles.parquet`, `tiles_meta.json`) of a build whose geometry package (manifest sha256) and `pack.tiling` /
  `tile_nm` are identical;
- the thinning pilot (yeast: `pack.thinning`, whole cell, its own pilot tiles packed in isolation as in a build);
- the pack of the subset tiles only (§3), `workers` ≤ 24;
- the metrics on the subset (§4) → `result.json`; with `--ref <run or build>` the comparison (§5) → `compare.json`.

Seeds: the config's own (Jurkat master 20261115, yeast master 0), never changed by the harness. The tile seed stream is
keyed by the tile id; a test's tile ids must equal the reference's for the same grid index (checked; the tile set does
not depend on the B-steps of Step 2). Hence with identical code and inputs a subset tile is the SAME tile as in the full
build (bitwise for Jurkat, §3).

## 3. Neighbour rule (reference context, the same for every setting of a cell)

Only the subset tiles are packed, in the full build's order (8 checkerboard colours; yeast: 6 size-class rounds with
alternating colour order). A neighbour (26-neighbourhood) that is itself in the subset is packed by the test and seen
as in a build. A neighbour that is NOT in the subset is taken from the declared **context build** of the cell:

- Jurkat: the accepted B3 build `jurkat/20261117-bc121a6-s20261115-b3` (`pack/tiles/<i>.npz`);
- yeast: the port s0 `yeast/20261116-aabe7b5-s0-x5repro` (`engine/work/yc3d/tiles_raw/<name>.npz`).

Its instances are re-indexed into the test's library by species name (instances of species not in the test's
`libs.json` are dropped and counted in `result.json`), and it is visible only when it would have been in a full build:
subset tile t in round r sees context tile j with the jobs of size class ≤ r if j's colour comes before t's colour in
round r's phase order, else with the jobs of size class < r (one round, Jurkat: the whole tile or nothing). Size class
of a job = the test's own round thresholds on the test's bounding radii (tile files keep their instances grouped by
job). Yeast approximation: the context tile is the final one, so the last round's squeeze moves are visible in earlier
rounds (squeeze changes < 1 % of instances). Pre-placed instances, RNA / DNA beads and exclusion zones are the test's
own (whole cell). Why this rule: (1) every setting of a cell sees the SAME context, so differences come from the test;
(2) the subset tiles see what they would see in a full build (no empty halo: an isolated subset over-places at its
faces); (3) with the context build's own code and inputs the subset reproduces the context build's tiles bitwise
(Jurkat), which validates the harness.

## 4. Metrics on the subset (`result.json`)

All on the subset tiles only, computed with ONE code version (the harness commit) for test and reference alike; the
reference is either another harness run or a full build read through the same subset tile list.

| Id | Metric | Definition on the subset |
|---|---|---|
| S-C10 | placed fractions | overall and **per region** (config region ids): placed / requested copies and mass (`job_n`, `job_placed`, library masses); placed of the unthinned request (`counts.npz` of the subset tiles) |
| S-P3 / S-P5 | per (type, region) | requested / placed copies per (library, region) pair (`pairs.tsv.gz`) |
| S-C7a | clashes | core C7(a) (`wholecell/validate/placement.py`): sphere overlap > 0.5 nm with the packer's collision trees inside each subset tile (packed-packed, packed-pre-placed, packed-RNA); per 10⁶ packed; max depth |
| S-C8 | containment | core C8(a) on the subset tiles, per 10⁶ packed, per region kind |
| S-C9 | orientation | core C9 on the membrane instances of the subset tiles: sign errors per 10⁶, median / p99 tilt |
| S-seam | seam index | suite definition (band 10 nm / interior ≥ 20 nm) over the subset tiles ≥ 95 % cytosol (and, Jurkat, nucleoplasm) |
| S-C11 | densities | per volume region: (packed mass placed in the subset tiles, shells and virtual regions moved to their volume as in `validate.densities` + pre-placed mass with centres in the subset tiles, assigned by `validate.densities.preplaced_volume`) / subset region volume; RNA excluded (not tile-resolved in every layout) |
| S-C12 | ribosomes | NN median of the 80S centres (packed + pre-placed ribosome libraries) inside the subset tiles, raw (no EM thinning) |
| S-stage | stage level, WHOLE cell | md5 of stage outputs (pack table, library manifest, pre-placed files, counts, pre-placed index); requested copies per region (whole-cell `counts.npz`); pre-placed instances per library; pack index drops; thinning factors; region volumes / areas of `tiles.parquet` |
| S-time | timing | wall seconds per stage and for the subset pack; tile CPU seconds (`pack_log.tsv`) |
| yeast | extras | `mem_PEROX` and `LD_peripheral` rows; FAA4 (library id `FAA4`, UniProt P47912) per region |

## 5. Comparison and pass rules (`compare.json`)

Tolerances are the plan §4 / suite ones; where a subset count is small a binomial term is added so that seed noise
alone does not fail.

| Id | Pass |
|---|---|
| S-C10 per region | test ≥ ref − max(0.003, 3·√(f(1−f)/n)) for every region with ≥ 100 requested copies on the subset (f = ref fraction, n = test requested copies); one-sided (a rise is reported). Overall placed mass fraction: same rule |
| S-P5 | X5 P5 rule: pairs with ref ≥ 10 copies within max(3√n, 1 % n); ≥ 99 % of pairs with ref < 10 within ±3. S-P3 (requested): max(3√n, 0.5 % n). Judged where the requests are made the same way (same request code); for core vs port (A2: the core split draws per tile differently) S-P3 / S-P5 are REPORT, and the judged item is the per-pair placed fraction |placed/requested(test) − placed/requested(ref)| ≤ max(0.01, 3·√(f(1−f)/n)) for ≥ 99 % of pairs with ≥ 10 requested copies |
| S-C7a | packed-packed ≤ 1 per 10⁶ packed and clash rate ≤ ref + max(10 %, 1 per 10⁶) |
| S-C8 | ≤ ref + max(10 %, 1 per 10⁶) |
| S-C9 | sign errors ≤ ref + max(10 %, 1 per 10⁶); median tilt within ±0.5° |
| S-seam | \|index − 1\| ≤ 0.1 (B1 K3) |
| S-C11 | per volume region with ≥ 0.01 µm³ on the subset: \|Δρ\| ≤ 2 mg/ml, or mass-consistent: \|Δρ − Δρ_requested\| ≤ 2 mg/ml (Δρ_requested from the requested mass on the subset) |
| S-C12 | \|Δ NN median\| ≤ 0.5 nm |
| bitwise | number of subset tiles with identical `types`, `pos`, `quat` (report; MUST be all tiles where a test declares bitwise equivalence, e.g. a code change designed to be bitwise) |
| S-stage | md5 equal / different per stage (report; MUST be equal for the stages a test declares unchanged) |

Verdict per test: PASS when every judged item passes; otherwise FAIL with the failing items (a failure explained by a
documented expected change of the step is listed as EXPECTED, as in the step's own criteria).

## 6. Criteria that cannot be judged on a subset (deferred to the combined build)

Whole-cell totals: C10 cell placed mass / copies fraction and P6 instance totals; C11 whole-cell densities, N/C ratio
and φ (crowding sample); C12 thinned-to-EM NN rows and ribosome concentration; B4 J3 / J4 polysome statistics vs EM
(whole-cell point clouds; the pre-placed stage output itself is whole cell and is compared at S-stage); suite rows
computed over seeded whole-cell tile samples (C7 / C8 / C9 / seam as the suite defines them); seed-spread rules
(X5 s0/s1/s2 intervals); compare_v7 / compare_x5 sections; localisation (C4/C5 loc_compare) and export; pack wall time
(B1 ≤ 2× reference); C13 / C16 build bookkeeping. Stage-level whole-cell outputs that the harness does compute
(pre-placed counts, index drops, thinning factors, region volumes) are compared at S-stage.
