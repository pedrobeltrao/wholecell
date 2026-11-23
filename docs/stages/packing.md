# Stage 7: packing

**Purpose.** Place every copy in the pack table inside its region, without overlaps, around the pre-placed
structures. The cell is cut into cubic tiles that are packed independently in phases, so neighbouring tiles never run
at the same time.

**Today.** Folder and CLI stage `pack`. Snakemake rules (`workflow/rules/pack.smk`): `pack_tiles` -> `pack_counts` ->
`pack_index` -> `pack_pack` -> `pack_summary`. With `pack.organisation` on, the step `organisation` runs after `tiles` and
before `counts` (rule `pack_organisation`, plan key `pack_organisation`).

**Inputs.** Geometry package, `inventory/pack_table.tsv`, the library manifest, `preplaced/preplaced_manifest.json`.

**Outputs** (`pack/`): `tiles.parquet`, `tiles_meta.json` (tile list and measures); `counts.npz`, `libs.json`
(expected copies per tile, species and region); `preplaced_{types,pos,quat,off}.npy`, `preplaced_libs.json` (pre-placed
index; `prepass_*` in builds before A0c); `tiles/<i>.npz` (placed instances per tile); `pack_log.tsv`; `pack_done.json`; `summary.json` (placed vs
expected mass, overall and per region).

**Config keys** (`pack:`): `tile_nm` (500 nm Jurkat and testcell, 200 nm yeast), `vox_nm`, `mesh_h_nm`, `meas_h_nm`,
`meas_vox_nm` (voxel size of the tile volume measures; default `meas_h_nm`), `tiling` (below), `requests` (below), `halo_nm`,
`max_tries`, `give_up`, `thin_try_f`, `thin_giveup_f`, `seed_prune`, `lumen_align`, `shell_area_factor`,
`thin_lumens`, `regions` (ordered; kinds `volume`, `virtual`, `membrane`, `shell`, `preplaced`), `virtual_providers` (Jurkat),
`region_map`, `region_map_default`, `mapping_rules`, `area_factors`, `open_face`, `index_drops`, `numba_threads`;
the algorithm switches below. Spec: [../reference/cell_config.md](../reference/cell_config.md) §3.5.

**Algorithm switches (A7 phases 1 and 2; yeast X5 algorithms ported into the core packer).** Defaults and meaning:
`wholecell/pack/switches.py`. With every switch at its default the packer runs the v7 code path (bit-identical; checked
on the Jurkat smoke tiles and the whole testcell).

| Switch | Default (v7) | What it does | Code |
|---|---|---|---|
| P1 `rounds: {n, split, alternate}` | n 1 | n visits per tile; non-membrane jobs split by bounding radius into n classes of equal requested sphere volume (largest first), membranes in round 0; phase order reversed in odd rounds. Intermediate state in `pack/rounds/<i>.npz` (removed when the tile is final) | `algorithms.round_edges`, `stage._pack_tile` |
| P2 `thinning: {enabled, pilot_tiles, ...}` | off | volume regions asking for > `min_rel` x the reference region's sphere volume fraction are pilot-packed (`pilot_tiles` non-adjacent tiles, no neighbours); a region that held < `keep_above` is thinned to `safety` x held (binomial thinning of its counts). Writes `thinning.json`, `counts_thinned.npz`; `summary.json` then also gives the fraction of the unthinned request | `stage._thinning`, `algorithms` |
| P3 `shell_cap` | off | a shell asks for at most the requested density of the bulk volume region on its side; the excess goes to that region (counts step, `libs.json: shell_cap`) | `algorithms.shell_cap` |
| P4 `budgets.per_region` | {} | per-region tries per copy / failures in a row (replace the global and thin-lumen budgets there) | `pack3d.pack_tile` |
| P5 `squeeze: {scope, threshold, ...}` | none | after the last round: insert missing copies, rigid-body overlap relaxation, exact RSA re-acceptance. `underfilled` = only label volume regions of a tile below `threshold` placed by mass; `all` = every tile with a deficit | `squeeze.py` |
| P7 `exclusion: [...]` | [] | copies of the listed species keep their centre `r_nm` away from pre-placed instances of the `around` libraries | `pack3d._in_exc`, `algorithms.exclusion_*` |
| P6 `organisation: {enabled, string, genes, subunits, pool_scores, ...}` | off | step `organisation`: STRING personalized-PageRank profiles per packed type (mean over subunits), field loadings (128 leading components of the profile correlation, float64) and partner CSR (loose partners r >= 0.4; contact partners from the AF3 pool only if `pool_scores` is given; off in published builds) -> `organisation/{network.npz, partners.npz, network.json}`. Counts: tile split x exp(field log weight) with the crowding correction (seed stream `pack.organisation.field`). Pack: partner attraction in the kernel (anchors = pre-pass instances of packed types, neighbour tiles, earlier rounds); shell jobs test a partner proposal against their band. Tiles get `att_stats`, `summary.json` an `organisation` entry | `organisation.py`, `stage._organisation_field`, `stage._partners` |
| P8 `shapes: library \| pack_shapes` | library | collision trees of every packed and pre-placed species from `library/pack_shapes/pack_trees_<per>.npz` (15 residues per sphere; `library --step pack_shapes`, also for frozen libraries) instead of the library trees; trees made from another library file (sig != manifest md5) are refused | `stage._shapes` |

**Presets and the common defaults (B1 mechanism, 2026-11-17).** `pack.preset` names one preset or a list merged in
order (CLI `preset=a,b`); presets come from `pack.presets` or are built in (`switches.builtin_presets`). A cell with no
`pack.preset` gets the built-in `common_defaults` = `switches.COMMON_DEFAULTS`, the ONE place that sets the common
packer defaults (algorithm keys only, `switches.ALGORITHM_KEYS`). It is empty until the B1 benchmark is read, so every
cell packs as before. Pins: `legacy_v7` (nothing merged: the cell's own keys with the v7 `DEFAULTS`; `jurkat_v7_port`),
`legacy_x5` (yeast, `yeast_x5_core`), `none` (`yeast_x5_port`); `jurkat` and `testcell` are unpinned. The yeast
`legacy_x5` preset holds only algorithm keys; its cell keys (thin lumens, bands tolerance, exclusion zones,
organisation layer, packing trees, pre-placed index, `thinning.resplit`) sit in its base `pack:` section (effective
section unchanged, md5-checked in `tests/test_pack_switches.py`). Benchmark settings `b1_s1` ... `b1_s5`, `b1_s4y`,
`b1_s4p` (`switches.B1_BENCHMARK`) and the benchmark: `scripts/b1_benchmark.py`, criteria
`docs/unification_criteria/B1_packer.md`; its `selection.json: common_defaults` is what goes into `COMMON_DEFAULTS`.

**Tiling and requests for every cell (A2).** `wholecell/pack/tiling.py` makes the tiles of any cell from its config
and geometry package (cellgeom_eval_v1): `pack.tiling.select` `centre_depth` (v7) or `enclosing_sphere` (yeast
`tiles3d.tile_list`), `box_margin_nm` / `pore_reach_nm` (evaluate a tile with only the prims and pores near it, yeast),
`labvol` (`tile_labvol.parquet`) and `meshes` (`tile_geom/<tile>.npz`: lo, measure, labvol, verts, faces, face_mem,
area; `mesh_membranes`). Region measures follow `pack.regions`: yeast bands are shells with `side: 0` and
`measure: volume`, the yeast wall depth is `sdf_lt`; no region ids are in code. `pack.requests.method: legacy_x5`
(`wholecell/pack/requests.py`) is the yeast `wc3d.distribute` split (stochastic rounding, multinomial with the
organisation log weights and crowding correction, thinning factors); its upstream inputs are declared until the core
makes them (per-type totals = X5 P1: A5/A6; library order and packing trees: A3; thinning pilot and organisation
loadings: A7). Check `tests/check_a2_yeast_tiles.py` (G4, G5, P3 vs the port build s0). Since A9 these inputs are core
outputs: totals and library order from `preplaced/serial/`, packing trees from the pre-placed context.

**A7 phase 2 checks.** Organisation layer on the yeast inputs vs the port build s0 (criterion P2 + addendum, core
environment, `scripts/organisation_a7_check.py`): pool on (X5 settings) passes at tolerance level: Z identical, the
loadings within 6e-8 after aligning one of 128 eigenvector signs (LAPACK of the core environment: the field
realisation of that component is mirrored), partner CSR identical up to one equal-strength tie at the 16-partner cut.
Pool off (published default) changes the partner rows of 74 % of the types (contact partners gone; 59 % of the
pool-on partner pairs kept). Testcell with the organisation layer on (human STRING v12, pool off): 45/45 types
with a profile, 232 loose links; 0.54 M of 2.92 M copies placed next to a partner; mass fraction unchanged (0.979).
Packing shapes on the testcell with the yeast switches on (rounds 6, thinning, shell cap, squeeze), 2 busiest colour-0
tiles: squeeze of the nucleus tile 636 s -> 203 s (whole run 673 s -> 237 s); the finer trees have a larger summed
sphere volume, so the pilot thins the nucleoplasm more (x 0.62 vs x 0.78). Tests: `tests/test_pack_organisation.py`. **A7 yeast side (2026-11-17).** The yeast cell packs through the core packer with `pack.preset: legacy_x5` (or
`preset=legacy_x5` on the command line; block `pack.presets.legacy_x5` of `config/cells/yeast.yaml`, one comment per
value with its X5 origin; presets are merged by `switches.effective`). Generic pieces: (1) region tables in any order:
the kernel indexes regions in its own order (volume-like, membranes, area shells; `Regions.KNAMES`, `k_of`, `c_of`;
counts and tiles keep config ids; identity for Jurkat / testcell); shells measured as volumes (yeast bands) are packed
as volume regions from their voxels (labels + band, the tile-measure rule `tiling.region_voxels`), with the band test
for partner proposals and squeezed copies (`volume_shell_tol_nm`); the pack evaluator uses the tiling box subset
(`pack.tiling.box_margin_nm`) and buffers for packages with more than 8 membranes. (2) Settings: rounds 6, pilot
thinning with `pilot_squeeze` and `resplit` (legacy_x5 requests re-split from the thinned totals), shell cap (also on
legacy_x5 totals, volume-measured shells by their measured volume), budgets 2000 / 1200 everywhere, squeeze `all` (every
label-volume molecule of a tile moves, yeast semantics), exclusion zones from a registered provider
(`around: {provider: yeast.endocytic_sites}`, `plugins/yeast/exclusion.py`), organisation with the AF3 pool,
`shapes: pack_shapes`. (3) `legacy_x5` requests read `pack/organisation/network.npz` (float64 loadings) when
`requests.organisation` is not given. (4) `pack.preplaced_index`: pre-placed instances from one grouped npz (yeast
`prepass.npz`, groups in `wc3d.PRE_GROUPS` order, beads as obstacles) -> the core index (+ `preplaced_groups.json`).
(5) Snakemake rule `pack_organisation`. Check `scripts/pack_a7_yeast_check.py` (P3 / P5 / P6 vs s0 with the
addendum seed spread, clash / containment / exclusion; subset or whole cell); the whole-cell run is prepared in
`_dev_unify/a7y/full_run.sh`. Subset (20 busy tiles: nucleus, cortex, mitochondria, ER, vacuole; s0 neighbours as
context, s0 pre-pass, core request split): placed / requested 0.9819 vs s0 0.9828; per-region fractions within 0.01
except mem_MITO_IM (0.924 vs 0.944); per-(type, region) shortfalls within max(3 sqrt n, 1 % n) for 5100 of 5106 pairs (the
P3 / P5 subset differences follow the requests, which the core split draws differently); clashes > 0.5 nm 0, centres
outside their region 0, ribosomes in exclusion zones 0 (s0 likewise); CPU 3344 s vs s0 2462 s. Jurkat smoke and testcell
(default switches) identical to v7repro2 / CI 1d6a575. Since A9 the yeast build packs through the core (`pack.preset: legacy_x5`); remaining: partners in the
narrow pore-cylinder shells still ignore partner proposals; Common defaults are set by the B1 benchmark.
through the core packer with Y settings against the X5 P3/P5/P6 criteria. Common defaults are set by the B1 benchmark.

**Checks.** `summary.json` mass fraction per region; validation steps `placement` and `densities`. Jurkat v7
reproduction B6-B11 (counts, index, per-tile and per-species placement, totals: 880,040,974 instances, 87.72 pg,
mass fraction 0.964). Testcell build `ci2`: 2,920,324 instances, mass fraction 0.979.

**Code.** `wholecell/pack/{stage.py, pack3d.py, regions.py, virtual.py, geomload.py, switches.py, algorithms.py, squeeze.py, requests.py, organisation.py}`; Jurkat virtual regions in
`plugins/jurkat`.

**Jurkat vs yeast today.** Yeast packs with the core packer, `pack.preset: legacy_x5` (A9): 200 nm tiles, six
size-class rounds, squeeze, pilot thinning, shell caps, organisation layer (**A7** ported these behind switches with
legacy defaults). `pack.engine` (`yeast_rsa_v1`, the vendored stage 3) remains only in the reproduction config
`config/cells/yeast_x5_port.yaml`. **B1** sets common defaults
after a benchmark; **B2** switches the organisation layer on (STRING; AF3 pool off). Until B2 is accepted its settings live in the overlays
`config/cells/jurkat_b2.yaml` (organisation on, STRING v12 9606, combined score >= 400) and `config/cells/yeast_b2.yaml`
(pool off; `yeast.yaml` keeps the X5 pool setting), built by the B2 queue job; check `scripts/b2_organisation_check.py`,
criteria `docs/unification_criteria/B2_organisation.md`.

**Methods text (draft).** The cell was divided into cubic tiles (500 nm for Jurkat). The expected copies of each
species, minus copies already used by pre-placed structures, were mapped to packing regions and split over tiles by a
multinomial draw proportional to each tile's region volume or membrane area. Tiles were packed in eight
checkerboard phases (2 x 2 x 2 colours) so that tiles packed at the same time never touch; each tile saw the
pre-placed structures and the already packed neighbour tiles within a 40 nm halo as obstacles. Within a tile,
molecules were placed by random sequential addition with overlap tests on sphere trees. Membrane proteins were
sampled on the membrane surface and aligned to the local membrane normal. The result does not depend on the number of
workers.
