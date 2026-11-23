# Cell config schema (`config/cells/<cell>.yaml`, schema_version 1)

Status: spec v1, 2026-11-15 (SPECS agent), written after reading the in-progress stage code (`wholecell/inventory/
stage.py`, `generic.py`, `context.py`, `wholecell/pack/regions.py`, `virtual.py`, `geomload.py`,
`plugins/jurkat/*`). Sections marked **(draft)** have no consumer code yet and are specified from ARCHITECTURE.md and
the audits; they will be aligned when their stage lands. Loader: `wholecell.config.load`. Audit references:
H = `data/derived/repro_audit/human/`, Y = `data/derived/repro_audit/yeast/`.

## 1. Principle

A cell = this file + a cellgeom_v1 package + adapters + plugins (ARCHITECTURE §1). Every value that differs between
cells lives here or in a file it names; code holds no cell, compartment, region or species list. Why: the audits list
~40 hard-coded human lists (H cell_specific_vs_general §7) and ~60 yeast ones (Y cell_specific_vs_general §8) that
had to be edited in code per cell or per version (`RNAMES_V5/6/7`, `jpack_test.SOL/TMR`, `NOT_IN_THIS_CELL`, `TGRID`).

## 2. Top-level keys

| Key | Type | Req. | Owner | Content |
|---|---|---|---|---|
| `schema_version` | int | yes | core | 1 |
| `cell_id` | str `[a-z0-9_]+` | yes | core | = file stem; folder name of builds and geometry packages |
| `name` | str | yes | core | display name (viewer card, stamp) |
| `organism` | {name, taxon, strain?, assembly?} | yes | core | taxon ids drive adapter file choices (yeast 559292/4932, Y §8 "organism") |
| `aliases` | {alias: text} | no | core | e.g. `v7`; build-id aliases (H merge_notes §8 "adopt build ids, keep v7 as alias") |
| `status` | `stub` \| `dev` \| `published` | no | core | a `stub` config may leave sections empty; builds refuse `stub` |
| `extends` | path | no | core | parent config (§5) |
| `plugins` | [str] | yes (may be []) | core | packages under `plugins/` imported before any stage (`registry.load_plugins`) |
| `seeds` | {master, legacy?} | yes | core | §7 |
| `build` | {allow_dirty, tag?} | no | core | §7 |
| `geometry` | map | yes | geometry | §3.1 |
| `inventory` | map | yes | inventory | §3.2 |
| `library` | map | yes | pack | §3.3 (draft) |
| `preplaced` | map | no | pack | §3.4 (draft) |
| `pack` | map | yes | pack | §3.5 |
| `validate` | map | no | validate | §3.6 (draft) |
| `export` | map | no | export | §3.7 (draft) |
| `frozen_inputs` | {name: spec} | yes (may be {}) | all | §4 |
| `notes` | str | no | — | free text |

No other top-level key is accepted (§6). Stage owners add keys only inside their section.

## 3. Sections

### 3.1 `geometry`

```yaml
geometry:
  toolkit: em                      # em | analytic | <plugin toolkit>  (who made the package; informational)
  package: "{geom_root}/jurkat/cg1-20261115-v7conv"   # immutable cellgeom_v1 directory
  manifest_sha256: "..."           # checked by the geometry stage (the package is itself a frozen input)
  on_mismatch: error               # error (default) | warn: record a manifest mismatch in geometry/package.json and go on
                                   # (locally generated packages, e.g. the testcell's; the paper cells use error)
  engine: sdf_sections_v1          # package layout: must equal the manifest's engine.name
  evaluator: cellgeom_eval_v1      # default and only core value (legacy evaluators: reproduction tests only)
  rules: "{repo}/config/geometry/rules/sdf_sections_v1.yaml"   # rule table of the layout (default when absent)
```
The package carries frame, handedness, labels, membranes, label rule, layers (`cellgeom_v1.md`). Development-only
fallback `pack.dev_legacy_geometry` (legacy npz) exists in `geomload.py` and is refused for `status: published`.

### 3.2 `inventory` (consumer: `wholecell/inventory/stage.py`)

```yaml
inventory:
  adapters:                        # kind -> registered component name (null/absent = skipped)
    rna: jurkat.rna_e61
    proteome: jurkat.geiger_ibaq_v0          # generic: proteome_table
    localisation: localisation        # core evidence table + combiner (preset); generic table: localisation_table
    secms: secms_toolkit                     # optional: SEC-MS toolkit -> inventory/secms/*
    assembly: complex_table                  # none | complex_table | assembly_steps (core step list + preset) | <plugin>
  total_protein_pg: 100            # REQUIRED, no default
  total_protein_source: "decision 2026-11-08 (H parameters.yaml inventory.total_protein; CLI arg of inventory_v0.py)"   # REQUIRED
  ingredients: generic.monomers    # inventory_builder
  rule_sets: [jurkat.v4, jurkat.v5, jurkat.v6, jurkat.v7]   # ordered; one step + changes.tsv each
  params:                          # per component; file-valued keys name frozen inputs (PARAM_FILES / FROZEN)
    jurkat.geiger_ibaq_v0: {proteome: inv_geiger_ibaq, uniprot: inv_uniprot_2026_11}
    proteome_table: {path: "{data_root}/...", source_id: paxdb_4932, unit: ppm}
    generic.signal_lumen: {threshold: 0.5}
  verify: {pack_table_md5: d1a0099caad792d0ba42cbc52b52f0a5}   # optional: compare with a legacy reference
```
| Rule | Why |
|---|---|
| multiple sources with weights are params of one adapter: `sources: [{source_id, input, weight}]`, recorded in `proteome_sources.tsv` / `localisation_sources.tsv` (inventory_schema.md) | human HPA main 1 / additional 0.3 and HPA:UniProt 50:50 (Y comparison §(a)) are weights, not code |
| `total_protein_pg` + `_source` required | human 100 pg decision vs yeast 3.5 pg measured (H merge_notes §8) |
| copy numbers come from one generic step: `inventory.copies` (sources with type / reader / aggregation / group split, combination and priority, imputation, absent lists, per-protein scaling, total protein fixed / anchors / source_sum, adjustments, anchors, future sources; inventory_schema.md §4b) with `adapters.proteome: copies` | 2026-11-16: the Jurkat and yeast methods are options of one step (`docs/reference/history/copy_numbers_report.md`) |
| curated lists (`ZERO`, `NOT_IN_THIS_CELL`, `CURATED_W`, `OBLIGATE`, ...) are TSV/YAML files with a source column, named in `params` and declared in `frozen_inputs` (`inv_*`) or shipped under `plugins/<cell>/data/` | H cell_specific §7; Y merge_notes §10.6 |
| pack-dependent rules (Jurkat v5 P, v6 R1) read only declared frozen tables | H rebuild_test, ARCHITECTURE §3 |
| future data slots (SEC-MS, surface proteomics) are adapters `adapter.secms`, `adapter.surface` | DECISION 2026-11-15 (new data expected) |
| SEC-MS: `adapters.secms: secms_toolkit` + `params.secms_toolkit` = the toolkit run config (`species: human\|yeast`, `preset`, `annotation`, `complexes`, `pair_candidates`, `datasets`; `toolkits/secms/README.md`); assembly rules that use SEC are options of the assembly preset (`params.assembly_steps.preset` / `overrides`; `complex_table`: `params.complex_table.assembly_options`, `wholecell/inventory/secms.py`) | 2026-11-16: same code and parameters for every cell; only dataset properties are per cell |

### 3.3 `library` (draft)

```yaml
library:
  mode: frozen                     # frozen | build
  frozen: lib_structures_v7        # frozen input: directory + manifest (16,246 npz, sha256 per file)
  builder: {name: stage0_default, params: {plddt_min: 50, sphere_tree: adaptive_2p5x}}   # mode build
  pack_tree_level: 15              # residues per sphere for packing trees (yeast PACK_PER), viewer gets library trees
```
Why: human library mutated in place by 5 scripts vs yeast scripted stage 0 (H merge_notes §2); phase 1 uses frozen.
Build mode (implemented; `wholecell/library/stage.py`): `species` (a TSV, or `{from: pack_table, accs: {table, gene,
acc}}`: written from this build's pack table), `steps` (incl. `fallback`: placeholder / experimental-chain structures,
assemblies split into monomers), `ledger: {max_missing_mass_frac}` (`library/species_ledger.tsv`; rebuild v9).

### 3.4 `preplaced` (draft)

Called `prepass` before A0c (2026-11-17); old keys still load (`wholecell/config.py` `LEGACY_CONFIG_KEYS`). Seed
stream names keep `prepass.` (the name is the seed).

```yaml
preplaced:
  modules:                         # order = placement order (legacy PRE_GROUPS order is a contract of the encoder)
    - {name: chromatin, plugin: jurkat.chromatin_v8, mode: frozen, input: pre_chromatin_v8_rep5, depends_on_geometry: cg1-...}
    - {name: lamina, plugin: jurkat.lamina, mode: build, seed_stream: prepass.lamina, params: {...}, compare_frozen: pre_lamina_v7}
    - {name: mrna, plugin: generic.mrna_threads, mode: frozen, input: pre_mrna, drop: pre_mrna_drop}
```
`mode: build` runs the registered module; `frozen` reads a declared input. `depends_on_geometry` must equal
`geometry.package` id, else the build fails (implemented in B6, 2026-11-17: `preplaced --step manifest`, suite C16
`depends_on_geometry_mismatch`; a module built in the build depends on its package by construction) (human chromatin/MTs built on older geometry, H merge_notes §7.7;
cellgeom_v1 check 10). Yeast lists its builders here (rabl_chromatin, spb, eisosomes, ...; Y merge_notes §10.7).

### 3.5 `pack` (consumer: `wholecell/pack/regions.py`, `virtual.py`)

```yaml
pack:
  tile_nm: 500                     # yeast 200
  halo_nm: 40
  phases: 8
  grid_spacing_nm: 6.0
  budgets: {max_tries: 200, give_up: 25, thin_lumen_try_factor: 10, thin_lumen_giveup_factor: 16}
  thin_lumens: {ER_lumen: ER, golgi_lumen: GOLGI}       # region -> engine membrane
  regions:                         # ORDERED: volume/virtual, then membranes (engine order), then shells; index = job_region
    - {name: cytosol, kind: volume, labels: [cytosol]}
    - {name: mt_lattice, kind: virtual, labels: [cytosol], provider: cyto_v7, select: mt_lattice}
    - {name: PM, kind: membrane, membrane: PM}
    - {name: sh_PM_cortex, kind: shell, membrane: PM, labels: [cytosol], side: 1, d0: 0, d1: 150, add_radius: true}
  virtual_providers: [{id: cyto_v7, plugin: jurkat.cytoskeleton_v7, params: {...}, inputs: {ne_mesh: geom_ne_hmff}}]
  region_map: {cytosol: {soluble: cytosol, tm: PM}, er: {soluble: ER_lumen, tm: ER}}    # localisation key -> region
  region_map_default: {soluble: cytosol, tm: PM}
  mapping_rules: {tm: [{rule: om_vs_im, frac: mito_om_frac}], soluble: [{rule: signal_lumen}]}
  area_factors: {GOLGI: 1.296, MITO_IM: 1.0}
  drop: [no_structure, tm_in_soluble_compartment]       # named drop rules (yeast drops these two names)
  index_drops: [{name, libs, lib_prefixes, mode, ..., reason, return_to_pool: {lib: species}}]   # implemented; reason required
  ledger: {max_unexplained_mass_frac: 1.0e-4, explain_no_structure: null}   # implemented (rebuild v9): pack/copy_ledger.tsv
```
Kinds: `volume` (labels; optional `sdf_lt: {membrane, nm}` keeps voxels with that SDF < nm, e.g. the yeast outer 60 nm
of wall), `virtual` (registered `virtual_region` provider), `membrane` (engine membrane), `shell` (membrane, labels,
side (+1, -1, or 0 = both sides, |SDF|), d0, d1, `d1_mem_half`, `add_radius`, `pore_limit_nm`, `measure`: `area`
(default; the membrane's area) or `volume` (the voxels of the band; yeast cytosolic bands)), `preplaced` (an id kept for
pre-placed instances; no measure, not packed). The list position is the region id. The canonical order (volume/virtual,
membranes in engine order, shells; no preplaced, two-sided or volume-measured shells) is required by the core counts (v7
method) and pack steps; the yeast table keeps its legacy ids 0-26 and is used by the tiles step and
`requests.method: legacy_x5` (A2). Why: replaces `jpack_test.RNAMES*`, `SOL/TMR`, area factors in `jwc.step_counts`
(H formats §1b) and yeast `tiles3d.REGION_NAMES`/`REG_MASK`/`BANDS`/`WALL_DEPTH` (yeast.yaml, one comment per region).

**Tiling and requests (A2).** Defaults = v7; keys and meaning in `wholecell/pack/tiling.py` and `requests.py`.

```yaml
pack:
  meas_vox_nm: 4.0                 # voxel size of tile volume measures (default meas_h_nm; yeast 4)
  tiling: {select: enclosing_sphere, sample_reach_nm: 88, box_margin_nm: 80, pore_reach_nm: 119,
           labvol: true, meshes: true, mesh_membranes: [PM, ER, NE, VAC, MITO_OM, MITO_IM, GOLGI, PEROX, LD]}  # yeast
  requests: {method: legacy_x5, totals: ..., lib_order: ..., shape_volumes: ..., thinning: ...,
             organisation: {loadings: ..., beta: 0.8, field_l_nm: 500, n_features: 128}}   # declared inputs (yeast)
```

**Packer switches (A7, yeast algorithms).** Defaults (= Jurkat v7, bit-identical) and the meaning of every key are in
ONE place, `wholecell/pack/switches.py` (`DEFAULTS`); a cell lists only what differs. Validated with the key name and
a nearest-name suggestion.

```yaml
pack:
  rounds: {n: 1, split: equal_volume, alternate: true}       # P1: n size-class visits per tile (v7: 1)
  thinning: {enabled: false, pilot_tiles: 16, min_rel: 1.1, safety: 0.97, keep_above: 0.98, reference_region: cytosol}  # P2
  shell_cap: false                                           # P3 (or {enabled: true})
  budgets: {per_region: {nucleolus: {max_tries: 2000, give_up: 1200}}}   # P4 (v7: {} = global + thin-lumen budgets)
  squeeze: {scope: none, threshold: 0.97, rounds: 4, iterations: 60, insert_frac: 0.15}   # P5: none | underfilled | all
  exclusion: [{name: ribo_endo, libs: [ribosome_80S], around: {lib_prefixes: [ACT]}, r_nm: 90}]   # P7 (v7: [])
  organisation:                                              # P6 (v7: enabled false); step `pack --step organisation`
    enabled: true
    string: {taxon: 9606, version: "12.0", links: <STRING links .txt.gz>, min_score: 400,
             ids: {aliases: <STRING aliases .txt.gz>, source: UniProt_AC}}   # yeast: ids: {table: inv_yeast_uniprot, columns: ["Gene Names (ordered locus)"]}
    genes: {table: <UniProt TSV>, columns: ["Gene Names (primary)"]}          # gene symbols of complexes -> UniProt
    subunits: {table: null, lib_map: null, ids: null}        # default: the pack table (ing_id, uniprot, genes)
    pool_scores: null                                        # AF3 pool pairs: off in published builds (STRING only)
    paralogs: null
    # damping 0.5, field_k 128, field_l_nm 500, beta 0.8, r_loose 0.4, max_partners 16, a_contact 0.6, a_loose 0.25,
    # k_half 8, p_max 0.6, attract_tries 40, loose_gap_nm [2, 15]: the yeast constants (switches.py)
  shapes: library                                            # P8: library (v7) | pack_shapes (library/pack_shapes/pack_trees_*.npz)
```
Paths: a frozen input name, an absolute path, or relative to the build directory.
New seed streams (legacy streams untouched): `pack.round` (rounds >= 1), `pack.thinning`, `pack.thinning.pilot`,
`pack.thinning.draw`, `pack.squeeze`, `pack.organisation.field` (the random fields of the tile split).

### 3.6 `validate` (draft)

```yaml
validate:
  references: [{name: opencell, adapter: localisation_table, input: val_opencell}, {name: hpa_jurkat, input: val_hpa_jurkat}]
  compartment_groups: {nucleus: [nucleoplasm, nucleolus, speckle], cytoplasm: [cytosol, mt_lattice, mt_plus, pcm, vim_cage]}
  heldout: {hpa_jurkat_half: val_heldout_hpa_20261110}  # frozen lists only; never regenerated in a build
  metrics: [densities, crowding, loc_compare, placement, ribosome_nn]
  compare_build: {reference: "{data_root}/data/derived/wc", criteria: "{repo}/docs/reference/history/v7_reproduction_criteria.md"}
```

### 3.7 `export` (draft; contracts in `tile_format_v2.md`, `viewer_multi_cell.md`)

```yaml
export:
  tiles: {tile_nm: 200, block_tiles: 5, q_org_nm: 80, q_range_nm: 360, proxies_per_tile: 128, shard_max_mb: 45}
  fields_lod_nm: {cell: 48, ne_mid: 48, default: 32, golgi_deep: 16}
  regions: {cytosol: {colour: "#d9d2c3", group: soluble}, mt_lattice: {display_as: cytosol}}
  label_style: {}                  # overrides of core label colours
  flag_bits: []                    # cell bits 4-13, e.g. yeast histone marks
  classes: "{repo}/plugins/jurkat/data/viewer_classes.tsv"
  bookmarks: "{repo}/plugins/jurkat/data/bookmarks.yaml"
  features: {open_face: {}, golgi_cisternae: {}}
  sim_defaults: {temperature_K: 310.15, viscosity_Pa_s: 6.9e-4, ionic_strength_M: 0.15}
  model_doc: "{repo}/docs/MODEL_jurkat.md"
  licence: {data: CC-BY-4.0, attribution: "{repo}/plugins/jurkat/data/attribution.md"}
  boxes: []                        # optional simbox exports (simulation_workflows.md §3)
```

## 4. `frozen_inputs`

```yaml
frozen_inputs:
  inv_geiger_ibaq: {path: "{data_root}/data/reference/...tsv", md5: "...", origin: "Expression Atlas E-PROT-?, 2026-11-07",
                    why: "Jurkat proteome (Geiger 2012)", licence: "EMBL-EBI terms", redistributable: true}
  lib_structures_v7: {path: "{data_root}/jurkat/yc/pipeline/lib", manifest: "{data_root}/.../lib_manifest.tsv",
                      md5: "<md5 of the manifest>", origin: "v7 library", why: "frozen for exact reproduction"}
```
| Field | Rule |
|---|---|
| name | prefixed by stage: `geom_`, `inv_`, `lib_`, `pre_`, `val_` (stage.py collects `inv_*`) |
| `path` | starts with a placeholder; directories need `manifest` (a file whose checksum stands for the folder) |
| `md5` or `sha256` | required, except `unchecked: true` + `why` (development only; refused for `published`) |
| `origin`, `why` | required (provenance; replaces hand-written parameter files and orphan intermediates, H merge_notes §7.6) |
| `licence`, `redistributable` | optional; feed the deposit/Zenodo split (H merge_notes §4-5) |

Verified by `build.verify_frozen` in the `frozen` stage before any stage reads them (`frozen_inputs.tsv` in the build).

## 5. Placeholders, site file and `extends`

Placeholders resolved in every string: `{data_root}`, `{results_root}`, `{geom_root}`, `{repo}`, `{yeast_root}` (`config.py`;
`{yeast_root}` = the READ-ONLY yeast_cell_viewer tree whose inputs the yeast config reads in place, 2026-11-16).
Absolute values live only in `config/site.yaml` (env overrides `WHOLECELL_DATA_ROOT` etc., `WHOLECELL_SITE=<file>`).
Why: the human tree has the workspace path in 300+ files and in data (`libs_v7.json`), H rebuild_test §6.

`extends: <path relative to this file>`: the parent is loaded (recursively), then this file is merged on top
(`config._merge`): **maps merge key by key, lists and scalars replace**. A child must set its own `cell_id`.
Example, the proteome swap of PLANNING_NOTES 2026-11-13 (same Jurkat geometry, primary-T-cell copy numbers):

```yaml
# config/cells/jurkat_activated.yaml
extends: jurkat.yaml
cell_id: jurkat_activated
name: Jurkat geometry, activated CD4 T-cell proteome (Rieckmann 2017)
seeds: {master: 20261201}          # replaces; `legacy` dropped by writing `seeds` without it
inventory:
  adapters: {proteome: proteome_table}
  params: {proteome_table: {path: "{data_root}/data/reference/rieckmann2017_cd4_activated.tsv", unit: copies,
                            source_id: rieckmann2017_act}}
  total_protein_source: "Rieckmann 2017 activated CD4, size-scaled to the Jurkat volume"
  rule_sets: [jurkat.v4, jurkat.v6]                     # list replaces the parent's
```
Caveat: lists replace wholesale, so changing one region means repeating `pack.regions`. The stamp records the
resolved config, so the effective content is always visible (`config.resolved.json`).

## 6. Validation of a config (`wholecell.config.validate`, run by `load` and by the build stage)

| Check | When |
|---|---|
| required keys present (§2); `schema_version` supported; `cell_id` = file stem and matches `[a-z0-9_]+` | load |
| no unknown top-level keys (typos such as `inventroy` fail) | load |
| every path string starts with a placeholder; no unresolved placeholder remains | load |
| `plugins` importable; every registered name used (`adapters`, `rule_sets`, `ingredients`, `virtual_providers.plugin`, `preplaced.modules.plugin`) exists in the registry | build start |
| frozen inputs: fields of §4; names used in `params`/`inputs` are declared | build start |
| region names unique; kinds ordered; membrane regions follow engine order; `region_map` targets exist (`regions.py`) | pack stage |
| every `labels`/`membrane` named in `pack.regions` exists in the geometry vocabulary (package manifest names or aliases) | geometry stage (needs the package) |
| `preplaced.modules.depends_on_geometry` = the package id; `export.regions` keys exist in `pack.regions` | build start |
| `status: stub` cannot be built; `published` forbids `unchecked` inputs, `dev_legacy_geometry`, `allow_dirty` | build start |

## 7. Build discipline (`wholecell/build.py`, `wholecell/seeds.py`)

| Item | Rule |
|---|---|
| build id | `<YYYYMMDD>-<commit7>[-dirty]-s<master>[-<tag>]` (`make_build_id`); folder `{results_root}/<cell_id>/<build_id>/` |
| stamp | `BUILD_STAMP.json` (commit, dirty, config sha256, seeds, PYTHONHASHSEED, host, platform, python) + `config.resolved.json`; an existing stamp with another commit/config is refused (no in-place rebuild) |
| dirty tree | refused unless `build.allow_dirty: true` (recorded) |
| write-once | stage outputs via `build.write_once` / `atomic_save`; a stage writes only in its folder; no live files |
| frozen inputs | verified before use (§4) |
| seeds | one `seeds.master`; named streams `seeds.rng(cfg, "pack.counts")` / `seed_int(cfg, "pack.tile", k)` (sha256 of the stream name, never `hash()`); `seeds.legacy: {stream: int or expression in k}` overrides a stream for exact reproduction (Jurkat v7: `pack.counts: 3`, `pack.tile: "k*7919+11"`, lamina 20261112) |
| environment | `PYTHONHASHSEED=0` enforced (`check_env`); numba cache per build outside the code tree (Y merge_notes §6.2) |

## 8. Adding a new cell type (checklist)

1. **Geometry**: produce a cellgeom_v1 package with a toolkit (`toolkits/em` for a segmented volume,
   `toolkits/analytic` for primitives) under `{geom_root}/<cell>/<id>/`; run the package validator.
2. **Vocabulary**: use core label/membrane names; declare extensions (ids ≥ 32 / ≥ 1000) in the package manifest.
3. **Adapters**: generic `proteome_table` / `localisation_table` / `complex_table` with TSVs, or a plugin adapter;
   declare every input in `frozen_inputs` with checksum, origin, licence.
4. **Regions**: write `pack.regions` + `region_map` from the cell's localisation keys; set `tile_nm`, budgets.
5. **Plugins** (only if needed): `plugins/<cell>/` registering rules, preplaced modules, virtual regions, viewer
   features (cell wall, eisosomes, Rabl chromatin, microvilli...). Generic plugins may be shared by several cells.
6. **Library and preplaced**: `library.mode: build` (or frozen) and the preplaced module list.
7. **Validate/export**: references, held-out lists (frozen before the first build), colours, bookmarks, MODEL doc,
   licence.
8. Copy `config/cells/testcell.yaml` as a template; run `python -m wholecell build --cell <cell>`; nothing in
   `wholecell/` should need editing. If it does, that is a bug in the core (a hidden cell-specific constant).

## Open questions

1. `config_hash` is computed on the **resolved** config, so the same config hashes differently on two sites
   (placeholders expanded). Hash the unresolved config + site separately?
2. List merge in `extends`: keep replace-only, or add keyed merge for `pack.regions` (by `name`)?
3. Should region ids be stable across builds of one cell (append-only `pack.regions`) so viewers and ledgers can
   compare builds by id, or is name-based comparison enough (regions.json per export)?
4. Where do cell physiology constants (temperature, ionic strength) live: `export.sim_defaults` or a top-level
   `physiology` block read by pack, export and applications?
5. `library`, `preplaced`, `validate`, `export` sections are drafts; to be aligned with their stage code.
