# wholecell architecture (phase 0-1)

Start with [README.md](README.md) for the reading order (overview, getting started, stages, adding a cell).

Shared pipeline for yeast, Jurkat and further cell types. Justification against the two audits is in each spec
(`docs/reference/`); this page fixes the module boundaries and contracts every component codes against.

## 1. A cell = four things, nothing else

| Part | Where | Contract |
|---|---|---|
| config | `config/cells/<cell>.yaml` | `docs/reference/cell_config.md`; loaded by `wholecell.config.load` |
| geometry package | `{geom_root}/<cell>/<geometry_id>/` (immutable, `cellgeom_v1`) | `docs/reference/cellgeom_v1.md`; produced by a toolkit (`toolkits/em`, `toolkits/analytic`) |
| data adapters | `wholecell/inventory/adapters/*` (generic) or plugins | registered `adapter.<kind>`; return tables in the inventory schema |
| plugins | `plugins/<name>/` listed in config `plugins:` | register named components (`wholecell.registry`) |

No module may contain a list of cells, compartments, labels, regions or species. The cell is chosen only by
`--cell <name>` (CLI) / `--config cell=<name>` (Snakemake). Constants that differ between cells are config values.

## 2. Repository layout

```
Snakefile                     entry: snakemake --config cell=<cell> [build_tag=...]
workflow/rules/*.smk          one file per stage; rules call `python -m wholecell <stage> --cell C --build B --step S`
config/site.yaml              site template: the only file with absolute paths (data_root, results_root, geom_root, python)
config/cells/*.yaml           jurkat.yaml (ported v7), testcell.yaml (synthetic CI cell), yeast.yaml (stub)
envs/                         conda specs per step family (yaml + `conda list --explicit` export)
wholecell/                    core package (stages below), config.py, build.py, seeds.py, registry.py, __main__.py
plugins/<cell or feature>/    cell-specific components registered by name (jurkat: rules v4-v7, virtual regions...)
toolkits/em/                  EM -> cellgeom_v1 generator (Jurkat: conversion of the v7 geometry + ported generator code)
toolkits/analytic/            analytic -> cellgeom_v1 generator (testcell; yeast-style primitive stub)
tests/                        unit tests + determinism test on testcell
docs/                         README.md (index), overview, getting_started, stages/, adding_a_cell, simulations
docs/reference/               the six specs; history/ = reproduction criteria, reports, unification plan (frozen)
```

## 3. Stages and their folders inside one build

`<results_root>/<cell_id>/<build_id>/` (build discipline: `wholecell/build.py`):

| Stage folder | Module | Reads | Writes |
|---|---|---|---|
| (root) | `build` | config | `BUILD_STAMP.json`, `config.resolved.json`, `frozen_inputs.tsv` |
| `geometry/` | `wholecell.geometry` | package path from config | `validation.json`, `measures.json` (volumes/areas), link to the package |
| `inventory/` | `wholecell.inventory` | adapters, frozen tables, geometry measures | `inventory_v0.tsv`, `assembly/*.tsv`, `rules/<step>/{ingredients.tsv,changes.tsv}`, `species.tsv`, `species_loc.tsv`, `field_provenance.tsv`, `pack_table.tsv` |
| `library/` | `wholecell.library` | frozen library (Jurkat) or library builder | `library_manifest.tsv` (species, file, sha256), `libs.json` |
| `preplaced/` | `wholecell.preplaced` | frozen pre-placed files or builders (registered `preplaced` modules) | per-module npz/json + `preplaced_manifest.json` |
| `pack/` | `wholecell.pack` | geometry, pack_table, library, preplaced | `tiles.parquet`, `tiles_meta.json`, `counts.npz`, `libs.json`, `preplaced_*.npy`, `tiles/<i>.npz`, `pack_log.tsv`, `summary.json` |
| `validate/` | `wholecell.validate` | pack outputs + references | `densities.json`, `crowding.json`, `loc_compare.tsv`/`.json`, `placement.json`, `MANIFEST.tsv` |
| `export/` | `wholecell.export` | everything above | `regions.json`, `cell.json` (cells.json entry), tile v2 data (optional target), sim manifest |

Rules: a stage writes only inside its folder; outputs are write-once (`build.write_once`, `atomic_save`); no
"live" files; every stochastic step draws from `wholecell.seeds` (named streams of the master seed; legacy seeds
only when declared in config for exact reproduction); `PYTHONHASHSEED=0` is enforced; frozen inputs are declared in
config with checksums and verified by the `frozen` stage before anything reads them; pack-dependent rules
(Jurkat v5 P, v6 R1) read only declared frozen tables.

## 4. Geometry contract (packer side)

- `wholecell.geometry.cellgeom.load(package_dir) -> Cell` with `.G` (float64 vector), `.pores` (n,6), `.pidx/.PT/.PQ/.PC`
  (empty prim arrays kept for the kernel signature), `.labels(P)`, `.mems(P)`, `.manifest`, `.vocab`.
- The package **layout** is named in the manifest (`engine: sdf_sections_v1` for Jurkat and testcell); the core evaluates
  it with the table-driven evaluator `cellgeom_eval_v1` (config `geometry.evaluator`) and the layout's rule table
  `config/geometry/rules/<layout>.yaml` (config `geometry.rules`; `cellgeom_v1.md` §4.3). An engine
  is a module exporting the numba interface the packer uses: `mem_sdfs, label_from, label_of, min_mem_dist, in_pore,
  face_dist, outer_depth, label_points, mem_points, field_sdf, n_fields` and the scalars `NMEM, HALF, M_PM, M_NE,
  M_OUTER, CYTO, PERI_W, PORE_R, MEM_HALF` plus engine label/membrane codes. `wholecell.geometry.engine` re-exports
  the one core evaluator (`cellgeom_eval_v1`; no environment switch) so `pack3d` can `from wholecell.geometry
  import engine as c3`; `sdf_sections_v1` stays importable for the legacy-reproduction test (B1).
- `sdf_sections_v1` = the Jurkat `cell3d.py` evaluator with every physical constant (PERI_W, IMS_W, PORE_R, MEM_HALF,
  cristae/Golgi constants) moved out of the code into reserved header slots of G, written from the manifest. It must
  reproduce the legacy evaluator exactly (labels identical, membrane SDFs |diff| < 1e-9 nm on >= 1e6 points).
- Engine label codes are compact (0..63, packer bitmasks are int64). The manifest maps engine codes to vocabulary
  names and to cellgeom ids (shared core ids < 1000, cell extensions >= 1000; no cap). Regions refer to names only.

## 5. Inventory contract (pack side)

The pack stage reads `inventory/pack_table.tsv`, the legacy-compatible view (columns of `ingredients_v7.tsv`:
`ing_id, kind, genes, stoich, copies, mass_da, mass_pg, comp_weights (JSON), transmembrane, n_tm, struct_src, uniprot,
signal_frac, mito_om_frac`). The canonical tables (`species.tsv`, `species_loc.tsv` long form with `rule_id`,
`field_provenance.tsv`, SEC-MS and surface-proteomics slots) follow `docs/reference/inventory_schema.md`. For Jurkat the
pack table must be byte-identical to the legacy `ingredients_v7.tsv` (md5 d1a0099caad792d0ba42cbc52b52f0a5).

## 6. Regions (pack)

Regions are config data (`pack.regions`, ordered; index = job_region id): `volume` (labels), `membrane` (engine
membrane), `shell` (membrane, labels, side, d0, d1 or `anchor`/`mem_half+X`, add-radius flag, pore-axis limit),
`virtual` (registered `virtual_region` plugin + params, packed with the shells). The localisation-key map
`pack.region_map: {key: {soluble: region, tm: region}}` replaces `jpack_test.SOL/TMR`; generic mapping rules
(signal-peptide lumen rule, OM vs IM, area factors) are core with parameters in config.

## 7. Seeds for the Jurkat reproduction

`config/cells/jurkat.yaml` declares `seeds.legacy` for the published v7 seeds (counts 3, tile `k*7919+11`, lamina
20261112, ...). With those, identical inputs and the same numba, the port is expected to be bit-identical to v7; a
master-seed build (no legacy block) measures the seed spread.
