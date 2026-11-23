# Overview

`wholecell` builds a molecular model of one whole cell: every protein complex, RNA and pre-placed structure is given a
position and orientation inside a 3D cell geometry. One pipeline serves several cell types (today: Jurkat T cell,
budding yeast, and a small synthetic `testcell`).

## What a build is

A **build** is one run of the pipeline for one cell. It lives in `<results_root>/<cell_id>/<build_id>/`
(`results_root` is set in `config/site.yaml`).

- **Build id**: `<date>-<commit7>-s<seed>[-<tag>]`, e.g. `20261115-ec0b403-s1-ci2` (testcell, commit ec0b403,
  master seed 1, tag `ci2`). It is computed by `python -m wholecell.plan`.
- **Write-once stage folders**: each stage writes only inside its own folder (`geometry/`, `inventory/`, ...) and
  never overwrites a file. A change means a new build id; old builds are never edited or deleted.
- **Frozen worktree**: `scripts/build_isolated.sh` refuses a dirty working tree, then runs the build from a detached
  `git worktree` of HEAD in `<build>/_code`. Later commits cannot change a running or finished build.
- **`BUILD_STAMP.json`** (build root): commit, dirty flag, config path and sha256, seeds, `PYTHONHASHSEED` (always 0),
  host, platform, Python version, start time. Next to it: `config.resolved.json` and `frozen_inputs.tsv` (every
  declared external input with its checksum, verified before any stage reads it).
- **Seeds**: every random step draws from a named stream of the master seed (`wholecell/seeds.py`). Jurkat also
  declares its legacy v7 seeds (`seeds.legacy`) to reproduce v7 exactly.

## A cell is four things

| Part | Where | Spec |
|---|---|---|
| config | `config/cells/<cell>.yaml` | [reference/cell_config.md](reference/cell_config.md) |
| geometry package | `{geom_root}/<cell>/<geometry_id>/`, immutable, made by a toolkit (`toolkits/em`, `toolkits/analytic`) | [reference/cellgeom_v1.md](reference/cellgeom_v1.md) |
| data adapters | `wholecell/inventory/adapters/*` (generic) or plugins | [reference/inventory_schema.md](reference/inventory_schema.md) |
| plugins | `plugins/<name>/`, listed in config `plugins:` | registered by name (`wholecell/registry.py`) |

No core module contains a list of cells, compartments or species. The cell is chosen only by
`--config cell=<name>`.

## Stages, in the order of the paper's Methods

| # | Stage | What it does | Today's folder / code |
|---|---|---|---|
| 1 | [geometry](stages/geometry.md) | load and check the cell geometry package; volumes and areas | `geometry/`, `wholecell/geometry` |
| 2 | [composition](stages/composition.md) | copy number of every protein | `inventory/` (steps `rna`, `copies`, `v0`) |
| 3 | [localisation](stages/localisation.md) | which compartment(s) each protein is in | `inventory/` (adapter + `rules/<set>/`) |
| 4 | [assembly](stages/assembly.md) | which proteins are in complexes, and how many copies | `inventory/` (steps `secms`, `assembly`, `ingredients`) |
| 5 | [structures](stages/structures.md) | one 3D structure per species (the library) | `library/` |
| 6 | [preplaced](stages/preplaced.md) | large structures placed before packing (chromatin, MTs, NPCs, ...) | `preplaced/`, `wholecell/preplaced` |
| 7 | [packing](stages/packing.md) | place every remaining copy, tile by tile | `pack/` |
| 8 | [validation](stages/validation.md) | compare the model with independent data | `validate/` |
| 9 | [export](stages/export.md) | viewer files, tiles, simulation boxes | `export/` |

Simulations (Brownian dynamics, synapse, condensates, ...) are separate **applications** that read a finished,
pinned build. They are not stages of the build. See [simulations.md](simulations.md).

## Current code vs plan

The stage names above are the target names from the unification plan
([reference/history/unification_plan.md](reference/history/unification_plan.md) §8). Today:

- One `inventory` stage covers composition, localisation and assembly. Its steps are
  `rna -> copies -> v0 -> secms -> assembly -> ingredients -> rules:<set> -> pack_table -> canonical`. The split into
  separate stages comes with plan steps A4 (assembly) and A5 (localisation).
- The pre-placed stage was called `prepass` until A0c (2026-11-17); folder, CLI stage, Snakemake rule file, registry
  kind, config section and file names are now `preplaced`. Builds made before the rename (the reference builds) keep
  `prepass/` and `pack/prepass_*.npy`; code that reads them goes through the alias map `LEGACY_NAMES` in
  `wholecell/build.py` (`ref_path`), and old configs and `python -m wholecell prepass` still work with a message.
- Yeast builds through the core stages with the `legacy_x5` presets (A9); yeast-only steps are plugin steps, some of
  which still call vendored yeast code (`plugins/yeast/engine`, `plugins/yeast/steps.yaml`). The whole-engine X5
  reproduction (`pipeline_engine: yeast_x5`) is kept only in `config/cells/yeast_x5_port.yaml`.
- Jurkat and yeast use one packer (A7; yeast with `pack.preset: legacy_x5`) and one geometry evaluator (A1);
  common defaults follow in B1.

`PORTING_STATUS.md` (repository root) lists what is ported and run, ported but not run, or still to do.
