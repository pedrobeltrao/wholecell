# wholecell

A pipeline that builds **whole-cell structural models**: every protein copy of a cell placed as a coarse 3D structure
inside a measured or analytic cell geometry, together with pre-placed large assemblies (chromatin, cytoskeleton, pores,
polysomes, ...). One shared pipeline builds different cell types from a cell config plus a cell plugin. The two cells of
the paper are a **Jurkat T cell** (human; geometry from volume EM) and a **budding yeast cell** (*S. cerevisiae*;
analytic geometry, yeast engine by G. Magnusson). A small synthetic **testcell** runs every stage in about 5 minutes and
is the CI target.

Status: **0.1.0-staging** (2026-11-23), staged from the development repository for the paper (Magnusson, Beltrao et
al., in preparation). See `docs/release_notes.md` for what was included, excluded or changed.

## Related repositories
- **wholecell** (this one): the build pipeline: stages geometry -> inventory -> library -> preplaced -> pack ->
  validate -> export.
- **wholecell-viewer**: one WebGPU viewer for the exported builds (cell picker; data-driven per
  `docs/reference/viewer_multi_cell.md`).
- **wholecell-paper**: the simulation applications run on the two paper builds (formerly `applications/` here), configs
  that pin the two builds, scripts for the main figures and tables, expected numbers.

Large data (inputs, builds, viewer tiles, simulation outputs) are deposited, not kept in git (see "Data" below).

## Install
Two conda environments (exact exports in `envs/`; details in `envs/README`), plus one for the yeast cell:

```
conda create -p <envs>/vcell --file envs/vcell.explicit.txt && <envs>/vcell/bin/pip install -r envs/vcell.pip.txt
conda env create -f envs/wholecell_smk.yaml -p <envs>/wholecell_smk
# yeast cell only: the yeast engine env (CPython 3.11; GPU with CUDA 12 for its validation step)
uv venv --python 3.11.16 <data_root>/conda_envs/yeast_x5 && uv pip install --python <data_root>/conda_envs/yeast_x5/bin/python -r envs/yeast_x5.pip.txt
```

Platform used for the paper builds: Linux x86_64 (AMD EPYC 7742), Python 3.12 (core), Snakemake 8.30.

## Configure the site
`config/site.yaml` is a template and the only file with absolute paths. Copy it **outside** the checkout, fill it in
and export `WHOLECELL_SITE` (builds refuse a dirty working tree and run from a frozen copy of the committed tree):

```
cp config/site.yaml ~/wholecell_site.yaml      # edit: data_root, results_root, geom_root, [yeast_root], [afdb_root], python, snakemake
export WHOLECELL_SITE=~/wholecell_site.yaml
```

| key | env override | needed for |
|---|---|---|
| `data_root` | `WHOLECELL_DATA_ROOT` | input data tree (Jurkat and yeast; layout in `docs/data_inputs.md`) |
| `results_root` | `WHOLECELL_RESULTS_ROOT` | where builds are written: `<results_root>/<cell>/<build_id>/` |
| `geom_root` | `WHOLECELL_GEOM_ROOT` | immutable geometry packages `<geom_root>/<cell>/<geometry_id>/` |
| `yeast_root` | `WHOLECELL_YEAST_ROOT` | optional; yeast model input tree (yeast only) |
| `afdb_root` | `WHOLECELL_AFDB_ROOT` | optional; folder with `UP000005640_9606_HUMAN_v6.tar` (Jurkat library) |
| `yeast_python` | `WHOLECELL_YEAST_PYTHON` | optional; yeast engine interpreter (default `<data_root>/conda_envs/yeast_x5/bin/python`) |

## Build the testcell (about 5 min on 8 cores)
```
PYTHONHASHSEED=0 <envs>/vcell/bin/python -m toolkits.analytic.generate toolkits/analytic/specs/testcell.yaml   # geometry, once
scripts/build_isolated.sh testcell <tag> 8      # prints the build folder; Snakemake runs detached
scripts/ci_testcell.sh                          # optional CI: builds with 8 and 4 workers must be identical
```
Walk-through and what to read in the outputs: `docs/getting_started.md`.

## Build the two paper cells
Step-by-step guide for outside users (inputs, environments, site file, commands, resources, checks): `docs/REBUILD.md`.
Inputs first (see "Data"); every declared input is checksum-verified by the first stage (`frozen_inputs.tsv` in the
build). Pins of the published builds (build ids, tags, manifests, config hashes): `config/paper_builds.yaml`.

```
export WHOLECELL_SITE=~/wholecell_site.yaml
scripts/build_isolated.sh jurkat paper 24       # -> <results_root>/jurkat/<date>-<commit7>-s20261115-paper
scripts/build_isolated.sh yeast  paper 24       # -> <results_root>/yeast/<date>-<commit7>-s0-paper
snakemake -n --config cell=jurkat               # dry run of the DAG (orchestration env)
```

Expected resources (paper builds, 24 cores on one node; stage times from the build logs):

| cell | wall time | notes |
|---|---|---|
| Jurkat | about 15 h | CPU only. Inventory ~10 min, library ~5.5 h, pre-placed ~1.2 h, pack setup ~0.7 h, pack ~5 h, validate ~0.8 h, export ~1.5 h. ~940 M placed instances. |
| yeast | about 7 h | Pack ~2 h; the vendored yeast validation (`validate_all`) needs 1 NVIDIA GPU and takes ~2.7 h of it. ~24 M placed instances. |

The published yeast build (`paper-yeast-u1`) is a downstream-only rerun (validate + export at commit 0a3ba79) that
reuses the upstream stages of `yeast/20261119-bedcea9-s0-paper`: `scripts/build_downstream_reuse.py --donor <build>`.
Freeze verdicts against the frozen criteria: `scripts/freeze_check.py` (criteria `docs/unification_criteria/freeze_v8.md`
and the two acceptance addenda).

## Add a cell
A new cell type = a cell config (`config/cells/<cell>.yaml`), a geometry package and, if needed, a plugin
(`plugins/<cell>/`). Worked example with the testcell: `docs/adding_a_cell.md`; stage reference: `docs/stages/`;
config spec: `docs/reference/cell_config.md`.

## Data
`docs/data_inputs.md` lists every input of the two paper builds (relative to `data_root`, `yeast_root`, `afdb_root`,
`geom_root`), with size, checksum, original source, licence and where to get it; machine-readable: `data_manifest.tsv`.
- Redistributable inputs and our derived intermediate files: Zenodo (https://zenodo.org/records/TODO-wholecell-inputs).
- Build outputs of the two paper builds: https://huggingface.co/datasets/TODO/wholecell-builds.
- Lab SEC-MS 2026 and Jurkat surface glyco-capture: processed tables with the inputs; raw MS files at PRIDE (TODO).
- Not redistributable (fetch from the source): CORUM, journal supplementary tables, HPA images, the unpublished yeast
  AlphaFold 3 pool predictions.

## Licence
- **Code**: PolyForm Noncommercial 1.0.0 (`LICENSE`).
- **Everything else we generated** (configs, curated tables, docs, figures, model outputs and deposits): CC BY-NC 4.0
  (`LICENSE-DATA`; files covered listed in `NOTICE.md`).
- **Commercial use** needs a separate licence: contact Pedro Beltrao (pbeltrao@ethz.ch); commercial licences are
  arranged through ETH Zurich technology transfer (`COMMERCIAL.md`).
- Third-party code and data keep their own terms (`NOTICE.md`). The `Required Notice` wording in `LICENSE` is a
  placeholder to be confirmed.

## Citation
See `CITATION.cff` (paper in preparation: Magnusson, Beltrao et al.).
