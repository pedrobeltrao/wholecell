# Rebuilding the two paper cells from this repository

This page lists the steps an outside user follows to rebuild the Jurkat T cell and the budding yeast cell of the paper
(Magnusson, Beltrao et al., in preparation) with this pipeline, and to check the result against our published builds.
Background: the input list is in [data_inputs.md](data_inputs.md) and `data_manifest.tsv`; the build stages are in
[stages/](stages/); the published builds are pinned in `config/paper_builds.yaml`.

**Inputs dataset.** The inputs are deposited on Hugging Face as the dataset `pedrobeltrao/wholecell-inputs`
(https://huggingface.co/datasets/pedrobeltrao/wholecell-inputs). `scripts/fetch_inputs.py` downloads it into the roots
of your site file. The build outputs go to a separate builds dataset (not yet deposited).

## 1. What you need

| | Jurkat (`config/cells/jurkat.yaml`) | yeast (`config/cells/yeast.yaml`) | testcell |
|---|---|---|---|
| CPU | one node, 24 cores (paper builds: AMD EPYC 7742) | one node, 24 cores | 8 cores |
| GPU | none | 1 NVIDIA GPU with CUDA 12, for the vendored validation step `validate_all` only | none |
| wall time | about 12.5-15 h at 24 cores (stage times below) | about 7-8 h (packing about 2.2 h; `validate_all` about 2.7 h on the GPU) | about 6 min, plus 2 min for its geometry |
| disk for the build | about 50 GB (the v8 build is 49 GB) | about 9 GB (the u1 donor build is 8.3 GB) | < 1 GB |
| inputs | about 22 GB under `data_root` (of which 15.3 GB is the local PDB assembly set; only the files listed in `<build>/library/work/pdb_assemblies_read.tsv` are read), 5.2 GB AlphaFold DB human tar, 1.4 GB geometry package | about 4 GB (yeast data subset 1.1 GB, AlphaFold DB yeast tar 1.0 GB, AFDB dimers 1.0 GB, PDB cache 0.6 GB, work seed and start library), 87 MB geometry package | none (in the repository) |

Jurkat stage times of the paper build (24 cores): inventory about 10 min, library about 5.5 h, pre-placed structures
about 1.2 h, pack setup about 0.7 h, packing about 5 h, validation about 0.8 h, export about 1.5 h. Peak memory was not
recorded; the paper builds ran on nodes with 512 GB of RAM, and we have not tested smaller machines.

## 2. Inputs

Every input of a cell is declared in its config under `frozen_inputs` (path relative to a site root, md5 or sha256).
`data_manifest.tsv` has one row per input with its source, licence and where to get it. There are four kinds:

1. **Third-party public data** (RNA-seq tables, Geiger 2012, UniProt, HPA tables, Complex Portal, SIFTS, PDB, STRING,
   AlphaFold DB, PaxDb, PRIDE SEC-MS, journal supplements). Where the licence allows it, the exact files the paper
   builds read are in the inputs dataset (their own licence, in its `cc_by/`, `cc0/`, `open/`, `cc_by_sa/` folders):
   exports change over time (UniProt REST answers, AFDB dimer fetches, the RCSB cache), so a new download would not
   match the checksums. Not in the dataset (its `EXCLUDED.tsv`): the AlphaFold DB v6 proteome tars (download from EBI)
   and the third-party tables whose terms do not clearly allow redistribution (fetch them from the journal or database,
   into the same path).
2. **Curated tables in this repository** (`plugins/*/data/`, `config/geometry/`, `tests/data/`): nothing to download.
3. **Deposited derived data and lab data**: the geometry packages (Jurkat from the jrc_jurkat-1 FIB-SEM segmentations,
   yeast analytic), the chromatin model (GPU molecular dynamics, not bitwise reproducible), the U-Net microtubule traces,
   the hold-out gene split, the lab SEC-MS and surface glyco-capture tables (raw files at PRIDE), the yeast start library
   (Boltz-2 tRNA models, builder templates) and work seed. In the inputs dataset (`ours/`, CC BY-NC 4.0).
4. **Inputs derived from earlier builds** (fixed inputs, in the inputs dataset `ours/` with their provenance). The pipeline cannot
   regenerate these from public data in one build, because they were made from the packs of earlier model versions:
   - `inv_v4pack_loc_compare`, `inv_v5pack_loc_compare`: per-gene placed compartment shares of the legacy v4 and v5
     Jurkat packs, read by the localisation rules v5 (C, N, P) and v6 (R2-R6a);
   - `inv_tier2_v5`, `inv_tier2_transfer`: imaging metrics of the v5 pack against HPA Jurkat images (rule v6
     calibration);
   - `pre_cytoskeleton_v7`: the centriole pair and the cortical F-actin of the v7 cytoskeleton (frozen obstacles; there
     is no core builder for the cortex and centrioles yet), with its summaries `inv_cytoskeleton_v7_json`,
     `inv_cytoskeleton_v6_json`;
   - `pre_mrna_v1`: the legacy mRNA beads (transcripts threaded by the polysome builder, and the nuclear mRNAs);
   - the legacy v7 structure library (`lib_structures_v7` / `inv_struct_library`) and its derived tables
     (`inv_membrane_lib_ids`, `inv_tm_route_libs`, `lib_afdb_fragment_index_tmroute`);
   - order and tie-break references `inv_legacy_ingredients_v7`, `inv_legacy_complexes_v0` (they fix row order and
     equal-coverage PDB picks, which the packer's random draws depend on);
   - yeast: `yeast_x5_lib_start` and `yeast_x5_work_seed` (from the earlier X5 run).

   Since v9 the Jurkat lamina consumption (rule L3) comes from the lamina this build makes (it used to be read from the
   legacy v7 lattice).

Optional inputs (a build runs without them): the archived v7 reference for `validate compare_v7` (skipped with a reason
when absent), the suite's no-regression reference builds (rows marked n/a), the legacy inventory tables of
`inventory.verify` (reported missing), and for yeast the 9.7 GB YeastRGB raw image archive
`validation/microscopy/yeastrgb_dubreuil2019/raw/RAW-C_SWAT_GFP-ORF.zip` (download it from YeastRGB,
Dubreuil et al. 2019 Nucleic Acids Res 47:D1245, into `{yeast_root}/yeast_cell_dev_bundle/pipeline/data/` at that path). The
yeast image N/C validation (`image_nc`, suite check C3) uses the per-ORF measurements cached in the work seed, so it runs
without the archive; `image_nc.json` then records `status: ok_from_cache_raw_archive_absent` and the log a WARNING. If
neither the archive nor the cache exists, the test does not run, `status: not_run_raw_archive_absent` is recorded and
the suite's C3 row carries the warning.

**Yeast data tree.** The yeast cell reads the yeast model's `pipeline/data` folder (`yeast_pipeline_data`). Only the
782 files matched by `plugins/yeast/data/pipeline_data_read_set.txt` (about 1.1 GB) are read; the rest of the 14.5 GB
tree (literature PDFs, raw microscopy, genome tracks) is not needed. The setup step checks the files against
the md5 manifest `{data_root}/data/derived/wholecell_frozen/yeast_x5/manifest_pipeline_data.tsv` (in the inputs
dataset). Files outside the read set may be absent. A file whose modification time differs from the manifest (as in any download) is
accepted when its md5 matches.

**Layout and download.** Every input sits under its root (`data_root`, `yeast_root`, `afdb_root`, `geom_root`) at the
path of `data_manifest.tsv`. With the site file in place:

```
python scripts/fetch_inputs.py --cell jurkat --list        # what comes from the dataset, what from elsewhere
python scripts/fetch_inputs.py --cell jurkat --source      # download (+ the AlphaFold DB tar from EBI)
python scripts/fetch_inputs.py --cell yeast --source
python scripts/fetch_inputs.py --cell yeast --group reference   # also the reference presets' inputs (optional)
```

The script checks every file against its md5 and unpacks the tar shards of large folders (the PDB assembly set, the
v7 structure library, the AFDB dimer models). It ends with the list of inputs to fetch by hand (third-party tables
whose terms do not allow redistribution), each with its source. Unpacked shards are recorded in
`<root>/.fetch_markers/`, outside the data folders (the folder listings are checked: the PDB assembly listing md5 in
`frozen` and in `inventory_assembly`); markers that an older copy of the script left inside a data folder
(`.fetched.*`) are moved there on the next run.

## 3. Software

```
conda create -p <envs>/vcell --file envs/vcell.explicit.txt && <envs>/vcell/bin/pip install -r envs/vcell.pip.txt
conda env create -f envs/wholecell_smk.yaml -p <envs>/wholecell_smk
# yeast only: the yeast engine environment (CPython 3.11; CUDA 12 for validate_all)
uv venv --python 3.11.16 <envs>/yeast_x5 && uv pip install --python <envs>/yeast_x5/bin/python -r envs/yeast_x5.pip.txt
```

Platform of the paper builds: Linux x86_64, Python 3.12 (core), Snakemake 8.30. The exact package versions are in
`envs/`.

## 4. Site file

Copy `config/site.yaml` outside the checkout, fill it in and export `WHOLECELL_SITE`:

```
data_root: /data/wholecell_inputs
results_root: /data/wholecell_builds
geom_root: /data/wholecell_geom
yeast_root: /data/yeast_cell_model           # yeast only
afdb_root: /data/alphafold_v6                # Jurkat only: holds UP000005640_9606_HUMAN_v6.tar
yeast_python: <envs>/yeast_x5/bin/python     # yeast only (default <data_root>/conda_envs/yeast_x5/bin/python)
python: <envs>/vcell/bin/python
snakemake: <envs>/wholecell_smk/bin/snakemake
max_workers: 24
```

Builds refuse a dirty working tree and run from a frozen `git worktree` of HEAD, so keep the site file outside the
checkout.

## 5. Commands

```
export WHOLECELL_SITE=/path/to/site.yaml
PY=<envs>/vcell/bin/python

# 1. testcell (no downloads): geometry once, then the build (about 6 min)
PYTHONHASHSEED=0 $PY -m toolkits.analytic.generate toolkits/analytic/specs/testcell.yaml
scripts/build_isolated.sh testcell mytest 8

# 2. check the inputs of a paper cell before a long build: writes <dir>/frozen_inputs.tsv and fails, naming every
#    input that is missing or has another checksum
PYTHONHASHSEED=0 $PY -m wholecell frozen --cell jurkat --build /tmp/check_jurkat
PYTHONHASHSEED=0 $PY -m wholecell frozen --cell yeast --build /tmp/check_yeast

# 3. dry run of the workflow (Jurkat 38 jobs, yeast 40); check that the exit code is 0
<envs>/wholecell_smk/bin/snakemake -n --cores 24 --config cell=jurkat build_tag=dry
<envs>/wholecell_smk/bin/snakemake -n --cores 24 --config cell=yeast build_tag=dry

# 4. the builds (each prints its build folder and runs detached; progress in <build>/snakemake.log)
scripts/build_isolated.sh jurkat rebuild 24
scripts/build_isolated.sh yeast rebuild 24
```

To resume a build that stopped, run the Snakemake command recorded in `<build>/launch.log` again from the build folder
(the launcher refuses an existing build folder).

The testcell geometry you generate does not match the pinned manifest checksum, because the manifest records the
generating commit and time and the last bits of two files depend on the code version. The testcell config sets
`geometry.on_mismatch: warn`, so the build records the mismatch in `geometry/package.json` and goes on. A testcell built
this way from a fresh clone was content-identical to our CI build (same tiles, counts and validation results; only
paths and timings differ). The paper cells keep the strict check, so use the deposited geometry packages for them.

## 6. Checking a rebuild against the published builds

1. **Inputs.** Compare your `<build>/frozen_inputs.tsv` with the one in the published build: same names and checksums.
2. **Code and configuration.** `<build>/BUILD_STAMP.json` has the commit and `config_sha256`; compare them with
   `config/paper_builds.yaml` (or the published build's stamp). A different commit or configuration explains any
   difference that follows.
3. **Content.** `python scripts/ci_compare_builds.py <your build> <published build> --out cmp.json` compares the two
   folders by content and ignores what legitimately differs (build id, paths, wall-clock fields, logs). Builds of the
   same commit and configuration with different worker counts are identical (testcell CI). We have not tested other
   CPUs or library versions; floating-point differences there could change the packing.
4. **Numbers.** Compare `pack/summary.json` (placed instances, placed mass fraction) and `validate/suite.json` (the
   validation checks C1-C17) with the published build.

## 7. Known limits

- The inputs of item 4 in section 2 come from earlier models; rebuilding them needs the legacy pipeline (in the inputs
  dataset, group `reference`, as `jurkat/yc/pipeline`) and the earlier packs. They are fixed inputs of this pipeline, not outputs.
- The yeast code falls back to downloading from RCSB when an entry is missing from the deposited PDB cache or
  `rcsb_cache.json`. With the deposited cache no download happens. Without it, the result depends on the RCSB state on
  the day of the build.
- If `genome/nucleosome_dyads.tsv.gz` is missing from the yeast data tree, the yeast engine rebuilds it and writes it
  into that tree. It is in the read set, so keep it.
- `inputs/yeast_organelle_meshes.zip` (yeast) is read by one non-fatal geometry report script. It is not a declared
  input (it is in the inputs dataset, group `optional`).
- The reference preset `yeast_x5_core.yaml` reads the per-pair PAE table of the lab's pooled AlphaFold 3 screen
  (`yeast_af3_pool_pairs`; source https://github.com/jurgjn/pooled-ppi-yeast). The public release v26.1 gives the same
  value for 7,082,835 of its 7,631,089 pairs but covers fewer pairs (10,264,748 in the table the build read), so the
  table itself is in the inputs dataset (`af3_pool/`, AlphaFold 3 output terms). `yeast.yaml` does not read it.
