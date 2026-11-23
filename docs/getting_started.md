# Getting started: build the testcell

The `testcell` is a small synthetic cell (analytic geometry, about 5.8 µm³; 45 structures copied from the Jurkat v7
library; synthetic proteome). It runs every stage on CPU and is the CI target. Its config is
`config/cells/testcell.yaml`.

## 1. Environments

Two conda environments are used (details and recreate commands in `envs/README`):

| Env | Used for | Files |
|---|---|---|
| core-cpu (`vcell`) | every `python -m wholecell <stage>` step | `envs/vcell.explicit.txt`, `envs/vcell.pip.txt` |
| orchestration (`wholecell_smk`, Snakemake 8.30) | the workflow only | `envs/wholecell_smk.yaml` |

```
conda create -p <path>/vcell --file envs/vcell.explicit.txt && <path>/vcell/bin/pip install -r envs/vcell.pip.txt
conda env create -f envs/wholecell_smk.yaml -p <path>/wholecell_smk
```

Then make a site file. `config/site.yaml` is a template and the only file with absolute paths: `data_root`,
`results_root`, `geom_root`, optional `yeast_root` and `afdb_root`, `python` (the vcell interpreter), `snakemake`,
`max_workers`. Copy it outside the checkout, fill it in and point `WHOLECELL_SITE=<path>` at it (builds refuse a dirty
working tree, and isolated builds run from a frozen copy of the committed tree, so an edited `config/site.yaml` would
not be seen). Single keys can also be set with `WHOLECELL_DATA_ROOT`, `WHOLECELL_RESULTS_ROOT`, `WHOLECELL_GEOM_ROOT`,
`WHOLECELL_YEAST_ROOT`, `WHOLECELL_AFDB_ROOT`. The testcell needs only `results_root`, `geom_root` and the two
interpreters.

The testcell geometry package must exist under `{geom_root}/testcell/analytic_v1`. If it does not, generate it
(about 2 min). Generation is deterministic for one code version and machine, but the manifest records the generating
commit and time, so a package you generate does not match the pinned `geometry.manifest_sha256`; the testcell config
sets `geometry.on_mismatch: warn`, so the build goes on and records the mismatch in `geometry/package.json`:

```
PYTHONHASHSEED=0 <vcell python> -m toolkits.analytic.generate toolkits/analytic/specs/testcell.yaml
```

## 2. Build

Recommended: an isolated build from a frozen worktree of HEAD (the tree must be committed):

```
scripts/build_isolated.sh testcell <tag> 8
```

It prints the build folder and returns at once; Snakemake runs detached. Progress is in `<build>/snakemake.log`, the
launch record in `<build>/launch.log`.

Plain form, from the repository checkout (uses the checkout itself, not a frozen copy):

```
snakemake --cores 8 --config cell=testcell [build_tag=<tag>]
```

Both refuse a dirty working tree (unless the config sets `build.allow_dirty: true`, which is recorded in the stamp).
Build ids are `<date>-<commit7>-s<seed>[-<tag>]`, so a second build of the same commit on the same day needs a new
tag.

**Time.** About 5 min on 8 cores: build `testcell/20261115-ec0b403-s1-ci2` ran from 16:41:57 (stamp) to 16:46:51
(last Snakemake job), 25 jobs.

**CI.** `scripts/ci_testcell.sh` is the testcell CI check (full build with 8 and 4 workers, outputs must be
identical). With `git config core.hooksPath scripts/hooks` it starts in the background after each commit and each merge
on `main` (hooks `post-commit`, `post-merge`), one run at a time; a commit that already has a CI result is skipped.

## 3. Where outputs land

`<results_root>/testcell/<build_id>/`:

```
BUILD_STAMP.json  config.resolved.json  frozen_inputs.tsv  launch.log  snakemake.log  logs/  _code/
geometry/   package.json, validation.json, measures.json
inventory/  v0/, assembly/, ingredients/, rules/<set>/, secms_*.tsv, pack_table.tsv, species.tsv, species_loc.tsv, ...
library/    library_manifest.tsv, library_meta.json
preplaced/  preplaced_manifest.json        (no pre-placed structures in the testcell; `prepass/` before A0c)
pack/       tiles.parquet, tiles/<i>.npz, counts.npz, preplaced_*.npy, pack_log.tsv, summary.json
validate/   densities.json, crowding.json, placement.json, MANIFEST.tsv, <step>.skipped.json
export/     regions.json, cell.json, web/ (YCT2 tiles), boxes/
```

## 4. What to read

1. `pack/summary.json`: placed vs expected mass. For `ci2`: 2,920,324 instances, mass fraction 0.979.
2. `validate/placement.json`: clashes and containment on sampled tiles (ci2: 104 tiles, clash rate 0 per 10⁶).
3. `validate/densities.json`: protein density per region in mg/ml (ci2: cytosol 84.4, nucleoplasm 81.7).
4. `validate/crowding.json`: crowding metrics in sampled boxes.
5. `validate/*.skipped.json`: checks that need reference data the testcell does not have (`loc_compare`,
   `ribosome_nn`, `chromatin`). Skipped is expected here.

Each stage is described in [stages/](stages/). A combined `validation/suite.json` and opening the testcell in the
viewer are planned (plan steps A8 and the viewer work), not available yet.
