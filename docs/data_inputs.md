# Data inputs of the paper builds

This page lists the input data that the two paper builds read: Jurkat `paper-jurkat-v8` and yeast
`paper-yeast-u1`, with the changes of the v9 / u2 configs noted (step-by-step rebuild: [REBUILD.md](REBUILD.md)). It also says where to get each input. The machine-readable list is
[`data_manifest.tsv`](../data_manifest.tsv), with one row per file or directory. That file gives the root, the relative
path, the size, the checksum, the role, the original source, the URL, the licence, whether it may be redistributed, and
where it will be deposited.

Deposit locations (placeholders until the deposits exist; nothing is uploaded yet):

| Location | What goes there |
|---|---|
| Hugging Face dataset `pedrobeltrao/wholecell-inputs` | the inputs of both cells: our derived data and lab tables (`ours/`, CC BY-NC 4.0), dated snapshots and copies of redistributable third-party data in their own licence folders (`cc_by/`, `cc0/`, `open/`, `cc_by_sa/`), the AF3 pool pairs of the yeast reference preset (`af3_pool/`). `scripts/fetch_inputs.py` downloads it into the site roots; staging: `scripts/deposit/stage_hf.py`, upload: `scripts/deposit/upload_hf.py` |
| `https://huggingface.co/datasets/TODO/wholecell-builds` | build outputs: the paper builds, the yeast donor pack, validation reference builds |
| `PRIDE TODO` | raw mass-spectrometry files of the lab SEC-MS 2026 runs (Jurkat NT / TCR-activated 24 h, yeast NaCl / YPD) and of the Jurkat surface glyco-capture time course |
| `fetch from source` | inputs not in the dataset (its `EXCLUDED.tsv`): the AlphaFold DB proteome tars and third-party tables whose terms do not clearly allow redistribution; `scripts/fetch_inputs.py` lists them with their source |

## 1. Roots

Every input path in a cell config is written relative to a configurable root, as a `{placeholder}`. The roots are set in
`config/site.yaml` or with environment variables:

| Placeholder | Site key / env | Content |
|---|---|---|
| `{data_root}` | `data_root` / `WHOLECELL_DATA_ROOT` | the input tree: `data/` (proteomics, transcriptomics, structures, reference, literature, derived) and `jurkat/yc/pipeline/` (legacy Jurkat structure library and tables) |
| `{yeast_root}` | `yeast_root` / `WHOLECELL_YEAST_ROOT` | the yeast model input tree (G. Magnusson): `inputs/`, `yeast_cell_dev_bundle/pipeline/data/`, `yeast_cell_dev_bundle/work/pdb/`, `runs/{X5_seed0,P2_seed0}/`. Only the yeast cell needs it. |
| `{afdb_root}` | `afdb_root` / `WHOLECELL_AFDB_ROOT` | the AlphaFold DB v6 proteome folder holding `UP000005640_9606_HUMAN_v6.tar` (https://ftp.ebi.ac.uk/pub/databases/alphafold/v6/); Jurkat only |
| `{geom_root}` | `geom_root` / `WHOLECELL_GEOM_ROOT` | immutable cellgeom_v1 geometry packages, as `<cell>/<geometry_id>/` |
| `{results_root}` | `results_root` / `WHOLECELL_RESULTS_ROOT` | builds, as `<cell>/<build_id>/`. It also holds the validation reference builds and the yeast donor build. |

## 2. Size overview (approximate)

| Root | Rows | Size | Largest items |
|---|---|---|---|
| data_root | 81 | ~21.6 GB | PDB assembly-1 mmCIF set (15.3 GB, 48,189 files; only the files the library build reads are needed: `<build>/library/work/pdb_assemblies_read.tsv`), chromatin v8 rep5 (2.6 GB), mRNA beads v1 (1.1 GB), v7 cytoskeleton (0.43 GB), Jurkat EM geometry v7 (0.56 GB) |
| yeast_root | 22 | ~18 GB distinct | `pipeline/data` (15.1 GB, of which ~15.0 GB is `validation/` literature PDFs and supplements the build does not read), AFDB yeast tar (1.0 GB), AFDB yeast dimers (0.97 GB), AF3 pool summary (0.44 GB), IntAct dump (0.84 GB, inside `inputs/`) |
| afdb_root | 1 | 5.2 GB | AlphaFold DB v6 human proteome tar |
| geom_root | 2 (+ testcell) | 1.47 GB | Jurkat `em_cell1_v8r` (1.38 GB), yeast `analytic_x5_s0_f5f3a5c723b2_memhalf1` (87 MB) |
| results_root (references) | 4 | ~54 GB as stored | yeast donor pack `20261119-bedcea9-s0-paper` (8.8 GB); the Jurkat Step-2 build (39 GB) and the yeast X5-repro build (6.1 GB) are validation references, and only their `validate/suite.json` plus a few files are needed |

## 3. Jurkat (`config/cells/jurkat.yaml`)

| Category | Inputs (manifest names) | Source | Redistribution |
|---|---|---|---|
| Proteomics | `inv_geiger2012_pg` | Geiger 2012, PXD002395 / Expression Atlas E-PROT-20 | fetch from source |
| SEC-MS | `inv_sec_fossati_features`, `inv_sec_fossati_fractions` | Fossati 2021, PXD018033 | fetch from source |
| SEC-MS (lab) | `inv_sec_lab2026_jurkat_nt`, `inv_sec_lab2026_jurkat_tcract24h` | Beltrao lab 2026 (CC BY-NC 4.0) | processed tables in HF inputs; raw to PRIDE (TODO) |
| Surface proteomics (lab) | `inv_surface_jurkat_glyco2026` (presence evidence + PM topology fix only) | Beltrao lab 2026 (CC BY-NC 4.0) | processed table in HF inputs; raw to PRIDE (TODO) |
| Transcriptomics | `inv_rna_emtab2706`, `inv_rna_gse145453_a1/a2`, `inv_rna_gse276689_rep1/2`, `inv_rna_gse301789`, `inv_rna_hpa` | Expression Atlas, GEO, HPA v25 | fetch from source (HPA: CC BY-SA 4.0) |
| Structures | `lib_afdb_human_tar` (`{afdb_root}`), `lib_afdb_human_tar_index`, `lib_afdb_complex_models`, `inv_afdb_*`, `inv_pdb_assemblies` / `lib_pdb_assemblies`, `inv_npc_7r5k`, `inv_sifts` | AlphaFold DB (CC BY 4.0), wwPDB (CC0), SIFTS | HF inputs (AFDB complexes: dated fetch; PDB set as tar shards; SIFTS snapshot); AlphaFold DB human tar: from EBI |
| Structure library | `inv_struct_library` / `lib_structures_v7` (`jurkat/yc/pipeline/lib`, 16,246 npz), `inv_tm_route_libs`, `lib_afdb_fragment_index_tmroute`, `inv_membrane_lib_ids` | this work (v7 library) | HF inputs |
| Assembly | `inv_complex_portal`, `inv_afdb_heterodimer_table`, `inv_afdb_homodimer_table`, `inv_sec_assembly_thresholds` | Complex Portal, AFDB/NVIDIA dimers | HF inputs |
| Interactions | `inv_string_links` (organisation layer; off in the published builds since v9/u2, so not read) | STRING v12.0 (CC BY 4.0) | fetch from source |
| Annotation | `inv_uniprot_full`, `inv_uniprot_reviewed`, `inv_uniprot_domains`, `inv_uniprot_topology`, `lib_tm_route_uniprot`, `lib_uniprot_human_table` | UniProt (CC BY 4.0) | HF inputs `cc_by/` (dated snapshots: the REST answers change over time) |
| Localisation | `inv_hpa_subcellular`, `inv_hpa_celllines`, `inv_hpa_jurkat_manifest`, `inv_itzhak2016` | HPA v25 (CC BY-SA 4.0; tables only, no images), Itzhak 2016 eLife (CC BY) | fetch from source / journal |
| Pre-placed structures (derived) | `pre_chromatin_v8_rep5(_json)`, `pre_cytoskeleton_v7(_json)`, `pre_lamina_v7(_json)`, `pre_chromatin_drop_lamina_v7(_npy)`, `pre_mrna_v1`, `pre_mrna_drop_v7`, `pre_mrna_nascent_v8`, `pre_npc_v3(_json)`, `pre_polysomes_v7(_json)`, `inv_*_json` summaries, `inv_mt_centerlines` | this work. Chromatin GPU MD and U-Net microtubule tracing are not bitwise reproducible, so these files are frozen. Since v9 `jurkat.yaml` no longer declares `pre_npc_v3(_json)`, `pre_polysomes_v7(_json)`, `pre_mrna_drop_v7`, `pre_cytoskeleton_v7_json` (v7 preset only) and rule L3 reads the build's own lamina (`pre_lamina_v7*` are compare references). | HF inputs |
| Order / regression references (derived) | `inv_legacy_ingredients_v7`, `inv_legacy_complexes_v0`, `inv_v4pack_loc_compare`, `inv_v5pack_loc_compare`, `inv_tier2_v5`, `inv_tier2_transfer` | this work (legacy v4–v7 models; derived from earlier packs, so fixed inputs: docs/REBUILD.md §2) | HF inputs |
| Hold-out split | `inv_fit_genes`, `inv_holdout_genes` | this work (seed 20261110) | HF inputs |
| Geometry | `geom_jgeom_v7` (source of the v7 package; declared by `jurkat_v7_port.yaml` only since v9), `geometry_jurkat` = `{geom_root}/jurkat/em_cell1_v8r` (manifest sha256 `6c86b168…`) | this work, from OpenOrganelle jrc_jurkat-1 FIB-SEM | HF inputs |
| Validation references | `val_chromatin_*`, `val_opencell`, `val_ribo_em`, `val_hpa_*`, `val_loc_*`, `val_uniprot_human`; also `ref_jurkat_v7_reference_validate`, `ref_jurkat_step2_build`, `hidden_v7_reference_tables` | OpenCell, OpenOrganelle, HPA, UniProt; this work | HF inputs (references: builds dataset) |

## 4. Yeast (`config/cells/yeast.yaml`)

| Category | Inputs | Source | Redistribution |
|---|---|---|---|
| Proteome | `inv_yeast_paxdb` | PaxDb integrated (id 4098344127) | fetch from source |
| SEC-MS | `inv_sec_swamy2022_yeast` (PXD028358 + Data S5), `inv_sec_pang2020_yeast` (PXD019513) | PRIDE | our re-exports in HF inputs (`ours/`) |
| SEC-MS (lab) | `inv_sec_lab2026_yeast_nacl`, `inv_sec_lab2026_yeast_ypd` | Beltrao lab 2026 (CC BY-NC 4.0) | processed tables in HF inputs; raw to PRIDE (TODO) |
| Localisation | `loc_yeast_cyclops_wt1..3` (CYCLoPs, Koh 2015), `loc_yeast_morgenstern2017`, `loc_yeast_vogtle2017`, `loc_yeast_uniprot_features` | CYCLoPs DB; journal supplements; UniProt | fetch from the source / journal |
| Structures | `yeast_afdb_tar`, `yeast_afdb_dimers`, `yeast_pdb_cache` | AlphaFold DB, wwPDB | dimers and PDB cache in HF inputs; AlphaFold DB yeast tar: from EBI |
| In-house AF3 pool | `yeast_af3_pool_summary` (in-house summary, 444 MB); `yeast_af3_pool_pairs` (cached per-pair PAE table, `yeast_x5_core.yaml` only) | Beltrao lab pooled AlphaFold 3 screen (https://github.com/jurgjn/pooled-ppi-yeast; public release v26.1 has fewer pairs) | the summary is not deposited (read only by `yeast_x5_port.yaml`); the pool pairs are in HF inputs `af3_pool/` (AlphaFold 3 output terms + the screen's usage request); `yeast.yaml` reads neither |
| Annotation / complexes | `inv_yeast_uniprot`, `inv_yeast_complex_portal`, `hidden_yeast_uniprot_topology`, `hidden_yeast_string_links` | UniProt, Complex Portal, STRING v12 | fetch from source |
| Yeast model data folder | `yeast_pipeline_data` (whole `pipeline/data`, manifest checksum) | yeast model | HF inputs hold the read subset (minus third-party tables with unclear terms, see `EXCLUDED.tsv`): `plugins/yeast/data/pipeline_data_read_set.txt` (782 of 2,631 files, ~1.1 GB of 14.5 GB; since v9 the engine setup check requires only these). Journal supplements in it: check terms |
| Frozen yeast model state | `yeast_x5_lib_start` (342 library entries: tRNA Boltz-2 models, builder templates), `yeast_x5_work_seed` (rcsb cache, organisation pool pairs, interaction pair summary, sweeps, validation tables) | yeast model (X5 run) | HF inputs (`ours/`; the pool pairs in `af3_pool/`, public source https://github.com/jurgjn/pooled-ppi-yeast; optional for `yeast.yaml`) |
| Size scaling | `inv_yeast_lanz2024_s3` | Lanz 2024 NSMB supplement | fetch from journal |
| Geometry | `geometry_yeast` = `{geom_root}/yeast/analytic_x5_s0_f5f3a5c723b2_memhalf1` (manifest sha256 `230b96d9…`) | this work; deterministic, from the vendored engine, but a regenerated package records its own commit and time in the manifest, so it does not match the pin | HF inputs |
| Validation references | `ref_yeast_x5repro_build` (suite no-regression reference, optional); `val_yeast_buildW_*`, `val_yeast_x5_build_*` are declared only by `yeast_x5_port.yaml` since v9 | yeast model runs; this work | HF |
| Donor pack | `donor_yeast_pack_bedcea9` = `{results_root}/yeast/20261119-bedcea9-s0-paper` | this work | HF. The paper yeast build is a downstream-only reuse of this pack (`scripts/build_downstream_reuse.py`; record in `DOWNSTREAM_REUSE.json`), so reproducing it exactly needs the donor pack. |

Not data, but required: two vendored yeast steps run with a separate interpreter, `{yeast_python}` (site key
`yeast_python`, default `{data_root}/conda_envs/yeast_x5/bin/python`; config `plugin_steps[].params.python`). Recreate that environment from
`envs/yeast_x5.pip.txt` (see `envs/README`).

## 5. Testcell

The testcell needs no download. Its inputs are in the repository (`tests/data/testcell/`, including the 45-structure
library), and its geometry package `{geom_root}/testcell/analytic_v1` is generated locally:
`python -m toolkits.analytic.generate toolkits/analytic/specs/testcell.yaml`.

## 6. Inputs not declared in `frozen_inputs` (hidden inputs)

These inputs are read by code or through config keys but are not checksum-verified by `check_frozen`:

- (v8 only) `{data_root}/data/reference/uniprot_glyco/uniprot_human_glyco_topo.tsv.gz` and
  `{data_root}/data/derived/afdb_plddt_summary.tsv.gz`, read by `wholecell/library/topology.py`: declared as frozen
  inputs `lib_uniprot_glyco_topo` and `lib_afdb_plddt_summary` since v9.
- Yeast: `{yeast_root}/inputs/` is passed whole to the vendored engine. The STRING links and `uniprot_topology.tsv` are
  read from `pipeline/data` (covered only by the `yeast_pipeline_data` manifest checksum).
- Jurkat regression and validation references: the `inventory.verify` tables and `validate.compare_reference` files
  under `{data_root}/data/derived/` (legacy v7 outputs, including `wc/tiles_v7/` at 28.6 GB), plus
  `validate.compare_reference.legacy_pipeline_dir` = `{data_root}/jurkat/yc/pipeline`. These are used for comparison
  only, never as model inputs.
- Suite references in `{results_root}`: the Jurkat Step-2 build and the yeast X5-repro build.
- The legacy library builders in `wholecell/library/legacy/*.py` read `{data_root}/data/derived/ingredients_v0.tsv` and
  `inventory_v0_complexed.tsv`, but only when the library is rebuilt from scratch. The paper builds use the frozen v7
  library, so they do not read these files (confidence medium).

The remaining items are comparison references: the build runs without them (the inventory compare reports them
missing, `validate compare_v7` is skipped with a reason, the suite marks no-regression rows n/a).

## 7. Redistribution limits (summary)

- **HPA** (CC BY-SA 4.0): the repository and the deposits contain no images, only tables. The HF inputs dataset holds the
  HPA tables the builds read and the tables derived from them under share-alike terms (`cc_by_sa/`).
- **Journal supplementary tables** (Lanz 2024, Morgenstern 2017, Vogtle 2017; Itzhak 2016 is CC BY): fetch from the
  journal.
- **CORUM** (academic use only): not read by the paper builds; if used, fetch it from the source.
- **AlphaFold 3 pool** (yeast): the in-house `summary_models.parquet` is not deposited; the per-pair PAE table read by
  `yeast_x5_core.yaml` is in `af3_pool/` under the AlphaFold 3 output terms and the usage request of
  https://github.com/jurgjn/pooled-ppi-yeast (no licence file there).
- **CYCLoPs**: database terms to be checked; fetch from the source.
- **Yeast `pipeline/data`**: contains third-party literature PDFs; only the read subset is deposited, without the
  third-party tables whose terms are unclear (`scripts/deposit/yeast_read_set_terms.tsv`).
- **Lab SEC-MS and surface glyco-capture**: no longer confidential (P. Beltrao, 2026-11-23). Processed tables go to
  the HF inputs dataset (CC BY-NC 4.0) and raw files to PRIDE (TODO). The config `licence` fields say so since v9.
- **Public databases** (AlphaFold DB CC BY 4.0, PDB CC0, UniProt CC BY 4.0, STRING CC BY 4.0, Complex Portal, SIFTS,
  GEO, Expression Atlas, PRIDE): attribute the source; dated snapshots of time-dependent exports are in the HF inputs dataset.
- **Our derived files** (pre-placed structures, libraries, geometry, references): CC BY-NC 4.0 (see `LICENSE-DATA`).

## 8. Verifying inputs

Each build runs the `check_frozen` stage (`wholecell/build.py::verify_frozen`) first. It checks every declared frozen
input against its md5 or sha256. For a directory, it checks the checksum of the manifest file named in the config;
if the config names no manifest, it checks only that the directory exists. The result is written to
`<build>/frozen_inputs.tsv`, and the build stops on any mismatch. Geometry packages are verified against
`geometry.manifest_sha256`. To compare your downloads with the paper builds, compare your `frozen_inputs.tsv` with
the one shipped in each paper build deposit.
