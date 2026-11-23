---
license: other
license_name: mixed (per folder, see below)
pretty_name: wholecell inputs (Jurkat T cell and budding yeast structural cell models)
tags:
- biology
- structural-biology
- proteomics
- whole-cell-model
---

# wholecell inputs

Input data of the two whole-cell structural models of Magnusson, Beltrao et al. (in preparation), a Jurkat T cell and a
budding yeast cell, as read by the pipeline https://github.com/TODO/wholecell (`config/cells/jurkat.yaml`,
`config/cells/yeast.yaml`). With this dataset, the AlphaFold DB proteome tars and the few third-party tables listed
in `EXCLUDED.tsv`, the two cells can be rebuilt (`docs/REBUILD.md` in the pipeline repository).

The repository id `pbeltrao/wholecell-inputs` is a placeholder until the deposit exists.

## Get the data
Use the pipeline's fetch script, which places every file where the cell configs expect it (under the roots of your
site file) and checks md5s:

```
export WHOLECELL_SITE=/path/to/site.yaml
python scripts/fetch_inputs.py --cell jurkat        # or yeast, all; --group reference for the reference presets
```

## Layout
`<licence folder>/<root>/<path>`. `<root>` is a key of the pipeline's site file (`data_root`, `yeast_root`,
`geom_root`) and `<path>` the path below it that the configs use (`data_manifest.tsv` of the pipeline). Folders with
more than 1,000 files are shipped as tar shards `<path>.partNNN.tar` (members in `TARS.tsv`), which the fetch script
unpacks into `<root>/<path>/`.

| file | content |
|---|---|
| `MANIFEST.tsv` | every file: path, kind (file or tar shard), bytes, md5, licence, target root and path, group (required / reference / optional), cell, input name of the pipeline |
| `TARS.tsv` | the members of each tar shard, with their md5 where known |
| `EXCLUDED.tsv` | inputs that are not here, why, and where to get them |
| `fetch_inputs.py` | copy of the pipeline's fetch script |

Groups: `required` = read by `jurkat.yaml` or `yeast.yaml` (including the B5 lamina comparison files);
`reference` = read only by the reference presets (`jurkat_v7_port.yaml`, `yeast_x5_core.yaml`, the legacy Jurkat
pipeline tree that made the earlier model versions); `optional` = not needed for a build.

## Licences (per folder)

| folder | licence | content |
|---|---|---|
| `ours/` | **CC BY-NC 4.0** | data made in this work (Beltrao lab, ETH Zurich; G. Magnusson): EM-derived and analytic geometry packages, chromatin model, microtubule traces, pre-placed structures, structure library v7, rule inputs derived from earlier model versions, hold-out split, lab SEC-MS and surface glyco-capture tables (raw files at PRIDE, accession TODO), the yeast model's start library, work seed, curated tables and Boltz-2 tRNA models, re-exports of public SEC-MS data |
| `cc_by/` | CC BY 4.0, third party | UniProt exports (dated), AlphaFold DB dimer tables and models (dated fetch), SIFTS, STRING v12, Itzhak et al. 2016, OpenCell, yeast annotation tables (SGD, Alliance, IntAct-derived), CZ cryoET portal and IDR tables, PLoS One supplements. Attribute the sources (listed in the pipeline's `data_manifest.tsv`) |
| `cc0/` | CC0 1.0, third party | wwPDB assembly files, Complex Portal, PRIDE tables (Fossati et al. 2021), GEO tracks, EMPIAR-11416 meshes |
| `open/` | terms of the source (public, no restriction stated) | Expression Atlas and GEO expression tables (Jurkat RNA, Geiger et al. 2012 via E-PROT-20) |
| `cc_by_sa/` | CC BY-SA 4.0 | Human Protein Atlas v25 tables and tables derived from them (no images) |
| `af3_pool/` | AlphaFold 3 output terms + usage request | per-pair best PAE of the Beltrao lab pooled AlphaFold 3 screen of the yeast proteome (https://github.com/jurgjn/pooled-ppi-yeast); read only by the reference preset `yeast_x5_core.yaml`. AlphaFold 3 outputs are subject to the AlphaFold 3 Output Terms of Use; the screen's authors ask that proteome-wide studies are not published from these predictions before their own publication |

## Not included
See `EXCLUDED.tsv`: the AlphaFold DB v6 proteome tars (download from EBI), journal supplementary tables and database
exports whose terms do not clearly allow redistribution (fetch from the journal or database into the same path), the
in-house AF3 pool summary of the yeast X5 run, and the comparison builds (separate builds dataset).

## Provenance and checksums
Every file has an md5 in `MANIFEST.tsv`; the pipeline's first build stage checks the declared inputs against the md5s in
the cell configs (`<build>/frozen_inputs.tsv`). Origins, versions and dates are in the configs (`frozen_inputs.*.origin`,
`version`) and in `data_manifest.tsv`. Several inputs were derived from earlier versions of the models (placed shares
of the v4/v5 Jurkat packs, the v7 cytoskeleton, the legacy mRNA beads and structure library, the yeast X5 run's
start library and work seed); they are fixed inputs of the pipeline, described in `docs/REBUILD.md`.

## Citation
Magnusson, Beltrao et al., in preparation. Please also cite the third-party sources you use.
