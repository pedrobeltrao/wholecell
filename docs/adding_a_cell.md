# Adding a cell

A new cell needs four things: a config, a geometry package, data adapters (only if a source format is new) and
plugins (only for biology the core lacks). Nothing in `wholecell/` should need editing; if it does, that is a hidden
cell-specific constant and a bug in the core.

The planned one-line start, `python -m wholecell new-cell mycell --from testcell` (plan §8), **does not exist yet**.
The manual steps below work today. The worked example is the testcell (`config/cells/testcell.yaml`).

## 1. Copy the config

```
cp config/cells/testcell.yaml config/cells/mycell.yaml
```

Set `cell_id: mycell` (must equal the file stem), `name`, `organism`, `seeds.master`. Keep `plugins: []` unless step 5
applies. Config sections follow the stage order (geometry, inventory, library, preplaced, pack, validate, export,
frozen_inputs). Full key list: [reference/cell_config.md](reference/cell_config.md). There is no
`config/defaults.yaml` yet, so every key the testcell sets must be set in your file too.

## 2. Make a geometry package

For a cell made of simple shapes, use the analytic toolkit. Copy the spec and edit shapes and sizes (nm, model frame
centred on the cell):

```
cp toolkits/analytic/specs/testcell.yaml toolkits/analytic/specs/mycell.yaml
# set cell_id: mycell, geometry_id (e.g. analytic_v1), seed; edit cell / nucleus / mito / er shapes
PYTHONHASHSEED=0 <vcell python> -m toolkits.analytic.generate toolkits/analytic/specs/mycell.yaml
```

It writes `{geom_root}/mycell/<geometry_id>/`, validates it and prints the `manifest_sha256`. Put the path and the
hash into `geometry.package` and `geometry.manifest_sha256`. A package is immutable: a changed shape means a new
`geometry_id`. For a segmented EM volume, `toolkits/em` is the generator (ported, run end to end only for the Jurkat
v7 conversion so far). Spec: [reference/cellgeom_v1.md](reference/cellgeom_v1.md).

## 3. Provide the data

The testcell uses the generic table adapters (`wholecell/inventory/generic.py`):

| Adapter | Config | File (testcell) | Columns |
|---|---|---|---|
| `proteome_table` | `inventory.adapters.proteome` | `tests/data/testcell/proteome.tsv` | gene, uniprot, mass_da, copies |
| `localisation_table` | `inventory.adapters.localisation` | `tests/data/testcell/localisation.tsv` | gene, compartments, transmembrane, n_tm, signal_peptide, mito_om |
| `complex_table` | `inventory.adapters.assembly` | `tests/data/testcell/complexes.tsv` | complex_id, subunits, copies, struct_src |

Declare every file under `frozen_inputs` with `path` (use `{repo}` or `{data_root}`), `md5`, `why`, and for
publication `origin`/`version`, `licence`, `redistributable`. The `frozen` stage checks them before anything reads
them. Set `inventory.total_protein_pg` (required, no default) and `total_protein_source`. For measured proteomics
use the generic copy-number step `inventory.copies` instead of a fixed table (see Jurkat and yeast configs and
[stages/composition.md](stages/composition.md)). A new source **format** needs a new registered adapter
(`adapter.<kind>`); a new file in a known format does not.

## 4. Library and regions

- `library`: `mode: frozen` (Jurkat, testcell) or `mode: build` (yeast: a scripted library from frozen start entries,
  AlphaFold DB models and curated PDB assemblies; [stages/structures.md](stages/structures.md)). With `frozen`, point it to a folder of species npz files with a committed md5
  manifest (`python -m wholecell.library.stage make-manifest <dir> <out.tsv.gz>`). Species in the pack table without a
  library entry are skipped by the packer (`wholecell/pack/stage.py` keeps only `ing_id`s found in the library), so
  check that every species you want placed has one. The testcell copies 45 entries from the Jurkat v7 library (`tests/data/testcell/lib`).
- `pack.regions`: one entry per packing region, using label and membrane names of your geometry package. Then
  `pack.region_map` maps each localisation key of your data to a soluble and a TM region.
- `preplaced.modules: []` unless you have pre-placed structures ([stages/preplaced.md](stages/preplaced.md)).

## 5. Plugins (only if needed)

A plugin is a package `plugins/<name>/` whose import registers components with `wholecell.registry.register(kind,
name)` (adapters, inventory rules, assembly steps / hooks, pre-placed builders, virtual regions). List it in `plugins:`.
Examples: `plugins/jurkat` (rules v4-v7, virtual regions), `plugins/yeast`.

A step the core stages do not have at all is a **plugin step** (kind `plugin_step`, `fn(cfg, bdir, params)`), listed
under `plugin_steps:` with the stage it belongs to, the build files it reads and writes, and the core rules that must
wait for it (`before:`); it runs as `python -m wholecell <stage> --step plugin:<id>` and becomes the rule
`<stage>_<id>` (`wholecell/plugin_steps.py`). The yeast config is the worked example (protein table, library species,
the vendored yeast validation scripts); `plugins/yeast/steps.yaml` records for every yeast step whether it is core,
toolkit or plugin. There is no other way to run cell code in the DAG (`pipeline_engine` exists only for the X5
reproduction config `config/cells/yeast_x5_port.yaml`).

## 6. Build and read the validation

```
git add config/cells/mycell.yaml toolkits/analytic/specs/mycell.yaml   # builds refuse a dirty tree
git commit -m "mycell: first config"
scripts/build_isolated.sh mycell first 8
```

Then read `pack/summary.json` (placed vs expected mass per region) and `validate/*.json`; see
[getting_started.md](getting_started.md) §4. Steps without references write `<step>.skipped.json` with the reason.
Add references under `validate:` when you have them.
