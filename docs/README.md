# wholecell docs

Read in this order. Each layer assumes the one before it.

1. [overview.md](overview.md): what a build is, the four parts of a cell, the stage list.
2. [getting_started.md](getting_started.md): build the synthetic testcell in about 5 minutes.
3. [stages/](stages/): one page per stage (purpose, inputs, outputs, config keys, checks, Methods text).
4. [adding_a_cell.md](adding_a_cell.md): add a new cell type, with the testcell as worked example.
5. [simulations.md](simulations.md): applications that run on a finished build.
6. [reference/](reference/): file formats and config specs; [reference/history/](reference/history/) holds frozen
   reproduction criteria, reports and the unification plan (never edited).

Also: [ARCHITECTURE.md](ARCHITECTURE.md) (module boundaries and contracts), `PORTING_STATUS.md` (what runs today),
[data_inputs.md](data_inputs.md) (input data, sources, licences, deposits; machine-readable `data_manifest.tsv`),
[release_notes.md](release_notes.md) (what the public release contains and how it was staged from the development
repository), `config/paper_builds.yaml` (pins of the two paper builds).
Renamed stages: builds made before a rename keep the old folder and file names (`prepass/` for `preplaced/`, A0c);
the one map of old names is `LEGACY_NAMES` in `wholecell/build.py`.
Config key reference: planned (generated from `config/defaults.yaml` comments once that file exists).
