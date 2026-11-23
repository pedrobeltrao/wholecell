# Notices: licences, provenance and third-party terms

## Licences of this repository

| What | Licence | File |
|---|---|---|
| Source code: `wholecell/`, `plugins/` (incl. the vendored yeast engine), `toolkits/`, `workflow/`, `Snakefile`, `scripts/`, `tests/*.py`, `envs/` | PolyForm Noncommercial 1.0.0 | `LICENSE` |
| Everything else we generated: configuration (`config/`), curated data tables (`plugins/*/data/`, `toolkits/*/data/`, `toolkits/*/datasets/`, `tests/data/`), docs (`docs/`, `*.md`), figures, the model outputs of the builds and the deposits listed in `data_manifest.tsv` | Creative Commons Attribution-NonCommercial 4.0 International (CC BY-NC 4.0) | `LICENSE-DATA` |

Commercial use (of the code or the data) needs a separate licence: see `COMMERCIAL.md`.
Decision: Pedro Beltrao, 2026-11-23 (supersedes the CC BY-NC 4.0-for-everything decision of 2026-11-16). G. Magnusson
agrees.

The `Required Notice` line at the end of `LICENSE` ("the wholecell authors (P. Beltrao lab, ETH Zurich; G. Magnusson)")
is a placeholder wording to be confirmed (the yeast code states no rights holder of its own).

## Provenance of copied code
- **Yeast engine** (`plugins/yeast/engine/`, plus yeast-derived modules elsewhere): written by Gudmundur Magnusson
  (yeast cell model, development repository `yeast_cell_dev_bundle`, commit e3b7448; archived public predecessor
  github.com/zoidpuff/yeast_cell_viewer). Files carry a header line `# origin: yeast yeast_cell_dev_bundle/<path> @ <commit>`.
  G. Magnusson agreed to the inclusion of this code under the licence of this repository (2026-11-23). The parenthetical
  "(no licence yet; local only, see wholecell/NOTICE.md)" in those header lines, and the same wording that
  `toolkits/analytic/yeast.py` writes into the provenance of new geometry packages, is historical: the headers are kept
  unchanged as provenance records and that wording no longer applies. Changes to the vendored files are listed in
  `plugins/yeast/engine/PORT_PATCHES.md`; original checksums in `plugins/yeast/engine/ORIGIN_sha256.json`.
- **Jurkat model code** (earlier Jurkat structural cell project of the Beltrao lab; archived public predecessor
  github.com/pedrobeltrao/jurkat-structural-cell): ported files carry `# origin: jurkat/<path> @ legacy workspace (2026-11-15)`.
- No other third-party source code is vendored (no foreign copyright headers, no minified or bundled libraries, no fonts
  or icons). The yeast viewer (`plugins/yeast/engine/web/app`) declares npm dependencies (molstar, MIT; dev tools
  typescript Apache-2.0, vite MIT, vite-plugin-singlefile MIT, @webgpu/types BSD-3-Clause) that are installed by npm,
  not included here.

## Third-party data in this repository
Small files only; each keeps its source's terms and is attributed here.
- `plugins/yeast/engine/pipeline/cell_meshes.npz`: cell and organelle meshes derived from the cryo-ET / FIB-SEM
  segmentations of EMPIAR-11416 (Dumoux et al. 2023; EMPIAR data CC0).
- `tests/data/testcell/lib/*.npz`: 45 coarse structures copied from the Jurkat v7 library, derived from AlphaFold DB
  (CC BY 4.0; Jumper et al. 2021, Varadi et al. 2024) and the PDB (CC0).
- `plugins/*/data/`: curated tables we derived (term maps, override rows, stoichiometries, anchors) that quote
  identifiers and values from UniProt (CC BY 4.0), Complex Portal (CC BY 4.0), PDB (CC0), HPA vocabulary terms (no HPA
  data or images) and the literature cited in each file.
- `toolkits/secms/data/`, `toolkits/secms/datasets/`: SEC-MS anchor sets and dataset descriptors (literature values
  with references in the files).

## Third-party data NOT in this repository (download; see `docs/data_inputs.md`, `data_manifest.tsv`)
- HPA images and tables (CC BY-SA 4.0): not included; fetch from proteinatlas.org (manifests only).
- Not redistributable: CORUM (academic use only), journal supplementary tables (fetch from the journal), the unpublished
  AlphaFold 3 pool predictions of the yeast project (their summary table `inputs/summary_models.parquet` is a
  checksum-verified input of the yeast config; the organisation layer of the paper builds runs without the AF3 pool).
- Public databases, attributed in `data_manifest.tsv`: AlphaFold DB (CC BY 4.0), PDB (CC0), UniProt (CC BY 4.0), Complex
  Portal (CC BY 4.0), STRING v12 (CC BY 4.0), PaxDb, OpenCell, GEO / ArrayExpress / Expression Atlas, PRIDE datasets.
- Lab data (Beltrao lab, ETH Zurich): SEC-MS 2026 (`lab2026_*`: Jurkat NT / TCR-activated 24 h, yeast YPD / NaCl) and
  the Jurkat surface glyco-capture time course. No longer confidential (Pedro Beltrao, 2026-11-23): processed tables are
  deposited with the inputs (CC BY-NC 4.0); raw MS files go to PRIDE (accession TODO).

## Dependencies (installed separately, not bundled)
Python and conda packages listed in `envs/` keep their own licences (mostly BSD / MIT / Apache-2.0; numba BSD-2,
gemmi MPL-2.0, openmm MIT/LGPL). Some conda builds pull GPL/LGPL system libraries (e.g. GNU readline, libgcc runtime
exceptions); they are installed by conda, not distributed here. The yeast engine environment export
(`envs/yeast_x5.pip.txt`) lists pymupdf / pymupdf4llm (AGPL-3.0): installed in that environment but not imported by
any code in this repository (can be left out when recreating it). Simulation tools with GPL licences used by the paper's
applications (CALVADOS 3: GPL-3.0; FreeDTS: GPL-2.0) are not used by this pipeline and not bundled.
