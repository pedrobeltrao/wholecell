# Copy-number step (`inventory/copies`): regression criteria (frozen before the first run)

Written 2026-11-16, before the generic copy-number step (`wholecell/inventory/copies.py`) was run on any real data.
It replaces the abundance half of the Jurkat proteome adapter `jurkat.geiger_ibaq_v0` (abundance -> copies) by one
cell-agnostic, config-driven step that the yeast method must also fit.

## Set-up (Jurkat)

- `config/cells/jurkat.yaml` with the legacy settings expressed in `inventory.copies`: one source (Geiger 2012 iBAQ ppb,
  MaxQuant protein groups, Jurkat columns, median over detected replicates, protein groups split equally among genes
  sharing the first gene's 3-letter prefix), gene universe = reviewed UniProt genes with a detection or an RNA call,
  RNA imputation (log-log ppb ~ TPM fit on detected genes, capped at the 5th percentile of detected ppb, `low` x 0.25,
  `off`/no RNA -> 0), absent list = `zero_genes` + `zero_patterns` of `plugins/jurkat/data/inventory_v0.yaml`,
  total protein fixed at 100 pg (`inventory.total_protein_pg`), copies rounded to integers. No per-protein scaling, no
  structure-count adjustment. Anchor checks run (diagnostic, they write their own tables only).
- `inventory.adapters.proteome: copies` (reads `inventory/copies/copies.tsv`); localisation and everything downstream
  unchanged.
- Run: the inventory stage only (`rna` ... `compare`), fresh dev folder `{results_root}/jurkat/_dev_copies/regress-<commit>`,
  PYTHONHASHSEED=0, CPU only. **No pack, no prepass, no GPU.**

References: legacy `{data_root}/data/derived/inventory_v0.tsv` (md5 6691469933fae5f72250b3ea03c048ab), port build
`jurkat/20261115-ec0b403-s20261115-v7repro2` (`inventory/v0/inventory_v0.tsv` md5 b40f1ba260f8a2460d62160bbbe68342,
`pack_table.tsv` md5 d1a0099caad792d0ba42cbc52b52f0a5 = legacy `ingredients_v7.tsv`).

## Criteria

| id | what | criterion | why |
|---|---|---|---|
| C1 | `v0/inventory_v0.tsv` vs the port v0 (v7repro2) | md5 identical (b40f1ba2...) | same float expressions in the same order; the port's v0 is the established reproduction |
| C2 | `v0/inventory_v0.tsv` vs legacy `inventory_v0.tsv` | `copies` column identical for every gene (exact integers, same gene set and order); other columns: content identical up to the tolerance already documented for the port (v7repro2 `compare_legacy.json` notes: `compartments` order-only cells, `ppb` of imputed genes rel. diff <= 1e-14) | the documented port tolerance |
| C3 | `copies/copies.tsv` | `copies` = v0 `copies` for every gene; `source` labels identical | the step is the only producer |
| C4 | `pack_table.tsv` | md5 = d1a0099caad792d0ba42cbc52b52f0a5 (byte-identical to v7 `ingredients_v7.tsv`) | the full inventory chain downstream of v0 is unchanged |
| C5 | `compare_legacy.json` | every table / ingredient entry has the same identical / content-identical status as in v7repro2; B3 all pass | as for the SEC-MS regression |
| C6 | legacy adapter equivalence (test) | in `tests/test_copies.py`, on the testcell's synthetic MaxQuant/RNA tables, the generic step with the Jurkat recipe gives the same table (md5 of the TSV) as the ported legacy adapter `jurkat.geiger_ibaq_v0`, which stays registered for this comparison | guards the recipe without the 1 GB Jurkat inputs |
| C7 | provenance | every gene of `copies.tsv` has `set_by`, and every rule that changed a value has a `copy_ledger.tsv` row with its `rule_id` | requirement of the step |

Not criteria (expected to change): canonical tables that print the proteome adapter name in provenance strings
(`species.tsv` `abundance_source`, `field_provenance.tsv` `source_id` for monomers: `v0:jurkat.geiger_ibaq_v0[...]` ->
`v0:copies[...]`); new tables under `inventory/copies/`.

Verdict: REPRODUCED if C1-C7 pass. If C1 fails but C2-C5 pass, the verdict is REPRODUCED (copies) with the byte
differences listed. Otherwise the failing rows are listed and Jurkat keeps `jurkat.geiger_ibaq_v0` until fixed.

## Yeast (feasibility, read-only)

Target: per-protein `copies_per_cell` of yeast build W (`runs/P2_seed0/work/features/protein_features.parquet`, written by
`features.py` from `inventory.build()`; commit 8f84387 = the inventory code of main 34d123f), recomputed by the generic step
from the yeast inputs read in place (`pipeline/data/paxdb_integrated.txt`, `uniprot_yeast.tsv`, Lanz 2024 Table S3) with
`config/cells/yeast.yaml` (`inventory.copies`). Nothing is written, run or locked in the yeast folder (no import of the
yeast `config.py`, which creates folders on import).

| id | criterion |
|---|---|
| Y1 | same protein set (ORF) as build W |
| Y2 | max relative difference of `copies_per_cell` <= 1e-9 over all proteins (float, the yeast method does not round); proteins with 0 copies in W have 0 |
| Y3 | step summaries equal build W's log: NOT_IN_THIS_CELL categories (counts and ppm shares), ribosome factor 1.185 and target 262,722, protein total 3.61 pg after the ribosome adjustment |
| Y4 | anchor check reproduces build W `work/validation/anchors.json`: nucleoporin median ratio and core histone ratios (rel. diff <= 1e-6) when given W's pore and nucleosome counts |
</content>
</invoke>
