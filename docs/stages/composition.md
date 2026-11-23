# Stage 2: composition (copy numbers)

**Purpose.** Give every protein a copy number per cell, from measured abundances scaled to a total protein mass, with
a ledger that records why each value was set.

**Today.** Part of the `inventory` stage (steps `rna`, `copies`, `v0`). The split into its own stage is part of
plan steps A4/A5.

**Inputs.** Proteome source tables and RNA tables declared in `frozen_inputs`, read by registered adapters
(`adapter.rna`, `adapter.proteome`). Jurkat: Geiger et al. 2012 iBAQ (Jurkat columns) plus RNA-based imputation for
undetected genes. Yeast: PaxDb ppm, no imputation. Testcell: a synthetic table.

**Outputs** (`inventory/`): `rna/rna.tsv` (when an RNA adapter is set); `copies/{copies.tsv, copy_ledger.tsv,
copy_sources.tsv, anchors.tsv, summary.json, report.txt}` (when `inventory.copies` is set); `v0/inventory_v0.tsv`
(per-gene copies plus first localisation weights); canonical `gene_abundance.tsv`, `proteome_sources.tsv`.

**Config keys** (`inventory:`): `total_protein_pg` (required, no default: Jurkat 100 pg, DECISION 2026-11-08; yeast
3.5 pg; testcell 0.4884 pg), `total_protein_source`, `adapters.{rna, proteome}`, `proteome_sources`, and `copies.*`:
`method_id`, `annotation`, `universe`, `value_column`, `sources`, `combine`, `imputation`, `absent`,
`per_protein_scaling`, `total_protein`, `adjustments`, `round_copies`, `anchor_check`, `anchors`, `future_sources`.
Spec: [../reference/inventory_schema.md](../reference/inventory_schema.md) §4b.

**Checks.** Anchor check (implied copies of known structures, e.g. ribosomes, nucleoporins, histones; diagnostic
only). Regression: `scripts/copies_regression_check.py` (Jurkat C1-C7, v0 and pack table md5-identical) and
`scripts/copies_yeast_check.py` (yeast Y1-Y4 vs build W), criteria in
[../reference/history/copies_regression_criteria.md](../reference/history/copies_regression_criteria.md).

**Code.** `wholecell/inventory/{copies.py, stage.py, generic.py}`, `wholecell/inventory/COPIES.md`; Jurkat RNA
adapter in `plugins/jurkat`.

**Jurkat vs yeast today.** The method is already common: yeast copies are reproduced by the same step (criterion I3,
max difference 2.6e-16); since A9 the yeast build uses them (the assembly's protein table takes copies from this step and
fails if they differ from the engine's). Declared, disabled
slots: SEC-MS abundance cross-check and surface-proteomics PM copies (planned surface adapter, Jurkat v8).

**Methods text (draft).** Protein copy numbers were derived from absolute or relative proteomics (Jurkat: iBAQ values
of Geiger et al. 2012; yeast: PaxDb). Relative values were converted to molar fractions and scaled so that the summed
protein mass equals the total protein per cell (Jurkat 100 pg; yeast 3.5 pg). Undetected Jurkat genes were imputed
from RNA expression by a log-log fit on detected genes, capped at the 5th percentile of detected values. Proteins
known to be absent were set to zero with a recorded reason, and every change was written to a ledger.
