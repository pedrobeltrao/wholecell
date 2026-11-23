# inventory step `copies`: abundance data -> copies per cell (any cell)

One step, one code path (`copies.py`) for every cell. A cell's copy-number method is a config block,
`inventory.copies`; the Jurkat (human) and yeast methods are two settings of the same options. Spec:
`docs/reference/inventory_schema.md` §4b. Reproduction: `docs/reference/history/copies_regression_criteria.md`. Report:
`docs/reference/history/copy_numbers_report.md`.

## Run

- In a build: set `inventory.copies` and `inventory.adapters.proteome: copies` -> step `copies` (after `rna`, before
  `v0`) writes `inventory/copies/*`; `v0` reads `copies.tsv` (columns `v0_columns`).
- Alone: `PYTHONHASHSEED=0 python -m wholecell inventory --cell <cell> --build <dir> --step copies`.
- Yeast (stub config, inputs read in place from `{yeast_root}`): `python scripts/dev/copies_yeast_check.py --out <new dir>`
  (runs only this step and compares with build W).

## Config (annotated; every key optional unless marked)

```yaml
inventory:
  total_protein_pg: 100              # REQUIRED; "anchors" when total_protein.mode is anchors
  adapters: {proteome: copies}
  copies:
    method_id: <label>
    avogadro: 6.02214076e23          # yeast legacy: 6.02214e23
    annotation: {reader: uniprot_reviewed_genes | uniprot_ordered_locus | annotation_table, path: <frozen name|path>, carry: [orf, gene_name]}   # REQUIRED
    universe: {require: [sources, rna], order: sorted | annotation}
    value_column: ppb                # name of the relative value used for the molar fractions
    sources:                         # REQUIRED (>= 1 enabled); priority = list order
      - {id, type: relative_ibaq|ppb|ppm|lfq|absolute_copies, unit, reader: maxquant_protein_groups|paxdb|table, path,
         column, n_rep_column, label, enabled: true, aggregate: median_detected|mean_detected|geomean_detected|single,
         group_split: {rule: family_prefix, n: 3} | {rule: first} | {rule: equal_all},          # maxquant_protein_groups
         sample_match | sample_columns, id_column,                                               # maxquant / table
         strip_prefix: "4932.", map_on: orf,                                                     # paxdb / table
         value_columns | value_regex}                                                            # table
    combine: {mode: priority | geometric_mean, order: [ids], harmonise: median_ratio}
    imputation: {method: none, missing_label: not_quantified}
             #  {method: rna_loglog, rna: rna/rna.tsv, tpm_column, call_column, floor_percentile: 5, low_factor: 0.25,
             #   zero_calls: ["off", no_rna], labels: {expressed, low, zero}, flag_detected_off: "..."}
    absent:
      - {id, genes: [..] | {gene: reason}, patterns: [regex], reason, label, note, match_on: key | gene_name}
      - {id, file: <yaml>, genes_key | patterns_key | categories_key}        # categories -> one list per category
    per_protein_scaling:
      - {rule_id, kind: size_slope, match_on: orf, ratio: {numerator, denominator},
         table: {path, sheet, id_column, slope_column, replicate_prefix, min_replicates, max_se, clip, default: median}}
      - {rule_id, kind: factor_table, path, id_column, factor_column}
    total_protein: {mode: fixed, value_pg: null}          # null -> inventory.total_protein_pg
              #    {mode: anchors, anchors: [ids], combine: geometric_mean|median|mean|min|max, round_to: 0.5}
              #    {mode: source_sum, source_id: <absolute source>}
    adjustments:
      - {rule_id, kind: family_target, structure, match_on, members: {regex, genes}, family_sub: [pattern, repl],
         target: {count} | {density_per_um3, volume_um3, source}}
    round_copies: true
    v0_columns: [..]
    anchor_check: {use: final | scaled}
    anchors:
      - {id, structure, statistic: median, match_on, filter_regex, exclude,
         families: {name: {stoich, genes, regex}} | stoich_map: {gene: n} or {file, key} | stoich_from_species: {table, species_id}
         | members: {regex, genes} + family_sub + stoich_overrides,
         counts: [{label, value, source, primary}]}
    future_sources:                  # declared, not enabled (enabled: true raises NotImplementedError)
      - {id, type: secms_abundance | surface_proteomics, role, notes}
```

## The two methods in this schema

| Part | Jurkat (`config/cells/jurkat.yaml`, legacy v0) | Yeast (`config/cells/yeast.yaml`, build W) |
|---|---|---|
| annotation / universe | reviewed UniProt by primary gene; genes with a detection or an RNA row, sorted | UniProt one row per ORF (gene, ORF fallback); all, file order |
| source | Geiger 2012 iBAQ ppb, MaxQuant groups, Jurkat columns, median over detected replicates, group split by 3-letter family prefix | PaxDb integrated ppm, summed per ORF |
| imputation | RNA log-log fit, capped at the 5th percentile of detected ppb; low x 0.25; off / no RNA -> 0 | none (not in PaxDb -> 0) |
| absent | 6 genes (verified losses, 9p21 deletion) + interferon-cluster patterns | NOT_IN_THIS_CELL: genotype, alpha-specific, meiosis/sporulation, glucose-repressed, absent in G1 |
| per-protein scaling | none | G1 size scaling (34/42)^slope, Lanz 2024 slopes |
| total protein | fixed 100 pg | fixed 3.5 pg (chosen between anchor-implied values in round R) |
| adjustment | none | ribosomal proteins (+ ASC1) scaled to 14,250 per um3 x 18.44 um3 cytosol |
| rounding | integers | none |
| anchors (check) | histones / 28.19 M nucleosomes; Nups (7R5K) / 300 EM or 1,000 model pores; r-proteins / 1.5 M EM-corrected ribosomes | Nups (Kim 2018) / 93 pores; histones / 69,714 nucleosomes; r-proteins / 2 x 10^5 |

## Adding a source type or reader

A reader is `@register("copies.reader", "<name>")`: annotation readers return a frame indexed by key with `uniprot,
mass_da` (+ carried columns); abundance readers return `(value Series, n_rep Series)` indexed by key, in the source's
unit. A new `type` needs an entry in `SOURCE_TYPES` and, if it is not molar, a conversion in the combine part (`lfq`:
value / mass). SEC-MS abundance and surface proteomics are declared slots: implement the reader (SEC-MS from
`inventory/secms/profiles.tsv`; surface from the `adapter.surface` table), set the role (cross-check: report only; PM
copies: an `adjustments` kind that sets PM-protein copies), then allow `enabled: true`.

Origins: Jurkat arithmetic from `plugins/jurkat/inventory/inventory_v0.py`; yeast parts follow yeast_cell_viewer
`pipeline/inventory.py`, `ingredients.py`, `validation/anchors.py` at 34d123f (G. Magnusson; no licence yet, local only).
