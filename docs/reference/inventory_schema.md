# Inventory schema (wholecell, phase 1)

What a cell's protein inventory is, how it is built, and which tables carry it. Code: `wholecell/inventory/`
(`schema.py` holds the column lists below, `stage.py` the steps, `canonical.py` the export), cell-specific parts in
plugins (`plugins/jurkat/inventory/`, curated lists in `plugins/jurkat/data/*.yaml`).

## 1. Why this shape (audits of 2026-11-15)

| Finding | Human (Jurkat v7) | Yeast | Decision |
|---|---|---|---|
| Levels | per gene `inventory_v0.tsv` -> assembly tables -> species `ingredients_v7.tsv` | per gene `inventory.build()` (in memory) -> `ingredients.csv` (+ `oligo_map.json`) | keep three levels, all written to disk, one folder per step |
| Species columns | `ing_id, kind, genes, stoich, copies, mass_da, mass_pg, comp_weights, transmembrane, n_tm, struct_src, uniprot, signal_frac, mito_om_frac` | `ing_id, label, kind, genes, copies_per_cell, comp_weights, mass_da, pdb, membrane, oligomer_n, n_tm` | one species table, human names (`copies`, `struct_src`); yeast maps 1:1 (`copies_per_cell`->`copies`, `pdb`->`struct_src`, oligomers become `kind=oligomer` rows instead of `oligo_map`) |
| Localisation | JSON `comp_weights` per species, rules rewrite it in place | JSON `comp_weights` (CYCLoPs keys) | pack view keeps the JSON; the canonical form is the long `species_loc.tsv` (one row per species x key) so a rule writes rows with its `rule_id` |
| Provenance | per-rule change tables `inventory_v{4..7}_changes.tsv` (best record, merge_notes) but 3 schemas (`reason` / `rule` / `rule_id`) and in-place rewrites | `loc_source` tags only | one `changes.tsv` schema for every rule step + `field_provenance.tsv`; no in-place rewrites |
| Ledger bug | `build_ledger._rules_table` looked for `rule`/`reason`, fell back to the last column: v6/v7 `rule_ids` = `mass_moved_pg` | - | `rule_id` is a required column of `changes.tsv`; `version_changes.tsv` reads it |
| Total protein | CLI argument (`inventory_v0.py 100`), stale default 80 pg | `TOTAL_PROTEIN_PG` 3.5 pg constant | `inventory.total_protein_pg` REQUIRED in config, no default; its source in `total_protein_source` |
| Pack dependence | v5 rule P, v6 R1 read earlier packs; v5 T reads the TM-route library list | none | such inputs are declared frozen tables (`frozen_inputs: inv_*` with md5) |
| SEC-MS, surface | Fossati 2021 SEC-SWATH used inside scripts; no surface data | none | explicit per-protein / per-complex SEC-MS tables and a surface-proteomics table (with condition / replicate), empty when no data; SEC-MS made by one toolkit for every cell (§4a, 2026-11-16) |

## 2. Species identifiers

`species_id` = the pack/library id (file `<species_id>.npz`), kept from the legacy grammar so libraries and packs stay
addressable: gene symbol (monomer), `CPX-<n>` (Complex Portal complex), `<GENE>_x<n>` (PDB homo-oligomer), `<GENE>_afx2`
(AFDB homodimer), `<A>-<B>` (AFDB heterodimer), curated names (`NUCLEOSOME`, `NPC`, `RIBO80S`). Yeast ids
(`CPX_554_<name>`, `<id>_oligo`, `<A>_<B>_het`) are valid species ids too.

`species_key` = stable, collision-free key (yeast merge note): `prot:<UniProt>`, `cpx:<CPX id>`, `homo:<UniProt>:x<n>`,
`afhomo:<UniProt>:x2`, `het:<UniProt A>+<UniProt B>`, `curated:<id>` (gene symbol when no accession).
`display_name` is free text. `kind` in {`protein`, `complex`, `oligomer`, `afdb_homodimer`, `afdb_heterodimer`}.

## 3. Pack view: `inventory/pack_table.tsv` (what the packer reads)

Columns exactly those of the legacy `ingredients_v7.tsv`, in this order (`schema.PACK_TABLE`):
`ing_id, kind, genes, stoich, copies, mass_da, comp_weights, transmembrane, n_tm, struct_src, uniprot, signal_frac,
mito_om_frac, mass_pg`. `genes`/`stoich` `;`-aligned; `copies` per cell (an assembly counts once); `mass_da` per copy,
`mass_pg` = copies x mass_da / 6.022e23 x 1e12; `comp_weights` JSON `{loc_key: weight}` summing to 1 (keys: config
vocabulary, mapped to pack regions by `pack.region_map`); `signal_frac` / `mito_om_frac` = subunit-mass fraction with a
signal peptide / outer-membrane annotation (packer: lumen only with a signal peptide; mito TM on OM vs IM).

## 4. Canonical tables (inventory/, written by step `canonical`)

All TSV with a header; `rule_id` values are `<version>.<rule>` (e.g. `v6.R1_nc`, `v4.f_kdel`) or a generic rule name.

| Table | Columns | Notes |
|---|---|---|
| `species.tsv` | `species_id, species_key, display_name, kind, genes, uniprots, stoich, n_chains, copies, mass_da, mass_pg, struct_src, tm, n_tm, signal_frac, om_frac, assembly_source, abundance_source, loc_sources, last_rule_id, build_step` | one row per species of the final step; `loc_sources` = base localisation source + every rule that set a key |
| `species_loc.tsv` | `species_id, loc_key, weight, rule_id, evidence, step` | long form of `comp_weights`; `rule_id` = last rule that changed that key's weight (from the change tables), else the base source (`v0:<adapter>[loc_source]` or `ingredients:<builder>[subunit mass-weighted]`) |
| `field_provenance.tsv` | `species_id, field, value, source_id, rule_id, step` | fields `copies, comp_weights, transmembrane, n_tm, mass_da, genes, stoich, signal_frac, mito_om_frac, struct_src`: last rule that set it, else the producing adapter/builder |
| `changes.tsv` (per rule step, `rules/<set>/`) and `changes_all.tsv` | `step, rule_set, rule_id, species_id, genes, field, old, new, evidence, check, mass_moved_pg` | one row per changed field; `field=row` (species deleted), `prepass/geometry/packer` (non-inventory history rows of v7) allowed. Jurkat rule steps also write `changes_legacy.tsv` in the legacy per-version format (byte-compatible) |
| `sources.tsv` | `source_id, kind, path, url, version, date, md5, licence, redistributable, used_by, description` | registry of every `inv_*` frozen input (+ config-only sources); `used_by` = steps that read it |
| `proteome_sources.tsv` | `source_id, weight, unit, normalisation, total_protein_pg, total_protein_source, condition, replicates, notes` | several sources with weights; the total-protein scale and its source on every row |
| `localisation_sources.tsv` | `source_id, weight, evidence_codes, term_map, notes` | weights, evidence codes (e.g. HPA reliability), term map file |
| `gene_abundance.tsv` | `gene, uniprot, source_id, value, unit, n_rep, used` | long, one row per gene x abundance source (config `inventory.canonical.abundance_columns`) |
| `secms_protein.tsv` | `uniprot, gene, source_id, condition, replicate, n_fractions, apex_fraction, apex_app_kda, mono_kda, assembled_frac, free_frac, threshold_k, calibration_id, notes` + §4a extras | per protein x condition x replicate x threshold (k x monomer MW); written by the `secms` step (toolkits/secms); Jurkat: Fossati 2021, k = 2/3/4/6 |
| `secms_complex.tsv` | `complex_id, species_id, source_id, condition, replicate, n_subunits_detected, n_subunits, coelution_peak, coelution_height, mean_corr, assembled_frac, free_frac, app_kda, expected_kda, notes` + §4a extras | per complex (Complex Portal, hu.MAP) co-elution and complex-level free fraction |
| `secms_complex_subunit.tsv` | §4a | per complex x subunit: in-complex / free fraction (the measured free pools) |
| `secms_pairs.tsv` | §4a | pairwise co-elution of candidate pairs (heterodimer evidence) |
| `secms_calibration.tsv` | §4a | calibration points (used + cross-check methods), fits, residuals, leave-one-out residuals |
| `surface.tsv` | `uniprot, gene, source_id, condition, replicate, pm_abundance, unit, pm_fraction, topology, n_term_side, evidence, notes` | surface proteomics: PM abundance, orientation (`topology` / `n_term_side` in {extracellular, cytoplasmic}), evidence; empty until an `adapter.surface` is configured |
| `version_changes.tsv` | `gene, transition, change, max_shift_target, tvd_target, shifts_target, copies_from, copies_to, reason_class, rule_ids, evidence` | gene-level target-share shifts between consecutive steps (coarse compartments via `inventory.canonical.key_compartment`), `rule_ids` from the `rule_id` column (fixed legacy ledger) |
| `secms_bound_fraction.tsv` (+ `secms_bound_fraction_summary.json`) | §4c | one row per protein (v0 gene): SEC-MS measured assembled fraction vs the modelled assembly state; written last by `canonical`, read-only on every other table (added 2026-11-19) |

### 4a. SEC-MS tables (toolkits/secms; step `secms`; added 2026-11-16)

One toolkit, one parameter set (`toolkits/secms/presets/default.yaml`) for every cell; species switch `human | yeast`
(Complex Portal taxon, native-mass anchor set). The `secms` step writes `inventory/secms/<table>.tsv` (+ `profiles.tsv`,
`secms_peaks.tsv`, `secms_calibrants.tsv`, `diagnostics.json`); `canonical` copies the first five to `inventory/`.
`replicate` is the replicate label of a run or `mean` (rows averaged over the replicates of a condition; SDs there).

| Table | Columns (in order; schema.py) | Meaning |
|---|---|---|
| `secms_protein` extras | `assembled_frac_lo, assembled_frac_hi, assembled_frac_sd, n_rep, n_peptides, void_frac, n_peaks, assembled_frac_peaks, tm, species, measurable, assembled_frac_mid` | lo/hi: apparent-MW scale shifted by -/+ 1 calibration residual SD; sd: over replicates (`mean` rows) or seeded bootstrap of a per-value standard error (matrix inputs); void_frac: intensity at fractions <= void (counted as assembled); assembled_frac_peaks: peak-level estimator (k = primary k); measurable: False = the column cannot resolve k x monomer (every a_k column NaN, "unmeasured"); assembled_frac_mid: share at >= sqrt(k (k - 1)) x monomer (midpoint between a (k-1)-mer and a k-mer; oligomer shares of the assembly rules with share source `secms_midpoint`) |
| `secms_complex` extras | `complex_source, name, stoich_known, measured, window_lo_kda, tm_any, calibrant, calibrant_rank, calibrant_kept, calibration_id, species` | measured = coverage >= min_subunits_detected and >= min_detected_frac; assembled_frac = median subunit in-complex fraction; calibrant columns = Complex Portal peak-level calibrant (diagnostic) |
| `secms_complex_subunit` | `complex_id, complex_source, uniprot, gene, stoich, source_id, condition, replicate, detected, mono_kda, in_complex_frac, free_frac, window_lo_kda, apex_fraction, apex_app_kda` | in-complex = intensity at app. MW >= max(window_lo x expected complex MW, k_floor x monomer), capped at the void's MW |
| `secms_pairs` | `uniprot_a, uniprot_b, source_id, condition, replicate, r, r_sd, n_rep, apex_delta, pair_source` | Pearson r of smoothed profiles (candidate orientation kept); `all_pairs` rows only with `pairs.all_pairs_min_r` |
| `secms_peaks` | `uniprot, source_id, condition, replicate, peak, apex_fraction, left, right, area_frac, app_kda, ratio_to_mono, cls` | cls in {monomer, assembled, void} |
| `secms_calibration` | `source_id, condition, replicate, calibration_id, method, role, kind, id, gene, fraction, kda, log10_fit, residual, kept, loo_residual` | role used / check; methods anchors (default), complex_portal, standards |
| `secms_calibrants` | `cpx, name, n, mw_kda, peak, co_height, mcorr, key, kept, source_id, condition, replicate` | Complex Portal peak-level calibrants (legacy `sec_calibrants_peak.tsv` columns first) |
| `profiles` | `uniprot, source_id, condition, replicate, kind, f1..fN` | kind smooth (used) / raw (top-N sums) |

Calibration and measurable range (2026-11-23, rebuild v9; `presets/default.yaml`): the Theil-Sen intercept is the median
anchor residual (`calibration.intercept: median_residual`; before: `medians` = median(y) - slope x median(x), which ran
1.1-1.8x above the median anchor in the v8 / u1 fits). a_k is computed only when the monomer mass is >= 10^(1 calibration
SD) x the apparent mass of the last fraction and k x monomer <= the apparent mass of the first non-void fraction /
10^(1 SD) (`assembled.range_margin_sd: 1.0`); otherwise `measurable` = False and a_k = NaN (before: a protein with
k x monomer below the last fraction scored 1 by construction). Presets `legacy_v7` and `v8_u1` keep the old rules
(`v8_u1` reproduces `inventory/secms/*` of the Jurkat v8 and yeast u1 builds). Analysis: `docs/methods_secms_bound_fraction.md` §5.

Assembly-state adapters read these tables through `wholecell.inventory.secms.SecView`; the rules that use them are
`assembly_options` (same module): `complex_copies: stoichiometry_cap | measured_free_pool`, homo-oligomer share,
AFDB homodimer share, heterodimer thresholds (defaults = Jurkat v7). Jurkat: toolkit preset `legacy_v7` reproduces the
v7 assembly tables exactly (`docs/reference/history/secms_regression_criteria.md`); `measured_free_pool` is OFF.

### 4b. Copy numbers (step `copies`; `wholecell/inventory/copies.py`; added 2026-11-16)

One cell-agnostic step turns abundance data into copies per cell for every cell (human Jurkat, yeast, others); everything
cell-specific is in `inventory.copies` of the cell config (step README: `wholecell/inventory/COPIES.md`). It runs after
`rna` and before `v0`; the proteome adapter `copies` hands its table to `v0` (`inventory.adapters.proteome: copies`).
Order of operations (each rule writes ledger rows with its `rule_id`):

| # | part | config key | options | rule ids |
|---|---|---|---|---|
| 1 | annotation (protein universe candidates) | `annotation` | reader `uniprot_reviewed_genes` (key = primary gene), `uniprot_ordered_locus` (key = gene, ORF fallback; carries `orf`, `gene_name`), `annotation_table` | - |
| 2 | typed abundance sources | `sources[]` | `type`: `relative_ibaq`, `ppb`, `ppm` (molar, as is), `lfq` (value / mass), `absolute_copies`; `reader`: `maxquant_protein_groups` (sample columns, `group_split`: `family_prefix` n / `first` / `equal_all`), `paxdb` (taxon prefix, ids summed), `table` (value columns = replicates); `aggregate`: `median_detected`, `mean_detected`, `geomean_detected`, `single`; `column`, `n_rep_column`, `label`, `enabled` | `source.<id>` |
| 3 | universe | `universe` | `require: [sources, rna]` (any of), `order: sorted \| annotation` | - |
| 4 | combination and priority | `combine` | `mode: priority` (first source with a value, list order or `order`) or `geometric_mean`; secondary sources harmonised to the first by `median_ratio` on shared proteins | `source.<id>`, `combine.geometric_mean` |
| 5 | imputation | `imputation` | `none` (`missing_label`, value 0) or `rna_loglog` (`rna`, `tpm_column`, `call_column`, `floor_percentile`, `low_factor`, `zero_calls`, `labels`, `flag_detected_off`) | `impute.rna_loglog`, `impute.rna_loglog_low`, `impute.rna_off`, `missing.zero` |
| 6 | absent lists with reasons | `absent[]` | `genes` (list or {gene: reason}), `patterns`, `categories` / `file` + `genes_key` / `patterns_key` / `categories_key`, `match_on`, `label`, `reason`, `note` | `absent.<id>[.<category>]` |
| 7 | per-protein scaling (relative values) | `per_protein_scaling[]` | `size_slope` (value x ratio^slope; table with replicate SE filter, clip, default median), `factor_table` | `scale.<rule_id>` |
| 8 | absolute scale | `total_protein` | `mode: fixed` (`value_pg`, default `inventory.total_protein_pg`), `anchors` (`anchors`, `combine`: geometric_mean / median / mean / min / max, `round_to`, `provisional_pg`; then `inventory.total_protein_pg: anchors`), `source_sum` (an absolute source keeps its mass); `avogadro` | `total.<mode>` (key `*`) |
| 9 | structure-count adjustments (after scaling) | `adjustments[]` | `family_target`: members (regex / genes), `family_sub`, target `count` or `density_per_um3` x `volume_um3`; the other proteins keep their copies | `adjust.<rule_id>` |
| 10 | rounding | `round_copies` | true (integers) / false | - |
| 11 | anchor check (diagnostic) | `anchors[]`, `anchor_check.use` | families: `families` {name: {stoich, genes / regex}}, `stoich_map` (inline or file), `stoich_from_species` (a pack-table row), `members` + `family_sub` + `stoich_overrides`; `filter_regex`, `exclude`, `match_on`; `counts[]` {label, value, source, primary}; `statistic: median` | - |
| 12 | future sources | `future_sources[]` | declared slots `secms_abundance` (cross-check), `surface_proteomics` (PM copies; data path = adapter.surface §4c, activated variant via `per_protein_scaling: factor_table`); `enabled: true` raises until implemented | - |
| - | v0 hand-over | `v0_columns` | columns of `copies.tsv` given to the v0 step (Jurkat: the legacy v0 columns) | - |

Outputs (`inventory/copies/`):

| Table | Columns | Notes |
|---|---|---|
| `copies.tsv` | key (`gene`), `uniprot, mass_da` [, carried annotation], one column per source (+ n_rep), imputation inputs (`tpm, rna_call`), `source` (label), value column (e.g. `ppb`, `ppm`), `note`, [`size_slope`, `scale_factor`], [`adjust_factor`], `copies`, `set_by` (source id or rule id that set the value), `rule_ids` (every rule that touched it, in order) | one row per protein of the universe |
| `copy_ledger.tsv` | `key, rule_id, field, old, new, factor, evidence` | one row per protein and rule; key `*` = global rows (total protein, N, mean mass) |
| `copy_sources.tsv` | `source_id, type, unit, reader, path, role, status, column, label, notes` | enabled / disabled sources and future slots (`declared, not enabled`) |
| `anchors.tsv` | `anchor_id, structure, count_label, count, family, members, stoich, copies, need, ratio` | ratio = family copies / (stoich x count) |
| `summary.json` | sources, imputation fit, absent lists (n, share of the value total), scaling rules, total protein (mode, pg, molecules, mean mass), adjustments, anchors (per anchor x count: median ratio, p10/p90, implied total = total / ratio), `anchor_summary` (implied totals of the primary counts, spread max/min, combined values), labels, protein pg | |

Implied total protein of an anchor = total protein used / median family ratio (the yeast convention: with a later
adjustment the ratio is taken on the final copies and divided into the scaling total, not the adjusted one).
Configs: Jurkat = legacy v0 arithmetic (`docs/reference/history/copies_regression_criteria.md` C1-C7), yeast = build W method (Y1-Y4).

Per-step intermediate tables keep their legacy formats: `rna/rna.tsv` (per gene TPM / zTPM / call), `v0/inventory_v0.tsv`
(`gene, uniprot, mass_da, <abundance columns>, source, copies, compartments ("key:w;key:w"), loc_source, transmembrane,
mass_pg`), `assembly/*` (adapter outputs; Jurkat: `complexes_v0, inventory_v0_complexed, sec_assembly, sec_calibrants_peak,
homo_oligomers_v0, afdb_homodimers_v0, afdb_heterodimers_v0, complex_copies_rule, assembly_options.json`; the two SEC
tables are views of `secms/*`), `ingredients/ingredients.tsv`
and `rules/<set>/ingredients.tsv` (pack-view columns).

### 4c. Per-protein SEC-MS bound fraction (`secms_bound_fraction.tsv`; `wholecell/inventory/secms_bound.py`; added 2026-11-19)

**Why.** SEC-MS measures much more assembly than the structural model represents: a protein can elute far above its
monomer mass although no structure of the assembly exists (yeast: EF-1A, EF3, Hsp70s, PAB1, clathrin; ~34 % of the
protein mass is kept monomeric with SEC a_2 ~ 0.75; `docs/unification_criteria/step3_yeast_secms_assembly_addendum_2026-11-19.md`).
DECISION (Pedro 2026-11-19): structure-less assemblies are not invented and the packing is unchanged; the gap is carried
as a per-protein attribute instead. Methods text and per-cell numbers: `docs/methods_secms_bound_fraction.md`.

**Definition.** For protein p with copies N_p (v0) and N_p^free copies in monomer species of the pack table:

- measured `secms_bound_fraction` = a_k(p), the share of the protein's SEC-MS elution profile (summed intensity) at
  apparent mass >= k x monomer mass (the toolkit's `assembled_frac` at `threshold_k`, void fractions counted as
  assembled); dimensionless, 0-1, a copy (intensity) fraction of the whole protein; NA ("unmeasured") when no dataset
  of the run can resolve k x monomer (toolkit `measurable`, §4a);
- modelled `modelled_assembled_fraction` = m_p = 1 - N_p^free / N_p (clipped to 0-1); monomer species = `kind protein`
  rows that are not library-level oligomers (`assembly/oligo_map.json` keys, yeast X5); copies in complexes, homo-oligomers,
  AF homo/heterodimers, family / pre-pass species and copies not placed all count as assembled (as in the addendum's m_p);
- `unrepresented_bound_fraction` = max(0, a_k - m_p): share of all copies that SEC-MS measures bound but the model keeps
  free; `unrepresented_share_of_free` = that gap / (1 - m_p) (0 when m_p = 1): the probability that a modelled monomer
  molecule of p is in an unmodelled assembly (the per-molecule quantity a simulation uses).

**Selection.** The SEC-MS run is the one the assembly rules read (`assembly/assembly_steps.json` options `secms`; for
other assembly adapters their `assembly_options`; else the `secms_toolkit.combine` label; else the only run): Jurkat
`jurkat_final` (lab 2026 NT + TCR-activated 24 h + Fossati 2021, combined), yeast `yeast_final` (lab 2026 NaCl + Swamy
2022 primary, lab YPD + Pang 2020 fill, combined), the v7 port Fossati 2021 only; no SEC-MS step (yeast_x5_core) = NA. k = the k of the rules (`k: rules`, default): for a protein
modelled in a homo-oligomer of state n, `k_by_state` of the homo-oligomer share (n <= 3 -> 2, n <= 5 -> 3, else 4), else
the rules' dimer k (`secms_fractions.threshold_k` / AF homodimer share k = 2). Optional config
`inventory.canonical.secms_bound_fraction: {source_id, condition, replicate, k: rules | <int>}`.

| Column | Unit / values | Meaning |
|---|---|---|
| `gene, uniprot, copies, mass_da, mass_pg` | -, -, copies per cell, Da, pg | v0 row (`mass_pg` = copies x mass) |
| `secms_source_id, secms_condition, secms_replicate` | - | SEC-MS run read (combined runs: condition = replicate = `combined`) |
| `datasets, n_datasets, confidence` | `source\|condition[;...]`, count, toolkit label | datasets behind the value (combined rows: the contributing conditions; one run: `source\|condition\|replicate`), combine confidence (`high`, `single`, `discordant`, `low`, `unmeasured`) |
| `threshold_k, k_basis` | integer x monomer MW, text | k used and why (`rules: default k 2`, `rules: homo-oligomer state n`, config) |
| `secms_bound_fraction` | 0-1 | a_k (above) |
| `secms_bound_fraction_sd` | 0-1 | SD over datasets (combined) or replicates (`mean` rows); NA for one replicate |
| `secms_bound_fraction_lo, _hi` | 0-1 | calibration band: a_k with the apparent-MW scale shifted by -/+ one calibration residual SD |
| `apex_app_kda, mono_kda` | kDa | apparent mass of the profile apex, monomer mass (an effective-size input) |
| `model_state` | `monomer_only`, `in_assembly`, `not_in_pack`, `absent` | named in no assembly species / in >= 1 assembly species with copies / named by no species (copies held by a family or pre-pass species, e.g. Jurkat histones in the curated nucleosome, or not placed) / no copies |
| `structure_represented` | bool (NA when absent) | assembly built as complex / oligomer / heterodimer in the model (`model_state` in_assembly or not_in_pack) vs kept monomeric for lack of a structure |
| `model_species, model_homo_state` | ids (top 5 by copies), n | assembly species that name the protein; largest homo-oligomer state |
| `modelled_free_copies, modelled_assembled_fraction` | copies, 0-1 | N_p^free, m_p |
| `unrepresented_bound_fraction, unrepresented_share_of_free` | 0-1 | gap (above); NA when not measured |
| `na_reason` | text | why the measured part is NA: `no SEC-MS step in this build`, `not quantified in SEC-MS run <id>`, `no UniProt accession`, `unmeasured at k ...` (k x monomer outside the calibrated range of every dataset), `a_k NaN`, `SEC-MS run not selected`, `no copies in this cell` |

`secms_bound_fraction_summary.json`: selection, k rules, total / SEC-covered mass, mass shares (structure represented,
monomer only, modelled assembled, measured bound of the covered mass, measured bound but unrepresented) and the top 20
contributors (gap x mass). Export: `wholecell/export/stage.py` `sim_manifest` copies the per-protein values onto the
monomer molecules of a simulation box (`simulation_workflows.md` §3); nothing else reads the table.

**Caveats.**
- Lysis, dilution and the column: SEC-MS runs on a diluted lysate (Jurkat Fossati 2021: 0.5 % NP-40). Weak or
  concentration-dependent assemblies dissociate (under-estimate of a_k: HSPD1 heptamer at 147-192 kDa in the lab 2026 runs,
  480 kDa in Fossati), while post-lysis aggregation, binding to released nucleic acids / membranes, detergent micelles
  (TM proteins) and co-elution of unrelated assemblies of similar size raise a_k (over-estimate). a_k is an apparent-mass
  share, not a proof of a specific partner.
- k threshold: a_k falls with k. With k = 2 the cut sits on the dimer's own mass, so a pure dimer whose apex is at ~2 x
  reads a_2 ~ 0.5 (on the median-residual calibration the dimeric anchors read 0.54 Jurkat / 0.45 yeast on average,
  the monomeric anchors 0.15 / 0.14; `docs/methods_secms_bound_fraction.md` §5). Monomers with an elongated shape or a large
  hydrodynamic radius can read a_2 > 0. The assembly rules therefore take oligomer shares from the midpoint share
  (`assembled_frac_mid`, share source `secms_midpoint`), while the bound fraction stays a_k. Compare values only at the
  stated `threshold_k` and estimator. Before 2026-11-23 the calibration intercept inflated apparent masses by 1.1-1.8x
  (dataset-dependent), which raised a_k of every protein (mass-weighted a_k of the covered mass 0.785 -> 0.666 Jurkat,
  0.809 -> 0.673 yeast with the fix).
- Filter vs no-filter datasets: lysate handling shifts a_k. A clarification spin / filter removes pelleted large
  assemblies (ribosomes, cytoskeleton, membranes; Swamy 2022: 15,000 g + 0.45 um filter, no concentration), a 100 kDa
  MWCO concentration step enriches complexes and biases a_k upwards (Pang 2020, low salt, one profile); the lab 2026 runs
  (Jurkat, yeast) have undocumented lysis / spin and void fraction. Void fractions count as assembled (Fossati void at
  fraction 8). The combined value mixes these (see `datasets`); yeast_final uses Pang and lab YPD only as fill.
- Calibration lower range: proteins whose monomer would elute within one calibration SD of the last fraction are
  unmeasured (NA), not 1 (lab yeast NaCl / YPD: monomers below ~43 kDa; Jurkat lab NT below ~21 kDa).
- Calibration upper range: apparent masses above the largest calibrant are extrapolated (yeast: top anchor 1,236 kDa,
  Pang only to 268 kDa); a_k of large k and apexes above ~1 MDa carry the calibration band `_lo`/`_hi`, not the SD.
- Condition dependence: values are for the measured conditions, combined with equal weight (Jurkat: lab 2026 resting
  NT + 24 h TCR activation + Fossati steady state; yeast: lab 2026 NaCl stress + Swamy 2022 SC log phase primary, lab YPD
  and Pang 2020 YEPD mid-log fill); `confidence = discordant` marks proteins whose datasets disagree by > 0.3
  (`discordant_range`). Activation / stress-specific assemblies enter the combined a_k of a resting / G1 model cell.
- Model side: m_p counts surplus copies not placed and family-held copies as assembled; a heterodimer with a much
  smaller partner is modelled assembled although SEC at k 2 could not see it (gap clipped at 0, never negative).
  `structure_represented = True` for `not_in_pack` assumes the copies sit in a family / pre-pass assembly (true for the
  Jurkat histones and yeast HHT/HHF; check `model_species` / `changes_all.tsv` for other cells).
- Pre-placed structures are not seen: copies consumed after the inventory by the pre-placed stage (Jurkat microtubules,
  spectrin / lamin / vimentin networks, H1 on chromatin) count as free in m_p, so their gap is over-stated (Jurkat: 2.5 of
  35.4 points; `docs/methods_secms_bound_fraction.md` §3). Being named in an assembly does not mean most copies are in it
  (`structure_represented` is a per-protein flag; use `modelled_assembled_fraction` for amounts).

### 4d. Surface proteomics (adapter.surface `surface_proteomics`; `wholecell/inventory/surface.py`; added 2026-11-19, OFF)

Generic for any cell. Not used by any build config: it runs only when a config sets
`inventory.adapters.surface: surface_proteomics` (Jurkat fragment: `config/fragments/jurkat_surface_proteomics.yaml`,
merged with `wholecell.config.load(cell, overrides=...)`). `copies.py` is unchanged; the copy-number slot
`future_sources: surface_proteomics` stays declared-not-enabled, and the supported uses go through existing hooks:

| use | how | effect on the build |
|---|---|---|
| presence / topology evidence (`mode: presence`, default) | canonical step writes `surface.tsv` (one row per protein x baseline condition: re-centred mean log2 intensity, `detected n/N`, topology class and `n_term_side` from UniProt with TM-alternation inference) | none (evidence table) |
| PM share (`mode: pm_share`) | `pm_fraction` from an offline fitted table (`pm_share_table`, header `# validated: true` required, else ValueError unless `allow_unvalidated`) | none until a localisation rule reads `surface.tsv` |
| activated variant (`mode: timepoint`, `factor_table()`) | robust per-protein surface changes of `timepoint.condition` vs `timepoint.reference` (detected >= 3 in both, \|log2\| >= 0.585, >= 3 paired replicates agree, optional q and PM-gene filter) -> table `gene, factor` for `copies.per_protein_scaling: [{kind: factor_table, ...}]` | scales those proteins' relative abundances in a variant build |

Reading: long table, one row per (protein, run[, comparison]); repeats per comparison collapse (identical measured value
required); rows with a condition outside `conditions` are dropped (counted). Re-centring (`recentre: complete_median`):
each run is shifted by the median, over proteins measured in every run, of its deviation from the protein's mean.
Intensities are relative within a protein; they are not copy numbers (peptide number / detectability differ between
proteins). Tests: `tests/test_surface.py` (synthetic tables only).

## 5. Components (registered by name, `wholecell.registry`)

| Kind | Call | Returns | Generic | Jurkat |
|---|---|---|---|---|
| `adapter.rna` | `f(ctx)` | DataFrame indexed by gene | - | `jurkat.rna_e61` |
| `adapter.proteome` | `f(ctx, total_protein_pg)` | indexed by gene: `uniprot, mass_da, ..., source, copies` | `proteome_table`, `copies` (§4b) | `copies` (legacy `jurkat.geiger_ibaq_v0` kept registered for the equivalence test) |
| `copies.reader` | annotation `f(ctx, spec)`; abundance `f(ctx, spec, annotation)` -> (value, n_rep) | readers of the copies step | `uniprot_reviewed_genes`, `uniprot_ordered_locus`, `annotation_table`, `maxquant_protein_groups`, `paxdb`, `table` | - |
| `adapter.localisation` | `f(ctx, proteome_df)` | indexed by gene: `compartments, loc_source, transmembrane` (+ optional `n_tm, signal_peptide, mito_om`) | `localisation` (core framework, preset), `localisation_table` | `jurkat.hpa_uniprot_itzhak_v0` (legacy reference) |
| `adapter.assembly` | `f(ctx)` | writes `assembly/*` (attribute `OUTPUTS`, or `outputs(cfg)`) | `none`, `complex_table`, `assembly_steps` (core step list + preset, `wholecell/inventory/assembly/`) | `assembly_steps`, preset `legacy_v7` |
| `adapter.secms` | `f(ctx)` | writes `secms/*` (attribute `OUTPUTS`) | `secms_toolkit` (toolkits/secms) | `secms_toolkit`, preset `legacy_v7` |
| `adapter.surface` | `f(ctx)` | DataFrame, `surface.tsv` columns | `surface_proteomics` (§4c) | - (fragment only, OFF) |
| `inventory_builder` | `f(ctx)` | pack-view DataFrame | `generic.monomers` | `jurkat.ingredients_v0` |
| `inventory_rule` | `f(ctx, src_table_path)` | `{"ingredients": df, "changes": df, "extra": {name: df or text}, "report": [..]}` | `generic.signal_lumen` | `jurkat.v4` ... `jurkat.v7` |

`ctx` (`wholecell/inventory/context.py`): `ctx.cfg`, `ctx.params(name)` (= `inventory.params.<name>`), `ctx.frozen(n)`,
`ctx.resolve(v)` (frozen name | `build:<rel>` | path), `ctx.need(rel)` (earlier output), `ctx.write_tsv(df, rel)` (write-once,
atomic). Components declare inputs for Snakemake: attributes `FROZEN`, `PREV` (inventory-relative outputs of earlier steps),
`PARAM_FILES` (params keys holding file references), `OUTPUTS`; every params string naming a frozen input is an input.
Core code holds no species, compartment or cell list; curated lists live in plugin data files with citations.

## 6. Steps (`python -m wholecell inventory --cell C --build B --step S`; `stage.steps(cfg)` gives inputs/outputs)

`rna` -> [`copies`] -> `v0` -> [`secms`] -> `assembly` -> `ingredients` -> `rules:<set>` (one per `inventory.rule_sets` entry, each reading the
previous step's `ingredients.tsv`) -> `pack_table` -> `canonical` -> `compare` (only with `inventory.verify`).

Config keys (`inventory:`): `total_protein_pg` (required; `anchors` when `copies.total_protein.mode: anchors`), `total_protein_source`, `copies` (§4b), `adapters.{rna, proteome, localisation,
secms, assembly, surface}`, `ingredients`, `rule_sets`, `proteome_sources`, `localisation_sources`, `params.<component>`,
`pack_table.order_reference` (optional, §8), `canonical.{key_compartment, abundance_columns, secms_bound_fraction}` (§4c), `verify.{tables, ingredients}`.

## 7. Rules: writing rows, not rewriting files

A rule receives the previous table, returns the new one and one `changes.tsv` row per changed field (old, new, evidence,
the validation check that motivated it). Pack-dependent rules read only declared frozen tables (Jurkat: v5 P/C/N =
`inv_v4pack_loc_compare`, v5 T = `inv_tm_route_libs`, v6 = `inv_v5pack_loc_compare` + `inv_tier2_v5`, hold-out lists).
Rules needing a pre-pass result take it as an explicit input: Jurkat v7 `params."jurkat.v7".lamina_consumption` =
`inv_lamina_consumption_v7` (legacy jlamina output) or `build:preplaced/.../lamina_cell1.json` (lamina rebuilt in the same
build). Geometry-derived constants are config values with their origin (v6 speckle volume 1.963 um3, v7 PCM sphere
cytosolic volume 0.2451 um3) until the geometry stage exports them.

## 8. Exact reproduction of a legacy build

Legacy outputs depend on two arbitrary orders: (a) the scripts ran with a random `PYTHONHASHSEED`, so the iteration order
of string sets decided the order of equal-weight keys in `compartments` / `comp_weights`, the order of some change rows
and evidence texts, and the pick among equal-coverage PDB candidates of a complex; (b) pandas/numpy quicksort is not
stable and was run on another CPU, so rows with equal sort keys came out in another order. Neither carries information.
The port runs with `PYTHONHASHSEED=0` (deterministic) and (Jurkat only) borrows from frozen legacy references what is
otherwise unrecoverable: the PDB pick among exact ties (`params.assembly_steps.inputs.pdb_tiebreak_reference`) and, for the
pack table only, the row order and the key order of comp_weights dicts whose content is exactly equal
(`inventory.pack_table.order_reference`, `legacy_order.py`; content differences are never overwritten and are reported in
`pack_table_order.json`). Result for Jurkat: `pack_table.tsv` md5 = `ingredients_v7.tsv` md5 (d1a0099c...).

## 9. Yeast mapping (for the shared pipeline)

`adapter.proteome` = `copies` (§4b; `config/cells/yeast.yaml`: PaxDb ppm, `NOT_IN_THIS_CELL` categories in
`plugins/yeast/data/inventory_copies.yaml`, G1 size scaling, ribosomes from the cryo-ET density, 3.5 pg fixed; reproduces
build W copies to 4e-16); `adapter.localisation` = UniProt CC -> GO CC ->
Vogtle/Morgenstern/CYCLoPs refinements (`loc_source` tags carry over); `NOT_IN_THIS_CELL` = curated zero list (plugin data);
`adapter.assembly` = Complex Portal + RCSB assemblies + AF2 heterodimers + `oligo_map`; yeast compartment keys stay as
`loc_key`s (`golgi_endosome`, `vacuole`, `cell_wall`, `mito_IM`, ...) mapped by the yeast `pack.region_map`.
