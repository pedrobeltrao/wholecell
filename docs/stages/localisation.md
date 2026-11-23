# Stage 3: localisation

**Purpose.** Give every species a set of compartment weights (`comp_weights`, e.g. 0.7 cytosol, 0.3 nucleus) that the
packer turns into regions.

**Today.** Part of the `inventory` stage. Initial weights are written in step `v0` by the core localisation framework
(adapter `localisation`, `wholecell/inventory/localisation/`, since A5 2026-11-16):

1. **Evidence.** Source adapters (`evidence.SOURCES`: `hpa_subcellular`, `uniprot_cc`, `fractionation_nc`,
   `localisation_table`; for yeast `uniprot_cc_statements` (statement order, Huh 2003 precedence), `go_cc_terms`,
   `class_columns` (Vögtle 2017), `proteome_classes` (Morgenstern 2017), `imaging_locscore` (CYCLoPs, with r-protein /
   prenylation exclusions)) write one long table `localisation/evidence.tsv` (`gene, uniprot, source, evidence_class,
   compartment, score, reliability, n_obs, note`); term maps are cell TSVs with a `reason` column.
2. **Combiner** (`combiner.STEPS`, one function per step, order from the preset): `drop_unreliable`,
   `restrict_membrane`, `split` (score | rank | equal), `combine` (blend | precedence by class), `new_from_imaging`,
   `round`, `nc_model` (fraction with background b | ratio_exp), `membrane_mechanism`; ordered-annotation steps (yeast):
   `ordered_annotation` (primary / fallback sources, per-gene key resolution, generic-key refinement, membrane
   fallbacks), `annotation_overrides`, `mito_fractionation` (measured sub-compartment, proteome veto, `mito_hc` flag),
   `rank_weights` (1, ½, ¼ by list order, secondary down-weighting rules, volume pool r = 1), `imaging_groups` (group
   shares from LOC-scores, annotation within groups, new groups from 0.5 with flags, N/C ratio exp(k(N − C)))
3. **Overrides.** Curated revisions are rows of a TSV (`rule_set, rule_id, key_type, key, field, op, value, reason,
   citation`) applied by one function, `overrides.apply` (field `annotation` = the ordered key list, used inside the
   combiner, e.g. yeast LOC_FIX).

Ordered rule sets (`rules:<set>`) then revise the weights. Each rule set writes a new folder and never edits the
previous one. The rule sets also make some non-localisation changes (see `changes.tsv`, column `rule_id`).

**Presets** (`params.localisation.preset`; `presets/<name>.yaml`): `legacy_v7` reproduces Jurkat v7 exactly
(inventory_v0.tsv, every rule table, pack_table md5 d1a0099c). Its `rules` section lists every step of the Jurkat rule
sets v4..v7 in order and what it is: `overrides` (v4.c, v4.f_curated, v4.f_drop_mito, v4.f_clic, v7.R4b, v7.F1: rows of
`plugins/jurkat/data/localisation/overrides_legacy_v7.tsv`), `frozen_overrides` (v5.P_nuclear_pool, v6.R1_nc: override
tables generated at build time from the frozen v4/v5 pack tables, written to `rules/<set>/overrides_<id>.tsv`; only in
this preset, plan §3b) or `legacy` (a rule function still in code, with a one-line reason). The rule-set code checks its
step order against this list. `legacy_x5` reproduces the yeast X5 comp_weights exactly (A5 yeast half, 2026-11-17:
6,067 proteins identical to the port build s0, `scripts/dev/localisation_x5_check.py`, results in
`yeast/_dev_unify/a5/`); its cell data are `plugins/yeast/data/localisation/` (term rules, SPB core, CYCLoPs map, LOC_FIX
override row); the `inventory.params.localisation` block and its frozen inputs are in `config/cells/yeast.yaml` (A9).

**Inputs.** Localisation sources declared in `frozen_inputs`. Jurkat: HPA subcellular, UniProt, Itzhak 2016
fractionation (adapter `jurkat.hpa_uniprot_itzhak_v0`), then rule sets `jurkat.v4` ... `jurkat.v7`. Some v5/v6 rules
read results of earlier packs; these are declared frozen tables. Testcell: a synthetic table plus
`generic.signal_lumen`.

**Outputs** (`inventory/`): `localisation/{evidence.tsv, gene_features.tsv, localisation_options.json}`;
`v0/inventory_v0.tsv`; `rules/<set>/{ingredients.tsv, changes.tsv, report.txt}`;
`changes_all.tsv`; canonical `species_loc.tsv` (long form with `rule_id`), `localisation_sources.tsv`,
`field_provenance.tsv`.

**Config keys** (`inventory:`): `adapters.localisation` (`localisation`), `params.localisation.{preset, overrides,
label, genes, sources, data}`, `localisation_sources`, `rule_sets`, `canonical.key_compartment`. The legacy adapter
`jurkat.hpa_uniprot_itzhak_v0` stays registered as the readable reference of the v0 arithmetic. The mapping of keys to
packing regions is in `pack.region_map` (see [packing.md](packing.md)).

**Checks.** Validation step `loc_compare` (placed shares vs OpenCell, HPA consensus, HPA Jurkat, plus a frozen
held-out half; [validation.md](validation.md)). Jurkat B3/B4: `pack_table.tsv` md5 equals the legacy
`ingredients_v7.tsv` (d1a0099c).

**Code.** `wholecell/inventory/localisation/{__init__, evidence, combiner, overrides}.py`, `presets/`;
`plugins/jurkat/inventory/rules_v4.py` ... `rules_v7.py`; tests `tests/test_localisation.py`; yeast check
`scripts/dev/localisation_x5_check.py`.

**Jurkat vs yeast today.** Yeast takes its localisation (UniProt/GO with Huh precedence, Vögtle/Morgenstern, CYCLoPs,
NUC_K) from the core with `legacy_x5` (A9; same weights as the engine, A5). Its evidence is keyed by ORF while the
copies table is keyed by gene name: `genes.key: orf` names the proteome-table column that holds the evidence key. The
CYCLoPs .xls screens are read with xlrd (in the core env vcell). **B5** switches Jurkat to common
defaults with no pack-dependent rules (plan §3).

**Methods text (draft).** Each protein received compartment weights from database annotations and proteomic
fractionation. For Jurkat, HPA and UniProt annotations were blended with equal weight where both existed, and the
nuclear/cytoplasmic split of soluble proteins was taken from fractionation data where at least three profiles were
available. Curated rules then revised specific groups; every revision was recorded with its evidence and rule
identifier.
