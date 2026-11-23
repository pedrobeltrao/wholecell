# Stage 4: assembly

**Purpose.** Decide which proteins are packed as complexes or oligomers, with which structure and how many copies,
and which remain monomers. The result is the species list the packer reads.

**Today.** Part of the `inventory` stage: steps `secms` (optional), `assembly`, `ingredients`, then `pack_table` and
`canonical` at the end of the inventory chain.

**Inputs.** `inventory_v0.tsv`; complex sources (Complex Portal, PDB assemblies via SIFTS, AFDB homo- and
heterodimer tables, STRING links for Jurkat; a complex table for the testcell); SEC-MS data processed by
`toolkits/secms` (Jurkat: Fossati 2021, 25 native-mass anchors; testcell: synthetic).

**Outputs** (`inventory/`): `secms_protein.tsv`, `secms_complex.tsv` (SEC-MS step); `assembly_steps`:
`assembly/{complexes_v0, inventory_v0_complexed, complex_copies_rule, sec_assembly, sec_calibrants_peak, homo_oligomers_v0,
afdb_homodimers_v0, afdb_heterodimers_v0}.tsv` + `assembly_options.json`, `assembly_steps.json` (generic adapters:
`assembly/complexes.tsv`, `assembly/complexed.tsv`); `ingredients/ingredients.tsv`; `pack_table.tsv` (the table the packer
reads); canonical `species.tsv`, `sources.tsv`, `version_changes.tsv`.

**Step list (core, `wholecell/inventory/assembly/`; A4).** `adapters.assembly: assembly_steps` runs one ordered list of
plain step functions (`steps.py`, registry `STEPS` = dict): `curated_machines` -> `complex_portal` -> `secms_fractions`
-> `afdb_homodimers` -> `afdb_heterodimers` -> `pdb_homo_oligomers` (`DEFAULT_ORDER`). The complex steps come first;
after them `complexes_v0.tsv`, `inventory_v0_complexed.tsv`, `complex_copies_rule.tsv` are written and the later steps
work on the residual pool. Behaviour (which steps, order, thresholds, allocation rules) comes from a named preset
(`presets/<name>.yaml`); cell data from the config: `inputs` (frozen input names), `data` (curated TSVs with a `reason`
column + scalars). Config errors name the key, the allowed values and the nearest valid name.

| Preset | Status | Notes |
|---|---|---|
| `legacy_v7` | used by Jurkat | v7 exactly (B3/B4 md5); runs `pdb_homo_oligomers` before the AFDB dimers, as v7 did |
| `legacy_x5` | used by yeast (A9) | X5 exactly: every s0 assembly table identical (`scripts/assembly_x5_check.py`); `pools: ingredients` (`ingredient_steps.py`), yeast lists in `plugins/yeast/data/assembly/*.tsv`, yeast plugin steps + structure hook (`plugins/yeast/assembly.py`) |

**Pools.** `legacy_v7` works on integer pools of `inventory_v0.tsv` (`steps.py`). `legacy_x5` (`pools: ingredients`,
`ingredient_steps.py`) works on one ingredient table with fractional copies: `curated_machines` (TSV machines with
paralog families; copy rules median / family / hook; Nup free-pool cap; legacy oligomers) -> plugin
`yeast_preribosomes` -> `complex_portal` (skip if a subunit is used except data.shared_subunits; mtDNA subunits out of
the copies rule; RCSB cache + local assembly1: >= 0.8 mass, >= 0.6 residues, < 60k residues, by resolution, preferred
entries, paralog templates; PDB chain-count stoichiometry; orientable TM) -> plugin steps (free r-proteins dropped, ER
ribosomes, mitoribosomes, III2IV2, Arp2/3, F-actin) -> `afdb_heterodimers` (file order, ipTM >= 0.5, pool >= 20,
paralog / Hsp70 / histone / compartment filters, proportional allocation, bound fraction 1.0 / 0.5) ->
`pdb_homo_oligomers` (`copies: none`: audit TSV, UniProt state, AFDB fallback). Outputs `assembly/{ingredients.csv,
oligo_map.json, complexes_tried.csv, complexes_added.csv, heterodimers.csv, oligomers.csv, f_actin.json}`. Plugin steps
are registry kind `assembly_step`; hooks (copy counts, the structure builder: beads, own RNA, membrane orientation)
kind `assembly_hook`, named in the preset options.

**Config keys** (`inventory:`): `adapters.{secms, assembly}`, `ingredients` (registered builder),
`params.secms_toolkit.{species, preset, annotation, complexes, pair_candidates, datasets}`,
`params.assembly_steps.{preset, overrides, inputs, data, workers}` (`overrides` may only use keys the preset has;
`complex_portal.copies: stoichiometry_cap | measured_free_pool`), for the generic `complex_table` adapter
`params.complex_table.assembly_options`, `pack_table.order_reference`, `canonical.*`. Spec: [../reference/inventory_schema.md](../reference/inventory_schema.md)
§3-§4a.

**Checks.** SEC-MS regression S1-S11 (`scripts/secms_regression_check.py`;
[../reference/history/secms_regression_criteria.md](../reference/history/secms_regression_criteria.md)); Jurkat
pack table md5 (B3). `measured_free_pool` stays off for Jurkat (see
[../reference/history/secms_toolkit_report.md](../reference/history/secms_toolkit_report.md)).

**Code.** `wholecell/inventory/assembly/{__init__.py (options, presets, runner), steps.py, presets/}`,
`wholecell/inventory/{secms.py, generic.py, canonical.py, legacy_order.py}`, `toolkits/secms/`; Jurkat curation
`plugins/jurkat/data/assembly/*.tsv`. Yeast: `plugins/yeast/assembly.py`, `plugins/yeast/data/assembly/*.tsv`. Test: `tests/test_assembly_steps.py`; check `scripts/assembly_x5_check.py`.

**Jurkat vs yeast today.** Jurkat runs the core step list (preset `legacy_v7`), yeast the same step list with
`legacy_x5` (A9): its stage-1 protein table is the core copies + localisation in the engine's column layout (plugin
step `yeast_protein_table`, `before: [inventory_assembly]`), the yeast plugin steps and structure hooks call the vendored
modules in the build's `engine/` layout (structure hooks add library entries to `engine/lib`, a writable copy of the
sealed library: plugin step `yeast_engine_lib`), and the builder `yeast.x5_ingredients` turns `assembly/ingredients.csv`
into the pack table. Exactness: `scripts/assembly_x5_check.py` (A4) and the A9 stage check. A9 added the plugin step
`yeast_barrel_monomers` (s1 rebuilt the outer-membrane beta-barrel monomers as oriented membrane library entries before
anything else; A4 compared tables only and missed it).

**Methods text (draft).** Proteins were grouped into complexes using Complex Portal entries with an experimental
structure, then AlphaFold-predicted homo- and heterodimers. The fraction of each subunit in assembled form was taken
from size-exclusion chromatography coupled to mass spectrometry (SEC-MS), calibrated on native-mass anchors.
Remaining copies were packed as monomers.
