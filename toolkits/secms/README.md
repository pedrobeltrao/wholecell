# toolkits/secms: SEC-MS -> standard assembly-state tables (human | yeast)

One code path and one parameter set (`presets/default.yaml`) turn any cell's SEC-MS data into the tables of
`docs/reference/inventory_schema.md` §4a. Per cell, only **dataset properties** change (format, paths, number of fractions,
void fraction, optional column standards, condition labels). Scope: **human and yeast** (species switch); see the end of
this file for what another species would need.

## Run

- In a build: `inventory.adapters.secms: secms_toolkit` and `inventory.params.secms_toolkit: <run config>` -> step
  `secms` -> `inventory/secms/*`; the assembly adapter reads them via `wholecell.inventory.secms.SecView`. Optional
  `inventory.params.secms_toolkit.combine: <combination spec>` (same keys as `--spec` below) appends the combined rows in
  the build (Step 3, 2026-11-19: config/cells/jurkat.yaml = final2026_jurkat + its combine spec, label `jurkat_final`, read
  by the assembly rules; config/cells/yeast.yaml = update2026_yeast + its combine spec, label `yeast_final`, tables only).
- Stand-alone: `PYTHONHASHSEED=0 python -m toolkits.secms run --config run.yaml --out <new dir>` (paths, placeholders
  `{data_root}` etc. allowed).

## Run config

```yaml
species: human                      # human | yeast  (species/<name>.yaml: Complex Portal taxon, anchor set, hu.MAP allowed)
preset: default                     # default | v8_u1 (v8 / u1 reproduction) | legacy_v7 (Jurkat v7 reproduction only)
overrides: {}                       # deep-merged on the preset; only with a DECISION in PLANNING_NOTES
annotation: {path: <UniProt TSV of the species>}        # columns default to Entry / Mass / Gene Names (primary) / Transmembrane
complexes:
  - {format: complex_portal, path: <ComplexTab export of the species>}      # taxon checked
  - {format: humap, path: <hu.MAP 3.0 complexes>, id_col: HuMAP3_ID, members_col: Uniprot_ACCs}   # human only, membership
pair_candidates: {path: <pairs TSV>, a: <col>, b: <col>, require_true: <bool column>}           # optional
datasets:
  - source_id: <label>
    format: openswath | peptide_long | protein_matrix | protein_long
    path: <file>
    fraction_map: {path: <sample sheet>, skiprows: 0, file_col: 0, fraction_col: 1, condition_col: null, replicate_col: null}  # openswath
    columns: {protein: .., peptide: .., fraction: .., intensity: .., qvalue: .., condition: .., replicate: ..}               # *_long
    value_regex: '^F(?P<fraction>\d+)\.(?P<condition>\w+)\.intensity\.mean$'   # protein_matrix (+ sem_regex)
    n_fractions: 60
    void_fraction: 8                # last fraction of the void volume (column property; null if unknown)
    standards: {path: <fraction, kda[, name]>}   # optional column standards (cross-check of the calibration)
    condition: <label>
    # protein_matrix: a protein with no signal in one run (condition x replicate) is dropped for that run.
    # protein-level options (dataset properties; io.py): duplicate_rows: sum | unique (precursor-level exports that repeat
    # the protein-group quantity, e.g. Spectronaut PG.IBAQ), exclude_ids_regex ('^NA$|^CON_'), id_groups {strip, sep,
    # mode: expand | first | drop} for '[P1;P2]' groups; calibration_from: <earlier source_id> (same column and fractions:
    # reuse its calibration; this dataset's own fits are reported as checks)
```

## Combining conditions (one estimate per cell)

`python -m toolkits.secms combine --run <run dir> --label <cell>_combined --out <dir>` (`combine.py`): per protein and k,
the mean of the condition values, condition SD and range as the uncertainty, one condition where only one scores
the protein; rows `condition = replicate = combined`, selectable by the assembly rules through
`assembly_options.secms`. `--spec <yaml>` (2026-11-17) adds per-dataset weights and roles (`primary` | `fill`: a fill
dataset is used only for proteins no primary dataset scores; confidence column high | discordant | single | low); final
per-cell specs: `datasets/final2026_{jurkat,yeast}_combine.yaml` (Jurkat: lab NT + TCRAct24h + Fossati 2021, equal;
yeast: lab NaCl + Pang 2020 PXD019513 equal, lab YPD fill). **Yeast update (2026-11-17, supersedes the yeast final spec): `datasets/update2026_yeast{,_combine}.yaml`** = lab NaCl + Swamy 2022 PXD028358 (unstressed, 120 mM KCl, no lysate concentration, Superose 6, 4 wild-type profiles; `replicates: split_weight` so the dataset counts once) primary with equal weight; Pang 2020 and lab YPD fill only. Lab datasets 2026 (confidential, local only): `datasets/lab2026_{jurkat,yeast}.yaml`;
analysis in `{data_root}/data/derived/secms_lab2026/`.

## Steps (per dataset, per run = condition x replicate)

| step | module | default (presets/default.yaml) |
|---|---|---|
| input | `io.py` | peptide or protein level; q-value < 0.01; shared (`/`) and decoy peptides removed |
| protein quantification | `quant.py` | top-3 peptides by summed intensity over all runs of the dataset (same peptides in every replicate), summed per fraction; >= 2 peptides; 3-fraction running mean |
| MW calibration | `calib.py` | **curated native-mass anchors of the species** (`data/anchors_native_mass_<species>.tsv`: apex fraction vs native mass), Theil-Sen on log10 MW, intercept = median residual (2026-11-23; before: median(y) - slope x median(x), 1.1-1.8x too high); QC gate n >= 8, Spearman <= -0.7, residual <= x2.5 (fail = stop, no silent fallback); diagnostics: residual SD, leave-one-out error, Complex Portal peak-level calibrants and column standards fitted as cross-checks |
| assembled fraction | `protein.py` | intensity at apparent MW >= k x monomer for k = 2, 3, 4, 6; lo/hi at -/+ 1 calibration SD; void fraction reported; unmeasured (NaN, `measurable` False) when the monomer is within 1 calibration SD of the last fraction's apparent mass or k x monomer above the first non-void fraction (2026-11-23); midpoint share `assembled_frac_mid` at sqrt(k (k - 1)) x monomer (oligomer shares of the assembly rules) |
| peaks | `protein.py` | find_peaks (height 5%, prominence 10% of max, distance 2); profile split at minima; peak class monomer / assembled / void; peak-level assembled fraction |
| complexes | `complexes.py` | co-elution (min of normalised subunit profiles, mean pairwise r); per subunit in-complex fraction at app. MW >= max(0.5 x complex MW, 2 x monomer), void-capped; measured if >= 2 subunits and >= 50% detected |
| pairs | `pairs.py` | Pearson r of smoothed profiles for candidate pairs (optional: all pairs above r) |
| replicates | `protein.py`, `pairs.py`, `pipeline.py` | `replicate = mean` rows with SD and n_rep; matrices with SEM: seeded bootstrap (30 draws) |

Preset `v8_u1` (2026-11-23) = default with the old calibration intercept and no measurable-range rule: it reproduces
`inventory/secms/*` of the Jurkat v8 and yeast u1 builds (checked on the stand-alone runs of the final dataset configs).
`legacy_v7` differs only in: single-peptide proteins kept, no calibration QC gate, no bootstrap, anchor set v1, old
intercept, no measurable-range rule. It reproduces
`sec_assembly.tsv`, `sec_calibrants_peak.tsv` and the orphan `sec_assembly_thresholds.tsv` (a_k) of Jurkat v7 exactly.

## Calibration choice (Pedro, 2026-11-16)

Primary = curated native-mass anchors (v7 approach, Jurkat rho -0.89, residual x1.85). Complex Portal peak-level
calibrants are a diagnostic only (Jurkat rho -0.28 to -0.35, too flat: partial co-elution in mild lysis). Column
standards, when a dataset has them, are a cross-check: on HeLa (Heusel 2020) the anchors put each fraction ~1.7x higher
in MW than the authors' standard calibration (see docs/reference/history/secms_toolkit_report.md).

## Species

| file | human | yeast |
|---|---|---|
| `species/<s>.yaml` | taxon 9606, hu.MAP allowed | taxon 559292, no hu.MAP |
| `data/anchors_native_mass_<s>.tsv` (anchor set v1) | 25 v7 anchors (native_kda); preset legacy_v7 | 26 anchors curated 2026-11-16 the same way (oligomer_n x UniProt mass for homomers; 20S / CCT literature masses); reviewed 2026-11-19 (v3) |
| `data/anchors_native_mass_<s>_v2.tsv` (anchor set v2, 2026-11-17) | v1 minus PFN1, HSP90AB1, YWHAZ, HSPA8 (21) | v1 minus COF1, HSP82, SSA1, BMH1 (22) |
| `data/anchors_native_mass_<s>_v3.tsv` (anchor set v3, 2026-11-19) | - (human stays v2) | v2 minus PDC1 (pH-dependent dimer/tetramer), EFT1 (EFT1/EFT2 one UniProt entry, never matched) (20); review with literature states and per-dataset residuals: `data/ANCHORS_yeast_v3_2026-11-19.md` |
| preset default `calibration.anchor_set` | v2 | v3 (`{human: v2, yeast: v3}`; a plain string applies to every species) |

Another species would need: a `species/<name>.yaml` (taxon, annotation columns), a Complex Portal (or equivalent)
export for its taxon, and a curated anchor table of >= ~20 soluble, stable assemblies from ~15 kDa to ~1 MDa with
citations. Nothing else in the code is species-specific; the species list is deliberately closed (`pipeline.SPECIES`).
