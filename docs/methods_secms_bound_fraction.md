# SEC-MS bound fraction: what the structural model represents and what SEC-MS measures

Attribute `secms_bound_fraction` (canonical inventory table `inventory/secms_bound_fraction.tsv`, definition, columns and
caveats in [reference/inventory_schema.md](reference/inventory_schema.md) §4c; code `wholecell/inventory/secms_bound.py`).
DECISION (Pedro, 2026-11-19, option 2): SEC-MS measures much more assembly than the structural model represents; the
model does not invent structure-less assemblies and the packing is unchanged; the measured-but-unrepresented share is
exported per protein for simulations to use.

## 1. Methods paragraph (paper draft)

**Assembly state: structural model vs size-exclusion chromatography.** The structural model places a protein in an
assembly only when a structure of that assembly is available or can be built: Complex Portal complexes with a PDB
template, homo-oligomers from PDB biological assemblies, and AlphaFold homo- and heterodimers that pass interface-
confidence gates; the share of each protein's copies put into an oligomer or dimer is taken from SEC-MS (the share of
its elution profile at apparent mass closer, on the log scale, to the oligomer than to the next smaller state:
>= sqrt(k (k - 1)) x monomer mass). All other copies
are packed as monomers. SEC-MS reports a different quantity: for each protein, the fraction of its elution profile at
apparent molecular mass at least k times its monomer mass (a_k; k = 2 for monomers, complex subunits and dimers, k = 3
or 4 for proteins modelled in larger homo-oligomers, as in the assembly rules). Apparent masses come from a robust
(Theil-Sen) line of log10 native mass against apex fraction through 20-21 curated anchor proteins per species, with the
intercept set so that the median anchor residual is zero; a_k is reported only where the column resolves it (monomer
mass at least one calibration residual SD above the apparent mass of the last fraction, k x monomer below the first
non-void fraction), and is otherwise "unmeasured" (NA). a_k is computed from calibrated SEC-MS runs combined
across datasets with equal weight (Jurkat: two unpublished lab runs, resting and 24 h TCR-activated, and Fossati et al.
2021; yeast: an unpublished lab run under NaCl stress and Swamy et al. 2022, with a lab YPD run and Pang et al. 2020 as
fill). SEC-MS sees every assembly that survives lysis and dilution, including those for which no structure exists
(translation factors on ribosomes, chaperones bound to clients, RNA-binding proteins on mRNPs), and it also sees
artefactual co-elution (aggregates, nucleic-acid or membrane fragments, detergent micelles). We therefore did not add
assemblies without structural evidence: the partners and stoichiometry of such assemblies are unknown, a placeholder body
of the apparent mass would fix a shape and a composition the data do not support, and SEC-MS cannot separate a specific
complex from lysis or column artefacts protein by protein. Instead, for every protein we report the measured assembled
fraction a_k with its uncertainty (dataset standard deviation and a calibration band obtained by shifting the
apparent-mass scale by one calibration residual standard deviation), the modelled assembled fraction m = 1 - (copies
packed as monomers) / (copies in the cell), whether the protein's assembly is represented by a structure, and the
unrepresented bound fraction u = max(0, a_k - m). [Numbers to be re-taken from the v9 / u2 builds' summary files: on
the v8 / u1 model states the 2026-11-23 estimator gives a mass-weighted a_k of 0.67 of the SEC-covered mass in both
cells (0.785 Jurkat / 0.809 yeast before), and 29 % (Jurkat) / 22 % (yeast) of the protein mass measured bound but not
represented (35 % / 28 % before); §5.] This attribute is exported with every
simulation box: each monomer molecule carries u / (1 - m), the probability that a free copy of that protein is, in the
cell, part of an assembly the model does not resolve.

## 2. How simulations may use the attribute (not implemented in the pipeline)

The box export (`export/boxes/<box>/molecules.parquet`, [reference/simulation_workflows.md](reference/simulation_workflows.md)
§3) adds two columns to every monomer molecule: `secms_bound_fraction` (a_k of its protein) and
`secms_unrepresented_share_of_free` (s = u / (1 - m)); assemblies, pre-pass structures and unmeasured proteins carry NaN.
`box.json` `secms_bound_fraction` records the table checksum, the SEC-MS run, the k values and the intended use. The
stub diffusion coefficients (`R_h_nm`, `D_t_um2_s`, `D_r_per_us`) do NOT use it. Intended uses, to be chosen and
pre-registered by each application:

- **Two-state mixture.** Treat each monomer molecule as free with probability 1 - s and bound with probability s
  (static assignment drawn per molecule with a recorded seed, or exchange between the states with a residence time when
  one is known). Do not change copy numbers or positions: the bound share stays at its packed position.
- **Effective hydrodynamic size of the bound share.** Use the SEC apex (`apex_app_kda` in the canonical table) as the
  apparent mass of the unresolved assembly: R_h of a globular particle of that mass (minimal sphere radius
  R_min = 0.066 x M^(1/3) nm with M in Da, Erickson 2009; Stokes radius ~1.2-1.3 x R_min) or a crowder bead of that radius tethered to the molecule; D_t, D_r by Stokes-Einstein. Bounded by the
  calibration band (`secms_bound_fraction_lo`, `_hi`), not by the dataset SD alone.
- **Reduced mobility.** For proteins whose partners are immobile on the simulation time scale (ribosome-bound
  translation factors, cytoskeleton or membrane binders), the bound share can be modelled as immobile or slowed by a
  factor, as a sensitivity arm rather than a default.
- Sensitivity: report results with s = 0 (the structural model as packed) and with s from the table, and, where it
  matters, at k = 2 vs the rules' k and with the calibration band, since a_k depends on k and on lysis / column
  conditions (§4c caveats).

## 3. Per-cell summary (inventory-only dry runs)

**Superseded numbers (2026-11-23).** The tables below are from the v8 / u1 estimator (calibration intercept
median(y) - slope x median(x), no measurable-range rule; toolkit preset `v8_u1`). The default since 2026-11-23 lowers the
measured a_k of the SEC-covered mass to 0.666 (Jurkat) and 0.673 (yeast) and leaves 104 Jurkat / 171 yeast proteins
unmeasured; analysis and new numbers in §5. Re-take this section from the v9 / u2 builds.

Inventory-only dry runs (CPU <= 6) of the branch code (commit a65b42d) into new folders
`wholecell_builds/{jurkat,yeast}/_dev_unify/secms_bound_fraction/new_a65b42d_<config>/` (Jurkat: full inventory chain,
geometry em_cell1_v8r copied from `jurkat/_dev_unify/step3_secms/new_a12bdac_wiring`; yeast: pseudo-build as
`scripts/dev/step3_yeast_secms_dryrun.py`, inventory/{yeast, v0, copies} + engine copied from build
`yeast/20261119-87a2018-s0-step2ds`, then secms -> assembly (preset x5_secms) -> ingredients -> pack_table -> canonical).
Numbers below come from `secms_bound_fraction.tsv` of these runs (mass = copies x monomer mass; shares of the cell's
protein mass unless stated); per-cell files `summary_<config>.json` and `top20_<config>.tsv` sit next to the runs.
The Jurkat v7 port (Fossati only) and yeast_x5_core (no SEC-MS step: every measured column NA, model columns filled)
were run for the unchanged check (`unchanged_check.json`).

| | Jurkat (`jurkat.yaml`, SEC-MS `jurkat_final`) | Yeast (`yeast.yaml`, SEC-MS `yeast_final`) |
|---|---|---|
| protein mass (pg) / proteins with copies | 100.00 / 13,804 | 3.61 / 5,759 |
| proteins measured by SEC-MS (share of mass) | 5,285 (91.9 %) | 4,355 (99.4 %) |
| mass with a modelled assembly (structure_represented) | 62.7 % | 67.4 % |
|   of which held by family / pre-pass species (not_in_pack) | 4.6 % | 0.0 % |
| modelled assembled share (sum m x mass) | 43.7 % | 59.1 % |
| mass packed only as monomers (monomer_only) | 37.3 % | 32.6 % |
|   SEC-measured a_k of that mass | 72.9 % (of 34.5 % of the cell measured) | 74.9 % (of 32.0 % of the cell measured) |
| SEC-covered mass: measured a_k vs modelled | 78.5 % vs 43.2 % | 80.9 % vs 59.4 % |
| **measured bound but unrepresented** (sum u x mass) | **35.4 %** of the cell (35.379 pg; 38.5 % of the SEC-covered mass) | **27.7 %** of the cell (0.997 pg; 27.8 % of the SEC-covered mass) |
|   of which in monomer-only proteins | 25.1 % | 24.0 % |
| k used (measured proteins) | k2: 5,168, k3: 87, k4: 30 | k2: 4,319, k3: 23, k4: 13 |

Top 20 contributors, Jurkat (u x mass; share of the cell's protein mass; measured a_k +/- SD [calibration band], modelled m, state):

| # | gene | mass (pg) | state | k | a_k +/- SD [lo, hi] | conf. | m | u | u x mass (% cell) | apex / mono (kDa) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | HSP90AB1 | 1.4725 | in_assembly | 2 | 0.96 +/- 0.05 [0.90, 0.99] | high | 0.00 | 0.96 | 1.41 | 836 / 83 |
| 2 | MYH9 | 0.6930 | monomer_only | 2 | 0.97 +/- 0.02 [0.94, 0.99] | high | 0.00 | 0.97 | 0.67 | 3420 / 227 |
| 3 | HSPA8 | 0.6620 | in_assembly | 2 | 0.84 +/- 0.07 [0.70, 0.92] | high | 0.06 | 0.77 | 0.51 | 251 / 71 |
| 4 | H1-2 | 0.4938 | monomer_only | 2 | 1.00 [1.00, 1.00] | single | 0.00 | 1.00 | 0.49 | 7168 / 21 |
| 5 | ACTB | 1.9932 | in_assembly | 2 | 0.63 +/- 0.14 [0.33, 0.88] | high | 0.39 | 0.24 | 0.48 | 87 / 42 |
| 6 | HMGB1 | 0.5043 | monomer_only | 2 | 0.92 +/- 0.06 [0.64, 0.98] | high | 0.00 | 0.92 | 0.46 | 178 / 25 |
| 7 | SPTAN1 | 0.3603 | monomer_only | 2 | 0.91 +/- 0.04 [0.87, 0.94] | high | 0.00 | 0.91 | 0.33 | 8503 / 285 |
| 8 | TUBB | 0.3637 | monomer_only | 2 | 0.89 +/- 0.10 [0.76, 0.96] | high | 0.00 | 0.89 | 0.32 | 289 / 50 |
| 9 | MSN | 0.4017 | monomer_only | 2 | 0.78 +/- 0.18 [0.53, 0.98] | discordant | 0.00 | 0.78 | 0.31 | 394 / 68 |
| 10 | H1-5 | 0.3182 | monomer_only | 2 | 0.98 +/- 0.01 [0.96, 0.99] | high | 0.00 | 0.98 | 0.31 | 4792 / 23 |
| 11 | TUBB4B | 0.9317 | in_assembly | 2 | 0.90 +/- 0.11 [0.71, 0.97] | high | 0.57 | 0.33 | 0.31 | 168 / 50 |
| 12 | HNRNPU | 0.3146 | monomer_only | 2 | 0.97 +/- 0.03 [0.93, 1.00] | high | 0.00 | 0.97 | 0.30 | 6570 / 91 |
| 13 | NUCLEOLIN | 0.3973 | monomer_only | 2 | 0.77 +/- 0.25 [0.49, 0.94] | discordant | 0.00 | 0.77 | 0.30 | 405 / 77 |
| 14 | HSP90AA1 | 0.5606 | in_assembly | 2 | 0.98 +/- 0.01 [0.97, 0.99] | high | 0.44 | 0.54 | 0.30 | 735 / 85 |
| 15 | HMGB2 | 0.3096 | monomer_only | 2 | 0.90 +/- 0.07 [0.60, 0.98] | high | 0.00 | 0.90 | 0.28 | 245 / 24 |
| 16 | HNRNPL | 0.2695 | monomer_only | 2 | 0.83 +/- 0.19 [0.63, 1.00] | discordant | 0.00 | 0.83 | 0.22 | 7830 / 64 |
| 17 | SPTBN1 | 0.3174 | monomer_only | 2 | 0.69 +/- 0.45 [0.17, 0.95] | discordant | 0.00 | 0.69 | 0.22 | 3103 / 275 |
| 18 | HNRNPA2B1 | 0.3979 | monomer_only | 2 | 0.54 +/- 0.37 [0.27, 0.95] | discordant | 0.00 | 0.54 | 0.21 | 189 / 37 |
| 19 | LMNB1 | 0.2293 | monomer_only | 2 | 0.91 +/- 0.11 [0.79, 1.00] | high | 0.00 | 0.91 | 0.21 | 2406 / 66 |
| 20 | MYO1G | 0.2181 | monomer_only | 2 | 0.96 +/- 0.03 [0.93, 0.98] | high | 0.00 | 0.96 | 0.21 | 17778 / 116 |

Top 20 contributors, Yeast (u x mass; share of the cell's protein mass; measured a_k +/- SD [calibration band], modelled m, state):

| # | gene | mass (pg) | state | k | a_k +/- SD [lo, hi] | conf. | m | u | u x mass (% cell) | apex / mono (kDa) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | YEF3 | 0.0687 | monomer_only | 2 | 0.79 +/- 0.20 [0.58, 0.94] | discordant | 0.00 | 0.79 | 1.51 | 1522 / 116 |
| 2 | TEF1; TEF2 | 0.0404 | in_assembly | 2 | 0.80 +/- 0.26 [0.49, 0.97] | discordant | 0.00 | 0.79 | 0.89 | 911 / 50 |
| 3 | SSA2 | 0.0370 | monomer_only | 2 | 0.84 +/- 0.05 [0.60, 0.96] | high | 0.00 | 0.84 | 0.86 | 240 / 69 |
| 4 | SSA1 | 0.0340 | monomer_only | 2 | 0.84 +/- 0.05 [0.59, 0.95] | high | 0.00 | 0.84 | 0.80 | 240 / 70 |
| 5 | PGK1 | 0.0565 | monomer_only | 2 | 0.40 +/- 0.15 [0.14, 0.86] | discordant | 0.00 | 0.40 | 0.63 | 82 / 45 |
| 6 | SSB1 | 0.0224 | monomer_only | 2 | 0.75 +/- 0.14 [0.57, 0.96] | discordant | 0.00 | 0.75 | 0.47 | 192 / 67 |
| 7 | SSB2 | 0.0208 | monomer_only | 2 | 0.77 +/- 0.12 [0.59, 0.97] | high | 0.00 | 0.77 | 0.45 | 209 / 67 |
| 8 | SSE1 | 0.0141 | monomer_only | 2 | 0.86 +/- 0.08 [0.70, 0.97] | high | 0.00 | 0.86 | 0.34 | 318 / 77 |
| 9 | PAB1 | 0.0112 | monomer_only | 2 | 0.97 +/- 0.00 [0.95, 0.98] | high | 0.00 | 0.97 | 0.30 | 1796 / 64 |
| 10 | URA2 | 0.0111 | monomer_only | 2 | 0.90 +/- 0.06 [0.82, 0.98] | high | 0.00 | 0.90 | 0.28 | 3570 / 245 |
| 11 | PSA1 | 0.0105 | monomer_only | 2 | 0.86 +/- 0.06 [0.70, 0.97] | high | 0.00 | 0.86 | 0.25 | 151 / 40 |
| 12 | KAR2 | 0.0100 | in_assembly | 2 | 0.86 +/- 0.03 [0.77, 0.93] | high | 0.03 | 0.84 | 0.23 | 1493 / 74 |
| 13 | SSC1 | 0.0090 | monomer_only | 2 | 0.81 +/- 0.05 [0.60, 0.92] | high | 0.00 | 0.81 | 0.20 | 230 / 71 |
| 14 | SAH1 | 0.0111 | monomer_only | 2 | 0.64 +/- 0.30 [0.40, 0.98] | discordant | 0.00 | 0.64 | 0.20 | 135 / 49 |
| 15 | CHC1 | 0.0072 | monomer_only | 2 | 0.99 +/- 0.00 [0.98, 1.00] | high | 0.00 | 0.99 | 0.20 | 9691 / 187 |
| 16 | HSP60 | 0.0067 | monomer_only | 2 | 0.95 +/- 0.02 [0.92, 0.98] | high | 0.00 | 0.95 | 0.18 | 1136 / 61 |
| 17 | EFT1; EFT2 | 0.0164 | monomer_only | 2 | 0.38 +/- 0.17 [0.18, 0.59] | discordant | 0.00 | 0.38 | 0.17 | 74 / 93 |
| 18 | ACS2 | 0.0074 | monomer_only | 2 | 0.83 +/- 0.10 [0.61, 0.97] | high | 0.00 | 0.83 | 0.17 | 255 / 75 |
| 19 | GDH1 | 0.0063 | monomer_only | 2 | 0.99 +/- 0.01 [0.95, 1.00] | high | 0.00 | 0.99 | 0.17 | 367 / 50 |
| 20 | RPN1 | 0.0080 | in_assembly | 2 | 0.96 +/- 0.02 [0.82, 1.00] | high | 0.23 | 0.74 | 0.16 | 720 / 109 |

Notes.
- "Mass with a modelled assembly" counts every protein named in at least one assembly species with copies (plus the
  family-held Jurkat histones); a protein can be named in an assembly and still have most copies free (HSP90AB1 m = 0.002,
  yeast TEF1/TEF2 m = 0.00): the copy-weighted modelled share is the row below it. The yeast monomer-only mass (32.6 %)
  is the addendum's "~34 % kept fully monomeric" under this naming rule; its SEC a_2 (74.9 %) matches the addendum's 0.75.
- Jurkat pre-placed structures: microtubules, spectrin / lamin / vimentin networks and H1 on chromatin consume copies at
  the pre-placed stage, after the inventory, so an inventory-only m counts them as free. Proteins of the Jurkat
  procedurally built polymer list (`assembly_steps.data.homo_oligomer_exclude_regex`) plus H1 variants: 60 proteins,
  5.2 % of the mass, 2.5 of the 35.4 points of unrepresented bound mass (H1-2, SPTAN1, TUBB, H1-5, TUBB4B, SPTBN1, LMNB1,
  VIM). The remaining ~33 points are chaperones (HSP90, HSC70), myosins, hnRNPs / HMGs, nucleolin and other proteins
  without a structure of their assembly.
- Jurkat SEC-MS covers 5,285 of 13,804 proteins with copies (91.9 % of the mass); unmeasured proteins are NA
  (`na_reason`), never 0, and are left out of the "of the SEC-covered mass" shares.
- Top contributors with `confidence = discordant` (datasets differ by > 0.3; e.g. MSN, NUCLEOLIN, hnRNPs, yeast YEF3,
  TEF1/2, PGK1) and wide calibration bands should be used with their bands, not the point value. PGK1 (yeast, a_2 0.40,
  apex 1.84 x monomer) is the k = 2 artefact of inventory_schema.md §4c (a monomer reading > 0 near the cut).
- Jurkat v7 port (Fossati 2021 only, preset legacy_v7): 3,040 proteins measured (78.6 % of the mass), modelled assembled
  42.9 %, measured bound but unrepresented 28.4 % of the protein mass (one dataset, mild lysis; not the reported value).

## 4. Unchanged check (branch base a77c683 vs a65b42d)

Same configs, same inputs, inventory-only, every file under `inventory/` md5-compared
(`_dev_unify/secms_bound_fraction/unchanged_check.json` per cell):

| config | files compared | identical | differ only by the run-folder path written inside | new files only | pack_table md5 (base = new) |
|---|---|---|---|---|---|
| jurkat | 72 | 67 | 5 (`compare_legacy.json`, `rules/jurkat.v{4..7}/report.txt`) | secms_bound_fraction.tsv, _summary.json | c073b868... |
| jurkat_v7_port | 72 | 67 | 5 (same) | same | d1a0099c... (= legacy ingredients_v7.tsv) |
| yeast | 44 | 44 | 0 | same | 7c143a7d... |
| yeast_x5_core | 35 | 35 | 0 | same (all measured columns NA: no SEC-MS step) | 09a4ad23... |

yeast_x5_port: the core inventory stage does not run for it (only the copies step is ported; `adapters.localisation:
null`); its vendored path never calls `canonical`, so it cannot change. Tests (scripts): test_secms_bound (new),
test_secms, test_assembly_steps, test_inventory, test_localisation, test_pack_switches, test_validate_suite: all pass.

## 5. Estimator revision (2026-11-23, rebuild v9)

Three changes to the SEC-MS assembled fraction (toolkits/secms; old rules: preset `v8_u1`, which reproduces
`inventory/secms/*` of the v8 / u1 builds exactly):

1. **Calibration intercept.** The Theil-Sen intercept was median(y) - slope x median(x) (scipy's default), which does
   not pass through the bulk of an unevenly spaced anchor set: in the v8 / u1 fits the median anchor residual was
   -0.05 to -0.26 log10, i.e. every apparent mass was 1.1-1.8x too high, by a different factor per dataset. The intercept
   is now median(y - slope x) (median anchor residual 0; `calibration.intercept: median_residual`). The slope is
   unchanged.
2. **Measurable range.** a_k is computed only when the monomer mass is at least 10^(1 calibration SD) x the apparent mass
   of the last fraction and k x monomer at most the first non-void fraction's apparent mass / 10^(1 SD)
   (`assembled.range_margin_sd: 1.0`); otherwise a_k is NA ("unmeasured"), not 1 by construction (before: every protein
   with 2 x monomer below the last fraction scored 1, e.g. yeast lab NaCl with 36.9 kDa at the last fraction).
3. **Midpoint share for oligomer shares.** a_k puts the boundary on the k-mer itself, so a pure dimer reads ~0.5 at
   k = 2 on an unbiased calibration. The assembly rules of `jurkat.yaml` and `yeast.yaml` take homo-oligomer and AF
   homodimer shares from `assembled_frac_mid` (intensity at >= sqrt(k (k - 1)) x monomer, share source
   `secms_midpoint`). The bound fraction attribute stays a_k.

Effect on the v8 / u1 SEC-MS inputs (stand-alone toolkit runs of the final dataset configs, combined with the build
specs; masses, copies and rule k from the builds' `secms_bound_fraction.tsv`):

| | Jurkat old | Jurkat new | yeast old | yeast new |
|---|---|---|---|---|
| proteins measured / unmeasured | 5,285 / 0 | 5,181 / 104 | 4,355 / 0 | 4,184 / 171 |
| a_k, mass-weighted (SEC-covered mass) | 0.785 | 0.666 | 0.809 | 0.673 |
| a_k, protein mean (median) | 0.830 (0.961) | 0.706 (0.873) | 0.876 (0.977) | 0.775 (0.924) |
| curated monomeric anchors, mean a_2 | 0.35 | 0.15 | 0.45 | 0.14 |
| cytosolic ribosomal proteins, mean a_k | 0.97 | 0.97 | 0.98 | 0.96 |

Apex-based (share of mass whose profile apex is >= k x monomer) and peak-based estimators on the new calibration give
0.66-0.69 in both cells; correcting the old or new values for what monomeric anchors read gives 0.61-0.67.

