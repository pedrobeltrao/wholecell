# Yeast SEC calibration anchors: review and set v3 (2026-11-19)

Review requested by Pedro (Step 3, 2026-11-19): for every anchor of set v2 (`anchors_native_mass_yeast_v2.tsv`, 22 rows),
the literature native mass and oligomeric state, the anchor's assigned mass, and its residual in each yeast dataset's
calibration fit. Result: **set v3 = v2 minus PDC1 and EFT1 (20 anchors)**, `anchors_native_mass_yeast_v3.tsv`, preset
`default` for yeast (`calibration.anchor_set: {human: v2, yeast: v3}`; human unchanged at v2). v1 and v2 files are kept
unchanged (provenance; v1 = preset legacy_v7, v2 selectable by `overrides: {calibration: {anchor_set: v2}}`).

## Method

- Datasets (`datasets/update2026_yeast.yaml`): lab 2026 NaCl (60 fr), Pang et al. 2020 PXD019513 (BioSep s4000, 70 fr,
  one profile), Swamy et al. 2022 PXD028358 (Superose 6, 27 fr, 4 wild-type profiles; mean shown). Lab YPD reuses the NaCl
  calibration (not a separate fit).
- Residual = log10(native kDa) - fitted log10(kDa) at the anchor's apex fraction (negative = elutes larger than its native
  mass), shown **relative to the median residual of the dataset's v3 anchors** (the Theil-Sen fit is not median-centred),
  each fit with the v3 set (rows removed in v2/v3 are scored against the v3 fit). Apex fractions are from the saved smoothed
  profiles of the 2026-11-17 update run (`secms_lab2026/final/yeast_update/runs/yeast/profiles.tsv`; the v2 QC numbers are
  reproduced exactly from them); the QC table below is the toolkit re-run with v3.
- Expected systematic pattern: small monomers (< 50 kDa) elute ~1.5x larger and ~150 kDa enzymes ~1.5-2x smaller than a
  single log-linear fit (curvature of the column calibration), in all three datasets. Only deviations beyond that,
  consistent across datasets and with a documented biochemical reason, justify removal (v2 rule, Pedro 2026-11-17).

## Per-anchor table

| anchor | assigned (kDa) | literature state, native mass | source | NaCl | Pang | Swamy | verdict |
|---|---|---|---|---|---|---|---|
| PFY1 | 1 x 13.7 | monomer, 13.7 | PDB 1YPR (S.c. profilin, Eads et al. 1998) | n.d. | -0.07 | -0.13 | keep |
| CPR1 | 1 x 17.4 | monomer, 17.4 | PDB 1IST (S.c. Cpr1); UniProt P14832 | -0.21 | -0.17 | -0.17 | keep (small-protein curvature; human PPIA kept) |
| SOD1 | 2 x 15.9 = 31.7 | homodimer, 31.7 | PDB 1SDY (Djinovic et al. 1992, yeast Cu,Zn SOD dimer) | -0.27 | -0.23 | -0.16 | keep (curvature; apo-Sod1/Ccs1 heterodimer minor) |
| PGK1 | 1 x 44.7 | monomer, 44.7 | PDB 3PGK (yeast PGK, Watson et al. 1982) | -0.18 | -0.17 | -0.04 | keep |
| TPI1 | 2 x 26.8 = 53.6 | homodimer, 53.6 | PDB 1YPI (yeast TIM) | -0.10 | -0.01 | 0.00 | keep |
| FBA1 | 2 x 39.6 = 79.2 | homodimer (class II aldolase), 79.2 | class II FBA dimer, PDB 1ZEN (E. coli homologue); UniProt P14540 homodimer | -0.14 | 0.00 | -0.03 | keep |
| MET6 | 1 x 85.9 | monomer, 85.9 | PDB 4QQU (C. albicans Met6, monomer; Ubhi & Robertus 2015) | 0.00 | -0.13 | 0.00 | keep |
| ENO2 | 2 x 46.9 = 93.8 | homodimer, 93.8 (ENO1/ENO2 heterodimers same mass) | PDB 1EBH / 2ONE (yeast enolase dimer) | +0.14 | +0.20 | +0.21 | keep (paralog mix does not change mass) |
| GPM1 | 4 x 27.6 = 110.4 | homotetramer (dimer of dimers), 110 | PDB 1QHF (Crowhurst et al. 1999); Rigden et al. 1998 (4PGM) | 0.00 | +0.14 | +0.11 | keep |
| PGI1 | 2 x 61.3 = 122.6 | homodimer, 122.6 | UniProt P12709 (homodimer; PGI family dimer, e.g. PDB 1HOX rabbit) | +0.05 | n.d. | +0.16 | keep |
| TDH3 | 4 x 35.7 = 143.0 | homotetramer, 143 (TDH1/2/3 mixed tetramers same mass) | PDB 3PYM (S.c. GAPDH3 tetramer) | +0.01 | +0.05 | +0.16 | keep |
| ADH1 | 4 x 36.8 = 147.4 | homotetramer, 147 (gel filtration MW marker, 150 kDa) | PDB 4W6Z (Raj et al. 2014, Biochemistry) | +0.39 | +0.47 | +0.17 | keep, **flag**: elutes 1.5-3x smaller (tetramer dissociation is known without stabilising Zn/Ca2+); same direction as TKL1 at the same mass (curvature), not removed |
| TKL1 | 2 x 73.8 = 147.6 | homodimer, 147.6 | PDB 1TRK (Lindqvist et al. 1992) | +0.23 | +0.39 | +0.24 | keep (curvature at ~150 kDa) |
| CDC19 | 4 x 54.5 = 218.2 | homotetramer (FBP-activated), 218 | PDB 1A3W (Jurica et al. 1998) | +0.03 | +0.23 | +0.13 | keep, flag: reversible aggregation under starvation/stress (Saad et al. 2017) - not in log phase |
| GLN1 | 10 x 41.8 = 417.7 | homodecamer, 418 | PDB 3FKY (He et al. 2009 Proteins; decamer); filaments only in starvation (Petrovska et al. 2014) | +0.15 | +0.43 | +0.27 | keep |
| CDC48 | 6 x 92.0 = 552 | homohexamer, 552 (+ cofactors Ufd1-Npl4, Shp1, Ubx: up to ~0.7 MDa) | PDB 6OA9 (Twomey et al. 2019, S.c. Cdc48-Ufd1-Npl4) | -0.21 | +0.19 | -0.06 | keep (cofactor-bound fraction; human VCP kept) |
| PRE1 | 700 | 20S CP 28 subunits, ~700-730 | PDB 1RYP (Groll et al. 1997) | -0.21 | n.d. | -0.09 | keep, **flag**: NaCl apex (fr 30) co-elutes with Blm10 (fr 29; Blm10-CP ~0.95 MDa) and is separate from 26S (RPT1/RPN2 fr 21); in Swamy every >= 1 MDa assembly (CP, 26S, CCT) peaks at fr 8-9 = top of the resolved range |
| PUP1 | 700 | as PRE1 | PDB 1RYP | -0.21 | n.d. | +0.18 | keep, as PRE1; one Swamy profile (LSCSC) apex at fr 18 (free precursor / noise; +0.96 in that run, absorbed by Theil-Sen) |
| CCT2 | 960 | TRiC 16-mer, ~960 | PDB 4V94 (Dekker et al. 2011, yeast TRiC) | +0.14 | -0.19 | -0.06 | keep, flag: in Pang and Swamy the apex is at fr 8 (void / top of range: censored) |
| TCP1 | 960 | as CCT2 | PDB 4V94 | +0.14 | -0.19 | -0.06 | keep, as CCT2 |
| **PDC1** | 4 x 61.5 = 246.0 | homotetramer **in equilibrium with dimer**: pH-dependent, midpoint ~pH 7.5, shifted by ThDP / Mg2+ / phosphate | Koenig et al. 1998 (SAXS, S.c. PDC); PDB 1PVD (Dyda et al. 1993) | +0.35 | +0.33 | +0.28 | **removed**: condition-dependent assembly; elutes ~2x smaller than the tetramer (= dimer) in all three datasets (Swamy buffer pH 7.4) |
| **EFT1** | 1 x 93.3 | eEF2 monomer, 93.3; partly ribosome-bound | PDB 1N0U (Joergensen et al. 2003) | +0.52 | +0.32 | 0.00 | **removed**: paralog ambiguity (EFT1/EFT2 encode the identical protein, one UniProt entry P32324 with gene field 'EFT1; EFT2': the row never matched any dataset in v1/v2, so removal changes no fit); forced to P32324 it is inconsistent (NaCl apex at the last fraction) |
| (COF1, BMH1, SSA1, HSP82) | - | removed in v2 | v2 note | -0.41/-0.37/-0.36/-0.42 | -0.41/-0.28/+0.03/+0.52 | -0.35/-0.33/-0.40/-0.38 | removed in v2 (confirmed) |

n.d. = not detected / not scored in that dataset (>= 2 peptides). Columns NaCl / Pang / Swamy = relative residual, log10.

Summary: 22 v2 anchors reviewed; 20 kept (5 flagged: ADH1, CDC19, PRE1/PUP1, CCT2/TCP1), 2 removed (PDC1, EFT1), 0 fixed
(assigned masses agree with the literature states for every kept anchor). Not changed here but noted: the void of Swamy
(fr ~8) and Pang (fr ~7-10) is not set in the dataset configs (`void_fraction: null`), so the >= 1 MDa anchors are
censored at the top of the range there; setting it is a dataset property, not an anchor-set change.

## Calibration QC with v3 (toolkit re-run, preset default, `datasets/update2026_yeast.yaml`)

QC gate: n >= 8, Spearman rho <= -0.7, residual SD fold <= 2.5.

| dataset / run | v2 n | v2 rho | v2 resid x | v3 n | v3 rho | v3 resid x | v3 LOO x | QC |
|---|---|---|---|---|---|---|---|---|
| lab NaCl | 20 | -0.852 | 1.552 | 19 | -0.852 | 1.510 | 1.714 | pass |
| Pang 2020 | 18 | -0.828 | 1.701 | 17 | -0.816 | 1.682 | 2.024 | pass |
| Swamy HSC16L | 21 | -0.967 | 1.406 | 20 | -0.966 | 1.366 | 1.707 | pass |
| Swamy HSC8p15L | 20 | -0.968 | 1.488 | 19 | -0.962 | 1.481 | 2.045 | pass |
| Swamy HSCSC | 21 | -0.959 | 1.399 | 20 | -0.957 | 1.369 | 1.639 | pass |
| Swamy LSCSC | 21 | -0.829 | 1.897 | 20 | -0.826 | 1.908 | 2.065 | pass |

Lab YPD: NaCl calibration (calibration_from), unchanged rule.

## Sources

- Koenig S, Svergun D, Koch MHJ, Huebner G, Schellenberger A (1998) Eur Biophys J / Biochemistry: pH- and cofactor-dependent
  dimer-tetramer equilibrium of yeast pyruvate decarboxylase (X-ray solution scattering; midpoint ~pH 7.5).
  https://www.yeastgenome.org/reference/S000064324 ; https://jp4.journaldephysique.org/10.1051/jp4:1994953
- RCSB PDB entries cited in the table (https://www.rcsb.org/structure/<id>): 1YPR, 1IST, 1SDY, 3PGK, 1YPI, 1ZEN, 4QQU,
  1EBH, 1QHF, 4PGM, 3PYM, 4W6Z, 1TRK, 1A3W, 3FKY, 6OA9, 1RYP, 4V94, 1PVD, 1N0U.
- He YX et al. (2009) Proteins 76:249 (Gln1 decamer, 3FKY); Petrovska I et al. (2014) eLife 3:e02409 (Gln1 filaments in
  starvation); Saad S et al. (2017) Nat Cell Biol 19:1202 (Cdc19 reversible aggregation); Raj SB et al. (2014)
  Biochemistry 53:5791 (ADH1 tetramer, 4W6Z); Ubhi DK, Robertus JD (2015) J Mol Biol 427:901 (Met6, monomer).
- UniProt P32324 (EF2_YEAST: gene names EFT1; EFT2, identical products).
- Data: per-anchor apex fractions and residuals, `secms_lab2026/final/yeast_update/runs/yeast/secms_calibration.tsv` (v2) and
  the v3 run `wholecell_builds/yeast/_dev_unify/step3_secms_assembly/secms_v3_run/secms_calibration.tsv`.
