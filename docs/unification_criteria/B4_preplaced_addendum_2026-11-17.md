# B4 pre-placed builders: addendum 2026-11-17 (per-cell polysome linker; J5b per pair class)

Addendum to `docs/unification_criteria/B4_preplaced.md` (frozen, md5 e715c2f90e17653369b7d5ac2c4d3d96, not edited).
**Written after the B4 pre-build check** (commit a182b2c; results `jurkat/_dev_unify/b4/measure_prebuild.json`), i.e.
post hoc: the frozen document predicted J3b 0.55 vs EM 0.697 and J4 32.4 nm vs the window 27.24-29.24 nm (with a random
complement for the packed free 80S 0.43 / 35.5 nm), and J5b per 10⁶ B4 instances moved for MT-MT when fewer
ribosomes were pre-placed. Decisions by the main session (2026-11-17, Pedro informed). Both the frozen values and the
addendum values are always reported (`scripts/b4_preplaced_check.py`: `measure.json` J5 `classes` = frozen,
`addendum` = this file; `verdict.json`: `items` / `verdict_frozen` = frozen, `items_addendum` / `verdict` = with
this addendum). Nothing in the frozen bounds of J3 / J4 changes.

## 1. Polysome linker: one common walk, a per-cell config parameter measured from EM

**What changes.** `core.polysomes` keeps one walk for both cells. The free-nucleotide linker between neighbouring
ribosomes (mean of the exponential draw per link) becomes the config parameter `params.linker_nt_mean`:
absent / null = per chain from the transcript (cds / ribosomes − footprint 28 nt; yeast, its literature gene table:
**unchanged**, the yeast pre-pass stays bitwise identical to the PATCH-4 reference, see checks); a number = the same
mean for every chain (Jurkat). Only the 3D walk reads it; ribosome loads per mRNA and mRNA lengths are unchanged.

**Why.** The B4 walk with the Jurkat transcript table (219 nt per ribosome, ~191 free nt per link) places
consecutive ribosomes 31.9 nm apart (median centre-to-centre), v7 had 23-29 nm. The spacing of ribosomes inside
polysomes is a structural property of the cell that EM measures; it should be a measured per-cell value, not the
consequence of a transcript-table density. **The parameter is measured per cell from EM data, not tuned to pass the
gate**: it is chosen on a calibration subset declared here, and the frozen gates J3b / J4 are then evaluated, unchanged,
on the whole cell.

**Data.** `val_ribo_em` (OpenOrganelle jrc_jurkat-1 `ribo_seg` instances, md5 ec43d48f…), cell 1 (`cell2 == 1`), class 2
(nuclear) excluded, voxel 4 × 4 × 3.44 nm, model frame (minus the geometry centre): 730,960 centres (J4 EM set) and
682,785 after merging pairs < 20 nm (J3 set).

**Calibration subset (declared before the scan).** 1-µm blocks of the model frame, i = floor(x / 1000 nm) (same j, k);
a point is in the subset iff (i + 2j + 3k) mod 5 == 0 (one block in five, spread over the cell; 140,574 of the 730,960
EM ribosomes). Statistics are taken over the points inside the subset (neighbours searched among all points). The
subset is declared separately from the evaluation; the evaluation (J3b / J4) is the whole cell as frozen; the
complement (4 blocks in 5, never used for the choice) is reported as well.

**Method** (`scripts/b4_polysome_calibration.py`). For each scanned value a full Jurkat walk (builders npc,
microtubules, polysomes of `config/cells/jurkat.yaml`, same seeds and EM priors) is run; the model ribosome set =
pre-placed + a uniform random cytosolic complement (pool − pre-placed, seed 20261117) standing in for the packed free
80S; it is thinned exactly as J3b / J4 (`b4_preplaced_check.ribosome_stats` streams). Score =
|NN_median − NN_EM| / 1 nm + |linked − linked_EM| / 0.05 (the frozen J4 / J3b tolerances as units) on the subset;
the value with the lowest score **among those that keep the polysome fraction J3a in its frozen window 0.68-0.72**
(a linker that makes the walk end chains short changes the EM-set polysome fraction) is taken.

Subset EM targets: NN median 27.20 nm, linked share 0.691; EM ribosomes in 34-nm groups: median NN inside groups
25.84 nm.

| linker_nt_mean | J3a (fraction) | consecutive centres, median nm | subset NN nm (model + complement) | subset linked | score |
|---|---|---|---|---|---|
| 5 | 0.285 | 28.94 | 37.67 | 0.399 | 16.30 |
| 20 | 0.476 | 29.15 | 35.00 | 0.456 | 12.49 |
| 60 | 0.651 | 30.20 | 34.31 | 0.473 | 11.46 (lowest, J3a fails) |
| 90 | 0.672 | 30.75 | 34.65 | 0.460 | 12.06 |
| 120 | 0.680 (0.6797) | 31.21 | 34.82 | 0.449 | 12.46 |
| **150** | **0.684** | 31.58 | 35.38 | 0.437 | **13.26 (chosen)** |
| per transcript (B4 default, ~191) | 0.684 | 31.90 | 35.70 | 0.426 | 13.80 |

**Value: `linker_nt_mean: 150.0`** in `config/cells/jurkat.yaml` (provenance comment there; calibration results
`jurkat/_dev_unify/b4/addendum_20261117/calib/calibration.json`; date 2026-11-17).

**Finding.** The linker moves the EM-comparable statistics only a little: the centre spacing of consecutive ribosomes
has a floor near 27 nm (the collided-disome pose, 27.0 nm; the walk's coarse 80S trees), shorter linkers make the walk
end chains early (fewer pre-placed ribosomes, J3a falls), and the thinned statistics are dominated by the thinning of a
model with 2.9× more ribosomes than EM detections. No scanned value reaches the EM targets on the subset (best NN
34.3 nm vs 27.2; best linked 0.47 vs 0.69). The gap is not a linker-length effect; what would close it (denser
packing of chains around the EM groups as in v7, finer walk trees, or the thinning model of J3b / J4) is a separate
decision.

**Evaluation on the whole Jurkat pre-placed stage** (`jurkat/_dev_unify/b4/addendum_20261117/stage/`, frozen J3b / J4,
not the subset): serial pass of all six builders (804 s,
47 GB), polysome arrays identical to the calibration run at 150 nt; 1,615,692 pre-placed RIBO80S, 490,136 left to the
packer; the packed free 80S are again represented by the random cytosolic complement (`measure_stage.json`).

| | frozen bound | B4 pre-build (default linker) | this stage (150 nt) | pass |
|---|---|---|---|---|
| J3a polysome fraction | 0.68-0.72 | 0.684 | 0.684 | yes |
| J3b linked share, pre-placed + complement, thinned | 0.697 ± 0.05 | 0.431 | **0.441** | **no** |
| J3b pre-placed only (reported) | | 0.552 | 0.568 | (no) |
| J4 NN median, pre-placed + complement, thinned | 27.24-29.24 nm | 35.54 | **35.22** | **no** |
| J4 pre-placed only (reported) | | 32.40 | 32.13 | (no) |

Evaluation-disjoint complement of the subset (4 blocks in 5): 35.19 nm / 0.442 (EM there 27.05 / 0.698). **J3b and
J4 fail on the full stage with the measured linker, as predicted by the pre-build check; the measured linker improves
them by 0.3 nm / 0.01 only.**

## 2. J5b: per pair class, with absolute counts

**What changes.** Frozen J5b normalises every pair class by the instances of *all* B4 modules on the sampled tiles, so
a class moves when another B4 module changes its instance count (pre-build: MT-MT +11 % per 10⁶ B4 instances with
unchanged MT-MT pairs, because fewer ribosomes are pre-placed). New rule:

- every pair class with at least one B4 module (npc, microtubules, polysomes, atp_synthase) is reported with its
  **absolute pair count** (sampled tiles, same tiles, trees and overlap definition as frozen J5) and its **rate =
  pairs per 10⁶ instances of the class's B4 module(s)** on those tiles (MT-MT per 10⁶ MT dimers; chromatin-NPC per 10⁶
  NPCs; NPC-polysome per 10⁶ NPCs + ribosomes);
- **gate: each class rate ≤ its reference (v7 layers) class rate + max(10 % rel., 1 per 10⁶)** (the frozen §4 rate
  rule); a class absent from the reference counts as rate 0 there; a class with a module absent from the reference
  (ATP synthase) is new and reported (no gate);
- the absolute total of pre-pass-pre-pass pairs and the frozen J5a (all pairs per 10⁶ pre-pass instances ≤ reference +
  max(10 %, 1)) are reported and J5a stays a gate unchanged.

The frozen J5b values are still computed and reported next to these.

**Pre-placed stage result** (v7 layers vs the stage, v7repro2 tile sample): `measure_stage.json` J5 (300 v7repro2 tiles; the lamina is the
frozen v7 lamina here, the build uses the rebuilt lattice of mode build). Absolute pre-pass-pre-pass pairs 2,131,007
(v7 layers 2,200,254); all pairs per 10⁶ pre-pass instances 822,622 vs 842,927.

| class | pairs stage / v7 | B4 instances stage / v7 | per 10⁶ stage | per 10⁶ v7 | limit | addendum | frozen J5b |
|---|---|---|---|---|---|---|---|
| microtubules-microtubules | 201,450 / 202,219 | 96,293 / 96,635 MT | 2,092,052 | 2,092,606 | 2,301,867 | pass | fail (1,213,920 vs 1,088,961 per 10⁶ B4) |
| polysomes-polysomes | 2,029 / 58,419 | 69,405 / 89,024 | 29,234 | 656,216 | 721,838 | pass | pass |
| chromatin-npc | 33 / 11 | 110 / 40 NPC | 300,000 | 275,000 | 302,500 | pass | fail (198.9 vs 59.2) |
| cytoskeleton_rest-npc | 28 / 0 | 110 / 40 | 254,545 | 0 | 1 | **fail** | pass (absent class treated as new) |
| lamina-npc (frozen lamina) | 2,511 / 0 | 110 / 40 | 22,827,273 | 0 | 1 | **fail** | pass (absent class treated as new) |
| npc-polysomes | 3 / 47 | 69,515 / 89,064 | 43 | 528 | 580 | pass | pass |
| cytoskeleton_rest-polysomes | 19 / 14,422 | 69,405 / 89,024 | 274 | 162,001 | 178,201 | pass | pass |
| cytoskeleton_rest-microtubules | 15 / 16 | 96,293 / 96,635 | 156 | 166 | 182 | pass | pass |
| microtubules-polysomes | 0 / 199 | 165,698 / 185,659 | 0 | 1,072 | 1,179 | pass | (not listed) |
| lamina-polysomes | 0 / 2 | 69,405 / 89,024 | 0 | 22 | 25 | pass | (not listed) |

J5b (addendum) on the stage: **fail** in two NPC classes: the 1,696 added pores overlap frozen centriole / cortical
actin instances (28 pairs; the NPC builder does not test the frozen obstacles) and the frozen v7 lamina (2,511 pairs;
the reason the lamina is in mode build; to be judged on the build with the rebuilt lamina). The frozen J5b fails in
MT-MT (normalisation artefact) and chromatin-NPC (2.7× more NPCs); both pass per class. J5a: 822,622 ≤ 842,927 + 10 %
(stage vs v7 layers, same tiles).

## Checks of this addendum

- Yeast (Y1): serial pre-pass of `config/cells/yeast.yaml` with this code (core env, 3 threads, 2,700 s;
  `yeast/_dev_unify/b4/run_vcell_b4_addenda/`): builders 1-18 identical to s0, 19-28 (polysomes, mrna included)
  identical to `ref_yeastenv_patched`, pre-pass log lines identical. Yeast does not set `linker_nt_mean`.
- Jurkat: calibration runs `jurkat/_dev_unify/b4/addendum_20261117/calib/` (7 walks, 7-9 min each, 1 thread);
  stage `addendum_20261117/stage/` (+ `measure_stage.json`, `run_stage.sh`).
- `tests/test_preplaced_b4.py` (new: linker config, J5b class rule), `tests/test_preplaced_serial.py` pass; snakemake
  dry runs jurkat (36 jobs), jurkat_v7_port (35), yeast (40), testcell (27) OK.
