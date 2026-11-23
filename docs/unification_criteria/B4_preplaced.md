# B4 common pre-placed builders: acceptance criteria (frozen before any B4 build)

Step B4 of the unification plan (`docs/reference/history/unification_plan.md` §1 rows "Pre-placed: cytoskeleton / NPC /
ribosomes-polysomes-mRNA / cell-type extras", §5b B4): the pre-placed builders that both cells need (microtubules,
nuclear pore complexes, polysomes and mRNA, mitochondrial machines) become ONE core module each, with their numerics
ported out of the vendored yeast engine (`# origin: yeast ... @ e3b7448`), driven by config. Jurkat switches from the
frozen v7 layers to these builders; yeast keeps reproducing its pre-pass. Written 2026-11-17 before any B4 code or
build existed; the only numbers used are those of the reference builds, the frozen v7 layers, the geometry package
em_cell1_v8 and the EM ribosome table (measurements listed with their source). The md5 of this file goes to
PLANNING_NOTES (coordinator). Later findings go to a dated addendum file, never into this one.

Plan §5b B4 acceptance: "ribosome NN median ± 1 nm; polysome fraction within literature (Jurkat EM, yeast 0.82-0.85
target); MT dimer count ± 1 %; pore count = config; pre-pass self-overlaps (C7a pre-pass-pre-pass) ≤ reference".
Made concrete below.

## References

| Cell | Reference build for no_regression (in this order of preference) | Pre-pass reference |
|---|---|---|
| Jurkat | the B3 build (path in `jurkat/_dev_unify/b3/full_build.status`, its `validate/suite.json`); if it has no suite.json: the B6 build (`jurkat/_dev_unify/b6/full_build.status`); else v7repro2 (`jurkat/_dev_unify/a8p2/20261115-ec0b403-s20261115-v7repro2/suite.json`) with the B6 and B3 EXPECTED lists added | the frozen v7 layers (`pre_polysomes_v7`, `pre_cytoskeleton_v7`, `pre_npc_v3`, `pre_mrna_v1`) as packed in v7repro2 |
| Yeast | the yeast B6 build (`yeast/_dev_unify/b6/full_build.status`, its `validate/suite.json`); else s0 (`yeast/_dev_unify/a8p2/20261116-aabe7b5-s0-x5repro/suite.json`) | port s0 + PORT PATCH 4 (+ PATCH 5): `yeast/_dev_unify/amphiphysin_fix/ref_yeastenv_patched/prepass/serial/` (core env: `.../run_vcell`, bitwise except the mRNA gene labels fixed by PATCH 5) |

Which reference was used is written into the verdict. No replicate-seed builds are run (Pedro); seed noise is judged
with the existing X5 s0/s1/s2 spread (yeast) and the §4 tolerances.

## What B4 builds (scope, fixed now)

ONE core module per builder in `wholecell/preplaced/builders/`, plain Python + numba kernels, no `ctx.mod` calls. The
geometry primitives a builder needs (label / membrane distances of a point, sphere-tree containment and overlap grid,
quaternions) come from the cell's evaluator through the context (`ctx.K`): yeast = its geometry evaluator of the X5
layout (vendored `cell3d` / `pack3d`, so its pre-pass stays bitwise), Jurkat = the core evaluator `cellgeom_eval_v1`
and core `pack3d`. Cell differences are config values and data tables only.

| Builder (core) | Yeast (unchanged behaviour) | Jurkat (B4: from frozen layer to builder) |
|---|---|---|
| `core.npc_pores` | one NPC (7N9F scaffold, lib `NPC`) per package pore, 93 pores | one NPC (7R5K, lib `NPC`) per package pore; pore count from the rule J1 (new geometry package) |
| `core.microtubules` (+ `core.microtubule_pools`) | track source `procedural`: nucleation from the SPB (MTOC) | track source `traced`: the EM polylines of the package layer `mt_traces` (variant b, 1,360 centerlines), the same 13_3 lattice (dimer 8.2 nm, phi -27.69 deg, dz 0.946 nm) as v7 |
| `core.polysomes` + `core.mrna_threads` | gene-specific transcripts, walk from random cytosolic starts, OM / cortical-ER tethers, collided disomes | the same walk, seeded by the EM point priors (`ribosome_priors` layer = `ribo_seg` cell-1 instances): every retained EM ribosome is placed at its EM position (± a seeded nudge) and starts or continues a chain; EM ribosomes of the ER class start ER-tethered chains; the remaining polysomal ribosomes extend the EM groups; free 80S are packed. mRNA threads through the chains (5 nt per bead, the v7 mRNA representation) for the mRNA set of the frozen transcript table (genes and lengths of `pre_mrna_v1`) |
| `core.atp_synthase_dimers` | rows on the rims of the crista prims | rows on the rims of the lamellar cristae (section `cristae`); dimer library entry built from the inventory monomer CPX-6151 (8H9T) with the yeast dimer geometry (monomer normals 84 deg apart, membrane centres on a circle of radius 8.8 nm; Davies et al. 2012 PNAS 109:13602) |
| `core.mtdna_nucleoids` | sites of the geometry, 0.45 genomes per fL | nucleoid sites sampled in the matrix (seeded, ≥ 500 nm apart, as the yeast fallback); genomes per cell = config (ASSUMPTION 1,000 mtDNA per cell, "hundreds to thousands", O'Hara et al. 2019 Genome Res 29:1878; 1.4 genomes per nucleoid, Kukat et al. 2011 PNAS 108:13534) |
| `core.npc_receptors`, `core.fg_mesh`, `core.mrnp_ends` | ported, unchanged | not enabled in B4 (Jurkat keeps the packer shells `sh_NPC`, `sh_NPC_fg`; PABPC1 / EIF4E packed free). Switching them on is a later decision |

Jurkat layers that stay frozen in B4: chromatin rep5, lamina (frozen copy), nascent RNA beads, and the non-MT part of
`pre_cytoskeleton_v7` (centrioles = kind 0, cortical F-actin = kind 2; selected from the frozen file by kind, not
copied). `config/cells/jurkat_v7_port.yaml` keeps every frozen v7 layer (the reproduction preset).

## Part J: Jurkat (measured pre-pack on the pre-placed stage, and again on the build)

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| J1 | **Pore count = config, from the literature density checked against the Nup anchor** (DECISION 2026-11-16, Pedro) | density rho = **11 NPC/µm²** (HeLa interphase steady state by EM, Otsuka et al. 2016 eLife 5:e19071; pore frequency nearly constant across cultured cell lines, Maul & Deaven 1977 J Cell Biol 73:748) × NE midsurface area of the geometry (em_cell1_v8 mesh: **245.1 µm²**) = **2,696 pores**. Anchor check: the nucleoporin anchor of the inventory (scaffold-Nup families vs 7R5K chains per pore; `docs/reference/history/copy_numbers_report.md` §4: median family ratio 2.24 at 1,000 pores, i.e. Nups for **2,240 pores**, recomputed on the build's inventory) must lie within a factor 2 of the count (absolute MS copy numbers are good to ~2-fold): 0.5 ≤ N_anchor / N ≤ 2. Here 2,240 / 2,696 = 0.83: passes. If the check failed, N = the nearer bound (2 N_anchor or N_anchor / 2). | `geometry.npc.n_pores` in config = 2,696; package pores = 2,696; placed NPC instances = 2,696; inventory NPC copies (assembly cap and rule v4 d `npc_copies`) = 2,696; anchor ratio in [0.5, 2] |
| J1b | Extra pores | a new geometry package (step from em_cell1_v8, never edited): the 1,000 pores of v8 (300 EM + 700 random) unchanged and first; 1,696 more at area-weighted random centres of NE-mesh triangles (the jgeom_build fill rule: numpy default_rng, seed recorded in the manifest), ≥ 150 nm from every other pore, axis = the triangle normal (outward), source code 0. Every other content file identical (sha256) to v8 | validator `pass`; `pores` check without waiver; min pore separation ≥ 150 nm; all package files except `pores.npz` and the manifest identical to v8 |
| J2 | **MT dimer count ± 1 %** | tubulin dimers placed on the traced MTs (lattice instances, excluding centrioles) | 2,325,171 ± 1 % = **2,301,919 - 2,348,423** (v7: kind 1 of `pre_cytoskeleton_v7`). Also reported: centerlines used (v7 1,360), total length (v7 1,456.2 µm), dimers dropped by clashes per cause |
| J3 | **Polysome fraction (Jurkat EM)** | (a) model bookkeeping: ribosomes in pre-placed chains of ≥ 2 / all cytosolic + ER RIBO80S copies of the inventory (pre-placed + packed). Target = config `polysome_fraction` = **0.70**, set from this cell's EM: the share of EM ribosomes (ribo_seg, cell 1, non-nuclear, pairs < 20 nm merged as duplicate detections: 682,785) in single-linkage groups of ≥ 2 at 34 nm (the v7 polysome-nucleus rule) = **0.697** (measured 2026-11-17, `jurkat/_dev_unify/b4/em_polysome_stats.json.txt`; 0.738 without the duplicate merge). (b) EM-comparable: the same linkage statistic on the model ribosome centres (pre-placed + packed RIBO80S) thinned at random (stream `validate.ribosome_nn`) to 682,785 | (a) 0.68-0.72; (b) within ± 0.05 of 0.697. v7 values reported (a: 0.945) |
| J4 | **Ribosome NN median ± 1 nm** | suite C12, `validate/ribosome_nn.json`: median NN distance of the model ribosomes thinned to the EM count (730,960). In B4 the free 80S are packed, so the row "model pre-pass + packed RIBO80S, thinned to EM count" is the one compared (v7: all ribosomes were pre-placed, its "thinned" row 28.24 nm) | **27.24-29.24 nm** (v7repro2 28.24; EM 27.05 reported). Pre-build: the pre-placed ribosomes thinned to the EM count, reported |
| J5 | **Pre-pass self-overlaps ≤ reference** | (a) suite C7a `prepass_prepass` pairs per 10⁶ pre-pass instances on the 300 sampled tiles (`validate/placement.json`; overlap > 0.5 nm of the collision sphere trees) vs the reference build (v7repro2: 838,914.7). (b) the same pair count split by module pair (each instance labelled by its pre-placed module), for every pair class with at least one B4 module (npc, microtubules, polysomes, atp_synthase, mtdna is beads: not counted), per 10⁶ instances of the B4 modules, vs the same classes in the reference (v7 polysomes, cytoskeleton MT part, npc). Pre-build: (b) on the reference build's tile sample (v7repro2 tiles, stream `validate.placement.tiles`) for the new and the v7 layers with the same code | (a) ≤ reference + max(10 % rel., 1 per 10⁶); (b) each class ≤ its reference value + max(10 % rel., 1 per 10⁶), or the class is new (ATP synthase) and reported |
| J6 | Mito machines (new layers) | ATP synthase dimers: count, share of the IM monomer pool used (cap: half, as yeast), every dimer's TM residues within 2.5 nm of the IM midsurface, no non-TM residue in a membrane beyond the yeast tolerance (2 %); mtDNA: nucleoids, genomes, Mb placed, all beads in `mito_matrix` | dimers > 0 and ≤ pool / 2; geometry tests 100 %; genomes = config; ≥ 95 % of the DNA placed |
| J7 | Consumption bookkeeping | module jsons `consumed` (RIBO80S, TUBA1B-TUBB4B, NPC, CPX-6151) = placed counts (two monomers per ATP dimer); pack counts = pack table − consumed; inventory rule v7 MT3b reads the new MT json | identical numbers |
| J8 | Geometry and frozen layers | `depends_on_geometry` of every frozen module = the new package; the B6 index drop `ne_lumen_b6` still applies to frozen modules; chromatin nucleosomes inside the INM of the new package ≥ 99.9 % (B6 C1; the INM is unchanged, only pores are added) | 0 mismatches; ≥ 99.9 % |
| J9 | Reference stays selectable | `config/cells/jurkat_v7_port.yaml` uses all frozen v7 layers and em_cell1_v7; its dry run lists the same preplaced jobs as at main 4b36303 | yes |

## Part Y: yeast unchanged

| # | Criterion | Pass |
|---|---|---|
| Y1 | Serial pre-pass on s0's inputs in the core env with the B4 code (`scripts/preplaced_a6_check.py` run, 3 threads) vs the PATCH-4 reference | every array of builders 1-28 identical to `ref_yeastenv_patched` (the core-env difference of PATCH 5 included: `mrna_gene` / `poly_chain_gene` now equal); pools after the pre-pass identical; log lines identical (timing removed). Any difference = NOT ACCEPTED unless it is float-summation order at the 1e-6 nm level (env tolerance, listed) |
| Y2 | Pore count = config | 93 = package pores = placed NPCs (unchanged) |
| Y3 | MT dimer count | identical to the reference (`n_tub`; layer arrays identical) |
| Y4 | Polysome fraction | suite C12 `polysome_in_chains_fraction` reported against the plan's 0.82-0.85 target (`POLY_FRACTION` 0.835: Arava et al. 2003 PNAS 100:3889, 85 ± 5 %; Metzl-Raz et al. 2017 eLife 6:e28034, 0.82). The reference has **0.784** (chains end short when the walk runs out of space); B4 does not change yeast, so this is reported, not a B4 failure. Changing it is a yeast model decision (G. Magnusson) |
| Y5 | Ribosome NN | C12 `nn_median_model` = reference (s0: 28.39 nm; CZII cryo-ET 29.38 reported) within the §4 tolerance (± 0.5 nm) |
| Y6 | Self-overlaps | C7a `prepass_prepass` per 10⁶ ≤ reference (s0: 642,511.1) + max(10 % rel., 1 per 10⁶) |

## Part B: build-level acceptance (one build per cell from the B4 commit, `scripts/build_isolated.sh <cell> b4 24`)

| # | Criterion | Pass |
|---|---|---|
| B1 | Build complete; pre-placed stage ran the B4 builders (Jurkat: `preplaced/serial/serial_manifest.json` lists them; yeast: same builders as A6) | yes |
| B2 | Common suite, `no_regression.tsv` vs the reference suite with the §4 tolerances | no regression outside the EXPECTED list below |
| B3 | J1-J8 recomputed on the build (`scripts/b4_preplaced_check.py verdict`); Y2-Y6 on the yeast build | the bounds above |

**EXPECTED changes (Jurkat; caused by B4, so not regressions within these bounds).**

- NE: 2.7× more NPC instances; C11 `density.perinuclear` (report) and perinuclear volume decrease by ≈ 1,696 × 2.54e-4 µm³ ≈ 0.43 µm³ (pore openings, B6 J4 arithmetic); free Nup pools decrease (more Nups in pores); C1/C2/C3 rows of Nups and NPC-shell species may move: placed metrics ≤ 0.01 (standard), loc_compare rows of Nups reported.
- Ribosomes: up to 30 % of RIBO80S now packed free (monosomes not at an EM position; v7: 0 % packed): C10 placed fractions ≥ reference − 0.003; C11 `density.cytosol` ± 3 mg/ml unless placed mass changed (then mass-consistent). C12 as J4.
- C7: (a) packed-packed / packed-pre-pass rates standard; pre-pass-pre-pass as J5. C7 `prepass_mrna_per_1e6` (report) changes with the rebuilt mRNA.
- MTs: same tracks; dimer phase and clash drops differ (J2).
- New matrix / IM content: ATP synthase dimers (IM) and mtDNA (matrix): C11 `density.matrix` ± 2 mg/ml or mass-consistent.
- C13 chromatin metrics unchanged (frozen file). Tile RNG trajectories diverge everywhere: pack metrics differ at seed-noise level (standard tolerances).
- `compare_v7` (still in the DAG) reports deviations in the pre-pass; informative only.

**EXPECTED changes (yeast).** None beyond the environment: Y1 identical arrays ⇒ the yeast build is the B6 yeast build's
content apart from code-path provenance; suite metrics equal to the reference within the §4 tolerances (tile packing
reads identical inputs).

**Verdict.** `<build>/validate/compare_B4/verdict.json` per cell (`scripts/b4_preplaced_check.py verdict`): `ACCEPTED`
if parts J (Jurkat) or Y (yeast) and B1-B3 pass and B2 has no regression outside the expected list; otherwise
`NOT ACCEPTED` with the failing items. Pre-build measurements: `{results_root}/<cell>/_dev_unify/b4/`.
