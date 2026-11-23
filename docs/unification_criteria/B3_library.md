# B3 Jurkat library rebuild: acceptance criteria (frozen before any B3 build)

Step B3 of the unification plan (`docs/reference/history/unification_plan.md` §1 "Structure library", §2b P8, §5b B3):
the Jurkat structure library is rebuilt by the core library stage (`library.mode: build`) with the Jurkat builder steps
ported from `wholecell/library/legacy/`, into a new canonical library per build (`<build>/library/lib/`, never edited in
place). The frozen v7 library (`lib_structures_v7`, `{data_root}/jurkat/yc/pipeline/lib`, manifest md5 `e4ba6bb3…`)
is never written; the reference preset `config/cells/jurkat_v7_port.yaml` keeps it (`mode: frozen`). Written 2026-11-17
before any B3 code or build existed; the only numbers used are those of the frozen v7 library and its legacy builder
logs. The md5 of this file goes to PLANNING_NOTES (coordinator). Later findings go to a dated addendum, never into this
file.

Plan §5b B3 acceptance: "per-species volume within 2 %, R_h within 2 % for ≥ 99 % of species; TM orientation:
hydrophobic belt centred within 0.5 nm; C9 ≤ reference". Made concrete below.

## What the rebuild is (scope, fixed now)

The v7 library was written in place by five legacy scripts over 2026-11-08/09 (PLANNING_NOTES §7, §7b, v5 R3):
`jlib_build` (assemblies + AFDB monomers) → `jlib_rebuild_tm` (TM assemblies again after the SIFTS fix, with
`tm_graft`) → `lib_refine_st` (adaptive sphere trees, ×3: soluble, membrane, tm_route ids) → `jlib_reorient_tm`
(1-2-pass TM monomers, incl. the "matrixfix" re-run) → `jlib_tm_route` (129 TM-routed ingredients). B3 runs the final
code of each script once, in a fixed order, as one core step each: `start` (frozen entries) → `jurkat_assemblies`
(option `tm_graft`) → `jurkat_monomers` → `reorient_tm` → `tm_route` → `adaptive_trees` → `seal` → `pack_shapes`.
Organism data are declared inputs (species table TSV in the repo; UniProt topology / REST cache, SIFTS, AFDB complex
models, AFDB human proteome tar and its fragment index, the RCSB assembly cache read offline: a missing file is an
error, nothing is downloaded).

- **Species table.** One row per ingredient of the legacy table the library was built from (16,279 rows; 16,242 have
  a v7 npz written by the scripts: 16,226 by `jlib_build`, 16 AFDB fragments by `tm_route`; 37 had no model): ing_id, kind, genes, accessions,
  stoich, transmembrane flag, signal fraction and mass of the legacy table the library was built from
  (`ingredients_v0_pre_v4.tsv`), and **the structure source recorded in the v7 npz (`src`)**. The library src differs
  from the inventory `struct_src` for 323 assemblies (legacy re-derivation of the complexes table after the library
  build); B3 keeps the structure that v7 packed. Moving the library to the inventory's struct_src is a separate,
  later decision (it would change ~2 % of species beyond these criteria).
- **Frozen start entries** (copied, md5-checked against the v7 manifest): the 2 pre-pass templates not built by these
  scripts (`F_ACTB`, `MT_TUBA1B-TUBB4B`) and their 2 json files; plus any entry the scripted steps cannot rebuild offline,
  each listed with its md5 and reason in the start table (expected: none beyond these 4 files).
- **Packing shapes.** The Jurkat pack keeps using the library's own adaptive sphere trees (`st_c`/`st_r`; `pack.shapes.
  source: library`, the v7 policy). The stage also writes the P8 packing shapes (`library/pack_shapes/`, 15 residues per
  sphere) for B1; they are not used by the pack in B3.

## Part L: library vs the frozen v7 library (measured on the built `library/lib/`, before the pack)

Measured by `scripts/library_b3_check.py` (written with the implementation), per species = per npz of the v7 manifest.

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| L1 | Species set | file names of `library/lib/` vs the v7 manifest (16,244 npz + 2 json) | identical sets; manifest `origin` column names the step that wrote each file; start entries = the declared list |
| L2 | Volume | V = union volume of the collision sphere tree (st_c, st_r): voxel count on a grid of spacing 0.25 nm aligned to the library origin (node (i,j,k)·h), × h³; same code for both libraries | \|V_new − V_v7\| / V_v7 ≤ 2 % for **≥ 99 %** of species |
| L3 | R_h | Kirkwood hydrodynamic radius of the residue beads, 1/R_h = ⟨1/r_ij⟩ over bead pairs i ≠ j; species with > 3,000 beads use a fixed subset of 3,000 beads (numpy default_rng(0).choice without replacement, sorted) | \|R_h,new − R_h,v7\| / R_h,v7 ≤ 2 % for **≥ 99 %** of species. Also reported: both L2 and L3 jointly; bead count identical; beads bitwise identical |
| L4 | TM set | species with the `membrane` key (v7: 3,530 TM + 1 NE) and their value; `tm_route` method of the 122 routed species | same species, same values; tm_route method identical for ≥ 95 % of the routed species (the rest listed) |
| L5 | Belt centre | z_belt = median z of the membrane residues: beads flagged in `tm_beads` where the key exists, else beads with \|z\| < 2.0 nm (the slab); library frame (membrane mid-plane z = 0) | \|z_belt\| ≤ 0.5 nm for ≥ 99 % of TM species, and the count failing is ≤ the v7 count + 1 % of TM species |
| L6 | Orientation vs v7 | for TM species with the same bead count: Kabsch superposition of the v7 beads onto the new beads; tilt = angle between the v7 normal (R·ẑ) and ẑ; shift = z of the v7 mid-plane origin in the new frame (t_z) | tilt ≤ 10° **and** \|shift\| ≤ 0.5 nm for ≥ 99 % of those species; no flipped species (tilt > 90°) without a listed cause |
| L7 | Sphere-tree sizes | spheres per species; adaptive rule (tree volume ≤ 2.5 × 1.21 nm³/kDa × mass, else the finest per = 5), species > 20,000 beads skipped as in v7 | sphere count identical for ≥ 95 % of species; total spheres within 2 %; share of non-skipped species over 2.5 × ≤ the v7 share + 0.5 % |
| L8 | Determinism | `seal` manifest md5s from two library builds of the same commit with different `nproc` (3 vs 24: the cheap check and the full build) | identical (the builders are seeded; KMeans random_state 0) |

Volumes and R_h use only the npz content (no inventory). A species failing L2/L3 is listed with its v7 builder log
note; > 1 % failing is NOT ACCEPTED unless a dated DECISION explains a systematic cause.

## Part Y: yeast unchanged

B3 does not change yeast: its library is already scripted (`library.mode: build`, A3/A9) and the yeast steps
(`complexes`, `af_models`, `pdb_oligomers`, `pdb_rna`, `seal`, `pack_shapes`) are not modified.

| # | Criterion | Pass |
|---|---|---|
| Y1 | `config/cells/yeast.yaml` `library:` section | byte-identical to main a4db4f1 except the `options` comment line (the Jurkat options become real steps) |
| Y2 | yeast library code path | `tests/test_library.py::test_build_mode_synthetic` passes unchanged; dry-run DAG of `yeast` lists the same jobs as at a4db4f1 |
| Y3 | yeast library output | the yeast library of the next yeast build from a B3-or-later commit has the same `library_manifest.tsv` (species, md5) as the yeast B6 build (checked when that build exists; no extra yeast build for B3) |

## Part B: build-level acceptance (one Jurkat build from the B3 commit, `scripts/build_isolated.sh jurkat b3 24`)

Reference: the **B6 Jurkat build** (`jurkat/20261117-04dec6f-s20261115-b6`, path from
`jurkat/_dev_unify/b6/full_build.status`) and its `validate/suite.json`, once its B6 verdict exists. If it does not
exist (or failed), the reference is v7repro2 with the B6 EXPECTED list of `B6_geometry.md` added to the B3 one.

| # | Criterion | Pass |
|---|---|---|
| B1 | Build complete; library stage mode `build`; `library/library_meta.json` lists the B3 steps; L1-L7 recomputed on the build's library | yes; same numbers as the pre-build check (L8) |
| B2 | Common suite, `no_regression.tsv` vs the reference suite with the §4 tolerances | no regression outside the EXPECTED list below |
| B3 | **C9 ≤ reference**: sign errors per 10⁶ membrane instances and median tilt | each ≤ the reference value + its §4 tolerance (rate: max(10 % rel., 1 per 10⁶)) |
| B4 | Library provenance | every library file has an `origin`; start entries md5 = v7 manifest; inputs read (AFDB tar, SIFTS, UniProt cache, fragment index, species table) recorded with md5 / size in `library_meta.json` |

**EXPECTED changes (Jurkat; caused by the rebuild, so not regressions within these bounds).** Library entries that are
not bitwise identical change the collision trees of their species, so tile RNG trajectories diverge: pack metrics
differ at the seed-noise level everywhere.

- C10 `placed_mass_fraction`, `placed_copies_fraction` ≥ reference − 0.003 (standard); C11 densities ± 2 mg/ml.
- C7 / C8: standard tolerances. C9 as B3 above.
- C1 / C2 / C3 placed metrics: ≤ 0.01 (standard).
- C12 / C14 (crowding, diffusion proxies): standard tolerances; library volume changes ≤ 2 % per species (L2) bound
  them.
- C13 frozen chromatin metrics unchanged (frozen file; pre-placed modules keep their frozen positions).
- The v7 reproduction comparison `compare_v7` (still in the DAG) reports deviations; informative only.
- Build wall time: + the library stage (expected ≤ 1.5 h on 24 cores; reported).

**Verdict.** `<build>/validate/compare_B3/verdict.json` (`scripts/library_b3_check.py verdict`): `ACCEPTED` if L1-L7,
Y1-Y2 and B1-B4 pass and B2 has no regression outside the expected list; otherwise `NOT ACCEPTED` with the failing
items. Pre-build measurements: `{results_root}/jurkat/_dev_unify/b3/`.
