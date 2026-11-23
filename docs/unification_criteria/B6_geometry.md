# B6 geometry fixes: acceptance criteria (frozen before any B6 build)

Step B6 of the unification plan (`docs/reference/history/unification_plan.md` §5b): Jurkat `ne_mid` fix in a new
geometry version, yeast per-membrane half thickness. Written 2026-11-17 before any B6 package or build existed. The only
B6-era numbers used here are analytic expectations and the already published v7 / X5 values. The md5 of this file goes
to PLANNING_NOTES (coordinator). Later findings go to a dated addendum file, never into this one.

References: Jurkat `jurkat/20261115-ec0b403-s20261115-v7repro2` (package `jurkat/em_cell1_v7`, suite
`jurkat/_dev_unify/a8p2/.../suite.json`); yeast port s0 `yeast/20261116-aabe7b5-s0-x5repro` (+ PORT PATCH 4; package
`yeast/analytic_x5_s0_f5f3a5c723b2`, suite `yeast/_dev_unify/a8p2/.../suite.json`). Monte Carlo (MC) volumes below use
2,000,000 points, seed 0, stream `geometry.mc_volumes` over the cell-field box. That is the validator's own measurement,
the same as in v7 `validation.json`. CI = the validator's 95 % half-width.

## Part J: Jurkat geometry (package `jurkat/em_cell1_v8`)

The v7 defect (v7 `validation_waivers.pores`): `ne_mid` is a true distance only on the nuclear side. Grid nodes outside
the HMFF midsurface mesh hold the 120 nm clamp, so the zero lies ~8-10 nm inside the mesh and the ONM (`ne_mid = +17`)
only a few nm outside that zero. Result: perinuclear volume 2.744 µm³ instead of ~8.3 µm³.

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| J1 | New package, never edit v7 | `toolkits/em` package step N -> N+1 (`pkgstep.open_pkg` / `publish`) from `em_cell1_v7`. Only `fields/ne_mid.npz` is recomputed: exact signed distance to the unchanged HMFF mesh (`meshes/ne_mid.ply`) on both sides, with the sign from ray-parity occupancy of the closed mesh. Same grid origin, brick layout and ±120 nm clamps. | every other content file (6 fields, 3 sections, pores, mesh, layer files) has the same sha256 as in v7; engine params unchanged; v7 folder untouched (manifest sha256 `eb16248b…` unchanged) |
| J2 | cellgeom validator | `cellgeom.validate` (written by `finalize`) | `pass` true. **`pores` passes with no waiver**: max \|ne_mid\| at pore centres < 4 nm (h/4), axis-vs-gradient median < 10°. `eikonal.ne_mid` median 0.95-1.05 **and** p95 ≤ 1.10 (v7 p95 7.09). `field_volumes_areas.ne_mid` passes against the mesh volume 271.2 µm³ and mesh area 245.1 µm². The only remaining waiver is `field_volumes_areas` for the unchanged legacy summary of er / endo / golgi. `label_rule` agreement 1.0. `layers`: every `depends_on` = `jurkat/em_cell1_v8`. |
| J3 | Field = true distance | 200,000 seeded points within 100 nm of the mesh: `ne_mid_v8` vs the exact signed distance to the mesh triangles | \|Δ\| p99 ≤ 1.0 nm; median \|ne_mid\| at pore centres ≤ 0.5 nm |
| J4 | Perinuclear volume | MC label volume `perinuclear_space` | expected **8.08 µm³** = 245.1 µm² × 34 nm (= 8.33, the plan's "~8.3") − 1,000 pores × 2.54e-4 µm³ (the opening per pore: a cylinder of r = PORE_R + PERI_W/2 = 62 nm × 34 nm, minus the inner half of the 17 nm torus tube). **Accept 7.7-8.5 µm³** (v7 2.744) |
| J5 | Other regions: direction and size of change (MC, v7 value in brackets) | MC label volumes v7 vs v8 | `cytosol` decreases by 4.0-7.0 µm³ (392.81; the ONM moves ~25 nm outward). `nucleoplasm` \|Δ\| ≤ 1.0 µm³ (258.87; the INM is unchanged because the nuclear side of the v7 field was exact). `nuclear_pore_channel` 0.30-0.60 µm³ (0.352; expected ~0.43). `er_lumen`, `mito_ims`, `mito_matrix`, `golgi_lumen`, `golgi_intercisternal`, `endolysosome_lumen` do not increase: each Δ ≤ +2·CI and ≥ −(2 % + 2·CI) (lumens within ~25 nm outside the mesh become NE lumen, since NE comes first in the label rule). `nucleolus`, `nuclear_speck`, `exterior`, `outside_imaged_volume` unchanged within 2·CI |
| J6 | Change confined to the NE | 2,000,000 uniform points, labels v7 vs v8 | 100 % identical where \|ne_mid_v8\| > 150 nm; every changed point lies within 150 nm of the NE midsurface |
| J7 | Areas | validator shell areas | NE midsurface 245.1 µm² (mesh) within 5 %; PM, ER, mito, Golgi, endo, nucleolus areas and volumes identical (identical field files) |
| J8 | NPC / pores unchanged | `pores.npz`; params PORE_R, PERI_W, PORE_ZWIN, PORE_LW | byte-identical: 1,000 pores (300 EM + 700 random), same params. The NPC count/density decision is NOT part of B6 (plan §5c item 1) |
| J9 | Reference stays selectable | config | `config/cells/jurkat.yaml` names `em_cell1_v8`; `config/cells/jurkat_v7_port.yaml` (`extends: jurkat.yaml`) names `em_cell1_v7` with its v7 declarations. The B1 check (`tests/check_b1_geometry.py`) on v7 still passes |

## Part Y: yeast geometry (new package, per-membrane half thickness)

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| Y1 | Cited half thicknesses | manifest `engine.params.MEM_HALF` = `membranes[].half_thickness_nm` = D_P-P / 2 from **Pogozheva et al. 2022 J Chem Inf Model 62:1036, SI Table S20** (fungal systems where they exist): PM = PMf 42.2 Å -> **2.110 nm**; ER, NE = ERf 37.9 Å -> **1.895**; GOLGI = GOLf 39.6 Å -> **1.980**; MITO_OM = MOM 39.5 Å -> **1.975**; MITO_IM = MIM 39.1 Å -> **1.955** (mammalian; the table has no fungal mitochondrial system). VAC, PEROX, LD, CELL_OUTER keep **2.0**: the table's VAC is the plant tonoplast, there is no peroxisome system, LD is a monolayer and CELL_OUTER is the cell-wall surface, not a bilayer. HALF stays 2.0 | exact values; schema check half_thickness_nm == MEM_HALF |
| Y2 | New package, never edit the old | package step from `analytic_x5_s0_f5f3a5c723b2`: only the manifest changes (params, membranes, generator.parent) | every content file has the same sha256 as the parent; the parent folder is untouched |
| Y3 | Packer clearance uses it | rule table `yeast_cell3d_v1.yaml` `clearance: thickness_corrected` (no wall exemption). A package whose MEM_HALF all equal HALF compiles to the legacy `plain` rule, so the X5 package keeps a uniform 2 nm, bitwise. The legacy core config `config/cells/yeast_x5_core.yaml` (`extends: yeast.yaml`) names the old package | old package: `min_mem_dist_g` identical to the A9 code on 1,000,000 points (bitwise). New package: \|Δclearance\| = \|ΔMEM_HALF\| of the nearest membrane, max ≤ 0.11 nm |
| Y4 | Validator | yeast package checks (`toolkits/analytic/yeast.py` validate: checksums, schema, native round trip, prims / fields / pores vs native) + schema half_thickness_nm == MEM_HALF | all pass |
| Y5 | Labels and volumes unchanged | 2,000,000 points, old vs new package, core evaluator | labels 100 % identical; MC volumes identical; band volumes (pack.regions shells, defined on midsurfaces) unchanged. The exclusion boundary of each membrane moves by \|ΔMEM_HALF\| ≤ 0.11 nm (plan: ≤ 1-2 nm) |

## Part C: chromatin and the other frozen pre-placed layers (Jurkat)

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| C1 | Frozen rep5 inside the new nucleus | `pre_chromatin_v8_rep5` (md5 cd4c3b4c…; 28,190,641 NUCLEOSOME of 51,387,965 instances): share of NUCLEOSOME centres whose v8 label is in {nucleoplasm, nucleolus, nuclear_speck} (inside the INM). Reported with the same share for v7, for all instance types, and the labels of the centres outside | **≥ 99.9 %** |
| C2 | Trim rule (applies only if C1 fails) | write a NEW layer file `chromatin_cell1_b6trim.npz` + json (new frozen input, md5 recorded; rep5 untouched): rep5 minus every instance of any type whose centre's v8 label is outside {nucleoplasm, nucleolus, nuclear_speck}, with all per-instance arrays filtered row-aligned and the consumed counts recomputed. Refit deferred (DECISION 2026-11-15) | trimmed file passes C1 at 100 %; the number trimmed is reported |
| C3 | Other frozen layers | polysomes, cytoskeleton, npc, lamina (frozen copy), mrna and mrna_nascent beads: v7 -> v8 label transitions of the centres | reported. **Rule:** pre-placed instance centres that fall in `perinuclear_space` under v8 (cytosolic items within ~25 nm outside the mesh) are dropped at indexing by a config index drop (`ne_lumen_b6`, mode `inside_labels`, all pre-placed libraries). The count is reported per module and expected ≤ 1 % of every module. RNA beads are only reported (obstacles) |
| C4 | `depends_on_geometry` updated | every frozen pre-placed module of `jurkat.yaml` declares `depends_on_geometry: jurkat/em_cell1_v8` (`jurkat_v7_port.yaml`: `jurkat/em_cell1_v7`); the preplaced stage fails the build on a mismatch; suite C16 `depends_on_geometry_mismatch` (was n/a) | 0 mismatches. Yeast: every module is built in the build on its own package, so 0 |

## Part B: build-level acceptance (one build per cell from the B6 commit, `scripts/build_isolated.sh <cell> b6 24`)

| # | Criterion | Pass |
|---|---|---|
| B1 | Build complete; geometry stage validation passes; `geometry/package.json` names the new package | yes |
| B2 | Common suite, `no_regression.tsv` vs the reference build's suite.json with the §4 tolerances | no regression outside the EXPECTED list below; any other regression needs a dated DECISION |
| B3 | J4-J6 / Y5 numbers recomputed on the build (`geometry/measures.json`, stream `geometry.measures`) | the same bounds as J4, J5 / Y5 |
| B4 | Chromatin (C1) and the index drop (C3) seen in the build: pack index log line `index drop 'ne_lumen_b6'` | count ≤ 1 % per module; C13 `prepass_md5_ok` = 1 (no trim) or the trimmed file's declared md5 |

**EXPECTED changes (Jurkat; caused by the fix, so not regressions within these bounds):**

- C11 `density.cytosol`: up, by at most +3.0 mg/ml. The same packed mass goes into ~1.4 % less cytosol (162.8 × 5.5/392.8 ≈ +2.3). Mass-consistency check: density × region volume within ±1 % of the reference.
- C11 `density.perinuclear`: the region triples in volume. Report only, since no species maps to it (reference 0).
- C11 `nuc_cyto_ratio`, `phi.cytosol`: report (down / up by ~1.4 %).
- C11 densities of ER_lumen, matrix, golgi_lumen, endo_lumen: within ±2 mg/ml, or mass-consistent with J5 volume losses ≤ 2 %.
- C1 / C2 / C3 placed metrics: ≤ 0.01 (standard tolerance; no exception expected).
- C10 `placed_mass_fraction` and `placed_copies_fraction` ≥ reference − 0.003 (standard; a denser cytosol may cost ≤ 0.002).
- C7 / C8 / C9 / C12: standard tolerances. The C3 index drop keeps pre-placed centres out of the new NE lumen.
- C13 frozen chromatin metrics unchanged (frozen file).
- C16 `depends_on_geometry_mismatch` = 0 (new row).
- The v7 reproduction comparison `compare_v7` (still in the DAG) is expected to report deviations around the NE. It is informative only, not a B6 criterion.

**EXPECTED changes (yeast).** Labels, volumes, band volumes, pre-placed builders and inputs are unchanged. Packing
differs at the seed-noise level everywhere: every tile with a membrane sees a clearance change of ≤ 0.11 nm, so the
tile RNG trajectories diverge. A §4 tolerance breach in a yeast metric is therefore **not** a regression if the B6 value
lies inside the range of the existing X5 port seed builds s0 / s1 / s2, widened by the §4 tolerance. Their suites are
computed into `_dev_unify/b6/seedspread_suite/` from the existing builds; no new replicate builds are run (Pedro). The
expected direction is for membrane-related C7 / C8 rows not to worsen beyond that band.

**Verdict.** `<build>/validate/compare_B6/verdict.json` per cell (`scripts/b6_geometry_check.py verdict`):
`ACCEPTED` if parts J/Y, C and B1-B4 pass and B2 has no regression outside the expected list. Otherwise
`NOT ACCEPTED` with the failing items. Pre-build measurements are in `{results_root}/<cell>/_dev_unify/b6/`.
