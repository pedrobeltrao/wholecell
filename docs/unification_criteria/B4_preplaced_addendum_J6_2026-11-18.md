# B4 pre-placed builders: addendum 2026-11-18 (J6: ATP synthase dimers on the rounded crista rims)

Addendum to `docs/unification_criteria/B4_preplaced.md` (frozen, md5 e715c2f90e17653369b7d5ac2c4d3d96, not edited) and
`B4_preplaced_addendum_2026-11-17.md` (md5 ad6dbb8bf4395efd251b6f8ff5ec6574, not edited). Geometry: Jurkat package
`em_cell1_v8r` (rounded crista rims, CR_RIM 10 nm; `B6b_crista_rims.md`, `B6b_crista_rims_addendum_2026-11-17.md`).
**Written after seeing the results** (B6b R9 / R10 and the measurements below), i.e. post hoc. Decision by the main
session 2026-11-18 (Pedro informed): (a) J6 is evaluated at a 4.5 nm transmembrane tolerance on the rounded rims;
(b) the rim sampling of the builder is fixed so that the dimer density reaches the literature range; no kinked-rim
geometry. Measurements: `{results_root}/jurkat/_dev_unify/b4/j6_fix/` (3 cores).

## 1. J6 transmembrane tolerance: 4.5 nm on the rounded rims

**What changes.** Frozen J6: "every dimer's TM residues within 2.5 nm of the IM midsurface". With this addendum the
bound is **4.5 nm** for the Jurkat lamellar cristae (em_cell1_v8r and later rounded-rim packages; config
`tm_tol_nm: 4.5`). The other J6 items are unchanged: dimers > 0 and ≤ pool / 2; no more than 2 % of the tested
(every third) soluble residues inside a bilayer (yeast tolerance); mtDNA items unchanged. Yeast is not affected (its
crista prims keep 2.5 nm). Reported with J6 from now on (`scripts/b4_j6_atp_check.py`): dimers per µm of rim (traced and
band rim length) against the literature 50-83 per µm, and the overlaps of the dimers with every pre-placed layer.

**Why.** The builder's transmembrane residues are the monomer beads with |z| < 2 nm of the library membrane plane: a
flat patch of ~9.4 × 9.8 nm per monomer. A rigid flat patch that wide cannot follow a smooth curved rim: fitting the
dimer's TM beads to an ideal cylinder (any axis offset, radii 10-17 nm, dimer angles 84-180°) leaves a best maximum
deviation of 2.5-2.7 nm, and real rim sites of em_cell1_v8r add 0.3-0.5 nm (best sites 3.0-3.3 nm; B6b addendum). In
cells the membrane bends sharply at the dimer interface and is flat under each monomer; a smooth half-cylinder rim is
a different shape. 2.5 nm is geometrically unreachable for this dimer model on any smooth rim.

**Measured counts by tolerance** (yeast dimer geometry 84°, 8.8 nm; ≤ 2 % soluble in membrane):

| TM tolerance | B6b: random sampling, cell-wide (400,000 samples) | B6b: random sampling, saturated, per µm (band) | rows sampling (§2), 3 mitochondria, per µm (traced 28.1 µm) |
|---|---|---|---|
| 2.5 nm (frozen) | 0 | 0.1 | 0.21 (6 dimers) |
| 3.0 nm | 14 | 1.1 | 4.8 (135) |
| 3.5 nm | 674 | 13.5 | 43.2 (1,214) |
| 4.5 nm (this addendum) | 1,809 | 26.2 | 56.7 (1,595) |

(B6b subset: lamellae 84, 35, 42 of em_cell1_v8r, the R10 sub-region; `j6_fix/tol_scan_subset.json`.) Only 4.5 nm
reaches the literature density; at 3.5 nm the rows sampling gives 43 per µm (below 50).

## 2. Rim sampling: dense rim rows (builder change)

**Before (B4 `random`).** `params.samples` (400,000) random points of the id-grid bounding box were projected onto
rim lines: ~3,950 rim points for a 4,040 µm rim, i.e. ~1 per µm; with a 13 nm exclusion around kept sites the
cell-wide density was 0.45 per µm (1,809 dimers), and even a saturated sub-region reached 26 per µm only.

**After (`rim_sampling: rows`, `core.atp_synthase_dimers`, lamellae only; config-driven, seeded):**

1. Seeds: `seeds_per_cell` (1.0) uniform seeds per enabled cristae id cell (not the bounding box); matrix seeds are
   projected onto their rim line (the B4 Newton steps on the boundary-IM distance within the slab mid-plane).
2. Rim lines: from each seed not on an already traced line, the whole rim line is traced in both directions every
   `trace_step_nm` (2 nm), each point projected back onto the rim; a line ends where the projection fails, the line
   leaves its lamella, jumps, or closes on itself. 1.0 seed per cell finds every rim line of the test sub-region
   (traced length identical at 0.8, 3.2 and 9.6 seeds per cell).
3. Rows: lines are walked in seeded random order; at each trace point the poses are tried in order
   (`tilt_deg` [0, -10, 10, -20, 20] × `rho_nm` [10, 9, 11, 8, 12, 7, 13, 6.5]); the first pose that does not overlap
   a placed dimer or an obstacle (`obstacles: [npc]`) and passes the unchanged site test is kept, then the walk moves
   `spacing_nm` = **12 nm** on (dimer spacing along the rows: ~12 nm, Strauss et al. 2008 EMBO J 27:1154; 12 to
   > 20 nm, Davies et al. 2011 PNAS 108:14121). Closed lines stop one spacing before their first dimer.
4. Tilt about the rim tangent (the dimer's ridge): the rims of lamellae that meet the boundary IM obliquely are
   sheared, so their cap normal leaves the slab plane. Measured on 300 rim points of the test sub-region: the in-plane
   pose (tilt 0) passes 4.5 nm at 63 %, tilts 0 / ±10 / ±20° at 94 %. Without tilts the rows sampling gives 37.5 per
   µm (13 nm spacing), with them 56.7.
5. Speed only (same decision): the site test runs in one numba kernel per point (`first_ok`: only the tested beads,
   adaptive bead order, early exits; pose chunks in parallel, the first passing pose index is returned, so the result
   does not depend on the thread count); verified identical to the per-pose Python loop on the test sub-region.

The B4 `random` mode stays the default of the lamellae path (unchanged code path); the yeast crista prims path
(`place_prims`, `_make_site_ok`) is byte-identical.

## 3. Results (Jurkat, whole cell)

Serial pre-placed stage of `config/cells/jurkat.yaml` (all six builders, em_cell1_v8r), 3 numba threads.
`stage/` = rows sampling without obstacles; `stage2/` = the committed config (+ `obstacles: [npc]`).
Measures: `scripts/b4_j6_atp_check.py` (`atp_check*.json`), `scripts/b4_preplaced_check.py measure`
(`measure_stage*.json`, J5b addendum per class).

| | B4 / B6b (`random`, 4.5 nm) | `stage` (rows) | `stage2` (rows + npc obstacles; committed) |
|---|---|---|---|
| dimers (monomers) | 1,809 (3,618) | 205,165 (410,330) | **205,152 (410,304)** |
| share of the CPX-6151 pool (868,688) | 0.4 % | 47.2 % | **47.2 %** (dimers ≤ pool / 2: yes) |
| rim lines / rim points / sites tried | - / 3,967 / 3,945 | 2,083 / 1,951,727 / 850,926 | 2,083 / 1,951,727 / 850,996 |
| per µm of rim, traced (3,902.3 µm) | 0.46 | 52.6 | **52.6** |
| per µm of rim, band (4,040 µm) | 0.45 | 50.8 | **50.8** |
| nearest dimer centre, median (p05-p95) | - | 15.1 nm (11.9-16.3) | 15.1 nm (11.9-16.3) |
| poses used: tilt 0 / ±10 / ±20° | 100 % / - / - | 45 / 29 / 25 % | 45 / 29 / 25 % |
| builder time (3 threads) | 52 s | 1,296 s | 1,278 s |
| serial stage (6 builders), peak RSS | - | 2,451 s, 48.9 GB | 2,283 s, 48.9 GB |

**Rim length.** Two measures: *traced* = the total length of the rim lines the builder traced (cap line, b = -CR_GAP,
2 nm steps), 3,902 µm in 2,083 lines; *band* = the B6b estimate (rim band volume / cross-section, 4,040 µm). Both
densities are reported; both are in the literature range 50-83 per µm (Strauss 2008; Davies 2011).

**Overlaps.** ATP-ATP: 0 monomer pairs between dimers with an overlap > 0.5 nm (J5 definition), 2 at > 0.3 nm (the
pack tolerance of the builder's test; float32 export rounding), 0 within a dimer. ATP-NPC: 23 pairs > 0.5 nm in
`stage` (pores next to mitochondria), hence `obstacles: [npc]`; `stage2`: 0 pairs > 0.3 / 0.5 nm (5 sphere pairs touching within the pack tolerance). ATP-polysomes: 0. Cytoskeleton, MT,
chromatin, lamina, mRNA: none within reach of a dimer. ATP-mtDNA: 2,752 dimer-sphere / bead pairs > 0.5 nm involving 785 monomers (0.19 %) in `stage2` (2,750 / 785 in `stage`): the mtDNA builder (unchanged) tests
its walk steps (3.4 nm) against the dimers but exports 2-bp beads interpolated along the steps, which can clip dimer
spheres; it is a property of the mtDNA walk, already present in B4 (8 bead pairs with 1,451 dimers), reported, not
changed here.

**Membrane distances** (seeded sample of 4,000 dimers, exported float32 poses, re-evaluated with the evaluator in an
independent kernel): transmembrane beads |IM| median 1.13 nm, p90 2.60, p99 3.58,
max 4.50; per-dimer maximum median 3.74 nm; dimers with all TM beads within 2.5 / 3.0 / 3.5 / 4.5 nm: 0 / 1.0 / 36.4 /
99.98 % (the one exception: 4.5000005 nm, float32 rounding of the exported pose; positions near 2 µm have a float32
step of 1.2e-4 nm); tested soluble beads within 0.6 nm of a bilayer: median 0.89 %, max 2.01 % (1 dimer: one bead
3e-5 nm inside the 0.6 nm threshold after float32 rounding; the builder's float64 test passes); soluble beads outside
the matrix / IMS: 0. J6 geometry tests at the addendum tolerance: 100 % (up to float32 export rounding).

**J5b (addendum 2026-11-17 rule, 300 v7repro2 tiles).** No overlapping pair involves the ATP synthase layer (18,638 ATP monomers on the sampled
tiles, 0 pairs in any class: `stage` and `stage2`). The other classes are unchanged from the 2026-11-17 addendum stage
(MT-MT 201,450 pairs, polysomes-polysomes 2,029, chromatin-NPC 33, ...); the only failing class is lamina-NPC (2,511
pairs with the frozen v7 lamina; the build uses the rebuilt lamina, unchanged by this addendum). All pre-pass pairs
2,130,904 (v7 layers 2,200,254), 816,761 vs 842,927 per 10⁶. J3a 0.684, J3b 0.441, J4 35.22 nm: as in the 2026-11-17
addendum (the polysome walk does not see the dimers). mtDNA: 479 nucleoids, 16.68 of 16.69 Mb placed (J6 ≥ 95 %).

## Checks of this addendum

- Jurkat: stage and stage2 runs, `j6_fix/run_stage.sh`, `run_stage2.sh` (+ code diff / head), measures listed above.
- Yeast (Y1): serial pre-pass of `config/cells/yeast.yaml` with this code (`yeast/_dev_unify/b4/run_vcell_j6_fix/`,
  `run_check_j6_fix.sh`): builders 1-18 identical to s0, 19-28 (atp_synthase included) identical to `ref_yeastenv_patched`, pre-pass log
lines identical (3 threads, 2,684 s; run with the code before `params.obstacles` was added, which touches the
lamellae path only). The lamellae `random` mode of the new code reproduces the B4 code exactly on
em_cell1_v8r (1,809 dimers, all arrays identical; whole cell, B4 params).
- `tests/test_preplaced_b4.py` (new: `test_atp_rows`), `tests/test_preplaced_serial.py` pass.
