# Yeast X5 reproduction: pre-declared criteria

Written 2026-11-16, **before** any wholecell yeast build was run (template: `docs/v7_reproduction_criteria.md`).
Committed to git before the run; its md5 is recorded in `PLANNING_NOTES.md` (DECISION 2026-11-16, yeast X5 port).
Changes after the run are allowed only as a dated addendum below the line at the end, never by editing the criteria.

## 0. What is compared

- **Reference (X5)**: yeast candidate build X5, run folder
  `{yeast_root}/runs/X5_seed0/` (READ-ONLY; owned by the yeast project, never written, run or locked):
  build_id `20261116-141711-e3b7448-s0`, commit `e3b74481d8` (branch round-x), `dirty: false`, master seed 0,
  host eu-g6-037, `YC_NPROC=40`, GPUs of that node; finished and validated 2026-11-16 16:33 ("no regressions against
  build W"). Compared files: `work/` (geometry, inventory, prepass, tiles, validation) and `lib/` (structure library).
- **Port (P)**: the wholecell build of `config/cells/yeast.yaml` (pipeline engine `yeast_x5`, packer `yeast_rsa_v1`,
  master seed 0) into a new folder `data/derived/wholecell_builds/yeast/<build_id>-x5repro/`. The yeast stage code is
  the vendored copy `plugins/yeast/engine/` (= `git archive e3b7448`, verified identical to `X5_seed0/tree` for
  pipeline/, scripts/, validation/, simgpu/, analysis/, web/, env/) with the documented port patches only
  (`plugins/yeast/engine/PORT_PATCHES.md`): geometry read from the cellgeom_v1 package, GPU count and lock from the
  environment, PDF step non-fatal. Python environment: a byte copy of the yeast bundle's `.venv` (CPython 3.11,
  numba 0.61.2, numpy 2.2.6, scipy 1.15.0, scikit-learn 1.9.1, cupy-cuda12x 14.2.0) at
  `conda_envs/yeast_x5`, run with a uv-managed CPython 3.11.
- **Inputs declared frozen** (checksums in the config / frozen manifests, verified before use):
  - yeast `pipeline/data` (the bundle's, read in place; md5 manifest of every file; files changed after X5 started,
    all under `validation/literature/endocytosis/`, are listed and are not read by the build);
  - `{yeast_root}/inputs` large inputs (AFDB v6 tar, AFDB dimers, AF3 pool summary, IntAct; md5);
  - the RCSB answer cache = X5's `work/rcsb_cache.json` as it was at the end of X5 (X5 started from the bundle cache and
    added its own queries at 14:21; using its final state makes stage 1 offline and identical by construction);
  - PDB assembly files of the bundle `work/pdb` (read in place through a link farm; md5 manifest);
  - `work/organisation/pool_pairs.parquet`, `work/interactions/pair_summary.parquet` (no generator in the yeast repo),
    `work/sweeps/` start copy (the bundle copy X5 started from), `work/validation/{image_nc_raw,proximity_sets}.tsv`;
  - **library start set**: the **342** entries of X5's `lib/` that X5 did not rewrite (modification time before the
    build stamp 14:17:11: copied in from the bundle library, never rebuilt). They include the 42 Boltz-2 tRNA species
    (manual GPU step, `trna_species.py`), the builder templates (`F_actin`, `TPM1_coil`, `Arp23`, `CP_cap`,
    `MT_tubulin`, `gTuRC_*`, `ATP_dimer`, `Pil1_dimer`, `Lsp1_dimer`, `Vac8_Nvj1`, `RVS161_RVS167_het`,
    `mitoribosome_74S_IM`, `preribo_*`, `ribosome_80S_{RAC,Sec61,Ssh1}` and their JSON/NPY side files) and stale
    leftovers. The port starts from exactly this set (copied, md5 manifest); every other library entry is rebuilt.
    (The audit of build S counted 305 such entries; X5 has 342.)
- **Hardware difference (stated in advance).** X5 ran on eu-g6-037 (40 workers, RTX 4090 GPUs per HANDOVER); the
  port runs on eu-g4-002 (AMD EPYC 7742, Zen 2, AVX2 only; one RTX 3090; ≤ 32 workers). The yeast bit-identity claim
  (rebuild of build R, 56 vs 64 workers) was established on one node class with one venv. numba kernels are compiled
  for the host CPU with `fastmath=True` and OpenBLAS selects CPU-specific kernels (KMeans sphere trees, PageRank), so
  bit-identity across CPU generations is **a test, not an assumption**. Tile seeds are independent of the worker
  count, so the worker count alone must not change any result.

## 1. Expected outcome and the two-level verdict

With identical code, inputs, library start set and seed, and the same package versions, the port is expected to be
bit-identical to X5 in every stage that has no CPU-dependent floating-point reduction; for the rest the exact level
is tested and the tolerance level applies.

- **Level E (exact)**: the artefact is identical (array equality / identical parsed content; npz zip timestamps,
  gzip headers, build ids, timings and absolute paths are not content).
- **Level T (tolerance)**: if not exact, the artefact is within the tolerances below.

Verdict: **REPRODUCED** if every *must* row passes at level E or T; **REPRODUCED WITH DEVIATIONS** if all *must* rows
pass but some *report* rows fail; **NOT REPRODUCED** if any *must* row fails. Each failure is explained. If the port
fails only at level E from the first numba/BLAS-dependent step on, the report states the first diverging artefact.

## 2. Geometry

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| G1 | `yc3d/geometry.json` (G vector, pores, PT/PQ/PC primitives, bounding spheres, npc_pos, mitochondria, consts) | parsed content identical | same keys, same primitive and pore counts; numeric \|Δ\| ≤ 0.01 nm (angles ≤ 1e-6 rad) | must |
| G2 | `yc3d/geometry_fields.npz` (real-shape SDF bricks: cell outline, vacuole, nucleus) | arrays identical | same brick index; values \|Δ\| ≤ 0.01 nm at the 99.9th percentile | must |
| G3 | cellgeom_v1 package (new; no reference): checksums, schema (core ids, extensions ≥ 1000 with a core parent, engine codes 0..63), native round trip (package → G, pores, PT/PQ/PC identical to the builder output), prims table ↔ native primitives identical, MC label volumes vs s2 log | all checks pass | — | must |
| G4 | tile list and per-tile measures: `tile_measure.parquet` (27 region columns), `tile_labvol.parquet` | identical | same tile set (5,912); per-region totals within 0.1 %; ≥ 99 % of tiles within 0.5 % of a tile volume per column | must |
| G5 | per-tile membrane meshes `yc3d/geom/*.npz` | arrays identical | membrane area per membrane within 0.1 % | report |

## 3. Inventory and copy numbers

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| I1 | `results/inventory/ingredients.csv` | identical content | same ingredient set; copies \|log2 ratio\| ≤ 0.01; compartment weights \|Δ\| ≤ 1e-6 | must |
| I2 | `features/protein_features.parquet` (6,067 proteins) | identical | `copies_per_cell` rel. Δ ≤ 1e-9; other numeric columns rel. Δ ≤ 1e-6; strings identical | must |
| I3 | shared core copy-number step (`wholecell.inventory.copies`, config `inventory.copies`) run in the build vs the port engine's `copies_per_cell` | max rel. Δ ≤ 1e-12, same 6,067 proteins, NOT_IN_THIS_CELL categories and shares as the X5 s1 log | — | must |
| I4 | s1 log numbers: NOT_IN_THIS_CELL counts/shares, ribosomes 262,722 (factor 1.185, protein 3.61 pg), curated complexes 11 / legacy homo-oligomers 22, mass in assemblies, pre-ribosomes 20,404 | identical lines | each number within 0.1 % | must |

## 4. Assembly state (stage 1)

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| A1 | `oligo_map.json`, `f_actin.json`, `complexes_tried.csv`, `complexes_added.csv`, `heterodimers.csv`, `oligomers.csv` | identical content | same rows/keys (sets identical); numeric rel. Δ ≤ 1e-9 | must |

## 5. Structure library

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| L1 | set of library files after the build vs X5 `lib/` (7,112 files) | identical name set | — | must |
| L2 | every entry rebuilt by the port (all except the 342 start entries): npz arrays, json, npy | arrays / content identical | same keys and shapes; beads \|Δ\| ≤ 1e-4 nm; sphere-tree centres \|Δ\| ≤ 0.05 nm and radii \|Δ\| ≤ 0.05 nm (KMeans on BLAS) for ≥ 99.9 % of entries, all entries same tree size | must |
| L3 | `results/library_af.json` (5,616 AF models written, 0 failed, 11 without a model) | identical | — | must |
| L4 | the 342 start entries unchanged at the end of the build (md5) | identical | — | report |
| L5 | `yc3d/lib_ids.npy` (type order), `libpack_meta.npz`, `pack_trees_15.npz` | identical | lib_ids identical (must); radii/centres \|Δ\| ≤ 0.05 nm | must |

## 6. Pre-placed structures (`yc3d/prepass.npz`, stage 3 serial pre-pass)

Groups: NPC (`npc`), Rabl chromatin (`nuc`, `seg_*`, `chain_target`), rDNA/nucleolus (`rpol`, `r90`, `rnuc`,
`rdna`, `prerrna`), microtubules/SPB/kinetochores (`mt`, `spb`, `kt`, `kto`), NPC FG mesh and receptors (`fg`,
`kap`), Pol II / origins / chromatin sites (`pol2*`, `mcm`, `orc`, `tf`, `lm`), eisosomes/MCC (`eis`, `mcc`), actin
cables and endocytic sites (`act`, `tpm`, `coat`, `rvs`, `pat`), cristae ATP synthase dimers and mtDNA (`atp`,
`mtdna`), vacuole Cvt/NVJ (`cvt`, `nvj`), polysomes/mRNA/mRNP (`poly`, `mrna`, `rnp`).

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| S1 | array key set (128 arrays) | identical | — | must |
| S2 | instance / bead counts per group (length of every `*_pos`, `*_pts`, `*_off` array) | identical | each within max(3, 1 %) | must |
| S3 | s3 log structure numbers (MTs 23 / 14,261 dimers; Rabl 33 segments, 64,487/64,487 arm nucleosomes; rDNA 75 genes, 3,750 Pol I, 5,227 rDNA nucleosomes; FG 11,904 domains; Pol II 515/654; origins 344/340; chromatin sites 20,230; eisosomes 22 furrows, 3,671 dimers; cables 79 filaments, 42,179 actin; patches 28; ATP synthase dimers 2,048; mtDNA 14 nucleoids; Cvt 259 Ape1; NVJ 518; polysomes 39,314 chains / 182,329 ribosomes; RAC 63,629; mRNP 42,766 Pab1) | identical lines | each count within 1 % | must |
| S4 | positions / orientations of every pre-pass array | arrays identical | per group: centroid within 20 nm, RMS radius about the cell centre within 1 % | report |

## 7. Packing (random sequential addition, squeeze, organisation layer)

| # | Artefact | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| P1 | `yc3d/expected_totals.parquet` (requested copies per library type × region) | identical | per pair rel. Δ ≤ 1e-6 | must |
| P2 | organisation layer `organisation/network.npz`, `partners.npz` | arrays identical | partner CSR identical; float arrays rel. Δ ≤ 1e-5 | must |
| P3 | `yc3d/tile_counts.parquet` (requested per tile × type × region) and `thinning.json` | identical | per (type, region) totals within max(3√n, 0.5 % n); thinning factors ± 0.005 | must |
| P4 | instances per tile (`pack_summary.parquet` placed, all 5,912 tiles) | identical | — | report |
| P5 | placed copies per (library type, region), whole cell (from `tiles_raw`) | identical | within max(3√n, 1 % n) for every pair with n_X5 ≥ 10; ≥ 99 % of pairs with n_X5 < 10 within ± 3 | must |
| P6 | totals: placed / requested (s3 log 25.17 M / 25.37 M = 0.9924), `instances_total` 25,664,134 | identical | placed fraction ± 0.003; instances ± 0.5 % | must |
| P7 | positions / orientations per tile (`tiles_raw/*.npz` types, pos, quat, job table) | arrays identical | not required | report |

## 8. Validation suite (the yeast code, run by the port on the port build)

| # | Metric | Level E | Level T tolerance | Class |
|---|---|---|---|---|
| V1 | stage-4 `results/wholecell3d_validation.json`: chromatin, polysomes, MTs, ATP synthase, eisosomes, clash rate, containment, orientation, TM residues, densities, mass by source, placed/expected per region | numeric content identical (excluding build id) | counts as S2/P6; clash and containment rates within max(25 % rel., 1e-5 abs.); orientation median ≤ 0.01°; densities ± 2 mg/ml; fractions ± 0.005 | must |
| V2 | CPU comparisons of `validate_all.sh` (`work/validation/*.json`, `validation/latest/*.json`): localisation benchmark (dev, test), proximity, PPiSeq, Singh 2020, kinase localisation, ribosome cryo-ET, chromatin stats, Hi-C, nucleosomes, mRNA stats/occupancy/localisation, tRNA-nucleolus, DNase, silent chromatin, image N/C, ER ribosome density, literature checks, copy numbers, **anchors**, placement checks, vacuole surface, membrane proteomics, metabolism | identical numeric content | every scalar within max(1 % relative, 0.01 absolute) (AUROC, Spearman, fractions: ± 0.01) | must |
| V3 | GPU steps (not bit-reproducible: float32 CuPy on a different GPU model, 1 GPU instead of 4): contact pass (`contacts/contacts.npz`), kinase encounters, diffusion sweeps of this build (`sweeps/{cyto120_mf,gem200_mf,nuc100_mf}.json`) and the diffusion benchmark | identical | contact pair count within 1 %; encounter statistics within 5 %; sweep volume fraction within 2 × the X5 box SD (0.010 / 0.013) and nuc100 φ ± 0.01; diffusion benchmark D/D0 per probe within 10 % | must |
| V4 | condensate negative control (`condensates/cond_test.log` log-phase and stress lines: n_species 2,285, n_pb 66, onset) and simulated microscopy counts | identical | counts within 1 %, onset ± 0.05 | report |
| V5 | `build_stamp.py check`: every result file is of this build, 0 older | 0 stale | — | must |

## 9. Shared-pipeline outputs (new; no X5 counterpart)

| # | Check | Pass | Class |
|---|---|---|---|
| X1 | YCT2 export of the port (tile_format_v2): decoded instance count = packed + pre-placed source count; position error ≤ q_range/131070; rotation error ≤ 1e-4 rad; quaternions w-first with w ≥ 0; region ids from `regions.json` (config `pack.regions` + pre-pass groups, > 32 regions allowed through the u16 palette); grid origin/extent from the geometry, not `TGRID = 26` | all pass | must |
| X2 | the yeast stage-5 viewer encoding (`yc3d/web/data`): decompressed tiles and library files vs X5 | identical | report |

## 10. Reporting

`validate/compare_x5/` in the port build holds `criteria_results.tsv` (row id, metric, X5, port, Δ, tolerance,
level E/T/fail, class) and `verdict.json`. The comparison script is `plugins/yeast/repro/compare_x5.py`; its
thresholds are constants that must match the tables above (reviewed at commit). Before the run it is tested on X5
against itself (every row must pass at level E).

---
(Addenda after the run go below this line, dated.)

**Addendum 2026-11-16 18:45 (during the run; criteria unchanged).**
1. Self-test before the run: `compare_x5.py --self-test` (X5 against itself) gave level E on every row it can
   evaluate (G1-G5, I1, I2, I4, A1, L1-L5, S1-S4, P1-P7, V1-V5).
2. Observed before packing: `geometry.json` differs from X5 by ≤ 1.5e-12 nm (PQ, bc, br, pores, mitochondria; ULPs of
   the CPU-specific numpy/numba code on Zen 2 vs the X5 node), while the SDF bricks, tile measures (27 × 5,912) and all
   5,912 mesh tiles are bit-identical; copies differ by ≤ 4e-16 relative. So level E is lost at stage 2 by hardware,
   as anticipated in §0; the RSA packing is chaotic, so the packed positions are expected to differ (level T).
3. Comparison-script fix (commit 8f1d0e1): I2 and A1 compared rows by position; ULP-level copy differences swap the
   order of tied rows, which is not content. Rows are now aligned by key (ORF; first column of each table). The
   workflow's own `validate/compare_x5/` (script of aabe7b5) is kept; the verdict is `validate/compare_x5_final/`
   (script of 8f1d0e1, run by the detached `validate/post_compare.sh` after the workflow ends).
