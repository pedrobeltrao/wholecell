# Unification plan (phase 2): one set of procedures for every cell

Status: design v1, 2026-11-16. Design only: nothing here has been built or run. Author: design agent for Pedro Beltrao.
Basis: DIRECTION 2026-11-16 in `PLANNING_NOTES.md` ("unify procedures; re-run all paper simulations on the unified
builds"), the two reproductions (Jurkat v7 **REPRODUCED**, 69/72 exact; yeast X5 port build
`yeast/20261116-aabe7b5-s0-x5repro` **running, verdict pending**), the audits in `data/derived/repro_audit/{human,yeast}/`,
`PORTING_STATUS.md`, `plugins/yeast/steps.yaml`, `plugins/yeast/engine/PORT_PATCHES.md` and `docs/specs/*`.

Abbreviations: **H** = Jurkat (human) pipeline as ported (`wholecell/` core + `plugins/jurkat`); **Y** = yeast X5
as vendored (`plugins/yeast/engine`, git archive e3b7448); **legacy** = Jurkat v7 / yeast X5 builds; **unified** =
builds made by the procedures proposed here.

---

## 0. Ground rules this plan follows

1. **One common procedure per stage.** A cell differs only in its *data*: config values, geometry package, source
   tables and curated lists. A cell-specific *plugin* is allowed only when the biology has no counterpart in the
   other cell, e.g. Rabl chromatin vs a genome-mapped polymer, the cell wall, or Jurkat virtual regions. Each
   plugin is listed with a reason (§1).
2. **Legacy reproductions are done once.** v7 is reproduced; X5 is pending. Unified builds are compared with
   legacy on **summary statistics** and with a **no-regression rule** on the common validation suite (§5). They are
   not compared on bit-identity.
3. **The paper uses one build per cell:** Jurkat v8 and the next yeast build. Every paper simulation is re-run on
   these builds as a workflow (§6). Legacy results stay as an archival record and are never deleted (provenance
   policy).
4. **Simplicity beats generality** wherever they conflict (Pedro, 2026-11-16; §8). Each such choice is marked
   **[simplicity]**.
5. **Licence.** Everything is CC BY-NC 4.0. Code derived from yeast (`plugins/yeast/engine`, and anything ported
   out of it into the core) stays **local** until the yeast licence is confirmed in the yeast repository. Every
   ported function keeps an `# origin: yeast ... @ e3b7448` line.

---

## 1. Stage-by-stage plan

Effort is in agent-days: one focused session, including its regression check. **No output change** means the
stage must reproduce the legacy artefact exactly, or at the X5 level T. **Changes outputs** means it is accepted on
summary statistics (§5).

### 1a. Summary table

| Stage | Jurkat (H) today | Yeast (Y) today | Proposed COMMON procedure | Stays a cell plugin (why) | Trade-off / expected effect on outputs | Effort |
|---|---|---|---|---|---|---|
| **Geometry generation** | EM toolkit (`toolkits/em`) → cellgeom_v1 `em_cell1_v7`; engine `sdf_sections_v1` (Jurkat cell3d with constants in G) | analytic toolkit (`toolkits/analytic/yeast.py`) → cellgeom_v1 package; engine `yeast_cell3d_v1` (vendored cell3d: prims + real-shape bricks) | **One generic evaluator engine `cellgeom_eval_v1`**: brick SDF fields + prim table + **table-driven label rule** + membrane table with per-membrane half-thickness (Pogozheva). It is the same evaluator in numba, TS and WGSL. Toolkits stay separate because their inputs differ (EM volume vs analytic or mesh). | Toolkits: EM (Jurkat, any EM cell) and analytic (yeast, synthetic cells). These are *generators*, not evaluators. The yeast wall, bud scar and eisosome furrows become prims plus core labels with extensions ≥ 1000. | ULP-level differences only. Acceptance as B1/G1: labels ≥ 99.99 % identical, SDF \|Δ\| < 0.01 nm. **v8 changes outputs on purpose**: `ne_mid` fix (perinuclear 2.7 → ~8.3 µm³), NPC count/density. Yeast adopts per-membrane half-thickness (was a uniform 2 nm), so its membrane bands move by ≤ 1–2 nm. | 5 |
| **Copy numbers** | core `inventory/copies` (Geiger iBAQ + RNA imputation, 100 pg) | the same core step reproduces build W/X5 copies (I3: 2.6e-16) | **Already common.** Add the declared slots: SEC-MS abundance cross-check, surface-proteomics PM copies (Jurkat v8). | None: sources and absent lists are config data. | None until new data are switched on. | 1 (surface adapter) |
| **Localisation** | v0 HPA/UniProt/Itzhak, then rules v4–v7. v5 P and v6 R1 read earlier packs. | UniProt/GO with Huh precedence → Vögtle/Morgenstern → CYCLoPs → NUC_K | **Evidence-weighting framework** (§3): source adapters write one long evidence table, then a common combiner, a common N/C model calibrated on imaging, and per-cell override tables with reasons. **No pack-dependent rules.** | Source adapters and curated override tables per cell. They are data, not code. | Jurkat changes: v8 imaging-fitted fractions, and v5 P / v6 R1 re-expressed without a pack. A yeast preset reproduces X5 weights exactly; the common defaults then change yeast slightly (§3d). | 6 |
| **Assembly state** | core assembly adapter + SEC-MS toolkit (Fossati 2021, 25 anchors); Complex Portal greedy; AFDB dimers (ipSAE / SEC r / STRING) | engine `ingredients/assemblies/heterodimers/...`: curated complexes, pre-ribosomes, CPX with RCSB ≥ 0.8 coverage, AFDB heterodimers with proportional allocation, PDB homo-oligomers + audit | **One ordered assembly step list in the core:** curated machines → Complex Portal with a structure → SEC-MS assembled fractions (toolkit, already common) → AFDB hetero/homo dimers → homo-oligomers of the residual pool. **Allocation logic from Y** (proportional heterodimer allocation, PDB audit). **Evidence from H** (SEC assembled fractions). | Curated machine lists (80S, NPC, nucleosomes, pre-ribosome stages for yeast) are config tables. The pre-ribosome builder is a yeast plugin until a human equivalent is wanted. | Jurkat: small shifts in dimer selection. Yeast: SEC-MS enters at the next build (anchor set still to review with G. Magnusson). A legacy preset per cell gives exact reproduction first. | 5 |
| **Structure library** | frozen v7 library (md5 manifest). Builders ported, not run: 5 in-place steps (jlib_build → rebuild_tm → refine_st ×3 → reorient_tm → tm_route). | `s0_library.py`: scripted, canonical library; separate 15-residue packing shapes per build; 342 frozen start entries (Boltz tRNAs, templates) | **Y policy:** a scripted library stage writing a **canonical library** (one new folder per build, never edited in place) plus **derived packing shapes**. H additions as options: TM grafting/routing (`tm_graft`, `jlib_tm_route`) and adaptive sphere trees (the ≤ 2.5× volume rule). | Curated PDB ids per organism (config). Frozen deposited entries (Boltz tRNAs, builder templates, and AFDB dimer models fetched on a date) are declared frozen inputs. | Jurkat library rebuilt: sphere trees and TM orientations differ slightly from the frozen v7 library. Accept on per-species volume, R_h and TM-orientation checks (§5). Yeast must reproduce X5 L2 at level T. | 4 |
| **Pre-placed: chromatin** | frozen v8 rep5 genome-mapped polymer (OpenMM, Hi-C/EM restrained, 28 M nucleosomes) | Rabl compact fibre + measured nucleosome map (in s3) | **Common output contract and validation, two builders.** Both write one `chromatin` layer (nucleosome poses, fibre polyline, genome coordinates, territory labels) and both run the same checks: P(s), nucleosome count vs genome/NRL, radial profile, and the EM/Hi-C test where data exist. | `jurkat.chromatin_genome_polymer` (frozen deposit, GPU, not bitwise) and `yeast.rabl_chromatin` with rDNA. **Why:** the two genomes are organised differently (Rabl vs territories) and the methods follow the data available. | None now. Jurkat chromatin stays as the frozen rep5 (DECISION 2026-11-15: full refit deferred), but must be re-checked against the v8 geometry (`ne_mid` fix): `depends_on_geometry` check, then trim or re-fit if nucleosomes fall outside the new nucleus. | 3 |
| **Pre-placed: cytoskeleton** | EM-traced MTs (U-Net, frozen), centrioles, cortical F-actin band 8–150 nm, lamina lattice | MTs from the SPB (16 KMT / 4 nuclear / 3 cytoplasmic), actin cables + patches, Arp2/3, capping | **One MT builder** (`microtubules`: nucleation from an MTOC layer, lattice of tubulin dimers) with two track sources: `procedural` (Y) or `traced` (an EM polyline layer, H). **One actin builder** for filaments along tracks or in a cortical band. | Yeast actin patches and endocytic coat (`actin_patches`, `amphiphysin`). The Jurkat vimentin cage is a virtual region (IF only in lymphocytes). | Jurkat: MTs re-instantiated from the traced polylines; same tracks, possibly different dimer phase. Accept on dimer count ± 1 %. Yeast: unchanged (level T). | 4 |
| **Pre-placed: NPC** | 7R5K in each of 1,000 pores (300 EM + 700 random); NPC channel / FG as packer shells | 7N9F + FG mesh + pore receptors in 93 pores | **One NPC builder:** pores come from the geometry package (`pores.npz`, source EM\|random); scaffold structure per organism (config); FG mesh + receptors (from Y) as an option. | Nothing beyond the structure id and the Nup gene list. | **Jurkat v8 changes:** the NPC density decision (3.3× the EM count; Nup anchors imply ~2.2× more pores than 1,000 → pore count from anchors or EM, Pedro). Moving FG from a shell to an FG-mesh pre-pass changes the sh_NPC_fg region (placed 0.771). | 3 |
| **Pre-placed: ribosomes / polysomes / mRNA** | EM ribo_seg (584 k) + EM-seeded extensions to 2.1 M; 250 k mRNA threads; nascent RNA | `polysome`, `mrna_genes`, `mrna_localization`, `rac`, `mrnp`, ER ribosomes (Sec61/Ssh1), mitoribosomes | **One polysome/mRNA builder** (Y code base: gene-specific mRNAs, polysome fraction, ER and mito ribosomes) with an **EM point-prior option** (`point_priors` layer; H) for seeding. | None: the EM seed positions are a geometry layer, not code. | Jurkat changes: polysome chains from the Y builder seeded at the EM positions. Ribosome NN median must stay within ± 1 nm of v7 (27.05 nm EM, 28.24 model). Yeast unchanged. | 5 |
| **Pre-placed: membranes and shells** | lamina (rods), cortical band, shells as packer regions | eisosomes/MCC, cortical ER (geometry), shells as packer bands, NVJ | **Shells are packer regions** (config `pack.regions`, kind `shell`; §2). Membrane-bound lattices (lamina, eisosome Pil1/Lsp1) are pre-pass builders of one type, "lattice on a membrane patch". | Lamina (Jurkat), eisosomes/MCC, NVJ (yeast). Each is a registered `prepass` module with its parameters in config. | None. | 2 |
| **Pre-placed: cell-type extras** | Golgi ribbon, speckles, open face, TCR pools (inventory) | SPB/kinetochores, rDNA, Pol II/origins/sites, Cvt, mtDNA, ATP-synthase dimers, cell wall | Generic mechanisms move to the core: ATP-synthase rows on crista rims, supercomplexes, mitoribosomes, mtDNA nucleoids. **Both cells have mitochondria**, so Jurkat gains these builders. | SPB/kinetochores, Cvt, rDNA trees, silent chromatin, Pol II sites (yeast); open face, Golgi ribbon (Jurkat geometry). | **Jurkat changes:** ATP-synthase dimers on cristae, an explicit mtDNA nucleoid layer. Small mass moves in the matrix and IM. | 4 |
| **Packing** | core v7 packer (500 nm tiles, 8 phases, single visit, shells, virtual regions, thin-lumen budget, open face, pre-pass index for 10⁹ instances) | `yeast_rsa_v1` (200 nm, 6 size-class rounds, squeeze, pilot thinning, shell cap, organisation layer, exclusion zones) | **ONE packer = the core packer, with the Y algorithms ported in** (§2). Tile size, rounds, squeeze scope, thinning, organisation and exclusion zones are config switches with common defaults. | Virtual-region providers (Jurkat), the organisation inputs (STRING taxon; AF3 pool optional and unpublished). These are data or plugins, not code paths in the packer. | Both cells change at summary-statistic level. Expected gains: Jurkat placed fraction 0.964 → ≥ 0.98 (rounds + thinning), seam check on 500 nm tiles. Yeast: 0.992 kept; clash rate should fall (H acceptance kernel). | 10 |
| **Validation suite** | core: densities, crowding, loc_compare, placement integrity, ribosome NN, frozen chromatin metrics, MANIFEST; Tier-2 imaging not ported | `validate_all.sh`: ~25 scripts (localisation, image N/C, proximity, Hi-C, cryo-ET, diffusion, anchors, placement, condensate test, ...) + build-stamp check | **One suite** (§4) with one metric definition per check, per-cell **reference adapters**, and a `no_regression` comparison against a declared reference build. | Reference datasets (HPA/OpenCell vs SWAT/CYCLoPs/hyperLOPIT, etc.) are adapters. Organism-only biology tests (yeast tRNA–nucleolus, silent chromatin; Jurkat Tier-2) stay as extra checks, labelled cell-specific. | No build output changes. Some legacy numbers change because defects are fixed (v7 densities with the v6 mRNA file; the yeast clash definition). | 8 |
| **Export (tiles)** | YCT2 encoder (testcell, 1 block of v7); legacy YCT1 shards | YCT2 export ported+run (42 regions via palette); YCT1 for comparison | **YCT2 only** (`tile_format_v2.md`): 1 µm shards + index, u16 regions via palette, w-first quaternions, grid from the geometry. | None. | No content change; legacy YCT1 readers kept for archived builds. | 2 |
| **Viewer** | Jurkat fork of the yeast WebGPU app (Jurkat WGSL sections) | yeast app (hand-mirrored cell3d in TS/WGSL, TGRID = 26) | **One app** per `viewer_multi_cell.md`: `cells.json` selection screen, table-driven geometry evaluator (the same label rule as §1 geometry), `features.json` toggles. | Procedural fillers (glucan/mannan, polyP, lipid tables) as per-label data; bookmarks per cell. | Visual parity check against both legacy viewers (screenshots of every bookmark). Yeast-derived viewer code stays local until the licence. | 8 |
| **Simulation export** | 7 per-application extractors, pickled radius libraries | one manifest schema (mmCIF + JSON) cut from the viewer data | **`wholecell.simbox/1`** (`simulation_workflows.md` §3), one exporter reading level-P tiles + cellgeom walls; temperature, viscosity and ionic strength from config. | None. | Boxes re-cut from the new builds. Legacy boxes are converted once for `compare_legacy`. | 3 |

The table sums to ≈ 73 agent-days. Scheduled in §5, the total is ≈ **77 agent-days**:

| Part | Agent-days |
|---|---|
| A0–A9 (no output change, incl. docs and getting started) | ≈ 48 |
| B1–B6 (change outputs, incl. the packer benchmark) | ≈ 20 |
| Viewer | ≈ 8 (in parallel) |
| Surface adapter | ≈ 1 |

Parallel agents can shorten the wall time to roughly 4–6 weeks.

### 1b. What disappears

- **`pipeline_engine: yeast_x5`** and the vendored engine. They are retired once every row of
  `plugins/yeast/steps.yaml` says `core` or `plugin`. The vendored tree is kept read-only under
  `plugins/yeast/engine/` as the X5 port record (provenance).
- **`pack.engine` choice.** One packer means there is nothing to select **[simplicity]**. `yeast_rsa_v1` stays
  registered only to re-run the X5 reproduction.
- **Geometry engine selection** by the env var `WHOLECELL_GEOM_ENGINE`. One evaluator reads the package
  **[simplicity]**. `sdf_sections_v1` and `yeast_cell3d_v1` stay only for the legacy-reproduction tests.
- **Hand-mirrored geometry** in TS/WGSL: replaced by the label-rule table. A golden test samples 10⁵ points and
  compares numba with TS (WGSL in headless Chromium where available).

---

## 2. The packer

### 2a. Comparison

Sources: packer agent report, 2026-11-16 (code-level read); `v7repro2/validate/placement.json`; X5
`work/results/wholecell3d_validation.json`; `runs/X5_seed0/work/s3_pack.log`; and the X5 port pack log (running).

| Aspect | H v7 packer (core) | Y seeded RSA (`yeast_rsa_v1`) |
|---|---|---|
| Kernel | numba RSA, sphere trees on a 6 nm grid, tol 0.3 nm, exact SDF region tests | same kernel lineage (common ancestor); extra arguments for per-region budgets, flat alignment, bands, exclusion zones |
| Tiles | 500 nm, halo 40 nm, 7,485 tiles | 200 nm, halo 36.7 nm, 5,912 tiles |
| Order | 8 checkerboard phases, **one visit per tile**; within a tile: membranes → shells/virtual → volumes, largest first | **6 size-class rounds** × 8 parity passes, order reversed on alternate rounds; heavy tiles first |
| Requests | multinomial split (tile × library × region) from tile measures; regions from config | `distribute` with organisation log-modulation, stochastic rounding; regions **hard-coded** (`tiles3d.REGION_NAMES`, `REG_MASK`, `BANDS`) |
| Determinism | named seed streams. **Shown:** the port (14 workers) reproduced v7 (22 workers) bit for bit on all checked artefacts; testcell identical with 8 and 4 workers. | master seed `k + 1000003·SEED`. **Shown:** rebuild of build R identical with 56 vs 64 workers. Not identical across CPU generations (X5 port: ULPs at stage 2, so packing is level T). |
| Speed | v7: ~3.3 h on 22 workers (legacy); 4.1 h on 14 workers (port, eu-g4). 880 M instances: **~12 M instances per core-hour**. | X5: 54 min on 40 workers for 25.2 M placed; rounds 0–4 took 370 s, **squeeze 2,862 s (89 %)**. ~4 M instances per core-hour without the squeeze, 0.7 M with it. |
| Membrane insertion | TM on tile meshes, +z along the gradient | same |
| Membrane shells | ✔ 19 shells (side, d0..d1, +radius, labels, pore limit) | ◐ 7 cytosolic bands (a shell with d0 = 0); shell cap ✔ |
| Virtual regions | ✔ (mt_lattice, mt_plus, pcm, vim_cage; plugin providers) | ✘ |
| Thin-lumen budgets | ✔ ×10 tries / ×16 give-up | ◐ per-region arrays exist, all 0 |
| Open face | ✔ face test, mesh clip, copy scaling, face-band drops | ✘ (whole cell) |
| Sharding / partial reruns | ✔ `tiles=`, `phases=`, `box=` subsets (viewer sharding is export) | ✘ |
| Organisation layer | ✘ | ✔ STRING PageRank field (128 components) + AF3 pool partners |
| Partner proximity | ◐ kernel supports it, v7 passes empty arrays | ✔ 4.59 M of 8.75 M eligible placed next to a partner |
| Pre-placed obstacles | separate pre-pass stage; **per-tile bucketed index** (51 M chromatin + 90.5 M mRNA beads) | built serially inside s3; one global array masked per tile (does not scale to 10⁸) |
| Squeeze / thinning | ✘ / ✘ | ✔ rigid-body relaxation with exact RSA re-acceptance / ✔ pilot-calibrated (nucleolus ×0.859) |
| Exclusion zones (type-specific) | ✘ | ✔ (ribosomes near endocytic sites) |
| Packing shapes | refined adaptive trees written into the library | separate 15-residue packing trees per build |
| Placed fraction | 0.964 by mass (0.983 by count). Worst: sh_NPC_fg 0.771, NE 0.802, sh_NPC 0.81, ER lumen 0.824 | 0.9924 of the thinned request (0.9895 of the inventory). Worst: nucleolus 0.917 of request, 0.787 of expected before thinning |
| Clash / containment | sphere-tree overlap > 0.5 nm on 300 tiles: **0 / 0** per 10⁶; orientation sign errors 0.36 per 10⁶. The pre-pass overlaps with itself (2.18 M pairs), which is a pre-pass problem. | bead overlap < 0.3 nm in 64 cubes: **2.55e-5** (1.0e-5 without actin), above its own 1e-5 limit; beads outside region 6.9e-5 |
| Cell-specificity | regions, region map, shells, thin lumens, virtual providers, budgets, seeds and open face are all config | region ids, bands, `FLAT_MEM`, `PRE_GROUPS`, taxon 4932 and the builders are in code |

The clash numbers above use two different definitions. A common metric is part of §4.

### 2b. Recommendation

**One packer: the core (H) packer, with the Y algorithms ported into it.** This reverses the yeast audit's
"yeast packer" verdict. The audit's reasons for preferring Y were its *algorithms*: rounds, squeeze, thinning,
determinism and the master seed. All of these port as contained changes into a kernel that shares the same ancestor.
The reverse port, adding H features into Y, would mean rewriting Y's hard-coded region tables, adding a per-tile
obstacle index for 10⁹ instances, and adding open-face, virtual-region and subset re-runs. In short:

- the H packer is already config-driven and plugin-based;
- it is proven at 10⁹ instances and reproduces v7 exactly;
- it does not depend on yeast code, which matters while the yeast licence is pending (§0.5).

**Features to port into the core packer.** Each is a config switch with a cell-independent default.

| # | Feature (from Y) | Config | Common default | Note |
|---|---|---|---|---|
| P1 | size-class rounds with alternating pass order | `pack.rounds: {n: 6, split: equal_volume, alternate: true}` | 6 | Fixes tile seams (yeast face density 0.49 vs 0.20 before). Seams in v7's 500 nm tiles were never measured: add a seam check (§4). Expected cost ×1–1.5 of a single visit (370 s for 25 M in Y). |
| P2 | pilot-calibrated thinning of over-full regions | `pack.thinning: {pilot_tiles: 16, enabled: true}` | on | Makes "requested = placeable" explicit and logged. Jurkat nucleolus, NPC shells and ER lumen are candidates. |
| P3 | shell cap (shell volume fraction ≤ cytosol's) | `pack.shell_cap: true` | on | |
| P4 | per-region try budgets (generalises thin_try_f) | `pack.budgets.per_region: {...}` | Y global 2000/1200 vs H 200/25 + ×10/×16: **benchmark first** (P8) | |
| P5 | squeeze relaxation with exact re-acceptance | `pack.squeeze: {scope: underfilled, threshold: 0.97}` | **targeted**: only (tile, region) pairs below 97 % placed after rounds | **[simplicity + cost]** A full squeeze is 89 % of Y pack time. At Jurkat scale (35× instances) a full squeeze extrapolates to ~28 h on 40 workers. Targeted squeeze keeps most of the gain. |
| P6 | organisation layer (STRING PPR field + partner CSR) | `pack.organisation: {string_taxon, pool_scores: optional}` | on; AF3 pool **off** in published builds unless deposited | The kernel already has partner attraction. Jurkat gains STRING-driven co-localisation. Output change: local neighbourhoods only, not compartment totals. |
| P7 | type-specific exclusion zones | `pack.exclusion: [{types, around: <prepass group>, r_nm}]` | none | Yeast endocytic sites; Jurkat could use it at the synapse face. |
| P8 | packing shapes separate from the library | library stage writes `pack_shapes/` (15 residues per sphere, adaptive refinement option) | on | Matches the library policy of §1. |
| P9 | master seed for everything | `seeds.master`; legacy streams only in `seeds.legacy` | on | Already core. |

**Kept from H as-is:** shells, virtual regions, open face, thin-lumen flat alignment, the pre-pass index, phases,
subset re-runs, and config regions.

**Tile size** is config: default 500 nm, yeast 200 nm. A common value is not required for correctness. The
benchmark (P8 below) decides whether to set 300 nm for both **[simplicity: one number in the docs]**.

**Benchmark before freezing defaults** (≈ 1 agent-day, CPU only): on a 500-tile Jurkat subset and on the whole yeast
cell, measure placed fraction, worst-region fraction, seam index, clash rate and wall time for:

1. v7 settings;
2. + rounds;
3. + rounds + thinning;
4. + targeted squeeze;
5. + full squeeze (yeast only).

Choose the cheapest setting that reaches placed ≥ 0.98 by mass in every volume region with ≥ 1 µm³ and seam
index ≤ 1.1.

**Acceptance of the unified packer:**

1. With the switches at legacy values (`rounds: 1`, no thinning, no squeeze, no organisation), Jurkat v7 is
   reproduced bit for bit (criteria B6–B11 level E). This is cheap and protects every H feature.
2. With yeast settings, the X5 summary statistics are met (P5/P6 tolerances of the X5 criteria).
3. With the common defaults, both cells pass §5b.

---

## 3. Localisation: a common evidence-weighting framework

### 3a. Shape

```
source adapters (per cell, data only)          common combiner                     per-cell overrides
  annotation:  UniProt CC, GO CC, HPA text  ─┐
  imaging:     HPA-Jurkat, OpenCell,         ├─> evidence.tsv ──> weights by class ──> nc_model ──> curated overrides ──> species_loc.tsv
               CYCLoPs, SWAT, Huh            │   (long form)       + precedence          (calibrated)   (reason + source per row)
  fractionation: Itzhak, Vögtle, Morgenstern │
  surface/proteomics: surface MS (v8)       ─┘
```

**`evidence.tsv`** (one row per protein × source × compartment) has these columns:

- `uniprot`, `source`, `evidence_class` (one of `imaging_cell` / `imaging_other` / `annotation_curated` /
  `annotation_ht` / `fractionation` / `proteomics_surface`);
- `compartment` (core vocabulary of `config/validate/loc_vocab_*.yaml`, mapped by each adapter);
- `score` (0–1);
- `reliability` (e.g. HPA Enhanced 1 / Supported 0.8 / Approved 0.5; CYCLoPs LOC-score);
- `n_obs`, `note`.

**Combiner (core, plain Python, one function per step):**

1. **Precedence by class** (R8d made general): same-cell imaging > curated annotation > other-cell imaging >
   high-throughput annotation > fractionation (fractionation is a tie-breaker only). Within a class, scores are
   reliability-weighted means.
   - Yeast: "keep Huh-supported terms only" becomes Huh as `imaging_other` outranking `annotation_ht`.
   - Jurkat: the 50:50 HPA/UniProt blend becomes a configurable `blend` between two classes.
2. **Multi-localisation:** compartments keep their relative scores. The split rule is a config choice: `rank`
   (Y: 1, ½, ¼) or `equal` (H). **Common default: `score`-proportional; equal when scores tie** **[simplicity: one
   rule]**.
3. **New compartment from imaging only** if score ≥ `new_min` (Y 0.5). Never nucleocytoplasm-only.
4. **N/C model** (`nc_model`, common): within the nucleus + cytosol pool, the concentration ratio is
   `r = f(imaging evidence)`. f is fitted per cell on the calibration data declared in config:
   - Y: `r = exp(k·(N − C LOC-score))`, k = 5 fitted on the SWAT dev half;
   - H v8: R8b category → ratio calibration on HPA-Jurkat and OpenCell images, replacing v6 R1 constants
     N 10 / Nc 3 / NC 1 / Cn 0.3 / C 0.02.

   Shares then follow from `fn = r·Vn / (r·Vn + Vc)`, using the **current build's** compartment volumes (geometry
   measures). This replaces v5 P, which raised the nuclear share to the volume ratio when the v4 pack placed none.
   Without imaging evidence the default is equal concentration (r = 1, same as Y), unless the annotation is
   nucleus-only or cytosol-only.
5. **Membrane moves need a mechanism.** v6 rule kept as a common rule: a move into a membrane or shell key needs a
   domain/GO feature or a second imaging source.
6. **Overrides** (per cell, a TSV with `uniprot, compartment shares, reason, citation, rule_id`): Y `LOC_FIX`, H v4
   curated list, v7 structural pools (TCR/CD3, PCM, MT bands, vim_cage, NPC receptors, lamins), and yeast cell-wall /
   prenylation exclusions. These are **data**. The rule code that applies them is one function.

### 3b. What stops being pack-dependent

- **v5 P and v6 R1** gate on the *placed* shares of the v4/v5 packs. In the framework they gate on the
  *intended* shares of the current inventory, plus a capacity check from the pilot thinning of §2 (P2). Nothing
  reads an earlier pack.
- The frozen tables `v4pack_loc_compare`, `v5pack_loc_compare` and `tier2_v5` stay declared frozen inputs **only**
  for the legacy-v7 reproduction preset.

### 3c. v8 plans folded in

- **Imaging-fitted fractions (R8a).** The 111 HPA-Jurkat-imaged genes get direct per-gene fractions: per-cell NNLS
  unmixing, probe-calibrated to 3D, then shrunk to the category prior. They enter as `imaging_cell` evidence with
  `score` = fraction.
  - Only nucleus, nucleolus, cytosol, ER+Golgi and PM are identifiable. Other compartments keep annotation
    evidence.
  - Yeast uses the same mechanism for `image_nc` (6,145 C-SWAT strains). Yeast then validates image N/C on a
    held-out split (`splits.py`, hash-based, build-independent), because its NUC_K is fitted on SWAT.
- **The 30 % HPA nuclear floor.** This is an *observation*, not a rule: cytoplasmic proteins show ~0.3 nuclear
  signal in HPA-Jurkat images (0.27 for MS-detected genes; nuclear proteins ~0.75). The framework gives every imaging
  source an explicit **background term** `b`: observed nuclear fraction = b + (1 − b)·true. `b` is estimated, not
  assumed, by the R8c tests:
  1. quantify ~300 OpenCell live HEK293T images with the same unmixing code. Live GFP has no antibody background,
     so ≈ 0.05 means the HPA floor is background and ≈ 0.3 means a real nuclear pool;
  2. use genes not expressed in Jurkat (RNA = 0) imaged by HPA as a direct background basis;
  3. compare antibodies for genes with ≥ 2 HPA antibodies.

  **Decision rule (proposal):**
  - If (1) gives ≤ 0.10 and (2) agrees within 0.05, set `b_HPA` = the median nuclear fraction of non-expressed
    genes, and do not apply it to OpenCell.
  - If (1) gives ≥ 0.2, treat it as a real pool: `b = 0`, and the nuclear shares of "cytoplasmic" proteins rise.
    This is a large mass move, so Pedro reviews it (> 2 pg rule).
  - In between: `b` is a fitted parameter with its CI propagated; report both extremes as a sensitivity.
- **Acceptance (v8 plan):**
  - changed-gene image agreement ≥ v7;
  - cross-validated image agreement > v7 held-out;
  - OpenCell agreement ≥ v7 (0.831);
  - no single rule moves > 2 pg without Pedro's review;
  - report agreement across versions v1 … v8 (DECISION 2026-11-14), not a held-out half.

### 3d. Effect on yeast

A `legacy_x5` preset (rank split, Huh precedence, NUC_K = 5, CYCLoPs thresholds 0.05/0.5) must give X5's
`comp_weights` within 1e-6. The common defaults (score-proportional split, precedence classes) then move some
multi-localised yeast proteins. This is accepted if the yeast localisation benchmark stays at or above X5:

- SWAT top-1 0.665 and hyperLOPIT 0.848, each ± 0.01;
- image N/C Spearman 0.81 ± 0.01.

---

## 4. Common validation suite

**One suite in `wholecell/validate/`.** It has:

- one metric definition per check;
- reference **adapters** per cell;
- an output `validate/suite.json` with, per metric: value, reference, direction (higher/lower better), tolerance
  and status;
- `validate/no_regression.tsv` against the declared reference build.

Every check is run on every cell that has the reference data. "n/a (no data)" is a valid status, never silently
skipped.

| # | Check | Metric (common definition) | Jurkat reference | Yeast reference | Direction | Origin |
|---|---|---|---|---|---|---|
| C1 | Localisation vs independent references | agree, agree_tvd, top-1, mass-weighted (intended and placed) | OpenCell (0.831), HPA consensus (0.832), HPA-Jurkat | SWAT N-GFP (0.665 top-1), SWAT NOP1, hyperLOPIT (0.848) | ↑ | H loc_compare + Y localization.py |
| C2 | Localisation across versions | the same metrics for every archived version (v1…v8; yeast S…next) | as C1 | as C1 | report | H imaging_across_versions |
| C3 | Image N/C | per-gene Spearman of model vs image N/C; AUROC nuclear vs not | HPA-Jurkat Tier-2 (0.48–0.53) via NNLS unmixing | YeastRGB C-SWAT (0.81, AUROC 0.98) | ↑ | Y image_nc + H tier2 (port) |
| C4 | Image compartment shares | 1 − TVD per imaged gene (5 resolvable patterns) | 111 HPA-Jurkat genes (~0.53) | n/a or CYCLoPs-derived (circular, report only) | ↑ | H v8 imaging |
| C5 | Copy-number anchors | ratio of anchor proteins to structural need (histones vs nucleosomes, r-proteins vs ribosomes, Nups vs pores) | 108 / 71 / 45 pg-equivalent | 0.90× / 1.36× | → 1 | core copies anchor check + Y anchors.py |
| C6 | Copy numbers vs an independent dataset | Spearman | a second Jurkat proteome (if available) | Das 2026 (0.81) | ↑ | Y copy_numbers.py |
| C7 | Clashes | **two levels, one definition each:** (a) sphere-tree overlap > 0.5 nm per 10⁶ instances on ≥ 300 seeded tiles; (b) bead overlap > 0.3 nm per bead pair in ≥ 64 seeded 100 nm cubes. Packed–packed, packed–pre-pass and pre-pass–pre-pass reported separately. | 0 (a) | 2.55e-5 (b) | ↓ | H placement + Y s4_validate |
| C8 | Containment | instances whose centre label is not allowed for their region; beads outside the region | 0 | 6.9e-5 (beads) | ↓ | both |
| C9 | Membrane orientation | sign errors per 10⁶; median tilt | 0.36 | 0° median | ↓ | both |
| C10 | Placement integrity | placed / requested per (species, region); species < 50 % placed; seam index = density in the 10 nm face band / tile interior | 0.964 mass | 0.992 | ↑ / seam → 1 | Y placement_checks + H summary |
| C11 | Densities and crowding | mg/ml per volume region (packed + pre-pass); excluded-volume fraction; N/C density ratio | cytosol 162.8, nucleoplasm 135.1 | cytosol 170 | within literature range (QPI/SMR for Jurkat; yeast literature checks) | both |
| C12 | Ribosomes | NN distance distribution vs cryo-ET/EM; concentration | EM median 27.05 nm vs 28.24 | CZII NN 29.4 vs 28.5; Delarue 23 µM vs 17.1 | → ref | both |
| C13 | Chromatin | P(s) slope; nucleosome count vs genome/NRL; radial profile; EM voxel r (where EM exists); Hi-C contact metrics | v8 metrics (P(s) −1.29 vs −1.13 ± 0.15 **fail**) | WLC R(s), Hi-C slope −1.01, nucleosome occupancy | report + fixed thresholds | both |
| C14 | Diffusion benchmark | probe-level \|log10(D_model / D_meas)\| on a frozen test set; within-compartment τ; inert probes | new prereg set (§6.1) | 93 literature D, DOI split (RMSE 0.21–0.23) | ↓ | H BD v2 + Y diffusion.py (one method, §6.1) |
| C15 | Proximity / interaction | AUROC of model neighbour enrichment vs PPI sets | BioPlex/OpenCell IP (adapter to write) | BioGRID PCA/BioID/XL-MS (0.85/0.81/0.99), PPiSeq | ↑ | Y proximity |
| C16 | Build integrity | every result file is from this build (MANIFEST, 0 stale); `depends_on_geometry` of every layer = this package | both | both | must | Y build_stamp + core MANIFEST |
| C17 | Cell-specific extras (labelled) | yeast: tRNA–nucleolus, silent chromatin, DNase, kinase localisation, ER ribosome density, mRNA localisation, condensate negative control; Jurkat: fluorescence panels, galleries | — | — | report | as now |

**No-regression rule (used by §5):**

- A metric **regresses** if it moves in the bad direction by more than its tolerance relative to the reference
  build: localisation/AUROC/Spearman ±0.01; rates max(10 % rel., 1 per 10⁶); densities ±2 mg/ml; NN median
  ±0.5 nm; placed fraction ±0.003.
- A build with any regression is **not accepted** unless a dated DECISION in `PLANNING_NOTES.md` explains it. An
  example is a defect fix that removes a flattering artefact.
- Improvements are reported and need no decision.

Effort: ≈ 8 agent-days. The work is porting the ~25 yeast validation scripts behind reference adapters, porting the
Jurkat Tier-2/fluorescence code, the C7 two-level clash metric, and the seam index.

---

## 5. Sequencing and acceptance

**Reference builds.**

| Cell | Reference build | Status |
|---|---|---|
| Jurkat | `jurkat/20261115-ec0b403-s20261115-v7repro2` | v7, REPRODUCED |
| Yeast | `yeast/20261116-aabe7b5-s0-x5repro` | once its verdict is REPRODUCED (or WITH DEVIATIONS, explained) |

If the X5 port fails, phase 2 for yeast waits. Phase 2 for Jurkat does not.

**Rule for every step.**

- Criteria are frozen first, in `docs/unification_criteria/<step>.md` with an md5 in PLANNING_NOTES (the same
  discipline as v7/X5).
- The step lands as one commit plus one build per cell.
- The verdict goes to `validate/compare_<step>/verdict.json`.

### 5a. Steps with no output change (refactors; acceptance = legacy artefacts identical or at X5 level T)

| # | Step | Acceptance | Effort |
|---|---|---|---|
| A0 | Docs + getting started (§8) and the testcell CI on every commit | testcell full build < 10 min on 8 cores; determinism with 8 vs 4 workers | 2 |
| A1 | Generic geometry evaluator `cellgeom_eval_v1` (table-driven label rule, prims, membrane table) | Jurkat B1 level E; yeast G1–G5 level E vs the port build; golden numba ↔ TS test | 5 |
| A2 | Yeast regions as config (`pack.regions`, no `REG_MASK`/`BANDS` in code); yeast tiling (`tiles3d`) moved to the core tiling step with tile size from config | yeast G4/G5 and P1/P3 identical to the port build | 3 |
| A3 | Library stage in the core (Y `s0` policy; H builders as options); Jurkat keeps the **frozen** v7 library in this step | yeast L1–L5 at level T vs the port build; Jurkat manifest md5 unchanged | 4 |
| A4 | Assembly steps in the core with presets `legacy_v7` / `legacy_x5` | Jurkat B3/B4 md5; yeast I1/A1 identical | 5 |
| A5 | Localisation framework with presets `legacy_v7` (frozen pack tables declared) and `legacy_x5` | Jurkat pack_table md5 d1a0099c; yeast comp_weights ≤ 1e-6 | 6 |
| A6 | Pre-pass builders as registered modules: yeast builders moved out of `s3_pack` unchanged; Jurkat frozen layers unchanged | yeast S1–S4 level E/T vs the port build; Jurkat B7 level E | 6 |
| A7 | Packer: port P1–P9 behind switches, legacy defaults | Jurkat B6–B11 level E; yeast pack through the core packer with Y settings meets X5 P3/P5/P6 at level T | 6 |
| A8 | Validation suite (§4), simbox exporter, YCT2 for both | suite runs on both reference builds; all values reproduce the legacy numbers where defined (code check), and defects are explained | 10 |
| A9 | Retire `pipeline_engine: yeast_x5` for building (it stays for reproduction only) | yeast config has no `pipeline_engine` block; X5 port re-run through the core reproduces at level T | 1 |

Subtotal ≈ **48 agent-days**. Several steps run in parallel: A1 ∥ A3 ∥ A4/A5, and A8 alongside A6–A7.

### 5b. Steps that change outputs (acceptance = summary statistics + no regression)

Each step is applied to **both cells** unless marked otherwise, starting from the A9 state.

| # | Step | Summary-statistic acceptance (vs reference build) | Effort |
|---|---|---|---|
| B1 | **Packer common defaults** after the benchmark (rounds, thinning, targeted squeeze, shell cap, budgets, packing shapes) | placed mass fraction ≥ reference − 0.003 overall, and ≥ 0.98 in every volume region ≥ 1 µm³ (else explained); seam index ≤ 1.1; C7(a) ≤ 1 per 10⁶ packed–packed and C7(b) ≤ 1e-5; C8 ≤ reference; densities per region ± 3 mg/ml unless placed mass changed (then mass-consistent); pack wall time ≤ 2× reference on the same node | 4 |
| B2 | **Organisation layer on** (STRING; AF3 pool off) | compartment totals unchanged (per (species, region) within max(3√n, 1 %)); C15 proximity AUROC ≥ reference; no regression elsewhere | 2 |
| B3 | **Library rebuilt for Jurkat** (scripted, TM routing + adaptive trees options) | per-species volume within 2 %, R_h within 2 % for ≥ 99 % of species; TM orientation: hydrophobic belt centred within 0.5 nm; C9 ≤ reference | 2 |
| B4 | **Common pre-pass builders** (polysome/mRNA, MT, NPC, mito machines) | ribosome NN median ± 1 nm; polysome fraction within literature (Jurkat EM, yeast 0.82–0.85 target); MT dimer count ± 1 %; pore count = config; pre-pass self-overlaps (C7a pre-pass–pre-pass) ≤ reference | 6 |
| B5 | **Localisation common defaults** (no pack-dependent rules; score-proportional split; N/C model with the current calibration) | C1 ≥ reference − 0.005 per reference; C3 ≥ reference − 0.01; mass moved by any single rule ≤ 2 pg without review | 3 |
| B6 | **Geometry fixes** (Jurkat `ne_mid`, new geometry version; yeast per-membrane half-thickness) | cellgeom validator passes; perinuclear volume as measured (~8.3 µm³); all `depends_on_geometry` updated; chromatin nucleosomes inside the nucleus ≥ 99.9 %; no regression | 3 |

Subtotal ≈ **20 agent-days**. Each step changes outputs on its own, so every regression has one cause. They run in the
order B6 → B3 → B4 → B1 → B2 → B5: geometry first, because everything depends on it.

### 5c. Freezing the paper builds

**Jurkat v8 = A9 + B1–B6 + the following, in this order:**

1. **Logged fixes:**
   - PM topology (42 units, 207 k copies on the wrong side: TFRC, ITGB1/ITGA4, TRBC1, CD7 …);
   - NPC density / pore count;
   - SOS copies and the PDGFRB 64 k check;
   - ledger `rule_id` column;
   - `ne_mid` (in B6).
2. **SEC-MS** (new Jurkat dataset through `toolkits/secms`, assembled fractions + per-subunit free pools;
   `measured_free_pool` stays OFF unless the gates in `secms_toolkit_report.md` pass).
3. **Surface proteomics** (`adapter.surface`: PM abundance and topology). It cross-checks or replaces the PM
   topology fix.
4. **Imaging localisation** (§3c), after the R8c decision on the nuclear floor.
5. Optional: TCR microvillus-tip enrichment (Pedro).

**Yeast next build = A9 + B1–B6 + yeast SEC-MS** (anchor set reviewed with G. Magnusson), plus whatever yeast
model rounds after X5 G. Magnusson wants included. Each such round must arrive as data or config changes, or as a
plugin through the unified pipeline, never as a new vendored engine.

**Freeze criteria (both cells):**

- every 5b step accepted;
- the validation suite has no unexplained regression vs the reference build;
- **3 seeds** built (seed spread reported for every suite metric; the paper build is seed 0, the others are
  replicates);
- a clean commit (not dirty), with the tag `paper-<cell>-v8` / `paper-yeast-<id>`;
- `export/MANIFEST.sha256` recorded in PLANNING_NOTES;
- the build folder made read-only.

After the freeze, any change means a new build id, and the simulations would have to be re-run. So freeze only when
the new data (SEC-MS, surface) are in.

**Gating items outside our control:** the new SEC-MS data (not local yet, both cells), Jurkat surface proteomics, the
R8c decision, and the yeast licence (needed for the public repo, not for the build).

---

## 6. Simulation re-run plan

Each application becomes `applications/<name>/` per `docs/specs/simulation_workflows.md`:

- standard rules `check_build → boxes → setup → run → analyse → prereg_check → evaluate → report`;
- inputs from one pinned build (`build_id` plus the export manifest sha256);
- boxes from `wholecell.simbox/1`.

GPU numbers below come from `bd_gpu_timing.md`, `bd_results.md`, `bd_v2_results.md`, PLANNING_NOTES, the yeast
`sweeps/*.json` wall times and the endocytosis `out/release_*.json` step rates. The legacy GPU was a shared RTX
3090/4090, so estimates are for **one unshared GPU of that class**.

**Pre-registration status:**

- **Reuse** means the md5-frozen file still defines the test, because it does not name build quantities. A dated
  addendum states the new build id before the run.
- **Re-register** means the definitions contain build quantities (copies, module compositions, boxes, θ), so a new
  prereg is frozen on the new build **before** any run.

### 6.1 BD diffusion v2, both cells (`bd_diffusion`)

**The two methods today.**

| Aspect | Human | Yeast |
|---|---|---|
| Engine | CuPy Euler–Maruyama, rigid sphere-tree bodies | CuPy Euler–Maruyama |
| dt | 0.05 ns | 0.2 ns, drift cap 0.15 nm |
| Temperature | 310 K | 303.15 K |
| Contact potential | soft harmonic core | krep 20 |
| Hydrodynamics | mean-field Tokuyama–Oppenheim H(φ) | none / mean-field Cichocki / RPY Ewald |
| Interactions | cosine wells (STRING K_d, quinary ε_ns 2.45 kT, chromatin) | none |
| Boxes | 500 nm periodic + walled lumen boxes (2 nm SDF) | 100–200 nm periodic, 8 pooled cytosol boxes |
| Run length | 30–100 µs | 2–10 µs |
| Output | per-probe D/D0 vs a 213-row benchmark | per-body D vs R_h → diffusion law → 93 literature D |

**Common method (proposal).**

- **One engine `wholecell/sim/bd`:**
  - the human engine as base: wells, walled boxes, long runs, sphere-tree bodies from `bodies.npz`;
  - the yeast drift cap;
  - hydrodynamics option `meanfield` (default; **one formula**, chosen by the benchmark below) and `rpy` as a
    check on a subset **[simplicity: one default]**;
  - temperature, viscosity and ionic strength from the build's `export.sim_defaults`.
- **One protocol:**
  - boxes: 8 periodic boxes per volume compartment (500 nm for cytosol/nucleus where the cell allows, else the
    largest cube that is ≥ 99 % one compartment; yeast cytosol 300 nm), plus walled lumen boxes (ER, matrix);
  - interactions: quinary + membrane + chromatin terms; STRING only in the size-neutral form (DECISION proposal
    2026-11-14, BD audit);
  - observables: both outputs, per-probe D/D0 against the literature benchmark and the D vs R_h law (which the
    yeast validation and the vesicle sims need).
- **Benchmark to choose the hydrodynamic default** (≈ 1 GPU-h): run both engines on the same yeast cytosol box and the
  same Jurkat cytosol box, with none / meanfield-TO / meanfield-Cichocki / RPY. Pick the mean-field form closest to
  RPY.

**Pre-registration.**

| File | Status |
|---|---|
| BD v2 held-out `bd_holdout.tsv` (md5 409929f8…) | used, so legacy only |
| BD v3 (671895d3…) | failed; legacy |
| Yeast DOI/sha split (`splits.py`, md5 b261ff88…) | hash-based and build-independent: **reuse** for yeast |
| Human | **re-register:** a fresh test set (new literature probes, or the frozen BD v3 test set ebba238b… as test rows, with mapping rules M1–M5) |

θ is either fitted on the training rows in a `fit` rule, or declared frozen (v2 tier-0 `theta_final.json`
20dbd5d5…). The prereg states which one.

**Cost.**

| Item | GPU-h |
|---|---|
| Jurkat: 8 + 8 periodic + 6 walled boxes, tier 0–2 at the v2 lengths | ~10 (v2 whole ~7, lumen ~1, more boxes) |
| Yeast: pooled sweeps (0.65 GPU-h on S3) | ~3 at longer length |
| Yeast: Gade test (2.3) | ~3 |
| Hydrodynamics benchmark | ~1 |
| **Total** | **≈ 17** (2 seeds) |

**Parallel units:** box × compartment × seed. One job per box, so 22 + 16 jobs fit on any number of GPUs.

### 6.2 Synapse S0–S2.1 with CRISPR and CD43 (`synapse`, Jurkat only; CPU)

**Inputs from the build:**

- PM surface (cellgeom field/mesh crop);
- PM units with heights (`species_loc` + library);
- cortex F-actin layer;
- ledger copies (SOS, RasGRP1, RASA2, CD43 …);
- θ from 6.1.

**Pre-registration.**

| Arm | Status | Reason |
|---|---|---|
| S0 test definitions f84599c8… and split 90220889… | reuse + addendum | |
| S1 017b1ec7… | reuse + addendum | it names `pm_units_v7` heights, so the addendum declares the v8 units *before* running |
| S2.0 be8f17b6… | reuse + addendum | parameters must become a frozen `params` rule output (they were computed live from the ledger) |
| S2.1 `S21_freeze.md5` + held-out definitions 88edf312… | literature held-out H1–H5 reused | the freeze covers code and params, so it is **re-frozen** as `S21v8_freeze.md5` with v8 copies |
| CRISPR f6090676… | **re-register** | it names v7 quantities: RASA2 10.4 %, RASGRP1 695 |
| Glyco/CD43 8fb8dcac… | **re-register** | CD43 68,627 copies "added to v7" |

The test *rules* are copied verbatim; only the build quantities change.

**Cost:** S2.0 34.8 + rerun 35.7 + S2.1 27.8 + CRISPR ~16 + glyco ~10 ≈ **125 CPU-h**. No GPU. The natural units are
runs (geometry × seed × dose × perturbation), which suit a CPU array job (≤ 1 h each). With 64 cores that is about
2 h wall.

### 6.3 Condensate residue gate (`condensates`)

**Module definitions from the new inventory.** The six level-2 modules (SG run-off/rest/no-RNA, PSP rest/no-RNA, PB
rest; md5s c8b07ca0…, 08ecd404…, 22145b40…, fd6d777a…, 0815b9b3…, eb1e89a6…) are species lists with copies drawn
from the v7 inventory. They are **regenerated by a `modules` rule** from the v8 inventory with the *same selection
code and rules*, then frozen in a new prereg.

**Pre-registration.**

| Item | Status |
|---|---|
| In-vitro panel (17 systems) and residue protocol `test_definitions_residue_frozen.md` 3b08bb5c… | build-independent: **reuse**. The panel does not need re-running unless the force field changes; its legacy results stand. |
| Gates on modules | **re-register** with the new modules |
| Cell-scale feasibility `criteria.md` | frozen md5 56a086a0… but the file now hashes 7f088438… because DEVIATION 1 was appended. `prereg_check` must hash the part above the deviation line, or the deviation goes to a separate addendum file. Fix this at import. |

**Cost.**

| Item | GPU-h |
|---|---|
| Modules only (6 systems × 2 seeds × 2 starts at 0.1 µs; residue chain stages modules_* ≈ 15 h) | ~15 |
| + cell-scale continuum on the new geometry (1.4 GPU-h used) | ~2 |
| **Total** | **≈ 17** |

Optional, not planned: the longer SG runs (~6 GPU-days) and in-situ nucleolus (30–165 GPU-h), each needing Pedro's
go.

**Parallel units:** system × seed × start. One CUDA process per GPU (OpenMM); up to 24 independent jobs.

### 6.4 Yeast vesicle / endocytosis (`vesicle_traffic`)

**Inputs from the build:**

- the deepest endocytic site (or all sites above a depth): crop + simbox with ~538 k molecules;
- the build's cytosol diffusion law (from 6.1, replacing `calibrate.py` against an old build);
- Golgi cisternae and endosomes for transit.

**Pre-registration: none exists, so register before the re-run.** It should cover:

- results 1–7 of the README as tests: vesicle area 6,380 ± 1,929 nm², centroid distance ~81 nm, rvs pinch
  28–50 nm, turgor bracket;
- that the actin-flow agreement is imposed (stated as such);
- ≥ 3 seeds;
- the T3 knock-outs (no helix, no Bzz1, no Syp1, no crowding) as pre-declared arms.

The sub-project code is on the local branch `sim-endocytosis` only. It must be pushed to a remote and frozen at a
commit (G. Magnusson).

**Cost.**

| Item | GPU-h |
|---|---|
| Release v2: 12 seeds × ~1.1 h (13.5 GPU-h) | ~14 |
| Release v3: 4 seeds × 0.9 h (3.5) | ~4 |
| Three more sites × 3 seeds | ~10 |
| **Total** | **≈ 28** |

Membrane stage, transit and the T3 scans are CPU. The legacy runs got ~3 runs per GPU at ~600 steps/s each, so the
units can be packed 3 per GPU.

### 6.5 NFAT / calcium A-series (`nfat_calcium`, Jurkat)

**Inputs:**

- whole-cell geometry at 100 nm (from the cellgeom package, not `jgeom_cell1_v7.npz`);
- pores (`pores.npz`, v8 count);
- ledger concentrations;
- D from 6.1 θ;
- chromatin layer (A2).

**Pre-registration:** `test/test_plan.json` 42e4e8ad… with a sha1 gene split (dev GSE85201 / held-out GSE178827).
It is build-independent: **reuse + addendum** naming the v8 build. Pore permeability is re-fitted to Kar & Parekh
2015 in a `fit` rule, declared in the addendum.

**Cost:** A1 7 scenarios × 0.25–0.8 GPU-h ≈ **5 GPU-h**; A2 16 runs, small (< 2 GPU-h). Total **≈ 7 GPU-h**.
Parallel units: scenario (A1), TF × start × sensitivity (A2).

### 6.6 Totals

| Workflow | GPU-h | CPU-h | Prereg action |
|---|---|---|---|
| BD diffusion (both) | ~17 | small | human re-register; yeast reuse |
| Synapse + CRISPR + CD43 | 0 | ~125 | reuse + addenda; CRISPR and glyco re-register; S2.1 re-freeze |
| Condensates (modules + cell-scale) | ~17 | small | panel reuse; modules re-register; fix the cellscale md5 handling |
| Yeast vesicle / endocytosis | ~28 | ~20 | new prereg |
| NFAT / Ca A-series | ~7 | small | reuse + addendum |
| Yeast condensate negative control (validation C17) | 0 | < 1 | reuse (`calibration_nucleus.json` frozen) |
| **Total (2 seeds where stated)** | **≈ 70** (+ 30 % contingency ≈ 90) | ≈ 150 | |

On 4 GPU nodes this is ≈ 1 day of wall time; on one GPU ≈ 4 days. Chromatin (deferred), in-situ condensates and the
long SG runs are excluded.

### 6.7 GPU job kit: we write the scripts, Pedro submits

We never submit jobs. Each application ships a **submit kit**:

```
applications/<name>/slurm/
  SUBMIT_README.md        what to submit, in which order, expected GPU-h / memory / wall per job, how to check
  jobs.tsv                one row per job: job_id, rule, params (box, seed, arm…), gpu (0/1), mem_gb, time_h, depends_on
  submit_<batch>.sh       parameterised sbatch array script: reads its row of jobs.tsv by $SLURM_ARRAY_TASK_ID
  collect_verify.py       after jobs finish: checks every expected output exists, run.json (seed, code md5, env lock,
                          GPU, wall) is complete, no NaN/blow-up flags, outputs in the right run folder; writes
                          collect_report.tsv + a one-line status; then runs the downstream CPU rules (analyse,
                          prereg_check, evaluate, report) locally
```

- **`jobs.tsv` is generated by Snakemake** (`snakemake --until run --dry-run` turned into a job table by a small
  script). The SLURM path and the local path run the *same* rule commands.
- **`submit_<batch>.sh` contains the following and nothing else:**
  - `#SBATCH --array=1-N%K --gpus=1 --mem=<from jobs.tsv max> --time=<from jobs.tsv max>`;
  - placeholder lines `#SBATCH --partition=TBD --account=TBD`, which Pedro fills in;
  - `conda run -p <env> python -m wholecell.sim run --job $ROW`.

  There are no node names or absolute paths beyond `config/site.yaml`.
- **Batches** follow the dependency order: (1) BD hydrodynamics benchmark → (2) BD boxes (both cells) →
  (3) synapse / A-series / vesicle runs that need θ and the diffusion law; condensates run independently. The
  `SUBMIT_README.md` of each application lists its batches, and a top-level `applications/SUBMIT_ORDER.md` gives
  the cross-application order with GPU-h totals (§6.6).
- **Re-runs:** `collect_verify.py --failed` writes `jobs_failed.tsv` and the matching `submit_rerun.sh`. A failed
  job is re-submitted with the same seed once. A second failure is reported and never silently re-seeded.
- **No agent monitors the jobs** (MEMORY "token economy"). `collect_verify.py` is run by Pedro, or by an agent when
  asked, after the queue is empty.

---

## 7. Risks and decisions for Pedro

**Decisions needed:**

1. **Packer base = the core (H) packer with the Y algorithms ported in.** This reverses the yeast audit's verdict;
   §2b gives the reasons. Agree?
2. **Squeeze scope:** targeted (only under-filled regions) as the common default, to keep Jurkat packing under
   ~2× the v7 time; full squeeze only if the benchmark shows a clear gain.
3. **AF3 pool scores** (unpublished) **off** in published builds unless they can be deposited. The yeast
   organisation layer then uses STRING only, which differs slightly from X5.
4. **Nuclear-floor rule** (§3c): estimate `b` from OpenCell live images and non-expressed genes, then apply the
   decision rule. Approve the ~300-image OpenCell quantification (CPU).
5. **NPC count for Jurkat v8:** EM count, Nup-anchor count (~2.2× the 1,000 pores), or the literature density?
6. **v8 waits for the SEC-MS and surface data?** Proposal: run A* and B* now, and freeze only when the data arrive.
   Otherwise freeze without them and accept a v9 later, which means re-running every simulation.
7. **Yeast intake:** the next yeast build starts from X5 + unified pipeline + SEC-MS. Later yeast rounds enter only
   through the unified pipeline. This needs G. Magnusson's agreement.
8. **BD:** a human re-registration with a fresh test set (BD v3 test rows) and one common BD engine (human base,
   mean-field default chosen by benchmark).
9. **Re-registration policy:** an "addendum before run" for build-independent preregs, and new preregs for those
   that name v7 quantities (CRISPR, glyco, condensate modules, S2.1 re-freeze).

**Risks:**

- **The X5 port verdict is still pending.** If it fails, the yeast side of phase 2 waits for a diagnosis. Jurkat is
  unaffected.
- **Jurkat-scale cost of the yeast algorithms** is extrapolated (rounds roughly the same as v7, full squeeze tens of
  hours). The benchmark in §2b must come before the defaults are fixed.
- **Chromatin vs the new geometry:** the `ne_mid` fix and the B6 check may push nucleosomes of the frozen rep5
  outside the nucleus. A trim is cheap; a re-fit (20–35 GPU-h) was deferred by DECISION.
- **Yeast licence:** none yet. The core packer port (P1–P9) carries yeast-origin code, so the public repo waits for
  it. CC BY-NC is intended.
- **External data timing** (SEC-MS, surface proteomics) sets the freeze date.
- **Prereg hashing:** files with appended deviations (cellscale `criteria.md`) break whole-file md5 checks. Import
  them with a "frozen part + addenda" convention.

---

## 8. Ease of understanding (primary design goal)

Other people will use and adapt this pipeline, so where generality and simplicity conflict, we choose simplicity and
say so.

**Few concepts.** A cell is **config + geometry package + data adapters + plugins** (ARCHITECTURE §1); there is
nothing else to learn. Phase 2 *removes* concepts:

- `pipeline_engine`;
- `pack.engine`;
- geometry-engine selection;
- legacy rule versions as code (v4–v7 become one combiner plus override tables);
- vendored yeast stages.

**Stage names match the paper's Methods**, in order:

`geometry` → `composition` (copy numbers) → `localisation` → `assembly` → `structures` (library) →
`preplaced` (today "prepass") → `packing` → `validation` → `export`; applications → `simulations/<name>`.

- The renames ("prepass" → `preplaced`, "inventory" split into `composition` / `localisation` / `assembly`) are
  done once in A0. Old folder names are kept as read-only aliases for the reference builds only.
- Each stage has one module folder, one doc page and one Snakemake rule file with the same name.

**Readable configs.**

- One YAML per cell, sections in stage order, and a comment on every non-obvious value with unit and source
  (citation or DECISION date).
- Defaults live in `config/defaults.yaml`, so a cell file lists only what differs. Today `jurkat.yaml` repeats
  packer budgets that will become defaults.
- Curated lists are TSV files with a `reason` column, not YAML literals.

**One obvious way to add a cell.**

```
python -m wholecell new-cell mycell --from testcell     # copies config + a tiny analytic geometry
# edit config/cells/mycell.yaml (sources, geometry package, regions); add data adapters only if a format is new
snakemake --config cell=mycell                          # builds; validation tells you what is missing
```

`docs/adding_a_cell.md` walks through it with the testcell as the worked example. Plugins are needed only for
biology the core lacks, such as a new pre-placed structure.

**Plain Python over clever abstractions.**

- One function per step with explicit inputs and outputs.
- The registry is a dict of names to functions.
- No metaclasses, no code generation for numba backends (one evaluator reads tables), and no environment-variable
  switches for behaviour.
- numba only in the kernels (packing, geometry), each with a pure-Python reference used in tests.

**Good error messages.**

- Config validation names the key, the expected type and the allowed values, and suggests the nearest valid name
  ("region `cytosl` not in geometry labels; did you mean `cytosol`?").
- Missing frozen inputs print the declared path and md5 and say how to obtain the file.
- A stage that finds a stale or other-geometry layer stops and names the layer and both geometry ids.
- Every check failure in the suite prints the metric, value, reference, tolerance and the doc anchor.

**Getting started in minutes.** `docs/getting_started.md`:

1. create the env (one command);
2. `snakemake --config cell=testcell` (all stages, about 5 min on 8 cores today);
3. open the testcell in the viewer (local server);
4. read `validation/suite.json`.

CI runs this on every commit.

**Docs structure.**

```
docs/
  overview.md            what a build is, the 4 parts of a cell, the stage diagram (1 page)
  getting_started.md     testcell in minutes
  stages/<stage>.md      one page per stage: purpose, inputs, outputs, config keys, checks, paper Methods text
  adding_a_cell.md       the walk-through above
  simulations.md         application workflows, prereg, GPU submit kits
  reference/             specs (cellgeom_v1, cell_config, inventory_schema, tile_format_v2, simbox, viewer),
                         config key reference (generated from defaults.yaml comments), criteria files
```

The present `docs/specs/*` move under `reference/`. Reproduction criteria and this plan go under
`reference/history/`.

**Where we chose simplicity over generality:**

| Choice | Instead of |
|---|---|
| one packer | selectable engines |
| one geometry evaluator | engine plug-ins |
| one multi-localisation split rule | per-cell split rules |
| one hydrodynamics default | three options |
| targeted squeeze | full |
| tile size as one documented number per cell | adaptive tiles |
| no code generation for TS/WGSL | a table-driven evaluator (slower to write, simpler to read) |
