# Porting status (phase 0-1, 2026-11-15; yeast X5 port 2026-11-16)

Legend: **ported+run** = in the pipeline and exercised by a build; **ported, not run** = code copied/parameterised
but frozen outputs are used instead; **stub** = interface only; **todo** = not started.

## Done

| Area | Status | Notes |
|---|---|---|
| Foundation | ported+run | `wholecell/{config,build,seeds,registry,__main__}.py`: site-only absolute paths, build stamp, write-once stage folders, frozen-input md5 verification, named seed streams + declared legacy seeds, PYTHONHASHSEED=0 enforced |
| Workflow | ported+run | `Snakefile`, `workflow/rules/*.smk`; DAG generated from the cell config (rule sets, pre-pass modules, validation steps); `scripts/build_isolated.sh` = frozen git worktree per build, detached |
| Envs | done | `envs/vcell.explicit.txt` + pip freeze (all CPU steps), `envs/wholecell_smk.yaml` |
| Specs | done | `docs/reference/{cellgeom_v1,cell_config,inventory_schema,tile_format_v2,viewer_multi_cell,simulation_workflows}.md`; `docs/ARCHITECTURE.md` |
| Geometry | ported+run | cellgeom_v1 writer/loader/validator; engine `sdf_sections_v1` (cell3d with constants in G); Jurkat package `em_cell1_v7` converted from the v7 npz (B1: 2.55 M points, labels and SDFs bit-identical); testcell package `analytic_v1` |
| EM toolkit | converter run; generator ported, not run | `toolkits/em/`: segmentation, refinement, HMFF NE, jgeom_build, cristae, face cap, speckles, Golgi; geometry steps restructured to write package N+1 |
| Analytic toolkit | ported+run (testcell); yeast stub | `toolkits/analytic/generate.py`, `yeast_stub.py` |
| Inventory (Jurkat) | ported+run | rna → v0 (100 pg from config, no default) → assembly → ingredients → rules v4..v7 (new file per step, `changes.tsv` with consistent `rule_id`) → pack_table (byte-identical to `ingredients_v7.tsv`) → canonical tables incl. SEC-MS per protein/complex and empty surface-proteomics slot; pack-dependent rule inputs are declared frozen tables |
| Inventory (generic) | ported+run | table adapters, complex table assembly, signal-peptide lumen rule (testcell) |
| Library | built per build (B3, 2026-11-17); frozen v7 kept for the reference preset | Jurkat `library.mode: build`: core stage with the Jurkat steps `jurkat_assemblies` (option `tm_graft`) → `jurkat_monomers` → `reorient_tm` → `tm_route` → `adaptive_trees` (`wholecell/library/jurkat.py`, legacy workers in `wholecell/library/legacy/`; organism tables `plugins/jurkat/data/library/`; offline, declared inputs); 4 frozen start entries (pre-pass templates) from the v7 library; `jurkat_v7_port.yaml` keeps `mode: frozen` (v7 manifest). Criteria `docs/unification_criteria/B3_library.md`, check `scripts/library_b3_check.py` |
| Pre-placed | common builders (B4, 2026-11-17); Jurkat frozen layers partly replaced | one core module per common builder (`wholecell/preplaced/builders/`: NPC, microtubules procedural|traced, polysomes with EM point priors, mRNA threads, mRNP ends, pore receptors, FG mesh, ATP synthase dimers, mtDNA; numerics ported from the vendored yeast modules, geometry through the cell's kernels `wholecell/preplaced/kern.py`); core context `core.cellgeom` with module exports. Jurkat (`jurkat.yaml`): geometry em_cell1_v8p (2,696 pores), NPC / traced MTs / polysomes / mRNA / ATP synthase / mtDNA built per build; chromatin v8 rep5, lamina (frozen copy), nascent RNA, centrioles + cortical F-actin and nuclear mRNAs stay frozen; `jurkat_v7_port.yaml` keeps every frozen v7 layer. Yeast: same builders, bitwise vs port s0 + PORT PATCH 4 (check: `scripts/preplaced_a6_check.py`); criteria `docs/unification_criteria/B4_preplaced.md` |
| Pack | ported+run | pack3d + jwc + jpack_test + jvirtual; regions/region_map/budgets/drops in config; tiles/counts/index and phase-0 tiles bit-identical to v7 in dev tests |
| Validation | ported+run | densities (VAL_GEOM/VAL_POLY defect fixed), crowding (stable sort), loc_compare + hold-out, ribosome NN, NEW placement integrity, frozen chromatin metrics, MANIFEST, `compare_builds.py` (criteria B1-B12, V1-V8) |
| Export | partial | `regions.json`, `cell.json` (cells.json entry), YCT2 encoder (tile_format_v2; testcell, 1-block v7 test), simulation-box manifest stub |
| SEC-MS toolkit | ported+run (2026-11-16) | `toolkits/secms/` (human \| yeast): inputs OpenSWATH / peptide long / protein matrix (+SEM) / protein long; top-N quantification, native-mass anchor calibration (+ Complex Portal / standards cross-checks, LOO), assembled fractions k = 2/3/4/6 with calibration lo/hi and replicate SD, peaks, per-complex co-elution and per-subunit free pools (Complex Portal, hu.MAP), candidate-pair co-elution; inventory step `secms` (`wholecell/inventory/secms.py`), canonical tables §4a; Jurkat assembly reads it (preset `legacy_v7`: v7 assembly tables md5-identical; orphan `sec_assembly_thresholds.tsv` regenerated); option `complex_copies: measured_free_pool` (OFF for Jurkat; dry run in `docs/reference/history/secms_toolkit_report.md`); testcell synthetic SEC-MS + `tests/test_secms.py` |
| Copy-number step | ported+run (2026-11-16) | `wholecell/inventory/copies.py` step `copies` (+ proteome adapter `copies`): typed sources, priority + `set_by`, RNA imputation, absent lists, per-protein scaling, total protein fixed / anchors / source_sum, structure-count adjustments, anchor check, copy ledger, future slots (SEC-MS abundance, surface proteomics); Jurkat runs it (legacy settings; v0 and pack_table md5-identical, `docs/reference/history/copies_regression_criteria.md`); yeast method expressed in `config/cells/yeast.yaml` (stub) and reproduces build W copies read-only (`scripts/copies_yeast_check.py`); report `docs/reference/history/copy_numbers_report.md` |
| Yeast X5 port | ported+run (2026-11-16) | Pipeline engine `yeast_x5` (`plugins/yeast/`: vendored yeast code `engine/` = git archive e3b7448 + 3 documented patches `engine/PORT_PATCHES.md`; `runner.py` = build_isolated.sh + run_all.sh semantics per step), packer `yeast_rsa_v1` selected by `pack.engine` (option a), geometry = ONE cellgeom_v1 package per build from `toolkits/analytic/yeast.py` (engine `yeast_cell3d_v1`, `plugins/yeast/geom_engine.py`), core steps reused in the same DAG: `inventory.copies`, export `regions` / `tiles_v2` (YCT2) / `cell_entry` via the adapter `plugins/yeast/pack_view.py`; DAG from `pipeline_engine.steps` (`workflow/rules/engine.smk`); step classification `plugins/yeast/steps.yaml`; criteria `docs/reference/history/yeast_x5_reproduction_criteria.md`; comparison `plugins/yeast/repro/compare_x5.py`. Verdict: see "Yeast X5 reproduction" below |
| Yeast through the core (A9) | wired (2026-11-17) | `config/cells/yeast.yaml` has no `pipeline_engine`: core geometry (named package), copies, localisation + assembly presets `legacy_x5`, library mode build, serial pre-placed builders, core packer preset `legacy_x5`, common suite, YCT2; yeast-only steps are plugin steps (`plugin_steps`, `wholecell/plugin_steps.py`, `plugins/yeast/bridge.py`; every row of `plugins/yeast/steps.yaml` is core / toolkit / plugin). The X5 reproduction (engine `yeast_x5`, packer `yeast_rsa_v1`) is `config/cells/yeast_x5_port.yaml` only |
| Reproducibility fixes | done | absolute paths → config; inventory_v0 total protein; hash() seeds → wholecell.seeds; ledger rule_id column; in-place chains → per-step files; scratch_checks code → toolkits/em; numba cache outside the legacy tree |

## Stubbed / partial

- Yeast (A9): the stages run in the core, but several core and plugin steps still call vendored engine modules in the
  build's `engine/` file layout (assembly hooks, the yeast pre-placed context, yeast-only builders and the library-side
  models of the common builders (MT lattice, ATP dimer entry, RAC, FG anchors, 80S channel caches), the yeast validation
  scripts); `plugins/yeast/steps.yaml` says which. `toolkits/analytic/yeast_stub.py` keeps only the primitive-table constants.
- Export: YCT2 ghosts, mRNA polylines, meshes, overview level; full Jurkat YCT2 encode is an optional target (not run).
- Simulation-box export: one-box stub; no application consumes it yet.
- Provenance ledger across builds (`wholecell/provenance/`): only the per-build `inventory/version_changes.tsv`.
- Geometry engine for analytic primitives: the yeast package carries `prims.npz` (§6.2) + baked fields + native files and
  is evaluated by `yeast_cell3d_v1` (the vendored cell3d); a generic `analytic_prims` engine, the declarative
  `label_rule` for yeast (waived in the package) and `cellgeom.load` support for prim packages are not done.

## Remaining

1. **Yeast, after X5 (phase 2)**: move engine steps into core stages one at a time, each checked against the X5 port
   build (not against X5 itself, which differs from any rebuild on this node at the ULP level): library builder, assembly
   state, tiling, shared pre-pass builders (polysomes/mRNA, MTs, NPC/FG, ATP synthase), validation suite; yeast-only
   builders become registered `preplaced` plugins. **Packer unification** (one engine for both cells) needs: the yeast
   rounds by size class + checkerboard phases + squeeze + pilot thinning + shell caps + organisation partners, and the
   human shells / virtual regions / thin-lumen budget / open face / 500 nm tiles / chromatin pre-pass index; a common
   request table (tile × type × region from config regions, no hard-coded ids), one geometry evaluator interface
   (`mems`, `labels`, `in_pore` from the cellgeom engine), tile size from config, seeds from `wholecell.seeds` with
   legacy streams declared per cell; acceptance = Jurkat v7 criteria and yeast X5 criteria at level T (E impossible
   across the two legacy RNG/ordering schemes). Copied yeast code must stay local until the yeast licence exists.
2. **New data adapters**: SEC-MS toolkit done (2026-11-16; new Jurkat / yeast datasets: checklist in
   `docs/reference/history/secms_toolkit_report.md`; yeast anchor set to review); Jurkat surface proteomics (PM abundance + topology); then
   re-express the pack-dependent rules v5 P / v6 R1 without a previous pack.
3. **Simulation workflows** (`applications/<name>/`, spec `simulation_workflows.md`): BD, synapse S0-S2.1, CRISPR/glycocalyx,
   condensates, NFAT/Ca, vesicle traffic (yeast); GPU profile to be defined with Pedro.
4. **Viewer**: implement `viewer_multi_cell.md` (cells.json, selection screen, table-driven geometry evaluator) and YCT2 decoding.
5. **Geometry**: the v7 `ne_mid` field is fixed in the new package `jurkat/em_cell1_v8` (B6, 2026-11-17: `toolkits/em/jne_fix.py`,
   cause = a uint8 `~occ` in `jgeom_build.mesh_sdf`; perinuclear 2.74 -> 8.02 µm³); yeast has per-membrane half thickness
   (package `analytic_x5_s0_f5f3a5c723b2_memhalf1`). Frozen pre-placed layers declare `depends_on_geometry` (checked against
   v8: chromatin rep5 inside the nucleus, criteria `docs/unification_criteria/B6_geometry.md`). Still open: run the restructured
   EM toolkit end to end; rebuild the chromatin/MT layers on the declared geometry (refit deferred, DECISION 2026-11-15).
   B6 full builds: scripts in `{results_root}/{jurkat,yeast}/_dev_unify/b6/full_build.sh`.
6. Library builder: run through the core stage for Jurkat (B3) and compared with the frozen v7 library (pre-build check
   `{results_root}/jurkat/_dev_unify/b3/`; full build script `.../b3/full_build.sh`); still open: moving the library
   structure source to the inventory's struct_src (321 assemblies differ), Tier-2 imaging and fluorescence validation.
7. Figures (`jurkat/paper`) not ported; `mkcache.py` and the figure caches were lost with the old /tmp scratchpad.

## Yeast X5 reproduction (2026-11-16)

Build `data/derived/wholecell_builds/yeast/20261116-aabe7b5-s0-x5repro` (launched 18:21, 24 workers, 1 RTX 3090;
expected end ~22:00). Verdict file: `validate/compare_x5_final/verdict.json` + `criteria_results.tsv` (written by the
detached `validate/post_compare.sh` after the workflow; the workflow's own `validate/compare_x5/` uses the pre-fix script).

Interim (stages 0-2 and the pre-pass log, checked 19:15):
- E: SDF fields (G2), tile measures/labvol (G4), 5,912 mesh tiles (G5), cellgeom package validation (G3), s1 log (I4),
  library file set (L1), 6,768/6,770 rebuilt library entries (L2: 2 within tolerance), library_af (L3), 342 start entries (L4);
  every pre-pass builder line checked so far (MTs, chromatin 64,487/88 fragments, rDNA, eisosomes, polysomes 39,314 /
  182,329 / 3,489 collided, RAC, mRNA) identical to X5.
- T: geometry.json (≤ 1.5e-12 nm, CPU ULPs: eu-g4 Zen 2 vs the X5 node), ingredients (4e-14 log2), protein_features
  (≤ 1.7e-12 rel.), core copies step vs engine (2.6e-16; I3), assembly tables (A1, after the row-order fix 8f1d0e1).
- Not rebuildable from scratch (frozen inputs): 42 Boltz-2 tRNA species + builder templates (342 lib entries), the
  AF3 pool pair summary (no generator), the RCSB answers (X5 final cache), the bundle's sweeps start copy.
