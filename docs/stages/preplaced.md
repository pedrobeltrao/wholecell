# Stage 6: preplaced structures

**Purpose.** Place large or connected structures before tile packing: chromatin, polysomes and mRNA, microtubules and
other cytoskeleton, nuclear pore complexes, lamina. Packing then treats them as fixed obstacles.

**Name.** `preplaced` everywhere: build folder `preplaced/`, `python -m wholecell preplaced`,
`workflow/rules/preplaced.smk` (rules `preplaced_<module>`, `preplaced_serial`, `preplaced_manifest`), code `wholecell/preplaced/`, registry
kinds `preplaced` / `preplaced_context`, config section `preplaced:`. Until A0c (2026-11-17) all of these were called
`prepass`. What still carries the old name, and why:
- builds made before A0c (the reference builds: Jurkat `v7repro2`, yeast `s0`, the testcell CI builds): `prepass/`,
  `pack/prepass_*.npy`. Readers find them through `LEGACY_NAMES` / `ref_path` in `wholecell/build.py`, the one map of
  old names; `scripts/ci_compare_builds.py` uses the same map to compare an old build with a new one;
- old configs: the keys listed in `wholecell/config.py` `LEGACY_CONFIG_KEYS` (`prepass:`, `pack.prepass_index`,
  `export.prepass_groups`, ..., region kind `prepass`) are renamed on load, with one message;
  `python -m wholecell prepass` runs the `preplaced` stage and prints "stage 'prepass' is now 'preplaced'";
- seed stream names (`prepass.<module>`, Jurkat `prepass.lamina`): the name is the seed, so it is kept;
- output content (not file names): JSON keys and values that readers and reference comparisons use, e.g.
  `pack/libs.json` key `prepass` (yeast: also `pack/summary.json`), `validate/placement.json` `n_prepass` and `clash_pairs`,
  `densities.json` `prepass_pg`, suite metric `prepass_md5_ok`, and the viewer format (tile_format_v2: region kind and
  group `prepass`, tile flag bit 2 `prepass`, box molecule `src` `prepass:<lib>`);
- the vendored yeast engine (`plugins/yeast/engine/`, its `yc3d/prepass.npz`, `wc3d.prepass`).

**Inputs.** Two kinds of builders, both registered as `preplaced` components:
- *module* builders (one independent layer, `--step module:<name>`; `frozen` = a declared input file). Jurkat:
  chromatin v8 rep5 (genome-mapped polymer), polysomes, traced microtubules, NPC (7R5K), mRNA and nascent RNA beads
  (all frozen), and the lamina (rebuilt and compared, frozen copy used).
- *serial* builders (`--step serial`, `wholecell/preplaced/serial.py`): an ordered list run in one process on one shared context
  (geometry, library pack, pools, obstacle sets). Each builder avoids what the ones before it placed and deducts its
  copies from the pools the next ones read, so the order is part of the result. Yeast (A6): the 28 builders of the
  vendored stage 3, in its order. Jurkat (B4): the common builders (below). Testcell: no modules.

**Outputs** (`preplaced/`): `preplaced_manifest.json` (what the pack index reads); per module builder
`<name>/selected.json`, `<name>/consumption.json`, `<name>/build/`, `<name>/compare_frozen.json`; serial run
`serial/preplaced.npz` (all arrays; yeast: same content as stage-3 `yc3d/prepass.npz`), `serial/layers/<name>.npz`
(the arrays each builder added or replaced), `serial/expected_after_preplaced.parquet` (pools after the pre-pass), `serial/lib_ids.npy` (type order of the
context's library pack), `serial/serial_manifest.json`. The pack stage indexes them into `pack/preplaced_{types,pos,quat,off}.npy`,
`preplaced_libs.json` (+ `preplaced_groups.json` for a grouped source).

**Config keys** (`preplaced:`): `numba_threads`, `modules` (ordered list; each with `name`, `kind` `instances|beads`,
`mode` `frozen|build`, `npz`/`json`, and for builders `builder`, `seed_stream`, `params`, `inputs`, `frozen`,
`rebuild_check`, `drops_target`), `drops`; `serial` (`context`: registered `preplaced_context`, `context_params`, `out`,
`builders`: ordered `{name, builder, params}`; yeast params = stream numbers and stoichiometries). Spec:
[../reference/cell_config.md](../reference/cell_config.md) §3.4 (draft).

**Checks.** Frozen files verified by checksum; builder vs frozen comparison (Jurkat B5, lamina; not identical, so the
frozen lamina is used); Jurkat B7 (pre-pass index identical, v7 port). Yeast serial run vs the port build s0: S1–S4 + pools
(`scripts/preplaced_a6_check.py`). B4: `scripts/b4_preplaced_check.py measure|verdict` (pores, MT dimers, polysome
fraction, ribosome NN, pre-pass self-overlaps by module pair, mito machines). Validation step `chromatin` reports frozen chromatin metrics.

**Code.** `wholecell/preplaced/{stage.py, serial.py, builders/}`; Jurkat-specific pieces in `plugins/jurkat`; yeast
serial builders and context in `plugins/yeast/preplaced/`.

**Common builders (B4, 2026-11-17).** The mechanisms both cells need are ONE core module each in
`wholecell/preplaced/builders/` (numerics ported from the vendored yeast modules, `# origin` lines; no engine module
calls): `npc_pores` (one scaffold per pore of the geometry package), `microtubules` (`core.microtubules`, track source
`procedural` = nucleation from an MTOC (yeast SPB) or `traced` = an EM polyline layer, Jurkat; `core.microtubules_mtoc` /
`core.microtubule_pools` keep the yeast names), `polysomes` (yeast walk with the 80S mRNA channel, single-stranded
linkers, collided disomes, membrane tethers; transcripts from a per-gene table (yeast) or a frozen per-mRNA table
(Jurkat); **point priors**: EM ribosome centres placed where measured and seeding the chains), `mrna_threads`
(threads through the chains + free mRNAs; `nt_per_bead` 1 (yeast) or 5 (Jurkat, the v7 representation)),
`mrnp_ends`, `npc_receptors`, `fg_mesh`, `atp_synthase_dimers` (rows on crista rims: yeast prims or the lamellar
cristae of a core package), `mtdna_nucleoids`. Geometry tests inside their numba kernels come from the cell's
evaluator through `ctx.K` (`wholecell/preplaced/kern.py`): the yeast X5 layout's own evaluator (vendored cell3d /
pack3d, so the yeast pre-pass stays bitwise) or the core evaluator `cellgeom_eval_v1` (Jurkat, new cells). Cell data
the builders read (yeast gene table, endocytic exclusion zones, OM/ER tethers; Jurkat mRNA table) are registered
`preplaced_data` providers (`plugins/yeast/preplaced/polysome_data.py`, `plugins/jurkat/preplaced.py`). Criteria:
`docs/unification_criteria/B4_preplaced.md`.

**Core context** (`core.cellgeom`, `wholecell/preplaced/context.py`): geometry from the build's package, library pack
of `context_params.libs` (library sphere trees), pools = pack-table copies, frozen instance layers as obstacles
(`context_params.obstacles`, selected by a key), named seed streams. Its `exports` write module files
(`preplaced/serial/modules/<module>.npz` + `.json` with `consumed`) that the manifest lists as `mode: serial`
modules, so the pack index and the consumption bookkeeping are unchanged; `frozen_subset` exports select rows of a
frozen file at build time (e.g. the centrioles and cortical actin of the v7 cytoskeleton, the nuclear mRNAs of v1).

**Jurkat vs yeast today.** Yeast serial builders (A6 order, unchanged): core `npc_pores`, `microtubules_mtoc`; yeast
`kinetochores`, `spb`, `rabl_chromatin`, `spindle_pools`; core `npc_receptors`, `fg_mesh`; yeast `rdna_transcription`,
`rdna_chromatin`, `pol2`, `origins`, `chromatin_sites`; core `microtubule_pools`; yeast `eisosomes`, `mcc`,
`actin_cables`, `endocytic_coat`, `amphiphysin`, `actin_patches`; core `atp_synthase_dimers`, `mtdna_nucleoids`;
yeast `cvt`, `nvj`; core `polysomes`; yeast `rac`; core `mrna_threads`, `mrnp_ends`. Every builder opens its own
generator from its yeast stream number. The yeast library-side models (MT lattice, ATP dimer entry, RAC, FG anchors,
the 80S channel and disome caches) are still made by the vendored modules called from the yeast context.
Jurkat (B4): geometry `em_cell1_v8p` (2,696 pores = 11 NPC/µm² × 245.1 µm², checked against the Nup anchor;
`toolkits/em/jpores.py`); one serial pass on the core context: `npc_pores` (7R5K), `microtubules` (traced),
`atp_synthase_dimers`, `mtdna_nucleoids`, `polysomes` (EM point priors, polysome fraction 0.70 from this cell's EM;
free 80S packed), `mrna_threads`; frozen: chromatin rep5, lamina (frozen copy, builder compared), nascent RNA, the
centrioles + cortical F-actin of the v7 cytoskeleton, the nuclear mRNAs of v1. FG mesh, pore receptors and mRNP ends
are not enabled for Jurkat (NPC shells stay packer regions). `config/cells/jurkat_v7_port.yaml` keeps every frozen
v7 layer (no serial pass, 1,000 pores). Jurkat chromatin stays the frozen rep5 and is re-checked against each geometry.

**Methods text (draft).** Large assemblies whose positions are constrained by data or by connectivity were placed
before the remaining proteins. These include the chromatin fibre (Jurkat: a genome-mapped polymer model; yeast: a
Rabl-configuration fibre with a measured nucleosome map), microtubules (Jurkat: traced from the EM volume; yeast:
nucleated from the spindle pole body), nuclear pore complexes at the pore positions of the geometry, polysomes and
mRNA. Each layer was checked against the geometry it was built on.
