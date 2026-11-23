# Stage 5: structures (library)

**Purpose.** Provide one 3D model per species: residue beads plus a sphere tree for collision tests, with
transmembrane species oriented in the membrane.

**Today.** Folder and CLI stage `library`, two modes. `frozen` (testcell; the Jurkat reference preset
`jurkat_v7_port`): import a frozen library, unchanged. `build` (yeast policy, A3, 2026-11-17; Jurkat since B3,
2026-11-17): a scripted, canonical library written to a new folder per build
(`library/lib/`, built in `lib.partial/`, renamed when complete and made read-only; never edited in place) from
declared frozen start entries plus scripted builders, and derived packing shapes. Yeast settings:
`library:` in `config/cells/yeast.yaml` (A9; the A3 fragment `library_x5.yaml` was merged into it).

**Inputs.** Jurkat: species table written by the stage from this build's pack table (`library.species: {from:
pack_table, accs: ...}` -> `library/species.tsv`; rebuild v9, 2026-11-23): one row per species the pack will request,
structure source = the inventory's `struct_src` (the current coverage rule). Until v8 it was
`plugins/jurkat/data/library/species_v7.tsv.gz` (16,279 rows: the legacy table the v7 library was built from, with the
v7 npz source for 321 complexes; kept in the repo as history): 117 pack-table species of v8 (72 monomers, 45
homo-oligomers, 0.51 pg) had no entry and were dropped by the pack without a record. TM-route ids, start entries (`start_v7.tsv`: the
pre-pass templates `F_ACTB`, `MT_TUBA1B-TUBB4B` + json, md5-checked against the v7 manifest), and frozen inputs:
AFDB v6 human tar (+ index, fragment index), AFDB complex models (manifest in the repo), UniProt topology / table /
REST cache, SIFTS, the RCSB assembly cache (read offline; md5 of every file read in `library/work/`). The frozen v7
library `lib_structures_v7` (16,246 files, md5 manifest in the repo) is the reference preset's library. Testcell: `lib_testcell` (45 npz copied from the Jurkat v7 library, in `tests/data/testcell/lib`).

**Outputs** (`library/`): `library_manifest.tsv` (species, path, md5; build mode adds `origin` = the step that wrote
the file), `library_meta.json`; build mode also `lib/`, `library_af.json`, `pack_shapes/`, `species_ledger.tsv`
(every species: `status` built | fallback | decomposed | missing, `source` pdb | afdb_complex | afdb_monomer |
afdb_fragment | pdb_chain | placeholder; the stage fails before sealing when `missing` species carry more than
`library.ledger.max_missing_mass_frac` of the requested mass, default 1e-4), `decomposed.tsv` (assemblies split into
subunit monomers; the pack returns their copies), `species.tsv` (when written from the pack table).

**Fallback structures** (step `fallback`, core, rebuild v9): a species still without a model after the builder steps
gets one by declared tiers, never silently none. Proteins: (1) `pdb_chain`: the chain of a cached RCSB assembly of the
same protein that maps (SIFTS) and shows the most residues, at least `min_cov` (0.5) of the UniProt length (soluble
proteins only; option `pdb_chain`); (2) `placeholder`: one bead per residue (UniProt length, else mass / 110 Da),
soluble proteins as a ball of the protein volume (1.21 nm3/kDa), transmembrane proteins with each UniProt TRANSMEM
segment as an ideal helix across the membrane (axes on a hexagonal lattice 1 nm apart, midplane z = 0) and the soluble
stretches as balls outside |z| >= 2.3 nm on their side (cytoplasmic side +z from TOPO_DOM; N-terminus cytoplasmic when
unannotated). Assemblies (option `decompose`): split into their subunit monomers (AlphaFold F1 by random access into the
AFDB tar, else the tiers above). Fallback npz carry `src` = `fallback:...` and `fallback` = the tier. Jurkat: after
`tm_route`, before `adaptive_trees`; yeast: last step (placeholders only, no SIFTS input).

**Build steps** (one function each, `wholecell/library/build.py:STEPS`, run in order): `start` (copy the frozen start
entries, md5-checked) → the builders of `library.steps`: `complexes` (curated complexes of the species table from
their RCSB assembly; membrane orientation method per complex from a TSV: `ring_axis | tm_slab | long_axis_cring`),
`af_models` (one AlphaFold DB model per protein, pLDDT < 50 left out, TM proteins oriented from UniProt topology;
writes `library_af.json`), `pdb_oligomers` (curated PDB homo-oligomers, TSV with reasons), `pdb_rna` (reference RNA:
tRNA-Phe 1EHZ) → `seal` → `pack_shapes`. Inputs read in place, offline (a missing PDB assembly is an error).
Jurkat steps (B3, `wholecell/library/jurkat.py`; the five legacy scripts that wrote v7 in place, final code, once each):
`jurkat_assemblies` (every non-protein species from `pdb:` / `afdb_complex:`; TM assemblies oriented by
`membrane_orient`, option `tm_graft`: AlphaFold grafting of missing TM segments; NPC ring axis) → `jurkat_monomers`
(AFDB F1 per protein, TM monomers oriented) → `reorient_tm` (1-2-pass TM monomers: soluble blocks out of the slab) →
`tm_route` (membrane libraries of the 129 TM-routed species) → `adaptive_trees` (per 30 → 15 / 8 / 5 until the tree
volume is ≤ 2.5× the protein volume, membrane species included; > 20,000 beads skipped). Logs in `library/work/`.
The species table (`ing_id, kind, uniprot, pdb`, ingredient order) comes from the inventory; yeast: the plugin step
`yeast_library_species` (`before: [library]`: vendored `ingredients.build_ingredients` on the build's protein table,
as `s0_library.py`).

**Packing shapes** (§2b P8): `library/pack_shapes/pack_trees_15.npz` (ids, sig = md5 of the library file, off, C, R;
15 residues per sphere, yeast `PACK_PER`; option `max_beads`: larger species keep their library tree, Jurkat 20,000
for the NPC) and `pack_meta.npz` (ids, bounding radius, structure source, PDB id,
membrane flag, bead count), separate from the library's own trees (30 / 40 residues per sphere).

**Config keys** (`library:`): `mode` (`frozen` | `build`); build: `start`, `species`, `pdb_dir`, `afdb_tar`,
`topology`, `steps`, `pack_shapes` (`per`, `ids`, `extra_dirs`), `options`, `nproc` (spec in `stage.py`); frozen: `frozen` (a frozen input with `path`,
`manifest`, `md5`), `copy` (false = reference the frozen folder, verified file by file; DECISION 2026-11-15).
Spec: [../reference/cell_config.md](../reference/cell_config.md) §3.3 (draft).

**Checks.** Every manifest file exists with the same md5, and the folder holds no extra file; otherwise the stage
fails.

**Code.** `wholecell/library/stage.py` (modes), `build.py` (build steps), `builders.py` (beads, sphere trees,
orientation; ported from yeast `structlib.py` / `build_structs.py` @ e3b7448), `pack_shapes.py`; check
`scripts/dev/library_a3_check.py` (vs the yeast port build s0); tests `tests/test_library.py`; Jurkat builders in `wholecell/library/legacy/` (`jlib_*`, `lib_refine_st`,
`tm_graft`, ...); `workflow/rules/library.smk`.

**Jurkat vs yeast today.** Yeast builds its library with the core stage, mode `build` (A9; policy of `s0_library.py`:
a canonical library plus separate 15-residue packing shapes, 342 frozen start entries: Boltz-2 tRNAs, builder
templates). Later yeast steps still write entries the vendored code builds (assembly structures, pre-placed models)
into `engine/lib`, a writable copy; the pack's shapes are the pre-placed context's `engine/work/yc3d/pack_trees_15.npz`
over the packed types. **A3** moved this policy into the core (Jurkat keeps the frozen v7 library); results: see the A3 commit message
(`_dev_unify/a3/run_*/check_report.json`). **B3** (2026-11-17) rebuilds the Jurkat library with the
scripted stage (TM routing, adaptive sphere trees); the Jurkat pack keeps the library's own adaptive trees
(`pack.shapes.source: library`) until B1. Check: `scripts/dev/library_b3_check.py` (criteria
[../unification_criteria/B3_library.md](../unification_criteria/B3_library.md)).

**Methods text (draft).** Each species was represented by a coarse-grained model with one bead per residue, derived
from experimental structures or AlphaFold models, and a hierarchy of bounding spheres for fast overlap tests.
Transmembrane proteins were oriented so that their hydrophobic belt lies in the bilayer. The library is a versioned,
checksummed set of files.
