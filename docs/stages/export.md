# Stage 9: export

**Purpose.** Write the files that other tools read: the viewer registry entry, region names and colours, streamed
tiles for the web viewer, and simulation boxes for applications.

**Today.** Folder and CLI stage `export`. Steps (`wholecell/export/stage.py`): `regions`, `cell_entry`,
`sim_manifest`, and `tiles_v2` (alias `web`) when `export.tiles.enabled`; `all` runs them in that order.

**Web tiles (YCT2).** Each 200 nm tile carries PALETTE + INST (the packed and pre-placed instances), and since
2026-11-24 also (`wholecell/export/meshpoly.py`):
- MESH: the lipid membranes of the build's geometry package (membrane SDFs of `cellgeom_eval_v1`; membranes whose
  `leaflets` is bilayer/monolayer), marching cubes at `mesh_h_nm` (8 nm) on the tile's own grid, wound so that face
  normals point to positive SDF (cytoplasmic / N side), clipped at an open face of the imaged volume (Jurkat). The viewer
  draws procedural lipids on these triangles.
- POLY: polymer threads of the pre-placed stage in pieces of up to 64 points, each in the tile of its middle point:
  mRNA (polysome mRNA, mRNA threads, Jurkat nuclear mRNA), DNA (linker DNA between consecutive nucleosomes, Hermite
  curves as the yeast engine drew them; rDNA, DNA through Pol II, mtDNA), extra with a kind (FG repeats, pre-rRNA,
  nascent RNA). Threads stored at 3 nm (Jurkat) are resampled to 0.6 nm.
- GHOST: neighbour instances of the lipid-mask regions (membrane regions + `export.tiles.lipid_mask_regions`) that
  reach the tile, so lipids are masked by proteins across tile borders.

PALETTE and INST stay byte-identical to an export without these sections (ghost-only palette entries are appended).
Blocks are encoded in parallel (`workers`, default 8; output independent of the worker count). Standalone re-export of
an existing build, without touching the build folder:

    python -m wholecell export --build <build_dir> --step web --out <dir> [workers=8] [mesh=0] [polylines=0] [ghosts=0]

(`--cell` defaults to the build's `config.resolved.json`; `<dir>` becomes the web folder). Check against the build's
own export and geometry: `scripts/dev/check_web_sections.py --build B --web DIR`.

**Inputs.** Pack outputs, the library, the geometry, the validation manifest. Yeast also writes the legacy YCT1
encoding (plugin step `yeast_s5_encode`, vendored `s5_encode.py`) for the X5 comparison and yeast validation scripts.

**Outputs** (`export/`): `regions.json`; `cell.json` (this build's entry for the viewer's `cells.json`); `web/`
(YCT2 tiles, `regions.json`, `version.json`, `MANIFEST.sha256`, encode checks); `boxes/<name>/` (simulation box,
schema `wholecell.simbox/1`; monomer molecules carry the per-protein SEC-MS bound-fraction attribute of
`inventory/secms_bound_fraction.tsv` as metadata columns, box.json `secms_bound_fraction` = provenance + intended use,
[../methods_secms_bound_fraction.md](../methods_secms_bound_fraction.md)).

**Config keys** (`export:`): `cell_type`, `state`, `status`, `data_base`, `licence`, `default_bookmark`, `steps`,
`tiles` (`enabled`, `tile_nm`, `block_tiles`, `q_org_nm`, `q_range_nm`, `shard_max_mb`, `mesh`, `mesh_h_nm`,
`mesh_skip`, `polylines`, `ghosts`, `lipid_mask_regions`, `workers`), `polylines` (`linker_dna`, `skip`),
`preplaced_groups`, `regions`,
`sim_defaults` (`temperature_K`, `viscosity_Pa_s`, `ionic_strength_M`). Specs:
[../reference/tile_format_v2.md](../reference/tile_format_v2.md),
[../reference/viewer_multi_cell.md](../reference/viewer_multi_cell.md),
[../reference/simulation_workflows.md](../reference/simulation_workflows.md) §3.

**Checks.** Encode round trip (`web/encode_check.json`, with section totals: mesh tiles, vertices, triangles, mesh
area per membrane, ghosts, polyline pieces / points per source); manifest of output hashes; unit test
`tests/test_export_meshpoly.py` (MESH / POLY / GHOST round trip, unchanged INST, determinism).

**Code.** `wholecell/export/{stage.py, yct2.py, meshpoly.py, quat.py}`; `workflow/rules/export.smk`.

**Jurkat vs yeast today.** Testcell and yeast export YCT2 tiles; for Jurkat the full YCT2 encode is an optional target
(`tiles.enabled: false`). Yeast reaches the core exporter through an adapter (`plugins/yeast/pack_view.py`). The
simulation-box exporter is a one-box stub. **A8** makes YCT2 the only format for both cells and adds the full simbox
exporter. The unified viewer (`viewer_multi_cell.md`) is separate work.

**Methods text (draft).** Each build was exported as streamed tiles (positions in nm, orientations as scalar-first
quaternions, region identifiers via a palette) with a checksum manifest, so that the viewer and simulations read the
exact build they cite.
