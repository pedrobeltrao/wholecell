# Stage 9: export

**Purpose.** Write the files that other tools read: the viewer registry entry, region names and colours, streamed
tiles for the web viewer, and simulation boxes for applications.

**Today.** Folder and CLI stage `export`. Steps (`wholecell/export/stage.py`): `regions`, `cell_entry`,
`sim_manifest`, and `tiles_v2` when `export.tiles.enabled`; `all` runs them in that order.

**Inputs.** Pack outputs, the library, the geometry, the validation manifest. Yeast also writes the legacy YCT1
encoding (plugin step `yeast_s5_encode`, vendored `s5_encode.py`) for the X5 comparison and yeast validation scripts.

**Outputs** (`export/`): `regions.json`; `cell.json` (this build's entry for the viewer's `cells.json`); `web/`
(YCT2 tiles, `regions.json`, `version.json`, `MANIFEST.sha256`, encode checks); `boxes/<name>/` (simulation box,
schema `wholecell.simbox/1`; monomer molecules carry the per-protein SEC-MS bound-fraction attribute of
`inventory/secms_bound_fraction.tsv` as metadata columns, box.json `secms_bound_fraction` = provenance + intended use,
[../methods_secms_bound_fraction.md](../methods_secms_bound_fraction.md)).

**Config keys** (`export:`): `cell_type`, `state`, `status`, `data_base`, `licence`, `default_bookmark`, `steps`,
`tiles` (`enabled`, `tile_nm`, `block_tiles`, `q_org_nm`, `q_range_nm`, `shard_max_mb`), `preplaced_groups`, `regions`,
`sim_defaults` (`temperature_K`, `viscosity_Pa_s`, `ionic_strength_M`). Specs:
[../reference/tile_format_v2.md](../reference/tile_format_v2.md),
[../reference/viewer_multi_cell.md](../reference/viewer_multi_cell.md),
[../reference/simulation_workflows.md](../reference/simulation_workflows.md) §3.

**Checks.** Encode round trip (`web/encode_check.json`); manifest of output hashes.

**Code.** `wholecell/export/{stage.py, yct2.py, quat.py}`; `workflow/rules/export.smk`.

**Jurkat vs yeast today.** Testcell and yeast export YCT2 tiles; for Jurkat the full YCT2 encode is an optional target
(`tiles.enabled: false`). Yeast reaches the core exporter through an adapter (`plugins/yeast/pack_view.py`). The
simulation-box exporter is a one-box stub. **A8** makes YCT2 the only format for both cells and adds the full simbox
exporter. The unified viewer (`viewer_multi_cell.md`) is separate work.

**Methods text (draft).** Each build was exported as streamed tiles (positions in nm, orientations as scalar-first
quaternions, region identifiers via a palette) with a checksum manifest, so that the viewer and simulations read the
exact build they cite.
