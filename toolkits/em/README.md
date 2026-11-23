# toolkits/em: EM segmentation -> cellgeom_v1 (Jurkat, jrc_jurkat-1 cell 1)

Paths come from `config/site.yaml` (`_paths.py`): inputs via `inp()` (the toolkit work folder
`{geom_root}/_work/jurkat_em/` if a regenerated file exists there, else the read-only legacy `{data_root}/data/derived`),
outputs via `out()` (work folder, write-once). Geometry steps write packages: `pkgstep.open_pkg(N)` then `pkgstep.publish(N+1)`.
Run as modules from the repo root: `PYTHONHASHSEED=0 python -m toolkits.em.<step> ...`.

## Run in wholecell

| module | what | result |
|---|---|---|
| `convert_legacy.py` | wraps the published v7 geometry (`jgeom_cell1_v7.npz`, NE PLY, MT traces, centrosome, pores, refs) | `{geom_root}/jurkat/em_cell1_v7` (validated, B1 exact) |
| `b1_check.py` | criterion B1: engine vs legacy `cell3d` (also `tests/check_b1_geometry.py`) | in `em_cell1_v7/validation.json` |
| `jne_fix.py` | B6 (2026-11-17): package step em_cell1_v7 -> **em_cell1_v8**, `ne_mid` recomputed as the true signed distance to the HMFF mesh on both sides (`jgeom_build.mesh_sdf` with the uint8 fix); every other file identical | `{geom_root}/jurkat/em_cell1_v8` (validated, `pores` passes without waiver) |

## Ported, in v7 order (NOT run in wholecell)

| # | module (origin) | status | in -> out |
|---|---|---|---|
| 0 | `dl_labels.sh` (scratch_checks) | ported, not restructured | OpenOrganelle labels mirror (needs `WHOLECELL_DATA_ROOT`) |
| 1 | `count_instances_per_cell.py` (scratch_checks) | ported, not restructured | ribo_seg / np_seg instances -> `*_instances.parquet` |
| 2 | `cells_s3.py` (jurkat/em) | ported, not restructured | 32 nm cell instances `cells_s3.npz` |
| 3 | `assign_instances.py` (jurkat/em) | ported; writes `*_instances_assigned.parquet` (legacy rewrote in place) | |
| 4 | `fetch_raw_cell1.py` (scratch_checks) | ported, not restructured | raw s1 crop zarr |
| 5 | `cell1_refine_s2.py` (scratch_checks; **the copy that made v7's `cell1_s2.npz`**) | ported, not restructured | `cell1_s2.npz` |
| 6 | `cell_refine_s2.py` (jurkat/em; run for cell 10 = Golgi template) | ported, not restructured | `cell10_s2.npz` |
| 7 | `build_labels.py` (jurkat/em) | ported, not restructured | `cell1_labels_s2.npz` (16 nm label grid) |
| 8 | `cell_compartments.py` (jurkat/em; paper fig 1 volumes) | ported, not restructured | `cell1_compartments.tsv` |
| 9 | `hmff_pilot_nucleus.py` + manual 20k-step FreeDTS run + `hmff_export.py` (scratch_checks) | ported, not restructured; FreeDTS (GPL-2.0, `{data_root}/external`) not bitwise reproducible -> the PLY is deposited | NE midsurface PLY |
| 10 | `jgeom_build.py` (+ library `cellshape.py`) | **restructured**: writes the first package `em_cell1_s1` (`pkgstep.publish_new`) | fields + pores |
| 11 | `jcristae.py` | **restructured** (package N -> N+1), seed 17 | cristae section |
| 12 | `jface_cap.py` | **restructured**: new package with updated summary + clip-plane f_vol / f_pm | open-face cap |
| 13 | `jspeckles.py` | **restructured** (args: src, id, chromatin npz), seed 5 | speckle section |
| 14 | `jgolgi.py` | **restructured**, seed 11 | Golgi cisternae section |
| 15 | `jgolgi_ribbon.py` | geometry **restructured** (blob replaced in slot F_GOLGI, cisternae rebuilt); pre-pass trim ported, writes copies to the work folder (belongs to the preplaced stage) | ribbon Golgi |
| 16 | `jspeckles.py` again with the v8 chromatin | as 13 | final v7 geometry |

`io.py` (jurkat/em) is the EM reader (un-flips `*_pred` in y; raw offset (0, 9, -2) s0). `jurkat/em/cytoskel/*`
(MT tracing, microvilli) is not ported here: its products enter as package layers (`mt_traces`, `microvilli`).

## Known issues before a rebuild

- **Fixed in B6 (em_cell1_v8, 2026-11-17)**: the v7 `ne_mid` field was a distance only on the nuclear side. Cause:
  `cellshape.voxelize` returns uint8 and `jgeom_build.mesh_sdf` took `~occ` (bitwise NOT: 254/255, never 0), so the
  outside EDT saw no inside voxel, every outside node stayed > 200 nm and never got the exact distance (the clamp).
  `mesh_sdf` now casts the occupancy to bool; `em_cell1_v7` keeps the defect (waived, for B1).
- `jface_cap` / `jgolgi_ribbon` / `jspeckles` depend on pre-pass inputs (chromatin, cytoskeleton): in the wholecell
  DAG they need the preplaced products of the *same* geometry (layer `depends_on`).
