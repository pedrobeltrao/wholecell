# Multi-cell viewer: data contract and design (no viewer code)

Status: design v1, 2026-11-15 (SPECS agent). Goal (DIRECTION 2026-11-15, PLANNING_NOTES): ONE WebGPU viewer, a
startup screen to pick a cell (yeast / Jurkat / future cells), every cell-specific fact read from data. Depends on
`tile_format_v2.md` (tiles, `regions.json`) and `cellgeom_v1.md` (geometry agent: membrane table, label rule,
prims, fields, layers). Audit references: H = `data/derived/repro_audit/human/`, Y = `data/derived/repro_audit/yeast/`.

## 1. Principle

A cell in the viewer = a data base URL. The app contains no list of cells, labels, membranes, regions, species,
colours or procedural shapes. Today each viewer hard-codes them in ≥ 4 places (yeast: `cell3d.py`, `shaders.ts`,
`geometry.ts`/`fields.ts`, pieces of `worker.ts` and `sim/params.ts`, Y merge_notes §6.1; human: Jurkat-only WGSL
sections for the open face, speckles, cisternal Golgi, GIC, `shaders.ts:672-831`, and region tables in `data.ts:115-141`).

## 2. `cells.json` (registry; one file per deployment)

Location: next to the app (`cells.json`), overridable by `?cells=<url>`. Written by `wholecell.export` as one
`cell.json` per build (`export/cell.json`, ARCHITECTURE §3) and merged by a deploy script; never edited by hand.

```json
{"schema": "wholecell.cells/1",
 "cells": [
  {"id": "jurkat", "name": "Jurkat T cell", "organism": {"name": "Homo sapiens", "taxon": 9606},
   "cell_type": "T lymphocyte (leukaemia line E6-1)", "state": "resting, cell 1 of jrc_jurkat-1",
   "data_base": "https://huggingface.co/datasets/<org>/<repo>/resolve/<rev>/jurkat/<build_id>/",
   "build": {"build_id": "20261116-abc1234-s20261115", "stamp": "BUILD_STAMP.json", "aliases": ["v7"]},
   "files": {"index": "tiles/index.json", "regions": "regions.json", "geometry": "geometry.json",
             "fields": "geometry_fields.bin", "meshes": "organelle_meshes.bin.gz", "bookmarks": "bookmarks.json",
             "model_doc": "MODEL.md", "features": "features.json", "library": "library/types.json"},
   "thumbnail": "thumb.webp", "stats": {"instances": 950527404, "species": 16174, "volume_um3": 712},
   "licence": {"data": "CC-BY-4.0", "attribution": "attribution.md"},
   "status": "published", "default_bookmark": "Whole cell"}]}
```

| Field | Required | Why |
|---|---|---|
| `id` | yes | deep links `?cell=<id>`; same as config `cell_id` |
| `name, organism, cell_type, state` | yes | selection card; organism mirrors config `organism` |
| `data_base` | yes | today: human `?data=` > `config.json` > `data/` (H formats §5), yeast fixed `data/` (yeast `data.ts:6`) |
| `build.build_id`, `stamp` | yes | provenance shown in About; cache key (replaces `version.json` timestamps, both viewers) |
| `files.*` | defaults as shown | relative to `data_base`; lets an old export keep its own names |
| `model_doc` | yes | both viewers already render a MODEL.md in-app (`about.ts`, `markdown.ts`) |
| `licence`, `attribution` | yes | licence blockers in H merge_notes §4 (HPA CC BY-SA, yeast repo unlicensed) must be visible per cell |
| `thumbnail`, `stats` | no | selection screen |
| `status` | no | `published` / `preview` / `archived` (legacy v3-v6 exports stay listable, provenance policy) |

## 3. Startup and deep links

| URL | Behaviour |
|---|---|
| `/` | selection screen: one card per `cells.json` entry (thumbnail, name, organism, state, instances, build id, licence) |
| `?cell=jurkat` | open that cell at `default_bookmark` |
| `?cell=jurkat&data=<url>` | cell metadata from the registry, data from `<url>` (preview a new build before registering it) |
| `?data=<url>` | legacy human override kept: no registry needed; the app reads `<url>/cell.json` if present, else treats it as a legacy (YCT1) export (§7) |
| `&view=<bookmark name>` or `&t=x,y,z&n=...&w=...` | bookmark or explicit camera, as today |
| `&build=<build_id>` | pick an archived build of the same cell (registry may list several builds per id) |

A "switch cell" button returns to the selection screen; GPU buffers of the old cell are released (one cell resident).

## 4. What the viewer reads per cell instead of hard-coding

| Today (hard-coded) | Per-cell data | Source in the pipeline |
|---|---|---|
| region names, membrane/soluble/shell id lists (`data.ts REGION_NAMES, MEM_REGIONS, SOLUBLE_REGIONS, SHELL_REGIONS`; yeast `data.ts:58-65`) | `regions.json` (`kind`, `display_as`, `colour`, `group`) | config `pack.regions` + `export.regions` |
| label names + label rule (human: free text in `geometry.json label_rule`; yeast/human: `label_from` re-implemented in WGSL and TS) | compiled label-rule table of cellgeom_v1 (`geometry.json label_rule_table`) | cellgeom_v1 manifest |
| membrane enum, half-thicknesses, face order (`shaders.ts:688`, `jenc MEMS`, yeast `HALF = 2.0` literal) | membrane table `{name, field, offset_nm, half_thickness_nm, n_side, leaflets}` | cellgeom_v1 manifest |
| field list and fixed field indices (`fields: 0 cell, 1 ne_mid, ...`, yeast `F_CELL, F_VAC, F_NUC`) | `geometry.json fields[]` by name + `geometry_fields.bin` (int16 bricks, LOD 48/32/16 nm as jenc_full) | cellgeom_v1 fields, re-bricked by export |
| procedural sections in G header slots (cristae, speckles, Golgi cisternae, open face; yeast prims `PT/PQ/PC`) | prim table + clip planes + layers of cellgeom_v1 | cellgeom_v1 `prims`, `domain.clip_planes`, `layers` |
| label / fill colours, palette indices (`palette[80+lab]`, `93u`, `92u`) | `geometry.json label_style[]` (fill colour, outline, display name) | config `export.label_style` (defaults by core label id) |
| species colour classes (`SPECIAL` by library id, yeast `CLASS_TYPES`) | `library/types.json` `class` + `features.json classes[]` | config `export.classes` |
| tile grid (`TGRID = 26`, `26*26*26`) | `tiles/index.json grid` | export |
| temperature, viscosity, ionic strength (yeast `T_K = 303.15` ×4) | `cell.json sim_defaults` | config `export.sim_defaults` (from the cell's physiology) |
| lipid compositions (`lipids.ts` Zinser 1991; human generic mammalian) | `features.json lipids` per membrane name | config `export.lipids` table with source |
| bookmarks (generated by encoders with cell-specific captions) | `bookmarks.json` (schema unchanged: name, caption, target, normal, up, width, slab, show) | config `export.bookmarks` + generated ones |

### 4a. One generic geometry evaluator (TS + WGSL)

Replace the hand-mirrored `cell3d` code with three generic pieces, written once and shared by both cells:

| Piece | Implements | Data |
|---|---|---|
| `fieldSdf(f, p)` | trilinear sampling of sparse int16 bricks (normative sampler of cellgeom_v1) | `geometry_fields.bin`, header per field |
| `primSdf(k, p)` | prim library (sphere, capsule, pill, ellipsoid, torus, slot, sheet, lamella, cisterna stack) with ops union / subtract / smooth-union k | prim table |
| `evalLabel(p)` / `evalMembranes(p)` | interpreter of the compiled label-rule table and membrane table (offsets of fields/prims, `n_side`) | label-rule table, membrane table, clip planes, pores |

The Python reference interpreter in `wholecell.geometry` and the TS/WGSL interpreter run the same table. A golden
test (export stage) samples 10⁵ seeded points per cell and compares Python vs TS (node) labels exactly and SDFs
within 0.01 nm; WGSL is checked headless when a GPU runner exists (Y merge_notes §6.1 asks for exactly this test).
Constants (`HALF`, `PERI_W`, `EIS_K`, `GOLGI_W/S`) come from the manifest; none appear as literals in shaders.
Interim: the human viewer keeps its v7 WGSL until the generic evaluator passes the golden test on Jurkat.

## 5. Per-cell plugin features, toggled from data (`features.json`)

A feature is a named renderer module compiled into the one app and **enabled only when the cell's `features.json`
lists it with its parameters**. Adding a cell never edits feature code; a new feature needs code once, then any
cell may enable it.

| Feature | Cell(s) | Data it needs | Today |
|---|---|---|---|
| `cell_wall` (procedural glucan, mannan, chitin at bud scars) | yeast | wall layer depths + budgets (cellgeom_v1 layer `cell_wall`, index `wall`) | `worker.ts:315-390`, `index.json wall` |
| `eisosomes` / `endocytic_sites` | yeast | site list layer, `EIS_K` | WGSL smooth unions `shaders.ts:752` |
| `polyp` | yeast | vacuole label name, density | yeast worker |
| `open_face` (cut cap + grey line, mesh clipping) | Jurkat (any EM crop) | `domain.clip_planes` | `faceDist`, `G1.x` |
| `microvilli` | Jurkat (planned, v8 TCR enrichment) | geometry layer of microvillus tubes (prims `slot`/capsule) | not implemented |
| `golgi_cisternae`, `speckles`, `cristae` | Jurkat; yeast cristae | prims/layers | Jurkat WGSL sections |
| `lipids` | all | per-membrane composition table | `lipids.ts` (yeast table) |
| `chromatin_annotation` (genome browser, marks) | both | `genome.json`, COLUMN `genome_bp`, flag bits | yeast `genome.json`, `annot/` |
| `sim_box` (in-browser BD + export) | both | `sim_defaults`, `library/{charges,contacts}.bin`, `diffusion.json` | yeast `sim/*`; human `dynamics.ts` |
| `motion_modes` | yeast | `library/modes.json` | yeast only (human lists it `absent`) |

`features.json`: `{"features": {"cell_wall": {"layer": "cell_wall", ...}, "open_face": {}}, "classes": [...],
"lipids": {...}}`. Missing optional files are declared (human `version.json absent[]` idea) so the app does not
request them.

## 6. Export side (what `wholecell.export` writes per build, `export/web/`)

`cell.json` (the registry entry of §2), `regions.json`, `tiles/` (YCT2), `overview/`, `library/`, `geometry.json`
(fields metadata, prims, membrane table, label-rule table, clip planes, pores, label_style), `geometry_fields.bin`,
`organelle_meshes.bin.gz` (u16 membrane ids), `bookmarks.json`, `features.json`, `MODEL.md` (copied from the
cell's doc, rendered in About), `attribution.md`, `thumb.webp`, `MANIFEST.sha256`, `encode_meta.json`. Hosting:
Hugging Face dataset, one folder per `<cell>/<build_id>` (immutable; a rebuild is a new folder); app on GitHub Pages.

## 7. Migration of the two current viewers

| Step | Human (`jurkat/viewer/app`) | Yeast (`web/app`, read-only source) |
|---|---|---|
| 0. Freeze | keep the deployed v7 site and `web_full_v7` unchanged | wait for the yeast freeze + licence (status_and_readiness §2, §6) |
| 1. One code base | start from the human fork (has shards, Range, overview LOD, `?data=`) and port yeast post-fork modules (`sim/*`, `lipids.ts`, wall, eisosomes, genome annotation) as features (Y comparison §(c)) | — |
| 2. Data readers | `regions.json`, `cells.json`, `features.json`; YCT2 decoder in worker + WGSL (`decodeInst` reads the palette) | same app |
| 3. Legacy readers | YCT1 reader kept: region names from legacy `index.json region_names`; a legacy export opens with an auto-generated `cell.json` (flags `legacy: true`) | per-tile `.bin.gz` without shards handled by the same reader (index without `shard` key) |
| 4. Geometry | generic evaluator behind a flag; golden test on Jurkat v7 geometry, then remove Jurkat WGSL sections | golden test on yeast build S geometry, then remove the hand-ported `cell3d` WGSL/TS |
| 5. Converted data | YCT1 → YCT2 converter (`tile_format_v2.md` §8) for archived human packs if wanted | converter for build S tiles; the original stays the archival record |
| 6. Deploy | one site, `cells.json` lists Jurkat (v7 build + aliases) and testcell; yeast added after freeze | — |

## Open questions

1. Hosting of the registry: one `cells.json` in the app repo (needs a site rebuild per new cell) or fetched from the
   data repo (no rebuild)? Proposed: app ships a default, `?cells=` overrides.
2. Several builds per cell in the selection screen (v7 and a later v8) or only the latest plus `&build=`?
3. Label-rule table format is owned by cellgeom_v1.md; the viewer needs it to be expressible in WGSL without
   recursion (flat ordered tests). To be confirmed with the geometry agent.
4. Should feature modules be lazily loaded bundles (smaller app) or compiled in (simpler)? Proposed: compiled in.
5. Who owns MODEL.md per cell (pipeline-generated parts such as counts vs hand-written text)?
