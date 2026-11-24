# Tile format v2 (`YCT2`): packed instances for the viewer and downstream consumers

Status: spec v2.0, 2026-11-15 (SPECS agent). Implementation status 2026-11-24: PALETTE, INST, GHOST, POLY (4, 5, 6) and
MESH are written by `wholecell.export` (`yct2.py`, `meshpoly.py`); COLUMN and LOD 1/2 are not yet written (index.json
`sections_written` / `not_written`). Details fixed by the implementation and the viewer reader are marked "(impl.)". Successor of `YCT1` (yeast `pipeline/webenc.py`, human
`jurkat/viewer/jenc.py` / `jenc_full.py`). Written by `wholecell.export` (stage folder `export/`, ARCHITECTURE §3).
Audit references: H = `data/derived/repro_audit/human/`, Y = `data/derived/repro_audit/yeast/`.

## 1. Two levels, one convention

| Level | File | Precision | Reader | Why |
|---|---|---|---|---|
| P (pack tile) | `pack/tiles/<i>.npz` (lossless, written by the packer) | float32 nm, model frame | box exporter, validation, encoder | both packers already write this schema (H formats §4, Y formats §4); simulations must not inherit viewer quantisation (Y formats §6c: `instances.npz` decoded from viewer tiles) |
| V (stream tile) | `export/web/tiles/s_*.bin` (`YCT2`, quantised, gzip, sharded) | ≤ 0.003 nm, ≤ 1e-4 rad | viewers, light analyses | one container for both cells: human shards scale to 10⁹ instances, yeast carries the richer per-instance semantics (Y comparison §b, "merge both") |

Both levels share: nm, the model frame of the cellgeom_v1 package (`docs/reference/cellgeom_v1.md`, frame + handedness
recorded there), world = R(q)·local + p (Y merge_notes §5), library +z = N-side of the membrane, and the
quaternion order of §2.

### 1a. Level P schema (`pack/tiles/<i>.npz`)

| key | dtype | meaning | change vs legacy |
|---|---|---|---|
| `types` | int32 (n) | species index into `pack/libs.json` | — |
| `pos` | float32 (n,3) | model-frame nm | — |
| `quat` | float32 (n,4) | **(w,x,y,z)**, unit, w ≥ 0 | sign canonicalised (new) |
| `region` | uint16 (n) | index into `pack.regions` (config order) | **new**: both legacy builds rebuild it from `np.repeat(job_region, job_placed)` (H formats §4, Y formats §4); stored so instance order may change |
| `flags` | uint16 (n) | core bits of §5 (bits 0-3) | new |
| `job_lib, job_region, job_n, job_placed, overflow` | as legacy | request table | — |
| `mrna_pts, mrna_off` | float32, int64 | in-tile polylines (optional) | — |

For the Jurkat reproduction (criteria B12) the port writes `region`/`flags` as extra keys; the legacy keys stay
byte-comparable to v7.

## 2. Quaternion order: (w, x, y, z), scalar first, on disk everywhere

What the code does today (checked in source):

| Component | Order | Evidence |
|---|---|---|
| pack3d `quat_mat`, `rand_quat`, `quat_z_to` (both forks) | w,x,y,z | `jurkat/yc/pipeline/pack3d.py:80` (`a,b,c,d = q[0..3]`, a = scalar) |
| pack tiles, pre-pass npz (both cells) | w,x,y,z | Y formats §0 "scalar-first in all packing arrays"; human builders convert scipy output (`jcytoskeleton.quats`, `actin_patches._quat`) |
| YCT1 encoders | w,x,y,z (smallest-three index 0 = w) | `jenc.quat_encode`, yeast `webenc.quat_encode` docstring "q (n,4) w,x,y,z" |
| TS decoders (`worker.ts instPose`, both viewers) | w,x,y,z | `jurkat/viewer/app/src/worker.ts:277` |
| WGSL `qrot` / sim `engine.ts` | x,y,z,w (vec4f, GPU-internal) | `shaders.ts:45`; `decodeInst` swizzles `o.q = q4.yzwx` |
| scipy `Rotation.from_quat/as_quat` | x,y,z,w | converted at every call site (`np.c_[q[:,1:], q[:,:1]]`) |
| legacy `jurkat/yc/pipeline/export.py placements_table` | x,y,z,w (columns `qx..qw`) | the only scalar-last file on disk; labelled columns |

Decision: **(w,x,y,z) for every file, array and wire format** (P, V, box exports, sim schema). It is what the
packers, all npz files of both builds, both YCT1 encoders and both TS decoders already use, so the choice needs no
data conversion and no packer change (v7 bit-identity, criteria B12). Scalar-last appears only inside GPU code and
third-party APIs, converted at the boundary with the helpers `wholecell.export.quat.{to_xyzw,from_xyzw}` (Python)
and `quat.ts` (TS). Canonical sign: w ≥ 0 (as `jcytoskeleton.quats`), so files are deterministic.

| Target | Native order | Conversion from wholecell (w,x,y,z) |
|---|---|---|
| numpy / pack3d / numba kernels | w,x,y,z | none |
| scipy `Rotation` | x,y,z,w | `Rotation.from_quat(q[:, [1,2,3,0]])`; back: `r.as_quat()[:, [3,0,1,2]]` (scipy ≥ 1.14 also accepts `scalar_first=True`) |
| three.js `Quaternion(x,y,z,w)` | x,y,z,w | `new Quaternion(q[1], q[2], q[3], q[0])` |
| gl-matrix `quat` | [x,y,z,w] | `quat.fromValues(q[1], q[2], q[3], q[0])` |
| WGSL (`qrot`) | vec4 (x,y,z,w) | swizzle once in `decodeInst` (as today) |
| Mol* `Quat` (linear-algebra) | [x,y,z,w] | `Quat.create(q[1], q[2], q[3], q[0])`; usually pass a `Mat4` instead (`Mat4.fromQuat`) |
| mmCIF export | no quaternion | explicit transformed coordinates (yeast box export), or `_pdbx_struct_oper_list.matrix[i][j]` + `vector[i]` from R(q) (row-major, Å) |
| OpenMM / HOOMD | rigid bodies: HOOMD `orientation` is (w,x,y,z) | none for HOOMD; OpenMM gets coordinates |

## 3. Grid, LOD and sharding (all values from config `export.tiles`, recorded in `index.json`)

| Item | Rule | Yeast | Jurkat | Why |
|---|---|---|---|---|
| `tile_nm` (viewer) | config; pack tiles re-binned by instance centre | 200 (= pack tile) | 200 (pack 500) | both viewers stream 200 nm today (H formats §5, Y formats §5) |
| `block_tiles` | block edge = `block_tiles·tile_nm`, must be a multiple of the pack tile | 5 (1 µm) | 5 (1 µm = 2×2×2 pack tiles) | jenc_full alignment of both grids (H formats §5) |
| grid origin / extent | `grid.origin_nm`, `grid.min_ijk`, `grid.max_ijk` in index | from geometry domain | from geometry domain | replaces yeast `TGRID = 26` / `26*26*26` (Y merge_notes §6.9) |
| tile key | `ijk = floor((p - origin)/tile_nm)` | — | — | as both encoders |
| LOD 0 | instances (this binary) | | | |
| LOD 1 | proxies: ≤ `proxies_per_tile` (128) stand-in spheres per tile, one gzip file per block `overview/BX_BY_BZ.bin.gz` (`YCO2` = YCO1 with u16 region) | | | human overview level (H formats §5) |
| LOD 2 (optional) | per-block region mass/volume table in `overview/blocks.json` | | | whole-cell colour fill before LOD 1 |

**Shards.** `tiles/s_BX_BY_BZ[.k].bin` = concatenated gzip members, one member per tile, members in (ix,iy,iz) order.
A block whose shard would exceed `shard_max_mb` (default 45) is split **at member boundaries** into `.0, .1, …`
(not byte-split as `split.json` does today, so Range offsets stay valid). Targets: < 10k files per folder and
< 100k per repo (Hugging Face guidance, H jenc_full docstring); Jurkat v7 ≈ 1,175 shards + 1,175 overview files.
Fetch: HTTP Range `[off, off+gz_bytes)`; fallback = whole shard (servers without Range). Large non-tile files
(beads, fields) may still be split by `split.json` (byte parts, concatenated by the reader).

**Determinism.** Every gzip member: `gzip.compress(blob, compresslevel=6, mtime=0)` (header OS byte fixed by
Python); fixes the human defect where file hashes differ per run (H stages §8, rebuild_test §5). zlib version is
recorded in `encode_meta.json`; equality is defined on decompressed members (crc32 in the index), file sha256 is
expected equal only with the same zlib.

## 4. Tile binary `YCT2` (little endian, sections 4-byte aligned)

Header: 32 × u32 (128 B).

| word | content |
|---|---|
| 0 | magic `YCT2` |
| 1 | format version: major << 16 \| minor (2.0) |
| 2-4 | ix, iy, iz (i32) |
| 5 | lod (0) |
| 6-8 | `tile_nm`, `q_org` (nm below tile origin), `q_range` (f32) |
| 9 | n_inst |
| 10 | n_ghost |
| 11 | n_palette |
| 12 | n_sections (entries of the section table) |
| 13 | flags (bit 0: u32 triangle indices) |
| 14 | `lib_tag` = first 4 bytes of sha256 of `library/types.json` ids list (refuses a tile/library mix) |
| 15 | `build_tag` = first 4 bytes of sha256(build_id) |
| 16-31 | reserved (0) |

Section table (n_sections × 4 u32: kind, byte offset from tile start, byte length, count), then the sections. Kinds:
`1 PALETTE, 2 INST, 3 GHOST, 4 POLY (mRNA), 5 POLY (linker), 6 POLY (extra), 7 MESH, 8 COLUMN`. Unknown kinds are
skipped by readers (forward compatibility; yeast v1 had to reuse "reserved" header words, Y formats §5).

**PALETTE** (n_palette × 2 u32): `species` (u32, index into `library/types.json`), `region` (u16) \| `pflags` (u16,
per-entry defaults). Instances are sorted by (species, region) as in YCT1, so the palette is the run table.

**INST / GHOST** (16 B per instance, `vec4u` on the GPU, same word layout as YCT1 except w0):

| field | encoding | YCT1 |
|---|---|---|
| w0 | palette index (16 bits) \| flags (14 bits) << 16 \| qidx (2 bits) << 30 | type 14 b, region 5 b, flags 11 b |
| w1 | px \| py << 16 | same |
| w2 | pz \| qa << 16 | same |
| w3 | qb \| qc << 16 | same |

Consequences: species and regions are unbounded (u32 / u16 via the palette) while the record stays 16 B. YCT1
caps were 16,383 types (Jurkat v7 uses 16,174, H `web_full_v7/encode_meta.json`) and 32 regions (Jurkat v7 needs 37
and aliased 5 of them away, H formats §5; yeast appended ids 24-26 at the cap, Y merge_notes §6.8). A tile with
> 65,535 palette entries is split by the encoder (never seen: a 200 nm tile holds < 2·10⁴ instances).

**Position**: `p = lo - q_org + u16 · q_range/65535`, `lo` = tile origin. Default `q_org` 80, `q_range` 360 nm
(window holds the tile and ghosts within 80 nm; step 0.0055 nm, max error 0.0027 nm). uint16 rather than float16:
float16 over a 360 nm window has 0.06-0.25 nm steps, non-uniform. **Orientation**: smallest-three over (w,x,y,z),
qidx = dropped largest component, others int16/32767 · 1/√2; per-component error ≤ 1.1e-5, rotation error
≤ 1e-4 rad (≤ 0.002 nm at a 20 nm radius). Ghosts: near-membrane instances of the 26 neighbours (lipid exclusion),
same encoding, flag bit 3 set.

**POLY sections** (YCT1 layouts kept): `u32 len[n]`, `i32 start[n×3]` in `poly_q` units from `lo`, `i8 delta[...]`;
extra polylines prefixed by `u16 kind[n]`. Kind names and `poly_q` per section come from `index.json`
`polyline_kinds` (yeast 4 = pre-rRNA, 5 = IDR; human mRNA 64-point pieces, nascent RNA), not from code.
(impl.) Section count = number of polylines; `kind[n]` is padded to 4 bytes, and each section is padded to 4 bytes.
`poly_q`: 4 mRNA 0.005 nm, 5 DNA 0.01 nm, 6 extra 0.005 nm. Pieces have at most 64 points and belong to the tile of
their middle point (`n_poly_beads` = points of all POLY sections of the tile). Extra kinds = viewer free-sphere kinds
(0 RNA, 4 pre-rRNA, 5 FG repeats; 2/3 reserved for the yeast wall polysaccharides).

**MESH**: `u16 vert[nv×3]` (position quantisation), `u16|u32 tri[nt×3]`, `u16 mem[nt]` (membrane id of the
cellgeom_v1 membrane table; YCT1 u8).
(impl.) Each array is padded to 4 bytes; u32 triangle indices iff header word 13 bit 0 (nv > 65535). Section count =
nt; nv is the index row `n_vert` (the viewer reads both from the index) and follows from the section length
(nv = floor((len - tri bytes - pad4(2 nt)) / 6)). `mem` = the membrane's `engine_code` (= `geometry.json`
`membranes[].code`, the column of `cellgeom_eval_v1` `mems`). Winding: face normal (b - a) × (c - a) points to
positive SDF (cytoplasmic / N side). Index rows carry `n_vert`, `n_tri`; index.json `mesh` records `mesh_h_nm` and
the meshed membranes.

**GHOST** (impl.): ghost words index the same PALETTE; (species, region) pairs used only by ghosts are appended after
the tile's own runs (`n_palette` counts both), so the INST words and the palette prefix they use are identical to a
tile without ghosts. Ghost source set = index.json `lipid_mask_regions` (membrane regions + config); flags ghost |
near_membrane. INST flags do not carry `near_membrane`; viewers treat instances of `lipid_mask_regions` as lipid
maskers.

**COLUMN** (typed per-instance extras, aligned with INST order): header `u16 name_id, u8 dtype, u8 width`, then
data. `name_id` indexes `index.json` `columns[]` (`{name, dtype, unit, description}`); dtypes `u8 u16 u32 i32 f16
f32`. Examples: `genome_bp` (u32, nucleosome dyad), `chain` (u32), `marks` (u8 bitmask, if not in flags), `plddt`
(u8). Readers ignore unknown columns.

## 5. Flags (14 bits per instance; names in `index.json` `flag_bits`)

| bit | name | meaning | source |
|---|---|---|---|
| 0 | `near_membrane` | centre within 3 nm + BR of a membrane | yeast bit 0 |
| 1 | `fixed` | anchored / fixed obstacle in simulations | yeast bit 1 (`simgpu FIXED_IDS`) |
| 2 | `prepass` | placed by a pre-pass module, not the packer (format name kept; the build stage is `preplaced` since A0c) | human `src` 0/1 of BD boxes (H formats §6) |
| 3 | `ghost` | copy from a neighbour tile | — |
| 4-13 | cell-declared | e.g. yeast histone marks (7 bits, order `genome_map.MARKS`); Jurkat: free | Y formats §5 (bits 2-8 undocumented in the docstring) |

Bits 0-3 are core and fixed; bits 4-13 are declared per export in `index.json` (`{bit, name, source}`), from
config `export.flag_bits`.

## 6. Dataset files written next to the tiles (`export/web/`)

| File | Content |
|---|---|
| `tiles/index.json` | `format: "YCT2"`, version, `build_id`, `tile_nm`, `block_tiles`, `grid`, `q_org`, `q_range`, `lod[]`, `flag_bits`, `columns`, `polyline_kinds`, `tile_columns` (named; count checked: yeast S had 13 values for 11 names, Y formats §5) and `tiles` rows `[ix,iy,iz,shard,off,gz_bytes,n_inst,n_ghost,n_poly_beads,n_vert,n_tri,label_mask,crc32]`; gzip-compressed copy `index.json.gz` |
| `regions.json` | per-export region table: `[{id, name, kind: volume|membrane|shell|virtual|prepass, geometry_refs (labels/membrane names of cellgeom_v1), display_as, colour, group}]`; ids = config `pack.regions` order, then pre-pass groups. Replaces `RNAMES_V5/6/7`, `jenc.region_map`, `PACK_ALIAS`, yeast `REGION_NAMES` and the viewer copies (`data.ts`) |
| `library/types.json`, `beads.bin`, `spheretree.bin` | as YCT1 (int16 beads at 0.005 nm); `types.json` gains `lib_sha256` per type (from `library/library_manifest.tsv`) |
| `MANIFEST.sha256` | sha256 of every file of the export (shards included) |
| `encode_meta.json` | build id, config sha256, encoder commit, zlib version, counts (instances per region = pack counts), coarsening |
| `version.json` | cache-buster `{v: build_id}` (both viewers) |

Geometry, bookmarks and cell metadata are specified in `viewer_multi_cell.md`.

## 7. Validation (export stage, `validate/` checks)

1. Decoded instance count per (species, region) = pack counts (H rebuild_test §8 "viewer data" row).
2. Round trip on a seeded 1 % of tiles: |Δp| ≤ q_range/131070, rotation error ≤ 1e-4 rad.
3. Every palette region id exists in `regions.json`; every species < len(types); `lib_tag` matches.
4. Shards: gzip members parse, crc32 matches the index, file sizes ≤ `shard_max_mb`.
5. Re-encode twice → identical `MANIFEST.sha256` (gzip determinism).

## 8. YCT1 → YCT2 converter (`python -m wholecell.export.yct1_to_yct2`)

Input: a legacy viewer data folder (yeast per-tile `tiles/<i>_<j>_<k>.bin.gz`, human shards + `index.json`, human
1.6 µm cube). Per tile: read the 16-word header; build the palette from (type 14 b, region 5 b) runs; region names
from the legacy `index.json region_names` (both builds write them) into `regions.json`, marking aliased human
regions (`pack_regions.viewer_id`) as unrecoverable; flags remapped (yeast bit 1 → `fixed`, bits 2-8 → declared
bits 4-10); w1-w3 copied verbatim (identical quantisation when q_org 80 / q_range 360, so poses are lossless);
polylines copied (yeast `hdr[13:15]` = nExtra/nExtraBeads); mesh `mem` u8 → u16; then shard. Check: decode v1
and v2 and compare poses bit-exactly. Legacy data stay untouched (provenance policy); converted copies are new
folders.

## Open questions

1. Container name: keep `YCT2` (lineage, yeast merge_notes §10.5) or a neutral `WCT2`? Assumed `YCT2`.
2. Human mRNA pieces vs yeast whole polylines + `annot/<tile>.json`: should polymer identities become COLUMNs on
   POLY sections in v2.0 or stay as side files?
3. Is 16 B per instance on the GPU still the binding constraint, or may w0 grow (20 B) to drop the palette?
4. Should LOD 1 proxies be per tile (as now) or per block octree for > 10 µm cells?
5. Index size for much larger cells (10⁶ tiles): switch `index.json` to a binary per-block index?
