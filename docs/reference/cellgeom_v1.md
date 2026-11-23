# cellgeom_v1: cell geometry package

Status: v1, 2026-11-15 (§4.3 rule table: 2026-11-16; yeast layout, prims: 2026-11-17). Code: `wholecell/geometry/cellgeom.py` (write / load /
validate), evaluator `wholecell/geometry/engines/cellgeom_eval_v1.py` with its rule tables (§4.3); legacy engine
`wholecell/geometry/engines/sdf_sections_v1.py` (reproduction tests only). Origin: proposal (e) of
`data/derived/repro_audit/yeast/comparison_with_human.md`; deviations are listed in §12.

A package is the **only** geometry input of a build. Toolkits write packages (`toolkits/em`, `toolkits/analytic`);
the pipeline, packer, viewer encoder and box exporters read them.

## 1. Folder, immutability, identity

```
{geom_root}/<cell_id>/<geometry_id>/      package_id = "<cell_id>/<geometry_id>"
  manifest.json                           everything below is declared here; sha256 + size per file
  fields/<name>.npz                       one sparse-brick SDF per field (§3)
  sections/<name>.npz                     procedural sections (§6); optional prims.npz (§6.2)
  pores.npz                               NPC pores (if the cell has a nucleus with pores)
  meshes/<name>.ply                       optional meshes (nm; frame stated in the layer entry)
  layers/<layer>/...                      optional layers (§8)
  validation.json                         validator output written at creation (not in manifest.files)
```

- Write-once: `PackageWriter` builds in `<dest>.tmp-<pid>`, writes `manifest.json`, validates, then renames;
  files are made read-only. An existing folder is never rewritten: **a new version = a new geometry_id folder**.
  A package that fails validation stays unpublished in the tmp folder (unless waived, §10).
- `load(dir)` verifies size + sha256 of every file in `manifest.files` (`verify="quick"`: sizes only).
- Builds pin a package by path and the **sha256 of manifest.json** (cell config `geometry.manifest_sha256`).
- Toolkit steps that refine a geometry read package N and write a new package N+1 (never in place).

## 2. manifest.json (top-level keys)

| key | content |
|---|---|
| `spec_version` | `"cellgeom_v1"` |
| `cell_id`, `geometry_id`, `package_id`, `created`, `description` | identity |
| `organism`, `cell_line`, `source` | dataset id/DOI, licence, voxel size, crop, legacy file + md5 |
| `generator` | toolkit, repo commit, dirty flag, method, params / seed (wholecell.seeds stream) |
| `engine` | `name` (evaluator module), `layout_id`, `params` (physical constants, §5), `field_slots` {field: slot} |
| `frame` | units, handedness, axes, origin, `source_affines` (§9) |
| `domain` | `enclosing_radius_nm`, `clip_planes` (§7) |
| `fields` | {name: slot, file, dtype, h_nm, brick, origin_nm, bricks, clamp_in/out_nm, bricks_stored, sign, provenance} or `{slot, empty: true}` |
| `labels` | vocabulary (§4) incl. engine codes |
| `membranes` | membrane table (§5) |
| `label_rule` | declarative rule (§4.2) |
| `sections` | {name: type, file, params, description, engine_params} (§6) |
| `pores` | file, n, columns, `source_codes`, `counts`, `min_sep_nm`, provenance |
| `layers` | {name: kind, schema_version, depends_on, files, refs, provenance} (§8) |
| `summary` | per-field `{vol_sdf_um3, area_um2, ...}` + `open_face: {f_vol, f_pm}` (§7); extra keys allowed |
| `validation_waivers` | {check: reason} documented accepted failures (§10) |
| `files` | {relpath: {sha256, bytes}} |

## 3. Sparse-brick SDF fields (shared yeast `cellshape.bricks` / Jurkat `jgeom_build.bricks` encoding)

One npz per field: `header` float64[12], `index` int32, `data` float32 (float64 if float32 is not lossless).

| header slot | 0-2 | 3 | 4-6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|
| value | origin xyz (nm, node 0) | spacing h | bricks nbx nby nbz | voxels per brick B (8) | inside clamp | outside clamp | 0 | len(index) |

- index (nbx·nby·nbz, C order x,y,z): `-1` brick fully inside (value -clamp_in), `-2` fully outside (+clamp_out),
  `k >= 0` stored brick k. data: (B+1)³ node values per stored brick (shared faces duplicated).
- Slots 10/11 are offsets relative to a concatenated [index, data] block; the loader rebases them into G.
- Sampler (normative): trilinear inside a stored brick, clamps elsewhere and outside the grid (+clamp_out).
- Sign: negative inside the volume the field encloses. Grid h <= 20 nm (Jurkat/testcell 16), B = 8. Band
  (exact distances) >= the largest membrane offset + 2h; record `band_nm` per field.
- **Absent field** (`empty: true`, no file): evaluates to +clamp everywhere = outside everywhere.

## 4. Labels

### 4.1 Vocabulary (no cap on the number of labels)

Each entry: `{id, name, go, description, engine_code, parent?}`. Shared core ids are < 1000 and stable
(`cellgeom.CORE_LABELS`); cell-specific **extension ids are >= 1000** and must name a core `parent`.

| id | name | GO | id | name | GO |
|---|---|---|---|---|---|
| 0 | exterior | GO:0005576 | 10 | golgi_lumen | GO:0005796 |
| 1 | cytosol | GO:0005829 | 11 | endolysosome_lumen | - |
| 2 | nucleoplasm | GO:0005654 | 12 | vacuole_lumen | GO:0005775 |
| 3 | nucleolus | GO:0005730 | 13 | peroxisome_matrix | GO:0005782 |
| 4 | perinuclear_space | GO:0005641 | 14 | lipid_droplet | GO:0005811 |
| 5 | er_lumen | GO:0005788 | 15 | cell_wall | GO:0005618 |
| 6 | nuclear_pore_channel | GO:0044613 | 16 | nuclear_speck | GO:0016607 |
| 7 | mito_ims | GO:0005758 | 17 | golgi_intercisternal | - |
| 8 | outside_imaged_volume | - | 18 | crista_lumen | GO:0061617 |
| 9 | mito_matrix | GO:0005759 | 19 | lysosome_lumen | GO:0043202 |

**Engine codes** are a compact per-package mapping (0..63, so packer label bitmasks fit in int64), declared per
label (`engine_code`). Labels without an engine code are vocabulary only (e.g. for layers). Regions and inventory
refer to labels by **name**, never by engine code or id.

### 4.2 Declarative label rule

Ordered list; first matching entry wins; an entry is `{if: [tests...], label: name}` or `{if: [...], then: [...]}`
(nested list, must end with an unconditional entry). Tests (all must hold):

| test | meaning |
|---|---|
| `{clip: name}` | point is on the removed side of clip plane `name` |
| `{mem: M, lt/gt: v}` | membrane SDF of membrane M (engine value) |
| `{field: F, lt/gt: v}` | raw field F (absent field = +inf) |
| `{var: matrix, lt/gt: v}` | matrix SDF (= -MITO_IM membrane SDF) |
| `{in_pore: rc, zwin: w}` | within rc of a pore axis, axial \|z\| <= w |
| `{in_section: S}` / `{has_section: S}` | inside a sphere of section S / section S exists |

`v`, `rc`, `w` are numbers or expressions over `engine.params` (e.g. `"-PERI_W-PORE_LW"`). The rule is for viewers
and future generic engines; the validator checks it reproduces the engine labels exactly on sampled points.

### 4.3 Label rule table (the evaluator's input; `cellgeom_eval_v1`)

The core evaluator `cellgeom_eval_v1` has no hand-coded cell logic. It reads a **rule table** (YAML, one per package
layout) with two parts, compiles it (`wholecell/geometry/rules.py`) into numeric rows appended to G, and its numba
kernels (`mem_sdfs`, `label_from`) loop over those rows. **One place:** `config/geometry/rules/<layout>.yaml`, where
layout = manifest `engine.name`; the cell config names it in `geometry.rules` (Jurkat and testcell:
`config/geometry/rules/sdf_sections_v1.yaml`). Packages are not changed; a future toolkit may write the table into new
packages.

**Membrane table** (`membranes:`): rows evaluated in order; every row except `smooth_union` defines the next engine
membrane (row k = code k, checked against the manifest; the number of membranes is the package's, `n_mems(G)` <= 16:
Jurkat/testcell 8, yeast 10). A row may use only membranes defined above it. Kinds:

| kind | keys | SDF |
|---|---|---|
| `field` | `field`, `sign` (±1), `offset` (value, or `[a, b]` = a − b) | sign × field − offset |
| `wall` | `of` | min(of, height above the open face) (domain wall) |
| `shell_pores` | `field`, `field_offset`, `width`, `pore_r`, `zwin` | \|f + field_offset\| − width/2; inside a pore cylinder: torus of radius pore_r + width/2 |
| `inner_cristae` | `of`, `offset` | matrix = of + offset carved by section `cristae`; SDF = −matrix |
| `golgi_cisternae` | `field`, `clear_of` [2 membranes] | blob field, or cisternal slabs of section `golgi_cis` inside it away from the two membranes |
| `prims` | `prims` (class), `base_field`, `ref`, `ref_radius` | min(base field, prims of the class) (no base field: 1e9); a cER-patch prim is a sheet at depth ref_radius − q5 below membrane `ref` |
| `prims_difference` | `prims`, `minus` (classes), `sign` | m = max(min prims, −min minus-prims); SDF = sign × m (yeast IM: matrix = IM capsules minus crista lumens) |
| `smooth_union` | `into` (membrane), `prims`, `k`, `side` (±1) | modifies `into`: x = side × into, d = min prims; if d < x + k: into = side × smin(x, d, k) (yeast eisosome furrows into the PM, bud scars on the outer wall) |

`fields:` gives analytic **fallbacks** for fields a package lacks (`sphere` {centre, radius}, `ellipsoid_approx`
{centre, axes, rot}: the yeast nucleus first-order distance); the fallback replaces the field everywhere (membrane
rows and `field:` terms). `shapes:` names analytic label shapes (`ball` {centre, radius}). `constants:` {file, key,
names} reads named scalars from a checksummed package file (yeast: the cell3d constants in `native/geometry.json`);
values may use them like `engine.params`. `clearance:` sets the packer clearance `min_mem_dist_g`: `thickness_corrected`
(min |d| − half thickness, + HALF; `wall` / `wall_exempt_for`: Jurkat CELL_OUTER is no wall for PM proteins) or `plain`
(yeast legacy: min |d| over all membranes, uniform 2 nm added by the caller).
Prims (`prims.npz`, §6.2) are loaded with the package; their types are mapped by name to the evaluator's types
(`PRIM_TYPES`: sphere, capsule, pill, cerpatch, slot, torus; formulas of the yeast cell3d) and classes are named by
`manifest.prims.classes`. The packer may pass a subset of prims (`pidx`, e.g. by bounding sphere).

Half thicknesses are package data (`membranes[].half_thickness_nm` = `engine.params.MEM_HALF`), not rule data.

**Label rule** (`labels:`): ordered rows, **first match wins**; a row is `{label, if: [conditions], requires: [...]}`
with a comment. `label` is a package vocabulary name (mapped to its `engine_code`). A condition is the string
`"<term> <op> <value>"`, op `<` or `>`, value a number or an expression over `engine.params` / constants; all
conditions must hold; no `if` = always (the last row must be unconditional). Terms: a membrane name (its SDF),
`field:<name>` (field or its fallback), `face` (height above the open face), `pore` (−1 inside the channel of a pore of
the `shell_pores` membrane, else +1), `speckle` (−1 inside a nuclear speckle sphere, else +1), `shape:<name>` (−1 inside
the named shape). `requires` lists package items (`field:<name>` non-empty, `section:<name>`, `clip:<name>`); a row
whose items are missing is dropped at compile time. Dropping is exact (such a row can never match) and lets one table
serve packages with fewer parts (testcell has no Golgi, endosomes, speckles, face).

| layout (manifest engine.name) | rule table | packages |
|---|---|---|
| `sdf_sections_v1` | `config/geometry/rules/sdf_sections_v1.yaml` (8 membranes, 14 labels) | Jurkat em_cell1_v7, testcell analytic_v1 |
| `yeast_cell3d_v1` | `config/geometry/rules/yeast_cell3d_v1.yaml` (10 membranes, 15 labels: two wall layers, cER lumen, vacuole, LD, peroxisome, nucleolus ball) | yeast analytic_x5_s* |

Exactness: the kernels perform the same float operations, in the same order, as the hand-coded `sdf_sections_v1` and
the vendored yeast `cell3d`, and the thresholds are IEEE-exact (e.g. −0.5·PERI_W). Two yeast forms follow what LLVM
fastmath makes of cell3d's constants: the PM offset is kept as two numbers (`dep − (CELL_R − PM_R)`, compiled as
`(dep − CELL_R) + PM_R`), and `smooth_union` uses smin = min(a, b) − t²·(0.25/k), t = max(k − |a − b|, 0), the form
cell3d's smin with a compile-time k folds into. Checked by Jurkat criterion B1 at level E on `em_cell1_v7` (labels and
the 8 membrane SDFs identical on 2.55 M points), by `tests/check_a1_yeast_geometry.py` on the yeast s0 port package
(labels, the 10 membrane SDFs and the clearance identical on 2.55 M points + box selections, under the numba of the
vendored env, 0.61.2; numba 0.68 changes cell3d's own ER values at the ulp level, see the A1 commit) and by
`tests/test_geometry.py` (python reference `wholecell/geometry/eval_ref.py`, which interprets the YAML directly, vs the
numba kernel on ~10⁴ points per package; yeast also vs the vendored cell3d and without fields vs cell3d's fallbacks).
The declarative `label_rule` of the manifest (§4.2) stays as package metadata and is still validated against the
evaluator.

Limits: at most 4 conditions per row; at most 16 membranes per package. The core packer still sizes its buffers with
the sdf_sections_v1 constant `NMEM = 8` (yeast packing in the core: A6/A7).

## 5. Engine, parameters, membranes

`engine.name` names the package **layout** (the engine it was written for). The core evaluates every package with
`cellgeom_eval_v1` (config `geometry.evaluator`; no environment switch) plus the rule table of its layout (§4.3);
`load(..., evaluator=engine.name)` still evaluates with the legacy engine itself (reproduction tests; the yeast layout's
own engine is the vendored cell3d, `plugins/yeast/geom_engine.load`).
`sdf_sections_v1` (Jurkat v7 `cell3d`, generalised; now the layout of Jurkat and testcell) has fixed
field slots (cell, ne_mid, mito, er, golgi, endo, nucleolus), 14 label codes and 8 membrane codes (PM, NE, ER,
MITO_OM, MITO_IM, GOLGI, ENDO, CELL_OUTER); `engine.params` fills header slots of G (slot map in the engine
docstring): PERI_W, IMS_W, PORE_R, HALF, CR_W, CR_GAP, CR_JW, CIS_CLR, PORE_ZWIN, PORE_LW, MEM_HALF[8], and the optional
CR_RIM (crista rim rounding radius, slot 40; absent = 0 = square lamella caps; B6b).

Membrane table entry: `{name, engine_code, field, derivation, offset_nm, half_thickness_nm, n_side, p_side,
leaflets, go}`. Sign convention: the engine membrane SDF is negative on `p_side` (lumenal / non-cytoplasmic) and its
gradient points to `n_side` (= +z of membrane-protein libraries). `half_thickness_nm` must equal
`params.MEM_HALF[engine_code]` (Jurkat: Pogozheva 2022 D_P-P/2). Example (Jurkat): NE = `|ne_mid| - PERI_W/2` with
pore tori, surfaces ONM / INM / pore membrane; MITO_IM = matrix boundary (`mito + IMS_W`, cristae section), n_side
mito_matrix; CELL_OUTER = domain wall `min(PM depth, z - z_face)`, leaflets none.

## 6. Procedural sections and analytic primitives

### 6.1 Sections (`sections/<name>.npz`, evaluated by the engine)

| type id | arrays | used by |
|---|---|---|
| `cristae_lamellae_v1` | `section`: n, id-grid origin xyz, spacing, dims xyz, n x 12 lamella params, dense id grid | MITO_IM, crista lumen |
| `speckle_spheres_v1` | `section`: n, n x (cx, cy, cz, r) | label nuclear_speck |
| `golgi_cisternae_v1` | `W`, `S`, `header`, `index`, `data` (deep Golgi SDF, §3 encoding) | GOLGI, golgi_intercisternal |

Section internals use offsets relative to the section start. An engine lists the types it evaluates
(`SECTION_TYPES`); `load` refuses unknown types.

### 6.2 Analytic primitive table (optional `prims.npz`, future engine `analytic_prims`)

Yeast-style: `type` int16, `params` float64 (n, 8), `comp` int16, `field_bind` int16, `op` int8 (union / subtract /
smooth-union, with `k`), `bsphere` (n, 4); type enum in `manifest.prims.types`, classes in `manifest.prims.classes`.
Written by `toolkits/analytic/yeast.py` (types sphere, capsule, pill, cerpatch, slot, torus); evaluated by
`cellgeom_eval_v1` through the rule table (§4.3: `prims`, `prims_difference`, `smooth_union` rows; the table, not
`op`, says how prims combine). Planned types (ellipsoid, sheet, lamella, cisterna_stack) and baked fields for
`sdf_sections_v1` consumers are not produced.

## 7. Boundary, summary

`domain.clip_planes`: `{name, normal, point, kept_side: "+normal", label: outside_imaged_volume, reason,
cap_method, f_vol, f_pm}`. Fields are extruded through an open face (no membrane on it); beyond the plane the label
is `outside_imaged_volume` and the face is a hard wall (CELL_OUTER). `sdf_sections_v1` supports one z plane.
`f_vol` / `f_pm` = modelled / (modelled + missing cap) cell volume / PM area, used to scale copy numbers.

`summary` (legacy-shaped, read by the packer): per field `vol_sdf_um3`, `area_um2` (+ `ne_mid.vol_inside_um3`), and
`open_face: {f_vol, f_pm}` (1.0 without a clip plane). Converters may keep a legacy summary verbatim (Jurkat).
Exposed as `cell.summary`.

## 8. Optional layers

`{kind, schema_version, depends_on, files, refs, provenance}`. `depends_on` must equal this `package_id` (checked);
`files` are package files (checksummed); `refs` are large external files `{path, md5, bytes, role}` (checked for
existence + md5). Layers carry no packing rules (those are cell config).

| kind | content (Jurkat em_cell1_v7) |
|---|---|
| `mesh_midsurface` | NE HMFF midsurface PLY (label frame nm; model = minus `frame.centre_label_nm`) |
| `cytoskeleton_tracks` | MT centrelines npz (`points_model_nm`, `offsets`, `flipped_for_lattice`, classes) |
| `centrosome` | centriole-pair centre, MTOC candidate (json) |
| `npc` | pore positions/counts (EM vs random) -> `pores.npz` |
| `point_priors` | ribosome instance centroids (ref, s0 voxels) |
| `microvilli` | finger density / lengths summary |
| `density_map`, `chromatin_territories` | EM density grids, territories (not in em_cell1_v7) |
| `reference_geometry` | superseded geometry kept for provenance (ref) |

## 9. Frame and units

nm, right-handed, model origin declared (`frame.origin`, Jurkat: bounding-box centre of the refined cell mask,
`centre_label_nm`). `frame.source_affines`: 4x4 maps from source index/frames to model nm, e.g.
`seg_s0_index_to_model`, `pred_s0_index_to_model` (OpenOrganelle `_pred` y-flip: y_seg = 2999 - y_pred at s0),
`raw_s0_index_to_model` (raw = label + (0, 9, -2) s0 voxels), `label_frame_nm_to_model`.

## 10. Validation (`cellgeom.validate`, run at creation and by the geometry stage)

| check | pass criterion |
|---|---|
| checksums | every file matches size + sha256 |
| schema | unique ids / engine codes; codes in 0..63; core ids match the table; extensions have a core parent; rule names known labels |
| label_rule | declarative rule == engine labels on 400k points (200k uniform, 200k within 50 nm of membranes): agreement 1.0 |
| eikonal | per field, median \|∇d\| (central differences h/2) within [0.95, 1.05] inside the band |
| nesting | >= 99 % of each organelle field's interior inside the cell field |
| volumes_mc | label volumes by seeded MC (reported) |
| field_volumes_areas | MC field volume / shell (\|d\| < 4 nm) area vs `summary` (vol within 4σ + 1 %, area within 5 % + 4σ) |
| clip_planes | every point beyond a plane has its declared label |
| pores | \|ne_mid\| at centres < h/4, axis vs ∇ne_mid median < 10°, spacing >= `min_sep_nm` |
| layers | depends_on == package_id, files present, refs exist with matching md5 |
| (toolkit extras) | e.g. `b1_legacy_equality` for the Jurkat conversion |

`validation_waivers` turn named failures into documented, reported waivers (`waived_failures`); they exist only for
converted legacy geometry whose defects must be reproduced exactly.

## 11. Consumers

`load(dir) -> Cell` (engine Cell subclass): `.G`, `.pores`, `.pidx/.PT/.PQ/.PC`, `.labels(P)`, `.mems(P)`,
`.mems64(P)`, `.manifest`, `.manifest_sha256`, `.vocab`, `.code_to_name`, `.name_to_code`, `.code_to_id`,
`.mem_code(name)`, `.summary`, `.layer_path(name)`, `.engine`, `.package_dir`, `.package_id`.

## 12. Deviations from proposal (e) and why

| proposal | v1 | why |
|---|---|---|
| core ids 0-31, extensions >= 32 | core < 1000, extensions >= 1000, engine codes separate (0..63) | no 32-region cap anywhere; packer bitmask limit lives only in the per-package engine mapping |
| label rule compiled to a generic interpreter replacing `label_from` | v1: declarative metadata. Since 2026-11-16 (A1): table-driven evaluator `cellgeom_eval_v1` with a repo rule table per layout (§4.3) | the rule table reproduces `sdf_sections_v1` exactly (B1 level E), so the hand-coded evaluator could be replaced |
| G sections become prims | kept as typed sections; prim table optional | bit-identical reproduction; no prim engine yet |
| header[12] with unused slots | slots 10/11 = relative [index, data] offsets | self-consistent header, rebased at load |
| data float32 | float32 when lossless, else float64 | cristae id grid / exactness |
| meshes + LODs, dense label grid required | optional | not needed by the packer; viewer encoder can derive them |
| validation: label grid vs EM >= 99.5 %, membrane separation, mesh checks | not in v1 | needs the source label grid; listed as v1.1 work |
| (not in proposal) | `validation_waivers`, legacy summary verbatim | converted legacy packages must reproduce known defects transparently |
