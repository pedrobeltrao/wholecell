# Stage 1: geometry

**Purpose.** Give the build one immutable cell geometry: membranes as signed distance fields, compartment labels,
nuclear pores and optional layers (EM traces, meshes). Every later stage asks the geometry "which label is this
point in" and "how far is the nearest membrane".

**Inputs.** A cellgeom_v1 package (`{geom_root}/<cell>/<geometry_id>/`, spec
[../reference/cellgeom_v1.md](../reference/cellgeom_v1.md)). The package is made by a **toolkit**, outside the build
for Jurkat and testcell:

- `toolkits/em` (Jurkat): `convert_legacy.py` wraps the published v7 geometry as package `em_cell1_v7` (the reference
  preset `config/cells/jurkat_v7_port.yaml`); the package step `jne_fix.py` writes **`em_cell1_v8`** (B6), the package
  of `config/cells/jurkat.yaml`. The generator from EM segmentations is ported but not run end to end.
- `toolkits/analytic` (testcell, yeast): `generate.py` builds a package from analytic shapes
  (`toolkits/analytic/specs/testcell.yaml`); `yeast.py` builds the yeast package.

**Outputs** (`geometry/`): `package.json` (path, id, manifest sha256, engine), `validation.json` (validator result),
`measures.json` (region volumes by seeded Monte Carlo, stream `geometry.measures`; membrane areas; field volumes).

**Config keys** (`geometry:`): `toolkit`, `package`, `manifest_sha256`, `engine`, `evaluator`, `rules`, `mc_points`
(every cell, yeast included since A9; the reproduction config `yeast_x5_port.yaml` keeps `generator`, `package_root`,
`seeds`).

**Checks.** Manifest sha256 equals the config value; engine equals the manifest's engine; package validator
(cellgeom_v1 §10; reused when fresh). Jurkat criterion B1 (v7 reproduction): labels and SDFs identical to the legacy
evaluator on 2,000,000 points.

**Code.** `wholecell/geometry/{stage.py, cellgeom.py, engine.py, engines/}`; `workflow/rules/geometry.smk`.

**Jurkat vs yeast today.** Jurkat and testcell use engine `sdf_sections_v1`. Yeast uses `yeast_cell3d_v1`
(`plugins/yeast/geom_engine.py`, the vendored yeast evaluator with analytic primitives), and its package is generated
inside the build by the engine step `geometry_package`. The active engine is chosen by the env var
`WHOLECELL_GEOM_ENGINE`. **A1** replaces both with one table-driven evaluator (`cellgeom_eval_v1`) and removes the env
var. Done for Jurkat/testcell (2026-11-16) and, as an evaluator, for yeast (2026-11-17): `cellgeom.load` evaluates the
yeast package with `config/geometry/rules/yeast_cell3d_v1.yaml` (prims, 10 membranes), identical to the vendored cell3d
(`tests/check_a1_yeast_geometry.py`). Since A9 the yeast config names its package like the other cells
(`analytic_x5_s0_f5f3a5c723b2`, written by `toolkits/analytic/yeast.py` in the port build s0) and the core stage loads
it with the core evaluator; only the reproduction config `yeast_x5_port.yaml` still generates a package per build.

**B6 (2026-11-17; criteria `docs/unification_criteria/B6_geometry.md`).** Two fixes that change outputs, each a new
package (the old ones are untouched and stay selectable):

- **Jurkat `em_cell1_v8`** (`toolkits/em/jne_fix.py`): `ne_mid` is now the true signed distance to the HMFF NE midsurface
  on both sides. In v7 every grid node outside the mesh held the 120 nm clamp: `jgeom_build.mesh_sdf` negated a uint8
  occupancy (`~occ` = 254/255), so the outside never got exact distances. The ONM therefore sat ~8 nm inside the
  midsurface. MC volumes v7 -> v8: perinuclear 2.744 -> 8.025 µm³ (expected 8.08 = 245.1 µm² × 34 nm minus 1,000 pore
  openings), cytosol 392.81 -> 387.79, ER lumen 17.95 -> 17.61, pore channel 0.352 -> 0.407; nucleoplasm and everything
  away from the NE unchanged. Pores, the mesh and the other fields are byte-identical. Frozen pre-placed modules now
  declare `depends_on_geometry` (checked by the preplaced manifest step and suite C16). Pre-placed centres that fall in
  the new NE lumen are dropped at indexing (`pack.index_drops` `ne_lumen_b6`). Reference preset:
  `config/cells/jurkat_v7_port.yaml`.
- **Yeast `analytic_x5_s0_f5f3a5c723b2_memhalf1`** (`toolkits/analytic/yeast_mem_half.py`): per-membrane half thickness
  = D_P-P/2 of Pogozheva et al. 2022 SI Table S20 (`config/geometry/membrane_half_thickness.tsv`: PM 2.11, ER/NE 1.895,
  Golgi 1.98, mito OM 1.975, IM 1.955 nm; vacuole, peroxisome, LD and wall keep 2.0). The rule table's clearance is
  `thickness_corrected`, and a uniform package compiles to the legacy plain rule (bitwise). Labels, SDFs and volumes
  are identical; the exclusion boundaries move by ≤ 0.11 nm. Reference preset: `config/cells/yeast_x5_core.yaml`.

**Methods text (draft).** The cell geometry was stored as an immutable package of sparse signed-distance fields, one
per membrane, with a declarative rule that assigns a compartment label to any point. Packages were generated from
segmented electron microscopy volumes (Jurkat) or from analytic shapes (yeast, synthetic test cell) and validated
before use (checksums, label coverage, distance-field gradient). Compartment volumes were estimated by Monte Carlo
sampling with a recorded seed.
