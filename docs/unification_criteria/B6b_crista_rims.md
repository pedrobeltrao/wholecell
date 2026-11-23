# B6b rounded crista rims (Jurkat): acceptance criteria (frozen before any B6b code or package)

Follow-up of B4 criterion J6 (`B4_preplaced.md`, md5 e715c2f9…): the common ATP-synthase dimer builder
(`core.atp_synthase_dimers`, rims `lamellae`) places 0 dimers on the cristae of `jurkat/em_cell1_v8p` at the frozen
2.5 nm TM tolerance, because the lamellar cristae (`crista_lumen` of the evaluator, section `cristae_lamellae_v1`,
made by `toolkits/em/jcristae.py`) end in a flat cap with square corners. In cells the crista rim is a tightly curved
ridge on which the dimer rows sit. MAIN-SESSION DECISION 2026-11-17 (Pedro informed): round the rims in the geometry if
the change is contained. Written 2026-11-17 before any B6b code or package existed; the only B6b-era numbers used are
measurements of the existing package `em_cell1_v8p` (listed below with their source). The md5 of this file is
reported to the coordinator. Later findings go to a dated addendum file, never into this one.

## Feasibility (decided before freezing)

Contained. The lumen of one lamella is `lam = max(slab, b + CR_GAP)` with `slab = |u| - CR_W` (u: distance to the
lamella mid-plane, b: matrix SDF of the boundary inner membrane), unioned with the junction ribbon. Rounding is one
rounded intersection in this one kernel (`crista_lumen`, three copies: `cellgeom_eval_v1`, legacy `sdf_sections_v1`,
python reference `geometry/eval_ref.py`) controlled by a NEW engine parameter; the package step changes only
`manifest.engine.params`. No section, field, rule-table or toolkit-geometry change. Estimate: ~0.5 day incl. checks.

## What changes (fixed now)

| Item | Specification |
|---|---|
| Engine parameter | `CR_RIM` (nm), header slot 40 of G (reserved range 26-62, unused by MEM_HALF ≤ 16 values at 18-33), in `PARAM_SLOTS` of both engines; NOT added to `DEFAULT_PARAMS` (packages without it get 0.0, so their manifests and G are unchanged). `CR_RIM = 0` = the current square cap, same float expressions (bitwise). |
| Kernel | `CR_RIM > 0`: `lam = min(-r, max(slab, e)) + sqrt(max(slab + r, 0)^2 + max(e + r, 0)^2)` with `e = b + CR_GAP`, `r = min(CR_RIM, CR_W)` (rounded intersection; the apex stays at `b = -CR_GAP`, the rim becomes a half-cylinder of radius r around the rim axis `b = -CR_GAP - r` in the mid-plane). The junction ribbon term and its corners are unchanged. |
| Value | **CR_RIM = 10.0 nm = CR_W**: the IM midsurface separation across the crista lumen is 20 nm (CR_W = 10; inter-cristae membrane distance 20 ± 4 nm, Davies et al. 2011 PNAS 108:14121; Mannella 2006), so a lamella whose lumen width is kept ends in a half-cylinder of midsurface radius 10 nm. This is the largest rounding possible without widening the lumen. Checks: the yeast dimer geometry used by the builder puts the dimer's membrane centres on a circle of 8.8 nm (Davies et al. 2012 PNAS 109:13602); dimer ribbons impose a "17-nm outer radius" (Strauss et al. 2008 EMBO J 27:1154; read as the outer membrane surface ≈ 15 nm midsurface, an upper bound, since a larger radius needs a bulged rim = wider lumen, out of scope). |
| Package | NEW `jurkat/em_cell1_v8r` from `em_cell1_v8p` by a `toolkits/em` package step (plain Python, `pkgstep.open_pkg` / `publish`) that sets `engine.params.CR_RIM = 10.0` and nothing else. `em_cell1_v8p` and all other packages are never written. |
| Config | `config/cells/jurkat.yaml` `geometry.package` -> `em_cell1_v8r` (+ manifest sha256); frozen pre-placed modules' `depends_on_geometry` -> `jurkat/em_cell1_v8r` (justified by R5: labels outside the mitochondria identical). The `atp_synthase` builder params return to the frozen B4/yeast values: `tm_tol_nm` 2.5, `dimer_angle_deg` 84, `rim_radius_nm` 8.8, max in-membrane 2 %, spacing 13 nm; only the rim search list `rho_nm` and `samples` may be set. |

## Pre-measurements on `em_cell1_v8p` (scratch script, 4,000,000 points uniform in the 699,279 cristae id-grid cells with an enabled lamella, 77.3 µm³, seed 0)

Rim band volume (|u| < CR_W, -CR_W < b + CR_GAP < 0, outside the junction ribbons) 0.806 µm³ -> rim line length
**L ≈ 4,031 µm** (band / (2 CR_W²)); matrix 28.94 µm³, IMS 11.10 µm³; IM area (|IM| < 1 nm shell / 2 nm) ≈ 1,265 µm².
Analytic expectations for r = 10 nm: volume moved from IMS (crista lumen) to matrix = L (2 - π/2) r² = 42.9 nm² per nm
of rim -> **+0.173 µm³ matrix / -0.173 µm³ IMS**; IM area change = -L (4 - π) r = -8.58 nm² per nm -> **-34.6 µm²
(-2.7 %)**. A full dimer row at the builder spacing (13 nm) on L would hold ~310,000 dimers (620,000 of the 868,688
CPX-6151 monomers).

## Criteria

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| R1 | Old packages unchanged | `CR_RIM` absent = 0: G of `em_cell1_v7`, `em_cell1_v8`, `em_cell1_v8p` built by the new code; label and membrane SDFs on 1,000,000 seeded points (core evaluator) vs the code at e03fffe; B1 check (`tests/check_b1_geometry.py --package {geom_root}/jurkat/em_cell1_v7`) | G identical (bitwise) to e03fffe; labels / SDFs bitwise identical; B1 passes |
| R2 | New package = parent + one parameter | `em_cell1_v8r` vs `em_cell1_v8p` | every content file (fields, sections, pores, meshes, layers) same sha256; manifest differs only in `engine.params.CR_RIM` (10.0), generator, package_id / geometry_id, layers' `depends_on` / `carried_from`, description, summary.cristae (rim note); v8p folder untouched (manifest sha256 4bfc54cb… unchanged) |
| R3 | Validator | `cellgeom.validate` (written by `finalize`) | `pass` true, the same checks and waiver as v8p (only `field_volumes_areas` legacy summary); `label_rule` agreement 1.0 |
| R4 | MC volumes (validator, 2,000,000 points, same seed / stream as v8p) | `volumes_mc` v8p vs v8r | `mito_matrix` up and `mito_ims` down by the same count of points; every other label volume **identical** (same points, same labels) |
| R5 | Labels change only at crista rims | 20,000,000 points uniform in the enabled cristae id cells (as the pre-measurement, seed 1) + 2,000,000 uniform points over the cell box (seed 2), v8p vs v8r | every changed point goes `mito_ims` -> `mito_matrix`, and lies in the rim corner zone \|u\| ≤ CR_W + 0.5, -CR_W - 0.5 ≤ b + CR_GAP ≤ 0.5 (within CR_W·√2 + 0.5 = 14.6 nm of the rim line); 0 changed points elsewhere |
| R6 | Volumes moved (paired, the 20 M points of R5) | ΔV = changed points × cell volume | **0.10-0.25 µm³** (expected 0.173; i.e. ≤ 2.3 % of IMS, ≤ 0.9 % of matrix); matrix + IMS identical |
| R7 | Other membranes unchanged | all membrane SDFs at the 20 M + 2 M points | PM, NE, ER, MITO_OM, GOLGI, ENDO, CELL_OUTER bitwise identical; MITO_IM identical outside the R5 zone |
| R8 | IM area | shell estimate (\|IM\| < 1 nm, paired points of R5) v8p vs v8r | ΔA between **-20 and -50 µm²** (expected -34.6) |
| R9 | ATP synthase dimers on the rims, frozen tolerance | builder `core.atp_synthase_dimers` (rims lamellae) on v8r with the B4 / yeast dimer (84°, 8.8 nm, TM tol **2.5 nm**, ≤ 2 % in membrane, spacing 13 nm), seed `preplaced.atp_synthase`, config `samples` | dimers > 0 and monomers used ≤ pool (868,688) / 2 (B4 J6, unchanged); count, sites tried, rho used reported |
| R10 | Density per µm of rim vs literature | the builder saturated on a sub-region (≥ 3 mitochondria with cristae, rim points ≥ 3 per 13 nm of rim), dimers / rim length of those lamellae (band method) | **38-83 dimers per µm** (literature: dimer spacing along the row ~12 nm, Strauss 2008; 12 to > 20 nm, Davies 2011 -> 50-83 per µm; lower bound 38 = half of a full 13 nm row). Cell-wide: count at the config sampling / L reported as rim coverage (sampling-limited, not a geometry criterion) |
| R11 | Tests and dry run | `tests/test_geometry.py`, `tests/test_preplaced_b4.py` as scripts; snakemake dry run of the Jurkat build | pass; dry run resolves with the new package |

No full build is part of B6b. Build-level acceptance (suite no_regression) is left to the next Jurkat build; EXPECTED
there: mito_matrix / mito_IMS densities change by ≤ 2.5 % (R6), IM-band metrics by the R8 area change, ATP synthase
monomers now placed as dimers at the frozen tolerance.
