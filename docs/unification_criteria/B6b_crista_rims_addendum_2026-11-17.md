# B6b rounded crista rims: addendum 2026-11-17 (findings after the frozen criteria)

Frozen criteria: `B6b_crista_rims.md` (md5 da317056…, commit 1d7b288). Package `jurkat/em_cell1_v8r` (manifest
sha256 6c86b168…), code 3254949. Measurements in `{results_root}/jurkat/_dev_unify/b6b/` (`geom_check.json`,
`atp_<package>_tol<t>.json`, `b1_v7.json`; `scripts/b6b_crista_rims_check.py`), 3 cores.

## Results

| # | Result | Verdict |
|---|---|---|
| R1 | G, labels and all membrane SDFs of em_cell1_v7 / v8 / v8p on 1,000,000 points identical (bitwise) to the code at e03fffe; B1 check on v7 PASS (2.55 M points, SDFs exactly equal) | PASS |
| R2 | 15/15 content files same sha256; engine params differ only in CR_RIM (absent -> 10.0); v8p manifest sha256 4bfc54cb… unchanged | PASS |
| R3 | validator pass; failed = waived = [field_volumes_areas] as v8p; label_rule 1.0 | PASS |
| R4 | validator MC (2 M): mito_matrix +0.168 µm³, mito_ims -0.168 µm³; every other label volume identical | PASS |
| R5 | 20 M points in the cristae id cells: 45,276 changed, all mito_ims -> mito_matrix, 0 outside the rim corner zone (max 14.12 nm from the rim axis); 2 M box points: 38 changed, all IMS -> matrix, in the zone | PASS |
| R6 | volume moved 0.175 µm³ (expected 0.173; rim length 4,040 µm) | PASS |
| R7 | PM, NE, ER, MITO_OM, GOLGI, ENDO, CELL_OUTER bitwise identical everywhere. MITO_IM: **values** also change outside the R5 zone (2.9 % of the 20 M points) | NOT PASSED as written, see below |
| R8 | IM area 1,264.9 -> 1,229.2 µm² (-35.7; expected -34.7) | PASS |
| R9 | dimers at tm_tol 2.5 nm: **0** (cell-wide, 400k samples) | **FAIL** |
| R10 | saturated density on 3 mitochondria (30.5 µm of rim, 4.2 rim points per 13 nm): 0.1 / µm at 2.5 nm | **FAIL** |
| R11 | test_geometry, test_preplaced_b4 (scripts) pass; snakemake dry run (cell=jurkat) resolves, 36 jobs; config loads em_cell1_v8r with the declared sha256 | PASS |

**R7 (criterion wording error).** The frozen text required MITO_IM to be identical outside the zone where labels can
change. That cannot hold: rounding a corner changes the *distance* to the membrane on the matrix side next to the
corner, without changing any label. Measured (4 M points): every IM value change outside the zone is on the matrix
side (IM > 0 before and after); |u| ≤ 32 nm, b + CR_GAP in [-10, +6] nm; changes > 0.01 nm lie ≤ 31.2 nm from the
rim axis; |Δ| median 0.79 nm, p99 4.8 nm, max 5.5 nm. This is the intended geometry change (the packer clearance
to the IM near the rims now follows the round rim). No other membrane changes.

## R9 / R10: the rounded rim is necessary but not sufficient for J6 at 2.5 nm

Dimer counts with the frozen yeast dimer (84°, 8.8 nm), ≤ 2 % in membrane, spacing 13 nm, config sampling
(400,000 samples, ~3,950 rim points) and saturated on the same 3 mitochondria (30.5 µm of rim):

| TM tolerance | v8r (round) cell-wide | v8r saturated / µm | v8p (square) cell-wide | v8p saturated / µm |
|---|---|---|---|---|
| 2.5 nm (frozen J6) | 0 | 0.1 | 0 | 0 |
| 3.0 nm | 14 | 1.1 | 0 | 0 |
| 3.5 nm | 674 | 13.5 | 1 | 0 |
| 4.5 nm | 1,809 | 26.2 | 21 | 1.2 |

(B4 placed 1,451 on v8p only with a non-yeast dimer opened to 14.9 nm at 4.5 nm.)
Literature: 50-83 dimers per µm of rim (spacing ~12 nm, Strauss 2008; 12 to > 20 nm, Davies 2011).

Why 2.5 nm cannot be met on any smooth rim: the builder's TM residues are the monomer beads with |z| < 2.0 nm of the
library membrane plane (906 beads, lateral extent 9.4 × 9.8 nm). A rigid flat TM patch that wide on a curved
membrane leaves its edges off the midsurface. Fitting the dimer's TM beads to an ideal cylinder (any axis
offset) gives a best maximum deviation of 2.5-2.7 nm for radii 10-17 nm and dimer angles 84-180°. Real rim sites
add ~0.3-0.5 nm (b is not an exact distance; the rim line curves within the plane), so best sites sit at 3.0-3.3 nm.
In cells the membrane bends sharply at the dimer interface (the dimer makes the kink) and is flat under each
monomer. A smooth half-cylinder is a different shape.

The cell-wide count is limited by sampling, not geometry: 400,000 samples over the whole id-grid box give ~3,950 rim
points for a 4,040 µm rim, which is 1 per µm. A full row would need ≥ 310 per µm.

## Open items (for the main session; not decided here)

1. **J6 tolerance (DECISION NEEDED).** Fixing the geometry does not make 2.5 nm reachable for this dimer model. Options:
   (a) evaluate J6 at 4.5 nm on the rounded rims (1,809 dimers at the config sampling; 26 per µm saturated), or
   (b) change the test, not the geometry: a 2.5 nm limit on the *median* TM deviation (best sites: median 1.3-1.4 nm),
   or TM = |z| < 1.5 nm, or
   (c) a kinked rim (two flat faces meeting at the dimer angle) instead of a half-cylinder. That is a new rim shape in
   `crista_lumen`: about one more day of work.
   `config/cells/jurkat.yaml` keeps tm_tol 4.5 nm, labelled DECISION NEEDED (same status as after B4).
2. **Rim coverage.** To approach the literature density cell-wide, the builder should draw its samples only in the
   enabled cristae id cells, not in the whole bounding box, and use more of them. This is a builder change (B4 code)
   and affects runtime. Not done.
3. The R7 wording: future geometry criteria should state label identity and *other-membrane* identity, plus a
   bounded neighbourhood for the changed membrane's SDF.
