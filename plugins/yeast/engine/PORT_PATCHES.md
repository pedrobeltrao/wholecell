# Vendored yeast code: origin and port patches

Origin: `git archive e3b7448` of the yeast development repository (branch round-x; worktree
`{yeast_root}/dev_polysomes`), extracted 2026-11-16; identical to `{yeast_root}/runs/X5_seed0/tree` for pipeline/,
scripts/, validation/, simgpu/, analysis/, web/, env/ (docs/figures differ: they are outputs of X5's validation).
Checksums of the extracted files: `ORIGIN_sha256.json` (before the header lines and patches below).
Every .py/.sh file carries a first (or second, after a shebang) line `# origin: yeast yeast_cell_dev_bundle/<path> @ e3b7448`.
**Licence (release staging 2026-11-23):** G. Magnusson agreed to the inclusion of this code in the wholecell repository
under the repository licence (PolyForm Noncommercial 1.0.0, `LICENSE`; see `NOTICE.md`). The words "no licence yet;
local only" in the `# origin:` header lines are historical (kept unchanged as provenance) and no longer apply.

Port patches (the only code changes; each marked in the file):

| id | file | change | why |
|---|---|---|---|
| P1 | scripts/s2_geometry.py | if `YC_GEOM_PACKAGE` is set, install `native/geometry.json` + `native/geometry_fields.npz` of that cellgeom_v1 package into `W3D` and load it, instead of calling `cellshape.build` / `cell3d.build` here | one geometry package per build (toolkits/analytic/yeast.py runs the same builder calls); unset = original behaviour |
| P2 | scripts/validate_all.sh | `sweeps_pooled.py --gpus ${YC_NGPU:-4}` | the script assumed 4 GPUs (jobs on devices 1-3 fail on a 1-GPU node) |
| P3 | scripts/validate_all.sh | PDF step `md_to_pdf.py ... \|\| echo` (non-fatal) | needs the bundle's headless-browser libraries (.chromelibs), not available; the PDF is not a result |
| P4 | pipeline/amphiphysin.py (`place`) | Rvs161-Rvs167 frame `W = np.c_[around, -u, -radial]` (was `[around, u, -radial]`) and library frame mid axis = short x long; **changes output** (2026-11-17, see below) | the old matrix was an improper rotation (det -1 = mirror image) |
| P5 | analysis/analyze.py, analysis/rcsb_coverage.py | author-machine absolute paths (`/Users/...`) replaced by env variables `YC_BUNDLE`, `YC_VIEWER_TYPES`, `YC_AFDB_DIMERS` (defaults relative to the working directory) | release staging 2026-11-23: no personal paths in the public tree; one-off analysis scripts, not used by the pipeline (no output change) |
| P6 | web/app/package.json, web/app/package-lock.json | `"license": "ISC"` (npm default) -> `"SEE LICENSE IN ../../../../../LICENSE"` | release staging 2026-11-23: the npm default contradicted the repository licence (no code change) |
| - | docs/model_validation_and_tuning.pdf | not included in the release tree (3.2 MB Chromium print of `model_validation_and_tuning.md`, same text and figures; rebuilt by `scripts/md_to_pdf.py`). `ORIGIN_sha256.json` still lists it (origin record, unchanged) | size; duplicate of the .md |

### P4 (2026-11-17): proper rotation for the Rvs161-Rvs167 amphiphysin placement - CHANGES OUTPUT

What: `pipeline/amphiphysin.py`, `place()` (e3b7448 lines 56 and 106; marked `PORT PATCH 4`). Each dimer's orientation
was `R = W @ Rl.T` with `W = np.c_[around, u, -radial]` and `around = np.cross(u, radial)`: W is left-handed (det -1).
`Rl` (the library entry's principal axes from `eigh`) is right-handed for RVS161_RVS167_het, so det R = -1 for every
dimer: an improper rotation, i.e. a mirror image. Consequences in the original code (port build s0 = the X5 reference):
all 303 placed dimers were affected - the residue beads tested for dimer-dimer clashes and membrane containment
(`Bz = (z - c0) @ R.T + ctr`) were the mirror image of the AlphaFold heterodimer, and the stored orientation
(`rvs_q`, used for the packing obstacles, the export and the viewer) was `Rotation.from_matrix(R)` of the improper
matrix, which scipy 1.15 silently turns into an unrelated proper rotation (median 142 deg from the intended
orientation; long axis around the tube and concave face on the membrane lost; stored dimers reach up to 2.8 nm into
the bilayer). The stored dimers themselves are proper (not mirrored); only their orientation is wrong.
Found because scipy 1.18 (core env) rejects a det -1 matrix in `Rotation.from_matrix` (unification A6).

Fix: `W = np.c_[around, -u, -radial]` (det +1) and `ax_mid = np.cross(ax_short, ax_long)` (library frame right-handed
whatever sign `eigh` returns). The library mid axis has no placement meaning (only long axis = around the tube and short
axis = concave side towards the membrane are intended, both unchanged), so flipping its world direction is free.
`Rl` is numerically unchanged for this entry (det was already +1), so the lattice (dimer length, width, thickness) is
unchanged. Check: `scripts/dev/check_amphiphysin_chirality.py`, results in
`data/derived/wholecell_builds/yeast/_dev_unify/amphiphysin_fix/`.

Changes output: the amphiphysin layer (`rvs_pos`, `rvs_q`; the count can change because the residue-level tests now
see the real dimer) and everything after it that depends on it: the later serial pre-pass builders (they avoid the Rvs
dimers' sphere trees and read the pools after the Rvs deduction) and the packing. The reference for the serial pre-pass
is now "port s0 + patch 4" (`_dev_unify/amphiphysin_fix/ref_yeastenv_patched`), not s0 itself (builders 1-18 are
unaffected and still identical to s0). Measured (yeast env, seed 0): amphiphysin 303 -> 293 dimers (all proper, long
axis exactly around the tube, concave face on the membrane, 2.3 nm clear of the bilayer); actin patches 45954 -> 45915
instances; polysomes 182329 -> 181965 ribosomes; RAC, mRNA threads and mRNP ends change accordingly; ATP synthase,
mtDNA, Cvt and NVJ layers identical to s0.

TO TELL THE YEAST TEAM (G. Magnusson): `pipeline/amphiphysin.py`, `place()`, the line `W = np.c_[around, u, -radial]`
(e3b7448 line 106) builds a left-handed frame, so `R = W @ Rl.T` has det -1 (mirror image of the Rvs161-Rvs167
dimer). With scipy <= 1.15 `Rotation.from_matrix(R)` (line 110) silently returns an unrelated rotation, with scipy 1.18
it raises "Non-positive determinant". Fix: `W = np.c_[around, -u, -radial]`, and make the library frame right-handed
(`ax_mid = np.cross(ax_short, ax_long)` before `Rl = np.c_[ax_long, ax_mid, ax_short]`). A scan of the other builders
(pipeline/*.py, scripts/*.py) found no other improper construction: every Kabsch has the det sign fix, every other
frame is `[a, b x a, b]`-type (right-handed) or checked (atp_dimers flips a column if det < 0).

## Patch 5 (2026-11-17): `pipeline/mrna_genes.py` gene_table(): missing gene names under pandas >= 3

Environment compatibility, NO output change against the yeast environment. `gene=p.iloc[:, 1].astype(str)` turns a
missing standard name into the string "nan" under pandas 2 (yeast env), but keeps it as NaN under pandas 3 (core env,
decision 2026-11-17 "one core environment"), so `g.gene.isin(["nan", ""])` missed it and genes without a standard name
(YDR154C, YAL037C-A, ...) were labelled "nan" instead of their ORF: in the core-env serial pre-pass, `mrna_gene`
9,338 / 90,243 and `poly_chain_gene` 3,760 / 39,314 labels (coordinates and counts unaffected). Fix: also test
`g.gene.isna()`. Check: `gene_table()` written by the yeast env (pandas 2) and the core env (pandas 3) with the patch:
byte-identical TSV, 5,083 genes, 0 named "nan". TO TELL THE YEAST TEAM (G. Magnusson): the same line breaks under
pandas 3; the one-line fix is above.

Not patched (handled by the environment, plugins/yeast/runner.py): paths (`YC_WORK`, `YC_LIB`, `YC_PDB`, `YC_INPUTS`),
seed and workers (`YC_SEED`, `YC_NPROC`), `PYTHONHASHSEED=0`, one BLAS thread, numba cache per build
(`NUMBA_CACHE_DIR`; run_all.sh deleted the in-tree caches), GPU lock per build (`YC_GPU_LOCK`), one GPU
(`CUDA_VISIBLE_DEVICES=0`). `.venv` and `pipeline/data` are links created at build time (ignored by git).
