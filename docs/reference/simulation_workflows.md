# Simulation workflows: applications that consume pipeline builds

Status: spec v1, 2026-11-15 (SPECS agent). Scope: how the simulations of both projects become reproducible
workflows reading **exports of a pipeline build** (ARCHITECTURE §3 `export/`). Audit references:
H = `data/derived/repro_audit/human/`, Y = `data/derived/repro_audit/yeast/`.

## 1. Principles

| Rule | Why (audit) |
|---|---|
| One Snakemake workflow per application under `applications/<name>/`, separate from the build Snakefile | human simulations are shell chains + manual launches (H merge_notes §1, stages §10); yeast `simgpu` has "no CLI orchestration" (Y validation_sims_notes §4) |
| An application reads only a **pinned build** (`build_id` + sha256 of its `export/MANIFEST.sha256`), never "the current pack" | human `wc/` = live pack; BD boxes cut from v4 then re-cut from v7 (H formats §6); yeast vesicle sims ran on a dirty mixed-commit build (Y status §3) |
| Boxes come from one exporter and one schema (§3), not per-application extractors | 7 human extractors (H stages §9) vs one yeast schema (Y comparison §(d): yeast better) |
| Pre-registration and held-out discipline are workflow rules (§4), not conventions | human strength (H merge_notes §2) but checked by hand; yeast `splits.py` frozen but unchecked |
| No pickles, no hand-written parameter files, no hard-coded PIDs or paths | BD `lib_<level>.npz` pickled dict; `theta_final.json` by hand; d2 waits on a PID (H stages §10, formats §6) |
| Legacy results are never moved or rewritten; re-runs are compared to them (§7) | provenance policy (paper planned, MEMORY "paper provenance") |

## 2. Layout and application config

```
applications/<name>/
  Snakefile                  includes rules/*.smk; `snakemake -d results/<run_id> --configfile config.yaml`
  config.yaml                schema below
  preregistration/           <date>_<test>.md + <date>_<test>.md5 + prereg.yaml (frozen; §4)
  legacy.yaml                archival record: legacy result files + md5 (read-only; §7)
  rules/*.smk  scripts/      thin entry points; engine code lives in wholecell/sim/ or plugins/<cell>/sim/
  envs/*.yaml                conda spec + `conda list --explicit` lock per engine
  results/<run_id>/          write-once (run_id = <date>-<app commit7>-<build_id short>-s<seed>)
```

```yaml
schema_version: 1
application: bd_diffusion
build: {cell: jurkat, build_id: 20261116-abc1234-s20261115, export_manifest_sha256: "..."}   # resolved via site.yaml results_root
seeds: {master: 7}                       # named streams via wholecell.seeds (same scheme as builds)
boxes: [{name: cyto_500, centre_nm: [..], edge_nm: 500, boundary: periodic, walls: false}, ...]
frozen_inputs: {bd_theta: {path: ..., md5: ..., origin: rule bd_v2.fit @ <run>, why: ...}}
gpu: {required: true, profile: slurm}    # §6
deposit: [run_bd]                        # non-deterministic rules (§5)
tolerances: {D_over_D0: {rel: 0.05}, heldout_mean_abs_log10: {abs: 0.02}}
```

Standard rule names (so a reader knows every workflow): `check_build` → `boxes` → `setup` → `run` → `analyse` →
`prereg_check` → `evaluate` → `report`; optional `fit`, `posthoc`, `rerun_compare`, `compare_legacy`.
`check_build` verifies the build stamp, `export/MANIFEST.sha256` and that the build's `regions.json` contains every
region name the application config uses (an application declares its cell requirements; no cell list in code).

### 2a. Application catalogue

| Application | Cells | Arms / tiers | Inputs from the build | Engine (env) | GPU | Legacy record (prereg md5) |
|---|---|---|---|---|---|---|
| `bd_diffusion` | any | v1 benchmark, v2 tiers 0-2 + lumen §6 (`*_v7sh`), v3 (audit), yeast crowding sweeps / encounters / Gade | boxes, library bodies, densities | CuPy BD (`bd`), yeast `simgpu` (RPY) | yes | BD v2 held-out 409929f8…, lock 2026-11-11; BD v3 671895d3… (failed) |
| `synapse` | jurkat | S0, S1, S2.0, S2 rerun, S2.1, crispr, glycocalyx (CD43 + ranking), external_evidence | PM surface mesh/field crop, membrane inventory (`species_loc`), ledger, BD theta | CPU (vcell) | no | S0 f84599c8…, S1 017b1ec7…, S2.0 be8f17b6…, S2.1 S21_freeze.md5, crispr f6090676…, glyco 8fb8dcac… |
| `condensates` | any (screen); jurkat (in-situ) | screen, d2 ladder/level-2, refit, residue (CALVADOS/Mpipi), in-situ nucleolus | inventory (species, copies, sequences), nascent RNA, chromatin crop, label grid | OpenMM 8.6 (`hic`), CALVADOS 0.9.1 (`calvados`), HOOMD 5.2 (`insitu`) | yes | 77d290fe…, 3b08bb5c…, in-situ ff22ce71… (parked) |
| `vesicle_traffic` | yeast | membrane (stage 0), release (1), transit (2), T3 knock-outs | endocytic-site crop, box, geometry | CuPy `mesobd` + CPU Helfrich solver | stage 1 | none yet (prototype on a dirty build, Y status §3): pre-register before the re-run |
| `nfat_calcium` | jurkat | A1, A2 + controls/sensitivity | geometry (NPCs, nucleus), ledger, BD theta, chromatin | CuPy (`bd`) | yes | test_plan 42e4e8ad… |
| `chromatin_validation` | any | human v8 phase 1 (EM resolution, PMA Hi-C, pilot); yeast Hi-C / nucleosome / DNase tests | chromatin pre-pass, EM crop, Hi-C | OpenMM (`hic`) | pilot | PREREG 94bdd94b… (H1-H4, E2 failed) |
| `yeast_condensate_test` | yeast | mean-field calibration + test | inventory, features | CPU | no | `calibration_nucleus.json` (frozen ratios) |

## 3. Simulation-box export schema (`wholecell.simbox/1`)

Written by the build's export stage (`export/boxes/<box>/`) when the build config lists boxes, or by an
application's `boxes` rule from the build (same code, `wholecell.export.simbox`). Reads **level P pack tiles**
(lossless; `tile_format_v2.md` §1), pre-pass, library, cellgeom_v1. Successor of the yeast manifest (Y formats §6a)
plus the human SDF walls (H formats §6).

| File | Content |
|---|---|
| `box.json` | `schema`, `build_id`, `cell_id`, `cellgeom` {id, sha256}, `library_manifest_sha256`, exporter commit, box {lo_nm, hi_nm, boundary: periodic\|walls\|embedded, shell_nm, tether}, frame, units (nm, ns, µm²/s), `temperature_K`, `viscosity_Pa_s`, `ionic_strength_M` (from config `export.sim_defaults`; replaces yeast 303.15 K ×4 and human 310 K literals), diffusion model {name, frozen input + md5}, counts, seed of the box choice, `secms_bound_fraction` {status, table, table_sha256, secms_source_id, threshold_k, columns, intended_use} |
| `molecules.parquet` | one row per molecule: `mol_id, species (library id), species_idx, region (name), region_id, x_nm, y_nm, z_nm, qw, qx, qy, qz, R_h_nm, D_t_um2_s, D_r_per_us, fixed (bool), tethered (bool), src (packed\|prepass:<module>), polymer_id, polymer_pos`, `secms_bound_fraction`, `secms_unrepresented_share_of_free` (monomer species: the protein's SEC-MS a_k and the share of its modelled free copies measured bound in unmodelled assemblies, `inventory_schema.md` §4c; NaN for assemblies, pre-pass and unmeasured proteins; metadata, the stub D_t / R_h do not use it) |
| `bodies.npz` | collision bodies per species (no pickle): `species` (str), `off` (CSR), `centres` float32 (n,3), `radii` float32 (n), `level` attr (single / r1.6 / r2.5 / r4 as BD v1), `rh, rg, vol, mass` per species; `np.load(allow_pickle=False)` must work |
| `beads.cif` | mmCIF, one pseudo-atom per residue (CA) / nucleotide (P), one chain per molecule, occupancy 1 mobile / 0 fixed, `_cell.length_*` in Å; membrane midplane beads as chain `ZZZZ` (yeast convention) |
| `walls.npz` | cellgeom_v1 crop: per field/membrane `sdf_<name>` float32 grid at `h_nm` (default 2, as human lumen boxes), `origin_nm`, `label_grid` uint8 at 10 nm + `label_names` (str array); signs/offsets exactly as cellgeom_v1 |
| `polymers.npz` | optional: mRNA / chromatin bead polylines (`pts`, `off`, `kind`) |
| `MANIFEST.sha256` | checksums of the files above |

Quaternions are (w,x,y,z) (`tile_format_v2.md` §2); the legacy yeast manifest key `quat_xyzw` is converted on
import. Converters: `legacy_bd_box.py` (human `box_<name>.npz` + pickled `lib_<level>.npz` → simbox, run once,
read-only on the legacy files) and the yeast browser export reader. `src`/`fixed` take the tile flags `prepass`
and `fixed`.

## 4. Frozen pre-registration as a workflow feature

```
preregistration/
  2026-11-20_bd_heldout.md        the plan: hypotheses, metrics, thresholds, arms, analysis code path
  2026-11-20_bd_heldout.md5       md5 of the .md (written by `python -m wholecell.prereg freeze`)
  prereg.yaml                     machine-readable part (below); its md5 is inside the .md5 file as well
```

```yaml
prereg: 2026-11-20_bd_heldout.md
build_id: 20261116-abc1234-s20261115
heldout: {path: preregistration/heldout_species.tsv, md5: ...}          # frozen split (seeded, recorded)
frozen_inputs: {bd_theta: md5..., benchmark: md5...}                   # names from config.frozen_inputs
evaluation_code: {path: scripts/evaluate.py, md5: ...}
tests: [{id: P1, metric: heldout_mean_abs_log10, op: "<=", value: 0.25}, ...]
```

| Rule | Does | Writes |
|---|---|---|
| `prereg_freeze` (manual target) | md5 of .md and prereg.yaml; refuses if `results/` of this application already contain `evaluate` outputs for the test; commits nothing itself (author commits) | `.md5` |
| `prereg_check` | verifies the .md5; every md5 in prereg.yaml against the files; `build_id` = config build; prereg file committed in git with commit time earlier than the first `run` output of the run (git log + output mtimes recorded in `run.json`) | `results/<run>/prereg_check.json` |
| `evaluate` | depends on `prereg_check`; computes only the declared tests; refuses if `eval.lock` exists | `eval.json`, `eval.lock` {time, prereg md5, input md5s, code commit, verdict per test} |
| `posthoc` | anything not pre-declared; may read `eval.json`; outputs labelled `posthoc_*` | never touches the lock |

Changes after a lock: a dated addendum file next to the prereg (as `docs/reference/history/v7_reproduction_criteria.md`), never an
edit. Existing legacy preregs (table §2a) are imported verbatim with their recorded md5s.

## 5. Parameters, determinism and deposits

**Hand-written / live-computed parameters become rule outputs or declared frozen inputs:**

| Legacy item | Problem (H stages §10, rebuild_test §7) | New form |
|---|---|---|
| BD v2 `theta_t0_*.json`, `theta_final.json` written by hand | not reproducible | `fit` rule output with the selection criterion in config; if a human choice remains, a `decision` frozen input {path, md5, origin: "hand-chosen 2026-11-10", why} |
| S2.0 parameters computed at run time from `ledger_long.parquet` | a ledger rebuild silently changes S2.0 | `params` rule writes `s2_params.json` from a pinned ledger (frozen input); legacy values exported once and checked equal |
| `BD_PROD` env from `logs/bd_prod_override.txt` | lumen runs lacked it | config key `run.prod_settings`; no env-based switches |
| epsmap seeds from `hash()` | not reproducible | `wholecell.seeds` named streams; `PYTHONHASHSEED=0` |
| d2 phase waiting on a PID; S2 jobs trimmed by hand | hidden dependency/selection | Snakemake dependencies; job tables as rule outputs with filters in config |

**Non-deterministic GPU steps.** Rules listed in `deposit:` (float32 CuPy BD, OpenMM CUDA, HOOMD) are run once;
their outputs are **deposited** (checksummed, then treated as frozen inputs by downstream rules). The target
`rerun_compare` re-runs them (new run_id) and compares against the deposit with the config tolerances, defaulting
to H rebuild_test §8:

| Artefact | Tolerance |
|---|---|
| BD D/D0 per probe | within run-to-run error (~3-5 %); held-out mean \|log10\| ± 0.02 |
| condensates | c_sat within replica spread |
| S1/S2 contact metrics, application verdicts | within seed spread; verdicts identical |
| chromatin | replica metrics within the 5-replica spread (EM voxel r 0.68 ± 0.02) |

Each run writes `run.json` (seed, code md5 per script, env lock hash, GPU model, driver, dt_final, wall time),
as the condensate driver already does (DECISION 2026-11-15 04:15).

## 6. GPU job submission (abstracted)

Rules declare resources only: `resources: gpu=1, gpu_mem_gb=8, runtime=600, cpus=4` and `conda: envs/<engine>.yaml`.
Execution is chosen by a Snakemake profile in `profiles/`: `local` (one shared GPU: `--resources gpu=1` acts as a
semaphore, one CUDA process at a time as in the residue chain), `slurm` (executor plugin; `slurm_extra` gres,
partition and account **left as placeholders until Pedro specifies the GPU nodes**). No node names, queues or
absolute paths in rules. Long chains run detached (`snakemake ... &` with a post-script writing a status file);
no agent monitors them (MEMORY "token economy").

## 7. Legacy results as the archival record; comparing a re-run

`legacy.yaml` lists, per legacy result used in the paper: path (under `{data_root}`, read-only), md5, producing
script + commit/date, build it came from (v4/v7 pack), known defects (H rebuild_test §7, e.g. Fig5 lumen files
superseded by `*_v7sh`). The `compare_legacy` rule reads the new run and the legacy files and writes
`compare_legacy.tsv` (metric, legacy, new, Δ, tolerance, pass) + `verdict.json` with the three-level verdict of
the v7 criteria (REPRODUCED / WITH DEVIATIONS / NOT REPRODUCED). Comparison criteria are written and committed
**before** the re-run (same discipline as §4). Legacy files are never copied into `results/`.

## Open questions

1. GPU nodes, partitions and accounts for the `slurm` profile (Pedro).
2. Should boxes be part of every build export (config `export.boxes`) or only cut by applications on demand?
3. Engines shared across applications (BD in `bd_diffusion` and `nfat_calcium`): one `wholecell/sim/bd` package,
   or keep the human CuPy engine and yeast `simgpu` side by side until a benchmark decides?
4. Prereg timestamp check relies on git commit times of a local repo; is a signed tag or an external timestamp
   (Zenodo/OSF) wanted for the paper?
5. Which legacy applications are re-run on the port build at all (vs kept only as archival record)?
