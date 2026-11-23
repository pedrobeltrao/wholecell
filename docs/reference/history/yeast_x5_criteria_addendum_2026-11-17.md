# Yeast X5 reproduction: criteria addendum 2026-11-17

Date: 2026-11-17. **Decided by Pedro 2026-11-17.** Addendum to the frozen criteria
`docs/reference/history/yeast_x5_reproduction_criteria.md` (not edited; frozen files and addenda are separate files).
Scored by `plugins/yeast/repro/seedspread_x5.py`, which reuses the extraction code and tolerance constants of
`plugins/yeast/repro/compare_x5.py` (version 8f1d0e1, which scored the s0 port).

**Frozen before any seed-1/2 output was read.** This file and its scoring code were committed while the port builds
with master seeds 1 and 2 (code aabe7b5, `.../wholecell_builds/yeast/20261117-aabe7b5-s{1,2}-x5repro`) were still
running and before any of their files had been opened. The only port build seen when writing it is s0
(`20261116-aabe7b5-s0-x5repro`, verdict NOT REPRODUCED; must rows failing L5, S3, P2, P3, V1, V2, V3; report rows
failing P4, P7, V4, X2).

## 1. Why

The s0 failures are of two kinds. (i) Non-content differences that the frozen criteria did not exclude explicitly:
the chromatin builder's wall time in the s3 log ("3s" vs "4s"); the `sig` strings of `pack_trees_15.npz` (C and R
identical); the order of equal-strength partners in the partner CSR (strengths differ by ≤ 1.4e-15; 13 of 69,656
indices swap among ties). (ii) Stochastic outcomes of the chaotic packing and the GPU steps, where an exact or
tight comparison with one X5 seed is not meaningful (e.g. P3: 18 (type, region) pairs requested 1 vs 0). For (ii) the
port is judged against its own seed spread (master seeds 0, 1, 2), and that spread must be consistent with X5 or with
the W4 seeds 0–2 (same yeast model, one round earlier).

## 2. Builds

- Port seeds: s0 = `20261116-aabe7b5-s0-x5repro`, s1/s2 = `20261117-aabe7b5-s{1,2}-x5repro` (same code, master seed 1/2).
- X5: `{yeast_root}/runs/X5_seed0/work` (READ-ONLY).
- W4: `{yeast_root}/runs/W4_seed{0,1,2}/work` (READ-ONLY).

## 3. Rules for S3, L5, P2 (s0 vs X5, as in the frozen criteria)

- **S3**: timing tokens — a number immediately followed by `s` (regex `(?<![\w.])\d+(?:\.\d+)?s\b`, e.g. `3s`,
  `12.4s`) — are removed before the numbers of a builder line are compared. All other numbers within 1 %, as before.
- **L5**: `lib_ids.npy` identical (unchanged); `libpack_meta.npz` and `pack_trees_15.npz` are compared on every key
  except the signature string `sig`: non-float keys identical, centres/radii (float keys) |Δ| ≤ 0.05 nm. With `sig`
  ignored, identical remaining keys count as level E.
- **P2**: tie-order permutations of equal-strength partners are equivalent. `network.npz`: as before (non-float keys
  identical, floats rel. Δ ≤ 1e-5). `partners.npz`: `ids` and `off` identical; each CSR row is compared as a set of
  (partner index, kind) — rows sorted by (index, kind, strength) — and strengths of matched elements within rel. 1e-5.

## 4. Seed-spread rule (P3, V1, V2, V3 must; P4, P7, V4, X2 report)

Each row is reduced to scalars (section 5). For one scalar with port values v0, v1, v2 (seeds 0, 1, 2), X5 value x
and W4 values w0..w2:

- **Margin** (fixed now): `margin = max(0.5 × (max(v) − min(v)), tol)`, where `tol` is the existing tolerance of the
  metric in the frozen criteria / compare script, evaluated at scale `n = max(|min v|, |max v|, |x|)` (so a count
  tolerance max(3√n, 0.5 % n) absorbs pairs like 1 vs 0: n = 1 gives a margin of 3). Interval
  `I = [min(v) − margin, max(v) + margin]`.
- **Reproduced** (scalar passes) iff (a) every port seed value lies in I (true by construction when all three seeds
  have the metric; a missing port value makes the scalar NA), AND (b) I contains x OR I overlaps
  [min W4, max W4] (over the W4 seeds that have the metric).
- Where a metric has no W4 counterpart (file or key missing in all three W4 seeds), (b) is tested against X5 only; the
  output records this per scalar (`ref = "X5 only (no W4 counterpart)"`) and counts it per row. If neither X5 nor
  W4 has the metric, the scalar is NA. An NA scalar never passes.
- Output per scalar: port min/max, margin, tolerance, interval, X5 value, W4 range, which reference matched
  (X5, W4, X5+W4, none).
- Level of a row scored this way: `SS` (passes the seed-spread rule), `fail`, or `NA` (nothing evaluable or an
  extraction error). NA counts as not passed (must and report rows alike).
- Booleans/flags and lists: flags are never scored; where the flag summarises a number (e.g.
  `acceptance.clash_rate_lt_1e-5` → the clash rate), that number is scored as its own scalar. Lists (per-bin vectors,
  `clash.top_type_pairs[*]` entries) and strings are report-only (not gating), as in compare_x5.

## 5. Metric definitions per row (scalars, tolerance, row-level pass rule)

| Row | Class | Scalars | tol (frozen criteria / compare_x5) | Row passes if |
|---|---|---|---|---|
| P3 | must | requested copies per (library type name, region), summed over tiles from `yc3d/tile_counts.parquet` (types keyed by name via `lib_ids.npy`; union of port and X5 pairs, a missing pair = 0; a W4 pair whose type is not in that W4 library has no W4 counterpart); every numeric leaf of `yc3d/thinning.json` | pairs max(3√n, 0.5 % n); thinning ± 0.005 | all scalars pass |
| V1 | must | every non-volatile numeric (non-bool) leaf of `results/wholecell3d_validation.json` (key union of X5 and port) | compare_x5 V1 leaf classes: integers max(3, 1 % n); clash/outside rates max(25 % n, 1e-5); median_angle 0.01°; density 2 mg/ml; values in [0, 1] ± 0.005; other max(1 % n, 0.01) | all scalars pass |
| V2 | must | per file (X5 `validation/*.json` + `validation/latest/*.json`, minus MANIFEST and diffusion_benchmark): every non-volatile numeric leaf; a file missing in a port seed is NA | max(1 % n, 0.01) | a file passes if all its scalars pass; the row passes if all files pass |
| V3 | must | contact pair count (sum of `C*` in `contacts/contacts.npz`); kinase encounters: row count and means of p_near_2, p_near_10, c_local_uM; φ of `sweeps/{cyto120_mf,gem200_mf,nuc100_mf}.json`; diffusion benchmark: the four **test** RMSE values (crowding abs/ratio = model `+ crowding(R_h) from BD`, binding abs/ratio = model `+ crowding + binding`; `rmse_log10_abs` / `rmse_log10_ratio`) | contacts 1 %; encounters 5 %; φ 2 × X5 box SD (cyto 0.026, gem 0.020), nuc100 ± 0.01; RMSE 10 % (the D/D0 10 % tolerance) | all scalars pass |
| P4 | report | placed instances per tile (`yc3d/pack_summary.parquet` `placed`, 5,912 tiles, keyed by `ti`); total placed | per tile max(3√n, 1 % n) (the P5 count tolerance); total 0.5 % (P6) | total passes and ≥ 99 % of tiles pass |
| P7 | report | number of `yc3d/tiles_raw/*.npz` files; total placed (sum of `job_placed`); placed per library type name (`job_lib`, `job_placed`) | file count exact (0); total 0.5 %; per type max(3√n, 1 % n) | file count and total pass and ≥ 99 % of types pass |
| V4 | report | every number (by position) of the `log phase` / `stress` lines of `condensates/cond_test.log`; each count of `microscopy/spot_counts.json` | lines max(1 % n, 0.05); spots max(1, 1 % n) | all scalars pass |
| X2 | report | number of viewer tile files (`yc3d/web/data/tiles/*`); number of viewer library files (`yc3d/web/data/library/*`) | exact (0) | both pass |

The ≥ 99 % share rule for P4 tiles and P7 types: thousands of per-tile / per-type tests at a ~3σ count margin expect a
few chance exceedances, which the frozen criteria already accept in P5 (≥ 99 % of small pairs within ± 3).

V3 diffusion benchmark: the test RMSE (log10) reference values below are used **in addition to** whatever the script
finds in the W4 files (W4 range = W4 file values ∪ table W4 rows; the X5 value is read from the X5 file and the table
X5 row is used only if the file lacks it). The other numbers of `diffusion_benchmark.json` (train RMSE, Spearman) are
not scored.

| Build | crowding abs | crowding ratio | binding abs | binding ratio |
|---|---|---|---|---|
| X5 | 0.207 | 0.091 | 0.241 | 0.163 |
| W4 s0 | 0.207 | 0.091 | 0.241 | 0.163 |
| W4 s1 | 0.212 | 0.095 | 0.248 | 0.170 |
| W4 s2 | 0.212 | 0.090 | 0.254 | 0.187 |
| port s0 | 0.218 | 0.093 | 0.263 | 0.196 |

## 6. All other rows

Every other row (G1–G5, I1–I4, A1, L1–L4, S1, S2, S4, P1, P5, P6, V2a, V5, X1) keeps its s0 result from
`<s0>/validate/compare_x5_final/` (E stays E, T stays T).

## 7. Verdict

- **REPRODUCED**: every must row passes under the frozen criteria + this addendum (level E, T or SS) and every
  report row passes.
- **REPRODUCED WITH DEVIATIONS**: every must row passes and at least one report row fails (or is NA).
- **NOT REPRODUCED**: otherwise (any must row fails or is NA).

## 8. Running

```
conda_envs/yeast_x5/bin/python plugins/yeast/repro/seedspread_x5.py --s0 <build s0> --s1 <build s1> --s2 <build s2> \
    --x5 {yeast_root}/runs/X5_seed0 --out <dir> [--w4 <W4_seed0/work> <W4_seed1/work> <W4_seed2/work>] [--nproc 4]
```

Writes `<out>/verdict.json` (verdict, counts, per-row level/pass and seed-spread summary with margin, port min/max,
X5 value, W4 range, matched reference), `<out>/criteria_results.tsv`, `<out>/details.json` (all non-passing scalars
per row) and `<out>/scalars.tsv.gz` (every scalar). Self-test before the seeds: `--s0 = --s1 = --s2 = s0` (degenerate
spread, margin = tolerance) must run end to end and give S3, L5 and P2 a pass on s0.
