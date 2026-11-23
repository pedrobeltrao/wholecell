# B4 pre-placed builders: addendum 2026-11-18 (polysome walk fix; J3b / J4 a documented model limitation for Jurkat)

Addendum to `docs/unification_criteria/B4_preplaced.md` (frozen, md5 e715c2f90e17653369b7d5ac2c4d3d96) and
`B4_preplaced_addendum_2026-11-17.md` (md5 ad6dbb8bf4395efd251b6f8ff5ec6574); neither is edited. **No criterion
changes**: the frozen J3a / J3b / J4 bounds and the 2026-11-17 J5b rule stand. Written after a timeboxed attempt (main
session decision 2026-11-18, Pedro informed: keep one common walk; the problem is chain growth, not the minimum
distance; make the walk continue chains at contact spacing with the measured linker, generically, yeast unchanged; if
J3b / J4 still fail, accept it as a documented model limitation; no EM-seeded spacing as a Jurkat-only case).

## 1. What changed in the walk (generic option, yeast default unchanged)

`core.polysomes` `params.walk` = {candidates, rounds, pick, dir_cone, orient_cone, priors_first}
(`wholecell/preplaced/builders/polysomes.py`, `WALK_SEARCH`).

- **candidates 0 (default, absent) = the legacy walk** (one pose per try, entry direction exactly against the link,
  free nt redrawn after 40 failed orientations, chain ends after 1,200 failures). Yeast does not set `walk`: its
  pre-pass is bitwise identical to the PATCH-4 reference (§4).
- **candidates K > 0**: each grown ribosome after the first (not a collided disome, not an EM prior) is chosen among
  rounds of K candidate poses on the shell the link can reach: entry = previous exit + L·D, L the end-to-end draw
  (Gaussian ssRNA chain) of the round's free nt (mean = `linker_nt_mean`, unchanged 150), D within `dir_cone` of the
  exit direction (default the legacy cone), entry direction within `orient_cone` of −D, random roll. Of the candidates
  that fit (no overlap, allowed labels, ER tether) in the first round with any fit, `pick` places the one nearest the
  previous centre (`nearest`), nearest the centroid of the chain so far (`compact`) or the first (`first`). Each round
  redraws the free nt (as the legacy walk does after 40 failures); **the chain ends only when rounds + 1 rounds give no
  fit**.
- `priors_first` (point priors only): all chains' EM priors are placed before any chain grows.

**Jurkat** (`config/cells/jurkat.yaml`): `walk: {candidates: 40, rounds: 29, pick: compact, orient_cone: 0.7854}`
(the legacy budget 40 × 30 = 1,200 tries per ribosome), `linker_nt_mean` unchanged at 150.

## 2. Calibration subset (addendum 2026-11-17 §1 subset and score; model = pre-placed + random cytosolic complement, thinned as J3b / J4)

`scripts/b4_polysome_calibration.py run --walk <json>` / `score`; results
`jurkat/_dev_unify/b4/walk_fix/calib/` (`score_*.json`; one full Jurkat walk per row, npc + microtubules +
polysomes, geometry em_cell1_v8r). Subset EM: NN 27.20 nm, linked 0.691.

| walk | linker nt | J3a | pre-placed | EM priors dropped | consecutive centres, median nm | subset NN nm | subset linked | score |
|---|---|---|---|---|---|---|---|---|
| legacy (npc_fix stage = before) | 150 | 0.684 | 1,615,692 | 130,940 | 31.58 | 35.38 | 0.437 | 13.26 |
| K40 R30, entry ±45°, nearest | 150 | 0.689 | 1,626,086 | 129,480 | 28.75 | 33.86 | 0.482 | 10.83 |
| **K40 R30, entry ±45°, compact (chosen)** | 150 | **0.688** | 1,623,965 | 126,561 | 31.08 | **33.40** | **0.499** | **10.03** |
| … + priors_first | 150 | 0.685 | 1,618,512 | 81,495 | 31.18 | 33.47 | 0.498 | 10.11 |
| … + priors_first + any direction | 150 | 0.697 | 1,644,241 | 81,495 | 31.65 | 33.74 | 0.486 | 10.63 |
| nearest, priors_first | 150 | 0.686 | 1,621,338 | 81,495 | 28.95 | 33.77 | 0.484 | 10.70 |
| compact | 60 | 0.663 (fails J3a) | 1,571,665 | 122,973 | 30.31 | 33.15 | 0.501 | 9.74 |
| compact | 20 | 0.538 (fails J3a) | 1,311,052 | 108,797 | 29.47 | 34.39 | 0.468 | 11.65 |

Not scored (development, `calib/superseded/`, `calib/killed/`): a first version that kept the free nt fixed per link
and drew fully random entry orientations ended chains early (1.13-1.33 M pre-placed); aligned entry orientation
(±0°) with the candidate search was > 3× slower and was stopped. With the walk fix, shorter linkers still end chains
(J3a fails at 60 nt), so the measured 150 nt stays. Choice = lowest subset score among J3a-eligible rows (the
2026-11-17 rule): `compact`.

## 3. Full Jurkat pre-placed stage (frozen J3 / J4; J5 with the 2026-11-17 J5b rule)

`jurkat/_dev_unify/b4/walk_fix/` (`run_stage.sh`, `stage/preplaced/serial`, `measure_stage.json`; serial 1,042 s,
47 GB, measure 1,097 s; polysome arrays identical to `calib/run_k40r30o45c`); packed free 80S = the addendum's random
cytosolic complement (pool − pre-placed = 481,863); lamina = the rebuilt lattice of `npc_fix` (NPCs unchanged).
Before = `npc_fix/measure_stage.json` (legacy walk, same geometry; = the 2026-11-17 stage values).

| | frozen bound | before (legacy walk) | walk fix | pass |
|---|---|---|---|---|
| J3a polysome fraction | 0.68-0.72 | 0.684 | **0.688** | yes |
| J3b linked share, pre-placed + complement, thinned | 0.697 ± 0.05 | 0.441 | **0.506** | **no** |
| J3b pre-placed only (reported) | | 0.568 | 0.659 | |
| J4 NN median, pre-placed + complement, thinned | 27.24-29.24 nm | 35.22 | **33.25** | **no** |
| J4 pre-placed only, thinned (reported) | | 32.13 | 30.89 | |
| J4 all model ribosomes, unthinned (reported) | | 29.09 | 28.64 | |
| J5a pre-pass pairs per 10⁶ (v7 layers 842,927 + 10 %) | ≤ 927,220 | 823,355 | 823,284 | yes |
| J5b polysomes-polysomes per 10⁶ polysome instances (limit 721,838) | | 29,234 (2,029 pairs) | 34,273 (2,402 pairs) | yes |
| J5b npc-polysomes per 10⁶ (limit 580) | | 43 (3) | 14 (1) | yes |
| J5b cytoskeleton_rest-polysomes per 10⁶ (limit 178,201) | | 274 (19) | 328 (23) | yes |
| J5b microtubules-polysomes, lamina-polysomes | | 0, 0 | 0, 0 | yes |
| J5b classes without polysomes (MT-MT, chromatin-NPC, cytoskeleton_rest-MT) | | unchanged | unchanged | yes |

J5b (2026-11-17 rule): **pass in every class** (no failing class; the NPC classes pass after the npc-overlap fix with
the rebuilt lamina). **J3b and J4 still fail** on the full stage: the walk fix moves them by +0.065 / −1.97 nm, about
a third of the gap to the window in both (J3b 0.506 vs ≥ 0.647; J4 33.25 vs ≤ 29.24 nm).

## 4. Accepted as a documented model limitation (Jurkat J3b / J4)

**Decision (main session, 2026-11-18):** the remaining J3b / J4 failure of the Jurkat pre-placed stage is accepted as a
model limitation; the bounds are not changed; v7's EM-seeded spacing is not used. `scripts/b4_preplaced_check.py`
keeps reporting J3b / J4 as fail; this addendum is the documented reason.

**What limits them** (measured, `walk_fix/calib`):

1. **The thinning model, not the chain geometry, sets the floor.** J3b / J4 thin the model (2,105,828 RIBO80S of the
   inventory) at random to the EM count (730,960; p = 0.347), i.e. they assume EM detects a random third of the
   ribosomes. The EM point set itself, thinned by the same ratio, gives NN median 39.4 nm and linked share 0.335 —
   far outside the windows. To pass, the model must have about 2.9× EM's local contact density: every ribosome needs
   on average ≈ 2.8 neighbours within 29 nm (EM: 0.96 at full count).
2. **Chain topology bounds the contacts.** A chain gives each ribosome at most two sequence neighbours; the compact
   walk reaches 0.87 neighbours within 29.24 nm per pre-placed ribosome (legacy 0.67; within 34 nm 2.53 vs 1.95),
   i.e. EM-like spacing at full count, but not the 2.9× clustering the thinned statistic asks for. Consecutive
   centres are already at the contact floor (q10 ≈ 26.3 nm; collided disome 27.0 nm; nearest pick median 28.7 nm).
3. **The packed free 80S** (23 % of the pool, here a uniform random complement) are mostly isolated after thinning and
   pull the median up (pre-placed only: 30.9 nm / 0.659 vs with complement 33.3 nm / 0.506).
4. EM priors that do not fit (126,561 of 682,785, 19 %; priors_first lowers this to 81,495 but does not improve the
   thinned statistics).

Closing the gap would need either a different reading of J3b / J4 (EM detection not a random subsample, or a smaller
ribosome pool) or ribosome clusters denser than the EM groups (v7) — both outside this walk fix and outside the frozen
criteria.

## Checks of this addendum

- Code: commit 9828d98 (runs: e89d2ba + its diff, `code.diff` / `stage_code.diff` in the run folders).
- Yeast (Y1): serial pre-pass of `config/cells/yeast.yaml` with this code (vcell env, 2 threads, 2,635 s;
  `yeast/_dev_unify/b4/run_vcell_walk_fix/`): builders 1-18 identical to s0, 19-28 (polysomes, mrna included)
  identical to `ref_yeastenv_patched` (`compare_ref.md`).
- Jurkat: calibration walks `jurkat/_dev_unify/b4/walk_fix/calib/` (1 thread each, 11-15 min); full stage
  `walk_fix/stage` + `measure_stage.json`.
- `tests/test_preplaced_b4.py` (new: `walk_option`, the walk default is legacy, Jurkat sets it),
  `tests/test_preplaced_serial.py` pass (as scripts).
