# Step-2 QUICK subsets (development checks), 2026-11-18

Pedro (via coordinator, 2026-11-18): the frozen Step-2 subsets (`step2_subsets.md`: Jurkat 56, yeast 58 tiles,
deliberately the busiest / membrane-richest tiles) stay for final acceptance only. For development iterations a
QUICK subset of ~20 representative, NOT busy tiles per cell. Same harness, same neighbour rule (§3), metrics (§4) and
comparisons (§5); a quick run compares only with a quick run or with a full build read on the quick tiles. No
acceptance decision is taken on a quick subset.

| Cell | File | Tiles | md5 (file) | md5 (`ref_tile`, int64) | Requested instances (frame build) |
|---|---|---|---|---|---|
| Jurkat | `config/subsets/jurkat_quick.tsv` | 19 | `3afafc1552604a463f91247b43946d51` | `0b4e174284050610db739649de115a19` | 2.27 M (frozen subset 9.83 M) |
| Yeast | `config/subsets/yeast_quick.tsv` | 19 | `2385e86f7d66b5cebd1b8e47db0f3a9d` | `51a1927c4cac25b275659e281277ceb6` | 0.11 M (frozen subset 0.50 M) |

Rule (`scripts/step2_subset_test.py select_quick`, frames as in `step2_subsets.md` §1): candidates = tiles with
requested instances ≤ the 75th percentile of the tiles with requests; region measure columns in order of their cell
total: the candidate closest to the median of the candidates' non-zero measure (17 columns at most); one cytosol tile
≥ 95 % (median requested); then, for any region with requests still absent, its median tile (Jurkat: pcm; yeast:
LD_peripheral). Every requested region of the frame build has requests on the quick subset. One tile overlaps the
frozen subset per cell.

Usage: `step2_subset_test.py run ... --subset quick [--reuse-auto <previous run or build>]`. `--reuse-auto` links every
upstream part (geometry, tiles, library, inventory, engine, _smk, pre-placed, pack organisation) whose code (git diff on
the part's paths) and config sections are unchanged and whose prerequisites are reused; the rest is rebuilt; the
pre-placed stage is also reused after a rebuilt inventory when the pack table, library manifest and geometry are
content-identical. Expected cost of a quick run when only the pack changes: Jurkat ~13 min library load + ~15 min pack;
yeast ~25-35 min pack (last-round squeeze phases bound by their slowest tile).
