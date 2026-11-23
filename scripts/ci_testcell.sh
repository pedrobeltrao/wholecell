#!/bin/bash
# Testcell CI: build the testcell of the current commit twice (8 and 4 workers) and check the two builds are identical.
#   scripts/ci_testcell.sh            run in the foreground (about 5-6 min; 12 cores)
#   scripts/ci_testcell.sh --detach   start in the background and return (the post-commit hook does this)
# Builds run one after the other (8 cores at a time; CI_SEQUENTIAL=0 runs them together, 12 cores). Set WHOLECELL_CI=0 to skip (e.g. while committing many small doc changes).
#
# Results: <results_root>/testcell/_ci/<build id of the 8-worker build>/
#   ci.json         status, wall time of each build, determinism result, comparison with the previous passing CI
#   determinism.json, vs_previous.json   file-level details (scripts/ci_compare_builds.py)
#   STATUS          one line: PASS | FAIL <reason>
# PASS needs: both builds complete, the 8-worker build under CI_MAX_SECONDS, and 8- vs 4-worker builds identical.
# The comparison with the previous passing CI build is reported, not judged: a refactor commit expects "identical",
# a commit that changes the model expects (and must explain) differences.
# Only one CI runs at a time (flock); a commit that already has a CI result is skipped.
set -uo pipefail
# The whole body is one function, read completely before it runs: editing this file during a CI run is safe.
main() {
  [ "${WHOLECELL_CI:-1}" = "0" ] && exit 0
  # a git hook exports GIT_DIR, GIT_INDEX_FILE, ... for the repository it runs in; they would break `git worktree add`
  unset $(env | sed -n 's/^\(GIT_[A-Z_]*\)=.*/\1/p')
  REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
  if [ "${1:-}" = "--detach" ]; then
    setsid nohup "$REPO/scripts/ci_testcell.sh" > /dev/null 2>&1 < /dev/null &
    exit 0
  fi
  SITE=${WHOLECELL_SITE:-$REPO/config/site.yaml}
  site() { grep -E "^$1:" "$SITE" | head -1 | sed -E "s/^$1:[[:space:]]*//; s/[[:space:]]+#.*$//"; }
  PY=$(site python); ROOT=$(site results_root)
  CI_MAX_SECONDS=${CI_MAX_SECONDS:-600}           # plan A0: testcell full build < 10 min on 8 cores
  CIROOT="$ROOT/testcell/_ci"; mkdir -p "$CIROOT"
  C7=$(git -C "$REPO" rev-parse --short=7 HEAD)

  exec 9> "$CIROOT/.lock"
  flock 9                                          # wait for any running CI

  # build ids: <date>-<commit7>-s<seed>-ci[N]-w8 / -w4 (scripts/build_isolated.sh refuses a dirty tree or an existing id).
  # A commit with a CI result (STATUS) is skipped; after a failed launch the next attempt uses ci2, ci3, ...
  if ls "$ROOT"/testcell/_ci/*-"$C7"-*/STATUS > /dev/null 2>&1 && ! grep -qs "^FAIL launch" "$ROOT"/testcell/_ci/*-"$C7"-*/STATUS; then
    exit 0
  fi
  N=""; while ls -d "$ROOT"/testcell/*-"$C7"-*-ci$N-w8 > /dev/null 2>&1; do N=$(( ${N:-1} + 1 )); done
  # Launch from a private clean worktree of the commit, so edits in the live checkout cannot interfere
  # (build_isolated.sh refuses a dirty tree; on 2026-11-16 a cherry-pick between the two launches did exactly that).
  SRC="$CIROOT/.src-$C7"
  git -C "$REPO" worktree remove --force "$SRC" > /dev/null 2>&1 || true
  git -C "$REPO" worktree add --detach "$SRC" "$C7" > /dev/null 2>&1
  trap 'git -C "$REPO" worktree remove --force "$SRC" > /dev/null 2>&1' EXIT
  launch() {  # tag cores -> build dir
    "$SRC/scripts/build_isolated.sh" testcell "$1" "$2" 2> "$CIROOT/.launch_$1.err" | tail -1
  }
  wait_build() {  # build dir -> seconds from start to end of snakemake
    local pid; pid=$(sed -n 's/^pid: //p' "$1/launch.log")
    while [ -e "/proc/$pid" ]; do sleep 10; done
    echo $(( $(date -r "$1/snakemake.log" +%s) - T0 ))
  }
  # CI_SEQUENTIAL=1: the 4-worker build starts after the 8-worker one (8 cores at a time instead of 12, for busy nodes)
  T0=$(date +%s)
  B8=$(launch ci$N-w8 8)
  if [ "${CI_SEQUENTIAL:-1}" = "1" ] && [ -d "$B8" ]; then S8=$(wait_build "$B8"); T0=$(date +%s); fi
  B4=$(launch ci$N-w4 4)
  OUT="$CIROOT/$(basename "${B8:-none-$C7}")"; mkdir -p "$OUT"
  fail() { echo "FAIL $*" > "$OUT/STATUS"; cat "$OUT/STATUS"; exit 1; }
  [ -d "$B8" ] && [ -d "$B4" ] || fail "launch: $(cat "$CIROOT"/.launch_ci$N-w*.err 2>/dev/null | tr '\n' ' ')"
  [ "${CI_SEQUENTIAL:-1}" = "1" ] || S8=$(wait_build "$B8")
  S4=$(wait_build "$B4")
  grep -q "(100%) done" "$B8/snakemake.log" || fail "8-worker build incomplete: $B8/snakemake.log"
  grep -q "(100%) done" "$B4/snakemake.log" || fail "4-worker build incomplete: $B4/snakemake.log"

  "$PY" "$REPO/scripts/ci_compare_builds.py" "$B8" "$B4" --out "$OUT/determinism.json" > "$OUT/determinism.txt"
  DET=$?
  PREV=$(cat "$CIROOT/LAST_PASS" 2>/dev/null || true)
  VSP="none"
  if [ -n "$PREV" ] && [ -d "$PREV" ]; then
    "$PY" "$REPO/scripts/ci_compare_builds.py" "$PREV" "$B8" --out "$OUT/vs_previous.json" > "$OUT/vs_previous.txt"
    [ $? -eq 0 ] && VSP="identical" || VSP="different"
  fi
  STATUS=PASS; REASON=""
  [ $DET -ne 0 ] && { STATUS=FAIL; REASON="8- vs 4-worker builds differ (determinism.txt)"; }
  [ "$S8" -gt "$CI_MAX_SECONDS" ] && { STATUS=FAIL; REASON="$REASON 8-worker build took ${S8}s > ${CI_MAX_SECONDS}s"; }
  cat > "$OUT/ci.json" <<EOF
{"commit": "$C7", "status": "$STATUS", "reason": "$REASON", "build_w8": "$B8", "build_w4": "$B4",
 "seconds_w8": $S8, "seconds_w4": $S4, "max_seconds": $CI_MAX_SECONDS, "determinism_identical": $([ $DET -eq 0 ] && echo true || echo false),
 "previous_pass": "$PREV", "vs_previous": "$VSP", "finished": "$(date -Is)"}
EOF
  echo "$STATUS $REASON" > "$OUT/STATUS"
  # LAST_PASS = the last passing build of a commit on main (the "previous" of the next CI); other branches never move it
  ONMAIN=0; git -C "$REPO" merge-base --is-ancestor "$C7" main 2> /dev/null && ONMAIN=1
  [ "$STATUS" = PASS ] && [ $ONMAIN = 1 ] && echo "$B8" > "$CIROOT/LAST_PASS"
  cat "$OUT/STATUS"
  [ "$STATUS" = PASS ]
}
main "$@"
