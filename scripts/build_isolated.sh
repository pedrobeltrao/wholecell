#!/bin/bash
# Launch an isolated, frozen build:  scripts/build_isolated.sh <cell> [tag] [cores]
#   1. refuses a dirty working tree (tracked files), resolves the build id (wholecell.plan: <date>-<commit7>-s<seed>[-tag])
#   2. creates <results_root>/<cell_id>/<build_id>/ and a frozen `git worktree` of HEAD in <build>/_code
#   3. runs snakemake from the worktree, detached (setsid nohup), --cores <cores> (default: site max_workers)
#   4. writes <build>/launch.log (command, commit, start time) ; snakemake output goes to <build>/snakemake.log
# Post-run steps (validation, compare with v7 for Jurkat) are rules of the DAG: nothing needs polling.
# Env: WHOLECELL_SITE (site file), SMK_EXTRA (extra snakemake arguments, e.g. "--until pack_index").
set -euo pipefail
CELL=${1:?usage: build_isolated.sh <cell> [tag] [cores]}
TAG=${2:-}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SITE=${WHOLECELL_SITE:-$REPO/config/site.yaml}
site() { grep -E "^$1:" "$SITE" | head -1 | sed -E "s/^$1:[[:space:]]*//; s/[[:space:]]+#.*$//"; }
PY=$(site python); SMK=$(site snakemake); MAXW=$(site max_workers)
CORES=${3:-${MAXW:-8}}
cd "$REPO"
if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
  echo "refusing: working tree has uncommitted changes (commit first; builds run from a frozen worktree of HEAD)" >&2
  exit 2
fi
COMMIT=$(git rev-parse HEAD)
PLAN=$(PYTHONHASHSEED=0 PYTHONPATH="$REPO" "$PY" -m wholecell.plan --cell "$CELL" ${TAG:+--build-tag "$TAG"} | grep '^WHOLECELL_PLAN ' | sed 's/^WHOLECELL_PLAN //')
read -r BUILD_ID BDIR < <(echo "$PLAN" | "$PY" -c 'import json,sys,os; p=json.load(sys.stdin); print(p["build_id"], os.path.join(p["results_root"], p["cell_id"], p["build_id"]))')
if [ -e "$BDIR" ]; then
  echo "refusing: $BDIR exists (no in-place rebuilds; use another tag)" >&2
  exit 3
fi
mkdir -p "$BDIR"
git worktree add --detach "$BDIR/_code" "$COMMIT" > /dev/null
{
  echo "cell: $CELL"; echo "build_id: $BUILD_ID"; echo "build_dir: $BDIR"; echo "commit: $COMMIT"
  echo "code: $BDIR/_code (git worktree, detached)"; echo "cores: $CORES"; echo "started: $(date -Iseconds)"
  echo "host: $(hostname)"
  echo "command: $SMK --snakefile $BDIR/_code/Snakefile --directory $BDIR --cores $CORES --config cell=$CELL build_id=$BUILD_ID build_dir=$BDIR --rerun-incomplete ${SMK_EXTRA:-}"
} > "$BDIR/launch.log"
cd "$BDIR"
PYTHONHASHSEED=0 setsid nohup "$SMK" --snakefile "$BDIR/_code/Snakefile" --directory "$BDIR" --cores "$CORES" \
  --config cell="$CELL" build_id="$BUILD_ID" build_dir="$BDIR" --rerun-incomplete ${SMK_EXTRA:-} \
  > "$BDIR/snakemake.log" 2>&1 < /dev/null &
echo "pid: $!" >> "$BDIR/launch.log"
echo "$BDIR"
