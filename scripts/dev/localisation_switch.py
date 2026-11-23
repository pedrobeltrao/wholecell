"""Switch a cell's localisation to preset `common` (step B5; docs/unification_criteria/B5_localisation.md §3).

    python scripts/localisation_switch.py [--repo DIR] [--check] jurkat yeast

APPLIED on main 2026-11-18 (Step-2 combined configuration; the reference configs jurkat_v7_port / yeast_x5_core pin the
legacy presets); kept for the record and tests. Originally: main keeps the legacy presets (legacy_v7 / legacy_x5) until the B5 queue job: that job checks out the merged main commit
in a launch worktree, runs this script there, commits the result ("switch commit") and builds it. After acceptance the
same commit is merged into main. Edits, per cell config (config/cells/<cell>.yaml), by exact line match (an unexpected
layout is an error, nothing is guessed):
  - inventory.params.localisation.preset  legacy_v7 | legacy_x5  ->  common
  - inventory.params.localisation.label   -> <cell>.localisation_common (provenance name in species_loc / field_provenance)
  - Jurkat inventory.localisation_sources  term_map strings -> the TSVs the adapters read (DECISION 2026-11-17: cleaned
    up in B5; they named yaml keys that no longer exist), the Itzhak row's role (N/C tie-breaker of the N/C model)
--check: report what would change, write nothing (exit 1 if a cell is not in its legacy state).
"""
import argparse
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
LEGACY = {"jurkat": "legacy_v7", "yeast": "legacy_x5"}
JURKAT_SOURCES = [
    ('term_map: "plugins/jurkat/data/inventory_v0.yaml:hpa_map"}',
     'term_map: "plugins/jurkat/data/localisation/hpa_term_map.tsv"}'),
    ('term_map: "plugins/jurkat/data/inventory_v0.yaml:uniprot_comp_rules"}',
     'term_map: "plugins/jurkat/data/localisation/uniprot_cc_rules.tsv"}'),
    ('term_map: "re-split of nucleus+cytosol+cytoskeleton share"}',
     'term_map: "N/C tie-breaker of the N/C model (HPA category NC; preset common)"}'),
    ("  localisation_sources:                 # -> canonical localisation_sources.tsv (term_map strings kept as in v7 for identical tables;",
     "  localisation_sources:                 # -> canonical localisation_sources.tsv (term_map: the curated TSV each adapter reads;"),
    ("                                        # the maps now live in plugins/jurkat/data/localisation/*.tsv)",
     "                                        # cleaned up in B5, DECISION 2026-11-17)"),
]


def _loc_block(lines):
    """line range of inventory.params.localisation (4-space key `localisation:` under `  params:`)"""
    start = None
    for i, l in enumerate(lines):
        if l.startswith("    localisation:") and (l[len("    localisation:"):].strip() + "#")[0] == "#":   # no value
            start = i
            break
    if start is None:
        raise SystemExit("no `    localisation:` block (inventory.params.localisation)")
    end = start + 1
    while end < len(lines) and (lines[end].startswith("      ") or not lines[end].strip() or lines[end].lstrip().startswith("#")):
        end += 1
    return start, end


def switch(cell, repo, check=False):
    p = os.path.join(repo, "config", "cells", f"{cell}.yaml")
    lines = open(p).read().split("\n")
    a, b = _loc_block(lines)
    changes = []
    for i in range(a, b):
        l = lines[i]
        if l.startswith("      preset: "):
            cur = l.split()[1]
            if cur == "common":
                return [f"{cell}: already preset common"], True
            if cur != LEGACY[cell]:
                raise SystemExit(f"{cell}: localisation preset is {cur!r}, expected {LEGACY[cell]!r}")
            lines[i] = "      preset: common                    # presets/common.yaml (B5 switch; was " + LEGACY[cell] + ")"
            changes.append(f"preset {cur} -> common")
        elif l.startswith("      label: "):
            lines[i] = f"      label: {cell}.localisation_common   # provenance name (B5 switch; was {l.split()[1]})"
            changes.append(f"label {l.split()[1]} -> {cell}.localisation_common")
    if not any(c.startswith("preset") for c in changes):
        raise SystemExit(f"{cell}: no `preset:` line in the localisation block")
    if cell == "jurkat":
        txt = "\n".join(lines)
        for old, new in JURKAT_SOURCES:
            if txt.count(old) != 1:
                raise SystemExit(f"jurkat: expected exactly one line containing {old!r}")
            txt = txt.replace(old, new)
            changes.append(f"localisation_sources: {old.strip()[:60]} -> {new.strip()[:60]}")
        lines = txt.split("\n")
    if not check:
        open(p, "w").write("\n".join(lines))
    return [f"{cell}: {c}" for c in changes], False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cells", nargs="+", choices=sorted(LEGACY))
    ap.add_argument("--repo", default=REPO)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    for c in a.cells:
        msgs, done = switch(c, a.repo, a.check)
        print("\n".join(msgs))
    return 0


if __name__ == "__main__":
    sys.exit(main())
