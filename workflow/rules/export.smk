# export (owner: validation/export agent, wholecell/export/stage.py). WIRED BUT NOT YET TESTED: runs the export stage's
# steps listed in config export.steps (default: [all]) after the pack summary and validation; sentinel per step.
if PLAN["export_steps"]:
    _EXP_STEPS = PLAN["export_steps"]
    _ED = lambda s: B("export", "_smk", f"{s}.done")
    _prev = [PACK_SUMMARY] + ([B("validate", "_smk", "manifest.done")] if _VAL_STAGE else [])
    for _s in _EXP_STEPS:
        rule:
            name: f"export_{SAFE(_s)}"
            input: list(_prev) + BEFORE(f"export_{SAFE(_s)}")
            output: _ED(_s)
            log: LOG(f"export_{SAFE(_s)}")
            params: step=_s
            shell: WC("export", "{params.step}") + " > {log} 2>&1 && touch {output}"
        _prev = [_ED(_s)]
    FINAL += [_ED(_EXP_STEPS[-1])]
