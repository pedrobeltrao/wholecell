# validation (owner: validation agent, wholecell/validate/stage.py): one rule per step of validate.stage.STEPS_ALL, then
# compare_v7 when the config declares validate.compare_reference. A step writes its result or <step>.skipped.json;
# the workflow tracks a sentinel validate/_smk/<step>.done written after the step succeeded.
_VAL_STAGE = bool(PLAN["validate_steps"])
if _VAL_STAGE:
    VAL_STEPS = [s for s in PLAN["validate_steps"] if s != "manifest"]
    _VD = lambda s: B("validate", "_smk", f"{s}.done")
    _VAL_IN = [PACK_SUMMARY, P_("preplaced_off.npy"), PACK_TABLE, GEOM]
    for _s in VAL_STEPS:
        rule:
            name: f"validate_{SAFE(_s)}"
            input: _VAL_IN + BEFORE(f"validate_{SAFE(_s)}")
            output: _VD(_s)
            threads: min(8, MAXW)
            log: LOG(f"validate_{SAFE(_s)}")
            params: step=_s
            shell: WC("validate", "{params.step}", "nproc={threads}") + " > {log} 2>&1 && touch {output}"

    rule validate_manifest:
        input: [_VD(s) for s in VAL_STEPS] + BEFORE("validate_manifest")
        output: _VD("manifest")
        log: LOG("validate_manifest")
        shell: WC("validate", "manifest") + " > {log} 2>&1 && touch {output}"

    FINAL += [_VD("manifest")]
    if PLAN["compare_reference"]:
        rule validate_compare_v7:
            input: _VD("manifest")
            output: _VD("compare_v7")
            log: LOG("validate_compare_v7")
            shell: WC("validate", "compare_v7") + " > {log} 2>&1 && touch {output}"

        FINAL += [_VD("compare_v7")]

    # common validation suite (plan s4, A8): after the manifest (and compare_v7), so neither lists its outputs; cell
    # checks without a common definition yet are plugin steps with `before: [validate_suite]` (the suite adapter reads them)
    rule validate_suite:
        input: [_VD("manifest")] + ([_VD("compare_v7")] if PLAN["compare_reference"] else []) + BEFORE("validate_suite")
        output: _VD("suite")
        threads: min(3, MAXW)
        log: LOG("validate_suite")
        shell: WC("validate", "suite", "nproc={threads}") + " > {log} 2>&1 && touch {output}"

    FINAL += [_VD("suite")]
