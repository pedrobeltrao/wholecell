# Whole-pipeline engine (config pipeline_engine; e.g. yeast_x5 from plugins/yeast). One rule per engine step; the step
# list and its dependencies come from the plan (engine.steps(cfg)), never from a cell list. stage `engine` runs
# `python -m wholecell engine --step <id>`; other stages (pack, export, inventory) run the core stage of that name
# (pack dispatches to config pack.engine). Each step leaves a sentinel engine/_smk/<id>.done; post-processing
# (validation suite, export, comparison with the reference) is part of the DAG, so nothing needs polling.
_ESTEPS = PLAN["engine_steps"]
_EDONE = lambda s: B("engine", "_smk", f"{s}.done")
for _s in _ESTEPS:
    rule:
        name: f"engine_{SAFE(_s['id'])}"
        input: [FROZEN_OK] + [_EDONE(a) for a in _s["after"]]
        output: _EDONE(_s["id"])
        log: LOG(f"engine_{SAFE(_s['id'])}")
        threads: MAXW if _s.get("heavy") else 1
        resources: gpu=1 if _s.get("gpu") else 0
        params: cmd=WC(_s["stage"], _s["step"])
        shell: "{params.cmd} > {log} 2>&1 && mkdir -p $(dirname {output}) && touch {output}"
_ELEAVES = [s["id"] for s in _ESTEPS if not any(s["id"] in t["after"] for t in _ESTEPS)]
FINAL += [_EDONE(s) for s in _ELEAVES]
