# plugin steps: cell-specific steps of the normal stage list (config plugin_steps; wholecell/plugin_steps.py). One rule
# <stage>_<id> per step: inputs = its declared build-relative files + the sentinels of the plugin steps in `after`;
# outputs = its declared files + the sentinel _smk/<stage>_<id>.done; the step gets nproc=<threads>. Core rules wait for it through BEFORE(<rule>).
_PIDS = {s["id"]: s for s in PLUG}
for _p in PLUG:
    rule:
        name: f"{_p['stage']}_{SAFE(_p['id'])}"
        input: [FROZEN_OK] + [B(x) for x in _p["inputs"]] + [PLUG_DONE(_PIDS[a]) for a in _p["after"]]
        output: [B(x) for x in _p["outputs"]] + [PLUG_DONE(_p)]
        log: LOG(f"{_p['stage']}_{SAFE(_p['id'])}")
        threads: MAXW if _p["threads"] == "max" else min(int(_p["threads"]), MAXW)
        resources: gpu=1 if _p["gpu"] else 0
        params: cmd=WC(_p["stage"], f"plugin:{_p['id']}"), done=PLUG_DONE(_p)
        shell: "{params.cmd} nproc={threads} > {log} 2>&1 && mkdir -p $(dirname {params.done}) && touch {params.done}"
FINAL += [PLUG_DONE(s) for s in PLUG]
