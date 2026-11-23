# inventory: one rule per step of wholecell.inventory.stage.steps(cfg) (v0 -> assembly -> ingredients -> rule sets
# -> pack_table -> canonical [-> compare]); inputs/outputs come from the stage's own plan, so pre-pass outputs a rule
# reads (e.g. preplaced/lamina/consumption.json for Jurkat rule v7) order the DAG automatically.
INV_PLAN = PLAN["inventory"]
for _step, _spec in INV_PLAN.items():
    rule:
        name: f"inventory_{SAFE(_step)}"
        input: [B(x) for x in _spec["inputs"]] + [FROZEN(n) for n in _spec["frozen"]] + [FROZEN_OK, GEOM] + BEFORE(f"inventory_{SAFE(_step)}")
        output: [B(x) for x in _spec["outputs"]]
        log: LOG(f"inventory_{SAFE(_step)}")
        params: step=_step
        shell: WC("inventory", "{params.step}") + " > {log} 2>&1"

PACK_TABLE = B("inventory", "pack_table.tsv")
FINAL += [B(x) for x in INV_PLAN[list(INV_PLAN)[-1]]["outputs"]]
