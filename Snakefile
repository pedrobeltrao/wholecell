# wholecell build workflow.  snakemake --cores N --config cell=<cell> [build_id=<id> | build_tag=<tag>]
# Usually launched by scripts/build_isolated.sh (frozen git worktree of HEAD under <build>/_code).
# The DAG is derived from the cell config (inventory rule sets, preplaced modules / serial builders, plugin steps,
# validation steps), never from a cell list.
# Every rule runs `python -m wholecell <stage> --cell C --build B --step S` with PYTHONHASHSEED=0.
include: "workflow/rules/common.smk"
include: "workflow/rules/build.smk"
if PLAN.get("pipeline_engine"):
    # reproduction configs only (config/cells/yeast_x5_port.yaml: the yeast X5 port run as a whole engine,
    # plugins/yeast/runner.py): DAG from the engine's step list
    include: "workflow/rules/engine.smk"
else:
    include: "workflow/rules/plugin_steps.smk"
    include: "workflow/rules/geometry.smk"
    include: "workflow/rules/inventory.smk"
    include: "workflow/rules/library.smk"
    include: "workflow/rules/preplaced.smk"
    include: "workflow/rules/pack.smk"
    include: "workflow/rules/validate.smk"
    include: "workflow/rules/export.smk"


rule all:
    input: FINAL
    default_target: True
