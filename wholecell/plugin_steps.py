"""Plugin steps: cell-specific steps inside the normal stage list (config `plugin_steps`).

A cell whose biology (or whose not-yet-ported legacy code) needs a step the core does not have lists it here; the
step is a registered function (kind `plugin_step`, fn(cfg, bdir, params)) of a plugin named in the config.

    plugin_steps:
      - id: yeast_protein_table             # unique; Snakemake rule <stage>_<id>
        stage: inventory                     # stage it belongs to (one of STAGES); docs/stages/<stage>.md lists it
        run: yeast.protein_table             # registered plugin_step
        inputs: [inventory/v0/inventory_v0.tsv]        # build-relative files it reads (outputs of other rules)
        outputs: [inventory/yeast/protein_table.parquet]   # build-relative files it writes
        before: [inventory_assembly]         # rules that must wait for it (their inputs get its sentinel)
        params: {...}                        # handed to the function
        threads: 1                           # or "max" (= the build's worker count); gpu: true reserves a GPU

Run: `python -m wholecell <stage> --cell C --build B --step plugin:<id>`. The workflow (workflow/rules/plugin_steps.smk)
writes the sentinel _smk/<stage>_<id>.done after the step succeeded; `before` hands that sentinel to the named rules.
"""
from __future__ import annotations

STAGES = ("geometry", "inventory", "library", "preplaced", "pack", "validate", "export")
KEYS = {"id", "stage", "run", "inputs", "outputs", "before", "after", "params", "threads", "gpu", "note"}


def specs(cfg) -> list[dict]:
    """the config's plugin steps, checked and normalised (order kept)"""
    out, seen = [], set()
    for i, s in enumerate(cfg.get("plugin_steps") or []):
        bad = set(s) - KEYS
        if bad:
            raise ValueError(f"plugin_steps[{i}]: unknown keys {sorted(bad)} (allowed: {sorted(KEYS)})")
        for k in ("id", "stage", "run"):
            if not s.get(k):
                raise ValueError(f"plugin_steps[{i}]: '{k}' is required")
        if s["stage"] not in STAGES:
            raise ValueError(f"plugin_steps[{i}] {s['id']}: stage '{s['stage']}' is not one of {STAGES}")
        if s["id"] in seen:
            raise ValueError(f"plugin_steps: id '{s['id']}' listed twice")
        seen.add(s["id"])
        th = s.get("threads", 1)
        out.append(dict(id=s["id"], stage=s["stage"], run=s["run"], inputs=list(s.get("inputs") or []),
                        outputs=list(s.get("outputs") or []), before=list(s.get("before") or []),
                        after=list(s.get("after") or []), threads=th if th == "max" else int(th), gpu=bool(s.get("gpu", False))))
    ids = {s["id"] for s in out}
    for s in out:
        for a in s["after"]:
            if a not in ids:
                raise ValueError(f"plugin_steps {s['id']}: after '{a}' is not a plugin step id")
    return out


def run_step(cfg, bdir, stage, sid, **kw):
    from . import registry
    spec = next((s for s in (cfg.get("plugin_steps") or []) if s.get("id") == sid), None)
    if spec is None:
        raise SystemExit(f"plugin step '{sid}' is not in the config plugin_steps")
    if spec["stage"] != stage:
        raise SystemExit(f"plugin step '{sid}' belongs to stage '{spec['stage']}', not '{stage}'")
    fn = registry.get("plugin_step", spec["run"])
    params = dict(spec.get("params") or {})
    params.update(kw)
    print(f"[{stage}] plugin step {sid} ({spec['run']})", flush=True)
    return fn(cfg, bdir, params)
