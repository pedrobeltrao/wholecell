"""python -m toolkits.secms run --config <run config yaml> --out <dir>
   python -m toolkits.secms combine --run <run dir> --label <cell label> --out <dir> [--spec <yaml>]   (conditions /
   datasets -> one estimate, combine.py; --spec: per-dataset weights and roles primary | fill; writes every table with
   the combined rows appended)
   (placeholders {data_root} ... resolved from the
site file; frozen-input names are not resolved here: give paths). Writes the tables + diagnostics.json + run_config.json."""
import argparse
import json
import os
import sys

import yaml


def main(argv=None):
    ap = argparse.ArgumentParser(prog="toolkits.secms")
    ap.add_argument("cmd", choices=["run", "combine"])
    ap.add_argument("--config")
    ap.add_argument("--run")
    ap.add_argument("--label")
    ap.add_argument("--min-fractions", type=int, default=0)
    ap.add_argument("--spec", help="combination spec YAML (label, datasets: [{source_id, weight, role}])")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if a.cmd == "combine":
        import pandas as pd
        from . import combine, pipeline
        if os.path.exists(a.out) and os.listdir(a.out):
            sys.exit(f"{a.out} exists and is not empty (write-once)")
        res = {k: pd.read_csv(os.path.join(a.run, f"{k}.tsv"), sep="\t", dtype={"replicate": str, "condition": str})
               for k in list(pipeline.TABLES) + ["profiles"]}
        spec = yaml.safe_load(open(a.spec)) if a.spec else None
        label = a.label or (spec or {}).get("label")
        if not label:
            sys.exit("--label (or label in --spec) required")
        out = combine.combine_tables(res, label, a.min_fractions, spec)
        pipeline.write(out, a.out)
        json.dump(dict(run=a.run, label=label, min_fractions=a.min_fractions, spec=spec),
                  open(os.path.join(a.out, "combine.json"), "w"), indent=1)
        print(f"wrote {a.out}")
        return
    from wholecell.config import load_site, _resolve
    from . import pipeline
    site = load_site()
    cfg = _resolve(yaml.safe_load(open(a.config)), site)
    if os.path.exists(a.out) and os.listdir(a.out):
        sys.exit(f"{a.out} exists and is not empty (write-once)")
    res = pipeline.run(cfg)
    pipeline.write(res, a.out)
    open(os.path.join(a.out, "diagnostics.json"), "w").write(pipeline.diagnostics_json(res["_diagnostics"]))
    json.dump(cfg, open(os.path.join(a.out, "run_config.json"), "w"), indent=1)
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
