"""Chromatin metrics of a frozen chromatin pre-pass (criteria V5): the metrics are not recomputed; they are copied
verbatim from the frozen validation outputs of the chromatin build after the pre-pass file checksum is verified.

Config validate.chromatin:
  preplaced_input frozen-input name (or path) of the chromatin pre-pass used by the pack (Jurkat: v8 rep5)
  expected_md5    its md5 (the copy is refused if it differs)
  metric_inputs   frozen-input names (or paths) of the metric files to copy (JSON parsed, TSV as row records)
  headline        {label: metric-row name in the metrics table} values pulled out for reports
Output: <build>/validate/chromatin_frozen.json
"""
from __future__ import annotations

import json
import os

import pandas as pd

from . import common


def _path(cfg, name):
    fi = cfg.get("frozen_inputs", {}) or {}
    return fi[name]["path"] if name in fi else name


def run(cfg, bdir, out_path, overwrite=False, preplaced_path=None):
    cc = (cfg.get("validate", {}) or {}).get("chromatin", {}) or {}
    p = preplaced_path or _path(cfg, cc["preplaced_input"])
    act = common.md5(p)
    exp = cc.get("expected_md5")
    if exp and act != exp:
        raise RuntimeError(f"chromatin pre-pass checksum {act} != expected {exp} ({p}); metrics not copied")
    files = {}
    rows = None
    for name in cc.get("metric_inputs", []) or []:
        f = _path(cfg, name)
        ent = dict(path=f, md5=common.md5(f))
        if f.endswith(".json"):
            ent["content"] = json.load(open(f))
        elif f.endswith((".tsv", ".txt")):
            T = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
            ent["content"] = T.to_dict("records")
            if "metric" in T.columns and rows is None:
                rows = {r["metric"]: r for r in ent["content"]}
        files[name] = ent
    head = {}
    for label, metric in (cc.get("headline", {}) or {}).items():
        if rows and metric in rows:
            head[label] = rows[metric]
    out = dict(prepass=dict(path=p, md5_expected=exp, md5_actual=act, ok=(exp is None or exp == act)), headline=head, files=files,
               note="copied verbatim from the frozen chromatin validation outputs; not recomputed (criteria V5)")
    common.write_json(out_path, out, overwrite)
    return out
