"""Per-cell localisation overrides: curated revisions as DATA (a TSV, one row per change, with reason and citation), applied
by one function.

Override table columns (OVERRIDE_COLUMNS; '#' lines are provenance comments):
  rule_set     the rule set that applies the row (e.g. jurkat.v4)
  rule_id      the change-table rule id the row produces (e.g. v4.f_curated)
  key_type     species (first species with that ing_id) | gene (every single-gene species of that gene)
  key          ing_id or gene
  field        comp_weights | transmembrane | n_tm | annotation (the ordered key list of the combiner's ordered-annotation
               steps, op set, value = JSON list; applied inside the combiner by step annotation_overrides, e.g. yeast LOC_FIX)
  op           set (value = the new shares / value) | drop (value = {"keys": [...], "fallback": {...}}: the keys removed,
               fallback shares when nothing is left) | set_keep_other (value = {"shares": {...}, "keep": f, "exclude":
               [...]}: shares x (1 - f) + f x the current weights outside `exclude`, renormalised among themselves)
  value        JSON
  reason       why (biology), shown in the change tables
  citation     literature / DECISION / data source

Consecutive rows with the same (rule_id, key) are one group: for each species the key resolves to, the group's rows are
applied in order (e.g. transmembrane, n_tm, then comp_weights of one protein).

Frozen-input override tables (generated at build time from declared frozen inputs, e.g. the v4/v5 pack tables of the
Jurkat legacy_v7 preset) have the same columns plus `genes` (the genes that selected the species); they are written next
to the rule set's change table and applied by the same function.
"""
from __future__ import annotations

import difflib
import json

import pandas as pd

OVERRIDE_COLUMNS = ["rule_set", "rule_id", "key_type", "key", "field", "op", "value", "reason", "citation"]
KEY_TYPES = ("species", "gene")
FIELDS = ("comp_weights", "transmembrane", "n_tm", "annotation")
OPS = {"comp_weights": ("set", "drop", "set_keep_other"), "transmembrane": ("set",), "n_tm": ("set",), "annotation": ("set",)}


class OverrideError(ValueError):
    pass


def _near(v, choices):
    m = difflib.get_close_matches(str(v), [str(c) for c in choices], n=1, cutoff=0.5)
    return f"; did you mean '{m[0]}'?" if m else ""


def check(T: pd.DataFrame, where: str) -> pd.DataFrame:
    miss = [c for c in OVERRIDE_COLUMNS if c not in T.columns]
    if miss:
        raise OverrideError(f"override table {where}: missing columns {miss} (required: {OVERRIDE_COLUMNS})")
    for n, r in enumerate(T.itertuples(index=False), start=1):
        at = f"override table {where} row {n} ({r.rule_id} {r.key})"
        if r.key_type not in KEY_TYPES:
            raise OverrideError(f"{at}: key_type '{r.key_type}' (allowed: {list(KEY_TYPES)}){_near(r.key_type, KEY_TYPES)}")
        if r.field not in FIELDS:
            raise OverrideError(f"{at}: field '{r.field}' (allowed: {list(FIELDS)}){_near(r.field, FIELDS)}")
        if r.op not in OPS[r.field]:
            raise OverrideError(f"{at}: op '{r.op}' for field {r.field} (allowed: {list(OPS[r.field])}){_near(r.op, OPS[r.field])}")
        try:
            v = json.loads(r.value)
        except ValueError as e:
            raise OverrideError(f"{at}: value is not JSON ({e}): {r.value!r}")
        if r.op == "set" and r.field == "comp_weights" and not (isinstance(v, dict) and v):
            raise OverrideError(f"{at}: comp_weights set needs a non-empty JSON object of shares")
        if r.field == "annotation" and not (isinstance(v, list) and v and all(isinstance(x, str) for x in v)):
            raise OverrideError(f"{at}: annotation set needs a non-empty JSON list of keys (in rank order)")
        if r.op == "drop" and not (isinstance(v, dict) and "keys" in v and "fallback" in v):
            raise OverrideError(f"{at}: drop needs {{\"keys\": [...], \"fallback\": {{...}}}}")
        if r.op == "set_keep_other" and not (isinstance(v, dict) and {"shares", "keep", "exclude"} <= set(v)):
            raise OverrideError(f"{at}: set_keep_other needs {{\"shares\": {{...}}, \"keep\": f, \"exclude\": [...]}}")
        if not str(r.reason).strip():
            raise OverrideError(f"{at}: reason is empty (every override needs a reason)")
    return T


def new_value(op: str, current, value):
    """the new field value of one override row (comp_weights: un-normalised dict; the rule set's setter normalises / rounds)"""
    if op == "set":
        return dict(value) if isinstance(value, dict) else value
    if op == "drop":
        drop = set(value["keys"])
        return {k: v for k, v in current.items() if k not in drop} or dict(value["fallback"])
    if op == "set_keep_other":
        keep, ex = value["keep"], set(value["exclude"])
        oth = {k: v for k, v in current.items() if k not in ex}
        if not oth or keep <= 0:
            return dict(value["shares"])
        so = sum(oth.values())
        return {**{k: s * (1.0 - keep) for k, s in value["shares"].items()}, **{k: keep * v / so for k, v in oth.items()}}
    raise OverrideError(f"op '{op}' (allowed: set, drop, set_keep_other)")


def groups(T: pd.DataFrame):
    """consecutive rows with the same (rule_id, key) -> lists of rows"""
    out, cur, last = [], [], None
    for r in T.itertuples(index=False):
        k = (r.rule_id, r.key_type, r.key)
        if cur and k != last:
            out.append(cur)
            cur = []
        cur.append(r)
        last = k
    if cur:
        out.append(cur)
    return out


def apply(T: pd.DataFrame, find, get, put) -> int:
    """apply an override table: THE one function for every override row.

    find(key_type, key) -> species row positions; get(i, field) -> the current value (comp_weights as a dict);
    put(i, row, new) records and sets the value (the rule set's change-table format). Returns the number of puts."""
    n = 0
    for grp in groups(T):
        r0 = grp[0]
        for i in find(r0.key_type, r0.key):
            for r in grp:
                put(i, r, new_value(r.op, get(i, r.field), json.loads(r.value)))
                n += 1
    return n


def select(T: pd.DataFrame, rule_set: str, rule_ids) -> pd.DataFrame:
    rule_ids = [rule_ids] if isinstance(rule_ids, str) else list(rule_ids)
    sub = T[(T.rule_set == rule_set) & T.rule_id.isin(rule_ids)]
    return sub.reset_index(drop=True)


def frame(rows) -> pd.DataFrame:
    """generated override rows (dicts) -> table with OVERRIDE_COLUMNS first"""
    D = pd.DataFrame(rows)
    if not len(D):
        return pd.DataFrame(columns=OVERRIDE_COLUMNS + ["genes"])
    return D[OVERRIDE_COLUMNS + [c for c in D.columns if c not in OVERRIDE_COLUMNS]]
