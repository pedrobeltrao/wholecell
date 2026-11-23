"""Present a computed species table in the row / key order of a reference table (exact reproduction of a legacy build).

Why: the legacy Jurkat scripts ran with a random PYTHONHASHSEED (set iteration order of strings -> order of equal-weight
keys in comp_weights JSON and of tied candidates) and sorted with a non-stable quicksort on another CPU (row order among
equal sort keys). Those orders carry no information but they change bytes, library indices and float summation order in
the packer. With `inventory.pack_table.order_reference: <frozen input>` the pack_table step keeps the port's CONTENT and
only borrows the reference's row order and, for dict cells whose parsed content is exactly equal, the reference's key
order. Anything that differs in content is kept as computed and reported (so a content difference can never be hidden).
"""
from __future__ import annotations

import json

import pandas as pd


def _d(s):
    try:
        return json.loads(s) if isinstance(s, str) and s.startswith("{") else None
    except ValueError:
        return None


def reorder(port_path: str, ref_path: str, key: str = "ing_id", dict_cols=("comp_weights",)):
    P = pd.read_csv(port_path, sep="\t", dtype=str, keep_default_na=False)
    R = pd.read_csv(ref_path, sep="\t", dtype=str, keep_default_na=False)
    rep = dict(rows=len(P), ref_rows=len(R), columns_equal=list(P.columns) == list(R.columns))
    sp, sr = set(P[key]), set(R[key])
    rep["only_port"], rep["only_ref"] = sorted(sp - sr), sorted(sr - sp)
    if sp != sr or P[key].duplicated().any() or R[key].duplicated().any():
        rep["applied"] = False
        return P, rep
    pos = {k: i for i, k in enumerate(R[key])}
    rep["rows_moved"] = int(sum(pos[k] != i for i, k in enumerate(P[key])))
    P = P.iloc[sorted(range(len(P)), key=lambda i: pos[P[key].iat[i]])].reset_index(drop=True)
    Rk = R.set_index(key)
    n_key, n_diff, diffs = 0, 0, []
    for c in P.columns:
        if c == key or c not in Rk.columns:
            continue
        rv = Rk[c].reindex(P[key]).values
        for i, (a, b) in enumerate(zip(P[c].values, rv)):
            if a == b:
                continue
            if c in dict_cols and _d(a) is not None and _d(a) == _d(b) and list(_d(a)) != list(_d(b)):
                P.iat[i, P.columns.get_loc(c)] = b
                n_key += 1
            else:
                n_diff += 1
                if len(diffs) < 20:
                    diffs.append(dict(species=P[key].iat[i], column=c, port=a, reference=b))
    rep.update(applied=True, dict_key_order_borrowed=n_key, content_cells_differing=n_diff, examples=diffs)
    return P, rep
