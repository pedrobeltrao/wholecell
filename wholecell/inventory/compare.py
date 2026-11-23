"""Port vs legacy comparison of inventory tables (criteria B3/B4 of docs/reference/history/v7_reproduction_criteria.md).

Config `inventory.verify`: {tables: {port_rel: legacy_path}, ingredients: {port_rel: legacy_path}} (port_rel relative to
the build's inventory/ folder). For every pair: md5 of both files and byte identity; for tables, row/column identity and
the first differing rows; for ingredient (species) tables additionally the B3 metrics: species set difference,
max |log2 copies ratio| (copies + 1 pseudo-count; zero-copy species flagged), max |comp_weights difference| over the
union of keys, and differences in the other fields. Result: inventory/compare_legacy.json.
"""
from __future__ import annotations

import json
import os
import re

import numpy as np
import pandas as pd

from .. import build

TOL_LOG2 = 0.01
TOL_W = 1e-6


def _table_diff(a, b, key=None, nshow=5):
    A = pd.read_csv(a, sep="\t", dtype=str, keep_default_na=False)
    B = pd.read_csv(b, sep="\t", dtype=str, keep_default_na=False)
    out = dict(rows_port=len(A), rows_legacy=len(B), columns_equal=list(A.columns) == list(B.columns))
    if not out["columns_equal"]:
        out["columns_port"], out["columns_legacy"] = list(A.columns), list(B.columns)
    if len(A) == len(B) and out["columns_equal"]:
        ne = (A.values != B.values)
        rows = np.nonzero(ne.any(1))[0]
        out["rows_differing"] = int(len(rows))
        out["cells_differing"] = int(ne.sum())
        out["columns_differing"] = {c: int(ne[:, j].sum()) for j, c in enumerate(A.columns) if ne[:, j].any()}
        out["examples"] = [{c: [A.iat[r, j], B.iat[r, j]] for j, c in enumerate(A.columns) if ne[r, j]} | {"_row": int(r)}
                           for r in rows[:nshow]]
    return out


_KV = re.compile(r"^[A-Za-z_]+:[-0-9.e]+(;[A-Za-z_]+:[-0-9.e]+)*$")


def _norm(v):
    """order-insensitive form of dict-like cells: JSON objects (comp_weights) and 'k:v;k:v' strings (compartments)"""
    if isinstance(v, str) and v.startswith("{") and v.endswith("}"):
        try:
            return json.dumps(sorted(json.loads(v).items()))
        except ValueError:
            return v
    if isinstance(v, str) and _KV.match(v):
        return ";".join(sorted(v.split(";")))
    return v


def _content_diff(a, b, rtol=1e-9):
    """content comparison: same rows/columns; dict-like cells compared order-insensitively; numeric cells with a relative
    tolerance. Separates 'order only' differences (set iteration order of the legacy process, random hash seed) and
    floating-point last-digit differences from real ones."""
    A = pd.read_csv(a, sep="\t", dtype=str, keep_default_na=False)
    B = pd.read_csv(b, sep="\t", dtype=str, keep_default_na=False)
    if A.shape != B.shape or list(A.columns) != list(B.columns):
        return dict(comparable=False)
    out = dict(comparable=True, order_only_cells={}, float_cells={}, max_rel_diff={}, real_cells={}, real_examples=[])
    for c in A.columns:
        x, y = A[c].values, B[c].values
        ne = x != y
        if not ne.any():
            continue
        xn = np.array([_norm(v) for v in x[ne]], dtype=object); yn = np.array([_norm(v) for v in y[ne]], dtype=object)
        oo = xn == yn
        if oo.any():
            out["order_only_cells"][c] = int(oo.sum())
        rest = np.nonzero(ne)[0][~oo]
        if not len(rest):
            continue
        xf = pd.to_numeric(pd.Series(x[rest]), errors="coerce").values
        yf = pd.to_numeric(pd.Series(y[rest]), errors="coerce").values
        num = ~np.isnan(xf) & ~np.isnan(yf)
        rel = np.abs(xf - yf) / np.maximum(np.abs(yf), 1e-300)
        fl = num & (rel <= rtol)
        if fl.any():
            out["float_cells"][c] = int(fl.sum()); out["max_rel_diff"][c] = float(rel[fl].max())
        real = rest[~fl]
        if len(real):
            out["real_cells"][c] = int(len(real))
            for r in real[:3]:
                out["real_examples"].append(dict(row=int(r), column=c, port=x[r][:300], legacy=y[r][:300]))
    out["identical_content"] = not out["real_cells"]
    if out["real_cells"]:
        # row-order-insensitive: compare the multisets of normalised rows (numbers rounded to 12 significant digits)
        def canon(T):
            def cell(v):
                v = _norm(v)
                try:
                    return f"{float(v):.12g}" if v not in ("", "True", "False") else v
                except (TypeError, ValueError):
                    return v
            return sorted("\t".join(cell(v) for v in row) for row in T.itertuples(index=False, name=None))
        ca, cb = canon(A), canon(B)
        out["rows_multiset_equal"] = ca == cb
        if not out["rows_multiset_equal"]:
            import collections
            d = collections.Counter(ca); d.subtract(collections.Counter(cb))
            out["rows_unmatched_port"] = sum(v for v in d.values() if v > 0)
            out["rows_unmatched_legacy"] = -sum(v for v in d.values() if v < 0)
            out["unmatched_examples_port"] = [k[:400] for k, v in d.items() if v > 0][:3]
            out["unmatched_examples_legacy"] = [k[:400] for k, v in d.items() if v < 0][:3]
        out["identical_content_any_row_order"] = out["rows_multiset_equal"]
    else:
        out["identical_content_any_row_order"] = True
    return out


def _w(s):
    return json.loads(s) if isinstance(s, str) and s else {}


def species_metrics(a, b):
    A = pd.read_csv(a, sep="\t").set_index("ing_id")
    B = pd.read_csv(b, sep="\t").set_index("ing_id")
    sa, sb = set(A.index), set(B.index)
    both = sorted(sa & sb)
    ca, cb = A.copies.reindex(both).astype(float), B.copies.reindex(both).astype(float)
    l2 = np.abs(np.log2((ca + 1) / (cb + 1)))
    wd = []
    for i in both:
        wa, wb = _w(A.at[i, "comp_weights"]), _w(B.at[i, "comp_weights"])
        wd.append(max([abs(wa.get(k, 0.0) - wb.get(k, 0.0)) for k in set(wa) | set(wb)] or [0.0]))
    wd = pd.Series(wd, index=both)
    other = {}
    for c in ("kind", "genes", "stoich", "transmembrane", "n_tm", "struct_src", "signal_frac", "mito_om_frac", "mass_da"):
        if c in A.columns and c in B.columns:
            x, y = A[c].reindex(both).astype(str), B[c].reindex(both).astype(str)
            other[c] = int((x != y).sum())
    res = dict(n_port=len(sa), n_legacy=len(sb), species_only_port=sorted(sa - sb)[:50], n_only_port=len(sa - sb),
               species_only_legacy=sorted(sb - sa)[:50], n_only_legacy=len(sb - sa),
               max_abs_log2_copies=float(l2.max()) if len(l2) else 0.0,
               n_log2_gt_tol=int((l2 > TOL_LOG2).sum()), worst_log2=l2.sort_values().tail(5).round(6).to_dict(),
               max_abs_comp_weights_diff=float(wd.max()) if len(wd) else 0.0, n_w_gt_tol=int((wd > TOL_W).sum()),
               worst_w=wd.sort_values().tail(5).to_dict(), other_field_differences=other,
               mass_pg_port=float(A.mass_pg.sum()), mass_pg_legacy=float(B.mass_pg.sum()))
    res["B3_pass"] = bool(res["n_only_port"] == 0 and res["n_only_legacy"] == 0 and res["max_abs_log2_copies"] <= TOL_LOG2
                          and res["max_abs_comp_weights_diff"] <= TOL_W)
    return res


def run(ctx):
    v = ctx.inv.get("verify") or {}
    out = dict(tables={}, ingredients={})
    for kind in ("tables", "ingredients"):
        for rel, leg in (v.get(kind) or {}).items():
            port = ctx.path(rel)
            r = dict(port=port, legacy=leg)
            if not os.path.exists(port) or not os.path.exists(leg):
                r["status"] = "missing " + ("port" if not os.path.exists(port) else "legacy")
                out[kind][rel] = r; continue
            r["md5_port"], r["md5_legacy"] = build.file_checksum(port), build.file_checksum(leg)
            r["identical"] = r["md5_port"] == r["md5_legacy"]
            if not r["identical"]:
                try:
                    r["diff"] = _table_diff(port, leg)
                    r["content"] = _content_diff(port, leg)
                except Exception as e:                     # noqa: BLE001
                    r["diff"] = f"unreadable: {e}"
            if kind == "ingredients":
                r["B3"] = species_metrics(port, leg)
            out[kind][rel] = r
            ctx.say(f"{rel}: {'IDENTICAL' if r['identical'] else 'differs'}"
                    + ("" if r["identical"] else f"; content identical {r.get('content', {}).get('identical_content')}"
                           f" (any row order: {r.get('content', {}).get('identical_content_any_row_order')})")
                    + (f"; B3 pass {r['B3']['B3_pass']}" if kind == "ingredients" else ""))
    out["summary"] = dict(
        identical=sorted(k for d in (out["tables"], out["ingredients"]) for k, r in d.items() if r.get("identical")),
        identical_content=sorted(k for d in (out["tables"], out["ingredients"]) for k, r in d.items()
                                 if r.get("identical") or r.get("content", {}).get("identical_content_any_row_order")),
        differing=sorted(k for d in (out["tables"], out["ingredients"]) for k, r in d.items() if r.get("identical") is False),
        B3_all_pass=all(r.get("B3", {}).get("B3_pass", False) for r in out["ingredients"].values()) if out["ingredients"] else None)
    b3 = [r["B3"] for r in out["ingredients"].values() if "B3" in r]
    if b3:
        out["summary"]["B3_max_abs_log2_copies"] = max(x["max_abs_log2_copies"] for x in b3)
        out["summary"]["B3_max_abs_comp_weights_diff"] = max(x["max_abs_comp_weights_diff"] for x in b3)
        out["summary"]["B3_species_set_diff"] = sum(x["n_only_port"] + x["n_only_legacy"] for x in b3)
    out["notes"] = v.get("notes", [])
    p = ctx.out("compare_legacy.json")
    build.atomic_save(p, lambda t: json.dump(out, open(t, "w"), indent=1, default=str))
