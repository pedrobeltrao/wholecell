"""Geometry rule tables: read a rule YAML (config/geometry/rules/<layout>.yaml) and compile it for cellgeom_eval_v1.

A rule table has a membrane table (how each membrane SDF is derived) and an ordered label rule (first matching row
wins); format and semantics: docs/reference/cellgeom_v1.md §4.3. compile_tables() resolves names (fields -> slots,
membranes / labels -> engine codes, param expressions -> numbers) and returns the numeric rows that the engine stores
in G. parse_condition() is shared with the pure-python reference evaluator (wholecell/geometry/eval_ref.py).
"""
from __future__ import annotations

import hashlib
import os
import re

import numpy as np
import yaml

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_COND = re.compile(r"^\s*(\S+)\s*([<>])\s*(.+?)\s*$")


def default_rules_path(layout):
    """the repo rule table of a package layout (manifest engine.name)"""
    return os.path.join(REPO, "config", "geometry", "rules", f"{layout}.yaml")


def load_rules(path):
    r = yaml.safe_load(open(path))
    r["_path"] = os.path.abspath(path)
    r["_sha256"] = hashlib.sha256(open(path, "rb").read()).hexdigest()
    return r


def value(v, params):
    """a number, or an expression over the scalar engine params (e.g. '-PERI_W-PORE_LW')"""
    if isinstance(v, (int, float)):
        return float(v)
    env = {k: float(x) for k, x in params.items() if not isinstance(x, (list, tuple, dict))}
    return float(eval(str(v), {"__builtins__": {}}, env))


def parse_condition(s):
    """'<term> <op> <value>' -> (term, op, value expression)"""
    m = _COND.match(s)
    if not m:
        raise ValueError(f"label rule condition {s!r} is not '<term> <op> <value>'")
    return m.group(1), m.group(2), m.group(3)


def active_rows(rules, present):
    """label rows kept for a package: rows whose `requires` are all present"""
    return [r for r in rules["labels"] if all(x in present for x in r.get("requires", []))]


def package_present(man):
    """items a package has: field:<name> (non-empty), section:<name>, clip:<name>"""
    out = {"field:" + n for n, f in man["fields"].items() if not f.get("empty")}
    out |= {"section:" + n for n in man.get("sections", {})}
    out |= {"clip:" + cp["name"] for cp in man.get("domain", {}).get("clip_planes", [])}
    return out


def _vec(v, params, n=None):
    out = [value(x, params) for x in (v if isinstance(v, (list, tuple)) else [v])]
    if n is not None and len(out) != n:
        raise ValueError(f"expected {n} values, got {len(out)}")
    return out


def analytic_data(spec, params):
    """analytic shape spec {kind, centre, radius | axes + rot} -> (kind code, data values) (engine A_* layout)"""
    from wholecell.geometry.engines import cellgeom_eval_v1 as ev
    kind, n = ev.ANALYTIC[spec["kind"]]
    if spec["kind"] in ("sphere", "ball"):
        d = _vec(spec.get("centre", [0, 0, 0]), params, 3) + [value(spec["radius"], params)]
    else:
        d = _vec(spec["centre"], params, 3) + _vec(spec["axes"], params, 3) + _vec(spec["rot"], params, 9)
    assert len(d) == n
    return kind, d


def membrane_names(rules):
    """membranes defined by the table, in engine-code order (smooth_union rows modify, they define none)"""
    return [m["name"] for m in rules["membranes"] if m["kind"] != "smooth_union"]


def compile_tables(rules, params, present, field_slots=None, label_codes=None, mem_codes=None, prim_classes=None):
    """-> engine Tables (membrane rows, label rows, field fallbacks, data, clearance) for cellgeom_eval_v1.build_G.
    field_slots {field: slot}, label_codes {vocabulary name: engine code}, mem_codes {membrane: engine code}: from the
    package manifest; None = the engine defaults (sdf_sections_v1 layout). prim_classes: [class name] (index = class
    code, manifest prims.classes). Iterating the result gives (mem_rows, label_rows)."""
    from wholecell.geometry.engines import cellgeom_eval_v1 as ev
    field_slots = field_slots or {n: k for k, n in enumerate(ev.FIELD_NAMES)}
    mem_codes = mem_codes or {n: k for k, n in enumerate(ev.MEM_NAMES)}
    label_codes = label_codes or DEFAULT_LABEL_CODES
    pcls = {n: k for k, n in enumerate(prim_classes or [])}
    mems = rules["membranes"]
    names = membrane_names(rules)
    for i, n in enumerate(names):
        if mem_codes.get(n) != i:
            raise ValueError(f"rule membrane {i} ({n}) != package engine code {mem_codes.get(n)}")
    if len(names) != len(mem_codes):
        raise ValueError(f"rule table defines {len(names)} membranes, package has {len(mem_codes)}")
    mrow = {n: i for i, n in enumerate(names)}
    data = []

    def put(vals):
        data.extend(vals)
        return len(data) - len(vals)

    def cls(name):
        if name not in pcls:
            raise ValueError(f"prim class {name!r} not in the package prims.classes {list(pcls)}")
        return pcls[name]

    M = np.zeros((len(mems), ev.MROW))
    defined = 0
    shell_row = {}
    for r, m in enumerate(mems):
        kind = ev.MEM_KINDS[m["kind"]]

        def earlier(name):
            if name not in mrow or mrow[name] >= defined:
                raise ValueError(f"membrane row {r} uses {name}: only membranes defined above may be used")
            return mrow[name]
        M[r, 0] = kind
        if kind == ev.K_FIELD:
            M[r, 1] = field_slots[m["field"]]; M[r, 2] = float(m.get("sign", 1))
            off = m.get("offset", 0)       # a value, or [a, b] = a - b evaluated by the kernel (cell3d float ops)
            M[r, 3:5] = [value(x, params) for x in off] if isinstance(off, list) else [value(off, params), 0.0]
        elif kind == ev.K_WALL:
            M[r, 1] = earlier(m["of"])
        elif kind == ev.K_SHELL_PORES:
            M[r, 1] = field_slots[m["field"]]
            M[r, 3:7] = [value(m["width"], params), value(m["pore_r"], params), value(m["zwin"], params),
                         value(m.get("field_offset", 0), params)]
            shell_row[m["name"]] = r
        elif kind == ev.K_INNER_CRISTAE:
            M[r, 1] = earlier(m["of"]); M[r, 3] = value(m["offset"], params)
        elif kind == ev.K_GOLGI_CIS:
            a, b = m["clear_of"]
            M[r, 1] = field_slots[m["field"]]; M[r, 3] = earlier(a); M[r, 4] = earlier(b)
        elif kind == ev.K_PRIMS:
            M[r, 1] = cls(m["prims"])
            base = m.get("base_field")
            M[r, 2] = field_slots[base] if base and ("field:" + base) in present else -1
            if m.get("ref"):
                M[r, 3] = earlier(m["ref"]); M[r, 4] = value(m["ref_radius"], params)
            else:
                M[r, 3] = -1
        elif kind == ev.K_PRIMS_DIFF:
            M[r, 1] = cls(m["prims"]); M[r, 2] = cls(m["minus"]); M[r, 3] = float(m.get("sign", 1))
        elif kind == ev.K_SMOOTH_UNION:
            M[r, 1] = earlier(m["into"]); M[r, 2] = cls(m["prims"])
            M[r, 3] = value(m["k"], params); M[r, 4] = float(m.get("side", 1)); M[r, 5] = 0.25 / M[r, 3]
        if kind != ev.K_SMOOTH_UNION:
            defined += 1
    # field fallbacks: analytic shapes for fields the package lacks
    FB = []
    for fname, spec in (rules.get("fields") or {}).items():
        fb = spec.get("fallback")
        if fb and ("field:" + fname) not in present:
            kind, d = analytic_data(fb, params)
            FB.append([field_slots[fname], kind, put(d), 0])
    shapes = {}
    for sname, spec in (rules.get("shapes") or {}).items():
        kind, d = analytic_data(spec, params)
        shapes[sname] = put([kind] + d)
    rows = active_rows(rules, present)
    L = np.zeros((len(rows), ev.RROW))
    for r, row in enumerate(rows):
        if row["label"] not in label_codes:
            raise ValueError(f"label rule names {row['label']!r}, which has no engine code in the package vocabulary")
        conds = row.get("if", [])
        if len(conds) > ev.NCOND:
            raise ValueError(f"label row {row['label']}: at most {ev.NCOND} conditions")
        L[r, 0] = label_codes[row["label"]]; L[r, 1] = len(conds)
        for c, s in enumerate(conds):
            term, op, v = parse_condition(s)
            if term in mrow:
                t, a = ev.T_MEM, mrow[term]
            elif term.startswith("field:"):
                t, a = ev.T_FIELD, field_slots[term[6:]]
            elif term == "face":
                t, a = ev.T_FACE, 0
            elif term == "pore":
                t, a = ev.T_PORE, next(iter(shell_row.values()))
            elif term == "speckle":
                t, a = ev.T_SPECKLE, 0
            elif term.startswith("shape:"):
                t, a = ev.T_SHAPE, shapes[term[6:]]
            else:
                raise ValueError(f"unknown label rule term {term!r}")
            L[r, 2 + 4 * c:6 + 4 * c] = [t, a, ev.OP_LT if op == "<" else ev.OP_GT, value(v, params)]
    if not rows or rows[-1].get("if"):
        raise ValueError("the label rule must end with an unconditional row")
    cl = rules.get("clearance") or {}
    mode = {"thickness_corrected": ev.CLEAR_THICKNESS, "plain": ev.CLEAR_PLAIN}[cl.get("mode", "thickness_corrected")]
    wall = mrow[cl["wall"]] if cl.get("wall") else -1
    exempt = mrow[cl["wall_exempt_for"]] if cl.get("wall_exempt_for") else -1
    if mode == ev.CLEAR_THICKNESS and wall < 0 and uniform_half(params):
        # min_m(|d_m| - h) + h == min_m |d_m| for one uniform half thickness h = HALF and no wall: compile the plain rule,
        # which gives the same numbers bit for bit as before B6 (yeast X5 packages, uniform 2 nm)
        mode = ev.CLEAR_PLAIN
    return ev.Tables(M, L, FB, data, n_mem=len(names), clearance=(mode, wall, exempt))


def uniform_half(params):
    """True if every MEM_HALF equals HALF (or the package declares no MEM_HALF)"""
    mh = params.get("MEM_HALF")
    return mh is None or all(float(x) == float(params.get("HALF", 2.0)) for x in np.atleast_1d(mh))


def package_constants(rules, pkg):
    """named scalars a rule table reads from a package file (`constants:` {file, key, names: {NAME: index}}); the
    file must be listed (checksummed) in manifest.files"""
    import json
    c = rules.get("constants")
    if not c:
        return {}
    if pkg is None:
        raise ValueError(f"rule table {rules.get('_path')} reads constants from {c['file']}: the package dir is needed")
    arr = json.load(open(os.path.join(pkg, c["file"])))[c["key"]]
    return {n: float(arr[i]) for n, i in c["names"].items()}


def rule_params(rules, man, pkg=None):
    """the value environment of a package: manifest engine.params + the rule table's package constants"""
    p = dict(man["engine"]["params"])
    p.update(package_constants(rules, pkg))
    return p


def compile_for_manifest(rules, man, pkg=None):
    """compile_tables with the names / codes / params / present items of a package manifest (+ its constants)"""
    label_codes = {L["name"]: L["engine_code"] for L in man["labels"] if L.get("engine_code") is not None}
    mem_codes = {m["name"]: m["engine_code"] for m in man["membranes"]}
    field_slots = {n: f["slot"] for n, f in man["fields"].items()}
    prim_classes = (man.get("prims") or {}).get("classes")
    return compile_tables(rules, rule_params(rules, man, pkg), package_present(man), field_slots, label_codes,
                          mem_codes, prim_classes)


# vocabulary name -> engine code of the sdf_sections_v1 layout (= em_cell1_v7 / analytic_v1 manifests)
DEFAULT_LABEL_CODES = {"exterior": 0, "cytosol": 1, "perinuclear_space": 2, "nucleoplasm": 3, "nucleolus": 4,
                       "mito_ims": 5, "mito_matrix": 6, "er_lumen": 7, "golgi_lumen": 8, "endolysosome_lumen": 9,
                       "nuclear_pore_channel": 10, "outside_imaged_volume": 11, "nuclear_speck": 12,
                       "golgi_intercisternal": 13}
