"""Pure-python reference of cellgeom_eval_v1 (no numba): membrane SDFs and labels of single points, interpreting the
rule YAML directly (not the compiled G tables). For tests only (slow: ~10^4 points); it cross-checks the numba kernels
and the table compiler. Spec: docs/reference/cellgeom_v1.md §4.3; kernels mirrored from engines/cellgeom_eval_v1.py.

    ref = RefCell(cell, rules)          # cell = cellgeom.load(...), rules = rules.load_rules(path)
    out, matrix = ref.mem_sdfs(p);  code = ref.label_of(p)
"""
from __future__ import annotations

import math

import numpy as np

from wholecell.geometry import rules as R
from wholecell.geometry.engines import cellgeom_eval_v1 as ev


def sample(h10, idx, data, p):
    """sparse-brick field (header10, index, data) at p: trilinear in a stored brick, clamps elsewhere"""
    h = h10[3]; B = int(h10[7])
    f = [(p[a] - h10[a]) / h for a in range(3)]
    nb = [int(h10[4]), int(h10[5]), int(h10[6])]
    if min(f) < 0.0 or any(f[a] >= nb[a] * B for a in range(3)):
        return h10[9]
    i = [int(x) for x in f]
    b = [x // B for x in i]
    bi = int(idx[(b[0] * nb[1] + b[1]) * nb[2] + b[2]])
    if bi == -1:
        return -h10[8]
    if bi == -2:
        return h10[9]
    n1 = B + 1
    t = [f[a] - i[a] for a in range(3)]
    q = bi * n1 ** 3 + ((i[0] - b[0] * B) * n1 + (i[1] - b[1] * B)) * n1 + (i[2] - b[2] * B)

    def v(dx, dy, dz):
        return data[q + dx * n1 * n1 + dy * n1 + dz]
    c00 = v(0, 0, 0) + (v(1, 0, 0) - v(0, 0, 0)) * t[0]
    c01 = v(0, 0, 1) + (v(1, 0, 1) - v(0, 0, 1)) * t[0]
    c10 = v(0, 1, 0) + (v(1, 1, 0) - v(0, 1, 0)) * t[0]
    c11 = v(0, 1, 1) + (v(1, 1, 1) - v(0, 1, 1)) * t[0]
    c0 = c00 + (c10 - c00) * t[1]; c1 = c01 + (c11 - c01) * t[1]
    return c0 + (c1 - c0) * t[2]


# ----------------------------------------------------------------- analytic primitives (cell3d formulas, plain python)
def prim_sdf(p, t, q):
    px, py, pz = p
    if t == ev.P_SPHERE:
        return math.sqrt((px - q[0]) ** 2 + (py - q[1]) ** 2 + (pz - q[2]) ** 2) - q[3]
    if t == ev.P_CAPSULE:
        a = p - q[0:3]; b = q[3:6] - q[0:3]
        h = min(max(float(a @ b) / float(b @ b), 0.0), 1.0)
        d = a - b * h
        return math.sqrt(float(d @ d)) - q[6]
    if t == ev.P_PILL:
        v = p - q[0:3]; z = float(v @ q[3:6]); r = v - z * q[3:6]
        s = max(math.sqrt(float(r @ r)) - q[6], 0.0)
        return math.sqrt(s * s + z * z) - q[7]
    if t == ev.P_TORUS:
        v = p - q[0:3]; h = float(v @ q[3:6])
        rho = math.sqrt(max(float(v @ v) - h * h, 0.0))
        return math.sqrt((rho - q[6]) ** 2 + h * h) - q[7]
    if t == ev.P_SLOT:
        m = 0.5 * (q[0:3] + q[3:6]); u = m / math.sqrt(float(m @ m))
        e = q[3:6] - q[0:3]; L = math.sqrt(float(e @ e)); e = e / L
        v = p - q[0:3]
        s = min(max(float(v @ e), 0.0), L); h = min(max(float(v @ u), 0.0), q[7])
        d = v - s * e - h * u
        return math.sqrt(float(d @ d)) - q[6]
    r = math.sqrt(float(p @ p))         # P_CERPATCH (plain shell patch)
    c = min(max(float(p @ q[0:3]) / max(r, 1e-9), -1.0), 1.0)
    s = max(math.acos(c) - q[4], 0.0) * r
    return math.sqrt(s * s + (r - q[5]) ** 2) - q[6]


def cer_sdf(p, q, pm, Rp):
    """cortical-ER sheet at depth Rp - q[5] below the reference membrane value pm (cell3d sdf_cer)"""
    r = math.sqrt(float(p @ p))
    c = min(1.0, max(-1.0, float(p @ q[0:3]) / max(r, 1e-9)))
    s = max(math.acos(c) - q[4], 0.0) * r
    dr = (Rp - q[5]) - pm
    return math.sqrt(s * s + dr * dr) - q[6]


def smin(a, b, k):
    h = max(k - abs(a - b), 0.0) / k
    return min(a, b) - h * h * k * 0.25


def analytic(spec, params, p):
    """field fallback shape (sphere, ellipsoid_approx) at p"""
    v = lambda x: R.value(x, params)  # noqa: E731
    c = np.array([v(x) for x in spec.get("centre", [0, 0, 0])])
    d = p - c
    if spec["kind"] == "sphere":
        return math.sqrt(float(d @ d)) - v(spec["radius"])
    ax = [v(x) for x in spec["axes"]]; rot = np.array([v(x) for x in spec["rot"]]).reshape(3, 3)
    l = rot @ d / np.array(ax)
    k0 = math.sqrt(float(l @ l)); k1 = math.sqrt(float(((l / np.array(ax)) ** 2).sum()))
    if k1 < 1e-12:
        return -min(ax)
    return k0 * (k0 - 1.0) / k1


class RefCell:
    def __init__(self, cell, rules):
        G = cell.G
        self.parts = ev.split_legacy_G(G, g_n=ev.G_N)
        self.params = ev.params_from_G(G)
        self.params.update(getattr(cell, "rule_params", None) or {})
        self.pores = np.asarray(cell.pores, np.float64)
        self.rules = rules
        self.PT, self.PQ, self.PC = cell.PT, cell.PQ, cell.PC
        man = getattr(cell, "manifest", None)
        if man is not None:
            self.slots = {n: f["slot"] for n, f in man["fields"].items()}
            self.codes = {L["name"]: L["engine_code"] for L in man["labels"] if L.get("engine_code") is not None}
            present = R.package_present(man)
            self.rows = R.active_rows(rules, present)
            self.classes = {n: k for k, n in enumerate((man.get("prims") or {}).get("classes", []))}
        else:
            self.slots = {n: k for k, n in enumerate(ev.FIELD_NAMES)}
            self.codes = R.DEFAULT_LABEL_CODES
            self.rows = rules["labels"]
            present = {"field:" + n for n in ev.FIELD_NAMES}
            self.classes = {}
        self.present = present
        self.fallback = {n: f["fallback"] for n, f in (rules.get("fields") or {}).items()
                         if f.get("fallback") and "field:" + n not in present}
        self.mem_index = {n: i for i, n in enumerate(R.membrane_names(rules))}

    # ------------------------------------------------------------- primitives
    def field(self, name, p):
        if name in self.fallback:
            return analytic(self.fallback[name], self.params, p)
        k = self.slots.get(name)
        if k is None or k >= len(self.parts["fields"]):
            return 1e9
        return sample(*self.parts["fields"][k], p)

    def prims(self, cls, p, ref=None, Rp=None):
        c = self.classes[cls]
        d = 1e9
        for i in np.flatnonzero(self.PC == c):
            t, q = int(self.PT[i]), self.PQ[i]
            d = min(d, cer_sdf(p, q, ref, Rp) if (t == ev.P_CERPATCH and ref is not None) else prim_sdf(p, t, q))
        return d

    def face(self, p):
        return p[2] - self.parts["face"] if self.parts["face"] is not None else 1e9

    def pore_index(self, p, rc, zwin):
        for k, q in enumerate(self.pores):
            v = p - q[:3]
            z = float(v @ q[3:])
            if z < -zwin or z > zwin:
                continue
            r = v - z * q[3:]
            if float(r @ r) < rc * rc:
                return k
        return -1

    def in_speckle(self, p):
        s = self.parts["speckle"]
        if s is None:
            return False
        for k in range(int(s[0])):
            c = s[1 + 4 * k:5 + 4 * k]
            if float((p - c[:3]) @ (p - c[:3])) < c[3] * c[3]:
                return True
        return False

    def crista_lumen(self, p, b):
        s = self.parts["cristae"]; P = self.params
        n = int(s[0]); h = s[4]; dims = [int(s[5]), int(s[6]), int(s[7])]
        ii = [int(math.floor((p[a] - s[1 + a]) / h)) for a in range(3)]
        if min(ii) < 0 or any(ii[a] >= dims[a] for a in range(3)):
            return 1e9
        k = int(s[8 + 12 * n + (ii[0] * dims[1] + ii[1]) * dims[2] + ii[2]])
        if k < 0:
            return 1e9
        q = s[8 + 12 * k:8 + 12 * k + 12]
        if q[11] == 0.0:
            return 1e9
        S = q[4]
        t = q[0] * p[0] + q[1] * p[1] + q[2] * p[2] - q[3]
        u = t - S * math.floor(t / S + 0.5)
        slab = abs(u) - P["CR_W"]
        lam = max(slab, b + P["CR_GAP"])
        if P.get("CR_RIM", 0.0) > 0.0:          # B6b rounded rim
            r = min(P["CR_RIM"], P["CR_W"]); e_ = b + P["CR_GAP"]
            lam = min(-r, lam) + math.sqrt(max(slab + r, 0.0) ** 2 + max(e_ + r, 0.0) ** 2)
        w = abs(q[5] * (p[0] - q[8]) + q[6] * (p[1] - q[9]) + q[7] * (p[2] - q[10])) - P["CR_JW"]
        return min(lam, max(slab, w))

    # ------------------------------------------------------------- membrane table
    def mem_sdfs(self, p):
        p = np.asarray(p, np.float64)
        out = np.zeros(len(self.mem_index)); matrix = 0.0
        val = lambda v: R.value(v, self.params)  # noqa: E731
        i = -1
        for m in self.rules["membranes"]:
            k = m["kind"]
            if k == "smooth_union":
                j = self.mem_index[m["into"]]; s = float(m.get("side", 1))
                d = self.prims(m["prims"], p); x = s * out[j]; kk = val(m["k"])
                if d < x + kk:
                    out[j] = s * smin(x, d, kk)
                continue
            i += 1
            if k == "field":
                off = m.get("offset", 0)
                off = val(off[0]) - val(off[1]) if isinstance(off, list) else val(off)
                out[i] = float(m.get("sign", 1)) * self.field(m["field"], p) - off
            elif k == "wall":
                out[i] = min(out[self.mem_index[m["of"]]], self.face(p))
            elif k == "shell_pores":
                pw, pr, zw = val(m["width"]), val(m["pore_r"]), val(m["zwin"])
                hp = self.field(m["field"], p) + val(m.get("field_offset", 0))
                d = abs(hp) - 0.5 * pw
                rc = pr + 0.5 * pw
                if abs(hp) < rc + pw:
                    j = self.pore_index(p, rc, zw)
                    if j >= 0:
                        q = self.pores[j]; v = p - q[:3]; z = float(v @ q[3:])
                        rho = math.sqrt(max(float(v @ v) - z * z, 0.0))
                        d = math.sqrt((rho - rc) ** 2 + hp * hp) - 0.5 * pw
                out[i] = d
            elif k == "inner_cristae":
                matrix = out[self.mem_index[m["of"]]] + val(m["offset"])
                if self.parts["cristae"] is not None and matrix < self.params["CR_W"]:
                    matrix = max(matrix, -self.crista_lumen(p, matrix))
                out[i] = -matrix
            elif k == "golgi_cisternae":
                g = self.field(m["field"], p)
                other = min(out[self.mem_index[c]] for c in m["clear_of"])
                gc = self.parts["golgi_cis"]
                if gc is not None and g < 0.0 and other > self.params["CIS_CLR"]:
                    W, S, parts = gc
                    v = -sample(*parts, p) - W
                    w = v - S * math.floor(v / S + 0.5)
                    g = max(g, abs(w) - W)
                out[i] = g
            elif k == "prims":
                base = m.get("base_field")
                v = self.field(base, p) if base and "field:" + base in self.present else 1e9
                if m.get("ref"):
                    out[i] = min(v, self.prims(m["prims"], p, out[self.mem_index[m["ref"]]], val(m["ref_radius"])))
                else:
                    out[i] = min(v, self.prims(m["prims"], p))
            elif k == "prims_difference":
                matrix = max(self.prims(m["prims"], p), -self.prims(m["minus"], p))
                out[i] = float(m.get("sign", 1)) * matrix
            else:
                raise ValueError(k)
        return out, matrix

    # ------------------------------------------------------------- label rule
    def term(self, t, p, out):
        if t in self.mem_index:
            return out[self.mem_index[t]]
        if t.startswith("field:"):
            return self.field(t[6:], p)
        if t == "face":
            return self.face(p)
        if t == "pore":
            m = next(m for m in self.rules["membranes"] if m["kind"] == "shell_pores")
            rc = R.value(m["pore_r"], self.params) + 0.5 * R.value(m["width"], self.params)
            return -1.0 if self.pore_index(p, rc, R.value(m["zwin"], self.params)) >= 0 else 1.0
        if t == "speckle":
            return -1.0 if self.in_speckle(p) else 1.0
        if t.startswith("shape:"):
            sh = self.rules["shapes"][t[6:]]
            c = np.array([R.value(x, self.params) for x in sh["centre"]]); r = R.value(sh["radius"], self.params)
            return -1.0 if float((p - c) @ (p - c)) < r * r else 1.0
        raise ValueError(t)

    def label_of(self, p):
        p = np.asarray(p, np.float64)
        out, _ = self.mem_sdfs(p)
        for row in self.rows:
            ok = True
            for s in row.get("if", []):
                t, op, v = R.parse_condition(s)
                x = self.term(t, p, out); v = R.value(v, self.params)
                if not (x < v if op == "<" else x > v):
                    ok = False
                    break
            if ok:
                return self.codes[row["label"]]
        return -1
