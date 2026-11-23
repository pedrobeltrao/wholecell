# new in wholecell (freeze v8 E7, 2026-11-19); orientation check ported from jurkat/applications/s0_geom.py::unit_geom
"""Plasma-membrane topology fix and check (freeze v8 criteria §4 E7 / §5d; DECISION 2026-11-19, D7).

Library convention: membrane normal +z = cytoplasm, mid-plane z = 0, so the ectodomain belongs on -z.

Library step `pm_topology_fix` runs after adaptive_trees, so the trees are rotated as well. It reads a TSV of units:
    lib  rule  genes  uniprot  n_PM  evidence
The two rules are:
  glyco_ectodomain_exoplasmic  The unit has a subunit captured by the Jurkat surface glyco-capture, so its ectodomain is
                               extracellular. The unit is rotated 180 deg about x: (x, y, z) -> (x, -y, -z). This is a
                               proper rotation, so chirality is kept. It is applied to `beads` and every `st_c*`
                               array; radii and bead kinds are unchanged.
  tm_alternation               The expected side comes from UniProt TOPO_DOM, with gaps filled by alternating sides
                               across the TRANSMEM helices. The unit is rotated only if the bead balance contradicts
                               that side; otherwise it is recorded as already correct (TRBC1: the constant domain is
                               N-terminal of its single TM, so it is extracellular).
Copies are not affected: no inventory table changes. Each rotated file gains `topology_fixed` (int8 1) and
`topology_rule`. A unit that is listed but missing, or already fixed, is an error.

check(): the E7 gate. It re-runs the S0 orientation rule on every PM unit of a library and compares it with a
reference library.
"""
from __future__ import annotations

import gzip
import json
import os
import re

import numpy as np

HALF = 2.0
ROT = np.array([1.0, -1.0, -1.0], np.float32)          # 180 deg about x
# read-only reference tables under the site data root (config/site.yaml data_root; docs/data_inputs.md)
from ..config import load_site as _load_site  # noqa: E402
VC = str(_load_site().get("data_root") or "{data_root}")
UNIPROT_GLYCO_TOPO = os.path.join(VC, "data", "reference", "uniprot_glyco", "uniprot_human_glyco_topo.tsv.gz")
AFDB_PLDDT = os.path.join(VC, "data", "derived", "afdb_plddt_summary.tsv.gz")
RULES = ("glyco_ectodomain_exoplasmic", "tm_alternation")


def ranges(s, key, note=None):
    out = []
    for m in re.finditer(rf"{key} (\d+)\.\.(\d+);(?: /note=\"([^\"]*)\")?", s or ""):
        if note is None or (m.group(3) and re.search(note, m.group(3))):
            out.append((int(m.group(1)), int(m.group(2)), m.group(3) or ""))
    return out


def rotate(z):
    """{key: array} with beads / st_c* rotated 180 deg about x (exact in float32: sign changes only)"""
    d = {k: z[k] for k in z.files}
    for k in d:
        if k == "beads" or k.startswith("st_c"):
            a = np.asarray(d[k])
            assert a.ndim == 2 and a.shape[1] == 3, (k, a.shape)
            d[k] = (a * ROT.astype(a.dtype)).astype(a.dtype)
    return d


def bead_balance(beads, half=HALF):
    z = np.asarray(beads, float)[:, 2]
    lo, hi = int((z < -half).sum()), int((z > half).sum())
    return lo, hi, (lo / (lo + hi) if lo + hi else np.nan)


def balance_flag(exp_frac, obs_frac):
    """S0 rule: flipped when the expected-larger side and the observed-larger side disagree, both clearly (> 0.15)"""
    if not (np.isfinite(exp_frac) and np.isfinite(obs_frac)):
        return False
    return (exp_frac - 0.5) * (obs_frac - 0.5) < 0 and abs(exp_frac - 0.5) > 0.15 and abs(obs_frac - 0.5) > 0.15


def alternation_counts(topo, tm, length):
    """(extracellular, cytoplasmic) residues: TOPO_DOM where annotated, gaps by alternation across TRANSMEM helices"""
    T = sorted((a, b) for a, b, _ in ranges(tm, "TRANSMEM"))
    if not T:
        return None
    bounds = [1] + [x for a, b in T for x in (a, b + 1)] + [length + 1]
    segs = [(bounds[2 * i], bounds[2 * i + 1] - 1) for i in range(len(T) + 1)]    # soluble segments between TMs
    ann = [(a, b, "e") for a, b, _ in ranges(topo, "TOPO_DOM", "xtracellular|umenal")] + \
          [(a, b, "c") for a, b, _ in ranges(topo, "TOPO_DOM", "ytoplasm")]
    side = [None] * len(segs)
    for i, (s0, s1) in enumerate(segs):
        v = {"e": 0, "c": 0}
        for a, b, k in ann:
            v[k] += max(0, min(b, s1) - max(a, s0) + 1)
        if v["e"] or v["c"]:
            side[i] = "e" if v["e"] >= v["c"] else "c"
    if all(s is None for s in side):
        return None
    k0 = next(i for i, s in enumerate(side) if s is not None)
    flip = {"e": "c", "c": "e"}
    for i in range(len(segs)):                                                   # alternate from the first annotated one
        side[i] = side[k0] if (i - k0) % 2 == 0 else flip[side[k0]]
    ex = sum(max(0, s1 - s0 + 1) for (s0, s1), s in zip(segs, side) if s == "e")
    cy = sum(max(0, s1 - s0 + 1) for (s0, s1), s in zip(segs, side) if s == "c")
    return ex, cy


# ------------------------------------------------------------------ annotations (read-only reference data)
_UP = None
_PL = None


def uniprot():
    global _UP
    if _UP is None:
        import pandas as pd
        u = pd.read_csv(UNIPROT_GLYCO_TOPO, sep="\t")
        u.columns = ["acc", "gene", "length", "topo", "tm", "intramem", "glyco", "compbias", "region", "lipid", "signal",
                     "domain", "repeat", "seq"]
        _UP = u.fillna("").drop_duplicates("gene").set_index("gene")
    return _UP


def plddt(acc):
    global _PL
    if _PL is None:
        _PL = {}
        with gzip.open(AFDB_PLDDT, "rt") as f:
            next(f)
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p[1] == "F1":
                    _PL[p[0]] = p[6]
    s = _PL.get(acc)
    return None if s is None else np.array([float(x) for x in s.split(",")])


def topo_counts(g):
    U = uniprot()
    if g not in U.index:
        return None
    r = U.loc[g]
    ex = sum(b - a + 1 for a, b, _ in ranges(r.topo, "TOPO_DOM", "xtracellular|umenal"))
    cy = sum(b - a + 1 for a, b, _ in ranges(r.topo, "TOPO_DOM", "ytoplasm"))
    return (ex, cy) if ex + cy else None


def orientation(path, genes, stoich, alternation=False):
    """S0 unit_geom orientation part: dict(method, flipped, exp_ecto_frac, obs_lo_frac). alternation=True fills
    TOPO_DOM gaps by TM alternation (balance method only)"""
    z = np.load(path, allow_pickle=True)
    b = z["beads"].astype(float)
    out = dict(method="", flipped=False, exp_ecto_frac=np.nan, obs_lo_frac=np.nan,
               topology_fixed=bool("topology_fixed" in z.files and int(z["topology_fixed"])))
    U = uniprot()
    if len(genes) == 1 and genes[0] in U.index and "plddt_frac_kept" in z.files:
        r = U.loc[genes[0]]
        pl = plddt(r.acc)
        if pl is not None and len(pl) == int(r.length):
            keep = pl >= 50
            seqn = np.arange(1, len(pl) + 1)
            for a, bb, _ in ranges(r.tm, "TRANSMEM"):
                keep |= (seqn >= a) & (seqn <= bb)
            if keep.sum() < 10:
                keep[:] = True
            if keep.sum() == len(b):
                kres = seqn[keep]
                isin = lambda rr, R: np.array([any(a <= x <= bb for a, bb, _ in R) for x in rr], bool)
                em = isin(kres, ranges(r.topo, "TOPO_DOM", "xtracellular"))
                cm = isin(kres, ranges(r.topo, "TOPO_DOM", "ytoplasm"))
                if em.sum() + cm.sum() > 0:
                    ze = b[em, 2].mean() if em.any() else np.nan
                    zc = b[cm, 2].mean() if cm.any() else np.nan
                    out.update(method="residue",
                               flipped=bool((np.isfinite(ze) and ze > 0.5) or (not np.isfinite(ze) and np.isfinite(zc) and zc < -0.5)))
                    return out
    ex = cy = 0
    for g, s in zip(genes, stoich):
        t = None
        if alternation and g in U.index:
            r = U.loc[g]
            t = alternation_counts(r.topo, r.tm, int(r.length))
        if t is None:
            t = topo_counts(g)
        if t:
            ex += s * t[0]
            cy += s * t[1]
    lo, hi, obs = bead_balance(b)
    if ex + cy > 0 and lo + hi > 0:
        exp = ex / (ex + cy)
        out.update(method="balance" + ("+alternation" if alternation else ""), flipped=bool(balance_flag(exp, obs)),
                   exp_ecto_frac=exp, obs_lo_frac=obs)
    else:
        out.update(method="none")
    return out


# ------------------------------------------------------------------ library step
def read_table(path):
    import pandas as pd
    T = pd.read_csv(path, sep="\t", dtype=str, comment="#").fillna("")
    bad = sorted(set(T["rule"]) - set(RULES))
    if bad:
        raise ValueError(f"pm_topology_fix {path}: unknown rules {bad} (known {RULES})")
    if T["lib"].duplicated().any():
        raise ValueError(f"pm_topology_fix {path}: duplicate lib ids")
    return T


def apply(lib_dir, table, log=print):
    """rotate the listed units in lib_dir (writable, in-progress library); returns the per-unit record"""
    rec = []
    for r in table.to_dict("records"):
        p = os.path.join(lib_dir, r["lib"] + ".npz")
        if not os.path.exists(p):
            raise FileNotFoundError(f"pm_topology_fix: {r['lib']} not in the library ({p})")
        z = np.load(p, allow_pickle=True)
        if "topology_fixed" in z.files:
            raise RuntimeError(f"pm_topology_fix: {r['lib']} already fixed")
        lo, hi, obs = bead_balance(z["beads"])
        genes = r["genes"].split(";")
        do = True
        why = r["rule"]
        if r["rule"] == "tm_alternation":
            st = [1] * len(genes)
            o = orientation(p, genes, st, alternation=True)
            do = bool(o["flipped"])
            why = (f"tm_alternation: expected ecto fraction {o['exp_ecto_frac']:.3f}, observed -z fraction {obs:.3f} -> "
                   + ("rotated" if do else "already ectodomain exoplasmic, not rotated"))
        if do:
            d = rotate(z)
            d["topology_fixed"] = np.int8(1)
            d["topology_rule"] = np.str_(r["rule"])
            tmp = p + f".tmp{os.getpid()}.npz"
            np.savez_compressed(tmp, **d)
            os.replace(tmp, p)
        lo2, hi2, obs2 = bead_balance(np.load(p)["beads"])
        rec.append(dict(lib=r["lib"], rule=r["rule"], rotated=do, n_PM=r.get("n_PM", ""), obs_lo_before=obs,
                        obs_lo_after=obs2, note=why))
        log(f"pm_topology_fix: {r['lib']} {'ROTATED' if do else 'unchanged'} (-z bead fraction {obs:.3f} -> {obs2:.3f}; {why})")
    return rec


def step_pm_topology_fix(X, p):
    from .build import _path, _md5, log
    path = _path(X.cfg, X.bdir, p["table"])
    T = read_table(path)
    if X.get("species_spec"):                           # species from the pack table: listed units that are not species
        have = T.lib.map(lambda i: os.path.exists(os.path.join(X.partial, i + ".npz")))   # of this build are skipped
        absent = sorted(T.lib[~have])
        if absent:
            log(f"pm_topology_fix: {len(absent)} listed units are not species of this build, skipped: {absent}")
        T = T[have]
    else:
        absent = []
    rec = apply(X.partial, T, log=log)
    for r in rec:
        if r["rotated"]:
            X.origin[r["lib"] + ".npz"] = "pm_topology_fix"
    X.results["pm_topology_fix"] = dict(table=p["table"], table_md5=_md5(path), units=len(rec),
                                        rotated=[r["lib"] for r in rec if r["rotated"]],
                                        not_rotated=[r["lib"] for r in rec if not r["rotated"]], absent=absent, records=rec)


# ------------------------------------------------------------------ E7 check
def pm_units(pack_table):
    """PM transmembrane units of a pack table: [(lib, genes, stoich, copies, pm_weight)]"""
    import pandas as pd
    P = pd.read_csv(pack_table, sep="\t")
    out = []
    for r in P.itertuples():
        if str(r.transmembrane) != "True":
            continue
        w = json.loads(r.comp_weights).get("plasma_membrane", 0.0) if isinstance(r.comp_weights, str) else 0.0
        if w <= 0:
            continue
        g = str(r.genes).split(";")
        s = [int(float(x)) if x not in ("", "nan") else 1 for x in str(r.stoich).split(";")]
        if len(s) != len(g):
            s = [1] * len(g)
        out.append((r.ing_id, g, s, int(r.copies), float(w)))
    return out


def _md5(p):
    import hashlib
    h = hashlib.md5()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 22), b""):
            h.update(c)
    return h.hexdigest()


def check(lib, ref_lib, pack_table, table):
    """E7: (1) library diff vs the reference library: only the rotated units differ, each = exact rotation of the
    reference; (2) S0 orientation rule on every PM unit: none of the glyco units flagged, the other reference flags
    unchanged (listed as known residuals), no new flag; tm_alternation units judged with alternation."""
    T = read_table(table) if isinstance(table, str) else table
    fix = dict(zip(T["lib"], T["rule"]))
    glyco = [l for l, r in fix.items() if r == "glyco_ectodomain_exoplasmic"]
    alt = [l for l, r in fix.items() if r == "tm_alternation"]
    files = sorted(f for f in os.listdir(lib) if f.endswith(".npz"))
    rfiles = sorted(f for f in os.listdir(ref_lib) if f.endswith(".npz"))
    changed, bad_rot = [], []
    for f in sorted(set(files) & set(rfiles)):
        a, b = os.path.join(lib, f), os.path.join(ref_lib, f)
        if _md5(a) == _md5(b):
            continue
        z, zr = np.load(a, allow_pickle=True), np.load(b, allow_pickle=True)
        lid = f[:-4]
        if "topology_fixed" in z.files:
            exp = rotate(zr)
            ok = all(np.array_equal(z[k], exp[k]) for k in exp) and set(z.files) - set(zr.files) == {"topology_fixed", "topology_rule"}
            changed.append(lid)
            if not ok:
                bad_rot.append(lid)
        else:
            same = set(z.files) == set(zr.files) and all(np.array_equal(z[k], zr[k]) for k in z.files)
            if not same:
                changed.append(lid)
                bad_rot.append(lid + " (content differs, not a topology fix)")
    rows = []
    for lid, g, s, cp, w in pm_units(pack_table):
        p, pr = os.path.join(lib, lid + ".npz"), os.path.join(ref_lib, lid + ".npz")
        if not os.path.exists(p):
            continue
        o = orientation(p, g, s)
        orr = orientation(pr, g, s) if os.path.exists(pr) else dict(flipped=None)
        oa = orientation(p, g, s, alternation=True) if lid in alt else None
        rows.append(dict(lib=lid, genes=";".join(g), copies=cp, pm_weight=w, method=o["method"], flagged=o["flipped"],
                         ref_flagged=orr["flipped"], flagged_alternation=None if oa is None else oa["flipped"],
                         fixed=o["topology_fixed"], rule=fix.get(lid, ""), exp_ecto_frac=o["exp_ecto_frac"],
                         obs_lo_frac=o["obs_lo_frac"]))
    R = {r["lib"]: r for r in rows}
    ref_flags = sorted(l for l, r in R.items() if r["ref_flagged"])
    new_flags = sorted(l for l, r in R.items() if r["flagged"])
    glyco_flagged = sorted(l for l in glyco if R.get(l, {}).get("flagged"))
    glyco_missing = sorted(l for l in glyco if l not in R)
    alt_flagged = sorted(l for l in alt if R.get(l, {}).get("flagged_alternation"))
    expect_rotated = sorted(set(glyco) | {l for l in alt if R.get(l, {}).get("fixed")})
    residual = sorted(set(ref_flags) - set(glyco))
    out = dict(
        pm_units=len(rows), ref_flagged=len(ref_flags), flagged=len(new_flags),
        changed_units=sorted(changed), rotation_errors=bad_rot,
        only_fix_units_changed=sorted(changed) == expect_rotated and not bad_rot,
        glyco_units=len(glyco), glyco_flagged=glyco_flagged, glyco_missing=glyco_missing,
        glyco_copies=int(sum(R[l]["copies"] for l in glyco if l in R)),
        alternation_units={l: dict(flagged_topodom=R.get(l, {}).get("flagged"), flagged_alternation=R.get(l, {}).get("flagged_alternation"),
                                   rotated=R.get(l, {}).get("fixed")) for l in alt},
        alternation_flagged=alt_flagged,
        residual_flags=residual, residual_copies=int(sum(R[l]["copies"] for l in residual)),
        new_flags=sorted(set(new_flags) - set(ref_flags)),
        residual_unchanged=sorted(set(new_flags)) == residual,
    )
    out["ok"] = bool(out["only_fix_units_changed"] and not glyco_flagged and not glyco_missing and not alt_flagged
                     and not out["new_flags"] and out["residual_unchanged"])
    return out, rows
