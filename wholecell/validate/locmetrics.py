"""One definition of the localisation agreement metrics (unification plan s4, check C1), used by both cells.

Per gene: a model distribution p over the cell's compartment classes (intended or placed copies, normalised to 1) and a
reference annotation q (class -> weight; weight >= 0.999 = a main class; a reference given as a plain set of classes
has every class at weight 1, i.e. every class is main).

  supported   classes the reference annotates (plus the classes a special term supports, Jurkat vocabulary only)
  unsup       model share in classes the reference does not annotate
  miss        share of the reference's main classes holding < 5 % of the model's copies
  agree       1 - (0.5 unsup + 0.5 miss)
  agree_tvd   1 - total-variation distance between p and the normalised reference weights
  tophit      1 if the model's largest class is supported (ties: first class in the cell's class order)

Aggregates over a gene set: plain mean, and mass-weighted mean (weight = model protein mass of the gene, copies x
mass), over the genes where the metric is defined.

Origins: Jurkat loc_compare.score (legacy jurkat/validate/loc_compare.py; numerics unchanged, it now calls score()
here) and yeast validation/localization.py:score @ e3b7448 ('agree' there is the same formula with every reference
class main; 'top1' is tophit). The yeast reference sets carry no weights, so for yeast agree_tvd compares p with the
uniform distribution over the reference classes.
"""
from __future__ import annotations

import numpy as np

MISS_SHARE = 0.05
MAIN_W = 0.999


def score(p, q, comps, special=None, main_needed=True):
    """p: {class: share}, q: {class: weight} (or a set of classes), comps: class order -> metric dict"""
    if not isinstance(q, dict):
        q = {c: 1.0 for c in q}
    sup = {c for c in q if c in comps}
    for sp, cs in (special or {}).items():
        if sp in q:
            sup |= set(cs)
    unsup = sum(v for c, v in p.items() if c not in sup)
    main = [c for c, w in q.items() if w >= MAIN_W and c in comps]
    miss = (sum(p.get(c, 0.0) < MISS_SHARE for c in main) / len(main)) if (main and main_needed) else 0.0
    qq = {c: w for c, w in q.items() if c in comps}
    s = sum(qq.values())
    tvd = 0.5 * sum(abs(p.get(c, 0.0) - qq.get(c, 0.0) / s) for c in comps) if s else np.nan
    top = max(p, key=p.get) if p else None
    return dict(unsup=unsup, miss=miss, score=0.5 * unsup + 0.5 * miss, tvd=tvd,
                tophit=float(top in sup) if top else np.nan)


def gene_row(p, q, comps, special=None):
    sc = score(p, q, comps, special)
    return dict(agree=1 - sc["score"], agree_tvd=1 - sc["tvd"], tophit=sc["tophit"], supported_mass=1 - sc["unsup"],
                main_present=1 - sc["miss"])


def aggregate(values, mass):
    """(plain mean, mass-weighted mean) over the genes where the value is defined; unrounded"""
    v = np.asarray(values, float); w = np.asarray(mass, float)
    m = np.isfinite(v)
    if not m.any():
        return None, None
    return float(v[m].mean()), float((v[m] * w[m]).sum() / w[m].sum()) if w[m].sum() > 0 else None


def score_set(model, ref, comps, special=None):
    """model: {gene: (p dict in class order, mass)}, ref: {gene: q} -> (summary dict, per-gene rows). Genes in both with
    a non-empty model distribution (sum > 0)."""
    rows = []
    for g, q in ref.items():
        if g not in model:
            continue
        p, mass = model[g]
        if not p or sum(p.values()) <= 0:
            continue
        rows.append(dict(gene=g, mass=mass, **gene_row(p, q, comps, special)))
    out = dict(n_genes=len(rows))
    w = [r["mass"] for r in rows]
    for k in ("agree", "agree_tvd", "tophit"):
        a, b = aggregate([r[k] for r in rows], w) if rows else (None, None)
        out[k] = dict(unweighted=a, mass_weighted=b)
    return out, rows
