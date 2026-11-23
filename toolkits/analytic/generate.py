"""Analytic cellgeom_v1 generator (engine sdf_sections_v1): bricked SDF fields computed from analytic shapes.

    PYTHONHASHSEED=0 python -m toolkits.analytic.generate [toolkits/analytic/specs/testcell.yaml] [--dest DIR]

Spec (YAML, see specs/testcell.yaml): cell (ellipsoid), nucleus (ellipsoid = perinuclear midsurface ne_mid) with pores,
optional nucleolus (sphere), mito (union of capsules) with optional procedural cristae, ER (union of tubes = capsules
and sheets = flat discs). Golgi and endo fields are absent (empty fields: outside everywhere). Grid h, brick B.
Deterministic: pores and cristae draw from wholecell.seeds streams 'geometry.analytic.pores' / '.cristae' of the spec
seed; everything else is analytic. Output: {geom_root}/<cell_id>/<geometry_id>/ (validated, immutable).

Ellipsoid SDF: the standard first-order approximation k0 (k0 - 1) / k1 (k0 = |p/r|, k1 = |p/r^2|), exact on the
surface, |grad| ~ 1 near it (checked by the validator's eikonal test).
"""
import argparse
import os
import sys

import numpy as np
import yaml
from skimage.measure import marching_cubes

from wholecell import seeds
from wholecell.build import git_info
from wholecell.config import REPO, load_site
from wholecell.geometry import cellgeom
from wholecell.geometry.engines import sdf_sections_v1 as e

DEFAULT_SPEC = os.path.join(REPO, "toolkits", "analytic", "specs", "testcell.yaml")


# ----------------------------------------------------------------- analytic SDFs (vectorised over P (n, 3))
def sdf_ellipsoid(P, centre, semi_axes):
    q = (P - np.asarray(centre, float)) / np.asarray(semi_axes, float)
    k0 = np.linalg.norm(q, axis=1)
    k1 = np.linalg.norm(q / np.asarray(semi_axes, float), axis=1)
    return np.where(k1 > 0, k0 * (k0 - 1.0) / np.maximum(k1, 1e-30), -min(semi_axes))


def grad_ellipsoid(P, centre, semi_axes, eps=0.5):
    g = np.zeros_like(P)
    for ax in range(3):
        d = np.zeros(3); d[ax] = eps
        g[:, ax] = (sdf_ellipsoid(P + d, centre, semi_axes) - sdf_ellipsoid(P - d, centre, semi_axes)) / (2 * eps)
    return g


def sdf_sphere(P, centre, radius):
    return np.linalg.norm(P - np.asarray(centre, float), axis=1) - radius


def sdf_capsule(P, a, b, r):
    a = np.asarray(a, float); b = np.asarray(b, float); ab = b - a
    t = np.clip(((P - a) @ ab) / (ab @ ab), 0.0, 1.0)
    return np.linalg.norm(P - (a + t[:, None] * ab), axis=1) - r


def sdf_disc(P, centre, normal, radius, half_thickness):
    n = np.asarray(normal, float); n /= np.linalg.norm(n)
    v = P - np.asarray(centre, float)
    h = v @ n
    rho = np.linalg.norm(v - h[:, None] * n, axis=1)
    dx = rho - radius; dz = np.abs(h) - half_thickness
    out = np.sqrt(np.maximum(dx, 0) ** 2 + np.maximum(dz, 0) ** 2)
    return out + np.minimum(np.maximum(dx, dz), 0.0)


def field_fn(name, spec):
    s = spec[name]
    if name in ("cell", "nucleus"):
        el = s["ellipsoid"]
        return lambda P: sdf_ellipsoid(P, el["centre"], el["semi_axes"])
    if name == "nucleolus":
        sp = s["sphere"]
        return lambda P: sdf_sphere(P, sp["centre"], sp["radius"])
    if name == "mito":
        return lambda P: np.min([sdf_capsule(P, c["a"], c["b"], c["r"]) for c in s["capsules"]], axis=0)
    if name == "er":
        def f(P):
            parts = [sdf_capsule(P, t["a"], t["b"], t["r"]) for t in s.get("tubes", [])]
            parts += [sdf_disc(P, d["centre"], d["normal"], d["radius"], d["half_thickness"]) for d in s.get("sheets", [])]
            return np.min(parts, axis=0)
        return f
    raise KeyError(name)


# ----------------------------------------------------------------- bricks (jgeom_build.bricks encoding)
def bricks(sdf, lo, h, B, band):
    """sparse brick encoding, numerically as jgeom_build.bricks / cellshape.bricks: (header10, index, data)"""
    bin_, bout = band
    n = np.array(sdf.shape); nb_ = -(-(n - 1) // B); pad = nb_ * B + 1 - n
    S = np.clip(np.pad(sdf, [(0, p) for p in pad], mode="edge"), -bin_, bout)
    idx = np.empty(tuple(nb_), np.int64); data = []
    for bx in range(nb_[0]):
        for by in range(nb_[1]):
            for bz in range(nb_[2]):
                blk = S[bx * B:bx * B + B + 1, by * B:by * B + B + 1, bz * B:bz * B + B + 1]
                if blk.max() <= -bin_: idx[bx, by, bz] = -1
                elif blk.min() >= bout: idx[bx, by, bz] = -2
                else: idx[bx, by, bz] = len(data); data.append(blk.ravel())
    head = np.array([lo[0], lo[1], lo[2], h, nb_[0], nb_[1], nb_[2], B, bin_, bout], np.float64)
    dat = np.concatenate(data) if data else np.zeros(0)
    return head, idx.ravel().astype(np.float64), dat.astype(np.float64)


def mc_area_um2(sdf, h):
    if not (sdf.min() < 0 < sdf.max()):
        return 0.0
    v, f, _, _ = marching_cubes(sdf, 0.0, spacing=(h, h, h))
    return float(np.linalg.norm(np.cross(v[f[:, 1]] - v[f[:, 0]], v[f[:, 2]] - v[f[:, 0]]), axis=1).sum() / 2 / 1e6)


# ----------------------------------------------------------------- procedural parts
def make_pores(spec, rng):
    nu = spec["nucleus"]; el = nu["ellipsoid"]; c = np.asarray(el["centre"], float); r = np.asarray(el["semi_axes"], float)
    n, sep = int(nu["pores"]["n"]), float(nu["pores"]["min_sep_nm"])
    out = []
    for _ in range(200 * n):
        u = rng.normal(size=3); u /= np.linalg.norm(u)
        p = (c + u * r)[None, :]
        for _ in range(20):                                      # Newton projection onto the zero level
            d = sdf_ellipsoid(p, c, r); g = grad_ellipsoid(p, c, r)
            p = p - (d / np.maximum((g * g).sum(1), 1e-12))[:, None] * g
        p = p[0]
        if all(np.linalg.norm(p - q[:3]) >= sep for q in out):
            g = grad_ellipsoid(p[None, :], c, r)[0]
            out.append(np.r_[p, g / np.linalg.norm(g)])
            if len(out) == n:
                break
    return np.array(out)


def make_cristae(spec, rng):
    """cristae_lamellae_v1 section for the capsule mitochondria (layout: sdf_sections_v1.crista_lumen)"""
    m = spec["mito"]; cs = m["cristae"]; caps = m["capsules"]
    idh, reach = float(cs["id_h_nm"]), float(cs["id_reach_nm"])
    par = []
    for c in caps:
        a = np.asarray(c["a"], float); b = np.asarray(c["b"], float); ax = (b - a) / np.linalg.norm(b - a)
        v = np.cross(ax, [0.0, 0.0, 1.0])
        if np.linalg.norm(v) < 1e-6:
            v = np.cross(ax, [1.0, 0.0, 0.0])
        v /= np.linalg.norm(v)
        S = rng.uniform(*cs["S_range"])
        par.append([*ax, rng.uniform(0, S), S, *v, *((a + b) / 2), 1.0])
    par = np.array(par)
    lo = np.min([np.minimum(c["a"], c["b"]) for c in caps], axis=0) - max(c["r"] for c in caps) - reach - idh
    hi = np.max([np.maximum(c["a"], c["b"]) for c in caps], axis=0) + max(c["r"] for c in caps) + reach + idh
    dims = np.ceil((hi - lo) / idh).astype(int)
    g = [lo[i] + (np.arange(dims[i]) + 0.5) * idh for i in range(3)]
    X = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3)
    D = np.stack([sdf_capsule(X, c["a"], c["b"], c["r"]) for c in caps], 1)
    ids = np.argmin(D, 1).astype(np.float64); ids[D.min(1) > reach] = -1.0
    return np.r_[len(caps), lo, idh, dims, par.ravel(), ids], par


def label_rule(fields, sections, clip):
    """declarative mirror of sdf_sections_v1.label_from restricted to the present fields / sections"""
    r = []
    if clip:
        r.append({"if": [{"clip": "open_face"}], "label": "outside_imaged_volume"})
    r += [{"if": [{"mem": "CELL_OUTER", "lt": 0}], "label": "exterior"}]
    if "ne_mid" in fields:
        nuc = [{"if": [{"field": "ne_mid", "gt": "-PERI_W-PORE_LW"}, {"in_pore": "PORE_R+0.5*PERI_W", "zwin": "PORE_ZWIN"}],
                "label": "nuclear_pore_channel"}]
        if "nucleolus" in fields:
            nuc.append({"if": [{"field": "nucleolus", "lt": 0}], "label": "nucleolus"})
        if "speckle" in sections:
            nuc.append({"if": [{"in_section": "speckle"}], "label": "nuclear_speck"})
        nuc.append({"label": "nucleoplasm"})
        r += [{"if": [{"mem": "NE", "lt": 0}], "label": "perinuclear_space"},
              {"if": [{"field": "ne_mid", "lt": "-0.5*PERI_W"}], "then": nuc},
              {"if": [{"field": "ne_mid", "lt": 0}], "label": "nuclear_pore_channel"}]
    if "mito" in fields:
        r.append({"if": [{"mem": "MITO_OM", "lt": 0}], "then": [{"if": [{"var": "matrix", "lt": 0}], "label": "mito_matrix"},
                                                                 {"label": "mito_ims"}]})
    for f, mem, lab in (("er", "ER", "er_lumen"), ("golgi", "GOLGI", "golgi_lumen"), ("endo", "ENDO", "endolysosome_lumen")):
        if f in fields:
            r.append({"if": [{"mem": mem, "lt": 0}], "label": lab})
    if "golgi_cis" in sections:
        r.append({"if": [{"has_section": "golgi_cis"}, {"field": "golgi", "lt": 0}], "label": "golgi_intercisternal"})
    r.append({"label": "cytosol"})
    return r


def _rule_labels(rules):
    for x in rules:
        if "label" in x:
            yield x["label"]
        else:
            yield from _rule_labels(x["then"])


def generate(spec, dest=None, validate_kw=None):
    from toolkits.em.convert_legacy import ENGINE_LABELS, labels_table, membranes_table
    site = load_site()
    cfg_seed = {"seeds": {"master": int(spec["seed"])}}
    h = float(spec["grid_h_nm"]); B = int(spec["brick"])
    dest = dest or os.path.join(site["geom_root"], spec["cell_id"], spec["geometry_id"])
    params = dict(e.DEFAULT_PARAMS); params.update(spec.get("engine_params") or {})
    el = spec["cell"]["ellipsoid"]; ca = np.asarray(el["semi_axes"], float); cc = np.asarray(el["centre"], float)
    pad = 160.0 + 4 * h
    lo = np.floor((cc - ca - pad) / h) * h; hi = np.ceil((cc + ca + pad) / h) * h
    n = np.round((hi - lo) / h).astype(int) + 1
    ax = [lo[i] + np.arange(n[i]) * h for i in range(3)]
    X = np.stack(np.meshgrid(*ax, indexing="ij"), -1).reshape(-1, 3)
    names = {"cell": "cell", "ne_mid": "nucleus", "mito": "mito", "er": "er", "nucleolus": "nucleolus"}
    grids, summ = {}, {}
    for fname, key in names.items():
        if key not in spec:
            continue
        grids[fname] = field_fn(key, spec)(X).reshape(tuple(n))
    # separations between surfaces (on the grid nodes): organelles inside the cell, outside the NE, apart
    sep = float(spec.get("min_separation_nm", 40.0)); hw = 0.5 * params["PERI_W"]
    chk = {}
    for f in ("mito", "er"):
        if f not in grids:
            continue
        inside = grids[f] < 0
        chk[f"{f}_pm_clearance_nm"] = float(-grids["cell"][inside].max())
        if "ne_mid" in grids:
            chk[f"{f}_ne_clearance_nm"] = float(grids["ne_mid"][inside].min() - hw)
    if "mito" in grids and "er" in grids:
        chk["mito_er_clearance_nm"] = float((grids["er"][grids["mito"] < 0]).min())
    if "nucleolus" in grids:
        chk["nucleolus_ne_clearance_nm"] = float(-(grids["ne_mid"][grids["nucleolus"] < 0]).max() - hw)
    bad = {k: v for k, v in chk.items() if v < sep}
    if bad:
        raise ValueError(f"spec violates min_separation_nm={sep}: {bad}")
    w = cellgeom.PackageWriter(dest)
    bands = {"cell": spec["cell"]["band_nm"], "ne_mid": spec["nucleus"]["band_nm"]}
    for f, key in names.items():
        if f in grids and f not in bands:
            bands[f] = spec[key]["band_nm"]
    for k, fname in enumerate(e.FIELD_NAMES):
        if fname not in grids:
            w.add_empty_field(fname, k)
            continue
        hd, ix, dt = bricks(grids[fname], lo, h, B, bands[fname])
        w.add_field(fname, k, hd, ix, dt, band_nm=list(map(float, bands[fname])), source="analytic",
                    shape=spec[names[fname]])
        g = grids[fname]
        if fname == "ne_mid":
            summ[fname] = {"vol_inside_um3": round(float((g < 0).sum() * h ** 3 / 1e9), 4), "area_um2": round(mc_area_um2(g, h), 3)}
        else:
            summ[fname] = {"vol_sdf_um3": round(float((g < 0).sum() * h ** 3 / 1e9), 4), "area_um2": round(mc_area_um2(g, h), 3)}
    for fname in ("golgi", "endo"):
        summ[fname] = {"vol_sdf_um3": 0.0, "area_um2": 0.0}
    summ["open_face"] = {"f_vol": 1.0, "f_pm": 1.0, "note": "no open face"}
    summ["checks"] = chk
    sections = {}
    if spec.get("mito", {}).get("cristae", {}).get("enabled"):
        sec, par = make_cristae(spec, seeds.rng(cfg_seed, "geometry.analytic.cristae"))
        w.add_section("cristae", "cristae_lamellae_v1", {"section": sec},
                      params={"n_mito": int(len(par)), "S_nm": [round(float(x), 2) for x in par[:, 4]],
                              "S_range": spec["mito"]["cristae"]["S_range"]},
                      description="lamellar cristae perpendicular to each capsule axis (sdf_sections_v1.crista_lumen)",
                      engine_params=["CR_W", "CR_GAP", "CR_JW", "IMS_W"])
        sections["cristae"] = 1
        summ["cristae"] = {"n_mito": int(len(par)), "with_cristae": int(len(par))}
    pores = make_pores(spec, seeds.rng(cfg_seed, "geometry.analytic.pores")) if "nucleus" in spec else np.zeros((0, 6))
    w.set_pores(pores, np.zeros(len(pores), np.int8), source_names=("procedural",), frame="model",
                min_sep_nm=float(spec["nucleus"]["pores"]["min_sep_nm"]),
                provenance="uniform directions on the nucleus ellipsoid, Newton-projected onto ne_mid = 0, normal = grad")
    rule = label_rule(set(grids), sections, clip=False)
    used = set(_rule_labels(rule))
    eng_labels = {c: nm for c, nm in ENGINE_LABELS.items() if nm in used}
    mems, _ = membranes_table(params["MEM_HALF"], "sdf_sections_v1 defaults (Pogozheva 2022 D_P-P/2)")
    R = float(np.linalg.norm(ca) + 100.0)
    gi = git_info()
    man = {"cell_id": spec["cell_id"], "geometry_id": spec["geometry_id"],
           "description": "synthetic analytic test cell (CI + demo)", "organism": {"name": "synthetic", "taxon": None},
           "source": {"dataset": "analytic", "spec": spec},
           "generator": {"toolkit": "toolkits/analytic/generate.py", "repo_commit": gi["commit"], "dirty": gi["dirty"],
                         "seed": int(spec["seed"]), "seed_streams": ["geometry.analytic.pores", "geometry.analytic.cristae"]},
           "engine": {"name": e.ENGINE_NAME, "layout_id": e.LAYOUT_ID, "params": params,
                      "field_slots": {nm: k for k, nm in enumerate(e.FIELD_NAMES)}},
           "frame": {"units": "nm", "handedness": "right", "axes": "x, y, z", "origin": "cell centre",
                     "source_affines": {}},
           "domain": {"enclosing_radius_nm": R, "clip_planes": []},
           "labels": labels_table(eng_labels), "membranes": mems, "label_rule": rule, "summary": summ}
    return w.finalize(man, validate_kw=validate_kw)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("spec", nargs="?", default=DEFAULT_SPEC)
    ap.add_argument("--dest", default=None)
    a = ap.parse_args(argv)
    spec = yaml.safe_load(open(a.spec))
    man, res = generate(spec, a.dest)
    print({"package": man["package_id"], "pass": res["pass"], "manifest_sha256": res["manifest_sha256"],
           "failed": res["failed_checks"]})
    return 0 if res["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
