"""This build's library against a reference build's library (report only; library.compare_reference).

Per species of this build's species table: structure source (npz `src`, else the species struct_src), chain set (the
chain record of the library chain rule; `-` for models without one), bead count, radius of gyration of the beads and, for
membrane species (npz `membrane` == TM) placed at the plasma membrane (pack-table comp_weights key plasma_membrane > 0),
the ectodomain height above the membrane. Changes above the thresholds are flagged; nothing fails.

Ectodomain height (simple structural definition, the structural part of the synapse kit's rule): library frame membrane
normal +z = cytoplasm, mid-plane z = 0, so the exoplasmic side is z < -2 nm (outer leaflet at -2 nm). The exit point is
the centre of the beads with -3.5 < z < -1.5 nm placed on z = -2 (else the centre of the exoplasmic beads at z = -2);
ecto_h_nm = the largest distance of an exoplasmic bead from the exit point (synapse kit h_up_nm), ecto_depth_nm = the
largest depth below the outer leaflet (-z_min - 2). The kit's full height adds disordered-tail extensions and a 2 nm
N-glycan cap and checks the side with the UniProt topology; that is not repeated here (a flipped model shows as a change).

Config: library.compare_reference: {build: <reference build directory ({results_root}/... or a path)>,
                                     thresholds: {height_nm: 1.0, rg_frac: 0.10, beads_frac: 0.10}}
Writes <library>/compare_reference.tsv and compare_reference.json (counts per flag).
"""
from __future__ import annotations

import json
import os

import numpy as np

HALF = 2.0
PM_KEYS = ("plasma_membrane",)
COLS = ["ing_id", "kind", "pm", "src_ref", "src_new", "src_changed", "chains_ref", "chains_new", "chains_changed",
        "beads_ref", "beads_new", "rg_ref", "rg_new", "rg_rel", "membrane_ref", "membrane_new", "ecto_h_ref", "ecto_h_new",
        "d_ecto_h", "ecto_depth_ref", "ecto_depth_new", "flags"]


def geom(path):
    """{src, chains, beads, rg, membrane, ecto_h, ecto_depth} of one library npz (None if absent)"""
    if not path or not os.path.exists(path):
        return None
    with np.load(path, allow_pickle=False) as z:
        b = z["beads"].astype(float)
        src = str(z["src"]) if "src" in z.files else ""
        mem = str(z["membrane"]) if "membrane" in z.files else ""
        ch = "-"
        if "chain_set" in z.files:
            rec = json.loads(str(z["chain_set"]))
            ch = ",".join(f"{d['chain']}:{d['acc'] or 'na'}" for d in rec["kept"])
    out = dict(src=src, chains=ch, beads=len(b), rg=float(np.sqrt(((b - b.mean(0)) ** 2).sum(1).mean())) if len(b) else 0.0,
               membrane=mem, ecto_h=np.nan, ecto_depth=np.nan)
    if mem == "TM":
        zz = b[:, 2]
        lo = zz < -HALF
        if lo.any():
            ex = b[(zz < -HALF + 0.5) & (zz > -HALF - 1.5)]
            ep = ex.mean(0) if len(ex) else np.array([b[lo, 0].mean(), b[lo, 1].mean(), -HALF])
            ep[2] = -HALF
            out["ecto_h"] = float(np.linalg.norm(b[lo] - ep, axis=1).max())
            out["ecto_depth"] = float(max(0.0, -zz.min() - HALF))
        else:
            out["ecto_h"] = out["ecto_depth"] = 0.0
    return out


def manifest(lib_out):
    """{species: npz path} of a library folder's library_manifest.tsv"""
    f = os.path.join(lib_out, "library_manifest.tsv")
    out = {}
    with open(f) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for ln in fh:
            r = dict(zip(head, ln.rstrip("\n").split("\t")))
            if r["path"].endswith(".npz"):
                out[r["species"]] = r["path"]
    return out


def pm_ids(pack_table):
    """species with a plasma-membrane share in the pack table (comp_weights JSON)"""
    import pandas as pd
    if not pack_table or not os.path.exists(pack_table):
        return set()
    T = pd.read_csv(pack_table, sep="\t", usecols=["ing_id", "comp_weights"], keep_default_na=False)
    out = set()
    for i, w in zip(T.ing_id, T.comp_weights):
        try:
            d = json.loads(w) if w else {}
        except ValueError:
            continue
        if any(float(d.get(k, 0) or 0) > 0 for k in PM_KEYS):
            out.add(i)
    return out


def compare(rows, new_path, ref_path, pm=(), thresholds=None, ref_src=None):
    """rows: [{ing_id, kind, struct_src}] of this build; new_path / ref_path: {id: npz}; ref_src: {id: struct_src} of the
    reference species table (when its npz carries no src). Returns a DataFrame (COLS)."""
    import pandas as pd
    t = dict(height_nm=1.0, rg_frac=0.10, beads_frac=0.10)
    t.update(thresholds or {})
    out = []
    for r in rows:
        i = r["ing_id"]
        g1, g0 = geom(new_path.get(i)), geom(ref_path.get(i))
        s1 = (g1 or {}).get("src") or r.get("struct_src", "")
        s0 = (g0 or {}).get("src") or (ref_src or {}).get(i, "")
        d = dict(ing_id=i, kind=r.get("kind", ""), pm=i in pm, src_ref=s0, src_new=s1)
        fl = []
        if g1 is None or g0 is None:
            fl.append("new_only" if g0 is None and g1 is not None else "ref_only" if g1 is None and g0 is not None else "none")
        for k, g in (("ref", g0), ("new", g1)):
            g = g or {}
            d.update({f"chains_{k}": g.get("chains", ""), f"beads_{k}": g.get("beads", np.nan), f"rg_{k}": g.get("rg", np.nan),
                      f"membrane_{k}": g.get("membrane", ""), f"ecto_h_{k}": g.get("ecto_h", np.nan),
                      f"ecto_depth_{k}": g.get("ecto_depth", np.nan)})
        d["src_changed"] = bool(s0 and s1 and s0 != s1)
        d["chains_changed"] = bool(g0 and g1 and d["chains_ref"] != "-" and d["chains_new"] != "-" and d["chains_ref"] != d["chains_new"])
        d["rg_rel"] = (d["rg_new"] - d["rg_ref"]) / d["rg_ref"] if g0 and g1 and d["rg_ref"] > 0 else np.nan
        d["d_ecto_h"] = d["ecto_h_new"] - d["ecto_h_ref"]
        if d["src_changed"]:
            fl.append("source")
        if d["chains_changed"]:
            fl.append("chains")
        if g0 and g1 and g0["beads"] and abs(g1["beads"] - g0["beads"]) / g0["beads"] > t["beads_frac"]:
            fl.append("beads")
        if np.isfinite(d["rg_rel"]) and abs(d["rg_rel"]) > t["rg_frac"]:
            fl.append("rg")
        if d["pm"] and np.isfinite(d["d_ecto_h"]) and abs(d["d_ecto_h"]) > t["height_nm"]:
            fl.append("pm_height")
        if g0 and g1 and (g0["membrane"] == "TM") != (g1["membrane"] == "TM"):
            fl.append("membrane")
        d["flags"] = ";".join(fl)
        out.append(d)
    return pd.DataFrame(out, columns=COLS)


def run(X):
    """library.compare_reference on a sealed library (build.step_seal)"""
    import pandas as pd
    from .build import _path, log, species
    c = X.lc["compare_reference"]
    ref = _path(X.cfg, X.bdir, c["build"])
    rlib = os.path.join(ref, "library")
    ref_path = manifest(rlib)
    new_path = {k: os.path.join(X.lib, os.path.basename(v)) for k, v in manifest(X.out).items()}
    ref_src = {}
    rs = os.path.join(rlib, "species.tsv")
    if os.path.exists(rs):
        T = pd.read_csv(rs, sep="\t", usecols=["ing_id", "struct_src"], dtype=str, keep_default_na=False)
        ref_src = dict(zip(T.ing_id, T.struct_src))
    sp = X.get("species_spec") or {}
    pt = os.path.join(X.bdir, sp.get("pack_table", "inventory/pack_table.tsv")) if sp else os.path.join(X.bdir, "inventory", "pack_table.tsv")
    R = compare(species(X).to_dict("records"), new_path, ref_path, pm_ids(pt), c.get("thresholds"), ref_src)
    R.to_csv(os.path.join(X.out, "compare_reference.tsv"), sep="\t", index=False)
    fl = R["flags"].str.split(";").explode()
    summ = dict(reference=ref, species=len(R), flagged=int((R["flags"] != "").sum()),
                by_flag={k: int(v) for k, v in fl[fl != ""].value_counts().items()},
                pm_height_changed=R.ing_id[R["flags"].str.contains("pm_height")].tolist()[:50],
                thresholds=dict(dict(height_nm=1.0, rg_frac=0.10, beads_frac=0.10), **(c.get("thresholds") or {})))
    json.dump(summ, open(os.path.join(X.out, "compare_reference.json"), "w"), indent=1)
    X.results["compare_reference"] = {k: summ[k] for k in ("reference", "species", "flagged", "by_flag")}
    log(f"library vs reference {ref}: {summ['flagged']} of {len(R)} species flagged {summ['by_flag']}")
    return R
