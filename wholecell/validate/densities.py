# origin: jurkat/validate/densities.py (legacy workspace, 2026-11-15)
# Changed: region list, shell->volume map, pre-pass library->region map, nucleus group and mRNA sources come from the
# config (validate.densities) and the build's region table instead of jpack_test / hard-coded dicts; pack files are the
# build's own (common.Inputs). VAL_GEOM/VAL_POLY defect fixed: the mRNA mass is taken from the build's own mRNA inputs
# (validate.inputs.<name> / preplaced manifest), never from a fixed legacy file. Library bead masses come from the
# validation library cache. Output values and rounding are those of the legacy script.
"""Macromolecular mass density per volume region of the packed cell.

mass = random-packed protein placed per region (tiles: job_lib x job_placed x library mass; shells added to the volume
region they occupy, virtual regions to their host volume) + pre-pass polymers (library bead mass, aa 110 Da / nt 330 Da,
assigned to a volume region by config) + RNA (nt x 330 Da, per configured source and region);
volume = region volumes from tiles.parquet. Output: <build>/validate/densities.json
"""
from __future__ import annotations

import fnmatch
import json
import os

import numpy as np
import pandas as pd

from . import common, regions as R_

AV = 6.022e23


def _tile_mass(args):
    files, mass, nreg = args
    out = np.zeros(nreg)
    for f in files:
        z = np.load(f)
        np.add.at(out, z["job_region"].astype(int), z["job_placed"] * mass[z["job_lib"]])
    return out


def _match(name, table, default):
    for pat, reg in (table or {}).items():
        if fnmatch.fnmatchcase(name, pat):
            return reg
    return default


def run(cfg, bdir, inp: common.Inputs, out_path, overwrite=False, nproc=8):
    dc = (cfg.get("validate", {}) or {}).get("densities", {}) or {}
    regs = R_.load_regions(cfg, inp)
    RN = R_.names(regs)
    vol_i = R_.by_kind(regs, "volume", "virtual")
    mem_i = R_.by_kind(regs, "membrane")
    sh_i = R_.by_kind(regs, "shell")
    default_region = dc.get("default_region", RN[vol_i[0]])
    mass, mass_src = common.lib_masses(inp, bdir, nproc=nproc)
    fs = inp.tile_files()
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(nproc) as ex:
        packed = sum(ex.map(_tile_mass, [(fs[i::64], mass, len(RN)) for i in range(min(64, len(fs)))]))
    M = pd.read_parquet(inp["tiles_parquet"])
    vol = {RN[i]: (float(M[RN[i]].sum()) / 1e9 if RN[i] in M.columns else 0.0) for i in vol_i}     # um3
    virt_to = dc.get("virtual_to", default_region)
    for i in R_.by_kind(regs, "virtual"):                 # virtual-region mass belongs to its host volume (volumes overlap)
        host = virt_to.get(RN[i], default_region) if isinstance(virt_to, dict) else virt_to
        packed[RN.index(host)] += packed[i]; packed[i] = 0.0; vol[RN[i]] = 0.0
    s2v = dc.get("shell_to_volume", {}) or {}
    split = dc.get("shell_split", {}) or {}
    for k in sh_i:
        if RN[k] in s2v:
            packed[RN.index(s2v[RN[k]])] += packed[k]
        for r_, f_ in (split.get(RN[k], {}) or {}).items():
            packed[RN.index(r_)] += f_ * packed[k]
    # pre-pass mass by library
    pre = {}
    pre_lib_mass = {}
    if inp.has("preplaced_libs") and inp.has("preplaced_types"):
        PL = json.load(open(inp["preplaced_libs"]))
        pt = np.load(inp["preplaced_types"], mmap_mode="r")
        cnt = np.bincount(np.asarray(pt), minlength=len(PL))
        # pre-pass entries by group (pack/preplaced_groups.json, parallel to preplaced_libs; yeast port): group -> region
        gv = dc.get("preplaced_group_volume") or {}
        grp = json.load(open(inp["preplaced_groups"])) if (gv and inp.has("preplaced_groups")) else None
        lc = common.lib_cache(inp, PL, common.cache_path(bdir, "lib_preplaced.npz"), nproc=min(nproc, 4))
        for i, (l, c) in enumerate(zip(PL, cnt)):
            m = float(lc.bead_mass[lc.idx[l]])
            k_ = f"{grp[i]}[{grp[:i].count(grp[i])}]" if grp is not None else None
            if grp is not None and (k_ in gv or grp[i] in gv):
                g_ = gv.get(k_, gv.get(grp[i]))
                split = g_ if isinstance(g_, dict) else {g_: 1.0}
            else:
                split = {_match(l, dc.get("preplaced_volume"), default_region): 1.0}
            for reg, f_ in split.items():
                pre[reg] = pre.get(reg, 0.0) + f_ * c * m
            key = l if l not in pre_lib_mass else f"{l}#{i}"
            pre_lib_mass[key] = dict(region=split if len(split) > 1 else next(iter(split)), copies=int(c), da_per_copy=m)
    # RNA (mRNA, nascent RNA): the build's own sources
    rna = {}
    rna_src = []
    for s in dc.get("mrna_sources", []) or []:
        p = inp.lookup_input(s["input"])
        if not p or not os.path.exists(p):
            raise FileNotFoundError(f"densities: RNA source '{s['input']}' not found ({p})")
        z = np.load(p, mmap_mode="r")
        nt = float(np.asarray(z["nt"]).sum())
        rna[s["region"]] = rna.get(s["region"], 0.0) + nt
        rna_src.append(dict(input=s["input"], path=p, md5=common.md5(p), nt=nt, region=s["region"]))
    if not rna_src and inp.has("mrna_pts") and float(dc.get("mrna_nt_per_bead", 0) or 0) > 0:
        nb_ = int(np.load(inp["mrna_off"])[-1]) if inp.has("mrna_off") else len(np.load(inp["mrna_pts"], mmap_mode="r"))
        reg = dc.get("mrna_region", default_region)
        rna[reg] = nb_ * float(dc["mrna_nt_per_bead"])
        rna_src.append(dict(input="mrna_pts (indexed beads x mrna_nt_per_bead)", nt=rna[reg], region=reg))
    rows = []
    for i in vol_i:
        r = RN[i]
        if vol[r] <= 0:
            continue
        prot = packed[i]; extra = pre.get(r, 0.0) + rna.get(r, 0.0) * 330.0
        rows.append(dict(region=r, volume_um3=round(vol[r], 1), packed_pg=round(prot / AV * 1e12, 2), prepass_pg=round(extra / AV * 1e12, 2),
                         mg_per_ml=round((prot + extra) / AV / (vol[r] * 1e-12) * 1e3, 1),
                         packed_only_mg_per_ml=round(prot / AV / (vol[r] * 1e-12) * 1e3, 1)))
    T = pd.DataFrame(rows)
    Ti = T.set_index("region")
    out = dict(rows=rows)
    den = dc.get("ratio_denominator", default_region)
    grp = [r_ for r_ in (dc.get("nucleus_group") or []) if r_ in Ti.index]
    if grp and den in Ti.index:
        c = Ti.loc[den, "mg_per_ml"]
        n = float((Ti.loc[grp, "packed_pg"] + Ti.loc[grp, "prepass_pg"]).sum() / Ti.loc[grp, "volume_um3"].sum() * 1e3)
        out["nuc_cyto_ratio"] = round(n / c, 3)
    out["mem_pg"] = {RN[i]: round(packed[i] / AV * 1e12, 3) for i in mem_i}
    out["shell_pg"] = {RN[i]: round(packed[i] / AV * 1e12, 3) for i in sh_i}
    out["mrna_pg"] = round(sum(rna.values()) * 330 / AV * 1e12, 2)
    out["mrna_pg_by_region"] = {k: round(v * 330 / AV * 1e12, 4) for k, v in rna.items()}
    out["packed_mass_source"] = mass_src
    out["inputs"] = dict(tiles_dir=inp["tiles_dir"], tiles_parquet=inp["tiles_parquet"], libs=inp["libs"],
                         prepass_libs=inp.get("preplaced_libs"), rna_sources=rna_src, prepass_libraries=pre_lib_mass)
    common.write_json(out_path, out, overwrite)
    return out
