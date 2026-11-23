"""Helpers shared by the Jurkat rule sets: data files, legacy change tables -> canonical changes.tsv rows."""
import os
import re

import pandas as pd
import yaml

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")
CH_COLS = ["rule_id", "species_id", "genes", "field", "old", "new", "evidence", "check", "mass_moved_pg"]


def data(name):
    return yaml.safe_load(open(os.path.join(DATA, name)))


def read_list(path):
    return {l.strip() for l in open(path) if l.strip() and not l.startswith("#")}


def normalise_changes(L: pd.DataFrame, version: str) -> pd.DataFrame:
    """legacy change table (v4: ing_id,gene,field,old,new,reason | v5: ...,rule | v6/v7: version,gene,ingredient,rule_id,
    field,old_value,new_value,evidence,motivating_check,mass_moved_pg) -> canonical rows with a consistent rule_id
    '<version>.<legacy rule>' (fixes the ledger bug: build_ledger looked for `rule`/`reason` and fell back to the last
    column, so v6/v7 rule_ids held mass_moved_pg)."""
    if len(L) == 0:
        return pd.DataFrame(columns=CH_COLS)
    out = pd.DataFrame(index=L.index)
    if "reason" in L.columns:                                   # v4
        pre = L.reason.astype(str).str.split(":").str[0].str.strip()
        out["rule_id"] = version + "." + pre.map(lambda s: re.sub(r"[^A-Za-z0-9]+", "_", s))
        out["species_id"], out["genes"] = L.ing_id, L.gene
        out["field"], out["old"], out["new"] = L.field, L.old, L.new
        out["evidence"], out["check"], out["mass_moved_pg"] = L.reason, "", ""
    elif "rule" in L.columns:                                   # v5
        out["rule_id"] = version + "." + L.rule.astype(str)
        out["species_id"], out["genes"] = L.ing_id, L.gene
        out["field"], out["old"], out["new"] = L.field, L.old, L.new
        out["evidence"], out["check"], out["mass_moved_pg"] = "", "", ""
    else:                                                       # v6 / v7
        out["rule_id"] = version + "." + L.rule_id.astype(str)
        out["species_id"], out["genes"] = L.ingredient, L.gene
        out["field"], out["old"], out["new"] = L.field, L.old_value, L.new_value
        out["evidence"], out["check"], out["mass_moved_pg"] = L.evidence, L.motivating_check, L.mass_moved_pg
    return out[CH_COLS].reset_index(drop=True)


# ---------------------------------------------------------------- preset common (B5): re-expressed rules, no pack read
NUC_KEYS = ("nucleus", "nucleolus", "speckle")              # soluble nuclear keys (capacity check)
CYT_KEYS = ("cytosol", "cytoskeleton", "unassigned")        # soluble cytoplasmic keys (capacity check)


def geometry_volumes(ctx):
    """(V nucleoplasm, V cytosol) in um3 of this build's geometry (<build>/geometry/measures.json, region_volumes_mc)"""
    import json
    p = os.path.join(ctx.bdir, "geometry", "measures.json")
    if not os.path.exists(p):
        raise FileNotFoundError(f"preset common: {p} needed (Vn / Vc of this build's geometry); run the geometry stage first")
    v = json.load(open(p))["region_volumes_mc"]["volume_um3"]
    return float(v["nucleoplasm"]), float(v["cytosol"])


def soluble_masses(I):
    """intended soluble nuclear / cytoplasmic protein mass (pg) of a species table (comp_weights x mass_pg)"""
    import json as _j
    mn = mc = 0.0
    for cw, m in zip(I.comp_weights, I.mass_pg):
        w = _j.loads(cw) if isinstance(cw, str) else {}
        mn += m * sum(w.get(k, 0.0) for k in NUC_KEYS)
        mc += m * sum(w.get(k, 0.0) for k in CYT_KEYS)
    return mn, mc


def capacity_scale(I, raises, vn, vc):
    """capacity check of a rule (criteria B5 §2): the rule may not raise the intended nucleoplasm protein density above
    the intended cytosol density (before the rule). raises = [(row, added nuclear share)] -> (scale in [0, 1], info)"""
    mn, mc = soluble_masses(I)
    add = sum(float(I.at[i, "mass_pg"]) * d for i, d in raises if d > 0)
    cap = mc / vc * vn - mn                     # nuclear mass that brings the nucleus to the cytosol density
    sc = 1.0 if add <= 0 or add <= cap else max(cap, 0.0) / add
    return sc, dict(nuclear_pg=mn, cytosol_pg=mc, vn_um3=vn, vc_um3=vc, nuclear_mg_ml=mn / vn * 1e3,
                    cytosol_mg_ml=mc / vc * 1e3, raise_pg=add, room_pg=cap, scale=sc)
