"""Generic (cell-agnostic) inventory components, registered by name:

  adapter.proteome     proteome_table      TSV: gene, uniprot, mass_da, abundance [, unit] (relative, e.g. ppm) or copies
  adapter.localisation localisation_table  TSV: gene, compartments ("key:w;key:w"), transmembrane [, n_tm, signal_peptide,
                                           mito_om] (bools)
  adapter.assembly     none                no assemblies (assembly/complexed.tsv = v0 with copies_in_complexes 0)
  adapter.assembly     complex_table       TSV: complex_id, subunits ("gene:stoich;..."), copies, struct_src; complexes
                                           consume subunit pools greedily in table order (clipped at the pool);
                                           assembly_options.complex_copies: measured_free_pool takes copies from SEC-MS
  inventory_builder    generic.monomers    assemblies (assembly/complexes.tsv) + free monomers -> pack-table species
  inventory_rule       generic.signal_lumen soluble species without a signal peptide (signal_frac < threshold) lose their
                                           lumenal keys to the cytosolic key (makes the packer's lumen rule explicit)

TSVs may start with '#' comment lines (README / provenance). Paths come from inventory.params.<component>.path.
"""
from __future__ import annotations

import json

import numpy as np
import pandas as pd

from ..registry import register
from . import schema

AV = schema.AVOGADRO


def _read(path):
    return pd.read_csv(path, sep="\t", comment="#")


def _parse_w(s):
    if not isinstance(s, str) or not s:
        return {}
    return {k: float(v) for k, v in (t.split(":") for t in s.split(";"))}


def _wjson(w, nd=4):
    return json.dumps({k: round(float(v), nd) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 10 ** -nd})


# ------------------------------------------------------------------------------------------------------- adapters
def proteome_table(ctx, total_protein_pg):
    p = ctx.params("proteome_table")
    T = _read(ctx.resolve(p["path"])).drop_duplicates("gene").set_index("gene").sort_index()
    inv = pd.DataFrame(index=T.index)
    inv["uniprot"] = T.uniprot
    inv["mass_da"] = T.mass_da.astype(float)
    src = p.get("source_id", "proteome_table")
    if "copies" in T.columns and T.copies.notna().all():
        inv["source"] = src
        inv["copies"] = T.copies.round().astype(np.int64)
        tot = (inv.copies * inv.mass_da).sum() / AV * 1e12
        if abs(tot - total_protein_pg) > 0.02 * total_protein_pg:
            raise ValueError(f"proteome_table: absolute copies sum to {tot:.4f} pg, inventory.total_protein_pg is "
                             f"{total_protein_pg} (> 2% apart)")
    else:
        frac = T.abundance / T.abundance.sum()
        mean_mass = (frac * inv.mass_da).sum()
        N = total_protein_pg * 1e-12 / (mean_mass / AV)
        inv["source"] = src
        inv["abundance"] = T.abundance
        inv["copies"] = (frac * N).round().astype(np.int64)
        inv = inv[["uniprot", "mass_da", "abundance", "source", "copies"]]
    ctx.say(f"proteome_table: {len(inv)} genes, {inv.copies.sum():,} copies, "
            f"{(inv.copies * inv.mass_da).sum() / AV * 1e12:.4f} pg (target {total_protein_pg} pg)")
    return inv


proteome_table.PARAM_FILES = ["path"]
register("adapter.proteome", "proteome_table")(proteome_table)


def localisation_table(ctx, prot):
    p = ctx.params("localisation_table")
    T = _read(ctx.resolve(p["path"])).drop_duplicates("gene").set_index("gene")
    out = pd.DataFrame(index=prot.index)
    out["compartments"] = T.compartments.reindex(out.index).fillna("unassigned:1.0")
    out["loc_source"] = p.get("source_id", "localisation_table")
    out["transmembrane"] = T.transmembrane.reindex(out.index).fillna(False).astype(bool).values
    for c in ("n_tm", "signal_peptide", "mito_om"):
        if c in T.columns:
            out[c] = T[c].reindex(out.index).fillna(0 if c == "n_tm" else False).values
    return out


localisation_table.PARAM_FILES = ["path"]
register("adapter.localisation", "localisation_table")(localisation_table)


def _complexed(inv, used):
    inv = inv.copy()
    inv["copies_in_complexes"] = pd.Series(used).reindex(inv.index).fillna(0).astype(np.int64)
    inv["frac_in_complexes"] = (inv.copies_in_complexes / inv.copies.replace(0, np.nan)).round(3)
    return inv


def assembly_none(ctx):
    inv = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t", index_col="gene")
    ctx.write_tsv(pd.DataFrame(columns=["complex_id", "genes", "stoich", "copies_built", "struct_src", "mass_da"]),
                  "assembly/complexes.tsv")
    ctx.write_tsv(_complexed(inv, {}), "assembly/complexed.tsv", index=True)


assembly_none.OUTPUTS = ["assembly/complexes.tsv", "assembly/complexed.tsv"]
register("adapter.assembly", "none")(assembly_none)


def assembly_complex_table(ctx):
    """params: path; assembly_options.complex_copies = stoichiometry_cap (default: the table's copies, clipped at the
    pools) | measured_free_pool (copies from the SEC-MS measured free pools when the complex is measured, see
    wholecell/inventory/secms.py; needs inventory.adapters.secms)"""
    from . import secms as S2
    p = ctx.params("complex_table")
    opts = S2.assembly_options(p)
    view = S2.SecView(ctx, opts["secms"]) if opts["complex_copies"] == "measured_free_pool" else None
    inv = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t", index_col="gene")
    C = _read(ctx.resolve(p["path"]))
    left = inv.copies.to_dict()
    pool = dict(left)
    acc_of = inv.uniprot.to_dict()
    used, rows = {}, []
    for r in C.itertuples(index=False):
        sub = {t.split(":")[0]: int(t.split(":")[1]) for t in r.subunits.split(";")}
        if any(g not in left for g in sub):
            ctx.say(f"complex {r.complex_id}: subunit missing from the inventory, skipped"); continue
        want, how = int(r.copies), "table"
        if view is not None:
            present = {g: s for g, s in sub.items() if pool[g] > 0}
            m, why = S2.measured_copies(view, r.complex_id, present, pool, acc_of, opts)
            if m is not None:
                want, how = int(m), f"measured_free_pool ({why})"
        n = int(min([want] + [left[g] // s for g, s in sub.items()]))
        for g, s in sub.items():
            left[g] -= n * s; used[g] = used.get(g, 0) + n * s
        rows.append(dict(complex_id=r.complex_id, genes=";".join(sub), stoich=";".join(str(s) for s in sub.values()),
                         copies_built=n, struct_src=r.struct_src,
                         mass_da=float(sum(inv.mass_da[g] * s for g, s in sub.items()))))
        ctx.say(f"complex {r.complex_id}: {n:,} copies (requested {want:,}, {how}; table {int(r.copies):,})")
    ctx.write_tsv(pd.DataFrame(rows, columns=["complex_id", "genes", "stoich", "copies_built", "struct_src", "mass_da"]),
                  "assembly/complexes.tsv")
    ctx.write_tsv(_complexed(inv, used), "assembly/complexed.tsv", index=True)


assembly_complex_table.OUTPUTS = ["assembly/complexes.tsv", "assembly/complexed.tsv"]
assembly_complex_table.PARAM_FILES = ["path"]
register("adapter.assembly", "complex_table")(assembly_complex_table)


# ------------------------------------------------------------------------------------------------------- builder
def monomers(ctx):
    """species = built assemblies + free monomers. comp_weights of an assembly = mass-weighted mean of its subunits;
    signal_frac / mito_om_frac = subunit-mass fraction flagged (localisation adapter columns signal_peptide / mito_om)."""
    inv = pd.read_csv(ctx.need("assembly/complexed.tsv"), sep="\t", index_col="gene")
    C = pd.read_csv(ctx.need("assembly/complexes.tsv"), sep="\t")
    W = {g: (_parse_w(s) or {"unassigned": 1.0}) for g, s in inv.compartments.items()}
    mass = inv.mass_da.fillna(0.0)
    tm = inv.transmembrane.fillna(False).astype(bool)
    ntm = inv.n_tm.fillna(0).astype(int) if "n_tm" in inv.columns else tm.astype(int)
    sig = inv.signal_peptide.fillna(False).astype(bool) if "signal_peptide" in inv.columns else pd.Series(False, inv.index)
    om = inv.mito_om.fillna(False).astype(bool) if "mito_om" in inv.columns else pd.Series(False, inv.index)
    rows = []
    for r in C.itertuples(index=False):
        if int(r.copies_built) <= 0:
            continue
        sub = dict(zip(r.genes.split(";"), map(int, str(r.stoich).split(";"))))
        tot = sum(mass[g] * s for g, s in sub.items()) or 1.0
        w = {}
        for g, s in sub.items():
            for k, v in W[g].items():
                w[k] = w.get(k, 0.0) + v * mass[g] * s / tot
        mf = lambda flag: round(sum(mass[g] * s for g, s in sub.items() if flag[g]) / tot, 3)
        rows.append(dict(ing_id=r.complex_id, kind="complex", genes=r.genes, stoich=str(r.stoich),
                         copies=int(r.copies_built), mass_da=float(tot), comp_weights=_wjson(w, 3),
                         transmembrane=bool(any(tm[g] for g in sub)), n_tm=int(sum(ntm[g] * s for g, s in sub.items())),
                         struct_src=r.struct_src, uniprot="", signal_frac=mf(sig), mito_om_frac=mf(om)))
    free = (inv.copies - inv.copies_in_complexes.fillna(0)).clip(lower=0).astype(np.int64)
    for g in inv.index[free > 0]:
        rows.append(dict(ing_id=g, kind="protein", genes=g, stoich="1", copies=int(free[g]), mass_da=float(mass[g]),
                         comp_weights=_wjson(W[g], 3), transmembrane=bool(tm[g]), n_tm=int(ntm[g]),
                         struct_src=f"afdb_monomer:{inv.uniprot[g]}", uniprot=inv.uniprot[g],
                         signal_frac=float(sig[g]), mito_om_frac=float(om[g])))
    I = pd.DataFrame(rows, columns=schema.PACK_TABLE[:-1])
    I["mass_pg"] = I.copies * I.mass_da / AV * 1e12
    ctx.say(f"generic.monomers: {len(I)} species, {I.copies.sum():,} copies, {I.mass_pg.sum():.4f} pg")
    return I


monomers.PREV = ["assembly/complexed.tsv", "assembly/complexes.tsv"]
register("inventory_builder", "generic.monomers")(monomers)


# ------------------------------------------------------------------------------------------------------- generic rule
def signal_lumen(ctx, src):
    """lumenal keys of soluble species without a signal peptide -> cytosolic key (params: lumen_keys, cytosol_key,
    threshold, rule_id)"""
    p = {"lumen_keys": ["er", "golgi", "endolysosomal", "secreted"], "cytosol_key": "cytosol", "threshold": 0.5,
         "rule_id": "generic.signal_lumen"} | ctx.params("generic.signal_lumen")
    I = pd.read_csv(src, sep="\t")
    log = []
    lk = set(p["lumen_keys"])
    for i in range(len(I)):
        if bool(I.at[i, "transmembrane"]) or float(I.at[i, "signal_frac"]) >= p["threshold"]:
            continue
        w = json.loads(I.at[i, "comp_weights"])
        mv = sum(v for k, v in w.items() if k in lk)
        if mv <= 0:
            continue
        nw = {k: v for k, v in w.items() if k not in lk}
        nw[p["cytosol_key"]] = nw.get(p["cytosol_key"], 0) + mv
        new = _wjson(nw, 4)
        log.append(dict(rule_id=p["rule_id"], species_id=I.at[i, "ing_id"], genes=I.at[i, "genes"], field="comp_weights",
                        old=I.at[i, "comp_weights"], new=new,
                        evidence=f"signal_frac {I.at[i, 'signal_frac']} < {p['threshold']}: lumenal share {mv:.3f} -> "
                                 f"{p['cytosol_key']}", check="", mass_moved_pg=round(mv * float(I.at[i, "mass_pg"]), 6)))
        I.at[i, "comp_weights"] = new
    ctx.say(f"generic.signal_lumen: {len(log)} species changed")
    return {"ingredients": I, "changes": pd.DataFrame(log, columns=[c for c in schema.CHANGES if c not in ("step", "rule_set")])}


register("inventory_rule", "generic.signal_lumen")(signal_lumen)
