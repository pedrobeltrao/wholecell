"""Canonical inventory tables (spec: docs/reference/inventory_schema.md) from the final pack table, the per-gene inventory,
the base species table and every rule step's changes.tsv. Cell-agnostic: everything cell-specific comes from config
(inventory.canonical, proteome_sources, localisation_sources, frozen_inputs metadata) or from adapter outputs.

Writes inventory/: species.tsv, species_loc.tsv, field_provenance.tsv, sources.tsv, proteome_sources.tsv,
localisation_sources.tsv, gene_abundance.tsv, secms_protein.tsv, secms_complex.tsv, surface.tsv, changes_all.tsv,
version_changes.tsv, secms_bound_fraction.tsv (+ secms_bound_fraction_summary.json; secms_bound.py, written last).
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from .. import registry
from . import schema

FIELDS = ["copies", "comp_weights", "transmembrane", "n_tm", "mass_da", "genes", "stoich", "signal_frac", "mito_om_frac",
          "struct_src"]


def _w(s):
    return json.loads(s) if isinstance(s, str) and s.startswith("{") else {}


def species_key(r, g2u):
    """stable key grammar: prot:<UniProt> | cpx:<CPX id> | curated:<id> | homo:<UniProt>:x<n> | afhomo:<UniProt>:x2 |
    het:<UniProt A>+<UniProt B> (fallback to gene symbols when no accession)"""
    gs = str(r.genes).split(";")
    acc = [g2u.get(g) if isinstance(g2u.get(g), str) else g for g in gs]
    if r.kind == "protein":
        return f"prot:{acc[0]}"
    if r.kind == "oligomer":
        return f"homo:{acc[0]}:x{str(r.stoich).split(';')[0]}"
    if r.kind == "afdb_homodimer":
        return f"afhomo:{acc[0]}:x2"
    if r.kind == "afdb_heterodimer":
        return "het:" + "+".join(acc)
    if str(r.ing_id).startswith("CPX-"):
        return f"cpx:{r.ing_id}"
    return f"curated:{r.ing_id}"


def _rule_changes(ctx):
    rs = ctx.inv.get("rule_sets", []) or []
    out = []
    for n in rs:
        p = ctx.path(f"rules/{n}/changes.tsv")
        if os.path.exists(p):
            out.append(pd.read_csv(p, sep="\t", dtype=str, keep_default_na=False))
    return pd.concat(out, ignore_index=True) if out else schema.empty(schema.CHANGES)


def _version_changes(tables, changes, key_comp):
    """gene-level target-share shifts between consecutive species tables (coarse compartments) with the rule ids taken
    from the `rule_id` column of the change tables (the legacy ledger read the wrong column for v6/v7)."""
    def gene_shares(T):
        acc = {}
        for r in T.itertuples():
            w = _w(r.comp_weights)
            for g, s in zip(str(r.genes).split(";"), str(r.stoich).split(";")):
                n = float(r.copies) * float(s)
                d = acc.setdefault(g, {"_n": 0.0})
                d["_n"] += n
                for k, v in w.items():
                    c = key_comp.get(k, k)
                    d[c] = d.get(c, 0.0) + n * v
        return acc
    rows = []
    prev_name, prev = None, None
    for name, T in tables:
        cur = gene_shares(T)
        if prev is not None:
            ch = changes[changes.step == name]
            by_gene = {}
            for r in ch.itertuples():
                for g in str(r.genes).split(";"):
                    if g:
                        by_gene.setdefault(g, []).append((r.rule_id, str(r.evidence)[:160]))
            for g in sorted(set(prev) | set(cur)):
                a, b = prev.get(g), cur.get(g)
                rr = by_gene.get(g)
                ids = ",".join(sorted({x for x, _ in rr})) if rr else ""
                ev = " | ".join(dict.fromkeys(e for _, e in rr if e))[:600] if rr else ""
                if a is None or b is None:
                    rows.append(dict(gene=g, transition=f"{prev_name}->{name}", change="added" if a is None else "removed",
                                     max_shift_target=np.nan, tvd_target=np.nan, shifts_target="",
                                     copies_from=a["_n"] if a else 0.0, copies_to=b["_n"] if b else 0.0,
                                     reason_class="inventory_rule" if rr else "gene_set", rule_ids=ids, evidence=ev))
                    continue
                na, nb = a["_n"], b["_n"]
                comps = (set(a) | set(b)) - {"_n"}
                sa = {c: a.get(c, 0.0) / na if na > 0 else 0.0 for c in comps}
                sb = {c: b.get(c, 0.0) / nb if nb > 0 else 0.0 for c in comps}
                d = {c: sb[c] - sa[c] for c in comps}
                mx = max([abs(x) for x in d.values()] or [0.0])
                if mx < 0.05 and not rr and na == nb:
                    continue
                rows.append(dict(gene=g, transition=f"{prev_name}->{name}",
                                 change="shift" if mx >= 0.1 else ("minor" if mx >= 0.05 else "no_target_shift"),
                                 max_shift_target=round(mx, 3), tvd_target=round(sum(abs(x) for x in d.values()) / 2, 3),
                                 shifts_target=";".join(f"{c}:{d[c]:+.2f}" for c in sorted(d, key=lambda c: (-abs(d[c]), c)) if abs(d[c]) >= 0.02),
                                 copies_from=round(na, 2), copies_to=round(nb, 2),
                                 reason_class="inventory_rule" if rr else "partner/assembly change", rule_ids=ids, evidence=ev))
        prev_name, prev = name, cur
    return pd.DataFrame(rows, columns=schema.VERSION_CHANGES)


def export(ctx):
    from .stage import steps
    cfg, inv = ctx.cfg, ctx.inv
    can = inv.get("canonical", {}) or {}
    key_comp = can.get("key_compartment", {}) or {}
    P = pd.read_csv(ctx.need("pack_table.tsv"), sep="\t")
    V0 = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t", dtype={"gene": str}).drop_duplicates("gene").set_index("gene")
    base = pd.read_csv(ctx.need("ingredients/ingredients.tsv"), sep="\t")
    ch = _rule_changes(ctx)
    ad = dict(inv.get("adapters", {}) or {})
    lobj = registry.get("adapter.localisation", ad["localisation"]) if ad.get("localisation") else None
    if hasattr(lobj, "provenance_label"):        # core localisation: provenance name from its params (label)
        ad["localisation"] = lobj.provenance_label(cfg)
    bname = inv.get("ingredients", "generic.monomers")
    g2u = V0.uniprot.to_dict() if "uniprot" in V0.columns else {}

    # ---- changes_all
    ctx.write_tsv(ch, "changes_all.tsv")
    sp_ch = ch[ch.field.isin(FIELDS)]
    last = {}                                                    # (species, field) -> last change row
    for r in sp_ch.itertuples():
        last[(r.species_id, r.field)] = r
    created = {}
    for r in ch[ch.field == "row"].itertuples():
        created[r.species_id] = r.step
    base_ids = set(base.ing_id)

    # ---- species_loc (one row per species x key; rule_id = last rule that changed that key's weight)
    key_rule = {}
    for r in ch[ch.field == "comp_weights"].itertuples():
        wo, wn = _w(r.old), _w(r.new)
        for k in set(wo) | set(wn):
            if abs(wo.get(k, 0.0) - wn.get(k, 0.0)) > 0:
                key_rule[(r.species_id, k)] = (r.rule_id, r.evidence, r.step)
    loc_rows, prov_rows, sp_rows = [], [], []
    for r in P.itertuples():
        w = _w(r.comp_weights)
        mono = r.kind == "protein"
        base_loc = (f"v0:{ad.get('localisation')}[{V0.loc_source.get(r.genes, '')}]" if mono and r.genes in V0.index
                    else f"ingredients:{bname}[subunit mass-weighted]")
        for k, v in w.items():
            rid, ev, st = key_rule.get((r.ing_id, k), (base_loc, "", "v0" if mono else "ingredients"))
            loc_rows.append(dict(species_id=r.ing_id, loc_key=k, weight=v, rule_id=rid, evidence=ev, step=st))
        for f in FIELDS:
            val = getattr(r, f)
            if (r.ing_id, f) in last:
                c = last[(r.ing_id, f)]
                prov_rows.append(dict(species_id=r.ing_id, field=f, value=val, source_id=f"rule:{c.rule_set}", rule_id=c.rule_id, step=c.step))
            else:
                if f == "copies":
                    src = (f"v0:{ad.get('proteome')}[{V0.source.get(r.genes, '')}]" if mono and r.genes in V0.index
                           else f"assembly:{ad.get('assembly') or 'none'}+ingredients:{bname}")
                elif f == "comp_weights":
                    src = base_loc
                elif f in ("transmembrane", "n_tm") and mono:
                    src = f"v0:{ad.get('localisation')}"
                else:
                    src = f"ingredients:{bname}"
                prov_rows.append(dict(species_id=r.ing_id, field=f, value=val, source_id=src, rule_id="", step="ingredients"
                                      if r.ing_id in base_ids else created.get(r.ing_id, "")))
        rules_here = sorted({c.rule_id for (s, f), c in last.items() if s == r.ing_id})
        gs = str(r.genes).split(";")
        sp_rows.append(dict(
            species_id=r.ing_id, species_key=species_key(r, g2u), display_name=r.ing_id, kind=r.kind, genes=r.genes,
            uniprots=";".join(g2u.get(g, "") if isinstance(g2u.get(g), str) else "" for g in gs), stoich=r.stoich,
            n_chains=int(sum(int(float(s)) for s in str(r.stoich).split(";"))), copies=r.copies, mass_da=r.mass_da,
            mass_pg=r.mass_pg, struct_src=r.struct_src, tm=r.transmembrane, n_tm=r.n_tm, signal_frac=r.signal_frac,
            om_frac=r.mito_om_frac,
            assembly_source="" if mono else f"{ad.get('assembly') or 'none'}/{bname}",
            abundance_source=(V0.source.get(r.genes, "") if mono else "subunit pools"),
            loc_sources=";".join([base_loc] + sorted({key_rule[(r.ing_id, k)][0] for k in w if (r.ing_id, k) in key_rule})),
            last_rule_id=rules_here[-1] if rules_here else "",
            build_step="ingredients" if r.ing_id in base_ids else created.get(r.ing_id, "")))
    ctx.write_tsv(pd.DataFrame(sp_rows, columns=schema.SPECIES), "species.tsv")
    ctx.write_tsv(pd.DataFrame(loc_rows, columns=schema.SPECIES_LOC), "species_loc.tsv")
    ctx.write_tsv(pd.DataFrame(prov_rows, columns=schema.FIELD_PROVENANCE), "field_provenance.tsv")

    # ---- sources registry (declared frozen inputs of the inventory + config-only sources)
    plan = steps(cfg)
    used = {}
    for s, d in plan.items():
        for n in d.get("frozen", []):
            used.setdefault(n, []).append(s)
    src_rows = []
    for n, spec in sorted((cfg.get("frozen_inputs") or {}).items()):
        if not n.startswith("inv_"):
            continue
        src_rows.append(dict(source_id=n, kind=spec.get("kind", ""), path=spec.get("path", ""), url=spec.get("url", ""),
                             version=spec.get("version", ""), date=spec.get("date", ""), md5=spec.get("md5", spec.get("sha256", "")),
                             licence=spec.get("licence", ""), redistributable=spec.get("redistributable", ""),
                             used_by=",".join(used.get(n, [])), description=spec.get("why", "")))
    known = {r["source_id"] for r in src_rows}
    for lst in ("proteome_sources", "localisation_sources"):
        for s in inv.get(lst, []) or []:
            if s["source_id"] not in known:
                src_rows.append(dict(source_id=s["source_id"], kind=lst.split("_")[0], path=s.get("path", ""), url=s.get("url", ""),
                                     version=s.get("version", ""), date="", md5="", licence=s.get("licence", ""),
                                     redistributable=s.get("redistributable", ""), used_by="v0", description=s.get("notes", "")))
                known.add(s["source_id"])
    ctx.write_tsv(pd.DataFrame(src_rows, columns=schema.SOURCES), "sources.tsv")
    ps = pd.DataFrame(inv.get("proteome_sources", []) or [], columns=schema.PROTEOME_SOURCES)
    ps["total_protein_pg"] = inv.get("total_protein_pg"); ps["total_protein_source"] = inv.get("total_protein_source", "")
    ctx.write_tsv(ps[schema.PROTEOME_SOURCES], "proteome_sources.tsv")
    ctx.write_tsv(pd.DataFrame(inv.get("localisation_sources", []) or [], columns=schema.LOCALISATION_SOURCES)[schema.LOCALISATION_SOURCES],
                  "localisation_sources.tsv")

    # ---- per-gene abundance (long), from the v0 columns named in inventory.canonical.abundance_columns
    ab = []
    for spec in can.get("abundance_columns", []) or []:
        c = spec["column"]
        if c not in V0.columns:
            continue
        sub = V0[V0[c].notna()]
        used_mask = sub.source.isin(spec.get("used_when_source", [])) if spec.get("used_when_source") else pd.Series(True, sub.index)
        ab.append(pd.DataFrame(dict(gene=sub.index, uniprot=sub.uniprot.values, source_id=spec["source_id"], value=sub[c].values,
                                    unit=spec.get("unit", ""), n_rep=sub[spec["n_rep_column"]].values if spec.get("n_rep_column") else "",
                                    used=used_mask.values)))
    ctx.write_tsv(pd.concat(ab, ignore_index=True)[schema.GENE_ABUNDANCE] if ab else schema.empty(schema.GENE_ABUNDANCE),
                  "gene_abundance.tsv")

    # ---- SEC-MS (per protein / per complex) and surface proteomics: adapter outputs or empty tables with the columns
    # SEC-MS: the secms step's tables (toolkits/secms) when inventory.adapters.secms is set, else an assembly adapter's
    # assembly/secms_*.tsv (older adapters), else empty tables with the columns (spec §4a)
    from . import secms as _secms
    for out, cols in _secms.canonical_tables().items():
        src = [ctx.path(f"secms/{out}")] if ad.get("secms") else [ctx.path(f"assembly/{out}")]
        T = pd.read_csv(src[0], sep="\t") if os.path.exists(src[0]) else schema.empty(cols)
        ctx.write_tsv(T.reindex(columns=cols), out)
    if ad.get("surface"):
        from . import surface as _surface  # noqa: F401  (registers adapter.surface surface_proteomics)
        S = registry.get("adapter.surface", ad["surface"])(ctx)
        ctx.write_tsv(S.reindex(columns=schema.SURFACE), "surface.tsv")
    else:
        ctx.write_tsv(schema.empty(schema.SURFACE), "surface.tsv")

    # ---- version changes (gene level, corrected rule ids)
    tabs = [("ingredients", base)]
    for n in inv.get("rule_sets", []) or []:
        tabs.append((f"rules:{n}", pd.read_csv(ctx.need(f"rules/{n}/ingredients.tsv"), sep="\t")))
    VC = _version_changes(tabs, ch, key_comp)
    ctx.write_tsv(VC, "version_changes.tsv")
    # ---- per-protein SEC-MS bound fraction vs the modelled assembly state (secms_bound.py; spec §4c). Written last and
    # read-only on every input: the tables above are unchanged by it.
    from . import secms_bound
    V0g = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t", dtype={"gene": str}).drop_duplicates("gene").set_index("gene")
    BF, bsum = secms_bound.compute(ctx, P, V0g)
    ctx.write_tsv(BF, "secms_bound_fraction.tsv")
    p = ctx.out("secms_bound_fraction_summary.json")
    from .. import build as _b
    _b.atomic_save(p, lambda t: open(t, "w").write(json.dumps(bsum, indent=1, default=float) + "\n"))
    ctx.say(f"canonical secms_bound_fraction: {bsum['n_measured']}/{bsum['n_proteins']} proteins measured "
            f"({bsum['secms_covered_frac']:.3f} of the mass); modelled assembled {bsum['modelled_assembled_frac']:.3f}, "
            f"measured bound but unrepresented {bsum['unrepresented_bound_frac']:.3f} of the protein mass")
    ctx.say(f"canonical: {len(sp_rows)} species, {len(loc_rows)} species_loc rows, {len(prov_rows)} provenance rows, "
            f"{len(src_rows)} sources, {len(VC)} version-change rows")
