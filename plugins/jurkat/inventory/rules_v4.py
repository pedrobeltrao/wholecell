# origin: jurkat/yc/pipeline/jinventory_v4.py (legacy workspace, 2026-11-15)
# changes: inventory_rule "jurkat.v4": reads the previous step's ingredients.tsv (no in-place rewrite, no _pre_v4 snapshot
# logic: always "apply"), returns the new table + the legacy change table (changes_legacy.tsv, byte-compatible) + the
# canonical changes (rule_id = v4.<letter>). Curated lists (CURATED_W, DROP_MITO, SOLUBLE_NOT_TM, OBLIGATE, KEEP_HET,
# NPC entity map, CD3 cap, PM curations, regexes) in plugins/jurkat/data/rules_v4.yaml; inputs are frozen (UniProt dump,
# topology, 7R5K mmCIF, membrane_lib_ids, library folder for the h-note). Numerics identical to legacy.
"""Inventory revision v4: curated overrides (single-pass type I/II TM -> PM, free CD3 chains, curated PM routing, NPC
stoichiometry from 7R5K, CPX-4747 cap, curated localisation fixes + KDEL, obligate OXPHOS complexes, cross-compartment AFDB
heterodimers dissolved). See the legacy docstring (origin file) for the full rule text."""
import json
import os
import re

import numpy as np
import pandas as pd

from wholecell.inventory import localisation as LOC
from wholecell.inventory.localisation import overrides as OV
from wholecell.registry import register
from .common import data, normalise_changes

# the step order of this rule set = the localisation preset's rules."jurkat.v4" list (checked at run time)
STEPS = [("v4.a", "legacy"), ("v4.b", "legacy"), ("v4.c", "overrides"), ("v4.d_e", "legacy"), ("v4.f_curated", "overrides"),
         ("v4.f_drop_mito", "overrides"), ("v4.f_clic", "overrides"), ("v4.f_kdel", "legacy"), ("v4.g", "legacy"),
         ("v4.h", "legacy")]

AV = 6.022e23


def wjson(w):
    return json.dumps({k: round(float(v), 3) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 0.001})


def pg(copies, mass):
    return copies * mass / AV * 1e12


def npc_7r5k(cif, entity_gene):
    import gemmi
    st = gemmi.read_structure(cif)
    blk = gemmi.cif.read(cif).sole_block()
    ent = {r[0]: (gemmi.cif.as_string(r[1]), float(r[2])) for r in blk.find("_entity.", ["id", "pdbx_description", "formula_weight"])}
    st.setup_entities()
    n = {}
    for ch in st[0]:
        sub = ch.get_polymer()
        if not len(sub): continue
        eid = st.get_entity_of(sub).name
        n[eid] = n.get(eid, 0) + 1
    rows = [(eid, ent[eid][0], entity_gene[ent[eid][0]], c, ent[eid][1]) for eid, c in n.items()]
    return pd.DataFrame(rows, columns=["entity", "description", "gene", "chains", "entity_fw"])


def apply_v4(ctx, src):
    p = ctx.params("jurkat.v4")
    K = data("rules_v4.yaml")
    NPC_COPIES, CD3_CAP, N_TERM_TM, MAJOR_F = K["npc_copies"], K["cd3_cap"], K["n_term_tm"], K["major_subunit_mass_fraction"]
    NPC_COPIES = int(p.get("npc_copies", NPC_COPIES))     # B4: the configured pore count (config params jurkat.v4.npc_copies)
    PM_CC = re.compile(K["pm_cc_regex"], re.I); PM_GO = re.compile(K["pm_go_regex"], re.I)
    NO_LOC = set(K["no_loc_keys"])
    OBLIGATE, KEEP_HET = K["obligate"], set(K["keep_het"])
    SOL, TMR = K["report_sol"], K["report_tmr"]
    LIB = ctx.frozen(p["struct_library"])
    CD3 = list(K["cd3_free_chains"])
    LOC.check_rule_steps(ctx, "jurkat.v4", STEPS)
    OV_C = LOC.override_rows(ctx, "jurkat.v4", "v4.c")                     # curated T-cell surface proteins
    OV_F = LOC.override_rows(ctx, "jurkat.v4", ["v4.f_curated", "v4.f_drop_mito", "v4.f_clic"])

    def region_mass(T, memlib):
        out = {}
        for r in T.itertuples():
            for comp, wt in json.loads(r.comp_weights).items():
                if r.transmembrane and r.ing_id in memlib:
                    reg = TMR.get(comp, "ER")
                    if reg == "MITO_IM" and r.mito_om_frac >= 0.5: reg = "MITO_OM"
                else:
                    reg = SOL.get(comp, "cytosol")
                    if reg in ("ER_lumen", "golgi_lumen", "endo_lumen") and r.signal_frac < 0.5: reg = "cytosol"
                out[reg] = out.get(reg, 0.0) + wt * r.mass_pg
        return out

    I0 = pd.read_csv(src, sep="\t")
    I = I0.copy()
    inv = pd.read_csv(ctx.need("assembly/inventory_v0_complexed.tsv"), sep="\t", index_col="gene")
    up = pd.read_csv(ctx.frozen(p["uniprot_full"]), sep="\t", usecols=["Entry", "Subcellular location [CC]", "Gene Ontology (cellular component)",
                                                                       "Signal peptide", "Length", "Sequence"]).set_index("Entry")
    topo = pd.read_csv(ctx.frozen(p["uniprot_topology"]), sep="\t").set_index("Entry")
    memlib = set(open(ctx.frozen(p["membrane_lib_ids"])).read().split())
    log = []
    say = ctx.say
    say(f"source table: {src}  ({len(I0)} rows, {I0.mass_pg.sum():.3f} pg)")

    def chg(i, field, new, reason):
        old = I.at[i, field]
        if str(old) == str(new): return
        log.append(dict(ing_id=I.at[i, "ing_id"], gene=I.at[i, "genes"], field=field, old=old, new=new, reason=reason))
        I.at[i, field] = new

    def row(ing):
        m = np.nonzero((I.ing_id == ing).values)[0]
        return int(m[0]) if len(m) else None

    # a. single-pass type I/II -> PM
    single = I.genes.map(lambda s: len(set(s.split(";"))) == 1) & I.transmembrane
    gene_up = inv.uniprot
    cand = []
    skip_a = set(CD3) | set(OV_C.key)
    for i in np.nonzero(single.values)[0]:
        r = I.iloc[i]; g = r.genes.split(";")[0]; w = json.loads(r.comp_weights)
        if g in skip_a: continue
        if not (w.get("unassigned", 0) >= 0.5 or set(w) <= NO_LOC): continue
        acc = r.uniprot if isinstance(r.uniprot, str) and r.uniprot else gene_up.get(g)
        if not isinstance(acc, str) or acc not in up.index: continue
        tmtxt = topo.Transmembrane.get(acc, np.nan)
        tms = re.findall(r"TRANSMEM (\d+)\.\.(\d+)", tmtxt) if isinstance(tmtxt, str) else []
        if len(tms) != 1: continue
        cc = str(up.at[acc, "Subcellular location [CC]"]); go = str(up.at[acc, "Gene Ontology (cellular component)"])
        sig = isinstance(up.at[acc, "Signal peptide"], str)
        tm0 = int(tms[0][0]); L = int(up.at[acc, "Length"])
        cc1, cc2 = "type I membrane" in cc, re.search(r"type II membrane|signal-anchor", cc, re.I) is not None
        if re.search(r"type III membrane|type IV membrane|tail-anchored", cc, re.I) and not (cc1 or cc2): continue
        t1 = cc1 or (sig and not cc2)
        t2 = cc2 or (not sig and tm0 <= N_TERM_TM)
        if not (t1 or t2): continue
        ev_cc = PM_CC.search(re.sub(r"\{[^}]*\}", "", cc)) is not None
        ev_go = PM_GO.search(go) is not None
        if not (ev_cc or ev_go): continue
        mv = sum(v for k, v in w.items() if k in NO_LOC)
        nw = {k: v for k, v in w.items() if k not in NO_LOC}
        nw["plasma_membrane"] = nw.get("plasma_membrane", 0) + mv
        why = (f"a: single-pass type {'I' if t1 else 'II'} (sig={int(sig)}, TM@{tm0}/{L}); PM support "
               f"{'+'.join(x for x, f in (('UniProtCC', ev_cc), ('GO_CC', ev_go)) if f)}; was {'unassigned' if 'unassigned' in w else 'cytosol/cytoskeleton'} -> ER default")
        cand.append(r.ing_id)
        chg(i, "comp_weights", wjson(nw), why)
    say(f"[a] single-pass type I/II TM proteins routed to PM: {len(cand)} ingredients")

    # b. free CD3 chains
    tcr = row(K["tcr_complex"]); ntcr = int(I.at[tcr, "copies"]); cap = int(round(CD3_CAP * ntcr))
    for g in CD3:
        i = row(g)
        if i is None: continue
        c0 = int(I.at[i, "copies"]); c1 = min(c0, cap)
        chg(i, "copies", c1, f"b: free {g} capped at {CD3_CAP} x TCR-CD3 ({ntcr}); unassembled chains are ER-retained/degraded; surplus removed")
        chg(i, "comp_weights", wjson({"er": 1.0}), "b: unassembled CD3 chains are ER-retained")

    def rows_of(key):
        if key.startswith("CPX-"): return [row(key)] if row(key) is not None else []
        return list(np.nonzero(I.genes.map(lambda x: set(x.split(";")) == {key}).values)[0])

    def find(kt, key):
        return ([row(key)] if row(key) is not None else []) if kt == "species" else rows_of(key)

    def get(i, field):
        return json.loads(I.at[i, field]) if field == "comp_weights" else I.at[i, field]

    # c. curated PM proteins (override rows v4.c; op set_keep_other on the rule-set input weights)
    def put_c(i, r, nw):
        chg(i, "comp_weights", wjson(nw), f"c: {r.key} is a {r.reason}; HPA/UniProt localisation was {I0.at[i, 'comp_weights']}")
    OV.apply(OV_C, find, lambda i, f: json.loads(I0.at[i, f]), put_c)

    # d. NPC stoichiometry (7R5K) + e. CPX-4747
    E = npc_7r5k(ctx.frozen(p["npc_cif"]), K["npc_entity_gene"])
    gm = inv.mass_da
    st = E.groupby("gene").chains.sum().to_dict()
    e98 = E[E.gene == "NUP98"]
    st["NUP98"] = int(round((e98.chains * e98.entity_fw).sum() / gm["NUP98"]))
    npc = row("NPC")
    old_st = dict(zip(I.at[npc, "genes"].split(";"), map(int, I.at[npc, "stoich"].split(";"))))
    genes = sorted(st, key=lambda g: -st[g] * gm[g])
    new_mass = float(sum(st[g] * gm[g] for g in genes))
    chg(npc, "genes", ";".join(genes), "d: NPC composition = PDB 7R5K (the placed library); basket/FG Nups not in 7R5K dropped")
    chg(npc, "stoich", ";".join(str(st[g]) for g in genes), "d: per-gene chain counts in 7R5K assembly 1 (NUP98 = Nup98+Nup96 precursor equivalents)")
    chg(npc, "mass_da", new_mass, "d: sum of stoich x full-length gene mass")
    upsig = up["Signal peptide"].notna(); uptm = topo.Transmembrane.notna()
    acc = {g: gm_acc for g, gm_acc in inv.uniprot.items()}
    sf = round(sum(st[g] * gm[g] for g in genes if upsig.get(acc.get(g), False)) / new_mass, 3)
    chg(npc, "signal_frac", sf, "d: recomputed for the 7R5K composition")
    chg(npc, "n_tm", int(sum(bool(uptm.get(acc.get(g), False)) for g in genes)), "d: recomputed for the 7R5K composition")
    chg(npc, "mass_pg", pg(I.at[npc, "copies"], new_mass), "d: mass_da changed")
    cx = row(K["ranbp2_complex"])
    if cx is not None and "RANBP2" in I.at[cx, "genes"].split(";"):
        avail = int(inv.copies["RANBP2"]) - NPC_COPIES * st.get("RANBP2", 0)
        c0 = int(I.at[cx, "copies"]); c1 = max(0, min(c0, avail))
        chg(cx, "copies", c1, f"e: RANBP2 inventory {int(inv.copies['RANBP2']):,} - {NPC_COPIES * st.get('RANBP2', 0):,} in NPCs; complex was overdrawn")
        chg(cx, "mass_pg", pg(c1, I.at[cx, "mass_da"]), "e: copies changed")
    touched = set(old_st) | set(st) | (set(I.at[cx, "genes"].split(";")) if cx is not None else set())
    other = {g: 0.0 for g in touched}
    for i, r in I.iterrows():
        if r.kind == "protein" or r.ing_id == "NUCLEOSOME": continue
        for g, s in zip(r.genes.split(";"), r.stoich.split(";")):
            if g in other: other[g] += r.copies * float(s)
    for g in sorted(touched, key=lambda g: -st.get(g, 0) * gm.get(g, 0)):
        tot = int(inv.copies.get(g, 0))
        i = row(g)
        f1 = int(max(0, round(tot - other[g])))
        if i is not None:
            chg(i, "copies", f1, f"d/e: free pool = inventory {tot:,} - assemblies {int(other[g]):,} (NPC 7R5K stoich {st.get(g, 0)}, was {old_st.get(g, 0)})")
            chg(i, "mass_pg", pg(f1, I.at[i, "mass_da"]), "d/e: copies changed")

    # f. curated localisation fixes (top-200 review)
    def setw(i, nw, why, tag):
        w0 = json.loads(I.at[i, "comp_weights"]); nw = {k: v / sum(nw.values()) for k, v in nw.items()}
        moved = 0.5 * sum(abs(nw.get(k, 0) - w0.get(k, 0)) for k in set(nw) | set(w0)) * I.at[i, "mass_pg"]
        if moved < 1e-12: return
        chg(i, "comp_weights", wjson(nw), f"f/{tag}: {why}")

    # override rows v4.f_curated / v4.f_drop_mito / v4.f_clic (comp_weights normalised by setw; fields via chg)
    def put_f(i, r, new):
        tag = r.rule_id.split(".f_", 1)[1]
        if r.field == "comp_weights": setw(i, new, r.reason, tag)
        else: chg(i, r.field, new, f"f/{tag}: {r.reason}")
    OV.apply(OV_F, find, get, put_f)
    kdel = re.compile(K["kdel_regex"])
    for i in np.nonzero(((~I.transmembrane) & (I.signal_frac >= 0.5)).values)[0]:
        gs = set(I.at[i, "genes"].split(";"))
        if len(gs) != 1: continue
        a = inv.uniprot.get(next(iter(gs)))
        sq = up.Sequence.get(a) if isinstance(a, str) else None
        if not isinstance(sq, str) or not kdel.search(sq): continue
        w = json.loads(I.at[i, "comp_weights"])
        if set(w) <= {"er"}: continue
        setw(i, {"er": 1.0}, f"ER-retrieval C-terminus {sq[-4:]} + signal peptide: ER lumen resident (non-lumen weights pack it in the cytosol)", "kdel")

    # g. obligate OXPHOS complexes absorb free core subunits
    nasm = {}
    for r in I[I.kind != "protein"].itertuples():
        for g in set(r.genes.split(";")): nasm[g] = nasm.get(g, 0) + 1
    gm_ = inv.mass_da
    for cid in OBLIGATE:
        ci = row(cid)
        if ci is None: continue
        gs = I.at[ci, "genes"].split(";"); ss = [float(x) for x in I.at[ci, "stoich"].split(";")]; N = int(I.at[ci, "copies"])
        fr = {g: (int(I.at[row(g), "copies"]) if row(g) is not None else 0) for g in gs}
        ded = [(g, s) for g, s in zip(gs, ss) if nasm.get(g, 0) == 1]
        if len(ded) < 2: continue
        imp = np.array([N + fr[g] / s for g, s in ded]); w = np.array([gm_[g] * s for g, s in ded])
        o = np.argsort(imp); Nmed = imp[o][np.searchsorted(np.cumsum(w[o]) / w.sum(), 0.5)]
        tot_m = sum(gm_[g] * s for g, s in zip(gs, ss))
        cap = min((N + fr[g] / s for g, s in zip(gs, ss) if gm_[g] * s >= MAJOR_F * tot_m), default=None)
        if cap is None: continue
        Ns = int(round(min(Nmed, cap)))
        if Ns <= N: continue
        chg(ci, "copies", Ns, f"g: obligate complex; N = min(subunit-mass-weighted median, lowest major subunit) of (N + free/stoich) over {len(ded)} dedicated subunits (was {N:,})")
        for g, s in zip(gs, ss):
            need = (Ns - N) * s; take = min(fr[g], need)
            if row(g) is not None and take > 0:
                chg(row(g), "copies", int(round(fr[g] - take)), f"g: free subunit absorbed into {cid} ({N:,} -> {Ns:,} copies)")
        say(f"[g] {cid}: {N:,} -> {Ns:,} copies")

    # h. dissolve cross-compartment AFDB heterodimers
    sigs = up["Signal peptide"].notna()

    def cls(g):
        if bool(topo.Transmembrane.notna().get(inv.uniprot.get(g), False)): return None
        if bool(sigs.get(inv.uniprot.get(g), False)): return "lumen"
        i = row(g)
        w = json.loads(I.at[i, "comp_weights"]) if i is not None and I.at[i, "kind"] == "protein" else {
            k: float(v) for k, v in (t.split(":") for t in str(inv.compartments.get(g, "cytosol:1")).split(";"))}
        tot = sum(w.values()) or 1.0
        mi = w.get("mitochondrion", 0) / tot
        cy = sum(v for k, v in w.items() if k in ("cytosol", "nucleus", "cytoskeleton", "unassigned", "plasma_membrane")) / tot
        return "mito" if mi > 0.5 else "cyto" if cy > 0.5 else None
    BAD = {frozenset(("lumen", "cyto")), frozenset(("lumen", "mito")), frozenset(("mito", "cyto"))}
    drop = []
    for i in np.nonzero((I.kind == "afdb_heterodimer").values)[0]:
        g1, g2 = I.at[i, "genes"].split(";")
        c1, c2 = cls(g1), cls(g2)
        if None in (c1, c2) or frozenset((c1, c2)) not in BAD or I.at[i, "ing_id"] in KEEP_HET: continue
        n = int(I.at[i, "copies"]); note = []
        for g, c in ((g1, c1), (g2, c2)):
            j = row(g)
            if j is not None and I.at[j, "kind"] == "protein":
                chg(j, "copies", int(I.at[j, "copies"]) + n, f"h: + {n:,} copies from dissolved heterodimer {I.at[i, 'ing_id']} ({c1} vs {c2})")
                if not os.path.exists(f"{LIB}/{g}.npz"): note.append(f"{g}: monomer row but no library")
            else:
                note.append(f"{g}: no monomer row, {n:,} copies dropped")
        log.append(dict(ing_id=I.at[i, "ing_id"], gene=I.at[i, "genes"], field="row", old=f"copies={n}", new="deleted",
                        reason=f"h: dissolved, partners in incompatible compartments ({g1}={c1}, {g2}={c2})" + (f"; {'; '.join(note)}" if note else "")))
        drop.append(i)
    say(f"[h] {len(drop)} heterodimers dissolved")

    # bookkeeping
    I["copies"] = I.copies.astype(np.int64)
    for i in set(np.nonzero((I.copies.values != I0.copies.values) | (I.mass_da.values != I0.mass_da.values))[0]):
        I.at[i, "mass_pg"] = pg(I.at[i, "copies"], I.at[i, "mass_da"])
    I = I.drop(index=I.index[drop]).reset_index(drop=True)
    rb, ra = region_mass(I0, memlib), region_mass(I, memlib)
    say("expected mass by pack region (pg, before -> after): " + "; ".join(
        f"{k} {rb.get(k, 0):.3f} -> {ra.get(k, 0):.3f}" for k in sorted(set(rb) | set(ra), key=lambda k: -ra.get(k, 0))))
    L = pd.DataFrame(log, columns=["ing_id", "gene", "field", "old", "new", "reason"])
    say(f"total mass {I0.mass_pg.sum():.4f} -> {I.mass_pg.sum():.4f} pg; {len(L)} field changes in {L.ing_id.nunique()} ingredients")
    for c in I0.columns:
        I[c] = I[c].astype(I0[c].dtype)
    assert list(I.columns) == list(I0.columns) and len(I) == len(I0) - len(drop)
    return {"ingredients": I, "changes": normalise_changes(L, "v4"), "extra": {"changes_legacy.tsv": L}}


apply_v4.OUTPUTS = ["changes_legacy.tsv"]
apply_v4.PREV = ["assembly/inventory_v0_complexed.tsv"]
register("inventory_rule", "jurkat.v4")(apply_v4)
