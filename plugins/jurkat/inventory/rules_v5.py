# origin: jurkat/yc/pipeline/jinventory_v5.py (legacy workspace, 2026-11-15)
# changes: inventory_rule "jurkat.v5": reads the previous step's table; the v4-curated gene set comes from the v4 step's
# changes_legacy.tsv of this build (was data/derived/inventory_v4_changes.tsv); pack-dependent inputs are declared frozen
# tables (v4-pack loc_compare.tsv, tm_route_libs.tsv, hpa_celllines.tsv); loc_compare.sources is the minimal port in
# locsrc.py; constants (VN/VC/VNUCL, NUCL_MAX, LAM, ...) in plugins/jurkat/data/rules_v5.yaml. Numerics identical.
"""Inventory revision v5: placement classes from the Tier-1 localisation review. Rules in order: T (TM routing for
ingredients with a membrane library), L (lamina), A (lipid anchors -> anchored shells), C (cortex), N (nucleolus, capped),
P (nuclear pool of small soluble proteins, v4-pack placed shares)."""
import json
import re

import numpy as np
import pandas as pd

from wholecell.inventory import localisation as LOC
from wholecell.inventory.localisation import overrides as OV
from wholecell.registry import register
from . import locsrc
from .common import capacity_scale, data, geometry_volumes, normalise_changes

# the step order of this rule set = the localisation preset's rules."jurkat.v5" list (checked at run time)
STEPS = [("v5.T_tm_routing", "legacy"), ("v5.L_lamina", "legacy"), ("v5.A_lipid", "legacy"), ("v5.C_cortex", "legacy"),
         ("v5.N_nucleolus", "legacy"), ("v5.P_nuclear_pool", ("frozen_overrides", "rule"))]   # legacy_v7 | common (B5)

AV = 6.02214076e23


def wjson(w):
    return json.dumps({k: round(float(v), 4) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 1e-4})


def curated_v4(ctx):
    ch = pd.read_csv(ctx.need("rules/jurkat.v4/changes_legacy.tsv"), sep="\t")
    cur = ch[ch.field.isin(["comp_weights", "transmembrane"]) & ch.reason.astype(str).str.match(r'^"?[bcf]')]
    return cur, {g for x in cur.gene.astype(str) for g in x.split(";")}


def apply_v5(ctx, src):
    p = ctx.params("jurkat.v5")
    K = data("rules_v5.yaml")
    VN, VC, VNUCL = K["volumes_um3"]["nucleoplasm"], K["volumes_um3"]["cytosol"], K["volumes_um3"]["nucleolus"]
    NUCL_MAX, MAX_ANCHOR_DA = K["nucleolus_cap_mg_ml"], K["max_anchor_da"]
    CYTO_KEYS = set(K["cyto_keys"]); ACT_TUB = re.compile(K["act_tub_regex"])
    LOC.check_rule_steps(ctx, "jurkat.v5", STEPS)
    I0 = pd.read_csv(src, sep="\t"); I = I0.copy()
    T = pd.read_csv(ctx.frozen(p["v4pack_loc_compare"]), sep="\t", index_col=0, low_memory=False)
    UPP = ctx.frozen(p["uniprot_full"])
    up = pd.read_csv(UPP, sep="\t", usecols=["Entry", "Gene Names (primary)", "Subcellular location [CC]", "Lipidation", "Transmembrane",
                                             "Intramembrane"]).set_index("Entry")
    C = pd.read_csv(ctx.frozen(p["hpa_celllines"]), sep="\t", keep_default_na=False)
    JUR = set(C[C.cell_line.str.upper() == K["hpa_jurkat_line"]].gene)
    _, CUR = curated_v4(ctx)
    inv = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t").drop_duplicates("gene").set_index("gene")
    gmass = inv.mass_da.to_dict(); gacc = inv.uniprot.to_dict()
    log = []
    say = ctx.say
    say(f"source table: {src} ({len(I0)} rows, {I0.mass_pg.sum():.3f} pg); HPA-JURKAT genes {len(JUR)}; v4-curated {len(CUR)}")

    def chg(i, field, new, rule):
        old = I.at[i, field]
        if str(old) == str(new): return
        log.append(dict(ing_id=I.at[i, "ing_id"], gene=I.at[i, "genes"], field=field, old=old, new=new, rule=rule))
        I.at[i, field] = new

    def W(i): return json.loads(I.at[i, "comp_weights"]) if isinstance(I.at[i, "comp_weights"], str) else {}
    def genes(i): return str(I.at[i, "genes"]).split(";")
    def stoich(i): return [float(x) for x in str(I.at[i, "stoich"]).split(";")]

    def frac_mass(i, S):
        m = [(g in S) * gmass.get(g, 0) * s for g, s in zip(genes(i), stoich(i))]; t = sum(gmass.get(g, 0) * s for g, s in zip(genes(i), stoich(i)))
        return sum(m) / t if t > 0 else 0.0

    def single(i): return len(set(genes(i))) == 1
    moved = {}

    def setw(i, nw, rule):
        ow = W(i)
        d = 0.5 * sum(abs(nw.get(k, 0) - ow.get(k, 0)) for k in set(nw) | set(ow))
        if d < 1e-6: return
        r = moved.setdefault(rule, [0.0, set(), set()])
        r[0] += d * I.at[i, "mass_pg"]; r[1] |= set(genes(i)); r[2].add(I.at[i, "ing_id"])
        chg(i, "comp_weights", wjson(nw), rule)

    # T. TM routing
    R = pd.read_csv(ctx.frozen(p["tm_route_libs"]), sep="\t")
    ok = set(R[R.method.astype(str) != "skipped"].ing_id)
    say(f"T: membrane libraries built for {len(ok)} of {len(R)} ingredients")
    for i in np.nonzero(I.ing_id.isin(ok).values)[0]:
        chg(i, "transmembrane", True, "T_tm_routing")
        if int(I.at[i, "n_tm"]) == 0:
            n = 0
            for g in set(genes(i)):
                a = gacc.get(g)
                if isinstance(a, str) and a in up.index:
                    n += sum(str(up.at[a, c]).count(k) for c, k in (("Transmembrane", "TRANSMEM"), ("Intramembrane", "INTRAMEM")))
            chg(i, "n_tm", max(n, 1), "T_tm_routing")
        w = W(i)
        if w.get("mitochondrion", 0) > 0:
            def _om(t): t = t.lower(); return "outer membrane" in t and "inner membrane" not in t and "matrix" not in t
            om = [_om(str(up.at[gacc[g], "Subcellular location [CC]"])) for g in set(genes(i))
                  if isinstance(gacc.get(g), str) and gacc[g] in up.index and "mitochondri" in str(up.at[gacc[g], "Subcellular location [CC]"]).lower()]
            if om: chg(i, "mito_om_frac", 1.0 if np.mean(om) >= 0.5 else 0.0, "T_tm_routing_OM")

    # L. lamina
    LAM = dict(K["lamina_share"])
    for i in range(len(I)):
        g = genes(i)
        if not single(i) or g[0] not in LAM or bool(I.at[i, "transmembrane"]): continue
        w = W(i); f = LAM[g[0]]; pool = w.get("nucleus", 0) + w.get("nuclear_envelope", 0)
        nw = {k: v for k, v in w.items() if k not in ("nucleus", "nuclear_envelope")}
        nw["lamina"] = f * pool
        if 1 - f > 0: nw["nucleus"] = (1 - f) * pool
        setw(i, nw, "L_lamina")

    # A. lipid anchors
    def lipid_kind(g):
        a = gacc.get(g)
        if not isinstance(a, str) or a not in up.index: return None
        t = str(up.at[a, "Lipidation"])
        if "GPI-anchor" in t: return "gpi"
        if re.search(r"myristoyl|palmitoyl|farnesyl|geranylgeranyl", t): return "cyto"
        return None
    LK = {g: lipid_kind(g) for g in inv.index}
    LK = {g: k for g, k in LK.items() if k}
    MAPA = dict(K["anchor_membranes"])
    lam_genes = set(LAM)
    for i in range(len(I)):
        if bool(I.at[i, "transmembrane"]) or I.at[i, "mass_da"] > MAX_ANCHOR_DA: continue
        ks = {LK[g] for g in set(genes(i)) if g in LK}
        if not ks or (single(i) and genes(i)[0] in lam_genes): continue
        side = "gpi" if ks == {"gpi"} else "cyto"
        w = W(i); nw = dict(w); hit = False
        for comp, pre in MAPA.items():
            if w.get(comp, 0) <= 0: continue
            if pre == "ne":
                if side == "gpi": continue
                key = "ne_anchor_nuc"
            else:
                key = f"{pre}_anchor_" + ({"pm": "exo"}.get(pre, "lumen") if side == "gpi" else "cyto")
            nw[key] = nw.get(key, 0) + nw.pop(comp); hit = True
        if w.get("mitochondrion", 0) > 0 and side == "cyto":
            om = any("mitochondrion outer membrane" in str(up.at[gacc[g], "Subcellular location [CC]"]).lower() for g in set(genes(i))
                     if g in LK and isinstance(gacc.get(g), str) and gacc[g] in up.index)
            if om: nw["om_anchor_cyto"] = nw.pop("mitochondrion"); hit = True
        if hit: setw(i, nw, f"A_lipid_{side}")

    # C. cortex (R6 selection as loc_compare.rules)
    cand = T.index[T.cause.isin(K["cortex_causes"]) & T.actin_cortex_annot.astype(bool) & (T.m_PM < 0.05)]
    HPAP = ctx.frozen(p["hpa_subcellular"])
    S = locsrc.sources(T.loc[cand], HPAP, UPP)
    R6 = set()
    for g in cand:
        if g in CUR: continue
        srcs = ("uniprot",) if g in JUR else ("hpa", "uniprot")
        if any(S[g].get(k, {}).get("q", {}).get("PM", 0) >= 0.999 and S[g].get(k, {}).get("w", 0) > 0 for k in srcs): R6.add(g)
    say(f"C: cortical genes {len(R6)} (of {len(cand)} R6 candidates)")
    for i in range(len(I)):
        if bool(I.at[i, "transmembrane"]): continue
        gs = set(genes(i)) & R6
        if not gs or not (single(i) or frac_mass(i, R6) >= 0.5): continue
        w = W(i); nw = {}; cx = 0.0
        for k, v in w.items():
            if k == "plasma_membrane": cx += v
            elif k in CYTO_KEYS: cx += 0.5 * v; nw[k] = 0.5 * v
            else: nw[k] = v
        if cx <= 0: continue
        nw["cortex"] = nw.get("cortex", 0) + cx
        setw(i, nw, "C_cortex")

    # N. nucleolus
    NF = K["nucleolar_share"]; OTH = K["nucleolar_other_nuclear_terms"]
    hm = T.hpa_main.fillna(""); rel = T.hpa_reliability.fillna("")
    NU = {}
    for g in T.index[hm.str.contains("Nucleoli") & rel.isin(["Enhanced", "Supported", "Approved"])]:
        if g in JUR or g in CUR: continue
        terms = [t.strip() for t in hm[g].split(";")]
        nuc_other = [t for t in terms if t in OTH]
        NU[g] = NF["enhanced_only_nucleoli"] if (rel[g] in ("Enhanced", "Supported") and not nuc_other) else NF["approved"]
    nset = set(NU); plan = []
    for i in range(len(I)):
        if bool(I.at[i, "transmembrane"]): continue
        gs = set(genes(i)) & nset
        if not gs: continue
        if not single(i) and frac_mass(i, nset) < 0.5: continue
        w = W(i)
        if w.get("nucleus", 0) <= 0: continue
        f = np.average([NU[g] for g in gs])
        plan.append((i, f, w))
    m_nucl = sum(f * w["nucleus"] * I.at[i, "mass_pg"] for i, f, w in plan)
    dens = m_nucl / VNUCL * 1e3
    sc = min(1.0, NUCL_MAX / dens) if dens > 0 else 1.0
    say(f"N: nucleolar genes {len(NU)}; nucleolar protein {m_nucl:.3f} pg = {dens:.0f} mg/ml before the cap; scale {sc:.3f}")
    for i, f, w in plan:
        nw = dict(w); x = sc * f * w["nucleus"]; nw["nucleus"] = w["nucleus"] - x; nw["nucleolus"] = nw.get("nucleolus", 0) + x
        setw(i, nw, "N_nucleolus")

    # P. nuclear pool (R2b)
    NP = K["nuclear_pool"]
    if LOC.step_kind(ctx, "jurkat.v5", "v5.P_nuclear_pool") == "rule":
        rowsP = _p_intended(ctx, I, NP, JUR, CUR, ACT_TUB, W, genes, single, say)
        OVP = OV.frame(rowsP)
        if I.ing_id.duplicated().any():
            raise ValueError("v5 P: species ids are not unique; a species-keyed override table is ambiguous")
        OV.apply(OVP, lambda kt, key: [int(np.nonzero((I.ing_id == key).values)[0][0])], lambda i, f: W(i),
                 lambda i, r, nw: setw(i, nw, r.rule_id.split(".", 1)[1]))
        return _finish(I0, I, log, moved, say, OVP)
    share = VN / (VN + VC)
    cand = T[(T.m_nucleus < NP["max_model_nucleus"]) & (T.m_cytosol >= NP["min_model_cytosol"]) & (T.unit_kda.fillna(T.mw_kda) < NP["max_unit_kda"])
             & (T.tm_ing_soluble_frac == 0) & ~(T.uniprot_tm.fillna(False).astype(bool)) & ~(T.uniprot_sp.fillna(False).astype(bool))].index
    S2 = locsrc.sources(T.loc[cand], HPAP, UPP)
    R2 = {g for g in cand if any(S2[g].get(k, {}).get("q", {}).get("nucleus", 0) > 0 for k in ("hpa", "uniprot"))
          and g not in JUR and g not in CUR and not ACT_TUB.match(g)}
    say(f"P: nuclear-pool genes {len(R2)} (R2b candidates {len(cand)}); nuclear share Vn/(Vn+Vc) = {share:.3f}")
    # frozen-input override table (preset legacy_v7 only): species selected from the v4 pack's placed shares; the target
    # shares are computed here and applied by the common override function (rows = overrides_P_nuclear_pool.tsv)
    rowsP = []
    for i in range(len(I)):
        if not single(i) or genes(i)[0] not in R2 or bool(I.at[i, "transmembrane"]) or I.at[i, "mass_da"] >= NP["max_mass_da"]: continue
        w = W(i); pool = w.get("cytosol", 0) + w.get("unassigned", 0) + w.get("nucleus", 0)
        tgt = share * pool
        if pool <= 0 or w.get("nucleus", 0) >= tgt: continue
        take = tgt - w.get("nucleus", 0); cy = w.get("cytosol", 0) + w.get("unassigned", 0)
        nw = dict(w)
        for k in ("cytosol", "unassigned"):
            if k in nw: nw[k] = nw[k] * (1 - take / cy)
        nw["nucleus"] = tgt
        g = genes(i)[0]
        rowsP.append(dict(rule_set="jurkat.v5", rule_id="v5.P_nuclear_pool", key_type="species", key=I.at[i, "ing_id"],
                          field="comp_weights", op="set", value=json.dumps(nw),
                          reason=f"small soluble protein with HPA/UniProt nuclear support placed almost only in the cytosol by "
                                 f"the v4 pack (placed nucleus {T.loc[[g], 'm_nucleus'].iloc[0]:.3f}, cytosol {T.loc[[g], 'm_cytosol'].iloc[0]:.3f}); "
                                 f"nucleus = Vn/(Vn+Vc) {share:.3f} of the cytosol+unassigned+nucleus pool",
                          citation="frozen v4-pack loc_compare (params jurkat.v5.v4pack_loc_compare); PLANNING_NOTES R2b",
                          genes=g))
    OVP = OV.frame(rowsP)
    if I.ing_id.duplicated().any():
        raise ValueError("v5 P: species ids are not unique; a species-keyed override table is ambiguous")
    OV.apply(OVP, lambda kt, key: [int(np.nonzero((I.ing_id == key).values)[0][0])], lambda i, f: W(i),
             lambda i, r, nw: setw(i, nw, r.rule_id.split(".", 1)[1]))

    return _finish(I0, I, log, moved, say, OVP)


def _finish(I0, I, log, moved, say, OVP):
    for k, (m, g, ii) in moved.items():
        say(f"  {k} | {len(ii)} | {len(g)} | {m:.3f}")
    L = pd.DataFrame(log, columns=["ing_id", "gene", "field", "old", "new", "rule"])
    say(f"total mass {I0.mass_pg.sum():.4f} -> {I.mass_pg.sum():.4f} pg; {len(L)} field changes in {L.ing_id.nunique()} ingredients")
    for c in I0.columns:
        if c != "transmembrane": I[c] = I[c].astype(I0[c].dtype)
    I["transmembrane"] = I.transmembrane.astype(bool)
    return {"ingredients": I, "changes": normalise_changes(L, "v5"),
            "extra": {"changes_legacy.tsv": L, "overrides_P_nuclear_pool.tsv": OVP}}


def _p_intended(ctx, I, NP, JUR, CUR, ACT_TUB, W, genes, single, say):
    """v5 P re-expressed (preset common, B5 criteria §2): no pack read. Single-gene soluble species (not TM, mass <
    max_mass_da, signal_frac < 0.5) of genes with nuclear support in this build's evidence table (an HPA or UniProt row
    with key nucleus, reliability not 0), not HPA-JURKAT-imaged, not v4-curated, not actin/tubulin, without an imaging
    N/C calibration at v0 (loc_source 'nc:hpa'), whose INTENDED nucleus share < max_model_nucleus and cytosol +
    unassigned >= min_model_cytosol: nucleus raised to Vn / (Vn + Vc) of the cytosol + unassigned + nucleus pool (this
    build's geometry), scaled by the capacity check. Returns override rows (applied by the common override function)."""
    from wholecell.inventory.localisation import evidence as EV
    VN, VC = geometry_volumes(ctx)
    share = VN / (VN + VC)
    E = EV.read(ctx.need("localisation/evidence.tsv"))
    NS = set(E.gene[E.source.isin(["hpa", "uniprot"]) & (E.compartment == "nucleus") & ~(E.reliability == 0)])
    v0 = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t", usecols=["gene", "loc_source"]).set_index("gene").loc_source
    NCIMG = {g for g, s in v0.items() if isinstance(s, str) and "nc:hpa" in s.split("+")}
    plan = []
    for i in range(len(I)):
        if not single(i) or bool(I.at[i, "transmembrane"]) or I.at[i, "mass_da"] >= NP["max_mass_da"]:
            continue
        if float(I.at[i, "signal_frac"]) >= 0.5:
            continue
        g = genes(i)[0]
        if g not in NS or g in JUR or g in CUR or g in NCIMG or ACT_TUB.match(g):
            continue
        w = W(i)
        tot = sum(w.values())
        n, cy = w.get("nucleus", 0.0), w.get("cytosol", 0.0) + w.get("unassigned", 0.0)
        pool = cy + n
        if tot <= 0 or pool <= 0 or n / tot >= NP["max_model_nucleus"] or cy / tot < NP["min_model_cytosol"]:
            continue
        tgt = share * pool
        if n >= tgt:
            continue
        plan.append((i, g, w, n, cy, tgt))
    sc, info = capacity_scale(I, [(i, tgt - n) for i, g, w, n, cy, tgt in plan], VN, VC)
    say(f"P (common, intended shares): nuclear-pool species {len(plan)}; Vn/(Vn+Vc) = {share:.3f} (Vn {VN:.1f}, Vc {VC:.1f} um3); "
        f"capacity {json.dumps({k: round(v, 3) for k, v in info.items()})}")
    rows = []
    for i, g, w, n, cy, tgt in plan:
        take = sc * (tgt - n)
        nw = dict(w)
        for k in ("cytosol", "unassigned"):
            if k in nw:
                nw[k] = nw[k] * (1 - take / cy)
        nw["nucleus"] = n + take
        rows.append(dict(rule_set="jurkat.v5", rule_id="v5.P_nuclear_pool", key_type="species", key=I.at[i, "ing_id"],
                         field="comp_weights", op="set", value=json.dumps(nw),
                         reason=f"small soluble protein with HPA/UniProt nuclear support and an intended nucleus share "
                                f"{n / sum(w.values()):.3f} (< {NP['max_model_nucleus']}); nucleus = Vn/(Vn+Vc) {share:.3f} of the "
                                f"cytosol+unassigned+nucleus pool x capacity scale {sc:.3f}",
                         citation="B5 re-expression of R2b on intended shares (docs/unification_criteria/B5_localisation.md §2)",
                         genes=g))
    return rows


apply_v5.OUTPUTS = ["changes_legacy.tsv", "overrides_P_nuclear_pool.tsv"]
apply_v5.PREV = ["rules/jurkat.v4/changes_legacy.tsv", "v0/inventory_v0.tsv", "localisation/evidence.tsv"]
register("inventory_rule", "jurkat.v5")(apply_v5)
