# origin: jurkat/yc/pipeline/jinventory_v6.py (legacy workspace, 2026-11-15)
# changes: inventory_rule "jurkat.v6": reads the previous step's table; v4-curated / v5-cortex gene sets from this build's
# v4 / v5 changes_legacy.tsv; pack-dependent inputs are declared frozen tables (v5-pack loc_compare_v5.tsv, tier2_v5.tsv,
# tier2_transfer.json, JURKAT image manifest, hold-out lists); the speckle volume (was read from jgeom_cell1.npz summary)
# is a config value (params speckle_volume_um3; legacy value 1.963 um3 of the v6 geometry). Constants and curated
# vocabularies in plugins/jurkat/data/rules_v6.yaml. Outputs: ingredients.tsv, changes (+ changes_legacy.tsv),
# rule_genes.tsv (= v6/inventory_v6_rule_genes.tsv), calibration.json (= v6/calibration_v6.json). NOT ported: the
# dry-run localisation score (dryrun_v6.json; it read the live wc/libs.json and the v5 pack cache; report-only).
# Numerics identical to legacy.
"""Inventory revision v6: systematic localisation rules from imaging categories plus mechanism features, in order
R4 (TM steady-state pools), R2 (actin binders -> cortex), R3 (peripheral membrane binders), R5 (Golgi binders), R6a
(NPC-proximal shell), R1 (N/C from imaging categories), R6b (speckles, capped). Hold-out: OpenCell never read; held-out
HPA-JURKAT genes use a non-JURKAT per-line consensus."""
import json
import re

import numpy as np
import pandas as pd

from wholecell.inventory import localisation as LOC
from wholecell.inventory.localisation import overrides as OV
from wholecell.registry import register
from .common import capacity_scale, data, geometry_volumes, normalise_changes, read_list
from .rules_v5 import curated_v4

VER = "v6"
# the step order of this rule set = the localisation preset's rules."jurkat.v6" list (checked at run time)
STEPS = [("v6.R4_tm_pool", "legacy"), ("v6.R2_cortex", "legacy"), ("v6.R3_peripheral", "legacy"), ("v6.R5_golgi", "legacy"),
         ("v6.R6a_npc", "legacy"), ("v6.R1_nc", ("frozen_overrides", "rule")), ("v6.R6b_speckle", "legacy")]   # R1: legacy_v7 | common (B5)


def apply_v6(ctx, src):
    p = ctx.params("jurkat.v6")
    K = data("rules_v6.yaml")
    VN, VC = K["volumes_um3"]["nucleoplasm"], K["volumes_um3"]["cytosol"]
    CYTO_KEYS, SP_KEYS, ANCHOR_KEYS = set(K["cyto_keys"]), set(K["sp_keys"]), set(K["anchor_keys"])
    ACT_TUB = re.compile(K["act_tub_regex"])
    NUC_T, SPERM, ORG_T, PM_T = set(K["nuc_terms"]), set(K["sperm_terms"]), set(K["organelle_terms"]), set(K["pm_terms"])
    R2_T = PM_T | {"Actin filaments"}
    MEMB_TGT = dict(K["membrane_targets"])
    R_CAT = dict(K["r_by_category"]); S_DEFAULT = K["default_share"]; MIN_CONTRAST = K["min_rim_contrast_log2"]
    SPK_MAX = K["speckle_cap_mg_ml"]; CL0, CL1 = K["share_clip"]
    MEMB_IPR, ACT_IPR = dict(K["memb_ipr"]), dict(K["act_ipr"])
    ACT_GO, ACT_CC = re.compile(K["act_go_regex"]), re.compile(K["act_cc_regex"])
    PERI_CC = re.compile(K["peri_cc_regex"], re.I)
    GOLGI_CC, GOLGI_BP = re.compile(K["golgi_cc_regex"]), re.compile(K["golgi_bp_regex"])
    NPC_CC, NPC_MF, SR_RE = re.compile(K["npc_cc_regex"]), re.compile(K["npc_mf_regex"]), re.compile(K["sr_regex"])
    CHECK = dict(K["check"])
    IZ = K["itzhak"]

    LOC.check_rule_steps(ctx, "jurkat.v6", STEPS)
    R1_RULE = LOC.step_kind(ctx, "jurkat.v6", "v6.R1_nc") == "rule"     # preset common: R1 on intended shares (B5)
    if R1_RULE:
        VN, VC = geometry_volumes(ctx)
    I0 = pd.read_csv(src, sep="\t"); I = I0.copy()
    T = pd.read_csv(ctx.frozen(p["v5pack_loc_compare"]), sep="\t", index_col=0, low_memory=False); T = T[~T.index.duplicated()]
    H = pd.read_csv(ctx.frozen(p["hpa_subcellular"]), sep="\t").drop_duplicates("Gene name").set_index("Gene name")
    C = pd.read_csv(ctx.frozen(p["hpa_celllines"]), sep="\t", keep_default_na=False)
    U = pd.read_csv(ctx.frozen(p["uniprot_domains"]), sep="\t").set_index("Entry")
    inv = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t").drop_duplicates("gene").set_index("gene")
    gmass = inv.mass_da.to_dict(); gacc = inv.uniprot.to_dict()
    HOLD, FIT = read_list(ctx.frozen(p["holdout_genes"])), read_list(ctx.frozen(p["fit_genes"]))
    cur, CUR = curated_v4(ctx)
    CUR_MITO = {g for x in cur[cur.reason.astype(str).str.contains("drop_mito")].gene.astype(str) for g in x.split(";")}
    c5 = pd.read_csv(ctx.need("rules/jurkat.v5/changes_legacy.tsv"), sep="\t")
    V5CX = {g for x in c5[c5.rule == "C_cortex"].gene.astype(str) for g in x.split(";")}
    log, moved = [], {}
    gene_rows = []
    say = ctx.say

    def hpa_cons(g):
        if g in HOLD:
            sub = C[(C.gene == g) & (C.cell_line.str.upper() != "JURKAT") & (C.ab_reliability != "uncertain") & (C.locations != "")]
            if not len(sub): return set(), set(), "", "non-JURKAT per-line consensus (none)"
            locs = [set(x.lower().split(";")) for x in sub.locations]
            fr = {t: sum(t in l for l in locs) / len(locs) for t in set().union(*locs)}
            fix = lambda t: NORM.get(t.lower(), t)
            main = {fix(t) for t, f in fr.items() if f >= 0.5}; add = {fix(t) for t, f in fr.items() if 0.25 <= f < 0.5}
            rel = H.at[g, "Reliability"] if g in H.index else ""
            return main - SPERM, add - SPERM, rel, f"non-JURKAT per-line consensus ({len(sub)} annotations)"
        if g not in H.index: return set(), set(), "", "none"
        h = H.loc[g]; sp = lambda x: set(str(x).split(";")) - {"nan", ""} - SPERM
        return sp(h["Main location"]), sp(h["Additional location"]), h["Reliability"], "HPA consensus"

    NORM = {t.lower(): t for col in ("Main location", "Additional location") for x in H[col].dropna() for t in x.split(";")}
    HC = {}

    def hc(g):
        if g not in HC: HC[g] = hpa_cons(g)
        return HC[g]

    def feat(g):
        a = gacc.get(g)
        if not isinstance(a, str) or a not in U.index: return None
        r = U.loc[a]
        if isinstance(r, pd.DataFrame): r = r.iloc[0]
        ipr = set(str(r["InterPro"]).rstrip(";").split(";")) if isinstance(r["InterPro"], str) else set()
        return dict(ipr=ipr, kw=str(r["Keywords"]), mf=str(r["Gene Ontology (molecular function)"]),
                    cc=str(r["Gene Ontology (cellular component)"]), bp=str(r["Gene Ontology (biological process)"]),
                    up=str(r["Subcellular location [CC]"]))
    FE = {}

    def fe(g):
        if g not in FE: FE[g] = feat(g)
        return FE[g]

    UPF = pd.read_csv(ctx.frozen(p["uniprot_full"]), sep="\t", usecols=["Entry", "Transmembrane", "Intramembrane"]).set_index("Entry")

    def is_tm(g):
        a = gacc.get(g)
        if not isinstance(a, str) or a not in UPF.index: return False
        r = UPF.loc[a]; return isinstance(r["Transmembrane"], str) and "TRANSMEM" in r["Transmembrane"] or \
            isinstance(r["Intramembrane"], str) and "INTRAMEM" in r["Intramembrane"]

    def W(i): return json.loads(I.at[i, "comp_weights"]) if isinstance(I.at[i, "comp_weights"], str) else {}
    def genes(i): return str(I.at[i, "genes"]).split(";")
    def stoich(i): return [float(x) for x in str(I.at[i, "stoich"]).split(";")]
    def single(i): return len(set(genes(i))) == 1

    def frac_mass(i, S):
        m = [(g in S) * gmass.get(g, 0) * s for g, s in zip(genes(i), stoich(i))]; t = sum(gmass.get(g, 0) * s for g, s in zip(genes(i), stoich(i)))
        return sum(m) / t if t > 0 else 0.0

    def cyto_keys(i, w):
        sp = I.at[i, "signal_frac"] >= 0.5
        return [k for k in w if k in CYTO_KEYS or (k in SP_KEYS and not sp)]

    def wjson(w): return json.dumps({k: round(float(v), 4) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 1e-4})
    EV = {}

    def setw(i, nw, rule, gs, ev=None):
        ow = W(i); d = 0.5 * sum(abs(nw.get(k, 0) - ow.get(k, 0)) for k in set(nw) | set(ow))
        if d < 1e-4: return False
        r = moved.setdefault(rule, [0.0, set(), set()]); r[0] += d * I.at[i, "mass_pg"]; r[1] |= set(gs); r[2].add(I.at[i, "ing_id"])
        if ev is None: ev = " | ".join(f"{g}: {EV.get((rule.split('_')[0], g), '')}" for g in sorted(gs))
        log.append(dict(version=VER, gene=I.at[i, "genes"], ingredient=I.at[i, "ing_id"], rule_id=rule, field="comp_weights",
                        old_value=I.at[i, "comp_weights"], new_value=wjson(nw), evidence=ev, motivating_check=CHECK[rule.split("_")[0]],
                        mass_moved_pg=round(d * I.at[i, "mass_pg"], 6)))
        I.at[i, "comp_weights"] = wjson(nw); return True

    def targets(S, allow_complex=True, allow_tm=False):
        for i in range(len(I)):
            if bool(I.at[i, "transmembrane"]) != allow_tm and not (allow_tm is None): continue
            gs = set(genes(i)) & S
            if not gs: continue
            if single(i) or (allow_complex and frac_mass(i, S) >= 0.5): yield i, gs

    excl = lambda g: g in CUR or ACT_TUB.match(g) is not None
    rel_ok = lambda r: r in ("Enhanced", "Supported", "Approved")
    model_genes = [g for g in T.index if isinstance(g, str)]
    say(f"source table: {src} ({len(I0)} rows, {I0.mass_pg.sum():.3f} pg); model genes {len(model_genes)}; v4-curated {len(CUR)}; "
        f"JURKAT hold-out {len(HOLD)} / fit {len(FIT)}")

    # calibration (Tier 2, JURKAT fit half only)
    T2 = pd.read_csv(ctx.frozen(p["tier2_v5"]), sep="\t", index_col=0)
    man = pd.read_csv(ctx.frozen(p["hpa_jurkat_manifest"]), sep="\t")
    ab_genes = man.groupby("antibody").gene.agg(set)
    shared = {g for gs in ab_genes for g in gs if len(gs) > 1 and (gs & HOLD) and (gs & FIT)}
    T2f = T2[(T2.split == "fit") & ~T2.index.isin(shared)]
    jl = lambda g: set(str(T2.at[g, "hpa_jurkat"]).split(";")) if g in T2.index and isinstance(T2.at[g, "hpa_jurkat"], str) else set()
    regs = lambda g: dict((k, float(v)) for k, v in (x.split(":") for x in str(T2.at[g, "model_regions"]).split(";") if ":" in x))
    cyto_ref = [g for g in T2f.index if jl(g) and jl(g) <= {"cytosol"} and regs(g).get("cytosol", 0) >= 0.9]
    pm_ref = [g for g in T2f.index if "plasma membrane" in jl(g) and regs(g).get("PM", 0) >= 0.8]
    Rc = float(np.median(2 ** T2f.loc[cyto_ref, "rim_hpa"])) if cyto_ref else np.nan
    Rp = float(np.median(2 ** T2f.loc[pm_ref, "rim_hpa"])) if pm_ref else np.nan
    cal = dict(note="Tier 2 HPA rim (linear) on JURKAT fit-half genes; share s = (R_gene - R_cytosol) / (R_PM_TM - R_cytosol), "
                    "clipped to 0.3-0.7; antibodies shared with held-out genes excluded",
               shared_antibody_genes_excluded=sorted(shared), cyto_ref=cyto_ref, pm_ref=pm_ref, R_cyto=Rc, R_pm=Rp)
    say(f"calibration: cytosol reference {len(cyto_ref)} genes (R {Rc:.2f}), PM-TM reference {len(pm_ref)} genes (R {Rp:.2f})")

    def share_from_rim(gs):
        gs = [g for g in gs if g in T2f.index and T2f.at[g, "rim_hpa"] == T2f.at[g, "rim_hpa"]]
        if not gs or not (np.log2(Rp / Rc) >= MIN_CONTRAST): return S_DEFAULT, gs, None
        Rx = float(np.median(2 ** T2f.loc[gs, "rim_hpa"])); s = (Rx - Rc) / (Rp - Rc)
        return float(np.clip(s, CL0, CL1)), gs, s

    tr = json.load(open(ctx.frozen(p["tier2_transfer"])))

    def jcat(g):
        t = {x.lower() for x in jl(g)}; nt = {x.lower() for x in NUC_T}; n = bool(t & nt); c = bool(t - nt - {"nuclear membrane"})
        return "N" if n and not c else "NC" if n and c else "C" if c else None
    catc = {}
    for g in T2f.index:
        k = jcat(g)
        if k and T2f.at[g, "nc_hpa"] == T2f.at[g, "nc_hpa"]: catc.setdefault(k, []).append(float(T2f.at[g, "nc_hpa"]))
    cal["nc_by_jurkat_category"] = {k: dict(n=len(v), median_hpa_nc=round(float(np.median(v)), 3)) for k, v in catc.items()}
    c0, d0 = K["nc_transfer_model"]["c0"], K["nc_transfer_model"]["d0"]
    for k, v in cal["nc_by_jurkat_category"].items():
        v["implied_log2_r"] = round(((v["median_hpa_nc"] - tr["nc"]["a"]) / tr["nc"]["b"] - c0) / d0, 2)
    cal["r_by_category"] = R_CAT

    # R4: TM steady-state pools
    R4 = {}
    for g in model_genes:
        if excl(g): continue
        m_, a_, rel, srcs = hc(g)
        if not rel_ok(rel) or not ({"Endoplasmic reticulum", "Golgi apparatus"} & m_): continue
        R4[g] = (m_, a_, rel, srcs)
    n4 = 0
    for i, gs in targets(set(R4), allow_complex=True, allow_tm=True):
        w = W(i); pm = sum(v for k, v in w.items() if k == "plasma_membrane" or k.startswith("pm_"))
        tot = sum(w.values())
        if tot <= 0 or pm / tot < K["r4_min_pm"]: continue
        g0 = sorted(gs)[0]; m_ = R4[g0][0]
        keep = K["r4_keep"]["pm_main"] if (m_ & PM_T) else K["r4_keep"]["other"]
        tg = {k: 1.0 for k, t in (("er", "Endoplasmic reticulum"), ("golgi", "Golgi apparatus")) if t in m_}
        mv = pm - keep * tot
        if mv <= 0: continue
        nw = {k: v * (1 - mv / pm) if (k == "plasma_membrane" or k.startswith("pm_")) else v for k, v in w.items()}
        for k in tg: nw[k] = nw.get(k, 0) + mv / len(tg)
        for g in gs: EV[("R4", g)] = f"TM; HPA main {';'.join(sorted(R4[g][0]))} ({R4[g][2]}, {R4[g][3]}); model PM {pm / tot:.2f} -> {keep}"
        n4 += setw(i, nw, "R4_tm_pool", gs)
    say(f"R4: candidates {len(R4)}; ingredients changed {n4}")

    # R2: actin-binding -> cortex
    def actin_mech(g):
        f = fe(g)
        if not f: return None
        ev = []
        if ACT_GO.search(f["mf"]): ev.append("GO MF " + ";".join(sorted(set(m.group(0) for m in ACT_GO.finditer(f["mf"])))))
        if ACT_CC.search(f["cc"]): ev.append("GO CC " + ";".join(sorted(set(m.group(0) for m in ACT_CC.finditer(f["cc"])))))
        if "Actin-binding" in f["kw"]: ev.append("UniProt KW Actin-binding")
        d = sorted(ACT_IPR[x] + f" ({x})" for x in f["ipr"] & set(ACT_IPR))
        if d: ev.append("InterPro " + ", ".join(d))
        strong = bool(d) or "Actin-binding" in f["kw"] or (ACT_GO.search(f["mf"]) and ACT_CC.search(f["cc"]))
        return "; ".join(ev) if strong else None
    R2 = {}
    for g in model_genes:
        if excl(g) or g in V5CX or is_tm(g): continue
        m_, a_, rel, srcs = hc(g)
        hit = (m_ | a_) & R2_T
        if not rel_ok(rel) or not hit: continue
        mech = actin_mech(g)
        if not mech: continue
        R2[g] = f"HPA {';'.join(sorted(hit))} ({'main' if m_ & R2_T else 'add'}, {rel}, {srcs}); mechanism {mech}"
    s2, s2g, s2raw = share_from_rim(set(R2) & FIT)
    cal["R2"] = dict(share=s2, fit_genes=s2g, raw=s2raw)
    say(f"R2: actin-binding genes {len(R2)}; cortex share {s2:.2f}")
    R2done = set()
    for i, gs in targets(set(R2)):
        w = W(i)
        if set(w) & ANCHOR_KEYS: continue
        ck = cyto_keys(i, w); pool = sum(w[k] for k in ck)
        if pool <= 0: continue
        nw = {k: (v * (1 - s2) if k in ck else v) for k, v in w.items()}; nw["cortex"] = nw.get("cortex", 0) + s2 * pool
        for g in gs: EV[("R2", g)] = R2[g] + f"; cortex = {s2:.2f} x cytoplasmic pool"
        if setw(i, nw, "R2_cortex", gs): R2done |= gs

    # R3: peripheral membrane binders
    def memtgt(txt):
        t = txt.lower()
        if re.search(r"golgi|copi", t): return "golgi"
        if re.search(r"endosom|lysosom|clathrin-coated vesicle|multivesicular", t): return "endo"
        if re.search(r"cell membrane|plasma membrane|coated pit|cell junction|ruffle|lamellipodium|cell projection|focal adhesion", t): return "pm"
        return None

    def memb_mech(g, main_t=None):
        f = fe(g)
        if not f: return None, set()
        d = sorted(MEMB_IPR[x] + f" ({x})" for x in f["ipr"] & set(MEMB_IPR))
        if d: return "InterPro " + ", ".join(d), {"pm", "endo", "golgi"}
        ev, ok = [], set()
        for m in PERI_CC.finditer(f["cc"]):
            t = memtgt(m.group(0))
            if t and main_t and t in main_t: ok.add(t); ev.append(f"GO CC {m.group(0)}")
        if "Chaperone" in f["kw"]: return None, set()
        upx = re.sub(r"\{[^}]*\}", "", f["up"]).replace("SUBCELLULAR LOCATION:", " ")
        for st in re.split(r"\.\s+", upx):
            if "Peripheral membrane protein" in st and "Cytoplasmic side" in st:
                t = memtgt(st.split(";")[0])
                if t and main_t and t in main_t: ok.add(t); ev.append(f"UniProt '{st.strip()}'")
        return ("; ".join(ev) if ok else None), ok
    R3 = {}
    for g in model_genes:
        if excl(g) or g in V5CX or g in R2 or is_tm(g): continue
        m_, a_, rel, srcs = hc(g)
        tg = {}
        for t in m_ | a_:
            if t in MEMB_TGT: tg[MEMB_TGT[t]] = max(tg.get(MEMB_TGT[t], 0), 1.0 if t in m_ else 0.3)
        if not rel_ok(rel) or not tg: continue
        mech, ok = memb_mech(g, {MEMB_TGT[t] for t in m_ if t in MEMB_TGT})
        tg = {k: v for k, v in tg.items() if k in ok}
        if not mech or not tg: continue
        R3[g] = (tg, f"HPA {';'.join(sorted((m_ | a_) & set(MEMB_TGT)))} ({rel}, {srcs}); mechanism {mech}")
    s3, s3g, s3raw = share_from_rim({g for g in set(R3) & FIT if "pm" in R3[g][0]})
    cal["R3"] = dict(share=s3, fit_genes=s3g, raw=s3raw)
    say(f"R3: peripheral binders {len(R3)}; membrane share {s3:.2f}")
    R3done = set()
    for i, gs in targets(set(R3)):
        w = W(i)
        if set(w) & ANCHOR_KEYS: continue
        ck = cyto_keys(i, w); pool = sum(w[k] for k in ck)
        if pool <= 0: continue
        tg = {}
        for g in gs:
            for k, v in R3[g][0].items(): tg[k] = tg.get(k, 0) + v
        st = sum(tg.values())
        nw = {k: (v * (1 - s3) if k in ck else v) for k, v in w.items()}
        for k, v in tg.items(): nw[f"{k}_periph_cyto"] = nw.get(f"{k}_periph_cyto", 0) + s3 * pool * v / st
        for g in gs: EV[("R3", g)] = R3[g][1] + f"; membrane shells {', '.join(f'{k} {v:.1f}' for k, v in R3[g][0].items())} x {s3:.2f}"
        if setw(i, nw, "R3_peripheral", gs): R3done |= gs

    # R5: Golgi binders
    R5 = {}
    for g in model_genes:
        if excl(g) or g in V5CX or g in R2 or g in R3 or is_tm(g): continue
        f = fe(g)
        if not f: continue
        mech = [m.group(0) for m in GOLGI_CC.finditer(f["cc"])] + [m.group(0) for m in GOLGI_BP.finditer(f["bp"])]
        if not mech: continue
        m_, a_, rel, srcs = hc(g)
        hg = "Golgi apparatus" in (m_ | a_) and rel_ok(rel)
        ig = 0.0
        if g in T.index:
            ig = float(str(T.at[g, "inv_compartments"]).count("golgi:") and
                       sum(float(x.split(":")[1]) for x in str(T.at[g, "inv_compartments"]).split(";") if x.startswith("golgi:")))
        if not (hg or ig >= 0.2): continue
        R5[g] = f"GO {';'.join(sorted(set(mech)))}; " + (f"HPA Golgi ({rel}, {srcs})" if hg else f"inventory golgi share {ig:.2f}")
    s5 = K["r5_share"]
    for i, gs in targets(set(R5)):
        w = W(i)
        if set(w) & ANCHOR_KEYS: continue
        ck = cyto_keys(i, w); pool = sum(w[k] for k in ck)
        if pool <= 0: continue
        nw = {k: (v * s5 if k in ck else v) for k, v in w.items()}; nw["golgi_periph_cyto"] = nw.get("golgi_periph_cyto", 0) + s5 * pool
        for g in gs: EV[("R5", g)] = R5[g] + "; golgi_periph_cyto = 0.5 x cytoplasmic pool"
        setw(i, nw, "R5_golgi", gs)

    # R6a: NPC-proximal shell
    R6a = {}
    sh6 = K["r6a_shares"]
    for g in model_genes:
        if excl(g) or is_tm(g): continue
        f = fe(g)
        if not f: continue
        cc = [m.group(0) for m in NPC_CC.finditer(f["cc"])]; mf = [m.group(0) for m in NPC_MF.finditer(f["mf"])]
        if cc: R6a[g] = (sh6["nuclear_pore_cc"], f"GO CC {';'.join(sorted(set(cc)))}")
        elif mf: R6a[g] = (sh6["transport_receptor_mf"], f"GO MF {';'.join(sorted(set(mf)))}")
    R6adone = set()
    for i, gs in targets(set(R6a)):
        w = W(i)
        ck = cyto_keys(i, w) + [k for k in ("nucleus", "nuclear_envelope") if k in w]; pool = sum(w[k] for k in ck)
        if pool <= 0: continue
        f = float(np.average([R6a[g][0] for g in gs]))
        nw = {k: (v * (1 - f) if k in ck else v) for k, v in w.items()}; nw["npc"] = nw.get("npc", 0) + f * pool
        for g in gs: EV[("R6a", g)] = R6a[g][1] + f"; npc = {f:.2f} x (nucleus + cytoplasmic) pool"
        if setw(i, nw, "R6a_npc", gs): R6adone |= gs

    # R1: nucleus / cytosol from imaging
    def r_itz():
        Z = pd.read_csv(ctx.frozen(p["itzhak"]), sep="\t")
        Z = Z[Z.primary_for_gene & (Z.n_good_mnc_profiles >= IZ["min_profiles"])].drop_duplicates("gene").set_index("gene")
        lo, hi = IZ["clip"]; base = IZ["nuclear_base"]
        out_ = {}
        for g, r in Z.iterrows():
            fn = max(r.frac_nuclear / max(r.frac_nuclear + r.frac_cytosol, 1e-9) - base, 0.0) / (1 - base)
            out_[g] = float(np.clip((fn / VN) / max((1 - fn) / VC, 1e-9), lo, hi)) if 0 < fn < 1 else (hi if fn >= 1 else lo)
        return out_
    ITZ = r_itz()
    # revision 2026-11-18 (preset common): the R1 target goes through the same N/C target function as the v0 N/C model
    # (combiner.nc_target: Itzhak prior, bound), with the per-gene prior / exclusivity inputs of v0 (localisation/
    # nc_model.tsv); chromatin-bound genes are skipped (their pool stays in the nucleus)
    NCO, NCT = None, None
    if R1_RULE:
        from wholecell.inventory.localisation import combiner as LCB
        NCO = LOC._resolve_data(LOC.options(ctx.params("localisation"))["combiner"], ctx.params("localisation")).get("nc_model") or {}
        if NCO.get("prior") or NCO.get("bound") or NCO.get("chromatin"):
            NCT = pd.read_csv(ctx.need("localisation/nc_model.tsv"), sep="\t").drop_duplicates("gene").set_index("gene")
        else:
            NCO = None

    def r1_target(g, cat, r):
        """(r, note) for preset common with the revision options; (r, '') otherwise; (None, note) for chromatin genes"""
        if NCO is None:
            return r, ""
        row = NCT.loc[g] if g in NCT.index else None
        if row is not None and bool(row.chromatin):
            return None, "chromatin-bound (pool stays nuclear)"
        fnf = float(row.fn_prior) if row is not None and pd.notna(row.fn_prior) else None
        ex = {"n": (int(row.ex_n_other) if row is not None and pd.notna(row.ex_n_other) else 0) + (cat == "N"),
              "c": (int(row.ex_c_other) if row is not None and pd.notna(row.ex_c_other) else 0) + (cat == "C")}
        rr, fn, ti = LCB.nc_target(r, fnf, ex, NCO, VN, VC, tie_break=(cat == "NC" and g in ITZ))
        note = (f"; Itzhak prior r {ti['r_prior']:.3g}" if "r_prior" in ti else "") + (f"; bounded {ti['bounded'][0]} -> {ti['bounded'][1]}" if ti.get("bounded") else "")
        return rr, note
    R1 = {}
    # universe: legacy = genes of the v5 pack (gate: placed nucleus + cytosol >= 0.8); common = genes of the current
    # table (gate on the species' intended shares below)
    r1_genes = sorted({g for x in I.genes.astype(str) for g in x.split(";")}) if R1_RULE else model_genes
    for g in r1_genes:
        if (excl(g) and g not in CUR_MITO) or ACT_TUB.match(g) or g in R2done | R3done | R6adone or is_tm(g): continue
        if not R1_RULE and (g not in T.index or T.at[g, "m_nucleus"] + T.at[g, "m_cytosol"] < 0.8): continue
        m_, a_, rel, srcs = hc(g)
        if not rel_ok(rel) or (m_ & ORG_T): continue
        nm, cm = bool(m_ & NUC_T), bool(m_ - NUC_T - ORG_T - PM_T)
        na, ca = bool(a_ & NUC_T), bool(a_ - NUC_T - ORG_T - PM_T)
        cat = "N" if nm and not cm and not ca else "Nc" if nm and not cm else "NC" if nm and cm else "Cn" if cm and na else "C" if cm else None
        if not cat: continue
        r = R_CAT[cat]; tb = ""
        if cat == "NC" and g in ITZ: r = ITZ[g]; tb = f"; Itzhak 2016 tie-breaker r {r:.2f}"
        r, note = r1_target(g, cat, r)
        if r is None: continue
        tb += note
        R1[g] = (cat, r, f"HPA main {';'.join(sorted(m_))} add {';'.join(sorted(a_)) or '-'} ({rel}, {srcs}); category {cat} -> target N/C {r:.3g}{tb}")
    nA = nB = 0
    rowsR1 = []
    for i, gs in targets(set(R1)):
        w = W(i)
        if set(w) & {"cortex", "lamina", "npc"}: continue
        ck = cyto_keys(i, w); wc = sum(w[k] for k in ck); wn = w.get("nucleus", 0.0); wl = w.get("nucleolus", 0.0) + w.get("speckle", 0.0)
        pool = wc + wn + wl
        if pool <= 0: continue
        if R1_RULE and pool < 0.8 * sum(w.values()): continue        # intended nucleus + cytoplasm share >= 0.8
        cats = {R1[g][0] for g in gs}
        fam = {"N": "n", "Nc": "n", "C": "c", "Cn": "c", "NC": "nc"}
        if len({fam[c] for c in cats}) > 1: continue
        lr = [np.log2(R1[g][1]) for g in gs]; r = float(2 ** np.average(lr, weights=[gmass.get(g, 1) for g in gs]))
        rm = ((wn + wl) / VN) / max(wc / VC, 1e-12)
        if abs(np.log2(max(rm, 1e-6) / r)) <= 1.0: continue
        if cats <= {"N", "Nc"} and rm > r: continue
        if cats <= {"C", "Cn"} and rm < r: continue
        fn = r * VN / (r * VN + VC); new_n = max(fn * pool - wl, 0.0); new_c = pool - new_n - wl
        nw = {k: v for k, v in w.items() if k not in ck and k != "nucleus"}
        if new_n > 0: nw["nucleus"] = new_n
        if wc > 0:
            for k in ck: nw[k] = w[k] * new_c / wc
        elif new_c > 0: nw["cytosol"] = new_c
        for g in gs: EV[("R1", g)] = R1[g][2] + f"; model N/C {rm:.3g} -> {r:.3g} (nuclear share of the soluble pool {(wn + wl) / pool:.2f} -> {fn:.2f})"
        rowsR1.append(dict(rule_set="jurkat.v6", rule_id="v6.R1_nc", key_type="species", key=I.at[i, "ing_id"],
                           field="comp_weights", op="set", value=json.dumps(nw),
                           reason=" | ".join(f"{g}: {EV.get(('R1', g), '')}" for g in sorted(gs)),
                           citation=("HPA imaging category -> target N/C (rules_v6.yaml r_by_category); gate on the INTENDED "
                                     "shares of the current table, this build's volumes (B5 criteria §2)") if R1_RULE else
                                    ("HPA imaging category -> target N/C (rules_v6.yaml r_by_category); gate on the frozen "
                                     "v5-pack placed shares (params jurkat.v6.v5pack_loc_compare)"),
                           genes=";".join(sorted(gs)), direction="raise" if r > rm else "lower"))
    if R1_RULE and rowsR1:                                          # capacity check (B5 criteria §2)
        pos = {I.at[i, "ing_id"]: i for i in range(len(I))}
        dl = []
        for rw in rowsR1:
            i = pos[rw["key"]]; w0 = W(i); nw = json.loads(rw["value"])
            dl.append((i, sum(nw.get(k, 0) for k in ("nucleus", "nucleolus", "speckle")) - sum(w0.get(k, 0) for k in ("nucleus", "nucleolus", "speckle"))))
        sc, info = capacity_scale(I, dl, VN, VC)
        say(f"R1 (common, intended shares): capacity {json.dumps({k: round(v, 3) for k, v in info.items()})}")
        if sc < 1.0:
            for rw, (i, d) in zip(rowsR1, dl):
                if d > 0:
                    w0 = W(i); nw = json.loads(rw["value"])
                    rw["value"] = json.dumps({k: w0.get(k, 0.0) + sc * (nw.get(k, 0.0) - w0.get(k, 0.0)) for k in set(w0) | set(nw)
                                              if w0.get(k, 0.0) + sc * (nw.get(k, 0.0) - w0.get(k, 0.0)) > 0})
    # frozen-input override table (preset legacy_v7 only): rows above, applied by the common override function
    OVR1 = OV.frame(rowsR1)
    if I.ing_id.duplicated().any():
        raise ValueError("v6 R1: species ids are not unique; a species-keyed override table is ambiguous")

    def put_r1(i, row, nw):
        nonlocal nA, nB
        if setw(i, nw, "R1_nc", set(row.genes.split(";")), ev=row.reason):
            if row.direction == "raise": nA += 1
            else: nB += 1
    OV.apply(OVR1, lambda kt, key: [int(np.nonzero((I.ing_id == key).values)[0][0])], lambda i, f: W(i), put_r1)
    say(f"R1: genes {len(R1)}; ingredients raised {nA}, lowered {nB}")

    # R6b: nuclear speckles
    R6b = {}
    s6b = K["r6b_shares"]
    for g in model_genes:
        if excl(g) or is_tm(g): continue
        m_, a_, rel, srcs = hc(g)
        if rel_ok(rel) and "Nuclear speckles" in m_:
            R6b[g] = (s6b["sr_type"] if SR_RE.match(g) else s6b["other"], f"HPA main Nuclear speckles ({rel}, {srcs})" + ("; SR-type" if SR_RE.match(g) else ""))
    vspk = float(p["speckle_volume_um3"])
    plan = []
    for i, gs in targets(set(R6b)):
        w = W(i)
        if w.get("nucleus", 0) <= 0: continue
        plan.append((i, gs, float(np.average([R6b[g][0] for g in gs])), w))
    m_spk = sum(f * w["nucleus"] * I.at[i, "mass_pg"] for i, gs, f, w in plan)
    dens = m_spk / vspk * 1e3; sc = min(1.0, SPK_MAX / dens) if dens > 0 else 1.0
    say(f"R6b: speckle genes {len(R6b)}; speckle volume {vspk:.2f} um3; requested {m_spk:.3f} pg = {dens:.0f} mg/ml -> scale {sc:.3f}")
    for i, gs, f, w in plan:
        x = sc * f * w["nucleus"]; nw = dict(w); nw["nucleus"] = w["nucleus"] - x; nw["speckle"] = nw.get("speckle", 0) + x
        for g in gs: EV[("R6b", g)] = R6b[g][1] + f"; speckle = {R6b[g][0]} x nuclear pool x {sc:.3f} (capacity)"
        setw(i, nw, "R6b_speckle", gs)
    cal["R6b"] = dict(volume_um3=vspk, requested_pg=m_spk, scale=sc, cap_mg_ml=SPK_MAX)

    for k, (m, g, ii) in moved.items(): say(f"  {k} | {len(ii)} | {len(g)} | {m:.3f}")
    L = pd.DataFrame(log, columns=["version", "gene", "ingredient", "rule_id", "field", "old_value", "new_value", "evidence",
                                   "motivating_check", "mass_moved_pg"])
    G = []
    for rule, dct in (("R1", {g: v[2] for g, v in R1.items()}), ("R2", R2), ("R3", {g: v[1] for g, v in R3.items()}),
                      ("R4", {g: f"HPA main {';'.join(sorted(v[0]))} ({v[2]}, {v[3]})" for g, v in R4.items()}), ("R5", R5),
                      ("R6a", {g: v[1] for g, v in R6a.items()}), ("R6b", {g: v[1] for g, v in R6b.items()})):
        changed = moved.get({"R1": "R1_nc", "R2": "R2_cortex", "R3": "R3_peripheral", "R4": "R4_tm_pool", "R5": "R5_golgi",
                             "R6a": "R6a_npc", "R6b": "R6b_speckle"}[rule], [0, set()])[1]
        for g, e in dct.items():
            G.append(dict(version=VER, gene=g, rule_id=rule, selected_evidence=e, changed=g in changed,
                          jurkat_half="holdout" if g in HOLD else "fit" if g in FIT else "", motivating_check=CHECK[rule]))
    say(f"total mass {I0.mass_pg.sum():.4f} -> {I.mass_pg.sum():.4f} pg; {len(L)} field changes in {L.ingredient.nunique()} ingredients")
    for c in I0.columns:
        if c != "transmembrane": I[c] = I[c].astype(I0[c].dtype)
    I["transmembrane"] = I.transmembrane.astype(bool)
    return {"ingredients": I, "changes": normalise_changes(L, VER),
            "extra": {"changes_legacy.tsv": L, "rule_genes.tsv": pd.DataFrame(G),
                      "calibration.json": json.dumps(cal, indent=1, default=float), "overrides_R1_nc.tsv": OVR1}}


apply_v6.OUTPUTS = ["changes_legacy.tsv", "rule_genes.tsv", "calibration.json", "overrides_R1_nc.tsv"]
apply_v6.PREV = ["rules/jurkat.v4/changes_legacy.tsv", "rules/jurkat.v5/changes_legacy.tsv", "v0/inventory_v0.tsv",
                 "localisation/nc_model.tsv"]
register("inventory_rule", "jurkat.v6")(apply_v6)
