# origin: jurkat/yc/pipeline/jinventory_v7.py (legacy workspace, 2026-11-15) + jvirtual.load_context/_segs (MT band lengths)
# changes: inventory_rule "jurkat.v7": reads the previous (v6) step's table and writes a new one (no snapshot copy);
# the lamina pre-pass consumption is an explicit input (params lamina_consumption: a frozen input name, "build:<path
# relative to the build dir>" for the rebuilt lamina of the preplaced stage, or a path; without it, rule L3 runs as the
# separate rule set jurkat.v7_lamina, v9); the PCM sphere volume (legacy:
# Monte Carlo on jgeom_cell1_v7.npz with cell3d, rng 3) is a config value (params pcm_volume_um3; recompute with the
# geometry engine for a new geometry); the MT band volumes are recomputed from the frozen centerlines (as legacy);
# V7_CHROMATIN_JSON env -> params chromatin_json. Constants in plugins/jurkat/data/rules_v7.yaml. Numerics identical.
"""Inventory revision v7: R4b (TCR/CD3 steady-state pools), R6c (PCM, capped), R7 (MAP +TIPs / lattice MAPs, capped),
R6d (vimentin cage), R6e (NPC channel receptors), L3 (lamin pool = lattice consumption), F1 (CPX-2511 nucleolar),
N8 (nucleosome copies = placed v8 chromatin) + non-inventory history rows (prepass / geometry / packer)."""
import json
import os
import re

import numpy as np
import pandas as pd

from wholecell.inventory import localisation as LOC
from wholecell.inventory.localisation import overrides as OV
from wholecell.registry import register
from .common import data, normalise_changes, read_list

VER = "v7"
# the step order of this rule set = the localisation preset's rules."jurkat.v7" list (checked at run time)
STEPS = [("v7.R4b", "overrides"), ("v7.R6c", "legacy"), ("v7.R7", "legacy"), ("v7.R6d", "legacy"), ("v7.R6e", "legacy"),
         ("v7.L3", "legacy"), ("v7.F1", "overrides"), ("v7.N8", "legacy"), ("v7.history", "legacy")]
AV = 6.022e23


def resolve_input(ctx, v):
    """frozen input name | 'build:<rel>' | path"""
    return ctx.resolve(v)


def mt_volumes(centerlines, band, plus_l, plus_min_l):
    """band volumes (um3) of mt_lattice / mt_plus from the centerline lengths (annulus band[0]-band[1] nm)"""
    z = np.load(centerlines); P = z["points_model_nm"].astype(np.float64); off = z["offsets"]; fl = z["flipped_for_lattice"]
    polys = [P[off[i]:off[i + 1]] for i in range(len(off) - 1)]
    plus = []
    for pp, f in zip(polys, fl):
        q = pp[::-1] if f else pp
        L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]
        if L[-1] < plus_min_l: continue
        plus.append(q[L >= L[-1] - plus_l])
    A = [pp[:-1] for pp in polys if len(pp) > 1]; B = [pp[1:] for pp in polys if len(pp) > 1]
    A, B = (np.concatenate(A), np.concatenate(B)) if A else (np.zeros((0, 3)), np.zeros((0, 3)))
    plus_um = float(sum(np.linalg.norm(np.diff(q, axis=0), axis=1).sum() for q in plus) / 1e3)
    ann = np.pi * (band[1] ** 2 - band[0] ** 2) / 1e6
    L_all = float(np.linalg.norm(B - A, axis=1).sum()) / 1e3
    return dict(lat_um3=ann * L_all, plus_um3=ann * plus_um, L_um=L_all, plus_um=plus_um, n_plus=len(plus))


def apply_v7(ctx, src):
    p = ctx.params("jurkat.v7")
    K = data("rules_v7.yaml")
    CYTO_KEYS, SP_KEYS, SPERM = set(K["cyto_keys"]), set(K["sp_keys"]), set(K["sperm_terms"])
    LOC.check_rule_steps(ctx, "jurkat.v7", STEPS)
    OV_R4B, OV_F1 = LOC.override_rows(ctx, "jurkat.v7", "v7.R4b"), LOC.override_rows(ctx, "jurkat.v7", "v7.F1")
    S_PCM, PCM_MAX = K["pcm"]["share"], K["pcm"]["cap_mg_ml"]
    MK = K["mt"]; S_PLUS, S_LAT, MT_MAX = MK["plus_share"], MK["lattice_share"], MK["cap_mg_ml"]
    S_VIM = K["vimentin_share"]; N_FG, N_PORES = K["npc_fg"]["receptors_per_pore"], K["npc_fg"]["n_pores"]
    PCM_T = set(K["pcm_terms"]); PCM_CC = re.compile(K["pcm_cc_regex"])
    PLUS_IPR, LAT_IPR = dict(K["plus_ipr"]), dict(K["lattice_ipr"]); MT_T = set(K["mt_terms"]); IF_IPR = K["if_ipr"]
    NPC_MF = re.compile(K["npc_mf_regex"]); CHECK = dict(K["check"])
    I = pd.read_csv(src, sep="\t"); I0 = I.copy()
    H = pd.read_csv(ctx.frozen(p["hpa_subcellular"]), sep="\t").drop_duplicates("Gene name").set_index("Gene name")
    C = pd.read_csv(ctx.frozen(p["hpa_celllines"]), sep="\t", keep_default_na=False)
    U = pd.read_csv(ctx.frozen(p["uniprot_domains"]), sep="\t").set_index("Entry")
    inv = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t").drop_duplicates("gene").set_index("gene")
    gmass = inv.mass_da.to_dict(); gacc = inv.uniprot.to_dict()
    HOLD = read_list(ctx.frozen(p["holdout_genes"]))
    ch = pd.read_csv(ctx.need("rules/jurkat.v4/changes_legacy.tsv"), sep="\t")
    cur = ch[ch.field.isin(["comp_weights", "transmembrane"]) & ch.reason.astype(str).str.match(r'^"?[bcf]')]
    CUR = {g for x in cur.gene.astype(str) for g in x.split(";")}
    NORM = {t.lower(): t for col in ("Main location", "Additional location") for x in H[col].dropna() for t in x.split(";")}
    log, moved, grow = [], {}, []
    say = ctx.say
    HC = {}

    def hc(g):
        if g in HC: return HC[g]
        if g in HOLD:
            sub = C[(C.gene == g) & (C.cell_line.str.upper() != "JURKAT") & (C.ab_reliability != "uncertain") & (C.locations != "")]
            if not len(sub): r = (set(), set(), "non-JURKAT per-line consensus (none)")
            else:
                locs = [set(x.lower().split(";")) for x in sub.locations]
                fr = {t: sum(t in l for l in locs) / len(locs) for t in set().union(*locs)}
                fx = lambda t: NORM.get(t.lower(), t)
                r = ({fx(t) for t, f in fr.items() if f >= 0.5} - SPERM, {fx(t) for t, f in fr.items() if 0.25 <= f < 0.5} - SPERM,
                     f"non-JURKAT per-line consensus ({len(sub)})")
        elif g in H.index:
            h = H.loc[g]; sp = lambda x: set(str(x).split(";")) - {"nan", ""} - SPERM
            r = (sp(h["Main location"]), sp(h["Additional location"]), "HPA consensus")
        else: r = (set(), set(), "none")
        HC[g] = r; return r

    FE = {}

    def fe(g):
        if g in FE: return FE[g]
        a = gacc.get(g); r = None
        if isinstance(a, str) and a in U.index:
            x = U.loc[a]; x = x.iloc[0] if isinstance(x, pd.DataFrame) else x
            r = dict(ipr=set(str(x["InterPro"]).rstrip(";").split(";")) if isinstance(x["InterPro"], str) else set(),
                     mf=str(x["Gene Ontology (molecular function)"]), cc=str(x["Gene Ontology (cellular component)"]), kw=str(x["Keywords"]))
        FE[g] = r; return r

    def Wt(i): return json.loads(I.at[i, "comp_weights"]) if isinstance(I.at[i, "comp_weights"], str) else {}
    def genes(i): return str(I.at[i, "genes"]).split(";")
    def stoich(i): return [float(x) for x in str(I.at[i, "stoich"]).split(";")]

    def frac_mass(i, S):
        gs, st = genes(i), stoich(i); t = sum(gmass.get(g, 0) * s for g, s in zip(gs, st))
        return sum((g in S) * gmass.get(g, 0) * s for g, s in zip(gs, st)) / t if t > 0 else 0.0

    def cyto_keys(i, w):
        sp = I.at[i, "signal_frac"] >= 0.5
        return [k for k in w if k in CYTO_KEYS or (k in SP_KEYS and not sp)]

    def wjson(w): return json.dumps({k: round(float(v), 4) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 1e-4})
    EV = {}

    def setw(i, nw, rule, gs):
        ow = Wt(i); d = 0.5 * sum(abs(nw.get(k, 0) - ow.get(k, 0)) for k in set(nw) | set(ow))
        if d < 1e-4: return False
        r = moved.setdefault(rule, [0.0, set(), set()]); r[0] += d * I.at[i, "mass_pg"]; r[1] |= set(gs); r[2].add(I.at[i, "ing_id"])
        log.append(dict(version=VER, gene=I.at[i, "genes"], ingredient=I.at[i, "ing_id"], rule_id=rule, field="comp_weights",
                        old_value=I.at[i, "comp_weights"], new_value=wjson(nw),
                        evidence=" | ".join(f"{g}: {EV.get((rule, g), '')}" for g in sorted(gs)),
                        motivating_check=CHECK[rule.split("_")[0]], mass_moved_pg=round(d * I.at[i, "mass_pg"], 6)))
        I.at[i, "comp_weights"] = wjson(nw); return True

    def targets(S):
        for i in range(len(I)):
            if bool(I.at[i, "transmembrane"]): continue
            gs = set(genes(i)) & S
            if not gs: continue
            if len(set(genes(i))) == 1 or frac_mass(i, S) >= 0.5: yield i, gs

    def move(i, key, share, rule, gs):
        w = Wt(i); ck = cyto_keys(i, w); pool = sum(w[k] for k in ck)
        if pool <= 0 or share <= 0: return 0.0
        x = share * pool; nw = dict(w)
        for k in ck: nw[k] = w[k] * (1 - share)
        nw[key] = nw.get(key, 0) + x
        return x if setw(i, nw, rule, gs) else 0.0

    model_genes = sorted({g for i in range(len(I)) for g in genes(i)})
    say(f"source {src} ({len(I)} rows, {I.mass_pg.sum():.3f} pg); JURKAT hold-out {len(HOLD)} genes (non-JURKAT consensus)")

    def find(kt, key):                      # species keys only (first species with that ing_id)
        k = I.index[I.ing_id == key]
        return [k[0]] if len(k) else []

    def get(i, field):
        return Wt(i) if field == "comp_weights" else I.at[i, field]

    # R4b TCR/CD3 (override rows v7.R4b)
    def put_r4b(i, r, nw):
        gs = set(genes(i))
        for g in gs: EV[("R4b", g)] = r.reason
        setw(i, nw, "R4b", gs)
    OV.apply(OV_R4B, find, get, put_r4b)

    # R6c PCM
    pcm = []
    for g in model_genes:
        if g in CUR: continue
        m, a, srcs = hc(g); f = fe(g)
        if not (m & PCM_T) or f is None or not PCM_CC.search(f["cc"]): continue
        pcm.append(g); EV[("R6c", g)] = f"HPA main {';'.join(sorted(m & PCM_T))} ({srcs}); GO CC {PCM_CC.search(f['cc']).group(0)}"
    req = {}
    for i, gs in targets(set(pcm)):
        w = Wt(i); pool = sum(w[k] for k in cyto_keys(i, w)); req[i] = (gs, S_PCM * pool * I.at[i, "copies"] * I.at[i, "mass_da"])
    vol_pcm = float(p["pcm_volume_um3"])
    mass_pcm = sum(v[1] for v in req.values()) / AV * 1e12
    cap = PCM_MAX * vol_pcm * 1e-3
    sc = min(1.0, cap / mass_pcm) if mass_pcm > 0 else 1.0
    for i, (gs, _) in req.items(): move(i, "pcm", S_PCM * sc, "R6c", gs)
    say(f"R6c PCM: {len(pcm)} genes, {len(req)} ingredients; requested {mass_pcm:.4f} pg in {vol_pcm:.3f} um3 -> scale {sc:.3f}")
    grow += [dict(gene=g, rule="R6c", evidence=EV[("R6c", g)]) for g in pcm]

    # R7 MAPs
    plus, latt = [], []
    MTX = re.compile(K["mt_exclude_regex"])
    for g in model_genes:
        if g in CUR or MTX.match(g): continue
        f = fe(g)
        if f is None: continue
        m, a, srcs = hc(g)
        loc_ok = bool((m | a) & MT_T)
        ii = [i for i in I.index[I.genes == g]]
        cs = max([Wt(i).get("cytoskeleton", 0) for i in ii] or [0])
        if not (loc_ok or cs >= 0.2): continue
        why = f"HPA {';'.join(sorted((m | a) & MT_T))} ({srcs})" if loc_ok else f"inventory cytoskeleton share {cs:.2f}"
        dp = [PLUS_IPR[x] for x in f["ipr"] if x in PLUS_IPR]; dl = [LAT_IPR[x] for x in f["ipr"] if x in LAT_IPR]
        if dp: plus.append(g); EV[("R7_plus", g)] = f"InterPro {', '.join(dp)}; {why}"
        elif dl: latt.append(g); EV[("R7_lattice", g)] = f"InterPro {', '.join(dl)}; {why}"
    ctx_mt = mt_volumes(ctx.frozen(p["mt_centerlines"]), tuple(MK["band_nm"]), MK["plus_length_nm"], MK["plus_min_mt_length_nm"])
    for key, gl, s0, rule, vol in (("mt_plus", plus, S_PLUS, "R7_plus", ctx_mt["plus_um3"]), ("mt_lattice", latt, S_LAT, "R7_lattice", ctx_mt["lat_um3"])):
        req = {}
        for i, gs in targets(set(gl)):
            w = Wt(i); pool = sum(w[k] for k in cyto_keys(i, w)); req[i] = (gs, s0 * pool * I.at[i, "copies"] * I.at[i, "mass_da"])
        mass = sum(v[1] for v in req.values()) / AV * 1e12; cap = MT_MAX * vol * 1e-3
        sc = min(1.0, cap / mass) if mass > 0 else 1.0
        for i, (gs, _) in req.items(): move(i, key, s0 * sc, rule, gs)
        say(f"{rule}: {len(gl)} genes, {len(req)} ingredients; requested {mass:.4f} pg in {vol:.3f} um3 -> scale {sc:.3f}")
        grow += [dict(gene=g, rule=rule, evidence=EV[(rule, g)]) for g in gl]

    # R6d vimentin
    vim = []
    for g in model_genes:
        f = fe(g)
        if f is None or IF_IPR not in f["ipr"] or g.startswith("LMN") or g in CUR: continue
        m, a, srcs = hc(g)
        if "Intermediate filaments" not in m: continue
        vim.append(g); EV[("R6d", g)] = f"InterPro IF rod; HPA main Intermediate filaments ({srcs})"
    for i, gs in targets(set(vim)): move(i, "vim_cage", S_VIM, "R6d", gs)
    say(f"R6d vimentin cage: {vim}")
    grow += [dict(gene=g, rule="R6d", evidence=EV[("R6d", g)]) for g in vim]

    # R6e NPC FG docking
    rec = [g for g in model_genes if fe(g) is not None and NPC_MF.search(fe(g)["mf"])]
    tl = [(i, gs) for i, gs in targets(set(rec)) if Wt(i).get("npc", 0) > 0]
    tot = sum(I.at[i, "copies"] for i, _ in tl); want = N_FG * N_PORES
    for i, gs in tl:
        w = Wt(i); x = min(w["npc"], want / tot)
        nw = dict(w); nw["npc"] = w["npc"] - x; nw["npc_fg"] = nw.get("npc_fg", 0) + x
        for g in gs: EV[("R6e", g)] = f"GO MF nuclear import/export receptor; {N_FG} receptors per pore in total, split by copies"
        setw(i, nw, "R6e", gs)
    grow += [dict(gene=g, rule="R6e", evidence=EV.get(("R6e", g), "")) for g in rec]

    # L3 nucleoplasmic lamin pool: lamina share = lattice copies / copies (explicit lamina pre-pass input). Without
    # params lamina_consumption the rule is left to the rule set jurkat.v7_lamina (apply_v7_lamina, run after this one).
    lam_path = resolve_input(ctx, p["lamina_consumption"]) if p.get("lamina_consumption") else None
    lj0 = json.load(open(lam_path)) if lam_path else {}
    if lam_path:
        l3_apply(I, lj0, Wt, setw, EV)
        say(f"L3 lamin pool ({lam_path}): {lj0['consumed']} in the lattice; rest -> nucleus")
    else:
        say("L3 lamin pool: deferred to rule set jurkat.v7_lamina (no params lamina_consumption)")

    # F1 CPX-2511 (SSU processome) compartment (override rows v7.F1: comp_weights, transmembrane, n_tm; v7 logged the
    # transmembrane change unconditionally and the n_tm change not at all - kept)
    def put_f1(i, r, new):
        if r.field == "comp_weights":
            EV[("F1", r.key)] = r.reason
            setw(i, new, "F1", {r.key})
        elif r.field == "transmembrane":
            oldtm = I.at[i, "transmembrane"]; I.at[i, "transmembrane"] = new
            log.append(dict(version=VER, gene=I.at[i, "genes"], ingredient=r.key, rule_id="F1", field="transmembrane", old_value=oldtm,
                            new_value=new, evidence=EV[("F1", r.key)], motivating_check=CHECK["F1"], mass_moved_pg=0.0))
        else:
            I.at[i, r.field] = new
    OV.apply(OV_F1, find, get, put_f1)

    # N8 nucleosome copies = placed genome-mapped chromatin
    cj8 = resolve_input(ctx, p["chromatin_json"])
    if os.path.exists(cj8):
        n8 = int(json.load(open(cj8))["n_nucleosomes"]); k = I.index[I.ing_id == "NUCLEOSOME"]
        if len(k) and n8 > 0:
            i = k[0]; old = int(I.at[i, "copies"])
            I.at[i, "copies"] = n8; I.at[i, "mass_pg"] = n8 * I.at[i, "mass_da"] / AV * 1e12
            log.append(dict(version=VER, gene=I.at[i, "genes"], ingredient="NUCLEOSOME", rule_id="N8", field="copies", old_value=old, new_value=n8,
                            evidence=f"genome-mapped chromatin v8 places {n8:,} nucleosomes (5.935 Gbp 2n with the Jurkat deletions, "
                                     f"210.5 bp/nucleosome); gap {old - n8:,} = 6.2 Gbp/200 bp estimate vs GRCh38 + NRL 203 + NFRs",
                            motivating_check="chromatin v8 report (inventory 31 M vs 28.19 M placed); coordinator DECISION",
                            mass_moved_pg=round((old - n8) * I.at[i, "mass_da"] / AV * 1e12, 4)))
            say(f"N8 NUCLEOSOME copies {old:,} -> {n8:,}")

    # non-inventory rows (per-gene history)
    nonint = []

    def nrow(gene, ing, rule, field, old, new, ev, chk, m=0.0):
        nonint.append(dict(version=VER, gene=gene, ingredient=ing, rule_id=rule, field=field, old_value=old, new_value=new,
                           evidence=ev, motivating_check=chk, mass_moved_pg=m))
    cjp = resolve_input(ctx, p["cytoskeleton_json"])
    cj = json.load(open(cjp)) if os.path.exists(cjp) else {}
    c6 = json.load(open(resolve_input(ctx, p["cytoskeleton_prev_json"])))
    if cj:
        nrow("TUBA1B;TUBB4B", "TUBA1B-TUBB4B", "MT3b", "prepass", f"consumed {c6['consumed']['TUBA1B-TUBB4B']:,} (procedural MTs + centrioles)",
             f"consumed {cj['consumed']['TUBA1B-TUBB4B']:,} (stage-3 b MTs + centrioles at the MTOC candidate)",
             "cytoskel/stage3 variant b (EM-scaffolded, 1.456 mm); DECISION stage 3 -> v7", "stage 3 report: old pre-pass radial, 2.3 mm; HPA aster ratio 3.14 vs 2.20")
        nrow("ACTB", "ACTB", "FA_trim", "prepass", f"consumed {c6['consumed']['ACTB']:,}", f"consumed {cj['consumed']['ACTB']:,}",
             "F-actin subunits clashing with the new MTs / centrioles or the v7 geometry removed", "stage 3 clash counts (16k subunits within 16 nm of the new axes)")
    if lj0:
        nonint += l2_rows(lj0)
    nrow("(Golgi TM + Golgi shells)", "*", "R5b_cisternae", "geometry", "EM blob, area factor 4.5",
         "juxta-centrosomal procedural Golgi (jgolgi_ribbon.py) with cisternal layers (jgolgi.py); factor = 4.5 x EM blob area / measured",
         "Ladinsky et al. 1999 cisternal spacing; DECISION Pedro 2026-11-10 (Golgi around the centrosome)",
         "v6 plan R5 + v7 prep (EM blob 6.4 um from the MTOC, overlapping mito/ER)")
    nrow("(ENDO anchored)", "*", "E_band", "packer", "sh_ENDO_cyto BR + 2 nm", "BR + 8 nm", "Rab hypervariable C-terminal linker",
         "v6 pack: ENDO-cyto 0.747 placed")

    L = pd.DataFrame(log + nonint)
    for r, (m, g, ing) in moved.items(): say(f"  {r}: {len(ing)} / {len(g)} / {m:.4f}")
    say(f"mass check: {I0.mass_pg.sum():.4f} -> {I.mass_pg.sum():.4f} pg; field changes {len(log)}; non-inventory rows {len(nonint)}")
    return {"ingredients": I, "changes": normalise_changes(L, VER),
            "extra": {"changes_legacy.tsv": L, "rule_genes.tsv": pd.DataFrame(grow),
                      "geom_context.json": json.dumps(dict(pcm_volume_um3=vol_pcm, **ctx_mt, lamina_consumption=lam_path), indent=1)}}


apply_v7.OUTPUTS = ["changes_legacy.tsv", "rule_genes.tsv", "geom_context.json"]
apply_v7.PREV = ["rules/jurkat.v4/changes_legacy.tsv", "v0/inventory_v0.tsv"]
apply_v7.PARAM_FILES = ["lamina_consumption", "chromatin_json", "cytoskeleton_json", "cytoskeleton_prev_json"]
register("inventory_rule", "jurkat.v7")(apply_v7)


# ---------------------------------------------------------------------------------------------- L3 as its own rule set
def l3_apply(I, lj, Wt, setw, EV):
    """rule L3 on table I: per lamin, lamina share = lattice copies / copies, the rest of the lamina share -> nucleus"""
    for g, P_ in lj["consumed"].items():
        k = I.index[I.ing_id == g]
        if not len(k): continue
        i = k[0]; w = Wt(i); C_ = float(I.at[i, "copies"]); wl = w.get("lamina", 0.0); s_ = min(wl, P_ / C_)
        nw = dict(w); nw["lamina"] = s_; nw["nucleus"] = nw.get("nucleus", 0.0) + (wl - s_)
        EV[("L3", g)] = f"lattice holds {P_:,} of {int(round(wl * C_)):,} lamina copies (jlamina.py); rest -> nucleoplasmic lamin pool"
        setw(i, nw, "L3", {g})


def l2_rows(lj, genes=("LMNB1", "LMNB2", "LMNA")):
    """non-inventory history rows L2_lattice (one per lamin)"""
    return [dict(version=VER, gene=g, ingredient=g, rule_id="L2_lattice", field="prepass", old_value="random lamina shell (v6 0.68 placed)",
                 new_value=f"lattice pre-pass {lj['consumed'][g]:,} of {lj['target'][g]:,} (rest: L3 nucleoplasm)",
                 evidence="jlamina.py two-layer rod lattice under the INM (Turgay et al. 2017)",
                 motivating_check="v6 pack: lamina shell 0.683 placed", mass_moved_pg=0.0) for g in genes]


def apply_v7_lamina(ctx, src):
    """Rule set jurkat.v7_lamina: rule L3 of v7 (+ the L2_lattice history rows) on the v7 table, as a step of its own.

    Why a separate step (v9, 2026-11-23): L3 reads the consumption of the lamina this build makes; the lamina builder
    needs the library (lamin and chromatin structures) and the library is built from this build's species table, so L3
    inside rules:jurkat.v7 made a cycle (pack_table -> library -> lamina -> rules:jurkat.v7 -> pack_table). L3 only
    changes the comp_weights of the lamins (lamina -> nucleus), which the library does not read: the library reads the
    v7 table (inventory.library_table, after: jurkat.v7) and the pack_table step checks that the two tables differ in
    comp_weights only. No other v7 rule reads or writes the lamin weights after L3, so the result equals L3 in place."""
    p = ctx.params("jurkat.v7_lamina")
    K = data("rules_v7.yaml"); CHECK = dict(K["check"])
    I = pd.read_csv(src, sep="\t"); I0 = I.copy()
    lam_path = resolve_input(ctx, p["lamina_consumption"])
    lj = json.load(open(lam_path))
    log, EV = [], {}

    def Wt(i): return json.loads(I.at[i, "comp_weights"]) if isinstance(I.at[i, "comp_weights"], str) else {}
    def wjson(w): return json.dumps({k: round(float(v), 4) for k, v in sorted(w.items(), key=lambda kv: -kv[1]) if v >= 1e-4})

    def setw(i, nw, rule, gs):
        ow = Wt(i); d = 0.5 * sum(abs(nw.get(k, 0) - ow.get(k, 0)) for k in set(nw) | set(ow))
        if d < 1e-4: return False
        log.append(dict(version=VER, gene=I.at[i, "genes"], ingredient=I.at[i, "ing_id"], rule_id=rule, field="comp_weights",
                        old_value=I.at[i, "comp_weights"], new_value=wjson(nw),
                        evidence=" | ".join(f"{g}: {EV.get((rule, g), '')}" for g in sorted(gs)),
                        motivating_check=CHECK[rule.split("_")[0]], mass_moved_pg=round(d * I.at[i, "mass_pg"], 6)))
        I.at[i, "comp_weights"] = wjson(nw); return True

    l3_apply(I, lj, Wt, setw, EV)
    ctx.say(f"L3 lamin pool ({lam_path}): {lj['consumed']} in the lattice; rest -> nucleus")
    L = pd.DataFrame(log + l2_rows(lj))
    ctx.say(f"mass check: {I0.mass_pg.sum():.4f} -> {I.mass_pg.sum():.4f} pg; field changes {len(log)}")
    return {"ingredients": I, "changes": normalise_changes(L, VER),
            "extra": {"changes_legacy.tsv": L,
                      "lamina_context.json": json.dumps(dict(lamina_consumption=lam_path, consumed=lj["consumed"],
                                                             target=lj.get("target")), indent=1)}}


apply_v7_lamina.OUTPUTS = ["changes_legacy.tsv", "lamina_context.json"]
apply_v7_lamina.PARAM_FILES = ["lamina_consumption"]
register("inventory_rule", "jurkat.v7_lamina")(apply_v7_lamina)
