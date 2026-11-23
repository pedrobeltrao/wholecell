# origin: jurkat/inventory/inventory_v0.py (legacy workspace, 2026-11-15)
# changes: split into adapter.proteome "jurkat.geiger_ibaq_v0" (abundance -> copies) and adapter.localisation
# "jurkat.hpa_uniprot_itzhak_v0" (compartment weights); the core v0 step joins them and adds mass_pg. TOTAL_PG comes from
# config inventory.total_protein_pg (required; the legacy default 80 pg is gone). ZERO / ZERO_PATTERNS moved to
# plugins/jurkat/data/inventory_v0.yaml; HPA_MAP / COMP / MEMB to plugins/jurkat/data/localisation/*.tsv (A5); paths are frozen inputs. Numerics identical to legacy.
"""Inventory v0 for Jurkat cell 1 (PLANNING_NOTES §4-results, §4b):
  abundance  = median Geiger 2012 iBAQ ppb (E-PROT-20 linear reprocessing) over detected replicates
  imputation = expressed (E6-1 RNA call) but undetected genes: log-log fit ppb ~ TPM on detected genes, capped at the
               5th percentile of detected ppb; 'low' genes x0.25
  zeros      = RNA 'off' and undetected; manual overrides (verified losses, homozygous deletions; data file)
  scale      = sum(copies * mass) = total_protein_pg (config)
  location   = HPA main 1 / additional 0.3, blended 50:50 with UniProt CC text; Itzhak 2016 N/C re-split (soluble)."""
import os
import re

import numpy as np
import pandas as pd
import yaml

from wholecell.registry import register

AVOGADRO = 6.02214076e23
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data")


def curated(name="inventory_v0.yaml"):
    return yaml.safe_load(open(os.path.join(DATA, name)))


def _uniprot(path):
    up = pd.read_csv(path, sep="\t")
    up = up[up.Reviewed == "reviewed"].rename(columns={"Gene Names (primary)": "gene"})
    up = up.dropna(subset=["gene"]).drop_duplicates("gene")
    up["gene"] = up.gene.str.split("; ").str[0]
    return up.drop_duplicates("gene").set_index("gene")


def geiger_ibaq_v0(ctx, total_pg):
    p = ctx.params("jurkat.geiger_ibaq_v0")
    cur = curated()
    ZERO = cur["zero_genes"]; ZERO_PATTERNS = cur["zero_patterns"]
    TOTAL_PG = float(total_pg)
    up = _uniprot(ctx.frozen(p["uniprot"]))
    # proteomics: MaxQuant protein groups (linear ppb); split a group's ppb equally among its genes that share the first
    # gene's 3-letter prefix (histone/tubulin families), ignore other members
    acc2gene = dict(zip(up.Entry, up.index))
    g = pd.read_csv(ctx.frozen(p["proteome"]), sep="\t")
    cols = [c for c in g.columns if p["sample_match"] in c]
    rows = []
    for ids, vals in zip(g["Majority.protein.IDs"], g[cols].apply(pd.to_numeric, errors="coerce").values):
        if np.all(np.isnan(vals)) or np.nansum(vals) <= 0: continue
        gs = list(dict.fromkeys(acc2gene[a.split("-")[0]] for a in str(ids).split(";") if a.split("-")[0] in acc2gene))
        if not gs: continue
        fam = [x for x in gs if x[:3] == gs[0][:3]]
        for x in fam: rows.append((x, *(vals / len(fam))))
    P = pd.DataFrame(rows, columns=["gene"] + cols).groupby("gene").sum(min_count=1)
    ppb = P.median(axis=1, skipna=True)
    nrep = P.notna().sum(axis=1)
    ppb = ppb[ppb > 0]
    ctx.say("protein groups mapped:", len(rows), "genes:", len(ppb))
    rna = pd.read_csv(ctx.need("rna/rna.tsv"), sep="\t", index_col="gene")
    genes = sorted(set(up.index) & (set(ppb.index) | set(rna.index)))
    inv = pd.DataFrame(index=pd.Index(genes, name="gene"))
    inv["uniprot"] = up.Entry.reindex(inv.index)
    inv["mass_da"] = up.Mass.reindex(inv.index).astype(float)
    inv["ppb_geiger"] = ppb.reindex(inv.index)
    inv["n_rep"] = nrep.reindex(inv.index).fillna(0).astype(int)
    inv["tpm"] = rna.tpm_median.reindex(inv.index)
    inv["rna_call"] = rna.call.reindex(inv.index).fillna("no_rna")
    d = inv.dropna(subset=["ppb_geiger"]).query("tpm > 0")
    b, a = np.polyfit(np.log10(d.tpm), np.log10(d.ppb_geiger), 1)
    floor = np.percentile(inv.ppb_geiger.dropna(), p["floor_percentile"])
    fit = 10 ** (a + b * np.log10(inv.tpm.clip(lower=1e-3)))
    inv["source"] = p["source_label"]
    inv["ppb"] = inv.ppb_geiger
    und = inv.ppb.isna()
    inv.loc[und & (inv.rna_call == "expressed"), "source"] = "imputed_rna"
    inv.loc[und & (inv.rna_call == "low"), "source"] = "imputed_rna_low"
    inv.loc[und & inv.rna_call.isin(["off", "no_rna"]), "source"] = "zero_off"
    m = inv.source == "imputed_rna"; inv.loc[m, "ppb"] = np.minimum(fit[m], floor)
    m = inv.source == "imputed_rna_low"; inv.loc[m, "ppb"] = np.minimum(fit[m], floor) * p["low_factor"]
    inv["ppb"] = inv.ppb.fillna(0.0)
    zero_hit = inv.index.isin(list(ZERO)) | inv.index.to_series().str.match("|".join(ZERO_PATTERNS)).values
    inv.loc[zero_hit, ["ppb", "source"]] = [0.0, "zero_override"]
    inv["note"] = inv.index.map(ZERO).fillna("")
    inv.loc[(inv.source == p["source_label"]) & (inv.rna_call == "off"), "note"] = "detected despite RNA off"
    frac = inv.ppb / inv.ppb.sum()
    mean_mass = (frac * inv.mass_da.fillna(inv.mass_da.median())).sum()
    N = TOTAL_PG * 1e-12 / (mean_mass / AVOGADRO)
    inv["copies"] = (frac * N).round().astype(np.int64)
    ctx.say(f"impute fit: log10 ppb = {a:.2f} + {b:.2f} log10 TPM; floor {floor:.0f} ppb; total {TOTAL_PG} pg -> {N:.3g} "
            f"molecules; mean mass {mean_mass/1e3:.1f} kDa")
    return inv


geiger_ibaq_v0.PREV = ["rna/rna.tsv"]
register("adapter.proteome", "jurkat.geiger_ibaq_v0")(geiger_ibaq_v0)


def hpa_uniprot_itzhak_v0(ctx, inv):
    p = ctx.params("jurkat.hpa_uniprot_itzhak_v0")
    from wholecell.inventory.localisation.evidence import read_curated
    up = _uniprot(ctx.frozen(p["uniprot"]))
    tmap = read_curated("plugins/jurkat/data/localisation/hpa_term_map.tsv")
    HPA_MAP = {}
    for t, c in zip(tmap.term, tmap.compartment):
        HPA_MAP.setdefault(c, []).append(t)
    wmain, wadd = p["hpa_main_weight"], p["hpa_additional_weight"]
    blend = p["hpa_blend"]
    term2comp = {t: c for c, ts in HPA_MAP.items() for t in ts}
    hpa = pd.read_csv(ctx.frozen(p["hpa"]), sep="\t").drop_duplicates("Gene name").set_index("Gene name")

    def hpa_w(gname):
        if gname not in hpa.index: return None
        w = {}
        for col, wt in (("Main location", wmain), ("Additional location", wadd)):
            for t in str(hpa.at[gname, col]).split(";"):
                c = term2comp.get(t.strip())
                if c: w[c] = max(w.get(c, 0), wt)
        return w or None
    rules = read_curated("plugins/jurkat/data/localisation/uniprot_cc_rules.tsv")
    COMP = list(zip(rules.compartment, rules.regex))
    loc_txt = up["Subcellular location [CC]"].reindex(inv.index).fillna("")
    loc_txt = loc_txt.str.replace(r"\{[^}]*\}", "", regex=True).str.replace("SUBCELLULAR LOCATION:", "", regex=False).str.split("Note=").str[0]

    def up_w(t):
        hit = [c for c, pat in COMP if re.search(pat, t.lower())]
        if hit and set(hit) - {"secreted"}: hit = [c for c in hit if c != "secreted"]
        return {c: 1.0 for c in hit} or {"unassigned": 1.0}
    tm = up.Transmembrane.reindex(inv.index).notna()
    MEMB = set(read_curated("plugins/jurkat/data/localisation/membrane_keys.tsv").compartment)
    W, SRC = {}, {}
    for gname in inv.index:
        uw = up_w(loc_txt[gname]); hw = hpa_w(gname)
        if hw is not None and hpa.at[gname, "Reliability"] == "Uncertain": hw = None
        if tm[gname]:
            uw = {k: v for k, v in uw.items() if k in MEMB} or uw
            if hw is not None: hw = {k: v for k, v in hw.items() if k in MEMB} or None
        norm = lambda w: {k: v / sum(w.values()) for k, v in w.items()}
        if hw is None: w, SRC[gname] = norm(uw), "uniprot"
        elif "unassigned" in uw: w, SRC[gname] = norm(hw), "hpa"
        else:
            hw, uw = norm(hw), norm(uw); w = {k: blend * hw.get(k, 0) + (1 - blend) * uw.get(k, 0) for k in set(hw) | set(uw)}
            SRC[gname] = "hpa+uniprot"
        W[gname] = {k: round(v, 3) for k, v in w.items() if v >= 0.01}
    ITZ = pd.read_csv(ctx.frozen(p["itzhak"]), sep="\t")
    ITZ = ITZ[ITZ.primary_for_gene & (ITZ.n_good_mnc_profiles >= p["itzhak_min_profiles"]) & ITZ.mnc_class.isin(["N", "C", "N/C"])].drop_duplicates("gene").set_index("gene")
    NUC_BASE = p["itzhak_nuclear_base"]; NC = ("nucleus", "cytosol", "cytoskeleton"); n_loc = 0
    for gname in inv.index:
        if tm[gname] or gname not in ITZ.index: continue
        nc, cy = ITZ.at[gname, "frac_nuclear"], ITZ.at[gname, "frac_cytosol"]
        fn = max(nc / (nc + cy) - NUC_BASE, 0.0) / (1 - NUC_BASE)
        w = dict(W[gname]); a_ = w.get("nucleus", 0.0); S = sum(w.get(k, 0.0) for k in NC)
        if S == 0 or (a_ == 0 and fn < 0.5): continue
        cyt = {k: w[k] for k in ("cytosol", "cytoskeleton") if w.get(k, 0) > 0} or {"cytosol": 1.0}; cs = sum(cyt.values())
        for k in NC: w.pop(k, None)
        for k, v in cyt.items(): w[k] = S * (1 - fn) * v / cs
        w["nucleus"] = S * fn
        W[gname] = {k: round(v, 3) for k, v in w.items() if v >= 0.01}; SRC[gname] += "+itzhak"; n_loc += 1
    ctx.say(f"nuclear share from Itzhak 2016 fractionation: {n_loc} soluble proteins")
    out = pd.DataFrame(index=inv.index)
    out["compartments"] = [";".join(f"{k}:{v}" for k, v in sorted(W[x].items(), key=lambda kv: -kv[1])) for x in inv.index]
    out["loc_source"] = pd.Series(SRC)
    out["transmembrane"] = up.Transmembrane.reindex(inv.index).notna().values
    return out


# Superseded 2026-11-16 (unification A5) by the core localisation framework (adapter `localisation`, preset legacy_v7,
# wholecell/inventory/localisation/), which reproduces this table exactly; kept registered for the equivalence test
# (tests/test_localisation.py) and as the readable reference of the v0 arithmetic.
register("adapter.localisation", "jurkat.hpa_uniprot_itzhak_v0")(hpa_uniprot_itzhak_v0)
