"""B5 localisation common defaults: inventory-only dry run and the build verdict (criteria
docs/unification_criteria/B5_localisation.md; frozen before this script).

    python scripts/b5_localisation_check.py dryrun  --cell jurkat|yeast --ref-build <complete build> --out <dir> [--nproc 3]
    python scripts/b5_localisation_check.py verdict --cell C --build <B5 build> --ref-build <reference build>
                                                    [--switch-commit SHA] [--complete 1]

dryrun: in <out>/{legacy,common}/ (pseudo-builds: geometry/ copied from the reference build, inventory/ written here)
  the inventory with the cell's legacy preset and with preset `common` from the same code and inputs (Jurkat: every
  inventory step but `compare`; yeast: copies + v0, the localisation is per protein), then
  - mass_moved.tsv     criteria §5: combiner ablation stages (a)-(f) on the common run's own evidence table, the rule-set
                       rules of both runs (Jurkat), and legacy -> common in total
  - metrics.json       intended-level metrics available without packing (criteria §5): Jurkat C1 intended (loc_compare on
                       the pseudo-build: tiles / library of the reference build are read for the placed level, which is
                       NOT meaningful and not reported) and C3i intended; yeast C1 intended (common locmetrics with the
                       reference sets of the reference suite run) and the intended image N/C on the test split
  - report.json        everything above + top-compartment changes
verdict: criteria §4 L1-L7 -> <build>/validate/compare_B5/{verdict.json, mass_moved.tsv, c3i.tsv}

Nothing is written into a reference build (read only). CPU: --nproc worker processes (default 3).
"""
import argparse
import copy
import glob
import json
import os
import shutil
import sys
import time

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from wholecell import registry  # noqa: E402
from wholecell.config import load as load_config  # noqa: E402

LEGACY = {"jurkat": "legacy_v7", "yeast": "legacy_x5"}
REVIEW_PG, REVIEW_FRAC = 2.0, 0.02            # criteria L5
C1_SETS = {"jurkat": ["opencell", "hpa_consensus_all", "hpa_jurkat", "hpa_jurkat_holdout_half"],
           "yeast": ["SWAT_native", "SWAT_NOP1", "hyperLOPIT"]}
C1_HEAD = {"jurkat": "agree", "yeast": "tophit"}                     # criteria L1 headline metric
YEAST_BENCH = {"SWAT_native": 0.655, "hyperLOPIT": 0.838, "image_nc": 0.80}   # criteria L4 (plan §3d - 0.01)
# yeast inventory compartments -> C1 classes (origin: yeast validation/localization.py:COMP_CLASS, CLASSES @ e3b7448)
Y_COMP_CLASS = {"cytosol": "cytosol", "nucleoplasm": "nucleus", "nucleolus": "nucleolus", "ER": "ER",
                "nuclear_envelope": "ER", "mito_IM": "mitochondrion", "mito_OM": "mitochondrion", "mito_IMS": "mitochondrion",
                "mito_matrix": "mitochondrion", "mitochondrion": "mitochondrion", "golgi_endosome": "golgi_endosome",
                "vacuole": "vacuole", "plasma_membrane": "periphery", "peroxisome": "peroxisome",
                "lipid_droplet": "lipid_droplet", "cell_wall": "secreted"}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def cfg_for(cell, preset):
    ov = {"inventory": {"params": {"localisation": {"preset": preset}}}}
    cfg = load_config(cell, overrides=ov)
    registry.load_plugins(cfg)
    return cfg


def tvd(a, b):
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def top(w):
    return max(sorted(w), key=lambda k: w[k]) if w else None


def parse_cw(s):
    if not isinstance(s, str) or not s:
        return {}
    if s.startswith("{"):
        return json.loads(s)
    return {k: float(v) for k, v in (x.rsplit(":", 1) for x in s.split(";") if x)}


# ------------------------------------------------------------------------------------------- combiner ablation (§5)
class _Ctx:
    def __init__(self, bdir):
        self.bdir = bdir

    def say(self, *a):
        pass


def _restrict_order(S, o, ctx=None):
    """ablation stage (a), yeast: a transmembrane protein's ordered key list keeps only membrane keys (when it has any)"""
    keys = set(o["keys"])
    for st in S.values():
        if st.tm and st.order:
            r = [k for k in st.order if k in keys]
            if r:
                st.order = r


def combiner_options(cfg):
    from wholecell.inventory import localisation as LOC
    p = cfg["inventory"]["params"]["localisation"]
    o = LOC.options(p)
    cb = LOC._resolve_data(o["combiner"], p)
    cal = (cb.get("nc_model") or {}).get("calibration") or {}
    if cal.get("exclude_genes_input"):
        f = cfg.frozen(cal["exclude_genes_input"])
        cal["exclude_genes"] = sorted({x.strip() for x in open(f) if x.strip() and not x.startswith("#")})
    return cb


REVISION = ("chromatin", "prior", "bound")      # nc_model options of the revision 2026-11-18


def _nc_unrevised(nc):
    """the common nc_model options without the revision 2026-11-18 options (the e094024 dry-run state)"""
    o = copy.deepcopy(nc)
    for k in REVISION:
        o[k] = None
    return o


def revision_stages(e, com, background):
    """ablation stages of the revision 2026-11-18 (docs/unification_criteria/B5_localisation_revision_2026-11-18.md), each on
    top of the previous: g imaging background (+ HPA keys all supported in new_from_imaging), h Itzhak prior, i bounded
    mixtures, j chromatin-bound pools; a stage whose option the cell does not set is a no-op (reported with 0)"""
    out = []
    g = copy.deepcopy(e)
    if background and com.get("imaging_background") is not None:
        i = g["steps"].index("combine")
        g["steps"] = g["steps"][:i] + ["imaging_background"] + g["steps"][i:]
        g["imaging_background"] = copy.deepcopy(com["imaging_background"])
        g["new_from_imaging"]["background_sources"] = com["new_from_imaging"].get("background_sources")
    b = (com.get("imaging_background") or {}).get("b")
    out.append(("g", f"imaging background b = {b} on HPA (+ background-corrected HPA blended as in v7: new_from_imaging guards off for it)", g))
    h = copy.deepcopy(g)
    h["nc_model"]["prior"] = copy.deepcopy(com["nc_model"].get("prior"))
    out.append(("h", "nc_model prior: Itzhak 2016 ratio blended with the imaging ratio (log space, weight 1)", h))
    i_ = copy.deepcopy(h)
    i_["nc_model"]["bound"] = copy.deepcopy(com["nc_model"].get("bound"))
    out.append(("i", "nc_model bound: nuclear share in [0.1, 0.9] unless >= 2 exclusive sources (categorical imaging calibration only)", i_))
    j = copy.deepcopy(i_)
    j["nc_model"]["chromatin"] = copy.deepcopy(com["nc_model"].get("chromatin"))
    out.append(("j", "nc_model chromatin: histone / chromatin-bound pools stay nuclear", j))
    return out


def ablation_stages(cell, leg, com):
    """[(stage id, description, combiner options)] (criteria §5): each stage switches one default on top of the previous"""
    st = [("0", f"legacy preset {LEGACY[cell]}", leg)]
    if cell == "jurkat":
        st.append(("a", "key selection / membrane restriction: Jurkat has no key selection; membrane restriction as v7 (no-op)", leg))
        st.append(("b", "split score: v7 already score-proportional (no-op)", leg))
        c = copy.deepcopy(leg)
        c["combine"] = copy.deepcopy(com["combine"])
        st.append(("c", "combine: precedence classes (annotation + imaging 0.5 : 0.5; Jurkat groups = keys)", c))
        d = copy.deepcopy(c)
        i = d["steps"].index("combine")
        d["steps"] = d["steps"][:i + 1] + ["new_from_imaging"] + d["steps"][i + 1:]
        d["new_from_imaging"] = copy.deepcopy(com["new_from_imaging"])
        st.append(("d", "new_from_imaging: imaging-only keys need score >= 0.5; never nucleus / cytosol", d))
        e = copy.deepcopy(d)
        e["nc_model"] = _nc_unrevised(com["nc_model"])
        st.append(("e", "nc_model: HPA category ratio / r = 1, this build's volumes (replaces the Itzhak re-split)", e))
        st += revision_stages(e, com, background=True)
    else:
        a = copy.deepcopy(leg)
        i = a["steps"].index("mito_fractionation")
        a["steps"] = a["steps"][:i + 1] + ["_restrict_order"] + a["steps"][i + 1:]
        a["_restrict_order"] = {"keys": copy.deepcopy(com["restrict_membrane"]["keys"])}
        st.append(("a", "key selection kept; membrane restriction of transmembrane proteins (new for yeast)", a))
        b = copy.deepcopy(a)
        b["rank_weights"]["base"] = 1.0
        b["rank_weights"]["secondary"] = []
        b["imaging_groups"]["rank"]["base"] = 1.0
        b["imaging_groups"]["rank"]["secondary"] = []
        st.append(("b", "split score: equal shares on ties (rank 1, 1/2, 1/4 and MITO_SECONDARY 0.1 dropped)", b))
        c = copy.deepcopy(com)
        c["steps"] = [s for s in c["steps"] if s not in ("new_from_imaging", "weight_overrides")]   # curated rows: stage f
        vol = leg["rank_weights"]["pool"]["volumes"]
        c["nc_model"] = dict(_nc_unrevised(com["nc_model"]), volumes={"nucleus": vol["nucleoplasm"], "cytoplasm": vol["cytosol"]},
                             add_min=1e9, new_min=1e9)
        c["round"] = copy.deepcopy(leg["round"])
        st.append(("c", "combine: precedence classes, annotation + imaging 0.5 : 0.5 at the CYCLoPs groups (X5: imaging "
                        "set the group shares); N/C as X5 (k = 5, X5 volumes, no pool side added)", c))
        d = copy.deepcopy(c)
        i = d["steps"].index("combine")
        d["steps"] = d["steps"][:i + 1] + ["new_from_imaging"] + d["steps"][i + 1:]
        st.append(("d", "new_from_imaging (common: group named only by imaging >= 0.5, never nucleocytoplasm, mito_hc)", d))
        e = copy.deepcopy(d)
        e["nc_model"] = _nc_unrevised(com["nc_model"])
        st.append(("e", "nc_model: this build's volumes, pool sides named by imaging may be added", e))
        st += revision_stages(e, com, background=False)
    st.append(("f", "round: 4 decimals, shares < 1e-4 dropped (= preset common; + curated weight_overrides rows of the cell, "
                    "2026-11-19)", com))
    return st


def run_ablation(cell, cfg_leg, cfg_com, bdir, inv_dir):
    """weights per gene for every ablation stage, on the evidence table of inv_dir"""
    from wholecell.inventory.localisation import combiner as CB, evidence as EVM
    CB.STEPS["_restrict_order"] = _restrict_order
    E = EVM.read(os.path.join(inv_dir, "localisation", "evidence.tsv"))
    F = pd.read_csv(os.path.join(inv_dir, "localisation", "gene_features.tsv"), sep="\t", keep_default_na=False,
                    na_values={"uniprot": [""]}).set_index("gene")
    F["transmembrane"] = F.transmembrane.astype(str).str.lower().eq("true")
    for c in ("n_tm",):
        if c in F.columns:
            F[c] = pd.to_numeric(F[c], errors="coerce").fillna(0).astype(int)
    genes = list(F.index)
    leg, com = combiner_options(cfg_leg), combiner_options(cfg_com)
    out = []
    for sid, desc, o in ablation_stages(cell, leg, com):
        S = CB.run(E, F, genes, o, _Ctx(bdir))
        out.append((sid, desc, {g: dict(S[g].w) for g in genes}))
        log(f"ablation stage {sid}: done")
    return out, genes


def gene_mass(inv_dir, cell):
    V = pd.read_csv(os.path.join(inv_dir, "v0", "inventory_v0.tsv"), sep="\t")
    key = "orf" if cell == "yeast" else "gene"
    return V.groupby(V[key].astype(str)).mass_pg.sum().to_dict(), float(V.mass_pg.sum())


def mass_rows(stages, mass, total):
    rows = []
    for (s0, _, W0), (s1, d1, W1) in zip(stages[:-1], stages[1:]):
        m = n = t = 0
        for g in W1:
            x = tvd(W0[g], W1[g])
            if x > 1e-9:
                n += 1
                m += x * mass.get(g, 0.0)
                t += top(W0[g]) != top(W1[g])
        rows.append(dict(rule=f"combiner.{s1}", level="combiner (per protein)", description=d1, changed=n, top_changed=t,
                         mass_moved_pg=round(m, 6), share_of_total=round(m / total, 6)))
    return rows


def top_movers(stages, mass, n=10):
    """per ablation stage: the n proteins that move the most mass (gene, pg, weights before -> after)"""
    out = {}
    for (s0, _, W0), (s1, _, W1) in zip(stages[:-1], stages[1:]):
        x = sorted(((tvd(W0[g], W1[g]) * mass.get(g, 0.0), g) for g in W1), reverse=True)[:n]
        out[s1] = [dict(gene=g, pg=round(m, 4), before={k: round(v, 3) for k, v in W0[g].items()},
                        after={k: round(v, 3) for k, v in W1[g].items()}) for m, g in x if m > 0]
    return out


def rule_rows(inv_dir, tag):
    """mass moved per rule of the Jurkat rule sets (changes.tsv comp_weights rows x species mass)"""
    rows = []
    for rs in sorted(glob.glob(os.path.join(inv_dir, "rules", "*"))):
        ch = os.path.join(rs, "changes.tsv")
        if not os.path.exists(ch):
            continue
        C = pd.read_csv(ch, sep="\t", dtype=str, keep_default_na=False)
        I = pd.read_csv(os.path.join(rs, "ingredients.tsv"), sep="\t").drop_duplicates("ing_id").set_index("ing_id").mass_pg
        C = C[C.field == "comp_weights"]
        agg = {}
        for r in C.itertuples(index=False):
            x = tvd(parse_cw(r.old), parse_cw(r.new)) * float(I.get(r.species_id, 0.0))
            a = agg.setdefault(r.rule_id, [0.0, set(), 0])
            a[0] += x; a[1].add(r.species_id); a[2] += top(parse_cw(r.old)) != top(parse_cw(r.new))
        for k, (m, sp, t) in agg.items():
            rows.append(dict(rule=k, level=f"rule set {os.path.basename(rs)} ({tag})", description="", changed=len(sp),
                             top_changed=t, mass_moved_pg=round(m, 6)))
    return rows


def flag(rows, total):
    for r in rows:
        r["share_of_total"] = round(r["mass_moved_pg"] / total, 6)
        r["review"] = bool(r["mass_moved_pg"] > REVIEW_PG or r["mass_moved_pg"] > REVIEW_FRAC * total)
    return rows


# ------------------------------------------------------------------------------------------------ metrics (intended)
def jurkat_c3i(cfg, loc_model, measures):
    """criteria L3 C3i: per Tier-2 gene log2 N/C from loc_model.parquet (intended and placed) vs the HPA image nc_hpa"""
    from scipy.stats import spearmanr
    c = cfg["validate"]["suite"]["image_nc"]
    T = pd.read_csv(c["table"], sep="\t")[["gene", c["measured"], c["split"]]].dropna()
    M = pd.read_parquet(loc_model, columns=["gene", "region", "intended", "placed"])
    v = measures["region_volumes_mc"]["volume_um3"]
    vn = v["nucleoplasm"] + v["nucleolus"] + v.get("nuclear_speck", 0.0)
    vc = v["cytosol"]
    nuc = M.region.isin(["nucleoplasm", "nucleolus", "speckle"])
    out, per = {}, {}
    for lvl in ("intended", "placed"):
        n = M[nuc].groupby("gene")[lvl].sum()
        cy = M[M.region == "cytosol"].groupby("gene")[lvl].sum()
        d = pd.DataFrame(dict(n=n, c=cy)).fillna(0.0)
        lo = 1e-6
        r = np.log2(np.maximum(d.n / vn, lo) / np.maximum(d.c / vc, lo)).clip(-8, 8)
        per[lvl] = r
        X = T.set_index("gene").join(r.rename("model"), how="inner").dropna()
        for sp in ("holdout", "all"):
            Y = X if sp == "all" else X[X[c["split"]] == sp]
            out[f"c3i.{lvl}.{sp}.spearman"] = float(spearmanr(Y[c["measured"]], Y.model).correlation) if len(Y) > 2 else None
            out[f"c3i.{lvl}.{sp}.n"] = int(len(Y))
    return out, pd.DataFrame(per)


def jurkat_loc_compare(cfg, pseudo, nproc):
    """loc_compare (model / compare / holdout) on a pseudo-build; returns loc_holdout.json"""
    from wholecell.validate import common, loc_compare as LC
    inp = common.Inputs(cfg, pseudo)
    C = LC.Ctx(cfg, pseudo, inp)
    LC.model(C, False, nproc=nproc)
    LC.compare(C, False)
    return LC.holdout(C, False)


def yeast_metrics(cell_cw, ref_inputs, image_tsv, measures):
    """yeast C1 intended (common locmetrics, reference sets of the reference suite run, test split) and the intended
    image N/C on the test split (model: (nucleoplasm + nucleolus) / V vs cytosol / V of the intended shares)"""
    from scipy.stats import spearmanr
    from wholecell.validate import locmetrics
    from wholecell.validate.suite import image_nc_metrics
    comps = json.load(open(os.path.join(ref_inputs, "classes.json")))["classes"]
    refs = json.load(open(os.path.join(ref_inputs, "refs_test.json")))
    Ms = pd.read_parquet(os.path.join(ref_inputs, "model_shares.parquet"))
    mass = (Ms["copies"] * Ms["mass_da"]).to_dict()
    model = {}
    for orf, w in cell_cw.items():
        s = {}
        for k, v in w.items():
            c = Y_COMP_CLASS.get(k)
            if c:
                s[c] = s.get(c, 0.0) + v
        z = sum(s.values()) or 1.0
        model[orf] = ({c: v / z for c, v in s.items()}, float(mass.get(orf, 0.0)))
    out = {}
    for name, ref in refs.items():
        summ, _ = locmetrics.score_set(model, ref, comps)
        for m in ("agree", "agree_tvd", "tophit"):
            for w in ("unweighted", "mass_weighted"):
                out[f"{name}.intended.{m}.{w}"] = summ[m][w]
        out[f"{name}.n_genes"] = summ["n_genes"]
    v = measures["region_volumes_mc"]["volume_um3"]
    vn, vc = v["nucleoplasm"] + v["nucleolus"], v["cytosol"]
    D = pd.read_csv(image_tsv, sep="\t", index_col=0)
    rows = []
    for orf, w in cell_cw.items():
        n = w.get("nucleoplasm", 0.0) + w.get("nucleolus", 0.0)
        c = w.get("cytosol", 0.0)
        if n + c >= 0.9 * sum(w.values()) and orf in D.index:
            lo = 2.0 ** -8 * 1e-3
            rows.append((orf, float(np.clip(np.log2(max(n / vn, lo) / max(c / vc, lo)), -8, 8))))
    R = pd.DataFrame(rows, columns=["orf", "model_int"]).set_index("orf").join(D, how="inner")
    for sp in ("test", "dev"):
        X = R[R.split == sp]
        r = image_nc_metrics(X.measured_log2.values, X.model_int.values)
        out[f"image_nc_intended.{sp}.spearman"] = r["spearman"]
        out[f"image_nc_intended.{sp}.auroc_nuclear"] = r["auroc_nuclear"]
        out[f"image_nc_intended.{sp}.n"] = r["n"]
    return out


# ------------------------------------------------------------------------------- revision 2026-11-18 report items
WATCH = {"jurkat": ["TKT", "PGAM1", "H2BC12", "H2BC17", "H2BC9", "H1-2", "H1-4", "H1-0", "H2AZ1", "H3-3A", "MACROH2A1", "HMGN2", "HMGB1"],
         "yeast": ["HTA1", "HTB1", "HHT1", "HHF1", "HTZ1", "HHO1"]}


def _species_shares(inv_dir, gene):
    """the final intended comp_weights of the single-gene species of a gene (pack table; None if absent)"""
    P = pd.read_csv(os.path.join(inv_dir, "pack_table.tsv"), sep="\t", usecols=["ing_id", "genes", "comp_weights", "mass_pg"])
    s = P[P.genes.astype(str) == gene]
    if not len(s):
        return None
    s = s.loc[s.mass_pg.idxmax()]
    return {k: round(v, 4) for k, v in parse_cw(s.comp_weights).items()}


def revision_report(cell, inv, stages, genes, mass, out):
    """report items of the revision (docs/unification_criteria/B5_localisation_revision_2026-11-18.md §5): watched genes'
    shares legacy -> common (v0 combiner and, Jurkat, the final species table), intended nucleoplasm / cytosol protein
    density (Jurkat: soluble nuclear / cytoplasmic keys of the pack table over Vnucleoplasm / Vcytosol), the top 20
    movers by mass legacy -> common (v0), proteins whose top compartment changes (top_changes.tsv)"""
    W0, Wf = stages[0][2], stages[-1][2]
    key = "orf" if cell == "yeast" else "gene"
    V = pd.read_csv(os.path.join(inv["common"], "v0", "inventory_v0.tsv"), sep="\t")
    name = dict(zip(V[key].astype(str), V["gene_name" if "gene_name" in V.columns else key].astype(str)))
    byname = {v: k for k, v in name.items()}
    rep = {"watch": {}}
    for g in WATCH[cell]:
        k = byname.get(g, g)
        if k not in W0:
            continue
        r = {"v0_legacy": {x: round(v, 4) for x, v in W0[k].items()}, "v0_common": {x: round(v, 4) for x, v in Wf[k].items()},
             "mass_pg": round(mass.get(k, 0.0), 5)}
        if cell == "jurkat":
            r["final_legacy"], r["final_common"] = _species_shares(inv["legacy"], g), _species_shares(inv["common"], g)
        rep["watch"][g] = r
    if cell == "jurkat":
        from plugins.jurkat.inventory.common import soluble_masses
        meas = json.load(open(os.path.join(out, "common", "geometry", "measures.json")))["region_volumes_mc"]["volume_um3"]
        rep["density_mg_ml"] = {}
        for t in ("legacy", "common"):
            P = pd.read_csv(os.path.join(inv[t], "pack_table.tsv"), sep="\t").drop_duplicates("ing_id")
            mn, mc = soluble_masses(P)
            rep["density_mg_ml"][t] = {"nucleoplasm": round(mn / meas["nucleoplasm"] * 1e3, 2),
                                       "cytosol": round(mc / meas["cytosol"] * 1e3, 2), "nuclear_pg": round(mn, 3),
                                       "cytoplasmic_pg": round(mc, 3)}
    mv = sorted(((tvd(W0[g], Wf[g]) * mass.get(g, 0.0), g) for g in genes), reverse=True)[:20]
    rep["top20_movers_v0"] = [dict(gene=name.get(g, g), pg=round(m, 4), legacy={k: round(v, 3) for k, v in W0[g].items()},
                                   common={k: round(v, 3) for k, v in Wf[g].items()}) for m, g in mv]
    rows = [dict(gene=name.get(g, g), mass_pg=mass.get(g, 0.0), top_legacy=top(W0[g]), top_common=top(Wf[g]))
            for g in genes if top(W0[g]) != top(Wf[g])]
    T = pd.DataFrame(rows, columns=["gene", "mass_pg", "top_legacy", "top_common"]).sort_values("mass_pg", ascending=False)
    T.to_csv(os.path.join(out, "top_changes.tsv"), sep="\t", index=False)
    rep["top_changed"] = {"n": int(len(T)), "mass_pg": round(float(T.mass_pg.sum()), 4),
                          "by_transition": {f"{a}->{b}": int(n) for (a, b), n in
                                            T.groupby(["top_legacy", "top_common"]).size().sort_values(ascending=False).head(12).items()}}
    return rep


# ------------------------------------------------------------------------------------------------------------ dryrun
def make_pseudo(d, ref, link_pack):
    os.makedirs(os.path.join(d, "geometry"), exist_ok=True)
    for f in ("measures.json", "package.json"):
        shutil.copyfile(os.path.join(ref, "geometry", f), os.path.join(d, "geometry", f))
    if link_pack:
        for x in ("pack", "library", "preplaced"):
            if os.path.exists(os.path.join(ref, x)) and not os.path.exists(os.path.join(d, x)):
                os.symlink(os.path.join(ref, x), os.path.join(d, x))
        cdir = os.path.join(d, "validate", "cache")
        os.makedirs(cdir, exist_ok=True)
        for f in glob.glob(os.path.join(ref, "validate", "cache", "*.npz")):
            if not os.path.exists(os.path.join(cdir, os.path.basename(f))):
                shutil.copyfile(f, os.path.join(cdir, os.path.basename(f)))


def run_inventory(cfg, d, cell):
    from wholecell.inventory import stage as INV
    plan = INV.steps(cfg)
    todo = [s for s in plan if s != "compare"] if cell == "jurkat" else [s for s in plan if s in ("rna", "copies", "v0")]
    for s in todo:
        t = time.time()
        INV.run(cfg, d, step=s)
        log(f"inventory {s}: {time.time() - t:.0f} s")


def cw_table(inv_dir, cell):
    V = pd.read_csv(os.path.join(inv_dir, "v0", "inventory_v0.tsv"), sep="\t")
    key = "orf" if cell == "yeast" else "gene"
    return {str(k): parse_cw(c) for k, c in zip(V[key], V.compartments)}


def dryrun(a):
    T0 = time.time()
    cell, out = a.cell, a.out
    os.makedirs(out, exist_ok=True)
    cfgs = {"legacy": cfg_for(cell, LEGACY[cell]), "common": cfg_for(cell, "common")}
    for tag, cfg in cfgs.items():
        d = os.path.join(out, tag)
        make_pseudo(d, a.ref_build, cell == "jurkat")
        log(f"{tag}: inventory ({cfg['inventory']['params']['localisation']['preset']}) -> {d}")
        run_inventory(cfg, d, cell)
    inv = {t: os.path.join(out, t, "inventory") for t in cfgs}
    mass, total = gene_mass(inv["common"], cell)
    stages, genes = run_ablation(cell, cfgs["legacy"], cfgs["common"], os.path.join(out, "common"), inv["common"])
    # sanity: stage 0 = the legacy run's v0 weights, stage f = the common run's (same evidence, same code)
    chk = {}
    for t, (sid, _, W) in (("legacy", stages[0]), ("common", stages[-1])):
        V = cw_table(inv[t], cell)
        chk[t] = max(tvd({k: round(v, 4) for k, v in W[g].items()}, V.get(g, {})) for g in genes)
    rows = mass_rows(stages, mass, total)
    movers = top_movers(stages, mass)
    W0, Wf = stages[0][2], stages[-1][2]
    tot = dict(rule="combiner.total", level="combiner (per protein)", description="legacy preset -> preset common (v0)",
               changed=sum(tvd(W0[g], Wf[g]) > 1e-9 for g in genes), top_changed=sum(top(W0[g]) != top(Wf[g]) for g in genes),
               mass_moved_pg=round(sum(tvd(W0[g], Wf[g]) * mass.get(g, 0.0) for g in genes), 6))
    rows.append(tot)
    metrics = {}
    if cell == "jurkat":
        rows += rule_rows(inv["legacy"], "legacy_v7") + rule_rows(inv["common"], "common")
        P = {t: pd.read_csv(os.path.join(inv[t], "pack_table.tsv"), sep="\t").drop_duplicates("ing_id").set_index("ing_id")
             for t in inv}
        common_ids = P["legacy"].index.intersection(P["common"].index)
        sp_tvd = [(tvd(parse_cw(P["legacy"].at[i, "comp_weights"]), parse_cw(P["common"].at[i, "comp_weights"])),
                   top(parse_cw(P["legacy"].at[i, "comp_weights"])) != top(parse_cw(P["common"].at[i, "comp_weights"])),
                   float(P["common"].at[i, "mass_pg"])) for i in common_ids]
        rows.append(dict(rule="pack_table.total", level="species (pack table)", description="legacy_v7 -> common, every step",
                         changed=sum(x > 1e-9 for x, _, _ in sp_tvd), top_changed=sum(t for _, t, _ in sp_tvd),
                         mass_moved_pg=round(sum(x * m for x, _, m in sp_tvd), 6)))
        meas = json.load(open(os.path.join(out, "common", "geometry", "measures.json")))
        for t, cfg in cfgs.items():
            d = os.path.join(out, t)
            log(f"{t}: loc_compare (intended level; nproc {a.nproc})")
            H = jurkat_loc_compare(cfg, d, a.nproc)
            m = {}
            for s in C1_SETS[cell]:
                if s in H["sets"]:
                    for k in ("agree", "agree_tvd", "tophit"):
                        for w in ("unweighted", "mass_weighted"):
                            m[f"{s}.intended.{k}.{w}"] = H["sets"][s]["intended"][k][w]
                    m[f"{s}.n_genes"] = H["sets"][s]["n_genes"]
            c3, per = jurkat_c3i(cfg, os.path.join(d, "validate", "loc_model.parquet"), meas)
            per.to_csv(os.path.join(out, f"c3i_{t}.tsv"), sep="\t")
            m.update({k: v for k, v in c3.items() if ".intended." in k})
            metrics[t] = m
    else:
        ref_inputs = a.ref_inputs
        image = a.image_nc
        meas = json.load(open(os.path.join(out, "common", "geometry", "measures.json")))
        for t in cfgs:
            metrics[t] = yeast_metrics(cw_table(inv[t], cell), ref_inputs, image, meas)
    flag(rows, total)
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(out, "mass_moved.tsv"), sep="\t", index=False)
    delta = {k: (round(metrics["common"][k] - metrics["legacy"][k], 4)
                 if isinstance(metrics["common"].get(k), (int, float)) and isinstance(metrics["legacy"].get(k), (int, float)) else None)
             for k in metrics.get("common", {})}
    revision = revision_report(cell, inv, stages, genes, mass, out)
    rep = dict(cell=cell, ref_build=a.ref_build, code=REPO, commit=os.popen(f"git -C {REPO} rev-parse HEAD").read().strip(),
               total_protein_pg=total, ablation_sanity_max_tvd=chk, metrics=metrics, delta_common_minus_legacy=delta,
               mass_moved=rows, review=[r["rule"] for r in rows if r["review"] and not r["rule"].endswith(".total")],
               top_movers=movers, revision_2026_11_18=revision,
               seconds=round(time.time() - T0, 1),
               notes=["placed-level metrics need packing (not computed)",
                      "Jurkat: loc_compare placed level reads the reference build's tiles and is not meaningful here"])
    json.dump(rep, open(os.path.join(out, "report.json"), "w"), indent=1, default=float)
    log(json.dumps(dict(metrics=metrics, delta=delta), indent=1, default=float))
    log(R.to_string())


# ----------------------------------------------------------------------------------------------------------- verdict
def _suite(b):
    p = os.path.join(b, "validate", "suite.json")
    return {r["metric"]: r for r in json.load(open(p))["metrics"]} if os.path.exists(p) else {}


def requested_density_change(cell, B, R, out, nproc):
    """{region: (requested mass in the region, B5 build - reference) / region volume} in mg/ml: requested instances per
    (library, region) over all tiles (validate cache tile_sums, else counted here into <out>) x library mass (libs.json
    mass_da, else the library bead mass), volumes from the B5 build's densities.json"""
    from wholecell.validate import common
    cfg = load_config(cell)
    res = {}
    for tag, b in (("build", B), ("reference", R)):
        common.CACHE_DIR = os.path.join(out, f"cache_{tag}")
        os.makedirs(common.CACHE_DIR, exist_ok=True)
        inp = common.Inputs(cfg, b)
        own = os.path.join(b, "validate", "cache", "tile_sums.npz")
        T = common.tile_sums(inp, own if os.path.exists(own) else os.path.join(common.CACHE_DIR, "tile_sums.npz"), nproc=nproc)
        L = json.load(open(inp["libs"]))
        m, _ = common.lib_masses(inp, b, nproc=nproc)
        rn = L.get("rnames") or []
        pg = (T["req"][:len(m)] * m[:, None]).sum(0) / 6.02214076e23 * 1e12
        res[tag] = {n: float(pg[i]) for i, n in enumerate(rn)}
    common.CACHE_DIR = None
    D = os.path.join(B, "validate", "densities.json")
    vol = {r["region"]: r["volume_um3"] for r in json.load(open(D))["rows"]} if os.path.exists(D) else {}
    return {k: (res["build"][k] - res["reference"].get(k, 0.0)) / vol[k] * 1e3 for k in res["build"] if vol.get(k)}


def verdict(a):
    from wholecell.validate.suite import TOL
    cell, B, R = a.cell, a.build, a.ref_build
    out = a.out or os.path.join(B, "validate", "compare_B5")      # --out: a test run that must not write into the build
    os.makedirs(out, exist_ok=True)
    S, SR = _suite(B), _suite(R)
    items = []

    def item(cid, name, value, ref, ok, note=""):
        items.append(dict(id=cid, check=name, value=value, reference=ref, status="pass" if ok else "FAIL", note=note))

    # L1 / L2
    for lvl in ("intended", "placed"):
        for s in C1_SETS[cell]:
            k = f"{s}.{lvl}.{C1_HEAD[cell]}.unweighted"
            v, r = (S.get(k) or {}).get("value"), (SR.get(k) or {}).get("value")
            item("L1", k, v, r, v is not None and r is not None and v >= r - 0.005, "C1 >= reference - 0.005")
    for k, r in SR.items():
        if r["check"] != "C1" or r["direction"] != "up" or k.endswith(f"{C1_HEAD[cell]}.unweighted") and k.split(".")[0] in C1_SETS[cell]:
            continue
        v = (S.get(k) or {}).get("value")
        if r.get("value") is None:
            continue
        item("L2", k, v, r["value"], v is not None and v >= r["value"] - 0.01, "no-regression +-0.01")
    # L3 / L4
    if cell == "yeast":
        for m in ("image_nc.test.spearman", "image_nc.test.auroc_nuclear"):
            v, r = (S.get(m) or {}).get("value"), (SR.get(m) or {}).get("value")
            item("L3", m, v, r, v is not None and r is not None and v >= r - 0.01, "C3 >= reference - 0.01")
        for s, thr in (("SWAT_native", YEAST_BENCH["SWAT_native"]), ("hyperLOPIT", YEAST_BENCH["hyperLOPIT"])):
            k = f"{s}.intended.tophit.unweighted"
            v = (S.get(k) or {}).get("value")
            item("L4", k, v, thr, v is not None and v >= thr, "plan §3d benchmark")
        v = (S.get("image_nc.test.spearman") or {}).get("value")
        item("L4", "image_nc.test.spearman", v, YEAST_BENCH["image_nc"], v is not None and v >= YEAST_BENCH["image_nc"], "plan §3d")
    else:
        cfg = load_config(cell)
        res = {}
        for tag, b in (("build", B), ("reference", R)):
            meas = json.load(open(os.path.join(b, "geometry", "measures.json")))
            lm = os.path.join(b, "validate", "loc_model.parquet")
            res[tag] = jurkat_c3i(cfg, lm, meas)[0] if os.path.exists(lm) else {}
        json.dump(res, open(os.path.join(out, "c3i.json"), "w"), indent=1)
        for lvl in ("intended", "placed"):
            k = f"c3i.{lvl}.holdout.spearman"
            v, r = res["build"].get(k), res["reference"].get(k)
            item("L3", k, v, r, v is not None and r is not None and v >= r - 0.01, "C3i >= reference - 0.01")
    # L5 mass moved (build's own evidence) + rule sets
    code = os.path.join(B, "_code")
    inv = os.path.join(B, "inventory")
    cfg_c, cfg_l = cfg_for(cell, "common"), cfg_for(cell, LEGACY[cell])
    mass, total = gene_mass(inv, cell)
    stages, genes = run_ablation(cell, cfg_l, cfg_c, B, inv)
    rows = mass_rows(stages, mass, total) + (rule_rows(inv, "common") if cell == "jurkat" else [])
    flag(rows, total)
    json.dump(top_movers(stages, mass), open(os.path.join(out, "top_movers.json"), "w"), indent=1)
    pd.DataFrame(rows).to_csv(os.path.join(out, "mass_moved.tsv"), sep="\t", index=False)
    review = [r["rule"] for r in rows if r["review"]]
    item("L5", "mass moved per rule", max(r["mass_moved_pg"] for r in rows), f"<= {REVIEW_PG} pg and <= {REVIEW_FRAC:.0%} of {total:.2f} pg",
         True, f"review: {review}" if review else "nothing to review")
    # top-compartment changes vs the reference build
    Wb, Wr = cw_table(inv, cell), cw_table(os.path.join(R, "inventory"), cell)
    tc = sum(top(Wb[g]) != top(Wr[g]) for g in Wb if g in Wr and Wb[g] and Wr[g])
    # L6 other checks: suite metrics C5-C17 against the reference build (plan s4 rule; a C11 density change is accepted
    # when the requested (intended) mass change of the region explains it: |d rho - d rho_intended| <= 2 mg/ml)
    dens_int = {}
    for k, r in SR.items():
        if r["check"] in ("C1", "C2", "C3", "C4") or r.get("value") is None:
            continue
        d = r["direction"]
        if not (d in ("up", "down", "stable") or str(d).startswith("target:")):
            continue
        v = (S.get(k) or {}).get("value")
        if v is None:
            item("L6", k, None, r["value"], False, "metric missing in the B5 build")
            continue
        tol = TOL[r["tol_kind"]][0](r["value"])
        if tol is None:
            continue
        if d == "up":
            bad = v < r["value"] - tol
        elif d == "down":
            bad = v > r["value"] + tol
        elif d == "stable":
            bad = abs(v - r["value"]) > tol
        else:
            x = float(str(d).split(":", 1)[1])
            bad = abs(v - x) > abs(r["value"] - x) + tol
        note = f"tolerance {tol} ({d})"
        if bad and r["check"] == "C11" and k.startswith("density."):
            if not dens_int:
                dens_int.update(requested_density_change(cell, B, R, out, a.nproc))
            reg = k.split(".", 1)[1]
            di = dens_int.get(reg)
            if di is not None:
                bad = abs((v - r["value"]) - di) > 2.0
                note += f"; requested-mass density change {di:+.2f} mg/ml explains it: {not bad} (criteria L6)"
            else:
                note += "; no requested-mass density for this region (not explained)"
        item("L6", k, v, r["value"], not bad, note)
    item("L7", "build complete", a.complete, 1, str(a.complete) == "1", "snakemake 100 %")
    st = {i["status"] for i in items}
    v = "NOT ACCEPTED" if "FAIL" in st else ("NEEDS REVIEW" if review else "ACCEPTED")
    J = dict(step="B5", cell=cell, build=B, reference=R, switch_commit=a.switch_commit, verdict=v, review_rules=review,
             top_compartment_changes_vs_reference=int(tc), total_protein_pg=total,
             criteria="docs/unification_criteria/B5_localisation.md", items=items, date=time.strftime("%Y-%m-%d %H:%M"))
    json.dump(J, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=float)
    log(f"B5 {cell}: {v}; {sum(i['status'] == 'FAIL' for i in items)} failing items; review {review}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["dryrun", "verdict"])
    ap.add_argument("--cell", required=True, choices=sorted(LEGACY))
    ap.add_argument("--ref-build", required=True)
    ap.add_argument("--out")
    ap.add_argument("--build")
    ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--ref-inputs", help="yeast: <reference suite run>/adapter/loc_inputs (reference sets, masses)")
    ap.add_argument("--image-nc", help="yeast: <reference build>/engine/work/validation/image_nc.tsv")
    ap.add_argument("--switch-commit", default="")
    ap.add_argument("--complete", default="1")
    a = ap.parse_args()
    from wholecell.validate import common
    common.MAX_PROCS = a.nproc
    dryrun(a) if a.mode == "dryrun" else verdict(a)


if __name__ == "__main__":
    main()
