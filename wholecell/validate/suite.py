"""Common validation suite (unification plan section 4, step A8): one metric definition per check C1-C17.

    python -m wholecell validate --cell C --build B --step suite [--out DIR] [ref_suite=<suite.json>] [strict=1]

What it does. Every check id has ONE function in CHECKS (a plain dict) that returns metric rows for a build. A check
reads the outputs of the core validation steps of this build (<build>/validate/*.json, written by the other steps of
this stage), computes a metric itself from the build's tiles (C10 placed fractions, the seam index), or, for a cell
whose validation still runs in its own code, asks the cell's reference adapter (registry kind `validate_adapter`,
named by config `validate.suite.adapter`; the yeast one is plugins/yeast/validate_adapter.py). Phase 2: where the plan
defines one metric, the suite computes it for every cell (C1 locmetrics.py, C3, C6, C10 here) from per-cell inputs that
the adapter supplies (model shares, reference sets, image and copy-number tables), and runs core steps the build lacks
(placement, orientation, crowding, densities: run_core_steps, into <out>/core). Rows still read verbatim through an
adapter carry `source: adapter (yeast code)`. A check with no data for this cell writes a row with status `n/a` and the
reason; nothing is skipped silently.

Outputs (default <build>/validate/, or --out DIR; the suite never writes into the build when --out is given):
  suite.json         per metric: check, metric, value, unit, reference (external reference value, when one exists),
                     direction (up | down | target:<x> | report), tolerance (no-regression tolerance), tolerance_rule,
                     status (ok | regress | improve | n/a), reason, source, file, legacy (the number the legacy code wrote
                     for this build, where defined), legacy_source, doc (anchor in docs/stages/validation.md)
  no_regression.tsv  every metric against the declared reference build's suite.json (config validate.suite.reference
                     {build, suite} or ref_suite=); rule of plan section 4 (TOL below)
  legacy_check.tsv   suite value vs legacy value per metric where a legacy number exists (the A8 code check)

No-regression rule (plan section 4): a metric regresses if it moves in its bad direction by more than its tolerance
relative to the reference build. Tolerances: localisation / AUROC / Spearman +-0.01; rates max(10 % rel., 1 per 1e6);
densities +-2 mg/ml; NN median +-0.5 nm; placed fraction +-0.003. Metrics the plan gives no tolerance for use the
suite defaults in TOL (marked "suite default"); `report` metrics are never regressions.

Seam index (C10, new in A8): tiles packed in parallel meet at faces; a packer artefact shows as a density step at the
faces. Definition: tiles whose volume is >= SEAM_MIN_FRAC in the seam region (config validate.suite.seam_region,
default cytosol; tile volumes from pack/tiles.parquet); a seeded sample of <= SEAM_N_TILES of them; packed instance
centres (all packed instances of the tile, count, not mass) per nm3 in the face band (centre < 10 nm from the nearest
tile face) divided by the same density in the tile interior (centre >= 20 nm from every face). 1 = no seam; > 1
over-filled faces; < 1 depleted faces. Reported overall and as the min/max over the 8 tile parity classes (the
checkerboard passes of a packer that packs in passes).
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pandas as pd

from . import common

DOC = "docs/stages/validation.md"
SEAM_BAND_NM = 10.0          # face band: centre closer than this to the nearest tile face
SEAM_INTERIOR_NM = 20.0      # interior: centre at least this far from every face
SEAM_MIN_FRAC = 0.95         # tile volume share of the seam region (Jurkat has no tile >= 0.98 cytosol: virtual regions)
SEAM_N_TILES = 300           # seeded tile sample
MIN_REQ = 100                # species with fewer requested copies are not counted in 'species below 50 % placed'
AVOGADRO = 6.02214076e23

# tolerance kinds: name -> (function of the reference value, rule text). Plan section 4 unless marked suite default.
TOL = {
    "loc": (lambda r: 0.01, "+-0.01 (plan s4: localisation/AUROC/Spearman)"),
    "rate1e6": (lambda r: max(0.1 * abs(r), 1.0), "max(10 % rel., 1 per 1e6) (plan s4: rates)"),
    "ratefrac": (lambda r: max(0.1 * abs(r), 1e-6), "max(10 % rel., 1e-6) (plan s4: rates, as a fraction)"),
    "density": (lambda r: 2.0, "+-2 mg/ml (plan s4: densities)"),
    "nn": (lambda r: 0.5, "+-0.5 nm (plan s4: NN median)"),
    "placed": (lambda r: 0.003, "+-0.003 (plan s4: placed fraction)"),
    "phi": (lambda r: 0.005, "+-0.005 (suite default; v7 criteria V2)"),
    "ratio": (lambda r: 0.01, "+-0.01 (suite default; v7 criteria V1 N/C ratio)"),
    "anchor": (lambda r: 0.05, "+-0.05 in |ratio - 1| (suite default)"),
    "seam": (lambda r: 0.02, "+-0.02 in |index - 1| (suite default)"),
    "rel10": (lambda r: 0.1 * abs(r), "10 % rel. (suite default: concentrations)"),
    "deg": (lambda r: 0.5, "+-0.5 degree (suite default: tilt)"),
    "count0": (lambda r: 0.0, "0 (must)"),
    "report": (lambda r: None, "report only"),
}
STATUS = ("ok", "regress", "improve", "n/a")


# ------------------------------------------------------------------ metric rows
def row(check, metric, value, *, unit="", direction="report", tol="report", reference=None, reference_source="",
        source="core", file="", legacy=None, legacy_source="", legacy_note="", reason="", note=""):
    """one metric row of suite.json (status is set later by the no-regression comparison)"""
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        value = None
    return dict(check=check, metric=metric, value=value, unit=unit, reference=reference, reference_source=reference_source,
                direction=direction, tol_kind=tol, tolerance=None, tolerance_rule=TOL[tol][1], status=None,
                reason=reason, source=source, file=file, legacy=legacy, legacy_source=legacy_source, legacy_note=legacy_note, note=note,
                doc=f"{DOC}#{check.lower()}")


def na(check, metric, reason, **kw):
    r = row(check, metric, None, **kw)
    r["reason"] = reason
    return r


SRC_CORE = "core"
SRC_CORE_SUITE = "core (step run by the suite on this build, A8 phase 2)"
SRC_NEW = "suite (new, wholecell/validate/suite.py)"
SRC_COMMON = "common (wholecell/validate/locmetrics.py)"
SRC_SUITE_COMMON = "common (wholecell/validate/suite.py)"
SRC_YEAST = "adapter (yeast code)"
SRC_JURKAT_LEGACY = "adapter (legacy Jurkat code)"


class Ctx:
    """one suite run: build dir, config, the core validation folder of the build, optional legacy folder and adapter"""

    def __init__(self, cfg, bdir, out, **kw):
        self.cfg, self.bdir, self.out, self.kw = cfg, bdir, out, kw
        self.vc = (cfg.get("validate") or {})
        self.sc = (self.vc.get("suite") or {})
        self.vdir = kw.get("validate_dir") or os.path.join(bdir, "validate")
        cr = self.vc.get("compare_reference") or {}
        # legacy validation folder (Jurkat: the ported validation run on the archived legacy v7 pack, itself code-checked
        # against the stored legacy files; code_check_vs_stored.tsv there)
        self.legacy_dir = kw.get("legacy_dir") or self.sc.get("legacy_validate_dir") or cr.get("validate_dir")
        self.adapter = None
        an = self.sc.get("adapter")
        if an:
            from .. import registry
            self.adapter = registry.get("validate_adapter", an)(self)
        self._json = {}

        self.core_dir = os.path.join(out, "core")      # core steps the suite ran itself (never inside the build)

    # core validation outputs of this build / of the legacy run; then those the suite ran into <out>/core
    def vjson(self, name, legacy=False):
        ds = [self.legacy_dir] if legacy else [self.vdir, self.core_dir]
        for d in ds:
            if not d:
                continue
            p = os.path.join(d, name)
            if p not in self._json:
                self._json[p] = json.load(open(p)) if os.path.exists(p) else None
            if self._json[p] is not None:
                return self._json[p]
        return None

    def vsource(self, name):
        """'core' if the build's own validate/ has the file, 'core (run by the suite)' if the suite ran the step"""
        return SRC_CORE if os.path.exists(os.path.join(self.vdir, name)) else SRC_CORE_SUITE

    def vfile(self, name):
        return os.path.join(self.vdir, name)

    def skipped(self, step):
        p = os.path.join(self.vdir, f"{step}.skipped.json")
        return json.load(open(p)).get("reason") if os.path.exists(p) else None

    def pack(self, name):
        return os.path.join(self.bdir, "pack", name)

    def adapt(self, check):
        """rows of the cell's reference adapter for this check, or None"""
        if self.adapter is None:
            return None
        f = self.adapter.get(check)
        return f() if f else None


def _get(d, *keys):
    for k in keys:
        if d is None:
            return None
        if isinstance(d, list):
            d = d[k] if isinstance(k, int) and k < len(d) else None
        else:
            d = d.get(k)
    return d


def _why_missing(ctx, step, fname):
    s = ctx.skipped(step)
    if s:
        return f"n/a (no data): validate step '{step}' skipped: {s}"
    return f"n/a (no data): {fname} not in {ctx.vdir}"


# ------------------------------------------------------------------ core steps run by the suite (A8 phase 2)
CORE_STEPS = ("placement", "orientation", "crowding", "densities")


def run_core_steps(ctx, steps):
    """run core validation steps that this build's validate/ lacks, into <out>/core/ (caches into <out>/cache/), with
    the inputs the cell's adapter names (core_inputs(): library folder, geometry package, ...). A reference build is
    read-only: nothing is written into it. 'orientation' = membrane signs + tilt on the placement tile sample, for a
    placement.json written before the tilt metric existed."""
    from . import placement as PLM
    os.makedirs(ctx.core_dir, exist_ok=True)
    extra = ctx.adapter.core_inputs() if ctx.adapter is not None and hasattr(ctx.adapter, "core_inputs") else {}
    inp = common.Inputs(ctx.cfg, ctx.bdir, **extra)
    nproc = min(3, common.nprocs(ctx.cfg, ctx.kw))
    done = {}
    for st in steps:
        f = {"orientation": "orientation.json"}.get(st, f"{st}.json")
        if ctx.vjson(f) is not None and not (st == "orientation" and "tilt_median_deg" not in (ctx.vjson(f) or {})):
            done[st] = "present"
            continue
        if st == "orientation":
            P = ctx.vjson("placement.json")
            if P is None:
                done[st] = "no placement.json"; continue
            if "tilt_median_deg" in P:
                done[st] = "in placement.json"; continue
            PLM.orientation_only(ctx.cfg, ctx.bdir, inp, P["settings"]["tiles"], os.path.join(ctx.core_dir, f))
        elif st == "placement":
            PLM.run(ctx.cfg, ctx.bdir, inp, os.path.join(ctx.core_dir, f), os.path.join(ctx.core_dir, "placement_tiles.tsv"),
                    True, nproc=nproc)
        elif st == "crowding":
            from . import crowding
            crowding.run(ctx.cfg, ctx.bdir, inp, os.path.join(ctx.core_dir, f), os.path.join(ctx.core_dir, "crowding.tsv"), True)
        elif st == "densities":
            from . import densities
            densities.run(ctx.cfg, ctx.bdir, inp, os.path.join(ctx.core_dir, f), True, nproc=nproc)
        else:
            raise ValueError(f"unknown core step {st} (suite core steps: {CORE_STEPS})")
        ctx._json.clear()
        done[st] = "run"
    return done


# ------------------------------------------------------------------ C1 localisation vs independent references
C1_SETS = {"opencell": "OpenCell (clean hold-out)", "hpa_consensus_all": "HPA consensus (all genes)",
           "hpa_jurkat": "HPA JURKAT per-line", "hpa_jurkat_holdout_half": "HPA JURKAT held-out half (clean)"}
C1_REF = {("opencell", "agree"): 0.831, ("hpa_consensus_all", "agree"): 0.832}      # plan s4 table (v7, intended)
C1_PLAN_YEAST = {("SWAT_native", "tophit"): 0.665, ("hyperLOPIT", "tophit"): 0.848}  # plan s4 table (X5, top-1)


def c1(ctx):
    rows = []
    H = ctx.vjson("loc_holdout.json")
    if H is not None:
        L = ctx.vjson("loc_holdout.json", legacy=True)
        for s in C1_SETS:
            if s not in H["sets"]:
                rows.append(na("C1", f"{s}", f"n/a (no data): set {s} not in loc_holdout.json"))
                continue
            for lvl in ("intended", "placed"):
                for m in ("agree", "agree_tvd", "tophit"):
                    for w in ("unweighted", "mass_weighted"):
                        v = H["sets"][s][lvl][m][w]
                        lg = _get(L, "sets", s, lvl, m, w)
                        rf = C1_REF.get((s, m)) if (lvl, w) == ("intended", "unweighted") else None
                        rows.append(row("C1", f"{s}.{lvl}.{m}.{w}", v, direction="up", tol="loc", reference=rf,
                                        reference_source="plan s4 (v7)" if rf else "", source=SRC_CORE + " (loc_compare)",
                                        file="loc_holdout.json", legacy=lg,
                                        legacy_source="legacy v7 run: loc_holdout.json" if lg is not None else "",
                                        note=f"n_genes {H['sets'][s]['n_genes']}"))
        return rows
    if ctx.adapter is not None and hasattr(ctx.adapter, "loc_inputs"):
        return c1_common(ctx)
    a = ctx.adapt("C1")
    if a is not None:
        return a
    return [na("C1", "localisation", _why_missing(ctx, "loc_compare", "loc_holdout.json"))]


# legacy yeast localization.py names of the common metrics (same formula; see locmetrics docstring)
C1_YEAST_LEGACY = {("agree", "unweighted"): "agree", ("agree", "mass_weighted"): "agree_mass_weighted",
                   ("tophit", "unweighted"): "top1"}


def c1_common(ctx):
    """C1 with the common metric (locmetrics) on the model shares and reference sets of the cell's adapter"""
    from . import locmetrics
    D = ctx.adapter.loc_inputs()
    rows = []
    pergene = []
    for s, ref in D["sets"].items():
        rep = s in D.get("report_only", ())
        for lvl in ("intended", "placed"):
            summ, per = locmetrics.score_set(D["model"][lvl], ref, D["comps"])
            pergene += [dict(set=s, level=lvl, **r) for r in per]
            lg_row = D["legacy"].get((s, lvl)) or {}
            for m in ("agree", "agree_tvd", "tophit"):
                for w in ("unweighted", "mass_weighted"):
                    lk = C1_YEAST_LEGACY.get((m, w))
                    lg = lg_row.get(lk) if lk else None
                    rf = C1_PLAN_YEAST.get((s, m)) if (lvl, w) == ("intended", "unweighted") else None
                    note = f"n_genes {summ['n_genes']}" + (f" (legacy n {lg_row.get('n')})" if lg_row else "")
                    if rep:
                        note += "; CYCLoPs is a model INPUT (fit, not a test): report only"
                    rows.append(row("C1", f"{s}.{lvl}.{m}.{w}", summ[m][w], direction="report" if rep else "up", tol="loc",
                                    reference=rf, reference_source="plan s4 (X5, top-1)" if rf else "",
                                    source=SRC_COMMON + " + adapter (yeast model shares, reference sets)",
                                    file="adapter/loc_inputs", legacy=lg,
                                    legacy_source=f"{D['legacy_file']} ({lk}, share mode, test split)" if lg is not None else "",
                                    legacy_note=("" if lg is not None else ""), note=note))
    if pergene:
        common.write_df(os.path.join(ctx.out, "c1_pergene.tsv"), pd.DataFrame(pergene), overwrite=True, index=False)
    if hasattr(ctx.adapter, "c1_extras"):
        rows += ctx.adapter.c1_extras()
    return rows


# ------------------------------------------------------------------ C2 localisation across versions
def c2(ctx):
    a = ctx.adapt("C2")
    if a is not None:
        return a
    p = ctx.sc.get("across_versions")
    if not p:
        return [na("C2", "localisation across versions",
                   "n/a (no data): no archived-version localisation table for this cell (validate.suite.across_versions)")]
    if not os.path.exists(p):
        return [na("C2", "localisation across versions", f"n/a (no data): {p} missing")]
    T = pd.read_csv(p, sep="\t")
    H = ctx.vjson("loc_holdout.json")
    rows = []
    want = [("hpa_consensus", "all", "hpa_consensus_all"), ("opencell", "all", "opencell"),
            ("hpa_jurkat", "all_excl_CD3", "hpa_jurkat"), ("hpa_jurkat", "holdout_half", "hpa_jurkat_holdout_half")]
    for ref, st, c1set in want:
        for _, r in T[(T.ref == ref) & (T.set == st)].iterrows():
            for m, c1m, w in (("agree", "agree", "unweighted"), ("agree_mass", "agree", "mass_weighted"), ("tophit", "tophit", "unweighted")):
                lg = None; ls = ""
                if r.version == ctx.sc.get("version") and H is not None and c1set in H["sets"]:
                    lg = H["sets"][c1set]["intended"][c1m][w]; ls = f"C1 {c1set}.intended.{c1m}.{w} of this build (cross-check)"
                ln = ("mass weights differ: imaging_across_versions.py weights by inventory_v0 mass_da x copies, loc_compare by the "
                      "model's per-gene mass (differences ~1e-4 expected)") if (lg is not None and m == "agree_mass") else ""
                rows.append(row("C2", f"{ref}.{st}.{r.version}.{m}", float(r[m]), direction="report", source=SRC_JURKAT_LEGACY,
                                file=p, legacy=lg, legacy_source=ls, legacy_note=ln, note=f"n {int(r.n)}; imaging_across_versions.py annot step"))
    return rows


# ------------------------------------------------------------------ rank statistics shared by C3, C6
def auroc(score, pos):
    """AUROC of `score` for positives `pos` (ties averaged; Mann-Whitney). origin: yeast validation/proximity_test.py:auroc
    @ e3b7448 (equal to sklearn roc_auc_score used by image_nc.py)"""
    from scipy.stats import rankdata
    pos = np.asarray(pos, bool)
    n1, n0 = int(pos.sum()), int((~pos).sum())
    if n1 == 0 or n0 == 0:
        return None
    r = rankdata(score)
    return float((r[pos].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def spearman(a, b):
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).correlation) if len(a) > 2 else None


# ------------------------------------------------------------------ C3 image N/C
C3_NUC, C3_CYT, C3_MIN_CLASS = 1.0, -1.0, 5


def image_nc_metrics(measured, model):
    """one definition (plan s4 C3; origin: yeast validation/image_nc.py:main @ e3b7448): per gene log2 N/C of the model
    and of the image; Spearman over the genes; AUROC of the image ratio for model-nuclear (model log2 N/C > 1) against
    model-cytosolic (< -1) genes; both classes need >= 5 genes"""
    measured = np.asarray(measured, float); model = np.asarray(model, float)
    nucl, cyto = model > C3_NUC, model < C3_CYT
    au = None
    if nucl.sum() >= C3_MIN_CLASS and cyto.sum() >= C3_MIN_CLASS:
        au = auroc(np.r_[measured[nucl], measured[cyto]], np.r_[np.ones(int(nucl.sum())), np.zeros(int(cyto.sum()))])
    return dict(n=int(len(measured)), spearman=spearman(measured, model), auroc_nuclear=au, model_nuclear=int(nucl.sum()),
                model_cytosolic=int(cyto.sum()))


def c3(ctx):
    T = None
    if ctx.adapter is not None and hasattr(ctx.adapter, "image_nc_table"):
        T = ctx.adapter.image_nc_table()
    elif ctx.sc.get("image_nc"):
        T = _image_nc_from_config(ctx)
    if T is None:
        return [na("C3", "image_nc.spearman", "n/a (no data): no per-gene image N/C table for this cell (validate.suite.image_nc)")]
    if T.get("missing"):            # the cell's image test did not run (e.g. yeast: raw image archive and cache absent)
        print(f"suite C3: {T['missing']}")
        return [na("C3", "image_nc.spearman", f"n/a (not run): {T['missing']}")]
    rows = []
    for sp in T["splits"]:
        d = T["df"] if sp == "all" else T["df"][T["df"].split == sp]
        r = image_nc_metrics(d.measured_log2.values, d.model_log2.values)
        lg = (T.get("legacy") or {}).get(sp, {})
        # a table whose model images come from another build (config image_nc.model_of_build) is reported, never scored
        # as this build's value (review B17: the Jurkat rows were the v7 pack's)
        other = T.get("of_build")
        main = sp == T["test"] and not other
        for m, ref in (("spearman", T.get("ref_spearman")), ("auroc_nuclear", T.get("ref_auroc"))):
            nm = f"image_nc.{sp}.{m}"
            if r[m] is None:
                rows.append(na("C3", nm, f"n/a (no data): {r['model_nuclear']} model-nuclear and {r['model_cytosolic']} "
                               f"model-cytosolic genes (needs >= {C3_MIN_CLASS} each); {T.get('why_auroc', '')}",
                               source=SRC_SUITE_COMMON + " + " + T["source"], file=T["file"]))
                continue
            rows.append(row("C3", nm, r[m], direction="up" if main else "report", tol="loc",
                            reference=ref if main else None, reference_source=T.get("ref_source", "") if (main and ref) else "",
                            source=SRC_SUITE_COMMON + " + " + T["source"], file=T["file"], legacy=lg.get(m),
                            legacy_source=T.get("legacy_source", "") if lg.get(m) is not None else "",
                            legacy_note=T.get("legacy_note", {}).get(m, ""),
                            note=(f"model images of build {other}, not of this build; " if other else "")
                                 + f"n {r['n']} genes ({sp}); model-nuclear {r['model_nuclear']}, model-cytosolic {r['model_cytosolic']}; "
                                 f"{T.get('note', '')}"))
    return rows


def _image_nc_from_config(ctx):
    """legacy per-gene image table declared in config (Jurkat: Tier-2 HPA-Jurkat table of the v7 pack)"""
    c = ctx.sc["image_nc"]
    p = c["table"]
    if not os.path.exists(p):
        return None
    D = pd.read_csv(p, sep="\t")
    df = pd.DataFrame(dict(gene=D[c.get("gene", "gene")], measured_log2=D[c["measured"]], model_log2=D[c["model"]],
                           split=D[c["split"]] if c.get("split") in D.columns else "all")).dropna(subset=["measured_log2", "model_log2"])
    leg = {}
    if c.get("legacy_json") and os.path.exists(c["legacy_json"]):
        J = json.load(open(c["legacy_json"]))
        for sp, path in (c.get("legacy_keys") or {}).items():
            v = _get(J, *path.split("/"))
            if v is not None:
                leg.setdefault(sp, {})["spearman"] = v
    return dict(df=df, splits=c.get("splits", ["all"]), test=c.get("test", "all"), source=SRC_JURKAT_LEGACY + f" ({os.path.basename(p)})",
                file=p, legacy=leg, legacy_source=c.get("legacy_json", ""),
                legacy_note={"spearman": c.get("legacy_note", "")}, note=c.get("note", ""), why_auroc=c.get("why_auroc", ""),
                ref_spearman=c.get("ref_spearman"), ref_source=c.get("ref_source", ""),
                of_build=c.get("model_of_build"))


# ------------------------------------------------------------------ C6 copy numbers vs an independent dataset
def copies_vs(model, ref):
    """one definition (plan s4 C6; origin: yeast validation/copy_numbers.py:main @ e3b7448): proteins with model > 0 and
    reference > 0; Spearman of copies, median model/reference ratio, share within 2-fold"""
    m = pd.DataFrame(dict(model=model, ref=ref)).dropna()
    m = m[(m.model > 0) & (m.ref > 0)]
    L = np.log10(m)
    return dict(n=int(len(m)), spearman=spearman(L.model, L.ref), median_ratio=float(10 ** np.median(L.model - L.ref)) if len(m) else None,
                within_2x=float((abs(L.model - L.ref) < np.log10(2)).mean()) if len(m) else None)


def c6(ctx):
    R = ctx.adapter.copies_reference() if ctx.adapter is not None and hasattr(ctx.adapter, "copies_reference") else None
    p = os.path.join(ctx.bdir, "inventory", "copies", "copies.tsv")
    if R is None:
        why = ctx.sc.get("copies_reference_na") or "no independent copy-number dataset declared for this cell"
        return [na("C6", "copies_vs_independent.spearman", f"n/a (no data): {why}")]
    if not os.path.exists(p):
        return [na("C6", "copies_vs_independent.spearman", f"n/a (no data): {p} missing (core copies step not run)")]
    C = pd.read_csv(p, sep="\t")
    key = R["key"]
    model = pd.Series(C["copies"].values, index=C[key].astype(str).str.upper()).groupby(level=0).sum()
    rows = []
    for name, ref in R["sets"].items():
        ref = ref.groupby(level=0).first()
        r = copies_vs(model.reindex(ref.index), ref)
        lg = R.get("legacy", {}).get(name, {})
        indep = name in R.get("independent", ())
        for m, tol, d in (("spearman", "loc", "up" if indep else "report"), ("median_ratio", "report", "report"), ("within_2x", "report", "report")):
            rows.append(row("C6", f"{name}.{m}", r[m], direction=d, tol=tol, reference=R.get("plan", {}).get((name, m)),
                            reference_source="plan s4 (X5)" if R.get("plan", {}).get((name, m)) else "",
                            source=SRC_SUITE_COMMON + " (core inventory/copies/copies.tsv) + adapter (reference table)",
                            file="inventory/copies/copies.tsv + " + R["file"], legacy=lg.get(m), legacy_source=R.get("legacy_source", "") if lg.get(m) is not None else "",
                            legacy_note=R.get("legacy_note", "") if lg.get(m) is not None else "",
                            note=f"n {r['n']} proteins (model > 0 and reference > 0)" + ("" if indep else "; not independent of the model input")))
    return rows


# ------------------------------------------------------------------ C5 copy-number anchors
def c5(ctx):
    p = os.path.join(ctx.bdir, "inventory", "copies", "summary.json")
    if os.path.exists(p):
        S = json.load(open(p))
        rows = []
        leg = ctx.adapter.legacy("C5") if ctx.adapter is not None and hasattr(ctx.adapter, "legacy") else {}
        for a, counts in (S.get("anchors") or {}).items():
            for lab, d in counts.items():
                if not d.get("primary"):
                    continue
                lg, ls = leg.get(a, (None, ""))
                rel = (leg.get("_notes") or {}).get(a, "")
                rows.append(row("C5", f"{a}.{lab}.ratio", d["ratio"], direction="target:1", tol="anchor", source=SRC_CORE + " (inventory copies anchor check)",
                                file=p, legacy=lg, legacy_source=ls,
                                note=f"median over {d.get('n_families')} families of copies / (stoichiometry x {d.get('count')} {d.get('structure')}); implied total {d.get('implied_total_pg', float('nan')):.3g} pg" + (f"; {rel}" if rel else "")))
                rows.append(row("C5", f"{a}.{lab}.implied_total_pg", d.get("implied_total_pg"), unit="pg", direction="report",
                                source=SRC_CORE + " (inventory copies anchor check)", file=p))
        if rows:
            return rows
    a = ctx.adapt("C5")
    if a is not None:
        return a
    return [na("C5", "anchors", "n/a (no data): no inventory/copies/summary.json with anchors (copy numbers not built by the core copies step)")]


# ------------------------------------------------------------------ C7 clashes
def c7(ctx):
    rows = []
    P = ctx.vjson("placement.json")
    if P is not None:
        L = ctx.vjson("placement.json", legacy=True)
        npk = P["n_packed"]; npre = P["n_prepass"]
        f = "placement.json"; src = ctx.vsource("placement.json") + " (placement)"
        ls = "legacy v7 run: placement.json"
        rows.append(row("C7", "a.clash_rate_per_1e6", P["clash_rate_per_1e6"], unit="per 1e6 packed", direction="down", tol="rate1e6",
                        reference=0.0, reference_source="plan s4", source=src, file=f, legacy=_get(L, "clash_rate_per_1e6"), legacy_source=ls,
                        note=f"sphere-tree overlap > {P['settings']['clash_nm']} nm, packed-packed + packed-pre-pass, {P['n_tiles_sampled']} seeded tiles; "
                             f"collision shapes: {os.path.basename(str(P['settings'].get('collision_shapes', 'library sphere trees (st_c, st_r)')))}"))
        for k, den, nm in (("packed_packed", npk, "packed"), ("packed_prepass", npk, "packed"), ("packed_mrna", npk, "packed")):
            v = P["clash_pairs"][k] / den * 1e6 if den else None
            lv = (_get(L, "clash_pairs", k) / _get(L, "n_packed") * 1e6) if L else None
            rows.append(row("C7", f"a.{k}_per_1e6", v, unit=f"pairs per 1e6 {nm}", direction="down", tol="rate1e6", source=src,
                            file=f, legacy=lv, legacy_source=ls if L else ""))
        for k in ("prepass_prepass", "prepass_mrna"):
            if not npre:
                rows.append(na("C7", f"a.{k}_per_1e6", "n/a (no data): no pre-pass instances in the sampled tiles", file=f))
                continue
            v = P["clash_pairs"][k] / npre * 1e6
            lv = (_get(L, "clash_pairs", k) / _get(L, "n_prepass") * 1e6) if L else None
            rows.append(row("C7", f"a.{k}_per_1e6", v, unit="pairs per 1e6 pre-pass", direction="report", source=src, file=f,
                            legacy=lv, legacy_source=ls if L else "",
                            note="pre-pass structures are built, not packed: their self-overlaps (polysome chains) are reported, not gated"))
        rows.append(row("C7", "a.max_clash_depth_nm", P["max_clash_depth_nm"], unit="nm", direction="report", source=src, file=f,
                        legacy=_get(L, "max_clash_depth_nm"), legacy_source=ls if L else ""))
    else:
        rows.append(na("C7", "a.clash_rate_per_1e6", _why_missing(ctx, "placement", "placement.json") +
                       ("; phase 2: run the core placement step on this cell's tiles" if ctx.adapter else "")))
    a = ctx.adapt("C7")
    if a is not None:
        rows += a
    else:
        rows.append(na("C7", "b.bead_clash_rate", "n/a (no data): bead-level clash metric not computed for this cell (phase 2: common bead-overlap code)"))
    return rows


# ------------------------------------------------------------------ C8 containment, C9 orientation
def c8(ctx):
    rows = []
    P = ctx.vjson("placement.json")
    if P is not None:
        L = ctx.vjson("placement.json", legacy=True)
        rows.append(row("C8", "a.containment_rate_per_1e6", P["containment_rate_per_1e6"], unit="per 1e6 packed", direction="down",
                        tol="rate1e6", reference=0.0, reference_source="plan s4", source=ctx.vsource("placement.json") + " (placement)", file="placement.json",
                        legacy=_get(L, "containment_rate_per_1e6"), legacy_source="legacy v7 run: placement.json" if L else "",
                        note="instance centres whose geometry label is not allowed for their region (volume), off-membrane > tol (membrane), outside the shell (shell)"))
        regs = {r["name"]: r["kind"] for r in __import__("wholecell.validate.regions", fromlist=["x"]).load_regions(ctx.cfg)} \
            if (ctx.cfg.get("pack") or {}).get("regions") else {}
        if regs:
            by = {}
            for name, d in (P.get("containment_by_region") or {}).items():
                k = regs.get(name, "other")
                b = by.setdefault(k, [0, 0]); b[0] += d["bad"]; b[1] += d["n"]
            for k, (bad, n) in sorted(by.items()):
                rows.append(row("C8", f"a.containment_rate_per_1e6.{k}", bad / n * 1e6 if n else None, unit=f"per 1e6 packed {k}-region instances",
                                direction="report", source=ctx.vsource("placement.json") + " (placement)", file="placement.json",
                                note=f"{bad} of {n} instances of {k} regions"))
    else:
        rows.append(na("C8", "a.containment_rate_per_1e6", _why_missing(ctx, "placement", "placement.json") +
                       ("; phase 2: run the core placement step on this cell's tiles" if ctx.adapter else "")))
    a = ctx.adapt("C8")
    if a is not None:
        rows += a
    return rows


def _legacy(ctx, check):
    return ctx.adapter.legacy(check) if ctx.adapter is not None and hasattr(ctx.adapter, "legacy") else {}


def c9(ctx):
    rows = []
    P = ctx.vjson("placement.json")
    leg = _legacy(ctx, "C9")
    if P is not None:
        L = ctx.vjson("placement.json", legacy=True)
        rows.append(row("C9", "sign_errors_per_1e6", P["orientation_error_rate_per_1e6"], unit="per 1e6 membrane instances",
                        direction="down", tol="rate1e6", reference=0.36 if ctx.cfg.get("cell_id") == "jurkat" else None,
                        reference_source="plan s4 (v7)" if ctx.cfg.get("cell_id") == "jurkat" else "", source=ctx.vsource("placement.json") + " (placement)",
                        file="placement.json", legacy=_get(L, "orientation_error_rate_per_1e6"),
                        legacy_source="legacy v7 run: placement.json" if L else "",
                        note=f"{P['n_membrane_instances']} membrane instances in the sampled tiles"))
        O = P if "tilt_median_deg" in P else ctx.vjson("orientation.json")
        if O is None or O.get("tilt_median_deg") is None:
            rows.append(na("C9", "median_tilt_deg", "n/a (no data): no membrane instances with tilt in the placement sample"))
        else:
            src = ctx.vsource("placement.json" if O is P else "orientation.json") + " (placement: orientation)"
            if O is not P and O.get("orientation_error_rate_per_1e6") != P["orientation_error_rate_per_1e6"]:
                print(f"suite C9: orientation re-run sign rate {O.get('orientation_error_rate_per_1e6')} != placement.json "
                      f"{P['orientation_error_rate_per_1e6']}")
            for m, unit, d, tol in (("median_tilt_deg", "degree", "down", "deg"), ("frac_tilt_lt_1deg", "", "report", "report"),
                                    ("p99_tilt_deg", "degree", "report", "report")):
                k = {"median_tilt_deg": "tilt_median_deg", "frac_tilt_lt_1deg": "tilt_frac_lt_1deg", "p99_tilt_deg": "tilt_p99_deg"}[m]
                lg, ls = leg.get(m, (None, ""))
                rows.append(row("C9", m, O.get(k), unit=unit, direction=d, tol=tol,
                                reference=0.0 if m == "median_tilt_deg" else None, reference_source="upright" if m == "median_tilt_deg" else "",
                                source=src, file="placement.json" if O is P else "orientation.json", legacy=lg, legacy_source=ls,
                                legacy_note=(leg.get("_note", "") if lg is not None else ""),
                                note=f"angle between the library +z axis and the membrane normal (SDF gradient) at the centre, "
                                     f"{O.get('tilt_n')} membrane instances of the placement sample"
                                     + (f"; signs re-checked: {O.get('orientation_error_rate_per_1e6')} per 1e6" if O is not P else "")))
    a = ctx.adapt("C9")
    if a:
        rows += a
    elif P is None:
        rows.append(na("C9", "sign_errors_per_1e6", _why_missing(ctx, "placement", "placement.json")))
    return rows


# ------------------------------------------------------------------ C10 placement integrity + seam index
def seam_counts(pos_local, tile_nm, band=SEAM_BAND_NM, interior=SEAM_INTERIOR_NM):
    """(n_band, n_interior) instance centres of ONE tile; pos_local = centres relative to the tile's low corner (nm)"""
    p = np.asarray(pos_local, np.float64)
    if len(p) == 0:
        return 0, 0
    d = np.minimum(p, tile_nm - p).min(1)
    return int((d < band).sum()), int((d >= interior).sum())


def seam_volumes(tile_nm, band=SEAM_BAND_NM, interior=SEAM_INTERIOR_NM):
    """(band volume, interior volume) of one cubic tile, nm3"""
    return tile_nm ** 3 - (tile_nm - 2 * band) ** 3, (tile_nm - 2 * interior) ** 3


def seam_index(n_band, n_interior, tile_nm, n_tiles=1, band=SEAM_BAND_NM, interior=SEAM_INTERIOR_NM):
    """density in the face band / density in the interior (summed over n_tiles tiles of the same size)"""
    vb, vi = seam_volumes(tile_nm, band, interior)
    if n_interior <= 0:
        return None
    return (n_band / (vb * n_tiles)) / (n_interior / (vi * n_tiles))


def seam(ctx, region=None, min_frac=SEAM_MIN_FRAC, n_tiles=SEAM_N_TILES, seed=None):
    """seam index of a build (definition in the module docstring) -> dict"""
    sc = ctx.sc.get("seam") or {}
    region = region or sc.get("region") or "cytosol"
    min_frac = float(sc.get("min_frac", min_frac)); n_tiles = int(sc.get("n_tiles", n_tiles))
    seed = int(sc.get("seed", 0) if seed is None else seed)
    meta = json.load(open(ctx.pack("tiles_meta.json")))
    T = float(meta["TILE"]); lo0 = np.asarray(meta.get("lo0", [0, 0, 0]), float)
    tp = pd.read_parquet(ctx.pack("tiles.parquet"))
    if region not in tp.columns:
        return dict(error=f"region '{region}' not a tile volume column of tiles.parquet")
    sel = tp[(tp[region] / T ** 3) >= min_frac]
    n_q = len(sel)
    if n_q == 0:
        return dict(error=f"no tile has >= {min_frac} of its volume in {region}", n_qualifying=0, tile_nm=T)
    rng = np.random.default_rng(seed)
    if len(sel) > n_tiles:
        sel = sel.iloc[np.sort(rng.choice(len(sel), n_tiles, replace=False))]
    nb = np.zeros(8); ni = np.zeros(8); nt = np.zeros(8)
    for _, r in sel.iterrows():
        f = ctx.pack(os.path.join("tiles", f"{int(r['tile'])}.npz"))
        if not os.path.exists(f):
            continue
        ijk = np.array([int(r.ti), int(r.tj), int(r.tk)])
        z = np.load(f)
        b, i = seam_counts(z["pos"].astype(np.float64) - (lo0 + ijk * T), T)
        ph = (ijk[0] & 1) + 2 * (ijk[1] & 1) + 4 * (ijk[2] & 1)
        nb[ph] += b; ni[ph] += i; nt[ph] += 1
    tot = seam_index(nb.sum(), ni.sum(), T, nt.sum())
    per = [seam_index(nb[k], ni[k], T, nt[k]) for k in range(8) if nt[k] > 0 and ni[k] > 0]
    return dict(index=tot, by_parity_min=min(per) if per else None, by_parity_max=max(per) if per else None,
                n_tiles=int(nt.sum()), n_qualifying=int(n_q), n_band=int(nb.sum()), n_interior=int(ni.sum()), tile_nm=T,
                region=region, min_frac=min_frac, seed=seed)


def c10(ctx):
    rows = []
    src = SRC_NEW + " from tile job tables (common.tile_sums)"
    if os.path.exists(ctx.pack("tiles.parquet")) and os.path.exists(ctx.pack("libs.json")):
        inp = common.Inputs(ctx.cfg, ctx.bdir)
        cache = os.path.join(ctx.bdir, "validate", "cache", "tile_sums.npz")
        if not os.path.exists(cache):
            cache = os.path.join(ctx.out, "cache", "tile_sums.npz")
        S = common.tile_sums(inp, cache, nproc=min(3, common.nprocs(ctx.cfg, ctx.kw)))
        req = S["req"]; pla = S["pla"]
        legacy = ctx.adapter.legacy("C10") if ctx.adapter is not None and hasattr(ctx.adapter, "legacy") else {}
        P = ctx.vjson("placement.json")
        lg = ls = None
        if P is not None and P.get("species_share"):
            ss = P["species_share"]
            lg = sum(v["placed"] for v in ss.values()) / sum(v["requested"] for v in ss.values())
            ls = "placement.json species_share (core placement step), summed over regions"
        elif "placed_copies_fraction" in legacy:
            lg, ls = legacy["placed_copies_fraction"]
        rows.append(row("C10", "placed_copies_fraction", float(pla.sum() / req.sum()), direction="up", tol="placed", source=src,
                        file="pack/tiles/*.npz", legacy=lg, legacy_source=ls or "",
                        note="packed instances placed / requested by the packer (pre-pass structures not included)"))
        extra = ctx.adapter.core_inputs() if ctx.adapter is not None and hasattr(ctx.adapter, "core_inputs") else {}
        m, msrc = common.lib_masses(common.Inputs(ctx.cfg, ctx.bdir, **extra), ctx.bdir, nproc=min(3, common.nprocs(ctx.cfg, ctx.kw)))
        if len(m) == req.shape[0]:
            sm = json.load(open(ctx.pack("summary.json"))) if os.path.exists(ctx.pack("summary.json")) else {}
            rows.append(row("C10", "placed_mass_fraction", float((pla.sum(1) * m).sum() / (req.sum(1) * m).sum()), direction="up",
                            tol="placed", reference=0.964 if ctx.cfg.get("cell_id") == "jurkat" else None,
                            reference_source="plan s4 (v7 pack summary)" if ctx.cfg.get("cell_id") == "jurkat" else "", source=src,
                            file=f"pack/tiles/*.npz + {msrc}", legacy=sm.get("mass_frac"),
                            legacy_source="pack/summary.json mass_frac (placed_pg / expected_pg; expected before packing)" if "mass_frac" in sm else "",
                            note=f"sum over (library, region) of placed x mass / requested x mass; mass per library entry = {msrc}"))
        else:
            rows.append(na("C10", "placed_mass_fraction", "n/a (no data): no mass per library entry (libs.json mass_da or library beads)"))
        tot_r = req.sum(1); tot_p = pla.sum(1); k = tot_r >= MIN_REQ
        nb50 = int((tot_p[k] / tot_r[k] < 0.5).sum())
        rel = ""
        if P is not None and P.get("species_share"):
            rel = (f"related legacy number (not the same definition): placement.json per-(species, region) pairs below 0.5, "
                   f"summed over regions = {sum(v['n_below_0_5'] for v in P['species_share'].values())}")
        elif "species_below_50pct" in legacy:
            rel = f"related legacy number (not the same definition): {legacy['species_below_50pct'][1]} = {legacy['species_below_50pct'][0]}"
        rows.append(row("C10", "species_below_50pct_placed", nb50, direction="report", source=src, file="pack/tiles/*.npz",
                        note=f"species (summed over regions) with >= {MIN_REQ} requested copies, of {int(k.sum())}; {rel}"))
        s = seam(ctx)
        if s.get("error"):
            rows.append(na("C10", "seam_index", f"n/a (no data): {s['error']}"))
        else:
            ln = legacy.get("seam_note", "")
            note = (f"{s['n_tiles']} of {s['n_qualifying']} tiles >= {s['min_frac']} {s['region']}, tile {s['tile_nm']:.0f} nm, "
                    f"{s['n_band']} centres in the {SEAM_BAND_NM:.0f} nm band, {s['n_interior']} in the interior (>= {SEAM_INTERIOR_NM:.0f} nm)")
            rows.append(row("C10", "seam_index", s["index"], direction="target:1", tol="seam", reference=1.0,
                            reference_source="no seam", source=SRC_NEW, file="pack/tiles.parquet + pack/tiles/*.npz",
                            legacy=None, legacy_source="none (new metric) " + ln, note=note))
            rows.append(row("C10", "seam_index_parity_min", s["by_parity_min"], direction="report", source=SRC_NEW,
                            note="min over the 8 tile parity classes"))
            rows.append(row("C10", "seam_index_parity_max", s["by_parity_max"], direction="report", source=SRC_NEW,
                            note="max over the 8 tile parity classes"))
    else:
        rows.append(na("C10", "placed_copies_fraction", "n/a (no data): no pack/tiles.parquet + pack/libs.json in the build"))
    a = ctx.adapt("C10")
    if a is not None:
        rows += a
    return rows


# ------------------------------------------------------------------ C11 densities and crowding
C11_REF = {"cytosol": (162.8, "plan s4 (v7)"), "nucleoplasm": (135.1, "plan s4 (v7)")}


def c11(ctx):
    rows = []
    D = ctx.vjson("densities.json")
    if D is not None:
        L = ctx.vjson("densities.json", legacy=True)
        lrows = {r["region"]: r for r in (L or {}).get("rows", [])}
        jur = ctx.cfg.get("cell_id") == "jurkat"
        leg = _legacy(ctx, "C11")              # cell legacy numbers of another definition (yeast stage 4), side by side
        refs = leg.get("_refs", {})
        for r in D["rows"]:
            rf = C11_REF.get(r["region"]) if jur else refs.get(r["region"])
            lg = lrows.get(r["region"], {}).get("mg_per_ml")
            ls = "legacy v7 run: densities.json" if lg is not None else ""
            ln = ""
            if lg is None and f"density.{r['region']}" in leg:
                lg, ls = leg[f"density.{r['region']}"]
                ln = leg.get("_density_note", "")
            rows.append(row("C11", f"density.{r['region']}", r["mg_per_ml"], unit="mg/ml", direction="report" if r["volume_um3"] < 1 else "target:range",
                            tol="density", reference=rf[0] if rf else None, reference_source=rf[1] if rf else "",
                            source=ctx.vsource("densities.json") + " (densities)", file="densities.json", legacy=lg,
                            legacy_source=ls, legacy_note=ln,
                            note=f"packed + pre-pass + RNA; volume {r['volume_um3']} um3; packed mass: {D.get('packed_mass_source', 'libs.json mass_da')}"))
        lg = (L or {}).get("nuc_cyto_ratio"); ls = "legacy v7 run: densities.json" if lg is not None else ""; ln = ""
        if lg is None and "nuc_cyto_ratio" in leg:
            (lg, ls), ln = leg["nuc_cyto_ratio"], leg.get("_ratio_note", "")
        rows.append(row("C11", "nuc_cyto_ratio", D.get("nuc_cyto_ratio"), direction="report", tol="ratio", source=ctx.vsource("densities.json") + " (densities)",
                        file="densities.json", legacy=lg, legacy_source=ls, legacy_note=ln,
                        note="nucleus group density / cytosol density"))
    C = ctx.vjson("crowding.json")
    if C is not None:
        L = ctx.vjson("crowding.json", legacy=True)
        lph = {r["region"]: r["phi"] for r in (L or {}).get("summary", [])}
        for r in C["summary"]:
            rows.append(row("C11", f"phi.{r['region']}", r["phi"], direction="report", tol="phi", source=ctx.vsource("crowding.json") + " (crowding)",
                            file="crowding.json", legacy=lph.get(r["region"]),
                            legacy_source="legacy v7 run: crowding.json" if r["region"] in lph else "",
                            note="excluded-volume fraction (bead model)"))
    if D is None and C is None:
        a = ctx.adapt("C11")
        if a is not None:
            return a
        return [na("C11", "density", _why_missing(ctx, "densities", "densities.json"))]
    if C is None:
        rows.append(na("C11", "phi", _why_missing(ctx, "crowding", "crowding.json")))
    # direction of a density: within the literature range; until ranges are declared (phase 2) the no-regression rule
    # treats it as two-sided (|change| > tolerance = regress)
    for r in rows:
        if r["direction"] == "target:range":
            r["direction"] = "stable"
    return rows


# ------------------------------------------------------------------ C12 ribosomes
def c12(ctx):
    R = ctx.vjson("ribosome_nn.json")
    if R is None:
        a = ctx.adapt("C12")
        if a is not None:
            return a
        return [na("C12", "ribosome_nn_median", _why_missing(ctx, "ribosome_nn", "ribosome_nn.json"))]
    L = ctx.vjson("ribosome_nn.json", legacy=True)
    by = {r["set"]: r for r in R["rows"]}; lby = {r["set"]: r for r in (L or {}).get("rows", [])}
    em = next((r for s, r in by.items() if s.startswith("EM")), None)
    rows = []
    ls = "legacy v7 run: ribosome_nn.json" if L else ""
    for key, name, ref in (("model pre-pass + packed RIBO80S", "nn_median_model", None),
                           ("model pre-pass + packed RIBO80S, thinned to EM count", "nn_median_model_thinned_to_EM", em)):
        r = by.get(key)
        if r is None:
            rows.append(na("C12", name, f"n/a (no data): row '{key}' not in ribosome_nn.json"))
            continue
        rows.append(row("C12", name, r["nn_median"], unit="nm", direction=f"target:{ref['nn_median']}" if ref else "report", tol="nn",
                        reference=ref["nn_median"] if ref else None, reference_source=ref["set"] if ref else "",
                        source=SRC_CORE + " (ribosome_nn)", file="ribosome_nn.json", legacy=_get(lby.get(key), "nn_median"),
                        legacy_source=ls, note=f"n {r['n']}"))
    D = ctx.vjson("densities.json")
    r = by.get("model pre-pass + packed RIBO80S")
    vol = next((x["volume_um3"] for x in (D or {}).get("rows", []) if x["region"] == "cytosol"), None)
    if r is not None and vol:
        c = r["n"] / (vol * 1e-15 * AVOGADRO) * 1e6
        rows.append(row("C12", "concentration_uM", c, unit="uM", direction="report", tol="rel10", source=SRC_NEW,
                        file="ribosome_nn.json + densities.json", note=f"{r['n']} ribosomes / cytosol {vol} um3 (all ribosomes, nucleus excluded by construction)"))
    return rows


# ------------------------------------------------------------------ C13 chromatin
def c13(ctx):
    F = ctx.vjson("chromatin_frozen.json")
    if F is None:
        a = ctx.adapt("C13")
        if a is not None:
            return a
        return [na("C13", "chromatin", _why_missing(ctx, "chromatin", "chromatin_frozen.json"))]
    L = ctx.vjson("chromatin_frozen.json", legacy=True)
    rows = [row("C13", "prepass_md5_ok", 1 if F["prepass"]["ok"] else 0, direction="up", tol="count0",
                source=SRC_CORE + " (chromatin, frozen)", file="chromatin_frozen.json", note=F["prepass"]["md5_actual"])]
    for k, h in F["headline"].items():
        def num(x):
            try:
                return float(x)
            except (TypeError, ValueError):
                return None
        lg = num(_get(L, "headline", k, "v8"))
        rows.append(row("C13", k, num(h["v8"]), direction="report", reference=num(h.get("target")),
                        reference_source=f"criterion {h['criterion']}; pass={h['pass_']}", source=SRC_CORE + " (chromatin, frozen v8 metrics)",
                        file="chromatin_frozen.json", legacy=lg, legacy_source="legacy v7 run: chromatin_frozen.json" if lg is not None else "",
                        note=h.get("use", "")))
    return rows


# ------------------------------------------------------------------ C16 build integrity
def c16(ctx):
    rows = []
    sp = os.path.join(ctx.bdir, "BUILD_STAMP.json")
    rows.append(row("C16", "build_stamp_present", 1 if os.path.exists(sp) else 0, direction="up", tol="count0",
                    source=SRC_CORE + " (BUILD_STAMP.json)", file="BUILD_STAMP.json"))
    mp = ctx.vfile("MANIFEST.tsv")
    if os.path.exists(mp):
        import hashlib
        stamp = None
        if os.path.exists(sp):
            stamp = "stamp:" + hashlib.sha256(open(sp, "rb").read()).hexdigest()[:16]
        M = pd.read_csv(mp, sep="\t", comment="#")
        stale = []
        for _, r in M.iterrows():
            p = ctx.vfile(r["file"])
            if not os.path.exists(p) or common.md5(p) != r["md5"] or (stamp and not str(r["build_stamp"]).startswith(stamp)):
                stale.append(r["file"])
        listed = set(M["file"])
        rows.append(row("C16", "manifest_stale_files", len(stale), direction="down", tol="count0", reference=0, reference_source="must",
                        source=SRC_CORE + " (validate MANIFEST.tsv)", file="MANIFEST.tsv",
                        note=f"{len(M)} files: md5 and build stamp re-checked; stale: {stale[:5]}"))
        extra = []
        for root, dirs, files in os.walk(ctx.vdir):
            dirs[:] = [d for d in dirs if d not in ("cache", "_smk")]
            for f in files:
                rel = os.path.relpath(os.path.join(root, f), ctx.vdir)
                if rel not in listed and rel != "MANIFEST.tsv" and ".tmp" not in rel:
                    extra.append(rel)
        rows.append(row("C16", "files_not_in_manifest", len(extra), direction="report", source=SRC_CORE + " (validate MANIFEST.tsv)",
                        file="MANIFEST.tsv", note=f"written after the manifest step (compare/suite outputs): {sorted(extra)[:6]}"))
    a = ctx.adapt("C16")
    if a is not None:
        rows += a
    elif not os.path.exists(mp):
        rows.append(na("C16", "manifest_stale_files", f"n/a (no data): no {mp}"))
    rows += _c16_depends_on_geometry(ctx)
    return rows


def _c16_depends_on_geometry(ctx):
    """frozen pre-placed modules whose depends_on_geometry != this build's geometry package (B6; contract in
    docs/reference/cell_config.md 3.4). Built modules depend on the build's package by construction."""
    from ..preplaced.stage import check_depends_on_geometry
    try:
        deps, bad = check_depends_on_geometry(ctx.cfg, ctx.bdir)
    except Exception as e:                          # noqa: BLE001  (no package / config: reported, not hidden)
        return [na("C16", "depends_on_geometry_mismatch", f"n/a (no data): {type(e).__name__}: {e}")]
    nf = sum(d["source"] == "frozen" for d in deps)
    serial = bool((ctx.cfg.get("preplaced") or {}).get("serial"))
    note = (f"{len(deps)} modules ({nf} frozen) vs package {deps[0]['package_id'] if deps else None}; "
            + (f"mismatches: {bad[:4]}" if bad else "all declared for this package")
            + ("; serial builders run in this build (depend on its package)" if serial else ""))
    return [row("C16", "depends_on_geometry_mismatch", len(bad), direction="down", tol="count0", reference=0,
                reference_source="must", source=SRC_CORE + " (config preplaced.modules[].depends_on_geometry)",
                file="preplaced/preplaced_manifest.json", note=note)]


# ------------------------------------------------------------------ phase-2 checks (n/a rows, never silent)
def _phase2(check, metric, what):
    def f(ctx):
        a = ctx.adapt(check)
        rows = [na(check, metric, f"n/a (phase 2): {what}")]
        if a:
            rows += a
        return rows
    return f


def c15(ctx):
    a = ctx.adapt("C15")
    if a:
        return a
    why = ctx.sc.get("proximity_na") or "no proximity / interaction reference wired for this cell"
    return [na("C15", "proximity_auroc", f"n/a (no data): {why}")]


def c17(ctx):
    rows = list(ctx.adapt("C17") or [])
    for e in ctx.sc.get("extras_na") or []:
        rows.append(na("C17", e["metric"], f"n/a ({e.get('kind', 'not ported')}): {e['reason']}"))
    if not rows:
        rows.append(na("C17", "cell_specific_extras", "n/a (no data): no cell-specific extras declared for this cell"))
    return rows


CHECKS = {
    "C1": c1,
    "C2": c2,
    "C3": c3,
    "C4": _phase2("C4", "image_compartment_shares", "1 - TVD per imaged gene (Jurkat v8 imaging) not wired yet"),
    "C5": c5,
    "C6": c6,
    "C7": c7,
    "C8": c8,
    "C9": c9,
    "C10": c10,
    "C11": c11,
    "C12": c12,
    "C13": c13,
    "C14": lambda ctx: [na("C14", "diffusion_benchmark", "n/a (BD v2): the diffusion benchmark waits for the common BD "
                           "engine (plan s6.1, one method for both cells); the yeast diffusion.py and Jurkat BD numbers are "
                           "different methods and are not reported as C14")],
    "C15": c15,
    "C16": c16,
    "C17": c17,
}


# ------------------------------------------------------------------ no-regression comparison
def compare(value, ref, direction, tol):
    """status of one metric against the reference build: ok | regress | improve | n/a (plan s4 rule)"""
    if value is None or ref is None:
        return "n/a"
    if direction in ("report", None) or tol is None:
        return "ok"
    d = value - ref
    if direction == "up":
        return "regress" if d < -tol else "improve" if d > tol else "ok"
    if direction == "down":
        return "regress" if d > tol else "improve" if d < -tol else "ok"
    if direction == "stable":
        return "regress" if abs(d) > tol else "ok"
    if direction.startswith("target:"):
        t = float(direction.split(":", 1)[1])
        dv, dr = abs(value - t), abs(ref - t)
        return "regress" if dv > dr + tol else "improve" if dv < dr - tol else "ok"
    raise ValueError(f"unknown direction {direction}")


def no_regression(rows, ref_rows, ref_name=""):
    """rows of this build vs the reference build's rows -> (table rows, rows with status filled)"""
    R = {(r["check"], r["metric"]): r for r in (ref_rows or [])}
    out = []
    for r in rows:
        rr = R.get((r["check"], r["metric"]))
        rv = rr["value"] if rr else None
        tol = TOL[r["tol_kind"]][0](rv) if rv is not None else None
        if r["value"] is None:
            st, why = "n/a", r["reason"] or f"n/a (no data): no value in {r['file'] or 'the inputs'}"
        elif ref_rows is None:
            st, why = "n/a", "n/a (no reference build declared: validate.suite.reference)"
        elif rr is None or rv is None:
            st, why = "n/a", f"n/a (metric has no value in the reference build {ref_name})"
        else:
            st, why = compare(r["value"], rv, r["direction"], tol), ""
        r["status"] = st
        r["tolerance"] = tol
        if why and not r["reason"]:
            r["reason"] = why
        out.append(dict(check=r["check"], metric=r["metric"], value=r["value"], reference_build_value=rv,
                        delta=(r["value"] - rv) if (r["value"] is not None and rv is not None) else None,
                        direction=r["direction"], tolerance=tol, tolerance_rule=r["tolerance_rule"], status=st,
                        reason=why, doc=r["doc"]))
    return out, rows


def _match(v, legacy):
    """identical, or equal at the precision the legacy file stores (legacy numbers rounded to n decimals)"""
    try:
        v, lg = float(v), float(legacy)
    except (TypeError, ValueError):
        return v == legacy
    txt = repr(legacy)
    dec = len(txt.split(".")[1]) if "." in txt and "e" not in txt.lower() else 12
    return abs(v - lg) <= max(0.5 * 10 ** -dec, 1e-9 * max(1.0, abs(lg))) + 1e-15


def legacy_table(rows):
    out = []
    for r in rows:
        if r["legacy"] is None or r["value"] is None:
            continue
        try:
            d = float(r["value"]) - float(r["legacy"])
        except (TypeError, ValueError):
            d = None
        out.append(dict(check=r["check"], metric=r["metric"], suite=r["value"], legacy=r["legacy"], delta=d,
                        match=_match(r["value"], r["legacy"]), source=r["source"], legacy_source=r["legacy_source"],
                        explanation=r.get("legacy_note", "")))
    return out


def _ref_rows(ctx, kw):
    """reference build's suite rows: ref_suite= > validate.suite.reference.suite; the reference build itself compares
    with its own rows"""
    rc = ctx.sc.get("reference") or {}
    rb = rc.get("build")
    if rb and os.path.realpath(rb.rstrip("/")) == os.path.realpath(ctx.bdir.rstrip("/")):
        return "self", "self (this build is the declared reference)"
    p = kw.get("ref_suite") or rc.get("suite")
    if not p:
        return None, "none declared"
    if not os.path.exists(p):
        print(f"suite: reference suite {p} not found (run the suite on the reference build first); no-regression n/a")
        return None, f"missing: {p}"
    return json.load(open(p))["metrics"], p


# ------------------------------------------------------------------ entry
def run(cfg, bdir, out=None, overwrite=True, checks=None, **kw):
    import time
    t0 = time.time()
    out = out or os.path.join(bdir, "validate")
    os.makedirs(out, exist_ok=True)
    if isinstance(checks, str):
        checks = [c.strip() for c in checks.split(",") if c.strip()]
    ctx = Ctx(cfg, bdir, out, **kw)
    common.CACHE_DIR = os.path.join(out, "cache") if os.path.realpath(out) != os.path.realpath(os.path.join(bdir, "validate")) else None
    cs = kw.get("core_steps")
    cs = ([c for c in cs.split(",") if c] if isinstance(cs, str) else cs) if cs is not None else \
        list(ctx.sc.get("core_steps") or []) + ["orientation"]
    ctx.core_run = run_core_steps(ctx, [c for c in CORE_STEPS if c in cs]) if cs else {}
    rows = []
    for cid, f in CHECKS.items():
        if checks and cid not in checks:
            continue
        try:
            got = f(ctx) or []
        except Exception as e:                      # a broken check must not hide the others; it is reported
            got = [na(cid, "error", f"n/a (error): {type(e).__name__}: {e}")]
            print(f"suite {cid}: ERROR {type(e).__name__}: {e}")
        rows += got
    for key, txt in (ctx.sc.get("notes") or {}).items():          # cell-specific explanations (config validate.suite.notes)
        c_, m_ = key.split("/", 1)
        for r in rows:
            if r["check"] == c_ and r["metric"] == m_:
                r["note"] = (r["note"] + "; " if r["note"] else "") + " ".join(str(txt).split())
    ref, ref_name = _ref_rows(ctx, kw)
    if ref == "self":
        ref = [dict(r) for r in rows]
    table, rows = no_regression(rows, ref, ref_name)
    for r in rows:
        if r["status"] == "regress":
            print(f"REGRESS {r['check']} {r['metric']}: value {r['value']}, reference build {next(t['reference_build_value'] for t in table if t['metric'] == r['metric'] and t['check'] == r['check'])}, "
                  f"tolerance {r['tolerance']} ({r['tolerance_rule']}), see {r['doc']}")
    by = {}
    for r in rows:
        b = by.setdefault(r["check"], dict(metrics=0, with_value=0, ok=0, regress=0, improve=0, **{"n/a": 0}))
        b["metrics"] += 1; b["with_value"] += r["value"] is not None; b[r["status"]] += 1
    stamp = os.path.join(bdir, "BUILD_STAMP.json")
    from ..build import git_info
    gi = git_info()
    res = dict(build=bdir, build_id=os.path.basename(bdir.rstrip("/")), cell=cfg.get("cell_id"),
               stamp_sha256=__import__("hashlib").sha256(open(stamp, "rb").read()).hexdigest() if os.path.exists(stamp) else None,
               suite_code=f"{gi['commit7']}{'-dirty' if gi['dirty'] else ''}", reference=ref_name,
               legacy_validate_dir=ctx.legacy_dir, adapter=ctx.sc.get("adapter"), core_steps_run=ctx.core_run,
               seconds=round(time.time() - t0, 1), summary=by, metrics=rows)
    common.write_json(os.path.join(out, "suite.json"), res, overwrite=True)
    common.write_df(os.path.join(out, "no_regression.tsv"), pd.DataFrame(table), overwrite=True, index=False)
    common.write_df(os.path.join(out, "legacy_check.tsv"), pd.DataFrame(legacy_table(rows)), overwrite=True, index=False)
    common.CACHE_DIR = None
    n = pd.Series([r["status"] for r in rows]).value_counts().to_dict()
    nv = sum(r["value"] is not None for r in rows)
    print(f"suite: {len(rows)} metrics in {len(by)} checks, {nv} with a value; status {n}; reference {ref_name}; {res['seconds']} s -> {out}")
    for cid, b in by.items():
        print(f"  {cid:4s} {b['with_value']}/{b['metrics']} with value; ok {b['ok']} regress {b['regress']} improve {b['improve']} n/a {b['n/a']}")
    return res
