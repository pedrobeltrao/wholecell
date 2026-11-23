"""SEC-MS in the inventory stage: the `secms` step (adapter.secms `secms_toolkit`) and the view that assembly-state
adapters use to read its tables (`SecView`) with the configurable assembly-state options (`assembly_options`).

Step `secms` (only when inventory.adapters.secms is set): runs toolkits/secms on inventory.params.secms_toolkit
(species human | yeast, preset, annotation, complexes, pair_candidates, datasets; every string naming a declared frozen
input is replaced by its path) and writes inventory/secms/{secms_protein, secms_complex, secms_complex_subunit,
secms_pairs, secms_peaks, secms_calibration, secms_calibrants, profiles}.tsv + diagnostics.json.
Optional `combine:` (2026-11-19, Step 3): a toolkits/secms/combine.py spec ({label, discordant_range, datasets: [{source_id,
weight, role primary | fill, replicates}], min_fractions}); the combined rows (source_id = label, condition = replicate =
'combined') are appended to every table, as `python -m toolkits.secms combine --spec`; the assembly rules then select
them with assembly_options.secms.source_id = label. Absent or null = no combination (the run tables only).

Assembly-state options (complex_table: inventory.params.complex_table.assembly_options; assembly_steps: the same
keys come from its preset, wholecell/inventory/assembly/steps.py legacy_options; defaults = the Jurkat v7 rules):
  secms: {source_id, condition, replicate}       which SEC-MS run feeds the rules (default: the only source /
                                                 condition; replicate 'mean' when present, else the first)
  complex_copies: stoichiometry_cap              copies of a complex = min(median, 2 x limiting) of subunit pools / stoich
                | measured_free_pool             copies = median over SEC-measured subunits of pool x in_complex_frac /
                                                 stoich (the measured free pool stays free), capped by what is left;
                                                 complexes without SEC coverage (or TM, see below) keep the cap rule
  measured_free_pool: {exclude_tm: true, min_subunits_detected: null, min_detected_frac: null (null = toolkit gate),
                       max_expected_kda: null (larger complexes keep the cap rule: pelleted by the pre-SEC clarification
                       spin or beyond the column), exclude_complex_ids: []}
  homo_oligomer_share: {source: secms | constant, k_by_state: [[3, 2], [5, 3], [null, 4]], default: 0.5, value: null,
                        apex_ratio: [0.4, 2.5], apex_check_min_k: 4, apex_check_min_frac: 0.2}
  afdb_homodimer_share: {source: secms | constant, k: 2, default: 0.3, value: null}
  heterodimer: {sec_r_min: 0.6, string_min: 700, bound_fraction: 0.5}
"""
from __future__ import annotations

import copy
import json
import os

import numpy as np
import pandas as pd

from ..registry import register
from . import schema

TABLES = ["secms_protein", "secms_complex", "secms_complex_subunit", "secms_pairs", "secms_peaks", "secms_calibration",
          "secms_calibrants", "profiles"]

ASSEMBLY_DEFAULTS = {
    "secms": {"source_id": None, "condition": None, "replicate": None},
    "complex_copies": "stoichiometry_cap",
    "measured_free_pool": {"exclude_tm": True, "min_subunits_detected": None, "min_detected_frac": None,
                           "max_expected_kda": None, "exclude_complex_ids": []},
    "homo_oligomer_share": {"source": "secms", "k_by_state": [[3, 2], [5, 3], [None, 4]], "default": 0.5, "value": None,
                            "apex_ratio": [0.4, 2.5], "apex_check_min_k": 4, "apex_check_min_frac": 0.2},
    "afdb_homodimer_share": {"source": "secms", "k": 2, "default": 0.3, "value": None},
    "heterodimer": {"sec_r_min": 0.6, "string_min": 700, "bound_fraction": 0.5},
}
COMPLEX_COPIES = ("stoichiometry_cap", "measured_free_pool")


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in (b or {}).items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else copy.deepcopy(v)
    return out


def assembly_options(params: dict) -> dict:
    o = _merge(ASSEMBLY_DEFAULTS, params.get("assembly_options") or {})
    if o["complex_copies"] not in COMPLEX_COPIES:
        raise ValueError(f"assembly_options.complex_copies must be one of {COMPLEX_COPIES}")
    return o


# ------------------------------------------------------------------------------------------------ step `secms`
def _frozen_paths(obj, ctx):
    names = set((ctx.cfg.get("frozen_inputs") or {}).keys())
    if isinstance(obj, str):
        head, _, rest = obj.partition("/")        # "<frozen name>/<path inside>" as the assembly inputs (2026-11-19)
        return ctx.frozen(obj) if obj in names else os.path.join(ctx.frozen(head), rest) if rest and head in names else obj
    if isinstance(obj, dict):                      # source_id stays a name (it labels the rows)
        return {k: (v if k == "source_id" else _frozen_paths(v, ctx)) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_frozen_paths(v, ctx) for v in obj]
    return obj


def secms_toolkit(ctx):
    from toolkits.secms import pipeline
    from .. import build
    p = ctx.params("secms_toolkit")
    comb = p.pop("combine", None)                  # not a toolkit run key: the combination spec (module doc)
    run_cfg = _frozen_paths(p, ctx)
    res = pipeline.run(run_cfg, resolve=ctx.resolve)
    if comb:
        from toolkits.secms import combine
        spec = {k: v for k, v in comb.items() if k != "min_fractions"}
        if not spec.get("label"):
            raise ValueError("inventory.params.secms_toolkit.combine.label is required")
        have = {str(d["source_id"]) for d in run_cfg["datasets"]}
        miss = [s for s in combine.spec_ids(spec) if s not in have]
        if miss:
            raise ValueError(f"secms_toolkit.combine: datasets {miss} are not in secms_toolkit.datasets {sorted(have)}")
        res = combine.combine_tables(res, str(spec["label"]), int(comb.get("min_fractions", 0) or 0), spec)
        res["_diagnostics"]["combine"] = comb
    for k in TABLES:
        ctx.write_tsv(res[k], f"secms/{k}.tsv")
    d = res["_diagnostics"]
    out = ctx.out("secms/diagnostics.json")
    build.atomic_save(out, lambda t: open(t, "w").write(pipeline.diagnostics_json(d)))
    for ds in d["datasets"]:
        for run, r in ds["runs"].items():
            c = r["calibration"]
            ctx.say(f"secms {ds['source_id']} [{run}]: calibration {c['used']} (n {c['n_kept']}, rho {c['spearman_rho']:.2f}, "
                    f"x{c['resid_fold']:.2f}); proteins scored {r['n_proteins_scored']}, complexes measured "
                    f"{r['n_complexes_measured']}, pairs {r['n_pairs']}")
    if comb:
        P = res["secms_protein"]
        P = P[(P.source_id == str(comb["label"])) & (P.threshold_k == P.threshold_k.min())]
        ctx.say(f"secms combined '{comb['label']}': {len(P)} proteins; confidence "
                + ", ".join(f"{k} {v}" for k, v in P.confidence.value_counts().items()))


secms_toolkit.OUTPUTS = [f"secms/{k}.tsv" for k in TABLES] + ["secms/diagnostics.json"]
register("adapter.secms", "secms_toolkit")(secms_toolkit)


# ------------------------------------------------------------------------------------------------ view for assembly
class SecView:
    """the SEC-MS tables of one run (source x condition x replicate) as lookups for assembly-state rules"""

    def __init__(self, ctx, sel: dict | None = None):
        sel = sel or {}
        P = pd.read_csv(ctx.need("secms/secms_protein.tsv"), sep="\t", dtype={"replicate": str, "condition": str})
        self.source_id = sel.get("source_id") or _only(P.source_id, "source_id")
        P = P[P.source_id == self.source_id]
        self.condition = sel.get("condition") or _only(P.condition, "condition")
        P = P[P.condition == self.condition]
        reps = list(dict.fromkeys(P.replicate))
        self.replicate = str(sel.get("replicate") or ("mean" if "mean" in reps else reps[0]))
        self._sel = lambda T: T[(T.source_id == self.source_id) & (T.condition.astype(str) == self.condition) &
                                (T.replicate.astype(str) == self.replicate)]
        self.protein = P[P.replicate == self.replicate]
        self.thresholds = self.protein.pivot(index="uniprot", columns="threshold_k", values="assembled_frac")
        self.ctx = ctx
        self._pairs = self._cpx = self._sub = None

    def a(self, acc, k, default=np.nan):
        return float(self.thresholds.at[acc, k]) if acc in self.thresholds.index and k in self.thresholds.columns else default

    def has(self, acc):
        return acc in self.thresholds.index

    def legacy_sec_assembly(self, k=2):
        """sec_assembly.tsv view (legacy columns): uniprot, gene, mono_kda, apex_frac, apex_app_kda, assembled_frac, tm"""
        p = self.protein[self.protein.threshold_k == k]
        return pd.DataFrame(dict(uniprot=p.uniprot, gene=p.gene, mono_kda=p.mono_kda, apex_frac=p.apex_fraction,
                                 apex_app_kda=p.apex_app_kda, assembled_frac=p.assembled_frac, tm=p.tm)).reset_index(drop=True)

    def calibrants(self):
        C = pd.read_csv(self.ctx.need("secms/secms_calibrants.tsv"), sep="\t", dtype={"replicate": str, "condition": str})
        return self._sel(C).reset_index(drop=True)

    def pair_r(self):
        if self._pairs is None:
            T = pd.read_csv(self.ctx.need("secms/secms_pairs.tsv"), sep="\t", dtype={"replicate": str, "condition": str})
            T = self._sel(T)
            self._pairs = dict(zip(zip(T.uniprot_a, T.uniprot_b), T.r))
        d = self._pairs

        def f(a, b):
            v = d.get((a, b))
            return float(v) if v is not None else float(d.get((b, a), np.nan))
        return f

    def complexes(self):
        if self._cpx is None:
            C = pd.read_csv(self.ctx.need("secms/secms_complex.tsv"), sep="\t", dtype={"replicate": str, "condition": str})
            S = pd.read_csv(self.ctx.need("secms/secms_complex_subunit.tsv"), sep="\t", dtype={"replicate": str, "condition": str})
            self._cpx = self._sel(C).drop_duplicates("complex_id").set_index("complex_id")
            S = self._sel(S)
            S = S[S.detected.astype(str) == "True"]
            self._sub = {c: dict(zip(g.uniprot, g.in_complex_frac)) for c, g in S.groupby("complex_id")}
        return self._cpx, self._sub


def _only(s, what):
    v = list(dict.fromkeys(s.astype(str)))
    if len(v) != 1:
        raise ValueError(f"SEC-MS tables hold several {what}s {v}: set assembly_options.secms.{what}")
    return v[0]


# ------------------------------------------------------------------------------------------------ rules
def homo_share(view: SecView | None, acc: str, state: int, o: dict):
    """(share, k used): share of the free pool put into the homo-oligomer"""
    h = o["homo_oligomer_share"]
    if h["source"] == "constant":
        return float(h["value"]), None
    k = None
    for lim, kk in h["k_by_state"]:
        if lim is None or state <= lim:
            k = kk
            break
    if view is not None and view.has(acc):
        return float(view.thresholds.at[acc, k]), k
    return float(h["default"]), k


def homo_apex_ok(acc: str, state: int, sa: pd.DataFrame, ak: float, o: dict) -> bool:
    """legacy check: a protein that is SEC-assembled (ak = a_k at apex_check_min_k > apex_check_min_frac) must have its
    SEC apex at state x monomer within apex_ratio"""
    h = o["homo_oligomer_share"]
    if acc not in sa.index:
        return True
    ratio = sa.apex_app_kda[acc] / (state * sa.mono_kda[acc])
    lo, hi = h["apex_ratio"]
    return not (not (lo <= ratio <= hi) and ak > h["apex_check_min_frac"])


def homodimer_share(view: SecView | None, acc: str, o: dict):
    h = o["afdb_homodimer_share"]
    if h["source"] == "constant":
        return float(h["value"])
    return float(view.thresholds.at[acc, h["k"]]) if view is not None and view.has(acc) else float(h["default"])


def measured_copies(view: SecView, cpx_id: str, present: dict, pool: dict, acc_of: dict, o: dict):
    """copies of a complex from SEC-measured free pools (None -> not measured: caller keeps the cap rule).
    present: {gene: stoich} of subunits with copies > 0; pool: gene -> copies; acc_of: gene -> UniProt"""
    m = o["measured_free_pool"]
    C, S = view.complexes()
    if cpx_id not in C.index or cpx_id not in S:
        return None, "no SEC row"
    row = C.loc[cpx_id]
    if cpx_id in (m.get("exclude_complex_ids") or []):
        return None, "excluded by config"
    if m.get("max_expected_kda") and float(row.expected_kda) > float(m["max_expected_kda"]):
        return None, f"expected {float(row.expected_kda):.0f} kDa > max_expected_kda"
    if m.get("exclude_tm", True) and str(row.tm_any) == "True":
        return None, "TM complex (detergent micelles; SEC under-recovers)"
    sub = S[cpx_id]
    nd = sum(1 for g in present if acc_of.get(g) in sub)
    md, mf = m.get("min_subunits_detected"), m.get("min_detected_frac")
    if md is None and mf is None:
        ok = str(row.measured) == "True"
    else:
        ok = nd >= int(md or 1) and nd / max(int(row.n_subunits), 1) >= float(mf or 0)
    if not ok:
        return None, f"coverage {nd}/{int(row.n_subunits)}"
    vals = [pool[g] * float(sub[acc_of[g]]) / s for g, s in present.items()
            if acc_of.get(g) in sub and np.isfinite(float(sub[acc_of[g]]))]
    if not vals:
        return None, "no finite subunit fraction"
    return float(np.median(vals)), f"measured ({len(vals)} subunits)"


def canonical_tables():
    """canonical SEC-MS outputs (inventory/<name>) and their columns"""
    return {"secms_protein.tsv": schema.SECMS_PROTEIN + schema.SECMS_PROTEIN_EXTRA,
            "secms_complex.tsv": schema.SECMS_COMPLEX + schema.SECMS_COMPLEX_EXTRA,
            "secms_complex_subunit.tsv": schema.SECMS_COMPLEX_SUBUNIT,
            "secms_pairs.tsv": schema.SECMS_PAIRS,
            "secms_calibration.tsv": schema.SECMS_CALIBRATION}


def to_json(o):
    return json.dumps(o, default=str)
