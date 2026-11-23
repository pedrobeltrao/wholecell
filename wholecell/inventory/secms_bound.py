"""Per-protein SEC-MS bound fraction vs the modelled assembly state (canonical table `secms_bound_fraction.tsv`;
spec: docs/reference/inventory_schema.md §4c; methods: docs/methods_secms_bound_fraction.md).

SEC-MS measures more assembly than the structural model represents: many proteins elute far above their monomer mass
(EF-1A, Hsp70s, PAB1, ...) but have no structure of the assembly, so the model keeps them monomeric. Structure-less
assemblies are NOT invented (packing unchanged, Pedro 2026-11-19); instead every protein carries, as an attribute,
  measured   secms_bound_fraction = a_k of the SEC-MS run the assembly rules read (fraction of the protein's profile at
             apparent mass >= k x monomer), with its uncertainty (dataset / replicate SD, calibration band) and datasets;
             NaN (na_reason 'unmeasured ...') when no dataset can resolve k x monomer (toolkits/secms/protein.py
             measurable range; 2026-11-23: such proteins previously scored 1 by construction)
  modelled   modelled_assembled_fraction = 1 - copies in monomer species / protein copies (v0);
  gap        unrepresented_bound_fraction = max(0, measured - modelled), and the same as a share of the modelled free
             (monomer) copies, unrepresented_share_of_free = gap / (1 - modelled), the quantity a simulation can apply to
             each monomer molecule (effective size / reduced mobility of a bound share; not implemented here).

Selection (defaults = what the assembly rules use; config inventory.canonical.secms_bound_fraction overrides):
  source_id / condition / replicate   the SEC-MS run of the assembly options (assembly/assembly_steps.json `secms`, else
                                      the assembly adapter's assembly_options, else the combine label, else the only run)
  k: rules | <int>                    rules (default): k_by_state of the homo-oligomer share for a protein modelled in a
                                      homo-oligomer of state n (n <= 3 -> 2, <= 5 -> 3, else 4 in the v7 rules), else the
                                      dimer k of the rules (secms_fractions.threshold_k / AF homodimer share k; 2)
Reads only build tables (pack_table, v0, assembly options, oligo_map.json when present, secms/secms_protein.tsv); writes
nothing but its own table (the caller writes it). Never changes copies, species or any other table.
"""
from __future__ import annotations

import json
import os

import numpy as np
import pandas as pd

from . import schema

AV = schema.AVOGADRO


def _tokens(g) -> list:
    return [t.strip() for t in str(g).split(";") if t.strip()]


def _key(g) -> str:
    return ";".join(sorted(_tokens(g)))


def _get(d, *path, default=None):
    for p in path:
        if not isinstance(d, dict) or p not in d:
            return default
        d = d[p]
    return d


def rule_options(ctx) -> dict:
    """SEC-MS run selection and k rules of the assembly step of this build: {sel, k_by_state, k_default, basis}"""
    from . import secms as S2
    p = ctx.path("assembly/assembly_steps.json")
    if os.path.exists(p):
        o = json.load(open(p)).get("options", {}) or {}
        sel = dict(o.get("secms") or {})
        kbs = (_get(o, "pdb_homo_oligomers", "share", "k_by_state")
               or S2.ASSEMBLY_DEFAULTS["homo_oligomer_share"]["k_by_state"])
        k0 = (_get(o, "secms_fractions", "threshold_k") or _get(o, "afdb_homodimers", "share", "k")
              or _get(o, "pdb_homo_oligomers", "af_share", "k") or S2.ASSEMBLY_DEFAULTS["afdb_homodimer_share"]["k"])
        basis = "assembly/assembly_steps.json"
    else:
        ad = (ctx.inv.get("adapters", {}) or {}).get("assembly") or ""
        o = S2.assembly_options(ctx.params(ad) if ad else {})
        sel = dict(o["secms"])
        kbs = o["homo_oligomer_share"]["k_by_state"]
        k0 = o["afdb_homodimer_share"]["k"]
        basis = f"assembly_options of adapter {ad or 'none'} (secms.ASSEMBLY_DEFAULTS merged)"
    return dict(sel=sel, k_by_state=kbs, k_default=int(k0), basis=basis)


def k_for_state(n, k_by_state, k_default):
    if n is None or not np.isfinite(n) or n < 2:
        return k_default
    for lim, kk in k_by_state:
        if lim is None or n <= lim:
            return int(kk)
    return k_default


def model_state(P: pd.DataFrame, oligo_map: dict | None = None) -> dict:
    """per gene token: monomer copies (keyed by the normalised gene label), assembly species naming it, homo state.
    Monomer species = kind protein that is not a library-level oligomer (oligo_map.json keys of ingredient-pool
    assembly presets, e.g. the yeast X5 `<gene>_oligo` / `<gene>__homo` rows)."""
    om = oligo_map or {}
    is_mono = (P.kind == "protein") & ~P.ing_id.isin(list(om))
    mono = P[is_mono].assign(_k=lambda d: d.genes.map(_key)).groupby("_k").copies.sum().to_dict()
    asm, homo = {}, {}
    for r in P[~is_mono].itertuples():
        toks = _tokens(r.genes)
        for t in toks:
            asm.setdefault(t, []).append((float(r.copies), r.ing_id))
        if r.ing_id in om:
            n = float(om[r.ing_id][1])
        elif len(set(toks)) == 1:
            try:
                n = float(sum(float(s) for s in str(r.stoich).split(";")))
            except ValueError:
                n = np.nan
        else:
            n = np.nan
        if len(set(toks)) == 1 and np.isfinite(n) and n > 1 and float(r.copies) > 0:
            homo[toks[0]] = max(homo.get(toks[0], 0.0), n)
    return dict(mono=mono, asm=asm, homo=homo)


def _select_run(ctx, sel: dict):
    """the SEC-MS rows of one run (SecView semantics); None + reason when unavailable"""
    from . import secms as S2
    p = ctx.path("secms/secms_protein.tsv")
    if not (ctx.inv.get("adapters", {}) or {}).get("secms") or not os.path.exists(p):
        return None, "no SEC-MS step in this build (inventory.adapters.secms not set)"
    sel = {k: v for k, v in (sel or {}).items() if v is not None}
    if not sel.get("source_id"):
        lab = (_get(ctx.params("secms_toolkit"), "combine", "label"))
        if lab:
            sel["source_id"] = str(lab)
    try:
        return S2.SecView(ctx, sel), ""
    except ValueError as e:                       # several runs and no selection: NA, never a failed canonical step
        return None, f"SEC-MS run not selected ({e}); set inventory.canonical.secms_bound_fraction.source_id"


def compute(ctx, P: pd.DataFrame, V0: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """the canonical secms_bound_fraction table (one row per v0 gene) and a summary dict"""
    can = (ctx.inv.get("canonical", {}) or {}).get("secms_bound_fraction", {}) or {}
    ro = rule_options(ctx)
    sel = dict(ro["sel"])
    sel.update({k: can[k] for k in ("source_id", "condition", "replicate") if can.get(k)})
    kmode = can.get("k", "rules")
    omp = ctx.path("assembly/oligo_map.json")
    om = json.load(open(omp)) if os.path.exists(omp) else {}
    ms = model_state(P, om)
    view, why = _select_run(ctx, sel)
    rows_by = {}
    if view is not None:
        T = view.protein.copy()
        T["threshold_k"] = T.threshold_k.astype(int)
        rows_by = {(u, int(k)): r for u, k, r in zip(T.uniprot, T.threshold_k, T.to_dict("records"))}
        have_u = set(T.uniprot)
        ks = sorted(set(T.threshold_k))
    out = []
    for g, r in V0.iterrows():
        g = str(g)
        u = r.get("uniprot")
        u = u if isinstance(u, str) and u else ""
        copies = float(r.copies) if pd.notna(r.copies) else 0.0
        mda = float(r.mass_da) if pd.notna(r.mass_da) else np.nan
        toks = _tokens(g)
        free = float(ms["mono"].get(_key(g), 0.0))
        asm = sorted({(c, i) for t in toks for c, i in ms["asm"].get(t, [])}, key=lambda x: (-x[0], x[1]))
        asm = [x for x in asm if x[0] > 0]
        n_homo = max([ms["homo"].get(t, 0.0) for t in toks] or [0.0])
        if copies <= 0:
            state, rep, m = "absent", None, np.nan
        elif asm:
            state, rep = "in_assembly", True
        elif free > 0:
            state, rep = "monomer_only", False
        else:
            state, rep = "not_in_pack", True           # copies held by a family / aggregate / pre-pass species
        if copies > 0:
            m = float(min(max(1.0 - free / copies, 0.0), 1.0))
        if kmode == "rules":
            k = k_for_state(n_homo if n_homo > 1 else None, ro["k_by_state"], ro["k_default"])
            kb = f"rules: homo-oligomer state {int(n_homo)}" if n_homo > 1 else f"rules: default k {ro['k_default']}"
        else:
            k, kb = int(kmode), "config inventory.canonical.secms_bound_fraction.k"
        d = dict(gene=g, uniprot=u, copies=copies, mass_da=mda, mass_pg=copies * mda / AV * 1e12 if np.isfinite(mda) else np.nan,
                 secms_source_id="", secms_condition="", secms_replicate="", datasets="", n_datasets=np.nan, confidence="",
                 threshold_k=k, k_basis=kb, secms_bound_fraction=np.nan, secms_bound_fraction_sd=np.nan,
                 secms_bound_fraction_lo=np.nan, secms_bound_fraction_hi=np.nan, apex_app_kda=np.nan, mono_kda=np.nan,
                 model_state=state, structure_represented=rep,
                 model_species=";".join(i for _, i in asm[:5]) + (f";+{len(asm) - 5} more" if len(asm) > 5 else ""),
                 model_homo_state=n_homo if n_homo > 1 else np.nan, modelled_free_copies=free,
                 modelled_assembled_fraction=m, unrepresented_bound_fraction=np.nan,
                 unrepresented_share_of_free=np.nan, na_reason="")
        if view is None:
            d["na_reason"] = why
        else:
            d.update(secms_source_id=view.source_id, secms_condition=view.condition, secms_replicate=view.replicate)
            if not u:
                d["na_reason"] = "no UniProt accession"
            elif u not in have_u:
                d["na_reason"] = f"not quantified in SEC-MS run {view.source_id}"
            elif (u, k) not in rows_by:
                d["na_reason"] = f"no a_k at k {k} (run thresholds {ks})"
            else:
                s = rows_by[(u, k)]
                conds = s.get("conditions")
                if isinstance(conds, str) and conds:
                    ds, nd = conds, s.get("n_conditions")
                else:
                    ds, nd = f"{view.source_id}|{view.condition}|{view.replicate}", s.get("n_rep", 1)
                a = float(s["assembled_frac"]) if pd.notna(s["assembled_frac"]) else np.nan
                d.update(datasets=ds, n_datasets=nd, confidence=s.get("confidence") if isinstance(s.get("confidence"), str) else "",
                         secms_bound_fraction=a, secms_bound_fraction_sd=s.get("assembled_frac_sd", np.nan),
                         secms_bound_fraction_lo=s.get("assembled_frac_lo", np.nan),
                         secms_bound_fraction_hi=s.get("assembled_frac_hi", np.nan),
                         apex_app_kda=s.get("apex_app_kda", np.nan), mono_kda=s.get("mono_kda", np.nan))
                if not np.isfinite(a):
                    nt = s.get("notes") if isinstance(s.get("notes"), str) else ""
                    d["na_reason"] = (f"unmeasured at k {k}: k x monomer outside the calibrated range of SEC-MS run "
                                      f"{view.source_id}" if "unmeasured" in nt else f"a_k NaN in SEC-MS run {view.source_id}")
                elif copies > 0:
                    gap = max(0.0, a - m)
                    d["unrepresented_bound_fraction"] = gap
                    d["unrepresented_share_of_free"] = min(gap / (1.0 - m), 1.0) if m < 1.0 else 0.0
        if copies <= 0 and not d["na_reason"]:
            d["na_reason"] = "no copies in this cell (model columns NA)"
        out.append(d)
    F = pd.DataFrame(out, columns=schema.SECMS_BOUND_FRACTION)
    return F, summary(F, ro, sel if view is None else dict(source_id=view.source_id, condition=view.condition,
                                                            replicate=view.replicate), kmode)


def summary(F: pd.DataFrame, ro=None, sel=None, kmode="rules", top=20) -> dict:
    """mass shares (of the cell's protein mass): modelled assembly, measured bound, measured bound but unrepresented"""
    X = F[F.copies > 0]
    mass = X.mass_pg.fillna(0.0)
    tot = float(mass.sum())
    cov = X.secms_bound_fraction.notna()
    mc = float(mass[cov].sum())
    unrep_pg = (X.unrepresented_bound_fraction.fillna(0.0) * mass)
    T = X.assign(unrep_pg=unrep_pg).sort_values(["unrep_pg", "gene"], ascending=[False, True]).head(top)
    w = lambda c, m=None: float((X[c].fillna(0.0) * mass)[m if m is not None else slice(None)].sum())  # noqa: E731
    return dict(
        selection=sel, k=kmode, rule_options=ro,
        total_pg=tot, secms_covered_pg=mc, secms_covered_frac=mc / tot if tot else np.nan,
        n_proteins=int(len(X)), n_measured=int(cov.sum()),
        n_unmeasured_out_of_range=int(X.na_reason.astype(str).str.startswith("unmeasured").sum()),
        measured_bound_protein_mean=float(X.secms_bound_fraction[cov].mean()) if cov.any() else np.nan,
        measured_bound_protein_median=float(X.secms_bound_fraction[cov].median()) if cov.any() else np.nan,
        mass_structure_represented_frac=float(mass[X.structure_represented == True].sum()) / tot if tot else np.nan,  # noqa: E712
        mass_monomer_only_frac=float(mass[X.model_state == "monomer_only"].sum()) / tot if tot else np.nan,
        modelled_assembled_frac=w("modelled_assembled_fraction") / tot if tot else np.nan,
        modelled_assembled_frac_of_covered=w("modelled_assembled_fraction", cov) / mc if mc else np.nan,
        measured_bound_frac_of_covered=w("secms_bound_fraction", cov) / mc if mc else np.nan,
        unrepresented_bound_frac=float(unrep_pg.sum()) / tot if tot else np.nan,
        unrepresented_bound_frac_of_covered=float(unrep_pg.sum()) / mc if mc else np.nan,
        unrepresented_bound_pg=float(unrep_pg.sum()),
        unrepresented_in_monomer_only_pg=float(unrep_pg[X.model_state == "monomer_only"].sum()),
        monomer_only_measured_bound_of_covered=(w("secms_bound_fraction", cov & (X.model_state == "monomer_only"))
                                                / float(mass[cov & (X.model_state == "monomer_only")].sum())
                                                if float(mass[cov & (X.model_state == "monomer_only")].sum()) else np.nan),
        top=[dict(gene=r.gene, uniprot=r.uniprot, mass_pg=round(float(r.mass_pg), 5), model_state=r.model_state,
                  measured=round(float(r.secms_bound_fraction), 3), k=int(r.threshold_k),
                  modelled=round(float(r.modelled_assembled_fraction), 3), unrep_pg=round(float(r.unrep_pg), 5),
                  unrep_frac_of_cell=round(float(r.unrep_pg) / tot, 5) if tot else None) for r in T.itertuples()
             if r.unrep_pg > 0])
