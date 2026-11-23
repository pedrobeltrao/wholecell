"""The assembly steps: one plain function per step, each `f(st)` on the shared State (pools, options, inputs).

origin: plugins/jurkat/inventory/assembly.py (port of legacy jurkat/inventory/{complexes_v0, sec_assembly,
homo_oligomers, afdb_dimers}.py, 2026-11-15/16); moved to the core 2026-11-16 (unification A4). The arithmetic is the
legacy one line for line (preset legacy_v7 reproduces the v7 tables byte for byte); the constants are preset options and
the Jurkat lists are cell data (inventory.params.assembly_steps.data). Later steps read the earlier steps' tables from
inventory/assembly/ as v7 did (same float round trips).
"""
from __future__ import annotations

import collections
import gzip
import hashlib
import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from .. import secms as S2


def read_curated(path: str) -> pd.DataFrame:
    """curated TSV (repo-relative or absolute path): lines starting with '#' are provenance comments"""
    from ...config import REPO
    p = path if os.path.isabs(path) else os.path.join(REPO, path)
    if not os.path.exists(p):
        raise FileNotFoundError(f"curated table {path} not found (looked at {p})")
    import io
    txt = "".join(line for line in open(p) if not line.startswith("#"))
    return pd.read_csv(io.StringIO(txt), sep="\t", dtype=str, keep_default_na=False)


class State:
    """what the steps share: ctx, options o, inputs / data of the cell, the SEC-MS view, and the pools of the complex
    steps (pool0 = input copies, pool = totals after curated rescaling, left = copies not in a complex)"""

    def __init__(self, ctx, params, o, view):
        self.ctx, self.o, self.view = ctx, o, view
        self.inputs = dict(params.get("inputs") or {})
        self.data = dict(params.get("data") or {})
        self.workers = int(params.get("workers", 8))
        self.inv = pd.read_csv(ctx.need("v0/inventory_v0.tsv"), sep="\t", index_col="gene")
        self.pool0 = self.inv.copies.to_dict()
        self.pool = dict(self.pool0)
        self.left = dict(self.pool0)
        self.complexes = None            # DataFrame of Complex Portal rows (complex_portal)
        self.curated_rows = []           # rows appended to complexes_v0 (curated machines)
        self.skip_name_regex = None      # Complex Portal names built by curated_machines instead
        self.decisions = []              # complex_copies_rule rows
        self.free = None                 # free pool of the dimer steps (pd.Series)
        self._pdb = None

    def frozen(self, key):
        return self.ctx.resolve(self.inputs[key])

    def pdb_local(self):
        """local PDB assembly1 folder and its entry ids; the folder listing is checked against inputs.pdb_listing_md5"""
        if self._pdb is None:
            A = self.frozen("pdb_assemblies")
            names = os.listdir(A)
            md5 = hashlib.md5("\n".join(sorted(names)).encode()).hexdigest()
            want = self.inputs.get("pdb_listing_md5")
            if want and md5 != want:
                raise RuntimeError(f"PDB assembly folder listing changed (md5 {md5} != inputs.pdb_listing_md5 {want}); the "
                                   f"assembly state would differ from the declared frozen set")
            self._pdb = (A, {f.split("-")[0] for f in names})
        return self._pdb


def legacy_options(o: dict) -> dict:
    """the SEC-MS rule options in the layout of wholecell/inventory/secms.py (ASSEMBLY_DEFAULTS); written to
    assembly/assembly_options.json and passed to the secms.py rule helpers"""
    cp, hx = o["complex_portal"], o["afdb_heterodimers"]
    return {"secms": o["secms"], "complex_copies": cp["copies"], "measured_free_pool": cp["measured_free_pool"],
            "homo_oligomer_share": o["pdb_homo_oligomers"]["share"],
            "afdb_homodimer_share": o["afdb_homodimers"]["share"],
            "heterodimer": {"sec_r_min": hx["evidence"]["sec_r_min"], "string_min": hx["evidence"]["string_min"],
                            "bound_fraction": hx["bound_fraction"]}}


# ------------------------------------------------------------------------------------------------ curated machines
def curated_machines(st):
    """nucleosome (option curated_machines.nucleosome = dna_count): N = genome_bp / bp_per_nucleosome; each histone
    family (data.histone_families) is rescaled to 2 N copies and consumed; Complex Portal rows whose name matches
    data.nucleosome.complex_name_regex are then skipped by complex_portal"""
    if st.o["curated_machines"]["nucleosome"] == "none":
        return
    ctx, pool, left = st.ctx, st.pool, st.left
    fam = dict(zip(*[read_curated(st.data["histone_families"])[c] for c in ("family", "gene_regex")]))
    fam_genes = {f: [g for g in pool if re.match(pat, g)] for f, pat in fam.items()}
    fam_tot = {f: sum(pool[g] for g in gs) for f, gs in fam_genes.items()}
    NUC_IBAQ = int(min(fam_tot.values()) // 2)
    nz = st.data["nucleosome"]
    NUC = int(nz["genome_bp"] / nz["bp_per_nucleosome"])
    for f, gs in fam_genes.items():
        k = 2 * NUC / fam_tot[f]
        for g in gs:
            pool[g] = int(round(pool[g] * k)); left[g] = pool[g]
        fam_tot[f] = 2 * NUC
    ctx.say(f"histone iBAQ supports {NUC_IBAQ:.3g} nucleosomes; DNA-based value used: {NUC:.3g}")
    for f, gs in fam_genes.items():
        for g in gs:
            left[g] = max(left[g] - round(2 * NUC * pool[g] / fam_tot[f]), 0)
    mk = nz["mass_kda"]
    st.curated_rows.append(dict(cpx="CURATED-NUCLEOSOME", name="Nucleosome core particle (generic, family totals)",
                                n_types=4, subunits="H2A:2;H2B:2;H3:2;H4:2", copies_rule=NUC, copies_built=NUC,
                                mass_kda=mk, n_zero_subunits=0, best_pdb=nz["pdb"], pdb_mass_cov=1.0,
                                mass_pg=NUC * mk * 1e3 / 6.022e23 * 1e12))
    st.skip_name_regex = nz["complex_name_regex"]


# ------------------------------------------------------------------------------------------------ Complex Portal
def _parse_complex_portal(st):
    """Complex Portal rows mapped to inventory genes, the copies rule on the input pools and the best local PDB entry"""
    ctx, inv = st.ctx, st.inv
    acc2gene = dict(zip(inv.uniprot, inv.index)); pool = st.pool0; mass = inv.mass_da.fillna(0).to_dict()
    so = st.o["complex_portal"]["structure"]
    cp = pd.read_csv(st.frozen("complex_portal"), sep="\t")
    sifts = pd.read_csv(st.frozen("sifts"), sep="\t", comment="#", usecols=["PDB", "CHAIN", "SP_PRIMARY"])
    _, local = st.pdb_local()
    sifts = sifts[sifts.PDB.isin(local)]
    pdb_up = sifts.groupby("PDB").SP_PRIMARY.agg(lambda s: set(s)).to_dict()
    up_pdb = defaultdict(set)
    for pp, ups in pdb_up.items():
        for u in ups:
            up_pdb[u].add(pp)
    # tie-break reference for equal-coverage PDB candidates (exact reproduction of an earlier build; else sorted order)
    tiebreak = {}
    if so["tie_break"] == "fewest_entities_then_reference" and st.inputs.get("pdb_tiebreak_reference"):
        tb = pd.read_csv(st.frozen("pdb_tiebreak_reference"), sep="\t", usecols=["cpx", "best_pdb"])
        tiebreak = dict(zip(tb.cpx, tb.best_pdb))
    n_tb = 0
    rows = []
    for r in cp.itertuples(index=False):
        parts = re.findall(r"([A-Z][A-Z0-9]{5,9})(?:-\w+)?\((\d+)\)", str(r[18]))
        sub, unmapped = {}, 0
        for a, s in parts:
            g = acc2gene.get(a)
            if g is None:
                unmapped += 1; continue
            sub[g] = sub.get(g, 0) + max(1, int(s))
        if len(sub) < 1 or (len(sub) == 1 and list(sub.values())[0] < 2):
            continue
        ratios = np.array([pool[g] / s for g, s in sub.items()])
        med = float(np.median(ratios))
        lim = float(np.percentile(ratios, 25)) if len(sub) >= 4 else float(ratios.min())
        copies = min(med, 2 * lim)
        m = sum(mass[g] * s for g, s in sub.items())
        accs = {inv.uniprot[g] for g in sub}
        cand = set().union(*[up_pdb.get(a, set()) for a in accs]) if accs else set()
        best, bcov = None, 0.0
        covs = {pp: sum(mass[g] * s for g, s in sub.items() if inv.uniprot[g] in pdb_up[pp]) / max(m, 1) for pp in cand}
        for pp in sorted(cand):                      # legacy iterated a set (random hash order): ties -> see below
            cov = covs[pp]
            if cov > bcov + 1e-9 or (abs(cov - bcov) < 1e-9 and best and len(pdb_up[pp]) < len(pdb_up[best])):
                best, bcov = pp, cov
        ref = tiebreak.get(r[0])
        if best is not None and isinstance(ref, str) and ref in covs and ref != best:
            # tied set = max coverage (1e-9) and fewest UniProt entries; the legacy pick among ties was arbitrary
            if abs(covs[ref] - bcov) < 1e-9 and len(pdb_up[ref]) == len(pdb_up[best]):
                best, bcov = ref, covs[ref]
                n_tb += 1
        if so["min_coverage"] is not None and (best is None or bcov < float(so["min_coverage"])):
            continue
        rows.append(dict(cpx=r[0], name=r[1], n_types=len(sub), unmapped=unmapped,
                         subunits=";".join(f"{g}:{s}" for g, s in sub.items()), copies_rule=round(copies),
                         limiting=min(sub, key=lambda g: pool[g] / sub[g]), mass_kda=round(m / 1e3, 1),
                         n_zero_subunits=int(sum(pool[g] == 0 for g in sub)), best_pdb=best, pdb_mass_cov=round(bcov, 2),
                         cp_pdb=";".join(re.findall(r"wwpdb:(\w{4})", str(r[8])))))
    C = pd.DataFrame(rows)
    ctx.say(f"PDB candidate ties resolved by the tie-break reference: {n_tb}")
    C["mass_pg"] = C.copies_rule * C.mass_kda * 1e3 / 6.022e23 * 1e12
    return C


def _cap_rule(vals, ntypes):
    vals = np.asarray(vals, float); med = np.median(vals)
    lim = np.percentile(vals, 25) if ntypes >= 4 else vals.min()
    return max(min(med, 2 * lim), 0)


def complex_portal(st):
    """greedy build: data.complex_priority first, then by mass; copies by option complex_portal.copies, capped by
    data.complex_copy_caps; complexes with >= max_zero_mass_frac of their mass in absent subunits are not built"""
    ctx, o = st.ctx, st.o["complex_portal"]
    lo = legacy_options(st.o)
    C = _parse_complex_portal(st)
    pool, left, mass = st.pool, st.left, st.inv.mass_da.fillna(0).to_dict()
    PRIORITY = list(read_curated(st.data["complex_priority"]).complex_id) if st.data.get("complex_priority") else []
    caps = {}
    if st.data.get("complex_copy_caps"):
        T = read_curated(st.data["complex_copy_caps"])
        caps = dict(zip(T.complex_id, T.copies.astype(int)))
    order = sorted(C.index, key=lambda i: (PRIORITY.index(C.cpx[i]) if C.cpx[i] in PRIORITY else len(PRIORITY), -C.mass_pg[i]))
    skip = C.name.str.contains(st.skip_name_regex, regex=True) if st.skip_name_regex else pd.Series(False, C.index)
    acc_of = st.inv.uniprot.to_dict()
    measured_on = o["copies"] == "measured_free_pool"
    zmax = float(o["max_zero_mass_frac"])
    built = {}
    for i in order:
        r = C.loc[i]; sub = {x.split(":")[0]: int(x.split(":")[1]) for x in r.subunits.split(";")}
        zero_mass = sum(mass[g] * s for g, s in sub.items() if pool[g] == 0) / max(sum(mass[g] * s for g, s in sub.items()), 1)
        present = {g: s for g, s in sub.items() if pool[g] > 0}
        if zero_mass >= zmax or not present or (len(sub) == 1 and list(sub.values())[0] < 2):
            built[i] = 0; continue
        if skip[i]:
            built[i] = 0; continue
        n_cap = int(_cap_rule([left[g] / s for g, s in present.items()], len(present)))
        n_meas, why = S2.measured_copies(st.view, r.cpx, present, pool, acc_of, lo) if st.view is not None else (None, "no SEC")
        avail = int(min(left[g] // s for g, s in present.items()))
        n = n_cap
        if measured_on and n_meas is not None:
            n = max(min(int(n_meas), avail), 0)
        if r.cpx in caps:
            n = min(n, caps[r.cpx])
        for g, s in present.items():
            left[g] = max(left[g] - n * s, 0)
        built[i] = n
        st.decisions.append(dict(cpx=r.cpx, n_cap=n_cap, n_measured=round(n_meas, 1) if n_meas is not None else None,
                                 available=avail, copies_built=n,
                                 rule="measured_free_pool" if (measured_on and n_meas is not None) else "stoichiometry_cap",
                                 measured_note=why))
    C["copies_built"] = pd.Series(built)
    st.complexes = C


def write_complex_tables(st):
    """complexes_v0.tsv (Complex Portal rows + curated machines), inventory_v0_complexed.tsv (pools after the complex
    steps), complex_copies_rule.tsv (cap vs SEC-measured copies per built complex)"""
    ctx, inv = st.ctx, st.inv.copy()
    C = st.complexes if st.complexes is not None else pd.DataFrame(columns=["cpx", "name", "n_types", "subunits",
                                                                            "copies_rule", "mass_kda", "copies_built"])
    if st.curated_rows:
        C = pd.concat([C, pd.DataFrame(st.curated_rows)], ignore_index=True)
    C["copies_built"] = C.copies_built.astype(int)
    C["mass_built_pg"] = C.copies_built * C.mass_kda * 1e3 / 6.022e23 * 1e12
    ctx.write_tsv(C.sort_values("mass_built_pg", ascending=False), "assembly/complexes_v0.tsv")
    inv["copies"] = pd.Series(st.pool).reindex(inv.index).fillna(inv.copies).astype(np.int64)
    inv["mass_pg"] = inv.copies * inv.mass_da.fillna(0) / 6.022e23 * 1e12
    inv["copies_in_complexes"] = inv.copies - pd.Series(st.left).reindex(inv.index).fillna(inv.copies)
    inv["frac_in_complexes"] = (inv.copies_in_complexes / inv.copies.replace(0, np.nan)).round(3)
    ctx.write_tsv(inv, "assembly/inventory_v0_complexed.tsv", index=True)
    D = pd.DataFrame(st.decisions, columns=["cpx", "n_cap", "n_measured", "available", "copies_built", "rule", "measured_note"])
    D["mass_kda"] = D.cpx.map(dict(zip(C.cpx, C.mass_kda)))
    ctx.write_tsv(D, "assembly/complex_copies_rule.tsv")
    ctx.say(f"complex copies: rule {st.o['complex_portal']['copies']}; SEC-measured complexes {D.n_measured.notna().sum()} "
            f"of {len(D)}; built by measured free pool {(D.rule == 'measured_free_pool').sum()}")
    tot = inv.mass_pg.sum(); inc = (inv.copies_in_complexes * inv.mass_da / 6.022e23 * 1e12).sum()
    ctx.say(f"complexes: parsed {len(C)}, built {(C.copies_built > 0).sum()}; protein mass in complexes {inc:.1f} of {tot:.1f} pg")


# ------------------------------------------------------------------------------------------------ SEC-MS fractions
def secms_fractions(st):
    """views of the SEC-MS toolkit tables: sec_assembly.tsv (assembled fraction at threshold_k) and the calibrant peaks"""
    ctx, view = st.ctx, st.view
    S = view.legacy_sec_assembly(int(st.o["secms_fractions"]["threshold_k"]))
    ctx.write_tsv(view.calibrants()[["cpx", "name", "n", "mw_kda", "peak", "co_height", "mcorr", "key", "kept"]],
                  "assembly/sec_calibrants_peak.tsv")
    ctx.write_tsv(S, "assembly/sec_assembly.tsv")
    ctx.say(f"SEC ({view.source_id}, {view.replicate}): proteins scored {len(S)}; median assembled fraction "
            f"{S.assembled_frac.median():.2f}")


# ------------------------------------------------------------------------------------------------ PDB homo-oligomers
_A = None
_MIN_RES = 30


def _n_chains(pid):
    import gemmi
    try:
        st = gemmi.read_structure(f"{_A}/{pid}-assembly1.cif.gz"); st.setup_entities()
        return sum(1 for ch in st[0] if ch.get_polymer().check_polymer_type() in (gemmi.PolymerType.PeptideL,)
                   and len(ch.get_polymer()) >= _MIN_RES)
    except Exception:
        return 0


def _init(a, min_res):
    global _A, _MIN_RES
    _A, _MIN_RES = a, min_res


def apex_check_exclude(st) -> dict:
    """data.apex_check_exclude (optional): proteins whose SEC-MS apex is a known artefact. List of {uniprot, gene,
    reason, share_source_id (optional)}; a listed protein skips the SEC apex check of pdb_homo_oligomers (it keeps its
    PDB state) and, when share_source_id is set, its oligomer share is a_k of that SEC-MS run (source_id; condition /
    replicate: the only ones, or 'mean') instead of the run that feeds the rules. Returns {uniprot: {"view": SecView | None}}.
    A reason is required (the list is curation, like the other data tables)"""
    rows = st.data.get("apex_check_exclude") or []
    if isinstance(rows, str):
        rows = read_curated(rows).to_dict("records")
    out = {}
    for r in rows:
        acc, why = str(r.get("uniprot") or "").strip(), str(r.get("reason") or "").strip()
        if not acc or not why:
            raise ValueError(f"assembly_steps.data.apex_check_exclude: every row needs uniprot and reason (got {r})")
        src = r.get("share_source_id") or None
        v = None
        if src:
            if st.view is None:
                raise ValueError(f"apex_check_exclude {acc}: share_source_id {src} needs the SEC-MS step (adapters.secms)")
            v = S2.SecView(st.ctx, {"source_id": str(src)})
            if not v.has(acc):
                raise ValueError(f"apex_check_exclude {acc}: no SEC-MS row in source {src}")
        out[acc] = {"view": v, "gene": r.get("gene"), "reason": why}
    return out


def pdb_homo_oligomers(st):
    """proteins mostly outside complexes whose single-UniProt PDB assemblies agree on a chain count >= 2: a SEC-MS
    share of the free pool (after complexes) becomes the oligomer"""
    ctx, h, view = st.ctx, st.o["pdb_homo_oligomers"], st.view
    lo = legacy_options(st.o)
    A, local = st.pdb_local()
    inv = pd.read_csv(ctx.need("assembly/inventory_v0_complexed.tsv"), sep="\t", index_col="gene")
    sifts = pd.read_csv(st.frozen("sifts"), sep="\t", comment="#", usecols=["PDB", "CHAIN", "SP_PRIMARY", "SP_BEG", "SP_END"])
    sifts = sifts[sifts.PDB.isin(local)]
    ent = sifts.groupby("PDB").SP_PRIMARY.agg(lambda s: set(s))
    homo = ent[ent.map(len) == 1].map(lambda s: next(iter(s)))
    L = pd.read_csv(st.frozen("uniprot_reviewed"), sep="\t", usecols=["Entry", "Length"]).set_index("Entry").Length
    sifts["span"] = pd.to_numeric(sifts.SP_END, errors="coerce") - pd.to_numeric(sifts.SP_BEG, errors="coerce") + 1
    cov = sifts.groupby(["PDB", "CHAIN", "SP_PRIMARY"]).span.sum().groupby(level=[0, 2]).max().reset_index()
    cov["frac"] = cov.span / cov.SP_PRIMARY.map(L)
    good = set(cov.loc[cov.frac >= h["min_chain_coverage"], "PDB"])
    homo = homo[homo.index.isin(good)]
    nmax = int(h["max_entries"])
    by_acc = homo.groupby(homo).apply(lambda s: sorted(s.index)[:nmax]).to_dict()
    sel = (inv.copies > 0) & (inv.frac_in_complexes.fillna(0) < h["max_frac_in_complexes"]) & inv.uniprot.isin(by_acc.keys())
    if st.data.get("homo_oligomer_exclude_regex"):
        sel &= ~inv.index.str.match(st.data["homo_oligomer_exclude_regex"])
    cand = inv[sel]
    try:
        sa = pd.read_csv(ctx.need("assembly/sec_assembly.tsv"), sep="\t").set_index("uniprot")
    except FileNotFoundError:
        sa = pd.DataFrame(columns=["apex_app_kda", "mono_kda"])
    excl = apex_check_exclude(st)
    jobs = {g: by_acc[a] for g, a in zip(cand.index, cand.uniprot)}
    allp = sorted({pp for v in jobs.values() for pp in v})
    with ProcessPoolExecutor(st.workers, initializer=_init, initargs=(A, int(h["min_chain_residues"]))) as ex:
        nc = dict(zip(allp, ex.map(_n_chains, allp, chunksize=50)))
    rows = []
    for g, pids in jobs.items():
        counts = [nc[pp] for pp in pids if nc[pp] >= 1]
        if not counts:
            continue
        state, support = collections.Counter(counts).most_common(1)[0]
        if state < 2 or support / len(counts) < h["min_support_frac"]:
            continue
        if state > h["large_state"] and support < h["large_state_min_support"]:
            continue
        a = inv.uniprot[g]
        a4 = view.a(a, h["share"]["apex_check_min_k"]) if view is not None else np.nan
        ex = excl.get(a)
        if ex is None and not S2.homo_apex_ok(a, state, sa, a4, lo):
            continue
        if ex is not None and ex["view"] is not None:          # share read from the named SEC-MS run instead
            a4 = ex["view"].a(a, h["share"]["apex_check_min_k"])
            share, _ = S2.homo_share(ex["view"], a, state, lo)
            ctx.say(f"homo-oligomers: {g} ({a}) apex check skipped (apex_check_exclude); share {share:.4f} from "
                    f"{ex['view'].source_id} / {ex['view'].condition} / {ex['view'].replicate}")
        else:
            share, _ = S2.homo_share(view, a, state, lo)
            if ex is not None:
                ctx.say(f"homo-oligomers: {g} ({a}) apex check skipped (apex_check_exclude); share {share:.4f}")
        free = inv.copies[g] - inv.copies_in_complexes.fillna(0)[g]
        n_olig = int(free * share // state)
        rows.append(dict(gene=g, uniprot=a, state=state, support=f"{support}/{len(pids)}", example_pdb=pids[0], sec_a4=a4,
                         share=share, n_oligomers=n_olig, copies_used=n_olig * state,
                         mass_pg=n_olig * state * inv.mass_da[g] / 6.022e23 * 1e12))
    H = pd.DataFrame(rows, columns=["gene", "uniprot", "state", "support", "example_pdb", "sec_a4", "share", "n_oligomers",
                                    "copies_used", "mass_pg"]).sort_values("mass_pg", ascending=False)
    ctx.write_tsv(H, "assembly/homo_oligomers_v0.tsv")
    ctx.say(f"homo-oligomers: candidates {len(cand)}, assigned {len(H)}")


# ------------------------------------------------------------------------------------------------ AFDB dimers
def _free_pool(st):
    """free copies after complexes and PDB homo-oligomers (when that step ran); shared by the AFDB dimer steps"""
    if st.free is None:
        ctx = st.ctx
        inv = pd.read_csv(ctx.need("assembly/inventory_v0_complexed.tsv"), sep="\t", index_col="gene")
        p = ctx.path("assembly/homo_oligomers_v0.tsv")
        ho = pd.read_csv(p, sep="\t").set_index("gene") if os.path.exists(p) else \
            pd.DataFrame(columns=["copies_used"], index=pd.Index([], name="gene"))
        used = inv.copies_in_complexes.fillna(0) + ho.copies_used.reindex(inv.index).fillna(0)
        st.free = (inv.copies - used).clip(lower=0)
        st.dimer_inv, st.homo = inv, ho
    return st.free


def afdb_homodimers(st):
    """AFDB homodimer models with max ipSAE > min_ipsae: a SEC-MS share of the free pool becomes dimers (best model per
    protein; genes with a PDB homo-oligomer skipped when skip_pdb_homo_oligomers)"""
    ctx, hd = st.ctx, st.o["afdb_homodimers"]
    lo = legacy_options(st.o)
    free = _free_pool(st)
    inv = st.dimer_inv
    MDA = inv.mass_da.fillna(0)
    acc2g = dict(zip(inv.uniprot, inv.index))
    skip = set(st.homo.index) if hd["skip_pdb_homo_oligomers"] else set()
    H = pd.read_csv(st.frozen("afdb_homodimers"), sep="\t")
    H["ipsae"] = H[["ipSAE_AB", "ipSAE_BA"]].max(axis=1)
    H = H[H.ipsae > hd["min_ipsae"]].sort_values("ipsae", ascending=False).drop_duplicates("uniprotAccession")
    rows = []
    for r in H.itertuples():
        g = acc2g.get(r.uniprotAccession)
        if g is None or g in skip or free[g] <= 0:
            continue
        share = S2.homodimer_share(st.view, r.uniprotAccession, lo)
        n = int(free[g] * share // 2)
        if n <= 0:
            continue
        free[g] -= 2 * n
        rows.append(dict(gene=g, uniprot=r.uniprotAccession, model=r.modelEntityId, ipsae=round(r.ipsae, 3),
                         share=round(share, 2), n_dimers=n, mass_pg=2 * n * MDA[g] / 6.022e23 * 1e12))
    HD = pd.DataFrame(rows, columns=["gene", "uniprot", "model", "ipsae", "share", "n_dimers", "mass_pg"]) \
        .sort_values("mass_pg", ascending=False)
    ctx.write_tsv(HD, "assembly/afdb_homodimers_v0.tsv")
    ctx.say(f"AFDB homodimers: {len(HD)} proteins, {HD.mass_pg.sum():.2f} pg")


def afdb_heterodimers(st):
    """AFDB heterodimer models passing the AFDB quality flag, supported by SEC-MS co-elution (r >= sec_r_min) or STRING
    (>= string_min), not already together in a built complex; greedy by ipSAE, n = bound_fraction x min(free a, free b)"""
    ctx, hx = st.ctx, st.o["afdb_heterodimers"]
    free = _free_pool(st)
    inv = st.dimer_inv
    MDA = inv.mass_da.fillna(0)
    acc2g = dict(zip(inv.uniprot, inv.index))
    X = pd.read_csv(st.frozen("afdb_heterodimers"), sep="\t")
    if hx["quality_column"]:
        X = X[X[hx["quality_column"]].astype(str).str.lower() == "true"]
    X = X[X.uniprot_ac_1.isin(acc2g) & X.uniprot_ac_2.isin(acc2g)]
    sec_r = st.view.pair_r() if st.view is not None else (lambda a, b: np.nan)
    ev = hx["evidence"]
    st_min = ev.get("string_min")
    st_score = {}
    acc2ensp = {}
    if st_min is not None:
        up = pd.read_csv(st.frozen("uniprot_full"), sep="\t", usecols=["Entry", "STRING"]).dropna()
        up["ensp"] = up.STRING.str.split(";").str[0]
        acc2ensp = dict(zip(up.Entry, up.ensp))
        need = {acc2ensp.get(a) for a in pd.concat([X.uniprot_ac_1, X.uniprot_ac_2]) if a in acc2ensp}
        with gzip.open(st.frozen("string_links"), "rt") as f:
            next(f)
            for line in f:
                p1, p2, sc = line.split()
                if p1 in need and p2 in need:
                    st_score[(p1, p2)] = int(sc)
    X["sec_r"] = [sec_r(a, b) for a, b in zip(X.uniprot_ac_1, X.uniprot_ac_2)]
    X["string"] = [st_score.get((acc2ensp.get(a), acc2ensp.get(b)), 0) for a, b in zip(X.uniprot_ac_1, X.uniprot_ac_2)]
    keep = X.sec_r >= ev["sec_r_min"] if ev.get("sec_r_min") is not None else pd.Series(False, X.index)
    if st_min is not None:
        keep = keep | (X.string >= st_min)
    X = X[keep].sort_values("max_ipSAE", ascending=False)
    if hx["skip_same_complex"]:
        C = pd.read_csv(ctx.need("assembly/complexes_v0.tsv"), sep="\t"); C = C[C.copies_built > 0]
        same = set()
        for sub in C.subunits.dropna():
            gs = [x.split(":")[0] for x in sub.split(";")]
            same |= {(a, b) for a in gs for b in gs if a != b}
        X = X[[(acc2g[a], acc2g[b]) not in same for a, b in zip(X.uniprot_ac_1, X.uniprot_ac_2)]]
    rows = []
    for r in X.itertuples():
        g1, g2 = acc2g[r.uniprot_ac_1], acc2g[r.uniprot_ac_2]
        n = int(hx["bound_fraction"] * min(free[g1], free[g2]))
        if n <= 0:
            continue
        free[g1] -= n; free[g2] -= n
        rows.append(dict(gene1=g1, gene2=g2, model=r.modelEntityId, ipsae=round(r.max_ipSAE, 3),
                         sec_r=round(r.sec_r, 2) if r.sec_r == r.sec_r else None, string=r.string, n=n,
                         mass_pg=n * MDA[g1] / 6.022e23 * 1e12 + n * MDA[g2] / 6.022e23 * 1e12))
    XD = pd.DataFrame(rows, columns=["gene1", "gene2", "model", "ipsae", "sec_r", "string", "n", "mass_pg"]) \
        .sort_values("mass_pg", ascending=False)
    ctx.write_tsv(XD, "assembly/afdb_heterodimers_v0.tsv")
    ctx.say(f"AFDB heterodimers: {len(XD)} pairs (of {len(X)} supported)")


STEPS = {"curated_machines": curated_machines, "complex_portal": complex_portal, "secms_fractions": secms_fractions,
         "pdb_homo_oligomers": pdb_homo_oligomers, "afdb_homodimers": afdb_homodimers,
         "afdb_heterodimers": afdb_heterodimers}
