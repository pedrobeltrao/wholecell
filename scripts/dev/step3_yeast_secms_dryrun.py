"""Step 3 yeast SEC-MS assembly: inventory-only dry run and the criteria metrics
(docs/unification_criteria/step3_yeast_secms_assembly.md; frozen before this script).

    python scripts/step3_yeast_secms_dryrun.py --ref-build <complete yeast build> --ref-inputs <suite>/adapter/loc_inputs --out <new dir>

Two pseudo-builds from the same code and inputs, <out>/legacy (assembly preset legacy_x5) and <out>/secms (preset
x5_secms of config/cells/yeast.yaml): inventory/{yeast, v0, copies} copied from the reference build, engine/ = a private
copy of the reference engine (lib and pdb copied, root linked, work = links to the reference work entries except the
RCSB cache, which the structure hooks rewrite, and numba_cache), then the inventory steps secms -> assembly ->
ingredients -> pack_table. Nothing is written into the reference build. Then report.json: S1 calibration QC, S2
comp_weights, S3 §3d benchmark (b5_localisation_check.yeast_metrics on inventory/v0), S4 agreement with SEC, S5 mass
moved per rule (+ top 20 changes), S6 complexes_added identical; legacy run vs the reference build's assembly tables.
"""
import argparse
import json
import os
import shutil
import sys
import time

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
sys.path[:0] = [os.path.join(REPO, "scripts", "dev"), os.path.join(REPO, "scripts")]

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

AV = 6.02214076e23
TAGS = {"legacy": {"inventory": {"params": {"assembly_steps": {"preset": "legacy_x5", "overrides": None}}}}, "secms": None}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def make_pseudo(d, ref):
    os.makedirs(os.path.join(d, "inventory"), exist_ok=True)
    for x in ("yeast", "v0", "copies"):
        shutil.copytree(os.path.join(ref, "inventory", x), os.path.join(d, "inventory", x))
    os.makedirs(os.path.join(d, "geometry"), exist_ok=True)
    for f in ("measures.json", "package.json"):
        shutil.copyfile(os.path.join(ref, "geometry", f), os.path.join(d, "geometry", f))
    e, re_ = os.path.join(d, "engine"), os.path.join(ref, "engine")
    os.makedirs(os.path.join(e, "work", "numba_cache"))
    os.symlink(os.path.join(re_, "root"), os.path.join(e, "root"))
    shutil.copytree(os.path.join(re_, "lib"), os.path.join(e, "lib"), symlinks=True)
    shutil.copytree(os.path.join(re_, "pdb"), os.path.join(e, "pdb"), symlinks=True)
    for f in os.listdir(os.path.join(re_, "work")):
        if f == "numba_cache":
            continue
        src, dst = os.path.join(re_, "work", f), os.path.join(e, "work", f)
        if f == "rcsb_cache.json":
            shutil.copyfile(src, dst)
        else:
            os.symlink(src, dst)


def run(tag, d, ref):
    from wholecell import registry
    from wholecell.config import load
    from wholecell.inventory import stage as INV
    cfg = load("yeast", overrides=TAGS[tag])
    registry.load_plugins(cfg)
    make_pseudo(d, ref)
    for s in ("secms", "assembly", "ingredients", "pack_table"):
        t = time.time()
        INV.run(cfg, d, step=s)
        log(f"{tag}: inventory {s} {time.time() - t:.0f} s")


# ------------------------------------------------------------------------------------------------ metrics
def model_share(inv):
    """per protein: m_p = 1 - copies in monomer ingredients / protein-table copies (criteria §2)"""
    T = pd.read_parquet(os.path.join(inv, "yeast", "protein_table.parquet"))
    T = T[T.copies_per_cell > 0].drop_duplicates("uniprot").set_index("uniprot")
    I = pd.read_csv(os.path.join(inv, "assembly", "ingredients.csv"))
    om = json.load(open(os.path.join(inv, "assembly", "oligo_map.json")))
    P = I[I.kind == "protein"]
    mono = P[~P.ing_id.isin(list(om))].groupby("uniprot").copies_per_cell.sum()
    homo = P[P.ing_id.isin(list(om))].groupby("uniprot").copies_per_cell.sum()
    M = pd.DataFrame(dict(gene=T.gene, copies=T.copies_per_cell, mass_da=T.mass_da))
    M["mono"] = mono.reindex(M.index).fillna(0.0)
    M["homo"] = homo.reindex(M.index).fillna(0.0)
    M["m"] = (1 - M.mono / M.copies).clip(0, 1)
    M["mass_pg"] = M.copies * M.mass_da / AV * 1e12
    return M


def same_definition(inv, M, k):
    """model share computed as SEC's (addendum 2026-11-19): per protein, copies in an assembly whose mass >= k x the
    monomer mass / copies (homo-oligomer: state n >= k; heterodimer: pair mass; complex: the largest complex ingredient
    naming the gene), copies not in any ingredient (surplus not placed) left out of the denominator"""
    I = pd.read_csv(os.path.join(inv, "assembly", "ingredients.csv"))
    om = json.load(open(os.path.join(inv, "assembly", "oligo_map.json")))
    P = I[I.kind == "protein"]
    homo = P[P.ing_id.isin(list(om))].assign(n=lambda d: d.ing_id.map(lambda i: om[i][1]))
    homo_k = homo[homo.n >= k].groupby("uniprot").copies_per_cell.sum()
    H = pd.read_csv(os.path.join(inv, "assembly", "heterodimers.csv"))
    H = H[H.copies > 0]
    mda = M.mass_da
    pm = H.acc_a.map(mda).fillna(0) + H.acc_b.map(mda).fillna(0)
    het_all = pd.concat([H.groupby("acc_a").copies.sum(), H.groupby("acc_b").copies.sum()]).groupby(level=0).sum()
    ka = H[pm >= k * H.acc_a.map(mda)].groupby("acc_a").copies.sum()
    kb = H[pm >= k * H.acc_b.map(mda)].groupby("acc_b").copies.sum()
    het_k = pd.concat([ka, kb]).groupby(level=0).sum()
    C = I[(I.kind == "complex") & ~I.ing_id.str.endswith("_het")]
    cm = {}
    for g, m in zip(C.genes, C.mass_da):
        for x in str(g).split(";"):
            x = x.strip()
            cm[x] = max(cm.get(x, 0.0), float(m))
    cmax = M.gene.map(lambda g: max([cm.get(x.strip(), 0.0) for x in str(g).split(";")] or [0.0]))
    other = (M.copies - M.mono - M.homo - het_all.reindex(M.index).fillna(0)).clip(lower=0)
    in_cpx = other.where(cmax > 0, 0.0)                     # copies with no complex ingredient: surplus not placed
    num = homo_k.reindex(M.index).fillna(0) + het_k.reindex(M.index).fillna(0) + in_cpx.where(cmax >= k * mda, 0.0)
    den = M.copies - other.where(cmax == 0, 0.0)
    return (num / den).clip(0, 1), float((other.where(cmax == 0, 0.0) * mda).sum() / AV * 1e12)


def agreement(M, S):
    out = {"total_pg": float(M.mass_pg.sum()), "model_share_all": float((M.m * M.mass_pg).sum() / M.mass_pg.sum())}
    for k in (2, 4):
        a = S[S.threshold_k == k].set_index("uniprot").assembled_frac
        m = M[M.index.isin(a.index)].copy()
        m["a"] = a.reindex(m.index)
        m = m.dropna(subset=["a"])
        w = m.mass_pg
        out[f"k{k}"] = dict(n=int(len(m)), covered_share=float(w.sum() / M.mass_pg.sum()),
                            model=float((m.m * w).sum() / w.sum()), sec=float((m.a * w).sum() / w.sum()),
                            spearman=float(m[["m", "a"]].corr("spearman").iloc[0, 1]),
                            mae_w=float(((m.m - m.a).abs() * w).sum() / w.sum()))
    return out


def sec_view(inv):
    S = pd.read_csv(os.path.join(inv, "secms", "secms_protein.tsv"), sep="\t", low_memory=False,
                    dtype={"replicate": str, "condition": str})
    return S[S.source_id == "yeast_final"]


def calib_qc(inv):
    d = json.load(open(os.path.join(inv, "secms", "diagnostics.json")))
    rows = []
    for ds in d["datasets"]:
        for run_, r in ds["runs"].items():
            c = r["calibration"]
            rows.append(dict(source=ds["source_id"], run=run_, used=c["used"], n=c["n_kept"], rho=round(c["spearman_rho"], 3),
                             resid_fold=round(c["resid_fold"], 3), qc_fail=c["qc_fail"]))
    return dict(anchor_set=d.get("anchor_set"), runs=rows, all_pass=all(not r["qc_fail"] for r in rows))


def comp_weights_check(L, N):
    out = {}
    for f, key in (("assembly/ingredients.csv", "ing_id"), ("ingredients/ingredients.tsv", None), ("pack_table.tsv", None)):
        pa, pb = os.path.join(L, f), os.path.join(N, f)
        sep = "," if f.endswith(".csv") else "\t"
        A, B = pd.read_csv(pa, sep=sep), pd.read_csv(pb, sep=sep)
        k = key or next(c for c in ("ing_id", "species_id", "id") if c in A.columns)
        cw = next((c for c in ("comp_weights", "compartments") if c in A.columns), None)
        if cw is None:
            out[f] = "no comp_weights column"
            continue
        a, b = A.set_index(k)[cw].astype(str), B.set_index(k)[cw].astype(str)
        common = a.index.intersection(b.index)
        diff = [i for i in common if a[i] != b[i]]
        new = [i for i in b.index.difference(a.index)]
        bad_new = [i for i in new if "__homo" in str(i) and str(i).split("__homo")[0] in b.index
                   and b[i] != b[str(i).split("__homo")[0]]]
        other_new = [i for i in new if "__homo" not in str(i)]
        out[f] = dict(key=k, n_common=len(common), n_diff=len(diff), diff=diff[:10], n_new=len(new),
                      n_new_bad=len(bad_new), new_other=other_new[:10], n_dropped=len(a.index.difference(b.index)),
                      dropped=list(a.index.difference(b.index))[:10])
    return out


def mass_moved(L, N, ML, MN):
    total = float(ML.mass_pg.sum())
    mda = ML.mass_da
    # R1 / R2: homo-oligomer copies per protein (by the legacy decision: PDB vs AF dimer)
    oL = pd.read_csv(os.path.join(L, "assembly", "oligomers.csv")).set_index("uniprot")
    oN = pd.read_csv(os.path.join(N, "assembly", "oligomers.csv")).set_index("uniprot")
    dh = (MN.homo - ML.homo).reindex(ML.index).fillna(0)
    dec = oL.decision.reindex(ML.index).fillna(oN.decision.reindex(ML.index)).fillna("")
    r2 = dec.str.startswith("AF dimer")
    rows = []
    for rule, sel in (("R1 homo-oligomer share (PDB / audit)", ~r2), ("R2 AF homodimer share", r2)):
        pg = float((dh[sel].abs() * mda[sel]).sum() / AV * 1e12)
        rows.append(dict(rule=rule, n_proteins=int((dh[sel].abs() >= 1).sum()), mass_pg=pg, frac=pg / total))
    # R3: heterodimer copies per pair
    hL = pd.read_csv(os.path.join(L, "assembly", "heterodimers.csv"))
    hN = pd.read_csv(os.path.join(N, "assembly", "heterodimers.csv"))
    key = ["acc_a", "acc_b"]
    H = hL[key + ["a", "b", "copies"]].merge(hN[key + ["copies", "decision"]], on=key, how="outer", suffixes=("_L", "_N")).fillna({"copies_L": 0, "copies_N": 0})
    H["mass_pair"] = H.acc_a.map(mda).fillna(0) + H.acc_b.map(mda).fillna(0)
    H["d_pg"] = (H.copies_N - H.copies_L) * H.mass_pair / AV * 1e12
    pg = float(H.d_pg.abs().sum())
    rows.append(dict(rule="R3 heterodimer evidence (SEC r / STRING)", n_proteins=int((H.d_pg.abs() > 0).sum()), mass_pg=pg,
                     frac=pg / total, pairs_legacy=int((H.copies_L > 0).sum()), pairs_new=int((H.copies_N > 0).sum())))
    # total state change per protein (monomer <-> assembled)
    dm = ((MN.m - ML.m) * ML.mass_pg).reindex(ML.index).fillna(0)
    rows.append(dict(rule="total (monomer <-> assembled, net per protein)", n_proteins=int((dm.abs() > 1e-9).sum()),
                     mass_pg=float(dm.abs().sum()), frac=float(dm.abs().sum() / total),
                     to_assembled_pg=float(dm.clip(lower=0).sum()), to_free_pg=float(-dm.clip(upper=0).sum())))
    for r in rows:
        r["review"] = bool(r["frac"] > 0.02)
    top = pd.DataFrame(dict(gene=ML.gene, mass_pg=ML.mass_pg, m_legacy=ML.m, m_secms=MN.m.reindex(ML.index), d_pg=dm,
                            legacy_decision=dec, new_decision=oN.decision.reindex(ML.index)))
    top = top.reindex(top.d_pg.abs().sort_values(ascending=False).index).head(20)
    return rows, top, H


def bench(inv, ref_inputs, image, meas):
    from b5_localisation_check import cw_table, yeast_metrics
    m = yeast_metrics(cw_table(inv, "yeast"), ref_inputs, image, meas)
    return {k: m[k] for k in ("SWAT_native.intended.tophit.unweighted", "hyperLOPIT.intended.tophit.unweighted",
                              "image_nc_intended.test.spearman")}


def aggregate_detail(L, N, ML, MN, S, H):
    """aggregate decomposition (addendum 2026-11-19): change of the mass-weighted model share (k2-covered proteins, % points)
    by rule / cause; the rule set (proteins R1/R2/R3 touch) model vs SEC; SEC calibration band; same-definition shares"""
    A = S.pivot_table(index="uniprot", columns="threshold_k", values="assembled_frac")
    S2 = S[S.threshold_k == 2].drop_duplicates("uniprot").set_index("uniprot")
    idx = ML.index.intersection(A[2].dropna().index)
    w = ML.mass_pg[idx]
    W = w.sum()
    oL = pd.read_csv(os.path.join(L, "assembly", "oligomers.csv")).drop_duplicates("uniprot").set_index("uniprot").decision
    oN = pd.read_csv(os.path.join(N, "assembly", "oligomers.csv")).drop_duplicates("uniprot").set_index("uniprot").decision
    dN, dL = oN.reindex(idx).fillna(""), oL.reindex(idx).fillna("")
    het_d = set(H[(H.copies_N - H.copies_L).abs() > 0.5][["acc_a", "acc_b"]].values.ravel())
    cat = pd.Series("unchanged", index=idx)
    af = lambda d: d.str.match(r"AF dimer AF-")           # accepted AF dimers (not 'AF dimer rejected / ... failed')
    cat[af(dN) & ~af(dL)] = "R2 AF dimer (new: ipSAE-only gate)"
    cat[af(dN) & af(dL)] = "R2 AF dimer share"
    cat[dN.str.startswith("PDB")] = "R1 PDB share"
    cat[dN.str.contains("apex above state")] = "R1 apex above state (kept)"
    cat[dN.str.contains("apex check: monomer")] = "R1 apex-check demotion"
    cat[idx.isin(list(het_d)) & (cat == "unchanged")] = "R3 heterodimer gate"
    dm = (MN.m - ML.m).reindex(idx).fillna(0)
    by = {c: round(float((dm[cat == c] * w[cat == c]).sum() / W * 100), 3) for c in sorted(set(cat))}
    rs = cat != "unchanged"
    mN, mL = MN.m.reindex(idx), ML.m.reindex(idx)
    out = dict(dm_pts_by_rule=by, dm_pts_total=round(float((dm * w).sum() / W * 100), 3),
               rule_set=dict(n=int(rs.sum()), mass_share=float(w[rs].sum() / W), legacy=float((mL[rs] * w[rs]).sum() / w[rs].sum()),
                             new=float((mN[rs] * w[rs]).sum() / w[rs].sum()),
                             sec_k2=float((A[2].reindex(idx)[rs] * w[rs]).sum() / w[rs].sum())),
               sec_k2_calibration_band=[float((S2.assembled_frac_lo.reindex(idx) * w).sum() / W),
                                        float((S2.assembled_frac_hi.reindex(idx) * w).sum() / W)])
    mono_only = (MN.mono / MN.copies).reindex(idx) > 0.99
    out["structureless_monomers"] = dict(mass_share=float(w[mono_only].sum() / W),
                                         sec_k2=float((A[2].reindex(idx)[mono_only] * w[mono_only]).sum() / w[mono_only].sum()),
                                         under_gap_pts=float(((A[2].reindex(idx) - mN).clip(lower=0)[mono_only] * w[mono_only]).sum() / W * 100))
    for k in (2, 4):
        a = A[k].reindex(idx)
        sd = {}
        for nm, inv, M in (("legacy", L, ML), ("new", N, MN)):
            m, unplaced = same_definition(inv, M, k)
            m = m.reindex(idx)
            sd[nm] = float((m * w).sum() / W)
            sd["unplaced_pg"] = unplaced
        sd["sec"] = float((a * w).sum() / W)
        out[f"same_definition_k{k}"] = sd
    return out


def report(out, ref, ref_suite):
    L, N = os.path.join(out, "legacy", "inventory"), os.path.join(out, "secms", "inventory")
    ML, MN = model_share(L), model_share(N)
    S = sec_view(N)
    agL, agN = agreement(ML, S), agreement(MN, S)
    rows, top, H = mass_moved(L, N, ML, MN)
    agg = aggregate_detail(L, N, ML, MN, S, H)
    top.to_csv(os.path.join(out, "top20_changes.tsv"), sep="\t", index=False, float_format="%.4g")
    H.to_csv(os.path.join(out, "heterodimer_changes.tsv"), sep="\t", index=False, float_format="%.4g")
    pd.DataFrame(rows).to_csv(os.path.join(out, "mass_moved.tsv"), sep="\t", index=False, float_format="%.4g")
    cw = comp_weights_check(L, N)
    rR = os.path.join(ref, "inventory", "assembly")
    same_ref = {}
    for f in ("ingredients.csv", "oligo_map.json", "heterodimers.csv", "oligomers.csv", "complexes_added.csv"):
        a, b = os.path.join(rR, f), os.path.join(L, "assembly", f)
        same_ref[f] = open(a, "rb").read() == open(b, "rb").read()
    cpx_same = open(os.path.join(L, "assembly", "complexes_added.csv"), "rb").read() == \
        open(os.path.join(N, "assembly", "complexes_added.csv"), "rb").read()
    v0_same = open(os.path.join(L, "v0", "inventory_v0.tsv"), "rb").read() == open(os.path.join(N, "v0", "inventory_v0.tsv"), "rb").read()
    ref_inputs = ref_suite
    image = os.path.join(ref, "engine", "work", "validation", "image_nc.tsv")
    meas = json.load(open(os.path.join(ref, "geometry", "measures.json")))
    bL, bN = bench(L, ref_inputs, image, meas), bench(N, ref_inputs, image, meas)
    R = dict(reference_build=ref, calibration=calib_qc(N), agreement=dict(legacy=agL, secms=agN), mass_moved=rows,
             comp_weights=cw, v0_identical=v0_same, complexes_added_identical=cpx_same, legacy_vs_reference_build=same_ref,
             bench=dict(legacy=bL, secms=bN), top20=top.round(4).to_dict("records"), aggregate=agg)
    s = R["calibration"]["all_pass"]
    cw_ok = v0_same and all(isinstance(v, dict) and v["n_diff"] == 0 and v["n_new_bad"] == 0 for v in cw.values())
    b_ok = all(abs(bL[k] - bN[k]) < 1e-12 for k in bL) and bN["SWAT_native.intended.tophit.unweighted"] >= 0.655 and \
        bN["hyperLOPIT.intended.tophit.unweighted"] >= 0.838 and bN["image_nc_intended.test.spearman"] >= 0.80
    a_ok = agN["k2"]["spearman"] > agL["k2"]["spearman"] and agN["k2"]["mae_w"] < agL["k2"]["mae_w"]
    review = [r["rule"] for r in rows if r["review"]]
    R["checks"] = dict(S1=s, S2=cw_ok, S3=b_ok, S4=a_ok, S5_review=review, S6=cpx_same)
    R["verdict"] = ("FAIL" if not (s and cw_ok and b_ok and a_ok and cpx_same) else "NEEDS REVIEW" if review else "PASS")
    json.dump(R, open(os.path.join(out, "report.json"), "w"), indent=1, default=str)
    print(json.dumps({k: R[k] for k in ("checks", "verdict")}, indent=1))
    print(json.dumps(R["agreement"], indent=1)); print(json.dumps(agg, indent=1)); print(pd.DataFrame(rows).to_string())
    print(top.round(3).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref-build", required=True)
    ap.add_argument("--ref-inputs", required=True, help="<reference suite run>/adapter/loc_inputs (reference sets, masses; "
                    "as b5_localisation_check.py --ref-inputs)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--report-only", action="store_true")
    ap.add_argument("--tag", choices=sorted(TAGS), help="run one pseudo-build in this process (internal)")
    a = ap.parse_args()
    if a.tag:                                      # the yeast engine is configured once per process (env YC_*)
        run(a.tag, os.path.join(a.out, a.tag), a.ref_build)
        return
    if not a.report_only:
        import subprocess
        for tag in TAGS:
            d = os.path.join(a.out, tag)
            if os.path.exists(d):
                sys.exit(f"{d} exists (write-once)")
            subprocess.run([sys.executable, os.path.abspath(__file__), "--ref-build", a.ref_build, "--out", a.out,
                            "--ref-inputs", a.ref_inputs, "--tag", tag], check=True)
    report(a.out, a.ref_build, a.ref_inputs)


if __name__ == "__main__":
    main()
