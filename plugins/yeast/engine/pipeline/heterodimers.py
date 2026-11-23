# origin: yeast yeast_cell_dev_bundle/pipeline/heterodimers.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Inventory v3 item A3: predicted AlphaFold DB heterodimers (YC_DIMERS/heterodimers) as packing ingredients.

A pair is used when
  - both partners are left in the monomer pool (>= MIN_COPIES each) after complexes were assigned,
  - ipTM >= MIN_IPTM (the dataset is already selected for confident interfaces, ipSAE >= 0.6),
  - the partners share a compartment (or a membrane protein faces its soluble partner's compartment),
  - neither partner is a homo-oligomer of more than 2 chains unless Complex Portal lists the pair together
    (TDH1-TDH3 would be a fragment of mixed GAPDH tetramers),
  - the partners are not both subunits of a complex already modeled (their leftovers are MS surplus),
  - they are not two Hsp70s (HSP70, 2026-11-16), nor paralogues without interaction evidence (2026-11-16).
Copies: each protein's pool is split evenly among its accepted partners; pairs supported by Complex Portal or IntAct
bind fully, others bind BOUND_FRACTION of what is available. IntAct evidence is recorded but not required.
"""
import re, json
import numpy as np, pandas as pd, gemmi
import build_structs as bs
import membrane_orient as mo
from config import DATA, DIMERS

MIN_IPTM = 0.5
MIN_COPIES = 20
BOUND_FRACTION = 0.5
# Hsp70 chaperones act as monomers. The AlphaFold pair Ssa1-Ssa2 was accepted as 294,000 heterodimers (1.9% of protein
# mass, audit of build W, 2026-11-16) because both belong to the seven-protein HMC complex of Complex Portal; pairs within
# the family are not built. (A general rule on paralogues was tried and dropped: name-prefix and Alliance paralogue flags
# also mark real heteromers such as NAC, Srp14-Srp21, Tps1-Tps2, Phb1-Phb2.)
HSP70 = {"SSA1", "SSA2", "SSA3", "SSA4", "SSB1", "SSB2", "SSE1", "SSE2"}
# Paralogues share their homodimer interface, so AlphaFold predicts a confident heterodimer for many pairs that are
# homodimers in the cell (Gpd1-Gpd2, Ald4-Ald5, Aro3-Aro4, Idp2-Idp3, the Pho3/5/11/12 acid phosphatases). A paralogue
# pair is built only with interaction evidence: Complex Portal, an IntAct publication or a direct interaction
# (2026-11-16). Paralogue: an Alliance paralogue pair, or sequences sharing at least PARALOG_K3 of the shorter one's
# distinct 3-mers (accepted pairs of build W: median 0.11; Gpd1-Gpd2 0.49, Ald4-Ald5 0.37, Aro3-Aro4 0.34). A shared
# gene-name prefix is not used (Spt3-Spt15 and Lst4-Lst7 interact). Paralogue heteromers with evidence (NAC,
# Srp14-Srp21, Phb1-Phb2, Tps1-Tps2, Bmh1-Bmh2, the Lsm rings) are kept.
PARALOG_K3 = 0.30
PARALOG_NEEDS_EVIDENCE = True
SEQ = {}
FACING = {"ER", "plasma_membrane", "golgi_endosome", "vacuole", "nuclear_envelope", "mito_OM", "peroxisome", "lipid_droplet"}


def _k3(a, b):
    """share of the shorter sequence's distinct 3-mers found in the other"""
    if len(a) < 3 or len(b) < 3:
        return 0.0
    A = {a[i:i + 3] for i in range(len(a) - 2)}; B = {b[i:i + 3] for i in range(len(b) - 2)}
    return len(A & B) / min(len(A), len(B))


def _family(a, b):
    ma, mb = re.match(r"^([A-Z]+)(\d+)$", a), re.match(r"^([A-Z]+)(\d+)$", b)
    return bool(ma and mb and ma.group(1) == mb.group(1))


def compatible(wa, wb, tma, tmb):
    """compartment weights of the pair, or None if the partners cannot meet"""
    shared = {c: min(wa.get(c, 0), wb.get(c, 0)) for c in set(wa) & set(wb)}
    if sum(shared.values()) >= 0.3:
        z = sum(shared.values()); return {c: v / z for c, v in shared.items() if v > 0}
    for wt, ws, tm_t, tm_s in ((wa, wb, tma, tmb), (wb, wa, tmb, tma)):
        if tm_t and not tm_s:      # membrane protein + soluble partner on the facing side
            mem = {c: v for c, v in wt.items() if c in FACING or c == "mito_IM"}
            side = ws.get("cytosol", 0) if "mito_IM" not in mem else ws.get("mito_matrix", 0) + ws.get("mito_IM", 0)
            if mem and side >= 0.3:
                z = sum(mem.values()); return {c: v / z for c, v in mem.items()}
    return None


def build(lib_id, model_id, acc_a, acc_b, tm):
    st = gemmi.read_structure(f"{DIMERS}/heterodimers/{model_id}.pdb")
    st.remove_hydrogens()
    chains = [ch.name for ch in st[0]]
    xyz, kind, accs, pos = mo.beads_with_positions(st, acc_of_chain={chains[0]: acc_a, chains[1]: acc_b})
    pl = np.array([a.b_iso for ch in st[0] for res in ch for a in [res.find_atom("CA", "*")] if a is not None])
    tmm, ns, ps = mo.masks(accs, pos, xyz)
    keep = (pl >= 50) | tmm
    if keep.sum() < 20:
        keep[:] = True
    extra = dict(src="AF-multimer", model=model_id, nchain=2, plddt=pl[keep].astype(np.float32), plddt_frac_kept=np.float32(keep.mean()))
    info = {}
    if tm:
        out, info = mo.orient(xyz[keep], tmm[keep], ns[keep], ps[keep], min_score=mo.min_score_for(accs))
        if out is None:
            return None, info
        extra.update(membrane="TM", tm=tmm[keep], memb_score=np.float32(info["score"]))
    else:
        out = xyz[keep] - xyz[keep].mean(0)
    bs.save(lib_id, out, np.zeros(len(out), np.int8), extra=extra)
    return len(out), info


# histone-fold proteins pair promiscuously in AlphaFold (H4-H2B, H2A.Z-TAF9 ...; review comment 6): a pair with a
# histone is accepted only if it is a native histone pair
HISTONES = {"HTA1": "H2A", "HTA2": "H2A", "HTB1": "H2B", "HTB2": "H2B", "HHT1": "H3", "HHT2": "H3", "HHF1": "H4",
            "HHF2": "H4", "HTZ1": "H2A.Z", "CSE4": "CenH3"}
HISTONE_PAIRS = {frozenset(("H2A", "H2B")), frozenset(("H3", "H4")), frozenset(("H2A.Z", "H2B")), frozenset(("CenH3", "H4"))}


def assign(ing, df, modeled_complexes, legacy_oligo, log=print):
    """ing: inventory after complexes. modeled_complexes: list of gene sets of modeled complexes.
    Returns (ing with heterodimer ingredients added and partners consumed, decision table)."""
    h = pd.read_csv(f"{DIMERS}/heterodimers/candidates.csv")
    ia = pd.read_csv(f"{DATA}/intact_pairs.tsv", sep="\t")
    IA = {frozenset((a, b)): (p, d) for a, b, p, d in zip(ia.acc_a, ia.acc_b, ia.n_publications, ia.direct)}
    global SEQ
    if PARALOG_NEEDS_EVIDENCE and not SEQ:
        from features import read_fasta
        SEQ = read_fasta()
    par = pd.read_csv(f"{DATA}/alliance_yeast_paralog_pairs.csv")
    PAR = {frozenset((a, b)) for a, b in zip(par.gene1_uniprot, par.gene2_uniprot)}
    cp = pd.read_csv(f"{DATA}/complexportal_559292.tsv", sep="\t")
    CP = set()
    for s in cp.iloc[:, 4]:
        accs = re.findall(r"([A-Z0-9]{6,10})(?:-\w+)?\(\d+\)", s)
        CP |= {frozenset((a, b)) for a in accs for b in accs if a != b}
    olig = dict(zip(df.uniprot, df.oligomer_n)); gene_of = dict(zip(df.uniprot, df.gene))
    prot = ing[ing.kind == "protein"].set_index("uniprot")
    in_complex = {}
    for gs in modeled_complexes:
        for g in gs:
            in_complex.setdefault(g, []).append(frozenset(gs))
    legacy_genes = set(legacy_oligo)

    rows, acc_ok = [], []
    for r in h.itertuples():
        a, b = r.uniprotAccession_1, r.uniprotAccession_2; k = frozenset((a, b))
        ga, gb = gene_of.get(a, r.gene_1), gene_of.get(b, r.gene_2)
        why = None
        ca = prot.copies_per_cell.get(a, 0) if a in prot.index else 0
        cb = prot.copies_per_cell.get(b, 0) if b in prot.index else 0
        pubs, direct = IA.get(k, (0, False))
        cpx = k in CP
        if min(ca, cb) < MIN_COPIES: why = "partner not in monomer pool"
        elif r.ipTM < MIN_IPTM: why = f"ipTM < {MIN_IPTM}"
        elif not cpx and (olig.get(a, 1) > 2 or olig.get(b, 1) > 2 or ga in legacy_genes or gb in legacy_genes): why = "partner is a larger homo-oligomer"
        elif set(in_complex.get(ga, [])) & set(in_complex.get(gb, [])): why = "both are subunits of a modeled complex"
        elif (ga in HISTONES or gb in HISTONES) and frozenset((HISTONES.get(ga, ga), HISTONES.get(gb, gb))) not in HISTONE_PAIRS:
            why = "non-native histone-fold pair"
        elif ga in HSP70 and gb in HSP70:
            why = "Hsp70 pair (Hsp70s act as monomers)"
        elif PARALOG_NEEDS_EVIDENCE and not (cpx or pubs > 0 or direct) and (k in PAR or _k3(SEQ.get(a, ""), SEQ.get(b, "")) >= PARALOG_K3):
            why = "paralogue pair without interaction evidence"
        w = None
        if why is None:
            w = compatible(json.loads(prot.loc[a, "comp_weights"]), json.loads(prot.loc[b, "comp_weights"]),
                           prot.n_tm.get(a, 0) > 0, prot.n_tm.get(b, 0) > 0)
            if w is None: why = "no shared compartment"
        rows.append(dict(model=r.modelEntityId, a=ga, b=gb, acc_a=a, acc_b=b, ipTM=r.ipTM, ipSAE=r.max_ipSAE, publications=pubs,
                         direct=direct, complex_portal=cpx, paralog=k in PAR or _family(ga, gb), decision=why or "accepted", copies=0.0))
        if why is None:
            acc_ok.append((len(rows) - 1, w))

    # allocation: each protein offers its copies to its accepted partners in proportion to the partners' copies, and a
    # pair forms the smaller of the two offers (2026-11-12; before, copies were split evenly over the partners, so a
    # rare paralog such as Btt1 took half of Egd2 and the obligate NAC heterodimer Egd1-Egd2 was half formed)
    partners = {}
    for i, _ in acc_ok:
        a_, b_ = rows[i]["acc_a"], rows[i]["acc_b"]
        partners.setdefault(a_, []).append(b_); partners.setdefault(b_, []).append(a_)
    cpc = prot.copies_per_cell
    offer = lambda x, y: cpc[x] * cpc[y] / sum(cpc[z] for z in partners[x])
    pool = prot.copies_per_cell.to_dict()
    new = []
    for i, w in acc_ok:
        r = rows[i]; a, b = r["acc_a"], r["acc_b"]
        frac = 1.0 if (r["complex_portal"] or r["publications"] > 0) else BOUND_FRACTION
        n = frac * min(offer(a, b), offer(b, a))
        tm = prot.n_tm.get(a, 0) > 0 or prot.n_tm.get(b, 0) > 0
        lib_id = re.sub(r"\W+", "_", f"{r['a']}-{r['b']}").replace("-", "_") + "_het"
        built, info = build(lib_id, r["model"], a, b, tm)
        if built is None:
            r["decision"] = f"membrane orientation failed ({info.get('fail')})"; continue
        r["copies"] = n
        pool[a] -= n; pool[b] -= n
        new.append(dict(ing_id=lib_id, label=f"{r['a']}-{r['b']} heterodimer", kind="complex", genes=f"{r['a']};{r['b']}",
                        n_subunit_types=2, copies_per_cell=n, comp_weights=json.dumps(w),
                        mass_da=float(df.set_index("uniprot").loc[[a, b], "mass_da"].sum()), pdb="", membrane="",
                        oligomer_n=1, uniprot="", n_tm=int(prot.n_tm.get(a, 0) + prot.n_tm.get(b, 0))))
    ing = ing.copy()
    m = ing.kind == "protein"
    ing.loc[m, "copies_per_cell"] = [max(0.0, pool.get(u, c)) for u, c in zip(ing.loc[m, "uniprot"], ing.loc[m, "copies_per_cell"])]
    ing = pd.concat([ing, pd.DataFrame(new)], ignore_index=True)
    ing = ing[(ing.kind != "protein") | (ing.copies_per_cell >= 1)].reset_index(drop=True)
    t = pd.DataFrame(rows)
    log(f"heterodimers: {len(new)} accepted; " + str(t.decision.str.replace(r" \(.*", "", regex=True).value_counts().to_dict()))
    return ing, t
