"""Chain rule of the library models built from an experimental assembly or an AFDB complex model (library rule
2026-11-25, Pedro Beltrao: "the monomer models should have only the chain of the protein of interest"; clarified the same
day: the chains of an assembly are its subunits and are not simply taken out - only chains positively identified as
other proteins are removed, ambiguous ones are kept and flagged for review).

Every chain of the structure is put in one of three classes:
  1 member   (kept)     its UniProt accession (SIFTS for PDB entries, per residue range of the cached RCSB assembly 1;
                        the model's own _struct_ref for AFDB complex models) is one of the species' accessions, or a
                        declared one (allow), or maps (UniProt primary / synonym gene names) to a gene of the species'
                        'genes' field (gene_match); a declared chain; a declared nucleic acid. All copies are kept
                        (homo-oligomer stoichiometry intact) - except for a protein species (kind protein, a monomer
                        model): one chain of the protein, the one with most residues.
  2 foreign  (removed)  positively identified as another protein: an antibody / Fab / nanobody / TCR-like chain (entity
                        description, foreign_desc regex; unless it is a member) or a chain mapped by SIFTS to a human
                        accession that is not a member, not a paralogue of a member (paralogs table) in a species whose
                        accession list is complete (partners, ligands, other complexes of the crystal). Residues of a
                        member chain mapped to such an accession (fusion partners: T4 lysozyme, BRIL) are removed too.
                        An unmapped or non-human chain whose entity description is the reviewed UniProt protein name /
                        primary gene symbol of human proteins (genes.columns.names; uL2-style ribosomal names expanded)
                        is judged by those accessions: a member's name -> member (ortholog), another protein -> foreign.
                        Nucleic acids (classify_nucleic): declared for the species (nucleic_acids: nucleosome DNA,
                        7SL / snRNA / snoRNA / telomerase / RNase P RNA, rRNA of ribosomal species) -> member; rRNA of a
                        species that does not declare it, and tRNA / mRNA in a structure holding rRNA (ribosome
                        context) -> foreign; any other undeclared nucleic acid -> flagged.
  3 flagged  (kept)     ambiguous: no SIFTS mapping (incl. short peptides such as an MHC-bound peptide), a nucleic acid
                        not declared, a non-human chain not identified by name (possible ortholog; SIFTS taxonomy), a
                        paralogue of a member, a
                        non-member chain in a species with an empty entry in its accession list (the missing subunit may
                        be it). Listed per species in the ledger for review; a declaration (allow / nucleic_acids) turns
                        it into a member, removing it needs a rule change.
Each model records the decision in its npz: `chain_set` (JSON: kept chains with class and reason, removed chains, declared
subunits without a chain), `n_struct_beads` (beads from the structure after the rule) and `beads_added` (AF grafts of
missing TM segments, built TM helices). The library stage then runs `check` on every structure-derived model and FAILS
when a kept chain is class 2 (classes recomputed from SIFTS / the AFDB model, independent of the builder's per-residue
mapping), when a model carries no record, or when it holds more beads than its kept chains plus the added ones.

Config (library.chain_filter):
    chain_filter: {enabled: true, gene_match: true,
                   genes: {path: <UniProt table>, columns: {acc: Entry, genes: [Gene Names (primary), Gene Names (synonym)],
                                                            names: Protein names, reviewed: Reviewed}},
                   nucleic_acids: <TSV ing_id, desc_regex, reason>, allow: <TSV ing_id, what (accession|chain), value, reason>,
                   sifts: <SIFTS pdb_chain_uniprot>, taxonomy: <SIFTS pdb_chain_taxonomy>,
                   paralogs: {path: <Ensembl paralogues TSV>, columns: {gene: Gene name, paralog: Human paralogue associated gene name}},
                   foreign_desc: <regex>, not_foreign_desc: <regex>}
"""
from __future__ import annotations

import json
import os
import re

import numpy as np

KEYS = ("enabled", "gene_match", "genes", "nucleic_acids", "allow", "sifts", "taxonomy", "paralogs", "foreign_desc",
        "not_foreign_desc")
STRUCT_PREFIXES = ("pdb:", "afdb_complex:", "fallback:pdb_chain:")
HUMAN = "9606"
# entity descriptions of chains that are not subunits of the species (antibody fragments, binders, TCRs)
FOREIGN_DESC = (r"antibod|\bfab\b|fab\s|fab'|\bf\(ab|nanobod|\bvhh\b|sybod|megabod|\bscfv\b|darpin|monobod|affibod|"
                r"immunoglobulin|\bigg\d?\b|\bigm\b|kappa|lambda light|light chain|heavy chain|\bmab\b|t[- ]?cell receptor|"
                r"\btcr\b|designed binder|binder protein")
NOT_FOREIGN_DESC = r"myosin|dynein|clathrin|ferritin|kinesin|neurofilament|spectrin|tubulin|immunoglobulin-like|ig-like"


class Policy:
    def __init__(self, gene_match=False, acc_genes=None, nucleic=None, allow_accs=None, allow_chains=None, paralogs=None,
                 acc_primary=None, foreign_desc=FOREIGN_DESC, not_foreign_desc=NOT_FOREIGN_DESC, gene_acc=None, name_accs=None):
        self.gene_match = bool(gene_match)
        self.acc_genes = acc_genes or {}             # accession -> frozenset of gene names (primary + synonyms)
        self.acc_primary = acc_primary or {}         # accession -> primary gene name
        self.gene_acc = gene_acc or {}               # primary gene name -> accession
        self.nucleic = nucleic or {}                 # ing_id -> regex on the entity description
        self.allow_accs = allow_accs or {}           # ing_id -> extra accessions (members)
        self.allow_chains = allow_chains or {}       # ing_id -> chain ids (members)
        self.paralogs = paralogs or {}               # gene -> set of paralogue genes
        self.name_accs = name_accs or {}             # normalised UniProt protein name / gene symbol -> human accessions
        self.foreign_rx = re.compile(foreign_desc, re.I) if foreign_desc else None
        self.taxonomy_path = None                    # SIFTS pdb_chain_taxonomy (chain_filter.taxonomy)
        self.taxonomy = {}                           # {pdb: {chain: set of taxonomy ids}} of the entries of a step
        self.sifts_path = None
        self.not_foreign_rx = re.compile(not_foreign_desc, re.I) if not_foreign_desc else None

    def nucleic_ok(self, ing, desc):
        rx = self.nucleic.get(ing)
        return rx is not None and re.search(rx, desc or "", re.I) is not None

    def accs_by_name(self, desc):
        """human accessions whose reviewed UniProt protein name or primary gene symbol is the entity description (an
        unmapped or non-human chain identified by name; ribosomal nomenclature uL2 / eS6 expanded)"""
        n = norm_name(desc)
        if not n:
            return set()
        m = RIBO_RX.match(n)
        if m:
            n = f"{'large' if m.group(2) == 'l' else 'small'} ribosomal subunit protein {m.group(1)}{m.group(2)}{m.group(3)}"
        return set(self.name_accs.get(n, ()))

    def foreign_by_desc(self, desc):
        d = desc or ""
        return bool(self.foreign_rx and self.foreign_rx.search(d) and not (self.not_foreign_rx and self.not_foreign_rx.search(d)))


RIBO_RX = re.compile(r"^([eub])([ls])(\d+[a-z]?)$")
MAX_NAME_ACCS = 10                                   # a name shared by more accessions identifies nothing


def norm_name(s):
    s = re.sub(r"\s+", " ", str(s or "").strip().strip("'\"").strip()).lower()
    s = re.sub(r",? (isoform|form) [\w-]+$", "", s)
    return "" if s in ("?", ".", "unknown", "uncharacterized protein", "protein") else s


def protein_names(text):
    """UniProt 'Protein names' -> [names]: the recommended name and every parenthesised alternative"""
    t = str(text or "").split(" [")[0]
    main = t.split(" (")[0]
    return [main] + re.findall(r"\(([^()]+)\)", t)


def _tsv(path):
    import pandas as pd
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, comment="#")


def load_policy(cf, resolve):
    """Policy from library.chain_filter (None when absent or enabled: false); resolve(value) -> path"""
    if not cf or cf.get("enabled", True) is False:
        return None
    bad = [k for k in cf if k not in KEYS]
    if bad:
        raise ValueError(f"library.chain_filter: unknown key(s) {bad} (known: {', '.join(KEYS)})")
    import pandas as pd
    acc_genes, acc_primary, gene_acc, name_accs = {}, {}, {}, {}
    if cf.get("genes"):
        g = cf["genes"]; c = g.get("columns", {})
        ca, cg = c.get("acc", "Entry"), list(c.get("genes", ["Gene Names (primary)", "Gene Names (synonym)"]))
        cn, crev = c.get("names"), c.get("reviewed")
        T = pd.read_csv(resolve(g["path"]), sep="\t", usecols=[ca] + cg + [x for x in (cn, crev) if x], dtype=str,
                        keep_default_na=False)
        if cn:                                       # reviewed entries: protein names and primary symbols -> accessions
            R = T[T[crev] == "reviewed"] if crev else T
            for a, nm, gp in zip(R[ca], R[cn], R[cg[0]]):
                for x in protein_names(nm) + [str(gp).split(";")[0]]:
                    k = norm_name(x)
                    if k:
                        name_accs.setdefault(k, set()).add(a)
            name_accs = {k: v for k, v in name_accs.items() if len(v) <= MAX_NAME_ACCS}
        T = T[[ca] + cg]
        for r in T.itertuples(index=False):
            names = set(" ".join(str(x) for x in r[1:]).replace(";", " ").split())
            if names:
                acc_genes[r[0]] = frozenset(names)
            prim = str(r[1]).split(";")[0].strip()
            if prim:
                acc_primary[r[0]] = prim
                gene_acc.setdefault(prim, r[0])
    nucleic, allow_accs, allow_chains = {}, {}, {}
    if cf.get("nucleic_acids"):
        for r in _tsv(resolve(cf["nucleic_acids"])).to_dict("records"):
            re.compile(r["desc_regex"])
            nucleic[r["ing_id"]] = r["desc_regex"]
    if cf.get("allow"):
        for r in _tsv(resolve(cf["allow"])).to_dict("records"):
            if r["what"] == "accession":
                allow_accs.setdefault(r["ing_id"], set()).add(r["value"])
            elif r["what"] == "chain":
                allow_chains.setdefault(r["ing_id"], set()).add(r["value"])
            else:
                raise ValueError(f"chain_filter allow: what must be accession | chain, got {r['what']!r} ({r})")
    par = {}
    if cf.get("paralogs"):
        p = cf["paralogs"]; c = p.get("columns", {})
        cg, cp = c.get("gene", "Gene name"), c.get("paralog", "Human paralogue associated gene name")
        T = pd.read_csv(resolve(p["path"]), sep="\t", usecols=[cg, cp], dtype=str, keep_default_na=False)
        T = T[(T[cg] != "") & (T[cp] != "")]
        for a, b in zip(T[cg], T[cp]):
            par.setdefault(a, set()).add(b)
    pol = Policy(cf.get("gene_match", False), acc_genes, nucleic, allow_accs, allow_chains, par, acc_primary,
                 cf.get("foreign_desc", FOREIGN_DESC), cf.get("not_foreign_desc", NOT_FOREIGN_DESC), gene_acc, name_accs)
    pol.taxonomy_path = resolve(cf["taxonomy"]) if cf.get("taxonomy") else None
    pol.sifts_path = resolve(cf["sifts"]) if cf.get("sifts") else None
    return pol


def load_taxonomy(path, pids):
    """{pdb: {chain: set of NCBI taxonomy ids}} of the SIFTS pdb_chain_taxonomy table for these entries"""
    import pandas as pd
    if not path:
        return {}
    pids = {str(p).lower() for p in pids}
    out = {}
    for T in pd.read_csv(path, sep="\t", comment="#", usecols=["PDB", "CHAIN", "TAX_ID"], dtype=str, chunksize=1 << 20):
        T = T[T.PDB.isin(pids)]
        for p, c, t in zip(T.PDB, T.CHAIN, T.TAX_ID):
            out.setdefault(p, {}).setdefault(c, set()).add(t)
    return out


def species_spec(row):
    """{ing_id, kind, genes (list), accs (set), acc_list (positional, with empty entries)} of a species-table row"""
    a = row.get("accs", "")
    acc_list = [str(x) for x in (a if isinstance(a, (list, tuple)) else str(a).split(";"))]
    accs = [x for x in acc_list if x]
    if row.get("uniprot"):
        accs.append(str(row["uniprot"]))
    g = row.get("genes", "")
    genes = [x for x in (g if isinstance(g, (list, tuple)) else str(g).split(";")) if x]
    return dict(ing_id=row["ing_id"], kind=row.get("kind", ""), genes=genes, accs=set(accs), acc_list=acc_list)


def chain_info(acc, chn, kind):
    """{chain: {n, n_na, accs: {acc: residues}, n_unmapped}} in order of first appearance"""
    out = {}
    for a, c, k in zip(acc, chn, kind):
        d = out.setdefault(str(c), dict(n=0, n_na=0, accs={}, n_unmapped=0))
        d["n"] += 1
        if int(k) == 1:
            d["n_na"] += 1
        elif a:
            d["accs"][str(a)] = d["accs"].get(str(a), 0) + 1
        else:
            d["n_unmapped"] += 1
    return out


def allowed_accs(spec, policy):
    return set(spec["accs"]) | (policy.allow_accs.get(spec["ing_id"], set()) if policy else set())


def _incomplete(spec):
    al = spec.get("acc_list") or []
    return [g for g, a in zip(spec["genes"], al) if not a] if len(al) == len(spec["genes"]) else []


def _member_genes(spec, policy):
    g = set(spec["genes"])
    g |= {policy.acc_primary[a] for a in allowed_accs(spec, policy) if a in policy.acc_primary}
    return g


def classify_protein(spec, policy, accs, desc="", taxa=None, n=0):
    """class of one protein chain: (cls in member | foreign | flagged, reason, member accessions hit)"""
    A = allowed_accs(spec, policy)
    m = set(a for a in accs if a)
    genes = set(spec["genes"])
    hit = m & A
    if hit:
        return "member", "accession", hit
    if m and policy.gene_match:
        g = {a for a in m if policy.acc_genes.get(a, frozenset()) & genes}
        if g:
            return "member", "gene_match", g
    if policy.foreign_by_desc(desc):
        return "foreign", "antibody_or_binder_by_description", set()
    non_human = bool(taxa) and HUMAN not in taxa
    named = policy.accs_by_name(desc) if (non_human or not m) else set()
    if named & A or (policy.gene_match and any(policy.acc_genes.get(a, frozenset()) & genes for a in named)):
        return "member", "name_match:" + ",".join(sorted(named))[:60], set()   # an ortholog / unmapped copy of a member
    cand = m if (m and not non_human) else named
    if not cand:
        if not m:
            return "flagged", ("unmapped_peptide" if 0 < n < 30 else "unmapped"), set()
        return "flagged", "non_human_possible_ortholog:" + ",".join(sorted(taxa)), set()
    mg = _member_genes(spec, policy)
    par = {policy.acc_primary.get(a, "") for a in cand} & {p for g in mg for p in policy.paralogs.get(g, ())}
    par.discard("")
    if par:
        return "flagged", "paralogue_of_member:" + ",".join(sorted(par)), set()
    inc = _incomplete(spec)
    if inc:
        return "flagged", "species_accessions_incomplete:" + ",".join(inc), set()
    return "foreign", "other_protein" if cand is m else "other_protein_by_name:" + ",".join(sorted(cand))[:60], set()


RRNA_RX = re.compile(r"\brrna\b|ribosomal rna|\b(5|5\.8|12|16|18|23|25|26|28|45|47)\s?s\b|\b[35]'?\s?ets\b|\bits\s?[12]\b",
                     re.I)
RIBO_CTX_RNA_RX = re.compile(r"\btrna|transfer rna|\bmrna\b|messenger rna|\bp-site|\be-site|\ba-site", re.I)


def ribosome_context(info, desc):
    """the structure holds ribosomal RNA (a ribosome / pre-ribosome in the crystal)"""
    return any(d["n_na"] == d["n"] and RRNA_RX.search(desc.get(c, "") or "") for c, d in info.items()) if desc else False


def classify_nucleic(spec, policy, desc, ribo_ctx):
    """class of one nucleic-acid chain (entity description): declared for the species (nucleic_acids) -> member;
    ribosomal RNA of a species that does not declare it -> foreign (a ribosome in the crystal); tRNA / mRNA in a
    structure holding rRNA -> foreign (ribosome context); any other undeclared nucleic acid -> flagged (kept)"""
    if policy.nucleic_ok(spec["ing_id"], desc):
        return "member", "declared_nucleic_acid"
    if RRNA_RX.search(desc or ""):
        return "foreign", "ribosomal_rna_not_member"
    if ribo_ctx and RIBO_CTX_RNA_RX.search(desc or ""):
        return "foreign", "ribosome_context_rna"
    return "flagged", "nucleic_acid_not_declared"


def decide(spec, info, policy, desc=None, tax=None):
    """classes of the chains of one structure. desc: {chain: entity description}; tax: {chain: set of taxonomy ids}.
    Returns {keep, na_keep, drop_acc, kept, removed, flagged, gene_matched, missing, unverifiable}"""
    ing, kind = spec["ing_id"], spec["kind"]
    chains_ok = policy.allow_chains.get(ing, set())
    desc, tax = desc or {}, tax or {}
    A = allowed_accs(spec, policy)
    ribo = ribosome_context(info, desc)
    keep, na_keep, drop_acc, kept, removed, gm = [], [], {}, [], [], []
    for c, d in info.items():
        n_p = d["n"] - d["n_na"]
        accs = d["accs"]
        ds = desc.get(c, "")
        base = dict(chain=c, desc=ds[:60])
        if c.split("-")[0] in chains_ok or c in chains_ok:
            keep.append(c); na_keep.append(c)
            kept.append(dict(base, acc=";".join(sorted(accs)) or "unmapped", n=d["n"], cls="member", why="declared_chain",
                             type="na" if n_p == 0 else "protein"))
            continue
        if n_p == 0:                                   # nucleic-acid chain (classify_nucleic)
            cls, why = classify_nucleic(spec, policy, ds, ribo)
            if cls == "foreign":
                removed.append(dict(base, acc="", n=d["n"], type="na", cls="foreign", why=why))
                continue
            keep.append(c); na_keep.append(c)
            kept.append(dict(base, acc="", n=d["n"], type="na", cls=cls, why=why))
            continue
        cls, why, hit = classify_protein(spec, policy, list(accs), ds, tax.get(c.split("-")[0]), n_p)
        if cls == "foreign":
            removed.append(dict(base, acc=";".join(sorted(accs)) or "unmapped", n=d["n"], cls="foreign", why=why,
                                type="protein"))
            continue
        keep.append(c); na_keep.append(c)              # nucleotides inside a kept protein chain stay with it
        if why == "gene_match":
            gm.extend(dict(chain=c, acc=a, genes=";".join(sorted(policy.acc_genes[a] & set(spec["genes"])))) for a in sorted(hit))
        fus = []
        if cls == "member" and not d.get("chain_level"):   # residues of a member chain mapped to a positively foreign protein
            for a in sorted(set(accs) - hit - A):
                c2, w2, _ = classify_protein(spec, policy, [a], "", tax.get(c.split("-")[0]), accs[a])
                if c2 == "foreign":
                    fus.append(a)
        if fus:
            drop_acc[c] = fus
            for a in fus:
                removed.append(dict(base, acc=a, n=accs[a], cls="foreign", why="fusion_residues", type="protein"))
        kept.append(dict(base, acc=";".join(sorted(hit)) if hit else (";".join(sorted(accs)) or "unmapped"),
                         n=d["n"] - sum(accs[a] for a in fus), cls=cls, why=why, type="protein"))
    if kind == "protein":                              # a monomer model: one chain of the protein
        prot = [k for k in kept if k["type"] == "protein" and k["cls"] == "member"]
        if len(prot) > 1:
            best = max(prot, key=lambda k: (k["n"], -prot.index(k)))
            for k in prot:
                if k is not best:
                    keep.remove(k["chain"]); kept.remove(k)
                    removed.append(dict(k, cls="extra_copy", why="extra_copy_monomer"))
    present = {a for k in kept if k["cls"] == "member" for a in k["acc"].split(";") if a and a != "unmapped"}
    missing = sorted(a for a in spec["accs"] if a not in present)
    return dict(keep=keep, na_keep=na_keep, drop_acc=drop_acc, kept=kept, removed=removed,
                flagged=[k for k in kept if k["cls"] == "flagged"], gene_matched=gm, missing=missing,
                unverifiable=_incomplete(spec))


def mask(acc, chn, kind, dec):
    keep, na_keep = set(dec["keep"]), set(dec["na_keep"])
    drop = dec["drop_acc"]
    m = np.zeros(len(chn), bool)
    for i, (a, c, k) in enumerate(zip(acc, chn, kind)):
        c = str(c)
        if c not in keep:
            continue
        if int(k) == 1:
            m[i] = c in na_keep
        else:
            m[i] = not (c in drop and str(a) in drop[c])
    return m


def record(dec, n_struct, source, n_in):
    """JSON record of the decision (npz key chain_set)"""
    return dict(source=source, n_in=int(n_in), n_struct=int(n_struct), kept=dec["kept"], removed=dec["removed"],
                gene_matched=dec["gene_matched"], missing=dec["missing"], unverifiable=dec["unverifiable"])


def has_member(rec):
    """a model needs at least one member chain (a structure of flagged chains only is not a model of the species)"""
    return any(k.get("cls", "member") == "member" for k in rec["kept"])


def npz_extra(rec, beads_added=0):
    return dict(chain_set=np.array(json.dumps(rec, sort_keys=True)), n_struct_beads=np.int32(rec["n_struct"]),
                beads_added=np.int32(beads_added))


def entity_descriptions(path):
    """{chain name: _entity.pdbx_description} of an mmCIF (assembly chain names as gemmi reads them)"""
    import gemmi
    blk = gemmi.cif.read(path).sole_block()
    desc = {r[0]: gemmi.cif.as_string(r[1]) for r in blk.find("_entity.", ["id", "pdbx_description"])}
    st = gemmi.make_structure_from_block(blk); st.setup_entities()
    sub2ent = {s: e.name for e in st.entities for s in e.subchains}
    out = {}
    for ch in st[0]:
        for res in ch:                                 # the polymer residues' entity (a chain also holds ligands, water)
            info = gemmi.find_tabulated_residue(res.name)
            if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                out[ch.name] = desc.get(sub2ent.get(res.subchain, ""), "")
                break
    return out


def afdb_chain_accs(path):
    """{chain id: UniProt accession} of an AFDB (ModelCIF) complex model from _struct_ref / _struct_ref_seq"""
    import gemmi
    blk = gemmi.cif.read(path).sole_block()
    ref = {r[0]: r[1] for r in blk.find("_struct_ref.", ["id", "pdbx_db_accession"])}
    return {r[1]: ref.get(r[0], "") for r in blk.find("_struct_ref_seq.", ["ref_id", "pdbx_strand_id"])}


def chain_level(info, chain_accs):
    """chains whose residues SIFTS maps to nothing (label_seq outside every SIFTS range: other numbering / isoform) take
    the chain-level SIFTS accessions (chain_accs {auth chain: set}); fusion trimming is not possible for them"""
    for c, d in info.items():
        n_p = d["n"] - d["n_na"]
        ca = (chain_accs or {}).get(c.split("-")[0])
        if n_p and not d["accs"] and ca:
            d["accs"] = {a: n_p for a in sorted(ca)}
            d["n_unmapped"] = 0
            d["chain_level"] = True
    return info


def apply(xyz, kind, acc, pos, chn, spec, policy, source, acc_decide=None, desc_fn=None, tax=None, chain_accs=None):
    """the rule on the bead arrays of one structure (per-bead chain chn; acc_decide: per-bead accession used for the
    decision when it differs from acc, e.g. the AFDB model's own chain accession; desc_fn() -> {chain: description},
    read only when a chain is not a member; tax: {chain: taxonomy ids}; chain_accs: chain-level SIFTS accessions,
    chain_level). Returns the filtered arrays and the record."""
    acc_d = acc if acc_decide is None else acc_decide
    info = chain_level(chain_info(acc_d, chn, kind), chain_accs)
    dec = decide(spec, info, policy, None, tax)
    if desc_fn is not None and any(k["cls"] != "member" for k in dec["kept"] + dec["removed"]):
        dec = decide(spec, info, policy, desc_fn(), tax)
    m = mask(acc_d, chn, kind, dec)
    rec = record(dec, int(m.sum()), source, len(m))
    return xyz[m], kind[m], acc[m], pos[m], chn[m], rec


def summary_cols(rec):
    """ledger columns of a chain record"""
    if not rec:
        return {}
    f = lambda d: f"{d['chain']}:{d['acc'] or 'na'}:{d['n']}"     # noqa: E731
    mem = [k for k in rec["kept"] if k.get("cls", "member") == "member"]
    fl = [k for k in rec["kept"] if k.get("cls") == "flagged"]
    return dict(chains_kept=",".join(f(d) for d in mem),
                chains_flagged=",".join(f"{f(d)}:{d['why']}" + (f"[{d['desc'][:30]}]" if d.get("desc") else "") for d in fl),
                chains_removed=",".join(f"{f(d)}:{d['why']}" + (f"[{d['desc'][:30]}]" if d.get("desc") else "") for d in rec["removed"]),
                subunits_missing=";".join(rec.get("missing", [])), accs_unverifiable=";".join(rec.get("unverifiable", [])),
                gene_matched=";".join(f"{d['chain']}:{d['acc']}:{d['genes']}" for d in rec.get("gene_matched", [])))


# ------------------------------------------------------------------ check (library stage fails on a violation)
def _sifts_accs(sifts_path, pids):
    import pandas as pd
    S = pd.read_csv(sifts_path, sep="\t", comment="#", usecols=["PDB", "CHAIN", "SP_PRIMARY"], dtype=str)
    S = S[S.PDB.isin(set(pids))]
    out = {}
    for p, c, a in zip(S.PDB, S.CHAIN, S.SP_PRIMARY):
        out.setdefault((p, c), set()).add(a)
    return out


def check(lib_dir, species_rows, policy, sifts_path, afdb_dir=None, taxonomy_path=None, skip_origins=("start",),
          origin=None, log=print):
    """every model of lib_dir built from a structure file (npz src pdb:<id> | afdb_complex:<model> |
    fallback:pdb_chain:<id>:<chain>): a chain record exists; no kept chain is class 2 (foreign) - its class recomputed
    from the SIFTS accessions / taxonomy of the entry (AFDB: the model's struct_ref), the recorded entity description and
    the policy; a monomer keeps one protein chain; the record adds up to n_struct_beads; the model holds no more beads
    than n_struct_beads + beads_added. Returns (rows, errors)."""
    origin = origin or {}
    rows, errors, todo = [], [], []
    for r in species_rows:
        f = os.path.join(lib_dir, f"{r['ing_id']}.npz")
        if not os.path.exists(f) or origin.get(f"{r['ing_id']}.npz") in skip_origins:
            continue
        with np.load(f, allow_pickle=False) as z:
            src = str(z["src"]) if "src" in z.files else ""
            if not src.startswith(STRUCT_PREFIXES):
                continue
            rec = json.loads(str(z["chain_set"])) if "chain_set" in z.files else None
            nb = len(z["beads"])
            ns = int(z["n_struct_beads"]) if "n_struct_beads" in z.files else None
            na = int(z["beads_added"]) if "beads_added" in z.files else 0
        todo.append((r, src, rec, nb, ns, na))

    def pid_of(s):
        return s.split(":")[1] if s.startswith("pdb:") else s.split(":")[2]
    pids = {pid_of(s) for _, s, _, _, _, _ in todo if s.startswith(("pdb:", "fallback:pdb_chain:"))}
    SI = _sifts_accs(sifts_path, pids) if pids and sifts_path else {}
    TX = load_taxonomy(taxonomy_path, pids) if pids and taxonomy_path else {}
    n_flag = 0
    for r, src, rec, nb, ns, na in todo:
        ing = r["ing_id"]
        spec = species_spec(r)
        bad = []
        if rec is None:
            errors.append((ing, src, "no chain record (chain_set) in a structure-derived model"))
            continue
        kept = rec["kept"]
        fa = None
        if src.startswith("afdb_complex:") and afdb_dir:
            fa = afdb_chain_accs(os.path.join(afdb_dir, src.split(":", 1)[1] + ".cif.gz"))
        for k in kept:
            c = k["chain"]
            if k.get("why") == "declared_chain":
                if c.split("-")[0] not in policy.allow_chains.get(ing, set()) and c not in policy.allow_chains.get(ing, set()):
                    bad.append(f"{c}: declared_chain not declared")
                continue
            if k["type"] == "na":                       # rRNA / ribosome-context RNA of a non-ribosome species: class 2
                ctx = any(RRNA_RX.search(x.get("desc", "") or "") for x in kept + rec["removed"] if x["type"] == "na")
                cls, why = classify_nucleic(spec, policy, k.get("desc", ""), ctx)
                if cls == "foreign":
                    bad.append(f"{c}: nucleic acid '{k.get('desc', '')}' {why} (class 2) still in the model")
                continue
            if src.startswith("afdb_complex:"):
                m, tx = ({fa.get(c, "")} if fa is not None else set(k["acc"].split(";"))), {HUMAN}
            else:
                key = (pid_of(src), c.split("-")[0])
                m, tx = SI.get(key, set()), TX.get(key[0], {}).get(key[1])
            cls, why, _ = classify_protein(spec, policy, sorted(m), k.get("desc", ""), tx, int(k["n"]))
            if cls == "foreign":
                bad.append(f"{c}: {sorted(m) or 'unmapped'} {why} (class 2) still in the model")
            n_flag += cls == "flagged"
        nprot = sum(1 for k in kept if k["type"] == "protein" and k.get("cls", "member") == "member")
        if spec["kind"] == "protein" and nprot != 1:
            bad.append(f"monomer model with {nprot} protein chains")
        if ns is None or sum(int(k["n"]) for k in kept) != ns or ns != int(rec["n_struct"]):
            bad.append(f"recorded chain residues {sum(int(k['n']) for k in kept)} != n_struct_beads {ns}")
        elif nb > ns + na:
            bad.append(f"{nb} beads > {ns} structure beads of the kept chains + {na} added beads")
        rows.append(dict(ing_id=ing, src=src, beads=nb, n_struct=ns, added=na, kept=len(kept),
                         flagged=sum(1 for k in kept if k.get("cls") == "flagged"), removed=len(rec["removed"]), ok=not bad))
        if bad:
            errors.append((ing, src, "; ".join(bad)))
    log(f"chain check: {len(rows)} structure-derived models, {len(errors)} violations ({n_flag} flagged chains kept)")
    return rows, errors
