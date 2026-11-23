"""Assembly steps on an ingredient table with fractional copies (option `pools: ingredients`; preset legacy_x5).

The same ordered step list as steps.py (curated_machines -> complex_portal -> afdb_heterodimers -> pdb_homo_oligomers),
here in the form the yeast X5 build used: the pools are the rows of an ingredient table (`ing`: one row per packed
species, kind protein | complex, float copies_per_cell), every step consumes copies from it and adds its own rows.
Plugin steps (registry kind `assembly_step`, e.g. the yeast pre-ribosomes) take and return the same table. Structure
building (beads, own RNA, membrane orientation) is not done here: a step that needs it calls the structure hook named
by its `structure.builder` option (registry kind `assembly_hook`); with builder null every assembly is accepted.

Ported from the vendored yeast engine (plugins/yeast/engine, X5 code e3b7448): every function carries its origin.
The arithmetic is the yeast one line for line (preset legacy_x5 reproduces the yeast port build s0 tables; check:
scripts/assembly_x5_check.py); constants are preset options and the yeast lists are cell data
(inventory.params.assembly_steps.data, TSVs with a reason column under plugins/yeast/data/assembly/).

Inputs (params.assembly_steps): inputs.protein_table (the stage-1 protein table: gene, uniprot, name, copies_per_cell,
mass_da, comp_weights, n_tm, oligomer_n), inputs.complex_portal (Complex Portal TSV, identifiers in column 5),
inputs.rcsb_cache (frozen RCSB responses, read-only: a miss is an error), inputs.pdb_assemblies (local assembly1
files), inputs.afdb_heterodimers / afdb_homodimers (AFDB candidates.csv), inputs.intact_pairs, inputs.paralog_pairs,
inputs.proteome_fasta.
"""
from __future__ import annotations

import collections
import io
import json
import os
import re

import numpy as np
import pandas as pd

from .steps import read_curated


def _lsum(xs):
    """left-to-right sum, as sum() of CPython <= 3.11 (the yeast build). CPython 3.12 sums floats with compensation,
    which moves results by an ULP (seen on Complex Portal coverages and compartment weights)"""
    t = 0
    for x in xs:
        t = t + x
    return t


# ------------------------------------------------------------------------------------------------ state
class IngredientState:
    """ctx, options o, params, the protein table df, the ingredient table ing, hooks, and what later steps read"""

    def __init__(self, ctx, params, o, df=None):
        self.ctx, self.o, self.params = ctx, o, params
        self.inputs = dict(params.get("inputs") or {})
        self.data = dict(params.get("data") or {})
        self.df = df if df is not None else read_protein_table(self.frozen("protein_table"))
        self.ing = None
        self.modeled = set()          # genes of the curated machines + legacy oligomers (complex_portal skips them)
        self.legacy = {}              # gene -> (lib_id, n): homo-oligomers kept as modeled
        self.oligo_map = None
        self.view = None              # SEC-MS view (wholecell.inventory.secms.SecView) when the preset selects one
        self._rcsb = None
        self._hooks = {}

    def frozen(self, key):
        """inputs.<key>: a declared frozen input name, 'name/relative path' inside one, 'build:...' or a path"""
        v = self.inputs.get(key)
        if v is None:
            raise KeyError(f"inventory.params.assembly_steps.inputs.{key} is required")
        head, _, rest = v.partition("/")
        if rest and head in (self.ctx.cfg.get("frozen_inputs") or {}):
            return os.path.join(self.ctx.frozen(head), rest)
        return self.ctx.resolve(v)

    def hook(self, name):
        if name not in self._hooks:
            from ... import registry
            h = registry.get("assembly_hook", name)
            self._hooks[name] = h(self) if isinstance(h, type) else h
        return self._hooks[name]

    def rcsb(self):
        if self._rcsb is None:
            self._rcsb = RcsbCache(self.frozen("rcsb_cache"))
        return self._rcsb

    def out(self, name):
        return self.ctx.out(f"assembly/{name}")

    def write_csv(self, df, name):
        from ... import build
        build.atomic_save(self.out(name), lambda t: df.to_csv(t, index=False))

    def write_json(self, obj, name, **kw):
        from ... import build
        build.atomic_save(self.out(name), lambda t: json.dump(obj, open(t, "w"), **kw))


def read_protein_table(path):
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path, sep="\t" if path.endswith(".tsv") else ",")


class RcsbCache:
    """frozen RCSB responses (the yeast rcsb_cache.json layout: 'acc:<UniProt>' -> entries, 'entry:<id>' -> {uniprot, res},
    'ent:<id>' -> {entity: [UniProt]}, 'homo:<acc>:<n>' -> entry or null); read-only, offline"""

    def __init__(self, path):
        self.path = path
        self.c = json.load(open(path))

    def get(self, key):
        if key not in self.c:
            raise KeyError(f"RCSB response '{key}' is not in the frozen cache {self.path} (offline build: refresh the "
                           f"cache with the yeast engine and declare the new file)")
        return self.c[key]

    def entries_with(self, acc):
        # origin: yeast pipeline/assemblies.py:entries_with @ e3b7448
        return self.get(f"acc:{acc}")

    def entry_info(self, ids):
        # origin: yeast pipeline/assemblies.py:entry_info @ e3b7448
        return {i: self.c.get(f"entry:{i}") for i in ids}

    def entity_uniprot(self, pid):
        # origin: yeast pipeline/assemblies.py:entity_uniprot @ e3b7448
        return {k: set(v) for k, v in self.get(f"ent:{pid}").items()}


def _mass_in_assemblies(ing, oligo_map):
    # origin: yeast scripts/s1_inventory.py:mass_in_assemblies @ e3b7448
    m = ing.copies_per_cell * ing.mass_da
    asmb = (ing.kind == "complex") | ing.ing_id.isin(list(oligo_map))
    return float(m[asmb].sum() / m.sum())


# ------------------------------------------------------------------------------------------------ curated machines
def _machine_families(spec_rows, cp, gene_of_acc, genes):
    """{family: (stoich, [paralog genes])} of one machine from its rows of data.curated_machine_subunits:
    members = 'cpx:<Complex Portal id>' (subunits and stoichiometry from Complex Portal), 're:<regex>' (every gene of
    the pool that matches, in pool order, family = gene with family_sub removed, stoich 1), or paralog genes
    separated by '|' (pooled into one family)
    origin: yeast pipeline/ingredients.py:complex_families @ e3b7448"""
    out = {}
    for r in spec_rows.itertuples(index=False):
        m = r.members
        if m.startswith("cpx:"):
            row = cp[cp["#Complex ac"] == m[4:]].iloc[0]
            for tok in row["Identifiers (and stoichiometry) of molecules in complex"].split("|"):
                mm = re.match(r"([A-Z0-9]+)(?:-\w+)?\((\d+)\)", tok)
                if mm and gene_of_acc.get(mm.group(1)):
                    g = gene_of_acc[mm.group(1)]
                    out[g] = (max(1, int(mm.group(2))), [g])
        elif m.startswith("re:"):
            pat, sub = re.compile(m[3:]), r.family_sub
            for g in genes:
                if pat.match(g):
                    f = re.sub(sub, "", g) if sub else g
                    if f in out:
                        out[f][1].append(g)
                    else:
                        out[f] = (int(r.stoich), [g])
        else:
            gs = m.split("|")
            if r.only_if_present == "yes":
                gs = [g for g in gs if g in genes]
                if not gs:
                    continue
            out[r.family or gs[0]] = (int(r.stoich), gs)
    return out


def _subunit_weights(df, genes, keep):
    """copy-weighted mean compartment weights of a complex's subunits, restricted to `keep` and renormalised
    origin: yeast pipeline/ingredients.py:subunit_weights @ e3b7448"""
    acc = {k: 0.0 for k in keep}
    for g in genes:
        r = df[df.gene == g]
        if len(r) == 0:
            continue
        w = json.loads(r.comp_weights.iloc[0]); n = float(r.copies_per_cell.iloc[0])
        for k in keep:
            acc[k] += n * w.get(k, 0.0)
    z = _lsum(acc.values())
    return {k: round(v / z, 4) for k, v in acc.items()} if z > 0 else {k: 1.0 / len(keep) for k in keep}


def curated_machines(st):
    """curated machines from data.curated_machines (one row per machine: id, label, PDB, compartment weights, membrane,
    copies rule) and data.curated_machine_subunits (paralog families and stoichiometry); copies by rule
    (median_families = median over families of family pool / stoich; family:<name>; hook:<assembly_hook>); each
    paralog consumes its share of the family pool. Then the residual monomer pool (free pool of the genes listed in
    free_pool_cap capped at max_frac_of_assembled of their copies in the machine; surplus not placed) and the legacy
    homo-oligomers (data.legacy_oligomers, kept as modeled when the gene is in the pool and the library has the entry)
    origin: yeast pipeline/ingredients.py:build_ingredients and scripts/s1_inventory.py:main @ e3b7448"""
    o, df = st.o["curated_machines"], st.df
    M = read_curated(st.data["curated_machines"])
    SUB = read_curated(st.data["curated_machine_subunits"])
    for c in ("family", "family_sub", "only_if_present"):
        if c not in SUB.columns:
            SUB[c] = ""
    cp = pd.read_csv(st.frozen("complex_portal"), sep="\t")
    gene_of_acc = dict(zip(df.uniprot, df.gene))
    pool = dict(zip(df.gene, df.copies_per_cell))
    mass = dict(zip(df.gene, df.mass_da))
    rows, consumed = [], {}
    for m in M.itertuples(index=False):
        fams = _machine_families(SUB[SUB.machine == m.machine], cp, gene_of_acc, list(pool))
        fam_copies = {f: _lsum(pool.get(g, 0) for g in gs) / s for f, (s, gs) in fams.items()}
        rule = m.copies
        if rule.startswith("hook:"):
            copies = st.hook(rule[5:])(st)
        elif rule.startswith("family:"):
            copies = fam_copies[rule[7:]]
        elif rule == "median_families":
            vals = [v for v in fam_copies.values() if v > 0]
            copies = float(np.median(vals)) if vals else 0.0
        else:
            raise ValueError(f"data.curated_machines: machine {m.machine}: copies '{rule}' (allowed: median_families, "
                             f"family:<family>, hook:<assembly_hook>)")
        m_tot = 0.0
        for f, (s, gs) in fams.items():
            tot = _lsum(pool.get(g, 0) for g in gs)
            m_tot += s * np.mean([mass.get(g, 0) for g in gs])
            for g in gs:
                share = pool.get(g, 0) / tot if tot > 0 else 1 / len(gs)
                consumed[g] = consumed.get(g, 0) + copies * s * share
        genes_all = sorted(g for _, gs in fams.values() for g in gs)
        cw = m.comp_weights
        if cw.startswith("from_subunits:"):
            cw = _subunit_weights(df, genes_all, tuple(cw[len("from_subunits:"):].split(",")))
        else:
            cw = json.loads(cw)
        rows.append(dict(ing_id=m.machine, label=m.label, kind="complex", genes=";".join(genes_all),
                         n_subunit_types=len(fams), copies_per_cell=copies, comp_weights=json.dumps(cw), mass_da=m_tot,
                         pdb=m.pdb, membrane=m.membrane or "", oligomer_n=1, uniprot=""))
    # residual monomer pool
    cap = o.get("free_pool_cap") or {}
    cap_genes, cap_re, cap_f = set(), None, None
    if cap.get("machine"):
        cap_genes = set(next(x["genes"] for x in rows if x["ing_id"] == cap["machine"]).split(";"))
        cap_re, cap_f = re.compile(cap["gene_regex"]), float(cap["max_frac_of_assembled"])
    surplus = [0, 0.0, 0.0]
    for r in df.itertuples():
        left = max(0.0, r.copies_per_cell - consumed.get(r.gene, 0.0))
        if r.gene in cap_genes and cap_re.match(str(r.gene)) and left > cap_f * consumed.get(r.gene, 0.0):
            c = cap_f * consumed.get(r.gene, 0.0)
            surplus[0] += 1; surplus[1] += left - c; surplus[2] += (left - c) * r.mass_da; left = c
        if left < 1:
            continue
        rows.append(dict(ing_id=r.gene, label=r.name, kind="protein", genes=r.gene, n_subunit_types=1,
                         copies_per_cell=left, comp_weights=r.comp_weights, mass_da=r.mass_da, pdb="",
                         membrane="", oligomer_n=r.oligomer_n, uniprot=r.uniprot, n_tm=r.n_tm))
    ing = pd.DataFrame(rows)
    ing["n_tm"] = ing.n_tm.fillna(0).astype(int)
    if cap:
        st.ctx.say(f"nucleoporin surplus not placed: {surplus[1]:.0f} copies of {surplus[0]} nucleoporins, "
                   f"{surplus[2] / 1e9:.2f} GDa (free pool capped at {cap_f:.0%} of the copies in pores)")
    # legacy homo-oligomers kept as modeled
    legacy = {}
    if st.data.get("legacy_oligomers"):
        i1 = ing.set_index("ing_id")
        builder = o.get("legacy_oligomers_require_library") and st.hook(o["legacy_oligomers_require_library"])
        for g in read_curated(st.data["legacy_oligomers"]).gene:
            lib_id = g + "_oligo"
            if g in i1.index and (not builder or builder.has_entry(lib_id)):
                legacy[g] = (lib_id, int(i1.loc[g, "oligomer_n"]))
    st.legacy = legacy
    st.modeled = set(g for gs in ing[ing.kind == "complex"].genes for g in gs.split(";")) | set(legacy)
    st.ing = ing
    st.ctx.say(f"curated complexes {int((ing.kind == 'complex').sum())}, legacy homo-oligomers {len(legacy)}; "
               f"mass in assemblies {_mass_in_assemblies(ing, legacy):.3f}")


# ------------------------------------------------------------------------------------------------ Complex Portal
RNA_TERMS = re.compile(r"URS[0-9A-F]+|CHEBI:13333[56]")          # RNAcentral ids; box C/D, H/ACA snoRNA classes
RNA_NAMES = re.compile(r"snoRNP|ribonucleoprotein|ribonuclease P|ribonuclease MRP|recognition particle|telomerase|snRNP|spliceosom", re.I)


def _complex_table(cp, df):
    """Complex Portal rows with >= 2 subunits mapped to the protein table; tm = any subunit with a TM segment, rna = the
    complex lists an RNA component (its own RNA chains are kept by the structure builder)
    origin: yeast pipeline/assemblies.py:complex_table @ e3b7448"""
    acc2gene = dict(zip(df.uniprot, df.gene)); ntm = dict(zip(df.gene, df.n_tm))
    rows = []
    for r in cp.itertuples(index=False):
        sub = {}
        for a, s in re.findall(r"([A-Z0-9]{6,10})(?:-\w+)?\((\d+)\)", r[4]):
            g = acc2gene.get(a)
            if g:
                sub[g] = max(1, int(s))
        if len(sub) < 2:
            continue
        rows.append(dict(cpx=r[0], name=r[1], sub=sub, tm=any(ntm.get(g, 0) > 0 for g in sub),
                         rna=bool(RNA_TERMS.search(r[4]) or RNA_NAMES.search(r[1]) or re.search(r"ribosomal subunit", r[1], re.I))))
    return rows


def _complex_copies(pool, sub, exclude):
    """median over subunits (except `exclude`, e.g. mtDNA-encoded ones that MS under-detects) of pool / stoichiometry,
    capped at 2x the limiting subunit (25th percentile for >= 4 subunit types)
    origin: yeast pipeline/assemblies.py:complex_copies @ e3b7448"""
    per = [pool.get(g, 0) / s for g, s in sub.items() if g not in exclude] or [pool.get(g, 0) / s for g, s in sub.items()]
    lim = min(per) if len(per) < 4 else float(np.percentile(per, 25))
    return float(min(np.median(per), 2 * lim))


def _assembly_coverage(fn, sub, acc_of, mass, ent_acc):
    """residue coverage (mass-weighted observed residues of the best chain / (mass / 110 Da)), resolution, polymer
    residues, the structure and the PDB chain-count stoichiometry {gene: chains} of an assembly1 file
    origin: yeast pipeline/extend_complexes.py:coverage @ e3b7448"""
    import gemmi
    st = gemmi.read_structure(fn); st.setup_entities()
    b = gemmi.cif.read(fn).sole_block()
    present, obs = {}, {}
    nres = 0
    for ch in st[0]:
        for span in ch.subchains():
            ent = st.get_entity_of(span)
            if ent is None or ent.entity_type != gemmi.EntityType.Polymer:
                continue
            nres += len(span)
            for a in ent_acc.get(ent.name, ()):
                present[a] = present.get(a, 0) + 1
                obs[a] = max(obs.get(a, 0), len(span))
    tot = _lsum(mass[g] for g in sub)
    cov = _lsum(mass[g] * min(1.0, obs.get(acc_of[g], 0) / (mass[g] / 110.0)) for g in sub) / tot
    pdb_sub = {g: present[acc_of[g]] for g in sub if present.get(acc_of[g], 0) > 0}
    res = None
    for tag in ["_refine.ls_d_res_high", "_em_3d_reconstruction.resolution", "_reflns.d_resolution_high"]:
        v = b.find_value(tag)
        try:
            res = float(v); break
        except (TypeError, ValueError):
            pass
    return cov, res, nres, st, pdb_sub


def _local_assembly(st, pid):
    fn = os.path.join(st.frozen("pdb_assemblies"), f"{pid}-assembly1.cif.gz")
    return fn if os.path.exists(fn) else None


def _best_assembly(st, sub, acc_of, mass, so, log, prefer=None):
    """RCSB entries holding the subunits (top `candidates` by subunit mass), kept when the subunits they contain make
    >= min_coverage of the complex mass; the best `verify_top` by resolution (a preferred entry first) are checked on
    the assembly1 file: subunits present >= min_coverage of the mass, residue coverage >= min_residue_coverage,
    < max_residues residues; the first passing wins. Returns (pid, res, residue_cov, structure, pdb_sub) or None.
    origin: yeast pipeline/assemblies.py:best_assembly @ e3b7448"""
    R = st.rcsb()
    hits = collections.Counter()
    for g in sub:
        for e in R.entries_with(acc_of[g]):
            hits[e] += mass[g] * sub[g]
    cand = [e for e, _ in hits.most_common(int(so["candidates"]))]
    info = R.entry_info(cand)
    tot = _lsum(mass[g] * sub[g] for g in sub)
    scored = []
    for e in cand:
        inf = info.get(e)
        if not inf: continue
        cov = _lsum(mass[g] * sub[g] for g in sub if acc_of[g] in inf["uniprot"]) / tot
        if cov >= so["min_coverage"]:
            scored.append((inf["res"], e))
    best = None
    order = sorted(scored)[:int(so["verify_top"])]
    if prefer:
        order = [(None, prefer)] + [x for x in order if x[1] != prefer]
    for res, pid in order:
        if res is None:
            res = (R.entry_info([pid]).get(pid) or {}).get("res")
        fn = _local_assembly(st, pid)
        if not fn: continue
        ent = R.entity_uniprot(pid)            # a cache miss is an error, not a skipped entry
        try:
            cov, r_, nres, s, psub = _assembly_coverage(fn, sub, acc_of, mass, ent)
        except Exception as ex:
            log(f"  {pid}: {ex}"); continue
        present = _lsum(mass[g] * sub[g] for g in psub) / tot
        if present >= so["min_coverage"] and cov >= so["min_residue_coverage"] and nres < so["max_residues"]:
            best = (pid, res, cov, s, psub); break
    return best


def complex_portal(st):
    """Complex Portal complexes, greedy by copies x mass: a complex sharing a subunit with an earlier one (or a curated
    machine / legacy oligomer) is skipped (overlap: skip_if_used; data.shared_subunits may be shared); copies by the
    stoichiometry cap rule (data.copies_exclude_subunits left out) on the PDB chain counts of the chosen RCSB
    assembly; membrane complexes need an orientable structure (structure builder). Writes complexes_tried.csv and
    complexes_added.csv; adds the complexes to the ingredient table
    origin: yeast pipeline/assemblies.py:extend_complexes, apply_complexes and scripts/s1_inventory.py:main @ e3b7448"""
    o, df, ing = st.o["complex_portal"], st.df, st.ing
    so = o["structure"]
    log = st.ctx.say
    cp = pd.read_csv(st.frozen("complex_portal"), sep="\t")
    acc_of = dict(zip(df.gene, df.uniprot)); mass = dict(zip(df.gene, df.mass_da)); ntm = dict(zip(df.gene, df.n_tm))
    pool = {r.ing_id: r.copies_per_cell for r in ing[ing.kind == "protein"].itertuples()}
    excl = set(read_curated(st.data["copies_exclude_subunits"]).gene) if st.data.get("copies_exclude_subunits") else set()
    shared = set(read_curated(st.data["shared_subunits"]).gene) if st.data.get("shared_subunits") else set()
    prefer, template = {}, {}
    if st.data.get("structure_overrides"):
        T = read_curated(st.data["structure_overrides"])
        for r in T.itertuples(index=False):
            if r.kind == "prefer_pdb":
                prefer[r.complex_id] = r.value
            elif r.kind == "template_paralog":
                g, par = r.value.split(":")
                template.setdefault(r.complex_id, {})[g] = par
            else:
                raise ValueError(f"data.structure_overrides: kind '{r.kind}' (allowed: prefer_pdb, template_paralog)")
    builder = st.hook(so["builder"]) if so.get("builder") else None
    rows = _complex_table(cp, df)
    for r in rows:
        r["copies"] = _complex_copies(pool, r["sub"], excl)
        r["mass_per_cell"] = r["copies"] * _lsum(mass.get(g, 0) * s for g, s in r["sub"].items())
    rows.sort(key=lambda r: -r["mass_per_cell"])
    chosen, tried, used = [], [], set(st.modeled)
    for r in rows:
        if (set(r["sub"]) & used) - shared:
            tried.append((r["cpx"], r["name"], "overlaps modeled complex", None, None)); continue
        if r["copies"] < 1:
            tried.append((r["cpx"], r["name"], "< 1 copy", None, None)); continue
        acc_c = dict(acc_of)
        for g, par in template.get(r["cpx"], {}).items():
            acc_c[g] = acc_of[par]
        b = _best_assembly(st, r["sub"], acc_c, mass, so, log, prefer=prefer.get(r["cpx"]))
        if b is None:
            tried.append((r["cpx"], r["name"], f"no assembly >= {so['min_coverage']:.0%} coverage", None, None)); continue
        pid, res, cov, s, psub = b
        copies = _complex_copies(pool, psub, excl)
        if copies < 1:
            tried.append((r["cpx"], r["name"], "< 1 copy", pid, cov)); continue
        ing_id = r["cpx"].replace("-", "_") + "_" + re.sub(r"\W+", "_", r["name"])[:32]
        tm = any(ntm.get(g, 0) > 0 for g in psub)
        if builder is not None:
            built, info = builder.complex(ing_id, r["name"], pid, s, psub, acc_c, r["rna"], tm and o["require_orientable_tm"])
        else:
            built, info = (None, 0, 0), {}
        if built is None:
            tried.append((r["cpx"], r["name"], f"membrane orientation failed ({info.get('fail')})", pid, cov)); continue
        nb_, nnt, n_rna = built
        for g, k in psub.items():
            pool[g] = max(0.0, pool.get(g, 0) - copies * k)
        chosen.append(dict(ing_id=ing_id, cpx=r["cpx"], label=r["name"], sub=psub, copies=copies,
                           mass_da=_lsum(mass[g] * k for g, k in psub.items()) + nnt * float(o["nucleotide_mass_da"]),
                           pdb=pid, resolution=res, coverage=cov, n_beads=nb_, n_nt=nnt,
                           n_tm=_lsum(ntm.get(g, 0) * k for g, k in psub.items()), memb_score=info.get("score")))
        used |= set(psub)
        tried.append((r["cpx"], r["name"], "added" + (f" (template: {template[r['cpx']]})" if r["cpx"] in template else ""), pid, cov))
        log(f"{r['cpx']:10s} {r['name'][:42]:42s} {pid} res={res} cov={cov:.2f} copies={copies:.0f} nt={nnt}"
            f"{' membrane ' + str(info) if tm else ''} sub={psub}")
    tried = pd.DataFrame(tried, columns=["complex_ac", "name", "outcome", "pdb", "coverage"])
    st.ing = _apply_complexes(st, ing, df, chosen, cp)
    st.write_csv(tried, "complexes_tried.csv")
    st.write_csv(pd.DataFrame(chosen).drop(columns=["sub"]).assign(subunits=[json.dumps(c["sub"]) for c in chosen]),
                 "complexes_added.csv")
    log(f"{len(chosen)} Complex Portal complexes added; outcomes {tried.outcome.value_counts().to_dict()}")


def _apply_complexes(st, ing, df, chosen, cp):
    """add the chosen complexes and subtract their subunits from the monomer pool; compartment weights = subunit
    weights, each subunit counting 1 / (number of Complex Portal complexes it is in); a generic compartment
    (option complex_portal.resolve_generic) is resolved to the sub-compartments the other subunits name
    origin: yeast pipeline/assemblies.py:apply_complexes @ e3b7448"""
    ing = ing.copy()
    cw_of = dict(zip(df.gene, df.comp_weights))
    n_cpx = collections.Counter(g for r in _complex_table(cp, df) for g in r["sub"])
    gen_k, gen_sub = None, ()
    rg = st.o["complex_portal"].get("resolve_generic")
    if rg:
        gen_k, gen_sub = rg["compartment"], tuple(rg["into"])
    new = []
    for c in chosen:
        w = {}
        spec = {g: 1.0 / max(1, n_cpx.get(g, 1)) for g in c["sub"]}
        z = _lsum(spec.values())
        for g in c["sub"]:
            for k, v in json.loads(cw_of[g]).items():
                w[k] = w.get(k, 0) + v * spec[g] / z
        spec = {k: v for k, v in w.items() if k in gen_sub}
        if gen_k in w and spec:
            gen = w.pop(gen_k); z = _lsum(spec.values())
            for k, v in spec.items():
                w[k] += gen * v / z
        for g, s in c["sub"].items():
            m = (ing.kind == "protein") & (ing.ing_id == g)
            ing.loc[m, "copies_per_cell"] = np.maximum(0, ing.loc[m, "copies_per_cell"] - c["copies"] * s)
        new.append(dict(ing_id=c["ing_id"], label=c["label"], kind="complex", genes=";".join(sorted(c["sub"])),
                        n_subunit_types=len(c["sub"]), copies_per_cell=c["copies"], comp_weights=json.dumps(w),
                        mass_da=c["mass_da"], pdb=c["pdb"], membrane="", oligomer_n=1, uniprot="", n_tm=int(c.get("n_tm", 0))))
    ing = pd.concat([ing, pd.DataFrame(new)], ignore_index=True)
    return ing[(ing.kind != "protein") | (ing.copies_per_cell >= 1)].reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ AFDB heterodimers
def read_fasta(path):
    # origin: yeast pipeline/features.py:read_fasta @ e3b7448
    seqs, acc, buf = {}, None, []
    for line in open(path):
        if line.startswith(">"):
            if acc:
                seqs[acc] = "".join(buf)
            acc = line.split("|")[1]; buf = []
        else:
            buf.append(line.strip())
    if acc:
        seqs[acc] = "".join(buf)
    return seqs


def _k3(a, b):
    """share of the shorter sequence's distinct 3-mers found in the other
    origin: yeast pipeline/heterodimers.py:_k3 @ e3b7448"""
    if len(a) < 3 or len(b) < 3:
        return 0.0
    A = {a[i:i + 3] for i in range(len(a) - 2)}; B = {b[i:i + 3] for i in range(len(b) - 2)}
    return len(A & B) / min(len(A), len(B))


def _family(a, b):
    # origin: yeast pipeline/heterodimers.py:_family @ e3b7448
    ma, mb = re.match(r"^([A-Z]+)(\d+)$", a), re.match(r"^([A-Z]+)(\d+)$", b)
    return bool(ma and mb and ma.group(1) == mb.group(1))


def _compatible(wa, wb, tma, tmb, min_shared, facing, inner, inner_side):
    """compartment weights of the pair, or None if the partners cannot meet: shared weight >= min_shared, or a
    membrane protein whose soluble partner sits on the facing side (cytosol; for the inner membrane `inner` the
    matrix side `inner_side`)
    origin: yeast pipeline/heterodimers.py:compatible @ e3b7448"""
    shared = {c: min(wa.get(c, 0), wb.get(c, 0)) for c in set(wa) & set(wb)}
    if _lsum(shared.values()) >= min_shared:
        z = _lsum(shared.values()); return {c: v / z for c, v in shared.items() if v > 0}
    for wt, ws, tm_t, tm_s in ((wa, wb, tma, tmb), (wb, wa, tmb, tma)):
        if tm_t and not tm_s:
            mem = {c: v for c, v in wt.items() if c in facing or c == inner}
            side = ws.get("cytosol", 0) if inner not in mem else _lsum(ws.get(c, 0) for c in inner_side)
            if mem and side >= min_shared:
                z = _lsum(mem.values()); return {c: v / z for c, v in mem.items()}
    return None


def _pair_evidence(st, ev, h):
    """heterodimer evidence gate as the inventory pools (steps.py afdb_heterodimers, Jurkat legacy_v7): a pair is kept
    when its SEC-MS co-elution r (secms_pairs of the selected run) >= sec_r_min OR its STRING combined score >=
    string_min. ev = 'none' / None: no gate. STRING ids: inputs.uniprot_full column `STRING` ('4932.YAL001C;'), else
    evidence.string_id {column, prefix} (yeast: UniProt 'Gene Names (ordered locus)', prefix '4932.')"""
    if not isinstance(ev, dict):
        return None
    import gzip
    sec = st.view.pair_r() if st.view is not None else (lambda a, b: np.nan)
    smin, rmin = ev.get("string_min"), ev.get("sec_r_min")
    score = {}
    ids = {}
    if smin is not None:
        sid = ev.get("string_id") or {}
        up = pd.read_csv(st.frozen("uniprot_full"), sep="\t")
        col, pre = sid.get("column", "STRING"), sid.get("prefix", "")
        for e, v in zip(up.Entry, up[col]):
            if isinstance(v, str):
                ids[e] = [pre + x.strip() for x in v.replace(";", " ").split() if x.strip()]
        need = {x for a in pd.concat([h.uniprotAccession_1, h.uniprotAccession_2]) for x in ids.get(a, [])}
        with gzip.open(st.frozen("string_links"), "rt") as f:
            next(f)
            for line in f:
                p1, p2, sc = line.split()
                if p1 in need and p2 in need:
                    score[(p1, p2)] = int(sc)

    def string(a, b):
        return max([score.get((x, y), 0) for x in ids.get(a, []) for y in ids.get(b, [])] or [0])

    def ok(a, b):
        r = sec(a, b)
        return bool((rmin is not None and r == r and r >= rmin) or (smin is not None and string(a, b) >= smin))
    ok.sec_r, ok.string = sec, string
    return ok


def afdb_heterodimers(st):
    """AFDB heterodimer models (file order) between proteins left in the monomer pool: gates min_pool and min_iptm;
    filters: no partner a larger homo-oligomer (UniProt n > 2 or a legacy oligomer) unless Complex Portal pairs them,
    not both in one modeled complex, native histone pairs only (data.histone_classes), not two members of
    data.monomeric_families (Hsp70), paralogues (data paralog_pairs or shared 3-mers >= paralog_k3) only with
    evidence (Complex Portal, IntAct), a shared compartment. Copies (allocation: proportional): each protein offers
    its copies to its accepted partners in proportion to their copies; a pair forms bound_fraction x the smaller
    offer (1.0 with Complex Portal / IntAct evidence). Writes heterodimers.csv
    origin: yeast pipeline/heterodimers.py:assign and scripts/s1_inventory.py:main @ e3b7448"""
    o, df, ing = st.o["afdb_heterodimers"], st.df, st.ing
    log = st.ctx.say
    builder = st.hook(o["builder"]) if o.get("builder") else None
    modeled_complexes = [{g.strip() for g in gs.split(";")} for gs in ing[ing.kind == "complex"].genes]
    legacy_oligo = st.legacy
    h = pd.read_csv(st.frozen("afdb_heterodimers"))
    ia = pd.read_csv(st.frozen("intact_pairs"), sep="\t")
    IA = {frozenset((a, b)): (p, d) for a, b, p, d in zip(ia.acc_a, ia.acc_b, ia.n_publications, ia.direct)}
    pf = o["paralog_filter"]
    SEQ = read_fasta(st.frozen("proteome_fasta")) if pf["needs_evidence"] else {}
    par = pd.read_csv(st.frozen("paralog_pairs"))
    PAR = {frozenset((a, b)) for a, b in zip(par.gene1_uniprot, par.gene2_uniprot)}
    cp = pd.read_csv(st.frozen("complex_portal"), sep="\t")
    CP = set()
    for s in cp.iloc[:, 4]:
        accs = re.findall(r"([A-Z0-9]{6,10})(?:-\w+)?\(\d+\)", s)
        CP |= {frozenset((a, b)) for a in accs for b in accs if a != b}
    HISTONES = dict(zip(*[read_curated(st.data["histone_classes"])[c] for c in ("gene", "histone_class")])) \
        if st.data.get("histone_classes") else {}
    HISTONE_PAIRS = {frozenset(p.split("+")) for p in o["native_histone_pairs"]}
    MONO = read_curated(st.data["monomeric_families"]) if st.data.get("monomeric_families") else pd.DataFrame(columns=["gene", "family"])
    mono_fam = dict(zip(MONO.gene, MONO.family))
    min_pool, min_iptm, max_n = o["min_pool"], o["min_iptm"], o["max_partner_oligomer_n"]
    olig = dict(zip(df.uniprot, df.oligomer_n)); gene_of = dict(zip(df.uniprot, df.gene))
    prot = ing[ing.kind == "protein"].set_index("uniprot")
    in_complex = {}
    for gs in modeled_complexes:
        for g in gs:
            in_complex.setdefault(g, []).append(frozenset(gs))
    legacy_genes = set(legacy_oligo)
    sc = o["shared_compartment"]
    ev_ok = _pair_evidence(st, o.get("evidence"), h)   # None = no gate (evidence: none, legacy_x5)
    rows, acc_ok = [], []
    for r in h.itertuples():
        a, b = r.uniprotAccession_1, r.uniprotAccession_2; k = frozenset((a, b))
        ga, gb = gene_of.get(a, r.gene_1), gene_of.get(b, r.gene_2)
        why = None
        ca = prot.copies_per_cell.get(a, 0) if a in prot.index else 0
        cb = prot.copies_per_cell.get(b, 0) if b in prot.index else 0
        pubs, direct = IA.get(k, (0, False))
        cpx = k in CP
        if min(ca, cb) < min_pool: why = "partner not in monomer pool"
        elif r.ipTM < min_iptm: why = f"ipTM < {min_iptm}"
        elif not cpx and (olig.get(a, 1) > max_n or olig.get(b, 1) > max_n or ga in legacy_genes or gb in legacy_genes): why = "partner is a larger homo-oligomer"
        elif set(in_complex.get(ga, [])) & set(in_complex.get(gb, [])): why = "both are subunits of a modeled complex"
        elif (ga in HISTONES or gb in HISTONES) and frozenset((HISTONES.get(ga, ga), HISTONES.get(gb, gb))) not in HISTONE_PAIRS:
            why = "non-native histone-fold pair"
        elif ga in mono_fam and gb in mono_fam and mono_fam[ga] == mono_fam[gb]:
            why = f"{mono_fam[ga]} pair ({mono_fam[ga]}s act as monomers)"
        elif pf["needs_evidence"] and not (cpx or pubs > 0 or direct) and (k in PAR or _k3(SEQ.get(a, ""), SEQ.get(b, "")) >= pf["k3"]):
            why = "paralogue pair without interaction evidence"
        if why is None and ev_ok is not None and not ev_ok(a, b):
            why = "no SEC co-elution / STRING evidence"
        w = None
        if why is None:
            w = _compatible(json.loads(prot.loc[a, "comp_weights"]), json.loads(prot.loc[b, "comp_weights"]),
                            prot.n_tm.get(a, 0) > 0, prot.n_tm.get(b, 0) > 0, sc["min_shared"], set(sc["facing"]),
                            sc["inner_membrane"], sc["inner_side"])
            if w is None: why = "no shared compartment"
        rows.append(dict(model=r.modelEntityId, a=ga, b=gb, acc_a=a, acc_b=b, ipTM=r.ipTM, ipSAE=r.max_ipSAE, publications=pubs,
                         direct=direct, complex_portal=cpx, paralog=k in PAR or _family(ga, gb), decision=why or "accepted", copies=0.0))
        if ev_ok is not None:
            rows[-1].update(sec_r=ev_ok.sec_r(a, b), string=ev_ok.string(a, b))
        if why is None:
            acc_ok.append((len(rows) - 1, w))
    if o["allocation"] != "proportional":
        raise ValueError(f"afdb_heterodimers.allocation '{o['allocation']}' is not available with pools: ingredients")
    partners = {}
    for i, _ in acc_ok:
        a_, b_ = rows[i]["acc_a"], rows[i]["acc_b"]
        partners.setdefault(a_, []).append(b_); partners.setdefault(b_, []).append(a_)
    cpc = prot.copies_per_cell
    offer = lambda x, y: cpc[x] * cpc[y] / _lsum(cpc[z] for z in partners[x])
    pool = prot.copies_per_cell.to_dict()
    bf = o["bound_fraction_by_evidence"]
    mass_by_acc = df.set_index("uniprot").mass_da
    new = []
    for i, w in acc_ok:
        r = rows[i]; a, b = r["acc_a"], r["acc_b"]
        frac = bf["complex_portal_or_intact"] if (r["complex_portal"] or r["publications"] > 0) else bf["other"]
        n = frac * min(offer(a, b), offer(b, a))
        tm = prot.n_tm.get(a, 0) > 0 or prot.n_tm.get(b, 0) > 0
        lib_id = re.sub(r"\W+", "_", f"{r['a']}-{r['b']}").replace("-", "_") + "_het"
        if builder is not None:
            built, info = builder.heterodimer(lib_id, r["model"], a, b, tm)
            if built is None:
                r["decision"] = f"membrane orientation failed ({info.get('fail')})"; continue
        r["copies"] = n
        pool[a] -= n; pool[b] -= n
        new.append(dict(ing_id=lib_id, label=f"{r['a']}-{r['b']} heterodimer", kind="complex", genes=f"{r['a']};{r['b']}",
                        n_subunit_types=2, copies_per_cell=n, comp_weights=json.dumps(w),
                        mass_da=float(mass_by_acc.loc[[a, b]].sum()), pdb="", membrane="",
                        oligomer_n=1, uniprot="", n_tm=int(prot.n_tm.get(a, 0) + prot.n_tm.get(b, 0))))
    ing = ing.copy()
    m = ing.kind == "protein"
    ing.loc[m, "copies_per_cell"] = [max(0.0, pool.get(u, c)) for u, c in zip(ing.loc[m, "uniprot"], ing.loc[m, "copies_per_cell"])]
    ing = pd.concat([ing, pd.DataFrame(new)], ignore_index=True)
    st.ing = ing[(ing.kind != "protein") | (ing.copies_per_cell >= 1)].reset_index(drop=True)
    t = pd.DataFrame(rows)
    log(f"heterodimers: {len(new)} accepted; " + str(t.decision.str.replace(r" \(.*", "", regex=True).value_counts().to_dict()))
    st.write_csv(t, "heterodimers.csv")


# ------------------------------------------------------------------------------------------------ homo-oligomers
def _homo_share_rule(st, o):
    """SEC-fed share of a homo-oligomer (as the inventory pools, steps.py pdb_homo_oligomers / afdb_homodimers, Jurkat
    legacy_v7; rules in wholecell/inventory/secms.py): option share {source: secms | constant, k_by_state, default,
    value, apex_ratio, apex_check_min_k, apex_check_min_frac} for PDB / audit states and af_share {source, k, default,
    value} for the AF homodimer fallback. The oligomer takes share x the ingredient's copies as a new ingredient row
    '<ing_id>__homo' (the parent's comp_weights, uniprot, ...; oligo_map -> the oligomer library entry); the parent row
    keeps the rest as monomers. SEC apex check (PDB states): a protein SEC-assembled at apex_check_min_k (a_k >
    apex_check_min_frac) whose SEC apex is outside apex_ratio x (state x monomer mass) stays a monomer
    (apex_check_side: both, as legacy_v7); apex_check_side: below (x5_secms, addendum 2026-11-19) demotes only an apex
    BELOW the window - an apex above it (a larger assembly, or beyond the top calibrant) is not evidence for a monomer:
    the oligomer keeps its SEC share. No share option (legacy_x5): None."""
    if o["copies"] != "share":
        return None
    from .. import secms as S2
    lo = {"homo_oligomer_share": o["share"], "afdb_homodimer_share": o.get("af_share") or
          {"source": "constant", "value": 1.0}}
    view = st.view
    if view is not None:
        sa = view.legacy_sec_assembly(int(o["share"].get("apex_k", 2))).set_index("uniprot")
    else:
        sa = pd.DataFrame(columns=["apex_app_kda", "mono_kda"])
    split = []
    side = o["share"].get("apex_check_side", "both")

    def f(r, oligo, row, from_af):
        lib_id, n = oligo[r.ing_id]
        if from_af:
            share, k = S2.homodimer_share(view, r.uniprot, lo), lo["afdb_homodimer_share"].get("k")
        else:
            a4 = view.a(r.uniprot, o["share"]["apex_check_min_k"]) if view is not None else np.nan
            if not S2.homo_apex_ok(r.uniprot, n, sa, a4, lo):
                above = sa.apex_app_kda[r.uniprot] / (n * sa.mono_kda[r.uniprot]) > o["share"]["apex_ratio"][1]
                if side == "below" and above:
                    row.update(decision=row["decision"] + "; SEC apex above state: kept")
                else:
                    del oligo[r.ing_id]
                    row.update(decision=row["decision"] + "; SEC apex check: monomer", lib_id="", share=0.0)
                    return
            share, k = S2.homo_share(view, r.uniprot, n, lo)
        if share != share:                                # SEC row without a value at this k: the rule's default
            share = float((lo["afdb_homodimer_share"] if from_af else lo["homo_oligomer_share"])["default"])
        share = float(min(max(share, 0.0), 1.0))
        row.update(share=share, sec_k=k, sec_row=bool(view is not None and view.has(r.uniprot)))
        split.append((r.ing_id, share))

    def apply(oligo):
        ing = st.ing.set_index("ing_id", drop=False)
        new = []
        for ing_id, share in split:
            lib_id, n = oligo.pop(ing_id)
            c = float(ing.at[ing_id, "copies_per_cell"])
            if share * c < 1:                                 # nothing assembled: monomer
                continue
            nid = f"{ing_id}__homo"
            x = ing.loc[ing_id].copy()
            x["ing_id"], x["copies_per_cell"] = nid, share * c
            x["label"] = f"{x['label']} homo-oligomer" if "label" in x.index and isinstance(x["label"], str) else nid
            new.append(x)
            ing.at[ing_id, "copies_per_cell"] = c - share * c
            oligo[nid] = (lib_id, n)
        ing = pd.concat([ing.reset_index(drop=True), pd.DataFrame(new)], ignore_index=True) if new else ing.reset_index(drop=True)
        st.ing = ing[(ing.kind != "protein") | (ing.copies_per_cell >= 1) | ing.ing_id.isin(list(oligo))].reset_index(drop=True)
    f.apply = apply
    return f


def pdb_homo_oligomers(st):
    """homo-oligomer structures of the residual monomer pool (copies: none = only the ingredient -> (structure, n) map,
    no copies move; copies: share = the SEC-fed share of _homo_share_rule becomes the oligomer, the rest monomers): proteins with >= min_copies, by copies descending; state from data.homo_oligomer_audit (PDB
    audit) or the UniProt subunit annotation (oligomer_n > 1: the best-resolution single-entity PDB entry whose
    assembly1 has n instances and >= 0.6 n x length modeled residues; frozen RCSB cache), else an AFDB homodimer
    (annotated: ipTM >= annotated_min_iptm; de novo: ipTM >= de_novo_min_iptm and ipSAE >= de_novo_min_ipsae; or, with
    ipsae_only_min set (x5_secms), the legacy_v7 gate max ipSAE > ipsae_only_min).
    Legacy oligomers are kept. Writes oligomers.csv, oligo_map.json (at the end of the step list)
    origin: yeast pipeline/assemblies.py:homo_oligomers and scripts/s1_inventory.py:main @ e3b7448"""
    o, df, ing = st.o["pdb_homo_oligomers"], st.df, st.ing
    log = st.ctx.say
    if o["copies"] not in ("none", "share"):
        raise ValueError(f"pdb_homo_oligomers.copies '{o['copies']}' is not available with pools: ingredients (allowed: none, share)")
    builder = st.hook(o["builder"]) if o.get("builder") else None
    R = st.rcsb()
    audit = {}
    if st.data.get("homo_oligomer_audit"):
        A = read_curated(st.data["homo_oligomer_audit"])
        audit = {g: (int(n), p) for g, n, p in zip(A.gene, A.chains, A.pdb)}
    afd = {r.uniprotAccession_1: r for r in pd.read_csv(st.frozen("afdb_homodimers")).itertuples()}
    fb = o["afdb_fallback"]
    sh = _homo_share_rule(st, o)                     # None = legacy_x5 (the whole ingredient is the oligomer)
    legacy = st.legacy
    oligo, rows = dict(legacy), []
    cand = ing[(ing.kind == "protein") & (ing.copies_per_cell >= o["min_copies"]) & ~ing.ing_id.isin(list(legacy))]
    for r in cand.sort_values("copies_per_cell", ascending=False).itertuples():
        n_uni = int(r.oligomer_n); lib_id = f"{r.ing_id}_oligo"; how = None; tm = r.n_tm > 0
        if n_uni > 1 or r.ing_id in audit:
            if r.ing_id in audit:
                n_uni, pid = audit[r.ing_id]
            else:
                pid = R.get(f"homo:{r.uniprot}:{n_uni}")
            if pid:
                if builder is not None:
                    built, info = builder.pdb_homo(lib_id, pid, n_uni, tm, only_acc=r.uniprot if r.ing_id in audit else None)
                else:
                    built, info = 0, {}
                if built is not None:
                    oligo[r.ing_id] = (lib_id, n_uni); how = f"PDB {pid} ({n_uni}-mer){' membrane' if tm else ''}"
                else:
                    how = f"PDB {pid} membrane orientation failed ({info.get('fail')})"
        if r.ing_id not in oligo and r.uniprot in afd:
            d = afd[r.uniprot]
            ok = d.ipTM >= fb["annotated_min_iptm"] if n_uni > 1 else (d.ipTM >= fb["de_novo_min_iptm"] and d.max_ipSAE >= fb["de_novo_min_ipsae"])
            if not ok and fb.get("ipsae_only_min") is not None:     # legacy_v7 afdb_homodimers gate (x5_secms addendum)
                ok = d.max_ipSAE > fb["ipsae_only_min"]
            if ok:
                if builder is not None:
                    built, info = builder.af_homodimer(lib_id, d.modelEntityId, r.uniprot, tm)
                else:
                    built, info = 0, {}
                if built is not None:
                    oligo[r.ing_id] = (lib_id, 2); how = f"AF dimer {d.modelEntityId} (ipTM {d.ipTM:.2f}){' membrane' if tm else ''}"
                else:
                    how = f"AF dimer membrane orientation failed ({info.get('fail')})"
            else:
                how = f"AF dimer rejected (ipTM {d.ipTM:.2f}, ipSAE {d.max_ipSAE:.2f})"
        row = dict(ing_id=r.ing_id, uniprot=r.uniprot, copies=r.copies_per_cell, uniprot_oligomer=n_uni,
                   decision=how or "monomer", lib_id=oligo.get(r.ing_id, ("", 0))[0])
        if sh is not None and r.ing_id in oligo:
            sh(r, oligo, row, from_af=how.startswith("AF dimer"))
        rows.append(row)
    if sh is not None:
        sh.apply(oligo)
    otab = pd.DataFrame(rows)
    st.write_csv(otab, "oligomers.csv")
    log("homo-oligomer decisions:", otab.decision.str.replace(r"\([^)]*\)|AF-\S+|\b[0-9][0-9A-Z]{3}\b", "", regex=True)
        .str.replace(r"\s+", " ", regex=True).str.strip().value_counts().to_dict())
    st.oligo_map = oligo


def finish(st):
    """after the step list: ingredients.csv (the species the packer reads) and oligo_map.json
    origin: yeast scripts/s1_inventory.py:main @ e3b7448"""
    om = st.oligo_map if st.oligo_map is not None else dict(st.legacy)
    st.ctx.say(f"mass in assemblies {_mass_in_assemblies(st.ing, om):.3f}")
    st.write_csv(st.ing, "ingredients.csv")
    st.write_json({k: list(v) for k, v in om.items()}, "oligo_map.json", indent=0)


STEPS = {"curated_machines": curated_machines, "complex_portal": complex_portal,
         "afdb_heterodimers": afdb_heterodimers, "pdb_homo_oligomers": pdb_homo_oligomers}
