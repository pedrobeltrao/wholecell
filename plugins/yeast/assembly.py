"""Yeast plugin parts of the core assembly step list (wholecell/inventory/assembly, preset legacy_x5; unification A4).

Registered here (kind `assembly_step`: a step the preset lists by name; kind `assembly_hook`: a function a core step
calls through an option value):

  assembly_step  yeast_barrel_monomers   outer-membrane beta-barrel monomers as oriented membrane entries (A9)
  assembly_step  yeast_preribosomes      pre-ribosome stages 90S / pre-60S x3 / pre-40S (engine pipeline/preribosomes.py)
  assembly_step  yeast_drop_free_rproteins  the free ribosomal-protein pool above the ribosome count is not packed
  assembly_step  yeast_er_ribosomes      80S on Ssh1 / Sec61 translocons (engine pipeline/er_ribosomes.py)
  assembly_step  yeast_mitoribosomes     74S mitoribosomes on the inner membrane (engine pipeline/mitoribosomes.py)
  assembly_step  yeast_supercomplexes    respiratory supercomplex III2IV2 (engine pipeline/supercomplexes.py)
  assembly_step  yeast_arp23             Arp2/3 as one 7-subunit complex (engine pipeline/actin_patches.add_arp23)
  assembly_step  yeast_f_actin           F-actin reservation (cables + patches) before the heterodimers
  assembly_hook  yeast_nucleosome_count  nucleosomes on the chromosome arms (chromatin map) + nucleosomal rDNA
  assembly_hook  yeast_npc_count         nuclear pores per cell (engine inventory.N_NPC, areal density x NE area)
  assembly_hook  yeast_structures        structure builder: beads, membrane orientation, library entry (engine
                                         assemblies.build_complex / build_pdb_homo / build_af_dimer, heterodimers.build)
  inventory_builder yeast.x5_ingredients the assembly ingredient table in the core pack-table layout (A9)

These steps are genuinely yeast biology (or the yeast structure-library builder, whose bead model and membrane
orientation belong to the structures stage) and call the vendored engine modules (plugins/yeast/engine, never edited)
until they move (unification A9 / structures stage). The engine modules read their paths from YC_* variables at import:
engine_env() sets them from params.assembly_steps.engine {root, work, lib, pdb, code, inputs} before the first import.
"""
from __future__ import annotations

import json
import os
import sys

from wholecell.registry import register

_ENV = None


def engine_env(E: dict):
    """point the vendored engine at a layout (must run before any engine module is imported; idempotent)"""
    global _ENV
    if _ENV is not None:
        if _ENV != E:
            raise RuntimeError(f"yeast engine already configured with {_ENV}; cannot switch to {E} in one process")
        return
    if os.environ.get("PYTHONHASHSEED") != "0":
        raise RuntimeError("the yeast engine steps need PYTHONHASHSEED=0 (set iteration order of the reference build)")
    for k in ("root", "work", "lib", "pdb", "code"):
        if not E.get(k):
            raise RuntimeError(f"inventory.params.assembly_steps.engine.{k} is required by the yeast plugin steps")
    os.environ.update(YC_ROOT=E["root"], YC_WORK=E["work"], YC_LIB=E["lib"], YC_PDB=E["pdb"], YC_SEED=str(int(E.get("seed", 0))),
                      YC_NPROC=str(E.get("nproc", 3)), OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    if E.get("inputs"):
        os.environ["YC_INPUTS"] = E["inputs"]
    os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(E["work"], "numba_cache"))   # never next to the vendored code
    p = os.path.join(E["code"], "pipeline")
    if p not in sys.path:
        sys.path.insert(0, p)
    _ENV = dict(E)


def resolve_layout(E: dict, bdir: str) -> dict:
    """engine layout values "build:<rel>" -> <bdir>/<rel> (config/cells/yeast.yaml: the build's engine/ folder)"""
    return {k: (os.path.join(bdir, v.split(":", 1)[1]) if isinstance(v, str) and v.startswith("build:") else v)
            for k, v in (E or {}).items()}


def _env(st):
    engine_env(resolve_layout(st.params.get("engine") or {}, st.ctx.bdir))


def _log(st):
    return lambda *a: st.ctx.say(" ".join(str(x) for x in a))


# ------------------------------------------------------------------ plugin steps (ingredient table in, out)
def barrel_monomers(st):
    """origin: yeast scripts/s1_inventory.py:main @ e3b7448 (first lines): log the absent / ribosome-density inventory
    numbers and rebuild the AlphaFold monomers of the outer-membrane beta-barrels (MDM10, MDM34, SAM50, TOM40) as
    oriented membrane entries over their soluble library entries (engine assemblies.build_barrel_monomers; writes
    <engine lib>/<gene>.npz, read by the pack). Found missing in the A9 stage check (L2: 4 entries)."""
    _env(st)
    import assemblies as A
    import inventory as inv
    say = _log(st)
    say("not in this cell (genotype, mating type, meiosis, glucose repression, G1; copies set to zero):",
        st.df.attrs.get("not_in_this_cell"))
    if not inv.RIBO_SCALE:                      # filled by inventory.build(); the protein table step ran it in another process
        inv.build()
    rs = inv.RIBO_SCALE
    say(f"ribosomes from the cryo-ET density: {rs['ribosomes']:,.0f} (PaxDb {rs['before']:,.0f}; ribosomal proteins x "
        f"{rs['factor']:.3f}); protein total {rs['protein_pg']:.2f} pg")
    A.build_barrel_monomers(st.df, log=say)


def preribosomes(st):
    _env(st)
    import preribosomes as prb
    st.ing, _ = prb.assign(st.ing, st.df, log=_log(st))


def drop_free_rproteins(st):
    """origin: yeast scripts/s1_inventory.py:main @ e3b7448 (the excess free r-protein pool, MS noise above the
    ribosome count, is not packed)"""
    ing = st.ing
    st.ing = ing[~((ing.kind == "protein") & ing.ing_id.str.match(r"^RP[LS]\d"))].reset_index(drop=True)


def er_ribosomes(st):
    _env(st)
    import er_ribosomes as er
    st.ing, _ = er.assign(st.ing, st.df, log=_log(st))


def mitoribosomes(st):
    _env(st)
    import mitoribosomes as mr
    st.ing, _ = mr.assign(st.ing, st.df, log=_log(st))


def supercomplexes(st):
    _env(st)
    import supercomplexes as sc
    st.ing = sc.add(st.ing, st.df, log=_log(st))


def arp23(st):
    _env(st)
    import actin_patches as ap
    st.ing = ap.add_arp23(st.ing, st.df, log=_log(st))


def f_actin(st):
    """origin: yeast scripts/s1_inventory.py:main @ e3b7448. F-actin (cables + patches) is reserved before the
    heterodimer step, which then pairs only monomeric actin; writes f_actin.json"""
    _env(st)
    import actin_patches as ap
    import cytoskeleton as cs
    df, ing2 = st.df, st.ing
    total = float(df.set_index("gene").loc["ACT1", "copies_per_cell"])
    r = (ing2.ing_id == "ACT1") & (ing2.kind == "protein"); free = float(ing2.loc[r, "copies_per_cell"].sum())
    want = {"cables": cs.CABLE_FRACTION * total, "patches": ap.PATCH_FRACTION * total}
    f = min(1.0, free / sum(want.values())); res = {k: v * f for k, v in want.items()}
    ing2.loc[r, "copies_per_cell"] = free - sum(res.values())
    st.write_json(dict(res, total_act1=total), "f_actin.json", indent=1)
    st.ctx.say(f"F-actin reserved: cables {res['cables']:.0f}, patches {res['patches']:.0f} of {total:.0f} ACT1; "
               f"monomer pool {free - sum(res.values()):.0f}")


# ------------------------------------------------------------------ hooks
def nucleosome_count(st):
    """nucleosomes on the chromosome arms (measured dyads, nucleosome_map.arm_count) + the nucleosomal rDNA"""
    _env(st)
    import nucleosome_map
    import rdna
    return nucleosome_map.arm_count() + rdna.nucleosomes_in_rdna()


def npc_count(st):
    _env(st)
    import inventory as inv
    return inv.N_NPC


class Structures:
    """structure builder hook: build the library entry (beads, own RNA, membrane orientation); returns what the
    vendored builders return ((n_beads, n_nt, n_rna) / n_beads, info) or (None, info) when a membrane assembly
    cannot be oriented"""

    def __init__(self, st):
        _env(st)
        import assemblies as A
        import heterodimers as H
        self.A, self.H = A, H

    def complex(self, ing_id, cname, pid, structure, psub, acc_of, keep_rna, tm):
        return self.A.build_complex(ing_id, cname, pid, structure, psub, acc_of, keep_rna, tm)

    def heterodimer(self, lib_id, model_id, acc_a, acc_b, tm):
        return self.H.build(lib_id, model_id, acc_a, acc_b, tm)

    def pdb_homo(self, lib_id, pid, n, tm, only_acc=None):
        return self.A.build_pdb_homo(lib_id, pid, n, tm, only_acc=only_acc)

    def af_homodimer(self, lib_id, model_id, acc, tm):
        return self.A.build_af_dimer(lib_id, model_id, acc, tm)

    def has_entry(self, lib_id):
        return os.path.exists(os.path.join(os.environ["YC_LIB"], f"{lib_id}.npz"))

    def fetch_assembly(self, pid):
        import extend_complexes as xc
        return xc.fetch_assembly(pid)

    def save_cache(self):
        self.A.save_cache()


for _n, _f in (("yeast_barrel_monomers", barrel_monomers), ("yeast_preribosomes", preribosomes), ("yeast_drop_free_rproteins", drop_free_rproteins),
               ("yeast_er_ribosomes", er_ribosomes), ("yeast_mitoribosomes", mitoribosomes),
               ("yeast_supercomplexes", supercomplexes), ("yeast_arp23", arp23), ("yeast_f_actin", f_actin)):
    register("assembly_step", _n)(_f)
f_actin.outputs = ["f_actin.json"]
register("assembly_hook", "yeast_nucleosome_count")(nucleosome_count)
register("assembly_hook", "yeast_npc_count")(npc_count)
register("assembly_hook", "yeast_structures")(Structures)


# ------------------------------------------------------------------ pack table (inventory_builder)
def x5_ingredients(ctx):
    """the assembly ingredient table (assembly/ingredients.csv: one row per packed species, fractional copies) in the
    core pack-table layout (schema.PACK_TABLE), row order kept. The yeast pack reads the ingredient table itself
    (pre-placed context: wc3d.expected_totals); this table is what the canonical tables, the suite and the export read.
    copies stay fractional (the yeast method); struct_src = the PDB entry or afdb_monomer:<UniProt>."""
    import numpy as np
    import pandas as pd
    from wholecell.inventory import schema
    I = pd.read_csv(ctx.need("assembly/ingredients.csv"))
    ntm = pd.to_numeric(I.get("n_tm"), errors="coerce").fillna(0).astype(int) if "n_tm" in I.columns else pd.Series(0, I.index)
    memb = I["membrane"].fillna("").astype(str) if "membrane" in I.columns else pd.Series("", I.index)
    pdb = I["pdb"].fillna("").astype(str) if "pdb" in I.columns else pd.Series("", I.index)
    uni = I["uniprot"].fillna("").astype(str) if "uniprot" in I.columns else pd.Series("", I.index)
    out = pd.DataFrame(dict(
        ing_id=I.ing_id, kind=I.kind, genes=I.genes.fillna(I.ing_id).astype(str),
        stoich=I.get("n_subunit_types", pd.Series("", I.index)).fillna("").astype(str),
        copies=I.copies_per_cell.astype(float), mass_da=I.mass_da.astype(float), comp_weights=I.comp_weights.astype(str),
        transmembrane=(ntm > 0) | (memb != ""), n_tm=ntm,
        struct_src=np.where(pdb != "", "pdb:" + pdb, "afdb_monomer:" + uni), uniprot=uni,
        signal_frac=0.0, mito_om_frac=0.0))
    out["mass_pg"] = out.copies * out.mass_da / schema.AVOGADRO * 1e12
    ctx.say(f"yeast.x5_ingredients: {len(out)} species, {out.copies.sum():,.0f} copies, {out.mass_pg.sum():.4f} pg")
    return out[schema.PACK_TABLE]


x5_ingredients.PREV = ["assembly/ingredients.csv"]
register("inventory_builder", "yeast.x5_ingredients")(x5_ingredients)
