"""Inventory stage: abundance -> localisation -> assembly state -> species (ingredients) -> rule sets -> pack table + canonical.

    python -m wholecell inventory --cell <cell> --build <bdir> [--step S]      (no step: every step whose outputs are missing)

Steps (in order; the list is built from the config, see `steps(cfg)`):
  rna           adapter.rna (inventory.adapters.rna; skipped when null)       -> inventory/rna/rna.tsv
  copies        generic copy-number step (only when inventory.copies is set; copies.py): typed abundance sources,
                priority/provenance, imputation, absent lists, per-protein scaling, total protein (fixed | anchors |
                source_sum), structure-count adjustments, anchor check, copy ledger -> inventory/copies/*; the
                proteome adapter `copies` hands its table to v0
  v0            adapter.proteome + adapter.localisation, joined per gene       -> inventory/v0/inventory_v0.tsv
                (core localisation = evidence table + combiner + preset, wholecell/inventory/localisation/
                                                                                -> inventory/localisation/*)
                (total protein: inventory.total_protein_pg, REQUIRED, no default; its source: total_protein_source)
  secms         adapter.secms (only when set; generic: secms_toolkit = toolkits/secms) -> inventory/secms/* (SEC-MS
                tables, docs/reference/inventory_schema.md §4a); independent of v0
  assembly      adapter.assembly (complexes, SEC-MS assembled fractions, oligomers, dimers) -> inventory/assembly/*
                (core: assembly_steps = ordered step list + preset, wholecell/inventory/assembly/)
  ingredients   inventory_builder (inventory.ingredients): per-gene inventory + assemblies -> packing species
                                                                                -> inventory/ingredients/ingredients.tsv
  rules:<name>  one step per entry of inventory.rule_sets (registered `inventory_rule`), each reading the previous
                step's ingredients.tsv -> inventory/rules/<name>/{ingredients.tsv, changes.tsv[, changes_legacy.tsv], report.txt}
  pack_table    the last step's ingredients.tsv (byte copy; with inventory.pack_table.order_reference = a frozen legacy
                table, rows / equal-content comp_weights key order follow the reference, see legacy_order.py)
                                                                                -> inventory/pack_table.tsv [, pack_table_order.json]
  canonical     canonical tables of docs/reference/inventory_schema.md             -> inventory/{species, species_loc,
                field_provenance, sources, proteome_sources, localisation_sources, gene_abundance, secms_protein,
                secms_complex, surface, changes_all, version_changes, secms_bound_fraction}.tsv
  compare       optional (inventory.verify present): port vs legacy reference  -> inventory/compare_legacy.json

Jurkat: rules:jurkat.v7 reads the lamina pre-pass consumption from params."jurkat.v7".lamina_consumption (frozen name
inv_lamina_consumption_v7 = legacy file, or "build:<rel path>" = the lamina rebuilt by the preplaced stage of this build).

`STEPS` documents the generic steps; `steps(cfg)` returns the concrete ordered dict
  {step: {"inputs": [paths], "frozen": [frozen input names], "outputs": [paths]}}
with paths relative to the build dir (Snakemake: inputs + cfg.frozen(n) for n in frozen). Registered components declare
their extra inputs with attributes FROZEN (frozen names), PREV (inventory-relative paths of earlier outputs),
PARAM_FILES (keys of their params block that hold file paths) and OUTPUTS (inventory-relative paths they write; or a
function outputs(cfg) when they depend on the config).
"""
from __future__ import annotations

import os
from collections import OrderedDict

from .. import build, registry
from . import schema
from .context import StepContext

STEPS = OrderedDict([
    ("rna", {"inputs": [], "outputs": ["inventory/rna/rna.tsv"]}),
    ("copies", {"inputs": ["inventory/rna/rna.tsv (if RNA imputation)"], "outputs": ["inventory/copies/*"]}),
    ("v0", {"inputs": ["inventory/rna/rna.tsv (if rna adapter)"], "outputs": ["inventory/v0/inventory_v0.tsv"]}),
    ("secms", {"inputs": [], "outputs": ["inventory/secms/* (adapter OUTPUTS)"]}),
    ("assembly", {"inputs": ["inventory/v0/inventory_v0.tsv", "inventory/secms/*"], "outputs": ["inventory/assembly/* (adapter OUTPUTS)"]}),
    ("ingredients", {"inputs": ["inventory/assembly/*"], "outputs": ["inventory/ingredients/ingredients.tsv"]}),
    ("rules:<name>", {"inputs": ["previous step ingredients.tsv"],
                      "outputs": ["inventory/rules/<name>/ingredients.tsv", "inventory/rules/<name>/changes.tsv"]}),
    ("pack_table", {"inputs": ["last ingredients.tsv"], "outputs": ["inventory/pack_table.tsv"]}),
    ("canonical", {"inputs": ["inventory/pack_table.tsv", "all changes.tsv"], "outputs": ["inventory/species.tsv", "..."]}),
    ("compare", {"inputs": ["all outputs"], "outputs": ["inventory/compare_legacy.json"]}),
])

CANONICAL_OUT = ["species.tsv", "species_loc.tsv", "field_provenance.tsv", "sources.tsv", "proteome_sources.tsv",
                 "localisation_sources.tsv", "gene_abundance.tsv", "secms_protein.tsv", "secms_complex.tsv",
                 "secms_complex_subunit.tsv", "secms_pairs.tsv", "secms_calibration.tsv",
                 "surface.tsv", "secms_bound_fraction.tsv", "secms_bound_fraction_summary.json", "changes_all.tsv", "version_changes.tsv"]


def _inv(cfg):
    return cfg.get("inventory", {}) or {}


def _comp(kind, name):
    return registry.get(kind, name)


def _scan_frozen(obj, names, out):
    """every string in a params block that names a declared frozen input"""
    if isinstance(obj, str):
        if obj in names:
            out.add(obj)
    elif isinstance(obj, dict):
        for v in obj.values():
            _scan_frozen(v, names, out)
    elif isinstance(obj, list):
        for v in obj:
            _scan_frozen(v, names, out)
    return out


def _attr_inputs(obj, cfg, pname):
    """declared inputs of a registered component: (frozen names, build-relative paths, other file paths). Frozen names =
    attribute FROZEN + every params string naming a frozen input; PARAM_FILES values of the form 'build:<rel>' are
    build-relative inputs (e.g. the rebuilt lamina of the preplaced stage)."""
    names = set((cfg.get("frozen_inputs") or {}).keys())
    params = ((_inv(cfg).get("params", {}) or {}).get(pname, {}) or {})
    fr = set(getattr(obj, "FROZEN", []) or []) | _scan_frozen(params, names, set())
    prev = ["inventory/" + p for p in (getattr(obj, "PREV", []) or [])]
    pf = []
    for k in getattr(obj, "PARAM_FILES", []) or []:
        v = params.get(k)
        if not isinstance(v, str) or v in names:
            continue
        pf.append(v[len("build:"):] if v.startswith("build:") else v)
    return sorted(fr), prev, pf


def rule_dir(name: str) -> str:
    return f"rules/{name}"


def _register_core():
    from . import assembly, copies, generic, localisation, secms  # noqa: F401  (core adapters/builders/rules/readers register themselves)


def steps(cfg) -> "OrderedDict[str, dict]":
    _register_core()
    inv = _inv(cfg)
    ad = inv.get("adapters", {}) or {}
    out = OrderedDict()

    def add(name, obj, pname, inputs, outputs):
        fr, prev, pf = _attr_inputs(obj, cfg, pname) if obj is not None else ([], [], [])
        out[name] = {"inputs": inputs + prev + pf, "frozen": fr, "outputs": outputs}

    if ad.get("rna"):
        add("rna", _comp("adapter.rna", ad["rna"]), ad["rna"], [], ["inventory/rna/rna.tsv"])
    if inv.get("copies"):
        from . import copies as CP
        out["copies"] = {"inputs": ["inventory/rna/rna.tsv"] if CP.needs_rna(cfg) else [], "frozen": CP.frozen_names(cfg),
                         "outputs": ["inventory/" + p for p in CP.OUTPUTS]}
    fr, prev, pf, extra = [], [], [], []
    for k in ("proteome", "localisation"):
        obj = _comp(f"adapter.{k}", ad[k])
        a, b, c = _attr_inputs(obj, cfg, ad[k])
        fr += a; prev += b; pf += c
        extra += ["inventory/" + p for p in getattr(obj, "OUTPUTS", []) or []]   # e.g. localisation/evidence.tsv
    out["v0"] = {"inputs": (["inventory/rna/rna.tsv"] if ad.get("rna") else []) + prev + pf, "frozen": sorted(set(fr)),
                 "outputs": ["inventory/v0/inventory_v0.tsv"] + extra}
    sec_outs = []
    if ad.get("secms"):
        sm = _comp("adapter.secms", ad["secms"])
        sec_outs = ["inventory/" + p for p in getattr(sm, "OUTPUTS", [])]
        add("secms", sm, ad["secms"], [], sec_outs)
    asm = _comp("adapter.assembly", ad.get("assembly") or "none")
    asm_outs = asm.outputs(cfg) if hasattr(asm, "outputs") else getattr(asm, "OUTPUTS", [])   # outputs(cfg): per config
    add("assembly", asm, ad.get("assembly") or "none", ["inventory/v0/inventory_v0.tsv"] + sec_outs,
        ["inventory/" + p for p in asm_outs])
    out["assembly"]["inputs"] = list(dict.fromkeys(out["assembly"]["inputs"]))
    bname = inv.get("ingredients", "generic.monomers")
    add("ingredients", _comp("inventory_builder", bname), bname, out["assembly"]["outputs"],
        ["inventory/ingredients/ingredients.tsv"])
    prev_tab = "inventory/ingredients/ingredients.tsv"
    for name in inv.get("rule_sets", []) or []:
        r = _comp("inventory_rule", name)
        d = f"inventory/{rule_dir(name)}"
        outs = [f"{d}/ingredients.tsv", f"{d}/changes.tsv", f"{d}/report.txt"]
        outs += [f"{d}/{x}" for x in (getattr(r, "OUTPUTS", []) or [])]
        add(f"rules:{name}", r, name, [prev_tab], outs)
        prev_tab = f"{d}/ingredients.tsv"
    ptref = (inv.get("pack_table", {}) or {}).get("order_reference")
    out["pack_table"] = {"inputs": [prev_tab], "frozen": [ptref] if ptref else [],
                         "outputs": ["inventory/pack_table.tsv"] + (["inventory/pack_table_order.json"] if ptref else [])}
    rule_changes = [f"inventory/{rule_dir(n)}/changes.tsv" for n in inv.get("rule_sets", []) or []]
    out["canonical"] = {"inputs": ["inventory/pack_table.tsv", "inventory/v0/inventory_v0.tsv",
                                   "inventory/ingredients/ingredients.tsv"] + rule_changes + sec_outs,
                        "frozen": sorted(n for n in (cfg.get("frozen_inputs") or {}) if n.startswith("inv_")),
                        "outputs": [f"inventory/{x}" for x in CANONICAL_OUT]}
    if inv.get("verify"):
        out["compare"] = {"inputs": out["canonical"]["outputs"], "frozen": [],
                          "outputs": ["inventory/compare_legacy.json"]}
    return out


def _last_table(cfg) -> str:
    rs = _inv(cfg).get("rule_sets", []) or []
    return f"{rule_dir(rs[-1])}/ingredients.tsv" if rs else "ingredients/ingredients.tsv"


# ----------------------------------------------------------------------------------------------- step implementations
def step_rna(ctx):
    name = ctx.inv["adapters"]["rna"]
    df = _comp("adapter.rna", name)(ctx)
    ctx.write_tsv(df, "rna/rna.tsv", index=True)


def step_v0(ctx):
    ad = ctx.inv["adapters"]
    if "total_protein_pg" not in ctx.inv or ctx.inv["total_protein_pg"] in (None, ""):
        raise ValueError("inventory.total_protein_pg is required (no default; legacy default 80 pg was a defect)")
    if ctx.inv["total_protein_pg"] == "anchors":     # set by the copies step (inventory.copies.total_protein.mode anchors)
        if ad["proteome"] != "copies" or ((ctx.inv.get("copies") or {}).get("total_protein") or {}).get("mode") != "anchors":
            raise ValueError("total_protein_pg: anchors needs adapters.proteome: copies with copies.total_protein.mode anchors")
        tot = None
    else:
        tot = float(ctx.inv["total_protein_pg"])
    prot = _comp("adapter.proteome", ad["proteome"])(ctx, tot)
    loc = _comp("adapter.localisation", ad["localisation"])(ctx, prot)
    schema.check(prot.reset_index(), schema.PROTEOME_OUT, "adapter.proteome")
    schema.check(loc.reset_index(), schema.LOCALISATION_OUT, "adapter.localisation")
    inv = prot.copy()
    for c in ("compartments", "loc_source", "transmembrane"):
        inv[c] = loc[c].reindex(inv.index).values if c != "loc_source" else loc[c].reindex(inv.index)
    inv["mass_pg"] = inv.copies * inv.mass_da / schema.AVOGADRO * 1e12
    for c in loc.columns:                       # extra per-gene annotations of the adapter (signal peptide, OM, ...)
        if c not in inv.columns and c != "gene":
            inv[c] = loc[c].reindex(inv.index).values
    inv.index.name = "gene"
    ctx.write_tsv(inv, "v0/inventory_v0.tsv", index=True)


def step_copies(ctx):
    from . import copies
    copies.step(ctx)


def step_secms(ctx):
    _comp("adapter.secms", ctx.inv["adapters"]["secms"])(ctx)


def step_assembly(ctx):
    name = ctx.inv.get("adapters", {}).get("assembly") or "none"
    _comp("adapter.assembly", name)(ctx)


def step_ingredients(ctx):
    name = ctx.inv.get("ingredients", "generic.monomers")
    I = _comp("inventory_builder", name)(ctx)
    schema.check(I, schema.PACK_TABLE, f"inventory_builder {name}")
    ctx.write_tsv(I, "ingredients/ingredients.tsv")


def step_rule(ctx, name: str, src_rel: str):
    import pandas as pd
    fn = _comp("inventory_rule", name)
    d = rule_dir(name)
    res = fn(ctx, ctx.need(src_rel))
    ctx.write_tsv(res["ingredients"], f"{d}/ingredients.tsv")
    ch = res.get("changes")
    if ch is None:
        ch = schema.empty(schema.CHANGES)
    ch = ch.copy()
    if "step" not in ch.columns:
        ch.insert(0, "step", f"rules:{name}")
    if "rule_set" not in ch.columns:
        ch.insert(1, "rule_set", name)
    schema.check(ch, schema.CHANGES, f"{name} changes")
    ctx.write_tsv(ch[schema.CHANGES + [c for c in ch.columns if c not in schema.CHANGES]], f"{d}/changes.tsv")
    for rel, df in (res.get("extra") or {}).items():
        if isinstance(df, pd.DataFrame):
            ctx.write_tsv(df, f"{d}/{rel}", index=False)
        else:
            p = ctx.out(f"{d}/{rel}")
            build.atomic_save(p, lambda t: open(t, "w").write(df))
    ctx.write_report(f"{d}/report.txt", res.get("report", ctx.report))


def step_pack_table(ctx):
    """pack_table.tsv = the last step's table; with inventory.pack_table.order_reference (a frozen input) the rows and
    the key order of equal comp_weights dicts follow the reference (legacy_order.py); report: pack_table_order.json"""
    import json
    src = ctx.need(_last_table(ctx.cfg))
    ref = (ctx.inv.get("pack_table", {}) or {}).get("order_reference")
    if not ref:
        ctx.copy(src, "pack_table.tsv")
        return
    from . import legacy_order
    P, rep = legacy_order.reorder(src, ctx.frozen(ref))
    rep["reference"] = ref
    ctx.say(f"pack_table order reference {ref}: {json.dumps({k: v for k, v in rep.items() if k != 'examples'})}")
    if not rep.get("applied"):
        ctx.copy(src, "pack_table.tsv")
    else:
        ctx.write_tsv(P, "pack_table.tsv")
    p = ctx.out("pack_table_order.json")
    build.atomic_save(p, lambda t: json.dump(rep, open(t, "w"), indent=1))


def step_canonical(ctx):
    from . import canonical
    canonical.export(ctx)


def step_compare(ctx):
    from . import compare
    compare.run(ctx)


def run(cfg, bdir, step=None, **kw):
    build.check_env()
    _register_core()
    plan = steps(cfg)
    todo = list(plan) if step in (None, "all") else [step]
    rs = _inv(cfg).get("rule_sets", []) or []
    for s in todo:
        if s not in plan:
            raise KeyError(f"unknown inventory step '{s}' (known: {list(plan)})")
        if step in (None, "all") and all(os.path.exists(os.path.join(bdir, o)) for o in plan[s]["outputs"]):
            continue
        ctx = StepContext(cfg, bdir, s)
        print(f"[inventory] step {s}", flush=True)
        if s == "rna":
            step_rna(ctx)
        elif s == "copies":
            step_copies(ctx)
        elif s == "v0":
            step_v0(ctx)
        elif s == "secms":
            step_secms(ctx)
        elif s == "assembly":
            step_assembly(ctx)
        elif s == "ingredients":
            step_ingredients(ctx)
        elif s.startswith("rules:"):
            name = s.split(":", 1)[1]
            i = rs.index(name)
            src = "ingredients/ingredients.tsv" if i == 0 else f"{rule_dir(rs[i - 1])}/ingredients.tsv"
            step_rule(ctx, name, src)
        elif s == "pack_table":
            step_pack_table(ctx)
        elif s == "canonical":
            step_canonical(ctx)
        elif s == "compare":
            step_compare(ctx)
