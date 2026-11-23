"""Assembly state in the core: one ordered list of plain step functions, behaviour chosen by a named preset.

    inventory:
      adapters: {assembly: assembly_steps}
      params:
        assembly_steps:
          preset: legacy_v7                  # wholecell/inventory/assembly/presets/<name>.yaml
          overrides: {secms: {source_id: inv_sec_fossati_features}}   # deep-merged on the preset (same keys only)
          inputs: {complex_portal: inv_complex_portal, sifts: inv_sifts, ...}   # frozen input names (see NEEDS)
          data: {complex_priority: <repo TSV>, histone_families: <repo TSV>, nucleosome: {...}, ...}   # cell curation
          workers: 8                         # processes for the PDB chain counts

Steps (STEPS, a dict name -> function in steps.py; DEFAULT_ORDER is the order of the unification plan):
  curated_machines     curated machines built before anything else (Jurkat: nucleosomes counted from the genome size)
  complex_portal       Complex Portal complexes, greedy over the pools, with the best local PDB entry as structure
  secms_fractions      SEC-MS assembled fractions (views of the `secms` step tables, toolkits/secms)
  afdb_homodimers      AFDB homodimers on the free pool (ipSAE gate, SEC-MS share)
  afdb_heterodimers    AFDB heterodimers on the free pool (SEC-MS r or STRING evidence)
  pdb_homo_oligomers   homo-oligomers of the residual pool from PDB assemblies (SEC-MS share)
The complex steps (curated_machines, complex_portal) come first; after them the runner writes
assembly/{complexes_v0, inventory_v0_complexed, complex_copies_rule}.tsv, which every later step reads. A preset may
reorder the later steps: legacy_v7 runs pdb_homo_oligomers before the AFDB dimers, as v7 did.

Outputs (inventory/assembly/): the tables above + sec_assembly.tsv, sec_calibrants_peak.tsv, homo_oligomers_v0.tsv,
afdb_homodimers_v0.tsv, afdb_heterodimers_v0.tsv (per step listed), assembly_options.json (the SEC-MS rule options in
the wholecell/inventory/secms.py layout) and assembly_steps.json (preset, steps, resolved options).

Pools (preset option `pools`): inventory = integer pools of inventory_v0.tsv and the v7 tables (steps.py; legacy_v7);
ingredients = an ingredient table with fractional copies that every step consumes from and adds rows to
(ingredient_steps.py; legacy_x5, the yeast X5 assembly state). With pools: ingredients a preset may also list plugin
steps (registry kind `assembly_step`, e.g. the yeast pre-ribosomes), and steps call hooks (kind `assembly_hook`: copy
counts, the structure builder) named in their options. Outputs then: assembly/ingredients.csv, oligo_map.json,
complexes_tried.csv, complexes_added.csv, heterodimers.csv, oligomers.csv (+ plugin step outputs, e.g. f_actin.json)
and assembly_steps.json.

Presets: legacy_v7 (Jurkat v7, exact: B3/B4 md5), legacy_x5 (yeast port build s0 assembly tables, exact:
scripts/assembly_x5_check.py), x5_secms (legacy_x5 + the SEC-fed rules of legacy_v7: homo-oligomer / AF homodimer
shares, SEC co-elution or STRING heterodimer evidence; Step 3, docs/unification_criteria/step3_yeast_secms_assembly.md).
run_standalone() runs a step list outside a build (checks).
"""
from __future__ import annotations

import copy
import difflib
import json
import os

import yaml

from ...registry import register

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = "inventory.params.assembly_steps"

DEFAULT_ORDER = ["curated_machines", "complex_portal", "secms_fractions", "afdb_homodimers", "afdb_heterodimers",
                 "pdb_homo_oligomers"]
COMPLEX_STEPS = ("curated_machines", "complex_portal")

# allowed values of the enumerated options ("<step>.<key>[.<key>]"); anything else is a config error
ENUMS = {
    "pools": ("inventory", "ingredients"),
    "curated_machines.nucleosome": ("dna_count", "none"),
    "curated_machines.table": ("data",),
    "complex_portal.order": ("priority_then_mass", "mass"),
    "complex_portal.overlap": ("consume_left", "skip_if_used"),
    "complex_portal.copies": ("stoichiometry_cap", "measured_free_pool"),
    "complex_portal.structure.source": ("sifts_local_pdb", "rcsb"),
    "complex_portal.structure.rank": ("resolution",),
    "complex_portal.structure.stoichiometry_from": ("pdb_chains",),
    "complex_portal.structure.coverage": ("uniprot_mass",),
    "complex_portal.structure.tie_break": ("fewest_entities_then_reference", "fewest_entities"),
    "pdb_homo_oligomers.share.source": ("secms", "constant"),
    "pdb_homo_oligomers.copies": ("none", "share"),
    "pdb_homo_oligomers.af_share.source": ("secms", "constant"),
    "pdb_homo_oligomers.share.apex_check_side": ("both", "below"),
    "pdb_homo_oligomers.order": ("copies_desc",),
    "pdb_homo_oligomers.state_source": ("uniprot_oligomer_n",),
    "afdb_homodimers.share.source": ("secms", "constant"),
    "afdb_heterodimers.order": ("max_ipsae", "file"),
    "afdb_heterodimers.allocation": ("bound_fraction_of_min_free", "proportional"),
    "afdb_heterodimers.evidence": ("none",),
}

# what each step needs: (frozen inputs required, inputs optional, data required, data optional)
NEEDS = {
    "curated_machines": ([], [], ["histone_families", "nucleosome"], []),
    "complex_portal": (["complex_portal", "sifts", "pdb_assemblies"], ["pdb_listing_md5", "pdb_tiebreak_reference"], [],
                       ["complex_priority", "complex_copy_caps"]),
    "secms_fractions": ([], [], [], []),
    "pdb_homo_oligomers": (["sifts", "pdb_assemblies", "uniprot_reviewed"], ["pdb_listing_md5"], [],
                           ["homo_oligomer_exclude_regex", "apex_check_exclude"]),
    "afdb_homodimers": (["afdb_homodimers"], [], [], []),
    "afdb_heterodimers": (["afdb_heterodimers", "uniprot_full", "string_links"], [], [], []),
}

# pools: ingredients (ingredient_steps.py): (inputs required, data required, data optional)
NEEDS_ING = {
    "curated_machines": (["protein_table", "complex_portal"], ["curated_machines", "curated_machine_subunits"],
                         ["legacy_oligomers"]),
    "complex_portal": (["complex_portal", "rcsb_cache", "pdb_assemblies"], [],
                       ["copies_exclude_subunits", "shared_subunits", "structure_overrides"]),
    "afdb_heterodimers": (["afdb_heterodimers", "intact_pairs", "paralog_pairs", "proteome_fasta", "complex_portal"], [],
                          ["histone_classes", "monomeric_families"]),
    "pdb_homo_oligomers": (["afdb_homodimers", "rcsb_cache"], [], ["homo_oligomer_audit"]),
}
OUTS_ING = {"complex_portal": ["complexes_tried.csv", "complexes_added.csv"], "afdb_heterodimers": ["heterodimers.csv"],
            "pdb_homo_oligomers": ["oligomers.csv"]}
ALWAYS_ING = ["ingredients.csv", "oligo_map.json", "assembly_steps.json"]

OUTS = {
    "secms_fractions": ["sec_assembly.tsv", "sec_calibrants_peak.tsv"],
    "pdb_homo_oligomers": ["homo_oligomers_v0.tsv"],
    "afdb_homodimers": ["afdb_homodimers_v0.tsv"],
    "afdb_heterodimers": ["afdb_heterodimers_v0.tsv"],
}
ALWAYS = ["complexes_v0.tsv", "inventory_v0_complexed.tsv", "complex_copies_rule.tsv", "assembly_options.json",
          "assembly_steps.json"]


class AssemblyConfigError(ValueError):
    pass


def _near(name, choices):
    m = difflib.get_close_matches(str(name), [str(c) for c in choices], n=1, cutoff=0.5)
    return f"; did you mean '{m[0]}'?" if m else ""


def presets() -> list[str]:
    return sorted(f[:-5] for f in os.listdir(os.path.join(HERE, "presets")) if f.endswith(".yaml") and ".draft" not in f)


def load_preset(name: str) -> dict:
    if name not in presets():
        draft = os.path.exists(os.path.join(HERE, "presets", f"{name}.draft.yaml"))
        raise AssemblyConfigError(
            f"{KEY}.preset = '{name}' is not a preset (known: {presets()})"
            + (f"; presets/{name}.draft.yaml is a commented draft, not wired yet" if draft else _near(name, presets())))
    return yaml.safe_load(open(os.path.join(HERE, "presets", f"{name}.yaml")))


def _merge_checked(base, over, path):
    """deep merge `over` on `base`; every key of `over` must exist in `base` (typos are errors, not new options)"""
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if k not in out:
            raise AssemblyConfigError(f"{KEY}.overrides.{path}{k}: unknown option (allowed: {sorted(out)}){_near(k, out)}")
        if isinstance(v, dict) and isinstance(out[k], dict):
            out[k] = _merge_checked(out[k], v, f"{path}{k}.")
        else:
            out[k] = copy.deepcopy(v)
    return out


def _get(o, dotted):
    for k in dotted.split("."):
        o = o.get(k) if isinstance(o, dict) else None
    return o


def options(params: dict) -> dict:
    """preset + overrides, validated: step names, enumerated values"""
    name = params.get("preset")
    if not name:
        raise AssemblyConfigError(f"{KEY}.preset is required (known: {presets()})")
    o = _merge_checked(load_preset(name), params.get("overrides") or {}, "")
    o["preset"] = name
    steps = list(o.get("steps") or [])
    if len(set(steps)) != len(steps):
        raise AssemblyConfigError(f"{KEY} (preset {name}) steps: a step is listed twice: {steps}")
    if o.get("pools", "inventory") == "ingredients":       # absent = inventory (legacy_v7 json unchanged)
        _check_ingredient_steps(o, steps, name)
        return o
    for s in steps:
        if s not in DEFAULT_ORDER:
            raise AssemblyConfigError(f"{KEY} (preset {name}) steps: unknown step '{s}' (known: {DEFAULT_ORDER}){_near(s, DEFAULT_ORDER)}")
    cx = [s for s in steps if s in COMPLEX_STEPS]
    if steps[:len(cx)] != cx or cx != [s for s in COMPLEX_STEPS if s in cx]:
        raise AssemblyConfigError(f"{KEY} (preset {name}) steps: the complex steps {list(COMPLEX_STEPS)} come first, in "
                                  f"that order (got {steps})")
    if "secms_fractions" in steps and "pdb_homo_oligomers" in steps and \
            steps.index("secms_fractions") > steps.index("pdb_homo_oligomers"):
        raise AssemblyConfigError(f"{KEY} (preset {name}) steps: secms_fractions must come before pdb_homo_oligomers "
                                  f"(its SEC apex check reads sec_assembly.tsv)")
    _check_enums(o, name)
    o["steps"] = steps
    return o


def _check_enums(o, name):
    for k, allowed in ENUMS.items():
        v = _get(o, k)
        if v is not None and not isinstance(v, dict) and v not in allowed:
            raise AssemblyConfigError(f"{KEY}.overrides.{k} = '{v}' (preset {name}): allowed {list(allowed)}{_near(v, allowed)}")


def _check_ingredient_steps(o, steps, name):
    """pools: ingredients: core steps of ingredient_steps.py + registered plugin steps; curated_machines first (it makes
    the ingredient table); hook names must be registered"""
    from ... import registry
    from . import ingredient_steps as IS
    plug = registry.names("assembly_step")
    known = list(IS.STEPS) + plug
    for s in steps:
        if s not in known:
            raise AssemblyConfigError(f"{KEY} (preset {name}) steps: unknown step '{s}' (core: {list(IS.STEPS)}; plugin "
                                      f"steps registered: {plug}; is the plugin listed in the cell config?){_near(s, known)}")
    if not steps or steps[0] != "curated_machines":
        raise AssemblyConfigError(f"{KEY} (preset {name}) steps: with pools: ingredients the list starts with "
                                  f"curated_machines (it makes the ingredient table); got {steps}")
    hooks = registry.names("assembly_hook")
    for k in ("curated_machines.legacy_oligomers_require_library", "complex_portal.structure.builder",
              "afdb_heterodimers.builder", "pdb_homo_oligomers.builder"):
        v = _get(o, k)
        if v and v not in hooks:
            raise AssemblyConfigError(f"{KEY}.overrides.{k} = '{v}' (preset {name}): not a registered assembly_hook "
                                      f"(registered: {hooks}){_near(v, hooks)}")
    _check_enums(o, name)
    o["steps"] = steps


def check_params(params: dict, o: dict, has_secms: bool):
    """every input / data entry the listed steps need is present"""
    inputs, data = params.get("inputs") or {}, params.get("data") or {}
    known_in = {x for s in NEEDS.values() for x in s[0] + s[1]} | {x for s in NEEDS_ING.values() for x in s[0]}
    known_data = {x for s in NEEDS.values() for x in s[2] + s[3]} | {x for s in NEEDS_ING.values() for x in s[1] + s[2]}
    if o.get("pools") == "ingredients":
        for s in o["steps"]:
            if s not in NEEDS_ING:
                continue                                   # plugin step: checks its own inputs
            req_in, req_data, _ = NEEDS_ING[s]
            for k in req_in:
                if not inputs.get(k):
                    raise AssemblyConfigError(f"{KEY}.inputs.{k} is required by step {s}{_near(k, inputs)}")
            for k in req_data:
                if data.get(k) in (None, ""):
                    raise AssemblyConfigError(f"{KEY}.data.{k} is required by step {s}{_near(k, data)}")
        _unknown(inputs, data, known_in, known_data)
        return
    for s in o["steps"]:
        req_in, _, req_data, _ = NEEDS[s]
        if s == "curated_machines" and o["curated_machines"]["nucleosome"] == "none":
            req_data = []
        if s == "afdb_heterodimers" and o["afdb_heterodimers"]["evidence"].get("string_min") is None:
            req_in = [x for x in req_in if x not in ("string_links", "uniprot_full")]
        for k in req_in:
            if not inputs.get(k):
                raise AssemblyConfigError(f"{KEY}.inputs.{k} is required by step {s} (a declared frozen input name)"
                                          f"{_near(k, inputs)}")
        for k in req_data:
            if data.get(k) in (None, ""):
                raise AssemblyConfigError(f"{KEY}.data.{k} is required by step {s}{_near(k, data)}")
    _unknown(inputs, data, known_in, known_data)
    if "secms_fractions" in o["steps"] and not has_secms:
        raise AssemblyConfigError(f"step secms_fractions needs inventory.adapters.secms (the SEC-MS toolkit step)")


def _unknown(inputs, data, known_in, known_data):
    for k in inputs:
        if k not in known_in:
            raise AssemblyConfigError(f"{KEY}.inputs.{k}: unknown input (known: {sorted(known_in)}){_near(k, known_in)}")
    for k in data:
        if k not in known_data:
            raise AssemblyConfigError(f"{KEY}.data.{k}: unknown entry (known: {sorted(known_data)}){_near(k, known_data)}")


def outputs(cfg) -> list[str]:
    p = (((cfg.get("inventory") or {}).get("params") or {}).get("assembly_steps") or {})
    o = options(p)
    if o.get("pools") == "ingredients":
        from ... import registry
        plug = set(registry.names("assembly_step"))
        extra = [f for s in o["steps"] for f in (OUTS_ING.get(s, []) if s not in plug else
                                                   getattr(registry.get("assembly_step", s), "outputs", []))]
        return [f"assembly/{x}" for x in ALWAYS_ING + extra]
    return [f"assembly/{x}" for x in ALWAYS + [f for s in o["steps"] for f in OUTS.get(s, [])]]


def assembly_steps(ctx):
    from .. import secms as S2
    from . import steps as S
    p = ctx.params("assembly_steps")
    o = options(p)
    has_secms = bool((ctx.inv.get("adapters") or {}).get("secms"))
    check_params(p, o, has_secms)
    if o.get("pools") == "ingredients":
        return run_ingredients(ctx, p, o, has_secms=has_secms)
    st = S.State(ctx, p, o, S2.SecView(ctx, o["secms"]) if has_secms else None)
    from ... import build
    legacy = S.legacy_options(o)
    sel = st.view
    info = dict(source_id=sel.source_id, condition=sel.condition, replicate=sel.replicate) if sel is not None else None
    build.atomic_save(ctx.out("assembly/assembly_options.json"),
                      lambda t: json.dump(dict(legacy, secms_used=info), open(t, "w"), indent=1))
    build.atomic_save(ctx.out("assembly/assembly_steps.json"),
                      lambda t: json.dump(dict(preset=o["preset"], steps=o["steps"], options=o,
                                               overrides=p.get("overrides") or {}), open(t, "w"), indent=1))
    ctx.say(f"assembly preset {o['preset']}: steps {' -> '.join(o['steps'])}")
    done_complexes = False
    for s in o["steps"]:
        if s not in COMPLEX_STEPS and not done_complexes:
            S.write_complex_tables(st); done_complexes = True
        print(f"[assembly] {s}", flush=True)
        S.STEPS[s](st)
    if not done_complexes:
        S.write_complex_tables(st)


def run_ingredients(ctx, p, o, df=None, has_secms=False):
    """pools: ingredients: the step list on one ingredient table (core steps + plugin steps), then ingredients.csv and
    oligo_map.json. A preset with a `secms` selection (x5_secms) reads the SEC-MS tables of the `secms` step"""
    from ... import build, registry
    from . import ingredient_steps as IS
    st = IS.IngredientState(ctx, p, o, df=df)
    if "secms" in o:
        if not has_secms:
            raise AssemblyConfigError(f"{KEY} (preset {o['preset']}): SEC-fed rules need inventory.adapters.secms")
        from .. import secms as S2
        st.view = S2.SecView(ctx, o["secms"])
        ctx.say(f"SEC-MS rules read {st.view.source_id} / {st.view.condition} / {st.view.replicate}")
    build.atomic_save(ctx.out("assembly/assembly_steps.json"),
                      lambda t: json.dump(dict(preset=o["preset"], steps=o["steps"], options=o,
                                               overrides=p.get("overrides") or {}), open(t, "w"), indent=1))
    ctx.say(f"assembly preset {o['preset']} (pools: ingredients): steps {' -> '.join(o['steps'])}")
    for s in o["steps"]:
        print(f"[assembly] {s}", flush=True)
        (IS.STEPS[s] if s in IS.STEPS else registry.get("assembly_step", s))(st)
    IS.finish(st)
    return st


def run_standalone(out_dir, params, cell="yeast", plugins=(), protein_table=None):
    """run a step list outside a build (checks, scripts/assembly_x5_check.py): outputs in <out_dir>/inventory/assembly;
    frozen input names resolve through config/cells/<cell>.yaml"""
    import importlib
    from ...config import load
    from ..context import StepContext
    for m in plugins:
        importlib.import_module(m)
    cfg = load(cell)
    params = copy.deepcopy(params)
    if protein_table:
        params.setdefault("inputs", {})["protein_table"] = protein_table
    o = options(params)
    check_params(params, o, False)
    os.makedirs(os.path.join(out_dir, "inventory", "assembly"), exist_ok=True)
    ctx = StepContext(cfg, out_dir, "assembly_steps")
    if o.get("pools") != "ingredients":
        raise AssemblyConfigError("run_standalone runs pools: ingredients presets (legacy_v7 runs in a build)")
    return run_ingredients(ctx, params, o)


assembly_steps.outputs = outputs
assembly_steps.PREV = []
register("adapter.assembly", "assembly_steps")(assembly_steps)
