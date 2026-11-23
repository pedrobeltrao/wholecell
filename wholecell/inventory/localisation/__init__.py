"""Localisation in the core: source adapters -> one evidence table -> common combiner -> per-cell override tables.

    inventory:
      adapters: {localisation: localisation}
      params:
        localisation:
          preset: legacy_v7                 # wholecell/inventory/localisation/presets/<name>.yaml
          overrides: {combiner: {combine: {weights: {...}}}}     # deep-merged on the preset (existing keys only)
          label: jurkat.hpa_uniprot_itzhak_v0   # provenance name in species_loc / field_provenance (default: localisation)
          genes: {input: inv_uniprot_full}  # gene -> UniProt accession + transmembrane feature (reviewed entries)
          sources:                          # evidence adapters (evidence.SOURCES), in evidence-table order
            - {source: hpa, adapter: hpa_subcellular, evidence_class: imaging_other, input: inv_hpa_subcellular,
               term_map: <repo TSV>, main_score: 1.0, additional_score: 0.3}
            - ...
          data: {membrane_keys: <repo TSV>, overrides: <repo TSV>}   # cell curation (TSVs with a reason column)

Step v0 calls the adapter: every source adapter writes rows of inventory/localisation/evidence.tsv (columns:
evidence.COLUMNS), the per-gene features go to gene_features.tsv, the combiner (combiner.STEPS, in the preset's
`combiner.steps` order) turns the evidence into weights (`compartments` of v0/inventory_v0.tsv) and
localisation_options.json records preset + resolved options.

Overrides (overrides.py): curated revisions are rows of a TSV (rule_set, rule_id, key, field, op, value, reason,
citation) applied by overrides.apply. The preset's `rules` section lists, per rule set, every step of the cell's revision
rules in order and says what each one is: `overrides` (rows of data.overrides), `frozen_overrides` (an override table
generated at build time from declared frozen inputs, e.g. earlier packs; legacy presets only) or `legacy` (a named rule
function still in code, with a one-line reason). The rule-set code checks its step order against this list
(check_rule_steps) and reads its override rows through override_rows().

Presets: common (step B5: one set of defaults for every cell; the cell data block `common:` of the localisation params
supplies vocabulary, key-selection data, imaging groups and the N/C calibration; docs/unification_criteria/
B5_localisation.md); legacy_v7 (Jurkat v7 exactly); legacy_x5 (yeast X5 / port build s0 comp_weights exactly: ordered-annotation
steps; cell data plugins/yeast/data/localisation/, check scripts/localisation_x5_check.py; not yet wired into the yeast
DAG, A9).
"""
from __future__ import annotations

import copy
import difflib
import json
import os

import yaml

from ...registry import register
from . import combiner, evidence, overrides

HERE = os.path.dirname(os.path.abspath(__file__))
KEY = "inventory.params.localisation"
OUTPUTS = ["localisation/evidence.tsv", "localisation/gene_features.tsv", "localisation/localisation_options.json",
           "localisation/nc_model.tsv"]     # nc_model.tsv: per-gene N/C model record (preset common; empty otherwise)
RULE_KINDS = ("overrides", "frozen_overrides", "legacy", "rule")
ENUMS = {
    "combiner.within_source": ("max", "sum"),
    "combiner.split": ("score", "rank", "equal"),
    "combiner.combine.mode": ("blend", "precedence", "classes"),
    "combiner.combine.key_order": ("set_union", "first_seen"),
    "combiner.nc_model.method": ("fraction", "ratio_exp", "ratio"),
}


class LocalisationConfigError(ValueError):
    pass


def _near(name, choices):
    m = difflib.get_close_matches(str(name), [str(c) for c in choices], n=1, cutoff=0.5)
    return f"; did you mean '{m[0]}'?" if m else ""


def presets() -> list[str]:
    return sorted(f[:-5] for f in os.listdir(os.path.join(HERE, "presets")) if f.endswith(".yaml") and ".draft" not in f)


def load_preset(name: str) -> dict:
    if name not in presets():
        draft = os.path.exists(os.path.join(HERE, "presets", f"{name}.draft.yaml"))
        raise LocalisationConfigError(
            f"{KEY}.preset = '{name}' is not a preset (known: {presets()})"
            + (f"; presets/{name}.draft.yaml is a commented draft, not wired yet" if draft else _near(name, presets())))
    return yaml.safe_load(open(os.path.join(HERE, "presets", f"{name}.yaml")))


def _merge_checked(base, over, path):
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if k not in out:
            raise LocalisationConfigError(f"{KEY}.overrides.{path}{k}: unknown option (allowed: {sorted(out)}){_near(k, out)}")
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
    """preset + overrides, validated (combiner step names, enumerated values, rule-step kinds)"""
    name = params.get("preset")
    if not name:
        raise LocalisationConfigError(f"{KEY}.preset is required (known: {presets()})")
    o = _merge_checked(load_preset(name), params.get("overrides") or {}, "")
    o["preset"] = name
    cb = o["combiner"]
    for s in cb["steps"]:
        if s not in combiner.STEPS:
            raise LocalisationConfigError(f"{KEY} (preset {name}) combiner.steps: unknown step '{s}' (known: "
                                          f"{list(combiner.STEPS)}){_near(s, combiner.STEPS)}")
        if s not in ("drop_unreliable",) + combiner.OPTIONAL and cb.get(s) is None:
            raise LocalisationConfigError(f"{KEY} (preset {name}) combiner.steps lists '{s}' but combiner.{s} is not set")
    for k, allowed in ENUMS.items():
        v = _get(o, k)
        if v is not None and v not in allowed:
            raise LocalisationConfigError(f"{KEY}.overrides.{k} = '{v}' (preset {name}): allowed {list(allowed)}{_near(v, allowed)}")
    for c in (cb.get("combine") or {}).get("classes", []):
        if c not in evidence.CLASSES:
            raise LocalisationConfigError(f"{KEY} (preset {name}) combiner.combine.classes: '{c}' is not an evidence class "
                                          f"(known: {list(evidence.CLASSES)}){_near(c, evidence.CLASSES)}")
    for rs, steps in (o.get("rules") or {}).items():
        for st in steps:
            if st.get("kind") not in RULE_KINDS:
                raise LocalisationConfigError(f"{KEY} (preset {name}) rules.{rs} step {st.get('id')}: kind '{st.get('kind')}' "
                                              f"(allowed: {list(RULE_KINDS)}){_near(st.get('kind'), RULE_KINDS)}")
            if st["kind"] in ("legacy", "rule") and not st.get("reason"):
                raise LocalisationConfigError(f"{KEY} (preset {name}) rules.{rs} step {st['id']}: a {st['kind']} step needs a reason")
    return o


def _data_path(params, key):
    v = (params.get("data") or {}).get(key)
    if not v:
        raise LocalisationConfigError(f"{KEY}.data.{key} is required{_near(key, params.get('data') or {})}")
    return v


def _resolve_data(o, params):
    """option values of the form 'data.<name>' -> the cell's data entry (a TSV: first column), 'rows.<name>' -> its rows
    (list of dicts, e.g. the override table), 'cell.<dotted path>' -> that value of the cell's localisation params
    (preset common: the cell data block `common:`; None when the cell has none, which skips an optional step); only the
    options of the combiner steps listed (an unused step needs no data)"""
    def walk(x):
        if isinstance(x, dict):
            return {k: walk(v) for k, v in x.items()}
        if isinstance(x, list):
            return [walk(v) for v in x]
        if isinstance(x, str) and x.startswith("cell."):        # preset common: the cell's own data block (or None)
            return walk(copy.deepcopy(_get(params, x[5:])))
        if isinstance(x, str) and x.startswith("data."):
            return list(evidence.read_curated(_data_path(params, x[5:])).iloc[:, 0])
        if isinstance(x, str) and x.startswith("rows."):
            return evidence.read_curated(_data_path(params, x[5:])).to_dict("records")
        return x
    return {k: (walk(v) if k in o["steps"] else v) for k, v in o.items()}


# ----------------------------------------------------------------------------------------------- the v0 adapter
class _SourceCtx:
    def __init__(self, ctx, up, uniprot_of):
        self.ctx, self.uniprot_table, self.uniprot_of = ctx, up, uniprot_of

    def resolve(self, v):
        return self.ctx.resolve(v)


def localisation(ctx, prot):
    import pandas as pd
    params = ctx.params("localisation")
    o = options(params)
    cb = _resolve_data(o["combiner"], params)
    cal = (cb.get("nc_model") or {}).get("calibration") or {}
    if cal.get("exclude_genes_input"):      # preset common: genes kept out of the imaging N/C calibration (hold-out list)
        f = ctx.resolve(cal["exclude_genes_input"])
        cal["exclude_genes"] = sorted({x.strip() for x in open(f) if x.strip() and not x.startswith("#")})
    gp = params.get("genes") or {}
    kc = gp.get("key")                      # column of the proteome table that holds the evidence key (yeast: orf);
    genes = list(prot[kc].astype(str)) if kc else list(prot.index)   # default: the proteome table's index
    if gp.get("reader") == "uniprot_ordered_locus":   # keyed by ordered locus (yeast ORF): uniprot, gene_name, n_tm
        up = None
        F = evidence.locus_features(ctx.resolve(gp["input"]), genes=genes)
    elif gp.get("input"):                   # UniProt dump: reviewed entries by gene (accession, Transmembrane feature)
        up = evidence.uniprot_genes(ctx.resolve(gp["input"]))
        F = evidence.gene_features(genes, up)
    elif gp.get("table"):                   # TSV gene, transmembrane (e.g. the testcell localisation table)
        T = pd.read_csv(ctx.resolve(gp["table"]), sep="\t", comment="#").drop_duplicates("gene").set_index("gene")
        up = None
        F = pd.DataFrame(index=pd.Index(genes, name="gene"))
        F["uniprot"] = prot["uniprot"].reindex(F.index) if "uniprot" in prot.columns else None
        F["transmembrane"] = T.transmembrane.reindex(F.index).fillna(False).astype(bool).values
    else:
        raise LocalisationConfigError(f"{KEY}.genes: set input (a UniProt dump: gene -> accession, transmembrane), "
                                      f"reader: uniprot_ordered_locus + input, or table (TSV gene, transmembrane)")
    sctx = _SourceCtx(ctx, up, F.uniprot.dropna().to_dict())
    rows = []
    seen = set()
    for sp in params.get("sources") or []:
        for k in ("source", "adapter", "evidence_class"):
            if not sp.get(k):
                raise LocalisationConfigError(f"{KEY}.sources[{sp.get('source', '?')}].{k} is required")
        if not sp.get("input") and not sp.get("inputs"):
            raise LocalisationConfigError(f"{KEY}.sources[{sp.get('source', '?')}].input is required (or inputs: a list)")
        if sp["adapter"] not in evidence.SOURCES:
            raise LocalisationConfigError(f"{KEY}.sources[{sp['source']}].adapter '{sp['adapter']}' (known: "
                                          f"{list(evidence.SOURCES)}){_near(sp['adapter'], evidence.SOURCES)}")
        if sp["evidence_class"] not in evidence.CLASSES:
            raise LocalisationConfigError(f"{KEY}.sources[{sp['source']}].evidence_class '{sp['evidence_class']}' (known: "
                                          f"{list(evidence.CLASSES)}){_near(sp['evidence_class'], evidence.CLASSES)}")
        if sp["adapter"] == "uniprot_cc" and up is None:
            raise LocalisationConfigError(f"{KEY}.sources[{sp['source']}]: adapter uniprot_cc needs genes.input (a UniProt dump)")
        if sp["source"] in seen:
            raise LocalisationConfigError(f"{KEY}.sources: source label '{sp['source']}' used twice")
        seen.add(sp["source"])
        r = evidence.SOURCES[sp["adapter"]](genes, sp, sctx)
        ctx.say(f"evidence {sp['source']} ({sp['adapter']}, {sp['evidence_class']}): {len(r)} rows, "
                f"{len({x[0] for x in r})} genes")
        rows += r
    ctx.write_tsv(evidence.to_frame(rows), "localisation/evidence.tsv")
    ctx.write_tsv(F, "localisation/gene_features.tsv", index=True)
    from ... import build
    build.atomic_save(ctx.out("localisation/localisation_options.json"),
                      lambda t: json.dump(dict(preset=o["preset"], combiner=cb, sources=params.get("sources"),
                                               overrides=params.get("overrides") or {}), open(t, "w"), indent=1))
    E = evidence.read(ctx.path("localisation/evidence.tsv"))              # the combiner reads the written table
    S = combiner.run(E, F, genes, cb, ctx)
    ctx.write_tsv(combiner.nc_table(S), "localisation/nc_model.tsv")
    out = pd.DataFrame(index=prot.index)
    out["compartments"] = [combiner.compartments_string(S[g].w) for g in genes]
    if kc:                                  # evidence keyed by another column: rows aligned by position
        out["loc_source"] = ["+".join(S[g].label) for g in genes]
        out["transmembrane"] = F.transmembrane.reindex(genes).values
    else:
        out["loc_source"] = pd.Series({g: "+".join(S[g].label) for g in genes})
        out["transmembrane"] = F.transmembrane.reindex(prot.index).values
    ctx.say(f"localisation preset {o['preset']}: combiner {' -> '.join(cb['steps'])}; loc_source "
            + ", ".join(f"{k} {v}" for k, v in out.loc_source.value_counts().items()))
    return out


localisation.OUTPUTS = OUTPUTS


def provenance_label(cfg) -> str:
    return ((((cfg.get("inventory") or {}).get("params") or {}).get("localisation") or {}).get("label") or "localisation")


localisation.provenance_label = provenance_label
register("adapter.localisation", "localisation")(localisation)


# -------------------------------------------------------------------------------- rule sets: plan + override rows
def rule_steps(ctx, rule_set: str) -> list[dict]:
    o = options(ctx.params("localisation"))
    rs = (o.get("rules") or {})
    if rule_set not in rs:
        raise LocalisationConfigError(f"{KEY} (preset {o['preset']}) rules: no entry for rule set '{rule_set}' (known: "
                                      f"{list(rs)}){_near(rule_set, rs)}")
    return rs[rule_set]


def check_rule_steps(ctx, rule_set: str, code_steps: list[tuple]) -> list[dict]:
    """the rule-set code's step order [(id, kind)] must equal the preset's rules.<rule_set> list; a code step may give a
    tuple of kinds it implements (e.g. ("frozen_overrides", "rule") for v5.P: legacy_v7 vs common)"""
    plan = rule_steps(ctx, rule_set)
    got = [(s["id"], s["kind"]) for s in plan]
    ok = len(got) == len(code_steps) and all(
        g[0] == c[0] and (g[1] in c[1] if isinstance(c[1], tuple) else g[1] == c[1]) for g, c in zip(got, code_steps))
    if not ok:
        raise LocalisationConfigError(f"{KEY} preset rules.{rule_set} = {got} but the rule-set code runs {list(code_steps)}; "
                                      f"the preset lists every step of the code in order (id, kind)")
    return plan


def step_kind(ctx, rule_set: str, rule_id: str) -> str:
    """the kind the preset gives one rule step (legacy_v7: frozen_overrides; common: rule)"""
    for st in rule_steps(ctx, rule_set):
        if st["id"] == rule_id:
            return st["kind"]
    raise LocalisationConfigError(f"{KEY}: rule {rule_set} {rule_id} is not in the preset's rules")


def override_table(ctx):
    params = ctx.params("localisation")
    p = _data_path(params, "overrides")
    return overrides.check(evidence.read_curated(p), p)


def override_rows(ctx, rule_set: str, rule_ids):
    """the data.overrides rows of these rule ids (each id must be an `overrides` step of the preset and have rows)"""
    rule_ids = [rule_ids] if isinstance(rule_ids, str) else list(rule_ids)
    plan = {s["id"]: s["kind"] for s in rule_steps(ctx, rule_set)}
    for r in rule_ids:
        if plan.get(r) != "overrides":
            raise LocalisationConfigError(f"{KEY}: rule {rule_set} {r} is not an `overrides` step of the preset "
                                          f"(steps: {plan}){_near(r, plan)}")
    T = overrides.select(override_table(ctx), rule_set, rule_ids)
    missing = [r for r in rule_ids if r not in set(T.rule_id)]
    if missing:
        raise LocalisationConfigError(f"{KEY}.data.overrides has no rows for {rule_set} {missing}")
    return T
