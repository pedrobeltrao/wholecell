"""Packer algorithm switches (unification plan A7, phase 1): the ONE place where their defaults live.

Each switch ports an algorithm of the yeast X5 packer (plugins/yeast/engine/pipeline/wc3d.py, squeeze.py @ e3b7448)
into the core packer. Every default below is the Jurkat v7 behaviour: with all defaults a build is bit-identical to
the v7 packer (docs/stages/packing.md). The cell config sets them under `pack:`; keys missing there take these values.

    from wholecell.pack import switches
    sw = switches.load(cfg["pack"], region_names)      # validated, defaults filled in; sw.rounds["n"], sw.squeeze, ...
"""
from __future__ import annotations

import copy
import difflib

DEFAULTS = {
    # P1 size-class rounds (origin: yeast wc3d.round_edges / job_rounds / s3_pack round loop).
    #   n: number of visits per tile (1 = v7: one visit, all jobs). With n > 1 the volume and shell jobs are split by
    #      library bounding radius into n classes (largest first) and round r packs class r in every tile;
    #      membrane jobs are packed in round 0.
    #   split: how the classes are cut. equal_volume = equal summed sphere volume (sum r^3 x copies) per class.
    #   alternate: reverse the order of the 8 checkerboard phases in every second round (each tile is as often ahead of
    #      a neighbour as behind it; fixes tile-face density seams).
    "rounds": {"n": 1, "split": "equal_volume", "alternate": True},
    # P2 pilot-calibrated thinning (origin: yeast wc3d.pilot_capacity / thin). Off in v7.
    #   pilot_tiles: tiles packed per crowded region (no neighbours, one visit) to measure what it holds.
    #   min_rel: a volume region is tested if its requested sphere volume fraction > min_rel x the reference region's.
    #   safety: request cut to safety x held (dimensionless). keep_above: regions that held >= this are left alone.
    #   reference_region: the region whose requested density is the yardstick (yeast: cytosol).
    #   pilot_squeeze: the pilot tiles are squeezed (P5 settings) before what they held is measured (yeast: the pilot
    #      ran the full tile packer, squeeze included). resplit: legacy_x5 requests are re-split from the thinned
    #      totals (yeast wc3d.distribute(thin(...)), seed stream pack.counts) instead of thinning each count binomially.
    "thinning": {"enabled": False, "pilot_tiles": 16, "min_rel": 1.1, "safety": 0.97, "keep_above": 0.98,
                 "reference_region": "cytosol", "pilot_squeeze": False, "resplit": False},
    # P3 shell cap (origin: yeast wc3d.shell_capacity). Off in v7.
    #   A shell (membrane face) asks for at most the requested sphere volume fraction of the bulk volume region on its
    #   side (the volume region holding the shell's first label); the excess copies of every species of the shell go
    #   to that region. Shell volume = membrane area x (d1 - d0, + mean bounding radius for add_radius shells).
    #   Pore-limited shells (pore_limit_nm > 0) are not capped (their volume is not area x width).
    "shell_cap": {"enabled": False},
    # P4 per-region try budgets (origin: yeast wc3d.thin_budget / pack3d reg_tries, reg_give_up).
    #   {region name: {max_tries: int, give_up: int}}: tries per copy and failures in a row for jobs of that region,
    #   replacing the global max_tries / give_up and the thin-lumen factors there. Empty = v7 budgets.
    "budgets": {"per_region": {}},
    # P5 squeeze relaxation with exact RSA re-acceptance (origin: yeast squeeze.py, wc3d.squeeze_tile). Off in v7.
    #   scope: none | underfilled (only volume regions of a tile below `threshold` placed by mass after the rounds;
    #      only their molecules move) | all (every tile with a deficit; all label-volume molecules of the tile move,
    #      missing copies are inserted where there is a deficit: the yeast squeeze).
    #   threshold: placed mass fraction below which a (tile, region) pair is squeezed (scope underfilled).
    #   rounds: insert-relax-accept rounds per tile; iterations: relaxation steps per round;
    #   insert_frac: missing copies inserted per round, as a share of the movable molecules (at least 50).
    #   Only label volume regions are squeezed (not shells, membranes or virtual regions).
    "squeeze": {"scope": "none", "threshold": 0.97, "rounds": 4, "iterations": 60, "insert_frac": 0.15},
    # P7 type-specific exclusion zones (origin: yeast wc3d.exclusion_zones / pack3d._in_exc). None in v7.
    #   list of {name, libs: [packed library ids], lib_prefixes: [...], around: {libs: [...], lib_prefixes: [...]},
    #   r_nm}: copies of the matching packed species (volume and shell jobs) may not have their centre within r_nm of
    #   a pre-placed instance of the `around` libraries (pre-pass index ids).
    #   around: {provider: <name>, params: {...}} instead: zone spheres (x, y, z, r; model nm) made once per build by a
    #   registered `exclusion_zones` provider (cell geometry features, e.g. the yeast endocytic sites); r_nm unused.
    "exclusion": [],
    # P6 organisation layer (origin: yeast pipeline/organisation.py @ e3b7448; code: wholecell/pack/organisation.py).
    #   Off in v7. Made by `pack --step organisation` (pack/organisation/network.npz, partners.npz), used by the counts
    #   step (tile split x exp(field log weight), crowding correction) and the pack step (partner attraction).
    #   string: STRING links file (taxon, version: recorded, checked against the file name), min combined score, and
    #      how STRING protein ids map to UniProt: ids.table + ids.columns (a UniProt table column holding the id after
    #      the taxon prefix, e.g. yeast ordered locus, or the full STRING id, e.g. the UniProt `STRING` cross-reference)
    #      or ids.aliases (STRING aliases file, first alias of source ids.source per STRING id).
    #   genes: UniProt table + gene-name columns: gene symbols of complexes -> UniProt (subunits without a uniprot field).
    #   subunits: table with ing_id, uniprot, genes (default: the inventory pack table) and optional lib_map (table with
    #      ing_id, lib_id: several ingredients per packed type, yeast); ids: npy of the network's type order (default:
    #      the pack-table ids in the library, or pack.requests.lib_order for legacy_x5).
    #   pool_scores: AlphaFold 3 pool pair scores (parquet a, b, pae, or the summary_models parquet) for contact
    #      partners, with paralogs (csv gene1_uniprot, gene2_uniprot). null = off (published builds: STRING only).
    #   The other keys are the yeast constants (organisation.py): PPR damping and profile, field components / length /
    #   Fourier features / beta, pool PAE cut, partner counts and strengths, attraction parameters.
    #   dense_max_nodes: networks up to this size use the dense PPR inverse (yeast), larger ones a Neumann series.
    #   field_cap_regions (freeze fix 2026-11-19): crowded volume regions in which the field may only lower a tile's
    #      requested mass below the unmodulated split, never raise it (requests.field_cap; the excess goes to the tiles
    #      below it, so the per-tile mass request equals the unmodulated one while the composition stays field-driven).
    #      Rule for the list: the regions whose placed mass fraction in the cell's reference build is < 0.95 (a field
    #      that piles mass into some tiles of a region the packer cannot fill loses it there: Jurkat nucleolus 0.864
    #      -> 0.807), except thin lumens (pack.thin_lumens: their deficit is membrane clearance, not crowding).
    #      [] = off (B2 as accepted).
    "organisation": {"enabled": False,
                     "string": {"taxon": None, "version": None, "links": None, "min_score": 400,
                                "ids": {"table": None, "columns": [], "value": "Entry", "aliases": None, "source": "UniProt_AC"}},
                     "genes": {"table": None, "columns": [], "value": "Entry"},
                     "subunits": {"table": None, "lib_map": None, "ids": None},
                     "pool_scores": None, "paralogs": None,
                     "damping": 0.5, "profile": "nodiag", "field_k": 128, "field_l_nm": 500.0, "n_features": 128,
                     "beta": 0.8, "pool_pae_max": 8.0, "max_partners": 16, "r_loose": 0.4, "a_contact": 0.6,
                     "a_loose": 0.25, "k_half": 8.0, "p_max": 0.6, "attract_tries": 40, "loose_gap_nm": [2.0, 15.0],
                     "dense_max_nodes": 8000, "field_cap_regions": []},
    # P8 packing shapes (origin: yeast wc3d.pack_trees, PACK_PER 15; written by the library stage, step pack_shapes).
    #   source: library (v7: the library's own sphere trees st_c / st_r) | pack_shapes (library/pack_shapes/
    #   pack_trees_<per>.npz, every packed and pre-placed species). file: explicit trees file (default: the one
    #   pack_trees_*.npz of <build>/library/pack_shapes). Shorthand `shapes: pack_shapes`. library: directory of the
    #   npz files the trees were made from when that is not the sealed library (yeast: engine/lib = sealed library +
    #   assembly products + pre-placed models, some sealed entries rewritten); its files override / extend the manifest
    #   entries and the stale-tree guard checks each tree's `sig` against them (md5, or the engine's size:mtime).
    "shapes": {"source": "library", "file": None},           # optional key `library` (not a default: unchanged records)
}

# ------------------------------------------------------------------ common defaults (B1)
# Keys a common default / benchmark setting may set: the packer algorithms (budgets included). Everything else in
# `pack:` (regions, tiles, requests, thin lumens, exclusion zones, organisation layer, pre-placed index, ...) is cell
# config. `shapes` is listed: P8 is decided by the B1 benchmark too.
ALGORITHM_KEYS = ("max_tries", "give_up", "thin_try_f", "thin_giveup_f", "rounds", "thinning", "shell_cap", "budgets",
                  "squeeze", "shapes")

# THE ONE PLACE that sets the common packer defaults (unification plan §2b, step B1; criteria
# docs/unification_criteria/B1_packer.md). Applied as the built-in preset `common_defaults` to every cell that does not
# pin a preset (pack.preset: legacy_v7 | legacy_x5 | none | a list). Empty = the v7 values (DEFAULTS above and the
# cell's own pack: keys). After the B1 benchmark: paste `common_defaults` of its selection.json here.
COMMON_DEFAULTS: dict = {}

# B1 benchmark settings (criteria B1_packer.md, Part A3). Every algorithm key is set explicitly, so a setting does not
# depend on DEFAULTS, on COMMON_DEFAULTS or on the cell's own values; merged after the cell's pinned preset (yeast:
# preset list [legacy_x5, b1_<s>]). `budgets.per_region: None` removes a preset's per-region budgets.
_SQ_OFF = {"scope": "none", "threshold": 0.97, "rounds": 4, "iterations": 60, "insert_frac": 0.15}
_TH_OFF = {"enabled": False, "pilot_tiles": 16, "min_rel": 1.1, "safety": 0.97, "keep_above": 0.98,
           "reference_region": "cytosol", "pilot_squeeze": False}
_V7_BUDGETS = {"max_tries": 200, "give_up": 25, "thin_try_f": 10.0, "thin_giveup_f": 16.0, "budgets": {"per_region": None}}


def _setting(rounds=1, thinning=False, squeeze="none", budgets="v7", shapes=None):
    s = copy.deepcopy(_V7_BUDGETS) if budgets == "v7" else {"max_tries": 2000, "give_up": 1200, "thin_try_f": 1.0,
                                                           "thin_giveup_f": 1.0, "budgets": {"per_region": None}}
    s["rounds"] = {"n": rounds, "split": "equal_volume", "alternate": True}
    s["thinning"] = dict(_TH_OFF, enabled=bool(thinning), pilot_squeeze=bool(thinning) and squeeze != "none")
    s["shell_cap"] = {"enabled": bool(thinning)}                 # P3 goes with P2 (criteria A3)
    s["squeeze"] = dict(_SQ_OFF, scope=squeeze)
    if shapes:
        s["shapes"] = {"source": shapes}
    return s


B1_BENCHMARK = {
    "b1_s1": _setting(),                                                        # v7 settings
    "b1_s2": _setting(rounds=6),                                                # + rounds
    "b1_s3": _setting(rounds=6, thinning=True),                                 # + thinning (+ shell cap)
    "b1_s4": _setting(rounds=6, thinning=True, squeeze="underfilled"),          # + targeted squeeze
    "b1_s5": _setting(rounds=6, thinning=True, squeeze="all"),                  # + full squeeze (yeast only)
    "b1_s4y": _setting(rounds=6, thinning=True, squeeze="underfilled", budgets="x5"),   # S4, 2000 / 1200 everywhere
    "b1_s4p": _setting(rounds=6, thinning=True, squeeze="underfilled", shapes="pack_shapes"),   # S4 + packing shapes
}

# Built-in presets (a cell's pack.presets entry of the same name wins): `common_defaults` (the default when nothing is
# pinned), `legacy_v7` (nothing merged: the cell's own pack: keys with the v7 DEFAULTS), the B1 benchmark settings.
DEFAULT_PRESET = "common_defaults"


def builtin_presets() -> dict:
    return {DEFAULT_PRESET: copy.deepcopy(COMMON_DEFAULTS), "legacy_v7": {}, **copy.deepcopy(B1_BENCHMARK)}


def preset_names(preset) -> list:
    """a preset selection (None, a name, 'a,b' or a list) -> list of names ([] = nothing merged)"""
    if preset is None:
        return []
    if isinstance(preset, str):
        preset = [x.strip() for x in preset.split(",")]
    names = [str(x) for x in preset if x not in (None, "")]
    return [] if names in ([], ["none"]) else names


def _merge(a, b):
    out = copy.deepcopy(a)
    for k, v in b.items():
        if v is None:
            out.pop(k, None)
        elif isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


# Presets (A7 yeast side; B1): `pack.presets: {<name>: {pack keys}}` + `pack.preset: <name>` (or a list of names merged
# in order, or the CLI argument preset=<name>[,<name>]) selects them; their keys are merged into `pack:` (mappings
# merged key by key, a null value removes the key, other values replace). No `pack.preset` -> `common_defaults`
# (COMMON_DEFAULTS); `none` -> nothing merged. The yeast config defines `legacy_x5` (the X5 packer settings).
def effective(pc: dict, preset=None) -> dict:
    """the `pack:` section with the selected presets merged in (`presets` / `preset` removed; `_preset` = the names)"""
    pc = dict(pc or {})
    sel = preset if preset not in (None, "") else pc.get("preset", DEFAULT_PRESET)
    presets = {**builtin_presets(), **(pc.pop("presets", None) or {})}
    pc.pop("preset", None)
    names = preset_names(sel)
    for name in names:
        if name not in presets:
            near = difflib.get_close_matches(str(name), list(presets), n=1)
            raise ValueError(f"config pack.preset: '{name}' is not one of pack.presets {sorted(presets)} or 'none'"
                             + (f"; did you mean '{near[0]}'?" if near else ""))
    if not names:
        return pc
    out = pc
    for name in names:
        out = _merge(out, presets[name])
    out["_preset"] = ",".join(names)
    return out


SPLITS = ("equal_volume",)
SCOPES = ("none", "underfilled", "all")
SHAPES = ("library", "pack_shapes")


class Switches:
    def __init__(self, d):
        self.__dict__.update(d)

    @property
    def legacy(self) -> bool:
        """True when every switch has its v7 value (the packer then runs the v7 code path)"""
        return (self.rounds["n"] == 1 and not self.thinning["enabled"] and not self.shell_cap["enabled"]
                and not self.budgets["per_region"] and self.squeeze["scope"] == "none" and not self.exclusion
                and not self.organisation["enabled"] and self.shapes["source"] == "library")

    def as_dict(self):
        return {k: copy.deepcopy(getattr(self, k)) for k in DEFAULTS}


def _err(path, msg):
    raise ValueError(f"config pack.{path}: {msg} (defaults and meaning: wholecell/pack/switches.py)")


def _check_keys(path, given, allowed):
    for k in given:
        if k not in allowed:
            near = difflib.get_close_matches(k, list(allowed), n=1)
            _err(f"{path}.{k}", "unknown key" + (f"; did you mean '{near[0]}'?" if near else f"; allowed: {sorted(allowed)}"))


def _region(path, name, region_names):
    if region_names is not None and name not in region_names:
        near = difflib.get_close_matches(str(name), list(region_names), n=1)
        _err(path, f"region '{name}' not in pack.regions" + (f"; did you mean '{near[0]}'?" if near else ""))


def load(pc: dict, region_names=None) -> Switches:
    """validated switches from the `pack:` config section (missing keys: DEFAULTS)"""
    pc = pc or {}
    out = copy.deepcopy(DEFAULTS)
    for key in ("rounds", "thinning", "squeeze", "budgets"):
        v = pc.get(key)
        if v is None:
            continue
        if not isinstance(v, dict):
            _err(key, f"expected a mapping, got {type(v).__name__} {v!r}")
        _check_keys(key, v, DEFAULTS[key])
        out[key].update(copy.deepcopy(v))
    sc = pc.get("shell_cap")
    if sc is not None:
        if isinstance(sc, bool):
            out["shell_cap"]["enabled"] = sc
        elif isinstance(sc, dict):
            _check_keys("shell_cap", sc, DEFAULTS["shell_cap"])
            out["shell_cap"].update(sc)
        else:
            _err("shell_cap", f"expected true/false or a mapping, got {sc!r}")
    og = pc.get("organisation")
    if og is not None:
        if isinstance(og, bool):
            og = {"enabled": og}
        if not isinstance(og, dict):
            _err("organisation", f"expected a mapping, got {og!r}")
        _check_keys("organisation", og, DEFAULTS["organisation"])
        for sub in ("string", "genes", "subunits"):
            if og.get(sub) is not None:
                if not isinstance(og[sub], dict):
                    _err(f"organisation.{sub}", f"expected a mapping, got {og[sub]!r}")
                _check_keys(f"organisation.{sub}", og[sub], DEFAULTS["organisation"][sub])
                if sub == "string" and og[sub].get("ids") is not None:
                    _check_keys("organisation.string.ids", og[sub]["ids"], DEFAULTS["organisation"]["string"]["ids"])
                    out["organisation"]["string"]["ids"].update(copy.deepcopy(og[sub]["ids"]))
                out["organisation"][sub].update({k: copy.deepcopy(v) for k, v in og[sub].items() if not (sub == "string" and k == "ids")})
        out["organisation"].update({k: copy.deepcopy(v) for k, v in og.items() if k not in ("string", "genes", "subunits")})
        o = out["organisation"]
        if o["enabled"]:
            if not o["string"]["links"]:
                _err("organisation.string.links", "the STRING links file is required when the organisation layer is on")
            if not (o["string"]["ids"]["table"] or o["string"]["ids"]["aliases"]):
                _err("organisation.string.ids", "needs table + columns (UniProt table) or aliases (STRING aliases file)")
            if o["profile"] not in ("nodiag", "raw"):
                _err("organisation.profile", f"'{o['profile']}' not one of ['nodiag', 'raw']")
            if not 0.0 < float(o["damping"]) < 1.0:
                _err("organisation.damping", f"expected a value in (0, 1), got {o['damping']!r}")
        if not isinstance(o["field_cap_regions"], (list, tuple)):
            _err("organisation.field_cap_regions", f"expected a list of volume region names, got {o['field_cap_regions']!r}")
        o["field_cap_regions"] = [str(x) for x in o["field_cap_regions"]]
        for name in o["field_cap_regions"]:
            _region("organisation.field_cap_regions", name, region_names)
    sh = pc.get("shapes")
    if sh is not None:
        if isinstance(sh, str):
            sh = {"source": sh}
        if not isinstance(sh, dict):
            _err("shapes", f"expected library | pack_shapes or a mapping, got {sh!r}")
        _check_keys("shapes", sh, list(DEFAULTS["shapes"]) + ["library"])
        out["shapes"].update(sh)
        if out["shapes"]["source"] not in SHAPES:
            near = difflib.get_close_matches(str(out["shapes"]["source"]), SHAPES, n=1)
            _err("shapes.source", f"'{out['shapes']['source']}' not one of {list(SHAPES)}" + (f"; did you mean '{near[0]}'?" if near else ""))
    if pc.get("exclusion") is not None:
        if not isinstance(pc["exclusion"], list):
            _err("exclusion", "expected a list of {name, libs, lib_prefixes, around: {libs, lib_prefixes}, r_nm}")
        out["exclusion"] = copy.deepcopy(pc["exclusion"])
    r = out["rounds"]
    if not isinstance(r["n"], int) or isinstance(r["n"], bool) or r["n"] < 1:
        _err("rounds.n", f"expected an integer >= 1, got {r['n']!r}")
    if r["split"] not in SPLITS:
        _err("rounds.split", f"'{r['split']}' not one of {list(SPLITS)}")
    th = out["thinning"]
    th["pilot_squeeze"] = bool(th["pilot_squeeze"]); th["resplit"] = bool(th["resplit"])
    if int(th["pilot_tiles"]) < 1:
        _err("thinning.pilot_tiles", f"expected >= 1, got {th['pilot_tiles']!r}")
    if th["enabled"]:
        _region("thinning.reference_region", th["reference_region"], region_names)
    sq = out["squeeze"]
    if sq["scope"] in (None, False):
        sq["scope"] = "none"
    if sq["scope"] not in SCOPES:
        near = difflib.get_close_matches(str(sq["scope"]), SCOPES, n=1)
        _err("squeeze.scope", f"'{sq['scope']}' not one of {list(SCOPES)}" + (f"; did you mean '{near[0]}'?" if near else ""))
    if not 0.0 < float(sq["threshold"]) <= 1.0:
        _err("squeeze.threshold", f"expected a fraction in (0, 1], got {sq['threshold']!r}")
    pr = out["budgets"]["per_region"] or {}
    if not isinstance(pr, dict):
        _err("budgets.per_region", "expected {region: {max_tries, give_up}}")
    for name, b in pr.items():
        _region(f"budgets.per_region.{name}", name, region_names)
        if not isinstance(b, dict):
            _err(f"budgets.per_region.{name}", "expected {max_tries: int, give_up: int}")
        _check_keys(f"budgets.per_region.{name}", b, {"max_tries": 0, "give_up": 0})
        for k, v in b.items():
            if not isinstance(v, int) or v < 1:
                _err(f"budgets.per_region.{name}.{k}", f"expected an integer >= 1, got {v!r}")
    out["budgets"]["per_region"] = pr
    for i, e in enumerate(out["exclusion"]):
        if not isinstance(e, dict):
            _err(f"exclusion[{i}]", "expected a mapping")
        _check_keys(f"exclusion[{i}]", e, {"name": 0, "libs": 0, "lib_prefixes": 0, "around": 0, "r_nm": 0})
        if not (e.get("libs") or e.get("lib_prefixes")):
            _err(f"exclusion[{i}]", "needs libs or lib_prefixes (the packed species that are excluded)")
        a = e.get("around")
        if not isinstance(a, dict) or not (a.get("libs") or a.get("lib_prefixes") or a.get("provider")):
            _err(f"exclusion[{i}].around", "needs {libs: [...]} or {lib_prefixes: [...]} (pre-placed library ids) or "
                 "{provider: <registered exclusion_zones name>, params: {...}}")
        _check_keys(f"exclusion[{i}].around", a, {"libs": 0, "lib_prefixes": 0, "provider": 0, "params": 0})
        if a.get("provider") and (a.get("libs") or a.get("lib_prefixes")):
            _err(f"exclusion[{i}].around", "either pre-placed libraries or a provider, not both")
        if not a.get("provider") and not float(e.get("r_nm", 0)) > 0:
            _err(f"exclusion[{i}].r_nm", f"expected a radius > 0 nm, got {e.get('r_nm')!r}")
    if len(out["exclusion"]) > 62:
        _err("exclusion", "at most 62 rules (one bit each)")
    return Switches(out)
