"""Validation stage: python -m wholecell validate --cell C --build B --step S [key=value ...]

All outputs go to <build>/validate/ (write-once; pass overwrite=1 to replace a file of an unfinished run). A step that is
disabled or not applicable for the cell (config validate.steps.<step>: false, or the configuration it needs is absent)
writes <build>/validate/<step>.skipped.json with the reason, so the workflow always has an output to depend on.

Steps (inputs -> outputs, all in <build>/validate/):
  densities     pack tiles (job tables), tiles.parquet, libs.json, preplaced_libs/types, library bead masses, the build's
                RNA sources (validate.densities.mrna_sources) -> densities.json
  crowding      tiles (instances), preplaced index, mRNA index, library beads -> crowding.json, crowding.tsv
  loc_compare   pack_table.tsv, inventory_v0.tsv, tiles (job tables), libs.json, preplaced index, reference tables
                (val_* frozen inputs named in validate.loc_compare.references) -> loc_model.parquet, loc_mapping.tsv,
                loc_compare.tsv, loc_holdout.json, loc_holdout_pergene.tsv        (sub-steps: loc_model, loc_scores,
                loc_holdout; 'loc_compare' runs all three). Skipped without validate.loc_compare.
  ribosome_nn   polysome pre-pass (the build's own), geometry centre, EM ribosome table (val_ribo_em) -> ribosome_nn.json
                Skipped without validate.ribosome_nn.
  placement     seeded tile sample: tiles, preplaced + mRNA index, library sphere trees, geometry, region table
                -> placement.json, placement_tiles.tsv
  chromatin     frozen chromatin pre-pass (md5 checked) + frozen chromatin metric files -> chromatin_frozen.json
                Skipped without validate.chromatin.
  manifest      every result file in validate/ -> MANIFEST.tsv (file, bytes, md5, build stamp id, validation commit)
  compare_v7    reference (validate.compare_reference: legacy paths + its validate/ dir) vs this build
                -> compare_v7/criteria_results.tsv, compare_v7/verdict.json  (wholecell/validate/compare_builds.py)
  suite         common validation suite (plan s4, wholecell/validate/suite.py): checks C1-C17 from the outputs above,
                the build's tiles and the cell's reference adapter (validate.suite.adapter) -> suite.json, no_regression.tsv,
                legacy_check.tsv (in validate/, or --out DIR; runs after manifest, so the manifest does not list them)
  code_check    (reference runs only) this run's outputs vs the stored legacy validation files
                -> code_check_vs_stored.tsv
  all           densities crowding loc_compare ribosome_nn placement chromatin manifest
Common arguments: pack_dir=, tiles_dir=, ref=<reference name> (run on the archived reference pack, files from
validate.compare_reference.files), geom=, lib_dir=, regions_json=, pack_table=, inventory_v0=, nproc= (<= 8),
overwrite=1. Step arguments: crowding n_boxes=; placement n_tiles=.

NOT ported in phase 1 (stated, criteria V8): Tier-2 imaging (tier2.py) and the fluorescence panels (fluor.py), the
galleries, the HTML page, and the loc_compare reference-building steps (fetch/celllines/jurkat) and rule-development
steps (rules/report).
"""
from __future__ import annotations

import os

os.environ.setdefault("NUMBA_NUM_THREADS", "8")
os.environ.setdefault("NUMBA_THREADING_LAYER", "workqueue")
os.environ.setdefault("NUMBA_CACHE_DIR", os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), ".numba_cache"))   # never into a legacy __pycache__

import json  # noqa: E402

from . import common  # noqa: E402

STEPS_ALL = ["densities", "crowding", "loc_compare", "ribosome_nn", "placement", "chromatin", "manifest"]


def _skip(bdir, step, reason, overwrite=True):
    p = os.path.join(common.outdir(bdir), f"{step}.skipped.json")
    common.write_json(p, dict(step=step, skipped=True, reason=reason), overwrite=True)
    print(f"validate {step}: skipped ({reason})")
    return None


def _enabled(cfg, step):
    st = ((cfg.get("validate", {}) or {}).get("steps", {}) or {})
    return st.get(step, True) is not False


def run(cfg, bdir, step=None, **kw):
    step = step or "all"
    if step == "all":
        for s in STEPS_ALL:
            run(cfg, bdir, s, **kw)
        return
    if step == "suite":                 # before outdir(): with --out the suite writes nothing into the build
        from . import suite
        return suite.run(cfg, bdir, out=kw.pop("out", None) or None, **kw)
    ow = common.flag(kw, "overwrite")
    vc = cfg.get("validate", {}) or {}
    np_ = common.nprocs(cfg, kw)
    out = common.outdir(bdir)
    if not _enabled(cfg, step):
        return _skip(bdir, step, f"disabled by config validate.steps.{step}")
    inp = common.Inputs(cfg, bdir, **kw)
    if step == "densities":
        from . import densities
        return densities.run(cfg, bdir, inp, os.path.join(out, "densities.json"), ow, nproc=np_)
    if step == "crowding":
        from . import crowding
        return crowding.run(cfg, bdir, inp, os.path.join(out, "crowding.json"), os.path.join(out, "crowding.tsv"), ow,
                            n_boxes=kw.get("n_boxes"))
    if step in ("loc_compare", "loc_model", "loc_scores", "loc_holdout"):
        if not vc.get("loc_compare"):
            return _skip(bdir, step, "no validate.loc_compare configuration (no localisation references for this cell)")
        from . import loc_compare as LC
        C = LC.Ctx(cfg, bdir, inp)
        if step in ("loc_compare", "loc_model"):
            LC.model(C, ow, nproc=np_)
        if step in ("loc_compare", "loc_scores"):
            LC.compare(C, ow)
        if step in ("loc_compare", "loc_holdout"):
            LC.holdout(C, ow)
        return
    if step == "ribosome_nn":
        if not vc.get("ribosome_nn"):
            return _skip(bdir, step, "no validate.ribosome_nn configuration (no EM ribosome reference)")
        from . import ribosome_nn
        from .geom import open_geometry
        c = kw.get("centre_nm")
        if c is None:
            geo = open_geometry(cfg, inp, bdir)
            c = geo.centre_nm
        else:
            c = [float(x) for x in str(c).split(",")]
        import numpy as np
        return ribosome_nn.run(cfg, bdir, inp, np.asarray(c, float), os.path.join(out, "ribosome_nn.json"), ow)
    if step == "placement":
        from . import placement
        return placement.run(cfg, bdir, inp, os.path.join(out, "placement.json"), os.path.join(out, "placement_tiles.tsv"), ow,
                             nproc=np_, n_tiles=kw.get("n_tiles"))
    if step == "chromatin":
        if not vc.get("chromatin"):
            return _skip(bdir, step, "no validate.chromatin configuration (no frozen chromatin pre-pass)")
        from . import chromatin
        return chromatin.run(cfg, bdir, os.path.join(out, "chromatin_frozen.json"), ow)
    if step == "manifest":
        return write_manifest(cfg, bdir, overwrite=True)
    if step == "compare_v7":
        from . import compare_builds
        return compare_builds.run(cfg, bdir, **kw)
    if step == "code_check":
        from . import code_check
        return code_check.run(cfg, bdir, **kw)
    raise ValueError(f"unknown validate step '{step}' (steps: {STEPS_ALL + ['compare_v7', 'suite', 'code_check', 'all']})")


def write_manifest(cfg, bdir, overwrite=True):
    """MANIFEST.tsv: every file under validate/ (except cache/ and the manifest) with bytes, md5, build stamp id"""
    import hashlib
    from ..build import git_info
    out = common.outdir(bdir)
    sp = os.path.join(bdir, "BUILD_STAMP.json")
    if os.path.exists(sp):
        stamp_id = "stamp:" + hashlib.sha256(open(sp, "rb").read()).hexdigest()[:16]
        st = json.load(open(sp))
        stamp_id += f" ({st.get('commit7', '')}{'-dirty' if st.get('dirty') else ''})"
    else:
        stamp_id = "no-build-stamp (reference run)"
    gi = git_info()
    rows = []
    for root, dirs, files in os.walk(out):
        dirs[:] = sorted(d for d in dirs if d not in ("cache",))
        for f in sorted(files):
            if f == "MANIFEST.tsv" or ".tmp" in f:
                continue
            p = os.path.join(root, f)
            rows.append((os.path.relpath(p, out), os.path.getsize(p), common.md5(p)))
    p = os.path.join(out, "MANIFEST.tsv")
    tmp = p + f".tmp{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write(f"# build_stamp\t{stamp_id}\n# validation_code\t{gi['commit7']}{'-dirty' if gi['dirty'] else ''}\n")
        fh.write("file\tbytes\tmd5\tbuild_stamp\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{stamp_id}\n")
    os.replace(tmp, p)
    print(f"MANIFEST.tsv: {len(rows)} files")
    return p
