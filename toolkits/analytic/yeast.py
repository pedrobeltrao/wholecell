"""Yeast analytic geometry generator -> ONE cellgeom_v1 package (engine yeast_cell3d_v1; spec docs/reference/cellgeom_v1.md).

Replaces the stub `yeast_stub.py` (kept for the primitive-table layout constants). The yeast model hand-copied its
geometry in four places (pipeline/cell3d.py, web geometry.ts, the WGSL shaders, the viewer section map) and rebuilt it
inside stage 2; here the generator runs ONCE per build, writes an immutable package, and every consumer of the build
(stage-2 tiling, packer, validation, encoder) reads that package (engine PORT_PATCHES P1). Viewers are meant to read
`prims.npz` + `fields/` + the manifest (docs/reference/viewer_multi_cell.md table-driven evaluator) instead of hand copies.

Generator = the yeast stage-2 builder calls, unchanged (vendored code, run in the engine environment):
    df, geo = inventory.build()                                   # volumes/areas from the inventory (stage 1)
    shapes  = cellshape.build(geo, np.random.default_rng(11))     # real outline / nucleus / vacuole (rng 11 fixed)
    cell    = cell3d.build(geo, seed=11, shapes=shapes, n_endo=round(patch actin / 1500), n_scars=BUD_SCARS)

    python toolkits/analytic/yeast.py --engine-code <plugins/yeast/engine> --dest-root <geom_root>/yeast --out <build>/geometry/package.json
(environment of the yeast engine: YC_WORK etc., set by plugins/yeast/runner.py)

Package layout (package_id yeast/analytic_x5_s<seed>_<sha12 of native/geometry.json>):
    manifest.json, validation.json, native/{geometry.json, geometry_fields.npz} (engine-native, read by cell3d.Cell.load),
    fields/{cell,vacuole,nucleus}.npz (§3 sparse bricks), prims.npz (§6.2), pores.npz, summary (MC label volumes, rng 0).
"""
import argparse
import datetime as _dt
import hashlib
import json
import os
import shutil
import stat
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, REPO)
from wholecell.geometry import cellgeom  # noqa: E402  (numpy-only helpers: savez_det, write_field_npz, sha256, core labels)
from toolkits.analytic.yeast_stub import PRIM_TYPES, OPS, YEAST_CLASSES  # noqa: E402

GEN_VERSION = "yeast_analytic_v1"
# yeast cell3d label code -> (vocabulary name, core id or extension id, parent core id, GO)
LABEL_MAP = [
    ("exterior", 0, None), ("wall_mannan", 1000, 15), ("wall_glucan", 1001, 15), ("cytosol", 1, None),
    ("er_lumen", 5, None), ("perinuclear_space", 4, None), ("nucleoplasm", 2, None), ("nucleolus", 3, None),
    ("vacuole_lumen", 12, None), ("mito_ims", 7, None), ("mito_matrix", 9, None), ("lipid_droplet", 14, None),
    ("golgi_endosome_lumen", 1002, 10), ("peroxisome_matrix", 13, None), ("nuclear_pore_channel", 6, None)]
EXT_DESC = {1000: ("outer mannoprotein layer of the cell wall", "GO:0009277"),
            1001: ("inner beta-glucan/chitin layer of the cell wall", "GO:0009277"),
            1002: ("lumen of Golgi cisternae and endosomes (one yeast class)", None)}
# yeast membrane code -> (field / derivation, n_side, p_side)
MEM_SIDES = {"PM": ("cell", "cytosol", "cell wall"), "ER": ("prims:cortical_er + NE", "cytosol", "er_lumen"),
             "NE": ("nucleus", "cytosol / nucleoplasm", "perinuclear_space"), "VAC": ("vacuole", "cytosol", "vacuole_lumen"),
             "MITO_OM": ("prims:mito_om", "cytosol", "mito_ims"), "MITO_IM": ("prims:mito_im + crista", "mito_matrix", "mito_ims"),
             "GOLGI": ("prims:golgi", "cytosol", "golgi_endosome_lumen"), "PEROX": ("prims:peroxisome", "cytosol", "peroxisome_matrix"),
             "LD": ("prims:lipid_droplet (monolayer)", "cytosol", "lipid_droplet"), "CELL_OUTER": ("cell (outer wall surface)", "exterior", "cell wall")}


def log(*a, fh=None):
    s = time.strftime("%H:%M:%S") + " " + " ".join(str(x) for x in a)
    print(s, flush=True)
    if fh:
        fh.write(s + "\n"); fh.flush()


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def build_native(engine_code, logf):
    """the yeast stage-2 builder calls (vendored code), returning (cell, geo, n_endo)"""
    sys.path.insert(0, os.path.join(engine_code, "pipeline"))
    sys.path.insert(0, os.path.join(engine_code, "scripts"))
    import _common as K  # noqa: F401  (paths: F_ACTIN_JSON; creates the inventory folder as the stage scripts do)
    import inventory as inv, cell3d as c3, cellshape, actin_patches
    L = lambda *a: log(*a, fh=logf)
    df, geo = inv.build()
    shapes = cellshape.build(geo, np.random.default_rng(11), log=L)
    fa = json.load(open(K.F_ACTIN_JSON))
    n_endo = int(round(fa["patches"] / actin_patches.ACTIN_PER_PATCH))
    cell = c3.build(geo, seed=11, log=L, shapes=shapes, n_endo=n_endo, n_scars=c3.BUD_SCARS)
    return cell, geo, n_endo, c3


def mc_volumes(cell, c3):
    """the stage-2 Monte-Carlo check (4M points in the cell sphere, rng 0), per label (um3)"""
    rng = np.random.default_rng(0); R0 = cell.G[0]
    P = rng.uniform(-R0, R0, (4_000_000, 3)); P = P[np.linalg.norm(P, axis=1) < R0]
    lab = c3.label_points(P, cell.G, cell.pores, np.arange(len(cell.PT)), cell.PT, cell.PQ, cell.PC)
    vol = np.bincount(lab, minlength=15) / len(P) * 4 / 3 * np.pi * (R0 / 1000) ** 3
    return {c3.LABEL_NAMES[i]: float(v) for i, v in enumerate(vol)}


def write_package(cell, geo, n_endo, c3, dest_root, engine_code, logf):
    tmp = os.path.join(dest_root, f".tmp-{os.getpid()}")
    if os.path.exists(tmp):
        shutil.rmtree(tmp)
    os.makedirs(os.path.join(tmp, "native")); os.makedirs(os.path.join(tmp, "fields"))
    # (a) engine-native files: exactly what stage 2 wrote before the port (cell.to_json)
    cell.to_json(os.path.join(tmp, "native", "geometry.json"))
    ff = os.path.join(tmp, "native", "geometry_fields.npz")      # re-save without zip timestamps (content-addressed package)
    if os.path.exists(ff):
        z = dict(np.load(ff)); os.remove(ff); cellgeom.savez_det(ff, **z)
    native_sha = cellgeom.sha256_file(os.path.join(tmp, "native", "geometry.json"))
    # (b) shared layout: fields (§3), prims (§6.2), pores
    fields_meta = {}
    names = ["cell", "vacuole", "nucleus"]
    for k, (h, idx, data) in enumerate(cell.fields or []):
        rel = f"fields/{names[k]}.npz"
        dt = cellgeom.write_field_npz(os.path.join(tmp, rel), h, idx, data)
        h = np.asarray(h, float)
        fields_meta[names[k]] = dict(slot=k, file=rel, dtype=dt, h_nm=float(h[3]), brick=int(h[7]), origin_nm=h[:3].tolist(),
                                     bricks=[int(x) for x in h[4:7]], clamp_in_nm=float(h[8]), clamp_out_nm=float(h[9]),
                                     bricks_stored=int((np.asarray(idx) >= 0).sum()), sign="negative inside",
                                     provenance="yeast cellshape.build (EMPIAR-11416 meshes, cell_meshes.npz), rng 11")
    for k in range(len(cell.fields or []), 3):
        fields_meta[names[k]] = dict(slot=k, empty=True)
    smooth_k = {"pm_invagination": float(c3.EIS_K), "bud_scar": float(c3.SCAR_K)}
    op = np.zeros(len(cell.PT), np.int8)
    op[np.isin(cell.PC, [c3.C_PMINV, c3.C_SCAR])] = OPS["smooth_union"]
    cellgeom.savez_det(os.path.join(tmp, "prims.npz"), type=cell.PT.astype(np.int16), params=cell.PQ.astype(np.float64),
                       comp=cell.PC.astype(np.int16), field_bind=np.full(len(cell.PT), -1, np.int16), op=op,
                       bsphere=np.c_[cell.bc, cell.br].astype(np.float64))
    cellgeom.savez_det(os.path.join(tmp, "pores.npz"), pores=cell.pores.astype(np.float64))
    # manifest
    labels = []
    for code, (name, lid, parent) in enumerate(LABEL_MAP):
        e = dict(id=lid, name=name, engine_code=code, engine_name=c3.LABEL_NAMES[code])
        if lid >= cellgeom.EXT_ID_MIN:
            e.update(parent=parent, description=EXT_DESC[lid][0], go=EXT_DESC[lid][1])
        else:
            core = {i: (n, g, d) for i, n, g, d in cellgeom.CORE_LABELS}[lid]
            e.update(go=core[1], description=core[2])
        labels.append(e)
    labels += [dict(id=i, name=n, go=g, description=d, engine_code=None) for i, n, g, d in cellgeom.CORE_LABELS
               if i not in {l["id"] for l in labels}]
    mems = [dict(name=m, engine_code=i, field=MEM_SIDES[m][0], derivation="cell3d.mem_sdfs", offset_nm=None,
                 half_thickness_nm=float(c3.HALF), n_side=MEM_SIDES[m][1], p_side=MEM_SIDES[m][2],
                 leaflets="monolayer" if m == "LD" else ("none" if m == "CELL_OUTER" else "bilayer"), go=None)
            for i, m in enumerate(c3.MEM_NAMES)]
    params = {k: (float(getattr(c3, k)) if np.isscalar(getattr(c3, k)) else list(np.atleast_1d(getattr(c3, k)).astype(float)))
              for k in ("HALF", "PORE_R", "PERI_W", "IMS_W", "CER_GAP", "CER_LUMEN", "EIS_K", "ENDO_R", "SCAR_R", "SCAR_r", "SCAR_K", "G_N")}
    params["MEM_HALF"] = [float(c3.HALF)] * len(c3.MEM_NAMES)
    vols = mc_volumes(cell, c3)
    man = dict(spec_version=cellgeom.SPEC_VERSION, cell_id="yeast", geometry_id=None, description=
               "S. cerevisiae BY4741 G1 cell: real outline/nucleus/vacuole SDFs (EMPIAR-11416) + analytic organelle primitives (yeast build X5 generator)",
               organism=dict(name="Saccharomyces cerevisiae", taxon=559292), cell_line="BY4741 (haploid, MATa)",
               source=dict(dataset="EMPIAR-11416 segmentations (Dumoux 2023) via yeast pipeline/cell_meshes.npz + organelle_stats.json",
                           licence="yeast code: no licence yet (intended CC BY-NC 4.0); local only"),
               generator=dict(toolkit="analytic", module="toolkits/analytic/yeast.py", version=GEN_VERSION,
                              yeast_code="plugins/yeast/engine @ e3b7448 (vendored)", seeds=dict(cellshape=11, cell3d="config.seed(11)", mc_volumes=0),
                              YC_SEED=int(os.environ.get("YC_SEED", "0")), n_endocytic_sites=n_endo, n_bud_scars=int(c3.BUD_SCARS),
                              inputs=dict(geo={k: (float(v) if np.isscalar(v) else str(v)) for k, v in geo.items() if np.isscalar(v)})),
               engine=dict(name="yeast_cell3d_v1", module="plugins.yeast.geom_engine", layout_id="cell3d_G23_fields3",
                           params=params, field_slots={n: fields_meta[n]["slot"] for n in names}),
               frame=dict(units="nm", handedness="right", origin="cell centre (yeast cell3d frame)", source_affines={}),
               domain=dict(enclosing_radius_nm=float(cell.G[0]), clip_planes=[]),
               fields=fields_meta, labels=labels, membranes=mems,
               label_rule=None,
               validation_waivers={"label_rule": "evaluator-defined (cell3d.label_from); the declarative rule for the yeast engine is not written yet (PORTING_STATUS)"},
               sections={},
               prims=dict(file="prims.npz", n=int(len(cell.PT)), types={k: v for k, v in PRIM_TYPES.items() if v < 6},
                          classes=YEAST_CLASSES, ops=OPS, smooth_k=smooth_k,
                          note="params = yeast PQ rows (8 floats per primitive, layout per type in toolkits/analytic/yeast_stub.PRIM_PARAMS)"),
               pores=dict(file="pores.npz", n=int(len(cell.pores)), columns=["cx", "cy", "cz", "nx", "ny", "nz"],
                          source_codes={"random": 0}, counts={"random": int(len(cell.pores))}, min_sep_nm=150.0,
                          provenance="cell3d.build: centres on the perinuclear midsurface, outside the nucleus-vacuole junction"),
               layers={},
               native=dict(files=["native/geometry.json", "native/geometry_fields.npz"], loader="cell3d.Cell.load (plugins.yeast.geom_engine.load)",
                           geometry_json_sha256=native_sha),
               summary=dict(mc_label_volumes_um3=vols, open_face=dict(f_vol=1.0, f_pm=1.0)))
    gid = f"analytic_x5_s{int(os.environ.get('YC_SEED', '0'))}_{native_sha[:12]}"
    man["geometry_id"] = gid; man["package_id"] = f"yeast/{gid}"
    man["created"] = _dt.datetime.now().isoformat(timespec="seconds")
    files = {}
    for root, _, fs in os.walk(tmp):
        for f in sorted(fs):
            p = os.path.join(root, f); rel = os.path.relpath(p, tmp)
            files[rel] = dict(sha256=cellgeom.sha256_file(p), bytes=os.path.getsize(p))
    man["files"] = dict(sorted(files.items()))
    json.dump(man, open(os.path.join(tmp, "manifest.json"), "w"), indent=1, default=float)
    res = validate(tmp, cell, c3, man)
    res["manifest_sha256"] = cellgeom.sha256_file(os.path.join(tmp, "manifest.json"))
    json.dump(res, open(os.path.join(tmp, "validation.json"), "w"), indent=1, default=float)
    if not res["pass"]:
        raise SystemExit(f"package validation failed; left in {tmp}: {[k for k, v in res['checks'].items() if not v['pass']]}")
    dest = os.path.join(dest_root, gid)
    if os.path.exists(dest):
        # immutable: an existing package with this id must have identical content files (manifest timestamps aside)
        old = json.load(open(os.path.join(dest, "manifest.json")))
        if old["files"] != man["files"]:
            raise SystemExit(f"{dest} exists with different content (immutable packages)")
        shutil.rmtree(tmp)
        log("package exists with identical content files:", dest, fh=logf)
    else:
        os.rename(tmp, dest)
        for root, _, fs in os.walk(dest):
            for f in fs:
                p = os.path.join(root, f)
                os.chmod(p, os.stat(p).st_mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))
        log("package written:", dest, fh=logf)
    return dest, man, res


def validate(pkg, cell, c3, man):
    """checks: checksums, schema, native round trip, prims/fields/pores vs native, MC volumes recorded"""
    chk = {}
    bad = [rel for rel, m in man["files"].items() if cellgeom.sha256_file(os.path.join(pkg, rel)) != m["sha256"]]
    chk["checksums"] = dict(pass_=not bad, bad=bad)
    ids = [l["id"] for l in man["labels"]]; codes = [l["engine_code"] for l in man["labels"] if l.get("engine_code") is not None]
    core = {i: n for i, n, _, _ in cellgeom.CORE_LABELS}
    errs = []
    if len(set(ids)) != len(ids): errs.append("duplicate label ids")
    if len(set(codes)) != len(codes): errs.append("duplicate engine codes")
    if any(c < 0 or c > cellgeom.MAX_ENGINE_CODE for c in codes): errs.append("engine code outside 0..63")
    for l in man["labels"]:
        if l["id"] < cellgeom.EXT_ID_MIN and core.get(l["id"]) != l["name"]: errs.append(f"core id {l['id']} name {l['name']}")
        if l["id"] >= cellgeom.EXT_ID_MIN and l.get("parent") not in core: errs.append(f"extension {l['name']} without core parent")
    if sorted(codes) != list(range(15)): errs.append("engine codes != cell3d labels 0..14")
    chk["schema"] = dict(pass_=not errs, errors=errs)
    L = c3.Cell.load(os.path.join(pkg, "native", "geometry.json"))
    rt = dict(G=bool(np.array_equal(L.G, cell.G)), pores=bool(np.array_equal(L.pores, cell.pores)), PT=bool(np.array_equal(L.PT, cell.PT)),
              PQ=bool(np.array_equal(L.PQ, cell.PQ)), PC=bool(np.array_equal(L.PC, cell.PC)), bc=bool(np.array_equal(L.bc, cell.bc)),
              br=bool(np.array_equal(L.br, cell.br)))
    chk["native_round_trip"] = dict(pass_=all(rt.values()), **rt)
    z = np.load(os.path.join(pkg, "prims.npz"))
    pr = dict(type=bool(np.array_equal(z["type"].astype(np.int64), L.PT)), params=bool(np.array_equal(z["params"], L.PQ)),
              comp=bool(np.array_equal(z["comp"].astype(np.int64), L.PC)), bsphere=bool(np.array_equal(z["bsphere"], np.c_[L.bc, L.br])))
    chk["prims_vs_native"] = dict(pass_=all(pr.values()), **pr)
    fl = {}
    names = ["cell", "vacuole", "nucleus"]
    for k, (h, idx, data) in enumerate(L.fields or []):
        hh, ii, dd = cellgeom.read_field_npz(os.path.join(pkg, f"fields/{names[k]}.npz"))
        fl[names[k]] = bool(np.array_equal(hh, np.asarray(h, float)) and np.array_equal(ii, np.asarray(idx, float))
                            and np.array_equal(dd, np.asarray(data, float)))
    G2 = c3.pack_fields(np.asarray(L.d["G"], float), [cellgeom.read_field_npz(os.path.join(pkg, f"fields/{n}.npz")) for n in names[:len(L.fields or [])]]) if L.fields else None
    fl["G_from_shared_fields"] = bool(G2 is None or np.array_equal(G2, L.G))
    chk["fields_vs_native"] = dict(pass_=all(fl.values()), **fl)
    pz = np.load(os.path.join(pkg, "pores.npz"))["pores"]
    chk["pores_vs_native"] = dict(pass_=bool(np.array_equal(pz, L.pores)), n=int(len(pz)))
    chk["mc_label_volumes"] = dict(pass_=True, um3=man["summary"]["mc_label_volumes_um3"], note="reported (stage-2 MC check, rng 0, 4M points)")
    checks = {k: {("pass" if kk == "pass_" else kk): vv for kk, vv in v.items()} for k, v in chk.items()}
    return dict(spec_version=cellgeom.SPEC_VERSION, package_id=man["package_id"], engine=man["engine"]["name"], checks=checks,
                waived=list(man.get("validation_waivers", {})), **{"pass": all(c["pass"] for c in checks.values())})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine-code", required=True)
    ap.add_argument("--dest-root", required=True)
    ap.add_argument("--out", required=True, help="<build>/geometry/package.json")
    ap.add_argument("--log", default=None)
    a = ap.parse_args()
    os.makedirs(a.dest_root, exist_ok=True)
    logf = open(a.log, "w") if a.log else None
    t0 = time.time()
    cell, geo, n_endo, c3 = build_native(a.engine_code, logf)
    log("prims", len(cell.PT), "pores", len(cell.pores), fh=logf)
    dest, man, res = write_package(cell, geo, n_endo, c3, a.dest_root, a.engine_code, logf)
    out = dict(package_dir=dest, package_id=man["package_id"], engine=man["engine"]["name"],
               manifest_sha256=cellgeom.sha256_file(os.path.join(dest, "manifest.json")),
               native_geometry_json_sha256=man["native"]["geometry_json_sha256"], validation_pass=res["pass"],
               linked_by="toolkits/analytic/yeast.py", seconds=round(time.time() - t0, 1))
    json.dump(out, open(a.out, "w"), indent=1)
    log(json.dumps(out), fh=logf)


if __name__ == "__main__":
    main()
