# new in wholecell (B4, 2026-11-17)
"""B4 common pre-placed builders: measurements and verdict (criteria docs/unification_criteria/B4_preplaced.md).

  measure   Jurkat pre-placed layers of a serial run (pre-build check or a build): J1/J1b pores, J2 MT dimers, J3 polysome
            fraction (a: chains of >= 2 / RIBO80S pool; b: EM-comparable linkage statistic on the ribosome centres thinned
            to the de-duplicated EM count), J4 ribosome NN median thinned to the EM count, J5b pre-pass self-overlaps by
            module pair on the reference's sampled tiles (new layers and the v7 layers, same code and trees), J6 mito.
            python scripts/b4_preplaced_check.py measure --serial <dir> --build <bdir: library, inventory> \
                   --ref-build <build with pack/tiles.parquet> --out <json> [--packed <ribosome positions .npy>]
  (addendum 2026-11-17: J5b also per pair class per 1e6 instances of the class's B4 modules, with absolute pair counts;
   the verdict reports the frozen rule set (verdict_frozen) and the addendum rule set (verdict); see
   docs/unification_criteria/B4_preplaced_addendum_2026-11-17.md)
  verdict   compare_B4/verdict.json of a build: measure on <build>/preplaced/serial (+ packed RIBO80S of the tiles),
            J5a from validate/placement.json vs the reference, suite no_regression vs --ref-suite with the EXPECTED
            exemptions; yeast: Y1 arrays of preplaced/serial vs the PATCH-4 reference, Y2-Y6 from the suite.
            python scripts/b4_preplaced_check.py verdict --cell jurkat|yeast --build <bdir> --ref-suite <suite.json> \
                   [--ref-build <bdir>] [--complete 1] [--extra-expected b6,b3]
"""
import argparse
import json
import os
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
from _devsite import VC as V, RESULTS  # noqa: E402,F401  (release: site paths)
EM_N, EM_DEDUP_N, EM_LINKED, EM_NN = 730960, 682785, 0.697, 27.05
V7 = dict(mt_dimers=2325171, nn_thinned=28.24, prepass_prepass_per_1e6=838914.7117814663, in_chains=0.9454945678979821)
YREF = dict(pores=93, nn=28.392134049682284, prepass_prepass_per_1e6=642511.1021396851, in_chains=0.7843018329010251)
YREF_DIR = f"{RESULTS}/yeast/_dev_unify/amphiphysin_fix/ref_yeastenv_patched/prepass/serial"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def linked_fraction(P, r=34.0):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree
    pr = cKDTree(P).query_pairs(r, output_type="ndarray")
    A = coo_matrix((np.ones(len(pr)), (pr[:, 0], pr[:, 1])), shape=(len(P), len(P)))
    _, lab = connected_components(A, directed=False)
    sz = np.bincount(lab)
    return float((sz[lab] >= 2).mean())


def nn_median(P):
    from scipy.spatial import cKDTree
    d, _ = cKDTree(P).query(P, k=2)
    return float(np.median(d[:, 1]))


def ribosome_stats(PM, seed=0):
    rng = np.random.default_rng(seed)
    thin = PM[rng.random(len(PM)) < EM_N / len(PM)]
    rng2 = np.random.default_rng(seed + 1)
    thin2 = PM[rng2.random(len(PM)) < EM_DEDUP_N / len(PM)]
    return dict(n=int(len(PM)), nn_median_all=round(nn_median(PM), 3), nn_median_thinned_to_EM=round(nn_median(thin), 3),
                linked_34nm_thinned_to_dedup_EM=round(linked_fraction(thin2), 4), em_nn=EM_NN, em_linked=EM_LINKED)


def _lib_trees(build):
    import pandas as pd
    L = pd.read_csv(os.path.join(build, "library", "library_manifest.tsv"), sep="\t")
    L = L[L.path.str.endswith(".npz")]
    return dict(zip(L.species, L.path))


def _layers_new(serial, cfg, lamina=None):
    """B4 layers of a serial run + chromatin rep5 + the lamina (frozen v7 copy, or `lamina`: e.g. the rebuilt lattice of
    mode build, <bdir>/preplaced/lamina/build/lamina_cell1.npz)"""
    from wholecell.preplaced.stage import _resolve
    out = {}
    for f in sorted(os.listdir(os.path.join(serial, "modules"))):
        if f.endswith(".npz"):
            z = np.load(os.path.join(serial, "modules", f))
            if "types" in z.files:
                out[f[:-4]] = (z["lib_ids"].astype(str), z["types"], z["pos"], z["quat"])
    for name in ("pre_chromatin_v8_rep5", "pre_lamina_v7"):
        z = np.load(lamina if (lamina and name == "pre_lamina_v7") else _resolve(cfg, None, name))
        out[name.replace("pre_", "").replace("_v8_rep5", "").replace("_v7", "")] = (z["lib_ids"].astype(str), z["types"], z["pos"], z["quat"])
    return out


def _layers_v7(cfg):
    from wholecell.preplaced.stage import _resolve
    out = {}
    for name, key in (("pre_chromatin_v8_rep5", "chromatin"), ("pre_polysomes_v7", "polysomes"), ("pre_cytoskeleton_v7", "cytoskeleton"),
                      ("pre_npc_v3", "npc"), ("pre_lamina_v7", "lamina")):
        z = np.load(_resolve(cfg, None, name))
        out[key] = (z["lib_ids"].astype(str), z["types"], z["pos"], z["quat"])
        if key == "cytoskeleton":                            # split into MT lattice (kind 1) and the rest (centrioles, actin)
            k = z["kind"]
            out["microtubules"] = (z["lib_ids"].astype(str), z["types"][k == 1], z["pos"][k == 1], z["quat"][k == 1])
            out["cytoskeleton_rest"] = (z["lib_ids"].astype(str), z["types"][k != 1], z["pos"][k != 1], z["quat"][k != 1])
            del out["cytoskeleton"]
    return out


def self_overlaps(layers, LIB, ref_build, cfg, b4_modules, n_tiles=300, tol=0.5, h=6.0):
    """pre-pass-pre-pass overlapping instance pairs on the sampled tiles (placement C7a definition: library sphere trees,
    overlap > tol, instances bucketed by centre tile, pairs within a tile), by module pair"""
    import pandas as pd
    from wholecell.validate import common
    from wholecell.validate.placement import _clash, _spheres
    M = pd.read_parquet(os.path.join(ref_build, "pack", "tiles.parquet"), columns=["tile", "ti", "tj", "tk"])
    meta = json.load(open(os.path.join(ref_build, "pack", "tiles_meta.json")))
    lo0, T = np.asarray(meta["lo0"]), float(meta["TILE"])
    tiles = common.sample_tiles(cfg, len(M), n_tiles)
    tij = M[["ti", "tj", "tk"]].values[tiles]
    want = {tuple(x): k for k, x in enumerate(tij)}
    names = list(layers)
    libs = sorted({str(i) for n in names for i in layers[n][0]})
    li = {s: k for k, s in enumerate(libs)}
    off = [0]; C, R = [], []
    for s in libs:
        z = np.load(LIB[s]); C.append(z["st_c"].astype(np.float64)); R.append(z["st_r"].astype(np.float64)); off.append(off[-1] + len(R[-1]))
    st_off, STC, STR = np.array(off, np.int64), np.concatenate(C), np.concatenate(R)
    per_tile = [[] for _ in tiles]
    for mi, n in enumerate(names):
        ids, t, P, Q = layers[n]
        g = np.array([li[str(x)] for x in ids], np.int64)[t]
        tt = np.floor((np.asarray(P, np.float64) - lo0) / T).astype(np.int64)
        key = (tt[:, 0] * 100003 + tt[:, 1]) * 100003 + tt[:, 2]
        kk = (tij[:, 0] * 100003 + tij[:, 1]) * 100003 + tij[:, 2]
        hit = np.isin(key, kk)
        idx = np.where(hit)[0]
        pos_k = {int(v): k for k, v in enumerate(kk)}
        for i in idx:
            per_tile[pos_k[int(key[i])]].append((mi, g[i], P[i], Q[i]))
    pairs = {}
    counts = {n: 0 for n in names}
    for k in range(len(tiles)):
        rows = per_tile[k]
        if len(rows) < 2:
            for r in rows: counts[names[r[0]]] += 1
            continue
        mod = np.array([r[0] for r in rows]); ty = np.array([r[1] for r in rows], np.int64)
        P = np.array([r[2] for r in rows], np.float64); Q = np.array([r[3] for r in rows], np.float64)
        for m in mod: counts[names[m]] += 1
        Cs, Rs, Is = _spheres(ty, P, Q, st_off, STC, STR, 0)
        cat = np.ones(len(rows), np.int8)
        a, b, d, tot = _clash(Cs, Rs, Is, cat, h, tol, 20_000_000)
        if len(a):
            key = np.unique(np.minimum(a, b) * (len(rows) + 1) + np.maximum(a, b))
            ia, ib = key // (len(rows) + 1), key % (len(rows) + 1)
            for x, y in zip(mod[ia], mod[ib]):
                c = "|".join(sorted((names[x], names[y])))
                pairs[c] = pairs.get(c, 0) + 1
    n_all = sum(counts.values()); n_b4 = sum(counts[n] for n in names if n in b4_modules)
    tot = sum(pairs.values())
    b4p = {c: v for c, v in pairs.items() if any(m in b4_modules for m in c.split("|"))}
    # B4 addendum 2026-11-17 (J5b per pair class): rate of a class = its pairs per 1e6 instances of the B4 module(s) of
    # that class on the sampled tiles (MT-MT per 1e6 MT dimers, NPC-chromatin per 1e6 NPCs, NPC-polysome per 1e6 NPCs +
    # ribosomes), so a class does not move when another B4 module changes its instance count
    cls_rate = {}
    for c, v in pairs.items():
        ms = sorted(set(m for m in c.split("|") if m in b4_modules))
        if ms:
            den = sum(counts.get(m, 0) for m in ms)
            cls_rate[c] = dict(pairs=int(v), den_modules=ms, den=int(den), per_1e6=round(v / max(den, 1) * 1e6, 3))
    return dict(n_tiles=int(len(tiles)), instances=counts, pairs=pairs, prepass_prepass_per_1e6=round(tot / max(n_all, 1) * 1e6, 3),
                prepass_prepass_pairs=int(tot), class_rates_addendum=cls_rate,
                b4_instances=n_b4, b4_pair_classes_per_1e6_b4={c: round(v / max(n_b4, 1) * 1e6, 3) for c, v in b4p.items()},
                b4_pairs_per_1e6_b4=round(sum(b4p.values()) / max(n_b4, 1) * 1e6, 3))


def j5_addendum(new, old, b4):
    """J5b rule of docs/unification_criteria/B4_preplaced_addendum_2026-11-17.md §2: every pair class with at least one
    B4 module, rate = pairs per 1e6 instances of the class's B4 modules (sampled tiles); pass iff rate_new <= rate_ref +
    max(0.1 rate_ref, 1). A class absent from the reference counts as rate 0 there, unless one of its modules does not
    exist in the reference layers (new module, e.g. ATP synthase): then it is reported (pass None). Absolute pair
    counts (per class and all pre-pass pairs) are reported next to the rates."""
    rows = {}
    for c in sorted(set(new["class_rates_addendum"]) | set(old["class_rates_addendum"])):
        n, r = new["class_rates_addendum"].get(c), old["class_rates_addendum"].get(c)
        mods = c.split("|")
        newmod = any(m not in old["instances"] for m in mods)
        if r is None:                       # rate 0 in the reference: denominator from the reference instance counts
            ms = sorted(set(m for m in mods if m in b4))
            r = dict(pairs=0, den_modules=ms, den=int(sum(old["instances"].get(m, 0) for m in ms)), per_1e6=0.0)
        if n is None:
            ms = sorted(set(m for m in mods if m in b4))
            n = dict(pairs=0, den_modules=ms, den=int(sum(new["instances"].get(m, 0) for m in ms)), per_1e6=0.0)
        lim = r["per_1e6"] + max(0.1 * r["per_1e6"], 1.0)
        rows[c] = dict(pairs_new=n["pairs"], pairs_ref=r["pairs"], den_new=n["den"], den_ref=r["den"], den_modules=n["den_modules"],
                       per_1e6_new=n["per_1e6"], per_1e6_ref=r["per_1e6"], limit=round(lim, 3),
                       pass_=None if newmod else bool(n["per_1e6"] <= lim), new_module=newmod)
    return dict(rule="addendum 2026-11-17 §2: per pair class, pairs per 1e6 instances of the class's B4 modules <= reference "
                     "+ max(10 % rel., 1 per 1e6)", classes=rows,
                pass_b=all(x["pass_"] is not False for x in rows.values()),
                failing=[c for c, x in rows.items() if x["pass_"] is False],
                absolute=dict(prepass_prepass_pairs_new=new.get("prepass_prepass_pairs"),
                              prepass_prepass_pairs_ref=old.get("prepass_prepass_pairs"),
                              per_1e6_all_new=new["prepass_prepass_per_1e6"], per_1e6_all_ref=old["prepass_prepass_per_1e6"],
                              instances_new=new["instances"], instances_ref=old["instances"]))


def measure(args):
    import pandas as pd
    from wholecell import config as C
    cfg = C.load(args.cell or "jurkat")
    serial = args.serial
    t0 = time.time()
    res = dict(criteria="docs/unification_criteria/B4_preplaced.md", serial=serial)
    man = json.load(open(os.path.join(serial, "serial_manifest.json")))
    res["builders"] = [dict(name=b["name"], seconds=b["seconds"]) for b in man["builders"]]
    mods = {f[:-4]: os.path.join(serial, "modules", f) for f in os.listdir(os.path.join(serial, "modules")) if f.endswith(".npz")}
    js = {k: json.load(open(v[:-4] + ".json")) for k, v in mods.items()}
    # J1 / J1b
    from wholecell.pack import geomload
    g = geomload.load(cfg, args.build if args.build and os.path.exists(os.path.join(args.build, "geometry", "package.json")) else None)
    pores = np.asarray(g.cell.pores)
    from scipy.spatial import cKDTree
    d, _ = cKDTree(pores[:, :3]).query(pores[:, :3], k=2)
    n_cfg = next(b["params"].get("n_pores") for b in cfg["preplaced"]["serial"]["builders"] if b["builder"] == "core.npc_pores")
    I = pd.read_csv(os.path.join(args.build, "inventory", "pack_table.tsv"), sep="\t")
    npc_inv = int(I.loc[I.ing_id == "NPC", "copies"].iloc[0]) if (I.ing_id == "NPC").any() else None
    anchor = 2240.0
    res["J1"] = dict(config=n_cfg, package=int(len(pores)), placed=js.get("npc", {}).get("n"), inventory_npc=npc_inv,
                     min_sep_nm=round(float(d[:, 1].min()), 2), nup_anchor_pores=anchor, anchor_ratio=round(anchor / max(n_cfg, 1), 3))
    res["J1"]["pass"] = bool(n_cfg == len(pores) == res["J1"]["placed"] and 0.5 <= res["J1"]["anchor_ratio"] <= 2
                             and res["J1"]["min_sep_nm"] >= 150.0 - 1e-6)
    res["J1"]["inventory_ok"] = npc_inv == n_cfg
    # J2
    mt = js.get("microtubules", {})
    n_mt = int(mt.get("n", 0))
    res["J2"] = dict(dimers=n_mt, v7=V7["mt_dimers"], rel=round(n_mt / V7["mt_dimers"] - 1, 5), json=mt, pass_=abs(n_mt / V7["mt_dimers"] - 1) <= 0.01)
    # J3 / J4
    z = np.load(mods["polysomes"])
    pool = float(I.loc[I.ing_id == "RIBO80S", "copies"].iloc[0])
    cnt = np.bincount(z["chain"])
    in2 = int((cnt[z["chain"]] >= 2).sum())
    PM = z["pos"].astype(np.float64)
    if args.packed:
        PK = np.load(args.packed); PM_all = np.concatenate([PM, PK])
    else:
        PK = None; PM_all = PM
    rs_pre = ribosome_stats(PM)
    rs_all = ribosome_stats(PM_all) if PK is not None else None
    res["J3"] = dict(pool=pool, preplaced=int(len(PM)), in_chains_ge2=in2, a=round(in2 / pool, 4), target=0.70,
                     b_preplaced_only=rs_pre["linked_34nm_thinned_to_dedup_EM"],
                     b_with_packed=rs_all["linked_34nm_thinned_to_dedup_EM"] if rs_all else None, em=EM_LINKED,
                     priors_placed=int((z["origin"] == 0).sum()), er_ribosomes=int(z["er"].sum()), v7_a=round(V7["in_chains"], 4))
    b = res["J3"]["b_with_packed"] if rs_all else None
    res["J3"]["pass_a"] = bool(0.68 <= res["J3"]["a"] <= 0.72)
    res["J3"]["pass_b"] = None if b is None else bool(abs(b - EM_LINKED) <= 0.05)
    res["J4"] = dict(preplaced=rs_pre, with_packed=rs_all, v7=V7["nn_thinned"], window=[27.24, 29.24])
    if rs_all:
        res["J4"]["pass"] = bool(27.24 <= rs_all["nn_median_thinned_to_EM"] <= 29.24)
    # J6
    atp = js.get("atp_synthase", {})
    res["J6"] = dict(atp_monomers=atp.get("n"), atp_consumed=atp.get("consumed"),
                     mtdna=json.load(open(os.path.join(serial, "modules", "mtdna.json"))) if "mtdna" in mods else None,
                     pool_cpx6151=float(I.loc[I.ing_id == "CPX-6151", "copies"].iloc[0]) if (I.ing_id == "CPX-6151").any() else None)
    # J7 consumption
    res["J7"] = {k: v.get("consumed") for k, v in js.items()}
    # J5b
    if args.ref_build and not args.no_overlaps:
        LIB = _lib_trees(args.build)
        b4 = ("polysomes", "microtubules", "npc", "atp_synthase")
        log("self-overlaps: new layers")
        new = self_overlaps(_layers_new(serial, cfg, args.lamina), LIB, args.ref_build, cfg, b4)
        res["lamina_layer"] = args.lamina or "frozen v7 (pre_lamina_v7)"
        log("self-overlaps: v7 layers")
        old = self_overlaps(_layers_v7(cfg), LIB, args.ref_build, cfg, b4)
        cls = {}
        for c, v in new["b4_pair_classes_per_1e6_b4"].items():
            r = old["b4_pair_classes_per_1e6_b4"].get(c)
            cls[c] = dict(new=v, ref=r, pass_=(r is None) or v <= r + max(0.1 * r, 1.0))
        res["J5"] = dict(new=new, v7=old, classes=cls, pass_b=all(x["pass_"] for x in cls.values()),
                         total_new_le_v7=bool(new["prepass_prepass_per_1e6"] <= old["prepass_prepass_per_1e6"] * 1.1 + 1))
        res["J5"]["addendum"] = j5_addendum(new, old, b4)
    res["seconds"] = round(time.time() - t0)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    json.dump(res, open(args.out, "w"), indent=1, default=str)
    log(json.dumps({k: res[k] for k in ("J1", "J2") if k in res}, default=str)[:600])
    log("J3", json.dumps(res["J3"]), "J4", json.dumps(res["J4"]["preplaced"]))
    if "J5" in res:
        log("J5", json.dumps(res["J5"]["classes"]), "total new/v7", res["J5"]["new"]["prepass_prepass_per_1e6"], res["J5"]["v7"]["prepass_prepass_per_1e6"])
        log("J5b addendum", json.dumps({c: (x["pairs_new"], x["pairs_ref"], x["per_1e6_new"], x["per_1e6_ref"], x["pass_"])
                                        for c, x in res["J5"]["addendum"]["classes"].items()}), "failing", res["J5"]["addendum"]["failing"])
    return res


def _selected_npz(bdir, name):
    """npz the build selected for a pre-placed module with a builder (preplaced/<name>/selected.json), else None"""
    p = os.path.join(bdir, "preplaced", name, "selected.json")
    return json.load(open(p))["npz"] if os.path.exists(p) else None


def _packed_ribosomes(bdir):
    from wholecell.validate import common
    from wholecell.validate.ribosome_nn import packed_positions
    inp = common.Inputs(None, bdir) if False else None
    import glob
    L = json.load(open(os.path.join(bdir, "pack", "libs.json")))["libs"]
    if "RIBO80S" not in L:
        return np.zeros((0, 3))
    li = L.index("RIBO80S"); out = []
    for f in sorted(glob.glob(os.path.join(bdir, "pack", "tiles", "*.npz"))):
        z = np.load(f)
        k = z["types"] == li
        if k.any():
            out.append(z["pos"][k].astype(np.float64))
    return np.concatenate(out) if out else np.zeros((0, 3))


def _suite(cell, bdir, ref_suite, out):
    import subprocess
    nr = os.path.join(out, "suite", "no_regression.tsv")
    if not os.path.exists(nr):
        subprocess.run([sys.executable, "-m", "wholecell", "validate", "--cell", cell, "--build", bdir, "--step", "suite",
                        "--out", os.path.join(out, "suite"), f"ref_suite={ref_suite}"], cwd=REPO, check=False)
    return nr


EXPECTED = {   # B4_preplaced.md EXPECTED (Jurkat): metric substrings exempted from B2 (they are judged by J1-J7 instead)
    "jurkat": ("perinuclear", "density.cytosol", "density.matrix", "prepass_mrna", "prepass_prepass", "nn_median", "c12",
               "nup", "npc"),
    "b6": ("density.cytosol", "perinuclear", "nuc_cyto_ratio", "phi.cytosol", "er_lumen", "matrix", "golgi_lumen", "endo",
           "depends_on_geometry"),
}


def verdict(args):
    import pandas as pd
    bdir = args.build
    out = os.path.join(bdir, "validate", "compare_B4")
    os.makedirs(out, exist_ok=True)
    items, extra, add = {}, {}, {}
    if args.cell == "jurkat":
        PK = _packed_ribosomes(bdir)
        np.save(os.path.join(out, "packed_ribo80s.npy"), PK)
        a = argparse.Namespace(cell="jurkat", serial=os.path.join(bdir, "preplaced", "serial"), build=bdir,
                               ref_build=args.ref_build, out=os.path.join(out, "measure.json"),
                               packed=os.path.join(out, "packed_ribo80s.npy"), no_overlaps=False,
                               lamina=_selected_npz(bdir, "lamina"))
        m = measure(a)
        items.update(J1=m["J1"]["pass"] and m["J1"]["inventory_ok"], J2=m["J2"]["pass_"], J3=m["J3"]["pass_a"] and bool(m["J3"]["pass_b"]),
                     J4=m["J4"].get("pass", False), J6=bool(m["J6"]["atp_monomers"]) and m["J6"]["mtdna"] is not None)
        pl = os.path.join(bdir, "validate", "placement.json")
        if os.path.exists(pl):
            P = json.load(open(pl)); cp = P["clash_pairs"]["prepass_prepass"]; npre = max(P.get("n_prepass", 1), 1)
            rate = cp / npre * 1e6
            ref = args.ref_prepass or V7["prepass_prepass_per_1e6"]
            extra["J5a"] = dict(rate=rate, ref=ref, pairs=cp, n_prepass=npre)
            items["J5"] = bool(rate <= ref + max(0.1 * ref, 1.0)) and m.get("J5", {}).get("pass_b", False)
            # addendum 2026-11-17 §2: J5a (all pre-pass pairs, frozen rule) and J5b per pair class (addendum rule)
            add["J5"] = bool(rate <= ref + max(0.1 * ref, 1.0)) and m.get("J5", {}).get("addendum", {}).get("pass_b", False)
        else:
            items["J5"] = False; add["J5"] = False
        extra["J5b_addendum"] = m.get("J5", {}).get("addendum")
        from wholecell import config as C
        extra["polysome_linker_nt_mean"] = next((b["params"].get("linker_nt_mean") for b in C.load("jurkat")["preplaced"]["serial"]["builders"]
                                                 if b["builder"] == "core.polysomes"), None)
        extra["measure"] = os.path.join(out, "measure.json")
    else:                                   # yeast: Y1 arrays vs the PATCH-4 reference
        ref = np.load(os.path.join(YREF_DIR, "prepass.npz"))
        new = np.load(os.path.join(bdir, "preplaced", "serial", "preplaced.npz"))
        diff = [k for k in ref.files if k not in new.files or not (ref[k].shape == new[k].shape and np.array_equal(ref[k], new[k]))]
        extra["Y1"] = dict(ref=YREF_DIR, arrays=len(ref.files), differing=diff)
        items["Y1"] = not [k for k in diff if k not in ("mrna_gene", "poly_chain_gene")]
        items["Y2"] = int(len(new["npc_pos"])) == YREF["pores"]
        items["Y3"] = bool(np.array_equal(ref["mt_pos"], new["mt_pos"]))
    nr = _suite(args.cell, bdir, args.ref_suite, out) if args.ref_suite else None
    regress, exempt = [], []
    if nr and os.path.exists(nr):
        N = pd.read_csv(nr, sep="\t", keep_default_na=False)
        regress = N[N.status == "regress"].to_dict("records")
        keys = EXPECTED["jurkat"] if args.cell == "jurkat" else ()
        for e in (args.extra_expected or "").split(","):
            keys = keys + EXPECTED.get(e, ())
        exempt = [r for r in regress if any(k in (str(r.get("check", "")) + "." + str(r.get("metric", ""))).lower() for k in keys)]
        regress = [r for r in regress if r not in exempt]
        if args.cell == "yeast":
            S = pd.read_csv(nr, sep="\t", keep_default_na=False)
            def val(c, mt):
                x = S[(S.check == c) & (S.metric == mt)]
                return float(x.value.iloc[0]) if len(x) and x.value.iloc[0] not in ("", "None") else None
            nn = val("C12", "nn_median_model"); pp = val("C7", "a.prepass_prepass_per_1e6"); ic = val("C12", "polysome_in_chains_fraction")
            extra.update(Y4=dict(in_chains=ic, target=[0.82, 0.85], ref=YREF["in_chains"]), Y5=dict(nn=nn, ref=YREF["nn"]),
                         Y6=dict(rate=pp, ref=YREF["prepass_prepass_per_1e6"]))
            items["Y5"] = nn is not None and abs(nn - YREF["nn"]) <= 0.5
            items["Y6"] = pp is not None and pp <= YREF["prepass_prepass_per_1e6"] * 1.1 + 1
    items["B1"] = bool(args.complete) and os.path.exists(os.path.join(bdir, "preplaced", "serial", "serial_manifest.json"))
    items["B2"] = bool(nr and os.path.exists(nr)) and not regress
    fail = [k for k, v in items.items() if not v]
    # both rule sets are always reported: items / verdict_frozen = B4_preplaced.md as frozen; items_addendum / verdict =
    # with the addendum 2026-11-17 (J5b per pair class; J3b / J4 bounds unchanged, the polysome linker is a config value)
    items_add = dict(items, **add)
    fail_add = [k for k, v in items_add.items() if not v]
    Vd = dict(step="B4", cell=args.cell, criteria="docs/unification_criteria/B4_preplaced.md",
              addendum="docs/unification_criteria/B4_preplaced_addendum_2026-11-17.md", build=bdir, ref_suite=args.ref_suite,
              no_regression=nr, items=items, items_addendum=items_add, regressions=regress[:50], exempt_expected=exempt, **extra,
              verdict_frozen="ACCEPTED" if not fail else "NOT ACCEPTED", failing_frozen=fail,
              verdict="ACCEPTED" if not fail_add else "NOT ACCEPTED", failing=fail_add)
    json.dump(Vd, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    log("verdict (with addendum)", Vd["verdict"], fail_add, "| frozen", Vd["verdict_frozen"], fail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["measure", "verdict"])
    ap.add_argument("--cell", default="jurkat"); ap.add_argument("--serial"); ap.add_argument("--build"); ap.add_argument("--ref-build")
    ap.add_argument("--out"); ap.add_argument("--packed"); ap.add_argument("--no-overlaps", action="store_true")
    ap.add_argument("--ref-suite"); ap.add_argument("--complete", type=int, default=1); ap.add_argument("--extra-expected", default="")
    ap.add_argument("--ref-prepass", type=float)
    ap.add_argument("--lamina", help="measure: lamina npz of the new layers (default: the frozen v7 lamina)")
    a = ap.parse_args()
    {"measure": measure, "verdict": verdict}[a.what](a)


if __name__ == "__main__":
    main()
