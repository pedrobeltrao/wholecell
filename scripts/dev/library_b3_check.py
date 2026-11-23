"""B3 (Jurkat library rebuilt through the core library stage): organism tables and the acceptance check.
Criteria: docs/unification_criteria/B3_library.md (frozen 2026-11-17, md5 b6c09336...) and its addendum
docs/unification_criteria/B3_library_addendum_2026-11-17.md (L5 belt measure redefined after the library-only result;
both L5 values are always reported: "L5 (frozen definition)" = key L5, "L5 (addendum 2026-11-17)" = key L5a, which is
the one the verdict uses).

  species   write plugins/jurkat/data/library/{species_v7.tsv.gz, start_v7.tsv, tm_route_ids.tsv} from the legacy
            tables the v7 library was built from (ingredients_v0_pre_v4.tsv, inventory_v0_complexed.tsv,
            data/derived/v5/tm_route_ids.tsv) and the `src` recorded in each v7 npz (run once; the outputs are committed)
  measure   per-species measures of a library folder -> TSV (volume of the sphere-tree union, Kirkwood R_h, beads,
            spheres, membrane, belt centre, ...)
  compare   L1-L8 of the criteria: a built library (<build>/library) vs the frozen v7 library -> JSON + per-species TSV
  verdict   compare_B3/verdict.json of a build: L1-L7 on the build's library + Y1/Y2 + B1-B4 (suite no_regression)
  rescore   L5 (frozen) and L5a (addendum) of an existing compare folder (per_species.tsv) -> <out>/l5_rescore.json

  python scripts/library_b3_check.py compare --build <bdir> --out <dir> [--nproc 3] [--other <bdir2>]
  python scripts/library_b3_check.py verdict --build <bdir> --ref-suite <suite.json> [--complete 1] [--b6-expected]
                                            [--cache <dir with measures_v7.tsv>] [--other <library run for L8>]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import time

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
H = 0.25                  # nm, voxel spacing of the tree volume (criteria L2)
RH_MAX = 3000             # beads; Kirkwood subset above this (criteria L3)
SLAB = 2.0                # nm, slab half thickness of the belt fallback (criteria L5)
BELT_MAX = 0.5            # nm, |z_belt| threshold (L5 and L5a)
SEG_CACHE = "l5a_tm_segments.pkl.gz"   # TRANSMEM bead indices per species (addendum L5a), cached next to measures_v7
TWO_MEMBRANES = ("GJA1_x12", "GJB2_x12", "VANGL1_x6", "VANGL2_x6")   # jlib_tm_route TWO_MEMB: slab from one layer
LIBDIR = "plugins/jurkat/data/library"


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def _md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------ species tables (run once, committed)
def make_species(args):
    import pandas as pd
    from wholecell.config import load
    cfg = load("jurkat")
    dr = cfg["_site"]["data_root"]
    lib = cfg["frozen_inputs"]["lib_structures_v7"]["path"]
    I = pd.read_csv(f"{dr}/data/derived/ingredients_v0_pre_v4.tsv", sep="\t")
    acc = pd.read_csv(f"{dr}/data/derived/inventory_v0_complexed.tsv", sep="\t", usecols=["gene", "uniprot"]).set_index("gene").uniprot
    src7, n_diff = {}, 0
    for i in I.ing_id:
        f = f"{lib}/{i}.npz"
        if os.path.exists(f):
            z = np.load(f)
            if "src" in z.files:
                src7[i] = str(z["src"])
    out = []
    for r in I.itertuples():
        genes = r.genes if isinstance(r.genes, str) else ""
        accs = ";".join(acc.get(g, "") if isinstance(acc.get(g, ""), str) else "" for g in genes.split(";")) if genes else ""
        s = r.struct_src
        lib_src = src7.get(r.ing_id, "")
        if lib_src and lib_src != s and not lib_src.startswith("afdb_monomer:"):
            n_diff += 1
            s_used, note = lib_src, f"v7 library src (inventory struct_src {r.struct_src})"
        elif lib_src.startswith("afdb_monomer:"):
            s_used, note = s, ""                              # tm_route fragment entries record the monomer
        else:
            s_used, note = s, ""
        out.append(dict(ing_id=r.ing_id, kind=r.kind, genes=genes, accs=accs, stoich=r.stoich if isinstance(r.stoich, str) else "",
                        struct_src=s_used, pdb=s_used.split(":", 1)[1] if s_used.startswith("pdb:") else "",
                        uniprot=r.uniprot if isinstance(r.uniprot, str) else "", transmembrane=bool(r.transmembrane),
                        signal_frac=r.signal_frac, mass_da=r.mass_da, src_note=note))
    S = pd.DataFrame(out)
    d = os.path.join(REPO, LIBDIR)
    os.makedirs(d, exist_ok=True)
    with gzip.GzipFile(os.path.join(d, "species_v7.tsv.gz"), "wb", mtime=0) as fh:
        fh.write(S.to_csv(sep="\t", index=False).encode())
    T = pd.read_csv(f"{dr}/data/derived/v5/tm_route_ids.tsv", sep="\t")
    T.to_csv(os.path.join(d, "tm_route_ids.tsv"), sep="\t", index=False)
    # start entries: v7 files no scripted step writes (not an ingredient of the species table)
    man = pd.read_csv(cfg["frozen_inputs"]["lib_structures_v7"]["manifest"], sep="\t")
    ids = set(S.ing_id)
    st = man[~man.file.str[:-4].isin(ids) | ~man.file.str.endswith(".npz")].copy()
    st["reason"] = st.file.map(lambda f: "pre-pass template (jcytoskeleton: helix 9GO5 / lattice 8QV3, 5W3F), not built by the library scripts")
    st[["file", "md5", "reason"]].to_csv(os.path.join(d, "start_v7.tsv"), sep="\t", index=False)
    log(f"species {len(S)} rows ({n_diff} with the v7 library src != inventory struct_src); start entries {len(st)}; "
        f"tm_route ids {len(T)}")


# ------------------------------------------------------------------ measures
def tree_volume(c, r, h=H):
    """union volume of spheres on the absolute grid (i, j, k) * h"""
    c = np.asarray(c, np.float64); r = np.asarray(r, np.float64)
    if not len(r):
        return 0.0
    lo = np.floor((c - r[:, None]).min(0) / h).astype(np.int64)
    hi = np.ceil((c + r[:, None]).max(0) / h).astype(np.int64) + 1
    occ = np.zeros(tuple(hi - lo), bool)
    for cc, rr in zip(c, r):
        a = np.floor((cc - rr) / h).astype(np.int64); b = np.ceil((cc + rr) / h).astype(np.int64) + 1
        g = [np.arange(a[k], b[k]) * h - cc[k] for k in range(3)]
        m = (g[0][:, None, None] ** 2 + g[1][None, :, None] ** 2 + g[2][None, None, :] ** 2) < rr * rr
        sl = tuple(slice(a[k] - lo[k], b[k] - lo[k]) for k in range(3))
        occ[sl] |= m
    return float(occ.sum()) * h ** 3


def kirkwood_rh(x):
    x = np.asarray(x, np.float64)
    n = len(x)
    if n < 2:
        return float("nan")
    if n > RH_MAX:
        x = x[np.sort(np.random.default_rng(0).choice(n, RH_MAX, replace=False))]
        n = RH_MAX
    s = 0.0
    for i in range(0, n, 500):
        d = np.sqrt(((x[i:i + 500, None, :] - x[None, :, :]) ** 2).sum(-1))
        d[d == 0] = np.inf
        s += (1.0 / d).sum()
    return float(n * (n - 1) / s)


def measure_one(args):
    sid, path = args
    try:
        z = np.load(path)
        b = z["beads"].astype(np.float64)
        rec = dict(id=sid, md5=_md5(path), nb=len(b), ns=len(z["st_r"]), vol=tree_volume(z["st_c"], z["st_r"]),
                   rh=kirkwood_rh(b), membrane=str(z["membrane"]) if "membrane" in z.files else "",
                   tm_route=str(z["tm_route"]) if "tm_route" in z.files else "",
                   beads_md5=hashlib.md5(np.ascontiguousarray(z["beads"]).tobytes()).hexdigest(),
                   tree_md5=hashlib.md5(np.ascontiguousarray(z["st_c"]).tobytes() + np.ascontiguousarray(z["st_r"]).tobytes()).hexdigest())
        if rec["membrane"]:
            if "tm_beads" in z.files and z["tm_beads"].astype(bool).any():
                zz = b[z["tm_beads"].astype(bool), 2]; how = "tm_beads"
            else:
                zz = b[np.abs(b[:, 2]) < SLAB, 2]; how = "slab"
            rec.update(z_belt=float(np.median(zz)) if len(zz) else float("nan"), belt_from=how)
        return rec
    except Exception as e:  # noqa: BLE001 - reported
        return dict(id=sid, error=repr(e)[:200])


def measure(entries, nproc, out_tsv):
    import pandas as pd
    from concurrent.futures import ProcessPoolExecutor
    if os.path.exists(out_tsv):
        return pd.read_csv(out_tsv, sep="\t", keep_default_na=False, na_values=[""])
    jobs = sorted(entries.items(), key=lambda kv: -os.path.getsize(kv[1]))          # big first
    with ProcessPoolExecutor(nproc) as ex:
        R = list(ex.map(measure_one, jobs, chunksize=4))
    T = pd.DataFrame(R).sort_values("id")
    T.to_csv(out_tsv, sep="\t", index=False)
    return pd.read_csv(out_tsv, sep="\t", keep_default_na=False, na_values=[""])


def _entries(folder):
    return {f[:-4]: os.path.join(folder, f) for f in os.listdir(folder) if f.endswith(".npz")}


def kabsch_tilt_shift(a, b):
    """superpose a (v7 beads) onto b (new beads): tilt of the v7 normal in the new frame (deg), z of the v7 origin"""
    pa, pb = a.mean(0), b.mean(0)
    Hm = (a - pa).T @ (b - pb)
    U, S, Vt = np.linalg.svd(Hm)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = pb - R @ pa
    n = R @ np.array([0, 0, 1.0])
    return float(np.degrees(np.arccos(np.clip(n[2], -1, 1)))), float(t[2])


# ------------------------------------------------------------------ L5a: TRANSMEM-segment belt (addendum 2026-11-17)
_SEG = {}


def _seg_one(row):
    """bead indices of the UniProt TRANSMEM residues (+ beta-barrel strands, as membrane_orient.masks) of one species, in
    the bead order of the library builder, re-derived from the species' structure source exactly as the builder read it;
    plus a sample of those beads' source coordinates (to verify the mapping on the library beads)"""
    L = _SEG["L"]
    jb, tg, mo = L["jb"], L["tg"], L["mo"]
    sid, src, accs, uni, tm = row
    try:
        kind_, ident = src.split(":", 1)
        if kind_ == "afdb_monomer":
            acc = ident
            m = tg.af_model(acc)
            if m is None:
                return sid, None, "no AFDB F1 model"
            pos, xyz, pl = m
            rng = mo.topology().get(acc, ([], [], []))[0]
            tmm = np.zeros(len(pos), bool)
            for a, b in rng:
                tmm |= (pos >= a) & (pos <= b)
            keep = pl >= 50                                   # build_structs.process_af / jlib_reorient_tm selection
            if rng:
                keep |= tmm
            if keep.sum() < 10:
                keep[:] = True
            xyz, tmm = xyz[keep], tmm[keep]
        elif kind_ == "pdb":
            xyz, kind, acc, pos, chn = jb.beads_pdb(ident, return_chain=True)
            if tm and sid != "NPC":
                xyz, kind, acc, pos, _ = tg.graft_missing_tm(xyz, kind, acc, pos, chn)
            tmm = mo.masks(acc, pos, xyz)[0]
        else:
            xyz, kind, acc, pos = jb.beads_afdb(ident, accs)
            tmm = mo.masks(acc, pos, xyz)[0]
        idx = np.flatnonzero(tmm).astype(np.int32)
        sub = idx[np.unique(np.linspace(0, len(idx) - 1, min(len(idx), 200)).astype(int))] if len(idx) else idx
        return sid, dict(nb=int(len(xyz)), idx=idx, sub=sub, ref=np.asarray(xyz[sub], np.float32)), ""
    except Exception as e:  # noqa: BLE001 - reported
        return sid, None, repr(e)[:160]


def tm_segment_masks(ids, nproc, cache):
    """{species: {nb, idx, sub, ref} | None} for the given species (TM species without `tm_beads`); cached (pickle)"""
    import pickle
    import tempfile
    import pandas as pd
    have = {}
    if os.path.exists(cache):
        with gzip.open(cache, "rb") as fh:
            have = pickle.load(fh)
    todo = sorted(set(ids) - set(have))
    if todo:
        from wholecell.config import load
        from wholecell.library.build import context
        from wholecell.library import jurkat as J
        cfg = load("jurkat")
        env0 = dict(os.environ)                                       # legacy() sets builder env vars: restored below
        tmp = tempfile.mkdtemp(prefix="l5a_", dir=os.path.dirname(os.path.abspath(cache)))
        X = context(cfg, tmp, nproc)
        X["partial"] = os.path.join(tmp, "lib.partial")               # legacy modules need a writable LIB; nothing is built
        os.makedirs(X.partial, exist_ok=True)
        _SEG["L"] = L = J.legacy(X)
        L["jb"].pdb_file = lambda pid: f"{X.pdb_dir}/{pid}-assembly1.cif.gz"   # read-only (no `used_pdb` bookkeeping)
        L["jb"].sifts(); L["mo"].topology(); L["mo"].barrels()               # loaded once, before the fork
        L["tg"]._IDX = json.load(open(L["tg"].IDX))
        sp = pd.read_csv(os.path.join(REPO, LIBDIR, "species_v7.tsv.gz"), sep="\t", keep_default_na=False).set_index("ing_id")
        rows = []
        for i in todo:
            r = sp.loc[i]
            accs = r.accs.split(";") if r.accs else [""] * len(r.genes.split(";"))
            rows.append((i, r.struct_src, accs, r.uniprot, str(r.transmembrane).lower() in ("true", "1")))
        log(f"L5a: TRANSMEM bead masks of {len(rows)} species ({nproc} processes)")
        with J._pool(nproc) as ex:
            for sid, rec, note in ex.map(_seg_one, rows, chunksize=4):
                have[sid] = rec if rec is not None else dict(error=note)
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)
        os.environ.clear(); os.environ.update(env0)
        with gzip.open(cache + ".tmp", "wb") as fh:
            pickle.dump(have, fh)
        os.replace(cache + ".tmp", cache)
    return {i: have.get(i) for i in ids}


def belt_addendum(lib_dir, ids, segs):
    """L5a z_belt per TM species of one library folder -> DataFrame (z_belt_a, belt_a_from, n_belt)
       tm_beads key (TM-routed species): beads flagged 1 (TRANSMEM + INTRAMEM, as stored);
       else: beads of the UniProt TRANSMEM residues (segs), when the mapping verifies on this library (same bead count,
       pairwise distances of the sampled TM beads within 0.02 nm of the source structure's);
       monotopic (`intramem_only`): tm_beads species whose accessions have no TRANSMEM feature -> not applicable;
       otherwise the frozen definition (`fallback_slab`)."""
    import pandas as pd
    from scipy.spatial.distance import pdist
    up = None
    rows = []
    for i in ids:
        z = np.load(os.path.join(lib_dir, i + ".npz"))
        b = z["beads"].astype(np.float64)
        rec = dict(id=i)
        if "tm_beads" in z.files and z["tm_beads"].astype(bool).any():
            m = z["tm_beads"].astype(bool)
            if up is None:
                up = _tm_route_tm_accs()
            how = "tm_beads" if up.get(i, True) else "intramem_only"
            zz = b[m, 2]
        else:
            s = segs.get(i) or {}
            ok = bool(s.get("nb") == len(b) and len(s.get("idx", [])) >= 1)
            if ok and len(s["sub"]) >= 2:
                ok = float(np.abs(pdist(b[s["sub"]]) - pdist(s["ref"].astype(np.float64))).max()) <= 0.02
            if ok:
                zz, how = b[s["idx"], 2], "tm_segments"
            else:
                zz, how = b[np.abs(b[:, 2]) < SLAB, 2], "fallback_slab"
                rec["map_note"] = s.get("error", "") or (f"bead count {s.get('nb')} vs {len(b)}" if s.get("nb") != len(b)
                                                         else "TRANSMEM bead distances differ / no TRANSMEM bead")
        rec.update(z_belt_a=float(np.median(zz)) if len(zz) else float("nan"), belt_a_from=how, n_belt=int(len(zz)))
        rows.append(rec)
    return pd.DataFrame(rows).set_index("id")


def _tm_route_tm_accs():
    """{ing_id: True if any accession of the TM-routed species has a UniProt TRANSMEM feature} (REST cache of tm_route)"""
    import pandas as pd
    from wholecell.config import load
    cfg = load("jurkat")
    path = cfg["frozen_inputs"]["lib_tm_route_uniprot"]["path"]
    U = json.load(open(path))
    sp = pd.read_csv(os.path.join(REPO, LIBDIR, "species_v7.tsv.gz"), sep="\t", keep_default_na=False).set_index("ing_id")
    out = {}
    for i, r in sp.iterrows():
        accs = [a for a in (r.accs.split(";") if r.accs else []) if a] or ([r.uniprot] if r.uniprot else [])
        out[i] = any(f[0] == "Transmembrane" for a in accs for f in U.get(a, {}).get("features", []))
    return out


def score_l5(C, n_tm_ids):
    """L5 (frozen definition) and L5a (addendum 2026-11-17) from a per-species table with z_belt_* / z_belt_a_* columns"""
    tm = C.loc[n_tm_ids]
    tm7 = C[C.membrane_v7.fillna("") != ""]                          # frozen L5 code: v7 count over the v7 TM set
    f_new = tm.z_belt_new.abs() > BELT_MAX
    f_v7 = tm7.z_belt_v7.abs() > BELT_MAX
    L5 = dict(label="L5 (frozen definition)", n_tm=int(len(tm)), fail_new=int(f_new.sum()), fail_v7=int(f_v7.sum()),
              frac_ok=round(1 - f_new.sum() / max(len(tm), 1), 5), frac_ok_v7=round(1 - f_v7.sum() / max(len(tm7), 1), 5),
              median_abs_new=float(tm.z_belt_new.abs().median()), fail_ids=tm.index[f_new][:30].tolist())
    L5["pass_"] = bool(L5["frac_ok"] >= 0.99 and L5["fail_new"] <= L5["fail_v7"] + 0.01 * len(tm))
    na = (tm.belt_a_from_new == "intramem_only")
    ap = tm[~na]
    g_new = ap.z_belt_a_new.abs() > BELT_MAX
    g_v7 = ap.z_belt_a_v7.abs() > BELT_MAX
    L5a = dict(label="L5 (addendum 2026-11-17)", definition="docs/unification_criteria/B3_library_addendum_2026-11-17.md",
               n_tm=int(len(tm)), not_applicable_intramem_only=int(na.sum()), n_scored=int(len(ap)),
               fail_new=int(g_new.sum()), fail_v7=int(g_v7.sum()),
               frac_ok=round(1 - g_new.sum() / max(len(ap), 1), 5), frac_ok_v7=round(1 - g_v7.sum() / max(len(ap), 1), 5),
               frac_ok_all_tm_counting_na_as_fail=round(1 - (g_new.sum() + na.sum()) / max(len(tm), 1), 5),
               belt_from_new=tm.belt_a_from_new.value_counts().to_dict(), belt_from_v7=tm.belt_a_from_v7.value_counts().to_dict(),
               median_abs_new=float(ap.z_belt_a_new.abs().median()), p99_abs_new=float(ap.z_belt_a_new.abs().quantile(0.99)),
               fail_ids=ap.index[g_new].tolist()[:60],
               fail_two_membranes=[i for i in ap.index[g_new] if i in TWO_MEMBRANES],
               intramem_only_ids=tm.index[na].tolist(),
               fallback_ids=tm.index[tm.belt_a_from_new == "fallback_slab"].tolist()[:60])
    L5a["pass_"] = bool(L5a["frac_ok"] >= 0.99 and L5a["fail_new"] <= L5a["fail_v7"] + 0.01 * len(ap))
    return L5, L5a


def add_belt_addendum(C, v7dir, lib, nproc, cache_dir):
    """adds z_belt_a_{v7,new}, belt_a_from_{v7,new}, n_belt_{v7,new} for the TM species of the new library"""
    mem = C.membrane_new.fillna("")
    ids = C.index[mem != ""].tolist()
    need = []
    for i in ids:
        if C.at[i, "belt_from_new"] != "tm_beads" or C.at[i, "belt_from_v7"] != "tm_beads":
            need.append(i)
    segs = tm_segment_masks(need, nproc, os.path.join(cache_dir, SEG_CACHE))
    for lab, d in (("v7", v7dir), ("new", lib)):
        B = belt_addendum(d, ids, segs)
        for k in B.columns:
            C.loc[ids, f"{k}_{lab}"] = B[k]
    return ids


def compare(args):
    import pandas as pd
    from wholecell.config import load
    cfg = load("jurkat")
    fi = dict(cfg["frozen_inputs"]["lib_structures_v7"])
    if args.v7dir:                                                    # tests: another reference folder + manifest
        fi.update(path=args.v7dir, manifest=args.v7man)
    v7dir = fi["path"]
    os.makedirs(args.out, exist_ok=True)
    lib = os.path.join(args.build, "library", "lib")
    t0 = time.time()
    man7 = pd.read_csv(fi["manifest"], sep="\t")
    have = sorted(os.listdir(lib))
    L1 = dict(v7_files=len(man7), new_files=len(have), missing=sorted(set(man7.file) - set(have)),
              extra=sorted(set(have) - set(man7.file)))
    L1["pass"] = not L1["missing"] and not L1["extra"]
    m7 = os.path.join(args.cache or args.out, "measures_v7.tsv")
    A = measure(_entries(v7dir), args.nproc, m7).set_index("id")
    log(f"v7 measured ({time.time() - t0:.0f} s)")
    B = measure(_entries(lib), args.nproc, os.path.join(args.out, "measures_new.tsv")).set_index("id")
    log(f"new measured ({time.time() - t0:.0f} s)")
    ids = sorted(set(A.index) & set(B.index))
    C = pd.DataFrame(index=ids)
    for k in ("nb", "ns", "vol", "rh", "membrane", "tm_route", "beads_md5", "tree_md5", "md5", "z_belt", "belt_from"):
        C[k + "_v7"] = A.loc[ids, k] if k in A else np.nan
        C[k + "_new"] = B.loc[ids, k] if k in B else np.nan
    C["dvol"] = (C.vol_new - C.vol_v7) / C.vol_v7
    C["drh"] = (C.rh_new - C.rh_v7) / C.rh_v7
    ok_v = C.dvol.abs() <= 0.02
    ok_r = C.drh.abs() <= 0.02
    n = len(C)
    L2 = dict(n=n, within_2pct=int(ok_v.sum()), frac=round(float(ok_v.mean()), 5), pass_=bool(ok_v.mean() >= 0.99),
              median_abs=float(C.dvol.abs().median()), p99_abs=float(C.dvol.abs().quantile(0.99)),
              worst=C.dvol.abs().sort_values(ascending=False).head(15).round(4).to_dict())
    L3 = dict(n=n, within_2pct=int(ok_r.sum()), frac=round(float(ok_r.mean()), 5), pass_=bool(ok_r.mean() >= 0.99),
              median_abs=float(C.drh.abs().median()), p99_abs=float(C.drh.abs().quantile(0.99)),
              joint_frac=round(float((ok_v & ok_r).mean()), 5),
              beads_identical=int((C.beads_md5_v7 == C.beads_md5_new).sum()), nb_identical=int((C.nb_v7 == C.nb_new).sum()),
              tree_identical=int((C.tree_md5_v7 == C.tree_md5_new).sum()), file_identical=int((C.md5_v7 == C.md5_new).sum()),
              worst=C.drh.abs().sort_values(ascending=False).head(15).round(4).to_dict())
    mem7, memn = C.membrane_v7.fillna(""), C.membrane_new.fillna("")
    tm = C[(mem7 != "") | (memn != "")]
    diff = tm[mem7[tm.index] != memn[tm.index]]
    r7 = C.tm_route_v7.fillna(""); rn = C.tm_route_new.fillna("")
    routed = C[r7 != ""]
    L4 = dict(tm_v7=int((mem7 != "").sum()), tm_new=int((memn != "").sum()), differ=int(len(diff)),
              differ_ids=diff.index[:40].tolist(), routed_v7=int(len(routed)),
              routed_same_method=int((r7[routed.index] == rn[routed.index]).sum()),
              routed_differ=routed.index[r7[routed.index] != rn[routed.index]].tolist())
    L4["pass_"] = bool(L4["differ"] == 0 and L4["routed_same_method"] >= 0.95 * max(L4["routed_v7"], 1))
    tmn = C[memn != ""]
    add_belt_addendum(C, v7dir, lib, args.nproc, args.cache or args.out)
    L5, L5a = score_l5(C, tmn.index)
    log(f"L5 (frozen definition): fail new {L5['fail_new']} / v7 {L5['fail_v7']} of {L5['n_tm']}; "
        f"L5 (addendum 2026-11-17): fail new {L5a['fail_new']} / v7 {L5a['fail_v7']} of {L5a['n_scored']}")
    rows = []
    for i in tmn.index[(C.loc[tmn.index, "nb_v7"] == C.loc[tmn.index, "nb_new"]).values]:
        a = np.load(os.path.join(v7dir, i + ".npz"))["beads"].astype(np.float64)
        b = np.load(os.path.join(lib, i + ".npz"))["beads"].astype(np.float64)
        tilt, sh = kabsch_tilt_shift(a, b)
        rows.append((i, tilt, sh))
    K = pd.DataFrame(rows, columns=["id", "tilt", "zshift"]).set_index("id")
    okk = (K.tilt <= 10) & (K.zshift.abs() <= 0.5)
    L6 = dict(n=int(len(K)), ok=int(okk.sum()), frac=round(float(okk.mean()), 5) if len(K) else None,
              flipped=K.index[K.tilt > 90].tolist(), median_tilt=float(K.tilt.median()) if len(K) else None,
              p99_tilt=float(K.tilt.quantile(0.99)) if len(K) else None, max_abs_shift=float(K.zshift.abs().max()) if len(K) else None,
              worst=K.tilt.sort_values(ascending=False).head(10).round(2).to_dict())
    L6["pass_"] = bool(len(K) and L6["frac"] >= 0.99 and not L6["flipped"])
    C = C.join(K)
    # L7: sphere counts and the adaptive rule
    sp = pd.read_csv(os.path.join(REPO, LIBDIR, "species_v7.tsv.gz"), sep="\t", keep_default_na=False).set_index("ing_id")
    kda = pd.to_numeric(sp.mass_da, errors="coerce").reindex(C.index) / 1e3
    vp = np.maximum(kda * 1.21, 1.0)
    big = C.nb_v7 > 20000
    rated = C[kda.notna() & ~big]
    over7 = float(((rated.vol_v7 / vp[rated.index]) > 2.5).mean()); overn = float(((rated.vol_new / vp[rated.index]) > 2.5).mean())
    L7 = dict(ns_identical=int((C.ns_v7 == C.ns_new).sum()), frac_identical=round(float((C.ns_v7 == C.ns_new).mean()), 5),
              spheres_v7=int(C.ns_v7.sum()), spheres_new=int(C.ns_new.sum()),
              d_total=round(float(C.ns_new.sum() / C.ns_v7.sum() - 1), 5), over_2p5_v7=round(over7, 5), over_2p5_new=round(overn, 5),
              note="tree volume here on the 0.25 nm grid (the builder's rule uses 0.5 nm), mass from the species table")
    L7["pass_"] = bool(L7["frac_identical"] >= 0.95 and abs(L7["d_total"]) <= 0.02 and overn <= over7 + 0.005)
    res = dict(criteria="docs/unification_criteria/B3_library.md", build=args.build, library=lib, v7=v7dir,
               addendum="docs/unification_criteria/B3_library_addendum_2026-11-17.md (L5a replaces L5 in the pass rule; "
                        "L5 = frozen definition, reported)",
               L1=L1, L2=L2, L3=L3, L4=L4, L5=L5, L5a=L5a, L6=L6, L7=L7, seconds=round(time.time() - t0))
    if args.other:                                                    # L8 determinism vs another build of the commit
        ma = {r.split("\t")[0]: r.split("\t")[2] for r in list(open(os.path.join(args.build, "library", "library_manifest.tsv")))[1:]}
        mb = {r.split("\t")[0]: r.split("\t")[2] for r in list(open(os.path.join(args.other, "library", "library_manifest.tsv")))[1:]}
        dif = sorted(k for k in set(ma) | set(mb) if ma.get(k) != mb.get(k))
        res["L8"] = dict(other=args.other, n=len(ma), differ=len(dif), differ_ids=dif[:30], pass_=not dif)
    C.to_csv(os.path.join(args.out, "per_species.tsv"), sep="\t")
    for k in ("L2", "L3", "L4", "L5", "L5a", "L6", "L7", "L8"):
        if k in res and "pass_" in res[k]:
            res[k]["pass"] = res[k].pop("pass_")
    res["pass_frozen_L5"] = all(res[k]["pass"] for k in ("L1", "L2", "L3", "L4", "L5", "L6", "L7"))
    res["pass"] = all(res[k]["pass"] for k in ("L1", "L2", "L3", "L4", "L5a", "L6", "L7"))   # addendum 2026-11-17
    json.dump(res, open(os.path.join(args.out, "library_compare.json"), "w"), indent=1, default=str)
    log(json.dumps({k: (v.get("pass") if isinstance(v, dict) else v) for k, v in res.items() if k.startswith("L")}))
    return res


def rescore(args):
    """L5 + L5a of an existing compare folder (its per_species.tsv: frozen-definition z_belt columns) -> <out>"""
    import pandas as pd
    from wholecell.config import load
    v7dir = args.v7dir or load("jurkat")["frozen_inputs"]["lib_structures_v7"]["path"]
    lib = os.path.join(args.build, "library", "lib")
    os.makedirs(args.out, exist_ok=True)
    t0 = time.time()
    C = pd.read_csv(os.path.join(args.compare, "per_species.tsv"), sep="\t", index_col=0, keep_default_na=False, na_values=[""])
    tmn = C[C.membrane_new.fillna("") != ""]
    add_belt_addendum(C, v7dir, lib, args.nproc, args.cache or args.out)
    L5, L5a = score_l5(C, tmn.index)
    for d in (L5, L5a):
        d["pass"] = d.pop("pass_")
    old = json.load(open(os.path.join(args.compare, "library_compare.json")))
    items = {k: old[k]["pass"] for k in ("L1", "L2", "L3", "L4", "L6", "L7")}
    res = dict(compare=args.compare, build=args.build, library=lib, v7=v7dir,
               addendum="docs/unification_criteria/B3_library_addendum_2026-11-17.md", L5=L5, L5a=L5a,
               other_items_from_compare=items, pass_frozen_L5=all(items.values()) and L5["pass"],
               pass_addendum=all(items.values()) and L5a["pass"], seconds=round(time.time() - t0))
    cols = [c for c in C.columns if c.startswith(("z_belt", "belt_from", "belt_a_from", "n_belt", "map_note"))] + ["membrane_new", "tm_route_new"]
    C.loc[tmn.index, cols].to_csv(os.path.join(args.out, "l5_per_species.tsv"), sep="\t")
    json.dump(res, open(os.path.join(args.out, "l5_rescore.json"), "w"), indent=1, default=str)
    log(json.dumps({k: res[k] for k in ("pass_frozen_L5", "pass_addendum")}), f"L5 fail {L5['fail_new']}/{L5['fail_v7']}",
        f"L5a fail {L5a['fail_new']}/{L5a['fail_v7']}")


# ------------------------------------------------------------------ verdict of a build
def verdict(args):
    import pandas as pd
    bdir = args.build
    out = os.path.join(bdir, "validate", "compare_B3")
    os.makedirs(out, exist_ok=True)
    args.out = out
    meta = json.load(open(os.path.join(bdir, "library", "library_meta.json")))
    lib = compare(args)
    items, fail = {}, []
    for k in ("L1", "L2", "L3", "L4", "L5a", "L6", "L7"):                  # L5a: addendum 2026-11-17
        items[k] = lib[k]["pass"]
    items_info = {"L5 (frozen definition)": lib["L5"]["pass"], "L5 (addendum 2026-11-17)": lib["L5a"]["pass"]}
    steps = meta.get("steps", [])
    items["B1"] = bool(args.complete) and meta.get("mode") == "build" and "tm_route" in steps and "adaptive_trees" in steps
    man = pd.read_csv(os.path.join(bdir, "library", "library_manifest.tsv"), sep="\t", keep_default_na=False)
    items["B4"] = bool((man.origin != "unknown").all() and meta.get("results", {}).get("inputs"))
    nr = os.path.join(bdir, "validate", "no_regression.tsv")
    if args.ref_suite:                                                # suite of this build against the B3 reference
        nr = os.path.join(out, "suite", "no_regression.tsv")
        if not os.path.exists(nr):
            import subprocess
            subprocess.run([sys.executable, "-m", "wholecell", "validate", "--cell", "jurkat", "--build", bdir, "--step", "suite",
                            "--out", os.path.join(out, "suite"), f"ref_suite={args.ref_suite}"], cwd=REPO, check=False)
    regress, c9 = [], []
    if os.path.exists(nr):
        N = pd.read_csv(nr, sep="\t", keep_default_na=False)
        st = "status" if "status" in N.columns else N.columns[-1]
        regress = N[N[st] == "regress"].to_dict("records")
        c9 = N[N.check == "C9"].to_dict("records") if "check" in N.columns else []
    exempt = []
    if args.b6_expected:                     # reference = v7repro2: the B6 EXPECTED changes (B6_geometry.md) are not B3's
        keys = ("density.cytosol", "perinuclear", "nuc_cyto_ratio", "phi.cytosol", "er_lumen", "matrix", "golgi_lumen",
                "endo", "depends_on_geometry")
        exempt = [r for r in regress if r.get("check") in ("C11", "C16") and any(k in str(r.get("metric", "")).lower() for k in keys)]
        regress = [r for r in regress if r not in exempt]
    items["B2"] = os.path.exists(nr) and not regress
    items["B3"] = os.path.exists(nr) and not [r for r in c9 if r.get(st) == "regress"]
    fail = [k for k, v in items.items() if not v]
    V = dict(step="B3", criteria="docs/unification_criteria/B3_library.md",
             addendum="docs/unification_criteria/B3_library_addendum_2026-11-17.md", l5_both=items_info, build=bdir, ref_suite=args.ref_suite,
             no_regression=nr, items=items, regressions=regress[:50], exempt_b6_expected=exempt, c9=c9, library=lib,
             verdict="ACCEPTED" if not fail else "NOT ACCEPTED", failing=fail,
             note="B2 regressions are listed for review against the EXPECTED list of the criteria (seed-level pack noise "
                  "within the standard tolerances is not a regression by construction)")
    json.dump(V, open(os.path.join(out, "verdict.json"), "w"), indent=1, default=str)
    log("verdict", V["verdict"], fail)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["species", "compare", "verdict", "rescore"])
    ap.add_argument("--compare", help="rescore: an existing compare folder (per_species.tsv, library_compare.json)")
    ap.add_argument("--build"); ap.add_argument("--out"); ap.add_argument("--other"); ap.add_argument("--cache")
    ap.add_argument("--ref-suite"); ap.add_argument("--complete", type=int, default=1)
    ap.add_argument("--nproc", type=int, default=3)
    ap.add_argument("--v7dir"); ap.add_argument("--v7man")
    ap.add_argument("--b6-expected", action="store_true", help="reference is v7repro2: exempt the B6 EXPECTED changes")
    a = ap.parse_args()
    {"species": make_species, "compare": compare, "verdict": verdict, "rescore": rescore}[a.what](a)


if __name__ == "__main__":
    main()
