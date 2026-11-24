"""Check a YCT2 web export with MESH / GHOST / POLY sections against its build (read-only on both).

usage: python scripts/dev/check_web_sections.py --build B --web DIR [--ref REF_WEB] [--sample N] [--workers 8] [--json OUT]

  identity  per tile of REF_WEB (default: <build>/export/web): the PALETTE entries the INST words use (prefix of the new
            palette) and the INST section are byte-identical; the tile sets of REF are all present
  mesh      on a deterministic sample of mesh tiles (all when --sample 0): |SDF| of every decoded vertex for its own
            membrane (max / median / p99), the share of faces whose normal points to positive SDF (cytoplasmic / N side),
            total area per membrane (all mesh tiles, from the decoded quantised vertices) vs geometry/measures.json
            (membrane_areas_summary, Jurkat) and the pack tiling's mem_<membrane> areas (pack/tiles.parquet)
  poly      decoded points per tile = index n_poly_beads; deltas fit int8 (decoder sanity)
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)
from wholecell.export import meshpoly as MP, yct2  # noqa: E402

S = {}


def _rows(web):
    j = json.load(open(os.path.join(web, "tiles", "index.json")))
    c = j["tile_columns"]
    return j, {tuple(r[:3]): dict(zip(c, r)) for r in j["tiles"]}


def _raw(web, r):
    with open(os.path.join(web, "tiles", r["shard"]), "rb") as f:
        f.seek(r["off"]); return gzip.decompress(f.read(r["gz_bytes"]))


def _sec(raw):
    h = np.frombuffer(raw[:128], "<u4"); ns = int(h[12])
    tab = np.frombuffer(raw[128:128 + 16 * ns], "<u4").reshape(ns, 4)
    return h, {int(k): (raw[o:o + l], int(c)) for k, o, l, c in tab}


def _identity(keys):
    bad = []; n = 0
    for k in keys:
        a = _raw(S["ref"], S["ref_rows"][k]); r = S["new_rows"].get(k)
        if r is None:
            bad.append((k, "missing")); continue
        b = _raw(S["web"], r)
        ha, sa = _sec(a); hb, sb = _sec(b)
        na = int(ha[11])
        if sa[2][0] != sb[2][0]:
            bad.append((k, "INST")); continue
        if sb[1][0][:8 * na] != sa[1][0]:
            bad.append((k, "PALETTE")); continue
        if int(hb[9]) != int(ha[9]):
            bad.append((k, "n_inst"))
        n += int(ha[9])
    return len(keys), n, bad[:20], len(bad)


def _mesh(keys):
    cell = S.get("cell")
    if cell is None:
        from wholecell.geometry import cellgeom
        cell = S["cell"] = cellgeom.load(S["pkg"], verify="quick")
    origin = np.asarray(S["origin"]); T = S["tile_nm"]; qo, qr = S["q_org"], S["q_range"]
    out = dict(area={}, sdf=[], normal_ok=0, normal_n=0, poly_ok=0, poly_bad=0, tiles=0)
    for k in keys:
        r = S["new_rows"][k]; raw = _raw(S["web"], r); h, sec = _sec(raw)
        lo = origin + np.array(k) * T
        if 4 in sec or 5 in sec or 6 in sec:
            npt = 0
            for kind in (4, 5, 6):
                if kind in sec:
                    res = MP.unpack_polylines(sec[kind][0], sec[kind][1], lo, MP.SECTIONS[kind]["poly_q"], extra=kind == 6)
                    pl = res[0] if kind == 6 else res
                    npt += sum(len(p) for p in pl)
            out["poly_ok" if npt == r["n_poly_beads"] else "poly_bad"] += 1
        if 7 not in sec:
            continue
        V, F, FM = MP.decode_mesh(sec[7][0], sec[7][1], bool(h[13] & 1), lo, qo, qr)
        assert len(V) == r["n_vert"] and len(F) == r["n_tri"], (k, len(V), r["n_vert"])
        out["tiles"] += 1
        tri = V[F]; n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]); ar = 0.5 * np.linalg.norm(n, axis=1)
        for c in np.unique(FM):
            out["area"][int(c)] = out["area"].get(int(c), 0.0) + float(ar[FM == c].sum())
        if not S["full"] and (hash(k) % S["every"]):           # SDF / normals on the sample only
            continue
        vm = np.zeros(len(V), np.int64); vm[F.ravel()] = np.repeat(FM, 3)
        sd = np.asarray(cell.mems64(V), np.float64)[np.arange(len(V)), vm]
        out["sdf"].append(np.abs(sd))
        good = ar > 1e-3
        cc = tri[good].mean(1); u = n[good] / (2 * ar[good])[:, None]; fm = FM[good].astype(np.int64)
        sp = np.asarray(cell.mems64(cc + 0.5 * u), np.float64)[np.arange(len(cc)), fm]
        sm = np.asarray(cell.mems64(cc - 0.5 * u), np.float64)[np.arange(len(cc)), fm]
        out["normal_ok"] += int((sp > sm).sum()); out["normal_n"] += len(cc)
    out["sdf"] = np.concatenate(out["sdf"]) if out["sdf"] else np.zeros(0)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", required=True); ap.add_argument("--web", required=True); ap.add_argument("--ref")
    ap.add_argument("--sample", type=int, default=2000, help="mesh tiles checked for SDF / normals (0 = all)")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--json")
    a = ap.parse_args()
    os.environ.setdefault("NUMBA_NUM_THREADS", "1")
    t0 = time.time()
    ref = a.ref or os.path.join(a.build, "export", "web")
    jn, new_rows = _rows(a.web); jr, ref_rows = _rows(ref)
    pj = json.load(open(os.path.join(a.build, "geometry", "package.json")))
    mesh_keys = [k for k, r in new_rows.items() if r["n_tri"] > 0]
    all_keys = sorted(new_rows)
    S.update(web=a.web, ref=ref, new_rows=new_rows, ref_rows=ref_rows, pkg=pj["package_dir"], origin=jn["grid"]["origin_nm"],
             tile_nm=jn["tile_nm"], q_org=jn["q_org"], q_range=jn["q_range"], full=a.sample == 0,
             every=max(1, len(mesh_keys) // max(a.sample, 1)))
    bt = jn.get("block_tiles", 5)
    nb = {tuple(x // bt for x in k) for k in new_rows}
    rk = sorted(k for k in ref_rows if tuple(x // bt for x in k) in nb)     # a partial export (blocks=) checks its blocks
    from wholecell.geometry import cellgeom
    S["cell"] = cellgeom.load(S["pkg"], verify="quick")       # before the fork: workers share the geometry (Jurkat ~5 GB)
    chunks = lambda L, n: [L[i::n] for i in range(n)]
    res = {}
    import multiprocessing as mp
    with ProcessPoolExecutor(a.workers, mp_context=mp.get_context("fork")) as ex:
        idn = list(ex.map(_identity, chunks(rk, 4 * a.workers)))
        ms = list(ex.map(_mesh, chunks(all_keys, 4 * a.workers)))
    res["identity"] = dict(ref_tiles=sum(x[0] for x in idn), ref_instances=sum(x[1] for x in idn),
                           mismatches=sum(x[3] for x in idn), examples=[e for x in idn for e in x[2]][:10],
                           new_tiles=len(new_rows), new_only_tiles=len(set(new_rows) - set(ref_rows)))
    sd = np.concatenate([m["sdf"] for m in ms]) if ms else np.zeros(0)
    area = {}
    for m in ms:
        for c, v in m["area"].items():
            area[c] = area.get(c, 0.0) + v
    names = {m["code"]: m["name"] for m in jn.get("mesh", {}).get("membranes", [])}
    meas = json.load(open(os.path.join(a.build, "geometry", "measures.json"))).get("membrane_areas_summary", {}) or {}
    tp = os.path.join(a.build, "pack", "tiles.parquet")
    tiling = {}
    if os.path.exists(tp):
        import pandas as pd
        M = pd.read_parquet(tp)
        tiling = {c[4:]: float(M[c].sum()) / 1e6 for c in M.columns if c.startswith("mem_")}
    ar = {}
    for c, v in sorted(area.items()):
        nm = names.get(c, str(c)); ref_a = (meas.get(nm) or {}).get("area_um2"); t_a = tiling.get(nm)
        ar[nm] = dict(mesh_um2=round(v / 1e6, 3), geometry_um2=ref_a, ratio_geometry=round(v / 1e6 / ref_a, 4) if ref_a else None,
                      pack_tiling_um2=round(t_a, 3) if t_a else None, ratio_tiling=round(v / 1e6 / t_a, 4) if t_a else None)
    nok = sum(m["normal_ok"] for m in ms); nn = sum(m["normal_n"] for m in ms)
    res["mesh"] = dict(mesh_tiles=sum(m["tiles"] for m in ms), vertices_checked=int(len(sd)),
                       abs_sdf_nm=dict(max=round(float(sd.max()), 4) if len(sd) else None,
                                       median=round(float(np.median(sd)), 4) if len(sd) else None,
                                       p99=round(float(np.percentile(sd, 99)), 4) if len(sd) else None),
                       normals_to_positive_sdf=dict(ok=nok, n=nn, share=round(nok / nn, 5) if nn else None), area=ar)
    res["poly"] = dict(tiles_ok=sum(m["poly_ok"] for m in ms), tiles_bad=sum(m["poly_bad"] for m in ms))
    res["seconds"] = round(time.time() - t0, 1)
    print(json.dumps(res, indent=1))
    if a.json:
        json.dump(res, open(a.json, "w"), indent=1)


if __name__ == "__main__":
    main()
