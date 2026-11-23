# origin: jurkat/yc/pipeline/jvirtual.py (legacy workspace, 2026-11-15)
# Changed: registered as the `virtual_region` provider "jurkat.cytoskeleton_v7" (wholecell.pack.virtual interface);
# module constants BAND / PLUS_L / PLUS_MIN_L / R_PCM / VIM_D became config params; inputs (MT centrelines,
# cytoskeleton json, NE mesh) come from config; PERI_W and MEM_HALF[NE] are read from the geometry (G) instead of
# cell3d constants; the base label (cytosol) is passed in by the core. Numerics unchanged.
"""Jurkat v7 virtual volume regions: cytosol voxels selected by distance rules.

  mt_lattice  centre band_nm[0]-band_nm[1] (14-26 nm) from an MT axis (MT outer radius 12.5 nm): lattice MAPs (R7)
  mt_plus     the same band, restricted to the plus-end zone: the last plus_l_nm of every MT centreline of length
              >= plus_min_l_nm (shorter segments are mostly trace fragments); plus end by the stage-3 lattice convention
  pcm         within r_pcm_nm of the centriole-pair centre (centrosome layer json `centre_nm`, or the legacy
              cytoskeleton pre-pass json log.centrioles.new_centre_nm)
  vim_cage    cytosol within vim_d_nm of the outer nuclear membrane surface (nearest NE-midsurface mesh vertex)
"""
import json

import numpy as np
from scipy.spatial import cKDTree

from wholecell.registry import register

DEFAULTS = dict(band_nm=[14.0, 26.0], plus_l_nm=400.0, plus_min_l_nm=1000.0, r_pcm_nm=400.0, vim_d_nm=300.0)


class Axes:
    def __init__(self, A, B):
        self.A, self.B = A, B; self.M = 0.5 * (A + B)
        self.hl = float(0.5 * np.linalg.norm(B - A, axis=1).max()) if len(A) else 0.0
        self.tree = cKDTree(self.M) if len(A) else None

    def dist(self, P, thr):
        out = np.full(len(P), np.inf)
        if self.tree is None or not len(P): return out
        d0, _ = self.tree.query(P, k=1, distance_upper_bound=thr + self.hl)
        cand = np.nonzero(np.isfinite(d0))[0]
        if not len(cand): return out
        lst = self.tree.query_ball_point(P[cand], thr + self.hl)
        for i, js in zip(cand, lst):
            js = np.asarray(js); a = self.A[js]; d = self.B[js] - a
            t = np.clip(np.einsum("ij,ij->i", P[i] - a, d) / np.maximum(np.einsum("ij,ij->i", d, d), 1e-9), 0, 1)
            out[i] = float(np.sqrt(((a + t[:, None] * d - P[i]) ** 2).sum(1)).min())
        return out


def _segs(polys):
    A = [p[:-1] for p in polys if len(p) > 1]; B = [p[1:] for p in polys if len(p) > 1]
    return (np.concatenate(A), np.concatenate(B)) if A else (np.zeros((0, 3)), np.zeros((0, 3)))


class CytoskeletonV7:
    """provider object (wholecell.pack.virtual interface)"""

    @staticmethod
    def load_context(params, inputs, geom):
        p = dict(DEFAULTS); p.update(params or {})
        z = np.load(inputs["mt_centerlines"]); P = z["points_model_nm"].astype(np.float64); off = z["offsets"]
        fl = z["flipped_for_lattice"]
        polys = [P[off[i]:off[i + 1]] for i in range(len(off) - 1)]
        plus = []
        for q0, f in zip(polys, fl):
            q = q0[::-1] if f else q0                       # lattice order: minus -> plus
            L = np.r_[0, np.cumsum(np.linalg.norm(np.diff(q, axis=0), axis=1))]
            if L[-1] < p["plus_min_l_nm"]: continue
            plus.append(q[L >= L[-1] - p["plus_l_nm"]])
        J = json.load(open(inputs["cytoskeleton_json"]))
        cen = np.array(J["centre_nm"] if "centre_nm" in J else J["log"]["centrioles"]["new_centre_nm"], np.float64)
        c3 = geom.engine; gp = geom.params
        ctx = dict(ax=Axes(*_segs(polys)), plus=Axes(*_segs(plus)), pcm=cen, n_plus=len(plus), p=p,
                   peri_w=float(gp["PERI_W"]), mh_ne=float(gp["MEM_HALF"][c3.M_NE]),
                   plus_um=float(sum(np.linalg.norm(np.diff(q, axis=0), axis=1).sum() for q in plus) / 1e3))
        if inputs.get("ne_mesh"):                            # NE midsurface vertices (56 nm mesh)
            import trimesh
            m = trimesh.load(inputs["ne_mesh"], process=False)
            ctx["ne_tree"] = cKDTree(np.asarray(m.vertices, np.float64) - geom.centre_nm)
        return ctx

    @staticmethod
    def masks(P, lab, ctx, names, base_codes):
        p = ctx["p"]; BAND = p["band_nm"]
        cyt = np.isin(lab, base_codes) if len(base_codes) > 1 else (lab == base_codes[0])
        out = {}
        idx = np.nonzero(cyt)[0]; Q = P[idx]
        if any(n in names for n in ("mt_lattice", "mt_plus")):
            d = ctx["ax"].dist(Q, BAND[1]); band = (d >= BAND[0]) & (d <= BAND[1])
        for n in names:
            m = np.zeros(len(P), bool)
            if n == "mt_lattice":
                m[idx[band]] = True
            elif n == "mt_plus":
                if band.any():
                    dp = ctx["plus"].dist(Q[band], BAND[1]); bi = idx[band]
                    m[bi[(dp >= BAND[0]) & (dp <= BAND[1])]] = True
            elif n == "pcm":
                m[idx[np.linalg.norm(Q - ctx["pcm"], axis=1) <= p["r_pcm_nm"]]] = True
            elif n == "vim_cage":
                # distance to the ONM surface ~ nearest NE-midsurface vertex - PERI_W/2 - MEM_HALF (vertex spacing ~56 nm)
                d_, _ = ctx["ne_tree"].query(Q, k=1, distance_upper_bound=p["vim_d_nm"] + 0.5 * ctx["peri_w"] + 50.0)
                ne = d_ - 0.5 * ctx["peri_w"] - ctx["mh_ne"]
                m[idx[ne <= p["vim_d_nm"]]] = True
            else:
                raise KeyError(f"jurkat.cytoskeleton_v7: unknown virtual region selector {n}")
            out[n] = m
        return out


register("virtual_region", "jurkat.cytoskeleton_v7")(CytoskeletonV7)
