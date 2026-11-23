# new in wholecell (unification B4, 2026-11-17)
"""Geometry primitives of one cell's evaluator, as used by the common pre-placed builders (`ctx.K`).

The builders (wholecell/preplaced/builders/) hold the placement algorithms; the geometry tests they make inside
their numba kernels come from the cell's evaluator through this namespace:

  K.p3      packing primitives: grid_new / g_insert / g_overlap (collision grid of spheres), quat_mat / rand_quat /
            quat_z_to / transform, check_volume (sphere tree inside allowed labels and clear of membranes),
            emit_spheres, arm_walks (bead walks), TOL
  K.c3      the evaluator: mem_sdfs (signed distance to every membrane), label_from (label code of a point)
  K.mmd     min_mem_dist(out, skip, G): smallest membrane clearance (yeast: uniform half thickness; core: per-membrane)
  K.half    half(G): the half thickness added by the callers of mmd
  K.nmem    length of the `out` buffer of mem_sdfs
  K.lab(name), K.mem(name)   label / membrane codes by common name (cytosol, nucleoplasm, mito_matrix, mito_IMS,
            nuclear_pore_channel; PM, NE, ER, MITO_OM, MITO_IM)

Two kernel sets exist: `yeast(p3, c3)` = the evaluator of the X5 geometry layout (the vendored cell3d / pack3d, so the
yeast pre-pass stays bitwise; retired with the vendored engine) and `core()` = the core evaluator cellgeom_eval_v1
and wholecell.pack.pack3d (Jurkat and every new cell). Numba kernels that use K are made by small factories
(`make(K)`) that close over these objects (numba freezes closure variables as constants); the factories cache one
compiled set per kernel set.
"""
from __future__ import annotations

import numba as nb

_CACHE = {}


class Kernels:
    def __init__(self, name, p3, c3, labels, mems, nmem, mmd, half):
        self.name, self.p3, self.c3 = name, p3, c3
        self._labels, self._mems = dict(labels), dict(mems)
        self.nmem, self.mmd, self.half = nmem, mmd, half

    def lab(self, name):
        if name not in self._labels:
            raise KeyError(f"label {name!r} not known to the {self.name} kernels (known: {sorted(self._labels)})")
        return int(self._labels[name])

    def mem(self, name):
        if name not in self._mems:
            raise KeyError(f"membrane {name!r} not known to the {self.name} kernels (known: {sorted(self._mems)})")
        return int(self._mems[name])

    def mask(self, *names):
        m = 0
        for n in names:
            m |= 1 << self.lab(n)
        return m


def factory(key, K, make):
    """compiled kernels of `make(K)`, one set per (builder key, kernel set)"""
    k = (key, K.name)
    if k not in _CACHE:
        _CACHE[k] = make(K)
    return _CACHE[k]


def yeast(p3, c3):
    """kernels of the yeast X5 geometry layout: the vendored cell3d / pack3d modules of the yeast context"""
    k = ("yeast", id(p3), id(c3))
    if k in _CACHE:
        return _CACHE[k]
    mmd_ = c3.min_mem_dist
    HALF = float(c3.HALF)

    @nb.njit(cache=False)
    def mmd(out, skip, G):
        return mmd_(out, skip)

    @nb.njit(cache=False)
    def half(G):
        return HALF

    labels = dict(cytosol=c3.CYTO, nucleoplasm=c3.NUCP, mito_matrix=c3.MATRIX, mito_IMS=c3.IMS, nuclear_pore_channel=c3.PORE)
    mems = dict(PM=c3.M_PM, ER=c3.M_ER, MITO_OM=c3.M_OM, MITO_IM=c3.M_IM)
    K = Kernels("yeast_x5", p3, c3, labels, mems, c3.NMEM, mmd, half)
    _CACHE[k] = K
    return K


def core():
    """kernels of the core evaluator (cellgeom_eval_v1) and the core packing primitives"""
    if ("core",) in _CACHE:
        return _CACHE[("core",)]
    from wholecell.geometry import engine as c3
    from wholecell.pack import pack3d as p3
    mmd_g, half_g = c3.min_mem_dist_g, c3.half

    @nb.njit(cache=False)
    def mmd(out, skip, G):
        return mmd_g(out, skip, G)

    @nb.njit(cache=False)
    def half(G):
        return half_g(G)

    labels = {n: i for i, n in enumerate(c3.LABEL_NAMES)}
    labels["nuclear_pore_channel"] = labels["pore_channel"]          # common name (yeast PORE)
    mems = {n: i for i, n in enumerate(c3.MEM_NAMES)}
    K = Kernels("core", p3, c3, labels, mems, c3.NMEM_MAX, mmd, half)
    _CACHE[("core",)] = K
    return K
