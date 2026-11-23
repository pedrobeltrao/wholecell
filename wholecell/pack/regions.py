# origin: jurkat/yc/pipeline/jpack_test.py (legacy workspace, 2026-11-15)
# Changed: the region tables (VOL_REG, MEM_REG, SHELL_REG, RNAMES_V7, SOL/TMR, LUMEN_MEM, area factors) are no longer
# module constants: they are built from the cell config `pack.regions` / `pack.region_map` / `pack.mapping_rules` /
# `pack.area_factors` / `pack.thin_lumens` and the geometry (label / membrane names -> engine codes, MEM_HALF from G).
# shell_arrays, regmask, regmem and qflat keep the legacy numerics (qflat copied verbatim).
"""Pack regions from config (docs/ARCHITECTURE.md section 6).

`pack.regions` is an ordered list; the index of a region is its job_region id in the tiles. Kinds, in this order:
  volume    {name, kind: volume, labels: [label names]}
  virtual   {name, kind: virtual, labels: [base label names], provider: <pack.virtual_providers id>, select: <name>}
            (a volume region whose voxels are chosen by a registered `virtual_region` plugin; packed with the shells)
  membrane  {name, kind: membrane, membrane: <engine membrane name>}   (all engine membranes except the outer wall, in
            engine order: the kernel indexes membrane regions by membrane code)
  shell     {name, kind: shell, membrane, labels, side (+1 cytoplasmic side of the SDF, -1 the other, 0 both sides:
            |SDF|), d0, d1, d1_mem_half (add MEM_HALF[membrane] to d1), add_radius (add the library bounding radius to
            d1), pore_limit_nm (centre within this distance of a pore axis; 0 = none), measure (tile measure:
            `area` = the membrane's area, default; `volume` = the voxels of `labels` with d0 <= side*SDF < d1)}
  preplaced {name, kind: preplaced} an id kept for pre-placed instances (no tile measure, never packed)
Volume regions may restrict their voxels by one membrane SDF: `sdf_lt: {membrane, nm}` (SDF < nm; e.g. the outer
60 nm of a cell wall).
Region ids: the position in `pack.regions` is the region id of the counts and of the packed tiles (job_region), in any
order (the yeast legacy ids are not grouped by kind). The packing kernel indexes regions in its own order ("kernel
order", A7): volume-like regions (volume, virtual, and shells measured as volumes, in config order), then the
membranes (engine order), then the area-measured shells; `k_of` / `c_of` map config ids <-> kernel ids (identity for
a table already in kernel order, e.g. Jurkat and testcell). A shell measured as a volume (yeast bands) is packed as a
volume region: seeds = its measured voxels (labels + band), molecule test = its labels, and a partner proposal must
lie within d1 + `pack.volume_shell_tol_nm` (default 4, yeast wc3d.SHELL_TOL) of the membrane (`band_arrays`).
Region mapping of an inventory localisation key: `pack.region_map: {key: {soluble: region, tm: region}}` with
`pack.region_map_default`, then the generic rules `pack.mapping_rules.{tm,soluble}` (e.g. OM vs IM by mito_om_frac,
signal-peptide lumen rule) applied in order.
"""
from __future__ import annotations

import numpy as np

KINDS = ("volume", "virtual", "membrane", "shell", "preplaced")
SHELL_MEASURES = ("area", "volume")


def _codes(names_list, manifest_map):
    """name -> engine code: the package manifest map first, then the engine's own name list"""
    out = {n: i for i, n in enumerate(names_list)}
    out.update(manifest_map or {})
    return out


def label_map(geom) -> dict:
    return _codes(geom.engine.LABEL_NAMES, geom.label_aliases)


def membrane_map(geom) -> dict:
    return _codes(geom.engine.MEM_NAMES, geom.membrane_aliases)


_REQUIRED = {"volume": ("labels",), "virtual": ("labels", "provider"), "membrane": ("membrane",),
             "shell": ("membrane", "labels", "d1"), "preplaced": ()}


def _check(r, lab, mem):
    """one pack.regions entry: kind, required keys, label / membrane names known to the geometry"""
    name = r.get("name", "?")
    if r.get("kind") not in KINDS:
        raise ValueError(f"pack.regions {name}: unknown kind {r.get('kind')!r} (kinds: {', '.join(KINDS)})")
    miss = [k for k in _REQUIRED[r["kind"]] if k not in r]
    if miss:
        raise ValueError(f"pack.regions {name} (kind {r['kind']}): missing {', '.join(miss)}")
    for l in r.get("labels", []) if r["kind"] != "preplaced" else []:
        if l not in lab:
            raise ValueError(f"pack.regions {name}: label {l!r} is not a label of the geometry package")
    ms = [r["membrane"]] if "membrane" in r else []
    if r.get("sdf_lt"):
        if r["kind"] != "volume":
            raise ValueError(f"pack.regions {name}: sdf_lt is for kind volume")
        ms.append(r["sdf_lt"]["membrane"])
    for m in ms:
        if m not in mem:
            raise ValueError(f"pack.regions {name}: membrane {m!r} is not a membrane of the geometry package")
    if r["kind"] == "shell":
        if float(r.get("side", 1.0)) not in (-1.0, 0.0, 1.0):
            raise ValueError(f"pack.regions {name}: side must be 1, -1 or 0 (both sides)")
        if r.get("measure", "area") not in SHELL_MEASURES:
            raise ValueError(f"pack.regions {name}: measure must be one of {SHELL_MEASURES}")


class Regions:
    """region tables for one (config, geometry) pair; attribute names follow the legacy jpack_test module"""

    def __init__(self, pcfg: dict, geom):
        regs = pcfg["regions"]
        c3 = geom.engine
        lab = label_map(geom); mem = membrane_map(geom)
        self.label_codes, self.mem_codes = lab, mem
        mh = np.asarray(geom.params["MEM_HALF"], np.float64)
        self._mh = mh
        for r in regs:
            _check(r, lab, mem)
        self.names = [r["name"] for r in regs]
        if len(set(self.names)) != len(self.names):
            raise ValueError("duplicate region names in pack.regions")
        self.spec = {r["name"]: r for r in regs}
        # kernel order (A7): volume-like regions, membranes (engine order), area-measured shells; preplaced ids: none
        self.not_packable = []
        if any(r["kind"] == "shell" and float(r.get("side", 1.0)) == 0.0 and r.get("measure", "area") == "area" for r in regs):
            self.not_packable.append("shells with side 0 measured by area (a two-sided band must be measured as a volume)")
        vshell = lambda r: r["kind"] == "shell" and r.get("measure", "area") == "volume"
        vol = [r for r in regs if r["kind"] in ("volume", "virtual") or vshell(r)]
        mems = [r for r in regs if r["kind"] == "membrane"]
        shells = [r for r in regs if r["kind"] == "shell" and not vshell(r)]
        self.VOL_REG = [(r["name"], tuple(int(lab[l]) for l in r["labels"])) for r in vol]
        self.MEM_REG = [(r["name"], int(mem[r["membrane"]])) for r in mems]
        if [m for _, m in self.MEM_REG] != list(range(len(self.MEM_REG))):
            raise ValueError("membrane regions must follow the engine membrane order (codes 0..n-1): "
                             f"{[(n, m) for n, m in self.MEM_REG]}")
        self.VIRTUAL = [r["name"] for r in vol if r["kind"] == "virtual"]
        self.virtual_spec = {r["name"]: r for r in vol if r["kind"] == "virtual"}
        self.SHELL_REG = []
        for r in shells:
            m = int(mem[r["membrane"]])
            d1 = float(r["d1"]) + (float(mh[m]) if r.get("d1_mem_half") else 0.0)
            rad = (1 if r.get("add_radius") else 0) | (int(r.get("pore_limit_nm", 0)) << 1)
            self.SHELL_REG.append((r["name"], m, tuple(int(lab[l]) for l in r["labels"]), float(r.get("side", 1.0)),
                                   float(r.get("d0", 0.0)), d1, rad))
        self.NV = len(self.VOL_REG)
        self.NSH0 = self.NV + len(self.MEM_REG)
        self.RNAMES = list(self.names)                     # region id = position in pack.regions
        self.KNAMES = [r[0] for r in self.VOL_REG] + [m[0] for m in self.MEM_REG] + [s_[0] for s_ in self.SHELL_REG]
        kpos = {n: i for i, n in enumerate(self.KNAMES)}
        self.k_of = np.array([kpos.get(n, -1) for n in self.names], np.int64)       # config id -> kernel id (-1: preplaced)
        self.c_of = np.array([self.names.index(n) for n in self.KNAMES], np.int64)  # kernel id -> config id
        self.identity = self.KNAMES == self.names
        self.SHELL_MEMNAME = {s_[0]: self.spec[s_[0]]["membrane"] for s_ in self.SHELL_REG}
        self.SHELL_MEMNAME.update({r["name"]: r["membrane"] for r in vol if vshell(r)})
        # pack order class per kernel region: membranes 0, area shells + virtual volumes 1, label volumes and
        # volume-measured shells 2
        self.pack_class = np.array([2 if (r["kind"] == "volume" or vshell(r)) else 1 if r["kind"] in ("virtual", "shell") else 0
                                    for r in (self.spec[n] for n in self.KNAMES)], np.int64)
        vtol = float(pcfg.get("volume_shell_tol_nm", 4.0))
        self.BAND_M = np.full(len(self.KNAMES), -1, np.int64); self.BAND_W = np.zeros(len(self.KNAMES))
        for i, r in enumerate(vol):
            if vshell(r):
                self.BAND_M[i] = int(mem[r["membrane"]]); self.BAND_W[i] = self.restriction(r["name"])[3] + vtol
        tl = pcfg.get("thin_lumens") or {}
        self.LUMEN_MEM = {k: int(mem[v]) for k, v in tl.items()}
        # region mapping
        self.region_map = pcfg.get("region_map") or {}
        self.map_default = pcfg.get("region_map_default") or {}
        self.rules = pcfg.get("mapping_rules") or {}
        for key, v in self.region_map.items():
            for side in ("soluble", "tm"):
                if side in v and v[side] not in self.spec:
                    raise ValueError(f"region_map {key}.{side}: unknown region {v[side]}")

    def require_core_packer(self, what="the core packer"):
        """the region table must be packable by the kernel (any order; see the module docstring)"""
        if self.not_packable:
            raise ValueError(f"pack.regions cannot run {what}: " + "; ".join(self.not_packable))

    def restriction(self, name):
        """voxel restriction of a volume-like region: (membrane code, side, d0, d1) of a volume-measured shell
        (d0 <= side*SDF < d1, |SDF| for side 0), (membrane code, 'lt', None, nm) of sdf_lt, or None"""
        r = self.spec[name]
        if r["kind"] == "shell":
            m = int(self.mem_codes[r["membrane"]])
            d1 = float(r["d1"]) + (float(self._mh[m]) if r.get("d1_mem_half") else 0.0)
            return m, float(r.get("side", 1.0)), float(r.get("d0", 0.0)), d1
        if r.get("sdf_lt"):
            return int(self.mem_codes[r["sdf_lt"]["membrane"]]), "lt", None, float(r["sdf_lt"]["nm"])
        return None

    def membrane_codes(self, G, c3):
        """pack3d.pack_tile MCODE: engine codes of the membranes named in the membrane-protein tests, from the package
        (not the sdf_sections_v1 constants c3.M_OUTER / M_PM / M_NE, which are other membranes in other packages, e.g.
        yeast code 7 = PEROX): [outer wall CELL_OUTER, the membrane exempt from it (the clearance rule's
        wall_exempt_for in the G tables; else PM when the package has no rule tables), pore membrane NE]; -1 = none"""
        outer = int(self.mem_codes.get("CELL_OUTER", -1))
        exempt = int(self.mem_codes.get("PM", -1))
        if hasattr(c3, "G_TABLES") and int(G[c3.G_TABLES]) > 0:
            o = int(G[c3.G_TABLES])
            exempt = int(G[o + 6]) if int(G[o + 5]) == outer else -1
        return np.array([outer, exempt, int(self.mem_codes.get("NE", -1))], np.int64)

    def band_arrays(self):
        """(BAND_M, BAND_W) per kernel region: membrane and largest |SDF| of a partner proposal in a volume-measured
        shell (-1: none); empty arrays when the table has no such shell (v7 kernel path)"""
        if (self.BAND_M < 0).all():
            return np.zeros(0, np.int64), np.zeros(0)
        return self.BAND_M, self.BAND_W

    # ---------------------------------------------------------------- tile measures
    def measures(self):
        """tiles.parquet measure columns, in region order: [(region index, column, kind of measure, spec)];
        kind 'volume' (voxels of labels, restrictions in spec), 'area' (membrane area; column mem_<membrane>, shared by
        a membrane region and its area-measured shells, listed once), 'zero' (preplaced ids)"""
        out, seen = [], set()
        for ri, r in enumerate(self.spec[n] for n in self.names):
            k = r["kind"]
            if k in ("volume", "virtual") or (k == "shell" and r.get("measure", "area") == "volume"):
                out.append((ri, r["name"], "volume", r))
            elif k == "membrane":
                c = "mem_" + r["membrane"]
                if c not in seen:
                    seen.add(c); out.append((ri, c, "area", r))
            elif k == "preplaced":
                out.append((ri, r["name"], "zero", r))
        return out

    # ---------------------------------------------------------------- legacy jpack_test helpers
    def shell_arrays(self):
        """pack3d.pack_tile shell arguments: n_shell0, SH_MEM, SH_MASK, SH_SIGN, SH_D0, SH_D1, SH_RAD"""
        msk = []
        for s_ in self.SHELL_REG:
            m = 0
            for l in s_[2]:
                m |= 1 << l
            msk.append(m)
        S = self.SHELL_REG
        return (self.NSH0, np.array([s_[1] for s_ in S], np.int64), np.array(msk, np.int64), np.array([s_[3] for s_ in S]),
                np.array([s_[4] for s_ in S]), np.array([s_[5] for s_ in S]), np.array([s_[6] for s_ in S], np.int64))

    def regmem(self):
        return np.array([self.LUMEN_MEM.get(r[0], -1) for r in self.VOL_REG], np.int64)

    def regmask(self):
        out = []
        for _, labs in self.VOL_REG:
            m = 0
            for l in labs:
                m |= 1 << l
            out.append(m)
        return np.array(out + [0] * (len(self.MEM_REG) + len(self.SHELL_REG)), np.int64)

    def tile_column(self, ri: int) -> str:
        """tiles.parquet column measuring region ri (volume, membrane area, or the shell's membrane area)"""
        r = self.spec[self.RNAMES[ri]]
        if r["kind"] in ("volume", "virtual", "preplaced") or (r["kind"] == "shell" and r.get("measure", "area") == "volume"):
            return r["name"]
        return "mem_" + r["membrane"]

    # ---------------------------------------------------------------- mapping of inventory keys
    def _apply_rules(self, rules, reg, row):
        for ru in rules:
            if "if_region" in ru and reg != ru["if_region"]:
                continue
            if "if_region_in" in ru and reg not in ru["if_region_in"]:
                continue
            x = getattr(row, ru["column"])
            op = ru["op"]; v = ru["value"]
            hit = (x >= v) if op == "ge" else (x < v) if op == "lt" else (x > v) if op == "gt" else (x <= v) if op == "le" else None
            if hit is None:
                raise ValueError(f"mapping rule op {op}")
            if hit:
                reg = ru["to"]
        return reg

    def map_tm(self, comp, row):
        reg = self.region_map.get(comp, {}).get("tm", self.map_default.get("tm"))
        return self._apply_rules(self.rules.get("tm", []), reg, row)

    def map_soluble(self, comp, row):
        reg = self.region_map.get(comp, {}).get("soluble", self.map_default.get("soluble"))
        return self._apply_rules(self.rules.get("soluble", []), reg, row)


def qflat(st_off, STC, STR):
    """per library: quaternion (w,x,y,z) rotating its smallest principal axis (sphere-tree centres weighted by r^3) to +z"""
    Q = np.zeros((len(st_off) - 1, 4)); Q[:, 0] = 1.0
    for i in range(len(st_off) - 1):
        c = STC[st_off[i]:st_off[i + 1]]; w = STR[st_off[i]:st_off[i + 1]] ** 3
        if len(c) < 3: continue
        m = (c * w[:, None]).sum(0) / w.sum(); X = c - m
        _, U = np.linalg.eigh((X * w[:, None]).T @ X)           # ascending: U[:, 0] = smallest axis
        R = np.stack([U[:, 1], U[:, 2], U[:, 0]])               # rows: e_mid -> x, e_big -> y, e_small -> z
        if np.linalg.det(R) < 0: R[0] *= -1
        tr = np.trace(R); w_ = np.sqrt(max(1.0 + tr, 1e-12)) / 2
        if w_ > 1e-3:
            Q[i] = [w_, (R[2, 1] - R[1, 2]) / (4 * w_), (R[0, 2] - R[2, 0]) / (4 * w_), (R[1, 0] - R[0, 1]) / (4 * w_)]
        else:
            from scipy.spatial.transform import Rotation
            x, y, z, ww = Rotation.from_matrix(R).as_quat(); Q[i] = [ww, x, y, z]
    return Q
