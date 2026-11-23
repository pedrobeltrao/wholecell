# origin: jurkat/yc/pipeline/cell3d.py (legacy workspace, 2026-11-15)
# Changed: every cell-specific physical constant moved from module constants into reserved header slots of G (slot map
# below); reserved header enlarged from 8 to 64 slots; G-reading numba accessors added for the packer; python helpers to
# assemble G from fields/sections (build_G) and to split a legacy (8-slot) G (split_legacy_G). Numerics identical.
"""sdf_sections_v1: sparse-brick signed-distance fields + procedural sections (Jurkat v7 evaluator, generalised).

Units nm; model frame (origin = cell centre). Membrane SDF sign: negative on the lumenal / non-cytoplasmic side, so the
gradient points to the cytoplasmic (N-) side (= +z of membrane-protein libraries). Mito IM: SDF_IM = -SDF_matrix.

G LAYOUT (float64 vector)
  G[0:G_N]           reserved header, G_N = 64 (legacy cell3d: 8)
  G[G_N]             number of fields nf
  G[G_N+1 : +12*nf]  per field (order F_*): origin x,y,z, spacing h, bricks nbx,nby,nbz, voxels per brick B,
                     inside clamp, outside clamp, ABSOLUTE offset of the brick index, ABSOLUTE offset of brick data
  then index/data blocks and procedural sections (absolute offsets in the header slots below).

HEADER SLOT MAP (G[slot]; 'legacy' = same slot in cell3d)
  0  G_R           radius of a sphere enclosing the cell (legacy)
  1  G_CRISTAE     offset of the cristae section, 0 = none (legacy)
  2  G_FACE        > 0: open boundary (truncation face) at z = G[G_ZFACE] (legacy)
  3  G_ZFACE       z of the open face (model nm) (legacy)
  4  G_SPECKLE     offset of the nuclear-speckle section, 0 = none (legacy)
  5  G_GOLGI_CIS   offset of the Golgi cisternae section, 0 = none (legacy)
  6  G_LAYOUT      layout id = LAYOUT_ID (7001.0) for this engine; 0 in legacy G (accessors fall back to nothing:
                   a legacy G must be converted with split_legacy_G + build_G before use)
  7  (reserved, 0)
  8  S_PERI_W      perinuclear width, ONM-INM midsurface distance            (cell3d PERI_W  = 34.0)
  9  S_IMS_W       mito OM-IM midsurface distance                            (cell3d IMS_W   = 8.0)
  10 S_PORE_R      pore-membrane radius at the NE midplane                   (cell3d PORE_R  = 45.0)
  11 S_HALF        reference bilayer half thickness (packer clearances)      (cell3d HALF    = 2.0)
  12 S_CR_W        crista lumen half separation                              (cell3d CR_W    = 10.0)
  13 S_CR_GAP      lamella end gap inside the boundary IM                    (cell3d CR_GAP  = 12.0)
  14 S_CR_JW       crista junction ribbon half width                         (cell3d CR_JW   = 10.0)
  15 S_CIS_CLR     no Golgi cisternae within this distance of mito/ER        (cell3d CIS_CLR = 10.0)
  16 S_PORE_ZWIN   pore axial half window in in_pore                         (cell3d literal 150.0)
  17 S_PORE_LW     PORE label window below the INM (hp > -PERI_W - this)     (cell3d literal 8.0)
  18..25 S_MEM_HALF  per-membrane exclusion half thickness, order M_*        (cell3d MEM_HALF[8])
  40 S_CR_RIM      crista rim rounding radius (B6b; optional param CR_RIM, 0 = square lamella cap)
  26..63           otherwise reserved (0)
  The defaults (DEFAULT_PARAMS) are the Jurkat v7 values. Packages write their values from the manifest
  (engine_params); the loader (wholecell.geometry.cellgeom) builds G with build_G.

FIELDS (engine slots F_*; a package maps its field names to these): cell (PM midsurface), ne_mid (perinuclear
midsurface), mito (OM), er, golgi (blob), endo, nucleolus. A missing field is an EMPTY field (0 bricks; evaluates to its
outside clamp everywhere = 'outside everywhere'); nucleolus may also be absent (nf = 6).

SECTIONS (formats unchanged from cell3d; offsets inside cristae/speckle sections are relative, the Golgi section holds
a 12-value field header with absolute offsets, rebased by build_G):
  cristae    n, id-grid origin xyz, spacing, dims xyz, n x 12 params (axis xyz, phase, S, junction dir xyz,
             centroid xyz, enabled), dense id grid (-1 none)
  speckle    n, n x (cx, cy, cz, r)
  golgi_cis  W, S, 12-value header of a deep Golgi SDF (absolute offsets), its index, its data

PACKER API (pack3d port). The legacy pack3d used c3.HALF, c3.PERI_W, c3.PORE_R, c3.MEM_HALF (via min_mem_dist) as
numba compile-time globals. In this engine they are per-package values in G. Use INSIDE numba code:
    half(G)  peri_w(G)  pore_r(G)  ims_w(G)  mem_half(G, m)  min_mem_dist_g(out, skip, G)
i.e. replace  c3.HALF -> c3.half(G),  c3.PERI_W -> c3.peri_w(G),  c3.PORE_R -> c3.pore_r(G),
              c3.min_mem_dist(out, s) -> c3.min_mem_dist_g(out, s, G).
Module-level scalars NMEM, M_* (membrane codes), label codes (EXT..GIC), F_* and G_* slots are true engine constants
and stay compile-time globals. The module scalars HALF, PERI_W, PORE_R, IMS_W, MEM_HALF and min_mem_dist(out, skip)
are kept ONLY as the Jurkat v7 defaults for legacy code paths; they are NOT the package values (numba freezes globals
at compile time and its cache does not see changes, so they are never re-set at import). For Jurkat v7 packages the
accessors return exactly these defaults.
Unchanged signatures: mem_sdfs(px,py,pz,G,pores,pidx,PT,PQ,PC,out), label_from(px,py,pz,G,pores,out,matrix),
label_of(...), in_pore(px,py,pz,pores,rc) [axial window = 150 nm default; in_pore_g(px,py,pz,pores,rc,G) uses the
package window], face_dist(pz,G), outer_depth(px,py,pz,G), field_sdf(px,py,pz,G,k), n_fields(G), label_points,
mem_points, Cell(G, pores).
"""
import numpy as np
import numba as nb

ENGINE_NAME = "sdf_sections_v1"
LAYOUT_ID = 7001.0

# ----------------------------------------------------------------- labels / membranes (engine codes, fixed)
EXT, CYTO, PERI, NUCP, NUCL, IMS, MATRIX, ERL, GOLGI, ENDO, PORE, OUTV, SPK, GIC = range(14)
LABEL_NAMES = ["exterior", "cytosol", "perinuclear", "nucleoplasm", "nucleolus", "mito_IMS", "mito_matrix",
               "ER_lumen", "golgi_lumen", "endolysosomal_lumen", "pore_channel", "outside_EM_volume", "speckle",
               "golgi_intercisternal"]
N_LABELS = 14
M_PM, M_NE, M_ER, M_OM, M_IM, M_GOLGI, M_ENDO, M_OUTER = range(8)
MEM_NAMES = ["PM", "NE", "ER", "MITO_OM", "MITO_IM", "GOLGI", "ENDO", "CELL_OUTER"]
NMEM = 8
F_CELL, F_NE, F_MITO, F_ER, F_GOLGI, F_ENDO, F_NUCLEOLUS = range(7)
FIELD_NAMES = ["cell", "ne_mid", "mito", "er", "golgi", "endo", "nucleolus"]
# cellgeom_v1 section type id -> build_G keyword (the procedural sections this engine evaluates)
SECTION_TYPES = {"cristae_lamellae_v1": "cristae", "speckle_spheres_v1": "speckle", "golgi_cisternae_v1": "golgi_cis"}

# ----------------------------------------------------------------- header slots
G_N = 64
G_R, G_CRISTAE, G_FACE, G_ZFACE, G_SPECKLE, G_GOLGI_CIS, G_LAYOUT = range(7)
S_PERI_W, S_IMS_W, S_PORE_R, S_HALF, S_CR_W, S_CR_GAP, S_CR_JW, S_CIS_CLR, S_PORE_ZWIN, S_PORE_LW = range(8, 18)
S_MEM_HALF = 18
S_CR_RIM = 40     # crista rim rounding radius (B6b, 2026-11-17); 0 / absent = square cap (unchanged)
PARAM_SLOTS = {"PERI_W": S_PERI_W, "IMS_W": S_IMS_W, "PORE_R": S_PORE_R, "HALF": S_HALF, "CR_W": S_CR_W,
               "CR_GAP": S_CR_GAP, "CR_JW": S_CR_JW, "CIS_CLR": S_CIS_CLR, "PORE_ZWIN": S_PORE_ZWIN,
               "PORE_LW": S_PORE_LW}

# ----------------------------------------------------------------- Jurkat v7 defaults (legacy code paths only)
HALF = 2.0
MEM_HALF = np.array([2.275, 2.01, 2.01, 1.975, 1.955, 2.095, 2.245, 2.275])
PORE_R = 45.0
PERI_W = 34.0
IMS_W = 8.0
CR_W = 10.0
CR_GAP = 12.0
CR_JW = 10.0
CIS_CLR = 10.0
DEFAULT_PARAMS = {"PERI_W": PERI_W, "IMS_W": IMS_W, "PORE_R": PORE_R, "HALF": HALF, "CR_W": CR_W, "CR_GAP": CR_GAP,
                  "CR_JW": CR_JW, "CIS_CLR": CIS_CLR, "PORE_ZWIN": 150.0, "PORE_LW": 8.0,
                  "MEM_HALF": [float(x) for x in MEM_HALF]}


# ----------------------------------------------------------------- parameter accessors (numba)
@nb.njit(cache=True, fastmath=True)
def half(G):
    return G[S_HALF]


@nb.njit(cache=True, fastmath=True)
def peri_w(G):
    return G[S_PERI_W]


@nb.njit(cache=True, fastmath=True)
def pore_r(G):
    return G[S_PORE_R]


@nb.njit(cache=True, fastmath=True)
def ims_w(G):
    return G[S_IMS_W]


@nb.njit(cache=True, fastmath=True)
def mem_half(G, m):
    return G[S_MEM_HALF + m]


# ----------------------------------------------------------------- numba evaluators
@nb.njit(cache=True, fastmath=True)
def n_fields(G):
    return int(G[G_N]) if G.shape[0] > G_N else 0


@nb.njit(cache=True, fastmath=True)
def field_sdf(px, py, pz, G, k):
    """signed distance field k (negative inside): trilinear within a stored brick, +-clamp elsewhere"""
    return field_sdf_o(px, py, pz, G, G_N + 1 + 12 * k)


@nb.njit(cache=True, fastmath=True)
def field_sdf_o(px, py, pz, G, o):
    """field_sdf for the 12-value header at offset o"""
    h = G[o + 3]; B = int(G[o + 7])
    fx = (px - G[o]) / h; fy = (py - G[o + 1]) / h; fz = (pz - G[o + 2]) / h
    nbx = int(G[o + 4]); nby = int(G[o + 5]); nbz = int(G[o + 6])
    if fx < 0.0 or fy < 0.0 or fz < 0.0 or fx >= nbx * B or fy >= nby * B or fz >= nbz * B:
        return G[o + 9]
    ix = int(fx); iy = int(fy); iz = int(fz)
    bx = ix // B; by = iy // B; bz = iz // B
    bi = int(G[int(G[o + 10]) + (bx * nby + by) * nbz + bz])
    if bi == -1:
        return -G[o + 8]
    if bi == -2:
        return G[o + 9]
    n1 = B + 1
    tx = fx - ix; ty = fy - iy; tz = fz - iz
    b = int(G[o + 11]) + bi * n1 * n1 * n1 + ((ix - bx * B) * n1 + (iy - by * B)) * n1 + (iz - bz * B)
    sx = n1 * n1; sy = n1
    c00 = G[b] + (G[b + sx] - G[b]) * tx
    c01 = G[b + 1] + (G[b + sx + 1] - G[b + 1]) * tx
    c10 = G[b + sy] + (G[b + sx + sy] - G[b + sy]) * tx
    c11 = G[b + sy + 1] + (G[b + sx + sy + 1] - G[b + sy + 1]) * tx
    c0 = c00 + (c10 - c00) * ty; c1 = c01 + (c11 - c01) * ty
    return c0 + (c1 - c0) * tz


@nb.njit(cache=True, fastmath=True)
def golgi_sdf(px, py, pz, G, other):
    """GOLGI membrane SDF: the blob field, or (G[G_GOLGI_CIS] > 0) the cisternal layers inside it. other = min of the
    mito-OM and ER SDFs at p: where the blob overlaps those organelles or lies within CIS_CLR of them, the blob field
    is kept (no cisternae cutting their membranes)"""
    g = field_sdf(px, py, pz, G, F_GOLGI)
    if G[G_GOLGI_CIS] > 0.0 and g < 0.0 and other > G[S_CIS_CLR]:
        o = int(G[G_GOLGI_CIS]); W = G[o]; S = G[o + 1]
        d = -field_sdf_o(px, py, pz, G, o + 2)
        v = d - W
        w = v - S * np.floor(v / S + 0.5)
        g = max(g, abs(w) - W)
    return g


@nb.njit(cache=True, fastmath=True)
def face_dist(pz, G):
    """height above the open (truncation) face; 1e9 when the geometry has no open face"""
    if G[G_FACE] > 0.0:
        return pz - G[G_ZFACE]
    return 1e9


@nb.njit(cache=True, fastmath=True)
def outer_depth(px, py, pz, G):
    """depth below the plasma membrane (negative outside the cell)"""
    return -field_sdf(px, py, pz, G, F_CELL)


@nb.njit(cache=True, fastmath=True)
def in_pore_w(px, py, pz, pores, rc, zwin):
    """index of the pore whose axis (through its centre along the NE normal, |axial| <= zwin) is within rc of p, or -1"""
    for k in range(pores.shape[0]):
        vx = px - pores[k, 0]; vy = py - pores[k, 1]; vz = pz - pores[k, 2]
        z = vx * pores[k, 3] + vy * pores[k, 4] + vz * pores[k, 5]
        if z < -zwin or z > zwin:
            continue
        rx = vx - z * pores[k, 3]; ry = vy - z * pores[k, 4]; rz = vz - z * pores[k, 5]
        if rx * rx + ry * ry + rz * rz < rc * rc:
            return k
    return -1


@nb.njit(cache=True, fastmath=True)
def in_pore(px, py, pz, pores, rc):
    """legacy signature (axial window = default 150 nm); pack3d.near_pore uses this"""
    return in_pore_w(px, py, pz, pores, rc, 150.0)


@nb.njit(cache=True, fastmath=True)
def in_pore_g(px, py, pz, pores, rc, G):
    """in_pore with the package's axial window G[S_PORE_ZWIN]"""
    return in_pore_w(px, py, pz, pores, rc, G[S_PORE_ZWIN])


@nb.njit(cache=True, fastmath=True)
def in_speckle(px, py, pz, G):
    """index of the nuclear speckle sphere containing p, or -1 (section at G[G_SPECKLE])"""
    o = int(G[G_SPECKLE]); n = int(G[o])
    for k in range(n):
        q = o + 1 + 4 * k
        dx = px - G[q]; dy = py - G[q + 1]; dz = pz - G[q + 2]
        if dx * dx + dy * dy + dz * dz < G[q + 3] * G[q + 3]:
            return k
    return -1


@nb.njit(cache=True, fastmath=True)
def crista_lumen(px, py, pz, G, b):
    """signed distance to the intracristal lumen (negative inside) of the mitochondrion at p; b = boundary-IM matrix SDF"""
    o = int(G[G_CRISTAE])
    n = int(G[o]); h = G[o + 4]; nx = int(G[o + 5]); ny = int(G[o + 6]); nz = int(G[o + 7])
    ix = int(np.floor((px - G[o + 1]) / h)); iy = int(np.floor((py - G[o + 2]) / h)); iz = int(np.floor((pz - G[o + 3]) / h))
    if ix < 0 or iy < 0 or iz < 0 or ix >= nx or iy >= ny or iz >= nz:
        return 1e9
    k = int(G[o + 8 + 12 * n + (ix * ny + iy) * nz + iz])
    if k < 0:
        return 1e9
    q = o + 8 + 12 * k
    if G[q + 11] == 0.0:
        return 1e9
    S = G[q + 4]
    t = G[q] * px + G[q + 1] * py + G[q + 2] * pz - G[q + 3]
    u = t - S * np.floor(t / S + 0.5)
    slab = abs(u) - G[S_CR_W]
    lam = max(slab, b + G[S_CR_GAP])
    if G[S_CR_RIM] > 0.0:                    # B6b: rounded rim (half-cylinder of radius r; apex unchanged)
        r = min(G[S_CR_RIM], G[S_CR_W]); e_ = b + G[S_CR_GAP]
        ra = max(slab + r, 0.0); rb = max(e_ + r, 0.0)
        lam = min(-r, lam) + np.sqrt(ra * ra + rb * rb)
    w = abs(G[q + 5] * (px - G[q + 8]) + G[q + 6] * (py - G[q + 9]) + G[q + 7] * (pz - G[q + 10])) - G[S_CR_JW]
    return min(lam, max(slab, w))


@nb.njit(cache=True, fastmath=True)
def mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    """Fill out[NMEM] with signed distances to each membrane (negative = lumenal side). Returns the matrix SDF.
    pidx/PT/PQ/PC are unused (kept for the yeast call signature)."""
    pw = G[S_PERI_W]
    dep = outer_depth(px, py, pz, G)
    out[M_PM] = dep
    out[M_OUTER] = min(dep, face_dist(pz, G))
    hp = field_sdf(px, py, pz, G, F_NE)
    d = abs(hp) - 0.5 * pw
    rc = G[S_PORE_R] + 0.5 * pw
    if abs(hp) < rc + pw:
        k = in_pore_w(px, py, pz, pores, rc, G[S_PORE_ZWIN])
        if k >= 0:
            vx = px - pores[k, 0]; vy = py - pores[k, 1]; vz = pz - pores[k, 2]
            z = vx * pores[k, 3] + vy * pores[k, 4] + vz * pores[k, 5]
            rho = np.sqrt(max(vx * vx + vy * vy + vz * vz - z * z, 0.0))
            d = np.sqrt((rho - rc) ** 2 + hp * hp) - 0.5 * pw
    out[M_NE] = d
    om = field_sdf(px, py, pz, G, F_MITO)
    out[M_OM] = om
    matrix = om + G[S_IMS_W]
    if G[G_CRISTAE] > 0.0 and matrix < G[S_CR_W]:
        matrix = max(matrix, -crista_lumen(px, py, pz, G, matrix))
    out[M_IM] = -matrix
    out[M_ER] = field_sdf(px, py, pz, G, F_ER)
    out[M_GOLGI] = golgi_sdf(px, py, pz, G, min(om, out[M_ER]))
    out[M_ENDO] = field_sdf(px, py, pz, G, F_ENDO)
    return matrix


@nb.njit(cache=True, fastmath=True)
def label_from(px, py, pz, G, pores, out, matrix):
    """Compartment label (engine code) given precomputed membrane SDFs (out) and matrix SDF."""
    if G[G_FACE] > 0.0 and pz < G[G_ZFACE]:
        return OUTV
    if out[M_OUTER] < 0:
        return EXT
    if out[M_NE] < 0:
        return PERI
    pw = G[S_PERI_W]
    hp = field_sdf(px, py, pz, G, F_NE)
    if hp < -0.5 * pw:
        if hp > -pw - G[S_PORE_LW] and in_pore_w(px, py, pz, pores, G[S_PORE_R] + 0.5 * pw, G[S_PORE_ZWIN]) >= 0:
            return PORE
        if n_fields(G) > F_NUCLEOLUS and field_sdf(px, py, pz, G, F_NUCLEOLUS) < 0:
            return NUCL
        if G[G_SPECKLE] > 0.0 and in_speckle(px, py, pz, G) >= 0:
            return SPK
        return NUCP
    if hp < 0:
        return PORE
    if out[M_OM] < 0:
        return MATRIX if matrix < 0 else IMS
    if out[M_ER] < 0:
        return ERL
    if out[M_GOLGI] < 0:
        return GOLGI
    if out[M_ENDO] < 0:
        return ENDO
    if G[G_GOLGI_CIS] > 0.0 and field_sdf(px, py, pz, G, F_GOLGI) < 0.0:
        return GIC
    return CYTO


@nb.njit(cache=True, fastmath=True)
def label_of(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    matrix = mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out)
    return label_from(px, py, pz, G, pores, out, matrix)


@nb.njit(cache=True, fastmath=True)
def min_mem_dist_g(out, skip, G):
    """smallest membrane clearance corrected for each membrane's thickness (package values):
    min_m(|d_m| - mem_half(G, m)) + half(G); so 'min_mem_dist_g(...) < r + half(G)' = closer than r to the surface.
    M_OUTER is a hard wall except for PM proteins (skip = M_PM)."""
    d = 1e9
    for m in range(NMEM):
        if m != skip and not (m == M_OUTER and skip == M_PM):
            d = min(d, abs(out[m]) - G[S_MEM_HALF + m])
    return d + G[S_HALF]


@nb.njit(cache=True, fastmath=True)
def min_mem_dist(out, skip):
    """LEGACY signature with the Jurkat v7 default MEM_HALF/HALF compiled in. Use min_mem_dist_g for packages."""
    d = 1e9
    for m in range(NMEM):
        if m != skip and not (m == M_OUTER and skip == M_PM):
            d = min(d, abs(out[m]) - MEM_HALF[m])
    return d + HALF


@nb.njit(cache=True, parallel=True)
def label_points(P, G, pores, pidx, PT, PQ, PC):
    lab = np.zeros(P.shape[0], np.int8)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM)
        lab[i] = label_of(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
    return lab


@nb.njit(cache=True, parallel=True)
def mem_points(P, G, pores, pidx, PT, PQ, PC):
    g = np.zeros((P.shape[0], NMEM), np.float32)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM)
        mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        for m in range(NMEM):
            g[i, m] = out[m]
    return g


@nb.njit(cache=True, parallel=True)
def mem_points64(P, G, pores, pidx, PT, PQ, PC):
    """mem_points in float64 (validation / equality tests)"""
    g = np.zeros((P.shape[0], NMEM), np.float64)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM)
        mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        for m in range(NMEM):
            g[i, m] = out[m]
    return g


@nb.njit(cache=True, parallel=True)
def field_points(P, G, k):
    """raw field k at points (float64)"""
    g = np.zeros(P.shape[0])
    for i in nb.prange(P.shape[0]):
        g[i] = field_sdf(P[i, 0], P[i, 1], P[i, 2], G, k)
    return g


@nb.njit(cache=True, parallel=True)
def pore_points(P, pores, rc, zwin):
    """in_pore_w at points (validation / declarative label rule)"""
    g = np.zeros(P.shape[0], np.int64)
    for i in nb.prange(P.shape[0]):
        g[i] = in_pore_w(P[i, 0], P[i, 1], P[i, 2], pores, rc, zwin)
    return g


@nb.njit(cache=True, parallel=True)
def speckle_points(P, G):
    g = np.full(P.shape[0], -1, np.int64)
    if G[G_SPECKLE] > 0.0:
        for i in nb.prange(P.shape[0]):
            g[i] = in_speckle(P[i, 0], P[i, 1], P[i, 2], G)
    return g


# ----------------------------------------------------------------- python helpers: assemble / split G
def empty_field(clamp_in=64.0, clamp_out=64.0):
    """a field with no bricks: evaluates to +clamp_out everywhere ('outside everywhere')"""
    return (np.array([0.0, 0.0, 0.0, 16.0, 0, 0, 0, 8, clamp_in, clamp_out], np.float64), np.zeros(0), np.zeros(0))


def _params_vector(params):
    p = dict(DEFAULT_PARAMS)
    p.update(params or {})
    v = np.zeros(G_N, np.float64)
    for k, s in PARAM_SLOTS.items():
        v[s] = float(p[k])
    v[S_CR_RIM] = float(p.get("CR_RIM", 0.0))         # optional (B6b); not in DEFAULT_PARAMS
    mh = np.asarray(p["MEM_HALF"], np.float64)
    assert mh.shape == (NMEM,), "MEM_HALF needs one value per engine membrane"
    v[S_MEM_HALF:S_MEM_HALF + NMEM] = mh
    return v


def params_from_G(G):
    out = {k: float(G[s]) for k, s in PARAM_SLOTS.items()}
    out["CR_RIM"] = float(G[S_CR_RIM])
    out["MEM_HALF"] = [float(x) for x in G[S_MEM_HALF:S_MEM_HALF + NMEM]]
    return out


def build_G(R, fields, params=None, face=None, cristae=None, speckle=None, golgi_cis=None):
    """assemble G. fields: list of (header10, index, data) in F_* order (header10 = origin xyz, h, nbx, nby, nbz, B,
    clamp_in, clamp_out). face: z of the open face or None. cristae / speckle: section arrays (relative offsets).
    golgi_cis: (W, S, (header10, index, data)). Sections are appended in the order cristae, speckle, golgi_cis."""
    head = _params_vector(params)
    head[G_R] = float(R)
    head[G_LAYOUT] = LAYOUT_ID
    if face is not None:
        head[G_FACE] = 1.0; head[G_ZFACE] = float(face)
    nf = len(fields)
    off = G_N + 1 + 12 * nf
    hdr, secs = [], []
    for h10, idx, data in fields:
        h10 = np.asarray(h10, np.float64); assert h10.shape == (10,)
        hdr.append(np.r_[h10, off, off + len(idx)])
        secs += [np.asarray(idx, np.float64), np.asarray(data, np.float64)]
        off += len(idx) + len(data)
    if cristae is not None:
        head[G_CRISTAE] = off; secs.append(np.asarray(cristae, np.float64)); off += len(cristae)
    if speckle is not None:
        head[G_SPECKLE] = off; secs.append(np.asarray(speckle, np.float64)); off += len(speckle)
    if golgi_cis is not None:
        W, S, (h10, idx, data) = golgi_cis
        head[G_GOLGI_CIS] = off
        o2 = off + 2 + 12
        secs.append(np.r_[W, S, np.asarray(h10, np.float64), o2, o2 + len(idx)])
        secs += [np.asarray(idx, np.float64), np.asarray(data, np.float64)]
        off = o2 + len(idx) + len(data)
    G = np.concatenate([head, np.array([nf], np.float64)] + hdr + secs)
    assert G.shape[0] == off
    return G


def pack_fields(R, fields, params=None):
    """legacy helper name: G with the enclosing radius R and the fields (list of (header10, index, data))"""
    return build_G(R, fields, params)


def _field_parts(G, o):
    h10 = G[o:o + 10].copy()
    nb_ = int(G[o + 4] * G[o + 5] * G[o + 6]); B = int(G[o + 7])
    io = int(G[o + 10]); do = int(G[o + 11])
    idx = G[io:io + nb_].copy()
    nst = int((idx >= 0).sum())
    assert nst == (int(idx.max()) + 1 if nb_ and idx.max() >= 0 else 0)
    data = G[do:do + nst * (B + 1) ** 3].copy()
    return h10, idx, data


def split_legacy_G(G, g_n=8):
    """split a legacy cell3d G (8-slot header) into build_G arguments: dict(R, fields, face, cristae, speckle,
    golgi_cis). Works also for an sdf_sections_v1 G with g_n=G_N."""
    G = np.asarray(G, np.float64)
    nf = int(G[g_n])
    fields = [_field_parts(G, g_n + 1 + 12 * k) for k in range(nf)]
    out = dict(R=float(G[0]), fields=fields, face=float(G[3]) if G[2] > 0 else None,
               cristae=None, speckle=None, golgi_cis=None)
    if G[1] > 0:
        o = int(G[1]); n = int(G[o]); nx, ny, nz = (int(G[o + 5]), int(G[o + 6]), int(G[o + 7]))
        out["cristae"] = G[o:o + 8 + 12 * n + nx * ny * nz].copy()
    if G[4] > 0:
        o = int(G[4]); n = int(G[o])
        out["speckle"] = G[o:o + 1 + 4 * n].copy()
    if G[5] > 0:
        o = int(G[5])
        out["golgi_cis"] = (float(G[o]), float(G[o + 1]), _field_parts(G, o + 2))
    return out


# ----------------------------------------------------------------- python container
class Cell:
    """geometry container: G, pores (centre on the perinuclear midsurface + outward normal), empty prims.
    cellgeom.load returns a subclass with .manifest, .vocab, .summary, .package_dir."""
    def __init__(self, G, pores=None):
        self.G = np.asarray(G, np.float64)
        if self.G.shape[0] <= G_N or self.G[G_LAYOUT] != LAYOUT_ID:
            raise ValueError("G is not an sdf_sections_v1 vector (legacy cell3d G? convert with split_legacy_G + build_G)")
        self.pores = np.zeros((0, 6)) if pores is None else np.asarray(pores, np.float64).reshape(-1, 6)
        self.PT = np.zeros(0, np.int64); self.PQ = np.zeros((0, 8)); self.PC = np.zeros(0, np.int64)
        self.pidx = np.zeros(0, np.int64)

    @property
    def params(self):
        return params_from_G(self.G)

    def labels(self, P):
        return label_points(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC)

    def mems(self, P):
        return mem_points(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC)

    def mems64(self, P):
        return mem_points64(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC)
