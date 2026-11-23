# origin: wholecell/geometry/engines/sdf_sections_v1.py (commit 087a66e). Field / section kernels copied verbatim;
# the hand-coded membrane derivations and label logic (mem_sdfs, label_from) replaced by two DATA tables stored in G
# (membrane table + label rule), read from config/geometry/rules/*.yaml by wholecell/geometry/rules.py.
# Analytic primitive kernels (sdf_prim, sdf_cer, smin, the nucleus ellipsoid) copied verbatim from the vendored yeast
# plugins/yeast/engine/pipeline/cell3d.py @ e3b7448 (A1 yeast half, 2026-11-17).
"""cellgeom_eval_v1: the generic, table-driven geometry evaluator (numba). Spec: docs/reference/cellgeom_v1.md §4.3, §5.

Evaluates sparse-brick SDF fields + procedural sections (G layout of sdf_sections_v1, see that module's docstring for
the header slot map) + analytic primitives (prims: PT type, PQ 8 params, PC class; pidx = the prims to visit) with a
MEMBRANE TABLE (how each membrane SDF is derived from fields / sections / prims / other membranes) and an ordered LABEL
RULE (first row whose conditions all hold gives the label). Both tables are compiled from YAML by
wholecell/geometry/rules.py and appended to G (header slot G_TABLES); a pure-python reference is geometry/eval_ref.py.

G = sdf_sections_v1 G  +  tables section at G[G_TABLES] (G_LAYOUT = LAYOUT_ID = 7101.0):
  TH header        nm (membranes), nrow (membrane-table rows), nr (label rows), nfb (field fallbacks),
                   clearance mode (0 thickness-corrected, 1 plain min |d|), wall membrane, its exempt membrane, 0
  nrow x MROW      kind, a, b, p0..p4      membrane rows (kinds K_* below), evaluated in order
  nr x RROW        code, nc, nc x (term, arg, op, value), zero padded to NCOND conditions
  nfb x 4          field slot, analytic kind (A_*), data offset, 0      (fields absent from the package)
  data             analytic parameters (fallbacks, label shapes), offsets relative to the data start
Membrane kinds. A row of any kind except K_SMOOTH_UNION defines the next membrane code (i = 0, 1, ...); K_SMOOTH_UNION
modifies an earlier membrane in place. "field k" below = sparse-brick field k, or its analytic fallback (nfb rows).
  K_FIELD           a = field slot, b = sign, offset p0 - p1    out = b * field - (p0 - p1)   (p0 = p1 = 0: none;
                                                                 two numbers: the same float ops as cell3d)
  K_WALL            a = membrane j                              out = min(out[j], height above the open face)
  K_SHELL_PORES     a = field slot; p0 width, p1 pore radius,  hp = field + p3; out = |hp| - width/2, a torus of
                    p2 axial window, p3 field offset             radius (pore_r + width/2) around each pore axis
  K_INNER_CRISTAE   a = membrane j, p0 offset                   matrix = out[j] + p0, carved by the cristae section;
                                                                 out = -matrix (mem_sdfs returns matrix)
  K_GOLGI_CIS       a = field slot, p0, p1 = membranes j, k     golgi blob, cisternae inside it away from out[j], out[k]
  K_PRIMS           a = prim class, b = base field slot (-1:    out = min(base, prims of class a); base = field b or
                    none), p0 = reference membrane j (-1), p1    1e9. A cER-patch prim is a sheet at depth p1 - q[5]
                    reference radius                             below membrane j (sdf_cer)
  K_PRIMS_DIFF      a, b = prim classes, p0 = sign              m = max(min prims a, -min prims b); out = p0 * m
                                                                 (p0 = -1: out = -m; mem_sdfs returns m)
  K_SMOOTH_UNION    a = membrane j, b = prim class, p0 = k,     x = p1 * out[j], d = min prims b; if d < x + k:
                    p1 = side (+1 / -1), p2 = 0.25 / k           out[j] = p1 * smin(x, d, k)  (smin_c)
Label terms: T_MEM (out[arg]), T_FIELD (field slot arg; slot beyond the stored fields = 1e9), T_FACE (height above
the open face, 1e9 without), T_PORE (-1 inside the channel of a pore of shell row arg, else +1), T_SPECKLE (-1
inside a speckle sphere, else +1), T_SHAPE (-1 inside the analytic shape at data offset arg, else +1).
op 0: term < value, op 1: term > value.
The operations reduce exactly to those of sdf_sections_v1 and of the vendored yeast cell3d (same float expressions in
the same order; Jurkat criterion B1 and the yeast check tests/check_a1_yeast_geometry.py, level E). Public interface
as sdf_sections_v1 (mem_sdfs, label_from, label_of, label_points, mem_points, ...); n_mems(G) membranes per package.
"""
import numpy as np
import numba as nb

ENGINE_NAME = "cellgeom_eval_v1"
LAYOUT_ID = 7101.0
LAYOUTS = ("sdf_sections_v1", "yeast_cell3d_v1")   # package layouts (manifest engine.name) read with a rule table

# ----------------------------------------------------------------- default code assignment (sdf_sections_v1 layout)
# Label / membrane codes of a package come from its manifest (labels[].engine_code, membranes[].engine_code). The
# constants below are the codes of the sdf_sections_v1 layout (Jurkat, testcell), used by the packer (c3.CYTO, c3.M_PM)
EXT, CYTO, PERI, NUCP, NUCL, IMS, MATRIX, ERL, GOLGI, ENDO, PORE, OUTV, SPK, GIC = range(14)
LABEL_NAMES = ["exterior", "cytosol", "perinuclear", "nucleoplasm", "nucleolus", "mito_IMS", "mito_matrix",
               "ER_lumen", "golgi_lumen", "endolysosomal_lumen", "pore_channel", "outside_EM_volume", "speckle",
               "golgi_intercisternal"]
N_LABELS = 14
M_PM, M_NE, M_ER, M_OM, M_IM, M_GOLGI, M_ENDO, M_OUTER = range(8)
MEM_NAMES = ["PM", "NE", "ER", "MITO_OM", "MITO_IM", "GOLGI", "ENDO", "CELL_OUTER"]
NMEM = 8          # membranes of the sdf_sections_v1 layout (core packer buffers until A6); a package has n_mems(G)
NMEM_MAX = 16     # evaluator buffers (out[]) hold up to NMEM_MAX membranes (yeast_cell3d_v1: 10)
F_CELL, F_NE, F_MITO, F_ER, F_GOLGI, F_ENDO, F_NUCLEOLUS = range(7)
FIELD_NAMES = ["cell", "ne_mid", "mito", "er", "golgi", "endo", "nucleolus"]
SECTION_TYPES = {"cristae_lamellae_v1": "cristae", "speckle_spheres_v1": "speckle", "golgi_cisternae_v1": "golgi_cis"}

# ----------------------------------------------------------------- header slots (as sdf_sections_v1, + G_TABLES)
G_N = 64
G_R, G_CRISTAE, G_FACE, G_ZFACE, G_SPECKLE, G_GOLGI_CIS, G_LAYOUT, G_TABLES = range(8)
S_PERI_W, S_IMS_W, S_PORE_R, S_HALF, S_CR_W, S_CR_GAP, S_CR_JW, S_CIS_CLR, S_PORE_ZWIN, S_PORE_LW = range(8, 18)
S_MEM_HALF = 18
S_CR_RIM = 40     # crista rim rounding radius (B6b, 2026-11-17); 0 / absent = square cap (unchanged)
S_N_MEM_HALF = G_N - 1   # number of MEM_HALF values (0 in G built before 2026-11-17 = NMEM)
PARAM_SLOTS = {"PERI_W": S_PERI_W, "IMS_W": S_IMS_W, "PORE_R": S_PORE_R, "HALF": S_HALF, "CR_W": S_CR_W,
               "CR_GAP": S_CR_GAP, "CR_JW": S_CR_JW, "CIS_CLR": S_CIS_CLR, "PORE_ZWIN": S_PORE_ZWIN,
               "PORE_LW": S_PORE_LW}

# ----------------------------------------------------------------- table encoding
K_FIELD, K_WALL, K_SHELL_PORES, K_INNER_CRISTAE, K_GOLGI_CIS, K_PRIMS, K_PRIMS_DIFF, K_SMOOTH_UNION = range(1, 9)
MEM_KINDS = {"field": K_FIELD, "wall": K_WALL, "shell_pores": K_SHELL_PORES, "inner_cristae": K_INNER_CRISTAE,
             "golgi_cisternae": K_GOLGI_CIS, "prims": K_PRIMS, "prims_difference": K_PRIMS_DIFF,
             "smooth_union": K_SMOOTH_UNION}
T_MEM, T_FIELD, T_FACE, T_PORE, T_SPECKLE, T_SHAPE = 1, 2, 3, 4, 5, 6
OP_LT, OP_GT = 0, 1
MROW = 8          # kind, a, b, p0..p4
NCOND = 4         # conditions per label row
RROW = 2 + 4 * NCOND
TH = 8            # tables header
FBROW = 4         # field fallback rows
CLEAR_THICKNESS, CLEAR_PLAIN = 0, 1
# analytic shapes (field fallbacks, label shapes): kind code, number of data values
A_SPHERE, A_ELLIPSOID, A_BALL = 1, 2, 3
ANALYTIC = {"sphere": (A_SPHERE, 4), "ellipsoid_approx": (A_ELLIPSOID, 15), "ball": (A_BALL, 4)}
# primitive types (prims.npz `type`, mapped by name from manifest prims.types at load); params = 8 floats (cell3d PQ)
P_SPHERE, P_CAPSULE, P_PILL, P_CERPATCH, P_SLOT, P_TORUS = range(6)
PRIM_TYPES = {"sphere": P_SPHERE, "capsule": P_CAPSULE, "pill": P_PILL, "cerpatch": P_CERPATCH, "slot": P_SLOT,
              "torus": P_TORUS}

# ----------------------------------------------------------------- Jurkat v7 defaults (legacy code paths only)
HALF = 2.0
MEM_HALF = np.array([2.275, 2.01, 2.01, 1.975, 1.955, 2.095, 2.245, 2.275])
PORE_R = 45.0
PERI_W = 34.0
IMS_W = 8.0
DEFAULT_PARAMS = {"PERI_W": PERI_W, "IMS_W": IMS_W, "PORE_R": PORE_R, "HALF": HALF, "CR_W": 10.0, "CR_GAP": 12.0,
                  "CR_JW": 10.0, "CIS_CLR": 10.0, "PORE_ZWIN": 150.0, "PORE_LW": 8.0,
                  "MEM_HALF": [float(x) for x in MEM_HALF]}


# ----------------------------------------------------------------- parameter accessors (numba; packer API)
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


# ----------------------------------------------------------------- field / section kernels (verbatim sdf_sections_v1)
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
def golgi_sdf(px, py, pz, G, k, other):
    """GOLGI membrane SDF: blob field k, or (G[G_GOLGI_CIS] > 0) the cisternal layers inside it. other = min of the
    membranes the cisternae keep clear of (mito OM, ER): within CIS_CLR of them the blob field is kept"""
    g = field_sdf(px, py, pz, G, k)
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
    """depth below the plasma membrane (negative outside the cell; field slot F_CELL)"""
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
def shell_pores(px, py, pz, G, pores, k, pw, pr, zwin, foff):
    """shell of width pw around midsurface field k + foff (|f + foff| - pw/2, negative inside the shell), opened by
    pore tori of radius pr + pw/2 around each pore axis (NE: perinuclear lumen, pores of the nuclear envelope).
    foff: yeast stores the ONM midsurface (nucleus field), the shell midsurface is nucleus + PERI_W/2"""
    hp = fval(px, py, pz, G, k)
    if foff != 0.0:
        hp = hp + foff
    d = abs(hp) - 0.5 * pw
    rc = pr + 0.5 * pw
    if abs(hp) < rc + pw:
        j = in_pore_w(px, py, pz, pores, rc, zwin)
        if j >= 0:
            vx = px - pores[j, 0]; vy = py - pores[j, 1]; vz = pz - pores[j, 2]
            z = vx * pores[j, 3] + vy * pores[j, 4] + vz * pores[j, 5]
            rho = np.sqrt(max(vx * vx + vy * vy + vz * vz - z * z, 0.0))
            d = np.sqrt((rho - rc) ** 2 + hp * hp) - 0.5 * pw
    return d


# ----------------------------------------------------------------- analytic primitives (verbatim yeast cell3d)
@nb.njit(cache=True, fastmath=True)
def sdf_cer(px, py, pz, q, pm, Rp):
    """cortical-ER patch: a sheet at constant depth (Rp - q[5]) below the plasma membrane (pm: signed PM distance,
    positive in the cytoplasm), within the angular cap q[0:3], cos q[3], q[4]; half-thickness q[6]"""
    r = np.sqrt(px * px + py * py + pz * pz)
    c = (px * q[0] + py * q[1] + pz * q[2]) / max(r, 1e-9)
    c = min(1.0, max(-1.0, c))
    s = max(np.arccos(c) - q[4], 0.0) * r
    dr = (Rp - q[5]) - pm
    return np.sqrt(s * s + dr * dr) - q[6]


@nb.njit(cache=True, fastmath=True)
def sdf_prim(px, py, pz, t, q):
    if t == P_SPHERE:
        dx = px - q[0]; dy = py - q[1]; dz = pz - q[2]
        return np.sqrt(dx * dx + dy * dy + dz * dz) - q[3]
    if t == P_CAPSULE:
        ax, ay, az = px - q[0], py - q[1], pz - q[2]
        bx, by, bz = q[3] - q[0], q[4] - q[1], q[5] - q[2]
        h = (ax * bx + ay * by + az * bz) / (bx * bx + by * by + bz * bz)
        h = min(max(h, 0.0), 1.0)
        dx = ax - bx * h; dy = ay - by * h; dz = az - bz * h
        return np.sqrt(dx * dx + dy * dy + dz * dz) - q[6]
    if t == P_PILL:   # points within q[7] of a flat disk (centre q0-2, normal q3-5, radius q6)
        vx, vy, vz = px - q[0], py - q[1], pz - q[2]
        z = vx * q[3] + vy * q[4] + vz * q[5]
        rx, ry, rz = vx - z * q[3], vy - z * q[4], vz - z * q[5]
        rho = np.sqrt(rx * rx + ry * ry + rz * rz)
        s = max(rho - q[6], 0.0)
        return np.sqrt(s * s + z * z) - q[7]
    if t == P_TORUS:  # ring: centre q0-2, axis q3-5, major radius q6, tube radius q7
        vx, vy, vz = px - q[0], py - q[1], pz - q[2]
        h = vx * q[3] + vy * q[4] + vz * q[5]
        rho = np.sqrt(max(vx * vx + vy * vy + vz * vz - h * h, 0.0))
        return np.sqrt((rho - q[6]) ** 2 + h * h) - q[7]
    if t == P_SLOT:   # furrow / tube: within q[6] of the rectangle swept from the bottom axis q0-2 -> q3-5 by
        # q[7] along the outward radial direction of its midpoint (axis perpendicular to it)
        mx, my, mz = 0.5 * (q[0] + q[3]), 0.5 * (q[1] + q[4]), 0.5 * (q[2] + q[5])
        m = np.sqrt(mx * mx + my * my + mz * mz); ux, uy, uz = mx / m, my / m, mz / m
        ex, ey, ez = q[3] - q[0], q[4] - q[1], q[5] - q[2]
        L = np.sqrt(ex * ex + ey * ey + ez * ez); ex /= L; ey /= L; ez /= L
        vx, vy, vz = px - q[0], py - q[1], pz - q[2]
        s = min(max(vx * ex + vy * ey + vz * ez, 0.0), L); h = min(max(vx * ux + vy * uy + vz * uz, 0.0), q[7])
        dx, dy, dz = vx - s * ex - h * ux, vy - s * ey - h * uy, vz - s * ez - h * uz
        return np.sqrt(dx * dx + dy * dy + dz * dz) - q[6]
    # P_CERPATCH: spherical shell patch (axis q0-2, cos half-angle q3, half-angle q4, mid radius q5, half thick q6)
    r = np.sqrt(px * px + py * py + pz * pz)
    c = (px * q[0] + py * q[1] + pz * q[2]) / max(r, 1e-9)
    c = min(max(c, -1.0), 1.0)
    th = np.arccos(c)
    s = max(th - q[4], 0.0) * r
    dr = r - q[5]
    return np.sqrt(s * s + dr * dr) - q[6]


@nb.njit(cache=True, fastmath=True)
def smin(a, b, k):
    """polynomial smooth minimum (rounds the union of two distance fields over ~k)"""
    h = max(k - abs(a - b), 0.0) / k
    return min(a, b) - h * h * k * 0.25


@nb.njit(cache=True, fastmath=True)
def smin_c(a, b, k, c):
    """smin with c = 0.25 / k given: min(a, b) - t^2 c, t = max(k - |a - b|, 0) (= smin: h = t / k, h^2 k / 4).
    This is the form LLVM (fastmath) folds cell3d's smin with a compile-time constant k into (t * t * fl(1/(4k)));
    the table reads k at run time, so the folded coefficient is stored in the row (exact reproduction)"""
    t = max(k - abs(a - b), 0.0)
    return min(a, b) - t * t * c


@nb.njit(cache=True, fastmath=True)
def prims_min(px, py, pz, cls, pidx, PT, PQ, PC, ref, refR):
    """min SDF over the prims of class cls among pidx (1e9 if none); cER-patch prims are sheets at depth refR - q[5]
    below the reference membrane value ref (refR < 0: the plain shell patch)"""
    d = 1e9
    for ii in range(pidx.shape[0]):
        i = pidx[ii]
        if PC[i] == cls:
            if PT[i] == P_CERPATCH and refR >= 0.0:
                v = sdf_cer(px, py, pz, PQ[i], ref, refR)
            else:
                v = sdf_prim(px, py, pz, PT[i], PQ[i])
            d = min(d, v)
    return d


# ----------------------------------------------------------------- analytic shapes (field fallbacks, label shapes)
@nb.njit(cache=True, fastmath=True)
def analytic_sdf(px, py, pz, G, kind, o):
    """A_SPHERE (data c xyz, r): |p - c| - r.  A_ELLIPSOID (c xyz, semi-axes abc, rotation rows 9; yeast cell3d
    nuc_ellipsoid, first-order distance k0 (k0 - 1) / k1)"""
    if kind == A_SPHERE:
        dx = px - G[o]; dy = py - G[o + 1]; dz = pz - G[o + 2]
        return np.sqrt(dx * dx + dy * dy + dz * dz) - G[o + 3]
    dx = px - G[o]; dy = py - G[o + 1]; dz = pz - G[o + 2]
    lx = (G[o + 6] * dx + G[o + 7] * dy + G[o + 8] * dz) / G[o + 3]
    ly = (G[o + 9] * dx + G[o + 10] * dy + G[o + 11] * dz) / G[o + 4]
    lz = (G[o + 12] * dx + G[o + 13] * dy + G[o + 14] * dz) / G[o + 5]
    k0 = np.sqrt(lx * lx + ly * ly + lz * lz)
    k1 = np.sqrt((lx / G[o + 3]) ** 2 + (ly / G[o + 4]) ** 2 + (lz / G[o + 5]) ** 2)
    if k1 < 1e-12:
        return -min(G[o + 3], min(G[o + 4], G[o + 5]))
    return k0 * (k0 - 1.0) / k1


@nb.njit(cache=True, fastmath=True)
def in_shape(px, py, pz, G, o):
    """label shape at data offset o: A_BALL (c xyz, r) -> inside if |p - c|^2 < r^2 (yeast nucleolus test)"""
    wx, wy, wz = px - G[o + 1], py - G[o + 2], pz - G[o + 3]
    return wx * wx + wy * wy + wz * wz < G[o + 4] * G[o + 4]


@nb.njit(cache=True, fastmath=True)
def _data0(G):
    o = int(G[G_TABLES])
    return o + TH + MROW * int(G[o + 1]) + RROW * int(G[o + 2]) + FBROW * int(G[o + 3])


@nb.njit(cache=True, fastmath=True)
def fval(px, py, pz, G, k):
    """field slot k: the sparse-brick field, or its analytic fallback when the package lacks the field"""
    o = int(G[G_TABLES])
    nfb = int(G[o + 3])
    if nfb > 0:
        f0 = o + TH + MROW * int(G[o + 1]) + RROW * int(G[o + 2])
        for j in range(nfb):
            q = f0 + FBROW * j
            if int(G[q]) == k:
                return analytic_sdf(px, py, pz, G, int(G[q + 1]), _data0(G) + int(G[q + 2]))
    return field_sdf(px, py, pz, G, k)


# ----------------------------------------------------------------- table-driven membranes and labels
@nb.njit(cache=True, fastmath=True)
def n_mems(G):
    """number of membranes (engine membrane codes) of the package"""
    return int(G[int(G[G_TABLES])])


@nb.njit(cache=True, fastmath=True)
def mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    """Fill out[:n_mems(G)] with the membrane SDFs (negative = lumenal side) by the membrane table; returns the matrix
    SDF of the K_INNER_CRISTAE / K_PRIMS_DIFF row (0.0 if none). pidx: prims to visit (all, or a box selection)."""
    o = int(G[G_TABLES])
    aux = 0.0
    i = 0
    for r in range(int(G[o + 1])):
        q = o + TH + MROW * r
        kind = int(G[q]); a = int(G[q + 1])
        if kind == K_FIELD:
            v = fval(px, py, pz, G, a)
            v = -v if G[q + 2] < 0.0 else v
            if G[q + 3] != 0.0 or G[q + 4] != 0.0:
                v = v - (G[q + 3] - G[q + 4])
            out[i] = v
        elif kind == K_WALL:
            out[i] = min(out[a], face_dist(pz, G))
        elif kind == K_SHELL_PORES:
            out[i] = shell_pores(px, py, pz, G, pores, a, G[q + 3], G[q + 4], G[q + 5], G[q + 6])
        elif kind == K_INNER_CRISTAE:
            matrix = out[a] + G[q + 3]
            if G[G_CRISTAE] > 0.0 and matrix < G[S_CR_W]:
                matrix = max(matrix, -crista_lumen(px, py, pz, G, matrix))
            out[i] = -matrix
            aux = matrix
        elif kind == K_GOLGI_CIS:
            out[i] = golgi_sdf(px, py, pz, G, a, min(out[int(G[q + 3])], out[int(G[q + 4])]))
        elif kind == K_PRIMS:
            b = int(G[q + 2])
            v = fval(px, py, pz, G, b) if b >= 0 else 1e9
            j = int(G[q + 3])
            ref = out[j] if j >= 0 else 0.0
            out[i] = min(v, prims_min(px, py, pz, a, pidx, PT, PQ, PC, ref, G[q + 4] if j >= 0 else -1.0))
        elif kind == K_PRIMS_DIFF:
            m = max(prims_min(px, py, pz, a, pidx, PT, PQ, PC, 0.0, -1.0),
                    -prims_min(px, py, pz, int(G[q + 2]), pidx, PT, PQ, PC, 0.0, -1.0))
            out[i] = -m if G[q + 3] < 0.0 else m
            aux = m
        if kind == K_SMOOTH_UNION:
            d = prims_min(px, py, pz, int(G[q + 2]), pidx, PT, PQ, PC, 0.0, -1.0)
            k = G[q + 3]
            x = -out[a] if G[q + 4] < 0.0 else out[a]
            if d < x + k:
                out[a] = -smin_c(x, d, k, G[q + 5]) if G[q + 4] < 0.0 else smin_c(x, d, k, G[q + 5])
        else:
            i += 1
    return aux


@nb.njit(cache=True, fastmath=True)
def label_from(px, py, pz, G, pores, out, matrix):
    """Label (engine code) by the label rule, given the membrane SDFs out (matrix: unused, kept for the signature).
    Rows are tried in order; the first row whose conditions all hold wins. The last raw field read is remembered."""
    o = int(G[G_TABLES])
    r0 = o + TH + MROW * int(G[o + 1])
    fk = -1; fv = 0.0
    for r in range(int(G[o + 2])):
        q = r0 + RROW * r
        ok = True
        for c in range(int(G[q + 1])):
            t = int(G[q + 2 + 4 * c]); a = int(G[q + 3 + 4 * c])
            if t == T_MEM:
                x = out[a]
            elif t == T_FIELD:
                if a != fk:
                    fk = a
                    fv = fval(px, py, pz, G, a) if n_fields(G) > a else 1e9
                x = fv
            elif t == T_FACE:
                x = face_dist(pz, G)
            elif t == T_PORE:
                m = o + TH + MROW * a
                x = -1.0 if in_pore_w(px, py, pz, pores, G[m + 4] + 0.5 * G[m + 3], G[m + 5]) >= 0 else 1.0
            elif t == T_SPECKLE:
                x = -1.0 if G[G_SPECKLE] > 0.0 and in_speckle(px, py, pz, G) >= 0 else 1.0
            else:   # T_SHAPE
                x = -1.0 if in_shape(px, py, pz, G, _data0(G) + a) else 1.0
            if G[q + 4 + 4 * c] == OP_LT:
                ok = x < G[q + 5 + 4 * c]
            else:
                ok = x > G[q + 5 + 4 * c]
            if not ok:
                break
        if ok:
            return int(G[q])
    return -1


@nb.njit(cache=True, fastmath=True)
def label_of(px, py, pz, G, pores, pidx, PT, PQ, PC, out):
    matrix = mem_sdfs(px, py, pz, G, pores, pidx, PT, PQ, PC, out)
    return label_from(px, py, pz, G, pores, out, matrix)


@nb.njit(cache=True, fastmath=True)
def min_mem_dist_g(out, skip, G):
    """smallest membrane clearance, package rule (tables header):
    thickness-corrected (Jurkat): min_m(|d_m| - mem_half(G, m)) + half(G); so 'min_mem_dist_g(...) < r + half(G)' =
      closer than r to the surface. The wall membrane is a hard wall except for proteins of its exempt membrane.
    plain (yeast legacy): min_m |d_m| (uniform half thickness, applied by the caller)."""
    o = int(G[G_TABLES])
    nm = int(G[o])
    d = 1e9
    if G[o + 4] == CLEAR_PLAIN:
        for m in range(nm):
            if m != skip:
                d = min(d, abs(out[m]))
        return d
    w = int(G[o + 5]); wx = int(G[o + 6])
    for m in range(nm):
        if m != skip and not (m == w and skip == wx):
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
        out = np.zeros(NMEM_MAX)
        lab[i] = label_of(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
    return lab


@nb.njit(cache=True, parallel=True)
def mem_points(P, G, pores, pidx, PT, PQ, PC):
    nm = n_mems(G)
    g = np.zeros((P.shape[0], nm), np.float32)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM_MAX)
        mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        for m in range(nm):
            g[i, m] = out[m]
    return g


@nb.njit(cache=True, parallel=True)
def mem_points64(P, G, pores, pidx, PT, PQ, PC):
    """mem_points in float64 (validation / equality tests)"""
    nm = n_mems(G)
    g = np.zeros((P.shape[0], nm), np.float64)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM_MAX)
        mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        for m in range(nm):
            g[i, m] = out[m]
    return g


@nb.njit(cache=True, parallel=True)
def clearance_points(P, G, pores, pidx, PT, PQ, PC, skip):
    """min_mem_dist_g at points (tests: the packer's membrane clearance)"""
    g = np.zeros(P.shape[0], np.float64)
    for i in nb.prange(P.shape[0]):
        out = np.zeros(NMEM_MAX)
        mem_sdfs(P[i, 0], P[i, 1], P[i, 2], G, pores, pidx, PT, PQ, PC, out)
        g[i] = min_mem_dist_g(out, skip, G)
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
    mh = np.asarray(p["MEM_HALF"], np.float64).ravel()
    assert 1 <= mh.shape[0] <= min(NMEM_MAX, S_N_MEM_HALF - S_MEM_HALF), "MEM_HALF: one value per engine membrane"
    v[S_MEM_HALF:S_MEM_HALF + len(mh)] = mh
    v[S_N_MEM_HALF] = len(mh)
    return v


def params_from_G(G):
    out = {k: float(G[s]) for k, s in PARAM_SLOTS.items()}
    out["CR_RIM"] = float(G[S_CR_RIM])
    n = int(G[S_N_MEM_HALF]) or NMEM
    out["MEM_HALF"] = [float(x) for x in G[S_MEM_HALF:S_MEM_HALF + n]]
    return out


class Tables:
    """compiled rule table (rules.compile_tables): membrane rows (nrow, MROW), label rows (nr, RROW), field fallback
    rows (nfb, FBROW), data blob, n membranes, clearance (mode, wall membrane, exempt membrane)"""
    def __init__(self, mem_rows, label_rows, fallbacks=None, data=None, n_mem=None, clearance=(CLEAR_THICKNESS, -1, -1)):
        self.mem_rows = np.asarray(mem_rows, np.float64).reshape(-1, MROW)
        self.label_rows = np.asarray(label_rows, np.float64).reshape(-1, RROW)
        self.fallbacks = np.zeros((0, FBROW)) if fallbacks is None else np.asarray(fallbacks, np.float64).reshape(-1, FBROW)
        self.data = np.zeros(0) if data is None else np.asarray(data, np.float64).ravel()
        self.n_mem = int(n_mem if n_mem is not None else (self.mem_rows[:, 0] != K_SMOOTH_UNION).sum())
        self.clearance = tuple(float(x) for x in clearance)

    def __iter__(self):          # M, L = compile_tables(...)
        return iter((self.mem_rows, self.label_rows))


def tables_vector(tables, label_rows=None):
    """encode the compiled tables (rules.compile_tables -> Tables; or legacy (mem_rows, label_rows)) as the G tables
    section"""
    if not isinstance(tables, Tables):
        tables = Tables(tables, label_rows, clearance=(CLEAR_THICKNESS, M_OUTER, M_PM)) if label_rows is not None \
            else Tables(*tables, clearance=(CLEAR_THICKNESS, M_OUTER, M_PM))
    t = tables
    assert 1 <= t.n_mem <= NMEM_MAX, f"1..{NMEM_MAX} membranes"
    head = np.zeros(TH)
    head[:7] = [t.n_mem, len(t.mem_rows), len(t.label_rows), len(t.fallbacks), *t.clearance]
    return np.r_[head, t.mem_rows.ravel(), t.label_rows.ravel(), t.fallbacks.ravel(), t.data]


def build_G(R, fields, params=None, face=None, cristae=None, speckle=None, golgi_cis=None, tables=None):
    """assemble G: the sdf_sections_v1 layout (fields in F_* order as (header10, index, data); sections cristae,
    speckle, golgi_cis) followed by the tables section. tables: Tables from rules.compile_tables; None
    = the default rule table (config/geometry/rules/sdf_sections_v1.yaml) with default codes and these params."""
    head = _params_vector(params)
    head[G_R] = float(R)
    head[G_LAYOUT] = LAYOUT_ID
    if face is not None:
        head[G_FACE] = 1.0; head[G_ZFACE] = float(face)
    if tables is None:
        from wholecell.geometry import rules
        present = {"field:" + FIELD_NAMES[k] for k, f in enumerate(fields) if int(np.prod(f[0][4:7])) > 0}
        present |= {"section:" + s for s, v in (("cristae", cristae), ("speckle", speckle), ("golgi_cis", golgi_cis))
                    if v is not None}
        if face is not None:
            present.add("clip:open_face")
        tables = rules.compile_tables(rules.load_rules(rules.default_rules_path("sdf_sections_v1")),
                                      params_from_G(head), present=present)
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
    tv = tables_vector(tables)
    head[G_TABLES] = off; secs.append(tv); off += len(tv)
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
    """split a legacy cell3d G (8-slot header) or an sdf_sections_v1 / cellgeom_eval_v1 G (g_n=G_N) into build_G
    arguments: dict(R, fields, face, cristae, speckle, golgi_cis) (the tables section is not returned)"""
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
    """geometry container: G (with tables), pores (centre on the perinuclear midsurface + outward normal), prims
    (empty unless set_prims). cellgeom.load returns a subclass with .manifest, .vocab, .summary, .package_dir, .rules."""
    def __init__(self, G, pores=None):
        self.G = np.asarray(G, np.float64)
        if self.G.shape[0] <= G_N or self.G[G_LAYOUT] != LAYOUT_ID or self.G[G_TABLES] <= 0:
            raise ValueError("G is not a cellgeom_eval_v1 vector (no rule tables; build it with build_G)")
        self.pores = np.zeros((0, 6)) if pores is None else np.asarray(pores, np.float64).reshape(-1, 6)
        self.PT = np.zeros(0, np.int64); self.PQ = np.zeros((0, 8)); self.PC = np.zeros(0, np.int64)
        self.pidx = np.zeros(0, np.int64)

    def set_prims(self, PT, PQ, PC):
        """analytic primitives (evaluator type codes PRIM_TYPES, 8 params, class codes); pidx = all of them"""
        self.PT = np.ascontiguousarray(PT, np.int64); self.PQ = np.ascontiguousarray(PQ, np.float64).reshape(-1, 8)
        self.PC = np.ascontiguousarray(PC, np.int64); self.pidx = np.arange(len(self.PT), dtype=np.int64)

    @property
    def n_mems(self):
        return n_mems(self.G)

    @property
    def params(self):
        return params_from_G(self.G)

    def labels(self, P):
        return label_points(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC)

    def mems(self, P):
        return mem_points(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC)

    def mems64(self, P):
        return mem_points64(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC)

    def clearance(self, P, skip=-1):
        return clearance_points(np.asarray(P, np.float64), self.G, self.pores, self.pidx, self.PT, self.PQ, self.PC,
                                int(skip))
