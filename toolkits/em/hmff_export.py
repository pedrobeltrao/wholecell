# origin: scratch_checks/hmff_export.py (legacy workspace, 2026-11-15)
# ported, not yet restructured, not run: absolute paths replaced by toolkits/em/_paths (inputs inp(), outputs out() in the work folder)
# Export the last HMFF trajectory frame as a PLY mesh in nm, in the label (s0 _seg) frame: crop-local nm + crop origin.
from toolkits.em._paths import inp, out  # noqa: F401  (wholecell port)
import os, sys, json, glob, numpy as np, open3d as o3d
run = sys.argv[1]; D = os.path.dirname(inp("hmff_pilot/inputs.npz"))  # noqa: E702
s = json.load(open(f"{D}/pilot_summary.json")); sf = s["scale"]; off = np.array([float(x) for x in s["offset"].split(",")])
origin = np.load(f"{D}/inputs.npz")["origin_nm"]
f = sorted(glob.glob(f"{run}/TrajTSI/*.tsi"), key=lambda p: int(''.join(c for c in p.split('/')[-1] if c.isdigit())))[-1]
L = open(f).read().split("\n")
iv = next(k for k, l in enumerate(L) if l.startswith("vertex")); nv = int(L[iv].split()[1])
it = next(k for k, l in enumerate(L) if l.startswith("triangle")); nt = int(L[it].split()[1])
V = np.array([[float(x) for x in L[iv + 1 + j].split()[1:4]] for j in range(nv)])
T = np.array([[int(x) for x in L[it + 1 + j].split()[1:4]] for j in range(nt)])
Vnm = (V - off) / sf + origin
m = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(Vnm), o3d.utility.Vector3iVector(T))
print(f, "verts", nv, "tris", nt, "area um2", round(m.get_surface_area() / 1e6, 1), "closed manifold", m.is_edge_manifold(allow_boundary_edges=False),
      "volume um3 (signed, divergence)", round(abs(np.einsum("ij,ij->i", Vnm[T[:, 0]], np.cross(Vnm[T[:, 1]], Vnm[T[:, 2]])).sum()) / 6e9, 1))
o3d.io.write_triangle_mesh(sys.argv[2], m)
