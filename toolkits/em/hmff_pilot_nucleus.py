# origin: scratch_checks/hmff_pilot_nucleus.py (legacy workspace, 2026-11-15)
# ported, not yet restructured, not run: absolute paths replaced by toolkits/em/_paths (inputs inp(), outputs out() in the work folder)
# HMFF pilot: fit a Helfrich membrane to the nuclear envelope of target cell 1 (jrc_jurkat-1), via Mosaic + FreeDTS v2.2.
# Initial mesh = marching cubes of the smoothed nucleus_seg mask (16 nm, isotropic); density = ne_pred (high = NE).
# Notes: Mosaic's equilibrate_fit is unusable here (500 Botsch iterations on the 1.4M-vertex MC mesh; Trimem not
# installed), so we decimate + remesh ourselves. FreeDTS v2.2 names the energy FreeDTS1.0_HMFF (Mosaic writes _MDFF).
# usage: hmff_pilot_nucleus.py STEPS [EDGE_nm] [XI]   (run with conda_envs/mosaic)
import os, re, sys, time, json, subprocess
import numpy as np, tensorstore as ts, mrcfile, open3d as o3d
from scipy import ndimage as ndi
from skimage import measure

from toolkits.em._paths import DATA_ROOT as ROOT, EM_LABELS, EXTERNAL, inp, out  # noqa: F401
OUT = os.path.dirname(out("hmff_pilot/inputs.npz"))
B = EM_LABELS
NUC_ID, VOX, PAD = 1, 16.0, 8      # nucleus instance of cell 1; isotropic working voxel (nm); zero-pad (voxels)
STEPS = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
EDGE = float(sys.argv[2]) if len(sys.argv) > 2 else 60.0
XI = float(sys.argv[3]) if len(sys.argv) > 3 else 5.0
THREADS = 16

from toolkits.em.io import read as em_read          # un-flips *_pred layers (they are y-flipped vs *_seg in the release)
def crop(layer, sl):
    return em_read(layer, "s2", sl)

t0 = time.time()
if not os.path.exists(f"{OUT}/inputs.npz"):
    # bbox of nucleus 1 from nucleus_seg s3, mapped to s2 (x2), padded 1 um, clipped to the volume
    nuc3 = ts.open({"driver": "n5", "kvstore": {"driver": "file", "path": B + "nucleus_seg/s3/"}}).result().read().result()
    shp2 = ts.open({"driver": "n5", "kvstore": {"driver": "file", "path": B + "nucleus_seg/s2/"}}).result().shape
    idx = np.argwhere(nuc3 == NUC_ID); lo, hi = idx.min(0) * 2 - 64, idx.max(0) * 2 + 66
    lo = np.maximum(lo, 0); hi = np.minimum(hi, shp2); sl = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    mask = crop("nucleus_seg", sl) == NUC_ID           # x,y,z at 16,16,13.76 nm
    dens = crop("ne_pred", sl).astype(np.float32) / 255.0
    zf = 13.76 / VOX                                   # resample z to isotropic 16 nm
    mask = ndi.zoom(mask.astype(np.float32), (1, 1, zf), order=1) > 0.5
    dens = ndi.zoom(dens, (1, 1, zf), order=1)
    lab, n = ndi.label(mask); mask = lab == (np.argmax(np.bincount(lab.ravel())[1:]) + 1)   # largest component
    mask = np.pad(mask, PAD); dens = np.pad(dens, PAD)  # closes surfaces cut by the crop/volume faces
    origin_nm = np.array([lo[0] * 16, lo[1] * 16, lo[2] * 13.76]) - PAD * VOX
    touches = [bool(lo[i] == 0 or hi[i] == shp2[i]) for i in range(3)]
    np.savez_compressed(f"{OUT}/inputs.npz", mask=mask, dens=dens.astype(np.float16), origin_nm=origin_nm, touches=touches)
    with mrcfile.new(f"{OUT}/ne_density.mrc", overwrite=True) as m:   # mrc axes z,y,x
        m.set_data(np.ascontiguousarray(dens.transpose(2, 1, 0)).astype(np.float32)); m.voxel_size = VOX
z = np.load(f"{OUT}/inputs.npz"); mask = z["mask"]; origin_nm = z["origin_nm"]
print("crop", mask.shape, "touches volume face (x,y,z):", z["touches"], "s", round(time.time() - t0, 1), flush=True)

# baseline: smoothed-mask surface (the mask-based SDF approach), closed by the zero padding
sm = ndi.gaussian_filter(mask.astype(np.float32), 2.0)
v, f, _, _ = measure.marching_cubes(sm, 0.5, spacing=(VOX, VOX, VOX))
base = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v), o3d.utility.Vector3iVector(f[:, ::-1].copy()))
base = base.remove_duplicated_vertices().remove_unreferenced_vertices().remove_degenerate_triangles()
area0 = base.get_surface_area(); vol0 = mask.sum() * VOX**3 * 1e-9
print(f"baseline: {len(base.vertices)} verts, area {area0/1e6:.1f} um2, mask vol {vol0:.1f} um3, "
      f"closed manifold {base.is_edge_manifold(allow_boundary_edges=False)}", flush=True)
o3d.io.write_triangle_mesh(f"{OUT}/nucleus_baseline.ply", base)

# decimate -> Botsch remesh to EDGE nm
from mosaic.meshing.utils import remesh, scale, center_mesh, compute_edge_lengths
from mosaic.formats.writer import write_topology_file
from mosaic.dts.screening import _build_dts_template, generate_screen
t1 = time.time()
# coarse surface for DTS: smooth the mask with sigma 4 voxels (64 nm; removes sub-100 nm pockets that a 60 nm-edge mesh
# cannot represent and that otherwise give folded faces), marching cubes at 32 nm, 20 Botsch iterations (no decimation:
# quadric decimation also introduces flips). Tested: sigma 2 -> 19-25 folded pairs; sigma 4 -> none, area -4%, volume same.
sm4 = ndi.zoom(ndi.gaussian_filter(mask.astype(np.float32), 4.0), 0.5, order=1)
v2, f2, _, _ = measure.marching_cubes(sm4, 0.5, spacing=(2 * VOX,) * 3)
m = o3d.geometry.TriangleMesh(o3d.utility.Vector3dVector(v2), o3d.utility.Vector3iVector(f2[:, ::-1].copy()))
m = m.remove_duplicated_vertices().remove_unreferenced_vertices().remove_degenerate_triangles()
m = remesh(m, EDGE, n_iter=20)
print("remesh done s", round(time.time() - t1, 1), flush=True)
tri = np.asarray(m.triangles); lab = np.asarray(m.cluster_connected_triangles()[0])
m.remove_triangles_by_mask(lab != np.argmax(np.bincount(lab))); m = m.remove_unreferenced_vertices()
e = compute_edge_lengths(m)
print(f"remeshed: {len(m.vertices)} verts, edge mean {e.mean():.1f} p0.1 {np.percentile(e,0.1):.1f} min {e.min():.1f} "
      f"p99.9 {np.percentile(e,99.9):.1f} max {e.max():.1f}; area {m.get_surface_area()/1e6:.1f} um2; "
      f"closed manifold {m.is_edge_manifold(allow_boundary_edges=False)}; s {time.time()-t1:.1f}", flush=True)
o3d.io.write_triangle_mesh(f"{OUT}/nucleus_remeshed.ply", m)

# DTS units: squared edge lengths must lie in [1, 5] (Mosaic's Min_Max_Lenghts); map the 0.1th edge percentile to 1.02
sf = 1.005 / e.min()                     # all squared edges >= 1 in DTS units
data, off = center_mesh(scale(m, sf))
offset = ",".join(str(-float(x)) for x in off)
qfile = f"{OUT}/nucleus_remeshed.q"; write_topology_file(file_path=qfile, data=data)
print("scale", sf, "offset", offset, "edge range DTS", e.min() * sf, e.max() * sf, flush=True)
cfg = _build_dts_template(f"{OUT}/ne_density.mrc", sf, offset,
                          {"kappa": 25.0, "kappa0": 0.0, "xi": XI, "invert_contrast": False},
                          {"steps": STEPS, "threads": THREADS, "output_period": max(STEPS // 10, 1)})
cfg = cfg.replace("FreeDTS1.0_MDFF", "FreeDTS1.0_HMFF")
maxsq = max(5.0, round((e.max() * sf) ** 2 * 1.1, 2))
cfg += "\nMinfaceAngle             = -0.5\n"
cfg = re.sub(r"Min_Max_Lenghts\s*=.*", f"Min_Max_Lenghts          = 1.0 {maxsq}", cfg)
os.system(f"rm -rf {OUT}/screen")
summ = generate_screen(f"{OUT}/screen", qfile, cfg)
run_dir = f"{OUT}/screen/run_1"
print(open(f"{run_dir}/input.dts").read(), flush=True)
env = dict(os.environ, DTS_CMD=f"{EXTERNAL}/FreeDTS/version_2.2/DTS", DTS_NT=str(THREADS))
t2 = time.time()
r = subprocess.run(["bash", f"{run_dir}/run.sh"], env=env, capture_output=True, text=True)
dt = time.time() - t2
print("DTS rc", r.returncode, "seconds", round(dt, 1), "s/1000 steps", round(dt / STEPS * 1000, 1), flush=True)
print(r.stdout[-3000:], r.stderr[-2000:])
json.dump(dict(steps=STEPS, threads=THREADS, edge_nm=EDGE, xi=XI, dts_seconds=dt, rc=r.returncode,
               baseline_area_um2=area0 / 1e6, mask_volume_um3=vol0, n_vertices=len(m.vertices),
               origin_nm=origin_nm.tolist(), scale=sf, offset=offset, qfile=qfile),
          open(f"{OUT}/pilot_summary.json", "w"), indent=1)
