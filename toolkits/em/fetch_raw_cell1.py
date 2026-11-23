# origin: scratch_checks/fetch_raw_cell1.py (legacy workspace, 2026-11-15)
# ported, not yet restructured, not run: absolute paths replaced by toolkits/em/_paths (inputs inp(), outputs out() in the work folder)
# Fetch the raw FIB-SEM (fibsem-uint16, s1 = 8x8x6.88 nm) crop around target cell 1 (+0.5 um margin) into a local zarr.
# Crop origin in s1 voxels is stored as attribute 'origin_s1_xyz'. Axes x,y,z (n5 order).
import numpy as np, tensorstore as ts, time, sys
from concurrent.futures import ThreadPoolExecutor
from toolkits.em._paths import DATA_ROOT as ROOT, EM_LABELS, EXTERNAL, inp, out  # noqa: F401
src = ts.open({"driver": "n5", "kvstore": {"driver": "http", "base_url":
      "https://janelia-cosem-datasets.s3.amazonaws.com/jrc_jurkat-1/jrc_jurkat-1.n5/em/fibsem-uint16/s1/"},
      "context": {"cache_pool": {"total_bytes_limit": 0}}}).result()
lab = np.load(inp("cells_s3.npz"))["labels"]
i = np.argwhere(lab == 1); m = np.array([16, 16, 18])            # 32 nm-grid voxels (~0.5 um margin)
lo = np.maximum((i.min(0) - m) * 4, 0); hi = np.minimum((i.max(0) + 1 + m) * 4, src.shape)
lo = lo // 64 * 64                                                 # chunk-aligned start
shape = tuple(int(x) for x in hi - lo)
print("crop s1 lo", lo, "hi", hi, "shape", shape, "GB", np.prod(shape) * 2 / 1e9, flush=True)
dst = ts.open({"driver": "zarr", "kvstore": {"driver": "file", "path": out("raw_cell1_s1.zarr")},
               "metadata": {"shape": shape, "chunks": [128, 128, 128], "dtype": "<u2",
                            "compressor": {"id": "zstd", "level": 3}}, "create": True, "open": True}).result()
open(out("raw_cell1_s1.origin.txt"), "w").write(" ".join(map(str, lo)) + "\n")
slabs = [(x, min(x + 128, shape[0])) for x in range(0, shape[0], 128)]
def job(s):
    a, b = s
    for t in range(5):
        try:
            d = src[lo[0] + a:lo[0] + b, lo[1]:hi[1], lo[2]:hi[2]].read().result()
            dst[a:b].write(d).result(); return a
        except Exception as e:
            print("retry", a, e, flush=True); time.sleep(10)
    raise RuntimeError(f"slab {a} failed")
t0 = time.time()
with ThreadPoolExecutor(6) as ex:
    for k, a in enumerate(ex.map(job, slabs)):
        print(f"slab {a} done ({k+1}/{len(slabs)}) {time.time()-t0:.0f}s", flush=True)
print("RAW_DONE", flush=True)
