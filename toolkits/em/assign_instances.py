# origin: jurkat/em/assign_instances.py (legacy workspace, 2026-11-15)
# ported, not yet restructured, not run: absolute paths replaced by toolkits/em/_paths (inputs inp(), outputs out() in the work folder)
# Assign instance centroids (ribo_seg, np_seg; s0 coordinates in the label frame) to cells of data/derived/cells_s3.npz:
# direct lookup, then nearest cell within 500 nm for centroids outside all cells. Rewrites columns cell, cell2.
from toolkits.em._paths import inp, out  # noqa: F401  (wholecell port)
import numpy as np, pandas as pd
from scipy import ndimage as ndi
lab = np.load(inp("cells_s3.npz"))["labels"]
dist, ind = ndi.distance_transform_edt(lab == 0, return_indices=True, sampling=(32, 32, 27.52))
for layer in ("ribo_seg", "np_seg"):
    r = pd.read_parquet(inp(f"{layer}_instances.parquet"))
    ix = np.clip((r[["x", "y", "z"]].values / 8).astype(int), 0, np.array(lab.shape) - 1)
    r["cell"] = lab[tuple(ix.T)]
    near = lab[ind[0][tuple(ix.T)], ind[1][tuple(ix.T)], ind[2][tuple(ix.T)]]; d = dist[tuple(ix.T)]
    r["cell2"] = np.where(r.cell > 0, r.cell, np.where(d <= 500, near, 0))
    r.to_parquet(out(f"{layer}_instances_assigned.parquet"))  # wholecell: new file instead of the legacy in-place rewrite
    print(layer, "total", len(r), "unassigned after 500 nm rule", int((r.cell2 == 0).sum()))
    if layer == "ribo_seg":
        names = {1: "cytosolic", 2: "nuclear", 3: "ER sheet", 4: "ER tube"}
        t = pd.crosstab(r.cell2, r.cls.map(names)); t["total"] = t.sum(axis=1)
        t["ER_bound_%"] = (100 * (t["ER sheet"] + t["ER tube"]) / t.total).round(1); print(t.to_string())
    else:
        print(r.groupby("cell2").size().to_string())
