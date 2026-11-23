# origin: scratch_checks/count_instances_per_cell.py (legacy workspace, 2026-11-15)
# ported, not yet restructured, not run: absolute paths replaced by toolkits/em/_paths (inputs inp(), outputs out() in the work folder)
# Per-cell instance counts (ribosomes with class; nuclear pores) from s0 instance layers, block-parallel from S3.
# Aggregates voxel count + coordinate sums per global instance id across 256^3 blocks, assigns by centroid to
# the coarse cell labels (data/derived/cells_s3.npz, 32 nm grid) and to the ribo_classified class at the centroid.
from toolkits.em._paths import inp, out  # noqa: F401  (wholecell port)
import tensorstore as ts, numpy as np, pandas as pd, itertools, sys
import multiprocessing as mp
B="https://janelia-cosem-datasets.s3.amazonaws.com/jrc_jurkat-1/jrc_jurkat-1.n5/labels/"
LAYER=sys.argv[1]; BS=256
op=lambda L: ts.open({'driver':'n5','kvstore':{'driver':'http','base_url':B+L+'/s0/'}},context=ts.Context({'cache_pool':{'total_bytes_limit':0}})).result()
dims=op(LAYER).shape
def block(ijk):
    i,j,k=ijk; sl=(slice(i,min(i+BS,dims[0])),slice(j,min(j+BS,dims[1])),slice(k,min(k+BS,dims[2])))
    a=op(LAYER)[sl].read().result()
    if not a.any(): return None
    m=a>0; ids=a[m]; c=np.argwhere(m)+[i,j,k]
    df=pd.DataFrame({"id":ids,"x":c[:,0],"y":c[:,1],"z":c[:,2]})
    g=df.groupby("id"); out=g.size().to_frame("n").join(g[["x","y","z"]].sum())
    if LAYER=="ribo_seg":
        cl=op("ribo_classified")[sl].read().result()[m]
        out=out.join(pd.Series(cl,index=ids).groupby(level=0).agg(lambda s: np.bincount(s,minlength=8).argmax()).rename("cls_mode"))
    return out.reset_index()
if __name__=="__main__":
    grid=list(itertools.product(range(0,dims[0],BS),range(0,dims[1],BS),range(0,dims[2],BS)))
    print(LAYER,"blocks",len(grid),flush=True)
    parts=[]
    with mp.get_context("spawn").Pool(24) as p:
        for n,r in enumerate(p.imap_unordered(block,grid,chunksize=4)):
            if r is not None: parts.append(r)
            if n%1000==0: print(n,flush=True)
    d=pd.concat(parts)
    agg=d.groupby("id").agg(n=("n","sum"),x=("x","sum"),y=("y","sum"),z=("z","sum"),**({"cls":("cls_mode","first")} if "cls_mode" in d else {}))
    for a in "xyz": agg[a]=agg[a]/agg.n
    cells=np.load(inp("cells_s3.npz"))["labels"]
    # s0 voxel (4,4,3.44 nm) -> 32 nm grid index: x/8, y/8, z*3.44/27.52 = z/8
    ix=np.clip((agg[["x","y","z"]].values/8).astype(int),0,np.array(cells.shape)-1)
    agg["cell"]=cells[ix[:,0],ix[:,1],ix[:,2]]
    agg.to_parquet(out(f"{LAYER}_instances.parquet"))
    print("instances",len(agg)); print(agg.groupby("cell").size().rename("count").to_string())
    if "cls" in agg: print(pd.crosstab(agg.cell,agg.cls).to_string())
