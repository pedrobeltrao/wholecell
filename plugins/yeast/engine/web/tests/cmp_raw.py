# origin: yeast yeast_cell_dev_bundle/web/tests/cmp_raw.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
import json, numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as Rot
J = json.load(open('/tmp/dec_js.json'))
res = {}
for key, d in J.items():
    ix, iy, iz = map(int, key.split('_'))
    raw = np.load(f'/workspace/yc3d/web_tmp/t_{ix:+03d}_{iy:+03d}_{iz:+03d}.npz')
    rows = np.array(d['rows']); t = rows[:, 0].astype(int); P = rows[:, 1:4]; M = rows[:, 4:13].reshape(-1, 3, 3)
    rt, rp, rq = raw['types'], raw['pos'], raw['quat']
    assert len(rt) == len(t), (key, len(rt), len(t))
    dp, da = [], []
    for ty in np.unique(rt):
        a = np.where(rt == ty)[0]; b = np.where(t == ty)[0]; assert len(a) == len(b)
        tr = cKDTree(P[b]); dist, j = tr.query(rp[a]); dp.append(dist)
        Rr = Rot.from_quat(np.c_[rq[a, 1:], rq[a, 0]]).as_matrix()   # scipy uses x,y,z,w
        Rj = M[b[j]]
        cosang = np.clip((np.einsum('nij,nij->n', Rr, Rj) - 1) / 2, -1, 1)
        da.append(np.degrees(np.arccos(cosang)))
    dp = np.concatenate(dp); da = np.concatenate(da)
    res[key] = dict(n=int(len(t)), pos_err_max_nm=float(dp.max()), rot_err_max_deg=float(da.max()), rot_err_p99_deg=float(np.percentile(da, 99)))
    print(key, res[key])
ok = all(r['pos_err_max_nm'] < 0.01 and r['rot_err_max_deg'] < 0.1 for r in res.values())
print('PASS' if ok else 'FAIL'); json.dump(dict(tiles=res, pass_=ok), open('/workspace/yc3d/web/tests/decoder_roundtrip.json', 'w'), indent=1)
