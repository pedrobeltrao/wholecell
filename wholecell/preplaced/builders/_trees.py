# origin: yeast pipeline/structlib.py:sphere_tree (+ BEAD_R) @ e3b7448 (unchanged)
"""Coarse sphere trees of residue beads (k-means clusters, 95th-percentile radius + bead radius), used by the
pre-placed builders for their own walk / obstacle shapes (e.g. the ribosome's ~30-sphere walk tree)."""
import numpy as np

BEAD_R = 0.42  # nm, residue bead radius used for rendering / sphere-tree padding


def sphere_tree(xyz, per=30, pad=BEAD_R, seed=0):
    from sklearn.cluster import KMeans
    n = len(xyz)
    if n == 0:
        return np.zeros((0, 3)), np.zeros(0)
    k = max(1, int(np.ceil(n / per)))
    if k == 1:
        c = xyz.mean(0, keepdims=True)
        r = np.array([np.percentile(np.linalg.norm(xyz - c, axis=1), 95) + pad])
        return c, r
    km = KMeans(n_clusters=k, n_init=1, random_state=seed, max_iter=50).fit(xyz)
    cs, rs = [], []
    for i in range(k):
        p = xyz[km.labels_ == i]
        if len(p) == 0:
            continue
        c = p.mean(0)
        cs.append(c); rs.append(np.percentile(np.linalg.norm(p - c, axis=1), 95) + pad)
    return np.array(cs), np.array(rs)
