"""Organisation layer (unification plan §2b P6; switch `pack.organisation`, defaults in wholecell/pack/switches.py).

Ported from the yeast packer (plugins/yeast/engine/pipeline/organisation.py @ e3b7448): a weak, data-driven hypothesis
on top of random packing. Compartments and cell totals are unchanged; only where in its regions a type goes and next to
whom.

  network   every protein of the STRING network (combined score >= min_score, edges weighted by score / 1000) gets its
            personalized PageRank profile (random walk with restart, continuation probability `damping`; profile
            `nodiag`: the seed's own entry dropped, renormalised). A packed type takes the mean profile of its subunits;
            the Pearson correlation r of two types' profiles is their network similarity (matrix R).
  field     (large scale, counts step) R ~ E E^T from its `field_k` leading eigenvectors (E = loadings). With field_k
            smooth unit-variance random fields F_k (random Fourier features, correlation length field_l_nm), a type's log
            weight in a tile is beta x E_i . F(tile centre); copies are split over tiles proportional to the region
            measure x exp(log weight), with a crowding correction (wholecell/pack/requests.py: tile_logmod,
            crowding_factors).
  partners  (small scale, pack step) contact partners from AlphaFold 3 pool pairs (optional, off in published builds)
            and loose partners (r >= r_loose); CSR arrays used by the kernel's partner attraction (pack3d.pack_tile):
            a copy is proposed next to an already placed partner with probability min(p_max, sum_p a_p),
            a_p = strength x A_kind x n_p / (k_half + n_p).

Step: python -m wholecell pack --cell C --build B --step organisation  ->  pack/organisation/{network.npz, partners.npz,
network.json}. network.npz: ids (network type order), tix (types with a profile), Z (float32, standardised profiles; for
the viewer and checks), E (float64 loadings: the float32 storage of the yeast build loses exactness, A2).
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os

import numpy as np

CONTACT, LOOSE = 0, 1


# ------------------------------------------------------------------ inputs
def read_table(path, **kw):
    """csv / tsv (optionally .gz) / parquet as a DataFrame"""
    import pandas as pd
    if path.endswith(".parquet"):
        return pd.read_parquet(path, **kw)
    base = path[:-3] if path.endswith(".gz") else path
    return pd.read_csv(path, sep="," if base.endswith(".csv") else "\t", **kw)


# origin: yeast organisation.py:_uniprot_maps and lib_subunits (gene -> UniProt) @ e3b7448
# Changed: one function for both maps; the columns are arguments.
def token_map(table, columns, value="Entry"):
    """{token: value} over the given columns of a UniProt table (tokens split on ';' and whitespace; first row wins)"""
    out = {}
    cols = [table[c].fillna("").astype(str).values for c in columns]
    for k, e in enumerate(table[value].values):
        for c in cols:
            for x in c[k].replace(";", " ").split():
                out.setdefault(x, e)
    return out


def aliases_map(path, source="UniProt_AC"):
    """{STRING protein id without the taxon prefix: first alias of `source`} from a STRING aliases file"""
    out = {}
    with gzip.open(path, "rt") if path.endswith(".gz") else open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            sid, alias, src = line.rstrip("\n").split("\t")[:3]
            if src == source:
                out.setdefault(sid.split(".", 1)[1], alias)
    return out


# origin: yeast organisation.py:string_edges @ e3b7448 (file, id map and score cut as arguments)
def string_edges(links, id_map, min_score):
    """STRING links as UniProt pairs (a < b), combined score >= min_score (max over STRING ids mapping to a pair)"""
    import pandas as pd
    s = pd.read_csv(links, sep=" ")
    s = s[s.combined_score >= min_score]
    a = s.protein1.str.split(".", n=1).str[1].map(id_map); b = s.protein2.str.split(".", n=1).str[1].map(id_map)
    d = pd.DataFrame(dict(a=a, b=b, score=s.combined_score.values)).dropna()
    d = d[d.a < d.b]
    return d.groupby(["a", "b"], as_index=False).score.max()


# origin: yeast organisation.py:pool_pairs @ e3b7448 (no cache; the cached pair table is accepted as input)
def pool_pairs(path):
    """best (lowest) chain-pair PAE min-min per UniProt pair: a parquet (a, b, pae) as cached by the yeast build, or
    the pool predictions' summary_models parquet (af3_id1, af3_id2, chain_pair_pae_min_min)"""
    import pyarrow.parquet as pq
    cols = pq.ParquetFile(path).schema_arrow.names
    if {"a", "b", "pae"} <= set(cols):
        return pq.read_table(path, columns=["a", "b", "pae"]).to_pandas()
    import pyarrow as pa
    import pyarrow.compute as pc
    t = pq.read_table(path, columns=["af3_id1", "af3_id2", "chain_pair_pae_min_min"])
    a, b = pc.utf8_upper(t["af3_id1"]), pc.utf8_upper(t["af3_id2"])
    lt = pc.less(a, b)
    t = pa.table(dict(a=pc.if_else(lt, a, b), b=pc.if_else(lt, b, a), pae=t["chain_pair_pae_min_min"]))
    g = t.group_by(["a", "b"]).aggregate([("pae", "min")]).to_pandas().rename(columns={"pae_min": "pae"})
    return g[g.a != g.b]


# origin: yeast organisation.py:pool_edges @ e3b7448 (pairs and paralog table as arguments)
def pool_edges(pairs, paralogs, pae_max):
    """confident pool pairs without paralog artefacts: drop (A, B) if a paralog A' of A scores better with B"""
    g = pairs[pairs.pae <= pae_max].copy()
    para = {}
    if paralogs is not None:
        for x, y in zip(paralogs.gene1_uniprot, paralogs.gene2_uniprot):
            para.setdefault(x, set()).add(y); para.setdefault(y, set()).add(x)
    best = {(a, b): p for a, b, p in zip(g.a, g.b, g.pae)}
    best.update({(b, a): p for (a, b), p in list(best.items())})
    keep = []
    for a, b, p in zip(g.a, g.b, g.pae):
        bad = any(best.get((a2, b), 99) < p for a2 in para.get(a, ())) or any(best.get((a, b2), 99) < p for b2 in para.get(b, ()))
        keep.append(not bad)
    return g[np.array(keep, bool)]


# origin: yeast organisation.py:lib_subunits @ e3b7448
# Changed: (ingredient, type) pairs and the gene map are arguments (yeast: expected_totals x ingredients.csv; core: the
# pack table, ingredient = type).
def lib_subunits(pairs, table, gene_map):
    """{type: set of UniProt accessions of its subunits}: the ingredient's uniprot field (';'-separated), else its
    genes mapped through gene_map. pairs: iterable of (ingredient id, type id); table: ing_id, uniprot, genes"""
    rows = table.set_index("ing_id")
    out = {}
    for ing_id, lib in pairs:
        if ing_id not in rows.index:
            continue
        acc, genes = rows.at[ing_id, "uniprot"], rows.at[ing_id, "genes"]
        if isinstance(acc, str) and acc:
            us = [x for x in acc.split(";") if x]
        elif isinstance(genes, str):
            us = [gene_map[g] for g in genes.split(";") if g in gene_map]
        else:
            us = []
        if us:
            out.setdefault(lib, set()).update(us)
    return out


# ------------------------------------------------------------------ network propagation
def _adjacency(edges, nodes):
    ix = {u: i for i, u in enumerate(nodes)}
    return edges.a.map(ix).values, edges.b.map(ix).values, edges.score.values / 1000.0


# origin: yeast organisation.py:ppr_profiles @ e3b7448
# Changed: rows for `need` only when the network has more than dense_max nodes (Neumann series on a sparse matrix:
# the dense inverse of the human network does not fit in memory); same profiles to ~1e-15.
def ppr_profiles(edges, damping, profile="nodiag", need=None, dense_max=8000):
    """personalized PageRank profiles: (nodes, X, row) with X[row[i]] the profile of node i (rows sum to 1)"""
    nodes = sorted(set(edges.a) | set(edges.b)); N = len(nodes)
    a, b, w = _adjacency(edges, nodes)
    if N <= dense_max:
        A = np.zeros((N, N))
        A[a, b] = w; A[b, a] = w
        W = A / A.sum(0)[None, :]                                  # column-stochastic transition matrix
        X = ((1 - damping) * np.linalg.inv(np.eye(N) - damping * W)).T   # row i: visits of walks restarting at i
        row = np.arange(N)
        seeds = np.arange(N)
    else:
        import scipy.sparse as sp
        ix = {u: i for i, u in enumerate(nodes)}
        seeds = np.array(sorted({ix[u] for u in (need or ()) if u in ix}), np.int64)
        A = sp.coo_matrix((np.r_[w, w], (np.r_[a, b], np.r_[b, a])), shape=(N, N)).tocsr()
        A.sum_duplicates()
        W = (A @ sp.diags(1.0 / np.asarray(A.sum(0)).ravel())).tocsr()
        term = np.zeros((N, len(seeds))); term[seeds, np.arange(len(seeds))] = 1.0 - damping
        acc = term.copy()
        while len(seeds) and term.max() > 1e-17:
            term = damping * (W @ term); acc += term
        X = acc.T.copy()
        row = -np.ones(N, np.int64); row[seeds] = np.arange(len(seeds))
    if profile == "nodiag":
        X[np.arange(len(seeds)), seeds] = 0.0; X /= X.sum(1, keepdims=True)
    return nodes, X, row


# origin: yeast organisation.py:type_profiles @ e3b7448 (row index of the profile matrix as argument)
def type_profiles(lib_ids, subunits, nodes, X, row):
    """mean profile over each type's subunits, standardised: (type indices tix, Z) with Z[k] . Z[l] = Pearson r"""
    ix = {u: i for i, u in enumerate(nodes)}
    tix, rows = [], []
    for i, lib in enumerate(lib_ids):
        us = [row[ix[u]] for u in subunits.get(lib, ()) if u in ix]
        if us:
            tix.append(i); rows.append(X[us].mean(0))
    Z = np.array(rows).reshape(len(rows), len(nodes))
    Z -= Z.mean(1, keepdims=True); Z /= np.linalg.norm(Z, axis=1, keepdims=True)
    return np.array(tix, np.int64), Z


# origin: yeast organisation.py:field_loadings @ e3b7448
def field_loadings(R, k):
    """E (types x k) with E E^T the best rank-k approximation of the correlation matrix R (positive part)"""
    lam, U = np.linalg.eigh(R)
    top = np.argsort(lam)[::-1][:k]
    lam = np.clip(lam[top], 0, None)
    return U[:, top] * np.sqrt(lam)[None, :], lam


# origin: yeast organisation.py:partners @ e3b7448 (constants as arguments)
def partners(lib_ids, subunits, pool, tix, R, pool_pae_max=8.0, r_loose=0.4, max_partners=16):
    """partner lists per type as CSR arrays over lib_ids indices: off (n+1), idx, strength (0..1), kind (CONTACT/LOOSE)"""
    lib_ids = [str(x) for x in lib_ids]
    by_u = {}
    for i, lib in enumerate(lib_ids):
        for u in subunits.get(lib, ()):
            by_u.setdefault(u, []).append(i)
    row = np.full(len(lib_ids), -1, np.int64); row[tix] = np.arange(len(tix))
    r_of = lambda i, j: float(R[row[i], row[j]]) if row[i] >= 0 and row[j] >= 0 else 0.0
    shares = lambda i, j: i == j or bool(subunits.get(lib_ids[i], set()) & subunits.get(lib_ids[j], set()))
    P = [dict() for _ in lib_ids]
    if pool is not None and len(pool):
        deg = {}
        for a, b in zip(pool.a, pool.b):
            deg[a] = deg.get(a, 0) + 1; deg[b] = deg.get(b, 0) + 1
        # contact: pool pairs, PAE term x network support, hubs of the pool graph pull less
        for a, b, p in zip(pool.a, pool.b, pool.pae):
            if a not in by_u or b not in by_u:
                continue
            s0 = float(np.clip((pool_pae_max - p) / (pool_pae_max - 2.0), 0.05, 1.0)) / np.sqrt(max(1.0, max(deg[a], deg[b]) / 20.0))
            for i in by_u[a]:
                for j in by_u[b]:
                    if shares(i, j):
                        continue
                    s = s0 * (1 + max(0.0, r_of(i, j))) / 2
                    for x, y in ((i, j), (j, i)):
                        if P[x].get(y, (CONTACT, -1))[1] < s:
                            P[x][y] = (CONTACT, s)
    # loose: the most correlated network neighbours
    Rm = R.copy(); np.fill_diagonal(Rm, -1)
    for k, i in enumerate(tix):
        js = np.nonzero(Rm[k] >= r_loose)[0]
        js = js[np.argsort(-Rm[k, js])][:2 * max_partners]
        for kk in js:
            j = int(tix[kk])
            if j in P[i] or shares(i, j):
                continue
            P[i][j] = (LOOSE, float(0.2 + 0.8 * (Rm[k, kk] - r_loose) / (1 - r_loose)))
    off = [0]; idx = []; st = []; kd = []
    for i in range(len(lib_ids)):
        top = sorted(P[i].items(), key=lambda kv: (kv[1][0], -kv[1][1]))[:max_partners]   # contact first, then stronger
        for j, (kind, s) in top:
            idx.append(j); st.append(s); kd.append(kind)
        off.append(len(idx))
    return dict(off=np.array(off, np.int64), idx=np.array(idx, np.int64), strength=np.array(st, np.float64),
                kind=np.array(kd, np.int64))


# origin: yeast organisation.py:params @ e3b7448
def params(o):
    """kernel attraction parameters prm[0:7]: A_contact, A_loose, K_half, P_max, attract tries, loose gap min, max (nm)"""
    return np.array([o["a_contact"], o["a_loose"], o["k_half"], o["p_max"], o["attract_tries"],
                     o["loose_gap_nm"][0], o["loose_gap_nm"][1]], np.float64)


def remap_partners(pr, net_ids, lib_ids):
    """partner CSR over the network's type order -> CSR over lib_ids (types and partners missing from lib_ids dropped;
    order within a row kept). Identity when the two orders are equal."""
    net_ids = [str(x) for x in net_ids]; lib_ids = [str(x) for x in lib_ids]
    if net_ids == lib_ids:
        return {k: pr[k] for k in ("off", "idx", "strength", "kind")}
    pos = {u: i for i, u in enumerate(lib_ids)}
    to_lib = np.array([pos.get(u, -1) for u in net_ids], np.int64)
    at = {u: k for k, u in enumerate(net_ids)}
    off = [0]; idx = []; st = []; kd = []
    for u in lib_ids:
        k = at.get(u)
        if k is not None:
            for e in range(pr["off"][k], pr["off"][k + 1]):
                j = to_lib[pr["idx"][e]]
                if j >= 0:
                    idx.append(j); st.append(pr["strength"][e]); kd.append(pr["kind"][e])
        off.append(len(idx))
    return dict(off=np.array(off, np.int64), idx=np.array(idx, np.int64), strength=np.array(st, np.float64),
                kind=np.array(kd, np.int64))


# ------------------------------------------------------------------ step
def _md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def resolve(cfg, bdir, v):
    """a path: frozen input name, absolute path, or relative to the build directory"""
    if v is None:
        return None
    fi = cfg.get("frozen_inputs") or {}
    if v in fi:
        return fi[v]["path"]
    return v if os.path.isabs(v) else os.path.join(bdir, v)


def network_ids(cfg, bdir, o, pc, lib_paths):
    """the network's type order: subunits.ids (npy), else pack.requests.lib_order (legacy_x5), else the pack-table
    ingredients that have a library entry (table order)"""
    if o["subunits"].get("ids"):
        return [str(x) for x in np.load(resolve(cfg, bdir, o["subunits"]["ids"]), allow_pickle=True)]
    rq = pc.get("requests") or {}
    if rq.get("method") == "legacy_x5":
        return [str(x) for x in np.load(resolve(cfg, bdir, rq["lib_order"]), allow_pickle=True)]
    t = read_table(lib_paths["pack_table"])
    return [str(i) for i in t.ing_id if i in lib_paths["library"]]


def build(cfg, bdir, o, ids, pack_table, out_dir, log=print):
    """network.npz, partners.npz and network.json for the network type order `ids`"""
    st_c = o["string"]
    links = resolve(cfg, bdir, st_c["links"])
    for key, tag in (("taxon", "{}."), ("version", "v{}")):
        if st_c.get(key) is not None and tag.format(st_c[key]) not in os.path.basename(links):
            raise ValueError(f"pack.organisation.string.{key} = {st_c[key]} does not match the links file {links}")
    ic = st_c["ids"]
    if ic.get("aliases"):
        id_map = aliases_map(resolve(cfg, bdir, ic["aliases"]), ic.get("source", "UniProt_AC"))
    else:
        id_map = {}                                       # a UniProt column of loci or of full STRING ids ("9606.ENSP...;")
        for k, v in token_map(read_table(resolve(cfg, bdir, ic["table"])), ic["columns"], ic.get("value", "Entry")).items():
            pre, _, rest = k.partition(".")
            id_map.setdefault(rest if rest and pre.isdigit() else k, v)
    gc = o["genes"]
    gene_map = token_map(read_table(resolve(cfg, bdir, gc["table"])), gc["columns"], gc.get("value", "Entry")) if gc.get("table") else {}
    sc = o["subunits"]
    table = read_table(resolve(cfg, bdir, sc["table"]) if sc.get("table") else pack_table)
    if sc.get("lib_map"):
        lm = read_table(resolve(cfg, bdir, sc["lib_map"]))
        pairs = list(lm[["ing_id", "lib_id"]].drop_duplicates().itertuples(index=False, name=None))
    else:
        pairs = [(i, i) for i in table.ing_id]
    subunits = lib_subunits(pairs, table, gene_map)
    edges = string_edges(links, id_map, int(st_c["min_score"]))
    need = sorted({u for lib in ids for u in subunits.get(lib, ())})
    nodes, X, row = ppr_profiles(edges, float(o["damping"]), o["profile"], need, int(o["dense_max_nodes"]))
    tix, Z = type_profiles(ids, subunits, nodes, X, row)
    del X
    R = Z @ Z.T
    E, lam = field_loadings(R, int(o["field_k"]))
    pool = None
    inputs = {"string_links": links}
    if o.get("pool_scores"):
        pp = resolve(cfg, bdir, o["pool_scores"])
        par = read_table(resolve(cfg, bdir, o["paralogs"]), usecols=["gene1_uniprot", "gene2_uniprot"]).dropna() if o.get("paralogs") else None
        pool = pool_edges(pool_pairs(pp), par, float(o["pool_pae_max"]))
        inputs.update(pool_scores=pp, paralogs=resolve(cfg, bdir, o["paralogs"]) if o.get("paralogs") else None)
    pr = partners(ids, subunits, pool, tix, R, float(o["pool_pae_max"]), float(o["r_loose"]), int(o["max_partners"]))
    os.makedirs(out_dir, exist_ok=True)
    from wholecell.build import atomic_save
    atomic_save(os.path.join(out_dir, "partners.npz"), lambda tmp: np.savez(tmp, ids=np.array(ids), **pr))
    atomic_save(os.path.join(out_dir, "network.npz"), lambda tmp: np.savez_compressed(
        tmp, ids=np.array(ids), tix=tix, Z=Z.astype(np.float32), E=E.astype(np.float64)))
    npart = np.diff(pr["off"])
    off_diag = R[np.triu_indices(len(R), 1)]
    stats = dict(string_links=int(len(edges)), proteins=len(nodes), types=len(ids), types_with_profile=int(len(tix)),
                 pairs_r_ge_loose=float((off_diag >= float(o["r_loose"])).mean()) if len(off_diag) else 0.0,
                 field_share=float(lam.sum() / max(len(R), 1)), pool_pairs=int(len(pool)) if pool is not None else 0,
                 types_with_partners=int((npart > 0).sum()), contact_links=int((pr["kind"] == CONTACT).sum()),
                 loose_links=int((pr["kind"] == LOOSE).sum()), ppr="dense" if len(nodes) <= int(o["dense_max_nodes"]) else "neumann")
    info = dict(params={k: v for k, v in o.items() if k not in ("string", "genes", "subunits")},
                string={k: v for k, v in st_c.items()}, inputs={k: (dict(path=v, md5=_md5(v)) if v else None) for k, v in inputs.items()},
                stats=stats)
    with open(os.path.join(out_dir, "network.json.tmp"), "w") as f:
        json.dump(info, f, indent=1, default=str)
    os.replace(os.path.join(out_dir, "network.json.tmp"), os.path.join(out_dir, "network.json"))
    log(f"organisation: {stats['string_links']} STRING links >= {st_c['min_score']} over {len(nodes)} proteins "
        f"(PPR {stats['ppr']}, damping {o['damping']}); {len(tix)} of {len(ids)} types with a profile (pairs with r >= "
        f"{o['r_loose']}: {stats['pairs_r_ge_loose']:.2%}); {E.shape[1]} field components hold {stats['field_share']:.0%} of the "
        f"correlation; pool pairs {stats['pool_pairs'] if pool is not None else 'off'}; {stats['types_with_partners']} types "
        f"with partners ({stats['contact_links']} contact, {stats['loose_links']} loose links)")
    return pr
