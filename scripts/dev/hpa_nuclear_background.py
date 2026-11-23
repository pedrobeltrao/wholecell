"""HPA antibody nuclear background (R8c / unification plan section 3c, step B5 input).

Question: cytoplasmic proteins show ~0.3 of their signal in the nucleus pattern of HPA-Jurkat ICC/IF images. Is this
antibody background (b > 0, to be subtracted: observed = b + (1 - b) * true) or a real nuclear pool (b = 0)?

Tests (plan section 3c):
  (1) OpenCell live HEK293T endogenous split-GFP images (no antibody, Cho et al. 2022 Science, public bucket
      s3://czb-opencell, CC BY 4.0) quantified with the SAME segmentation / per-cell metrics as the HPA images;
  (2) HPA-imaged genes with (near) zero Jurkat RNA as a direct background basis;
  (3) genes imaged with >= 2 HPA antibodies.

Image code: the segmentation + per-cell metrics are a verbatim copy of seg_unmix() in
jurkat/validate/imaging_across_versions.py (source sha256 recorded in the outputs), parameterised by pixel size and
background fallback, with a REDUCED unmixing (nucleus = Hoechst/DAPI, nucleolus, cytoplasm = uniform cell-minus-nucleus,
PM ring; all PSF-blurred) that needs no tubulin / ER channel, so it runs identically on HPA and OpenCell. The full
5-basis unmixing is recomputed on HPA fields and checked against data/derived/v8/imaging_eval/hpa_cells.tsv.

Steps (conda_envs/vcell/bin/python scripts/hpa_nuclear_background.py <step>):
  select    stratified OpenCell gene list           -> data/derived/v8/opencell_background/gene_list.tsv
  fetch     best-focus Hoechst+GFP plane + max projection per FOV (HTTP range requests, no full stacks)
            -> data/derived/v8/opencell_background/images/*.npz, fovs.tsv
  quant_oc  segment + quantify OpenCell planes       -> opencell_background/oc_cells.tsv, oc_genes.tsv
  quant_hpa segment + quantify HPA-Jurkat fields     -> hpa_background/hpa_cells_reduced.tsv, hpa_genes_reduced.tsv
  estimate  b with bootstrap CI, decision rule       -> hpa_background/estimate.json, per-gene tables
  all
CPU: NPROC <= 3 (BG_NPROC)."""
import os, sys, io, json, time, hashlib
for _k in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"): os.environ.setdefault(_k, "1")
import numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")

from _devsite import VC as V, RESULTS  # noqa: E402,F401  (release: site paths)
D = f"{V}/data/derived"; REF = f"{V}/data/reference"
OCD = f"{D}/v8/opencell_background"; OUT = f"{D}/v8/hpa_background"; EVAL = f"{D}/v8/imaging_eval"
IMG = f"{REF}/hpa_subcellular/jurkat_images"
SRC = f"{V}/jurkat/validate/imaging_across_versions.py"
BUCKET = "https://czb-opencell.s3.amazonaws.com"
NPROC = min(int(os.environ.get("BG_NPROC", 3)), 3)
PX_HPA, PX_OC = 80.0, 103.0          # nm/pixel: HPA confocal 63x; OpenCell 63x/1.47 + 6.5 um sCMOS -> 103 nm
FWHM = 250.0
BASES = ["nucleus", "nucleolus", "cytosol", "ER", "PM"]
RBASES = ["nucleus", "nucleolus", "cyto", "PM"]
NUC_T = {"nucleoplasm", "nuclear_punctae", "chromatin", "nucleolus_gc", "nucleolus_fc_dfc", "nuclear_membrane"}
ORG_T = {"er", "golgi", "vesicles", "membrane", "mitochondria", "cell_contact", "focal_adhesions", "cytoskeleton",
         "centrosome"}
N_STRATA = dict(cyto_only=120, nonnuclear_organelle=60, nuclear_only=70, mixed=50)
FOV_PER_GENE = 2
SEED = 0


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def sha(f):
    return hashlib.sha256(open(f, "rb").read()).hexdigest()


# ================================================================= segmentation (copy of imaging_across_versions.seg_unmix)
def seg_quant(blue, red, green, yellow, px, full=True, min_nuc_um2=12.0, max_nuc_um2=160.0, drop_border=True, terr_um=4.0):
    """Verbatim seg_unmix segmentation + metrics; adds the reduced unmixing (r_*) and a background fallback for fields
    with < 2 % background pixels (confluent OpenCell fields). full=False skips the 5-basis fit (no tubulin/ER)."""
    from scipy import ndimage as ndi
    from scipy.optimize import nnls
    from skimage.filters import threshold_otsu
    from skimage.segmentation import watershed, expand_labels
    from skimage.feature import peak_local_max
    from skimage.morphology import remove_small_objects, white_tophat, disk
    um = 1000.0 / px; psf = FWHM / 2.3548 / px
    b = ndi.gaussian_filter(blue.astype(np.float32), 1.0)
    if b.max() <= 0: return []
    nuc = ndi.binary_fill_holes(b > threshold_otsu(b))
    nuc = remove_small_objects(nuc, int(8 * um * um))
    dt = ndi.distance_transform_edt(nuc)
    pk = peak_local_max(dt, min_distance=int(1.2 * um), labels=ndi.label(nuc)[0], exclude_border=False)
    mk = np.zeros(nuc.shape, np.int32); mk[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    nl = watershed(-dt, mk, mask=nuc)
    area = np.bincount(nl.ravel()) / um ** 2
    ok = (area >= min_nuc_um2) & (area <= max_nuc_um2); ok[0] = False
    nl = np.where(ok[nl], nl, 0)
    if nl.max() == 0: return []
    sc = lambda x: x.astype(np.float32) / max(np.percentile(x, 99.5), 1e-6)
    s = ndi.gaussian_filter(sc(red) + sc(green), 2.0)
    cm = ndi.binary_fill_holes((s > threshold_otsu(s)) | (nl > 0))
    terr = expand_labels(nl, distance=int(terr_um * um)); cl = np.where(cm, terr, 0)
    bgm = ~ndi.binary_dilation(cm, iterations=int(1 * um))
    few_bg = bgm.sum() < 0.02 * bgm.size
    bgv = lambda x: (np.percentile(x, 1.0) if few_bg else np.median(x[bgm]))
    g = green.astype(np.float32) - bgv(green)
    th = white_tophat(ndi.gaussian_filter(green.astype(np.float32), 0.7), footprint=disk(int(0.6 * um)))
    bsub = lambda x: np.clip(x.astype(np.float32) - bgv(x), 0, None)
    B_ = bsub(blue)
    if full: R_, Y_ = bsub(red), bsub(yellow)
    out = []
    H_, W_ = g.shape
    objs = ndi.find_objects(cl)
    r_rim, r_n = 0.4 * um, 0.24 * um
    for lab, sl in enumerate(objs, start=1):
        if sl is None: continue
        if drop_border and (sl[0].start == 0 or sl[1].start == 0 or sl[0].stop == H_ or sl[1].stop == W_): continue
        p0 = [max(sl[0].start - 3, 0), max(sl[1].start - 3, 0)]
        ss = (slice(p0[0], min(sl[0].stop + 3, H_)), slice(p0[1], min(sl[1].stop + 3, W_)))
        cell = cl[ss] == lab; nucm = (nl[ss] == lab) & cell
        if nucm.sum() < 10 or cell.sum() < nucm.sum() + 20: continue
        depth = ndi.distance_transform_edt(cell)
        own = (terr[ss] == lab) & ~cell
        outd = ndi.distance_transform_edt(~cell)
        nd_out = ndi.distance_transform_edt(~nucm); nd_in = ndi.distance_transform_edt(nucm)
        rim = cell & (depth <= r_rim) & ~nucm
        n_in = nucm & (nd_in > r_n)
        c_in = cell & (depth > r_rim) & (nd_out > r_n)
        gg = g[ss]; tot = gg[cell].sum()
        p99 = np.percentile(gg[cell], 99); eps = 0.02 * max(p99, 1e-6)
        r = dict(cell_area_um2=cell.sum() / um ** 2, nuc_area_um2=nucm.sum() / um ** 2, n_cyto_px=int(c_in.sum()),
                 nuc_frac=float(gg[nucm].sum() / tot) if tot > 0 else np.nan, few_bg=bool(few_bg))
        if c_in.sum() >= 20 and n_in.sum() >= 20:
            cy = gg[c_in].mean()
            r["nc"] = float(np.clip(np.log2((gg[n_in].mean() + eps) / (cy + eps)), -6, 6))
            r["nc_raw_ratio"] = float(gg[n_in].mean() / cy) if cy > 0 else np.nan     # no eps / clip
            r["cyto_mean"] = float(cy); r["nuc_mean"] = float(gg[n_in].mean())
        else:
            r["nc"] = r["nc_raw_ratio"] = np.nan
        t = th[ss]; v_ = t[cell]; med = np.median(v_); mad = 1.4826 * np.median(np.abs(v_ - med)) + 1e-9
        sp, k = ndi.label((t > med + 4 * mad) & cell)
        if k:
            a = np.bincount(sp.ravel()); keep = a >= 3; keep[0] = False
            r["punct"] = float(np.clip(t[keep[sp]].sum(), 0, None) / max(np.clip(gg[cell], 0, None).sum(), 1e-9))
        else:
            r["punct"] = 0.0
        dom = cell | (own & (outd <= 3))
        bb = B_[ss] * dom
        nucleo = n_in & (ndi.gaussian_filter(B_[ss], 1.0) < 0.6 * np.median(B_[ss][n_in])) if n_in.sum() >= 20 else np.zeros_like(cell)
        ring = cell & ~ndi.binary_erosion(cell)
        y = gg[dom].astype(np.float64)
        common = {"nucleus": bb, "nucleolus": ndi.gaussian_filter(nucleo.astype(np.float32), psf) * dom,
                  "PM": ndi.gaussian_filter(ring.astype(np.float32), psf) * dom}
        fits = [("r", dict(common, cyto=ndi.gaussian_filter((cell & ~nucm).astype(np.float32), psf) * dom), RBASES)]
        if full:
            fits.append(("u", dict(common, cytosol=R_[ss] * dom, ER=Y_[ss] * dom), BASES))
        for pre, pats, BS in fits:
            names = [n for n in BS if pats[n].sum() > 0]
            A = np.stack([pats[n][dom] / pats[n][dom].sum() for n in names], 1)
            if y.sum() > 0 and len(names) >= 3:
                sc_ = y.sum(); coef, res = nnls(A.astype(np.float64), y / sc_)
                tt = coef.sum()
                for n in BS: r[f"{pre}_{n}"] = float(coef[names.index(n)] / tt) if (n in names and tt > 0) else 0.0
                yy = y / sc_; r[f"{pre}_r2"] = float(1 - res ** 2 / max(((yy - yy.mean()) ** 2).sum(), 1e-30))
        out.append(r)
    return out


# ================================================================= (1) OpenCell
def select():
    os.makedirs(OCD, exist_ok=True)
    O = pd.read_csv(f"{REF}/opencell/opencell-localization-annotations.csv")
    sp = lambda v: set(v.split(";")) if isinstance(v, str) else set()
    O["g3"] = O.annotations_grade_3.map(sp)
    O["all"] = [a | sp(b) | sp(c) for a, b, c in zip(O.g3, O.annotations_grade_2, O.annotations_grade_1)]
    def stratum(g3, al):
        if "cytoplasmic" in g3 and not (al & NUC_T): return "cyto_only"
        if (g3 & ORG_T) and not (al & NUC_T) and "cytoplasmic" not in al: return "nonnuclear_organelle"
        if (g3 & {"nucleoplasm", "chromatin"}) and al <= NUC_T: return "nuclear_only"
        if "cytoplasmic" in g3 and "nucleoplasm" in g3: return "mixed"
        return None
    O["stratum"] = [stratum(a, b) for a, b in zip(O.g3, O["all"])]
    H = pd.read_csv(f"{EVAL}/hpa_genes.tsv", sep="\t", index_col=0)
    O["hpa_jurkat_overlap"] = O.target_name.isin(H.index)
    rng = np.random.default_rng(SEED); pick = []
    for st, n in N_STRATA.items():
        c = O[O.stratum == st].sort_values("target_name")
        forced = c[c.hpa_jurkat_overlap]
        rest = c[~c.hpa_jurkat_overlap]
        k = max(0, min(n - len(forced), len(rest)))
        pick.append(pd.concat([forced, rest.iloc[np.sort(rng.choice(len(rest), k, replace=False))]]))
    forced_other = O[O.hpa_jurkat_overlap & O.stratum.isna()].assign(stratum="other_hpa_overlap")
    S = pd.concat(pick + [forced_other])
    S["annot_all"] = S["all"].map(lambda s: ";".join(sorted(s)))
    S = S[["target_name", "ensg_id", "stratum", "annotations_grade_3", "annotations_grade_2", "annotations_grade_1", "annot_all", "hpa_jurkat_overlap"]]
    S.to_csv(f"{OCD}/gene_list.tsv", sep="\t", index=False)
    log("selected", len(S), S.stratum.value_counts().to_dict(), "overlap", S.hpa_jurkat_overlap.sum())


def _get(sess, url, rng_=None, tries=4):
    for k in range(tries):
        try:
            r = sess.get(url, headers={"Range": f"bytes={rng_[0]}-{rng_[1]}"} if rng_ else {}, timeout=60)
            if r.status_code in (200, 206): return r.content
        except Exception:
            pass
        time.sleep(2 * (k + 1))
    raise RuntimeError(url)


def _list(sess, prefix):
    import re
    x = _get(sess, f"{BUCKET}/?list-type=2&max-keys=1000&prefix={prefix}").decode()
    return re.findall(r"<Key>([^<]+)</Key>", x)


def _fetch_gene(row):
    import requests, tifffile
    from scipy import ndimage as ndi
    gene, ensg = row
    sess = requests.Session(); out = []
    try:
        keys = [k for k in _list(sess, f"microscopy/raw/{gene}_{ensg}/") if k.endswith("_stack.tif")]
    except Exception as e:
        return [dict(gene=gene, status=f"list_failed:{e}")]
    if not keys: return [dict(gene=gene, status="no_fov")]
    keys = sorted(keys, key=lambda k: hashlib.md5(k.encode()).hexdigest())[:FOV_PER_GENE]   # deterministic pseudo-random
    for key in keys:
        fid = key.split("_")[-2]; fn = f"{OCD}/images/{gene}_{fid}.npz"
        if os.path.exists(fn):
            z = np.load(fn); out.append(dict(gene=gene, fov=fid, key=key, z_best=int(z["z_best"]), nz=int(z["nz"]), status="ok")); continue
        try:
            url = f"{BUCKET}/{key}"
            t = tifffile.TiffFile(io.BytesIO(_get(sess, url, (0, 65535))))
            p = t.pages[0]; md = t.imagej_metadata; nz, nc_ = int(md.get("slices", 1)), int(md.get("channels", 2))
            Hh, Ww = p.shape; nb = Hh * Ww * 2; off = p.dataoffsets[0]; dt_ = np.dtype((">" if t.byteorder == ">" else "<") + "u2")
            plane = lambda z, c: np.frombuffer(_get(sess, url, (off + (z * nc_ + c) * nb, off + (z * nc_ + c + 1) * nb - 1)), dt_).reshape(Hh, Ww).astype(np.uint16)
            focus = lambda a: float(ndi.gaussian_laplace(a.astype(np.float32), 2.0).std())   # in-focus nuclear texture
            sc = {z: focus(plane(z, 0)) for z in range(0, nz, 3)}
            z0 = max(sc, key=sc.get)
            for z in range(max(z0 - 2, 0), min(z0 + 3, nz)):
                if z not in sc: sc[z] = focus(plane(z, 0))
            zb = max(sc, key=sc.get)
            hoe, gfp = plane(zb, 0), plane(zb, 1)
            pj = tifffile.imread(io.BytesIO(_get(sess, url.replace("_stack.tif", "_proj.tif"))))
            np.savez_compressed(fn, hoechst=hoe, gfp=gfp, proj=pj.astype(np.uint16), z_best=zb, nz=nz,
                                focus_z=np.array(sorted(sc)), focus_v=np.array([sc[k] for k in sorted(sc)]))
            out.append(dict(gene=gene, fov=fid, key=key, z_best=zb, nz=nz, status="ok"))
        except Exception as e:
            out.append(dict(gene=gene, fov=fid, key=key, status=f"failed:{e}"))
    return out


def fetch():
    from concurrent.futures import ThreadPoolExecutor
    os.makedirs(f"{OCD}/images", exist_ok=True)
    S = pd.read_csv(f"{OCD}/gene_list.tsv", sep="\t")
    rows = []
    with ThreadPoolExecutor(6) as ex:            # I/O-bound threads; numpy work is tiny
        for i, r in enumerate(ex.map(_fetch_gene, zip(S.target_name, S.ensg_id))):
            rows += r
            if i % 20 == 0: log("fetched", i, "/", len(S))
    F = pd.DataFrame(rows); F.to_csv(f"{OCD}/fovs.tsv", sep="\t", index=False)
    log(F.status.str.split(":").str[0].value_counts().to_dict())


def _oc_fov(args):
    gene, fov = args
    z = np.load(f"{OCD}/images/{gene}_{fov}.npz")
    res = []
    for mode, hoe, gfp in (("mid", z["hoechst"], z["gfp"]), ("proj", z["proj"][0], z["proj"][1])):
        for c in seg_quant(hoe, gfp, gfp, None, PX_OC, full=False):
            res.append(dict(gene=gene, fov=fov, plane=mode, **c))
    return res


def quant_oc():
    import multiprocessing as mp
    F = pd.read_csv(f"{OCD}/fovs.tsv", sep="\t"); F = F[F.status == "ok"]
    with mp.get_context("fork").Pool(NPROC) as pool:
        rows = [r for rr in pool.imap_unordered(_oc_fov, list(zip(F.gene, F.fov.astype(str))), chunksize=2) for r in rr]
    C = pd.DataFrame(rows); C.to_csv(f"{OCD}/oc_cells.tsv", sep="\t", index=False)
    log("OpenCell cells", C.groupby("plane").size().to_dict())


# ================================================================= HPA side, same code
def _hpa_field(base):
    from PIL import Image
    ch = {c: np.asarray(Image.open(f"{IMG}/{base}_{c}.jpg").convert("RGB")).astype(np.float32) for c in ("blue", "red", "green", "yellow")}
    y = ch["yellow"]; y = 0.5 * (y[..., 0] + y[..., 1])
    return base, seg_quant(ch["blue"][..., 2], ch["red"][..., 0], ch["green"][..., 1], y, PX_HPA, full=True)


def quant_hpa():
    import multiprocessing as mp
    os.makedirs(OUT, exist_ok=True)
    M = pd.read_csv(f"{IMG}/manifest.tsv", sep="\t")
    G_ = M[(M.channel == "green") & M.present].copy()
    G_["base"] = G_.file.str.replace("_green.jpg", "", regex=False)
    have_y = set(M[(M.channel == "yellow") & M.present].file.str.replace("_yellow.jpg", "", regex=False))
    G_ = G_[G_.base.isin(have_y)]
    with mp.get_context("fork").Pool(NPROC) as pool:
        res = dict(pool.imap_unordered(_hpa_field, list(G_.base), chunksize=2))
    rows = []
    for r in G_.itertuples():
        for k, c in enumerate(res.get(r.base, [])):
            rows.append(dict(gene=r.gene, antibody=r.antibody, field=r.base, ab_reliability=r.ab_reliability, cell=k, **c))
    C = pd.DataFrame(rows); C.to_csv(f"{OUT}/hpa_cells_reduced.tsv", sep="\t", index=False)
    T = pd.read_csv(f"{EVAL}/hpa_cells.tsv", sep="\t")
    J = C.merge(T, on=["field", "cell"], suffixes=("", "_ref"))
    chk = dict(n_cells=len(C), n_ref=len(T), matched=len(J), few_bg_fields=int(C.groupby("field").few_bg.first().sum()))
    for m in ["nuc_frac", "nc", "punct"] + [f"u_{b}" for b in BASES]:
        chk[f"max_abs_diff_{m}"] = float((J[m] - J[f"{m}_ref" if f"{m}_ref" in J else m]).abs().max())
    json.dump(chk, open(f"{OUT}/hpa_reproduction_check.json", "w"), indent=1)
    log("HPA reproduction", chk)


# ================================================================= estimate
def _gene_agg(C, keys):
    g = C.groupby(keys)
    G = pd.DataFrame(dict(n_cells=g.size()))
    for m in ("nuc_frac", "nc", "punct", "r_r2", "u_r2", "cell_area_um2", "nuc_area_um2"):
        if m in C: G[m] = g[m].median()
    for pre, BS in (("r", RBASES), ("u", BASES)):
        cols = [f"{pre}_{b}" for b in BS]
        if all(c in C for c in cols):
            U = g[cols].mean(); U = U.div(U.sum(1), axis=0); G = G.join(U)
            G[f"{pre}_nucshare"] = G[f"{pre}_nucleus"] + G[f"{pre}_nucleolus"]
    return G


def _logistic(x, a, k, lo, hi):
    return lo + (hi - lo) / (1 + np.exp(-(a + k * np.asarray(x, float))))


def _fit_logistic(x, y):
    from scipy.optimize import curve_fit
    try:
        p, _ = curve_fit(_logistic, x, y, p0=[0.5, 1.0, 0.05, 0.85], bounds=([-10, 0, 0, 0.3], [10, 10, 0.5, 1.0]), maxfev=20000)
    except Exception:
        p = np.array([0.5, 1.0, 0.05, 0.85])
    return p


def _boot_median(x, n=2000, seed=SEED):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) == 0: return [np.nan, np.nan]
    rng = np.random.default_rng(seed)
    m = np.median(rng.choice(x, (n, len(x))), axis=1)
    return [float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))]


def estimate():
    from scipy.stats import spearmanr
    from sklearn.isotonic import IsotonicRegression
    os.makedirs(OUT, exist_ok=True)
    res = dict(task="HPA-Jurkat antibody nuclear background b (observed = b + (1-b)*true)", date=time.strftime("%Y-%m-%d"),
               script="wholecell/scripts/hpa_nuclear_background.py", seg_source=SRC.replace(V + "/", ""), seg_source_sha256=sha(SRC))
    # ---------------- HPA per gene (full = original 5-basis reference values; reduced = same code, no tubulin/ER)
    Hc = pd.read_csv(f"{OUT}/hpa_cells_reduced.tsv", sep="\t")
    Hall = _gene_agg(Hc, ["gene"])                                       # all antibodies (incl. uncertain)
    Hg = _gene_agg(Hc[Hc.ab_reliability != "uncertain"], ["gene"])        # as imaging_eval/hpa_genes.tsv
    Href = pd.read_csv(f"{EVAL}/hpa_genes.tsv", sep="\t", index_col=0)
    Hg["s_full_ref"] = (Href.u_nucleus + Href.u_nucleolus).reindex(Hg.index)
    F = pd.read_csv(f"{EVAL}/model_fractions.tsv", sep="\t"); F = F[F.version == "v7"].set_index("gene")
    Hg["model_v7_nuc"] = (F.nucleus + F.nucleolus).reindex(Hg.index)
    V0 = pd.read_csv(f"{D}/inventory_v0.tsv", sep="\t", usecols=["gene", "source"]).drop_duplicates("gene").set_index("gene").source
    Hg["ms_detected"] = Hg.index.isin(V0[V0 == "geiger"].index)
    R = pd.read_csv(f"{D}/rna_e61.tsv", sep="\t", index_col=0)
    for c in ("tpm_median", "HPA_E6.1_tpm", "call"): Hg[c] = R[c].reindex(Hg.index)
    M = pd.read_csv(f"{IMG}/manifest.tsv", sep="\t").drop_duplicates("gene").set_index("gene").locations
    Hg["hpa_locations"] = M.reindex(Hg.index)
    nucw = ("nucleoplasm", "nuclear", "nucleoli", "nucleus")
    Hg["hpa_annot_nuclear"] = Hg.hpa_locations.fillna("").map(lambda s: any(w in s for w in nucw))
    cyto = Hg[Hg.model_v7_nuc == 0]; nucl = Hg[Hg.model_v7_nuc >= 0.9]
    rf = lambda x: round(float(x), 4)
    hpa = dict(
        n_genes=len(Hg), n_cells=int(Hg.n_cells.sum()),
        floor_full_cyto_model0=dict(n=len(cyto), median=rf(cyto.s_full_ref.median()), ci95=_boot_median(cyto.s_full_ref)),
        floor_full_cyto_model0_ms=dict(n=int(cyto.ms_detected.sum()), median=rf(cyto[cyto.ms_detected].s_full_ref.median())),
        floor_full_cyto_hpa_annot_nonnuclear=dict(n=int((~Hg.hpa_annot_nuclear).sum()), median=rf(Hg[~Hg.hpa_annot_nuclear].s_full_ref.median())),
        nuclear_full_model_ge0p9=dict(n=len(nucl), median=rf(nucl.s_full_ref.median())),
        floor_reduced_cyto_model0=dict(median=rf(cyto.r_nucshare.median()), ci95=_boot_median(cyto.r_nucshare)),
        nuclear_reduced=rf(nucl.r_nucshare.median()),
        nc_log2_cyto_model0=dict(median=rf(cyto.nc.median()), ci95=_boot_median(cyto.nc)),
        nc_log2_nuclear=rf(nucl.nc.median()),
        nuc_frac_mask_cyto_model0=rf(cyto.nuc_frac.median()),
        rho_full_vs_reduced_nucshare=rf(spearmanr(Hg.s_full_ref, Hg.r_nucshare, nan_policy="omit")[0]),
        median_nuc_over_cell_area=rf((Hc.nuc_area_um2 / Hc.cell_area_um2).median()))
    res["hpa_jurkat"] = hpa
    # ---------------- (1) OpenCell
    S = pd.read_csv(f"{OCD}/gene_list.tsv", sep="\t").set_index("target_name")
    Oc = pd.read_csv(f"{OCD}/oc_cells.tsv", sep="\t")
    Fv = pd.read_csv(f"{OCD}/fovs.tsv", sep="\t")
    Og = {}
    for pl in ("mid", "proj"):
        G = _gene_agg(Oc[Oc.plane == pl], ["gene"]); G["stratum"] = S.stratum.reindex(G.index); Og[pl] = G
    # map the OpenCell nuclear/cytoplasmic concentration ratio into HPA-Jurkat geometry: isotonic fit of the HPA full
    # nuclear share (reference unmixing) on the HPA per-gene nc over all imaged genes, applied to OpenCell nc.
    hh = Hg.dropna(subset=["nc", "s_full_ref"])
    iso = IsotonicRegression(out_of_bounds="clip").fit(hh.nc, hh.s_full_ref)
    lgp = _fit_logistic(hh.nc.values, hh.s_full_ref.values)              # smooth monotone map (primary)
    # analytic alternative: per Jurkat cell, s = r*An / (r*An + Ac) in the cell's own section areas (geometry only)
    an = Hc.nuc_area_um2.values; ac = (Hc.cell_area_um2 - Hc.nuc_area_um2).clip(lower=1e-3).values
    geo = lambda r: float(np.nanmedian(r * an / (r * an + ac)))
    for pl, G in Og.items():
        G["s_jurkat_iso"] = iso.predict(G.nc.fillna(G.nc.median()))
        G["s_jurkat"] = _logistic(G.nc.values, *lgp)
        G["s_jurkat_geo"] = [geo(2 ** x) if np.isfinite(x) else np.nan for x in G.nc]
    om = Og["mid"]
    oc = dict(n_genes_selected=len(S), n_fov_ok=int((Fv.status == "ok").sum()), n_fov_failed=int((Fv.status != "ok").sum()),
              n_genes_quantified=int(om.index.nunique()), n_cells_mid=int((Oc.plane == "mid").sum()),
              pixel_nm=PX_OC, few_bg_fraction=rf(Oc.few_bg.mean()))
    for pl, G in Og.items():
        d = {}
        for st in list(N_STRATA) + ["other_hpa_overlap"]:
            x = G[G.stratum == st]
            if not len(x): continue
            d[st] = dict(n_genes=len(x), n_cells=int(x.n_cells.sum()),
                         reduced_nucshare_native=rf(x.r_nucshare.median()), reduced_nucshare_native_ci95=_boot_median(x.r_nucshare),
                         nc_log2=rf(x.nc.median()), nc_log2_ci95=_boot_median(x.nc), nuc_frac_mask=rf(x.nuc_frac.median()),
                         s_jurkat=rf(np.nanmedian(x.s_jurkat)), s_jurkat_ci95=_boot_median(x.s_jurkat),
                         s_jurkat_iso=rf(x.s_jurkat_iso.median()), s_jurkat_iso_ci95=_boot_median(x.s_jurkat_iso),
                         s_jurkat_geo=rf(x.s_jurkat_geo.median()))
        oc[f"plane_{pl}"] = d
    res["opencell"] = oc
    cyO = om[om.stratum.isin(["cyto_only", "nonnuclear_organelle"])]
    cyO1 = om[om.stratum == "cyto_only"]
    # ---------------- (2) non-expressed HPA genes (no RNA = 0 gene exists in the set; use near-zero thresholds)
    Ha = Hall.copy(); Ha = Ha.join(Hg[["s_full_ref", "tpm_median", "HPA_E6.1_tpm", "ms_detected", "model_v7_nuc", "hpa_locations"]], how="left")
    Ha["tpm_median"] = R.tpm_median.reindex(Ha.index); Ha["ms_detected"] = Ha.index.isin(V0[V0 == "geiger"].index)
    # full-unmixing share per gene from the reference per-cell table, ALL antibodies (non-expressed: reliability moot)
    T = pd.read_csv(f"{EVAL}/hpa_cells.tsv", sep="\t")
    Tg = T.groupby("gene")[["u_nucleus", "u_nucleolus", "u_cytosol", "u_ER", "u_PM"]].mean()
    Ha["s_full_allab"] = (Tg.u_nucleus + Tg.u_nucleolus) / Tg.sum(1)
    ne = {}
    for name, sel in (("rna_lt0p5", (Ha.tpm_median < 0.5) & ~Ha.ms_detected), ("rna_lt1", (Ha.tpm_median < 1) & ~Ha.ms_detected),
                      ("rna_lt1_model_nonnuclear", (Ha.tpm_median < 1) & ~Ha.ms_detected & (Ha.model_v7_nuc == 0))):
        x = Ha[sel]
        ne[name] = dict(n=int(len(x)), genes=sorted(x.index), median_s_full=rf(x.s_full_allab.median()),
                        ci95=_boot_median(x.s_full_allab), median_s_reduced=rf(x.r_nucshare.median()), median_nc=rf(x.nc.median()))
    ex = Ha[(Ha.tpm_median >= 10) & Ha.ms_detected & (Ha.model_v7_nuc == 0)]
    ne["contrast_expressed_ms_model_nonnuclear_tpm_ge10"] = dict(n=len(ex), median_s_full=rf(ex.s_full_allab.median()))
    ne["rho_tpm_vs_s_full_model_nonnuclear"] = rf(spearmanr(Ha[Ha.model_v7_nuc == 0].tpm_median, Ha[Ha.model_v7_nuc == 0].s_full_allab, nan_policy="omit")[0])
    res["hpa_nonexpressed"] = ne
    # ---------------- (3) multi-antibody genes
    Tab = T.groupby(["gene", "antibody"]).agg(n=("cell", "size"), rel=("ab_reliability", "first"),
                                               un=("u_nucleus", "mean"), uo=("u_nucleolus", "mean"))
    Tab["s_full"] = Tab.un + Tab.uo
    multi = Tab[Tab.index.get_level_values(0).isin([g for g, k in Tab.groupby(level=0).size().items() if k >= 2])]
    mr = []
    for g, x in multi.groupby(level=0):
        mr.append(dict(gene=g, n_ab=len(x), s_min=rf(x.s_full.min()), s_max=rf(x.s_full.max()), delta=rf(x.s_full.max() - x.s_full.min()),
                       antibodies=";".join(f"{a}:{r_}:{s:.3f}" for (_, a), r_, s in zip(x.index, x.rel, x.s_full)),
                       model_v7_nuc=rf(Hg.model_v7_nuc.get(g, np.nan))))
    MR = pd.DataFrame(mr)
    res["hpa_multi_antibody"] = dict(n_genes=len(MR), median_abs_delta=rf(MR.delta.median()) if len(MR) else None,
                                     min_s_over_model_nonnuclear=rf(MR[MR.model_v7_nuc == 0].s_min.median()) if len(MR) else None,
                                     rows=mr, note="only 5 genes; a between-antibody spread of this size bounds the antibody-specific part of the floor")
    # ---------------- decision
    t1 = float(np.nanmedian(cyO.s_jurkat)); t1_ci = _boot_median(cyO.s_jurkat)
    t1_native = float(cyO.r_nucshare.median())
    floor = float(cyto.s_full_ref.median())
    t2 = ne["rna_lt1"]["median_s_full"]
    # bootstrap: resample HPA genes (floor AND the nc -> share map) and OpenCell non-nuclear genes
    rng = np.random.default_rng(SEED); bs = []; t1s = []
    a_ = cyto.s_full_ref.dropna().values; onc = cyO.nc.dropna().values; hx, hy = hh.nc.values, hh.s_full_ref.values
    for _ in range(2000):
        fl = np.median(rng.choice(a_, len(a_)))
        ih = rng.integers(0, len(hx), len(hx)); p_ = _fit_logistic(hx[ih], hy[ih])
        fr = float(np.median(_logistic(rng.choice(onc, len(onc)), *p_)))
        bs.append((fl - fr) / (1 - fr)); t1s.append(fr)
    t1_ci = [float(np.percentile(t1s, 2.5)), float(np.percentile(t1s, 97.5))]
    b_fit = (floor - t1) / (1 - t1)
    if t1 <= 0.10 and abs(t2 - floor) <= 0.05:
        branch = "background: b_HPA = median nuclear fraction of non-expressed genes; not applied to OpenCell"
        b_est = t2; b_ci = ne["rna_lt1"]["ci95"]
    elif t1 >= 0.20:
        branch = "real pool: b = 0 (large mass move, Pedro reviews)"
        b_est = 0.0; b_ci = [0.0, 0.0]
    else:
        branch = "intermediate: b fitted with CI; report both extremes (b=0 and b=HPA floor)"
        b_est = b_fit; b_ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    if t1 <= 0.10 and abs(t2 - floor) > 0.05:
        branch = "background by test (1) but test (2) disagrees with the floor by > 0.05: b fitted from (1) with CI, (2) reported"
        b_est = b_fit; b_ci = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
    res["decision"] = dict(
        test1_opencell_nonnuclear_median_nuclear_fraction_jurkat_geometry=rf(t1), test1_ci95=t1_ci,
        test1_opencell_cyto_only_jurkat_geometry=rf(np.nanmedian(cyO1.s_jurkat)),
        test1_isotonic_map=rf(cyO.s_jurkat_iso.median()), test1_analytic_area_map_reduced_scale=rf(cyO.s_jurkat_geo.median()),
        test1_projection_plane_sensitivity=rf(np.nanmedian(Og["proj"][Og["proj"].stratum.isin(["cyto_only", "nonnuclear_organelle"])].s_jurkat)),
        nc_log2_opencell_mid_nonnuclear=rf(cyO.nc.median()), nc_log2_hpa_floor_genes=rf(cyto.nc.median()),
        conc_ratio_nuc_over_cyto_opencell=rf(2 ** cyO.nc.median()), conc_ratio_nuc_over_cyto_hpa=rf(2 ** cyto.nc.median()),
        test2_variants=dict(rna_lt0p5=ne["rna_lt0p5"]["median_s_full"], rna_lt1_model_nonnuclear=ne["rna_lt1_model_nonnuclear"]["median_s_full"]),
        test1_native_hek_reduced_unmixing=rf(t1_native),
        hpa_floor_model_nonnuclear=rf(floor), test2_nonexpressed_median=t2, test2_minus_floor=rf(t2 - floor),
        b_fit_from_test1=rf(b_fit), b_fit_ci95=[rf(x) for x in np.percentile(bs, [2.5, 97.5])],
        branch=branch, b_HPA=rf(b_est), b_HPA_ci95=[rf(x) for x in b_ci], b_OpenCell=0.0,
        sensitivity_extremes=dict(b_zero=0.0, b_full_floor=rf(floor)),
        formula="observed_nuclear_fraction = b + (1 - b) * true  ->  true = (observed - b) / (1 - b), clipped to [0, 1]")
    # ---------------- tables
    Hg.round(4).to_csv(f"{OUT}/hpa_genes_background.tsv", sep="\t")
    Ha.round(4).to_csv(f"{OUT}/hpa_genes_allantibodies.tsv", sep="\t")
    MR.to_csv(f"{OUT}/hpa_multi_antibody.tsv", sep="\t", index=False)
    pd.concat({k: v for k, v in Og.items()}, names=["plane", "gene"]).round(4).to_csv(f"{OUT}/opencell_genes.tsv", sep="\t")
    pd.concat({k: v for k, v in Og.items()}, names=["plane", "gene"]).round(4).to_csv(f"{OCD}/oc_genes.tsv", sep="\t")
    res["outputs"] = sorted(os.listdir(OUT)) + [f"../opencell_background/{f}" for f in ("gene_list.tsv", "fovs.tsv", "oc_cells.tsv", "oc_genes.tsv", "images/")]
    json.dump(res, open(f"{OUT}/estimate.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    log(json.dumps(res["decision"], indent=1))
    return res


if __name__ == "__main__":
    steps = sys.argv[1:] or ["all"]
    order = ["select", "fetch", "quant_oc", "quant_hpa", "estimate"]
    for s in (order if "all" in steps else steps):
        t0 = time.time(); globals()[s](); log(s, f"{time.time() - t0:.0f}s")
