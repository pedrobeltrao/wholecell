# origin: jurkat/inventory/rna_e61.py (legacy workspace, 2026-11-15)
# changes: wrapped as adapter.rna "jurkat.rna_e61"; file paths come from frozen inputs named in config
# (inventory.params."jurkat.rna_e61".sources, ordered = output column order) instead of relative paths + globs (the glob
# matches of 2026-11-15 are listed explicitly, in sorted order); Hart zTPM / call thresholds are config values.
# Numerics identical to legacy.
"""E6-1 RNA consensus: per-gene TPM from several bulk sources, renormalised to protein-coding genes, plus Hart-style
zTPM per source and an expression call. Returns a DataFrame indexed by HGNC symbol (index name `gene`)."""
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from wholecell.registry import register


def _files(ctx, spec):
    fr = spec["frozen"]
    return [ctx.frozen(f) for f in (fr if isinstance(fr, list) else [fr])]


def rna_e61(ctx):
    p = ctx.params("jurkat.rna_e61")
    srcs = p["sources"]
    hspec = srcs[p["ens2sym_from"]]
    h = pd.read_csv(_files(ctx, hspec)[0], sep="\t")
    ens2sym = dict(zip(h.Gene, h["Gene name"]))
    pc = set(pd.read_csv(ctx.frozen(p["protein_coding_frozen"]), sep="\t", usecols=["Gene Names (primary)"]).iloc[:, 0].dropna())
    rna = {}

    def rsem(f):
        x = pd.read_csv(f, sep="\t"); x["sym"] = x.gene_id.str.split(".").str[0].map(ens2sym); return x.groupby("sym").TPM.sum()

    def fc(f):
        x = pd.read_csv(f, sep="\t", comment="#"); c = x.columns[-1]; rpk = x[c] / x.Length
        return pd.Series((rpk / rpk.sum() * 1e6).values, index=x.Geneid.map(ens2sym)).groupby(level=0).sum()

    for name, spec in srcs.items():
        fs = _files(ctx, spec); rd = spec["reader"]
        if rd == "hpa_celline":
            rna[name] = h[h["Cell line"] == spec["cell_line"]].groupby("Gene name").TPM.max()
        elif rd == "expression_atlas":
            e = pd.read_csv(fs[0], sep="\t")
            rna[name] = e.groupby("Gene Name")[e.columns[-1]].max()
        elif rd == "rsem_mean":
            rna[name] = pd.concat([rsem(f) for f in fs], axis=1).mean(axis=1)
        elif rd == "featurecounts_mean":
            rna[name] = pd.concat([fc(f) for f in fs], axis=1).mean(axis=1)
        elif rd == "deseq2_tpm_wt_ct":
            x = pd.read_csv(fs[0], sep="\t")
            symcol = [c for c in x.columns if c.lower() in ("gene_name", "genename", "symbol", "gene_symbol", "gene")] or \
                     [c for c in x.columns if "name" in c.lower() or "symbol" in c.lower()][:1]
            ct = [c for c in x.columns if c.endswith("_tpm") and "WT" in c and "_CT" in c and c[-5].isdigit()]
            rna[name] = x.groupby(symcol[0])[ct].sum().mean(axis=1)
        else:
            raise ValueError(f"unknown RNA reader {rd}")
    R = pd.DataFrame(rna)
    R = R.loc[R.index.isin(pc)].fillna(0)
    R = R / R.sum() * 1e6                                  # protein-coding TPM
    bw = p.get("kde_bw", 0.1)

    def hart(s):                                           # Hart 2013 zTPM: KDE mode, right-half sd, on log2 TPM>0
        v = np.log2(s[s > 0].values); g = np.linspace(-8, 16, 2000); k = gaussian_kde(v, bw)(g)
        mu = g[k.argmax()]; sd = np.sqrt(np.mean((v[v > mu] - mu) ** 2))
        return (np.log2(s.clip(lower=1e-6)) - mu) / sd
    Z = R.apply(hart)
    out = R.add_suffix("_tpm").join(Z.add_suffix("_z"))
    out["tpm_median"] = R.median(axis=1); out["z_median"] = Z.median(axis=1)
    out["n_expr"] = (Z >= p["z_expressed"]).sum(axis=1); out["n_off"] = (Z < p["z_off"]).sum(axis=1)
    out["call"] = np.where(out.n_expr >= p["min_sources_expressed"], "expressed", np.where(out.n_off == R.shape[1], "off", "low"))
    out.index.name = "gene"
    ctx.say(str(R.shape), out.call.value_counts().to_dict())
    return out


rna_e61.PARAM_FILES = []
register("adapter.rna", "jurkat.rna_e61")(rna_e61)
