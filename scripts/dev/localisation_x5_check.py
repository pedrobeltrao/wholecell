#!/usr/bin/env python
"""Check the core localisation framework with preset legacy_x5 against the yeast port build s0 (unification plan §3d /
A5, yeast half): the yeast source adapters write evidence.tsv, the combiner (legacy_x5) gives per-protein comp_weights,
compared with s0's engine/work/features/protein_features.parquet (one row per UniProt entry, comp_weights = the
engine's inventory.build(), rounded to 4 decimals). Target: identical; acceptance |delta| <= 1e-6 per protein x key.

    python scripts/localisation_x5_check.py [--s0 <build dir>] [--out <results dir>]

Runs the real v0 adapter (wholecell.inventory.localisation.localisation) with config/cells/yeast.yaml (its
inventory.params.localisation block; the A5 fragment was merged in A9), writing under <out>/run_<stamp>/inventory/. The
CYCLoPs screens are .xls (BIFF): when this env has no xlrd, they are converted once to TSV (lossless float repr) by the
yeast engine env (pipeline_engine.python_env) into <out>/cyclops_tsv/ and read from there.
Results: <out>/check_report.json, <out>/comp_weights_diff.tsv (proteins with any difference).
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

os.environ.setdefault("PYTHONHASHSEED", "0")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

S0 = ("{results_root}/yeast/20261116-aabe7b5-s0-x5repro")
OUT = "{results_root}/yeast/_dev_unify/a5"
CYCLOPS = ["loc_yeast_cyclops_wt1", "loc_yeast_cyclops_wt2", "loc_yeast_cyclops_wt3"]


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def convert(src, dst, header):
    """(run in an env with xlrd) Excel -> TSV, the table below the header row; floats written as repr (round trip)"""
    import pandas as pd
    pd.read_excel(src, header=header).to_csv(dst, sep="\t", index=False)


def main():
    import pandas as pd
    import yaml
    from wholecell import config as C
    from wholecell.inventory.context import StepContext
    from wholecell.inventory import localisation as L
    from wholecell.inventory.localisation import evidence as EV

    ap = argparse.ArgumentParser()
    ap.add_argument("--s0", default=None)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    t0 = time.time()
    site = C.load_site()
    s0 = a.s0 or S0.replace("{results_root}", site["results_root"])
    out = a.out or OUT.replace("{results_root}", site["results_root"])
    os.makedirs(out, exist_ok=True)
    cfg = C.load("yeast")
    raw = C.load_raw("yeast")
    frag = {"inventory": {"params": {"localisation": raw["inventory"]["params"]["localisation"]}},
            "frozen_inputs": {k: v for k, v in raw["frozen_inputs"].items() if k.startswith("loc_yeast_")}}
    fi = cfg["frozen_inputs"]
    used = ["inv_yeast_uniprot"] + [k for k in frag["frozen_inputs"]]
    md5s = {}
    for k in used:
        got = md5(fi[k]["path"])
        md5s[k] = got
        if got != fi[k]["md5"]:
            raise SystemExit(f"frozen input {k}: md5 {got} != declared {fi[k]['md5']} ({fi[k]['path']})")
    conv = {}
    try:
        import xlrd  # noqa: F401
    except ImportError:
        cdir = os.path.join(out, "cyclops_tsv")
        os.makedirs(cdir, exist_ok=True)
        py = os.path.join(site["data_root"], "conda_envs", "yeast_x5", "bin", "python")   # an env with xlrd
        for k in CYCLOPS:
            src = fi[k]["path"]
            dst = os.path.join(cdir, os.path.basename(src).replace(".xls", ".tsv"))
            if not os.path.exists(dst):
                subprocess.run([py, os.path.abspath(__file__), "--convert", src, dst, "3"], check=True)
            conv[k] = {"from": src, "from_md5": md5s[k], "tsv": dst, "tsv_md5": md5(dst), "converter": py}
        over = dict(frag)
        over["frozen_inputs"] = {k: dict(v) for k, v in frag["frozen_inputs"].items()}
        for k in CYCLOPS:
            over["frozen_inputs"][k]["path"] = conv[k]["tsv"]
        over["inventory"]["params"]["localisation"]["sources"][-1]["header"] = 0
        cfg = C.load("yeast", overrides=over)
    stamp = time.strftime("%Y%m%d-%H%M%S")
    bdir = os.path.join(out, f"run_{stamp}")
    ctx = StepContext(cfg, bdir, "v0")
    uni = cfg.frozen("inv_yeast_uniprot")
    genes = list(EV.locus_table(uni).index)
    prot = pd.DataFrame(index=pd.Index(genes, name="gene"))
    prot["orf"] = genes                                     # config genes.key: orf (evidence keyed by ORF)
    res = L.localisation(ctx, prot)
    F = pd.read_csv(os.path.join(bdir, "inventory", "localisation", "gene_features.tsv"), sep="\t", keep_default_na=False)
    acc = dict(zip(F.gene, F.uniprot))

    def parse(s):
        return {k: float(v) for k, v in (kv.rsplit(":", 1) for kv in s.split(";"))} if s else {}

    mine = {acc[g]: parse(s) for g, s in res.compartments.items()}
    ref = pd.read_parquet(os.path.join(s0, "engine", "work", "features", "protein_features.parquet"))
    rows = []
    maxd = 0.0
    nkey = 0
    for u in ref.index:
        r = json.loads(ref.at[u, "comp_weights"])
        m = mine.get(u)
        if m is None:
            rows.append(dict(uniprot=u, orf=ref.at[u, "orf"], gene=ref.at[u, "gene"], max_abs_delta=float("inf"),
                             s0=json.dumps(r), core="MISSING", s0_loc_source=ref.at[u, "loc_source"]))
            continue
        d = max(abs(r.get(k, 0.0) - m.get(k, 0.0)) for k in set(r) | set(m))
        ks = set(r) != set(m)
        nkey += ks
        if d > 0 or ks:
            rows.append(dict(uniprot=u, orf=ref.at[u, "orf"], gene=ref.at[u, "gene"], max_abs_delta=d,
                             s0=json.dumps(r), core=json.dumps(m), s0_loc_source=ref.at[u, "loc_source"]))
        maxd = max(maxd, d)
    D = pd.DataFrame(rows, columns=["uniprot", "orf", "gene", "max_abs_delta", "s0", "core", "s0_loc_source"])
    D = D.sort_values("max_abs_delta", ascending=False)
    D.to_csv(os.path.join(out, "comp_weights_diff.tsv"), sep="\t", index=False)
    n_diff = int((D.max_abs_delta > 0).sum())
    n_tol = int((D.max_abs_delta > 1e-6).sum())
    rep = dict(
        check="localisation legacy_x5 vs yeast port build s0 comp_weights (plan §3d, A5 yeast half)",
        s0=s0, reference="engine/work/features/protein_features.parquet (comp_weights)", run=bdir,
        n_s0=int(len(ref)), n_core=len(mine), n_compared=int(len(ref)), max_abs_delta=maxd,
        n_proteins_differing=n_diff, n_proteins_over_1e6=n_tol, n_key_set_differences=int(nkey),
        identical=(n_diff == 0 and nkey == 0 and len(D) == 0), pass_1e6=(n_tol == 0 and maxd <= 1e-6),
        worst=D.head(20).to_dict("records"), frozen_md5=md5s, cyclops_conversion=conv,
        python=sys.executable, wall_s=round(time.time() - t0, 1))
    json.dump(rep, open(os.path.join(out, "check_report.json"), "w"), indent=1, default=str)
    print(f"legacy_x5 vs s0: {len(ref)} proteins, max |delta| {maxd:.3g}, differing {n_diff}, > 1e-6 {n_tol}, "
          f"key-set differences {nkey}; identical={rep['identical']} pass={rep['pass_1e6']} ({rep['wall_s']} s)")
    if len(D):
        print(D.head(10).to_string(index=False))
    return 0 if rep["pass_1e6"] else 1


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--convert":
        convert(sys.argv[2], sys.argv[3], int(sys.argv[4]))
        sys.exit(0)
    sys.exit(main())
