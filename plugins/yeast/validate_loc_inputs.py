"""Yeast C1 reference adapter, model + reference side (A8 phase 2): dump what the common localisation metric needs.

    python plugins/yeast/validate_loc_inputs.py --bdir <build> --out <dir>

Runs the build's own snapshot of the vendored yeast code (<build>/_code/plugins/yeast/engine, read-only; the code that
wrote the build's localization_test.json) with YC_WORK=<build>/engine/work, and calls (no numerics of its own):
  localization.model_shares()            per ORF: assigned_<class> (inventory comp_weights; a build with a localisation
                                         preset other than legacy_x5: from its core protein table), placed_<class>
                                         (types.json placed copies per region), copies, mass_da
  localization.load_swat (native, NOP1), load_hyperlopit, load_cyclops_fit   ORF -> set of classes
  splits.protein_split                   held-out 'test' half (fixed 2026-11-11)
Writes <out>/model_shares.parquet and <out>/refs_test.json ({set: {orf: [classes]}}, test split) and
<out>/classes.json (the class order, localization.CLASSES). Nothing is written into the build (imports only read).
"""
import argparse
import json
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bdir", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    code = os.path.join(a.bdir, "_code", "plugins", "yeast", "engine")
    if not os.path.isdir(os.path.join(code, "pipeline")):
        here = os.path.dirname(os.path.abspath(__file__))
        code = os.path.join(here, "engine")
    work = os.path.join(a.bdir, "engine", "work")
    os.environ.update(YC_WORK=work, YC_LIB=os.path.join(a.bdir, "engine", "lib"), YC_PDB=os.path.join(a.bdir, "engine", "pdb"),
                      YC_SEED="0", YC_NPROC="1")
    sys.path[:0] = [os.path.join(code, "pipeline"), os.path.join(code, "validation")]
    import localization as LZ  # noqa: E402
    from splits import protein_split  # noqa: E402
    os.makedirs(a.out, exist_ok=True)
    M = LZ.model_shares()
    # assigned (intended) shares: the engine's inventory.build() recomputes the X5 localisation; a build made with another
    # localisation preset (B5: common) takes them from its own core protein table instead (same class map)
    cr = os.path.join(a.bdir, "config.resolved.json")
    preset = (((json.load(open(cr)).get("inventory") or {}).get("params") or {}).get("localisation") or {}).get("preset") \
        if os.path.exists(cr) else None
    pt = os.path.join(a.bdir, "inventory", "yeast", "protein_table.parquet")
    if preset not in (None, "legacy_x5") and os.path.exists(pt):
        import pandas as pd
        T = pd.read_parquet(pt, columns=["orf", "comp_weights"])
        A = {}
        for orf, cw in zip(T.orf, T.comp_weights):
            s_ = {}
            for c, v in json.loads(cw).items():
                k = LZ.COMP_CLASS.get(c)
                if k:
                    s_[k] = s_.get(k, 0.0) + v
            z = sum(s_.values()) or 1.0
            A[orf] = {f"assigned_{k}": v / z for k, v in s_.items()}
        A = pd.DataFrame.from_dict(A, orient="index").reindex(M.index).fillna(0.0)
        for c in [c for c in M.columns if c.startswith("assigned_")]:
            M[c] = 0.0
        for c in A.columns:
            M[c] = A[c].values
        print(f"yeast loc inputs: assigned shares from {pt} (localisation preset {preset})")
    M.to_parquet(os.path.join(a.out, "model_shares.parquet"))
    refs = {"SWAT_native": LZ.load_swat(), "SWAT_NOP1": LZ.load_swat("NOP1pr-GFP localization "),
            "hyperLOPIT": LZ.load_hyperlopit(), "CYCLoPs": LZ.load_cyclops_fit()}
    out = {k: {o: sorted(v) for o, v in r.items() if protein_split(o) == "test"} for k, r in refs.items()}
    json.dump(out, open(os.path.join(a.out, "refs_test.json"), "w"))
    json.dump(dict(classes=LZ.CLASSES, code=code), open(os.path.join(a.out, "classes.json"), "w"))
    print("yeast loc inputs:", len(M), "ORFs;", {k: len(v) for k, v in out.items()})
    # C6 reference tables (copy_numbers.load, reading only): Das 2026 (independent) and Ho 2018 (not independent; the
    # legacy comparison set is restricted to proteins in all three)
    # origin: yeast validation/copy_numbers.py:load @ e3b7448 (reading of the two tables, unchanged)
    import pandas as pd
    from config import DATA
    LIT = f"{DATA}/validation/literature"
    d = pd.read_excel(f"{LIT}/das2026_peroxisome/418_2026_2458_MOESM2_ESM.xlsx", sheet_name="Table S1")
    idc = [c for c in d.columns if "Majority" in str(c)][0]
    cc = [c for c in d.columns if "Mean copy" in str(c) and "Glc" in str(c)][0]
    d = pd.DataFrame(dict(orf=d[idc].astype(str).str.split(";").str[0].str.strip().str.upper(),
                          das=pd.to_numeric(d[cc], errors="coerce")))
    h = pd.read_excel(f"{LIT}/ho2018/mmc5.xlsx", header=2)
    h = pd.DataFrame(dict(orf=h["Systematic Name"].astype(str).str.upper(),
                          ho=pd.to_numeric(h["Mean molecules per cell"], errors="coerce")))
    d.to_csv(os.path.join(a.out, "copies_das2026.tsv"), sep="\t", index=False)
    h.to_csv(os.path.join(a.out, "copies_ho2018.tsv"), sep="\t", index=False)


if __name__ == "__main__":
    main()
