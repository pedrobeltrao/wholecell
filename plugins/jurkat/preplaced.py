# new in wholecell (B4, 2026-11-17)
"""Jurkat data of the common pre-placed builders (registered `preplaced_data` providers).

jurkat.mrna_table: the cytosolic mRNA molecules of the frozen v7 mRNA layer (`pre_mrna_v1`, jmrna.py v1: 250,000
  mRNAs of E6-1 TPM x total, DepMap RSEM effective length + 350 nt; kinds 0 polysomal, 1 free cytosolic, 2 nuclear):
  one row per molecule of the selected kinds (gene, length nt). The common polysome builder decides which of them are
  translated and how many ribosomes each carries; the nuclear ones stay a frozen bead subset (config exports).
"""
from wholecell.registry import register


@register("preplaced_data", "jurkat.mrna_table")
def mrna_table(ctx, params):
    import numpy as np
    import pandas as pd
    from wholecell.preplaced.stage import _resolve
    z = np.load(_resolve(ctx.cfg, ctx.bdir, params.get("npz", "pre_mrna_v1")))
    k = np.isin(z["kind"], params.get("kinds", [0, 1]))
    T = pd.DataFrame(dict(gene=z["gene"][k].astype(str), length=z["nt"][k].astype(np.int64)))
    ctx.log(f"mRNA table {params.get('npz', 'pre_mrna_v1')}: {len(T):,} cytosolic mRNAs of {T.gene.nunique():,} genes, "
            f"median {int(T.length.median())} nt")
    return T
