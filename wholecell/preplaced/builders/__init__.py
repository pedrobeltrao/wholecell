"""Core preplaced builders (registered as `preplaced` components).

Module builders (one independent layer, preplaced --step module:<name>): lamina_lattice.
Serial builders (fn(ctx, params), run in order on one shared context by preplaced --step serial; wholecell/preplaced/
serial.py): npc_pores, microtubules (core.microtubules procedural|traced; core.microtubules_mtoc + microtubule_pools for yeast), npc_receptors, fg_mesh, atp_synthase_dimers,
mtdna_nucleoids, polysomes, mrna_threads, mrnp_ends. Since B4 their numerics are core code; geometry comes from the cell's kernels (ctx.K, kern.py).
"""
from . import lamina_lattice  # noqa: F401
from . import (atp_synthase_dimers, fg_mesh, microtubules, mrna_threads, mrnp_ends, mtdna_nucleoids,  # noqa: F401
               npc_pores, npc_receptors, polysomes)
