# origin: yeast pipeline/wc3d.py:prepass (NPC part, lines "li_npc = ..." to "npc_q = np.array(npc_q)") @ e3b7448
# Changed: own registered serial builder (was the first half of wc3d.prepass); library id and stream number are params;
# B4: geometry primitives through ctx.K (no engine module), the pore count checked against params.n_pores, the
# placed NPCs recorded as an instance group (module export, consumption).
# npc-overlap fix (2026-11-17): params.drop_frozen (core context): frozen pre-placed instances overlapping an NPC are dropped.
"""Nuclear pore complexes at the geometry's pore positions (serial pre-placed builder `core.npc_pores`).

One NPC per pore of the cell geometry (cell.npc_pos / cell.pores normals; the pores come from the geometry package,
`pores.npz`, EM detections and random fill), library +z (cytoplasmic face) along the pore normal, random spin about
it (stream `seed`). params: lib_id (scaffold structure of the organism: yeast 7N9F `NPC`, Jurkat 7R5K `NPC`),
seed, n_pores (optional: the configured pore count; a package with another count is an error), ingredient (pool
consumed, default = lib_id), drop_frozen (optional, core context: obstacle sets of frozen pre-placed layers, e.g.
[centrioles, actin]), drop_tol_nm (default 0). Adds pre: npc_pos, npc_q, li_npc; obstacles: obs["npc"] (sphere trees
of all NPCs).

Pre-placed obstacles: the NPC positions are the geometry's pores (J1: one NPC per pore) and cannot move to avoid an
obstacle. Rule: an instance of a frozen pre-placed layer (context_params.obstacles) that overlaps an NPC (library
sphere trees, overlap > drop_tol_nm) is dropped from that layer (logged; left out of the obstacle sets the later
builders see and of the frozen_subset exports, so its copies return to the pack). Frozen layers built for another
geometry version (the v7 centrioles / cortical actin did not know the 1,696 pores added in em_cell1_v8p/v8r) can
then not overlap the pores' NPCs. Builders that run before core.npc_pores and place movable instances must list
"npc" among their own obstacles instead (none does: the NPCs are placed first).
"""
import numpy as np

from wholecell.registry import register


@register("preplaced", "core.npc_pores")
def build(ctx, params):
    p3 = ctx.K.p3
    lp, cell = ctx.lp, ctx.cell
    li_npc = int(np.where(lp["ids"] == params.get("lib_id", "NPC"))[0][0])
    npc_pos = cell.npc_pos
    if params.get("n_pores") is not None and len(npc_pos) != int(params["n_pores"]):
        raise ValueError(f"core.npc_pores: the geometry has {len(npc_pos)} pores, config n_pores = {params['n_pores']}")
    q = np.empty(4); npc_q = []
    rng = np.random.default_rng(ctx.seed(params.get("seed", 5)))
    for u in cell.pores[:, 3:6]:          # NE normals
        p3.quat_z_to(u.astype(np.float64), rng.uniform(0, 2 * np.pi), q); npc_q.append(q.copy())
    npc_q = np.array(npc_q)
    types = np.full(len(npc_pos), li_npc, np.int32)
    if params.get("drop_frozen"):          # core context: frozen layers give way to the NPCs at the geometry's pores
        from wholecell.preplaced.context import drop_frozen_overlaps
        ctx.log(f"NPCs vs frozen obstacles {list(params['drop_frozen'])}:")
        drop_frozen_overlaps(ctx, types, npc_pos, npc_q.astype(np.float32), set(params["drop_frozen"]), float(params.get("drop_tol_nm", 0.0)),
                             why="overlap with an NPC at a geometry pore (core.npc_pores drop_frozen)")
    ctx.pre.update(npc_pos=npc_pos.astype(np.float32), npc_q=npc_q.astype(np.float32), li_npc=li_npc)
    ctx.obs["npc"] = ctx.emit(types, npc_pos.astype(np.float32), npc_q.astype(np.float32))
    if hasattr(ctx, "groups"):            # core context: module export + consumption
        from wholecell.preplaced.context import add_group
        add_group(ctx, "npc", types, npc_pos, npc_q, pore=np.arange(len(npc_pos), dtype=np.int32),
                  consumed={params.get("ingredient", params.get("lib_id", "NPC")): len(npc_pos)})
        ctx.log(f"NPCs: {len(npc_pos)} at the geometry's pores")
