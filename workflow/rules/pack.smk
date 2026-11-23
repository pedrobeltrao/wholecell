# pack: tiles -> [organisation] -> index -> counts (index drops returned to the pool, copy ledger) -> pack (8 checkerboard phases, resumable) -> summary
P_ = lambda *a: B("pack", *a)
PACK_SUMMARY = P_("summary.json")
# pack.organisation on (plan key pack_organisation, preset merged): STRING network + partners before the counts
PACK_ORG = [P_("organisation", "network.npz"), P_("organisation", "partners.npz")] if PLAN.get("pack_organisation") else []

if PACK_ORG:
    rule pack_organisation:
        input: [FROZEN_OK, PACK_TABLE, LIB_MANIFEST] + PRE_SERIAL + BEFORE("pack_organisation")
        output: PACK_ORG, P_("organisation", "network.json")
        log: LOG("pack_organisation")
        shell: WC("pack", "organisation") + " > {log} 2>&1"

rule pack_tiles:
    input: [FROZEN_OK, GEOM] + BEFORE("pack_tiles")
    output: P_("tiles.parquet"), P_("tiles_meta.json")
    threads: MAXW
    log: LOG("pack_tiles")
    shell: WC("pack", "tiles", "workers={threads}") + " > {log} 2>&1"

rule pack_counts:
    input: [P_("tiles.parquet"), PACK_TABLE, LIB_MANIFEST, PRE_MANIFEST, PACK_ORG, P_("index_drops.tsv")] + BEFORE("pack_counts")
    output: P_("counts.npz"), P_("libs.json")
    log: LOG("pack_counts")
    shell: WC("pack", "counts") + " > {log} 2>&1"

rule pack_index:
    input: [P_("tiles.parquet"), P_("tiles_meta.json"), PRE_MANIFEST] + BEFORE("pack_index")
    output: P_("preplaced_types.npy"), P_("preplaced_pos.npy"), P_("preplaced_quat.npy"), P_("preplaced_off.npy"), P_("preplaced_libs.json"), P_("index_drops.tsv")
    log: LOG("pack_index")
    shell: WC("pack", "index") + " > {log} 2>&1"

rule pack_pack:
    input: [P_("counts.npz"), P_("libs.json"), P_("preplaced_off.npy"), P_("preplaced_libs.json"), LIB_MANIFEST, PACK_ORG] + BEFORE("pack_pack")
    output: P_("pack_done.json")
    threads: MAXW
    log: LOG("pack_pack")
    shell: WC("pack", "pack", "workers={threads}") + " >> {log} 2>&1"

rule pack_summary:
    input: [P_("pack_done.json")] + BEFORE("pack_summary")
    output: PACK_SUMMARY
    log: LOG("pack_summary")
    shell: WC("pack", "summary") + " > {log} 2>&1"

FINAL += [PACK_SUMMARY]
