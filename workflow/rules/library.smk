# structure library: library.mode frozen (import, md5-verified per file) or build (scripted canonical library + packing
# shapes; the species table may come from a plugin step that lists `before: [library]`, or from this build's inventory
# pack table: library.species {from: pack_table}, plan key library_inputs)
LIB_MANIFEST = B("library", "library_manifest.tsv")

rule library:
    input: [FROZEN_OK] + [B(x) for x in PLAN.get("library_inputs") or []] + BEFORE("library")
    output: LIB_MANIFEST
    log: LOG("library")
    threads: min(int(PLAN.get("library_nproc") or 1), MAXW)
    shell: WC("library", "import") + " > {log} 2>&1"
