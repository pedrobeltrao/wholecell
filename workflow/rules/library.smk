# structure library: library.mode frozen (import, md5-verified per file) or build (scripted canonical library + packing
# shapes; the species table may come from a plugin step that lists `before: [library]`)
LIB_MANIFEST = B("library", "library_manifest.tsv")

rule library:
    input: [FROZEN_OK] + BEFORE("library")
    output: LIB_MANIFEST
    log: LOG("library")
    threads: min(int(PLAN.get("library_nproc") or 1), MAXW)
    shell: WC("library", "import") + " > {log} 2>&1"
