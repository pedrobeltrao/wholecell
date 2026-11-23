# pre-placed structures: modules with a builder get one rule each (build + compare with the frozen reference +
# selection); frozen modules are read directly; serial builders (preplaced.serial: one ordered list on a shared context,
# e.g. yeast) are one rule `preplaced_serial`; the manifest lists the files the pack uses (config order).
PRE_MODS = PLAN["preplaced_modules"]
PRE_BUILT = {m["name"]: m for m in PRE_MODS if m.get("builder")}
PRE_MANIFEST = B("preplaced", "preplaced_manifest.json")


def _pre_input(v):
    if not isinstance(v, str) or v.startswith("layer:"):
        return []                                   # package layers: covered by geometry/package.json
    if v.startswith("build:"):
        return [B(v.split(":", 1)[1])]
    if v.startswith("module:"):
        n = v.split(":", 1)[1]
        if n in PRE_BUILT:
            return [B("preplaced", n, "selected.json")]
        m = next(x for x in PRE_MODS if x["name"] == n)
        return _pre_input(m["npz"])
    if v.startswith("frozen:"):
        return [FROZEN(v.split(":", 1)[1])]
    if v in PLAN["frozen"]:
        return [FROZEN(v)]
    return [v]


for _n, _m in PRE_BUILT.items():
    rule:
        name: f"preplaced_{SAFE(_n)}"
        input: [FROZEN_OK, GEOM, LIB_MANIFEST] + [p for v in (_m.get("inputs") or {}).values() for p in _pre_input(v)] + BEFORE(f"preplaced_{SAFE(_n)}")
        output: B("preplaced", _n, "selected.json"), B("preplaced", _n, "consumption.json")
        log: LOG(f"preplaced_{SAFE(_n)}")
        params: step=f"module:{_n}"
        shell: WC("preplaced", "{params.step}") + " > {log} 2>&1"

PRE_SERIAL = []
if PLAN.get("preplaced_serial"):
    # serial builders on one context (wholecell/preplaced/serial.py); the context reads the inventory (pack table) and
    # whatever plugin steps list `before: [preplaced_serial]` (yeast: the engine view of this build's tiles)
    _SO = PLAN["preplaced_serial"]["out"]
    PRE_SERIAL = [B(_SO, f) for f in ("preplaced.npz", "expected_after_preplaced.parquet", "lib_ids.npy", "serial_manifest.json")]
    rule preplaced_serial:
        input: [FROZEN_OK, GEOM, LIB_MANIFEST, PACK_TABLE] + BEFORE("preplaced_serial")
        output: PRE_SERIAL
        threads: min(8, MAXW)
        log: LOG("preplaced_serial")
        shell: WC("preplaced", "serial") + " > {log} 2>&1"

rule preplaced_manifest:
    input: [FROZEN_OK, GEOM] + [B("preplaced", n, "selected.json") for n in PRE_BUILT] + PRE_SERIAL + BEFORE("preplaced_manifest")
    output: PRE_MANIFEST
    log: LOG("preplaced_manifest")
    shell: WC("preplaced", "manifest") + " > {log} 2>&1"
