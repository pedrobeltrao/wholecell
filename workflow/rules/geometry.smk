# geometry: link the immutable cellgeom_v1 package named in the config into the build (geometry/package.json).
# If the geometry stage module exists (wholecell/geometry/stage.py, owner: geometry) it does this and validates;
# otherwise a minimal link is written here (package path from config geometry.package).
_GEOM_STAGE = PLAN["geometry_stage"]

rule geometry:
    input: [FROZEN_OK] + BEFORE("geometry")
    output: GEOM
    log: LOG("geometry")
    params:
        pkg = PLAN["geometry_package"]
    run:
        if _GEOM_STAGE:
            shell(WC("geometry") + " > {log} 2>&1")
        else:
            import json
            if not params.pkg or not os.path.isdir(params.pkg):
                raise ValueError(f"config geometry.package missing or not a directory: {params.pkg!r}")
            os.makedirs(os.path.dirname(output[0]), exist_ok=True)
            json.dump({"package_dir": params.pkg, "linked_by": "workflow/rules/geometry.smk"}, open(output[0], "w"), indent=1)
