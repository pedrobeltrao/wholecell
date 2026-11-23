# new in wholecell (B6b, 2026-11-17): no legacy origin
"""Rounded crista rims for a Jurkat package (package N -> N+1; B6b, criteria docs/unification_criteria/B6b_crista_rims.md).

The lamellar cristae (section cristae_lamellae_v1, toolkits/em/jcristae.py) end CR_GAP inside the boundary inner
membrane in a flat cap with square corners. In cells the crista rim is a tightly curved ridge that carries the ATP
synthase dimer rows (Davies et al. 2011 PNAS 108:14121; Strauss et al. 2008 EMBO J 27:1154). The evaluator rounds the
cap when the engine parameter CR_RIM > 0 (crista_lumen: rounded intersection of the slab and the cap plane, radius
r = min(CR_RIM, CR_W); apex unchanged, rim = half-cylinder of radius r). This step sets CR_RIM and changes nothing else:
fields, sections, pores, meshes and layers are carried byte-identical.
CR_RIM = CR_W = 10 nm: half the 20 nm IM midsurface separation across the crista lumen (inter-cristae membrane
distance 20 +- 4 nm, Davies 2011), i.e. the largest rounding that keeps the lumen width.

usage: NUMBA_NUM_THREADS=3 python -m toolkits.em.jrims [src={geom_root}/jurkat/em_cell1_v8p] [geometry_id=em_cell1_v8r] [cr_rim=10]
"""
import copy
import hashlib
import os
import sys
import time

CODE = ("toolkits/em/jrims.py", "wholecell/geometry/engines/cellgeom_eval_v1.py")
CR_RIM = 10.0


def _sha(rel):
    from wholecell.config import REPO
    return hashlib.sha256(open(os.path.join(REPO, rel), "rb").read()).hexdigest()


def main(src="{geom_root}/jurkat/em_cell1_v8p", geometry_id="em_cell1_v8r", cr_rim=CR_RIM):
    from toolkits.em import pkgstep
    t0 = time.time()
    cr_rim = float(cr_rim)
    cell, parts = pkgstep.open_pkg(src)
    man = cell.manifest
    assert float(man["engine"]["params"].get("CR_RIM", 0.0)) == 0.0, "parent already has rounded rims"
    assert cr_rim <= float(man["engine"]["params"]["CR_W"]), "CR_RIM > CR_W would need a wider lumen"
    summ = pkgstep.summary_copy(cell)
    summ.setdefault("cristae", {})["CR_RIM"] = cr_rim
    summ["cristae"]["rims"] = ("rounded (B6b): half-cylinder of radius CR_RIM at each lamella end; the sample IM/OM "
                               "ratios above were measured on the square caps (rounding lowers IM area ~2.7 %)")
    cell.manifest = man = copy.deepcopy(man)
    man["description"] = man.get("description", "") + "; crista rims rounded (CR_RIM, B6b)"
    man_new, res = pkgstep.publish(cell, parts, geometry_id, step="jrims", summary=summ, params={"CR_RIM": cr_rim},
                                   generator_extra={"method": f"engine.params.CR_RIM = {cr_rim}; all content files carried",
                                                    "criteria": "docs/unification_criteria/B6b_crista_rims.md",
                                                    "code_sha256": {f: _sha(f) for f in CODE}})
    print(f"{man_new['package_id']}: validation pass {res['pass']}, failed {res.get('failed_checks')}, "
          f"waived {res.get('waived_failures')} ({time.time() - t0:.0f} s)", flush=True)
    return man_new, res


if __name__ == "__main__":
    main(*sys.argv[1:4])
