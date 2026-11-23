# origin: yeast yeast_cell_dev_bundle/scripts/vacuole_template.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""One real segmented vacuole as the shape of the model's single G1 vacuole -> pipeline/vacuole_template.npz (kept in the
repository).

Why one lobe (2026-11-14; IMPROVEMENTS round W): in unbudded G1 cells the vacuole is one compact body. Cryo soft X-ray
tomography of 289 wild-type G1 cells (Chen et al. 2025 Mol Biol Cell 36:ar132, Supplementary file S8) gives a median
vacuole area 1.03-1.13 times that of a sphere of the same volume; the vacuoles of all five segmented FIB-SEM cells of
EMPIAR-11416 are fragmented (11-38 pieces, 2.3-3.4 times), and those cells are 44-74 um3, larger than G1 cells (about
25 um3), as expected after vacuoles fragment at bud emergence. Single segmented lobes are as compact as G1 vacuoles
(1.03-1.23). The template is the largest lobe that is fully captured (z extent close to its in-plane extent: the
slice-wise segmentation loses the caps of flatter ones) and compact: TEMPLATE below.
Usage: YC_ORGANELLE_MESHES=.../yeast_organelle_meshes python scripts/vacuole_template.py
"""
from _common import *
from organelle_stats import load_obj, moments
SRC = os.environ.get("YC_ORGANELLE_MESHES", "")
TEMPLATE = "cell_06/vacuole_01"     # 1.41 um3, area 1.06 x the equal-volume sphere, z / in-plane extent 1.09
OUT = f"{C.ROOT}/vacuole_template.npz"


def main():
    if not SRC:
        raise SystemExit("set YC_ORGANELLE_MESHES to the yeast_organelle_meshes folder")
    man = json.load(open(f"{SRC}/outputs/manifest.json"))
    p = next(p for p in man["parts"] if p["id"] == TEMPLATE)
    V, F = load_obj(f"{SRC}/outputs/{p['obj']['high']}")
    V = V + np.array(p["obj_offset_um"]["high"])
    vol, cen, ax = moments(V, F)
    np.savez_compressed(OUT, V=(V - cen).astype(np.float32), F=F.astype(np.int32),
                        source=np.array(f"EMPIAR-11416 {TEMPLATE}, yeast_organelle_meshes bundle (high LOD), um, centred"))
    log(f"{TEMPLATE}: {vol:.2f} um3, area {p['high_area_um2']:.2f} um2, semi-axes {np.round(ax, 2)} -> {OUT}")


if __name__ == "__main__":
    main()
