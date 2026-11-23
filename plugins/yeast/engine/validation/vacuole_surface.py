# origin: yeast yeast_cell_dev_bundle/validation/vacuole_surface.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Is the model's vacuole the right shape? Membrane area against volume, from three independent kinds of data (2026-11-14).

  1. Cryo soft X-ray tomography of unbudded early-G1 BY4741 in YPD (Chen et al. 2025 Mol Biol Cell 36:ar132,
     Supplementary file S8, wild type): per cell, vacuole volume and area. Area = area of the equal-volume sphere /
     corrected sphericity (SpherCorr; equals the corrected surface of FIJI 3D Suite, checked on the Cell and Nucleus
     sheets: literature/vacuole_surface/README.md)
  2. The segmented cryo-FIB/SEM cells of EMPIAR-11416 (yeast_organelle_meshes bundle, manifest): all vacuole pieces
     per cell, and single pieces above 0.3 um3
  3. Lipids: phospholipid per protein of isolated vacuoles (Zinser et al. 1991 J Bacteriol 173:2026, Table 4:
     0.51 mg/mg; Schneiter et al. 1999 J Cell Biol 146:741, Table 3: 1/0.9) times the model's vacuolar protein (membrane
     and lumen, expected copies), 760 g/mol per phospholipid (assumption) and 0.674 nm2 per lipid (DOPC at 30 C), two
     leaflets: the bilayer area the lipids of the vacuole can make
The model's own vacuole: volume and area from the build's stage-2 log.
Output: work/validation/vacuole_surface.json, docs/figures/vacuole_surface.png
"""
import os, sys, re, json, math, glob, zipfile
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import DATA, WORK, W3D, ROOT  # noqa: E402

LIT = f"{DATA}/validation/literature"
S8 = f"{LIT}/chen2025_sxt/mbc-36-ar132-s008.xlsx"
MESHES = os.path.join(os.environ.get("YC_INPUTS", os.path.join(ROOT, "..", "..", "inputs")), "yeast_organelle_meshes.zip")
PL_PER_PROTEIN = {"Zinser 1991": 0.51, "Schneiter 1999": 1 / 0.9}
PL_MW, A_LIPID = 760.0, 0.674          # g/mol (assumption), nm2 per lipid
sphere_area = lambda v: (36 * math.pi * v * v) ** (1 / 3)


def chen_g1():
    C = pd.read_excel(S8, sheet_name="Cell"); V = pd.read_excel(S8, sheet_name="Vacuole")
    strain = dict(zip(C["Cell ID"].str.rsplit("_", n=1).str[0], C["Strain"]))
    V = V[V["Cell ID"].str.rsplit("_", n=1).str[0].map(strain) == "WT"]
    px = V["pix size"].astype(str).str.split().str[0].astype(float)
    vol = V["Volume(Pix)"] * px ** 3
    area = vol.map(sphere_area) / V["SpherCorr(Pix)"]
    return vol.values, area.values


def fibsem():
    z = zipfile.ZipFile(MESHES); m = json.loads(z.read("yeast_organelle_meshes/outputs/manifest.json"))
    cells, pieces = {}, []
    for p in m["parts"]:
        if p["organelle"] == "vacuole" and p.get("provenance") == "real":
            c = cells.setdefault(p["cell_id"], [0.0, 0.0, 0]); c[0] += p["high_volume_um3"]; c[1] += p["high_area_um2"]; c[2] += 1
            if p["high_volume_um3"] > 0.3:
                pieces.append((p["high_volume_um3"], p["high_area_um2"]))
    return cells, np.array(pieces)


def model():
    log = open(f"{WORK}/s2_geometry.log").read()
    m = re.findall(r"vacuoles?: (?:\d+ segmented lobes.*?-> )?(\d+) compartments?, ([\d.]+) um3, surface ([\d.]+) um2", log)
    n, vol, area = m[-1]
    e = pd.read_parquet(f"{W3D}/expected_totals.parquet")
    ing = pd.read_csv(f"{WORK}/results/inventory/ingredients.csv").drop_duplicates("ing_id").set_index("ing_id")
    g = e[e.region.isin(["mem_VAC", "vacuole"])]
    prot_fg = float((g.expected * g.ing_id.map(ing.mass_da).fillna(0.0)).sum() / 6.022e23 * 1e15)
    lipid_area = {k: prot_fg * r * 1e-15 / PL_MW * 6.022e23 * A_LIPID * 1e-6 / 2 for k, r in PL_PER_PROTEIN.items()}
    return dict(compartments=int(n), volume_um3=float(vol), area_um2=float(area), vacuolar_protein_fg=prot_fg, lipid_area_um2=lipid_area)


def main():
    v, a = chen_g1(); cells, pieces = fibsem(); mo = model()
    r = a / np.array([sphere_area(x) for x in v])
    res = dict(chen2025_g1=dict(n=int(len(v)), volume_median=float(np.median(v)), area_median=float(np.median(a)),
                                area_over_sphere_median=float(np.median(r)), area_over_sphere_p95=float(np.percentile(r, 95)),
                                share_above_1p5=float((r > 1.5).mean())),
               fibsem_cells={c: dict(vacuoles=n, volume_um3=round(V_, 2), area_um2=round(A_, 1), area_over_sphere=round(A_ / sphere_area(V_), 2))
                             for c, (V_, A_, n) in sorted(cells.items())},
               fibsem_pieces_area_over_sphere_median=float(np.median(pieces[:, 1] / np.array([sphere_area(x) for x in pieces[:, 0]]))),
               model=mo, model_area_over_sphere=mo["area_um2"] / sphere_area(mo["volume_um3"]))
    os.makedirs(f"{WORK}/validation", exist_ok=True)
    json.dump(res, open(f"{WORK}/validation/vacuole_surface.json", "w"), indent=1)
    print(json.dumps(res, indent=1))
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8.2, 5.6))
    xs = np.linspace(0.05, 11, 200); ax.plot(xs, [sphere_area(x) for x in xs], color="#6b7280", lw=1.2, ls="--", zorder=1)
    ax.text(10.4, sphere_area(10.4) - 2.5, "one sphere", fontsize=8.5, color="#6b7280", ha="right")
    ax.scatter(v, a, s=12, color="#2a78b5", alpha=0.55, lw=0, label=f"G1 vacuoles, soft X-ray tomography (Chen 2025; n = {len(v)})", zorder=2)
    ax.scatter(pieces[:, 0], pieces[:, 1], s=22, facecolor="white", edgecolor="#7b9e3c", lw=1.3, label="single FIB-SEM vacuole pieces > 0.3 um3", zorder=3)
    for c, (V_, A_, n) in cells.items():
        ax.scatter(V_, A_, s=60, marker="s", color="#d95f02", zorder=4)
        off = {"cell_04": (6, 4), "cell_05": (-6, -13)}.get(c, (5, 3)); ha = "right" if off[0] < 0 else "left"
        ax.annotate(f"{c.replace('cell_', 'cell ')} ({n} pieces)", (V_, A_), xytext=off, textcoords="offset points", fontsize=8, color="#1f2328", ha=ha)
    ax.scatter([], [], s=60, marker="s", color="#d95f02", label="whole FIB-SEM cells, all pieces (EMPIAR-11416)")
    lo_, hi_ = min(mo["lipid_area_um2"].values()), max(mo["lipid_area_um2"].values())
    ax.plot([mo["volume_um3"]] * 2, [lo_, hi_], color="#8e5aa8", lw=6, alpha=0.5, solid_capstyle="butt", zorder=2,
            label="membrane the vacuole's lipids can make (Zinser 1991, Schneiter 1999)")
    ax.scatter(mo["volume_um3"], mo["area_um2"], s=150, marker="*", color="#1f2328", zorder=5, label=f"this build: {mo['area_um2']:.1f} um2")
    ax.set_xlim(0, 11); ax.set_ylim(0, 80); ax.grid(color="#e5e7eb", lw=0.6); ax.set_axisbelow(True)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    ax.set_xlabel("vacuole volume (um3)", color="#6b7280"); ax.set_ylabel("vacuolar membrane area (um2)", color="#6b7280")
    ax.set_title("Vacuole membrane area against volume", loc="left", fontsize=11.5, color="#1f2328")
    ax.legend(loc="upper left", fontsize=8, frameon=False)
    fig.tight_layout(); os.makedirs(os.path.join(HERE, "..", "docs", "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "..", "docs", "figures", "vacuole_surface.png"), dpi=110)


if __name__ == "__main__":
    main()
