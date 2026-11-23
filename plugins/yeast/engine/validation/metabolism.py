# origin: yeast yeast_cell_dev_bundle/validation/metabolism.py @ e3b7448 (no licence yet; local only, see wholecell/NOTICE.md)
"""Metabolism against the structural model (2026-11-16): held-out data, nothing here is tuned.

Data: pipeline/data/validation/literature/metabolism/ (README there; every number read in its source).
1. Glycolytic capacity. For each glycolytic step, the flux the model cell must carry (Park et al. 2016 Nat Chem Biol
   12:482, net fluxes in mM/min per cell water, minimal glucose batch) divided by the model's copies of the step's
   isoenzymes gives the turnover each enzyme molecule would need (k_app). Chen & Nielsen 2021 (PNAS 118:e2108391118)
   measured k_app in glucose batch from their own fluxes and proteomics, and k_max, the highest k_app over 26
   conditions. Model k_app near the measured one: the enzyme copies match the flux. Above k_max: too few copies.
2. Metabolite census. Park 2016's absolute concentrations (88 yeast entries; Pi, mitochondrial CoA and CO2 are not
   measured there and are kept apart) as molecules, mass and volume per model cell.
3. Osmotic balance. Measured metabolites plus K+ (Masaryk & Sychrova 2022) against the medium give an ideal
   van 't Hoff turgor, to set beside the 0.2-0.6 MPa used by the endocytosis sub-project.
The model cell's water is taken as WATER_FRAC of its volume (Park's 56 fL of water in a cell with 19 pg dry mass
gives about 0.8; derived, not measured for this cell).
Output: work/validation/metabolism.json
"""
import os, sys, json, re
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "pipeline"))
from config import DATA, WORK  # noqa: E402

D = f"{DATA}/validation/literature/metabolism"
NA = 6.02214e23
CELL_VOL = 34.0           # um3 = fL (inventory.CELL_VOL)
WATER_FRAC = 0.80         # derived from Park 2016's conversion factors (56 fL water, 19 pg dry mass at 1.35 g/mL)
RT = 2.519                # MPa L / mol at 30 C (van 't Hoff: Pi = c R T)
# Park 2016 flux rows -> isoenzymes of the step (Chen & Nielsen's grouping)
STEPS = [("HXK", "glucose uptake (net)", ["HXK1", "HXK2", "GLK1", "YLR446W", "YDR516C"]),
         ("PGI", "PGI net", ["PGI1"]), ("PFK", "PFK net", ["PFK1", "PFK2"]), ("FBA", "FBA net", ["FBA1"]),
         ("TPI", "TPI net", ["TPI1"]), ("GAPDH", "GAPD net", ["TDH1", "TDH2", "TDH3"]), ("PGK", "PGK net", ["PGK1"]),
         ("GPM", "PGM net", ["GPM1", "YOR283W"]), ("ENO", "ENO net", ["ENO1", "ENO2", "YPL281C", "YMR323W", "YOR393W"]),
         ("PYK", "PYK net", ["CDC19", "PYK2"]), ("PDC", "PDC net", ["PDC1", "PDC5", "PDC6"])]
CHEN_KEY = {"HXK": "HXK2", "PGI": "PGI1", "PFK": "PFK", "FBA": "FBA1", "TPI": "TPI1", "GAPDH": "TDH3", "PGK": "PGK1",
            "GPM": "GPM1", "ENO": "ENO1", "PYK": "CDC19", "PDC": "PDC1", "ADH": "ADH1"}
# molar masses (g/mol) of the metabolites that make up most of the pool; the rest are counted at DEFAULT_MW
MW = {"glutamate": 147.1, "glutamine": 146.1, "citrulline": 175.2, "alanine": 89.1, "arginine": 174.2, "pyruvate": 88.1,
      "trehalose": 342.3, "threonine": 119.1, "aspartate": 133.1, "asparagine": 132.1, "lysine": 146.2,
      "ornithine": 132.2, "glutathione": 307.3, "serine": 105.1, "fructose-1,6-bisphosphate": 340.1,
      "glucose-6-phosphate": 260.1, "fructose-6-phosphate": 260.1, "valine": 117.1, "NAD+": 663.4, "ATP": 507.2,
      "citrate": 192.1, "proline": 115.1, "UDP-N-acetyl-glucosamine": 607.4}
DEFAULT_MW = 250.0
VBAR = 0.65               # mL/g, partial specific volume of small organic solutes (approximate; unchecked)
BASIC = ("arginine", "lysine", "ornithine", "citrulline")   # mostly vacuolar in yeast (not resolved by Park's whole-cell values)


def model_copies():
    f = pd.read_parquet(f"{WORK}/features/protein_features.parquet")
    c = {}
    for g, o, n in zip(f.gene.astype(str), f.orf.astype(str), f.copies_per_cell):
        for k in re.split(r"[;\s]+", g) + [o]:
            if k:
                c[k.upper()] = c.get(k.upper(), 0.0) + float(n)
    return c


def chen(kind):
    k = pd.read_csv(f"{D}/kinetics.tsv", sep="\t")
    k = k[k.parameter.str.startswith(kind) & k.source.str.startswith("Chen")]
    out = {}
    for genes, v in zip(k.gene.astype(str), k.value):
        gs = genes.split()
        key = next((s for s in CHEN_KEY.values() if s in gs), None) or ("PFK" if "PFK1" in gs else None)
        if key:
            out[key] = float(v)
    return out


def glycolysis(copies, water_L):
    fl = pd.read_csv(f"{D}/fluxes.tsv", sep="\t")
    park = {q: float(v) for q, v, s in zip(fl.quantity, fl.value, fl.source) if s.startswith("Park")}
    kapp, kmax = chen("kapp in vivo"), chen("kmax in vivo")
    rows = []
    for step, q, genes in STEPS:
        flux = park[q] * 1e-3 / 60 * water_L * NA            # molecules/s per model cell
        n = sum(copies.get(g, 0.0) for g in genes)
        key = CHEN_KEY[step]
        k_model = flux / n if n else float("nan")
        rows.append(dict(step=step, flux_per_s=flux, enzyme_copies=n, k_app_model=k_model,
                         k_app_measured=kapp.get(key), k_max_measured=kmax.get(key),
                         model_over_measured=(k_model / kapp[key]) if key in kapp else None,
                         capacity_over_flux=(n * kmax[key] / flux) if key in kmax else None))
    return rows


def census(water_L):
    m = pd.read_csv(f"{D}/metabolites.tsv", sep="\t")
    p = m[m.source.str.startswith("Park") & (m.unit == "M")].copy()
    p["value"] = p.value.astype(float)
    p["name"] = p.metabolite.str.replace(r"\[.\]\**|\*+", "", regex=True).str.strip()
    not_measured = p.note.astype(str).str.startswith("NOT MEASURED")
    drop = p.name.str.contains(r"assumed 1/2|hexose-phosphate", regex=True) | not_measured
    meas = p[~drop].copy()
    whole = meas[meas.compartment.str.startswith("whole")]
    meas = pd.concat([whole, meas[~meas.compartment.str.startswith("whole") & ~meas.name.isin(set(whole.name))]])
    meas = meas.drop_duplicates("name")
    meas["mw"] = meas.name.map(MW).fillna(DEFAULT_MW)
    tot_M = float(meas.value.sum())
    g_per_L = float((meas.value * meas.mw).sum())
    basic = float(meas[meas.name.isin(BASIC)].value.sum())
    top = meas.sort_values("value", ascending=False).head(12)
    return dict(n_metabolites=len(meas), total_M=tot_M, molecules_per_cell=tot_M * water_L * NA,
                mass_g_per_L_water=g_per_L, mass_pg_per_cell=g_per_L * water_L * 1e12,
                volume_fraction_of_water=g_per_L * VBAR * 1e-3, mw_known_share=float(meas[meas.name.isin(MW)].value.sum() / tot_M),
                basic_amino_acids_M=basic, amino_acids_M=float(meas[meas.name.str.match(r"^(alanine|arginine|asparagine|aspartate|glutamate|glutamine|isoleucine|lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine)")].value.sum()),
                atp_molecules=float(meas.set_index("name").value.get("ATP", 0)) * water_L * NA,
                pi_literature_input_M=float(p[p.name.str.startswith("phosphate")].value.sum()),
                largest=[(n, round(v * 1e3, 1)) for n, v in zip(top.name, top.value)])


def osmotic(cen, water_L):
    m = pd.read_csv(f"{D}/metabolites.tsv", sep="\t")
    k = m[m.metabolite.str.startswith("K+")].iloc[0]
    conv = {r.metabolite: float(r.value) for r in m.itertuples() if str(r.metabolite).startswith("[conversion]") and str(r.source).startswith("Park")}
    L_per_gDW = conv["[conversion] aqueous volume per cell"] / conv["[conversion] dry weight per cell"]
    k_M = float(k.value) * 1e-9 / 1e-3 / L_per_gDW                  # nmol/mg DW -> mol/L water
    # YNB without amino acids + 2% glucose (Difco formulation, salts in g/L: (NH4)2SO4 5, KH2PO4 1, MgSO4 0.5, NaCl 0.1,
    # CaCl2 0.1; unchecked), counted as fully dissociated ions
    medium = 0.111 + 3 * 5 / 132.1 + 2 * 1 / 136.1 + 2 * 0.5 / 120.4 + 2 * 0.1 / 58.4 + 3 * 0.1 / 111.0
    inside = dict(metabolites=cen["total_M"], potassium=k_M, phosphate_input=cen["pi_literature_input_M"])
    tot = sum(inside.values())
    return dict(inside_osmolar=inside, inside_total=tot, medium_osmolar=medium,
                turgor_ideal_MPa=(tot - medium) * RT,
                turgor_without_Pi_MPa=(tot - inside["phosphate_input"] - medium) * RT,
                note="ideal solutions; K+ measured in 100 mM KCl medium; Pi is a literature input in Park 2016; "
                     "the vacuole's own osmolytes (polyphosphate, basic amino acids) are not separated")


def main():
    water_L = CELL_VOL * WATER_FRAC * 1e-15
    copies = model_copies()
    gl = glycolysis(copies, water_L)
    cen = census(water_L)
    osm = osmotic(cen, water_L)
    prot_g_per_L = 3.5 / (CELL_VOL * WATER_FRAC) * 1e3              # pg / fL = g/mL; x 1000 -> g/L (inventory TOTAL_PROTEIN_PG)
    res = dict(water_fL=water_L * 1e15, glycolysis=gl, census=cen, osmotic=osm,
               protein_g_per_L_water=prot_g_per_L)
    print(f"{'step':6s} {'flux/s':>10s} {'copies':>10s} {'k_model':>8s} {'k_meas':>8s} {'k_max':>8s} {'model/meas':>10s} {'cap/flux':>8s}")
    for r in gl:
        f = lambda v, p=1: f"{v:8.{p}f}" if isinstance(v, float) and np.isfinite(v) else f"{'-':>8s}"
        print(f"{r['step']:6s} {r['flux_per_s']:10.3g} {r['enzyme_copies']:10.3g} {f(r['k_app_model'])} {f(r['k_app_measured'])} "
              f"{f(r['k_max_measured'])} {f(r['model_over_measured'], 2):>10s} {f(r['capacity_over_flux'], 2)}")
    print(f"metabolites: {cen['n_metabolites']} measured, {cen['total_M'] * 1e3:.0f} mM, {cen['molecules_per_cell']:.3g} molecules, "
          f"{cen['mass_pg_per_cell']:.2f} pg ({cen['mass_g_per_L_water']:.0f} g/L water; protein {prot_g_per_L:.0f} g/L), "
          f"{cen['volume_fraction_of_water'] * 100:.1f}% of the water by volume; amino acids {cen['amino_acids_M'] * 1e3:.0f} mM "
          f"(basic {cen['basic_amino_acids_M'] * 1e3:.0f} mM); ATP {cen['atp_molecules']:.3g} molecules")
    print(f"osmotic: inside {osm['inside_total']:.3f} Osm {json.dumps({k: round(v, 3) for k, v in osm['inside_osmolar'].items()})}, "
          f"medium {osm['medium_osmolar']:.3f} -> ideal turgor {osm['turgor_ideal_MPa']:.2f} MPa ({osm['turgor_without_Pi_MPa']:.2f} without Pi)")
    os.makedirs(f"{WORK}/validation", exist_ok=True)
    json.dump(res, open(f"{WORK}/validation/metabolism.json", "w"), indent=1, default=float)
    return res


if __name__ == "__main__":
    main()
