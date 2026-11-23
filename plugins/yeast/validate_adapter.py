"""Yeast reference adapter of the common validation suite (wholecell/validate/suite.py, plan section 4, A8 phase 1).

The yeast build still validates with its own code (vendored engine: validation/*.py and stage 4 s4_validate, run by
the engine steps s4_validate and validate_all). This adapter READS those numbers from <build>/engine/work/ and maps
them onto the suite's check ids. Rows it returns carry `source: adapter (yeast code)`, and their legacy value is the
same number (read verbatim). Phase 2 (A8): for checks with one common definition it supplies INPUTS instead
(loc_inputs: model shares + reference sets via validate_loc_inputs.py; image_nc_table; copies_reference; core_inputs
for the core placement / crowding / densities steps) and legacy numbers of other definitions (legacy(check)) shown
side by side with the reason.

Files read (all under <build>/engine/work/): results/wholecell3d_validation.json (stage 4), validation/latest/
localization_test.json, validation/anchors.json, validation/placement_checks.json, validation/ribosome_cryoet.json,
validation/literature_checks.json, validation/chromatin_stats.json, validation/hic_test.json,
validation/nucleosome_test.json, validation/MANIFEST.json (build_stamp.py); <build>/logs/engine_geometry_package.log.
Known port defect: validation/vacuole_surface.json is never written by the port (C17 row n/a (port defect)).
"""
from __future__ import annotations

import datetime
import json
import os
import re

from wholecell.registry import register
from wholecell.validate.suite import SRC_YEAST, na, row

VERBATIM = "same file (the adapter reads the legacy number verbatim)"
C1_REFS = {"SWAT_native": "SWAT N-GFP", "SWAT_NOP1": "SWAT NOP1 promoter", "hyperLOPIT": "hyperLOPIT",
           "CYCLoPs (input: fit, not a test)": "CYCLoPs (fit input: report only)"}
C1_PLAN = {("SWAT_native", "top1"): 0.665, ("hyperLOPIT", "top1"): 0.848}     # plan s4 table (X5)


class YeastX5SuiteAdapter:
    def __init__(self, ctx):
        self.ctx = ctx
        self.W = os.path.join(ctx.bdir, "engine", "work")
        self._c = {}

    def j(self, rel):
        p = os.path.join(self.W, rel)
        if p not in self._c:
            self._c[p] = json.load(open(p)) if os.path.exists(p) else None
        return self._c[p]

    def core_inputs(self):
        """inputs of the core validation steps the port build does not hold in the core layout (A8 phase 2)"""
        out = {"lib_dir": os.path.join(self.ctx.bdir, "engine", "lib"),
               # collision shapes of the yeast packer (wc3d.PACK_PER = 15 residues per sphere): C7(a) is measured on them
               "pack_trees": os.path.join(self.ctx.bdir, "engine", "work", "yc3d", "pack_trees_15.npz")}
        pj = os.path.join(self.ctx.bdir, "geometry", "package.json")
        if os.path.exists(pj):
            out["geom"] = json.load(open(pj))["package_dir"]
        return out

    def get(self, check):
        return getattr(self, check.lower(), None)

    def _r(self, check, metric, value, rel, **kw):
        kw.setdefault("legacy", value)
        kw.setdefault("legacy_source", VERBATIM)
        return row(check, metric, value, source=SRC_YEAST, file=f"engine/work/{rel}", **kw)

    def _missing(self, check, metric, rel):
        return na(check, metric, f"n/a (no data): engine/work/{rel} not written by this build", source=SRC_YEAST)

    # --- C1 inputs for the COMMON metric (A8 phase 2): per-ORF model class shares + reference class sets (test split)
    def loc_inputs(self):
        """{comps, sets: {name: {orf: [classes]}}, model: {intended|placed: {orf: (p, mass)}}, legacy: {(set, level): row},
        report_only: {set names}} from plugins/yeast/validate_loc_inputs.py (the build's own yeast code, run in a
        subprocess; outputs into the suite folder)"""
        import subprocess
        import sys
        import pandas as pd
        d = os.path.join(self.ctx.out, "adapter", "loc_inputs")
        if not os.path.exists(os.path.join(d, "classes.json")):
            here = os.path.dirname(os.path.abspath(__file__))
            env = dict(os.environ, OPENBLAS_NUM_THREADS="1", PYTHONDONTWRITEBYTECODE="1")
            subprocess.run([sys.executable, os.path.join(here, "validate_loc_inputs.py"), "--bdir", self.ctx.bdir, "--out", d],
                           check=True, env=env)
        comps = json.load(open(os.path.join(d, "classes.json")))["classes"]
        refs = json.load(open(os.path.join(d, "refs_test.json")))
        M = pd.read_parquet(os.path.join(d, "model_shares.parquet"))
        mass = (M["copies"] * M["mass_da"]).to_dict() if "mass_da" in M.columns else {}
        model = {}
        for lvl, pre in (("intended", "assigned"), ("placed", "placed")):
            cols = [c for c in comps if f"{pre}_{c}" in M.columns]
            A = M[[f"{pre}_{c}" for c in cols]].to_numpy()
            model[lvl] = {o: ({c: float(v) for c, v in zip(cols, A[i])}, float(mass.get(o, 0.0)))
                          for i, o in enumerate(M.index)}
        L = self.j("validation/latest/localization_test.json") or []
        name = {"CYCLoPs (input: fit, not a test)": "CYCLoPs"}
        legacy = {(name.get(r["name"], r["name"]), "intended" if r["prefix"] == "assigned" else "placed"): r
                  for r in L if r["mode"] == "share"}
        return dict(comps=comps, sets=refs, model=model, legacy=legacy, report_only={"CYCLoPs"},
                    legacy_file="engine/work/validation/latest/localization_test.json")

    def c1_extras(self):
        """yeast-only C1 numbers kept as labelled extras (macro F1 of the multi-label class calls, yeast definition)"""
        rel = "validation/latest/localization_test.json"
        L = self.j(rel) or []
        rows = []
        for r in L:
            if r["mode"] != "share" or r["name"] not in C1_REFS:
                continue
            lvl = "intended" if r["prefix"] == "assigned" else "placed"
            rows.append(self._r("C1", f"{r['name'].split(' ')[0]}.{lvl}.yeast_macro_f1", r["macro_f1"], rel, direction="report",
                                note=f"yeast-only extra: macro F1 of class calls (share >= 0.2) over classes with >= 10 reference proteins; n {r['n']}"))
        return rows

    # --- C1 localisation (yeast localization.py, test split, 'share' mode; assigned = intended): phase-1 read-through
    def c1(self):
        rel = "validation/latest/localization_test.json"
        L = self.j(rel)
        if L is None:
            return [self._missing("C1", "localisation", rel)]
        rows = []
        for r in L:
            if r["mode"] != "share" or r["name"] not in C1_REFS:
                continue
            lvl = "intended" if r["prefix"] == "assigned" else "placed"
            fit = r["name"].startswith("CYCLoPs")
            for m in ("top1", "agree", "agree_mass_weighted", "macro_f1"):
                rf = C1_PLAN.get((r["name"], m))
                rows.append(self._r("C1", f"{r['name'].split(' ')[0]}.{lvl}.{m}", r[m], rel, direction="report" if fit else "up",
                                    tol="loc", reference=rf, reference_source="plan s4 (X5)" if rf else "",
                                    note=f"yeast definitions (localization.py), not yet the common agree/agree_tvd/tophit; n {r['n']}"))
        return rows

    def c2(self):
        return [na("C2", "localisation across versions", "n/a (no data): no archived yeast versions evaluated with one code "
                   "(phase 2: run C1 on the archived yeast builds S...X5)", source=SRC_YEAST)]

    # --- C5: legacy numbers for the core anchor table (validation/anchors.py)
    def legacy(self, check):
        if check == "C5":
            A = self.j("validation/anchors.json") or {}
            out = {}
            if "nucleoporins_per_pore_need" in A:
                out["nucleoporins"] = (A["nucleoporins_per_pore_need"]["median"], "validation/anchors.json nucleoporins_per_pore_need.median")
            if "core_histones_median" in A:
                out["core_histones"] = (A["core_histones_median"], "validation/anchors.json core_histones_median")
            if "ribosomes" in A:
                out["_notes"] = {"ribosomal_proteins": (
                    f"related legacy number (not the same definition): validation/anchors.json ribosomes.ratio = model 80S count / "
                    f"literature ribosomes = {A['ribosomes']['ratio'][0]:.3f}-{A['ribosomes']['ratio'][1]:.3f}")}

            return out
        if check == "C9":
            S = self.j("results/wholecell3d_validation.json") or {}
            m = S.get("membrane_orientation") or {}
            out = {"_note": ("yeast stage 4: 10,000 seeded membrane instances of the whole cell, angle to the normal at the "
                             "projected membrane point; suite: every membrane instance of the 300 placement tiles, normal = SDF "
                             "gradient at the centre (centres are on the membrane, < 0.1 nm)")}
            if "median_angle_deg" in m:
                out["median_tilt_deg"] = (m["median_angle_deg"], "results/wholecell3d_validation.json membrane_orientation.median_angle_deg (stage 4)")
                out["frac_tilt_lt_1deg"] = (m["frac_angle_lt_1deg"], "results/wholecell3d_validation.json membrane_orientation.frac_angle_lt_1deg (stage 4)")
            return out
        if check == "C11":
            S = self.j("results/wholecell3d_validation.json") or {}
            out = {f"density.{k}": (v, "results/wholecell3d_validation.json density_mg_per_ml (stage 4)")
                   for k, v in (S.get("density_mg_per_ml") or {}).items()}
            out["_refs"] = {"cytosol": (170.0, "plan s4 (X5)")}
            out["_density_note"] = (
                "definition change: stage 4 books every pre-pass thread too (DNA of transcribed genes, rDNA + pre-rRNA, mtDNA, "
                "2 bp or nt per bead), takes volumes from the label measure (tile_labvol; cell wall = mannan + glucan labels) "
                "and leaves NPCs, ATP synthase dimers and eisosomes out; the core step counts packed + indexed pre-pass "
                "instances (library bead mass) per region volume of tiles.parquet (cell wall = the outer 60 nm region), "
                "microtubules in the cytosol (nuclear MTs: their gTuRC only)")
            Lc = self.j("validation/literature_checks.json")
            nc = ((Lc or {}).get("organelles") or {}).get("model_density_ratio_to_cytosol", {}).get("nucleoplasm")
            if nc is not None:
                out["nuc_cyto_ratio"] = (nc, "validation/literature_checks.json organelles.model_density_ratio_to_cytosol.nucleoplasm")
                out["_ratio_note"] = ("definition change: legacy = nucleoplasm / cytosol (stage-4 densities); common = nucleus group "
                                      "(nucleoplasm + nucleolus) / cytosol, as Jurkat")
            return out
        if check == "C10":
            S = self.j("results/wholecell3d_validation.json") or {}
            P = self.j("validation/placement_checks.json") or {}
            out = {}
            if "placed_fraction_total" in S:
                out["placed_copies_fraction"] = (S["placed_fraction_total"], "results/wholecell3d_validation.json placed_fraction_total (stage 4)")
            if P.get("species"):
                out["species_below_50pct"] = (P["species"]["summary"]["species_below_50pct"],
                                              "validation/placement_checks.json (against the inventory, not the packing request)")
            if P.get("seams"):
                r = P["seams"]["face_over_interior_range"]
                out["seam_note"] = (f"(yeast placement_checks: sphere-volume fraction 0-5 nm / 15-35 nm from the face, "
                                    f"range over passes {r[0]:.3f}-{r[1]:.3f}; a different definition)")
            return out
        return {}

    # --- C7 (b) bead overlap (stage 4)
    def c7(self):
        rel = "results/wholecell3d_validation.json"
        S = self.j(rel)
        if S is None:
            return [self._missing("C7", "b.bead_clash_rate", rel)]
        c = S["clash"]
        return [
            self._r("C7", "b.bead_clash_rate", c["clash_rate"], rel, unit="clashing beads per sampled bead", direction="down", tol="ratefrac",
                    reference=2.55e-5, reference_source="plan s4 (X5)",
                    note=f"bead overlap > 0.3 nm, {c['cubes']} seeded cubes, {c['beads_sampled']} beads (per bead, not per bead pair)"),
            self._r("C7", "b.bead_clash_rate_without_actin_structures", c["clash_rate_without_actin_structures"], rel, direction="down", tol="ratefrac"),
            self._r("C7", "b.cross_tile_clashing_beads", c["cross_tile_clashing_beads"], rel, direction="down", tol="count0"),
        ]

    # --- C8 containment (stage 4 + placement_checks shells)
    def c8(self):
        rel = "results/wholecell3d_validation.json"
        S = self.j(rel)
        if S is None:
            return [self._missing("C8", "b.beads_outside_region", rel)]
        c = S["containment"]
        rows = [self._r("C8", "b.beads_outside_region_fraction", c["frac_beads_outside_region"], rel, direction="down", tol="ratefrac",
                        reference=6.9e-5, reference_source="plan s4 (X5)", note=f"{c['volume_instances_sampled']} volume instances, {c['beads']} beads"),
                self._r("C8", "b.instances_any_bead_outside_fraction", c["frac_instances_any_bead_outside"], rel, direction="down", tol="ratefrac")]
        P = self.j("validation/placement_checks.json")
        if P and P.get("shells"):
            worst = max(v["outside"] for v in P["shells"].values())
            rows.append(self._r("C8", "shell_outside_fraction_max", worst, "validation/placement_checks.json", direction="report",
                                note="membrane-face proteins further than shell + 4 nm from their membrane, worst shell"))
        return rows

    # --- C9 orientation: the common rows (core placement on the yeast tiles) carry the stage-4 numbers as legacy
    def c9(self):
        return []

    # --- C10 extra yeast numbers (report)
    def c10(self):
        rel = "results/wholecell3d_validation.json"
        S = self.j(rel)
        rows = []
        if S is not None:
            rows.append(self._r("C10", "yeast.placed_fraction_of_inventory", S["placed_fraction_of_inventory"], rel, direction="report"))
            rows.append(self._r("C10", "yeast.placed_fraction_outside_mitochondria", S["placed_fraction_outside_mitochondria"], rel, direction="report"))
        P = self.j("validation/placement_checks.json")
        if P and P.get("species"):
            rows.append(self._r("C10", "yeast.placed_share_by_volume", P["species"]["summary"]["placed_share_by_volume"],
                                "validation/placement_checks.json", direction="report"))
        if P and P.get("seams"):
            r = P["seams"]["face_over_interior_range"]
            rows.append(self._r("C10", "yeast.seam_volume_ratio_0_5nm_min", r[0], "validation/placement_checks.json", direction="report",
                                note="yeast definition: sphere-volume fraction 0-5 nm / 15-35 nm from the tile face, min over passes"))
            rows.append(self._r("C10", "yeast.seam_volume_ratio_0_5nm_max", r[1], "validation/placement_checks.json", direction="report"))
        return rows

    # --- C11 densities (stage 4) and N/C ratio (literature_checks)
    def c11(self):
        rel = "results/wholecell3d_validation.json"
        S = self.j(rel)
        if S is None:
            return [self._missing("C11", "density", rel)]
        rows = [self._r("C11", f"density.{k}", v, rel, unit="mg/ml", direction="stable", tol="density",
                        reference=170.0 if k == "cytosol" else None, reference_source="plan s4 (X5)" if k == "cytosol" else "")
                for k, v in S["density_mg_per_ml"].items()]
        L = self.j("validation/literature_checks.json")
        nc = ((L or {}).get("organelles") or {}).get("model_density_ratio_to_cytosol", {}).get("nucleoplasm")
        if nc is not None:
            rows.append(self._r("C11", "nuc_cyto_ratio", nc, "validation/literature_checks.json", direction="report", tol="ratio",
                                note="yeast definition: nucleoplasm (not the nucleus group) / cytosol"))
        rows.append(na("C11", "phi", "n/a (no data): excluded-volume fraction not computed by the yeast suite with the core "
                       "bead model (phase 2: core crowding step on the yeast tiles)", source=SRC_YEAST))
        return rows

    # --- C12 ribosomes (ribosome_cryoet.py)
    def c12(self):
        rel = "validation/ribosome_cryoet.json"
        R = self.j(rel)
        if R is None:
            return [self._missing("C12", "nn_median_model", rel)]
        return [
            self._r("C12", "nn_median_model", R["model_nn_median_nm"], rel, unit="nm", direction=f"target:{R['cryoet_nn_median_nm']}",
                    tol="nn", reference=R["cryoet_nn_median_nm"], reference_source=f"cryo-ET picks (CZII), {R['cryoet_tomograms']} tomograms"),
            self._r("C12", "nn_median_model_thinned", R["model_thinned_to_5200_per_um3_nn_median_nm"], rel, unit="nm", direction="report"),
            self._r("C12", "concentration_uM", R["model_uM"], rel, unit="uM", direction=f"target:{R['delarue2018_uM']}", tol="rel10",
                    reference=R["delarue2018_uM"], reference_source="Delarue 2018"),
            self._r("C12", "polysome_in_chains_fraction", R["polysomes"]["in_chains_fraction"], rel, direction="report"),
        ]

    # --- C13 chromatin (hic_test, chromatin_stats, nucleosome_test, stage 4)
    def c13(self):
        rows = []
        H = self.j("validation/hic_test.json")
        if H:
            rows.append(self._r("C13", "hic_ps_slope", H["model"]["ps_slope"], "validation/hic_test.json", direction="report",
                                reference=H.get("lazarstefanita2017_G1", {}).get("ps_slope"), reference_source="Lazar-Stefanita 2017 G1 Hi-C"))
            rows.append(self._r("C13", "hic_cis_fraction", H["model"]["cis_fraction"], "validation/hic_test.json", direction="report",
                                reference=H.get("lazarstefanita2017_G1", {}).get("cis_fraction"), reference_source="Lazar-Stefanita 2017 G1 Hi-C"))
        C = self.j("validation/chromatin_stats.json")
        if C:
            rows.append(self._r("C13", "ps_slope_10_200kb", C["ps_slope_10_200kb"], "validation/chromatin_stats.json", direction="report"))
            if "50" in C.get("rms_distance_nm", {}):
                rg = C.get("bystricky2004_wlc_rms_range_nm", {}).get("50")
                rows.append(self._r("C13", "wlc_rms_50kb_nm", C["rms_distance_nm"]["50"], "validation/chromatin_stats.json", unit="nm",
                                    direction="report", reference=rg, reference_source="Bystricky 2004 WLC range"))
        N = self.j("validation/nucleosome_test.json")
        if N:
            rows.append(self._r("C13", "nucleosome_count", N["count"], "validation/nucleosome_test.json", direction="report",
                                reference=N.get("count_oberbeckmann_range"), reference_source="Oberbeckmann range"))
            rows.append(self._r("C13", "nucleosome_occupancy_auroc", N["occupancy"]["auroc"], "validation/nucleosome_test.json",
                                direction="up", tol="loc"))
        if not rows:
            rows.append(self._missing("C13", "chromatin", "validation/{hic_test,chromatin_stats,nucleosome_test}.json"))
        return rows

    # --- C3 image N/C: the per-ORF table of validation/image_nc.py (model log2 N/C from placed copies, YeastRGB C-SWAT
    #     measured log2 N/C, split); the metric is computed by the suite (common definition)
    def image_nc_table(self):
        import pandas as pd
        rel = "validation/image_nc.tsv"
        p = os.path.join(self.W, rel)
        if not os.path.exists(p):
            return None
        J = self.j("validation/image_nc.json") or {}
        if str(J.get("status", "ok")).startswith("not_run"):      # port P7: absence of the raw archive is reported
            return dict(missing=f"WARNING: {J.get('warning', J['status'])} (status {J['status']}, engine/work/validation/image_nc.json)")
        d = pd.read_csv(p, sep="\t", index_col=0)
        leg = {sp: {"spearman": J[sp]["spearman"], "auroc_nuclear": J[sp]["auroc_nuclear_vs_cytosolic"]} for sp in ("dev", "test") if sp in J}
        return dict(df=d.reset_index().rename(columns={"orf": "gene"}), splits=["test", "dev"], test="test",
                    source="adapter (yeast image table: YeastRGB C-SWAT, image_nc.py)", file=f"engine/work/{rel}", legacy=leg,
                    legacy_source="engine/work/validation/image_nc.json", ref_spearman=0.81, ref_auroc=0.98,
                    ref_source="plan s4 (X5)", note="model: placed copies (nucleoplasm + nucleolus vs cytosol concentration); "
                    "image: YeastRGB C-SWAT mNeonGreen, Nup49-mCherry nuclear mask"
                    + (f"; status {J['status']} (raw archive absent: measured ratios from the work-seed cache)"
                       if J.get("status", "ok") != "ok" else ""))

    # --- C6 reference copy numbers (Das 2026 independent; Ho 2018 not independent), from the C1 input dump
    def copies_reference(self):
        import pandas as pd
        d = os.path.join(self.ctx.out, "adapter", "loc_inputs")
        if not os.path.exists(os.path.join(d, "copies_das2026.tsv")):
            self.loc_inputs()
        D = pd.read_csv(os.path.join(d, "copies_das2026.tsv"), sep="\t")
        H = pd.read_csv(os.path.join(d, "copies_ho2018.tsv"), sep="\t")
        J = self.j("validation/copy_numbers.json") or {}
        pr = J.get("pairs", {})
        leg = {"das2026": {k: pr.get("model_vs_das", {}).get(k) for k in ("spearman", "median_ratio", "within_2x")},
               "ho2018": {k: pr.get("model_vs_ho", {}).get(k) for k in ("spearman", "median_ratio", "within_2x")}}
        return dict(key="orf", sets={"das2026": pd.Series(D.das.values, index=D.orf.astype(str)),
                                     "ho2018": pd.Series(H.ho.values, index=H.orf.astype(str))},
                    independent={"das2026"}, plan={("das2026", "spearman"): 0.81}, legacy=leg,
                    legacy_source=f"engine/work/validation/copy_numbers.json ({J.get('proteins')} proteins)",
                    legacy_note=("definition change: legacy compares on the proteins in model, Das AND Ho (all > 0, n "
                                 f"{J.get('proteins')}); common: model and the one reference (> 0). Model copies are the core "
                                 "inventory/copies step = engine features copies_per_cell (identical, criterion I3)"),
                    file="adapter/loc_inputs/copies_{das2026,ho2018}.tsv (yeast copy_numbers.py tables)")

    # --- C15 proximity AUROC (proximity_test.py, gap 10 nm; contact matrix of the yeast GPU contacts step): read
    def c15(self):
        rel = "validation/proximity_test_g10.json"
        J = self.j(rel)
        if J is None:
            return [self._missing("C15", "proximity_auroc", rel)]
        rows = []
        plan = {"pca": 0.85, "prox_label": 0.81, "xlms": 0.99}
        for name, v in J["sets"].items():
            for k, lab in (("snapshot+assembly", "snapshot_assembly"), ("snapshot", "snapshot"), ("compartment", "compartment"),
                           ("abundance", "abundance")):
                x = v["median_auroc"].get(k)
                main = k == "snapshot+assembly"
                rows.append(self._r("C15", f"{name}.median_auroc.{lab}", x, rel, direction="up" if main else "report", tol="loc",
                                    reference=plan.get(name) if main else None, reference_source="plan s4 (X5)" if (main and name in plan) else "",
                                    note=f"median per-bait AUROC over {v['baits_ge5']} baits with >= 5 positives, gap {J['gap_nm']} nm; "
                                         "BioGRID 5.0.263 sets; contacts from the yeast GPU contacts step (common neighbour code: with BD v2)"))
        return rows

    # --- C16 build stamp (build_stamp.py MANIFEST.json), re-checked against file mtimes
    def c16(self):
        rel = "validation/MANIFEST.json"
        M = self.j(rel)
        if M is None:
            return [self._missing("C16", "yeast.stale_result_files", rel)]
        stale = [f for f, v in M["files"].items() if not v.get("of_this_build")]
        moved = []
        for f, v in M["files"].items():
            p = os.path.join(self.W, f)
            if not os.path.exists(p):
                moved.append(f)
                continue
            mt = datetime.datetime.fromtimestamp(os.path.getmtime(p)).isoformat(timespec="seconds")
            if mt != v.get("modified"):
                moved.append(f)
        return [self._r("C16", "yeast.stale_result_files", len(stale), rel, direction="down", tol="count0", reference=0, reference_source="must",
                        note=f"{len(M['files'])} files; build_stamp mixed={M.get('mixed')}"),
                row("C16", "yeast.result_files_changed_since_stamp_check", len(moved), source=SRC_YEAST + " + suite mtime re-check",
                    file=f"engine/work/{rel}", direction="down", tol="count0", note=f"missing or modified after the check: {moved[:5]}")]

    # --- C17 cell-specific extras: vacuole surface (port defect) + model vacuole from the geometry log
    def c17(self):
        rows = [na("C17", "yeast.vacuole_surface", "n/a (port defect): the yeast port never wrote validation/vacuole_surface.json (known "
                   "port defect); the model vacuole numbers are only in logs/engine_geometry_package.log (rows below)",
                   source=SRC_YEAST)]
        X = [("trna_nucleolus", "validation/trna_nucleolus.json", ("model_trna_within_overlap", "model_sup53_nm")),
             ("silent_chromatin", "validation/silent_chromatin.json", ("sets/sir_bound (model)/h4k16ac_hypoacetylated_fraction",
                                                                     "sets/sir_bound (model)/sir4_enriched_3x_fraction")),
             ("dnase", "validation/dnase_test.json", ("core_over_linker", "nucleosomes_protected_fraction", "model_occupancy_in_dhs")),
             ("kinase_localization", "validation/kinase_localization.json", ("positive_passed", "positive_total", "negative_passed", "negative_total")),
             ("er_ribosome_density", "validation/er_ribosome_density.json", ("per_um2",)),
             ("mrna_localization", "validation/mrna_localization.json", ("spearman_mlr", "auroc_class1_vs_3_4", "er_spearman")),
             ("condensate_negative_control", "condensates/test_cytosol.json", ("log phase/unstable_at_nucleus_scale", "log phase/onset",
                                                                              "log phase/auroc_pb", "log phase/auroc_sg_only"))]
        for name, rel, keys in X:
            J = self.j(rel)
            if J is None:
                rows.append(self._missing("C17", f"yeast.{name}", rel)); continue
            for k in keys:
                v = J
                for part in k.split("/"):
                    v = v.get(part) if isinstance(v, dict) else None
                if isinstance(v, bool):
                    v = int(v)
                rows.append(self._r("C17", f"yeast.{name}.{k.split('/')[-1]}", v, rel, direction="report",
                                    note="cell-specific extra (labelled, report only; yeast code)"))
        lp = os.path.join(self.ctx.bdir, "logs", "engine_geometry_package.log")
        if os.path.exists(lp):
            for line in open(lp):
                m = re.search(r"vacuole: 1 compartment, ([\d.]+) um3, surface ([\d.]+) um2, sphericity ([\d.]+)", line)
                if m:
                    for k, v, u in (("volume_um3", m.group(1), "um3"), ("surface_um2", m.group(2), "um2"), ("sphericity", m.group(3), "")):
                        rows.append(row("C17", f"yeast.model_vacuole_{k}", float(v), unit=u, direction="report", source=SRC_YEAST + " (geometry log)",
                                        file="logs/engine_geometry_package.log", legacy=float(v), legacy_source=VERBATIM))
                    break
        return rows


register("validate_adapter", "yeast_x5")(YeastX5SuiteAdapter)
