"""Inventory table schemas (spec: docs/reference/inventory_schema.md). Column lists only; no cell, compartment or species lists.

Every table is a TSV with a header row. Adapters and rules return pandas DataFrames with (at least) these columns;
`check(df, name)` raises on missing columns. Extra columns are allowed and carried along.
"""
from __future__ import annotations

AVOGADRO = 6.02214076e23

# --- pack view (legacy-compatible; = columns of the Jurkat ingredients_v7.tsv, in this order) ---
PACK_TABLE = ["ing_id", "kind", "genes", "stoich", "copies", "mass_da", "comp_weights", "transmembrane", "n_tm",
              "struct_src", "uniprot", "signal_frac", "mito_om_frac", "mass_pg"]

# --- per-gene inventory (level 1; legacy inventory_v0.tsv columns, `gene` = index label) ---
INVENTORY_V0 = ["gene", "uniprot", "mass_da", "ppb_geiger", "n_rep", "tpm", "rna_call", "source", "ppb", "note", "copies",
                "compartments", "loc_source", "transmembrane", "mass_pg"]
# adapter.proteome returns the abundance part, adapter.localisation the location part (joined on `gene`)
PROTEOME_OUT = ["gene", "uniprot", "mass_da", "source", "copies"]
LOCALISATION_OUT = ["gene", "compartments", "loc_source", "transmembrane"]

# --- canonical tables ---
SPECIES = ["species_id", "species_key", "display_name", "kind", "genes", "uniprots", "stoich", "n_chains", "copies",
           "mass_da", "mass_pg", "struct_src", "tm", "n_tm", "signal_frac", "om_frac", "assembly_source",
           "abundance_source", "loc_sources", "last_rule_id", "build_step"]
SPECIES_LOC = ["species_id", "loc_key", "weight", "rule_id", "evidence", "step"]
FIELD_PROVENANCE = ["species_id", "field", "value", "source_id", "rule_id", "step"]
CHANGES = ["step", "rule_set", "rule_id", "species_id", "genes", "field", "old", "new", "evidence", "check",
           "mass_moved_pg"]
SOURCES = ["source_id", "kind", "path", "url", "version", "date", "md5", "licence", "redistributable", "used_by",
           "description"]
PROTEOME_SOURCES = ["source_id", "weight", "unit", "normalisation", "total_protein_pg", "total_protein_source",
                    "condition", "replicates", "notes"]
LOCALISATION_SOURCES = ["source_id", "weight", "evidence_codes", "term_map", "notes"]
GENE_ABUNDANCE = ["gene", "uniprot", "source_id", "value", "unit", "n_rep", "used"]          # long: one row per gene x source
SECMS_PROTEIN = ["uniprot", "gene", "source_id", "condition", "replicate", "n_fractions", "apex_fraction",
                 "apex_app_kda", "mono_kda", "assembled_frac", "free_frac", "threshold_k", "calibration_id", "notes"]
SECMS_COMPLEX = ["complex_id", "species_id", "source_id", "condition", "replicate", "n_subunits_detected",
                 "n_subunits", "coelution_peak", "coelution_height", "mean_corr", "assembled_frac", "free_frac",
                 "app_kda", "expected_kda", "notes"]
# SEC-MS extensions (toolkits/secms; docs/reference/inventory_schema.md §4a): extra columns appended to the two tables above
SECMS_PROTEIN_EXTRA = ["assembled_frac_lo", "assembled_frac_hi", "assembled_frac_sd", "n_rep", "n_peptides", "void_frac",
                       "n_peaks", "assembled_frac_peaks", "tm", "species", "measurable",
                       "assembled_frac_mid"]
SECMS_COMPLEX_EXTRA = ["complex_source", "name", "stoich_known", "measured", "window_lo_kda", "tm_any", "calibrant",
                       "calibrant_rank", "calibrant_kept", "calibration_id", "species"]
SECMS_COMPLEX_SUBUNIT = ["complex_id", "complex_source", "uniprot", "gene", "stoich", "source_id", "condition", "replicate",
                         "detected", "mono_kda", "in_complex_frac", "free_frac", "window_lo_kda", "apex_fraction",
                         "apex_app_kda"]
SECMS_PAIRS = ["uniprot_a", "uniprot_b", "source_id", "condition", "replicate", "r", "r_sd", "n_rep", "apex_delta",
               "pair_source"]
SECMS_PEAKS = ["uniprot", "source_id", "condition", "replicate", "peak", "apex_fraction", "left", "right", "area_frac",
               "app_kda", "ratio_to_mono", "cls"]
SECMS_CALIBRATION = ["source_id", "condition", "replicate", "calibration_id", "method", "role", "kind", "id", "gene",
                     "fraction", "kda", "log10_fit", "residual", "kept", "loo_residual"]
SECMS_CALIBRANTS = ["cpx", "name", "n", "mw_kda", "peak", "co_height", "mcorr", "key", "kept", "source_id", "condition",
                    "replicate"]
SECMS_PROFILE_ID = ["uniprot", "source_id", "condition", "replicate", "kind"]          # + one column per fraction: f1..fN
# per-protein SEC-MS bound fraction vs the modelled assembly state (secms_bound.py; docs/reference/inventory_schema.md §4c)
SECMS_BOUND_FRACTION = ["gene", "uniprot", "copies", "mass_da", "mass_pg", "secms_source_id", "secms_condition",
                        "secms_replicate", "datasets", "n_datasets", "confidence", "threshold_k", "k_basis",
                        "secms_bound_fraction", "secms_bound_fraction_sd", "secms_bound_fraction_lo",
                        "secms_bound_fraction_hi", "apex_app_kda", "mono_kda", "model_state", "structure_represented",
                        "model_species", "model_homo_state", "modelled_free_copies", "modelled_assembled_fraction",
                        "unrepresented_bound_fraction", "unrepresented_share_of_free", "na_reason"]
SURFACE = ["uniprot", "gene", "source_id", "condition", "replicate", "pm_abundance", "unit", "pm_fraction",
           "topology", "n_term_side", "evidence", "notes"]
VERSION_CHANGES = ["gene", "transition", "change", "max_shift_target", "tvd_target", "shifts_target", "copies_from",
                   "copies_to", "reason_class", "rule_ids", "evidence"]


def check(df, cols, name="table"):
    miss = [c for c in cols if c not in df.columns and not (c == df.index.name)]
    if miss:
        raise ValueError(f"{name}: missing columns {miss}")
    return df


def empty(cols):
    import pandas as pd
    return pd.DataFrame({c: pd.Series(dtype=object) for c in cols})
