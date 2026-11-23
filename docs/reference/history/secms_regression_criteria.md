# SEC-MS toolkit + assembly adapter: regression criteria vs Jurkat v7 (frozen before the first run)

Written 2026-11-16, before the toolkit (`toolkits/secms/`) or the new assembly adapter was run on Jurkat data.
Reference: the legacy v7 tables in `{data_root}/data/derived/` (identical by md5 to the port build
`jurkat/20261115-ec0b403-s20261115-v7repro2`, `inventory/compare_legacy.json`).

## Set-up

- Jurkat config with its legacy settings: `inventory.adapters.secms: secms_toolkit`, `params.secms_toolkit.preset:
  legacy_v7` (top-3 peptides, m_score < 0.01, no shared/decoy peptides, 3-fraction running mean, curated native-mass
  anchors without outlier rejection, assembled at app. MW >= k x monomer on fraction centres, k = 2/3/4/6);
  assembly options at their legacy values (`complex_copies: stoichiometry_cap`, homo-oligomer share = SEC a_k by state
  with default 0.5, AFDB homodimer share = SEC a2 with default 0.3, heterodimers SEC r >= 0.6 or STRING >= 700, bound
  0.5 of the limiting free pool). The new option `complex_copies: measured_free_pool` is OFF.
- The legacy code was deterministic for every table below (the two arbitrary orders of the legacy run - hash seed and
  unstable sort - are already handled by the port and gave md5 identity in v7repro2). Hence **exact** criteria.
- Run: inventory steps `secms`, `assembly` (then the rest of the inventory chain through `compare`) in a fresh dev build
  folder under `{results_root}/jurkat/_dev_secms/`, PYTHONHASHSEED=0, CPU only. No packing.

## Criteria

| id | what | criterion | why exact |
|---|---|---|---|
| S1 | `assembly/sec_assembly.tsv` | md5 = legacy (12785dddc48fdb917f19a0eeeb190180) | deterministic legacy script; same float expressions |
| S2 | `assembly/sec_calibrants_peak.tsv` | md5 = legacy (8d2bc00fbc9d75db9dcfdeed10f0e628) | idem |
| S3 | `assembly/homo_oligomers_v0.tsv` (shares, n_oligomers) | md5 = legacy (47c8a0a62da9cf22776e23689fe72ac1) | idem (port already identical) |
| S4 | `assembly/afdb_homodimers_v0.tsv` | md5 = legacy (205a99f82885b97c0decf9a22929c068) | idem |
| S5 | `assembly/afdb_heterodimers_v0.tsv` (sec_r, n) | md5 = legacy (2beee72eb8410c71990d7d9b04ce0161) | idem |
| S6 | `assembly/complexes_v0.tsv` | md5 = legacy (0c95b67246b26e85407f5e6989d2317e) | stoichiometry cap unchanged |
| S7 | `assembly/inventory_v0_complexed.tsv` | content identical to legacy; only the known order-only `compartments` cells (2,218, legacy hash order) may differ, as in v7repro2 | as in the port |
| S8 | toolkit per-protein a_k (k = 2, 3, 4, 6) vs the orphan `sec_assembly_thresholds.tsv` | every protein of the orphan table present; a_k equal (abs diff 0, tolerance 1e-12 only for float formatting); `rule` column equal to `frac_in_complexes` (NaN -> 0); the orphan's `mass_pg` column is NOT a criterion (pre-check 2026-11-16: it is 0.8 x the 100 pg inventory, i.e. it was written from an 80 pg inventory) | the orphan has no generator; this criterion is its generator |
| S9 | `pack_table.tsv` | md5 = d1a0099caad792d0ba42cbc52b52f0a5 (= ingredients_v7.tsv) | whole chain downstream of assembly unchanged |
| S10 | `compare_legacy.json` | every table/ingredient entry has the same identical / content status as in v7repro2 | |
| S11 | new option off | Jurkat resolved config has `complex_copies: stoichiometry_cap`; a test asserts it | Pedro: not enabled for Jurkat yet |

Verdict: REPRODUCED if S1-S11 all pass; otherwise the failing rows are listed with the first differing cells, and the
legacy inline SEC code is NOT removed until fixed. Diagnostics (toolkit defaults vs legacy, measured free pools dry
run, HeLa format test) are reported in `docs/secms_toolkit_report.md` and are not criteria.
