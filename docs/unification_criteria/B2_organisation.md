# B2 organisation layer on (STRING; AF3 pool off): settings and acceptance criteria (frozen before any B2 run)

Step B2 of the unification plan (`docs/reference/history/unification_plan.md` §2b P6, §4 C15, §5, §5b B2): the
organisation layer ported in A7 phase 2 (`wholecell/pack/organisation.py`; switch `pack.organisation`, defaults in
`wholecell/pack/switches.py`) is switched on for both cells with the published-build setting: STRING only, AlphaFold 3
pool scores OFF. Order of steps: B6 → B3 → B4 → B1 → **B2** → B5; B2 is applied on the B1 state. Written 2026-11-18
before any B2 code, cheap check or build existed. The only numbers used are those of existing builds and checks (named
with their paths). The md5 of this file goes to PLANNING_NOTES (coordinator). Later findings go to a dated addendum
file, never into this one.

Pedro's decision that binds this step: **AF3 pool scores OFF in published builds** (the pool is in-house and not
deposited); the organisation layer uses STRING only. No replicate-seed builds are run; seed noise is judged with the §4
tolerances (and, for yeast, the existing X5 s0/s1/s2 spread).

Plan §2b P6: "organisation layer (STRING PPR field + partner CSR) | `pack.organisation: {string_taxon, pool_scores:
optional}` | on; AF3 pool off in published builds unless deposited | ... Output change: local neighbourhoods only, not
compartment totals." Plan §5b B2: "compartment totals unchanged (per (species, region) within max(3√n, 1 %)); C15
proximity AUROC ≥ reference; no regression elsewhere". Made concrete below.

## 1. What the layer does (A7 phase 2 code; unchanged by B2)

Every protein of the STRING network (combined score ≥ min_score, edge weight score / 1000) gets a personalized
PageRank profile (restart 0.5, own entry dropped); a packed type takes the mean profile of its subunits; R = Pearson
correlation of the type profiles. **Large scale (counts step):** R ≈ E Eᵀ from its 128 leading eigenvectors; 128 smooth
random fields (correlation length 500 nm, seed stream `pack.organisation.field`) give each type a log weight β E_i·F(x)
(β = 0.8) per tile; copies of a (type, region) are split over tiles ∝ region measure × exp(log weight) with a crowding
correction. The per (type, region) totals are the same split totals as without the layer. **Small scale (pack step):**
loose partners = the ≤ 16 network neighbours with r ≥ 0.4 (strength 0.2 + 0.8 (r − 0.4) / 0.6); contact partners only
from AF3 pool pairs (**none in B2**). The kernel proposes a copy next to an already placed partner (gap 2–15 nm for
loose partners) with probability min(0.6, Σ_p 0.25 × strength × n_p / (8 + n_p)), 40 tries; every proposal passes the
same clash / containment tests as a random one.

## 2. B2 settings per cell

| Key (`pack.organisation`) | Jurkat B2 | Yeast B2 |
|---|---|---|
| enabled | **true** (B1 state: false) | true (unchanged) |
| string.taxon / version | **9606** / "12.0" | 4932 / "12.0" (unchanged) |
| string.links | frozen input `inv_string_links` (`data/interactions/string/9606.protein.links.v12.0.txt.gz`, md5 5b7c91d02926e723ad30e81ab9d009a4, already declared) | `{yeast_root}/yeast_cell_dev_bundle/pipeline/data/string/4932.protein.links.v12.0.txt.gz` (unchanged) |
| string.min_score | **400** (combined score; = yeast X5 `organisation.STRING_MIN`, = the A7 phase 2 testcell run) | 400 (unchanged) |
| string.ids | `inv_uniprot_full` (20,415 reviewed human entries), column `STRING` (full STRING id "9606.ENSP…") | `inv_yeast_uniprot`, column "Gene Names (ordered locus)" (unchanged) |
| genes | `inv_uniprot_full`, column "Gene Names (primary)" (complex subunits without a UniProt field) | as now |
| subunits | the inventory pack table (`uniprot`, else `genes` mapped); types = the pack-table ingredients with a library entry (default of `organisation.network_ids`) | as now (assembly ingredients, lib_map, packed library order) |
| **pool_scores / paralogs** | **null / null** | **null / null** (B1 state: the cached AF3 pool table `wholecell_frozen/yeast_x5/work_seed/organisation/pool_pairs.parquet` + Alliance paralogs) |
| other constants | `switches.DEFAULTS` (= the yeast X5 constants: damping 0.5, profile nodiag, field_k 128, field_l_nm 500, n_features 128, beta 0.8, r_loose 0.4, max_partners 16, a_loose 0.25, k_half 8, p_max 0.6, attract_tries 40, loose_gap_nm [2, 15], dense_max_nodes 8000 → the human network uses the Neumann-series PPR) | unchanged (the same values) |
| field seed | stream `pack.organisation.field` from the Jurkat master seed | legacy stream value 11 (unchanged) |

**Yeast: B2 is "pool off", not "organisation on".** The yeast cell config has carried the organisation layer with the
AF3 pool ON since A7 (the X5 setting, needed to reproduce X5). In B2 the yeast build drops the pool; the X5 reproduction
configs (`yeast_x5_core.yaml`, `yeast_x5_port.yaml`, preset `legacy_x5`) keep it (reproduction only, never published).
Dropping the pool changes the partner sets: A7 phase 2 measured on the s0 inputs (`scripts/organisation_a7_check.py`)
pool off vs on **62,813 vs 69,656 partner links, 74 % of the type rows change** (contact partners gone, 59 % of the
pool-on partner pairs kept). The STRING profiles, R and the field loadings E do not depend on the pool, so the
large-scale field is unchanged; only the small-scale partner attraction changes.

**Switch mechanism (no effect on B1).** The B2 settings live in two config overlays, `config/cells/jurkat_b2.yaml` and
`config/cells/yeast_b2.yaml` (`extends:` the cell config, `cell_id` unchanged: builds land under `jurkat/` and `yeast/`,
the same pattern as `jurkat_v7_port.yaml` / `yeast_x5_core.yaml`), each holding only the `pack.organisation` keys of the
table above. `config/cells/{jurkat,yeast}.yaml` and `switches.COMMON_DEFAULTS` are not touched, so no build of `jurkat`
or `yeast` (B1 included) changes, whichever commit it is built from. The B2 queue job builds
`scripts/build_isolated.sh jurkat_b2 b2 24` / `yeast_b2 b2 24` from a main commit that contains both B1 (common
defaults) and the B2 overlays. After acceptance the overlay keys are folded into the base cell configs and the overlays
removed (one commit; the effective config equals the accepted build's, checked by comparing `config.resolved.json`).

## 3. References

| Cell | Reference build (no_regression and every comparison) |
|---|---|
| Jurkat | the B1 build: the last `RUNNING build <dir>` line of `jurkat/_dev_unify/b1/full_build.status`, with its `validate/suite.json` |
| Yeast | the B1 build: the same from `yeast/_dev_unify/b1/full_build.status` (organisation on, AF3 pool ON) |

No fallback: B2 is defined against the B1 state. Without a complete B1 build (snakemake "100%", `validate/suite.json`)
the B2 verdict is not computed (`NO REFERENCE`). The B1 build's own verdict (ACCEPTED or not) is recorded.

## 4. Acceptance items (verdict `<build>/validate/compare_B2/verdict.json`)

| # | Criterion | Measurement | Pass |
|---|---|---|---|
| O0 | Settings applied | `config.resolved.json` of the B2 build vs the B1 build: `pack` sections differ only in `pack.organisation` (Jurkat: enabled + the STRING wiring; yeast: pool_scores, paralogs), and nothing outside `pack` differs except cell name / build id / paths of the build; `pack/organisation/network.json`: `stats.pool_pairs` 0, `stats.contact_links` 0, `string.taxon` 9606 / 4932, `string.version` "12.0", `string.min_score` 400, links md5 recorded (Jurkat = the frozen md5) | all |
| O1 | Compartment totals | placed copies per (species, region) = Σ over `pack/tiles/*.npz` of `job_placed` by (`libs.json` library id of `job_lib`, region name of `job_region`), B2 vs reference; n = reference copies; pair within tolerance when \|Δ\| ≤ max(3√n, 0.01 n) (n = 0: \|Δ\| ≤ 3) | ≥ 99.5 % of the pairs within tolerance AND the out-of-tolerance pairs' summed \|Δ mass\| ≤ 0.1 % of the placed packed mass; every out-of-tolerance pair listed (worst 100 by \|Δ\| / tolerance) |
| O1r | Requested totals (report) | the same with `job_n` (the packer's request per (species, region)) | report; identical when no thinning is active; with thinning the pilot factors may move (factors from `pack/thinning.json` reported side by side) |
| O2 | C15 proximity | yeast: suite C15 rows `<set>.median_auroc.snapshot_assembly` (adapter: `engine/work/validation/proximity_test_g10.json`, per-bait AUROC of P(a copy of A has a copy of B within 10 nm) in the packed cell + same complex, median over baits with ≥ 5 positives; BioGRID 5.0.263 sets) for every set in the build (pca, prox_label, xlms, kid_*, kinase_*), B2 vs reference | each ≥ reference − 0.01 (§4 AUROC tolerance); the strict "≥ reference" of the plan is reported per set for pca / prox_label / xlms (plan values 0.85 / 0.81 / 0.99; s0 0.815 / 0.844 / 0.994). Jurkat: **n/a** (below) |
| O3 | Partner placement (report + activity gate) | `pack/summary.json` `organisation`: **placed_next_to_partner / instances** (the fraction of packed copies placed by a partner proposal next to a STRING partner) and copies_with_partner_proposals / instances, B2 and reference; types with partners / with a profile (network.json) | report; gate: Jurkat fraction > 0 (the layer acts). Yeast: report the change pool on (B1) → off (B2) |
| O4 | No regression elsewhere | suite re-run with `ref_suite=` the reference's `validate/suite.json` into `compare_B2/suite/`, `no_regression.tsv` | no regression outside the EXPECTED list (§5) |
| O5 | Build complete | snakemake "steps (100%) done"; C16 build integrity (MANIFEST, 0 stale; `depends_on_geometry` = this package) | yes |

Verdict: `ACCEPTED` if O0, O1, O2 (yeast; n/a for Jurkat), O3 gate, O4, O5 pass; otherwise `NOT ACCEPTED` with the
failing items. A failing item can be accepted only with a dated DECISION in PLANNING_NOTES (§4 no-regression rule).

**C15 for Jurkat: n/a (stated, not skipped silently).** On disk under `virtualcell/data`: no BioPlex, OpenCell IP-MS or
HuRI table of their own (`reference/opencell/` holds only `opencell-localization-annotations.csv`); BioGRID 5.0.263
human (`reference/tcell_crispr_screens/raw/networks/BIOGRID-ORGANISM-Homo_sapiens-5.0.263.tab3.txt`) contains the
BioPlex AP-MS records (Huttlin 2015 / 2017 / 2021: 23,712 / 56,473 / 169,808 rows) and the BioID cell map (Go 2021,
33,974 rows); IntAct human (`interactions/intact/Human.zip`) and `reference/complexes/omnipath_corum_complexes.tsv` are
there too. But there is no Jurkat model-neighbour code (the yeast C15 reads the yeast GPU contacts step; the common
contact counting comes with BD v2, plan §6.1; suite `proximity_na`). So the Jurkat B2 gate is O0, O1, O3, O4, O5 only:
**compartment totals + partner activity + no regression**. A Jurkat C15 adapter (BioGRID human, BioPlex/BioID subsets,
via the common contact counting) is future work, not part of B2.

**Circularity (both cells, reported in the verdict note).** STRING's "experiments" and "database" channels import
BioGRID / IntAct records, so a C15 gain from STRING-driven placement is partly circular (the reference sets feed the
network). The yeast B2 vs B1 comparison is nevertheless the clean measurement of the pool's contribution (STRING
unchanged; only the AF3 contact partners removed), and the AF3 pool is not derived from BioGRID.

## 5. EXPECTED changes (caused by B2; not regressions within these bounds)

- Tile RNG trajectories diverge everywhere (the counts split and the partner proposals change the draws): pack metrics
  differ at seed-noise level; judged by the §4 tolerances.
- Jurkat: per-tile copy numbers of each species follow the field (large-scale co-localisation of network neighbours,
  correlation length 500 nm); region totals unchanged (O1). Copies are placed next to STRING partners (O3).
- Yeast: contact partners (AF3 pool) are gone, loose partners kept; the partner fraction (O3) drops; C15 rows may move
  either way (judged by O2). The field (STRING only) is unchanged.
- C10 placed fractions, C7(a) clashes, C8(a) containment, C9, C11 densities, C12 ribosome NN, C1/C3/C6: the standard §4
  tolerances apply, no exemption (the layer changes where copies go inside a region, not how many are requested).
- Pack wall time (partner proposals: up to 40 extra tries per proposal; the counts step evaluates 128 fields per tile;
  the `pack_organisation` step itself): reported (B2 vs B1 `pack_pack`, `pack_counts`, `pack_organisation` rule times,
  same node / cores), not gated (the plan sets no time bound for B2).
