# B3 Jurkat library rebuild: addendum 2026-11-17 (criterion L5, belt measure)

Addendum to `docs/unification_criteria/B3_library.md` (frozen, md5 `b6c0933613d977a7770f109ca7d2c914`, not edited).
Decided by the main session 2026-11-17 (Pedro informed). **Written AFTER seeing the B3 library-only result**, so it
is post hoc: every B3 report gives both values, "L5 (frozen definition)" and "L5 (addendum 2026-11-17)". Only the
belt measure of L5 changes. The threshold (0.5 nm), the pass rule, and every other criterion (L1-L4, L6-L8, Y1-Y3,
B1-B4, the EXPECTED list, the verdict rule) stay as frozen.

## Why

Library-only build `jurkat/_dev_unify/b3/run_59634d8` (`compare/library_compare.json`): L1-L4, L6 and L7 pass. L5
fails: 340 of 3,531 TM species (9.6 %) have |z_belt| > 0.5 nm. The frozen v7 library gives **the same 340 species**
(all 3,531 TM species have bitwise-identical beads in both libraries), so the failure is in the measurement, not in
the rebuild:

- 3,409 TM species have no `tm_beads` key, so the frozen belt is "all beads with |z| < 2.0 nm". For asymmetric
  proteins (a large soluble domain on one side) the soluble beads that reach into that slab pull the median off
  centre: 309 of the 340 failures. The builders put the median (assemblies, `membrane_orient.orient`) or mean
  (monomers, `build_structs.process_af`) of the UniProt TRANSMEM residues at z = 0. They never centre a 2 nm slab of
  all beads.
- 31 failures have `tm_beads`: 27 are INTRAMEM-only (monotopic) species, whose hairpins sit in one leaflet by
  design (`jlib_tm_route` anchored, z 0.3-1.9 nm). The other 4 are two-membrane assemblies (GJA1_x12, GJB2_x12,
  VANGL1_x6, VANGL2_x6), where the slab comes from one layer of chains.
- Frozen-definition |z_belt| on both libraries: median 0.13 nm, 90 % 0.49, 95 % 0.78, 99 % 1.30 nm.

The plan says "hydrophobic belt centred within 0.5 nm". The hydrophobic belt means the membrane-spanning segments,
not every bead that happens to lie near the mid-plane.

## New definition: L5 (addendum 2026-11-17)

z_belt = median z (library frame, mid-plane z = 0) of the **membrane-segment beads** of each TM species (the
`membrane` key is present):

1. **`tm_beads` key present** (TM-routed species): the beads flagged 1, as stored (TRANSMEM + INTRAMEM, structural
   or built).
2. **Otherwise**: the beads of the UniProt TRANSMEM residues, plus beta-barrel strands as in
   `membrane_orient.masks`. The bead-to-residue map is re-derived from the species' structure source, read exactly
   as the builder read it:
   - `afdb_monomer`: the AFDB v6 F1 model through the tar index, keeping pLDDT ≥ 50 or TRANSMEM residues (all
     residues if fewer than 10 are kept);
   - `pdb`: `jlib_build.beads_pdb` (SIFTS), with `tm_graft` for TM assemblies;
   - `afdb_complex`: `jlib_build.beads_afdb`;
   - TM ranges come from `uniprot_topology.tsv`.

   The map is accepted on a given library only if it has the same bead count and the pairwise distances of up to
   200 sampled TRANSMEM beads match the source within 0.02 nm. If not, the species falls back to the frozen
   definition (`fallback_slab`, listed).
3. **Not applicable**: species with `tm_beads` where no accession has a UniProt TRANSMEM feature (INTRAMEM-only /
   monotopic, from the tm_route UniProt REST cache). These are listed and excluded from the denominator. Their belt
   is one leaflet by construction. Also reported: the fraction when they are counted as failures.

**Pass (unchanged):** |z_belt| ≤ 0.5 nm for ≥ 99 % of the scored TM species, and the count failing ≤ the v7 count
(same measure on the frozen v7 library) + 1 % of the scored TM species.

Note: the builders centre exactly these residues, so the measure mostly checks that the stored frame is the one the
TM annotation defines. It cannot detect a wrong annotation. Code: `scripts/library_b3_check.py` (`L5a`;
`rescore` sub-command). `compare` and `verdict` report both L5 and L5a, and the pass rule and the verdict use L5a.
The TRANSMEM map is cached as `l5a_tm_segments.pkl.gz` next to `measures_v7.tsv`.

## Re-score of the library-only result (`jurkat/_dev_unify/b3/rescore_addendum/`)

| | rebuilt (run_59634d8) | frozen v7 |
|---|---|---|
| L5 (frozen definition): failing / 3,531, fraction ok | 340, 0.9037 (FAIL) | 340, 0.9037 |
| L5 (addendum 2026-11-17): failing / 3,504 scored, fraction ok | 21, 0.9940 (PASS) | 21, 0.9940 |

- Belt source on both libraries: 3,409 `tm_segments` (all maps verified, 0 fallbacks), 95 `tm_beads`, 27 not
  applicable (INTRAMEM-only).
- Counting the 27 not-applicable species as failures: 0.9864.
- z_belt is identical per species on both libraries. Median |z_belt| 0.009 nm.
- The 21 failures: the 4 two-membrane assemblies (by design), NPC (−1.06 nm, NE ring axis), MCUR1 (−5.0 nm), and
  15 monomers at 0.5-1.2 nm. In those monomers, process_af centres the TRANSMEM mean, not the median, and the TM
  annotation is asymmetric or partly peripheral.

Library-level B3 (L1-L7) under the addendum: **pass** (L1-L4, L6, L7 as in `library_compare.json`; L5a passes).
L8 and Part B are decided by the full build.
