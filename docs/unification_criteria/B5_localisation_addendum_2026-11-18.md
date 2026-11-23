# B5 localisation: criteria addendum 2026-11-18 (gate on independent references)

Date: 2026-11-18. **Post-hoc change, made after seeing the B5 revision dry-run results, approved by Pedro on
2026-11-18.** The frozen criteria `B5_localisation.md` (md5 49a43723bfafb90ef5092b96f42bbf20) are not edited; this file
amends how their C1 gate is applied. Revision of the defaults: `B5_localisation_revision_2026-11-18.md`
(md5 98e239af5bcd1b47f3868900c2fe1b19).

## Change

The B5 localisation gate (C1 agree ≥ reference − 0.005; image N/C as frozen) applies to the **independent** references
only:

| Reference | Role | Why |
|---|---|---|
| OpenCell (live HEK293T, no antibody) | **gated** | independent of every localisation input |
| HPA-Jurkat **hold-out half** | **gated** | held out from the N/C calibration |
| image N/C correlation (C3, hold-out) | **gated** | as frozen |
| HPA consensus | reported, **not gated** | partly circular: HPA calls are inputs of the localisation |
| HPA-Jurkat full set | reported, **not gated** | partly circular (contains the fit half); the HPA antibody-background correction (b = 0.19, `data/derived/v8/hpa_background/estimate.json`) deliberately moves the model away from raw HPA calls |

All four C1 sets are always reported, with their gate status.

## Result under this addendum (dry run `jurkat/_dev_unify/b5/dryrun_7b9880f/`)

OpenCell +0.0055 (pass), HPA-Jurkat hold-out −0.0036 (pass), image N/C hold-out 0.606 (pass); not gated: HPA consensus
−0.0081, HPA-Jurkat full −0.0073. Yeast §3d passes (SWAT 0.6647, hyperLOPIT 0.8483, image N/C 0.833). **B5 accepted
(Pedro, 2026-11-18)** for the combined Step-2 build; the build-level check uses this addendum.
