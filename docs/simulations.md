# Simulations

> Release note (2026-11-23): the simulation applications run on the paper builds live in the separate `wholecell-paper`
> repository; this page is the original design note.

Simulations are **applications** that read a finished build. They are not stages of the build, and a build never
depends on them.

Nothing in this page is implemented in this repository yet; there is no `applications/` folder. The design is:

- Spec: [reference/simulation_workflows.md](reference/simulation_workflows.md). One Snakemake workflow per
  application under `applications/<name>/`, standard rules `check_build -> boxes -> setup -> run -> analyse ->
  prereg_check -> evaluate -> report`. Each application reads one pinned build (`build_id` plus the sha256 of its
  export manifest) and simulation boxes in the `wholecell.simbox/1` schema (written by the export stage, today a
  one-box stub).
- Pre-registration: test plans are frozen with an md5 before a run, and the `prereg_check` rule enforces it
  (spec §4).
- Re-run plan for the paper: [reference/history/unification_plan.md](reference/history/unification_plan.md) §6
  (BD diffusion, synapse, condensates, vesicle traffic, NFAT/calcium; GPU estimates in §6.6).
- GPU jobs (plan §6.7): each application ships a submit kit in `applications/<name>/slurm/` (`SUBMIT_README.md`,
  `jobs.tsv`, `submit_<batch>.sh`, `collect_verify.py`). We write the scripts; Pedro submits them. No agent monitors
  the queue.

Legacy simulation results stay where they are as the archival record and are compared with re-runs, never
overwritten.
