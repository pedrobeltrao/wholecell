"""Yeast plugin (port of yeast build X5, 2026-11-16; docs/reference/history/yeast_x5_reproduction_criteria.md).

Registers:
  plugin_step      yeast.*         bridge.py: the yeast steps of the core stage list (config plugin_steps): the engine
                                   file layout inside a core build, the protein table, library species, the vendored
                                   validation scripts (A9)
  inventory_builder yeast.x5_ingredients  assembly.py: the assembly ingredient table as the pack table (A9)
  pipeline_engine  yeast_x5        (reproduction only: config/cells/yeast_x5_port.yaml) runner.YeastX5Engine: the vendored yeast stage code (engine/, git archive e3b7448 +
                                   PORT_PATCHES.md) run step by step inside a wholecell build (workflow/rules/engine.smk)
  pack_engine      yeast_rsa_v1    (reproduction only) runner.pack_yeast_rsa: the yeast seeded packer (stage 3), selected
                                   by pack.engine in yeast_x5_port.yaml
  geometry_engine  yeast_cell3d_v1 geom_engine: the yeast analytic-primitive + SDF-brick evaluator behind the
                                   cellgeom_v1 package written by toolkits/analytic/yeast.py
  validate_adapter yeast_x5        validate_adapter.YeastX5SuiteAdapter: reads the yeast validation outputs for the
                                   common suite (wholecell/validate/suite.py)
  preplaced        yeast.*         preplaced/: the yeast serial pre-placed builders (+ preplaced_context yeast.s3)
  exclusion_zones  yeast.endocytic_sites  exclusion.py: ribosome exclusion zones of the endocytic sites (core pack.exclusion)
Which yeast step is yeast-specific (plugin), generic (core candidate) or a toolkit is tabulated in steps.yaml.
Copied yeast code carries `# origin: yeast ... @ e3b7448`; no yeast licence yet (intended CC BY-NC 4.0): local only.
"""
from wholecell.registry import register

from . import assembly, bridge, exclusion, preplaced, runner, validate_adapter  # noqa: F401  (assembly: plugin steps + hooks of the core assembly preset legacy_x5;
#   preplaced: serial pre-placed builders + context yeast.s3; validate_adapter: suite adapter yeast_x5)

register("pipeline_engine", "yeast_x5")(runner.YeastX5Engine)
register("pack_engine", "yeast_rsa_v1")(runner.pack_yeast_rsa)
register("geometry_engine", "yeast_cell3d_v1")("plugins.yeast.geom_engine")
