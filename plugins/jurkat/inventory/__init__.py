"""Jurkat inventory components (registered on import): adapter.rna jurkat.rna_e61, adapter.proteome jurkat.geiger_ibaq_v0,
adapter.localisation jurkat.hpa_uniprot_itzhak_v0, inventory_builder
jurkat.ingredients_v0, inventory_rule jurkat.v4 / jurkat.v5 / jurkat.v6 / jurkat.v7. Curated lists: plugins/jurkat/data/.
The assembly state is the core step list (wholecell/inventory/assembly, preset legacy_v7; since 2026-11-16 A4); the
Jurkat curation it reads is plugins/jurkat/data/assembly/*.tsv."""
from . import rna_e61, inventory_v0, ingredients, rules_v4, rules_v5, rules_v6, rules_v7  # noqa: F401
