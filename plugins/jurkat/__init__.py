"""Jurkat plugin: cell-specific components (inventory rules v4-v7, virtual regions, curated lists). Importing this
package registers them (wholecell.registry)."""

from . import virtual_regions  # noqa: F401  (virtual_region jurkat.cytoskeleton_v7; owner: pack)
from . import inventory  # noqa: F401  (inventory adapters, builder, rules v4-v7; owner: inventory)
from . import preplaced  # noqa: F401  (preplaced_data providers of the common pre-placed builders; B4)
