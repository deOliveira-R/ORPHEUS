r"""Energy group-structure instances — the canonical grids for energy condensation.

A growing family of multigroup :class:`~orpheus.data.energy_grid.EnergyGrid` instances
(the grid *data*; the TYPE + within-group flux models stay in
:mod:`orpheus.data.energy_grid`). The condensation source is the library fine grid
(:data:`~orpheus.data.group_structures.ornl.ORNL_421`); the targets are the WIMS-D
few-group structures (:mod:`~orpheus.data.group_structures.wims`).
"""

from __future__ import annotations

from orpheus.data.group_structures.ornl import ORNL_421
from orpheus.data.group_structures.wims import (
    CONDENSE_172_TO_69,
    WIMS_69,
    WIMS_172,
)

__all__ = ["ORNL_421", "WIMS_69", "WIMS_172", "CONDENSE_172_TO_69"]
