"""Machine-readable benchmark catalogue for LA-13511 (1999).

Each :class:`La13511Case` carries:

* Problem identifier in Sood naming convention
  (``<Material>-<Groups>-<Scattering>-<Geometry>``).
* Problem number from LA-13511 Tables 2-4.
* Geometry, group count, scattering order.
* Cross sections in **ORPHEUS convention**: ``g=0`` fast, ``g=N-1``
  slow. The conversion from Sood's reverse convention (``g=N`` fast,
  ``g=1`` slow) is handled by the catalogue at construction time —
  consumers see ORPHEUS-ordered arrays directly.
* Reference values: :math:`r_c` (mfp + cm), :math:`k_\\infty`,
  scalar-flux ratios at the 4 published sample points
  :math:`r/r_c \\in \\{0.25, 0.5, 0.75, 1.0\\}`.
* Provenance: which LA-13511 table holds the case and which primary
  reference the values are sourced from.
"""
from __future__ import annotations

from .la13511 import (
    La13511Case,
    PUA_1_0_IN,
    PU_2_0_IN,
    UA_1_0_SL_STUB,
    UA_1_0_CY_STUB,
    UA_1_0_SP_STUB,
    ALL_FIRST_SLICE,
)

__all__ = [
    "La13511Case",
    "PUA_1_0_IN",
    "PU_2_0_IN",
    "UA_1_0_SL_STUB",
    "UA_1_0_CY_STUB",
    "UA_1_0_SP_STUB",
    "ALL_FIRST_SLICE",
]
