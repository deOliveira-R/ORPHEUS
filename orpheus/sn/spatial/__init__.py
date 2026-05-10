"""Per-cell spatial-update strategies for the SN sweep.

This package hosts the strategy contract and the concrete
strategies — Diamond Difference (here) and (Wave C-extension)
Linear Discontinuous, Exponential Characteristic, Step — that a
1-D SN sweep uses to march through a cell and produce its
average flux + outgoing states.

The contract itself lives in :mod:`orpheus.sn.spatial.cell_update`:

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate` — the
  ``@runtime_checkable Protocol``.
* :class:`~orpheus.sn.spatial.cell_update.UpstreamState` —
  per-cell input state (spatial face flux + optional angular
  half-flux).
* :class:`~orpheus.sn.spatial.cell_update.CellResult` — per-cell
  output state (average flux, outgoing spatial flux, outgoing
  angular state).

Round 1 (Issue #157) shipped the contract.
Round 2 (Issue #158, this file's ``DiamondDifference`` re-export)
ships the first concrete strategy as a bit-identical extraction of
the existing inlined sweep math at :mod:`orpheus.sn.sweep`.  Wave
C-extension will add Linear Discontinuous, Step, and Exponential
Characteristic.

See the SN reshape campaign plan at
``.claude/plans/sn_reshape.md`` and the Wave C plan at
``.claude/plans/mossy-mapping-pine.md`` for context.
"""

from .boundary_face_flux import (
    BoundaryFaceFlux,
    BoundaryFaceFluxBase,
    CellCenter,
    DDExtrapolation,
)
from .cell_update import CellResult, CellUpdate, UpstreamState
from .diamond import DiamondDifference

__all__ = [
    "BoundaryFaceFlux",
    "BoundaryFaceFluxBase",
    "CellCenter",
    "CellResult",
    "CellUpdate",
    "DDExtrapolation",
    "DiamondDifference",
    "UpstreamState",
]
