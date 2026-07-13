"""Method-generic spatial discretization schemes — the per-cell closure layer.

This package hosts the spatial trial-space/closure contract and the
concrete schemes — Diamond Difference and Linear Discontinuous — that a
transport method consumes to march through a cell.  The layer is
METHOD-GENERIC by construction: geometry arrives as data
(:class:`~orpheus.geometry.reduced_operator.StreamingTerms`), angular
closure constants arrive as plain floats on the
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet (#236's
spatial ⊗ angular separation), and nothing here imports from a method
package.  The discrete-ordinates sweep was the first consumer, not the
owner — the layer was promoted from ``orpheus.sn.sweep`` at #272
(2026-07-04) so the assembly consumption mode can serve every method.

The contract lives in :mod:`orpheus.transport.spatial.scheme`:

* :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` — the
  ``@runtime_checkable Protocol``.
* :class:`~orpheus.transport.spatial.scheme.UpstreamState` — per-cell
  input state (spatial face flux + optional angular half-flux).
* :class:`~orpheus.transport.spatial.scheme.CellResult` — per-cell
  output state (average flux, outgoing spatial flux, outgoing angular
  state).

History: Round 1 (Issue #157) shipped the contract; Round 2 (Issue
#158) shipped ``DiamondDifference`` as a bit-identical extraction of
the then-inlined sweep math; the UBLD kernels + Linear Discontinuous
followed (#240 family); #272 promoted the layer here.  The open #158
arms (Step, Exponential Characteristic, the curvilinear LD derivation)
land in this package.
"""

from .scheme import CellResult, DiscretizationScheme, UpstreamState
from .diamond import DiamondDifference
from .linear_discontinuous import LinearDiscontinuous

__all__ = [
    "CellResult",
    "DiamondDifference",
    "DiscretizationScheme",
    "LinearDiscontinuous",
    "UpstreamState",
]
