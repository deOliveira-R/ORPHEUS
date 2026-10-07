r"""The characteristic reference: transport integrated along the lines of a concentric body.

One reference posed from the specification (P1 of
``.claude/plans/characteristic_reference_architecture.md``), built on the
geometric kernel (:mod:`orpheus.geometry.chord`) and the reference
numerics kernel (:mod:`orpheus.derivations.common.dense_pencil`):

- :mod:`~orpheus.derivations.continuous.characteristic.walls` — each
  boundary point with what its law returns, read from the law's factors;
- :mod:`~orpheus.derivations.continuous.characteristic.closure` — the
  period of each line's unfolded path and the least solution of its cycle,
  the line part of the boundary resolvent;
- :mod:`~orpheus.derivations.continuous.characteristic.basis` — the panel
  basis the emission density and the flux are represented in;
- :mod:`~orpheus.derivations.continuous.characteristic.transport` — the
  transport along each line on that basis: the traversals' source
  integrals, the Volterra block and the angular flux;
- :mod:`~orpheus.derivations.continuous.characteristic.assembly` — the
  Galerkin assembly over lines: one group's transport block, its line part
  and its diffuse walls' coupling.
"""

from .assembly import GroupTransport, LineRule
from .basis import PanelBasis
from .closure import LinePeriod, TrappedSource, WallCoupling
from .transport import TraversalRule
from .walls import Wall, Walls

__all__ = [
    "GroupTransport",
    "LinePeriod",
    "LineRule",
    "PanelBasis",
    "TrappedSource",
    "TraversalRule",
    "Wall",
    "WallCoupling",
    "Walls",
]
