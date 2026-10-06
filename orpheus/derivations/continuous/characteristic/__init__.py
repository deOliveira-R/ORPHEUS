r"""The characteristic reference: transport integrated along the lines of a concentric body.

One reference posed from the specification (P1 of
``.claude/plans/characteristic_reference_architecture.md``), built on the
geometric kernel (:mod:`orpheus.geometry.chord`) and the reference
numerics kernel (:mod:`orpheus.derivations.common.dense_pencil`):

- :mod:`~orpheus.derivations.continuous.characteristic.walls` — each
  boundary point with what its law returns, read from the law's factors;
- :mod:`~orpheus.derivations.continuous.characteristic.closure` — the
  period of each line's unfolded path and the least solution of its cycle,
  the line part of the boundary resolvent.
"""

from .closure import LinePeriod, TrappedSource
from .walls import Wall, Walls

__all__ = ["LinePeriod", "TrappedSource", "Wall", "Walls"]
