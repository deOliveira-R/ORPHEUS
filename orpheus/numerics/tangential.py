r"""The half-range of a direction cosine: the one selector of outflow, inflow and tangential ordinates.

For a face whose outward normal is the axis ``axis_index`` taken with a ``sign``, an ordinate :math:`n` is *outflow* when
:math:`\mathrm{sign}\,\mu_n > 0`, *inflow* when it is :math:`< 0` and *tangential* when it is exactly :math:`0`. The comparisons are
exact: no tolerance, band or epsilon decides a class. The tangential set and the kernel of the form :math:`\mathrm{Diag}(|\mu_n| w_n)`
are one datum on every shipped rule, because a rule carries no zero weight on a non-tangential ordinate, and a tolerance on the cosine
would split them. Every shipped tangential cosine is exactly ``0.0`` (the azimuths are roots of unity); a generator that emitted
round-off cosines would have those ordinates read as genuine outflow or inflow, and the gate of
:func:`half_range` on a rotated rule is the control that catches it.

A leaf of ``numerics`` whose only intra-package import is the declaration module, so that the angular trace spaces, the quadrature
certification (:func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing`) and the octant partition of
:class:`~orpheus.numerics.quadrature.Quadrature` all read it without ``quadrature`` importing ``spaces`` (the module-cycle gate counts
every import kind). :func:`half_range` types its argument by the protocol :class:`AxisCosineSource` for the same reason.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from orpheus.numerics.declaration import Coord

__all__ = ["AxisCosineSource", "HalfRange", "half_range"]


class AxisCosineSource(Protocol):  # concept: quadrature | role: interface | kind: element
    r"""Anything that answers the orbit-mean direction cosine of every ordinate along an axis."""

    def mean_axis_cosine(self, axis_index: int) -> np.ndarray:
        r"""The ``(N,)`` array of cosines :math:`\langle \Omega_i \rangle`, one per ordinate."""
        ...


@dataclass(frozen=True)
class HalfRange:  # concept: trace | role: part | kind: set
    r"""A partition of the ordinate index set ``range(n)`` into outflow, inflow and tangential ordinates.

    Each part is a :class:`~orpheus.numerics.declaration.Coord` on the ordinate set. The parts are pairwise disjoint and cover
    ``range(n)``; a value that does not is refused at construction.
    """

    outflow: Coord
    inflow: Coord
    tangential: Coord
    n: int

    def __post_init__(self) -> None:
        for name in ("outflow", "inflow", "tangential"):
            part = getattr(self, name)
            if part.n != self.n:
                raise ValueError(f"the {name} part is not on the ordinate set: Coord.n = {part.n}, n = {self.n}")
        out, inn, tan = (set(p.indices) for p in (self.outflow, self.inflow, self.tangential))
        if out & inn or out & tan or inn & tan:
            raise ValueError("the outflow, inflow and tangential parts are not disjoint")
        if len(out) + len(inn) + len(tan) != self.n:
            raise ValueError("the outflow, inflow and tangential parts do not cover the ordinates")

    def sign(self) -> np.ndarray:
        r"""The ``(n,)`` int8 array that is ``+1`` on ``outflow``, ``-1`` on ``inflow`` and ``0`` on ``tangential``."""
        out = np.zeros(self.n, dtype=np.int8)
        out[list(self.outflow.indices)] = 1
        out[list(self.inflow.indices)] = -1
        return out


def half_range(quadrature: AxisCosineSource, axis_index: int, *, sign: int = +1) -> HalfRange:
    r"""The half-range of ``quadrature`` along axis ``axis_index``, outward being the direction of ``sign``.

    Reads ``sign * quadrature.mean_axis_cosine(axis_index)``. Outflow is the set where that product is strictly positive, inflow
    where it is strictly negative, tangential where it is exactly zero.

    Raises
    ------
    ValueError
        If ``sign`` is not ``+1`` or ``-1``, if the cosines are not a one-dimensional array, or if a cosine is not finite.
    """
    if sign not in (+1, -1) or isinstance(sign, bool):
        raise ValueError(f"sign must be +1 or -1, got {sign!r}")
    c = sign * np.asarray(quadrature.mean_axis_cosine(axis_index), dtype=np.float64)
    if c.ndim != 1:
        raise ValueError(f"the direction cosines must be a one-dimensional array, got shape {c.shape}")
    if not np.all(np.isfinite(c)):
        raise ValueError(f"non-finite direction cosine along axis {axis_index}: ordinates {np.flatnonzero(~np.isfinite(c)).tolist()}")
    n = int(c.size)

    def part(mask: np.ndarray) -> Coord:
        return Coord(n, tuple(int(i) for i in np.flatnonzero(mask)))

    return HalfRange(part(c > 0.0), part(c < 0.0), part(c == 0.0), n)
