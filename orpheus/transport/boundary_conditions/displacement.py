r"""Combinations outside the cone (``displacement.py`` of the conditions package).

**What a displacement is.** The cone :math:`K` of walls (``s, d >= 0``, ``s + d <= 1``, every inflow weight ``>= 0``)
is a convex set; a combination of conditions that leaves it (``white(0.3) - specular(0.5)``, ``0.6 * specular +
0.6 * white``, a negative multiple) is a point of the AFFINE HULL of :math:`K` times the source space: a direction in
which a parameter may move, the sensitivity coordinates of a perturbation study, not a wall. :class:`WallDisplacement`
is that value. It is NOT a :class:`~.condition.BoundaryCondition`, so it is never bindable: ``BoundaryConditions``
refuses it with a ``TypeError`` naming the type, and that refusal is the one narrowing site (no consumer branches on
the sort).

**Why a separate type and not a flag.** The cone is a value predicate (decision 50): the result of arithmetic is a
``WallCondition`` iff its coordinates lie in :math:`K`. Weights of a displacement are signed, which a ``WallCondition``
cannot hold, so the two cannot share a type without a runtime flag that every consumer would have to read.

**Same operators, one exit.** A displacement defines ``+``, ``-``, unary ``-`` and scalar ``*`` through the same
``_PairArithmetic`` as the condition, exiting through ``_combine``: ``WallCondition + WallDisplacement`` is the
condition of the summed pair iff the sum re-enters :math:`K` (a finite-difference step and its inverse returns to a
condition), else a displacement.

**Gotchas.**

* :meth:`WallDisplacement.of` ALWAYS builds a displacement, even when the combination lies in :math:`K` (then it is
  still unbindable and unequal to the condition of the same coordinates, since the types differ). It is the
  constructor tests use when the operand signs must not decide the type.
* The network and the inflow of a displacement hold the same drop rule as a condition's (a zero weight is not
  stored, ``-0.0`` is ``0.0``).

Layer: ``transport`` (L2). Imports ``values`` and the sibling modules only.
"""

from __future__ import annotations

from dataclasses import dataclass

from orpheus.transport.boundary_conditions.condition import WallCondition, _PairArithmetic
from orpheus.transport.boundary_conditions.inflow import InflowTerms
from orpheus.transport.boundary_conditions.network import ResponseNetwork
from orpheus.values.content import ContentIdentity, FrozenMapping

__all__ = ["WallDisplacement"]


@dataclass(frozen=True, eq=False)
class WallDisplacement(_PairArithmetic, ContentIdentity):
    """A combination outside the cone: a point of the affine hull of the wall cone times the source space. Never bindable."""

    network: ResponseNetwork
    inflow: InflowTerms

    @classmethod
    def of(cls, *terms: tuple[float, WallCondition]) -> WallDisplacement:
        """The signed combination ``sum_i c_i x_i`` of conditions, ALWAYS a displacement."""
        network = ResponseNetwork(FrozenMapping({}))
        inflow = InflowTerms.empty()
        for coefficient, condition in terms:
            network = network + condition.network * coefficient
            inflow = inflow + condition.inflow * coefficient
        return cls(network, inflow)
