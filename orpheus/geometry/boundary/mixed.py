r"""Mixed (rank-N) boundary condition.

See :class:`MixedBoundaryOperator` for the algebraic definition.
This composer is scheduled for removal in Wave 11 (replaced by Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` algebra).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, ClassVar, Sequence

import numpy as np

from orpheus.numerics.operator import CAP_APPLY

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["MixedBoundaryOperator"]


@dataclass(frozen=True)
class MixedBoundaryOperator(BoundaryTraceLaw, key="mixed"):
    r"""Linear combination of rank-1 BC primitives.

    Realises a rank-N tensor decomposition

    .. math::

        R = \sum_\alpha c_\alpha \, G_\alpha \otimes A_\alpha

    by storing a list of ``(coefficient, primitive)`` pairs and
    summing ``coefficient * primitive.apply(...)``. The
    coefficients :math:`c_\alpha` are typically convex (sum to 1) for
    a partial-current Marshak boundary (Bell & Glasstone 1970 §1.5:
    ``MixedBoundaryOperator([(0.3, SpecularBoundaryOperator()),
    (0.7, WhiteBoundaryOperator())])`` is "30% specular, 70% diffuse");
    the linear-combination interface does not enforce this so other
    use cases (asymmetric weights, gain media) can also be expressed.

    Parameters
    ----------
    components : sequence of (coefficient, BoundaryOperator)
        The rank-N decomposition. Each component contributes
        ``coefficient * primitive.apply(...)`` to the incoming flux.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    components: tuple[tuple[float, BoundaryTraceLaw], ...] = field(
        default_factory=tuple
    )

    def __init__(
        self,
        components: Sequence[tuple[float, BoundaryTraceLaw]],
    ) -> None:
        # Frozen-dataclass-with-Sequence-arg pattern: take a Sequence,
        # store as a tuple. ``object.__setattr__`` to bypass the frozen
        # guard during construction.
        object.__setattr__(self, "components", tuple(components))

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        result = np.zeros_like(psi_out)
        for coeff, primitive in self.components:
            result = result + coeff * primitive.apply(psi_out, quadrature)
        return result
