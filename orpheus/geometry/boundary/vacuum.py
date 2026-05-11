r"""Vacuum boundary condition.

See :class:`VacuumBoundaryOperator` for the algebraic definition.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY

from ._base import BoundaryOperator

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["VacuumBoundaryOperator"]


@dataclass(frozen=True)
class VacuumBoundaryOperator(BoundaryOperator, key="vacuum"):
    r"""Vacuum boundary: :math:`R = 0`.

    The empty sum in the tensor decomposition: no incoming flux,
    irrespective of what leaks out. Algebraically a rank-0 case of
    :eq:`bc-tensor-decomposition`.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    #: String tag for legacy string-kind comparisons. The Wave B
    #: refactor preserves the existing BC-resolution test contract
    #: (``sn_mesh.bc_right == "vacuum"`` continues to evaluate True)
    #: while consumers transition to direct
    #: :meth:`apply` calls.
    kind: str = "vacuum"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, VacuumBoundaryOperator):
            return True
        return NotImplemented

    def __hash__(self) -> int:
        return hash(("VacuumBoundaryOperator",))

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        return np.zeros_like(psi_out)
