r"""Vacuum (prescribed-zero inflow) boundary law.

See :class:`VacuumInflow` for the algebraic definition. The legacy
``VacuumBoundaryOperator`` name is re-exported as a deprecated alias
from the package ``__init__.py`` (Wave 7 rename per Grand Report v3
vocabulary).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["VacuumInflow"]


@dataclass(frozen=True)
class VacuumInflow(BoundaryTraceLaw, key="vacuum"):
    r"""Vacuum boundary: :math:`R = 0`, :math:`q = 0`.

    The empty sum in the tensor decomposition: no incoming flux,
    irrespective of what leaks out. Algebraically a rank-0 case of
    :eq:`bc-tensor-decomposition`.

    Wave-7 rename note
    ------------------
    Previously named ``VacuumBoundaryOperator``. The Grand Report v3
    naming convention (``VacuumInflow`` — "the prescribed-zero inflow
    law") is now canonical. The legacy name is preserved as a
    deprecated alias in ``orpheus.geometry.boundary``; remove the
    alias in a future cleanup wave.

    The :attr:`kind` attribute stays ``"vacuum"`` (the registry key
    under which this class is indexed) for backward compat with
    string-kind comparisons (``sn_mesh.bc_right == "vacuum"``).
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
        if isinstance(other, VacuumInflow):
            return True
        return NotImplemented

    def __hash__(self) -> int:
        # Hash on the canonical (post-rename) class name; the legacy
        # ``VacuumBoundaryOperator`` alias resolves to this class so
        # ``hash(VacuumBoundaryOperator()) == hash(VacuumInflow())``.
        return hash(("VacuumInflow",))

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Wave 7 design note: the 2-arg legacy path keeps zeros-all
        # for backward compat — production sweeps still call
        # ``bc.apply(psi_out, quadrature)`` and the §16A.5
        # inflow-only-mask requires per-face inflow indices that this
        # signature can't deliver. The §16A.5 trace-correct
        # behaviour activates when consumers route through
        # :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`
        # (which receives an :class:`SNMethodSpace` carrying the
        # inflow indices) — see the realizer's vacuum dispatch
        # branch. Wave 8's SNMesh wiring will hoist this delegation
        # to the realizer, at which point this 2-arg body becomes
        # dead code and is dropped in Wave 11.
        return np.zeros_like(psi_out)
