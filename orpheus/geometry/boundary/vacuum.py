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
        quadrature: "AngularQuadrature | None" = None,
    ) -> np.ndarray:
        r"""Return zeros (backward-compat fallback).

        Issue #176 / C176.3: ``quadrature`` is now optional. The
        §16A.5-correct inflow-only-mask requires per-face inflow
        indices that this signature cannot deliver — for the
        production-correct path, route through
        :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`
        with an :class:`~orpheus.sn.method_space.SNMethodSpace`
        carrying the face's inflow indices:

        .. code-block:: python

            from orpheus.sn.boundary_realizer import SNBoundaryRealizer
            from orpheus.sn.method_space import SNMethodSpace
            method_space = SNMethodSpace(
                quadrature=quad, face="right",
                inflow_indices=np.flatnonzero(quad.mu_x < 0),
            )
            op = SNBoundaryRealizer().realize(VacuumInflow(), method_space)
            psi_in = op.apply(psi_out)

        Direct calls keep the legacy zeros-all body. The two paths
        agree on the inflow rows (production-relevant subset for SN
        sweeps); they diverge on the outflow rows. ``quadrature`` is
        accepted for backward-compat with the legacy 2-arg form but
        is unused — the zeros-all body needs no angular information.
        """
        return np.zeros_like(psi_out)
