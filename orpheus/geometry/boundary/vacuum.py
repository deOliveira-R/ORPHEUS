r"""Vacuum (prescribed-zero inflow) boundary law.

See :class:`VacuumInflow` for the algebraic definition. This class was
previously named ``VacuumBoundaryOperator``; the legacy alias was
retired in Wave O step O.4a.1. ``VacuumInflow`` (the Grand Report v3
vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._base import BoundaryTraceLaw


__all__ = ["VacuumInflow"]


@dataclass(frozen=True)
class VacuumInflow(BoundaryTraceLaw, key="vacuum"):
    r"""Vacuum boundary: :math:`R = 0`, :math:`q = 0`.

    The empty sum in the tensor decomposition: no incoming flux,
    irrespective of what leaks out. Algebraically a rank-0 case of
    :eq:`bc-tensor-decomposition`.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. The SN realisation is an
    :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
    that zeroes the per-face inflow ordinates only (the §16A.5
    trace-correct semantics). Realise via:

    .. code-block:: python

        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer, SNMethodSpace,
        )
        op = SNBoundaryRealizer().realize(
            VacuumInflow(),
            SNMethodSpace(
                quadrature=quad, face="right",
                inflow_indices=np.flatnonzero(quad.mu_x < 0),
            ),
        )
        psi_in = op.apply(psi_out)

    Rename history
    --------------
    Previously named ``VacuumBoundaryOperator``. The Grand Report v3
    naming convention (``VacuumInflow`` — "the prescribed-zero inflow
    law") is now canonical. The legacy alias was retired in Wave O
    step O.4a.1 — ``VacuumInflow`` is the sole importable name.

    The :attr:`kind` attribute stays ``"vacuum"`` (the registry key
    under which this class is indexed) for backward compat with
    string-kind comparisons (``sn_mesh.bc_right == "vacuum"``).
    """

    #: String tag for legacy string-kind comparisons. Preserved
    #: across the Issue #186 descriptor cleanup so the SN-side
    #: ``sn_mesh.bc_right == "vacuum"`` test contract still holds.
    kind: str = "vacuum"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, VacuumInflow):
            return True
        return NotImplemented

    def __hash__(self) -> int:
        # Hash on the canonical class name. All instances are
        # value-equal (the descriptor is stateless), so they share
        # one hash bucket.
        return hash(("VacuumInflow",))
