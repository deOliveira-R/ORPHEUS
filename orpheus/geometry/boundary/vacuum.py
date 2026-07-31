r"""Vacuum (prescribed-zero inflow) boundary law.

See :class:`VacuumInflow` for the algebraic definition. This class was
previously named ``VacuumBoundaryOperator``; the legacy alias was
retired in Wave O step O.4a.1. ``VacuumInflow`` (the Grand Report v3
vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._base import BoundaryTraceLaw
from ._factors import NullMap, ScalarResponse


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

        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        op = SNBoundaryRealizer().realize(
            VacuumInflow(),
            SNMethodSpace(
                quadrature=quad, face="xmax",
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

    :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.kind` is ``"vacuum"`` —
    inherited from the base, which derives it from the registry key. It was a
    mutable dataclass FIELD here until the B0 cleanup, which meant
    ``VacuumInflow(kind="banana")`` constructed a law whose tag matched no
    registry entry; deriving it removes that state. The string-comparison
    contract (``sn_mesh.bc["xmax"] == "vacuum"``) is unchanged.
    """

    # ── The affine form's two factors (B1) ──────────────────────────────
    @property
    def geometry_map(self) -> "NullMap":
        r""":math:`G = 0` — no outgoing flux returns.

        Vacuum is the rank-0 corner of :math:`\gamma_-\psi = R\,G\,\gamma_+\psi
        + q`: geometry, response AND source all vanish, so
        :math:`\gamma_-\psi = 0`. That the realized SN operator is a *projector*
        rather than the zero map is a representation choice at the operator
        layer (it preserves the outflow rows), not a different law — and phase
        **B3** narrows the domain so the two agree.
        """
        return NullMap()

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = 0` — the defining property of a vacuum face."""
        return ScalarResponse(0.0)

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
