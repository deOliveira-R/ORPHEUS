r"""Pure albedo boundary condition.

See :class:`AlbedoBoundary` for the algebraic definition. The legacy
``AlbedoBoundaryOperator`` name is re-exported as a deprecated alias
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


__all__ = ["AlbedoBoundary"]


@dataclass(frozen=True)
class AlbedoBoundary(BoundaryTraceLaw, key="albedo"):
    r"""Pure albedo boundary: scalar multiple of the outgoing flux.

    Tensor decomposition :math:`(I, \alpha)` where :math:`I` is the
    angular identity and :math:`\alpha \in [0, 1]` is the albedo:

    .. math::

        \psi_{\text{in}}(\Omega) = \alpha \, \psi_{\text{out}}(\Omega).

    No angular redistribution. Useful as a *building block* for
    Wave-0 ``OperatorSum``-algebra mixed boundaries (e.g.
    ``c1 * AlbedoBoundary(0.5).realize(ms) + c2 * other.realize(ms)``,
    where albedo and specular shares are independent parameters),
    and as a stand-alone primitive when the boundary is a pure
    attenuator with no angular structure.

    Wave-7 rename note
    ------------------
    Previously named ``AlbedoBoundaryOperator``. The legacy name is
    preserved as a deprecated alias.

    Parameters
    ----------
    albedo : float
        Albedo coefficient. ``0`` is vacuum, ``1`` is perfect
        same-direction return.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    albedo: float = 0.0

    # ------------------------------------------------------------------
    # §16A.12 universal invariants — Wave 7 / C7.6 overrides.
    # ------------------------------------------------------------------

    def assert_response_positive_if_declared(self) -> None:
        r"""Albedo coefficient must be non-negative.

        Raises
        ------
        BoundaryResponseNotPositiveError
            When ``self.albedo < 0``.
        """
        if self.albedo < 0.0:
            from ._errors import BoundaryResponseNotPositiveError
            raise BoundaryResponseNotPositiveError(
                f"Albedo BC albedo={self.albedo} < 0",
                law="albedo",
            )

    def assert_submarkov(self) -> None:
        r"""Albedo satisfies the sub-Markov bound :math:`\alpha \le 1`.

        Raises
        ------
        SubmarkovViolationError
            When ``self.albedo > 1``.
        """
        if self.albedo > 1.0:
            from ._errors import SubmarkovViolationError
            raise SubmarkovViolationError(
                f"Albedo BC albedo={self.albedo} > 1",
                law="albedo",
            )

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Wave 7 (C7.3): delegate to the Wave-5 SNBoundaryRealizer.
        # The realizer dispatches on isinstance(self, AlbedoBoundary)
        # and returns one of:
        #
        # * ``ZeroOperator()`` if ``albedo == 0.0`` (vacuum-equivalent
        #   in the bare-attenuator algebra, but distinct from the
        #   §16A.5 inflow-only-mask of ``VacuumInflow``).
        # * ``IdentityOperator()`` if ``albedo == 1.0``.
        # * ``ScaledOperator(albedo, IdentityOperator())`` for
        #   ``albedo`` in (0, 1).
        #
        # The output is bit-equivalent to the pre-Wave-7 body
        # ``self.albedo * psi_out`` (legacy multiplies on input;
        # ScaledOperator multiplies the identity-passed-through input
        # — same single multiplication).
        #
        # Local imports break the cycle ``orpheus.geometry.boundary``
        # → ``orpheus.sn.boundary_realizer`` → ``orpheus.geometry.boundary``.
        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer,
            SNMethodSpace,
        )
        op = SNBoundaryRealizer().realize(
            self, SNMethodSpace.minimal(quadrature),
        )
        return op.apply(psi_out)
