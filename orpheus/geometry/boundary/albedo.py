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
        quadrature: "AngularQuadrature | None" = None,
    ) -> np.ndarray:
        r"""Scale the outgoing flux by :attr:`albedo`.

        No angular structure — :math:`\psi_{\text{in}}(\Omega) =
        \alpha \, \psi_{\text{out}}(\Omega)` for every ordinate. The
        ``quadrature`` argument is accepted for backward-compat with
        the legacy 2-arg form but is unused; the scalar
        multiplication needs no angular information.

        Issue #176 / C176.3: ``quadrature`` is now optional and
        ignored. The body is implemented inline (rather than
        delegating to the realizer) because the Albedo branch's
        realized output is just ``ScaledOperator(albedo,
        IdentityOperator())`` whose ``apply(psi)`` is ``albedo *
        psi`` — identical to inlining here, and avoids an import
        cycle in test paths that exercise direct-call.
        """
        return self.albedo * psi_out
