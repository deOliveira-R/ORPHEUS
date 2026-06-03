r"""Pure albedo boundary condition.

See :class:`AlbedoBoundary` for the algebraic definition. This class
was previously named ``AlbedoBoundaryOperator``; the legacy alias was
retired in Wave O step O.4a.1. ``AlbedoBoundary`` (the Grand Report
v3 vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass

from ._base import BoundaryTraceLaw


__all__ = ["AlbedoBoundary"]


@dataclass(frozen=True)
class AlbedoBoundary(BoundaryTraceLaw, key="albedo"):
    r"""Pure albedo boundary: scalar multiple of the outgoing flux.

    Tensor decomposition :math:`(I, \alpha)` where :math:`I` is the
    angular identity and :math:`\alpha \in [0, 1]` is the albedo:

    .. math::

        \psi_{\text{in}}(\Omega) = \alpha \, \psi_{\text{out}}(\Omega).

    No angular redistribution. Useful as a *building block* for
    descriptor-tree composition (e.g.
    ``0.5 * AlbedoBoundary(1.0) + 0.5 * ReflectiveBoundary(axis="x")``)
    and as a stand-alone primitive when the boundary is a pure
    attenuator with no angular structure.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. The SN realisation collapses to
    :class:`~orpheus.numerics.operator.ZeroOperator` (α=0),
    :class:`~orpheus.numerics.operator.IdentityOperator` (α=1), or
    ``ScaledOperator(α, IdentityOperator)`` (α ∉ {0, 1}). Realise via
    :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`.

    Rename history
    --------------
    Previously named ``AlbedoBoundaryOperator``. The legacy alias was
    retired in Wave O step O.4a.1 — ``AlbedoBoundary`` is the sole
    importable name.

    Parameters
    ----------
    albedo : float
        Albedo coefficient. ``0`` is vacuum, ``1`` is perfect
        same-direction return.
    """

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
