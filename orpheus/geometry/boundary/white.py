r"""White (Lambertian) boundary condition.

See :class:`WhiteBoundary` for the algebraic definition. This class
was previously named ``WhiteBoundaryOperator``; the legacy alias was
retired in Wave O step O.4a.1. ``WhiteBoundary`` (the Grand Report v3
vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    import numpy as np

    from orpheus.numerics.quadrature import Quadrature


__all__ = ["WhiteBoundary"]


@dataclass(frozen=True)
class WhiteBoundary(BoundaryTraceLaw, key="white"):
    r"""White (Lambertian) boundary with optional albedo.

    Tensor decomposition :math:`(G_{\text{diff}}, \alpha)` where
    :math:`G_{\text{diff}}` is the cosine-weighted angular average
    over the outgoing hemisphere, broadcast isotropically over the
    incoming hemisphere:

    .. math::

        \psi_{\text{in}}(\Omega) = \frac{\alpha}{\pi}
            \sum_{\Omega' :\, \Omega' \cdot \hat{n} > 0}
            w(\Omega')\,|\Omega' \cdot \hat{n}|\,
            \psi_{\text{out}}(\Omega'),

    where the result is independent of the incoming :math:`\Omega`
    direction (Lambertian / cosine emission). The factor :math:`1/\pi`
    is the canonical Lambertian normalisation used in radiative
    transfer; the implementation here normalises by the outgoing
    cosine-weighted weight sum so the BC is conservative for any
    quadrature: the total returned current equals the incoming
    current (times :math:`\alpha`), which is the property the consumer
    actually needs -- see Bell & Glasstone 1970 §1.5.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` method. The SN realisation is
    :class:`~orpheus.sn.boundary.angular.AngularAverageOperator`
    (α=1 fast path) or ``ScaledOperator(α, AngularAverageOperator)``
    (α ≠ 1). Realise via:

    .. code-block:: python

        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        law = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.8)
        op = SNBoundaryRealizer().realize(law, SNMethodSpace.minimal(quad))
        psi_in = op.apply(psi_out)

    Rename history
    --------------
    Previously named ``WhiteBoundaryOperator``. The legacy alias was
    retired in Wave O step O.4a.1 — ``WhiteBoundary`` is the sole
    importable name.

    Parameters
    ----------
    axis : str
        Boundary normal axis: ``"x"``, ``"y"``, or ``"z"``.
    outward_sign : int
        Sign of the outward normal along ``axis``: ``+1`` for the
        upper face (``xmax`` / ``ymax``) and ``-1`` for the lower face
        (``xmin`` / ``ymin``). Selects which ordinates are *outgoing*
        at this face.
    albedo : float
        Diffuse albedo. Defaults to 1 (perfectly reflecting).
    """

    axis: str = "x"
    outward_sign: int = +1
    albedo: float = 1.0

    # ------------------------------------------------------------------
    # §16A.12 universal invariants — Wave 7 / C7.6 overrides.
    # ------------------------------------------------------------------

    def assert_response_positive_if_declared(self) -> None:
        r"""White albedo must be non-negative.

        Raises
        ------
        BoundaryResponseNotPositiveError
            When ``self.albedo < 0``.
        """
        if self.albedo < 0.0:
            from ._errors import BoundaryResponseNotPositiveError
            raise BoundaryResponseNotPositiveError(
                f"White BC albedo={self.albedo} < 0",
                law="white",
            )

    def assert_submarkov(self) -> None:
        r"""White albedo satisfies the sub-Markov bound
        :math:`\alpha \le 1`.

        Raises
        ------
        SubmarkovViolationError
            When ``self.albedo > 1``.
        """
        if self.albedo > 1.0:
            from ._errors import SubmarkovViolationError
            raise SubmarkovViolationError(
                f"White BC albedo={self.albedo} > 1",
                law="white",
            )

    def assert_realizable(
        self,
        quadrature: "Quadrature",
        *,
        inflow_indices: "Optional[np.ndarray]" = None,
    ) -> None:
        r"""Universal invariants + the sub-Markov bound (ERR-046).

        ``assert_submarkov`` is white-specific (not one of the §16A.12
        universal five), so it joins the certification here — the same
        extension pattern as :class:`ReflectiveBoundary`'s table
        checks.
        """
        super().assert_realizable(quadrature, inflow_indices=inflow_indices)
        self.assert_submarkov()
