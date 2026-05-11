r"""White (Lambertian) boundary condition.

See :class:`WhiteBoundary` for the algebraic definition. The legacy
``WhiteBoundaryOperator`` name is re-exported as a deprecated alias
from the package ``__init__.py`` (Wave 7 rename per Grand Report v3
vocabulary).
"""

from __future__ import annotations

from dataclasses import dataclass

from ._base import BoundaryTraceLaw


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
    :class:`~orpheus.sn.angular_operator.AngularAverageOperator`
    (α=1 fast path) or ``ScaledOperator(α, AngularAverageOperator)``
    (α ≠ 1). Realise via:

    .. code-block:: python

        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer, SNMethodSpace,
        )
        law = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.8)
        op = SNBoundaryRealizer().realize(law, SNMethodSpace.minimal(quad))
        psi_in = op.apply(psi_out)

    Wave-7 rename note
    ------------------
    Previously named ``WhiteBoundaryOperator``. The legacy name is
    preserved as a deprecated alias.

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
