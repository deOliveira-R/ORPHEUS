r"""White (Lambertian) boundary condition.

See :class:`WhiteBoundary` for the algebraic definition. The legacy
``WhiteBoundaryOperator`` name is re-exported as a deprecated alias
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

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    axis: str = "x"
    outward_sign: int = +1
    albedo: float = 1.0

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        # Direction cosine along the boundary normal axis.
        if self.axis == "x":
            mu_n = quadrature.mu_x
        elif self.axis == "y":
            mu_n = quadrature.mu_y
        elif self.axis == "z":
            mu_n = getattr(quadrature, "mu_z", None)
            if mu_n is None:
                raise ValueError(
                    "WhiteBoundary(axis='z') requires a quadrature with mu_z "
                    "(2-D / 3-D adapters: Lebedev, level-symmetric, "
                    "product). The 1-D Gauss-Legendre adapter has no "
                    "mu_z attribute."
                )
        else:
            raise ValueError(f"Unknown axis: {self.axis!r}")

        weights = quadrature.weights
        # Outgoing ordinates at this face: those whose direction cosine
        # along the outward normal is positive.
        outgoing_mask = (self.outward_sign * mu_n) > 0.0
        cos_w = weights * (self.outward_sign * mu_n)
        # Cosine-weighted outgoing-current normalisation. ``np.where``
        # zeroes contributions from non-outgoing ordinates.
        cos_w = np.where(outgoing_mask, cos_w, 0.0)

        norm = cos_w.sum()
        if norm <= 0.0:
            # Degenerate quadrature -- no outgoing ordinates. Return
            # zero rather than producing a NaN.
            return np.zeros_like(psi_out)

        # Cosine-weighted average of the outgoing flux.
        # Shape (N_ord, ...) -> broadcast (..., ) average.
        psi_avg = (
            cos_w.reshape((-1,) + (1,) * (psi_out.ndim - 1))
            * psi_out
        ).sum(axis=0) / norm

        # Broadcast over all ordinates; sweeps consume only entries
        # whose direction is *incoming* at the face, but it is cheap
        # and conventional to fill the whole array uniformly.
        result = np.broadcast_to(
            psi_avg[None, ...] * self.albedo,
            psi_out.shape,
        ).copy()
        return result
