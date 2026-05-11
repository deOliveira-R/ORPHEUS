r"""SN-specific angular operator primitives for the boundary realizer.

Hosts angular primitives that are SN-specific (i.e. depend on the
:class:`~orpheus.sn.quadrature.AngularQuadrature` direction-cosine /
weight arrays) but are independent of any single boundary condition.
They are consumed by the Wave-5
:class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` to realize
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances as
single-arg :class:`~orpheus.numerics.operator.LinearOperator` objects.

Promotion criterion: when a second non-SN consumer arrives, an
operator here should be promoted to ``orpheus/numerics/`` (per
Cardinal Rule 2 — no shared abstraction with only one consumer).
Today: ``AngularAverageOperator`` is SN-only (the white BC consumes
it via the SN realizer).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import (
    CAP_APPLY,
    LinearOperatorMixin,
)

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature

__all__ = ["AngularAverageOperator"]


class AngularAverageOperator(LinearOperatorMixin):
    r"""Cosine-weighted Lambertian average over an outgoing hemisphere.

    Realizes the geometric projection :math:`G_{\text{diff}}` in the
    §15.2 sum-of-tensor-products form of the white boundary law
    :math:`R_{\text{white}} = G_{\text{diff}} \otimes \alpha`. The
    action on an angular flux array :math:`\psi(\Omega, \ldots)` is

    .. math::

        (G_{\text{diff}}\,\psi)(\Omega, \ldots) \;=\;
        \frac{\sum_{\Omega' :\, \mathrm{sign}\,(\Omega' \cdot \hat n) > 0}
              w(\Omega')\,|\Omega' \cdot \hat n|\,\psi(\Omega', \ldots)}
             {\sum_{\Omega' :\, \mathrm{sign}\,(\Omega' \cdot \hat n) > 0}
              w(\Omega')\,|\Omega' \cdot \hat n|}.

    The result is independent of the output :math:`\Omega` direction
    (Lambertian / cosine emission), so the apply broadcasts the scalar
    angular average back over all ordinates uniformly. The denominator
    is the outgoing cosine-weighted weight sum — this normalization
    makes the operator conservative for any quadrature: the cosine-
    weighted outgoing current equals the cosine-weighted incoming
    current (per Bell & Glasstone 1970 §1.5).

    The operator stores defensive copies of the cosine-weighted
    outgoing mask values and the normalization scalar at construction
    time — no reference to the source quadrature is retained. This
    decouples the operator's lifetime from the quadrature object and
    makes the operator safely picklable / hash-stable.

    Parameters
    ----------
    cos_w : numpy.ndarray
        Shape ``(N_ord,)``. Cosine-weighted ordinate weights, zeroed
        on the non-outgoing hemisphere:
        ``np.where(outgoing_mask, weights * (outward_sign * mu_n), 0.0)``.
        Built once at construction; defensively copied. Always
        non-negative.
    norm : float
        ``cos_w.sum()``. The cosine-weighted outgoing-current
        normalization. Computed once at construction. MUST be
        strictly positive — a zero ``norm`` indicates a degenerate
        quadrature (no outgoing ordinates on the chosen face) and
        the constructor raises ``ValueError``.

    Capabilities
    ------------
    ``{CAP_APPLY}``. The white BC is not self-adjoint in the
    unweighted inner product; it IS self-adjoint under the cosine-
    weighted inner product. ``apply_transpose`` is NOT advertised
    because the unweighted transpose semantics differ from the
    weighted-adjoint that the BC analytically possesses — exposing
    ``apply_transpose`` here would invite the harmful confusion of
    two different ``.T`` semantics. The Hilbert adjoint via ``.H``
    on a weighted ``FunctionSpace`` is the right channel; that's
    deferred until a consumer needs it (file follow-up issue if/when
    it does).

    See Also
    --------
    :class:`~orpheus.geometry.boundary.WhiteBoundaryOperator`
        Legacy 2-arg BC; Wave 7 will rename to ``WhiteBoundary`` and
        delegate to ``AngularAverageOperator`` via the SN realizer.
    """

    capabilities: ClassVar[frozenset[str]] = frozenset({CAP_APPLY})

    def __init__(self, cos_w: np.ndarray, norm: float) -> None:
        cos_w = np.asarray(cos_w, dtype=float)
        if cos_w.ndim != 1:
            raise ValueError(
                f"AngularAverageOperator cos_w must be 1-D; "
                f"got shape {cos_w.shape}"
            )
        if (cos_w < 0).any():
            raise ValueError(
                "AngularAverageOperator cos_w must be non-negative "
                "(zeroed on non-outgoing ordinates)"
            )
        norm = float(norm)
        if norm <= 0.0:
            raise ValueError(
                f"AngularAverageOperator norm must be positive; got "
                f"{norm}. A zero norm indicates a degenerate quadrature "
                f"(no outgoing ordinates on the chosen face)."
            )
        # Defensive copy — no reference to caller's array retained.
        self._cos_w = cos_w.copy()
        self._norm = norm
        self.n_ordinates = cos_w.size

    @classmethod
    def from_quadrature(
        cls,
        quadrature: "AngularQuadrature",
        axis: str,
        outward_sign: int,
    ) -> "AngularAverageOperator":
        r"""Construct from a quadrature + face orientation.

        Parameters
        ----------
        quadrature : AngularQuadrature
            Quadrature whose ``weights`` and ``mu_x`` / ``mu_y`` /
            ``mu_z`` arrays drive the cosine-weighted average. Only
            slices of those arrays are retained (copied) — the
            quadrature reference is NOT held by the operator.
        axis : {"x", "y", "z"}
            Boundary-normal axis. ``"z"`` requires a quadrature with
            ``mu_z`` (Lebedev, level-symmetric, product). The 1-D
            Gauss-Legendre adapter has no ``mu_z`` and raises
            ``ValueError`` if ``axis="z"``.
        outward_sign : {+1, -1}
            Sign of the outward normal along ``axis``. ``+1`` for
            upper faces (``xmax``, ``ymax``); ``-1`` for lower faces
            (``xmin``, ``ymin``). Selects which ordinates are
            *outgoing* at this face.

        Returns
        -------
        AngularAverageOperator
            Self-contained operator. The quadrature is not retained.
        """
        if axis == "x":
            mu_n = quadrature.mu_x
        elif axis == "y":
            mu_n = quadrature.mu_y
        elif axis == "z":
            mu_n = getattr(quadrature, "mu_z", None)
            if mu_n is None:
                raise ValueError(
                    "AngularAverageOperator.from_quadrature(axis='z') "
                    "requires a quadrature with mu_z (2-D / 3-D "
                    "adapters: Lebedev, level-symmetric, product). The "
                    "1-D Gauss-Legendre adapter has no mu_z attribute."
                )
        else:
            raise ValueError(f"Unknown axis: {axis!r}")

        if outward_sign not in (+1, -1):
            raise ValueError(
                f"outward_sign must be +1 or -1; got {outward_sign}"
            )

        weights = quadrature.weights
        # Outgoing ordinates: positive direction-cosine along the outward normal.
        outgoing_mask = (outward_sign * mu_n) > 0.0
        cos_w = weights * (outward_sign * mu_n)
        # Zero contributions from non-outgoing ordinates.
        cos_w = np.where(outgoing_mask, cos_w, 0.0)
        norm = float(cos_w.sum())
        if norm <= 0.0:
            raise ValueError(
                f"AngularAverageOperator.from_quadrature: no outgoing "
                f"ordinates on axis={axis!r}, outward_sign={outward_sign}. "
                f"Quadrature is degenerate for this face."
            )
        return cls(cos_w=cos_w, norm=norm)

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Cosine-weighted average over the outgoing hemisphere.

        Parameters
        ----------
        psi : numpy.ndarray
            Angular flux. Shape ``(N_ord, ...)`` where the leading
            axis is the ordinate axis and trailing axes are spatial /
            group / etc. The operator broadcasts the cosine-weighted
            average over all leading-axis entries (Lambertian
            emission), so the output shape equals ``psi.shape``.

        Returns
        -------
        numpy.ndarray
            Same shape as ``psi``. The output's leading-axis slices
            are all equal (the broadcast average) and the output is
            a fresh array (calling code may mutate freely).
        """
        if psi.shape[0] != self.n_ordinates:
            raise ValueError(
                f"AngularAverageOperator.apply: psi.shape[0] = "
                f"{psi.shape[0]}, expected {self.n_ordinates}"
            )
        # Cosine-weighted scalar average over the outgoing hemisphere.
        psi_avg = (
            self._cos_w.reshape((-1,) + (1,) * (psi.ndim - 1)) * psi
        ).sum(axis=0) / self._norm
        # Broadcast the (..., ) average back over every ordinate.
        return np.broadcast_to(
            psi_avg[None, ...], psi.shape,
        ).copy()
