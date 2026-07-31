r"""SN-specific angular operator primitives for the boundary realizer.

Hosts angular primitives that are SN-specific (i.e. depend on the
:class:`~orpheus.numerics.quadrature.Quadrature` direction-cosine /
weight arrays) but are independent of any single boundary condition.
They are consumed by the Wave-5
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` to realize
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances as
single-arg :class:`~orpheus.numerics.operator.LinearOperator` objects.

Promotion criterion: when a second non-SN consumer arrives, an
operator here should be promoted to ``orpheus/numerics/`` (per
Cardinal Rule 2 — no shared abstraction with only one consumer).
Today: ``AngularAverageOperator`` is SN-only (the white BC consumes
it via the SN realizer); :class:`IncomingSourceOperator` is SN-only
(the prescribed-inflow BC consumes it via the SN realizer).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import (
    LinearOperator,
)

if TYPE_CHECKING:
    from orpheus.geometry.boundary._source import InflowSourceSpec
    from orpheus.numerics.quadrature import Quadrature

__all__ = ["AngularAverageOperator", "IncomingSourceOperator"]


class AngularAverageOperator(LinearOperator):
    r"""Cosine-weighted Lambertian average over an outgoing hemisphere.

    Realizes the geometric projection :math:`G_{\text{diff}}` in the
    §15.2 sum-of-tensor-products form of the white boundary law
    :math:`B_{\text{white}} = G_{\text{diff}} \otimes \alpha`.

    .. note::

       :math:`B` is the **realized composite**, matching the sibling
       reflective arm's ``B = G_patch ⊗ K_omega ⊗ K_g``
       (:mod:`orpheus.sn.boundary.realizer`). This spelling was
       :math:`R_{\text{white}}` until B1, which collided with the affine
       form's :math:`R` — the **response factor alone**, the crossing
       :math:`\Gamma_+ \to \Gamma_-` that
       :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.response_kernel`
       now carries as a type. Here that factor is the scalar
       :math:`\alpha`; :math:`G_{\text{diff}}` is the geometry.

    The action on an angular flux array :math:`\psi(\Omega, \ldots)` is

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
    apply-only (structurally non-invertible, non-adjointable as
    posed). The white BC is not self-adjoint in the
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
    :class:`~orpheus.geometry.boundary.WhiteBoundary`
        The white-reflection BC; delegates to
        ``AngularAverageOperator`` via the SN realizer.
    """

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
        quadrature: "Quadrature",
        axis: str,
        outward_sign: int,
    ) -> "AngularAverageOperator":
        r"""Construct from a quadrature + face orientation.

        Parameters
        ----------
        quadrature : Quadrature
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


class IncomingSourceOperator(LinearOperator):
    r"""Prescribed inflow source — returns the source value, ignores input.

    Realises the :math:`q` term in the §16A.1 affine BC form

    .. math::

        \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

    for the rank-0 case where :math:`R = G = 0` and only the source
    matters: :math:`\gamma_- \psi = q`. The :meth:`apply` IGNORES its
    input and returns the source evaluated on a probe inflow trace
    matching the input shape. Used by the SN realizer for
    :class:`~orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow`.

    The operator stores the
    :class:`~orpheus.geometry.boundary._source.InflowSourceSpec` and —
    when the realizing method space carries them (#52, ERR-047) — the
    face's inflow-ordinate indices. An
    :class:`~orpheus.geometry.boundary._source.InflowSourceSpec` fills
    whatever block *shape* it is handed (it carries no trace
    knowledge), so restricting the delivered :math:`q` to
    :math:`\Gamma_-` is THIS operator's contract: with a mask, the
    evaluation is zeroed on every non-inflow ordinate (outflow AND
    tangential — neither is in :math:`\Gamma_-`), making the affine
    form's :math:`q \in \Gamma_-` hold by construction. Without a
    mask, the raw evaluation passes through — legal only for
    :math:`q \equiv 0` sources on the realize path, which
    :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
    certifies BEFORE this operator is built.

    Apply-only. The operator is rank-0 in the
    input (every input maps to the same source value); it is NOT
    invertible and NOT naturally self-adjoint. ``solve`` and
    ``apply_transpose`` are deliberately NOT advertised.

    Parameters
    ----------
    source : InflowSourceSpec
        The :math:`q` source generator. Typically
        :class:`~orpheus.geometry.boundary._source.NoSource` (no
        inflow, equivalent to vacuum) or
        :class:`~orpheus.geometry.boundary._source.ConstantInflowSource`
        for a uniform inflow level. Custom
        :class:`~orpheus.geometry.boundary._source.InflowSourceSpec`
        implementations may inject spatially / energy- / angularly-
        varying inflow.
    inflow_indices : np.ndarray, optional
        The face's inflow-ordinate indices (axis 0 of the boundary
        block). When given, ``n_ordinates`` is required and the
        delivered source is masked to these slots.
    n_ordinates : int, optional
        Length of the ordinate axis; validates the mask range
        (mirrors :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`'s
        construction guards) and pins the expected ``apply`` shape.
    """

    def __init__(
        self,
        source: "InflowSourceSpec",
        *,
        inflow_indices: "np.ndarray | None" = None,
        n_ordinates: "int | None" = None,
    ) -> None:
        self.source = source
        if inflow_indices is None:
            self._inflow_mask: "np.ndarray | None" = None
            self.n_ordinates: "int | None" = (
                int(n_ordinates) if n_ordinates is not None else None
            )
            return
        if n_ordinates is None:
            raise ValueError(
                "IncomingSourceOperator: inflow_indices requires "
                "n_ordinates (the mask needs the ordinate-axis length)."
            )
        indices = np.asarray(inflow_indices, dtype=np.intp)
        if indices.ndim != 1:
            raise ValueError(
                f"IncomingSourceOperator inflow_indices must be 1-D; "
                f"got shape {indices.shape}"
            )
        if indices.size > 0 and (
            indices.min() < 0 or indices.max() >= int(n_ordinates)
        ):
            raise ValueError(
                f"IncomingSourceOperator inflow_indices out of range "
                f"[0, {n_ordinates}); got min={int(indices.min())}, "
                f"max={int(indices.max())}"
            )
        mask = np.zeros(int(n_ordinates), dtype=float)
        mask[indices] = 1.0
        self._inflow_mask = mask
        self.n_ordinates = int(n_ordinates)

    def apply(self, psi_out: np.ndarray) -> np.ndarray:
        r"""Return the source evaluated at ``psi_out.shape``, masked to
        :math:`\Gamma_-` when inflow indices were wired. ``psi_out``
        itself is IGNORED (the affine source does not depend on the
        outgoing flux).

        The source is asked to fill an array of the incoming-face
        shape; the mask (ordinate axis 0, broadcast over trailing
        group / spatial axes) then zeroes every non-inflow slot so the
        delivered :math:`q` lives on the incoming trace by
        construction (ERR-047).
        """
        q = self.source.evaluate(tuple(int(s) for s in psi_out.shape))
        if self._inflow_mask is None:
            return q
        if q.shape[0] != self._inflow_mask.shape[0]:
            raise ValueError(
                f"IncomingSourceOperator: ordinate axis mismatch — the "
                f"inflow mask covers {self._inflow_mask.shape[0]} "
                f"ordinates but the requested block has axis-0 length "
                f"{q.shape[0]}."
            )
        return q * self._inflow_mask.reshape(
            (-1,) + (1,) * (q.ndim - 1)
        )
