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
Today: :class:`PartialCurrentOperator` and
:class:`IsotropicEmissionOperator` — the two links of the factored Lambertian
(white / diffuse-albedo, via the SN realizer) — are both SN-only.

``IncomingSourceOperator`` lived here until **P3** of the affine-boundary-source
campaign (2026-08-05). It realized prescribed inflow as a rank-0 map whose
``apply`` ignored its input and returned :math:`q` — i.e. an **affine** map
declared as a :class:`~orpheus.numerics.operator.LinearOperator`, measured
``A(0) = 2.5`` and ``A(2x) − 2A(x) = 2.5`` at :math:`q = 2.5`. It was retired
onto the **zero morphism** every other zero-response law already returns
(:func:`~orpheus.sn.boundary.realizer._narrowed_zero_operator`): the affine law
is :math:`\gamma_-\psi = L\,\gamma_+\psi + q`, this tier realizes :math:`L`, and
for prescribed inflow :math:`L = 0`. The source travels the boundary-source
channel instead (:ref:`bc-affine-source-channel`). Retiring it removed a
DOUBLE delivery, not merely an untidy type — see the realizer's own arm for the
measurement.

``AngularAverageOperator`` lived here until **G6.3 step 3b** (2026-08-04) as the
WELDED Lambertian: one operator :math:`\Gamma_+ \to \Gamma_-` that withheld its
transpose because its Euclidean and cosine-weighted readings differ. It was
retired onto ``IsotropicEmissionOperator @ PartialCurrentOperator``, which
removes that ambiguity rather than choosing between the readings — each link has
one honest transpose — and names the intermediate (the outgoing partial current)
that its body carried as an anonymous local.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    LinearOperator,
    checked_space_extent,
)

if TYPE_CHECKING:
    from orpheus.numerics.space import FunctionSpace

__all__ = [
    "IsotropicEmissionOperator",
    "PartialCurrentOperator",
]


# ─────────────────────────────────────────────────────────────────────
# The Lambertian, FACTORED — G6.3 step 3 (#330)
# ─────────────────────────────────────────────────────────────────────


class PartialCurrentOperator(LinearOperator):
    r""":math:`C : \Gamma_+(f) \to S(f)` — the cosine-weighted partial current.

    .. math::

        (C\psi)_g \;=\; \sum_{\Omega \in \Gamma_+}
                        w(\Omega)\,|\Omega\cdot\hat n|\;\psi_g(\Omega)
                \;=\; J^+_g,

    the **outgoing partial current** — a named physical quantity with units,
    which until G6.3 existed only as the anonymous local ``psi_avg`` inside
    the retired welded ``AngularAverageOperator``'s apply (and there already
    divided by the
    normalisation, so not even a current).

    The first link of the factored Lambertian
    :math:`\Gamma_+(f) \to S(f) \to \Gamma_-(f)`; its partner is
    :class:`IsotropicEmissionOperator`. Splitting here rather than anywhere else
    is what makes the intermediate NAMEABLE: this is precisely where the angular
    axis is integrated out, and the result is the quantity an albedo law scales
    (:math:`J^- = \alpha J^+`, the albedo family in its own native form).

    **Its transpose is a theorem, not a hand-roll.** :math:`C` is the row vector
    :math:`\cos w`, so :math:`C^{\mathsf T}` is the outer product with it — one
    line, no convention to get wrong. That is the whole reason the composite's
    adjoint stops being deferred.

    Parameters
    ----------
    cos_w : np.ndarray
        The outflow cosine weights :math:`w\,|\Omega\cdot\hat n|` restricted to
        :math:`\Gamma_+`, strictly positive since **B3.4a** (the operator is
        HANDED its half-trace and classifies nothing).
    domain, codomain : FunctionSpace, optional
        :math:`\Gamma_+(f)` and :math:`S(f)`, from
        :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_space`
        / :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.current_space`.
    """

    def __init__(
        self,
        cos_w: np.ndarray,
        *,
        domain: "FunctionSpace | None" = None,
        codomain: "FunctionSpace | None" = None,
    ) -> None:
        cos_w = np.asarray(cos_w, dtype=float)
        if cos_w.ndim != 1:
            raise ValueError(
                f"PartialCurrentOperator cos_w must be 1-D; got shape "
                f"{cos_w.shape}"
            )
        if cos_w.size == 0:
            raise ValueError(
                "PartialCurrentOperator cos_w is empty — a face with no "
                "outflow ordinates has no partial current to contract."
            )
        if not np.all(cos_w > 0.0):
            raise ValueError(
                "PartialCurrentOperator cos_w must be STRICTLY positive: its "
                "domain is Γ₊, whose tangential rows are excluded by "
                "construction since B3.4a. A zero entry means the caller "
                "passed a full-face weight vector instead of the outflow "
                f"restriction. Got min={cos_w.min():.3e}."
            )
        self._cos_w = cos_w.copy()
        self._cos_w.flags.writeable = False
        # ⭐ The weights and the domain both state |Γ₊|; a disagreement is
        # SILENT at apply-time (the arrays broadcast) and is the B3.4a twin
        # classifier's last hiding place — a caller classifying the retired
        # `> 0.0` way on BOTH the weights and the input would otherwise get a
        # wrong average with no error anywhere.
        self._domain = checked_space_extent(
            domain, cos_w.size,
            owner="PartialCurrentOperator", role="domain",
        )
        self._codomain = codomain

    @property
    def domain(self) -> "FunctionSpace | None":
        return self._domain

    @property
    def codomain(self) -> "FunctionSpace | None":
        return self._codomain

    @property
    def is_adjointable(self) -> bool:
        """``True`` — the transpose is an outer product, with no ambiguity."""
        return True

    @property
    def cos_w(self) -> np.ndarray:
        """The outflow cosine weights (read-only view)."""
        return self._cos_w

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Contract :math:`\Gamma_+` to the partial current: ``(|Γ₊|, …) → (…)``."""
        psi = np.asarray(psi)
        if psi.shape[0] != self._cos_w.size:
            raise ValueError(
                f"PartialCurrentOperator.apply: psi has {psi.shape[0]} rows "
                f"along the leading axis, expected |Γ₊| = {self._cos_w.size}. "
                f"Its DOMAIN is the outflow half-trace — restrict with the "
                f"trace space's γ₊ first."
            )
        return (
            self._cos_w.reshape((-1,) + (1,) * (psi.ndim - 1)) * psi
        ).sum(axis=0)

    def apply_transpose(self, current: np.ndarray) -> np.ndarray:
        r"""``Cᵀ(s) = cos_w ⊗ s``: ``(…) → (|Γ₊|, …)``.

        The Euclidean transpose. The *Hilbert* adjoint additionally carries the
        two spaces' metrics and is reached through ``.H``; both are honest
        here, which is why this one can be advertised (contrast
        the retired ``AngularAverageOperator``, which withheld its transpose precisely
        because the composite's two readings differed).
        """
        current = np.asarray(current)
        return (
            self._cos_w.reshape((-1,) + (1,) * current.ndim)
            * current[None, ...]
        )


class IsotropicEmissionOperator(LinearOperator):
    r""":math:`B : S(f) \to \Gamma_-(f)` — a current, re-emitted isotropically.

    .. math::

        (B J)_g(\Omega) \;=\; \frac{J_g}{\textstyle\sum_{\Gamma_+}
                              w\,|\Omega'\cdot\hat n|},
        \qquad \Omega \in \Gamma_-,

    independent of :math:`\Omega` — Lambertian / cosine emission. The second
    link of the factored Lambertian; its partner is
    :class:`PartialCurrentOperator`.

    **The normalisation belongs HERE, not in the contraction**, and that is
    physics rather than bookkeeping: :math:`C` produces a *current* and
    :math:`B` must produce an *intensity*, so the division is exactly the
    unit-changing step. Placing it here leaves :math:`S(f)` carrying an honest
    :math:`J^+`, which in turn lets an albedo enter as the pure scalar law
    :math:`J^- = \alpha J^+` rather than as a factor smeared across two
    operators.

    The denominator is the OUTGOING cosine-weighted weight sum, which makes the
    composite conservative. The measured discussion of that choice, and of the
    tangential-band hazard the B3.4a narrowing dissolved, moved to
    :class:`PartialCurrentOperator`'s strictly-positive ``cos_w`` guard when
    the welded ``AngularAverageOperator`` was retired at G6.3 step 3b.

    Parameters
    ----------
    norm : float
        :math:`\sum_{\Gamma_+} w\,|\Omega\cdot\hat n|`. MUST be strictly
        positive.
    n_inflow : int
        :math:`|\Gamma_-|`, the codomain extent. Required rather than inferred:
        the operator maps BETWEEN two half-traces, and the zero of the space it
        lands in is not deducible from the one it came from.
    domain, codomain : FunctionSpace, optional
        :math:`S(f)` and :math:`\Gamma_-(f)`.
    """

    def __init__(
        self,
        norm: float,
        n_inflow: int,
        *,
        domain: "FunctionSpace | None" = None,
        codomain: "FunctionSpace | None" = None,
    ) -> None:
        norm = float(norm)
        if not norm > 0.0:
            raise ValueError(
                f"IsotropicEmissionOperator norm must be strictly positive; "
                f"got {norm:g}. It is the outgoing cosine-weighted weight sum, "
                f"which vanishes only for an empty Γ₊."
            )
        n_inflow = int(n_inflow)
        if n_inflow <= 0:
            raise ValueError(
                f"IsotropicEmissionOperator n_inflow must be positive; got "
                f"{n_inflow}. A face with no inflow ordinates has nothing to "
                f"re-emit into."
            )
        self._norm = norm
        self._n_inflow = n_inflow
        self._domain = domain
        self._codomain = checked_space_extent(
            codomain, n_inflow,
            owner="IsotropicEmissionOperator", role="codomain",
        )

    @property
    def domain(self) -> "FunctionSpace | None":
        return self._domain

    @property
    def codomain(self) -> "FunctionSpace | None":
        return self._codomain

    @property
    def is_adjointable(self) -> bool:
        """``True`` — the transpose is a sum over :math:`\\Gamma_-`."""
        return True

    @property
    def norm(self) -> float:
        """The outgoing cosine-weighted weight sum."""
        return self._norm

    @property
    def n_inflow(self) -> int:
        r""":math:`|\Gamma_-|`."""
        return self._n_inflow

    def apply(self, current: np.ndarray) -> np.ndarray:
        r"""Broadcast the normalised current over :math:`\Gamma_-`: ``(…) → (|Γ₋|, …)``.

        The array is fresh (callers may mutate freely).
        """
        current = np.asarray(current)
        return np.broadcast_to(
            (current / self._norm)[None, ...],
            (self._n_inflow,) + current.shape,
        ).copy()

    def apply_transpose(self, phi: np.ndarray) -> np.ndarray:
        r"""``Bᵀ(φ) = Σ_{Γ₋} φ / norm``: ``(|Γ₋|, …) → (…)``."""
        phi = np.asarray(phi)
        if phi.shape[0] != self._n_inflow:
            raise ValueError(
                f"IsotropicEmissionOperator.apply_transpose: phi has "
                f"{phi.shape[0]} rows along the leading axis, expected "
                f"|Γ₋| = {self._n_inflow}."
            )
        return phi.sum(axis=0) / self._norm
