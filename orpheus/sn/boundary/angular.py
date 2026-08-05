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
(white / diffuse-albedo, via the SN realizer) — and
:class:`IncomingSourceOperator` (prescribed inflow) are all SN-only.

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
from orpheus.numerics.face_layout import AXIS_NAMES, face_name
from orpheus.numerics.spaces.angular_trace_space import (
    TANGENTIAL_EPS,
    build_omega_dot_n,
)

if TYPE_CHECKING:
    from orpheus.geometry.boundary._source import InflowSourceSpec
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.numerics.space import FunctionSpace

__all__ = [
    "IncomingSourceOperator",
    "IsotropicEmissionOperator",
    "PartialCurrentOperator",
]


class IncomingSourceOperator(LinearOperator):
    r"""Prescribed inflow source — returns the source value, ignores input.

    Realises the :math:`q` term in the §16A.1 affine BC form

    .. math::

        \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q

    for the rank-0 case where :math:`R = 0` and only the source
    matters: :math:`\gamma_- \psi = q`. (Campaign phase **B3.0**
    corrected the older spelling ":math:`R = G = 0`" — it is the
    RESPONSE that vanishes; :math:`G` is the identity deck element,
    because the zero map is not a bijection and so cannot be a geometry
    map at all. Writing both as zero spelled one vanishing twice, once
    in the wrong tier.) The :meth:`apply` IGNORES its
    input and returns the source evaluated on a probe inflow trace
    matching the input shape. Used by the SN realizer for
    :class:`~orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow`.

    **B3.4a — the inflow mask dissolved.** An
    :class:`~orpheus.geometry.boundary._source.InflowSourceSpec` fills
    whatever block *shape* it is handed (it carries no trace knowledge),
    so pre-B3.4a this operator emitted a FULL-FACE block and then zeroed
    every non-inflow ordinate — outflow AND tangential — to make the
    affine form's :math:`q \in \Gamma_-` hold. With the codomain narrowed
    to :math:`\Gamma_-` the operator simply asks the spec to fill a
    :math:`|\Gamma_-|`-row block: the rows the mask used to zero are no
    longer in the codomain to be emitted on, so ERR-047 is closed by the
    TYPE rather than by an erasure. That is the same dissolution the
    vacuum projector underwent at B3.2 — the mask was never load-bearing
    physics, only the cost of a codomain that was too big.

    The masked branch had a companion: an *un*masked fallback for a method
    space carrying no inflow indices, legal only for :math:`q \equiv 0`.
    It is retired with the mask. Post-B3.2 every realization needs
    :math:`\gamma_+` too, so a method space with no face data cannot reach
    this operator at all — the fallback was already unreachable on the
    realize path.

    Apply-only. The operator is rank-0 in the input (every input maps to
    the same source value); it is NOT invertible and NOT naturally
    self-adjoint. ``solve`` and ``apply_transpose`` are deliberately NOT
    advertised.

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
    n_inflow : int
        :math:`|\Gamma_-|` — the codomain size, i.e. how many inflow
        ordinates this face carries. Required: the source is asked to
        fill the CODOMAIN's shape, which the domain does not determine.
    """

    def __init__(
        self,
        source: "InflowSourceSpec",
        *,
        n_inflow: int,
    ) -> None:
        self.source = source
        n_inflow = int(n_inflow)
        if n_inflow < 0:
            raise ValueError(
                f"IncomingSourceOperator n_inflow must be non-negative; "
                f"got {n_inflow}."
            )
        #: :math:`|\Gamma_-|` — the codomain size the source is asked to fill.
        self.n_inflow = n_inflow

    def apply(self, psi_out: np.ndarray) -> np.ndarray:
        r"""Return :math:`q` on :math:`\Gamma_-`. ``psi_out`` is IGNORED
        (an affine source does not depend on the outgoing flux) — only its
        TRAILING axes are read, to size the group / spatial block.

        Since **B3.4a** the source is asked to fill
        ``(|Γ₋|,) + psi_out.shape[1:]`` directly, so the affine form's
        :math:`q \in \Gamma_-` holds because there are no other rows to
        write, not because a mask erased them.
        """
        psi_out = np.asarray(psi_out)
        return self.source.evaluate(
            (self.n_inflow,) + tuple(int(s) for s in psi_out.shape[1:])
        )


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
        :class:`AngularAverageOperator`, which withheld its transpose precisely
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
