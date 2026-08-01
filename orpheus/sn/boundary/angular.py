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

from typing import TYPE_CHECKING

import numpy as np

from orpheus.numerics.operator import (
    LinearOperator,
)
from orpheus.numerics.face_layout import AXIS_NAMES
from orpheus.numerics.spaces.angular_trace_space import (
    TANGENTIAL_EPS,
    build_omega_dot_n,
)

if TYPE_CHECKING:
    from orpheus.geometry.boundary._source import InflowSourceSpec
    from orpheus.numerics.quadrature import Quadrature

__all__ = ["AngularAverageOperator", "IncomingSourceOperator"]


class AngularAverageOperator(LinearOperator):
    r"""Cosine-weighted Lambertian re-emission, :math:`\Gamma_+ \to \Gamma_-`.

    Realizes the **response** factor :math:`R_{\text{diff}}` of the white
    boundary law — the spec
    :class:`~orpheus.geometry.boundary.LambertianReemission`. White's
    geometry factor is :class:`~orpheus.geometry.boundary.IdentityMap`: a
    white face fixes no geometry, and all of its content is constitutive.

    **Not a geometry, despite the older name.** Campaign phase **B3.0**
    moved this out of the geometry tier: a geometry map is the composition
    operator of a measure-preserving bijection, hence **multiplicative**
    (:math:`G(\psi\varphi) = (G\psi)(G\varphi)`), which an average never
    is. The misassignment survived because a rank-one response annihilates
    :math:`G` entirely (:ref:`bc-factor-roles`), so the slot it occupied
    had no observable consequence.

    The action, on the **outflow half-trace** :math:`\Gamma_+`, is

    .. math::

        (R_{\text{diff}}\,\psi)(\Omega) \;=\;
        \frac{\sum_{\Omega' \in \Gamma_+}
              w(\Omega')\,|\Omega' \cdot \hat n|\,\psi(\Omega')}
             {\sum_{\Omega' \in \Gamma_+}
              w(\Omega')\,|\Omega' \cdot \hat n|},
        \qquad \Omega \in \Gamma_-.

    The result is independent of the output :math:`\Omega` (Lambertian /
    cosine emission), so the apply broadcasts one scalar over every
    **inflow** ordinate. The denominator is the outgoing cosine-weighted
    weight sum — the normalization that makes the operator conservative
    for any quadrature: the cosine-weighted outgoing current equals the
    cosine-weighted incoming current (Bell & Glasstone 1970 §1.5).

    Narrowed at **B3.4a — and the narrowing dissolved a real bug.**
    Pre-B3.4a this was a full-face endomorphism whose ``cos_w`` carried ``N``
    entries zeroed off the outgoing hemisphere by its OWN test,
    ``(outward_sign * mu_n) > 0.0`` — a *second* outflow classifier, which
    disagrees with the trace space's ``TANGENTIAL_EPS`` one wherever an
    ordinate sits in the tangential band on the positive side of exact zero.
    With the domain narrowed the operator no longer classifies anything — it
    is *handed* :math:`\Gamma_+` — so the twin is not fixed but
    **unspellable**, and every ``cos_w`` entry is strictly positive by
    construction rather than mostly zero.

    Measured 2026-08-01, over every production quadrature × face:

    .. list-table::
       :header-rows: 1
       :widths: 34 22 44

       * - quadrature
         - tangential / face
         - MIS-ADMITTED by ``> 0.0``
       * - ``gauss_legendre(8)``
         - 0
         - 0
       * - ``product(2,4)`` / ``product(3,4)``
         - 4 / 6
         - **2 / 3**, on ``xmin`` ``xmax`` ``ymax`` — *none* on ``ymin``
       * - ``level_symmetric(6)``
         - 0
         - 0
       * - ``lebedev(17)``
         - 12
         - 0

    So the defect is the **``product`` family only**, and there only on three
    of four faces — the ``ymin`` asymmetry is floating-point sign noise about
    exact zero, so it is not even stable under a face flip. ``lebedev`` carries
    the most tangential ordinates of any production quadrature and mis-admits
    none of them.

    .. warning::

       **Do not justify this narrowing with a ULP table.** On ``product(2,4)``
       ``xmax`` the mis-admitted rows carry ``cos_w ≈ 7.85e-17`` against a norm
       of ``2.5650996603237282`` — below its ULP — so ``Δnorm`` is **exactly
       0.0** and the entire discrepancy sits ψ-weighted in the NUMERATOR. It
       therefore scales with the flux ratio between the tangential rows and
       :math:`\Gamma_+` and is **unbounded by floating point**: measured
       5.9e-14 at a 1e3 ratio, 6.0e-11 at 1e6, and **6.2e-05 at 1e12** — the
       last reproducing the 6.1e-05 recorded independently in
       :class:`~orpheus.geometry.boundary.LambertianReemission` from a cylinder
       under ``product(2, 4)``.

       Two mechanically DIFFERENT effects therefore both measure ≤ 1 ULP on an
       O(1) probe: a bounded reduction-order change (``lebedev``,
       ``level_symmetric`` — they mis-admit nothing) and this unbounded
       mis-admission. A magnitude table cannot tell them apart. The
       justification is **structural — one classifier, not two**; the numbers
       are evidence of what the second classifier was doing, never the
       argument.

    Narrowed, this operator IS the rank-one kernel
    :math:`u \otimes v` with :math:`u = \mathbf{1}_{\Gamma_-}` and
    :math:`v = \cos\!w/\!\operatorname{norm}` on :math:`\Gamma_+`. Phase
    **B5** types it as such (which is what makes its adjoint structurally
    available); B3.4a stops short of that, deliberately — the
    Euclidean-vs-cosine-metric transpose question below is B5's to answer.

    The operator stores defensive copies of the cosine weights and the
    normalization scalar — no reference to the source quadrature is
    retained, so its lifetime is decoupled from the quadrature object and
    it stays picklable / hash-stable.

    Parameters
    ----------
    cos_w : numpy.ndarray
        Shape ``(|Γ₊|,)`` — the cosine-weighted ordinate weights **on the
        outflow half-trace**, ``weights[Γ₊] * |Ω·n̂|[Γ₊]``. Every entry is
        strictly positive: :math:`\Gamma_+` is exactly where
        :math:`\Omega\cdot\hat n > +\varepsilon_{\text{tan}}`, so a
        non-positive entry means a non-outflow ordinate leaked into the
        domain and the constructor refuses it.
    norm : float
        ``cos_w.sum()`` — the cosine-weighted outgoing-current
        normalization. MUST be strictly positive.
    n_inflow : int
        :math:`|\Gamma_-|`, the codomain size. The average is broadcast
        over exactly this many rows. Required: the operator maps BETWEEN
        two half-traces, and the zero of the space it lands in is not
        deducible from the space it came from.

    Capabilities
    ------------
    apply-only (structurally non-invertible, non-adjointable as posed).
    The white BC is not self-adjoint in the unweighted inner product; it
    IS self-adjoint under the cosine-weighted one. ``apply_transpose`` is
    NOT advertised because the unweighted transpose semantics differ from
    the weighted adjoint the BC analytically possesses — exposing it here
    would invite two different ``.T`` semantics. The Hilbert adjoint via
    ``.H`` on a weighted ``FunctionSpace`` is the right channel.

    See Also
    --------
    :class:`~orpheus.geometry.boundary.WhiteBoundary`
        The white-reflection BC; delegates here via the SN realizer.
    :class:`~orpheus.geometry.boundary.LambertianReemission`
        The response-kernel SPEC this operator realizes.
    """

    def __init__(
        self, cos_w: np.ndarray, norm: float, n_inflow: int,
    ) -> None:
        cos_w = np.asarray(cos_w, dtype=float)
        if cos_w.ndim != 1:
            raise ValueError(
                f"AngularAverageOperator cos_w must be 1-D; "
                f"got shape {cos_w.shape}"
            )
        # STRICTLY positive since B3.4a, where the old ``>= 0`` was right:
        # ``cos_w`` used to be a full-face array zeroed off the outgoing
        # hemisphere, so zeros were the normal case. Now it is the
        # RESTRICTION to Γ₊, where |Ω·n̂| > TANGENTIAL_EPS and the weights
        # are positive — so a non-positive entry means a tangential or
        # inflow ordinate leaked into the domain. That is the failure the
        # narrowing exists to make unspellable; refuse it loudly rather
        # than average over a row that contributes nothing.
        if (cos_w <= 0).any():
            raise ValueError(
                f"AngularAverageOperator cos_w must be strictly positive "
                f"on Γ₊ (the outflow half-trace, where Ω·n̂ > "
                f"TANGENTIAL_EPS); got {int((cos_w <= 0).sum())} "
                f"non-positive of {cos_w.size}. A zero entry means a "
                f"non-outflow ordinate is in the domain — pass the "
                f"RESTRICTION to Γ₊, not a masked full-face array."
            )
        norm = float(norm)
        if norm <= 0.0:
            raise ValueError(
                f"AngularAverageOperator norm must be positive; got "
                f"{norm}. A zero norm indicates a degenerate quadrature "
                f"(no outgoing ordinates on the chosen face)."
            )
        n_inflow = int(n_inflow)
        if n_inflow <= 0:
            raise ValueError(
                f"AngularAverageOperator n_inflow must be positive; got "
                f"{n_inflow}. A face with no inflow ordinates has no Γ₋ "
                f"for the re-emitted flux to land on."
            )
        # Defensive copy — no reference to caller's array retained.
        self._cos_w = cos_w.copy()
        self._norm = norm
        #: :math:`|\Gamma_+|` — the DOMAIN size (outflow ordinates read).
        self.n_outflow = int(cos_w.size)
        #: :math:`|\Gamma_-|` — the CODOMAIN size (inflow ordinates written).
        self.n_inflow = n_inflow

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

        Notes
        -----
        Since **B3.4a** the two half-traces come from
        :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`
        — the codebase's SINGLE face-name → signed-projection primitive —
        classified against ``TANGENTIAL_EPS``, exactly as the trace space
        does. Pre-B3.4a this method carried its own ``> 0.0`` test, which
        is the producer-side twin the narrowing retires: ``(axis,
        outward_sign)`` and a face NAME are the same datum, so there is no
        reason for two classifiers, and the strict compare is wrong
        wherever a quadrature carries tangential ordinates (see the class
        docstring's measurement).
        """
        if axis not in AXIS_NAMES:
            raise ValueError(
                f"Unknown axis: {axis!r}; expected one of {AXIS_NAMES}"
            )
        if outward_sign not in (+1, -1):
            raise ValueError(
                f"outward_sign must be +1 or -1; got {outward_sign}"
            )
        # (axis, outward_sign) IS a face normal — render it as the face
        # NAME the one primitive speaks, rather than re-deriving mu_axis
        # and its sign here. A quadrature with no genuine cosines on this
        # axis is refused there, naming the rank mismatch between the face
        # and the cubature.
        face = f"{axis}{'max' if outward_sign == +1 else 'min'}"
        omega_dot_n = build_omega_dot_n(quadrature, (face,))[0]
        outflow = np.flatnonzero(omega_dot_n > +TANGENTIAL_EPS)
        inflow = np.flatnonzero(omega_dot_n < -TANGENTIAL_EPS)
        if outflow.size == 0:
            raise ValueError(
                f"AngularAverageOperator.from_quadrature: no outgoing "
                f"ordinates on axis={axis!r}, outward_sign={outward_sign}. "
                f"Quadrature is degenerate for this face."
            )
        if inflow.size == 0:
            raise ValueError(
                f"AngularAverageOperator.from_quadrature: no incoming "
                f"ordinates on axis={axis!r}, outward_sign={outward_sign} "
                f"— the re-emitted flux has no Γ₋ to land on. Quadrature "
                f"is degenerate for this face."
            )
        # No mask, no ``where``: every retained entry is strictly positive
        # because Γ₊ IS the set where Ω·n̂ exceeds the tangential band.
        cos_w = np.asarray(quadrature.weights, dtype=float)[outflow] * (
            omega_dot_n[outflow]
        )
        return cls(
            cos_w=cos_w, norm=float(cos_w.sum()), n_inflow=int(inflow.size),
        )

    def apply(self, psi: np.ndarray) -> np.ndarray:
        r"""Cosine-weighted average of :math:`\Gamma_+`, re-emitted on
        :math:`\Gamma_-`.

        Parameters
        ----------
        psi : numpy.ndarray
            The **outflow** half-trace. Shape ``(|Γ₊|, ...)`` where the
            leading axis is the restricted ordinate axis and trailing axes
            are group / spatial / etc.

        Returns
        -------
        numpy.ndarray
            Shape ``(|Γ₋|, ...)`` — the leading axis is the INFLOW
            half-trace. Its slices are all equal (the broadcast Lambertian
            average) and the array is fresh (callers may mutate freely).
        """
        psi = np.asarray(psi)
        if psi.shape[0] != self.n_outflow:
            raise ValueError(
                f"AngularAverageOperator.apply: psi.shape[0] = "
                f"{psi.shape[0]}, expected |Γ₊| = {self.n_outflow}. Since "
                f"B3.4a this operator's DOMAIN is the outflow half-trace, "
                f"not the whole face slot — restrict with the trace "
                f"space's γ₊ before applying."
            )
        # Cosine-weighted scalar average over the outflow half-trace.
        psi_avg = (
            self._cos_w.reshape((-1,) + (1,) * (psi.ndim - 1)) * psi
        ).sum(axis=0) / self._norm
        # Broadcast that one number over every INFLOW ordinate — the
        # codomain, which is a different space from the domain even when
        # (as on every quadrature in the tree) the two happen to have
        # equal size.
        return np.broadcast_to(
            psi_avg[None, ...], (self.n_inflow,) + psi.shape[1:],
        ).copy()


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
