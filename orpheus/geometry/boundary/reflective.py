r"""Specular (reflective) boundary condition.

See :class:`ReflectiveBoundary` for the algebraic definition. This
class was previously named ``SpecularBoundaryOperator``; the legacy
alias was retired in Wave O step O.4a.1. ``ReflectiveBoundary`` (the
Grand Report v3 vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Optional

import numpy as np

from orpheus.numerics.spaces.angular_trace_space import TANGENTIAL_EPS

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.numerics.quadrature import Quadrature


__all__ = ["ReflectiveBoundary"]


#: Axis tag → ``Quadrature.axis_cosines`` column. The law's ``axis``
#: field speaks the SN slab tags; the quadrature's canonical per-axis
#: accessor speaks column indices.
_AXIS_INDEX: dict[str, int] = {"x": 0, "y": 1, "z": 2}

#: Relative tolerance for the ERR-042 measure-preservation check.
#: Reflection partners are constructed from the same symmetric node/
#: weight sets, so a correct table agrees to a few ULP; a wrong table
#: (mispaired weight classes) violates at O(1). Twelve orders of
#: headroom on both sides.
_MEASURE_RTOL: float = 1e-12


@dataclass(frozen=True)
class ReflectiveBoundary(BoundaryTraceLaw, key="reflective"):
    r"""Specular reflection with optional albedo.

    Tensor decomposition :math:`(G_{\text{refl}}, \alpha)` where
    :math:`G_{\text{refl}}` is the index permutation realising
    :math:`\Omega \mapsto \Omega - 2(\Omega \cdot \hat{n}) \hat{n}` on
    the quadrature ordinates and :math:`\alpha \in [0, 1]` is the
    specular albedo (1 = perfect reflection; the standard ``BC.reflective``
    case). The same operator is the
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward` of
    the angular measure under the reflection map, with the Jacobian
    convention ``|R| = 1`` since reflections are isometries.

    This is a **pure descriptor** (Issue #186 / B3 + β2) — it carries
    no ``apply`` / ``apply_transpose`` methods. Realise via
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` to obtain
    the 1-arg :class:`~orpheus.numerics.operator.PermutationOperator`
    (α=1 fast path) or ``ScaledOperator(α, PermutationOperator)``
    (α ≠ 1):

    .. code-block:: python

        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        law = ReflectiveBoundary(axis="x", albedo=0.7)
        op = SNBoundaryRealizer().realize(
            law, SNMethodSpace.minimal(quad),
        )
        psi_in = op.apply(psi_out)        # forward
        # The realised operator is adjointable as well (a working
        # apply_transpose), consumed by the sensitivity-analysis adjoint
        # pipeline:
        phi_out = op.apply_transpose(phi_in)

    For axis reflections, the index permutation is its own inverse:
    applying the reflection twice returns each ordinate to itself.
    This makes :math:`G_{\text{refl}}^T = G_{\text{refl}}` and so the
    transpose action is identical to the forward action.

    Rename history
    --------------
    Previously named ``SpecularBoundaryOperator``. The legacy alias
    was retired in Wave O step O.4a.1 — ``ReflectiveBoundary`` is the
    sole importable name.

    Parameters
    ----------
    axis : str
        Axis of reflection: ``"x"``, ``"y"``, or ``"z"``. The
        :meth:`~orpheus.numerics.quadrature.Quadrature.reflection_index`
        method maps each ordinate to its reflected partner under this
        axis.
    albedo : float
        Specular albedo. Defaults to 1 (perfect reflection).
    """

    axis: str = "x"
    albedo: float = 1.0

    #: String tag for legacy string-kind comparisons. The default
    #: ``albedo == 1.0`` ReflectiveBoundary (the standard
    #: ``BC.reflective`` case) compares equal to the string
    #: ``"reflective"``; tagged ReflectiveBoundary instances
    #: with ``albedo != 1`` compare equal to ``"partial"`` instead.
    @property
    def kind(self) -> str:
        return "reflective" if self.albedo == 1.0 else "partial"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, ReflectiveBoundary):
            return (
                self.axis == other.axis and self.albedo == other.albedo
            )
        return NotImplemented

    def __hash__(self) -> int:
        # Hash on the canonical (post-rename) class name.
        return hash(("ReflectiveBoundary", self.axis, self.albedo))

    # ------------------------------------------------------------------
    # §16A.12 universal invariants — Wave 7 / C7.6 overrides.
    # ------------------------------------------------------------------

    def assert_is_involutive(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""Reflection permutation must be an involution
        (:math:`\pi \circ \pi = \mathrm{id}`).

        Raises
        ------
        ReflectionNotInvolutiveError
            When ``reflection_index(self.axis)`` is not its own
            inverse on the quadrature's ordinate set.
        """
        ref = quadrature.reflection_index(self.axis)
        if not np.array_equal(ref[ref], np.arange(quadrature.N)):
            from ._errors import ReflectionNotInvolutiveError
            raise ReflectionNotInvolutiveError(
                f"reflection_index({self.axis!r}) is not an involution",
                law="reflective",
            )

    def assert_geometry_map_measure_preserving(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""The reflection table preserves the direction-cosine measure
        :math:`w(\Omega)\,|\Omega\cdot\hat n|` (ERR-042).

        Specular reflection is an isometry, so the pushforward of the
        discrete face measure under the index map :math:`\pi` must be
        the measure itself: with :math:`m_n = w_n\,|\mu_{a,n}|` on the
        reflection axis :math:`a`,

        .. math::

            m_{\pi(n)} \;=\; m_n \qquad \forall\, n .

        Violation means the reflected current at the face differs from
        the incident current — a persistent phantom net leakage at a
        purely reflective boundary (amplified by heterogeneous
        :math:`\Sigma_t` contrast).

        This is checked DIRECTLY, independent of the involution
        property (:meth:`assert_is_involutive`, ERR-044): an involutive
        table that pairs ordinates from different weight classes passes
        the involution check while breaking the measure — exactly the
        hole the pre-#52 delegation ("weight equality is implied by
        construction") left open. The catalog's lesson: the inflow
        partition, the involution, and the measure are INDEPENDENT
        invariants; all must hold.

        Raises
        ------
        BoundaryGeometryMapNotMeasurePreservingError
            When any ordinate and its reflection partner carry
            different direction-cosine measures (relative tolerance
            :data:`_MEASURE_RTOL`; tangential ordinates carry zero
            measure, so a tangential↔non-tangential mispairing is
            caught by the same comparison).
        """
        perm = quadrature.reflection_index(self.axis)
        mu_axis = quadrature.axis_cosines(_AXIS_INDEX[self.axis])
        cosine_measure = quadrature.weights * np.abs(mu_axis)
        partner_measure = cosine_measure[perm]
        if not np.allclose(
            partner_measure, cosine_measure, rtol=_MEASURE_RTOL, atol=0.0
        ):
            from ._errors import (
                BoundaryGeometryMapNotMeasurePreservingError,
            )

            worst = int(
                np.argmax(np.abs(partner_measure - cosine_measure))
            )
            raise BoundaryGeometryMapNotMeasurePreservingError(
                f"reflection_index({self.axis!r}) does not preserve the "
                f"direction-cosine measure w·|μ_{self.axis}|: ordinate "
                f"{worst} carries m={cosine_measure[worst]:.6e} but its "
                f"partner {int(perm[worst])} carries "
                f"m={partner_measure[worst]:.6e} (ERR-042: wrong "
                f"reflection-index table, or a quadrature whose weights "
                f"are inconsistent with its nodes).",
                law="reflective",
            )

    def assert_reflection_maps_inflow_to_outflow(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""Every non-tangential ordinate reflects to the OPPOSITE
        sign class on the reflection axis (ERR-045).

        Specular reflection's geometric contract at a face with normal
        along axis :math:`a`: every inflow ordinate has an outflow
        partner with :math:`\mu_{a,\pi(n)} = -\mu_{a,n}` (and vice
        versa — the sign classes swap, so the check is face-side-free).
        A table that maps an inflow ordinate to itself (or to any
        same-sign partner) plants a length-1 self-loop in the sweep
        dependency graph that degenerates the sweep convergence theory
        — while passing BOTH the involution check (a self-map is
        trivially involutive) and the measure check (a self-map
        trivially preserves :math:`m_n`). Third independent invariant
        of the reflection contract.

        Tangential ordinates (:math:`|\mu_a| \le` the shared
        :data:`~orpheus.numerics.spaces.angular_trace_space.TANGENTIAL_EPS`)
        are exempt: they belong to neither trace, and a tangential
        self-map is the geometrically correct image.

        Raises
        ------
        ReflectionDidNotMapInflowToOutflowError
            When any non-tangential ordinate's partner is tangential
            or lies in the same sign class.
        """
        perm = quadrature.reflection_index(self.axis)
        mu_axis = quadrature.axis_cosines(_AXIS_INDEX[self.axis])
        partner_mu = mu_axis[perm]
        active = np.abs(mu_axis) > TANGENTIAL_EPS
        partner_opposite = (np.abs(partner_mu) > TANGENTIAL_EPS) & (
            np.sign(partner_mu) != np.sign(mu_axis)
        )
        bad = active & ~partner_opposite
        if np.any(bad):
            from ._errors import ReflectionDidNotMapInflowToOutflowError

            worst = int(np.flatnonzero(bad)[0])
            raise ReflectionDidNotMapInflowToOutflowError(
                f"reflection_index({self.axis!r}) maps ordinate {worst} "
                f"(μ_{self.axis}={mu_axis[worst]:+.6e}) to ordinate "
                f"{int(perm[worst])} "
                f"(μ_{self.axis}={partner_mu[worst]:+.6e}) — same sign "
                f"class instead of the outflow partner (ERR-045: wrong "
                f"reflection-index table, or a non-axis-aligned "
                f"reflection that needs a different BC type).",
                law="reflective",
            )

    def assert_realizable(
        self,
        quadrature: "Quadrature",
        *,
        inflow_indices: Optional[np.ndarray] = None,
    ) -> None:
        r"""Universal invariants + the three reflection-table checks.

        The catalog's ERR-045 lesson verbatim: *"the inflow partition,
        the involution property, and the inflow → outflow image are
        three independent invariants. All three must hold; checking
        only one or two leaves a hole."* The base template fires the
        measure check (via the universal five); this override adds the
        involution (ERR-044) and the inflow→outflow image (ERR-045).
        """
        super().assert_realizable(quadrature, inflow_indices=inflow_indices)
        self.assert_is_involutive(quadrature)
        self.assert_reflection_maps_inflow_to_outflow(quadrature)
