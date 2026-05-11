r"""Specular (reflective) boundary condition.

See :class:`ReflectiveBoundary` for the algebraic definition. The
legacy ``SpecularBoundaryOperator`` name is re-exported as a
deprecated alias from the package ``__init__.py`` (Wave 7 rename per
Grand Report v3 vocabulary).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY, CAP_APPLY_TRANSPOSE

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["ReflectiveBoundary"]


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

    Wave-7 rename note
    ------------------
    Previously named ``SpecularBoundaryOperator``. The legacy name is
    preserved as a deprecated alias in
    ``orpheus.geometry.boundary``.

    Transpose
    ---------

    For axis reflections, the index permutation is its own inverse:
    applying the reflection twice returns each ordinate to itself.
    This makes :math:`G_{\text{refl}}^T = G_{\text{refl}}` and so the
    transpose action is identical to :meth:`apply`. The
    :pydata:`capabilities` set advertises ``apply_transpose`` --
    consumed by the sensitivity-analysis adjoint pipeline (which
    needs to assemble :math:`A^* \, y` for an operator whose
    boundary slot is reflective).

    Parameters
    ----------
    axis : str
        Axis of reflection: ``"x"``, ``"y"``, or ``"z"``. The
        :meth:`~orpheus.sn.quadrature.AngularQuadrature.reflection_index`
        method maps each ordinate to its reflected partner under this
        axis.
    albedo : float
        Specular albedo. Defaults to 1 (perfect reflection).
    """

    capabilities: ClassVar[frozenset[str]] = frozenset(
        {CAP_APPLY, CAP_APPLY_TRANSPOSE}
    )

    #: Wave-7 sweep-cycle signal (§15A.2). A reflective face maps
    #: outgoing ordinates back into the domain as inflow ordinates,
    #: creating a dependency cycle in the SN sweep DAG. The cycle
    #: detector reads this ClassVar to identify which faces need
    #: Krylov closure (rather than a single sweep).
    creates_sweep_cycle: ClassVar[bool] = True

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
        self, quadrature: "AngularQuadrature"
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
        self, quadrature: "AngularQuadrature"
    ) -> None:
        r"""Reflection is an isometry; measure preservation follows
        from involution-ness combined with per-ordinate weight
        equality between reflected partners.

        The involution check is the load-bearing piece; the
        per-ordinate weight equality is implied by the quadrature's
        construction (its weights array is symmetric under the
        reflection-index permutation by construction in every
        adapter).
        """
        self.assert_is_involutive(quadrature)

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature | None" = None,
    ) -> np.ndarray:
        r"""Apply specular reflection. ``quadrature`` is **required**.

        The geometric operator is an index permutation built from
        :meth:`~orpheus.sn.quadrature.AngularQuadrature.reflection_index`
        — without the quadrature there is no way to construct the
        permutation. Issue #176 / C176.3 makes ``quadrature``
        keyword-optional for Liskov compatibility with the abstract
        :meth:`BoundaryTraceLaw.apply` (1-arg), but ``None`` raises
        a :class:`BoundaryError` rather than silently returning an
        incorrect identity.

        For modern usage, prefer routing through
        :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`
        which captures the quadrature at realization time and
        produces a 1-arg :class:`PermutationOperator`.

        Raises
        ------
        BoundaryError
            If ``quadrature`` is ``None``.
        """
        if quadrature is None:
            from ._errors import BoundaryError
            raise BoundaryError(
                "ReflectiveBoundary.apply requires a quadrature "
                "argument (the reflection-index lookup is "
                "load-bearing). Modern usage: route through "
                "SNBoundaryRealizer.realize() with an "
                "SNMethodSpace, which captures the quadrature at "
                "construction time and produces a 1-arg "
                "PermutationOperator.",
                law="reflective",
            )
        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer,
            SNMethodSpace,
        )
        op = SNBoundaryRealizer().realize(
            self, SNMethodSpace.minimal(quadrature),
        )
        return op.apply(psi_out)

    def apply_transpose(
        self,
        psi_in: np.ndarray,
        quadrature: "AngularQuadrature | None" = None,
    ) -> np.ndarray:
        r"""Transpose action of specular reflection.

        Index permutations under axis reflection are involutions:
        applying the reflection twice maps each ordinate back to
        itself. Hence :math:`G_{\text{refl}}^T = G_{\text{refl}}`
        and the transpose acts as the same permutation scaled by
        :math:`\alpha`. Verified by the reciprocity test in
        ``tests/geometry/test_boundary.py``:

        .. math::

           \langle B(\psi_{\text{out}}), \varphi_{\text{in}} \rangle
           \;=\; \langle \psi_{\text{out}}, B^T(\varphi_{\text{in}}) \rangle.

        The reciprocity holds for the Euclidean inner product on
        the trace space (and, with quadrature-weight metadata on the
        :class:`FunctionSpace`, for the cosine-weighted inner
        product too -- the permutation commutes with diagonal
        reweighting along the ordinate axis).

        Wave 7 (C7.3): delegates to the realizer's
        :class:`PermutationOperator` (which advertises
        ``apply_transpose`` natively, scaled by albedo via the same
        ``ScaledOperator`` wrap as the forward apply).

        Issue #176 / C176.3: ``quadrature`` keyword-optional for
        Liskov compatibility; ``None`` raises :class:`BoundaryError`
        as in :meth:`apply`.
        """
        if quadrature is None:
            from ._errors import BoundaryError
            raise BoundaryError(
                "ReflectiveBoundary.apply_transpose requires a "
                "quadrature argument.",
                law="reflective",
            )
        from orpheus.sn.boundary_realizer import (
            SNBoundaryRealizer,
            SNMethodSpace,
        )
        op = SNBoundaryRealizer().realize(
            self, SNMethodSpace.minimal(quadrature),
        )
        return op.apply_transpose(psi_in)
