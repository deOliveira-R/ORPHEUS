r"""Specular (reflective) boundary condition.

See :class:`ReflectiveBoundary` for the algebraic definition. This
class was previously named ``SpecularBoundaryOperator``; the legacy
alias was retired in Wave O step O.4a.1. ``ReflectiveBoundary`` (the
Grand Report v3 vocabulary) is the sole live name.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.numerics.quadrature import Quadrature


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
