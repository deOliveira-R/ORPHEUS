r"""Specular (reflective) boundary condition.

See :class:`SpecularBoundaryOperator` for the algebraic definition.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, ClassVar

import numpy as np

from orpheus.numerics.operator import CAP_APPLY, CAP_APPLY_TRANSPOSE

from ._base import BoundaryTraceLaw

if TYPE_CHECKING:
    from orpheus.sn.quadrature import AngularQuadrature


__all__ = ["SpecularBoundaryOperator"]


@dataclass(frozen=True)
class SpecularBoundaryOperator(BoundaryTraceLaw, key="reflective"):
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

    axis: str = "x"
    albedo: float = 1.0

    #: String tag for legacy string-kind comparisons. The default
    #: ``albedo == 1.0`` SpecularBoundaryOperator (the standard
    #: ``BC.reflective`` case) compares equal to the string
    #: ``"reflective"``; tagged SpecularBoundaryOperator instances
    #: with ``albedo != 1`` compare equal to ``"partial"`` instead.
    @property
    def kind(self) -> str:
        return "reflective" if self.albedo == 1.0 else "partial"

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, SpecularBoundaryOperator):
            return (
                self.axis == other.axis and self.albedo == other.albedo
            )
        return NotImplemented

    def __hash__(self) -> int:
        return hash(("SpecularBoundaryOperator", self.axis, self.albedo))

    def apply(
        self,
        psi_out: np.ndarray,
        quadrature: "AngularQuadrature",
    ) -> np.ndarray:
        ref = quadrature.reflection_index(self.axis)
        return self.albedo * psi_out[ref]

    def apply_transpose(
        self,
        psi_in: np.ndarray,
        quadrature: "AngularQuadrature",
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
        """
        perm = quadrature.reflection_index(self.axis)
        return self.albedo * psi_in[perm]
