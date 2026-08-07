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

from ._base import BoundaryTraceLaw
from ._factors import ScalarResponse, SelfPairedDeck
from ._specular import (
    assert_specular_pairing_involutive,
    assert_specular_pairing_maps_inflow_to_outflow,
    assert_specular_pairing_measure_preserving,
)

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
        Axis of reflection: ``"x"``, ``"y"``, or ``"z"`` — the mirror
        plane's NORMAL. The ordinate pairing is derived from the mirror
        motion via
        :meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`.
    albedo : float
        Specular albedo. Defaults to 1 (perfect reflection).
    """

    axis: str = "x"
    albedo: float = 1.0

    # ── The affine form's two factors (B1) ──────────────────────────────
    @property
    def geometry_map(self) -> "SelfPairedDeck":
        r""":math:`G = G_{\text{refl}}`, the mirror about :attr:`axis`."""
        return SelfPairedDeck.mirror(axis=self.axis)

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = \alpha`, the specular albedo.

        The SAME number the SN realizer already multiplies by —
        ``float(law.albedo) * base`` — now named rather than reached for.
        """
        return ScalarResponse(self.albedo)

    @property
    def kind(self) -> str:
        r"""``"reflective"`` at :math:`\alpha = 1`, else ``"partial"``.

        The ONE law that legitimately overrides the base's registry-key
        derivation, because the *declaration* vocabulary distinguishes the two:
        :meth:`~orpheus.geometry.mesh.BC.to_alpha` maps ``BC("partial",
        albedo=…)`` to its specular albedo, so a partially-reflecting face is
        declared as ``"partial"`` and must report itself the same way.

        .. warning::

           ``"partial"`` is **not** a law-registry key and appears in **no**
           method's ``BOUNDARY_OPERATOR_REGISTRY`` — so a face reporting it
           matches no admission entry, and an exact float compare
           (``albedo == 1.0``) decides which name a caller sees. That is a
           genuine semantic wrinkle, NOT an oversight, and it is deliberately
           left alone by the B0 cleanup: campaign phase **B2** removes the
           string dispatch that reads this tag at all, at which point the
           question dissolves rather than needing an answer.
        """
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
        (:math:`\pi \circ \pi = \mathrm{id}`) — ERR-044.

        Delegates to :func:`~orpheus.geometry.boundary._specular.assert_specular_pairing_involutive`;
        see that module for why the specular pairing's invariants stopped being
        this law's methods at campaign phase **B3.4b**.
        """
        assert_specular_pairing_involutive(
            quadrature, self.axis, law_key="reflective",
        )

    def assert_geometry_map_measure_preserving(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""The reflection table preserves the direction-cosine measure
        :math:`w(\Omega)\,|\Omega\cdot\hat n|` — ERR-042.

        This is the polymorphic hook the base template fires as one of the
        universal five, which is why it stays a *method* while its body lives
        with the other two pairing invariants in
        :mod:`~orpheus.geometry.boundary._specular`.
        """
        assert_specular_pairing_measure_preserving(
            quadrature, self.axis, law_key="reflective",
        )

    def assert_reflection_maps_inflow_to_outflow(
        self, quadrature: "Quadrature"
    ) -> None:
        r"""Every non-tangential ordinate reflects to the OPPOSITE sign class
        on the reflection axis — ERR-045.

        Delegates to
        :func:`~orpheus.geometry.boundary._specular.assert_specular_pairing_maps_inflow_to_outflow`.
        """
        assert_specular_pairing_maps_inflow_to_outflow(
            quadrature, self.axis, law_key="reflective",
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

        Since **B3.4b** the same three certify
        :class:`~orpheus.geometry.boundary.AlbedoBoundary` with a
        :class:`~orpheus.geometry.boundary.SpecularReturn` closure, which
        realizes through the same mirror MOTION (one spelling,
        ``_mirror_motion``; the deck kernel derives the permutation from it
        since G6.3 step 7) with the pairing
        in :math:`R` instead of :math:`G`. One certification, two laws —
        :mod:`~orpheus.geometry.boundary._specular`.
        """
        super().assert_realizable(quadrature, inflow_indices=inflow_indices)
        self.assert_is_involutive(quadrature)
        self.assert_reflection_maps_inflow_to_outflow(quadrature)
