r"""The mirror symmetry plane: the reflective deck law.

A partially or perfectly specular wall is a different law, a response:
:class:`~orpheus.geometry.boundary.AlbedoBoundary` with a
:class:`~orpheus.geometry.boundary.SpecularReturn` closure.

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
    r"""A symmetry plane of the domain: the deck mirror about ``axis``.

    Tensor decomposition :math:`(G_{\text{refl}}, I)`: :math:`G_{\text{refl}}`
    is the Koopman operator of the mirror motion
    :math:`\Omega \mapsto \Omega - 2(\Omega \cdot \hat{n}) \hat{n}`, realized
    as the index permutation it induces on the quadrature ordinates, and the
    response is the identity. The same operator is the
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward` of the
    angular measure under the reflection map, with the Jacobian convention
    ``|R| = 1`` since reflections are isometries.

    **A symmetry carries no amplitude.** The law states that the domain is a
    fundamental domain of the group the mirror generates, so the solution on
    the far side is the mirror image of the solution on this side; nothing is
    absorbed and nothing is emitted, and the only datum is which plane. A
    surface that returns a fraction :math:`\alpha` of its outflow specularly
    is a response, :class:`~orpheus.geometry.boundary.AlbedoBoundary` with a
    :class:`~orpheus.geometry.boundary.SpecularReturn` closure, and is a
    different value even at :math:`\alpha = 1`, where the two realize to one
    matrix. Until 2026-10-01 this law took an ``albedo`` and so spelled the
    partial wall a second time with both factors non-trivial (the fossil of
    the 2026-05 tensor-decomposition framing, and the root of ERR-094); for
    the same reason a deck law cannot be scaled or mixed
    (:mod:`~orpheus.geometry.boundary._composition`).

    This is a **pure descriptor** (Issue #186 / B3 + β2): it carries no
    ``apply`` / ``apply_transpose`` methods. Realize it via
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` to obtain the
    :class:`~orpheus.numerics.operator.PermutationOperator`:

    .. code-block:: python

        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace
        law = ReflectiveBoundary(axis="x")
        op = SNBoundaryRealizer().realize(
            law, SNMethodSpace.minimal(quad),
        )
        psi_in = op.apply(psi_out)        # forward
        # The realized operator is adjointable as well (a working
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
    """

    axis: str = "x"

    # ── The affine form's two factors (B1) ──────────────────────────────
    @property
    def geometry_map(self) -> "SelfPairedDeck":
        r""":math:`G = G_{\text{refl}}`, the mirror about :attr:`axis`."""
        return SelfPairedDeck.mirror(axis=self.axis)

    @property
    def response_kernel(self) -> "ScalarResponse":
        r""":math:`R = I`, the unit scalar response: a symmetry adds no
        physics. The same factor :class:`~orpheus.geometry.boundary.PeriodicBoundary`,
        the other deck law, declares."""
        return ScalarResponse(1.0)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, str):
            return other == self.kind
        if isinstance(other, ReflectiveBoundary):
            return self.axis == other.axis
        return NotImplemented

    def __hash__(self) -> int:
        # Hash on the canonical (post-rename) class name.
        return hash(("ReflectiveBoundary", self.axis))

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
