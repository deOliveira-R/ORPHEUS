r"""Universal ``assert_*`` invariant tests (Wave 7 / C7.6).

Each concrete :class:`BoundaryTraceLaw` overrides the relevant
``assert_*`` invariants from the §16A.12 universal-invariants
catalog (Grand Report v3). This file pins the invariant contract
PER CONCRETE BC × PER APPLICABLE INVARIANT.

V&V tags
========

* ``@pytest.mark.l1`` — these tests verify a numerical contract
  (the invariant holds on physically valid inputs and raises on
  physically invalid ones), not a software invariant. They pair
  with the matching ERR-NNN error in
  ``.claude/skills/vv-principles/error_catalog.md``.
* ``@pytest.mark.catches("ERR-NNN")`` decorators wire each negative
  test to its catalogue entry per the ERR-040..ERR-047 mapping
  shipped in Wave 3.

Invariant → ERR-NNN map
=======================

* :meth:`ReflectiveBoundary.assert_is_involutive` → ERR-044
  (``ReflectionNotInvolutiveError``)
* :meth:`ReflectiveBoundary.assert_geometry_map_measure_preserving`
  → ERR-042 (``BoundaryGeometryMapNotMeasurePreservingError``) — the
  DIRECT ``w·|μ|`` pushforward check (#52; the pre-#52 delegation to
  the involution check left the unequal-weight-pairing hole open)
* :meth:`ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow`
  → ERR-045 (``ReflectionDidNotMapInflowToOutflowError``) — the
  third independent reflection-table invariant (#52)
* :meth:`BoundaryTraceLaw.assert_source_lives_on_incoming_trace` →
  ERR-047 (``BoundarySourceNotOnIncomingTraceError``) — the
  universal q ∈ Γ_- certification (#52; real body on the ABC)
* :meth:`WhiteBoundary.assert_response_positive_if_declared` →
  ERR-043 (``BoundaryResponseNotPositiveError``)
* :meth:`WhiteBoundary.assert_submarkov` → ERR-046
  (``SubmarkovViolationError``)
* :meth:`AlbedoBoundary.assert_response_positive_if_declared` →
  ERR-043
* :meth:`AlbedoBoundary.assert_submarkov` → ERR-046

Production wiring (#52): every invariant above fires at REALIZE time
through :meth:`BoundaryTraceLaw.assert_realizable` (called by
:meth:`SNBoundaryRealizer.realize` at entry), so the negative legs
below include realize-time catchers — the invariants are production
guards, not test-only helpers. The ERR-041 realizer-seam guard
(vacuum trace orientation) lives with the realizer's own tests in
``tests/sn/operators/test_sn_boundary_realizer.py``.

Plan: ``.claude/plans/transient-giggling-cake.md`` Wave 7 / C7.6;
#52 (the four unbuilt BC-layer invariants, 2026-07-12).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    PeriodicBoundary,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectiveBoundary,
    SubmarkovViolationError,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.geometry.boundary._errors import ReflectionNotInvolutiveError
from orpheus.numerics.face_layout import AXIS_NAMES
from orpheus.numerics.quadrature import Quadrature
# B3.4a: the prescribed-inflow legs below realize through the SN realizer
# (they already imported ``SNBoundaryRealizer`` in-body), and a narrowed law
# needs BOTH half-traces — which ``SNMethodSpace.minimal`` cannot name. This
# is the ONE place the face-ful method space is stood up (Pattern 2).
from tests.sn._test_helpers import face_method_space


# ─────────────────────────────────────────────────────────────────────
# ReflectiveBoundary.assert_is_involutive (ERR-044)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-044")
class TestReflectiveInvolutionInvariant:
    """:meth:`ReflectiveBoundary.assert_is_involutive` correctness."""

    def test_passes_for_default_axis_reflection_lebedev(self) -> None:
        """Lebedev-17's x-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="x")
        quad = Quadrature.lebedev(17)
        # No raise.
        bc.assert_is_involutive(quad)

    def test_passes_for_level_symmetric_y_axis(self) -> None:
        """Level-symmetric S6 y-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="y")
        quad = Quadrature.level_symmetric(sn_order=6)
        bc.assert_is_involutive(quad)

    def test_passes_for_gauss_legendre_x_axis(self) -> None:
        """1-D Gauss-Legendre x-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="x")
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        bc.assert_is_involutive(quad)

    def test_raises_for_non_involutive_perm(self) -> None:
        """A monkey-patched non-involutive ``reflection_index`` raises.

        Construct a fake quadrature whose ``reflection_index`` returns
        a 1-cycle rotation (NOT its own inverse). The invariant must
        catch this.
        """
        quad = Quadrature.lebedev(17)
        N = quad.N

        # Wrap with a fake reflection_index that returns a rotation
        # (1-cycle: i → (i+1) % N), which is NOT an involution for N ≠ 2.
        class _FakeQuad:
            N = quad.N

            @staticmethod
            def reflection_index(axis: str) -> np.ndarray:
                return np.roll(np.arange(N), 1)

        bc = ReflectiveBoundary(axis="x")
        with pytest.raises(ReflectionNotInvolutiveError):
            bc.assert_is_involutive(_FakeQuad())  # type: ignore[arg-type]


class _TableQuad:
    """Quadrature stand-in: REAL nodes/weights, injectable reflection
    table.

    The negative legs below need tables that are wrong in exactly ONE
    invariant while every other datum stays production-real — the
    surgical way to prove the three reflection-table invariants are
    independent (the ERR-045 catalog lesson: "checking only one or two
    leaves a hole").
    """

    def __init__(self, base: Quadrature, table: np.ndarray) -> None:
        self.N = base.N
        self.weights = base.weights
        self._base = base
        self._table = np.asarray(table)

    def axis_cosines(self, axis_index: int) -> np.ndarray:
        return self._base.axis_cosines(axis_index)

    def reflection_index(self, axis: str) -> np.ndarray:
        return self._table


@pytest.mark.l1
@pytest.mark.catches("ERR-042")
class TestReflectiveMeasurePreservingInvariant:
    """:meth:`ReflectiveBoundary.assert_geometry_map_measure_preserving`
    checks the ``w·|μ|`` pushforward DIRECTLY (#52).

    The pre-#52 body delegated to the involution check on the claim
    that weight symmetry is "implied by construction" — leaving open
    the exact ERR-042 hole: an involutive table pairing ordinates from
    DIFFERENT weight classes. The neighbor-pair mutant below is that
    hole made flesh; its raise is the invariant's tooth, and its clean
    pass through the involution check is the independence proof.
    """

    def test_passes_for_real_quadratures(self) -> None:
        for quad in (
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
            Quadrature.gauss_legendre(n_ordinates=8),
            Quadrature.product(n_mu=8, n_phi=8),
        ):
            for axis in ("x", "y", "z"):
                ReflectiveBoundary(axis=axis).assert_geometry_map_measure_preserving(quad)

    def test_raises_for_involutive_weight_class_mispairing(self) -> None:
        """Neighbor-pair table on GL-8: a perfect involution that
        pairs ordinates with UNEQUAL weights and |μ| — the measure
        check must red where the involution check stays green."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.array([1, 0, 3, 2, 5, 4, 7, 6]))
        bc = ReflectiveBoundary(axis="x")

        # Independence: the mutant sails through the ERR-044 check…
        bc.assert_is_involutive(mutant)  # type: ignore[arg-type]

        # …and reds the ERR-042 check.
        with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError):
            bc.assert_geometry_map_measure_preserving(mutant)  # type: ignore[arg-type]

    def test_raises_at_realize_time(self, monkeypatch) -> None:
        """The invariant is a PRODUCTION guard: the same mispaired
        table poisoning a real quadrature's precomputed partners
        reddens ``SNBoundaryRealizer.realize`` itself (via the
        ``assert_realizable`` certification at entry)."""
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        monkeypatch.setitem(
            quad.reflection_partners, 0,
            np.array([1, 0, 3, 2, 5, 4, 7, 6]),
        )
        with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError):
            SNBoundaryRealizer().realize(
                ReflectiveBoundary(axis="x"),
                SNMethodSpace.minimal(quad),
            )


# ─────────────────────────────────────────────────────────────────────
# ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow (ERR-045)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-045")
class TestReflectiveInflowToOutflowInvariant:
    """:meth:`ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow`
    (#52) — the THIRD independent reflection-table invariant.

    The identity table is the canonical ERR-045 mutant: every ordinate
    self-maps, so it is trivially involutive (ERR-044 green) and
    trivially measure-preserving (ERR-042 green — ``m[id] ≡ m``), yet
    every non-tangential inflow ordinate maps to ITSELF instead of its
    outflow partner — a length-1 self-loop in the sweep dependency
    graph. Only this check can see it.
    """

    def test_passes_for_real_quadratures(self) -> None:
        for quad in (
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
            Quadrature.gauss_legendre(n_ordinates=8),
            Quadrature.product(n_mu=8, n_phi=8),
        ):
            for axis in ("x", "y", "z"):
                ReflectiveBoundary(axis=axis).assert_reflection_maps_inflow_to_outflow(quad)

    def test_identity_table_raises_while_sibling_invariants_pass(self) -> None:
        """The independence triad on ONE mutant: identity passes
        involution AND measure, reds ONLY inflow→outflow."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.arange(quad.N))
        bc = ReflectiveBoundary(axis="x")

        bc.assert_is_involutive(mutant)  # type: ignore[arg-type]
        bc.assert_geometry_map_measure_preserving(mutant)  # type: ignore[arg-type]

        with pytest.raises(ReflectionDidNotMapInflowToOutflowError):
            bc.assert_reflection_maps_inflow_to_outflow(mutant)  # type: ignore[arg-type]

    def test_tangential_self_map_is_exempt(self) -> None:
        """A 1-D quadrature's y-axis table is the identity over
        ALL-tangential ordinates (``mu_y ≡ 0``) — the geometrically
        correct image. The exemption must keep it green (an
        over-strict check would red every degenerate axis)."""
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        np.testing.assert_array_equal(
            quad.reflection_index("y"), np.arange(quad.N),
        )
        ReflectiveBoundary(axis="y").assert_reflection_maps_inflow_to_outflow(quad)

    def test_raises_at_realize_time(self, monkeypatch) -> None:
        """Production wiring: the identity-poisoned table reddens
        ``SNBoundaryRealizer.realize`` itself."""
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        monkeypatch.setitem(
            quad.reflection_partners, 0, np.arange(quad.N),
        )
        with pytest.raises(ReflectionDidNotMapInflowToOutflowError):
            SNBoundaryRealizer().realize(
                ReflectiveBoundary(axis="x"),
                SNMethodSpace.minimal(quad),
            )


# ─────────────────────────────────────────────────────────────────────
# WhiteBoundary.assert_response_positive_if_declared (ERR-043)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-043")
class TestWhiteResponsePositiveInvariant:
    """:meth:`WhiteBoundary.assert_response_positive_if_declared`."""

    def test_passes_for_zero_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.0).assert_response_positive_if_declared()

    def test_passes_for_unit_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0).assert_response_positive_if_declared()

    def test_passes_for_partial_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.5).assert_response_positive_if_declared()

    def test_raises_for_negative_albedo(self) -> None:
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=-0.1)
        with pytest.raises(BoundaryResponseNotPositiveError):
            bc.assert_response_positive_if_declared()


# ─────────────────────────────────────────────────────────────────────
# WhiteBoundary.assert_submarkov (ERR-046)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-046")
class TestWhiteSubmarkovInvariant:
    """:meth:`WhiteBoundary.assert_submarkov`."""

    def test_passes_for_zero_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.0).assert_submarkov()

    def test_passes_for_unit_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0).assert_submarkov()

    def test_passes_for_partial_albedo(self) -> None:
        WhiteBoundary(axis="x", outward_sign=+1, albedo=0.5).assert_submarkov()

    def test_raises_for_supercritical_albedo(self) -> None:
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.2)
        with pytest.raises(SubmarkovViolationError):
            bc.assert_submarkov()


# ─────────────────────────────────────────────────────────────────────
# AlbedoBoundary.assert_response_positive_if_declared (ERR-043)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-043")
class TestAlbedoResponsePositiveInvariant:
    """:meth:`AlbedoBoundary.assert_response_positive_if_declared`."""

    def test_passes_for_zero_albedo(self) -> None:
        AlbedoBoundary(albedo=0.0).assert_response_positive_if_declared()

    def test_passes_for_unit_albedo(self) -> None:
        AlbedoBoundary(albedo=1.0).assert_response_positive_if_declared()

    def test_passes_for_partial_albedo(self) -> None:
        AlbedoBoundary(albedo=0.5).assert_response_positive_if_declared()

    def test_raises_for_negative_albedo(self) -> None:
        with pytest.raises(BoundaryResponseNotPositiveError):
            AlbedoBoundary(albedo=-0.1).assert_response_positive_if_declared()


# ─────────────────────────────────────────────────────────────────────
# AlbedoBoundary.assert_submarkov (ERR-046)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-046")
class TestAlbedoSubmarkovInvariant:
    """:meth:`AlbedoBoundary.assert_submarkov`."""

    def test_passes_for_zero_albedo(self) -> None:
        AlbedoBoundary(albedo=0.0).assert_submarkov()

    def test_passes_for_unit_albedo(self) -> None:
        AlbedoBoundary(albedo=1.0).assert_submarkov()

    def test_passes_for_partial_albedo(self) -> None:
        AlbedoBoundary(albedo=0.5).assert_submarkov()

    def test_raises_for_supercritical_albedo(self) -> None:
        with pytest.raises(SubmarkovViolationError):
            AlbedoBoundary(albedo=1.2).assert_submarkov()


# ─────────────────────────────────────────────────────────────────────
# Wave-7 sweep-cycle flag (§15A.2)
# ─────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────
# PrescribedInflow (Wave 7 / C7.5)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
class TestPrescribedInflowConstruction:
    """:class:`PrescribedInflow` construction + registry + sweep-cycle
    classification."""

    def test_default_source_is_no_source(self) -> None:
        from orpheus.geometry.boundary import NoSource, PrescribedInflow

        bc = PrescribedInflow()
        assert isinstance(bc.source, NoSource)

    def test_constant_inflow_source_attaches(self) -> None:
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        bc = PrescribedInflow(source=ConstantInflowSource(value=2.5))
        assert isinstance(bc.source, ConstantInflowSource)
        assert bc.source.value == 2.5

    def test_kind_is_prescribed_inflow(self) -> None:
        from orpheus.geometry.boundary import PrescribedInflow

        assert PrescribedInflow().kind == "prescribed_inflow"

    def test_registered_under_prescribed_inflow_key(self) -> None:
        from orpheus.geometry.boundary import (
            BoundaryTraceLaw,
            PrescribedInflow,
        )

        assert (
            BoundaryTraceLaw.registry["prescribed_inflow"] is PrescribedInflow
        )

    def test_create_via_registry(self) -> None:
        from orpheus.geometry.boundary import (
            BoundaryTraceLaw,
            PrescribedInflow,
        )

        bc = BoundaryTraceLaw.create("prescribed_inflow")
        assert isinstance(bc, PrescribedInflow)


@pytest.mark.l1
class TestPrescribedInflowApply:
    r""":class:`PrescribedInflow` realised-op semantics.

    Issue #186 (B3 + β2): descriptors are no longer callable.
    Realisation through
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` produces
    an :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`
    whose :meth:`apply` carries the rank-0 contract: input is ignored,
    output depends only on the source.

    MIGRATED at **B3.4a**, where the operator's CODOMAIN narrowed from the
    whole face slot to :math:`\Gamma_-`. The delivered-:math:`q` claim — "the
    source lands on the inflow trace and nowhere else" — is unchanged, but
    its mechanism moved from an ERASURE (emit ``N`` rows, then zero every
    non-inflow one through a mask) to an ABSENCE (emit :math:`|\Gamma_-|`
    rows; the others are not in the codomain to be written). The assertions
    follow: where they used to index the delivered array at ``off_trace`` and
    find zeros, they now assert the array HAS no such row.
    """

    def test_realized_apply_with_no_source_returns_zeros(self) -> None:
        r"""Default ``NoSource``: realized ``apply`` returns zeros on
        :math:`\Gamma_-`."""
        from orpheus.geometry.boundary import PrescribedInflow
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        bc = PrescribedInflow()
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        psi_out = np.random.default_rng(0).standard_normal(
            (space.outflow_indices.size, 3),
        )
        psi_in = op.apply(psi_out)
        assert psi_in.shape == (space.inflow_indices.size, 3)
        np.testing.assert_array_equal(psi_in, 0.0)

    def test_realized_apply_delivers_v_on_gamma_minus_and_has_no_other_row(
        self,
    ) -> None:
        r"""``ConstantInflowSource(v)``: the realized ``apply`` delivers ``v``
        on :math:`\Gamma_-` — and there is no row of the output that is not an
        inflow row.

        RE-POSED at **B3.4a** from ``..._is_masked_to_inflow``. The pre-B3.4a
        assertion read ``psi_in[off_trace] == 0`` — the mask's erasure. With
        the codomain narrowed the erasure is unspellable, so the claim it
        protected (:math:`q \in \Gamma_-`, the #52 / ERR-047 contract) is
        stated as the absence it now is: exactly :math:`|\Gamma_-|` rows,
        all carrying ``v``.

        `[M]` ``|Γ₊| == |Γ₋|`` on every quadrature × face pair in the tree,
        so a row count alone cannot distinguish :math:`\Gamma_-` from
        :math:`\Gamma_+`; what it CAN distinguish — and what the ERR-047
        hazard actually was — is a codomain that still contains the outflow
        and tangential rows, i.e. ``N``. The strict ``< quad.N`` leg is that
        discriminator, and the second half of the test makes it independent
        of the input's own leading axis.
        """
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        bc = PrescribedInflow(source=ConstantInflowSource(value=3.7))
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        inflow = space.inflow_indices
        op = SNBoundaryRealizer().realize(bc, space)
        assert op.n_inflow == inflow.size < quad.N, (
            f"the delivered q spans {op.n_inflow} rows; Γ₋ has "
            f"{inflow.size} of {quad.N} — a codomain of {quad.N} is the "
            f"pre-#52 unmasked full-face shape (ERR-047)."
        )
        psi_out = np.random.default_rng(1).standard_normal(
            (space.outflow_indices.size, 2),
        )
        psi_in = op.apply(psi_out)
        assert psi_in.shape == (inflow.size, 2)
        np.testing.assert_array_equal(psi_in, 3.7)
        # …and the row count is the OPERATOR's, not the input's: handing it
        # the whole face slot still yields |Γ₋| rows, so no caller can
        # persuade it to emit an off-trace row.
        np.testing.assert_array_equal(
            op.apply(np.ones((quad.N, 2))), np.full((inflow.size, 2), 3.7),
        )

    def test_realized_apply_ignores_psi_out(self) -> None:
        """The rank-0 contract: input is ignored, output depends only on source."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        bc = PrescribedInflow(source=ConstantInflowSource(value=1.0))
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        n_out = space.outflow_indices.size
        psi_out_a = np.random.default_rng(2).standard_normal((n_out, 2))
        psi_out_b = 1000.0 * np.ones((n_out, 2))
        # Same source → same output, regardless of psi_out.
        np.testing.assert_array_equal(
            op.apply(psi_out_a),
            op.apply(psi_out_b),
        )


@pytest.mark.l1
class TestIncomingSourceOperator:
    r""":class:`IncomingSourceOperator` standalone (independent of
    :class:`PrescribedInflow`).

    Every leg builds with ``n_inflow`` DIFFERENT from the probe's leading
    axis. `[M]` on every reachable face ``|Γ₊| == |Γ₋|``, so a fixture where
    the two agree cannot tell "the source fills the CODOMAIN" from "the
    source echoes the input's shape" — the error class sits inside the shape
    functional's invariance group (``vv`` Mode 12). Unequal sizes are the
    only way to see it, and this operator is hand-constructible, so they cost
    nothing here.
    """

    def test_apply_fills_the_codomain_not_the_input_shape(self) -> None:
        from orpheus.geometry.boundary import ConstantInflowSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        op = IncomingSourceOperator(ConstantInflowSource(value=2.5), n_inflow=5)
        psi_out = np.zeros((8, 3))  # 8 ≠ 5: the leading axes must not agree
        result = op.apply(psi_out)
        assert result.shape == (5, 3), (
            f"emitted {result.shape} for an 8-row probe with |Γ₋| = 5 — the "
            f"codomain is echoing the domain, not filling Γ₋."
        )
        np.testing.assert_array_equal(result, 2.5)

    def test_apply_ignores_input(self) -> None:
        from orpheus.geometry.boundary import ConstantInflowSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        op = IncomingSourceOperator(ConstantInflowSource(value=1.0), n_inflow=3)
        result_a = op.apply(np.zeros((6, 2)))
        result_b = op.apply(99.0 * np.ones((6, 2)))
        np.testing.assert_array_equal(result_a, result_b)
        assert result_a.shape == (3, 2)

    def test_predicates_are_apply_only(self) -> None:
        from orpheus.geometry.boundary import NoSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        op = IncomingSourceOperator(NoSource(), n_inflow=3)
        assert callable(getattr(op, "apply", None))
        # rank-0 / non-invertible — neither structural axis advertised.
        assert not op.is_invertible
        assert not op.is_adjointable


@pytest.mark.l1
class TestSNRealizerPrescribedInflowDispatch:
    """:class:`SNBoundaryRealizer` dispatches
    :class:`PrescribedInflow` to :class:`IncomingSourceOperator`."""

    def test_realize_returns_incoming_source_operator(self) -> None:
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.angular import IncomingSourceOperator
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        bc = PrescribedInflow(source=ConstantInflowSource(value=1.5))
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        assert isinstance(op, IncomingSourceOperator)
        # The source is the same one we passed in.
        assert op.source is bc.source
        # …and the realizer sized the codomain from the face's Γ₋.
        assert op.n_inflow == space.inflow_indices.size

    def test_realize_with_default_no_source(self) -> None:
        from orpheus.geometry.boundary import NoSource, PrescribedInflow
        from orpheus.sn.boundary.angular import IncomingSourceOperator
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        bc = PrescribedInflow()
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        assert isinstance(op, IncomingSourceOperator)
        assert isinstance(op.source, NoSource)

    def test_prescribed_inflow_on_a_faceless_space_is_refused(self) -> None:
        r"""B3.4a negative: even the rank-0 law needs a face.

        NET-NEW. Pre-B3.4a ``SNMethodSpace.minimal`` realized prescribed
        inflow happily (an unmasked full-face ``q`` — the ERR-047 shape). The
        law's codomain is now :math:`\Gamma_-`, which a quadrature alone
        cannot name, so the realizer refuses; the guard shipped with no
        negative test.

        The refusal is attributed to the DOMAIN guard (``outflow_indices``),
        which fires first: prescribed inflow ignores its input, but a law
        that cannot name its own domain is not realized, it is guessed.
        """
        from orpheus.geometry.boundary import BoundaryError, PrescribedInflow
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        quad = Quadrature.lebedev(17)
        with pytest.raises(BoundaryError, match="outflow_indices") as excinfo:
            SNBoundaryRealizer().realize(
                PrescribedInflow(), SNMethodSpace.minimal(quad),
            )
        assert excinfo.value.law == "prescribed_inflow"


# ─────────────────────────────────────────────────────────────────────
# BoundaryTraceLaw.assert_source_lives_on_incoming_trace (ERR-047)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-047")
class TestSourceLivesOnIncomingTraceInvariant:
    """:meth:`BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
    (#52) — the real universal body replacing the no-op default.

    The affine form's :math:`q` lives on :math:`\\Gamma_-`. An
    :class:`InflowSourceSpec` fills whatever shape it is handed, so
    the delivered-:math:`q` guarantee needs the realizer's inflow
    mask; the invariant certifies the mask EXISTS whenever
    :math:`q \\not\\equiv 0`, and the delivered-q leg pins the masked
    postcondition end-to-end (the ERR-047 mechanism: an unmasked
    constant writes into outflow slots the sweep silently discards —
    the total inflow lands SHORT of intent)."""

    def test_nonzero_source_without_inflow_set_raises(self) -> None:
        """Law-level negative: the exact catalog catcher —
        ``ConstantInflowSource(2.0)`` on a mixed inflow/outflow
        quadrature, no outflow mask available."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        quad = Quadrature.lebedev(17)
        with pytest.raises(BoundarySourceNotOnIncomingTraceError):
            law.assert_source_lives_on_incoming_trace(quad, None)

    def test_nonzero_source_with_inflow_set_passes(self) -> None:
        """Positive: with the face's inflow indices supplied, the
        realizer masks and the invariant certifies."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        quad = Quadrature.lebedev(17)
        law.assert_source_lives_on_incoming_trace(
            quad, np.flatnonzero(quad.mu_x < 0),
        )

    def test_homogeneous_laws_certify_masklessly(self) -> None:
        """Positive sweep: every homogeneous law's ``NoSource`` is
        :math:`q \\equiv 0` — trivially on :math:`\\Gamma_-`, with or
        without a mask. The universal body must not demand trace data
        it does not need."""
        quad = Quadrature.lebedev(17)
        for law in (
            VacuumInflow(),
            ReflectiveBoundary(axis="x"),
            WhiteBoundary(axis="x", outward_sign=+1),
            AlbedoBoundary(albedo=0.5),
            PeriodicBoundary(),
        ):
            law.assert_source_lives_on_incoming_trace(quad, None)

    def test_zero_valued_constant_certifies_masklessly(self) -> None:
        """``ConstantInflowSource(0.0)`` has zero support — same
        trivial certification as ``NoSource``."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )

        law = PrescribedInflow(source=ConstantInflowSource(value=0.0))
        law.assert_source_lives_on_incoming_trace(
            Quadrature.lebedev(17), None,
        )

    def test_raises_at_realize_time(self) -> None:
        """Production wiring: realizing the nonzero source on a
        quadrature-only method space (no inflow indices anywhere)
        reddens ``SNBoundaryRealizer.realize`` itself — 'the realiser
        asked to apply the source without an outflow mask' (the
        catalog's catcher spec, verbatim)."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        quad = Quadrature.lebedev(17)
        with pytest.raises(BoundarySourceNotOnIncomingTraceError):
            SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))

    def test_delivered_q_has_no_row_off_the_incoming_trace(self) -> None:
        r"""The end-to-end postcondition the invariant's mask-exists
        certification rides on: the REALIZED operator delivers the source
        value on :math:`\Gamma_-`, and there is NO row of its output that is
        an outflow or tangential slot.

        RE-POSED at **B3.4a** from ``..._vanishes_off_the_incoming_trace``.
        The pre-B3.4a assertion indexed the delivered array at ``off_trace``
        and found zeros — the mask's erasure, which is the postcondition
        ERR-047 names. With the codomain narrowed to :math:`\Gamma_-` those
        rows are not emitted at all, so the ERASURE has become an ABSENCE and
        this leg asserts the absence: the delivered block has exactly
        :math:`|\Gamma_-|` rows out of ``N``, every one carrying ``v``.

        The ``catches("ERR-047")`` claim stays LIVE, and the FULL-FACE probe
        is what keeps it live. `[M]` ``|Γ₊| == |Γ₋|`` on every reachable face,
        so a :math:`\Gamma_+`-sized probe cannot tell "``q`` fills the
        codomain" from "``q`` echoes whatever it is handed" — the very
        regression that IS ERR-047 (an unmasked ``q`` spanning the whole face
        slot, whose outflow entries the sweep then discards, leaving the
        total inflow SHORT). Feeding the whole slot and requiring
        :math:`|\Gamma_-|` rows back is the discriminator; without it this
        leg sits inside the shape functional's invariance group (``vv``
        Mode 12) and the marker would be a phantom.
        """
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        bc = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        inflow = space.inflow_indices
        # The face partition is THREE-way — inflow ⊔ outflow ⊔ tangential —
        # so "no off-trace row" is counted against N, never against |Γ₊|.
        off_trace = np.setdiff1d(np.arange(quad.N), inflow)
        assert off_trace.size > 0  # activation: there IS an off-trace to miss
        op = SNBoundaryRealizer().realize(bc, space)
        delivered = op.apply(np.ones((space.outflow_indices.size, 3)))
        assert delivered.shape == (inflow.size, 3), (
            f"delivered q spans {delivered.shape[0]} rows; Γ₋ has "
            f"{inflow.size} and the face slot has {quad.N}. Any row beyond "
            f"Γ₋ is an off-trace slot the sweep would discard (ERR-047)."
        )
        np.testing.assert_array_equal(delivered, 2.0)
        # The ERR-047 discriminator: the whole face slot in, |Γ₋| rows out.
        full_face = op.apply(np.ones((quad.N, 3)))
        assert full_face.shape == (inflow.size, 3), (
            f"a full-face probe produced {full_face.shape[0]} rows — the "
            f"delivered q is sized from its INPUT, so it spans the whole "
            f"face slot again and its off-trace entries are the ERR-047 "
            f"source the sweep discards."
        )
        np.testing.assert_array_equal(full_face, 2.0)

    def test_codomain_size_is_the_construction_contract(self) -> None:
        r"""The operator's construction contract: :math:`|\Gamma_-|` is the
        one datum it needs, and it must be a real codomain size.

        RE-POSED at **B3.4a** from ``test_mask_construction_guards``. The
        mask plumbing whose guards that test pinned (``inflow_indices`` +
        ``n_ordinates``, with rank and out-of-range checks) is RETIRED — the
        codomain no longer contains the rows a mask would zero. The claim
        those guards protected, ":math:`q` lands on :math:`\Gamma_-` and
        nowhere else", is now carried by the codomain size itself, so that is
        what this leg guards: a negative row count is refused, and a valid
        one is exactly the number of rows :meth:`apply` emits — independent
        of the probe's own leading axis (``vv`` Mode 12: on every reachable
        face ``|Γ₊| == |Γ₋|``, so an equal-sized probe could not tell the two
        apart).
        """
        from orpheus.geometry.boundary import ConstantInflowSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        src = ConstantInflowSource(value=1.0)
        with pytest.raises(ValueError, match="n_inflow must be non-negative"):
            IncomingSourceOperator(src, n_inflow=-1)
        # Positive control: a real codomain size constructs, and the emitted
        # row count IS that size (7 ≠ 3, so the probe cannot supply it).
        op = IncomingSourceOperator(src, n_inflow=3)
        assert op.n_inflow == 3
        assert op.apply(np.ones((7, 2))).shape == (3, 2)


# ─────────────────────────────────────────────────────────────────────
# B3.4b — ONE certification, TWO carriers of the specular pairing
# ─────────────────────────────────────────────────────────────────────


def _non_involutive_but_otherwise_valid(
    quadrature: Quadrature, *, axis: str,
) -> np.ndarray:
    r"""A table that breaks ONLY the involution (ERR-044).

    Why not the obvious ``np.roll(arange(N), 1)``: an n-cycle rotation breaks
    all three invariants at once (it pairs unequal weight classes, so ERR-042
    fires FIRST in the aggregate and the row proves nothing about ERR-044's
    independence). Worse, **no pure 3-cycle can ever be the ERR-044-only
    mutant**: a cycle :math:`a \to b \to c \to a` needs
    :math:`\mathrm{sign}(b) = -\mathrm{sign}(a)`,
    :math:`\mathrm{sign}(c) = -\mathrm{sign}(b) = \mathrm{sign}(a)` and then
    :math:`\mathrm{sign}(a) = -\mathrm{sign}(c) = -\mathrm{sign}(a)` — a
    contradiction, so every odd cycle breaks ERR-045 too.

    The construction that works is :math:`\pi \circ \sigma`, where :math:`\pi`
    is the TRUE mirror and :math:`\sigma` swaps two ordinates of the SAME sign
    and the SAME measure :math:`w|\mu_a|`:

    * **measure** — :math:`m[\pi[\sigma[i]]] = m[\sigma[i]] = m[i]`, since
      :math:`\sigma` stays inside a measure class and :math:`\pi` preserves it;
    * **sign** — :math:`\mathrm{sign}(\pi[\sigma[i]]) = -\mathrm{sign}(i)`,
      since :math:`\sigma` preserves sign and :math:`\pi` flips it;
    * **involution** — BROKEN: for the swapped pair :math:`\{a, b\}`,
      :math:`T[T[a]] = \pi[\sigma[\pi[b]]] = \pi[\pi[b]] = b \neq a`
      (:math:`\sigma` does not touch :math:`\pi[b]`, which has the opposite
      sign).
    """
    perm = quadrature.reflection_index(axis)
    mu = quadrature.axis_cosines(AXIS_NAMES.index(axis))
    measure = quadrature.weights * np.abs(mu)
    # Group same-sign ordinates by measure; the first class with ≥2 members
    # supplies the swap. Rounded because the class is a float coincidence of
    # a symmetric node set, not an exact-compare question.
    classes: dict[tuple[int, float], list[int]] = {}
    for n in range(int(quadrature.N)):
        classes.setdefault(
            (int(np.sign(mu[n])), round(float(measure[n]), 12)), []
        ).append(n)
    for members in classes.values():
        if len(members) >= 2:
            a, b = members[0], members[1]
            sigma = np.arange(int(quadrature.N))
            sigma[a], sigma[b] = b, a
            return perm[sigma]
    raise AssertionError(
        f"{quadrature!r} has no same-sign measure class with two members, so "
        f"the ERR-044-only mutant is not constructible here — pick a "
        f"quadrature with a degenerate measure class."
    )


@pytest.mark.l1
@pytest.mark.catches("ERR-042", "ERR-044", "ERR-045")
class TestSpecularPairingCertifiedOnBothCarriers:
    r"""The three table invariants fire for BOTH laws that stand on the pairing.

    Campaign phase **B3.4b** moved the specular pairing's certification out of
    :class:`ReflectiveBoundary`'s methods and into
    :mod:`~orpheus.geometry.boundary._specular`, because the 2026-08-01 ruling
    gave :class:`AlbedoBoundary` a :class:`SpecularReturn` closure that realizes
    through the SAME ``quadrature.reflection_index(axis)`` table — with the
    pairing in :math:`R` instead of :math:`G`.

    That module's central claim is *"one certification, two laws"*, and until
    this class it was a **docstring claim with no test**: MEASURED
    2026-08-01, deleting the
    ``isinstance(self.reemission, SpecularReturn)`` clause from
    ``AlbedoBoundary.assert_realizable`` left the entire B3.4b gate file
    (``test_reemission_closure.py``, 178 tests) GREEN. A wrong table would have
    been caught on the reflective route and **silently realized** on the albedo
    one — exactly the twin-path failure this campaign exists to remove.

    Structure: each broken table is wrong in exactly ONE invariant (the
    independence the ERR-045 lesson demands), and each row asserts BOTH
    carriers raise the SAME error type with their OWN ``law=`` attribution.
    The attribution leg is what makes this more than a duplicate of the three
    classes above — it proves the shared module reports the caller, not a
    hard-coded ``"reflective"``.
    """

    #: ``(id, quadrature_factory, table_factory, error_type)``. Each table is
    #: production-real except for the ONE property it breaks.
    #:
    #: The QUADRATURE is per-row, not shared, and that is load-bearing: the
    #: ERR-044-only mutant needs a **degenerate measure class** (two same-sign
    #: ordinates carrying equal ``w·|μ_x|``) to swap inside, and MEASURED,
    #: ``gauss_legendre(8)`` has none — its four per-side measures are all
    #: distinct. ``level_symmetric(6)`` supplies the degeneracy.
    _BROKEN = [
        # Neighbour-pair on GL-8: a perfect involution pairing UNEQUAL weight
        # classes — ERR-042 only (the involution and sign checks stay green).
        ("weight_class_mispairing",
         lambda: Quadrature.gauss_legendre(n_ordinates=8),
         lambda q: np.array([1, 0, 3, 2, 5, 4, 7, 6]),
         BoundaryGeometryMapNotMeasurePreservingError),
        # The identity table: every ordinate self-maps — trivially involutive
        # and trivially measure-preserving, but every inflow ordinate maps to
        # ITSELF instead of its outflow partner. ERR-045 only.
        ("self_map_identity",
         lambda: Quadrature.gauss_legendre(n_ordinates=8),
         lambda q: np.arange(q.N),
         ReflectionDidNotMapInflowToOutflowError),
        # The true mirror composed with a swap of two SAME-SIGN, SAME-MEASURE
        # ordinates: measure-preserving and sign-crossing, but NOT its own
        # inverse. ERR-044 only. See :func:`_non_involutive_but_otherwise_valid`
        # for why the obvious n-cycle will not do.
        ("mirror_after_same_class_swap",
         lambda: Quadrature.level_symmetric(sn_order=6),
         lambda q: _non_involutive_but_otherwise_valid(q, axis="x"),
         ReflectionNotInvolutiveError),
    ]

    #: The two carriers, with the tier each holds the pairing in.
    def _carriers(self):
        from orpheus.geometry.boundary import SpecularReturn

        return [
            ("reflective", ReflectiveBoundary(axis="x", albedo=1.0), "G"),
            ("albedo", AlbedoBoundary(0.5, SpecularReturn(axis="x")), "R"),
        ]

    @pytest.mark.parametrize(
        "table_id,build_quad,build_table,error_type", _BROKEN,
        ids=[row[0] for row in _BROKEN],
    )
    def test_both_carriers_reject_the_same_broken_table(
        self, table_id, build_quad, build_table, error_type,
    ) -> None:
        """A table broken in ONE invariant reddens BOTH laws' certification."""
        quad = build_quad()
        mutant = _TableQuad(quad, build_table(quad))
        seen = {}
        for law_id, law, tier in self._carriers():
            with pytest.raises(error_type) as exc:
                law.assert_realizable(mutant)  # type: ignore[arg-type]
            seen[law_id] = exc.value.law
        assert seen == {"reflective": "reflective", "albedo": "albedo"}, (
            f"{table_id}: each carrier must attribute the failure to ITSELF; "
            f"got {seen}. A hard-coded law_key in "
            f"orpheus.geometry.boundary._specular would blame the wrong law "
            f"and send a reader to the wrong file."
        )

    def test_the_real_table_certifies_on_both_carriers(self) -> None:
        """POSITIVE control (``vv-principles`` anti-#11).

        Without it the negatives above validate only the RAISING behaviour,
        not the invariant claim: a certification that raised on every table
        would pass all three rows.
        """
        for quad in (
            Quadrature.gauss_legendre(n_ordinates=8),
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
        ):
            for law_id, law, tier in self._carriers():
                law.assert_realizable(quad)

    def test_the_albedo_carrier_leaves_the_G_HOOK_a_no_op(self) -> None:
        r"""The TIER claim, not just the check-fires claim.

        ``AlbedoBoundary``'s :math:`G` **is** the identity map and **does**
        preserve the measure, so the base template's polymorphic
        ``assert_geometry_map_measure_preserving`` hook is correctly a no-op
        here even on a table that is catastrophically broken — the pairing is
        certified as a law-specific extension because it lives in :math:`R`.

        This is the gate that reds if someone "closes the hole" by overriding
        the :math:`G` hook instead of extending ``assert_realizable``: that
        would re-assert the very conflation the 2026-08-01 ruling forbids
        (a wall's constitutive return posing as a symmetry of the domain).
        """
        from orpheus.geometry.boundary import SpecularReturn

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.array([1, 0, 3, 2, 5, 4, 7, 6]))
        law = AlbedoBoundary(0.5, SpecularReturn(axis="x"))
        # The G hook is silent — G is the identity deck element, which fixes
        # no geometry.
        law.assert_geometry_map_measure_preserving(mutant)  # type: ignore[arg-type]
        # …while the aggregate still reds, via the R-tier extension.
        with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError):
            law.assert_realizable(mutant)  # type: ignore[arg-type]

    def test_a_diffuse_closure_does_NOT_fire_the_pairing_checks(self) -> None:
        """SCOPE control: only a SPECULAR closure stands on the table.

        An ``IsotropicReturn`` closure realizes through the Lambertian average
        and never touches ``reflection_index``, so a broken table must leave it
        untouched. Without this leg, an ``assert_realizable`` that fired the
        pairing checks unconditionally would pass every row above while
        refusing perfectly good diffuse laws.
        """
        from orpheus.geometry.boundary import IsotropicReturn

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        mutant = _TableQuad(quad, np.roll(np.arange(quad.N), 1))
        AlbedoBoundary(
            0.5, IsotropicReturn(axis="x", outward_sign=+1),
        ).assert_realizable(mutant)  # type: ignore[arg-type]
        # …and the bare law likewise (no closure, no pairing).
        AlbedoBoundary(0.5).assert_realizable(mutant)  # type: ignore[arg-type]

    def test_raises_at_realize_time_on_the_albedo_route(self, monkeypatch) -> None:
        """PRODUCTION-seam leg: the certification is wired at
        ``SNBoundaryRealizer.realize``, not only reachable by hand.

        The sibling ``TestReflectiveMeasurePreservingInvariant`` proves this
        for the reflective route; B3.4b's new route needs its own, because the
        realizer reaches ``assert_realizable`` through a different dispatch arm.
        """
        from orpheus.geometry.boundary import SpecularReturn
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        quad = Quadrature.gauss_legendre(n_ordinates=8)
        monkeypatch.setitem(
            quad.reflection_partners, 0,
            np.array([1, 0, 3, 2, 5, 4, 7, 6]),
        )
        with pytest.raises(BoundaryGeometryMapNotMeasurePreservingError) as exc:
            SNBoundaryRealizer().realize(
                AlbedoBoundary(0.5, SpecularReturn(axis="x")),
                face_method_space(quad, face="xmax"),
            )
        assert exc.value.law == "albedo"


# ─────────────────────────────────────────────────────────────────────
# B3.4c — the factor tier's adjointability claim must match what the
#         realizer actually builds
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestFactorAdjointabilityMatchesTheRealizedOperator:
    r"""``G.is_adjointable and R.is_adjointable == realized.is_adjointable``.

    Two independent encodings of one claim, compared — the ERR-041 discipline
    applied to adjointability. The factor tier DECLARES whether the map it
    names exposes an honest transpose; the realized operator ANSWERS from its
    own composition tree. Nothing compared them before **B3.4c**, and they had
    silently drifted apart: :class:`SpatialWrap` declared ``False`` while the
    operator realizing it answered ``True``, so a consumer got opposite
    answers depending on which one it asked. (The declaration was reporting an
    implementation gap — #183, periodic's unbuilt partner channel — in a slot
    whose contract is a property of the MAP. B3.4c built the channel, the
    declaration became true, and this gate is what keeps the two tiers from
    parting again.)

    The conjunction is the right shape because the realized operator is
    :math:`R \circ G`: a composition is Euclidean-adjointable exactly when both
    factors are. It bites in both directions — a factor that over-claims
    (``True`` with no realized transpose) and one that under-claims (the
    B3.4c-era ``SpatialWrap``) are both red.
    """

    #: ``(law_id, law)`` over every law with a realizable SN arm. The albedo
    #: rows carry both closures because since B3.4b an albedo face answers by
    #: its CLOSURE, not its class — a specular closure is adjointable and a
    #: diffuse one is not, so a single albedo row would pin only one arm.
    def _laws(self):
        from orpheus.geometry.boundary import IsotropicReturn, SpecularReturn

        return [
            ("vacuum", VacuumInflow()),
            ("reflective", ReflectiveBoundary(axis="x", albedo=1.0)),
            ("reflective_partial", ReflectiveBoundary(axis="x", albedo=0.5)),
            ("white", WhiteBoundary(axis="x", outward_sign=-1, albedo=1.0)),
            ("albedo_specular",
             AlbedoBoundary(0.5, SpecularReturn(axis="x"))),
            ("albedo_isotropic",
             AlbedoBoundary(0.5, IsotropicReturn(axis="x", outward_sign=-1))),
            ("periodic", PeriodicBoundary(axis="x")),
        ]

    def test_every_law_agrees_with_its_realization(self) -> None:
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        quad = Quadrature.level_symmetric(sn_order=6)
        space = face_method_space(
            quad, face="xmin", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        realizer = SNBoundaryRealizer()
        disagreed = {}
        for law_id, law in self._laws():
            declared = (
                law.geometry_map.is_adjointable
                and law.response_kernel.is_adjointable
            )
            realized = realizer.realize(law, space).is_adjointable
            if declared != realized:
                disagreed[law_id] = (declared, realized)
        assert not disagreed, (
            f"the factor tier and the realized operator disagree on "
            f"adjointability for {disagreed} (declared, realized). Whichever "
            f"is wrong, a consumer reading one gets an answer the other "
            f"contradicts — which is the two-sources-of-truth B3.4c closed "
            f"for periodic."
        )

    def test_both_answers_are_actually_exercised(self) -> None:
        """The row set must contain a ``True`` AND a ``False``.

        Without this the agreement above is satisfied by a set that is all one
        value, and a mutation flipping every declaration at once would pass.
        White and the diffuse albedo closure supply the ``False`` (the
        Lambertian is self-adjoint only under the cosine-weighted metric).
        """
        answers = {
            law.geometry_map.is_adjointable
            and law.response_kernel.is_adjointable
            for _law_id, law in self._laws()
        }
        assert answers == {True, False}, (
            f"the law set exercises only {answers}; a one-valued set makes "
            f"the agreement gate blind to a uniform flip."
        )

    def test_prescribed_inflow_is_the_documented_exception(self) -> None:
        r"""An AFFINE law is the one place the conjunction does NOT hold — and
        it is a statement about the algebra, not a gap.

        :class:`PrescribedInflow` carries both factors trivially adjointable
        (:math:`G = \mathrm{id}`, :math:`R = 0`) while its realized
        :class:`IncomingSourceOperator` declines a transpose. Nothing is
        inconsistent: the factor pair describes the LINEAR part
        :math:`R G \gamma_+\psi`, and this law's entire content is the affine
        term :math:`q`, which the tier does not carry. Pinned as a row rather
        than left out of the set above, so "prescribed inflow is missing" can
        never be an oversight that hides a real drift.
        """
        from orpheus.geometry.boundary import PrescribedInflow
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer

        law = PrescribedInflow()
        assert law.geometry_map.is_adjointable
        assert law.response_kernel.is_adjointable
        quad = Quadrature.level_symmetric(sn_order=6)
        realized = SNBoundaryRealizer().realize(
            law,
            face_method_space(
                quad, face="xmin", faces=("xmin", "xmax", "ymin", "ymax"),
            ),
        )
        assert not realized.is_adjointable, (
            "IncomingSourceOperator advertised a transpose — if the affine "
            "source became adjointable, this law joins the conjunction gate "
            "above instead of standing outside it."
        )
