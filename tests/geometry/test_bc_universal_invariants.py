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
from orpheus.numerics.quadrature import Quadrature


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


@pytest.mark.foundation
class TestCreatesSweepCycleFlag:
    """The :pydata:`creates_sweep_cycle` ClassVar is the §15A.2
    signal consumed by the SN sweep planner.

    Reflective + Periodic create cycles; the rest do not.
    """

    def test_reflective_creates_sweep_cycle(self) -> None:
        assert ReflectiveBoundary.creates_sweep_cycle is True

    def test_periodic_creates_sweep_cycle(self) -> None:
        assert PeriodicBoundary.creates_sweep_cycle is True

    def test_vacuum_does_not_create_sweep_cycle(self) -> None:
        assert VacuumInflow.creates_sweep_cycle is False

    def test_white_does_not_create_sweep_cycle(self) -> None:
        assert WhiteBoundary.creates_sweep_cycle is False

    def test_albedo_does_not_create_sweep_cycle(self) -> None:
        assert AlbedoBoundary.creates_sweep_cycle is False


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

    def test_does_not_create_sweep_cycle(self) -> None:
        from orpheus.geometry.boundary import PrescribedInflow

        assert PrescribedInflow.creates_sweep_cycle is False

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
    """:class:`PrescribedInflow` realised-op semantics.

    Issue #186 (B3 + β2): descriptors are no longer callable.
    Realisation through
    :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` produces
    an :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`
    whose :meth:`apply` carries the rank-0 contract: input is ignored,
    output depends only on the source.
    """

    def test_realized_apply_with_no_source_returns_zeros(self) -> None:
        """Default ``NoSource``: realized ``apply`` returns zeros."""
        from orpheus.geometry.boundary import PrescribedInflow
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow()
        quad = Quadrature.lebedev(17)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        psi_out = np.random.default_rng(0).standard_normal((quad.N, 3))
        psi_in = op.apply(psi_out)
        assert psi_in.shape == psi_out.shape
        np.testing.assert_array_equal(psi_in, np.zeros_like(psi_out))

    def test_realized_apply_with_constant_source_is_masked_to_inflow(self) -> None:
        """``ConstantInflowSource(v)``: the realized ``apply`` delivers
        ``v`` on the face's inflow slots and ZERO everywhere else — the
        #52 / ERR-047 delivered-q contract (the pre-#52 unmasked
        full-array ``v`` was the ERR-047 hazard this file's pin used
        to bless)."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow(source=ConstantInflowSource(value=3.7))
        quad = Quadrature.lebedev(17)
        inflow = np.flatnonzero(quad.mu_x < 0)
        op = SNBoundaryRealizer().realize(
            bc,
            SNMethodSpace(
                quadrature=quad, face="xmax", inflow_indices=inflow,
            ),
        )
        psi_out = np.random.default_rng(1).standard_normal((quad.N, 2))
        psi_in = op.apply(psi_out)
        assert psi_in.shape == psi_out.shape
        np.testing.assert_array_equal(
            psi_in[inflow], np.full((inflow.size, 2), 3.7),
        )
        off_trace = np.setdiff1d(np.arange(quad.N), inflow)
        np.testing.assert_array_equal(
            psi_in[off_trace], np.zeros((off_trace.size, 2)),
        )

    def test_realized_apply_ignores_psi_out(self) -> None:
        """The rank-0 contract: input is ignored, output depends only on source."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow(source=ConstantInflowSource(value=1.0))
        quad = Quadrature.lebedev(17)
        op = SNBoundaryRealizer().realize(
            bc,
            SNMethodSpace(
                quadrature=quad, face="xmax",
                inflow_indices=np.flatnonzero(quad.mu_x < 0),
            ),
        )
        psi_out_a = np.random.default_rng(2).standard_normal((quad.N, 2))
        psi_out_b = 1000.0 * np.ones((quad.N, 2))
        # Same source → same output, regardless of psi_out.
        np.testing.assert_array_equal(
            op.apply(psi_out_a),
            op.apply(psi_out_b),
        )


@pytest.mark.l1
class TestIncomingSourceOperator:
    """:class:`IncomingSourceOperator` standalone (independent of
    :class:`PrescribedInflow`)."""

    def test_apply_returns_source_evaluation(self) -> None:
        from orpheus.geometry.boundary import ConstantInflowSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        op = IncomingSourceOperator(ConstantInflowSource(value=2.5))
        psi_out = np.zeros((8, 3))
        result = op.apply(psi_out)
        assert result.shape == (8, 3)
        np.testing.assert_array_equal(result, np.full((8, 3), 2.5))

    def test_apply_ignores_input(self) -> None:
        from orpheus.geometry.boundary import ConstantInflowSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        op = IncomingSourceOperator(ConstantInflowSource(value=1.0))
        result_a = op.apply(np.zeros((6, 2)))
        result_b = op.apply(99.0 * np.ones((6, 2)))
        np.testing.assert_array_equal(result_a, result_b)

    def test_predicates_are_apply_only(self) -> None:
        from orpheus.geometry.boundary import NoSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        op = IncomingSourceOperator(NoSource())
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
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow(source=ConstantInflowSource(value=1.5))
        quad = Quadrature.lebedev(17)
        op = SNBoundaryRealizer().realize(
            bc,
            SNMethodSpace(
                quadrature=quad, face="xmax",
                inflow_indices=np.flatnonzero(quad.mu_x < 0),
            ),
        )
        assert isinstance(op, IncomingSourceOperator)
        # The source is the same one we passed in.
        assert op.source is bc.source

    def test_realize_with_default_no_source(self) -> None:
        from orpheus.geometry.boundary import NoSource, PrescribedInflow
        from orpheus.sn.boundary.angular import IncomingSourceOperator
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow()
        quad = Quadrature.lebedev(17)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, IncomingSourceOperator)
        assert isinstance(op.source, NoSource)


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

    def test_delivered_q_vanishes_off_the_incoming_trace(self) -> None:
        """The end-to-end postcondition the invariant's mask-exists
        certification rides on: the REALIZED operator's output is the
        source value on :math:`\\Gamma_-` and exactly zero on outflow
        AND tangential slots (neither is in :math:`\\Gamma_-`)."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow(source=ConstantInflowSource(value=2.0))
        quad = Quadrature.lebedev(17)
        inflow = np.flatnonzero(quad.mu_x < 0)
        op = SNBoundaryRealizer().realize(
            bc,
            SNMethodSpace(
                quadrature=quad, face="xmax", inflow_indices=inflow,
            ),
        )
        delivered = op.apply(np.ones((quad.N, 3)))
        np.testing.assert_array_equal(
            delivered[inflow], np.full((inflow.size, 3), 2.0),
        )
        off_trace = np.setdiff1d(np.arange(quad.N), inflow)
        np.testing.assert_array_equal(
            delivered[off_trace], np.zeros((off_trace.size, 3)),
        )

    def test_mask_construction_guards(self) -> None:
        """The mask plumbing's own contract: indices without a length
        are refused; out-of-range indices are refused (mirrors
        :class:`IncomingOrdinateMaskTensor`'s construction guards)."""
        from orpheus.geometry.boundary import ConstantInflowSource
        from orpheus.sn.boundary.angular import IncomingSourceOperator

        src = ConstantInflowSource(value=1.0)
        with pytest.raises(ValueError):
            IncomingSourceOperator(
                src, inflow_indices=np.array([0, 1]),
            )
        with pytest.raises(ValueError):
            IncomingSourceOperator(
                src, inflow_indices=np.array([0, 9]), n_ordinates=4,
            )
        with pytest.raises(ValueError):
            IncomingSourceOperator(
                src,
                inflow_indices=np.array([[0], [1]]),
                n_ordinates=4,
            )
