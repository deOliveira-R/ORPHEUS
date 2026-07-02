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
  → ERR-042 / ERR-044 (delegates to involution check)
* :meth:`WhiteBoundary.assert_response_positive_if_declared` →
  ERR-043 (``BoundaryResponseNotPositiveError``)
* :meth:`WhiteBoundary.assert_submarkov` → ERR-046
  (``SubmarkovViolationError``)
* :meth:`AlbedoBoundary.assert_response_positive_if_declared` →
  ERR-043
* :meth:`AlbedoBoundary.assert_submarkov` → ERR-046

Plan: ``.claude/plans/transient-giggling-cake.md`` Wave 7 / C7.6.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundary,
    BoundaryResponseNotPositiveError,
    PeriodicBoundary,
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


@pytest.mark.l1
class TestReflectiveMeasurePreservingInvariant:
    """:meth:`ReflectiveBoundary.assert_geometry_map_measure_preserving`
    delegates to the involution check."""

    def test_passes_for_real_quadratures(self) -> None:
        for quad in (
            Quadrature.lebedev(17),
            Quadrature.level_symmetric(sn_order=6),
            Quadrature.gauss_legendre(n_ordinates=8),
        ):
            ReflectiveBoundary(axis="x").assert_geometry_map_measure_preserving(quad)

    def test_raises_when_involution_fails(self) -> None:
        """Non-involutive perm fails the measure-preserving check
        via the delegated involution test."""
        quad = Quadrature.lebedev(17)
        N = quad.N

        class _FakeQuad:
            N = quad.N

            @staticmethod
            def reflection_index(axis: str) -> np.ndarray:
                return np.roll(np.arange(N), 1)

        with pytest.raises(ReflectionNotInvolutiveError):
            ReflectiveBoundary(axis="x").assert_geometry_map_measure_preserving(
                _FakeQuad(),  # type: ignore[arg-type]
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

    def test_realized_apply_with_constant_source_returns_constant(self) -> None:
        """``ConstantInflowSource(v)``: realized ``apply`` returns ``v``."""
        from orpheus.geometry.boundary import (
            ConstantInflowSource,
            PrescribedInflow,
        )
        from orpheus.sn.boundary.realizer import SNBoundaryRealizer
        from orpheus.sn.mesh.method_space import SNMethodSpace

        bc = PrescribedInflow(source=ConstantInflowSource(value=3.7))
        quad = Quadrature.lebedev(17)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        psi_out = np.random.default_rng(1).standard_normal((quad.N, 2))
        psi_in = op.apply(psi_out)
        assert psi_in.shape == psi_out.shape
        np.testing.assert_array_equal(psi_in, np.full_like(psi_out, 3.7))

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
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
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
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
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
