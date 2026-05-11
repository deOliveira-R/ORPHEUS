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
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
)


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
        quad = LebedevSphere.create(17)
        # No raise.
        bc.assert_is_involutive(quad)

    def test_passes_for_level_symmetric_y_axis(self) -> None:
        """Level-symmetric S6 y-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="y")
        quad = LevelSymmetricSN.create(sn_order=6)
        bc.assert_is_involutive(quad)

    def test_passes_for_gauss_legendre_x_axis(self) -> None:
        """1-D Gauss-Legendre x-axis reflection is a true involution."""
        bc = ReflectiveBoundary(axis="x")
        quad = GaussLegendre1D.create(n_ordinates=8)
        bc.assert_is_involutive(quad)

    def test_raises_for_non_involutive_perm(self) -> None:
        """A monkey-patched non-involutive ``reflection_index`` raises.

        Construct a fake quadrature whose ``reflection_index`` returns
        a 1-cycle rotation (NOT its own inverse). The invariant must
        catch this.
        """
        quad = LebedevSphere.create(17)
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
            LebedevSphere.create(17),
            LevelSymmetricSN.create(sn_order=6),
            GaussLegendre1D.create(n_ordinates=8),
        ):
            ReflectiveBoundary(axis="x").assert_geometry_map_measure_preserving(quad)

    def test_raises_when_involution_fails(self) -> None:
        """Non-involutive perm fails the measure-preserving check
        via the delegated involution test."""
        quad = LebedevSphere.create(17)
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
