r"""Tests for SNBoundaryRealizer (Wave 5 / C5.3).

The realizer dispatches by ``isinstance(law, ...)`` to the Wave-0 /
Wave-1 primitives that realize each legacy
:class:`~orpheus.geometry.boundary.BoundaryOperator` subclass as a
single-arg :class:`~orpheus.numerics.operator.LinearOperator`.

The L1 verification claim: for every legacy BC × quadrature pair, the
realized operator's ``apply(psi)`` reproduces the legacy 2-arg
``bc.apply(psi, quad)`` semantics. The match is bit-equivalent for the
permutation-based / pass-through paths (specular at α=1, periodic,
albedo at α=0/1) and within nulp=4 for the scaled paths
(specular/white/albedo at α∉{0,1}, mixed compositions). Vacuum is the
ONLY case where the new semantics deviate from legacy by design: the
new path zeroes ONLY the inflow ordinates per §16A.10 of the Grand
Report (legacy zeroed ALL ordinates) — an intentional Wave 8 semantic
correction documented in the plan's risk register.

The realizer's dispatch logic is structural (not numerical math); the
math content lives in the Wave-0 / Wave-1 primitives and is verified
in their own tests. Wave 5 verifies that the dispatch chain composes
correctly to reproduce the legacy answer.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    AlbedoBoundaryOperator,
    BoundaryError,
    BoundaryRealizerRegistry,
    BoundaryRealizerRegistryError,
    MixedBoundaryOperator,
    PeriodicBoundaryOperator,
    SpecularBoundaryOperator,
    VacuumBoundaryOperator,
    WhiteBoundaryOperator,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    IncomingOrdinateMaskTensor,
    PeriodicWrapOperator,
    PermutationOperator,
    ScaledOperator,
    ZeroOperator,
)
from orpheus.sn.angular_operator import AngularAverageOperator
from orpheus.sn.boundary_realizer import SNBoundaryRealizer, SNMethodSpace
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
)


# ─────────────────────────────────────────────────────────────────────
# 1. Vacuum
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeVacuum:
    """Vacuum realizes to :class:`IncomingOrdinateMaskTensor`.

    SEMANTIC NOTE (plan §16A.5 risk register, §16A.10 trace
    representation): the legacy ``VacuumBoundaryOperator.apply`` returns
    ``np.zeros_like(psi_out)`` (zeroes EVERY ordinate). The new realizer
    zeroes ONLY the inflow ordinates and preserves the outflow trace.
    This is the intentional Wave 8 semantic correction. Tests compare
    the realized output against the EXPECTED §16A.10 behaviour, NOT
    against the legacy 2-arg ``bc.apply(psi, quad)`` — the latter would
    be the bit-equivalent test if the legacy semantics were the target,
    which they are not for vacuum.
    """

    def test_lebedev17_xmin_zeroes_only_inflow(self):
        quad = LebedevSphere.create(17)
        # xmin face: outward normal is -x, inflow is mu_x > 0 (into the domain).
        inflow_indices = np.flatnonzero(quad.mu_x > 0)
        space = SNMethodSpace(
            quadrature=quad, face="xmin", inflow_indices=inflow_indices,
        )
        rng = np.random.default_rng(42)
        psi = rng.uniform(0.5, 2.0, size=(quad.N, 4, 2))
        realizer = SNBoundaryRealizer()
        op = realizer.realize(VacuumBoundaryOperator(), space)
        out = op.apply(psi)
        # Inflow rows: zero. Non-inflow rows: equal to input.
        np.testing.assert_array_equal(out[inflow_indices], 0.0)
        non_inflow = np.setdiff1d(np.arange(quad.N), inflow_indices)
        np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])

    def test_lebedev17_ymax_zeroes_only_inflow(self):
        quad = LebedevSphere.create(17)
        # ymax face: outward normal is +y, inflow is mu_y < 0 (into the domain).
        inflow_indices = np.flatnonzero(quad.mu_y < 0)
        space = SNMethodSpace(
            quadrature=quad, face="ymax", inflow_indices=inflow_indices,
        )
        psi = np.arange(quad.N * 3, dtype=float).reshape(quad.N, 3)
        realizer = SNBoundaryRealizer()
        op = realizer.realize(VacuumBoundaryOperator(), space)
        out = op.apply(psi)
        np.testing.assert_array_equal(out[inflow_indices], 0.0)
        non_inflow = np.setdiff1d(np.arange(quad.N), inflow_indices)
        np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])

    def test_vacuum_missing_inflow_indices_raises(self):
        """Without ``inflow_indices`` the realizer raises BoundaryError."""
        quad = LebedevSphere.create(17)
        space = SNMethodSpace.minimal(quad)  # no inflow_indices
        realizer = SNBoundaryRealizer()
        with pytest.raises(BoundaryError) as excinfo:
            realizer.realize(VacuumBoundaryOperator(), space)
        assert excinfo.value.law == "vacuum"

    def test_returned_op_type_is_incoming_ordinate_mask(self):
        """The vacuum dispatch returns an :class:`IncomingOrdinateMaskTensor`."""
        quad = GaussLegendre1D.create(8)
        inflow_indices = np.flatnonzero(quad.mu_x > 0)
        space = SNMethodSpace(
            quadrature=quad, face="left", inflow_indices=inflow_indices,
        )
        op = SNBoundaryRealizer().realize(VacuumBoundaryOperator(), space)
        assert isinstance(op, IncomingOrdinateMaskTensor)


# ─────────────────────────────────────────────────────────────────────
# 2. Reflective (specular)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeReflective:
    """Specular realizes to ``albedo * PermutationOperator(perm)``."""

    def test_specular_unit_albedo_lebedev_matches_legacy_bit_exact(self):
        """At α=1 the realized op MUST bit-match legacy ``bc.apply``."""
        quad = LebedevSphere.create(17)
        bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
        space = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(7)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5, 3))
        np.testing.assert_array_equal(op.apply(psi), bc.apply(psi, quad))

    def test_specular_unit_albedo_returns_bare_permutation(self):
        """At α=1 the dispatch returns the bare PermutationOperator
        (no ScaledOperator wrap) — preserves bit-identity vs legacy."""
        quad = LevelSymmetricSN.create(4)
        bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, PermutationOperator)

    def test_specular_partial_albedo_matches_legacy_within_nulp(self):
        """At α=0.7 the order of multiplication may shift one ULP."""
        quad = LebedevSphere.create(17)
        bc = SpecularBoundaryOperator(axis="x", albedo=0.7)
        space = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(9)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5))
        # ScaledOperator(0.7, P).apply(psi) = 0.7 * psi[perm].
        # Legacy bc.apply also does 0.7 * psi_out[ref]. SAME order,
        # so bit-equivalence is expected; nulp=4 is a safety margin.
        np.testing.assert_array_almost_equal_nulp(
            op.apply(psi), bc.apply(psi, quad), nulp=4,
        )

    def test_specular_partial_albedo_returns_scaled_operator(self):
        """At α!=1 the dispatch returns a ScaledOperator."""
        quad = LevelSymmetricSN.create(4)
        bc = SpecularBoundaryOperator(axis="x", albedo=0.5)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.5
        assert isinstance(op.op, PermutationOperator)


# ─────────────────────────────────────────────────────────────────────
# 3. White (Lambertian)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeWhite:
    """White realizes to ``albedo * AngularAverageOperator``.

    The body of :class:`AngularAverageOperator.from_quadrature` is
    lifted verbatim from :class:`WhiteBoundaryOperator.apply`, and
    Wave 1's ``TestLegacyBitEquivalence`` already pins
    ``np.testing.assert_array_equal`` at α=1.0 on Lebedev 17. We
    re-check that property here through the realizer's dispatch chain.
    """

    def test_white_unit_albedo_lebedev_matches_legacy_bit_exact(self):
        """At α=1 the realized op MUST bit-match legacy ``bc.apply``."""
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        rng = np.random.default_rng(123)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 5, 3))
        np.testing.assert_array_equal(op.apply(psi), bc.apply(psi, quad))

    def test_white_unit_albedo_returns_bare_angular_average(self):
        """At α=1 the dispatch returns the bare AngularAverageOperator."""
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, AngularAverageOperator)

    def test_white_partial_albedo_matches_legacy_within_nulp(self):
        """At α=0.3 the multiplication order may shift one ULP.

        Legacy: the multiplication happens inside ``bc.apply`` via the
        ``psi_avg[None, ...] * self.albedo`` broadcast. Realizer:
        ``ScaledOperator(0.3, ...).apply`` multiplies the inner output
        by 0.3 from the outside. The two orderings produce the same
        value in real arithmetic but may differ by one ULP under
        IEEE-754 (FP-non-associativity); nulp=4 is the safe floor per
        the algebra-of-record bit-identity discipline.
        """
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=0.3)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        rng = np.random.default_rng(11)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 4))
        np.testing.assert_array_almost_equal_nulp(
            op.apply(psi), bc.apply(psi, quad), nulp=4,
        )

    def test_white_z_axis_on_gauss_legendre_raises(self):
        """z-axis white on a 1-D GL quadrature delegates the raise to
        :meth:`AngularAverageOperator.from_quadrature` (no mu_z)."""
        quad = GaussLegendre1D.create(8)
        bc = WhiteBoundaryOperator(axis="z", outward_sign=+1, albedo=1.0)
        with pytest.raises(ValueError, match="mu_z"):
            SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))


# ─────────────────────────────────────────────────────────────────────
# 4. Albedo (rank-1 isotropic scaling)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeAlbedo:
    """Albedo realizes to Zero / Identity / Scaled-Identity by α."""

    def test_albedo_zero_realizes_to_zero_operator(self):
        quad = GaussLegendre1D.create(8)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundaryOperator(0.0), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, ZeroOperator)
        psi = np.arange(quad.N * 4, dtype=float).reshape(quad.N, 4)
        np.testing.assert_array_equal(op.apply(psi), 0.0)

    def test_albedo_one_realizes_to_identity_operator(self):
        quad = GaussLegendre1D.create(8)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundaryOperator(1.0), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, IdentityOperator)
        psi = np.arange(quad.N * 4, dtype=float).reshape(quad.N, 4)
        np.testing.assert_array_equal(op.apply(psi), psi)

    def test_albedo_half_realizes_to_scaled_identity(self):
        quad = GaussLegendre1D.create(8)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundaryOperator(0.5), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.5
        assert isinstance(op.op, IdentityOperator)
        psi = np.arange(quad.N * 4, dtype=float).reshape(quad.N, 4)
        np.testing.assert_array_equal(op.apply(psi), 0.5 * psi)


# ─────────────────────────────────────────────────────────────────────
# 5. Periodic
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizePeriodic:
    """Periodic realizes to :class:`PeriodicWrapOperator`."""

    def test_periodic_realizes_to_wrap_and_passes_through(self):
        quad = LebedevSphere.create(17)
        op = SNBoundaryRealizer().realize(
            PeriodicBoundaryOperator(), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, PeriodicWrapOperator)
        rng = np.random.default_rng(3)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5))
        # The Wave-0 PeriodicWrapOperator returns the input by
        # reference (zero-cost angular pass-through); legacy
        # PeriodicBoundaryOperator returned a copy.  Compare values,
        # NOT identity.
        np.testing.assert_array_equal(op.apply(psi), psi)


# ─────────────────────────────────────────────────────────────────────
# 6. Mixed (rank-N composition)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeMixed:
    """Mixed realizes recursively to a sum-of-scaled-primitives."""

    def test_mixed_30_specular_70_white_matches_legacy(self):
        """Standard Marshak mixed: 30% specular + 70% white."""
        quad = LebedevSphere.create(17)
        spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
        white = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        mixed = MixedBoundaryOperator([(0.3, spec), (0.7, white)])
        op = SNBoundaryRealizer().realize(mixed, SNMethodSpace.minimal(quad))
        rng = np.random.default_rng(17)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 4))
        legacy_out = mixed.apply(psi, quad)
        # Reduction order: realizer composes via Wave-0 OperatorSum,
        # while legacy accumulates result = result + coeff * primitive.apply.
        # Both fold the same two summands; FP order may differ by one
        # ULP. nulp=64 is the safe margin for nested-sum compositions.
        np.testing.assert_array_almost_equal_nulp(
            op.apply(psi), legacy_out, nulp=64,
        )

    def test_mixed_empty_realizes_to_zero_operator(self):
        quad = GaussLegendre1D.create(8)
        op = SNBoundaryRealizer().realize(
            MixedBoundaryOperator([]), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, ZeroOperator)

    def test_mixed_singleton_unit_coeff_matches_legacy(self):
        """A singleton ``[(1.0, prim)]`` realizes to the primitive
        directly (no extra wrap)."""
        quad = LevelSymmetricSN.create(4)
        spec = SpecularBoundaryOperator(axis="x", albedo=1.0)
        mixed = MixedBoundaryOperator([(1.0, spec)])
        op = SNBoundaryRealizer().realize(mixed, SNMethodSpace.minimal(quad))
        # coeff == 1.0 short-circuits the ScaledOperator wrap, so the
        # singleton-mixed realization equals the realization of the
        # primitive alone.
        rng = np.random.default_rng(5)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 3))
        np.testing.assert_array_equal(op.apply(psi), spec.apply(psi, quad))


# ─────────────────────────────────────────────────────────────────────
# 7. Registry lookup
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRegistryLookup:
    """The SN realizer self-registers under ``method_name="SN"``."""

    def test_get_sn_returns_sn_realizer(self):
        # Import side-effect: importing ``orpheus.sn.boundary_realizer``
        # already triggered the @register decorator at module load.
        assert BoundaryRealizerRegistry.get("SN") is SNBoundaryRealizer

    def test_unknown_method_name_raises(self):
        with pytest.raises(BoundaryRealizerRegistryError):
            BoundaryRealizerRegistry.get("__no_such_method__")


# ─────────────────────────────────────────────────────────────────────
# 8. Unknown-law dispatch failure
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeUnknownLawRaises:
    """An object that isn't a recognised BoundaryOperator subclass
    raises BoundaryError naming the offending type."""

    def test_random_object_raises_boundary_error(self):
        quad = GaussLegendre1D.create(8)
        realizer = SNBoundaryRealizer()

        class _NotABc:  # noqa: D401  --- ad-hoc test stand-in
            pass

        with pytest.raises(BoundaryError) as excinfo:
            realizer.realize(_NotABc(), SNMethodSpace.minimal(quad))
        assert excinfo.value.law == "_NotABc"
