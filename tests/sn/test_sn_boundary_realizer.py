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

    def test_specular_unit_albedo_lebedev_matches_hand_computed(self):
        """At α=1 the realized op MUST bit-match the hand-computed
        ``psi[reflection_index]`` gather.

        Issue #186 (B3 + β2) rewrite: descriptors are no longer
        callable, so the previous ``bc.apply(psi, quad)`` half cannot
        be evaluated. The hand-computed gather is strictly stronger:
        we assert the realised op's output against the structural
        definition of specular reflection, not against another
        implementation.
        """
        quad = LebedevSphere.create(17)
        bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
        space = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(7)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5, 3))
        expected = psi[quad.reflection_index("x")]
        np.testing.assert_array_equal(op.apply(psi), expected)

    def test_specular_unit_albedo_returns_bare_permutation(self):
        """At α=1 the dispatch returns the bare PermutationOperator
        (no ScaledOperator wrap) — preserves bit-identity vs legacy."""
        quad = LevelSymmetricSN.create(4)
        bc = SpecularBoundaryOperator(axis="x", albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, PermutationOperator)

    def test_specular_partial_albedo_matches_hand_computed(self):
        """At α=0.7 the realized op output equals ``0.7 * psi[reflection_index]``
        bit-exactly.

        Issue #186 (B3 + β2) rewrite: hand-computed RHS instead of
        the legacy ``bc.apply``. ``ScaledOperator(0.7, P).apply(psi)``
        first gathers ``psi[perm]`` (via the inner PermutationOperator)
        then multiplies by 0.7 — identical to ``0.7 * psi[perm]``.
        """
        quad = LebedevSphere.create(17)
        bc = SpecularBoundaryOperator(axis="x", albedo=0.7)
        space = SNMethodSpace.minimal(quad)
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(9)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5))
        expected = 0.7 * psi[quad.reflection_index("x")]
        np.testing.assert_array_equal(op.apply(psi), expected)

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

    def test_white_unit_albedo_lebedev_matches_angular_average_operator(self):
        """At α=1 the realized op MUST bit-match a bare
        :class:`AngularAverageOperator.from_quadrature` output.

        Issue #186 (B3 + β2) rewrite: the realiser's fast path returns
        the bare :class:`AngularAverageOperator` at α=1, so this test
        pins the dispatch chain rather than legacy equivalence. The
        :class:`AngularAverageOperator` body is itself verified by the
        Wave-1 cosine-weighted-current and self-adjointness tests in
        ``tests/sn/test_angular_average_operator.py`` — those are the
        structurally-independent references.
        """
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        ref = AngularAverageOperator.from_quadrature(
            quad, axis="x", outward_sign=+1,
        )
        rng = np.random.default_rng(123)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 5, 3))
        np.testing.assert_array_equal(op.apply(psi), ref.apply(psi))

    def test_white_unit_albedo_returns_bare_angular_average(self):
        """At α=1 the dispatch returns the bare AngularAverageOperator."""
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, AngularAverageOperator)

    def test_white_partial_albedo_matches_albedo_times_angular_average(self):
        """At α=0.3 the realized op output equals
        ``0.3 * AngularAverageOperator(...).apply(psi)`` within
        ``nulp=4``.

        Issue #186 (B3 + β2) rewrite: hand-computed RHS via the
        Wave-1 :class:`AngularAverageOperator` primitive (the
        structurally-independent reference; its body is verified in
        ``tests/sn/test_angular_average_operator.py``).
        ``ScaledOperator(0.3, ...).apply`` multiplies the inner
        output by 0.3 from the outside; the hand-computed
        ``0.3 * ref.apply(psi)`` may differ by one ULP under
        IEEE-754 FP-non-associativity. ``nulp=4`` is the safe floor
        per the ``algebra-of-record`` bit-identity discipline.
        """
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=0.3)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        ref = AngularAverageOperator.from_quadrature(
            quad, axis="x", outward_sign=+1,
        )
        rng = np.random.default_rng(11)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 4))
        expected = 0.3 * ref.apply(psi)
        np.testing.assert_array_almost_equal_nulp(
            op.apply(psi), expected, nulp=4,
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
# 6. Mixed (rank-N composition) — Wave 11 removal
#
# The ``MixedBoundaryOperator`` composer was removed in Wave 11; the
# realizer no longer dispatches on a "mixed" type.  Rank-N compositions
# are now built via Wave-0 ``OperatorSum``/``ScaledOperator`` algebra
# over already-realised leaves; the tree-walking helper
# :func:`orpheus.sn.boundary_realize.realize_recursively` realises a
# ``BoundaryTraceLaw``-rooted expression by recursing through the
# Wave-0 composers (its tests live in
# ``tests/sn/test_boundary_realize.py``).
# ─────────────────────────────────────────────────────────────────────


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
