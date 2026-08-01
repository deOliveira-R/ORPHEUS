r"""Tests for SNBoundaryRealizer (Wave 5 / C5.3).

The realizer dispatches by ``isinstance(law, ...)`` to the Wave-0 /
Wave-1 primitives that realize each legacy
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` subclass as a
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
    AlbedoBoundary,
    BoundaryError,
    BoundaryRealizer,
    PeriodicBoundary,
    ReflectiveBoundary,
    VacuumAppliedToOutgoingTraceError,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    PeriodicWrapOperator,
    PermutationOperator,
    ScaledOperator,
    TensorProductOperator,
    ZeroOperator,
)
from orpheus.sn.boundary.angular import AngularAverageOperator
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import face_method_space, face_trace


# ─────────────────────────────────────────────────────────────────────
# 1. Vacuum
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeVacuum:
    r"""Vacuum realizes to the ZERO MAP :math:`\Gamma_+ \to \Gamma_-`.

    RE-POSED at campaign phase **B3.2** (C-1). The pre-B3.2 realization was an
    :class:`IncomingOrdinateMaskTensor` — a FULL-FACE projector that zeroed the
    inflow rows and *preserved* the outflow ones — and these tests asserted
    that pass-through. Two campaign phases documented the preserved rows as
    having "no consumer today"; B3.2 removes the question instead of answering
    it, because with the law's domain narrowed to :math:`\Gamma_+` those rows
    are not in the operator's domain at all.

    Vacuum's entire content is :math:`R = 0`, so the honest object is the zero
    map between the two half-traces. The assertions below therefore state:
    the image is zero, it lives on :math:`\Gamma_-`, and the operator is NOT an
    endomorphism of the face slot.
    """

    def test_lebedev17_xmin_is_the_zero_map_onto_gamma_minus(self):
        quad = Quadrature.lebedev(17)
        # xmin face: outward normal is -x, inflow is mu_x > 0 (into the domain).
        space = face_method_space(quad, face="xmin")
        inflow_indices = space.inflow_indices
        outflow_indices = space.outflow_indices
        rng = np.random.default_rng(42)
        psi = rng.uniform(0.5, 2.0, size=(quad.N, 4, 2))
        realizer = SNBoundaryRealizer()
        op = realizer.realize(VacuumInflow(), space)
        out = op.apply(psi[outflow_indices])
        # The whole image is zero, and it is Γ_- shaped.
        assert out.shape == (inflow_indices.size, 4, 2), (
            f"vacuum emitted {out.shape}; the narrowed codomain is Γ₋, i.e. "
            f"{(inflow_indices.size, 4, 2)}."
        )
        np.testing.assert_array_equal(out, 0.0)
        # …and it is NOT an endomorphism of the face slot. |Γ₊| == |Γ₋| here
        # (measured true for EVERY quadrature × face in the tree), so the shape
        # check above cannot distinguish Γ₊ → Γ₋ from Γ₊ → Γ₊ — vv Mode 12.
        # A full-face input must not produce a full-face image.
        assert op.apply(psi).shape[0] != quad.N, (
            "vacuum emitted a full-face image — it is still an endomorphism "
            "of the whole face slot, not the zero map Γ₊ → Γ₋."
        )

    def test_lebedev17_ymax_is_the_zero_map_onto_gamma_minus(self):
        quad = Quadrature.lebedev(17)
        # ymax face: outward normal is +y, inflow is mu_y < 0 (into the domain).
        space = face_method_space(
            quad, face="ymax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        # Independent spelling of the face's orientation: Ω·n_ymax = +mu_y, so
        # inflow ⟺ mu_y < 0. Pinning it here keeps the fixture honest about
        # WHICH face it is testing, rather than trusting the helper.
        np.testing.assert_array_equal(
            space.inflow_indices, np.flatnonzero(quad.mu_y < -1e-12),
        )
        psi = np.arange(quad.N * 3, dtype=float).reshape(quad.N, 3)
        op = SNBoundaryRealizer().realize(VacuumInflow(), space)
        out = op.apply(psi[space.outflow_indices])
        assert out.shape == (space.inflow_indices.size, 3)
        np.testing.assert_array_equal(out, 0.0)

    def test_vacuum_missing_inflow_indices_raises(self):
        """Without ``inflow_indices`` the realizer raises BoundaryError."""
        quad = Quadrature.lebedev(17)
        space = SNMethodSpace.minimal(quad)  # no inflow_indices
        realizer = SNBoundaryRealizer()
        with pytest.raises(BoundaryError) as excinfo:
            realizer.realize(VacuumInflow(), space)
        assert excinfo.value.law == "vacuum"

    def test_vacuum_missing_outflow_indices_raises(self):
        r"""B3.2 negative: a law's DOMAIN is :math:`\Gamma_+`, so a method
        space that cannot name the face's outflow ordinates cannot realize one.

        The sibling of the inflow negative above, and the reason
        ``SNMethodSpace.minimal`` no longer suffices for vacuum or reflective:
        face orientation is not derivable from a quadrature.
        """
        quad = Quadrature.lebedev(17)
        inflow_only = SNMethodSpace(
            quadrature=quad, face="xmin",
            inflow_indices=np.flatnonzero(quad.mu_x > 1e-12),
        )
        with pytest.raises(BoundaryError, match="outflow_indices"):
            SNBoundaryRealizer().realize(VacuumInflow(), inflow_only)

    def test_vacuum_returns_zero_operator_with_both_space_hooks(self):
        r"""The vacuum dispatch returns a :class:`ZeroOperator` carrying BOTH
        space hooks — forward emits the zero of :math:`\Gamma_-`, transpose the
        zero of :math:`\Gamma_+`.

        RE-POSED (C-1) from the pre-B3.2 assertion
        ``isinstance(op.ops[0], IncomingOrdinateMaskTensor)``. The type moved
        because the CONTRACT moved; the leg is kept (not deleted) because
        "which object realizes vacuum" is exactly what a future refactor could
        silently regress.

        The two hooks are asserted SEPARATELY and with different lengths where
        possible: a zero map between two DIFFERENT spaces must emit the zero of
        the space it lands in, and relying on the endomorphic ``0.0 * x`` echo
        would be wrong in principle and merely lucky in practice.
        """
        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(VacuumInflow(), space)
        assert isinstance(op, ZeroOperator)
        probe_plus = np.ones((space.outflow_indices.size, 3))
        probe_minus = np.ones((space.inflow_indices.size, 3))
        np.testing.assert_array_equal(
            op.apply(probe_plus), np.zeros((space.inflow_indices.size, 3)),
        )
        np.testing.assert_array_equal(
            op.apply_transpose(probe_minus),
            np.zeros((space.outflow_indices.size, 3)),
        )


# ─────────────────────────────────────────────────────────────────────
# 2. Reflective (specular)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeReflective:
    """Specular realizes to ``albedo * PermutationOperator(perm)``."""

    def test_specular_unit_albedo_lebedev_matches_hand_computed(self):
        r"""At α=1 the realized op MUST bit-match the hand-computed narrowed
        gather ``psi[reflection_index][inflow]``.

        RE-POSED at **B3.2** (C-1): the pre-B3.2 reference was the FULL-face
        gather ``psi[reflection_index]``. The narrowed reference is that SAME
        retired expression, RESTRICTED to :math:`\Gamma_-` — which is precisely
        the bit-identity claim the phase makes, so this leg doubles as the
        law-level half of it.

        Note what the reference does NOT use: ``to_local``. The expected value
        is built from ``reflection_index`` and the inflow index set alone, so
        an error in the production remap cannot cancel against the reference
        (they share no code above the numpy line — ``algebra-of-record``
        structural independence).
        """
        quad = Quadrature.lebedev(17)
        bc = ReflectiveBoundary(axis="x", albedo=1.0)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(7)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5, 3))
        expected = psi[quad.reflection_index("x")][space.inflow_indices]
        np.testing.assert_array_equal(
            op.apply(psi[space.outflow_indices]), expected,
        )

    def test_specular_unit_albedo_returns_tensor_product(self):
        r"""At α=1 the dispatch returns a :class:`TensorProductOperator`
        whose factors are ``(PermutationOperator, IdentityOperator)`` — with
        the permutation now on the **reduced** ordinate axis.

        Depth B step D-B+1 (2026-05-27): the first production tensor-network
        instance in ORPHEUS. Pre-D-B+1 the realiser returned a bare
        :class:`PermutationOperator(perm, axis=0)`; the implicit numpy
        broadcast across the group axis played the role of ``I_group``.
        Promoting to a typed TP makes the §16A.10
        ``B = G_patch ⊗ K_omega ⊗ K_g`` algebra type-visible.

        **B3.2** adds the size leg: the permutation's length is
        :math:`|\Gamma_+|`, not ``quad.N``. That single assertion is what makes
        the narrowing TYPE-VISIBLE on this object — without it, a revert to the
        full-``N`` permutation would keep every structural claim here green.
        """
        quad = Quadrature.level_symmetric(4)
        bc = ReflectiveBoundary(axis="x", albedo=1.0)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(bc, space)
        assert isinstance(op, TensorProductOperator)
        assert len(op.ops) == 2
        assert isinstance(op.ops[0], PermutationOperator)
        assert isinstance(op.ops[1], IdentityOperator)
        assert op.ops[0].perm.size == space.outflow_indices.size < quad.N, (
            f"the realized permutation has {op.ops[0].perm.size} entries; the "
            f"narrowed law acts on Γ₊ ({space.outflow_indices.size} of "
            f"{quad.N} ordinates)."
        )

    def test_specular_partial_albedo_matches_hand_computed(self):
        r"""At α=0.7 the realized op output equals
        ``0.7 * psi[reflection_index][inflow]`` bit-exactly.

        RE-POSED at B3.2 alongside its α=1 sibling — same retired reference,
        restricted to :math:`\Gamma_-`.
        """
        quad = Quadrature.lebedev(17)
        bc = ReflectiveBoundary(axis="x", albedo=0.7)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(9)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5))
        expected = 0.7 * psi[quad.reflection_index("x")][space.inflow_indices]
        np.testing.assert_array_equal(
            op.apply(psi[space.outflow_indices]), expected,
        )

    def test_specular_partial_albedo_returns_scaled_tensor_product(self):
        """At α≠1 the dispatch returns ``ScaledOperator(α, TP)`` where
        TP is the 2-factor :class:`TensorProductOperator` from D-B+1."""
        quad = Quadrature.level_symmetric(4)
        bc = ReflectiveBoundary(axis="x", albedo=0.5)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(bc, space)
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.5
        assert isinstance(op.op, TensorProductOperator)
        assert len(op.op.ops) == 2
        assert isinstance(op.op.ops[0], PermutationOperator)
        assert isinstance(op.op.ops[1], IdentityOperator)
        assert op.op.ops[0].perm.size == space.outflow_indices.size < quad.N

    def test_specular_narrowing_needs_a_bijection_between_half_traces(self):
        r"""B3.2 negative: a face where :math:`|\Gamma_-| \neq |\Gamma_+|` has
        NO specular realization, and says so.

        A specular mirror is a bijection between the two half-traces, so an
        asymmetric face is not a narrowing edge case — it is an ill-posed
        request. **[M]** No real quadrature × face pair in the tree produces
        one (every one measured ``|I| == |O|``), so the guard is unreachable
        from a mesh; a hand-built method space is the only way to exercise it,
        and that is exactly why this negative exists — otherwise the branch
        would ship with no test at all.
        """
        quad = Quadrature.gauss_legendre(8)
        trace = face_trace(quad)
        lopsided = SNMethodSpace(
            quadrature=quad, face="xmax",
            inflow_indices=trace.inflow_indices_for_face("xmax"),
            outflow_indices=trace.outflow_indices_for_face("xmax")[:-1],
        )
        with pytest.raises(BoundaryError, match=r"BIJECTION"):
            SNBoundaryRealizer().realize(
                ReflectiveBoundary(axis="x", albedo=1.0), lopsided,
            )


# ─────────────────────────────────────────────────────────────────────
# 3. White (Lambertian)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeWhite:
    """White realizes to ``albedo * AngularAverageOperator``.

    The body of :class:`AngularAverageOperator.from_quadrature` is
    lifted verbatim from :class:`WhiteBoundary.apply`, and
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
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        ref = AngularAverageOperator.from_quadrature(
            quad, axis="x", outward_sign=+1,
        )
        rng = np.random.default_rng(123)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 5, 3))
        np.testing.assert_array_equal(op.apply(psi), ref.apply(psi))

    def test_white_unit_albedo_returns_tensor_product(self):
        """At α=1 the dispatch returns a 2-factor
        :class:`TensorProductOperator` ``(AngularAverageOperator, IdentityOperator)``.

        Wave T step T.1 (2026-05-30): white BC lifted from a bare
        :class:`AngularAverageOperator` to the 2-factor TP shape
        introduced by D-B+1 for specular reflection.  The TP fold
        reduces to the bare ``AngularAverageOperator.apply`` (the
        :class:`IdentityOperator` factor returns ``x`` unchanged), so
        bit-identity at the apply level is preserved.
        """
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, TensorProductOperator)
        assert len(op.ops) == 2
        assert isinstance(op.ops[0], AngularAverageOperator)
        assert isinstance(op.ops[1], IdentityOperator)

    def test_white_partial_albedo_returns_scaled_tensor_product(self):
        """At α≠1 the dispatch returns ``ScaledOperator(α, TP)`` where
        TP is the 2-factor :class:`TensorProductOperator`."""
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3)
        op = SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.3
        assert isinstance(op.op, TensorProductOperator)
        assert len(op.op.ops) == 2
        assert isinstance(op.op.ops[0], AngularAverageOperator)
        assert isinstance(op.op.ops[1], IdentityOperator)

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
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3)
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
        :meth:`AngularAverageOperator.from_quadrature`.

        A 1-D GL quadrature has no outgoing ordinates on the z-axis, so
        the white-BC average operator cannot be built. The Wave-T
        angular-average lift replaced the legacy ``mu_z`` wording with a
        degenerate-face diagnostic; the test asserts that message.
        """
        quad = Quadrature.gauss_legendre(8)
        bc = WhiteBoundary(axis="z", outward_sign=+1, albedo=1.0)
        with pytest.raises(ValueError, match="degenerate for this face"):
            SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))


# ─────────────────────────────────────────────────────────────────────
# 4. Albedo (rank-1 isotropic scaling)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeAlbedo:
    """Albedo realizes to Zero / Identity / Scaled-Identity by α."""

    def test_albedo_zero_realizes_to_zero_operator(self):
        quad = Quadrature.gauss_legendre(8)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(0.0), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, ZeroOperator)
        psi = np.arange(quad.N * 4, dtype=float).reshape(quad.N, 4)
        np.testing.assert_array_equal(op.apply(psi), 0.0)

    def test_albedo_one_realizes_to_identity_operator(self):
        quad = Quadrature.gauss_legendre(8)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(1.0), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, IdentityOperator)
        psi = np.arange(quad.N * 4, dtype=float).reshape(quad.N, 4)
        np.testing.assert_array_equal(op.apply(psi), psi)

    def test_albedo_half_realizes_to_scaled_tensor_product(self):
        """At α=0.5 the dispatch returns
        ``ScaledOperator(0.5, IdentityOperator() & IdentityOperator())``.

        Wave T step T.1 (2026-05-30): the inner identity lifts to a
        2-factor TP (I & I), making the §16A.10 algebra type-visible
        while the apply remains a no-op (both factors return ``x``
        unchanged; the ``ScaledOperator`` wrapper supplies the α
        multiplication).
        """
        quad = Quadrature.gauss_legendre(8)
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(0.5), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.5
        assert isinstance(op.op, TensorProductOperator)
        assert len(op.op.ops) == 2
        assert isinstance(op.op.ops[0], IdentityOperator)
        assert isinstance(op.op.ops[1], IdentityOperator)
        psi = np.arange(quad.N * 4, dtype=float).reshape(quad.N, 4)
        np.testing.assert_array_equal(op.apply(psi), 0.5 * psi)


# ─────────────────────────────────────────────────────────────────────
# 5. Periodic
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizePeriodic:
    """Periodic realizes to :class:`PeriodicWrapOperator`."""

    def test_periodic_returns_tensor_product_and_passes_through(self):
        """Periodic realizes to a 2-factor :class:`TensorProductOperator`
        ``(PeriodicWrapOperator, IdentityOperator)`` whose apply passes
        through values unchanged.

        Wave T step T.1 (2026-05-30): the bare
        :class:`PeriodicWrapOperator` (identity-with-copy body) is now
        lifted to a 2-factor TP, mirroring the D-B+1 specular pattern.
        Compare values, NOT identity: the TP fold returns a copy from
        the first factor and the :class:`IdentityOperator` second
        factor returns it unchanged.
        """
        quad = Quadrature.lebedev(17)
        op = SNBoundaryRealizer().realize(
            PeriodicBoundary(), SNMethodSpace.minimal(quad),
        )
        assert isinstance(op, TensorProductOperator)
        assert len(op.ops) == 2
        assert isinstance(op.ops[0], PeriodicWrapOperator)
        assert isinstance(op.ops[1], IdentityOperator)
        rng = np.random.default_rng(3)
        psi = rng.uniform(-1.0, 1.0, size=(quad.N, 5))
        np.testing.assert_array_equal(op.apply(psi), psi)


# ─────────────────────────────────────────────────────────────────────
# 6. Mixed (rank-N composition) — Wave 11 removal
#
# The ``MixedBoundaryOperator`` composer was removed in Wave 11; the
# realizer no longer dispatches on a "mixed" type.  Rank-N compositions
# are now built via Wave-0 ``OperatorSum``/``ScaledOperator`` algebra
# over already-realised leaves; the tree-walking helper
# :func:`orpheus.geometry.boundary.realize_recursively` realises a
# ``BoundaryTraceLaw``-rooted expression by recursing through the
# Wave-0 composers (its tests live in
# ``tests/geometry/test_law_composition.py``).
# ─────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────
# 7. Realizer identity
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizerIdentity:
    """The structural pins that replaced the registry-lookup pins when
    #290 P7b dissolved ``BoundaryRealizerRegistry``."""

    def test_method_name_attribute(self):
        assert SNBoundaryRealizer.method_name == "SN"

    def test_conforms_to_the_boundary_realizer_protocol(self):
        # The Protocol the walker and ``SNMesh.realize_boundary_law``
        # dispatch through.
        assert isinstance(SNBoundaryRealizer(), BoundaryRealizer)


# ─────────────────────────────────────────────────────────────────────
# 8. Unknown-law dispatch failure
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeUnknownLawRaises:
    """An object that isn't a recognised BoundaryTraceLaw subclass
    raises BoundaryError naming the offending type."""

    def test_random_object_raises_boundary_error(self):
        quad = Quadrature.gauss_legendre(8)
        realizer = SNBoundaryRealizer()

        class _NotABc:  # noqa: D401  --- ad-hoc test stand-in
            pass

        with pytest.raises(BoundaryError) as excinfo:
            realizer.realize(_NotABc(), SNMethodSpace.minimal(quad))
        assert excinfo.value.law == "_NotABc"


# ─────────────────────────────────────────────────────────────────────
# 9. Vacuum trace-orientation guard (#52, ERR-041)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-041")
class TestVacuumTraceOrientationGuard:
    """The realizer-seam ERR-041 guard: claimed ``inflow_indices``
    cross-checked against the orientation the FACE NAME alone implies.

    The catalog's threat model verbatim: "a realiser swapped the
    face's (inflow, outflow) annotation by mistake" — the resulting
    vacuum operator zeroes the OUTGOING flux the sweep just computed,
    while every shape-level test stays green (``apply`` returns a
    same-shaped array regardless of orientation). The face name and
    the index array are two independently-supplied encodings of one
    orientation; the guard reads the face's signed projection
    Ω·n̂ through the single
    :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`
    primitive and refuses any claimed inflow index that is outgoing.
    """

    def test_swapped_annotation_raises(self):
        """The full annotation swap: the OUTGOING set of xmax
        (μ_x > 0) claimed as its inflow set."""
        quad = Quadrature.gauss_legendre(8)
        swapped = SNMethodSpace(
            quadrature=quad, face="xmax",
            inflow_indices=np.flatnonzero(quad.mu_x > 0),
        )
        with pytest.raises(VacuumAppliedToOutgoingTraceError):
            SNBoundaryRealizer().realize(VacuumInflow(), swapped)

    def test_single_contaminated_index_raises(self):
        """The sharper form: ONE outgoing ordinate hidden inside an
        otherwise-correct inflow set (a per-ordinate indexing slip,
        not a wholesale swap)."""
        quad = Quadrature.gauss_legendre(8)
        correct = np.flatnonzero(quad.mu_x < 0)
        one_outgoing = np.flatnonzero(quad.mu_x > 0)[:1]
        contaminated = np.concatenate([correct, one_outgoing])
        space = SNMethodSpace(
            quadrature=quad, face="xmax", inflow_indices=contaminated,
        )
        with pytest.raises(VacuumAppliedToOutgoingTraceError):
            SNBoundaryRealizer().realize(VacuumInflow(), space)

    def test_correct_orientation_realizes(self):
        """Positive control: the true inflow set of each face realizes
        clean through the SAME guarded path, both face sides.

        B3.2: the realized object is now the ``Γ₊ → Γ₋`` zero map, so the
        type assertion moves from ``TensorProductOperator`` to
        :class:`ZeroOperator`. The CLAIM — "the correctly-oriented face
        realizes silently" — is unchanged, which is the point of the leg.
        """
        quad = Quadrature.gauss_legendre(8)
        # The face-side pin the pre-B3.2 spelling carried inline: xmax has
        # outward normal +x̂ so Ω·n = +μ_x and inflow ⟺ μ_x < 0; xmin is the
        # mirror image. Written out per face rather than as a signed formula —
        # a sign convention compressed into an expression is exactly the kind
        # of thing this test exists to catch.
        for face, inflow in (
            ("xmax", np.flatnonzero(quad.mu_x < -1e-12)),
            ("xmin", np.flatnonzero(quad.mu_x > +1e-12)),
        ):
            space = face_method_space(quad, face=face)
            np.testing.assert_array_equal(space.inflow_indices, inflow)
            op = SNBoundaryRealizer().realize(VacuumInflow(), space)
            assert isinstance(op, ZeroOperator)

    def test_faceless_space_carries_no_orientation_truth(self):
        r"""A method space with hand-supplied half-traces and NO face name has
        no independent orientation encoding — the guard cannot fire there
        (documented escape; the canonical ``SNMesh.realize_boundary_law`` path
        always carries a face). The realize must succeed, wrong indices and
        all: the caller owns the claim.

        B3.2 sharpens the fixture rather than the claim. Since a law's DOMAIN
        is :math:`\Gamma_+`, a faceless space must now supply BOTH index sets
        by hand — but supplying them is not supplying an *orientation truth*,
        because nothing cross-checks them against a face normal. So the escape
        survives verbatim: the deliberately-swapped inflow set (the OUTGOING
        ordinates of a nominal xmax face) still realizes silently.
        """
        quad = Quadrature.gauss_legendre(8)
        faceless = SNMethodSpace(
            quadrature=quad,
            inflow_indices=np.flatnonzero(quad.mu_x > 0),
            outflow_indices=np.flatnonzero(quad.mu_x < 0),
        )
        op = SNBoundaryRealizer().realize(VacuumInflow(), faceless)
        assert isinstance(op, ZeroOperator)

    def test_unknown_face_name_fails_loud(self):
        """A face name outside the ``{axis}{min|max}`` convention is
        uncertifiable orientation data — the guard's primitive raises
        a loud ``ValueError`` naming the valid faces rather than
        silently skipping the check (the pre-#52 legacy strings
        ``"left"``/``"right"`` migrate, not linger)."""
        quad = Quadrature.gauss_legendre(8)
        legacy = SNMethodSpace(
            quadrature=quad, face="right",
            inflow_indices=np.flatnonzero(quad.mu_x < 0),
        )
        with pytest.raises(ValueError, match="Unknown face name"):
            SNBoundaryRealizer().realize(VacuumInflow(), legacy)
