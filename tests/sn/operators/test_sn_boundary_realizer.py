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
    IsotropicReturn,
    PeriodicBoundary,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumAppliedToOutgoingTraceError,
    VacuumInflow,
    WhiteBoundary,
)
from orpheus.numerics.operator import (
    IdentityOperator,
    PermutationOperator,
    ScaledOperator,
    TensorProductOperator,
    ZeroOperator,
)
from orpheus.sn.boundary.realizer import SNBoundaryRealizer
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import face_method_space, face_trace


# ─────────────────────────────────────────────────────────────────────
# 1. Vacuum
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1

def assert_lambertian_kernel(kernel, space, quad) -> None:
    r"""The surviving claims about a realized Lambertian, posed BEHAVIOURALLY.

    Re-posed at G6.3 step 3b (#330), where the welded ``AngularAverageOperator``
    was retired onto the chain ``IsotropicEmission @ PartialCurrent``. The old
    spelling pinned the concrete CLASS and read ``n_outflow`` / ``n_inflow`` off
    it; neither survives a factorisation, and neither was the claim.

    ⭐ **The load-bearing leg is the SIZE one**, and B3.4a's docstring says why:
    *"without it a revert to the full-face kernel would keep every structural
    claim here green."* It is re-posed as a **behavioural** pair — the kernel
    ACCEPTS ``|Γ₊|`` rows and REFUSES ``quad.N`` — which is strictly stronger
    than reading an attribute, because it holds for any future spelling and
    cannot be satisfied by an operator that merely reports the right number.
    """
    n_out = int(space.outflow_indices.size)
    n_in = int(space.inflow_indices.size)
    assert n_out < quad.N, (
        "fixture no longer exercises the narrowing — Γ₊ is the whole face"
    )

    out = np.asarray(kernel.apply(np.ones((n_out,))))
    assert out.shape == (n_in,), (
        f"the realized Lambertian emits {out.shape} on Γ₋ of size {n_in}"
    )
    with pytest.raises(ValueError):
        kernel.apply(np.ones((quad.N,)))


def lambertian_reference(quad, axis, sign, psi_out, space):
    r"""An INDEPENDENT closed form for the Lambertian, in explicit numpy.

    .. math::

        (R\psi)(\Omega) = \frac{\sum_{\Gamma_+} w\,|\Omega'\cdot\hat n|\,
        \psi(\Omega')}{\sum_{\Gamma_+} w\,|\Omega'\cdot\hat n|},
        \qquad \Omega \in \Gamma_-.

    ⭐ **Why this exists rather than a production operator.** These gates used
    ``AngularAverageOperator.from_quadrature`` as their reference. G6.3 step 3b
    retires that class onto the chain ``IsotropicEmission @ PartialCurrent`` —
    which is *what the realizer now builds*. Re-pointing the reference at the
    chain would compare the realizer's output **with itself through a wrapper**:
    green forever, still named a bit-identity check, and unable to detect the
    drift its docstring advertises. That is the documented rewire-demotion trap.

    So the reference is rebuilt from the quadrature's own ``weights`` and
    direction cosines with no production operator in the path — strictly MORE
    independent than the class it replaces, which was itself production code.
    """
    mu_n = np.asarray(getattr(quad, f"mu_{axis}"), dtype=float) * sign
    out_i = np.asarray(space.outflow_indices)
    cos_w = np.asarray(quad.weights, dtype=float)[out_i] * mu_n[out_i]
    psi_out = np.asarray(psi_out)
    contracted = (
        cos_w.reshape((-1,) + (1,) * (psi_out.ndim - 1)) * psi_out
    ).sum(axis=0)
    n_in = int(np.asarray(space.inflow_indices).size)
    return np.broadcast_to(
        (contracted / cos_w.sum())[None, ...],
        (n_in,) + contracted.shape,
    ).copy()

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
    r"""White realizes to ``albedo * AngularAverageOperator``, narrowed at
    **B3.4a** to :math:`\Gamma_+ \to \Gamma_-`.

    MIGRATED at B3.4a. Every leg below used to build its method space with
    ``SNMethodSpace.minimal(quad)`` — legal while white was a full-face
    endomorphism that read only the quadrature. Since the law's DOMAIN is the
    outflow half-trace of a PARTICULAR face, a faceless space can no longer
    name it and the realizer refuses (see
    :meth:`test_white_on_a_faceless_space_is_refused`); the fixture is now
    :func:`~tests.sn._test_helpers.face_method_space`, whose face must MATCH
    the law's declared orientation (:class:`TestWhiteOrientationGuard`).

    The :class:`AngularAverageOperator` body is verified independently by the
    current-conservation and reciprocity gates in
    ``tests/sn/operators/test_angular_average_operator.py``; these legs pin
    the realizer's dispatch chain against it.
    """

    def test_white_unit_albedo_lebedev_matches_angular_average_operator(self):
        r"""At α=1 the realized op MUST bit-match a bare
        :class:`AngularAverageOperator.from_quadrature` output.

        Issue #186 (B3 + β2) rewrite: the realiser's fast path returns
        the bare :class:`AngularAverageOperator` at α=1, so this test
        pins the dispatch chain rather than legacy equivalence. **B3.4a**
        restricts the probe to :math:`\Gamma_+` — the claim (dispatch adds
        nothing) is unchanged; only the domain it is stated on moved.
        """
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(123)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 5, 3))
        psi_out = psi[space.outflow_indices]
        got = op.apply(psi_out)
        assert got.shape == (space.inflow_indices.size, 5, 3), (
            f"white emitted {got.shape}; the narrowed codomain is Γ₋, i.e. "
            f"{(space.inflow_indices.size, 5, 3)}."
        )
        np.testing.assert_array_equal(
            got, lambertian_reference(quad, "x", +1, psi_out, space),
        )

    def test_white_unit_albedo_returns_tensor_product(self):
        r"""At α=1 the dispatch returns a 2-factor
        :class:`TensorProductOperator` ``(the Lambertian chain, IdentityOperator)``.

        Wave T step T.1 (2026-05-30): white BC lifted from a bare
        :class:`AngularAverageOperator` to the 2-factor TP shape
        introduced by D-B+1 for specular reflection.  The TP fold
        reduces to the bare ``AngularAverageOperator.apply`` (the
        :class:`IdentityOperator` factor returns ``x`` unchanged), so
        bit-identity at the apply level is preserved.

        **B3.4a** adds the size leg, mirroring the specular sibling: the
        kernel's DOMAIN is :math:`|\Gamma_+|`, not ``quad.N``. Without it a
        revert to the full-face kernel would keep every structural claim
        here green.
        """
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        assert isinstance(op, TensorProductOperator)
        assert len(op.ops) == 2
        assert isinstance(op.ops[1], IdentityOperator)
        assert_lambertian_kernel(op.ops[0], space, quad)

    def test_white_partial_albedo_returns_scaled_tensor_product(self):
        """At α≠1 the dispatch returns ``ScaledOperator(α, TP)`` where
        TP is the 2-factor :class:`TensorProductOperator`."""
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.3
        assert isinstance(op.op, TensorProductOperator)
        assert len(op.op.ops) == 2
        assert isinstance(op.op.ops[1], IdentityOperator)
        assert_lambertian_kernel(op.op.ops[0], space, quad)

    def test_white_partial_albedo_matches_albedo_times_angular_average(self):
        r"""At α=0.3 the realized op output equals
        ``0.3 * AngularAverageOperator(...).apply(γ₊ψ)`` within ``nulp=4``.

        Issue #186 (B3 + β2) rewrite: hand-computed RHS via the
        Wave-1 :class:`AngularAverageOperator` primitive (the
        structurally-independent reference; its body is verified in
        ``tests/sn/operators/test_angular_average_operator.py``).
        ``ScaledOperator(0.3, ...).apply`` multiplies the inner
        output by 0.3 from the outside; the hand-computed
        ``0.3 * ref.apply(psi)`` may differ by one ULP under
        IEEE-754 FP-non-associativity. ``nulp=4`` is the safe floor
        per the ``algebra-of-record`` bit-identity discipline.
        """
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=0.3)
        space = face_method_space(quad, face="xmax")
        op = SNBoundaryRealizer().realize(bc, space)
        rng = np.random.default_rng(11)
        psi = rng.uniform(0.0, 2.0, size=(quad.N, 4))
        psi_out = psi[space.outflow_indices]
        expected = 0.3 * lambertian_reference(quad, "x", +1, psi_out, space)
        np.testing.assert_array_almost_equal_nulp(
            op.apply(psi_out), expected, nulp=4,
        )

    def test_white_z_axis_on_gauss_legendre_raises(self):
        r"""A ``z``-axis white law on a 1-D GL quadrature is a RANK MISMATCH
        between the face the law names and the cubature, and says so.

        RE-POSED at **B3.4a**. Pre-B3.4a the raise came from
        :meth:`AngularAverageOperator.from_quadrature`'s own "no outgoing
        ordinates / degenerate for this face" guard, reached because the
        operator classified its own hemisphere with a strict compare against
        ``mu_z == zeros(N)``. It now comes from
        :func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`
        — the single face-name → signed-projection primitive, which the
        orientation cross-check consults first — and names the defect
        directly. The message is strictly more specific, so the pattern
        tightens.

        The method space is ``xmax`` (a face this quadrature genuinely has):
        the defect under test is the LAW's degenerate axis, and a ``zmax``
        method space could not be stood up on a GL quadrature at all.
        """
        quad = Quadrature.gauss_legendre(8)
        bc = WhiteBoundary(axis="z", outward_sign=+1, albedo=1.0)
        with pytest.raises(ValueError, match="requires genuine mu_z"):
            SNBoundaryRealizer().realize(bc, face_method_space(quad, face="xmax"))

    def test_white_on_a_faceless_space_is_refused(self):
        r"""B3.4a negative: white's DOMAIN is :math:`\Gamma_+`, so a
        quadrature-only method space cannot realize it.

        The sibling of the vacuum/reflective negatives above, and the reason
        every leg in this class migrated off ``SNMethodSpace.minimal``.
        NET-NEW: the guard reached white for the first time at B3.4a and
        shipped with no negative.
        """
        quad = Quadrature.lebedev(17)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
        with pytest.raises(BoundaryError, match="outflow_indices") as excinfo:
            SNBoundaryRealizer().realize(bc, SNMethodSpace.minimal(quad))
        assert excinfo.value.law == "white"


# ─────────────────────────────────────────────────────────────────────
# 3b. White orientation cross-check (B3.4a, the ERR-041 pattern)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestWhiteOrientationGuard:
    r"""The law's declared ``axis``/``outward_sign`` is cross-checked against
    the face it is installed on.

    NET-NEW at **B3.4a**. :class:`WhiteBoundary` carries its OWN orientation
    while the method space independently names the face — two encodings of
    one datum, and until B3.4a nothing compared them. A white law declared
    for ``xmax`` and installed on ``xmin`` averaged over that face's INFLOW
    hemisphere and reported nothing: the Lambertian's shape is
    orientation-blind, so every shape- and type-level assertion stayed green.
    Same threat model as ERR-041 (the vacuum trace-orientation swap), same
    guard shape; the comparison is on index SETS rather than sizes because
    `[M]` ``|Γ₊| == |Γ₋|`` on every quadrature × face pair in the tree, which
    makes a size comparison Mode-12 blind.

    The guard is green by construction on the canonical
    ``SNMesh.realize_boundary_law`` path (both encodings derive from the same
    face label); it bites on hand-built method spaces and on a mis-declared
    law — which is exactly where the tree's white fixtures live.
    """

    def test_law_declared_for_the_opposite_face_raises(self):
        quad = Quadrature.gauss_legendre(8)
        bc = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)  # xmax
        with pytest.raises(BoundaryError, match=r"Γ₊") as excinfo:
            SNBoundaryRealizer().realize(bc, face_method_space(quad, face="xmin"))
        assert excinfo.value.law == "white"

    def test_law_declared_for_a_different_axis_raises(self):
        quad = Quadrature.level_symmetric(4)
        bc = WhiteBoundary(axis="y", outward_sign=+1, albedo=1.0)  # ymax
        space = face_method_space(
            quad, face="xmax", faces=("xmin", "xmax", "ymin", "ymax"),
        )
        with pytest.raises(BoundaryError, match=r"Γ₊") as excinfo:
            SNBoundaryRealizer().realize(bc, space)
        assert excinfo.value.law == "white"

    def test_matching_orientation_realizes_on_every_face(self):
        r"""Positive control: each of the four 2-D faces realizes silently
        when the law names it, and the realized kernel's domain is that
        face's own :math:`\Gamma_+`.

        Without this leg the negatives above are unattributable — a realize
        that raised for ANY reason would satisfy them.
        """
        quad = Quadrature.level_symmetric(4)
        faces = ("xmin", "xmax", "ymin", "ymax")
        for face in faces:
            axis, sign = face[0], (+1 if face.endswith("max") else -1)
            space = face_method_space(quad, face=face, faces=faces)
            op = SNBoundaryRealizer().realize(
                WhiteBoundary(axis=axis, outward_sign=sign, albedo=1.0), space,
            )
            assert isinstance(op, TensorProductOperator), face
            assert_lambertian_kernel(op.ops[0], space, quad)


# ─────────────────────────────────────────────────────────────────────
# 4. Albedo (rank-1 isotropic scaling)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeAlbedo:
    r"""Albedo realizes through its re-emission CLOSURE — or is refused.

    **RE-POSED at campaign phase B3.4b.** These rows asserted the pre-B3.4b
    dispatch: a BARE ``AlbedoBoundary(α)`` on a FACELESS ``SNMethodSpace``
    returning ``ZeroOperator`` / ``IdentityOperator`` / ``ScaledOperator(α,
    I & I)`` — three full-face ENDOMORPHISMS. The composite ``ι₋ ∘ law ∘ γ₊``
    then read them positionally: MEASURED, inflow row ``j`` received ``α``
    times outflow row ``j``, a pairing between two independently-sorted index
    sets carrying no geometry at all. (And the coincidence was
    quadrature-dependent — on ``product(2,4)`` and ``level_symmetric(6)`` that
    positional pairing happens to EQUAL the specular one, on
    ``gauss_legendre`` and ``lebedev(17)`` it does not.)

    So the bare spelling is now REFUSED, and the three shape claims re-pose
    onto the COMPLETED law, whose closure names the pairing.

    Division of labour with ``tests/geometry/test_reemission_closure.py``
    --------------------------------------------------------------------

    That module owns the refusal's MESSAGE contract (both completions named,
    the array-position defect named, attribution) and the ``≡`` theorems.
    This class deliberately asserts something it does not: that the refusal
    is reached **through the realizer's albedo arm specifically**, i.e. before
    and independently of the FACE-data guard — which is the one thing a
    faceless-vs-face-ful distinction can pin, and the reason the old rows
    lived here in the first place.
    """

    #: A face-ful space; since B3.2 every narrowed law needs both half-traces.
    @staticmethod
    def _space(n_ord: int = 8):
        return face_method_space(Quadrature.gauss_legendre(n_ord), face="xmax")

    def test_albedo_zero_realizes_to_the_narrowed_zero_map(self):
        r"""α=0 → ``ZeroOperator``, but now :math:`\Gamma_+ \to \Gamma_-`.

        The claim SURVIVES the narrowing and changes shape: pre-B3.4b this was
        an endomorphism of the whole face slot; now it emits ``|Γ₋|`` rows.
        Note this branch could not be reached AT ALL before B3.4b for the
        sibling laws — ``ScaledOperator`` refuses a zero scalar, so
        ``ReflectiveBoundary(axis, 0.0)`` and ``WhiteBoundary(…, 0.0)`` were
        legal laws that died in the numerics layer. The fold into
        ``_attenuated_kernel_operator`` fixed all four routes at once.
        """
        space = self._space()
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(0.0, SpecularReturn(axis="x")), space,
        )
        assert isinstance(op, ZeroOperator)
        psi = np.arange(space.outflow_indices.size * 4, dtype=float).reshape(
            space.outflow_indices.size, 4,
        )
        image = op.apply(psi)
        assert image.shape == (space.inflow_indices.size, 4)
        np.testing.assert_array_equal(image, 0.0)

    def test_albedo_one_realizes_to_the_bare_kernel_tensor_product(self):
        r"""α=1 → the bare 2-factor TP of the closure's kernel.

        Pre-B3.4b this returned a bare ``IdentityOperator`` — the spelling
        that made the law an endomorphism. The α=1 fast path survives (no
        ``ScaledOperator`` wrapper) but the inner factor is now the closure's
        angular kernel, so the operator carries the pairing instead of
        asserting there is none.
        """
        space = self._space()
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(1.0, SpecularReturn(axis="x")), space,
        )
        assert isinstance(op, TensorProductOperator)
        assert len(op.ops) == 2
        assert isinstance(op.ops[0], PermutationOperator)
        assert isinstance(op.ops[1], IdentityOperator)
        psi = np.arange(space.outflow_indices.size * 4, dtype=float).reshape(
            space.outflow_indices.size, 4,
        )
        assert op.apply(psi).shape == (space.inflow_indices.size, 4)

    def test_albedo_half_realizes_to_scaled_tensor_product(self):
        r"""α ∉ {0,1} → ``ScaledOperator(α, kernel & I)``.

        The SHAPE claim the pre-B3.4b row made is intact — a ``ScaledOperator``
        wrapping a 2-factor TP, the §16A.10 ``B = G_patch ⊗ K_ω ⊗ K_g``
        algebra kept type-visible. What changed is the inner factor: an
        ``IdentityOperator`` (no pairing) became the closure's angular kernel.
        """
        space = self._space()
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(0.5, SpecularReturn(axis="x")), space,
        )
        assert isinstance(op, ScaledOperator)
        assert op.scalar == 0.5
        assert isinstance(op.op, TensorProductOperator)
        assert len(op.op.ops) == 2
        assert isinstance(op.op.ops[0], PermutationOperator)
        assert isinstance(op.op.ops[1], IdentityOperator)

    def test_the_diffuse_closure_realizes_to_the_lambertian_kernel(self):
        """The OTHER closure reaches the OTHER shared body.

        Without this row the class would pin only the specular route and a
        dispatch that ignored ``IsotropicReturn`` entirely would pass.
        """
        space = self._space()
        op = SNBoundaryRealizer().realize(
            AlbedoBoundary(0.5, IsotropicReturn(axis="x", outward_sign=+1)),
            space,
        )
        assert isinstance(op, ScaledOperator)
        assert isinstance(op.op, TensorProductOperator)

    def test_the_bare_spelling_is_refused_by_the_ALBEDO_ARM_not_the_face_guard(
        self,
    ):
        r"""The refusal is reached on a FACE-FUL space — so it is the albedo
        arm's own, not a fallout of the missing-face guard.

        This is the leg that justifies keeping a refusal test HERE alongside
        the message-contract one next door (Pattern 2: two gates, two
        different claims). ``_outflow_restriction`` also raises
        ``BoundaryError`` with ``law="albedo"``, so a refusal probed on a
        FACELESS space cannot tell the two apart — and that weaker probe is
        exactly what a reader migrating the old rows would reach for, since
        the old rows used ``SNMethodSpace.minimal``.
        """
        space = self._space()          # face-ful: the face guard CANNOT fire
        with pytest.raises(BoundaryError) as exc:
            SNBoundaryRealizer().realize(AlbedoBoundary(0.5), space)
        assert exc.value.law == "albedo"
        assert "outflow_indices" not in str(exc.value), (
            "the raise came from the missing-face guard, not the "
            "angular-resolution refusal — this space HAS both half-traces, so "
            "that guard firing would mean the arm never ran."
        )
        # CONTROL: the faceless space raises too, but for the OTHER reason —
        # which is what makes the assertion above discriminating rather than
        # decorative.
        with pytest.raises(BoundaryError) as faceless:
            SNBoundaryRealizer().realize(
                AlbedoBoundary(0.5, SpecularReturn(axis="x")),
                SNMethodSpace.minimal(Quadrature.gauss_legendre(8)),
            )
        assert "outflow_indices" in str(faceless.value)


# ─────────────────────────────────────────────────────────────────────
# 5. Periodic
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizePeriodic:
    """Periodic realizes to the IDENTITY between the partner's Γ₊ and this
    face's Γ₋ — the crossing lives in the channel, not in the operator."""

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

        **Re-posed at B3.4c** — the same migration B3.2 made for vacuum /
        reflective and B3.4a for white / prescribed inflow, now that periodic
        is the LAST law to narrow. The probe used to run on
        ``SNMethodSpace.minimal(quad)``, a quadrature-only space with no face;
        periodic's domain is a *different face's* :math:`\\Gamma_+`, so a
        space that cannot name the installation face cannot name the partner
        either, and the realizer refuses instead of guessing.

        The pass-through claim is unchanged and is exactly what the identity
        body means — but it is now the identity on the **narrowed** trace, so
        the probe is :math:`\\Gamma_+`-shaped rather than full-face. That the
        two shapes differ (49 of 74 ordinates on ``lebedev(17)``) is what
        makes this a statement about the narrowed law rather than about an
        endomorphism that happens to be the identity.
        """
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmin")
        op = SNBoundaryRealizer().realize(PeriodicBoundary(axis="x"), space)
        assert isinstance(op, TensorProductOperator)
        assert len(op.ops) == 2
        assert all(isinstance(f, IdentityOperator) for f in op.ops)
        n_domain = np.asarray(space.outflow_indices).size
        assert n_domain != quad.N, (
            "the probe must be narrower than the full face, or the shape "
            "cannot distinguish a narrowed law from an endomorphism"
        )
        rng = np.random.default_rng(3)
        psi = rng.uniform(-1.0, 1.0, size=(n_domain, 5))
        np.testing.assert_array_equal(op.apply(psi), psi)

    def test_periodic_without_a_face_is_refused(self):
        r"""A quadrature-only space cannot name the partner, so it is refused.

        Periodic is the ONLY law whose domain is a different face's
        :math:`\Gamma_+`, which makes the face strictly more load-bearing here
        than for any other law: vacuum needs it to size :math:`\Gamma_-`,
        reflective to certify the pairing, and periodic to know *which face it
        is reading at all*. Guessing would silently realize the pre-B3.4c
        defect — a face returning its own outflow as its inflow.
        """
        with pytest.raises(BoundaryError, match="without a face"):
            SNBoundaryRealizer().realize(
                PeriodicBoundary(axis="x"),
                SNMethodSpace.minimal(Quadrature.lebedev(17)),
            )

    def test_periodic_refuses_a_face_off_its_declared_axis(self):
        r"""A wrap along ``y`` installed on an ``x`` face identifies nothing.

        The ERR-041 mis-declaration class, in the geometry tier: the
        translation :math:`x \mapsto x + L_y` does not carry an x-face
        anywhere, so there is no partner to name. Answering with the
        installation face would realize it as a bare identity on the wrong
        half-trace — silently, since the shapes agree.
        """
        quad = Quadrature.level_symmetric(sn_order=6)
        with pytest.raises(BoundaryError, match="installed on face 'xmin'"):
            SNBoundaryRealizer().realize(
                PeriodicBoundary(axis="y"),
                face_method_space(quad, face="xmin"),
            )


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
