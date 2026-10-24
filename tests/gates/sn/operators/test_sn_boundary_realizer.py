r"""Tests of the SN realizer table ``SN_REALIZERS`` and the typed face operators it returns (wave 2, slice S2b re-posing).

This file tested ``SNBoundaryRealizer`` (Wave 5 / C5.3), an ``isinstance`` ladder over the string-tagged boundary laws. The ladder
is gone: SN admits and realizes through ONE table keyed by unit (a chain of the closed sum, ``ZeroNetwork``, ``InflowTerms``),
whose images are the typed operators of :mod:`orpheus.sn.boundary.realized_operators`. The verification claims of the old file
that SURVIVE are posed here on the new objects, each through the door a problem uses
(``SN_REALIZERS.realize(declared, site, method_space)``), with the independent references the old file built (a numpy gather for
the mirror, a closed form for the Lambertian, the weighted-adjoint law of ERR-076) kept as they were:

* the zero map :math:`\Gamma_+ \to \Gamma_-` (vacuum, ``specular(0)``, ``white(0)``);
* the scaled ordinate permutation (the mirror: hand-computed gather, bit-exact, unit and partial);
* the factored Lambertian (closed-form reference);
* the translation wrap (an identity relabelling between two distinct faces);
* the metric correctness of every realized face (ERR-076, the weighted Hilbert adjoint);
* the zero map realizing on the trace ladder (ERR-041's invariant is witnessed on the ladder);
* the face PACKING ORDER is book-keeping (G6.5).

**Retired with the ladder, each by its structural reason** (the demotion check of ``retirement-audit``: the claim is unspellable now,
or a gate of the new tree owns it):

* *the white orientation guard* (``TestWhiteOrientationGuard``: a law declared for ``xmax`` installed on ``xmin``): a shared wall
  holds no axis and no outward sign, and the hemisphere of the Lambertian IS ``trace.outflow_space(face)``, derived from the face
  name once, so the mismatch is unspellable. The replacement witness is law L-E (naturality of the mirror,
  ``tests/gates/sn/architecture/test_sn_realizer_table.py::test_mirror_naturality``), which catches a wrong-hemisphere realizer
  without reading ``face_normal``; the positive control (each of four 2-D faces realizes the Lambertian over its own
  :math:`\Gamma_+`) is kept below.
* *the closure-free albedo refusal* and *the periodic wrap off its declared axis*: a shape-less albedo and a wrap with an axis
  letter are unspellable (``WallCondition`` is shaped by its chain; a deck's partner is derived from the face by ``pairing_of``).
* *``TestRealizerIdentity`` and ``TestRealizeUnknownLawRaises``*: there is no realizer class to name and no law to be
  unrecognised; an unregistered unit meets the table's door (``test_the_unregistered_door_names_the_state_cone_...``).
"""

from __future__ import annotations

from fractions import Fraction

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.geometry.faces import face_normal
from orpheus.geometry.transformation import Permutation
from orpheus.numerics.elementary_operators import PermutationOperator
from orpheus.numerics.operator import IdentityOperator, OperatorSum, TensorProductOperator, ZeroMorphism
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.boundary.realized_operators import RankOneReemission, ScaledOrdinatePermutation, ZeroMap, _crossing
from orpheus.sn.boundary.realizer import SN_REALIZERS
from orpheus.sn.mesh.method_space import SNMethodSpace
from orpheus.transport.boundary_conditions import BoundaryError, BoundaryRefused, FaceSite, RefusalAxis, WallCondition
from tests._harness.float_bounds import gamma
from tests._harness.references import mirror_partner_indices
from tests.gates.sn._test_helpers import face_method_space, face_trace


def _realize(declared, space, face: str | None = None):
    """The face operator of ``declared`` on ``face`` (default: the method space's own), through the door a problem uses."""
    face = face if face is not None else space.face
    site = FaceSite(face, CoordSystem.CARTESIAN, 1 + face_normal(face)[0])
    return SN_REALIZERS.realize(declared, site, space)


def _inflow(space):
    """The face's :math:`\\Gamma_-` ordinate indices, read off the trace ladder."""
    return space.inflow_indices_for_face(space.face)


def _outflow(space):
    """The face's :math:`\\Gamma_+` ordinate indices, read off the trace ladder."""
    return space.trace.outflow_indices_for_face(space.face)


def _scaled_parts(op):
    """The ``(alpha, word)`` of a scaled operator: the one-term ``OperatorSum``."""
    assert isinstance(op, OperatorSum), f"expected alpha * word, got {type(op).__name__}"
    [(alpha, word)] = op.terms
    return alpha, word


def assert_lambertian_kernel(kernel, space, quad) -> None:
    r"""The surviving claims about a realized Lambertian, posed BEHAVIOURALLY.

    ⭐ **The load-bearing leg is the SIZE one**: the kernel ACCEPTS :math:`|\Gamma_+|` rows and REFUSES ``quad.n_ordinates``. It is
    strictly stronger than reading an attribute, because it holds for any future spelling and cannot be satisfied by an operator
    that merely reports the right number; without it a revert to the full-face kernel would keep every structural claim green.
    """
    n_out = int(_outflow(space).size)
    n_in = int(_inflow(space).size)
    assert n_out < quad.n_ordinates, "fixture no longer exercises the narrowing: Gamma_+ is the whole face"
    out = np.asarray(kernel.apply(np.ones((n_out,))))
    assert out.shape == (n_in,), f"the realized Lambertian emits {out.shape} on Gamma_- of size {n_in}"
    with pytest.raises(ValueError):
        kernel.apply(np.ones((quad.n_ordinates,)))


def lambertian_reference(quad, axis, sign, psi_out, space):
    r"""An INDEPENDENT closed form for the Lambertian, in explicit numpy.

    .. math::

        (R\psi)(\Omega) = \frac{\sum_{\Gamma_+} w\,|\Omega'\cdot\hat n|\,\psi(\Omega')}{\sum_{\Gamma_+} w\,|\Omega'\cdot\hat n|},
        \qquad \Omega \in \Gamma_-.

    Rebuilt from the quadrature's own ``weights`` and direction cosines with no production operator in the path, so a rewire of the
    realizer cannot compare its output with itself through a wrapper.
    """
    mu_n = np.asarray(getattr(quad, f"mu_{axis}"), dtype=float) * sign
    out_i = np.asarray(_outflow(space))
    cos_w = np.asarray(quad.weights, dtype=float)[out_i] * mu_n[out_i]
    psi_out = np.asarray(psi_out)
    contracted = (cos_w.reshape((-1,) + (1,) * (psi_out.ndim - 1)) * psi_out).sum(axis=0)
    n_in = int(np.asarray(_inflow(space)).size)
    return np.broadcast_to((contracted / cos_w.sum())[None, ...], (n_in,) + contracted.shape).copy()


# ─────────────────────────────────────────────────────────────────────
# 1. Vacuum
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeVacuum:
    r"""The empty network realizes to the ZERO MAP :math:`\Gamma_+ \to \Gamma_-`.

    Vacuum's entire content is :math:`R = 0`, so the honest object is the zero map between the two half-traces. The assertions
    state: the image is zero, it lives on :math:`\Gamma_-`, and the operator is NOT an endomorphism of the face slot (with the
    domain narrowed to :math:`\Gamma_+` the outflow rows are simply not in the operator's domain).
    """

    def test_lebedev17_xmin_is_the_zero_map_onto_gamma_minus(self):
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmin")  # xmin: outward normal -x, inflow is mu_x > 0
        rng = np.random.default_rng(42)
        psi = rng.uniform(0.5, 2.0, size=(quad.n_ordinates, 4, 2))
        op = _realize(WallCondition.vacuum(), space)
        assert isinstance(op, ZeroMap) and op.is_vacuum
        out = op.apply(psi[_outflow(space)])
        assert out.shape == (_inflow(space).size, 4, 2)
        np.testing.assert_array_equal(out, 0.0)
        # |Gamma_+| == |Gamma_-| on every quadrature in the tree, so the shape check above cannot tell Gamma_+ -> Gamma_- from
        # Gamma_+ -> Gamma_+ (vv Mode 12): a full-face input must not produce a full-face image.
        assert op.apply(psi).shape[0] != quad.n_ordinates

    def test_lebedev17_ymax_is_the_zero_map_onto_gamma_minus(self):
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="ymax", faces=("xmin", "xmax", "ymin", "ymax"))
        # An independent spelling of the face's orientation: Omega.n_ymax = +mu_y, so inflow <=> mu_y < 0.
        np.testing.assert_array_equal(_inflow(space), np.flatnonzero(quad.mu_y < -1e-12))
        psi = np.arange(quad.n_ordinates * 3, dtype=float).reshape(quad.n_ordinates, 3)
        out = _realize(WallCondition.vacuum(), space).apply(psi[_outflow(space)])
        assert out.shape == (_inflow(space).size, 3)
        np.testing.assert_array_equal(out, 0.0)

    def test_vacuum_without_the_trace_ladder_is_refused(self):
        """A quadrature-only method space carries no trace ladder, so the zero map cannot be bound to its two spaces."""
        quad = Quadrature.lebedev(17)
        with pytest.raises(BoundaryError, match="trace ladder"):
            _realize(WallCondition.vacuum(), SNMethodSpace.minimal(quad), face="xmin")

    def test_the_zero_map_is_a_bound_zero_morphism_with_both_space_hooks(self):
        r"""The derived operator is a bound :class:`ZeroMorphism`: forward emits the zero of :math:`\Gamma_-`, transpose the zero of
        :math:`\Gamma_+`, both derived from the DECLARED pair. Relying on the endomorphic ``0.0 * x`` echo would be wrong in
        principle and merely lucky in practice (the shapes coincide on every reachable fixture)."""
        quad = Quadrature.gauss_legendre(8)
        space = face_method_space(quad, face="xmin")
        op = _realize(WallCondition.vacuum(), space)
        assert isinstance(op.inner, ZeroMorphism)
        np.testing.assert_array_equal(op.apply(np.ones((_outflow(space).size, 3))), np.zeros((_inflow(space).size, 3)))
        np.testing.assert_array_equal(op.apply_transpose(np.ones((_inflow(space).size, 3))), np.zeros((_outflow(space).size, 3)))

    @pytest.mark.parametrize("spelling", [WallCondition.specular(0.0), WallCondition.white(0.0)], ids=["specular(0)", "white(0)"])
    def test_a_wall_of_zero_weight_is_the_same_zero_map(self, spelling):
        """``specular(0.0)`` and ``white(0.0)`` ARE ``vacuum()`` (one value), so they realize the ``ZeroMap``, not a zero-scaled kernel."""
        space = face_method_space(Quadrature.gauss_legendre(8), face="xmax")
        assert _realize(spelling, space) == _realize(WallCondition.vacuum(), space)


# ─────────────────────────────────────────────────────────────────────
# 2. Specular (the mirror)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeSpecular:
    r"""The specular chain (and the mirror deck) realize to ``alpha * PermutationOperator(perm)``."""

    def test_specular_unit_weight_lebedev_matches_hand_computed(self):
        r"""At weight 1 the realized op MUST bit-match the hand-computed narrowed gather ``psi[mirror_partner][inflow]``.

        The expected value is built from :func:`~tests._harness.references.mirror_partner_indices` (pure numpy over the direction
        cosines) and the inflow index set alone, so neither a remap error nor a pairing-derivation error in production can cancel
        against the reference (``algebra-of-record`` structural independence).
        """
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmin")
        rng = np.random.default_rng(7)
        psi = rng.uniform(-1.0, 1.0, size=(quad.n_ordinates, 5, 3))
        expected = psi[mirror_partner_indices(quad, "x")][_inflow(space)]
        for declared in (WallCondition.specular(1.0), MirrorDeck()):
            np.testing.assert_array_equal(_realize(declared, space).apply(psi[_outflow(space)]), expected)

    def test_specular_unit_weight_returns_tensor_product(self):
        r"""At weight 1 the derived operator is a :class:`TensorProductOperator` whose factors are ``(PermutationOperator,
        IdentityOperator)``, with the permutation on the **reduced** ordinate axis.

        The permutation's length is :math:`|\Gamma_+|`, not ``quad.n_ordinates``: that single assertion is what makes the narrowing
        TYPE-VISIBLE on this object, since a revert to the full-``N`` permutation would keep every structural claim green.
        """
        quad = Quadrature.level_symmetric(4)
        space = face_method_space(quad, face="xmin")
        op = _realize(WallCondition.specular(1.0), space)
        assert isinstance(op, ScaledOrdinatePermutation) and op.permutes_ordinates and op.return_fraction == 1.0
        assert isinstance(op.inner, TensorProductOperator) and len(op.inner.ops) == 2
        assert isinstance(op.inner.ops[0], PermutationOperator) and isinstance(op.inner.ops[1], IdentityOperator)
        assert op.inner.ops[0].perm.size == _outflow(space).size < quad.n_ordinates

    def test_specular_partial_weight_matches_hand_computed(self):
        """At weight 0.7 the realized op output equals ``0.7 * psi[mirror_partner][inflow]`` bit-exactly."""
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmin")
        rng = np.random.default_rng(9)
        psi = rng.uniform(-1.0, 1.0, size=(quad.n_ordinates, 5))
        expected = 0.7 * psi[mirror_partner_indices(quad, "x")][_inflow(space)]
        np.testing.assert_array_equal(_realize(WallCondition.specular(0.7), space).apply(psi[_outflow(space)]), expected)

    def test_specular_partial_weight_returns_scaled_tensor_product(self):
        """At weight != 1 the specular wall's operator is ``alpha * TP`` (the one-term combination), TP the 2-factor product."""
        quad = Quadrature.level_symmetric(4)
        space = face_method_space(quad, face="xmin")
        alpha, tp = _scaled_parts(_realize(WallCondition.specular(0.5), space).inner)
        assert alpha == 0.5
        assert isinstance(tp, TensorProductOperator) and len(tp.ops) == 2
        assert isinstance(tp.ops[0], PermutationOperator) and isinstance(tp.ops[1], IdentityOperator)
        assert tp.ops[0].perm.size == _outflow(space).size < quad.n_ordinates

    def test_the_scaled_mirror_composed_with_its_transpose_is_alpha_squared(self):
        r"""``Bᵀ(B(x)) == alpha² x`` on :math:`\Gamma_+` (migrated from ``tests/gates/geometry/test_boundary.py``, wave 2 S2b D.14).

        The narrowed ``B`` is ``alpha`` times a bijection :math:`\Gamma_+ \to \Gamma_-`, so ``Bᵀ B = alpha² I`` on :math:`\Gamma_+`:
        a mirror is a bijection and not merely an involution of a larger space. (A :math:`\Gamma_+ \to \Gamma_-` map cannot be composed
        with itself; the old spelling ``B(B(x))`` worked only while the law was a square endomorphism of the whole face.)
        """
        quad = Quadrature.gauss_legendre(n_ordinates=8)
        space = face_method_space(quad, face="xmax")
        x = np.random.default_rng(7).standard_normal((int(_outflow(space).size), 2))
        op = _realize(WallCondition.specular(0.7), space)
        np.testing.assert_allclose(op.apply_transpose(op.apply(x)), 0.7 ** 2 * x, rtol=1e-13, atol=1e-14)

    @pytest.mark.catches("ERR-074")  # migrated from test_deck_kernel.py (the retired ``_deck_kernel`` refusal row), wave 2 S2b D.14
    def test_an_uncertified_mirror_is_refused_on_the_angular_resolution_axis(self):
        r"""A rule NOT closed under the mirror is refused at ``realize()`` with the pairing's own diagnosis.

        ``product(4, 5)`` has no x-mirror closure (odd :math:`n_\varphi`, the ERR-074 fixture): ``ordinate_pairing`` raises the
        not-closed error, which the table's ``admit`` turns into ``Refused(ANGULAR_RESOLUTION)`` and the one door into
        ``BoundaryRefused`` BEFORE any realizer body runs. Pinned on the distinctive message fragment: ``BoundaryRefused`` is a
        ``ValueError`` and a type-only pin could not tell this from a lookup miss.
        """
        quad = Quadrature.product(n_mu=4, n_phi=5)
        for declared in (WallCondition.specular(1.0), MirrorDeck()):
            with pytest.raises(BoundaryRefused, match="no ordinate pairing about axis") as info:
                _realize(declared, face_method_space(quad, face="xmin"))
            assert info.value.admission.axis is RefusalAxis.ANGULAR_RESOLUTION

    def test_a_pairing_needs_a_bijection_between_the_half_traces(self):
        r"""A face where :math:`|\Gamma_-| \neq |\Gamma_+|` has NO realization, and says so.

        A pairing is a bijection between the two half-traces, so an asymmetric face is an ill-posed request. No real quadrature x
        face pair in the tree produces one (every one measured :math:`|I| = |O|`), so the guard is unreachable from a mesh; it is
        exercised on the crossing check itself with two half-traces of different rules.
        """
        small, large = Quadrature.gauss_legendre(4), Quadrature.gauss_legendre(8)
        outflow, inflow = face_trace(small).outflow_space("xmax"), face_trace(large).inflow_space("xmax")
        with pytest.raises(ValueError, match=r"BIJECTION"):
            _crossing(Permutation.identity(8), outflow, inflow, "probe")

    def test_a_spaceless_method_space_cannot_pose_the_pairing(self):
        r"""The pairing is returned fully bound :math:`\Gamma_+(f') \to \Gamma_-(f)`: a method space that stands up no half-trace
        spaces is refused (before the wave-2 carve this path realized a silently UNBOUND permutation whose ``.H`` was the bare
        Euclidean transpose, the ERR-076 shape)."""
        quad = Quadrature.gauss_legendre(8)
        spaceless = SNMethodSpace(quadrature=quad, face="xmax")
        with pytest.raises(BoundaryError, match="trace ladder"):
            _realize(WallCondition.specular(1.0), spaceless)


# ─────────────────────────────────────────────────────────────────────
# 3. White (Lambertian)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeWhite:
    r"""The diffuse chain realizes to ``alpha * (emission @ contraction)``, :math:`\Gamma_+ \to \Gamma_-`.

    The kernel's body is the factored Lambertian of :mod:`orpheus.sn.boundary.angular`, verified independently by the
    current-conservation and reciprocity gates; these legs pin the realizer's dispatch chain against an INDEPENDENT closed form.
    """

    def test_white_unit_weight_lebedev_matches_the_closed_form(self):
        r"""At weight 1 the realized op bit-matches :func:`lambertian_reference` (the closed form, divisor the inflow mass: bit-equal
        to the outflow sum on every registry rule, the census of ``test_realized_boundary.py``)."""
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        rng = np.random.default_rng(123)
        psi_out = rng.uniform(0.0, 2.0, size=(quad.n_ordinates, 5, 3))[_outflow(space)]
        got = _realize(WallCondition.white(1.0), space).apply(psi_out)
        assert got.shape == (_inflow(space).size, 5, 3)
        np.testing.assert_array_equal(got, lambertian_reference(quad, "x", +1, psi_out, space))

    def test_white_unit_weight_returns_tensor_product(self):
        """At weight 1 the derived operator is a 2-factor product ``(the Lambertian chain, IdentityOperator)``."""
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        op = _realize(WallCondition.white(1.0), space)
        assert isinstance(op, RankOneReemission) and not op.permutes_ordinates and op.return_fraction == 1.0
        assert isinstance(op.inner, TensorProductOperator) and len(op.inner.ops) == 2
        assert isinstance(op.inner.ops[1], IdentityOperator)
        assert_lambertian_kernel(op.inner.ops[0], space, quad)

    def test_white_partial_weight_returns_scaled_tensor_product(self):
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        alpha, tp = _scaled_parts(_realize(WallCondition.white(0.3), space).inner)
        assert alpha == 0.3
        assert isinstance(tp, TensorProductOperator) and len(tp.ops) == 2 and isinstance(tp.ops[1], IdentityOperator)
        assert_lambertian_kernel(tp.ops[0], space, quad)

    def test_white_partial_weight_matches_weight_times_the_closed_form(self):
        r"""At weight 0.3 the output equals ``0.3 * lambertian_reference`` within ``nulp = 4``: ``(0.3 * TP).apply`` scales the
        inner output from the outside, the hand-computed product may differ by one ULP under IEEE-754 non-associativity."""
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmax")
        psi_out = np.random.default_rng(11).uniform(0.0, 2.0, size=(quad.n_ordinates, 4))[_outflow(space)]
        np.testing.assert_array_almost_equal_nulp(
            _realize(WallCondition.white(0.3), space).apply(psi_out), 0.3 * lambertian_reference(quad, "x", +1, psi_out, space), nulp=4,
        )

    def test_a_face_the_quadrature_does_not_carry_has_no_trace(self):
        r"""A ``z`` face on a 1-D Gauss-Legendre rule is a RANK MISMATCH between the face and the cubature, and says so where the
        trace is built (the single face-name to signed-projection primitive, ``build_omega_dot_n``), before any realizer runs."""
        with pytest.raises(ValueError, match="requires genuine mu_z"):
            face_method_space(Quadrature.gauss_legendre(8), face="zmax", faces=("zmin", "zmax"))

    def test_white_without_the_trace_ladder_is_refused(self):
        r"""White's DOMAIN is :math:`\Gamma_+`, so a quadrature-only method space cannot realize it (every leg migrated off
        ``SNMethodSpace.minimal``)."""
        with pytest.raises(BoundaryError, match="trace ladder"):
            _realize(WallCondition.white(1.0), SNMethodSpace.minimal(Quadrature.lebedev(17)), face="xmax")

    def test_every_face_of_a_box_realizes_over_its_own_outflow(self):
        r"""Positive control for the retired orientation guard: each of the four 2-D faces realizes the Lambertian, and the kernel's
        domain is THAT face's :math:`\Gamma_+` (the hemisphere is derived from the face name, so no mismatch is spellable)."""
        quad = Quadrature.level_symmetric(4)
        faces = ("xmin", "xmax", "ymin", "ymax")
        for face in faces:
            space = face_method_space(quad, face=face, faces=faces)
            op = _realize(WallCondition.white(1.0), space)
            assert isinstance(op, RankOneReemission), face
            assert op.domain == space.trace.outflow_space(face) and op.codomain == space.trace.inflow_space(face), face
            assert_lambertian_kernel(op.inner.ops[0], space, quad)


# ─────────────────────────────────────────────────────────────────────
# 4. The weighted walls and the mixed wall
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeTheWeightedWalls:
    r"""Each chain of a wall realizes at its own weight, and a wall of two chains realizes to the sum of the images.

    Replaces ``TestRealizeAlbedo``: a shape-less albedo is unspellable, and ``specular(a)`` and ``white(a)`` are the two shapes. The
    zero weight is the ``ZeroMap`` (``TestRealizeVacuum``); the unit and partial weights are pinned per chain above and below.
    """

    @staticmethod
    def _space(n_ord: int = 8):
        return face_method_space(Quadrature.gauss_legendre(n_ord), face="xmax")

    def test_weight_one_is_the_bare_kernel_product_and_a_partial_weight_the_scaled_one(self):
        space = self._space()
        unit = _realize(WallCondition.specular(1.0), space)
        assert isinstance(unit.inner, TensorProductOperator) and isinstance(unit.inner.ops[0], PermutationOperator)
        half = _realize(WallCondition.specular(0.5), space)
        alpha, tp = _scaled_parts(half.inner)
        assert alpha == 0.5 and isinstance(tp, TensorProductOperator)

    def test_a_mixed_wall_is_the_sum_of_the_images_in_the_canonical_order(self):
        space = self._space()
        op = _realize(WallCondition.from_amplitudes(0.3, 0.5), space)
        assert op == _realize(WallCondition.specular(0.3), space) + _realize(WallCondition.white(0.5), space)
        assert op.return_fraction == 0.8 and not op.is_vacuum

    def test_the_diffuse_chain_reaches_the_lambertian_kernel(self):
        """The OTHER chain reaches the OTHER body: a dispatch that ignored the diffuse chain would pass the specular rows alone."""
        _, tp = _scaled_parts(_realize(WallCondition.white(0.5), self._space()).inner)
        assert isinstance(tp, TensorProductOperator) and not isinstance(tp.ops[0], PermutationOperator)


# ─────────────────────────────────────────────────────────────────────
# 5. Periodic (the translation wrap)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestRealizeTranslation:
    r"""The translation deck realizes to the identity RELABELLING between the partner's :math:`\Gamma_+` and this face's
    :math:`\Gamma_-`: a typed arrow between two distinct spaces, derived from the wrap's motion."""

    def test_the_wrap_returns_tensor_product_and_passes_through(self):
        r"""The derived operator is a 2-factor product ``(PermutationOperator, IdentityOperator)``, the angular factor a BOUND identity
        relabelling whose apply passes values through unchanged. Compare values, NOT identity: the product fold returns a copy.

        The probe is :math:`\Gamma_+`-shaped rather than full-face (49 of 74 ordinates on ``lebedev(17)``): that the two shapes
        differ is what makes this a statement about the narrowed law rather than about an endomorphism that happens to be the
        identity.
        """
        quad = Quadrature.lebedev(17)
        space = face_method_space(quad, face="xmin")
        op = _realize(TranslationDeck(), space)
        assert isinstance(op, ScaledOrdinatePermutation) and not op.permutes_ordinates and op.domain_face == "xmax"
        angular, group = op.inner.ops
        assert isinstance(op.inner, TensorProductOperator) and isinstance(angular, PermutationOperator) and isinstance(group, IdentityOperator)
        trace = space.trace
        assert op.domain == trace.outflow_space("xmax") and op.codomain == trace.inflow_space("xmin")
        n_domain = np.asarray(trace.outflow_indices_for_face("xmax")).size
        assert n_domain != quad.n_ordinates
        psi = np.random.default_rng(3).uniform(-1.0, 1.0, size=(n_domain, 5))
        np.testing.assert_array_equal(op.apply(psi), psi)
        np.testing.assert_array_equal(op.apply_transpose(psi), psi)

    def test_a_translation_without_the_trace_ladder_is_refused(self):
        r"""A quadrature-only space cannot name the partner's half-trace, so it is refused: periodic is the ONLY realization whose
        domain is a different face's :math:`\Gamma_+`, which makes the trace ladder strictly more load-bearing here."""
        with pytest.raises(BoundaryError, match="trace ladder"):
            _realize(TranslationDeck(), SNMethodSpace.minimal(Quadrature.lebedev(17)), face="xmin")


# ─────────────────────────────────────────────────────────────────────
# 6. The realized face is METRIC-correct (ERR-076)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
@pytest.mark.catches("ERR-076")
class TestTheRealizedFaceIsMETRICCorrect:
    r"""⭐ The consequence of binding every operator to its half-trace spaces, gated at the tier the physics lives in.

    :class:`~orpheus.numerics.operator.AdjointOperator` reads the operator's spaces to apply the metrics, so an UNBOUND operator's
    ``.H`` is not a weaker Hilbert adjoint: it is the **Euclidean transpose**, silently. `[M]` on a 3-group ``gauss_legendre(8)``
    xmin face the two agree to **0.0** for the specular mirror (whose metric cancels exactly, :math:`G_{\Gamma_-} = G_{\Gamma_+}\circ\pi`)
    and differ by **87 % relative** for the Lambertian. So a reflective-only fixture cannot see this, which is why the diffuse wall is
    parametrized here and why the twin below asserts the fixture can still tell the two apart.
    """

    _DECLARED = {
        "specular": lambda: WallCondition.specular(1.0),
        "lambertian": lambda: WallCondition.white(1.0),
        "specular_scaled": lambda: WallCondition.specular(0.3),
        "lambertian_scaled": lambda: WallCondition.white(0.4),
        "translation": lambda: TranslationDeck(),
        "mirror_deck": lambda: MirrorDeck(),
    }

    @staticmethod
    def _realized(key, n_ordinates=8):
        quad = Quadrature.gauss_legendre(n_ordinates=n_ordinates)
        space = face_method_space(quad, face="xmax")
        return _realize(TestTheRealizedFaceIsMETRICCorrect._DECLARED[key](), space), quad

    @pytest.mark.parametrize("key", sorted(_DECLARED))
    def test_the_realized_operator_carries_its_spaces(self, key):
        """The precondition: without both ends, ``.H`` cannot be weighted."""
        op, _ = self._realized(key)
        assert op.domain is not None, f"{key} lost its domain"
        assert op.codomain is not None, f"{key} lost its codomain"

    @pytest.mark.parametrize("key", sorted(_DECLARED))
    def test_the_realized_operator_satisfies_the_WEIGHTED_adjoint_law(self, key):
        r""":math:`\langle Rx, y\rangle_{\Gamma_-} = \langle x, R^{*}y\rangle_{\Gamma_+}` under the partial-current metric.

        ⚠ **What this row gates is that the metric is REACHED**, not the adjoint formula: given :math:`R^{*} = G_{\Gamma_+}^{-1}
        R^{\mathsf T} G_{\Gamma_-}` the identity is algebra for any :math:`R`. It discriminates because an operator with no spaces never
        gets that far: :class:`AdjointOperator` applies no metric at all and :math:`R^{*}` collapses to the bare :math:`R^{\mathsf T}`,
        which does NOT satisfy the weighted law.
        """
        op, _ = self._realized(key)
        rng = np.random.default_rng(20260806)
        x = rng.standard_normal(op.domain.shape)
        y = rng.standard_normal(op.codomain.shape)
        assert op.codomain.inner_product(op.apply(x), y) == pytest.approx(op.domain.inner_product(x, op.H.apply(y)), rel=1e-13)

    def test_the_TRANSLATION_metric_agrees_across_the_pair(self):
        r""":math:`G_{\Gamma_+(x_{\min})} = G_{\Gamma_-(x_{\max})}`: the opposite-normals theorem at the METRIC tier, across two DIFFERENT
        faces. The agreement is bit-exact (:math:`|\Omega\cdot\hat n|` is the modulus of a sign flip and :math:`w_n` is one array); the
        operator-level coincidence of the weighted ``.H`` and the bare transpose is exact only analytically (1 ULP measured), so it
        is stated at ``nulp = 4``: the FP-non-associativity class, not a relaxed method tolerance (the Lambertian control sits 15
        orders away)."""
        op, _ = self._realized("translation")
        np.testing.assert_array_equal(op.domain.inner_product_weights, op.codomain.inner_product_weights)
        y = np.random.default_rng(20260807).standard_normal(op.codomain.shape)
        np.testing.assert_array_almost_equal_nulp(np.asarray(op.H.apply(y)), np.asarray(op.apply_transpose(y)), nulp=4)

    def test_the_LAMBERTIAN_fixture_can_tell_a_weighted_adjoint_from_a_bare_one(self):
        """The discrimination check: a metric that cancels makes the weighted law hold for the UNWEIGHTED operator too (the mirror).
        This asserts the Lambertian row is not in that regime, so its green is evidence."""
        op, _ = self._realized("lambertian")
        y = np.random.default_rng(20260806).standard_normal(op.codomain.shape)
        hilbert, euclidean = np.asarray(op.H.apply(y)), np.asarray(op.apply_transpose(y))
        relative = np.abs(hilbert - euclidean).max() / np.abs(hilbert).max()
        assert relative > 0.1, f"the Hilbert adjoint and the Euclidean transpose differ by only {relative:.2e} relative on this fixture"

    def test_the_SPECULAR_mirror_is_the_documented_blind_case(self):
        r"""And the negative leg: on the mirror the weighted adjoint and the bare transpose coincide (:math:`G_{\Gamma_-}\circ\pi =
        G_{\Gamma_+}` because a mirror preserves :math:`|\Omega\cdot\hat n|\,w_n`). Two legs, each at the tier where it is exact.

        **(1) The metric tier, bit-exact, a float theorem.** With :math:`R^{\mathsf T}` a permutation (asserted), metric cancellation is
        the identity :math:`R^{\mathsf T}(G_{\Gamma_-} y) = G_{\Gamma_+}(R^{\mathsf T} y)`: each entry is the product of the SAME two floats.

        **(2) The operator tier, within a derived bound.** :math:`R^{*}y` evaluates entry :math:`i` as :math:`\mathrm{fl}(\mathrm{fl}(g\,y)/g)`,
        two roundings around an exact permutation, so :math:`|(R^{*}y)_i - (R^{\mathsf T}y)_i| \le \gamma_2\,|(R^{\mathsf T}y)_i|` (Higham
        Lemma 3.1); ``array_equal`` would be :math:`(a y)/a = y`, not an IEEE identity.
        """
        op, _ = self._realized("specular")
        y = np.random.default_rng(20260806).standard_normal(op.codomain.shape)
        transposed = np.asarray(op.apply_transpose(y))
        if not np.array_equal(np.sort(transposed.ravel()), np.sort(y.ravel())):
            pytest.fail("the specular transpose is not a permutation of its input")
        np.testing.assert_array_equal(
            np.asarray(op.apply_transpose(op.codomain.apply_metric(y))), np.asarray(op.domain.apply_metric(transposed)),
            err_msg="G_{Gamma-} o pi != G_{Gamma+}: the mirror does not preserve |Omega.n|*w",
        )
        hilbert = np.asarray(op.H.apply(y)).ravel()
        bound = gamma(2)
        for h, t in zip(hilbert, transposed.ravel()):
            if abs(Fraction(float(h)) - Fraction(float(t))) > bound * abs(Fraction(float(t))):
                pytest.fail(f"R* y = {h!r} vs R^T y = {t!r}: outside gamma_2, the metric does not cancel through the specular adjoint")


# ─────────────────────────────────────────────────────────────────────
# 7. The zero map reads its half-space off the trace ladder (ERR-041's invariant is witnessed there)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestVacuumRealizesOnTheTrace:
    r"""The zero map of a vacuum face is built from the trace ladder alone.

    ERR-041 (the zero map fills only the incoming trace) used to be a realizer guard that cross-checked a hand-supplied inflow index
    set against the face name. The method space holds no index set any more, so the guard and its tests retired (WF2.4); the invariant
    itself, that ``trace.inflow_space(face)`` holds exactly the ordinates with :math:`\Omega\cdot n_f < 0`, is witnessed on the ladder in
    ``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates`` (the ``catches`` marker lives
    there). What stays here is what the realizer still owes: a face-ful method space realizes, a faceless or unknown-face one refuses.
    """

    def test_correct_orientation_realizes(self):
        """Positive control: each face realizes through the realizer, and the half-space it reads is the written-out per-face predicate."""
        quad = Quadrature.gauss_legendre(8)
        for face, inflow in (("xmax", np.flatnonzero(quad.mu_x < -1e-12)), ("xmin", np.flatnonzero(quad.mu_x > +1e-12))):
            space = face_method_space(quad, face=face)
            np.testing.assert_array_equal(_inflow(space), inflow)
            assert isinstance(_realize(WallCondition.vacuum(), space), ZeroMap)

    def test_a_faceless_space_carries_no_trace(self):
        r"""A method space with no trace ladder has no :math:`\Gamma_-` to bind, so the realize REFUSES loudly instead of building an
        unbindable zero map."""
        faceless = SNMethodSpace(quadrature=Quadrature.gauss_legendre(8))
        with pytest.raises(BoundaryError, match="trace ladder"):
            _realize(WallCondition.vacuum(), faceless, face="xmax")

    def test_unknown_face_name_fails_loud_on_the_ladder(self):
        """A face name outside the layout is refused loudly, by the trace ladder every half-space is read from (a ``ValueError`` naming the
        available faces). The realizer itself no longer validates the name: the guard that did (it read the signed projection through
        ``build_omega_dot_n``) retired with the index sets, and the name now reaches the ladder at the first read."""
        space = face_method_space(Quadrature.gauss_legendre(8), face="xmax")
        assert _inflow(space).size > 0  # positive control: a known face reads
        with pytest.raises(ValueError, match="Unknown face"):
            space.inflow_indices_for_face("right")


# ─────────────────────────────────────────────────────────────────────
# 8. G6.5 acceptance: the face PACKING ORDER is book-keeping
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestFacePackingOrderIsBookkeeping:
    r"""⭐ G6.5's acceptance gate: permute the face packing order, bit-identity.

    The face packing order moves only which flat offset each face's slot occupies, book-keeping with no arithmetic content. So the
    same declaration realized on two traces differing ONLY in layout order must produce per-face outputs identical to the BIT.

    ⭐ **The 3-D arm carries the y-faces, and they are load-bearing**: the ONLY flat-offset read in the whole realization path is the
    per-face metric slice, and with both faces on ONE axis the two slots' metric slices are bit-identical, so a wrong-slot metric
    read is INVISIBLE on a same-axis fixture (Mode 12 at the fixture tier). ``gauss_legendre`` structurally cannot carry a y-face, so
    its arm stays 2-face and gates the permutation class only. Per ``vv-principles`` #23 the valid positive control for an
    INVARIANCE gate is knob-dependent (an ignored-order mutation reddening the activation rows).
    """

    _QUADS = {
        "gauss_legendre(8)": (lambda: Quadrature.gauss_legendre(n_ordinates=8), ("xmin", "xmax"), ("xmax", "xmin")),
        "product(4,4)": (lambda: Quadrature.product(n_mu=4, n_phi=4), ("xmin", "xmax", "ymin", "ymax"), ("ymin", "xmax", "xmin", "ymax")),
    }
    _DECLARED = {
        "specular": lambda: WallCondition.specular(1.0),
        "lambertian": lambda: WallCondition.white(1.0),
        "translation": lambda: TranslationDeck(),
        "vacuum": lambda: WallCondition.vacuum(),
    }

    @pytest.mark.parametrize("quad_name", sorted(_QUADS))
    def test_the_packing_order_genuinely_moves_the_buffer(self, quad_name):
        """Activation: without this, the class could green on an ignored order. The two orderings place at least one face at a
        different flat offset, while the per-face spaces compare equal and carry bit-identical metrics (the FACE does not move)."""
        factory, forward, reversed_ = self._QUADS[quad_name]
        quad = factory()
        trace_a, trace_b = face_trace(quad, forward), face_trace(quad, reversed_)
        offsets_a = {f: trace_a.layout.faces[f].offset for f in forward}
        offsets_b = {f: trace_b.layout.faces[f].offset for f in forward}
        assert offsets_a != offsets_b, "the two orderings packed every face identically: the fixture is vacuous"
        for face in forward:
            assert trace_a.outflow_space(face) == trace_b.outflow_space(face)
            np.testing.assert_array_equal(
                np.asarray(trace_a.outflow_space(face).inner_product_weights), np.asarray(trace_b.outflow_space(face).inner_product_weights),
            )
            np.testing.assert_array_equal(
                np.asarray(trace_a.outflow_space(face).ordinate_indices), np.asarray(trace_b.outflow_space(face).ordinate_indices),
            )

    @pytest.mark.parametrize("quad_name", sorted(_QUADS))
    @pytest.mark.parametrize("key", sorted(_DECLARED))
    def test_realized_face_is_bit_identical_across_packing_order(self, quad_name, key):
        """``apply``, ``apply_transpose`` and the weighted ``.H``, all exact: ``array_equal``, not ``allclose``, because nothing
        arithmetic changed, so any drift at all is a packing dependence leaking into the math."""
        factory, forward, reversed_ = self._QUADS[quad_name]
        quad = factory()
        op_a, op_b = (_realize(self._DECLARED[key](), face_method_space(quad, face="xmax", faces=faces)) for faces in (forward, reversed_))
        rng = np.random.default_rng(20260807)
        x = rng.standard_normal(op_a.domain.shape)
        y = rng.standard_normal(op_a.codomain.shape)
        np.testing.assert_array_equal(op_a.apply(x), op_b.apply(x))
        np.testing.assert_array_equal(np.asarray(op_a.apply_transpose(y)), np.asarray(op_b.apply_transpose(y)))
        np.testing.assert_array_equal(np.asarray(op_a.H.apply(y)), np.asarray(op_b.H.apply(y)))
