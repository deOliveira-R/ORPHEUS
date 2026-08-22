"""Foundation tests for Issue #197 PR-TYPED-3 typed source types.

Pins the structural contract of :class:`ScalarSourceSink` and
:class:`AngularSourceSink` — shape validation, dunder algebra
(within-type AND the load-bearing cross-type), and the bit-identity
of the new dunder path with the procedural
``np.broadcast_to(...).copy() + Q_aniso`` pattern it replaces.

The foundation mark: these tests verify software invariants
(constructor shape check, frozen-ness, dunder algebra) rather than
physics claims.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.source_sinks import (
    ScalarSourceSink,
    AngularSourceSink,
    HarmonicMomentSourceSink,
)
from orpheus.transport.fields.harmonic_moment_flux import HarmonicMomentFlux

from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


# ── Fixtures ─────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Build a small slab :class:`SNMesh` for unit testing."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _stretched_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
    """Same shape as ``_slab_mesh``, doubled width — the cell VOLUMES differ,
    so the carrier mints an UNEQUAL space (the F2 content discriminator)."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, nx + 1),
        mat_ids=np.zeros(nx, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 3, ng: int = 1) -> SNMesh:
    """Build a small 2-D Cartesian :class:`SNMesh`."""
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ════════════════════════════════════════════════════════════════════
# ScalarSourceSink — construction + within-type algebra
# ════════════════════════════════════════════════════════════════════


class TestScalarSource:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        Q = ScalarSourceSink.zeros_on(m)
        assert Q.values.shape == (m.ng, *m.spatial_shape)
        assert np.all(Q.values == 0.0)
        assert isinstance(Q, ScalarSourceSink)
        assert Q.mesh is m

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="ScalarSourceSink.*does not match"):
            ScalarSourceSink.from_mesh(np.zeros((m.ng + 1, *m.spatial_shape)), m)

    def test_within_type_add(self) -> None:
        m = _slab_mesh()
        a = ScalarSourceSink.from_mesh(np.ones((m.ng, *m.spatial_shape)), m)
        b = ScalarSourceSink.from_mesh(2.0 * np.ones((m.ng, *m.spatial_shape)), m)
        c = a + b
        assert isinstance(c, ScalarSourceSink)
        assert np.all(c.values == 3.0)
        # Originals untouched (frozen)
        assert np.all(a.values == 1.0)
        assert np.all(b.values == 2.0)

    def test_within_type_sub(self) -> None:
        m = _slab_mesh()
        a = ScalarSourceSink.from_mesh(3.0 * np.ones((m.ng, *m.spatial_shape)), m)
        b = ScalarSourceSink.from_mesh(np.ones((m.ng, *m.spatial_shape)), m)
        c = a - b
        assert isinstance(c, ScalarSourceSink)
        assert np.all(c.values == 2.0)

    def test_scalar_mul_left_and_right(self) -> None:
        m = _slab_mesh()
        a = ScalarSourceSink.from_mesh(np.ones((m.ng, *m.spatial_shape)), m)
        left = 3.0 * a
        right = a * 3.0
        assert isinstance(left, ScalarSourceSink)
        assert isinstance(right, ScalarSourceSink)
        np.testing.assert_array_equal(left.values, right.values)
        assert np.all(left.values == 3.0)

    def test_space_content_binding_check(self) -> None:
        """CS4b S3 (F2): twin carriers mix; differing volumes refuse."""
        a = ScalarSourceSink.from_mesh(np.ones((2, 4)), _slab_mesh())
        b = ScalarSourceSink.from_mesh(np.ones((2, 4)), _slab_mesh())
        _ = a + b  # twin content — legal since the F2 re-key
        c = ScalarSourceSink.from_mesh(np.ones((2, 4)), _stretched_mesh())
        with pytest.raises(ValueError, match="equal space"):
            _ = a + c

    def test_linf_and_copy(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(42)
        a = ScalarSourceSink.from_mesh(rng.standard_normal((m.ng, *m.spatial_shape)), m)
        assert a.linf == pytest.approx(float(np.abs(a.values).max()))
        b = a.copy()
        assert b is not a
        assert b.values is not a.values
        np.testing.assert_array_equal(b.values, a.values)


# ════════════════════════════════════════════════════════════════════
# AngularSourceSink — construction + within-type algebra
# ════════════════════════════════════════════════════════════════════


class TestAngularSource:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        Qa = AngularSourceSink.zeros_on(m)
        assert Qa.values.shape == (m.quad.N, m.ng, *m.spatial_shape)
        assert np.all(Qa.values == 0.0)
        assert isinstance(Qa, AngularSourceSink)

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="AngularSourceSink.*does not match"):
            AngularSourceSink.from_mesh(
                np.zeros((m.quad.N + 1, m.ng, *m.spatial_shape)), m,
            )

    def test_within_type_add(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, *m.spatial_shape)
        a = AngularSourceSink.from_mesh(np.ones(shape), m)
        b = AngularSourceSink.from_mesh(2.0 * np.ones(shape), m)
        c = a + b
        assert isinstance(c, AngularSourceSink)
        assert np.all(c.values == 3.0)

    def test_scalar_mul(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, *m.spatial_shape)
        a = AngularSourceSink.from_mesh(np.ones(shape), m)
        c = 2.5 * a
        assert isinstance(c, AngularSourceSink)
        assert np.all(c.values == 2.5)

    def test_at_ordinate_selector(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, *m.spatial_shape)
        rng = np.random.default_rng(7)
        values = rng.standard_normal(shape)
        a = AngularSourceSink.from_mesh(values, m)
        np.testing.assert_array_equal(a.at_ordinate(2), values[2])


# ════════════════════════════════════════════════════════════════════
# Cross-class arithmetic via canonical subspace containment (refined
# Issue #207 principle): ScalarSourceSink ⊂ AngularSourceSink via the
# broadcast injection iso → 1 ⊗ iso. The dunder applies the injection
# inside the operation; result type is the LARGER (containing) space.
# ════════════════════════════════════════════════════════════════════


class TestCrossClassDunder:
    def test_iso_plus_aniso_returns_per_ordinate(self) -> None:
        """``ScalarSourceSink + AngularSourceSink → AngularSourceSink``.

        The canonical subspace-containment dunder: iso is broadcast
        across the Ω axis, then added to the per-ordinate operand.
        Reads as the math ``Q_total = Q_iso + Q_aniso``.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(11)
        iso_values = rng.standard_normal((m.ng, *m.spatial_shape))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape))
        iso = ScalarSourceSink.from_mesh(iso_values, m)
        aniso = AngularSourceSink.from_mesh(aniso_values, m)

        combined = iso + aniso
        assert isinstance(combined, AngularSourceSink)
        # Broadcast: per-ordinate the combined value is iso + aniso[n].
        expected = iso_values[None] + aniso_values
        np.testing.assert_array_equal(combined.values, expected)
        assert combined.mesh is m

    def test_aniso_plus_iso_commutative(self) -> None:
        """``AngularSourceSink + ScalarSourceSink → AngularSourceSink``
        (commutative — delegates to ScalarSourceSink for Pattern 2)."""
        m = _slab_mesh()
        rng = np.random.default_rng(13)
        iso_values = rng.standard_normal((m.ng, *m.spatial_shape))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape))
        iso = ScalarSourceSink.from_mesh(iso_values, m)
        aniso = AngularSourceSink.from_mesh(aniso_values, m)

        forward = iso + aniso
        reverse = aniso + iso
        assert isinstance(reverse, AngularSourceSink)
        np.testing.assert_array_equal(forward.values, reverse.values)

    def test_iso_plus_aniso_coherence_is_the_marginal(self) -> None:
        """CS4b S3 (F2): the cross-class dunder's coherence is the SPACE
        relation — the iso operand must be the angular operand's
        non-angular marginal. Twin carriers satisfy it (the EQUAL leg);
        a volumes-differing carrier breaks it (the UNEQUAL leg)."""
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # distinct instance, same content
        iso = ScalarSourceSink.from_mesh(np.ones((m1.ng, *m1.spatial_shape)), m1)
        aniso = AngularSourceSink.from_mesh(
            np.ones((m2.quad.N, m2.ng, *m2.spatial_shape)), m2,
        )
        _ = iso + aniso  # twin content — legal since the F2 re-key
        m3 = _stretched_mesh()
        aniso3 = AngularSourceSink.from_mesh(
            np.ones((m3.quad.N, m3.ng, *m3.spatial_shape)), m3,
        )
        with pytest.raises(ValueError, match="marginal"):
            _ = iso + aniso3

    def test_pattern_7_normalised_iso_plus_aniso(self) -> None:
        r"""Pattern 7 producer-side normalisation: divide by sum_w
        BEFORE the cross-class add. Mirrors the ScatteringOperator
        call site at ``scattering.py:942``.

        ``(iso / sum_w) + aniso`` reads as the math
        :math:`Q_n = Q_{\text{iso}} / W + Q_{\text{aniso},n}` —
        the per-ordinate magnitude after the producer's dimensional
        normalisation.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(17)
        iso_values = rng.standard_normal((m.ng, *m.spatial_shape))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape))
        iso = ScalarSourceSink.from_mesh(iso_values, m)
        aniso = AngularSourceSink.from_mesh(aniso_values, m)
        sum_w = float(m.quad.weights.sum())

        # Pattern 7 entry point — caller-applied /sum_w.
        combined = (iso / sum_w) + aniso

        # Bit-identical reference: broadcast iso/sum_w into per-ordinate,
        # add aniso element-wise.
        Q_reference = np.broadcast_to(
            (iso_values / sum_w)[None], aniso_values.shape,
        ).copy()
        Q_reference += aniso_values
        np.testing.assert_array_equal(combined.values, Q_reference)

    def test_per_ordinate_from_isotropic_alternative_factory(self) -> None:
        r"""Named factory :meth:`AngularSourceSink.from_isotropic` is
        the alternative to caller-applied /sum_w + cross-class dunder.

        It bakes the ``/sum_w`` normalisation into the factory,
        equivalent to ``(iso / sum_w).as_per_ordinate()``. Caller
        chooses between this and the dunder form based on which
        better reads at the call site.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(23)
        iso_values = rng.standard_normal((m.ng, *m.spatial_shape))
        iso = ScalarSourceSink.from_mesh(iso_values, m)

        via_factory = AngularSourceSink.from_isotropic(iso.values, m)
        via_dunder = (iso / float(m.quad.weights.sum())).as_per_ordinate()
        np.testing.assert_array_equal(via_factory.values, via_dunder.values)

    def test_iso_plus_iso_stays_isotropic(self) -> None:
        """``ScalarSourceSink + ScalarSourceSink → ScalarSourceSink``.

        Within-type addition must STAY within type — the cross-type
        path only fires when the partner is :class:`AngularSourceSink`.
        """
        m = _slab_mesh()
        a = ScalarSourceSink.from_mesh(np.ones((m.ng, *m.spatial_shape)), m)
        b = ScalarSourceSink.from_mesh(2.0 * np.ones((m.ng, *m.spatial_shape)), m)
        c = a + b
        assert isinstance(c, ScalarSourceSink)
        assert not isinstance(c, AngularSourceSink)

    def test_aniso_plus_aniso_within_type(self) -> None:
        """``AngularSourceSink + AngularSourceSink → AngularSourceSink``."""
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, *m.spatial_shape)
        a = AngularSourceSink.from_mesh(np.ones(shape), m)
        b = AngularSourceSink.from_mesh(2.0 * np.ones(shape), m)
        c = a + b
        assert isinstance(c, AngularSourceSink)
        assert np.all(c.values == 3.0)


# ════════════════════════════════════════════════════════════════════
# ScalarSourceSink.as_per_ordinate() — explicit broadcast conversion
# ════════════════════════════════════════════════════════════════════


class TestAsPerOrdinate:
    def test_round_trip_shape_and_values(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(19)
        iso_values = rng.standard_normal((m.ng, *m.spatial_shape))
        iso = ScalarSourceSink.from_mesh(iso_values, m)

        aniso = iso.as_per_ordinate()
        assert isinstance(aniso, AngularSourceSink)
        assert aniso.values.shape == (m.quad.N, m.ng, *m.spatial_shape)
        # Every ordinate slice IS the isotropic source.
        for n in range(m.quad.N):
            np.testing.assert_array_equal(aniso.values[n], iso_values)

    def test_as_per_ordinate_owns_data(self) -> None:
        """The returned array must own its data (writable) — the
        ``.copy()`` in :meth:`ScalarSourceSink.as_per_ordinate` is the
        load-bearing piece that prevents readonly-broadcast-view
        surprises at the caller."""
        m = _slab_mesh()
        iso = ScalarSourceSink.from_mesh(np.ones((m.ng, *m.spatial_shape)), m)
        aniso = iso.as_per_ordinate()
        assert aniso.values.flags.writeable

    def test_as_per_ordinate_equivalent_to_dunder_with_zero(self) -> None:
        """``iso.as_per_ordinate() == iso + AngularSourceSink.zeros_on()``.

        The explicit conversion and the dunder-against-zero produce
        equivalent results (modulo the type of the zero partner).
        """
        m = _slab_mesh()
        rng = np.random.default_rng(23)
        iso = ScalarSourceSink.from_mesh(
            rng.standard_normal((m.ng, *m.spatial_shape)), m,
        )
        via_conv = iso.as_per_ordinate()
        via_dunder = iso + AngularSourceSink.zeros_on(m)
        np.testing.assert_array_equal(via_conv.values, via_dunder.values)


# ════════════════════════════════════════════════════════════════════
# 2-D mesh smoke — typed sources work for 2-D as well
# ════════════════════════════════════════════════════════════════════


class Test2DTypedSources:
    def test_factory_shapes_2d(self) -> None:
        m = _2d_mesh()
        iso = ScalarSourceSink.zeros_on(m)
        aniso = AngularSourceSink.zeros_on(m)
        assert iso.values.shape == (m.ng, *m.spatial_shape)
        assert aniso.values.shape == (m.quad.N, m.ng, *m.spatial_shape)

    def test_cross_type_add_2d(self) -> None:
        m = _2d_mesh()
        rng = np.random.default_rng(29)
        iso = ScalarSourceSink.from_mesh(
            rng.standard_normal((m.ng, *m.spatial_shape)), m,
        )
        aniso = AngularSourceSink.from_mesh(
            rng.standard_normal((m.quad.N, m.ng, *m.spatial_shape)), m,
        )
        combined = iso + aniso
        expected = iso.values[None] + aniso.values
        np.testing.assert_array_equal(combined.values, expected)


# ════════════════════════════════════════════════════════════════════
# HarmonicMomentSourceSink — the moment-space source/sink leaf (P4).
#
# Intrinsic-property test (feedback_test_intrinsic_properties): the
# leaf's DEFINING law is the bare-MomentField *vector-space* algebra —
# `source ± source` is CLOSED (a rate density adds vectorially). (Until
# campaign 1 CS3, 2026-08-19, this stood in deliberate contrast to the
# then-FluxRole sibling HarmonicMomentFlux, whose `flux + flux` raised
# and whose `flux − flux` minted a displacement; since the cone carve
# both families carry the same plain V algebra and the contrast is
# history — the cross-CLASS gates below are what still discriminate.)
# Positive AND negative cases (vv-principles L11 / Mode
# 11 — never negative-only): the positives pin the closed algebra, the
# negatives pin the class-identity / L / mesh gates.
# ════════════════════════════════════════════════════════════════════


def _moment_shape(m: SNMesh, L: int) -> tuple[int, ...]:
    """The ``(L+1, 2L+1, ng, *spatial)`` moment layout for mesh ``m``."""
    return (L + 1, 2 * L + 1, m.ng, *m.spatial_shape)


class TestHarmonicMomentSourceSink:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        q = HarmonicMomentSourceSink.zeros_for_mesh_and_L(m, L=2)
        assert q.values.shape == _moment_shape(m, 2)
        assert np.all(q.values == 0.0)
        assert isinstance(q, HarmonicMomentSourceSink)
        assert q.mesh is m
        assert q.L == 2

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        bad = np.zeros((2 + 1, 2 * 2 + 1, m.ng + 1, *m.spatial_shape))
        with pytest.raises(
            ValueError, match="HarmonicMomentSourceSink.*does not match",
        ):
            HarmonicMomentSourceSink.from_mesh_and_L(bad, m, L=2)

    def test_space_is_shared_across_the_role_leaves(self) -> None:
        """CS4b S2 (F1-sub): the two role leaves compose the SAME space —
        ``SphericalHarmonicSpace(L) * mesh.bulk_space``, the carrier's one
        cached cell-group mint — so role lives ONLY in the class, and the
        class arm is the sole role guard
        (``test_cross_class_with_flux_rejected`` below pins its raise).
        Until CS4b each leaf minted a role-named cell-group tag and this
        test asserted the spaces UNEQUAL; that defensive duplication of
        the class gate retired with ``_CELL_GROUP_NAME``."""
        m = _slab_mesh()
        q = HarmonicMomentSourceSink.zeros_for_mesh_and_L(m, L=2)
        phi = HarmonicMomentFlux.zeros_for_mesh_and_L(m, L=2)
        assert q.values.shape == phi.values.shape  # same layout
        assert q.space == phi.space  # ONE space; role is class identity

    # ── The DEFINING law: closed vector-space algebra (positive) ─────

    def test_within_class_add_is_closed(self) -> None:
        """``source + source → source`` — CLOSED (a rate density adds
        vectorially). The load-bearing contrast with the flux leaf."""
        m = _slab_mesh()
        shape = _moment_shape(m, 2)
        a = HarmonicMomentSourceSink.from_mesh_and_L(np.ones(shape), m, L=2)
        b = HarmonicMomentSourceSink.from_mesh_and_L(2.0 * np.ones(shape), m, L=2)
        c = a + b
        assert isinstance(c, HarmonicMomentSourceSink)
        assert np.all(c.values == 3.0)
        # frozen — operands untouched
        assert np.all(a.values == 1.0)

    def test_within_class_sub_is_closed(self) -> None:
        """``source − source → source`` — a plain vector difference (since
        the CS3 cone carve the flux sibling behaves the same way; the
        surviving distinction is the CROSS-CLASS gate below)."""
        m = _slab_mesh()
        shape = _moment_shape(m, 2)
        a = HarmonicMomentSourceSink.from_mesh_and_L(3.0 * np.ones(shape), m, L=2)
        b = HarmonicMomentSourceSink.from_mesh_and_L(np.ones(shape), m, L=2)
        c = a - b
        assert isinstance(c, HarmonicMomentSourceSink)
        assert np.all(c.values == 2.0)

    def test_scalar_mul_left_and_right(self) -> None:
        m = _slab_mesh()
        shape = _moment_shape(m, 1)
        a = HarmonicMomentSourceSink.from_mesh_and_L(np.ones(shape), m, L=1)
        left, right = 3.0 * a, a * 3.0
        assert isinstance(left, HarmonicMomentSourceSink)
        np.testing.assert_array_equal(left.values, right.values)
        assert np.all(left.values == 3.0)

    def test_flux_vs_source_roles_share_algebra_but_never_mix(self) -> None:
        """Re-posed at the CS3 cone carve (2026-08-19): both role families
        now carry the SAME vector additive algebra (flux lives in V — the
        old contrast row asserted the flux side FORBIDS ``+``), so the
        surviving load-bearing distinction is the CROSS-CLASS gate: a state
        and a rate density never mix, even at equal shape and units
        (``test_cross_class_with_flux_rejected`` below is the sharp form)."""
        m = _slab_mesh()
        shape = _moment_shape(m, 2)
        phi_a = HarmonicMomentFlux.from_mesh_and_L(np.ones(shape), m, L=2)
        phi_b = HarmonicMomentFlux.from_mesh_and_L(np.ones(shape), m, L=2)
        s = phi_a + phi_b
        if type(s) is not HarmonicMomentFlux:
            raise AssertionError("flux + flux left the leaf type")
        q_a = HarmonicMomentSourceSink.from_mesh_and_L(np.ones(shape), m, L=2)
        q_b = HarmonicMomentSourceSink.from_mesh_and_L(np.ones(shape), m, L=2)
        if not isinstance(q_a + q_b, HarmonicMomentSourceSink):
            raise AssertionError("source + source left the leaf type")
        with pytest.raises(TypeError):  # the surviving distinction
            _ = phi_a + q_a  # type: ignore[operator]

    # ── The gates: class identity, L-match, mesh-binding (negative) ──

    def test_cross_class_with_flux_rejected(self) -> None:
        """``source + flux`` (same shape, different role) → TypeError via
        Field's Layer-1 class-identity gate."""
        m = _slab_mesh()
        shape = _moment_shape(m, 2)
        q = HarmonicMomentSourceSink.from_mesh_and_L(np.ones(shape), m, L=2)
        phi = HarmonicMomentFlux.from_mesh_and_L(np.ones(shape), m, L=2)
        with pytest.raises(TypeError, match="same-class partner"):
            _ = q + phi

    def test_cross_class_with_angular_source_rejected(self) -> None:
        """``source + AngularSourceSink`` (different storage family) →
        TypeError; there is no containment injection for the moment leaf."""
        m = _slab_mesh()
        q = HarmonicMomentSourceSink.from_mesh_and_L(
            np.ones(_moment_shape(m, 2)), m, L=2,
        )
        ang = AngularSourceSink.from_mesh(
            np.ones((m.quad.N, m.ng, *m.spatial_shape)), m,
        )
        with pytest.raises(TypeError, match="same-class partner"):
            _ = q + ang

    def test_cross_class_with_scalar_source_rejected(self) -> None:
        """``source + ScalarSourceSink`` → TypeError — unlike the
        Scalar↔Angular pair, the moment leaf carries no subspace-
        containment dunder (no consumer yet; build-the-seam-not-spec)."""
        m = _slab_mesh()
        q = HarmonicMomentSourceSink.from_mesh_and_L(
            np.ones(_moment_shape(m, 2)), m, L=2,
        )
        sca = ScalarSourceSink.from_mesh(np.ones((m.ng, *m.spatial_shape)), m)
        with pytest.raises(TypeError, match="same-class partner"):
            _ = q + sca

    def test_L_mismatch_rejected(self) -> None:
        """``source(L=1) + source(L=2)`` → ValueError.

        Different ``L`` produces different ``SphericalHarmonicSpace`` shapes,
        so Field's space-equality check fires first with the more-general
        "requires equal space" message (the explicit L-match in
        ``MomentField._check_partner`` is the documented belt-and-suspenders
        secondary gate — same as the flux leaf, see
        ``test_harmonic_moment_flux.py``). Both gate the same invariant.
        """
        m = _slab_mesh()
        a = HarmonicMomentSourceSink.from_mesh_and_L(
            np.ones(_moment_shape(m, 1)), m, L=1,
        )
        b = HarmonicMomentSourceSink.from_mesh_and_L(
            np.ones(_moment_shape(m, 2)), m, L=2,
        )
        with pytest.raises(ValueError, match="equal space"):
            _ = a + b

    def test_space_content_binding_rejected(self) -> None:
        """CS4b S3 (F2): twin carriers mint EQUAL moment spaces (the SH
        factor and the shared scalar-bulk cell-group agree); a carrier
        with differing volumes mints an UNEQUAL cell-group factor and the
        sum refuses on space content."""
        m1, m2 = _slab_mesh(), _stretched_mesh()
        a = HarmonicMomentSourceSink.from_mesh_and_L(
            np.ones(_moment_shape(m1, 2)), m1, L=2,
        )
        b = HarmonicMomentSourceSink.from_mesh_and_L(
            np.ones(_moment_shape(m2, 2)), m2, L=2,
        )
        with pytest.raises(ValueError, match="equal space"):
            _ = a + b

    def test_2d_construction(self) -> None:
        m = _2d_mesh()
        q = HarmonicMomentSourceSink.zeros_for_mesh_and_L(m, L=2)
        assert q.values.shape == _moment_shape(m, 2)
        assert isinstance(q, HarmonicMomentSourceSink)
