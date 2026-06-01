"""Foundation tests for Issue #197 PR-TYPED-3 typed source types.

Pins the structural contract of :class:`IsotropicSource` and
:class:`PerOrdinateSource` — shape validation, dunder algebra
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
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from orpheus.transport.sources import IsotropicSource, PerOrdinateSource

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
# IsotropicSource — construction + within-type algebra
# ════════════════════════════════════════════════════════════════════


class TestIsotropicSource:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        Q = m.zeros_isotropic_source()
        assert Q.values.shape == (m.ng, m.nx, m.ny)
        assert np.all(Q.values == 0.0)
        assert isinstance(Q, IsotropicSource)
        assert Q.mesh is m

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="IsotropicSource.*does not match"):
            IsotropicSource.from_mesh(np.zeros((m.ng + 1, m.nx, m.ny)), m)

    def test_within_type_add(self) -> None:
        m = _slab_mesh()
        a = IsotropicSource.from_mesh(np.ones((m.ng, m.nx, m.ny)), m)
        b = IsotropicSource.from_mesh(2.0 * np.ones((m.ng, m.nx, m.ny)), m)
        c = a + b
        assert isinstance(c, IsotropicSource)
        assert np.all(c.values == 3.0)
        # Originals untouched (frozen)
        assert np.all(a.values == 1.0)
        assert np.all(b.values == 2.0)

    def test_within_type_sub(self) -> None:
        m = _slab_mesh()
        a = IsotropicSource.from_mesh(3.0 * np.ones((m.ng, m.nx, m.ny)), m)
        b = IsotropicSource.from_mesh(np.ones((m.ng, m.nx, m.ny)), m)
        c = a - b
        assert isinstance(c, IsotropicSource)
        assert np.all(c.values == 2.0)

    def test_scalar_mul_left_and_right(self) -> None:
        m = _slab_mesh()
        a = IsotropicSource.from_mesh(np.ones((m.ng, m.nx, m.ny)), m)
        left = 3.0 * a
        right = a * 3.0
        assert isinstance(left, IsotropicSource)
        assert isinstance(right, IsotropicSource)
        np.testing.assert_array_equal(left.values, right.values)
        assert np.all(left.values == 3.0)

    def test_mesh_binding_check(self) -> None:
        a = IsotropicSource.from_mesh(np.ones((2, 4, 1)), _slab_mesh())
        b = IsotropicSource.from_mesh(np.ones((2, 4, 1)), _slab_mesh())  # different mesh instance
        with pytest.raises(ValueError, match="distinct SNMesh"):
            _ = a + b

    def test_linf_and_copy(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(42)
        a = IsotropicSource.from_mesh(rng.standard_normal((m.ng, m.nx, m.ny)), m)
        assert a.linf == pytest.approx(float(np.abs(a.values).max()))
        b = a.copy()
        assert b is not a
        assert b.values is not a.values
        np.testing.assert_array_equal(b.values, a.values)


# ════════════════════════════════════════════════════════════════════
# PerOrdinateSource — construction + within-type algebra
# ════════════════════════════════════════════════════════════════════


class TestPerOrdinateSource:
    def test_construct_from_factory(self) -> None:
        m = _slab_mesh()
        Qa = m.zeros_per_ordinate_source()
        assert Qa.values.shape == (m.quad.N, m.ng, m.nx, m.ny)
        assert np.all(Qa.values == 0.0)
        assert isinstance(Qa, PerOrdinateSource)

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="PerOrdinateSource.*does not match"):
            PerOrdinateSource.from_mesh(
                np.zeros((m.quad.N + 1, m.ng, m.nx, m.ny)), m,
            )

    def test_within_type_add(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        a = PerOrdinateSource.from_mesh(np.ones(shape), m)
        b = PerOrdinateSource.from_mesh(2.0 * np.ones(shape), m)
        c = a + b
        assert isinstance(c, PerOrdinateSource)
        assert np.all(c.values == 3.0)

    def test_scalar_mul(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        a = PerOrdinateSource.from_mesh(np.ones(shape), m)
        c = 2.5 * a
        assert isinstance(c, PerOrdinateSource)
        assert np.all(c.values == 2.5)

    def test_at_ordinate_selector(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        rng = np.random.default_rng(7)
        values = rng.standard_normal(shape)
        a = PerOrdinateSource.from_mesh(values, m)
        np.testing.assert_array_equal(a.at_ordinate(2), values[2])


# ════════════════════════════════════════════════════════════════════
# Cross-class arithmetic via canonical subspace containment (refined
# Issue #207 principle): IsotropicSource ⊂ PerOrdinateSource via the
# broadcast injection iso → 1 ⊗ iso. The dunder applies the injection
# inside the operation; result type is the LARGER (containing) space.
# ════════════════════════════════════════════════════════════════════


class TestCrossClassDunder:
    def test_iso_plus_aniso_returns_per_ordinate(self) -> None:
        """``IsotropicSource + PerOrdinateSource → PerOrdinateSource``.

        The canonical subspace-containment dunder: iso is broadcast
        across the Ω axis, then added to the per-ordinate operand.
        Reads as the math ``Q_total = Q_iso + Q_aniso``.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(11)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny))
        iso = IsotropicSource.from_mesh(iso_values, m)
        aniso = PerOrdinateSource.from_mesh(aniso_values, m)

        combined = iso + aniso
        assert isinstance(combined, PerOrdinateSource)
        # Broadcast: per-ordinate the combined value is iso + aniso[n].
        expected = iso_values[None, :, :, :] + aniso_values
        np.testing.assert_array_equal(combined.values, expected)
        assert combined.mesh is m

    def test_aniso_plus_iso_commutative(self) -> None:
        """``PerOrdinateSource + IsotropicSource → PerOrdinateSource``
        (commutative — delegates to IsotropicSource for Pattern 2)."""
        m = _slab_mesh()
        rng = np.random.default_rng(13)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny))
        iso = IsotropicSource.from_mesh(iso_values, m)
        aniso = PerOrdinateSource.from_mesh(aniso_values, m)

        forward = iso + aniso
        reverse = aniso + iso
        assert isinstance(reverse, PerOrdinateSource)
        np.testing.assert_array_equal(forward.values, reverse.values)

    def test_iso_plus_aniso_rejects_distinct_meshes(self) -> None:
        """The cross-class dunder still enforces mesh-identity."""
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # distinct instance, same shape
        iso = IsotropicSource.from_mesh(np.ones((m1.ng, m1.nx, m1.ny)), m1)
        aniso = PerOrdinateSource.from_mesh(
            np.ones((m2.quad.N, m2.ng, m2.nx, m2.ny)), m2,
        )
        with pytest.raises(ValueError, match="mesh-bound"):
            _ = iso + aniso

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
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny))
        iso = IsotropicSource.from_mesh(iso_values, m)
        aniso = PerOrdinateSource.from_mesh(aniso_values, m)
        sum_w = float(m.quad.weights.sum())

        # Pattern 7 entry point — caller-applied /sum_w.
        combined = (iso / sum_w) + aniso

        # Bit-identical reference: broadcast iso/sum_w into per-ordinate,
        # add aniso element-wise.
        Q_reference = np.broadcast_to(
            (iso_values / sum_w)[None, :, :, :], aniso_values.shape,
        ).copy()
        Q_reference += aniso_values
        np.testing.assert_array_equal(combined.values, Q_reference)

    def test_per_ordinate_from_isotropic_alternative_factory(self) -> None:
        r"""Named factory :meth:`PerOrdinateSource.from_isotropic` is
        the alternative to caller-applied /sum_w + cross-class dunder.

        It bakes the ``/sum_w`` normalisation into the factory,
        equivalent to ``(iso / sum_w).as_per_ordinate()``. Caller
        chooses between this and the dunder form based on which
        better reads at the call site.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(23)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        iso = IsotropicSource.from_mesh(iso_values, m)

        via_factory = PerOrdinateSource.from_isotropic(iso.values, m)
        via_dunder = (iso / float(m.quad.weights.sum())).as_per_ordinate()
        np.testing.assert_array_equal(via_factory.values, via_dunder.values)

    def test_iso_plus_iso_stays_isotropic(self) -> None:
        """``IsotropicSource + IsotropicSource → IsotropicSource``.

        Within-type addition must STAY within type — the cross-type
        path only fires when the partner is :class:`PerOrdinateSource`.
        """
        m = _slab_mesh()
        a = IsotropicSource.from_mesh(np.ones((m.ng, m.nx, m.ny)), m)
        b = IsotropicSource.from_mesh(2.0 * np.ones((m.ng, m.nx, m.ny)), m)
        c = a + b
        assert isinstance(c, IsotropicSource)
        assert not isinstance(c, PerOrdinateSource)

    def test_aniso_plus_aniso_within_type(self) -> None:
        """``PerOrdinateSource + PerOrdinateSource → PerOrdinateSource``."""
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        a = PerOrdinateSource.from_mesh(np.ones(shape), m)
        b = PerOrdinateSource.from_mesh(2.0 * np.ones(shape), m)
        c = a + b
        assert isinstance(c, PerOrdinateSource)
        assert np.all(c.values == 3.0)


# ════════════════════════════════════════════════════════════════════
# IsotropicSource.as_per_ordinate() — explicit broadcast conversion
# ════════════════════════════════════════════════════════════════════


class TestAsPerOrdinate:
    def test_round_trip_shape_and_values(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(19)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        iso = IsotropicSource.from_mesh(iso_values, m)

        aniso = iso.as_per_ordinate()
        assert isinstance(aniso, PerOrdinateSource)
        assert aniso.values.shape == (m.quad.N, m.ng, m.nx, m.ny)
        # Every ordinate slice IS the isotropic source.
        for n in range(m.quad.N):
            np.testing.assert_array_equal(aniso.values[n], iso_values)

    def test_as_per_ordinate_owns_data(self) -> None:
        """The returned array must own its data (writable) — the
        ``.copy()`` in :meth:`IsotropicSource.as_per_ordinate` is the
        load-bearing piece that prevents readonly-broadcast-view
        surprises at the caller."""
        m = _slab_mesh()
        iso = IsotropicSource.from_mesh(np.ones((m.ng, m.nx, m.ny)), m)
        aniso = iso.as_per_ordinate()
        assert aniso.values.flags.writeable

    def test_as_per_ordinate_equivalent_to_dunder_with_zero(self) -> None:
        """``iso.as_per_ordinate() == iso + zeros_per_ordinate_source()``.

        The explicit conversion and the dunder-against-zero produce
        equivalent results (modulo the type of the zero partner).
        """
        m = _slab_mesh()
        rng = np.random.default_rng(23)
        iso = IsotropicSource.from_mesh(
            rng.standard_normal((m.ng, m.nx, m.ny)), m,
        )
        via_conv = iso.as_per_ordinate()
        via_dunder = iso + m.zeros_per_ordinate_source()
        np.testing.assert_array_equal(via_conv.values, via_dunder.values)


# ════════════════════════════════════════════════════════════════════
# 2-D mesh smoke — typed sources work for 2-D as well
# ════════════════════════════════════════════════════════════════════


class Test2DTypedSources:
    def test_factory_shapes_2d(self) -> None:
        m = _2d_mesh()
        iso = m.zeros_isotropic_source()
        aniso = m.zeros_per_ordinate_source()
        assert iso.values.shape == (m.ng, m.nx, m.ny)
        assert aniso.values.shape == (m.quad.N, m.ng, m.nx, m.ny)

    def test_cross_type_add_2d(self) -> None:
        m = _2d_mesh()
        rng = np.random.default_rng(29)
        iso = IsotropicSource.from_mesh(
            rng.standard_normal((m.ng, m.nx, m.ny)), m,
        )
        aniso = PerOrdinateSource.from_mesh(
            rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny)), m,
        )
        combined = iso + aniso
        expected = iso.values[None, :, :, :] + aniso.values
        np.testing.assert_array_equal(combined.values, expected)
