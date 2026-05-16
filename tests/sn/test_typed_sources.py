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
from orpheus.sn.quadrature import GaussLegendre1D, LevelSymmetricSN
from orpheus.sn.sources import IsotropicSource, PerOrdinateSource

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
    quad = GaussLegendre1D.create(n_ordinates=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _2d_mesh(nx: int = 3, ny: int = 3, ng: int = 1) -> SNMesh:
    """Build a small 2-D Cartesian :class:`SNMesh`."""
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
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
        with pytest.raises(ValueError, match="IsotropicSource expects"):
            IsotropicSource(np.zeros((m.ng + 1, m.nx, m.ny)), m)

    def test_within_type_add(self) -> None:
        m = _slab_mesh()
        a = IsotropicSource(np.ones((m.ng, m.nx, m.ny)), m)
        b = IsotropicSource(2.0 * np.ones((m.ng, m.nx, m.ny)), m)
        c = a + b
        assert isinstance(c, IsotropicSource)
        assert np.all(c.values == 3.0)
        # Originals untouched (frozen)
        assert np.all(a.values == 1.0)
        assert np.all(b.values == 2.0)

    def test_within_type_sub(self) -> None:
        m = _slab_mesh()
        a = IsotropicSource(3.0 * np.ones((m.ng, m.nx, m.ny)), m)
        b = IsotropicSource(np.ones((m.ng, m.nx, m.ny)), m)
        c = a - b
        assert isinstance(c, IsotropicSource)
        assert np.all(c.values == 2.0)

    def test_scalar_mul_left_and_right(self) -> None:
        m = _slab_mesh()
        a = IsotropicSource(np.ones((m.ng, m.nx, m.ny)), m)
        left = 3.0 * a
        right = a * 3.0
        assert isinstance(left, IsotropicSource)
        assert isinstance(right, IsotropicSource)
        np.testing.assert_array_equal(left.values, right.values)
        assert np.all(left.values == 3.0)

    def test_mesh_binding_check(self) -> None:
        a = IsotropicSource(np.ones((2, 4, 1)), _slab_mesh())
        b = IsotropicSource(np.ones((2, 4, 1)), _slab_mesh())  # different mesh instance
        with pytest.raises(ValueError, match="distinct SNMesh"):
            _ = a + b

    def test_linf_and_copy(self) -> None:
        m = _slab_mesh()
        rng = np.random.default_rng(42)
        a = IsotropicSource(rng.standard_normal((m.ng, m.nx, m.ny)), m)
        assert a.linf() == pytest.approx(float(np.abs(a.values).max()))
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
        with pytest.raises(ValueError, match="PerOrdinateSource expects"):
            PerOrdinateSource(
                np.zeros((m.quad.N + 1, m.ng, m.nx, m.ny)), m,
            )

    def test_within_type_add(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        a = PerOrdinateSource(np.ones(shape), m)
        b = PerOrdinateSource(2.0 * np.ones(shape), m)
        c = a + b
        assert isinstance(c, PerOrdinateSource)
        assert np.all(c.values == 3.0)

    def test_scalar_mul(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        a = PerOrdinateSource(np.ones(shape), m)
        c = 2.5 * a
        assert isinstance(c, PerOrdinateSource)
        assert np.all(c.values == 2.5)

    def test_at_ordinate_selector(self) -> None:
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        rng = np.random.default_rng(7)
        values = rng.standard_normal(shape)
        a = PerOrdinateSource(values, m)
        np.testing.assert_array_equal(a.at_ordinate(2), values[2])


# ════════════════════════════════════════════════════════════════════
# Cross-type dunder algebra (the load-bearing PR-TYPED-3 pattern)
# ════════════════════════════════════════════════════════════════════


class TestCrossTypeDunder:
    def test_iso_plus_aniso_returns_per_ordinate(self) -> None:
        """``IsotropicSource + PerOrdinateSource → PerOrdinateSource``.

        The load-bearing case for PR-TYPED-3 — replaces the procedural
        ``np.broadcast_to(...).copy(); Q += Q_aniso`` pattern with one
        line of typed algebra.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(11)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny))
        iso = IsotropicSource(iso_values, m)
        aniso = PerOrdinateSource(aniso_values, m)

        combined = iso + aniso
        assert isinstance(combined, PerOrdinateSource)
        # Broadcast: per-ordinate the combined value is iso + aniso[n].
        expected = iso_values[None, :, :, :] + aniso_values
        np.testing.assert_array_equal(combined.values, expected)
        # Mesh propagates through the dunder.
        assert combined.mesh is m

    def test_aniso_plus_iso_commutative(self) -> None:
        """``PerOrdinateSource + IsotropicSource → PerOrdinateSource``.

        The reverse-order case must yield the same result (algebra is
        commutative by construction — the dunder delegates to the
        IsotropicSource path).
        """
        m = _slab_mesh()
        rng = np.random.default_rng(13)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny))
        iso = IsotropicSource(iso_values, m)
        aniso = PerOrdinateSource(aniso_values, m)

        forward = iso + aniso
        reverse = aniso + iso
        assert isinstance(reverse, PerOrdinateSource)
        np.testing.assert_array_equal(forward.values, reverse.values)

    def test_bit_identity_with_legacy_broadcast_pattern(self) -> None:
        r"""``(iso + aniso).values == np.broadcast_to(iso.values[None, ...], aniso.values.shape).copy() + aniso.values``.

        Proves the dunder replicates the procedural pattern at FP
        bit-identity (the legacy ``ScatteringOperator.apply`` pattern
        ``Q = np.broadcast_to(...).copy(); Q += Q_aniso``).  Both
        produce a newly-allocated ``(N, ng, nx, ny)`` array; numpy's
        broadcast add yields the same byte-for-byte result.
        """
        m = _slab_mesh()
        rng = np.random.default_rng(17)
        iso_values = rng.standard_normal((m.ng, m.nx, m.ny))
        aniso_values = rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny))
        iso = IsotropicSource(iso_values, m)
        aniso = PerOrdinateSource(aniso_values, m)

        # Legacy pattern (from scattering.py:799 pre-PR-TYPED-3).
        Q_legacy = np.broadcast_to(
            iso_values[None, :, :, :], aniso_values.shape,
        ).copy()
        Q_legacy += aniso_values

        # New typed dunder.
        Q_typed = (iso + aniso).values

        np.testing.assert_array_equal(Q_typed, Q_legacy)

    def test_iso_plus_iso_stays_isotropic(self) -> None:
        """``IsotropicSource + IsotropicSource → IsotropicSource``.

        Within-type addition must STAY within type — the cross-type
        path only fires when the partner is :class:`PerOrdinateSource`.
        """
        m = _slab_mesh()
        a = IsotropicSource(np.ones((m.ng, m.nx, m.ny)), m)
        b = IsotropicSource(2.0 * np.ones((m.ng, m.nx, m.ny)), m)
        c = a + b
        assert isinstance(c, IsotropicSource)
        assert not isinstance(c, PerOrdinateSource)

    def test_aniso_plus_aniso_within_type(self) -> None:
        """``PerOrdinateSource + PerOrdinateSource → PerOrdinateSource``."""
        m = _slab_mesh()
        shape = (m.quad.N, m.ng, m.nx, m.ny)
        a = PerOrdinateSource(np.ones(shape), m)
        b = PerOrdinateSource(2.0 * np.ones(shape), m)
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
        iso = IsotropicSource(iso_values, m)

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
        iso = IsotropicSource(np.ones((m.ng, m.nx, m.ny)), m)
        aniso = iso.as_per_ordinate()
        assert aniso.values.flags.writeable

    def test_as_per_ordinate_equivalent_to_dunder_with_zero(self) -> None:
        """``iso.as_per_ordinate() == iso + zeros_per_ordinate_source()``.

        The explicit conversion and the dunder-against-zero produce
        equivalent results (modulo the type of the zero partner).
        """
        m = _slab_mesh()
        rng = np.random.default_rng(23)
        iso = IsotropicSource(
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
        iso = IsotropicSource(
            rng.standard_normal((m.ng, m.nx, m.ny)), m,
        )
        aniso = PerOrdinateSource(
            rng.standard_normal((m.quad.N, m.ng, m.nx, m.ny)), m,
        )
        combined = iso + aniso
        expected = iso.values[None, :, :, :] + aniso.values
        np.testing.assert_array_equal(combined.values, expected)
