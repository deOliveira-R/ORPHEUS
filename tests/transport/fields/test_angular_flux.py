r"""L0/L2 — pure-Field :class:`AngularFlux` foundation tests.

Pins the post-D-H.1 contract for
:class:`orpheus.transport.fields.angular_flux.AngularFlux`:

* Inherits :class:`~orpheus.numerics.field.Field`; no hand-coded
  dunders.
* Pure Field — NO ``boundary`` attribute, NO ``_history`` attribute.
* Mesh-binding rejection (cross-mesh arithmetic raises).
* ``from_mesh`` / ``from_ndarray`` / ``zeros_on`` factories.
* Frozen contract.
* ``integrate_angular`` reduction to ScalarFlux.

Symmetric with :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
algebra (no special-case for boundary or history).

References
----------

* Depth B plan §3.3 (pure AngularFlux at L2).
* Pre-D-H.1 ``orpheus/sn/angular_flux.py:AngularFlux`` (legacy with
  ``boundary`` + ``_history`` fields — retired in D-H.2).
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.field import Field
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux

from tests.sn._test_helpers import placeholder_materials


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Mesh fixtures
# ───────────────────────────────────────────────────────────────────────


def _slab_mesh(nx: int = 4, ng: int = 2) -> SNMesh:
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


def _cartesian_2d_mesh(nx: int = 3, ny: int = 2, ng: int = 2) -> SNMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# A. Field-ABC inheritance
# ───────────────────────────────────────────────────────────────────────


class TestFieldInheritance:
    def test_inherits_field(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        assert isinstance(psi, Field)

    def test_pure_field_no_boundary_attribute(self) -> None:
        """Critical: post-D-H.1 AngularFlux has NO .boundary field."""
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        assert not hasattr(psi, "boundary")

    def test_pure_field_no_history_attribute(self) -> None:
        """Critical: post-D-H.1 AngularFlux has NO ._history field."""
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        assert not hasattr(psi, "_history")

    def test_pure_field_no_history_depth_attribute(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        assert not hasattr(psi, "history_depth")


# ───────────────────────────────────────────────────────────────────────
# B. Algebra (inherited from Field)
# ───────────────────────────────────────────────────────────────────────


class TestAlgebra:
    def test_flux_add_flux_legal_and_update_round_trip(self) -> None:
        """Flux lives in V (campaign 1 CS3): ψ + ψ' is the plain vector sum
        in the same leaf type, and the update step ψ + (ψ' − ψ) ≈ ψ' is
        ordinary arithmetic (until 2026-08-19 the #208 affine gate raised on
        the first and typed the second through a displacement mint)."""
        m = _slab_mesh()
        a = AngularFlux.from_mesh(np.ones((m.quad.N, m.ng, *m.spatial_shape)), m)
        b = AngularFlux.from_mesh(2.0 * np.ones((m.quad.N, m.ng, *m.spatial_shape)), m)
        s = a + b
        if type(s) is not AngularFlux:
            raise AssertionError("flux + flux left the leaf type")
        np.testing.assert_array_equal(s.values, 3.0)
        out = a + (b - a)  # the update step, plain V arithmetic
        if type(out) is not AngularFlux:
            raise AssertionError("the update step left the leaf type")
        np.testing.assert_array_almost_equal_nulp(out.values, b.values, nulp=4)

    def test_sub(self) -> None:
        m = _slab_mesh()
        a = AngularFlux.from_mesh(np.full((m.quad.N, m.ng, *m.spatial_shape), 5.0), m)
        b = AngularFlux.from_mesh(np.full((m.quad.N, m.ng, *m.spatial_shape), 2.0), m)
        c = a - b
        np.testing.assert_array_equal(c.values, 3.0)

    def test_scalar_mul(self) -> None:
        m = _slab_mesh()
        a = AngularFlux.from_mesh(np.full((m.quad.N, m.ng, *m.spatial_shape), 2.0), m)
        c = a * 3.0
        np.testing.assert_array_equal(c.values, 6.0)

    def test_scalar_div(self) -> None:
        m = _slab_mesh()
        a = AngularFlux.from_mesh(np.full((m.quad.N, m.ng, *m.spatial_shape), 8.0), m)
        c = a / 4.0
        np.testing.assert_array_equal(c.values, 2.0)

    def test_neg(self) -> None:
        m = _slab_mesh()
        a = AngularFlux.from_mesh(np.full((m.quad.N, m.ng, *m.spatial_shape), 1.5), m)
        c = -a
        np.testing.assert_array_equal(c.values, -1.5)


# ───────────────────────────────────────────────────────────────────────
# C. Mesh-binding rejection
# ───────────────────────────────────────────────────────────────────────


class TestMeshBinding:
    def test_cross_carrier_discrimination_is_space_content(self) -> None:
        """CS4b S3 (F2): twin carriers mint EQUAL spaces and mix; a carrier
        with different cell volumes mints an UNEQUAL space and refuses."""
        m1 = _slab_mesh()
        m2 = _slab_mesh()  # different instance, same content
        a = AngularFlux.zeros_on(m1)
        b = AngularFlux.zeros_on(m2)
        _ = a - b  # twin content — legal since the F2 re-key
        _ = a + b
        c = AngularFlux.zeros_on(_stretched_mesh())
        with pytest.raises(ValueError, match="equal space"):
            a - c
        with pytest.raises(ValueError, match="equal space"):
            a + c

    def test_cross_class_rejected(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        phi = ScalarFlux.from_mesh(np.zeros((m.ng, *m.spatial_shape)), m)
        with pytest.raises(TypeError, match="same-class"):
            psi + phi  # type: ignore[operator]


# ───────────────────────────────────────────────────────────────────────
# D. Construction factories
# ───────────────────────────────────────────────────────────────────────


class TestConstruction:
    def test_from_mesh(self) -> None:
        m = _slab_mesh()
        arr = np.ones((m.quad.N, m.ng, *m.spatial_shape))
        psi = AngularFlux.from_mesh(arr, m)
        assert psi.values.shape == (m.quad.N, m.ng, *m.spatial_shape)
        assert psi.mesh is m

    def test_from_ndarray_alias(self) -> None:
        m = _slab_mesh()
        arr = np.ones((m.quad.N, m.ng, *m.spatial_shape))
        psi = AngularFlux.from_ndarray(arr, m)
        assert isinstance(psi, AngularFlux)

    def test_zeros_on(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        np.testing.assert_array_equal(psi.values, 0.0)

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = _slab_mesh()
        with pytest.raises(ValueError, match="AngularFlux.*does not match"):
            AngularFlux.from_mesh(
                np.zeros((m.quad.N, m.ng, m.nx + 1)), m,
            )

    def test_2d_construction(self) -> None:
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        N = m.quad.N
        psi = AngularFlux.from_mesh(np.zeros((N, 2, 3, 2)), m)
        assert psi.values.shape == (N, 2, 3, 2)


# ───────────────────────────────────────────────────────────────────────
# E. Frozen contract
# ───────────────────────────────────────────────────────────────────────


class TestFrozen:
    def test_assign_values_raises(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        with pytest.raises(FrozenInstanceError):
            psi.values = np.zeros(psi.values.shape)  # type: ignore[misc]

    def test_assign_mesh_raises(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        with pytest.raises(FrozenInstanceError):
            psi.mesh = m  # type: ignore[misc]


# ───────────────────────────────────────────────────────────────────────
# F. Angular integration (reduction to ScalarFlux)
# ───────────────────────────────────────────────────────────────────────


class TestIntegrateAngular:
    def test_returns_scalar_flux(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        phi = psi.integrate_angular()
        assert isinstance(phi, ScalarFlux)
        assert phi.values.shape == (m.ng, *m.spatial_shape)

    def test_uniform_flux_gives_sum_of_weights(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.from_mesh(
            np.ones((m.quad.N, m.ng, *m.spatial_shape)), m,
        )
        phi = psi.integrate_angular()
        # φ = Σ_n w_n ψ_n; uniform ψ_n = 1 → φ = Σ_n w_n
        expected = m.quad.weights.sum()
        np.testing.assert_allclose(phi.values, expected, rtol=1e-15)

    def test_isotropic_flux_recovers_scalar(self) -> None:
        """For isotropic flux ψ_n(r,g) = φ(r,g) / Σw, integrate gives φ."""
        m = _slab_mesh()
        sum_w = m.quad.weights.sum()
        phi_target = 5.0
        psi_values = np.full((m.quad.N, m.ng, *m.spatial_shape), phi_target / sum_w)
        psi = AngularFlux.from_mesh(psi_values, m)
        phi = psi.integrate_angular()
        np.testing.assert_allclose(phi.values, phi_target, rtol=1e-15)


# ───────────────────────────────────────────────────────────────────────
# G. Metadata read-throughs
# ───────────────────────────────────────────────────────────────────────


class TestMetadata:
    def test_N_property(self) -> None:
        m = _slab_mesh()
        psi = AngularFlux.zeros_on(m)
        assert psi.N == m.quad.N

    def test_ng_property(self) -> None:
        m = _slab_mesh(ng=3)
        psi = AngularFlux.zeros_on(m)
        assert psi.ng == 3

    def test_spatial_shape_via_mesh(self) -> None:
        """C5.2 (#225): spatial reads are rank-generic via the mesh.

        The ``nx``/``ny`` field read-throughs are RETIRED — an
        ``(nx, ny)``-keyed read silently truncates a 3-D tensor.
        """
        m = _cartesian_2d_mesh(nx=3, ny=5)
        psi = AngularFlux.zeros_on(m)
        np.testing.assert_equal(psi.mesh.spatial_shape, (3, 5))
        with pytest.raises(AttributeError):
            psi.nx
        with pytest.raises(AttributeError):
            psi.ny
