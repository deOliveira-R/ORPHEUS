r"""L0/L2 — pure-Field :class:`AngularFlux` foundation tests.

Pins the post-D-H.1 contract for
:class:`orpheus.transport.fields.angular_flux.AngularFlux`:

* Inherits :class:`~orpheus.numerics.field.Field`; no hand-coded
  dunders.
* Pure Field — NO ``boundary`` attribute, NO ``_history`` attribute.
* Mesh-binding rejection (cross-mesh arithmetic raises).
* space-primary construction on the carrier's cached mint (the
  mesh-keyed sugar retired at CS4b S5).
* Frozen contract.
* ``integrate_angular`` reduction to ScalarFlux.

Symmetric with :class:`~orpheus.transport.fields.scalar_flux.ScalarFlux`
algebra (no special-case for boundary or history).

References
----------

* Depth B plan §3.3 (pure AngularFlux at L2).
* Pre-D-H.1 ``orpheus/sn/angular_flux.py:AngularFlux`` (legacy with
  ``boundary`` + ``_history`` fields — retired in D-H.2).

Layer DAG Phase 4 (2026-11-01): the file builds material meshes and spells each space it reads from the raw
constructors of orpheus.transport.phase_spaces (the SN problem's and the diffusion mesh's cached mints are those
constructors on the problem's own ingredients, gated in
tests/gates/sn/mesh/test_the_mints_are_the_raw_constructors.py and
tests/gates/diffusion/test_scalar_trace_is_the_raw_constructor.py).
"""
from __future__ import annotations

from dataclasses import FrozenInstanceError

import numpy as np
import pytest

from orpheus.geometry import StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh2D, Mesher
from orpheus.numerics.field import Field
from orpheus.transport.mesh.material_mesh import MaterialMesh
from orpheus.transport.phase_spaces import angular_bulk_space
from orpheus.transport.fields.angular_flux import AngularFlux
from orpheus.transport.fields.scalar_flux import ScalarFlux
from orpheus.numerics.axis_marginal import AxisRetractionOperator

from tests.gates.transport._materials import placeholder_materials
from orpheus.geometry.deck import MirrorDeck
from tests.gates.transport._phase_space_cases import posing_quadrature
from tests.gates.transport._phase_space_cases import slab_mesh


pytestmark = [pytest.mark.foundation]


# ───────────────────────────────────────────────────────────────────────
# Mesh fixtures
# ───────────────────────────────────────────────────────────────────────


def _stretched_mesh(nx: int = 4, ng: int = 2) -> MaterialMesh:
    """Same shape as ``_slab_mesh``, doubled width — the cell VOLUMES differ,
    so the carrier mints an UNEQUAL space (the F2 content discriminator)."""
    mesh = Mesher(StructuredGeometry.from_homogeneous(2.0, 0)).partition(
        CellsByCount.uniform_width(nx),
    ).mesh
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


def _cartesian_2d_mesh(nx: int = 3, ny: int = 2, ng: int = 2) -> MaterialMesh:
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, nx + 1),
        edges_y=np.linspace(0, 1, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        face_roles={"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": MirrorDeck(), "ymax": MirrorDeck()},
    )
    return MaterialMesh(mesh, placeholder_materials(ng=ng))


# ───────────────────────────────────────────────────────────────────────
# A. Field-ABC inheritance
# ───────────────────────────────────────────────────────────────────────


class TestFieldInheritance:
    def test_inherits_field(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        assert isinstance(psi, Field)

    def test_pure_field_no_boundary_attribute(self) -> None:
        """Critical: post-D-H.1 AngularFlux has NO .boundary field."""
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        assert not hasattr(psi, "boundary")

    def test_pure_field_no_history_attribute(self) -> None:
        """Critical: post-D-H.1 AngularFlux has NO ._history field."""
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        assert not hasattr(psi, "_history")

    def test_pure_field_no_history_depth_attribute(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
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
        m = slab_mesh()
        a = AngularFlux(values=np.ones((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape)), space=angular_bulk_space(m, posing_quadrature(m)))
        b = AngularFlux(values=2.0 * np.ones((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape)), space=angular_bulk_space(m, posing_quadrature(m)))
        s = a + b
        if type(s) is not AngularFlux:
            raise AssertionError("flux + flux left the leaf type")
        np.testing.assert_array_equal(s.values, 3.0)
        out = a + (b - a)  # the update step, plain V arithmetic
        if type(out) is not AngularFlux:
            raise AssertionError("the update step left the leaf type")
        np.testing.assert_array_almost_equal_nulp(out.values, b.values, nulp=4)

    def test_sub(self) -> None:
        m = slab_mesh()
        a = AngularFlux(values=np.full((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape), 5.0), space=angular_bulk_space(m, posing_quadrature(m)))
        b = AngularFlux(values=np.full((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape), 2.0), space=angular_bulk_space(m, posing_quadrature(m)))
        c = a - b
        np.testing.assert_array_equal(c.values, 3.0)

    def test_scalar_mul(self) -> None:
        m = slab_mesh()
        a = AngularFlux(values=np.full((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape), 2.0), space=angular_bulk_space(m, posing_quadrature(m)))
        c = a * 3.0
        np.testing.assert_array_equal(c.values, 6.0)

    def test_scalar_div(self) -> None:
        m = slab_mesh()
        a = AngularFlux(values=np.full((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape), 8.0), space=angular_bulk_space(m, posing_quadrature(m)))
        c = a / 4.0
        np.testing.assert_array_equal(c.values, 2.0)

    def test_neg(self) -> None:
        m = slab_mesh()
        a = AngularFlux(values=np.full((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape), 1.5), space=angular_bulk_space(m, posing_quadrature(m)))
        c = -a
        np.testing.assert_array_equal(c.values, -1.5)


# ───────────────────────────────────────────────────────────────────────
# C. Mesh-binding rejection
# ───────────────────────────────────────────────────────────────────────


class TestMeshBinding:
    def test_cross_carrier_discrimination_is_space_content(self) -> None:
        """CS4b S3 (F2): twin carriers mint EQUAL spaces and mix; a carrier
        with different cell volumes mints an UNEQUAL space and refuses."""
        m1 = slab_mesh()
        m2 = slab_mesh()  # different instance, same content
        a = AngularFlux.zeros(angular_bulk_space(m1, posing_quadrature(m1)))
        b = AngularFlux.zeros(angular_bulk_space(m2, posing_quadrature(m2)))
        _ = a - b  # twin content — legal since the F2 re-key
        _ = a + b
        c = AngularFlux.zeros(angular_bulk_space(_stretched_mesh(), posing_quadrature(_stretched_mesh())))
        with pytest.raises(ValueError, match="equal space"):
            a - c
        with pytest.raises(ValueError, match="equal space"):
            a + c

    def test_cross_class_rejected(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        phi = ScalarFlux(values=np.zeros((m.ng, *m.spatial_shape)), space=m.bulk_space)
        with pytest.raises(TypeError, match="same-class"):
            psi + phi  # type: ignore[operator]


# ───────────────────────────────────────────────────────────────────────
# D. Construction factories
# ───────────────────────────────────────────────────────────────────────


class TestConstruction:
    def test_ctor_on_the_carrier_mint(self) -> None:
        m = slab_mesh()
        space = angular_bulk_space(m, posing_quadrature(m))
        arr = np.ones((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape))
        psi = AngularFlux(values=arr, space=space)
        assert psi.values.shape == (posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape)
        assert psi.space is space

    def test_zeros_on_the_carrier_mint(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        np.testing.assert_array_equal(psi.values, 0.0)

    def test_shape_validation_rejects_wrong_shape(self) -> None:
        m = slab_mesh()
        with pytest.raises(ValueError, match="AngularFlux.*does not match"):
            AngularFlux(values=np.zeros((posing_quadrature(m).n_ordinates, m.ng, m.nx + 1)), space=angular_bulk_space(m, posing_quadrature(m)))

    def test_2d_construction(self) -> None:
        m = _cartesian_2d_mesh(nx=3, ny=2, ng=2)
        N = posing_quadrature(m).n_ordinates
        psi = AngularFlux(values=np.zeros((N, 2, 3, 2)), space=angular_bulk_space(m, posing_quadrature(m)))
        assert psi.values.shape == (N, 2, 3, 2)


# ───────────────────────────────────────────────────────────────────────
# E. Frozen contract
# ───────────────────────────────────────────────────────────────────────


class TestFrozen:
    def test_assign_values_raises(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        with pytest.raises(FrozenInstanceError):
            psi.values = np.zeros(psi.values.shape)  # type: ignore[misc]

    def test_assign_mesh_raises(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        with pytest.raises(FrozenInstanceError):
            psi.mesh = m  # type: ignore[misc]


# ───────────────────────────────────────────────────────────────────────
# F. Angular integration (reduction to ScalarFlux)
# ───────────────────────────────────────────────────────────────────────


class TestIntegrateAngular:
    def test_returns_scalar_flux(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        phi = psi.integrate_angular()
        assert isinstance(phi, ScalarFlux)
        assert phi.values.shape == (m.ng, *m.spatial_shape)

    def test_uniform_flux_gives_sum_of_weights(self) -> None:
        m = slab_mesh()
        psi = AngularFlux(values=np.ones((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape)), space=angular_bulk_space(m, posing_quadrature(m)))
        phi = psi.integrate_angular()
        # φ = Σ_n w_n ψ_n; uniform ψ_n = 1 → φ = Σ_n w_n
        expected = posing_quadrature(m).weights.sum()
        np.testing.assert_allclose(phi.values, expected, rtol=1e-15)

    def test_isotropic_flux_recovers_scalar(self) -> None:
        """For isotropic flux ψ_n(r,g) = φ(r,g) / Σw, integrate gives φ."""
        m = slab_mesh()
        sum_w = posing_quadrature(m).weights.sum()
        phi_target = 5.0
        psi_values = np.full((posing_quadrature(m).n_ordinates, m.ng, *m.spatial_shape), phi_target / sum_w)
        psi = AngularFlux(values=psi_values, space=angular_bulk_space(m, posing_quadrature(m)))
        phi = psi.integrate_angular()
        np.testing.assert_allclose(phi.values, phi_target, rtol=1e-15)


# ───────────────────────────────────────────────────────────────────────
# G. Metadata read-throughs
# ───────────────────────────────────────────────────────────────────────


class TestMetadata:
    def test_N_property(self) -> None:
        m = slab_mesh()
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        assert psi.N == posing_quadrature(m).n_ordinates

    def test_ng_property(self) -> None:
        m = slab_mesh(ng=3)
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        assert psi.ng == 3

    def test_spatial_shape_via_mesh(self) -> None:
        """C5.2 (#225): spatial reads are rank-generic via the mesh.

        The ``nx``/``ny`` field read-throughs are RETIRED — an
        ``(nx, ny)``-keyed read silently truncates a 3-D tensor.
        """
        m = _cartesian_2d_mesh(nx=3, ny=5)
        psi = AngularFlux.zeros(angular_bulk_space(m, posing_quadrature(m)))
        np.testing.assert_equal(psi.values.shape[2:], (3, 5))
        with pytest.raises(AttributeError):
            psi.nx
        with pytest.raises(AttributeError):
            psi.ny


class TestIntegrateAngularRidesTheAxisMarginal:
    """CS4b S6.2: both halves of the reduction ride the angular axis's ONE
    memoized marginal (:mod:`orpheus.numerics.axis_marginal`)."""

    def test_space_is_the_memoized_marginal(self) -> None:
        m = slab_mesh()
        space = angular_bulk_space(m, posing_quadrature(m))
        psi = AngularFlux.zeros(space)
        phi = psi.integrate_angular()
        R = AxisRetractionOperator.of(space, "angular")
        assert phi.space is R.codomain
        assert phi.space == m.bulk_space
