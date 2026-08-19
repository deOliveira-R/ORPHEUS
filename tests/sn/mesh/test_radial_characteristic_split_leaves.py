r"""B.1c — the split ψ½ leaf types (System B's interior ⊕ boundary flux pair).

Intrinsic-property gates for the Phase-B (coupled-block campaign) split-locus leaf
types: the two field bases
:class:`~orpheus.transport.fields._bases.RadialCharacteristicInteriorField` /
:class:`~orpheus.transport.fields._bases.RadialCharacteristicBoundaryField` and
their flux role leaves. These are the DATA carriers System B's composite (B.1d)
pairs; here we pin their construction, views, presence, and the per-locus V
flux algebra ([[feedback-test-intrinsic-properties]]).

The LOAD-BEARING gate is **blockwise algebra over DISTINCT reps**: the split
flux leaves carry distinct Field bases, so arithmetic stays per-locus (interior
with interior, boundary with boundary — cross-locus refuses at Layer 1), never
a collision onto the unified ψ½ or the wrong locus. (Until campaign 1 CS3,
2026-08-19, this was spelled through the affine-torsor mint: ``flux ⊖ flux``
resolved via ``Displacement.sibling_of`` to per-locus displacement siblings;
the displacement family is retired — flux lives in V.)

vv Mode-8 discipline: ``np.testing.assert_*`` / ``pytest.raises`` only.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.units import ANGULAR_FLUX_UNITS, ANGULAR_RATE_UNITS
from orpheus.transport.fields.radial_characteristic_boundary_flux import (
    RadialCharacteristicBoundaryFlux,
)
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

_NG, _NX = 2, 4


def _mesh_1d(coord: CoordSystem, quad, *, nx: int = _NX, ng: int = _NG) -> SNMesh:
    edges = np.linspace(0.0, 1.0, nx + 1)
    mat_ids = np.zeros(nx, dtype=int)
    if coord is CoordSystem.CARTESIAN:
        mesh = Mesh1D(
            edges=edges, mat_ids=mat_ids, coord=coord,
            bc_left=BC("reflective"), bc_right=BC("reflective"),
        )
    else:
        mesh = Mesh1D(
            edges=edges, mat_ids=mat_ids, coord=coord, bc_right=BC("reflective"),
        )
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _sphere() -> SNMesh:
    return _mesh_1d(CoordSystem.SPHERICAL, Quadrature.gauss_legendre(4))


def _slab() -> SNMesh:
    # Since Q5.6.3 the slab is the only admitted non-carrying 1-D geometry
    # (the file's former `_cyl_product` non-carrying fixture became
    # unconstructible at the admission flip; its refusal is gated in
    # ``test_cylindrical_quadrature_admission.py``).
    return _mesh_1d(CoordSystem.CARTESIAN, Quadrature.gauss_legendre(4))


def _rand(cls, mesh, seed: int):
    """A random flux leaf on ``mesh`` (distinct nonzero values)."""
    space = cls._face_space_of(mesh)
    rng = np.random.default_rng(seed)
    return cls.from_mesh(rng.random(space.shape[0]) + 0.5, mesh)


# ── Construction + views + presence ──


class TestSplitLeafConstruction:
    def test_interior_flux_shape_space_and_cells_view(self) -> None:
        sn = _sphere()
        f = RadialCharacteristicInteriorFlux.zeros_on(sn)
        if f.space.name != "radial_characteristic_interior":
            pytest.fail(f"interior flux space is {f.space.name!r}")
        np.testing.assert_array_equal(f.cells(0, -1).shape, (_NG, _NX))
        np.testing.assert_array_equal(f.cells(0, +1).shape, (_NG, _NX))

    def test_boundary_flux_shape_space_and_corner_view(self) -> None:
        sn = _sphere()
        f = RadialCharacteristicBoundaryFlux.zeros_on(sn)
        if f.space.name != "radial_characteristic_boundary":
            pytest.fail(f"boundary flux space is {f.space.name!r}")
        np.testing.assert_array_equal(f.corner(0, -1).shape, (_NG,))

    def test_from_mesh_round_trips_values(self) -> None:
        sn = _sphere()
        space = RadialCharacteristicInteriorFlux._face_space_of(sn)
        vals = np.arange(space.shape[0], dtype=float)
        f = RadialCharacteristicInteriorFlux.from_mesh(vals, sn)
        np.testing.assert_array_equal(f.values, vals)

    def test_view_writes_back_through_the_flat_buffer(self) -> None:
        sn = _sphere()
        f = RadialCharacteristicInteriorFlux.zeros_on(sn)
        f.cells(0, -1)[...] = 7.0
        # the (0,-1) cells block is the first (ng*nx) of the flat buffer
        np.testing.assert_array_equal(f.values[: _NG * _NX], np.full(_NG * _NX, 7.0))

    @pytest.mark.parametrize("mesh_fn", [_slab])
    @pytest.mark.parametrize(
        "cls",
        [
            RadialCharacteristicInteriorFlux,
            RadialCharacteristicBoundaryFlux,
            RadialCharacteristicInteriorSourceSink,
            RadialCharacteristicBoundarySourceSink,
        ],
    )
    def test_presence_noncarrying_mesh_raises(self, cls, mesh_fn) -> None:
        sn = mesh_fn()
        with pytest.raises(ValueError, match="carries no radial_characteristic"):
            cls.zeros_on(sn)


# ── The V flux algebra (per locus, over DISTINCT reps — CS3 cone carve) ──


class TestSplitLeafVectorAlgebra:
    def test_interior_subtraction_stays_interior_typed(self) -> None:
        """(Until CS3 this minted RadialCharacteristicInteriorDisplacement;
        the per-locus DISTINCTNESS claim now lives on the flux pair below.)"""
        sn = _sphere()
        a = _rand(RadialCharacteristicInteriorFlux, sn, 1)
        b = _rand(RadialCharacteristicInteriorFlux, sn, 2)
        d = a - b
        if type(d) is not RadialCharacteristicInteriorFlux:
            pytest.fail(f"interior − became {type(d).__name__}")
        np.testing.assert_array_equal(d.values, a.values - b.values)

    def test_boundary_subtraction_stays_boundary_typed(self) -> None:
        sn = _sphere()
        a = _rand(RadialCharacteristicBoundaryFlux, sn, 3)
        b = _rand(RadialCharacteristicBoundaryFlux, sn, 4)
        d = a - b
        if type(d) is not RadialCharacteristicBoundaryFlux:
            pytest.fail(f"boundary − became {type(d).__name__}")

    def test_the_two_locus_leaves_are_distinct_types(self) -> None:
        # Distinct Field bases ⟹ the two loci never blur (blockwise algebra).
        if RadialCharacteristicInteriorFlux is RadialCharacteristicBoundaryFlux:
            pytest.fail("interior/boundary leaves collapsed to one type")
        with pytest.raises(TypeError):  # cross-locus arithmetic refuses
            _ = _rand(RadialCharacteristicInteriorFlux, _sphere(), 21) + _rand(  # type: ignore[operator]
                RadialCharacteristicBoundaryFlux, _sphere(), 22,
            )

    @pytest.mark.parametrize(
        "cls", [RadialCharacteristicInteriorFlux, RadialCharacteristicBoundaryFlux],
    )
    def test_update_recovery_a_plus_b_minus_a_is_b(self, cls) -> None:
        sn = _sphere()
        a, b = _rand(cls, sn, 5), _rand(cls, sn, 6)
        recovered = a + (b - a)  # plain V arithmetic
        if type(recovered) is not cls:
            pytest.fail(f"torsor action returned {type(recovered).__name__}")
        np.testing.assert_allclose(recovered.values, b.values)

    @pytest.mark.parametrize(
        "cls", [RadialCharacteristicInteriorFlux, RadialCharacteristicBoundaryFlux],
    )
    def test_flux_plus_flux_is_legal_same_typed(self, cls) -> None:
        """The CS3 inversion of the retired affine-gate row."""
        sn = _sphere()
        a, b = _rand(cls, sn, 7), _rand(cls, sn, 8)
        s = a + b
        if type(s) is not cls:
            pytest.fail(f"flux + flux returned {type(s).__name__}")
        np.testing.assert_array_equal(s.values, a.values + b.values)

    @pytest.mark.parametrize(
        "cls", [RadialCharacteristicInteriorFlux, RadialCharacteristicBoundaryFlux],
    )
    def test_scalar_scaling_stays_a_flux(self, cls) -> None:
        sn = _sphere()
        a = _rand(cls, sn, 9)
        out = 2.0 * a
        if type(out) is not cls:
            pytest.fail(f"scalar * flux returned {type(out).__name__}")
        np.testing.assert_allclose(out.values, 2.0 * a.values)
        np.testing.assert_allclose((a / 2.0).values, 0.5 * a.values)

    def test_cross_locus_subtraction_is_forbidden(self) -> None:
        sn = _sphere()
        i = _rand(RadialCharacteristicInteriorFlux, sn, 10)
        b = _rand(RadialCharacteristicBoundaryFlux, sn, 11)
        with pytest.raises(TypeError, match="same-class"):
            _ = i - b  # type: ignore[operator]

    def test_cross_mesh_arithmetic_is_forbidden(self) -> None:
        a = _rand(RadialCharacteristicInteriorFlux, _sphere(), 12)
        b = _rand(RadialCharacteristicInteriorFlux, _sphere(), 12)  # a DIFFERENT mesh
        with pytest.raises(ValueError, match="mesh"):
            _ = a - b


# ── The split SourceSink leaves (B.2b b1 — plain vector role per locus) ──


class TestSplitSourceSinkLeaves:
    r"""The q½ source split: plain vector-space leaves over the locus bases.

    Unlike the flux leaves (affine points, ``⊖`` mints displacements), source
    sums are CLOSED — no role mixin. The Field class-identity gate is what
    keeps source, flux, and displacement arithmetic from mixing (all roles per
    locus share the SAME split space, so the class is the gate). The units row
    pins the B.2b ruling that the split DISSOLVES the unified leaf's
    documented corner-units deviation into two honest per-locus identities.
    """

    @pytest.mark.parametrize(
        "cls",
        [RadialCharacteristicInteriorSourceSink, RadialCharacteristicBoundarySourceSink],
    )
    def test_source_plus_source_stays_source(self, cls) -> None:
        sn = _sphere()
        a, b = _rand(cls, sn, 20), _rand(cls, sn, 21)
        out = a + b
        if type(out) is not cls:
            pytest.fail(f"source + source returned {type(out).__name__}")
        np.testing.assert_array_equal(out.values, a.values + b.values)

    def test_cross_role_sum_is_forbidden(self) -> None:
        # Source and flux share the SAME interior space — the CLASS identity
        # (not the space) is the arithmetic gate rejecting the cross-role sum.
        sn = _sphere()
        src = _rand(RadialCharacteristicInteriorSourceSink, sn, 22)
        flx = _rand(RadialCharacteristicInteriorFlux, sn, 23)
        with pytest.raises(TypeError):
            _ = src + flx  # type: ignore[operator]

    def test_units_dissolve_the_unified_corner_deviation(self) -> None:
        r"""Interior cells = volumetric rate (like ``AngularSourceSink``);
        boundary corner = trace-like flux datum (like
        ``AngularBoundarySourceSink``) — two honest per-locus identities where
        the unified leaf carried one identity plus a documented deviation."""
        if RadialCharacteristicInteriorSourceSink.UNITS is not ANGULAR_RATE_UNITS:
            pytest.fail("interior q½ cells must carry ANGULAR_RATE_UNITS")
        if RadialCharacteristicBoundarySourceSink.UNITS is not ANGULAR_FLUX_UNITS:
            pytest.fail("boundary q½ corner must carry ANGULAR_FLUX_UNITS "
                        "(the trace convention — no volumetric cm⁻³)")
