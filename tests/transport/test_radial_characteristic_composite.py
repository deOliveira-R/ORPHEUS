r"""B.1d — RadialCharacteristicComposite (System B as an independent composite).

Intrinsic-property + split-fidelity gates for the Phase-B (coupled-block campaign)
ψ½ composite
``RadialCharacteristicComposite = Composite[RadialCharacteristicInteriorFlux,
RadialCharacteristicBoundaryFlux]`` — System B's own ``interior ⊕ boundary`` field,
parallel to System A's :class:`~orpheus.transport.full_field.FullField`. It is a
TRIVIAL subclass of the relaxed :class:`~orpheus.transport.full_field.Composite`
base (B.1a — no hook overrides), so this file pins BOTH:

* **The intrinsic composite algebra on the RAY carrier** — the multi-instantiation
  crux (test-architect ``coupled_operator_b1_split_verification`` G-C1): System B
  is the THIRD ``Composite`` instantiation and the FIRST whose leaves are neither
  ``BulkField`` nor ``BoundaryField``. Running the inherited generic algebra on
  these leaves would RED if the generic body hardcoded an ``AngularFlux`` /
  ``BulkField`` assumption, while System A (``AngularFlux``) + the scalar
  ``test_composite`` instantiation stay green — the asymmetry is the evidence.
* **The split-fidelity bridge** (:meth:`from_unified` / :meth:`to_unified`) — the
  "split IS the unified, re-typed" contract: the additive-migration seam +
  retirement licence for the unified ψ½ leaf.

vv Mode-8 discipline: ``np.testing.assert_*`` / ``pytest.raises`` only.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
import orpheus.transport.displacements  # noqa: F401  eager _BY_REP registry
from orpheus.transport.displacements.radial_characteristic_boundary_displacement import (
    RadialCharacteristicBoundaryDisplacement,
)
from orpheus.transport.displacements.radial_characteristic_interior_displacement import (
    RadialCharacteristicInteriorDisplacement,
)
from orpheus.transport.fields.radial_characteristic_boundary_flux import (
    RadialCharacteristicBoundaryFlux,
)
from orpheus.transport.fields.radial_characteristic_flux import RadialCharacteristicFlux
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.full_field import Composite, FullField
from orpheus.transport.radial_characteristic_composite import (
    RadialCharacteristicComposite,
)
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

_SIGNS = (-1, +1)


def _mesh(coord: CoordSystem, **bc) -> SNMesh:
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 5), mat_ids=np.zeros(4, dtype=int),
        coord=coord, **bc,
    )
    return SNMesh(mesh, Quadrature.gauss_legendre(4), placeholder_materials(ng=2))


def _sphere() -> SNMesh:
    return _mesh(CoordSystem.SPHERICAL, bc_right=BC("reflective"))


def _slab() -> SNMesh:
    return _mesh(
        CoordSystem.CARTESIAN, bc_left=BC("reflective"), bc_right=BC("reflective"),
    )


def _rand_unified(sn: SNMesh, seed: int) -> RadialCharacteristicFlux:
    n = RadialCharacteristicFlux.zeros_on(sn).values.size
    rng = np.random.default_rng(seed)
    return RadialCharacteristicFlux.from_mesh(rng.random(n) + 0.5, sn)


def _rand_composite(sn: SNMesh, seed: int) -> RadialCharacteristicComposite:
    return RadialCharacteristicComposite.from_unified(_rand_unified(sn, seed))


# ── Construction (the multi-instantiation structural crux) ──


class TestCompositeConstruction:
    def test_is_a_composite_but_not_a_fullfield(self) -> None:
        c = RadialCharacteristicComposite.from_mesh(_sphere())
        if not isinstance(c, Composite):
            pytest.fail("System B is not a Composite")
        if isinstance(c, FullField):
            pytest.fail("System B must NOT be a FullField (a distinct composite)")

    def test_has_no_radial_characteristic_third_slot(self) -> None:
        # A pure 2-block composite — the ψ½ third block is FullField's, not this.
        c = RadialCharacteristicComposite.from_mesh(_sphere())
        if hasattr(c, "radial_characteristic"):
            pytest.fail("System B is pure 2-block; it carries no ψ½ third slot")

    def test_leaf_types_are_the_ray_split_loci(self) -> None:
        # The 3rd Composite instantiation, FIRST with non-Bulk/Boundary leaves.
        c = RadialCharacteristicComposite.from_mesh(_sphere())
        if type(c.interior) is not RadialCharacteristicInteriorFlux:
            pytest.fail(f"interior is {type(c.interior).__name__}")
        if type(c.boundary) is not RadialCharacteristicBoundaryFlux:
            pytest.fail(f"boundary is {type(c.boundary).__name__}")

    def test_mesh_reads_off_the_leaves(self) -> None:
        sn = _sphere()
        c = RadialCharacteristicComposite.from_mesh(sn)
        if c.mesh is not sn:
            pytest.fail("composite mesh must be the leaves' mesh")


# ── The split-fidelity bridge (retirement licence) ──


class TestSplitFidelityBridge:
    def test_from_unified_maps_cells_to_interior_and_corner_to_boundary(self) -> None:
        sn = _sphere()
        u = _rand_unified(sn, 1)
        c = RadialCharacteristicComposite.from_unified(u)
        for level in u.levels:
            for sign in _SIGNS:
                np.testing.assert_array_equal(
                    c.interior.cells(level, sign), u.cells(level, sign),
                )
                np.testing.assert_array_equal(
                    c.boundary.corner(level, sign), u.corner(level, sign),
                )

    def test_to_unified_after_from_unified_is_exact(self) -> None:
        """The retirement licence: to_unified(from_unified(u)) == u, bitwise."""
        sn = _sphere()
        u = _rand_unified(sn, 2)
        back = RadialCharacteristicComposite.from_unified(u).to_unified()
        if type(back) is not RadialCharacteristicFlux:
            pytest.fail(f"to_unified returned {type(back).__name__}")
        np.testing.assert_array_equal(back.values, u.values)

    def test_from_unified_after_to_unified_is_exact(self) -> None:
        sn = _sphere()
        c = _rand_composite(sn, 3)
        rt = RadialCharacteristicComposite.from_unified(c.to_unified())
        np.testing.assert_array_equal(rt.interior.values, c.interior.values)
        np.testing.assert_array_equal(rt.boundary.values, c.boundary.values)


# ── Intrinsic composite algebra on the ray carrier (multi-instantiation) ──


class TestCompositeAlgebra:
    def test_flux_plus_flux_is_the_affine_gate(self) -> None:
        sn = _sphere()
        a, b = _rand_composite(sn, 4), _rand_composite(sn, 5)
        with pytest.raises(TypeError, match="affine"):
            _ = a + b

    def test_subtraction_mints_a_displacement_composite_per_block(self) -> None:
        sn = _sphere()
        a, b = _rand_composite(sn, 6), _rand_composite(sn, 7)
        d = a - b
        if type(d.interior) is not RadialCharacteristicInteriorDisplacement:
            pytest.fail(f"interior block minted {type(d.interior).__name__}")
        if type(d.boundary) is not RadialCharacteristicBoundaryDisplacement:
            pytest.fail(f"boundary block minted {type(d.boundary).__name__}")
        np.testing.assert_array_equal(
            d.interior.values, a.interior.values - b.interior.values,
        )

    def test_torsor_add_displacement_recovers_the_point(self) -> None:
        sn = _sphere()
        a, b = _rand_composite(sn, 8), _rand_composite(sn, 9)
        recovered = a + (b - a)  # flux ⊕ displacement → flux
        if type(recovered) is not RadialCharacteristicComposite:
            pytest.fail(f"torsor returned {type(recovered).__name__}")
        np.testing.assert_allclose(recovered.interior.values, b.interior.values)
        np.testing.assert_allclose(recovered.boundary.values, b.boundary.values)

    def test_scalar_scales_both_blocks(self) -> None:
        sn = _sphere()
        a = _rand_composite(sn, 10)
        out = 2.0 * a
        np.testing.assert_allclose(out.interior.values, 2.0 * a.interior.values)
        np.testing.assert_allclose(out.boundary.values, 2.0 * a.boundary.values)
        np.testing.assert_allclose((a / 2.0).boundary.values, 0.5 * a.boundary.values)

    def test_to_flat_from_flat_round_trip_preserves_type(self) -> None:
        sn = _sphere()
        a = _rand_composite(sn, 11)
        expected = np.concatenate([a.interior.values.ravel(), a.boundary.values])
        np.testing.assert_array_equal(a.to_flat(), expected)
        back = Composite.from_flat(a.to_flat(), a)
        if type(back) is not RadialCharacteristicComposite:
            pytest.fail(f"from_flat returned {type(back).__name__}")
        np.testing.assert_array_equal(back.interior.values, a.interior.values)
        np.testing.assert_array_equal(back.boundary.values, a.boundary.values)

    def test_from_flat_rejects_wrong_size(self) -> None:
        sn = _sphere()
        a = _rand_composite(sn, 12)
        with pytest.raises(ValueError, match="flat.size"):
            Composite.from_flat(a.to_flat()[:-1], a)

    def test_copy_is_value_independent(self) -> None:
        sn = _sphere()
        a = _rand_composite(sn, 13)
        c = a.copy()
        np.testing.assert_array_equal(c.interior.values, a.interior.values)
        c.interior.cells(0, -1)[0, 0] += 1.0
        if a.interior.values[0] == c.interior.values[0]:
            pytest.fail("copy aliases the original interior buffer")

    def test_rejects_non_composite_partner(self) -> None:
        a = _rand_composite(_sphere(), 14)
        with pytest.raises(TypeError, match="same-class partner"):
            _ = a + 42  # type: ignore[operator]


# ── Presence = ray-carrying ──


class TestPresence:
    def test_from_mesh_raises_on_a_noncarrying_mesh(self) -> None:
        with pytest.raises(ValueError, match="carries no radial_characteristic"):
            RadialCharacteristicComposite.from_mesh(_slab())
