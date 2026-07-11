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

from dataclasses import dataclass

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
import orpheus.transport.displacements  # noqa: F401  eager _BY_REP registry
from orpheus.transport.displacements.radial_characteristic_boundary_displacement import (
    RadialCharacteristicBoundaryDisplacement,
)
from orpheus.transport.displacements.radial_characteristic_displacement import (
    RadialCharacteristicDisplacement,
)
from orpheus.transport.displacements.radial_characteristic_interior_displacement import (
    RadialCharacteristicInteriorDisplacement,
)
from orpheus.transport.fields._bases import (
    RadialCharacteristicField,
    RadialCharacteristicInteriorField,
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
from orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)
from orpheus.transport.source_sinks.radial_characteristic_source_sink import (
    RadialCharacteristicSourceSink,
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


def _rand_unified_of(cls, sn: SNMesh, seed: int):
    """A random unified ψ½ leaf of the given ROLE class (bridge role rows)."""
    n = cls.zeros_on(sn).values.size
    rng = np.random.default_rng(seed)
    return cls.from_mesh(rng.random(n) + 0.5, sn)


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

    def test_source_members_construct_role_erased(self) -> None:
        r"""B.2b DP2: the slots bind the FIELD BASES (the FullField precedent),
        so an operator emission — the (interior, boundary) SOURCE pair that
        A_BA / B_b write — rides the SAME composite class as the flux state."""
        sn = _sphere()
        comp = RadialCharacteristicComposite(
            interior=RadialCharacteristicInteriorSourceSink.zeros_on(sn),
            boundary=RadialCharacteristicBoundarySourceSink.zeros_on(sn),
        )
        if type(comp.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"interior member is {type(comp.interior).__name__}")
        if type(comp.boundary) is not RadialCharacteristicBoundarySourceSink:
            pytest.fail(f"boundary member is {type(comp.boundary).__name__}")

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

    # ── Role preservation (B.2b DP2 — one bridge body, exact-class table) ──

    @pytest.mark.parametrize(
        "unified_cls,interior_cls,boundary_cls",
        [
            (RadialCharacteristicSourceSink,
             RadialCharacteristicInteriorSourceSink,
             RadialCharacteristicBoundarySourceSink),
            (RadialCharacteristicDisplacement,
             RadialCharacteristicInteriorDisplacement,
             RadialCharacteristicBoundaryDisplacement),
        ],
    )
    def test_bridge_preserves_the_member_role(
        self, unified_cls, interior_cls, boundary_cls,
    ) -> None:
        r"""from_unified follows the unified leaf's ROLE: a source leaf splits
        into the source pair (the A_BA / B_b emission carrier), a displacement
        into the displacement pair — values mapped bitwise per leg."""
        sn = _sphere()
        u = _rand_unified_of(unified_cls, sn, 31)
        comp = RadialCharacteristicComposite.from_unified(u)
        if type(comp.interior) is not interior_cls:
            pytest.fail(f"interior member is {type(comp.interior).__name__}, "
                        f"expected {interior_cls.__name__}")
        if type(comp.boundary) is not boundary_cls:
            pytest.fail(f"boundary member is {type(comp.boundary).__name__}, "
                        f"expected {boundary_cls.__name__}")
        for level in u.levels:
            for sign in _SIGNS:
                np.testing.assert_array_equal(
                    comp.interior.cells(level, sign), u.cells(level, sign),
                )
                np.testing.assert_array_equal(
                    comp.boundary.corner(level, sign), u.corner(level, sign),
                )

    @pytest.mark.parametrize(
        "unified_cls",
        [RadialCharacteristicSourceSink, RadialCharacteristicDisplacement],
    )
    def test_round_trip_is_exact_per_role(self, unified_cls) -> None:
        r"""The retirement licence holds PER ROLE: to∘from == id bitwise, and
        the round-tripped leaf carries the original role class."""
        sn = _sphere()
        u = _rand_unified_of(unified_cls, sn, 32)
        back = RadialCharacteristicComposite.from_unified(u).to_unified()
        if type(back) is not unified_cls:
            pytest.fail(f"round trip returned {type(back).__name__}, "
                        f"expected {unified_cls.__name__}")
        np.testing.assert_array_equal(back.values, u.values)

    def test_from_unified_rejects_an_unknown_role(self) -> None:
        r"""An off-table unified role refuses loudly (a new role joins by a
        table row, never a silent flux-pair default)."""
        @dataclass(frozen=True, eq=False, kw_only=True, repr=False)
        class _AlienUnified(RadialCharacteristicField):
            pass

        sn = _sphere()
        with pytest.raises(TypeError, match="role-preserving"):
            RadialCharacteristicComposite.from_unified(_AlienUnified.zeros_on(sn))

    def test_to_unified_rejects_a_mixed_role_composite(self) -> None:
        r"""Role-erased slots make a mixed-role composite CONSTRUCTABLE, but it
        has no unified counterpart — the bridge is where role identity is
        parsed (the consumption-site F2 discipline)."""
        sn = _sphere()
        comp = RadialCharacteristicComposite(
            interior=RadialCharacteristicInteriorFlux.zeros_on(sn),
            boundary=RadialCharacteristicBoundarySourceSink.zeros_on(sn),
        )
        with pytest.raises(TypeError, match="mixed-role"):
            comp.to_unified()

    def test_to_unified_rejects_an_unknown_interior_role(self) -> None:
        @dataclass(frozen=True, eq=False, kw_only=True, repr=False)
        class _AlienInterior(RadialCharacteristicInteriorField):
            pass

        sn = _sphere()
        comp = RadialCharacteristicComposite(
            interior=_AlienInterior.zeros_on(sn),
            boundary=RadialCharacteristicBoundaryFlux.zeros_on(sn),
        )
        with pytest.raises(TypeError, match="no unified role"):
            comp.to_unified()


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


# ── System B's member space (B.2b b2 — the family-blind FullFieldSpace reuse) ──


class TestRadialCharacteristicCompositeSpace:
    r"""``SNMesh.radial_characteristic_composite_space`` — System B's member space.

    The DP1 ruling realized: the SAME family-blind ``FullFieldSpace`` class
    System A uses, instantiated over the two split ψ½ spaces under the
    identity name ``"radial_characteristic"``. The metric rows exercise the
    ``_rebuild`` presence-dispatch 2-BLOCK arm on the REAL composite (before
    B.2b this path raised — ``_rebuild`` always passed the seed kwarg to a
    hook that has none); the FullField row pins the 3-slot arm byte-identical
    (test-architect G-b2: exercise BOTH arms). Space ``==`` is
    ``(name, shape)``-only (Mode-12 note), so member identity is asserted
    with ``is``, and every value row is ``array_equal`` (bit-id — a re-label
    moves no arithmetic).
    """

    def test_identity_name_and_shape(self) -> None:
        sn = _sphere()
        space = sn.radial_characteristic_composite_space
        if space is None:
            pytest.fail("sphere (carrying) must mint the composite space")
        n_i = sn.radial_characteristic_interior_space.shape[0]
        n_b = sn.radial_characteristic_boundary_space.shape[0]
        if space.name != "radial_characteristic":
            pytest.fail(f"space name is {space.name!r}")
        np.testing.assert_array_equal(space.shape, (n_i + n_b,))

    def test_member_spaces_are_the_split_instances(self) -> None:
        # Space == is (name, shape)-only — assert the MEMBER identity with
        # ``is`` (the cached split spaces, one instance per mesh).
        sn = _sphere()
        space = sn.radial_characteristic_composite_space
        if space.interior_space is not sn.radial_characteristic_interior_space:
            pytest.fail("interior member is not THE cached interior split space")
        if space.trace_space is not sn.radial_characteristic_boundary_space:
            pytest.fail("trace member is not THE cached boundary split space")
        if space.radial_characteristic_space is not None:
            pytest.fail("System B's own space must carry NO third slot")

    def test_presence_and_caching(self) -> None:
        if _slab().radial_characteristic_composite_space is not None:
            pytest.fail("a non-carrying mesh has no System B, hence no space")
        sn = _sphere()
        if (sn.radial_characteristic_composite_space
                is not sn.radial_characteristic_composite_space):
            pytest.fail("the composite space must be cached (one per mesh)")

    def test_metric_trio_dispatches_member_wise_on_the_real_composite(self) -> None:
        r"""apply_metric / apply_inverse_metric / inner_product on the REAL
        2-block composite ≡ direct split-space application per member,
        bitwise — the presence-dispatch 2-block arm, type-preserving. This
        row is the TOOTH for a kwarg-inversion mutation (always passing the
        seed kwarg REDs it with a TypeError)."""
        sn = _sphere()
        space = sn.radial_characteristic_composite_space
        x = _rand_composite(sn, 51)
        y = _rand_composite(sn, 52)

        gx = space.apply_metric(x)
        if type(gx) is not RadialCharacteristicComposite:
            pytest.fail(f"apply_metric returned {type(gx).__name__}")
        np.testing.assert_array_equal(
            gx.interior.values,
            sn.radial_characteristic_interior_space.apply_metric(x.interior.values),
        )
        np.testing.assert_array_equal(
            gx.boundary.values,
            sn.radial_characteristic_boundary_space.apply_metric(x.boundary.values),
        )

        ginv = space.apply_inverse_metric(x)
        np.testing.assert_array_equal(
            ginv.interior.values,
            sn.radial_characteristic_interior_space.apply_inverse_metric(
                x.interior.values,
            ),
        )
        np.testing.assert_array_equal(
            ginv.boundary.values,
            sn.radial_characteristic_boundary_space.apply_inverse_metric(
                x.boundary.values,
            ),
        )

        # inner_product is a SCALAR functional (Mode-12-blind to per-member
        # swaps on its own) — paired with the array rows above per the memo.
        np.testing.assert_array_equal(
            space.inner_product(x, y),
            sn.radial_characteristic_interior_space.inner_product(
                x.interior.values, y.interior.values,
            )
            + sn.radial_characteristic_boundary_space.inner_product(
                x.boundary.values, y.boundary.values,
            ),
        )

    def test_metric_preserves_the_member_role(self) -> None:
        r"""The rebuilt members keep their leaf classes (replace-based
        rebuild) — a source composite through G stays a source composite."""
        sn = _sphere()
        space = sn.radial_characteristic_composite_space
        src = RadialCharacteristicComposite(
            interior=RadialCharacteristicInteriorSourceSink.zeros_on(sn),
            boundary=RadialCharacteristicBoundarySourceSink.zeros_on(sn),
        )
        out = space.apply_metric(src)
        if type(out.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"interior became {type(out.interior).__name__}")
        if type(out.boundary) is not RadialCharacteristicBoundarySourceSink:
            pytest.fail(f"boundary became {type(out.boundary).__name__}")

    def test_full_field_three_slot_arm_is_byte_identical(self) -> None:
        r"""The OTHER presence-dispatch arm: a seed-carrying FullField still
        routes its ψ½ block through the seed space (byte-identity vs the
        direct leaf-space application — the pre-B.2b behavior)."""
        from orpheus.transport.fields.angular_boundary_flux import (
            AngularBoundaryFlux,
        )
        from orpheus.transport.fields.angular_flux import AngularFlux

        sn = _sphere()
        ff = FullField.zeros(
            interior=AngularFlux, boundary=AngularBoundaryFlux, mesh=sn,
            radial_characteristic=RadialCharacteristicFlux,
        )
        rng = np.random.default_rng(53)
        ff.interior.values[...] = rng.random(ff.interior.values.shape)
        ff.boundary.values[...] = rng.random(ff.boundary.values.shape)
        if ff.radial_characteristic is None:
            pytest.fail("sphere FullField must carry the ψ½ block (R12a)")
        ff.radial_characteristic.values[...] = rng.random(
            ff.radial_characteristic.values.shape,
        )
        out = sn.full_field_space.apply_metric(ff)
        np.testing.assert_array_equal(
            out.radial_characteristic.values,
            sn.radial_characteristic_space.apply_metric(
                ff.radial_characteristic.values,
            ),
        )

    def test_flat_layout_is_coextensive_with_the_space_shape(self) -> None:
        r"""The elegance-carried B.2a reminder, pinned on the REAL members:
        ``to_flat().size == prod(space.shape)`` AND the flat ORDER is
        ``concat(interior, boundary)`` — the three offset spellings
        (``to_flat``, the space shape, the future ``system_slices``) name
        ONE layout."""
        sn = _sphere()
        space = sn.radial_characteristic_composite_space
        c = _rand_composite(sn, 54)
        flat = c.to_flat()
        np.testing.assert_array_equal(flat.size, int(np.prod(space.shape)))
        np.testing.assert_array_equal(
            flat,
            np.concatenate(
                [c.interior.values.ravel(), c.boundary.values.ravel()],
            ),
        )
