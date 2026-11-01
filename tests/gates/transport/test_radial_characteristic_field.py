r"""B.1d — RadialCharacteristicField (System B as an independent composite).

Intrinsic-property + split-fidelity gates for the Phase-B (coupled-block campaign)
ψ½ composite
``RadialCharacteristicField = Composite[RadialCharacteristicInteriorFlux,
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

Layer DAG Phase 4 (2026-11-01): the composite's spaces are posed by hand from the numerics split-space classes, and
the rows whose subject is the SN problem's mint (content equals the mint, members are its cached split spaces,
presence and caching, the absent space refused) moved to
tests/gates/sn/mesh/test_radial_characteristic_split_spaces.py.
"""

from __future__ import annotations


import numpy as np
import pytest

from orpheus.numerics.quadrature import Quadrature
from typing import NamedTuple

from orpheus.numerics.spaces.full_field_space import FullFieldSpace
from orpheus.numerics.spaces.radial_characteristic_space import (
    RadialCharacteristicBoundarySpace,
    RadialCharacteristicInteriorSpace,
)
from tests.gates.transport import _phase_space_cases as cases
from tests.gates.transport._phase_space_cases import angular_full_field_space
from orpheus.transport.fields.radial_characteristic_boundary_flux import (
    RadialCharacteristicBoundaryFlux,
)
from orpheus.transport.fields.radial_characteristic_interior_flux import (
    RadialCharacteristicInteriorFlux,
)
from orpheus.transport.full_field import Composite, FullField
from orpheus.transport.radial_characteristic_field import (
    RadialCharacteristicField,
)
from orpheus.transport.fields.radial_characteristic_boundary_source_sink import (
    RadialCharacteristicBoundarySourceSink,
)
from orpheus.transport.fields.radial_characteristic_interior_source_sink import (
    RadialCharacteristicInteriorSourceSink,
)

pytestmark = pytest.mark.foundation


class _PsiHalfSpaces(NamedTuple):
    """The ψ½ (radial-characteristic) spaces of one solid sphere, posed by hand (layer DAG Phase 4)."""

    interior: RadialCharacteristicInteriorSpace
    boundary: RadialCharacteristicBoundarySpace
    composite: FullFieldSpace


def _sphere() -> _PsiHalfSpaces:
    """System B's spaces on a 4-cell solid sphere of radius 1, 2 groups, GL4: the march levels are ``(0,)`` (one
    level carries; ``[M]`` the SN problem's ``radial_characteristic_levels`` on this mesh) and the cell volumes
    ``4/3 pi (r_o^3 - r_i^3)`` from the edges. The SN problem's mint of the same spaces is gated in
    ``tests/gates/sn/mesh/test_radial_characteristic_split_spaces.py``; this file's subject is the composite FIELD."""
    edges = np.linspace(0.0, 1.0, 5)
    volumes = 4.0 * np.pi / 3.0 * np.diff(edges ** 3)
    interior = RadialCharacteristicInteriorSpace.for_levels((0,), ng=2, nx=4, cell_volumes=volumes)
    boundary = RadialCharacteristicBoundarySpace.for_levels((0,), ng=2, nx=4, cell_volumes=volumes)
    return _PsiHalfSpaces(interior, boundary, FullFieldSpace.from_blocks(interior, boundary))


def _rand_composite(sn: _PsiHalfSpaces, seed: int) -> RadialCharacteristicField:
    """A random ψ½ FLUX composite (4e native — built via ``from_flat`` over a
    zero flux composite; the retired unified-leaf bridge is gone)."""
    rng = np.random.default_rng(seed)
    template = RadialCharacteristicField.flux_zeros(sn.composite)
    n = template.to_flat().size
    return RadialCharacteristicField.from_flat(rng.random(n) + 0.5, template)


# ── Construction (the multi-instantiation structural crux) ──


class TestCompositeConstruction:
    def test_is_a_composite_but_not_a_fullfield(self) -> None:
        c = RadialCharacteristicField.flux_zeros(_sphere().composite)
        if not isinstance(c, Composite):
            pytest.fail("System B is not a Composite")
        if isinstance(c, FullField):
            pytest.fail("System B must NOT be a FullField (a distinct composite)")

    def test_has_no_radial_characteristic_third_slot(self) -> None:
        # A pure 2-block composite — the ψ½ third block is FullField's, not this.
        c = RadialCharacteristicField.flux_zeros(_sphere().composite)
        if hasattr(c, "radial_characteristic"):
            pytest.fail("System B is pure 2-block; it carries no ψ½ third slot")

    def test_leaf_types_are_the_ray_split_loci(self) -> None:
        # The 3rd Composite instantiation, FIRST with non-Bulk/Boundary leaves.
        c = RadialCharacteristicField.flux_zeros(_sphere().composite)
        if type(c.interior) is not RadialCharacteristicInteriorFlux:
            pytest.fail(f"interior is {type(c.interior).__name__}")
        if type(c.boundary) is not RadialCharacteristicBoundaryFlux:
            pytest.fail(f"boundary is {type(c.boundary).__name__}")

    def test_source_members_construct_role_erased(self) -> None:
        r"""B.2b DP2: the slots bind the FIELD BASES (the FullField precedent),
        so an operator emission — the (interior, boundary) SOURCE pair that
        A_BA / B_b write — rides the SAME composite class as the flux state."""
        sn = _sphere()
        comp = RadialCharacteristicField(
            interior=RadialCharacteristicInteriorSourceSink.zeros(sn.interior),
            boundary=RadialCharacteristicBoundarySourceSink.zeros(sn.boundary),
        )
        if type(comp.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"interior member is {type(comp.interior).__name__}")
        if type(comp.boundary) is not RadialCharacteristicBoundarySourceSink:
            pytest.fail(f"boundary member is {type(comp.boundary).__name__}")


# ── The split-fidelity bridge (TestSplitFidelityBridge) RETIRED at Phase C 4e ──
#
# The whole class tested the ``from_unified`` / ``to_unified`` bridge between the
# retired unified single-buffer leaf (whose ``RadialCharacteristicField`` name
# the composite has since taken, 4e-e1b) and the composite — the
# "demotion licence" (``to_unified ∘ from_unified == id`` + role-preserving
# splits + the mixed/unknown-role bridge guards). The bridge and the unified
# leaf are GONE (the composite is the native representation), so the licence has
# no referent — the walk-slot rewrite landing bit-identically on the frozen
# baselines IS the discharged licence. Successors carrying the surviving claims:
#   * value fidelity / round-trip → ``TestCompositeAlgebra::
#     test_to_flat_from_flat_round_trip_preserves_type`` (the native flat
#     round-trip) + the ``_rand_composite`` ``from_flat`` build;
#   * role identity (flux / source members) → the member-type
#     rows in ``TestCompositeConstruction`` + ``TestCompositeAlgebra::
#     test_subtraction_is_same_typed_per_block`` (its pre-CS3 name was
#     ``test_subtraction_mints_a_displacement_composite_per_block``);
#   * the mixed-role "guard" DISSOLVED — role-erased slots make a mixed-role
#     composite LEGAL by construction (B.2b DP2); there is no bridge to refuse it.


# ── Intrinsic composite algebra on the ray carrier (multi-instantiation) ──


class TestCompositeAlgebra:
    def test_flux_plus_flux_propagates_per_block(self) -> None:
        """flux + flux is legal member-wise (flux lives in V — campaign 1
        CS3; the affine gate that used to refuse this here is retired)."""
        sn = _sphere()
        a, b = _rand_composite(sn, 4), _rand_composite(sn, 5)
        s = a + b
        np.testing.assert_array_equal(
            s.interior.values, a.interior.values + b.interior.values,
        )
        np.testing.assert_array_equal(
            s.boundary.values, a.boundary.values + b.boundary.values,
        )

    def test_subtraction_is_same_typed_per_block(self) -> None:
        """The difference keeps the SPLIT leaf types per block (until CS3 it
        minted the two RC displacement siblings — the per-block mint the
        composite torsor forced)."""
        sn = _sphere()
        a, b = _rand_composite(sn, 6), _rand_composite(sn, 7)
        d = a - b
        if type(d.interior) is not RadialCharacteristicInteriorFlux:
            pytest.fail(f"interior block became {type(d.interior).__name__}")
        if type(d.boundary) is not RadialCharacteristicBoundaryFlux:
            pytest.fail(f"boundary block became {type(d.boundary).__name__}")
        np.testing.assert_array_equal(
            d.interior.values, a.interior.values - b.interior.values,
        )

    def test_update_step_recovers_the_point(self) -> None:
        sn = _sphere()
        a, b = _rand_composite(sn, 8), _rand_composite(sn, 9)
        recovered = a + (b - a)  # plain V arithmetic, blockwise
        if type(recovered) is not RadialCharacteristicField:
            pytest.fail(f"vector recovery returned {type(recovered).__name__}")
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
        if type(back) is not RadialCharacteristicField:
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


# ── System B's member space (B.2b b2 — the family-blind FullFieldSpace reuse) ──


class TestRadialCharacteristicFieldSpace:
    r"""System B's member space (the SN problem's ``radial_characteristic_field_space``; here posed by hand).

    The DP1 ruling realized: the SAME family-blind ``FullFieldSpace`` class
    System A uses, instantiated over the two split ψ½ spaces. The identity
    name is the family rule (CS4b S4): ``full_field#<digest>`` folded from
    the MEMBER spaces' content — the pre-S4 ``"radial_characteristic"``
    role tag left the space name (G2.3: role is class identity; the role
    lives on ``RadialCharacteristicField``). The metric rows exercise the
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
        space = sn.composite
        if space is None:
            pytest.fail("sphere (carrying) must mint the composite space")
        n_i = sn.interior.shape[0]
        n_b = sn.boundary.shape[0]
        # The family prefix + a member-content digest (never pin the full
        # literal — R4); System B is DISCRIMINATED from System A by member
        # content, not by a role tag.
        if not space.name.startswith("full_field#"):
            pytest.fail(f"space name is {space.name!r}")
        system_a = angular_full_field_space(cases.sphere(nx=4, ng=2), Quadrature.gauss_legendre(4))
        if space == system_a:
            pytest.fail("System B's composite must differ from System A's")
        np.testing.assert_array_equal(space.shape, (n_i + n_b,))


    def test_metric_trio_dispatches_member_wise_on_the_real_composite(self) -> None:
        r"""apply_metric / apply_inverse_metric / inner_product on the REAL
        2-block composite ≡ direct split-space application per member,
        bitwise — the presence-dispatch 2-block arm, type-preserving. This
        row is the TOOTH for a kwarg-inversion mutation (always passing the
        seed kwarg REDs it with a TypeError)."""
        sn = _sphere()
        space = sn.composite
        x = _rand_composite(sn, 51)
        y = _rand_composite(sn, 52)

        gx = space.apply_metric(x)
        if type(gx) is not RadialCharacteristicField:
            pytest.fail(f"apply_metric returned {type(gx).__name__}")
        np.testing.assert_array_equal(
            gx.interior.values,
            sn.interior.apply_metric(x.interior.values),
        )
        np.testing.assert_array_equal(
            gx.boundary.values,
            sn.boundary.apply_metric(x.boundary.values),
        )

        ginv = space.apply_inverse_metric(x)
        np.testing.assert_array_equal(
            ginv.interior.values,
            sn.interior.apply_inverse_metric(
                x.interior.values,
            ),
        )
        np.testing.assert_array_equal(
            ginv.boundary.values,
            sn.boundary.apply_inverse_metric(
                x.boundary.values,
            ),
        )

        # inner_product is a SCALAR functional (Mode-12-blind to per-member
        # swaps on its own) — paired with the array rows above per the memo.
        np.testing.assert_array_equal(
            space.inner_product(x, y),
            sn.interior.inner_product(
                x.interior.values, y.interior.values,
            )
            + sn.boundary.inner_product(
                x.boundary.values, y.boundary.values,
            ),
        )

    def test_metric_preserves_the_member_role(self) -> None:
        r"""The rebuilt members keep their leaf classes (replace-based
        rebuild) — a source composite through G stays a source composite."""
        sn = _sphere()
        space = sn.composite
        src = RadialCharacteristicField(
            interior=RadialCharacteristicInteriorSourceSink.zeros(sn.interior),
            boundary=RadialCharacteristicBoundarySourceSink.zeros(sn.boundary),
        )
        out = space.apply_metric(src)
        if type(out.interior) is not RadialCharacteristicInteriorSourceSink:
            pytest.fail(f"interior became {type(out.interior).__name__}")
        if type(out.boundary) is not RadialCharacteristicBoundarySourceSink:
            pytest.fail(f"boundary became {type(out.boundary).__name__}")

    # ``test_full_field_three_slot_arm_is_byte_identical`` RETIRED at B.2d
    # d2 with the arm itself: a seed-carrying 3-slot FullField is
    # unrepresentable (the eviction), so the presence-dispatch arm it pinned
    # dissolved — the 2-block dispatch above is the ONLY arm, and System B's
    # metric rides THIS composite space (the member-wise trio rows).

    def test_flat_layout_is_coextensive_with_the_space_shape(self) -> None:
        r"""The elegance-carried B.2a reminder, pinned on the REAL members:
        ``to_flat().size == prod(space.shape)`` AND the flat ORDER is
        ``concat(interior, boundary)`` — the three offset spellings
        (``to_flat``, the space shape, the future ``system_slices``) name
        ONE layout."""
        sn = _sphere()
        space = sn.composite
        c = _rand_composite(sn, 54)
        flat = c.to_flat()
        np.testing.assert_array_equal(flat.size, int(np.prod(space.shape)))
        np.testing.assert_array_equal(
            flat,
            np.concatenate(
                [c.interior.values.ravel(), c.boundary.values.ravel()],
            ),
        )
