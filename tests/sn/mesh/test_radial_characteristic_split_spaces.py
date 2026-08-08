r"""B.1b — the split ψ½ ray spaces (System B's interior ⊕ boundary layout).

Space-level gates for the Phase-B (coupled-block campaign) split of the unified
``RadialCharacteristicSpace`` (keyed ``(level, sign, part)``; retired BY the
split, hence a literal and not a cross-reference) into the two per-``(level,
sign)`` spaces
:class:`RadialCharacteristicInteriorSpace` (the ``cells`` legs, ``G_sd = V_cell``)
and :class:`RadialCharacteristicBoundarySpace` (the ``corner`` legs, ``G = V(r =
R)``) — the substrate for posing the ψ½ ray as its own composite.

The split-fidelity oracle is a **hand-known arange fingerprint** (NOT
``split == unified`` — the split spaces and the unified space read the ONE
shared ``_radial_characteristic_legs`` walk, so a split-vs-unified check was
self-referential to a walk bug even while the unified space still existed; the
arange pins the layout order/offsets against a value computed independently of
production).

⚠ **Mandatory config-blindness fixture (test-architect
``coupled_operator_b1_split_verification.md``):** the sphere carries ONE level,
so a per-level-position offset bug (``2*pos → pos``) is INVISIBLE on it
(``pos ∈ {0}`` ⇒ ``2·0 = 0``). A MULTI-LEVEL synthetic ``for_levels((0, 2, 5))``
fixture (non-contiguous levels) is required — such a bug reds it while staying
green on the sphere, and the asymmetry IS the evidence.

vv Mode-8 discipline: ``np.testing.assert_*`` / ``pytest.raises`` only (no bare
``assert`` on numerical claims — the canonical invocation is ``python -O``).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.numerics.spaces.radial_characteristic_space import (
    RadialCharacteristicBoundarySpace as Boundary,
    RadialCharacteristicInteriorSpace as Interior,
    _RadialCharacteristicSubSpace,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation

#: A GRADED (non-uniform) cell-volume vector — so ``tile(V, ng)`` (interior) and
#: ``full(*, V[-1])`` (boundary) differ bit-exactly and the metric-partition
#: mutation bites (a uniform V would make the two metrics coincide).
_CV = np.array([0.065, 0.458, 1.244, 2.422])
_NG, _NX = 2, 4
_SPHERE_LEVELS: tuple[int, ...] = (0,)  # the production carrying instance
_MULTI_LEVELS: tuple[int, ...] = (0, 2, 5)  # the config-blindness fixture


def _legs(levels: tuple[int, ...]) -> list[tuple[int, int]]:
    """The canonical ``(level, sign)`` leg order: level ascending, sign ∈ (-1, +1)."""
    return [(lvl, s) for lvl in levels for s in (-1, +1)]


# ── Mesh builders (replicated from test_radial_characteristic_carrier) ──


def _mesh_1d(coord: CoordSystem, quad, *, nx: int = 4, ng: int = 2) -> SNMesh:
    edges = np.linspace(0.0, 1.0, nx + 1)
    mat_ids = np.zeros(nx, dtype=int)
    # Cartesian gets a reflective LEFT edge too; curvilinear leaves the r = 0
    # pole/axis to the Mesh1D default (the symmetry condition).
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
    """Sphere-GL: τ_raw,0 ≈ 0.42 ∈ (0, 1) — the CARRYING instance (1 level)."""
    return _mesh_1d(CoordSystem.SPHERICAL, Quadrature.gauss_legendre(4))


def _cyl_folded() -> SNMesh:
    """Cylinder on the admitted folded family — carrying on EVERY level.

    (Until Q5.6.3 this file held the two non-carrying cylinder fixtures;
    the admission flip made them unconstructible — their refusals are
    gated in ``test_cylindrical_quadrature_admission.py``.)"""
    return _mesh_1d(
        CoordSystem.CYLINDRICAL, Quadrature.folded_product(n_mu=2, n_phi=4),
    )


def _slab() -> SNMesh:
    """Cartesian: never carries."""
    return _mesh_1d(CoordSystem.CARTESIAN, Quadrature.gauss_legendre(4))


# ── Construction + identity ──


class TestSplitSpaceConstruction:
    def test_split_shapes_sum_to_the_combined_layout(self) -> None:
        i = Interior.for_levels(_MULTI_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        b = Boundary.for_levels(_MULTI_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        n_legs = len(_MULTI_LEVELS) * 2
        np.testing.assert_array_equal(i.shape, (n_legs * _NG * _NX,))
        np.testing.assert_array_equal(b.shape, (n_legs * _NG,))
        # The split shapes sum to the hand-known combined ψ½ layout size (per
        # leg: cells ng·nx + corner ng) — the successor of the retired unified
        # space's total (4e; the composite space is (ni + nb,)).
        np.testing.assert_equal(
            i.shape[0] + b.shape[0], n_legs * (_NG * _NX + _NG),
        )

    def test_identity_is_name_shape(self) -> None:
        i1 = Interior.for_levels(_SPHERE_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        i2 = Interior.for_levels(_SPHERE_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        b = Boundary.for_levels(_SPHERE_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        if i1 != i2 or hash(i1) != hash(i2):
            pytest.fail("(name, shape) identity: two equal-input builds differ")
        if i1 == b:
            pytest.fail("interior and boundary spaces must have distinct names")
        if i1.name != "radial_characteristic_interior":
            pytest.fail(f"interior name is {i1.name!r}")
        if b.name != "radial_characteristic_boundary":
            pytest.fail(f"boundary name is {b.name!r}")

    def test_empty_levels_rejected(self) -> None:
        with pytest.raises(ValueError, match="levels is empty"):
            Interior.for_levels((), ng=_NG, nx=_NX, cell_volumes=_CV)

    def test_nonascending_levels_rejected(self) -> None:
        with pytest.raises(ValueError, match="strictly ascending"):
            Boundary.for_levels((2, 0), ng=_NG, nx=_NX, cell_volumes=_CV)

    def test_cell_volumes_must_be_spd(self) -> None:
        with pytest.raises(ValueError, match="strictly positive"):
            Interior.for_levels(
                _SPHERE_LEVELS, ng=_NG, nx=_NX,
                cell_volumes=np.array([0.0, 1.0, 1.0, 1.0]),
            )

    def test_subspace_base_is_abstract(self) -> None:
        # The base carries no _PART; for_levels on it hits the abstract ClassVar.
        with pytest.raises(AttributeError):
            _RadialCharacteristicSubSpace.for_levels(
                _SPHERE_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV,
            )


# ── G-B1: split-fidelity via HAND-KNOWN arange fingerprint ──


class TestSplitFidelityArangeFingerprint:
    """Each split ``(level, sign)`` slot maps to a CONTIGUOUS, hand-known arange
    range (legs in level-ascending, sign ∈ (-1, +1) order). Independent of the
    unified space, so a shared-walk bug cannot hide. The multi-level fixture is
    the load-bearing case (the sphere is offset-blind at one level)."""

    @pytest.mark.parametrize("levels", [_SPHERE_LEVELS, _MULTI_LEVELS])
    def test_interior_slots_are_hand_known_aranges(self, levels) -> None:
        i = Interior.for_levels(levels, ng=_NG, nx=_NX, cell_volumes=_CV)
        buf = np.arange(i.shape[0], dtype=float)
        block = _NG * _NX
        for leg_idx, (level, sign) in enumerate(_legs(levels)):
            expected = np.arange(
                leg_idx * block, (leg_idx + 1) * block, dtype=float,
            ).reshape(_NG, _NX)
            np.testing.assert_array_equal(i.slot_view(buf, level, sign), expected)

    @pytest.mark.parametrize("levels", [_SPHERE_LEVELS, _MULTI_LEVELS])
    def test_boundary_slots_are_hand_known_aranges(self, levels) -> None:
        b = Boundary.for_levels(levels, ng=_NG, nx=_NX, cell_volumes=_CV)
        buf = np.arange(b.shape[0], dtype=float)
        block = _NG
        for leg_idx, (level, sign) in enumerate(_legs(levels)):
            expected = np.arange(
                leg_idx * block, (leg_idx + 1) * block, dtype=float,
            )
            np.testing.assert_array_equal(b.slot_view(buf, level, sign), expected)

    def test_unknown_level_and_sign_reject(self) -> None:
        i = Interior.for_levels(_SPHERE_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        buf = np.zeros(i.shape[0])
        with pytest.raises(KeyError, match="carries no"):
            i.slot_view(buf, 7, -1)  # level 7 not carrying
        with pytest.raises(ValueError, match="sign must be"):
            i.slot_view(buf, 0, 0)  # sign 0 illegal


# ── G-B2: recompose exact — the "split IS the unified, re-typed" contract ──


# The G-B3 "split ↔ unified buffer" recompose licence (``test_split_then_
# recompose_is_the_unified_buffer``) RETIRED at Phase C 4e with its referent:
# the unified ``RadialCharacteristicSpace`` buffer is gone. Its "the split
# loses nothing" claim is now carried by the successors that do NOT reference
# the unified: the arange-fingerprint slots (``TestSplitFidelityArangeFingerprint``
# — the exact per-slot layout) + the composite ``to_flat``/``from_flat``
# round-trip in the composite carrier suite (``test_radial_characteristic_field``).


# ── G-B3: metric partition + conservation (ERR-067-adjacent gauge) ──


class TestMetricPartition:
    def test_interior_metric_is_tile_v_cell(self) -> None:
        i = Interior.for_levels(_MULTI_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        expected = np.concatenate([np.tile(_CV, _NG)] * (len(_MULTI_LEVELS) * 2))
        np.testing.assert_array_equal(i.inner_product_weights, expected)

    def test_boundary_metric_is_v_at_outer_radius(self) -> None:
        b = Boundary.for_levels(_MULTI_LEVELS, ng=_NG, nx=_NX, cell_volumes=_CV)
        expected = np.full(len(_MULTI_LEVELS) * 2 * _NG, _CV[-1])
        np.testing.assert_array_equal(b.inner_product_weights, expected)

    # ``test_unified_metric_is_the_interleave_of_the_split`` RETIRED at 4e with
    # its referent (the unified metric is gone). The split gauges are pinned
    # directly, without the unified, by ``test_interior_metric_is_tile_v_cell``
    # + ``test_boundary_metric_is_v_at_outer_radius`` above (single-sourced from
    # ``_radial_characteristic_legs``, so no cross-space drift is spellable).


# ── G-C2: mesh presence = ray-carrying (the new SNMesh properties) ──


class TestMeshPresence:
    def test_sphere_carries_both_split_spaces_aligned_with_composite(self) -> None:
        sn = _sphere()
        i = sn.radial_characteristic_interior_space
        b = sn.radial_characteristic_boundary_space
        c = sn.radial_characteristic_field_space
        if i is None or b is None or c is None:
            pytest.fail("sphere must carry the ψ½ spaces (R12a)")
        if not (i.levels == b.levels):
            pytest.fail(f"level mismatch: {i.levels} / {b.levels}")
        # The split shapes sum to the composite space's flat size (4e — the
        # composite is the successor of the retired unified space).
        np.testing.assert_array_equal(
            (i.shape[0] + b.shape[0],), c.shape,
        )

    def test_folded_cylinder_carries_both_split_spaces_aligned_with_composite(
        self,
    ) -> None:
        """Q5.6.3 — the multi-level twin of the sphere presence gate."""
        sn = _cyl_folded()
        i = sn.radial_characteristic_interior_space
        b = sn.radial_characteristic_boundary_space
        c = sn.radial_characteristic_field_space
        if i is None or b is None or c is None:
            pytest.fail("folded cylinder must carry the ψ½ spaces (R12a)")
        n_levels = len(sn.quad.level_indices)
        if i.levels != tuple(range(n_levels)) or b.levels != i.levels:
            pytest.fail(f"level mismatch: {i.levels} / {b.levels}")
        np.testing.assert_array_equal(
            (i.shape[0] + b.shape[0],), c.shape,
        )

    @pytest.mark.parametrize("mesh_fn", [_slab])
    def test_noncarrying_meshes_have_no_split_spaces(self, mesh_fn) -> None:
        sn = mesh_fn()
        # The presence trichotomy is single-sourced: no composite → no split.
        if sn.radial_characteristic_field_space is not None:
            pytest.skip("mesh unexpectedly carries a seed level")
        if sn.radial_characteristic_interior_space is not None:
            pytest.fail("non-carrying mesh has an interior split space")
        if sn.radial_characteristic_boundary_space is not None:
            pytest.fail("non-carrying mesh has a boundary split space")

    def test_split_spaces_are_cached(self) -> None:
        sn = _sphere()
        if sn.radial_characteristic_interior_space is not sn.radial_characteristic_interior_space:
            pytest.fail("interior space is not cached (one instance per mesh)")
        if sn.radial_characteristic_boundary_space is not sn.radial_characteristic_boundary_space:
            pytest.fail("boundary space is not cached")
