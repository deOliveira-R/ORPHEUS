r"""The ψ½ direct-march slot-coordination invariant: ``slot_view`` is keyed by
the level POSITION (``p_idx``), and the gate that admits a ``p_idx`` shares its
source with the space's key validator.

**The invariant (why the direct ψ½ solve is correct by construction).** The
direct-march block in
:mod:`orpheus.sn.loss_representation` (the ``for p_idx, level in enumerate(levels)``
solve) passes ``p_idx`` — the enumerate POSITION — to
:meth:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicInteriorSpace.slot_view`
``(buffer, level, sign)``, whose ``_slot_key`` raises unless the argument is a
member of ``self.levels``. This is **correct, not a coincidence**:

* the space's ``levels`` are level POSITIONS (``enumerate(raw)`` in
  :meth:`SNMesh.radial_characteristic_levels`), the same coordinate that keys the
  space slots AND indexes ``level_ordinates_list[p_idx]`` in the sweep;
* the gate ``seed_levels = frozenset(mesh.radial_characteristic_levels)`` and the
  space ``for_levels(mesh.radial_characteristic_levels)`` read the SAME tuple, so
  the gate ``p_idx in seed_levels`` admits ONLY a valid slot key — ``_slot_key``
  can never raise, for any carrying configuration (contiguous, sparse, or
  multi-level).

The proposed "fix" (pass ``level`` not ``p_idx``) is the actual bug: ``level``
is the walk's ``mu_level_idx`` — ``None`` on the sphere — so it would raise
``None not in (0,)``. ``p_idx`` is deliberate.

**Why this gate exists (net-new coverage).** The carrier gates
(:mod:`tests.sn.mesh.test_radial_characteristic_carrier`) pin *which* levels
carry; this module pins the p_idx↔level *coordination* — the one property that
guarantees the direct-march writes each level's seed into that level's slot
without a crash or a silent-wrong-slot.

**Q5.6.3 — the multi-level carrier arrived.** This module's original premise
("no production curvilinear mesh has a multi-level carrier; if a future
quadrature ever carries >1 level, this reddens — the trigger to re-examine")
FIRED exactly as designed: the folded cylinder (``Quadrature.folded_product``)
carries EVERY μ-level, so the multi-level case is now the admitted cylinder
NORM, not a hypothetical. The battery below therefore splits by family:
sphere → single slot ``(0,)``; folded cylinder → all level positions.  The
full-product / LS families REFUSE at construction since the flip landed —
their negatives live in ``test_cylindrical_quadrature_admission.py``.
The multi-level coordination itself is verified
DIRECTLY (distinct per-(p_idx, sign) markers — chosen asymmetric across
levels, because folded per-level data are bit-palindromic under the ξ-mirror
and a p ↔ mirror(p) slot swap would be invisible to symmetric markers).

Promoted from the numerics-investigator's p_idx-vs-level investigation (2026-07-07,
`refactor/sn-walk-unification`) — the verdict was NOT-A-BUG, and this pins the
invariant that makes it so. All gates raise via a ``_require`` helper (a raise
statement fires under the canonical ``python -O``; a bare ``assert`` would be
stripped — vv Mode 8).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import placeholder_materials

pytestmark = pytest.mark.foundation


def _require(cond: bool, msg: str) -> None:
    """``assert`` replacement that FIRES under ``python -O`` (a raise statement is
    not stripped, unlike the ``assert`` statement — vv Mode 8)."""
    if not cond:
        raise AssertionError(msg)


def _mesh_1d(coord: CoordSystem, quad, *, nx: int = 4, ng: int = 2) -> SNMesh:
    edges = np.linspace(0.0, 1.0, nx + 1)
    mat_ids = np.zeros(nx, dtype=int)
    # Cartesian carries a left (inner) BC; curvilinear's inner edge is the pole
    # (a regularity condition, not a BC face). Explicit branches — a **kwargs
    # dict spread confuses the type-checker's positional-arg inference.
    if coord is CoordSystem.CARTESIAN:
        mesh = Mesh1D(edges=edges, mat_ids=mat_ids, coord=coord,
                      bc_right=BC("reflective"), bc_left=BC("reflective"))
    else:
        mesh = Mesh1D(edges=edges, mat_ids=mat_ids, coord=coord,
                      bc_right=BC("reflective"))
    return SNMesh(mesh, quad, placeholder_materials(ng=ng))


def _n_angular_levels(sn: SNMesh) -> int:
    """The count of angular levels the sweep iterates (its ``levels`` variable).

    Mirrors loss_representation ~L3971-3979: sphere is the single-level ``[None]``
    case; every other curvilinear quadrature iterates
    ``range(len(quad.level_indices))``.
    """
    # (the retired ``sn.curvature is not None`` conjunct was redundant:
    # ``coord is SPHERICAL`` already implies non-Cartesian.  The chart is
    # read off the mesh since P4.1a retired the bundle's copy of it.)
    if sn.reduced is not None and sn.coord is CoordSystem.SPHERICAL:
        return 1
    li = getattr(sn.quad, "level_indices", None)
    return len(li) if li is not None else 1


# A broad battery covering every production curvilinear quadrature family.
# Each row carries its EXPECTED carrier as a function of the built mesh —
# the R12a predicate's per-family answer, asserted exactly.
_SPHERES = [("sphere_GL%d" % n, CoordSystem.SPHERICAL,
             Quadrature.gauss_legendre(n), lambda sn: (0,))
            for n in (2, 4, 6, 8, 10, 12)]
# Q5.6.3: the ADMITTED cylinder family — every μ-level carries (folded
# march starts are neither edge-node nor degenerate).  folded(4,6) puts
# the surviving bit-exact μ_r = 0 pure-azimuthal ordinate inside a
# carrying level (parent n_φ ≡ 2 mod 4).
_CYL_FOLDED = [("cyl_folded_%dx%d" % (m, p), CoordSystem.CYLINDRICAL,
                Quadrature.folded_product(n_mu=m, n_phi=p),
                lambda sn: tuple(range(len(sn.quad.level_indices))))
               for m, p in ((2, 4), (4, 8), (4, 6), (6, 12))]
# The pre-flip dead-seed families (full products, level-symmetric) held
# ``()`` rows here until the Q5.6.3 admission flip made their SNMesh
# construction REFUSE outright — their negatives now live one tier up,
# in ``test_cylindrical_quadrature_admission.py`` (single source: the
# refusal module owns "which rules are refused and why"; this battery
# owns "what an ADMITTED mesh's carrier looks like").
_ALL = _SPHERES + _CYL_FOLDED


@pytest.mark.parametrize(
    "name,coord,quad,expected_fn", _ALL, ids=[c[0] for c in _ALL],
)
def test_carrying_levels_match_the_family_prediction(name, coord, quad,
                                                     expected_fn):
    """The carrier set is EXACTLY the R12a predicate's per-family answer.

    Sphere → ``(0,)`` (one M-M level, τ_raw,0 ∈ (0,1)); folded cylinder →
    EVERY level position (the Q5.6.3 admitted family — the multi-level
    carrier this module's original premise called "a future quadrature");
    full-product / LS cylinder → ``()`` (τ_raw,0 ∈ {0,1} on every level,
    dead seed) until the flip refuses them outright.
    """
    sn = _mesh_1d(coord, quad)
    carrying = sn.radial_characteristic_levels

    # (a) carrying set is always a subset of the sweep's level positions
    n_levels = _n_angular_levels(sn)
    _require(set(carrying) <= set(range(n_levels)),
             f"{name}: carrying levels {carrying} escape the level-position "
             f"range [0, {n_levels}) — the space key / sweep p_idx coordinate.")

    # (b) the per-family exact carrier (multi-level is the folded NORM)
    expected = expected_fn(sn)
    _require(carrying == expected,
             f"{name}: carrier {carrying} != the family prediction {expected}")


@pytest.mark.parametrize(
    "name,coord,quad,expected_fn", _ALL, ids=[c[0] for c in _ALL],
)
def test_gate_source_is_the_space_key_source(name, coord, quad, expected_fn):
    """THE INVARIANT: the gate ``p_idx in seed_levels`` and the key validator
    ``level in space.levels`` read the SAME ``radial_characteristic_levels``.

    ``seed_levels = frozenset(mesh.radial_characteristic_levels)`` and
    ``space = mesh.radial_characteristic_space`` built ``for_levels(mesh.
    radial_characteristic_levels)``. Therefore the gate ONLY admits ``p_idx ∈
    space.levels`` — so ``slot_view(buf, p_idx, sign)``'s ``_slot_key`` guard can
    NEVER raise, regardless of contiguity/multiplicity.  Since Q5.6.3 the
    folded rows exercise this on GENUINE multi-level carriers, not only the
    sphere's single slot.
    """
    del expected_fn  # the invariant is family-independent
    sn = _mesh_1d(coord, quad)
    seed_levels = frozenset(sn.radial_characteristic_levels)
    space = sn.radial_characteristic_interior_space

    if not seed_levels:
        # Absence spelled None (never a zero-DOF space) — nothing to key.
        _require(space is None, f"{name}: empty carrier but space is not None")
        return

    _require(space is not None, f"{name}: carrying mesh but space is None")
    # Same source ⟹ gate ⊆ key-set: every admitted p_idx is a valid slot key.
    _require(seed_levels == frozenset(space.levels),
             f"{name}: gate source {seed_levels} != space key set "
             f"{frozenset(space.levels)} — a p_idx could pass the gate yet "
             f"crash cells_view, or vice-versa.")


def test_sphere_p_idx_is_the_correct_slot_key_and_level_is_None():
    """Sphere: the ``level is None`` carrying geometry. ``p_idx == 0`` is the
    correct ``slot_view`` key; the loop's ``level`` is ``None`` (the walk's
    ``mu_level_idx``), so the proposed "fix" (pass ``level`` not ``p_idx``)
    would CRASH here — proving ``p_idx`` is deliberate, not a typo.  (The
    folded cylinder — carrying since Q5.6.3 — has integer levels where the
    two coordinates AGREE numerically; the sphere is where they differ, so
    it remains the discriminating instance.)
    """
    sn = _mesh_1d(CoordSystem.SPHERICAL, Quadrature.gauss_legendre(4))
    space = sn.radial_characteristic_interior_space
    _require(space is not None and space.levels == (0,),
             f"sphere space.levels expected (0,), got "
             f"{None if space is None else space.levels}")
    assert space is not None  # narrow for the type-checker

    # The sweep's sphere branch is levels=[None]: enumerate([None]) →
    # (p_idx=0, level=None). p_idx=0 is the slot key.
    buf = np.zeros(space.shape)
    # p_idx = 0 is a valid key (no raise) and returns the (ng, nx) cells view.
    view = space.slot_view(buf, 0, -1)
    _require(view.shape == (space.ng, space.nx),
             f"cells_view(0,-1) shape {view.shape} != {(space.ng, space.nx)}")

    # The proposed "fix" passes ``level`` (== None for the sphere) — a crash.
    with pytest.raises((KeyError, TypeError, ValueError)):
        space.slot_view(buf, None, -1)  # type: ignore[arg-type]


def test_slot_coordination_writes_the_processed_levels_seed():
    """No silent-wrong-slot: the seed for the level at position ``p_idx`` (its
    ordinates are ``level_ordinates_list[p_idx]``) is written into the space slot
    keyed by that SAME position ``p_idx``. Verified structurally on the sphere
    (1 level ↔ 1 slot) — the coordinate is the level POSITION for the space slots,
    the sweep iteration, AND ``level_ordinates_list``.
    """
    sn = _mesh_1d(CoordSystem.SPHERICAL, Quadrature.gauss_legendre(4))
    space = sn.radial_characteristic_interior_space
    _require(space is not None, "sphere carrier space is None")
    assert space is not None  # narrow for the type-checker

    # Sphere: one M-M level (all N ordinates), one space slot keyed 0.
    _require(space.n_levels == 1 and space.levels == (0,),
             f"sphere expected 1 level (0,), got n_levels={space.n_levels} "
             f"levels={space.levels}")

    # Write a distinct marker through slot key p_idx=0 and read it back — the
    # -1 (inward) and +1 (outward) cells views are disjoint memory regions.
    buf = np.zeros(space.shape)
    space.slot_view(buf, 0, -1)[...] = 1.0
    space.slot_view(buf, 0, +1)[...] = 2.0
    np.testing.assert_array_equal(space.slot_view(buf, 0, -1), 1.0)
    np.testing.assert_array_equal(space.slot_view(buf, 0, +1), 2.0)


@pytest.mark.parametrize("n_mu,n_phi", [(2, 4), (4, 8), (4, 6)])
def test_folded_cylinder_multilevel_slots_are_disjoint_and_level_keyed(
    n_mu, n_phi,
):
    """Q5.6.3 — the multi-level coordination, verified DIRECTLY.

    The folded cylinder carries every μ-level, so ``slot_view`` must key
    ``n_mu × 2`` disjoint regions by ``(p_idx, sign)``.  Markers are chosen
    ASYMMETRIC across levels (``10·p_idx ± 1``): folded per-level data are
    bit-palindromic under the ξ-mirror (only signed μ_z distinguishes a
    level from its mirror partner — a Mode-12 hazard), so a p ↔ mirror(p)
    slot-permutation bug would be INVISIBLE to level-symmetric markers.
    Distinct per-level values make any permutation, overlap, or dropped
    write red here.
    """
    sn = _mesh_1d(CoordSystem.CYLINDRICAL,
                  Quadrature.folded_product(n_mu=n_mu, n_phi=n_phi))
    space = sn.radial_characteristic_interior_space
    _require(space is not None,
             f"folded({n_mu},{n_phi}) cylinder must carry — space is None")
    assert space is not None  # narrow for the type-checker
    _require(space.levels == tuple(range(n_mu)),
             f"folded({n_mu},{n_phi}) space.levels {space.levels} != "
             f"{tuple(range(n_mu))}")

    buf = np.zeros(space.shape)
    for p in space.levels:
        space.slot_view(buf, p, -1)[...] = 10.0 * p + 1.0
        space.slot_view(buf, p, +1)[...] = 10.0 * p + 2.0
    # Read back AFTER all writes: any slot overlap / permutation would have
    # clobbered an earlier marker.
    for p in space.levels:
        np.testing.assert_array_equal(
            space.slot_view(buf, p, -1), 10.0 * p + 1.0,
            err_msg=f"level {p} inward slot lost/permuted its marker",
        )
        np.testing.assert_array_equal(
            space.slot_view(buf, p, +1), 10.0 * p + 2.0,
            err_msg=f"level {p} outward slot lost/permuted its marker",
        )
