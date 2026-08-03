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

The proposed "fix" (pass ``level`` not ``p_idx``) is the actual bug: ``level`` is
the walk's ``mu_level_idx`` — ``None`` on the sphere, the ONLY carrying geometry —
so it would raise ``None not in (0,)``. ``p_idx`` is deliberate.

**Why this gate exists (net-new coverage).** The carrier gates
(:mod:`tests.sn.mesh.test_radial_characteristic_carrier`) pin *which* levels carry
(sphere → 1, cyl/cart → 0); this module pins the p_idx↔level *coordination* — the
one property that guarantees the direct-march writes each level's seed into that
level's slot without a crash or a silent-wrong-slot. It reddens if a future
curvilinear quadrature ever carries a multi-level / sparse-value set (the trigger
to re-examine the p_idx/level naming — still safe by the invariant, but worth
revisiting), or if the gate source and the space-key source ever diverge.

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
    if sn.curvature is not None and sn.reduced is not None \
            and sn.reduced.coord is CoordSystem.SPHERICAL:
        return 1
    li = getattr(sn.quad, "level_indices", None)
    return len(li) if li is not None else 1


# A broad battery covering every production curvilinear quadrature family.
_SPHERES = [("sphere_GL%d" % n, CoordSystem.SPHERICAL, Quadrature.gauss_legendre(n))
            for n in (2, 4, 6, 8, 10, 12)]
_CYL_PRODUCT = [("cyl_product_%dx%d" % (m, p), CoordSystem.CYLINDRICAL,
                 Quadrature.product(n_mu=m, n_phi=p))
                for m in (2, 4, 6) for p in (4, 8)]
_CYL_LS = [("cyl_LS_S%d" % s, CoordSystem.CYLINDRICAL, Quadrature.level_symmetric(s))
           for s in (2, 4, 8, 16)]
_ALL = _SPHERES + _CYL_PRODUCT + _CYL_LS


@pytest.mark.parametrize("name,coord,quad", _ALL, ids=[c[0] for c in _ALL])
def test_carrying_levels_are_never_sparse_or_multilevel(name, coord, quad):
    """No production curvilinear mesh has a multi-level / sparse-value carrier.

    Sphere → exactly ``(0,)`` (one M-M level, τ_raw,0 ∈ (0,1)); cylinder (both
    product and level-symmetric) → ``()`` (τ_raw,0 ∈ {0,1} on every level, dead
    seed). If a future quadrature ever carries >1 level, this reddens — the
    trigger to re-examine the p_idx/level coordination (still safe by the
    invariant below, but the naming should be revisited).
    """
    sn = _mesh_1d(coord, quad)
    carrying = sn.radial_characteristic_levels

    # (a) carrying set is always a subset of the sweep's level positions
    n_levels = _n_angular_levels(sn)
    _require(set(carrying) <= set(range(n_levels)),
             f"{name}: carrying levels {carrying} escape the level-position "
             f"range [0, {n_levels}) — the space key / sweep p_idx coordinate.")

    # (b) it is NEVER multi-level in production (the concern's premise)
    _require(len(carrying) <= 1,
             f"{name}: MULTI-LEVEL carrier {carrying} — the concern's premise "
             f"is now reachable; re-audit the p_idx/level coordination.")

    # (c) per the R12a trichotomy: sphere carries {0}, everything else empty
    if coord is CoordSystem.SPHERICAL:
        _require(carrying == (0,), f"{name}: sphere expected (0,), got {carrying}")
    else:
        _require(carrying == (), f"{name}: {coord} expected (), got {carrying}")


@pytest.mark.parametrize("name,coord,quad", _ALL, ids=[c[0] for c in _ALL])
def test_gate_source_is_the_space_key_source(name, coord, quad):
    """THE INVARIANT: the gate ``p_idx in seed_levels`` and the key validator
    ``level in space.levels`` read the SAME ``radial_characteristic_levels``.

    ``seed_levels = frozenset(mesh.radial_characteristic_levels)`` and
    ``space = mesh.radial_characteristic_space`` built ``for_levels(mesh.
    radial_characteristic_levels)``. Therefore the gate ONLY admits ``p_idx ∈
    space.levels`` — so ``slot_view(buf, p_idx, sign)``'s ``_slot_key`` guard can
    NEVER raise, regardless of contiguity/multiplicity.
    """
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
    """Sphere: the ONLY carrying geometry. ``p_idx == 0`` is the correct
    ``slot_view`` key; the loop's ``level`` is ``None`` (the walk's
    ``mu_level_idx``), so the proposed "fix" (pass ``level`` not ``p_idx``) would
    CRASH here — proving ``p_idx`` is deliberate, not a typo.
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
