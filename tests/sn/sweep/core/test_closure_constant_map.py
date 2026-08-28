r"""Production catcher for the closure's per-ordinate ``(N,)`` constant map.

P4.9a re-derivation of the retired ``CellVisit`` c/τ-stamp catcher
(Issue #236 Phase 2 B2/B3/C).  The stamp died with the un-weld: the mesh
no longer copies closure data onto visits — the SN walk assembles the
scheme's angular arguments from the closure's OWN per-global-ordinate
accessors (``c_{in,out}_per_ordinate`` / ``tau_per_ordinate``), and the
march indexes ``tau_per_ordinate`` at the global ordinate.  So the
ordinate-MAP hazard this file existed for did not disappear, it moved
one producer up: a wrong per-level→global GATHER inside the closure's
construction (a ``c_in``↔``c_out`` swap, a wrong ``tau`` index, a
cylinder level block scattered to the wrong ordinates) would now ship
through every consumer of the ``(N,)`` accessors — the cache populator,
the walk's degenerate assembly, and the march itself.

This module is the structurally-independent catcher for that map: for
EVERY (cell, global ordinate) of a sphere, a MULTI-LEVEL cylinder, and a
Cartesian slab, the closure's accessor values equal the constants
recomputed via the INDEPENDENT surrogate — τ from
``contamination.morel_montry_weights`` (a DIFFERENT code path to the
same BMC-2010 weight; vv L11 — NOT the closure's own
``morel_montry_tau_per_level``, which would be tautological), α from the
operator's surviving dome — all at 0-ULP.  The per-cell loop also pins
the surrogate's cell-independence (the accessors are per-ordinate only).

The Leg-1 producer-equivalence gate (``test_tau_producer_equivalence.py``)
pins the closure τ PRODUCER against the same independent reference; THIS
gate pins the per-ordinate map every consumer reads.  Assertions use
``np.testing.assert_array_equal`` / ``pytest.fail`` (function calls,
fire under ``-O``; the canonical ORPHEUS invocation strips bare
``assert`` — vv Mode 8).

Mutation-self-check: swap the c-gather (``c_in``↔``c_out``) or perturb
the τ gather (``×1.1`` / wrong ordinate) in the closure's construction
and the sphere / cylinder cases go RED; restore by re-edit (NEVER
``git checkout``, L28).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.sn.mesh.reduced_operator import (
    cylindrical_streaming,
    slab_streaming,
    spherical_streaming,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import placeholder_materials
from tests.sn.sweep.core._c_surrogate import (
    c_from_constants,
    mm_constants_for_ordinate,
)


def _assert_closure_map_matches_inline(
    closure,
    op,
    cell_idx: int,
    global_ordinate: int,
    direction_idx: int,
    *,
    mu_level_idx=None,
    where: str,
) -> None:
    """Closure accessors at ``global_ordinate`` equal the independent
    values for this (cell, ordinate) at 0-ULP.

    ``direction_idx`` is what the surrogate keys on (within-level for the
    cylinder); ``global_ordinate`` is what the ``(N,)`` accessors key on —
    asserting across both pins the per-level→global gather.
    """
    tau_ref, alpha_in, alpha_out = mm_constants_for_ordinate(
        op, cell_idx, direction_idx, mu_level_idx=mu_level_idx,
    )
    c_in_inline, c_out_inline = c_from_constants(tau_ref, alpha_in, alpha_out)
    np.testing.assert_array_equal(
        np.asarray(closure.c_out_per_ordinate[global_ordinate]),
        np.asarray(c_out_inline),
        err_msg=f"{where}: c_out accessor != α_out/τ independent recompute",
    )
    np.testing.assert_array_equal(
        np.asarray(closure.c_in_per_ordinate[global_ordinate]),
        np.asarray(c_in_inline),
        err_msg=(
            f"{where}: c_in accessor != (1-τ)/τ·α_out+α_in independent "
            f"recompute"
        ),
    )
    np.testing.assert_array_equal(
        np.asarray(closure.tau_per_ordinate[global_ordinate]),
        np.asarray(tau_ref),
        err_msg=f"{where}: tau accessor != independent τ reference (map)",
    )


@pytest.mark.foundation
def test_sphere_closure_map_matches_inline():
    """Every sphere (cell, ordinate) pair: the closure ``(N,)`` accessors
    equal the independent surrogate at 0-ULP.

    Single-level (sphere: ``direction_idx`` IS the global ordinate);
    loops all ``N`` ordinates × all cells so the full accessor and the
    surrogate's cell-independence are both exercised.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 2.0, 11),
        mat_ids=np.zeros(10, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(8)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    op = spherical_streaming(mesh, quad)
    closure = sn_mesh.pole_angular_closure

    for n in range(quad.N):
        for cell in range(sn_mesh.nx):
            _assert_closure_map_matches_inline(
                closure, op, cell, n, n,
                where=f"sphere ordinate {n} cell {cell}",
            )


@pytest.mark.foundation
def test_multilevel_cylinder_closure_map_matches_inline():
    """Every multi-level cylinder (cell, ordinate): accessors at 0-ULP.

    The cylinder global ordinate is the per-level→global permutation
    ``level_indices[mu_level_idx][m]``; a multi-level quadrature
    (2 μ-levels) exercises the gather across DISTINCT level blocks — a
    wrong block scatter would mis-map a whole level and reds here.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.5, 9),
        mat_ids=np.zeros(8, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.folded_product(n_mu=2, n_phi=4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    op = cylindrical_streaming(mesh, quad)
    closure = sn_mesh.pole_angular_closure

    level_indices = quad.level_indices
    if len(level_indices) < 2:
        pytest.fail(
            f"expected a multi-level cylinder quadrature; got "
            f"{len(level_indices)} μ-level(s)"
        )

    n_pairs = 0
    for level_p in range(len(level_indices)):
        block = np.asarray(level_indices[level_p])
        for within in range(int(block.size)):
            global_n = int(block[within])
            for cell in range(sn_mesh.nx):
                _assert_closure_map_matches_inline(
                    closure, op, cell, global_n, within,
                    mu_level_idx=level_p,
                    where=(
                        f"cyl level {level_p} within {within} "
                        f"(global {global_n}) cell {cell}"
                    ),
                )
            n_pairs += 1
    if n_pairs != quad.N:
        pytest.fail(
            f"cylinder level blocks covered {n_pairs} ordinates, "
            f"expected {quad.N}"
        )


@pytest.mark.foundation
def test_slab_closure_map_is_neutral():
    """The Cartesian identity closure's map is the neutral element.

    ``c_in == c_out == 0.0`` and ``τ == 1.0`` at every ordinate —
    verified via the inline recompute (which also yields the neutral
    values for slab) AND an explicit literal so a non-zero leak reds.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("vacuum"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.gauss_legendre(4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())
    op = slab_streaming(mesh, quad)
    closure = sn_mesh.pole_angular_closure

    for n in range(quad.N):
        _assert_closure_map_matches_inline(
            closure, op, 0, n, n,
            where=f"slab ordinate {n}",
        )
    np.testing.assert_array_equal(
        closure.c_in_per_ordinate, np.zeros(quad.N),
        err_msg="slab c_in not neutral zero",
    )
    np.testing.assert_array_equal(
        closure.c_out_per_ordinate, np.zeros(quad.N),
        err_msg="slab c_out not neutral zero",
    )
    np.testing.assert_array_equal(
        closure.tau_per_ordinate, np.ones(quad.N),
        err_msg="slab tau not neutral one",
    )
