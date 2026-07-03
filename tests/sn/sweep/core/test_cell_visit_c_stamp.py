r"""Production catcher for the ``CellVisit`` c- AND τ-stamp (Issue #236 Phase 2 B2/B3).

The Morel--Montry angular weight ``τ`` and its derived constants
``c_out = α_out/τ`` and ``c_in = (1-τ)/τ·α_out + α_in`` are an
angular-closure property.  B2 single-sourced the c-constants onto the
angular closure; B3 sourced ``τ`` itself off the closure for the LIVE
sweep (the ``DiamondDifference.update`` angular recurrence reads
``CellVisit.tau``).  All three are stamped onto every :class:`CellVisit`
at ONE production site, :meth:`SNMesh._make_cell_visit`, which reads the
closure's per-global-ordinate ``(N,)`` accessors
(``c_{in,out}_per_ordinate`` AND ``tau_per_ordinate``).
``global_ordinate`` is resolved differently per geometry:

* slab / sphere — the global ordinate IS ``direction_idx``;
* cylinder — ``level_indices[mu_level_idx][m]`` (the per-level→global
  permutation the closure also keys on).

A WRONG global-ordinate map in ``_make_cell_visit`` (e.g. a
``c_in``↔``c_out`` swap, a wrong ``tau_per_ordinate`` index, or a
cylinder level-index that scatters the wrong ``(M_p,)`` block) would
ship SILENTLY: the live matvec reads the closure's ``cell_contribution``
directly (NOT ``DD.residual``, the only reader of ``CellVisit.c_*``),
and the diamond / cell-balance fixtures build visits with a SURROGATE
stamp, so neither exercises the production stamp.  The ``τ`` stamp is
the live ``DD.update`` reader's sole input yet is read from a SEPARATE
accessor than the c-stamps — so a τ-only corruption (``×1.1`` / wrong
ordinate) leaves the c-assertions green and is invisible without a
dedicated τ arm (vv L11 Mode 11 — the named twins never call the
rewired reader; qa B3 mutation-1 proved a 0.2% cylinder scalar drift
ships with no committed catcher).

This module is the structurally-independent catcher: it walks a REAL
production ``dag_walk`` (sphere + a MULTI-LEVEL cylinder + a Cartesian
slab), and for EVERY visit asserts the stamped ``visit.c_in`` /
``visit.c_out`` equal the constants recomputed INLINE from that visit's
OWN ``streaming_terms`` (``α_in`` / ``α_out`` / ``τ_mm``), AND that the
stamped ``visit.tau`` equals that visit's OWN ``streaming_terms.tau_mm``
— all at 0-ULP.

The c reference is the HAND-TRANSCRIBED surrogate
``c_from_streaming_terms`` (NOT the production closure — importing the
closure's ``c`` to compare against itself would be tautological,
vv L11).  The τ reference is the GEOMETRY-produced ``st.tau_mm``, a
structurally-independent code path to the same BMC-2010 weight as the
closure's ``morel_montry_tau_per_level`` (the Leg-1 producer-equivalence
gate pins the two producers equal; THIS gate pins the stamp's
ordinate map, so a τ-stamp failure localizes to ``_make_cell_visit``).
Assertions use ``np.testing.assert_array_equal`` / ``pytest.fail``
(function calls, fire under ``-O``; the canonical ORPHEUS invocation
strips bare ``assert`` — vv Mode 8).

Mutation-self-check: swap ``c_in``↔``c_out`` (B2) OR perturb the τ
stamp (``tau=...*1.1``, B3) in ``_make_cell_visit`` and the sphere /
cylinder cases go RED; restore by re-edit (NEVER ``git checkout``, L28).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.mesh.augmented_mesh import SNMesh
from tests.sn._test_helpers import placeholder_materials
from tests.sn.sweep.core._c_surrogate import c_from_streaming_terms


def _assert_visit_stamp_matches_inline(visit, *, where: str) -> None:
    """Stamped ``visit.c_*`` / ``visit.tau`` equal the inline ``st`` values (0-ULP).

    Uses ``assert_array_equal`` (a function call — fires under ``-O``;
    vv Mode 8) on scalars wrapped in ``np.asarray`` so the same bit
    test applies to the slab neutral (``c = 0.0``, ``τ = 1.0``) and the
    curvilinear constants alike.

    The c reference is the hand-transcribed surrogate; the τ reference
    is the GEOMETRY-produced ``st.tau_mm`` (Issue #236 Phase 2 B3) — a
    structurally-independent producer to the closure's
    ``tau_per_ordinate`` (Leg-1 pins the producers equal; this pins the
    stamp's ordinate map).  ``visit.tau`` is read from a SEPARATE closure
    accessor than the c-stamps, so a τ-only corruption is invisible to
    the c-arm — this τ-arm is its dedicated catcher.
    """
    st = visit.streaming_terms
    c_in_inline, c_out_inline = c_from_streaming_terms(st)
    np.testing.assert_array_equal(
        np.asarray(visit.c_out), np.asarray(c_out_inline),
        err_msg=f"{where}: visit.c_out != α_out/τ inline recompute",
    )
    np.testing.assert_array_equal(
        np.asarray(visit.c_in), np.asarray(c_in_inline),
        err_msg=f"{where}: visit.c_in != (1-τ)/τ·α_out+α_in inline recompute",
    )
    np.testing.assert_array_equal(
        np.asarray(visit.tau), np.asarray(st.tau_mm),
        err_msg=f"{where}: visit.tau != streaming_terms.tau_mm (stamp map)",
    )


@pytest.mark.foundation
def test_sphere_production_stamp_matches_inline():
    """Every sphere ``dag_walk`` visit carries the closure's c at 0-ULP.

    Single-level (sphere: ``direction_idx`` IS the global ordinate);
    walks all ``N`` ordinates so the full closure ``(N,)`` accessor is
    exercised through the production stamp.
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

    n_visits = 0
    for n in range(quad.N):
        for visit in sn_mesh.dag_walk(ordinate_idx=n):
            _assert_visit_stamp_matches_inline(
                visit, where=f"sphere ordinate {n} cell {visit.cell_idx}",
            )
            n_visits += 1
    if n_visits != quad.N * sn_mesh.nx:
        pytest.fail(
            f"sphere walk produced {n_visits} visits, "
            f"expected {quad.N * sn_mesh.nx}"
        )


@pytest.mark.foundation
def test_multilevel_cylinder_production_stamp_matches_inline():
    """Every multi-level cylinder ``dag_walk`` visit carries c at 0-ULP.

    The cylinder global ordinate is the per-level→global permutation
    ``level_indices[mu_level_idx][m]``; a multi-level quadrature
    (``product(n_mu=2, ...)`` → 2 μ-levels) exercises the gather across
    DISTINCT level blocks — a wrong block scatter would mis-stamp a
    whole level and reds here.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.0, 1.5, 9),
        mat_ids=np.zeros(8, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"),
        bc_right=BC("vacuum"),
    )
    quad = Quadrature.product(n_mu=2, n_phi=4)
    sn_mesh = SNMesh(mesh, quad, placeholder_materials())

    level_indices = quad.level_indices
    if len(level_indices) < 2:
        pytest.fail(
            f"expected a multi-level cylinder quadrature; got "
            f"{len(level_indices)} μ-level(s)"
        )

    n_visits = 0
    for level_p in range(len(level_indices)):
        m_count = int(np.asarray(level_indices[level_p]).size)
        for within in range(m_count):
            for visit in sn_mesh.dag_walk(
                ordinate_idx=within, mu_level_idx=level_p,
            ):
                _assert_visit_stamp_matches_inline(
                    visit,
                    where=(
                        f"cyl level {level_p} within {within} "
                        f"cell {visit.cell_idx}"
                    ),
                )
                n_visits += 1
    if n_visits != quad.N * sn_mesh.nx:
        pytest.fail(
            f"cylinder walk produced {n_visits} visits, "
            f"expected {quad.N * sn_mesh.nx}"
        )


@pytest.mark.foundation
def test_slab_production_stamp_is_neutral_zero():
    """Every slab (Identity-closure) ``dag_walk`` visit has c_in == c_out == 0.0.

    Cartesian carries no angular redistribution (α ≡ 0, τ ≡ 1), so the
    identity closure stamps the neutral ``0.0`` — verified at 0-ULP via
    the inline recompute (which also yields ``0.0`` for slab ``st``) AND
    an explicit ``== 0.0`` so a non-zero leak reds.
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

    n_visits = 0
    for n in range(quad.N):
        for visit in sn_mesh.dag_walk(ordinate_idx=n):
            _assert_visit_stamp_matches_inline(
                visit, where=f"slab ordinate {n} cell {visit.cell_idx}",
            )
            np.testing.assert_array_equal(
                np.asarray([visit.c_in, visit.c_out]),
                np.zeros(2),
                err_msg=(
                    f"slab ordinate {n} cell {visit.cell_idx}: "
                    f"c not neutral zero"
                ),
            )
            n_visits += 1
    if n_visits != quad.N * sn_mesh.nx:
        pytest.fail(
            f"slab walk produced {n_visits} visits, "
            f"expected {quad.N * sn_mesh.nx}"
        )
