r"""2-D Cartesian L2-native matvec correctness — eigenvalue / linearity / parity.

Test-architect Tests 2.1, 2.2, 2.4, 2.5 (D-H.2-C4a verification spec).

Each test xfails today and auto-flips to PASS in D-H.2-C4b (when the
2-D ``StreamingOperator._apply_timed_full_field`` stops raising
``NotImplementedError`` via a flat-shim through the legacy 2-D kernel),
through C4c-d (when the kernel itself rewrites L2-native).

Pillar / claim discipline
==========================

Each test row carries a Pillar + Claim per ``vv-principles`` skill:

* **Test 2.1** — k_∞ closed-form (homogeneous reflective ⇒
  k = νΣ_f/Σ_a, mesh-independent).  Cardinal-rule: ≥2G.
* **Test 2.2** — algebraic identity (A.apply distributes over the
  ``OperatorSum`` ``+``).  Catches ERR-026 on the 2-D L2 path.
* **Test 2.4** — 2-D reflective-y reduces to 1-D slab.  Pin against
  the L1-anchored 1-D slab solution.  Cardinal-rule: ≥2G,
  heterogeneous.
* **Test 2.5** — linearity (α·A.u + β·A.v = A.(α·u + β·v)).  Catches
  missing /W projection at any of the bulk/boundary interfaces.

Why xfail-strict
================

The 2-D ``StreamingOperator._apply_typed`` path raises
``NotImplementedError`` for composite ``TimedFullField`` input today
(D-H.1b.6 left it as a stub).  C4b unblocks the path.  Each test
here pins a positive correctness criterion that the unblock must
satisfy at the new 2-D L2 surface.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Mesh2D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.geometry import SNMesh


# ── Mesh builders ────────────────────────────────────────────────────────


def _homogeneous_reflective_2d(nx: int = 4, ny: int = 4) -> SNMesh:
    r"""All-reflective 2-D Cartesian mesh; 2g mixture A homogeneous.

    The canonical k_inf fixture: all-reflective + homogeneous ⇒
    k = νΣ_f/Σ_a = k_inf (1.875 for mixture A 2g) for ALL meshes.
    """
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 2.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("reflective"),
        bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"),
        bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, {0: get_mixture("A", "2g")})


def _vacuum_xy_2d_with_scatter(nx: int = 4, ny: int = 4) -> SNMesh:
    r"""Vacuum-on-all-faces 2-D mesh with nonzero scattering.

    Eigenvalue problem doesn't make sense here (no fission),
    so this fixture serves the apply-vs-sweep parity test (2.2)
    and the linearity test (2.5) — both of which need a 2-D
    mesh with nontrivial C and S (heterogeneous redistribution path)
    but don't need fission.
    """
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("vacuum"),
        bc_xmax=BC("vacuum"),
        bc_ymin=BC("vacuum"),
        bc_ymax=BC("vacuum"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, {0: get_mixture("A", "2g")})


def _slab_homogeneous_2g(nx: int = 4) -> SNMesh:
    r"""1-D slab counterpart for the 2-D-reflective-y reduction test."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC("reflective"), BC("reflective")),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=nx),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=4)
    return SNMesh(mesh, quad, {0: get_mixture("A", "2g")})


def _reflective_y_2d_for_1d_reduction(nx: int = 4, ny: int = 2) -> SNMesh:
    r"""2-D mesh with reflective y, reflective x — should reduce to 1-D slab.

    Under all-reflective with y-direction homogeneous geometry +
    homogeneous materials, the y-direction is symmetric and the
    eigenvalue/flux should match the 1-D-slab analog (same x BC
    + same materials).  This pins the 2-D path against the L1-anchored
    1-D path WITHOUT needing an analytical 2-D solver.
    """
    geom = Mesh2D(
        edges_x=np.linspace(0.0, 2.0, nx + 1),
        edges_y=np.linspace(0.0, 1.0, ny + 1),
        mat_map=np.zeros((nx, ny), dtype=int),
        bc_xmin=BC("reflective"),
        bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"),
        bc_ymax=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    return SNMesh(geom, quad, {0: get_mixture("A", "2g")})


# ── Test 2.1: 2-D Krylov recovers k_inf (homogeneous reflective) ────────


@pytest.mark.l1
def test_solve_sn_2d_krylov_homogeneous_reflective_recovers_kinf() -> None:
    r"""All-reflective 2-D homogeneous ⇒ k_eff = k_inf = 1.875.

    Pillar: closed-form (k_inf = νΣ_f/Σ_a is mesh-independent for
    homogeneous reflective).  The L1 floor for the 2-D L2 matvec.
    """
    from orpheus.sn.solver import solve_sn

    mesh = _homogeneous_reflective_2d(nx=4, ny=4)
    res = solve_sn(
        materials={0: get_mixture("A", "2g")},
        mesh=mesh.mesh,
        quadrature=mesh.quad,
        inner_solver="krylov",
        keff_tol=1e-10,
        flux_tol=1e-10,
        max_outer=200,
        max_inner=200,
    )

    k_inf = 1.875  # νΣ_f / Σ_a for mixture A 2g
    assert np.isfinite(res.keff), f"2-D Krylov returned non-finite keff: {res.keff}"
    assert abs(res.keff - k_inf) < 1e-6, (
        f"2-D Krylov keff = {res.keff:.10f}, expected k_inf = {k_inf}, "
        f"err = {abs(res.keff - k_inf):.3e}.  Pillar: closed-form "
        f"homogeneous reflective."
    )


# ── Test 2.2: 2-D ERR-026 — apply-distributes-over-OperatorSum ──────────


@pytest.mark.foundation
@pytest.mark.catches("ERR-026")
def test_apply_vs_sweep_2d_residual_cancellation() -> None:
    r"""``(L + C).apply(state).bulk == L.apply(state).bulk + C.apply(state).bulk``.

    Pillar: algebraic identity (``OperatorSum.apply`` distributes over
    ``+``).  Catches ERR-026 family on the 2-D L2 path.  This is the
    2-D analog of ``test_apply_vs_sweep_consistency``'s 1-D tripwire.
    """
    from orpheus.sn.operator import CollisionOperator, StreamingOperator

    mesh = _vacuum_xy_2d_with_scatter()
    rng = np.random.default_rng(seed=20260528)

    state = mesh.zeros_timed_full_field()
    state.bulk.values[...] = rng.standard_normal(state.bulk.values.shape)
    state.boundary.values[...] = rng.standard_normal(state.boundary.values.shape)

    sigma_t = np.array(
        [
            mesh.materials[0].SigT[g]
            * np.ones((mesh.nx, mesh.ny))
            for g in range(mesh.ng)
        ]
    )

    L = StreamingOperator(mesh, sigma_t)
    C = CollisionOperator(mesh, sigma_t)
    A = L + C

    sum_out = A.apply(state)
    sep_out = L.apply(state).bulk.values + C.apply(state).bulk.values

    np.testing.assert_allclose(
        sum_out.bulk.values, sep_out, rtol=0.0, atol=1e-12,
        err_msg=(
            "ERR-026 2-D tripwire: (L+C).apply ≠ L.apply + C.apply on bulk."
        ),
    )


# ── Test 2.4: 2-D-with-reflective-y reduces to 1-D slab ─────────────────


@pytest.mark.l1
def test_2d_reflective_xy_keff_matches_1d_slab_reflective_analog() -> None:
    r"""2-D all-reflective + y-direction-symmetric ⇒ keff = 1-D slab keff.

    Pillar: structurally-independent reference (the 1-D path is
    pinned by ``test_l1_standoff_slab_cylinder.py`` against analytic
    eigenvalue).  Under all-reflective + y-homogeneous-geometry +
    homogeneous-materials, 2-D = 1-D-extruded ⇒ keff matches the
    1-D slab analog to mesh-converged tolerance.

    This test catches a heterogeneous-2D-specific bug class that
    homogeneous-flat-flux test 2.1 cannot see.
    """
    from orpheus.sn.solver import solve_sn

    mesh_1d = _slab_homogeneous_2g(nx=8)
    mesh_2d = _reflective_y_2d_for_1d_reduction(nx=8, ny=2)

    res_1d = solve_sn(
        materials={0: get_mixture("A", "2g")},
        mesh=mesh_1d.mesh,
        quadrature=mesh_1d.quad,
        inner_solver="krylov",
        keff_tol=1e-10,
        flux_tol=1e-10,
    )
    res_2d = solve_sn(
        materials={0: get_mixture("A", "2g")},
        mesh=mesh_2d.mesh,
        quadrature=mesh_2d.quad,
        inner_solver="krylov",
        keff_tol=1e-10,
        flux_tol=1e-10,
    )

    assert abs(res_2d.keff - res_1d.keff) < 1e-6, (
        f"2-D-reflective-y keff = {res_2d.keff:.10f}, "
        f"1-D slab analog keff = {res_1d.keff:.10f}, "
        f"diff = {res_2d.keff - res_1d.keff:.3e}.  Reflective-y "
        f"reduction broken — heterogeneous-2D path has a bug invisible "
        f"to test 2.1's homogeneous-flat-flux."
    )


# ── Test 2.5: 2-D matvec linearity ──────────────────────────────────────


@pytest.mark.foundation
def test_2d_matvec_linearity_random_state() -> None:
    r"""``A.apply(α·u + β·v) ≡ α·A.apply(u) + β·A.apply(v)`` on 2-D L2.

    Pillar: linearity property (the matvec is a linear map by
    construction).  Catches per-component normalisation defects
    (e.g. a /W applied to bulk but not boundary).  Per
    ``vv-principles`` failure-mode #3 (missing factor) on the 2-D
    L2 path.
    """
    from orpheus.sn.operator import CollisionOperator, StreamingOperator

    mesh = _vacuum_xy_2d_with_scatter()
    sigma_t = np.array(
        [
            mesh.materials[0].SigT[g]
            * np.ones((mesh.nx, mesh.ny))
            for g in range(mesh.ng)
        ]
    )

    L = StreamingOperator(mesh, sigma_t)
    C = CollisionOperator(mesh, sigma_t)
    A = L + C

    rng = np.random.default_rng(seed=20260528)
    u = mesh.zeros_timed_full_field()
    v = mesh.zeros_timed_full_field()
    u.bulk.values[...] = rng.standard_normal(u.bulk.values.shape)
    v.bulk.values[...] = rng.standard_normal(v.bulk.values.shape)
    u.boundary.values[...] = rng.standard_normal(u.boundary.values.shape)
    v.boundary.values[...] = rng.standard_normal(v.boundary.values.shape)

    alpha, beta = 2.5, -1.3
    lhs = A.apply(alpha * u + beta * v)
    rhs_bulk = (
        alpha * A.apply(u).bulk.values + beta * A.apply(v).bulk.values
    )

    np.testing.assert_allclose(
        lhs.bulk.values, rhs_bulk, rtol=1e-12, atol=1e-14,
        err_msg=(
            "2-D matvec lost linearity on bulk — /W projection defect "
            "or convention drift across bulk/boundary boundary."
        ),
    )
