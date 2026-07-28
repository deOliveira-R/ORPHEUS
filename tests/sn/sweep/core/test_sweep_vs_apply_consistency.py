r"""Apply-vs-sweep consistency invariants for the starting-direction march (Issue #168).

Phase D fixed the apply-matvec path's Carlson coupled-pole seed; Phase F
backported the same fix to the SI/sweep path via
:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`.
The two paths share the Hébert §3.9.4 Eqs. (3.432)-(3.435) starting-direction
math, so they solve the SAME curvilinear fixed point.

#282 route (a) (#280 Phase 2.5d, 2026-07-04) retired the seed-strategy zoo:
the apply-path strategy object is gone, so the former apply-strategy ≡
sweep-helper equivalence test is vacuous and was dropped (Issue #248-style
retirement cleanup).  What SURVIVES is the direct
source-driven march :func:`carlson_inward_sweep_from_source` (now returning the
``(phi_cells, phi_face_final)`` tuple) and the end-to-end SI-vs-Krylov
eigenvalue consistency — the two genuine behavioural invariants this file pins.

This is a foundation-tagged test file: pure software-contract checks, no
theory-page equation labels.  ERR-026 tripwire — these tests stay green when
the architectural fix holds; if they break, the apply matvec and the SI sweep
have drifted from their shared starting-direction seed math.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.sweep.psi_half_angle_seed import (
    carlson_inward_sweep_from_source,
)


# ═══════════════════════════════════════════════════════════════════════
# Linearity of the source-driven starting-direction march
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.sentinel
@pytest.mark.catches("ERR-026")
def test_carlson_seed_helper_is_linear_in_Q_bar():
    r"""``carlson_inward_sweep_from_source`` is linear in ``Q_bar``.

    Required for the operator to remain linear: the seed propagates
    through the M-M angular recurrence as the half-angle face flux
    initial condition.  A non-linear seed would break
    :meth:`StreamingCollisionOperator.apply_transpose` / dense matrix probing.
    """
    nx = 8
    sigma_t_gx = np.full((1, nx), 0.5)
    dr = np.full(nx, 0.125)
    bc_outer = np.full((1,), 0.0)

    rng = np.random.default_rng(seed=20260512)
    Q_a = rng.standard_normal((1, nx))
    Q_b = rng.standard_normal((1, nx))
    alpha, beta = 2.3, -0.7

    # New API (#282 route (a)): the march returns (phi_cells, phi_face_final).
    seed_a, _ = carlson_inward_sweep_from_source(Q_a, sigma_t_gx, dr, bc_outer)
    seed_b, _ = carlson_inward_sweep_from_source(Q_b, sigma_t_gx, dr, bc_outer)
    seed_lin, _ = carlson_inward_sweep_from_source(
        alpha * Q_a + beta * Q_b, sigma_t_gx, dr, bc_outer,
    )

    np.testing.assert_allclose(
        seed_lin, alpha * seed_a + beta * seed_b,
        rtol=1e-13, atol=1e-13,
    )


@pytest.mark.foundation
@pytest.mark.catches("ERR-026")
def test_carlson_seed_helper_is_linear_in_bc():
    r"""``carlson_inward_sweep_from_source`` is linear in ``bc_outer_value``."""
    nx = 8
    sigma_t_gx = np.full((1, nx), 0.5)
    dr = np.full(nx, 0.125)
    Q_bar = np.zeros((1, nx))

    rng = np.random.default_rng(seed=20260512)
    bc_a = rng.standard_normal((1,))
    bc_b = rng.standard_normal((1,))
    alpha, beta = 1.7, 0.3

    seed_a, _ = carlson_inward_sweep_from_source(Q_bar, sigma_t_gx, dr, bc_a)
    seed_b, _ = carlson_inward_sweep_from_source(Q_bar, sigma_t_gx, dr, bc_b)
    seed_lin, _ = carlson_inward_sweep_from_source(
        Q_bar, sigma_t_gx, dr, alpha * bc_a + beta * bc_b,
    )

    np.testing.assert_allclose(
        seed_lin, alpha * seed_a + beta * seed_b,
        rtol=1e-13, atol=1e-13,
    )


# ═══════════════════════════════════════════════════════════════════════
# Apply ↔ sweep eigenvalue / scalar-flux consistency on
# homogeneous reflective curvilinear (the canonical Phase F probe)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.catches("ERR-026")
@pytest.mark.catches("ERR-053")
def test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere():
    r"""SI and Krylov inner solvers agree on homogeneous reflective sphere.

    The eigenvalue ``k = νΣ_f / Σ_a`` is shape-independent for
    homogeneous reflective; both paths must produce identical k_eff
    to numerical precision.  This is a foundation-level pin: any
    drift would indicate the apply-matvec and the SI sweep are no
    longer solving the same equation.

    Note: this is the WEAK form of apply-vs-sweep consistency.
    A future stronger form (post-Phase-F-extension) could pin
    eigenvalue agreement on HETEROGENEOUS MR — pre-Phase-F such a
    pin would have caught ERR-026 manifestation #6.  See the Phase F
    closeout memo §9.
    """
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.solver import solve_sn

    fuel = get_mixture("A", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom,
        region_meshes=(RegionMesh(n_cells=20),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)

    res_si = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver="source_iteration",
        keff_tol=1e-12, flux_tol=1e-10,
    )
    res_kr = solve_sn(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        inner_solver="krylov",
        keff_tol=1e-12, flux_tol=1e-10,
    )

    # Eigenvalue agreement: homogeneous reflective is the easy case
    # (flat eigenmode, k = νΣ_f/Σ_a), so this should hold to ~1e-6.
    assert abs(res_si.keff - res_kr.keff) < 1e-6, (
        f"SI keff = {res_si.keff:.10f}, Krylov keff = {res_kr.keff:.10f}, "
        f"diff = {res_si.keff - res_kr.keff:.3e}"
    )
