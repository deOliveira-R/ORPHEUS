r"""Phase F apply-vs-sweep consistency invariants (Issue #168).

Phase D fixed the apply-matvec path's Carlson coupled-pole seed
(:func:`~orpheus.sn.operator.transport_operator_matvec_spherical`).
Phase F backports the same fix to the SI/sweep path
(:func:`~orpheus.sn.sweep._sweep_1d_spherical`,
:func:`~orpheus.sn.sweep._sweep_1d_cylindrical`) via
:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`.
The two paths share the same Hébert §3.9.4 Eqs. (3.432)-(3.435) seed
math — they ought to produce identical seeds on identical inputs.

These tests pin that structural-alignment invariant.  A future
regression where one path drifts from the canonical seed math but
the other does not should be caught here.

This is a foundation-tagged test file: pure software-contract checks,
no theory-page equation labels.  ERR-026 tripwire — these tests stay
green when the architectural fix lands; if they break, the apply
matvec and the SI sweep have drifted from their shared Carlson seed.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.spatial.psi_half_angle_seed import (
    CarlsonInwardSweep,
    CarlsonSweepContext,
    carlson_inward_sweep_from_source,
)


# ═══════════════════════════════════════════════════════════════════════
# Phase F apply ↔ sweep Carlson seed structural alignment
# ═══════════════════════════════════════════════════════════════════════


def _build_inputs(nx, sigma_t_value, bc_outer_value):
    """Build matching apply-path and sweep-path inputs."""
    sigma_t_gx = np.full((1, nx), sigma_t_value)
    dr = np.full(nx, 1.0 / nx)
    bc_outer = np.full((1,), bc_outer_value)
    return sigma_t_gx, dr, bc_outer


@pytest.mark.foundation
@pytest.mark.catches("ERR-026")
@pytest.mark.parametrize("nx", [4, 8, 16])
@pytest.mark.parametrize("sigma_t_value", [0.1, 0.5, 1.5])
@pytest.mark.parametrize("bc_outer_value", [0.0, 0.3, 1.0])
@pytest.mark.parametrize("psi_const", [1.0, 2.5])
def test_carlson_seed_apply_sweep_equivalence_flat_psi(
    nx, sigma_t_value, bc_outer_value, psi_const,
):
    r"""apply-path ``CarlsonInwardSweep`` ≡ sweep-path helper on flat ψ.

    Builds matching inputs for both paths.  The apply path consumes
    a per-ordinate ψ field of shape ``(ng, M, nx)`` and folds to the
    L=0 moment ``φ_0 = Σw · ψ``; the sweep path consumes ``Q̄ = 0.5·
    Σ_t · φ_0`` directly (the within-group source the SI loop
    builds).  For identical underlying flat-ψ scalar flux, both
    paths MUST produce bit-identical (up to FP-non-associativity)
    seed values.
    """
    sigma_t_gx, dr, bc_outer = _build_inputs(
        nx, sigma_t_value, bc_outer_value,
    )

    # Apply path: ψ_n = ψ_const for all (n, i) ⇒ φ_0 = Σw · ψ_const.
    # Use 2-ordinate quadrature with weights summing to 2 (canonical
    # Hébert (3.432) convention for L=0 isotropic).
    M = 2
    weights = np.array([1.0, 1.0])  # Σw = 2
    mu = np.array([-1.0, 1.0])
    psi_level = np.full((1, M, nx), psi_const)
    ctx = CarlsonSweepContext(
        sigma_t=sigma_t_gx, dr=dr,
        mu_quad=mu, weights=weights,
        bc_outer_value=bc_outer,
    )
    seed_apply = CarlsonInwardSweep()(psi_level, ctx)

    # Sweep path: feeds Q̄ = 0.5 · Σ_t · φ_0 = 0.5 · Σ_t · 2 · ψ_const
    #            = Σ_t · ψ_const directly.
    Q_bar = np.full((1, nx), sigma_t_value * psi_const)
    seed_sweep = carlson_inward_sweep_from_source(
        Q_bar=Q_bar,
        sigma_t=sigma_t_gx,
        dr=dr,
        bc_outer_value=bc_outer,
    )

    np.testing.assert_allclose(
        seed_sweep, seed_apply,
        rtol=1e-14, atol=1e-14,
        err_msg=(
            f"Apply-vs-sweep Carlson seed drift "
            f"(nx={nx}, Σ_t={sigma_t_value}, bc_outer="
            f"{bc_outer_value}, ψ_const={psi_const})"
        ),
    )


@pytest.mark.foundation
@pytest.mark.sentinel
@pytest.mark.catches("ERR-026")
def test_carlson_seed_helper_is_linear_in_Q_bar():
    r"""``carlson_inward_sweep_from_source`` is linear in ``Q_bar``.

    Required for the operator to remain linear: the seed propagates
    through the M-M angular recurrence as the half-angle face flux
    initial condition.  A non-linear seed would break
    :meth:`InvertibleOperator.apply_transpose` / dense matrix probing.
    """
    nx = 8
    sigma_t_gx = np.full((1, nx), 0.5)
    dr = np.full(nx, 0.125)
    bc_outer = np.full((1,), 0.0)

    rng = np.random.default_rng(seed=20260512)
    Q_a = rng.standard_normal((1, nx))
    Q_b = rng.standard_normal((1, nx))
    alpha, beta = 2.3, -0.7

    seed_a = carlson_inward_sweep_from_source(Q_a, sigma_t_gx, dr, bc_outer)
    seed_b = carlson_inward_sweep_from_source(Q_b, sigma_t_gx, dr, bc_outer)
    seed_lin = carlson_inward_sweep_from_source(
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

    seed_a = carlson_inward_sweep_from_source(Q_bar, sigma_t_gx, dr, bc_a)
    seed_b = carlson_inward_sweep_from_source(Q_bar, sigma_t_gx, dr, bc_b)
    seed_lin = carlson_inward_sweep_from_source(
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
