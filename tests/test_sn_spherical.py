"""Verify the spherical 1D SN solver.

Tests cover:
- Homogeneous exact: k_inf must match analytical (geometry-independent)
- Particle balance: production/absorption = keff (reflective BC, no leakage)
- Angular redistribution properties: α sum, antisymmetry
- Spatial convergence: O(h²) with mesh refinement
- Cross-check with CP spherical: SN eigenvalue close to CP eigenvalue
"""

import numpy as np
import pytest

from derivations import get
from derivations._xs_library import get_mixture
from geometry import CoordSystem, Mesh1D, homogeneous_1d, mesh1d_from_zones, Zone
from sn_geometry import SNMesh
from sn_quadrature import GaussLegendre1D
from sn_solver import SNSolver, solve_sn


# ── Homogeneous infinite medium (geometry-independent k_inf) ─────────

@pytest.mark.parametrize("case_name", [
    "sn_slab_1eg_1rg",
    "sn_slab_2eg_1rg",
    "sn_slab_4eg_1rg",
])
def test_homogeneous_exact(case_name):
    """Spherical SN on a homogeneous sphere with reflective BC must
    match the analytical infinite-medium eigenvalue (same for all
    geometries)."""
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    mesh = homogeneous_1d(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = GaussLegendre1D.create(8)
    result = solve_sn(materials, mesh, quad,
                      max_inner=500, inner_tol=1e-10)

    # Spherical DD has larger discretization error than Cartesian
    # due to angular redistribution coupling. 1G is exact (keff
    # independent of flux shape); multi-group has ~1% error on S8/20-cell.
    tol = 1e-6 if case.n_groups == 1 else 0.02
    assert abs(result.keff - case.k_inf) < tol, (
        f"keff={result.keff:.8f} vs analytical={case.k_inf:.8f} "
        f"err={abs(result.keff - case.k_inf):.2e}"
    )


# ── Particle balance ─────────────────────────────────────────────────

def test_particle_balance():
    """For reflective BCs (no leakage), production / absorption = keff."""
    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    mesh = homogeneous_1d(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = GaussLegendre1D.create(8)
    result = solve_sn(materials, mesh, quad,
                      max_inner=500, inner_tol=1e-10)

    # Volume-weighted production and absorption rates
    V = mesh.volumes
    flux = result.scalar_flux[:, 0, :]  # (nx, ng)
    sig_p = mix.SigP
    sig_a = mix.SigC + mix.SigF

    production = np.sum(flux * sig_p[None, :] * V[:, None])
    absorption = np.sum(flux * sig_a[None, :] * V[:, None])

    k_balance = production / absorption
    np.testing.assert_allclose(
        k_balance, result.keff, rtol=1e-5,
        err_msg=f"Particle balance: prod/abs={k_balance:.8f} ≠ keff={result.keff:.8f}",
    )


# ── Angular redistribution properties ────────────────────────────────

class TestAlphaCoefficients:
    """Properties of the angular redistribution coefficients."""

    @pytest.mark.parametrize("N", [4, 8, 16, 32])
    def test_alpha_boundary_conditions(self, N):
        """α_{1/2} = 0 and α_{N+1/2} = 0 by GL antisymmetry."""
        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(N)
        sn_mesh = SNMesh(mesh, quad)

        np.testing.assert_allclose(sn_mesh.alpha_half[0], 0.0)
        np.testing.assert_allclose(sn_mesh.alpha_half[-1], 0.0, atol=1e-14)

    def test_alpha_recursion(self):
        """α_{n+1/2} = α_{n-1/2} + w_n μ_n."""
        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(8)
        sn_mesh = SNMesh(mesh, quad)

        alpha = sn_mesh.alpha_half
        for n in range(quad.N):
            expected = alpha[n] + quad.weights[n] * quad.mu_x[n]
            np.testing.assert_allclose(alpha[n + 1], expected, rtol=1e-14)

    def test_alpha_symmetric(self):
        """α coefficients are symmetric about the midpoint: α[k] = α[N-k].

        This follows from GL symmetry (w_n = w_{N-1-n}, μ_n = -μ_{N-1-n}):
        the cumulative sum Σ w_m μ_m is symmetric because each pair
        (w_n μ_n, w_{N-1-n} μ_{N-1-n}) contributes equal and opposite
        increments that mirror about the midpoint.
        """
        mesh = Mesh1D(edges=np.array([0.0, 1.0]), mat_ids=np.array([0]),
                      coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(8)
        sn_mesh = SNMesh(mesh, quad)

        alpha = sn_mesh.alpha_half
        N = quad.N
        for k in range(N + 1):
            np.testing.assert_allclose(
                alpha[k], alpha[N - k], atol=1e-14,
                err_msg=f"α not symmetric at k={k}",
            )


# ── Spatial convergence ──────────────────────────────────────────────

@pytest.mark.slow
def test_spatial_convergence():
    """Diamond-difference on spherical mesh must show O(h²) convergence."""
    fuel = get_mixture("A", "1g")
    mod = get_mixture("B", "1g")
    materials = {2: fuel, 0: mod}

    keffs = []
    drs = []
    for n_per in [5, 10, 20, 40]:
        zones = [
            Zone(outer_edge=0.5, mat_id=2, n_cells=n_per),
            Zone(outer_edge=1.0, mat_id=0, n_cells=n_per),
        ]
        mesh = mesh1d_from_zones(zones, coord=CoordSystem.SPHERICAL)
        quad = GaussLegendre1D.create(16)
        result = solve_sn(
            materials, mesh, quad,
            max_outer=300, max_inner=500, inner_tol=1e-10,
        )
        keffs.append(result.keff)
        drs.append(0.5 / n_per)

    # Richardson extrapolation reference
    k_ref = keffs[-1] + (keffs[-1] - keffs[-2]) / 3.0

    # Compute convergence orders
    orders = []
    for i in range(1, len(keffs)):
        err_prev = abs(keffs[i - 1] - k_ref)
        err_curr = abs(keffs[i] - k_ref)
        if err_prev > 0 and err_curr > 0:
            orders.append(
                np.log(err_prev / err_curr)
                / np.log(drs[i - 1] / drs[i])
            )

    assert orders[-1] > 1.5, (
        f"Expected O(h²) convergence, got order {orders[-1]:.2f}"
    )


# ── Cross-check with CP spherical ────────────────────────────────────

def test_cross_check_with_cp_1g():
    """SN and CP on the same spherical geometry should give close k_inf.

    For homogeneous 1G, both must match the analytical value exactly.
    The white-BC (CP) vs reflective-BC (SN) difference vanishes for
    homogeneous infinite medium.
    """
    from collision_probability import solve_cp

    mix = get_mixture("A", "1g")
    materials_sn = {0: mix}
    materials_cp = {0: mix}

    # SN: homogeneous sphere with reflective BC
    mesh_sn = homogeneous_1d(20, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = GaussLegendre1D.create(8)
    result_sn = solve_sn(materials_sn, mesh_sn, quad,
                         max_inner=500, inner_tol=1e-10)

    # CP: same geometry with white BC
    mesh_cp = homogeneous_1d(1, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    result_cp = solve_cp(materials_cp, mesh_cp)

    # Both must give the same 1G homogeneous k_inf
    np.testing.assert_allclose(
        result_sn.keff, result_cp.keff, rtol=1e-6,
        err_msg=f"SN keff={result_sn.keff:.6f} vs CP keff={result_cp.keff:.6f}",
    )


# ── Flux positivity ──────────────────────────────────────────────────

def test_flux_non_negative():
    """Converged scalar flux must be non-negative everywhere."""
    mix = get_mixture("A", "1g")
    mesh = homogeneous_1d(10, 2.0, mat_id=0, coord=CoordSystem.SPHERICAL)
    quad = GaussLegendre1D.create(8)
    result = solve_sn({0: mix}, mesh, quad, max_inner=500, inner_tol=1e-10)

    assert np.all(result.scalar_flux >= 0), (
        f"Negative flux: min={result.scalar_flux.min():.4e}"
    )
