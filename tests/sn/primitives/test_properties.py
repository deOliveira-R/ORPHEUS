"""Unit tests for SN 1D solver properties.

Tests structural properties that the SN solution must satisfy,
independent of the reference eigenvalue:
- Gauss-Legendre quadrature weights sum to 2
- Flux symmetry for symmetric geometry with reflective BCs
- Particle balance: production / absorption = keff (no leakage with reflective BCs)
"""

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.geometry import (
    BC,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.solver import solve_sn


def _homogeneous_slab_mesh(n_cells: int, total_width: float, mat_id: int = 0) -> Mesh1D:
    """Single-region Cartesian mesh helper."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(Region(mat_id=mat_id, outer_thickness_cm=total_width),),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=n_cells),))


def _slab_fuel_moderator_mesh(
    n_fuel: int, n_mod: int, t_fuel: float, t_mod: float,
) -> Mesh1D:
    """Two-region fuel-moderator slab mesh."""
    geom = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=2, outer_thickness_cm=t_fuel),
            Region(mat_id=0, outer_thickness_cm=t_mod),
        ),
        bcs=(BC.reflective, BC.reflective),
    )
    return Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=n_fuel),
        RegionMesh(n_cells=n_mod),
    ))

pytestmark = pytest.mark.l0  # SN property checks (quadrature weights, symmetry, balance)


def test_gl_weights_sum():
    """Gauss-Legendre weights on [-1,1] must sum to 2."""
    for N in [4, 8, 16, 32]:
        quad = Quadrature.gauss_legendre(N)
        np.testing.assert_allclose(
            quad.weights.sum(), 2.0, atol=1e-14,
            err_msg=f"GL({N}) weights sum to {quad.weights.sum()}, expected 2.0",
        )


def test_gl_symmetry():
    """GL quadrature points must be symmetric: μ[i] = -μ[N-1-i]."""
    quad = Quadrature.gauss_legendre(16)
    np.testing.assert_allclose(
        quad.mu_x, -quad.mu_x[::-1], atol=1e-14,
    )


def test_flux_symmetry():
    """For symmetric geometry, the scalar flux must be symmetric about the center."""
    case = get("sn_slab_1eg_1rg")
    mix = next(iter(case.materials.values()))

    # Build a symmetric 2-region slab: fuel | moderator | moderator | fuel
    from orpheus.derivations.common.xs_library import get_mixture
    fuel = get_mixture("A", "1g")
    mod = get_mixture("B", "1g")
    materials = {2: fuel, 0: mod}

    # Symmetric layout: 10 fuel | 10 mod (half-cell with reflective BCs)
    mesh = _slab_fuel_moderator_mesh(
        n_fuel=10, n_mod=10, t_fuel=0.5, t_mod=0.5,
    )
    quad = Quadrature.gauss_legendre(8)
    result = solve_sn(materials, mesh, quad, max_outer=200,
                      max_inner=500, inner_tol=1e-10)

    # With reflective BCs at both ends, a half-cell geometry is symmetric
    # about its midpoint only if the materials are arranged symmetrically.
    # Here fuel|mod is NOT symmetric about the center, but the flux
    # should still be smooth and monotonic from fuel to moderator.
    # A stronger test: a homogeneous slab must have exactly flat flux.
    mesh_homo = _homogeneous_slab_mesh(20, 2.0, mat_id=0)
    result_homo = solve_sn({0: mix}, mesh_homo, quad, max_outer=200,
                           max_inner=500, inner_tol=1e-10)
    # PR-INDEX-5: scalar_flux principled (ng, nx, ny) — group-0 radial slice.
    flux = result_homo.scalar_flux.values[0, :]  # (nx,) for group 0
    np.testing.assert_allclose(
        flux, flux[0], rtol=1e-6,
        err_msg="Homogeneous slab flux is not flat",
    )


def test_particle_balance():
    """For reflective BCs (no leakage), production / absorption = keff."""
    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    mesh = _homogeneous_slab_mesh(20, 2.0, mat_id=0)
    quad = Quadrature.gauss_legendre(8)
    result = solve_sn(materials, mesh, quad,
                      max_inner=500, inner_tol=1e-10)

    # Volume-weighted production and absorption rates.
    # PR-INDEX-5: scalar_flux principled (ng, nx, ny=1) → (ng, nx).
    dx = mesh.widths
    flux = result.scalar_flux.values  # (ng, nx)
    sig_p = mix.SigP
    sig_a = mix.SigC + mix.SigF

    production = np.sum(flux * sig_p[:, None] * dx[None, :])
    absorption = np.sum(flux * sig_a[:, None] * dx[None, :])

    k_balance = production / absorption
    np.testing.assert_allclose(
        k_balance, result.keff, rtol=1e-6,
        err_msg=f"Particle balance: prod/abs={k_balance:.8f} ≠ keff={result.keff:.8f}",
    )
