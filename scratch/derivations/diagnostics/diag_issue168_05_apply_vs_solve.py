"""Diagnostic: are apply and solve the same discrete operator?

If apply truly inverts solve, then apply(solve(q)) == q to round-off.
If they implement different discrete equations, the residual is O(1).

Tests on a slab (Cartesian) and on a sphere with simple Q.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import build_spherical_mms_case
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import SNStreamingOperator
from orpheus.sn.scattering import ScatteringOperator
from orpheus.sn.fission import FissionOperator
from orpheus.sn.solver import SNSolver


def test_apply_solve_consistency_sphere():
    """For solve(q): produces ψ_solve. Apply(ψ_solve) should give q if
    apply == L = solve^{-1}. Quantify the residual."""
    case = build_spherical_mms_case(n_ordinates=4)
    nc = 10
    mesh = case.build_mesh(nc)
    Q = case.external_source(mesh)

    # Build the SNStreamingOperator directly via SNSolver
    sn_mesh = SNMesh(mesh, case.quadrature)
    solver = SNSolver(case.materials, sn_mesh)
    L = solver.L
    eq_map = L._ensure_eq_map()
    n = eq_map.n_unknowns

    # 1. Build a packed RHS analogous to what _solve_krylov uses
    from orpheus.sn.solver import _build_rhs_spherical
    sum_w = float(case.quadrature.weights.sum())
    fission_src_norm = np.zeros((mesh.nx, 1, 1))  # no fission for fixed-source MMS
    phi_init = np.zeros((mesh.nx, 1, 1))
    rhs = _build_rhs_spherical(
        fission_src_norm, phi_init, eq_map, case.quadrature,
        solver.sig_s, solver.sig2, mesh.mat_map,
        mesh.nx, 1,
    )
    # Add the external source contribution (Q is already in the structured form)
    # Actually rhs is computed from scattering+fission only, plus the per-cell source.
    # For fixed-source MMS the exterior source is the dominant RHS — let me just
    # use the krylov solver directly.

    # Solve via krylov to get a converged solution
    result_krylov = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q,
        inner_solver="krylov", max_inner=100, inner_tol=1e-12,
    )
    phi_krylov = result_krylov.scalar_flux[:, 0, 0]

    # Solve via source iteration to get the WDD sweep solution
    result_si = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q,
        inner_solver="source_iteration", max_inner=100, inner_tol=1e-12,
    )
    phi_si = result_si.scalar_flux[:, 0, 0]

    print()
    print(f"=== Apply vs solve consistency on sphere MMS, nc={nc} ===")
    phi_ref = case.phi_exact(mesh.centers)
    print(f"  phi_ref      : {phi_ref}")
    print(f"  phi_krylov   : {phi_krylov}")
    print(f"  phi_si (WDD) : {phi_si}")
    print()
    print(f"  ||krylov - WDD||_inf = {np.max(np.abs(phi_krylov - phi_si)):.4e}")
    print(f"  ||krylov - exact||_inf = {np.max(np.abs(phi_krylov - phi_ref)):.4e}")
    print(f"  ||WDD    - exact||_inf = {np.max(np.abs(phi_si    - phi_ref)):.4e}")
    print()


def test_outer_cell_balance_check():
    """Check the discrete balance equation at the OUTER cell explicitly.
    For exact phi_exact and using the symmetric-closure operator, what is
    the residual at cell nx-1?"""
    case = build_spherical_mms_case(n_ordinates=4)
    print()
    print(f"=== Truncation error vs h (sphere, isotropic MMS) ===")
    print(f"  nc      total     inner       outer       interior")
    last = None
    for nc in (10, 20, 40, 80):
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)

        sn_mesh = SNMesh(mesh, case.quadrature)
        solver = SNSolver(case.materials, sn_mesh)
        L = solver.L
        eq_map = L._ensure_eq_map()

        quad = case.quadrature
        W = quad.weights.sum()
        psi_exact_packed = np.zeros((1, eq_map.n_eq))
        A_at_centers = case.phi_exact(mesh.centers)
        for k in range(eq_map.n_eq):
            n = eq_map.ordinate[k]
            i = eq_map.ix[k]
            psi_exact_packed[0, k] = A_at_centers[i] / W
        psi_exact_packed = psi_exact_packed.ravel(order='F')

        L_psi = L.apply(psi_exact_packed)

        Q_packed = np.zeros((1, eq_map.n_eq))
        for k in range(eq_map.n_eq):
            n = eq_map.ordinate[k]
            i = eq_map.ix[k]
            Q_packed[0, k] = Q[n, i, 0, 0]
        Q_packed = Q_packed.ravel(order='F')

        sigma_s = case.sigma_s
        A_phi = A_at_centers
        residual_packed = L_psi.copy()
        for k in range(eq_map.n_eq):
            n = eq_map.ordinate[k]
            i = eq_map.ix[k]
            residual_packed[k] -= (sigma_s / W) * A_phi[i]
            residual_packed[k] -= Q_packed[k] / W

        # Separate inner-cell, outer-cell, interior residuals
        i0_mask = np.zeros(eq_map.n_eq, dtype=bool)
        ifin_mask = np.zeros(eq_map.n_eq, dtype=bool)
        for k in range(eq_map.n_eq):
            i0_mask[k] = (eq_map.ix[k] == 0)
            ifin_mask[k] = (eq_map.ix[k] == nc - 1)
        interior_mask = ~(i0_mask | ifin_mask)

        r_inf = np.max(np.abs(residual_packed))
        r_inner = np.max(np.abs(residual_packed[i0_mask])) if i0_mask.any() else 0
        r_outer = np.max(np.abs(residual_packed[ifin_mask])) if ifin_mask.any() else 0
        r_interior = np.max(np.abs(residual_packed[interior_mask])) if interior_mask.any() else 0
        # Find which cell has the max interior residual
        interior_idx = np.where(interior_mask)[0]
        if len(interior_idx) > 0:
            k_max = interior_idx[np.argmax(np.abs(residual_packed[interior_mask]))]
            i_max = eq_map.ix[k_max]
        else:
            i_max = -1
        print(f"  {nc:3d}  inf={r_inf:.4e}  i=0:{r_inner:.4e}  i=N-1:{r_outer:.4e}  interior:{r_interior:.4e} (at i={i_max})")
        last = r_inf

    # Print residual per cell, per ordinate
    # (per-cell breakdown removed; only summary stats remain above)
