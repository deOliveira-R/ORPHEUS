"""Diagnostic: very fine-grained residual analysis.

Goal: pin down WHY the truncation residual stays at ~0.4 — is it the
outer boundary, the inner boundary, the interior, the redistribution
term, or what?
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import build_spherical_mms_case
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver


def test_residual_per_cell_per_ordinate_n40():
    """At nc=40, list residuals per cell, per ordinate. Look for the
    structure of the error — is it spread over all cells, or focused
    at a few?"""
    case = build_spherical_mms_case(n_ordinates=4)
    nc = 80
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
    residual_packed = L_psi.copy()
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        residual_packed[k] -= (sigma_s / W) * A_at_centers[i]
        residual_packed[k] -= Q_packed[k] / W   # FIX: Q is un-divided, divide by W

    # Reshape into (n_ord, n_cells)
    res_arr = np.zeros((4, nc))
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        res_arr[n, i] = residual_packed[k]

    print()
    print(f"=== Residual structure at nc={nc} ===")
    print(f"  Max |res| per ordinate:")
    for n in range(4):
        ix_max = np.argmax(np.abs(res_arr[n]))
        print(f"    n={n}  mu={quad.mu_x[n]:+.3f}  max|res|={np.max(np.abs(res_arr[n])):.4e}  "
              f"at i={ix_max}, r={mesh.centers[ix_max]:.3f}")

    # First few and last few cells per ordinate
    print(f"\n  Residual values at outer boundary (i = {nc-3}..{nc-1}):")
    print(f"   i      r       n=0       n=1       n=2       n=3")
    for i in range(nc-3, nc):
        print(f"   {i:2d}  {mesh.centers[i]:.4f}   "
              f"{res_arr[0,i]:+.4e}  {res_arr[1,i]:+.4e}  {res_arr[2,i]:+.4e}  {res_arr[3,i]:+.4e}")

    print(f"\n  Residual values at inner boundary (i = 0..2):")
    print(f"   i      r       n=0       n=1       n=2       n=3")
    for i in range(3):
        print(f"   {i:2d}  {mesh.centers[i]:.4f}   "
              f"{res_arr[0,i]:+.4e}  {res_arr[1,i]:+.4e}  {res_arr[2,i]:+.4e}  {res_arr[3,i]:+.4e}")

    print(f"\n  Residual values at mid-domain (i = nx/2 - 2..nx/2+2):")
    print(f"   i      r       n=0       n=1       n=2       n=3")
    for i in range(nc//2 - 2, nc//2 + 3):
        print(f"   {i:2d}  {mesh.centers[i]:.4f}   "
              f"{res_arr[0,i]:+.4e}  {res_arr[1,i]:+.4e}  {res_arr[2,i]:+.4e}  {res_arr[3,i]:+.4e}")


def test_residual_at_n0_only_check_n4():
    """For each ordinate separately, check the residual at i=0 (innermost).
    The redistribution at i=0 has alpha boundary contributions; check whether
    these cause a problem."""
    pass
