"""Diagnostic: same truncation residual analysis but using transport_sweep
output (WDD asymmetric closure) to see if THAT operator's residual is O(1)
or O(h^p).

If the WDD sweep's discrete equations also have O(1) residual on the exact
ψ, then BOTH operators are inconsistent on this MMS — meaning the MMS
itself has a problem (e.g., source manufactured with one convention,
operator using another).

If the WDD sweep's residual is O(h) or O(h^2), then the symmetric closure
is the broken one.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import build_spherical_mms_case
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import SNStreamingOperator
from orpheus.sn.solver import SNSolver
from orpheus.sn.sweep import transport_sweep


def test_wdd_sweep_residual_scaling():
    """Take exact phi_exact, apply WDD sweep, compare with input source.
    The sweep computes: psi_out = sweep(Q_total) where Q_total includes scattering.
    If the sweep is consistent, sweep(Q_built_from_phi_exact) should give phi_exact.
    Quantify how close.
    """
    case = build_spherical_mms_case(n_ordinates=4)
    print()
    print(f"=== Sweep self-consistency on phi_exact (sphere MMS) ===")
    print(f"  nc   ||phi_sweep - phi_exact||_inf   ratio")
    last_err = None
    for nc in (10, 20, 40, 80):
        mesh = case.build_mesh(nc)
        Q_ext = case.external_source(mesh)
        sn_mesh = SNMesh(mesh, case.quadrature)

        # Exact scalar flux at cell centers
        phi_exact = case.phi_exact(mesh.centers)

        # Build a one-shot sweep using exact phi as scattering input.
        # In source iteration: Q_total = Q_ext + (Sigma_s/W) * phi_exact (per ordinate, isotropic source)
        sigma_s = case.sigma_s
        W = case.quadrature.weights.sum()
        # Q_total (n, nx, 1, 1) — anisotropic source slot for sweep
        Q_total = Q_ext.copy()  # already (N, nx, 1, 1) of mu·dphi/dr + (Σt-Σs)·phi
        # Add (Σs/W)·phi for the isotropic ansatz reconstruction
        for n in range(case.quadrature.N):
            Q_total[n, :, 0, 0] += (sigma_s / W) * phi_exact

        # Q_iso_for_sweep is the isotropic external source = sum_n w_n Q_total[n] / sum_n w_n? No.
        # transport_sweep takes: Q (nx,1,ng) iso + Q_aniso (n,nx,1,ng) per-ordinate
        # The isotropic part is constructed by the sweep as Q_iso/W + Q_aniso for each ordinate.
        # Path: for our purposes, set Q_iso = 0 and put everything in Q_aniso.
        nx_local = sn_mesh.nx
        Q_iso = np.zeros((nx_local, 1, 1))
        Q_aniso = Q_total

        # sig_t for sweep is per-cell
        from orpheus.data.macro_xs.cell_xs import assemble_cell_xs
        xs = assemble_cell_xs(case.materials, sn_mesh.mat_map)
        sig_t = xs.sig_t.reshape(nx_local, 1, 1)

        psi_bc = {}
        angular_flux, scalar_flux = transport_sweep(
            Q_iso, sig_t, sn_mesh, psi_bc, Q_aniso=Q_aniso,
        )
        phi_sweep = scalar_flux[:, 0, 0]

        err = np.max(np.abs(phi_sweep - phi_exact))
        ratio = last_err / err if last_err is not None else float("nan")
        print(f"  {nc:3d}  {err:.4e}                       {ratio:.3f}")
        last_err = err


def test_apply_residual_normalized():
    """Check the symmetric-closure operator's residual but NORMALIZED by ψ.
    The raw residual ~0.4 may be normal if ψ scale grows with nc; rescale
    by ||ψ_exact||_inf and check if the normalized residual decays."""
    case = build_spherical_mms_case(n_ordinates=4)
    print()
    print(f"=== Normalized truncation residual ||L psi - rhs|| / ||psi||_inf ===")
    print(f"  nc   ||resid||_inf   ||psi||_inf   normalized")
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
        residual_packed = L_psi.copy()
        for k in range(eq_map.n_eq):
            n = eq_map.ordinate[k]
            i = eq_map.ix[k]
            residual_packed[k] -= (sigma_s / W) * A_at_centers[i]
            residual_packed[k] -= Q_packed[k]

        r_inf = np.max(np.abs(residual_packed))
        psi_inf = np.max(np.abs(psi_exact_packed))
        print(f"  {nc:3d}  {r_inf:.4e}    {psi_inf:.4e}    {r_inf/psi_inf:.4e}")
