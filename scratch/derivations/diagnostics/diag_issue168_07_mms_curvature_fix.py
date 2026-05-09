"""Diagnostic: test the hypothesis that the MMS source is missing the
spherical-curvature term (2μ/r)·A.

For an isotropic trial ψ_n(r) = A_phi(r)/W on a sphere, the continuous
spherical SN streaming term is:

   μ ∂ψ/∂r + (2μ/r) ψ

The slab streaming is just μ ∂ψ/∂r. The current MMS source only
includes the slab streaming, which mismatches the actual spherical
operator by the geometric (2μ/r) ψ term.

Patch: add the (2μ/r) ψ term to the MMS source and rerun.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import build_spherical_mms_case
from orpheus.sn import solve_sn_fixed_source
from orpheus.sn.geometry import SNMesh
from orpheus.sn.solver import SNSolver


def make_correct_spherical_mms_source(case, mesh):
    """Compute the spherical-correct MMS source.

    Continuous spherical 1D SN with isotropic ψ_n(r) = A_phi(r)/W:
       μ ∂ψ/∂r + (2μ/r) ψ + Σ_t ψ - (Σ_s/W) φ = q_n^ext / W
    where φ = A_phi.
    Solving for q_n^ext:
       q_n^ext = μ A_phi' + (2μ/r) A_phi + (Σ_t - Σ_s) A_phi
    """
    r = mesh.centers
    A = case.phi_exact(r)
    Ap = case.dphi_exact(r)
    mu = case.quadrature.mu_x

    streaming = mu[:, None] * Ap[None, :]
    curvature = 2.0 * mu[:, None] * A[None, :] / r[None, :]   # NEW term
    removal = (case.sigma_t - case.sigma_s) * A[None, :]
    Q = streaming + curvature + removal
    return Q[:, :, None, None]


def _l2_1d(phi_num, phi_ref, volumes):
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


def test_spherical_mms_with_curvature_correction():
    """If the MMS is missing the (2μ/r)·A term, adding it should make
    Krylov-on-apply converge at O(h^2) without any operator changes.
    """
    case = build_spherical_mms_case(n_ordinates=4)
    n_cells = [10, 20]

    print()
    print("=== Spherical MMS with curvature correction (2μ/r)·A in source ===")
    errs = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q_corrected = make_correct_spherical_mms_source(case, mesh)

        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q_corrected,
            inner_solver="krylov", max_inner=200, inner_tol=1e-9,
        )
        phi_num = result.scalar_flux[:, 0, 0]
        phi_ref = case.phi_exact(mesh.centers)
        err = _l2_1d(phi_num, phi_ref, mesh.volumes)
        max_err = np.max(np.abs(phi_num - phi_ref))
        errs.append(err)
        print(f"  nc={nc:3d}  L2 err = {err:.4e}  max err = {max_err:.4e}")

    if len(errs) > 1:
        orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
        print(f"  Orders: {orders}")


def test_residual_with_curvature_correction():
    """Apply L to the exact ψ, with the corrected MMS source. Residual
    should now decay as O(h²) — confirming the hypothesis."""
    case = build_spherical_mms_case(n_ordinates=4)
    print()
    print(f"=== Truncation residual with curvature-corrected MMS source ===")
    print(f"  nc   ||resid||_inf   ratio")
    last = None
    for nc in (10, 20, 40, 80):
        mesh = case.build_mesh(nc)
        Q = make_correct_spherical_mms_source(case, mesh)
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

        residual_packed = L_psi.copy()
        for k in range(eq_map.n_eq):
            n = eq_map.ordinate[k]
            i = eq_map.ix[k]
            residual_packed[k] -= (case.sigma_s / W) * A_at_centers[i]
            residual_packed[k] -= Q_packed[k]

        r_inf = np.max(np.abs(residual_packed))
        ratio = last / r_inf if last is not None else float("nan")
        print(f"  {nc:3d}  {r_inf:.4e}    {ratio:.3f}")
        last = r_inf
