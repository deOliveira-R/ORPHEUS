"""Diagnostic: Option A' — read the OUTGOING face from the BC.

The cleaner statement: at the outer boundary, the face flux for outgoing
mu>0 is the value the BC operator returns for incoming directions, on the
ORDINATE-MIRRORED ordinate. For vacuum BC: the outgoing face flux is the
value extrapolated FROM the cell center to the face. For reflective BC:
the outgoing flux is symmetric. For white BC: the flux is the angle-
weighted incoming average, etc.

Actually, on closer inspection, the OUTER face for MU>0 ordinates is the
OUTGOING face — NOT a BC-controlled face. The BC controls only what flows
INWARD. Outgoing flux at r=R is simply ψ(R, μ>0), which is what the
operator must approximate.

For the MMS with vacuum BC at r=R: phi_exact(R) = sin(pi) = 0, so the
EXACT psi_face_out = 0 too (isotropic ansatz). Try setting psi_right = 0
explicitly and see if the order improves.

This isn't a "fix" — it's a probe to see if the boundary value alone is
the bottleneck.
"""
from __future__ import annotations

import numpy as np
import pytest

import orpheus.sn.operator as opmod
from orpheus.derivations.continuous.mms.sn import build_spherical_mms_case
from orpheus.sn import solve_sn_fixed_source


def _patched_spherical_matvec_zero_outer(
    solution, eq_map, quad, sig_t,
    nx, ng,
    face_areas, volumes,
    alpha_half, redist_dAw, tau_mm,
    *,
    bc_outer=None,
):
    """psi_right = 0 explicitly at outer face for outgoing mu>0."""
    fi = opmod.solution_to_angular_flux_spherical(
        solution, eq_map, quad, nx, ng, bc_outer=bc_outer
    )
    ref_x = quad.reflection_index("x")
    A = face_areas
    V = volumes[:, 0]
    dAw = redist_dAw
    alpha = alpha_half
    N = quad.N
    mu = quad.mu_x

    lhs = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        psi_ni = fi[:, n, i, 0]

        if i < nx - 1:
            psi_right = 0.5 * (fi[:, n, i, 0] + fi[:, n, i + 1, 0])
        else:
            if mu[n] > 1e-15:
                psi_right = np.zeros(ng)  # FORCED ZERO for vacuum BC
            else:
                psi_right = fi[:, ref_x[n], i, 0]

        if i > 0:
            psi_left = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
        else:
            psi_left = 0.0

        streaming = mu[n] * (A[i + 1] * psi_right - A[i] * psi_left) / V[i]
        dA_w = dAw[i, n]
        tau_n = tau_mm[n]

        if n < N - 1:
            psi_angle_right = tau_n * fi[:, n + 1, i, 0] + (1.0 - tau_n) * fi[:, n, i, 0]
        else:
            psi_angle_right = fi[:, n, i, 0]
        if n > 0:
            psi_angle_left = tau_mm[n - 1] * fi[:, n, i, 0] + (1.0 - tau_mm[n - 1]) * fi[:, n - 1, i, 0]
        else:
            psi_angle_left = fi[:, n, i, 0]
        redistribution = dA_w * (alpha[n + 1] * psi_angle_right
                                 - alpha[n] * psi_angle_left) / V[i]
        collision = sig_t[i, 0, :] * psi_ni
        lhs[:, k] = streaming + redistribution + collision

    return lhs.ravel(order='F')


def _l2_1d(phi_num, phi_ref, volumes):
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


def test_spherical_zero_outer_face():
    """If outer face = 0 explicitly, what is the order?"""
    case = build_spherical_mms_case(n_ordinates=4)
    n_cells = [10, 20, 40]

    original = opmod.transport_operator_matvec_spherical
    opmod.transport_operator_matvec_spherical = _patched_spherical_matvec_zero_outer
    try:
        errs = []
        for nc in n_cells:
            mesh = case.build_mesh(nc)
            Q = case.external_source(mesh)
            result = solve_sn_fixed_source(
                case.materials, mesh, case.quadrature, Q,
                inner_solver="krylov", max_inner=200, inner_tol=1e-9,
            )
            phi_num = result.scalar_flux[:, 0, 0]
            phi_ref = case.phi_exact(mesh.centers)
            errs.append(_l2_1d(phi_num, phi_ref, mesh.volumes))
    finally:
        opmod.transport_operator_matvec_spherical = original

    errs = np.asarray(errs)
    orders = np.log2(errs[:-1] / errs[1:])
    print()
    print("=== Probe: psi_right = 0 explicitly at outer face for mu>0 ===")
    for nc, e in zip(n_cells, errs):
        print(f"    nc={nc:3d}  L2 err = {e:.4e}")
    print(f"    Orders: {orders}")
    print()
