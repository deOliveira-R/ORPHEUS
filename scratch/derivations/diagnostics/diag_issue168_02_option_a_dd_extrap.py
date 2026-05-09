"""Diagnostic: Option A — DD diamond extrapolation at outer face.

Monkey-patches transport_operator_matvec_spherical and re-runs the MMS
to measure the convergence-order improvement.

The DD relation states:
    psi_cell = 0.5 * (psi_face_in + psi_face_out)
=>  psi_face_out = 2 * psi_cell - psi_face_in

At the outer cell i = nx-1 for outgoing mu>0:
- psi_face_in (the LEFT face of the outer cell) is computed from
  arithmetic average of cell centers as everywhere else:
    psi_face_in = 0.5 * (fi[:, n, nx-2, 0] + fi[:, n, nx-1, 0])
- psi_face_out is then:
    psi_face_out = 2 * fi[:, n, nx-1, 0] - psi_face_in
                 = 2*fi[:,nx-1] - 0.5*(fi[:,nx-2] + fi[:,nx-1])
                 = 1.5*fi[:,nx-1] - 0.5*fi[:,nx-2]
This is a one-sided second-order extrapolation through the last two
cell centers.

Created by numerics-investigator on 2026-05-09.
"""
from __future__ import annotations

import numpy as np
import pytest

import orpheus.sn.operator as opmod
from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case,
    build_cylindrical_mms_case,
)
from orpheus.sn import solve_sn_fixed_source


def _patched_spherical_matvec(
    solution, eq_map, quad, sig_t,
    nx, ng,
    face_areas, volumes,
    alpha_half, redist_dAw, tau_mm,
    *,
    bc_outer=None,
):
    """Identical to transport_operator_matvec_spherical EXCEPT the outer-face
    cell-center substitution is replaced with DD diamond extrapolation."""
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
                # === OPTION A: DD diamond extrapolation ===
                if i > 0:
                    psi_face_in = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
                else:
                    psi_face_in = 0.0
                psi_right = 2.0 * fi[:, n, i, 0] - psi_face_in
                # === end Option A ===
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


def _patched_cylindrical_matvec(
    solution, eq_map, quad, sig_t,
    nx, ng,
    face_areas, volumes,
    alpha_per_level, redist_dAw_per_level, tau_mm_per_level,
    *,
    bc_outer=None,
):
    """Same DD diamond extrapolation patch for the cylindrical operator."""
    fi = opmod.solution_to_angular_flux_cylindrical(
        solution, eq_map, quad, nx, ng, bc_outer=bc_outer
    )
    ref_x = quad.reflection_index("x")
    A = face_areas
    V = volumes[:, 0]
    N = quad.N
    mu = quad.mu_x

    ord_to_level = np.empty(N, dtype=int)
    ord_to_local = np.empty(N, dtype=int)
    for p, level_idx in enumerate(quad.level_indices):
        for m_local, n in enumerate(level_idx):
            ord_to_level[n] = p
            ord_to_local[n] = m_local

    lhs = np.empty((ng, eq_map.n_eq))
    for k in range(eq_map.n_eq):
        n = eq_map.ordinate[k]
        i = eq_map.ix[k]
        psi_ni = fi[:, n, i, 0]

        p = ord_to_level[n]
        m_local = ord_to_local[n]
        alpha = alpha_per_level[p]
        dAw = redist_dAw_per_level[p]
        tau_level = tau_mm_per_level[p]
        level_idx = quad.level_indices[p]
        M = len(level_idx)

        if i < nx - 1:
            psi_right = 0.5 * (fi[:, n, i, 0] + fi[:, n, i + 1, 0])
        else:
            if mu[n] > 1e-15:
                if i > 0:
                    psi_face_in = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
                else:
                    psi_face_in = 0.0
                psi_right = 2.0 * fi[:, n, i, 0] - psi_face_in
            else:
                psi_right = fi[:, ref_x[n], i, 0]

        if i > 0:
            psi_left = 0.5 * (fi[:, n, i - 1, 0] + fi[:, n, i, 0])
        else:
            psi_left = 0.0

        streaming = mu[n] * (A[i + 1] * psi_right - A[i] * psi_left) / V[i]

        dA_w = dAw[i, m_local]
        tau_m = tau_level[m_local]

        if m_local < M - 1:
            n_next = level_idx[m_local + 1]
            psi_angle_right = tau_m * fi[:, n_next, i, 0] + (1.0 - tau_m) * fi[:, n, i, 0]
        else:
            psi_angle_right = fi[:, n, i, 0]

        if m_local > 0:
            n_prev = level_idx[m_local - 1]
            tau_prev = tau_level[m_local - 1]
            psi_angle_left = tau_prev * fi[:, n, i, 0] + (1.0 - tau_prev) * fi[:, n_prev, i, 0]
        else:
            psi_angle_left = fi[:, n, i, 0]

        redistribution = dA_w * (alpha[m_local + 1] * psi_angle_right
                                 - alpha[m_local] * psi_angle_left) / V[i]

        collision = sig_t[i, 0, :] * psi_ni
        lhs[:, k] = streaming + redistribution + collision

    return lhs.ravel(order='F')


def _l2_1d(phi_num, phi_ref, volumes):
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


def _convergence_table_with_patch(case, n_cells, patched_matvec, attr_name):
    """Run convergence with monkey-patched matvec."""
    original = getattr(opmod, attr_name)
    setattr(opmod, attr_name, patched_matvec)
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
        errs = np.asarray(errs)
        orders = np.log2(errs[:-1] / errs[1:])
        return errs, orders
    finally:
        setattr(opmod, attr_name, original)


def test_spherical_option_a_dd_extrapolation():
    """With DD diamond extrapolation at the outer face,
    spherical MMS should converge at >= 1.9 (= O(h^2))."""
    case = build_spherical_mms_case(n_ordinates=4)
    n_cells = [20, 40, 80]
    print()
    print("=== Option A: DD diamond extrapolation at outer face (sphere) ===")
    errs, orders = _convergence_table_with_patch(
        case, n_cells, _patched_spherical_matvec,
        "transport_operator_matvec_spherical",
    )
    for nc, e in zip(n_cells, errs):
        print(f"    nc={nc:3d}  L2 err = {e:.4e}")
    print(f"    Orders (n -> 2n): {orders}")
    print()


def test_cylindrical_option_a_dd_extrapolation():
    """Same Option A patch for cylindrical."""
    case = build_cylindrical_mms_case()
    n_cells = [10, 20]
    print()
    print("=== Option A: DD diamond extrapolation at outer face (cyl) ===")
    errs, orders = _convergence_table_with_patch(
        case, n_cells, _patched_cylindrical_matvec,
        "transport_operator_matvec_cylindrical",
    )
    for nc, e in zip(n_cells, errs):
        print(f"    nc={nc:3d}  L2 err = {e:.4e}")
    print(f"    Orders: {orders}")
    print()
