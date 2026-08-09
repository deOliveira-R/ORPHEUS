"""Symbolic per-cell algebra walkthrough for the 2-cell × 2-ordinate
minimal reproducer.

Problem:
    - n_cells = 2, GL-2 (μ_0 = -1/√3, μ_1 = +1/√3, w_0 = w_1 = 1)
    - Homogeneous mixture B 1G: Σ_t = 2, Σ_s = 1.9 (c_self = 0.95)
    - Sphere R = 2 cm; reflective outer BC
    - Isotropic source Q = 1 → analytical streaming-equilibrium φ = 10
      everywhere → ψ_iso = 5 for every ordinate at every cell

Empirical baselines (refl, 2x2, mixture B):
    ψ_SI:  [[4.65, 6.06], [3.12, 6.16]]   off by up to 2.0
    ψ_K:   [[5.00, 5.00], [5.00, 5.00]]   exact

This script:

1. Imports the production sweep + apply-matvec.
2. Runs each to fixed point at scattering convergence.
3. Symbolically substitutes ψ_SI* into BOTH the SI per-cell equation
   and the apply-matvec per-cell residual; confirms SI satisfies its
   own equation to ~1e-14 and FAILS to satisfy the apply-matvec
   residual.
4. Same for ψ_K* — should satisfy apply-matvec residual to ~1e-14
   and fail SI's.
5. Algebraically subtracts the two equations at fixed point and
   identifies the single differing term.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.solver import solve_sn_fixed_source
from orpheus.sn.sweep.psi_half_angle_seed import (
    CarlsonSweepContext, carlson_inward_sweep_from_source,
)


def build_problem():
    fuel = get_mixture("B", "1g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=2),))
    quad = GaussLegendre1D.create(n_ordinates=2)
    sn_mesh = SNMesh(mesh, quad)
    return fuel, mesh, quad, sn_mesh


def get_fixed_points():
    """Run both solvers; return ψ_SI*, ψ_K*."""
    fuel, mesh, quad, sn_mesh = build_problem()
    Q = np.ones((quad.N, 2, 1, 1))  # external isotropic source
    res_si = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        external_source=Q, boundary_condition=BC.reflective,
        inner_solver="source_iteration",
        max_inner=10000, inner_tol=1e-14,
    )
    res_kr = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad,
        external_source=Q, boundary_condition=BC.reflective,
        inner_solver="krylov",
        max_inner=2000, inner_tol=1e-14,
    )
    psi_si = res_si.angular_flux[:, :, 0, 0]  # (N, nx)
    psi_kr = res_kr.angular_flux[:, :, 0, 0]
    return psi_si, psi_kr, sn_mesh, quad, fuel


# ─────────────────────────────────────────────────────────────────────
# Symbolic derivation: per-cell equation at the fixed point
# ─────────────────────────────────────────────────────────────────────

# Coefficients (all values from sn_mesh.reduced for this 2x2 case):
#   A     = [0, 31.665, 50.265]                  (face areas, nx+1)
#   V     = [16.755, 16.755]                      (cell volumes)
#   ΔA    = [31.665, 18.600]                      (per-cell, A[i+1]-A[i])
#   μ     = [-0.5774, +0.5774]
#   w     = [1, 1]
#   α_h   = [0, 0.5774, 0]                        (per-half-ordinate, N+1)
#   τ_mm  = [0.5, 0.5774]                         (per-ordinate)
#   c_in[n]  = (1-τ_n)/τ_n · α_{n+1/2} + α_{n-1/2}
#   c_out[n] = α_{n+1/2} / τ_n
#   dAw[i, n] = ΔA[i] / w[n]
#
# Total cross sections (1G mixture B):
#   Σ_t = 2.0, Σ_s = 1.9 (self-scatter), no fission
#   In within-group iteration form: Q_within[i] = Σ_s · φ_0[i] + Q_ext


def per_cell_si_equation(
    psi_avg, n, i, *, sn_mesh, quad, sigma_t, source_per_ord_per_cell,
    psi_spatial_in, psi_angle_in,
):
    """Return the SI per-cell residual at (n, i) for trial psi_avg.

    The SI per-cell equation (in residual form, to make verification
    by substitution explicit) is:

        denom[n, i] · ψ_avg[n, i]
          - source[n, i]·V[i]/Σw  (i.e., already-weight-normalised volumetric source)
          - |μ_n|·(A_in + A_out)·ψ_spatial_in[n, i]
          - dAw[i, n] · c_in[n] · ψ_angle_in[n, i] = 0

    where
        denom[n, i] = 2|μ_n|·A_downstream[n, i] + dAw[i, n]·c_out[n] + Σ_t[i]·V[i]
        A_downstream = A_outer (i+1) for μ ≥ 0, A_inner (i) for μ < 0
        A_in + A_out = A[i] + A[i+1]  (symmetric)
        c_in = (1-τ)/τ·α_out + α_in
        c_out = α_out / τ

    ψ_spatial_in and ψ_angle_in are SUPPLIED by the caller (different
    between SI and apply-matvec — that's the whole point of this audit).
    """
    red = sn_mesh.reduced
    A = red.face_areas
    V = sn_mesh.volumes[:, 0]
    mu = quad.mu_x
    w = quad.weights
    alpha_h = red.alpha_half          # (N+1,) per half-ordinate
    tau = red.tau_mm                  # (N,)
    dAw = red.redist_dAw              # (nx, N)

    abs_mu = abs(mu[n])
    A_total = A[i] + A[i + 1]
    A_downstream = A[i + 1] if mu[n] >= 0 else A[i]
    c_in = (1.0 - tau[n]) / tau[n] * alpha_h[n + 1] + alpha_h[n]
    c_out = alpha_h[n + 1] / tau[n]

    denom = 2.0 * abs_mu * A_downstream + dAw[i, n] * c_out + sigma_t * V[i]
    # source carries Q·V/Σw on the strategy contract
    source_VW = source_per_ord_per_cell  # already-normalised

    residual = (
        denom * psi_avg
        - source_VW
        - abs_mu * A_total * psi_spatial_in
        - dAw[i, n] * c_in * psi_angle_in
    )
    return residual


# ─────────────────────────────────────────────────────────────────────
# Build the per-cell upstream state for each method
# ─────────────────────────────────────────────────────────────────────

def trace_si_upstream(psi, sn_mesh, quad, sigma_t, Q_within_iso, bc_outer_face):
    """Reproduce the SI sweep's per-cell upstream state at convergence.

    Given a candidate fixed-point ψ_SI* (shape (N, nx)) and the
    weight-normalised per-cell-per-ordinate source Q_within (per cell
    only, since iso) shape (nx, ng), trace the sweep and report:
        - psi_spatial_in (the spatial-upstream face flux at this cell)
        - psi_angle_in (the M-M ψ_{n-1/2} face value at this cell)
        - psi_avg (the per-cell-traced cell average)
    """
    N = quad.N
    nx = sn_mesh.nx
    mu = quad.mu_x
    w = quad.weights
    sum_w = w.sum()

    # Carlson seed for psi_angle (uses Q_within_iso=full source, NOT psi field)
    red = sn_mesh.reduced
    bc_outer_obj = sn_mesh.bc_right
    inflow_full = bc_outer_obj.apply(bc_outer_face)  # (N, ng)
    most_inward_idx = int(np.argmin(mu))
    bc_outer_value = inflow_full[most_inward_idx, :]  # (ng,)

    # SI Carlson seed: Q_bar = 0.5 · Q_total (the full within-group source)
    Q_bar = 0.5 * Q_within_iso.T  # (ng, nx)
    sigma_t_gx = np.array([[sigma_t]*nx])
    dr = sn_mesh.dx
    phi_aux = carlson_inward_sweep_from_source(
        Q_bar=Q_bar, sigma_t=sigma_t_gx, dr=dr, bc_outer_value=bc_outer_value,
    )  # (ng, nx) = (1, nx)
    psi_angle_1d = phi_aux.T[:, 0].copy()

    # Build QV_iso per cell
    V = sn_mesh.volumes[:, 0]
    QV_iso = Q_within_iso[:, 0] * V / sum_w  # (nx,) for 1G

    upstream_info = []  # list of (n, i, psi_spatial_in, psi_angle_in, psi_avg_si, psi_spatial_out, psi_angle_out)

    bc_outer = bc_outer_face.copy()

    for n in range(N):
        mu_n = mu[n]
        if mu_n < 0:
            psi_in_full = bc_outer_obj.apply(bc_outer)
            psi_spatial_in = psi_in_full[n, 0]
        else:
            psi_spatial_in = 0.0
        for visit in sn_mesh.iter_cell_visits(ordinate_idx=n):
            i = visit.cell_idx
            # Build per-cell terms from sn_mesh.reduced
            A = red.face_areas
            alpha_h = red.alpha_half
            tau_mm = red.tau_mm
            dAw = red.redist_dAw

            abs_mu = abs(mu_n)
            A_total = A[i] + A[i + 1]
            A_downstream = A[i + 1] if mu_n >= 0 else A[i]
            c_in = (1.0 - tau_mm[n]) / tau_mm[n] * alpha_h[n + 1] + alpha_h[n]
            c_out = alpha_h[n + 1] / tau_mm[n]

            denom = 2.0 * abs_mu * A_downstream + dAw[i, n] * c_out + sigma_t * V[i]
            numer = (QV_iso[i]
                     + abs_mu * A_total * psi_spatial_in
                     + dAw[i, n] * c_in * psi_angle_1d[i])
            psi_avg_si = numer / denom

            psi_spatial_out = 2.0 * psi_avg_si - psi_spatial_in
            psi_angle_out = (psi_avg_si - (1.0 - tau_mm[n]) * psi_angle_1d[i]) / tau_mm[n]

            upstream_info.append({
                "n": n, "i": i, "mu_n": mu_n,
                "psi_spatial_in": psi_spatial_in,
                "psi_angle_in": psi_angle_1d[i],
                "psi_avg_si": psi_avg_si,
                "psi_spatial_out": psi_spatial_out,
                "psi_angle_out": psi_angle_out,
                "denom": denom, "c_in": c_in, "c_out": c_out,
                "A_total": A_total, "A_down": A_downstream,
            })
            # State updates (the SI's defining property)
            psi_angle_1d[i] = psi_angle_out
            psi_spatial_in = psi_spatial_out
        if mu_n >= 0:
            bc_outer[n, 0] = psi_spatial_in
    return upstream_info, psi_angle_1d, bc_outer


def trace_apply_upstream(psi, sn_mesh, quad, sigma_t):
    """Reproduce the apply-matvec's per-cell upstream state.

    Given input ψ field (N, nx), evaluate the SAME per-cell algebra
    structure but with the apply-matvec's input-driven (not in-sweep-
    propagated) upstream state.

    Specifically:
        - psi_angle_in[n, i] = ψ_{n-1/2,i} from M-M recurrence
          starting at the Carlson seed phi_aux[i] from INPUT ψ
        - psi_spatial_in[n, i] = WDD-propagated face flux from
          cell i-1 (outward) or i+1 (inward), using INPUT ψ_cell[n, *]
          for the recurrence ψ_face_out = 2·ψ_cell - ψ_face_in.
    """
    N = quad.N
    nx = sn_mesh.nx
    mu = quad.mu_x
    w = quad.weights

    red = sn_mesh.reduced
    A = red.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha_h = red.alpha_half
    tau_mm = red.tau_mm
    dAw = red.redist_dAw

    bc_outer_obj = sn_mesh.bc_right

    # Phase D Carlson seed: uses INPUT psi field for phi_0
    # phi_0[i] = Σ_n w_n · ψ[n, i]
    phi_0 = (w[:, None] * psi).sum(axis=0)   # (nx,)
    Q_bar = 0.5 * sigma_t * phi_0[None, :]    # (ng=1, nx)
    sigma_t_gx = np.array([[sigma_t]*nx])
    dr = sn_mesh.dx
    # bc_outer_value: BC applied to cell-centre ψ at outermost cell
    # (the apply-matvec's `fi[:, :, -1, 0]` proxy)
    psi_outermost_cell = psi[:, -1].reshape(N, 1)  # (N, ng=1)
    outer_inflow = bc_outer_obj.apply(psi_outermost_cell)  # (N, ng=1)
    most_inward = int(np.argmin(mu))
    bc_outer_value = outer_inflow[most_inward, :]  # (ng=1,)
    phi_aux = carlson_inward_sweep_from_source(
        Q_bar=Q_bar, sigma_t=sigma_t_gx, dr=dr, bc_outer_value=bc_outer_value,
    )  # (ng=1, nx)

    # M-M recurrence (over ordinates m=0..N-1, vectorised over cells)
    # ψ_half_left starts at phi_aux[0, :]; for each m:
    #   ψ_half_right = (ψ[m, :] - (1-τ_m)·ψ_half_left) / τ_m
    # Store ψ_angle_in[n, i] = ψ_half_left at the START of step m=n
    psi_angle_in_apply = np.zeros((N, nx))
    psi_half_left = phi_aux[0, :].copy()
    for m in range(N):
        psi_angle_in_apply[m, :] = psi_half_left.copy()
        psi_half_right = (psi[m, :] - (1.0 - tau_mm[m]) * psi_half_left) / tau_mm[m]
        psi_half_left = psi_half_right

    # WDD spatial: outward pass i=0..nx-1 with ψ_face_in at i=0 set to ψ_cell[0]
    # (the "pole face = cell-centre" choice in operator.py:781-786).
    psi_spatial_in_apply = np.zeros((N, nx))
    outflow_at_boundary = np.zeros(N)

    # outward
    out_mask = mu >= 0
    for n in range(N):
        if not out_mask[n]:
            continue
        psi_face_in = psi[n, 0]  # pole-face init = cell-centre value at i=0
        for i in range(nx):
            psi_spatial_in_apply[n, i] = psi_face_in
            psi_cell = psi[n, i]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            psi_face_in = psi_face_out
        outflow_at_boundary[n] = psi_face_in  # outermost outflow

    # BC apply on outflow
    inflow_full = bc_outer_obj.apply(outflow_at_boundary.reshape(N, 1))  # (N, 1)
    inflow_vec = inflow_full[:, 0]

    # inward
    for n in range(N):
        if out_mask[n]:
            continue
        psi_face_in = inflow_vec[n]
        for i in range(nx - 1, -1, -1):
            psi_spatial_in_apply[n, i] = psi_face_in
            psi_cell = psi[n, i]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            psi_face_in = psi_face_out

    return psi_angle_in_apply, psi_spatial_in_apply, phi_aux


def verify_si_residual(psi_si, sn_mesh, quad, sigma_t, Q_within_iso, bc_outer_face):
    """Show SI fixed point satisfies the SI per-cell equation."""
    info, _, _ = trace_si_upstream(psi_si.copy(), sn_mesh, quad, sigma_t,
                                    Q_within_iso, bc_outer_face)
    print("\nSI per-cell equation check (residual = denom·ψ_avg_si − RHS):")
    print(f"  Visit order:  n, i, psi_avg_si (computed), psi_si[n,i] (state), |diff|")
    max_diff = 0.0
    for d in info:
        diff = abs(d["psi_avg_si"] - psi_si[d["n"], d["i"]])
        max_diff = max(max_diff, diff)
        print(f"  n={d['n']} i={d['i']}  μ={d['mu_n']:+.4f}  "
              f"psi_avg_si={d['psi_avg_si']:.6f}  "
              f"psi_si[n,i]={psi_si[d['n'], d['i']]:.6f}  "
              f"|diff|={diff:.2e}")
    return max_diff, info


def verify_apply_residual(psi, sn_mesh, quad, sigma_t, Q_within_iso):
    """Show that ψ_K satisfies T·ψ = q_total at fixed point (residual ≈ 0)."""
    red = sn_mesh.reduced
    A = red.face_areas
    V = sn_mesh.volumes[:, 0]
    mu = quad.mu_x
    w = quad.weights
    sum_w = w.sum()
    alpha_h = red.alpha_half
    tau_mm = red.tau_mm
    dAw = red.redist_dAw
    N = quad.N
    nx = sn_mesh.nx

    psi_angle_in_apply, psi_spatial_in_apply, _ = trace_apply_upstream(
        psi, sn_mesh, quad, sigma_t,
    )

    # The apply-matvec's LHS for (n, i): streaming + redist + collision = T·ψ
    # The fixed-point equation: T·ψ = (Q_within_iso·V/Σw)/V ?  Actually
    # apply-matvec's lhs is the BALANCE EQUATION DIVIDED BY V (since it
    # uses streaming = mu·(A[i+1]·ψ_out − A[i]·ψ_in) / V[i]).
    # The fixed-point condition: (T·ψ)[n, i] = q[n, i] = Q_within_iso[i] / Σw.

    print("\nApply-matvec per-cell residual check:")
    print(f"  n, i, μ, ψ[n,i], T·ψ[n,i], q[n,i], |residual|")
    max_resid = 0.0
    info = []
    for n in range(N):
        for i in range(nx):
            mu_n = mu[n]
            psi_face_in = psi_spatial_in_apply[n, i]
            psi_cell = psi[n, i]
            psi_face_out = 2.0 * psi_cell - psi_face_in
            # Streaming: mu_n · (A[i+1]·ψ_face_out − A[i]·ψ_face_in) / V[i]  (signed mu_n)
            streaming = mu_n * (A[i + 1] * psi_face_out - A[i] * psi_face_in) / V[i]
            # Redistribution: (dAw / V) · (α_{n+1/2}·ψ_{n+1/2} − α_{n-1/2}·ψ_{n-1/2})
            # where ψ_{n+1/2} is M-M-derived from ψ_{n-1/2} and ψ_cell:
            #   ψ_{n+1/2} = (ψ_cell − (1−τ_n)·ψ_{n-1/2}) / τ_n
            psi_half_left = psi_angle_in_apply[n, i]
            psi_half_right = (psi_cell - (1.0 - tau_mm[n]) * psi_half_left) / tau_mm[n]
            redist = (dAw[i, n] / V[i]) * (
                alpha_h[n + 1] * psi_half_right - alpha_h[n] * psi_half_left
            )
            collision = sigma_t * psi_cell
            T_psi = streaming + redist + collision
            q_ni = Q_within_iso[i, 0] / sum_w  # source / Σw
            residual = T_psi - q_ni
            max_resid = max(max_resid, abs(residual))
            print(f"  n={n} i={i}  μ={mu_n:+.4f}  ψ={psi_cell:.6f}  "
                  f"T·ψ={T_psi:.6f}  q={q_ni:.6f}  |r|={abs(residual):.3e}")
            info.append({
                "n": n, "i": i, "psi_face_in_apply": psi_face_in,
                "psi_face_out_apply": psi_face_out,
                "psi_angle_in_apply": psi_half_left,
                "psi_angle_out_apply": psi_half_right,
                "T_psi": T_psi, "residual": residual,
            })
    return max_resid, info


def cross_verify(psi_si, psi_kr, sn_mesh, quad, sigma_t, Q_within_iso, bc_outer_face):
    """Show ψ_SI* fails apply-matvec residual; ψ_K* fails SI."""
    # Substitute ψ_SI into apply-matvec residual
    print("\n=== Substituting ψ_SI* into apply-matvec residual ===")
    max_resid_si_in_apply, info_si_in_apply = verify_apply_residual(
        psi_si, sn_mesh, quad, sigma_t, Q_within_iso,
    )
    print(f"  max |T·ψ_SI − q| = {max_resid_si_in_apply:.4e}")

    print("\n=== Substituting ψ_K* into apply-matvec residual ===")
    max_resid_kr_in_apply, info_kr_in_apply = verify_apply_residual(
        psi_kr, sn_mesh, quad, sigma_t, Q_within_iso,
    )
    print(f"  max |T·ψ_K − q| = {max_resid_kr_in_apply:.4e}")

    print("\n=== Substituting ψ_SI* into SI per-cell equation ===")
    max_resid_si_in_si, info_si_in_si = verify_si_residual(
        psi_si, sn_mesh, quad, sigma_t, Q_within_iso, bc_outer_face,
    )
    print(f"  max |ψ_avg_si(traced) − ψ_si| = {max_resid_si_in_si:.4e}")

    print("\n=== Substituting ψ_K* into SI sweep trace ===")
    max_resid_kr_in_si, info_kr_in_si = verify_si_residual(
        psi_kr, sn_mesh, quad, sigma_t, Q_within_iso, bc_outer_face,
    )
    print(f"  max |ψ_avg_si(traced from ψ_K) − ψ_K| = {max_resid_kr_in_si:.4e}")

    return info_si_in_apply, info_kr_in_apply, info_si_in_si, info_kr_in_si


def main():
    psi_si, psi_kr, sn_mesh, quad, fuel = get_fixed_points()
    print(f"ψ_SI: {psi_si}")
    print(f"ψ_K:  {psi_kr}")

    # Within-group iteration form: Q_within = Σ_s · φ_0 + Q_ext  (at fixed point)
    sigma_t = float(fuel.SigT[0])
    sigma_s = float(fuel.SigS[0].toarray()[0, 0])
    # external Q = 1 isotropic
    Q_ext_iso = 1.0
    # At fixed point: φ_0 = Σ_n w_n · ψ[n, :]
    w = quad.weights
    phi_0_si = (w[:, None] * psi_si).sum(axis=0)  # (nx,)
    phi_0_kr = (w[:, None] * psi_kr).sum(axis=0)
    print(f"\nφ_0 (SI):    {phi_0_si}")
    print(f"φ_0 (Krylov): {phi_0_kr}")
    Q_within_si = (sigma_s * phi_0_si + Q_ext_iso).reshape(-1, 1)  # (nx, ng=1)
    Q_within_kr = (sigma_s * phi_0_kr + Q_ext_iso).reshape(-1, 1)
    print(f"Q_within_si: {Q_within_si.ravel()}")
    print(f"Q_within_kr: {Q_within_kr.ravel()}")

    # bc_outer_face: at convergence, contains the outgoing face flux at outermost cell
    # for positive-μ ordinates. Trace 100 times to reach SI's actual converged bc_outer.
    def iterate_bc(psi_state, Q_within, n_iter=200):
        bc_outer = np.zeros((quad.N, 1))
        for _ in range(n_iter):
            _, _, bc_outer = trace_si_upstream(
                psi_state.copy(), sn_mesh, quad, sigma_t, Q_within, bc_outer,
            )
        return bc_outer
    bc_outer_face_si = iterate_bc(psi_si, Q_within_si)
    bc_outer_face_kr = iterate_bc(psi_kr, Q_within_kr)
    print(f"\nConverged bc_outer for SI: {bc_outer_face_si.ravel()}")
    print(f"Converged bc_outer for K:  {bc_outer_face_kr.ravel()}")

    # Substitute ψ_SI* into apply-matvec (with Q_within_si)
    print("\n" + "=" * 70)
    print("CROSS-VERIFICATION at the two fixed points")
    print("=" * 70)
    print("\n=== ψ_SI* in apply-matvec residual (with Q_within_si) ===")
    max_resid_si_in_apply, info_si_in_apply = verify_apply_residual(
        psi_si, sn_mesh, quad, sigma_t, Q_within_si,
    )
    print(f"  max |T·ψ_SI − q| = {max_resid_si_in_apply:.4e}")

    print("\n=== ψ_K* in apply-matvec residual (with Q_within_kr) ===")
    max_resid_kr_in_apply, info_kr_in_apply = verify_apply_residual(
        psi_kr, sn_mesh, quad, sigma_t, Q_within_kr,
    )
    print(f"  max |T·ψ_K − q| = {max_resid_kr_in_apply:.4e}")

    print("\n=== ψ_SI* in SI sweep trace (with Q_within_si) ===")
    max_resid_si_in_si, info_si_in_si = verify_si_residual(
        psi_si, sn_mesh, quad, sigma_t, Q_within_si, bc_outer_face_si,
    )
    print(f"  max |ψ_avg_si(traced) − ψ_si| = {max_resid_si_in_si:.4e}")

    print("\n=== ψ_K* in SI sweep trace (with Q_within_kr) ===")
    max_resid_kr_in_si, info_kr_in_si = verify_si_residual(
        psi_kr, sn_mesh, quad, sigma_t, Q_within_kr, bc_outer_face_kr,
    )
    print(f"  max |ψ_avg_si(traced from ψ_K) − ψ_K| = {max_resid_kr_in_si:.4e}")

    # Now compute the TERM-BY-TERM difference of the two operators at (n=1, i=0)
    # (the pole on the outward sweep)
    print("\n" + "=" * 70)
    print("Per-cell term comparison at (n=1, i=0) [pole, outward ordinate]")
    print("=" * 70)
    for d_si, d_apply in zip(info_si_in_si, info_si_in_apply):
        if d_si["n"] == 1 and d_si["i"] == 0:
            print(f"\nSI sweep trace at (n=1, i=0) using ψ_SI*:")
            print(f"  psi_spatial_in (SI)   = {d_si['psi_spatial_in']:.6f}")
            print(f"  psi_angle_in   (SI)   = {d_si['psi_angle_in']:.6f}")
            print(f"  psi_avg_si (traced)  = {d_si['psi_avg_si']:.6f}")
            print(f"\nApply-matvec trace at (n=1, i=0) using ψ_SI*:")
            print(f"  psi_face_in_apply  = {d_apply['psi_face_in_apply']:.6f}")
            print(f"  psi_angle_in_apply = {d_apply['psi_angle_in_apply']:.6f}")
            print(f"  T·ψ_SI − q (residual) = {d_apply['residual']:.4e}")
            print(f"\nINPUT-STATE difference at (n=1, i=0):")
            print(f"  Δ(psi_spatial_in)  = {d_apply['psi_face_in_apply'] - d_si['psi_spatial_in']:.6e}")
            print(f"  Δ(psi_angle_in)    = {d_apply['psi_angle_in_apply'] - d_si['psi_angle_in']:.6e}")


if __name__ == "__main__":
    main()
