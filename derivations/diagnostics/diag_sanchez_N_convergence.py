"""Convergence of the Sanchez µ-ortho recipe with N.

Created by numerics-investigator on 2026-04-21 for Issue #119.

Best recipe so far (from diag_sanchez_fractional_scan.py):
    f_basis = µ-orthonormal (f^0=√2, f^1=6µ-4, ...)
    α_P = α_G = α_We = α_Wa = 1 (µ-weight on all integrands)
    No Jacobian
    Model A (both P and G skip blocked rays)
    No Gelbard

At N=2: 1.43% error
At N=1: what does it give? (should be some finite value ~ F.4)
At N=3, 4: does it decrease toward 0?

If yes, the recipe IS converging but slowly. We may need N=4-5 to hit 0.1%.
If no (plateau or increase), there's still a bug.
"""
import numpy as np
import pytest
import mpmath

from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry, composite_gl_r, build_volume_kernel,
    gl_float,
)


def f_mu_orthonormal(n):
    if n == 0:
        return lambda mu: np.sqrt(2.0)
    elif n == 1:
        return lambda mu: 6.0 * mu - 4.0
    elif n == 2:
        return lambda mu: np.sqrt(6.0) * (10.0 * mu**2 - 12.0 * mu + 3.0)
    elif n == 3:
        return lambda mu: np.sqrt(2.0) * (70.0 * mu**3 - 120.0 * mu**2 + 60.0 * mu - 8.0)
    raise ValueError(f"n_mode {n} not implemented (need more polynomials)")


def compute_P_esc(geometry, r_nodes, radii, sig_t, n_mode, face, *,
                  alpha_mu=1.0, n_angular=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    r_0 = float(geometry.inner_radius)
    omega_low, omega_high = geometry.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geometry.ray_direction_cosine(omega_pts)
    angular_factor = geometry.angular_weight(omega_pts)
    pref = geometry.prefactor
    N_r = len(r_nodes)
    P = np.zeros(N_r)
    f_n = f_mu_orthonormal(n_mode)

    for i in range(N_r):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            if face == "outer":
                rho_out = geometry.rho_max(r_i, cos_om, R)
                if rho_out <= 0.0:
                    continue
                rho_in_minus, _ = geometry.rho_inner_intersections(r_i, cos_om)
                if rho_in_minus is not None and rho_in_minus < rho_out:
                    continue
                tau = geometry.optical_depth_along_ray(
                    r_i, cos_om, rho_out, radii, sig_t,
                )
                K_esc = geometry.escape_kernel_mp(tau, dps)
                mu_s = (rho_out + r_i * cos_om) / R
            elif face == "inner":
                rho_in_minus, _ = geometry.rho_inner_intersections(r_i, cos_om)
                if rho_in_minus is None:
                    continue
                tau = geometry.optical_depth_along_ray(
                    r_i, cos_om, rho_in_minus, radii, sig_t,
                )
                K_esc = geometry.escape_kernel_mp(tau, dps)
                sin_om = np.sqrt(max(0.0, 1.0 - cos_om * cos_om))
                h_sq = r_i * r_i * sin_om * sin_om
                mu_s_sq = max(0.0, (r_0 * r_0 - h_sq) / (r_0 * r_0))
                mu_s = float(np.sqrt(mu_s_sq))
            total += (
                omega_wts[k] * angular_factor[k]
                * mu_s * f_n(mu_s) * K_esc
            )
        P[i] = pref * total
    return P


def compute_G_bc(geometry, r_nodes, radii, sig_t, n_mode, face, *,
                 n_surf_quad=32, dps=25):
    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    r_0 = float(geometry.inner_radius)
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N_r = len(r_nodes)
    G = np.zeros(N_r)
    f_n = f_mu_orthonormal(n_mode)

    for i in range(N_r):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            if face == "outer":
                rho_out = geometry.rho_max(r_i, ct, R)
                if rho_out <= 0.0:
                    continue
                rho_in_minus, _ = geometry.rho_inner_intersections(r_i, ct)
                if rho_in_minus is not None and rho_in_minus < rho_out:
                    continue
                if len(radii) == 1:
                    tau = sig_t[0] * rho_out
                else:
                    tau = geometry.optical_depth_along_ray(
                        r_i, ct, rho_out, radii, sig_t,
                    )
                mu_s = (rho_out + r_i * ct) / R
            elif face == "inner":
                rho_in_minus, _ = geometry.rho_inner_intersections(r_i, ct)
                if rho_in_minus is None:
                    continue
                if len(radii) == 1:
                    tau = sig_t[0] * rho_in_minus
                else:
                    tau = geometry.optical_depth_along_ray(
                        r_i, ct, rho_in_minus, radii, sig_t,
                    )
                sin_om = float(np.sqrt(max(0.0, 1.0 - ct * ct)))
                h_sq = r_i * r_i * sin_om * sin_om
                mu_s_sq = max(0.0, (r_0 * r_0 - h_sq) / (r_0 * r_0))
                mu_s = float(np.sqrt(mu_s_sq))
            total += (
                theta_wts[k] * st * mu_s * f_n(mu_s) * float(np.exp(-tau))
            )
        G[i] = 2.0 * total
    return G


def compute_W_rank_n(geometry, radii, sig_t, N, dps=25):
    R = float(radii[-1])
    r_0 = float(geometry.inner_radius)
    theta_c = float(np.arcsin(r_0 / R))
    sig_t_val = float(sig_t[0])
    dim = 2 * N
    W = np.zeros((dim, dim))

    with mpmath.workdps(dps):
        for m in range(N):
            for n in range(N):
                f_m = f_mu_orthonormal(m)
                f_n = f_mu_orthonormal(n)
                def integrand(th, f_m=f_m, f_n=f_n):
                    cos_th = float(mpmath.cos(th))
                    sin_th = float(mpmath.sin(th))
                    return (
                        cos_th * sin_th
                        * cos_th * cos_th  # µ_e · µ_a weights
                        * f_n(cos_th) * f_m(cos_th)
                        * mpmath.exp(-2.0 * sig_t_val * R * cos_th)
                    )
                val = 2.0 * float(mpmath.quad(
                    integrand, [mpmath.mpf(theta_c), mpmath.pi / 2],
                ))
                W[m, n] = val

        def _chord(t):
            h_sq = R * R * float(mpmath.sin(t)) ** 2
            return R * float(mpmath.cos(t)) - float(
                mpmath.sqrt(mpmath.mpf(r_0 * r_0 - h_sq))
            )
        for m in range(N):
            for n in range(N):
                f_m = f_mu_orthonormal(m)
                f_n = f_mu_orthonormal(n)
                def integrand(th, f_m=f_m, f_n=f_n):
                    cos_th = float(mpmath.cos(th))
                    sin_th = float(mpmath.sin(th))
                    c_in_sq = 1.0 - (R * sin_th / r_0) ** 2
                    if c_in_sq < 0.0:
                        c_in_sq = 0.0
                    c_in = float(mpmath.sqrt(mpmath.mpf(c_in_sq)))
                    return (
                        cos_th * sin_th
                        * cos_th * c_in  # µ_e emit, µ_a arrive
                        * f_n(cos_th) * f_m(c_in)
                        * mpmath.exp(-sig_t_val * _chord(th))
                    )
                val = 2.0 * float(mpmath.quad(
                    integrand, [mpmath.mpf(0.0), mpmath.mpf(theta_c)],
                ))
                W[N + m, n] = val

        area_ratio = (R / r_0) ** 2
        for m in range(N):
            for n in range(N):
                W[m, N + n] = area_ratio * W[N + n, m]
    return W


def build_K_bc(geom, r_nodes, r_wts, radii, sig_t, N, *,
               n_angular=24, n_surf_quad=24, dps=15):
    N_r = len(r_nodes)
    P = np.zeros((2 * N, N_r))
    G = np.zeros((N_r, 2 * N))

    R_out = float(radii[-1])
    r_in = float(geom.inner_radius)
    div_outer = R_out * R_out
    div_inner = r_in * r_in

    sig_t_n = np.array([sig_t[geom.which_annulus(ri, radii)] for ri in r_nodes])
    rv = np.array([geom.radial_volume_weight(rj) for rj in r_nodes])

    for n in range(N):
        for face_idx, face in enumerate(["outer", "inner"]):
            P_arr = compute_P_esc(
                geom, r_nodes, radii, sig_t, n, face,
                n_angular=n_angular, dps=dps,
            )
            G_arr = compute_G_bc(
                geom, r_nodes, radii, sig_t, n, face,
                n_surf_quad=n_surf_quad, dps=dps,
            )
            row = face_idx * N + n
            P[row, :] = rv * r_wts * P_arr
            if face_idx == 0:
                G[:, row] = sig_t_n * G_arr / div_outer
            else:
                G[:, row] = sig_t_n * G_arr / div_inner

    W = compute_W_rank_n(geom, radii, sig_t, N, dps=dps)
    R_mat = np.linalg.inv(np.eye(2 * N) - W)
    return G @ R_mat @ P


def solve_k_eff(K_bc, K_vol, sig_t_v, sig_s_v, nu_sig_f_v):
    import scipy.linalg as sla
    K = K_vol + K_bc
    N_r = K.shape[0]
    A = np.diag(sig_t_v * np.ones(N_r)) - K * sig_s_v
    B = K * nu_sig_f_v
    eigvals = sla.eigvals(np.linalg.solve(A, B))
    return float(np.max(np.real(eigvals)))


def test_sanchez_mu_ortho_plateaus_above_0p1_percent():
    """REGRESSION GATE: the Sanchez µ-weighted orthonormal recipe for
    hollow-sphere rank-N per-face white-BC closure plateaus at ~1.43 %
    k_eff residual for N ≥ 2 and does NOT reach the 0.1 % target.

    This test PASSES as long as the recipe is still INSUFFICIENT. If a
    future fix closes the gap (N=2 error < 0.1 %), this test will FAIL
    and should be retired in favour of a direct success gate.

    Setup: R=5, r_0/R=0.3, Σ_t=1, Σ_s=0.5, νΣ_f=0.75 (k_inf=1.5).
    """
    sig_t_v, sig_s_v, nu_sig_f_v = 1.0, 0.5, 0.75
    k_inf = nu_sig_f_v / (sig_t_v - sig_s_v)
    R = 5.0
    r_0 = 1.5
    geom = CurvilinearGeometry(kind='sphere-1d', inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([sig_t_v])
    r_nodes, r_wts, panels = composite_gl_r(radii, 2, 4, dps=15, inner_radius=r_0)
    K_vol = build_volume_kernel(
        geom, r_nodes, panels, radii, sig_t,
        n_angular=24, n_rho=24, dps=15,
    )

    errs = {}
    for N in [1, 2, 3]:
        K_bc = build_K_bc(geom, r_nodes, r_wts, radii, sig_t, N)
        k_eff = solve_k_eff(K_bc, K_vol, sig_t_v, sig_s_v, nu_sig_f_v)
        errs[N] = abs(k_eff - k_inf) / k_inf * 100

    # N=1 ~ 2.55 %, N=2 ~ 1.42 %, N=3 ~ 1.42 % (plateau)
    assert errs[1] > 0.1, (
        f"Sanchez µ-ortho N=1 error {errs[1]:.4f}% now below 0.1% — "
        f"this regression gate must be revisited."
    )
    assert errs[2] > 0.1, (
        f"Sanchez µ-ortho N=2 error {errs[2]:.4f}% now below 0.1% — "
        f"Issue #119 fix landed! Retire this gate."
    )
    # Confirm plateau: N=3 within 10 % of N=2
    assert abs(errs[3] - errs[2]) / errs[2] < 0.10, (
        f"N=2 err {errs[2]:.4f}% vs N=3 err {errs[3]:.4f}% — "
        f"recipe NOT plateauing, re-examine."
    )
    # Confirm recipe at 1.4 % order of magnitude
    assert 1.0 < errs[2] < 2.0, (
        f"N=2 err {errs[2]:.4f}% outside expected 1-2% plateau band. "
        f"If much lower, investigate recipe improvement. If higher, "
        f"investigate regression."
    )


def main():
    """Print convergence table for manual inspection."""
    sig_t_v, sig_s_v, nu_sig_f_v = 1.0, 0.5, 0.75
    k_inf = nu_sig_f_v / (sig_t_v - sig_s_v)
    R = 5.0
    r_0 = 1.5
    geom = CurvilinearGeometry(kind='sphere-1d', inner_radius=r_0)
    radii = np.array([R])
    sig_t = np.array([sig_t_v])
    r_nodes, r_wts, panels = composite_gl_r(radii, 2, 4, dps=15, inner_radius=r_0)
    K_vol = build_volume_kernel(
        geom, r_nodes, panels, radii, sig_t,
        n_angular=24, n_rho=24, dps=15,
    )

    print(f"\nk_inf = {k_inf}, R = {R}, r_0 = {r_0}, σ_t R = {sig_t_v*R}")
    print(f"{'N':>3s} {'k_eff':>12s} {'err%':>10s}")
    print("-" * 30)

    for N in [1, 2, 3, 4]:
        try:
            K_bc = build_K_bc(geom, r_nodes, r_wts, radii, sig_t, N)
            k_eff = solve_k_eff(K_bc, K_vol, sig_t_v, sig_s_v, nu_sig_f_v)
            err = abs(k_eff - k_inf) / k_inf * 100
            print(f"{N:3d} {k_eff:12.6f} {err:9.4f}%")
        except Exception as e:
            print(f"{N:3d} ERROR: {e}")


if __name__ == "__main__":
    main()
