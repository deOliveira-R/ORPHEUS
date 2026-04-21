"""Fractional µ-weight scan: where do the missing factors of µ go?

Created by numerics-investigator on 2026-04-21 for Issue #119.

At the Sanchez µ-ortho recipe we get 1.43% at N=2 (r_0/R=0.3, σ_t R=5).
The two closest prior-state recipes (µ-ortho all-µ, f=const Pµ=Gµ=1) each
gave ~1.4-10.9% residual. Something in the µ-weight placement is off.

Hypothesis: We're double-counting or missing a µ factor somewhere. Let me
scan with fractional µ exponents on {P, G, W} and find the minimum.

Integrand:
  P^n(r_i) ∝ ∫ W_Ω · µ^α_P · f^n(µ) · K_esc(τ) dΩ
  G^n(r_i) ∝ ∫ sin θ · µ^α_G · f^n(µ) · e^{-τ} dθ
  W^mn ∝ ∫ cos θ · sin θ · (µ_e · µ_a)^α_W · f^m(µ_a) · f^n(µ_e) · e^{-τ} dθ

Best-so-far at N=2: Sanchez basis with α_P=α_G=1, α_W=0, no Jacobian → 1.6%.

We also saw: µ-ortho basis, α_P=α_G=α_W=1, no Jacobian → 1.43%.

Let me expand the scan to cover combinations of (α_P, α_G, α_W) ∈ {0,1}×{0,1}×{0,1}.
Also swap the sign of the mode-1 P/G terms (could be a sign convention).
"""
import numpy as np
import pytest
import mpmath

from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry, composite_gl_r, build_volume_kernel,
    gl_float, _shifted_legendre_eval,
    BoundaryClosureOperator,
)


def f_mu_orthonormal(n):
    if n == 0:
        return lambda mu: np.sqrt(2.0)
    elif n == 1:
        return lambda mu: 6.0 * mu - 4.0
    elif n == 2:
        return lambda mu: np.sqrt(6.0) * (10.0 * mu**2 - 12.0 * mu + 3.0)
    raise ValueError


def compute_P_esc_face_mode(
    geometry, r_nodes, radii, sig_t, n_mode, face, *,
    f_basis, alpha_mu,
    n_angular=32, dps=25,
):
    """General P_esc with µ^alpha_mu weight in the integrand."""
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
    f_n = f_basis(n_mode)

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

            weight = mu_s ** alpha_mu if alpha_mu > 0 else 1.0
            total += (
                omega_wts[k] * angular_factor[k]
                * weight * f_n(mu_s) * K_esc
            )
        P[i] = pref * total
    return P


def compute_G_bc_face_mode(
    geometry, r_nodes, radii, sig_t, n_mode, face, *,
    f_basis, alpha_mu,
    n_surf_quad=32, dps=25,
):
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
    f_n = f_basis(n_mode)

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
            weight = mu_s ** alpha_mu if alpha_mu > 0 else 1.0
            total += (
                theta_wts[k] * st * weight * f_n(mu_s) * float(np.exp(-tau))
            )
        G[i] = 2.0 * total
    return G


def compute_W_rank_n(
    geometry, radii, sig_t, N, *, f_basis, alpha_W_emit, alpha_W_arrive, dps=25,
):
    """W with fractional µ weights on emit and arrive sides."""
    R = float(radii[-1])
    r_0 = float(geometry.inner_radius)
    theta_c = float(np.arcsin(r_0 / R))
    sig_t_val = float(sig_t[0])
    dim = 2 * N
    W = np.zeros((dim, dim))

    with mpmath.workdps(dps):
        for m in range(N):
            for n in range(N):
                f_m = f_basis(m)
                f_n = f_basis(n)
                def integrand(th, f_m=f_m, f_n=f_n):
                    cos_th = float(mpmath.cos(th))
                    sin_th = float(mpmath.sin(th))
                    w_e = cos_th ** alpha_W_emit if alpha_W_emit > 0 else 1.0
                    w_a = cos_th ** alpha_W_arrive if alpha_W_arrive > 0 else 1.0
                    return (
                        cos_th * sin_th * w_e * w_a
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
                f_m = f_basis(m)
                f_n = f_basis(n)
                def integrand(th, f_m=f_m, f_n=f_n):
                    cos_th = float(mpmath.cos(th))
                    sin_th = float(mpmath.sin(th))
                    c_in_sq = 1.0 - (R * sin_th / r_0) ** 2
                    if c_in_sq < 0.0:
                        c_in_sq = 0.0
                    c_in = float(mpmath.sqrt(mpmath.mpf(c_in_sq)))
                    w_e = cos_th ** alpha_W_emit if alpha_W_emit > 0 else 1.0
                    w_a = c_in ** alpha_W_arrive if alpha_W_arrive > 0 else 1.0
                    return (
                        cos_th * sin_th * w_e * w_a
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
               f_basis, alpha_P, alpha_G, alpha_W_e, alpha_W_a,
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
            P_arr = compute_P_esc_face_mode(
                geom, r_nodes, radii, sig_t, n, face,
                f_basis=f_basis, alpha_mu=alpha_P,
                n_angular=n_angular, dps=dps,
            )
            G_arr = compute_G_bc_face_mode(
                geom, r_nodes, radii, sig_t, n, face,
                f_basis=f_basis, alpha_mu=alpha_G,
                n_surf_quad=n_surf_quad, dps=dps,
            )
            row = face_idx * N + n
            P[row, :] = rv * r_wts * P_arr
            if face_idx == 0:
                G[:, row] = sig_t_n * G_arr / div_outer
            else:
                G[:, row] = sig_t_n * G_arr / div_inner

    W = compute_W_rank_n(
        geom, radii, sig_t, N,
        f_basis=f_basis,
        alpha_W_emit=alpha_W_e, alpha_W_arrive=alpha_W_a,
        dps=dps,
    )
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


def test_full_scan():
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

    print(f"\nk_inf = {k_inf}")
    print(f"{'(αP, αG, αWe, αWa)':<25s} {'k_eff':>10s} {'err%':>10s}")
    print("-" * 55)

    # Scan α ∈ {0, 1} — binary exponents
    best_err = float('inf')
    best_combo = None
    for a_P in [0, 1]:
        for a_G in [0, 1]:
            for a_We in [0, 1]:
                for a_Wa in [0, 1]:
                    try:
                        K_bc = build_K_bc(
                            geom, r_nodes, r_wts, radii, sig_t, 2,
                            f_basis=f_mu_orthonormal,
                            alpha_P=a_P, alpha_G=a_G,
                            alpha_W_e=a_We, alpha_W_a=a_Wa,
                        )
                        k_eff = solve_k_eff(K_bc, K_vol, sig_t_v, sig_s_v, nu_sig_f_v)
                        err = abs(k_eff - k_inf) / k_inf * 100
                        combo = (a_P, a_G, a_We, a_Wa)
                        mark = ""
                        if err < best_err:
                            best_err = err
                            best_combo = combo
                            mark = " *BEST*"
                        print(f"{str(combo):<25s} {k_eff:10.6f} {err:9.3f}%{mark}")
                    except Exception as e:
                        print(f"{str((a_P, a_G, a_We, a_Wa)):<25s} ERROR: {e}")
    print(f"\nBest recipe: {best_combo} with err {best_err:.4f}%")


if __name__ == "__main__":
    test_full_scan()
