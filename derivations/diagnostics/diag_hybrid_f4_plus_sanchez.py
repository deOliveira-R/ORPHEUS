"""Hybrid: F.4 mode 0 + Sanchez µ-weighted mode n≥1.

Created by numerics-investigator on 2026-04-21 for Issue #119.

At N=1, pure Sanchez µ-ortho gives 2.55% error — that's 30× worse than
F.4's 0.077%. The Sanchez recipe's mode-0 is different from F.4's.

Strategy: use F.4's compute_P_esc_outer / compute_G_bc_outer for mode 0
(so N=1 closure is bit-exact with F.4 at 0.077%) and stack µ-weighted
Sanchez primitives for n ≥ 1.

Consistency concern: this mixes bases (F.4's mode 0 is Lambert angular-
flux scalar, Sanchez's mode n≥1 is µ-weighted partial-current moment).
But if the W matrix treats the blocks accordingly, it may still close.

For W:
  W_{00}: F.4 scalar transmission (no µ)
  W_{mn} for m,n ≥ 1: Sanchez µ-weighted moment (µ in integrand)
  W_{0n}, W_{n0} for n ≥ 1: mixed — emit is F.4 (no µ), arrive is
    Sanchez (µ on arrival), OR vice versa.
"""
import numpy as np
import mpmath

from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry, composite_gl_r, build_volume_kernel,
    gl_float, _shifted_legendre_eval,
    compute_P_esc_outer, compute_P_esc_inner,
    compute_G_bc_outer, compute_G_bc_inner,
    compute_hollow_sph_transmission,
)


def f_mu_orthonormal(n):
    if n == 0:
        return lambda mu: np.sqrt(2.0)
    elif n == 1:
        return lambda mu: 6.0 * mu - 4.0
    elif n == 2:
        return lambda mu: np.sqrt(6.0) * (10.0 * mu**2 - 12.0 * mu + 3.0)
    raise ValueError


def compute_P_esc_mode_n1_plus(
    geometry, r_nodes, radii, sig_t, n_mode, face, *,
    n_angular=32, dps=25,
):
    """Sanchez µ-weighted primitive for n ≥ 1. Drop-in replacement for
    modes beyond 0."""
    if n_mode == 0:
        raise ValueError("Use F.4 compute_P_esc_outer/inner for mode 0")

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
            total += omega_wts[k] * angular_factor[k] * mu_s * f_n(mu_s) * K_esc
        P[i] = pref * total
    return P


def compute_G_bc_mode_n1_plus(
    geometry, r_nodes, radii, sig_t, n_mode, face, *,
    n_surf_quad=32, dps=25,
):
    if n_mode == 0:
        raise ValueError("Use F.4 compute_G_bc_outer/inner for mode 0")

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
            total += theta_wts[k] * st * mu_s * f_n(mu_s) * float(np.exp(-tau))
        G[i] = 2.0 * total
    return G


def compute_hybrid_W(geom, radii, sig_t, N, dps=25):
    """W with mode 0 = F.4 scalar, modes 1..N-1 = µ-weighted.

    W has blocks for [outer_0, outer_1, ..., outer_{N-1}, inner_0, ..., inner_{N-1}].
    For mode 0, use the F.4 scalar compute_hollow_sph_transmission block.
    For mode n ≥ 1, use µ-weighted integrand. Cross-blocks are µ on arrival only
    or µ on emit only, depending on which mode is non-zero.
    """
    R = float(radii[-1])
    r_0 = float(geom.inner_radius)
    theta_c = float(np.arcsin(r_0 / R))
    sig_t_val = float(sig_t[0])
    dim = 2 * N
    W = np.zeros((dim, dim))

    # F.4 scalar W (2×2)
    W_scalar = compute_hollow_sph_transmission(r_0, R, radii, sig_t, dps=dps)
    # Place F.4 W into the (0,0) "mode-0 both faces" block
    # Mode-0 outer is index 0, mode-0 inner is index N.
    W[0, 0] = W_scalar[0, 0]
    W[0, N] = W_scalar[0, 1]
    W[N, 0] = W_scalar[1, 0]
    W[N, N] = W_scalar[1, 1]

    # For pure mode n ≥ 1 entries and mode 0 × mode n ≥ 1 cross entries,
    # use the µ-weighted integrand consistent with Sanchez.
    with mpmath.workdps(dps):
        # W_oo block (emission at outer → arrival at outer)
        for m in range(N):
            for n in range(N):
                if m == 0 and n == 0:
                    continue  # already set from F.4
                f_m = f_mu_orthonormal(m)
                f_n = f_mu_orthonormal(n)
                # Use µ-weight consistent with higher modes.
                # mode 0 basis: f^0 = √2 (NOT F.4's constant 1 here)
                # This is a basis mismatch — but the F.4 W_00 is already loaded.
                # Cross entries must use the SAME basis as their mode ≥ 1
                # function. So f_n for n=0 here must be √2 (not 1) for the
                # cross-terms — else the inner product isn't well-defined.
                def integrand(th, f_m=f_m, f_n=f_n):
                    cos_th = float(mpmath.cos(th))
                    sin_th = float(mpmath.sin(th))
                    # µ on both sides per Sanchez
                    return (
                        cos_th * sin_th * cos_th * cos_th
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
                if m == 0 and n == 0:
                    continue  # already set
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
                        cos_th * sin_th * cos_th * c_in
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
                if m == 0 and n == 0:
                    continue
                W[m, N + n] = area_ratio * W[N + n, m]
    return W


def build_K_bc_hybrid(geom, r_nodes, r_wts, radii, sig_t, N, *,
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

    # Mode 0: F.4 scalar primitives
    P_out_0 = compute_P_esc_outer(geom, r_nodes, radii, sig_t, n_angular=n_angular, dps=dps)
    P_in_0 = compute_P_esc_inner(geom, r_nodes, radii, sig_t, n_angular=n_angular, dps=dps)
    G_out_0 = compute_G_bc_outer(geom, r_nodes, radii, sig_t, n_surf_quad=n_surf_quad, dps=dps)
    G_in_0 = compute_G_bc_inner(geom, r_nodes, radii, sig_t, n_surf_quad=n_surf_quad, dps=dps)
    P[0, :] = rv * r_wts * P_out_0
    P[N, :] = rv * r_wts * P_in_0
    G[:, 0] = sig_t_n * G_out_0 / div_outer
    G[:, N] = sig_t_n * G_in_0 / div_inner

    # Modes n ≥ 1: Sanchez µ-weighted
    for n in range(1, N):
        P_out_n = compute_P_esc_mode_n1_plus(
            geom, r_nodes, radii, sig_t, n, "outer",
            n_angular=n_angular, dps=dps,
        )
        P_in_n = compute_P_esc_mode_n1_plus(
            geom, r_nodes, radii, sig_t, n, "inner",
            n_angular=n_angular, dps=dps,
        )
        G_out_n = compute_G_bc_mode_n1_plus(
            geom, r_nodes, radii, sig_t, n, "outer",
            n_surf_quad=n_surf_quad, dps=dps,
        )
        G_in_n = compute_G_bc_mode_n1_plus(
            geom, r_nodes, radii, sig_t, n, "inner",
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[n, :] = rv * r_wts * P_out_n
        P[N + n, :] = rv * r_wts * P_in_n
        G[:, n] = sig_t_n * G_out_n / div_outer
        G[:, N + n] = sig_t_n * G_in_n / div_inner

    W = compute_hybrid_W(geom, radii, sig_t, N, dps=dps)
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


def test_hybrid_n_convergence():
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

    print(f"\nHybrid F.4 + Sanchez: k_inf={k_inf}, R={R}, r_0={r_0}")
    print(f"{'N':>3s} {'k_eff':>12s} {'err%':>10s}")
    print("-" * 30)
    for N in [1, 2, 3]:
        try:
            K_bc = build_K_bc_hybrid(geom, r_nodes, r_wts, radii, sig_t, N)
            k_eff = solve_k_eff(K_bc, K_vol, sig_t_v, sig_s_v, nu_sig_f_v)
            err = abs(k_eff - k_inf) / k_inf * 100
            print(f"{N:3d} {k_eff:12.6f} {err:9.4f}%")
        except Exception as e:
            print(f"{N:3d} ERROR: {e}")


if __name__ == "__main__":
    test_hybrid_n_convergence()
