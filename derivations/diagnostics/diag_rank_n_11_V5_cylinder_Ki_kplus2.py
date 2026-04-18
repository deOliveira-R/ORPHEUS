"""V5: cylinder with Ki_{k+2} polynomial expansion (Phase C) + R = B^{-1}."""
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")
import numpy as np
from numpy.polynomial import Polynomial, Legendre
from orpheus.derivations import peierls_geometry as pg
from orpheus.derivations._kernels import _shifted_legendre_eval, ki_n_mp
from orpheus.derivations.peierls_geometry import (
    SPHERE_1D, CYLINDER_1D, gl_float,
    compute_P_esc_mode, compute_G_bc_mode,
    BoundaryClosureOperator,
)


def shifted_legendre_coeffs(n_max):
    """p_{n,k} = coeff of x^k in P̃_n(x) = P_n(2x-1)."""
    coeffs = np.zeros((n_max + 1, n_max + 1))
    for n in range(n_max + 1):
        leg = np.zeros(n + 1); leg[n] = 1.0
        P_poly = Legendre(leg).convert(kind=Polynomial)
        shifted = Polynomial([0.0])
        for k, c in enumerate(P_poly.coef):
            shifted = shifted + c * Polynomial([-1.0, 2.0]) ** k
        coefs = shifted.coef
        coeffs[n, :len(coefs)] = coefs
    return coeffs


LEG_COEFFS = shifted_legendre_coeffs(10)


def compute_B(N):
    mu, w = gl_float(40, 0.0, 1.0, dps=30)
    B = np.zeros((N, N))
    for n in range(N):
        pn = _shifted_legendre_eval(n, mu)
        for m in range(N):
            pm = _shifted_legendre_eval(m, mu)
            B[n, m] = float(np.sum(w * mu * pn * pm))
    return B


def compute_P_cyl_v5(geometry, r_nodes, radii, sig_t, n_mode, n_beta=32, dps=25):
    """Cylinder P_esc_mode via Ki_{k+2} polynomial expansion + (ρ_2D/R)² Jacobian."""
    R_cell = float(radii[-1])
    inv_R2 = 1.0 / (R_cell * R_cell)
    beta_pts, beta_wts = gl_float(n_beta, 0.0, np.pi, dps)
    cos_betas = np.cos(beta_pts)
    N_r = len(r_nodes)
    P = np.zeros(N_r)
    p_n = LEG_COEFFS[n_mode, :n_mode + 1]
    for i in range(N_r):
        r_i = r_nodes[i]
        total = 0.0
        for k_b in range(n_beta):
            cb = cos_betas[k_b]
            rho_2d = geometry.rho_max(r_i, cb, R_cell)
            if rho_2d <= 0.0:
                continue
            if len(radii) == 1:
                tau = sig_t[0] * rho_2d
            else:
                tau = geometry.optical_depth_along_ray(
                    r_i, cb, rho_2d, np.asarray(radii), np.asarray(sig_t))
            mu_2d = (rho_2d + r_i * cb) / R_cell
            jacobian = rho_2d * rho_2d * inv_R2
            sum_ki = 0.0
            for k in range(n_mode + 1):
                ki_kp2 = float(ki_n_mp(k + 2, float(tau), dps))
                sum_ki += p_n[k] * (mu_2d ** k) * ki_kp2
            total += beta_wts[k_b] * jacobian * sum_ki
        P[i] = (1.0 / np.pi) * total
    return P


def compute_G_cyl_v5(geometry, r_nodes, radii, sig_t, n_mode, n_beta=32, dps=25):
    """Cylinder G_bc_mode via Ki_{k+2} observer-centered, NO Jacobian."""
    R_cell = float(radii[-1])
    beta_pts, beta_wts = gl_float(n_beta, 0.0, np.pi, dps)
    cos_betas = np.cos(beta_pts)
    N_r = len(r_nodes)
    G = np.zeros(N_r)
    p_n = LEG_COEFFS[n_mode, :n_mode + 1]
    for i in range(N_r):
        r_i = r_nodes[i]
        total = 0.0
        for k_b in range(n_beta):
            cb = cos_betas[k_b]
            rho_2d = geometry.rho_max(r_i, cb, R_cell)
            if rho_2d <= 0.0:
                continue
            if len(radii) == 1:
                tau = sig_t[0] * rho_2d
            else:
                tau = geometry.optical_depth_along_ray(
                    r_i, cb, rho_2d, np.asarray(radii), np.asarray(sig_t))
            mu_2d = (rho_2d + r_i * cb) / R_cell
            sum_ki = 0.0
            for k in range(n_mode + 1):
                ki_kp2 = float(ki_n_mp(k + 2, float(tau), dps))
                sum_ki += p_n[k] * (mu_2d ** k) * ki_kp2
            total += beta_wts[k_b] * sum_ki
        # Prefactor: cylinder observer-centered form gives 4·Ki_2 at r=0 for mode 0.
        # (1/π) · ∫_0^π Ki_2(τ) dβ · 4 = 4·Ki_2(τ) for homogeneous at r=0.
        G[i] = (4.0 / np.pi) * total
    return G


SIG_T, SIG_S, NU_SIG_F = np.array([1.0]), np.array([0.5]), np.array([0.75])
K_INF = 1.5


def build_v5_cylinder(geometry, r_nodes, r_wts, radii, sig_t, N):
    """V5: canonical Marshak for cylinder with Ki_{k+2} 3D reformulation."""
    R_cell = float(radii[-1])
    N_r = len(r_nodes)
    sig_t_n = np.array([sig_t[geometry.which_annulus(ri, radii)] for ri in r_nodes])
    rv = np.array([geometry.radial_volume_weight(rj) for rj in r_nodes])
    divisor = geometry.rank1_surface_divisor(R_cell)
    P = np.zeros((N, N_r))
    G = np.zeros((N_r, N))
    for n in range(N):
        P_esc_n = compute_P_cyl_v5(geometry, r_nodes, radii, sig_t, n,
                                   n_beta=32, dps=25)
        G_bc_n = compute_G_cyl_v5(geometry, r_nodes, radii, sig_t, n,
                                  n_beta=32, dps=25)
        P[n, :] = rv * r_wts * P_esc_n
        G[:, n] = sig_t_n * G_bc_n / divisor
    B = compute_B(N)
    return BoundaryClosureOperator(P=P, G=G, R=np.linalg.inv(B))


def _solve(build_fn, geom, R, N):
    def patched(*args, **kw):
        return build_fn(geom, args[1], args[2], args[3], args[4], N)
    orig = pg.build_closure_operator
    pg.build_closure_operator = patched
    try:
        sol = pg.solve_peierls_1g(
            geom, np.array([R]), SIG_T, SIG_S, NU_SIG_F,
            boundary="white",
            n_panels_per_region=2, p_order=5,
            n_angular=24, n_rho=24, n_surf_quad=24, dps=25,
            n_bc_modes=N,
        )
        return sol.k_eff
    finally:
        pg.build_closure_operator = orig


print("V5: CYLINDER canonical Marshak (Ki_{k+2} 3D, Jacobian, R=B^-1)")
print("="*70)
for R in [1.0, 2.0, 5.0, 10.0]:
    line = f"  R={R:5.1f}: "
    for N in [1, 2, 3, 5, 8]:
        k = _solve(build_v5_cylinder, CYLINDER_1D, R, N)
        err = abs(k - K_INF) / K_INF * 100
        line += f"N={N}:{err:6.2f}%  "
    print(line, flush=True)
