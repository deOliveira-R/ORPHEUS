"""V6: corrected integrand with |μ_s| cosine weight (canonical cosine moment)."""
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")
import numpy as np
from orpheus.derivations import peierls_geometry as pg
from orpheus.derivations._kernels import _shifted_legendre_eval
from orpheus.derivations.peierls_geometry import (
    SPHERE_1D, CYLINDER_1D, gl_float,
    compute_P_esc, compute_G_bc, compute_G_bc_mode,
    BoundaryClosureOperator,
)


def compute_P_esc_mode_v6(geometry, r_nodes, radii, sig_t, n_mode,
                          n_angular=32, dps=25):
    """Corrected: integrand = W_Ω · μ_exit · P̃_n(μ_exit) · K_esc (cosine-weighted)."""
    R_cell = float(radii[-1])
    omega_pts, omega_wts = gl_float(n_angular, 0.0, np.pi, dps)
    cos_omegas = np.cos(omega_pts)
    angular_factor = geometry.angular_weight(omega_pts)
    pref = geometry.prefactor
    N_r = len(r_nodes)
    P = np.zeros(N_r)
    for i in range(N_r):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_max_val = geometry.rho_max(r_i, cos_om, R_cell)
            if rho_max_val <= 0.0:
                continue
            tau = geometry.optical_depth_along_ray(
                r_i, cos_om, rho_max_val,
                np.asarray(radii), np.asarray(sig_t))
            K_esc = geometry.escape_kernel_mp(tau, dps)
            mu_exit = (rho_max_val + r_i * cos_om) / R_cell
            p_tilde = float(_shifted_legendre_eval(n_mode, np.array([mu_exit]))[0])
            # Cosine-weighted integrand (canonical Marshak moment)
            total += omega_wts[k] * angular_factor[k] * mu_exit * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_B(N):
    mu, w = gl_float(40, 0.0, 1.0, dps=30)
    B = np.zeros((N, N))
    for n in range(N):
        pn = _shifted_legendre_eval(n, mu)
        for m in range(N):
            pm = _shifted_legendre_eval(m, mu)
            B[n, m] = float(np.sum(w * mu * pn * pm))
    return B


def build_v6(geometry, r_nodes, r_wts, radii, sig_t, N):
    """V6: canonical with μ_s cosine weight inside integrand + R = B^{-1}."""
    R_cell = float(radii[-1])
    N_r = len(r_nodes)
    sig_t_n = np.array([sig_t[geometry.which_annulus(ri, radii)] for ri in r_nodes])
    rv = np.array([geometry.radial_volume_weight(rj) for rj in r_nodes])
    divisor = geometry.rank1_surface_divisor(R_cell)
    P = np.zeros((N, N_r))
    G = np.zeros((N_r, N))
    for n in range(N):
        P_esc_n = compute_P_esc_mode_v6(geometry, r_nodes, radii, sig_t, n,
                                        n_angular=32, dps=25)
        G_bc_n = compute_G_bc_mode(geometry, r_nodes, radii, sig_t, n,
                                    n_surf_quad=32, dps=25)
        P[n, :] = rv * r_wts * P_esc_n
        G[:, n] = sig_t_n * G_bc_n / divisor
    B = compute_B(N)
    return BoundaryClosureOperator(P=P, G=G, R=np.linalg.inv(B))


SIG_T, SIG_S, NU_SIG_F = np.array([1.0]), np.array([0.5]), np.array([0.75])
K_INF = 1.5


def _solve(geom, R, N):
    def patched(*args, **kw):
        return build_v6(geom, args[1], args[2], args[3], args[4], N)
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


print("V6: canonical Marshak with |μ_s| cosine weight + R = B^{-1}")
print("="*70)
for geom_name, geom in [("SPHERE", SPHERE_1D), ("CYLINDER", CYLINDER_1D)]:
    print(f"\n{geom_name}:")
    for R in [1.0, 2.0, 5.0, 10.0]:
        line = f"  R={R:5.1f}: "
        for N in [1, 2, 3, 5, 8]:
            k = _solve(geom, R, N)
            err = abs(k - K_INF) / K_INF * 100
            line += f"N={N}:{err:6.2f}%  "
        print(line, flush=True)
