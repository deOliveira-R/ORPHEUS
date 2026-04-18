"""Phase A / C diagnostic: scan variants of P/G/R to find the right fix.

Variants tested:

  V1 (BASELINE): current code
      Mode 0: legacy P_esc, G_bc.  Modes n≥1: compute_*_mode with (ρ/R)² jacobian.
      R = diag(1, 3, 5, ...).

  V2 (CONSISTENT): use compute_*_mode (with Jacobian) for ALL modes n≥0.
      Question: does mode-0 with Jacobian give a better baseline, curing plateau?

  V3 (CYLINDER OBSERVER): replace cylinder compute_G_bc_mode with observer-
      centered form (Ki_{k+2} polynomial expansion) for modes n≥1.
      Question: does cylinder stop diverging at high N?

  V4 (CYLINDER JACOBIAN): cylinder compute_P_esc_mode uses (ρ_2D/R)² jacobian
      with Ki_{k+2} polynomial expansion instead of Ki_2 alone.

Expected outcomes per Part IV §29:
  Phase A (sphere plateau):  V2 should resolve the sphere plateau at N≥3
  Phase C (cylinder divergence): V3+V4 should stop cylinder divergence at high N
"""
from __future__ import annotations
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

import numpy as np
from orpheus.derivations import peierls_geometry as pg
from orpheus.derivations._kernels import _shifted_legendre_eval, ki_n_mp
from orpheus.derivations.peierls_geometry import (
    CYLINDER_1D, SPHERE_1D, gl_float, composite_gl_r,
    compute_P_esc, compute_G_bc,
    compute_P_esc_mode, compute_G_bc_mode,
    BoundaryClosureOperator, reflection_marshak,
)


SIG_T_ONE = np.array([1.0])
SIG_T, SIG_S, NU_SIG_F = SIG_T_ONE, np.array([0.5]), np.array([0.75])
K_INF = 1.5


def _shifted_legendre_coeffs(n_max: int) -> np.ndarray:
    """Return a (n_max+1, n_max+1) matrix of Taylor coefficients of
    shifted Legendre: row n is [p_{n,0}, p_{n,1}, ...] in
    P̃_n(x) = Σ_k p_{n,k} x^k."""
    coeffs = np.zeros((n_max + 1, n_max + 1))
    coeffs[0, 0] = 1.0  # P̃_0 = 1
    if n_max == 0:
        return coeffs
    coeffs[1, :] = [-1.0, 2.0, *([0.0] * (n_max - 1))]  # 2x - 1
    for n in range(1, n_max):
        # Bonnet recurrence for P_n(2x - 1):
        # (n+1) P̃_{n+1}(x) = (2n+1)(2x - 1) P̃_n(x) - n P̃_{n-1}(x)
        # Rewrite as polynomial-in-x recurrence.
        # Let y = 2x - 1, so P̃_n(x) = P_n(y).
        # (n+1) P_{n+1}(y) = (2n+1) y P_n(y) - n P_{n-1}(y)
        # In x coords: multiply by y = 2x - 1 = 2·(x-coeff_shift).
        # Easiest: build symbolically via polynomial multiplication.
        pass
    # Use numpy.polynomial as simpler:
    from numpy.polynomial import Polynomial
    coeffs = np.zeros((n_max + 1, n_max + 1))
    for n in range(n_max + 1):
        # Shifted Legendre P̃_n(x) = P_n(2x - 1).
        # Build the Legendre poly P_n and compose with 2x - 1.
        leg = np.zeros(n + 1)
        leg[n] = 1.0
        from numpy.polynomial import Legendre
        L_n = Legendre(leg)  # P_n in Legendre basis
        P_poly = L_n.convert(kind=Polynomial)  # P_n in monomial basis
        # Substitute y = 2x - 1.  If P_n(y) = sum c_k y^k, then
        # P̃_n(x) = sum c_k (2x - 1)^k.
        shifted = Polynomial([0.0])
        for k, c in enumerate(P_poly.coef):
            # (2x - 1)^k via Polynomial arithmetic.
            term = Polynomial([-1.0, 2.0]) ** k
            shifted = shifted + c * term
        coefs = shifted.coef
        coeffs[n, :len(coefs)] = coefs
    return coeffs


# Pre-compute up to n = 10
LEG_COEFFS = _shifted_legendre_coeffs(10)


# ═════════════════════════════════════════════════════════════════════
# Variants
# ═════════════════════════════════════════════════════════════════════

_ORIG_BUILD_CLOSURE = pg.build_closure_operator


def build_closure_V1_baseline(geometry, r_nodes, r_wts, radii, sig_t, N):
    """V1: current code = legacy mode 0 + compute_*_mode for n≥1."""
    return _ORIG_BUILD_CLOSURE(
        geometry, r_nodes, r_wts, radii, sig_t,
        n_angular=32, n_surf_quad=32, dps=25,
        n_bc_modes=N, reflection="marshak",
    )


def build_closure_V2_consistent(geometry, r_nodes, r_wts, radii, sig_t, N):
    """V2: use compute_*_mode for ALL modes n ≥ 0.

    Mode 0 thus gets the (ρ_max/R)² Jacobian as well. Rank-1 will
    NOT match legacy bit-exactly.
    """
    R_cell = float(radii[-1])
    N_r = len(r_nodes)
    sig_t_n = np.array([sig_t[geometry.which_annulus(ri, radii)]
                        for ri in r_nodes])
    rv = np.array([geometry.radial_volume_weight(rj) for rj in r_nodes])
    divisor = geometry.rank1_surface_divisor(R_cell)

    P = np.zeros((N, N_r))
    G = np.zeros((N_r, N))
    for n in range(N):
        P_esc_n = compute_P_esc_mode(
            geometry, r_nodes, radii, sig_t, n,
            n_angular=32, dps=25,
        )
        G_bc_n = compute_G_bc_mode(
            geometry, r_nodes, radii, sig_t, n,
            n_surf_quad=32, dps=25,
        )
        P[n, :] = rv * r_wts * P_esc_n
        G[:, n] = sig_t_n * G_bc_n / divisor
    return BoundaryClosureOperator(P=P, G=G, R=reflection_marshak(N))


def _compute_P_esc_mode_cyl_observer(
    geometry, r_nodes, radii, sig_t, n_mode, n_beta=32, dps=25,
):
    """V3/V4 cylinder observer-centered: Ki_{k+2} polynomial expansion
    with (ρ_2D/R)² Jacobian.
    """
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
                    r_i, cb, rho_2d,
                    np.asarray(radii), np.asarray(sig_t),
                )
            mu_2d = (rho_2d + r_i * cb) / R_cell
            jacobian = rho_2d * rho_2d * inv_R2
            sum_ki = 0.0
            for k in range(n_mode + 1):
                ki_kp2 = float(ki_n_mp(k + 2, float(tau), dps))
                sum_ki += p_n[k] * (mu_2d ** k) * ki_kp2
            total += beta_wts[k_b] * jacobian * sum_ki
        P[i] = (1.0 / np.pi) * total
    return P


def _compute_G_bc_mode_cyl_observer(
    geometry, r_nodes, radii, sig_t, n_mode, n_beta=32, dps=25,
):
    """Cylinder observer-centered G_bc: Ki_{k+2} expansion WITHOUT
    Jacobian (G_bc is a response, not a moment).
    """
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
                    r_i, cb, rho_2d,
                    np.asarray(radii), np.asarray(sig_t),
                )
            mu_2d = (rho_2d + r_i * cb) / R_cell
            sum_ki = 0.0
            for k in range(n_mode + 1):
                ki_kp2 = float(ki_n_mp(k + 2, float(tau), dps))
                sum_ki += p_n[k] * (mu_2d ** k) * ki_kp2
            total += beta_wts[k_b] * sum_ki
        # Normalization: at r=0 for n=0 this should give 2·Ki_2·something
        # to match legacy → tune to 4/π so G_bc_0(0) = 4·Ki_2
        # (the observer-centered analogue of sphere's 2·exp(-τ)).
        G[i] = (2.0 / np.pi) * total * 2.0  # prefactor 4/π per sphere analogy
    return G


def build_closure_V3_cyl_observer(geometry, r_nodes, r_wts, radii, sig_t, N):
    """V3 cylinder: observer-centered Ki_{k+2} expansion for n ≥ 1.

    Mode 0 stays legacy for bit-exact recovery.
    """
    if geometry.kind != "cylinder-1d":
        raise ValueError("V3 applies only to cylinder")

    R_cell = float(radii[-1])
    N_r = len(r_nodes)
    sig_t_n = np.array([sig_t[geometry.which_annulus(ri, radii)]
                        for ri in r_nodes])
    rv = np.array([geometry.radial_volume_weight(rj) for rj in r_nodes])
    divisor = geometry.rank1_surface_divisor(R_cell)

    P = np.zeros((N, N_r))
    G = np.zeros((N_r, N))

    # Mode 0: legacy
    P_esc_0 = compute_P_esc(
        geometry, r_nodes, radii, sig_t, n_angular=32, dps=25,
    )
    G_bc_0 = compute_G_bc(
        geometry, r_nodes, radii, sig_t, n_surf_quad=32, dps=25,
    )
    P[0, :] = rv * r_wts * P_esc_0
    G[:, 0] = sig_t_n * G_bc_0 / divisor

    # Modes n ≥ 1: observer-centered Ki_{k+2} polynomial expansion
    for n in range(1, N):
        P_esc_n = _compute_P_esc_mode_cyl_observer(
            geometry, r_nodes, radii, sig_t, n, n_beta=32, dps=25,
        )
        G_bc_n = _compute_G_bc_mode_cyl_observer(
            geometry, r_nodes, radii, sig_t, n, n_beta=32, dps=25,
        )
        P[n, :] = rv * r_wts * P_esc_n
        G[:, n] = sig_t_n * G_bc_n / divisor

    return BoundaryClosureOperator(P=P, G=G, R=reflection_marshak(N))


# ═════════════════════════════════════════════════════════════════════
# Solver via build_closure_operator replacement
# ═════════════════════════════════════════════════════════════════════

def _solve_with_closure(closure_fn, geometry, R, N):
    """Run solve_peierls_1g with a monkey-patched build_closure_operator."""
    orig = pg.build_closure_operator
    def patched(*args, **kwargs):
        return closure_fn(geometry, args[1], args[2], args[3], args[4], N)
    pg.build_closure_operator = patched
    try:
        sol = pg.solve_peierls_1g(
            geometry, np.array([R]), SIG_T, SIG_S, NU_SIG_F,
            boundary="white",
            n_panels_per_region=2, p_order=5,
            n_angular=32, n_rho=32, n_surf_quad=32, dps=25,
            n_bc_modes=N,
        )
        return sol.k_eff
    finally:
        pg.build_closure_operator = orig


def run_scan(name, closure_fn, geometries, R_values, N_values):
    print(f"\n{'='*70}")
    print(f"  {name}")
    print(f"{'='*70}")
    for geom_name, geom in geometries:
        if closure_fn is build_closure_V3_cyl_observer and geom_name != "cyl":
            continue
        print(f"\n  {geom_name}:")
        for R in R_values:
            line = f"    R={R:5.1f}: "
            for N in N_values:
                k = _solve_with_closure(closure_fn, geom, R, N)
                err = abs(k - K_INF) / K_INF * 100
                line += f"N={N}:{err:6.2f}%  "
            print(line, flush=True)


if __name__ == "__main__":
    geometries = [("sph", SPHERE_1D), ("cyl", CYLINDER_1D)]
    R_values = [1.0, 2.0, 5.0, 10.0]
    N_values = [1, 2, 3, 5, 8]

    run_scan("V1 BASELINE (current code)", build_closure_V1_baseline,
             geometries, R_values, N_values)

    run_scan("V2 CONSISTENT (compute_*_mode for mode 0 too)",
             build_closure_V2_consistent,
             geometries, R_values, N_values)

    # V3 is cylinder-only
    run_scan("V3 CYL OBSERVER (Ki_{k+2} polynomial expansion for n ≥ 1)",
             build_closure_V3_cyl_observer,
             geometries, R_values, N_values)
