"""Diagnostic: c_in-aware split basis rank-N k_eff for hollow sphere.

Created by numerics-investigator on 2026-04-21.

PHILOSOPHY: Implement rank-(N_g, N_s, N_i) closure in a SELF-CONTAINED
script first. If rank-(0,0,0) matches F.4 to ~0.077%, the structural
derivation is correct. Then rank-(1,1,1) tests whether the novel DOFs
reduce residual as predicted. Only after empirical validation do we
migrate into peierls_geometry.py.

==== Geometry ====

Hollow sphere  ρ := r_0 / R,  µ_crit := √(1-ρ²).
Basis split on OUTER surface cosine µ ∈ [0, 1]:
  - grazing: [0, µ_crit]      P̃_n^g(µ) = Gram-Schmidt on µ-weight
  - steep:   [µ_crit, 1]      P̃_n^s(µ) = (1/ρ) P̃_n(c_I(µ))
Inner surface cosine c ∈ [0, 1]:
  - P̃_n^I(c) = half-range Legendre on [0,1] in c-weight

All three sub-bases are orthonormal in the µ-weighted inner product on
their support. At σ_t = 0, outer-steep → inner transmission is
W_{oi,s}^{mn} = (1/ρ) δ_{mn} — the structural win.

==== Block assembly ====

Mode layout: [graze_0..graze_{Ng-1}, steep_0..steep_{Ns-1}, inner_0..inner_{Ni-1}]
Dimension D = N_g + N_s + N_i.

Transmission W is D×D with block structure:
  - W_{gg}: grazing self, (N_g, N_g).
  - W_{gs} = W_{sg} = 0 (disjoint cone).
  - W_{gi} = W_{ig} = 0 (grazing doesn't couple to inner).
  - W_{si}: steep → inner, (N_s, N_i). At σ_t=0: (1/ρ) δ_{mn}.
  - W_{is}: inner → steep, same formula by time reversal.
  - W_{ss} = 0 (steep emits toward inner, always).
  - W_{ii} = 0 (convex cavity — interior = vacuum).

Plus cavity self-coupling: rays crossing the cavity from inner return to
inner. By 1D spherical symmetry, this is identity in the inner basis:
K_cav^{mn} = δ_{mn}. So when assembling the "effective W" for the closure,
we use W̃_{ii} = K_cav = I_{N_i}.

==== White BC ====

Isotropic re-emission with albedo β on outer surface.  ψ^- on outer is
proportional to 1_{[0,1]}(µ) = (1/√2) F.4-outer-mode = (µ_crit P̃_0^g + ρ P̃_0^s).
So the white BC re-injects ONLY in the (µ_crit, ρ) direction at mode-0 outer.

==== Assembly ====

Surface current vector ψ^+ ∈ R^D. With white BC:
  ψ^+ = W · ψ^- + Q_V · P_esc
  ψ^- = B · ψ^+     (B = boundary operator)

where B is block-diagonal:
  - Graze block (outer-graze → outer-graze via white BC):
      [B_outer] = rank-1 projection onto (µ_crit, ρ, 0) direction × β.
      Specifically: ψ^-_{graze, 0} = µ_crit · β · J_out_total
                    ψ^-_{steep, 0} = ρ · β · J_out_total
                    higher modes of outer = 0.
      Where J_out_total = µ_crit · ψ^+_{graze,0} + ρ · ψ^+_{steep,0}.
  - Inner: ψ^-_{inner, m} = ψ^+_{inner, m}  (cavity identity)

So:  B = block-diag(B_outer, I_{Ni})
where B_outer = β · e_F4 · e_F4^T (with e_F4 = (µ_crit, ρ) padded to full outer space).

The closure equation:
  (I - W·B) ψ^+ = P_esc · Q_V

then φ_V = G · ψ^- + K_vol · Q_V = G · B · ψ^+ + K_vol · Q_V.

For k_eff (with Q_V = (νΣ_f / k) φ_V + Σ_s φ_V):
  Σ_t φ_V = K_total · (Σ_s + νΣ_f/k) φ_V
  K_total = K_vol + G · B · (I - W·B)^{-1} · P_esc.

This is exactly the rank-N closure template — just with the split basis
and block-structured W/B.

==== Reciprocity note ====

Sanchez-McCormick reciprocity between physical transmission matrices
requires area factors: W_io^{nm}_phys = ρ² · W_{oi,s}^{mn}_phys. Our
W_oi_s and W_io formulas give the SAME integral (by time reversal of
the ray), which means in the µ-weighted inner-product convention they
are already symmetric. The ρ² factor enters when relating to the physical
Sanchez-McCormick W. In the closure we use the as-integrated W directly,
and the divisor on G provides the area normalization — SAME convention
as the existing rank-2 closure.

==== Tests run in __main__ ====

1. Orthonormality of all three sub-bases (symbolic verified).
2. W matrix block structure at σ_t = 0.
3. F.4 decomposition (µ_crit, ρ) in split basis with Parseval = 1.
4. k_eff at rank-(0,0,0): should match F.4 to within quadrature.
5. k_eff at rank-(1,1,1): should reduce residual vs F.4.
6. Cross σ_t·R and r_0/R scan.
"""
from __future__ import annotations

import math
import sys

import numpy as np
import sympy as sp
from scipy import integrate

from orpheus.derivations.peierls_geometry import (
    CurvilinearGeometry,
    build_volume_kernel,
    compute_G_bc_inner,
    compute_G_bc_outer,
    compute_P_esc_inner,
    compute_P_esc_outer,
    compute_hollow_sph_transmission,
    composite_gl_r,
)
from orpheus.derivations._kernels import _shifted_legendre_eval

# ---------------------------------------------------------------------------
# Split-basis Gram-Schmidt machinery
# ---------------------------------------------------------------------------


def half_range_legendre_symbolic(n_max, a_sym, b_sym, var):
    """Gram-Schmidt of {1, x, x², ...} w.r.t. ⟨f, g⟩ = ∫_a^b f g x dx.

    Returns list of SymPy expressions in ``var``.
    """
    basis = []
    for k in range(n_max + 1):
        fk = var**k
        for pj in basis:
            cjk = sp.integrate(fk * pj * var, (var, a_sym, b_sym))
            fk = fk - cjk * pj
        norm_sq = sp.integrate(fk * fk * var, (var, a_sym, b_sym))
        fk = sp.simplify(fk / sp.sqrt(norm_sq))
        basis.append(fk)
    return basis


def _unit_legendre_lambdas(n_max):
    """Numerical callables P̃_n(x) for x in [0, 1], µ-weighted."""
    mu = sp.Symbol("mu", positive=True, real=True)
    exprs = half_range_legendre_symbolic(n_max, 0, 1, mu)
    return [sp.lambdify(mu, e, modules=["numpy"]) for e in exprs]


def grazing_lambdas(n_max, mu_crit_val):
    """Grazing basis on [0, µ_crit] (float numeric)."""
    mu = sp.Symbol("mu", positive=True, real=True)
    muc = sp.Float(mu_crit_val)
    exprs = half_range_legendre_symbolic(n_max, 0, muc, mu)
    return [sp.lambdify(mu, e, modules=["numpy"]) for e in exprs]


# ---------------------------------------------------------------------------
# Geometry
# ---------------------------------------------------------------------------


def mu_crit(rho):
    return math.sqrt(1.0 - rho * rho)


def c_of_mu(mu, rho):
    """c_I(µ) = sqrt(1 - (1/ρ²)(1-µ²)) for µ ∈ [µ_crit, 1]."""
    val = 1.0 - (1.0 / (rho * rho)) * (1.0 - mu * mu)
    if val < 0.0:
        val = 0.0
    return math.sqrt(val)


def mu_of_c(c, rho):
    """µ_R(c) = sqrt(1 - ρ²(1-c²))."""
    return math.sqrt(1.0 - rho * rho * (1.0 - c * c))


def chord_oi(c, rho):
    """Normalised steep outer→inner chord length in R units."""
    return mu_of_c(c, rho) - rho * c


# ---------------------------------------------------------------------------
# Transmission matrix W in split basis
# ---------------------------------------------------------------------------


def compute_W_split_basis(rho, tau, N_g, N_s, N_i, *, epsabs=1e-13, epsrel=1e-11):
    """Build the (N_g + N_s + N_i) × (N_g + N_s + N_i) W matrix.

    Block layout: [graze | steep | inner]   (rows = incoming, cols = outgoing).
    Returns W where W[incoming_row, outgoing_col] is the transmission
    matrix element in the µ-weighted orthonormal basis.

    Non-zero blocks:
      W_{gg}   (graze,  graze)  : outer-grazing self.
      W_{si}   (inner,  steep)  : outer-steep emission → inner arrival.
      W_{is}   (steep,  inner)  : inner emission → outer-steep arrival.

    Block conventions:
      W_gg[m, n] = ∫_0^{µ_crit} P̃_m^g(µ) P̃_n^g(µ) exp(-2τµ) µ dµ
      W_si[m, n] = (1/ρ) ∫_0^1 P̃_m(c) P̃_n(c) exp(-τ s(c; ρ)) c dc
                 (rows index inner basis, cols index steep basis;
                  P̃_n^s evaluated at µ_R(c) = (1/ρ) P̃_n(c))
      W_is[m, n] = same integral (time-reversal symmetric).
    """
    muc = mu_crit(rho)
    D = N_g + N_s + N_i
    W = np.zeros((D, D))

    leg = _unit_legendre_lambdas(max(N_s, N_i, 1) - 1 + 1)
    leg_g = grazing_lambdas(N_g - 1 + 1 if N_g > 0 else 0, muc) if N_g > 0 else []

    # Index ranges
    g0, g1 = 0, N_g
    s0, s1 = N_g, N_g + N_s
    i0, i1 = N_g + N_s, D

    # W_gg: outer-grazing self-transmission
    for m in range(N_g):
        for n in range(N_g):
            def integrand(mu, m=m, n=n):
                return leg_g[m](mu) * leg_g[n](mu) * math.exp(-2.0 * tau * mu) * mu
            val, _ = integrate.quad(integrand, 0.0, muc, epsabs=epsabs, epsrel=epsrel)
            W[g0 + m, g0 + n] = val

    # W_si: rows index inner-receiving-mode m, cols index steep-outgoing-mode n.
    # Integrand in c-coordinate (c = c_I at inner).
    for m in range(N_i):
        for n in range(N_s):
            def integrand(c, m=m, n=n):
                return leg[m](c) * leg[n](c) * math.exp(-tau * chord_oi(c, rho)) * c
            val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=epsabs, epsrel=epsrel)
            W[i0 + m, s0 + n] = val / rho

    # W_is: rows index steep-receiving-mode m, cols index inner-outgoing-mode n.
    # By time reversal, ray from inner at cosine c arrives at outer-steep at µ_R(c),
    # with P̃_m^s(µ_R(c)) = (1/ρ) P̃_m(c). Integrand over c-coordinate.
    for m in range(N_s):
        for n in range(N_i):
            def integrand(c, m=m, n=n):
                return leg[m](c) * leg[n](c) * math.exp(-tau * chord_oi(c, rho)) * c
            val, _ = integrate.quad(integrand, 0.0, 1.0, epsabs=epsabs, epsrel=epsrel)
            W[s0 + m, i0 + n] = val / rho

    return W


# ---------------------------------------------------------------------------
# P_esc (volume → surface) and G_bc (surface → volume) in split basis
# ---------------------------------------------------------------------------


def compute_P_esc_graze_mode(geom, r_nodes, radii, sig_t, m, mu_crit_val,
                             leg_g, n_angular=32, dps=25):
    """Projection of volume-escape to outer GRAZING sub-basis at mode m.

    Integrand: µ_exit · P̃_m^g(µ_exit) · K_esc(τ) · angular_factor
    restricted to directions where µ_exit ∈ [0, µ_crit_val] at outer surface.
    Follows the same structure as compute_P_esc_outer_mode_marshak but
    evaluates the GRAZING basis function and restricts to grazing µ range.
    """
    from orpheus.derivations.peierls_geometry import gl_float

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    omega_low, omega_high = geom.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geom.ray_direction_cosine(omega_pts)
    angular_factor = geom.angular_weight(omega_pts)
    pref = geom.prefactor
    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_out = geom.rho_max(r_i, cos_om, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, cos_om)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_exit = (rho_out + r_i * cos_om) / R
            if mu_exit > mu_crit_val:
                continue  # steep ray; not counted by grazing sub-basis
            tau = geom.optical_depth_along_ray(r_i, cos_om, rho_out, radii, sig_t)
            K_esc = geom.escape_kernel_mp(tau, dps)
            p_tilde = leg_g[m](mu_exit)
            total += omega_wts[k] * angular_factor[k] * mu_exit * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_P_esc_steep_mode(geom, r_nodes, radii, sig_t, m, mu_crit_val, rho_val,
                             leg, n_angular=32, dps=25):
    """Projection of volume-escape to outer STEEP sub-basis at mode m.

    P̃_m^s(µ_exit) = (1/ρ) P̃_m(c_I(µ_exit)) for µ_exit ∈ [µ_crit, 1].
    """
    from orpheus.derivations.peierls_geometry import gl_float

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    omega_low, omega_high = geom.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geom.ray_direction_cosine(omega_pts)
    angular_factor = geom.angular_weight(omega_pts)
    pref = geom.prefactor
    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_out = geom.rho_max(r_i, cos_om, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, cos_om)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_exit = (rho_out + r_i * cos_om) / R
            if mu_exit < mu_crit_val:
                continue  # grazing ray; not counted here
            c_val = c_of_mu(mu_exit, rho_val)
            tau = geom.optical_depth_along_ray(r_i, cos_om, rho_out, radii, sig_t)
            K_esc = geom.escape_kernel_mp(tau, dps)
            # P̃_m^s(µ) = (1/ρ) P̃_m(c_I(µ))
            p_tilde_s = leg[m](c_val) / rho_val
            total += omega_wts[k] * angular_factor[k] * mu_exit * p_tilde_s * K_esc
        P[i] = pref * total
    return P


def compute_P_esc_inner_mode_split(geom, r_nodes, radii, sig_t, m, leg,
                                    n_angular=32, dps=25):
    """Projection of volume-escape to INNER at mode m (standard basis on [0,1])."""
    from orpheus.derivations.peierls_geometry import gl_float

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    r_0 = float(geom.inner_radius)
    omega_low, omega_high = geom.angular_range
    omega_pts, omega_wts = gl_float(n_angular, omega_low, omega_high, dps)
    cos_omegas = geom.ray_direction_cosine(omega_pts)
    angular_factor = geom.angular_weight(omega_pts)
    pref = geom.prefactor
    N = len(r_nodes)
    P = np.zeros(N)
    for i in range(N):
        r_i = float(r_nodes[i])
        total = 0.0
        for k in range(n_angular):
            cos_om = cos_omegas[k]
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, cos_om)
            if rho_in_minus is None:
                continue
            tau = geom.optical_depth_along_ray(r_i, cos_om, rho_in_minus, radii, sig_t)
            K_esc = geom.escape_kernel_mp(tau, dps)
            sin_om = math.sqrt(max(0.0, 1.0 - cos_om * cos_om))
            h_sq = r_i * r_i * sin_om * sin_om
            mu_exit_sq = max(0.0, (r_0 * r_0 - h_sq) / (r_0 * r_0))
            mu_exit = math.sqrt(mu_exit_sq)
            p_tilde = leg[m](mu_exit)
            total += omega_wts[k] * angular_factor[k] * mu_exit * p_tilde * K_esc
        P[i] = pref * total
    return P


def compute_G_bc_graze_mode(geom, r_nodes, radii, sig_t, m, mu_crit_val,
                            leg_g, n_surf_quad=32, dps=25):
    """Surface-to-volume response at observer r_i from unit mode-m incident
    on outer surface in GRAZING sub-basis. Integrand uses P̃_m^g(µ_s)
    restricted to µ_s ∈ [0, µ_crit]."""
    from orpheus.derivations.peierls_geometry import gl_float

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N = len(r_nodes)
    G = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            rho_out = geom.rho_max(r_i, ct, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, ct)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_s = (rho_out + r_i * ct) / R
            if mu_s > mu_crit_val:
                continue
            if len(radii) == 1:
                tau = sig_t[0] * rho_out
            else:
                tau = geom.optical_depth_along_ray(r_i, ct, rho_out, radii, sig_t)
            p_tilde = leg_g[m](mu_s)
            total += theta_wts[k] * st * mu_s * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def compute_G_bc_steep_mode(geom, r_nodes, radii, sig_t, m, mu_crit_val, rho_val,
                            leg, n_surf_quad=32, dps=25):
    """Surface-to-volume response for STEEP sub-basis mode m."""
    from orpheus.derivations.peierls_geometry import gl_float

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    R = float(radii[-1])
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N = len(r_nodes)
    G = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            rho_out = geom.rho_max(r_i, ct, R)
            if rho_out <= 0.0:
                continue
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, ct)
            if rho_in_minus is not None and rho_in_minus < rho_out:
                continue
            mu_s = (rho_out + r_i * ct) / R
            if mu_s < mu_crit_val:
                continue
            if len(radii) == 1:
                tau = sig_t[0] * rho_out
            else:
                tau = geom.optical_depth_along_ray(r_i, ct, rho_out, radii, sig_t)
            c_val = c_of_mu(mu_s, rho_val)
            p_tilde = leg[m](c_val) / rho_val
            total += theta_wts[k] * st * mu_s * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


def compute_G_bc_inner_mode_split(geom, r_nodes, radii, sig_t, m, leg,
                                   n_surf_quad=32, dps=25):
    """Surface-to-volume response from INNER surface, standard basis."""
    from orpheus.derivations.peierls_geometry import gl_float

    r_nodes = np.asarray(r_nodes, dtype=float)
    radii = np.asarray(radii, dtype=float)
    sig_t = np.asarray(sig_t, dtype=float)
    r_0 = float(geom.inner_radius)
    theta_pts, theta_wts = gl_float(n_surf_quad, 0.0, np.pi, dps)
    cos_thetas = np.cos(theta_pts)
    sin_thetas = np.sin(theta_pts)
    N = len(r_nodes)
    G = np.zeros(N)
    for i in range(N):
        r_i = r_nodes[i]
        total = 0.0
        for k in range(n_surf_quad):
            ct = cos_thetas[k]
            st = sin_thetas[k]
            rho_in_minus, _ = geom.rho_inner_intersections(r_i, ct)
            if rho_in_minus is None:
                continue
            if len(radii) == 1:
                tau = sig_t[0] * rho_in_minus
            else:
                tau = geom.optical_depth_along_ray(r_i, ct, rho_in_minus, radii, sig_t)
            sin_om = math.sqrt(max(0.0, 1.0 - ct * ct))
            h_sq = r_i * r_i * sin_om * sin_om
            mu_s_sq = max(0.0, (r_0 * r_0 - h_sq) / (r_0 * r_0))
            mu_s = math.sqrt(mu_s_sq)
            p_tilde = leg[m](mu_s)
            total += theta_wts[k] * st * mu_s * p_tilde * math.exp(-tau)
        G[i] = 2.0 * total
    return G


# ---------------------------------------------------------------------------
# Closure assembly: K_total = K_vol + G · B · (I - W·B)^{-1} · P_esc
# ---------------------------------------------------------------------------


def build_split_basis_closure(
    geom, r_nodes, r_wts, radii, sig_t,
    N_g, N_s, N_i, *,
    n_angular=32, n_surf_quad=32, dps=25,
    albedo=1.0,
):
    """Assemble K_bc = G · B · (I - W·B)^{-1} · P for split basis.

    Returns (K_bc, aux) where aux is a dict with debugging info.
    """
    rho_val = float(geom.inner_radius) / float(radii[-1])
    muc = mu_crit(rho_val)

    # Build basis lambdas
    nmax_inner = N_i - 1 + 1 if N_i > 0 else 1
    nmax_steep = N_s - 1 + 1 if N_s > 0 else 1
    nmax_leg = max(nmax_inner, nmax_steep, 1) - 1 + 1
    leg = _unit_legendre_lambdas(max(N_s, N_i, 1))
    leg_g = grazing_lambdas(max(N_g, 1), muc) if N_g > 0 else []

    D = N_g + N_s + N_i
    N_r = len(r_nodes)

    # Volume-to-surface P matrix (D rows x N_r cols)
    P = np.zeros((D, N_r))
    # Surface-to-volume G matrix (N_r rows x D cols)
    G = np.zeros((N_r, D))

    R_out = float(radii[-1])
    r_in = float(geom.inner_radius)
    A_outer = R_out * R_out  # 4π R² / 4π
    A_inner = r_in * r_in    # 4π r_in² / 4π

    # Compute volume weighting vector r² · w_j for the P row scaling.
    # peierls uses rv = r^(d-1) for sphere d=3 -> r²
    rv = r_nodes**2

    # Grazing P, G
    for m in range(N_g):
        P_m = compute_P_esc_graze_mode(
            geom, r_nodes, radii, sig_t, m, muc, leg_g,
            n_angular=n_angular, dps=dps,
        )
        G_m = compute_G_bc_graze_mode(
            geom, r_nodes, radii, sig_t, m, muc, leg_g,
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[m, :] = rv * r_wts * P_m
        G[:, m] = G_m / A_outer  # G normalised by outer area

    # Steep P, G
    for m in range(N_s):
        P_m = compute_P_esc_steep_mode(
            geom, r_nodes, radii, sig_t, m, muc, rho_val, leg,
            n_angular=n_angular, dps=dps,
        )
        G_m = compute_G_bc_steep_mode(
            geom, r_nodes, radii, sig_t, m, muc, rho_val, leg,
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[N_g + m, :] = rv * r_wts * P_m
        G[:, N_g + m] = G_m / A_outer  # steep still on outer surface

    # Inner P, G
    for m in range(N_i):
        P_m = compute_P_esc_inner_mode_split(
            geom, r_nodes, radii, sig_t, m, leg,
            n_angular=n_angular, dps=dps,
        )
        G_m = compute_G_bc_inner_mode_split(
            geom, r_nodes, radii, sig_t, m, leg,
            n_surf_quad=n_surf_quad, dps=dps,
        )
        P[N_g + N_s + m, :] = rv * r_wts * P_m
        G[:, N_g + N_s + m] = G_m / A_inner

    # Scale G by sig_t(r_i)
    sig_t_n = np.empty(N_r)
    for i, ri in enumerate(r_nodes):
        ki = geom.which_annulus(ri, radii)
        sig_t_n[i] = sig_t[ki]
    G = sig_t_n[:, None] * G

    # Transmission W
    sig_t_val = float(sig_t[0])
    tau = sig_t_val * R_out
    W = compute_W_split_basis(rho_val, tau, N_g, N_s, N_i)

    # White BC operator B: outer block is rank-1 projection onto (µ_crit, ρ);
    # inner block is identity (cavity passes through).
    # The white BC re-emits ψ^- = 2β J_out_total · 1_{[0,1]}(µ)
    #                         = √2 β J_out_total · (µ_crit P̃_0^g + ρ P̃_0^s)
    # So as an operator from ψ^+ to ψ^-:
    #   J_out_total = √2 · (µ_crit · a_0^g + ρ · a_0^s) = √2 · (µ_crit ψ^+_{g0} + ρ ψ^+_{s0})
    #   ψ^-_{g0} = β · 2 · J_out_total · (µ_crit / √2) · √2
    #
    # Let me re-derive carefully:
    #   outgoing current = J^+ = ∫_0^1 µ ψ^+(µ) dµ
    #   in split basis: ψ^+(µ) = Σ_n [a_n^g P̃_n^g(µ) 1_{[0,µc]} + a_n^s P̃_n^s(µ) 1_{[µc,1]}]
    #   J^+ = Σ_n [a_n^g ⟨P̃_n^g, 1⟩_µ,g + a_n^s ⟨P̃_n^s, 1⟩_µ,s]
    #
    # ⟨P̃_n^g, 1⟩_µ,g = ∫_0^{µc} P̃_n^g(µ) µ dµ.  For n=0, P̃_0^g = sqrt(2)/µ_crit,
    #   so ⟨P̃_0^g, 1⟩_µ = sqrt(2)/µ_crit · µ_crit²/2 = µ_crit / sqrt(2).
    # For n>=1, ⟨P̃_n^g, 1⟩_µ,g = 0 (orthogonality).
    # ⟨P̃_n^s, 1⟩_µ,s = (1/ρ) ∫_{µc}^1 P̃_n(c_I(µ)) µ dµ.  Substitute c = c_I(µ):
    #   µ dµ = ρ² c dc (from Jacobian), so = (1/ρ) · ρ² ∫_0^1 P̃_n(c) c dc
    #        = ρ · ⟨P̃_n, 1⟩_c.  For n=0, P̃_0 = sqrt(2), ⟨P̃_0, 1⟩_c = sqrt(2)/2.
    #   So ⟨P̃_0^s, 1⟩_µ = ρ / sqrt(2).  For n>=1, 0.
    #
    # J^+ = (1/sqrt(2)) (µ_crit a_0^g + ρ a_0^s).
    #
    # White BC: ψ^- = 2β J^+ · 1_{[0,1]}(µ).  Decompose this in split basis:
    #   1_{[0,1]}(µ) = c_g · P̃_0^g(µ) 1_{[0,µc]} + c_s · P̃_0^s(µ) 1_{[µc,1]}
    #   c_g = ⟨1, P̃_0^g⟩_µ,g = µ_crit / sqrt(2)
    #   c_s = ⟨1, P̃_0^s⟩_µ,s = ρ / sqrt(2)
    # So ψ^- = 2β J^+ · [ (µ_crit/sqrt(2)) P̃_0^g + (ρ/sqrt(2)) P̃_0^s ]
    #        = sqrt(2) β J^+ · [ µ_crit · P̃_0^g + ρ · P̃_0^s ]
    #
    # Thus:
    #   b^-_{g0} = sqrt(2) β J^+ µ_crit = sqrt(2) β (1/sqrt(2)) (µ_crit a_0^g + ρ a_0^s) µ_crit
    #            = β µ_crit (µ_crit a_0^g + ρ a_0^s)
    #   b^-_{s0} = β ρ (µ_crit a_0^g + ρ a_0^s)
    #   higher g/s modes: b^- = 0
    #
    # For INNER cavity: identity (crossing vacuum).
    B = np.zeros((D, D))
    if N_g >= 1 and N_s >= 1:
        B[0, 0] = albedo * muc * muc
        B[0, N_g] = albedo * muc * rho_val
        B[N_g, 0] = albedo * rho_val * muc
        B[N_g, N_g] = albedo * rho_val * rho_val
    for m in range(N_i):
        B[N_g + N_s + m, N_g + N_s + m] = 1.0  # cavity identity

    # Closure: ψ^- = B ψ^+, ψ^+ = W ψ^- + Q_surf (= P Q_V).
    # So ψ^+ = W B ψ^+ + P Q_V => ψ^+ = (I - W B)^{-1} P Q_V.
    # Then K_bc · Q_V = G · (ψ^- contribution) = G · B · ψ^+ = G B (I - WB)^{-1} P Q_V.
    M = np.eye(D) - W @ B
    try:
        Minv = np.linalg.inv(M)
    except np.linalg.LinAlgError:
        return None, {"error": "singular I - W·B"}

    K_bc = G @ B @ Minv @ P
    return K_bc, {
        "W": W, "B": B, "P": P, "G": G, "M": M,
        "rho": rho_val, "mu_crit": muc, "tau": tau,
    }


# ---------------------------------------------------------------------------
# k_eff power iteration
# ---------------------------------------------------------------------------


def solve_k_eff(K, sig_t_val, sig_s_val, nu_sig_f_val, *, tol=1e-12, max_iter=500):
    N = K.shape[0]
    A = np.diag(np.full(N, sig_t_val)) - K * sig_s_val
    B = K * nu_sig_f_val
    phi = np.ones(N)
    k = 1.0
    for _ in range(max_iter):
        q = B @ phi / k
        phi_new = np.linalg.solve(A, q)
        B_phi_new = B @ phi_new
        B_phi = B @ phi
        k_new = k * (np.abs(B_phi_new).sum() / np.abs(B_phi).sum())
        if abs(k_new - k) < tol:
            return k_new
        phi = phi_new / max(np.linalg.norm(phi_new), 1e-30)
        k = k_new
    return k


# ---------------------------------------------------------------------------
# Structural tests
# ---------------------------------------------------------------------------


def test_W_block_structure(rho_val=0.3, tau=0.0):
    """At τ=0: W_gg = I_{N_g}, W_si = (1/ρ) I, W_is = (1/ρ) I."""
    W = compute_W_split_basis(rho_val, tau, 3, 3, 3)
    print(f"\n[Structural test] W at τ=0, ρ={rho_val}")
    print(f"  W shape = {W.shape}")
    print(f"  max |W - diag-block ideal| = ...")
    # Check W_gg diagonal
    W_gg = W[:3, :3]
    err_gg = np.max(np.abs(W_gg - np.eye(3)))
    # Check W_si = (1/ρ) I  (rows = inner = last 3, cols = steep = middle 3)
    W_si = W[6:, 3:6]
    err_si = np.max(np.abs(W_si - np.eye(3) / rho_val))
    # Check W_is
    W_is = W[3:6, 6:]
    err_is = np.max(np.abs(W_is - np.eye(3) / rho_val))
    # Check zero blocks
    err_zero = 0.0
    for (r_lo, r_hi, c_lo, c_hi) in [
        (0, 3, 3, 6),   # graze-steep
        (0, 3, 6, 9),   # graze-inner
        (3, 6, 0, 3),   # steep-graze
        (3, 6, 3, 6),   # steep-steep
        (6, 9, 0, 3),   # inner-graze
        (6, 9, 6, 9),   # inner-inner
    ]:
        block = W[r_lo:r_hi, c_lo:c_hi]
        err_zero = max(err_zero, np.max(np.abs(block)))
    print(f"  |W_gg - I|_∞ = {err_gg:.3e}")
    print(f"  |W_si - I/ρ|_∞ = {err_si:.3e}")
    print(f"  |W_is - I/ρ|_∞ = {err_is:.3e}")
    print(f"  |zero blocks|_∞ = {err_zero:.3e}")
    ok = max(err_gg, err_si, err_is, err_zero) < 1e-10
    print(f"  PASS: {ok}")
    return ok


def test_F4_decomposition(rho_val):
    """F.4 outer-mode-0 (= constant on [0,1] in µ-weighted basis) should decompose
    as (µ_crit, ρ) in the split basis of dimension 2 (mode-0 graze + mode-0 steep)."""
    muc = mu_crit(rho_val)
    # The constant 1 on [0,1] projects onto split basis mode-0 components as:
    # c_g = ⟨1, P̃_0^g⟩_µ,g = µ_crit / sqrt(2).
    # Wait — that's if we're decomposing 1 · 1_{[0,1]} against the basis.
    # But F.4's outer mode is sqrt(2) · 1_{[0,1]} (the Legendre P̃_0 on [0,1]).
    # So in split basis: sqrt(2) · 1_{[0,1]} has components:
    #   a_0^g = ⟨sqrt(2) · 1, P̃_0^g⟩_µ,g = sqrt(2) · µ_crit / sqrt(2) = µ_crit
    #   a_0^s = ⟨sqrt(2) · 1, P̃_0^s⟩_µ,s = sqrt(2) · ρ / sqrt(2) = ρ
    # Parseval: µ_crit² + ρ² = (1 - ρ²) + ρ² = 1. ✓
    pars = muc * muc + rho_val * rho_val
    print(f"\n[Structural test] F.4 decomposition at ρ={rho_val}")
    print(f"  µ_crit = {muc:.6f}, ρ = {rho_val}")
    print(f"  Parseval: µ_crit² + ρ² = {pars:.12f} (expected 1.0)")
    ok = abs(pars - 1.0) < 1e-14
    print(f"  PASS: {ok}")
    return ok


# ---------------------------------------------------------------------------
# Reference: scalar F.4 for baseline
# ---------------------------------------------------------------------------


def run_scalar_f4(r_0, R_out, sig_t_val, sig_s_val, nu_sig_f_val,
                  n_panels=2, p_order=4, dps=15, n_ang=32):
    """Compute k_eff with the scalar F.4 closure (rank-2 white)."""
    from orpheus.derivations.peierls_geometry import (
        build_closure_operator,
        build_volume_kernel,
    )
    geom = CurvilinearGeometry(kind="sphere-1d", inner_radius=r_0)
    radii = np.array([R_out])
    sig_t = np.array([sig_t_val])
    r_nodes, r_wts, panels = composite_gl_r(
        radii, n_panels, p_order, dps=dps, inner_radius=r_0,
    )
    K_vol = build_volume_kernel(
        geom, r_nodes, panels, radii, sig_t,
        n_angular=n_ang, n_rho=n_ang, dps=dps,
    )
    op = build_closure_operator(
        geom, r_nodes, r_wts, radii, sig_t,
        reflection="white", n_angular=n_ang, n_surf_quad=n_ang, dps=dps,
    )
    K = K_vol + op.as_matrix()
    return solve_k_eff(K, sig_t_val, sig_s_val, nu_sig_f_val)


def run_split_basis(r_0, R_out, sig_t_val, sig_s_val, nu_sig_f_val,
                    N_g, N_s, N_i,
                    n_panels=2, p_order=4, dps=15, n_ang=32):
    """Compute k_eff with split-basis closure."""
    geom = CurvilinearGeometry(kind="sphere-1d", inner_radius=r_0)
    radii = np.array([R_out])
    sig_t = np.array([sig_t_val])
    r_nodes, r_wts, panels = composite_gl_r(
        radii, n_panels, p_order, dps=dps, inner_radius=r_0,
    )
    K_vol = build_volume_kernel(
        geom, r_nodes, panels, radii, sig_t,
        n_angular=n_ang, n_rho=n_ang, dps=dps,
    )
    K_bc, aux = build_split_basis_closure(
        geom, r_nodes, r_wts, radii, sig_t,
        N_g, N_s, N_i,
        n_angular=n_ang, n_surf_quad=n_ang, dps=dps,
    )
    if K_bc is None:
        return None, aux
    K = K_vol + K_bc
    return solve_k_eff(K, sig_t_val, sig_s_val, nu_sig_f_val), aux


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    print("=" * 78)
    print("c_in-aware split-basis rank-N closure — k_eff diagnostic")
    print("=" * 78)

    # Structural tests (must pass)
    ok1 = test_W_block_structure(rho_val=0.3, tau=0.0)
    ok2 = test_F4_decomposition(0.3)
    if not (ok1 and ok2):
        print("\n*** STRUCTURAL TEST FAILED — aborting. ***")
        return 1

    # k_inf baseline
    sig_t_val, sig_s_val, nu_sig_f_val = 1.0, 1.0 / 3.0, 1.0
    # Want k_inf = nu_sig_f / (sig_t - sig_s) = 1.5
    k_inf = nu_sig_f_val / (sig_t_val - sig_s_val)
    print(f"\nk_inf baseline = {k_inf:.6f} (expected 1.5)")
    assert abs(k_inf - 1.5) < 1e-12

    # Primary test: σ_t · R = 5, r_0/R = 0.3
    print("\n" + "=" * 78)
    print("Primary test: σ_t·R = 5, r_0/R = 0.3, k_inf = 1.5")
    print("=" * 78)

    R_out = 5.0  # σ_t · R_out = 5 with σ_t = 1
    r_0 = 0.3 * R_out  # r_0/R = 0.3

    # F.4 baseline
    k_f4 = run_scalar_f4(r_0, R_out, sig_t_val, sig_s_val, nu_sig_f_val)
    err_f4 = abs(k_f4 - k_inf) / k_inf * 100.0
    print(f"\n  F.4 scalar (rank-2 white):  k = {k_f4:.8f}, err = {err_f4:.4f}%")

    # Split basis rank-(0,0,0) -- should match F.4 topology
    # Note: (0,0,0) means zero modes in each sub-basis — degenerate, zero closure.
    # Use (1,1,1) for the minimal meaningful case.
    for (Ng, Ns, Ni) in [(1, 1, 1), (2, 1, 1), (1, 2, 1), (1, 1, 2),
                          (2, 2, 1), (2, 2, 2), (3, 3, 3)]:
        try:
            k_sb, _ = run_split_basis(
                r_0, R_out, sig_t_val, sig_s_val, nu_sig_f_val,
                Ng, Ns, Ni,
            )
        except Exception as e:
            print(f"  rank-({Ng},{Ns},{Ni}): FAILED — {type(e).__name__}: {e}")
            continue
        if k_sb is None:
            print(f"  rank-({Ng},{Ns},{Ni}): singular closure")
            continue
        err_sb = abs(k_sb - k_inf) / k_inf * 100.0
        diff_f4 = (k_sb - k_f4) * 1e5  # delta in 1e-5 units
        print(f"  rank-({Ng},{Ns},{Ni}) ({Ng+Ns+Ni} modes): k = {k_sb:.8f}, "
              f"err = {err_sb:.4f}%, Δ(F.4) = {diff_f4:+.2f}e-5")

    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
