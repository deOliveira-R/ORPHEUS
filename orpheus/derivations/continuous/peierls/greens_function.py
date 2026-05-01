r"""Plan-2 **Variant α** Green's function reference for sphere homogeneous
specular — research-grade prototype.

Iterates the angle-resolved Green's function :math:`\tilde t(r' \to r,
\mu)` (Sanchez 1986 [SanchezTTSP1986]_ Eq. A5 specular leg, :math:`\alpha
= 1`) along **bouncing characteristics** rather than assembling the
angle-integrated Peierls kernel :math:`g_\alpha(\rho' \to \rho)` (the
hypersingular form that killed Phase 5; see
:file:`/.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`).

Architecture (per :file:`/.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`)
---------------------------------------------------------------------

For each phase-space grid point :math:`(r_i, \mu_q)`:

1. **First-leg trajectory integral** along the chord from :math:`r_i`
   in direction :math:`-\Omega_\mu` to surface entry:

   .. math::

       F(r_i, \mu_q) \;=\; \int_0^{L_0(r_i, \mu_q)}
           q(\bigl|r_i - s\,\Omega_\mu\bigr|)\,e^{-\Sigma_t s}\,\mathrm d s.

2. **Bounce-period integral** along the antipodal chord at impact
   parameter :math:`h(r_i, \mu_q) = R\sqrt{1 - \mu_{\rm surf}^2}`:

   .. math::

       B(\mu_{\rm surf}) \;=\; \int_0^{L_p}
           q(\bigl|r_{\rm chord}(s)\bigr|)\,e^{-\Sigma_t s}\,\mathrm d s,
       \qquad L_p = 2 R \mu_{\rm surf},
       \qquad |r_{\rm chord}(s)| = \sqrt{h^2 + (s - L_p/2)^2}.

3. **Geometric bounce sum** (closed form):

   .. math::

       \psi_{\rm surf}(\mu_{\rm surf}) \;=\;
           T(\mu_{\rm surf})\cdot B(\mu_{\rm surf}),
       \qquad T(\mu) = \frac{1}{1 - e^{-\Sigma_t L_p}}.

4. **Total angular flux**:

   .. math::

       \psi(r_i, \mu_q) \;=\; F(r_i, \mu_q) +
           e^{-\Sigma_t L_0}\,\psi_{\rm surf}(\mu_{\rm surf}).

Each per-pair sample is a 1-D quadrature along a single characteristic
— pointwise finite by V_α1 (operator action on constant trial = ω₀
algebraically; closed form derived via the surface fixed-point
:math:`\psi_{\rm surf} = q/\Sigma_t`). The :math:`\mu \to 0`
geometric-series factor :math:`T(\mu) \sim 1/(2\Sigma_t R \mu)` is
finite per :math:`(r, \mu)` pair; only when computing the scalar flux
:math:`\phi(r) = 2\pi \int \psi(r, \mu)\,\mathrm d\mu` does the
integrable :math:`1/\mu` singularity emerge — handled here by Gauss-
Legendre on :math:`\mu \in (0, 1]` with sufficient nodes (Gauss-Jacobi
or :math:`u^2 = \mu` substitution available as future optimisations).

For the first prototype, the assumptions are:

- Homogeneous sphere (single :math:`\Sigma_t`, :math:`\Sigma_s`,
  :math:`\nu\Sigma_f`); multi-region deferred to Plan 2 follow-on.
- Isotropic scattering (:math:`\omega_1 = 0`).
- Perfect specular BC (:math:`\alpha = 1`).
- :math:`\psi(r, \mu)` symmetric under :math:`\mu \to -\mu`, so we
  discretise :math:`\mu \in (0, 1]` only and double-count for the full
  angular integral. This holds for the rank-1 isotropic eigenmode of
  the closed-sphere problem (V_α1).

Validation status (Plan 2 B4 smoke test)
-----------------------------------------

V_α1 algebraically proves :math:`(K \cdot \mathrm{const})(r, \mu) =
\omega_0 \cdot \mathrm{const}` for closed homogeneous sphere, giving
:math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`. The B4 prototype
verifies the **numerical implementation** matches this algebraic
identity to within quadrature error.

Smoke test target (per Plan 2): thin homogeneous sphere :math:`R = 5`,
:math:`\Sigma_t = 0.5`, fuel-A-like XS, :math:`k_\infty = 0.20833`.
Convergence to within 0.05 % of :math:`k_\infty` at moderate
(:math:`n_r, n_\mu`) is the implementation-correctness gate.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
  DOI: 10.1080/00411458608210456.
- :mod:`orpheus.derivations.continuous.peierls.origins.specular.greens_function`
  (V_α1, V_α2, V_α3 SymPy verifications).
- :mod:`orpheus.derivations.continuous.peierls.origins.specular.continuous_mu`
  (V1-V4 kernel-form identities, predecessor).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline


@dataclass(frozen=True)
class GreensFunctionResult:
    """Result of Variant α power iteration on homogeneous sphere with
    specular BC."""

    k_eff: float
    psi: np.ndarray  # (n_r, n_mu) angular flux on the grid
    phi: np.ndarray  # (n_r,) scalar flux on the radial grid
    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


def _apply_operator(
    psi: np.ndarray,
    r_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    mu_weights: np.ndarray,
    R: float,
    sigma_t: float,
    source_coeff: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Apply the Variant α operator: :math:`\psi^{(n+1)} = K[\psi^{(n)}]`.

    Computes the angle-resolved scattering kernel along bouncing
    characteristics for each grid point :math:`(r_i, \mu_q)`. The
    source ``q(r) = source_coeff * φ(r)`` is built via cubic spline
    interpolation of the scalar flux :math:`\phi(r) = 4\pi \int_0^1
    \psi(r, \mu)\,\mathrm d\mu` (the factor 4π = 2 × 2π combines the
    azimuthal :math:`2\pi` with the :math:`\mu`-symmetry doubling).

    Parameters
    ----------
    psi : (n_r, n_mu) ndarray
        Current angular flux iterate.
    r_nodes : (n_r,) ndarray
        Radial grid nodes on :math:`(0, R]`.
    mu_nodes : (n_mu,) ndarray
        Angular grid on :math:`(0, 1]` (outward direction cosines).
    mu_weights : (n_mu,) ndarray
        Quadrature weights matching ``mu_nodes`` for the full
        :math:`\int_0^1 \cdot\,\mathrm d\mu` integral.
    R : float
        Outer radius.
    sigma_t : float
        Total cross section.
    source_coeff : float
        Combined source factor :math:`(\Sigma_s + \nu\Sigma_f/k) / 4\pi`
        — multiplied by :math:`\phi(r)` to give the isotropic source
        per steradian per unit volume.
    n_traj_quad : int
        Gauss-Legendre order for the trajectory + bounce-chord
        integrals.

    Returns
    -------
    (n_r, n_mu) ndarray
        Updated angular flux :math:`\psi^{(n+1)}`.
    """
    # Scalar flux φ(r_i) = ∫ ψ(r_i, Ω) dΩ = 2 · 2π · ∫_0^1 ψ(r_i, µ) dµ.
    # The factor 2 absorbs the µ → -µ symmetry (we only discretise µ > 0).
    phi = 4.0 * np.pi * np.sum(psi * mu_weights[None, :], axis=1)

    # Cubic spline interpolant for φ(r) — smooth source between grid nodes.
    # Trajectory points satisfy 0 ≤ |r(s)| ≤ R (geometry confines them);
    # GL nodes are open (excludes endpoints), so r_chord can fall slightly
    # below r_nodes[0] (at chord midpoint with small impact parameter)
    # or slightly above r_nodes[-1] (near the surface). Cubic extrapolation
    # is well-behaved for the physically smooth φ(r).
    phi_interp = CubicSpline(r_nodes, phi, extrapolate=True)

    # Gauss-Legendre on s ∈ [0, 1] (we'll rescale per integral).
    s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
    s_unit = 0.5 * (s_quad_raw + 1.0)
    w_unit = 0.5 * w_quad_raw

    n_r = len(r_nodes)
    n_mu = len(mu_nodes)
    psi_new = np.zeros((n_r, n_mu))

    for i in range(n_r):
        r = r_nodes[i]
        for q_idx in range(n_mu):
            mu = mu_nodes[q_idx]

            # Geometry: first-leg distance + surface cosine.
            # Line from r in direction -Ω_µ exits at distance L_0.
            # Discriminant: R² − r²(1 − µ²) — always positive for r < R.
            disc = R * R - r * r * (1.0 - mu * mu)
            sqrt_disc = np.sqrt(max(disc, 0.0))
            L_0 = sqrt_disc - r * mu  # > 0 for any r < R
            mu_surf = sqrt_disc / R   # surface cosine of trajectory
            L_p = 2.0 * R * mu_surf   # antipodal chord = bounce period

            # First-leg trajectory integral.
            #   F = ∫_0^L_0 q(|r(s)|) e^{-σ_t s} ds, q = source_coeff · φ.
            # Trajectory: r(s) = r + s · (-Ω_µ), |r(s)|² = r² − 2rsµ + s².
            s_pts = s_unit * L_0
            r_traj_sq = r * r - 2.0 * r * mu * s_pts + s_pts * s_pts
            r_traj = np.sqrt(np.clip(r_traj_sq, 0.0, R * R))
            phi_traj = phi_interp(r_traj)
            integrand_F = (
                source_coeff * phi_traj * np.exp(-sigma_t * s_pts)
            )
            F = L_0 * np.sum(w_unit * integrand_F)

            # Bounce-period integral on antipodal chord.
            #   B = ∫_0^L_p q(|r_chord(s)|) e^{-σ_t s} ds.
            # Chord parametrization: at impact parameter h, position at
            # arc length s ∈ [0, L_p] is |r_chord(s)|² = h² + (s − L_p/2)².
            s_pts_p = s_unit * L_p
            h_sq = R * R * (1.0 - mu_surf * mu_surf)
            r_chord_sq = h_sq + (s_pts_p - 0.5 * L_p) ** 2
            r_chord = np.sqrt(np.clip(r_chord_sq, 0.0, R * R))
            phi_chord = phi_interp(r_chord)
            integrand_B = (
                source_coeff * phi_chord * np.exp(-sigma_t * s_pts_p)
            )
            B = L_p * np.sum(w_unit * integrand_B)

            # Geometric bounce-sum closure: ψ_surf = T(µ_surf) · B.
            # T(µ) = 1/(1 − e^{-σ_t L_p}) with L_p = 2R·µ_surf.
            denom = 1.0 - np.exp(-sigma_t * L_p)
            T_factor = 1.0 / denom if denom > 1e-15 else 1.0 / (sigma_t * L_p)
            psi_surf = T_factor * B

            # Total angular flux at (r_i, µ_q).
            psi_new[i, q_idx] = F + np.exp(-sigma_t * L_0) * psi_surf

    return psi_new


def solve_greens_function_specular_sphere(
    R: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> GreensFunctionResult:
    r"""Power iteration on the Variant α operator for homogeneous
    sphere with perfect specular BC and isotropic scattering.

    Solves the k-eigenvalue problem :math:`(I - K_{\rm scat})\,\psi =
    (\nu\Sigma_f/k)\,K_{\rm transport}\,\psi` via fission-source
    iteration with the Variant α operator (angle-resolved Green's
    function with specular bounces summed analytically).

    For closed sphere (specular BC, no leakage) the rank-1 isotropic
    mode is the unique eigenmode with :math:`k_{\rm eff} = k_\infty =
    \nu\Sigma_f/\Sigma_a` (V_α1 algebraic identity). This solver
    verifies the **numerical implementation** matches that identity
    within quadrature error — it does not extract physically-different
    information from the specular-BC closed sphere (which is trivially
    :math:`k_\infty`).

    Parameters
    ----------
    R : float
        Outer radius.
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G; multi-group / multi-region deferred).
    n_r, n_mu : int
        Radial and angular quadrature orders.
    n_traj_quad : int
        Trajectory + bounce-chord quadrature order.
    max_iter, tol : int, float
        Power-iteration convergence parameters.
    initial_psi : (n_r, n_mu) ndarray, optional
        Initial angular flux. Default: uniform constant
        (rank-1 isotropic — the algebraic eigenmode).

    Returns
    -------
    :class:`GreensFunctionResult`
        Converged k_eff, ψ, φ, grid info, iteration count.
    """
    # Radial grid: Gauss-Legendre on (0, R) — open quadrature avoids
    # the boundary point r = R exactly (where the Phase-5 hypersingular
    # spike lives in the angle-integrated kernel; here the angle-resolved
    # kernel is finite at r = R but we don't need that node).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    r_weights = R * 0.5 * r_quad_wts

    # Angular grid: GL on (0, 1] (outward µ only; inward µ < 0 absorbed
    # via the µ → -µ symmetry of the homogeneous closed-sphere mode).
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = 0.5 * (mu_quad_pts + 1.0)
    mu_weights = 0.5 * mu_quad_wts

    # Initial guess
    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
    else:
        psi = np.ones((n_r, n_mu))

    # Initial k_eff from k_inf
    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; "
            "k_inf undefined for non-absorbing medium."
        )
    k_eff = nu_sigma_f / sigma_a

    # Volume-integrated fission rate for Rayleigh quotient.
    # Sphere volume integral: ∫_0^R f(r) · 4π r² dr.
    def fission_rate(phi: np.ndarray) -> float:
        return float(
            nu_sigma_f * 4.0 * np.pi
            * np.sum(phi * r_nodes ** 2 * r_weights)
        )

    iterations = 0
    converged = False

    for it in range(max_iter):
        iterations = it + 1
        # Source coefficient: (σ_s + νσ_f/k) / 4π.
        source_coeff = (sigma_s + nu_sigma_f / k_eff) / (4.0 * np.pi)

        psi_new = _apply_operator(
            psi, r_nodes, mu_nodes, mu_weights, R, sigma_t,
            source_coeff, n_traj_quad=n_traj_quad,
        )

        # Scalar flux for Rayleigh quotient.
        phi_old = 4.0 * np.pi * np.sum(psi * mu_weights[None, :], axis=1)
        phi_new = 4.0 * np.pi * np.sum(
            psi_new * mu_weights[None, :], axis=1
        )

        # k_new = k * <νΣ_f, φ_new> / <νΣ_f, φ_old>.
        # For homogeneous medium νΣ_f factors out.
        k_new = k_eff * fission_rate(phi_new) / fission_rate(phi_old)

        # Normalisation: ψ → ψ_new / fission_rate(φ_new) (keeps the
        # iterate from drifting in magnitude).
        norm = fission_rate(phi_new)
        psi_normed = psi_new / norm

        # Convergence on k.
        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    phi = 4.0 * np.pi * np.sum(psi * mu_weights[None, :], axis=1)

    return GreensFunctionResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )
