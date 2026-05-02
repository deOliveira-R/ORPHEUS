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


def _apply_operator_with_source_profile(
    psi: np.ndarray,
    source_profile: np.ndarray,
    r_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    R: float,
    sigma_t: float,
    alpha: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group Variant α operator.

    Evaluates :math:`\psi_g^{(n+1)}(r,\mu) = \int_0^{L_{\rm back}}
    q_g(|r(s)|)\,e^{-\Sigma_{t,g} s}\,\mathrm d s + \alpha\,
    e^{-\Sigma_{t,g} L_{\rm back}}\,T_\alpha(\mu_{\rm surf})\,B_g`
    along bouncing characteristics, taking the source profile
    :math:`q_g(r_i)` (per steradian per unit volume) as input rather
    than computing it from a scalar coefficient times the in-group
    flux. Multi-group friendly: the caller supplies whatever cross-
    group + fission terms apply.

    Parameters
    ----------
    psi : (n_r, n_mu) ndarray
        Current group's angular flux iterate (used only in vacuum
        branch where ``psi`` doesn't enter — kept for shape/symmetry).
    source_profile : (n_r,) ndarray
        :math:`q_g(r_i) / (4\pi)` — already-divided isotropic source
        per steradian per unit volume. Cubic-spline-interpolated for
        evaluation along trajectory points.
    r_nodes, mu_nodes : ndarray
        Radial GL grid on :math:`(0, R)`, angular GL grid on
        :math:`[-1, 1]`.
    R, sigma_t : float
        Outer radius and per-group total cross section.
    alpha : float
        Surface reflectivity in :math:`[0, 1]`.
    n_traj_quad : int
        Trajectory + bounce-period quadrature order.

    Returns
    -------
    (n_r, n_mu) ndarray
        Updated angular flux for this group.
    """
    source_interp = CubicSpline(r_nodes, source_profile, extrapolate=True)

    # Gauss-Legendre on s ∈ [0, 1] (rescaled per integral).
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

            # Backward-leg geometry (see _apply_operator docstring).
            disc = R * R - r * r * (1.0 - mu * mu)
            sqrt_disc = np.sqrt(max(disc, 0.0))
            L_back = r * mu + sqrt_disc
            mu_surf = sqrt_disc / R
            L_p = 2.0 * R * mu_surf

            # First-leg trajectory integral.
            s_pts = s_unit * L_back
            r_traj_sq = r * r - 2.0 * r * mu * s_pts + s_pts * s_pts
            r_traj = np.sqrt(np.clip(r_traj_sq, 0.0, R * R))
            integrand_F = source_interp(r_traj) * np.exp(-sigma_t * s_pts)
            F = L_back * np.sum(w_unit * integrand_F)

            # Vacuum branch.
            if alpha == 0.0:
                psi_new[i, q_idx] = F
                continue

            atten_first_leg = np.exp(-sigma_t * L_back)

            # Bounce-period integral on antipodal chord.
            s_pts_p = s_unit * L_p
            h_sq = R * R * (1.0 - mu_surf * mu_surf)
            r_chord_sq = h_sq + (s_pts_p - 0.5 * L_p) ** 2
            r_chord = np.sqrt(np.clip(r_chord_sq, 0.0, R * R))
            integrand_B = source_interp(r_chord) * np.exp(-sigma_t * s_pts_p)
            B = L_p * np.sum(w_unit * integrand_B)

            # Geometric bounce-sum closure.
            denom = 1.0 - alpha * np.exp(-sigma_t * L_p)
            psi_surf = alpha * B / denom

            psi_new[i, q_idx] = F + atten_first_leg * psi_surf

    return psi_new


def _apply_operator(
    psi: np.ndarray,
    r_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    mu_weights: np.ndarray,
    R: float,
    sigma_t: float,
    source_coeff: float,
    alpha: float,
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

    Boundary condition is parametrised by :math:`\alpha \in [0, 1]`
    (Sanchez 1986 Eq. A3.a with :math:`\beta = 0`):

    - :math:`\alpha = 1`: perfect specular reflection (closed sphere)
    - :math:`\alpha = 0`: vacuum BC (no surface re-entry)
    - :math:`0 < \alpha < 1`: partial-reflection albedo

    The bounce-sum closure carries :math:`\alpha`:

    .. math::

        \psi_{\rm surf} \;=\; \frac{\alpha\,B(\mu_{\rm surf})}
                                   {1 - \alpha\,e^{-\Sigma_t L_p}}.

    For :math:`\alpha = 0` this collapses to :math:`\psi_{\rm surf} = 0`,
    so the total flux reduces to the first-leg integral
    :math:`\psi(r, \mu) = F(r, \mu)`. The bounce-period integral
    :math:`B` is not computed in this branch (efficiency).

    Parameters
    ----------
    psi : (n_r, n_mu) ndarray
        Current angular flux iterate.
    r_nodes : (n_r,) ndarray
        Radial grid nodes on :math:`(0, R]`.
    mu_nodes : (n_mu,) ndarray
        Angular grid on :math:`[-1, 1]` (cosine of :math:`\Omega` with
        outward radial). Vacuum BC breaks the :math:`\mu \to -\mu`
        symmetry — the full range must be discretised. For closed sphere
        the eigenmode is :math:`\mu`-isotropic anyway, so this incurs
        only redundant compute, not wrong answers.
    mu_weights : (n_mu,) ndarray
        Quadrature weights matching ``mu_nodes`` for the full
        :math:`\int_{-1}^{1} \cdot\,\mathrm d\mu` integral.
    R : float
        Outer radius.
    sigma_t : float
        Total cross section.
    source_coeff : float
        Combined source factor :math:`(\Sigma_s + \nu\Sigma_f/k) / 4\pi`
        — multiplied by :math:`\phi(r)` to give the isotropic source
        per steradian per unit volume.
    alpha : float
        Surface reflectivity in :math:`[0, 1]`. 1 = perfect specular,
        0 = vacuum, intermediate = partial-reflection albedo.
    n_traj_quad : int
        Gauss-Legendre order for the trajectory + bounce-chord
        integrals.

    Returns
    -------
    (n_r, n_mu) ndarray
        Updated angular flux :math:`\psi^{(n+1)}`.
    """
    # Scalar flux φ(r_i) = 2π · ∫_{-1}^{1} ψ(r_i, µ) dµ.
    phi = 2.0 * np.pi * np.sum(psi * mu_weights[None, :], axis=1)
    # 1G source profile: q(r)/(4π) = source_coeff · φ(r). The source_coeff
    # already absorbs the 1/(4π) factor (it's (Σ_s + νΣ_f/k) / (4π)).
    source_profile = source_coeff * phi
    return _apply_operator_with_source_profile(
        psi, source_profile, r_nodes, mu_nodes, R, sigma_t, alpha,
        n_traj_quad=n_traj_quad,
    )


def solve_greens_function_sphere(
    R: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha: float = 1.0,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> GreensFunctionResult:
    r"""Power iteration on the Variant α operator for homogeneous sphere
    with isotropic scattering, parametrised by surface reflectivity α.

    Solves the k-eigenvalue problem :math:`(I - K_{\rm scat})\,\psi =
    (\nu\Sigma_f/k)\,K_{\rm transport}\,\psi` via fission-source
    iteration with the Variant α operator (angle-resolved Green's
    function with bounces summed analytically — Sanchez 1986 Eq. A5
    specular leg).

    Boundary condition is parametrised by :math:`\alpha \in [0, 1]`
    (Sanchez 1986 Eq. A3.a with :math:`\beta = 0`):

    - :math:`\alpha = 1`: perfect specular reflection. Closed sphere,
      no leakage, rank-1 isotropic eigenmode with
      :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a` (V_α1).
      Trivially exact at machine precision (V_α1.numerical).
    - :math:`\alpha = 0`: vacuum BC, no surface re-entry. Spatial
      eigenmode is non-trivial (peaked at center, depleted at surface
      due to leakage). :math:`k_{\rm eff} < k_\infty`. **This case
      stress-tests the trajectory machinery** — closed sphere does not.
      Cross-check against Pomraning-Siewert 1982 closed form.
    - :math:`0 < \alpha < 1`: partial-reflection albedo. Intermediate
      between the two extremes.

    Parameters
    ----------
    R : float
        Outer radius.
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G; multi-group deferred to A3).
    alpha : float, default 1.0
        Surface reflectivity. 1 = closed sphere (specular), 0 = vacuum.
    n_r, n_mu : int
        Radial and angular quadrature orders.
    n_traj_quad : int
        Trajectory + bounce-chord quadrature order.
    max_iter, tol : int, float
        Power-iteration convergence parameters. Vacuum BC needs
        :math:`O(10^2)` iterations from a uniform initial guess; closed
        sphere converges in 1.
    initial_psi : (n_r, n_mu) ndarray, optional
        Initial angular flux. Default: uniform constant.

    Returns
    -------
    :class:`GreensFunctionResult`
        Converged k_eff, ψ, φ, grid info, iteration count.
    """
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha = {alpha} must lie in [0, 1]")

    # Radial grid: Gauss-Legendre on (0, R) — open quadrature avoids
    # the boundary point r = R exactly (where the Phase-5 hypersingular
    # spike lives in the angle-integrated kernel; here the angle-resolved
    # kernel is finite at r = R but we don't need that node).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    r_weights = R * 0.5 * r_quad_wts

    # Angular grid: GL on [-1, 1] (full µ range). For vacuum BC the
    # µ → -µ symmetry is broken (different trajectories for outward vs
    # inward Ω), so both signs are needed. For closed-sphere specular
    # the eigenmode is µ-isotropic anyway, so the extra inward nodes
    # are redundant but harmless.
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    # Initial guess
    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
    else:
        psi = np.ones((n_r, n_mu))

    # Initial k_eff: use k_inf as starting guess where defined.
    # For α < 1 the converged k_eff < k_inf due to leakage; the
    # Rayleigh quotient walks it down during iteration.
    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; non-absorbing "
            "medium not supported (k_inf undefined)."
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
            source_coeff, alpha, n_traj_quad=n_traj_quad,
        )

        # Scalar flux for Rayleigh quotient.
        phi_old = 2.0 * np.pi * np.sum(psi * mu_weights[None, :], axis=1)
        phi_new = 2.0 * np.pi * np.sum(
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

    phi = 2.0 * np.pi * np.sum(psi * mu_weights[None, :], axis=1)

    return GreensFunctionResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group extension (A3) — homogeneous sphere, isotropic scattering
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class GreensFunctionMGResult:
    """Result of multi-group Variant α power iteration."""

    k_eff: float
    psi_g: np.ndarray   # (G, n_r, n_mu) — angular flux per group
    phi_g: np.ndarray   # (G, n_r) — scalar flux per group
    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


def solve_greens_function_sphere_mg(
    R: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), sigma_s[g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,), defaults to (1, 0, ..., 0)
    *,
    alpha: float = 1.0,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> GreensFunctionMGResult:
    r"""Multi-group Variant α power iteration on homogeneous sphere
    with isotropic scattering, parametrised by surface reflectivity α.

    Solves the multi-group k-eigenvalue problem

    .. math::

       \Sigma_{t,g}\,\psi_g(r,\Omega) \;=\;
            \frac{1}{4\pi}\Bigl[
                \sum_{g'}\Sigma_{s,g'\to g}\,\phi_{g'}(r)
                + \frac{\chi_g}{k}\sum_{g'}\nu\Sigma_{f,g'}\,\phi_{g'}(r)
            \Bigr] + \text{(transport operator)}

    via fission-source iteration. At each outer iteration:

    1. Compute scalar fluxes :math:`\phi_g(r)` from the current
       :math:`\psi_g(r,\mu)`.
    2. Compute fission rate :math:`F(r) = \sum_{g'} \nu\Sigma_{f,g'}\,
       \phi_{g'}(r)`.
    3. For each group :math:`g`, build the source profile
       :math:`q_g(r)/(4\pi) = \frac{1}{4\pi}\bigl[\sum_{g'}
       \Sigma_{s,g'\to g}\,\phi_{g'}(r) + (\chi_g/k)\,F(r)\bigr]` and
       apply the per-group operator.
    4. Update :math:`k_{\rm eff}` via Rayleigh quotient on
       volume-integrated fission rate.
    5. Normalise; check convergence.

    Convention: ``sigma_s[g_from, g_to]`` = scatter cross section from
    group ``g_from`` to group ``g_to``. This matches
    :func:`orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
    where the matrix-form transport operator uses ``sigma_s.T``.

    For Option A scope: handles arbitrary G with arbitrary scattering
    matrix (downscatter + upscatter both supported algorithmically;
    upscatter just slows convergence). Cross-check at closed sphere
    (α=1) reduces to ``kinf_homogeneous`` exactly.

    Parameters
    ----------
    R : float
        Outer radius (cm).
    sigma_t : (G,) ndarray
        Per-group total cross section.
    sigma_s : (G, G) ndarray
        Scattering matrix; ``sigma_s[g_from, g_to]``.
    nu_sigma_f : (G,) ndarray
        Per-group production cross section :math:`\nu\Sigma_f`.
    chi : (G,) ndarray, optional
        Fission spectrum, sums to 1. Default: all-fast emission
        (``chi = [1, 0, ..., 0]``).
    alpha : float, default 1.0
        Surface reflectivity. 1 = closed sphere, 0 = vacuum.
    n_r, n_mu, n_traj_quad : int
        Quadrature orders.
    max_iter, tol : int, float
        Power-iteration parameters.
    initial_psi : (G, n_r, n_mu) ndarray, optional
        Initial angular flux. Default: uniform constant per group.
    initial_k : float, optional
        Initial k_eff guess. Default: ``kinf_homogeneous`` from XS.

    Returns
    -------
    :class:`GreensFunctionMGResult`
    """
    sigma_t = np.asarray(sigma_t, dtype=float)
    sigma_s = np.asarray(sigma_s, dtype=float)
    nu_sigma_f = np.asarray(nu_sigma_f, dtype=float)
    G = len(sigma_t)
    if sigma_s.shape != (G, G):
        raise ValueError(
            f"sigma_s shape must be ({G}, {G}); got {sigma_s.shape}"
        )
    if len(nu_sigma_f) != G:
        raise ValueError(
            f"nu_sigma_f length must be {G}; got {len(nu_sigma_f)}"
        )
    if chi is None:
        chi = np.zeros(G)
        chi[0] = 1.0
    else:
        chi = np.asarray(chi, dtype=float)
        if len(chi) != G:
            raise ValueError(
                f"chi length must be {G}; got {len(chi)}"
            )
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha = {alpha} must lie in [0, 1]")

    # Radial grid and angular grid (same as 1G).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    r_weights = R * 0.5 * r_quad_wts

    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    # Initial guess
    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (G, n_r, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({G}, {n_r}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((G, n_r, n_mu))

    # Initial k_eff: use k_inf from infinite-medium balance.
    if initial_k is not None:
        k_eff = float(initial_k)
    else:
        # Diagonal of A = diag(σ_t) - σ_s.T; F = outer(χ, νΣ_f).
        # k_inf = largest eigenvalue of A^{-1} F.
        A = np.diag(sigma_t) - sigma_s.T
        F = np.outer(chi, nu_sigma_f)
        eig = np.linalg.eigvals(np.linalg.solve(A, F))
        k_eff = float(np.real(eig).max())
        if k_eff <= 0:
            raise ValueError(
                f"k_inf estimate non-positive ({k_eff}); XS likely "
                "degenerate."
            )

    # Volume-integrated fission rate for Rayleigh quotient.
    def total_fission_rate(phi_g: np.ndarray) -> float:
        # ∫_0^R [Σ_g νΣ_f,g · φ_g(r)] · 4π r² dr
        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)
        return float(4.0 * np.pi * np.sum(F_r * r_nodes ** 2 * r_weights))

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        # Current scalar fluxes per group: (G, n_r).
        phi_g = 2.0 * np.pi * np.sum(
            psi * mu_weights[None, None, :], axis=2,
        )

        # Fission rate F(r) = Σ_g' νΣ_f,g' · φ_g'(r). Shape (n_r,).
        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)

        # Per-group source profiles q_g(r)/(4π) at r-nodes:
        #   = inv_4pi · [Σ_g' σ_s[g',g] · φ_g'(r) + (χ_g/k) · F(r)]
        # Scatter contribution: einsum 'sg,sr->gr' over s=source group.
        scatter_source = np.einsum('sg,sr->gr', sigma_s, phi_g)
        fission_source = (chi[:, None] / k_eff) * F_r[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        # Apply per-group operator.
        psi_new = np.zeros_like(psi)
        for g in range(G):
            psi_new[g] = _apply_operator_with_source_profile(
                psi[g], source_profile_g[g], r_nodes, mu_nodes, R,
                float(sigma_t[g]), alpha, n_traj_quad=n_traj_quad,
            )

        # Update k_eff via Rayleigh quotient on fission rate.
        phi_g_new = 2.0 * np.pi * np.sum(
            psi_new * mu_weights[None, None, :], axis=2,
        )
        Frate_old = total_fission_rate(phi_g)
        Frate_new = total_fission_rate(phi_g_new)
        if Frate_old < 1e-30:
            raise RuntimeError(
                f"Fission rate vanished at iter {iterations}; XS likely "
                "non-multiplying."
            )
        k_new = k_eff * Frate_new / Frate_old

        # Normalise so the fission rate stays O(1) — prevents drift.
        norm = Frate_new
        psi_normed = psi_new / norm

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    phi_g = 2.0 * np.pi * np.sum(
        psi * mu_weights[None, None, :], axis=2,
    )

    return GreensFunctionMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )
