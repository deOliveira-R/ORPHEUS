r"""Phase-3A slab **Variant α** Green's function reference (homogeneous
1D slab, symmetric reflective specular BC :math:`\alpha_{\rm left} =
\alpha_{\rm right} = \alpha`) — research-grade prototype.

Mirrors the sphere implementation in :mod:`.greens_function` and the
cylinder implementation in :mod:`.greens_function_cylinder`. Iterates
the angle-resolved Green's function along **bouncing characteristics**
in a 2D phase-space :math:`(x, \mu)` with :math:`x \in [0, L]` and
signed :math:`\mu \in [-1, 1]`.

Structural difference from sphere/cylinder: **2 bounces per period**
---------------------------------------------------------------------

A specular trajectory at fixed :math:`|\mu|` alternates between the two
parallel walls. One full period contains TWO surface reflections —
versus ONE for sphere/cylinder. Consequences:

- Per-period chord: :math:`L_{\rm period} = 2L/|\mu|` (full transit out
  + reverse transit back).
- Per-period reflection product: :math:`\alpha^2` (two reflections per
  period). This enters the geometric resolvent
  :math:`T(\alpha^2, \tau_{\rm period}) = 1/(1 - \alpha^2
  e^{-\Sigma_t L_{\rm period}})`.
- The leading :math:`\alpha` factor in :math:`\psi_{\rm surf} = \alpha
  B T` remains :math:`\alpha^1` — single reflection at the FIRST
  surface arrival.

This is captured back-compatibly by the
:func:`apply_variant_alpha_closure` API extension: slab passes
``alpha_per_period=alpha**2``; sphere/cylinder use the default
(``alpha_per_period = alpha``).

Slab Variant α architecture
----------------------------

For each phase-space grid point :math:`(x_i, \mu_q)`:

1. **First-leg trajectory integral** along the backward chord from
   :math:`x_i` to the first surface arrival. The chord length depends
   on the sign of :math:`\mu`:

   .. math::

      L_{\rm first}(x, \mu) = \begin{cases}
          x / \mu       & \mu > 0 \text{ (came from } x = 0 \text{)} \\
          (L-x) / |\mu| & \mu < 0 \text{ (came from } x = L \text{)}
      \end{cases}.

   The integrand :math:`q(x'(s))\,e^{-\Sigma_t s}` is parametrised by
   the 1D arclength :math:`s \in [0, L_{\rm first}]`; the slab
   position along the chord is :math:`x'(s) = x - s\,\mathrm{sgn}(\mu)`.

2. **Bounce-period integral** along the full out-and-back trajectory at
   conserved :math:`|\mu|`:

   .. math::

      L_{\rm period}(\mu) = \frac{2L}{|\mu|}, \qquad
      x_{\rm chord}(s) = \begin{cases}
          x_{\rm entry} + s\,\mathrm{sgn}(\mu_{\rm out})  & 0 \le s \le L/|\mu| \\
          x_{\rm reflected} - (s - L/|\mu|)\,\mathrm{sgn}(\mu_{\rm out})
                                                            & L/|\mu| \le s \le 2L/|\mu|
      \end{cases}

   where the trajectory enters the slab at the first-arrival wall in
   the original outgoing direction, traverses to the opposite wall
   (full transit at chord :math:`L/|\mu|`), reflects to the antipodal
   :math:`\mu`, and returns (full transit at chord :math:`L/|\mu|`).

3. **Geometric bounce sum** (closed form, 2-bounce-per-period):

   .. math::

      \psi_{\rm surf}(\mu) = \frac{\alpha\,B(\mu)}
                                  {1 - \alpha^2\,e^{-\Sigma_t L_{\rm period}}}.

4. **Total angular flux**: :math:`\psi(x_i, \mu) = F + e^{-\Sigma_t
   L_{\rm first}}\,\psi_{\rm surf}`.

5. **Scalar flux**: :math:`\phi(x) = 2\pi \int_{-1}^{1} \psi(x, \mu)\,
   \mathrm d\mu`.

Symmetry exploited: for the homogeneous closed slab the eigenmode is
:math:`x`-uniform and :math:`\mu`-isotropic. We discretise the full
:math:`\mu` range to support vacuum BC where the symmetry is partially
broken (still :math:`x \to L-x` symmetric, but not :math:`\mu \to -\mu`
when the trajectory machinery treats inward/outward differently).

Grazing-ray pathology
---------------------

At :math:`\mu \to 0` the period collapses to :math:`L_{\rm period} \to
\infty`; :math:`e^{-\Sigma_t L_{\rm period}} \to 0` and the resolvent
:math:`T \to 1`. The bounce-period integral :math:`B \to L_{\rm period}
\cdot (q/\Sigma_t)` for constant source — but the leading
:math:`e^{-\Sigma_t L_{\rm first}}` attenuation in the closure
(:math:`L_{\rm first} \to \infty` as well) drives :math:`\psi_{\rm surf}
\to 0` exponentially. **Slab is structurally immune to the sphere's
:math:`\mu \to 0` Hadamard divergence** — the planar geometry ensures
:math:`L_{\rm first}` and :math:`L_{\rm period}` co-diverge.

We use Gauss-Legendre on :math:`\mu \in [-1, 1]` (open at endpoints)
which avoids the singular :math:`\mu = 0` point exactly.

Assumptions for Phase 3A prototype
-----------------------------------

- Homogeneous slab (single :math:`\Sigma_t`, :math:`\Sigma_s`,
  :math:`\nu\Sigma_f`); multi-region deferred.
- Isotropic scattering.
- **Symmetric** reflective specular BC :math:`\alpha_{\rm left} =
  \alpha_{\rm right} = \alpha \in [0, 1]`. Asymmetric BC (rank-2
  resolvent) is Phase 3B.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- Hébert, A. (2009). *Applied Reactor Physics* §3.8.5 — slab :math:`E_n`
  forms.
- :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_slab`
  — V_α1_slab/V_α2_slab/V_α3_slab SymPy verifications.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3A
  slab Variant α plan.
- Sphere reference solver:
  :mod:`orpheus.derivations.continuous.peierls_greens_function.greens_function`.
- Cylinder reference solver:
  :mod:`orpheus.derivations.continuous.peierls_greens_function.greens_function_cylinder`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from orpheus.derivations.continuous.peierls_greens_function.variant_alpha_core import (
    apply_variant_alpha_closure,
)


@dataclass(frozen=True)
class SlabGreensResult:
    """Result of Variant α power iteration on homogeneous slab with
    symmetric reflective specular BC."""

    k_eff: float
    psi: np.ndarray  # (n_x, n_mu) angular flux on the grid
    phi: np.ndarray  # (n_x,) scalar flux on the spatial grid
    x_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


@dataclass(frozen=True)
class SlabGreensMGResult:
    """Result of multi-group Variant α power iteration on slab."""

    k_eff: float
    psi_g: np.ndarray   # (G, n_x, n_mu)
    phi_g: np.ndarray   # (G, n_x)
    x_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


# ═══════════════════════════════════════════════════════════════════════
# Geometry helpers — slab chord arithmetic
# ═══════════════════════════════════════════════════════════════════════


def _first_leg_chord_slab(x: float, mu: float, L: float) -> float:
    r"""Backward chord length from interior point :math:`x` (with
    direction-cosine :math:`\mu`) to the first surface arrival on the
    slab boundary.

    .. math::

        L_{\rm first}(x, \mu) = \begin{cases}
            x / \mu       & \mu > 0 \\
            (L-x) / |\mu| & \mu < 0
        \end{cases}.

    At :math:`\mu = 0` the chord is undefined (returns ``np.inf``);
    callers must avoid this point via open quadrature on :math:`(-1, 1)`.

    Parameters
    ----------
    x : float
        Interior position in :math:`[0, L]`.
    mu : float
        Signed direction-cosine in :math:`[-1, 1]`.
    L : float
        Slab width.
    """
    if mu > 0:
        return x / mu
    if mu < 0:
        return (L - x) / abs(mu)
    return np.inf


def _bounce_period_chord_slab(mu: float, L: float) -> float:
    r"""Full bounce-period chord (out + back transit) at :math:`|\mu|`.

    :math:`L_{\rm period}(\mu) = 2 L / |\mu|`. At :math:`\mu = 0`
    returns ``np.inf`` (degenerate grazing).
    """
    if mu == 0.0:
        return np.inf
    return 2.0 * L / abs(mu)


# ═══════════════════════════════════════════════════════════════════════
# Per-group Variant α operator (slab)
# ═══════════════════════════════════════════════════════════════════════


def _apply_operator_slab(
    source_profile: np.ndarray,
    x_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    L: float,
    sigma_t: float,
    alpha: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group slab Variant α operator.

    Evaluates :math:`\psi^{(n+1)}_g(x, \mu) = F + e^{-\Sigma_t L_{\rm
    first}}\,\alpha B / (1 - \alpha^2 e^{-\Sigma_t L_{\rm period}})`
    for every phase-space grid point :math:`(x_i, \mu_q)` along the
    bouncing characteristic.

    The 2-bounce-per-period structure enters via
    ``alpha_per_period = alpha**2`` in the call to
    :func:`apply_variant_alpha_closure`.

    Parameters
    ----------
    source_profile : (n_x,) ndarray
        :math:`q_g(x_i)/(4\pi)` — already-divided isotropic source
        per steradian per unit volume.
    x_nodes : (n_x,) ndarray
        Spatial Gauss-Legendre nodes on :math:`(0, L)`.
    mu_nodes : (n_mu,) ndarray
        Direction-cosine Gauss-Legendre nodes on :math:`(-1, 1)`
        (open quadrature avoids :math:`\mu = 0` exactly).
    L, sigma_t : float
        Slab width and per-group total cross section.
    alpha : float
        Surface reflectivity (symmetric, both walls). In :math:`[0, 1]`.
    n_traj_quad : int
        Trajectory + bounce-period quadrature order.

    Returns
    -------
    (n_x, n_mu) ndarray
        Updated angular flux.
    """
    source_interp = CubicSpline(x_nodes, source_profile, extrapolate=True)

    s_quad_raw, w_quad_raw = np.polynomial.legendre.leggauss(n_traj_quad)
    s_unit = 0.5 * (s_quad_raw + 1.0)
    w_unit = 0.5 * w_quad_raw

    n_x = len(x_nodes)
    n_mu = len(mu_nodes)
    psi_new = np.zeros((n_x, n_mu))

    for i in range(n_x):
        x = x_nodes[i]
        for q_idx in range(n_mu):
            mu = mu_nodes[q_idx]
            abs_mu = abs(mu)

            # First-leg backward chord and direction.
            # mu > 0: came from x = 0, backward direction is -x_hat,
            #         so x'(s) = x - s; entry wall at x = 0.
            # mu < 0: came from x = L, backward direction is +x_hat,
            #         so x'(s) = x + s; entry wall at x = L.
            if mu > 0:
                L_first = x / mu
                # Backward chord parametrisation: x' = x - s, s ∈ [0, L_first].
                # Entry point at s = L_first is x = 0 (matches x'(L_first) = 0).
                s_pts_first = s_unit * L_first
                x_traj = x - s_pts_first
                # Bounce-period trajectory: starts at entry wall (x=0)
                # in direction +mu_hat (mu > 0); traverses full slab
                # to x = L (s_2 = L/mu = L/|mu|), reflects to -mu,
                # traverses back to x = 0 (s_2 = 2L/|mu|).
                # Full period s ∈ [0, 2L/|mu|]:
                # x_chord(s) = s · |mu| / 1 if s ≤ L/|mu| (going +x)
                #            = 2L - s · |mu| if s ≥ L/|mu| (going -x).
                # More compactly: x_chord(s) = L - |L - s · |mu||.
                L_half = L / abs_mu  # = L/|mu|, half-period chord arclength
                L_period = 2.0 * L_half
                s_pts_period = s_unit * L_period
                # x_chord ∈ [0, L]; folded distance from x=0:
                x_chord = L - np.abs(L - s_pts_period * abs_mu)
            else:
                # mu < 0
                L_first = (L - x) / abs_mu
                s_pts_first = s_unit * L_first
                # Backward direction is +x_hat: x' = x + s.
                # Entry point at s = L_first is x = L (matches x'(L_first) = L).
                x_traj = x + s_pts_first
                # Bounce period: enters at x = L in direction -mu_hat
                # (mu < 0); traverses to x = 0 (s_2 = L/|mu|), reflects,
                # traverses back to x = L.
                L_half = L / abs_mu
                L_period = 2.0 * L_half
                s_pts_period = s_unit * L_period
                # Mirror image: x_chord starts at L, decreases linearly,
                # bounces at 0, climbs back. Folded distance from x=L:
                x_chord = np.abs(L - s_pts_period * abs_mu)

            x_traj = np.clip(x_traj, 0.0, L)
            x_chord = np.clip(x_chord, 0.0, L)

            # First-leg integral: ∫_0^{L_first} q(x'(s)) e^{-Σ_t s} ds.
            integrand_F = (
                source_interp(x_traj) * np.exp(-sigma_t * s_pts_first)
            )
            F = L_first * np.sum(w_unit * integrand_F)

            # Vacuum branch — short-circuit before computing B.
            if alpha == 0.0:
                psi_new[i, q_idx] = F
                continue

            # Bounce-period integral: ∫_0^{L_period} q(x_chord(s)) e^{-Σ_t s} ds.
            integrand_B = (
                source_interp(x_chord) * np.exp(-sigma_t * s_pts_period)
            )
            B = L_period * np.sum(w_unit * integrand_B)

            # Shared Variant α closure with slab's 2-bounce-per-period
            # reflection product (alpha_per_period = alpha**2).
            psi_new[i, q_idx] = apply_variant_alpha_closure(
                F=F, B=B,
                tau_first_leg=sigma_t * L_first,
                tau_period=sigma_t * L_period,
                alpha=alpha,
                alpha_per_period=alpha * alpha,
            )

    return psi_new


# ═══════════════════════════════════════════════════════════════════════
# Scalar-flux extractor
# ═══════════════════════════════════════════════════════════════════════


def _scalar_flux_from_psi(
    psi: np.ndarray,
    mu_weights: np.ndarray,
) -> np.ndarray:
    r"""Reduce angular flux to scalar flux on the spatial grid.

    .. math::

        \phi(x) = 2\pi \int_{-1}^{1} \psi(x, \mu)\,\mathrm d\mu.

    The factor :math:`2\pi` converts from the 1D angular distribution
    in :math:`\mu` to the full :math:`\Omega`-integrated scalar flux
    (azimuthal :math:`2\pi`).

    Parameters
    ----------
    psi : (..., n_x, n_mu) ndarray
    mu_weights : (n_mu,) ndarray

    Returns
    -------
    (..., n_x) ndarray
    """
    return 2.0 * np.pi * np.einsum('...xm,m->...x', psi, mu_weights)


# ═══════════════════════════════════════════════════════════════════════
# 1G homogeneous solver
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_slab(
    L: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha: float = 1.0,
    n_x: int = 16,
    n_mu: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> SlabGreensResult:
    r"""Power iteration on the slab Variant α operator (homogeneous,
    isotropic scattering, symmetric reflective specular BC).

    Solves the k-eigenvalue problem on a homogeneous 1D slab of width
    :math:`L` via fission-source iteration with the slab Variant α
    operator (angle-resolved Green's function with bouncing
    characteristics summed analytically — 2 bounces per period).

    Boundary condition is parametrised by :math:`\alpha \in [0, 1]`
    applied SYMMETRICALLY to both walls:

    - :math:`\alpha = 1`: closed slab (perfect specular both walls,
      no leakage). :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`
      EXACTLY by V_α1_slab. Load-bearing Phase-3A acceptance test.
    - :math:`\alpha = 0`: vacuum on both walls. Spatial eigenmode
      peaked at center, depleted at walls. :math:`k_{\rm eff} <
      k_\infty`.
    - :math:`0 < \alpha < 1`: partial-reflection albedo on both walls.

    Asymmetric BC (:math:`\alpha_{\rm left} \ne \alpha_{\rm right}`)
    requires the rank-2 resolvent and is deferred to Phase 3B.

    Parameters
    ----------
    L : float
        Slab width (cm).
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G).
    alpha : float, default 1.0
        Symmetric wall reflectivity.
    n_x : int
        Spatial Gauss-Legendre quadrature order on :math:`(0, L)`.
    n_mu : int
        Direction-cosine Gauss-Legendre order on :math:`(-1, 1)`.
    n_traj_quad : int
        Trajectory + bounce-period quadrature order.
    max_iter, tol : int, float
        Power-iteration parameters.
    initial_psi : (n_x, n_mu) ndarray, optional

    Returns
    -------
    :class:`SlabGreensResult`
    """
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha = {alpha} must lie in [0, 1]")

    # Spatial grid: GL on (0, L).
    x_quad_pts, x_quad_wts = np.polynomial.legendre.leggauss(n_x)
    x_nodes = L * 0.5 * (x_quad_pts + 1.0)
    x_weights = L * 0.5 * x_quad_wts

    # Direction-cosine grid: GL on [-1, 1] (open avoids µ = 0 exactly).
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (n_x, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({n_x}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((n_x, n_mu))

    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; non-absorbing "
            "medium not supported."
        )
    k_eff = nu_sigma_f / sigma_a

    # Volume-integrated fission rate per unit transverse area:
    # ∫_0^L νΣ_f φ(x) dx.
    def fission_rate(phi_x: np.ndarray) -> float:
        return float(nu_sigma_f * np.sum(phi_x * x_weights))

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_x = _scalar_flux_from_psi(psi, mu_weights)
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_eff) * phi_x

        psi_new = _apply_operator_slab(
            source_profile, x_nodes, mu_nodes, L,
            sigma_t, alpha, n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(psi_new, mu_weights)

        Frate_old = fission_rate(phi_x)
        Frate_new = fission_rate(phi_new)
        if Frate_old < 1e-30:
            raise RuntimeError(
                f"Fission rate vanished at iter {iterations}; "
                "non-multiplying medium."
            )
        k_new = k_eff * Frate_new / Frate_old
        psi_normed = psi_new / Frate_new

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    phi_x = _scalar_flux_from_psi(psi, mu_weights)

    return SlabGreensResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi_x,
        x_nodes=x_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_slab_mg(
    L: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), [g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,)
    *,
    alpha: float = 1.0,
    n_x: int = 16,
    n_mu: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> SlabGreensMGResult:
    r"""Multi-group slab Variant α power iteration (homogeneous,
    isotropic scattering, symmetric reflective specular BC).

    Multi-group analog of :func:`solve_greens_function_slab`.
    Convention: ``sigma_s[g_from, g_to]``.

    At :math:`\alpha = 1` (closed slab) reduces exactly to
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` —
    the load-bearing MG verification.

    Parameters
    ----------
    L : float
        Slab width (cm).
    sigma_t : (G,) ndarray
    sigma_s : (G, G) ndarray
    nu_sigma_f : (G,) ndarray
    chi : (G,) ndarray, optional
        Fission spectrum. Default: all-fast emission.
    alpha : float, default 1.0
    n_x, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (G, n_x, n_mu) ndarray, optional
    initial_k : float, optional

    Returns
    -------
    :class:`SlabGreensMGResult`
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
            raise ValueError(f"chi length must be {G}; got {len(chi)}")
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha = {alpha} must lie in [0, 1]")

    x_quad_pts, x_quad_wts = np.polynomial.legendre.leggauss(n_x)
    x_nodes = L * 0.5 * (x_quad_pts + 1.0)
    x_weights = L * 0.5 * x_quad_wts

    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (G, n_x, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({G}, {n_x}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((G, n_x, n_mu))

    if initial_k is not None:
        k_eff = float(initial_k)
    else:
        # Use kinf as initial guess (largest eigenvalue of A^{-1} F).
        A = np.diag(sigma_t) - sigma_s.T
        F = np.outer(chi, nu_sigma_f)
        eig = np.linalg.eigvals(np.linalg.solve(A, F))
        k_eff = float(np.real(eig).max())
        if k_eff <= 0:
            raise ValueError(f"k_inf estimate non-positive ({k_eff}).")

    def total_fission_rate(phi_g_x: np.ndarray) -> float:
        F_x = np.einsum('g,gx->x', nu_sigma_f, phi_g_x)
        return float(np.sum(F_x * x_weights))

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_g = 2.0 * np.pi * np.einsum(
            'gxm,m->gx', psi, mu_weights,
        )  # (G, n_x)

        F_x = np.einsum('g,gx->x', nu_sigma_f, phi_g)  # (n_x,)
        scatter_source = np.einsum(
            'sg,sx->gx', sigma_s, phi_g,
        )  # (G, n_x)
        fission_source = (chi[:, None] / k_eff) * F_x[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi)
        for g in range(G):
            psi_new[g] = _apply_operator_slab(
                source_profile_g[g], x_nodes, mu_nodes, L,
                float(sigma_t[g]), alpha, n_traj_quad=n_traj_quad,
            )

        phi_g_new = 2.0 * np.pi * np.einsum(
            'gxm,m->gx', psi_new, mu_weights,
        )
        Frate_old = total_fission_rate(phi_g)
        Frate_new = total_fission_rate(phi_g_new)
        if Frate_old < 1e-30:
            raise RuntimeError(
                f"Fission rate vanished at iter {iterations}."
            )
        k_new = k_eff * Frate_new / Frate_old
        psi_normed = psi_new / Frate_new

        rel_dk = abs(k_new - k_eff) / max(abs(k_eff), 1e-30)
        psi = psi_normed
        k_eff = k_new

        if rel_dk < tol:
            converged = True
            break

    phi_g = 2.0 * np.pi * np.einsum('gxm,m->gx', psi, mu_weights)

    return SlabGreensMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        x_nodes=x_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )
