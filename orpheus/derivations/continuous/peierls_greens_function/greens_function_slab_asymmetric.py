r"""Phase-3B asymmetric slab **Variant α** Green's function reference
(homogeneous 1D slab, **independent** specular reflectivities
:math:`\alpha_L \in [0, 1]` and :math:`\alpha_R \in [0, 1]`) — research-
grade prototype.

Standalone module, parallel to :mod:`.greens_function_slab` (Phase 3A
symmetric :math:`\alpha_L = \alpha_R`). Implements the **rank-2
boundary-to-boundary scattering resolvent** :math:`T = (I - S)^{-1}`
for the asymmetric-slab monodromy

.. math::

    S(\alpha_L, \alpha_R, \tau) = \begin{pmatrix}
        0                              & \alpha_L\,e^{-\tau} \\
        \alpha_R\,e^{-\tau}            & 0
    \end{pmatrix}, \qquad
    \tau = \Sigma_t\,L /|\mu|

(single-transit optical depth — NOT the full out+back period).

Why a separate module
---------------------

The Phase-3A symmetric slab module (:mod:`.greens_function_slab`) uses
the rank-1 scalar resolvent :math:`T(\alpha^2, 2\tau) = 1/(1 - \alpha^2
e^{-2\tau})` with full bounce-period (out + back) source integrals.
The Phase-3B asymmetric slab uses the rank-2 matrix resolvent with
**single-transit** source integrals — a structurally different
arithmetic form that is mathematically equivalent at :math:`\alpha_L =
\alpha_R = \alpha` but NOT bit-equal at the IEEE-754 level. Keeping
them as parallel modules:

- preserves the Phase-3A bit-equal invariant for sphere/cylinder/
  symmetric-slab regression (rank-1 callers do not pay the matrix
  inversion cost),
- isolates the rank-2 algebra at one call site, and
- gives the cross-domain frame's predicted "rank-1 vs rank-2 split"
  a clean implementation boundary.

The reduce-to-symmetric consistency check (rank-2 with :math:`\alpha_L
= \alpha_R = \alpha` agrees with Phase-3A rank-1 to ≤ 1e-10) is the
solver-level L1 gate verifying the two arithmetic forms match.

Asymmetric slab Variant α architecture
---------------------------------------

Phase-space: :math:`(x, \mu)` with :math:`x \in [0, L]`, signed
:math:`\mu \in [-1, 1]`. For each interior point :math:`(x_i, \mu_q)`:

1. **First-leg trajectory integral** (sign-dependent backward chord):

   .. math::

      L_{\rm first}(x, \mu) = \begin{cases}
          x / \mu       & \mu > 0 \text{ (came from } x = 0 \text{)} \\
          (L-x) / |\mu| & \mu < 0 \text{ (came from } x = L \text{)}
      \end{cases}.

2. **Single-transit B integrals**: source contributions along the
   single chord from each wall in the trajectory direction at conserved
   :math:`|\mu|`. Both are evaluated at every grid point because both
   feed the rank-2 closure regardless of which wall is the first
   surface arrival:

   .. math::

      B_{LR}(\mu) &= \int_0^{L/|\mu|}
                       q(0 + s\cdot|\mu|)\,e^{-\Sigma_t s}\,
                       \mathrm d s \\
      B_{RL}(\mu) &= \int_0^{L/|\mu|}
                       q(L - s\cdot|\mu|)\,e^{-\Sigma_t s}\,
                       \mathrm d s.

3. **Rank-2 surface closure**: for :math:`\mu > 0` the first surface
   arrival is the LEFT wall, so the interior reconstruction uses
   :math:`\psi_L^+(\mu)`. For :math:`\mu < 0` the first arrival is the
   RIGHT wall, so :math:`\psi_R^-(\mu)` is used. The rank-2 matrix
   inversion is encoded in :func:`apply_variant_alpha_closure_rank2`
   (shared in :mod:`.variant_alpha_core`).

4. **Total angular flux**: :math:`\psi(x_i, \mu) = F + e^{-\Sigma_t
   L_{\rm first}}\,\psi_{\rm surface}`, where :math:`\psi_{\rm
   surface}` is the appropriate component of :math:`\psi_{\rm surf}`.

5. **Scalar flux**: :math:`\phi(x) = 2\pi \int_{-1}^{1} \psi(x, \mu)\,
   \mathrm d\mu`.

Method-of-images symmetry test (load-bearing)
----------------------------------------------

The user-pinned acceptance test for Phase 3B: the slab :math:`[0, L=1]`
with :math:`\alpha_L = 1, \alpha_R = 0` (left reflective, right vacuum)
**MUST give the same** :math:`k_{\rm eff}` **and the same** :math:`\phi(x
\in [0, 1])` **as a symmetric vacuum slab** :math:`[-1, 1]` **with**
:math:`\alpha_L = \alpha_R = 0` **on both ends.** Reasoning: the
reflective BC at :math:`x = 0` enforces the even-symmetry plane
satisfied by the fundamental eigenmode of the symmetric vacuum slab —
by method of images, the half-domain solution is identical to the
right half of the doubled-domain solution.

This test is the canonical structural-independence check for the
rank-2 prototype: success implies the rank-2 algebra correctly
reproduces the method-of-images equivalence built into the underlying
transport equation.

Assumptions for Phase 3B prototype
-----------------------------------

- Homogeneous slab (single :math:`\Sigma_t`, :math:`\Sigma_s`,
  :math:`\nu\Sigma_f`); multi-region deferred.
- Isotropic scattering.
- **Independent** reflective specular BCs :math:`\alpha_L \in [0, 1]`
  and :math:`\alpha_R \in [0, 1]` — the rank-2 resolvent handles any
  combination including symmetric corner cases.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- Hébert, A. (2009). *Applied Reactor Physics* §3.8.5 — slab :math:`E_n`
  forms and rank-1 white-BC closure.
- :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_slab_asymmetric`
  — V_α1_slab_asym/V_α2_slab_asym/V_α3_slab_asym SymPy verifications.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3B
  asymmetric slab plan.
- :mod:`.greens_function_slab` — Phase-3A symmetric-slab reference
  solver (the structural template this module mirrors with the rank-2
  resolvent).
- :file:`/.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match memo predicting the rank-2 generalisation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline

from orpheus.derivations.continuous.peierls_greens_function.variant_alpha_core import (
    apply_variant_alpha_closure_rank2,
)


@dataclass(frozen=True)
class SlabAsymmetricGreensResult:
    """Result of Variant α power iteration on homogeneous slab with
    INDEPENDENT specular reflectivities :math:`\\alpha_L, \\alpha_R`."""

    k_eff: float
    psi: np.ndarray  # (n_x, n_mu) angular flux on the grid
    phi: np.ndarray  # (n_x,) scalar flux on the spatial grid
    x_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


@dataclass(frozen=True)
class SlabAsymmetricGreensMGResult:
    """Result of multi-group Variant α power iteration on asymmetric
    slab."""

    k_eff: float
    psi_g: np.ndarray   # (G, n_x, n_mu)
    phi_g: np.ndarray   # (G, n_x)
    x_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


# ═══════════════════════════════════════════════════════════════════════
# Per-group Variant α operator (asymmetric slab, rank-2 resolvent)
# ═══════════════════════════════════════════════════════════════════════


def _apply_operator_slab_asymmetric(
    source_profile: np.ndarray,
    x_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    L: float,
    sigma_t: float,
    alpha_left: float,
    alpha_right: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group asymmetric-slab Variant α operator with rank-2 resolvent.

    For each :math:`(x_i, \mu_q)` evaluates:

    1. The first-leg :math:`F` integral along the backward chord to
       the first surface arrival.
    2. The two single-transit :math:`B` integrals: :math:`B_{LR}` from
       :math:`x = 0` to :math:`x = L` and :math:`B_{RL}` from :math:`x =
       L` to :math:`x = 0`, both at conserved :math:`|\mu|`.
    3. The rank-2 surface closure
       :func:`apply_variant_alpha_closure_rank2` selecting the
       appropriate surface flux component (:math:`\psi_L^+` for
       :math:`\mu > 0`, :math:`\psi_R^-` for :math:`\mu < 0`).
    4. The interior reconstruction :math:`F + e^{-\tau_{\rm first\,leg}}
       \,\psi_{\rm surface}`.

    Parameters
    ----------
    source_profile : (n_x,) ndarray
        :math:`q_g(x_i)/(4\pi)` — already-divided isotropic source.
    x_nodes : (n_x,) ndarray
    mu_nodes : (n_mu,) ndarray
    L, sigma_t : float
    alpha_left, alpha_right : float
        Per-wall reflectivities in :math:`[0, 1]`.
    n_traj_quad : int
        Trajectory + single-transit quadrature order.

    Returns
    -------
    (n_x, n_mu) ndarray
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

            # First-leg backward chord. The trajectory direction is +μ,
            # so the backward chord goes in -μ direction. Per the
            # transport equation's formal solution along characteristics,
            # at arclength s back from (x, μ), the position is
            # x'(s) = x - μ·s (with the appropriate sign of μ for the
            # 1D slab — μ > 0 case has x' = x - μ·s; μ < 0 case has
            # x' = x - μ·s = x + |μ|·s).
            if mu > 0:
                L_first = x / mu
                s_pts_first = s_unit * L_first
                x_traj = x - mu * s_pts_first
                surface = "left"
            else:
                # mu < 0
                L_first = (L - x) / abs_mu
                s_pts_first = s_unit * L_first
                x_traj = x + abs_mu * s_pts_first
                surface = "right"

            x_traj = np.clip(x_traj, 0.0, L)

            # Single-transit chord arclength.
            L_transit = L / abs_mu
            s_pts_transit = s_unit * L_transit

            # B_LR: chord from x = 0 to x = L in the +x direction at
            # conserved |mu|. Position along chord:
            #   x_chord(s) = 0 + |mu|·s, s ∈ [0, L_transit].
            x_chord_LR = np.clip(abs_mu * s_pts_transit, 0.0, L)

            # B_RL: chord from x = L to x = 0 in the -x direction.
            #   x_chord(s) = L - |mu|·s, s ∈ [0, L_transit].
            x_chord_RL = np.clip(L - abs_mu * s_pts_transit, 0.0, L)

            # First-leg integral: F = ∫_0^{L_first} q(x'(s)) e^{-Σ_t s} ds.
            integrand_F = (
                source_interp(x_traj) * np.exp(-sigma_t * s_pts_first)
            )
            F = L_first * np.sum(w_unit * integrand_F)

            # Vacuum-vacuum branch — both walls absorbing, ψ_surf = 0.
            if alpha_left == 0.0 and alpha_right == 0.0:
                psi_new[i, q_idx] = F
                continue

            # Single-transit B integrals.
            integrand_B_LR = (
                source_interp(x_chord_LR) * np.exp(-sigma_t * s_pts_transit)
            )
            B_LR = L_transit * np.sum(w_unit * integrand_B_LR)

            integrand_B_RL = (
                source_interp(x_chord_RL) * np.exp(-sigma_t * s_pts_transit)
            )
            B_RL = L_transit * np.sum(w_unit * integrand_B_RL)

            psi_new[i, q_idx] = apply_variant_alpha_closure_rank2(
                F=F, B_RL=B_RL, B_LR=B_LR,
                tau_first_leg=sigma_t * L_first,
                tau_single_transit=sigma_t * L_transit,
                alpha_left=alpha_left,
                alpha_right=alpha_right,
                surface=surface,
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
    """
    return 2.0 * np.pi * np.einsum('...xm,m->...x', psi, mu_weights)


# ═══════════════════════════════════════════════════════════════════════
# 1G homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_slab_asymmetric(
    L: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha_left: float = 1.0,
    alpha_right: float = 1.0,
    n_x: int = 16,
    n_mu: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> SlabAsymmetricGreensResult:
    r"""Power iteration on the asymmetric-slab Variant α operator
    (homogeneous, isotropic scattering, rank-2 resolvent).

    Solves the k-eigenvalue problem on a homogeneous 1D slab of width
    :math:`L` with INDEPENDENT specular reflectivities on the two
    walls. Uses the rank-2 boundary-to-boundary scattering resolvent
    :math:`T = (I - S)^{-1}`.

    Parameters
    ----------
    L : float
        Slab width (cm).
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G).
    alpha_left : float, default 1.0
        Left-wall (:math:`x = 0`) specular reflectivity in :math:`[0, 1]`.
    alpha_right : float, default 1.0
        Right-wall (:math:`x = L`) specular reflectivity in :math:`[0, 1]`.
    n_x, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (n_x, n_mu) ndarray, optional

    Returns
    -------
    :class:`SlabAsymmetricGreensResult`
    """
    if not (0.0 <= alpha_left <= 1.0):
        raise ValueError(f"alpha_left = {alpha_left} must lie in [0, 1]")
    if not (0.0 <= alpha_right <= 1.0):
        raise ValueError(f"alpha_right = {alpha_right} must lie in [0, 1]")

    # Spatial grid: GL on (0, L).
    x_quad_pts, x_quad_wts = np.polynomial.legendre.leggauss(n_x)
    x_nodes = L * 0.5 * (x_quad_pts + 1.0)
    x_weights = L * 0.5 * x_quad_wts

    # Direction-cosine grid: GL on [-1, 1].
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

    def fission_rate(phi_x: np.ndarray) -> float:
        return float(nu_sigma_f * np.sum(phi_x * x_weights))

    iterations = 0
    converged = False
    inv_4pi = 1.0 / (4.0 * np.pi)

    for it in range(max_iter):
        iterations = it + 1

        phi_x = _scalar_flux_from_psi(psi, mu_weights)
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_eff) * phi_x

        psi_new = _apply_operator_slab_asymmetric(
            source_profile, x_nodes, mu_nodes, L,
            sigma_t, alpha_left, alpha_right,
            n_traj_quad=n_traj_quad,
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

    return SlabAsymmetricGreensResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi_x,
        x_nodes=x_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_slab_asymmetric_mg(
    L: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), [g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,)
    *,
    alpha_left: float = 1.0,
    alpha_right: float = 1.0,
    n_x: int = 16,
    n_mu: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> SlabAsymmetricGreensMGResult:
    r"""Multi-group asymmetric-slab Variant α power iteration
    (homogeneous, isotropic scattering, rank-2 resolvent).

    Multi-group analog of :func:`solve_greens_function_slab_asymmetric`.
    Convention: ``sigma_s[g_from, g_to]``. At :math:`\alpha_L =
    \alpha_R = 1` reduces exactly to
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`.

    Parameters
    ----------
    L : float
    sigma_t : (G,) ndarray
    sigma_s : (G, G) ndarray
    nu_sigma_f : (G,) ndarray
    chi : (G,) ndarray, optional
        Fission spectrum. Default: all-fast emission.
    alpha_left, alpha_right : float
    n_x, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (G, n_x, n_mu) ndarray, optional
    initial_k : float, optional

    Returns
    -------
    :class:`SlabAsymmetricGreensMGResult`
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
    if not (0.0 <= alpha_left <= 1.0):
        raise ValueError(f"alpha_left = {alpha_left} must lie in [0, 1]")
    if not (0.0 <= alpha_right <= 1.0):
        raise ValueError(f"alpha_right = {alpha_right} must lie in [0, 1]")

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
            psi_new[g] = _apply_operator_slab_asymmetric(
                source_profile_g[g], x_nodes, mu_nodes, L,
                float(sigma_t[g]), alpha_left, alpha_right,
                n_traj_quad=n_traj_quad,
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

    return SlabAsymmetricGreensMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        x_nodes=x_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )
