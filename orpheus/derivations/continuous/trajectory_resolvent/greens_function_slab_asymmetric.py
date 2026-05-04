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
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric`
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

from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    SlabAsymmetricChordOracle,
)
from orpheus.derivations.continuous.trajectory_resolvent.power_iteration import (
    power_iterate_variant_alpha,
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

    Thin facade over :class:`SlabAsymmetricChordOracle.apply_operator`
    (the R3 ChordOracle Protocol). The chord arithmetic + rank-2
    closure live in the oracle; this function preserves the legacy
    call signature for back-compatibility of the public solver API.

    Bit-equal with the pre-R3 inlined body — the oracle was extracted
    verbatim and this facade preserves every FP operation in the same
    order.
    """
    oracle = SlabAsymmetricChordOracle(
        x_nodes=x_nodes, mu_nodes=mu_nodes, L=L,
        alpha_left=alpha_left, alpha_right=alpha_right,
    )
    return oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj_quad,
    )


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

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_x = _scalar_flux_from_psi(psi_iter, mu_weights)
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_iter) * phi_x

        psi_new = _apply_operator_slab_asymmetric(
            source_profile, x_nodes, mu_nodes, L,
            sigma_t, alpha_left, alpha_right,
            n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(psi_new, mu_weights)

        return psi_new, fission_rate(phi_x), fission_rate(phi_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

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

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_g = 2.0 * np.pi * np.einsum(
            'gxm,m->gx', psi_iter, mu_weights,
        )  # (G, n_x)

        F_x = np.einsum('g,gx->x', nu_sigma_f, phi_g)  # (n_x,)
        scatter_source = np.einsum(
            'sg,sx->gx', sigma_s, phi_g,
        )  # (G, n_x)
        fission_source = (chi[:, None] / k_iter) * F_x[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi_iter)
        for g in range(G):
            psi_new[g] = _apply_operator_slab_asymmetric(
                source_profile_g[g], x_nodes, mu_nodes, L,
                float(sigma_t[g]), alpha_left, alpha_right,
                n_traj_quad=n_traj_quad,
            )

        phi_g_new = 2.0 * np.pi * np.einsum(
            'gxm,m->gx', psi_new, mu_weights,
        )
        return psi_new, total_fission_rate(phi_g), total_fission_rate(phi_g_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

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
