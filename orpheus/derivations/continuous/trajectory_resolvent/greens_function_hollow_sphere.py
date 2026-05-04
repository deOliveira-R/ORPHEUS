r"""Phase-3C-1 hollow sphere **Variant α** Green's function reference
(homogeneous shell :math:`r \in [R_{\rm in}, R_{\rm out}]`, **independent**
specular reflectivities :math:`\alpha_{\rm in} \in [0, 1]` on the inner
cavity surface and :math:`\alpha_{\rm out} \in [0, 1]` on the outer
surface) — research-grade prototype.

Standalone module, parallel to :mod:`.greens_function` (solid sphere)
and :mod:`.greens_function_slab_asymmetric` (asymmetric slab). Implements
the **rank-2 boundary-to-boundary scattering resolvent** on the through-
ray phase-space subset (impact parameter :math:`b \le R_{\rm in}`) and
the **rank-1 outer-only resolvent** on the outer-only subset
(:math:`b > R_{\rm in}`).

Phase-space partition by impact parameter
------------------------------------------

For each interior point :math:`(r, \mu)` with :math:`r \in [R_{\rm in},
R_{\rm out}]` and :math:`\mu \in [-1, 1]`, the **impact parameter**
:math:`b = r\sqrt{1 - \mu^2}` is conserved along straight-line travel.
Two cases:

- :math:`b > R_{\rm in}` — outer-only rays. The trajectory line does not
  intersect the inner sphere; the particle bounces between two points
  on the OUTER surface only. **Closure rank-1** (orbit-space class:
  one-surface-compact), structurally identical to a solid-sphere
  ray at the same outer radius.
- :math:`b \le R_{\rm in}` — through-rays. The trajectory line crosses
  the inner cavity. Under interpretation (A) — see plan §"Hollow
  sphere geometry" — the inner surface is a specular reflector with
  reflectivity :math:`\alpha_{\rm in}`, so the particle bounces
  alternately inner ↔ outer. **Closure rank-2** (orbit-space class:
  two-surface), structurally identical to the asymmetric slab with
  curvilinear chord algebra: both have an M/G interval with two
  distinct BC endpoints.

Interpretation (A): :math:`\alpha_{\rm in} = 0` means "particle is
lost to the cavity" (perfect inner absorber); :math:`\alpha_{\rm in}
= 1` means "perfect inner reflector". This makes the hollow sphere a
clean **orbit-space M/G analog** of the asymmetric slab — the
through-ray subset of phase-space carries the same 1-D-interval-with-
two-distinct-BC-endpoints orbit-space structure as slab asymmetric,
only the chord algebra and the G-equivariant lift back to the
higher-dim M differ.

Architecture
------------

Phase-space: :math:`(r, \mu)` with :math:`r \in [R_{\rm in}, R_{\rm
out}]`, signed :math:`\mu \in [-1, 1]`. For each interior point
:math:`(r_i, \mu_q)`:

1. **Compute** :math:`b = r_i\sqrt{1 - \mu_q^2}` and partition the
   trajectory into outer-only (:math:`b > R_{\rm in}`) or through-ray
   (:math:`b \le R_{\rm in}`).

2. **First-leg trajectory integral** along the backward chord to the
   first surface arrival:

   - Outer-only: :math:`L_{\rm first} = r\mu + \sqrt{R_{\rm out}^2 -
     b^2}`, surface = outer (rank-1 sphere formula).
   - Through-ray with :math:`\mu > 0`: backward goes inward; first
     arrival is INNER at :math:`L_{\rm first} = r\mu - \sqrt{R_{\rm
     in}^2 - b^2}`, surface = inner.
   - Through-ray with :math:`\mu < 0`: backward goes outward; first
     arrival is OUTER at :math:`L_{\rm first} = r\mu + \sqrt{R_{\rm
     out}^2 - b^2}`, surface = outer.

3. **Outer-only branch** (:math:`b > R_{\rm in}`): rank-1 closure on
   bounce-period chord :math:`L_p = 2\sqrt{R_{\rm out}^2 - b^2}`
   identical to solid sphere.

4. **Through-ray branch** (:math:`b \le R_{\rm in}`): rank-2 closure
   with single-transit (shell-traversal) optical depth
   :math:`\tau_{\rm step} = \Sigma_t\,(\sqrt{R_{\rm out}^2 - b^2}
   - \sqrt{R_{\rm in}^2 - b^2})`. The two single-transit B integrals
   :math:`B_{\rm in}` (outer→inner shell chord) and :math:`B_{\rm out}`
   (inner→outer shell chord) feed the rank-2 resolvent.

5. **Total angular flux**:
   :math:`\psi(r, \mu) = F + e^{-\Sigma_t L_{\rm first}}\,\psi_{\rm
   surface}` with the appropriate :math:`\psi_{\rm surface}` component.

6. **Scalar flux**: :math:`\phi(r) = 2\pi \int_{-1}^{1} \psi(r, \mu)\,
   \mathrm d\mu`.

Assumptions for Phase 3C-1 prototype
-------------------------------------

- Homogeneous shell (single :math:`\Sigma_t, \Sigma_s, \nu\Sigma_f`);
  multi-region deferred (Phase 3C-2 annulus is the next priority).
- Isotropic scattering.
- Independent reflective specular BCs :math:`\alpha_{\rm in} \in
  [0, 1]` and :math:`\alpha_{\rm out} \in [0, 1]`.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase
  3C-1 hollow sphere plan.
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_hollow_sphere`
  — V_α1_hollow_sph/V_α2_hollow_sph/V_α3_hollow_sph SymPy verifications.
- :mod:`.greens_function_slab_asymmetric` — Phase-3B asymmetric slab
  reference solver (the rank-2 template lifted here to curvilinear
  2-surface).
- :mod:`.greens_function` — solid-sphere reference solver (the
  rank-1 template the outer-only branch reuses).
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: rank-2 BIE block resolvent at 2-surface topologies.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    HollowSphereChordOracle,
)
from orpheus.derivations.continuous.trajectory_resolvent.power_iteration import (
    power_iterate_variant_alpha,
)


@dataclass(frozen=True)
class HollowSphereGreensResult:
    """Result of Variant α power iteration on homogeneous hollow sphere
    with INDEPENDENT specular reflectivities
    :math:`\\alpha_{\\rm in}, \\alpha_{\\rm out}`."""

    k_eff: float
    psi: np.ndarray  # (n_r, n_mu) angular flux on the grid
    phi: np.ndarray  # (n_r,) scalar flux on the radial grid
    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


@dataclass(frozen=True)
class HollowSphereGreensMGResult:
    """Result of multi-group Variant α power iteration on homogeneous
    hollow sphere with INDEPENDENT specular BCs."""

    k_eff: float
    psi_g: np.ndarray   # (G, n_r, n_mu)
    phi_g: np.ndarray   # (G, n_r)
    r_nodes: np.ndarray
    mu_nodes: np.ndarray
    iterations: int
    converged: bool


# ═══════════════════════════════════════════════════════════════════════
# Per-group Variant α operator (hollow sphere, rank-2 + impact-parameter
# partition)
# ═══════════════════════════════════════════════════════════════════════


def _apply_operator_hollow_sphere(
    source_profile: np.ndarray,
    r_nodes: np.ndarray,
    mu_nodes: np.ndarray,
    R_in: float,
    R_out: float,
    sigma_t: float,
    alpha_in: float,
    alpha_out: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group hollow-sphere Variant α operator with rank-2 + rank-1
    closure.

    Thin facade over :class:`HollowSphereChordOracle.apply_operator`
    (the R3 ChordOracle Protocol). The :math:`b`-partition + dual
    rank-1 (outer-only) / rank-2 (through-ray) closure routing live in
    the oracle; this function preserves the legacy call signature.

    Bit-equal with the pre-R3 inlined body.
    """
    oracle = HollowSphereChordOracle(
        r_nodes=r_nodes, mu_nodes=mu_nodes,
        R_in=R_in, R_out=R_out,
        alpha_in=alpha_in, alpha_out=alpha_out,
    )
    return oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj_quad,
    )


# ═══════════════════════════════════════════════════════════════════════
# Scalar-flux extractor — shared with sphere
# ═══════════════════════════════════════════════════════════════════════


def _scalar_flux_from_psi(
    psi: np.ndarray,
    mu_weights: np.ndarray,
) -> np.ndarray:
    r"""Reduce angular flux to scalar flux on the radial grid.

    .. math::

        \phi(r) = 2\pi \int_{-1}^{1} \psi(r, \mu)\,\mathrm d\mu.
    """
    return 2.0 * np.pi * np.einsum('...rm,m->...r', psi, mu_weights)


# ═══════════════════════════════════════════════════════════════════════
# 1G homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_hollow_sphere(
    R_in: float,
    R_out: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha_in: float = 1.0,
    alpha_out: float = 1.0,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> HollowSphereGreensResult:
    r"""Power iteration on the hollow-sphere Variant α operator
    (homogeneous shell, isotropic scattering, rank-2 + outer-only
    closure).

    Solves the k-eigenvalue problem on a homogeneous shell with
    INDEPENDENT specular reflectivities on the inner and outer
    surfaces. Phase-space is partitioned by impact parameter
    :math:`b = r\sqrt{1 - \mu^2}` into outer-only (:math:`b > R_{\rm
    in}`, rank-1 closure) and through-ray (:math:`b \le R_{\rm in}`,
    rank-2 closure) subsets.

    Parameters
    ----------
    R_in : float
        Inner cavity radius (cm). Must be > 0 and < R_out.
    R_out : float
        Outer surface radius (cm).
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G). All in the shell material.
    alpha_in : float, default 1.0
        Inner-surface specular reflectivity in :math:`[0, 1]`.
    alpha_out : float, default 1.0
        Outer-surface specular reflectivity in :math:`[0, 1]`.
    n_r, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (n_r, n_mu) ndarray, optional

    Returns
    -------
    :class:`HollowSphereGreensResult`
    """
    if not (0.0 <= alpha_in <= 1.0):
        raise ValueError(f"alpha_in = {alpha_in} must lie in [0, 1]")
    if not (0.0 <= alpha_out <= 1.0):
        raise ValueError(f"alpha_out = {alpha_out} must lie in [0, 1]")
    if not (0.0 < R_in < R_out):
        raise ValueError(
            f"Need 0 < R_in < R_out; got R_in={R_in}, R_out={R_out}"
        )

    # Radial grid: GL on [R_in, R_out] (open quadrature avoids
    # endpoints which are surfaces).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    r_weights = (R_out - R_in) * 0.5 * r_quad_wts

    # Direction-cosine grid: GL on [-1, 1].
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (n_r, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({n_r}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((n_r, n_mu))

    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; non-absorbing "
            "medium not supported."
        )
    k_eff = nu_sigma_f / sigma_a

    # Volume-integrated fission rate for Rayleigh quotient.
    # Shell volume integral: ∫_{R_in}^{R_out} f(r) · 4π r² dr.
    def fission_rate(phi_r: np.ndarray) -> float:
        return float(
            nu_sigma_f * 4.0 * np.pi
            * np.sum(phi_r * r_nodes ** 2 * r_weights)
        )

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_r = _scalar_flux_from_psi(psi_iter, mu_weights)
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_iter) * phi_r

        psi_new = _apply_operator_hollow_sphere(
            source_profile, r_nodes, mu_nodes, R_in, R_out,
            sigma_t, alpha_in, alpha_out,
            n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(psi_new, mu_weights)

        return psi_new, fission_rate(phi_r), fission_rate(phi_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

    phi_r = _scalar_flux_from_psi(psi, mu_weights)

    return HollowSphereGreensResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi_r,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver (asymmetric BC)
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_hollow_sphere_mg(
    R_in: float,
    R_out: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), [g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,)
    *,
    alpha_in: float = 1.0,
    alpha_out: float = 1.0,
    n_r: int = 24,
    n_mu: int = 24,
    n_traj_quad: int = 64,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> HollowSphereGreensMGResult:
    r"""Multi-group hollow-sphere Variant α power iteration (homogeneous
    shell, isotropic scattering, rank-2 + outer-only closure).

    Multi-group analog of :func:`solve_greens_function_hollow_sphere`.
    Convention: ``sigma_s[g_from, g_to]``. At :math:`\alpha_{\rm in} =
    \alpha_{\rm out} = 1` reduces exactly to
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`.

    Parameters
    ----------
    R_in, R_out : float
    sigma_t : (G,) ndarray
    sigma_s : (G, G) ndarray
    nu_sigma_f : (G,) ndarray
    chi : (G,) ndarray, optional
        Fission spectrum. Default: all-fast emission.
    alpha_in, alpha_out : float
    n_r, n_mu, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (G, n_r, n_mu) ndarray, optional
    initial_k : float, optional

    Returns
    -------
    :class:`HollowSphereGreensMGResult`
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
    if not (0.0 <= alpha_in <= 1.0):
        raise ValueError(f"alpha_in = {alpha_in} must lie in [0, 1]")
    if not (0.0 <= alpha_out <= 1.0):
        raise ValueError(f"alpha_out = {alpha_out} must lie in [0, 1]")
    if not (0.0 < R_in < R_out):
        raise ValueError(
            f"Need 0 < R_in < R_out; got R_in={R_in}, R_out={R_out}"
        )

    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R_in + (R_out - R_in) * 0.5 * (r_quad_pts + 1.0)
    r_weights = (R_out - R_in) * 0.5 * r_quad_wts

    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = mu_quad_pts
    mu_weights = mu_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (G, n_r, n_mu):
            raise ValueError(
                f"initial_psi shape must be ({G}, {n_r}, {n_mu}); "
                f"got {psi.shape}"
            )
    else:
        psi = np.ones((G, n_r, n_mu))

    if initial_k is not None:
        k_eff = float(initial_k)
    else:
        A = np.diag(sigma_t) - sigma_s.T
        F_mat = np.outer(chi, nu_sigma_f)
        eig = np.linalg.eigvals(np.linalg.solve(A, F_mat))
        k_eff = float(np.real(eig).max())
        if k_eff <= 0:
            raise ValueError(f"k_inf estimate non-positive ({k_eff}).")

    def total_fission_rate(phi_g_r: np.ndarray) -> float:
        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g_r)
        return float(
            4.0 * np.pi * np.sum(F_r * r_nodes ** 2 * r_weights)
        )

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_g = 2.0 * np.pi * np.einsum(
            'grm,m->gr', psi_iter, mu_weights,
        )

        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)
        scatter_source = np.einsum('sg,sr->gr', sigma_s, phi_g)
        fission_source = (chi[:, None] / k_iter) * F_r[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi_iter)
        for g in range(G):
            psi_new[g] = _apply_operator_hollow_sphere(
                source_profile_g[g], r_nodes, mu_nodes, R_in, R_out,
                float(sigma_t[g]), alpha_in, alpha_out,
                n_traj_quad=n_traj_quad,
            )

        phi_g_new = 2.0 * np.pi * np.einsum(
            'grm,m->gr', psi_new, mu_weights,
        )
        return psi_new, total_fission_rate(phi_g), total_fission_rate(phi_g_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

    phi_g = 2.0 * np.pi * np.einsum('grm,m->gr', psi, mu_weights)

    return HollowSphereGreensMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        r_nodes=r_nodes,
        mu_nodes=mu_nodes,
        iterations=iterations,
        converged=converged,
    )
