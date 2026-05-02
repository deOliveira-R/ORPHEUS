r"""Phase-3A slab **Variant α** Green's function reference (homogeneous
1D slab, symmetric reflective specular BC :math:`\alpha_{\rm left} =
\alpha_{\rm right} = \alpha`) — research-grade prototype.

Mirrors the sphere implementation in :mod:`.greens_function` and the
cylinder implementation in :mod:`.greens_function_cylinder`. Iterates
the angle-resolved Green's function along **bouncing characteristics**
in a 2D phase-space :math:`(x, \mu)` with :math:`x \in [0, L]` and
signed :math:`\mu \in [-1, 1]`.

ERR-035 fix (2026-05-02) — first-principles rank-2 closure
-----------------------------------------------------------

The original Phase-3A symmetric closure used a heuristic
generalisation of the sphere/cylinder rank-1 form:

    ψ_surf = α · B_period / (1 - α² · e^{-2τ})

(with `B_period` being the full out-and-back bounce-period source
integral and `α²` = per-period reflection product for slab's 2-bounce
geometry). This formula coincides with the first-principles closure
ONLY at α∈{0, 1} corners. At intermediate α it under-estimated
:math:`\psi_{\rm surf}` by ~1.3e-4 relative on fuel-A τ_L=5 (caught
by the Phase-3B reduce-to-symmetric consistency check at α=0.5).

The first-principles rank-2 closure at :math:`\alpha_L = \alpha_R =
\alpha` is

    ψ_surf = α · B_single / (1 - α · e^{-τ})

(single-transit B integral, single-transit denominator). This module
now delegates the operator application to the Phase-3B asymmetric
solver in :mod:`.greens_function_slab_asymmetric` with :math:`\alpha_L
= \alpha_R = \alpha`. The two formulations are mathematically
equivalent on uniform-source / closed-slab cases (V_α1_slab, vacuum
α=0), and the rank-2 path is correct at intermediate α where the
heuristic was not.

Architectural consequence: this module exposes the same public API
(``solve_greens_function_slab``, ``solve_greens_function_slab_mg``,
``SlabGreensResult``, ``SlabGreensMGResult``) for back-compatibility
of all existing callers and tests. The internals are thin wrappers
over the rank-2 path; the heuristic ``_apply_operator_slab`` and the
``_bounce_period_chord_slab`` / ``_first_leg_chord_slab`` helpers
have been removed.

Slab Variant α architecture (rank-2 framing)
---------------------------------------------

For each phase-space grid point :math:`(x_i, \mu_q)` the rank-2
asymmetric solver evaluates:

1. **First-leg trajectory integral** along the backward chord from
   :math:`x_i` to the first surface arrival.

2. **Single-transit B integrals** (rank-2 monodromy):

   .. math::

      B_{LR}(\mu) &= \int_0^{L/|\mu|}
                       q(0 + s\cdot|\mu|)\,e^{-\Sigma_t s}\,
                       \mathrm d s \\
      B_{RL}(\mu) &= \int_0^{L/|\mu|}
                       q(L - s\cdot|\mu|)\,e^{-\Sigma_t s}\,
                       \mathrm d s.

3. **Rank-2 surface closure** at :math:`\alpha_L = \alpha_R = \alpha`:

   .. math::

      \psi_{\rm surf}(\mu) = \frac{\alpha\,B(\mu)}
                                  {1 - \alpha\,e^{-\Sigma_t L/|\mu|}}

   (single-transit denominator, structurally equivalent to the rank-2
   matrix inversion at the symmetric corner).

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

At :math:`\mu \to 0` the transit chord :math:`L_{\rm transit} \to
\infty`; :math:`e^{-\Sigma_t L_{\rm transit}} \to 0` and the rank-2
denominator :math:`(1 - \alpha\,e^{-\Sigma_t L_{\rm transit}}) \to 1`.
The single-transit B integrand picks up the integrable :math:`1/|\mu|`
weight from the chord length but the leading :math:`e^{-\Sigma_t
L_{\rm first}}` attenuation in the closure (:math:`L_{\rm first} \to
\infty` as well) drives :math:`\psi_{\rm surf} \to 0` exponentially.
**Slab is structurally immune to the sphere's :math:`\mu \to 0`
Hadamard divergence** — the planar geometry ensures :math:`L_{\rm
first}` and :math:`L_{\rm transit}` co-diverge.

We use Gauss-Legendre on :math:`\mu \in [-1, 1]` (open at endpoints)
which avoids the singular :math:`\mu = 0` point exactly.

Assumptions for Phase 3A prototype
-----------------------------------

- Homogeneous slab (single :math:`\Sigma_t`, :math:`\Sigma_s`,
  :math:`\nu\Sigma_f`); multi-region deferred.
- Isotropic scattering.
- **Symmetric** reflective specular BC :math:`\alpha_{\rm left} =
  \alpha_{\rm right} = \alpha \in [0, 1]`. Asymmetric BC is provided
  directly by the rank-2 module
  :mod:`.greens_function_slab_asymmetric`.

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
- ERR-035 closeout memo:
  ``.claude/agent-memory/method-implementer/err035_phase3a_delegation_fix.md``.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.derivations.continuous.peierls_greens_function.greens_function_slab_asymmetric import (
    solve_greens_function_slab_asymmetric,
    solve_greens_function_slab_asymmetric_mg,
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
# 1G homogeneous solver — delegates to rank-2 at α_L = α_R = α
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
    characteristics summed analytically).

    Boundary condition is parametrised by :math:`\alpha \in [0, 1]`
    applied SYMMETRICALLY to both walls:

    - :math:`\alpha = 1`: closed slab (perfect specular both walls,
      no leakage). :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`
      EXACTLY by V_α1_slab. Load-bearing Phase-3A acceptance test.
    - :math:`\alpha = 0`: vacuum on both walls. Spatial eigenmode
      peaked at center, depleted at walls. :math:`k_{\rm eff} <
      k_\infty`.
    - :math:`0 < \alpha < 1`: partial-reflection albedo on both walls.

    Implementation note (ERR-035 fix, 2026-05-02): this function
    delegates to :func:`solve_greens_function_slab_asymmetric` with
    :math:`\alpha_L = \alpha_R = \alpha`, replacing the original
    Phase-3A heuristic closure (which was wrong at intermediate α —
    see module docstring). The delegated rank-2 closure is the
    first-principles correct form.

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
    res_asym = solve_greens_function_slab_asymmetric(
        L=L,
        sigma_t=sigma_t,
        sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f,
        alpha_left=alpha,
        alpha_right=alpha,
        n_x=n_x,
        n_mu=n_mu,
        n_traj_quad=n_traj_quad,
        max_iter=max_iter,
        tol=tol,
        initial_psi=initial_psi,
    )
    return SlabGreensResult(
        k_eff=res_asym.k_eff,
        psi=res_asym.psi,
        phi=res_asym.phi,
        x_nodes=res_asym.x_nodes,
        mu_nodes=res_asym.mu_nodes,
        iterations=res_asym.iterations,
        converged=res_asym.converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver — delegates to rank-2 at α_L = α_R = α
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

    Implementation note (ERR-035 fix, 2026-05-02): delegates to
    :func:`solve_greens_function_slab_asymmetric_mg` with
    :math:`\alpha_L = \alpha_R = \alpha`.

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
    res_asym = solve_greens_function_slab_asymmetric_mg(
        L=L,
        sigma_t=sigma_t,
        sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f,
        chi=chi,
        alpha_left=alpha,
        alpha_right=alpha,
        n_x=n_x,
        n_mu=n_mu,
        n_traj_quad=n_traj_quad,
        max_iter=max_iter,
        tol=tol,
        initial_psi=initial_psi,
        initial_k=initial_k,
    )
    return SlabGreensMGResult(
        k_eff=res_asym.k_eff,
        psi_g=res_asym.psi_g,
        phi_g=res_asym.phi_g,
        x_nodes=res_asym.x_nodes,
        mu_nodes=res_asym.mu_nodes,
        iterations=res_asym.iterations,
        converged=res_asym.converged,
    )
