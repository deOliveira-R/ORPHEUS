r"""Phase-1 cylinder **Variant α** Green's function reference (homogeneous
infinite cylinder, specular/vacuum/partial-reflection BC) — research-grade
prototype.

Mirrors the sphere implementation in :mod:`.greens_function`. Iterates the
angle-resolved Green's function along **bouncing characteristics** in a
3D phase-space :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`, with the
axial direction kept **explicit** (Issue #129 discipline) — the cylinder
Nyström kernel pre-integrates axial which causes the 22 % planar-limit
mismatch; Variant α avoids this by retaining the full angular
distribution.

Cylinder Variant α architecture
--------------------------------

For each phase-space grid point :math:`(r_i, \mu_{q,\rm axial},
\varphi_{p,\rm az})`:

1. **First-leg trajectory integral** along the in-plane backward chord
   from :math:`r_i` to first surface arrival, with 3D path length
   :math:`L_0 = L_{\rm 2D, first} / s_{\rm in\!-\!plane}`,
   :math:`s_{\rm in\!-\!plane} = \sqrt{1 - \mu_{\rm axial}^2}`:

   .. math::

       L_{\rm 2D, first}(r, \varphi_{\rm az}) =
           r\,\cos\varphi_{\rm az} +
           \sqrt{R^2 - r^2 \sin^2\varphi_{\rm az}}.

   The integrand :math:`q(r'(s_{\rm 2D}))\,e^{-\Sigma_t s_{\rm 3D}}`
   is parametrised by the 2D arclength :math:`s_{\rm 2D} \in [0,
   L_{\rm 2D, first}]`; the radial position along the chord is
   :math:`r'(s_{\rm 2D}) = \sqrt{r^2 - 2 r s_{\rm 2D} \cos\varphi_{\rm
   az} + s_{\rm 2D}^2}` and :math:`s_{\rm 3D} = s_{\rm 2D}/s_{\rm
   in\!-\!plane}`.

2. **Bounce-period integral** along the antipodal in-plane chord at
   conserved impact parameter :math:`b = r\,|\sin\varphi_{\rm az}|`:

   .. math::

       L_{\rm 2D, period} = 2\sqrt{R^2 - b^2}, \qquad
       L_{\rm period} = L_{\rm 2D, period}/s_{\rm in\!-\!plane},

   .. math::

       r_{\rm chord}(s_{\rm 2D}) =
           \sqrt{b^2 + (s_{\rm 2D} - L_{\rm 2D, period}/2)^2}.

3. **Geometric bounce sum** (closed form):

   .. math::

       \psi_{\rm surf} = \frac{\alpha\,B}
                              {1 - \alpha\,e^{-\Sigma_t L_{\rm period}}}.

4. **Total angular flux**: :math:`\psi(r_i, \mu_{\rm axial},
   \varphi_{\rm az}) = F + e^{-\Sigma_t L_0}\,\psi_{\rm surf}`.

5. **Scalar flux**: :math:`\phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm
   axial}\!\int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\,\psi(r, \mu_{\rm
   axial}, \varphi_{\rm az})`.

Symmetry exploited: the closed-cylinder isotropic eigenmode is
symmetric under :math:`\varphi_{\rm az} \to -\varphi_{\rm az}` and
:math:`\mu_{\rm axial} \to -\mu_{\rm axial}`. We discretise both
ranges fully — the redundancy is harmless and supports vacuum BC
where the symmetry is lost.

Grazing-ray pathology
---------------------

Two grazing-ray loci exist for cylinder Variant α:

- **Axial grazing**: :math:`\mu_{\rm axial} \to \pm 1` makes
  :math:`s_{\rm in\!-\!plane} \to 0` and :math:`L_{\rm period} \to
  \infty`. The bounce-sum geometric series :math:`1 / (1 - \alpha
  e^{-\Sigma_t L_{\rm period}}) \to 1` finite; :math:`B \to 0` faster
  via the :math:`L_{\rm 2D, period}/s_{\rm in\!-\!plane}` weight on
  the integrand (whose magnitude is bounded by :math:`(q/\Sigma_t) \cdot
  s_{\rm in\!-\!plane}` after the :math:`\Sigma_t s_{\rm 3D}` exponential
  collapses the integration support). Net contribution at fixed
  :math:`\varphi_{\rm az}` is finite. We use Gauss-Legendre on
  :math:`\mu_{\rm axial} \in [-1, 1]` which avoids the endpoints
  exactly.

- **Tangential in-plane** (:math:`\varphi_{\rm az} \to \pm\pi/2`):
  trajectory grazes the surface at large :math:`b \to R`, so
  :math:`L_{\rm 2D, period} \to 0` and :math:`L_{\rm period} \to 0`.
  The bounce-sum factor :math:`1/(1 - \alpha e^{-\Sigma_t L_{\rm
  period}}) \to 1/(1 - \alpha)` which **diverges at** :math:`\alpha =
  1` (closed cylinder). Mitigated by the analogous :math:`B \to 0`
  scaling — :math:`\psi_{\rm surf} = \alpha B / (1 - \alpha
  e^{-\Sigma_t L_{\rm period}}) \to 0/0` form, with V_α1 ensuring the
  ratio limits to :math:`q/\Sigma_t` finite. We use Gauss-Legendre on
  :math:`\varphi_{\rm az} \in [0, 2\pi)` which never lands on
  :math:`\pm\pi/2` exactly.

Assumptions for Phase 1 prototype
----------------------------------

- Homogeneous cylinder (single :math:`\Sigma_t`, :math:`\Sigma_s`,
  :math:`\nu\Sigma_f`); multi-region deferred to Phase 1b.
- Isotropic scattering.
- Specular BC parametrised by :math:`\alpha \in [0, 1]`.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder`
  — V_α1_cyl/V_α2_cyl/V_α3_cyl SymPy verifications.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 1
  cylinder Variant α plan.
- Sphere reference solver:
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.greens_function`
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from orpheus.derivations.continuous.trajectory_resolvent.chord_oracle import (
    CylinderChordOracle,
)
from orpheus.derivations.continuous.trajectory_resolvent.power_iteration import (
    power_iterate_variant_alpha,
)


@dataclass(frozen=True)
class CylinderGreensResult:
    """Result of Variant α power iteration on homogeneous cylinder."""

    k_eff: float
    psi: np.ndarray  # (n_r, n_mu_axial, n_phi_az)
    phi: np.ndarray  # (n_r,)
    r_nodes: np.ndarray
    mu_axial_nodes: np.ndarray
    phi_az_nodes: np.ndarray
    iterations: int
    converged: bool


@dataclass(frozen=True)
class CylinderGreensMGResult:
    """Result of multi-group Variant α power iteration on cylinder."""

    k_eff: float
    psi_g: np.ndarray  # (G, n_r, n_mu_axial, n_phi_az)
    phi_g: np.ndarray  # (G, n_r)
    r_nodes: np.ndarray
    mu_axial_nodes: np.ndarray
    phi_az_nodes: np.ndarray
    iterations: int
    converged: bool


# ═══════════════════════════════════════════════════════════════════════
# Geometry helpers (in-plane chord arithmetic)
# ═══════════════════════════════════════════════════════════════════════


def _first_leg_2d_chord(
    r: float, phi_az: float, R: float,
) -> float:
    r"""Compute the 2D in-plane backward chord length from interior
    point :math:`r` (with azimuth :math:`\varphi_{\rm az}`) to the
    cylinder surface at :math:`r = R`.

    .. math::

        L_{\rm 2D, first}(r, \varphi_{\rm az}) =
            r\,\cos\varphi_{\rm az} +
            \sqrt{R^2 - r^2 \sin^2\varphi_{\rm az}}.

    Parameters
    ----------
    r : float
        Interior radial position in :math:`[0, R]`.
    phi_az : float
        Azimuth of the in-plane velocity component measured from
        :math:`\hat r`.
    R : float
        Cylinder outer radius.
    """
    sin_phi = np.sin(phi_az)
    cos_phi = np.cos(phi_az)
    disc = R * R - r * r * sin_phi * sin_phi
    return r * cos_phi + np.sqrt(max(disc, 0.0))


def _impact_parameter(r: float, phi_az: float) -> float:
    r"""Conserved impact parameter for cylinder specular trajectories.

    :math:`b = r\,|\sin\varphi_{\rm az}|`.
    """
    return r * abs(np.sin(phi_az))


def _bounce_period_2d_chord(b: float, R: float) -> float:
    r"""2D in-plane bounce-period chord (full antipodal chord on
    in-plane circle).

    :math:`L_{\rm 2D, period}(b) = 2\sqrt{R^2 - b^2}`. At :math:`b
    \ge R` returns 0 (trajectory never re-enters).
    """
    if b >= R:
        return 0.0
    return 2.0 * np.sqrt(R * R - b * b)


# ═══════════════════════════════════════════════════════════════════════
# Per-group Variant α operator (cylinder)
# ═══════════════════════════════════════════════════════════════════════


def _apply_operator_cylinder(
    source_profile: np.ndarray,
    r_nodes: np.ndarray,
    mu_axial_nodes: np.ndarray,
    phi_az_nodes: np.ndarray,
    R: float,
    sigma_t: float,
    alpha: float,
    *,
    n_traj_quad: int,
) -> np.ndarray:
    r"""Per-group cylinder Variant α operator.

    Thin facade over :class:`CylinderChordOracle.apply_operator` (the
    R3 ChordOracle Protocol). The 2D chord arithmetic + 3D arclength
    lift + rank-1 closure live in the oracle; this function preserves
    the legacy call signature for back-compatibility.

    Bit-equal with the pre-R3 inlined body — every FP operation runs
    in the same order.
    """
    oracle = CylinderChordOracle(
        r_nodes=r_nodes,
        mu_axial_nodes=mu_axial_nodes,
        phi_az_nodes=phi_az_nodes,
        R=R, alpha=alpha,
    )
    return oracle.apply_operator(
        source_profile, sigma_t=sigma_t, n_traj_quad=n_traj_quad,
    )


# ═══════════════════════════════════════════════════════════════════════
# Scalar-flux extractor
# ═══════════════════════════════════════════════════════════════════════


def _scalar_flux_from_psi(
    psi: np.ndarray,
    mu_axial_weights: np.ndarray,
    phi_az_weights: np.ndarray,
) -> np.ndarray:
    r"""Reduce angular flux to scalar flux on the radial grid.

    .. math::

        \phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm axial}\!
                    \int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\,
                    \psi(r, \mu_{\rm axial}, \varphi_{\rm az}).

    Parameters
    ----------
    psi : (..., n_r, n_mu, n_phi) ndarray
        Angular flux. Leading axes are passed through.
    mu_axial_weights : (n_mu,) ndarray
    phi_az_weights : (n_phi,) ndarray

    Returns
    -------
    (..., n_r) ndarray
    """
    return np.einsum(
        '...rmp,m,p->...r', psi, mu_axial_weights, phi_az_weights,
    )


# ═══════════════════════════════════════════════════════════════════════
# 1G homogeneous solver
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_cylinder(
    R: float,
    sigma_t: float,
    sigma_s: float,
    nu_sigma_f: float,
    *,
    alpha: float = 1.0,
    n_r: int = 16,
    n_mu_axial: int = 12,
    n_phi_az: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 200,
    tol: float = 1e-10,
    initial_psi: np.ndarray | None = None,
) -> CylinderGreensResult:
    r"""Power iteration on the cylinder Variant α operator (homogeneous,
    isotropic scattering, specular/vacuum/partial-reflection BC).

    Solves the k-eigenvalue problem on an infinite homogeneous cylinder
    of radius :math:`R` via fission-source iteration with the cylinder
    Variant α operator (angle-resolved Green's function with bounces
    summed analytically).

    Boundary condition is parametrised by :math:`\alpha \in [0, 1]`:

    - :math:`\alpha = 1`: perfect specular (closed cylinder, no
      leakage). :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`
      EXACTLY by V_α1_cyl. The load-bearing Phase-1 acceptance test.
    - :math:`\alpha = 0`: vacuum. Spatial eigenmode peaked at center,
      depleted at surface. :math:`k_{\rm eff} < k_\infty`.
    - :math:`0 < \alpha < 1`: partial-reflection albedo.

    Parameters
    ----------
    R : float
        Outer radius.
    sigma_t, sigma_s, nu_sigma_f : float
        Group cross sections (1G).
    alpha : float, default 1.0
        Surface reflectivity. 1 = closed cylinder (specular), 0 = vacuum.
    n_r : int
        Radial Gauss-Legendre quadrature order on :math:`(0, R)`.
    n_mu_axial : int
        Axial-cosine Gauss-Legendre order on :math:`[-1, 1]`.
    n_phi_az : int
        Azimuthal Gauss-Legendre order on :math:`[0, 2\pi)`.
    n_traj_quad : int
        Trajectory + bounce-period quadrature order.
    max_iter, tol : int, float
        Power-iteration parameters.
    initial_psi : (n_r, n_mu_axial, n_phi_az) ndarray, optional

    Returns
    -------
    :class:`CylinderGreensResult`
    """
    if not (0.0 <= alpha <= 1.0):
        raise ValueError(f"alpha = {alpha} must lie in [0, 1]")

    # Radial grid: GL on (0, R).
    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    r_weights = R * 0.5 * r_quad_wts

    # Axial-cosine grid: GL on [-1, 1].
    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu_axial)
    mu_axial_nodes = mu_quad_pts
    mu_axial_weights = mu_quad_wts

    # Azimuthal grid: GL on [0, 2π) — open quadrature avoids the
    # tangential ±π/2 grazing rays.
    phi_quad_pts, phi_quad_wts = np.polynomial.legendre.leggauss(n_phi_az)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)
    phi_az_weights = np.pi * phi_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (n_r, n_mu_axial, n_phi_az):
            raise ValueError(
                f"initial_psi shape must be ({n_r}, {n_mu_axial}, "
                f"{n_phi_az}); got {psi.shape}"
            )
    else:
        psi = np.ones((n_r, n_mu_axial, n_phi_az))

    # Initial k_eff.
    sigma_a = sigma_t - sigma_s
    if sigma_a <= 0:
        raise ValueError(
            f"sigma_a = sigma_t - sigma_s = {sigma_a} ≤ 0; non-absorbing "
            "medium not supported."
        )
    k_eff = nu_sigma_f / sigma_a

    # Volume-integrated fission rate per unit cylinder length:
    # ∫_0^R νΣ_f φ(r) · 2π r dr.
    def fission_rate(phi_r: np.ndarray) -> float:
        return float(
            nu_sigma_f * 2.0 * np.pi
            * np.sum(phi_r * r_nodes * r_weights)
        )

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        # Source profile q(r)/(4π) = (σ_s + νσ_f/k)/(4π) · φ(r).
        phi_r = _scalar_flux_from_psi(
            psi_iter, mu_axial_weights, phi_az_weights,
        )
        source_profile = inv_4pi * (sigma_s + nu_sigma_f / k_iter) * phi_r

        psi_new = _apply_operator_cylinder(
            source_profile, r_nodes, mu_axial_nodes, phi_az_nodes, R,
            sigma_t, alpha, n_traj_quad=n_traj_quad,
        )

        phi_new = _scalar_flux_from_psi(
            psi_new, mu_axial_weights, phi_az_weights,
        )

        return psi_new, fission_rate(phi_r), fission_rate(phi_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

    phi_r = _scalar_flux_from_psi(psi, mu_axial_weights, phi_az_weights)

    return CylinderGreensResult(
        k_eff=float(k_eff),
        psi=psi,
        phi=phi_r,
        r_nodes=r_nodes,
        mu_axial_nodes=mu_axial_nodes,
        phi_az_nodes=phi_az_nodes,
        iterations=iterations,
        converged=converged,
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group homogeneous solver
# ═══════════════════════════════════════════════════════════════════════


def solve_greens_function_cylinder_mg(
    R: float,
    sigma_t: np.ndarray,        # (G,)
    sigma_s: np.ndarray,        # (G, G), [g_from, g_to]
    nu_sigma_f: np.ndarray,     # (G,)
    chi: np.ndarray | None = None,  # (G,)
    *,
    alpha: float = 1.0,
    n_r: int = 16,
    n_mu_axial: int = 12,
    n_phi_az: int = 24,
    n_traj_quad: int = 32,
    max_iter: int = 300,
    tol: float = 1e-9,
    initial_psi: np.ndarray | None = None,
    initial_k: float | None = None,
) -> CylinderGreensMGResult:
    r"""Multi-group cylinder Variant α power iteration (homogeneous,
    isotropic scattering, parametrised by :math:`\alpha`).

    Multi-group analog of :func:`solve_greens_function_cylinder`.
    Convention: ``sigma_s[g_from, g_to]``.

    At :math:`\alpha = 1` (closed cylinder) reduces exactly to
    :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` —
    this is the load-bearing MG verification.

    Parameters
    ----------
    R : float
        Outer radius.
    sigma_t : (G,) ndarray
    sigma_s : (G, G) ndarray
    nu_sigma_f : (G,) ndarray
    chi : (G,) ndarray, optional
        Fission spectrum. Default: all-fast emission.
    alpha : float, default 1.0
    n_r, n_mu_axial, n_phi_az, n_traj_quad : int
    max_iter, tol : int, float
    initial_psi : (G, n_r, n_mu_axial, n_phi_az) ndarray, optional
    initial_k : float, optional

    Returns
    -------
    :class:`CylinderGreensMGResult`
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

    r_quad_pts, r_quad_wts = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_quad_pts + 1.0)
    r_weights = R * 0.5 * r_quad_wts

    mu_quad_pts, mu_quad_wts = np.polynomial.legendre.leggauss(n_mu_axial)
    mu_axial_nodes = mu_quad_pts
    mu_axial_weights = mu_quad_wts

    phi_quad_pts, phi_quad_wts = np.polynomial.legendre.leggauss(n_phi_az)
    phi_az_nodes = np.pi * (phi_quad_pts + 1.0)
    phi_az_weights = np.pi * phi_quad_wts

    if initial_psi is not None:
        psi = np.asarray(initial_psi, dtype=float).copy()
        if psi.shape != (G, n_r, n_mu_axial, n_phi_az):
            raise ValueError(
                f"initial_psi shape must be ({G}, {n_r}, {n_mu_axial}, "
                f"{n_phi_az}); got {psi.shape}"
            )
    else:
        psi = np.ones((G, n_r, n_mu_axial, n_phi_az))

    if initial_k is not None:
        k_eff = float(initial_k)
    else:
        # Use kinf_homogeneous as initial guess.
        A = np.diag(sigma_t) - sigma_s.T
        F = np.outer(chi, nu_sigma_f)
        eig = np.linalg.eigvals(np.linalg.solve(A, F))
        k_eff = float(np.real(eig).max())
        if k_eff <= 0:
            raise ValueError(
                f"k_inf estimate non-positive ({k_eff})."
            )

    def total_fission_rate(phi_g_r: np.ndarray) -> float:
        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g_r)
        return float(2.0 * np.pi * np.sum(F_r * r_nodes * r_weights))

    inv_4pi = 1.0 / (4.0 * np.pi)

    def _step(psi_iter, k_iter):
        phi_g = np.einsum(
            'grmp,m,p->gr', psi_iter, mu_axial_weights, phi_az_weights,
        )  # (G, n_r)

        F_r = np.einsum('g,gr->r', nu_sigma_f, phi_g)  # (n_r,)
        scatter_source = np.einsum('sg,sr->gr', sigma_s, phi_g)  # (G, n_r)
        fission_source = (chi[:, None] / k_iter) * F_r[None, :]
        source_profile_g = inv_4pi * (scatter_source + fission_source)

        psi_new = np.zeros_like(psi_iter)
        for g in range(G):
            psi_new[g] = _apply_operator_cylinder(
                source_profile_g[g], r_nodes, mu_axial_nodes,
                phi_az_nodes, R, float(sigma_t[g]), alpha,
                n_traj_quad=n_traj_quad,
            )

        phi_g_new = np.einsum(
            'grmp,m,p->gr', psi_new, mu_axial_weights, phi_az_weights,
        )
        return psi_new, total_fission_rate(phi_g), total_fission_rate(phi_g_new)

    pi_result = power_iterate_variant_alpha(
        _step, psi, initial_k=k_eff, max_iter=max_iter, tol=tol,
    )
    psi = pi_result.psi
    k_eff = pi_result.k_eff
    iterations = pi_result.iterations
    converged = pi_result.converged

    phi_g = np.einsum(
        'grmp,m,p->gr', psi, mu_axial_weights, phi_az_weights,
    )

    return CylinderGreensMGResult(
        k_eff=float(k_eff),
        psi_g=psi,
        phi_g=phi_g,
        r_nodes=r_nodes,
        mu_axial_nodes=mu_axial_nodes,
        phi_az_nodes=phi_az_nodes,
        iterations=iterations,
        converged=converged,
    )
