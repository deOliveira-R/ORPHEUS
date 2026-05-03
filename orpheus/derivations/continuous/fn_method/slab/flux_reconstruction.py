r"""Interior flux reconstruction for the bare-critical slab via the
Kaper-Lindeman-Leaf (KLL) 1974 Fredholm-iteration recipe.

This module ships the **rich-machinery extension** of the F_N
boundary-only solver (:mod:`.one_group`):

* :func:`solve_kll_slab_continuum_coefficient` — solves the Fredholm
  integral equation for :math:`A(\nu)` via successive iteration
  (KLL Eqs. 6 + 7).
* :func:`slab_scalar_flux_kll` — evaluates the converged scalar flux
  :math:`\phi(z)` interior to the slab via KLL Eq. 7.
* :func:`slab_angular_flux_from_scalar` — reconstructs the angular
  flux :math:`\psi(z, \mu)` via characteristic integration of the BTE
  with the converged :math:`\phi(z)` as the source.

Why KLL (and not just F_N continued)
------------------------------------

The F_N method gives the boundary angular flux :math:`\psi(\pm a,
\mp\mu) = \sum_\alpha a_\alpha \mu^\alpha` directly. Reconstructing
the **interior** flux from the boundary in the F_N framework requires
inverting the Case full-range expansion, which boils down to the same
Fredholm integral equation KLL solves. Hence we just call KLL
directly: it is the canonical published reference for slab interior
flux, structurally independent of F_N collocation, with
benchmark-quality reference values (KLL Tables III, IV).

The F_N solver still drives this module: it provides :math:`r_c` (the
critical half-thickness) which is the fixed parameter of the KLL
Fredholm equation. In the verification chain:

1. F_N gives :math:`r_c` (boundary representation, half-range basis).
2. KLL gives :math:`\phi(z)` at that :math:`r_c` (interior
   representation, full-range Case expansion).
3. The agreement of both with the analytical critical condition
   :math:`r_c \cdot \nu_0\,\mathrm{trig}(...) = 0` is the structural
   independence proof for the medium-only scalar :math:`r_c`.
4. The flux ratios :math:`\phi(z)/\phi(0)` from KLL agree with the
   benchmark tables KLL III + Sood Table 14 to ≤ 1e-5.

References
----------

* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94, Eqs. 5-7.
* Mitsis 1963, *Nucl. Sci. Eng.* **17**, 55 + ANL-6787 (Eqs. 13-22 of
  the Argonne report — the integral equation derivation).
* Case & Zweifel 1967, *Linear Transport Theory*, §6.6.A.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

from ..core.dispersion import case_nu0


# ─────────────────────────────────────────────────────────────────────
# Wiener-Hopf X-function and related medium-only quantities
# ─────────────────────────────────────────────────────────────────────


def _theta_kll(mu: float, c: float) -> float:
    r"""KLL Eq. 5 (paragraph after): :math:`\theta(\mu)` for :math:`\mu
    \in (0, 1)`.

    .. math::

       \theta(\mu) = \arg\!\left(1 - \tfrac{1}{2}\,c\,\mu\,
       \log\frac{1 + \mu}{1 - \mu} + \tfrac{1}{2} i \pi c \mu\right).

    For :math:`c < 1` this is the argument of the Case
    :math:`\Lambda^{+}(\mu)` (limiting value of the dispersion function
    on the upper side of the cut). For :math:`c > 1` (multiplying
    medium) the same expression continues to define :math:`\theta`
    used in :math:`\Gamma(z)` and :math:`X(z)`.
    """
    if mu <= 0.0 or mu >= 1.0:
        raise ValueError(f"mu must be in (0, 1), got {mu}")
    real_part = 1.0 - 0.5 * c * mu * math.log((1.0 + mu) / (1.0 - mu))
    imag_part = 0.5 * math.pi * c * mu
    return math.atan2(imag_part, real_part)


def _Gamma_real_imag(c: float, u0: float) -> tuple[float, float]:
    r"""Compute :math:`\Gamma(i u_0) = \Gamma_{10} + i \Gamma_{20}`
    per KLL Eqs. 7-8 (after the introduction of :math:`u_0 = |\nu_0|`).

    .. math::

       \Gamma_{10} = \frac{1}{\pi} \int_0^1 \frac{\mu}{\mu^2 + u_0^2}\,
       \theta(\mu)\,d\mu, \qquad
       \Gamma_{20} = \frac{1}{\pi} \int_0^1 \frac{u_0}{\mu^2 + u_0^2}\,
       \theta(\mu)\,d\mu .

    These are the real and imaginary parts of :math:`\Gamma(z)` at
    :math:`z = i u_0`; together they give
    :math:`X(i u_0) = X_{10} + i X_{20}` via KLL.
    """

    def integrand_re(mu: float) -> float:
        return (mu / (mu * mu + u0 * u0)) * _theta_kll(mu, c)

    def integrand_im(mu: float) -> float:
        return (u0 / (mu * mu + u0 * u0)) * _theta_kll(mu, c)

    Gamma_10, _ = quad(integrand_re, 0.0, 1.0, limit=200)
    Gamma_20, _ = quad(integrand_im, 0.0, 1.0, limit=200)
    return Gamma_10 / math.pi, Gamma_20 / math.pi


def _X10_X20(c: float, u0: float) -> tuple[float, float]:
    r"""Compute :math:`X_{10}, X_{20}` (real and imaginary parts of
    :math:`X(i u_0)`) per KLL.

    .. math::

       X_{10} = \frac{e^{\Gamma_{10}}}{1 + u_0^2}\,
       (\cos\Gamma_{20} - u_0 \sin\Gamma_{20}), \qquad
       X_{20} = \frac{e^{\Gamma_{10}}}{1 + u_0^2}\,
       (\sin\Gamma_{20} + u_0 \cos\Gamma_{20}) .

    These factor :math:`X(i u_0)` from the dispersion-function argument
    integral.
    """
    Gamma_10, Gamma_20 = _Gamma_real_imag(c, u0)
    pref = math.exp(Gamma_10) / (1.0 + u0 * u0)
    X10 = pref * (math.cos(Gamma_20) - u0 * math.sin(Gamma_20))
    X20 = pref * (math.sin(Gamma_20) + u0 * math.cos(Gamma_20))
    return X10, X20


def _X_negative_real_axis(mu: float, c: float, u0: float) -> float:
    r"""Compute :math:`X(-\mu)` for :math:`\mu \in (0, 1)`.

    KLL: :math:`X(-\mu) = (1 + \mu)^{-1} \exp[\Gamma(-\mu)]`, where

    .. math::

       \Gamma(-\mu) = \frac{1}{\pi}\int_0^1 \theta(\tau)\,
       \frac{d\tau}{\tau + \mu} .

    This is :math:`X(z)` evaluated below the cut. The value is real
    and positive for :math:`\mu \in (0, 1)` and :math:`c > 1`.
    """

    def integrand(tau: float) -> float:
        return _theta_kll(tau, c) / (tau + mu)

    integral, _ = quad(integrand, 0.0, 1.0, limit=200)
    Gamma_minus_mu = integral / math.pi
    return math.exp(Gamma_minus_mu) / (1.0 + mu)


def _gamma_kll(mu: float, c: float, u0: float) -> float:
    r"""KLL: :math:`\gamma(\mu) = \tfrac{1}{2} c (c - 1) \mu X(-\mu) (\mu^2 +
    u_0^2)\,g(\mu)` for :math:`c > 1`.

    Where :math:`g(\mu) = 1/[(1 - \tfrac{1}{2} c \mu \log((1+\mu)/(1-\mu)))^2
    + (\tfrac{\pi c \mu}{2})^2]` is :math:`|\Lambda^{+}(\mu)|^{-2}`.
    """
    if mu <= 0.0 or mu >= 1.0:
        raise ValueError(f"mu must be in (0, 1), got {mu}")
    real_part = 1.0 - 0.5 * c * mu * math.log((1.0 + mu) / (1.0 - mu))
    imag_part = 0.5 * math.pi * c * mu
    g = 1.0 / (real_part * real_part + imag_part * imag_part)
    X_neg_mu = _X_negative_real_axis(mu, c, u0)
    return 0.5 * c * (c - 1.0) * mu * X_neg_mu * (mu * mu + u0 * u0) * g


# ─────────────────────────────────────────────────────────────────────
# KLL Fredholm iteration for A(ν)
# ─────────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class KLLSlabFluxResult:
    r"""Converged KLL Fredholm-iteration result for the slab.

    Attributes
    ----------
    b : float
        Slab half-thickness (mfp). Same as the F_N
        :attr:`SlabFNResult.a_critical_mfp`.
    c : float
        Mean number of secondaries per collision.
    u0 : float
        Magnitude :math:`u_0 = |\nu_0|` of the imaginary Case discrete
        eigenvalue (:math:`c > 1`).
    nu_grid : np.ndarray, shape (n_nodes,)
        Quadrature nodes :math:`\nu \in (0, 1)` used for the
        Fredholm-equation solution.
    nu_weights : np.ndarray, shape (n_nodes,)
        Quadrature weights matching ``nu_grid``.
    A_nu : np.ndarray, shape (n_nodes,)
        Converged continuum coefficient :math:`A(\nu)` evaluated on
        ``nu_grid``.
    iterations : int
        Number of KLL iterations performed.
    converged : bool
        Whether the iteration converged within ``tol``.
    sup_diff_history : list[float]
        Sup-norm of :math:`A_n - A_{n-1}` per iteration.

    Notes
    -----
    The result is suitable for evaluating
    :func:`slab_scalar_flux_kll` and
    :func:`slab_angular_flux_from_scalar`. The :class:`SlabFNResult`
    +  :class:`KLLSlabFluxResult` together form the rich-machinery
    pair: F_N gives :math:`r_c`, KLL gives interior :math:`\phi(z)`.
    """

    b: float
    c: float
    u0: float
    nu_grid: np.ndarray
    nu_weights: np.ndarray
    A_nu: np.ndarray
    iterations: int
    converged: bool
    sup_diff_history: list[float]


def _build_KLL_slab_kernel(
    nu_grid: np.ndarray,
    b: float,
    c: float,
    u0: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    r"""Build the KLL slab Fredholm kernel matrix and source vector.

    KLL Eq. 6 (rearranged):

    .. math::

       \mu A(\mu) = \mu A_1(\mu) - \int_0^1 K_A(\mu, \nu)\,\nu A(\nu)\,
       d\nu, \quad \mu \in (0, 1)

    where :math:`K_A(\mu, \nu) = \gamma(\mu)/(\mu + \nu) \cdot
    e^{-2b/\nu}\,X(-\nu)` (kernel) and :math:`\mu A_1(\mu)` is the
    inhomogeneous source (KLL):

    .. math::

       \mu A_1(\mu) = -\frac{2 u_0\,\gamma(\mu)}{\mu^2 + u_0^2}\,
       \big\{u_0 [X_{10} \cos(b/u_0) + X_{20} \sin(b/u_0)]
       - \mu [X_{10} \sin(b/u_0) - X_{20} \cos(b/u_0)]\big\} .

    Returns ``(K_matrix, mu_A1_source, integration_weights)`` so the
    iteration step is :math:`\mu A_n = \mu A_1 - K \cdot (\nu A_{n-1})
    \cdot W`.
    """
    n_nodes = len(nu_grid)
    X10, X20 = _X10_X20(c, u0)

    # Build the source vector μ A_1(μ) at each grid point.
    mu_A1 = np.empty(n_nodes)
    for i, mu in enumerate(nu_grid):
        gamma_mu = _gamma_kll(mu, c, u0)
        prefactor = -2.0 * u0 * gamma_mu / (mu * mu + u0 * u0)
        cos_part = X10 * math.cos(b / u0) + X20 * math.sin(b / u0)
        sin_part = X10 * math.sin(b / u0) - X20 * math.cos(b / u0)
        mu_A1[i] = prefactor * (u0 * cos_part - mu * sin_part)

    # Build kernel K(μ, ν) = γ(μ) / (μ + ν) · X(-ν) · exp(-2b/ν).
    # Cache γ(μ_i), X(-ν_j), exp(-2b/ν_j) once.
    gamma_at_grid = np.array([_gamma_kll(mu, c, u0) for mu in nu_grid])
    X_neg_at_grid = np.array(
        [_X_negative_real_axis(nu, c, u0) for nu in nu_grid]
    )
    exp_factor_at_grid = np.exp(-2.0 * b / nu_grid)

    K_matrix = np.empty((n_nodes, n_nodes))
    for i, mu in enumerate(nu_grid):
        for j, nu in enumerate(nu_grid):
            K_matrix[i, j] = (
                gamma_at_grid[i]
                / (mu + nu)
                * X_neg_at_grid[j]
                * exp_factor_at_grid[j]
            )

    return K_matrix, mu_A1, X_neg_at_grid


def solve_kll_slab_continuum_coefficient(
    b: float,
    c: float,
    *,
    n_nodes: int = 64,
    tol: float = 1e-12,
    max_iter: int = 200,
) -> KLLSlabFluxResult:
    r"""Solve the KLL Fredholm integral equation Eq. 6 for
    :math:`A(\nu)` via successive iteration.

    Parameters
    ----------
    b : float
        Slab half-thickness in mean-free-paths. Typically the F_N
        result :attr:`SlabFNResult.a_critical_mfp`.
    c : float
        Mean number of secondaries per collision (:math:`> 1`).
    n_nodes : int, default 64
        Number of Gauss-Legendre quadrature nodes on :math:`(0, 1)`.
        Convergence is geometric in :math:`n_{\rm nodes}`; 64 nodes
        gives ~1e-10 accuracy on flux ratios for typical :math:`b`.
    tol : float, default 1e-12
        Convergence tolerance on sup-norm of :math:`A_n - A_{n-1}`.
    max_iter : int, default 200
        Maximum number of iterations.

    Returns
    -------
    KLLSlabFluxResult

    Notes
    -----
    The iteration converges geometrically with rate
    :math:`\frac{c}{2} e^{-2b}` (KLL p. 96 — see "convergence criteria"
    paragraph, with their factor). For :math:`c < 2` (i.e. all
    physical multiplying media), the rate is :math:`< 1` and
    convergence is guaranteed.
    """
    if c <= 1.0:
        raise ValueError(f"KLL slab requires c > 1, got c = {c}")
    if b <= 0.0:
        raise ValueError(f"b must be positive, got {b}")

    u0, is_imag = case_nu0(c)
    if not is_imag:
        raise RuntimeError(
            "Internal: expected supercritical c>1 with imaginary nu_0"
        )

    # Gauss-Legendre nodes on (0, 1).
    nu_grid_unit, w_unit = np.polynomial.legendre.leggauss(n_nodes)
    # Map [-1, 1] → (0, 1).
    nu_grid = 0.5 * (nu_grid_unit + 1.0)
    nu_weights = 0.5 * w_unit

    K_matrix, mu_A1, _ = _build_KLL_slab_kernel(nu_grid, b, c, u0)

    # Iterate: μ A_n(μ) = μ A_1(μ) - K · (ν A_{n-1}) · W.
    # Where ν A is the unknown (not A directly — we store μA on the grid).
    # Initial guess: μ A_0 = μ A_1 (the inhomogeneous source).
    mu_A_curr = mu_A1.copy()
    sup_diff_history: list[float] = []

    for iteration in range(max_iter):
        # ν A_{n-1}(ν): we have μ A on the grid, ν is the same grid.
        # Compute (K · (νA · W))_i at row i.
        nu_A_prev = mu_A_curr  # because grid is shared
        # The Fredholm equation has ∫ K(μ, ν) ν A(ν) dν → discretize as
        # K @ (ν A · W) where W is the diagonal of weights.
        weighted_nu_A = nu_A_prev * nu_weights
        rhs_correction = K_matrix @ weighted_nu_A
        mu_A_next = mu_A1 - rhs_correction

        sup_diff = float(np.max(np.abs(mu_A_next - mu_A_curr)))
        sup_diff_history.append(sup_diff)
        mu_A_curr = mu_A_next

        if sup_diff < tol:
            converged = True
            break
    else:
        converged = False
        iteration = max_iter - 1

    # A_ν = (μ A) / μ where μ = ν on the grid.
    A_nu = mu_A_curr / nu_grid

    return KLLSlabFluxResult(
        b=b,
        c=c,
        u0=u0,
        nu_grid=nu_grid,
        nu_weights=nu_weights,
        A_nu=A_nu,
        iterations=iteration + 1,
        converged=converged,
        sup_diff_history=sup_diff_history,
    )


# ─────────────────────────────────────────────────────────────────────
# Scalar-flux evaluation via KLL Eq. 7
# ─────────────────────────────────────────────────────────────────────


def slab_scalar_flux_kll(
    kll_result: KLLSlabFluxResult,
    z: float | np.ndarray,
) -> float | np.ndarray:
    r"""Evaluate KLL Eq. 7 for the bare-critical slab scalar flux.

    .. math::

       \phi(z) = a\,\Big[\cos(z/u_0)
       + \int_0^1 A(\nu)\,e^{-b/\nu}\,\cosh(z/\nu)\,d\nu\Big]

    The normalisation :math:`a = 1` is used here. The published
    benchmark (KLL Tables III, IV; Sood Table 14) tabulates
    :math:`\phi(z)/\phi(0)` which is :math:`a`-independent — use
    :func:`slab_scalar_flux_ratio` for that.

    Parameters
    ----------
    kll_result : :class:`KLLSlabFluxResult`
        Converged Fredholm result.
    z : float or array
        Position(s) in :math:`[-b, b]` (mfp). Symmetric: :math:`\phi(-z)
        = \phi(z)`.

    Returns
    -------
    Same shape as ``z``.
    """
    z_arr = np.asarray(z, dtype=float)
    u0 = kll_result.u0
    b = kll_result.b
    nu_grid = kll_result.nu_grid
    nu_weights = kll_result.nu_weights
    A_nu = kll_result.A_nu

    # cos(z/u0) term.
    cos_term = np.cos(z_arr / u0)

    # ∫ A(ν) exp(-b/ν) cosh(z/ν) dν via Gauss-Legendre quadrature.
    # Use the numerically-stable form
    #   exp(-b/ν) cosh(z/ν) = (exp(-(b-z)/ν) + exp(-(b+z)/ν)) / 2
    # which avoids overflow at small ν (cosh(z/ν) → ∞ but
    # exp(-b/ν) → 0 faster — the product is bounded by 1/2 since
    # |z| ≤ b).
    if z_arr.ndim == 0:
        z_val = float(z_arr)
        if abs(z_val) > b + 1e-12:
            raise ValueError(
                f"|z| = {abs(z_val)} exceeds slab half-thickness b = {b}"
            )
        # exp_factor[ν] = exp(-(b-|z|)/ν) + exp(-(b+|z|)/ν), both ≤ 1.
        # Note z² is even so cosh(z/ν) = cosh(|z|/ν).
        z_abs = abs(z_val)
        exp1 = np.exp(-(b - z_abs) / nu_grid)
        exp2 = np.exp(-(b + z_abs) / nu_grid)
        kernel = 0.5 * (exp1 + exp2)
        integrand = A_nu * kernel
        integral = float(np.sum(integrand * nu_weights))
        phi = cos_term + integral
        if np.isscalar(z):
            return float(phi)
        return phi
    else:
        # Vectorised for array z.
        z_abs = np.abs(z_arr)
        if np.any(z_abs > b + 1e-12):
            raise ValueError(
                f"max |z| = {np.max(z_abs)} exceeds slab half-thickness "
                f"b = {b}"
            )
        z_col = z_abs[..., None]
        exp1 = np.exp(-(b - z_col) / nu_grid[None, ...])
        exp2 = np.exp(-(b + z_col) / nu_grid[None, ...])
        kernel = 0.5 * (exp1 + exp2)  # (..., n_nodes)
        integrand = A_nu[None, ...] * kernel
        integral = np.sum(integrand * nu_weights[None, ...], axis=-1)
        return cos_term + integral


def slab_scalar_flux_ratio(
    kll_result: KLLSlabFluxResult,
    z_over_b: float | np.ndarray,
) -> float | np.ndarray:
    r"""Compute :math:`\phi(z) / \phi(0)` at fractional position
    :math:`z/b`.

    This is the quantity tabulated in KLL Table III + Sood Table 14.

    Parameters
    ----------
    kll_result : :class:`KLLSlabFluxResult`
    z_over_b : float or array
        Fractional position(s) :math:`z/b \in [-1, 1]`.

    Returns
    -------
    Same shape as ``z_over_b``.
    """
    z_arr = np.asarray(z_over_b, dtype=float) * kll_result.b
    phi_z = slab_scalar_flux_kll(kll_result, z_arr)
    phi_0 = slab_scalar_flux_kll(kll_result, 0.0)
    return phi_z / phi_0


# ─────────────────────────────────────────────────────────────────────
# Angular-flux reconstruction via characteristic integration
# ─────────────────────────────────────────────────────────────────────


def slab_angular_flux_from_scalar(
    kll_result: KLLSlabFluxResult,
    z: float | np.ndarray,
    mu: float | np.ndarray,
    *,
    n_quad: int = 64,
) -> np.ndarray | float:
    r"""Reconstruct interior angular flux :math:`\psi(z, \mu)` from the
    converged scalar flux :math:`\phi(z)` via characteristic
    integration of the BTE.

    For the bare-critical slab on :math:`[-b, b]` with vacuum BCs and
    isotropic scattering, the angular flux satisfies

    .. math::

       \mu \frac{\partial \psi}{\partial z} + \psi = \frac{c}{2} \phi(z)

    with :math:`\psi(-b, \mu) = 0` for :math:`\mu > 0` and
    :math:`\psi(b, \mu) = 0` for :math:`\mu < 0`. Integration along
    the characteristic gives the closed-form

    .. math::

       \psi(z, \mu > 0) &= \frac{c}{2 \mu}\int_{-b}^{z}\phi(z')\,
       e^{-(z-z')/\mu}\,dz', \\
       \psi(z, \mu < 0) &= \frac{c}{2|\mu|}\int_{z}^{b}\phi(z')\,
       e^{-(z'-z)/|\mu|}\,dz' .

    The integration is performed via Gauss-Legendre quadrature with
    ``n_quad`` nodes mapped to the characteristic-segment interval.

    Parameters
    ----------
    kll_result : :class:`KLLSlabFluxResult`
        Converged scalar-flux result (used to evaluate :math:`\phi(z')`
        on the integration nodes).
    z : float or array
        Position(s) in :math:`(-b, b)` mfp.
    mu : float or array
        Direction cosine(s) in :math:`[-1, 1] \setminus \{0\}`. The
        :math:`\mu = 0` direction is structural-singular for the slab
        (perpendicular streaming); reject it.
    n_quad : int, default 64
        Number of Gauss-Legendre nodes for the characteristic-segment
        integral. 64 nodes gives ~1e-12 accuracy on the integral itself.

    Returns
    -------
    Array of shape ``np.broadcast(z, mu).shape`` with the angular flux.

    Notes
    -----
    The normalisation here is :math:`a = 1` (matching
    :func:`slab_scalar_flux_kll`). Use :math:`\phi = \int_{-1}^{1}
    \psi\,d\mu` as a closure-consistency check (see
    :func:`slab_scalar_flux_from_angular_quadrature` for that
    self-consistency).
    """
    z_arr = np.asarray(z, dtype=float)
    mu_arr = np.asarray(mu, dtype=float)
    if np.any(mu_arr == 0.0):
        raise ValueError(
            "slab_angular_flux_from_scalar: mu = 0 is structural-"
            "singular (perpendicular streaming); choose μ ≠ 0."
        )

    # Broadcast z and mu to common shape.
    z_b, mu_b = np.broadcast_arrays(z_arr, mu_arr)
    out = np.empty_like(z_b)

    b = kll_result.b
    c = kll_result.c

    # Gauss-Legendre quadrature on [-1, 1].
    nodes_unit, weights_unit = np.polynomial.legendre.leggauss(n_quad)

    flat_z = z_b.ravel()
    flat_mu = mu_b.ravel()
    flat_out = out.ravel()

    for i, (zi, mui) in enumerate(zip(flat_z, flat_mu, strict=True)):
        if mui > 0:
            # ψ(z, μ > 0) = (c/(2μ)) ∫_(-b)^z φ(z') exp(-(z-z')/μ) dz'.
            lo, hi = -b, zi
            if hi <= lo:
                # z < -b → outside the slab (or at the boundary).
                flat_out[i] = 0.0
                continue
            mid = 0.5 * (lo + hi)
            half = 0.5 * (hi - lo)
            z_prime = mid + half * nodes_unit
            phi_at_z_prime = slab_scalar_flux_kll(kll_result, z_prime)
            integrand = phi_at_z_prime * np.exp(-(zi - z_prime) / mui)
            integral = half * np.sum(integrand * weights_unit)
            flat_out[i] = (c / (2.0 * mui)) * integral
        else:
            # μ < 0.
            mui_abs = -mui
            lo, hi = zi, b
            if hi <= lo:
                flat_out[i] = 0.0
                continue
            mid = 0.5 * (lo + hi)
            half = 0.5 * (hi - lo)
            z_prime = mid + half * nodes_unit
            phi_at_z_prime = slab_scalar_flux_kll(kll_result, z_prime)
            integrand = phi_at_z_prime * np.exp(-(z_prime - zi) / mui_abs)
            integral = half * np.sum(integrand * weights_unit)
            flat_out[i] = (c / (2.0 * mui_abs)) * integral

    if z_b.ndim == 0 and mu_b.ndim == 0:
        return float(flat_out[0])
    return out


def slab_scalar_flux_from_angular_quadrature(
    kll_result: KLLSlabFluxResult,
    z: float | np.ndarray,
    *,
    n_mu: int = 96,
    n_quad: int = 64,
) -> float | np.ndarray:
    r"""Compute :math:`\phi(z) = \int_{-1}^{1} \psi(z, \mu)\,d\mu` by
    direct numerical quadrature.

    This is the closure-consistency test: if
    :func:`slab_angular_flux_from_scalar` is correct, then
    :math:`\phi^{\rm angular-quad}(z)` should equal :math:`\phi^{\rm KLL}(z)`
    to within the quadrature accuracy.

    Parameters
    ----------
    kll_result : :class:`KLLSlabFluxResult`
    z : float or array
        Position(s) in :math:`(-b, b)`.
    n_mu : int, default 96
        Number of Gauss-Legendre :math:`\mu` quadrature nodes.
    n_quad : int, default 64
        Number of nodes for the inner characteristic integral
        (passed through to :func:`slab_angular_flux_from_scalar`).

    Returns
    -------
    Same shape as ``z``.
    """
    z_arr = np.asarray(z, dtype=float)
    nodes_mu, w_mu = np.polynomial.legendre.leggauss(n_mu)
    # Skip μ = 0 (none in Gauss-Legendre nodes for n_mu even).
    if z_arr.ndim == 0:
        psi_at_z = slab_angular_flux_from_scalar(
            kll_result, z_arr, nodes_mu, n_quad=n_quad
        )
        return float(np.sum(psi_at_z * w_mu))
    out = np.empty_like(z_arr)
    for i, zi in enumerate(z_arr.ravel()):
        psi_at_z = slab_angular_flux_from_scalar(
            kll_result, zi, nodes_mu, n_quad=n_quad
        )
        out.ravel()[i] = np.sum(psi_at_z * w_mu)
    return out


# ─────────────────────────────────────────────────────────────────────
# Surface angular flux from F_N coefficients (closed-form)
# ─────────────────────────────────────────────────────────────────────


def slab_surface_angular_flux_fn(
    coefficients_a: np.ndarray,
    mu: float | np.ndarray,
) -> float | np.ndarray:
    r"""Evaluate the F_N surface angular flux representation.

    Per Siewert-Benoist Part I Eq. 45:

    .. math::

       \psi(-b, -\mu) = \sum_{\alpha=0}^{N} a_\alpha \mu^\alpha,
       \qquad \mu > 0 .

    By symmetry of the bare-critical slab eigenmode, :math:`\psi(b,
    \mu) = \psi(-b, -\mu)` for :math:`\mu > 0`. So this entry point
    gives the **outgoing** angular flux at both surfaces.

    Parameters
    ----------
    coefficients_a : np.ndarray, shape (N+1,)
        F_N expansion coefficients :math:`a_\alpha` from
        :attr:`SlabFNResult.coefficients_a`. Complex-valued in general
        but real for the bare-critical problem (imaginary parts at the
        machine-noise floor).
    mu : float or array
        Outgoing direction cosine in :math:`(0, 1]`.

    Returns
    -------
    Real angular flux value(s). Imaginary part of the F_N coefficients
    is dropped (it is at the noise floor).
    """
    mu_arr = np.asarray(mu, dtype=float)
    coeffs_real = np.asarray(coefficients_a).real
    # Σ_α a_α μ^α: numpy.polynomial.polynomial.polyval uses ascending
    # power order (a_0 + a_1 μ + a_2 μ² + ...), which matches F_N.
    return np.polynomial.polynomial.polyval(mu_arr, coeffs_real)
