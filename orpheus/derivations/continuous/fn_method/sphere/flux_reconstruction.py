r"""Interior flux reconstruction for the bare-critical sphere via the
Kaper-Lindeman-Leaf (KLL) 1974 Fredholm-iteration recipe.

Sphere analog of :mod:`..slab.flux_reconstruction`. The structural
specification is identical (KLL Eqs. 13-15 mirror Eqs. 5-7 with sign
flips and :math:`\sinh \leftrightarrow \cosh` due to the
anti-symmetric :math:`r\Phi(r)` reformulation), so we reuse the
medium-only quantities (:math:`X(z)`, :math:`\theta(\mu)`,
:math:`\gamma(\mu)`) imported from the slab module.

Public API:

* :func:`solve_kll_sphere_continuum_coefficient` — solves the sphere
  Fredholm integral equation Eq. 14 for :math:`B(\nu)` via successive
  iteration.
* :func:`sphere_scalar_flux_kll` — evaluates KLL Eq. 15 for the
  bare-critical sphere scalar flux :math:`\phi(r)`.
* :func:`sphere_scalar_flux_ratio` — :math:`\phi(r)/\phi(0)` as
  tabulated in KLL Table VII.
* :func:`sphere_angular_flux_from_scalar` — interior angular flux
  :math:`\psi(r, \mu)` via chord-length characteristic integration.

References
----------

* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94, Eqs. 13-15.
* Mitsis 1963, ANL-6787.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

from ..core.dispersion import case_nu0
from ..slab.flux_reconstruction import (
    _gamma_kll,
    _X_negative_real_axis,
    _X10_X20,
)


@dataclass(frozen=True)
class KLLSphereFluxResult:
    r"""Converged KLL Fredholm-iteration result for the sphere.

    Mirrors :class:`...slab.flux_reconstruction.KLLSlabFluxResult`,
    with ``R`` replacing ``b`` (radius vs half-thickness) and ``B_nu``
    the sphere continuum coefficient.

    Attributes
    ----------
    R : float
        Sphere radius (mfp). Same as the F_N
        :attr:`SphereFNResult.R_critical_mfp`.
    c : float
        Mean number of secondaries per collision.
    u0 : float
        Magnitude :math:`u_0 = |\nu_0|` of the imaginary Case discrete
        eigenvalue.
    nu_grid : np.ndarray
        Quadrature nodes :math:`\nu \in (0, 1)`.
    nu_weights : np.ndarray
        Quadrature weights.
    B_nu : np.ndarray
        Converged continuum coefficient :math:`B(\nu)` on
        ``nu_grid``.
    iterations : int
        Number of iterations used.
    converged : bool
        Whether the iteration converged within ``tol``.
    sup_diff_history : list[float]
        Sup-norm of :math:`B_n - B_{n-1}` per iteration.
    """

    R: float
    c: float
    u0: float
    nu_grid: np.ndarray
    nu_weights: np.ndarray
    B_nu: np.ndarray
    iterations: int
    converged: bool
    sup_diff_history: list[float]


def _build_KLL_sphere_kernel(
    nu_grid: np.ndarray,
    R: float,
    c: float,
    u0: float,
) -> tuple[np.ndarray, np.ndarray]:
    r"""Build the KLL sphere Fredholm kernel matrix and source vector.

    KLL Eq. 14:

    .. math::

       \mu B(\mu) = \mu B_1(\mu) + \int_0^1 K_B(\mu, \nu)\,\nu B(\nu)\,
       d\nu, \quad \mu \in (0, 1)

    where (KLL):

    .. math::

       \mu B_1(\mu) = \frac{2 |u_0|^2 \gamma(\mu)}{\mu^2 + u_0^2}\,
       \big\{|u_0| [X_{10} \sin(R/u_0) - X_{20} \cos(R/u_0)]
       + \mu [X_{10} \cos(R/u_0) + X_{20} \sin(R/u_0)]\big\}

    and :math:`K_B(\mu, \nu) = \gamma(\mu)/(\mu + \nu) \cdot
    e^{-2R/\nu} X(-\nu)`. The kernel has the SAME form as the slab —
    only the sign on the inhomogeneous term changes between the slab
    and sphere Fredholm equations (and the sin/cos pattern
    rearranges, reflecting :math:`\sinh \leftrightarrow \cosh` in the
    flux formulae).

    Returns ``(K_matrix, mu_B1_source)``.
    """
    n_nodes = len(nu_grid)
    X10, X20 = _X10_X20(c, u0)

    # Build μ B_1(μ) at each grid point.
    mu_B1 = np.empty(n_nodes)
    for i, mu in enumerate(nu_grid):
        gamma_mu = _gamma_kll(mu, c, u0)
        prefactor = 2.0 * u0 * u0 * gamma_mu / (mu * mu + u0 * u0)
        sin_minus_part = X10 * math.sin(R / u0) - X20 * math.cos(R / u0)
        cos_plus_part = X10 * math.cos(R / u0) + X20 * math.sin(R / u0)
        mu_B1[i] = prefactor * (u0 * sin_minus_part + mu * cos_plus_part)

    # Build kernel K(μ, ν) = γ(μ)/(μ + ν) · X(-ν) · exp(-2R/ν) — same
    # form as slab.
    gamma_at_grid = np.array([_gamma_kll(mu, c, u0) for mu in nu_grid])
    X_neg_at_grid = np.array(
        [_X_negative_real_axis(nu, c, u0) for nu in nu_grid]
    )
    exp_factor_at_grid = np.exp(-2.0 * R / nu_grid)

    K_matrix = np.empty((n_nodes, n_nodes))
    for i, mu in enumerate(nu_grid):
        for j, nu in enumerate(nu_grid):
            K_matrix[i, j] = (
                gamma_at_grid[i]
                / (mu + nu)
                * X_neg_at_grid[j]
                * exp_factor_at_grid[j]
            )

    return K_matrix, mu_B1


def solve_kll_sphere_continuum_coefficient(
    R: float,
    c: float,
    *,
    n_nodes: int = 64,
    tol: float = 1e-12,
    max_iter: int = 200,
) -> KLLSphereFluxResult:
    r"""Solve the KLL Fredholm integral equation Eq. 14 for
    :math:`B(\nu)` via successive iteration.

    Parameters
    ----------
    R : float
        Sphere radius in mean-free-paths.
    c : float
        Mean number of secondaries per collision (:math:`> 1`).
    n_nodes : int, default 64
        Number of Gauss-Legendre quadrature nodes on :math:`(0, 1)`.
    tol : float, default 1e-12
        Convergence tolerance.
    max_iter : int, default 200
        Maximum iterations.

    Returns
    -------
    KLLSphereFluxResult

    Notes
    -----
    The iteration converges with rate :math:`(c/2)\,e^{-2R}`, so for
    physical :math:`R_c \gtrsim 0.5` mfp and :math:`c < 2` the
    convergence is geometric and rapid (~5 iterations for 1e-10).
    """
    if c <= 1.0:
        raise ValueError(f"KLL sphere requires c > 1, got c = {c}")
    if R <= 0.0:
        raise ValueError(f"R must be positive, got {R}")

    u0, is_imag = case_nu0(c)
    if not is_imag:
        raise RuntimeError(
            "Internal: expected supercritical c>1 with imaginary nu_0"
        )

    nu_grid_unit, w_unit = np.polynomial.legendre.leggauss(n_nodes)
    nu_grid = 0.5 * (nu_grid_unit + 1.0)
    nu_weights = 0.5 * w_unit

    K_matrix, mu_B1 = _build_KLL_sphere_kernel(nu_grid, R, c, u0)

    mu_B_curr = mu_B1.copy()
    sup_diff_history: list[float] = []

    for iteration in range(max_iter):
        # Sphere Fredholm has + sign on integral (vs slab -).
        nu_B_prev = mu_B_curr
        weighted_nu_B = nu_B_prev * nu_weights
        rhs_correction = K_matrix @ weighted_nu_B
        mu_B_next = mu_B1 + rhs_correction

        sup_diff = float(np.max(np.abs(mu_B_next - mu_B_curr)))
        sup_diff_history.append(sup_diff)
        mu_B_curr = mu_B_next

        if sup_diff < tol:
            converged = True
            break
    else:
        converged = False
        iteration = max_iter - 1

    B_nu = mu_B_curr / nu_grid

    return KLLSphereFluxResult(
        R=R,
        c=c,
        u0=u0,
        nu_grid=nu_grid,
        nu_weights=nu_weights,
        B_nu=B_nu,
        iterations=iteration + 1,
        converged=converged,
        sup_diff_history=sup_diff_history,
    )


# ─────────────────────────────────────────────────────────────────────
# Scalar-flux evaluation via KLL Eq. 15
# ─────────────────────────────────────────────────────────────────────


def sphere_scalar_flux_kll(
    kll_result: KLLSphereFluxResult,
    r: float | np.ndarray,
) -> float | np.ndarray:
    r"""Evaluate KLL Eq. 15 for the bare-critical sphere scalar flux.

    .. math::

       \phi(r) = \frac{2 b}{r}\Big[\frac{u_0\,\sin(r/u_0)}{1}
       - \int_0^1 B(\nu)\,e^{-R/\nu}\,\sinh(r/\nu)\,d\nu\Big]

    where :math:`b` is the (arbitrary) normalisation constant. We use
    :math:`b = 1`. The :math:`r \to 0` limit is well-defined:

    .. math::

       \phi(0) = 2 \Big[1 - \int_0^1 \frac{B(\nu)}{\nu}\,e^{-R/\nu}\,
       d\nu\Big]

    using :math:`\lim_{r \to 0} \sin(r/u_0)/(r/u_0) = 1` and
    :math:`\lim_{r \to 0} \sinh(r/\nu)/r = 1/\nu`. (Note: KLL Eq. 15
    uses the prefactor :math:`2 b /r` and the full inner expression
    has units that combine to a finite scalar flux at :math:`r = 0`;
    the ``b = 1`` normalisation amounts to a multiplicative factor
    that cancels in :math:`\phi(r)/\phi(0)`.)

    Numerical-stability note: like the slab analog, we rewrite
    :math:`e^{-R/\nu} \sinh(r/\nu) = (e^{-(R-r)/\nu} -
    e^{-(R+r)/\nu})/2` to avoid overflow in :math:`\sinh` at small
    :math:`\nu`.

    Parameters
    ----------
    kll_result : :class:`KLLSphereFluxResult`
    r : float or array
        Radial position(s) in :math:`[0, R]` mfp.

    Returns
    -------
    Same shape as ``r``.
    """
    r_arr = np.asarray(r, dtype=float)
    u0 = kll_result.u0
    R = kll_result.R
    nu_grid = kll_result.nu_grid
    nu_weights = kll_result.nu_weights
    B_nu = kll_result.B_nu

    if np.any(r_arr < -1e-12) or np.any(r_arr > R + 1e-12):
        raise ValueError(f"r must be in [0, R={R}], got {r}")

    # The form is (2/r) [u0 sin(r/u0) - ∫ B(ν) e^(-R/ν) sinh(r/ν) dν].
    # Note: u0 sin(r/u0) = sinc-like; at r → 0, → r (since sin(x)/x → 1
    # times x → 0 ... actually u0 · sin(r/u0) = u0 · (r/u0) = r at r=0,
    # so the prefactor (2/r) gives finite 2 at r = 0).
    # And ∫ B(ν) e^(-R/ν) sinh(r/ν) dν → ∫ B(ν) e^(-R/ν) (r/ν) dν, so
    # (2/r) · that → 2 ∫ B(ν) e^(-R/ν) / ν dν, finite.

    if r_arr.ndim == 0:
        r_val = float(r_arr)
        if abs(r_val) < 1e-15:
            # r → 0 limit.
            sin_term = 1.0  # sin(r/u0)/(r/u0) → 1, times u0 ... actually
            # u0 · sin(r/u0) ≈ r so (2/r) · u0 · sin(r/u0) → 2.
            # Use Taylor: u0 sin(r/u0) ≈ r - r^3/(6 u0^2) + ...
            # Use the limit form directly.
            sinh_integrand_limit = B_nu * np.exp(-R / nu_grid) / nu_grid
            integral_lim = float(np.sum(sinh_integrand_limit * nu_weights))
            return 2.0 * (1.0 - integral_lim)

        # exp(-R/ν) sinh(r/ν) = (exp(-(R-r)/ν) - exp(-(R+r)/ν)) / 2.
        if r_val > R:
            raise ValueError(f"r = {r_val} exceeds R = {R}")
        exp1 = np.exp(-(R - r_val) / nu_grid)
        exp2 = np.exp(-(R + r_val) / nu_grid)
        kernel = 0.5 * (exp1 - exp2)
        integrand = B_nu * kernel
        integral = float(np.sum(integrand * nu_weights))

        sin_term_val = u0 * math.sin(r_val / u0)
        phi = (2.0 / r_val) * (sin_term_val - integral)
        return float(phi) if np.isscalar(r) else phi

    # Vectorised.
    out = np.empty_like(r_arr)
    flat = out.ravel()
    flat_r = r_arr.ravel()
    for i, r_val in enumerate(flat_r):
        flat[i] = sphere_scalar_flux_kll(kll_result, float(r_val))
    return out


def sphere_scalar_flux_ratio(
    kll_result: KLLSphereFluxResult,
    r_over_R: float | np.ndarray,
) -> float | np.ndarray:
    r"""Compute :math:`\phi(r) / \phi(0)` at fractional position
    :math:`r/R`.

    This is the quantity tabulated in KLL Table VII + Sood Table 8.

    Parameters
    ----------
    kll_result : :class:`KLLSphereFluxResult`
    r_over_R : float or array
        Fractional position(s) in :math:`[0, 1]`.

    Returns
    -------
    Same shape as ``r_over_R``.
    """
    r_arr = np.asarray(r_over_R, dtype=float) * kll_result.R
    phi_r = sphere_scalar_flux_kll(kll_result, r_arr)
    phi_0 = sphere_scalar_flux_kll(kll_result, 0.0)
    return phi_r / phi_0


# ─────────────────────────────────────────────────────────────────────
# Angular-flux reconstruction via chord-length characteristic
# ─────────────────────────────────────────────────────────────────────


def sphere_angular_flux_from_scalar(
    kll_result: KLLSphereFluxResult,
    r: float | np.ndarray,
    mu: float | np.ndarray,
    *,
    n_quad: int = 64,
) -> np.ndarray | float:
    r"""Reconstruct interior angular flux :math:`\psi(r, \mu)` from the
    converged scalar flux :math:`\phi(r)` via chord-length
    characteristic integration.

    For the bare-critical sphere of radius :math:`R` with vacuum BC,
    the angular flux at :math:`(r, \mu)` (where :math:`\mu = \cos
    \theta` is the cosine of the angle between the direction of motion
    and the outward radial vector) follows from the BTE along the
    chord:

    .. math::

       \psi(r, \mu) = \frac{c}{2}\int_0^{s_{\rm in}}
       \phi\!\left(\sqrt{r^2 - 2 r s' \mu + s'^2}\right)\,e^{-s'}\,ds'

    where :math:`s_{\rm in}(r, \mu)` is the chord length back to the
    surface (vacuum inflow). For :math:`\mu \in [-1, 1]`,

    .. math::

       s_{\rm in}(r, \mu) = r\mu + \sqrt{R^2 - r^2(1 - \mu^2)} .

    (Sign convention: :math:`\mu > 0` means motion outward; the chord
    enters the sphere at distance :math:`s_{\rm in}` BEHIND the current
    point, where the chord crosses the surface and vacuum inflow is
    zero.)

    Parameters
    ----------
    kll_result : :class:`KLLSphereFluxResult`
    r : float or array
        Radial position(s) in :math:`[0, R]`.
    mu : float or array
        Direction cosine(s) in :math:`[-1, 1]`. The :math:`\mu = 0`
        case is the perpendicular-streaming chord — handled correctly
        by the formula (chord length is :math:`\sqrt{R^2 - r^2}`).
    n_quad : int, default 64
        Gauss-Legendre nodes for the characteristic integral.

    Returns
    -------
    Array shape ``np.broadcast(r, mu).shape``.

    Notes
    -----
    For :math:`r = R, \mu < 0` (incoming particles at the surface),
    :math:`s_{\rm in} = 0` and :math:`\psi = 0` (vacuum inflow).
    For :math:`r = R, \mu > 0` (outgoing), :math:`s_{\rm in} = 2 R \mu`
    (the diameter chord projected on the radial direction).
    """
    r_arr = np.asarray(r, dtype=float)
    mu_arr = np.asarray(mu, dtype=float)

    R = kll_result.R
    c = kll_result.c

    if np.any(r_arr < -1e-12) or np.any(r_arr > R + 1e-12):
        raise ValueError(f"r must be in [0, R={R}], got max {np.max(r_arr)}")

    r_b, mu_b = np.broadcast_arrays(r_arr, mu_arr)
    out = np.empty_like(r_b, dtype=float)

    nodes_unit, weights_unit = np.polynomial.legendre.leggauss(n_quad)

    flat_r = r_b.ravel()
    flat_mu = mu_b.ravel()
    flat_out = out.ravel()

    for i, (ri, mui) in enumerate(zip(flat_r, flat_mu, strict=True)):
        # Chord length back to surface.
        # s_in(r, μ) = r μ + sqrt(R² - r²(1 - μ²)) per characteristic
        # geometry. For r = R, μ = -1: s_in = -R + R = 0 (vacuum
        # incoming). For r = 0: s_in = R (radius).
        discriminant = R * R - ri * ri * (1.0 - mui * mui)
        if discriminant < 0:
            # Should not happen for r ∈ [0, R].
            flat_out[i] = 0.0
            continue
        s_in = ri * mui + math.sqrt(discriminant)
        if s_in <= 0.0:
            # Particle just left the boundary with no path inside.
            flat_out[i] = 0.0
            continue

        # Map [-1, 1] → [0, s_in].
        mid = 0.5 * s_in
        half = 0.5 * s_in
        s_prime = mid + half * nodes_unit
        # Position along chord: r(s') = sqrt(r² - 2 r s' μ + s'²).
        # When the chord starts at the inner point at s'=0 and goes
        # back, the distance from origin is sqrt(r² - 2 r s' cos(θ_in) +
        # s'²) where cos(θ_in) is the angle to the inner-pointing
        # direction. Since μ is outward-pointing, the inward s'
        # corresponds to going opposite to μ, hence the (-μ) sign:
        r_prime = np.sqrt(np.maximum(
            ri * ri - 2.0 * ri * s_prime * mui + s_prime * s_prime,
            0.0,
        ))
        phi_at_r_prime = sphere_scalar_flux_kll(kll_result, r_prime)
        integrand = phi_at_r_prime * np.exp(-s_prime)
        integral = half * np.sum(integrand * weights_unit)
        flat_out[i] = (c / 2.0) * integral

    if r_b.ndim == 0 and mu_b.ndim == 0:
        return float(flat_out[0])
    return out


def sphere_scalar_flux_from_angular_quadrature(
    kll_result: KLLSphereFluxResult,
    r: float | np.ndarray,
    *,
    n_mu: int = 96,
    n_quad: int = 64,
) -> float | np.ndarray:
    r"""Compute :math:`\phi(r) = \int_{-1}^{1} \psi(r, \mu)\,d\mu` by
    direct numerical quadrature over :math:`\mu`.

    Closure-consistency test: should equal
    :func:`sphere_scalar_flux_kll` to within quadrature accuracy.
    """
    r_arr = np.asarray(r, dtype=float)
    nodes_mu, w_mu = np.polynomial.legendre.leggauss(n_mu)
    if r_arr.ndim == 0:
        psi_at_r = sphere_angular_flux_from_scalar(
            kll_result, float(r_arr), nodes_mu, n_quad=n_quad
        )
        return float(np.sum(psi_at_r * w_mu))
    out = np.empty_like(r_arr)
    for i, ri in enumerate(r_arr.ravel()):
        psi_at_r = sphere_angular_flux_from_scalar(
            kll_result, float(ri), nodes_mu, n_quad=n_quad
        )
        out.ravel()[i] = np.sum(psi_at_r * w_mu)
    return out


# ─────────────────────────────────────────────────────────────────────
# Surface angular flux from F_N coefficients
# ─────────────────────────────────────────────────────────────────────


def sphere_surface_angular_flux_fn(
    coefficients_a: np.ndarray,
    mu: float | np.ndarray,
) -> float | np.ndarray:
    r"""Evaluate the F_N surface angular flux for the sphere.

    Per Siewert-Thomas 1986 Eq. 25 (specialised to 1G):

    .. math::

       \psi(R, \mu) = \sum_{\alpha=0}^{N} a_\alpha\,\mu^\alpha,
       \qquad \mu \in [0, 1] (\text{outgoing})

    By the sphere BC :math:`\Psi(-R, \mu) = -\Psi(R, \mu)` (Siewert-
    Thomas Eq. 46), the inner and outer surface angular fluxes are
    related by sign — but since we only have a positive-r physical
    surface (the sphere geometry only has one boundary), the
    representation above is the complete surface angular flux.

    Parameters
    ----------
    coefficients_a : np.ndarray, shape (N+1,)
        F_N expansion coefficients :math:`a_\alpha` from
        :attr:`SphereFNResult.coefficients_a`.
    mu : float or array
        Outgoing direction cosine in :math:`[0, 1]`.

    Returns
    -------
    Real angular flux value(s).
    """
    mu_arr = np.asarray(mu, dtype=float)
    coeffs_real = np.asarray(coefficients_a).real
    return np.polynomial.polynomial.polyval(mu_arr, coeffs_real)
