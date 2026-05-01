r"""Special-function kernels for continuous reference solutions.

This module provides the two kernel families used across the
flat-source CP and pointwise Peierls references:

- **Exponential integrals** :math:`E_n(x)`: double-precision
  :func:`e3` / :func:`e3_vec` thin wrappers over
  :func:`scipy.special.expn`, plus arbitrary-precision :func:`e_n`
  / :func:`e_n_mp` via :func:`mpmath.expint`.

- **Bickley-Naylor functions** :math:`\mathrm{Ki}_n(x)`: arbitrary-
  precision :func:`ki_n` / :func:`ki_n_mp` via :func:`mpmath.quad` on
  the canonical A&S 11.2 definition. Double-precision fast paths live
  with their consumers (e.g. :func:`~.cp_geometry._ki3_mp` builds a
  Chebyshev interpolant of :math:`e^{\tau}\,\mathrm{Ki}_3(\tau)` from
  :func:`ki_n_mp` at module load).

Both families satisfy the differential identities

.. math::

   E_n'(x) = -E_{n-1}(x), \qquad
   \mathrm{Ki}_n'(x) = -\mathrm{Ki}_{n-1}(x),

exposed as :func:`e_n_derivative` / :func:`ki_n_derivative` and
verified term-by-term in ``tests/derivations/test_kernels.py`` (L0).

.. note::

   The legacy :class:`BickleyTables` tabulation (20 000-point
   :math:`\mathrm{Ki}_3` lookup with ~:math:`10^{-3}` absolute
   accuracy and a naming discrepancy under the A&S convention —
   Issue #94) was retired in Phase B.4. Every former consumer now
   uses the mpmath-built Chebyshev interpolant in
   :mod:`~.cp_geometry`.
"""

from __future__ import annotations

import math
from typing import Callable

import mpmath
import numpy as np
from scipy.special import expn


# ═══════════════════════════════════════════════════════════════════════
# Chord geometry — concentric annular / spherical regions
# ═══════════════════════════════════════════════════════════════════════

def chord_half_lengths(radii: np.ndarray, y_pts: np.ndarray) -> np.ndarray:
    r"""Half-chord lengths :math:`\ell_k(y)` for concentric annuli / shells.

    For a cylinder or sphere partitioned into ``N`` regions with outer
    radii ``radii = [r_1, r_2, ..., r_N]`` (``0 < r_1 < ... < r_N = R``)
    and a chord of impact parameter :math:`y`, the **half-chord length**
    through region :math:`k` (inner radius :math:`r_{k-1}`, outer radius
    :math:`r_k`; :math:`r_0 \equiv 0`) is

    .. math::

        \ell_k(y) = \begin{cases}
            0, & y \ge r_k,\\
            \sqrt{r_k^{2}-y^{2}}, & r_{k-1} \le y < r_k,\\
            \sqrt{r_k^{2}-y^{2}}-\sqrt{r_{k-1}^{2}-y^{2}}, & y < r_{k-1}.
        \end{cases}

    This is the primitive consumed by the CP cylinder / sphere
    derivations (:mod:`~orpheus.derivations.continuous.flat_source_cp.cylinder`,
    :mod:`~orpheus.derivations.continuous.flat_source_cp.sphere`) and by the Phase-4 Peierls
    cylinder reference (:mod:`~orpheus.derivations.continuous.peierls.cylinder`).
    Tested at L0 in ``tests/derivations/test_kernels.py``.

    Parameters
    ----------
    radii : np.ndarray, shape (N,)
        Outer radii of the ``N`` concentric regions, strictly increasing.
    y_pts : np.ndarray, shape (n_y,)
        Impact parameters at which the chord lengths are evaluated.

    Returns
    -------
    np.ndarray, shape (N, n_y)
        ``chords[k, m] = ℓ_k(y_pts[m])``.
    """
    N = len(radii)
    chords = np.zeros((N, len(y_pts)))
    r_inner = np.zeros(N)
    r_inner[1:] = radii[:-1]
    y2 = y_pts ** 2

    for k in range(N):
        r_out, r_in = radii[k], r_inner[k]
        outer = np.sqrt(np.maximum(r_out ** 2 - y2, 0.0))
        if r_in > 0:
            inner = np.sqrt(np.maximum(r_in ** 2 - y2, 0.0))
            mask = y_pts < r_in
            chords[k, mask] = outer[mask] - inner[mask]
        mask_p = (y_pts >= r_in) & (y_pts < r_out)
        chords[k, mask_p] = outer[mask_p]

    return chords


# ═══════════════════════════════════════════════════════════════════════
# E₃ kernel (slab geometry)
# ═══════════════════════════════════════════════════════════════════════

def e3(x: float) -> float:
    """Third-order exponential integral E₃(x)."""
    return float(expn(3, max(x, 0.0)))


def e3_vec(x: np.ndarray) -> np.ndarray:
    """Vectorised E₃."""
    return expn(3, np.maximum(x, 0.0))


# ═══════════════════════════════════════════════════════════════════════
# High-precision tier — mpmath-backed E_n and Ki_n
# ═══════════════════════════════════════════════════════════════════════
#
# These wrap :func:`mpmath.expint` (for :math:`E_n`) and
# :func:`mpmath.quad` (for :math:`\mathrm{Ki}_n`) at a caller-chosen
# working precision. They are the primitives consumed by the
# Phase-4 Peierls Nyström references; the Phase-0 kernel-identity
# tests exercise them directly.
#
# The Bickley function is defined as
#
#     Ki_n(x) = ∫_0^{π/2} cos^{n-1}(θ) exp(-x/cos θ) dθ
#
# (Bickley & Naylor 1935, Abramowitz & Stegun 11.2). The integrand
# is smooth on (0, π/2) but develops an exp(-∞) at θ = π/2 when x > 0
# and a mild endpoint behaviour at θ = 0 when x → 0. ``mpmath.quad``
# handles both correctly via adaptive Gauss–Legendre.


def e_n(n: int, x: float, precision_digits: int = 50) -> float:
    r"""Exponential integral :math:`E_n(x)` at arbitrary precision.

    Wraps :func:`mpmath.expint` and returns a Python ``float``
    rounded from the mpmath result at the requested working precision.
    For double-precision callers the default of 50 digits is
    massively over-specified but keeps all consumers consistent;
    Phase-4 Peierls callers request 50 explicitly and store the
    mpmath result unrounded via :func:`e_n_mp`.

    Parameters
    ----------
    n : int
        Order of the exponential integral. ``n >= 1``.
    x : float
        Argument. ``x >= 0``.
    precision_digits : int
        mpmath working precision in decimal digits.
    """
    with mpmath.workdps(precision_digits):
        return float(mpmath.expint(n, x))


def e_n_mp(n: int, x, precision_digits: int = 50):
    """Exponential integral :math:`E_n(x)` returning an mpmath ``mpf``.

    Use this when the caller will chain further mpmath arithmetic
    and must preserve precision.
    """
    with mpmath.workdps(precision_digits):
        return mpmath.expint(n, mpmath.mpf(x))


def e_n_derivative(n: int, x: float, precision_digits: int = 50) -> float:
    r"""Analytical derivative :math:`E_n'(x) = -E_{n-1}(x)`.

    For ``n == 1`` we fall back to the closed form
    :math:`E_1'(x) = -e^{-x}/x` since :math:`E_0(x) = e^{-x}/x`
    diverges at :math:`x = 0`.

    This function exists primarily as the RHS of the kernel-identity
    tests — the identity is *defined* at the mpmath layer, and the
    test verifies that the numerical derivative of ``e_n(n, x)``
    matches ``e_n_derivative(n, x)``.
    """
    if n < 1:
        raise ValueError(f"E_n derivative requires n >= 1, got {n}")
    if n == 1:
        with mpmath.workdps(precision_digits):
            if x == 0:
                raise ValueError("E_1'(0) diverges")
            return float(-mpmath.exp(-mpmath.mpf(x)) / mpmath.mpf(x))
    with mpmath.workdps(precision_digits):
        return float(-mpmath.expint(n - 1, x))


def _ki_integrand(n: int, x):
    r"""Integrand for :math:`\mathrm{Ki}_n(x)` on :math:`[0, \infty)`.

    Uses the substitution :math:`u = \tan\theta`, :math:`du = \sec^{2}\theta\,d\theta`,
    :math:`\cos\theta = 1/\sqrt{1+u^{2}}`, to map
    :math:`\theta \in [0, \pi/2]` to :math:`u \in [0, \infty)` and remove
    the essential singularity at :math:`\theta = \pi/2`. The transformed
    integrand

    .. math::

        \mathrm{Ki}_n(x) \;=\; \int_{0}^{\infty}
            (1+u^{2})^{-(n+1)/2}\,
            \exp\!\bigl(-x\,\sqrt{1+u^{2}}\bigr)\,du

    is smooth and bounded on :math:`[0, \infty)` for all :math:`x \ge 0`
    and decays like :math:`u^{-(n+1)}\,e^{-x u}` at infinity, which
    :func:`mpmath.quad` handles uniformly.
    """
    n_mp = mpmath.mpf(n)
    x_mp = mpmath.mpf(x)
    exponent = -(n_mp + 1)

    def integrand(u):
        s = mpmath.sqrt(1 + u * u)
        return s ** exponent * mpmath.exp(-x_mp * s)

    return integrand


def ki_n(n: int, x: float, precision_digits: int = 50) -> float:
    r"""Bickley–Naylor function :math:`\mathrm{Ki}_n(x)` at arbitrary precision.

    Canonical definition:

    .. math::

        \mathrm{Ki}_n(x) \;=\; \int_{0}^{\pi/2}
            \cos^{n-1}\theta\;
            \exp\!\left(-\tfrac{x}{\cos\theta}\right)\,d\theta.

    Internally evaluated via the substitution :math:`u = \tan\theta`,
    which removes the essential singularity at :math:`\theta = \pi/2`
    and gives a smooth, exponentially decaying integrand on
    :math:`[0, \infty)`. :func:`mpmath.quad` converges to the target
    precision for every :math:`(n, x)` with :math:`n \ge 1` and
    :math:`x \ge 0`.

    Cross-checked against the Wallis closed form
    :math:`\mathrm{Ki}_n(0)` (see :func:`ki_n_at_zero`) in the
    L0 kernel-identity tests.

    Parameters
    ----------
    n : int
        Order. Must be :math:`\ge 1`.
    x : float
        Argument. :math:`x \ge 0`.
    precision_digits : int
        mpmath working precision in decimal digits.
    """
    if n < 1:
        raise ValueError(f"Ki_n requires n >= 1, got {n}")
    if x < 0:
        raise ValueError(f"Ki_n requires x >= 0, got {x}")
    with mpmath.workdps(precision_digits):
        val = mpmath.quad(_ki_integrand(n, x), [0, mpmath.inf])
        return float(val)


def ki_n_mp(n: int, x, precision_digits: int = 50):
    """Bickley–Naylor :math:`\\mathrm{Ki}_n` returning an mpmath ``mpf``."""
    if n < 1:
        raise ValueError(f"Ki_n requires n >= 1, got {n}")
    with mpmath.workdps(precision_digits):
        return mpmath.quad(_ki_integrand(n, float(x)), [0, mpmath.inf])


def ki_n_float(n: int, x: float) -> float:
    r"""Fast float-precision Bickley–Naylor :math:`\mathrm{Ki}_n(x)`.

    Uses :func:`scipy.integrate.quad` on the tanh-substituted integrand

    .. math::

       \mathrm{Ki}_n(x) \;=\; \int_{0}^{\infty}
           (1+u^{2})^{-(n+1)/2}\,\exp\!\bigl(-x\,\sqrt{1+u^{2}}\bigr)\,
           \mathrm d u

    which is smooth and exponentially decaying. ~100× faster than
    :func:`ki_n_mp` at ``precision_digits=20``, with ~1e-14 relative
    accuracy — sufficient for every float-valued Peierls kernel
    assembly (the output is always cast to float before use in the
    Nyström matrix). Designed for use inside hot loops of
    :func:`~orpheus.derivations.continuous.peierls.geometry.build_volume_kernel`.
    """
    if n < 1:
        raise ValueError(f"Ki_n requires n >= 1, got {n}")
    if x < 0:
        raise ValueError(f"Ki_n requires x >= 0, got {x}")
    # scipy.integrate.quad handles the exponentially-decaying integrand
    # on [0, ∞) via its internal tanh-sinh / Gauss-Kronrod fallback.
    # Integrand is s^{-(n+1)} · exp(-x·s) with s = √(1+u²) — identical
    # to the mpmath _ki_integrand so the two routines agree to 1e-14.
    from scipy.integrate import quad
    exponent = -(n + 1)

    def integrand(u: float) -> float:
        s = math.sqrt(1.0 + u * u)
        return s ** exponent * math.exp(-x * s)

    val, _err = quad(integrand, 0.0, np.inf, epsabs=1e-15, epsrel=1e-13)
    return float(val)


def ki_n_derivative(n: int, x: float, precision_digits: int = 50) -> float:
    r"""Analytical derivative :math:`\mathrm{Ki}_n'(x) = -\mathrm{Ki}_{n-1}(x)`.

    For :math:`n = 1` the limit :math:`\mathrm{Ki}_0(x) = K_0(x)`
    (modified Bessel function) is used via ``mpmath.besselk``.
    Reference: Abramowitz & Stegun 11.2.11.
    """
    if n < 1:
        raise ValueError(f"Ki_n' requires n >= 1, got {n}")
    if n == 1:
        # Ki_1'(x) = -Ki_0(x) = -K_0(x) (modified Bessel)
        with mpmath.workdps(precision_digits):
            if x == 0:
                raise ValueError("Ki_1'(0) diverges (K_0 has log singularity)")
            return float(-mpmath.besselk(0, x))
    return -ki_n(n - 1, x, precision_digits)


# ── Special values (used by kernel-identity tests) ───────────────────

def ki_n_at_zero(n: int) -> float:
    r"""Exact closed form :math:`\mathrm{Ki}_n(0)` for ``n >= 1``.

    Evaluates :math:`\int_0^{\pi/2} \cos^{n-1}\theta\,d\theta` in
    closed form via the Wallis integral. Used by kernel-identity
    tests (and as a sanity initialiser for the Bickley table).

    .. math::

        \mathrm{Ki}_1(0) = \tfrac{\pi}{2},\quad
        \mathrm{Ki}_2(0) = 1,\quad
        \mathrm{Ki}_3(0) = \tfrac{\pi}{4},\quad
        \mathrm{Ki}_4(0) = \tfrac{2}{3}.
    """
    with mpmath.workdps(50):
        # Wallis: ∫_0^{π/2} cos^k θ dθ for k = n-1
        k = n - 1
        if k == 0:
            return float(mpmath.pi / 2)
        if k % 2 == 0:  # even
            # (π/2) * (k-1)!! / k!!
            num = mpmath.mpf(1)
            den = mpmath.mpf(1)
            for i in range(1, k, 2):
                num *= i
            for i in range(2, k + 1, 2):
                den *= i
            return float(mpmath.pi / 2 * num / den)
        # odd k
        num = mpmath.mpf(1)
        den = mpmath.mpf(1)
        for i in range(2, k, 2):
            num *= i
        for i in range(1, k + 1, 2):
            den *= i
        return float(num / den)


def e_n_at_zero(n: int) -> float:
    r"""Exact closed form :math:`E_n(0) = 1/(n-1)` for :math:`n > 1`.

    :math:`E_1(0)` diverges logarithmically — excluded here.
    """
    if n <= 1:
        raise ValueError(f"E_n(0) requires n > 1 (E_1 diverges), got {n}")
    return 1.0 / (n - 1)


# ═══════════════════════════════════════════════════════════════════════
# Shifted Legendre polynomials (Gelbard DP_N half-range basis)
# ═══════════════════════════════════════════════════════════════════════

def _shifted_legendre_eval(n: int, mu: np.ndarray) -> np.ndarray:
    r"""Evaluate the shifted Legendre polynomial :math:`\tilde P_n(\mu)
    = P_n(2\mu - 1)` on :math:`\mu \in [0, 1]`.

    This is the Gelbard half-range / double-P\ :sub:`N` basis used by
    the rank-N Marshak white-boundary closure for the pointwise Peierls
    Nyström operator (see :ref:`theory-peierls-unified` §8 and Sanchez
    & McCormick 1982 §III.F.1 Eqs. 165–169). The orthogonality relation

    .. math::

       \int_0^1 \tilde P_n(\mu)\,\tilde P_m(\mu)\,\mathrm d\mu
         \;=\; \frac{\delta_{nm}}{2n+1}

    is the reason this specific rescaling is preferred over the Case–
    Zweifel half-range polynomials (which carry a material-dependent
    weight) for CP / Peierls work.

    Evaluated via the three-term Bonnet recurrence

    .. math::

       (k+1)\,P_{k+1}(x) = (2k+1)\,x\,P_k(x) - k\,P_{k-1}(x),
       \qquad P_0 = 1,\; P_1 = x,

    on the affine-transformed argument :math:`x = 2\mu - 1 \in [-1, 1]`.
    This is stable for all :math:`n \lesssim 100` in double precision,
    far beyond the :math:`n \le 8` range of practical rank-N closures.

    Parameters
    ----------
    n : int
        Polynomial order, :math:`n \ge 0`.
    mu : np.ndarray
        Evaluation points, assumed to lie in :math:`[0, 1]`.

    Returns
    -------
    np.ndarray
        :math:`\tilde P_n(\mu)` at each input point, same shape as ``mu``.
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    mu = np.asarray(mu, dtype=float)
    x = 2.0 * mu - 1.0
    if n == 0:
        return np.ones_like(x)
    if n == 1:
        return x
    P_km1 = np.ones_like(x)
    P_k = x
    for k in range(1, n):
        P_kp1 = ((2 * k + 1) * x * P_k - k * P_km1) / (k + 1)
        P_km1, P_k = P_k, P_kp1
    return P_k
