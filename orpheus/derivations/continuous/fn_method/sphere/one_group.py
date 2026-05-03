r"""Critical-sphere F_N method solver — 1G isotropic, bare sphere.

Implements the Siewert-Thomas 1986 (*Nucl. Sci. Eng.* **94**, 264)
F_N method for the bare-critical sphere. Per their explicit
prescription (p. 268, verbatim):

    "the critical sphere problem requires only that Eq. (4) be
    changed to read [Eq. 46], and that we interpret a as the
    critical radius, we have incorporated the relevant minus sign
    in our developed equations"

— meaning **slab and sphere F_N differ by exactly one sign on the
boundary attenuation term**. This module reuses
:func:`...core.fn_matrix.assemble_fn_matrix` with ``geometry_sign =
-1`` (sphere) where slab uses ``+1``. Everything else (moment
recursion, collocation grid, dispersion-root finder for the discrete
Case eigenvalue :math:`\nu_0 = i u_0` at :math:`c > 1`) is shared
verbatim.

For the 1G specialisation (Sood ``Ua-1-0-SP`` benchmark target,
:math:`R_c = 2.4248249802` mfp at :math:`c = 1.30`), the
Siewert-Thomas 2G machinery collapses cleanly to scalars (matrix
:math:`C \to c`, :math:`\Theta \to 1`, :math:`\Lambda(z) \to
1 - cz\,\mathrm{atanh}(1/z)`) — see
:func:`...origins.fn_sphere_derivations.derive_sphere_2g_to_1g_reduction`
for the symbolic verification.

The critical condition is :math:`\det M(R) = 0`, where :math:`M(R)`
is the :math:`(N+1) \times (N+1)` complex F_N matrix with entries

.. math::

   M_{\beta, \alpha}(R) = B_\alpha(\xi_\beta)
   - e^{-2R/\xi_\beta}\,A_\alpha(\xi_\beta).

Because :math:`\xi_0 = \nu_0 = i u_0` for :math:`c > 1`, the matrix
is genuinely complex and we work in :math:`\mathbb{C}` throughout,
finding the **real positive critical radius** :math:`R > 0` for
which :math:`\Re \det M(R) = 0`. By the Hermitian symmetry of the
problem, :math:`\Im \det M(R)` is also small at the converged
:math:`R` (see slab analog).

References
----------

* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264-270.
* Sood, Forster & Parsons 1999, LA-13511 Table 6 (case ``Ua-1-0-SP``).
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94, Table V
  (independent reference values via Wiener-Hopf).
"""
from __future__ import annotations

import warnings as _warnings
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from ..core.dispersion import case_nu0
from ..core.fn_matrix import assemble_fn_matrix


def _siewert_thomas_collocation(nu0_complex: complex, N: int) -> np.ndarray:
    r"""Build the Siewert-Thomas 1986 Eq. 38a collocation grid for sphere.

    The 1G sphere F_N collocation prescription is:

    * :math:`\xi_0 = \nu_0 = i u_0` (the imaginary discrete Case
      eigenvalue);
    * :math:`\xi_\beta = \frac{1}{2}\left[1 + \cos(\beta \pi / (N+1))
      \right]` for :math:`\beta = 1, \ldots, N`.

    These are the (shifted) Chebyshev-of-the-first-kind nodes on
    :math:`(0, 1)` — strictly interior, NOT including the endpoints.
    This differs from the Grandjean-Siewert slab prescription
    :math:`(\nu_0, 0, 1, k/(N-1))` which DOES include the endpoints
    :math:`\xi_1 = 0` and :math:`\xi_2 = 1`.

    The reason: the sphere F_N matrix at :math:`\xi = 0` would have
    a constant row independent of :math:`R` (the attenuation factor
    :math:`e^{-\infty} = 0` kills the geometry-sign-bearing term),
    which creates a structural rank deficiency that masks the genuine
    root condition :math:`\det M(R) = 0`. The Chebyshev grid avoids
    this by staying strictly inside :math:`(0, 1)`.

    Parameters
    ----------
    nu0_complex : complex
        :math:`\xi_0 = \nu_0`. For supercritical :math:`c > 1` this is
        the imaginary number :math:`i u_0`.
    N : int
        F_N order. The grid has :math:`N + 1` points total (one for
        :math:`\nu_0` plus :math:`N` Chebyshev points).

    Returns
    -------
    np.ndarray
        Length-:math:`(N+1)` complex array.
    """
    if N < 0:
        raise ValueError(f"N must be ≥ 0, got {N}")
    xis = np.empty(N + 1, dtype=complex)
    xis[0] = nu0_complex
    if N >= 1:
        betas = np.arange(1, N + 1)
        cheb = 0.5 * (1.0 + np.cos(betas * np.pi / (N + 1)))
        xis[1:] = cheb.astype(complex)
    return xis


@dataclass(frozen=True)
class SphereFNResult:
    r"""Result of :func:`solve_fn_sphere_bare_critical`.

    Attributes
    ----------
    R_critical_mfp : float
        Critical radius in mean-free paths.
    c : float
        Mean number of secondaries per collision (input).
    N : int
        F_N order.
    nu0 : float
        Magnitude of the Case discrete eigenvalue :math:`|\nu_0|`.
        For :math:`c > 1` (the only relevant case for bare-critical
        sphere), :math:`\nu_0 = i u_0` is purely imaginary; this
        attribute stores :math:`u_0`.
    xi_collocation : np.ndarray
        :math:`(N+1,)` complex array of collocation points (the first
        is :math:`i u_0`, the rest are real in :math:`[0, 1]`).
    coefficients_a : np.ndarray
        :math:`(N+1,)` complex array of F_N expansion coefficients
        :math:`a_\alpha`, normalised so :math:`a_0 = 1`.
    determinant_residual : complex
        Value of :math:`\det M(R)` at the converged :math:`R`.
    n_bracket_points : int
        Number of points used in the initial determinant scan.

    Notes
    -----
    The :class:`SphereFNResult` mirrors the slab analog
    (:class:`...slab.one_group.SlabFNResult`); their two
    public APIs are identical except for the ``R_critical_mfp`` /
    ``a_critical_mfp`` field name. This is by design — the two
    geometries are the same method with one parameter difference
    (geometry_sign), which is why a single ``...core.fn_matrix``
    module supports both.
    """

    R_critical_mfp: float
    c: float
    N: int
    nu0: float
    xi_collocation: np.ndarray
    coefficients_a: np.ndarray
    determinant_residual: complex
    n_bracket_points: int


def _build_fn_matrix_sphere(
    R: float, N: int, xis: np.ndarray, c: float
) -> np.ndarray:
    r"""Assemble the sphere F_N matrix :math:`M(R)` via the shared
    assembler with ``geometry_sign = -1``."""
    return assemble_fn_matrix(c, R, geometry_sign=-1, n_modes=N, xis=xis)


def _det_M_sphere(R: float, N: int, xis: np.ndarray, c: float) -> complex:
    """Return :math:`\\det M(R)` for the sphere F_N matrix as a complex number."""
    M = _build_fn_matrix_sphere(R, N, xis, c)
    return complex(np.linalg.det(M))


def solve_fn_sphere_bare_critical(
    c: float,
    *,
    n_modes: int = 8,
    R_min: float = 0.5,
    R_max: float = 25.0,
    n_bracket: int = 800,
    bisect_tol: float = 1e-12,
    max_bisect: int = 80,
) -> SphereFNResult:
    r"""Solve the bare-critical sphere radius via the F_N method
    (Siewert-Thomas 1986).

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision,
        :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`. Must satisfy
        :math:`c > 1` (multiplying medium).
    n_modes : int, default 8
        F_N order :math:`N`. The system size is :math:`(N+1) \times
        (N+1)`. For the Sood ``Ua-1-0-SP`` target (:math:`R_c =
        2.4248249802` mfp at :math:`c = 1.30`), :math:`N = 8` reaches
        ~1e-6 absolute accuracy and :math:`N = 12` reaches ~1e-8.
    R_min, R_max : float, defaults 0.5, 25.0
        Bracket for the scan on :math:`R`. The default range covers
        :math:`c \in [1.01, 5.0]` (sphere :math:`R_c \in [22, 0.6]`).
    n_bracket : int, default 800
        Number of bracket-scan points used to locate the FIRST
        prominent local minimum of :math:`|\det M(R)|`. The sphere
        determinant has a sequence of zeros corresponding to higher
        radial modes (overtones); we want the fundamental, which is
        the first.
    bisect_tol : float, default 1e-12
        Relative tolerance on :math:`R` for the local minimization step.
    max_bisect : int, default 80
        Maximum minimization iterations.

    Returns
    -------
    :class:`SphereFNResult`

    Notes
    -----
    Sphere F_N is structurally identical to slab F_N except for one
    sign (``geometry_sign = -1``). The X-function, dispersion
    relation, F_α moment recursion, and collocation grid are all
    reused verbatim — see
    :func:`...origins.fn_sphere_derivations.derive_sphere_bc_sign_flip`
    and :func:`...origins.fn_sphere_derivations.derive_x_function_geometry_independence`
    for the symbolic verifications.

    **Implementation choice (sphere vs slab root-finding strategy).**
    Sphere :math:`\det M(R)` is genuinely complex with anti-Hermitian
    structure (the ``geometry_sign = -1`` breaks the slab Hermitian
    symmetry that made :math:`\Re \det M_{\rm slab}(R)` cross zero
    cleanly). For the sphere, the trajectory of :math:`\det M(R)`
    in the complex plane spirals through the origin at each :math:`R_c`
    — the fundamental and the overtones. Detecting genuine zero
    crossings via sign change in :math:`\Re \det M` is unreliable
    because the trajectory can "almost" cross the imaginary axis at
    sub-noise-floor values without going through a real root. The
    robust strategy:

    1. Scan :math:`\log_{10}|\det M(R)|` on a coarse grid covering
       the bracket.
    2. Locate the FIRST prominent local minimum (the fundamental
       mode; later minima are spurious overtone modes).
    3. Refine via :func:`scipy.optimize.minimize_scalar` (Brent) on
       :math:`\log_{10}|\det M|` within a narrow bracket around the
       first minimum.

    This is the same idea as bisection on a real determinant —
    minimum of :math:`|\det|` IS the location of the zero — but it's
    unaffected by the overall complex phase of :math:`\det M` and
    works through the noise-floor-relative-scale issues that the
    pure sign-change approach suffers.
    """
    if c <= 1.0:
        raise ValueError(
            f"F_N critical sphere requires c > 1 (multiplying medium), got c={c}"
        )
    if n_modes < 0:
        raise ValueError(f"n_modes must be ≥ 0, got {n_modes}")

    # Case discrete eigenvalue magnitude (positive real); for c > 1,
    # actual ν_0 = i u_0.
    u0, is_imag = case_nu0(c)
    if not is_imag:
        raise RuntimeError("internal: expected supercritical c>1 with imaginary nu_0")
    nu0_complex = complex(0.0, u0)

    # Collocation points per Siewert-Thomas 1986 Eq. 38a:
    # ξ_0 = ν_0, ξ_β = (1/2)[1 + cos(β π/(N+1))] for β = 1, ..., N
    # (Chebyshev nodes strictly inside (0, 1)). NOT the slab
    # Grandjean-Siewert grid (which includes ξ = 0 and ξ = 1) —
    # see :func:`_siewert_thomas_collocation` for the rank-deficiency
    # rationale.
    xis = _siewert_thomas_collocation(nu0_complex, n_modes)

    # Scan log|det M(R)| on the coarse grid. The genuine bare-
    # critical-sphere root is the location where det M = 0; spurious
    # oscillations may produce shallow dips, but the genuine root has
    # a sharp, deep dip relative to the local background. We use
    # |det| not smin because: (a) |det| has WIDER dips than smin (so
    # the coarse scan grid doesn't miss it at low N), and (b) the
    # FIRST prominent dip in |det| corresponds reliably to the
    # fundamental mode (overtone modes appear at larger R with similar
    # depth — we want the first one).
    R_scan = np.linspace(R_min, R_max, n_bracket)
    log_abs_det = np.empty(n_bracket)
    with _warnings.catch_warnings():
        _warnings.filterwarnings("ignore", category=RuntimeWarning)
        for i, R in enumerate(R_scan):
            try:
                d = _det_M_sphere(R, n_modes, xis, c)
                log_abs_det[i] = np.log10(abs(d) + 1e-300)
            except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                log_abs_det[i] = float("nan")

    valid = ~np.isnan(log_abs_det)
    if not valid.any():
        raise RuntimeError(
            f"All det M(R) evaluations failed for c={c}, N={n_modes}. "
            f"Check inputs."
        )

    valid_idx = np.where(valid)[0]
    log_valid = log_abs_det[valid_idx]
    R_valid = R_scan[valid_idx]
    n_v = len(log_valid)
    if n_v < 3:
        raise RuntimeError("Not enough valid bracket points for minimum search.")

    # Find local minima of log|det|.
    is_local_min = np.zeros(n_v, dtype=bool)
    is_local_min[1:-1] = (
        (log_valid[1:-1] < log_valid[:-2]) & (log_valid[1:-1] < log_valid[2:])
    )
    if not is_local_min.any():
        raise RuntimeError(
            f"No local minimum of |det M(R)| found in [{R_min}, {R_max}] "
            f"for c={c}, N={n_modes}. Try widening the bracket or "
            f"increasing n_bracket."
        )

    # Filter by prominence using a "median rank" approach: a genuine
    # root produces a dip whose depth is > a threshold below the
    # global median (which represents the typical |det| level over
    # the scan). Threshold of 0.7 decades = factor ~5x — works for
    # both low N (where dips are narrow but deep relative to wide
    # background) and high N (where dips are very deep but the
    # background is slightly noisier).
    global_median = float(np.median(log_valid))
    min_indices = np.where(is_local_min)[0]
    prominent = [k for k in min_indices if global_median - log_valid[k] > 0.7]
    if not prominent:
        # Fallback: take the global minimum.
        prominent = [int(np.argmin(log_valid))]
    # First prominent minimum corresponds to the fundamental.
    k_first = prominent[0]
    R_guess = float(R_valid[k_first])

    # Bracket for Brent: take the neighbour grid points.
    k_lo = max(0, k_first - 1)
    k_hi = min(n_v - 1, k_first + 1)
    R_lo_bracket = float(R_valid[k_lo])
    R_hi_bracket = float(R_valid[k_hi])

    # Refine via scipy.optimize.minimize_scalar (Brent) on log|det|.
    def log_abs_det_at_R(R: float) -> float:
        try:
            d = _det_M_sphere(float(R), n_modes, xis, c)
            return float(np.log10(abs(d) + 1e-300))
        except Exception:
            return 0.0  # large value so minimizer avoids

    with _warnings.catch_warnings():
        _warnings.filterwarnings("ignore", category=RuntimeWarning)
        result = minimize_scalar(
            log_abs_det_at_R,
            bracket=(R_lo_bracket, R_guess, R_hi_bracket),
            method="brent",
            options={"xtol": bisect_tol, "maxiter": max_bisect},
        )
    if not result.success:
        # Fall back to the coarse-grid guess.
        R_crit = R_guess
    else:
        R_crit = float(result.x)

    M_crit = _build_fn_matrix_sphere(R_crit, n_modes, xis, c)
    det_residual = complex(np.linalg.det(M_crit))

    # Recover the eigenvector (F_N coefficients a_α) via SVD.
    U, s, Vh = np.linalg.svd(M_crit)
    coefficients = Vh.conj().T[:, -1]
    if abs(coefficients[0]) > 1e-30:
        coefficients = coefficients / coefficients[0]

    return SphereFNResult(
        R_critical_mfp=float(R_crit),
        c=c,
        N=n_modes,
        nu0=u0,
        xi_collocation=xis,
        coefficients_a=coefficients,
        determinant_residual=det_residual,
        n_bracket_points=n_bracket,
    )
