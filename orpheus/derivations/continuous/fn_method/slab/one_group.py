r"""Critical-slab F_N method solver — 1G isotropic, bare slab.

Implements the Siewert-Benoist 1979 (Part I, Eq. 49) / Grandjean-Siewert
1979 (Part II, Eq. 25) F_N collocation system for the critical
half-thickness of a bare homogeneous slab with isotropic scattering:

.. math::

   \sum_{\alpha=0}^{N} a_\alpha\,
   \big[B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta}\,A_\alpha(\xi_\beta)\big]
   = 0, \qquad \beta = 0, 1, \ldots, N

with collocation points :math:`\xi_0 = \nu_0,\, \xi_1 = 0,\, \xi_2 = 1`,
and the remaining :math:`N - 2` interior points equally spaced in
:math:`(0, 1)` (per the Grandjean-Siewert prescription).

For multiplying media (:math:`c > 1`), :math:`\nu_0` is purely
imaginary, :math:`\nu_0 = i u_0` with :math:`u_0` the real positive
root of :math:`1 - c u\,\mathrm{atan}(1/u) = 0`. The first
collocation row is therefore complex; we work in complex linear
algebra throughout and find the **real positive critical
half-thickness** :math:`a > 0` for which the complex determinant
:math:`\det M(a)` has a real-part zero crossing. Symmetry of the
problem ensures the imaginary part is also small at the converged
:math:`a`.

The :math:`\xi = 0` collocation point requires special handling:
:math:`\exp(-2a/0) = 0` (since :math:`a > 0`), and
:math:`\xi B_{\alpha-1}(\xi) \to 0` as :math:`\xi \to 0`, so:

* :math:`B_0(0) = 2/c - 1` (the :math:`\xi \log(1+1/\xi) \to 0` limit).
* :math:`B_\alpha(0) = -1/(\alpha+1)` for :math:`\alpha \ge 1`.
* :math:`A_\alpha(0)` does not contribute (the :math:`\exp` factor is 0).

References
----------

* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156-160. Part I.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161-168. Part II,
  Section III, Eqs. 25-26 + Table XI ("The Critical Thickness").
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94 (Table I).
* Sood, Forster & Parsons 1999, LANL LA-13511 (Table 4).
"""
from __future__ import annotations

import cmath
from dataclasses import dataclass

import numpy as np

from ..core.dispersion import case_nu0
from ..core.fn_matrix import assemble_fn_matrix
from ..core.moments import collocation_points


@dataclass(frozen=True)
class SlabFNResult:
    r"""Result of :func:`solve_fn_slab_bare_critical`.

    Attributes
    ----------
    a_critical_mfp : float
        Critical half-thickness in mean free paths.
    c : float
        Mean number of secondaries per collision.
    N : int
        F_N order.
    nu0 : float
        Magnitude of the Case discrete eigenvalue
        :math:`|\nu_0|`. Real for :math:`c < 1`; magnitude of
        :math:`i u_0` for :math:`c > 1`.
    xi_collocation : np.ndarray
        :math:`(N+1,)` complex array of collocation points
        (real for :math:`\xi \in [0, 1]`, imaginary for the
        :math:`\xi_0 = \nu_0` entry when :math:`c > 1`).
    coefficients_a : np.ndarray
        :math:`(N+1,)` complex array of F_N expansion coefficients
        :math:`a_\alpha`, normalised so that :math:`a_0 = 1`. These
        define :math:`\psi(-a, -\mu) = \sum_\alpha a_\alpha \mu^\alpha`
        and (by symmetry of the critical problem)
        :math:`\psi(a, \mu) = \psi(-a, -\mu) = \sum a_\alpha \mu^\alpha`.
    determinant_residual : complex
        Value of :math:`\det M(a)` at the converged :math:`a`. Should
        be small (complex magnitude :math:`\lesssim 10^{-10}`).
    n_bracket_points : int
        Number of points used in the initial determinant scan.
    converged : bool
        Whether the bisection met ``bisect_tol`` — **required, no
        default** (#340).

        ⛔ The consumer at
        :meth:`orpheus.derivations.continuous.fn_method.moment_space`
        asserted ``converged=True`` here until 2026-08-09, justified by
        the comment *"bisection always converges given bracket"*.  That
        is false: ``max_bisect`` is a caller-tunable parameter and is
        forwarded, so the loop can exhaust.  The bracket width is the
        witness — not ``determinant_residual``, whose spec is only
        "should be small".
    """

    a_critical_mfp: float
    c: float
    N: int
    nu0: float
    xi_collocation: np.ndarray
    coefficients_a: np.ndarray
    determinant_residual: complex
    n_bracket_points: int
    converged: bool


def _build_fn_matrix(a: float, N: int, xis: np.ndarray, c: float) -> np.ndarray:
    r"""Assemble the slab F_N matrix :math:`M(a)` via the shared assembler.

    Row :math:`\beta` (collocation point :math:`\xi_\beta`) has entries
    :math:`M_{\beta\alpha} = B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta}
    A_\alpha(\xi_\beta)` for :math:`\alpha = 0, \ldots, N`. The
    :math:`\xi = 0` row uses :math:`e^{-\infty} = 0`, leaving only the
    :math:`B_\alpha(0)` part.

    Implemented as a thin wrapper over :func:`...core.fn_matrix.assemble_fn_matrix`
    with ``geometry_sign = +1`` (slab BC, Siewert-Benoist Eq. 4 / Grandjean-Siewert
    Eq. 25). Sphere reuses the same assembler with ``geometry_sign = -1``.
    """
    return assemble_fn_matrix(c, a, geometry_sign=+1, n_modes=N, xis=xis)


def _det_M(a: float, N: int, xis: np.ndarray, c: float) -> complex:
    """Return :math:`\\det M(a)` as a complex number."""
    M = _build_fn_matrix(a, N, xis, c)
    return complex(np.linalg.det(M))


def solve_fn_slab_bare_critical(
    c: float,
    *,
    n_modes: int = 5,
    a_min: float = 0.05,
    a_max: float = 20.0,
    n_bracket: int = 400,
    bisect_tol: float = 1e-12,
    max_bisect: int = 80,
) -> SlabFNResult:
    r"""Solve the critical half-thickness of a bare 1G slab via the
    F_N method.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision,
        :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`. Must satisfy
        :math:`c > 1` (multiplying medium); the F_N method is also
        defined for :math:`c < 1` half-space problems but those are
        not critical-slab problems.
    n_modes : int, default 5
        F_N order, :math:`N`. The system size is :math:`(N+1) \times
        (N+1)`. :math:`N = 5` (F_5) gives 4-5 sig figs accuracy across
        the Grandjean-Siewert critical-thickness table; :math:`N = 8`
        gets to ~1e-5 absolute on the half-thickness.
    a_min, a_max : float, default 0.05, 20.0
        Bracket for the bisection on :math:`a`. The bracket is
        sampled at :math:`n_{\rm bracket}` points and the first
        sign change of :math:`\Re \det M(a)` is bisected.
    n_bracket : int, default 400
        Number of bracket-scan points.
    bisect_tol : float, default 1e-12
        Absolute tolerance on :math:`a` for the bisection step.
    max_bisect : int, default 80
        Maximum bisection iterations.

    Returns
    -------
    :class:`SlabFNResult`

    Notes
    -----
    The collocation prescription follows Grandjean-Siewert Part II
    Section III: :math:`\xi_0 = \nu_0`, :math:`\xi_1 = 0`,
    :math:`\xi_2 = 1`, and the remaining :math:`N - 2` points equally
    spaced strictly inside :math:`(0, 1)`. Different prescriptions
    give different convergence rates but converge to the same limit.
    """
    if c <= 1.0:
        raise ValueError(
            f"F_N critical slab requires c > 1 (multiplying medium), got c={c}"
        )
    if n_modes < 0:
        raise ValueError(f"n_modes must be ≥ 0, got {n_modes}")

    # Case discrete eigenvalue magnitude (positive real); for c > 1,
    # actual ν_0 = i u_0.
    u0, is_imag = case_nu0(c)
    if not is_imag:
        # Should not happen given c > 1 check, but guard anyway.
        raise RuntimeError("internal: expected supercritical c>1 with imaginary nu_0")
    nu0_complex = complex(0.0, u0)

    # Collocation points: ξ_0 = ν_0 = i u_0 (complex), ξ_1 = 0, ξ_2 = 1,
    # remaining equally spaced in (0, 1).
    real_xis = collocation_points(0.0, n_modes)  # placeholder for xi_0
    xis = np.empty(n_modes + 1, dtype=complex)
    xis[0] = nu0_complex
    xis[1:] = real_xis[1:].astype(complex)

    # Bracket scan on the real part of det M(a). Some intermediate
    # `a` values give matrices with vanishing determinants (genuine
    # zeros that we WANT to bracket) or near-singular matrices that
    # numpy reports as division warnings. Suppress the warnings — the
    # zeros are what we're looking for.
    import warnings as _warnings
    a_scan = np.linspace(a_min, a_max, n_bracket)
    real_det_scan = np.empty(n_bracket)
    with _warnings.catch_warnings():
        _warnings.filterwarnings("ignore", category=RuntimeWarning)
        for i, a in enumerate(a_scan):
            try:
                real_det_scan[i] = _det_M(a, n_modes, xis, c).real
            except (np.linalg.LinAlgError, ValueError, FloatingPointError):
                real_det_scan[i] = float("nan")

    # Find the first sign change.
    valid = ~np.isnan(real_det_scan)
    sgn = np.sign(real_det_scan[valid])
    if np.all(sgn == sgn[0]):
        raise RuntimeError(
            f"No sign change of det M(a) found in [{a_min}, {a_max}] "
            f"with N={n_modes}, c={c}. Try widening the bracket or "
            f"increasing n_bracket."
        )
    # Find first crossing index relative to the valid mask.
    valid_idx = np.where(valid)[0]
    diffs = np.diff(np.sign(real_det_scan[valid_idx]))
    first_crossing_in_valid = int(np.where(diffs != 0)[0][0])
    i0 = valid_idx[first_crossing_in_valid]
    i1 = valid_idx[first_crossing_in_valid + 1]
    a_lo, a_hi = a_scan[i0], a_scan[i1]
    f_lo, f_hi = real_det_scan[i0], real_det_scan[i1]

    # Bisection on a.
    for _ in range(max_bisect):
        if a_hi - a_lo < bisect_tol * max(1.0, a_lo):
            break
        a_mid = 0.5 * (a_lo + a_hi)
        f_mid = _det_M(a_mid, n_modes, xis, c).real
        if f_mid * f_lo < 0:
            a_hi, f_hi = a_mid, f_mid
        else:
            a_lo, f_lo = a_mid, f_mid

    a_crit = 0.5 * (a_lo + a_hi)
    M_crit = _build_fn_matrix(a_crit, n_modes, xis, c)
    det_residual = complex(np.linalg.det(M_crit))

    # Recover the eigenvector (F_N coefficients a_α) via SVD.
    # The null vector of M is the F_N expansion of ψ(-a, -μ).
    # Numerical null space: smallest singular value's right vector.
    U, s, Vh = np.linalg.svd(M_crit)
    coefficients = Vh.conj().T[:, -1]
    # Normalize so a_0 = 1.
    if abs(coefficients[0]) > 1e-30:
        coefficients = coefficients / coefficients[0]

    # The loop's OWN stop predicate on the final bracket (#340) — not
    # "did we break", and not an appeal to bisection's guaranteed
    # asymptotic convergence, which says nothing about a FINITE
    # ``max_bisect``.
    converged = bool(a_hi - a_lo < bisect_tol * max(1.0, a_lo))

    return SlabFNResult(
        a_critical_mfp=float(a_crit),
        c=c,
        N=n_modes,
        nu0=u0,
        xi_collocation=xis,
        coefficients_a=coefficients,
        determinant_residual=det_residual,
        n_bracket_points=n_bracket,
        converged=converged,
    )


def fn_slab_flux_at_x_cosine_only(
    result: SlabFNResult,
    x: float | np.ndarray,
) -> float | np.ndarray:
    r"""**Approximate** scalar flux using only the discrete Case mode
    (Mitsis 1963 first-order asymptotic).

    For the supercritical bare slab the discrete-mode contribution is

    .. math::

       \phi^{(0)}(x) = \cos(x/u_0)

    where :math:`u_0 = |\nu_0|` is the magnitude of the imaginary-
    axis Case discrete eigenvalue. **This omits the continuum
    contribution** :math:`\int_0^1 A(\nu) e^{-b/\nu}\cosh(x/\nu)\,d\nu`
    (KLL Eq. 7), which adds :math:`\mathcal{O}(10\%)` to :math:`\phi(x)`
    near the slab surface. Hence this entry point is exposed as
    diagnostic only — full flux reconstruction within the F_N
    framework is deferred (see GitHub issue tracking; the angular flux
    half-range :math:`\psi(\pm a, \mp\mu) = \sum a_\alpha \mu^\alpha`
    representation lives in :attr:`SlabFNResult.coefficients_a` and
    a future helper can integrate it self-consistently).

    Use case for tests: agreement with the cosine-only term is
    expected to ~1% near the slab center and degrade to ~20% at the
    edge — useful as a sanity check, not as an L1 verification.

    Parameters
    ----------
    result : :class:`SlabFNResult`
    x : float or array
        Position(s) in :math:`[-a, a]` (mfp). Symmetric in ``x``.

    Returns
    -------
    Same shape as ``x``: :math:`\cos(x/u_0)`.
    """
    u0 = result.nu0
    x_arr = np.asarray(x, dtype=float)
    phi = np.cos(x_arr / u0)
    if np.isscalar(x):
        return float(phi)
    return phi
