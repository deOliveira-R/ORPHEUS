r"""Case dispersion-relation roots for the F_N method.

The Case-Zweifel singular eigenfunction expansion of the 1G isotropic
transport equation has continuum eigenvalues :math:`\nu \in (-1, 1)`
and a pair of discrete eigenvalues :math:`\pm\nu_0` (real, with
:math:`|\nu_0| > 1` for :math:`c \neq 1`) determined by the dispersion
relation

.. math::

   \Lambda(\nu_0) = 1 - c\,\nu_0\,\mathrm{atanh}(1/\nu_0) = 0
   \qquad (c < 1)

or equivalently for :math:`c > 1` with :math:`\nu_0 = i|\nu_0|`
purely imaginary,

.. math::

   1 - c\,|\nu_0|\,\mathrm{atan}(1/|\nu_0|) = 0 .

Both branches are handled by the same scalar findroot below; we work
with the real positive variable :math:`\nu_0` for :math:`c < 1` and
the real positive variable :math:`u_0 = |\nu_0|` for :math:`c > 1`,
keeping the API uniform.

References
----------

* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156, Eq. 10-11.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161, Eq. 6-7.
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94 — gives
  the form :math:`1 - c\,\nu \tan^{-1}(1/\nu) = 0` for the supercritical
  case using :math:`|\nu_0|`. Same equation, different symbol.
* Case 1960, *Ann. Phys.* **9**, 1.
"""
from __future__ import annotations

import math


def case_dispersion_root_subcritical(c: float, *, tol: float = 1e-15) -> float:
    r"""Solve :math:`1 - c\,\nu_0\,\mathrm{atanh}(1/\nu_0) = 0` for
    :math:`\nu_0 > 1`, valid for :math:`0 < c < 1`.

    Uses bisection on a guaranteed bracket :math:`[1+\epsilon, \nu_*]`
    where :math:`\nu_*` is grown until :math:`\Lambda(\nu_*) > 0`.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision, :math:`(\Sigma_s +
        \nu\Sigma_f)/\Sigma_t`. Must satisfy :math:`0 < c < 1`.
    tol : float
        Absolute tolerance on :math:`\nu_0` (default 1e-15 — bisection
        gets to machine precision in ~50 iterations).

    Returns
    -------
    float
        The unique :math:`\nu_0 > 1` satisfying the dispersion relation.

    Notes
    -----
    For :math:`c \to 1^-` the root approaches :math:`\nu_0 \to \infty`;
    for :math:`c \to 0^+` it approaches :math:`\nu_0 \to 1^+`. The
    function is smooth and monotone on the bracket so bisection is
    sufficient (Newton would also work; bisection is chosen for
    robustness without derivative tuning).
    """
    if not (0.0 < c < 1.0):
        raise ValueError(f"subcritical c required, got c={c}")

    def Lam(nu: float) -> float:
        return 1.0 - c * nu * math.atanh(1.0 / nu)

    # Bracket: at nu = 1 + eps, atanh(1/nu) → +∞ so Λ → -∞.
    lo = 1.0 + 1e-12
    hi = 2.0
    while Lam(hi) < 0.0:
        hi *= 2.0
        if hi > 1.0e12:
            raise RuntimeError(
                f"Could not bracket dispersion root for c={c}"
            )
    # Bisection.
    while hi - lo > tol * max(1.0, lo):
        mid = 0.5 * (lo + hi)
        if Lam(mid) < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def case_dispersion_root_supercritical(c: float, *, tol: float = 1e-15) -> float:
    r"""Solve :math:`1 - c\,u_0\,\mathrm{atan}(1/u_0) = 0` for :math:`u_0 > 0`,
    where :math:`u_0 = |\nu_0|` and :math:`\nu_0 = i u_0` is purely
    imaginary, valid for :math:`c > 1`.

    Uses bisection on :math:`[\epsilon, u_*]` where :math:`u_*` is
    grown until :math:`\Lambda(u_*) > 0`, exploiting the monotonicity
    of :math:`u \mapsto u\,\mathrm{atan}(1/u)`.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision. Must satisfy :math:`c > 1`.
    tol : float
        Absolute tolerance on :math:`u_0`.

    Returns
    -------
    float
        :math:`u_0 = |\nu_0|`. The Case discrete root is :math:`\nu_0 = i u_0`.

    Notes
    -----
    For :math:`c \to 1^+`, :math:`u_0 \to \infty`; for :math:`c \to \infty`,
    :math:`u_0 \to 0^+`. Some F_N references work with :math:`\nu_0` complex
    directly; we keep it real positive for numerical convenience.
    """
    if c <= 1.0:
        raise ValueError(f"supercritical c required (c>1), got c={c}")

    def Lam(u: float) -> float:
        return 1.0 - c * u * math.atan(1.0 / u)

    # As u → 0+: u·atan(1/u) → u·(π/2) → 0, so Λ → 1 (positive).
    # As u → ∞: u·atan(1/u) → u·(1/u - 1/(3u^3)+...) → 1, so Λ → 1 - c < 0.
    # Hence Λ has exactly one zero on (0, ∞) for c > 1.
    lo = 1e-12
    hi = 1.0
    while Lam(hi) > 0.0:
        hi *= 2.0
        if hi > 1.0e12:
            raise RuntimeError(
                f"Could not bracket supercritical dispersion root for c={c}"
            )
    while hi - lo > tol * max(1.0, lo):
        mid = 0.5 * (lo + hi)
        if Lam(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def case_nu0(c: float) -> tuple[float, bool]:
    r"""Unified API for the Case discrete eigenvalue magnitude.

    Returns ``(nu0, is_imag)`` where:

    * ``is_imag = False`` (subcritical, :math:`c < 1`): :math:`\nu_0`
      is real and :math:`\nu_0 > 1`.
    * ``is_imag = True`` (supercritical, :math:`c > 1`): the returned
      number is :math:`|\nu_0|` and the actual Case eigenvalue is
      :math:`\nu_0 = i|\nu_0|`.

    For :math:`c = 1` the dispersion relation is degenerate; raises.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision.

    Returns
    -------
    tuple[float, bool]
    """
    if c == 1.0:
        raise ValueError(
            "c = 1 is the degenerate case (no leakage limit); "
            "the dispersion relation has no finite root."
        )
    if c < 1.0:
        return case_dispersion_root_subcritical(c), False
    return case_dispersion_root_supercritical(c), True
