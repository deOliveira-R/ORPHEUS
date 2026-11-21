r"""The Bickley-Naylor function :math:`\mathrm{Ki}_3`, the kernel of cylindrical collision probability.

.. math::

    \mathrm{Ki}_3(\tau) \;=\; \int_0^{\pi/2} \cos^2\theta\, e^{-\tau/\cos\theta}\, d\theta
    \;=\; \int_0^\infty (1 + u^2)^{-2}\, e^{-\tau\sqrt{1 + u^2}}\, du ,

the second form being the substitution :math:`u = \tan\theta`, which removes the essential singularity at :math:`\theta = \pi/2`.

:func:`ki3` evaluates :math:`e^{-\tau}` times the degree-63 Chebyshev interpolant of the scaled kernel :math:`f(\tau) = e^{\tau}\,\mathrm{Ki}_3(\tau)`
on :math:`[0, 50]` (the scaling turns the exponentially decaying tail into a slowly varying function that 64 points fit to about
:math:`3\times 10^{-6}` relative), and clamps to 0 above 50 (:math:`\mathrm{Ki}_3(50) \approx 3\times 10^{-23}`). The interpolant is rebuilt from
the integral at 30 decimal digits, lazily and once, never stored: a stored table of coefficients is an untraceable constant.

This is a special function, not a space, field or operator. It is the same interpolant that the reference collision-probability derivation holds
(a twin, equal bit for bit, pinned by ``tests/gates/derivations/test_cp_ki3_copy.py``), so a defect of the shared approximation is invisible to that
cross-check.
"""

from __future__ import annotations

import functools

import mpmath
import numpy as np
from numpy.typing import ArrayLike, NDArray

__all__ = ["ki3", "ki3_scaled_chebyshev"]

_KI3_TAU_MAX = 50.0
_KI3_DEG = 63  # 64 Chebyshev points of the first kind
_KI3_DPS = 30


def _ki3_integral(tau: float) -> mpmath.mpf:
    r""":math:`\mathrm{Ki}_3(\tau)` as :math:`\int_0^\infty (1+u^2)^{-2} e^{-\tau\sqrt{1+u^2}}\,du` at ``_KI3_DPS`` digits."""
    with mpmath.workdps(_KI3_DPS):
        tau_mp = mpmath.mpf(float(tau))

        def integrand(u):
            s = mpmath.sqrt(1 + u * u)
            return s ** -4 * mpmath.exp(-tau_mp * s)

        return mpmath.quad(integrand, [0, mpmath.inf])


@functools.lru_cache(maxsize=1)
def ki3_scaled_chebyshev() -> np.polynomial.Chebyshev:
    r"""The Chebyshev interpolant of :math:`e^{\tau}\,\mathrm{Ki}_3(\tau)` on :math:`[0, 50]`, degree 63, built once (about 0.3 s).

    Its coefficients are read-only: the object is cached and shared, so a write would change every later kernel call."""

    def scaled(tau: np.ndarray) -> np.ndarray:
        return np.array([float(_ki3_integral(float(t))) * float(np.exp(t)) for t in tau])

    poly = np.polynomial.Chebyshev.interpolate(scaled, deg=_KI3_DEG, domain=[0.0, _KI3_TAU_MAX])
    poly.coef.flags.writeable = False
    return poly


def ki3(tau: ArrayLike) -> NDArray[np.float64]:
    r""":math:`\mathrm{Ki}_3(\tau)` for optical thicknesses ``tau >= 0``: the output has the input's shape; 0 exactly for ``tau > 50``."""
    poly = ki3_scaled_chebyshev()
    tau = np.asarray(tau, dtype=float)
    out = np.zeros_like(tau, dtype=float)
    mask = tau <= _KI3_TAU_MAX
    if np.any(mask):
        out[mask] = poly(tau[mask]) * np.exp(-tau[mask])
    return out
