r"""Closed-form moment integrals :math:`B_\alpha(\xi)` and
:math:`A_\alpha(\xi)` for the slab F_N method.

Per Grandjean-Siewert 1979 (Part II) Eqs 10 + 20 (and equivalently
Siewert-Benoist 1979 (Part I) Eqs 29 + 48), the slab F_N matrix
entries are

.. math::

   B_\alpha(\xi) &= \frac{2}{c\xi}\int_0^1 \mu^{\alpha+1}\phi(\xi,\mu)\,d\mu,
   \\
   A_\alpha(\xi) &= \frac{2}{c\xi}\int_0^1 \mu^{\alpha+1}\phi(-\xi,\mu)\,d\mu .

Both reduce to closed-form recursions in :math:`\alpha`:

.. math::

   B_0(\xi) &= \frac{2}{c} - 1 - \xi\,\log\!\left(1 + \frac{1}{\xi}\right), \\
   B_\alpha(\xi) &= \xi B_{\alpha-1}(\xi) - \frac{1}{\alpha+1}, \qquad \alpha \ge 1, \\
   A_0(\xi) &= 1 - \xi\,\log\!\left(1 + \frac{1}{\xi}\right), \\
   A_\alpha(\xi) &= -\xi A_{\alpha-1}(\xi) + \frac{1}{\alpha+1}, \qquad \alpha \ge 1 .

The :math:`\log(1+1/\xi)` factor is benign for :math:`\xi \in [0, 1]`
plus :math:`\xi = \nu_0 > 1`: at :math:`\xi = 0` it diverges as
:math:`\log(1/\xi)` but is multiplied by :math:`\xi`, giving a finite
limit :math:`\xi \log(1+1/\xi) \to 0`. We special-case :math:`\xi = 0`
to avoid evaluating ``log(inf)``.

References
----------

* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156, Eqs. 29+48.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161, Eqs. 10+20.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np


def _xi_log_one_plus_inv_xi(xi: float) -> float:
    r"""Return :math:`\xi\,\log(1 + 1/\xi)` with the :math:`\xi = 0`
    limit handled (returns 0 there)."""
    if xi == 0.0:
        return 0.0
    return xi * math.log1p(1.0 / xi)


def B_alpha(alpha: int, xi: float, c: float) -> float:
    r"""Compute :math:`B_\alpha(\xi)` via the closed-form recursion.

    Parameters
    ----------
    alpha : int
        Moment order, :math:`\alpha \ge 0`.
    xi : float
        Collocation point, :math:`\xi \in [0, 1]` or :math:`\xi = \nu_0`.
    c : float
        Mean number of secondaries per collision.

    Returns
    -------
    float
    """
    if alpha < 0:
        raise ValueError(f"alpha must be ≥ 0, got {alpha}")
    B = (2.0 / c) - 1.0 - _xi_log_one_plus_inv_xi(xi)
    for n in range(1, alpha + 1):
        B = xi * B - 1.0 / (n + 1)
    return B


def A_alpha(alpha: int, xi: float, c: float) -> float:  # noqa: ARG001
    r"""Compute :math:`A_\alpha(\xi)` via the closed-form recursion.

    Note that :math:`A_\alpha` does NOT depend on :math:`c` — the
    parameter is kept in the signature for API symmetry with
    :func:`B_alpha`.

    Parameters
    ----------
    alpha : int
        Moment order, :math:`\alpha \ge 0`.
    xi : float
        Collocation point, :math:`\xi \in [0, 1]` or :math:`\xi = \nu_0`.
    c : float
        Unused but kept for symmetric signature with :func:`B_alpha`.

    Returns
    -------
    float
    """
    if alpha < 0:
        raise ValueError(f"alpha must be ≥ 0, got {alpha}")
    A = 1.0 - _xi_log_one_plus_inv_xi(xi)
    for n in range(1, alpha + 1):
        A = -xi * A + 1.0 / (n + 1)
    return A


def B_alpha_array(N: int, xi: float, c: float) -> np.ndarray:
    r"""Compute :math:`(B_0, B_1, \ldots, B_N)` at :math:`\xi` for fixed
    :math:`c`. More efficient than calling :func:`B_alpha` N+1 times."""
    out = np.empty(N + 1, dtype=float)
    out[0] = (2.0 / c) - 1.0 - _xi_log_one_plus_inv_xi(xi)
    for n in range(1, N + 1):
        out[n] = xi * out[n - 1] - 1.0 / (n + 1)
    return out


def A_alpha_array(N: int, xi: float, c: float) -> np.ndarray:  # noqa: ARG001
    r"""Compute :math:`(A_0, A_1, \ldots, A_N)` at :math:`\xi`.
    :math:`A_\alpha` does not depend on :math:`c`."""
    out = np.empty(N + 1, dtype=float)
    out[0] = 1.0 - _xi_log_one_plus_inv_xi(xi)
    for n in range(1, N + 1):
        out[n] = -xi * out[n - 1] + 1.0 / (n + 1)
    return out


def collocation_points(nu0: float, N: int) -> np.ndarray:
    r"""Equally-spaced F_N collocation points per Grandjean-Siewert
    Part II "we let :math:`\xi_0 = \nu_0`, :math:`\xi_1 = 0`,
    :math:`\xi_2 = 1`, and the remaining :math:`\xi_\beta` are spaced
    equally in the interval :math:`[0, 1]`."

    Returns an ``(N+1,)`` array.

    For :math:`N = 0`: just ``(nu0,)``.
    For :math:`N = 1`: ``(nu0, 0)``.
    For :math:`N = 2`: ``(nu0, 0, 1)``.
    For :math:`N \ge 3`: ``(nu0, 0, 1)`` followed by :math:`N-2`
    interior points equally spaced strictly inside :math:`(0, 1)`.

    The exact "remaining ξ_β are spaced equally in [0, 1]" prescription
    is read here as: the remaining :math:`N - 2` interior points are
    placed at :math:`(k+1)/(N-1)` for :math:`k = 0, \ldots, N-3`,
    sandwiched between the two endpoints already at 0 and 1. This
    matches what Grandjean-Siewert 1979 Table III-X imply (their
    F_5 = 6-point grid for [0, 1] with endpoints at 0 and 1 gives
    points at 0, 0.25, 0.5, 0.75, 1.0 plus :math:`\nu_0`).
    """
    if N < 0:
        raise ValueError(f"N must be ≥ 0, got {N}")
    pts = np.empty(N + 1, dtype=float)
    pts[0] = nu0
    if N >= 1:
        pts[1] = 0.0
    if N >= 2:
        pts[2] = 1.0
    if N >= 3:
        # The remaining N-2 collocation points fill in [0,1] equally
        # between the existing endpoints at 0 (ξ_1) and 1 (ξ_2). For
        # F_5 (N=5) Grandjean-Siewert Table III implies points at
        # {0, 0.25, 0.50, 0.75, 1.0}; that's 0, 1, plus 3 interior
        # at multiples of 1/4. General: divide [0,1] into N-1 equal
        # subintervals, giving N points; we already have 0 and 1,
        # so the N-2 interior are at k/(N-1) for k=1, ..., N-2.
        pts[3:] = np.array(
            [k / (N - 1) for k in range(1, N - 1)]
        )
    return pts


def collocation_points_explicit(nu0: float, N: int) -> Sequence[float]:
    """Same as :func:`collocation_points` but as a list, for diagnostics."""
    return list(collocation_points(nu0, N))
