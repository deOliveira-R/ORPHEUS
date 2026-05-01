r"""Shifted-Legendre monomial coefficients (canonical implementation).

Provides the single source of truth for the monomial-basis coefficients
:math:`(c_n^0, c_n^1, \ldots, c_n^n)` of the shifted Legendre polynomial

.. math::

   \tilde P_n(\mu) \;=\; P_n(2\mu - 1)
                   \;=\; \sum_{k=0}^{n} c_n^k\,\mu^k.

These coefficients are used by:

- :func:`peierls_geometry.compute_P_esc_*_mode` and friends — the
  rank-:math:`N` Knyazev :math:`\mathrm{Ki}_{2+k}` and :math:`E_{n+2}`
  expansions in cylinder and slab geometries.
- The lifted SymPy derivation modules
  :mod:`orpheus.derivations.continuous.peierls.origins.cylinder_knyazev` and
  :mod:`orpheus.derivations.continuous.peierls.origins.specular.slab` consume this same
  canonical coefficient source (the workbench duplicates were
  removed when those scripts were lifted into the package).

Implementation
--------------

The canonical implementation uses :func:`numpy.polynomial.legendre.leg2poly`
to convert standard-Legendre coefficients to monomial coefficients, then
expands :math:`(2\mu - 1)^j` binomially. This avoids SymPy at runtime
(production callers like the rank-:math:`N` cylinder primitive call this
inside the angular-quadrature loop) while remaining bit-exact for
:math:`n \lesssim 20` in float64.

Historical note: prior to the workbench-to-package lift, the
``peierls_specular_slab.py`` and ``peierls_cylinder_3d_mode_n.py``
derivation scripts in ``scratch/derivations/`` carried duplicate
SymPy-based copies of these coefficients (computed via
:func:`sympy.legendre` + :func:`sympy.Poly`). Those duplicates were
removed when the scripts were lifted to
:mod:`orpheus.derivations.continuous.peierls.origins.specular.slab` and
:mod:`orpheus.derivations.continuous.peierls.origins.cylinder_knyazev` respectively;
the lifted modules now import this canonical implementation directly.
"""

from __future__ import annotations

import functools
import math

import numpy as np


@functools.lru_cache(maxsize=64)
def shifted_legendre_monomial_coefs(n: int) -> tuple[float, ...]:
    r"""Monomial-basis coefficients :math:`(c_n^0, c_n^1, \ldots, c_n^n)`
    of :math:`\tilde P_n(\mu) = P_n(2\mu - 1) = \sum_k c_n^k\,\mu^k`.

    Computed once and cached per polynomial order. Returns coefficients
    in ascending monomial order (constant term first).

    Examples
    --------
    >>> shifted_legendre_monomial_coefs(0)
    (1.0,)
    >>> shifted_legendre_monomial_coefs(1)
    (-1.0, 2.0)
    >>> shifted_legendre_monomial_coefs(2)
    (1.0, -6.0, 6.0)
    """
    if n < 0:
        raise ValueError(f"n must be non-negative, got {n}")
    # Build P_n in standard form via Bonnet, then substitute x = 2µ - 1.
    # Use numpy.polynomial: leg2poly converts Legendre coefficients to
    # monomial coefficients in x, then expand binomially.
    # SymPy is overkill at runtime; do this with numpy's Polynomial class.
    leg_coefs = np.zeros(n + 1)
    leg_coefs[n] = 1.0
    # P_n(x) coefficients in monomial basis (x ascending order)
    px_ascending = np.polynomial.legendre.leg2poly(leg_coefs)
    # Substitute x = 2µ - 1: expand each x^j as (2µ - 1)^j
    coefs = np.zeros(n + 1)
    for j, aj in enumerate(px_ascending):
        if aj == 0.0:
            continue
        # (2µ - 1)^j = sum_{m=0}^j C(j,m) (2µ)^m (-1)^(j-m)
        for m in range(j + 1):
            binom = math.comb(j, m)
            coefs[m] += aj * binom * (2.0**m) * ((-1.0) ** (j - m))
    return tuple(float(c) for c in coefs)
