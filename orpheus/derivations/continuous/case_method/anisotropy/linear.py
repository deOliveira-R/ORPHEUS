r"""Linearly-anisotropic scattering utilities.

The Atalay 1997 framework is restricted to linearly-anisotropic
scattering (P_1) in the **single-pair-of-discrete-modes** regime,
parameterised by the mean cosine :math:`f_1 \in [0, 1)`.

Validity bound (Atalay Eq 5):

.. math::

   c \le 1 + \frac{1}{3 f_1}

Outside this band, complex-conjugate eigenvalue pairs appear in the
dispersion relation (Dahl-Sjöstrand 1979; Kohut 1993), and the
first-order Fredholm iteration breaks. For :math:`f_1 = 0`
(isotropic) the bound is trivially satisfied for all :math:`c`.
"""
from __future__ import annotations

import math


def atalay_validity_max_c(f1: float) -> float:
    r"""Maximum allowed :math:`c` for the Atalay 1997 framework.

    Per Atalay Eq 5: :math:`c_{\rm max} = 1 + 1/(3 f_1)`.

    For :math:`f_1 = 0`, returns :math:`+\infty`.
    """
    if f1 < 0:
        raise ValueError(f"f1 must be ≥ 0, got {f1}")
    if f1 == 0.0:
        return math.inf
    return 1.0 + 1.0 / (3.0 * f1)


def check_atalay_validity(c: float, f1: float) -> None:
    r"""Raise ``ValueError`` if (c, f_1) violate Atalay Eq 5.

    Convenience pre-condition for the slab/sphere solvers.
    """
    if c <= 1.0:
        raise ValueError(
            f"c = {c} ≤ 1: Atalay framework requires multiplying medium c > 1."
        )
    c_max = atalay_validity_max_c(f1)
    if c > c_max:
        raise ValueError(
            f"c = {c} violates Atalay Eq 5 validity bound "
            f"c ≤ 1 + 1/(3 f_1) = {c_max} for f_1 = {f1}. "
            f"Outside this band complex eigenvalues appear "
            f"(Dahl-Sjöstrand 1979 / Kohut 1993)."
        )
