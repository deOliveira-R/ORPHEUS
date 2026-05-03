r"""Core machinery for the Case singular-eigenfunction method (Atalay 1997).

Modules
-------

* :mod:`.dispersion` — Atalay Eq. 12 quadratic in c; finds the purely
  imaginary discrete Case eigenvalue :math:`\nu_0 = i u_0`.
* :mod:`.x_function` — Atalay Eq. 40 X-function via mpmath quadrature
  with successive subdivision near the right endpoint.
* :mod:`.extrapolated_endpoint` — Atalay Eq. 42 :math:`z_0(c, f_1)`.
* :mod:`.half_range` — Atalay Eqs. 38 (slab) and 51 (sphere)
  :math:`K_j(c, R)` and :math:`L_j(c, R)` integral moments.
"""
from __future__ import annotations

from .dispersion import (
    case_atalay_dispersion_quadratic_coeffs,
    case_atalay_u0,
    nu_bar_atalay,
)
from .extrapolated_endpoint import atalay_z0
from .half_range import atalay_K_moments, atalay_L_moments, T_function, T1_function
from .x_function import atalay_X_function

__all__ = [
    # dispersion
    "case_atalay_u0",
    "case_atalay_dispersion_quadratic_coeffs",
    "nu_bar_atalay",
    # x function
    "atalay_X_function",
    # extrapolated endpoint
    "atalay_z0",
    # half range moments
    "atalay_K_moments",
    "atalay_L_moments",
    "T_function",
    "T1_function",
]
