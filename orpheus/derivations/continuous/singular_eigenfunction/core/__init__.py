r"""Shared primitives for the singular-eigenfunction expansion family.

The dispersion-function root primitive
:func:`...fn_method.core.dispersion.case_nu0` is shared from the F_N
package; this is acceptable cross-package reuse because the dispersion
function :math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu) = 0`
is a **medium property** common to every singular-eigenfunction
treatment of the 1G isotropic transport equation, not an
algorithmic-pillar primitive of either package.

This module also hosts the **linearly-anisotropic** Atalay 1997 primitives
that the slab + sphere solvers depend on:

* :mod:`.dispersion` — Atalay Eq. 12 quadratic in c; finds the purely
  imaginary discrete Case eigenvalue :math:`\nu_0 = i u_0`.
* :mod:`.x_function` — Atalay Eq. 40 X-function via mpmath quadrature
  with successive subdivision near the right endpoint.
* :mod:`.extrapolated_endpoint` — Atalay Eq. 42 :math:`z_0(c, f_1)`.
* :mod:`.half_range` — Atalay Eqs. 38 (slab) and 51 (sphere)
  :math:`K_j(c, R)` and :math:`L_j(c, R)` integral moments.

The cylinder Branch-2 solver (Westfall-Metcalf 1973) currently does not
share these primitives — its kernel reduces to modified Bessel functions
:math:`K_0(t/\mu)\,I_0(r/\mu)` rather than exponentials, so a separate
Bessel-Wronskian primitive lives inline in :mod:`.cylinder.one_group`.
A future ``bessel_radial.py`` core module is the natural reservation
point if cylinder + slab/sphere ever need to share radial machinery.
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
