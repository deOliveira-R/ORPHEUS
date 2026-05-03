r"""Branch-1 SymPy algebra-of-record for the singular-eigenfunction
expansion family.

Modules
-------

* :mod:`.cylinder_derivations` — Westfall-Metcalf 1973 cylinder
  derivations: bare-cylinder integral-equation reduction, dispersion
  function (re-derivation matching Case 1960 / Westfall-Metcalf Eq. 18),
  Bessel-Wronskian identity used in the integrodifferential reduction
  (Eq. 9), pseudo-eigenfunction structure (Eq. 17 + 19), bare-cylinder
  closure (Eq. 32 + 30 simplification).

* :mod:`.slab_sphere_derivations` — Atalay 1997 slab + sphere derivations:
  dispersion relation reduction Eq. 11→12 (linearly anisotropic), parity
  flip Eqs. 47-49 mapping slab-odd → sphere, half-range relations
  Eqs. 28-31, Fredholm-form prefactor Eq. 32, criticality conditions
  Eqs. 46 (slab) / 54 (sphere), X-function Eq. 40, extrapolated endpoint
  Eq. 42, validity bound Eq. 5.

The algebra-of-record discipline (skill ``algebra-of-record``) requires
each verifiable claim to be one of:

* **State 1A** — closed-form SymPy identity (algebraic verification).
* **State 1B** — semi-analytical (single integral, mpmath/scipy quadrature).
* **State 1C** — MMS (manufactured solution).

Atalay's machinery is mostly **State 1A** at the algebra layer (dispersion
quadratic, parity flip, half-range relations, Fredholm prefactor,
criticality conditions). The X(μ) function and extrapolated endpoint
:math:`z_0` are **State 1B** (single mpmath integrals).
"""
from __future__ import annotations

from .slab_sphere_derivations import (
    derive_atalay_critical_slab_eq46,
    derive_atalay_critical_sphere_eq54_via_parity_flip,
    derive_atalay_dispersion_linear_anisotropic,
    derive_atalay_extrapolated_endpoint_eq42,
    derive_atalay_fredholm_form_eq27_to_eq32,
    derive_atalay_half_range_eqs28_to_31,
    derive_atalay_symmetry_conditions_eq13_14_47_to_49,
    derive_atalay_validity_bound,
    derive_atalay_validity_bound_eq5,
    derive_atalay_x_function_eq40,
)

__all__ = [
    "derive_atalay_dispersion_linear_anisotropic",
    "derive_atalay_symmetry_conditions_eq13_14_47_to_49",
    "derive_atalay_half_range_eqs28_to_31",
    "derive_atalay_fredholm_form_eq27_to_eq32",
    "derive_atalay_critical_slab_eq46",
    "derive_atalay_critical_sphere_eq54_via_parity_flip",
    "derive_atalay_extrapolated_endpoint_eq42",
    "derive_atalay_x_function_eq40",
    "derive_atalay_validity_bound",
    "derive_atalay_validity_bound_eq5",
]
