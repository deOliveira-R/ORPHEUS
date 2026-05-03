r"""Branch-1 SymPy algebra-of-record for the Case singular-eigenfunction
method (Atalay 1997).

Each ``derive_*()`` function returns a dict with at minimum a
``"name"`` field, the symbolic expressions verified, and a
``"pass": bool`` flag. The :mod:`tests.derivations.test_case_method_*`
test gates pin one foundation-tagged test per ``derive_*()``.

The algebra-of-record discipline (skill ``algebra-of-record``) requires
each verifiable claim to be one of:

* **State 1A** — closed-form SymPy identity (algebraic verification).
* **State 1B** — semi-analytical (single integral, mpmath/scipy quadrature).
* **State 1C** — MMS (manufactured solution).

Atalay's machinery is mostly **State 1A** at the algebra layer:

* Dispersion relation Eq. 11 reduces to quadratic in c (Eq. 12) — closed-form.
* Parity flip Eqs. 47-49 → odd-mode boundary conditions for sphere — algebraic identity.
* Half-range relations Eqs. 28-31 share Wiener-Hopf X-function structure — algebraic.
* Eq. 32 first-order Fredholm form prefactor + bracket structure — algebraic identity.
* Criticality conditions Eqs. 46 / 54 — algebraic identity (taking
  log of complex-conjugate quotient = i·(±π/2 - 2·arg)).
* Geometry parity-flip equivalence (sphere = slab with sign + sin↔cos shuffle).

The X(μ) function (Eq. 40) is **State 1B** — single integral, evaluated
by mpmath quadrature with successive subdivision near the right endpoint
where ``ln(ν - μ)`` is improper. The z_0 extrapolated endpoint (Eq. 42)
is also State 1B. Both reductions are SymPy-symbolic (the integrand is
a SymPy expression); only the quadrature evaluation is numerical.

The validity bound c ≤ 1 + 1/(3 f_1) (Eq. 5) is documented but not
algebraically derivable from Eq. 12 in closed form (it's a transcendental
condition on the dispersion-root reality); the SymPy module
verifies the limit f_1 → 0 reduces to "all c" trivially.
"""
from __future__ import annotations

from .derivations import (
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
