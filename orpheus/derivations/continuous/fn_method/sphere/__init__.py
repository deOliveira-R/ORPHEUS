r"""Sphere bare-critical 1G F_N solver.

Public entry point:

* :func:`solve_fn_sphere_bare_critical` — given the mean number of
  secondaries per collision :math:`c > 1`, returns the bare-critical
  sphere radius :math:`R_c` (mfp) via the F_N method (Siewert-Thomas
  1986).

The sphere F_N is the slab F_N with one sign flip on the boundary
attenuation block (Siewert-Thomas Eq. 46 vs Siewert-Benoist Eq. 4).
Both share the same moment-integral recursion, collocation grid,
dispersion-root finder, and X-function machinery — the literature-
researcher-validated structural fact verified symbolically in
:mod:`..origins.fn_sphere_derivations`.

References
----------

* Siewert & Thomas 1986, *Nucl. Sci. Eng.* **94**, 264.
* Sood, Forster & Parsons 1999, LA-13511 Table 6 (case ``Ua-1-0-SP``).
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94 (Table V).
"""
from __future__ import annotations

from .flux_reconstruction import (
    KLLSphereFluxResult,
    solve_kll_sphere_continuum_coefficient,
    sphere_angular_flux_from_scalar,
    sphere_scalar_flux_from_angular_quadrature,
    sphere_scalar_flux_kll,
    sphere_scalar_flux_ratio,
    sphere_surface_angular_flux_fn,
)
from .one_group import (
    SphereFNResult,
    solve_fn_sphere_bare_critical,
)

__all__ = [
    "SphereFNResult",
    "solve_fn_sphere_bare_critical",
    # Rich-machinery flux reconstruction (KLL Fredholm + characteristics).
    "KLLSphereFluxResult",
    "solve_kll_sphere_continuum_coefficient",
    "sphere_scalar_flux_kll",
    "sphere_scalar_flux_ratio",
    "sphere_angular_flux_from_scalar",
    "sphere_scalar_flux_from_angular_quadrature",
    "sphere_surface_angular_flux_fn",
]
