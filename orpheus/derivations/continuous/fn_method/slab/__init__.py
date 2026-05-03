"""Slab F_N benchmarks.

Implements the bare-critical-slab F_N method per Siewert-Benoist 1979
(Part I) Section IV + Grandjean-Siewert 1979 (Part II) Section III.

Public entry point:

* :func:`solve_fn_slab_bare_critical` — given the mean number of
  secondaries per collision :math:`c > 1`, returns the critical
  half-thickness :math:`a` (mfp) and the F_N expansion coefficients.
* :func:`fn_slab_flux_at_x_cosine_only` — diagnostic discrete-mode
  flux approximation (NOT the full F_N reconstruction; see docstring).

Verified against:

* LA-13511 ``Ua-1-0-SL`` (Sood problem 12): :math:`r_c = 0.93772556`
  mfp at :math:`c = 1.30` — exact from KLL 1974 Table I.
* Grandjean-Siewert 1979 Table XI: critical thickness for
  :math:`c \\in \\{1.10, 1.30, 1.50, 1.70, 1.90\\}`.

References
----------

* Siewert & Benoist 1979, *Nucl. Sci. Eng.* **69**, 156-160.
* Grandjean & Siewert 1979, *Nucl. Sci. Eng.* **69**, 161-168.
* Kaper, Lindeman & Leaf 1974, *Nucl. Sci. Eng.* **54**, 94.
"""
from __future__ import annotations

from .flux_reconstruction import (
    KLLSlabFluxResult,
    slab_angular_flux_from_scalar,
    slab_scalar_flux_from_angular_quadrature,
    slab_scalar_flux_fn_projection,
    slab_scalar_flux_fn_projection_ratio,
    slab_scalar_flux_kll,
    slab_scalar_flux_ratio,
    slab_surface_angular_flux_fn,
    solve_kll_slab_continuum_coefficient,
)
from .one_group import (
    SlabFNResult,
    fn_slab_flux_at_x_cosine_only,
    solve_fn_slab_bare_critical,
)
from .reflected import (
    SlabReflectedFNResult,
    solve_fn_slab_reflected_critical,
)

__all__ = [
    "SlabFNResult",
    "solve_fn_slab_bare_critical",
    "fn_slab_flux_at_x_cosine_only",
    # Reflected-slab F_N (Neshat-Maiorino 1980).
    "SlabReflectedFNResult",
    "solve_fn_slab_reflected_critical",
    # Rich-machinery flux reconstruction (KLL Fredholm + characteristics).
    "KLLSlabFluxResult",
    "solve_kll_slab_continuum_coefficient",
    "slab_scalar_flux_kll",
    "slab_scalar_flux_ratio",
    "slab_angular_flux_from_scalar",
    "slab_scalar_flux_from_angular_quadrature",
    "slab_surface_angular_flux_fn",
    # Path A.i — F_N projection flux extraction (Peierls iteration).
    "slab_scalar_flux_fn_projection",
    "slab_scalar_flux_fn_projection_ratio",
]
