r"""Singular eigenfunction expansion (Mitsis-Westfall-Metcalf family).

This package collects ORPHEUS Branch-2 production solvers built on the
**Mitsis-style direct singular-eigenfunction expansion** of the
monoenergetic 1-G isotropic-scattering integral transport equation
(Eq. 1 of Westfall & Metcalf 1973):

.. math::

   \rho(\vec r) = \frac{1}{4\pi} \int_{\vec r'}
   \frac{c(\vec r')\,\rho(\vec r')\,e^{-|\vec r - \vec r'|}}
        {|\vec r - \vec r'|^2}\,d^3 r' .

Why a *separate* package from :mod:`...fn_method`
==================================================

Both packages solve the same physical equation via singular-eigen-
function technology — but the **structural mathematics is different**:

* :mod:`...fn_method` (Siewert et al. 1979-1986) builds an
  :math:`(N+1) \times (N+1)` complex matrix :math:`M` from
  Case-eigenfunction collocation against a Wiener-Hopf
  X-function, then solves :math:`\det M = 0` for the critical
  dimension. Slab and sphere are unified by a geometry sign
  :math:`s = \pm 1`.

* :mod:`.cylinder` (Westfall-Metcalf 1973) follows the
  **Mitsis 1963 path**: the cylindrical kernel reduces to modified
  Bessel functions :math:`K_0(t/\mu)\,I_0(r/\mu)` via the addition
  theorem (Eq. 5-6), the singular eigenfunctions take the same
  Case form (Eq. 17 + 19), but the radial expansion
  :math:`R_{l\nu}(r) = \alpha_l I_0(r/\nu) + \beta_l K_0(r/\nu)`
  (Eq. 13) is **Bessel-rather-than-exponential**. There is no
  Wiener-Hopf X-function in this formulation. The criticality
  condition reduces to a **Fredholm integral equation** in the
  half-range moment :math:`A'(\nu)` (Eq. 30 + 32), not a finite
  determinantal equation.

Sharing the same dispersion function
:math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu)` (Eq. 18)
is a *medium* property — both packages reuse
:func:`...fn_method.core.dispersion.case_nu0` for the dominant
Case eigenvalue. Below the trusted-library line, primitive sharing
is fine (per ``algebra-of-record`` § "Structural independence
applies above the trusted-library line"). The structural
difference lies in everything *above* that line: collocation
grid, matrix assembly, root-finding strategy.

Discipline
==========

* :mod:`.origins` — **Branch-1 SymPy** algebra-of-record:
  state-1A closed-form derivations of the dispersion function,
  pseudo-eigenfunction structure, and Bessel-Wronskian identities;
  state-1B mpmath integral evaluations where SymPy chokes (e.g.
  the Fredholm sequence for :math:`A'(\nu)`).

* :mod:`.cylinder` — **Branch-2** numpy/scipy/mpmath production
  solver for the bare-critical infinite cylinder. Reuses the
  shared dispersion-root primitive from :mod:`...fn_method.core.dispersion`
  (trusted-library line, see above) but the criticality search,
  matrix machinery, and flux reconstruction are independent.

* :mod:`.core` — primitives shared across geometries within this
  package family. Currently empty (cylinder is the only
  inhabitant); slab/sphere via Mitsis 1963 / Smith 1969 may
  populate later if a benchmark need surfaces.

The Sphinx narrative lives in :doc:`/theory/singular_eigenfunction`.
The current page is **stub-grade**; the archivist agent expands
each TODO marker into the full derivation + numerical evidence
+ literature backstory after method-implementer ships the
prototype.

References
----------

.. [WestfallMetcalf1973] Westfall, R. M. & Metcalf, D. R. (1973).
   "Singular Eigenfunction Solution of the Monoenergetic Neutron
   Transport Equation for Finite Radially Reflected Critical
   Cylinders." *Nuclear Science and Engineering* **52**, 1-11.
   DOI 10.13182/NSE73-A23285.

.. [Mitsis1963] Mitsis, G. F. (1963). "Transport Solutions to the
   Monoenergetic Critical Problems." Argonne National Laboratory
   report ANL-6787.

.. [Case1960] Case, K. M. (1960). "Elementary Solutions of the
   Transport Equation and Their Applications." *Annals of Physics*
   **9**, 1-23.

.. [Sood1999] Sood, A., Forster, R. A., Parsons, D. K. (1999).
   "Analytical Benchmark Test Set for Criticality Code Verification."
   Los Alamos National Laboratory report LA-13511, Table 13
   (case ``Ua-1-0-CY``).
"""
from __future__ import annotations

from .cylinder import (
    CylinderSingularEigenfunctionResult,
    solve_singular_eigenfunction_cylinder_bare_critical,
)

__all__ = [
    "CylinderSingularEigenfunctionResult",
    "solve_singular_eigenfunction_cylinder_bare_critical",
]
