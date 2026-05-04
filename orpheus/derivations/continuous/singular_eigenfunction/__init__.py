r"""Case singular eigenfunction expansion family.

This package collects ORPHEUS Branch-2 production solvers built on
**Case's singular eigenfunction expansion** [Case1960]_ — the
foundational analytic technique for the one-speed (1G) integral
transport equation. The package is organized by geometry and
scattering anisotropy:

* :mod:`.slab` — Atalay 1997 reflected-slab criticality with
  **linearly anisotropic scattering** (P_1 expansion).
* :mod:`.sphere` — Atalay 1997 reflected-sphere criticality (linearly
  anisotropic) via the slab-odd-mode parity flip.
* :mod:`.cylinder` — Westfall-Metcalf 1973 bare radially-reflected
  cylinder, **isotropic scattering**, modified-Bessel-K radial kernels.
* :mod:`.anisotropy` — linearly-anisotropic scattering primitives
  (validity bounds, half-range relation specialisations).
* :mod:`.core` — shared above-trusted-library primitives (Atalay
  X-function, dispersion roots, half-range moment integrals).
* :mod:`.origins` — Branch-1 SymPy algebra-of-record (slab-sphere
  Atalay derivations, cylinder Westfall-Metcalf derivations).

Why this consolidation
======================

This package merges what were historically separate ``case_method/``
(Atalay 1997 slab + sphere) and ``singular_eigenfunction/cylinder/``
(Westfall-Metcalf 1973) folders. The merge is justified because **all
primary authors call this technique "singular eigenfunction expansion"**,
not "Case's method" — and they all explicitly cite Case 1960 as the
foundational source:

* **Westfall-Metcalf 1973** (introduction): "*Since the introduction of
  the singular eigenfunction expansion technique by Case [Ref. 1] in
  1960, a wide variety of transport problems have been treated by this
  method.*"
* **Mitsis 1963** (abstract): "*Transport solutions ... developed by
  the method of singular expansion modes.*"
* **Atalay 1997** (abstract): "*Case's singular eigenfunction method
  is used to formulate the criticality conditions.*"

The two ORPHEUS folders represented **two parametric variations on the
same method**, not two methods. Different geometry (slab/sphere vs
cylinder), different scattering anisotropy (linear vs isotropic),
different kernel reduction (exponentials vs Bessel-K via the addition
theorem), but **identical mathematical machinery** above the trusted-
library line: discrete eigenvalue :math:`\nu_0` (same dispersion
function), continuum modes on :math:`(-1, 1)`, half-range
orthogonality, Fredholm reduction of the boundary condition. The
folder rename also cures the previous "case_method" folder's
violation of the no-author-names project rule.

Distinction from F_N (mod:`...fn_method`)
==========================================

Both packages solve the same physical equation via singular-eigen-
function technology — but the **structural mathematics is different**:

* :mod:`...fn_method` (Siewert et al. 1979-1986) builds an
  :math:`(N+1) \times (N+1)` complex matrix :math:`M` from
  Case-eigenfunction collocation against a Wiener-Hopf
  X-function, then solves :math:`\det M = 0` for the critical
  dimension. Slab and sphere are unified by a geometry sign
  :math:`s = \pm 1`.

* This package follows the **Mitsis 1963 path**: direct singular
  eigenfunction expansion with Atalay's first-order Fredholm
  iteration for slab/sphere with linear anisotropy, and the
  Westfall-Metcalf Bessel-K cylinder reduction for the radially
  reflected isotropic cylinder.

Sharing the same dispersion function
:math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu)` is a *medium*
property — both packages reuse
:func:`...fn_method.core.dispersion.case_nu0` for the dominant
Case eigenvalue. Below the trusted-library line, primitive sharing
is fine (per ``algebra-of-record`` § "Structural independence
applies above the trusted-library line"). The structural
difference lies in everything *above* that line: collocation
grid, matrix assembly, root-finding strategy.

Validity range — slab + sphere (Atalay)
========================================

The bi-orthogonality relations are valid for :math:`c \le 1 + 1/(3 f_1)`
(Atalay Eq. 5), the range where the transport operator has only one
pair of discrete modes :math:`\pm \nu_0`. For :math:`f_1 = 0` (isotropic
limit) the bound is trivial (all :math:`c`); for :math:`f_1 = 0.30` the
upper limit is :math:`c \le 2.111`. Outside this band complex
eigenvalues appear (Dahl-Sjöstrand 1979, Kohut 1993) — first-order
Atalay does not detect them.

Discipline
==========

* :mod:`.origins` — **Branch-1 SymPy** algebra-of-record:
  closed-form derivations + mpmath integral evaluations for cases
  SymPy chokes on.
* :mod:`.slab`, :mod:`.sphere`, :mod:`.cylinder` — **Branch-2**
  numpy/scipy/mpmath production solvers.
* :mod:`.core` — primitives shared across geometries within this
  package family. Reuses :func:`...fn_method.core.dispersion.case_nu0`
  (trusted-library line) but the criticality search, matrix
  machinery, and flux reconstruction are independent.

The Sphinx narrative lives in :doc:`/theory/singular_eigenfunction`.

References
----------

.. [Case1960] Case, K. M. (1960). "Elementary Solutions of the
   Transport Equation and Their Applications." *Annals of Physics*
   **9**, 1-23.

.. [Mitsis1963] Mitsis, G. F. (1963). "Transport Solutions to the
   Monoenergetic Critical Problems." Argonne National Laboratory
   report ANL-6787.

.. [WestfallMetcalf1973] Westfall, R. M. & Metcalf, D. R. (1973).
   "Singular Eigenfunction Solution of the Monoenergetic Neutron
   Transport Equation for Finite Radially Reflected Critical
   Cylinders." *Nuclear Science and Engineering* **52**, 1-11.
   DOI 10.13182/NSE73-A23285.

.. [Atalay1997] Atalay, M.A. (1997).
   "The reflected slab and sphere criticality problem with anisotropic
   scattering in one-speed neutron transport theory."
   *Progress in Nuclear Energy* **31**(3), 229-252.
   DOI: 10.1016/0149-1970(95)00094-1.

.. [BurkartIshiguroSiewert1976] Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976).
   "Neutron transport in two dissimilar media with anisotropic scattering."
   *Nuclear Science and Engineering* **61**, 72-81.
   (Used as cross-check for Atalay vacuum-slab + anisotropic results.)

.. [McCormickKuscer1965] McCormick, N.J., Kušcer, I. (1965).
   "Bi-orthogonality relations for solving half-space transport problems."
   *J. Math. Phys.* **6**, 1939.

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
from .slab import (
    CaseMethodSlabResult,
    solve_case_method_slab_critical,
)
from .sphere import (
    CaseMethodSphereResult,
    solve_case_method_sphere_critical,
)
from .spectrum import Spectrum

__all__ = [
    # cylinder (Westfall-Metcalf 1973)
    "CylinderSingularEigenfunctionResult",
    "solve_singular_eigenfunction_cylinder_bare_critical",
    # slab (Atalay 1997)
    "CaseMethodSlabResult",
    "solve_case_method_slab_critical",
    # sphere (Atalay 1997)
    "CaseMethodSphereResult",
    "solve_case_method_sphere_critical",
    # math-heart class (3rd instance of the pattern)
    "Spectrum",
]
