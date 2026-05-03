r"""F_N method analytical benchmark family — Sood/Forster/Parsons LA-13511 (1999).

This package is a structurally-independent verification ground for
ORPHEUS transport solvers, mounted on the published Sood/Forster/Parsons
analytical benchmark test set [SoodLA13511_1999]_. The benchmark
catalogue contains 75 critical configurations (43 1G, 30 2G, 1 3G, 1 6G;
24 infinite, 24 slabs, 9 cylinders, 14 spheres, 4 ISLC) whose critical
radii, :math:`k_\infty`, and selected scalar-flux ratios were computed
in the peer-reviewed transport-theory literature using **Case singular
eigenfunctions, F_N method, B_N method, and Green's-function method**.

The package follows the algebra-of-record bifurcation discipline:

* :mod:`.origins` — Branch 1 SymPy: closed-form derivations from the
  transport equation. These verify the *algebra* — that the published
  closed forms follow from the reduction chain documented in LA-13511
  Section V + Appendix A.
* :mod:`.multi_group` — Branch 2 numpy/scipy: production code that
  evaluates the closed forms on the catalogued cases. These verify the
  *numerics* — that the algebra evaluates to the published reference
  values.
* :mod:`.benchmarks` — machine-readable catalogue of LA-13511 cases.

Subpackage stubs (:mod:`.core`, :mod:`.slab`, :mod:`.sphere`,
:mod:`.cylinder`) are placeholders for the F_N machinery from
Kaper-Lindeman-Leaf 1974 (slab + sphere) and Westfall-Metcalf 1973
(cylinder). They will be populated when the structurally-independent
F_N reference solvers are added; the LA-13511 report does not contain
the F_N method specification (only the truth set), so populating those
folders requires acquiring the cited journal references.

This first slice ships **only the k_inf cases** (:math:`PUa-1-0-IN`
and :math:`PU-2-0-IN`), which need only LA-13511 Appendix A — pure
rational algebra in cross sections, no special functions. They prove
the Branch-1/Branch-2 pattern works at this complexity, and serve as
the V&V ground floor for the multi-group infrastructure that the F_N
method extensions will need.

References
----------

.. [SoodLA13511_1999] Sood, A., Forster, R.A., Parsons, D.K. (1999).
   "Analytical Benchmark Test Set for Criticality Code Verification."
   Los Alamos National Laboratory report LA-13511. 75 problems, 67
   tables. PNE-2003 is the journal-published condensation.
"""
from __future__ import annotations
