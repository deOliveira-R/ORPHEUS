.. _theory-bn-method:

==========================================================================
B_N method — boundary-collocated F_N variant (reserved, lower priority)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred (lower
   priority).** Reserves ``:label:`` anchors for the B_N
   reference-solver theory. Python folder:
   :mod:`orpheus.derivations.continuous.bn_method` (empty; README
   pins canonical references).


Why this folder exists
======================

The B_N method is the **boundary-collocated sister** of the F_N
method (:ref:`theory-fn-method`): it collocates the F_N equations
on the boundary points instead of interior collocation points. Used
historically for diffuse-reflection benchmarks and for problems where
the boundary-layer structure dominates the angular response.

In the ORPHEUS verification stack, F_N already covers the slab and
sphere benchmark surface, so B_N is a niche addition rather than a
foundational pillar. It is reserved for taxonomic completeness and as
a possible future extension when diffuse-reflection benchmarks
require additional cross-checking.


Canonical references
====================

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Brockmann, R. (1981).** B_N method for anisotropic-scattering
       treatment.
     - ✓
     - Foundational treatment of the B_N collocation scheme.
   * - **Sood, A., Forster, R.A., Parsons, D.K. (2003).** LA-13511
       benchmark catalogue.
     - ✓
     - Cites the B_N method in the diffuse-reflection benchmark
       context.
   * - **Garibba & Rojas (1980s).** Technical reports.
     - n/a (NOT local)
     - Original diffuse-reflection benchmark formulations cited by
       LA-13511.


Cross-references
================

- :ref:`theory-fn-method` — sister method (interior collocation).
- :ref:`theory-singular-eigenfunction` — underlying singular-
  eigenfunction expansion both methods rest on.
- ``orpheus.derivations.continuous.bn_method`` — Python package
  (currently empty; implementation pending).
