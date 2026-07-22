.. _theory-galerkin-sn-hybrid:

==========================================================================
Hybrid collocation-Galerkin-S_N (reserved, lower priority)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred (lower
   priority).** Reserves ``:label:`` anchors for the hybrid
   Galerkin-S_N reference-solver theory. Python folder:
   ``orpheus.derivations.continuous.galerkin_sn_hybrid`` (empty;
   README pins canonical references).


Why this folder exists
======================

The hybrid collocation-Galerkin-S_N angular discretisation
(Morel 1989) sits between modal expansion (P_N) and discrete-
ordinates (S_N): it keeps discrete-ordinate :term:`sweeps <sweep>` as the spatial
mechanism while applying Galerkin projection in the angular variable.

The hybrid is useful for problems where neither P_N nor S_N alone
behaves well — typically highly anisotropic scattering combined with
strong streaming, where ray effects in S_N and ringing in P_N both
bite the converged solution.

This is research-grade with a narrow application surface, hence the
lower-priority placeholder status.


Canonical references
====================

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Morel, J.E. (1989).** "A Hybrid Collocation-Galerkin-Sn
       Method for Solving the Boltzmann Transport Equation."
     - ✓
     - Foundational hybrid construction.


Cross-references
================

- :ref:`theory-pn-method` — modal expansion route (one parent).
- :ref:`theory-discrete-ordinates` — discrete-ordinates route
  (other parent).
- ``orpheus.derivations.continuous.galerkin_sn_hybrid`` — Python
  package (currently empty; implementation pending).
