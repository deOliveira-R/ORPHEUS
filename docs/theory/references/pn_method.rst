.. _theory-pn-method:

==========================================================================
P_N method — spherical-harmonics expansion (reserved)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred.** Reserves
   ``:label:`` anchors for the P_N reference-solver theory. Python
   folder: :mod:`orpheus.derivations.continuous.pn_method` (empty;
   README pins canonical references).


Why this folder exists
======================

The P_N method expands the :term:`angular flux` in Legendre / spherical-
harmonic polynomials. It is the historical alternative to the
:term:`discrete-ordinates <ordinate>` S_N family, with a different error structure
(modal truncation vs. ordinate-set discretisation).

Within the ORPHEUS verification stack, P_N is the **natural
structurally-independent cross-check** for multi-region sphere
results from :ref:`theory-trajectory-resolvent`. ORPHEUS already
cross-checks against Garcia 2021 multi-region sphere truth values
inside ``tests/derivations/test_peierls_greens_function_garcia2021.py``;
having ``pn_method/`` reserved makes this asymmetry explicit:
today Garcia 2021 is a *truth set without a method-of-record*, but
once this folder is populated the cross-check becomes
**method-vs-method** instead of **method-vs-tabulated-numbers**.


Canonical references
====================

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Garcia, R.D.M. & Siewert, C.E. (2017).**
     - ✓
     - P_N particular solution for radiative transfer in spherical
       geometry.
   * - **Garcia, R.D.M. et al. (2019).**
     - ✓
     - Spherical-shell P_N.
   * - **Garcia, R.D.M. et al. (2021).** "Accurate spherical
       harmonics solutions ... multi-region spherical geometry."
     - ✓
     - Multi-region sphere P_N for fixed-source. The primary truth
       set ORPHEUS cross-checks against.
   * - **Davison, B. (1957).** *Neutron Transport Theory.* Oxford,
       chs 9–11.
     - n/a
     - Foundational textbook treatment of the P_N method.


Cross-references
================

- :ref:`reference-solvers-three-meanings` — taxonomy positioning.
- :ref:`theory-trajectory-resolvent` — current cross-check target.
- ``orpheus.derivations.continuous.pn_method`` — Python package
  (currently empty; implementation pending).
