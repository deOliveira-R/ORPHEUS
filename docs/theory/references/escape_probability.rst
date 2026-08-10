.. _theory-escape-probability:

==========================================================================
Escape-probability methods (reserved)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred.** Reserves
   ``:label:`` anchors for the escape-probability reference-solver
   theory. Python folder ``orpheus.derivations.continuous.escape_probability``
   is **reserved, not yet implemented** — a README pinning canonical
   references and no module, so it is named as a literal rather than a
   ``:mod:`` role.


Why this folder exists
======================

Escape-probability methods are the **resonance-self-shielding
ancestor** of collision-probability methods. The Wigner-Seitz cell
construction (cylindrical fuel pin in a moderator-equivalent annulus)
is the classical lattice-physics tool for resonance treatment and
multi-group cross-section condensation in thermal reactors.

Structurally distinct from but mathematically related to
:mod:`orpheus.derivations.continuous.flat_source_cp` (flat-source-per-
region collision probability). The escape-probability formulation
emphasises the *first-flight escape* as the primary geometric
quantity; collision probability emphasises the
*collision-rate matrix* as the primary algebraic object. Both reduce
to the same E_3 / Ki_3 kernel evaluations in homogeneous geometries.


Canonical references
====================

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Carlvik, I. (1967).** "A method for calculating collision
       probabilities in general cylindrical geometry and applications
       to flux distributions and Dancoff factors." *NSE* **30**.
     - ✓
     - Finite-cylinder and cuboid collision probabilities via
       :math:`E_3` kernels.
   * - **Carlvik, I. (1965).** Geneva Conference Vol. 2 p. 225.
     - n/a (NOT local)
     - Foundational paper introducing the escape-probability
       formulation. Cited by Sanchez 1977 Ref. 27.
   * - **Benoist, P. (1981).** "Integral Transport Theory Formalism
       for Diffusion Coefficient Calculations in Wigner-Seitz Cell."
     - ✓
     - Diffusion-coefficient corrections from the escape-probability
       formulation in heterogeneous lattices.
   * - **Stammler, R.J.J. & Abbate, M.J. (1983).** *Methods of Steady
       State Reactor Physics in Nuclear Design.* Chs 4 + 6.
     - ✓
     - Textbook treatment of the Wigner-Seitz cell + lattice
       collision probabilities.
   * - **Hébert, A. (2009).** *Applied Reactor Physics.* Ch. 3.
     - ✓
     - Modern textbook treatment of integral transport, Peierls form,
       and the escape-probability connection.


Cross-references
================

- :ref:`theory-collision-probability` — sister method (flat-source
  CP is the production-grade specialisation).
- :ref:`theory-peierls-nystrom` — production-grade Peierls reference
  consuming the same kernel families.
- ``orpheus.derivations.continuous.escape_probability`` — Python
  package (currently empty; implementation pending).
