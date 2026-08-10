.. _theory-spectral-collocation:

==========================================================================
Spectral collocation — Chebyshev / Gegenbauer / Jacobi (reserved)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred (no concrete
   reference selected yet).** Reserves ``:label:`` anchors for a
   future spectral-collocation reference-solver theory. Python folder
   ``orpheus.derivations.continuous.spectral_collocation`` is **reserved,
   not yet implemented** — a README pinning candidate references and no
   module, so it is named as a literal rather than a ``:mod:`` role.


Why this folder exists
======================

Pure spectral-collocation methods (Chebyshev, Gegenbauer, Jacobi
polynomial families) for the transport equation. Reserved for
**taxonomic completeness** — when a concrete reference is selected
and a benchmark identified, this page will be upgraded with the
canonical-reference table and the implementation will be queued.

Spectral collocation is most useful as a Branch-1 (symbolic +
spectrally-converged numerical) verification reference rather than
a Branch-2 production solver. The high spectral accuracy is wasted
on production geometries with discontinuous material interfaces;
its sweet spot is the smooth-coefficient verification problem where
exponential convergence yields machine-precision references at
modest expansion order.


Candidate references
====================

(Not all locally available — concrete selection pending.)

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Garcia, R.D.M., Siewert, C.E. (1979 / 1980s).** Chebyshev /
       Gegenbauer collocation work for transport theory.
     - partial
     - Some of this family lives on HAL; not all PDFs are local.
   * - **Hauser, M., Pomraning, G.C. (1980s).** Spectral collocation
       for transport problems with high-anisotropy scattering.
     - n/a
     - Candidate reference; needs concrete selection.


Cross-references
================

- :ref:`theory-galerkin-spectral` — Galerkin-spectral sibling (a
  closely related, already-implemented method).
- :ref:`theory-fn-method` — F_N method as a sibling collocation
  scheme.
- ``orpheus.derivations.continuous.spectral_collocation`` — Python
  package (currently empty; implementation pending).
