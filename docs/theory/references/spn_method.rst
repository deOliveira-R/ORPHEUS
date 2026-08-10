.. _theory-spn-method:

==========================================================================
SP_N method — Simplified P_N (reserved)
==========================================================================

.. note::

   **Stub-grade theory page — implementation deferred.** Reserves
   ``:label:`` anchors for the SP_N reference-solver theory. Python
   folder ``orpheus.derivations.continuous.spn_method`` is **reserved, not
   yet implemented** — a README pinning canonical references and no
   module, so it is named as a literal rather than a ``:mod:`` role.


Why this folder exists
======================

The Simplified P_N (SP_N) method is the asymptotic reduction of the
full P_N equations in slab geometry, generalised to multi-D. SP_N
occupies the **bridge between diffusion and full transport** — its
operator is a small extension of the diffusion equation but its
accuracy is closer to true transport.

For ORPHEUS this is pedagogically valuable (the "teaching tool"
mission emphasises the spectrum from diffusion to transport) and
practically useful as an intermediate-fidelity reference for the
production diffusion solver.


Canonical references
====================

.. list-table::
   :header-rows: 1
   :widths: 32 18 50

   * - Reference
     - Local PDF
     - Provides
   * - **Makine, J. et al. (2018).** "Exact transport representations
       of the classical and nonclassical simplified P_N equations."
     - ✓
     - The headline reference for this folder. Exact transport
       representations of SP_N derived rigorously.
   * - **Brantley, P.S., Larsen, E.W. (2000).** "The Simplified P_3
       Approximation." *NSE* **134**, 1–21.
     - n/a
     - Foundational SP_3 formulation.


Cross-references
================

- :ref:`theory-diffusion-1d` — the operator SP_N reduces to in the
  diffusion limit.
- :ref:`theory-pn-method` — the parent operator SP_N is reduced from.
- ``orpheus.derivations.continuous.spn_method`` — Python package
  (currently empty; implementation pending).
