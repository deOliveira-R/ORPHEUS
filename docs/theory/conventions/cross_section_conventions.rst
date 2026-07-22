.. _theory-conventions-cross-sections:

=========================
Cross-Section Conventions
=========================

**Which way the scattering arrow points, and where cross sections live.**
This page carries the SN cross-section conventions --- the
``SigS[g_from, g_to]`` storage law with its transpose acting form
``Q_scatter = SigS^T @ phi`` and the per-material versus per-cell split.
It *links* rather than restates the two single sources of truth: rows
1--2 of :ref:`notation-crosswalk` fix the arrow direction, and
:eq:`sigs-convention` is the governing storage identity. The array
layout beneath is the sibling :doc:`indexing_and_layout` page.

Cross-section convention
========================

Cross sections follow the same priority as the scalar flux:
:math:`g` first, then spatial.  Per-cell cross sections are stored as
``(ng, nx, ny)`` numpy arrays on
:class:`~orpheus.sn.solver.SNSolver`:

.. code-block:: python

   class SNSolver:
       sig_t: np.ndarray   # (ng, nx, ny) total
       sig_a: np.ndarray   # (ng, nx, ny) absorption
       sig_p: np.ndarray   # (ng, nx, ny) production (νΣ_f)
       chi:   np.ndarray   # (ng, nx, ny) fission spectrum

The producer :func:`~orpheus.data.macro_xs.assemble_cell_xs`
emits the flat ``(N_cells, ng)`` shape (CP also consumes that flat
shape --- the producer is *unchanged* by the SN migration).  The
``.T.reshape(ng, nx, ny)`` bridge lives at exactly one site,
:meth:`SNSolver.__init__`, and the
:ref:`cell-flattening invariant <sn-cell-flattening-invariant>`
asserts the round-trip in ``__debug__`` builds.


SigS scattering convention --- still ``[g_from, g_to]``
=======================================================

The :ref:`scattering-matrix-convention` is unchanged by the migration.
:attr:`Mixture.SigS` matrices are stored as
``SigS[l][g_from, g_to]``; the in-scatter source uses the transpose:
``Q_scatter = SigS^T @ phi``.

The layout migration affects the **storage** of the resulting flux
arrays, not the **convention** of the cross-section matrices.


Per-material vs per-cell cross sections
=======================================

:attr:`Mixture.SigT`, :attr:`Mixture.NuSigF`, etc. are stored per-mixture
as shape ``(ng,)`` (group-only).  The per-cell arrays on
:class:`SNSolver` (``sig_t``, ``sig_a``, ``sig_p``, ``chi``) are
shape ``(ng, nx, ny)`` --- group first, then spatial.  The bridge is
:func:`~orpheus.data.macro_xs.assemble_cell_xs`, which lifts the
per-mixture group-only array to the per-cell flat shape
``(N_cells, ng)``; :meth:`SNSolver.__init__` then transposes and
reshapes to the principled per-cell ``(ng, nx, ny)``.  CP consumes
the producer's flat shape directly --- CP is unaffected by the SN
migration.
