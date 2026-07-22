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

Cross sections follow the same priority as the :term:`scalar flux`:
:math:`g` first, then spatial.  The canonical per-cell XS state is
**one attribute** --- ``SNSolver.mat_xs``, a
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
wrapping the per-material
:class:`~orpheus.data.macro_xs.mixture.Mixture` data together with
the per-cell typed views (Issue #197 PR-TYPED-1; the earlier
per-attribute surface ``sig_t`` / ``sig_a`` / ``sig_p`` / ``chi``
was retired by PR-TYPED-2).  Every operator (:math:`L, C, S, F`)
reads cross sections through this single source; its accessors
(e.g. ``mat_xs.total_cross_section``) yield the principled
``(ng, nx, ny)`` layout.

The producer :func:`~orpheus.data.macro_xs.assemble_cell_xs`
emits the flat ``(N_cells, ng)`` shape (CP also consumes that flat
shape --- the producer is *unchanged* by the SN migration).  The
:ref:`cell-flattening invariant <sn-cell-flattening-invariant>` is
asserted in ``__debug__`` builds at construction, exercised through
the ``mat_xs.total_cross_section`` accessor against the producer's
legacy flat shape.


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

:attr:`Mixture.SigT`, :attr:`Mixture.SigP`, etc. are stored
per-mixture as shape ``(ng,)`` (group-only).  The per-cell views on
``SNSolver.mat_xs`` are shape ``(ng, nx, ny)`` --- group first, then
spatial.  The bridge is
:func:`~orpheus.data.macro_xs.assemble_cell_xs`, which lifts the
per-mixture group-only array to the per-cell flat shape
``(N_cells, ng)``; the
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
accessors expose the principled per-cell ``(ng, nx, ny)``.  CP
consumes the producer's flat shape directly --- CP is unaffected by
the SN migration.
