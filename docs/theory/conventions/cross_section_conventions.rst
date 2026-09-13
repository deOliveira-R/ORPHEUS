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
**one field per Problem** --- ``mat_xs`` on the carrier hub
(:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` and its
method subclasses), a
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
wrapping the per-material
:class:`~orpheus.data.macro_xs.mixture.Mixture` data together with
the per-cell typed views (Issue #197 PR-TYPED-1; the earlier
per-attribute surface ``sig_t`` / ``sig_a`` / ``sig_p`` / ``chi``
was retired by PR-TYPED-2).  Every operator (:math:`L, C, S, F`)
reads cross sections through this single source; its accessors
(e.g. ``mat_xs.total_cross_section``) yield the principled
``(ng, nx, ny)`` layout.

.. note:: **It was a solver attribute until 2026-09-13.**

   The solver used to build its own field in ``__init__`` and hold it as
   an attribute, and every operator read ``self.mat_xs``.  At the
   consumers campaign's step 2 the field became a
   :func:`~functools.cached_property` of the **hub** — one
   :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField` per
   Problem, shared by every consumer of that Problem — and
   ``SNSolver.mat_xs`` was **deleted**; the solver reads
   ``self.sn_mesh.mat_xs``.  The move is what makes
   :math:`\sigma_t` a Problem *datum* rather than solver state: see
   :ref:`sn-sigma-is-a-problem-datum`.  (``DiffusionSolver.mat_xs``
   survives as a solver-side *read* of the same hub property — a
   homonym, not the retired attribute.)

The producer :func:`~orpheus.data.macro_xs.cell_xs.assemble_cell_xs`
emits the flat ``(N_cells, ng)`` shape (CP also consumes that flat
shape --- the producer is *unchanged* by the SN migration).  The
:ref:`cell-flattening invariant <sn-cell-flattening-invariant>` is the
contract between that flat shape and the principled one.  It was an
``assert`` inside ``if __debug__:`` at solver construction until
2026-09-13 --- inert under the canonical ``python -O`` runner, which
strips ``assert`` statements --- and is now carried by two
``foundation`` gates, one of them the
``@pytest.mark.verifies("sn-cell-flatten-roundtrip")`` witness on the
hub's own :math:`\sigma_t` datum.


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
per-mixture as shape ``(ng,)`` (group-only).  The per-cell views on the
hub's ``mat_xs`` are shape ``(ng, nx, ny)`` --- group first, then
spatial.  The bridge is
:func:`~orpheus.data.macro_xs.cell_xs.assemble_cell_xs`, which lifts the
per-mixture group-only array to the per-cell flat shape
``(N_cells, ng)``; the
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`
accessors expose the principled per-cell ``(ng, nx, ny)``.  CP
consumes the producer's flat shape directly --- CP is unaffected by
the SN migration.

One of the four per-cell views is **not** a gather at read time.  Since
2026-09-13 :math:`\sigma_t` is stored on the hub as
``sigma_t_cell``, assembled once at construction in exactly the layout
above, and ``mat_xs.total_cross_section`` *reads that datum*; the other
three (:math:`\sigma_a`, :math:`\nu\Sigma_f`, :math:`\chi`) are still
gathered lazily from the materials.  The asymmetry is deliberate and is
the whole content of :ref:`sn-sigma-is-a-problem-datum` --- a
:math:`\sigma_t` that is not a pure function of the materials is what
makes a depletion or thermal-feedback step expressible as *another
Problem* instead of a mutation of a live solver.
