Data Package (``data``)
=======================

Mixture
-------

.. automodule:: orpheus.data.macro_xs.mixture
   :members:
   :undoc-members:
   :show-inheritance:

Cell Cross Sections
-------------------

.. automodule:: orpheus.data.macro_xs.cell_xs
   :members:
   :undoc-members:
   :show-inheritance:

Energy-Group Structure and Condensation
---------------------------------------

The energy-axis value object for spectrum-weighted energy condensation
(the energy-axis transpose of spatial homogenization). An
:class:`~orpheus.data.energy_grid.EnergyGrid` is the energy analogue of a
coarse :class:`~orpheus.geometry.mesh.Mesh1D`: a strictly descending
boundary array (the canonical fast-first convention) that yields **both**
halves of a discrete frame — :meth:`~orpheus.data.energy_grid.EnergyGrid.as_measure`
(the source view, project *from*) and
:meth:`~orpheus.data.energy_grid.EnergyGrid.as_basis` (the nested target
view, project *to*) — plus the irreducibly binary
:meth:`~orpheus.data.energy_grid.EnergyGrid.overlap_to` mismatch factory
(the non-nested fractional :class:`~orpheus.numerics.basis.OverlapBasis`).
The condensation is **data-native**: the grid, the overlap factory, and
the Petrov-Galerkin frame all live in ``data`` / ``numerics``, with no
transport dependency. There is no ``GroupCondensation`` /
``CondensationFrame`` type — the collapse identity lives in
:meth:`~orpheus.data.energy_grid.EnergyGrid.overlap_to` (the trial) and
the :meth:`Mixture.condense <orpheus.data.macro_xs.mixture.Mixture.condense>`
verb (above). See :ref:`sn-condensation-no-frame-subclass` for that
rationale.

.. automodule:: orpheus.data.energy_grid
   :members:
   :undoc-members:
   :show-inheritance:

See :ref:`sn-energy-condensation` for the theory: the rate-preserving
collapse, the fractional-overlap re-binning that lifts the nesting
requirement, the Petrov-Galerkin discipline (flux = test weight), and the
downsampling-only semantic. The per-material collapse verb is
:meth:`Mixture.condense <orpheus.data.macro_xs.mixture.Mixture.condense>`
(above); the orchestration is :meth:`Solution.condense
<orpheus.sn.solution.Solution.condense>`.

Energy-group structures
~~~~~~~~~~~~~~~~~~~~~~~~~

The named :class:`~orpheus.data.energy_grid.EnergyGrid` instances (the
*data*; the type + within-group flux models stay in
:mod:`orpheus.data.energy_grid` above), in
:mod:`orpheus.data.group_structures`:

* :data:`~orpheus.data.group_structures.ornl.ORNL_421` — the canonical
  421-group fine grid the production library is generated on (the
  condensation source; GENERATED from the live library ``eg``,
  L4-self-validated against ``pwr_like_mix().eg``, superseding the
  retired ``EGB421.txt``).
* ``orpheus.data.group_structures.wims`` — the WIMS-D/IAEA few-group
  condensation targets ``WIMS_69`` / ``WIMS_172`` and the Table 11.3
  derivation-validation oracle ``CONDENSE_172_TO_69`` (IAEA-TECDOC /
  STI/Pub/1264, *WIMS-D Library Update*, Tables 11.1–11.3).

.. automodule:: orpheus.data.group_structures.ornl
   :members:

.. note::

   ``orpheus.data.group_structures.wims`` is cited in prose (not
   ``automodule``\ d) because its module docstring carries a plain-text
   (non-RST) transcription-defect note whose 2-space-indented numbered
   list trips ``-W`` (``ERROR: Unexpected indentation`` /
   ``WARNING: Block quote ends without a blank line``). Re-flowing that
   note to valid RST (blank line before the numbered list; consistent
   indentation) unblocks an ``.. automodule:: orpheus.data.group_structures.wims``
   block here — a docs-only follow-up that does not touch the grid data.

.. note::

   Verification benchmarks have moved to the :doc:`derivations <derivations>`
   package. See :doc:`../theory/verification` for the full reference case table.
