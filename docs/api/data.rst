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

The energy-axis value objects for spectrum-weighted energy condensation
(the energy-axis transpose of spatial homogenization), in
:mod:`orpheus.data.energy_grid`:

* :class:`~orpheus.data.energy_grid.EnergyGrid` — a multigroup energy
  structure (strictly descending boundaries, the canonical fast-first
  convention; the energy analogue of a coarse
  :class:`~orpheus.geometry.mesh.Mesh1D`).
* :class:`~orpheus.data.energy_grid.GroupCondensation` — the fine→coarse
  partition map, built by conservative fractional re-binning (a
  partition-of-unity overlap table; nested → one-hot degenerate). Carries
  the upscaling guard and the
  :attr:`~orpheus.data.energy_grid.GroupCondensation.locally_interpolated`
  provenance report.
* :class:`~orpheus.data.energy_grid.WithinGroupSpectrum` — the selectable
  within-group flux model :math:`w(E)` for straddle apportionment, with
  :class:`~orpheus.data.energy_grid.InverseEnergySpectrum` (1/E, the
  default) as the first realisation.

See :ref:`sn-energy-condensation` for the theory: the rate-preserving
collapse, the fractional-overlap re-binning that lifts the nesting
requirement, the Petrov-Galerkin discipline (flux = test weight), and the
downsampling-only semantic. The per-material collapse verb is
:meth:`Mixture.condense <orpheus.data.macro_xs.mixture.Mixture.condense>`
(above); the orchestration is :meth:`Solution.condense
<orpheus.sn.solution.Solution.condense>`.

.. note::

   Per-symbol autodoc for :mod:`orpheus.data.energy_grid` is pending a
   one-line docstring fix in ``GroupCondensation.from_grids`` (a
   ``:class:`EnergyGrid`s`` role with a trailing ``s`` that renders as an
   unterminated inline-literal under ``-W``). Once corrected to
   ``:class:`EnergyGrid`\\ s`` (or rephrased), an ``.. automodule::
   orpheus.data.energy_grid`` block surfaces the full per-symbol API
   here.

.. note::

   Verification benchmarks have moved to the :doc:`derivations <derivations>`
   package. See :doc:`../theory/verification` for the full reference case table.
