Method of Characteristics Solvers
==================================

The :mod:`orpheus.moc` package provides a 2-D Method of
Characteristics (MOC) solver for PWR-style pin cells. MOC solves
the same multi-group transport equation as SN but trades the
discrete ordinates representation for a ray-traced characteristic
sweep: neutrons are propagated along straight-line tracks through
a flat-source region mesh, and the source iteration converges the
in-scatter and fission contributions.

.. contents::
   :local:
   :depth: 2

.. seealso::

   :ref:`theory-method-of-characteristics` — derivation of the
   characteristic form, the exponential attenuation kernel, the
   track / segment data structures, and the MOC approximations
   (flat source, isotropic scattering, cyclic tracks).


Package Layout
--------------

The MOC implementation is split into four submodules with clear
responsibilities:

* :mod:`orpheus.moc.geometry` — track generation and segment
  data structures. Defines :class:`~orpheus.moc.geometry.Track`
  (a single ray and its list of
  :class:`~orpheus.moc.geometry.Segment` intersections with the
  flat-source region mesh) and :class:`~orpheus.moc.geometry.MOCMesh`
  (the full set of tracks with azimuthal and polar weights).
  Ray tracing routines
  (:func:`~orpheus.moc.geometry._ray_circle_intersections`,
  :func:`~orpheus.moc.geometry._ray_box_intersections`,
  :func:`~orpheus.moc.geometry._trace_single_ray`) live here.
* :mod:`orpheus.moc.quadrature` — 2-D angular quadrature.
  :class:`~orpheus.moc.quadrature.MOCQuadrature` combines an
  azimuthal quadrature (equispaced in
  :math:`[0, 2\pi)` with cyclic track spacing) with a polar
  quadrature (Tabuchi–Yamamoto or equivalent) and normalises
  track weights so that the integrated flat source recovers the
  cell volume exactly.
* :mod:`orpheus.moc.core` — :class:`~orpheus.moc.core.MOCSolver`,
  the :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver`
  protocol implementation. Holds the flux arrays, the
  exponential look-up table, and the characteristic sweep
  routine.
* :mod:`orpheus.moc.solver` — the convenience entry point
  :func:`~orpheus.moc.solver.solve_moc` and the
  :class:`~orpheus.moc.solver.MoCResult` container.


Geometry Note: Cyclic Tracks
----------------------------

MOC relies on **cyclic** (also called *modular* or *periodic*)
track spacing: the track that exits the right boundary re-enters
the left boundary at the mirror-image offset, and the track that
exits the top re-enters the bottom. Cyclic tracks let the
reflective boundary condition be implemented by a simple index
permutation on the angular flux array, with no interpolation —
which is the trick that keeps the characteristic sweep linear in
the number of segments.

The shipped implementation reaches that permutation from the *ray
spacing* side rather than the angle side.
:meth:`MOCQuadrature.create
<orpheus.moc.quadrature.MOCQuadrature.create>` lays down
:attr:`~orpheus.moc.quadrature.MOCQuadrature.n_azi` **unadjusted**
midpoint-even azimuthal angles on :math:`[0, \pi)`,
:math:`\varphi_a = \pi(a + \tfrac12)/n_{\rm azi}` — a set closed
under :math:`\varphi \mapsto \pi - \varphi`, which is what makes the
reflected direction always exist in the quadrature. Per angle,
:class:`~orpheus.moc.geometry.MOCMesh` then rounds the requested
perpendicular ray spacing to an integer number of rays spanning the
cell, so the *effective* spacing (queryable via
:meth:`~orpheus.moc.geometry.MOCMesh.effective_spacing`) differs
slightly per angle. Tracks live in the ``tracks`` list on
:class:`~orpheus.moc.geometry.MOCMesh`, grouped by azimuthal index in
``tracks_per_azi``; the reflective links (``fwd_link`` / ``bwd_link``
on each :class:`~orpheus.moc.geometry.Track`) are resolved by
nearest-endpoint match against the reflected angle's track family.


Integration with the Eigenvalue Loop
------------------------------------

:class:`~orpheus.moc.core.MOCSolver` implements the
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` protocol,
so :func:`~orpheus.numerics.eigenvalue.power_iteration` drives
the outer eigenvalue convergence without any MOC-specific logic.
The solver's ``solve_fixed_source`` performs one characteristic
sweep over all tracks and all azimuthal / polar directions,
accumulating a flat-source contribution per region weighted by
the exponential attenuation
:math:`(1 - e^{-\tau_k})/\Sigma_t` — see
:ref:`theory-method-of-characteristics` for the full derivation.

Geometry construction currently reuses
:class:`~orpheus.geometry.mesh.Mesh1D` for the underlying radial
discretisation of concentric pin-cell regions: the MOC mesh is
built by tracking rays through the Wigner–Seitz cell that
:meth:`StructuredGeometry.wigner_seitz_pin_cell
<orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell>`
declares and :meth:`Mesh1D.from_geometry
<orpheus.geometry.mesh.Mesh1D.from_geometry>` discretises. 2-D
Cartesian assemblies are not yet supported; see the open MOC issues
for the roadmap.


API Reference
-------------

Solver entry point
~~~~~~~~~~~~~~~~~~

.. automodule:: orpheus.moc.solver
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Solver core
~~~~~~~~~~~

.. automodule:: orpheus.moc.core
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Track / segment geometry
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: orpheus.moc.geometry
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Angular quadrature
~~~~~~~~~~~~~~~~~~

.. automodule:: orpheus.moc.quadrature
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
