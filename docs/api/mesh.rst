Mesh
====

The :mod:`orpheus.mesh` package is the mesh: the discretisation overlay
on a geometry. A geometry (:doc:`/api/geometry`) gives the shape of a
problem, its regions and its boundary declarations; a mesh divides the
geometry's intervals into cells and carries the quantities derived from
that division: cell edges, the material ID of each cell, cell volumes and
face areas, and the boundary declaration on each face. Solvers receive a
mesh and build mutable, method-specific state on top of it.

The package holds three modules:

* :mod:`orpheus.mesh.structured`: the structured meshes
  :class:`~orpheus.mesh.structured.Mesh1D` and
  :class:`~orpheus.mesh.structured.Mesh2D`, and the per-region
  discretisation descriptor :class:`~orpheus.mesh.structured.RegionMesh`.
* :mod:`orpheus.mesh.factories`: the 2-D pin-cell factory
  :func:`~orpheus.mesh.factories.pwr_pin_2d` and the equal-volume
  subdivision of one interval.
* :mod:`orpheus.mesh.axis`: the per-axis primitives
  (:class:`~orpheus.mesh.axis.AxisMesh`,
  :class:`~orpheus.mesh.axis.RadialAxisMesh`,
  :class:`~orpheus.mesh.axis.AxisCoord`,
  :class:`~orpheus.mesh.axis.FaceLabel`) whose tensor product is a
  structured phase-space mesh, and the pure shape functions on tuples of
  axes.

Every public name of the three modules is exported from
:mod:`orpheus.mesh` itself (``from orpheus.mesh import Mesh1D,
RegionMesh``); :mod:`orpheus.geometry` and :mod:`orpheus.transport.mesh`
do not export them.

.. contents::
   :local:
   :depth: 2


Where the mesh sits in the layers
---------------------------------

Posing a problem is a chain of overlays: materials, then geometry, then
mesh. A mesh is a complex enough overlay on the geometry to be its own
package rather than a module of :mod:`orpheus.geometry`: its vocabulary
(cells, subdivision schemes, per-axis primitives, face inventories) is
not the geometry's, and it is the home for every kind of mesh, of which
the structured ones are the kinds that exist today.

The package belongs to the input layer, one step above the geometry
(:ref:`architecture-layering`):

* :mod:`orpheus.mesh` may import :mod:`orpheus.geometry` (shapes,
  coordinate systems, the boundary tag and laws),
  :mod:`orpheus.numerics` and :mod:`orpheus.data`, and never
  :mod:`orpheus.transport` or a method package.
* :mod:`orpheus.geometry`, :mod:`orpheus.data` and
  :mod:`orpheus.numerics` never import :mod:`orpheus.mesh`.

:file:`tests/gates/test_layer_imports.py` enforces both directions.
Binding cross sections to the cells is not the mesh's job: that is
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, the
method-agnostic material mesh at L2 in :mod:`orpheus.transport.mesh`,
which the S\ :sub:`N` hub :class:`~orpheus.sn.problem.SNProblem` and the
diffusion hub :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`
subclass.


Design Principles
-----------------

**Frozen dataclasses.**
Both :class:`~orpheus.mesh.structured.Mesh1D` and
:class:`~orpheus.mesh.structured.Mesh2D` are
``@dataclass(frozen=True)``. Once constructed, their fields cannot
be reassigned. This turns every solver entry point into a pure
function of its inputs and prevents whole classes of bugs where a
downstream routine accidentally mutates mesh state shared across
iterations.

**Equal-volume subdivision.**
Curvilinear zones (cylindrical, spherical) are subdivided by default
into **equal-volume** annuli or shells rather than equal-width cells.
With equal widths the cell volume grows with radius
(:math:`V \propto r\,\Delta r` for an annulus,
:math:`V \propto r^2\,\Delta r` for a shell), so the innermost cells
hold a small fraction of the zone's volume; equal volumes give every
cell of the zone the same share. In Cartesian geometry the two schemes
coincide.

**Precomputed volumes — the ULP escape hatch.**
:class:`~orpheus.mesh.structured.Mesh1D` accepts an optional
``precomputed_volumes`` override, and
:meth:`~orpheus.mesh.structured.Mesh1D.from_geometry` always sets it.
For an ``"equal-volume"`` region the volumes come from
:func:`~orpheus.mesh.factories._subdivide_zone`, which returns the
*algebraic* cell volume (e.g.
:math:`V_{\rm cell} = \pi(r_{\rm out}^2 - r_{\rm in}^2)/n` in the
cylindrical case) broadcast as a scalar to every cell in the region.
Deriving those volumes from the *edges* after the fact via
:func:`~orpheus.geometry.coord.compute_volumes_1d` would pass
through a ``sqrt → **2`` or ``cbrt → **3`` round trip that loses
roughly one ULP per cell and breaks the invariant "every cell in
an equal-volume region is bit-identical" at ``rtol=1e-14`` (ERR-020).
For a ``"uniform"`` region, ``from_geometry`` computes the region's
volumes from its own edges with
:func:`~orpheus.geometry.coord.compute_volumes_1d`. A mesh constructed
directly, without ``precomputed_volumes``, derives every volume from
its edges.

**Boundary declarations on the faces.**
Each face of a structured mesh carries a boundary declaration: a
:class:`~orpheus.geometry.boundary.BC` tag, an already-typed boundary
law, or ``None`` for the method's default.
:class:`~orpheus.mesh.structured.Mesh1D` has ``bc_left`` and
``bc_right``; :class:`~orpheus.mesh.structured.Mesh2D` has ``bc_xmin``,
``bc_xmax``, ``bc_ymin`` and ``bc_ymax``. The tag, the laws and the
deferred resolution of a tag by each method's hub are geometry-layer
concepts, documented on :doc:`/api/geometry`.


Construction — the geometry to mesh path
----------------------------------------

The recommended 1-D construction path is **two-layered**: declare a
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
(pure shape — a coordinate system, breakpoints, one material id per
interval and one :class:`~orpheus.geometry.boundary.BC` per boundary
point), then discretize it with :meth:`Mesh1D.from_geometry
<orpheus.mesh.structured.Mesh1D.from_geometry>` by supplying one
:class:`~orpheus.mesh.structured.RegionMesh` per interval. The mesh
starts at the first breakpoint and ends at the last, and the laws reach
``bc_left`` / ``bc_right`` by the boundary point they belong to: a slab
or a hollow cylinder or sphere gives ``(inner, outer)``, a solid
cylinder or sphere gives its one law to ``bc_right`` and leaves
``bc_left`` ``None``. The geometry
carries **no** cell counts; the discretization description enters at
exactly one point, the ``from_geometry`` call.

.. code-block:: python

   from orpheus.geometry import BC, CoordSystem, StructuredGeometry
   from orpheus.mesh import Mesh1D, RegionMesh

   geom = StructuredGeometry(
       coord=CoordSystem.SPHERICAL,
       breakpoints=(0.0, 5.0),
       mat_ids=(0,),
       boundaries=(BC.vacuum,),
   )
   mesh = Mesh1D.from_geometry(
       geom, region_meshes=(RegionMesh(n_cells=64),),
   )

This split is what lets **reference** solution generators (``Billiard``,
``MomentSpace``, ``Spectrum``, ``BasisSpace``) consume the
:class:`StructuredGeometry` directly — they need no mesh — while
**discrete production** solvers (``solve_sn`` / ``solve_cp`` /
``solve_moc`` / ``solve_mc``) consume the resulting
:class:`~orpheus.mesh.structured.Mesh1D`. The two conventional PWR
shapes that ship as :class:`StructuredGeometry` classmethods, and the
material ID convention, are on :doc:`/api/geometry`.

Per-region subdivision
~~~~~~~~~~~~~~~~~~~~~~

:class:`~orpheus.mesh.structured.RegionMesh` selects the scheme per
region. ``"equal-volume"`` (the default) routes through
:func:`~orpheus.mesh.factories._subdivide_zone`, which carries the
three coordinate-system invariants:

* **Cartesian** — equal-width cells
  :math:`x_k = x_0 + (k/n)\,(x_n - x_0)`.
* **Cylindrical** — equal-volume annuli
  :math:`r_k = \sqrt{r_0^2 + (k/n)\,(r_n^2 - r_0^2)}`.
* **Spherical** — equal-volume shells
  :math:`r_k = \sqrt[3]{r_0^3 + (k/n)\,(r_n^3 - r_0^3)}`.

It returns both the edges **and** the exact per-cell volume (a
broadcast scalar), which ``from_geometry`` passes to the frozen
:class:`Mesh1D` as ``precomputed_volumes`` — see the design principle
above. ``"uniform"`` instead lays down equal radial extents and derives
the region's volumes from its edges via
:func:`~orpheus.geometry.coord.compute_volumes_1d`.

.. note:: **Retired 1-D factory surface.**

   Phase F retired the free functions ``Zone``, ``mesh1d_from_zones``,
   ``pwr_pin_equivalent``, ``pwr_slab_half_cell``, ``homogeneous_1d``
   and ``slab_fuel_moderator`` from the factories module (then
   ``orpheus.geometry.factories``, now :mod:`orpheus.mesh.factories`).
   Their jobs are now split between the geometry layer (the
   :class:`StructuredGeometry` classmethods) and the mesh layer
   (:meth:`Mesh1D.from_geometry
   <orpheus.mesh.structured.Mesh1D.from_geometry>` +
   :class:`~orpheus.mesh.structured.RegionMesh`) — which is what removed
   the old "a factory both shapes AND meshes the problem" conflation.
   A homogeneous or two-region slab is now a one-liner
   ``StructuredGeometry`` literal, so no dedicated convenience
   function survives for it.


Structured meshes
-----------------

.. automodule:: orpheus.mesh.structured
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Two-dimensional factories
-------------------------

There is no 2-D :class:`StructuredGeometry` yet, so the 2-D Cartesian
path keeps a standalone factory:
:func:`~orpheus.mesh.factories.pwr_pin_2d` builds a
:class:`~orpheus.mesh.structured.Mesh2D` on a uniform grid with material
IDs assigned by radial distance from the pin centre (by default
``2 = fuel``, ``1 = clad``, ``0 = coolant``).

.. automodule:: orpheus.mesh.factories
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Axis primitives
---------------

The S\ :sub:`N` phase space factors as a tensor product of per-axis 1-D
meshes. :mod:`orpheus.mesh.axis` declares the per-axis primitive:
:class:`~orpheus.mesh.axis.Axis1D` is the protocol every 1-D axis
satisfies; :class:`~orpheus.mesh.axis.AxisMesh` is a Cartesian axis with
two boundary-bearing endpoints (``min`` and ``max``);
:class:`~orpheus.mesh.axis.RadialAxisMesh` is a solid radial axis
(sphere or cylinder) with one endpoint (``outer``), because the pole at
:math:`r = 0` is a coordinate singularity and not a boundary face.
:class:`~orpheus.mesh.axis.AxisCoord` is the coordinate system of one
axis, distinct from the whole-mesh
:class:`~orpheus.geometry.coord.CoordSystem` because an
:math:`(r, z)` mesh mixes a radial axis and a Cartesian one.

The shape algebra is a set of pure functions on a tuple of axes
(:func:`~orpheus.mesh.axis.spatial_shape`,
:func:`~orpheus.mesh.axis.face_labels`,
:func:`~orpheus.mesh.axis.face_shape`,
:func:`~orpheus.mesh.axis.face_outflow_ordinates`,
:func:`~orpheus.mesh.axis.n_unknowns_flat`,
:func:`~orpheus.mesh.axis.coord_system`), so a gate can exercise it on a
synthetic axis tuple without building a full phase space.
:func:`~orpheus.mesh.axis.axes_from_legacy_mesh` and
:func:`~orpheus.mesh.axis.legacy_mesh_from_axes` convert between an
axis tuple and a :class:`~orpheus.mesh.structured.Mesh1D` or
:class:`~orpheus.mesh.structured.Mesh2D`. The architectural narrative,
the axis tuple as the key of the dimension-agnostic face inventory, is
in :doc:`/theory/foundations/boundary_conditions`.

.. automodule:: orpheus.mesh.axis
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
