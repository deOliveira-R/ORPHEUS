Mesh
====

The :mod:`orpheus.mesh` package is the mesh: the discretisation overlay
on a geometry. A geometry (:doc:`/api/geometry`) gives the shape of a
problem, its regions and its boundary declarations; a mesh divides the
geometry's intervals into cells and carries the quantities derived from
that division: cell edges, a region label on each cell (with the region →
material map, from which the material ID of each cell is derived), cell
volumes and face areas, and the boundary declaration on each face. Solvers receive a
mesh and build mutable, method-specific state on top of it.

The package holds six modules:

* :mod:`orpheus.mesh.structured`: the structured meshes
  :class:`~orpheus.mesh.structured.Mesh1D` and
  :class:`~orpheus.mesh.structured.Mesh2D`.
* :mod:`orpheus.mesh.face_laws`: the face inventory rule
  :func:`~orpheus.mesh.face_laws.face_inventory` and the value
  :class:`~orpheus.mesh.face_laws.FaceLaws` both meshes store their
  boundary laws in.
* :mod:`orpheus.mesh.mesher`: the :class:`~orpheus.mesh.mesher.Mesher`,
  the meshing session that builds a ``Mesh1D`` from a geometry and
  interval rules, and the only way one is built.
* :mod:`orpheus.mesh.partition`: the interval rules
  (:class:`~orpheus.mesh.partition.CellsByCount`,
  :class:`~orpheus.mesh.partition.CellsByMaxWidth`,
  :class:`~orpheus.mesh.partition.Refined`,
  :class:`~orpheus.mesh.partition.CellEdges`) and the spacing rules
  (:class:`~orpheus.mesh.partition.EqualWidth`,
  :class:`~orpheus.mesh.partition.EqualVolume`) that place the cells of
  one interval.
* :mod:`orpheus.mesh.factories`: the 2-D pin-cell factory
  :func:`~orpheus.mesh.factories.pwr_pin_2d`.
* :mod:`orpheus.mesh.axis`: the per-axis primitives
  (:class:`~orpheus.mesh.axis.AxisMesh`,
  :class:`~orpheus.mesh.axis.RadialAxisMesh`,
  :class:`~orpheus.mesh.axis.AxisCoord`,
  :class:`~orpheus.mesh.axis.FaceLabel`) whose tensor product is a
  structured phase-space mesh, and the pure shape functions on tuples of
  axes.

Every public name of the six modules is exported from
:mod:`orpheus.mesh` itself (``from orpheus.mesh import Mesher,
CellsByCount``); :mod:`orpheus.geometry` and :mod:`orpheus.transport.mesh`
do not export them.

.. contents::
   :local:
   :depth: 2


Where the mesh sits in the layers
---------------------------------

Posing a problem is a chain of overlays: materials, then geometry, then
mesh. A mesh is a complex enough overlay on the geometry to be its own
package rather than a module of :mod:`orpheus.geometry`: its vocabulary
(cells, interval rules, per-axis primitives, face inventories) is
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

**Two spacings, and no default.**
Cells are placed in an interval at equal steps of a measure coordinate
:math:`T(r) = r^{p}`. :class:`~orpheus.mesh.partition.EqualWidth` steps
in :math:`r` (:math:`p = 1`) in every coordinate system;
:class:`~orpheus.mesh.partition.EqualVolume` steps in the coordinate
system's own measure coordinate (:math:`r`, :math:`r^2`, :math:`r^3`),
so every cell of the interval holds the same share of its measure. With
equal widths on a curved body the cell volume grows with radius
(:math:`V \propto r\,\Delta r` for an annulus,
:math:`V \propto r^2\,\Delta r` for a shell), so the innermost cells
hold a small fraction of the interval's volume. On a slab the two
spacings are one body and give one mesh. A counted rule states its
spacing: there is no default (``CellsByCount.uniform_width(n)`` and
``CellsByCount.uniform_volume(n)`` are the short spellings).

**Stored volumes, checked against the one measure.**
:class:`~orpheus.mesh.structured.Mesh1D` stores its cell volumes; it
does not recompute them from its edges. An equal-volume cell's volume
is the equal share :math:`m/n` of its interval's measure, one scalar
broadcast over the interval, because re-deriving it from the edges goes
through a ``sqrt → **2`` or ``cbrt → **3`` round trip that moves it by
about one unit in the last place (ulp) per cell and breaks the invariant
"every cell of an equal-volume interval is bit-identical" (ERR-020). Any
other cell's volume is the geometry's measure of the realised cell. The
constructor checks every stored volume against the coordinate system's
measure of its cell, :meth:`~orpheus.geometry.coord.CoordSystem.measure`,
within a band of :math:`2p + 5` ulp derived from the rounding of the
edge placement (7 on a slab, 9 on a cylinder, 11 on a sphere), and
refuses a non-positive volume. The derivation and its measurements are
on :ref:`structured-geometry-mesh`.

**Boundary laws on the faces.**
Each boundary face of a structured mesh carries a boundary law: a
:class:`~orpheus.geometry.boundary.BC` tag or an already-typed law.
Both :class:`~orpheus.mesh.structured.Mesh1D` and
:class:`~orpheus.mesh.structured.Mesh2D` store their laws as one value,
``face_laws``, a :class:`~orpheus.mesh.face_laws.FaceLaws`: a frozen,
ordered, picklable mapping from face name to law over exactly the faces
that :func:`~orpheus.mesh.face_laws.face_inventory` derives. A
``Mesh1D`` has ``xmin`` and ``xmax`` on a slab or a hollow body and
``xmax`` alone on a solid cylinder or sphere, whose centre is an
interior point; ``boundary_points`` gives the faces' positions and
``outer_law`` reads ``face_laws["xmax"]``. A ``Mesh2D`` adds ``ymin``
and ``ymax``; a solid :math:`(r, z)` mesh has no ``xmin``, its axis
:math:`r = 0` being interior. The axis primitives declare a law on each
endpoint too. ``None`` is not a law on
any of them, as on a geometry: every declaration passes one element
parser, and no consumer supplies a default law for an unstated face
(:ref:`structured-geometry-no-default-law`). The tag, the laws and the
deferred resolution of a tag by each method's hub are geometry-layer
concepts, documented on :doc:`/api/geometry`.


Construction — the geometry to mesh path
----------------------------------------

The 1-D construction path is **two-layered**: declare a
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
(pure shape: a coordinate system, breakpoints, one material id per
interval and one :class:`~orpheus.geometry.boundary.BC` per boundary
point), then mesh it with a :class:`~orpheus.mesh.mesher.Mesher`:
``partition`` takes one interval rule for every interval, or a tuple of
one rule per interval, builds the mesh and returns the mesher;
``refine(k)`` rebuilds it with ``k`` times the cells of every rule
(``k`` a power of two); ``mesh`` returns the current mesh. The mesh
starts at the first breakpoint and ends at the last, every breakpoint
is a cell edge, each cell takes the material of its interval, and the
face laws are the geometry's boundary laws. The geometry carries
**no** cell counts; the discretisation enters at exactly one point,
the rules handed to the mesher.

.. code-block:: python

   from orpheus.geometry import BC, StructuredGeometry
   from orpheus.mesh import CellsByCount, Mesher

   geom = StructuredGeometry.sphere((0.0, 5.0), (0,), outer=BC.vacuum)
   mesher = Mesher(geom).partition(CellsByCount.uniform_volume(64))
   mesh = mesher.mesh
   finer = mesher.refine(2).mesh
   assert (mesh.N, finer.N) == (64, 128)
   assert dict(mesh.face_laws) == {"xmax": BC.vacuum}

The general constructor
``Mesh1D(coord, edges, volumes, region_ids, region_materials, face_laws)`` is what the
mesher calls, and what a relabelling
(:meth:`~orpheus.mesh.structured.Mesh1D.with_distinct_cell_ids`) or the
axis adapter calls; a test or a script builds its mesh through the
mesher.

This split is what lets **reference** solution generators (``CharacteristicDerivation``,
``MomentSpace``, ``Spectrum``, ``BasisSpace``) consume the
:class:`StructuredGeometry` directly — they need no mesh — while
**discrete production** solvers (``solve_sn`` / ``solve_cp`` /
``solve_moc`` / ``solve_mc``) consume the resulting
:class:`~orpheus.mesh.structured.Mesh1D`. The two conventional PWR
shapes that ship as :class:`StructuredGeometry` classmethods, and the
material ID convention, are on :doc:`/api/geometry`.

Interval rules
~~~~~~~~~~~~~~

An interval rule says how one interval :math:`[a, b]` is divided:

* :class:`~orpheus.mesh.partition.CellsByCount` — ``n`` cells placed by
  a spacing rule;
* :class:`~orpheus.mesh.partition.CellsByMaxWidth` — the fewest cells
  whose nominal widest cell is no wider than a bound;
* :class:`~orpheus.mesh.partition.Refined` — ``k * rule``, ``k`` times
  the cells of a counted rule by the same spacing, ``k`` a power of two
  so every coarse edge stays a fine edge bit for bit;
* :class:`~orpheus.mesh.partition.CellEdges` — the edges written out,
  for a grid whose irregularity is the point; it cannot be refined.

A spacing rule places ``n`` cells at
:math:`r_j = T^{-1}\bigl(T(a) + (j/n)\,(T(b) - T(a))\bigr)`, with both
end edges pinned to the breakpoints. For
:class:`~orpheus.mesh.partition.EqualVolume` this is, per coordinate
system:

* **Cartesian** — equal-width cells
  :math:`x_k = x_0 + (k/n)\,(x_n - x_0)`.
* **Cylindrical** — equal-volume annuli
  :math:`r_k = \sqrt{r_0^2 + (k/n)\,(r_n^2 - r_0^2)}`.
* **Spherical** — equal-volume shells
  :math:`r_k = \sqrt[3]{r_0^3 + (k/n)\,(r_n^3 - r_0^3)}`.

The theory (the one measure, the volume band, the Mesher's design and
the candidates it replaced) is :ref:`structured-geometry-mesh`.

.. note:: **Retired 1-D factory and mesh surfaces.**

   Phase F retired the free functions ``Zone``, ``mesh1d_from_zones``,
   ``pwr_pin_equivalent``, ``pwr_slab_half_cell``, ``homogeneous_1d``
   and ``slab_fuel_moderator`` from the factories module (then
   ``orpheus.geometry.factories``, now :mod:`orpheus.mesh.factories`).
   Their jobs were split between the geometry layer (the
   :class:`StructuredGeometry` classmethods) and the mesh layer, then
   ``Mesh1D.from_geometry`` with one ``RegionMesh`` per interval, which
   removed the old "a factory both shapes AND meshes the problem"
   conflation. P1 step 3 of #405 (2026-09-29) retired
   ``Mesh1D.from_geometry``, ``RegionMesh``, the ``_subdivide_zone``
   helper, the ``precomputed_volumes`` field and the ``bc_left`` /
   ``bc_right`` fields in favour of the Mesher, the interval rules and
   ``face_laws``.


Structured meshes
-----------------

.. automodule:: orpheus.mesh.structured
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


The Mesher
----------

.. automodule:: orpheus.mesh.mesher
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Interval rules and spacing rules
--------------------------------

.. automodule:: orpheus.mesh.partition
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


The 2-D mesh and its factory
----------------------------

There is no 2-D :class:`StructuredGeometry` yet, so a
:class:`~orpheus.mesh.structured.Mesh2D` is constructed directly from
its edges, its material map and its laws:
``Mesh2D(edges_x, edges_y, mat_map, *, face_laws, coord=CARTESIAN)``.
``face_laws`` is keyword-only and required: any mapping from face name
to law whose keys are exactly the faces ``face_inventory`` derives
(the rule ``Mesh1D`` uses, with both ends of the second axis added).
The faces are named by :class:`~orpheus.mesh.axis.FaceLabel`, the names
S\ :sub:`N`'s resolved ``bc`` table uses. A missing face, an extra face
(``xmin`` on a solid :math:`(r, z)` mesh among them), ``None``, a
non-law object and a non-mapping are refused. The laws are stored as a
``FaceLaws`` in inventory order, so ``tuple(mesh.face_laws)`` is the
inventory:

.. code-block:: python

   import numpy as np
   from orpheus.geometry import BC
   from orpheus.mesh import Mesh2D

   mesh = Mesh2D(
       [0.0, 1.0, 2.0], [0.0, 1.0], np.zeros((2, 1), dtype=int),
       face_laws={"xmin": BC.reflective, "xmax": BC.vacuum,
                  "ymin": BC.reflective, "ymax": BC.reflective},
   )
   assert tuple(mesh.face_laws) == ("xmin", "xmax", "ymin", "ymax")

The design, the inventory table and the candidates it replaced are
:ref:`structured-geometry-face-laws`.

The one 2-D factory, :func:`~orpheus.mesh.factories.pwr_pin_2d`, builds
a Cartesian ``Mesh2D`` on a uniform grid with material IDs assigned by
radial distance from the pin centre (by default ``2 = fuel``,
``1 = clad``, ``0 = coolant``). Its keyword ``law``, the law on all four
faces, is required: a factory does not choose a physical boundary on
the caller's behalf.

.. automodule:: orpheus.mesh.factories
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Face laws
---------

.. automodule:: orpheus.mesh.face_laws
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
:math:`r = 0` is a coordinate singularity and not a boundary face. Each
endpoint's law is a required constructor argument (``bc_low`` and
``bc_high``, or ``bc_outer``), parsed like every other declaration, so
an axis's ``bc`` table holds a law on every endpoint and never ``None``.
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
