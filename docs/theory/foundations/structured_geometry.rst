.. _theory-structured-geometry:

================================================
Structured Geometry — the geometry/mesh contract
================================================

Key facts
=========

* **Two roles, one input axis.** ORPHEUS solvers split into
  *discrete production* (CP, SN, MOC, MC, ``solve_homogeneous_infinite``)
  and *continuous reference* (Billiard, MomentSpace, Spectrum,
  BasisSpace). Both consume the **same geometry layer** —
  :class:`~orpheus.geometry.structured_geometry.StructuredGeometry` +
  ``materials: dict[int, Mixture]`` — but diverge on whether they
  want a discrete mesh.
* :class:`StructuredGeometry` is **pure shape**, a frozen value with
  four keyword-only fields: the coordinate system ``coord`` (a
  :class:`~orpheus.geometry.coord.CoordSystem` member), the
  ``breakpoints`` :math:`r_0 < r_1 < \dots < r_R`, one material id per
  interval :math:`[r_k, r_{k+1}]` in ``mat_ids``, and one boundary law
  per boundary point in ``boundaries`` (a
  :class:`~orpheus.geometry.boundary.BC` tag or a typed
  ``BoundaryTraceLaw``). **No cell counts, no critical-dimension
  scalars, no energy-group count, no infinite-medium kind.** The field
  set is :ref:`structured-geometry-value`.
* **The boundary is derived, not declared.** The boundary points are
  the topological boundary of :math:`[r_0, r_R]` in the coordinate
  system, read by the derived property ``boundary_points``: two on a
  slab, one (the outer surface) on a solid cylinder or sphere, whose
  centre :math:`r = 0` is an interior point and carries no law, and two
  (inner, outer) on a hollow cylinder or sphere (:math:`r_0 > 0`). The
  constructor refuses a law count that disagrees, and refuses ``None``
  as a law (:ref:`structured-geometry-derived-boundary`).
* **Breakpoints are stored bit for bit, never re-derived.** Re-adding a
  stack's widths loses bits (``(0, .17, .45, .62, 1.0)`` comes back as
  ``(0.0, 0.17, 0.45000000000000007, 0.6200000000000001, 1.0)``), so a
  geometry states positions. A registry that publishes thicknesses uses
  :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.from_thicknesses`,
  the left fold :math:`r_{k+1} = r_k + t_k`
  (:ref:`structured-geometry-stored-breakpoints`).
* **Hollow cylinders and spheres are declarable, and every method
  refuses a declared law it would drop.** S\ :sub:`N` and diffusion
  admit only a reflective inner law on a hollow body (#511), which is
  verified to be the void cavity it models; CP admits a slab only with
  equal left and right laws, and no inner law (#513); MoC admits only a
  solid cylinder (#514); MC admits only a slab or a solid cylinder, and
  only a ``periodic`` left law (#513).
  Each refusal is a ``NotImplementedError`` at the method's own door
  (:ref:`structured-geometry-hollow-inner-law`).
* The geometry → mesh transition is **the single explicit point**
  where discretization information enters the pipeline. A
  :class:`~orpheus.mesh.mesher.Mesher` loads the geometry, divides
  each interval by an interval rule (:mod:`orpheus.mesh.partition`:
  ``CellsByCount.uniform_width(n)``, ``CellsByCount.uniform_volume(n)``,
  ``CellsByMaxWidth``, ``CellEdges``, and ``k * rule`` to refine), and
  returns a :class:`~orpheus.mesh.structured.Mesh1D`,
  ``Mesher(geom).partition(rule).mesh``, whose first edge is
  :math:`r_0` and whose last is :math:`r_R`. The mesh holds the
  coordinate system, the cells, a material per cell and a law per
  boundary face, and no geometry; the geometry owns the measure, which
  has one definition (:ref:`structured-geometry-mesh`).
* Reference solvers (``Billiard``, ``MomentSpace``, ``Spectrum``,
  ``BasisSpace``) take ``(geometry: StructuredGeometry, materials,
  **method_kwargs)`` directly via ``__init__``. They never see a
  mesh. They never see ``n_cells``. All four read the body they were
  handed through one function,
  :func:`~orpheus.derivations.common.reference_body.reference_body`,
  which classifies every geometry as exactly one of four shapes (a
  homogeneous body, a hollow body of one material, a symmetric
  reflected slab, a layered body) and knows no solver; the boundary
  laws are read by one function too,
  :func:`~orpheus.derivations.common.reference_body.specular_albedo`,
  as one specular albedo per boundary point. Each generator serves the
  shapes and the laws its solvers solve and refuses the rest through
  one door,
  :func:`~orpheus.derivations.common.reference_body.refuse_unserved`
  (the table of which generator serves which shape under which laws:
  :ref:`structured-geometry-reference-body`).
* Discrete production solvers take ``(materials, mesh, params)``
  where ``mesh`` is built by a ``Mesher``.
* Slab convention: :attr:`StructuredGeometry.domain_extent_cm` is
  :math:`r_R - r_0`, the **full slab width** on a slab (end to end).
  F_N's natural half-thickness ``a = L / 2`` is recovered inside
  :class:`MomentSpace`. On a solid cylinder or sphere the same
  property is the outer radius, and on a hollow one the shell
  thickness.
* The Sood case registry adapter is
  :meth:`La13511Case.to_geometry()
  <orpheus.derivations.continuous.sood_registry.case.La13511Case.to_geometry>`,
  which materialises a :class:`StructuredGeometry` from the case's
  ``geometry_kind`` tag (mapped to a
  :class:`~orpheus.geometry.coord.CoordSystem` member) and
  ``truth.critical_dimension_mfp`` (cm = mfp / Σ_t).
  Infinite-medium cases raise — for ``k_\infty`` use
  :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
  or :meth:`MomentSpace.solve_kinf`.
* Two non-trivial classmethods earn their keep on
  :class:`StructuredGeometry`:
  :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
  (a solid cylinder with the ``r_cell = pitch / √π`` equal-area
  transformation, its radii stored as literal breakpoints)
  and
  :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.pwr_slab_half_cell`
  (a Cartesian half-cell from the reflective fuel-centre symmetry
  plane, built by ``from_thicknesses``).


Architectural role
==================

Before Phase F the geometry layer was conflated with the registry-
truth layer (the legacy ``GeometrySpec`` carried
``critical_dimension_mfp`` / ``critical_dimension_cm`` / ``n_groups``,
which are method-of-evaluation artefacts, not geometric properties)
and with the mesh layer (the same ``GeometrySpec`` had a ``build()``
method that took a cell count). Both conflations leaked solver-tuning
parameters into reference-solver call sites that have no use for
them — a reference solver that solves a Sood-Pu sphere needs the
coordinate system and the radius in cm, it does not need a cell count
and it does not need the published critical dimension's name.

Phase F separates the three concerns into three layers:

1. **Geometry layer** —
   :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`.
   Pure shape and boundary laws: a coordinate system, breakpoints, a
   material id per interval, a law per boundary point. No cell counts.
   No scalars from a published table.
2. **Mesh layer** —
   :class:`~orpheus.mesh.structured.Mesh1D`, the
   :class:`~orpheus.mesh.mesher.Mesher` and the interval rules of
   :mod:`orpheus.mesh.partition`. Discrete representation.
   Discretization is supplied by the interval rules a Mesher applies,
   not pinned to the geometry. The mesh layer is its
   own package, :mod:`orpheus.mesh` (:doc:`/api/mesh`), which imports
   :mod:`orpheus.geometry` and is never imported by it.
3. **Registry layer** —
   :class:`~orpheus.derivations.continuous.sood_registry.case.La13511Case`,
   :class:`~orpheus.derivations.continuous.sood_registry.case.La13511Truth`.
   Published reference values (``k_eff_or_kinf``,
   ``critical_dimension_mfp``, flux ratios). The case's
   ``to_geometry()`` adapter materialises a
   :class:`StructuredGeometry` for solvers that want one.

The layer boundaries are load-bearing: they are what makes a
reference-solver call site read

.. code-block:: python

   moment = MomentSpace(
       geometry=case.to_geometry(),
       materials=case.materials,
       fn_order=10,
   )

instead of carrying a cell count it never uses.


.. _structured-geometry-value:

The value: coordinate system, breakpoints, materials, laws
==========================================================

A :class:`StructuredGeometry` is an interval of positions in one
coordinate system, cut into material intervals, with a boundary law at
every point of its boundary. Its four fields are keyword-only, so every
construction names each of them, and the dataclass is frozen:

.. list-table:: The fields of :class:`StructuredGeometry`
   :header-rows: 1
   :widths: 16 34 50

   * - Field
     - Type and constraint
     - What it carries, and why it is spelled this way
   * - ``coord``
     - a :class:`~orpheus.geometry.coord.CoordSystem` member
       (``CARTESIAN``, ``CYLINDRICAL``, ``SPHERICAL``); any other type,
       a string included, is a ``TypeError``
     - The measure of an interval (length, annulus area, shell volume)
       and the topology of its boundary. It is the same enum every mesh,
       axis and volume formula dispatches on, so the geometry and the
       mesh built from it cannot disagree on the chart.
   * - ``breakpoints``
     - at least two real numbers, finite, strictly increasing;
       :math:`r_0 \ge 0` on a cylinder or a sphere (a radius); a
       ``bool`` is a ``TypeError``; stored as a tuple of ``float``, bit
       for bit as given, except that ``-0.0`` is stored as ``+0.0``
     - The positions :math:`r_0 < r_1 < \dots < r_R` where one material
       interval ends and the next begins, the first and the last being
       the ends of the region. A slab admits any :math:`r_0`, since a
       position on a line has no preferred origin.
   * - ``mat_ids``
     - one ``int`` per interval (``R`` of them); a ``bool`` or a
       ``float`` is a ``TypeError``, a wrong count a ``ValueError``
     - The key of the material filling :math:`[r_k, r_{k+1}]`, into the
       ``materials: dict[int, Mixture]`` that consumers receive beside
       the geometry. Inside-out on a cylinder or a sphere, left to
       right on a slab. Adjacent intervals may share a material (an
       interval boundary need not be a material boundary).
   * - ``boundaries``
     - a sequence, stored as a ``tuple``, one entry per boundary point,
       each a :class:`~orpheus.geometry.boundary.BC` tag or a typed
       ``BoundaryTraceLaw``; a non-sequence (a string included) and a
       ``None`` entry are a ``TypeError``
     - The law each boundary point carries, paired one to one with
       ``boundary_points`` in the order (inner, outer). A typed law is
       admitted beside the tag because a tag cannot carry a function (a
       prescribed inflow whose source is a manufactured solution), and
       declaring such a law on the geometry is what makes it survive the
       method-mesh rebuild every public solver entry performs.

Three properties are derived, never stored: ``is_hollow`` and
``boundary_points`` (the next section), and
:attr:`StructuredGeometry.domain_extent_cm`, the width
:math:`r_R - r_0` of the interval of positions. The three sequence
fields are parsed alike: any sequence is canonicalised to a tuple (of
``float``, of ``int``, of laws), and a non-sequence, a string included,
is refused with a ``TypeError`` naming the field.

**The geometry never interprets a law.** It checks that each entry is a
law and that there is one per boundary point; what a ``BC`` tag means
(what ``"white"`` does to the returning flux, whether a method supports
``"albedo"``) is resolved by each method's mesh through its own
admission table, at solver construction.

**What is not here, and where it lives instead.** No cell counts and no
discretisation rule: those are the mesh layer's, the interval rules a
``Mesher`` applies. No critical dimension: that is a registry's
truth record. No group count: that is the materials'. No infinite
medium: an infinite medium is either a problem with no geometry
(``solve_homogeneous_infinite``, ``MomentSpace.solve_kinf``) or a finite
domain with reflective laws.

A declaration for each shape of boundary:

.. code-block:: python

   from orpheus.geometry import BC, CoordSystem, StructuredGeometry

   # A bare sphere: one interval, one law at the outer surface.
   sphere = StructuredGeometry(
       coord=CoordSystem.SPHERICAL,
       breakpoints=(0.0, 2.872),
       mat_ids=(0,),
       boundaries=(BC.vacuum,),
   )
   assert sphere.boundary_points == (2.872,)

   # A reflected slab (reflector | core | reflector): two laws.
   slab = StructuredGeometry(
       coord=CoordSystem.CARTESIAN,
       breakpoints=(0.0, 0.5, 2.5, 3.0),
       mat_ids=(1, 0, 1),
       boundaries=(BC.vacuum, BC.vacuum),
   )
   assert slab.boundary_points == (0.0, 3.0)

   # A hollow sphere: an inner and an outer law.
   shell = StructuredGeometry(
       coord=CoordSystem.SPHERICAL,
       breakpoints=(0.5, 1.0, 2.0),
       mat_ids=(1, 0),
       boundaries=(BC.reflective, BC.vacuum),
   )
   assert shell.boundary_points == (0.5, 2.0)
   assert shell.domain_extent_cm == 1.5


.. _structured-geometry-derived-boundary:

The boundary is derived, not declared
=====================================

The number of laws a geometry takes is not a property of its
coordinate system alone. The boundary of the region :math:`[r_0, r_R]`
is its **topological boundary** in its coordinate system, and that
depends on :math:`r_0`:

.. list-table:: Boundary points, and so laws, per coordinate system
   :header-rows: 1
   :widths: 28 36 36

   * - Coordinate system
     - :math:`r_0 = 0`
     - :math:`r_0 > 0`
   * - Cartesian (slab)
     - 2: :math:`(r_0, r_R)`, left and right
     - 2: :math:`(r_0, r_R)`, left and right
   * - cylindrical, spherical
     - 1: :math:`(r_R,)`, the outer surface
     - 2: :math:`(r_0, r_R)`, inner and outer

On a slab :math:`[r_0, r_R]` is an interval of a line and its boundary
is its two ends, wherever :math:`r_0` sits. On a cylinder or a sphere
the coordinate :math:`r` is a radius and the region is a disk or ball
of radius :math:`r_R` with, when :math:`r_0 > 0`, the concentric disk
or ball of radius :math:`r_0` removed. When :math:`r_0 = 0` nothing is
removed: the centre :math:`r = 0` is a point of the region's
**interior**, every neighbourhood of it lies inside the body, and a
point of the interior carries no boundary law. The boundary is the
outer surface alone. When :math:`r_0 > 0` the removed cavity has a
surface, and that surface is a boundary point of the radial interval
like any other, with its own law.

**Why the centre takes no law, not a reflective one.** In the radial
chart the point :math:`r = 0` is where the chart itself degenerates
(the areas :math:`2\pi r` and :math:`4\pi r^2` vanish there), and what a
solution must satisfy at it is regularity: a finite flux, with the
symmetry the chart imposes. That is a property of the **coordinate
chart**, which the methods' curvilinear machinery carries (an
S\ :sub:`N` radial axis, for one, has one law slot, the outer surface,
and treats the pole as a coordinate singularity rather than an
endpoint), not a choice the user makes. Declaring a law there would state a second,
possibly contradictory, condition at a point where the chart already
fixes one; the constructor therefore refuses it.

Whether the centre is in the region is decided in one place,
:meth:`CoordSystem.boundary_points <orpheus.geometry.coord.CoordSystem.boundary_points>`,
which lists the positions of the boundary of :math:`[r_0, r_R]`, inner
first: both ends on a slab and on a cylinder or a sphere with
:math:`r_0 > 0`, the outer end alone on a solid one. The geometry's
``boundary_points`` and the mesh's ``boundary_faces`` both read it; the
derived property ``is_hollow`` is true exactly for a cylinder or a
sphere with :math:`r_0 > 0`, and never for a slab, which has no centre.
One check, ``parse_boundary_laws``, requires one law per boundary point,
for the geometry's ``boundaries`` and a mesh's ``face_laws`` alike. The refusals are keyed to the three ways a law
count can be wrong, each naming the reason in its message:

* **a law at the centre** — two laws on a solid cylinder or sphere:
  *"the centre r = 0 of a solid … body is an interior point and
  carries no law"*;
* **a hollow body missing its inner law** — one law with
  :math:`r_0 > 0`: *"a hollow … body … has an inner surface, which
  needs its own law"*;
* **a slab with one law**: *"a slab has two boundary points (left,
  right)"*;

and, before the count is read, ``None`` in any position: *"None is not
a boundary law: declare the law the boundary point carries."* ``None``
is refused because it means nothing there: a geometry declares the
problem, and a default is a method's, not the problem's. Until P1 step
3b a ``None`` face on a :class:`~orpheus.mesh.structured.Mesh1D` meant
"use the solver's default"; since then the mesh refuses it with the same
message, so every face of a 1-D mesh carries a declared law (``Mesh2D``
and the S\ :sub:`N` axis tuples still admit ``None`` until step 3c).

The laws are indexed by boundary point, (inner, outer), never by a
coordinate-specific name (``left``, ``centreline``, ``outer``), so the
Mesher carries them onto the mesh unchanged for every coordinate
system: the mesh's ``face_laws`` are the geometry's ``boundaries``,
paired one to one with its ``boundary_faces``. Two boundary points give
two faces (inner, outer); one point gives one face, since the centre is
an interior point. The call site names each law through the named-face
constructors, which build the indexed tuple:
``StructuredGeometry.slab(breakpoints, mat_ids, left=, right=)``;
``cylinder(…)`` and ``sphere(…)`` with ``outer=`` and ``inner=``, the
second given exactly when the body is hollow;
``uniform_boundary(coord, breakpoints, mat_ids, law)``, one law on every
boundary point in any coordinate system; and
``from_homogeneous(width, boundary)``, a slab :math:`[0, w]` of
material 0 with one law on both faces.

The gates of these laws are the foundation rows
``tests/gates/geometry/test_structured_geometry.py::TestTheBoundaryIsDerived``
(the six cells coordinate × {:math:`r_0 = 0`, :math:`r_0 > 0`} and the
keyed refusals); the mesh routing is ``TestMeshingAGeometry`` in the
same file, whose ``test_a_hollow_body_propagates_its_inner_law`` pins the
two laws of a hollow body onto ``face_laws``, the faces onto
``boundary_faces`` :math:`= (r_0, r_R)`, and the first edge onto
:math:`r_0`.


.. _structured-geometry-stored-breakpoints:

Breakpoints are stored, never re-derived
========================================

A geometry states the **positions** of its interval boundaries, and
stores them exactly as given. The one exception is ``-0.0``, which is
stored as ``+0.0``: the two compare equal, so they are one breakpoint,
and a digest taken over the stored bits must see one value. It does not store thicknesses and add
them up, because floating-point addition does not return the positions
a thickness list was taken from:

.. code-block:: python

   import itertools
   import numpy as np

   E = (0.0, 0.17, 0.45, 0.62, 1.0)
   widths = np.diff(E).tolist()
   refolded = tuple(itertools.accumulate(widths, initial=0.0))
   assert refolded == (0.0, 0.17, 0.45000000000000007, 0.6200000000000001, 1.0)
   assert refolded != E

Two of the five positions come back one unit in the last place (ULP)
away from the numbers that were written. This is the one-ULP site the
census of the reference-solution campaign found in the test corpus
(``test_g_adjoint_reciprocity.py``, a slab whose interfaces are stated
as positions), and it is the reason for the rule: an interface position
moved by one ULP moves the cell edges and volumes of every mesh built
on it, and a reference or a snapshot pinned bitwise against the
written positions then disagrees with a mesh that was supposed to be
the same. Storing the breakpoints makes the geometry the single source
of its own interface positions.

**The thickness constructor.** Registries and factories that publish a
layered configuration as thicknesses (reflector, core, reflector; fuel
half-width, cladding, coolant) build it with
:meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.from_thicknesses`:

.. math::

   r_0 = r_0^{\rm given}, \qquad r_{k+1} = r_k + t_k,

evaluated left to right, which is ``itertools.accumulate(thicknesses,
initial=r_0)``. This is the same sequential sum the mesh construction
of the time (``Mesh1D.from_geometry``, since retired) evaluated while a
geometry stored thicknesses, so a registry's geometry
keeps its bits across the change to stored breakpoints. A thickness :math:`t_k \le 0` is refused as a
non-increasing breakpoint pair. The association is load-bearing: a
pairwise or compensated (``math.fsum``) cumulative sum differs from
the left fold on ``(0.1, 0.2, 0.3, 0.4, 0.5)``, which is the input the
gate uses to show it can tell them apart.

**The Wigner–Seitz factory states radii.** A pin cell is published as
radii (fuel, cladding, cell), so
:meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
stores ``(0, r_fuel, r_clad, r_cell)`` as literal breakpoints, with
:math:`r_{\rm cell} = {\rm pitch}/\sqrt{\pi}`. `[M]` 2026-09-29, at the
default arguments ``(0.9, 1.1, 3.6)``: the literal radii and the
thickness fold of the earlier spelling agree bit for bit,
``(0.0, 0.9, 1.1, 2.0310825007719226)``, so no mesh built on the
default pin cell moves. The half-cell slab
(:meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.pwr_slab_half_cell`)
is published as thicknesses and uses ``from_thicknesses``.

The gates are ``TestBreakpointsAreStored`` (the round trip of
``(0, .17, .45, .62, 1.0)``, with a positive control that re-adding
the widths misses) and ``TestFromThicknesses`` (the left fold at three
origins, the association control, and the refusal of a non-positive
thickness) in ``tests/gates/geometry/test_structured_geometry.py``.


.. _structured-geometry-mesh:

The mesh refines the geometry: the one measure, the Mesher, the interval rules
==============================================================================

The step from a geometry to a mesh is the one point where discretisation
enters the pipeline. Four objects take part, each in its own place:

* the **geometry**, which owns the measure of its coordinate system
  (lengths on a slab, areas per unit height on a cylinder, volumes on a
  sphere);
* the **interval rules** (:mod:`orpheus.mesh.partition`), each of which
  says how one material interval :math:`[r_k, r_{k+1}]` is divided into
  cells;
* the **Mesher** (:mod:`orpheus.mesh.mesher`), a meshing session on one
  geometry, which applies the rules, holds the mesh it built, and refines
  it;
* the **mesh**, :class:`~orpheus.mesh.structured.Mesh1D`, which holds the
  coordinate system, the cells, the material of each cell and the law on
  each boundary face, and no geometry.

A session reads like this:

.. code-block:: python

   from orpheus.geometry import BC, StructuredGeometry
   from orpheus.mesh import CellEdges, CellsByCount, CellsByMaxWidth, EqualWidth, Mesher

   # The CP and MoC default pin cell: equal-volume cells, 10 / 3 / 7 per region.
   pin = StructuredGeometry.wigner_seitz_pin_cell(r_fuel=0.9, r_clad=1.1, pitch=3.6)
   mesher = Mesher(pin).partition((
       CellsByCount.uniform_volume(10),
       CellsByCount.uniform_volume(3),
       CellsByCount.uniform_volume(7),
   ))
   coarse = mesher.mesh
   fine = mesher.refine(2).mesh           # twice the cells of every rule
   assert (coarse.N, fine.N) == (20, 40)
   assert set(coarse.edges) <= set(fine.edges)   # every coarse edge stays, bit for bit

   # One rule per interval, and the spacings may differ between intervals.
   slab = StructuredGeometry.slab(
       (0.0, 0.5, 2.5, 3.0), (1, 0, 1), left=BC.vacuum, right=BC.reflective,
   )
   mesh = Mesher(slab).partition((
       CellsByCount.uniform_width(2),
       CellsByMaxWidth(0.25, EqualWidth()),
       CellEdges([2.5, 2.6, 3.0]),
   )).mesh
   assert mesh.N == 12
   assert mesh.mat_ids.tolist() == [1, 1] + [0] * 8 + [1, 1]
   assert mesh.face_laws == (BC.vacuum, BC.reflective)
   assert mesh.boundary_faces == (0.0, 3.0)


.. _structured-geometry-one-measure:

The geometry owns the measure, and it has one definition
--------------------------------------------------------

The measure of the cells between edges :math:`r_0 < r_1 < \dots` is

.. math::

   m_j = c\,\bigl(T(r_{j+1}) - T(r_j)\bigr), \qquad T(r) = r^{d},

with :math:`d` the dimension the position sweeps (1 on a slab, 2 on a
cylinder, 3 on a sphere) and :math:`c` equal to 1, :math:`\pi` and
:math:`\tfrac43\pi`: a slab's length per unit transverse area, a
cylinder's area per unit height, a sphere's volume. :math:`T` is the
coordinate in which the measure is uniform, the
:class:`~orpheus.geometry.coord.MeasureCoordinate` of the coordinate
system. The coordinate system owns the definition
(:meth:`~orpheus.geometry.coord.CoordSystem.measure`, with
``measure_constant`` for :math:`c` and ``measure_coordinate`` for
:math:`T`), the geometry answers it through
:meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.measure`,
and the measure of an interval is the one-cell case,
``measure([r_k, r_{k+1}])``. ``compute_volumes_1d`` is the same call.
There is no second spelling anywhere in the mesh layer.

:math:`T` is evaluated on numpy arrays, whose power is correctly
rounded. Python's scalar ``b**2`` calls the C library's ``pow``, which
is not: `[M]` 2026-09-29 (the elegance review of P1 step 3a, macOS libm),
22 of 20 000 random squares differ from the correctly rounded value,
against 0 for numpy's. Adopting the one definition moved the volumes of
0 of the 414 equal-volume intervals the pre-carve capture recorded, so
bit identity with the retired subdivision helper was dropped as a
constraint rather than kept by a second spelling.

The measure belongs to the geometry because a measure means nothing
without the coordinate system that defines it. The first design of this
step stored measures in a free-standing value and was refuted on exactly
that point (:ref:`structured-geometry-mesh-refuted`).


The mesh: cells, a material per cell, a law per boundary face
-------------------------------------------------------------

:class:`~orpheus.mesh.structured.Mesh1D` is constructed from exactly
what a 1-D discretisation is:
``Mesh1D(coord, edges, volumes, mat_ids, face_laws)``. It holds no
geometry. Its construction laws, each a refusal with a keyed message:

* ``edges`` are at least two finite, strictly increasing positions, and
  :math:`r_0 \ge 0` on a cylinder or a sphere;
* ``volumes`` are positive, one per cell, and each is the coordinate
  system's measure of its cell (the formula above)
  within a band of :math:`2p + 5` units in the last place (ulp) of
  :math:`c\,T(r_{j+1})`, :math:`p` the exponent of :math:`T`;
* ``mat_ids`` are one integer per cell;
* ``face_laws`` are one law per boundary point of
  :math:`[r_0, r_R]`, inner first, each a :class:`~orpheus.geometry.boundary.BC`
  tag or a typed ``BoundaryTraceLaw``. ``None`` is refused, with the
  same message the geometry gives (*"None is not a boundary law: declare
  the law the boundary point carries."*), because both declarations go
  through one check, ``parse_boundary_laws``.

The derived quantities are ``widths``, ``centers``, ``areas``,
``boundary_faces`` (the positions paired with ``face_laws``) and
``outer_law`` (the law on :math:`r = r_R`, a slab's right face, which
is the one law collision probability, characteristics and Monte Carlo
read). Equality is bitwise over the five fields. A mesh has no hash
yet: its content identity, and the discretisation digest that keys on
it, land in step 5 of the reference-solution plan together with the
shared encoder.

**Why the volumes are stored, and why they are checked.** An
equal-volume cell's volume is stored as the equal share :math:`m/n` of
its interval's measure, because the shell between the realised edges is
not exact: the ``sqrt`` or ``cbrt`` that places an edge and the power
that re-evaluates :math:`T` there do not round-trip (ERR-020). `[M]` on
a sphere :math:`[0, 3]` in 12 equal-volume cells, the 12 stored shares
are one value and the 12 shells recomputed from the edges are 10
distinct values. A stored volume that is not checked could be anything,
though, so the constructor checks each against the shell of its cell.
The band is derived, not chosen. Each realised edge carries at most 1
ulp from forming :math:`t_j = T(a) + f_j\,\Delta T`, at most 1 ulp from
the root (allowing a :math:`\sqrt[3]{\cdot}` that is not correctly
rounded) amplified by :math:`p` when :math:`T` is re-evaluated, and
½ ulp from that evaluation: :math:`p + 1.5` ulp per edge. Two edges,
plus the subtraction, the constant and the share (about 2 ulp together),
give :math:`2p + 5`: 7 on a slab, 9 on a cylinder, 11 on a sphere.
`[M]` The largest gap over about 46 000 random legal equal-volume
partitions is 3.05, 6.53 and 7.62 ulp (the elegance review of P1 step
3b, macOS libm), while a volume from the wrong coordinate system is
:math:`10^{15}` to :math:`10^{16}` ulp off (`[M]` 6e15 in that review;
1.9e16 for a cylinder's annuli handed to a slab mesh on :math:`[0, 2]`,
2026-09-29). The first version used one band of 8 ulp
for every coordinate system; the sphere's measured 7.62 left it 0.38 ulp
of margin, which is why the band became a function of :math:`p`. A
volume that comes with its edges (a :class:`~orpheus.mesh.partition.CellEdges`
rule, or a width-spaced cell on a curved body) is the measure itself
and sits 0 ulp from it.

**Why the mesh holds no geometry.** Before settling the constructor the
user asked what a ``Mesh1D`` needs the geometry for. The measured
answer: the coordinate system and the breakpoints only; the materials
and the laws rode along, with about 30 and about 29 production reads
going through the mesh. So the mesh stores the coordinate system, the
cells, the material of each cell and the law of each face, and the
geometry stays one layer down. The laws are **per face**, not per
coordinate-specific name: in 1-D each boundary point is one face, and
the per-face form is the seed for 2-D, where one boundary is several
faces that may carry different laws. A specialised
``(geometry, edges, volumes)`` constructor was considered and not built:
its one consumer would be the Mesher's own lift, while a relabelling
(``with_distinct_cell_ids``), an adaptation or the axis adapter starts
from a mesh or from axes.


The Mesher: a meshing session
-----------------------------

:class:`~orpheus.mesh.mesher.Mesher` loads one
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry` and
holds its current mesh, the way a meshing program holds the mesh it
previews:

* ``partition(rule)`` applies one interval rule to every interval, and
  ``partition((rule_0, …, rule_{R-1}))`` one rule per interval; it builds
  the current mesh and returns the mesher, so a session chains;
* ``refine(k)`` rebuilds the mesh from ``k * rule`` for each rule of the
  last partition (:math:`k` a power of two);
* ``mesh`` returns the current mesh.

Building the mesh is the **lift** of the geometry onto the cells: each
cell takes the material of the interval it lies in, and the face laws are
the geometry's boundary laws, one per boundary point. The mesher
asserts that every breakpoint is a cell edge, bit for bit, because an
interval rule is an open protocol and a rule that misplaced an end edge
would produce cells straddling two materials. So every cell lies in
exactly one interval: the mesh **refines** the geometry's region
partition, which is the gate the reference-solution plan requires of
this step.

The Mesher is the only construction path, in production and in the
tests: ``Mesher(g).partition(CellsByCount.uniform_width(8)).mesh``.
There is no ``Mesh1D(geometry, rule)`` shorthand. Quality measures, a
preview, adaptive refinement on a field (a flux gradient) and a protocol
for external meshers are the next members of the session, filed as
#539; the session exists now because they need it.


Interval rules and spacing rules
--------------------------------

.. list-table:: The interval rules of :mod:`orpheus.mesh.partition`
   :header-rows: 1
   :widths: 26 44 30

   * - Rule
     - Cells it places in :math:`[a, b]`
     - Refines as
   * - ``CellsByCount(n, spacing)``, with ``CellsByCount.uniform_width(n)``
       and ``CellsByCount.uniform_volume(n)``
     - :math:`n` cells, placed by the spacing rule
     - ``k * rule`` (``Refined``)
   * - ``CellsByMaxWidth(h, spacing)``
     - the fewest cells whose nominal widest cell is no wider than
       :math:`h`
     - ``k * rule`` (``Refined``)
   * - ``Refined(rule, k)``, spelled ``k * rule``
     - :math:`k` times the cells of ``rule``, by the same spacing;
       ``k * (m * rule)`` is ``(k m) * rule``
     - ``k * rule``
   * - ``CellEdges(edges)``
     - the edges written out, the end edges being the interval's
       breakpoints bit for bit; the geometry gives the measures
     - refused (``TypeError``): it has no spacing rule to place new edges

A spacing rule is the measure coordinate whose equal steps place the
cells. With :math:`T(r) = r^{p}` it places :math:`n` cells at

.. math::

   r_j = T^{-1}\!\left(T(a) + \tfrac{j}{n}\bigl(T(b) - T(a)\bigr)\right),
   \qquad j = 0, \dots, n,

with both end edges pinned to the breakpoints. ``EqualWidth`` steps in
:math:`T(r) = r` in every coordinate system; ``EqualVolume`` steps in the
coordinate system's own measure coordinate (:math:`r`, :math:`r^2`,
:math:`r^3`), where the measure is uniform. The two are one body with
one parameter, the exponent :math:`p`. When the spacing's coordinate is
the measure coordinate, the cells are equal shares of the interval's
measure and each is stored as the share :math:`m/n`, which is ERR-020's
fix in its present spelling; otherwise each volume is the geometry's
measure of the realised cell. There is no default spacing: a count
without one is refused (the user's ruling of 2026-09-25), so
``CellsByCount(8)`` does not silently mean equal volume.

**The refinement factor is a power of two.** The fine fractions are
:math:`j\,\mathrm{fl}(1/(kn))`; for :math:`k = 2^m` that equals
:math:`(j/k)\,\mathrm{fl}(1/n)` exactly, a power-of-two scaling, so every
coarse edge is a fine edge bit for bit. For any other :math:`k` they
differ by up to 2 ulp (`[M]` 3909 of 7176 cases at :math:`k = 3`, the
docstring of ``_parse_factor``), and cells that do not nest are not a
refinement. The same reason makes ``2 * CellsByMaxWidth(h, s)`` the
refinement and ``CellsByMaxWidth(h / 2, s)`` not: halving the bound can
give an odd count.

**The width bound is nominal.** ``CellsByMaxWidth`` bounds the width of
the widest cell as the spacing places it: stepping in :math:`r`, every
cell's nominal width is :math:`\mathrm{fl}((b - a)/n)`, and a realised
width differs from it by up to 2 ulp of :math:`b`; stepping in
:math:`r^2` or :math:`r^3` the cells narrow outward and the first cell
is the widest.


.. _structured-geometry-issue-495:

#495, fixed at its root
-----------------------

The retired per-region descriptor had two spellings of an equal-width
slab: ``"equal-volume"`` stored the exact share :math:`(b - a)/n`, while
``"uniform"`` placed its edges with ``np.linspace`` and recomputed each
volume from them, so the volumes of one slab interval were not equal to
each other: `[M]` on :math:`[0, 3]`, 3, 5, 5 and 3 distinct volumes at
:math:`n = 5, 7, 9, 11` (the pre-carve reading recorded in the gate's
docstring). The two spellings described one mesh, and one of them was
wrong. The repair is not a patch to ``"uniform"``: on a slab
:math:`T(r) = r` is both the width coordinate and the measure
coordinate, so ``EqualWidth`` and ``EqualVolume`` are one body with
:math:`p = 1`, both store the share, and they give one mesh. The gate
``TestIssue495Mesh`` in ``tests/gates/mesh/test_mesher.py`` asserts that
``uniform_width(n)`` and ``uniform_volume(n)`` build equal meshes on
four slabs at every :math:`n` from 1 to 64 and at 100, 127, 255 and 1000,
and that they differ at every :math:`n \ge 2` on a cylinder and a sphere,
where the two coordinates differ.


.. _structured-geometry-mesh-refuted:

What was tried, and why it was refuted
--------------------------------------

The design went through four candidates before the one described above;
each is recorded with the structural reason it failed, so that a later
session does not re-attempt it.

1. **A free-standing measured partition.** Step 3a (``2c62f9f2``)
   shipped a ``Partition`` value holding the edges and the measures of
   each interval. `[REFUTED 2026-09-29]` by the review of that commit:
   the value held measures without the coordinate system that gives
   them meaning. `[M]` a cylinder's equal-volume partition was accepted
   by a slab with the same breakpoints, and
   ``Partition(([0, 1],), ([42.0],))`` was accepted. The fact that
   survives: no measure exists apart from its geometry. ``Partition``
   retired; the mesh stores the cells.
2. **A three-stage Mesher**, ``Mesher(geometry)``,
   ``.partition(rule)``, ``.mesh(cells)``, the user's first answer to
   that review. `[REFUTED 2026-09-29]` as a separate object in 1-D: it
   held nothing but the geometry, and the partition *is* the mesh, so
   the three stages collapse to one function of (geometry, rule).
3. **The Mesher as a session**, the user's second answer: the mesher
   builds a mesh it can check for quality, refine, and later adapt to a
   field, and ``.mesh`` returns it; it is the seam for a front end. This
   one landed. What changed is that the object now holds state (the
   current mesh) and mesh-to-mesh operations (refine, adapt) that are
   functions of neither the geometry nor the rule alone, which is what
   the second candidate lacked.
4. **A mesh that stores its geometry**, ``Mesh1D(geometry, partition)``,
   the specification's shape for this step. Superseded by the bare
   constructor, for the reason measured above: the mesh needs the
   coordinate system and the breakpoints, not the geometry.

Two smaller decisions went the same way. Bit identity of the
equal-volume edges with the retired subdivision helper was dropped in
favour of the one measure (0 of 414 captured intervals moved). The fixed
8-ulp volume band became :math:`2p + 5` once the sphere's measured
worst case sat within 0.38 ulp of it.


What this step does not do
--------------------------

* ``None`` as a boundary declaration is gone from ``Mesh1D`` and
  ``StructuredGeometry``; it survives on ``Mesh2D``'s four faces, on the
  axis tuples, and as S\ :sub:`N`'s ``boundary_condition=`` parameters.
  Their retirement, with the consumers' defaults, is the next sub-step
  (3c) of #405. Until then the axis adapter gives the adapter mesh of an
  axis tuple the reflective law S\ :sub:`N` resolves for an undeclared
  axis law (``ELEGANCE-DEBT[guard]``, #405) and for the inner face of a
  hollow radial axis, which has no law slot (``SCOPE-BOUNDARY[guard]``,
  #511).
* The discretisation digest and the mesh's hash are step 5's.
* Quality, preview, adaptation and external meshers are #539.


The gates
---------

``tests/gates/mesh/test_mesh1d.py`` carries the constructor's laws
(``TestConstructionLaws``, ``TestTheVolumeLaw`` with the band and the
wrong-coordinate control, ``TestTheValue`` for equality,
``TestTheRetirements`` for the retired fields and constructors, and
``TestOneBoundaryLawCheck`` for the shared law check).
``tests/gates/mesh/test_mesher.py`` carries the session (``TestTheLift``:
every cell in exactly one interval, one rule meaning that rule on every
interval; ``TestRefine``; ``TestEveryBreakpointIsAnEdge``, whose stub
rules miss their interval and are refused; ``TestSessionRefusals``;
``TestIssue495Mesh``). ``tests/gates/mesh/test_partition.py`` carries
the rules and the measure (``TestTheOneMeasure``, ``TestNesting``,
``TestSumLaw``, ``TestEqualVolume``, ``TestEqualWidth``,
``TestIssue495``, ``TestRefinement``, ``TestRefinementRefusals``,
``TestMaxWidthCount``, ``TestRuleRefusals``, ``TestCellEdges``).
``TestMeshingAGeometry`` in
``tests/gates/geometry/test_structured_geometry.py`` meshes named
geometries the way production does, and holds four of the six tests
that catch ERR-020; the other two are
``TestEqualVolume::test_measures_are_the_interval_measure_over_n`` and
``TestEqualVolume::test_each_interval_takes_its_own_share`` in
``test_partition.py``.


.. _structured-geometry-hollow-inner-law:

Hollow bodies, and the laws each method reads (#511, #513, #514)
================================================================

A hollow cylinder or sphere (:math:`r_0 > 0`) is a legal geometry:
it declares two laws, and a ``Mesher`` builds a mesh whose first edge
is :math:`r_0` and whose first face law, ``face_laws[0]``, is the inner
law. The same routing puts a slab's left law in ``face_laws[0]``.

The geometry never interprets a law, so what happens to the law on the
first face is each method's. **No method in the tree today reads every
law a geometry can declare**: S\ :sub:`N` and diffusion have no
inner-surface slot on a curvilinear axis, collision probability (CP)
and Monte Carlo (MC) read only the outer law (``outer_law``), and the method of
characteristics (MoC) reads any mesh as a solid pin cell. Before P1
step 2 each dropped the law it did not read, silently; since step 2,
each **refuses a declared law it would drop, at its own door**, with
``NotImplementedError`` naming the issue that tracks the missing
machinery. Every such refusal is a **declared scope boundary**
(``SCOPE-BOUNDARY[guard]`` in the code, with the machinery that would
retire it, the ruling, and what would overturn it): the edge of what
the method computes, not a defect in the declaration. In CP and MC
the law the method reads is resolved first against its registry, and
only then is the law it would drop refused, so an unsupported read law
still reports the registry's own ``ValueError``; MoC's refusal is on the
mesh itself and comes before any law is read.

.. list-table:: What each method reads, and what it refuses
   :header-rows: 1
   :widths: 14 30 30 26

   * - Method
     - Laws it reads
     - Refused (``NotImplementedError``)
     - Guard, issue
   * - S\ :sub:`N`, diffusion
     - slab: left and right; solid cylinder or sphere: the outer law;
       hollow: the outer law, and the cavity is reflective
     - a hollow cylinder or sphere whose inner law is not
       ``reflective``
     - ``_refuse_an_inner_law_the_radial_axis_drops``
       (``orpheus/mesh/axis.py``), #511
   * - CP
     - the outer law only
     - a slab whose left law differs from its right law; any inner law
       on a hollow cylinder or sphere
     - ``_refuse_a_law_cp_drops`` (``orpheus/cp/solver.py``), #513
   * - MoC
     - the outer law of a solid cylinder, read as the Wigner–Seitz
       cylinder of a square pin cell
     - any mesh that is not a solid cylinder: a slab, a sphere, a hollow
       cylinder
     - ``_refuse_a_mesh_moc_misreads`` (``orpheus/moc/geometry.py``),
       #514
   * - MC
     - the outer law only; ``periodic``, its registry's one kind,
       applied to every face of the unit cell
     - a hollow cylinder; a declared left or inner law that is not
       ``periodic``
     - ``_refuse_a_mesh_mc_misreads`` (``orpheus/mc/solver.py``), #513

Every guard reads declared laws only: since P1 step 3b a
:class:`~orpheus.mesh.structured.Mesh1D` refuses ``None`` as a face law,
as the geometry does, so a mesh built directly from edges declares its
laws too. Until then an undeclared law was admitted by every guard.

**S**\ :sub:`N` **and diffusion (#511).** Both build their spatial axis
through ``axes_from_legacy_mesh``, and the curvilinear axis it builds,
:class:`~orpheus.mesh.axis.RadialAxisMesh`, has one law slot,
``bc_outer``: it has no inner-surface trace. `[M]` 2026-09-25 (the P1
specification's probe ``hollow_inner_law.py``): before the guard, on a
cylinder and a sphere with :math:`r_0 = 0.5`, an inner ``vacuum`` law
and an inner ``reflective`` law gave bit-identical S\ :sub:`N` fluxes,
while the slab control, whose left law has a slot, moved. The guard
sits at the one adapter, so it covers both routes to a hollow mesh (the
geometry and a mesh built directly from edges) and both methods. It
retires when the radial axis gains an inner-surface trace, which is the
curvilinear S\ :sub:`N` build the user has deferred.

**The admitted inner law is a void cavity, verified as if derived.**
Admitting ``reflective`` is correct only if what the methods compute
for it is a physical answer, and it is: in one-dimensional cylindrical
or spherical symmetry, a ray that enters a void core crosses it along
a chord and leaves at the same impact parameter with its radial
direction cosine reversed, which is exactly specular reflection at
:math:`r_0`. So a hollow body with a reflective inner law must equal a
**solid** body whose core :math:`[0, r_0]` is empty. The gate
``test_a_reflective_inner_law_is_a_void_cavity`` makes the core a pure
absorber of total cross section :math:`\varepsilon` (no source in the
core), measures the maximum relative gap between the two bodies' shell
scalar fluxes at :math:`\varepsilon = 10^{-6}` and :math:`10^{-4}`, and
extrapolates it linearly to :math:`\varepsilon = 0`. `[M]` 2026-09-29,
the gate's own fixture (a 2-group material on the shell
:math:`[0.5, 2.0]` meshed 8 cells, the core 4 cells, vacuum outer law;
Gauss–Legendre 8 on the sphere, ``folded_product(n_mu=4, n_phi=8)`` on
the cylinder):

.. list-table:: Hollow reflective body against a solid body with an absorbing core
   :header-rows: 1
   :widths: 16 18 18 22 26

   * - Geometry
     - gap / ε at ε = 1e-6
     - gap / ε at ε = 1e-4
     - gap extrapolated to ε = 0
     - control: black core, σ = 50
   * - sphere
     - 0.3263
     - 0.3263
     - :math:`2.0\times10^{-11}`
     - 0.36
   * - cylinder
     - 0.9983
     - 0.9982
     - :math:`1.4\times10^{-10}`
     - 0.55

The gap is linear in :math:`\varepsilon` (its ratio to
:math:`\varepsilon` is constant to four digits across two decades) and
extrapolates to zero at round-off, so the discrete void-cavity answer
is the hollow reflective answer; a black core, which is not a void,
differs at :math:`O(1)`. The gate asserts the growth (the larger gap
exceeds ten times the smaller), the intercept (below :math:`10^{-3}`
of the larger gap) and the control (above 0.1). The gate's docstring
records the same linearity on the sphere under Gauss–Legendre 8 and 16
and core meshes of 2 and 8 cells.

**CP (#513).** CP reads only the outer law, the outer cell surface.
`[M]` 2026-09-25 (the specification's probe ``cp_slab_left_law.py``):
a white, a vacuum and an undeclared left law gave the same
:math:`k = 1.8749980808246423` to all sixteen digits. On a slab the
only declaration that means what CP computes is therefore left law
equal to right law, and that is the one admitted. On a hollow cylinder
or sphere what CP realises at the inner surface is not established, so
no inner law is admitted. Which law CP realises on a slab's left face
is itself open, and recorded on #513: `[M]` 2026-09-25 (probe
``cp_mirror.py``), a fuel | moderator slab with both faces white and its
mirror image give :math:`k` differing by :math:`5.0\times10^{-5}`
relative, and neither equals the mirrored double slab, so the left face
is neither the right law nor a mirror.

**MoC (#514).** ``MOCMesh`` reads a :class:`~orpheus.mesh.structured.Mesh1D`
as the Wigner–Seitz cylinder of a square pin cell: region 0 is a disk
and the pitch is ``edges[-1] * sqrt(pi)``. `[M]` 2026-09-25 (the probe
``partial_law_readers.py``): a Cartesian mesh with the cylinder's edges
gave exactly the cylinder's :math:`k`, and a hollow cylinder with
:math:`r_0 = 0.1` was realised inconsistently (the region areas read
``edges[0]`` while the tracks treat region 0 as a disk). The refusal is
on the mesh, before any law is read, and it retires with a 2-D
geometry value for MoC to read.

**MC (#513).** ``MCMesh`` reads only the outer law and applies its
registry's one kind, ``periodic``, to every face of the unit cell.
`[M]` 2026-09-25 (the same probe): a slab with a ``vacuum`` left law and
a ``periodic`` right law constructed and ran as periodic. A declared
left or inner law other than ``periodic`` is refused. So is a hollow
cylinder: MC's material lookup clamps a radius below the first edge to
region 0, so the cavity would be filled with the innermost material
(`[M]` 2026-09-29, qa: breakpoints ``(0.5, 1.0, 2.0)`` returned
material 7 at :math:`r = 0`).

The gates are ``tests/gates/mesh/test_hollow_inner_law.py`` (the
S\ :sub:`N` refusal for the inner laws ``vacuum``, ``white`` and
``albedo(0.5)`` on a cylinder and a sphere; the direct-mesh route;
diffusion; ``test_an_undeclared_inner_law_is_the_reflective_one``,
which until P1 step 3b recorded that an undeclared inner law and a
declared ``reflective`` one gave the same S\ :sub:`N` flux bit for bit,
and now asserts that ``None`` is refused, which leaves
``inner=BC.reflective`` as the one spelling of the cavity; and the
void-cavity witness above) and
``tests/gates/mesh/test_dropped_laws_are_refused.py`` (CP: a slab whose
laws differ refused both ways, an equal left law built, a hollow
cylinder and sphere with an inner law refused; MoC: a slab and a hollow
cylinder refused, the solid pin cell built; MC: a ``vacuum`` or
``reflective`` left law refused, a ``periodic`` one built, a hollow
cylinder refused).


.. _structured-geometry-reference-body:

Reference generators read one of four body shapes, and its laws
================================================================

Four continuous reference generators take a geometry directly:
``Spectrum`` (singular eigenfunctions, :doc:`/theory/references/singular_eigenfunction`),
``MomentSpace`` (F\ :sub:`N`, :doc:`/theory/references/fn_method`),
``BasisSpace`` (Galerkin spectral, :doc:`/theory/references/galerkin_spectral`)
and ``Billiard`` (trajectory resolvent,
:doc:`/theory/references/trajectory_resolvent`). Each wraps a set of
solver functions, and each solver function solves a few body shapes:
the shapes the reference literature states its benchmarks on. A
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry` can
describe more than any one generator solves, in its body and in its
boundary laws, so a generator has to read which shape it was handed
under which laws, route what it serves to the right solver, and refuse
the rest with a stated scope. Read naively, a
geometry of several materials would be solved as its first material
alone, and a hollow cylinder or sphere, whose width :math:`r_R - r_0` is
a shell thickness, as a solid of the wrong radius, both without an
error.

**Why four shapes, and why these.** The benchmark specification the
reference generators are measured against is Sood, Forster and
Parsons :cite:`SoodForsterParsons2003`, the edition the corpus cites;
the 1999 report :cite:`SoodLA13511_1999` numbers the same problems the
same way (:ref:`sood-registry-editions`). It
specifies 75 problems, and 20 of them are multi-media: symmetric
three-region slabs (problems 4, 25, 26), a slab reflected on one side
only (3), reflected cylinders (9, 10, 27, 28), reflected spheres (16,
18, 20), a four-region Fe / U / Fe / Na slab (30), two-group reflected
slabs (58 to 61) and infinite slab lattice cells (63 to 66). `[M]`
2026-09-29, the explorer's census of the report, recorded in #536 with
the table number of each problem; the registry
(:doc:`/theory/references/sood_registry`) holds 0 of the 20. The four
shapes are the distinctions the solvers under these generators draw
over that population:

* a **homogeneous body**, which every family solves;
* a **hollow body of one material**, which the trajectory resolvent
  solves (``solve_greens_function_hollow_sphere`` and
  ``solve_greens_function_annulus``);
* a **symmetric reflected slab**, which the F\ :sub:`N` reflected-slab
  solver of Neshat and Maiorino (1980) solves in one group
  (:func:`~orpheus.derivations.continuous.fn_method.slab.reflected.solve_fn_slab_reflected_critical`,
  :ref:`the F_N page's reflected-slab section <fn-method-reflected-slab>`);
  Sood's problems 4, 25 and 26 have this shape;
* a **layered body**, everything else: the reflected cylinders and
  spheres, which the trajectory resolvent's multi-region solvers
  solve when the body is solid, and the one-sided and four-region
  slabs (problems 3 and 30), which no solver in the tree solves.

The reflected slab is a shape of its own, rather than a layered slab
with a predicate beside it, because it is the one layered slab a
solver takes: its input is exactly one core, one reflector material and
one reflector thickness. The discriminating neighbours are problems 3
and 30, which hold two materials in a slab and are not this shape.

The classification
------------------

:func:`~orpheus.derivations.common.reference_body.reference_body` is
the one place a geometry is read as a body. It is **total**: every
geometry is exactly one ``ReferenceBody`` (the union of the four
shape classes below), and it raises for none. First, adjacent intervals of one material are
merged into one material **run**, because an interior breakpoint
between two intervals of the same material is not a material
interface (a sphere of intervals ``(0, 1)`` and ``(1, 2.5)`` of one
material is one body of radius 2.5). Then, with :math:`n` the number
of runs:

.. list-table::
   :header-rows: 1
   :widths: 24 40 36

   * - Shape
     - The geometry
     - Fields
   * - :class:`~orpheus.derivations.common.reference_body.HomogeneousBody`
     - :math:`n = 1`; a slab, or a solid cylinder or sphere
       (:math:`r_0 = 0`)
     - ``coord``, ``extent_cm`` (the full slab width, or the radius),
       ``mat_id``
   * - :class:`~orpheus.derivations.common.reference_body.HollowBody`
     - :math:`n = 1`; a hollow cylinder (an annulus) or a hollow
       sphere (:math:`r_0 > 0`)
     - ``coord``, ``inner_radius_cm``, ``outer_radius_cm``, ``mat_id``
   * - :class:`~orpheus.derivations.common.reference_body.ReflectedSlab`
     - a slab with :math:`n = 3` whose outer runs share one material
       and one width
     - ``core_width_cm``, ``reflector_width_cm``, ``core_mat_id``,
       ``reflector_mat_id``; ``coord`` is Cartesian
   * - :class:`~orpheus.derivations.common.reference_body.LayeredBody`
     - every other geometry of :math:`n \ge 2` runs: a layered
       cylinder or sphere, solid or hollow, or a slab that is not a
       symmetric reflected slab
     - ``coord``, ``breakpoints`` (the run boundaries
       :math:`r_0 < \dots < r_n`), ``mat_ids`` (one per run);
       the property ``is_hollow``

A slab has no centre, so its :math:`r_0` is a translation: a
one-material slab on :math:`[3, 5]` is the homogeneous body of width 2,
and a slab is never hollow
(:attr:`~orpheus.geometry.structured_geometry.StructuredGeometry.is_hollow`).

**Equal widths, up to the breakpoints' rounding.** A registry states a
layered slab as thicknesses, and
:meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.from_thicknesses`
stores the breakpoints as a sequential sum, so the right reflector's
width :math:`r_3 - r_2` is a difference of two rounded sums and is not,
in general, bit-equal to the left one (`[M]` 2026-09-29: the
thicknesses ``(0.1, 0.5, 0.1, 0.2)`` give a last breakpoint of
``0.8999999999999999``). Each of the four breakpoints carries at most
one unit in the last place of the largest magnitude
:math:`s = \max(|r_0|, |r_3|)`, so the two widths are compared with an
absolute tolerance of :math:`4\,\mathrm{ulp}(s)`; a larger difference
is a real asymmetry, and the slab is a layered body. The gate is a
Sood-style slab built from the thicknesses ``(0.5, 0.8603, 0.5)``,
read as a reflected slab.

**Why the classification knows no solver.** The shape is a fact about
the geometry; which shapes a family serves is a fact about its
solvers. Kept apart, a solver that lands changes one owner's match,
never the classification, and the breakpoints and materials are read
in one place, so "hollow", "layered" and "symmetric" cannot be derived
four slightly different ways. The alternative the ruling rejected was
a reader per generator with no shared type: each of the four would
re-derive those three predicates from breakpoints and material ids,
which is the twin the shared reading exists to prevent.

The boundary laws: one specular albedo per boundary point
---------------------------------------------------------

The body shape is half of what a generator is handed; the other half is
the law at each boundary point of the geometry. Every solver under these
four generators parametrises a boundary the same way, by one **specular
albedo** :math:`\alpha \in [0, 1]`: the fraction of the arriving
angular flux returned into the mirror direction, :math:`\alpha = 0` for
vacuum and :math:`\alpha = 1` for a perfect mirror. The singular
eigenfunction solvers call it :math:`R` (Atalay 1997), the trajectory
resolvent :math:`\alpha`. A law is read as that albedo in one place,
:func:`~orpheus.derivations.common.reference_body.specular_albedo`
(and :func:`~orpheus.derivations.common.reference_body.specular_albedos`
for every boundary point of a geometry, inner first):

.. list-table::
   :header-rows: 1
   :widths: 55 45

   * - The declared law
     - Its specular albedo
   * - ``BC.vacuum``, :class:`~orpheus.geometry.boundary.VacuumInflow`
     - 0
   * - ``BC.reflective``
     - 1
   * - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - its ``albedo``
   * - ``BC("partial", {"albedo": a})``
     - ``a``
   * - :class:`~orpheus.geometry.boundary.AlbedoBoundary` whose
       re-emission is
       :class:`~orpheus.geometry.boundary.SpecularReturn`
     - its ``albedo``
   * - any other law: white, an ``AlbedoBoundary`` returning
       isotropically or with no stated re-emission, periodic, a
       prescribed inflow
     - refused through the door

A refused law returns neutrons in an angular shape a specular albedo
cannot express, so no value of :math:`\alpha` would be an honest
reading of it. `[M]` 2026-09-29, by calling the reader on each law:
``ReflectiveBoundary(albedo=0.7)`` reads 0.7,
``AlbedoBoundary(0.3, SpecularReturn())`` reads 0.3, and ``BC("white")``,
``BC("periodic")`` and ``AlbedoBoundary(0.3, IsotropicReturn())`` are
refused.

**Why the laws are part of what a generator serves.** Before the laws
were read, two generators answered a question they had not been asked,
and neither said so. ``MomentSpace`` on a reflected slab whose outer
faces were declared mirrors returned the vacuum core half-thickness, bit
for bit, because its solver has vacuum outer faces and nothing read the
declared law (the elegance review of 2026-09-29; the gate
``test_moment_space_refuses_a_reflected_slab_with_mirror_faces`` is its
witness). ``Billiard`` took its albedos from a separate ``alpha``
constructor parameter, a second declaration of the boundary law that
won over the geometry's own, so a hollow sphere declared
(reflective, vacuum) was solved with whatever ``alpha`` said (the qa
review of 2026-09-29). The ruling (the user, 2026-09-29) made the
geometry's laws the only declaration: ``Billiard``'s ``alpha``
parameter and its ``with_alpha`` copy method retired, and each
generator's served pattern names the laws as well as the shapes.

Which generator serves which shape, under which laws
----------------------------------------------------

Each generator matches on the shape and the albedos and either routes
them or refuses them. `[M]` 2026-09-29, by constructing each generator
on each shape and law set (one group, :math:`\Sigma_t = 1`; the solid
bodies of radius or width 2, the hollow ones on :math:`[0.5, 2]`; the
laws vacuum, mirror, partial with :math:`\alpha = 0.5`, and a slab with
unequal faces):

.. list-table::
   :header-rows: 1
   :widths: 18 20 14 20 28

   * - Shape
     - ``Spectrum``
     - ``BasisSpace``
     - ``MomentSpace``
     - ``Billiard``
   * - homogeneous slab
     - served with one albedo :math:`R` on both faces; unequal faces
       refused
     - vacuum faces only
     - vacuum faces only
     - any albedos: equal faces are ``slab``, unequal faces
       ``slab_asymmetric`` (the two-surface billiard)
   * - homogeneous solid sphere
     - served with the outer albedo :math:`R`
     - vacuum only
     - vacuum only
     - any albedo (``sphere``)
   * - homogeneous solid cylinder
     - vacuum only (bare); a reflected cylinder refused
     - refused (out of pillar)
     - refused (out of pillar)
     - any albedo (``cylinder``)
   * - hollow sphere, annulus
     - refused
     - refused
     - refused
     - ``hollow_sphere``, ``annulus``, with
       :math:`(\alpha_{\rm in}, \alpha_{\rm out})` from the two laws
   * - symmetric reflected slab
     - refused
     - refused
     - vacuum outer faces only; one group with one :math:`\Sigma_t`;
       critical core half-thickness only
     - refused
   * - layered solid sphere, cylinder
     - refused
     - refused
     - refused
     - any outer albedo (``sphere_mr``, ``cylinder_mr``)
   * - layered slab, hollow layered body
     - refused
     - refused
     - refused
     - refused

In every column a law with no specular albedo (white, periodic, an
isotropic return) is refused. The names in parentheses are
``Billiard.geometry_kind``, the key of its dispatch onto the
``solve_greens_function_*`` functions.

**Billiard.** The albedos are the geometry's laws and nothing else:
``Billiard`` has no albedo parameter, and its derived
``alpha_payload`` is built from
:func:`~orpheus.derivations.common.reference_body.specular_albedos`.
The shape and the albedos are read in one match. A homogeneous slab
whose two faces declare one albedo is the one-surface billiard
(``{"alpha"}``), and one whose faces differ is the two-surface billiard
(``{"alpha_left", "alpha_right"}``, closure rank 2); a solid sphere or
cylinder reads its one outer albedo; a hollow body reads
``{"alpha_in", "alpha_out"}`` from its inner and outer laws. A solid
layered sphere routes to ``solve_greens_function_sphere_mr`` and a
solid layered cylinder to ``solve_greens_function_cylinder_mr``; the
cylinder arm is new with this change, the solver it reaches is not
(:ref:`peierls-greens-cylinder-mr`). A layered body's cross sections are
stacked one mixture per run, ``sigma_t`` and ``nu_sigma_f`` of shape
``(n_runs, G)`` and ``sigma_s`` of shape ``(n_runs, G, G)``, and its
geometry payload is the outer radius of each run. A hollow layered body
and every layered or reflected slab are refused: no trajectory-resolvent
solver takes them. ``solve_fixed_source`` is built for ``sphere_mr``
only; it returns the total scalar flux (the sum over groups) and the
per-group fluxes in its metadata (ERR-091 records the arm's earlier
defect, :doc:`/theory/verification/error_catalog`).

**MomentSpace.** Its bare slab and sphere solvers and the reflected-slab
solver all have vacuum outer faces, so every law must read as albedo 0
(:func:`~orpheus.derivations.common.reference_body.require_vacuum`). A
symmetric reflected slab routes to
:func:`~orpheus.derivations.continuous.fn_method.slab.reflected.solve_fn_slab_reflected_critical`.
That solver is one-group and puts the core and the reflector on one
mean-free-path scale, so ``MomentSpace`` refuses, at construction, a
reflected slab of more than one group or whose two media have
different total cross sections. Sood states 11 of its 12 one-group
multi-media problems this way (its Tables 2, 9 and 13; problem 30 is
the exception; the census is in #536). The question the solver answers
is the critical core half-thickness in mean free paths, so, as for the
bare slab, the geometry's core width is the unknown and is not read;
only the reflector width is, converted to mean free paths by the shared
:math:`\Sigma_t`. Each medium's :math:`c` is
:attr:`~orpheus.data.macro_xs.mixture.Mixture.scattering_ratio`, the
one definition of the Case–Zweifel secondaries per collision
:math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`, which ``MomentSpace.c``
also reads. The result is a
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
with ``parameter_kind="core_half_thickness_mfp"``. Flux reconstruction
of a reflected slab (the two-media Peierls integral) is not built and
is refused. A cylinder is refused: the F\ :sub:`N` cylinder is out of
pillar (Westfall and Metcalf 1972), and the refusal points to the
singular-eigenfunction cylinder,
:mod:`orpheus.derivations.continuous.singular_eigenfunction.cylinder`,
and to #170.

**BasisSpace** serves a homogeneous slab or sphere with vacuum
boundaries: the Galerkin-spectral method as shipped is bare-critical.
A cylinder is refused for the same reason as in ``MomentSpace``.

**Spectrum** serves a homogeneous body, and reads the reflection
coefficient :math:`R` of Atalay (1997) from the laws. Atalay's slab
solver puts one :math:`R` on both faces, so a slab whose two faces
declare different albedos is refused; the sphere solver takes :math:`R`
on the outer surface; the cylinder solver is bare (Westfall and Metcalf
1972), so a cylinder with a nonzero albedo is refused (the reflected
cylinder of Westfall and Metcalf 1973 :cite:`WestfallMetcalf1973` is
not built). :math:`R = 1` is admitted at construction and refused by
the slab and sphere solvers when ``solve_critical`` runs: under perfect
reflection the size drops out of the criticality condition, and Atalay
omits :math:`R = 1` from his tables. The laws are read, and refused, at
construction. Reflected slabs and spheres as two media, and the reflected
cylinder, are not built under ``Spectrum``; the refusal names them.

The one door
------------

Every refusal of a shape or a law goes through
:func:`~orpheus.derivations.common.reference_body.refuse_unserved`
``(what, owner=, missing=)``, which raises ``NotImplementedError``
worded *"<owner> does not solve a <what>: <missing> (#536)."* ``what``
describes the refused configuration: a body shape
(:func:`~orpheus.derivations.common.reference_body.describe`), a
boundary law, or a precondition (for instance *"Billiard does not solve
a layered cartesian body of 2 material runs (0, 1): a
trajectory-resolvent solver for a layered or reflected slab or a hollow
layered body (#536)."*, or *"MomentSpace does not solve a body with
reflecting boundaries (specular albedos (1.0, 1.0)): …"*). The
``MomentSpace`` and ``BasisSpace`` cylinder refusals, the reflected
slab's preconditions (one group, one :math:`\Sigma_t`) and its flux
reconstruction go through the door as well. Refusals of a material
property (a multi-group mixture where a solver is one-group, an
anisotropy order out of pillar) and the solve-time scope of an arm
(``solve_fixed_source`` outside ``sphere_mr``) stay their own
``NotImplementedError``, raised where that property is read.

One door gives the refusals one spelling and one entry in the guard
ledger. The door is a ``SCOPE-BOUNDARY[guard]`` (``coding-standards``,
"A guard is elegance debt"): the machinery it stands in for is the
missing solver, which ``missing`` names per owner; the ruling is the
user's, of 2026-09-29; it moves when a solver for the refused
configuration lands in the owner's family, and #536 lists the Sood
problems that would then be reachable. The refusal is
``NotImplementedError``, not ``ValueError``, because the geometry is
well formed and the scope is the generator's.

Gates and evidence
------------------

``tests/gates/derivations/test_reference_body.py`` (59 rows, all
``foundation``; `[M]` 2026-09-29, ``python -O -m pytest``, 59 passed):

* **The classification** (13 rows): the three homogeneous solids; the
  merge of intervals of one material; the two hollow bodies; the
  reflected slab built from thicknesses (the rounding case above);
  four slabs that are layered and not reflected (unequal widths,
  different reflector materials, one-sided as in problem 3, four runs
  as in problem 30); a layered solid sphere; a layered hollow
  cylinder.
* **The shape refusals** (17 rows): every owner on every shape it does
  not serve, each asserting the refusal starts with the owner's name
  and cites #536.
* **The law reader** (8 rows): the six specular spellings read as their
  albedo; white and an ``AlbedoBoundary`` with no stated re-emission
  refused.
* **The served laws** (9 rows): ``MomentSpace`` refuses a reflected
  slab with mirror faces (the witness of the silent vacuum answer
  above); ``MomentSpace`` and ``BasisSpace`` refuse a partially
  reflecting sphere; ``Spectrum`` refuses a slab with unequal faces and
  a reflected cylinder; ``Billiard`` reads a slab with unequal faces as
  the two-surface billiard and a hollow sphere's two laws as
  :math:`(\alpha_{\rm in}, \alpha_{\rm out})`; ``Spectrum`` and
  ``Billiard`` refuse a white law.
* **The reflected slab's preconditions and route** (5 rows): unequal
  :math:`\Sigma_t`, two groups and flux reconstruction, each refused;
  the symmetry tolerance on the thicknesses ``(0.3, 1.1, 0.3)``, whose
  reflector widths differ by a fraction of an ulp, with a 5-ulp
  asymmetry as the control; and the cm-to-mean-free-path conversion,
  where :math:`\Sigma_t = 2` with a 0.25 cm reflector gives the
  :math:`\Sigma_t = 1`, 0.5 cm answer bit for bit.
* **The routes** (5 rows): ``MomentSpace`` on a reflected slab returns
  the bare solver's :math:`\tau_c` bit for bit; ``Billiard`` on a
  two-group layered solid sphere and cylinder returns the bare
  multi-region solver's :math:`k` bit for bit (the done-when of #190);
  ``Billiard`` on a hollow sphere and an annulus picks the hollow arm
  and its payload (#421).
* **The body's material** (1 row): ``Billiard`` takes its cross
  sections from the body's material id, not from key 0 of the
  materials dict.
* **ERR-091** (1 row): the multi-region sphere's fixed-source arm
  reports two groups on a two-group source and returns their sum as the
  scalar flux.

The route rows are routing claims: they establish that the facade
reaches the solver unchanged, not that the solver is right. Each
solver's own verification is on its family's page. The published
values for the reflected slab pass through the facade:
``tests/gates/cross_method/test_eigenvalue.py::test_fn_reflected_slab_matches_truth``
drives its adapter through ``MomentSpace`` over four cases, Sood's
problem 4 (:math:`\tau_c = 0.43015` mfp, the 2003 edition, tolerance
:math:`10^{-3}`) and three Neshat–Maiorino cases (`[M]` 2026-09-29: 4
passed). That tolerance cannot tell the 2003 value from the 1999
edition's 0.43014; the gate that can is on the F\ :sub:`N` page
(:ref:`fn-method-reflected-slab-facade`).

The gates were mutation-tested (`[M]` 2026-09-29, the implementation
session's battery, before the laws joined the served pattern): the
positive control, the old shared refusal of every multi-material or
hollow geometry put back, reddened all 5 route rows and 34 rows in all,
and each of the battery's five other arms reddened the row it was aimed
at.

**What stays out of scope**: transcribing the 20 multi-media problems
into the registry (#536); the F\ :sub:`N` cylinder (#170); a solver
for problems 3, 30, 58 to 61 and 63 to 66, which no family has.


End-state spot checks
=====================

These are the canonical call sites the reset enables.

Production user — no registry, no truth, no critical anything::

    geom = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=3.6,
    )
    mesh = Mesher(geom).partition((
        CellsByCount.uniform_volume(10),
        CellsByCount.uniform_volume(3),
        CellsByCount.uniform_volume(7),
    )).mesh
    materials: Materials = {2: uo2_fuel(), 1: zircaloy(), 0: borated_water()}
    result = solve_cp(materials, mesh, CPParams())
    print(result.keff)

Registry consumer — Sood case::

    case = SOOD2003_CASES["Ua-1-0-SP"]
    geom = case.to_geometry()
    mesh = Mesher(geom).partition(CellsByCount.uniform_volume(64)).mesh
    result = solve_cp(case.materials, mesh, CPParams())

Reference solver direct — NO mesh, NO Mesher::

    moment = MomentSpace(
        geometry=case.to_geometry(),
        materials=case.materials,
        fn_order=10,
    )
    ref_sol = moment.solve_critical()

Infinite-medium ``k_inf`` — no geometry at all::

    mix = case.materials[0]
    result_inf = solve_homogeneous_infinite(mix)
    print(result_inf.k_inf)

F_N ``k_inf`` — same shape, just a Mixture::

    k_inf = MomentSpace.solve_kinf(mix)


Design rationale and references
===============================

The full architectural reset rationale, including the locked design
decisions, the per-solver migration record, and the inventory of
removed surfaces, lives in the implementation plan
``.claude/plans/dazzling-cuddling-boot.md``. Specifically:

* Locked decision 1 (``Materials`` type alias) — zero migration cost,
  matches every solver.
* Locked decision 2 (``Mixture.eg`` is ``Optional``) — synthetic
  Sood-style XS no longer fabricates a fake energy grid.
* Locked decision 3 (slab convention) — full slab width, accepting
  ULP drift from the F_N half-thickness path inside ``MomentSpace``.
* Locked decision 5 (the geometry carries no discretisation; spelled at
  the time as "``Region`` is geometry-only") — no ``n_cells`` on the
  geometry; cell counts lived on ``RegionMesh`` at the mesh layer, one
  per interval. Since P1 step 3 (2026-09-29) they are the interval
  rules a ``Mesher`` applies (:ref:`structured-geometry-mesh`).
* Locked decision 6 (``Mesh1D.from_geometry``) — the single explicit
  point where discretization enters. The point is unchanged; since P1
  step 3 it is ``Mesher(geometry).partition(rules)``, and
  ``from_geometry`` is retired.

Phase F replaces the now-deleted ``transport_solver_protocol.rst`` (the
Phase D casualty — that protocol conflated discrete and reference
solver roles) as the documentation entry point for this architectural
contract.


Connection coefficients (reduced streaming operator)
====================================================

Connection coefficients are **differential-geometric data of the
coordinate chart**, not solver-specific.  In SO(3)-charts language,
the spherical redistribution term :math:`(1-\mu^2)/r\,\partial_\mu`
and the cylindrical redistribution term
:math:`-(1/r)\,\partial_\varphi(\xi\,\cdot)` are the **same
connection-coefficient operator** viewed in two coordinate charts.
A curvilinear S\ :sub:`N` :term:`sweep <sweep>` marches through this data:

* **chord lengths**: cell radial widths
  (:attr:`~orpheus.mesh.structured.Mesh1D.widths`),
* **face areas**: :math:`A_{i+1/2} = 4\pi r_{i+1/2}^2` (sphere) or
  :math:`2\pi r_{i+1/2}` (cylinder),
* **the geometry factor** :math:`\Delta A_i / w_n` that ensures
  :term:`per-ordinate <ordinate>` flat-flux consistency,
* **the** :math:`\alpha` **dome recursion**, and
* **the Morel--Montry angular closure** :math:`\tau_{mm}`
  (one geometry-free formula :eq:`morel-montry-closure`, read against a
  per-geometry cell partition :eq:`angular-cell-partition` — see below).

Per Cardinal Rule 2 (architecture is critical), the same data **MUST
NOT** be duplicated across the solvers that need it.
:class:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator`
lifts the math into **one** primitive rather than a per-solver one.

⛔ That sentence read *"in* ``orpheus.geometry.reduced_operator``
*lifts the math into a* **geometry-layer** *primitive rather than a
solver-side one"* until 2026-08-28.  The single-sourcing half stands and
is the Cardinal-Rule-2 point; the **layer** half was refuted by the
un-weld arc's P4.4.  `[M]` the primitive holds no geometry — its one
geometric datum, ``face_areas``, was a verbatim copy of
:attr:`~orpheus.mesh.structured.Mesh1D.areas`, already single-sourced in
:func:`~orpheus.geometry.coord.compute_areas_1d`, while ``delta_A`` has
zero non-S\ :sub:`N` consumers and
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` carries a
:math:`\Delta A` divided by a *quadrature weight*.  It was also an
island in its own package: every genuine geometry primitive had 1-4
intra-``geometry/`` consumers and this one had **0**.  The layer test it
failed — *a datum belongs to the layer that can define it without naming
a method; everything else is posing, and posing belongs to the method
head that poses it* — moved it to
:mod:`orpheus.sn.mesh.reduced_operator`.

.. important:: **Which solvers is "the solvers that need it"?**  Until
   2026-08-27 this page — and the module's own docstring — read *"SN,
   MoC, and CP curvilinear sweeps all march through the same data"* and
   *"downstream consumers (SN, MoC, CP) share"* it.  That is
   **measurably false, and structurally so**: `[M]` no file under
   ``orpheus/moc/``, ``orpheus/cp/`` or ``orpheus/mc/`` names this
   primitive under **any** of the eight spellings the census at
   :ref:`connection-coefficient-census` enumerates — while both of that
   census's positive controls name every one of them.  That is not a
   migration which has yet to happen, but a term those methods never
   form.  The chart data is still correctly geometry-layer; what was
   wrong is the consumer list.  The reason is worked in
   :ref:`who-needs-a-connection-coefficient` below, and the
   re-runnable predicate — deliberately a predicate rather than a table
   of counts — is at :ref:`connection-coefficient-census`.

The :math:`\alpha` dome recursion (sphere) — Hébert (2009) §3.9.4
Eqs. 3.423-3.424, after Lathrop, K., & Carlson, B. (1966),
*J. Comp. Phys.* 1:173, in the ORPHEUS factor-of-2-absorbed
normalization:

.. math::
   :label: alpha-dome-recursion

   \alpha_{n+\tfrac12} = \alpha_{n-\tfrac12} - w_n\,\mu_n,
   \qquad \alpha_{\tfrac12} = 0.

.. vv-status: alpha-dome-recursion documented
   Rationale: this is the literature-transcribed definition of the
   :math:`\alpha` recursion (Hébert §3.9.4 Eqs. 3.423-3.424, after
   Lathrop & Carlson 1966), a representational identity rather than a
   solver claim.  ⛔ The label read ``bailey-dome-recursion`` until
   2026-08-27; that name encoded the wrong-paper attribution retracted
   at Issue #168 Phase B (see
   :ref:`sn-citation-corrections`).  The verifiable content is the
   dome-closure contract — ``tests/gates/geometry/test_reduced_operator.py``
   (``test_every_shipped_gauss_legendre_dome_closes``,
   ``test_every_shipped_folded_product_dome_closes_on_every_level``,
   and the negative control ``test_a_dome_that_does_not_close_is_refused``)
   plus ``tests/gates/sn/sweep/curvilinear/test_alpha_closed_form.py``
   (``test_production_alpha_is_a_non_negative_closing_dome``).
   ⚠ The SAME recursion is also stated on the S\ :sub:`N` methods page as
   :eq:`alpha-recursion`, which is the label the ``verifies`` markers
   target; see the de-duplication note at the end of this section.

For Gauss--Legendre :term:`quadrature` with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome
(:math:`\alpha_{1/2} = 0 \to \text{peak} \to \alpha_{N+1/2} = 0`)
that closes back to zero at the upper boundary by GL antisymmetry.
The cylindrical analog runs **per-**\ :math:`\mu`\ **-level**: each level
:math:`p` carries its own :math:`(M+1)`-tuple of
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine and :math:`M` is the
number of azimuthal ordinates on that level — **and it closes on that arm
too**, :math:`\alpha^{(p)}_{M+1/2} = 0` on every level, for the same reason
and by the same telescoping.

⭐ **Both ends are zero, and only one of them is an axiom.**  The recursion
is strictly one-sided — it is seeded at :math:`\alpha_{1/2} = 0` and never
consults the far end — so telescoping it over the level gives
:math:`\alpha_{M+1/2} = -\sum_m w_m c_m` in the level's marching cosine
:math:`c` (:math:`\mu` sphere, :math:`\eta` cylinder).  The far endpoint
therefore vanishes **iff the measure's first moment in the marching
coordinate does**, which makes it a *property of the quadrature* — a real
admission contract a bad rule can violate — rather than a property of the
recursion.  One body computes the dome
(:func:`~orpheus.sn.angular.redistribution.alpha_dome`, called by both
curvilinear factories, with the derivations-side name delegating to it) and
one guard refuses a non-closing measure
(``_assert_alpha_dome_closes``, per level on the cylinder so the offending
level is named).  Until ``bea6a367`` (2026-08-12) the contract was a bare
``assert`` on the sphere arm and *nothing* on the cylinder — and a bare
``assert`` is stripped by the canonical ``python -O`` runner, so it did not
run at all.  Full account, including why the fix had to start with
de-duplicating three copies of the recursion:
:ref:`sn-alpha-dome-closes`.

The Morel--Montry closure weight is the **barycentric coordinate of the
ordinate between the two edges of its own angular cell** — predicate
**P2**, :cite:`BaileyMorelChang2010` Eq. 43 = Lathrop 2000 Eq. 23 —
equivalently the UNIQUE closure weight exact for an angular flux affine
in the radial cosine:

.. math::
   :label: morel-montry-closure

   \tau_m
     \;=\; \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}} .

.. vv-status: morel-montry-closure documented
   Rationale: this is the literature-transcribed definition of the M-M
   angular closure weight (Bailey-Morel-Chang 2010 Eq. 43 = Lathrop 2000
   Eq. 23); it is a representational identity, not a solver claim. The
   verifiable content is the producer-equivalence gate
   ``tests/gates/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` (both
   arms now compare against HAND-AUTHORED references — the analytic arc
   closed form on the cylinder, an inline cumulative-weight expression on
   the sphere) plus the ν-closure and P3 gates in
   ``tests/gates/sn/sweep/test_tau_arc_wellposedness.py``.  τ is closure-owned,
   NOT a reduced-operator field — see the τ-ownership note below.

⭐ **There is no "raw" and no "clamped" τ** (Q5.6.4, 2026-08-11): the
:math:`[\tfrac12, 1]` absorber RETIRED, ``morel_montry_tau_raw_per_level``
retired with it, and :eq:`morel-montry-closure` carries **no geometry**.
One generic body serves both arms — the geometry lives entirely in the
**cell partition**, which is a separate object with its own producer.

.. _angular-cell-partition-section:

The angular cell partition — where the geometry actually lives
--------------------------------------------------------------

A level's ordinates each own an angular *cell*; the partition is the
:math:`(M+1)` cell edges in the radial direction cosine
(:math:`\mu` sphere, :math:`\eta` cylinder), and it is the object
:eq:`morel-montry-closure` reads:

.. math::
   :label: angular-cell-partition

   \mu_{m+1/2} \;=\;
   \begin{cases}
     \mu_{m-1/2} + w_m, \qquad \mu_{1/2} = -1
       & \text{sphere: cumulative WEIGHT} \\[8pt]
     \sin\theta\,\cos\omega_{m+1/2},\quad
       \omega_{m+1/2} = \tfrac12\bigl(\omega_m + \omega_{m+1}\bigr),
       \;\; \omega_{1/2} = \pi,\;\; \omega_{M+1/2} = 0
       & \text{cylinder: MIDPOINT in }\omega
   \end{cases}

.. vv-status: angular-cell-partition documented
   Rationale: a geometry-of-the-rule construction, not a physics-equation
   claim with an L0..L3 ladder slot — the partition is a property of the
   quadrature, produced solve-free.  The verifiable content is
   ``tests/gates/sn/sweep/test_angular_cell_partition.py`` — the **direct**
   value gate on the producer, both arms, added 2026-08-11: a
   hand-written cumulative-weight reference (sphere) and the analytic
   equispaced-arc closed form :math:`e_k = \sin\theta\cos(\pi - k\Delta
   \omega)` (cylinder), each with a negative control (the uniform
   partition; the retired chord partition) per ``vv-principles`` #19,
   plus the closing identities, the march-orientation sign law, and a
   labelled control recording that :math:`M = 2` — i.e. every
   ``folded_product(·, 4)`` fixture — is structurally BLIND to the
   partition choice.  Then
   ``tests/gates/sn/sweep/curvilinear/test_tau_producer_equivalence.py``
   (:math:`\tau` = P2 applied to the partition, same two references),
   ``tests/gates/sn/sweep/test_tau_arc_wellposedness.py`` (the P3 theorem and
   the attainable closed endpoint) and
   ``tests/gates/sn/verification/mms/test_mms_ordering_blindness.py``
   ``::test_the_full_circle_double_cover_is_REFUSED_by_the_cell_partition``
   (the non-monotone-arc refusal).  All are ``foundation`` gates —
   software/structural invariants of a discrete construction.

   ⚠ Until 2026-08-11 the partition producer had **no value gate at
   all**: every listed test read :math:`\tau`, which is P2 *applied* to
   the partition, so a wrong partition that kept :math:`\tau` inside
   :math:`[0,1]` was visible only to the two cylinder :math:`\tau` rows.
   The recurrence's worst partial amplification
   :math:`A(M) = \max_m \prod_{k\le m}(1-\tau_k)/\tau_k` — the number
   quoted as "recurrence error-amplification" in the Q5.6.4
   adjudication — is likewise now committed, in
   ``tests/gates/sn/sweep/curvilinear/test_psi_half_positivity.py``.

Both branches are **derived, not conventional**, and they are derived
from *different* facts:

* **Sphere** — a Gauss--Legendre weight *is* the cell's
  :math:`\mu`-measure, so accumulating weights from :math:`\mu_{1/2} = -1`
  partitions :math:`[-1, 1]` exactly.  This is
  :cite:`BaileyMorelChang2010` Eq. (12) **verbatim**, corroborated
  independently by Lathrop 2000 p. 249 (:math:`\sum \Delta\mu_m = 2`).
  Unchanged by Q5.6.4.
* **Cylinder** — the azimuthal march is a march in :math:`\omega`,
  **arc by arc**, so the cell boundary is the midpoint *in the variable
  the march marches in*.  Taking the midpoint of the **stored**
  :math:`\omega` values keeps it correct for any monotone arc; on an
  equispaced-:math:`\omega` rule it is exactly the half-angle boundary
  :math:`\omega_m \pm \Delta\omega/2`.

Feeding the equispaced-:math:`\omega` case of
:eq:`angular-cell-partition` through :eq:`morel-montry-closure` gives a
closed form, which is the cylinder arm's **structurally independent**
reference in the producer-equivalence gate (it shares no code path with
the producer):

.. math::

   \tau_m \;=\; \tfrac12
     \;+\; \tfrac12\,\cot\omega_m\,\tan\!\bigl(\Delta\omega/4\bigr) .

`[M]` (2026-08-11) agreement with the producer on
``folded_product(n_mu=4, n_phi)``, maximised over all four levels:
:math:`1.1\mathrm{e}{-16}` / :math:`2.2\mathrm{e}{-16}` /
:math:`7.8\mathrm{e}{-16}` / :math:`7.4\mathrm{e}{-15}` /
:math:`2.3\mathrm{e}{-14}` at :math:`n_\varphi = 4/8/16/32/64`.  ⚠ It
**degrades with refinement** (:math:`\arctan2`/:math:`\cos` round-off in
both paths), so the gate asserts ``atol=1e-13`` rather than a
machine-epsilon bound; a row beyond :math:`n_\varphi = 64` must widen it
(``vv-principles`` #16 — never assert tighter than the producer
achieves).

.. warning:: **Do NOT "unify"** :math:`\alpha` **and** :math:`\tau`
   **onto one expression.**

   Both reference the same partition — that is the point of
   :eq:`angular-cell-partition`, and deriving it twice is exactly the
   defect Q5.6.4 fixed.  But they impose **two different conditions** on
   it: :math:`\tau` the **zeroth** moment (P2, above), :math:`\alpha` the
   **first**-moment conservation recursion (P4,
   :eq:`alpha-dome-recursion`; Hébert 3.397--3.399, after Alcouffe &
   O'Dell).  Forcing :math:`\alpha` to equal the geometric tangential
   cosine at these edges silently drives Lathrop's defect
   :math:`\delta \to 0`, i.e. :math:`\tau \to \tfrac12` — the angular
   *diamond* scheme (Hébert 3.406/3.431), a **different method** with a
   different diffusion limit.  ⚠ Hébert's own
   :math:`\eta_{p,q\pm1/2}` is a constant-flux conservation recursion,
   **not** a trig evaluation at a bisected :math:`\omega`; the closed
   form above is a theorem about *our* equispaced-:math:`\omega` rule,
   not the literature's definition of the partition.

.. _tau-ownership-note:

.. note:: **Where τ lives.**

   :math:`\tau_m` is **not** a
   :class:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator`
   field.  The geometry-layer primitive carries the SPATIAL curvature
   coefficients (``face_areas``, ``delta_A``), a reference to the
   ANGULAR factor
   (:class:`~orpheus.sn.angular.redistribution.AngularRedistribution`,
   which owns the :math:`\alpha`-dome and :math:`\mu_{\rm start}` as of
   the 2026-08-26 un-weld) and — since the P4-remainder, 2026-08-29 —
   the generator-stamped angular SPACE factor
   :attr:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.angular_axis`,
   which is how it reaches the quadrature at all.  Three fields, three
   different objects: the spatial chart, the :math:`\alpha` data, and
   the space.  The M-M closure weight is produced by
   :func:`~orpheus.sn.angular.closure.morel_montry_tau_per_level`
   reading the single partition producer
   :func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`.
   ⛔ This note said that producer was called *"by* :class:`SNProblem`
   *against the quadrature and its own* ``self.coord``\ *"* until
   2026-08-29.  `[M]` its one production caller is
   :class:`~orpheus.sn.angular.closure.MorelMontryAngularSweep`'s
   constructor, against ``angular.coord`` and the quadrature it recovers
   from its ``angular_axis`` operand — which is the same ruling stated
   correctly (τ is the closure's, not the hub's), and is why the closure
   needed the space factor handed to it in the first place.
   The split is deliberate: τ is a property of
   the *angular closure scheme*, selectable per mesh, while the curvature
   coefficients are a property of the *geometry*.  Statements elsewhere
   in the corpus describing ``tau_mm`` / ``tau_mm_per_level`` as factory
   outputs predate that move, and any naming
   ``morel_montry_tau_raw_per_level`` predate Q5.6.4.

The admissible range is :math:`\tau \in [0, 1]` — predicate **P3**, an
ordinate lies inside its own angular cell — **enforced since Q5.5
(2026-08-07)**: the producer RAISES on :math:`\tau \notin [0, 1]`.  On a
well-posed monotone march membership certifies the march; a value outside
certifies an ILL-POSED march.  Both realized cases were caught by
measurement: (a) mis-ordered members — T22's ω-ordered mis-ordering
measured :math:`\tau = 1.079` at the producer before surfacing as a NaN
400 lines downstream, the pre-Q5.5 absorption silently laundering it into
a finite wrong answer; (b) a quadrature incompatible with the arm's edge
convention — a raw 3-D ``level_symmetric(4)`` rule fed to the 1-D
spherical arm (24 unsorted ``mu_x`` with duplicates, weights summing to
:math:`4\pi`) measured :math:`\tau \in [-20.3,\, 1.13]` with 23 of 24
ordinates outside, consumed SILENTLY by the *unclamped* sphere closure
until the guard landed — seven operator-equivalence tests ran exactly
this configuration and stayed green because both compared spellings share
the :math:`\tau` (the Mode-12 annihilation).  Issue #336 tracks the
refuse-or-reduce design for ``SNProblem`` on a spherical mesh with a
non-μ-line rule.  The closed endpoints are legal march starts — ``0`` is
an edge-node start and ``1`` an η-degenerate tie
(:func:`~orpheus.sn.angular.closure.march_start_structure_per_level`).
The guard does NOT catch the double cover: a full-circle level's
:math:`[0, 1, 0, 1, \ldots]` fingerprint is entirely inside
:math:`[0, 1]`; that detector is the singular set :math:`\Sigma`, and
since Q5.6.4 the non-monotone arc is refused one frame earlier still, by
:eq:`angular-cell-partition`'s own producer (a full-circle level carries
:math:`\omega` of both signs, so "the midpoint in :math:`\omega`" is not
defined for it).

⭐ **On the cylinder, P3 is now a THEOREM.**  With the partition taken as
the ω-midpoint, a strictly ω-monotone level has
:math:`\omega_{m-1/2} > \omega_m > \omega_{m+1/2}` by construction, and
:math:`\cos` is monotone on :math:`(0, \pi)`, so
:math:`\eta_{m-1/2} < \eta_m < \eta_{m+1/2}` — :math:`\tau \in (0, 1)` is
**forced** (`[M]` 4000 random monotone arcs: :math:`\min\tau =
4.739\mathrm{e}{-7}`, :math:`\min(1-\tau) = 7.599\mathrm{e}{-10}`, never
outside).  Its only equality case is a node ON an arc endpoint, i.e. a
node on :math:`\Sigma` — so **cylinder-P3 reduces to the fold criterion**
:math:`\Sigma = \emptyset`.  P3 keeps its teeth on the **sphere**, where
cumulative-weight edges genuinely need not bracket their nodes (case (b)
above).  Stated plainly so no audit reads cylinder-P3 as live coverage it
is not.

.. _sn-tau-absorber-retirement:

The retired :math:`[\tfrac12, 1]` absorber — what it was compensating for
-------------------------------------------------------------------------

**The clamp was TWO objects welded together** (T27, adjudicated
2026-08-02; membership guard landed Q5.5, 2026-08-07; absorber retired
Q5.6.4, 2026-08-11).  The cylinder-only expression
``max(0.5, min(1.0, τ_raw))`` fused the :math:`[0, 1]` *membership* — now
the P3 guard above — to a :math:`[\tfrac12, 1]` *absorption* whose stated
purpose was blocking an edge-node division.  The
:math:`[\tfrac12, 1]` box was **never** the admissible range of
:math:`\tau`: the sphere ran outside it, unclamped and correct, and `[M]`
at :math:`S_8` Gauss--Legendre **four of eight** M-M τ sit below
:math:`\tfrac12`.  :cite:`BaileyMorelChang2010`'s own :math:`S_2` example
gives :math:`\tau_1 = \mu_1 + 1 = 1 - 1/\sqrt3 \approx 0.4226 < \tfrac12`
(their Eq. 47).  **No source prescribes any limiter on** :math:`\tau`.

.. _sn-tau-absorber-provenance:

Where the number :math:`[\tfrac12, 1]` actually came from
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"No source prescribes it" is the weaker statement.  The stronger one,
established 2026-08-11: the interval is real, it is **Grant's**, and it
is on a **different parameter**.

:cite:`ReedLathrop1970` carry two independent weighted-diamond
parameters — a **spatial** weight :math:`a_{i+1/2}` (their Eqs. 5a/5b
and 11) and an **angular** weight :math:`\tau_m` (their Eqs. 6a/6b).
Their footnote 8 (printed p. 239) discusses Grant's choice for the
*spatial* one: Grant lets that weight depend on the sign of its
argument, and Reed & Lathrop note this "is necessary only to keep
[it] between :math:`\tfrac12` and 1" — then add, decisively, that
**Grant does not determine angular weights at all**.

So :math:`[\tfrac12, 1]` is a bound on the SPATIAL weighted-diamond
parameter of Grant, I. P. (1968), *J. Comp. Phys.* 2(4):381-402,
doi:10.1016/0021-9991(68)90044-2.  ⚠ **Grant 1968 is not in the local
library and has not been read**; the attribution rests on Reed &
Lathrop's footnote, read on the rendered page.  Transplanting that
interval onto the angular :math:`\tau` is exactly what the retired
cylinder absorber did — the number was inherited across a parameter
boundary, which is why no search of the angular literature could ever
find its source.

⭐ **What it was actually compensating for was a WRONG PARTITION, and
that is why retiring it alone made things worse.**  Until 2026-08-11 the
cylinder edges were taken at the midpoint of consecutive :math:`\eta`
values — the **chord** midpoint — with the endpoints pinned at
:math:`\mp\sin\theta`.  That partition is :eq:`angular-cell-partition`
*with its end cells stretched*.  The identity

.. math::

   \tfrac12\bigl[\cos\omega_a + \cos\omega_b\bigr]
     \;=\; \cos\!\Bigl(\tfrac{\omega_a+\omega_b}{2}\Bigr)
            \cos\!\Bigl(\tfrac{\Delta\omega}{2}\Bigr)

(the :math:`\kappa` prefactor's sibling) makes every *interior* chord
edge exactly :math:`\cos(\Delta\omega/2) \times` the arc edge (`[M]`
agreement :math:`10^{-16}`), while the two endpoints stay **unscaled** —
so the outermost cells stretch to absorb the shrink.  The :math:`\eta`
error vanishes as :math:`\Delta\omega \to 0`, but the implied
:math:`\omega`-width spread does **not**: it converges to
:math:`\approx 17.45\,\%` (`[M]` 18.71 / 17.59 / 17.48 / 17.46 % at
:math:`n_\varphi = 8/16/32/64`) against a quadrature whose own cells are
bit-exactly equal.  That :math:`O(1)` inconsistency — one object, the
"boundary between azimuthal cell :math:`m` and :math:`m+1`", derived
independently by :math:`\alpha` (at the real half-angle
:math:`\omega_{m-1/2}`) and by :math:`\tau` (at the chord midpoint), in
disagreement — is what the absorber hid.

**The absorber is condemned on its own terms, with no MMS involved.**  The
ν-closure diagnostic marches the level *implied by* :math:`\tau`
(:math:`\nu_{1/2} = -\sin\theta`,
:math:`\nu_{m+1/2} = (\eta_m - (1-\tau_m)\nu_{m-1/2})/\tau_m`) and asks
whether it lands on :math:`+\sin\theta`.  It is solve-free and it
separates a derived τ from a fabricated one — `[M]`
:math:`\nu/\sin\theta` at close:

.. list-table:: ν-closure: does the march implied BY τ close the level?
   :header-rows: 1
   :widths: 16 21 21 21 21

   * - :math:`n_\varphi`
     - arc ω (production)
     - chord (retired)
     - clamped (retired)
     - :math:`\tau \equiv \tfrac12`
   * - 8
     - ``1.000000``
     - ``1.000000``
     - **``1.016389``**
     - **``1.164784``**
   * - 16
     - ``1.000000``
     - ``1.000000``
     - ``1.001930``
     - ``1.039182``
   * - 32
     - ``1.000000``
     - ``1.000000``
     - ``1.000238``
     - ``1.009677``
   * - 64
     - ``1.000000``
     - ``1.000000``
     - ``1.000030``
     - ``1.002412``

⟹ the clamped values and :math:`\tau \equiv \tfrac12` correspond to **no
partition of the level at all** — their implied march overshoots the
level's own endpoint by 1.6 % and 16.5 % respectively.  (:math:`\tau
\equiv \tfrac12` is Hébert's angular *diamond* scheme; it is listed here
because it is the tempting rescue, and because it optimises truncation
order while breaking the diffusion limit :math:`\tau` exists to fix.)

**Why the sphere's convention cannot simply be transplanted.**
Accumulating weights in :math:`\eta` — even correctly renormalised to
:math:`\sum \bar w = 2\sin\theta`, :cite:`BaileyMorelChang2010` Eq. 52 —
violates **P3** on our rule and **worsens with refinement**: `[M]`
ordinates outside their own cell go 0/4 → 4/8 → 12/16 → 28/32 at
:math:`n_\varphi = 8/16/32/64` **per level** (`[M]` 0/16, 16/32, 48/64,
112/128 over the four levels of ``folded_product(n_mu=4, n_phi)``, with
:math:`\tau` ranging out to :math:`[-2.86,\, 3.86]` at
:math:`n_\varphi = 64`), and the solve diverges (NaN) from
:math:`n_\varphi \ge 16`.  The reason is structural: an arc cell's
:math:`\eta`-measure is
:math:`2\sin\theta\,\sin\omega_m\,\sin(\Delta\omega/2) \propto
\sin\omega_m`, **not** constant — `[M]` at :math:`n_\varphi = 16` it
spans :math:`0.30`--:math:`1.53 \times` the uniform width across one
level (and :math:`0.08`--:math:`1.57\times` at
:math:`n_\varphi = 64`) — while a trapezoid weight is.  ⟹ Eq. (52) is
not a law; it is the
statement that *in their* quadrature the weight equals the cell's
:math:`\eta`-measure.  Ours does not, so we satisfy the same **predicate**
by a different partition.  (This diagnosis is not new — it was written
into the original M-M implementation as *"weights are uniform in*
:math:`\varphi` *-space, not* :math:`\eta` *-space"*, and is preserved
verbatim in
:mod:`orpheus.derivations.discrete.sn.angular_differencing`.  The
diagnosis was right; the chord-midpoint *substitute* was never checked
against a partition predicate.)

**On the σ_y-folded arc the absorption's reason is structurally gone.**
With midpoint nodes the smallest-η ordinate sits at
:math:`\omega_0 = \pi - \Delta\omega/2`, so the closed form above gives

.. math::
   :label: morel-montry-folded-arc

   \tau_0
     \;=\; \tfrac12 - \tfrac12\cot\!\bigl(\Delta\omega/2\bigr)
            \tan\!\bigl(\Delta\omega/4\bigr)
   \;\xrightarrow[\;n_\varphi \to \infty\;]{}\; \tfrac14
   \;\;\text{from inside},
   \qquad
   \tau_m \in \Bigl[\tfrac14,\, \tfrac34\Bigr],
   \qquad
   \tau_m + \tau_{M-1-m} \;=\; 1 ,

since :math:`\cot(\Delta\omega/2)\tan(\Delta\omega/4) \to \tfrac12`.
`[M]` (2026-08-11, :math:`n_\mu = 4`, folded staggered), with
:math:`|\Sigma| = 0` throughout — :math:`\tau \in \{0, 1\}` is
**structurally unreachable** on the fold:

.. list-table:: Folded-arc τ box and the reversal identity
   :header-rows: 1
   :widths: 20 45 35

   * - :math:`n_\varphi`
     - τ range
     - reversal residual
   * - 4
     - ``[0.292893, 0.707107]``
     - 0.5 ULP
   * - 8
     - ``[0.259892, 0.740108]``
     - 0.5 ULP
   * - 16
     - ``[0.252425, 0.747575]``
     - 2.0 ULP
   * - 32
     - ``[0.250603, 0.749397]``
     - 7.0 ULP
   * - 64
     - ``[0.250151, 0.749849]``
     - 12.0 ULP

.. note:: **Retraction (2026-08-11, Q5.6.4).**  Until Q5.6.4 this
   equation read :math:`\tau_{{\rm raw},0} \to \tfrac15`,
   :math:`\tau_{{\rm raw},m} \in [\tfrac15, \tfrac45]`, with the
   reversal identity holding **bit-exactly** (residual :math:`0.0` at
   every :math:`n_\varphi`; measured :math:`[0.2195, 0.7805]` at
   :math:`n_\varphi = 8` falling to :math:`[0.200289, 0.799711]` at 64).
   Those numbers are correct — *for the retired chord partition they were
   measured on*.  The box is now :math:`[\tfrac14, \tfrac34]`.

   ⚠ **And the reversal identity is no longer bit-exact: that is a trade,
   not a regression.**  The chord partition's reversal symmetry was exact
   *because* both end cells were stretched symmetrically — the 17.5 %
   ω-width defect cancelled itself under
   :math:`\omega \to \pi - \omega`.  The ω partition has the correct cells
   and pays 0.5--12 ULP of :math:`\arctan2` / :math:`\cos` round-off.  The
   gate asserts 64 ULP, because the residual grows with arc refinement and
   a bit-exact assertion would be a latent false red
   (``vv-principles`` #16 — never assert tighter than the producer
   achieves).

.. (vv-status rationale) morel-montry-folded-arc: Verified by the Q5.5
   mechanism gates in ``tests/gates/sn/sweep/test_tau_arc_wellposedness.py``,
   re-posed at Q5.6.4 —
   ``test_the_fold_mechanism_is_an_empty_singular_set`` asserts the
   MECHANISM (Σ = ∅, computed via ``singular_set_under``, never declared) and
   ``test_the_folded_tau_is_bounded_with_the_reversal_identity`` the
   CONSEQUENCE (τ ⊂ [1/4, 3/4] per level plus the reversal identity at
   64 ULP), each at n_φ ∈ {8, 16, 32, 64} on the folded staggered product
   (n_μ = 4).  Teeth measured 2026-08-07: reverting Q5.2's offset to
   δ = 0 reds BOTH legs at every n_φ (8/10 red), which attributes the
   pass to the mechanism; the [0, 1] guard's negative companion reds
   alone under a no-opped guard (1/10 red).  A geometry-of-the-rule
   invariant, not a physics-equation claim with an L0..L3 ladder slot.
.. vv-status: morel-montry-folded-arc documented

**The honest accuracy cost, ratified rather than hidden.**  The
principled partition is not uniformly more accurate on the one
manufactured fixture available: `[M]` on the anisotropic cylindrical MMS
at :math:`n_x = 320` it is BETTER at :math:`n_\varphi = 8`
(:math:`3.128\mathrm{e}{-3}` vs :math:`3.511\mathrm{e}{-3}`) and
:math:`\sim 1.8`--:math:`2\times` WORSE at :math:`n_\varphi = 16/32/64`.
Principled :math:`\ne` more accurate: a scheme satisfying P2/P3 wins over
one with a smaller number on a single MMS, and the L2 norm measures
truncation order — exactly what :math:`\tau \equiv \tfrac12` optimises
and exactly the quantity that is blind to the diffusion limit
:math:`\tau` exists to fix.  The full ladder, and why the MMS is the
wrong instrument for this decision, is at
:ref:`sn-cylinder-angular-floor` in :doc:`/theory/verification/sn`.

The predicate ladder (P0--P4) these decisions are made against, the
:math:`\tau`/:math:`\beta` nomenclature (both letters are overloaded, and
both collisions have cost real time), and a written record of **which
diagnostics are blind on which rules** live in
:mod:`orpheus.derivations.discrete.sn.angular_differencing`.

API surface
-----------

The geometry-layer primitive is built by three factory functions, one
per coordinate system:

* :func:`~orpheus.sn.mesh.reduced_operator.slab_streaming(mesh, ang)
  <orpheus.sn.mesh.reduced_operator.slab_streaming>` — Cartesian 1-D;
  no curvature math.  Both its factors are present and **neutral**: the
  ANGULAR one carries a zero dome and the diameter-ray start, and the
  SPATIAL one carries a unit cross-section with **zero area change**
  (``face_areas == ones(nx+1)``, ``delta_A == zeros(nx)``) — because a
  slab having "no curvature" IS its faces not changing area.

  ⛔ Until P4.1b (2026-08-27) this read *"its spatial arrays remain
  ``None``"*, and they were stored fields.  They are **derived** from
  the mesh now (``face_areas`` is ``mesh.areas`` itself, ``delta_A`` its
  difference), so no factory computes them and the per-coordinate
  ``Optional`` union is dead on **both** factors rather than one.  That
  is what let :meth:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.streaming_terms`
  collapse from three chart arms to one shared body: the retired
  CARTESIAN arm was the spherical body with ``1.0`` / ``1.0`` / ``0.0``
  written out by hand.
* :func:`~orpheus.sn.mesh.reduced_operator.spherical_streaming(mesh, ang)
  <orpheus.sn.mesh.reduced_operator.spherical_streaming>` — 1-D spherical
  with the dome recursion :eq:`alpha-dome-recursion` and Morel--Montry
  closure :eq:`morel-montry-closure`.
* :func:`~orpheus.sn.mesh.reduced_operator.cylindrical_streaming(mesh, ang)
  <orpheus.sn.mesh.reduced_operator.cylindrical_streaming>` — 1-D
  cylindrical with **per-**\ :math:`\mu`\ **-level** :math:`\alpha`,
  :math:`\Delta A/w` and :math:`\mu_{\rm start}` lists (τ is
  closure-owned — see the :ref:`τ-ownership note <tau-ownership-note>`).  Requires the
  angular measure to expose ``level_indices`` (e.g., a
  :class:`~orpheus.numerics.quadrature.Quadrature` built from
  :meth:`Quadrature.level_symmetric
  <orpheus.numerics.quadrature.Quadrature.level_symmetric>` or
  :meth:`Quadrature.product
  <orpheus.numerics.quadrature.Quadrature.product>`).

The per-cell, per-direction inputs needed by a sweep cell update are
extracted via
:meth:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.streaming_terms`,
which returns a
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` dataclass —
`[M]` today exactly four fields, ``face_area_inner``,
``face_area_outer``, ``volume`` and ``abs_mu``, **populated on every
chart**.

``volume`` is the per-cell volume :math:`V_i`; ``abs_mu`` is the
absolute primary direction cosine :math:`|\mu|` (sphere) /
:math:`|\eta|` (cylinder, radial) / :math:`|\mu_x|` (slab).  All four
are populated by all three factories, so a downstream sweep cell update
— see :doc:`/theory/methods/sn/index`, "Cell update strategies (the
strategy contract)" — receives a self-contained per-cell,
per-direction packet and need not reach back into ``SNProblem``.

.. warning:: ⛔ **Two claims in this paragraph were retired, and the
   second is the one that matters.**

   It read *"whose populated fields are geometry-dependent (slab is
   minimal; sphere/cylinder carry the full curvature-coefficient
   bundle)"* and *"the* ``alpha_in is None`` *test discriminates slab
   from curvilinear inside cell-update strategies"*.  Both were true of
   a packet that has since lost every field they named.  Issue #196
   Phase G Step 2.5 gave the slab *neutral* curvature rather than
   ``None``\ s; Issue #236 Step C deleted the Morel–Montry
   ``alpha_in`` / ``alpha_out`` / ``tau_mm`` fields outright (τ is
   closure-owned — see the :ref:`τ-ownership note <tau-ownership-note>`
   above); and P4.7 (2026-08-29) shed the last three, ``mu``,
   ``chord_length`` and ``delta_A_over_w``.  **A spatial scheme no
   longer discriminates slab from curvilinear at all**, because the
   curvature data reaches it already reduced to numbers whose slab
   values are the neutral element of the arithmetic they enter.

   ⚠ Note also what the surviving fields are *not*.  ``abs_mu`` is the
   **ordinate's**, not the geometry's, so the packet is not "purely
   geometric" — a reading refuted 2026-08-28.  It is the evaluation
   point of a spatial closure for a *directional* method, which is why
   it lives beside the scheme contract in
   :mod:`orpheus.transport.spatial.scheme` rather than in
   :mod:`orpheus.geometry`.  The degenerate cylindrical cell is
   signalled by the geometric ``visit.face_area_downstream == 0.0``,
   never by a numerical threshold on :math:`|\mu|`.

Geometric labels, not flow-direction labels
-------------------------------------------

The two face-area fields on
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` are
**purely geometric**: ``face_area_inner`` is :math:`A_{i-1/2}` (the
face closer to :math:`r=0`), ``face_area_outer`` is
:math:`A_{i+1/2}` (the face farther from :math:`r=0`).  These labels
are independent of the sweep's marching direction.  For an outward
sweep (centre :math:`\to` boundary) the inner face is upstream; for
an inward sweep (boundary :math:`\to` centre) the outer face is
upstream.  But that resolution is **SN-specific** — the SN sweep is
a topological sort of a directed cell graph for a given ordinate,
where edges are oriented by
:math:`\mathrm{sign}(\Omega \cdot \hat n_{\text{face}})`.  MoC uses
a different mathematical structure (fiber bundles + solution
sheaves), CP / diffusion / MC do not have a sweep at all.

Per Cardinal Rule 2, the geometry layer therefore stays geometric.
Sweep-direction resolution lives in the SN module:
:class:`~orpheus.transport.spatial.scheme.CellVisit` is the
SN-specific per-visit packet that composes the geometric
:class:`StreamingTerms` together with the sweep-resolved
``face_area_downstream``.  The SN sweep iterates
:meth:`~orpheus.sn.problem.SNProblem.dag_walk`, which encodes
the inward / outward branching, the cylindrical per-level
traversal, and the pure-azimuthal degenerate handling — yielding
one :class:`CellVisit` per cell in DAG-topological order.  The
cell-update strategy then sees only resolved data; no
sign-of-:math:`\mu` branching inside the strategy.

Likewise, the signed primary direction cosine ``mu`` is read from
the **global ordinate index** for all three coordinate systems:
slab and sphere have ``direction_idx`` :math:`=` global ordinate;
cylindrical resolves the global index through
``level_indices[mu_level_idx][direction_idx]`` because cylindrical
``direction_idx`` is the within-level azimuthal index
:math:`m \in [0, M)`.  ``abs_mu`` follows the same convention.

Bit-identical contract — and what it pins TODAY
-----------------------------------------------

When the lift landed, the factories were required to produce arrays
bit-identical to the historical inline implementations
``SNProblem._setup_spherical`` and ``SNProblem._setup_cylindrical``.  Hash
equality — :func:`numpy.array_equal`, never ``np.allclose`` — was
enforced at test time by ``tests/gates/geometry/test_reduced_operator.py``
(``foundation``-tagged), so the two paths had to share every
floating-point bit.  That is what made the lift safe at the time: the
then-consumers (the SN sweep and the curvilinear Krylov operator) were
unaffected because the two paths computed the same data.

.. warning:: **That gate has since been DEMOTED — read its green
   accordingly.**

   The two legacy setup methods no longer exist (see
   :ref:`snmesh-as-router` below).  ``SNProblem.__init__`` now calls
   :func:`~orpheus.sn.mesh.reduced_operator.spherical_streaming` /
   :func:`~orpheus.sn.mesh.reduced_operator.cylindrical_streaming`
   itself, so the surviving hash-equality legs compare a fresh factory
   call against ``problem.reduced`` — *the value that same factory
   produced*, routed through the Problem's constructor.

   ⛔ This paragraph used to add *"and the two ``SNProblem.face_areas`` /
   ``SNProblem.delta_A`` legs are deprecated read-throughs to that same
   object"*.  Those accessors **retired at P4.1c** (2026-08-27) — `[M]`
   11 readers, **0** of them in ``orpheus/``, and every one of the
   tests that read them existed to verify the shims themselves.  The
   legs now read ``problem.reduced.*`` directly.  The gate therefore pins the
   **wiring** (the constructor really does route to the geometry-layer
   primitive, for every geometry and every quadrature order in its
   parametrization), not the **math**: no structurally-independent
   second implementation remains on the other side of the comparison.

   Measured 2026-08-03: garbaging every array the factories emit
   leaves **all 47** tests in that file green.

   The mathematical content is pinned elsewhere.  These are the gates
   to cite for a correctness claim — identified by that same mutation,
   each one **structurally independent** (a closed form), with the SN
   curvilinear regression snapshots
   (``tests/gates/sn/regression/test_dd_regression.py``) corroborating but
   nowhere the sole evidence:

   * ``delta_A`` — the closed-form L0 term check
     ``TestL0TermVerification::test_delta_A_magnitude`` in
     ``tests/gates/sn/primitives/test_quadrature.py``, against
     :math:`4\pi\,\Delta(r^2)` / :math:`2\pi\,\Delta r`.
     ⛔ This entry read "**sole catcher**; the snapshots are blind here,
     and correctly so — ``delta_A`` has no production consumer."  True
     until 2026-08-26, **false now**: retiring the fused ``redist_dAw``
     cache made ``delta_A`` the spatial factor that BOTH
     ``streaming_terms`` and the angular closure read, so every
     curvilinear snapshot rides on it.
   * ``angular.alpha_per_level`` (was ``alpha_half`` /
     ``alpha_per_level``, until the 2026-08-26 un-weld moved the dome to
     the angular factor) — the L0 per-ordinate flat-flux identity
     ``test_per_ordinate_flat_flux_consistency`` on both arms
     (``catches("ERR-006", "ERR-007")``);
     ``tests/gates/sn/sweep/curvilinear/test_alpha_closed_form.py`` (the
     Dirichlet-kernel closed form; **cylindrical α only** — every
     fixture there is ``CoordSystem.CYLINDRICAL``); plus both snapshot
     families.
   * ``redist_dAw`` / ``redist_dAw_per_level`` — **RETIRED 2026-08-26**
     as a fused product neither of its two consumers owned.  Its catcher
     — ``tests/gates/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py``,
     the L0 closed-form :math:`\varphi = Q/(\Sigma_t(1-c))` identity —
     still covers the QUANTITY, now formed at each consumer from
     ``delta_A`` and the weights.  ⭐ And the historical note that the
     flat-flux identity "recomputes :math:`\Delta A / w` rather than
     reading the production array" is exactly why that gate needed no
     migration: it was already forming the product, not reading the
     cache.
   * ``face_areas`` — ``tests/gates/geometry/test_geometry.py`` pins the
     producer :func:`~orpheus.geometry.coord.compute_areas_1d`
     against its closed form; the snapshots pin the forwarding.

   ``tests/gates/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` is
   **not** among them, despite an earlier revision of this warning
   naming it.  #236 Step C moved :math:`\tau` to the angular closure
   (see the :ref:`τ-ownership note <tau-ownership-note>` above), which
   derives it from :math:`(\mu, w)` alone — so that gate passes
   untouched (5 passed, 0.03 s) under fully-garbaged factories.  It
   remains the right gate to cite for :math:`\tau` **itself**; it is
   simply blind to the reduced-operator arrays.

   None of these is evidence that the *lift* preserved bits, which is
   now unfalsifiable by construction.

The forward-looking half of the original rationale still holds
unchanged: new consumers bind to the geometry-layer primitive instead
of duplicating the curvature math.

.. _snmesh-as-router:

SNProblem as router
-------------------

After Round 1.1 of Wave D of the SN reshape campaign, :class:`SNProblem`
**routes** to :class:`ReducedStreamingOperator` rather than computing
the connection coefficients itself.  The :meth:`SNProblem.__init__`
ladder calls :func:`slab_streaming` / :func:`spherical_streaming` /
:func:`cylindrical_streaming` directly; the historical
``SNProblem._setup_spherical`` and ``SNProblem._setup_cylindrical`` methods
no longer exist.  ``self.reduced`` is the new canonical accessor
every downstream consumer should bind to::

    problem.reduced.streaming_terms(cell_idx, dir_idx, mu_level_idx)

returns the per-(cell, direction) packet a sweep cell update needs —
no more reaching into ``SNProblem`` for a half-dozen separate arrays.

That migration has since **completed**.  Of the eight legacy attribute
names the lift originally re-exposed as :class:`DeprecationWarning`
``@property`` accessors, only **two** survive on :class:`SNProblem` today
— ``face_areas`` and ``delta_A``, still read-throughs to the matching
field on ``self.reduced``.  The other six (``alpha_half``,
``redist_dAw``, ``alpha_per_level``, ``redist_dAw_per_level``,
``tau_mm``, ``tau_mm_per_level``) are gone — and since 2026-08-26 the
first four have no ``self.reduced`` field left to route to either, the
α-dome having moved to the angular factor and ``redist_dAw`` having
retired as a fused product.  Consumers bind to
``streaming_terms(...)`` or to ``problem.reduced.*`` directly, and the
two ``tau_mm`` names have no ``self.reduced`` field left to route to at
all — τ is closure-owned now, not a factory output (see the
:ref:`τ-ownership note <tau-ownership-note>` above).

The Cartesian path is unchanged: ``SNProblem._setup_cartesian`` still
populates the :math:`2|\mu|/\Delta x` and :math:`2|\mu_y|/\Delta y`
streaming stencils used by the DD-denominator precomputation in the
Cartesian sweep (these are SN-specific and not represented in
:class:`ReducedStreamingOperator`).  Slab geometry additionally gets
a slab :class:`ReducedStreamingOperator` for completeness so
``problem.reduced`` is always populated.

.. _who-needs-a-connection-coefficient:

Who needs a connection coefficient — and who does not
-----------------------------------------------------

The primitive above is **structurally S**\ :sub:`N`\ **-only**, and that
is a statement about the mathematics, not about how far a migration has
got.  The distinction matters because the two readings license opposite
work: *"MoC and CP have not migrated yet"* invites someone to go and
wire them up, while *"MoC and CP never form this term"* says the
primitive is correctly placed and correctly consumed by exactly one
solver family — S\ :sub:`N`.  (`[M]` the whole
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` /
:class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` chain is
referenced only from ``orpheus/sn/`` and from inside
``orpheus/transport/`` itself, and by **no** file under
``orpheus/moc/``, ``orpheus/cp/``, ``orpheus/mc/`` or
``orpheus/diffusion/``.  As above, the count is deliberately not frozen
here — re-run the predicate at
:ref:`connection-coefficient-census` with
``r"transport\.spatial|DiscretizationScheme|StreamingTerms|CellVisit"``
substituted for its pattern set; the two controls stay non-zero and the
four subjects stay at zero.)

A solver family needs an :math:`\alpha` dome **iff all three** of the
following hold.

1. **It carries an angular unknown** — a :math:`\psi` that survives
   discretisation still wearing a direction index.
2. **That index is read in a local, rotating frame** — a basis that
   turns as the spatial point moves, so a particle streaming in a fixed
   *physical* direction changes its *coordinate* label as it travels.
3. **The resulting angular derivative is discretised by collocation**
   on that index — :math:`\partial_\mu` (sphere) or
   :math:`\partial_\varphi` (cylinder) approximated as a difference
   between neighbouring ordinates, rather than by an expansion in a
   basis that differentiates exactly.

Condition 2 is what mints the term at all.  In a curvilinear chart the
direction cosines are measured against the *local* radial and azimuthal
axes, so straight-line streaming is a continuous relabelling of
:math:`(\mu, \varphi)`; the redistribution term
:math:`(1-\mu^2)/r\,\partial_\mu` (sphere) or
:math:`-(1/r)\,\partial_\varphi(\xi\,\cdot)` (cylinder) is the
bookkeeping for that relabelling.  Condition 3 is what turns its weight
into a *recursion*: collocation supplies no exact derivative, so the
weights must instead be built to preserve the one property that has to
survive — a spatially flat angular flux must redistribute to zero, per
ordinate — and :eq:`alpha-dome-recursion` together with its closure
contract (:ref:`sn-alpha-dome-closes`) is precisely the construction
that delivers it.

Adjudicating the shipped families against those three conditions:

.. list-table:: Which families satisfy the three conditions
   :header-rows: 1
   :widths: 15 19 22 20 24

   * - Family
     - (1) angular unknown
     - (2) local rotating frame
     - (3) collocated angular derivative
     - Consequence
   * - S\ :sub:`N`, curvilinear
     - yes — :math:`\psi_{n,i}`
     - yes — :math:`(\eta, \xi, \mu)` on the local radial frame
     - yes — the half-angle recursion
     - **needs the dome**
   * - MoC
     - yes — :math:`\psi` per track
     - **no** — :math:`\Omega` is fixed in the GLOBAL frame
     - n/a
     - term relocates into track segmentation
   * - CP
     - **no** — angle is integrated out before discretisation
     - n/a
     - n/a
     - term never appears
   * - MC
     - **no** — directions are sampled, not indexed
     - n/a
     - n/a
     - term never appears

**MoC fails condition 2.**  The method of characteristics is *defined*
by choosing the global frame in which :math:`\Omega` is constant along a
track, so :math:`\Omega \cdot \nabla \psi = \mathrm{d}\psi/\mathrm{d}s`
is chart-free and there is no angular derivative left to discretise.
Curvature does not disappear — it moves into *segmentation*, the
ray-region intersection that produces the chord lengths.  `[M]` the
shipped inner loop in :mod:`orpheus.moc.core` forms
:math:`\tau = \Sigma_t \, \ell_{\rm seg} / \sin\theta_p` per segment and
applies plain exponential attenuation
:math:`\Delta\psi = (\psi - Q/\Sigma_t)\,(1 - e^{-\tau})`; no ordinate
couples to its neighbour anywhere in the sweep.  This is also why
:mod:`orpheus.sn.loss_representation.sweep_graph` records that MoC will
define a *per-ray traversal* analog rather than reuse the
S\ :sub:`N` sweep graph.

**CP fails condition 1.**  Collision probability integrates the angular
variable analytically *before* anything is discretised — the transport
kernel is already an angle-integrated function of optical path, so no
angular unknown, and therefore no angular index, ever exists.  `[M]`
:class:`orpheus.cp.solver.CPSolver` dispatches on
:class:`~orpheus.geometry.coord.CoordSystem` to three real setups —
slab, cylinder, and **sphere**, i.e. the curvilinear cases the false
claim was about — and each installs a scalar kernel:
:math:`F(\tau) = e^{-\tau}` with a :math:`y`-weighted quadrature on the
sphere, the :math:`\mathrm{Ki}_3` kernel on the cylinder, :math:`E_3`
on the slab.  There is no :math:`\alpha`, no :math:`\Delta A / w`, and
nothing for either to act on.

**MC fails condition 1 as well, for a different reason.**  Monte Carlo
samples directions from a continuous distribution rather than indexing
a fixed set, so there is no neighbouring ordinate to redistribute *to*.
`[M]` :class:`orpheus.mc.solver.MCMesh` admits ``CARTESIAN`` or
``CYLINDRICAL`` — so it, too, solves a curvilinear problem with no
:math:`\alpha` anywhere.

⭐ **The curvilinear counter-examples are the load-bearing evidence.**
"MoC and CP have not migrated yet" predicts that neither has a
curvilinear capability to migrate.  Both do.  `[M]`
:class:`orpheus.moc.geometry.MOCMesh` wraps a **cylindrical**
:class:`~orpheus.mesh.structured.Mesh1D` and ray-traces **concentric
annuli** (``_ray_circle_intersections``); :class:`orpheus.cp.solver.CPSolver`
ships a real **sphere**; and, as a third witness,
:class:`orpheus.mc.solver.MCMesh` ships a real **cylinder**.  All three
solve curved geometry, and none carries one line of redistribution
machinery.  A capability that exists *and* declines the primitive
refutes the migration reading in a way that an absent capability never
could — which is exactly why the claim survived unchallenged for so
long: nobody had looked at what those packages already do.

.. _connection-coefficient-census:

Reproducing the census — the predicate, not a table of counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The census behind the `[M]` claims above is published as the **recipe,
not as a table of counts**.  Eight independent spellings of the concept
— three symbol spellings, three concept spellings, an attribute-access
spelling and the prose paraphrase — counted as regex *occurrences* over
``*.py`` under each package root:

.. code-block:: python

   import pathlib, re

   PATTERNS = {
       "reduced_operator":         r"reduced_operator",
       "ReducedStreamingOperator": r"\bReducedStreamingOperator\b",
       "AngularRedistribution":    r"\bAngularRedistribution\b",
       "alpha":                    r"\balpha\b",
       "delta_A / face_areas":     r"\b(?:delta_A|face_areas)\b",
       "redistribut*":             r"\bredistribut\w*",
       ".reduced":                 r"\.reduced\b",
       "connection coefficient":   r"connection[ -]coefficient",
   }
   SUBJECTS = ["orpheus/moc", "orpheus/cp", "orpheus/mc"]
   CONTROLS = ["orpheus/sn", "orpheus/geometry"]

   def count(pkg, pat):
       return sum(len(re.findall(pat, f.read_text(encoding="utf-8")))
                  for f in pathlib.Path(pkg).rglob("*.py"))

   for name, pat in PATTERNS.items():
       # POSITIVE CONTROL: a zero here means the filter is broken,
       # not that the tree is clean.
       assert all(count(p, pat) > 0 for p in CONTROLS), name
       # THE FINDING:
       assert all(count(p, pat) == 0 for p in SUBJECTS), name

**Run as written on 2026-08-27 it passes: every one of the eight
patterns is non-zero in both controls and exactly zero in all three
subjects.**  The zeros are the finding and they are falsifiable — the
day one of them stops being zero, the claim on this page is refuted and
should be re-argued rather than patched.

.. note:: **Why the control COUNTS are deliberately not printed here.**
   A positive control has to be **non-zero**; its particular value
   carries no part of the argument.  Freezing it would put a number in
   the corpus that moves under any edit to ``orpheus/sn`` or
   ``orpheus/geometry`` — including the ⛔ correction blocks this very
   pass added, which name the module and so raise several of the control
   counts.  ⛔ An earlier revision of this section did print them, and
   did so wrongly in three independent ways at once: the counts were
   taken **before** this pass's own edits; the column set was a
   **different** partition of "six spellings" from the one the prose
   beside it named (so two of the six numbers belonged to spellings the
   prose never mentioned); and the ``redistribut`` column was
   case-insensitive and unanchored, which silently absorbed every
   ``AngularRedistribution``.  All three are the failure mode this page
   exists to document, so the section now points at a re-runnable
   predicate instead (`plan-authoring` §9: never copy a number the tree
   can re-measure; `plan-authoring` §2: a number without its predicate
   is not re-runnable).

The set of files that name the primitive at all is small enough to
enumerate, and an enumeration — unlike a count — can be checked by
reading it.  `[M]` 2026-08-28, re-run after P4.2 + P4.3 completed the
un-weld (an earlier revision of this list, re-run after P4.4 only, counted
**15** and named ``geometry/reduced_operator.py`` a definer): **14** files
under ``orpheus/`` name any of ``ReducedStreamingOperator``, the three
``*_streaming`` factories, ``StreamingTerms``, ``AngularRedistribution``,
``angular_redistribution``, ``alpha_dome`` or ``AngularMeasure``.  Three
are the definers (``sn/mesh/reduced_operator.py`` — the connection
operator and factories; ``transport/spatial/scheme.py`` —
``StreamingTerms``, beside the contract that consumes it;
``sn/angular/redistribution.py`` — the :math:`\alpha` cluster and
``AngularMeasure``).  The rest: in ``orpheus/sn/``,
``angular/__init__.py``, ``angular/closure.py``, ``problem.py``
(``mesh/augmented_mesh.py`` until #412, 2026-09-18),
``operators/radial_characteristic.py``, ``solver.py`` and
``sweep/cache.py``; in ``orpheus/transport/spatial/``, ``__init__.py``,
``cell_balance.py``, ``diamond.py`` and ``linear_discontinuous.py``; and
``orpheus/derivations/discrete/sn/angular_differencing.py``.  Not one of
them is under ``orpheus/moc/``, ``orpheus/cp/`` or ``orpheus/mc/`` — the
load-bearing half, and it is unchanged.  ⭐ And for the first time not
one is under ``orpheus/geometry/`` either: the un-weld's own done-when,
now a property of this enumeration.

⚠ ``transport/spatial/linear_discontinuous.py`` is **new to this list and
was missing from it before P4.4**, not added by it: `[M]` it named
``StreamingTerms`` at the previous commit too.  An enumeration is only
checkable by re-running its own predicate, which is how the gap was
found — see `plan-authoring` §2 on a universal owing its denominator.

.. note:: **What WOULD change this answer.**  The three conditions are
   the claim, so the honest way to falsify it is to break one.  A
   discrete-ordinates scheme that expanded the angular flux in a basis
   which differentiates exactly — spherical harmonics, or a
   discontinuous-Galerkin / finite-element discretisation *in angle* —
   would satisfy 1 and 2 and fail 3, and would then need a different
   object entirely (a mass/stiffness pair in :math:`\mu`), not this
   recursion.  No such scheme exists in this codebase.  Conversely,
   neither MoC nor CP is likely to acquire condition 2 or condition 1
   respectively without ceasing to be MoC or CP.

.. note:: **Two labels, one recursion.**  :eq:`alpha-dome-recursion`
   here and :eq:`alpha-recursion` on
   :doc:`/theory/methods/sn/curvilinear_one_group` state the *same*
   recurrence.  This page carries it in the **geometry-primitive**
   register — what the chart owes a consumer, seeded at
   :math:`\alpha_{1/2} = 0` — and the methods page in the
   **discretisation** register, where it enters the cell update with
   :eq:`alpha-cylindrical` as the per-level arm.  Only the methods-page
   label is a ``verifies`` target.  The duplication is a genuine
   single-source-of-truth smell; collapsing it is **recommended and
   deliberately not done here**, because it would move a generated
   V&V-matrix row and re-point ``verifies`` markers that a
   documentation pass may not edit.

Migration roadmap
-----------------

This primitive is the foundation for several follow-on issues in the
SN reshape campaign (``.claude/plans/sn_reshape.md``):

* **Issue 10 (Wave D Round 1.1) — DONE**: :class:`SNProblem` consumes
  :class:`ReducedStreamingOperator` via the dispatch ladder above.
  The connection-coefficient math no longer lives in :class:`SNProblem`.
* **SN operator algebra (Depth B, 2026-05)** —
  :class:`~orpheus.sn.operators.streaming.StreamingOperator` /
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` consume the
  primitive (as ``SNProblem.reduced``) through the loss-representation walk:
  :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` reads the
  connection coefficients off ``self.mesh.coord`` inside the walk.
  (Depth B consumed them through the per-geometry
  ``transport_operator_matvec_*`` matvecs; that family and its unified
  successor were deleted in the typed-field (#197) and walk-unification
  (#280 campaigns) refactors — the primitive itself is unchanged.)
* **MoC and CP campaigns (post-Wave-1)** — ⛔ **retracted 2026-08-27.**
  This item read *"reuse the same primitive with their own consumption
  patterns (track-segment chord march for MoC; ray-traced chord-length
  integrals for CP)"*.  It was never a description of the tree, and it
  is not a migration still owed: neither method forms an angular
  redistribution term at all, for two independent structural reasons
  worked in :ref:`who-needs-a-connection-coefficient` above.  The
  roadmap item is closed as **not applicable**, not as pending.  What
  MoC and CP *do* share with S\ :sub:`N` is the L2 transport layer —
  fields, sources, cross-section data, the scattering kernel and the
  eigenvalue driver (:mod:`orpheus.transport`) — which is a different
  and genuinely satisfied sharing claim.


.. _structured-geometry-history:

Development history
===================

Reverse-chronological (latest first) changelog of this page's subject,
the geometry value. Entries marked *(in development)* live on an
unmerged feature branch and have no landed merge-to-``main`` hash yet;
trust ``git`` over this table for merge status.

.. list-table::
   :header-rows: 1
   :widths: 10 54 10 26

   * - When
     - Milestone
     - Issue
     - Where
   * - 2026-09-29
     - **The mesh refines the geometry, and a Mesher builds it.** The
       measure got one definition, ``CoordSystem.measure``, asked
       through the geometry. The interval rules (``CellsByCount``,
       ``CellsByMaxWidth``, ``Refined`` as ``k * rule`` for :math:`k` a
       power of two, ``CellEdges``) and the spacing rules
       (``EqualWidth``, ``EqualVolume``, one body in the measure
       coordinate) replaced ``RegionMesh``; one rule or one per interval,
       so mixed spacing became spellable. ``Mesh1D`` became
       ``(coord, edges, volumes, mat_ids, face_laws)`` with no geometry,
       stored volumes checked against the measure within :math:`2p + 5`
       ulp, and one declared law per boundary face; ``from_geometry``,
       ``precomputed_volumes``, ``bc_left`` / ``bc_right``, the subdivision
       helper and ``None`` as a ``Mesh1D`` face law retired. The
       ``Mesher`` session (``partition``, ``refine``, ``mesh``) became the
       only construction path. #495 was fixed at its root (on a slab the
       two spacings are one body), and ERR-020's fix became the stored
       share :math:`m/n`. The named-face geometry constructors
       (``slab``, ``cylinder``, ``sphere``, ``uniform_boundary``,
       ``from_homogeneous``) landed with it. Rulings: the user,
       2026-09-29, P1 step 3 of ``.claude/plans/reference_cache.md``;
       what was refuted on the way: :ref:`structured-geometry-mesh-refuted`.
     - #405, #495, #539
     - *(in development)* branch ``refactor/reference-specification``
   * - 2026-09-29
     - **Each reference generator serves the bodies its solvers
       solve.** The ``homogeneous_body`` reading, which refused every
       multi-material or hollow geometry for all four generators,
       retired in favour of
       :func:`~orpheus.derivations.common.reference_body.reference_body`,
       a total classification into four shapes, and one refusal door,
       :func:`~orpheus.derivations.common.reference_body.refuse_unserved`.
       ``Billiard`` gained reachable multi-region sphere, hollow
       sphere and annulus arms, and a new multi-region cylinder arm;
       ``MomentSpace`` gained the symmetric reflected slab, and the
       cross-method adapter for Sood's problem 4 went through it
       instead of calling the solver. The boundary laws joined what a
       generator serves: one reader,
       :func:`~orpheus.derivations.common.reference_body.specular_albedo`,
       turns a law into a specular albedo; ``Billiard``'s ``alpha``
       parameter and ``with_alpha`` retired, and ``Spectrum``'s
       ``_extract_R_refl`` with them. The multi-region sphere's
       fixed-source arm lost a defect on the way (ERR-091). The step
       closed #190 and #421.
       Rulings: the user, 2026-09-29, P1 step 2b of
       ``.claude/plans/reference_cache.md``.
     - #405, #190, #421, #536
     - ``1aab17a6``
   * - 2026-09-29
     - **The geometry value is (coordinate system, breakpoints,
       material per interval, law per boundary point).** The string
       kind tag (``"SLB"`` / ``"CYL"`` / ``"SPH"``), the ``Region``
       class, and the ``regions``, ``bcs`` and ``n_endpoints`` fields
       retired, with the tag maps ``_GEOMETRY_TO_COORD`` and
       ``_GEOMETRY_TO_N_ENDPOINTS``. The law count became derived from
       the topological boundary (so hollow cylinders and spheres became
       declarable, and the S\ :sub:`N` / diffusion inner-law guard of
       #511 landed with them); breakpoints became stored rather than
       re-added from thicknesses; ``Mesh1D.from_geometry`` lost its
       ``origin=`` argument (the first breakpoint is the origin); the
       four one-material reference generators gained the shared
       ``homogeneous_body`` reading, which refused the multi-material
       geometry they had read as its first region's material. CP,
       MoC and MC gained refusals of the laws they drop (#513, #514),
       and the admitted S\ :sub:`N` inner law gained its void-cavity
       witness. Rulings: the
       user, 2026-09-29, P1 step 2 of ``.claude/plans/reference_cache.md``;
       gates S2.1 to S2.5 of ``.claude/plans/reference_p1_spec.md``.
     - #405, #511, #513, #514
     - ``b91f1591``
   * - 2026-05-04
     - **Phase F: the geometry layer separates from the registry and
       mesh layers.** ``StructuredGeometry`` was then a kind tag, a
       tuple of ``Region(mat_id, outer_thickness_cm)`` and a tuple of
       endpoint ``bcs`` whose length the tag fixed (``SLB`` 2, ``CYL``
       and ``SPH`` 1, the centreline described as "implicit
       reflective"). The endpoint count was justified by the billiard's
       orbit-space rank, and hollow bodies were anticipated as new tags
       (``HSPH``, ``ANN``). Both readings were superseded on
       2026-09-29: the count follows from :math:`r_0`, not from the
       coordinate system alone, so no new tag is needed, and the centre
       of a solid body is an interior point, not a reflecting surface.
     - —
     - plan ``.claude/plans/dazzling-cuddling-boot.md``
