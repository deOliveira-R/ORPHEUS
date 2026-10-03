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
* **Equality and hash are content, through one encoder.** The geometry,
  its boundary laws, the ``BC`` tags, the mesh, its ``FaceLaws``, the
  mixtures and the ``Materials`` declaration all derive ``==`` and
  ``hash`` from one content digest,
  :func:`~orpheus.numerics.content.content_digest` (blake2b-256 of a
  type-tagged, length-prefixed encoding), so two values that are the
  same physics are equal and hash alike in every process and on every
  platform, which is what a persistent cache key needs. A real scalar is
  its value (``1 == 1.0``, ``-0.0`` is ``+0.0``), NaN is refused when a
  value is constructed, and an object's digest covers its class's
  schema, so an entry written under an older schema misses
  (:ref:`structured-geometry-content-identity`).
* **A source and a detector are stated before any mesh.**
  :class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant` is one
  real value per (region, group), a function on the angle-integrated
  space whose regions are the geometry's interval indices;
  :class:`~orpheus.numerics.mesh_free_function.Symbolic` is one SymPy
  expression :math:`q_g(r,\mu,\varphi)` per group, a function on phase
  space stored as ``srepr`` text with the SymPy version as content.
  Neither carries a role or a density: a table enters phase space through
  the angular section as a source and through the retraction's adjoint as
  a detector, and the measure's mass between them is never typed. Each
  coordinate system declares the angular chart :math:`(\mu,\varphi)` is
  read in, and the sphere declares no azimuth reference
  (:ref:`structured-geometry-mesh-free-functions`).
* **What is asked is a value with no physics in it.**
  :class:`~orpheus.numerics.question.Eigen` ``(parameter, point, mode)``
  asks where, along one direction of the system's parameter space, the
  system is singular; :class:`~orpheus.numerics.question.FixedSource`
  ``(source, point)`` asks for the flux a source drives and
  :class:`~orpheus.numerics.question.Response` ``(detector, point)`` for a
  detector's importance. The parameter and the point's keys are opaque
  keys a specification (later a system) resolves; the point is a frozen
  mapping of offsets from the physical value, empty by default; the mode
  is ``Fundamental()`` or ``Nearest(tau)``. No value carries an adjoint
  flag: the question's type is its role, and the eigen adjoint belongs to
  the answer. Nothing behind the values (the pencil a parameter derives,
  the mode law, the system) exists yet: #529
  (:ref:`structured-geometry-question-values`).
* **Hollow cylinders and spheres are declarable, and every method
  refuses a declared law it would drop.** S\ :sub:`N` and diffusion
  admit only a reflective inner law on a hollow body (#511), which is
  verified to be the void cavity it models; CP admits a slab only with
  a reflective left law, the mirror its slab kernel computes at the left
  face whatever is declared there, and no inner law (#513); MoC admits only a
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
  plane, built by ``from_thicknesses``). Each carries its model's laws
  and takes no law argument: white on the Wigner–Seitz cell's outer
  surface, reflective on both faces of the half-cell (candidate 9 of
  :ref:`structured-geometry-mesh-refuted`).


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
``boundary_points``, the mesh's ``boundary_points`` and the face
inventory of both meshes (:ref:`structured-geometry-face-laws`) read it; the
derived property ``is_hollow`` is true exactly for a cylinder or a
sphere with :math:`r_0 > 0`, and never for a slab, which has no centre.
One check, ``parse_boundary_laws``, requires one law per boundary point
of the geometry's ``boundaries`` (a mesh's ``face_laws`` are checked
against the face inventory, with the same reasons in its messages). The
refusals are keyed to the three ways a law count can be wrong, each
naming the reason in its message:

* **a law at the centre** — two laws on a solid cylinder or sphere:
  *"the centre r = 0 of a solid … body is an interior point and
  carries no law"*;
* **a hollow body missing its inner law** — one law with
  :math:`r_0 > 0`: *"a hollow … body … has an inner surface, which
  needs its own law"*;
* **a slab with one law**: *"a slab has two boundary points (left,
  right)"*;

and, before the count is read, ``None`` in any position, refused by
the element parser ``parse_boundary_law`` that ``parse_boundary_laws``
calls on each law: *"StructuredGeometry.boundaries[0] is None, and None
is not a boundary law: declare the law the boundary carries."* ``None``
is refused because it means nothing there: a geometry declares the
problem, and a default is a method's, not the problem's. Every other
boundary declaration in the tree (the face laws of
:class:`~orpheus.mesh.structured.Mesh1D` and
:class:`~orpheus.mesh.structured.Mesh2D`, and the endpoint laws of the
S\ :sub:`N` axes) passes through the same parser and is refused the
same way (:ref:`structured-geometry-no-default-law`).

The geometry's laws are indexed by boundary point, (inner, outer),
never by a coordinate-specific name (``left``, ``centreline``,
``outer``), and the Mesher lifts them onto the mesh's named faces
through the face inventory, pairing the tuple with the inventory in
order: two boundary points give the faces ``xmin`` and ``xmax``; one
point gives ``xmax`` alone, since the centre is an interior point. The
geometry keeps its positional tuple, so one declaration has two
spellings across the lift, a tuple on the geometry and a
:class:`~orpheus.mesh.face_laws.FaceLaws` on the mesh. The call site names each law through the named-face
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
two laws of a hollow body onto ``face_laws["xmin"]`` and
``face_laws["xmax"]``, the positions onto ``boundary_points``
:math:`= (r_0, r_R)`, and the first edge onto :math:`r_0`.


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
   assert dict(mesh.face_laws) == {"xmin": BC.vacuum, "xmax": BC.reflective}
   assert mesh.boundary_points == (0.0, 3.0)


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
* ``face_laws`` are any mapping from face name to law over exactly the
  mesh's face inventory: ``xmin`` and ``xmax`` on a slab or a hollow
  cylinder or sphere, ``xmax`` alone on a solid one. Each law is a
  :class:`~orpheus.geometry.boundary.BC` tag or a typed
  ``BoundaryTraceLaw``; ``None`` is refused with the message the
  geometry gives (*"Mesh1D.face_laws['xmin'] is None, and None is not a
  boundary law: declare the law the boundary carries."*), because every
  declaration goes through one element parser, ``parse_boundary_law``.
  The mesh stores a :class:`~orpheus.mesh.face_laws.FaceLaws`
  (:ref:`structured-geometry-face-laws`).

The derived quantities are ``widths``, ``centers``, ``areas``,
``boundary_points`` (the positions of the boundary faces, inner first,
in the order of ``face_laws``) and ``outer_law``, which reads
``face_laws["xmax"]`` (the law on :math:`r = r_R`, a slab's right face,
the one law collision probability, characteristics and Monte Carlo
read). Equality and hash are the content identity of the five fields
(:ref:`structured-geometry-content-identity`): two meshes with the same
coordinate system, edges, volumes, material ids and face laws are equal
and hash alike in every process, whatever objects built them, and the
derived ``widths``, ``centers`` and ``areas`` are not content
(``compare=False``), because they are functions of the five. The
mesh's digest is the content key of the discretisation, and it reads no
material data, only the material ids.

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
geometry stays one layer down. The laws are **per face**, keyed by the
face's name (:ref:`structured-geometry-face-laws`): in 1-D each
boundary point is one face, and the per-face form is the seed for 2-D,
where one side of the boundary may be several faces carrying different
laws. A specialised
``(geometry, edges, volumes)`` constructor was considered and not built:
its one consumer would be the Mesher's own lift, while a relabelling
(``with_distinct_cell_ids``), an adaptation or the axis adapter starts
from a mesh or from axes.


.. _structured-geometry-face-laws:

.. _structured-geometry-mesh2d:

The face laws: one inventory rule, one value, both meshes
---------------------------------------------------------

A mesh's boundary is a finite set of faces, each carrying one law. Both
meshes, :class:`~orpheus.mesh.structured.Mesh1D` and
:class:`~orpheus.mesh.structured.Mesh2D`, derive the set by one rule
and store the laws in one value, both in :mod:`orpheus.mesh.face_laws`.

**The inventory rule.**
``face_inventory(coord, first_axis_edges, dimension)`` names the
boundary faces of a mesh, axis by axis, inner face first:

* along the first axis (:math:`x`, or :math:`r`), the coordinate
  system's boundary points of :math:`[x_0, x_{N}]`,
  :meth:`~orpheus.geometry.coord.CoordSystem.boundary_points`, the rule
  the geometry uses: both ends on a slab and on a hollow cylinder or
  sphere, the outer surface alone on a solid one, whose centre or axis
  :math:`r = 0` is interior and carries no law;
* along every further axis (the :math:`y` of :math:`(x, y)`, the
  :math:`z` of :math:`(r, z)`), both ends;
* each face is named by the crosswalk
  :class:`~orpheus.mesh.axis.FaceLabel` (axis index and endpoint give
  ``face_name``), a solid radial axis's outer surface being ``xmax``.

.. list-table:: The face inventory, by mesh
   :header-rows: 1
   :widths: 34 30 36

   * - Mesh
     - Faces
     - Why
   * - ``Mesh1D``, slab
     - ``xmin``, ``xmax``
     - both ends bound the slab
   * - ``Mesh1D``, solid cylinder or sphere (:math:`r_0 = 0`)
     - ``xmax``
     - the centre is an interior point
   * - ``Mesh1D``, hollow cylinder or sphere (:math:`r_0 > 0`)
     - ``xmin``, ``xmax``
     - the inner surface :math:`r = r_0` bounds the body
   * - ``Mesh2D``, :math:`(x, y)`
     - ``xmin``, ``xmax``, ``ymin``, ``ymax``
     - both ends of both axes bound the rectangle
   * - ``Mesh2D``, solid :math:`(r, z)`
     - ``xmax``, ``ymin``, ``ymax``
     - the axis :math:`r = 0` is an interior line and carries no law
   * - ``Mesh2D``, hollow :math:`(r, z)`
     - ``xmin``, ``xmax``, ``ymin``, ``ymax``
     - the inner surface bounds the body

The names are the ones S\ :sub:`N`'s resolved boundary table
(``SNProblem.bc``) is keyed by, and the ones the axis adapter reads and
writes (``axes_from_legacy_mesh`` reads ``face_laws["xmin"]`` and so on
into the axes' endpoint slots; ``legacy_mesh_from_axes`` builds the
mapping from ``face_labels(axes)``), so a face has one name from its
declaration to its realised operator. A consumer that asks whether a
body has an inner face asks ``"xmin" in mesh.face_laws``.

**The value.** :class:`~orpheus.mesh.face_laws.FaceLaws` is a frozen,
ordered, picklable mapping from face name to law, a subclass of
:class:`~orpheus.numerics.content.FrozenMapping`. A mesh is given its
laws as any mapping keyed by face name, and builds the value with
``FaceLaws.over(inventory, laws, where, coord)``, which refuses, each
with a keyed message:

* a declaration whose key set is not the inventory exactly, a missing
  face and an extra face alike; the message names the inventory, and
  adds the reason when the difference is ``xmin``: the centre of a solid
  body carries no law, or a hollow body's inner surface needs its own;
* ``None``, or any object that is not a ``BC`` tag or a
  ``BoundaryTraceLaw``, on a face, through the one element parser
  (*"Mesh2D.face_laws['xmin'] is None, and None is not a boundary law:
  declare the law the boundary carries."*); a string such as
  ``"vacuum"`` is not a tag and is refused too;
* a declaration that is not a mapping, the retired positional tuple
  among them.

The value iterates the face names in inventory order, whatever the
order of the caller's mapping, so ``tuple(mesh.face_laws)`` is the
inventory; assigning to it raises. Its equality and hash are content
identity (:ref:`structured-geometry-content-identity`): the same faces
carrying equal laws, in any order, are one value. A ``FaceLaws`` is
equal only to another ``FaceLaws``, never to a plain ``dict`` of its
items; to compare its items with a dictionary, compare
``dict(mesh.face_laws)``, as the examples below do.
It is picklable because the reference-solution cache of #405
pickles its entries, and a mesh is part of what it stores: `[M]` the
read-only view the first 2-D spelling stored, a ``MappingProxyType``,
raises ``TypeError: cannot pickle 'mappingproxy' object``.

``Mesh1D(coord, edges, volumes, mat_ids, face_laws)`` passes
``dimension=1`` to the rule, and
``Mesh2D(edges_x, edges_y, mat_map, *, face_laws, coord=CARTESIAN)``,
whose ``face_laws`` is keyword-only with no default, passes
``dimension=2``. ``Mesh1D.boundary_points`` gives the positions of its
faces (the vocabulary of ``CoordSystem``), and ``Mesh1D.outer_law``
reads ``face_laws["xmax"]``.

.. code-block:: python

   import numpy as np
   from orpheus.geometry import BC, StructuredGeometry
   from orpheus.geometry.coord import CoordSystem
   from orpheus.mesh import CellsByCount, Mesh2D, Mesher

   rod1d = Mesher(StructuredGeometry.cylinder((0.0, 1.0), (0,), outer=BC.white)).partition(
       CellsByCount.uniform_volume(4)).mesh
   assert tuple(rod1d.face_laws) == ("xmax",)
   assert rod1d.outer_law == BC.white

   box = Mesh2D(
       [0.0, 1.0, 2.0], [0.0, 1.0], np.zeros((2, 1), dtype=int),
       face_laws={"ymax": BC.vacuum, "xmin": BC.reflective,
                  "xmax": BC.vacuum, "ymin": BC.reflective},
   )
   assert tuple(box.face_laws) == ("xmin", "xmax", "ymin", "ymax")
   assert dict(box.face_laws) == {"xmin": BC.reflective, "xmax": BC.vacuum,
                                  "ymin": BC.reflective, "ymax": BC.vacuum}

   # A solid (r, z) mesh: the axis r = 0 carries no law.
   rod = Mesh2D(
       [0.0, 0.5, 1.0], [0.0, 2.0], np.zeros((2, 1), dtype=int),
       face_laws={"xmax": BC.vacuum, "ymin": BC.reflective, "ymax": BC.reflective},
       coord=CoordSystem.CYLINDRICAL,
   )
   assert tuple(rod.face_laws) == ("xmax", "ymin", "ymax")

   import pickle
   assert pickle.loads(pickle.dumps(rod.face_laws)) == rod.face_laws

A ``Mesh2D`` has no geometry value to be built from (there is no 2-D
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`), so
it is constructed directly.

The one 2-D factory, :func:`~orpheus.mesh.factories.pwr_pin_2d`, takes
a required keyword ``law`` for its four faces, because a factory does
not choose a physical boundary on the caller's behalf: a lattice cell
is reflective or periodic, an isolated cell vacuum.


.. _structured-geometry-no-default-law:

No boundary law is undeclared
-----------------------------

Every boundary declaration in the tree passes through one element
parser, ``parse_boundary_law`` in
:mod:`orpheus.geometry.structured_geometry`, which accepts a
:class:`~orpheus.geometry.boundary.BC` tag or a ``BoundaryTraceLaw`` and
refuses ``None`` and every other object. Its five callers are the five
places a law is declared:

.. list-table:: Where a boundary law is declared
   :header-rows: 1
   :widths: 34 66

   * - Declaration
     - Laws
   * - ``StructuredGeometry.boundaries``
     - one per boundary point, through ``parse_boundary_laws``
   * - ``Mesh1D.face_laws``
     - one per boundary face, keyed by face name, through
       ``FaceLaws.over``
   * - ``Mesh2D.face_laws``
     - one per boundary face, keyed by face name, through
       ``FaceLaws.over``
   * - ``AxisMesh.bc_low``, ``AxisMesh.bc_high``
     - both required, parsed at construction
   * - ``RadialAxisMesh.bc_outer``
     - required, parsed at construction

The consumers take the declaration verbatim. The shared resolution,
:func:`~orpheus.transport.method.resolve_boundary_conditions`, reads the
law each axis declares for each face and parses the tag into its typed
law, with no default to fall back on; the S\ :sub:`N` entry points
(``solve_sn``, ``solve_sn_adjoint``, ``solve_sn_fixed_source``,
``solve_sn_adjoint_fixed_source``, ``solve_sn_multiplying_source``)
have no parameter that supplies a law, and nothing fills a face. The law
on a boundary is read from the problem's declaration alone.

The one law the tree supplies is not a default: the inner surface of a
hollow radial axis. :class:`~orpheus.mesh.axis.RadialAxisMesh` has no
slot for it, so S\ :sub:`N` realises the cavity as reflective, and the
geometry and the mesh refuse every other inner law at the axis adapter
(:ref:`structured-geometry-hollow-inner-law`, #511). A declared
``reflective`` inner law is therefore the only one that reaches it.

**Why there is no default.** A ``None`` that a consumer read as a law
was a behavioural default, and a behavioural default hides a
precondition no type states (lesson L19 of
:doc:`/development/evidence/lessons`: type the precondition or force an
explicit choice). Here the precondition was *which problem the caller
meant*, and the answer depended on the door the mesh went through. The
fixed-source entries filled an all-undeclared ``Mesh2D`` or axis tuple
with vacuum (their ``boundary_condition="vacuum"`` parameter); the
eigenvalue entry and a directly constructed ``SNProblem`` resolved the
same undeclared faces as reflective (the shared resolver's default);
and the fill was all-or-nothing, so a partly declared mesh kept
reflective faces even under a fixed-source entry. One mesh was two
problems, and nothing in the mesh said which. The operator algebra met
the same pattern in an operator's optional spaces, where ``None``
silently meant "Euclidean" (:ref:`bound-operator`), and removed it the
same way: the value must say what it is.

`[M]` The exposure before the retirement (the post-step-3b runtime
capture, ``-m "not slow"``, serial; the census
``scratch/reference_architecture/p1step3c/blast_set.md``): 272 tests in
36 files resolved an undeclared face through ``SNProblem``, as
face-resolutions ``xmin`` 293, ``xmax`` 295, ``ymin`` 282, ``ymax``
282, ``zmin`` 33 and ``zmax`` 33; 6 tests received an entry's fill; the
diffusion hub resolved none. Statically, 147 ``Mesh2D`` constructions
under ``tests/`` (68 of them missing a law) and 88 axis constructions
(56 missing a law, a floor: a law forwarded through a test helper's
default is invisible to the static count) had to state their laws.

Deciding what an undeclared face means is not the mesh's job, and not a
method's either: a problem whose boundary is unstated is not yet a
problem. The retirement therefore has no successor default anywhere.


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

The design went through the candidates below before the one described
above (candidates 1 to 4 for the 1-D mesh, 5 to 9 for the declared laws);
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
5. **Four required fields on** ``Mesh2D``, ``bc_xmin``, ``bc_xmax``,
   ``bc_ymin`` and ``bc_ymax`` with their ``None`` defaults removed.
   `[REFUTED 2026-09-29]` by the user's ruling for step 3c: a solid
   :math:`(r, z)` mesh would then have to declare a law on its axis
   :math:`r = 0`, which is not a boundary, or keep ``None`` there
   meaning "no face here", a second meaning of ``None`` beside the
   retired "use the default". The per-face mapping has neither problem:
   the axis is simply not a key. The fact that survives: in 2-D, as in
   1-D, the set of faces is derived from the coordinate system, so the
   declaration is keyed by the derived faces, not by fixed slots.
6. **Keeping** ``None`` **on** ``Mesh2D`` **and on the axes**, with only the
   entries' ``boundary_condition=`` fill retired and the resolver's
   reflective default kept as the one convention. This was the scope
   question of step 3c. `[REFUTED 2026-09-29]` by the user's ruling
   that ``None`` retires as a boundary declaration everywhere: under
   this candidate the method, not the declaration, decides what an
   unstated face is, so a mesh cannot be read as a problem on its own
   (:ref:`structured-geometry-no-default-law`).
7. **A default law on** ``pwr_pin_2d`` (the factory used to build its
   mesh with every face ``None``, which the fixed-source entries read as
   vacuum and the eigenvalue entry as reflective). Refused in the design
   of step 3c: the factory has no knowledge of whether its cell is a
   lattice member or isolated, so its ``law`` keyword is required.
8. **A local rename of the two face-law spellings, with their
   unification filed as an issue.** Before the review of step 3c,
   ``Mesh1D.face_laws`` was a positional tuple (inner face first) and
   ``Mesh2D.face_laws`` a read-only mapping keyed by face name, and
   ``boundary_faces`` meant positions on one and names on the other
   (`[M]` the elegance review: ``zip(m2.boundary_faces, m2.face_laws)``
   paired names with names, in silence). `[REFUTED 2026-09-30]` by the
   user's ruling: a rename would leave one concept in two spellings,
   the stop signal of Cardinal Rule 2, so both meshes store one
   ``FaceLaws`` derived by one inventory rule. The 2-D mapping was also
   unpicklable, which the step-5 cache could not have stored.
9. **A law parameter on the named cells.**
   ``wigner_seitz_pin_cell`` and ``pwr_slab_half_cell`` took an
   overridable ``boundaries=`` (default white, and reflective on both
   faces). `[REFUTED 2026-09-30]` by the user's ruling that a named
   cell's law is part of its model: the Wigner–Seitz cell is the
   cylindricalised lattice cell with isotropic (white) re-entry, and
   the half-cell's two faces are symmetry planes of the lattice
   (reflective). Two alternatives were refused with it: requiring the
   law (the name already fixes it) and keeping the overridable default
   (a default a caller can inherit without seeing it). Another law is a
   different body, built with ``StructuredGeometry.cylinder`` or
   ``StructuredGeometry.slab``; `[M]` 2 test sites overrode the law
   before the ruling.

Two smaller decisions went the same way. Bit identity of the
equal-volume edges with the retired subdivision helper was dropped in
favour of the one measure (0 of 414 captured intervals moved). The fixed
8-ulp volume band became :math:`2p + 5` once the sphere's measured
worst case sat within 0.38 ulp of it.


What the mesh layer does not carry
----------------------------------

* The inner surface of a hollow radial axis has no law slot on
  :class:`~orpheus.mesh.axis.RadialAxisMesh`. The axis adapter gives the
  adapter mesh of such an axis the reflective cavity S\ :sub:`N`
  computes there (``SCOPE-BOUNDARY[guard]``, #511); the geometry and the
  mesh refuse any other inner law before it reaches the axis.
* A ``Mesh2D`` has no geometry value to be built from, so no Mesher
  builds it.
* A cache keyed on the mesh's digest. The digest exists
  (:ref:`structured-geometry-content-identity`); no code in the tree
  stores an entry under it, and the specification that composes it with
  the materials, the question and the source is #405's.
* Quality, preview, adaptation and external meshers are #539.


The gates
---------

``tests/gates/mesh/test_mesh1d.py`` carries the constructor's laws
(``TestConstructionLaws``, ``TestTheVolumeLaw`` with the band and the
wrong-coordinate control, ``TestTheValue`` for equality,
``TestTheRetirements`` for the retired fields and constructors, and
``TestTheDefaultsAreGone``: no S\ :sub:`N` entry takes a
``boundary_condition``, the fill helper is gone, and the resolver
neither completes an undeclared face nor constructs a law).
``tests/gates/mesh/test_mesh2d_face_laws.py`` carries the face laws
(``test_the_inventory`` against hand-written 2-D inventories,
``TestRefusals``, the storage in inventory order,
``TestTheOneElementParser``, a route gate that every declaration of the
five classes calls ``parse_boundary_law``; ``TestTheOneInventoryRule``,
that both meshes name the same first-axis faces and both constructors
ask ``face_inventory``; ``TestFaceLaws`` for the value's order,
immutability, content equality (equal to another ``FaceLaws`` of the
same items in any order, and not to a plain ``dict``) and pickle round
trips of the value and
of both meshes; and ``TestTheNamedCells``, that the two named cells take
no law keyword and carry their model's laws).
``tests/gates/mesh/test_axis_adapter_laws.py`` carries the axes
(``TestTheAxisLaws``: an omitted law and a non-law refused on both
axis classes, the default-filling verb gone, no ``None`` in a ``bc``
table; ``TestTheRoundTrip`` through the adapter in one and two
dimensions), and ``test_axes_declared_laws_are_taken_verbatim`` in
``tests/gates/sn/solve/test_d3_admission.py`` checks, face by face,
that the S\ :sub:`N` entry realises a mixed declaration on a 3-D axis
tuple as declared.
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
is :math:`r_0` and whose face law ``face_laws["xmin"]`` is the inner
law. The same routing puts a slab's left law in ``face_laws["xmin"]``.

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
     - the outer law; on a slab the left face is always a mirror (the
       kernel's image term), which a reflective left law declares
     - a slab whose left law is not ``reflective``; any inner law on a
       hollow cylinder or sphere
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
On a slab its kernel computes a mirror at the left face whatever is
declared there. The slab's reduced collision probabilities carry, beside
the direct path between cells :math:`i` and :math:`j`, a reflected path
through the plane :math:`x = r_0`, whose optical length is the sum of
the two cells' optical distances from that plane (the reflected path
:eq:`dc-slab` of the slab :math:`E_3` kernel,
:doc:`/theory/methods/collision_probability`; ``gap_c = bnd_pos[i] +
bnd_pos[j]`` in ``CPMesh._compute_slab_rcp``). So the slab CP
computes is half of a slab symmetric about its left face, and a declared
left law has no slot to enter. `[M]` 2026-09-25 (the specification's
probe ``cp_slab_left_law.py``): a white, a vacuum and an undeclared left
law gave the same :math:`k = 1.8749980808246423` to all sixteen digits.
`[M]` 2026-09-30 (Census A's probe ``probe_cp_left.py``, re-measured on
the fixed tree): a one-group fuel | moderator slab on :math:`[0, 2]`
(the fuel 0.5 thick, 40 equal cells in each material; S\ :sub:`N` with Gauss–Legendre
:math:`S_{32}`), white on the right face; "flipped" puts the moderator
on the left:

.. list-table:: Which law CP computes at a slab's left face
   :header-rows: 1
   :widths: 40 20 20

   * - Solve
     - fuel on the left
     - flipped
   * - CP, left law ``reflective`` (bit-identical to the ``white`` left
       law CP was given before 2026-09-30)
     - 1.215476
     - 1.212883
   * - S\ :sub:`N`, reflective left, white right
     - 1.215494
     - 1.212884
   * - S\ :sub:`N`, white on both faces
     - 1.212537
     - 1.212537

CP agrees with the S\ :sub:`N` mirror to :math:`1.5\times10^{-5}` and
:math:`8\times10^{-7}` relative, and differs from the white left face by
:math:`2.4\times10^{-3}` and :math:`2.9\times10^{-4}`; and CP's value
changes when the slab is flipped, which a white left face, symmetric
under the flip, would not do. A reflective left law is therefore the
only declaration that means what CP computes, and it is the one
admitted; any other left law, an equal right law included, is refused.
On a hollow cylinder or sphere what CP realises at the inner surface is
not established, so no inner law is admitted.

.. dropdown:: What was admitted first, and why it was wrong
   :color: muted

   From P1 step 2 (2026-09-29) until 2026-09-30 the guard admitted a
   left law EQUAL to the right one, on the reading that equal laws are
   what CP computes, and this page recorded the left face's law as
   open, verbatim: "a fuel | moderator slab with both faces white and
   its mirror image give :math:`k` differing by :math:`5.0\times10^{-5}`
   relative, and neither equals the mirrored double slab, so the left
   face is neither the right law nor a mirror" (probe ``cp_mirror.py``).
   The inference does
   not follow. Under the mirror reading the two orientations are two
   different bodies (a mirror, then fuel, moderator and a white face;
   a mirror, then moderator, fuel and a white face), so their
   eigenvalues must differ; and a doubled slab solved by CP is itself
   mirrored at its own left face, so it is a third body, not the
   reference the comparison needed. The table above compares CP with
   S\ :sub:`N` solved under each candidate left law instead, and CP
   matches the mirror. The equal-law guard was therefore admitting a declaration
   CP replaced: a white left law was solved as a mirror (Census A, row
   G1; the guard and this page corrected 2026-09-30).

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
left law is not the mirror refused, an equal white left law included,
a reflective left law built under a white and a vacuum right law, a
hollow cylinder and sphere with an inner law refused; MoC: a slab and a hollow
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
   * - ``BC.reflective``, :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - 1 (a mirror is a symmetry and has no amplitude)
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
reading of it. `[M]`, by calling the reader on each law:
``ReflectiveBoundary`` reads 1 on each of the axes ``x``, ``y`` and ``z``
and ``AlbedoBoundary(0.3, SpecularReturn())`` reads 0.3 (2026-10-01);
``BC("white")``, ``BC("periodic")`` and
``AlbedoBoundary(0.3, IsotropicReturn())`` are refused (2026-09-29). A
partially specular face is declared as
``AlbedoBoundary(α, SpecularReturn(axis))``; the mirror has no albedo to
read (:ref:`bc-deck-law-not-an-operand`).

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


.. _structured-geometry-content-identity:

Content identity: one encoder for every value a key covers
==========================================================

Every value a reference question is posed from (the mixtures and the
``Materials`` declaration, the geometry with its boundary laws, the
mesh with its face laws) derives its equality and its hash from one
content digest. The digest is computed by one encoder,
:mod:`orpheus.numerics.content`, whose module docstring is the design
record; this section states the problem it solves, the canonical form
of each kind of value with its reason, the schema tag, the two
refusals, what the encoder replaced and the gates that hold it. It
landed as step 5 of the first phase of the reference-solution campaign
(#405), on 2026-10-02.

The axis of a function space uses the same encoder: an
:class:`~orpheus.numerics.axis.Axis` is a content-identity value whose
parts exclude its generator, and every derived space name is a digest
of its content (:ref:`spaces-identity-bridge`,
:ref:`spaces-generator-identity-exclusion`).

.. _structured-geometry-content-identity-problem:

Why one content identity, and why it is the digest
--------------------------------------------------

The reference-solution cache of #405 stores the answer to a question
under a key computed from the question's content: the materials, the
geometry with its laws, and the discretisation. The key must have two
properties at once. Two values that are the same physics must give the
same key in every process and on every platform, or a stored answer is
never found again. Two values that differ in anything a consumer reads
must give different keys, or the cache serves one problem's answer to
another.

Python's ``hash`` has neither property across processes. The hash of a
``str`` or ``bytes`` is salted per interpreter (``PYTHONHASHSEED``), so
a tuple key built from them changes from one run to the next. `[M]` at
``1dc31163``, the hash of the 2-group mixture A of
``orpheus.derivations.common.xs_library`` was 2365733758199073423 under
``PYTHONHASHSEED=1`` and -2362898479209405550 under ``PYTHONHASHSEED=2``
(the test-architect's probe, ``.claude/plans/reference_p1_spec.md``
§1.5, correction 7).

A per-type equality cannot provide the key either. Each type that
defines its own ``__eq__`` defines its own "same value", and a key folded
from several such types inherits every disagreement among them: one type
separates ``-0.0`` from ``+0.0`` and another does not, one admits NaN,
one compares by identity, one cannot be hashed. So there is one
definition of "the same content", the encoder, and every type's ``==``
and ``hash`` is derived from its digest (instrument doctrine X4, one
definition per quantity): equality and the cache key cannot drift
apart, because they are one computation.

.. dropdown:: What each type's equality was before the encoder (`[M]` at ``1dc31163``)
   :color: muted

   The explorer's census measured each type on the tree the step
   started from (``scratch/reference_architecture/p1step5/census.md``
   §2 and §4):

   * a ``Mixture`` compared its arrays by their bytes, so ``-0.0`` and
     ``+0.0`` gave two unequal mixtures, and so did one sparse matrix with
     and without an explicitly stored zero; a NaN cross section was
     admitted;
   * ``Materials`` compared by identity (``eq=False``: two declarations of
     the same mixtures were two values) and could not be pickled, because
     it held a read-only ``mappingproxy``, so no ``MaterialMesh`` could be
     pickled either;
   * a ``BC`` tag could not be hashed (its ``params`` was a ``dict``),
     that dict stayed mutable after construction (the shared constants
     ``BC.vacuum`` and the others included), its values were not checked
     (a ``str`` was admitted), and a NaN parameter made two equal-looking
     tags unequal;
   * a geometry holding a ``BC`` tag could not be hashed;
   * a ``FaceLaws`` compared equal to a plain ``dict`` of its items, and
     could not be hashed (``Mapping`` defines ``__eq__``, which removes the
     inherited hash);
   * a ``Mesh1D`` compared bitwise through a hand-written ``__eq__`` and
     could not be hashed;
   * a ``Mesh2D`` raised ``ValueError`` on ``==`` between any two distinct
     meshes (the generated ``__eq__`` compared arrays), aliased the
     caller's writeable arrays, and stored ``-0.0``;
   * ``VacuumInflow`` and ``ReflectiveBoundary`` hand-wrote ``__eq__`` and
     ``__hash__``.

.. _structured-geometry-content-identity-forms:

The canonical forms, each with its reason
-----------------------------------------

:func:`~orpheus.numerics.content.encode` turns a value into bytes,
recursively. Every chunk carries a one-byte type tag and an 8-byte
length, so no two different values share a byte stream: the encoding is
injective on the values it admits, and ``("a", "b")`` and ``("ab",)``
cannot collide. :func:`~orpheus.numerics.content.content_digest` is the
blake2b digest of those bytes, 32 bytes (256 bits), wide enough that a
collision between two specifications is not a failure mode the cache
considers. The canonical form of each kind of value follows Python's
``==`` wherever ``==`` is defined on it, because that is the user's
ruling of 2026-10-02: the digest follows ``==``.

.. list-table:: The canonical forms
   :header-rows: 1
   :widths: 20 40 40

   * - Value
     - Canonical form
     - Reason
   * - real scalar (``bool``, ``int``, ``float``, the numpy real
       scalars)
     - one IEEE-754 double, little-endian; ``-0.0`` becomes ``+0.0``;
       NaN is refused (``ValueError``, naming the path to it); an
       integer of magnitude above :math:`2^{53}` is refused
     - ``True == 1 == 1.0`` in Python, so ``BC("albedo", {"albedo": 1})``
       and ``{"albedo": 1.0}`` are one value; ``-0.0 == 0.0``; NaN is
       not equal to itself, so a value holding it has no equality to
       encode; an integer beyond :math:`2^{53}` has no exact double, and
       ``==`` would separate it from the double it rounds to
   * - boolean, integer or real array
     - its shape, then its entries as canonical doubles (``-0.0``
       becomes ``+0.0``, NaN and wide integers refused, the width
       checked on the integers themselves so an unsigned 64-bit entry
       cannot wrap); a complex array, and any other dtype, is refused
     - ``np.array_equal`` already identifies an integer material map
       with its float twin: the dtype is storage, not content. Rounding
       a ``float64`` value to ``float32`` is a change of value, and
       moves the digest
   * - sparse matrix
     - the matrix in compressed sparse row form, with duplicate entries
       summed, explicitly stored zeros eliminated and indices sorted,
       then its shape, row pointers, column indices and values
     - a sparse matrix is the matrix, not its storage: an explicitly
       stored zero equals its absence, and ``int64`` indices equal
       ``int32`` ones. Its stored values are checked for NaN and width
       before the cast to doubles; a complex matrix is refused
   * - ``str``, ``bytes``, ``None``
     - their bytes (UTF-8 for a string), each kind with its own tag
     - ``"1"`` is not ``1`` and ``b"a"`` is not ``"a"``
   * - ``tuple`` and ``list``
     - their elements in order, under two different tags
     - ``(1,) != [1]`` in Python
   * - mapping
     - its (key, value) pairs, ordered by the encoded key
     - the order of insertion is not content, so two declarations in
       different orders are one value; a ``str`` key and an ``int`` key
       stay distinct
   * - ``set``, ``frozenset``
     - the encoded elements, sorted
     - iteration order is not content
   * - ``Enum`` member
     - its class's ``module.qualname`` and its member name (checked
       before the scalars, since an ``IntEnum`` member is an ``int``)
     - two enumerations with a member of one name are different values;
       ``CoordSystem.SPHERICAL`` is not the string ``"spherical"``
   * - a :class:`~orpheus.numerics.content.ContentIdentity` value, or
       a frozen dataclass whose own equality is by value (``eq=True``)
     - the schema tag (next section), then its parts in schema order:
       the dataclass fields with ``compare=True`` unless the class says
       otherwise
     - a value is its class's schema and its parts; a field declared
       ``field(compare=False)`` is not content, the one spelling of
       that (an axis's generator, a mesh's derived widths)
   * - a mutable part inside an object: a ``list``, a ``dict``, a
       ``set``, a writeable array, a sparse matrix with writeable arrays
     - refused with
       :class:`~orpheus.numerics.content.ContentlessError`
     - the part could change after its owner was keyed, under a cached
       digest; a ``tuple``, a ``frozenset``, a
       :class:`~orpheus.numerics.content.FrozenMapping` and a read-only
       array are admitted. A value handed to ``content_digest`` directly
       is not a part and is not held to this
   * - anything else
     - refused with
       :class:`~orpheus.numerics.content.ContentlessError`, naming the
       path from the root value
     - a function, a plain object, a dataclass compared by identity
       (``eq=False``), a mutable dataclass or an object array has no
       content a persistent key can carry

**One limitation, by design.** The schema tag names a class by its module
and qualified name, so two classes defined under one name (inside one
function, on two calls) share a tag. A class defined inside a function is
therefore not a persistent type; nothing in the package defines one.

.. code-block:: python

   import pickle

   import numpy as np

   from orpheus.geometry import BC, CoordSystem, StructuredGeometry
   from orpheus.geometry.boundary import AlbedoBoundary, VacuumInflow
   from orpheus.numerics.content import ContentlessError, content_digest, encode

   # A real scalar is its value: the digest follows ==.
   assert BC("albedo", {"albedo": 1}) == BC("albedo", {"albedo": 1.0})
   assert hash(BC("albedo", {"albedo": 1})) == hash(BC("albedo", {"albedo": 1.0}))
   assert encode(-0.0) == encode(0.0) and encode(True) == encode(1)

   # Containers keep the distinctions == keeps; a mapping's order is not content.
   assert encode((1,)) != encode([1])
   assert encode({"a": 1, "b": 2}) == encode({"b": 2, "a": 1})
   assert encode({"1": 0}) != encode({1: 0})

   # An array is its values on its shape; the dtype is storage.
   assert encode(np.array([1, 2])) == encode(np.array([1.0, 2.0]))
   assert encode(np.zeros((2, 1))) != encode(np.zeros((1, 2)))

   # NaN is refused where the value is made, naming the parameter.
   try:
       AlbedoBoundary(float("nan"))
   except ValueError as err:
       assert "AlbedoBoundary.albedo" in str(err)
   else:
       raise AssertionError("a NaN albedo was admitted")

   # A tag and the typed law it resolves to are two declarations.
   assert BC.vacuum != VacuumInflow()

   # Two geometries built apart are one value, in every process.
   def sphere():
       return StructuredGeometry(
           coord=CoordSystem.SPHERICAL, breakpoints=(0.25, 0.5, 1.0, 2.0),
           mat_ids=(0, 1, 0), boundaries=(AlbedoBoundary(0.5), BC.vacuum),
       )

   a, b = sphere(), sphere()
   assert a is not b and a == b and hash(a) == hash(b)
   assert len(content_digest(a)) == 32
   assert pickle.loads(pickle.dumps(a)) == a

   # Anything else has no content.
   try:
       encode(lambda x: x)
   except ContentlessError as err:
       assert "has no content identity" in str(err)
   else:
       raise AssertionError("a function was encoded")

.. _structured-geometry-content-identity-schema:

The schema tag: a key covers the class's schema
-----------------------------------------------

An object encodes as its schema tag followed by its parts. The tag is
the text ``<module>.<qualname>|v<version>|<part names>``: the class's
module and qualified name, its ``__content_version__`` (a class
attribute, 1 unless the class raises it) and the names of its parts in
order. `[M]` the tag of a ``BC`` reads
``orpheus.geometry.boundary._tag.BC|v1|kind,params``. Adding, removing,
renaming or reordering a part, raising the version, and moving or
renaming the class each change every digest of that class.

This implements the user's ruling of 2026-10-01 on #405: the cache key
covers the schema of every persisted class, so an entry written under an
older schema MISSES. The lookup computes the key under the current schema
and finds nothing, and the stale object is never read. A ``__setstate__``
refusal on each class was declined, because it adds a guard per class and
per carve, while a schema-covering key handles every class and every
future carve once.

**Pickling goes through the constructor.** A default unpickle restores an
object's fields without running its ``__post_init__``: its arrays come
back writeable and its admission laws never run, so a pickle written
before a carve loads silently into the class as it is after the carve
(the reflective cleanup's qa review measured a
``ReflectiveBoundary('x', 0.7)`` pickled before the albedo was removed
loading as a mirror that ignores its 0.7, and a pickled
``0.7*R + 0.3*W`` loading as an admitted ``LawSum``;
``scratch/boundary_ontology/reflective_qa_review.md``, G1). So
``ContentIdentity.__reduce__`` pickles a dataclass value as its class and
its ``init`` fields and rebuilds it by calling the constructor: every law
re-runs on load, the arrays are read-only copies again, and a pickle
written under an older schema (a removed field) fails to load with a
``TypeError`` naming the field instead of loading as a value it is not.
``FrozenMapping`` pickles as its items, and ``PrescribedInflow`` has its
own ``__reduce__``. The key and the pickle therefore guard the schema
twice, on the two routes a stored value can take.

Three consequences follow, each a rule for whoever changes a class.

* **A change in what a part means, with its name unchanged, raises
  the class's content version.** Set ``__content_version__`` one higher:
  the names and the order are in the tag already, but a new unit or a
  new convention for an existing part is not, and only the version says
  it.
* **Moving or renaming a content class invalidates every cached entry
  that holds it.** The tag carries the module path (``_tag`` in the
  example above, a private module), so even a pure move such as P1
  step 1's changes the digests. The result is a cache miss and a
  recomputation, never a wrong answer.
* **The tag separates classes with no parts.** ``VacuumInflow`` and
  ``ZeroFluxBoundary`` have no fields, so their parts are the same empty
  sequence and only the class name distinguishes them; `[M]` the
  test-architect's battery arm that drops the tag reddens the
  ``LawSum`` row whose operand changes from one to the other (spec
  §1.5, correction 10). The tag is load-bearing for these types, not
  only for schema evolution.

.. _structured-geometry-content-identity-mixin:

Equality and hash are the digest
--------------------------------

:class:`~orpheus.numerics.content.ContentIdentity` is the mixin that
turns the digest into ``==`` and ``hash``. A class takes it in one of
two ways: as a frozen dataclass declared ``eq=False`` (so the generated
``__eq__`` and ``__hash__`` cannot shadow the mixin's), whose parts are
by default its fields with ``compare=True``; or by overriding
``content_parts()`` to return the named parts itself, which only
``FrozenMapping`` does. The mixin declares ``__slots__ = ()``, so it adds
no instance storage to a slotted class.

* ``a == b`` holds when ``a is b``, or when ``type(a) is type(b)`` and
  their digests agree. A subclass value is never equal to a value of
  its base class: an ``EnergyAxis`` and a generic ``Axis`` with the
  same fields are two values.
* ``hash(a)`` is the first 8 bytes of the digest read as a signed
  integer, so it is the same in every process.
* The digest is computed once per object and cached in a module-level
  table keyed by the object's ``id``, with an entry that is removed
  when the object dies. It is not stored on the instance: a digest in
  the instance's ``__dict__`` would be pickled with it and read back
  under a later schema, which is the stale key the schema tag exists to
  prevent. The cache is sound only because the values are frozen and
  their arrays are read-only copies.

`[M]` 22 concrete classes carry the mixin on this tree (23 with the
abstract ``BoundaryTraceLaw``), walked through ``__subclasses__`` at
runtime:

.. list-table:: The content classes, and their parts
   :header-rows: 1
   :widths: 30 70

   * - Class
     - Parts
   * - ``Mixture``
     - every cross-section field, ``chi`` and the energy grid ``eg``
       (``None`` when the mixture has no grid); each sparse Legendre
       block as a matrix
   * - ``Materials``
     - one part, ``mixtures``: a ``FrozenMapping`` from material id to
       ``Mixture``. The ids are coerced to ``int`` at admission (an
       ``np.int64`` id is the same id); the declared order is kept for
       iteration and is not content
   * - ``BC``
     - ``kind`` and ``params``, a ``FrozenMapping`` from parameter name
       to real number; ``bc.params`` is never equal to a plain ``dict``,
       so compare ``dict(bc.params)``
   * - the seven registered boundary laws (``VacuumInflow``,
       ``ReflectiveBoundary``, ``WhiteBoundary``, ``AlbedoBoundary``,
       ``PeriodicBoundary``, ``PrescribedInflow``,
       ``ZeroFluxBoundary``), and the law algebra's ``LawSum`` and
       ``LawScaled``
     - their dataclass fields; the parts of a law (``SpecularReturn``,
       ``IsotropicReturn``, ``NoSource``, ``ConstantInflowSource``) are
       frozen dataclasses, encoded by the same rule without carrying
       the mixin
   * - ``StructuredGeometry``
     - ``coord``, ``breakpoints``, ``mat_ids``, ``boundaries``
   * - ``Mesh1D``
     - ``coord``, ``edges``, ``volumes``, ``mat_ids``, ``face_laws``; the
       derived ``widths``, ``centers`` and ``areas`` are
       ``compare=False``
   * - ``Mesh2D``
     - ``edges_x``, ``edges_y``, ``mat_map``, ``face_laws``, ``coord``,
       stored as read-only copies with ``-0.0`` canonicalised
   * - ``FaceLaws``
     - a ``FrozenMapping`` subclass: one part, ``items``, the mapping
       from face name to law, so the order of the faces is not content
   * - ``CellEdges``
     - ``edges``
   * - ``Axis``, ``EnergyAxis``, ``HarmonicAxis``, ``LegendreAxis``
     - their dataclass fields: ``label``, ``shape``, ``weights``,
       ``kind``; ``EnergyAxis`` adds ``edges`` and ``LegendreAxis`` adds
       ``spent_axis``. The ``generator`` is declared
       ``field(compare=False)``: provenance, not content
       (:ref:`spaces-generator-identity-exclusion`)
   * - ``FrozenMapping``
     - one part, ``items``: its (key, value) pairs, order-free. It keeps
       the declared order for iteration, is frozen and picklable, and is
       equal only to a mapping of its own type, never to a plain
       ``dict``. It is the one spelling of a frozen mapping in the
       content types

**Two exceptions, both kept by ruling.** ``VacuumInflow`` and
``ReflectiveBoundary`` keep a string arm in ``__eq__``:
``VacuumInflow() == "vacuum"`` and ``ReflectiveBoundary("x") ==
"reflective"`` are ``True``, and every other comparison is the mixin's.
Retiring string equality is its own cleanup (the user, 2026-10-01).
The cost is stated here so nobody relies on the opposite: `[M]`
``hash(VacuumInflow()) != hash("vacuum")``, so the rule "equal values
hash alike" does not hold across the string arm, and a law and its
kind string must not be mixed as keys of one ``dict`` or members of one
``set``.

**A tag is not its typed law.** ``BC.vacuum`` and ``VacuumInflow()``
are different declarations and unequal (`[M]`), so two geometries that
mean the same face but spell it differently have different digests. The
failure this produces is a cache miss, never a wrong hit.

**Where the in-process keys read it.**
``MaterialMesh._contractibility_key`` folds each mixture's
``content_digest``, and ``MaterialMesh``'s boundary-law
key is now the law itself for a ``BC`` tag as for a typed law, since
both hash by content; only a law with no content keeps the key
``(qualified class name, id(law))``, which is honest in process because
a callable has no content to compare. These keys are hashed by Python
and live inside one process; the persistent key is the digest.

.. _structured-geometry-content-identity-refusals:

The two refusals: no content, and NaN
-------------------------------------

**A value with no content is refused, never guessed.** The encoder
raises :class:`~orpheus.numerics.content.ContentlessError` (a
``TypeError``) on a part it cannot encode, and the message names the
path from the root value, for example
``_WithGeneratorInKey.generator: Quadrature is a mutable dataclass,
whose content can change after it is keyed`` (`[M]`, the recipe of
:ref:`spaces-generator-identity-third-answer`). A content-identity value
holding such a part has no digest. It is then equal only to itself
(identity is the honest equality of a value with no content) and it is
unhashable, because an identity hash would let it into a persistent
key. The case in the tree is a ``PrescribedInflow`` whose source is a
plain object or a function, such as the manufactured inflow of
``tests/gates/sn/verification/analytical/test_mms_declared_inflow.py``:
two laws over one such source object are two values, and a geometry
holding one inherits the refusal, with the path through
``boundaries``.

⚠ The directional ``Quadrature`` is such a part today: it is a mutable
dataclass (``frozen=False``) with a hand-written ``_identity_key``, so
it compares and hashes by content in process but has no content
encoding (`[M]` ``encode(Quadrature.gauss_legendre(4))`` raises
``ContentlessError``). Nothing at this step keys on a quadrature; a
persistent key that must cover an S\ :sub:`N` angular discretisation
needs the quadrature moved onto the encoder first.

**NaN is refused at construction, by parsing at the boundary.** A value
holding NaN would have no equality to encode, so the refusal is placed
where the value is made, and ``==`` and ``hash`` never meet a NaN on a
value that could be constructed. The parsers, each a ``ValueError``
naming the field or parameter:

* :func:`~orpheus.numerics.scalars.parse_real`, the one definition of
  "a real number" at L1 for every layer that admits numbers
  (:ref:`structured-geometry-one-real-parser`): it refuses ``bool`` and
  every non-real type (``TypeError``), NaN (``ValueError``, its message
  fragment ``is NaN, which is not a number``) and an integer a double
  cannot carry, passes infinities (a caller needing a finite value asks
  :func:`~orpheus.numerics.scalars.parse_finite_real`), and returns
  ``+0.0`` for ``-0.0``;
* ``BC`` parameters: each name must be a ``str`` and each value goes
  through ``parse_real``, so a ``str`` or ``bool`` value is a
  ``TypeError`` and NaN a ``ValueError``; the parsed parameters are
  stored as a ``FrozenMapping``;
* ``AlbedoBoundary.albedo``, ``WhiteBoundary.albedo``,
  ``ConstantInflowSource.value``, ``LawScaled.scalar``, and the
  ``alpha`` of ``ScalarResponse``, ``LambertianReemission`` and
  ``SpecularReemission``, through ``parse_real``; the ``outward_sign`` of
  ``WhiteBoundary``, ``LambertianReemission`` and ``IsotropicReturn`` is
  an ``int``;
* ``Mesh2D``'s edges through ``parse_positions`` (an infinite or a
  ``bool`` edge is refused) and its material map through
  ``parse_integer``: a material id is an ``int``, so a float map, even
  ``[[1.0]]``, is refused, as in ``Mesh1D``;
* every array of a ``Mixture``: the dense fields, the stored values of
  every sparse Legendre block and the energy grid through
  :func:`~orpheus.numerics.scalars.canonical_reals`, the message naming
  the field and the entry (since ``1fde7b59``; a loop of its own before),
  and ``chi`` through the emission-spectrum law.

The encoder refuses NaN too, as the backstop for a bare value handed to
``encode``, and it does so through the same function, so the message is
the parsers' message. Signed zero is canonicalised at both places for the
same reason: the parsers store ``+0.0`` where ``-0.0`` was given, and the
encoder maps ``-0.0`` to ``+0.0`` whatever it receives; both spell the
fold once, in :func:`~orpheus.numerics.scalars.canonical_real`.

.. _structured-geometry-one-real-parser:

One definition of a real number, at L1
--------------------------------------

Every value a key covers admits numbers: a boundary law's albedo, a
geometry's breakpoints, a mesh's edges, a mixture's cross sections, a
table's entries, a question's offsets, and the encoder itself. Each
admission must agree with the encoder on what a number *is*, or a value
that constructs could carry bits the digest does not describe. That
agreement is one module, :mod:`orpheus.numerics.scalars`, which imports
nothing from ``orpheus`` and sits below the encoder that uses it (#559).
It moved from ``orpheus/geometry/scalars.py`` to ``numerics``, the
lowest layer whose vocabulary suffices, because the encoder (L1) needed
it and may not import the input layer; the geometry and the mesh import
L1, which the layer rule allows (:ref:`architecture-layering`).

**The rule** is written once, in
:func:`~orpheus.numerics.scalars.canonical_real`: NaN is refused (it is
not a number, and not equal to itself, so a value holding one has no
content equality), and ``-0.0`` becomes ``+0.0`` (they compare equal, so
they are one value). **The conversion** every admitting site makes is
:func:`~orpheus.numerics.scalars.exact_double`: the double that carries
a real scalar exactly as ``==`` sees it, which is the rule plus two
refusals,

- **an integer beyond** :math:`2^{53}` is refused. A double carries every
  integer up to :math:`2^{53}` exactly and no further: `[M]`
  ``float(2**53 + 1) == float(2**53)`` is ``True``, so admitting the
  integer would make two unequal integers one value, and ``==`` would
  separate the integer from the double stored for it. Before
  ``8122abc6`` the parsers and the encoder disagreed here (qa, `[M]`):
  ``parse_real(2**53 + 1)`` stored :math:`2^{53}` silently while the
  encoder refused the same integer;
- **a magnitude beyond the range of a double** is refused with its key.
  ``float()`` raises an unkeyed ``OverflowError`` on such a value; it is
  caught and re-raised as a ``ValueError`` naming the field. `[M]`
  ``parse_real(Fraction(10**400, 3), "x")`` reads ``x: Fraction(...) lies
  beyond the range of a double``. An integer that large is caught one
  check earlier, by the :math:`2^{53}` refusal.

Infinities are real and pass both functions; a caller that needs a
finite value asks for it (below). The functions, each adding one thing to
the one beneath it:

.. list-table:: The parsers of :mod:`orpheus.numerics.scalars`
   :header-rows: 1
   :widths: 26 74

   * - Function
     - What it adds
   * - ``canonical_real(value, where)``
     - the rule: NaN refused with the fragment ``which is not a number``,
       ``-0.0`` made ``+0.0``
   * - ``exact_double(value, where)``
     - the :math:`2^{53}` and overflow refusals, then ``canonical_real``;
       admits ``bool``
   * - ``canonical_reals(values, where)``
     - ``exact_double`` entrywise on an array (the integer extremes
       checked on the integers themselves, a refused entry named by its
       index): the encoder's array and sparse paths and ``Mixture``'s
       fields
   * - ``parse_real(value, where)``
     - refuses ``bool`` and every non-``Real`` type (``TypeError``), then
       ``exact_double``
   * - ``parse_finite_real(value, where)``
     - ``parse_real``, then refuses an infinity: the one spelling of
       "finite", used by the question's offsets and :math:`\tau`, the
       geometry's breakpoints and ``parse_positions``
   * - ``parse_finite_reals(value, where)``
     - ``parse_finite_real`` for every entry of an array of any rank, so a
       ``bool`` hidden in a list is refused as the scalar is; a read-only
       ``float`` copy; ``RegionwiseConstant``'s table
   * - ``parse_positive_real``, ``parse_integer``,
       ``parse_positive_integer``, ``parse_entries``, ``parse_positions``
     - the geometry's and the mesh's other inputs, each over the parsers
       above

`[M]` this pass, 13 production modules import the module (a ``git
grep`` of ``orpheus.numerics.scalars import`` under ``orpheus/``, its
known member ``content.py`` found): six boundary-law modules and
``structured_geometry.py`` in ``geometry``; ``partition.py`` and
``structured.py`` in ``mesh``; ``mixture.py`` in ``data``; and
``content.py``, ``mesh_free_function.py`` and ``question.py`` in
``numerics``.

**Why the encoder admits** ``bool`` **and the parsers refuse it.** The
two answer different questions. The encoder digests a value that already
exists, and its rule is that the digest follows ``==`` (the user's ruling
of 2026-10-02, :ref:`structured-geometry-content-identity-forms`):
``True == 1 == 1.0`` in Python, so the three encode alike, and
``Eigen(True) == Eigen(1)``. A parser decides what an *input* means at
the boundary where a value is made, and a ``bool`` given where a number
is expected (an albedo of ``True``, an offset of ``True``) is a mistake in
the call, not a number the caller meant; admitting it would turn the
mistake into a value. So ``exact_double`` admits ``bool`` and
``parse_real`` checks the type first.

**What it replaced.** Before #559 the rule had three spellings: the
encoder's private ``_real``, the geometry layer's ``parse_real``, and
``RegionwiseConstant``'s inline NaN, infinity and dtype check. The first
commit (``6b179059``) moved the module and routed the first and third
through it, and its claim that the rule was then written once was
premature: the encoder's array and sparse paths kept their own NaN
refusal and ``-0.0`` fold (qa's finding), which ``8122abc6`` routed through
``canonical_reals``. The elegance review found two more pairs in the
module itself, fixed in the same commit: "a real entry" was decided by
the coerced dtype in ``parse_finite_reals`` and by ``isinstance`` in
``parse_real`` (`[M]` ``parse_finite_reals([1.0, True])`` returned
``[1., 1.]`` while ``parse_positions`` refused the same list), and
"finite" was spelled four times. The archivist's pass found a fourth
NaN refusal, ``Mixture``'s own loop ("holds NaN, which is not a cross
section"), outside the witness's vocabulary; ``1fde7b59`` passed every
field of a mixture through ``canonical_reals`` and retired the loop.

**The witness**, ``tests/gates/numerics/test_one_real_parser.py``,
keeps the copies from coming back (X4: the mechanism that keeps several
admitting sites equal needs a test that fails when one diverges). `[M]`
this pass, 20 rows, 20 passed under ``-O``:

- by AST over every production file, the NaN fragment ``which is not a
  number`` is written at exactly 1 site, with the activation leg that the
  site found is ``scalars.py``; first red, on the tree before #559, 2
  sites wrote it and the encoder refused NaN with a message of its own;
- a second filter in another vocabulary: by AST, every ``raise`` under
  ``orpheus/`` whose message text says ``NaN`` sits in ``scalars.py``
  (k = 1 of 1), its first red the old ``Mixture`` raise;
- every admitting site (``canonical_real``, ``parse_real``,
  ``parse_finite_real``, ``parse_finite_reals``, ``encode``,
  ``RegionwiseConstant``, ``Mixture``) refuses NaN with the one message;
- the sites agree on ``-0.0``; a finite parse refuses an infinity; a
  ``bool`` entry is refused by every parse, the list spelling included;
- the parsers and the encoder agree on what a double carries
  (:math:`2^{53} + 1` refused by ``parse_real``, by ``encode`` and by an
  encoded array; ``10**400`` and an overflowing ``Fraction`` refused with
  their key); a rank-0 array parses.

.. warning::

   **Open: the frozen encoder admits a read-only view of a writeable
   array (#561).** The encoder refuses a writeable array as a part,
   because it could change after its owner was keyed, but it checks only
   the array's own ``writeable`` flag, and a read-only VIEW over a
   writeable base passes. `[M]` this pass: a frozen content value holding
   ``base.view()`` with the view's flag cleared is admitted; after
   ``base[0] = 1.0`` the value keeps its cached digest, and a fresh value
   over the same view no longer equals it. A content type that copies its
   arrays at construction (``Mixture``, ``Mesh2D``,
   ``RegionwiseConstant``) is safe by that copy; the hole is for a type
   that stores an array uncopied. The mapping analogue was closed in
   ``8122abc6``: a ``MappingProxyType`` part is refused, since it is a
   view over a ``dict`` its caller may still hold, and
   :class:`~orpheus.numerics.content.FrozenMapping` encodes its own items
   through a private view over storage only it holds, byte for byte as
   before (the RECORD pins did not move). The strict array check (a view
   refused unless it owns its data) reddens at least 19 owners, measured
   by the orchestrator, which is why it is an issue and not an edit.

.. _structured-geometry-content-identity-replaced:

What goes through the encoder, and what it replaced
---------------------------------------------------

Every space-name digest goes through the encoder: ``FunctionSpace.of_axes``, the
angular and scalar trace spaces, the two radial characteristic spaces
and ``FullFieldSpace.from_blocks`` each pass a tuple of their parts to
:func:`~orpheus.numerics.content.name_digest`, the one 8-byte space-name
digest (the first 8 bytes of ``content_digest``, in hexadecimal), and the
three trace spaces read their face layout's content through
``FaceLayout.structure``, the layout's ``(key, offset, size)`` triples; the
route gate S5.6 (below) reds any of them that hashes its own bytes.

**What is outside the encoder, and why.** The Sood registry's result cache
(``orpheus/derivations/continuous/sood_registry/cache.py``) keys its
entries on a SHA-256 of sorted JSON; `[M]` in the census that encoding
separates ``1`` from ``1.0`` and ``-0.0`` from ``0.0``, admits NaN,
identifies ``{1: x}`` with ``{"1": x}`` and a tuple with a list, and has
no production consumer (one test file reads it), so it is outside the
encoder. The in-process keys that compare by content within
one process and are never persisted are outside it too:
``MaterialMesh._identity_key`` and ``_contractibility_key``, the
``SNProblem`` extension of them, ``Quadrature._identity_key`` and
``RigidMotion._exact_key``.

.. dropdown:: The encoders and keys the encoder replaced, and who read equality (`[M]` at ``1dc31163``)
   :color: muted

   The explorer's census of the existing identity machinery
   (``scratch/reference_architecture/p1step5/census.md`` §1) used the
   predicate ``git grep -nE
   "_structural_bytes|_identity_key|blake2b|hashlib|sha256|digest|_law_key|_contractibility_key|_canonicalize"``
   over ``orpheus``, ``tests``, ``derivations`` and ``tools``: 415 lines,
   308 of them one JSON ledger, with ``Axis._structural_bytes`` as the
   positive control. It found five byte encoders that produced a
   cross-process digest. All five moved onto the encoder:

   .. list-table:: The five byte encoders, `[M]` at ``1dc31163``
      :header-rows: 1
      :widths: 26 28 46

      * - Encoder
        - What it encoded
        - What it lacked
      * - ``Axis._structural_bytes``, read by ``FunctionSpace.of_axes``
          for the space name
        - the class's ``__qualname__``, then label, shape, kind and the
          weights' bytes (``EnergyAxis`` added the edges,
          ``LegendreAxis`` the spent axis)
        - it was type-tagged and length-prefixed, but encoded a scalar by
          its ``repr``, so ``1`` and ``1.0`` gave different bytes, and its
          class tag had no module and no version
      * - ``AngularTraceSpace.for_layout``
        - the ``repr`` of the face layout, then the raw bytes of
          :math:`\Omega\cdot\hat n`, the quadrature weights and the nodes
        - no type tag, no length prefix, no shape or dtype, ``-0.0`` not
          canonicalised, NaN not refused
      * - ``ScalarTraceSpace``'s name mint
        - the ``repr`` of the layout, then the raw weight bytes
        - the same; its comment said it mirrored the angular mint
      * - ``RadialCharacteristicSpace``'s name mint
        - the ``repr`` of the layout and of the levels, then the raw metric
          bytes
        - the same; the three trace mints were one encoder written three
          times
      * - ``FullFieldSpace.from_blocks``
        - the text ``"<name>:<shape>|<name>:<shape>"`` of its two blocks
        - a text join with no tags

   The in-process keys that duplicated a type's equality retired with them:
   ``Axis._identity_key`` and ``Axis``'s hand-written ``__eq__`` and
   ``__hash__`` (the dataclass fields replace them, with the generator
   declared ``compare=False``),
   ``Mixture._identity_key`` with its ``__eq__`` and ``__hash__``,
   ``CellEdges``' bytes hash, ``Mesh1D``'s hand-written ``__eq__``, and the
   ``BC`` branch of ``MaterialMesh``'s law key, which keyed a tag by
   ``("BC", kind, sorted params)`` because a tag could not be hashed.

   **Who reads equality and hash in production**, from the census's
   runtime spy (§3a: every ``__eq__`` and ``__hash__`` of the affected
   classes wrapped, 44 wrappers on 23 classes, over ``tests/gates`` with
   ``-m "not slow"``, 23 of 24 shards; ``tests/gates/derivations`` was
   excluded because it had not finished): the S\ :sub:`N` geometry intern
   (keyed on the problem and its contractibility key; 1984 hash and 1538
   equality calls, of which 44 reached a law), ``MaterialMesh``'s law key
   (44), collision probability's law refusal (109 tuple-membership
   comparisons), ``Mesh1D`` equality through ``FaceLaws`` (3) and
   ``HomogeneousProblem`` equality over its mixture (13). No production
   code used a ``Mixture``, ``Materials``, ``BC``, ``StructuredGeometry``,
   ``Mesh1D`` or ``Mesh2D`` as a dictionary key. Two production keys
   changed meaning with the move: ``MaterialMesh``'s contractibility key
   folds each mixture's canonical digest, so ``-0.0`` and an explicitly
   stored zero no longer separate two mixtures there, and its law key
   keys a ``BC`` tag by the tag itself.


.. _structured-geometry-content-identity-gates:

The gates and the fingerprint
-----------------------------

The specification of the gates is
``.claude/plans/reference_p1_spec.md`` §1.5 (re-specified 2026-10-02 by
the test-architect, with its mutation battery). Every gate is a
``foundation`` test (a software invariant of the encoder and the types,
with no theory-page label), and each carries ``rests_on``. `[M]`
2026-10-02, before the step's review fixes, ``.venv/bin/python -O -m
pytest`` over the five files: 253 collected, 253 passed; the review
fixes added rows, and their run of record is the full suite.

.. list-table:: The content-identity gates
   :header-rows: 1
   :widths: 40 10 50

   * - File
     - Rows at the first pass
     - What it holds
   * - ``tests/gates/numerics/test_content_identity.py``
     - 45
     - the encoder and the type population: seed stability (S5.1:
       digests, hashes and a space name printed by two interpreters
       under ``PYTHONHASHSEED`` 1 and 2 agree, with a ``str`` hash as
       the control that the harness sees salting); the canonical forms
       (S5.4); the route gate (S5.6: a decoy rebound over the
       encoder's one recursive entry must move the digest and the hash
       of every roster type and every derived space name, so a type or
       a space still hashing its own bytes reds as a second encoder);
       the fingerprint (S5.7); the schema tag (S5.8: a field added,
       renamed or reordered, a version raised, a class moved or
       renamed each move the digest); conformance (S5.9: every
       content class takes ``__hash__`` from the mixin and ``__eq__``
       from it or a declared override, a dataclass is ``eq=False`` and
       frozen, and every concrete content class has a roster entry)
   * - ``tests/gates/numerics/test_content_identity_axis.py``
     - 32
     - the ``Axis`` family: equal content is one value, a moved part
       moves the digest, pickle round trips
   * - ``tests/gates/data/test_content_identity_data.py``
     - 38
     - ``Mixture`` and ``Materials``: equal content (S5.2), each part
       perturbed (S5.3, the population read off the production
       instance), ``-0.0`` and NaN (S5.4), ``Materials`` key coercion
       and order, a pickled ``MaterialMesh``
   * - ``tests/gates/geometry/test_content_identity_geometry.py``
     - 101
     - ``BC``, the seven registered laws and their parts, ``LawSum``,
       ``LawScaled`` and ``StructuredGeometry``: S5.2 to S5.4, the
       read-only ``BC.params``, the tag that is not its law, the kept
       string arm, the NaN and non-real refusals at construction, and
       the contentless law (S5.5)
   * - ``tests/gates/mesh/test_content_identity_mesh.py``
     - 37
     - ``FaceLaws``, ``CellEdges``, ``Mesh1D`` and ``Mesh2D``: S5.2,
       S5.3, ``FaceLaws`` is not a ``dict``, and ``Mesh2D``'s read-only
       copies, canonical signed zero and equality (S5.10)

The population is quantified, not listed: each file declares a roster
of its types with the expected part names, a perturbation per part and
pairs of equal content. The population row of S5.3 reads the part names
off the production instance and requires the roster to perturb exactly
those, so a part added later reddens until it is perturbed; S5.9 closes
the population from the other side, so a new content class reddens
until it joins a roster. The battery of the specification ran 41
textual arms against the encoder with a positive control (every chunk
emptied: 48 new reds), and read each arm's red set against its target
row.

**The fingerprint.** S5.7 is a RECORD: it pins the hexadecimal digests
of three values built from literals, so any edit to the encoder or to
the schema of these classes reddens it on purpose. It was re-pinned once
within the step: ``FrozenMapping`` changed the schema of the geometry's
``BC`` parameters and of ``Materials``, so those two digests moved, while
mixture A's did not. The three values are
the 2-group mixture A of ``xs_library``, a hollow sphere of three
intervals with an albedo specular inner law and a vacuum outer law, and
a ``Materials`` of the 2-group mixtures A and B. Every cache key of the
cache to come is a digest of this kind, so a moved digest invalidates
every stored entry; the message says so, and a re-pin carries its
reason in the commit message. The inputs are literals and IEEE-exact
arithmetic, with no output of libm or LAPACK, so the pinned bytes are
the same on every platform (V&V anti-pattern 38, never pin a platform
library's float output). The ``gates`` CI workflow runs the five files
under ``-O`` on its Linux runner, a step named "The content digests are
the same bytes on this platform", so a platform difference in the
encoding reddens there against the pins measured on the macOS host.


.. _structured-geometry-mesh-free-functions:

Mesh-free functions: the source and the detector a specification states
=======================================================================

A reference question (the reference-solution cache of #405) is stated
before any mesh exists: the materials, the geometry with its laws, and
the source and the detector of a fixed-source problem. The source and
the detector therefore cannot be arrays on a discretised space; they are
functions stated *intensionally*, by a rule a later projection evaluates
on whatever unknowns a method chooses. The module
:mod:`orpheus.numerics.mesh_free_function` holds the two forms such a
function takes, and it landed as step 6 of the campaign's first phase,
on 2026-10-02 (merged to ``main`` with the step's close-out,
``61f82a18``):

- :class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant`, one
  real value per (region, energy group): a function on the
  **angle-integrated** space;
- :class:`~orpheus.numerics.mesh_free_function.Symbolic`, one SymPy
  expression :math:`q_g(r,\mu,\varphi)` per energy group: a function on
  **phase space**, stored as text.

Both are content-identity values
(:ref:`structured-geometry-content-identity`), so a cache key built
from a question that holds one is stable across processes. Neither
carries a role (source or detector) nor a density convention; the field
of the specification that holds the value is the role, and the role
picks how the value enters phase space.

**Why this page is the home.** The two types are mathematics and live
in ``numerics`` (L1; a region is an index into a partition, a group an
index into an axis, and the measure on the direction sphere is
mathematics), but what they *mean* is defined here: a
``RegionwiseConstant``'s regions are this page's interval indices
(:ref:`structured-geometry-value`), a ``Symbolic`` reads its direction
in the chart each :class:`~orpheus.geometry.coord.CoordSystem` declares
(below), and both enter the question whose content identity is the
chapter above. The arrows that carry them into phase space are the
space layer's, and are derived there
(:ref:`spaces-collapse-pair-two-lifts`); this chapter points at them and
does not re-derive them. The specification of the step, with its gates
and its measured first reds, is ``.claude/plans/reference_p1_spec.md``
§1.6; the census that preceded it is
``scratch/reference_architecture/p1step6/census.md``.

.. _structured-geometry-mesh-free-two-types:

Two types, because a per-region table is not on phase space
-----------------------------------------------------------

The first specification (2026-09-25) made ``RegionwiseConstant`` *the
isotropic, piecewise-constant special case of* ``Symbolic``: a table
value :math:`Q` would lower to the phase-space function
:math:`Q/4\pi`, a density per steradian fixed on the type. That design
is refuted (`[REFUTED 2026-10-02]` for the question *"what is a
per-region table, as a function?"*), for two reasons that are one
reason.

**The** :math:`4\pi` **is the measure, not a convention.** The user's
ruling of 2026-10-02 (``.claude/plans/reference_cache.md``, "P1 step 6
opened", ruling 1): "integrated over all directions" and "divided by
:math:`4\pi`" are both statements about the angular measure, and the
operators that perform them already exist. A per-region table is a
function on the space the angular retraction :math:`R` maps *to*, the
angle-integrated space (region, group). It enters phase space by one of
two arrows: a source rate through the section :math:`E`, so that
:math:`R(EQ) = Q`; a detector through the adjoint :math:`R^{\dagger}`,
the Riesz representative of :math:`\psi \mapsto \langle\Sigma_d,
R\psi\rangle`. The two differ by the mass of the measure,
:math:`R\circ R^{\dagger}`, which is :math:`4\pi` on the sphere and 2 on
the orbit space of a one-dimensional rule, and which is therefore never
written into a type.

**One lowering would be wrong for one of the two roles.** `[M]` the
census and the test-architect's SymPy probe
(``scratch/reference_architecture/p1step6/ta/symbolic_probe.out``, S4):
a detector lowered by the source's rule, :math:`\Sigma_d/4\pi`, gives
:math:`1/(4\pi)` of the true response on the continuous sphere, and the
same rule applied with the slab's weight sum 2 gives
:math:`R(Q/2) = 2\pi Q` for a source. A type fixing one density would
have made every ``RegionwiseConstant`` used as a detector wrong by a
factor, and a hand-written :math:`1/4\pi` wrong on the slab and the
sphere, where the discrete measure's mass is 2.

Three consequences are structural:

- **No role and no density is a field** of either type (gate S6.11:
  ``dataclasses.fields(RegionwiseConstant)`` is exactly ``values``).
  The alternative the census offered, "the value type declares which
  density it carries", is the spelling this rules out.
- **The two types are never equal to each other.** Two content types
  are never equal across types (:ref:`structured-geometry-content-identity`),
  so the step-8 leg "a ``RegionwiseConstant`` equals the ``Symbolic`` it
  lowers to" cannot hold as written; it is re-posed as "both lift to the
  same function" (spec §1.8).
- **A** ``Symbolic`` **needs no lift onto the sphere**, because it is
  already a density with respect to :math:`\mathrm d\Omega`; onto a rule
  whose ordinates are orbits it needs a pushforward (below).

.. _structured-geometry-regionwise-constant:

``RegionwiseConstant``: one value per region and group
------------------------------------------------------

``values`` is a real table of shape ``(regions, groups)``, at least one
of each, stored as a **read-only float copy**: the caller's array can
change afterwards without moving the value or its digest (gate S6.8).
Region :math:`k` is the geometry's interval :math:`[r_k, r_{k+1}]`, a
positional index, ``0 … len(mat_ids) − 1``; a region is not a material
(two regions holding one material are two rows). The reference
specification checks the row count against the geometry's intervals and
the column count against the materials' group count (S6.12;
:ref:`structured-geometry-specification-datum`).

The constructor refuses, each with its own message: a table of rank
other than 2 (naming the rank), an empty region or group axis (naming
which), a NaN entry (naming its index, the encoder's rule), an infinite
entry (naming its index: a rate or a response of infinite size is not a
function value), and a non-real or boolean dtype. Content identity
follows the encoder's canonical forms (gate S6.9): ``-0.0`` is ``0.0``,
an integer table equals its float twin, and a ``(2, 3)`` and a
``(3, 2)`` table with the same bytes are different values, because the
shape is part of the content. ``n_regions`` and ``n_groups`` read the
shape.

The finite-real check is :func:`~orpheus.numerics.scalars.parse_finite_reals`,
the one definition of "a finite, canonical real" at L1 (#559), which the
geometry, the mesh, the data layer and the content encoder share
(:ref:`structured-geometry-one-real-parser`).

.. _structured-geometry-symbolic:

``Symbolic``: a function on phase space, stored as text
--------------------------------------------------------

**The coordinates are owned.** ``Symbolic.r``, ``Symbolic.mu`` and
``Symbolic.phi`` are class-level ``Symbol(name, real=True)``: the
position, the cosine to the chart's polar axis and the azimuth about it
(the chart is the next section). The position coordinate is :math:`r`
on every chart, the slab's :math:`x` included. ``real=True`` and not
``nonnegative`` for :math:`r`, because a slab admits any :math:`r_0`.
The symbols are built on first access through a descriptor,
``_OwnedSymbol``, which takes its name from ``__set_name__``, so
importing the module builds no SymPy object.

An expression with any other free symbol is refused, and the message
names the symbol **and its assumptions** (gate S6.13). A symbol whose
name is an owned name but whose assumptions differ is a different symbol
to SymPy (``Symbol('mu', positive=True, real=True) == Symbol('mu',
real=True)`` is ``False``), and it gets its own message, *"same name,
different assumptions: use Symbolic.mu (real=True)"*. That case is not
hypothetical: the tree's own manufactured-solution builder declares
:math:`r` and :math:`\mu` with ``positive=True``
(``orpheus/derivations/continuous/mms/sn.py``), and the transport
equation's derivation module declares :math:`r` ``nonnegative``. A check
by name would have admitted them and produced a function in coordinates
the type does not own; one Branch-1 coordinate vocabulary is part of
#557.

**Storage is text, and identity is by spelling.** ``srepr`` holds one
canonical ``sympy.srepr`` string per group, and the content is that
text (gate S6.16: a live SymPy expression as a content part is refused
by the encoder). Two spellings of one function are two values: `[M]`
``Symbolic.of((r + 1)**2)`` and ``Symbolic.of(expand((r + 1)**2))`` are
unequal, while ``simplify`` reduces their difference to 0. That costs
a cache miss and never a wrong hit, which is the only direction a cache
key may err in; deciding equality of two expressions in general is not
something a key can afford. `[M]` the round trip
``Symbolic.from_srepr(s.srepr) == s`` holds, with assumptions kept, on
the gate's four-member population (a polynomial; :math:`\sin(\pi r)
e^{-\mu} + \cos\varphi\,\sqrt{1-\mu^2}\,r`; a ``Piecewise`` with a
30-digit ``Float`` and a ``Rational``; a 12-term sum), and the digest is
the same in two interpreters under ``PYTHONHASHSEED`` 0 and 1 (S6.14).

**The SymPy version is content** (gate S6.17), because the two ways a
new SymPy can disagree with a stored string differ in kind. A SymPy
that *writes* a different ``srepr`` for the same expression gives a
different key: a miss, which costs a recomputation. A SymPy that
*parses* a stored string into a different expression gives the same key
for a different function: a wrong hit, which returns a wrong answer. A
pinned record of the text detects both, but only when the suite runs
under the new version; a user who upgraded SymPy and ran nothing would
be unprotected against the wrong hit. With ``sympy_version`` a content
part, every read compares versions, so the wrong hit is impossible, at
the price of a miss on every SymPy upgrade. The full version string is
used, not major.minor, because nothing measured bounds a patch
release's effect on ``srepr``. The RECORD pin stays beside it as the
early notice that the format moved (its message: *"SymPy changed its
srepr format: every Symbolic cache key is invalidated (the version part
already forces the miss); re-pin with the version"*).

**Stored text is parsed through a whitelist, never by** ``eval``
**alone.** The first ``from_srepr`` evaluated stored text with
``sympify``, which calls ``eval``; the qa review (2026-10-02) wrote a
file with a stored string. A cache is a store of text a later process
reads back, so the parse is an input boundary. The text is parsed to a
Python AST first, and only these nodes are admitted: a call whose callee
is the name of a SymPy class (every subclass of ``sympy.Basic`` by class
name, since ``ExprCondPair`` is not a top-level export) or a SymPy
singleton (``pi``, ``true``, ``oo``); numeric, string and boolean
literals; a unary sign on a literal; keywords and tuples. Anything else
is refused before evaluation, naming the node or the name (an attribute,
a subscript, a lambda, ``getattr``, ``Matrix``, a syntax error; gate
S6.16's six rows), and the admitted tree is evaluated with no builtins.
`[M]` ``Symbolic.from_srepr(("__import__('os')",))`` is refused:
the message names ``'__import__'`` as "not a SymPy class or
constant".

**Only real scalar functions are admitted** (the qa review found each
of these admitted by the first build). Refused, each with its own
message: an object that is not a scalar expression (a relation
``r > 1``; a matrix); ``nan``, ``zoo``, ``oo`` or ``-oo`` anywhere in the
expression; the imaginary unit; an undefined function such as
:math:`f(r)`; a ``Piecewise`` with no otherwise branch, which has no
value outside its conditions. None of these is a real function value
everywhere on phase space.

**Isotropy is decided by substitution.** ``is_isotropic`` asks whether,
in every group,

.. math::

   q_g(r,\mu,\varphi) - q_g(r,\mu',\varphi') \;\overset{\text{simplify}}{=}\; 0

with :math:`\mu',\varphi'` fresh real symbols: the value does not change
when the direction does. A difference that ``simplify`` cannot reduce
counts as a dependence, so the undecided case falls on the anisotropic
side, which is the side a consumer refusing anisotropy refuses. Two
earlier predicates failed, each on a row the gate now carries (S6.15,
ten rows):

- *by free symbols* (anisotropic iff :math:`\mu` or :math:`\varphi`
  appears): wrong on :math:`\sin^2\varphi + \cos^2\varphi` and on
  :math:`\mu^2 + (1-\mu^2)\cos^2\varphi + (1-\mu^2)\sin^2\varphi`
  (:math:`\Omega\cdot\Omega`), which are constant;
- *by derivatives* (isotropic iff :math:`\partial_\mu q` and
  :math:`\partial_\varphi q` simplify to 0): wrong on a **step** in the
  direction, ``Piecewise((1, mu > 0), (0, True))``, whose derivative is
  0 wherever it is defined (the qa review). Substitution sees the step,
  because :math:`q(\mu=0.5) - q(\mu'=-0.5) = 1`.

`[M]` this pass: :math:`1+\mu`, :math:`\cos\varphi` and the two steps in
:math:`\mu` and :math:`\varphi` are anisotropic; :math:`2`,
:math:`r^2+1`, :math:`\mu-\mu`, the two trigonometric traps and a step in
:math:`r` are isotropic.

.. code-block:: python

   import sympy as sp
   from orpheus.numerics.mesh_free_function import RegionwiseConstant, Symbolic

   r, mu, phi = Symbolic.r, Symbolic.mu, Symbolic.phi
   q = Symbolic.of(sp.exp(-r) * (1 + mu), sp.Integer(2))   # two groups
   assert q.n_groups == 2 and not q.is_isotropic
   assert Symbolic.from_srepr(q.srepr, q.sympy_version) == q
   assert Symbolic.of(sp.sin(phi)**2 + sp.cos(phi)**2).is_isotropic

   table = RegionwiseConstant([[1.0, 0.5], [0.0, 2.0], [1.0, 0.5]])  # 3 regions, 2 groups
   assert (table.n_regions, table.n_groups) == (3, 2)

**From a density on the sphere to a rule's ordinates: the pushforward
(not built; phase P4).** A ``Symbolic`` is a density with respect to
:math:`\mathrm d\Omega = \mathrm d\mu\,\mathrm d\varphi`. On a rule
whose ordinates are points of the whole sphere (``folded_product``,
Lebedev) its value at an ordinate is :math:`q` itself. A one-dimensional
rule is different: its ordinates are points of the **orbit space** of
:math:`\mu`, the sphere quotiented by the rotations and reflections
about the polar axis (:ref:`manifold-orbit-space`), and the measure on
that space has mass 2. The per-ordinate value there is the integral of
:math:`q` over each orbit, the pushforward along the quotient map,

.. math::

   q^{\flat}(r, \mu) \;=\; \int_0^{2\pi} q(r, \mu, \varphi)\,\mathrm d\varphi ,

which is :math:`2\pi q` when :math:`q` does not depend on
:math:`\varphi`. That factor is the ratio of the two measures' masses,
:math:`4\pi/2`, and it is derived from the two measures when the arrow
is built, never typed. The arrow belongs to the projection of a
specification onto a method's unknowns, phase P4
(``.claude/plans/reference_cache.md``, "P1 step 6 built", the P4
obligations). The elegance review of 2026-10-02 named the hazard in the
first draft of the module docstring, which said a ``Symbolic`` "takes no
lift" and called the mass "4π on the sphere" beside a spherical S\
:sub:`N` problem whose :math:`\Sigma w` is 2: the ERR-004 / ERR-051
class, one step ahead of its first consumer.

**SymPy is imported inside the functions that need it.** `[M]` in a
fresh interpreter, importing the module and building a
``RegionwiseConstant`` leaves ``sympy`` out of ``sys.modules`` (gate
S6.21, with the positive control that building a ``Symbolic`` loads
it). SymPy moved from the ``test`` and ``docs`` extras into
``[project].dependencies`` in the same step: `[M]` (the test-architect's
census, ``scratch/reference_architecture/p1step6/ta/deps_census.out``) of the 8 third-party
top-level modules imported anywhere under ``orpheus/`` (366 files, by
AST), 7 were declared and SymPy was not, although
``Quadrature.gauss_legendre`` imports it at run time (gate S6.20, which
prints its input count and controls on ``numpy``).

.. _structured-geometry-angular-chart:

The angular chart each coordinate system declares
-------------------------------------------------

A function :math:`q(r,\mu,\varphi)` reads its direction in a chart, and
the chart is the coordinate system's local frame, declared once by the
coordinate system (the user's ruling of 2026-10-02, ruling 2). The
declaration is :attr:`CoordSystem.angular_chart
<orpheus.geometry.coord.CoordSystem.angular_chart>`, an
:class:`~orpheus.geometry.coord.AngularChart` over the three columns of
the local orthonormal frame at the position :math:`r`, the same columns
a quadrature's ordinates carry (``Quadrature.axis_cosines(k)``):

.. math::

   \Omega \;=\; \mu\,\hat e_{\rm polar}
   + \sqrt{1-\mu^2}\,\bigl(\sin\varphi\,\hat e_\perp
   + \cos\varphi\,\hat e_{\rm ref}\bigr),

so :math:`\mu = \Omega\cdot\hat e_{\rm polar}` and :math:`\varphi` is
the azimuth about the polar axis, measured from :math:`\hat e_{\rm ref}`
towards the remaining column :math:`\hat e_\perp`. The polar axis is the
one direction a one-dimensional position distinguishes.

.. list-table:: The declared charts (3 of 3 members; gate S6.18)
   :header-rows: 1
   :widths: 16 26 26 32

   * - Coordinate system
     - Polar axis (column)
     - Azimuth reference (column)
     - The column left perpendicular
   * - ``CARTESIAN`` (slab)
     - :math:`\hat e_x` (0)
     - :math:`\hat e_z` (2)
     - :math:`\hat e_y` (1); a 1-D rule has no such column
   * - ``CYLINDRICAL``
     - :math:`\hat e_r` (0)
     - :math:`\hat e_z` (2), the cylinder's axis
     - :math:`\hat e_\varphi` (1), the column the quadrature names
       :math:`\xi = \Omega\cdot\hat e_\varphi`
   * - ``SPHERICAL``
     - :math:`\hat e_r` (0)
     - none
     - —

**The sphere declares no reference, and cannot.** A reference
direction perpendicular to :math:`\hat e_r`, chosen at every position
of the sphere, would be a continuous tangent vector field on the
sphere that vanishes nowhere, and the hairy-ball theorem says no such
field exists: every choice is singular somewhere. So
``azimuth_reference`` is ``None`` on the sphere, and a function that
depends on :math:`\varphi` has no well-defined value beside a spherical
geometry. The orchestrator's ruling on the specification's third open
question: such a ``Symbolic`` is refused when a specification pairs it
with a spherical geometry, keyed on the dependence predicate restricted to
:math:`\varphi`, ``Symbolic.depends_on(phi)``. That refusal belongs to the
specification, which reads the chart's ``azimuth_reference`` for it
(:ref:`structured-geometry-specification-datum`). On the slab and the
cylinder a :math:`\varphi`-dependent function is admitted.

**Whether a problem can see the azimuth is not the chart's to say.**
The first build gave the chart a third field, ``azimuth_observable``,
true on the cylinder only. The qa review retired it the same day, for
two reasons: it was false for a two-dimensional Cartesian mesh, which
uses ``CoordSystem.CARTESIAN`` and whose rules (``level_symmetric(4)``,
``product(4, 8)``) carry azimuthal information; and nothing read it.
Observability is a property of the problem's symmetry, and a quadrature
already states it as the orbit space its ordinates live on (a slab's
:math:`\mu` rule is declared on the sphere quotiented by
:math:`O(2)` about the polar axis). The projection (P4) reads it from
there.

**The declaration and the quadrature are one definition** (X4, gate
S6.19): the columns the sweeps read as :math:`\mu` (``mu_x`` on the slab
and the sphere, ``eta`` on the cylinder) are ``array_equal`` to
``axis_cosines(chart.polar_axis)`` on one S\ :sub:`N`-admitted rule per
coordinate system (Gauss–Legendre 8; Gauss–Legendre 8;
``folded_product(4, 8)``), and on the cylinder the column the chart
leaves perpendicular is the one the rule names :math:`\xi`. `[M]` this
pass, the cylinder's ``folded_product(4, 8)`` has :math:`\xi > 0` on 16
of 16 ordinates and :math:`\varphi = \operatorname{atan2}(\xi, \Omega
\cdot\hat e_z) \in (0.222, 2.919)`: the rule is folded about the plane
of :math:`\hat e_r` and :math:`\hat e_z`, so on the cylinder the azimuth
is observable modulo the reflection :math:`\varphi\mapsto-\varphi`.

.. note::

   **A gate that was designed green, and its re-posing.** The
   specification's S6.19(b) re-synthesised the cylinder's two off-polar
   columns from :math:`(\mu, \varphi)`, with :math:`\varphi` measured
   from the declared reference, and asserted they matched to 2 ULP.
   `[REFUTED 2026-10-02]` for the question *"does this pin the declared
   reference?"*: :math:`\varphi = \operatorname{atan2}(\Omega\cdot\hat
   e_\perp, \Omega\cdot\hat e_{\rm ref})` rebuilds the two columns it
   was computed from for *any* choice of reference, so the mutation
   "reference = column 1" left the row green (the step's mutation
   battery, ``scratch/reference_architecture/p1step6/impl/mut_6c.py``).
   The row now asserts the quadrature's own naming, that the
   perpendicular column is ``Quadrature.xi``, and the same mutation
   reddens it.

.. _structured-geometry-two-lifts-branch-1:

Branch 1: the continuous measure and its two lifts
--------------------------------------------------

The discrete arrows are the space layer's
(:class:`~orpheus.numerics.operator.AxisSectionOperator`,
:class:`~orpheus.numerics.operator.AxisPullbackOperator`). Their
continuous counterpart, written so that the measure's mass is derived
and never typed, is the algebra of record
:mod:`orpheus.derivations.common.angular_measure` (Branch 1, closed-form
SymPy). It reads no quadrature and no production operator, so it is
structurally independent of the discrete arrows (X4).

The direction sphere :math:`S^2` carries
:math:`\mathrm d\Omega = \mathrm d\mu\,\mathrm d\varphi` over
:math:`\mu\in[-1,1]`, :math:`\varphi\in[0,2\pi)`. The integration
domain is written once, as the constant ``SPHERE``, which is the only
place the module spells :math:`\pi`. On it the module defines

.. math::

   R\,q = \int_{S^2} q\,\mathrm d\Omega,
   \qquad
   m = R\,1,
   \qquad
   E\,Q = \frac{Q}{m},
   \qquad
   R^{\dagger}\Sigma = \Sigma ,

the retraction, the mass of the measure, the section and the pullback.
SymPy computes :math:`m = \int_{-1}^{1}\int_0^{2\pi}\mathrm d\varphi\,
\mathrm d\mu = 4\pi`. Two verification functions prove the two lifts,
each pinned by a ``foundation`` test in
``tests/gates/derivations/test_angular_measure_symbolic.py``:

- **V_E, the section keeps the rate**
  (``derive_section_identity``): :math:`R(E\,Q) - Q` simplifies to 0 for
  every entry of a three-region, two-group table of ``Rational`` values,
  :math:`R(Q/m) = (Q/m)\,R\,1 = Q`.
- **V_R†, the pullback is the adjoint of the retraction**
  (``derive_adjoint_identity``), with :math:`\Sigma(r) = 1 + r` and a
  :math:`\psi` anisotropic in both angles,
  :math:`\psi = (1 + \mu + \mu^2\cos\varphi)\,e^{-r}`, on
  :math:`r\in[0,1]`:

  .. math::

     R\psi = e^{-r}\!\int_{-1}^{1}\!\!\int_0^{2\pi}
             (1 + \mu + \mu^2\cos\varphi)\,\mathrm d\varphi\,\mathrm d\mu
           = e^{-r}\,(4\pi + 0 + 0) = 4\pi e^{-r},

  because :math:`\int\mu\,\mathrm d\mu` and
  :math:`\int\cos\varphi\,\mathrm d\varphi` vanish over the domain; then

  .. math::

     \int_0^1\!\Sigma\,R\psi\,\mathrm dr
     = 4\pi\!\int_0^1\!(1+r)\,e^{-r}\,\mathrm dr
     = 4\pi\Bigl[-(2+r)\,e^{-r}\Bigr]_0^1
     = 8\pi - \frac{12\pi}{e},

  and the left side, :math:`\int_0^1\!\int_{S^2}(R^{\dagger}\Sigma)\,
  \psi\,\mathrm d\Omega\,\mathrm dr`, evaluates to the same value
  (`[M]` this pass, SymPy 1.14.0: both sides
  :math:`8\pi - 12\pi e^{-1}`). The function also returns the ratio a
  detector lifted by the **wrong** arrow would read,
  :math:`\int\!\!\int (E\Sigma)\,\psi / \int\Sigma\,R\psi = 1/m =
  1/(4\pi)`, and the test asserts it: the wrong arrow is measured, not
  assumed.

Two further rows: the mass is :math:`4\pi` (the comparison constant is
written in the test, as the independent value), and an AST pass over the
module finds ``pi`` only inside ``SPHERE`` (a typed ``4*pi`` or
``1/(4*pi)`` would add a site, and the step's mutation battery reddened
the row with one). Branch 1 here is State 1A of the
algebra-of-record discipline: the identities close in elementary
functions, so no ``mpmath`` stage is needed.

**The seed of #557.** The manufactured-solution builders divide by a
hand-read ``quadrature.weights.sum()`` (12 of 12 S\ :sub:`N` cases) and
the Green's-function multi-region sphere divides by a typed
:math:`4\pi`; #557 asks them to lift through one Branch-1 measure object
instead, and this module is that object's first instance. The elegance
review's direction for it: a measure *value* parameterised by its
domain, with the sphere and the orbit space of :math:`\mu` (mass 2) as
its two instances, rather than free functions over the one ``SPHERE``
constant; and one coordinate vocabulary shared with ``Symbolic``.

.. _structured-geometry-mesh-free-gates:

The gates, and what was refuted on the way
-------------------------------------------

`[M]` 2026-10-02, ``pytest --collect-only`` over the step's files: 110
rows.

.. list-table::
   :header-rows: 1
   :widths: 44 8 48

   * - File
     - Rows
     - What it holds (spec §1.6 row ids)
   * - ``tests/gates/numerics/test_mesh_free_function.py``
     - 53
     - S6.8 construction and refusals; S6.11 no role or density; S6.12
       the group count; S6.13 the owned symbols and the stray-symbol
       refusals; S6.14 the round trip and process stability; S6.15 the
       isotropy predicate; S6.16 storage as text, identity by spelling,
       the non-values, real scalar admission, the parse whitelist; S6.17
       the version as content and the RECORD pin; S6.21 the lazy import
       and no geometry import (by AST, with a positive control on the
       module's own ``orpheus.numerics.content`` import)
   * - ``tests/gates/numerics/test_content_identity_mesh_free.py``
     - 14
     - S6.9 and S6.16: both types join the content rosters (equal content
       is one value, a moved part moves the digest, pickle round trips),
       so the encoder's conformance gate (S5.9) sees them
   * - ``tests/gates/geometry/test_angular_chart.py``
     - 8
     - S6.18 the declaration on 3 of 3 members and the chart's two
       refusals; S6.19 the tie to the quadrature's columns
   * - ``tests/gates/derivations/test_angular_measure_symbolic.py``
     - 4
     - S6.10 the two lifts in Branch 1, the derived mass, ``pi`` only in
       the domain
   * - ``tests/gates/test_dependencies_declared.py``
     - 1
     - S6.20 every third-party import under ``orpheus/`` is declared
   * - ``tests/gates/numerics/test_retraction_adjoint_is_the_pullback.py``
     - 28
     - S6.1 to S6.4 and the axis verbs
       (:ref:`spaces-collapse-pair-pullback`)
   * - ``tests/gates/sn/solve/test_detector_lift_is_the_retraction_adjoint.py``
     - 2
     - S6.6 the S\ :sub:`N` detector lift's route (:ref:`sn-adjoint-dual-lift`)

The mutation batteries ran in process, one arm per "mutation witness"
of the specification, each read against its target row
(``scratch/reference_architecture/p1step6/impl/mut_6*.py``): SymPy
undeclared; the detector lifted by the section; a typed ``4*pi``; stray
symbols checked by name; isotropy by free symbols; the version dropped
from the content (step ``f7309da3``); then, after the review, the
derivative predicate (the two step rows red), a bare ``sympify`` (the
six whitelist rows red), the real-scalar checks removed (11 rows red).

**Refuted on the way, each with its structural reason** (the
specification's "Corrections found while building" and the review
reports carry the measurements):

- ``RegionwiseConstant`` as ``Symbolic``'s isotropic case, lowered by
  :math:`Q/4\pi`: a detector would read :math:`1/(4\pi)` of its response,
  and a typed :math:`4\pi` is wrong where the discrete mass is 2
  (above).
- A density convention on the type: the role, not the value, picks the
  arrow (S6.11).
- Isotropy by free symbols, then by derivatives: the trigonometric
  traps, then the steps (above).
- ``from_srepr`` through ``sympify``: stored text ran code.
- The stray-symbol check by name: the tree's own ``positive=True``
  coordinates passed it.
- ``azimuth_observable`` on the chart: false for a 2-D Cartesian mesh,
  and unread.
- S6.19(b), re-synthesising the cylinder's columns: designed green.
- The module name ``phase_space_function``: a per-region table is not a
  function on phase space, so the module is named after what both types
  share, that they are stated before a mesh (renamed before any page
  referred to it).

**Not built, and where it lands.** The pushforward onto an orbit space is
P4. The specification's fields, with the refusal of a
:math:`\varphi`-dependent ``Symbolic`` beside a sphere, landed with step 8
(:ref:`structured-geometry-specification`).
The MoC solver's isotropic source lift is a hand-written
:math:`1/(4\pi)`, not the angular section (#556). The derivations' typed
measure masses are #557. A composite holding a retraction daggers to the
generic sandwich instead of the leaf pullback (#558). The three
finite-real parsers were #559, made one at L1 by step 7
(:ref:`structured-geometry-one-real-parser`).


.. _structured-geometry-question-values:

The question values: what is asked of a system, with no physics in it
=====================================================================

A reference question (#405) is the materials, the geometry with its
laws, and *what is asked*: a criticality eigenvalue, the flux a source
drives, or the importance of a detector. The first two fields are the
chapters above; this chapter is the third. The module
:mod:`orpheus.numerics.question` holds three question values and two
mode values, and it landed as step 7 of the campaign's first phase on
2026-10-02, with its prerequisite, one definition of a real number at L1
(#559, :ref:`structured-geometry-one-real-parser`):

- :class:`~orpheus.numerics.question.Eigen` ``(parameter, point, mode)``,
  where the system is singular along one direction;
- :class:`~orpheus.numerics.question.FixedSource` ``(source, point)``,
  the flux a given source drives;
- :class:`~orpheus.numerics.question.Response` ``(detector, point)``,
  the importance of a given detector;
- the modes :class:`~orpheus.numerics.question.Fundamental` ``()`` and
  :class:`~orpheus.numerics.question.Nearest` ``(tau)``, which say which
  pole an ``Eigen`` asks for.

The union aliases ``Question = Eigen | FixedSource | Response`` and
``Mode = Fundamental | Nearest`` are the closed sets, and every value is
a content-identity value (:ref:`structured-geometry-content-identity`),
so a cache key built from a question is stable across processes.

**Why this page is the home.** The values are mathematics and live in
``numerics`` (L1), but what they are *for* is the reference question
whose other fields this page defines: the source and the detector a
question holds are the mesh-free functions of the chapter above, and
the question is keyed by the encoder of the chapter before it. The
ontology the values follow is the posing sequence's layer 2,
"the problem (a question bound to a system)"
(``.claude/plans/posing_sequence.md``); the posings the solvers consume
today, :class:`~orpheus.numerics.posing.EigenPosing` and
:class:`~orpheus.numerics.posing.SourcePosing` on the operator pencil,
are documented at :ref:`eigenvalue-posing`. The specification of the
step, with its rulings, its gates and its mutation battery, is
``.claude/plans/reference_p1_spec.md`` §1.7; the census that preceded
it is ``scratch/reference_architecture/p1step7/census.md``, and the two
reviews are ``qa_report.md`` and ``elegance_report.md`` beside it.

.. _structured-geometry-question-values-three:

Three questions over one family of operators
--------------------------------------------

The posing sequence states every question on one object, a family of
operators :math:`E(p)` over a **parameter space** :math:`P`. A point
:math:`p \in P` fixes every coefficient the system declares as variable
(a set of reaction cells scaled together, a geometric extent, a nuclide
density), and :math:`E(p)` is the balance operator of the system at
that point, loss minus production. The **physical point**
:math:`p_\star` is the system as specified. A question is posed at a
**base point** :math:`p_0 = p_\star + \sum_j \delta_j\,e_j`, where each
:math:`e_j` is a declared direction and :math:`\delta_j` its offset; the
offsets are the question's ``point`` field.

**The eigen question** asks where, on the line through :math:`p_0` along
one direction :math:`e_d` (the ``parameter``), the family is singular:

.. math::

   \text{find } \sigma \text{ and } \psi \ne 0 \text{ with }
   E(p_0 + \sigma e_d)\,\psi = 0 .

The solutions :math:`\sigma` are the **poles** of the resolvent
:math:`E(p_0 + \sigma e_d)^{-1}` along that line, and the ``mode`` says
which pole is wanted. When the direction is affine in :math:`\sigma`,
:math:`E(p_0 + \sigma e_d) = A - \sigma T_d` is the operator pencil of
:ref:`the-operator-pencil` restricted to the line. The k-eigenvalue is
the case :math:`e_d` = the fission-emission cells, :math:`T_d = F`:
:math:`A\psi = \sigma F\psi`, and :math:`k = 1/\sigma` is the spectral map
:data:`~orpheus.numerics.posing.K_MAP` reads. The classical
c-eigenvalue (secondaries per collision, the Sood benchmarks) is the
direction that scales every emission cell, scattering and fission
together; a boron search scales one absorption cell; a critical size is
a geometric-extent direction. These are four values of one ``Eigen``,
differing only in the key.

**The fixed-source question** asks for :math:`\psi = E(p_0)^{-1} q`, the
response to a source :math:`q` that does not depend on the unknown. A
subcritical multiplying system driven by :math:`q` is ``FixedSource(q)``
at the physical point; the same source with fission switched off is
``FixedSource(q, point)`` whose point moves the fission-emission
direction to its chart zero (where on the key's chart that zero lies is
the coordinate's chart, which the mode law of #529 reads; the
specification resolves the key and leaves its chart there, by ruling).

**The response question** asks for :math:`\psi^\dagger_R =
E(p_0)^{-\dagger} R`, the importance of a detector :math:`R`: the
adjoint fixed-source question. Its answer serves every source at once,
because the detector's reading of the flux a source drives is

.. math::

   \langle R,\, E(p_0)^{-1} q \rangle
   \;=\; \langle E(p_0)^{-\dagger} R,\, q \rangle
   \;=\; \langle \psi^\dagger_R,\, q \rangle .

.. list-table:: The three questions
   :header-rows: 1
   :widths: 18 30 22 30

   * - Value
     - Asks for
     - Datum
     - Answer
   * - ``Eigen(parameter, point, mode)``
     - the pole :math:`\sigma` of :math:`E(p_0 + \sigma e_d)^{-1}` that
       ``mode`` selects
     - none: the direction is a key
     - the pole and its mode :math:`\psi`; the adjoint mode
       :math:`\psi^\dagger` belongs to this answer, not to the question
   * - ``FixedSource(source, point)``
     - :math:`E(p_0)^{-1} q`
     - the source :math:`q`, a mesh-free function
     - the flux
   * - ``Response(detector, point)``
     - :math:`E(p_0)^{-\dagger} R`
     - the detector :math:`R`, a mesh-free function
     - the importance :math:`\psi^\dagger_R`

Nothing behind the questions exists yet: no system declares a parameter
space, no pencil or spectral map is derived from a parameter, and no mode
law finds a pole. Those are the posing sequence's unit 6 (#529), which
binds a question to a system. Step 7 mints the values, so that the
reference specification (step 8) can hold one and a cache can key on it.

.. _structured-geometry-question-values-physics-free:

Physics-free: the parameter is an opaque key, resolved elsewhere
----------------------------------------------------------------

``Eigen.parameter`` and every key of a point are **opaque keys**:
``numerics`` never reads what a key names. That is the posing sequence's
ruling of 2026-09-27 ("questions are physics-free"), and it is why the
values live in ``numerics``. A question that named a reaction, a
material or an extent would need the vocabulary of the input layer or of
the transport layer, and a module's home is the lowest-knowledge layer
whose vocabulary suffices (:ref:`architecture-layering`). Gate S7.13
holds the layer twice: by AST, every ``orpheus`` import of the module is
under ``orpheus.numerics`` (with the activation leg that its known
import of ``content`` is seen); and in a fresh interpreter, importing the
module and building one value of each kind with a table datum loads only
the ``orpheus`` sub-packages ``numerics`` and ``geometry`` and no SymPy.
``geometry`` is admitted there because a cold ``import orpheus.numerics``
already loads it (``numerics.invariance`` imports
``geometry.transformation``, an exception the layer gate lists); the
positive control is that a cold ``orpheus.sn.problem`` loads
``transport``.

**What a key means is declared by whoever resolves it** (the orchestrator's
ruling 5 of 2026-10-02). A key resolves to a **coordinate declaration**,
one of three kinds:

- ``CellCoefficient``, a *set* of reaction-grid cells scaled together
  (the user's ruling 1 of 2026-10-02). k is the set of fission-emission
  cells; the classical c is every emission cell, scattering and fission;
  a boron search is one absorption cell (the absorption cells arrive with
  #526; today a cell's channel is one of the three emissions,
  :ref:`structured-geometry-specification-coordinates`). The posing
  plan's "one cell's
  coefficient" is the one-element set, which is how the ruling reconciled
  the two meanings c had carried (all collision emission in the
  specification, one cell's coefficient in the posing plan);
- ``GeometryExtent``, a geometric degree of freedom (a critical size);
- ``NuclideDensity(nuclide, regions)``, a number density (deferred until
  the materials keep number densities).

In phase P1 the reference specification resolves every key against its
materials and geometry and refuses one that is not a coordinate of its
problem (:ref:`structured-geometry-specification-canonical`); #529 later
declares the coordinates on the system. Two of the three kinds exist, and
neither lives in ``numerics``: :class:`~orpheus.data.cells.CellCoefficient`
in ``orpheus/data`` and :class:`~orpheus.geometry.extent.GeometryExtent` in
``orpheus/geometry`` (:ref:`structured-geometry-specification-coordinates`).
Numerics therefore names no default parameter: ``Eigen()`` does not
construct (S7.3), because "k" is the system's name for its
fission-emission direction, not a word numerics knows. The specification
derives no default either: its question is required.

**A key must be hashable and have content.** The parameter is parsed by
``_admit_key``: it must be hashable, and the question's digest is then
taken at construction, which refuses a key with no content (a function,
a dataclass compared by identity) with
:class:`~orpheus.numerics.content.ContentlessError` naming
``Eigen.parameter``. A string, an integer, a tuple and a
:class:`~orpheus.numerics.content.FrozenMapping` are admitted (S7.6). The
hashability check was added by the review: on ``31a2dc46`` the parameter
was checked only by taking the digest, and the encoder admitted a
``MappingProxyType`` and a read-only view of a writeable array, so
``Eigen(MappingProxyType(d))`` constructed and, after the caller wrote
``d["cells"] = 2.0``, still equalled the question built from the old
content and kept its cached digest: a wrong cache hit, measured by both
reviews. Both are unhashable, so ``hash(key)`` enforces the ``Hashable``
the annotation already declared; `[M]` this pass, both are refused with
the message "an unhashable mappingproxy is not a key (it is mutable, so
its content is not fixed)" (and "ndarray" for the view). A point key
needs no separate check: the ``dict`` a point is built from refuses an
unhashable key before the question sees it.

``Eigen(None)`` constructs (`[M]`): ``None`` is a hashable,
digestable key. It names no coordinate, and the specification's key
resolution refuses it (S8.1 (h1): "the parameter: None (a NoneType) is
not a coordinate").

.. _structured-geometry-question-values-point:

The point: offsets by key, and why a zero offset is kept
--------------------------------------------------------

``point`` is a :class:`~orpheus.numerics.content.FrozenMapping` from a
parameter key to its **offset from the physical value** (the user's
ruling 2 of 2026-10-02). The empty mapping is the physical point and the
default of all three kinds, and a caller's ``dict`` is frozen at the
boundary: writing the ``dict`` afterwards moves neither the stored point
nor the digest, and the stored point refuses item assignment (S7.8). The
point is the existing frozen mapping, not a ``Point`` class: the name is
a type alias in ``orpheus/derivations/discrete/sn/face_transmission.py``,
and one frozen mapping is the content types' rule. The elegance review
recorded when a point type is owed (a qualified name such as
``ParameterPoint``): when a fourth question holds a point (the evolution
question), or when the point gains a verb. The verb in view is
``moved(key, offset)``: the answer of ``Eigen(p, x, mode)`` is a pole
:math:`\sigma`, the critical point is :math:`x` moved by :math:`\sigma`
along :math:`p`, and the derived question "the source problem at the
critical point" is ``FixedSource(q, x.moved(p, σ))``.

Every offset is parsed by
:func:`~orpheus.numerics.scalars.parse_finite_real`, so the point obeys
the one real-number rule: a NaN or an infinite offset is a
``ValueError`` and a string, a complex or a ``bool`` offset a
``TypeError``, each message naming the key (``the offset of 'boron' is
NaN, which is not a number``); ``-0.0`` is stored as ``+0.0``, and an
offset of ``1`` and one of ``1.0`` are one question. The order in which
keys were given is not content.

**A zero offset is kept, not canonicalised away.** ``Eigen("b", {"b":
0.0})`` and ``Eigen("b")`` are different questions with different digests
(`[M]`, and gate S7.8). A zero offset names a key that the specification
still sees and resolves: a key that is not a coordinate of its problem is
refused there, and dropping it here would hide the refusal. The cost is two cache
keys for one physical question, which is a miss and never a wrong hit,
the only direction a key may err in.

**The base point is not the line's origin.** :math:`\sigma = 0` on the
eigen line is the base point, and the coordinate's *chart zero* (where
its coefficient vanishes) is a different point, declared with the
coordinate. The posing sequence measured the difference: re-posing from
the k pole along c with the parameter's own component removed finds the
c pole, not offset 0.

**The parameter's own key may carry an offset.**
``Eigen("b", point={"b": 0.3})`` is admitted, and the specification does
not refuse a coordinate that is both the parameter and a point key (the
orchestrator's ruling 3 on the step-7 NEEDS). The posing
ontology rules the answer *invariant* under the point's own component
along the question's direction: moving the base point along the line
re-labels the poles, it does not move them. The question is therefore
well posed, a refusal would contradict that law, and its cost is again
a second cache key for one answer. The invariance gate is #529's.

.. _structured-geometry-question-values-modes:

The modes, and what is deferred until its trigger
-------------------------------------------------

``mode`` is a value, ``Fundamental()`` by default (S7.8):

- ``Fundamental()`` asks for the first pole met from the
  removal-dominated end of the direction, on the bulk cone; that end is
  the boundary of the coordinate's declared admissible range. A
  direction whose sign is indefinite (fuel plus moderator) refuses it.
  ``Fundamental`` has no field, so any two are equal.
- ``Nearest(tau)`` asks for the pole nearest :math:`\tau`, read in the
  parameter's chart. ``tau`` is parsed by ``parse_finite_real``: a NaN,
  an infinite or a complex :math:`\tau` is refused, naming ``tau``
  (S7.7).

A string (``"fundamental"``), a bare float and ``None`` are refused as a
mode with a ``TypeError`` (S7.7): a stringly-typed selector is the
missing type the closed set exists to replace. What each mode
*selects* is the mode law of #529, not a property of the value, and no
step-7 gate can see it. An ordinal ``Index(n)`` was refused by the
posing sequence's review: ordering by real part, by modulus and by the
shifted spectrum gives three orders, so "the n-th pole" depends on the
Strategy that found it, and a question may not.

Three members of the ontology are deliberately absent, each with the
event that mints it (rulings 4 and 6). Gate S7.2 asserts the module has
no attribute ``Enclosed``, ``Evolution``, ``Alpha``, ``Time``, ``Index``
or ``Pseudospectrum``, so adding one is an edit of the gate on purpose.

.. list-table:: Deferred, and the trigger that mints each
   :header-rows: 1
   :widths: 24 42 34

   * - Member
     - Why it is not in step 7
     - Trigger
   * - ``Enclosed(region)``, the poles inside a region of the plane
     - no P1 reference asks for more than one mode; beyond the
       direction's continuum edge the request is refused, and the
       information a flagged answer would carry is a question of its own,
       the pseudospectrum over the region
     - the pseudospectrum and spectrum work of #529 and #531
   * - the time direction (the α-eigenvalue, ``Eigen`` along the
       Laplace direction)
     - time is not one of the three declared coordinate kinds; no P1
       reference asks for α. Adding it extends the closed set, it does
       not edit a case
     - the transient question
   * - ``Evolution(initial, source, f, point)``, the initial-value
       question
     - its answer is a function of the generator applied to an initial
       state; nothing in P1 asks it
     - the transient question

.. _structured-geometry-question-values-role:

The role is the type: no adjoint flag
-------------------------------------

No question value carries a forward/adjoint flag, and none has an ``H``,
``adjoint``, ``transpose``, ``dagger`` or ``is_adjoint`` attribute (S7.3).
The role is the question's TYPE (the user's ruling 3 of 2026-10-02):
``FixedSource`` holds a source and ``Response`` a detector, so neither
value has an exclusive-or of fields, and the states a flag would make
spellable cannot be written. `[M]` ``FixedSource(detector=t)`` and
``Response(source=t)`` raise ``TypeError`` (an unexpected keyword
argument), and a datum that is not a mesh-free function is refused with a
message naming the role and the received type (``the source is a
mesh-free function (RegionwiseConstant or Symbolic), got a
numpy.ndarray``; S7.4, eight refusal rows and four positive rows, either
type in either role). The member list in that message is read from the
alias
``MeshFreeFunction = RegionwiseConstant | Symbolic``, which the review
moved into :mod:`orpheus.numerics.mesh_free_function` so that a third
function type is added in one place.

**This is step 6's ruling one level up.** Step 6 ruled that a mesh-free
function carries no role and no density, and that the field holding it is
the role, because the role picks the arrow into phase space: a source
rate enters through the section :math:`E`, a detector through the
retraction's adjoint :math:`R^{\dagger}`, and the two differ by the
measure's mass (:ref:`structured-geometry-mesh-free-two-types`). The
field that picks the arrow is now ``FixedSource.source`` or
``Response.detector``, and the question's type is the role.

**The eigen adjoint belongs to the answer.** ``Eigen`` has no adjoint
either. The adjoint mode :math:`\psi^\dagger` of an eigen answer is the
null vector of :math:`E(\text{pole})^{\dagger}`, the same pole, so
"the adjoint eigen question" asks nothing the eigen question does not;
it is a second component of one answer (posing ruling of 2026-09-28).
The response question is different: its datum is a detector, which the
forward question does not have, so it is a question of its own.

**A source and a detector of one function are two questions.**
``FixedSource(q, p)`` and ``Response(q, p)``, the same function at the
same point, are unequal both ways, have two digests, are two members of
a set, and neither equals its datum (S7.10, over a table and a
``Symbolic``, at the physical point and at an offset point). A cache
therefore never serves a flux for an importance. The schema tag carries
both the class name and the part names (``source`` against
``detector``), so either alone separates the two digests.

.. warning::

   **What no step-7 gate can see: the right function in the wrong
   type.** ``FixedSource(R)`` written where ``Response(R)`` was meant
   constructs, by ruling, because the source and the detector are one
   function type with no units. The digest separates the two (S7.10), so
   the confusion is never a wrong cache hit, but the answer is the wrong
   one. Its catcher is a solver-level value gate owed by phase P4: the
   response :math:`\langle R, A^{-1}q\rangle` against an independent
   reference computed with the roles as intended, on a fixture where
   :math:`A` is not self-adjoint (streaming with vacuum faces,
   anisotropic scattering or upscatter) and :math:`q \ne R` with
   different supports. The adjoint identity
   :math:`\langle R, A^{-1}q\rangle = \langle A^{-\dagger}R, q\rangle`
   cannot see the swap: it holds for either assignment of the two
   functions, so the swap lies inside its invariance group (V&V failure
   mode 12), and on a self-adjoint fixture the swapped response is even
   equal.

The elegance review raised, and did not grade, the name ``Response``:
what the question asks for is the importance :math:`\psi^\dagger_R`,
and the response :math:`\langle R, \psi\rangle` is a scalar derived from
it. The name was ruled (ruling 3) and is kept.

.. _structured-geometry-question-values-identity:

Content identity, admitted eagerly
----------------------------------

Every value is a frozen dataclass declared ``eq=False`` that takes
``==`` and ``hash`` from :class:`~orpheus.numerics.content.ContentIdentity`
(:ref:`structured-geometry-content-identity-mixin`). The parts are the
fields: ``Eigen`` (``parameter``, ``point``, ``mode``), ``FixedSource``
(``source``, ``point``), ``Response`` (``detector``, ``point``),
``Nearest`` (``tau``) and ``Fundamental`` (none).

**Admission is eager.** Each ``__post_init__`` parses its fields and then
takes the value's digest, so a question that cannot be keyed raises at
construction (``ContentlessError``, ``ValueError`` or ``TypeError``,
naming the field), never first when a cache calls ``hash``. A value that
exists can always be keyed. Gate S7.5 builds a point whose key has no
content (a function; a dataclass compared by identity) for each kind and
requires the refusal at construction; the battery's arm that made
admission lazy reddened those rows.

The content rows are the encoder's pattern of step 5, applied to five
new types, which join the content roster that gate S5.9 closes (S5.9
walks ``orpheus.numerics`` with ``pkgutil``, and `[M]` on the
specification's prototype it reddened with "content types with no
roster entry: ['Eigen', 'FixedSource', 'Fundamental', 'Nearest',
'Response']" until the union was amended). Per type, the population of
parts is read off ``content_parts`` and matched to the roster; each part
moved by one leg (one ULP of one offset, a key added or relabelled, the
mode and :math:`\tau`, one coefficient of the datum) moves the digest,
``==`` and set membership while no other part moves (S7.9, 22 legs);
equal content is one value (the point's insertion order, a ``dict``
against the ``FrozenMapping``, ``-0.0`` against ``0.0``, ``1`` against
``1.0``); pickle round trips. S7.11 prints six values' digests and
hashes in two interpreters under ``PYTHONHASHSEED`` 1 and 2, with a
``str`` hash as the control that must differ. S7.12 is the RECORD: the
hexadecimal digests of four canonical questions (``Eigen("fission-emission")``;
an ``Eigen`` with a tuple key, a two-key point and ``Nearest(0.5)``;
``FixedSource(t, {"fission-emission": -1.0})``; ``Response(t)``, with
``t`` a 2 × 2 table) pinned as literals. A moved digest invalidates every
cache key that holds a question, the message says so, and a re-pin
carries its reason; qa's arm that bumped ``Eigen.__content_version__``
reddened it alone.

.. _structured-geometry-question-values-refuted:

The rulings, and the specification they replaced
------------------------------------------------

The first specification of the step (2026-09-25; in git at
``403f357d``, ``.claude/plans/reference_p1_spec.md`` §1.7) was a closed
sum of four cases, ``Eigen(k)``, ``Eigen(c)``, ``FixedSource(σ)`` and
``CriticalParameter``, each with a forward and an adjoint, and no
``Eigen(α)``. The step's census found it stale on every case: it
predated the posing sequence's rulings of 2026-09-27 to 2026-09-29.
`[REFUTED 2026-10-02]` for the question *"what are the values a question
is spelled with?"*; each part was retired for a structural reason.

.. list-table:: The 2026-09-25 cases, and why each was retired
   :header-rows: 1
   :widths: 24 46 30

   * - Retired
     - Structural reason
     - Replaced by
   * - ``Eigen(k)`` and ``Eigen(c)`` as two cases
     - a case per eigenvalue name enumerates directions inside the type,
       so a boron search or a critical size is a new case and an edit of
       every match; the direction is data. Its spectral map was to be a
       field (``Eigen(k)``'s map *is* ``K_MAP``), but a
       ``SpectralMap`` holds three lambdas, which have no content, so the
       old gate could only digest the map's ``repr``. And "c" meant two
       things in two plans
     - one ``Eigen(parameter, point, mode)``; the spectral map is derived
       from the parameter by #529, not stored; c is a set of cells
       (ruling 1)
   * - ``FixedSource(σ)``, σ = 1 multiplying and σ = 0 fission off
     - σ is a coordinate of the base point along one direction, so a
       field for it welds one direction into the fixed-source type and
       leaves boron or an extent unspellable
     - the point (ruling 2): the multiplying question is the physical
       point, fission off is the point that moves the fission-emission
       key to its chart zero
   * - ``CriticalParameter``
     - a critical size is a pole along a geometric direction, which is an
       ``Eigen`` whose key is a ``GeometryExtent``; a separate type is
       ``Eigen`` with the direction welded in. The posing ontology says
       it "does not exist in the tree, and nothing is dissolved"
     - ``Eigen`` over a ``GeometryExtent`` key, resolved by step 8
       (S8.1 (h3) refuses the extent without a geometry). S7.2's AST
       census of class definitions under ``orpheus/`` (more than 100
       files, positive control ``SpectralMap`` found) reds the day the
       name is re-minted
   * - a forward and an adjoint per case (a boolean flag, or an ``.H``)
     - a flag makes the detector-less adjoint and the forward question
       given a detector spellable, and both must then be refused; for an
       eigen question the adjoint is part of the answer, and for a fixed
       source the adjoint has a different datum
     - ``Response(detector, point)``, its own type (ruling 3); no flag on
       any value
   * - "no ``Eigen(α)``", with ``ALPHA_MAP`` refused
     - α is the time direction, which is absent until its trigger, not
       forbidden; refusing a spectral map at the question would have
       placed the map in the question again
     - deferred to the transient question (ruling 4)
   * - the specification's ``source`` field beside the question
     - the question now holds its datum, so a source on the
       specification as well would be a second definition of one datum,
       and every refusal that kept the two consistent would guard a state
       the design should make unspellable
     - no ``source`` field: the specification's fields are a material and
       a question, or the materials, a geometry and a question
       (:ref:`structured-geometry-specification-layer`), and the old
       refusals S8.1 (d) to (f) are structural legs (the orchestrator's
       ruling 1 on the step-7 NEEDS)

**The rulings of record**, all of 2026-10-02:

1. The user: a parameter's direction is a set of reaction-grid cells
   scaled together.
2. The user: the point is a frozen mapping from parameter key to an
   offset from the physical value; the empty mapping is the physical
   point and the default.
3. The user: the adjoint fixed-source question is its own question,
   ``Response(detector, point)``; no value carries a forward/adjoint
   flag, and the eigen adjoint belongs to the eigen answer.
4. The orchestrator: α is not in step 7.
5. The orchestrator: the parameter is an opaque key in ``numerics``; the
   coordinate kinds live with whoever resolves the key.
6. The orchestrator: the modes are ``Fundamental`` and ``Nearest(τ)``.

And on the specification's NEEDS: the specification's ``source`` field
retires; the coordinate declaration's shape is step 8's first design
item; ``Eigen("b", point={"b": 0.3})`` is admitted.

**Found by the reviews, and fixed in** ``8122abc6``: the parameter's
missing hashability check (above); the frozen encoder's admission of a
``MappingProxyType`` part, now refused
(:ref:`structured-geometry-one-real-parser` has the array analogue,
#561); the four spellings of "finite" and the two of "a real entry"
(the next section); ``MeshFreeFunction`` defined in its consumer, moved
to its module. Withdrawn by the elegance review on its second pass, each
with its reason: a shared base for ``FixedSource`` and ``Response`` (the
field name IS the role, and a base would read the field by a string); a
generic admit-by-annotation hook on the mixin (three instances whose
parses differ); ``Eigen(True) == Eigen(1)`` as a defect (it follows
Python's ``==`` in the ``dict`` and in the digest alike).

.. _structured-geometry-question-values-gates:

The gates, and the battery
--------------------------

Every gate is a ``foundation`` test declaring ``rests_on``. `[M]`
2026-10-02, this pass, ``.venv/bin/python -O -m pytest`` over the two
files: 116 rows, 116 passed.

.. list-table::
   :header-rows: 1
   :widths: 9 55 36

   * - Id
     - What it asserts
     - Mutation witness (the battery's arm)
   * - S7.1
     - the closed set: ``get_args(Question)`` and ``get_args(Mode)`` are
       exactly the three and the two; the module's content classes are
       exactly those five; a ``match`` with ``assert_never`` names each
       once, and a foreign object reaches ``case _``
     - an ``Enclosed`` class minted in the module
   * - S7.2
     - the struck and deferred names are absent; no class
       ``CriticalParameter`` anywhere under ``orpheus/`` (AST census, input
       count printed, positive control)
     - the same arm
   * - S7.3
     - the fields are exactly the roles, by name; no field is a
       ``bool``; no adjoint attribute; the wrong-role spellings and
       ``Eigen()`` do not construct; the signatures' defaults; two
       detectors refused
     - a boolean ``adjoint`` field on ``FixedSource``; a default
       ``parameter``
   * - S7.4
     - a datum that is not a mesh-free function is refused, naming the
       role; both types admitted in both roles
     - the datum check removed (8 of 8 refusal rows red)
   * - S7.5
     - a non-finite or non-real offset refused, naming the key; a
       contentless point key refused at construction
     - the finiteness check removed; admission made lazy
   * - S7.6
     - a contentless or unhashable parameter refused (a function, an
       identity-compared dataclass, a ``list``, a ``MappingProxyType``, a
       read-only view); a ``str``, ``int``, ``tuple`` and
       ``FrozenMapping`` admitted
     - admission made lazy
   * - S7.7
     - :math:`\tau` finite and real; a string, a float and ``None``
       refused as a mode
     - the finiteness check removed; the mode check removed
   * - S7.8
     - the default point is an empty ``FrozenMapping`` and the default
       mode ``Fundamental()``; ``point={}`` and ``point=FrozenMapping()``
       are one value; a caller's ``dict`` is frozen at the boundary; a
       zero offset is not the empty point
     - a non-empty default point; the point stored as a read-only view of
       the caller's ``dict``
   * - S7.9
     - every part is content (population, 22 perturbation legs, equal
       pairs, pickle)
     - the positive control (below)
   * - S7.10
     - a source and a detector of one function are two questions
     - the schema tag reduced to a constant (the digest leg; the ``==``
       leg is a second, independent tooth)
   * - S7.11
     - the same digests and hashes in two processes
     - a ``str`` part encoded through the salted ``hash()``
   * - S7.12
     - RECORD: four digests pinned
     - qa's ``__content_version__`` bump
   * - S7.13
     - the layer, by subprocess and by AST
     - the subprocess also importing ``orpheus.transport``; the module
       importing ``orpheus.mesh``

The specification's battery (``scratch/reference_architecture/p1step7/ta/battery/``,
a ``-p`` plugin that rebinds in every ``sys.modules`` binding at each
test's setup and refuses an arm that rebinds nothing) ran 14 arms and a
positive control, first against a prototype and then against the real
module at ``31a2dc46`` (114 rows then; logs ``real/arm_*.log``). `[M]`
the unmutated run 0 red; the positive control, every content part
emptied, 33 red; each arm reddened its target rows (red counts 1 to 17),
and the red SETS, read against the target rows, are the record, since
the counts moved between the prototype and the module (A3 9 rows, not 12:
the encoder's NaN refusal now goes through the one parser and names the
key). The last leg of the specification's S7.13 puts
``"orpheus.numerics.question"`` (with ``mesh_free_function``, owed by
S6.21, and ``scalars``) on the entry-point list of
``tests/gates/test_layer_imports.py``, so each imports cleanly from a cold
interpreter.

**Landed since, and not built.** The coordinates and the key resolution,
with the refusal of ``Eigen(None)`` and of a key that is not a coordinate
of the problem, landed with step 8 (S8.1 (h1) to (h4);
:ref:`structured-geometry-specification`). Not built: the pencil
and spectral map a resolved parameter derives, the mode law, the
binding to a system, and the invariance gate under the point's own
component: #529. ``Enclosed`` and the pseudospectrum question: #529 and
#531. The time direction and the evolution question: the transient
question. The role-confusion value gate: P4.


.. _structured-geometry-specification:

The reference specification: a question with its materials, keyed
=================================================================

A reference (#405) answers one question about one system, and the
**reference specification** is that question written down: the materials,
the geometry with its boundary laws when the problem has one, and what is
asked (:ref:`structured-geometry-question-values`). It is never an answer.
It is the key a reference cache stores an answer under, so it is a
content value (:ref:`structured-geometry-content-identity`), admitted in a
canonical form at construction: two specifications that ask the same
question of the same system are one value with one digest. It landed as
step 8 of the campaign's first phase on 2026-10-02, in a new package,
:mod:`orpheus.specification`, with the two coordinates its questions are
keyed by, :class:`~orpheus.data.cells.CellCoefficient` and
:class:`~orpheus.geometry.extent.GeometryExtent`.

.. code-block:: python

   from orpheus.data.cells import CellCoefficient, Channel
   from orpheus.data.materials import Materials
   from orpheus.derivations.continuous.sood_registry import SOOD2003_CASES
   from orpheus.geometry import GeometryExtent
   from orpheus.numerics.question import Eigen
   from orpheus.specification import GeometrySpecification, InfiniteMediumSpecification

   k = Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))

   # k-infinity: one material, posed on energy alone (there is no geometry field)
   infinite = SOOD2003_CASES["PUa-1-0-IN"]
   (material_id, mixture), = infinite.materials.items()
   k_inf = InfiniteMediumSpecification(material_id, mixture, k)

   # a critical size: the width of interval 0 of a bare sphere
   sphere = SOOD2003_CASES["Ua-1-0-SP"]
   critical_radius = GeometrySpecification(
       Materials(sphere.materials), sphere.to_geometry(), Eigen(GeometryExtent(0)),
   )

   k_inf.question.parameter   # resolved: the one cell (0, FISSION_EMISSION), no quantifier
   k_inf.question == k        # False: the stored question is the canonical one
   critical_radius.unreadable # (Symbolic.phi,): a sphere has no azimuth reference

**Why this page is the home.** A specification composes the three values
the chapters above define: the geometry with its laws, the content
encoder, the mesh-free functions and the question values. It is the
consumer for which each of them was built. The specification of the step,
with its rulings, its gates and its mutation battery, is
``.claude/plans/reference_p1_spec.md`` §1.8 (the blocks "The two types"
and "After the elegance re-review" record the final shape) and the
sections "The step-8 rulings", "Step 8 built" and "The infinite medium is
the point in phase space" of the same file; the census that preceded it
is ``scratch/reference_architecture/p1step8/census.md``, and the four
reviews are ``qa_report.md``, ``qa_report_2.md``, ``elegance_report.md``
and ``elegance_report_2.md`` beside it.

.. _structured-geometry-specification-layer:

The layer a question is posed at is the type
--------------------------------------------

The posing filtration commits phase space in stages: the materials, then
the geometry with its laws, then the mesh, then the method
(:doc:`/architecture/conceptual_view`). A reference question is posed at
one of two of those stages, and each stage is a type:

.. list-table:: The two specification types
   :header-rows: 1
   :widths: 24 38 38

   * -
     - :class:`~orpheus.specification.specification.InfiniteMediumSpecification`
     - :class:`~orpheus.specification.specification.GeometrySpecification`
   * - Fields
     - ``material_id``, ``mixture``, ``question``
     - ``materials``, ``geometry``, ``question``
   * - Posed on
     - energy alone: the infinite medium, the point in position and in
       direction
     - a finite :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
       with its boundary laws
   * - ``materials``
     - ``Materials({material_id: mixture})``, built on read
     - the given ``Materials``, restricted to the ids the geometry assigns
   * - ``n_groups``
     - the mixture's ``ng``
     - :meth:`Materials.uniform_group_count
       <orpheus.data.materials.Materials.uniform_group_count>` over the
       materials kept, read once at construction
   * - ``n_regions``
     - 1 (a class constant)
     - ``len(geometry.intervals)``, one region per interval
   * - ``unreadable``
     - ``(r, mu, phi)``: no position and no direction chart
     - ``(phi,)`` when the chart has no azimuth reference (the sphere),
       otherwise ``()``
   * - ``resolve_extent(e)``
     - refuses: "names an extent, and the infinite medium has no geometry"
     - :meth:`GeometryExtent.resolve
       <orpheus.geometry.extent.GeometryExtent.resolve>` against the
       geometry

``Specification = InfiniteMediumSpecification | GeometrySpecification``
is the closed union, exported with ``Coordinate = CellCoefficient |
GeometryExtent``. Each type answers the five members of that surface for
itself, and the admission (the canonical question, the datum's fit, the
digest) is written once, as free functions over the union that read only
the surface. No code path asks which layer it is on: `[M]` this pass, ``git
grep`` finds 0 ``isinstance`` tests on either type in ``orpheus/`` and
``tests/`` (the pattern's control, ``isinstance(`` followed by
``StructuredGeometry``, finds the one known site in the module). Gate
S8.1's structural row ``test_s8_1_the_layer_is_the_type`` asserts the two
field tuples, the union and the infinite medium's ``materials``,
``n_regions`` and ``unreadable``.

**Why the infinite medium has no geometry field.** The infinite medium is
posed on energy alone, and it holds a material and nothing else: admitting
a geometry would admit spatial dimension, and with it more than one
material and a direction chart, none of which it has. The physics, and
why the definition is exact for every question it admits, is
:ref:`infinite-medium-definition` (the user's ruling of 2026-10-02). In
the type this is two rules made unspellable rather than checked: "the
infinite medium has one material" is a type with one ``mixture`` field,
and "the infinite medium has no geometry" is a type with no ``geometry``
field. ``GeometrySpecification`` refuses ``geometry=None`` with a
``TypeError`` ("the geometry is a StructuredGeometry, got a NoneType").

.. dropdown:: First got wrong: the infinite medium as ``geometry=None``, and then as a geometry value
   :color: muted

   **What was tried.** The first build of the step (``29fe4266``) was one
   type, ``Specification(materials, geometry | None, question)``, with
   ``None`` standing for the infinite medium. Its elegance review measured
   the cost: ``geometry is None`` was tested 5 times in 3 functions (the
   materials rule, the extent resolution, the region count and the chart
   read). Its finding F5 named the pattern, a repeated conditional is a
   missing type, and proposed the destination an earlier review had
   proposed too (finding J of the W5 review of 2026-09-25, "the infinite
   medium is a geometry value, not ``None``"): a geometry value standing
   for the infinite medium, with one region, no extent and no chart, from
   which the 5 branches would fall out. That value was written in the
   working tree and never committed.

   **Why it failed.** A geometry with no placement and no finiteness is
   not a geometry: it defines nothing a geometry defines, so it carries no
   mathematical information, and it misdirects, because whatever accepts
   a geometry admits spatial dimension (the user, 2026-10-02). The
   finding's FACT stands, and it is what the design rests on: repeated
   ``is None`` discrimination is a missing type. The missing type is the
   layer, not a degenerate geometry.

   **What replaced it.** Two types, one per layer (``eaa74163``). `[M]`
   the second elegance review found 0 ``geometry is None`` tests left; the
   only ``is None`` in the module reads the chart's azimuth reference.

.. _structured-geometry-specification-canonical:

The canonical form: one question, one key
-----------------------------------------

Every specification is stored in a canonical form, computed at
construction, so that the digest, ``==`` and ``hash`` see the question and
not its spelling.

- **A spectator is not part of the key** (the orchestrator's ruling 3 on
  the step-8 NEEDS, the leak principle). ``GeometrySpecification`` keeps
  only the materials its geometry assigns,
  ``materials.restrict(set(geometry.mat_ids))``: a declared material that
  no interval uses changes no answer, so carrying a channel or not, of the
  same group count or not, it changes neither the specification nor its
  digest. ``restrict`` also refuses an assigned id the declaration lacks,
  with the pinned fragment ``references material ids [1]``.
- **Every key is resolved.** The ``Eigen`` parameter and every key of the
  question's point must be a coordinate of this problem, and each is
  replaced by its resolved form: a
  :class:`~orpheus.data.cells.CellCoefficient` by the explicit non-zero
  cells it scales in the kept materials, a
  :class:`~orpheus.geometry.extent.GeometryExtent` by itself once the
  geometry is shown to have its interval. A key that is not a coordinate
  is a ``TypeError`` naming where it sits ("the parameter: 'fission-emission'
  (a str) is not a coordinate"; the same for "the point key"), so
  ``Eigen(None)``, which numerics admits as a hashable key, is refused
  here (S8.1 (h1)).
- **Two spellings of one direction are one specification.**
  ``every(FISSION_EMISSION)`` and the explicit set of the fissile
  materials' fission cells resolve to one ``CellCoefficient``; an explicit
  zero cell beside a non-zero one is dropped, in the parameter and in a
  point key alike (qa's finding F4 found the point-key case unpinned; the
  row ``test_s8_10_a_zero_cell_in_a_point_key_is_dropped`` now pins it).
  Two point keys that resolve to one coordinate are refused, naming both
  (S8.1 (h5)).
- **The question is required.** A specification derives no default
  question: a physics default ("k") would be an optimistic default, and
  on a problem with no fissile material it would be refused anyway.

``spec.question`` is therefore the CANONICAL question, which may be
unequal to the value the caller passed (`[M]` ``k_inf.question == k`` is
``False`` in the example above). It follows that a specification is
re-posed from the caller's question, never from ``spec.question``. qa's
finding F2 measured why: ``dataclasses.replace(spec, materials=other)``
re-runs the admission on the stored, already-resolved question, so a
question written ``every(FISSION_EMISSION)`` would then name only the
materials that were fissile before, a different question posed silently
(a cache miss, never a wrong hit, since the digests differ). Storing the
caller's question beside the canonical one was considered and refuted by
the second elegance review: two specifications whose canonical forms are
equal would then compare equal and behave differently under ``replace``,
so equality would stop being a congruence. With the canonical form stored
alone, equal values behave alike, and the module docstring states the
re-posing rule.

**The order of admission.** ``GeometrySpecification.__post_init__`` checks
the field types, restricts the materials, reads the group count, then
admits the question (its type, its datum, its keys) and takes the digest;
``InfiniteMediumSpecification`` parses its id and its mixture first. Two
wiring rows of S8.1 pin the order where it decides a message: the group
count is read before the keys, and the composition before the keys.

.. _structured-geometry-specification-coordinates:

The coordinates: a set of cells, and the width of an interval
-------------------------------------------------------------

The question values hold opaque keys
(:ref:`structured-geometry-question-values-physics-free`); a specification
resolves them, and a key IS the coordinate value: there is no table from
names to coordinates, so nothing is defined twice and nothing is declared
implicitly (the user's ruling 1 of 2026-10-02). Two kinds of coordinate
exist, each in the lowest layer whose vocabulary defines it.

**The channels and the cell coefficient** (:mod:`orpheus.data.cells`). A
material's cross sections form a grid of cells; a cell is the pair
``(material id, Channel)``, and :class:`~orpheus.data.cells.Channel` is a
closed enum of the three **emission** channels (the user's ruling 2):

.. list-table:: The channels, and what carrying one means
   :header-rows: 1
   :widths: 26 30 44

   * - ``Channel``
     - The operator it scales
     - A mixture carries the cell when
   * - ``FISSION_EMISSION``
     - the fission emission :math:`\chi \otimes \nu\Sigma_f`
     - it is producing (:attr:`Mixture.is_producing
       <orpheus.data.macro_xs.mixture.Mixture.is_producing>`,
       :math:`\nu\Sigma_f > 0`); a producing mixture's :math:`\chi` is a
       probability simplex, which ``Mixture`` enforces, so the emission is
       non-zero
   * - ``SCATTERING_EMISSION``
     - the scattering transfer :math:`\Sigma_s`
     - some block of its Legendre stack ``SigS`` has a non-zero entry (a
       stack whose only non-zero block is :math:`P_1` is carried)
   * - ``N2N_EMISSION``
     - the :math:`(n,2n)` transfer :math:`\Sigma_{2n}`
     - some block of ``Sig2`` has a non-zero entry

"Carries" means **the cell is non-zero** (the orchestrator's ruling 5):
every channel field exists on every ``Mixture``, so a predicate on the
field's existence would refuse nothing. The predicate is
:meth:`Channel.is_carried_by <orpheus.data.cells.Channel.is_carried_by>`,
one exhaustive ``match`` with an arm per member, so a fourth member
without an arm reddens pyright through the declared ``-> bool`` (`[M]` the
second elegance review's mutant ``m2``).

**Why only the emission channels.** Every method's fission operator is the
fission emission alone, and the scattering and :math:`(n,2n)` emissions
are the other gains, so scaling an emission cell moves exactly the
operator it names. A removal channel (a capture cell, an absorber search)
is not a cell yet, because :math:`\Sigma_t` is stored on the ``Mixture``
beside its parts: scaling a capture cell would not move the collision
operator, which reads the stored total (the census, §1; `[M]` at
``11a3f058``). The removal cells arrive with the reaction grid of posing
unit 3 (#526), which derives every total where it is used.

A :class:`~orpheus.data.cells.CellCoefficient` is a DIRECTION in the
system's parameter space, the set of cells scaled together. It has **two
fields**: ``cells``, the explicit pairs, and
``channels_in_every_material``, the channels named in every material that
carries them, spelled :meth:`CellCoefficient.every
<orpheus.data.cells.CellCoefficient.every>`. Both are frozen sets, so the
order and repetition of the input are not content, and a direction names
at least one cell or channel. k is ``every(FISSION_EMISSION)``; the
classical c (secondaries per collision) is every emission cell; a single
cell is a one-element set.
:meth:`CellCoefficient.resolve <orpheus.data.cells.CellCoefficient.resolve>`
turns a direction into the explicit non-zero cells it scales in a given
``Materials``: a channel named in every material becomes the cells of
every material that carries it, an explicit cell the material does not
carry is dropped (it scales a zero), a cell on a material outside the
problem is refused ("material 9 is not among the problem's materials
(ids: [0, 1])"), and a direction with no non-zero cell left is refused as
a "zero direction" (a zero direction has no pole to find). The resolved
key has an empty second field, so "a cache key holds no quantifier" is a
property of the stored value that a gate asserts, not a promise of the
resolving code (the first elegance review's finding F7 replaced a
quantifier token stored inside the cell set, ``EVERY_MATERIAL``, with the
second field). Resolving a resolved key returns it.

**The geometric extent** (:mod:`orpheus.geometry.extent`).
:class:`~orpheus.geometry.extent.GeometryExtent` ``(interval)`` is the
width of one interval of a ``StructuredGeometry``, in cm, with every
interval outside it translated outward (the user's ruling 3). On a
one-interval body it is the width of that interval: the critical radius
of a solid sphere, or the full width of a bare slab in the Sood
benchmarks. On a reflected body it is the core grown under a reflector of
fixed thickness, which is the coordinate the Neshat–Maiorino reflected
slab varies. It covers the 31 critical-extent references of the Sood and
Atalay registries (`[M]` the census at ``11a3f058``, of 53 cases: 12
one-group slabs, 6 Sood and 6 Atalay, 5 two-group slabs, 11 spheres and 3
cylinders; the other 22 ask :math:`k_\infty` with no geometry). The index counts from the innermost
interval and is non-negative: ``-1`` is refused rather than read as "the
last interval", which would be a second spelling of one coordinate whose
meaning moves with the interval count. :meth:`GeometryExtent.resolve
<orpheus.geometry.extent.GeometryExtent.resolve>` refuses an index the
geometry does not have.

**The homes.** ``Channel`` and ``CellCoefficient`` live in
``orpheus/data``, because they name a ``Mixture``'s channels;
``GeometryExtent`` lives in ``orpheus/geometry``. A single protocol
``resolve(spec)`` for both was rejected by the first elegance review for
the layering: one signature would have to take both the materials and the
geometry, or a specification, so ``data`` would read geometry or both
packages would import the layer above them. The one ``match`` over
``Coordinate`` sits in the specification, the one place that knows both,
and is the boundary dispatch over a closed set.

.. _structured-geometry-specification-datum:

The datum fits the problem
--------------------------

A ``FixedSource`` holds a source and a ``Response`` a detector, each a
mesh-free function (:ref:`structured-geometry-mesh-free-functions`). The
specification admits the datum against its problem, with one ``match``
over the two function types that ends in ``assert_never`` (the second
elegance review's finding G1: a third function type would otherwise be
admitted unchecked):

- **groups**: the datum's group count is ``spec.n_groups``;
- **regions**: a ``RegionwiseConstant`` has ``spec.n_regions`` regions,
  one per interval of the geometry. The count is the intervals', not the
  distinct materials': on a slab whose intervals hold the materials
  ``(1, 0, 1)`` there are 3 regions. qa's finding F1 found every region
  and extent fixture assigning one material per interval, where the two
  counts agree, so a mutant counting distinct materials was green on 165
  of 165 rows; the fixture ``slab3_repeated()`` now separates them;
- **coordinates**: a ``Symbolic`` depends on no coordinate in
  ``spec.unreadable``. The infinite medium reads none of :math:`r`,
  :math:`\mu`, :math:`\varphi`; a geometry reads all three, except the
  azimuth :math:`\varphi` on a chart with no azimuth reference. That is
  the sphere, where a reference perpendicular to :math:`\hat e_r` at every
  position would be a continuous nowhere-zero tangent field on the
  sphere, which the hairy-ball theorem forbids
  (:ref:`structured-geometry-angular-chart`). The rule reads the chart's
  ``azimuth_reference``, not the coordinate-system tag (the first
  elegance review's finding F2: the chart is where the fact is declared,
  and this was its first consumer).

The dependence test is :meth:`Symbolic.depends_on
<orpheus.numerics.mesh_free_function.Symbolic.depends_on>`, the one
predicate the isotropy test also reads (``is_isotropic`` is ``not
depends_on(mu, phi)``). A function does NOT depend on the named
coordinates iff, in every group, ``simplify`` reduces
:math:`q_g - q_g|_{c \to c'}` to 0, every named :math:`c` replaced by a
fresh real symbol at once. An undecided difference counts as a
dependence, the side a consumer refusing the dependence refuses; so
``sin(φ)**2 + cos(φ)**2`` is independent of :math:`\varphi` and admitted
beside a sphere, which a free-symbols test would get wrong. The query
refuses an empty argument list and a symbol the function does not own
(the first elegance review's finding F4: a caller's ``Symbol("phi")``
without ``real=True`` answered "independent", the optimistic answer).

.. _structured-geometry-specification-eager:

Admission is eager, and what a key cannot hold
----------------------------------------------

Each type takes its digest at the end of ``__post_init__``, so a
specification that cannot be keyed is refused at construction, never
first when a cache calls ``hash``. The encoder is the one definition of
"keyable", so a separate parse would be its twin (the first elegance
review withdrew its own objection to digest-as-admission on that ground,
and the digest is memoised, so the eager call costs nothing extra). What
the earlier parses leave for the digest to refuse is the boundary-law
payload: a geometry whose law holds a function, such as
:class:`~orpheus.geometry.boundary.prescribed_inflow.PrescribedInflow`
built from a callable, is refused with
:class:`~orpheus.numerics.content.ContentlessError` naming
``GeometrySpecification.geometry`` (S8.1 ``contentless_geometry``).

.. warning::

   **An obligation for phase P4, not met in P1.** The non-vacuum
   manufactured-solution reference
   (``orpheus/derivations/continuous/mms/sn.py``, its prescribed-inflow
   branch) declares a callable-bearing inflow law, so a specification of
   it cannot be built today. The census counts 13 manufactured-solution
   case types among the fixed-source producers; before P4 stores their
   references, the inflow needs a content-bearing law, a ``Symbolic``
   boundary inflow on the trace. How many of the 13 carry a non-vacuum
   face is not counted (`[R]`).

.. _structured-geometry-specification-group-count:

The group-count rule has one home
---------------------------------

A problem is one energy discretisation, so its materials must agree on the
group count. The rule lives once, in the input layer:
:meth:`Materials.uniform_group_count
<orpheus.data.materials.Materials.uniform_group_count>` returns the one
count of the declared mixtures or raises
:exc:`~orpheus.data.materials.InconsistentMaterialsError` (a
``ValueError``, message fragment "uniform ng"). The exception moved from
``orpheus/transport/mesh/material_mesh.py`` to
:mod:`orpheus.data.materials` with every importer re-pointed and no shim;
``MaterialMesh.ng`` is the call.

The two callers read two different sets, by ruling. ``MaterialMesh``
reads every declared material, its rule since before step 8, now stated:
a material mesh is built from a declaration and refuses a declaration of
mixed group counts. The specification reads the materials its geometry
assigns, because it restricts first (ruling 3), so a spectator of another
group count is dropped, not refused (the row
``test_s8_1_a_spectator_of_another_group_count_is_dropped``). The
infinite medium reads its one mixture's ``ng`` and never consults the
rule. Gate S8.9 holds the home: the class is defined in
``orpheus.data.materials`` and stays a ``ValueError``; an AST census of
``orpheus/`` and ``tests/`` finds every importer reading that home; and a
route row installs a decoy rule raising a sentinel and requires both
``MaterialMesh`` and the specification to raise it.

.. _structured-geometry-specification-deferred:

What is not here, and the event that brings each
------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - Absent
     - Why
     - Trigger
   * - ``NuclideDensity(nuclide, regions)``, the third coordinate kind
     - ``Materials`` holds macroscopic cross sections only: number
       densities are summed away in ``compute_macro_xs``, so no fixture
       can declare a nuclide, and no P1 or P4 reference asks a
       composition search (the census: 0 boron references)
     - ``Materials`` keeping number densities (the user's ruling 4)
   * - the removal cells (capture, absorption)
     - :math:`\Sigma_t` is stored beside its parts, so scaling a part
       would not move the collision operator
     - posing unit 3 (#526), which un-stores the totals
   * - a coordinate's chart: its zero, its physical value, its admissible
       range
     - read only by the mode law, which finds a pole
     - #529 (the orchestrator's ruling 6)
   * - binding a specification to a method's problem
     - P1 mints the key; the test of S8.4 lifts a specification by hand
       (``MaterialMesh(Mesher(spec.geometry)...mesh, spec.materials)``)
     - phase P4
   * - the direction-dependent problem uniform in space
     - it is the spatial marginal of a body problem, the retraction of a
       geometry problem's solution along space
       (:ref:`infinite-medium-spatial-marginal`)
     - none: it needs no specification type

Two consequences are accepted, each a cache miss and never a wrong hit:
the infinite medium of one mixture under two material ids is two keys,
because the id is inside the resolved cells and ``CellCoefficient`` names
ids by ruling; and an infinite medium and a reflective slab of the same
mixture ask one physical question and are two keys (the row
``test_s8_10_the_two_layers_are_two_keys``). The two extent refusals of
the infinite medium, as the parameter and as a point key, print one text,
accepted by the orchestrator: it is the coordinate's resolution refusal,
worded like ``CellCoefficient.resolve``'s and ``GeometryExtent.resolve``'s,
neither of which names where the key sits.

.. _structured-geometry-specification-gates:

The gates, and the battery
--------------------------

Every gate is a ``foundation`` test except the assembly rung (``l1``).
`[M]` 2026-10-02, this pass, at ``65b533fe``,
``.venv/bin/python -O -m pytest`` over the seven files: 204 rows, 204
passed (``tests/gates/specification/`` 142, ``tests/gates/data/test_cells.py``
40, ``tests/gates/geometry/test_geometry_extent.py`` 8,
``tests/gates/numerics/test_mesh_free_depends_on.py`` 14).

.. list-table::
   :header-rows: 1
   :widths: 8 56 36

   * - Id
     - What it asserts
     - Mutation witness (the battery's arm)
   * - S8.1
     - the refusals, each keyed on a message fragment the triggering
       argument determines, and the fragments pairwise disjoint; the
       admission of each rule's legal inputs; the wiring order; the layer
       is the type; the azimuth rule per coordinate system; the 3 × 3
       channel-carriage matrix through the specification; the fixture
       mixtures' balance
     - one swallowing arm per refusal (``Xg_regions``, ``Xg_groups``,
       ``Xg_coordinates``, ``Xh12``, ``Xh3``, ``Xh5``, ``Xcontent``);
       ``B1``, ``C1``; ``N1`` (regions by distinct materials); ``U1`` to
       ``U3`` (the readable coordinates); ``R1``; ``F2`` (a ``geometry``
       field on the infinite medium); ``SPECTATOR`` (spectators kept)
   * - S8.2
     - digest: a roster of the four content classes through the content
       helpers (population, a moved part moves the digest, equal content
       is one value, pickle); a table and the ``Symbolic`` it lowers to
       are one right-hand side and two keys; seed-stable under
       ``PYTHONHASHSEED`` 1 and 2 with a ``str`` control
     - ``CI0`` (content collapsed to the type tag); ``D1`` (the written
       question stored); ``D2`` (the question dropped from the parts)
   * - S8.3
     - the layer: building two specifications in a fresh interpreter loads
       nothing above the input tier; an AST census admits only ``data``,
       ``geometry``, ``numerics`` and the package; the layer table's
       rows in ``tests/gates/test_layer_imports.py``
     - ``L1``, a file edit importing ``orpheus.transport.mesh.material_mesh``
   * - S8.4
     - the assembly rung (``l1``): a two-interval specification with a
       spectator builds a ``MaterialMesh``; the infinite medium's k from
       the 0-D solve and from an S\ :sub:`N` reflective slab equals the
       dense pencil :math:`\rho(A^{-1}F)` assembled in the test from the
       raw arrays, :math:`A = \mathrm{diag}(\Sigma_t) - \Sigma_{s0}^{\mathsf
       T} - 2\Sigma_2^{\mathsf T}`, on a fixture with upscatter and
       :math:`(n,2n)`
     - ``M1`` (``N2N_MULTIPLICITY`` = 1) reds both k rows; the mesh row is
       reddened by no arm, declared: a break there is the material mesh's
   * - S8.5
     - RECORD: the digests of two canonical specifications (k on the
       infinite medium; a fixed source on the slab) pinned as literals; a
       moved digest invalidates every cache key holding a specification
     - any encoder or schema edit
   * - S8.6
     - the channels: exactly three members of a plain ``Enum``; the
       carriage predicate over 4 fixtures × 3 channels; a block above
       :math:`P_0` alone is carried
     - ``K1`` (two channels swapped), ``K2`` (always carried), ``K5``
       (never carried); the member list by a file edit
   * - S8.7
     - the cell coefficient: a content value (order, repetition and numpy
       ids are not content); 5 construction refusals; the two fields and
       ``every`` equal to the field spelling; resolution per channel;
       idempotence; 4 resolution refusals; ``every(F)`` on materials
       carrying none is a zero direction
     - ``V1`` (no parse), ``V3`` (``every`` ignores its channels), ``K3``
       (``resolve`` the identity), ``K4`` (zero cells kept); idempotence by
       no arm, declared
   * - S8.8
     - the extent: a content value keyed by its index; 5 construction
       refusals (``-1``, a ``bool``, a ``float``, a ``str``, ``None``);
       resolution; the range is the intervals', not the distinct
       materials'
     - ``V2`` (no parse), ``E1`` (no range check), ``N2``
   * - S8.9
     - the group-count rule's one home (above)
     - ``C1`` (the rule reads nothing) reds the existing ``MaterialMesh``
       rows too, which proves the mesh routes through it; ``C2`` (the mesh
       reads the reachable materials only) reds the spectator row; the
       home and census rows were red on the pre-carve tree
   * - S8.10
     - canonicalisation: ``every`` resolves to the explicit non-zero cells
       in the parameter and in a point key; the stored question is unequal
       to the written one; four spellings of one direction are one
       specification; idempotent under re-posing and pickle; the materials
       are a ``Materials`` value; the two layers are two keys
     - ``D1``, ``K3``, ``V3``; ``K6`` (a point key keeps its zero cells)
       reds exactly the point-key row
   * - S8.11
     - ``Symbolic.depends_on``: 9 expressions × :math:`(r, \mu, \varphi)`,
       the joint predicate, agreement with ``is_isotropic``, one dependent
       group, and the refusal of an empty or non-owned query
     - ``Z1`` (free symbols), ``Z2`` (arguments ignored), ``Z3`` and
       ``Z4`` (constant answers), ``Z5`` (no owned-symbol refusal)

**The battery** (``scratch/reference_architecture/p1step8/ta/battery/final/``:
``step8_battery.py``, a ``-p`` plugin that installs one arm per test and
restores it at teardown, ``run_battery.sh``, ``logs/``, ``redsets.txt``)
ran on the real modules over 330 rows: the 204 step-8 rows and the
neighbouring files the arms reach (``test_mesh_free_function.py``,
``test_content_identity.py``, ``test_materials.py``,
``test_material_mesh.py``, ``test_snmesh_materials_pr_typed_0.py``).
`[M]` (``logs/none.log``, ``logs/P0.log``): the unmutated run 330 passed;
the positive control ``P0``, neither type checking anything, 65 failed.
Each of 39 in-process arms first proved it bites on a probe (0
uninstallable) and reddened its target rows, the red sets read against
the targets; ``L1``, the file-edit arm, reddened the two S8.3 rows and the
layer-table row for the package. The rows reddened by no arm are
declared, 15 of the 204: the ``CellCoefficient`` and ``GeometryExtent``
equal pairs and pickles, the seed row (whose teeth are its own ``str``
control, since a subprocess sees no plugin), the four fixture-balance
rows, the ``Enum`` member list, and the S8.9 home and census rows.


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
   * - 2026-10-02
     - **The reference specification: a question with its materials,
       keyed, and the layer it is posed at is its type.**
       :mod:`orpheus.specification` landed with
       ``InfiniteMediumSpecification(material_id, mixture, question)``,
       posed on energy alone, and ``GeometrySpecification(materials,
       geometry, question)``, admitted in a canonical form (spectators
       dropped, every key resolved to explicit non-zero cells, the datum
       fitted to the problem's groups, regions and readable coordinates),
       with the coordinates :class:`~orpheus.data.cells.Channel` and
       :class:`~orpheus.data.cells.CellCoefficient` in ``data`` and
       :class:`~orpheus.geometry.extent.GeometryExtent` in ``geometry``,
       :meth:`Symbolic.depends_on
       <orpheus.numerics.mesh_free_function.Symbolic.depends_on>`, and the
       group-count rule's one home, :meth:`Materials.uniform_group_count
       <orpheus.data.materials.Materials.uniform_group_count>`, with
       ``InconsistentMaterialsError`` moved from ``transport`` to
       ``data`` (no shim). Ruled the same day (the user): the cells are the
       three emission channels, an extent is one interval's width, the
       nuclide density is deferred, and the infinite medium is the point
       in phase space, so the first build's ``geometry | None`` and the
       proposed geometry value for the infinite medium were both retired
       for a second type. Record:
       :ref:`structured-geometry-specification`.
     - #405
     - ``29fe4266`` (the specification), ``eaa74163`` (the two types),
       ``2bc2e8c2`` (the re-review fixes), ``65b533fe`` (the physics docs)
   * - 2026-10-02
     - **What is asked is a physics-free value, and a real number has one
       definition.** :mod:`orpheus.numerics.question` landed with
       ``Eigen(parameter, point, mode)``, ``FixedSource(source, point)``,
       ``Response(detector, point)`` and the modes ``Fundamental`` and
       ``Nearest``, content-identity values admitted at construction.
       Ruled the same day (the user): a parameter's direction is a set of
       reaction-grid cells scaled together, the point is a frozen mapping
       of offsets from the physical value, and the adjoint fixed-source
       question is its own type, with no adjoint flag on any value. The
       2026-09-25 cases ``Eigen(k)``, ``Eigen(c)``, ``FixedSource(σ)`` and
       ``CriticalParameter``, each with a forward and an adjoint, were
       retired before any was built. ``orpheus/geometry/scalars.py``
       moved to :mod:`orpheus.numerics.scalars` (no shim) and became the
       one conversion the parsers, the encoder and ``Mixture`` make; the
       frozen encoder refuses a ``MappingProxyType`` part. Records:
       :ref:`structured-geometry-question-values`,
       :ref:`structured-geometry-one-real-parser`.
     - #405, #559, #561
     - ``6b179059``, ``31a2dc46``, ``8122abc6``, ``1fde7b59`` on ``main``
   * - 2026-10-02
     - **The source and the detector a specification states are
       mesh-free functions.** :mod:`orpheus.numerics.mesh_free_function`
       landed with ``RegionwiseConstant`` (a per-(region, group) table on
       the angle-integrated space) and ``Symbolic`` (one SymPy expression
       per group on phase space, stored as ``srepr`` text, the SymPy
       version a content part, the text parsed through an AST whitelist,
       only real scalar functions admitted, isotropy decided by
       substitution). Ruled the same day (the user): the :math:`4\pi` is
       the measure, so a table's role picks its arrow into phase space,
       the section for a source and the retraction's adjoint for a
       detector, and neither type carries a density.
       :attr:`CoordSystem.angular_chart
       <orpheus.geometry.coord.CoordSystem.angular_chart>` declares the
       chart per coordinate system (no azimuth reference on the sphere).
       Branch 1: :mod:`orpheus.derivations.common.angular_measure`. SymPy
       became a core dependency. Record:
       :ref:`structured-geometry-mesh-free-functions`.
     - #405
     - ``61f82a18`` on ``main``
   * - 2026-10-02
     - **Equality and hash are content, through one encoder.**
       :mod:`orpheus.numerics.content` landed: ``encode``,
       ``content_digest`` (blake2b-256), ``ContentlessError`` and the
       ``ContentIdentity`` mixin, with a schema tag per class (the
       user's ruling of 2026-10-01: a key covers the schema, so an
       entry written under an older schema misses) and canonical forms
       that follow ``==`` (the user's ruling of 2026-10-02). ``Mixture``,
       ``Materials`` (no longer compared by identity; ids coerced to
       ``int``; picklable), ``BC`` (hashable, read-only real-valued
       ``params``), every boundary law with ``LawSum`` and
       ``LawScaled``, ``StructuredGeometry``, ``FaceLaws`` (no longer
       equal to a plain ``dict``), ``CellEdges``, ``Mesh1D`` (hashable)
       and ``Mesh2D`` (read-only copied arrays, ``==`` works) moved onto
       it, with the ``Axis`` family and the five space-name digests.
       NaN is refused when a value is constructed. ``FrozenMapping``
       became the one frozen mapping (``Materials.mixtures``,
       ``BC.params``, and the base of ``FaceLaws``), ``name_digest`` the
       one space-name digest, and pickling goes through the constructor.
       The encoder refuses a mutable part, an identity-compared
       dataclass and a complex array. Retired:
       ``Axis._identity_key``, ``Axis._structural_bytes``, the axis's
       ``content_parts`` overrides, ``Mixture._identity_key``,
       ``Mesh1D``'s hand-written ``__eq__``, the ``MappingProxyType``
       stores and ``__reduce__`` methods of ``Materials`` and ``BC``,
       and the ``BC`` branch of ``MaterialMesh``'s law key. Gates:
       spec §1.5 of ``.claude/plans/reference_p1_spec.md``; record:
       :ref:`structured-geometry-content-identity`.
     - #405
     - ``a5113ac0`` on ``main``
   * - 2026-09-30
     - **CP's slab left face is the mirror its kernel computes.** The
       guard ``_refuse_a_law_cp_drops`` admitted a slab whose left law
       equalled its right law, on the claim that equal laws are what CP
       computes; Census A measured a white left law solved as a mirror
       (CP 1.212883 against S\ :sub:`N` reflective | white 1.212884 and
       white | white 1.212537), so the guard now admits only a reflective
       left law, and the CP gates declare it (the test helper
       ``cp_body``). Ruling: the user, 2026-09-30,
       ``.claude/plans/boundary_law_ontology.md``, "Third exchange".
     - #513
     - ``8f9300b2`` on ``main``
   * - 2026-09-30
     - **No boundary law is undeclared.** ``None`` retired as a boundary
       declaration everywhere: ``Mesh2D``'s four ``bc_*`` fields gave
       way to ``face_laws``, a mapping from face name to law over the
       inventory the 1-D topology law derives axis by axis; the axis
       classes' laws became required and parsed; one element parser,
       ``parse_boundary_law``, replaced ``_check_boundary_declaration``
       (which admitted ``None``). The S\ :sub:`N` entries'
       ``boundary_condition=`` parameter, the ``_apply_default_bcs``
       fill, the axes' ``with_uniform_bc`` and the shared resolver's
       reflective default retired; ``pwr_pin_2d`` takes a required
       ``law``. After review, both meshes store one ``FaceLaws`` (a
       frozen, ordered, picklable mapping from face name to law) over
       one inventory rule, ``face_inventory``: ``Mesh1D``'s positional
       tuple and ``boundary_faces`` (renamed ``boundary_points``) and
       ``Mesh2D.boundary_faces`` retired; ``wigner_seitz_pin_cell`` and
       ``pwr_slab_half_cell`` lost ``boundaries=``, their laws being
       part of the model. Rulings: the user, 2026-09-29 and 2026-09-30,
       P1 step 3c of ``.claude/plans/reference_cache.md``; the refused
       candidates are 5 to 9 of :ref:`structured-geometry-mesh-refuted`.
     - #405
     - ``4318adc5``
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
     - ``3a6468e6``
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
