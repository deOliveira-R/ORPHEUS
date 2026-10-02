Geometry Infrastructure
========================

The :mod:`orpheus.geometry` package is the geometry layer: the shapes,
the coordinate systems and the boundary declarations of a problem, with
no discretisation. It holds four things:

* :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`: the
  1-D layered geometry (a coordinate system, breakpoints, one material id
  per interval and one boundary law per boundary point), a pure shape and
  boundary description that carries no cell counts. Reference solvers
  consume it directly.
* :class:`~orpheus.geometry.coord.CoordSystem` and the volume and area
  formulas of each coordinate system (:mod:`orpheus.geometry.coord`).
* :mod:`orpheus.geometry.boundary`: the boundary-condition tag
  :class:`~orpheus.geometry.boundary.BC` and the typed boundary laws.
* :mod:`orpheus.geometry.transformation`: rigid motions of
  :math:`\mathbb{R}^d` and permutations.

The mesh, which divides a geometry's intervals into cells, is an overlay
on the geometry and lives in its own package, :mod:`orpheus.mesh`,
documented on :doc:`/api/mesh`. The dependency runs one way:
:mod:`orpheus.mesh` imports :mod:`orpheus.geometry`, and nothing in
:mod:`orpheus.geometry` imports :mod:`orpheus.mesh`
(:file:`tests/gates/test_layer_imports.py` enforces it; the layer
assignment is :ref:`architecture-layering`).

.. contents::
   :local:
   :depth: 2


Design Principles
-----------------

**Coordinate-aware volumes and surfaces.**
All geometric quantities route through
:mod:`orpheus.geometry.coord`, which dispatches on a
:class:`~orpheus.geometry.coord.CoordSystem` enum
(``CARTESIAN``, ``CYLINDRICAL``, ``SPHERICAL``). This keeps the
physics solvers coordinate-agnostic — the same
:func:`~orpheus.sn.solver.solve_sn` entry point handles slab,
cylinder, and sphere without branching on geometry.

**Boundary condition declaration and deferred resolution.**
Boundary conditions follow a two-phase pattern: *declare* on the
geometry, *resolve* at solver construction. The geometry layer
provides :class:`~orpheus.geometry.boundary.BC`, a frozen dataclass
carrying a ``kind`` string (e.g. ``"vacuum"``, ``"reflective"``,
``"white"``) and an optional ``params`` dict for numeric parameters
(e.g. ``{"albedo": 0.7}``). The meshes of :doc:`/api/mesh` store a
``BC`` tag or an already-typed boundary law on each boundary face:
:class:`~orpheus.mesh.structured.Mesh1D` has ``face_laws``, one per
boundary face, keyed by face name (``xmin`` and ``xmax``, or ``xmax``
alone on a solid cylinder or sphere);
:class:`~orpheus.mesh.structured.Mesh2D` has ``face_laws`` over its
faces ``xmin``, ``xmax``, ``ymin``, ``ymax`` (no ``xmin`` on a solid
:math:`(r, z)` mesh). Both are one value,
:class:`~orpheus.mesh.face_laws.FaceLaws`. A
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`, both
meshes and the S\ :sub:`N` axis primitives refuse ``None``, through one
element parser,
:func:`~orpheus.geometry.structured_geometry.parse_boundary_law`,
because a declaration states the problem and a default is a method's.
No method fills an unstated face: the shared resolution of
S\ :sub:`N` and diffusion reads the declared law and nothing else
(:ref:`structured-geometry-no-default-law`).

The geometry module makes **no assumptions** about what a given
``kind`` means physically. Semantics are resolved by each method's own
hub (``SNProblem`` — the SN **Problem**, ``SNMesh`` until #412,
2026-09-18 — ``DiffusionMesh``, and the future
``CPMesh`` / ``MOCMesh`` / ``MCMesh``) at construction time via a
class-level ``BOUNDARY_OPERATOR_REGISTRY: dict[str,
type[BoundaryTraceLaw]]`` mapping kind strings to typed boundary
LAWS, whose realized per-face operators land in the face-name-keyed
``bc`` dict each hub carries (``SNProblem.bc`` /
``DiffusionMesh.bc`` — #290 P7a moved the diffusion resolution off
the solver onto the phase space, the SN pattern). The realization
translates the abstract declaration into method-specific operator
state (e.g. zeroing the incoming angular flux for SN vacuum, or the
scalar albedo row :math:`J^- = \mathcal{A} J^+` for diffusion). If a
mesh carries a ``kind`` that the method does not support,
construction raises ``ValueError`` listing the supported kinds.

This pattern has three advantages:

1. **Solver-agnostic problem setup.** The same ``Mesh1D`` with a
   vacuum outer law can be passed to SN, CP, or diffusion
   solvers without modification — each method-mesh resolves the tag
   through its own registry.
2. **Extensibility.** Adding a new BC type (e.g. albedo, periodic)
   requires only a boundary-law class and a one-line addition to the
   method-mesh's ``BOUNDARY_OPERATOR_REGISTRY``. No geometry code
   changes.
3. **Discoverability.** Each boundary-law class carries a docstring
   that serves as a human-readable description, queryable at
   runtime via ``{k: v.__doc__ for k, v in
   SNProblem.BOUNDARY_OPERATOR_REGISTRY.items()}``.

The current registry contents per method are:

.. list-table:: Supported boundary conditions by solver
   :header-rows: 1
   :widths: 20 40 40

   * - Solver
     - Supported kinds
     - Default (when ``None``)
   * - SN
     - ``vacuum``, ``reflective``
     - ``reflective``
   * - CP
     - ``white``, ``vacuum``
     - ``white``
   * - MOC
     - ``reflective``
     - ``reflective``
   * - MC
     - ``periodic``
     - ``periodic``
   * - Diffusion (``DiffusionMesh``)
     - ``vacuum``, ``reflective``, ``albedo``, ``zero_flux``
     - ``reflective``

Boundary Conditions
-------------------

The :class:`~orpheus.geometry.boundary.BC` dataclass is the single type
used to declare boundary conditions on geometry surfaces. It is
defined in :mod:`orpheus.geometry.boundary` and exported from
:mod:`orpheus.geometry` for convenience:

.. code-block:: python

   from orpheus.geometry import BC, StructuredGeometry
   from orpheus.mesh import CellsByCount, Mesher

   # Pre-built convenience instances (tab-completable)
   bc_v = BC.vacuum       # BC("vacuum")
   bc_r = BC.reflective   # BC("reflective")
   bc_w = BC.white        # BC("white")

   # Custom BC with parameters
   bc_a = BC("albedo", params={"albedo": 0.7})

   # Declare each boundary point's law on the geometry; the mesher
   # carries the laws onto the mesh's named faces.
   geom = StructuredGeometry.slab((0.0, 10.0), (0,), left=BC.reflective, right=BC.vacuum)
   mesh = Mesher(geom).partition(CellsByCount.uniform_width(20)).mesh
   assert dict(mesh.face_laws) == {"xmin": BC.reflective, "xmax": BC.vacuum}
   assert mesh.outer_law == BC.vacuum

Three convenience class-level instances are pre-defined:
:obj:`BC.vacuum <orpheus.geometry.boundary.BC.vacuum>`,
:obj:`BC.reflective <orpheus.geometry.boundary.BC.reflective>`, and
:obj:`BC.white <orpheus.geometry.boundary.BC.white>`.
These are ordinary ``BC`` instances, not subclasses — they exist
solely to avoid spelling out ``BC("vacuum")`` at every call site.

.. Deliberately INDEXED (no ``:noindex:``, unlike its neighbours): this is
   the one directive that registers ``BC`` and its three tag constants in
   the Python domain, so ``:obj:`BC.vacuum <orpheus.geometry.boundary.BC.vacuum>```
   above — and every other qualified reference to them — resolves to a link
   instead of rendering as plain text.  The ``automodule`` below keeps its
   ``:noindex:``, so it re-renders ``BC`` without claiming the entry, and
   there is no duplicate-object warning.  (#302's general case is unfixed:
   `[M]` 2026-08-10 the whole built inventory is 1014 entries because nearly
   every api page is ``:noindex:``.)

.. autoclass:: orpheus.geometry.boundary.BC
   :members:
   :undoc-members:
   :show-inheritance:


Coordinate Systems
------------------

:mod:`orpheus.geometry.coord` defines the
:class:`~orpheus.geometry.coord.CoordSystem` enum and the
coordinate-aware volume / surface primitives:

* ``compute_volumes_1d(coord, edges)``
* ``compute_areas_1d(coord, edges)``
* ``compute_volumes_2d(coord, edges_x, edges_y)``

All three dispatch on ``coord`` and return NumPy arrays sized to
match the mesh. The 1-D spherical volume formula,

.. math::

   V_i = \frac{4\pi}{3}\bigl(r_{i+1}^3 - r_i^3\bigr),

and the cylindrical formula,

.. math::

   V_i = \pi\bigl(r_{i+1}^2 - r_i^2\bigr),

are the standard shell / annulus expressions. The surface arrays
return :math:`4\pi r^2` (spherical) or :math:`2\pi r` (cylindrical,
per unit height) at each edge — these drive the :math:`\Delta A /
w_m` redistribution factor in the curvilinear SN sweeps (see
:ref:`theory-discrete-ordinates`).

Each coordinate system also declares the angular chart
:math:`(\mu, \varphi)` a function on phase space reads its direction in,
:attr:`CoordSystem.angular_chart
<orpheus.geometry.coord.CoordSystem.angular_chart>`, an
:class:`~orpheus.geometry.coord.AngularChart` over the columns of the local
direction frame (the theory, with why the sphere declares no azimuth
reference: :ref:`structured-geometry-angular-chart`).

.. automodule:: orpheus.geometry.coord
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Structured geometry
-------------------

A :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
is the pure shape of a 1-D problem: its coordinate system ``coord``
(a :class:`~orpheus.geometry.coord.CoordSystem` member), its
``breakpoints`` :math:`r_0 < \dots < r_R` (stored bit for bit as given),
one material id per interval in ``mat_ids``, and one law per boundary
point in ``boundaries`` (a :class:`~orpheus.geometry.boundary.BC` tag or
a typed law; ``None`` is refused). The boundary points are derived, not
declared: two on a slab, the outer surface alone on a solid cylinder or
sphere (whose centre is an interior point and carries no law), inner
and outer on a hollow one. It carries **no** cell counts. The reasons
for each of these choices are on
:doc:`/theory/foundations/structured_geometry`. Reference solution generators (``Billiard``, ``MomentSpace``,
``Spectrum``, ``BasisSpace``) consume it directly, because they need no
mesh; the discrete production solvers consume the
:class:`~orpheus.mesh.structured.Mesh1D` that a
:class:`~orpheus.mesh.mesher.Mesher` builds from it (the construction
path is on :doc:`/api/mesh`). The named-face constructors
``StructuredGeometry.slab(breakpoints, mat_ids, left=, right=)``,
``cylinder`` and ``sphere`` (``outer=``, and ``inner=`` exactly when the
body is hollow), ``uniform_boundary(coord, breakpoints, mat_ids, law)``
and ``from_homogeneous(width, boundary)`` name each boundary point's law
at the call site.

A configuration published as thicknesses is built with
:meth:`StructuredGeometry.from_thicknesses
<orpheus.geometry.structured_geometry.StructuredGeometry.from_thicknesses>`,
the left fold :math:`r_{k+1} = r_k + t_k`. Two conventional PWR shapes
ship as :class:`StructuredGeometry` classmethods:

* :meth:`StructuredGeometry.pwr_slab_half_cell
  <orpheus.geometry.structured_geometry.StructuredGeometry.pwr_slab_half_cell>`
  — Cartesian 3-region (fuel / clad / coolant) half-cell starting at
  the reflective symmetry plane :math:`x = 0`. Both faces are
  reflective, and that is part of the model (both are symmetry planes
  of the lattice), not a parameter.
* :meth:`StructuredGeometry.wigner_seitz_pin_cell
  <orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell>`
  — cylindrical Wigner--Seitz equivalent pin cell. The square unit
  cell of side *pitch* is replaced by a cylinder of equal
  cross-sectional area, :math:`r_{\rm cell} = {\rm pitch} /
  \sqrt{\pi}`. The outer law is white, and that is part of the model
  (isotropic re-entry is what maps the lattice to one cylindrical
  cell), not a parameter.

Neither named cell takes a law: the same stack under another law is a
different body, built with ``StructuredGeometry.slab`` or
``StructuredGeometry.cylinder``.

**Material ID convention:**
``2 = fuel``, ``1 = clad``, ``0 = coolant / moderator``. This
ordering matches the synthetic cross-section library used by the
L0 / L1 verification suites.

.. automodule:: orpheus.geometry.structured_geometry
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
