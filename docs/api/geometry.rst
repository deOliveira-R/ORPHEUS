Geometry Infrastructure
========================

The :mod:`orpheus.geometry` package provides the spatial data
structures that every deterministic solver (SN, CP, MOC, diffusion)
consumes. A *mesh* is an immutable description of a domain — cell
edges, material assignments, and derived quantities such as
volumes and surfaces. Solvers receive a mesh and build mutable,
solver-specific state on top of it.

.. contents::
   :local:
   :depth: 2


Design Principles
-----------------

**Frozen dataclasses.**
Both :class:`~orpheus.geometry.mesh.Mesh1D` and
:class:`~orpheus.geometry.mesh.Mesh2D` are
``@dataclass(frozen=True)``. Once constructed, their fields cannot
be reassigned. This turns every solver entry point into a pure
function of its inputs and prevents whole classes of bugs where a
downstream routine accidentally mutates mesh state shared across
iterations.

**Coordinate-aware volumes and surfaces.**
All geometric quantities route through
:mod:`orpheus.geometry.coord`, which dispatches on a
:class:`~orpheus.geometry.coord.CoordSystem` enum
(``CARTESIAN``, ``CYLINDRICAL``, ``SPHERICAL``). This keeps the
physics solvers coordinate-agnostic — the same
:func:`~orpheus.sn.solver.solve_sn` entry point handles slab,
cylinder, and sphere without branching on geometry.

**Equal-volume subdivision.**
Curvilinear zones (cylindrical, spherical) are subdivided into
**equal-volume** annuli / shells rather than equal-width cells.
This gives uniform statistical weighting across the zone and
avoids skinny inner cells that would dominate the CFL-like step
limits of explicit sweeps.

**Precomputed volumes — the ULP escape hatch.**
:class:`~orpheus.geometry.mesh.Mesh1D` accepts an optional
``precomputed_volumes`` override.
:meth:`~orpheus.geometry.mesh.Mesh1D.from_geometry` sets it for every
``"equal-volume"`` region by routing through
:func:`~orpheus.geometry.factories._subdivide_zone`, which returns the
*algebraic* cell volume (e.g.
:math:`V_{\rm cell} = \pi(r_{\rm out}^2 - r_{\rm in}^2)/n` in the
cylindrical case) broadcast as a scalar to every cell in the region.
Deriving volumes from the *edges* after the fact via
:func:`~orpheus.geometry.coord.compute_volumes_1d` would pass
through a ``sqrt → **2`` or ``cbrt → **3`` round trip that loses
roughly one ULP per cell and breaks the invariant "every cell in
an equal-volume region is bit-identical" at ``rtol=1e-14`` (ERR-020).
``"uniform"`` regions and manually constructed meshes with arbitrary
edges still derive volumes from edges — the override only kicks in on
the equal-volume subdivision path.

**Boundary condition declaration and deferred resolution.**
Boundary conditions follow a two-phase pattern: *declare* on the
geometry, *resolve* at solver construction. The geometry layer
provides :class:`~orpheus.geometry.mesh.BC`, a frozen dataclass
carrying a ``kind`` string (e.g. ``"vacuum"``, ``"reflective"``,
``"white"``) and an optional ``params`` dict for numeric parameters
(e.g. ``{"albedo": 0.7}``). The meshes store ``BC | None`` on each
boundary face — :class:`~orpheus.geometry.mesh.Mesh1D` has
``bc_left`` and ``bc_right``;
:class:`~orpheus.geometry.mesh.Mesh2D` has ``bc_xmin``,
``bc_xmax``, ``bc_ymin``, and ``bc_ymax``. A value of ``None``
means "use the solver's default," which varies by method (e.g.
reflective for SN eigenvalue, white for CP).

The geometry module makes **no assumptions** about what a given
``kind`` means physically. Semantics are resolved by each method's
augmented mesh (``SNMesh``, ``DiffusionMesh``, and the future
``CPMesh`` / ``MOCMesh`` / ``MCMesh``) at construction time via a
class-level ``BOUNDARY_OPERATOR_REGISTRY: dict[str,
type[BoundaryTraceLaw]]`` mapping kind strings to typed boundary
LAWS, whose realized per-face operators land in the face-name-keyed
``bc`` dict each method-mesh carries (``SNMesh.bc`` /
``DiffusionMesh.bc`` — #290 P7a moved the diffusion resolution off
the solver onto the phase space, the SN pattern). The realization
translates the abstract declaration into method-specific operator
state (e.g. zeroing the incoming angular flux for SN vacuum, or the
scalar albedo row :math:`J^- = \mathcal{A} J^+` for diffusion). If a
mesh carries a ``kind`` that the method does not support,
construction raises ``ValueError`` listing the supported kinds.

This pattern has three advantages:

1. **Solver-agnostic problem setup.** The same ``Mesh1D`` with
   ``bc_right=BC.vacuum`` can be passed to SN, CP, or diffusion
   solvers without modification — each method-mesh resolves the tag
   through its own registry.
2. **Extensibility.** Adding a new BC type (e.g. albedo, periodic)
   requires only a boundary-law class and a one-line addition to the
   method-mesh's ``BOUNDARY_OPERATOR_REGISTRY``. No geometry code
   changes.
3. **Discoverability.** Each boundary-law class carries a docstring
   that serves as a human-readable description, queryable at
   runtime via ``{k: v.__doc__ for k, v in
   SNMesh.BOUNDARY_OPERATOR_REGISTRY.items()}``.

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

The :class:`~orpheus.geometry.mesh.BC` dataclass is the single type
used to declare boundary conditions on geometry surfaces. It is
exported from :mod:`orpheus.geometry` for convenience:

.. code-block:: python

   from orpheus.geometry import BC, Mesh1D, CoordSystem
   import numpy as np

   # Pre-built convenience instances (tab-completable)
   bc_v = BC.vacuum       # BC("vacuum")
   bc_r = BC.reflective   # BC("reflective")
   bc_w = BC.white        # BC("white")

   # Custom BC with parameters
   bc_a = BC("albedo", params={"albedo": 0.7})

   # Attach to mesh faces — None means "use solver default"
   mesh = Mesh1D(
       edges=np.linspace(0, 10, 21),
       mat_ids=np.zeros(20, dtype=int),
       coord=CoordSystem.CARTESIAN,
       bc_left=BC.reflective,
       bc_right=BC.vacuum,
   )

Three convenience class-level instances are pre-defined:
:obj:`BC.vacuum <orpheus.geometry.mesh.BC.vacuum>`,
:obj:`BC.reflective <orpheus.geometry.mesh.BC.reflective>`, and
:obj:`BC.white <orpheus.geometry.mesh.BC.white>`.
These are ordinary ``BC`` instances, not subclasses — they exist
solely to avoid spelling out ``BC("vacuum")`` at every call site.

.. Deliberately INDEXED (no ``:noindex:``, unlike its neighbours): this is
   the one directive that registers ``BC`` and its three tag constants in
   the Python domain, so ``:obj:`BC.vacuum <orpheus.geometry.mesh.BC.vacuum>```
   above — and every other qualified reference to them — resolves to a link
   instead of rendering as plain text.  The ``automodule`` below keeps its
   ``:noindex:``, so it re-renders ``BC`` without claiming the entry, and
   there is no duplicate-object warning.  (#302's general case is unfixed:
   `[M]` 2026-08-10 the whole built inventory is 1014 entries because nearly
   every api page is ``:noindex:``.)

.. autoclass:: orpheus.geometry.mesh.BC
   :members:
   :undoc-members:
   :show-inheritance:


Mesh
----

.. automodule:: orpheus.geometry.mesh
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Coordinate Systems
------------------

:mod:`orpheus.geometry.coord` defines the
:class:`~orpheus.geometry.coord.CoordSystem` enum and the
coordinate-aware volume / surface primitives:

* ``compute_volumes_1d(coord, edges)``
* ``compute_surfaces_1d(coord, edges)``
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

.. automodule:: orpheus.geometry.coord
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


Construction — the geometry to mesh path
----------------------------------------

The recommended 1-D construction path is **two-layered**: declare a
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
(pure shape — geometry kind, an ordered tuple of
:class:`~orpheus.geometry.structured_geometry.Region` layers, and the
endpoint :class:`~orpheus.geometry.mesh.BC`\ s), then discretize it
with :meth:`Mesh1D.from_geometry
<orpheus.geometry.mesh.Mesh1D.from_geometry>` by supplying one
:class:`~orpheus.geometry.mesh.RegionMesh` per region. The geometry
carries **no** cell counts; the discretization description enters at
exactly one point, the ``from_geometry`` call.

.. code-block:: python

   from orpheus.geometry import BC, Mesh1D, RegionMesh
   from orpheus.geometry.structured_geometry import (
       Region, StructuredGeometry,
   )

   geom = StructuredGeometry(
       geometry="SPH",
       regions=(Region(mat_id=0, outer_thickness_cm=5.0),),
       bcs=(BC.vacuum,),
   )
   mesh = Mesh1D.from_geometry(
       geom, region_meshes=(RegionMesh(n_cells=64),),
   )

This split is what lets **reference** solution generators (``Billiard``,
``MomentSpace``, ``Spectrum``, ``BasisSpace``) consume the
:class:`StructuredGeometry` directly — they need no mesh — while
**discrete production** solvers (``solve_sn`` / ``solve_cp`` /
``solve_moc`` / ``solve_mc``) consume the resulting
:class:`~orpheus.geometry.mesh.Mesh1D`.

Two conventional PWR shapes ship as
:class:`StructuredGeometry` classmethods:

* :meth:`StructuredGeometry.pwr_slab_half_cell
  <orpheus.geometry.structured_geometry.StructuredGeometry.pwr_slab_half_cell>`
  — Cartesian 3-region (fuel / clad / coolant) half-cell starting at
  the reflective symmetry plane :math:`x = 0`.
* :meth:`StructuredGeometry.wigner_seitz_pin_cell
  <orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell>`
  — cylindrical Wigner--Seitz equivalent pin cell. The square unit
  cell of side *pitch* is replaced by a cylinder of equal
  cross-sectional area, :math:`r_{\rm cell} = {\rm pitch} /
  \sqrt{\pi}`.

Per-region subdivision
~~~~~~~~~~~~~~~~~~~~~~

:class:`~orpheus.geometry.mesh.RegionMesh` selects the scheme per
region. ``"equal-volume"`` (the default) routes through
:func:`~orpheus.geometry.factories._subdivide_zone`, which carries the
three coordinate-system invariants:

* **Cartesian** — equal-width cells
  :math:`x_k = x_0 + (k/n)\,(x_n - x_0)`.
* **Cylindrical** — equal-volume annuli
  :math:`r_k = \sqrt{r_0^2 + (k/n)\,(r_n^2 - r_0^2)}`.
* **Spherical** — equal-volume shells
  :math:`r_k = \sqrt[3]{r_0^3 + (k/n)\,(r_n^3 - r_0^3)}`.

It returns both the edges **and** the exact per-cell volume (a
broadcast scalar) so the frozen :class:`Mesh1D` is built with
``precomputed_volumes`` set — see the design principle above.
``"uniform"`` instead lays down equal radial extents and derives
volumes from the edges via
:func:`~orpheus.geometry.coord.compute_volumes_1d`.

**Material ID convention:**
``2 = fuel``, ``1 = clad``, ``0 = coolant / moderator``. This
ordering matches the synthetic cross-section library used by the
L0 / L1 verification suites.

.. note:: **Retired 1-D factory surface.**

   Phase F retired the free functions ``Zone``, ``mesh1d_from_zones``,
   ``pwr_pin_equivalent``, ``pwr_slab_half_cell``, ``homogeneous_1d``
   and ``slab_fuel_moderator`` from
   :mod:`orpheus.geometry.factories`. Their jobs are now split between
   the geometry layer (the two classmethods above) and the mesh layer
   (:meth:`Mesh1D.from_geometry
   <orpheus.geometry.mesh.Mesh1D.from_geometry>` +
   :class:`~orpheus.geometry.mesh.RegionMesh`) — which is what removed
   the old "a factory both shapes AND meshes the problem" conflation.
   A homogeneous or two-region slab is now a one-liner
   ``StructuredGeometry`` literal, so no dedicated convenience
   function survives for it.

Structured geometry
~~~~~~~~~~~~~~~~~~~

.. automodule:: orpheus.geometry.structured_geometry
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:

Two-dimensional factories
~~~~~~~~~~~~~~~~~~~~~~~~~

There is no 2-D :class:`StructuredGeometry` yet, so the 2-D Cartesian
path keeps a standalone factory:
:func:`~orpheus.geometry.factories.pwr_pin_2d` builds a
:class:`~orpheus.geometry.mesh.Mesh2D` on a uniform grid with material
IDs assigned by radial distance from the pin centre.

.. automodule:: orpheus.geometry.factories
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
