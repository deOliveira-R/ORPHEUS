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
``precomputed_volumes`` override. The
:func:`~orpheus.geometry.factories.mesh1d_from_zones` factory sets
it by computing the *algebraic* cell volume (e.g.
:math:`V_{\rm cell} = \pi(r_{\rm out}^2 - r_{\rm in}^2)/n` in the
cylindrical case) and broadcasting that scalar to every cell in
the zone. Deriving volumes from the *edges* after the fact via
:func:`~orpheus.geometry.coord.compute_volumes_1d` would pass
through a ``sqrt → **2`` or ``cbrt → **3`` round trip that loses
roughly one ULP per cell and breaks the invariant "every cell in
an equal-volume zone is bit-identical" at ``rtol=1e-14``. Manually
constructed meshes with arbitrary edges still derive volumes from
edges as before — the override only kicks in on the factory path.

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
``kind`` means physically. Semantics are resolved by each solver's
augmented mesh (``SNMesh``, ``CPMesh``, ``MOCMesh``, ``MCMesh``,
``DiffusionSolver``) at construction time via a class-level
``BC_REGISTRY: dict[str, Callable]``. This registry maps kind
strings to factory functions that translate the abstract declaration
into solver-specific internal state (e.g. setting incoming angular
flux to zero for SN vacuum, or choosing a collision-probability
transform for CP white). If a mesh carries a ``kind`` that the
solver does not support, construction raises ``ValueError`` listing
the supported kinds.

This pattern has three advantages:

1. **Solver-agnostic problem setup.** The same ``Mesh1D`` with
   ``bc_right=BC.vacuum`` can be passed to SN, CP, or diffusion
   solvers without modification — each resolves the tag through
   its own registry.
2. **Extensibility.** Adding a new BC type (e.g. albedo, periodic)
   requires only a new factory function and a one-line addition
   to the solver's ``BC_REGISTRY``. No geometry code changes.
3. **Discoverability.** Each factory function carries a docstring
   that serves as a human-readable description, queryable at
   runtime via
   ``{k: v.__doc__ for k, v in SolverMesh.BC_REGISTRY.items()}``.

The current ``BC_REGISTRY`` contents per solver are:

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
   * - Diffusion
     - ``vacuum``, ``reflective``
     - ``vacuum``


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
:obj:`BC.vacuum`, :obj:`BC.reflective`, and :obj:`BC.white`.
These are ordinary ``BC`` instances, not subclasses — they exist
solely to avoid spelling out ``BC("vacuum")`` at every call site.

.. autoclass:: orpheus.geometry.mesh.BC
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:


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


Factories
---------

The factory layer is the recommended construction path.
:class:`~orpheus.geometry.factories.Zone` describes one material
region by its outer boundary and cell count;
:func:`~orpheus.geometry.factories.mesh1d_from_zones` builds a
coordinate-aware :class:`~orpheus.geometry.mesh.Mesh1D` from a
list of zones. Three subdivision strategies are baked in:

* **Cartesian** — equal-width cells
  :math:`x_k = x_0 + (k/n)\,(x_n - x_0)`.
* **Cylindrical** — equal-volume annuli
  :math:`r_k = \sqrt{r_0^2 + (k/n)\,(r_n^2 - r_0^2)}`.
* **Spherical** — equal-volume shells
  :math:`r_k = \sqrt[3]{r_0^3 + (k/n)\,(r_n^3 - r_0^3)}`.

Each returns both the edges and the exact per-cell volume (a
broadcast scalar) so the frozen :class:`Mesh1D` can be built with
``precomputed_volumes`` set — see the design principle above.

**Convenience constructors:**

* :func:`~orpheus.geometry.factories.pwr_slab_half_cell` — Cartesian
  3-zone (fuel / clad / coolant) half-cell with a reflective
  symmetry plane at :math:`x = 0`.
* :func:`~orpheus.geometry.factories.pwr_pin_equivalent` — cylindrical
  Wigner--Seitz equivalent pin cell. The square unit cell of side
  *pitch* is replaced by a cylinder of equal cross-sectional area,
  :math:`r_{\rm cell} = {\rm pitch} / \sqrt{\pi}`.
* :func:`~orpheus.geometry.factories.homogeneous_1d` — single-material
  uniform mesh for homogeneous-medium tests and analytical
  benchmarks.
* :func:`~orpheus.geometry.factories.slab_fuel_moderator` — 2-zone
  Cartesian slab (fuel / moderator) for classic L1 verification
  problems.
* :func:`~orpheus.geometry.factories.pwr_pin_2d` — 2-D Cartesian mesh
  with material IDs assigned by radial distance from the pin centre.

**Material ID convention:**
``2 = fuel``, ``1 = clad``, ``0 = coolant / moderator``. This
ordering matches the synthetic cross-section library used by the
L0 / L1 verification suites.

.. automodule:: orpheus.geometry.factories
   :members:
   :undoc-members:
   :show-inheritance:
   :noindex:
