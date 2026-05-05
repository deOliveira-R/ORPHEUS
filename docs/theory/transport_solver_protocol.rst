.. _theory-transport-solver-protocol:

==========================
Transport Solver Protocol
==========================

The :class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol is the structural unification of the project's
**math-heart pattern** for analytical reference solvers
(``Billiard``, ``MomentSpace``, the parallel ``Spectrum`` and
``BasisSpace``). They answer the same question — *given materials
and a geometry, what is the critical configuration?* — and the
Protocol gives them a shared attribute / method surface so
cross-method comparators can substitute one for another without
case-by-case ``isinstance`` branching.

Production discrete-mesh solvers (:class:`~orpheus.cp.solver.CPSolver`,
:class:`~orpheus.sn.solver.SNSolver`) deliberately stay **off** the
Protocol. They consume ``(materials, mesh, params)`` directly via the
canonical :func:`~orpheus.cp.solver.solve_cp` /
:func:`~orpheus.sn.solver.solve_sn` entry points. Forcing them under
the Protocol facade hid solver tuning behind a thin wrapper and
bundled cross-cutting parameters into a redundant struct; production
callers chew through 1000-group industrial data and need direct
access to the solver-tuning surface. The geometry → mesh transition
(:meth:`~orpheus.geometry.mesh.Mesh1D.from_geometry`) is an explicit
caller-side step, not a hidden one inside the solver.

.. contents:: On this page
   :local:
   :depth: 2


Three-axis carve
================

The Protocol's design follows a deliberate three-axis carving
of the transport problem:

1. **What** — the materials (``dict[int, Mixture]``) and the
   geometry (``GeometrySpec``). These are *method-agnostic* — the
   same problem description is served to every Protocol-conformant
   reference solver (F_N moment space, trajectory resolvent, ...).
2. **How** — the method-specific computational specialisation. F_N
   commits to a moment basis order :math:`N`; trajectory_resolvent
   commits to a billiard table classified by orbit-space rank.
   Each Protocol-conforming class owns its own method kwargs;
   the Protocol does NOT dictate them.
3. **What is asked** — :meth:`solve_critical` returns a
   :class:`~orpheus.derivations.common.solution_types.CriticalSolution`.
   This is method-agnostic: the same dataclass is populated by
   every conforming class.

Why this carving matters: the *what* and *what-is-asked* axes
are the **substrate** the Protocol lives on; the *how* axis
is each class's private business. Forcing the *how* into a
common surface would couple unrelated methods (F_N's ``n_modes``
has no meaning for a billiard rank). Keeping it private respects
each method's mathematical structure.


Protocol surface (ABNF)
========================

The full surface every conforming class must expose:

.. code-block:: abnf

   TransportSolver = materials geometry_spec method_name solve_critical

   materials      = "materials" ":" Mapping<int, Mixture>
   geometry_spec  = "geometry_spec" ":" GeometrySpec
   method_name    = "method_name" ":" string
   solve_critical = "def solve_critical(self) -> CriticalSolution"

That is the *entire* contract. Method-specific kwargs
(``fn_order``, ``quadrature``, ``flux_reconstruction``) are carried on
each class as their own fields; the Protocol does not see them. This
deliberate minimalism follows from the *unify-after-two-instances*
memory (``feedback_unify_after_two_instances.md``): the Protocol was
designed only after both ``Billiard`` and ``MomentSpace`` had shipped,
so its surface matches what they already exposed.


Why a Protocol, not an ABC
==========================

PEP-544 :class:`typing.Protocol` is **structural typing**. Any
class with the four named members above automatically conforms,
without inheriting from a base. This matches the project's
reality:

* ``Billiard`` and ``MomentSpace`` were designed independently
  before the Protocol; an ABC would have required retroactive
  subclassing.
* The Protocol is :func:`runtime_checkable`, so
  ``isinstance(x, TransportSolver)`` works for the schema-gate
  tests in :mod:`tests.derivations.test_transport_solver_protocol`
  without a registration step.

Production discrete-mesh solvers (CPSolver, SNSolver) deliberately
stay off the Protocol — see § *Production solvers stay direct*
below.


Protocol vs adapter layer
==========================

The cross-method gates in
:mod:`tests.cross_method.test_eigenvalue` already use a
**per-method adapter layer** (``FNSlabAdapter``,
``TrajectoryResolventSphereAdapter``, ...) that handles unit
conversions (mfp ↔ cm, half-thickness ↔ full slab) and
parameter selection. This adapter layer pre-dates the Protocol
and **stays in place** — the test bodies in
``test_eigenvalue.py`` continue to consume adapters, the 84
test IDs are preserved (CI / pytest-xdist contract).

What the Protocol adds is a *new* dispatch path:
:func:`Billiard.from_problem` / :func:`MomentSpace.from_problem`
accept the same production-protocol input pair
``(materials, geometry_spec)`` and return Protocol-conforming
objects directly, without going through an adapter.

The two layers serve different audiences:

* **Adapter layer** — cross-method test bodies that need
  unit conversion + parameter selection from a
  :class:`CrossMethodCase`. Stays as-is.
* **Protocol layer** — downstream consumers that already have
  ``(materials, geometry_spec)`` in production form.

The :mod:`tests.cross_method.test_polymorphism` foundation
gate pins that the two layers agree on the canonical case
sets — adapter ↔ Protocol routes converge on the same
function-level call and produce the same float result.


Concrete realisations
=====================

Four Protocol-conforming continuous-reference classes ship today:

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Class
     - Pillar
     - Notes
   * - :class:`~orpheus.derivations.continuous.trajectory_resolvent.Billiard`
     - trajectory_resolvent
     - Birkhoff transfer-operator resolvent on a billiard. Six
       geometries (slab, sphere, cylinder, slab-asymmetric,
       hollow_sphere, annulus) at two orbit-space classes.
   * - :class:`~orpheus.derivations.continuous.fn_method.moment_space.MomentSpace`
     - fn_method
     - Galerkin half-range projection on F_N moments. Slab and
       sphere (cylinder out of pillar — see
       :mod:`...singular_eigenfunction.cylinder`).
   * - ``Spectrum``
     - singular_eigenfunction
     - Case singular-eigenfunction expansion (parallel agent —
       see :mod:`orpheus.derivations.continuous.singular_eigenfunction`).
   * - ``BasisSpace``
     - galerkin_spectral
     - Dahl-Sjöstrand Legendre-Galerkin block-matrix linearisation
       (parallel agent — see
       :mod:`orpheus.derivations.continuous.galerkin_spectral`).

The math-rich theory for the continuous-reference classes lives
on each method's existing theory page; the
:doc:`reference_solvers` § "Three meanings of the Green's
function" places ``Billiard`` (path-integral) and
``MomentSpace`` (spectral) within the same Green's-function
landscape.


Production solvers stay direct
==============================

Production discrete-mesh solvers (:class:`~orpheus.cp.solver.CPSolver`,
:class:`~orpheus.sn.solver.SNSolver`) deliberately stay **off** the
Protocol. The canonical entry points are the free functions
:func:`~orpheus.cp.solver.solve_cp` /
:func:`~orpheus.sn.solver.solve_sn`, which consume
``(materials, mesh, params)`` directly:

.. code-block:: python

    from orpheus.cp.solver import CPParams, solve_cp
    from orpheus.geometry import Mesh1D
    from orpheus.geometry.structured_geometry import StructuredGeometry

    # Build the mesh explicitly. Production callers control
    # cell counts per region — multi-region geometries via
    # Mesh1D.from_geometry, single-region via the geometry
    # factory helpers.
    geom = StructuredGeometry.sphere(R_cm=5.0, mat_id=0, bc="white")
    mesh = Mesh1D.from_geometry(geom, region_meshes=({"n_cells": 50},))

    # Solver tuning lives on the production-side params struct.
    params = CPParams(n_quad_y=64, max_outer=500, keff_tol=1e-6)

    result = solve_cp({0: pu_mixture}, mesh=mesh, params=params)
    # result.keff, result.flux, result.keff_history, ...

Why CP / SN stay off the Protocol:

* Production solvers chew through 1000-group industrial data;
  bundling their cross-cutting parameters
  (``n_cells``, ``n_quad_y``, ``max_outer``, ``keff_tol``,
  ``solver_mode``) behind a Protocol facade hid solver tuning
  behind a thin wrapper.
* The geometry → mesh transition is an explicit caller-side
  step (:meth:`~orpheus.geometry.mesh.Mesh1D.from_geometry`).
  Production solvers stay focused on the solve, not on geometry
  management.
* Mismatches between geometry and angular quadrature
  (e.g. cylindrical SN requires level-symmetric / product
  quadrature, not 1-D Gauss-Legendre) are surfaced at the
  caller — a Protocol facade silently picked Gauss-Legendre
  regardless of geometry, masking the bug.


Foundation-tier verification
============================

Conformance to the Protocol is a **software invariant**, not
a physics claim. The L1 chains stay where they were
(``test_fn_la13511_*.py``,
``test_peierls_greens_function_*.py``,
:mod:`tests.cross_method.test_eigenvalue`). The Protocol unifies
the *shape* those L1-verified components must expose.

Foundation-tier tests:

* :mod:`tests.derivations.test_transport_solver_protocol` —
  conformance + ``KNOWN_TRANSPORT_SOLVERS`` registry +
  attribute-surface drift detection. Covers ``Billiard`` and
  ``MomentSpace`` (the four continuous-reference pillars). The
  registry-completeness test also pins that production
  ``cp`` / ``sn`` tags are NOT in the registry — production
  solvers stay off the Protocol.
* :mod:`tests.derivations.test_trajectory_resolvent_billiard` —
  11 tests. Bit-equality between :meth:`Billiard.solve_critical`
  and the underlying ``solve_greens_function_*`` entry points
  via :func:`np.array_equal` (NEVER :func:`np.allclose`). Covers
  sphere/cylinder/slab/slab_asymmetric in 1G + sphere MG.
* :mod:`tests.derivations.test_moment_space_from_problem_unified`
  — 6 tests. Bit-equality with the function-level F_N API.
* :mod:`tests.cross_method.test_polymorphism` — 6 tests.
  Adapter ↔ Protocol agreement on canonical cases.

The legacy↔new factory bit-equality test
(``test_billiard_from_problem_unified``) was removed when the
input-cleanup track sunset the dual-factory dispatch on
:meth:`Billiard.from_problem` — with one path remaining, the file
became a self-comparison. The structural-independence coverage
that remains lives in
:mod:`tests.derivations.test_trajectory_resolvent_billiard` (facade
vs underlying solver) and :mod:`tests.cross_method.test_polymorphism`
(Protocol vs adapter layer).


Structural-independence preservation
=====================================

The Protocol does **not** introduce new in-house cross-method
dependencies. Per :doc:`/skills/algebra-of-record` § "Structural
independence applies above the trusted-library line":

* Continuous-reference solvers (``Billiard``, ``MomentSpace``,
  ``Spectrum``, ``BasisSpace``) consume scipy / numpy primitives
  but **do not share** any in-house module above the
  trusted-library line.
* Production discrete-mesh solvers (``CPSolver``, ``SNSolver``)
  consume their own scipy / numpy primitives via separate code
  paths and do not share in-house modules with the
  continuous-reference family above the trusted-library line.

The Protocol is **above** the trusted-library line — but it
carries only attribute / method shapes, not numerical primitives.
Bit-equality preservation through the Protocol facade does
not propagate any computation; it merely keeps the unified
result-type semantics intact. Branch-1 ⊥ Branch-2 cross-checks
remain valid by construction.


Cross-references
================

* :doc:`reference_solvers` § "Three meanings of the Green's
  function" — the broader Green's-function landscape that
  ``Billiard`` and ``MomentSpace`` realise differently.
* :doc:`/skills/algebra-of-record` § "Branch 1 / Branch 2
  separation" — the Protocol lives at the math layer.
* :doc:`/skills/vv-principles` § "Hierarchical claim taxonomy"
  — software-invariant tier vs L0/L1 physics claims.

.. todo::
   Archivist expansion: add a worked example of substituting one
   pillar for another in a cross-method comparator (e.g.,
   ``MomentSpace`` swapped in for ``Billiard`` on a closed-sphere
   k_inf case). The example would walk through the structural-
   typing substitution with annotated code, and would highlight
   how the unified ``CriticalSolution`` field set lets the
   downstream comparator stay pillar-blind.
