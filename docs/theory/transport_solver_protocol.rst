.. _theory-transport-solver-protocol:

==========================
Transport Solver Protocol
==========================

The :class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol is the structural unification of the project's
**math-heart pattern** with its **production discrete-mesh
solvers**. Both classes of object — the analytical reference
solvers (``Billiard``, ``MomentSpace``, the parallel ``Spectrum``
and ``BasisSpace``) and the production CP/SN adapters
(``CPMeshAdapter``, ``SNMeshAdapter``) — answer the same question:
*given materials and a geometry, what is the critical configuration?*
The Protocol gives them a shared attribute / method surface so
cross-method comparators can substitute one for another without
case-by-case ``isinstance`` branching.

.. contents:: On this page
   :local:
   :depth: 2


Three-axis carve
================

The Protocol's design follows a deliberate three-axis carving
of the transport problem:

1. **What** — the materials (``dict[int, Mixture]``) and the
   geometry (``GeometrySpec``). These are *method-agnostic* — the
   same problem description is served to a CP solver, an SN solver,
   and an F_N moment space.
2. **How** — the method-specific computational specialisation. F_N
   commits to a moment basis order :math:`N`; trajectory_resolvent
   commits to a billiard table classified by orbit-space rank;
   discrete CP/SN commit to a discretisation
   (:class:`~orpheus.derivations.common.discretization_spec.DiscretizationSpec`).
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
has no meaning for a CP solver; CP's ``n_chord_quad`` has no
meaning for F_N). Keeping it private respects each method's
mathematical structure.


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
(``fn_order``, ``quadrature``, :class:`DiscretizationSpec`,
``flux_reconstruction``) are carried on each class as their own
fields; the Protocol does not see them. This deliberate
minimalism follows from the *unify-after-two-instances* memory
(``feedback_unify_after_two_instances.md``): the Protocol was
designed only after both ``Billiard`` and ``MomentSpace`` had
shipped, so its surface matches what they already exposed.


Why a Protocol, not an ABC
==========================

PEP-544 :class:`typing.Protocol` is **structural typing**. Any
class with the four named members above automatically conforms,
without inheriting from a base. This matches the project's
reality:

* ``Billiard`` and ``MomentSpace`` were designed independently
  before the Protocol; an ABC would have required retroactive
  subclassing.
* Discrete-mesh adapters are thin facades over production
  solvers (``solve_cp``, ``solve_sn``); the Protocol lets them
  conform without modifying the production solver's API.
* The Protocol is :func:`runtime_checkable`, so
  ``isinstance(x, TransportSolver)`` works for the schema-gate
  tests in :mod:`tests.derivations.test_transport_solver_protocol`
  without a registration step.


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

Six Protocol-conforming classes / adapters ship today:

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
   * - :class:`~tests.cross_method.discrete_adapters.CPMeshAdapter`
     - cp
     - Production CP solver behind the Protocol. ``solve_cp``
       wraps + :class:`CPParams` parameterisation.
   * - :class:`~tests.cross_method.discrete_adapters.SNMeshAdapter`
     - sn
     - Production SN solver behind the Protocol. 1-D
       Gauss-Legendre only today.

The math-rich theory for the continuous-reference classes lives
on each method's existing theory page; the
:doc:`reference_solvers` § "Three meanings of the Green's
function" places ``Billiard`` (path-integral) and
``MomentSpace`` (spectral) within the same Green's-function
landscape.


Discretisation parameters for discrete adapters
=================================================

Production solvers (CP, SN) need *additional* parameters that
continuous-reference solvers don't: cell counts, angular
quadrature orders, chord-integration orders. The
:class:`~orpheus.derivations.common.discretization_spec.DiscretizationSpec`
dataclass packs these into a single immutable object that
``CPMeshAdapter.from_problem`` and ``SNMeshAdapter.from_problem``
consume.

.. code-block:: python

    from orpheus.derivations.common import (
        DiscretizationSpec, GeometrySpec,
    )
    from tests.cross_method.discrete_adapters import (
        CPMeshAdapter, SNMeshAdapter,
    )

    spec = GeometrySpec(geometry="sphere",
        critical_dimension_mfp=2.5, critical_dimension_cm=5.0,
        n_groups=1)

    cp_adapter = CPMeshAdapter.from_problem(
        materials={0: pu_mixture},
        geometry_spec=spec,
        discretization=DiscretizationSpec(
            n_cells=50, n_chord_quad=64,
        ),
    )
    sol = cp_adapter.solve_critical()  # CriticalSolution

The fields of :class:`DiscretizationSpec`:

* ``n_cells`` — production mesh refinement (passed to
  :meth:`GeometrySpec.build`).
* ``n_angular`` — angular quadrature order (SN only; CP
  integrates analytically over angle).
* ``n_chord_quad`` — chord integration order (CP curvilinear
  only; slab CP is closed-form in :math:`E_3`).

Defaults match production-side conventions
(``CPParams.n_quad_y = 64``,
``GaussLegendre1D.create(n_ordinates=16)``).


Foundation-tier verification
============================

Conformance to the Protocol is a **software invariant**, not
a physics claim. The L1 chains stay where they were
(``test_fn_la13511_*.py``,
``test_peierls_greens_function_*.py``,
:mod:`tests.cross_method.test_eigenvalue`). The Protocol unifies
the *shape* those L1-verified components must expose.

Foundation-tier tests in this scaffold:

* :mod:`tests.derivations.test_transport_solver_protocol` — 8
  tests. Conformance + ``KNOWN_TRANSPORT_SOLVERS`` registry +
  attribute-surface drift detection.
* :mod:`tests.derivations.test_billiard_from_problem_unified` —
  11 tests. Bit-equality between legacy and new Billiard
  factory paths via :func:`np.array_equal` (NEVER
  :func:`np.allclose`). Asymmetric inputs catch numerical-bug
  signatures Sig 3 (scattering transpose), Sig 4 (BC swap),
  Sig 5 (factor-of-2).
* :mod:`tests.derivations.test_moment_space_from_problem_unified`
  — 6 tests. Bit-equality with the function-level F_N API.
* :mod:`tests.cross_method.test_polymorphism` — 6 tests.
  Adapter ↔ Protocol agreement on canonical cases.
* :mod:`tests.cross_method.test_discrete_adapters_smoke` — 7
  tests (1 skipped). Constructibility, non-crash solve,
  Protocol conformance for ``CPMeshAdapter`` / ``SNMeshAdapter``.

Total: 38 foundation-tier tests + 1 skipped placeholder for
the L4 cross-check.


Structural-independence preservation
=====================================

The Protocol does **not** introduce new in-house cross-method
dependencies. Per :doc:`/skills/algebra-of-record` § "Structural
independence applies above the trusted-library line":

* Continuous-reference solvers (``Billiard``, ``MomentSpace``,
  ``Spectrum``, ``BasisSpace``) consume scipy / numpy primitives
  but **do not share** any in-house module above the
  trusted-library line.
* Discrete-mesh adapters (``CPMeshAdapter``, ``SNMeshAdapter``)
  call into production CP / SN code paths which themselves do
  not share in-house modules with the continuous-reference
  family above the trusted-library line.

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
