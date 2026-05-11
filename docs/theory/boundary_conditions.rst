.. _theory-boundary-conditions:

============================================
Boundary Conditions — Trace-Law Architecture
============================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before touching anything in** :mod:`orpheus.geometry.boundary`
**or any** ``BoundaryRealizer``.

- A boundary condition is a **method-agnostic affine map** on the
  transport equation's boundary trace:
  :math:`\gamma_- \psi = R\,G\,\gamma_+ \psi + q`, where
  :math:`\gamma_\pm` are the inflow / outflow trace operators,
  :math:`G` is a geometric map (permutation, pushforward, angular
  average), :math:`R` is a response amplitude (albedo, sub-Markov
  kernel), and :math:`q` is an optional prescribed inflow source.
  See :eq:`affine-bc-form`.
- The architecture has **three concrete layers**, defined in three
  files and connected by two registries:

  +-------+-----------------------+---------------------------------------------+
  | Layer | What                  | Where                                       |
  +=======+=======================+=============================================+
  | 1     | Trace structure       | :mod:`orpheus.numerics.trace_space`         |
  |       | (Γ\_-, Γ\_+ + mask)   | (Cartesian only; curvilinear deferred)      |
  +-------+-----------------------+---------------------------------------------+
  | 2     | Boundary law          | :mod:`orpheus.geometry.boundary` (ABC +     |
  |       | (method-agnostic)     | 6 concrete laws, dual-registry: kind-keyed) |
  +-------+-----------------------+---------------------------------------------+
  | 3     | Method realizer       | :mod:`orpheus.sn.boundary_realizer`         |
  |       | (per-method strategy) | (SN functional) + 4 method stubs            |
  +-------+-----------------------+---------------------------------------------+

- Rank-N boundary conditions (Marshak, partial-current mixes) are
  expressed via Wave-0 operator algebra over realized leaves:
  ``0.3 * spec_realized + 0.7 * white_realized`` produces an
  :class:`~orpheus.numerics.operator.OperatorSum` of
  :class:`~orpheus.numerics.operator.ScaledOperator` wrappers. There
  is no dedicated ``MixedBoundaryOperator`` class (retired Wave 11).
- The :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` is the
  **only** functional realizer today. ``MoCBoundaryRealizer``,
  ``MCBoundaryRealizer``, ``CPBoundaryRealizer``, and
  ``DiffusionBoundaryRealizer`` are stubs that self-register at
  import time and raise :class:`NotImplementedError` from
  ``realize()``. The stubs hold the dispatch architecture in place
  so a future MoC / MC / CP / diffusion modernization session can
  bolt on a functional body without touching the SN side.
- The :attr:`creates_sweep_cycle` ``ClassVar`` flag on each
  :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` subclass
  signals to the SN sweep planner (§15A.2) which boundary types
  introduce cycles in the directed cell-visit graph. ``True`` on
  :class:`~orpheus.geometry.boundary.ReflectiveBoundary` and
  :class:`~orpheus.geometry.boundary.PeriodicBoundary`; ``False``
  on all other laws.
- The eight typed errors :class:`~orpheus.geometry.boundary.IncomingOutgoingTraceClassificationError`
  through :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`
  (ERR-040..ERR-047 in the V&V error catalog at
  ``.claude/skills/vv-principles/error_catalog.md``) replace the
  pre-refactor generic :class:`ValueError` raises; every one is
  pinned by a ``@pytest.mark.catches("ERR-NNN")`` decorator on the
  test that fires it.
- **Vacuum semantic correction (§16A.5).** ``VacuumInflow`` realizes
  to :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`,
  which zeros **only the inflow ordinates** and preserves the
  outflow trace. The pre-refactor ``VacuumBoundaryOperator.apply``
  used ``np.zeros_like(psi_out)`` and zeroed everything. Both are
  consistent at every existing 2-arg ``bc.apply(psi, quad)`` call
  site (all 13 read inflow rows only — audit in Wave 8 closeout),
  but the §16A.10 inflow-only mask is the trace-correct
  representation. The new behavior is the right one; the old
  behavior is preserved only on the curvilinear path until
  Issue #176 lands the curvilinear realizer.

.. admonition:: V&V status

   This page is **L4-informational** with respect to correctness.
   The architecture documented here is structural — it makes the
   code understandable and composable but does not by itself verify
   any equation. The verification load is carried by:

   - L0 foundation tests on individual primitives
     (:mod:`tests.numerics`,
     :mod:`tests.geometry.test_boundary_trace_law`,
     :mod:`tests.geometry.test_bc_errors`).
   - L1 equivalence-snapshot tests (the Wave-6 harness at
     :mod:`tests.geometry.test_bc_equivalence_snapshot`) that pin
     realizer-vs-legacy bit-equivalence (``nulp ≤ 4`` for non-vacuum
     BCs; intentional semantic divergence for vacuum per §16A.5
     above).
   - L1 universal-invariant tests
     (:mod:`tests.geometry.test_bc_universal_invariants`) that fire
     ERR-043 / ERR-044 / ERR-046 under fault-injection.

   No equation on this page makes a claim that requires a closed-form
   or MMS reference; all equations are **definitional** or
   **structural-architecture** statements drawn from Grand Report v3
   §16, §16A, §16A.10.


.. _bc-overview-three-layers:

The §16A three-layer decomposition
==================================

A boundary condition in transport-theory codes is, in the
discrete-form-typical mathematical sense, a **single linear operator**
that takes the outgoing angular flux at a face and returns the
incoming angular flux. In ORPHEUS we explicitly factor that single
operator into three layers because each layer carries different
mathematical, physical, and architectural responsibilities. The split
is taken verbatim from Grand Report v3 §16A.3 and the source-of-record
plan ``.claude/plans/transient-giggling-cake.md``.

.. _affine-bc-form:

Layer 2 — the affine law on the boundary trace
----------------------------------------------

The full mathematical form of every boundary law in this codebase is
the **affine map**

.. math::
   :label: affine-bc-form

   \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q,

where:

* :math:`\gamma_\pm` are the **trace operators** that restrict the
  angular flux :math:`\psi(\mathbf{r}, \Omega)` from the volumetric
  function space to the inflow / outflow boundary trace spaces
  :math:`\Gamma_\pm` (see :ref:`bc-trace-structure` below for the
  formal definition).
* :math:`G : \Gamma_+ \to \Gamma_+` is the **geometric map** — a
  measure-preserving permutation, pushforward, spatial wrap-around,
  or hemispheric cosine-weighted average. It carries pure geometry
  (it changes nothing about the physical interaction at the
  boundary; it just relabels the angular fluxes that meet there).
* :math:`R : \Gamma_+ \to \Gamma_-` is the **response kernel** — a
  scalar amplitude in :math:`[0, 1]` for the standard sub-Markov BCs
  (albedo, white, partial-current) or a full angular kernel in
  general weak-form BCs (deferred; see the
  :class:`~orpheus.geometry.boundary.BoundaryError` catalog and
  Issue #175 close-out follow-ups).
* :math:`q \in \Gamma_-` is the **prescribed inflow source** — a
  vector-valued quantity on :math:`\Gamma_-` only. The empty case
  :math:`q \equiv 0` is the homogeneous BC; the inhomogeneous case
  :math:`q \neq 0` is the rank-0 affine BC
  :class:`~orpheus.geometry.boundary.PrescribedInflow`.

Three remarks make this form load-bearing:

1. **Method-agnostic.** Nothing in :eq:`affine-bc-form` is SN-specific.
   The same affine map describes how MoC track-bundles, MC particle
   histories, CP boundary-to-region coupling matrices, and diffusion
   bilinear-form weak BCs all interact with the geometry. Each
   method's *realization* of the operators :math:`G`, :math:`R`,
   :math:`q` differs (see :ref:`bc-realizer-layer`); the algebraic
   shape of the law itself is universal.
2. **Affine, not linear.** The :math:`q` term is what makes the map
   affine. Most published transport-theory references treat the
   homogeneous case (:math:`q \equiv 0`) and never give the affine
   form an explicit name; ORPHEUS does because two distinct rank-0
   cases (:class:`~orpheus.geometry.boundary.VacuumInflow` with
   :math:`R = G = q = 0` and
   :class:`~orpheus.geometry.boundary.PrescribedInflow` with
   :math:`R = G = 0` but :math:`q \neq 0`) need a single uniform
   contract.
3. **The three operators are first-class.** The
   :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` ABC exposes
   :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.geometry_map`,
   :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.response_kernel`,
   and :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.source`
   as Python properties on every concrete subclass. The properties
   default to ``None / None / NoSource()``; concrete laws override
   when applicable. The split lets cross-method realizers introspect
   the law's geometric and response components separately — the SN
   realizer dispatches on the law's class today, but a future
   weak-form realizer might dispatch on the geometry / response /
   source independently.


Layer 1 — trace structure
-------------------------

The trace operators :math:`\gamma_\pm` carry their domain
information in two typed :class:`~orpheus.numerics.space.FunctionSpace`
subclasses:

* :class:`~orpheus.numerics.trace_space.InflowTraceSpace` represents
  :math:`\Gamma_- = \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^2
  : \Omega \cdot \hat n(\mathbf{r}) < 0\}` — the per-face directional
  half of the boundary on which the incoming angular flux is
  constrained by the law.
* :class:`~orpheus.numerics.trace_space.OutflowTraceSpace` represents
  :math:`\Gamma_+` symmetrically — the boundary half on which the
  outgoing flux is *not* constrained by the BC but is *consumed* by
  it (as :math:`\gamma_+ \psi`).

Both carry a **per-face directional mask** that discretizes the
sign predicate :math:`\mathrm{sign}(\Omega \cdot \hat n_f)` for each
face :math:`f` of the spatial mesh — see
:ref:`bc-trace-structure` for the geometric convention. The mask
is what the :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`
reads to build the sparse vacuum-mask operator
(:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`)
that zeros precisely the inflow ordinates at one face.


.. _bc-realizer-layer:

Layer 3 — the method realizer
-----------------------------

A single :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
describes the physics at the boundary but is **not** by itself
ready for consumption by a transport sweep. The conversion from
method-agnostic law to method-specific
:class:`~orpheus.numerics.operator.LinearOperator` is the job of a
:class:`~orpheus.geometry.boundary.BoundaryRealizer`. For each
transport method (``"SN"``, ``"MoC"``, ``"MC"``, ``"CP"``,
``"diffusion"``) there is one realizer class, registered in the
:class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry` under
the method name.

The realizer's :meth:`realize` method takes the law plus a
**method space** — a method-specific container holding the
quadrature, mesh, trace masks, and any other discretization
metadata the realizer needs — and returns a 1-arg
:class:`~orpheus.numerics.operator.LinearOperator` whose
:meth:`apply` carries the method-specific realization of the
affine BC :eq:`affine-bc-form`.

Why this third layer? Because the same affine law is realized by
**structurally different** linear operators in each transport
method:

* SN realizes vacuum as a sparse per-ordinate
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor` on
  the inflow ordinates of the affected face.
* MoC realizes vacuum by zeroing the entering boundary fluxes of
  every track that intersects the face.
* MC realizes vacuum by killing particle histories at the face.
* CP realizes vacuum as zero rows in the boundary-to-region
  coupling matrix.

Splitting the realizer out of the law makes each piece independently
testable and gives every method a single bolt-in point for its BC
treatment — see :ref:`bc-cross-method-stubs` for the four stubs
that hold the architecture today.


.. _bc-trace-structure:

Trace structure (Γ\_-, Γ\_+)
============================

The transport equation lives on a phase space :math:`\Omega \times
S^d` where :math:`\Omega \subset \mathbb{R}^d` is the spatial domain
and :math:`S^d` is the unit sphere of directions. The boundary
:math:`\partial\Omega` carries an outward unit normal :math:`\hat
n(\mathbf{r})` at every regular point. For an angular flux
:math:`\psi(\mathbf{r}, \Omega)` defined on the full phase space, the
**boundary trace** splits naturally into two pieces by the sign of
:math:`\Omega \cdot \hat n`:

.. math::
   :label: trace-sign-predicate

   \Gamma_- \;=\; \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^d
                  : \Omega \cdot \hat n(\mathbf{r}) < 0\},
   \qquad
   \Gamma_+ \;=\; \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^d
                  : \Omega \cdot \hat n(\mathbf{r}) > 0\}.

Points with :math:`\Omega \cdot \hat n = 0` are **tangential** —
they belong to neither half. For axis-aligned ordinates on
axis-aligned faces these arise exactly (no round-off) for face
normals perpendicular to the ordinate's direction cosine; for
general curvilinear faces or generic ordinates they arise only
at a measure-zero subset that the discrete representation
identifies via a small tolerance (``_TANGENTIAL_EPS = 1e-12``).

In the discrete setting, the spatial boundary is a union of finite
faces :math:`\{f_1, \ldots, f_F\}` and the angular variable is a
finite ordinate set :math:`\{\Omega_n : n = 1, \ldots, N\}`. The
sign predicate :eq:`trace-sign-predicate` then collapses to a
**per-face boolean mask** of shape :math:`(F, N)`:

.. math::
   :label: inflow-mask-discrete

   \mathrm{inflow\_mask}[f, n]
   \;=\; \bigl(\Omega_n \cdot \hat n_f < -\epsilon\bigr),
   \qquad
   \mathrm{outflow\_mask}[f, n]
   \;=\; \bigl(\Omega_n \cdot \hat n_f > +\epsilon\bigr).

This mask is the discrete realization of :math:`\Gamma_\pm`. It is
the load-bearing primitive that downstream consumers need:

* The SN realizer's vacuum branch reads
  ``inflow_mask[f]`` for the specific face :math:`f` and converts
  it to an integer array of ordinate indices via
  :meth:`~orpheus.numerics.trace_space.InflowTraceSpace.inflow_indices_for_face`.
  Those indices are the constructor argument to
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`.
* The universal invariant
  :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
  uses the inflow mask to validate that a
  :class:`~orpheus.geometry.boundary.BoundarySource` has no nonzero
  entries on outflow ordinates (per
  :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`,
  ERR-047).
* Future curvilinear Krylov solvers (deferred — see Issue #176)
  will read the inflow mask to construct the boundary-condition
  projection matrix-free, instead of routing through the legacy
  2-arg ``law.apply(psi, quad)`` shim.

The class :class:`~orpheus.numerics.trace_space.InflowTraceSpace`
carries the mask as an ``Optional[np.ndarray]`` field excluded from
:meth:`__eq__` and :meth:`__hash__` — preserving the
:class:`~orpheus.numerics.space.FunctionSpace` identity convention
``(name, shape)``. Construction goes through the classmethod
factory :meth:`InflowTraceSpace.from_mesh_and_quadrature`; the bare
dataclass constructor is reserved for trace spaces whose mask will
be populated later (or never).

**Curvilinear support is deferred.** The factory raises
:class:`NotImplementedError` for ``Mesh1D`` and ``Mesh2D`` with
``coord ∈ {SPHERICAL, CYLINDRICAL}`` because no curvilinear-Krylov
consumer needs the mask yet. The deferral is grep-able for the
future consumer; see Issue #176 and the
"BC: curvilinear InflowTraceSpace.from_mesh_and_quadrature
implementation" follow-up tracked under ``module:geometry``.


.. _bc-law-layer:

Boundary law (``BoundaryTraceLaw`` ABC + concretes)
===================================================

The base class :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
is an ``abc.ABC`` that combines
:class:`~orpheus.numerics.operator.LinearOperatorMixin` (for
operator-algebra dunders like ``+``, ``*``, ``@``) and
:class:`~orpheus.numerics.registry.RegistryMixin` (so each concrete
subclass self-registers under its ``key=`` class-creation kwarg).
The ABC ships:

1. Three first-class properties carrying the
   :eq:`affine-bc-form` operators: ``geometry_map``,
   ``response_kernel``, ``source`` (defaulting to ``None``, ``None``,
   :class:`~orpheus.geometry.boundary.NoSource`).
2. The :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.creates_sweep_cycle`
   ``ClassVar`` boolean signal — see :ref:`bc-sweep-cycle`.
3. Five **universal** ``assert_*`` invariants (no-op defaults;
   concrete laws override the relevant ones) and four **specific**
   invariants on the BCs that need them — see
   :ref:`bc-universal-invariants`.
4. A :meth:`realize` hook that defers to the
   :class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry` —
   see :ref:`bc-realizer-layer-detail`.
5. An abstract :meth:`apply` that accepts
   ``(psi_out, *args, **kwargs)`` — transitional contract
   accommodating both the legacy 2-arg
   ``apply(psi_out, quadrature)`` from the curvilinear path and
   the 1-arg ``apply(psi)`` that resolved realizer output uses.
   The two paths converge once Issue #176 ships the curvilinear
   realizer.

The six concrete laws ship under :mod:`orpheus.geometry.boundary`,
one per submodule. The Grand Report v3 vocabulary is used verbatim
for the class names; pre-refactor names are retained as deprecated
aliases in the package ``__init__`` (see :ref:`bc-naming-audit`).

.. list-table:: Concrete ``BoundaryTraceLaw`` subclasses
   :header-rows: 1
   :widths: 18 16 13 13 13 27

   * - Class
     - Registry key
     - :math:`G_\alpha`
     - :math:`R_\alpha`
     - :math:`q`
     - Sweep-cycle flag
   * - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - ``"vacuum"``
     - rank-0 (none)
     - 0
     - 0
     - ``False``
   * - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - ``"reflective"``
     - axis-reflection permutation
     - albedo
     - 0
     - **``True``** (signals §15A.2)
   * - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - ``"white"``
     - cosine-weighted hemispheric average (Lambertian)
     - albedo
     - 0
     - ``False``
   * - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - ``"periodic"``
     - spatial pushforward (caller-supplied)
     - 1
     - 0
     - **``True``** (signals §15A.2)
   * - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - ``"albedo"``
     - identity in angle
     - albedo
     - 0
     - ``False``
   * - :class:`~orpheus.geometry.boundary.PrescribedInflow`
     - ``"prescribed_inflow"``
     - 0
     - 0
     - :math:`q \in \Gamma_-`
     - ``False``

Two are rank-1 in space (periodic), four are rank-1 in angle
(reflective, white, albedo) and write the same incoming flux to
every inflow ordinate; one is rank-0 (vacuum / prescribed inflow);
and Marshak / partial-current boundaries are rank-N via Wave-0
algebra over realized leaves — see :ref:`bc-rank-n-algebra` below.


.. _bc-naming-audit:

Naming audit: pre-refactor vs Grand Report v3 vocabulary
--------------------------------------------------------

Wave 7 of the refactor renamed every concrete BC to match the
Grand Report v3 vocabulary verbatim. The pre-refactor names are
preserved as deprecated re-exports in
:mod:`orpheus.geometry.boundary.__init__` and the package's
public ``__all__`` so existing import sites keep working unchanged
during the deprecation window:

.. list-table:: Wave 7 BC renames
   :header-rows: 1
   :widths: 35 35 30

   * - Pre-refactor name
     - Wave 7 canonical name
     - Why renamed
   * - ``VacuumBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - Emphasizes "inflow set to zero", not "operator that vacuums";
       distinguishes from the rank-N case
       :class:`PrescribedInflow` which also writes only the inflow
       trace.
   * - ``SpecularBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - "Specular" is one specific axis-aligned reflection;
       "Reflective" is the family name that the Grand Report uses.
       A future ``SymmetryBoundary`` (deferred) will share the
       reflective-family base but apply on non-physical octant
       boundaries.
   * - ``WhiteBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - Drops the redundant "Operator" suffix that pre-dated the
       law / realizer split. The law is no longer "the operator";
       it's the abstract physical statement that gets *realized*
       to an operator.
   * - ``PeriodicBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - Same rationale: pre-refactor "Operator" suffix is
       structurally misleading.
   * - ``AlbedoBoundaryOperator``
     - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - Same rationale.
   * - ``MixedBoundaryOperator``
     - **retired Wave 11** — see :ref:`bc-rank-n-algebra`.
     - Replaced by Wave-0 ``OperatorSum`` algebra over realized
       leaves; the dedicated class added no value over the
       inherited dunders.


.. _bc-sweep-cycle:

The ``creates_sweep_cycle`` signal
----------------------------------

The SN sweep is a topological sort of the directed cell-visit graph
where edges are oriented by :math:`\mathrm{sign}(\Omega_n \cdot
\hat n_f)`. For most BCs the boundary is the *root* of the sort —
inflow values come from the BC, get propagated through the cells,
and exit as outflow values that the BC consumes but doesn't feed
back. For two BC families this is no longer true:

* **Reflective.** The outflow flux at a face is mapped to an inflow
  flux at the **same** face (under the reflection permutation). The
  cell visits that produce the outflow are predecessors of the cell
  visits that consume the reflected inflow — but those are *the
  same cells*. The sweep DAG acquires a cycle and a self-consistent
  fixed-point iteration must converge the reflective inflow against
  the outflow.
* **Periodic.** Same situation, but the cycle spans two different
  faces (outflow at face A maps to inflow at face B, and vice
  versa).

The :attr:`~orpheus.geometry.boundary.BoundaryTraceLaw.creates_sweep_cycle`
``ClassVar`` is the structural signal that lets the §15A.2
sweep-cycle detector identify these BCs without inspecting the
realization. It is ``True`` on :class:`ReflectiveBoundary` and
:class:`PeriodicBoundary`, ``False`` on every other law.

Vacuum, white, albedo, and prescribed-inflow are all cycle-free:
they consume :math:`\gamma_+ \psi` and produce a function of it
that only depends on :math:`\gamma_+` (or, for prescribed inflow,
nothing at all). No fixed point is needed.


.. _bc-realizer-layer-detail:

Realizer (``BoundaryRealizer`` Protocol + ``SNBoundaryRealizer``)
=================================================================

The :class:`~orpheus.geometry.boundary.BoundaryRealizer` Protocol is
``@runtime_checkable`` and lives at
:mod:`orpheus.geometry.boundary._realizer`. Its contract is one
attribute and one method:

.. code-block:: python

   @runtime_checkable
   class BoundaryRealizer(Protocol):
       method_name: str
       def realize(
           self,
           law: BoundaryTraceLaw,
           method_space: Any,
       ) -> LinearOperator: ...

The Protocol intentionally does *not* prescribe how
:meth:`realize` dispatches over law types — different methods will
have different optimal dispatch strategies. The SN realizer uses
``isinstance`` because the law set is small and stable; a future
realizer that needs runtime extension might use the
:class:`~orpheus.numerics.registry.RegistryMixin` machinery instead.

The Wave 5 SN dispatch table is the documented standard:

.. list-table:: SN realization map (law → Wave-0 / Wave-1 primitive)
   :header-rows: 1
   :widths: 24 38 38

   * - Law
     - Realized representation (α = 1 fast path)
     - Realized representation (α ∉ {0, 1})
   * - :class:`VacuumInflow`
     - :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
       with the per-face ``inflow_indices`` from the method space's
       :class:`~orpheus.numerics.trace_space.InflowTraceSpace`.
     - n/a (vacuum has no α parameter)
   * - :class:`ReflectiveBoundary(axis, α)`
     - bare
       :class:`~orpheus.numerics.operator.PermutationOperator`
       with ``perm = quadrature.reflection_index(axis)``
     - ``ScaledOperator(α, PermutationOperator(...))``
   * - :class:`WhiteBoundary(axis, outward_sign, α)`
     - bare
       :class:`~orpheus.sn.angular_operator.AngularAverageOperator.from_quadrature(quad, axis, outward_sign)`
     - ``ScaledOperator(α, AngularAverageOperator(...))``
   * - :class:`AlbedoBoundary(α)` with α=0
     - :class:`~orpheus.numerics.operator.ZeroOperator`
     -
   * - :class:`AlbedoBoundary(α)` with α=1
     - :class:`~orpheus.numerics.operator.IdentityOperator`
     -
   * - :class:`AlbedoBoundary(α)` with α ∉ {0, 1}
     -
     - ``ScaledOperator(α, IdentityOperator())``
   * - :class:`PeriodicBoundary`
     - :class:`~orpheus.numerics.operator.PeriodicWrapOperator`
       (today an angular identity; spatial-pushforward extension
       pending — see "BC: PeriodicWrapOperator spatial-pushforward
       implementation" follow-up, ``module:sn``)
     - n/a (periodic has no α parameter)
   * - :class:`PrescribedInflow(source)`
     - :class:`~orpheus.sn.angular_operator.IncomingSourceOperator(source)`
       — :meth:`apply` ignores the outgoing flux and returns
       ``source.evaluate(probe_inflow_trace)``
     - n/a

The α = 1.0 fast paths return the **bare** primitive (no
``ScaledOperator`` wrap). This is load-bearing for bit-identity:
without it, the legacy "perfect reflection" case
``SpecularBoundaryOperator(axis="x", albedo=1.0)`` would shift by
one ULP under the realizer relative to its pre-refactor
``np.take(psi_out, reflection_index, axis=0)`` body — see the
Wave 6 snapshot harness for the bit-equivalence pin.

The :class:`~orpheus.sn.method_space.SNMethodSpace` dataclass is the
realizer's second argument. It carries:

* :attr:`~orpheus.sn.method_space.SNMethodSpace.quadrature` — the
  angular quadrature (mandatory).
* :attr:`~orpheus.sn.method_space.SNMethodSpace.face` — the face
  label (``"left"``, ``"xmin"``, etc.) so the vacuum branch can
  look up the right inflow indices.
* :attr:`~orpheus.sn.method_space.SNMethodSpace.inflow_indices` —
  the per-face inflow ordinate indices for the vacuum branch.
* :attr:`~orpheus.sn.method_space.SNMethodSpace.mesh`,
  :attr:`~orpheus.sn.method_space.SNMethodSpace.inflow_trace`,
  :attr:`~orpheus.sn.method_space.SNMethodSpace.outflow_trace` —
  full mesh + trace metadata for any future realizer branch that
  needs more than the per-face index list.

The :meth:`SNMethodSpace.for_face` factory is the standard
construction site at
:meth:`orpheus.sn.geometry.SNMesh._resolve_bcs` time; the
:meth:`SNMethodSpace.minimal` factory returns a quadrature-only
method space for unit tests that don't need mesh + face metadata.


.. _bc-dual-registry:

Dual registry (``BoundaryLawRegistry`` + ``BoundaryRealizerRegistry``)
======================================================================

Two registries serve two independent lookup keys:

1. **Law registry** — keyed by ``BC.kind`` string
   (``"vacuum"``, ``"reflective"``, ``"white"``, ``"periodic"``,
   ``"albedo"``, ``"prescribed_inflow"``). The registry IS
   :attr:`BoundaryTraceLaw.registry` (a class-level dict
   maintained by :class:`~orpheus.numerics.registry.RegistryMixin`).
   Concrete laws self-register at module import time via the
   ``key=`` class-creation kwarg:

   .. code-block:: python

      class VacuumInflow(BoundaryTraceLaw, key="vacuum"):
          ...

   Lookup is :meth:`BoundaryTraceLaw.create("vacuum")` or direct
   dictionary access ``BoundaryTraceLaw.registry["vacuum"]``.
   :meth:`orpheus.sn.geometry.SNMesh._resolve_one` uses the
   second form to recover the law class from a mesh-declared
   :class:`~orpheus.geometry.mesh.BC`.

2. **Realizer registry** — keyed by method name (``"SN"``,
   ``"MoC"``, ``"MC"``, ``"CP"``, ``"diffusion"``). The registry is
   :class:`~orpheus.geometry.boundary.BoundaryRealizerRegistry`, a
   stand-alone class with a class-level ``_registry`` dict. The
   :meth:`~orpheus.geometry.boundary.BoundaryRealizerRegistry.register`
   decorator and :meth:`~orpheus.geometry.boundary.BoundaryRealizerRegistry.get`
   classmethod are the canonical add / lookup surface.

The two registries are **disjoint by design** (§16A.11 lines
3252–3257 of the Grand Report). They could not collide on a
shared key — laws are keyed by physical-name strings, realizers by
transport-method-name strings — but more importantly they describe
two orthogonal extension axes:

* Adding a new BC type means adding one
  :class:`BoundaryTraceLaw` subclass and registering it under one
  ``BC.kind`` string. **N** new realizer cases need to be added
  (one per method), but no existing realizer needs to change.
* Adding a new transport method means adding one
  :class:`BoundaryRealizer` and registering it under one method
  name. **M** new realizer branches need to be implemented (one
  per existing BC), but no existing law needs to change.

The realizer registry uses a **stand-alone** registry class instead
of :class:`RegistryMixin` because realizers are independent
strategies keyed by string — not a class hierarchy with a shared
base. Collisions at
:meth:`BoundaryRealizerRegistry.register` raise
:class:`~orpheus.geometry.boundary.BoundaryRealizerRegistryError`
at import time rather than silently overriding; this surfaces
duplicate-registration bugs during package import.


.. _bc-cross-method-stubs:

Cross-method realizer stubs
---------------------------

.. note::

   The cross-method realizers are **scaffolding for future
   modernization**, not a mature shared abstraction. The
   :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` is the
   single functional realizer today. ``MoCBoundaryRealizer``,
   ``MCBoundaryRealizer``, ``CPBoundaryRealizer``, and
   ``DiffusionBoundaryRealizer`` are stubs that hold the dispatch
   table in place — they self-register at module import time so
   :meth:`BoundaryRealizerRegistry.get("MoC")` returns the class,
   but the class's :meth:`realize` raises :class:`NotImplementedError`
   with a grep-able message pointing at the follow-up issue.

   When the second functional realizer ships (whichever of MoC / MC
   / CP / diffusion adopts the unified BC architecture first), the
   shared abstraction that the realizers form will become
   discoverable; until then, building a ``BoundaryRealizerBase``
   ABC would violate the "Unify after two instances" architectural
   discipline (Cardinal Rule 2 generalization). Each stub is its
   own ~30-line class with explicit ``method_name`` and ``realize``
   method.

The four stub files are:

* :mod:`orpheus.moc.boundary_realizer` — ``MoCBoundaryRealizer``.
  Follow-up: "BC: MoCBoundaryRealizer functional implementation"
  (label ``module:moc, type:feature``).
* :mod:`orpheus.mc.boundary_realizer` — ``MCBoundaryRealizer``.
  Follow-up: "BC: MCBoundaryRealizer functional implementation"
  (label ``module:mc, type:feature``).
* :mod:`orpheus.cp.boundary_realizer` — ``CPBoundaryRealizer``.
  Follow-up: "BC: CPBoundaryRealizer functional implementation"
  (label ``module:cp, type:feature``).
* :mod:`orpheus.diffusion.boundary_realizer` —
  ``DiffusionBoundaryRealizer``. Follow-up:
  "BC: DiffusionBoundaryRealizer functional implementation"
  (label ``module:diffusion, type:feature``).

Each method's ``__init__.py`` imports the boundary-realizer
submodule unconditionally so plain ``import orpheus.<method>``
auto-registers the stub. The SN realizer is **not**
auto-imported by ``orpheus.sn.__init__`` (it's a heavy module
that every SN consumer pays for); the
:class:`~orpheus.sn.geometry.SNMesh` constructor imports it
explicitly when it needs it.


.. _bc-worked-example:

Worked example — end to end
===========================

The following walks the
``BC("vacuum") → VacuumInflow → SNBoundaryRealizer.realize →
IncomingOrdinateMaskTensor`` chain that
:meth:`orpheus.sn.geometry.SNMesh._resolve_bcs` performs at SNMesh
construction time. The example uses a 1-D Cartesian slab; the same
chain runs on Mesh2D with face labels ``xmin`` / ``xmax`` /
``ymin`` / ``ymax``.

Step 1 — declaration on the mesh
--------------------------------

The user declares the vacuum BC on the mesh's left face:

.. code-block:: python

   from orpheus.geometry.mesh import Mesh1D, BC

   mesh = Mesh1D(
       edges=np.linspace(0.0, 1.0, 11),
       mat_ids=np.zeros(10, dtype=int),
       coord=CoordSystem.CARTESIAN,
       bc_left=BC("vacuum"),
       bc_right=BC("reflective"),
   )

The :class:`~orpheus.geometry.mesh.BC` dataclass is a thin wrapper
``BC(kind: str, params: dict)`` with no SN-specific knowledge. The
mesh is method-agnostic.

Step 2 — law resolution (in ``SNMesh.__init__``)
------------------------------------------------

When :class:`~orpheus.sn.geometry.SNMesh` is constructed against the
mesh, its :meth:`_resolve_bcs` method walks the four (1-D: two)
faces:

.. code-block:: python

   law_cls = SNMesh.BOUNDARY_OPERATOR_REGISTRY["vacuum"]
   # law_cls is VacuumInflow  (registry key -> law class lookup)
   law = law_cls()
   # law is a zero-arg instance: VacuumInflow has no parameters

The :attr:`SNMesh.BOUNDARY_OPERATOR_REGISTRY` is the SN-side view
of the law registry; today it carries only ``"vacuum"`` and
``"reflective"`` because those are the only kinds the SN sweep
pipeline has been wired for in production (the other three —
white, periodic, albedo — are realizable but require sweep-side
plumbing tracked in separate issues). Adding a new kind is one
dict-entry edit.

Step 3 — method space construction
----------------------------------

For the **Cartesian** path, the SNMesh has precomputed the
``InflowTraceSpace`` once for the whole mesh at construction time:

.. code-block:: python

   from orpheus.numerics.trace_space import InflowTraceSpace

   self._inflow_trace = InflowTraceSpace.from_mesh_and_quadrature(
       mesh, self.quad, faces=("left", "right"),
   )

The per-face inflow mask is a ``(2, N)`` boolean array indexed by
``(face_idx, ordinate_idx)``. For the left face (``axis=0,
sign=-1``), the inflow predicate is ``-mu_x[n] < 0``, i.e.
``mu_x[n] > 0`` — the rightward-pointing ordinates are inflow at
the left face.

The :meth:`SNMethodSpace.for_face` factory carves out the per-face
slice:

.. code-block:: python

   from orpheus.sn.method_space import SNMethodSpace

   method_space = SNMethodSpace.for_face(
       mesh=self.mesh,
       quadrature=self.quad,
       face="left",
       inflow_trace=self._inflow_trace,
       outflow_trace=self._outflow_trace,
   )
   # method_space.inflow_indices is a 1-D int array:
   # [n for n in range(N) if mu_x[n] > 1e-12]

Step 4 — realization
--------------------

The :class:`SNBoundaryRealizer` is now invoked. Looking up by
method name produces the realizer class; instantiation is
stateless:

.. code-block:: python

   from orpheus.geometry.boundary import BoundaryRealizerRegistry

   realizer_cls = BoundaryRealizerRegistry.get("SN")
   # realizer_cls is SNBoundaryRealizer
   realizer = realizer_cls()
   realized = realizer.realize(law, method_space)

The realizer's vacuum branch fires:

.. code-block:: python

   # Inside SNBoundaryRealizer.realize:
   if isinstance(law, VacuumInflow):
       return IncomingOrdinateMaskTensor(
           inflow_indices=method_space.inflow_indices,
           n_ordinates=quad.N,
           axis=0,
       )

The returned ``realized`` is a
:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`,
which is a self-adjoint, idempotent Wave-0 primitive with
``capabilities = frozenset({CAP_APPLY, CAP_APPLY_TRANSPOSE})``.

Step 5 — bind for the legacy contract
-------------------------------------

Until Issue #176 lands the curvilinear realizer, the
``SNMesh.bc_*`` attributes must respond to both the 1-arg
``bc.apply(psi)`` contract (post-Wave-9 production call sites) and
the 2-arg legacy ``bc.apply(psi, quad)`` contract (a small number
of test sites and the curvilinear sweep where the realizer cannot
yet ship). The
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim is the bridge:

.. code-block:: python

   from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator

   return _BoundBoundaryOperator(realized, kind=bc.kind)
   # In curvilinear (no inflow_trace) the shim wraps the bare law
   # with a bound quadrature so the 2-arg apply body still runs:
   # _BoundBoundaryOperator(law, quadrature=self.quad, kind=bc.kind)

The shim is **internal** to the package (not in
:attr:`__all__`) — a test pins its private status. Issue #176
will delete the shim once the curvilinear realizer ships and every
consumer routes through the 1-arg contract.

Step 6 — consumption by the sweep
---------------------------------

At each sweep call site the resolved operator is applied with
the uniform 1-arg interface:

.. code-block:: python

   psi_in = self.bc_left.apply(psi_out_full)
   # Shim forwards to IncomingOrdinateMaskTensor.apply, which zeros
   # only the inflow ordinates; outflow rows pass through unchanged.

The downstream sweep reads ``psi_in[inflow_ord]`` only — every
production call site was audited in Wave 8 (see the §16A.5
compatibility-audit table in the Wave 8 close-out memo and the
:ref:`bc-vacuum-semantic-correction` section below).


.. _bc-universal-invariants:

Universal ``assert_*`` invariants
=================================

The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` ABC ships
**five universal** assertion methods (no-op defaults; concrete laws
override) plus **four specific** assertions on the BCs that need
them. Together they form the structural verification surface that
the :mod:`tests.geometry.test_bc_universal_invariants` suite
exercises (~30 L1 tests).

The split between "universal" and "specific" follows the
Grand Report v3 §16A.12 + §27.6 catalog. Universal invariants are
properties the affine law :eq:`affine-bc-form` should satisfy for
**any** physically meaningful BC; specific invariants pin properties
that only a subset of laws claim (e.g. involution is meaningful for
reflective but not for white).

.. list-table:: Universal invariants
   :header-rows: 1
   :widths: 30 25 45

   * - Method
     - Pinned error
     - What it asserts
   * - ``assert_inflow_outflow_classification``
     - :class:`~orpheus.geometry.boundary.IncomingOutgoingTraceClassificationError`
       (ERR-040)
     - Every ordinate at the face is either inflow or outflow (no
       tangential ordinates allowed by the law's contract). Default:
       no-op; reflective overrides to require strict partition.
   * - ``assert_outgoing_leakage_unconstrained``
     - n/a (architectural contract)
     - The outgoing trace flux is not constrained by the BC.
       Default: no-op. Future Dirichlet-outflow / prescribed
       cell-edge interface laws would override.
   * - ``assert_geometry_map_measure_preserving``
     - :class:`~orpheus.geometry.boundary.BoundaryGeometryMapNotMeasurePreservingError`
       (ERR-042)
     - The geometric map :math:`G` preserves the angular measure
       :math:`w(\Omega)\,|\Omega \cdot \hat n|`. Default: no-op.
       Reflective overrides (delegates to the involution check).
   * - ``assert_response_positive_if_declared``
     - :class:`~orpheus.geometry.boundary.BoundaryResponseNotPositiveError`
       (ERR-043)
     - If a response kernel is declared, it produces non-negative
       output on the inflow trace. Default: no-op. Albedo / white
       override.
   * - ``assert_source_lives_on_incoming_trace``
     - :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`
       (ERR-047)
     - The source :math:`q` is nonzero only on :math:`\Gamma_-`.
       Default: no-op (the :class:`NoSource` default trivially
       satisfies). Prescribed-inflow + future user-source classes
       override.

.. list-table:: Specific invariants
   :header-rows: 1
   :widths: 40 25 35

   * - Method
     - Pinned error
     - Where it's defined
   * - :meth:`ReflectiveBoundary.assert_is_involutive`
     - :class:`~orpheus.geometry.boundary.ReflectionNotInvolutiveError`
       (ERR-044)
     - :class:`ReflectiveBoundary` (reflection-index table must
       satisfy ``perm[perm] == arange``).
   * - :meth:`ReflectiveBoundary.assert_maps_inflow_to_outflow`
     - :class:`~orpheus.geometry.boundary.ReflectionDidNotMapInflowToOutflowError`
       (ERR-045)
     - :class:`ReflectiveBoundary` (every inflow ordinate maps to an
       outflow ordinate under the reflection).
   * - :meth:`ReflectiveBoundary.assert_direction_norm_preserved`
     - inherits ERR-042 via measure-preservation
     - :class:`ReflectiveBoundary` (delegates to involution check).
   * - :meth:`WhiteBoundary.assert_submarkov` /
       :meth:`AlbedoBoundary.assert_submarkov`
     - :class:`~orpheus.geometry.boundary.SubmarkovViolationError`
       (ERR-046)
     - :class:`WhiteBoundary` + :class:`AlbedoBoundary` (the sub-Markov
       kernel constraint :math:`\int R\,\mathrm{d}y \leq 1` per row;
       albedo :math:`> 1` violates this physically).


.. _bc-named-error-catalog:

Named-error catalog (ERR-040..ERR-047)
======================================

Per Grand Report v3 §26A.4 and the `vv-principles` skill's
"Log every caught bug" directive, every typed error has a matching
``@pytest.mark.catches("ERR-NNN")`` decorator on the test that
proves it fires under the right fault-injection. The eight errors
shipped under :mod:`orpheus.geometry.boundary._errors` are:

.. list-table::
   :header-rows: 1
   :widths: 18 26 12 44

   * - Error class
     - Trigger condition
     - Mode
     - Mechanism
   * - :class:`~orpheus.geometry.boundary.IncomingOutgoingTraceClassificationError`
       (ERR-040)
     - Tangential ordinate at a face that requires strict partition.
     - #5 (index)
     - ``assert_inflow_outflow_classification`` finds
       ``|Ω · n| ≤ ε`` on a face where the law's contract
       forbids it.
   * - :class:`~orpheus.geometry.boundary.VacuumAppliedToOutgoingTraceError`
       (ERR-041)
     - Vacuum law applied to an outgoing trace.
     - #6 (convention)
     - Vacuum sets only :math:`\gamma_- \psi = 0`; applying it on
       :math:`\Gamma_+` is geometrically meaningless and typically
       indicates a wrong face annotation.
   * - :class:`~orpheus.geometry.boundary.BoundaryGeometryMapNotMeasurePreservingError`
       (ERR-042)
     - Geometric map :math:`G` does not preserve
       :math:`w(\Omega) |\Omega \cdot \hat n|`.
     - #5 + #6
     - Wrong reflection-index table or inconsistent quadrature
       :math:`\mu_n` / weights.
   * - :class:`~orpheus.geometry.boundary.BoundaryResponseNotPositiveError`
       (ERR-043)
     - Response kernel produces negative output.
     - #1 (sign)
     - Sign-flipped kernel construction (e.g. accidental
       :math:`-\alpha`).
   * - :class:`~orpheus.geometry.boundary.ReflectionNotInvolutiveError`
       (ERR-044)
     - Reflection permutation is not its own inverse:
       :math:`\pi \circ \pi \neq \mathrm{id}`.
     - #5 (index)
     - Wrong reflection axis or partial permutation in the
       :meth:`quadrature.reflection_index` table.
   * - :class:`~orpheus.geometry.boundary.ReflectionDidNotMapInflowToOutflowError`
       (ERR-045)
     - Reflection maps an inflow ordinate to itself.
     - #5 (index)
     - Non-axis-aligned reflection mislabeled as ``ReflectiveBoundary``
       (the right BC family would be a future ``SymmetryBoundary``).
   * - :class:`~orpheus.geometry.boundary.SubmarkovViolationError`
       (ERR-046)
     - Sub-Markov BC with :math:`\alpha > 1`.
     - #4 (factor)
     - Albedo / white kernel scalar exceeds 1.0 — physically this
       would imply a source on the boundary surface.
   * - :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`
       (ERR-047)
     - Boundary source :math:`q` has nonzero outflow entries.
     - #6 (convention)
     - User-supplied :class:`BoundarySource` has nonzero entries on
       :math:`\Gamma_+`; geometrically meaningless and indicates a
       wrong source-shape contract.

All eight extend :class:`ValueError` (via the
:class:`~orpheus.geometry.boundary.BoundaryError` base) so existing
``except ValueError`` consumers from the pre-refactor code keep
working. Every catch site can additionally pattern-match on the
typed subclass to recover the offending law name from
:attr:`BoundaryError.law`.


.. _bc-rank-n-algebra:

Wave-0 algebra for rank-N boundaries
====================================

Rank-N (Marshak, partial-current) boundary conditions are
**not** a special class. They are expressed via the
:class:`~orpheus.numerics.operator.OperatorSum` /
:class:`~orpheus.numerics.operator.ScaledOperator` algebra
that :class:`~orpheus.numerics.operator.LinearOperatorMixin`
provides to every :class:`BoundaryTraceLaw` subclass through
operator dunders (``+``, ``-``, ``*``).

The standard Marshak boundary (Bell & Glasstone 1970 §1.5) — a
mix of specular reflection (weight :math:`c_1`) and diffuse
white reflection (weight :math:`c_2`) — is:

.. code-block:: python

   from orpheus.geometry.boundary import (
       ReflectiveBoundary, WhiteBoundary,
   )
   from orpheus.sn.boundary_realizer import (
       SNBoundaryRealizer, SNMethodSpace,
   )

   # Realize each leaf against the SN method space.
   ms = SNMethodSpace.for_face(...)
   spec_realized = SNBoundaryRealizer().realize(
       ReflectiveBoundary(axis="x"), ms,
   )
   white_realized = SNBoundaryRealizer().realize(
       WhiteBoundary(axis="x", outward_sign=+1), ms,
   )

   # Wave-0 algebra over realized leaves.
   marshak = 0.3 * spec_realized + 0.7 * white_realized
   # marshak is OperatorSum(
   #     ScaledOperator(0.3, PermutationOperator(...)),
   #     ScaledOperator(0.7, AngularAverageOperator(...))
   # )

The result is a Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator`-wrapped Wave-0
primitives, consumable directly by the SN sweep / Krylov path
via its 1-arg :meth:`apply`.

The pre-refactor implementation
---------------------------------

The pre-Wave-11 code carried a dedicated
``MixedBoundaryOperator(components: list[tuple[float,
BoundaryOperator]])`` class whose :meth:`apply` body looped over
``components`` and summed ``coeff * primitive.apply(psi, quad)``.
The SN realizer dispatched on it via an
``isinstance(law, MixedBoundaryOperator)`` branch that ran the same
loop with ``coeff * realize(primitive, ms)`` summed via
:class:`OperatorSum`.

Wave 11 deleted this code because:

1. **Closure-by-inheritance.** Every
   :class:`BoundaryTraceLaw` already inherits ``+`` /
   ``__rmul__`` from
   :class:`~orpheus.numerics.operator.LinearOperatorMixin`. Writing
   ``0.3 * spec + 0.7 * white`` directly produces the same
   :class:`OperatorSum` shape — no extra class needed.
2. **The class added no semantic value.** The pre-refactor
   ``MixedBoundaryOperator`` was a *container* that delayed the
   realization until ``apply`` time. After the realizer split
   that delay no longer fit: every BC must be realized *before*
   the SN sweep starts (the resolved operator carries the
   per-face inflow indices, which the bare law has no access to).
3. **Structural simplicity.** The pre-refactor mixed-BC class
   carried its own ``apply`` body, ``apply_transpose`` body,
   registry entry, and the SN realizer's special branch. Wave 11
   deleted ~340 lines (the class + the dispatch branch + four
   class-specific tests) and added ~210 lines for the
   :func:`~orpheus.sn.boundary_realize.realize_recursively`
   walker that handles the **general** case of any Wave-0
   composer tree (not just rank-N sums).

The Wave 6 snapshot harness verified that the bit-identity of the
Marshak case (``0.3 * specular + 0.7 * white``) is **preserved** by
the Wave 11 transition: the realized
``OperatorSum(ScaledOperator, ScaledOperator)`` produces the same
floating-point output as the deleted ``MixedBoundaryOperator.apply``
because the reduction tree is identical.


.. _bc-realize-recursively:

The ``realize_recursively`` walker
==================================

Users may want to build a boundary law as an expression tree:

.. code-block:: python

   law = (
       0.3 * ReflectiveBoundary(axis="x")
       + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
   )

Note: the multiplication is *before* realization. ``law`` is an
:class:`OperatorSum` containing :class:`ScaledOperator` wrappers
around the **leaf** :class:`BoundaryTraceLaw` instances — neither
of those is itself a :class:`BoundaryTraceLaw`, so passing this
``law`` directly to :meth:`SNBoundaryRealizer.realize` would raise
:class:`BoundaryError`.

:func:`~orpheus.sn.boundary_realize.realize_recursively` walks the
expression tree, realizing every :class:`BoundaryTraceLaw` leaf
via the SN realizer and reassembling the result through the same
Wave-0 composers:

.. code-block:: python

   from orpheus.sn.boundary_realize import realize_recursively

   realized = realize_recursively(law, method_space)
   # realized is:
   #   OperatorSum(
   #       ScaledOperator(0.3, PermutationOperator(...)),
   #       ScaledOperator(0.7, AngularAverageOperator(...)),
   #   )

The walker recognizes the five Wave-0 composers
(:class:`OperatorSum`, :class:`OperatorProduct`,
:class:`ScaledOperator`, :class:`TensorProductOperator`,
:class:`SumOfTensorProductsOperator`); leaves must be
:class:`BoundaryTraceLaw` instances. Unknown node types raise
:class:`BoundaryError` with the offending type in
:attr:`law` so the failure is debuggable.

Per the "Unify after two instances" architectural rule, the walker
currently lives at :mod:`orpheus.sn.boundary_realize` —
SN-specific. When a second functional realizer ships, the walker
will move to ``orpheus.geometry.boundary.realize`` and become
method-agnostic via
:meth:`BoundaryRealizerRegistry.get(method_name)` lookup.


.. _bc-vacuum-semantic-correction:

The vacuum semantic correction (§16A.5)
=======================================

The most subtle design decision of the refactor concerns vacuum.
Pre-Wave-7 the legacy ``VacuumBoundaryOperator.apply`` body was:

.. code-block:: python

   def apply(self, psi_out: np.ndarray, quadrature) -> np.ndarray:
       return np.zeros_like(psi_out)

This returns **all zeros**, including the outflow ordinates that
the BC has no physical interpretation for (vacuum says nothing
about :math:`\gamma_+ \psi`; it only sets :math:`\gamma_- \psi = 0`).

The post-Wave-8 SN realizer's vacuum branch returns
:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor` —
a sparse mask that zeros **only the inflow ordinates** and
preserves the outflow trace. This is the §16A.10 trace-correct
representation: the operator's apply is a projector onto the
inflow ordinate subspace, which is the right algebraic object
for the affine law :eq:`affine-bc-form` to read.

Why this matters
----------------

Three downstream consequences make the inflow-only mask the right
contract:

1. **Sensitivity adjoints.** A future adjoint sensitivity path
   needs the outflow trace preserved to compute the response
   of an outgoing-current functional to the inflow BC. The
   zeros-all body would silently lose the gradient at the
   boundary.
2. **Compositional clarity.** The realized vacuum operator is
   now self-adjoint and idempotent (it's a projector). The
   zeros-all operator is also idempotent but is the
   ``ZeroOperator`` projector, which is the wrong type tag
   for "inflow-mask only".
3. **Algebraic uniformity.** Every other rank-1 law (reflective,
   white, albedo) acts on the **inflow ordinates only** and
   leaves the outflow rows untouched. The legacy vacuum was
   the asymmetric special case; the post-Wave-8 vacuum joins
   the family.

The §16A.5 compatibility audit
------------------------------

The semantic correction is **observably equivalent** at every
existing production call site because every site reads inflow
rows only. The Wave 8 close-out audited all 13 ``bc.apply(...)``
sites in :mod:`orpheus.sn.sweep` and :mod:`orpheus.sn.operator`:

* ``sweep.py:334,351`` (1-D slab) read
  ``psi_face_left_in[n_half + n]`` for positive-μ ordinates
  only.
* ``sweep.py:508,654`` (spherical / cylindrical) gate the read
  by ``mu_n < 0`` / ``eta_n < 0`` (inflow at the outer face).
* ``sweep.py:843-854`` (2-D wavefront) reads
  ``full_face_x[oct_idx]`` only, where ``oct_idx ↔ inflow
  ordinates`` by construction.
* ``operator.py:230-256`` (2-D FD matvec) has explicit
  ``mu_x[n] > 1e-15`` / ``mu_y[n] < -1e-15`` gates per face.
* ``operator.py:530`` (spherical FD matvec) gates on
  ``quad.mu_x[n] < -1e-15``.

No call site reads the outflow rows of the BC result, so the
semantic correction has no observable effect on existing tests.
The Wave 6 snapshot harness gates this explicitly: the
``vacuum_lebedev17`` case compares **inflow-row outputs only**,
documenting the intentional semantic divergence in a comment
adjacent to the test case.

Why "Option a" (Wave 7) over the alternatives
---------------------------------------------

The Wave 7 brief considered three migration strategies for the
2-arg legacy path:

* **(a) Vacuum-stays-legacy.** Vacuum's standalone
  :meth:`apply(psi_out, quad)` preserves the pre-refactor
  zeros-all body until the §16A.5 inflow-only mask activates
  via the realizer path. ``VacuumInflow.apply`` body keeps the
  legacy ``np.zeros_like(psi_out)`` form for backward
  compatibility with tests that pin standalone-apply behavior.
* **(b) Face-aware BC.** Add a ``face`` constructor argument
  to ``VacuumInflow`` so the standalone apply knows which
  ordinates to mask. Would require every existing test that
  instantiates ``VacuumBoundaryOperator()`` to pass a face
  argument.
* **(c) Combined ABC merge.** Move the apply body into a
  shared base class so both rank-1 vacuum and rank-N mixed BCs
  share one code path. Would not address the face-awareness
  question (no face → still zeros-all).

Option (a) was chosen because:

* The realizer path is the architecturally correct one; the
  standalone-apply path is a transitional artifact. Distorting
  the law's interface (Option b) to better serve the
  transitional path would inflate every user-facing test
  signature for one wave's duration.
* The compatibility audit above confirmed no production
  consumer reads outflow rows, so the legacy zeros-all on the
  standalone path produces identical results to the inflow-only
  mask on the realizer path **at every observable point**.
* The Wave 6 snapshot harness already excluded vacuum from the
  bit-identity contract (intentional documented semantic
  divergence). No new test failures.

Option (a) is documented in the
:class:`VacuumInflow` docstring and the Wave 7 close-out memo;
the corresponding test sites at
``tests/geometry/test_boundary.py:87,281,282`` carry inline
comments marking the retained legacy semantics.


.. _bc-numerical-evidence:

Numerical evidence
==================

The Wave 6 snapshot harness
(:mod:`tests.geometry.test_bc_equivalence_snapshot`) is the
load-bearing equivalence pin for the realizer-vs-legacy migration.
Eight snapshot cases pin pre-refactor outputs as ``.npz`` files;
post-refactor outputs are compared at ``nulp ≤ 4`` (or ``nulp = 64``
for cases where the reduction tree intentionally changed). The
cases:

.. list-table:: Wave 6 BC equivalence snapshot cases
   :header-rows: 1
   :widths: 28 28 18 26

   * - Snapshot
     - BC
     - Quadrature
     - Tolerance
   * - ``vacuum_lebedev17``
     - ``VacuumInflow()``
     - Lebedev order 17
     - inflow-rows-only
       (intentional §16A.5 semantic divergence on outflow)
   * - ``albedo_05_lebedev17``
     - ``AlbedoBoundary(0.5)``
     - Lebedev order 17
     - ``nulp ≤ 4``
   * - ``specular_x_lebedev17``
     - ``ReflectiveBoundary(axis="x", albedo=1.0)``
     - Lebedev order 17
     - bit-identical (α=1 fast path)
   * - ``specular_y_partial_07_LS6``
     - ``ReflectiveBoundary(axis="y", albedo=0.7)``
     - LevelSymmetricSN(6)
     - ``nulp ≤ 4``
   * - ``white_xmax_LS4``
     - ``WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)``
     - LevelSymmetricSN(4)
     - bit-identical (α=1 fast path)
   * - ``white_xmin_partial_03_GL``
     - ``WhiteBoundary(axis="x", outward_sign=-1, albedo=0.3)``
     - GaussLegendre1D(8)
     - ``nulp ≤ 4``
   * - ``periodic_lebedev17``
     - ``PeriodicBoundary()``
     - Lebedev order 17
     - bit-identical (angular identity body)
   * - ``mixed_30spec_70white_LS4``
     - ``0.3 * spec + 0.7 * white`` (Wave-0 algebra)
     - LevelSymmetricSN(4)
     - ``nulp ≤ 64``

The 1-ULP tolerance for α=1 specular / white reflects the
α=1.0 fast path returning the bare primitive (no
``ScaledOperator`` wrap → same FP reduction tree as the legacy
body). The 4-ULP tolerance for partial-albedo cases reflects the
``ScaledOperator(α, primitive).apply(x) = α * primitive.apply(x)``
multiplication being the single new floating-point operation per
ordinate (one extra rounding step). The 64-ULP tolerance for the
Marshak case reflects the
:class:`OperatorSum` reduction order (the legacy ``MixedBoundaryOperator``
loop summed left-to-right via accumulation; the post-Wave-11
:class:`OperatorSum` binary tree adds in a different but
mathematically equivalent order — see the
``vv-principles`` skill's "bit-identity vs principled-equivalence"
discussion for the general framing of this kind of FP-non-
associativity drift).

The snapshot ``.npz`` files live at
``tests/geometry/snapshots/bc_equivalence_*.npz`` and are committed
to the repository. The generator script
``tests/geometry/_generate_bc_equivalence_snapshots.py`` regenerates
them on the legacy commit; the comparison test enforces
equivalence on every subsequent commit. Wave 11 (``MixedBoundaryOperator``
removal) regenerated the mixed snapshot — the new payload was
**bit-identical** to the pre-Wave-11 payload because the realizer's
deleted internal mixed-BC path was already composing via
``OperatorSum``-of-``ScaledOperator``.


.. _bc-cartesian-curvilinear-split:

Cartesian / curvilinear split
=============================

For technical reasons :class:`~orpheus.numerics.trace_space.InflowTraceSpace`
is implemented for **Cartesian** meshes only (Mesh1D with
``coord=CARTESIAN`` and Mesh2D with ``coord=CARTESIAN``). The
curvilinear case (Mesh1D ``SPHERICAL`` / ``CYLINDRICAL``, Mesh2D
``CYLINDRICAL``) raises :class:`NotImplementedError` from the
:meth:`InflowTraceSpace.from_mesh_and_quadrature` factory.

Why? Because no current curvilinear-Krylov consumer needs the
per-face inflow mask — the curvilinear sweep paths compute the
inflow / outflow predicate on the fly inside the inner loop. The
deferral is grep-able for the future consumer that does need the
sparse representation.

This forces :meth:`SNMesh._resolve_one` to take two branches:

.. list-table:: Resolution path by mesh coord
   :header-rows: 1
   :widths: 28 36 36

   * - Mesh coord
     - Path
     - Resolved BC type
   * - ``CARTESIAN`` (1-D slab + 2-D Cartesian)
     - ``SNBoundaryRealizer().realize(law, method_space)``
       wrapped in
       :class:`_BoundBoundaryOperator(realized, kind=bc.kind)`
     - 1-arg Wave-0 ``LinearOperator``
   * - ``SPHERICAL`` / ``CYLINDRICAL``
     - bare law (e.g. ``ReflectiveBoundary(axis="x")``) wrapped in
       :class:`_BoundBoundaryOperator(law, quadrature=self.quad, kind=bc.kind)`
     - 2-arg legacy law; shim's ``apply(psi)`` forwards
       ``inner.apply(psi, bound_quad)``

Both paths present a uniform 1-arg ``apply(psi)`` contract to the
sweep call sites (Wave 9 migrated 13 sites from 2-arg to 1-arg);
they differ only in what the shim wraps internally.

The curvilinear path is **bit-identical** to the pre-Wave-9 direct
``bc.apply(psi, quad)`` because the shim's bound quadrature is the
same :class:`~orpheus.sn.quadrature.AngularQuadrature` object the
call sites used to pass. Same object, same FP reduction tree, same
output.

The deferred work is tracked at:

* Issue #176 — "BC refactor: drop 2-arg apply + delete
  ``_BoundBoundaryOperator`` once curvilinear realizer ships"
  (``module:sn, module:geometry, type:improvement``).
* "BC: curvilinear InflowTraceSpace.from_mesh_and_quadrature
  implementation" follow-up (``module:geometry, type:feature``).

When both ship, the curvilinear path will route through the
realizer, the shim is deleted, and every consumer becomes
``LinearOperator``-only.


Anti-pattern catalog
====================

A short list of patterns the refactor's authors considered and
rejected, with the reasoning preserved so future sessions don't
re-attempt them:

1. **Single ``BoundaryOperator`` ABC carrying both law and
   realizer responsibilities.** This is the pre-refactor shape.
   Rejected because the law / realizer split is the architectural
   point of the refactor: laws are method-agnostic, realizers are
   method-specific. Keeping them in one class would force every
   law to know about the SN sweep's inflow-mask plumbing, which is
   precisely what Wave 0 / Wave 1 / Wave 2 was supposed to abstract
   away.
2. **Dedicated ``MixedBoundaryOperator`` class.** Pre-Wave-11
   shape. Rejected (deleted Wave 11) — the Wave-0 algebra dunders
   on every :class:`BoundaryTraceLaw` already produce
   :class:`OperatorSum` shapes; the dedicated class added a
   second realization path with no semantic difference. See
   :ref:`bc-rank-n-algebra`.
3. **Shared ``BoundaryRealizerBase`` ABC for cross-method
   realizers.** Considered Wave 5. Rejected per the "Unify after
   two instances" architectural discipline: only one functional
   realizer ships today. Building the abstraction on a single
   instance would force a particular shape on every future method
   based on SN's current dispatch idiom (``isinstance``); MoC
   might want a different dispatch shape (e.g. registry-keyed by
   law class name) that the ABC would foreclose.
4. **Adding ``face`` to ``VacuumInflow``'s constructor for
   semantic correctness on the standalone-apply path.** Option (b)
   in :ref:`bc-vacuum-semantic-correction`. Rejected because it
   would have inflated every test signature for one wave to fix a
   path that the refactor was retiring anyway. The transitional
   legacy-vacuum body is the right cost.
5. **Auto-importing every cross-method realizer at
   :mod:`orpheus.geometry.boundary` import time.** Considered to
   make :meth:`BoundaryRealizerRegistry.get("MoC")` work without
   the caller having to ``import orpheus.moc`` first. Rejected
   because :mod:`orpheus.sn` is a heavy module that every consumer
   of the boundary package would then pay for. The current rule
   (each method's :mod:`__init__.py` imports its own realizer; the
   SN realizer is imported explicitly by
   :class:`SNMesh.__init__`) is the right cost-localization.


References
==========

* Grand Report v3 §16, §16A.1–5 (affine boundary form + trace
  structure), §16A.10 (sparse trace primitives), §16A.11 (dual
  registry), §16A.12 + §27.6 (universal invariants), §26A.4
  (named-error catalog). Source: ``.claude/plans/neutron_transport_grand_report_v3.md``.
* The 12-wave refactor plan:
  ``.claude/plans/transient-giggling-cake.md``.
* Lewis, E. E. & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*. American Nuclear Society. §3.4 (boundary
  conditions in transport).
* Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.5 (albedo, white, and Marshak
  boundary conditions).
* The tensor decomposition equation :eq:`bc-tensor-decomposition`
  at :ref:`bc-tensor-decompositions` (in
  :doc:`discrete_ordinates`) shows the algebra
  :math:`R = \sum_\alpha G_\alpha \otimes A_\alpha` that this page
  refines into the affine form :eq:`affine-bc-form`.
* :ref:`operator-algebra` for the Wave-0 primitives the realized
  BCs decompose into.
* The V&V error catalog in the ``vv-principles`` skill
  (``.claude/skills/vv-principles/error_catalog.md``) carries
  the ERR-040..ERR-047 entries in canonical form.
