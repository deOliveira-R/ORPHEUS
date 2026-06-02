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
  | 1     | Trace structure       | :mod:`orpheus.numerics.spaces.trace_space`  |
  |       | (Γ\_-, Γ\_+ + mask)   | (all Mesh1D coord systems + 2-D Cartesian;  |
  |       |                       | 2-D cylindrical Mesh2D deferred)            |
  +-------+-----------------------+---------------------------------------------+
  | 2     | Boundary law          | :mod:`orpheus.geometry.boundary` (ABC +     |
  |       | (method-agnostic)     | 6 concrete laws, dual-registry: kind-keyed) |
  +-------+-----------------------+---------------------------------------------+
  | 3     | Method realizer       | :mod:`orpheus.sn.boundary_realizer`         |
  |       | (per-method strategy) | (SN functional) + 4 method stubs            |
  +-------+-----------------------+---------------------------------------------+

- Rank-N boundary conditions (Marshak, partial-current mixes) are
  expressed via a **descriptor-tree algebra** on the unrealised laws
  themselves. The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
  algebra dunders (``+``, ``-``, ``*``, ``/``, unary ``-``) return
  :class:`~orpheus.geometry.boundary.LawSum` /
  :class:`~orpheus.geometry.boundary.LawScaled` nodes — a closed
  algebra over ``BoundaryTraceLaw | LawSum | LawScaled``. The tree is
  a **pure descriptor** with no ``apply`` method; the
  :func:`~orpheus.sn.boundary_realize.realize_recursively` type
  transformer walks it once per face and emits an operator tree of
  :class:`~orpheus.numerics.operator.OperatorSum` /
  :class:`~orpheus.numerics.operator.ScaledOperator` composers around
  realised 1-arg leaves. See :ref:`bc-trace-law-descriptor-model` and
  :ref:`bc-rank-n-algebra`. There is no dedicated
  ``MixedBoundaryOperator`` class (retired Wave 11); there is also no
  ``apply`` method on the raw law (retired Issue #186, B3 + β2).
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
- **Vacuum semantic correction (§16A.5).** ``VacuumInflow`` realises
  to :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`,
  which zeroes **only the inflow ordinates** and preserves the
  outflow trace. The pre-refactor ``VacuumBoundaryOperator.apply``
  used ``np.zeros_like(psi_out)`` and zeroed everything; the §16A.10
  inflow-only mask is the trace-correct representation.
  **Post Issue #186 (2026-05-11)** the realizer-routed
  inflow-only mask is the **sole** path to vacuum action — the
  zeros-all body has been deleted along with every other
  :meth:`apply` method on every concrete BC. The §16A.5
  two-paths-divergence is therefore eliminated by design (no
  second path remains). The realizer-routed mask is uniform across
  every supported mesh (1-D Cartesian / spherical / cylindrical +
  2-D Cartesian) since Issue #188 lifted the curvilinear
  :class:`InflowTraceSpace` deferral.

.. admonition:: V&V status

   This page is **L4-informational** with respect to correctness.
   The architecture documented here is structural — it makes the
   code understandable and composable but does not by itself verify
   any equation. The verification load is carried by:

   - L0 foundation tests on individual primitives
     (:mod:`tests.numerics`,
     :mod:`tests.geometry.test_boundary_trace_law`,
     :mod:`tests.geometry.test_bc_errors`).
   - L1 realiser-output snapshot tests (the Wave-6 harness at
     :mod:`tests.geometry.test_bc_equivalence_snapshot`, with the
     legacy halves dropped post Issue #186 / C-B3.7). The
     surviving ``test_realizer_*`` halves pin the realised-operator
     output against committed ``.npz`` snapshots at ``nulp ≤ 4``
     for non-vacuum BCs (intentional semantic capture of the
     §16A.5 inflow-only mask for vacuum).
   - L1 descriptor-tree algebra tests
     (:mod:`tests.geometry.test_law_composition`) pinning the
     :class:`LawSum` / :class:`LawScaled` closed-algebra contract
     (18 tests: foundation + L1 coverage).
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

.. note::

   **SN apply matvec honours the affine BC contract (Issue #168
   Phase C, 2026-05-12).** The
   :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
   and ``_cylindrical`` matvecs were rewritten as one sweep
   iteration semantically: the BC trace law is applied **at least
   once** per matvec at the boundary edge on the WDD-propagated
   outflow face values (:math:`\gamma_+ \psi`), not on cell-centre
   approximations. The pre-Phase-C cell-centre-as-face-value
   contamination — and the Phase A ``BoundaryFaceFlux`` Protocol
   that patched it — both retire in Phase C. See
   :ref:`bc-trace-contract-respected-by-matvec` for the
   verification gate that pins this contract, and
   :ref:`bc-phase-d-two-bc-applies-per-matvec` for the Phase D
   strengthening that audits **two** BC apply calls per matvec
   (Phase D Carlson context + Phase C trace law).


Layer 1 — trace structure
-------------------------

The trace operators :math:`\gamma_\pm` carry their domain
information in two typed :class:`~orpheus.numerics.space.FunctionSpace`
subclasses:

* :class:`~orpheus.numerics.spaces.trace_space.InflowTraceSpace` represents
  :math:`\Gamma_- = \{(\mathbf{r}, \Omega) \in \partial\Omega \times S^2
  : \Omega \cdot \hat n(\mathbf{r}) < 0\}` — the per-face directional
  half of the boundary on which the incoming angular flux is
  constrained by the law.
* :class:`~orpheus.numerics.spaces.trace_space.OutflowTraceSpace` represents
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
  :meth:`~orpheus.numerics.spaces.trace_space.InflowTraceSpace.inflow_indices_for_face`.
  Those indices are the constructor argument to
  :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`.
* The universal invariant
  :meth:`~orpheus.geometry.boundary.BoundaryTraceLaw.assert_source_lives_on_incoming_trace`
  uses the inflow mask to validate that a
  :class:`~orpheus.geometry.boundary.InflowSourceSpec` has no nonzero
  entries on outflow ordinates (per
  :class:`~orpheus.geometry.boundary.BoundarySourceNotOnIncomingTraceError`,
  ERR-047).
* The SN curvilinear sweep (1-D spherical / cylindrical) consumes
  the same realizer-routed mask as the Cartesian path — Issue #188
  (C188.1+C188.2 in :mod:`orpheus.numerics.spaces.trace_space`, C188.3 in
  :mod:`orpheus.sn.geometry`) lifted the curvilinear deferral and
  Issue #176 then dropped the legacy 2-arg shim that existed only
  to bridge that deferral.

The class :class:`~orpheus.numerics.spaces.trace_space.InflowTraceSpace`
carries the mask as an ``Optional[np.ndarray]`` field excluded from
:meth:`__eq__` and :meth:`__hash__` — preserving the
:class:`~orpheus.numerics.space.FunctionSpace` identity convention
``(name, shape)``. Construction goes through the classmethod
factory :meth:`InflowTraceSpace.from_mesh_and_quadrature`; the bare
dataclass constructor is reserved for trace spaces whose mask will
be populated later (or never).

**Coord-system coverage (post Issue #188).** The factory supports
every :class:`~orpheus.geometry.mesh.Mesh1D` coord system
(``CARTESIAN`` / ``SPHERICAL`` / ``CYLINDRICAL``) because they all
share the same ``("left", "right")`` face structure and the
:class:`~orpheus.sn.quadrature.GaussLegendre1D` adapter's
:math:`\mu_x` is the direction cosine along the chosen axis (the
Cartesian-x axis for slab; the radial axis for curvilinear). The
:class:`~orpheus.geometry.mesh.Mesh2D` factory supports
``coord=CARTESIAN`` only — 2-D cylindrical (axisymmetric
:math:`(r, z)`) ``Mesh2D`` raises :class:`NotImplementedError` and
will continue to do so until ORPHEUS gains a 2-D cylindrical SN
sweep (a face-normal table and azimuthal-averaging mask predicate
will then be needed).


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
5. **No ``apply`` method at all** (Issue #186 / B3 + β2,
   2026-05-11). The descriptor model that survived the C176.3
   Option A interim was retired in favour of a pure-descriptor
   contract: :class:`BoundaryTraceLaw` is no longer a
   :class:`~orpheus.numerics.operator.LinearOperatorMixin` subclass,
   no concrete law carries ``apply`` / ``apply_transpose``
   methods, and no ``capabilities`` ``ClassVar`` is defined.
   The §16A.3 three-layer split (descriptor / realizer / operator)
   is now enforced by the **type system**, not by convention —
   ``law.apply(psi)`` on a raw law is an ``AttributeError`` at
   runtime and a static-type error at the linter level. The
   :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` is the
   sole bridge from descriptor to callable; see
   :ref:`bc-trace-law-descriptor-model` for the design rationale
   and the predecessor approaches that were tried and rejected.

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
and Marshak / partial-current boundaries are rank-N via the
**descriptor-tree algebra** on the unrealised laws (:class:`LawSum`
/ :class:`LawScaled` over :class:`BoundaryTraceLaw` leaves) —
realised once per face by
:func:`~orpheus.sn.boundary_realize.realize_recursively`. See
:ref:`bc-rank-n-algebra` below.


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
     - Replaced by the descriptor-tree algebra
       (:class:`LawSum` / :class:`LawScaled` over
       :class:`BoundaryTraceLaw` leaves; Issue #186 / B3 + β2);
       the dedicated class added no value over the inherited
       algebra dunders.


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
       :class:`~orpheus.numerics.spaces.trace_space.InflowTraceSpace`.
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

   from orpheus.numerics.spaces.trace_space import InflowTraceSpace

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

Step 5 — wrap with a kind tag
-----------------------------

Every ``SNMesh.bc_*`` attribute carries a uniform 1-arg
``apply(psi)`` contract (Wave 9 migrated 13 production sites from
2-arg to 1-arg). Post Issue #186 / C-B3.4 the
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim is a **strict 1-arg passthrough** that adds two pieces of
metadata to the realized operator:

* a free-form ``kind`` string tag carrying the originating
  :class:`~orpheus.geometry.mesh.BC` kind (``"vacuum"`` /
  ``"reflective"`` / ...) — load-bearing for the
  ``sn_mesh.bc_left == "vacuum"`` string-equality surface that
  several legacy SN tests and the BC-resolution diagnostic rely
  on;
* :meth:`capabilities` delegation to the wrapped inner operator
  so consumers composing the shim with other Wave-0 primitives
  inherit the right feature set.

The shim's ``apply`` / ``apply_transpose`` signatures are strict
1-arg ``(self, psi)`` — extra positional or keyword arguments
raise :class:`TypeError`. The pre-Issue-#186 affordance that
swallowed ``*_extra, **_kw`` was the last remnant of the 2-arg
legacy era; it was dropped in C-B3.4 alongside the descriptor
cleanup because every production and test call site is now
strict 1-arg.

.. code-block:: python

   from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator

   return _BoundBoundaryOperator(realized, kind=bc.kind)

The shim is **internal** to the package (not in :attr:`__all__`)
— a test pins its private status.

**Historical note (pre Issue #176 / #186).** The Wave-8/9
implementation carried an optional ``quadrature=`` kwarg that,
when non-``None``, bound an
:class:`~orpheus.sn.quadrature.AngularQuadrature` and forwarded
``inner.apply(psi, bound_quad)`` to a legacy 2-arg
:class:`BoundaryTraceLaw` body. That dual-mode existed ONLY
because :meth:`InflowTraceSpace.from_mesh_and_quadrature` raised
:class:`NotImplementedError` for curvilinear ``Mesh1D``, which
forced :meth:`SNMesh._resolve_one` to bypass the realizer for
spherical / cylindrical meshes. Issue #188 lifted that
deferral; Issue #176 (C176.1) dropped the bound-quadrature mode
here because no production-issued shim carried
``_quadrature is not None`` after C188.3. Issue #186 (C-B3.4)
then dropped the residual ``*_extra, **_kw`` argument-swallow
because, with concrete-BC ``apply`` methods retired and all
production / test sites strict 1-arg, the defensive net was
dead code.

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
     - User-supplied :class:`InflowSourceSpec` has nonzero entries on
       :math:`\Gamma_+`; geometrically meaningless and indicates a
       wrong source-shape contract.

All eight extend :class:`ValueError` (via the
:class:`~orpheus.geometry.boundary.BoundaryError` base) so existing
``except ValueError`` consumers from the pre-refactor code keep
working. Every catch site can additionally pattern-match on the
typed subclass to recover the offending law name from
:attr:`BoundaryError.law`.


.. _bc-rank-n-algebra:

Descriptor-tree algebra for rank-N boundaries
=============================================

Rank-N (Marshak, partial-current) boundary conditions are
**not** a special class. They are expressed via a closed
**descriptor-tree algebra** over
``BoundaryTraceLaw | LawSum | LawScaled`` nodes — pure
declarative structure with **no** ``apply`` method on any node.
The :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` algebra
dunders (``+``, ``-``, ``*``, ``/``, unary ``-``) return
:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled` instances, never
operators. The :func:`~orpheus.sn.boundary_realize.realize_recursively`
type transformer is the **sole** path from descriptor tree to
operator tree.

The §15.2 sum-of-tensor-products form

.. math::
   :label: bc-rank-n-tensor-decomposition

   R \;=\; \sum_{\alpha} c_{\alpha}\, G_{\alpha},
   \qquad c_{\alpha} \in \mathbb{R},
   \quad G_{\alpha} \in
   \{\text{permutation, average, mask, wrap, identity, source}\},

maps onto the LawXxx algebra as ``c1 * law_1 + c2 * law_2 + ...``,
where each ``c_i * law_i`` term is a :class:`LawScaled` node and
the sum is a :class:`LawSum` node.

The standard Marshak boundary (Bell & Glasstone 1970 §1.5) — a
mix of specular reflection (weight :math:`c_1`) and diffuse
white reflection (weight :math:`c_2`) — is:

.. code-block:: python

   from orpheus.geometry.boundary import (
       LawScaled, LawSum,
       ReflectiveBoundary, WhiteBoundary,
   )
   from orpheus.sn.boundary_realize import realize_recursively
   from orpheus.sn.boundary_realizer import SNMethodSpace

   # Build the descriptor tree — no realization yet.
   spec = ReflectiveBoundary(axis="x", albedo=1.0)
   white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
   marshak_law = 0.3 * spec + 0.7 * white
   # marshak_law is:
   #   LawSum(
   #       LawScaled(0.3, ReflectiveBoundary(axis="x", albedo=1.0)),
   #       LawScaled(0.7, WhiteBoundary(axis="x", outward_sign=+1,
   #                                    albedo=1.0)),
   #   )
   # NOT callable — no .apply method on LawSum or its children.
   assert not hasattr(marshak_law, "apply")

   # Realize the tree at one face. realize_recursively walks
   # LawSum / LawScaled / leaf-law nodes and emits the matching
   # Wave-0 operator-tree composers around realized 1-arg leaves.
   ms = SNMethodSpace.for_face(...)
   marshak_op = realize_recursively(marshak_law, ms)
   # marshak_op is:
   #   OperatorSum(
   #       ScaledOperator(0.3, PermutationOperator(...)),
   #       ScaledOperator(0.7, AngularAverageOperator(...))
   #   )
   psi_in = marshak_op.apply(psi_out)   # 1-arg LinearOperator

The output is a Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator`-wrapped Wave-0
primitives, consumable by the SN sweep / Krylov path via the
uniform 1-arg :meth:`apply`. The descriptor-tree algebra and the
operator-tree algebra are **separate type families**: the
descriptor tree is built with ``LawXxx`` nodes that have **no**
``apply``; the operator tree is built with ``OperatorXxx`` nodes
that **do** have ``apply``. The two families never inter-compose —
mixing a :class:`LawNode` with an already-realized
:class:`~orpheus.numerics.operator.LinearOperator` via ``+`` is a
type error; the user MUST call :func:`realize_recursively` first.

Closed-algebra guarantees
-------------------------

* **Constant folding on scalars.**
  ``LawScaled(α, LawScaled(β, x))`` collapses to
  ``LawScaled(α * β, x)`` at construction time. The intermediate
  ``LawScaled`` nesting never appears at rest, which keeps the tree
  shallow under repeated scalar multiplication.
* **No associativity flattening on sums.** ``(a + b) + c`` is
  :class:`LawSum(LawSum(a, b), c)`, distinct from
  :class:`LawSum(a, LawSum(b, c))`. The walker treats both shapes
  identically — the realized output is the same Wave-0 operator
  algebra value up to floating-point non-associativity in the
  final sum.
* **Subtraction rewrites via :class:`LawScaled(-1, ...)`.** The
  unary ``-`` operator and the binary ``-`` operator both produce
  trees containing only :class:`LawSum` / :class:`LawScaled`
  nodes — there is no dedicated ``LawDifference`` type.
* **Division rewrites via :class:`LawScaled(1/α, ...)`.** Pure
  syntactic sugar for ``LawScaled(α, ...).__truediv__``.

The pre-refactor implementations
--------------------------------

Two prior approaches converged on the present descriptor-tree
design through empirical falsification:

**Wave 11 (~2026-03)** — ``MixedBoundaryOperator(components:
list[tuple[float, BoundaryOperator]])`` class whose
:meth:`apply` body looped over ``components`` and summed
``coeff * primitive.apply(psi, quad)``. The SN realizer
dispatched on it via an ``isinstance(law, MixedBoundaryOperator)``
branch that ran the same loop with
``coeff * realize(primitive, ms)`` summed via
:class:`OperatorSum`. Wave 11 deleted this code because the
delayed-realization-by-container pattern broke down once vacuum
needed per-face inflow indices that the bare-law container had
no access to.

**β1 interim landing (Issue #186 / B3, ~2026-04)** — every
:class:`BoundaryTraceLaw` inherited the Wave-0 operator-algebra
dunders from :class:`~orpheus.numerics.operator.LinearOperatorMixin`,
so writing ``0.3 * spec + 0.7 * white`` directly produced an
:class:`OperatorSum` of :class:`ScaledOperator`-wrapped raw
:class:`BoundaryTraceLaw` leaves (NOT realized). The
:func:`realize_recursively` walker then traversed the Wave-0
composer tree, realized each leaf, and emitted a parallel
operator tree. This achieved the right algebraic shape but
**violated the type system**: the resulting tree was an
:class:`OperatorSum` instance (a :class:`LinearOperator`!) whose
:meth:`apply` could not actually be called before realization —
calling it raised :class:`BoundaryError` at apply-time because the
leaves were laws, not operators. The convention "you must realize
this OperatorSum before calling apply" was a runtime contract that
the type system did nothing to enforce. β1 retained the
ergonomic of "the same ``+`` operator before and after
realization" but at the cost of conflating two type families.

**β2 (this scope, Issue #186 / B3 + β2)** — separates the two
type families explicitly: :class:`LawSum` / :class:`LawScaled` for
the descriptor tree (no :meth:`apply`); :class:`OperatorSum` /
:class:`ScaledOperator` for the operator tree (with
:meth:`apply`). The static type system enforces "you cannot call
this until it's realized" — :class:`LawSum` has no :meth:`apply`
method on the class, so the linter flags ``tree.apply(...)``
without running the program. The ergonomic of "the same ``+``"
survives because both the law-tree and the operator-tree use the
same Python ``+`` syntax; the runtime dispatch on type tells the
reader (and the type checker) which algebra is in effect.

The Wave 6 snapshot harness verifies that the bit-identity of the
Marshak case ``0.3 * spec + 0.7 * white`` is **preserved** by the
β2 transition: the realized
``OperatorSum(ScaledOperator, ScaledOperator)`` reduction tree
matches the β1-era output to the same ULP tolerance, because the
operator-tree shape after :func:`realize_recursively` is
algebraically identical.


.. _bc-realize-recursively:

The ``realize_recursively`` walker — a descriptor → operator type transformer
=============================================================================

:func:`~orpheus.sn.boundary_realize.realize_recursively` is the
**type transformer** from the descriptor-tree algebra
(``BoundaryTraceLaw | LawSum | LawScaled``) to the operator-tree
algebra (``LinearOperator`` with
:class:`~orpheus.numerics.operator.OperatorSum` /
:class:`~orpheus.numerics.operator.ScaledOperator` composers
around realized 1-arg leaves). Calling it is the **only** path
from a non-callable descriptor to a callable operator.

The dispatch is exhaustive on the descriptor-tree node types:

.. code-block:: python

   def realize_recursively(
       node: BoundaryTraceLaw | LawSum | LawScaled,
       method_space: SNMethodSpace,
   ) -> LinearOperator:
       if isinstance(node, BoundaryTraceLaw):
           # Leaf: realize via the SN realizer registry.
           return SNBoundaryRealizer().realize(node, method_space)
       if isinstance(node, LawScaled):
           # Scalar-times-law: wrap the realized inner in ScaledOperator.
           inner_op = realize_recursively(node.inner, method_space)
           return ScaledOperator(node.scalar, inner_op)
       if isinstance(node, LawSum):
           # Sum: realize each side, wrap in OperatorSum.
           a_op = realize_recursively(node.a, method_space)
           b_op = realize_recursively(node.b, method_space)
           return OperatorSum(a_op, b_op)
       raise TypeError(
           f"realize_recursively expected BoundaryTraceLaw | LawSum | "
           f"LawScaled, got {type(node).__name__}."
       )

Usage on the descriptor tree:

.. code-block:: python

   from orpheus.geometry.boundary import (
       ReflectiveBoundary, WhiteBoundary,
   )
   from orpheus.sn.boundary_realize import realize_recursively

   # Build the descriptor tree.
   law = (
       0.3 * ReflectiveBoundary(axis="x")
       + 0.7 * WhiteBoundary(axis="x", outward_sign=+1)
   )
   # law is LawSum(LawScaled(0.3, ...), LawScaled(0.7, ...)).

   # Realize once, at face resolution time.
   realized = realize_recursively(law, method_space)
   # realized is:
   #   OperatorSum(
   #       ScaledOperator(0.3, PermutationOperator(...)),
   #       ScaledOperator(0.7, AngularAverageOperator(...)),
   #   )
   psi_in = realized.apply(psi_out)   # 1-arg LinearOperator

Type-system contract
--------------------

The walker's input is intentionally narrow:

* :class:`BoundaryTraceLaw` instances and the two descriptor-tree
  composer dataclasses (:class:`LawSum`, :class:`LawScaled`) are
  the only valid node shapes.
* Wave-0 operator-tree composers
  (:class:`~orpheus.numerics.operator.OperatorProduct`,
  :class:`~orpheus.numerics.operator.TensorProductOperator`,
  :class:`~orpheus.numerics.operator.SumOfTensorProductsOperator`,
  :class:`~orpheus.numerics.operator.OperatorSum`,
  :class:`~orpheus.numerics.operator.ScaledOperator`) are **not**
  recognized — they belong to the operator tree, not the
  descriptor tree, so they should never appear in the realizer's
  input.
* Unknown nodes raise :class:`TypeError` (not
  :class:`BoundaryError`) with the offending type name in the
  message, because this is a **typing** failure (caller passed
  the wrong kind of object), not a BC-domain failure.

The β1 → β2 transition (see :ref:`bc-rank-n-algebra`) sharpened
the walker's type signature from "any Wave-0 composer tree with
:class:`BoundaryTraceLaw` leaves" to "any descriptor-tree node".
The dispatch table shrank from five Wave-0 composers + leaf to
three descriptor types + leaf (counting the leaf as the same
:class:`BoundaryTraceLaw` branch). The eliminated branches
(:class:`OperatorProduct`, :class:`TensorProductOperator`,
:class:`SumOfTensorProductsOperator`) handle operator composition
patterns that have no descriptor-tree analog — they were dead
dispatch paths once :class:`LawSum` / :class:`LawScaled` replaced
the in-tree Wave-0 algebra. Removing them clarified the walker's
role: it is **exactly** the type transformer between the two
algebras, nothing more.

Per the "Unify after two instances" architectural rule, the walker
currently lives at :mod:`orpheus.sn.boundary_realize` —
SN-specific. When a second functional realizer ships, the walker
will move to ``orpheus.geometry.boundary.realize`` and become
method-agnostic via
:meth:`BoundaryRealizerRegistry.get(method_name)` lookup. The leaf
dispatch (:class:`BoundaryTraceLaw` → realized op) changes; the
composer dispatch (:class:`LawSum` / :class:`LawScaled` → operator
composers) is method-agnostic and stays as-is.


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

The algebra: where the two semantics agree and where they diverge
-------------------------------------------------------------------

Decompose the angular ordinate set at a given face :math:`f` as

.. math::
   :label: ordinate-partition-inflow-outflow

   \{1, \ldots, N\} \;=\; I_f \,\sqcup\, O_f \,\sqcup\, T_f,

with :math:`I_f = \{n : \Omega_n \cdot \hat n_f < -\epsilon\}` the
inflow set, :math:`O_f = \{n : \Omega_n \cdot \hat n_f > +\epsilon\}`
the outflow set, and :math:`T_f` the (measure-zero in the continuum,
:math:`\epsilon`-band in the discretisation) tangential set. For
:math:`\psi_{\text{out}} \in \mathbb{R}^N` representing the trace
ordinate values at face :math:`f`, the two vacuum representations
are:

.. math::
   :label: vacuum-legacy-vs-trace-correct

   \mathrm{zeros\_all}(\psi_{\text{out}})[n] &= 0
       \qquad \forall\, n \in \{1, \ldots, N\}, \\[2pt]
   \mathrm{inflow\_mask}(\psi_{\text{out}})[n] &=
       \begin{cases}
           0 & n \in I_f, \\[2pt]
           \psi_{\text{out}}[n] & n \in O_f \cup T_f.
       \end{cases}

The two functions agree **on** :math:`I_f` (both give 0) and
**diverge** on :math:`O_f` (legacy gives 0, trace-correct gives
:math:`\psi_{\text{out}}[n]`). They diverge on :math:`T_f` too
in principle, but ORPHEUS's quadrature adapters carry every
tangential ordinate at :math:`\mu = 0` so :math:`\psi_{\text{out}}[n] = 0`
on :math:`T_f` for a properly-initialised flux — making the
divergence physically restricted to :math:`O_f`.

The §16A.5 production-relevant subset is **the inflow rows**.
Every SN sweep call site reads :math:`\psi_{\text{in}}[n]` only for
:math:`n \in I_f` — outflow rows are never consumed downstream.
The Wave 8 close-out audited all 13 ``bc.apply(...)`` sites in
:mod:`orpheus.sn.sweep` and :mod:`orpheus.sn.operator`:

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

No call site reads :math:`\psi_{\text{in}}[n]` for
:math:`n \in O_f`, so the two semantics produce **bit-identical
observable output** for every existing production consumer.
This is what makes the §16A.5 trace-correct representation a
safe semantic upgrade: it strictly extends the correctness
boundary (outflow rows are now correct under the new algebra
where they were ill-defined under the old one) without changing
any observable downstream value.

Post Issue #188 + #176 (2026-05-11) the **realizer path is
uniform** across every supported mesh — 1-D Cartesian, 1-D
spherical, 1-D cylindrical, and 2-D Cartesian. The curvilinear
sweeps now consume the realizer-routed
:class:`IncomingOrdinateMaskTensor` exactly like the slab sweep
does. **Empirical confirmation**: spherical 26/26 + cylindrical
25/25 + MMS curvilinear 2/2 xfail (pre-existing ERR-026) green
on C188.3 — the curvilinear sweeps were previously consuming the
zeros-all legacy body via the bound-quadrature shim, and now
consume the inflow-only-mask realizer output. Production result
is bit-identical on inflow rows (the only rows the sweep reads),
confirming the Wave 8 call-site audit empirically.

The Wave 6 snapshot harness gates this explicitly: the
``vacuum_lebedev17`` case compares **inflow-row outputs only**,
documenting the intentional semantic divergence in a comment
adjacent to the test case.

"Option a" (Wave 7) — historical context, retired Issue #186
-------------------------------------------------------------

The Wave 7 brief considered three migration strategies for the
2-arg legacy path. **Option (a)** ("vacuum-stays-legacy") landed:
:class:`VacuumInflow` carried a standalone
:meth:`apply(psi_out, quad)` whose body was
``np.zeros_like(psi_out)`` (the pre-§16A.5 zeros-all form), and
the realizer path produced the inflow-only-mask form via
:class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`.
The two paths agreed on inflow rows (the production-relevant
subset) and diverged on outflow rows.

Options (b) ("face-aware BC" — add ``face`` to the constructor)
and (c) ("combined ABC merge") were rejected because they would
have distorted the law's interface to better serve a transitional
path the refactor was retiring anyway.

**Status after Issue #186 (2026-05-11):** Option (a)'s standalone
``apply`` body has been **deleted**. :class:`VacuumInflow` (like
every other concrete BC) is a pure descriptor with no ``apply``
method; the only path to vacuum action is
:func:`realize_recursively` or :class:`SNBoundaryRealizer`
producing :class:`IncomingOrdinateMaskTensor` output. The
**two-paths-divergence is therefore eliminated by design** — there
is no longer a "second path" that could disagree with the realizer
path. The §16A.5 inflow-only-mask body is the **unique** vacuum
semantics in the post-#186 codebase.

This is the load-bearing architectural payoff of B3 + β2: the
documentation no longer needs to caveat which path you're on,
because there is only one path. See
:ref:`bc-trace-law-descriptor-model` for the design rationale.


.. _bc-trace-law-descriptor-model:

The trace-law descriptor model (Issue #186 / B3 + β2)
=====================================================

Issue #186 / Scope B3 + β2 (landed 2026-05-11 on branch
``feature/bc-curvilinear-realizer-cleanup``) is the architectural
**closure** of the BC trace-law refactor. It collapses the
remaining 2-arg ``apply`` affordance from the Wave-8/9 era into a
**pure-descriptor** contract:

* :class:`BoundaryTraceLaw` no longer inherits
  :class:`~orpheus.numerics.operator.LinearOperatorMixin`.
* The abstract :meth:`apply` method that the mixin used to provide
  is gone. So is ``apply_transpose``. So is the
  ``capabilities: ClassVar[frozenset[str]]`` advertisement.
* Every concrete BC (vacuum / reflective / white / albedo /
  periodic / prescribed-inflow) is now a **frozen dataclass**
  carrying only its parameters (axis, albedo, source, ...), its
  :attr:`kind` tag, its :attr:`creates_sweep_cycle` ``ClassVar``,
  its :attr:`geometry_map` / :attr:`response_kernel` /
  :attr:`source` property overrides, and the relevant
  :meth:`assert_*` invariants. **No** ``apply`` method on any
  concrete BC.
* The base class :class:`BoundaryTraceLaw` carries a **minimal
  algebra** that returns :class:`LawSum` / :class:`LawScaled`
  nodes — the descriptor-tree composition algebra documented at
  :ref:`bc-rank-n-algebra`. The dunders are: ``+``, ``-``, ``*``,
  ``/``, unary ``-``, plus their reflected variants. Each returns
  a new descriptor-tree node; none returns an operator.
* :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` is the
  **sole** bridge from descriptor to callable. There is no
  alternative path. Calling ``law.apply(psi)`` raises
  :class:`AttributeError` at runtime; a static type checker flags
  the call without running it.

The §16A.3 three-layer architecture (descriptor / realizer /
operator) is now **enforced by the type system**, not by convention.

What was tried and rejected before B3 + β2 landed
-------------------------------------------------

This subsection preserves the design history because the rejected
paths are the load-bearing intellectual content of the close-out
— future sessions asking "why does the realizer exist?" need to
see why every alternative failed.

**Option A** (Issue #176 / C176.3, ~2026-04). Concrete BC
:meth:`apply` methods kept a keyword-optional
``quadrature: AngularQuadrature | None = None`` parameter with
defensive errors:

* :class:`ReflectiveBoundary.apply` / :class:`WhiteBoundary.apply`
  raised :class:`BoundaryError` when ``quadrature is None``
  because their geometric / response operators needed the
  quadrature to construct themselves.
* :class:`VacuumInflow` / :class:`AlbedoBoundary` /
  :class:`PeriodicBoundary` / :class:`PrescribedInflow` accepted
  and ignored the ``quadrature`` parameter.

Option A was the **interim** landing — it preserved the
direct-call convenience pattern ("sketching code can write
``ReflectiveBoundary(axis='x').apply(psi, quad)``") while routing
production through the realizer. The C176.3 audit identified
three architectural costs that made Option A unsustainable:

1. **Asymmetric semantics on ``quadrature=None``.** Two BCs
   raised; four accepted-and-ignored. The behaviour was
   inconsistent across the law family and required per-BC
   documentation of "when is this method usable".
2. **Vacuum two-paths-divergence.** Direct ``VacuumInflow.apply(psi)``
   returned ``np.zeros_like(psi)`` (the pre-§16A.5 zeros-all body).
   The realizer-routed path returned
   :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
   output (the §16A.5 inflow-only mask). The two paths agreed on
   inflow rows (the production-relevant subset; see
   :ref:`bc-vacuum-semantic-correction`) but diverged on outflow
   rows. The divergence was harmless at every existing production
   call site but was a documentation-burden landmine for future
   adjoint-sensitivity consumers that read outflow rows.
3. **Liskov violation.** The abstract
   :meth:`BoundaryTraceLaw.apply(self, psi_out)` was strict 1-arg
   (post Issue #176 / C176.4). The concrete
   :meth:`apply(self, psi_out, quadrature=None)` was technically
   Liskov-substitutable (the optional parameter has a default), but
   calling ``bc.apply(psi)`` polymorphically on a
   :class:`BoundaryTraceLaw`-typed parameter could fail at runtime
   for Reflective/White — the type signature said "this works" but
   the runtime behaviour said "this raises". The static type
   system could not catch the failure because the contract was
   carried only in the docstring.

The B3 audit cataloged every remaining 2-arg ``bc.apply(psi,
quad)`` call site in production AND tests; none was load-bearing
for correctness. The Wave-6 snapshot regression carried legacy
halves (regenerated through the realizer path), the
:class:`PrescribedInflow` invariant tests ignored the
``quadrature`` argument anyway, and the realizer-vs-legacy
equivalence assertions could be replaced by hand-computed
expressions strictly stronger than legacy-agreement. The B3 sweep
rewrote every such site in one commit cycle (C-B3.7 through
C-B3.12) and deleted the ``apply`` methods.

**β1 interim** (sub-option within Issue #186 B3, considered but
not landed as a final state). Keep
:class:`LinearOperatorMixin` inheritance on
:class:`BoundaryTraceLaw` and drop only the abstract
:meth:`apply`. Rank-N composition would then build an
:class:`OperatorSum` of :class:`ScaledOperator`-wrapped *unrealized*
laws, and :func:`realize_recursively` would walk the resulting
Wave-0 composer tree. This β1 form was **algebraically equivalent**
to β2 but conflated two type families: the
:class:`OperatorSum` instance representing a not-yet-realized
descriptor tree was structurally identical to the
:class:`OperatorSum` representing an actual operator composition,
and only runtime inspection of the leaves could tell which was
which. β2 was preferred because the type system can then enforce
"this is a law, that is an operator" by static inspection — a
type checker rejects ``law_tree.apply(...)`` at the linter level
without ever running the code. β2 is more verbose (two new
dataclasses) but is the architecturally-checkable form.

**The vacuum two-paths-divergence is eliminated by design.** Under
B3 + β2 there is no "direct path" any more — every consumer must
realize the law before applying it, and realize-then-apply goes
through the inflow-only-mask path. The Wave-6 snapshot harness's
``vacuum_lebedev17`` case still pins inflow rows only (the
intentional §16A.5 divergence is still there in the algebra), but
no caller can accidentally invoke the pre-§16A.5 zeros-all path
because that path no longer exists in the code.

Empirical justification
-----------------------

The 18-test :mod:`tests.geometry.test_law_composition` suite pins
the descriptor-tree contract (foundation + L1 tests):

* Algebra closure on every dunder for every node-type pairing
  (laws × LawSum × LawScaled).
* :class:`LawScaled` constant folding
  (``2 * (3 * spec) == LawScaled(6, spec)``).
* The walker's exhaustive dispatch on the three node types and
  its :class:`TypeError` on unknown nodes.
* Walker value-correctness against hand-composed expectation
  (``realize_recursively(law_tree).apply(psi)``
  ``== 0.3 * realize(spec).apply(psi) + 0.7 * realize(white).apply(psi)``).
* The absence of ``apply`` on any descriptor-tree node
  (``not hasattr(tree, "apply")``).

The Wave-6 snapshot harness (now realizer-only — the legacy
halves are gone) verifies that the realized output is
bit-identical to the pre-B3 realizer-path output on every
non-mixed case, and identical up to the documented ULP tolerance
on the Marshak mixed case (the operator tree is structurally the
same; only the route to it changed).

Call-site contract
------------------

There is **one** way to call a boundary law's ``apply``:

.. code-block:: python

   from orpheus.geometry.boundary import ReflectiveBoundary
   from orpheus.sn.boundary_realizer import (
       SNBoundaryRealizer, SNMethodSpace,
   )

   law = ReflectiveBoundary(axis="x", albedo=0.5)
   ms = SNMethodSpace.minimal(quad)   # or .for_face(...) in production
   op = SNBoundaryRealizer().realize(law, ms)
   psi_in = op.apply(psi_out)

For descriptor-tree composition:

.. code-block:: python

   from orpheus.sn.boundary_realize import realize_recursively

   tree = 0.3 * ReflectiveBoundary(axis="x") + 0.7 * WhiteBoundary(
       axis="x", outward_sign=+1,
   )
   op_tree = realize_recursively(tree, ms)
   psi_in = op_tree.apply(psi_out)

No call site routes through a putative ``law.apply(psi)`` — that
method does not exist.


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


.. _bc-phase-d-two-bc-applies-per-matvec:

Phase D extension — two BC apply calls per curvilinear matvec
==============================================================

Phase C (:ref:`bc-trace-contract-respected-by-matvec`) established
that the SN curvilinear matvec applies the BC trace law **once per
matvec** at the boundary edge, consuming the WDD-propagated outflow
trace :math:`\gamma_+\psi`.  Phase D (Issue #168 Phase D, 2026-05-12,
landed on ``refactor/sn-operator-algebra``) extends the matvec with
a **second** BC apply call, used to build the Carlson coupled-pole
seed's ``bc_outer_value`` (see
:ref:`sn-phase-d-carlson-coupled-pole-sweep` in
:doc:`discrete_ordinates` for the math).  The §16A.3 affine trace
law contract is therefore exercised **twice per matvec** in the
post-Phase-D code path:

.. list-table:: BC apply call sequence inside one curvilinear matvec
   :header-rows: 1
   :widths: 12 28 30 30

   * - Order
     - Caller / purpose
     - Input shape & meaning
     - Output use
   * - **#1**
     - Phase D Carlson context build
       (:func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
       / ``_cylindrical`` early in the call)
     - ``(N, ng)`` — outer-cell cell-centred :math:`\psi` (NOT the
       face trace; a first-order proxy used only to construct the
       linear-in-:math:`\psi` ``bc_outer_value`` scalar at
       :math:`\mu = -1`)
     - Extract the most-inward-ordinate row; scalar feeds into
       :class:`CarlsonSweepContext.bc_outer_value`
   * - **#2**
     - Phase C BC trace law application (at the boundary edge after
       the WDD sweep completes)
     - ``(N, ng)`` — WDD-propagated outflow face values
       :math:`\gamma_+\psi` (the §16A.3 contract input)
     - Fill the inflow rows used as
       :math:`\psi^{\text{face}}_{\text{in}}` for the inward sweep
       phase

The two calls are **structurally distinct**:

* **Call #1** is a Phase D-specific use of the BC operator as a
  *linear-in-ψ* construction of the inward-zero-weight ordinate's
  outer-face flux.  For vacuum BC the
  :class:`IncomingOrdinateMaskTensor` zeroes the inflow ordinate,
  giving ``bc_outer_value = 0``; for reflective BC the
  :class:`PermutationOperator` mirrors outgoing :math:`\leftrightarrow`
  incoming, giving ``bc_outer_value = ψ_cell[N-1]`` (i.e. the
  cell-centred outer-cell value).  Both behaviours preserve operator
  linearity in the input :math:`\psi`.

  The input shape ``(N, ng)`` for cell-centres is a structural
  proxy: the BC trace law expects a trace ``(N, ng)`` shape, and
  feeding it the outer cell-centre array IS the right shape for
  Call #1's linear extraction.  The §16A.3 contract is not
  literally honoured here in the *interpretation* sense (the input
  is not a face trace), but the resulting scalar is linear in the
  matvec's input :math:`\psi` and gives the correct value at
  :math:`\mu = -1` on the only configurations the apply matvec
  cares about (reflective + flat :math:`\psi` :math:`\rightarrow` C,
  vacuum :math:`\rightarrow` 0).  This is the **principled
  shortcut**: a linearly-compatible scalar extraction whose values
  match the canonical inward-zero-weight ordinate's flux on the
  load-bearing test configurations.

* **Call #2** is the canonical Phase C use — the BC trace law
  consumes the **actual** WDD-propagated outflow face trace and
  produces the inflow trace per the §16A.3 affine-bc-form
  contract.  This is the call the Phase C
  :ref:`Gate 1.5 <bc-trace-contract-respected-by-matvec>` test
  pins.

Capture-and-compare Gate 1.5 strengthening
------------------------------------------

The pre-Phase-D Gate 1.5 test was a "round-trip" check: invoke
``bc.realize().apply(...)`` independently and compare against the
matvec's observable output.  Phase D **strengthens** the gate to a
capture-and-compare check that audits the *exact* value the matvec
passes into the BC trace law — necessary because the matvec now
calls :meth:`bc.apply` twice and the test must locate Call #2 (the
§16A.3 call) unambiguously.

The Phase D test
:func:`tests.sn.test_phase_c_gates.test_bc_trace_contract_capture_and_compare_sphere`
(parametrised over ``vacuum`` and ``reflective``):

#. Monkey-patches :meth:`sn_mesh.bc_right.apply` with a recorder
   wrapper that appends every input array passed to it during one
   matvec call.
#. Independently reconstructs the WDD-propagated outflow trace via
   a reference implementation
   (``_outflow_at_boundary_for_sphere(sn_mesh, sig_t, psi_input)``).
#. **Locates Call #2** by matching shape ``(N, ng)`` AND content
   (the captured input that bit-matches the independent reference
   IS the Phase C call).
#. Asserts the located input matches the reference to
   ``rtol=1e-14`` — bit-equal up to FP non-associativity.

The test is **foundation-tagged** (``@pytest.mark.foundation``)
because it pins a software invariant about the matvec's
two-application sequence, not a math claim about the BC operator.
The matching strategy by *both* shape and content protects against a
future regression that adds a third BC apply with the right shape
but wrong content — the test would still locate the canonical Phase
C call provided its content matches the WDD reference.

Both ``vacuum`` and ``reflective`` parametrised cases pass.  The
``vacuum`` case is the load-bearing check because Call #1 produces
non-trivial output under vacuum (the
:class:`IncomingOrdinateMaskTensor` zeroes inflow ordinates, so the
extracted ``bc_outer_value`` is zero — but the **input** to Call #1
is the outer cell-centre value which is **not** zero on a non-trivial
:math:`\psi`).  Locating Call #2 unambiguously requires content
matching, not just shape matching.


.. _bc-phase-f-three-bc-applies-per-sweep-iteration:

Phase F extension — BC applies in the SI sweep path
====================================================

Phase D (Issue #168 Phase D, :ref:`bc-phase-d-two-bc-applies-per-matvec`
above) instituted the *two BC apply calls per curvilinear matvec*
contract on the apply-matvec path (the within-group operator,
:class:`~orpheus.sn.operator.InvertibleOperator` post-Depth-B; the
matvec lives at :func:`~orpheus.sn.operator.transport_operator_matvec_unified`).  Phase F
(Issue #168 Phase F, 2026-05-12, also landed on
``refactor/sn-operator-algebra``) propagates the same pattern to the
**SI/sweep path** (:func:`~orpheus.sn.sweep.transport_sweep` →
:func:`~orpheus.sn.sweep._sweep_1d_spherical` /
``_sweep_1d_cylindrical``).  See
:ref:`sn-phase-f-carlson-sweep-path-backport` in
:doc:`discrete_ordinates` for the math and the
twin-path-fix-incompleteness anti-pattern that motivated the
backport.

The post-Phase-F SI sweep iteration applies the BC operator in
**three** distinct invocations per sweep call — one for the Phase F
Carlson seed, plus :math:`N_{\text{inward ordinates}}` legacy inflow
applications inside the per-ordinate loop (which are not
fundamentally new — they predated Phase D and Phase F).  The
load-bearing addition is the **Phase F seed call**:

.. list-table:: BC apply call sequence inside one SI sphere sweep
   :header-rows: 1
   :widths: 12 28 30 30

   * - Order
     - Caller / purpose
     - Input shape & meaning
     - Output use
   * - **#1**
     - Phase F Carlson seed
       (:func:`~orpheus.sn.sweep._sweep_1d_spherical` early in
       the call, before the per-ordinate loop)
     - ``(N, ng)`` — persistent outer-face outflow buffer
       ``bc_outer`` carrying the previous outward sweep's
       outgoing flux per ordinate (zero on the first SI
       iteration)
     - Extract the most-inward-ordinate row; scalar feeds into
       :func:`carlson_inward_sweep_from_source` as
       ``bc_outer_value`` (the seed for Hébert (3.434)–(3.435)
       at the outer face)
   * - **#2 … #N₋**
     - Per-inward-ordinate inflow read inside the
       per-ordinate sweep loop
     - ``(N, ng)`` — same ``bc_outer`` buffer (re-read each
       iteration, since intervening outward ordinates may have
       updated it)
     - Read the inflow row ``psi_in_full[n]`` for the current
       inward ordinate :math:`\mu_n < 0` as the spatial
       sweep's incoming flux

**Comparison with the apply-matvec twin** (per
:ref:`bc-phase-d-two-bc-applies-per-matvec`):

* The apply matvec consolidates its inflow logic into the **single
  Phase C trace law call** at the boundary edge — the BC operator
  is invoked once on the WDD-propagated outflow trace, producing
  the inflow trace per the §16A.3 affine-bc-form contract for ALL
  ordinates simultaneously.
* The SI sweep, by contrast, processes ordinates **sequentially**;
  each inward ordinate independently reads its inflow row from
  the persistent ``bc_outer`` buffer.  The per-ordinate apply
  calls are **not** a §16A.3 trace-law application of the same
  semantic kind — they are *consumer reads* of an already-updated
  buffer.  The Phase F seed call (Call #1) is the only call that
  semantically mirrors Phase D's apply-matvec Call #1 (a
  Phase-specific use of the BC operator as a linear-in-:math:`\psi`
  construction of the inward-zero-weight ordinate's outer-face
  flux).

The Phase F seed call's role is exactly analogous to the
apply-matvec's Phase D Call #1 (per the
:ref:`bc-phase-d-two-bc-applies-per-matvec` table): a
linear-in-:math:`\psi` extraction of the inward-zero-weight
ordinate's outer-face flux.  Vacuum BC zeros it; reflective BC
yields the most-inward ordinate's mirrored outflow.  In both
contexts the BC operator's **linearity** is the load-bearing
contract — the Phase F helper
:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`
must remain linear in the input to preserve the SI loop's
fixed-point convergence properties.

The cylindrical path has the analogous structure with a
**per-:math:`\mu`-level** seed call: one BC apply per level,
each invocation extracting the level-specific most-inward
ordinate's row from the same persistent ``bc_outer_cyl``
buffer.  The
:func:`~orpheus.sn.sweep._sweep_1d_cylindrical` body invokes
the BC operator ``len(quad.level_indices) + N_inward`` times
total per sweep — once per :math:`\mu`-level for the Phase F
seed, plus once per inward ordinate inside each level's
azimuthal loop.

Phase F leaves the §16A.3 BC trace contract semantics
**unchanged**: the SI sweep path does not call the §16A.3 trace
law application that the apply matvec uses at the boundary
edge — the SI sweep updates ``bc_outer`` directly from each
outward ordinate's last visit (line ~593 of
:file:`orpheus/sn/sweep.py`), then the next inward ordinate
reads it via ``apply``.  The semantic contract is the same
(BC operator maps outflow trace :math:`\gamma_+\psi` to inflow
trace :math:`\gamma_-\psi`), but the *invocation pattern* is
per-ordinate sequential rather than once-per-matvec
collective.

No new Gate 1.5 test variant is needed for the Phase F seed
call.  The Phase F bit-identity test
:func:`tests.sn.spatial.test_sweep_vs_apply_consistency`
(57 foundation tests) pins that the sweep-path's
``bc_outer_value`` extraction matches the apply-path's Phase D
Call #1 result on every test configuration — the structural
invariant that the two paths' Carlson seeds agree on matching
inputs subsumes the BC-apply-input pinning.


.. _bc-curvilinear-realizer-unification:

Curvilinear realizer unification (Issue #188 + #176 close-out)
==============================================================

The pre-cleanup architecture carried a **Cartesian / curvilinear
split** at :meth:`SNMesh._resolve_one`: the slab and 2-D Cartesian
paths constructed an :class:`InflowTraceSpace` and routed the BC
through :class:`SNBoundaryRealizer`, while spherical and
cylindrical ``Mesh1D`` bypassed the realizer entirely because
:meth:`InflowTraceSpace.from_mesh_and_quadrature` raised
:class:`NotImplementedError` on those coord systems. The bypass
wrapped the bare law instance in
:class:`_BoundBoundaryOperator(law, quadrature=self.quad)`, a
dual-mode shim whose ``apply(psi)`` forwarded
``law.apply(psi, bound_quad)`` to the legacy 2-arg body.

**Why the split existed.** The Wave 2 ``InflowTraceSpace``
factory deferred curvilinear support because no curvilinear-Krylov
consumer at the time needed the per-face mask — the curvilinear
sweep paths computed the inflow / outflow predicate on the fly
inside the inner loop. The deferral was load-bearing for the
shim's dual-mode design: ``quadrature=None`` (Cartesian, wraps a
realized 1-arg op) and ``quadrature=<set>`` (curvilinear, wraps a
legacy 2-arg BC).

**Why the split has been removed.** Issue #188 / C188.1+C188.2
discovered that all three :class:`Mesh1D` coord systems
(``CARTESIAN`` / ``SPHERICAL`` / ``CYLINDRICAL``) share the same
``("left", "right")`` face structure with the radial axis as the
outward normal and the :class:`~orpheus.sn.quadrature.GaussLegendre1D`
adapter's :math:`\mu_x` as the direction cosine along that axis —
identically to the slab case. The mask predicate
:math:`\Omega \cdot \hat n < -\epsilon` therefore applies
unchanged. The factory's curvilinear guard was lifted; only 2-D
cylindrical :class:`Mesh2D` (which has no SN sweep in ORPHEUS
today) continues to raise :class:`NotImplementedError`.

Issue #188 / C188.3 then collapsed :meth:`SNMesh._resolve_one`
to a single path: every supported mesh (1-D Cartesian / spherical
/ cylindrical + 2-D Cartesian) builds an
:class:`SNMethodSpace.for_face` and routes through
:class:`SNBoundaryRealizer`. Issue #176 / C176.1 then dropped the
``quadrature=`` kwarg from the shim because no
production-issued shim carried ``_quadrature is not None`` after
C188.3. Issue #176 / C176.3+C176.4 then trimmed the concrete BC
``apply`` signatures to the Option-A interim (keyword-optional
``quadrature=None`` with defensive errors). Issue #186 / B3 + β2
then retired Option A entirely — every concrete BC ``apply``
method was deleted; see
:ref:`bc-trace-law-descriptor-model` for the retrospective.

The architectural sequence is therefore:

* **Issue #188 unblocks Issue #176.** The shim's dual mode existed
  ONLY because curvilinear :class:`InflowTraceSpace` could not be
  constructed. Once #188 lifted that, #176's "drop the 2-arg form"
  cleanup became possible without breaking curvilinear sweeps.
* **Issue #176 unblocks Issue #186.** The Option-A interim was
  necessary because dropping ``apply`` outright before the
  curvilinear sweeps consumed realizer output (#188) and the test
  fleet migrated to the realizer-routed contract (#176 / C176.5
  cleanup commits) would have broken curvilinear regression.
  Once those landed, the descriptor cleanup (#186 / B3 + β2)
  became the next step on the architectural ladder.

.. _bc-1d-y-placeholder-design:

1-D y-face placeholders
-----------------------

A :class:`Mesh1D` mesh has no physical ``y`` faces, but
:class:`SNMesh` exposes ``bc_ymin`` / ``bc_ymax`` attributes
uniformly with the 2-D path so cross-dimensional code can read
them without coord-system gating. After Issue #188, those
placeholders route through the realizer with
:meth:`SNMethodSpace.minimal(quad)` and a
:class:`ReflectiveBoundary(axis="y")` law. Three pre-investigated
facts make this safe:

1. **``GaussLegendre1D.reflection_index("y")`` returns the
   identity permutation.** Every ordinate is its own ``y``-reflected
   partner because :math:`\mu_y \equiv 0` on the 1-D quadrature
   (only ``mu_x`` is populated). The realized
   :class:`PermutationOperator` is therefore a no-op.
2. **The SN realizer's reflective branch reads only ``law.axis``
   and ``quad.reflection_index(axis)``** — never the per-face
   :attr:`inflow_indices`. The minimal method space carries the
   quadrature only; the absence of an :class:`InflowTraceSpace`
   for the ``y`` faces is not a problem for this branch.
3. **The 1-D trace space has
   ``face_names == ("left", "right")``** — calling
   :meth:`SNMethodSpace.for_face(face="ymin", inflow_trace=...)`
   would raise. ``.minimal(quad)`` is therefore required.

Production consumers in
:meth:`solution_to_angular_flux` only invoke ``bc_ymin.apply(...)``
on 1-D meshes when ``curvature is None`` (the slab case), and the
inner ``if mu_y[n] > 1e-15`` filter discards every ordinate on
GL1D — so the no-op behaviour is observably correct.

Closure
-------

Closed by branch ``feature/bc-curvilinear-realizer-cleanup``
(2026-05-11). Three GitHub issues converged on this branch:

* **Issue #188** — curvilinear :class:`InflowTraceSpace` support
  (commits ``9cf2b0a`` + ``17067d5``). Lifted the
  :class:`NotImplementedError` guard on spherical / cylindrical
  Mesh1D so every supported mesh can build a per-face inflow mask.
* **Issue #176** — drop 2-arg ``apply`` + simplify shim (commits
  ``cf29ce4`` + ``a4a43c2`` + ``913e501`` + ``188bf9a``). Collapsed
  the dual-mode shim into a strict 1-arg passthrough; landed the
  Option-A interim with keyword-optional ``quadrature=None`` on
  the concrete laws.
* **Issue #186 (B3 + β2)** — pure-descriptor cleanup (commits
  ``f71a32c`` + ``da414eb`` + ``89d09a4`` + ``633cc69`` +
  ``bb674da`` + the test-migration trail). Retired the Option-A
  ``apply`` methods, dropped :class:`LinearOperatorMixin`
  inheritance from :class:`BoundaryTraceLaw`, and formalised the
  descriptor-tree composition algebra via the new
  :class:`LawSum` / :class:`LawScaled` types. The architectural
  sequence is therefore Issue #188 → #176 → #186: each step
  unblocked the next.

The :class:`_BoundBoundaryOperator` shim survives because the
``kind``-string tag is load-bearing for the BC-resolution
diagnostic and several legacy ``sn_mesh.bc_left ==
"vacuum"``-style test sites; the dual-mode bound-quadrature
backing is gone (#176), and the ``*_extra, **_kw`` swallow is
gone (#186 / C-B3.4). Every supported mesh consumes a strict
1-arg :class:`LinearOperator` produced by
:class:`SNBoundaryRealizer` for single BCs, or by
:func:`realize_recursively` for rank-N descriptor trees.

Plan documents:

* ``.claude/plans/transient-giggling-cake.md`` — the foundational
  12-wave BC trace-law refactor plan (Waves 0–12 close-out
  documented at :ref:`theory-boundary-conditions`).
* ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` —
  the #188 + #176 cleanup plan (Option-A landing).
* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` — the
  Issue #186 B3 + β2 cleanup plan (descriptor-model landing).

Grand Report v3 §16A.3 (the three-layer architecture) is now
**enforced by the type system** — descriptors have no ``apply``,
operators do. Grand Report v3 §16A.5 (the trace-correct vacuum
representation) is uniform across coord systems and the legacy
zeros-all path no longer exists.


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
6. **Cartesian-vs-curvilinear bypass in
   :meth:`SNMesh._resolve_one` + dual-mode shim.** Pre Issue #188
   shape: curvilinear ``Mesh1D`` bypassed the realizer and wrapped
   the bare 2-arg law in
   ``_BoundBoundaryOperator(law, quadrature=self.quad)``, while
   Cartesian routed through the realizer with
   ``_BoundBoundaryOperator(realized)``. **Retired Issue #188 +
   #176**; documented here because the seductive trap is to
   "preserve flexibility" by keeping the dual mode after the
   curvilinear deferral lifts. The right move is to delete the
   bypass and consolidate on one path — see
   :ref:`bc-curvilinear-realizer-unification`.
7. **Option A (keyword-optional ``quadrature=None`` on the
   concrete laws).** Landed Issue #176 / C176.3 as the interim
   form; **retired Issue #186 / B3 + β2** in favour of the
   pure-descriptor model (no ``apply`` on any law). The
   architectural costs of Option A (asymmetric semantics on
   ``quadrature=None``, vacuum two-paths-divergence, Liskov
   violation under polymorphic typing) made it unsustainable as
   the long-term contract; the interim was kept only long enough
   to land curvilinear realizer unification (Issue #188 first)
   before the descriptor cleanup could ship. See
   :ref:`bc-trace-law-descriptor-model` for the full retrospective.
8. **Calling ``apply`` on a raw BC descriptor.** Under the
   pre-#186 contract this either worked (with surprising
   semantics — see Option A entry above) or raised
   :class:`BoundaryError`. Under post-#186 it's a **static type
   error** — :class:`BoundaryTraceLaw` has no :meth:`apply`
   method on the class, and neither do :class:`LawSum` /
   :class:`LawScaled`. The correct contract is
   ``SNBoundaryRealizer().realize(law, ms).apply(psi)`` for a
   single BC, or ``realize_recursively(tree, ms).apply(psi)`` for
   a descriptor tree. The realizer is the **sole** bridge; the
   §16A.3 three-layer split is enforced by the type system.
9. **In-tree Wave-0 operator algebra over unrealized
   :class:`BoundaryTraceLaw` instances (β1 form).** Considered as
   the rank-N composition mechanism during Issue #186 B3
   exploration. Rejected in favour of the separate-type-family
   approach (β2 / :class:`LawSum` / :class:`LawScaled`). β1
   produced :class:`OperatorSum` trees whose leaves were laws,
   not operators — the type checker could not distinguish a
   not-yet-realized "operator" from a real operator, and calling
   :meth:`apply` on the β1 tree raised at the leaf realization
   step. β2 makes the law-vs-operator distinction inspectable
   statically: :class:`LawSum` has no :meth:`apply` method, full
   stop. See :ref:`bc-rank-n-algebra` for the detailed
   comparison.


References
==========

* Grand Report v3 §16, §16A.1–5 (affine boundary form + trace
  structure), §16A.10 (sparse trace primitives), §16A.11 (dual
  registry), §16A.12 + §27.6 (universal invariants), §26A.4
  (named-error catalog). Source: ``.claude/plans/neutron_transport_grand_report_v3.md``.
* The 12-wave refactor plan:
  ``.claude/plans/transient-giggling-cake.md``.
* The post-Wave-12 cleanup plan (Issue #188 + #176):
  ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md``.
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
