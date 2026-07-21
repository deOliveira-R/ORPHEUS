.. _theory-discrete-ordinates:

==========================================
Discrete Ordinates Method (S\ :sub:`N`)
==========================================

.. contents:: Contents
   :local:
   :depth: 3


.. Machine header — the ``nexus-meta`` schema for this page.  Ingestion is
.. PENDING nexus#1 Phase 2: the ``nexus-meta`` directive is NOT yet
.. registered, so the schema is rendered here as a collapsed sphinx-design
.. dropdown and machine-consumed later.

.. dropdown:: Machine header — ``nexus-meta`` schema (module · operators · conventions · invariants)
   :color: muted

   .. code-block:: yaml

      module: sn
      method: discrete-ordinates
      aliases: [SN, discrete ordinates, Sₙ, transport sweep, ordinate method]
      governing_equation: "A ψ = (1/k) F ψ  [eigenvalue] ;  A ψ = q  [fixed source]"
      operators:
        L: streaming, BULK only (Ω·∇) — the boundary law is the sibling B, NOT folded into L
        C: collision / removal (Σ_t)
        S: scattering in-scatter gain (Σ_s0ᵀ φ + anisotropic moments)
        B: boundary law as a first-class SIBLING operator (reflective / vacuum / white trace), every geometry
        F: fission production (χ ⊗ νΣ_f, rank-1 dyad)
      composites:
        A: "L + C - S - B — the within-group loss operator; the Krylov driver applies it"
        (L+C): "lower-triangular under the upwind cell ordering; (L+C)⁻¹ IS the transport sweep"
      key_types: [AngularFlux, SNMesh, HarmonicMomentField, SweepDependencyGraph]
      entry_points:                    # qualnames; Nexus links via implements edges
        - orpheus.sn.solver.solve_sn
        - orpheus.sn.solver.SNSolver
      conventions:
        sign: "μ>0 is +x / outward-radial outflow; inward inflow at the left / outer boundary"
        scattering: "Mixture.SigS[l][g_from, g_to]; the in-scatter source uses the TRANSPOSE  Q = SigSᵀ @ φ"
        diamond_difference: "ψᵃ = (1+β)ψ_out − β ψ_in; Morel–Montry sets β = 0 (Bailey–Morel–Chang 2010 Eq. 43, unique exact-on-linear-in-μ)"
        quadrature_norm: "GL weights sum to 2; Lebedev / level-symmetric / product sum to 4π; moments carry NO 4π prefactor (Σw normalisation)"
        layout: "angular flux (N, ng, nx, ny); scalar flux and per-cell XS (ng, nx, ny); 1-D keeps ny=1 (singleton, not squeezed)"
        group_ordering: "fast → thermal; downscatter makes SigS upper-triangular"
        starting_direction: "curvilinear half-angle seed ψ_{1/2} is first-class typed state (System B), marched directly (Issue #282 route (a)); only levels with first-ordinate raw τ ∈ (0,1) carry the block (R12a)"
      invariants:
        - "particle balance PER ORDINATE (flat-flux residual = 0) — the strong check, NOT the telescoped scalar balance"
        - "sweep ≡ matvec (one loss representation, two applications: solve vs residual)"
        - "α redistribution dome ≥ 0 (negative → NaN / overflow)"
      depends_on: [transport_methods, operator_algebra, spherical_harmonics, frame]
      verification: [L0, L1, L2]       # authored claim; cross-checked vs the Verification slice (§ below)


.. _sn-synopsis:

Synopsis
========

The discrete ordinates (S\ :sub:`N`) method solves the
:ref:`multi-group eigenvalue problem <mg-eigenvalue-problem>` in
integro-differential form by discretising the direction variable
:math:`\hat{\Omega}` into a finite ordinate set
:math:`\{(\hat{\Omega}_m, w_m)\}`, **retaining the angular flux**
:math:`\psi(\mathbf{r}, \hat{\Omega}, E)` rather than collapsing to the
scalar flux (contrast the collision-probability integral form).  It resolves
streaming, anisotropic scattering, and interface angular current directly.
ORPHEUS supports three coordinate systems under one balance framework:
**Cartesian** (slab / 2-D, no inter-ordinate coupling), **spherical** (1-D
radial, a single :math:`\alpha`-redistribution dome coupling all ordinates in
:math:`\mu`), and **cylindrical** (1-D radial, an independent :math:`\alpha`
dome per :math:`\mu`-level).  All three share a geometry factor
:math:`\Delta A / w` that guarantees per-ordinate flat-flux consistency; the
curvilinear formulation follows [BaileyMorelChang2010]_ (Eq. 43, the
Morel–Montry angular-closure weight — unique exact-on-linear-in-:math:`\mu`),
the general framework [LewisMiller1984]_, and the angular discretisation
[CaseZweifel1967]_ / [Hebert2009]_ (§3.9.4).

The solver is posed as an **operator algebra** over five operators: streaming
:math:`L` (bulk :math:`\hat{\Omega}\cdot\nabla`), collision / removal
:math:`C`, the scattering gain :math:`S`, the boundary law :math:`B` — a
first-class **sibling** operator, *not* folded into :math:`L` — and the rank-1
fission dyad :math:`F`.  They compose the within-group loss operator
:math:`A = L + C - S - B`, so the eigenvalue problem is
:math:`A\,\psi = \tfrac{1}{k}\,F\,\psi` (fixed source: :math:`A\,\psi = q`).
The sub-composite :math:`(L+C)` is lower-triangular under the upwind cell
ordering, which is exactly why :math:`(L+C)^{-1}` **is** the transport sweep
(:doc:`/theory/methods/sn/loss_representation`).  :class:`SNSolver` satisfies
the :class:`~numerics.eigenvalue.EigenvalueSolver` protocol and
:func:`solve_sn` returns an :class:`SNResult`.  Because the protocol places the
scattering source *inside* ``solve_fixed_source``, the inner source iteration
(in-scatter + anisotropic convergence) stays encapsulated in the SN sweep,
while the outer :func:`~numerics.eigenvalue.power_iteration` loop is the one
shared by CP, MoC, diffusion, and the homogeneous solver (see
:doc:`/api/numerics` for the protocol contract).

The spatial closure is **diamond difference** with the Morel–Montry weight
(:math:`\beta = 0`); the per-ordinate discrete transport is the **sweep**,
which is byte-identical to the loss-operator **matvec** — one loss
representation, two applications (``solve`` vs residual).  Both the 1-D scan
(:meth:`~orpheus.sn.loss_representation.CumprodScan.sweep`) and the 2-D
wavefront sweep
(:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph`,
per-octant batched dispatch over a mesh-time-precomputed DAG) are **bare**:
the reflective coupling :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}`
rides as a sibling :math:`-B` source term rather than a re-applied boundary
condition (:ref:`bare-sweep-extraction`, and the canonical algebra
:ref:`bc-extraction` in :doc:`/theory/foundations/boundary_conditions`).  2-D Cartesian eigenvalue
problems solve through **both** inner drivers —
:class:`~orpheus.numerics.iteration.SourceIteration` (the geometry-agnostic
default) and :class:`~orpheus.numerics.iteration.KrylovAcceleration` —
verified SI ≡ Krylov ≡ closed-form :math:`k_\infty`
(:ref:`bc-extraction-2d-si-krylov-twin`).  Interior cell-face fluxes are typed
as an interior 1-cochain :math:`C^1_{\rm int}` carried in the rolling front
(:ref:`wavefront-flux-cochain`).

Curvilinear redistribution — the geometric :math:`\alpha`-dome, distinct from
Legendre :math:`P_1^+` scattering anisotropy — and its half-angle
**starting-direction seed** are first-class typed state.  The Issue #282 route
(a) design marches :math:`\psi_{1/2}` directly from the true :math:`q_{1/2}`
source through System B's named resolvent
:meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
(a single-pass exact inverse on the *full* Legendre fold), only on levels
whose first-ordinate raw :math:`\tau \in (0,1)` (**R12a**); see
:ref:`sn-282-direct-starting-direction-solve`.  The "#229 floor" resolution
— three distinct curvilinear errors separated by a
volume-weighted-:math:`L_2`-vs-:math:`L_\infty` norm difference — is
:ref:`sn-curvilinear-aniso-norm-reconciliation`.  The storage / index /
scattering / closure **conventions** and the load-bearing **invariants** are
captured structurally in the machine header above and cross-linked below;
verification is L0–L2 plus semi-analytical (Sood, Case singular-eigenfunction)
benchmarks; the traps that hide solver bugs behind green tests are collected in
:ref:`sn-gotchas`.

.. admonition:: Conventions

   - Scattering matrix: :ref:`scattering-matrix-convention` — ``SigS[g_from, g_to]``, source uses transpose
   - **Storage layout**: :ref:`theory-sn-index-convention` — ``(N, ng, nx, ny)`` for ψ, ``(ng, nx, ny)`` for φ / σ
   - Multi-group balance: :eq:`mg-balance` in :ref:`theory-homogeneous`
   - Cross sections: :ref:`theory-cross-section-data`
   - Verification: :ref:`synthetic-xs-library` — regions A/B/C/D
   - Eigenvalue: :ref:`power-iteration-algorithm` shared with all deterministic solvers


Architecture
============

.. note:: **Implementation map — automation pending.**  The auto-generated
   Nexus filtered flow-graph figure (root symbol + traversal depth →
   graphviz) that will head this section is blocked on the nexus#20
   flow-graph directive; until it ships, the architecture below is
   **hand-authored**.  See :doc:`/api/numerics` for the live
   operator-protocol surface and :doc:`/theory/foundations/operator_algebra` for the algebra.

Two-Layer Mesh Pattern
----------------------

The S\ :sub:`N` solver follows the same two-layer pattern as the CP
solver.  This pattern (base :class:`~geometry.mesh.Mesh1D` + augmented
mesh) is shared with :ref:`theory-collision-probability` and
:ref:`theory-method-of-characteristics`.

1. **Base geometry** --- :class:`~geometry.mesh.Mesh1D` or
   :class:`~geometry.mesh.Mesh2D` stores cell edges, material IDs,
   coordinate system, and **boundary condition declarations**.
   Each face carries an optional :class:`~geometry.mesh.BC` field
   (``bc_left``/``bc_right`` for 1-D;
   ``bc_xmin``/``bc_xmax``/``bc_ymin``/``bc_ymax`` for 2-D).
   When ``None`` (the default), the solver applies its own default
   --- for the SN solver, that default is reflective.
   See :ref:`boundary-conditions` for details.

2. **Augmented geometry** --- :class:`SNMesh` pairs the spatial mesh
   with an angular quadrature, precomputing the coordinate-specific
   streaming stencil.  Its **primary representation is the per-axis
   tuple** :attr:`SNMesh.axes` (the SN phase space factors as a tensor
   product of per-axis 1-D meshes): a legacy ``Mesh1D`` / ``Mesh2D`` is
   converted to axes **once** at the inbound boundary, and
   :meth:`SNMesh.from_axes` stores the caller's tuple verbatim. After
   C5 (:ref:`sn-axis-primary-c5`) :attr:`SNMesh.mesh` is *inbound
   provenance only* — ``None`` for an axis-native :math:`d \ge 3` mesh,
   which carries no legacy mesh at all.  It also **resolves boundary
   conditions**: each ``BC`` tag
   on the mesh is looked up in :attr:`SNMesh.BOUNDARY_OPERATOR_REGISTRY` and converted
   to a validated kind string (``"vacuum"`` or ``"reflective"``)
   stored in the face-name-keyed :attr:`SNMesh.bc` dict
   (``sn_mesh.bc["xmin"]``, ``sn_mesh.bc["xmax"]``, ... — the dict
   keys are the mesh's true boundary faces; see
   :ref:`bc-face-name-carve`).
   The sweep reads these resolved strings directly --- it never
   inspects the raw :class:`~geometry.mesh.BC` objects.  Precomputed
   stencil contents per coordinate system:

   - **Cartesian**: one per-axis array ``streaming(a)[n,i] =
     2|mu_a|/da[i]`` for every axis ``a < ndim`` (built over
     ``range(ndim)`` from ``quad.axis_cosines(a)`` since C3.6 ---
     no hand-listed x/y pair, no phantom axis on a slab) --- the
     diamond-difference denominator terms, precomputed to avoid
     per-cell division in the sweep hot loop.
   - **Spherical**: ``face_areas`` (:math:`4\pi r^2`), ``delta_A``,
     ``alpha_half`` (angular redistribution dome), and
     ``redist_dAw`` (:math:`\Delta A_i / w_n`, shape ``(nx, N)``).
   - **Cylindrical**: ``face_areas`` (:math:`2\pi r`), ``delta_A``,
     ``alpha_per_level`` (per-level redistribution domes), and
     ``redist_dAw_per_level`` (list of ``(nx, M)`` arrays).

   The Morel--Montry angular weight :math:`\tau` is **not** a factory
   output: it is owned by the angular closure
   (:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`),
   since the geometry-side producer was retired in Issue #236 Phase 2
   Step C (see :ref:`sn-tau-c-on-cellvisit-live`).  The geometry
   factories carry **geometry only** — face areas, the
   :math:`\alpha`-dome, the redistribution factor :math:`\Delta A/w`,
   and the level starting-direction edge :math:`\mu_{1/2}`.

3. **Solver** --- :func:`solve_sn` creates an ``SNMesh``, builds the
   ``SNSolver``, and runs power iteration. At :math:`d \le 2` the input
   is a ``Mesh1D`` / ``Mesh2D``; at :math:`d = 3` the input is the
   **axes tuple** itself (the only 3-D entry — there is no ``Mesh3D``;
   see :ref:`sn-c5-3d-admission`).

.. code-block:: text

   Mesh1D / Mesh2D (d<=2)   OR   axes tuple (d=3, axis-native)
       |                              |
       |  axes_from_legacy_mesh       |  (stored verbatim)
       +------------+-----------------+
                    v
   SNMesh.axes  (PRIMARY)  -->  SNMesh (stencil + quadrature
                                + alpha coefficients + resolved BCs)
                    |
                    v
   solve_sn() --> SNResult

Quadrature Dispatch
-------------------

The geometry-and-quadrature dispatch is a first-class polymorphism:
:func:`~orpheus.sn.loss_representation.default_for` selects the
:class:`~orpheus.sn.loss_representation.LossRepresentation` whose declared
``supports`` predicate matches the mesh — the 1-D chain scan
(:class:`~orpheus.sn.loss_representation.CumprodScan`, any geometry) or the
multi-D anti-hyperplane wavefront — and the operator then calls it
branchlessly.  (This replaced the pre-carve procedural branch on the
``SNMesh.curvature`` string tags and the since-retired operator-free
``transport_sweep`` entry; the full carve is documented in
:doc:`/theory/methods/sn/loss_representation`.)  Boundary conditions are **not** passed as
a parameter to the sweep --- it reads the resolved BC kind strings
directly from the face-name-keyed :attr:`SNMesh.bc` dict
(``sn_mesh.bc["xmin"]``, ``sn_mesh.bc["xmax"]``, ...; see
:ref:`bc-face-name-carve`).

For 1D meshes (``ny=1``):

- **Gauss--Legendre** quadrature takes the fast 1-D chain-scan
  (:class:`~orpheus.sn.loss_representation.CumprodScan`) path (all
  :math:`\mu_y = 0`, so no y-streaming).
- **Lebedev** quadrature falls through to the 2D wavefront sweep.
  Ordinates with :math:`\mu_x \neq 0` stream along *x*; the
  *y*-streaming terms cancel via reflective BCs on the single-cell
  *y*-dimension.  Ordinates with :math:`\mu_x = \mu_y = 0`
  (z-directed) reduce to pure collision:
  :math:`\psi = Q \cdot w_{\text{norm}} / \Sigt{}`.

Both quadratures recover the analytical eigenvalue exactly on
homogeneous problems (verified to machine precision for 1G/2G/4G).


The Transport Equation
======================

The 1-D slab form — the base of the broadening progression — is
posed in :doc:`slab_one_group`; the multi-group energy extension in
:doc:`slab_multigroup`; the 2-D Cartesian form in
:doc:`cartesian_multid`. The curvilinear geometries below extend
them.

Spherical 1D
-------------

In spherical coordinates the transport equation acquires an **angular
redistribution term** that couples ordinates:

.. math::
   :label: transport-spherical

   \mu \frac{\partial \psi}{\partial r}
   + \frac{1 - \mu^2}{r} \frac{\partial \psi}{\partial \mu}
   + \Sigt{} \psi = \frac{Q}{W}

The curvature term :math:`(1 - \mu^2)/r \cdot \partial\psi/\partial\mu`
arises because a neutron streaming radially at angle :math:`\mu` *rotates*
its direction cosine as it moves to a different radius.  Discretising this
term requires diamond difference in **both space and angle**.

Cylindrical 1D
---------------

For an infinitely long cylinder with azimuthal symmetry, the transport
equation in the radial variable :math:`r` is:

.. math::
   :label: transport-cylindrical

   \frac{\eta}{r} \frac{\partial(r\psi)}{\partial r}
   - \frac{1}{r} \frac{\partial(\xi\psi)}{\partial\varphi}
   + \Sigt{} \psi = \frac{Q}{W}

where the direction cosines are:

- :math:`\eta = \sin\theta\cos\varphi` --- radial projection (streaming)
- :math:`\xi = \sin\theta\sin\varphi` --- azimuthal component
- :math:`\mu = \cos\theta` --- axial component

The constraint :math:`\eta^2 + \xi^2 + \mu^2 = 1` holds.  The azimuthal
redistribution :math:`-\partial(\xi\psi)/\partial\varphi` couples ordinates
on each :math:`\mu`-level.

The Discrete Balance Equation
=============================

This is the core of the S\ :sub:`N` method.  The balance equations are
presented from simplest to most complex, in the progression chapters:
Cartesian geometries have no angular redistribution — the 1-D slab
balance is derived in :doc:`slab_one_group`, the 2-D Cartesian balance
in :doc:`cartesian_multid`; curvilinear geometries add :math:`\alpha`
coupling and a geometry factor :math:`\Delta A/w` —
:doc:`curvilinear_one_group`.

.. _cell-update-strategies:

Cell update strategies (the strategy contract)
==============================================

The discrete balance equation (the slab DD :eq:`dd-cartesian-1d` of
:doc:`slab_one_group`, the M-M-closed curvilinear update :eq:`dd-solve`
of :doc:`curvilinear_one_group`) yields, for
each cell, a small algebraic system: combine the upstream face flux
(and, for sphere/cylinder, the upstream angular half-flux) with a
local source and the cell's total cross section to produce the
cell-average flux plus the downstream states.  The closure relating
:math:`\overline{\psi}_i` to :math:`\psi_{i-1/2}` and
:math:`\psi_{i+1/2}` is **not unique** — Diamond Difference (DD),
weighted DD, Linear Discontinuous (LD), Step, and Exponential
Characteristic (EC) are all valid choices, each with different
truncation error, positivity, and cost.  Per Cardinal Rule 2
(architecture), the cell-update math is **the same algebra** in slab,
sphere, and cylindrical 1-D — only the populated fields of the
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` packet
change.  Lifting the closure into a strategy contract makes the
sweep driver thin and lets each closure be unit-tested in isolation.

The strategy contract is owned by
:mod:`orpheus.transport.spatial.scheme`.

The Protocol
------------

The contract is a ``@runtime_checkable`` ``typing.Protocol`` —
satisfied by structural typing, not inheritance — exposing two class-
level traits and a single :meth:`update` method:

* :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme`

  - ``is_linear: bool`` — whether the closure is linear in its inputs.
    Diamond Difference is linear; Step's positivity-fixup, weighted-DD
    with a flux-dependent weight, and EC with a clipped argument are
    not.
  - ``is_positivity_preserving: bool`` — whether non-negative inputs
    guarantee non-negative outputs.  Diamond Difference is **not**
    positivity preserving (Lewis & Miller §5.3, where DD's tendency
    to produce negative cell-edge fluxes is exhibited and motivates
    the choice of Step or weighted-DD in stiff cells); Step is.
  - ``update(visit, total_xs, source, upstream_state) ->
    CellResult`` — the cell update itself.  ``visit`` is a
    :class:`~orpheus.transport.spatial.scheme.CellVisit` packet (see
    next subsection) that combines the geometric
    :class:`~orpheus.geometry.reduced_operator.StreamingTerms` with
    sweep-direction-resolved data.

The two helper dataclasses (frozen, slotted) carry the per-cell
state:

* :class:`~orpheus.transport.spatial.scheme.UpstreamState`

  - ``spatial_upstream: np.ndarray`` — shape ``(ng,)``.  Face flux
    entering the cell from the upstream face (always populated).
  - ``angular_upstream: np.ndarray | None`` — shape ``(ng,)`` for
    sphere/cylinder; ``None`` for slab.  :math:`\psi_{n-1/2,\,i}`,
    the half-flux at the upstream half-angle.

* :class:`~orpheus.transport.spatial.scheme.CellResult`

  - ``cell_average_flux: np.ndarray`` — shape ``(ng,)``.  The cell-
    average flux :math:`\overline{\psi}_i = \mathrm{numer}/\mathrm{denom}`
    from the closure's algebraic solve.
  - ``outgoing_spatial_flux: np.ndarray | None`` — shape ``(ng,)`` in
    the typical case; ``None`` for the cylindrical pure-azimuthal
    degenerate case where the cell has no radial face flow (see
    below).
  - ``outgoing_angular_state: np.ndarray | None`` — shape ``(ng,)``
    for sphere/cylinder; ``None`` for slab.  :math:`\psi_{n+1/2,\,i}`
    from the Morel--Montry closure.

The SN sweep DAG and ``CellVisit``
-----------------------------------

The SN sweep is a **topological sort of a directed cell graph**.
For a given ordinate :math:`\Omega_n`, every face :math:`f` of the
mesh is oriented by the sign of :math:`\Omega_n \cdot \hat n_f` — an
edge from cell :math:`A` to cell :math:`B` if :math:`\Omega_n` points
from :math:`A` into :math:`B` across that face.  The sweep walks
cells in a topological order over this DAG so that, when each cell is
visited, all its upstream face fluxes are already known.  This is
the SN-specific graph-theoretic concept; MoC uses a different
mathematical structure (fiber bundles + solution sheaves over
characteristic curves), and CP / diffusion / MC have no sweep at
all.  Per Cardinal Rule 2 (architecture), no shared
``SweepGraph`` Protocol is hoisted across solvers — the sweep DAG
lives in :mod:`orpheus.sn`.

The contract's :meth:`update` consumes a
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet rather
than a raw
:class:`~orpheus.geometry.reduced_operator.StreamingTerms`.
The :class:`CellVisit` composes:

* ``cell_idx: int`` — the cell being visited.
* ``streaming_terms: StreamingTerms`` — the **purely geometric**
  primitive (``face_area_inner`` / ``face_area_outer`` are
  geometric labels — inner = closer to :math:`r=0`, outer =
  farther — independent of sweep direction).
* ``face_area_downstream: float | None`` — **sweep-direction-
  resolved**.  For an outward sphere or cylinder sweep
  (:math:`\mu \ge 0`) it equals ``streaming_terms.face_area_outer``;
  for an inward sweep (:math:`\mu < 0`) it equals
  ``streaming_terms.face_area_inner``.  ``None`` for slab (slab DD
  does not read face areas) and for the cylindrical pure-azimuthal
  degenerate case (no spatial flow).

The :class:`CellVisit` packets are produced by
:meth:`SNMesh.dag_walk(*, ordinate_idx=..., direction_sign=...,
mu_level_idx=None)
<orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk>` — a generator that
yields cells in DAG-topological order.  The method takes EXACTLY ONE
of ``ordinate_idx`` (single-ordinate visits, used by the sweep
driver) or ``direction_sign`` (direction-keyed visits, used by the
apply matvec) as a keyword-only argument.  Both invocation modes
encapsulate the inward / outward branching, the cylindrical
per-:math:`\mu`-level traversal, and the pure-azimuthal degenerate
handling.  The sweep at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) consumes this
generator::

    for visit in sn_mesh.dag_walk(ordinate_idx=n):
        upstream = UpstreamState(
            spatial_upstream=psi_face,
            angular_upstream=psi_angle[visit.cell_idx],
        )
        result = scheme.update(
            visit, total_xs, source, upstream,
        )
        ...

The cell-update strategy receives only **resolved** data — no
sign-of-:math:`\mu` branching inside the strategy.  This pattern
moves the graph-theoretic concept to where it belongs (the SN
module) and keeps the geometry-layer
:class:`~orpheus.geometry.reduced_operator.StreamingTerms`
geometry-only and reusable by future MoC / CP / diffusion modules
that have different mathematical structures.

Slab vs curvilinear discrimination
-----------------------------------

.. note:: **Superseded mechanism (Issue #196 Phase G Step 2.5 → Issue
   #236 Step C).** The ``alpha_in is None`` slab/curvilinear
   discrimination described below was **retired**: Issue #196 Phase G
   Step 2.5 gave slab *neutral* curvature (``face_area_inner =
   face_area_outer = 1.0``, ``delta_A_over_w = 0.0``) so the unified
   cell-balance helper consumes the same packet regardless of geometry,
   and Issue #236 Step C removed the Morel--Montry ``alpha_in`` /
   ``alpha_out`` / ``tau_mm`` fields from
   :class:`~orpheus.geometry.reduced_operator.StreamingTerms` entirely
   (it is now **purely geometric**; :math:`\tau` is closure-owned — see
   :ref:`sn-tau-c-on-cellvisit-live`).  Slab is now distinguished at the
   sweep level by ``upstream_state.angular_upstream is None``, not by a
   ``StreamingTerms`` field test.  The prose below records the historical
   pre-Step-2.5 convention and is pending a dedicated rewrite to the
   neutral-curvature mechanism; the authoritative current description is
   the :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
   docstring.

A strategy distinguishes slab from curvilinear by a single field test
on the visit's streaming terms:

* **Slab** — ``visit.streaming_terms.alpha_in is None`` (and the rest
  of the curvature bundle, ``alpha_out``, ``delta_A_over_w``,
  ``tau_mm``, ``face_area_inner`` / ``face_area_outer``, are all
  ``None``).  ``upstream_state.angular_upstream is None``.  The
  strategy returns ``CellResult(outgoing_angular_state=None, ...)``.
* **Sphere or cylinder** —
  ``visit.streaming_terms.alpha_in is not None``; the full curvature
  bundle is populated.  ``upstream_state.angular_upstream`` carries
  :math:`\psi_{n-1/2,\,i}`.  The strategy returns the M-M-closed
  ``outgoing_angular_state``.

This single-field discrimination convention is locked in by
foundation-tier protocol-conformance tests in
``tests/sn/sweep/core/test_discretization_scheme_protocol.py``; concrete
strategies and the Wave D sweep rewrite read this same field to
dispatch.

Cylindrical pure-azimuthal degenerate case
-------------------------------------------

For cylindrical 1-D radial sweeps with a product or level-symmetric
quadrature, ordinates with axial direction cosine
:math:`|\mu_z| \to 1` have radial direction cosine
:math:`|\eta| = \sqrt{1 - \mu_z^2} \to 0`.  In this limit the cell
has **no radial face flow** — the streaming term
:math:`\mu_x \cdot \partial_r` vanishes — and the cell-update
algebra collapses to the redistribution-only form

.. math::

   \mathrm{denom} = (\Delta A / w)\,c_{\rm out} + \Sigma_t\,V_i,
   \qquad
   \mathrm{numer} = Q_i\,V_i + (\Delta A/w)\,c_{\rm in}\,\psi_{n-1/2,\,i},

with no spatial-flux contribution.  The strategy contract signals
this case by setting ``CellResult.outgoing_spatial_flux = None``: the
sweep driver, on receiving ``None``, skips the face-flux update for
that cell.  The angular M-M closure remains active — angular
redistribution physics is still present.

The numerical threshold is ``streaming_terms.abs_mu < 1e-15``, with
``abs_mu`` populated from the **global ordinate**
:math:`|\eta|` on the streaming-terms packet (resolved through
``level_indices`` for cylindrical geometry — see
:doc:`/theory/foundations/structured_geometry`, "Connection coefficients (reduced
streaming operator)").  In this case
:meth:`SNMesh.dag_walk
<orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk>` yields visits with
``face_area_downstream = 0.0`` to signal "no spatial flow" to the
strategy (Issue #196 Step 2.5 retired the ``None`` sentinel — the
slab carries ``1.0`` and degenerate cylindrical carries ``0.0`` so
the cell-balance helper consumes one geometry-blind number).

The DD recurrence
------------------

For non-degenerate cells, the closure relation reduces — for slab
geometry, the cell-update math is the DD recurrence
:eq:`dd-recurrence` (derived at :ref:`sweep-cumprod` in
:doc:`slab_one_group`); for curvilinear
geometry, it is the curvilinear DD form combining the
:math:`\Delta A/w` redistribution with the M-M angular closure.  The
sweep driver inlines this math today; Wave D (Issue #159) will
rewrite the driver to dispatch through a
:class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` strategy.

The first concrete strategy — :class:`DiamondDifference` — is shipped
in Round 2 of Wave C (Issue #158) as a bit-identical extraction of
the existing inlined sweep math.  Linear Discontinuous (Lewis &
Miller §5.3 — preview), Exponential Characteristic, and Step
strategies are deferred to a Wave C-extension session, each with its
own MMS spatial-convergence verification.

Diamond Difference
------------------

The first concrete strategy is
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`
(:mod:`orpheus.transport.spatial.diamond`).  It implements the **same**
algebra as the existing inlined sweep — Round 2 of Wave C is a
bit-identical extraction, gated by ``np.array_equal`` hand-calc tests
in ``tests/sn/sweep/core/test_diamond.py`` against the sweep's scalar
formulas at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``).

Per Wave C decision **D5** (one geometry-polymorphic class), the
strategy is a single :class:`DiamondDifference` that handles slab,
sphere, and cylinder by branching on two
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` fields:
``alpha_in is None`` (slab vs curvilinear) and ``abs_mu < 1e-15``
(cylindrical pure-azimuthal degenerate vs not).

**Slab branch** (``streaming_terms.alpha_in is None``).  The flat /
Cartesian DD recurrence reduces to the per-cell scalar form of
:eq:`dd-recurrence`:

.. math::
   :label: dd-slab-scalar

   \overline{\psi}_i \;=\; \tfrac12\bigl(\psi_{i-1/2}
                                         + \psi_{i+1/2}\bigr),
   \qquad
   \psi_{i+1/2} \;=\; \frac{2|\mu_n| - \Delta x_i\,\Sigma_t}
                              {2|\mu_n| + \Delta x_i\,\Sigma_t}\,
                       \psi_{i-1/2}
                   \;+\; \frac{2\,Q_i\,\Delta x_i / W}
                              {2|\mu_n| + \Delta x_i\,\Sigma_t},

with ``W = Σ_n w_n`` the quadrature weight sum, mirroring the
sweep's vectorised cumprod path
(``_sweep_1d_cumprod`` (the dissolved ``sweep.py``) lines 117–123) and per-
cell solver (``_solve_recurrence`` (the dissolved ``sweep.py``) lines 208–
222) at the operation level.  Per the strategy contract, ``source``
arrives at the cell update **already weight-normalised** by the
sweep — for slab, ``source = Q · Δx / W`` (and slab cell volume is
``V = Δx``).  For slab, the strategy sets
:attr:`~orpheus.transport.spatial.scheme.CellResult.outgoing_angular_state`
to ``None`` — slab geometry has no angular redistribution.

**Curvilinear branch** (``streaming_terms.alpha_in is not None`` and
``abs_mu ≥ 1e-15``).  Sphere or cylinder, away from the cylindrical
pure-azimuthal degenerate case.  The strategy couples the M-M
angular closure :eq:`mm-weights` to the WDD spatial closure
:eq:`wdd-closure`, with the redistribution constants

.. math::
   :label: dd-mm-closure-constants

   c_{\rm out} \;=\; \alpha_{n+\tfrac12}/\tau_n,
   \qquad
   c_{\rm in}  \;=\; \frac{1 - \tau_n}{\tau_n}\,\alpha_{n+\tfrac12}
                       + \alpha_{n-\tfrac12},

built from the Bailey 2009 dome :eq:`alpha-recursion` and the
Morel–Montry weight :eq:`mm-weights`.  The cell-update is then

.. math::
   :label: dd-curvilinear-scalar

   \overline{\psi}_{n,i} \;=\;
       \frac{Q_i\,V_i / W
             + |\mu_n|\,(A_{i-1/2} + A_{i+1/2})\,\psi^s_{n,\,{\rm in}}
             + (\Delta A_i / w_n)\,c_{\rm in}\,\psi_{n-\tfrac12,\,i}}
            {2|\mu_n|\,A^s_{\rm out}
             + (\Delta A_i / w_n)\,c_{\rm out}
             + \Sigma_t\,V_i},

mirroring ``_sweep_1d_spherical`` (the dissolved ``sweep.py``) lines
350–355 (and the structurally identical cylindrical branches at
sweep.py:511–531 / sweep.py:548–575) verbatim, with closures

.. math::

   \psi^s_{\rm out} \;=\; 2\overline{\psi}_{n,i}
                        - \psi^s_{n,\,{\rm in}},
   \qquad
   \psi_{n+\tfrac12,\,i} \;=\;
       (\overline{\psi}_{n,i}
         - (1 - \tau_n)\,\psi_{n-\tfrac12,\,i})/\tau_n.

**Cylindrical pure-azimuthal degenerate branch**
(``streaming_terms.alpha_in is not None`` and ``abs_mu < 1e-15``).
For a level whose axial direction cosine :math:`|\mu_z| \to 1`, the
radial direction cosine :math:`|\eta| \to 0` and the cell has no
radial face flow — the :math:`2|\mu| A_{\rm out}` and
:math:`|\mu|(A_{\rm in} + A_{\rm out})\,\psi^s_{\rm in}`
contributions drop out:

.. math::

   \mathrm{denom} \;=\; (\Delta A / w)\,c_{\rm out} + \Sigma_t\,V_i,
   \qquad
   \mathrm{numer} \;=\; Q_i\,V_i / W
                       + (\Delta A / w)\,c_{\rm in}\,
                          \psi_{n-\tfrac12,\,i},

mirroring ``_sweep_1d_cylindrical`` (the dissolved ``sweep.py``) lines
533–543 verbatim.  The strategy returns
:attr:`~orpheus.transport.spatial.scheme.CellResult.outgoing_spatial_flux`
``= None`` to signal "no face-flux write" to the sweep driver; the
M-M angular closure remains active.

**Traits and forward references.**  Diamond Difference has

* :attr:`~orpheus.transport.spatial.diamond.DiamondDifference.is_linear`
  ``= True`` — the cell average and downstream states are affine
  combinations of ``source`` and ``upstream_state``;
* :attr:`~orpheus.transport.spatial.diamond.DiamondDifference.is_positivity_preserving`
  ``= False`` — Lewis & Miller §5.3 exhibits the canonical thin-
  cell / large-source counter-example where DD's
  :math:`\psi_{\rm out} = 2\overline{\psi} - \psi_{\rm in}` produces
  negative outgoing flux from positive inputs.

Wave C-extension and Wave D will ship
:class:`Step` (positivity-preserving, :math:`\mathcal{O}(\Delta x)`),
:class:`LinearDiscontinuous` (:math:`\mathcal{O}(\Delta x^2)`,
better robustness in optically-thick cells), and
:class:`ExponentialCharacteristic` (positivity-preserving by
construction) as alternatives, each with its own MMS spatial-
convergence verification.

References
----------

* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §5.3 covers Diamond Difference, weighted-DD,
  Step, and Linear Discontinuous closures and their positivity /
  truncation properties; §4.5 covers the Morel--Montry angular
  closure used for :math:`\psi_{n+1/2,\,i}`.
* Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
  *A piecewise linear finite element discretization of the
  diffusion equation for arbitrary polyhedral grids*.
  JCP 227, 3738--3757.  Eq. 50 (dome recursion), Eq. 74
  (Morel--Montry).

See also
--------

* :mod:`orpheus.transport.spatial.scheme` — the contract module.
* :mod:`orpheus.transport.spatial.diamond` — the
  :class:`~orpheus.transport.spatial.diamond.DiamondDifference` concrete
  strategy.
* :doc:`/theory/foundations/structured_geometry`, "Connection coefficients (reduced
  streaming operator)" — the upstream side of the contract: where the
  per-cell, per-direction streaming-terms packet is built.


.. _sweep-algorithm:

Sweep Algorithm
===============

Because each cell's outgoing flux becomes the next cell's incoming flux,
the equations must be solved in the direction of neutron travel --- this
is called a **transport sweep**.

The 1-D slab sweep — the cumprod recurrence and the generic affine
outflow reconstruction — is derived in :doc:`slab_one_group`; the 2-D
Cartesian wavefront, its octant dependency graph, and the
multi-dimensional LD (UBLD) system in :doc:`cartesian_multid`; the 2-D
LD stress MMS in :doc:`verification`; the curvilinear machinery — the
sequential ordinate sweep, the pole angular closure (#168 Phase B),
the sweep-frame apply matvec (Phase C), and the direct :math:`\psi_{1/2}`
starting-direction solve (#282 route (a)) — in
:doc:`curvilinear_one_group`. The sections below preserve the #168 /
#195 / #282 campaign record (the dispatch consolidation, the Carlson
phases D–F, ERR-058).

.. _unified-sweep-dispatch:

Unified sweep dispatch
-----------------------

.. note::

   **Superseded (coupled-block campaign step 6, R-6.1, 2026-07-12).**
   This section records the Wave-D Round-2 consolidation (#161) — the
   *first* unification of the four sweep paths under a single entry, the
   operator-free ``transport_sweep``.  That entry was itself retired at
   step 6; the dispatch is now the first-class ``LossRepresentation``
   polymorphism (:func:`~orpheus.sn.loss_representation.default_for` selects
   :class:`~orpheus.sn.loss_representation.CumprodScan` for the 1-D scan,
   the anti-hyperplane wavefront for multi-D), documented in
   :doc:`/theory/methods/sn/loss_representation`.  The Wave-D narrative below is preserved as
   the origin of that unification: read ``transport_sweep`` and the
   ``ReducedStreamingOperator`` boolean-dispatch code as the then-current
   implementation, not today's.

Wave D Round 2 of the SN reshape campaign (Issue #161) consolidated
the four pre-existing sweep paths (1-D cumprod / 2-D wavefront /
spherical / cylindrical) under one operator-free ``transport_sweep``
entry point that branched
on a single boolean from the
:class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
primitive (Wave B Issue #6 / Wave D Round 1):

.. code-block:: python

   def transport_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso=None):
       reduced = sn_mesh.reduced
       if reduced is not None and reduced.requires_upstream_angular_state:
           return _curvilinear_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
       return _cartesian_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)

The pre-Wave-D dispatch did string-equality on
``sn_mesh.curvature == "spherical"`` / ``"cylindrical"`` /
``None``.  The new dispatch reads the canonical geometry-layer
property
:attr:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.requires_upstream_angular_state`
— ``False`` for slab + 2-D Cartesian (no angular redistribution
between successive half-angles), ``True`` for spherical +
cylindrical.  Two-D Cartesian sets ``sn_mesh.reduced is None``
(no curvilinear math is needed), and the dispatch falls through
to the Cartesian path.  Why this matters:

* The :class:`ReducedStreamingOperator` is the primitive that
  already encodes "does this geometry need angular
  redistribution to march the sweep?", so the dispatch reads its
  property directly instead of round-tripping through a
  string tag — Cardinal Rule 2 (architecture).  Consumers
  outside the SN sweep (MoC, CP) read the same property when
  they need the same dispatch.
* The dispatch surface shrinks from four string-equality checks
  to one boolean — a structural simplification that makes the
  control flow easier to reason about and to extend with
  additional cell-update strategies (Wave C-extension).

Cell update strategy parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The curvilinear sweep dispatches per-cell to
:meth:`~orpheus.transport.spatial.scheme.DiscretizationScheme.update`:

.. code-block:: python

   for cell_idx in march_order:
       st = reduced.streaming_terms(cell_idx, dir_idx, mu_level_idx=p)
       upstream = UpstreamState(
           spatial_upstream=psi_face,
           angular_upstream=psi_angle[cell_idx],
       )
       result = sn_mesh.scheme.update(
           streaming_terms=st,
           total_xs=sig_t[cell_idx],
           source=QV[cell_idx],
           upstream_state=upstream,
       )
       psi_face = result.outgoing_spatial_flux  # may be None for cylindrical degenerate
       psi_angle[cell_idx] = result.outgoing_angular_state

The cell-update strategy lives on
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.scheme` (introduced in
this round as a constructor argument with default
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`).  The
default reproduces the inlined sweep math bit-identically — every
regression snapshot at ``tests/sn/regression/snapshots/`` was
generated with DD and continues to match bit-for-bit when the
unified sweep dispatches via ``scheme.update(...)``.  See
:ref:`cell-update-strategies` for the strategy contract and
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` for the
DD scalar form.

Wave C-extension will ship :class:`Step`, :class:`LinearDiscontinuous`,
and :class:`ExponentialCharacteristic` as positivity-preserving /
higher-order alternatives; the unified dispatch infrastructure is
in place to receive them — users will pass
``scheme=LinearDiscontinuous()`` etc. at
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` construction.

The 1-D cumprod fast path (DD-only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Cartesian dispatch checks three preconditions before
selecting the 1-D cumprod fast path:

#. ``scheme is DiamondDifference`` — the cumprod recurrence
   :eq:`dd-recurrence` is a DD-specific algebraic identity
   (Lewis & Miller §5.3); LD / EC / Step closures do not admit
   the same recurrence.
#. Quadrature is GL1D (``ny == 1`` and all ``mu_y`` vanish).
#. Source is isotropic (``Q_aniso is None``).

If any precondition fails, the Cartesian dispatch routes
through the 2-D wavefront sweep (which handles 1-D as a
special case).  Preserving the cumprod fast path inside the
unified algorithm is required to keep the 1-D regression
snapshots bit-identical and to retain the historical
sub-millisecond sweep time for typical 1-D problems.

The 2-D wavefront sweep dispatches its DD per-cell algebra
through the storage-free kernel pair
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
/ :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
on the strategy attached to the
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (Wave 2 of the SN
performance plan; closes Issue #4 — see
:ref:`sweep-octant-dependency-graph` for the full architecture and
:ref:`sweep-dispatch-relayering` for the S6.4(e) re-layering).
The "inlined DD math" formerly carried inside
:func:`~orpheus.sn.loss_representation._sweep_jacobi` was lifted into
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` as a single
bit-identical extraction, vectorised across the
``(N_oct, n_diag, ng)`` slice — the ordinate axis, anti-diagonal
axis, and group axis simultaneously.  Wave C-extension's LD / EC
/ Step closures override the kernel pair and become drop-in
alternatives at SNMesh construction time:
``SNMesh(mesh, quad, scheme=LinearDiscontinuous())``.  The
open design point of "how to parameterise the 2-D wavefront
without breaking anti-diagonal vectorisation" is now closed: the
storage walk is the scheduler, the level operation owns the
direction fork, and the kernel pair is the closure — the contract
is per-level batched evaluation.

ERR-026 closure status (partial through Wave E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Superseded (2026-06-12, Issue #195).**

   This subsection records the Wave-E-era reading: ERR-026 PARTIAL,
   the curvilinear ``"krylov"`` default "would regress MMS to
   :math:`\mathcal{O}(h)`", the open second-order follow-up.  That
   reading was the best available *then* and is preserved as history.
   It is now **superseded**: ERR-058 (#195) showed the wrong fixed
   point was the *closure-seed* family, not a boundary-truncation
   order; the curvilinear default returned to ``"source_iteration"``
   (SI :math:`\equiv` Krylov bit-identical post-unification); and the
   isotropic MMS is :math:`\mathcal{O}(h^2)`-consistent.  See
   :ref:`sn-err-058-closure-seed-closeout` for the mechanism,
   structural obstruction, and production decision.  The numerical
   values below stay as bug-era evidence; their *interpretation* is
   carried by the close-out.

The curvilinear sweep's one-directional WDD closure
:math:`\psi_{n+1/2} = (\overline{\psi} - (1 - \tau_{mm})\,
\psi_{n-1/2})/\tau_{mm}` is preserved bit-identically by
:class:`DiamondDifference` (Wave C extracted it from the
inlined sweep verbatim).  ERR-026 (catalogued in
:doc:`/development` and the V&V matrix at
:doc:`/verification/matrix`) lives in this closure: the
solver's source-iteration path converges to a non-flat
fixed point even though the matrix-free ``apply`` path with
the symmetric closure is exact for **constant** sources.

Wave D's gating contract was bit-identity for ``scheme =
DiamondDifference`` — the bug is preserved by construction so
the regression snapshots stay green.  Wave E (Issues #98 #99 #164)
took two passes at the closure:

* **Wave E Round 2** wired
  :func:`~orpheus.sn.solver.solve_sn_fixed_source` to route
  through Krylov-on-:meth:`InvertibleOperator.apply` (the
  symmetric closure) with the sweep-as-``solve`` as preconditioner.
  This closes ERR-026 on **constant-source reflective-BC
  problems** — the canonical
  :file:`tests/sn/test_sweep_operator_inconsistency.py` regression
  suite confirms the krylov path gives the analytical flat flux
  to round-off where the sweep does not.
* **Wave E Round 3** (Issue #98 follow-up) closed the BC-faithfulness
  gap that Round 2 identified: the FD operator's
  then-``solution_to_angular_flux*`` codec and the
  matvec helpers consumed the
  :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances on
  the :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (Wave B Issue 7
  tensor-decomposed BC algebra), dispatching boundary fills via the
  realiser-routed 1-arg :meth:`apply` on the resolved
  :class:`~orpheus.numerics.operator.LinearOperator`. Vacuum,
  reflective, white, periodic, albedo, and mixed BCs are now
  plumbed uniformly through the FD operator; bit-identity to the
  pre-Round 3 hard-coded reflective fill is preserved for
  :class:`ReflectiveBoundary(axis=…, albedo=1.0)` (the standard
  ``BC.reflective`` case), which is the load-bearing condition for
  the 11 frozen regression snapshots to stay green. (Note:
  post Issue #186 / B3 + β2 the law itself is a pure descriptor;
  ``BoundaryTraceLaw.apply`` no longer exists. The realiser
  produces the 1-arg :class:`LinearOperator` whose :meth:`apply`
  the matvec calls; the Wave-E Round-3 prose describes the
  contract as it existed at the time, but the architectural
  conclusion — uniform BC consumption through a 1-arg
  ``apply`` — is the same.)

What is **still** open after Round 3: empirically the symmetric-
closure FD operator at the curvilinear outer face uses cell-center
as a face-flux approximation (``psi_right = fi[:, n, i, 0]`` at
``i = nx-1`` for outgoing :math:`\mu > 0`).  This is exact for
constant solutions but only first-order accurate on non-constant
solutions like the manufactured ``A(r) = sin(πr/R)`` ansatz used
by the curvilinear MMS test suite.  Switching the
``solve_sn_fixed_source`` curvilinear default from
``"source_iteration"`` to ``"krylov"`` would *regress* the MMS
convergence rate from the WDD sweep's
~:math:`\mathcal{O}(h^{1.3})` (ERR-026-affected, but a benign
volumetric-error mode for these MMS) to
~:math:`\mathcal{O}(h^{1})` (FD operator's boundary truncation).
Round 3 therefore *keeps* ``inner_solver="source_iteration"`` as
the default for all geometries; ``"krylov"`` is opt-in and
correct for constant-source problems but not the right default
for MMS.

The two ``xfail-strict`` tripwires at
``tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py``
remain ``xfail`` through Round 3 with updated reason strings
reflecting the partial closure.  Full ERR-026 closure on MMS
depends on a follow-up that extrapolates the curvilinear
outer-face flux at second order (DD diamond relation at the
boundary, or analogous ghost-cell technique).

Adams & Larsen 2002 §III.B's "preconditioner correctness vs
operator correctness" frame is the right lens: the sweep's WDD
fixed-point bias is the wrong answer for a *primary solve*, but
as a *preconditioner* the same fixed point is just an effective
scaling of the residual — it does not poison the converged
solution determined by the operator.  The operator must be
correct *and* second-order-accurate; Round 3 closed the
correctness piece (BC-faithfulness), the second-order piece is
the open follow-up.

.. _sn-boundary-face-flux-protocol:

Boundary face-flux strategies — Phase A (RETIRED in Phase C)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Issue #168 empirical investigation (recorded at
``.claude/agent-memory/numerics-investigator/issue_168_three_defects.md``)
found **three independent O(1) boundary truncation defects** in the
historical curvilinear FD operator.  Phase A (2026-05-10) addressed
**Defects 1 + 2** via a ``BoundaryFaceFlux`` strategy Protocol — a
one-sided second-order DD diamond extrapolation
:math:`\psi^{\text{face}}_{N-1/2} = \tfrac{3}{2}\,\psi_{N-1} -
\tfrac{1}{2}\,\psi_{N-2}` plus a structural decoupling of
cell-centre storage from BC face-value storage in
the then-``solution_to_angular_flux_spherical`` codec
(returning a ``(fi, boundary_face_flux)`` tuple where ``fi`` was
pure cell-centre storage and the BC face flux lived in its own
companion array).

**Phase C retired the Phase A Protocol entirely.** The sweep-frame
apply matvec subsumes the boundary-face closure into the WDD
propagation chain — the BC trace law owns the boundary edge per
the §16A.3 contract, no separate algebraic extrapolation is needed.
The retired symbols are:

* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`` (Protocol)
* ``orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`` (default)
* ``orpheus.sn.spatial.boundary_face_flux.CellCenter`` (ablation)
* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFluxBase`` (ABC)
* The ``boundary_face_flux`` field on
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`
* The 21 foundation tests at
  :file:`tests/sn/sweep/test_boundary_face_flux.py`

See :ref:`phase-c-sweep-frame-matvec` for the replacement
architecture. The Phase A subsection is preserved as historical
context for the empirical-defects-investigation reasoning chain.

.. _sn-phase-d-carlson-coupled-pole-sweep:

Phase D Carlson coupled-pole sweep (Issue #168 Phase D)
-------------------------------------------------------

.. attention:: **Superseded by Issue #282 route (a) (2026-07-04).**

   The swappable ``PsiHalfAngleSeed`` strategy family
   (``ZeroSeed`` / ``CarlsonInwardSweep`` / ``AngularEdgeExtrapolation``)
   whose design this section — and the Phase F and ERR-058 sections that
   follow — build up was **retired** by Issue #282 route (a).  The
   starting-direction half-angle flux :math:`\psi_{1/2}` is now
   **first-class typed state** the sweep marches *directly* from the true
   within-group source, not a functional of the previous iterate.  Any
   "current default / retained strategy / the seed lives as a strategy
   field" claim in the three sections below is **historical** — read them
   for the *why* (what was tried and the diagnoses that narrowed the
   defect), but for the CURRENT design see
   :ref:`sn-282-direct-starting-direction-solve`.  In particular the
   ``AngularEdgeExtrapolation`` "iterate extrapolation" seed those
   sections land on was itself the #282 walk-order back edge that route
   (a) removes.

.. admonition:: Key Facts
   :class: important

   * Phase D (commit landed 2026-05-12 on
     ``refactor/sn-operator-algebra``) closes the structural bug
     in :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
     by replacing the hardcoded ``psi_half_left = 0`` seed with
     the canonical Hébert §3.9.4 Eqs. (3.432)–(3.435) inward
     :math:`\mu = -1` sweep output.
   * The seed lives in the **M-M angular recurrence**
     (:func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
     in ``orpheus/sn/sweep/pole_angular_closure.py``), **NOT**
     in the WDD spatial pole-face initial condition the
     :ref:`Phase C plan <sn-curvilinear-trajectory-resolvent-crosscheck-section>`
     proposed.  The diagnostic memo at
     ``.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md``
     empirically falsified intervention ``[A]`` (WDD pole-face
     replacement) and confirmed intervention ``[B]`` (M-M
     half-angle seed replacement).
   * Architectural choice is **Option α (composition)**: the
     seed lives as a
     ``PsiHalfAngleSeed``
     strategy field on :class:`MorelMontryAngularSweep`, not as a
     sibling Protocol on :class:`SNMesh`. The Legacy / Bailey
     closures have no ``psi_half_left`` variable to seed; a
     sibling Protocol would force every consumer to handle an
     irrelevant Protocol.
   * The **L = 0 isotropic-only** assumption is load-bearing: the
     apply matvec's :math:`L` operator currently carries only
     :math:`\Sigma_t \psi` (scattering is composed externally
     via a separate operator).  A future refactor that moves
     scattering INTO :math:`L` MUST extend the moment-folded
     source in :eq:`hebert-3-432-source` to include
     :math:`\ell \ge 1` terms.

The Hébert §3.9.4 equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hébert §3.9.4 (pp. 141–144 of [Hebert2009]_) opens the sphere
difference relations at Eq. (3.418) (angularly-integrated
divergence form), introduces the :math:`\alpha`-recursion
:eq:`bailey-dome-recursion` and the cell-balance with
redistribution divisor :math:`\Delta S_i / (2\,\mathcal{W}_n)` at
Eq. (3.428), and then specialises to the auxiliary starting
direction :math:`\mu = -1`.  At this direction the
angular-redistribution coefficient :math:`(1 - \mu^2)` is
**identically zero**, so the streaming–collision balance
decouples from the :math:`\alpha`-cascade and reduces to a plain
DD inward recurrence in radius.

The continuous form at :math:`\mu = -1` (Hébert Eq. (3.432)) is

.. math::
   :label: hebert-3-432

   -\frac{\partial}{\partial r}\,\phi_{-1/2}(r)
   \;+\; \Sigma(r)\,\phi_{-1/2}(r)
   \;=\; \sum_{\ell=0}^{L}
         \frac{2\ell + 1}{2}\,Q_\ell(r)\,P_\ell(-1).

The subscript :math:`-1/2` is Hébert's half-integer index for the
auxiliary starting ordinate — it labels the **inward
zero-weight** direction that sits one half-step above
:math:`\mu = -1` in the :math:`\alpha`-cascade, not a physical
ordinate at :math:`\mu = -0.5`.  The right-hand side is the
Legendre expansion of the scattering source :math:`Q` evaluated
at :math:`\mu = -1`, where :math:`P_\ell(-1) = (-1)^\ell`.

For an **isotropic** operator (:math:`L = 0`, the current ORPHEUS
apply matvec) the source collapses to

.. math::
   :label: hebert-3-432-source

   \bar Q_i \;=\; \sum_\ell \frac{2\ell+1}{2}\,Q_\ell(r_i)\,(-1)^\ell
   \;\;\xrightarrow{L=0}\;\;
   \frac{1}{2}\,\Sigma_t(r_i)\,\phi_0(r_i),

where :math:`\phi_0` is the scalar-flux Legendre :math:`\ell = 0`
moment of the input :math:`\psi`.  The :math:`L = 0` collapse is the
*apply matvec's* reach (isotropic scattering); the **source** side is
NOT collapsed — Issue #282 route (a) folds **all** Legendre moments of
the true within-group source, because streaming manufactures angular
structure an isotropic flux does not have (an :math:`\ell = 0`-only
fold floored the anisotropic curvilinear MMS).  See
:ref:`sn-282-source-fold` for the full fold and the load-bearing
:math:`\ell = 1` term.

Discretising Eq. :eq:`hebert-3-432` on a sub-mesh of cell width
:math:`\Delta r_i` gives the DD cell-balance Hébert Eq. (3.433):

.. math::

   -\bigl(\bar\phi_{i+1/2} - \bar\phi_{i-1/2}\bigr)
   \;+\; \Delta r_i \cdot \Sigma_i \cdot \bar\phi_i
   \;=\; \Delta r_i \cdot \bar Q_i,

with Hébert's typographic conventions

.. math::

   \bar\phi_i \;\equiv\; \phi_{1/2,\,i}, \qquad
   \bar Q_i \;\equiv\; Q_{1/2,\,i}, \qquad
   \Delta r_i \;=\; r_{i+1/2} - r_{i-1/2}.

The negative sign on the streaming jump comes from
:math:`\mu = -1 < 0` — particles travel **inward**, so the
discrete jump is :math:`-(\phi_{i+1/2} - \phi_{i-1/2})`.
**Critically**, no :math:`\alpha`-redistribution divisor appears
in this balance because :math:`(1 - \mu^2) = 0` at the endpoint.
This is the entire reason Hébert can solve the :math:`\mu = -1`
sweep in closed form with a plain DD recurrence: the coupled
angular cascade is decoupled at the starting direction.

Combining the DD auxiliary relation
:math:`\phi_{n,i} = \frac{1}{2}(\phi_{n,i-1/2} + \phi_{n,i+1/2})`
specialised to the :math:`-1/2` ordinate with the balance and
solving for :math:`\bar\phi_i` in terms of the known
outgoing-face value :math:`\bar\phi_{i+1/2}` (further from the
centre — known because we sweep **inward** from the outer BC)
yields Hébert Eq. (3.434):

.. math::
   :label: hebert-3-434

   \bar\phi_i \;=\; \frac{\Delta r_i \cdot \bar Q_i
                            + 2 \cdot \bar\phi_{i+1/2}}
                          {\Delta r_i \cdot \Sigma_i + 2}.

Stepping inward to the next face uses the textbook DD auxiliary
relation rearranged (Hébert Eq. (3.435)):

.. math::
   :label: hebert-3-435

   \bar\phi_{i-1/2} \;=\; 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}.

The pair :eq:`hebert-3-434`–:eq:`hebert-3-435` IS the spatial
recurrence.  Together they realise a tridiagonal-style inward
sweep on the radial mesh: outer face :math:`\rightarrow` cell
centre :math:`\rightarrow` inner face :math:`\rightarrow` next
cell centre :math:`\rightarrow \ldots \rightarrow` pole face
:math:`\bar\phi_{1/2}` at :math:`r = 0`.

.. note::

   The three labels :eq:`hebert-3-432-source`,
   :eq:`hebert-3-434`, :eq:`hebert-3-435` are also declared in the
   :mod:`~orpheus.sn.sweep.psi_half_angle_seed` module docstring
   (the canonical algebra-of-record).  Each :math:`:label:` is
   unique across the documentation graph; the Sphinx page is the
   **presentation layer** for the equations the code module owns
   as source-of-truth.  The
   ``@pytest.mark.verifies("hebert-3-43X")`` wiring on the L0
   algebraic-identity tests in
   :file:`tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py` is tracked
   at Issue #194; without that wiring the labels appear in the V&V
   audit as "documented but not tested" (orphan labels).

Why :math:`\mu = -1` is the natural starting direction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The M-M angular closure on sphere is a per-cell
:math:`\alpha`-cascade that **couples** the angular flux across
ordinates within one spatial cell: ordinate :math:`n` reads
:math:`\alpha_{n-1/2}` from the previous (more-inward-:math:`\mu`)
ordinate.  To start the cascade at the smallest-:math:`\mu`
ordinate, one needs a value for :math:`\alpha_{1/2}` AND for the
angular-edge flux :math:`\phi_{1/2,i}` at that seed half-integer.

The :math:`\alpha_{1/2} = 0` seed is **free**: it comes from
:math:`1 - \mu^2` evaluated at :math:`\mu = -1`, i.e.,
"The first value :math:`\alpha` is equal to :math:`1 - (-1)^2 = 0`"
(text below Hébert Eq. (3.422)).  That handles the *angular*
half of the problem.

The flux value :math:`\phi_{1/2,i}`, however, is NOT free.  It
is the **spatial flux profile** at the auxiliary starting
direction, and it must be solved for as a function of position
:math:`i` along the radial mesh.  Eqs. :eq:`hebert-3-432` through
:eq:`hebert-3-435` provide exactly that spatial solve.

At :math:`\mu = -1` the sphere streaming operator collapses to
pure radial divergence **without** the angular-redistribution
coupling.  As Hébert writes (p. 143):

   *"We observe that these directions correspond to particles
   entering the external surface and moving toward the central
   axis with* :math:`\mu = -1`. *The angular redistribution term
   vanishes on these points so that Eq. (3.164) simplifies to
   [Eq. (3.432)]."*

This is the **only** direction on the unit interval :math:`[-1, 1]`
where the spatial 1D-sphere problem reduces to a closed-form
linear recurrence in radius alone, without an inner angular
solve.  Picking any intermediate :math:`\mu` would leave the
coupling term active and re-introduce the cascade
chicken-and-egg.  See also
:ref:`sn-phase-d-pomraning-structural-singularity` for the
deeper structural reason :math:`\mu = \pm 1` is the only
admissible starting direction in any curvilinear geometry.

Why "zero-weight"
~~~~~~~~~~~~~~~~~

In an :math:`N`-point Gauss–Legendre quadrature on :math:`[-1, 1]`
the endpoints :math:`\mu = \pm 1` are **not** base points (the
polynomial is approximated by interior nodes only).  They have no
quadrature weight, hence "zero-weight" — the flux value at
:math:`\mu = -1` does NOT contribute to any
:math:`\sum_n \mathcal{W}_n \phi_n` integral that builds the scalar
flux moments.

The :math:`\mu = -1` ordinate is therefore a **purely auxiliary
numerical construct**: its flux values exist for the sole purpose
of seeding the :math:`\alpha`-cascade for the finite-weight
ordinates that follow.  After the cascade is initialised, the
angular-edge values :math:`\bar\phi_{i\pm 1/2}` are discarded;
only the **cell-centred values** :math:`\bar\phi_i \equiv
\phi_{1/2,i}` are kept (Hébert, p. 143, between Eqs. (3.435) and
(3.436)).  Those cell-centred values feed the finite-weight
ordinates' cell-balance Eq. (3.436) via the
:math:`(\alpha_{n-1/2} + \alpha_{n+1/2})\,\phi_{n-1/2,i} /
(2\,\mathcal{W}_n)` redistribution term, with
:math:`\phi_{n-1/2,i} = \phi_{1/2,i}` at the first
finite-weight ordinate :math:`n = 1`.

The flat-:math:`\psi` algebraic verification trace
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Phase D hypothesis is: *for a flat angular flux*
:math:`\psi_{\text{cell}} = C` *across all cells, the
inward sweep returns* :math:`\bar\phi_{1/2} = C`.  The algebra
verifies this in closed form.

Take a homogeneous problem with constant :math:`\Sigma_t = \Sigma`
and source :math:`\bar Q_i` constructed so the consistent fixed
point is :math:`\bar\phi_i = C` everywhere.  Specialising
Eq. :eq:`hebert-3-432-source` to :math:`L = 0` and applying the
flat-:math:`\psi` ansatz gives the consistent source
:math:`\bar Q = \frac{1}{2} \Sigma \cdot 2C = \Sigma \cdot C` (the
:math:`\phi_0` integral over flat unit-:math:`\psi` against GL
weights summing to 2 returns :math:`2C`; lumped into the discrete
:math:`\bar Q_i = \Sigma \cdot C`).

Substituting into Eq. :eq:`hebert-3-434` with inductive hypothesis
:math:`\bar\phi_{i+1/2} = C`:

.. math::

   \bar\phi_i
   \;=\; \frac{\Delta r \cdot \Sigma \cdot C + 2C}
              {\Delta r \cdot \Sigma + 2}
   \;=\; C \cdot \frac{\Delta r \cdot \Sigma + 2}
                     {\Delta r \cdot \Sigma + 2}
   \;=\; C.

Eq. :eq:`hebert-3-435` then gives
:math:`\bar\phi_{i-1/2} = 2C - C = C`.  The recurrence is
self-similar: every face and cell value stays at :math:`C`.  Hence
:math:`\bar\phi_{1/2}(r = 0) = C` for flat :math:`\psi` on the
consistent flat source — **the hypothesis holds**.

This trace establishes the Phase D fix as a **closed-form
analytical reference** in the
:doc:`algebra-of-record </development>` State-1A pillar sense: the
identity :math:`(L \cdot \psi_{\text{flat}})_{n,i,g} = \Sigma_t
\cdot \psi_{n,i,g}` is verifiable by exact algebra on the discrete
operator, no numerical quadrature required.  The L0 foundation test
:func:`tests.sn.sweep.curvilinear.test_psi_half_angle_seed.TestCarlsonFlatPsiAlgebraicIdentity.test_carlson_flat_psi_identity_reflective`
pins this identity at machine precision (``rtol=1e-13``).

The corrected injection-point story
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The single largest **architectural correction** of Phase D is
where the canonical inward-sweep output is injected.  The Phase D
plan (and the literature memo's §7 implementation note) routed
the inward-sweep result :math:`\bar\phi_i` into the **WDD
spatial pole-face initial condition** at the then-production
``transport_operator_matvec_spherical`` matvec's (since deleted)
``psi_face_in`` initialisation — the very same site the
:ref:`sn-curvilinear-trajectory-resolvent-crosscheck-section` discussion
identified as the Phase C Carlson seed location.

The numerics-investigator diagnostic
(:file:`tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`)
falsified that hypothesis empirically.  Four interventions tested
against the M-M failing configuration on the flat-:math:`\psi`
probe:

.. list-table:: Phase D injection-point intervention sweep (Σ_t = 0.5)
   :header-rows: 1
   :widths: 8 40 30 22

   * - Probe
     - What it changes
     - Site
     - max\|residual\|
   * - ``[A]``
     - Carlson seed for WDD ``psi_face_in``
     - ``operator.py:738``
     - **1.89e+01 FAIL** (unchanged)
   * - ``[B]``
     - Carlson seed for M-M half-angle ``ψ_{1/2,i}``
     - ``pole_angular_closure.py:411``
     - **1.78e-15 PASS**
   * - ``[C]``
     - BOTH ``[A]`` + ``[B]``
     - both
     - **1.78e-15 PASS** (no extra effect)
   * - ``[D]``
     - M-M half-angle ``ψ_{1/2,i}`` = cell-centre value
     - ``pole_angular_closure.py:411``
     - **1.78e-15 PASS** (degenerate)

Reading the table:

* ``[A]`` confirms the WDD spatial pole-face IC is **not** what's
  wrong.  The Phase C
  ``psi_face_in = fi[:, outgoing_mask, 0, 0]`` Lewis–Miller
  cell-centre seed is already structurally equivalent to the
  Carlson inward-sweep output **on flat ψ** — both equal
  :math:`\psi_{\text{cell}}[0]` in that limit.  Replacing the
  WDD seed changes nothing.
* ``[B]`` is the canonical Carlson intervention: feeding
  :math:`\bar\phi_i` into the M-M recurrence's ``psi_half_left``
  closes the residual to machine precision.
* ``[D]`` is the **falsification check**: on the flat-:math:`\psi`
  reflective probe ``[B]`` and ``[D]`` coincide because the
  inward sweep returns :math:`\bar\phi_i \equiv \psi_{\text{cell}}`
  exactly.  The probe **cannot distinguish** the two.

To prove the Carlson seed is canonical (not merely coincidentally
correct), the diagnostic includes a vacuum-BC structural
independence cross-check.  On a vacuum-BC probe the inward sweep
returns a non-trivial spatial profile

.. math::

   \bar\phi_i \;=\; (0.613, 0.572, 0.527, 0.478, 0.423, 0.362,
                     0.295, 0.220, 0.138, 0.048),

distinctly **not** equal to the cell-centred flat
:math:`\psi_{\text{cell}} = \mathbf{1}`.  The two seeds differ by
up to 0.95 in absolute value, and the resulting operator residuals
differ by max-abs 7.31 — the Carlson seed ``[B]`` is mathematically
distinct and quantitatively superior to the degenerate
broadcast-cell-centre seed ``[D]``.  This is the
**structural-independence evidence** that pins the Phase D fix as
canonical, not as a coincidental match on a degenerate probe.

The pinning test for this structural distinction is
:func:`tests.sn.sweep.curvilinear.test_psi_half_angle_seed.TestCarlsonFlatPsiAlgebraicIdentity.test_carlson_vacuum_BC_flat_source_nx_3`
— a vacuum-BC hand calculation on the Carlson inward sweep
(``rtol=1e-13``) whose values are distinct from the degenerate
broadcast-cell-centre seed.  Without this test a future
regression that replaced the Carlson sweep with a naive
broadcast-cell-centre would pass every flat-:math:`\psi`
reflective test silently.

The bug Phase B baked in
~~~~~~~~~~~~~~~~~~~~~~~~

The pre-Phase-D production code at
``orpheus/sn/sweep/pole_angular_closure.py:411`` carried the
hardcoded zero seed:

.. code-block:: python

   psi_half_left = np.zeros((ng, nx), dtype=psi_level.dtype)
   for m in range(M):
       tau_m = tau_level[m]
       psi_half_right = (
           psi_level[:, m, :] - (1.0 - tau_m) * psi_half_left
       ) / tau_m
       redist[:, m, :] = (
           dAw_level[:, m].reshape(1, nx)
           * (alpha_level[m + 1] * psi_half_right
              - alpha_level[m] * psi_half_left)
           / volume.reshape(1, nx)
       )
       psi_half_left = psi_half_right

The Phase B docstring justified the zero seed as: *"for the
forward apply matvec we adopt* :math:`\phi_{1/2,i} = 0`, *the
unique choice that makes the recursion's seed consistent with*
:math:`\alpha_{1/2} = 0` *and that the sweep converges to under
fixed-point iteration."*  This reasoning is wrong — the
:math:`\alpha_{1/2} \psi_{1/2}` product vanishes regardless of
:math:`\psi_{1/2}` because :math:`\alpha_{1/2} = 0`, but the seed
ALSO enters the **denominator-propagation chain**: every
subsequent half-angle face flux
:math:`\psi_{m+1/2,i,g}` depends on :math:`\psi_{m-1/2,i,g}`
recursively, and the chain inherits the seed through the M-M
weighting :math:`(1 - \tau_m)`.  Setting the seed to zero when
Hébert's structural form says :math:`\psi_{1/2,i,g} =
\bar\phi_{1/2,i}` (the inward-sweep output) is a **wrong term
initialisation** — Mode 3 in the
``vv-principles`` 6-failure-mode taxonomy (see ``error_catalog.md``
ERR-026 entry).

How the wrong seed survived Phase B
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The zero seed survived Phase B's L1 flat-flux-identity test
(:file:`tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py`)
because that test compared the three closures (Legacy / BFF /
M-M) **against each other on flat ψ**, NOT against the closed-form
fixed-point identity :math:`L \cdot \psi = \Sigma_t \cdot \psi`.
All three closures collapse to the same wrong-but-internally-
consistent value on flat :math:`\psi`, so cross-comparison passes
while the absolute closed-form check would have caught it
immediately.

The cylindrical case ALSO carries the zero seed in production but
:ref:`Cylindrical Gate 1.1 <sn-phase-d-gate-1-1-empirical>`
**passes** empirically.  The mechanism is the **dead first-ordinate
seed** of the level-symmetric quadrature exercised here: the
first-swept ordinate's seed weight is zero
(:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`), so
the wrong ``psi_half_left = 0`` seed is annihilated at source per
level.  (This was originally read as per-:math:`\mu`-level
:math:`\alpha`-dome telescoping "cancelling" the seed; #280 Phase
2.5b corrected that — it is a dead weight, level-symmetric-only, and
**false for a product quadrature**, where the seed is a live
self-coupling and the cold solve was seed-lagged until the
direct-seed fold.)  The sphere cascade has no equivalent dead-seed
weight — a wrong seed propagates directly to a wrong fixed point.  Phase D's fix
updates the cylindrical path too for **structural alignment with
the canonical form** (architectural correctness), but cylindrical
behaviour is empirically a regression-stability check, not a new
PASS.

.. _sn-phase-d-pomraning-structural-singularity:

Pomraning structural-singularity cross-reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pomraning (1989) [Pomraning1989]_ frames the curvilinear pole
problem as **geometric**: :math:`r = 0` is structurally singular
in any curvilinear streaming operator because the
angular-derivative coefficients in the streaming term (the
:math:`(1 - \mu^2)/r` factor in the sphere streaming operator)
contain :math:`1/r`.  At :math:`r = 0` the coefficient diverges;
the natural discretisation must somehow handle this.  In his words
(p. 339, right column):

   *"It was pointed out that if the bounding surface of the
   system is used as one of the coordinate surfaces and one
   considers a family of nonintersecting surfaces that starts
   with the bounding surface and progresses inward to fill the
   system, then these surfaces will eventually shrink to a
   surface with a zero area, namely a line or a point. ... A
   special case of this elliptical example is a sphere, where
   the innermost surface is simply a point.  Hence, in general
   there will exist points on the innermost surface where the
   coefficients of the angular derivatives in the streaming term
   are infinite, since these coefficients contain the reciprocal
   of the radii of curvature ... Prime examples of such singular
   points are found in the usual spherical and cylindrical
   geometry formulations where* :math:`1/r` *terms are extant and
   the attendant difficulties are well known, particularly in
   numerical treatments."*

The naive engineering response would be **extrapolation**: pick
:math:`\psi_{\text{face}}(r = 0)` by fitting a polynomial in
:math:`r` through nearby interior cells.  This is what an
incautious starting heuristic does; it is also what produces
the M-M wrong fixed point ERR-026 diagnoses.

The Carlson coupled-pole response is **canonical** because it
sidesteps the singularity entirely: at the auxiliary direction
:math:`\mu = -1` the singular :math:`(1 - \mu^2)/r` term is
**identically zero** (the numerator vanishes), so the spatial
sweep at this direction sees **no singularity at all**.  The
equation tells the discretisation what
:math:`\bar\phi_{1/2}(r = 0)` should be — there is no need to
guess.  The cost is that the :math:`\mu = -1` sweep must be
solved first, then its result used as the seed for the cascade
at finite-weight ordinates (where :math:`(1 - \mu^2) > 0` and
the singularity would otherwise be felt).  This is exactly the
price Pomraning warns about: "difficulties must be dealt with".
The Carlson construction deals with it by **exploiting** the
singularity's vanishing at :math:`\mu = \pm 1` rather than
trying to regularise it at intermediate :math:`\mu`.

Option α: composition over sibling Protocol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The seed is **M-M-specific**: only
:class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
carries a ``psi_half_left`` variable to seed.  The Legacy and
Bailey closures don't have one — their half-angle face flux
evaluation collapses to cell-centre values unconditionally.  Two
architectures were considered:

* **Option α (composition, shipped)** — the seed strategy lives as
  a ``PsiHalfAngleSeed``
  field on
  :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`.
  The abstraction stays local to the closure that consumes it.

* **Option B (sibling Protocol on SNMesh, rejected)** — the seed
  would be a separate Protocol attribute on
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`, applied by the matvec
  before calling the pole closure.  This would force every
  consumer (Legacy / BFF / M-M) to handle a Protocol that is a
  **no-op** for the non-M-M strategies, violating the
  single-responsibility principle and forcing unrelated tests to
  thread the Protocol through call signatures.

The
``CarlsonSweepContext``
dataclass bundles the four inputs the Carlson sweep needs that
are NOT in the
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
strategy's ordinary per-cell call signature (``sigma_t``, ``dr``,
``mu_quad``, ``weights``, ``bc_outer_value``), keeping the
call-signature expansion to a single new optional keyword — a minimal
blast-radius extension that Legacy and Bailey closures ignore by
documented closure contract.

Linear-operator preservation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both seed strategies — ``ZeroSeed`` and
``CarlsonInwardSweep`` — are **linear in the input** ``psi_cells``
(verified by the ``is_linear: ClassVar[bool] = True`` trait, pinned
by foundation tests).  Linearity is the load-bearing property:
the apply matvec must be a linear operator, otherwise the
operator-algebra operations of
:class:`~orpheus.sn.operators.streaming.InvertibleOperator`
(apply, apply_transpose, dense matrix probing) break.  The
``CarlsonInwardSweep`` is linear because:

* The :math:`\phi_0` moment is a linear projection of input
  :math:`\psi` (Legendre integration is linear).
* :math:`\bar Q = \frac{1}{2} \Sigma_t \cdot \phi_0` is linear
  in :math:`\psi` (:math:`\Sigma_t` is constant).
* The recurrence Eqs. :eq:`hebert-3-434`–:eq:`hebert-3-435` is
  an affine function of :math:`(\bar Q, \bar\phi_{i+1/2})` with
  constant coefficients depending only on
  :math:`(\Sigma_t, \Delta r)`.
* The ``bc_outer_value`` is constructed in the matvec by applying
  the realised BC operator to the cell-centred outer-cell
  :math:`\psi`, then extracting the most-inward ordinate's value
  — both operations are linear in the input :math:`\psi`.

The foundation test
:func:`tests.sn.sweep.curvilinear.test_psi_half_angle_seed.TestSeedLinearity.test_carlson_inward_sweep_is_linear`
pins the linearity directly; the operator-level linearity gate
in :file:`tests/sn/test_streaming_operator.py` pins it transitively
at the matvec boundary (``rtol=1e-12`` — relaxed from the
pre-Phase-D ``rtol=1e-13`` to absorb ~10×ULP non-associativity
drift, justified by the three principled-relaxation criteria of
the ``vv-principles`` bit-identity-vs-principled-equivalence
framework).

The L = 0 isotropic-only limitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current
``CarlsonInwardSweep``
evaluates only the :math:`\ell = 0` (isotropic) Legendre moment
when building the moment-folded source in
:eq:`hebert-3-432-source`.  This is **consistent with the apply
matvec's structure**: the
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` apply matvec
carries only an isotropic collision term :math:`\Sigma_t \psi`;
anisotropic scattering (P\ :sub:`1`\ +) is composed externally via
a separate scattering operator, not included in :math:`L`.

.. warning::

   **The L = 0 isotropic-only assumption is load-bearing for the
   Phase D fix.**  If a future refactor moves scattering INTO
   :math:`L` (e.g., to enable a "monolithic" SN apply that
   includes within-group scattering), the Carlson seed becomes
   WRONG: the source at :math:`\mu = -1` (Eq. :eq:`hebert-3-432`)
   needs the full Legendre-moment sum

   .. math::

      \bar Q_i \;=\; \sum_\ell \frac{2\ell+1}{2}\,Q_\ell(r)\,(-1)^\ell,

   not just :math:`\Sigma_t \phi_0`.  This is a Mode-6
   convention-drift risk per the ``vv-principles`` skill (the
   definition-site assumption disagreeing with the usage-site
   intention).  A foundation test pinning the isotropic-only
   assumption (e.g., asserting the apply matvec does NOT couple
   to ``self_scattering``) would catch a future drift; in its
   absence, this WARNING block and the module docstring's
   matching admonition are the only safeguards.  Track the
   future-refactor case under a fresh GitHub issue when the
   monolithic apply work is scheduled.

.. _sn-phase-d-default-flips:

Default flips
~~~~~~~~~~~~~

Phase D ships **two default flips** that activate the full
canonical curvilinear closure path:

#. :attr:`SNMesh.pole_angular_closure
   <orpheus.sn.mesh.augmented_mesh.SNMesh.pole_angular_closure>` default
   flipped from
   ``LegacyTauSymmetricInterpolation``
   to
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`.
   :class:`MorelMontryAngularSweep`'s own constructor default for
   ``psi_half_seed`` is
   ``CarlsonInwardSweep``,
   so the single :class:`SNMesh` flip activates the full Phase D
   fix (canonical M-M closure + canonical Carlson seed) without
   requiring downstream call sites to thread the new strategy
   explicitly.

#. :class:`~orpheus.sn.solver.SNSolver`'s ``inner_solver`` default
   flipped from ``"source_iteration"`` to ``"krylov"`` for
   **curvilinear geometries** (spherical, cylindrical); Cartesian
   stays at ``"source_iteration"``.  The rationale: the Phase D
   fix lives in the apply matvec, and the Krylov path is the one
   that uses
   :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply`.  The
   sweep path (``"source_iteration"``) uses the spatial WDD
   recurrence and is unaffected by the Phase D fix — leaving its
   ERR-026-affected curvilinear behaviour in place would be wrong
   for the production default.

   .. note:: **Reverted (2026-06-12, Issue #195).**

      This curvilinear ``"krylov"`` default was undone by the ERR-058
      fix.  The premise — that the sweep path was ERR-026-affected
      while Krylov-on-apply was not — held only because the sweep and
      matvec were *distinct* discrete systems at Phase D time.  After
      the Depth-B/Wave-T unification they are ONE system, and the
      ERR-058 closure-seed fix makes that system O(h²)-consistent: SI
      :math:`\equiv` Krylov bit-identical, SI :math:`\sim 10^2\times`
      faster.  The curvilinear default returned to
      ``"source_iteration"``.  See
      :ref:`sn-err-058-closure-seed-closeout`.

.. _sn-phase-d-gate-1-1-empirical:

Empirical Gate 1.1 outcome (Phase D — full 12-cell crosstab)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Phase D acceptance gate is Gate 1.1 on **all three** pole
closures across both curvilinear geometries and both :math:`\Sigma_t`
values.  The parametrised test
:func:`tests.sn.sweep.core.test_phase_c_gates.test_apply_curvilinear_per_ordinate_flat_flux_residual`
produces the 12-cell crosstab:

.. list-table:: Gate 1.1 outcome under Phase D Carlson seed (2026-05-12)
   :header-rows: 1
   :widths: 18 30 26 26

   * - Geometry
     - Pole closure
     - :math:`\Sigma_t = 0`
     - :math:`\Sigma_t = 0.5`
   * - Sphere
     - ``LegacyTauSymmetricInterpolation``
     - PASS
     - PASS
   * - Sphere
     - ``BaileyFlatFluxRedist``
     - PASS
     - PASS
   * - Sphere
     - ``MorelMontryAngularSweep``
     - **XPASS**
     - **XPASS**
   * - Cylinder
     - ``LegacyTauSymmetricInterpolation``
     - PASS
     - PASS
   * - Cylinder
     - ``BaileyFlatFluxRedist``
     - PASS
     - PASS
   * - Cylinder
     - ``MorelMontryAngularSweep``
     - **XPASS**
     - **XPASS**

All 12 cells PASS or XPASS.  The 4 XPASS cells under M-M closure
are the ERR-026 markers — they now flip from FAIL to XPASS on
xfail-strict=False, unblocking the marker-removal commit that
Phase D Step 5 will execute (deferred per the closeout memo's
acceptance gate item 6).

This is the load-bearing **empirical evidence** for the
ERR-026 identity-and-rate scope closure.  The asymmetry between
the Phase C (cylinder PASS / sphere FAIL) and Phase D
(both PASS) crosstabs is the diagnostic mark of the Phase D
intervention: the sphere case required the Carlson seed because
its single-cascade structure has no telescoping; the cylinder case
already passed under Phase C because — for the **level-symmetric**
quadrature exercised here — the first-swept ordinate's seed weight
is exactly zero (:math:`c_{\rm in}[m_0] = (1-\tau)/\tau = 0` at raw
:math:`\tau = 1`), so the zero-seed inconsistency was annihilated at
source (a **dead** seed), not "absorbed" by any telescoping of the
solve.

.. note::

   **Correction (#280 Phase 2.5b, 2026-07-05).**  An earlier reading
   of this crosstab attributed the cylinder's Phase-C pass to
   ":math:`\alpha`-dome telescoping absorbing the zero-seed
   inconsistency" and generalised it to "the cylinder solve is
   seed-insensitive / was already exact."  That is a **level-symmetric-
   only** artefact and is **false for a product quadrature**: there the
   starting direction coincides with the first-swept ordinate
   (:math:`\mu_{\rm start} \equiv \mu_{m_0}`, :math:`t = 0`, #229), so
   :math:`c_{\rm in}[m_0] \ne 0` and the seed is a **live per-ordinate
   self-coupling** that contributes :math:`O(1)` to the :math:`m_0`
   cell diagonal.  The product-cylinder cold ``(L+C).solve`` was in
   fact seed-**lagged** (cold error :math:`\approx 0.57`) until the
   #280 2.5b direct-seed fold folded that self-coupling into the
   :math:`m_0` diagonal (:math:`c_{\rm out} \to c_{\rm out} -
   c_{\rm in}`), making the cold solve a single-pass direct inverse.
   The augmented :math:`(L+C)` is block-lower-triangular because the
   seed contribution lands **on the block diagonal** (forward
   substitution resolves it) — *not* because the seed "telescopes
   away."  Distinct claim, still valid: the :math:`\alpha`-dome
   telescopes under the angular weight sum
   :math:`\sum_n w_n \psi_n`, which is why **scalar / balance** V&V
   gates are blind to a wrong per-ordinate seed (anti-pattern #8) —
   that blindness statement is unaffected by this correction.

.. _sn-phase-d-gate-1-5-capture-and-compare:

Gate 1.5 strengthening — capture-and-compare BC apply input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase C's Gate 1.5
(:ref:`bc-trace-contract-respected-by-matvec`) was a "round-trip"
check: invoke ``bc.realize().apply(...)`` independently and
compare against the matvec's observable output.  Phase D
strengthens this to a **capture-and-compare** check that pins the
exact value the matvec passes into the BC trace law:

#. Patch ``sn_mesh.bc["xmax"].apply`` (the outer radial face —
   a sphere's ``"outer"`` endpoint renders as ``"xmax"`` since
   C4 / #220, see :ref:`bc-face-name-carve`) to capture every input
   array passed to it during one matvec call.
#. Independently reconstruct the WDD-propagated outflow trace via
   a reference implementation (:func:`_outflow_at_boundary_for_sphere`).
#. Assert the captured BC apply input matches the reference to
   ``rtol=1e-14`` — exactly bit-equal up to FP non-associativity.

The strengthening matters because the Phase D matvec now calls
``bc["xmax"].apply`` **twice** per matvec:

#. **Phase D Carlson context call** — applied to cell-centred
   outer-cell :math:`\psi` to build ``bc_outer_value`` for the
   ``CarlsonSweepContext``.  See the BC companion section
   :ref:`bc-phase-d-two-bc-applies-per-matvec`.
#. **Phase C BC trace law call** — applied to the WDD-propagated
   outflow face value at the boundary edge, per the
   :ref:`affine-bc-form` contract.

The capture-and-compare test
:func:`tests.sn.sweep.core.test_phase_c_gates.test_bc_trace_contract_capture_and_compare_sphere`
(parametrised over ``vacuum`` and ``reflective``) **locates the
Phase C call by shape and content matching**: of the two captured
inputs, the one whose shape matches ``(N, ng)`` and whose values
match the independent reference is the Phase C trace law call.
Both vacuum and reflective parametrised cases pass; the test is
foundation-tagged because it pins a software invariant (the
matvec's two-application sequence) rather than a math claim.

.. _sn-phase-d-err-026-closure-narrative:

ERR-026 PARTIAL → PARTIAL (narrowed scope)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Retraction (2026-06-12, Issue #195).**

   The sub-claim table below classified the residual MMS magnitude as
   a benign "pre-asymptotic transient" that finer :math:`n_x` would
   clear (status OPEN, tracked at #195).  **That classification is
   wrong.**  The curvilinear isotropic MMS error did not shrink under
   refinement — it PLATEAUED mesh-independently (orders :math:`\to 0`),
   because the dominant defect was the *angular* closure seed (the
   Carlson proxy source), not a spatial-truncation constant.  The
   "rate :math:`[3.33, 2.46]` already correct, only the constant is
   large" reading was an artefact of the
   :math:`\alpha`-dome-telescoping blindness of the scalar residual.
   ERR-058 (#195) replaced both seeds; the isotropic MMS is now a
   clean :math:`\mathcal{O}(h^2)` ladder.  The table is preserved as
   bug-era evidence — its STATUS / Closed-by interpretation is
   superseded by :ref:`sn-err-058-closure-seed-closeout`.

Phase D **narrows** ERR-026's open scope.  The bug ERR-026
originally diagnosed — *"curvilinear sweep WDD angular closure
converges to wrong fixed-source solution"* — had three
sub-claims, each addressed by a different Wave:

.. list-table:: ERR-026 sub-claim closure tracking
   :header-rows: 1
   :widths: 35 35 30

   * - Sub-claim
     - Status
     - Closed by
   * - Operator identity:
       :math:`(L \cdot \psi_{\text{flat}})_{n,i,g} = \Sigma_t \cdot \psi_{n,i,g}`
       on per-ordinate flat-flux probe
     - **CLOSED**
     - Phase D Carlson seed (Gate 1.1 XPASS)
   * - Convergence rate:
       :math:`\mathcal{O}(h^2)` MMS rate at fixed :math:`N`
     - **CLOSED (rate)**
     - Phase D Carlson seed (empirical rate [3.33, 2.46] across
       refinements; both above the L1 acceptance floor of 1.9)
   * - Convergence magnitude: pre-asymptotic absolute MMS error
       below quadrature floor at practical ``nx`` (:math:`\le 160`)
     - **OPEN**
     - Tracked at `Issue #195
       <https://github.com/deOliveira-R/ORPHEUS/issues/195>`_;
       requires either finer ``nx`` or a higher-order spatial
       closure refinement to fully close

The convergence-rate evidence ``[3.33, 2.46]`` is the slope
sequence measured at successive refinement levels; both values are
above the L1 acceptance floor of 1.9 (second-order accuracy
demonstrated robustly), satisfying the rate sub-claim.  However,
the **absolute magnitude** at the largest tested ``nx`` (=160)
remains above the L1 tolerance that the test architect specified
for full closure on the pre-asymptotic regime.  This is **NOT** a
violation of the Phase D fix — the rate is correct, the asymptotic
regime is the right shape, but the **constant-coefficient** in
front of the :math:`\mathcal{O}(h^2)` term is larger than the
test's pre-asymptotic-magnitude budget at practical mesh
resolutions.

The pre-asymptotic regime is the consequence of the Carlson
sweep's L0-truncated source: at coarse :math:`nx` the Legendre
:math:`\phi_0` moment is computed from the cell-centred input
:math:`\psi` against the GL quadrature — an integration whose own
truncation contributes to the constant in
:math:`\bar Q = \frac{1}{2} \Sigma_t \phi_0`.  Refining
:math:`nx` reduces this contribution, but the rate at which it
reduces is set by the WDD spatial closure's own truncation order,
not by the Carlson sweep itself.  See Issue #195 for the candidate
follow-up paths (higher-order pole-face spatial closure, or a
:math:`\phi_0` recomputation that uses the M-M angular recurrence
output rather than the cell-centred input).

The 4 ``xfail-strict`` ERR-026 tripwires
(:file:`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`,
sphere + cylinder × isotropic + anisotropic ansatz) therefore
**stay xfail** through Phase D Step 3.  They will ``xpass`` under
the Phase D defaults (which is what triggers the deferred Step 5
marker-removal commit); the pre-asymptotic-magnitude regression
that prevents `strict=True` flipping is Issue #195's domain.  The
narrative for ``error_catalog.md`` therefore reads:

   ERR-026 status: **PARTIAL CLOSURE** (was PARTIAL through Phase
   C, narrowed scope through Phase D).  The structural bug (M-M
   recurrence hardcoded ψ\ :sub:`1/2,i` = 0 seed) is closed by the
   Phase D Carlson coupled-pole sweep; Gate 1.1 sphere MMS PASS
   confirms the operator identity and the second-order
   convergence rate is recovered.  The pre-asymptotic-magnitude
   open question (Issue #195) is what keeps the status at PARTIAL
   rather than CLOSED.

.. _sn-phase-d-files-touched:

Files touched by Phase D
~~~~~~~~~~~~~~~~~~~~~~~~

The full Phase D footprint (per the closeout memo at
``.claude/agent-memory/method-implementer/issue_168_phase_d_step3_closeout.md``):

**New modules**

* :mod:`orpheus.sn.sweep.psi_half_angle_seed` — Protocol family
  + ABC + 2 strategies (``ZeroSeed`` + ``CarlsonInwardSweep``)
  + ``CarlsonSweepContext`` dataclass.
* :file:`tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py` — 25
  foundation + L0 + L1 tests covering Protocol conformance,
  registry/self-registration, immutability, shape contract,
  bit-identity for ``ZeroSeed``, L0 algebraic identities
  (flat-:math:`\psi` at varying C, vacuum-BC nx=3 hand
  calculation, multi-region :math:`\Sigma_t` step), linearity, and
  L1 structural-independence (Carlson vs Zero on vacuum-BC probe).

**Modified files**

* :mod:`orpheus.sn.sweep.pole_angular_closure` —
  :class:`MorelMontryAngularSweep` gains a
  ``psi_half_seed: PsiHalfAngleSeed`` field; the per-level M-M
  recurrence (then ``_mm_weighted_angular_recurrence_single_level``,
  now :func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`)
  accepts an
  optional ``psi_half_seed`` array; Protocol signatures extended
  with an optional ``carlson_context`` kwarg (Legacy + Bailey
  ignore it).
* :mod:`orpheus.sn.sweep` ``__init__`` re-exports the new
  symbols.
* :mod:`orpheus.sn.operators.streaming` — spherical + cylindrical matvecs
  build the
  ``CarlsonSweepContext``
  before calling ``pole_angular_closure``.
* :mod:`orpheus.sn.mesh.augmented_mesh` — :class:`SNMesh` default flipped to
  :class:`MorelMontryAngularSweep`.
* :mod:`orpheus.sn.solver` — curvilinear default ``inner_solver``
  flipped to ``"krylov"``.
* :file:`tests/sn/test_phase_c_gates.py` — Gate 1.5 strengthened
  with capture-and-compare.
* :file:`tests/sn/test_streaming_operator.py` (post-D-K successor
  to the retired ``test_snstreamingoperator.py``) — 3 tests updated
  (one test docstring rewritten to pin the Phase D fix; two
  bit-identity tests threaded with ``sn_mesh.pole_angular_closure``;
  one linearity tolerance relaxed ``rtol=1e-13 → 1e-12``).

The agent-memory trail for Phase D session reproducibility:

* Literature memo:
  ``.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md``
  — Hébert §3.9.4 derivation + flat-:math:`\psi` algebra +
  architecture-shape correction + open questions.
* Diagnostic memo:
  ``.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md``
  — empirical evidence + 4 plan corrections + structural-
  independence cross-check.
* Step 3 closeout:
  ``.claude/agent-memory/method-implementer/issue_168_phase_d_step3_closeout.md``
  — what shipped + 3 deviations + V&V evidence chain.
* Diagnostic script:
  :file:`tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`
  — self-contained CLI probe reproducing the diagnostic table.


.. _sn-phase-f-carlson-sweep-path-backport:

Phase F Carlson seed sweep-path backport (Issue #168 Phase F)
-------------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   * Phase F (commit chain landed 2026-05-12 on
     ``refactor/sn-operator-algebra``, atop Phase E ``6708a4a``)
     backports the Phase D Carlson coupled-pole seed
     (``CarlsonInwardSweep``,
     Hébert §3.9.4 Eqs. (3.432)–(3.435)) from the apply-matvec path
     (the then-production ``transport_operator_matvec_spherical``
     / ``_cylindrical`` matvec — since deleted, #197 / #280 campaigns —
     fixed in Phase D Step 3) into the SI/sweep
     path
     (``_sweep_1d_spherical`` (the dissolved ``sweep.py``) and
     ``_sweep_1d_cylindrical``).
   * The bug is the **structural twin** of the Phase D defect: the
     SI loop in :file:`orpheus/sn/sweep.py` initialised
     ``psi_angle = np.zeros((nx, ng))`` at the spherical sweep
     entry (line 474, pre-Phase-F) and at the cylindrical per-level
     loop entry (line 634, pre-Phase-F) — the same hardcoded
     zero seed Phase D diagnosed as wrong-term-initialization on
     the apply-matvec twin
     (``orpheus/sn/sweep/pole_angular_closure.py:411``).
   * Phase F factors a **NEW free function**
     :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
     ``(Q_bar, sigma_t, dr, bc_outer_value) -> (ng, nx)`` that runs
     :eq:`hebert-3-434`–:eq:`hebert-3-435` driven by the SI
     within-group source ``Q_1d`` rather than by an apply-path
     :math:`\psi` Legendre fold.
     ``CarlsonInwardSweep.__call__`` is refactored to delegate
     to the same helper after folding ``psi_level → Q̄ = 0.5 ·
     Σ_t · φ_0``.  **One helper, two consumers** — Cardinal Rule 2
     (architecture) enforced via reuse without duplication.
   * Empirical result on ``sphere_2g_3reg`` n=40
     (heterogeneous A|B|A reflective 2-group sphere):
     ``sf[0]/sf[1]`` ratio at the pole was **0.522** (DIVERGING
     to **0.473** under refinement to n=320); post-Phase-F it is
     **0.778** and STABLE under refinement (still 0.777 at n=320).
     The outer-cell reflective-face defect ``sf[-1]/sf[-2]`` was
     **0.887** → **0.997** (essentially CLOSED).
     :math:`\psi(r=0)` quasi-isotropy: ``cv(ψ@i=0)``
     **0.520** → **0.404**, ``max/min(ψ@i=0)`` **6.4×** →
     **1.16×** (Pomraning 1989 prediction substantially approached).
   * **What was logged as open** *(now CLOSED, #196)*: the residual
     O(h) per-cell WDD spatial-closure asymmetry between SI and Krylov
     paths was logged as **manifestation #7 of ERR-026**.  It is now
     **CLOSED** — ERR-058 (#195) showed the gap was a shared
     closure-seed defect, not a discretisation asymmetry; #196 verified
     SI :math:`\equiv` Krylov to the iteration floor on the
     heterogeneous eigenvalue path and added the permanent regression
     gate (see :ref:`sn-phase-f-residual-o-h-open` and
     :ref:`sn-issue-196-eigenvalue-equivalence`).  The Phase E
     flux-shape sentinel
     (:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`)
     **no longer xfails** — it runs as a plain L1 test, the
     structurally-independent Variant-α anchor.

The twin-path bug Phase D left open
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase D's fix lived entirely in the **apply-matvec path**.  The
Phase D Carlson seed is invoked by
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` via the
``MorelMontryAngularSweep.psi_half_seed`` composition; that
covers every Krylov-driven call.  But ORPHEUS's curvilinear
production default is **source iteration**, which (pre-Phase-F) dispatched
through the then-production ``transport_sweep`` entry rather than
through the apply matvec, and the two paths ran **different code**
to seed the M-M half-angle recurrence:

.. list-table:: Apply vs SI/sweep dispatch divergence (pre-Phase-F)
   :header-rows: 1
   :widths: 24 38 38

   * - Path
     - Carlson seed site
     - Pre-Phase-F state
   * - Apply matvec (Krylov)
     - ``_mm_weighted_angular_recurrence_single_level``
       :math:`\to`
       ``CarlsonInwardSweep``
       via ``MorelMontryAngularSweep.psi_half_seed``
     - **CORRECT** — Phase D Carlson seed installed; Gate 1.1
       XPASS on the per-ordinate flat-flux residual probe
       (residual :math:`\le 10^{-15}`).
   * - SI/sweep (Source Iteration)
     - ``_sweep_1d_spherical`` (the dissolved ``sweep.py``) line 474
       (spherical) and ``_sweep_1d_cylindrical`` line 634
       (cylindrical per-:math:`\mu`-level loop)
     - **WRONG** — hardcoded ``psi_angle = np.zeros((nx, ng))``,
       the very same Phase B zero seed Phase D diagnosed and
       replaced on the apply-matvec twin.  The bug survived
       Phase D's regression suite untouched.

The cylindrical site has its own per-:math:`\mu`-level twin —
each level's azimuthal recurrence enters with the same hardcoded
zero.  Cylindrical Gate 1.1 passed empirically pre-Phase-D
because — for the **level-symmetric** quadrature exercised here — the
first-swept ordinate's seed weight is the **dead** first-ordinate weight
(:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`), which
annihilates the wrong seed at source.  (This is **not**
:math:`\alpha`-dome ':math:`\alpha=0`' level-edge cancellation
*absorbing* the seed — a level-symmetric-only reading, false for a
product quadrature; see the #280 Phase 2.5b correction at
:ref:`sn-phase-d-gate-1-1-empirical`.)  Cardinal Rule 2 (architecture)
nonetheless demands the structural fix on the sister path even when the
empirical signature is invisible there.

Phase F-Step-2 mesh-refinement evidence (sphere)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Step 2 numerics-investigator probe ran SN on
``sphere_2g_3reg`` (A|B|A reflective 2-group sphere, R=2.0 cm,
GL-8) at :math:`n_{\text{total}} \in \{40, 80, 160, 320\}` and
Variant α (composite-GL trajectory-resolvent reference) at
:math:`n_r \in \{24, 36, 48, 72, 96\}` matching effective
refinements.  The full table from
``.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md``:

.. list-table:: SN sphere pre-Phase-F mesh refinement (g=0 ratios)
   :header-rows: 1
   :widths: 10 18 18 18 18 18

   * - :math:`n_{\text{total}}`
     - :math:`k_{\text{eff}}`
     - :math:`\bigl|\text{sf}[0]/\text{sf}[1] - 1\bigr|`
       (pole)
     - :math:`\bigl|\text{sf}[N{-}1]/\text{sf}[N{-}2] - 1\bigr|`
       (outer)
     - log-log slope (pole)
     - log-log slope (outer)
   * - 40
     - 1.3578153066
     - 4.78e-01
     - 1.13e-01
     - —
     - —
   * - 80
     - 1.3576649296
     - 4.94e-01
     - 9.75e-02
     - −0.049 (**DIV**)
     - +0.21
   * - 160
     - 1.3576295736
     - 5.11e-01
     - 6.59e-02
     - −0.049 (**DIV**)
     - +0.57
   * - 320
     - 1.3576226569
     - 5.27e-01
     - 3.88e-02
     - −0.043 (**DIV**)
     - +0.76

A linear-in-:math:`h` extrapolation of the pole ratio gave
``ratio = 0.473 + 1.06·h`` — a **fixed structural asymptote
at 0.473**, not 1.  The outer ratio converged toward 1 at
:math:`\sim \mathcal{O}(h^{3/4})`, slower than the
:math:`\mathcal{O}(h^2)` DD interior, consistent with a
first-order BC-trace truncation that *vanishes* under
refinement.  Variant α at all five refinements gave inner /
outer ratios **monotonically → 1** (1.001949 → 1.000010 inner;
1.027508 → 1.001004 outer), confirming SN as the outlier and
ruling out the BC-interpretation alternative.

Per the
``vv-principles`` Step 2 decision matrix, the pole cell fires
**Branch 3 (DIVERGENT, high urgency)** and the outer cell fires
**Branch 1 (O(h^p), file follow-up)**.  This made the dispatch
to Step 3 (deep diagnostic) mandatory.

Phase F-Step-3 isolation: SI vs Krylov split
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Step 3 diagnostic ran the **same problem** through
:func:`~orpheus.sn.solver.solve_sn` with the only variable
changed being the ``inner_solver`` kwarg:

.. list-table:: SI vs Krylov on ``sphere_2g_3reg`` n=40 (pre-Phase-F)
   :header-rows: 1
   :widths: 26 18 18 19 19

   * - Inner solver
     - :math:`k_{\text{eff}}`
     - :math:`\text{sf}[0]/\text{sf}[1]`
     - :math:`\text{sf}[N{-}1]/\text{sf}[N{-}2]`
     - cv(ψ\@i=0)
   * - ``"source_iteration"`` (sweep)
     - 1.38069560
     - **0.5223**
     - 0.8871
     - **0.520**
   * - ``"krylov"`` (apply matvec)
     - 1.38464040
     - **1.0288**
     - 0.9745
     - 0.445

The Krylov path **eliminates the pole anomaly entirely** at
n=40, and Krylov's pole ratio converges to 1 cleanly under
refinement (1.029 at n=40 → 1.002 at n=80 → 1.0018 at n=160 —
:math:`\sim\mathcal{O}(h^2)` consistent with second-order DD).
Same materials, same quadrature, same mesh, same
:class:`MorelMontryAngularSweep` pole closure with the Phase D
Carlson seed installed on its ``psi_half_seed`` field — the
**only** difference is which inner-solver dispatch is used.
The Krylov path went through
:meth:`InvertibleOperator.apply` (which consumes the Phase D
Carlson seed correctly); the SI path went through the then-production
``transport_sweep`` entry (which carried the **legacy zero
seed** untouched by Phase D).

This split is the smoking gun that pins the bug to
:file:`orpheus/sn/sweep.py:474` (and :file:`orpheus/sn/sweep.py:634`
for the cylindrical per-level twin).  See
``.claude/agent-memory/numerics-investigator/phase_f_step3_diagnostic.md``
for the full empirical trail.

Source-driven Hébert (3.434)–(3.435) — the math
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The apply-matvec's Phase D Carlson seed consumes a
:math:`\psi`-current array shaped ``(ng, M, nx)`` and folds it
to :math:`\bar Q = \frac{1}{2} \Sigma_t \phi_0` via the
Legendre :math:`\ell = 0` projection (Eq. :eq:`hebert-3-432-source`
with :math:`P_0(-1) = 1`).  The SI/sweep path has **no such
current :math:`\psi` array at sweep start** — the entire point
of one SI iteration is to *produce* the updated angular flux.
What the SI loop **does** carry at sweep start is the
within-group source

.. math::
   :label: phase-f-q-1d-decomposition

   Q_{\text{1d}}(i, g) \;\equiv\;
   Q^{\text{scatt}}_{\text{within}}(i, g)
   \;+\; \frac{1}{k_{\text{eff}}}\,Q^{\text{fiss}}(i, g)
   \;+\; Q^{\text{ext}}(i, g),

i.e. the **isotropic** within-group source from the previous
power-iteration's scalar flux + fission moment + external
source.

On the fixed-point solution of the SI loop, the operator
identity :math:`L \cdot \psi = Q_{\text{1d}}` is satisfied
ordinate-by-ordinate, with :math:`L` carrying only
:math:`\Sigma_t \psi` for the current isotropic ORPHEUS scope.
The scalar-flux Legendre moment satisfies
:math:`\phi_0 = \sum_n \mathcal{W}_n \psi_n`.  Combining
gives, on the fixed point,

.. math::
   :label: phase-f-source-eq-sigt-phi0

   \Sigma_t(r) \cdot \phi_0(r) \;=\; Q_{\text{1d}}(r),

so the cell-averaged source at :math:`\mu = -1` (Eq.
:eq:`hebert-3-432-source` collapsed to :math:`L = 0`,
isotropic) admits two equivalent expressions:

.. math::
   :label: phase-f-q-bar-twin-forms

   \bar Q_i
   \;=\; \tfrac{1}{2}\,\Sigma_t(r_i) \cdot \phi_0(r_i)
   \quad\text{(apply path: builds }\phi_0\text{ from input }\psi\text{)}

   \bar Q_i
   \;=\; \tfrac{1}{2}\,Q_{\text{1d}}(r_i)
   \quad\text{(sweep path: takes }Q_{\text{1d}}\text{ directly).}

**The two are identical on the fixed point** by Eq.
:eq:`phase-f-source-eq-sigt-phi0`.  Off the fixed point they
differ by the SI residual :math:`r_k = Q_{\text{1d}} -
\Sigma_t \phi_0^{(k)}`, which vanishes as the SI loop
converges.  The sweep path's source-driven Carlson seed is
therefore the **canonically equivalent** invocation of the
same Hébert §3.9.4 math, packaged for a code path that has
:math:`Q_{\text{1d}}` available but not the per-ordinate
:math:`\psi`.

The factor :math:`\tfrac{1}{2}` is the Legendre fold weight
:math:`(2\ell + 1)/2` at :math:`\ell = 0` times
:math:`P_0(-1) = 1`.  For an :math:`L \ge 1` anisotropic
operator (not currently in scope for ORPHEUS's apply
matvec, but flagged in the
``CarlsonInwardSweep``
class docstring's L=0 WARNING block), additional terms
:math:`(2\ell + 1) Q_\ell \cdot (-1)^\ell / 2` for
:math:`\ell \ge 1` would enter — the source-driven helper
would need a moment vector ``Q_ell[ell, i, g]`` rather than
the present ``Q_bar[i, g]`` to recover the canonical
construction.

With :math:`\bar Q_i` from either formula, the inward DD
recurrence Eqs. :eq:`hebert-3-434`–:eq:`hebert-3-435`
proceed identically to the apply path:

.. math::
   :label: phase-f-carlson-seed-source-driven

   \bar\phi_i \;=\;
   \frac{\Delta r_i \cdot \tfrac{1}{2}\,Q_{\text{1d}}(r_i)
          \;+\; 2 \cdot \bar\phi_{i+1/2}}
        {\Delta r_i \cdot \Sigma_t(r_i) \;+\; 2},
   \qquad
   \bar\phi_{i-1/2} \;=\; 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}

(sequential in cells from :math:`i = nx - 1` inward to
:math:`i = 0`, vectorised across groups).  The
``bc_outer_value`` at :math:`\bar\phi_{nx+1/2}` is the
outer-face angular flux at :math:`\mu = -1`, realised through
the BC operator on the persistent outflow buffer
``bc_outer``.

Equivalence on the converged eigenmode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Foundation test
:func:`tests.sn.sweep.core.test_sweep_vs_apply_consistency`
pins the source-vs-:math:`\psi` equivalence directly: for any
flat-:math:`\psi` field ``ψ_const`` with ``bc_outer_value =
ψ_const`` (reflective) and ``Q_1d = Σ_t · Σw · ψ_const`` (the
within-group source built by SI from
``φ_0 = Σw · ψ_const``), the two helpers return
**bit-identical seeds** (up to FP non-associativity).
Apply-path:
``CarlsonInwardSweep``
``(psi_level=ψ_const·ones, ctx)`` produces ``Q̄ = 0.5 · Σ_t · Σw
· ψ_const``; sweep-path:
:func:`carlson_inward_sweep_from_source`
``(Q_bar=0.5·Q_1d, ...)`` produces the same ``Q̄`` — the
recurrence is identical, the bit-equal result is the
**single-invariant property** the test pins.

The architectural choice: one helper, two consumers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase F's structural choice is **factor the helper, delegate
from the strategy** — the Cardinal Rule 2 (architecture)
imperative.  The pre-Phase-F implementation had Eqs.
:eq:`hebert-3-434`–:eq:`hebert-3-435` open-coded inside
``CarlsonInwardSweep.__call__``.  Naive options for the
backport:

* **Option 1 (REJECTED) — duplicate the recurrence loop in
  the sweep path.**  Two copies of the inward DD recurrence,
  one driven by ``Q̄ = 0.5 · Σ_t · φ_0`` (apply path), one by
  ``Q̄ = 0.5 · Q_1d`` (sweep path).  Equivalent at the
  algorithmic level but a Cardinal-Rule-2 architecture
  violation: a future bug fix to one copy would need to be
  audited against the sister — exactly the failure mode that
  produced the Phase F bug in the first place.
* **Option 2 (REJECTED) — invoke**
  ``CarlsonInwardSweep``
  **directly from the sweep path with a synthesized**
  ``psi_level``
  **array.**  The strategy's ``__call__(psi_level,
  context)`` Protocol signature takes ``(ng, M, nx)`` — the
  sweep would have to allocate a flat-:math:`\psi` proxy of
  the right shape just to feed it through the Legendre fold
  that would extract ``φ_0`` from the proxy.  Mathematically
  equivalent but wasteful and obscures intent.
* **Option 3 (SHIPPED) — factor**
  :func:`carlson_inward_sweep_from_source`
  **as a free function that takes** ``Q̄``
  **directly, and have the strategy delegate.**

``CarlsonInwardSweep.__call__`` now reads (in essence):

.. code-block:: python

   def __call__(self, psi_level, context):
       # ψ -> φ_0 -> Q̄ fold (apply-path-specific)
       phi_0 = np.einsum("gmi,m->gi", psi_level, context.weights)
       Q_bar = 0.5 * context.sigma_t * phi_0
       # Delegate to the source-driven recurrence
       return carlson_inward_sweep_from_source(
           Q_bar=Q_bar,
           sigma_t=context.sigma_t,
           dr=context.dr,
           bc_outer_value=context.bc_outer_value,
       )

The sweep path consumes the helper directly with ``Q_bar = 0.5
· Q_1d.T``.  **Single source of truth, two structurally
equivalent invocation points.**  A future bug fix to the
recurrence (e.g., an :math:`L \ge 1` anisotropic extension)
lands in **one** place; both consumers inherit it
automatically.

Why the cylindrical site needed the fix too
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cylindrical Gate 1.1 passes empirically pre-Phase-F (for the
**level-symmetric** quadrature exercised here the first-swept
ordinate's seed weight is the **dead** first-ordinate weight
:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`, which
annihilates the wrong zero seed at source — **not** :math:`\alpha`-dome
telescoping "absorbing" it, a level-symmetric-only reading false for a
product quadrature; see the :ref:`sn-phase-d-gate-1-1-empirical`
discussion and its #280 Phase 2.5b correction for the
sphere-vs-cylinder asymmetry).  Phase F nonetheless fixes
both sites for two reasons:

#. **Cardinal Rule 2 (architecture)**: structural alignment of
   the canonical math at both sites prevents a future
   refactor from introducing an asymmetric bug that only the
   sphere catches.  The sweep-path helper is the same code
   regardless of geometry; consuming it consistently from
   both geometries is the architecturally clean choice.
#. **Defense in depth against future stress probes**: on any
   cylinder rule where the first-ordinate seed weight is **live**
   (a **product** quadrature already is — :math:`c_{\rm in}[m_0]\ne 0`,
   #280 Phase 2.5b), the wrong zero seed enters the fixed point.
   Fixing both sites now is cheap insurance.

The cylindrical fix sits inside the per-:math:`\mu`-level
loop (lines 678–714 of :file:`orpheus/sn/sweep.py`).  The
helper is invoked **once per level** with the level-specific
``bc_outer_value`` extracted from the persistent outflow
buffer at the most-inward ordinate of the level.  The
linearity of the helper in ``Q_bar`` and ``bc_outer_value``
ensures the per-level invocations remain commutative with
the outer-loop level iteration.

Phase F empirical evidence (post-fix)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The post-Phase-F state recovers the canonical SN behaviour on
the smoking-gun case:

.. list-table:: ``sphere_2g_3reg`` n=40 — pre/post Phase F
   :header-rows: 1
   :widths: 36 22 22 20

   * - Diagnostic
     - Pre-Phase-F (SI)
     - Post-Phase-F (SI)
     - Krylov (reference)
   * - :math:`\text{sf}[0]/\text{sf}[1]`
       (pole ratio, target ~1)
     - **0.522**
     - **0.778**
     - 1.029
   * - :math:`\text{sf}[N{-}1]/\text{sf}[N{-}2]`
       (outer ratio, target ~1)
     - 0.887
     - **0.997**
     - 0.974
   * - :math:`\text{cv}(\psi@i=0)`
       (Pomraning isotropy, target ~0)
     - 0.520
     - **0.404**
     - 0.445
   * - :math:`\max/\min(\psi@i=0)`
       (target ~1)
     - **6.4×**
     - **1.16×**
     - 1.18×
   * - :math:`k_{\text{eff}}`
     - 1.38069560
     - 1.38069560
     - 1.38464040

The pole ratio jumps from a structural plateau at 0.473–0.522
(divergent under refinement) up to a stable ~0.778 plateau
that holds at n=320 — the **structural divergence is
gone**.  ``sf[0]/sf[1] = 0.778`` is not yet ``1`` because the
SI fixed point still differs from the Krylov fixed point by
the residual O(h) WDD spatial-closure asymmetry (see
:ref:`sn-phase-f-residual-o-h-open` below), but the
**diverging-vs-refinement** signature that made the Phase E
flux-shape sentinel xfail-strict is closed.

.. note:: **Retraction (2026-06-13, Issue #196).**

   The table below logs a residual SI :math:`\neq` Krylov
   :math:`\mathcal{O}(h)` gap (pole ratio 0.778 vs Krylov 1.029;
   :math:`\Delta k` 0.286 % at n=40 halving per mesh doubling) and
   reads it as a benign discretisation artefact of "two methods now
   solving the same equation".  **That interpretation is wrong.**  The
   methods did NOT yet solve the same equation at Phase F: the
   *shared* closure seeds were still O(1)-wrong on non-flat fields
   (ERR-058).  After ERR-058 (#195) fixed the seeds, the
   :math:`\mathcal{O}(h)` gap **collapsed to the iteration floor** —
   SI :math:`\equiv` Krylov to :math:`|\Delta k|\approx
   1.9\mathrm{e}{-11}` and L∞ flux-shape :math:`\approx
   2.4\mathrm{e}{-10}` on ``sphere_2g_3reg`` n=40 (from a bug-era
   3.9e-3 / ~30 %); the pole ratio reaches 1 to that floor.  The
   measured numbers stay below as bug-era evidence; the
   production-decision record and post-fix evidence are
   :ref:`sn-issue-196-eigenvalue-equivalence`.

Mesh-refinement convergence (SI vs Krylov, post-Phase-F):

.. list-table:: Post-Phase-F SI-vs-Krylov convergence on ``sphere_2g_3reg`` (bug-era — gap closed by ERR-058/#196)
   :header-rows: 1
   :widths: 10 22 18 22 18 14

   * - :math:`n`
     - :math:`k_{\text{eff}}` (SI)
     - :math:`\text{sf}[0]/\text{sf}[1]` (SI)
     - :math:`k_{\text{eff}}` (Kr)
     - :math:`\text{sf}[0]/\text{sf}[1]` (Kr)
     - :math:`\Delta k`
   * - 40
     - 1.38069560
     - 0.7776
     - 1.38464040
     - 1.0288
     - 0.286 %
   * - 80
     - 1.38075258
     - 0.7771
     - 1.38261730
     - 1.0125
     - 0.135 %
   * - 160
     - 1.38078077
     - 0.7771
     - 1.38167934
     - 1.0018
     - 0.065 %

*(Bug-era reading, retracted — see the note above.)* The
:math:`k_{\text{eff}}` gap between SI and Krylov drops
by a factor of 2 per mesh doubling — apparent clean
:math:`\mathcal{O}(h)` convergence to a shared limit.
Pre-Phase-F the SI sat on the wrong structural fixed point
(0.473–0.522 ratio asymptote diverging from 1) while Krylov
converged to ~1 — the two methods **solved different
equations**, and refinement made it worse for SI.  Phase F
removed the *divergent* signature but, as ERR-058 later
showed, left a *shared* O(1)-on-non-flat seed defect: the two
paths still did not solve the same equation, and the residual
:math:`\mathcal{O}(h)` gap above is the slow trace of that
shared defect, not a discretisation artefact.  ERR-058 (#195)
fixed the seeds and the gap collapsed to the iteration floor
(:ref:`sn-issue-196-eigenvalue-equivalence`); there is no
residual O(h) gap in production.

Files touched by Phase F
~~~~~~~~~~~~~~~~~~~~~~~~

**Modified production code**

* :mod:`orpheus.sn.sweep.psi_half_angle_seed` — NEW free
  function :func:`carlson_inward_sweep_from_source` (lines
  358–419 of the module); ``CarlsonInwardSweep.__call__``
  refactored to delegate after folding ``psi_level → Q̄``;
  ``__all__`` extended.
* ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) —
  :func:`_sweep_1d_spherical` line ≈ 472–530: replaces the
  legacy ``psi_angle = np.zeros((nx, ng))`` with the Phase F
  Carlson seed call (uses ``bc_outer_obj.apply(bc_outer)`` to
  derive ``bc_outer_value`` at the most-inward ordinate, mirror
  of the apply-path's Phase D logic);
  :func:`_sweep_1d_cylindrical` lines ≈ 678–714: per-level
  Carlson seed inside the :math:`\mu`-level loop, replaces
  the inline level-zero init.

**Tests added**

* :func:`tests.sn.sweep.core.test_phase_c_gates.test_sweep_curvilinear_per_ordinate_flat_flux_residual`
  — **Gate 1.6**, the dual of Gate 1.1 for the SI/sweep path.
  Parametrised over geometry (sphere × cylinder) and
  :math:`\Sigma_t \in \{0.5, 1.5\}`.  Pins
  apply-path-vs-sweep-path bit-identity on the helper output
  AND the flat-:math:`\psi` algebraic identity at :math:`\Sigma
  w = 2` (Hébert convention).  Carries
  ``@pytest.mark.verifies("dd-curvilinear-scalar")`` and
  ``@pytest.mark.catches("ERR-026")`` — see
  :ref:`sn-phase-f-test-wiring` for the proposed extension to
  the Phase F equation labels.
* :file:`tests/sn/sweep/core/test_sweep_vs_apply_consistency.py` —
  NEW file, **57 foundation tests** pinning:

  #. Apply-path vs sweep-path Carlson seed bit-equivalence on
     matching ``Q̄`` (the load-bearing structural invariant).
  #. Linearity of
     :func:`carlson_inward_sweep_from_source` in ``Q_bar`` and
     ``bc_outer_value`` independently (Protocol-shape contract
     preservation).
  #. SI-vs-Krylov :math:`k_{\text{eff}}` agreement on
     homogeneous reflective spheres (the degenerate case
     where the Phase F fix is invariant — same eigenvalue
     pre- and post-fix).

**Updated tests**

* :func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`
  — *(Phase-F action, since superseded.)* Phase F updated the
  ``xfail-strict`` reason string from *"UNRESOLVED structural
  discrepancy with hypothesised pole issue"* to *"Phase F closed
  gross divergence; residual O(h) drift awaits further work"*, on
  the expectation a future tightening would self-enforce removal.
  **The xfail was removed by ERR-058 (#195):** the canary now runs
  as a plain L1 test (the structurally-independent Variant-α anchor;
  see :ref:`sn-issue-196-eigenvalue-equivalence`).

**Snapshot regeneration**

* 6 curvilinear regression snapshots regenerated under the
  Phase F fix:

  * ``tests/sn/regression/snapshots/sphere_2g_homogeneous_dd_n20.npz``
  * ``tests/sn/regression/snapshots/sphere_2g_3reg_dd_n40.npz``
  * ``tests/sn/regression/snapshots/sphere_2g_p1_aniso_dd_n20.npz``
  * ``tests/sn/regression/snapshots/cyl_1g_homogeneous_LS4_dd_n20.npz``
  * ``tests/sn/regression/snapshots/cyl_1g_homogeneous_product_dd_n20.npz``
  * ``tests/sn/regression/snapshots/cyl_2g_3reg_LS4_dd_n40.npz``

  Bit-identity break is principled per the
  ``vv-principles`` *"Bit-identity vs principled-equivalence"*
  framework: the new seed is the canonical Hébert value
  (replaces the diagnosed wrong zero); the
  structurally-independent verification reference is
  Variant α (composite-GL trajectory-resolvent, accessed via
  Gate 4.2); the drift is algorithmic (intended) and
  well-defined.  All 5 Gate 4.2 snapshots still PASS at the
  Phase E tightened tolerances (sphere
  :math:`r_{\text{tol}} = 2 \times 10^{-2}`, cylinder
  :math:`3 \times 10^{-2}`).

.. _sn-phase-f-residual-o-h-open:

ERR-026 manifestation #7 — CLOSED by ERR-058 (#195), verified + pinned by #196
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: Status — manifestation #7 is CLOSED (Issue #196, 2026-06-13)
   :class: important

   **This was the LAST open manifestation of ERR-026; closing it
   formally retires the curvilinear-SN wrong-fixed-point family.**  The
   "residual :math:`\mathcal{O}(h)` SI-vs-Krylov gap" reading below —
   including the SI :math:`\neq` Krylov tables above (pole ratio 0.778
   vs Krylov 1.029; :math:`\Delta k` converging :math:`\mathcal{O}(h)`)
   and Options (a)/(b)/(c) — was the *two-distinct-systems* picture and
   is **bug-era history**.  The gap was NOT a discretisation artefact;
   it was the shared closure-seed defect (ERR-058), manifest
   differently on the two then-distinct paths.

   * **ERR-048** (Phase G Step 2, 2026-05-13) closed only the **L0
     flat-field** twin-agreement: it patched the SI sweep to MATCH the
     apply-matvec conventions on the homogeneous streaming-equilibrium
     gauntlet (pole-face WDD IC mirror + Carlson seed normalisation).
     The **L1 heterogeneous eigenvalue** :math:`\mathcal{O}(h)`
     asymmetry that manifestation #7 names **PERSISTED** — which is
     exactly why #196 stayed OPEN — because the shared closure seeds
     were still *exact-on-flat / O(1)-wrong-on-non-flat* (the ERR-058
     defect).
   * **ERR-058** (Issue #195, 2026-06-12) was the TERMINAL fix: it
     replaced the shared closure seeds with correct ones — the
     coupled-pole spatial seed :math:`\psi(0,+\mu)=\psi(0,-\mu)` and the
     ``AngularEdgeExtrapolation``
     half-angle seed — so BOTH inner solvers operate on the SAME correct
     discrete operator.
   * **#196** (2026-06-13) VERIFIED the eigenvalue-path equivalence and
     added the permanent regression gate.  See
     :ref:`sn-issue-196-eigenvalue-equivalence` for the measured
     evidence (sphere :math:`|\Delta k|=4.68\mathrm{e}{-12}`, cylinder
     :math:`1.91\mathrm{e}{-11}` on the bug-era snapshot cases) and the
     gate description.

   Option (c) (keep SI, accept an O(h) gap) is moot — there is no gap;
   Option (b) (flip to Krylov) is the opposite of what landed (SI is
   restored as the faster default).  The full production-decision
   record is :ref:`sn-issue-196-eigenvalue-equivalence`.

The bug-era reading (preserved as history) ran: Phase F closed the
**structural** pole defect (the divergent ratio at the pole cell on
heterogeneous MR) and the **outer-cell** defect (sf[-1]/sf[-2]
essentially reaches 1); what was thought to remain was a milder
**convergence-rate** gap between SI and Krylov on heterogeneous MR
snapshots (at n=40 per-cell shape differing by ~5 %, apparently
converging :math:`\mathcal{O}(h)` toward zero under refinement),
logged in ``error_catalog.md`` as **ERR-026 manifestation #7**:

   *"SI-vs-Krylov per-cell agreement (residual O(h) WDD
   asymmetry) — OPEN, new follow-up after Phase F."*

That row now reads **CLOSED by ERR-058 (#195), verified + pinned by
#196**.  The Phase E flux-shape sentinel
:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`
**no longer xfails** — it runs as a plain L1 test (the
structurally-independent Variant-α anchor; see
:ref:`sn-issue-196-eigenvalue-equivalence`).  The two viable
closures that were tracked as Phase F-extensions are recorded here
**only as bug-era history** — neither was taken, because both
presupposed the shared fixed point was correct and only the
arithmetic differed, which the terminal diagnosis refuted:

* **Option (a) — Sweep WDD-closure refinement** *(bug-era,
  not taken).*  Investigate the per-cell WDD recurrence
  :math:`\psi_{n+1/2} = (\psi_n - (1-\tau)\psi_{n-1/2})/\tau`
  in ``_sweep_1d_spherical`` to identify the residual
  numerical asymmetry vs the apply matvec's symmetric closure
  :math:`\psi_{n\pm 1/2} = \tau \psi_{\text{next}} +
  (1-\tau) \psi_{\text{this}}`.  This presumed the seed was
  correct and only the spatial closure differed — false: the
  seed itself was the defect.
* **Option (b) — Flip curvilinear ``inner_solver`` default to
  Krylov** *(bug-era, not taken — in fact reverted).*
  :func:`solve_sn` for spherical / cylindrical would route
  through the Krylov inner (which carried the Phase D Carlson
  seed and produced the cleanly-converging fixed point).  This
  was the Phase-D-era flip; ERR-058 made the SI sweep correct, so
  the default was **reverted to ``source_iteration``** (SI is
  :math:`\sim 10^2\times` faster and now equivalent).

Phase F shipped **option (c)** at the time (keep SI default, achieve
structural alignment of the seed math, accept a residual O(h) gap)
on the reasoning that the methods "now solve the same equation".
The terminal diagnosis (ERR-058) showed they did **not** yet solve
the same equation — the *shared* fixed point was itself wrong on
non-flat fields.  Once the seeds were fixed there is no residual gap
to accept: SI and Krylov agree to the iteration floor at the
eigenvalue level (see :ref:`sn-issue-196-eigenvalue-equivalence`),
and bit-identically at the fixed-source level.

The anti-pattern Phase F surfaced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Twin-path fix incompleteness** is a Mode-3
(missing-factor / wrong-term-initialization) anti-pattern.
The Phase D fix was scoped to the apply-matvec path because
Gate 1.1 runs through the apply-matvec; the SI/sweep path's
zero seed was untouched.  The bug survived Phase D's entire
regression suite untouched because:

* Phase D's Gate 1.1 MMS xfail-strict marker is on the
  apply-matvec path; the SI/sweep path didn't run that probe.
* The 6 curvilinear regression snapshots were SI-generated
  under the wrong seed; the snapshots **encoded the bug
  bit-identically** and "passed" by tautology.
* Homogeneous degenerate cases (1G, flat-flux reflective) gave
  k = νΣ_f / Σ_a independent of the flux shape, masking the
  structural divergence on the eigenvalue side.
* The heterogeneous-MR case was marked ``xfail`` for **flux
  shape** (Phase E), not for **eigenvalue**, so the
  shape-sentinel signal was deliberately not enforced.

The lesson, **proposed for addition to**
``vv-principles/SKILL.md`` *§ Anti-patterns* (per the Phase F
closeout memo §"Lessons (proposed for skill catalogue)"):

   *Whenever a fix is applied to one of two structurally-
   mirrored production paths (apply-matvec vs SI/sweep,
   prepass vs postpass, etc.), MUST audit the OTHER path for
   the same defect.  Mode 3 wrong-term-initialization
   defects often appear in pairs; fixing one path without
   auditing its sister is a Cardinal Rule 2 (architecture)
   violation that ERR-026 instantiated twice.*

.. _sn-phase-f-test-wiring:

Test wiring proposal — Phase F equation labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase F declares three new equation labels:
:eq:`phase-f-q-1d-decomposition`,
:eq:`phase-f-source-eq-sigt-phi0`,
:eq:`phase-f-q-bar-twin-forms`, and
:eq:`phase-f-carlson-seed-source-driven`.  The
:eq:`phase-f-carlson-seed-source-driven` label is the
canonical Hébert (3.434)–(3.435) recurrence in
source-driven form — semantically the **same recurrence**
as :eq:`hebert-3-434` and :eq:`hebert-3-435` but with the
sweep-path source substitution made explicit.

The new Gate 1.6 test
:func:`tests.sn.sweep.core.test_phase_c_gates.test_sweep_curvilinear_per_ordinate_flat_flux_residual`
already carries
``@pytest.mark.verifies("dd-curvilinear-scalar")`` and
``@pytest.mark.catches("ERR-026")``.  Per the project's
V&V harness wiring, the test SHOULD additionally declare:

* ``@pytest.mark.verifies("phase-f-carlson-seed-source-driven")``
  — pins the source-driven recurrence on the bit-identity
  helper-vs-strategy probe.
* ``@pytest.mark.verifies("phase-f-q-bar-twin-forms")``
  — pins the apply-vs-sweep source-equivalence identity (the
  load-bearing structural invariant of Phase F).

The other two labels
(:eq:`phase-f-q-1d-decomposition` and
:eq:`phase-f-source-eq-sigt-phi0`) document the
*decomposition* of the SI source and the *fixed-point
identity* :math:`\Sigma_t \phi_0 = Q_{\text{1d}}`; both are
verified transitively by the existing SI convergence
infrastructure (the SI inner-tolerance is the gate that
enforces the fixed-point identity to machine precision).
The proposed wiring is tracked as a follow-up to the V&V
audit harness (see Issue #194 for the sister case of
``hebert-3-43X`` labels — same pattern, same fix).

Pointers
~~~~~~~~

* **Phase F plan**:
  ``.claude/plans/issue_168_phase_f_curvilinear_boundary_eigenvector.md``
  — context, hypothesis, three-step structure, sub-agent
  dispatch chain.
* **Step 2 numerics memo**:
  ``.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md``
  — mesh-refinement convergence study, SN-vs-Variant-α
  outlier identification, Step 2 branch-3 decision.
* **Step 3 diagnostic memo**:
  ``.claude/agent-memory/numerics-investigator/phase_f_step3_diagnostic.md``
  — fix-site identification (the smoking gun), SI-vs-Krylov
  isolation, Option-A-vs-B implementation analysis.
* **Phase F closeout memo**:
  ``.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md``
  — what shipped, the empirical evidence tables, files
  touched, residual-open items.
* **ERR-026 catalogue narrative**:
  ``.claude/skills/vv-principles/error_catalog.md`` (§ ERR-026,
  *"What Wave H Phase F added"*) — manifestation table update
  #6 CLOSED, #7 (new) OPEN.
* **Sister section on the BC apply call sequence**:
  :ref:`bc-phase-f-three-bc-applies-per-sweep-iteration` in
  :doc:`/theory/foundations/boundary_conditions` — extends the Phase D
  two-BC-applies-per-matvec narrative to cover the SI sweep's
  Phase F invocation.


.. _sn-err-058-closure-seed-closeout:

ERR-058 — the curvilinear closure-seed fix (Issue #195 CLOSED)
--------------------------------------------------------------

.. admonition:: Status banner
   :class: important

   **Issue #195 CLOSED 2026-06-12.**  ERR-058 closes the curvilinear
   *wrong-fixed-point* family — the open loop the Phase A–F narrative
   tracked under the name "ERR-026 PARTIAL CLOSURE".  Two
   independent closure SEEDS in the curvilinear within-group operator
   were wrong on every non-flat field; both are now replaced.  In
   production:

   * The **half-angle thread seed** is
     ``AngularEdgeExtrapolation``
     (the new
     :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
     ``psi_half_seed`` default).  It replaces
     ``CarlsonInwardSweep``,
     whose proxy source :math:`\bar Q = \Sigma_t\phi_0/\!\sum w` was the
     dominant defect.
   * The **spatial pole-face seed** of the outward (:math:`\mu>0`)
     sweep is the *mirror inward sweep's pole-face outflow* — the
     Carlson coupled-pole continuity :math:`\psi(0,+\mu)=\psi(0,-\mu)`
     — replacing the historical innermost-cell-centre read
     :math:`\psi(\Delta r/2)`.
   * The curvilinear inner default returned from ``"krylov"`` to
     ``"source_iteration"`` (both
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` and the
     eigenvalue :func:`~orpheus.sn.solver.solve_sn`): post-unification
     the sweep and matvec are ONE discrete system, so SI
     :math:`\equiv` Krylov **bit-identical for fixed-source** and to
     the **iteration floor for eigenvalue** (the eigenvalue solve wraps
     the inner in power iteration — see
     :ref:`sn-issue-196-bit-identical-vs-floor`).  SI is
     :math:`\sim 10^2\times` faster than GMRES (no restart, ERR-053).

   ``CarlsonInwardSweep`` is **retained** (not deleted) as the
   registered host of the canonical Hébert §3.9.4 recurrence
   (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`),
   reachable only by explicit opt-in.  Its class docstring carries a
   ``.. warning::`` block recording the proxy-source caveat by design,
   so a future session cannot re-activate it as a default unaware of
   the falsification.

   .. note:: **Retraction (2026-07-04, Issue #282 route (a)).**  The
      ``CarlsonInwardSweep`` *strategy class* was NOT ultimately
      retained — route (a) deleted the whole ``PsiHalfAngleSeed``
      family.  What survives is the pure Hébert recurrence
      :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
      (a free function, now the SOLVE engine driven by the **true** q½
      source rather than the falsified proxy, no opt-in), plus the
      inlined
      :meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`
      for non-carrying cylinder levels.  See
      :ref:`sn-282-seed-strategy-zoo`.

   The **anisotropic** curvilinear MMS gates improved :math:`\sim 50\times`
   and are now limited by a *fixed-quadrature angular floor* of the
   half-angle thread interpolation — a test-design retune tracked at
   `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_,
   **not** a residual instance of this wrong-fixed-point class.

Motivation preserved — what the Phase A–F loop was chasing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Phase A–F sections — A, D, and F above; B and C with the
production machinery in :doc:`curvilinear_one_group` — are preserved
verbatim as the *investigation history*.  Their reasoning was sound at the time and is
pedagogically load-bearing — a future reader asking "why did anyone
try a Carlson inward sweep, an apply-vs-sweep twin audit, a
Krylov-default flip?" must find the answer there.  This subsection only
flips the tenses on the *terminal* claims those sections reached:

* Phase D **was expected to** close ERR-026 once the apply-matvec
  half-angle seed was made non-zero; it closed the per-ordinate
  *flat-flux identity* but left the assembled operator wrong on
  non-flat profiles (the Carlson proxy source).  The "PARTIAL CLOSURE /
  pre-asymptotic transient" framing it shipped **was** the best
  available reading of the evidence then; it is now superseded.
* Phase F **was expected to** close the SI-vs-Krylov gap by backporting
  the same Carlson seed to the sweep.  It did make the two paths share
  the *seed strategy*, but both still drove it with the wrong proxy
  source, so the residual "manifestation #7 O(h) gap" it logged was a
  symptom of the shared defect, not a discretisation artefact.  After
  the Depth-B/Wave-T matvec unification, the sweep and matvec became
  ONE discrete system; the gap vanished by construction once the seed
  was fixed.

The premise the *issue itself* carried — a benign "pre-asymptotic
transient" that finer meshes would clear — was **empirically refuted**:
on ``main`` the isotropic curvilinear MMS error PLATEAUS
mesh-independently (sphere :math:`\sim 0.0413`, cylinder
:math:`\sim 0.0494`, :math:`n_x` 20 :math:`\to` 640, orders
:math:`\to 0`), with SI :math:`\equiv` Krylov bit-identical.  No
refinement helps a plateau.

The two manifestations — one class, both flat-field-exact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERR-058 is a **closure-seed inconsistency**: two discrete seeds inside
the curvilinear within-group operator were constructed so that they are
*exact* on spatially/angularly flat fields and *O(1)-wrong* on every
other field.  Because a discrete closure seed is part of the operator,
each seed had to be verified per-ordinate on a NON-flat field — and was
not.

.. _sn-err-058-manifestation-a:

Manifestation (a) — the spatial pole-face seed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The outward (:math:`\mu>0`) radial sweep needs an inflow value at the
pole face :math:`r=0`.  The historical matvec
(now the fused
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`) and
the sweep both read the innermost CELL-CENTRE flux :math:`\psi(\Delta r/2)`
as if it were the pole-FACE value — a half-cell offset.  On a flat
radial profile :math:`\psi(\Delta r/2)=\psi(0)`, so the read is exact;
on the manufactured :math:`A(r)=\sin(\pi r/R)` ansatz it is
:math:`\mathcal{O}(h)`-wrong.  The DD face chain propagates that seed
error as an *undamped odd–even alternation*, and the area-weighted
streaming amplifies it by :math:`\sim A/V \sim 1/r` near the pole.

**The fix — Carlson coupled-pole continuity.**  At :math:`r=0` the
outward characteristic is the *continuation* of the inward one: a
neutron travelling inward along :math:`-\mu` that reaches the centre
emerges travelling outward along :math:`+\mu`, so

.. math::
   :label: sn-err-058-coupled-pole-continuity

   \psi(0,\,+\mu) \;=\; \psi(0,\,-\mu).

.. (vv-status rationale) Representational identity: the r=0
   pole-continuity boundary condition coupling the mirror ordinates.
   Not a solver claim (no eigenvalue / flux value). The verifiable
   content is the per-ordinate operator-admission gate
   (test_curvilinear_operator_admits_mms, catches ERR-058) + the
   strategy-owned seed-adjoint bit-identity (test_g_adjoint_reciprocity).
.. vv-status: sn-err-058-coupled-pole-continuity documented

The :math:`-1` (inward) sweep is therefore run FIRST.  Its pole-face
outflow, read at the *mirror* ordinate
``quad.reflection_index("x")``, seeds the :math:`+1` (outward) sweep.
This is **data** — propagated from the outer boundary, lower-triangular
in cell-visit order — not a self-reference.  It is the
"inward-determines-outward" pole condition deferred at Phase C
(`Issue #192 <https://github.com/deOliveira-R/ORPHEUS/issues/192>`_),
now landed.  The seed is exact on flat :math:`\psi` (so every flat-flux
gate is untouched), lower-triangular (so the operator stays
forward-substitutable), and the matvec and sweep capture/consume it
identically (so the pair stays ONE discrete system).

The **adjoint** routes the :math:`+1` seed cotangent into the
:math:`-1` reversal's initial outflow cotangent at the mirrored
ordinates (see
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action_transpose`,
pinned by the dense-probe transpose oracle
``derivations/diagnostics/diag_p42_adjoint_oracle.py``).

.. _sn-coupled-pole-mu-level-invariant:

The μ-level-preservation invariant the mirror seed relies on
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The coupled-pole continuity :eq:`sn-err-058-coupled-pole-continuity`
seeds the outward (:math:`+\mu`) pole face from the inward (:math:`-\mu`)
sweep's pole outflow read **at the mirror ordinate** — concretely,
``pole_face_seed = outflow_at_inner.T[quad.reflection_index("x")]`` in
the fused matvec
(:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk` and
its adjoint partner) and ``psi_in = pole_outflow[mirror[global_n]]`` in
the SI sweep twin (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.sweep`).
That single index — ``reflection_index("x")[n]`` — is what makes the
seed correct, and it carries a load-bearing assumption that is invisible
in the code but essential to the physics.

**The invariant.**  For the mirror seed to realise
:math:`\psi(0,+\mu_r)=\psi(0,-\mu_r)`, the partner
``reflection_index("x")[n]`` MUST be the *intra-level sign-flip partner*
of ordinate :math:`n`: the ordinate in the **same** :math:`\mu`-level
(same axial cosine :math:`\mu_z` — the level index) with the radial
cosine :math:`\mu_x` negated and :math:`\mu_y,\mu_z` held.

.. math::
   :label: sn-coupled-pole-mu-level-invariant-eq

   m \;=\; \mathrm{reflection\_index}(\text{"x"})[n]
   \;\Longrightarrow\;
   \mu_x[m] = -\,\mu_x[n],\quad
   \mu_y[m] = \mu_y[n],\quad
   \mu_z[m] = \mu_z[n].

.. (vv-status rationale) Structural / representational identity: the
   defining property the x-mirror partner MUST satisfy for the
   coupled-pole seed to be physically correct (intra-level μ_x
   sign-flip). Not a solver claim (no eigenvalue / flux value). The
   verifiable content is the foundation gate
   test_x_reflection_is_intra_level_signflip_partner (asserts the three
   equalities over gauss_legendre/level_symmetric/product cubatures) +
   its involution sibling; documented-only.
.. vv-status: sn-coupled-pole-mu-level-invariant-eq documented

**Why the physics demands it.**  The pole continuity is a statement at a
*fixed axial direction*: a neutron travelling inward along :math:`-\mu_x`
at axial cosine :math:`\mu_z` reaches the centre and emerges travelling
outward along :math:`+\mu_x` at the **same** :math:`\mu_z`.  The axial
component does not turn at the pole — only the radial one reverses.  So
the reflected partner must stay in the same :math:`\mu`-level; a
cross-level partner would couple two different axial directions and seed
the outward sweep with a value from the *wrong* characteristic.

**Why it holds by construction today.**  Two facts conspire.  First,
``reflection_index("x")`` resolves to
``reflection_partners[0] = _find_reflections(-mu_x, mu_y, mu_z, ...)``
(:func:`orpheus.numerics.quadrature.directional._compute_sphere_reflection_partners`),
which nearest-neighbour-matches each node against its image with
**only** :math:`\mu_x` negated — :math:`\mu_y,\mu_z` are passed through
unchanged.  Second, the cylinder/sphere level is grouped on the
**axial** cosine: the level factories key ``level_indices`` on
:math:`|\mu_z|` (sphere / level-symmetric — ``rules_sphere.py``) or hold
:math:`\mu_z=\mu_{\rm GL}` fixed per level (product — ``rules_product.py``),
never on :math:`\mu_x`.  Because the x-mirror holds :math:`\mu_z` and the
level is indexed by :math:`\mu_z`, the x-mirror provably maps an ordinate
to a partner in its own level.  The two facts are *independent* code
sites, so the invariant is an emergent property of their agreement — not
something either site enforces alone.

.. warning::

   **This is a silent-corruption surface — a Mode-7 blind spot at the
   operator-internals level.**  If ``reflection_index("x")`` ever
   returned a **cross-level** partner — a future cubature whose
   reflection table is built differently, or a refactor of
   :func:`~orpheus.numerics.quadrature.directional._compute_sphere_reflection_partners`
   that no longer holds :math:`\mu_z` — then
   ``pole_outflow[mirror[n]]`` would read a *different axial direction's*
   pole value, and the break would be **completely silent under the
   existing solver suite**.  The reason is the same blindness that hid
   ERR-058 itself: on a spatially/angularly **flat** :math:`\psi` field
   the mirror partner's pole value equals the ordinate's own value, so
   the seed is exact regardless of *which* ordinate it reads.  Every
   flat-flux gate, every streaming-equilibrium L0, every reflective
   :math:`k_\infty` check would stay green while the operator quietly
   coupled the wrong characteristics on any non-flat field (``vv-principles``
   Mode 7 — the ansatz-simplification blindness — operating on the
   operator's own internals, exactly the ERR-058 class).  A scalar /
   particle-balance residual cannot see it either, because the
   :math:`\alpha`-dome telescoping (above) sums away per-ordinate seed
   errors.

**Why it now has its own foundation gate.**  Because the solver tests are
structurally blind to a cross-level regression, the invariant is pinned
*directly* — not through any flux or eigenvalue, but as a property of the
quadrature's reflection table itself — by the foundation test
:func:`tests.sn.sweep.curvilinear.test_coupled_pole_mu_level_invariant.test_x_reflection_is_intra_level_signflip_partner`.
It asserts all three equalities of
:eq:`sn-coupled-pole-mu-level-invariant-eq` (intra-level membership,
:math:`\mu_x` sign-flip, :math:`\mu_y,\mu_z` held) over the
``gauss_legendre`` / ``level_symmetric`` / ``product`` cubatures the
curvilinear sweep actually uses; the sibling
:func:`~tests.sn.sweep.curvilinear.test_coupled_pole_mu_level_invariant.test_x_reflection_is_an_involution`
pins ``mirror ∘ mirror = id`` (the partner relation is symmetric — a
necessary corollary of the sign-flip).  Both are
``@pytest.mark.foundation``; the first carries
``@pytest.mark.verifies("sn-err-058-coupled-pole-continuity")``, tying
the table-level invariant to the continuity equation it underwrites.
This gate is the regression tripwire that turns the silent-corruption
surface into a loud one: a cross-level reflection table fails it
immediately, *before* any solver ever runs.

.. note::

   **Re-scope of Issue #193.**  This invariant is what
   `Issue #193 <https://github.com/deOliveira-R/ORPHEUS/issues/193>`_
   now pins.  The issue *originally* targeted a different "level-locality"
   concern — that the cylindrical matvec's
   ``bc_outer.apply``-once-then-per-level-extract pattern stayed correct.
   That concern **dissolved**: Wave O O.4a.2 removed the ``bc_outer.apply``
   keystone from the matvec entirely (the reflective coupling :math:`B`
   moved *outside* the bare sweep as a first-class sibling — see the
   boundary-condition extraction record at :ref:`bc-extraction`), and the
   surviving SI-sweep seed reads the **raw** inflow trace with no
   ``apply`` at all, so the
   apply/restrict commutativity the original test would have guarded is
   now vacuous.  The genuinely load-bearing :math:`\mu`-level invariant
   *moved* to the coupled-pole seed mirror documented here, and that is
   what #193 was re-scoped to gate.

.. _sn-err-058-manifestation-b:

Manifestation (b) — the angular half-angle thread seed (the dominant defect)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Morel–Montry angular recurrence (Hébert §3.9.4 Eqs. 3.437/3.439)
threads the half-angle face fluxes
:math:`\psi_{m\pm 1/2,i}` across a :math:`\mu`-level and needs a starting
seed :math:`\psi_{1/2,i}` at the level's most-inward angular edge.  The
Phase D / Phase F ``CarlsonInwardSweep`` solved the canonical
*sweep-side* starting-direction ODE (Hébert Eqs. 3.432–3.435) for that
seed, but drove it with the **proxy source**

.. math::
   :label: sn-err-058-proxy-source

   \bar Q_i \;=\; \frac{\Sigma_{t,i}\,\phi_{0,i}}{\sum_n w_n},
   \qquad \phi_{0,i} = \sum_n w_n\,\psi_{n,i},

.. (vv-status rationale) Literature-transcribed definition of the
   falsified proxy source (the CarlsonInwardSweep half-angle seed).
   Recorded as the diagnosed defect, not a solver claim. Its falsity
   is what the per-ordinate operator-admission gate (catches ERR-058)
   detects; documented-only.
.. vv-status: sn-err-058-proxy-source documented

which equals the true within-group source ONLY at the flat-flux
equilibrium :math:`\Sigma_t\phi_0 = \bar Q`.  On any non-equilibrium
field — every MMS reference, every vacuum or heterogeneous problem —
the seed solves the *wrong* starting-direction ODE.  The measured
consequence on the isotropic MMS input
(:math:`\psi_n = A(r)/W`, scalar value :math:`0.5`): the seed returns
:math:`\bar\phi = 0.5777` where the correct angle-flat value is
:math:`0.5000`, and the per-ordinate redistribution residual reaches
:math:`\pm 55` at the pole, :math:`\pm 13` in the bulk, against a
continuous streaming of :math:`\pm 0.31`.  **This was the dominant
defect.**

**The fix —**
``AngularEdgeExtrapolation``.
For the *operator* (matvec) to be consistent, the seed must approximate
the *input field's* own value at the level edge :math:`\mu_{\rm start}`
— a pure angular-extrapolation problem with NO dependence on
:math:`\Sigma_t`, the source, or the boundary trace.  The new strategy
extrapolates linearly in :math:`\mu` through the level's two most-inward
distinct-:math:`\mu` ordinates :math:`(m_0, m_1)`:

.. math::
   :label: sn-err-058-edge-extrapolation

   \psi_{1/2,i} \;=\; (1-t)\,\psi_{m_0,i} + t\,\psi_{m_1,i},
   \qquad
   t \;=\; \frac{\mu_{\rm start} - \mu_{m_0}}{\mu_{m_1} - \mu_{m_0}}.

.. (vv-status rationale) Representational identity: the
   operator-consistent half-angle thread seed (AngularEdgeExtrapolation,
   the new psi_half_seed default) as a fixed linear map. Not a solver
   claim. The verifiable content is the per-ordinate operator-admission
   gate (catches ERR-058), the isotropic MMS L1 ladders, and the
   strategy-owned seed-adjoint bit-identity; documented-only.
.. vv-status: sn-err-058-edge-extrapolation documented

The starting-direction edge :math:`\mu_{\rm start}` (sphere
:math:`-1`; cylinder :math:`-\sqrt{1-\xi_p^2}`, the level's most-inward
azimuthal edge) is single-sourced from the SAME
:math:`\alpha`/:math:`\tau` construction site as the
:math:`\alpha`-dome (``orpheus.geometry.reduced_operator``) and threaded
to the strategy via the REQUIRED
``CarlsonSweepContext.mu_start``
field — no default, so a forgotten cylinder site cannot silently fall
back to the sphere value.

**Exactness ladder.**  The extrapolation is

* **exact on angle-flat fields**, because the barycentric weights sum
  to one: :math:`(1-t)+t=1`.  Every per-ordinate flat-flux identity
  gate is therefore untouched.
* **exact on linear-in-:math:`\mu` fields**: write the level's input as
  :math:`\psi_{m,i}=a_i+b_i\,\mu_m`.  Then

  .. math::

     \psi_{1/2,i}
       &= (1-t)(a_i + b_i\mu_{m_0}) + t(a_i + b_i\mu_{m_1}) \\
       &= a_i + b_i\bigl[(1-t)\mu_{m_0} + t\mu_{m_1}\bigr]
        = a_i + b_i\,\mu_{\rm start},

  the last bracket collapsing to :math:`\mu_{\rm start}` by the
  definition of :math:`t` in :eq:`sn-err-058-edge-extrapolation`.  The
  M-M recurrence is itself a Möbius/affine map in :math:`\mu`; seeded
  with :math:`\psi_{1/2}=a+b\,\mu_{1/2}` it threads the ENTIRE
  half-angle grid exactly as :math:`\psi_{m+1/2}=a+b\,\mu_{m+1/2}` (for
  *unclamped* :math:`\tau`).  Hence the P1-class anisotropic MMS
  references — whose ansatz is exactly :math:`(A(r)+B(r)\mu)/W` — are
  *admitted* by the operator.
* **O(\Delta\mu^2)-consistent** on general smooth angular profiles —
  the same order as the angular discretisation itself.
* **linear in the input**, so the operator-algebra operations
  (:meth:`apply`, :meth:`apply_transpose`, dense probing) are
  preserved.  The strategy OWNS its adjoint
  (``PsiHalfAngleSeedBase.seed_adjoint``),
  a fixed linear scatter of the seed cotangent onto the two stencil
  ordinates, so a strategy swap on
  :class:`MorelMontryAngularSweep` swaps both the forward and reverse
  maps at once.

.. note::

   The :math:`\tau`-clamp
   (:math:`\tau \to \max(0.5,\min(1.0,\tau_{\rm raw}))`,
   Bailey–Morel–Chang) breaks the *exact* linear-in-:math:`\mu`
   threading wherever it is active.  This is part of the residual
   anisotropic angular floor quantified at :ref:`Issue #229
   <sn-err-058-aniso-floor>` below — NOT a wrong-fixed-point defect.

Why every gate stayed green — the blindness analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both seeds hid behind a regime in which they are exact, and the V&V
suite sat *entirely inside* that regime.  This is
``vv-principles`` Mode 7 (MMS / ansatz simplification bias) operating
not on a manufactured solution but on the *operator's own internals*.

.. list-table:: Which fields each closure seed is exact on (the blind regime)
   :header-rows: 1
   :widths: 30 34 36

   * - Closure seed
     - Exact on
     - Gate that sat in the blind regime
   * - Spatial pole-face (a)
     - flat radial :math:`\psi`
       (:math:`\psi(\Delta r/2)=\psi(0)`)
     - streaming-equilibrium L0; reflective-equilibrium k\ :sub:`∞`
   * - Angular thread (b)
     - flat-flux equilibrium
       (:math:`\Sigma_t\phi_0 = \bar Q`)
     - per-ordinate flat-flux identity (Gate 1.1); homogeneous
       reflective

**The :math:`\alpha`-dome telescoping made scalar checks blind to
(b).**  The M-M redistribution coefficients form a dome that telescopes
under the angular weight sum: :math:`\sum_n w_n\,(\alpha_{n+1/2} -
\alpha_{n-1/2}) = 0` REGARDLESS of the half-angle thread values.  Any
weight-summed (scalar-flux / particle-balance) residual therefore cannot
see a wrong half-angle thread — ``vv-principles`` anti-pattern #8
("NEVER accept particle balance as L0 evidence; require per-ordinate
residual") instantiated *inside a diagnostic*.  During the #195
investigation this telescoping made the scalar residual go
:math:`\mathcal{O}(h^2)` after fixing only (a), while the per-ordinate
residual was still :math:`\mathcal{O}(10)` — which mis-supported a
"near-singular operator / two-solutions gauge mode" hypothesis until a
dense SVD showed :math:`\sigma_{\min}\approx 0.9` (never near-singular)
and the *per-ordinate* check named the real defect.

**Historical compensation explains the Phase-D-era O(h²) reading.**  At
Phase D time the apply path measured :math:`\mathcal{O}(h^2)` under
Krylov (the premise this issue inherited), because its closure internals
compensated differently from the sweep.  The Depth-B/Wave-T matvec
rebuild changed the redistribution assembly and surfaced the latent seed
inconsistency; the SWEEP, by contrast, had ALWAYS plateaued (#195's own
SI data :math:`[0.083, 0.095, 0.098]`).  The wrong fixed point was the
sweep's all along — the same class as #98's original 35 %-at-:math:`r=0`
finding.

The three refuted intermediate hypotheses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Recording the dead paths so a future session does not re-run them
(Sphinx is the brain):

#. **"Pre-asymptotic transient"** (the issue's own premise).  Refuted by
   the :math:`n_x` 20 :math:`\to` 640 plateau — orders :math:`\to 0`, no
   refinement helps.
#. **A pure :math:`r=0`-regularity extrapolation spatial seed**
   (:math:`1.5\,\psi_0 - 0.5\,\psi_1`).  Implemented; it drove the
   *scalar* residual to :math:`\mathcal{O}(h^2)` but the solution still
   plateaued, because the dominant defect was the *angular* seed (b),
   invisible to the scalar residual by the telescoping above.  Superseded
   by the coupled-pole seed :eq:`sn-err-058-coupled-pole-continuity` for
   (a) — which is *data* rather than a one-sided stencil — once (b) was
   diagnosed.
#. **A "near-null gauge mode" theory** (apparent two-solutions paradox).
   Falsified by a dense SVD: :math:`\sigma_{\min}\approx 0.9`, the
   operator is well-conditioned.  The paradox was an artefact of the
   scalar-blind diagnostic, not a property of the operator.

Production closure decision — post-fix evidence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Post-fix (measured 2026-06-12), the isotropic curvilinear MMS solution
error collapses into a clean second-order ladder, with SI and Krylov
bit-identical:

.. list-table:: Post-fix isotropic curvilinear MMS L2 ladders (SI ≡ Krylov)
   :header-rows: 1
   :widths: 16 16 16 16 16 16

   * - :math:`n_x`
     - 20
     - 40
     - 80
     - 160
     - 320
   * - sphere :math:`\|\phi_h-A\|_{L^2}`
     - 1.49e-2
     - 3.73e-3
     - 9.28e-4
     - 2.31e-4
     - 5.74e-5
   * - sphere order
     -
     - 2.00
     - 2.01
     - 2.01
     - 2.01
   * - cylinder :math:`\|\phi_h-A\|_{L^2}`
     - 2.16e-3
     - 5.39e-4
     - 1.35e-4
     - 3.37e-5
     -
   * - cylinder order
     -
     - 2.00
     - 2.00
     - 2.00
     -

The magnitude band :math:`10^{-8} < {\rm err} < 10^{-3}` is satisfied
(sphere :math:`n_x \ge 80`, cylinder :math:`n_x \ge 40`).  SI converges
:math:`\sim 10^2\times` faster than GMRES here (sphere :math:`n_x=160`:
:math:`\sim 0.11\,\mathrm{s}` SI vs :math:`\sim 69\,\mathrm{s}` Krylov),
which is why the curvilinear default returned to source iteration.

The decisive *structural* gate is the **per-ordinate, volume-weighted**
operator-admission residual of :math:`\psi_{\rm ref}` (the scalar
residual is blind, per the telescoping above):

.. list-table:: Per-ordinate volume-weighted residual of ψ_ref under (L+C−S−B), post-fix
   :header-rows: 1
   :widths: 25 25 25 25

   * - Geometry
     - :math:`n_c=40`
     - :math:`n_c=80`
     - measured order
   * - sphere
     - 1.97e-3
     - 9.7e-4
     - :math:`\approx 1.5` (pole-adjacent bounded band under
       the :math:`r^2\,dr` weight)
   * - cylinder
     - 5.50e-5
     - 1.37e-5
     - :math:`\approx 2.0` (pointwise :math:`\mathcal{O}(h^2)`
       everywhere)

The sphere's sub-quadratic residual order is benign: the
pole-adjacent cells legitimately carry a bounded non-decaying
*pointwise* residual where the closure truncation meets the
:math:`\Delta A/V \sim 1/h` geometry factor on cells whose volume
vanishes as :math:`r^2\,dr`; the solution-error ladder above proves it
harmless.  **Bug-era** values for this gate were :math:`\mathcal{O}(10^{-1})`-class
(per-ordinate pointwise up to :math:`\pm 55` at the pole) — three-plus
orders of magnitude above the post-fix bounds, which is the margin the
ERR-058 catcher asserts.

The quadrature/truncation floor is the radial DD closure order itself;
the post-fix sphere/cylinder ladders sit at the DD design order
(2.00–2.01), so "have you tried finer quadrature?" is pre-empted — the
solution-error *is* second-order, and the only residual non-quadratic
behaviour is the volume-weighted pole band, which the solution error
does not inherit.

.. _sn-err-058-aniso-floor:

The anisotropic angular floor — deferred to Issue #229
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Resolved (2026-06-13, W1–W5).**

   This deferral is now fully treated at
   :ref:`sn-curvilinear-aniso-norm-reconciliation`.  The W1–W5
   root-cause program found the single "floor" sketched below was
   **three distinct errors** (a sphere pole-cell spatial closure, a
   sphere angular :math:`\tau`-clamp floor, and a cylinder angular
   floor), separated by a norm difference (volume-weighted L2 vs
   pointwise / L∞).  In particular, the "Open research paths" item
   below — "Unclamped-:math:`\tau` threading on a linear-in-:math:`\mu`
   shell" — was **executed** (W1): the sphere clamp was removed (it was
   mis-cited and 100 % spurious; see
   :ref:`sn-tau-clamp-vindication`), which cleaned the coarse rate but,
   surprisingly, *raised* the S16 fine floor (the prior lower floor was
   a fortuitous cancellation, not a gain).  The numbers preserved below
   are correct historical evidence; the comprehensive treatment and the
   per-error production decisions are at the reconciliation section.

The **anisotropic** curvilinear MMS (ansatz :math:`(A(r)+B(r)\mu)/W`)
dropped :math:`\sim 50\times` under the fix and is now limited by a
*fixed-quadrature angular floor*, NOT by a residual wrong-fixed-point
defect.  The mechanism: the aniso MMS imposes :math:`\psi_n` per
ordinate, so there is no angular error *at the imposed ordinates* — but
the M-M redistribution consumes half-angle THREAD values
:math:`\psi_{m\pm 1/2}` that the recurrence *interpolates*.  On an
angle-varying ansatz the thread's interpolation error is an
angular-quadrature-resolution effect: under spatial refinement at fixed
quadrature the solution converges to an angular floor, and the
pure-spatial rate + magnitude assertions cannot both hold once the
spatial error drops below it.  The floor *scales with quadrature order*
in both geometries — the structural signature confirming the
angular-thread attribution:

.. list-table:: Anisotropic angular floor vs quadrature order (post-ERR-058, SI inner)
   :header-rows: 1
   :widths: 22 24 54

   * - Case
     - Quadrature
     - Floor behaviour
   * - sphere aniso
     - S16 (shipped)
     - :math:`n_x` 10→160: ``[5.9e-2, 1.5e-2, 4.0e-3, 1.15e-3, 7.3e-4]``;
       floor :math:`\approx 7\mathrm{e}{-4}`
   * - sphere aniso
     - S32
     - err@80 = 9.5e-4, err@160 = 2.9e-4 (floor drops :math:`\sim 2.5\times`)
   * - cylinder aniso
     - :math:`n_\mu{=}4` (shipped)
     - :math:`n_x` 40→80: ``1.91e-2 → 1.90e-2`` (hard floor 1.9e-2)
   * - cylinder aniso
     - :math:`n_\mu{=}8`
     - :math:`n_x` 40→80: ``7.50e-3 → 7.39e-3``
       (floor drops :math:`\sim 2.6\times` per :math:`n_\mu` doubling)

The :math:`\tau`-clamp (above) contributes a constant to this floor by
breaking the exact linear-in-:math:`\mu` threading where active.  The
quadrature-aware retune (raise the case quadrature, or split the claim
into a pre-floor spatial-O(h²) segment + a separate
angular-convergence assertion) is `Issue #229
<https://github.com/deOliveira-R/ORPHEUS/issues/229>`_.

Infrastructure retained
~~~~~~~~~~~~~~~~~~~~~~~~~

Per the aggressive-retirement *exception* (a correct primitive that
would be needed if the obstruction is ever bypassed is kept as an
oracle), ERR-058 deletes no correct machinery:

.. list-table:: Curvilinear closure-seed primitives — status after ERR-058
   :header-rows: 1
   :widths: 38 18 44

   * - Primitive
     - Status
     - Why kept
   * - ``AngularEdgeExtrapolation``
     - **production**
     - The new ``psi_half_seed`` default; operator-consistent on
       non-flat fields.
   * - ``CarlsonInwardSweep``
       + :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
     - retained, opt-in
     - Correct *source-driven* Hébert §3.9.4 recurrence; would seed a
       future TRUE-source sweep-side closure.  Proxy-source caveat
       pinned in its docstring ``warning`` block.
   * - ``ZeroSeed``
     - retained, ablation
     - Reproduces the Phase B behaviour for A/B regression-safety
       comparison.
   * - Coupled-pole spatial seed in
       :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
       / :meth:`~orpheus.sn.loss_representation._OneDimScanWalk.sweep`
     - **production**
     - The :math:`\psi(0,+\mu)=\psi(0,-\mu)` continuity; matvec + sweep
       share it (one discrete system).

.. note:: **Retraction (2026-07-04, Issue #282 route (a)).**  The three
   ``PsiHalfAngleSeed`` strategy rows above (``AngularEdgeExtrapolation``
   / ``CarlsonInwardSweep`` / ``ZeroSeed``) are superseded: route (a)
   **deleted** the whole strategy family — including the
   ``AngularEdgeExtrapolation`` "production default", which was itself
   the #282 walk-order back edge.  What is genuinely retained is the
   free-function Hébert engine
   :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
   (now the SOLVE driver, on the **true** q½ source) and the inlined
   :meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`
   (non-carrying cylinder levels).  The coupled-pole spatial seed row is
   unaffected.  See :ref:`sn-282-seed-strategy-zoo`.

Open research paths
~~~~~~~~~~~~~~~~~~~~~

Two paths could lift the anisotropic angular floor without changing the
isotropic O(h²) result:

#. **TRUE-source-driven sweep-side seed** — **LANDED as Issue #282 route
   (a) (2026-07-04).**  This path predicted the resolution exactly:
   replace the ``AngularEdgeExtrapolation`` *input-field* extrapolation
   with the canonical Hébert recurrence
   :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
   driven by the genuine within-group source
   :math:`\bar Q_i = \sum_\ell \tfrac{2\ell+1}{2}Q_\ell(r_i)(-1)^\ell`
   (the full Legendre fold, not the :math:`\Sigma_t\phi_0` proxy), making
   the *starting-direction transport* exact rather than the *input-field
   value* exact.  Route (a) shipped precisely this — and the "likely
   diagnostic probe" proposed here (holding spatial mesh fixed and
   sweeping quadrature order) is exactly the **angular-order N-sweep**
   that certified the re-pose principled.  The one refinement over the
   prediction: the seed also had to become **first-class typed state**
   (not just a better strategy) to kill the walk-order back edge that
   made the *solve* non-direct.  See
   :ref:`sn-282-direct-starting-direction-solve`.
#. **Unclamped-:math:`\tau` threading on a linear-in-:math:`\mu` shell.**
   The exact linear-:math:`\mu` threading (above) holds only for
   unclamped :math:`\tau`; quantify the clamp's contribution to the
   floor and, where the cell is well-resolved
   (:math:`\tau_{\rm raw}\in[0.5,1.0]`), thread unclamped to recover the
   exact P1 admission.  Likely probe: the
   :ref:`Issue #229 <sn-err-058-aniso-floor>` floor table with the clamp
   disabled on resolved cells.

Session trail (V&V audit trail)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **ERR-058 catalogue narrative**:
  ``.claude/skills/vv-principles/error_catalog.md`` (§ ERR-058) — the
  authoritative two-manifestation mechanism + post-fix evidence.
* **Re-scope record**: `Issue #195
  <https://github.com/deOliveira-R/ORPHEUS/issues/195>`_ comments
  (2026-06-12) — the premise refutation and the decisive probe-3
  residual evidence.
* **Diagnostics**: ``derivations/diagnostics/diag_195_probe{1,2,3}_*.py``
  (the plateau / error-profile / operator-admission probes), promoted to
  the gate ``tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py``.
* **Investigator memo**:
  ``.claude/agent-memory/numerics-investigator/issue_195_root_cause_2_pole_closure.md``.

Verification chain
~~~~~~~~~~~~~~~~~~~

The ERR-058 fix is pinned by, in order of structural decisiveness:

#. :func:`tests.sn.verification.mms.test_curvilinear_operator_admits_mms.test_operator_admits_isotropic_mms_per_ordinate`
   (``@pytest.mark.l1`` + ``catches("ERR-058")``) — the fast
   per-ordinate volume-weighted operator-admission gate (the structurally
   decisive check, immune to the telescoping blindness).
#. :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_spherical_mms_converges_second_order`
   and
   :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_cylindrical_mms_converges_second_order`
   (``catches("ERR-058")``) — the end-to-end L1 ladders whose
   ``xfail`` markers came off with this fix; they ``verifies`` the
   :eq:`sn-mms-spherical-psi` / :eq:`sn-mms-spherical-qext` /
   :eq:`sn-mms-cylindrical-psi` / :eq:`sn-mms-cylindrical-qext` labels.
#. The flat-flux and streaming-equilibrium gates pin the flat-field
   exactness BOTH fixes preserve (so they did not regress).
#. :func:`tests.sn.operators.test_g_adjoint_reciprocity` — pins the
   strategy-owned seed adjoints.

.. note::

   **vv-status (eq-labels added by this section).**  The labels
   :eq:`sn-err-058-coupled-pole-continuity`,
   :eq:`sn-err-058-proxy-source`, and
   :eq:`sn-err-058-edge-extrapolation` are *structural / representational*
   identities (the pole-continuity boundary condition; the falsified
   proxy-source definition; the operator-consistent edge-extrapolation
   map).  They are NOT solver claims (no eigenvalue / flux-value claim).
   Per the vv-status discipline they are ``documented`` — the
   verifiable content is the per-ordinate operator-admission gate
   (``catches("ERR-058")``) plus the strategy-owned adjoint
   bit-identity, named in the verification chain above.

.. _sn-issue-196-eigenvalue-equivalence:

Issue #196 — eigenvalue-path SI≡Krylov verification and the permanent gate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: Status — manifestation #7 verified and pinned (Issue #196 CLOSED, 2026-06-13)
   :class: important

   ERR-058 (#195, above) replaced the wrong shared closure seeds; **#196
   is the verification and regression-gate step** that confirms the
   replacement closed ERR-026 manifestation #7 at the *eigenvalue*
   level and locks the closure against re-introduction.  This was the
   LAST open manifestation of ERR-026 — with #196 closed, the
   curvilinear-SN wrong-fixed-point family is **formally retired**.

The two-layer history (why the L0 close did not suffice)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Manifestation #7 names a specific defect: a residual :math:`\mathcal{O}(h)`
SI-vs-Krylov WDD spatial-closure asymmetry on curvilinear SN.  Pre-fix,
the source-iteration inner (drives the curvilinear sweep) and the Krylov
inner (drives the apply-matvec) produced eigenvalues differing at
:math:`\mathcal{O}(h)`: **0.286 %** on ``sphere_2g_3reg`` at n=40, **~30
%** per-cell on eigenvector shape, the gap halving under mesh refinement.

Two closures touched this defect, and the honest history distinguishes
them:

* **ERR-048** (Phase G Step 2, 2026-05-13) closed only the **L0
  flat-field** twin-agreement.  It patched the SI sweep
  (``_sweep_1d_spherical`` / ``_sweep_1d_cylindrical``) to MATCH the
  apply-matvec conventions on the homogeneous streaming-equilibrium
  gauntlet (pole-face WDD IC mirror + Carlson seed :math:`\bar Q`
  normalisation).  The **L1 heterogeneous eigenvalue**
  :math:`\mathcal{O}(h)` asymmetry that manifestation #7 names
  **PERSISTED** — which is exactly why #196 stayed OPEN — because the
  shared closure seeds were still *exact-on-flat /
  O(1)-wrong-on-non-flat* (the ERR-058 defect).
* **ERR-058** (Issue #195, 2026-06-12) was the TERMINAL fix.  It
  replaced the shared closure seeds with correct ones (the coupled-pole
  spatial seed :math:`\psi(0,+\mu)=\psi(0,-\mu)` and the
  ``AngularEdgeExtrapolation``
  half-angle seed), making BOTH inner solvers operate on the SAME
  correct discrete operator.

.. _sn-issue-196-bit-identical-vs-floor:

Bit-identical (fixed-source) vs floor-equivalent (eigenvalue)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The post-fix SI :math:`\equiv` Krylov agreement is **not the same kind
of agreement** on the two solver entry points, and the distinction is
load-bearing:

* **Fixed-source** (:func:`~orpheus.sn.solver.solve_sn_fixed_source` on
  the curvilinear MMS ladders): SI :math:`\equiv` Krylov is
  **BIT-IDENTICAL**.  Post-unification the sweep and the matvec are ONE
  discrete operator on one quadrature; the within-group inner
  (``A.solve`` vs Krylov-on-``apply``) realises the *same* :math:`A^{-1}`
  arithmetic, so the two paths return the same bits.
* **Eigenvalue** (:func:`~orpheus.sn.solver.solve_sn` with
  ``inner_solver="source_iteration"`` vs ``"krylov"``): SI
  :math:`\equiv` Krylov to the **ITERATION FLOOR** (~:math:`1.9\mathrm{e}{-11}`
  in :math:`k_{\text{eff}}`, ~:math:`2.4\mathrm{e}{-10}` in flux shape),
  **NOT bit-identical**.  The eigenvalue solve wraps the inner in power
  iteration; SI and Krylov are *different iteration schemes* that
  converge to the **same correct fixed point** only to ~``inner_tol``.
  Same physics, not the same arithmetic.

Confusing the two would mis-state the verification claim.  The earlier
close-out's "SI :math:`\equiv` Krylov bit-identical on the curvilinear
MMS ladders" is correct **for the fixed-source ladders specifically**;
the eigenvalue verification below is *floor-equivalence*, which is the
right and sufficient claim for an eigenvalue solve.

Measured eigenvalue-path equivalence (Issue #196)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All values measured 2026-06-12 under tight iteration tolerances for
BOTH inner solvers (``keff_tol=1e-12``, ``flux_tol=1e-10``,
``inner_tol=1e-10``).  The eigenvalue snapshot cases are the exact
acceptance cases of #196:

.. list-table:: SI≡Krylov eigenvalue equivalence — the bug-era heterogeneous snapshot cases
   :header-rows: 1
   :widths: 30 22 22 26

   * - Case
     - :math:`|\Delta k|` (post-fix)
     - max :math:`|\Delta\varphi|_{L\infty}` (post-fix)
     - Bug-era (pre-ERR-058)
   * - ``sphere_2g_3reg_dd_n40``
     - :math:`4.68\mathrm{e}{-12}`
     - :math:`5.88\mathrm{e}{-11}`
     - :math:`|\Delta k|=3.9\mathrm{e}{-3}` (0.286 %), shape ~30 %
   * - ``cyl_2g_3reg_LS4_dd_n40``
     - :math:`1.91\mathrm{e}{-11}`
     - :math:`4.32\mathrm{e}{-11}`
     - same :math:`\mathcal{O}(h)` family, gap halving under refinement

The homogeneous (k_inf-degenerate, flat-flux) curvilinear snapshots
agree at the rounding floor — as expected, since on a flat eigenmode the
redistribution terms null and SI/Krylov differ only by accumulated
round-off:

.. list-table:: SI≡Krylov eigenvalue equivalence — homogeneous (flat-flux) snapshots
   :header-rows: 1
   :widths: 38 26 26

   * - Case
     - :math:`|\Delta k|`
     - relative :math:`\varphi` diff
   * - ``sphere_2g_homogeneous_dd_n20``
     - :math:`6.92\mathrm{e}{-13}`
     - :math:`2.15\mathrm{e}{-10}`
   * - ``cyl_1g_homogeneous_LS4_dd_n20``
     - :math:`2.22\mathrm{e}{-16}`
     - :math:`2.27\mathrm{e}{-11}`
   * - ``cyl_1g_homogeneous_product_dd_n20``
     - :math:`6.66\mathrm{e}{-16}`
     - :math:`1.10\mathrm{e}{-10}`

.. note::

   The homogeneous cases agree to the floor but, on their own, supply
   **no** evidence for the curvilinear closure — a flat eigenmode is
   degenerate (``flat = flat``; 1-group :math:`k=\nu\Sigma_f/\Sigma_a`
   is flux-shape independent, vv-principles anti-patterns #3/#4).  The
   load-bearing evidence is the **heterogeneous 2-group** cases above,
   where the flux is genuinely non-flat and the angular-redistribution
   terms are exercised.

The permanent regression gate (Issue #196)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The manifestation-#7 catcher is
:func:`tests.sn.eigenvalue.test_keff_curvilinear.test_si_krylov_eigenvalue_equivalence_sphere`
and
:func:`tests.sn.eigenvalue.test_keff_curvilinear.test_si_krylov_eigenvalue_equivalence_cylinder`,
each carrying ``@pytest.mark.catches("ERR-026")``.  Configuration:
heterogeneous 2-group fuel|moderator (region A inner, region B outer,
n=10+10), solved twice under ``inner_solver="source_iteration"`` and
``inner_solver="krylov"`` at the tight tolerances above.  The gate
asserts:

.. list-table:: Manifestation-#7 gate thresholds vs measured vs bug-era
   :header-rows: 1
   :widths: 26 24 24 26

   * - Quantity
     - Asserted bound
     - Measured (post-fix)
     - Bug-era (would trip)
   * - sphere :math:`|\Delta k|`
     - :math:`< 1\mathrm{e}{-7}`
     - :math:`1.9\mathrm{e}{-11}`
     - :math:`3.9\mathrm{e}{-3}`
   * - sphere per-group :math:`|\Delta\varphi|_{L\infty}`
     - :math:`< 1\mathrm{e}{-6}`
     - :math:`2.4\mathrm{e}{-10}`
     - ~30 %
   * - cylinder :math:`|\Delta k|`
     - :math:`< 1\mathrm{e}{-7}`
     - :math:`1.1\mathrm{e}{-11}`
     - same family
   * - cylinder per-group :math:`|\Delta\varphi|_{L\infty}`
     - :math:`< 1\mathrm{e}{-6}`
     - :math:`2.6\mathrm{e}{-11}`
     - ~30 %

A **non-flat-flux guard** (group-0 radial ``max/min > 1.2``) precedes
the equivalence assertion so the test cannot pass vacuously on a
degenerate flat mode — sphere radial ``max/min`` = 3.34, cylinder =
1.67, both well above the guard.  The bug-era values (3.9e-3 / ~30 %)
exceed the asserted bounds by **4–5 orders of magnitude**, so the gate
would have tripped loudly on the pre-fix code.

.. note::

   **Runtime-mode discipline (vv-principles anti-pattern #8).**  The
   canonical ORPHEUS invocation is ``python -O``, under which bare
   ``assert`` statements are stripped to no-ops.  These gates assert via
   bare ``assert`` inside the *collected test module*, which pytest
   rewrites at collection time so the asserts fire under ``-O``.  This
   was confirmed empirically: a negative control with a
   :math:`1\mathrm{e}{-15}` tolerance failed as required under ``-O``.

Structural-independence — why SI≡Krylov alone is not the whole proof
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SI :math:`\equiv` Krylov is **twin-path agreement** — necessary but not
sufficient (vv-principles L11: two implementations agreeing is
cross-implementation evidence, not correctness evidence).  Both inner
solvers could in principle converge to the same *wrong* fixed point.
The independent ground that makes the closure a *correctness* claim, not
merely a *consistency* claim, comes from two structurally-independent
legs:

* The **k_inf homogeneous legs** — on a uniform reflective infinite
  medium :math:`k_\infty=\nu\Sigma_f/\Sigma_a` is an analytical
  (closed-form) eigenvalue the SN snapshots must reproduce.
* The **Variant-α Green's-function cross-check**
  (:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`),
  now a plain L1 test (xfail removed), which compares the SN flux-shape
  snapshot against the composite-GL trajectory-resolvent reference
  within 8 % (sphere) / 12 % (cylinder).  This reference is a
  semi-analytical pillar structurally independent of the SN sweep, so
  agreement pins the *converged-to value*, not just twin-path
  consistency.

Production-decision record — curvilinear default reverted to SI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The curvilinear :attr:`~orpheus.sn.solver.SNSolver.inner_solver` default
is now ``"source_iteration"``, **reverted from the Phase-D Krylov
flip**.  The Phase-D flip existed ONLY because the sweep's fixed point
was wrong; ERR-058 made it correct, so SI is restored as the default —
it is :math:`\sim 10^2\times` faster than GMRES (no restart) and now
equivalent (bit-identical fixed-source / floor-equivalent eigenvalue).

Crucially, **neither of the old Phase-F closures was taken**:

* *Option (a) — make SI bit-identical to Krylov by refining the WDD
  closure* presupposed the seed was correct and only the spatial
  arithmetic differed.
* *Option (b) — flip the default to Krylov* presupposed the Krylov fixed
  point was the correct one and SI's was a discretisation artefact.

Both presupposed the *shared* fixed point was correct and only the
arithmetic differed.  The terminal diagnosis (ERR-058) showed the shared
fixed point **itself** was wrong on non-flat fields; once the seeds were
fixed, both inner solvers reach the same *correct* fixed point and the
"choose between them" framing dissolves — SI is restored on speed alone.


.. _sn-curvilinear-aniso-norm-reconciliation:

The curvilinear anisotropic-MMS "floor", reconciled (W1–W5)
===========================================================

.. admonition:: Key Facts — curvilinear anisotropic SN
   :class: important

   - **The single "#229 floor" was three distinct errors**, separated
     by a *norm difference*: the production gates measure a
     **volume-weighted L2** :math:`\sqrt{\sum_i V_i\,\Delta_i^2}`; the
     root-cause probes measured **pointwise / L∞**.  An error
     concentrated at the :math:`r \to 0` pole cell is loud in L∞ and
     near-silent in L2.
   - **(a) Sphere central-cell** :math:`\mathcal{O}(h)` **spatial
     closure (#233)** — L∞-only; :math:`\sim 75\,\%` an MMS
     midpoint-vs-shell-volume-average comparison artifact +
     :math:`\sim 25\,\%` genuine but **inherent** first order.  WONTFIX
     for diamond difference (Hébert §3.9.4 / Stacey §9.9 use plain
     diamond at the central cell).  See
     :ref:`sn-pole-cell-spatial-closure`.
   - **(b) Sphere angular** :math:`\tau`-**clamp floor** — fixed by W1
     (the clamp was mis-cited and 100 % spurious on physical fields).
     The sphere now uses the raw Bailey-Morel-Chang 2010 Eq. 43 weight.
     See :ref:`sn-tau-clamp-vindication`.
   - **(c) Cylinder angular floor (#229)** — the half-angle-thread
     INTERPOLATION floor; scales with the **azimuthal** quadrature
     :math:`n_\varphi`, structurally blocked (needs a 2-D
     :math:`(\eta, \varphi)` closure).  The sphere has a pre-floor
     :math:`\mathcal{O}(h^2)` window (clean at S32); the cylinder does
     **not** at any practical quadrature.  See
     :ref:`sn-cylinder-angular-floor`.
   - **Two unrelated "anisotropic" paths** (Issue #9): Path-(I) =
     geometric angular redistribution :math:`(1-\mu^2)/r\,\partial_\mu`
     (the M-M :math:`\alpha`-dome, P0-only — what #229 concerns);
     Path-(II) = :math:`P_1{+}` Legendre SCATTERING
     :math:`R\,\Lambda\,M` (geometry-agnostic).  See
     :ref:`sn-p1-scattering-curvilinear`.
   - **Norm gotcha**: a convergence-rate gate on a volume-weighted L2
     norm CANNOT see a pole-cell defect (the pole sits at one cell of
     :math:`V \sim h^3` → :math:`\sqrt{V} \sim h^{1.5}` →
     :math:`\sim h^{2.5}` contribution → subdominant).  A pointwise /
     L∞ probe is required to surface it.

This section closes the curvilinear-anisotropic-SN investigation
program (W1–W5, branch ``fix/curvilinear-aniso-pole-and-clamp``,
2026-06-13).  It is the sequel to the ERR-058 / #195 / #196 curvilinear
*isotropic* closure-seed family above; that family fixed the
wrong-fixed-point class (now formally retired), and what remained was
the *anisotropic* floor — which this program resolved into three
distinct, separately-actionable errors.

The headline — one floor was three errors, settled by a norm difference
-----------------------------------------------------------------------

The ERR-058 close-out (above) deferred a residual "anisotropic angular
floor" to Issue #229, citing a single
:ref:`floor table <sn-err-058-aniso-floor>`.  The W1–W5 root-cause
program found that the apparent single floor was **three structurally
distinct errors**, and the reason they had been conflated is a **norm
difference** in how two independent investigations measured the same
solves:

* The verification gates (test-architect) measure the **volume-weighted
  L2** norm :math:`\|\Delta\|_{2,V} = \sqrt{\sum_i V_i\,\Delta_i^2}` —
  the natural norm for a finite-volume scheme whose unknown is a
  cell-volume average.
* The diagnostic probes (numerics-investigator) measured **pointwise /
  L∞** — :math:`\max_i |\Delta_i|`.

The two norms weight the :math:`r \to 0` pole cell completely
differently.  Under the spherical volume weight :math:`V \sim h^3`, a
fixed pointwise error at the single pole cell contributes
:math:`\sqrt{V} \sim h^{1.5}` to the L2 sum, so an L∞-:math:`\mathcal{O}(h)`
pole error appears as :math:`\sim h^{2.5}` in L2 — **subdominant to the
interior** :math:`\mathcal{O}(h^2)`, hence invisible.  This is exactly
why the production L2 gate stayed green throughout while a pointwise
probe found a first-order pole cell.

.. list-table:: The three errors behind the "#229 floor"
   :header-rows: 1
   :widths: 6 38 22 16 18

   * - #
     - Error
     - Dominant norm
     - Quadrature-scaling?
     - Status
   * - (a)
     - Sphere pole-cell spatial closure :math:`\mathcal{O}(h)` at
       :math:`r \to 0`
     - L∞ / pointwise central flux (diluted in L2 by :math:`V \propto
       r^2`; invisible in :math:`k_{\rm eff}`)
     - no (spatial)
     - **#233 — documented inherent limitation (ERR-059, WONTFIX-for-DD)**
   * - (b)
     - Sphere angular :math:`\tau`-clamp floor (:math:`\sim 7\mathrm{e}{-4}`
       @ S16)
     - volume-weighted L2 at fine mesh
     - yes (angular)
     - **fixed (W1 unclamp)**
   * - (c)
     - Cylinder angular floor
     - both
     - yes (azimuthal :math:`n_\varphi`)
     - **structurally blocked (#229)**

The remainder of this section treats each error in turn, then the two
unrelated anisotropic paths (Issue #9).

.. _sn-tau-clamp-vindication:

(b) The :math:`\tau`-clamp vindication (W1)
-------------------------------------------

The spherical Morel--Montry weighted-diamond weight is

.. math::
   :label: sn-tau-mm-raw

   \tau_n^{\rm raw}
       \;=\; \frac{\mu_n - \mu_{n-1/2}}{\mu_{n+1/2} - \mu_{n-1/2}}
       \;\in\; [0, 1],

the **unique** weight exact for an angular flux linear in :math:`\mu`
([BaileyMorelChang2010]_ Eq. 43; the same object as
:eq:`mm-weights`).  The production code had wrapped it in a
:math:`[\tfrac12, 1]` clamp,
:math:`\tau_n = \mathrm{clip}(\tau_n^{\rm raw}, \tfrac12, 1)`, cited
to Lewis & Miller §4.5.

W1 established, by three independent lines of evidence, that the clamp
is **mis-cited and 100 % spurious on physical fields**:

#. **Literature.** [BaileyMorelChang2010]_ state the admissible range
   is :math:`\tau \in [0, 1]` and recommend *exactly* the unclamped
   :math:`\tau^{\rm raw}` (their Eq. 43) as the unique exact-on-linear
   weight; Hébert §3.9.4 uses pure diamond (:math:`\tau = \tfrac12`),
   no clamp.  Lewis & Miller §4.5 does **not** prescribe the
   :math:`[\tfrac12, 1]` clamp — the citation was wrong.
#. **Positivity is never needed.** On every realistic converged solve
   (smooth MMS, homogeneous eigenvalue :math:`k_{\rm eff} = 1`, thick
   absorber) there are ZERO negative half-angle fluxes, clamped or
   unclamped — every clamp activation is spurious (measured: 160 / 320
   / 80 / 240 activations across stress configs, 0 protective).  The
   half-flux negativity that *does* transiently appear in early SI
   iterates is inherited from a negative *input* :math:`\psi` and the
   clamp barely reduces it.  On Gauss--Legendre quadrature
   :math:`\tau^{\rm raw} \in [0.39, 0.61]` (never 0), so the unclamped
   weight is always interior to :math:`[0, 1]`.
#. **Stability without it.** Unclamped sphere source iteration
   converges with strictly positive, finite scalar flux on every
   stress config (thick absorber, near-vacuum, :math:`c = 0.999`, S64);
   the clamp costs a few SI iterations on low-scattering problems but is
   dispensable for stability.

.. admonition:: The architectural reason the static removal is correct
   :class: note

   A *dynamic* negative-flux fixup (where :math:`\tau` depends on the
   iterate :math:`\psi`) would make the streaming operator
   **nonlinear**, breaking the linear-Krylov matvec and the SI ≡ Krylov
   twin identity (Pattern-2 discipline, Cardinal Rule 2).  Because the
   fixup is *never needed* on physical fields, the principled W1 fix is
   to **drop the clamp** (a config-time, static change) and use the
   linear unclamped :math:`\tau^{\rm raw}`.  The weight :math:`\tau` is
   single-sourced in the pole-angular closure
   (:func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`,
   since Issue #236 Step C — see :ref:`sn-tau-c-on-cellvisit-live`) and
   inherited by every consumer (the SI sweep and the Krylov matvec
   both), so both twins stay linear and stay identical.

**Geometry split.**  W1 removed the clamp for the **sphere only**.  The
**cylinder keeps** it: product / level-symmetric quadratures place the
most-inward azimuthal ordinate exactly on :math:`\eta = -\sin\theta`,
giving :math:`\tau^{\rm raw} = 0` **exactly** (bit-exact, not "near
zero") — an unclamped recurrence divides by zero there.  This is a
genuine *structural* singularity the sphere provably lacks; the
cylinder's real fix is a 2-D :math:`(\eta, \varphi)` closure
(:ref:`sn-cylinder-angular-floor`), not unclamping.  See
:eq:`morel-montry-clamp` in :doc:`/theory/foundations/structured_geometry` for the
equation-of-record carrying both branches.

**Mixed accuracy signature (the gotcha).**  Unclamping does NOT
uniformly improve the anisotropic solve.  It *cleans the coarse
convergence rate* (sphere S16 coarse orders 1.978 → 1.995) but *raises
the S16 fine-mesh floor* (:math:`7.3\mathrm{e}{-4} \to 1.2\mathrm{e}{-3}`):

.. list-table:: W1 sphere aniso MMS, matched-quadrature S16 (volume-weighted L2)
   :header-rows: 1
   :widths: 20 40 40

   * - :math:`n_x`
     - Clamped (pre-W1)
     - Unclamped (post-W1)
   * - 10 → 40
     - coarse orders 1.979 / 1.978
     - coarse orders 1.995 / 1.999
   * - 80
     - 1.16e-3
     - 1.40e-3
   * - 160 (floor)
     - **7.3e-4**
     - **1.2e-3**

The lower *clamped* floor was a **fortuitous cancellation**, not a
genuine accuracy gain — the clamp's constant bias happened to partly
offset the angular-thread interpolation floor at S16.  Removing it
exposes the true #229 floor (next subsection), which is what the
unclamped weight should converge to.  Iso solves are unchanged in real
arithmetic (the clamp is silent on flat-in-:math:`\mu` fields) but
**not bit-identical** at IEEE-754: the closure
:math:`(\overline\psi - (1-\tau)\psi_{\rm in})/\tau` returns
:math:`\psi` exactly only :math:`\sim 81\,\%` of the time and within
1 ULP otherwise (reduction-order non-associativity), so the converged
homogeneous-reflective sphere drifts :math:`|\Delta k| = 2.3\mathrm{e}{-13}`,
:math:`\max|\Delta\phi| = 4.4\mathrm{e}{-13}` — an FP-tail, anchored to
the closed-form :math:`k_\infty = 1.875`.  One snapshot
(``sphere_2g_3reg_dd_n40``, genuinely non-flat) was regenerated
(:math:`k\;1.380766 \to 1.381001`); the two flat snapshots drift only
in the FP tail and were not regenerated.

**W1 gates.** ``tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py``
(closure-unit :math:`\tau`-independence on flat fields; converged
homogeneous-reflective iso anchored to :math:`k_\infty`; unclamped
positivity) + the W1 ``@slow`` aniso gates appended to
``tests/sn/verification/mms/test_curvilinear_aniso_convergence.py``
(the S32 clean-:math:`\mathcal{O}(h^2)` full-ladder claim; the S16
coarse-rate-cleaner-unclamped discriminator; the floor-scales-with-
quadrature pin).  Landed in commit ``b2d8a6d``.

.. _sn-cylinder-angular-floor:

(c) The cylinder angular floor (#229) — structurally blocked
------------------------------------------------------------

The anisotropic ansatz :math:`\psi_{\rm chosen} = (A(r) + B(r)\,\mu)/W`
imposes :math:`\psi_n` per ordinate, so there is **no angular error at
the imposed ordinates**.  But the M-M redistribution consumes
half-angle THREAD values :math:`\psi_{m\pm 1/2}` that the recurrence
**interpolates** (they are not imposed).  On an angle-varying ansatz the
thread's interpolation error is an angular-quadrature-resolution effect:
under spatial refinement at fixed quadrature the solution converges to
an angular floor, and a pure-spatial-rate assertion cannot hold once
the spatial error drops below it.

**The cylinder floor scales with the AZIMUTHAL quadrature**
:math:`n_\varphi`, **not the polar** :math:`n_\mu`.  This is the
load-bearing physical fact (and a correction to an earlier mislabel):
the radial direction cosine is :math:`\eta = \sin\theta\,\cos\varphi`,
so the M-M thread marches in azimuth :math:`\varphi` *per polar
:math:`\mu`-level*.  Measured at :math:`n_x = 80`:

.. list-table:: Cylinder aniso floor vs azimuthal quadrature (:math:`n_x = 80`, volume-weighted L2)
   :header-rows: 1
   :widths: 25 25 50

   * - :math:`n_\varphi`
     - L2 error
     - Behaviour
   * - 8
     - 1.90e-2
     - hard floor
   * - 16
     - 7.37e-3
     - drops :math:`2.58\times`
   * - 32
     - 3.10e-3
     - drops :math:`2.38\times`

while :math:`n_\mu` (polar) refinement at fixed :math:`n_\varphi`
leaves the floor **flat**.

**Why it is structurally blocked.**  Product and level-symmetric
quadratures carry **duplicate azimuthal** :math:`\eta`: ordinates come
in :math:`\pm\varphi` symmetry pairs with the same :math:`|\eta|` but
opposite :math:`\xi` (e.g. :math:`\varphi = \pi/4` and
:math:`\varphi = 7\pi/4` both give
:math:`\eta = \sin\theta/\sqrt 2`).  The M-M thread marches in
:math:`\eta` alone, so a field whose true variation is in the full
:math:`(\eta, \varphi)` plane is **not threadable exactly** by a 1-D
:math:`\eta`-march — a structural mismatch, not a tuning problem.  No
partition (midpoint / cumulative-weight / ordinate-interior) gives
:math:`\tau^{\rm raw} \in [\tfrac12, 1]` with bounded edges; the
cumulative-weight partition is exact on level-symmetric but needs
:math:`\tau^{\rm raw} \in [-4.5, 5.5]` (edges outside the level).
Closing the cylinder floor requires a genuine 2-D
:math:`(\eta, \varphi)` angular closure — **out of scope**, tracked by
`Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_.

**The sphere–cylinder asymmetry.**  The sphere DOES have a pre-floor
:math:`\mathcal{O}(h^2)` window: at S16 the coarse orders clear 1.99 and
the floor (:math:`\approx 2.9\mathrm{e}{-4}` at S32, :math:`n_x = 160`)
sits below the segment's finest spatial error, so the clean
second-order window extends to :math:`n_x = 80` at S32.  The cylinder
has **no** such window at *any* practical quadrature — even
:math:`n_\mu = 16` (:math:`N = 512`) reaches only order 1.80 on the
coarsest :math:`\{5, 10, 20\}` segment before the angular floor
dominates.  The mathematics, not runtime, is the blocker.

**W3 gate retune (Issue #229).**  Per the vv-principles anti-pattern
"a claim that cannot hold MUST NOT be asserted; pin what IS true
instead", W3 removed all five aniso xfail markers and migrated the six
equation labels to green tests:

* **Sphere** ``test_sn_spherical_aniso_mms_converges_second_order`` →
  coarse-segment ``orders[:2] > 1.9`` + magnitude band
  :math:`1\mathrm{e}{-8} < \mathrm{err}[-1] < 5\mathrm{e}{-3}`
  (loosened from :math:`1\mathrm{e}{-3}` because the W1 unclamp removed
  the fortuitous-cancellation lower floor) + ``catches("ERR-026")``.
* **Cylinder** ``test_sn_cylindrical_aniso_mms_converges_second_order``
  → floor band :math:`1\mathrm{e}{-3} < \mathrm{err}[-1] <
  5\mathrm{e}{-2}`, **no rate claim** (the floor dominates).  The
  cylinder phase-C spatial test was **repurposed** into
  ``test_cyl_aniso_floor_scales_with_quadrature``
  (:math:`\mathrm{err}(n_\varphi{=}16) < \mathrm{err}(n_\varphi{=}8)/2`
  — the verified-floor second claim that pins the angular attribution).
* The sphere prescribed-inflow redistribution test dropped its
  strict-xfail and rate gate (band :math:`1\mathrm{e}{-8} <
  \mathrm{err} < 5\mathrm{e}{-3}` + a kept converged-value
  ``assert_allclose(2e-2)``).

Landed in commit ``679a1e6`` (audit exit 0; the
:eq:`sn-mms-spherical-psi` / ``-qext`` / :eq:`sn-mms-cylindrical-psi` /
``-qext`` labels and the two spatial-convergence labels are now all
green-tested).

.. _sn-pole-cell-spatial-closure:

(a) The sphere/cylinder pole-cell :math:`\mathcal{O}(h)` spatial closure (#233, ERR-059)
----------------------------------------------------------------------------------------

This is the **new** discovery of the program and the one *surviving*
manifestation in the curvilinear-SN family.  The curvilinear scalar
flux is **first-order** :math:`\mathcal{O}(h)` at the :math:`r \to 0`
central cell in the pointwise / L∞ norm — distinct from #168 (outer
face, CLOSED), ERR-058 (the closure seed, CLOSED), and the angular
floors above.  It decomposes into three parts, **none** of which
warrants a code fix.

Decomposition
~~~~~~~~~~~~~

**Part 1 — :math:`\sim 75\,\%` MMS comparison artifact (not a solver
bug at all).**  The production spherical MMS evaluates the source at the
cell MIDPOINT ``mesh.centers`` and compares
:math:`\phi_{\rm solver}` against :math:`\phi_{\rm exact}(\text{midpoint})`.
But the spherical DD discrete unknown **IS** the cell-volume average

.. math::
   :label: sn-pole-cell-shell-average

   \overline{\phi}_{n,i}
       \;=\; \frac{4\pi}{V_i}\int_{r_{i-1/2}}^{r_{i+1/2}} r^2\,\phi_n(r)\,dr

([Hebert2009]_ Eq. 3.430 — the unknown is *defined* as the shell
average, not a point value; the diamond relation Eq. 3.431 relates it to
the face fluxes).  Under :math:`r^2\,dr` weighting the volume-average and
the midpoint point-value differ by :math:`\mathcal{O}(h)` at the pole
cell, because :math:`r_{\rm lo} = 0` maximally skews the weight (the
volume-centroid sits at :math:`\tfrac34 h`, not :math:`\tfrac12 h`).
Using the *shell-averaged* source AND comparing against the
*shell-volume-average* drops the pole error :math:`\sim 4\times`
(:math:`0.0212 \to 0.00497`) — confirming the bulk of the apparent
error is a comparison subtlety, not solver truncation.

**Part 2 — :math:`\sim 25\,\%` genuine but LITERATURE-ACCEPTED INHERENT
first order.**  Even the fully consistent finite-volume MMS (shell-avg
source + shell-avg reference) leaves the pole at clean
:math:`\mathcal{O}(h^{1.00})`.  The root cause: at :math:`r_{\rm lo} = 0`
the inner face area :math:`A(0) = 0`, so the diamond closure
:math:`\overline\psi = \tfrac12(\psi_{\rm in} + \psi_{\rm out})` gives
:math:`\psi_{\rm out} = 2\overline\psi`, **over-predicting the pole
outer face by exactly +50 %** (mesh-independent rel. error 0.5000), while
the true face is :math:`A(h)` and :math:`2\langle A\rangle_{\rm vol} =
2\cdot\tfrac34 A(h) = 1.5\,A(h)`.  Deeper still, the conservative
*balance itself* is inconsistent at the pole: fed the EXACT cell average
and EXACT inflow it solves for an outer face :math:`-46\,\%` wrong, and
the residual-per-volume plateaus mesh-independently — because
:math:`A_{\rm in} = 0` degenerates the streaming surface integral while
:math:`V \sim h^3`.

[Hebert2009]_ §3.9.4 and Stacey §9.9 **both** use exactly this plain
diamond + Carlson-starting-direction + symmetry scheme at the central
cell with **no special** :math:`\mathcal{O}(h^2)` **closure, and
neither flags reduced order there**.  First-order at the single pole
cell is the accepted, unflagged behaviour of the standard scheme.

**Part 3 — NOT cleanly fixable by a local closure.**  W2 tested the
volume-weighted linear reconstruction
:math:`\overline\psi = \beta\,\psi_{\rm out} + (1-\beta)\,\psi_{\rm in}`
with :math:`\beta = \tfrac34` at the pole (the value that makes
:math:`\overline\psi` :math:`\mathcal{O}(h^3)`-consistent against
:math:`\langle A\rangle_{\rm vol}` at exact faces).  Validated
end-to-end with a faithful production-sweep monkeypatch (and a
:math:`\beta = \tfrac12`-identity regression guard verified to
:math:`3\mathrm{e}{-16}`): :math:`\beta = \tfrac34` does **NOT** restore
:math:`\mathcal{O}(h^2)` — the pole stays :math:`\mathcal{O}(h)`,
magnitude slightly *worse* (:math:`0.0050 \to 0.0106`), and a full-mesh
:math:`\beta` degrades the interior.  Closure-consistency at exact faces
:math:`\neq` fixed-point accuracy: the propagated face flux couples back
through the balance.  A genuine fix needs a non-local higher-order
central-cell reconstruction the canon does not provide — a linear-
discontinuous (Issue #6), cell-update (#158), or nodal scheme
([WuXieFischer1999]_ NSE 133).

Why it is invisible to L2 and to :math:`k_{\rm eff}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The production ``test_sn_spherical_mms_converges_second_order`` uses the
volume-weighted L2 norm.  The pole :math:`\mathcal{O}(h)` at one cell of
:math:`V \sim h^3` contributes :math:`\sqrt V \sim h^{1.5}` →
:math:`h^{2.5}` to L2 — subdominant.  Both midpoint and volume-average
L2 references converge clean :math:`\mathcal{O}(h^{2.00})`; only the L∞
(pole) is :math:`\mathcal{O}(h)`.  For :math:`k_{\rm eff}`: a reflective
sphere recovers :math:`k_\infty = 1.875` exactly, mesh-independent; a
vacuum sphere converges monotone to :math:`\sim 1.78590` at
:math:`\mathcal{O}(h^{1.48})` (combined pole + outer-face first order;
increments :math:`2\mathrm{e}{-5}` at :math:`n_x = 160`, far below
engineering tolerance).  **This is why #233 needed an L∞ / per-cell
probe to surface** — no L2 or eigenvalue gate could see it.

The cylinder shares the same defect, masked
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cylinder pole vs. **midpoint** is :math:`\mathcal{O}(h^2)`
(1.94 / 1.97 / 1.98) but vs. the **volume average** is
:math:`\mathcal{O}(h)` (0.99 / 0.99 / 1.00) — the SAME diamond
inconsistency, masked by the midpoint comparison: the cylinder's
:math:`r\,dr` (linear) weight puts the volume-centroid at
:math:`\tfrac23 h` while diamond's :math:`\tfrac12 A(h)` happens to
:math:`\approx` the midpoint :math:`A(h/2)`, so the midpoint comparison
the gate uses is *accidentally* :math:`\mathcal{O}(h^2)` for the
cylinder.  The cylinder pole is therefore **not** "clean
:math:`\mathcal{O}(h^2)`" — it is the same :math:`\mathcal{O}(h)`
volume-average defect, hidden by the comparison choice.  Cylinder global
L2 is also clean :math:`\mathcal{O}(h^{2.00})`.

The characterization gate (W2, #233)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per the vv-principles "pin what is TRUE + protect the floor WITHOUT
calcifying the limitation" discipline, W2 ships a **characterization**
gate, not a fix gate
(``tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py``,
commit ``255eba4``):

* **Guarantee tests** (carry ``verifies("dd-curvilinear-scalar", ...)``,
  the :eq:`dd-curvilinear-scalar` cell-update label): global
  volume-weighted L2 is :math:`\mathcal{O}(h^2)` (``orders > 1.9``).
  The sphere is asserted under **both** references — midpoint AND the
  Hébert-3.430 shell-volume-average :eq:`sn-pole-cell-shell-average`,
  built from ``scipy.integrate.quad`` (a trusted-library integrator,
  structurally independent of the solver).  Agreement on the order
  across two structurally-different references proves the L2 order is
  REAL, not a midpoint artifact.
* **Characterization tests** (NO ``verifies`` — they pin a *limitation*,
  not a correctness claim): the pole L∞ order is **lower-bounded only**
  (:math:`> 0.8` — "at least first order, does not regress"), the pole
  is the L∞-dominant cell (fraction :math:`> 0.99`), and the interior
  is clean :math:`\mathcal{O}(h^2)` (:math:`> 1.8`).  **No upper bound**
  on the pole order, so a future LD / nodal scheme that lifts the pole
  to :math:`\mathcal{O}(h^2)` keeps the gate green
  (:math:`2.0 > 0.8`) — the characterization gate pins what is true and
  the regression floor without blocking a legitimate improvement
  (vv-principles anti-patterns #5 / #17).

Measured (sphere :math:`n_{\rm ord} = 16`, ladder
:math:`[40, 80, 160, 320]`): L2 midpoint orders :math:`2.01\times3`; L2
shell-average :math:`2.00\times3`; L∞ (pole) :math:`0.91 / 0.95 / 0.97`;
interior :math:`1.84 / 1.92 / 1.96`; pole fraction :math:`1.00` every
mesh.  Cylinder pole-vs-midpoint :math:`1.94 / 1.97 / 1.98`;
pole-vs-volavg :math:`0.99 / 0.99 / 1.00`.

.. _sn-p1-scattering-curvilinear:

Issue #9 — :math:`P_1` anisotropic SCATTERING in curvilinear (the two unrelated paths)
--------------------------------------------------------------------------------------

A persistent source of confusion in this cluster is that "anisotropic"
names **two structurally unrelated** things in a curvilinear SN solve.
Issue #9 is about the *second*; everything above (#229, the
:math:`\alpha`-dome, the :math:`\tau`-clamp) is about the *first*.

* **Path-(I) — geometric angular redistribution.**  The
  :math:`(1-\mu^2)/r\,\partial_\mu\psi` term (sphere) /
  :math:`\xi^2 B / r` (cylinder), threaded by the Morel--Montry Carlson
  :math:`\alpha`-dome.  This is **P0-only**; the "anisotropy" lives in
  the *angular-flux ansatz*, not in the scattering kernel.  The
  existing curvilinear aniso MMS cases
  (:math:`\psi = (A + \zeta B)/W`) exercise ONLY this path.  #229 is a
  Path-(I) test-design floor.
* **Path-(II) — Legendre SCATTERING moments.**  The
  :math:`P_1{+}` scattering source :math:`R\,\Lambda\,M`
  (``scattering.build_aniso_source``, ``scattering_order ≥ 1``),
  geometry-**agnostic**, wired identically for all geometries through
  :func:`~orpheus.sn.coupled_system.build_within_group_system` (the
  :math:`S` gain of the :math:`(L+C),\,S,\,B` decomposition carries
  :math:`P_1` when
  ``scattering_order = 1``).  No curvilinear test exercised Path-(II)
  before #9 — it is NEW coverage of an existing capability (NO
  ``orpheus/`` change; Path-(II) works as-is).

L0 — the operator-admits trick
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than derive a costly symbolic :math:`P_1`-source MMS, W4 feeds a
*known* anisotropic angular flux :math:`\psi_{{\rm ref},n} = (A + \zeta
B)/W` to the within-group :math:`S` operator at
``scattering_order = 1`` and isolates the :math:`P_1` contribution as
:math:`S_1.\mathrm{apply}(\psi) - S_0.\mathrm{apply}(\psi)`, asserted
**per ordinate** (NOT weight-summed — the :math:`\alpha`-dome
telescopes, vv anti-pattern #8) against a structurally-independent
hand-reference:

* **Sphere** (fully SH-table-INDEPENDENT — :math:`P_1(\mu) = \mu`
  directly):

  .. math::
     :label: sn-p1-sphere-hand-ref

     q_n^{P_1} \;=\; \frac{1}{W}\,3\,\mu_n\,\Sigma_{s1}\,\phi_1,
     \qquad
     \phi_1 \;=\; B(r)\,\frac{\sum_n w_n\,\mu_n^2}{W}.

* **Cylinder** (explicit :math:`Y_1^m` moment-sum, independent of the
  production :math:`R\,\Lambda\,M` einsum):

  .. math::
     :label: sn-p1-cylinder-hand-ref

     q_n^{P_1} \;=\; \frac{1}{W}\,3\,\Sigma_{s1}
                  \sum_m Y_1^m(\Omega_n)\,\phi_1^m,
     \qquad
     \phi_1^m \;=\; \sum_n w_n\,Y_1^m(\Omega_n)\,\psi_n.

Both agree at machine precision (rel. :math:`4.7\mathrm{e}{-15}` sphere /
:math:`5.6\mathrm{e}{-15}` cylinder), with a
``max|S₁−S₀| > 1e-6`` negative control (vv anti-pattern #11 — a dropped
:math:`P_1` makes :math:`S_1 - S_0 \equiv 0` and fails the non-zero
hand-ref match).  **1-group is legitimate here**: this is a
flux-shape / OPERATOR claim (the per-ordinate :math:`P_1` source reads
:math:`\phi_1`, flux-shape-dependent by construction), NOT an eigenvalue
claim — the 1-group-degeneracy rule applies only to *eigenvalue*
verification.

L1 — the directional eigenvalue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Forward-peaked :math:`P_1` scattering (:math:`\bar\mu > 0`) **lowers**
:math:`k_{\rm eff}` versus :math:`P_0`.  The physics: positive
:math:`\bar\mu` preserves the forward direction, so in a finite,
vacuum-bounded sphere forward-preserved scattered neutrons are more
likely to cross the outer boundary → **enhanced leakage** → lower
:math:`k_{\rm eff}`.  This requires a **vacuum** outer BC (a reflective
sphere has no leakage → :math:`P_0 \equiv P_1`).  Validated robust:

* Homogeneous vacuum sphere :math:`R = 4 / 10 / 25`:
  :math:`\Delta = k_{\rm eff}^{P_1} - k_{\rm eff}^{P_0} =
  -3.76\mathrm{e}{-3} / -1.32\mathrm{e}{-3} / -2.88\mathrm{e}{-4}` — sign
  always negative, :math:`|\Delta|` **grows as the sphere shrinks**
  (the leakage-monotone signature, the structural negative control a
  sign-flipped or absorption-mimicking :math:`P_1` would violate).
* Heterogeneous fuel-core(:math:`r < 5`)+moderator-shell vacuum sphere
  :math:`R = 10`: :math:`\Delta = -1.40\mathrm{e}{-2}` (:math:`140\times`
  the :math:`1\mathrm{e}{-3}` detection bar), with materials
  ``get_mixture("A","2g")`` (the only fissile 2-group mixture;
  asymmetric downscatter-only P0 avoids the 1-group degeneracy) and
  ``get_mixture("C","2g")``.

Two L1 rows pin this: a heterogeneous-sphere
:math:`k_{\rm eff}^{P_1} < k_{\rm eff}^{P_0}` AND
:math:`1\mathrm{e}{-3} < (P_0 - P_1) < 5\mathrm{e}{-2}`; and a
leakage-monotone control
:math:`(P_0 - P_1)|_{R=4} > (P_0 - P_1)|_{R=25} > 0` (the mechanism
pin).  These are the **first curvilinear exercise** of the
geometry-agnostic ``pn-scatter`` / ``flux-moments`` labels (prior tests
were 2-D Cartesian only).  L0 lands in
``tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py``,
L1 in ``tests/sn/eigenvalue/test_keff_curvilinear.py::TestSphereP1DirectionalEigenvalue``
(commit ``d5878e9``).  L2 is deferred (subsumed by L0+L1; a
:math:`P_1`-convergence L2 needs the :math:`\sigma_{s1}`-MMS source and
rides the same #229 floor).

Infrastructure retained
-----------------------

Per the aggressive-retirement exception, the program deletes no correct
machinery:

.. list-table:: Curvilinear aniso program — primitives status
   :header-rows: 1
   :widths: 34 18 48

   * - Primitive
     - Status
     - Why kept / what changed
   * - Spherical :math:`\tau_m^{\rm raw}` (unclamped)
     - **production**
     - W1: the unique exact-on-linear weight; single-sourced in
       :func:`~orpheus.geometry.reduced_operator.spherical_streaming`,
       inherited by SI sweep + Krylov matvec.
   * - Cylindrical :math:`\tau_m` clamp
     - **production**
     - Retained — the :math:`\tau^{\rm raw} = 0` structural
       :math:`\div 0` block; removing it needs a 2-D
       :math:`(\eta,\varphi)` closure (#229).
   * - Pole-cell characterization gate
     - **regression net**
     - Pins the inherent :math:`\mathcal{O}(h)` pole limitation
       (lower-bounded, not calcified) + the global :math:`\mathcal{O}(h^2)`
       guarantee under two independent references.
   * - Shell-volume-average reference :eq:`sn-pole-cell-shell-average`
     - **oracle**
     - The Hébert-3.430 finite-volume unknown; the principled MMS
       reference that removes the :math:`\sim 75\,\%` comparison
       artifact.  Built from ``scipy.integrate.quad``.

Open research paths (research-tag, not production-blocking)
-----------------------------------------------------------

#. **Higher-order central-cell spatial scheme** (lifts the #233 pole
   :math:`\mathcal{O}(h)`).  The canon ([Hebert2009]_ §3.9.4, Stacey
   §9.9) provides no drop-in :math:`\mathcal{O}(h^2)` central-cell
   diamond closure; the documented route is a non-local higher-order
   *spatial* scheme — linear-discontinuous (Issue #6), step-
   characteristic, or the Green's-function nodal method of
   [WuXieFischer1999]_ (NSE 133, "very high precision on coarse meshes
   relative to standard fine-mesh DD").  Likely diagnostic probe: the
   pole-cell per-cell rate under the shell-average reference, holding
   quadrature fixed.
#. **2-D** :math:`(\eta, \varphi)` **cylinder angular closure** (lifts
   the #229 cylinder floor).  The 1-D :math:`\eta`-thread cannot
   represent the duplicate-azimuthal-:math:`\eta` variation of product /
   level-symmetric quadratures; a genuine 2-D angular closure (or a
   Gauss-type azimuthal quadrature with distinct :math:`\eta` values,
   GitHub Issue #1) is required.  Likely probe: the floor-scaling table
   above with the azimuthal quadrature replaced by a distinct-:math:`\eta`
   set.

Session trail (V&V audit trail)
-------------------------------

* **Commits** (branch ``fix/curvilinear-aniso-pole-and-clamp``,
  2026-06-13): ``b2d8a6d`` (W1 sphere unclamp), ``255eba4`` (W2 #233
  pole-cell characterization gate), ``d5878e9`` (W4 #9 :math:`P_1`
  scattering coverage), ``679a1e6`` (W3 #229 gate retune).
* **Diagnostics**: the W1–W4 ``diag_01..31`` decomposition probes (the
  decisive ones: the
  :math:`E_{\rm test} = E_{\rm artifact}(\text{midpoint} -
  \text{volavg}) + E_{\rm true}(\text{solver} - \text{volavg})`
  decomposition; the discrete-balance residual fed exact fields; the
  faithful production-sweep monkeypatch with a :math:`\beta = \tfrac12`
  identity guard).
* **Literature**: [Hebert2009]_ §3.9.4, Stacey §9.9 (plain diamond at
  the central cell, no special closure), [BaileyMorelChang2010]_ Eq. 43
  (the exact-on-linear weight), [WuXieFischer1999]_ (the nodal route to
  :math:`\mathcal{O}(h^2)` at the origin).
* **vv catalogue**: ``error_catalog.md`` — ERR-059 (the pole-cell
  inherent limitation) + the :math:`\tau`-clamp mis-citation finding +
  the ERR-026 surviving-manifestation note.
* **Issues**: #229 (cylinder floor + sphere gate retune), #9
  (:math:`P_1` curvilinear scattering), #233 (pole-cell, stays OPEN to
  track the future higher-order scheme).

.. note:: **vv-status (eq-labels added by this section).**  The labels
   :eq:`sn-tau-mm-raw`, :eq:`sn-pole-cell-shell-average`,
   :eq:`sn-p1-sphere-hand-ref`, and :eq:`sn-p1-cylinder-hand-ref` are
   *structural / representational* identities (the literature-
   transcribed M-M weight; the Hébert-3.430 finite-volume unknown; the
   structurally-independent :math:`P_1` hand-references).  They are NOT
   solver claims.  The :math:`\tau`-clamp / pole-cell / :math:`P_1`
   *verifiable* content is the W1 clamp-silence + positivity gates, the
   W2 ``verifies("dd-curvilinear-scalar")`` guarantee tests, and the W4
   ``verifies("pn-scatter","flux-moments")`` per-ordinate operator-
   admission gates named above — so these eq-labels are ``documented``.





InvertibleOperator: the streaming-collision operator algebra
==============================================================

The streaming-collision loss composite :math:`L+C` — its operator
surface (``apply`` / ``solve`` / ``apply_transpose``), the one-walk
matvec ≡ sweep fact, reciprocity, and the typed carrier — is
documented in :doc:`slab_one_group` (:ref:`sn-streaming-operator`).
This section retains the historical record of the superseded
two-closure design (kept for the ERR-026 closure-bias reasoning it
records) and the scattering adjoint.

Apply and solve use **different** closures by design (historical)
-----------------------------------------------------------------

.. note:: **Superseded architecture (Wave D / early Wave E).**  This
   subsection describes a design in which ``apply`` (a separate
   finite-difference operator) and ``solve`` (the WDD sweep) were
   **two distinct discretisations** of :math:`A`.  That split was
   dissolved by the #206 Phase C **matvec ≡ sweep** unification: today
   ``apply`` and ``solve`` run the **one** loss-representation walk
   (:meth:`~orpheus.sn.loss_representation.LossRepresentation.loss_action`
   for the matvec, :meth:`~orpheus.sn.loss_representation.LossRepresentation.sweep`
   for the inverse — "one walk", a code fact per L21), so there is no
   longer a separate FD operator and no by-design bit-difference between
   the two directions.  The FD-operator family
   (``transport_operator_matvec*`` and its unified successor) was deleted
   in the typed-field (#197) and walk-unification (#280 campaigns)
   refactors.  The historical two-closure narrative below is retained for
   the ERR-026 closure-bias reasoning it records.

This was the load-bearing architectural fact about
:class:`InvertibleOperator`, and the reason the operator's
:meth:`apply` was **not** bit-identical to its :meth:`solve`.

The Wave-D-era SN dispatch in ORPHEUS shipped **two distinct
discretisations** of the same continuous operator
:math:`A = \Omega\cdot\nabla + \Sigma_t`, built at different times
for different consumers:

* The **finite-difference operator** (the ``transport_operator_matvec_*``
  family in :mod:`orpheus.sn.operators.streaming`, since deleted) was
  built for the BiCGSTAB inner
  solver path (:meth:`SNSolver._solve_bicgstab_*`).  It used
  upwind cell-center FD on Cartesian and arithmetic face averages
  with **τ-symmetric Morel-Montry angular interpolation** on
  curvilinear (see the "Explicit Transport Operator" subsection of
  the BiCGSTAB Alternative above and the warning at the head of
  :mod:`orpheus.sn.operators.streaming`).

* The **sweep operator**
  (:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`,
  dispatching its selected representation through the
  :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme`
  Protocol with :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
  as the default strategy) uses the **WDD asymmetric closure**
  :math:`\psi_{n+1/2} = (\overline\psi - (1-\tau)\,\psi_{n-1/2})/\tau`.
  This is the historical SN sweep's closure: the upper-triangular
  forward-substitution form that lets a sweep run in
  :math:`O(N\cdot N_{\rm cells})` work.

In the fine-mesh limit both discretisations converge to the same
continuous operator.  On coarse meshes they differ:

* For Cartesian the difference is :math:`O(h)` (upwind FD has
  the same first-order consistency as DD on uniform meshes;
  divergence appears on non-uniform meshes — see the warning at
  the head of :mod:`orpheus.sn.operators.streaming`).
* For curvilinear the WDD asymmetric closure has a closure-bias-
  driven self-consistent fixed point that is **not** the
  fine-mesh-limit transport solution (ERR-026).

The reconciliation is **Wave E Issue 15**: the SN solver's inner
loop migrates from "sweep with WDD closure" to "Krylov on apply,
sweep as preconditioner".  Krylov on :meth:`apply` uses the
symmetric closure (the one that agrees with analytical references)
as the system to solve; the WDD sweep is invoked only as a
preconditioner that accelerates the Krylov iteration without
poisoning the converged solution with its closure bias.  ERR-026
closes when Wave E lands; the 2 xfail-strict tripwires at
:file:`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
are the gating bug-catchers for that closure.

.. _sn-scattering-adjoint:

The scattering adjoint, free from the harmonic frame
----------------------------------------------------------

The streaming operator's analytic adjoint is hard (the subsection above):
sign-flipping the upwind direction, transposing the M–M closure, re-deriving
the per-level azimuthal redistribution — each an AI-failure-mode trap — so
:math:`A^*` is taken by the dense-transpose fallback.  The **scattering**
operator :math:`S` is the counterexample: campaign **#276 P3** (commit
``15185e5``, closes
`#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_) made
:math:`S^{T}` fall out **for free**, because :math:`S` is already written as
a harmonic-frame conjugation.

The modernised in-scatter source is ONE frame-conjugated operator
(:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`):

.. math::

   \mathrm{full\_scatter\_kernel}
   \;=\; R \circ (\Lambda_{\ell\ge 0} + N_{2n}) \circ M ,

where :math:`M` / :math:`R` are the angular frame's analysis /
reconstruction faces, :math:`\Lambda_{\ell\ge 0}` is the per-:math:`\ell`
moment-space group transfer
(:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`),
and :math:`N_{2n}` is the distinct :math:`(n,2n)` multiplication channel
(:class:`~orpheus.transport.operators.scattering.N2NMomentOperator`) —
summed with :math:`\Lambda` in moment space and conjugated by the frame
*together* (one analysis, one reconstruction) for the WHOLE
P0 + :math:`\ell\ge1` + :math:`(n,2n)` source.  Its transpose is therefore
the product transpose

.. math::

   \mathrm{full\_scatter\_kernel}^{T}
   \;=\; M^{T} \circ (\Lambda + N_{2n})^{T} \circ R^{T},

which :meth:`OperatorProduct.apply_transpose
<orpheus.numerics.operator.OperatorProduct.apply_transpose>` assembles from
the leaf transposes — the frame's :math:`M^{T}` / :math:`R^{T}` faces (landed
in the Frame/Basis carve), the per-:math:`\ell` group transpose
:math:`\Lambda^{T}`, and :math:`N_{2n}^{T}` — with **no per-geometry
derivation to verify** (the trap the streaming adjoint above could not
avoid).  The per-ordinate adjoint scattering source is then

.. math::

   S^{T}\chi \;=\; \tfrac{1}{W}\,\mathrm{full\_scatter\_kernel}^{T}\,\chi ,

the producer-side :math:`1/W` transposing as the scalar it is
(:math:`(A/W)^{T} = A^{T}/W`).
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` now
reports ``is_adjointable = True`` (it has a working ``apply_transpose``),
and the old "no ``apply_transpose``" class-docstring confession is
retired.

**Forward fast-path, adjoint frame-path — and why the asymmetry is
principled.**  The production FORWARD source keeps the scalar fast-path
(:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_kernel`
for P0 + :math:`(n,2n)`, and the per-:math:`\ell` ``build_aniso_source``)
for SI-sweep performance; the **adjoint** — not the hot path — rides the
validated frame form instead.  The two are thus structurally *different*
representations of the same operator, which is exactly what makes the
verification a genuine cross-check rather than a tautology: the per-group
Euclidean reciprocity
:math:`\langle S\psi, \chi\rangle = \langle\psi, S^{T}\chi\rangle`
(``tests/sn/operators/test_scattering_adjoint.py``,
``TestFullScatterKernel::test_S_euclidean_reciprocity``) pins the frame-form
:math:`S^{T}` against the *independent* scalar fast-path :math:`S`, and the
forward equivalence
:math:`(1/W)\,\mathrm{full\_scatter\_kernel}.\mathrm{apply} \equiv
S.\mathrm{apply}` holds to :math:`\sim 10^{-12}`.

.. note::

   This :math:`S^{T}` is the **Euclidean** transpose (the plain
   group-and-angle matvec adjoint), NOT the metric Hilbert adjoint
   :math:`S^{\dagger} = G^{-1}S^{T}G` — that angular-Gram weighting is the
   :attr:`~orpheus.numerics.operator.LinearOperator.H` wrapper's job.  The
   campaign and commit name it "S†" colloquially; the precise object the
   operator computes is the transpose.

This is the discrete scattering adjoint the SN adjoint chain builds on: the
adjoint flux :math:`\psi^{*}` solving :math:`(L+C-S)^{T}\psi^{*} = q^{*}`,
adjoint-weighted homogenisation, perturbation theory, and detector
sensitivity all need :math:`S^{T}`.  Its companion forward step (campaign
**#276 P2**, commit ``dcea43a``) routes the SN forward *isotropic* source
through the same model-shared
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
(:math:`\Sigma_{s0}`) and
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
(:math:`2\Sigma_{2n}`) operators (0-ULP bit-identical), so the
:math:`K_\mathrm{iso}` energy operator — which also assembles the
infinite-medium loss matrix (:ref:`direct-eigensolve-assembly`) — is one
cross-model source.  These model-shared operators live in
:mod:`orpheus.transport.operators`.


Wave E and beyond — landed and forward
--------------------------------------

* **Wave E Issue 14 lifted** the iteration primitives (power
  iteration, source iteration) into stand-alone operator-algebra
  consumers in :mod:`orpheus.numerics.iteration`, decoupling them
  from :class:`SNSolver`.
* **Wave E Issue 15 wired** :class:`SNSolver` to
  :class:`InvertibleOperator`: the inner solve becomes Krylov on
  :meth:`apply` with sweep as preconditioner, which removes the WDD
  asymmetric closure from the converged-solution path for the
  reflective-BC eigenvalue case.  Wave E Round 3 will extend the
  equation-map layout to the vacuum-BC path that closes ERR-026
  for fixed-source curvilinear MMS.
* **Forward**: when production reciprocity becomes performance-
  critical, an :math:`O(n)` analytic-adjoint matvec replaces the
  dense-transpose fallback; the new path is gated by the
  reciprocity tests in :file:`tests/sn/test_streaming_operator.py`
  (post-D-K successor) and :file:`tests/sn/test_phase_c_gates.py`
  Gate 1.3 (xfail pending Wave O).

References
----------

* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport.*  §10 (adjoint transport, reciprocity, and
  perturbation theory).
* Adams, M. L., & Larsen, E. W. (2002).  *Fast iterative methods
  for discrete-ordinates particle transport calculations.*
  Progress in Nuclear Energy 40 (1), 3–159.  Reviews
  Krylov-on-apply with sweep preconditioning.
* Trefethen, L. N., & Bau, D. (1997).  *Numerical Linear Algebra.*
  §3.2 (matrix-free Krylov view).


.. _sn-solver-operator-algebra-coordinator:

SNSolver as an operator-algebra coordinator
============================================

:class:`~orpheus.sn.solver.SNSolver` consumes the operator triple
:math:`(A, S, F)` directly.  At construction time, the solver builds:

* :attr:`SNSolver.L` — :class:`InvertibleOperator` carrying the
  symmetric-closure streaming-collision operator.
* :attr:`SNSolver.S` — :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
  carrying the P0 + (n,2n) + Pℓ Galerkin reconstruction (Wave D
  Issue 13).
* :attr:`SNSolver.F` — :class:`~orpheus.transport.operators.fission.FissionOperator`
  carrying the rank-1-in-energy fission emission (Wave D Issue 13).

Each of these is a :class:`~orpheus.numerics.operator.LinearOperator`
in the Wave A operator-algebra sense: predicate-typed, composable
under :class:`~orpheus.numerics.operator.OperatorSum` and
:class:`~orpheus.numerics.operator.OperatorProduct`, and protocol-
conforming so the iteration primitives in
:mod:`orpheus.numerics.iteration` consume them without SN-specific
plumbing.  The within-group inner solve is built once from a single
source of truth — the :func:`~orpheus.sn.coupled_system.build_within_group_system`
builder assembles the :class:`~orpheus.sn.coupled_system.WithinGroupSystem`
record, the honest within-group decomposition :math:`(L+C,\ S,\ B)` as a
named splitting :math:`A = M - N`: the invertible resolvent
:math:`M = (L+C)` plus its two lagged coupling gains :math:`N = (S,\ B_a)`
(the bulk scattering :math:`S` and the trace boundary reflection :math:`B`;
zero within-group fission), handed to the **variadic** driver
:math:`\text{Driver}(L_{\rm resolvent},\,*\text{gains})` (Wave O step
O.2a — the transitional :math:`S + B` fold is retired; see
:ref:`bc-extraction-variadic-driver` in :doc:`/theory/foundations/boundary_conditions`).
:func:`_within_group_krylov` wraps the matching
:class:`~orpheus.numerics.iteration.KrylovAcceleration` — and the
decomposition is shared verbatim across the eigenvalue source-iteration
inner (:meth:`SNSolver._solve_source_iteration`), the eigenvalue Krylov
inner (:meth:`SNSolver._solve_krylov`), and both fixed-source paths.

The within-group inner solve consumes the primitives directly
-------------------------------------------------------------

:class:`SNSolver`'s within-group inner solve **is** the
:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KrylovAcceleration` primitive — not
a verbatim replica of its loop.
:meth:`SNSolver._solve_source_iteration` constructs a
:class:`SourceIteration` from the :func:`~orpheus.sn.coupled_system.build_within_group_system` SSoT and
runs it; :meth:`SNSolver._solve_krylov` constructs a
:class:`KrylovAcceleration` from :func:`_within_group_krylov` and runs
that.  The Layer-3 resolvent of the SN row in the
:ref:`eigenvalue-posing` architecture is exactly these primitive
instances.

The primitive is **type-agnostic and angular-capable**: it operates on
the typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
composite, which carries the full angular flux on its bulk.  Pℓ
anisotropic scattering therefore rides the angular bulk with no special
plumbing — :meth:`ScatteringOperator.apply` on the timeless
:class:`~orpheus.transport.full_field.FullField` operator carrier (the
driver's :class:`~orpheus.transport.timed_full_field.TimedFullField` iterate
reaches it via MRO) reads the angular moments off the composite and builds
the anisotropic source via :meth:`ScatteringOperator.build_aniso_source`,
all inside the primitive's normal RHS path.  There is **no scalar-flux
limitation** and **no pending "Approach A" cleanup**: the earlier
framing — that :class:`SourceIteration` carried only scalar flux and SN
had to replicate the loop verbatim until the angular state could be
threaded through — was a property of an interim scalar-only carrier
that the typed composite retired.  The
``.claude/skills/algebra-of-record`` "Branch 2 implements the same
operator algebra" discipline is satisfied: SN is the discretized
Branch-2 consumer of the shared primitive, not a parallel loop.

The (L − S − F)·ψ = (1/k)·F·ψ framing at the solver level
-----------------------------------------------------------

Beyond driving the within-group inner solve, the :math:`(A, S, F)`
framing organises the solver's outer API surface:

* :meth:`SNSolver.compute_fission_source` returns
  :math:`F\,\phi/k` — a thin delegator to :meth:`F.apply` with the
  :math:`1/k` outer-loop scaling applied at the solver level.
* :meth:`SNSolver.solve_fixed_source` solves
  :math:`(A - S)\,\psi = q_{\rm ext}` (with :math:`q_{\rm ext}` the
  fission source built by ``compute_fission_source``).  Two paths:

  * ``inner_solver="source_iteration"`` — sweep-driven fixed-point
    iteration; :math:`A^{-1}` is the WDD asymmetric sweep.  ERR-026-
    affected for curvilinear vacuum-BC cases.
  * ``inner_solver="krylov"`` — GMRES on :meth:`A.apply` with the
    sweep as preconditioner.  Reflective-BC equation map only
    (Round 3 owns the vacuum-BC extension).

* :meth:`SNSolver.compute_keff` returns **fission production over net
  removal**, :eq:`sn-keff-update` — the volume-weighted method-layer
  functional :math:`R_{\nu\Sigma_f}(\phi) / (R_{\Sigma_a}(\phi) + L -
  E_{2n}(\phi))`, derived in :ref:`sn-keff-estimator` below.  The
  SN-specific volume weighting lives here (in the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  fields) — one honest realization of the same discipline the
  operator-form :meth:`KEigenvalue.compute_keff` spells with the
  measure absorbed into the operators' action.  (Pre-#291 this method
  returned the leakage-blind :math:`\sum F\phi V / \sum \Sigma_a\phi V`
  ratio; see :ref:`sn-keff-estimator` for why that was a
  non-eigenvalue on any vacuum-bounded problem.)

The solver-level :math:`1/k` scaling (in
:meth:`~SNSolver.compute_fission_source`) and the volume-weighted
eigenvalue estimate (in :meth:`~SNSolver.compute_keff`) are exactly the
points where SN's specifics live; the rest of the solver is
operator-algebra coordination over the canonical
:func:`~orpheus.numerics.eigenvalue.power_iteration` boundary.  These
two K-specific hooks are also precisely *why* the Layer-4 loop is not
yet literally K/α-agnostic — relocating the eigenvalue scaling into the
algorithm is the first step of the α-wave (see the honest-scope caveat
in :ref:`eigenvalue-posing`).

The eigenvalue :math:`\keff` is determined by **power iteration**: an
outer loop updates :math:`k` from the net-removal balance
:eq:`sn-keff-update` (fission production over absorption + leakage −
:math:`(n,2n)` emission), with an inner loop that solves the
within-group scattering problem.

.. _sn-keff-estimator:

The reported eigenvalue: fission production over net removal
------------------------------------------------------------

:meth:`~orpheus.sn.solver.SNSolver.compute_keff` reports the eigenvalue
of the problem the inner solve **actually poses**.  This is the SN
symptom (#291) and the MoC/CP/homogeneous root (#259) of a single
discipline: *the reported* :math:`k` *must be the eigenvalue of the
fixed-source map every method scales only fission by* :math:`1/k`
*through* — scattering and the :math:`(n,2n)` emission are plain gains
assembled **inside** :meth:`~orpheus.sn.solver.SNSolver.solve_fixed_source`,
never scaled by :math:`1/k`.  An estimator that disagrees with its own
posed problem converges cleanly and silently to a **non-eigenvalue
ratio**.

.. math::
   :label: sn-keff-update

   k \;=\; \frac{R_{\nu\Sigma_f}(\phi)}
                {R_{\Sigma_a}(\phi) \;+\; L \;-\; E_{2n}(\phi)}

.. (vv-status rationale) Governing/definitional identity: the reported k
   IS the eigenvalue of the posed fixed-source map, not a solver
   eigenvalue-correctness claim against an external analytical reference
   (that rests on the multi-group heterogeneous L1/L2 references
   elsewhere on this page). The verifiable content is the cross-engine
   consistency gate tests/sn/eigenvalue/test_keff_estimator_gate.py
   (reported k == the converged fixed-point map ratio k* = P(Mφ*)/P(φ*),
   map-ratio ground-truth noise ≤ 2e-11) with in-file mutation teeth.
   Wiring @pytest.mark.verifies("sn-keff-update") onto that gate is a
   test-side follow-up (this docs pass could not touch tests/).
.. vv-status: sn-keff-update documented

The three terms are typed volume-integrated reaction-rate functionals
and one boundary functional:

* **Numerator** :math:`R_{\nu\Sigma_f}(\phi) = \int_V \nu\Sigma_f\,\phi\,dV`
  — the fission production, the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  over :math:`\nu\Sigma_f` (the :math:`\phi^\dagger\!=\!1` degenerate of
  the homogenization Petrov–Galerkin bilinear).  The :math:`(n,2n)`
  emission is **not** production here — contrast
  :meth:`~orpheus.sn.solver.SNSolver.compute_production_rate`, the
  ERR-052 renormalisation scale anchor, which keeps *total* physical
  production (fission **plus** the :math:`(n,2n)` emission).  The role
  split is the load-bearing #259 correction: the same physical
  :math:`(n,2n)` neutrons are a **scale** quantity for the outer
  renormalisation but a **removal-reduction** in the eigenvalue balance.
* **Absorption** :math:`R_{\Sigma_a}(\phi) = \int_V \Sigma_a\,\phi\,dV`,
  with :math:`\Sigma_a = \Sigma_f + \Sigma_c + \Sigma_L +
  \sum_{g'}\Sigma_{2,g\to g'}` — i.e. ``absorption_xs`` counts the
  :math:`(n,2n)` **collision once** (the neutron is removed from its
  incident group by the collision).  See
  :attr:`~orpheus.data.macro_xs.mixture.Mixture.absorption_xs`.
* **Leakage** :math:`L` — the net vacuum-boundary outflow (below).  On a
  reflective (lattice) problem it is a **structural zero**.
* **Emission** :math:`E_{2n}(\phi) = \int_V \sum_{g,g'} 2\,\Sigma_{2,g'\to
  g}\,\phi_{g'}\,dV` — the :math:`(n,2n)` **emission** (two neutrons out
  per collision; the factor 2).  A gain, so it **reduces** net removal.

The net :math:`(n,2n)` effect on removal is therefore
:math:`\underbrace{\sum_{g'}\Sigma_{2,g\to g'}}_{\text{collision, in }\Sigma_a}
- \underbrace{2\Sigma_2}_{E_{2n}} = -\Sigma_2` — **one extra neutron
gained** per collision, exactly the physics of a neutron-doubling
reaction.

**Balance identity (divergence-telescoping).**  The angle- and
group-summed discrete cell balance for cell :math:`i` in the posed
eigenproblem is

.. math::

   \underbrace{\sum_{f\in\partial i}\!\bigl(\textstyle\sum_g J_g\bigr)\,\Delta A_f}
              _{\text{net face flow}}
   \;+\; \Sigma_{t,i}\,\phi_i\,V_i
   \;=\; \frac{1}{k}\,R_{\nu\Sigma_f,i}
        \;+\; \Sigma_{s,i}\,\phi_i\,V_i
        \;+\; E_{2n,i}

(streaming + total collision on the left; the isotropic fission source
scaled by :math:`1/k`, plus the *unscaled* scatter and :math:`(n,2n)`
gains, on the right).  Summing over **all** cells, every interior face
is shared by two cells with opposite outward normals and equal current
(continuity), so the interior face-flow terms **telescope to zero** —
only the domain-boundary faces survive, and their sum is the net
leakage :math:`L`.  With :math:`\Sigma_t - \Sigma_s = \Sigma_a` this
collapses to

.. math::

   \frac{R_{\nu\Sigma_f}(\phi)}{k} \;=\; R_{\Sigma_a}(\phi) \;+\; L
                                        \;-\; E_{2n}(\phi),

which is :eq:`sn-keff-update` rearranged.  This is the same discrete
divergence discipline the diffusion page states as
:math:`\mathbf 1^{\mathsf T}(C-S)=\Sigma_a` with interior-face
telescoping (see :ref:`diffusion-leakage-boundary-leaves`); SN and
diffusion report the *same* balance-law eigenvalue, differing only in
how the streaming operator is discretised.

The leakage functional
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   L \;=\; \sum_{f\,\in\,\text{vacuum}} \oint_{f} dA\,
           \sum_g J_g(\mathbf{r}_f)\,,
   \qquad
   J_g \;=\; \sum_m (\Omega_m\cdot\hat n_f)\, w_m\, \psi_{m,g}

is the face-area integral of the boundary trace's **net outward
current**.  The angular-to-scalar reduction :math:`J_g` is
:meth:`AngularBoundaryFlux.net_current
<orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current>`
— the single source of the :math:`\Omega\cdot\hat n\,w` contraction, the
angular sibling of the scalar trace's :math:`J = J^+ - J^-`
(:meth:`ScalarBoundaryFlux.net_current
<orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux.net_current>`).
It is spelled through the trace space's own atoms — the signed
projection table
:attr:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.omega_dot_n`
and the :math:`|\Omega\cdot\hat n|\odot w` partial-current metric (using
the identity :math:`\operatorname{sign}(\Omega\cdot\hat n)\cdot
|\Omega\cdot\hat n|\,w = \Omega\cdot\hat n\,w`) — so no consumer
re-derives the cosine weighting.  Tangential ordinates carry zero weight
and drop out.

The face measure :math:`dA` is supplied by
:meth:`SNSolver._face_area_of`, matching the cell
``volume_measure`` exactly so the balance identity closes:

.. list-table:: Boundary-face measure by geometry
   :header-rows: 1
   :widths: 30 30 40

   * - Geometry
     - Face measure :math:`\Delta A`
     - Source
   * - 1-D slab
     - :math:`1` (per unit cross-section)
     - :attr:`MaterialMesh.areas <orpheus.transport.mesh.material_mesh.MaterialMesh.areas>`
   * - 1-D cylinder
     - :math:`2\pi R` (per unit height)
     - ``MaterialMesh.areas``
   * - 1-D sphere
     - :math:`4\pi R^2`
     - ``MaterialMesh.areas``
   * - 2-D Cartesian
     - transverse edge-cell widths (unit depth)
     - ``mesh.axes`` transverse extent
   * - 3-D Cartesian
     - :math:`\Delta A_{\mathbf c} = \prod_{j\ne a}\Delta_j[c_j]`
       (transverse-area outer product)
     - ``mesh.axes`` transverse extents

The :math:`d \ge 2` Cartesian arms are ONE generic body: the outer
product of the *other* axes' edge widths in **ascending axis order** —
the same codimension-1 enumeration as
:func:`~orpheus.transport.mesh.axis.face_shape`, so the measure array
broadcasts cell-for-cell against the ``(ng, *face_spatial)`` net
current, and the 2-D width vector is just the single-transverse-axis
degenerate (bit-identical to the pre-3-D spelling).

The 3-D arm originally shipped as a **typed refusal**
(``NotImplementedError``): guessing the transverse product's cell
ordering would silently mis-weight the leakage sum, and Cardinal Rule 1
forbids returning a wrong-but-plausible number.  The wire landed
2026-07-13 when the first 3-D vacuum eigenvalue consumer arrived (the
d=3 Mode-9 G-S≡Jacobi gate), with the ordering pinned twice in
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: an **object-level
pin** (face measure ≡ the boundary layer's ``volumes / Δ_axis``, the
mesh's own ascending-axis enumeration — vv Mode-12 discipline: pin the
object, not only the k functional) and the **k* map-ratio gate** on a
Mode-2 asymmetric all-vacuum box, whose teeth are proven by permanent
in-process mutants — a reversed transverse enumeration moves the
reported k by a measured **13.9 %** against the estimator-independent
:math:`k^*` (clean agreement :math:`6\times10^{-10}`), and a transposed
enumeration crash-REDs on the broadcast.  A trace carrying a
``#251`` transverse face-moment tail is refused loudly at the
consumption site (the face integral must consume ONLY the
transverse-average moment — higher Legendre face moments integrate to
zero over each face cell — and that slot-0 read has no consumer yet).

Reflective faces are a structural zero
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reflective law equates a face's inflow to its reflected outflow
**exactly**, so the net current there vanishes *by construction*.
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` therefore
**skips** reflective faces — it never accumulates them, rather than
accumulating a value that ought to be zero but carries
:math:`\pm`-cancelling angular-sum floating-point noise.

This is a deliberate design choice with a bit-level payoff.  On an
all-reflective (lattice) problem :math:`L` is a structural ``0.0``, and
on a :math:`\Sigma_2`-free mixture :math:`E_{2n}` is exactly ``0.0`` (the
per-material :math:`(n,2n)` loop adds nothing), so
:eq:`sn-keff-update` reduces **bit-identically** to the historical
lattice functional ``production / absorption``.  Every pre-existing
reflective eigenvalue anchor is preserved to the last ULP — the
unification adds terms that vanish structurally, not numerically, on the
cases it must not perturb.

The scale bridge: trace of the last inner solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The leakage term reads the **typed** boundary trace of the last inner
solve (``self._psi_typed.boundary``), whereas the numerator/denominator
reaction rates consume the bare-array flux :math:`\phi` the estimator is
handed.  These two representations can be at **different scales**:
:func:`~orpheus.numerics.eigenvalue.power_iteration` renormalises
:math:`\phi` to unit production rate *between* the inner solve and the
:math:`k`-update (ERR-052), so the stored trace belongs to the
**un-renormalised** last iterate while the estimator's :math:`\phi` is
its renormalised multiple.

Leakage is degree-1 homogeneous in :math:`\psi`, so the fix is a single
rescale by the fission-production ratio of the two fluxes
(``self._phi_of_trace``, stored alongside the trace at **both**
inner-path returns) — exactly ``1.0`` when the caller passes the
returned flux itself.  The **contract** is therefore: the flux handed to
:meth:`~orpheus.sn.solver.SNSolver.compute_keff` must be a scalar
multiple of the last inner solve's flux (true for ``power_iteration`` and
for every manual solve-then-estimate loop).

If a vacuum face exists but **no** inner solve has stored a trace,
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` raises
``RuntimeError`` — the leakage cannot be answered honestly, and silently
returning it as zero would *reproduce the #291 omission*.  Fail loud;
never return a non-eigenvalue.

The R7 :math:`(n,2n)` convention fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The historical spelling put the :math:`(n,2n)` emission in the
**numerator** as production,

.. math::

   k_{\text{old}} \;=\; \frac{R_{\nu\Sigma_f} + E_{2n}}{R_{\Sigma_a}},

which is a **non-eigenvalue** of the posed map whenever
:math:`\Sigma_2 \neq 0` *and* :math:`k \neq 1`.  The reason is exactly
the posing asymmetry: the inner solve scales **only** fission by
:math:`1/k`; the :math:`(n,2n)` emission is an *unscaled* gain in the
sweep source.  So the eigenvalue of that map is
:math:`k^\star = R_{\nu\Sigma_f}/(R_{\Sigma_a} - E_{2n})` (reflective,
:math:`L=0`), and putting the *unscaled* emission in the numerator does
not recover it.  Writing :math:`f = R_{\nu\Sigma_f}`,
:math:`a = R_{\Sigma_a}`, :math:`e = E_{2n} = 2s_2` and substituting
:math:`f = k^\star (a - e)`:

.. math::

   k_{\text{old}}
   \;=\; \frac{k^\star (a - e) + e}{a}
   \;=\; k^\star \;+\; \frac{2 s_2\,(1 - k^\star)}{a}.

The two agree only when :math:`s_2 = 0` (no :math:`(n,2n)`) or
:math:`k^\star = 1` (critical).  For a supercritical
:math:`k^\star > 1` the correction is negative
(:math:`k_{\text{old}} < k^\star`); for subcritical, positive.  The MoC
and CP pages carry the same fork
(:eq:`moc-keff-update`, :eq:`cp-keff-update`); CP was the one member
already spelled on net removal.

What was tried and found
~~~~~~~~~~~~~~~~~~~~~~~~~~

The #291 bias was characterised pre-fix (commit ``d1daaac``, Gauss–
Legendre :math:`n=8`, map-ratio ground truth noise :math:`\le 2\times
10^{-11}`) across the five gate configurations:

.. list-table:: Pre-fix reported :math:`k` vs the posed-problem eigenvalue :math:`k^\star`
   :header-rows: 1
   :widths: 40 18 18 24

   * - Configuration
     - Pre-fix reported
     - Posed :math:`k^\star`
     - Bias
   * - homog. 2G vacuum slab (width 8)
     - 1.83767525
     - 0.98163269
     - :math:`+87.2\%` (:math:`L/A = 0.872`)
   * - het. vacuum sphere P\ :sub:`0`
     - 0.86484694
     - 0.70601977
     - :math:`+22.5\%`
   * - het. vacuum sphere P\ :sub:`1`
     - 0.85080423
     - 0.67876772
     - :math:`+25.3\%`
   * - reflective control (:math:`\Sigma_2=0`)
     - 1.87500000
     - 1.87500000
     - :math:`\equiv` (bias :math:`1.2\times 10^{-10}`)
   * - reflective :math:`\Sigma_2\neq 0`
     - 1.92857143
     - 2.61278195
     - :math:`-26.2\%` (R7 defect)

The two failure classes are visible in one table: the vacuum rows are
the **leakage omission** (#291) — the reported :math:`k` overshoots by
the leakage-to-absorption ratio :math:`L/A`; the last row is the **R7
:math:`(n,2n)` convention** — zero leakage, yet a
:math:`-26.2\%` error because the emission was mis-posed as production.
The reflective-control row is exactly the bit-identity guarantee above.
The exact check on the R7 row is
:math:`0.78/(0.5185 - 0.2200) = 2.61278`, and
:math:`(0.78 + 0.2200)/0.5185 = 1.92857` reproduces the old value —
matching :math:`k_{\text{old}} = k^\star + 2s_2(1-k^\star)/a` term for
term.

Post-fix, reported :math:`k` and the map-ratio :math:`k^\star` agree to
:math:`\le 6\times 10^{-10}` on all five.  The P\ :sub:`0`\ –P\ :sub:`1`
sphere gap :math:`\Delta` roughly **doubled** (:math:`1.404\times
10^{-2} \to 2.725\times 10^{-2}`) but stays inside the diagnostic
:math:`(10^{-3}, 5\times 10^{-2})` band — the P\ :sub:`1` anisotropic
correction is now measured against the correct eigenvalue on both
solves.

The V&V decision was a **principled re-baseline** (per ``vv-principles``
bit-identity-vs-principled-equivalence): the old reported :math:`k` was a
*different functional* from the posed problem's eigenvalue, so the new
value is not a regression to be tolerance-matched but a correction to be
verified against a structurally-independent reference (the fixed-point
map ratio).

Verification
~~~~~~~~~~~~

The permanent gate is
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: it asserts the
reported :math:`k` equals the converged fixed-point map ratio
:math:`k^\star = P(M\phi^\star)/P(\phi^\star)` across the four physics
regimes — {vacuum slab, vacuum sphere (pinning the :math:`4\pi R^2`
face-area convention), reflective bitwise-degenerate, reflective
:math:`\Sigma_2\neq 0`} — with **in-file mutation teeth**: a
leakage-drop mutation reds the vacuum legs while staying bitwise-green
on reflective; a leakage sign-flip crash-reds through the scale-bridge
guard; and the old :math:`(n,2n)`-in-numerator convention reds the
:math:`\Sigma_2\neq 0` leg.

This is a **consistency** gate: the map ratio is the structurally-
independent ground truth for "does the estimator return the eigenvalue
of its own posed map", and is blind by construction to *which*
eigenvalue that is.  The SN solver's eigenvalue **correctness** — that
the posed map's eigenvalue is the *physically right* :math:`k` — rests
on the multi-group heterogeneous L1/L2 references elsewhere on this
page, not on this gate.

Two Inner Solvers
-----------------

**Source iteration (sweep-based):**

- Operator: :math:`T^{-1}` (diamond-difference sweep)
- Solution variable: scalar flux :math:`\phi(x, y, g)`
- Fixed-point: :math:`\phi^{(k+1)} = T^{-1}(S \cdot \phi^{(k)} + Q_f)`
- Convergence rate: spectral radius of :math:`T^{-1}S` (~0.97 for 421
  groups)
- Cost per iteration: one transport sweep
- Works for all geometries

**Krylov (direct operator):**

- Operator: :math:`L + C` applied matrix-free
  (:meth:`InvertibleOperator.apply` — the same one-walk
  discretization as the sweep; L21 matvec ≡ sweep)
- Solution variable: angular flux :math:`\psi(x, y, n, g)` (much
  larger than scalar flux)
- System: :math:`(L+C)\psi = b` where :math:`b` = fission + scattering
- Convergence: GMRES with sweep preconditioner, typically ~100
  Krylov iterations at ``tol=1e-4`` (always converges)
- Available for all geometries (Cartesian, spherical, cylindrical)

Wave E Round 2 (Issue #164) replaced the legacy BiCGSTAB FD-operator
path with this Krylov path.  See the Krylov alternative in
:doc:`slab_one_group` for the full discussion.

The two paths share the **one** loss-representation discretization
(matvec ≡ sweep, #206 Phase C), so they converge to the same fixed
point; the Wave-D-era design in which they carried different spatial
closures — and disagreed on coarse-mesh :math:`\keff` — is recorded
in the streaming-collision history section above.


.. _sn-consuming-the-frame:

Consuming the frame in SN
=========================

Spatial homogenisation and energy condensation are **discrete-frame
projections** — the Petrov-Galerkin coefficient extraction
:math:`G^{-1}M` of a flux- (or spectrum-) weighted frame. All of that
theory — rate preservation, the source-group / sink-sum matrix rules, the
metric-fold-vs-bilinear adjoint argument, fractional-overlap re-binning,
the condense/homogenize asymmetry law, and the verification gates — is
the frame page's headline **Petrov-Galerkin** consumer; see
:ref:`sn-spatial-homogenization` and :ref:`sn-energy-condensation`
(:doc:`/theory/foundations/frame`). This section keeps only the **SN-layer orchestration**:
how the SN :class:`~orpheus.sn.solution.Solution` drives that machinery
from a converged flux.

Homogenisation: the solve → homogenize → re-solve loop
------------------------------------------------------

:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
takes a coarse mesh (:class:`~orpheus.geometry.mesh.Mesh1D` or
:class:`~orpheus.geometry.mesh.Mesh2D`) and returns a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` — the coarse
geometry already carrying one freshly-homogenised effective
:class:`~orpheus.data.macro_xs.mixture.Mixture` per coarse cell. The SN
:class:`~orpheus.sn.solution.Solution` owns the converged flux, so the SN
layer is what builds the flux-weighted **test** basis the frame consumes;
the frame itself, and the rate-preservation theory that *forces* the flux
weighting (rather than a plain volume average), live in
:ref:`sn-spatial-homogenization` (:doc:`/theory/foundations/frame`). The returned
``MaterialMesh`` is re-promoted to a solvable phase space by
:meth:`SNMesh.from_material_mesh
<orpheus.sn.mesh.augmented_mesh.SNMesh.from_material_mesh>`, closing the
**solve → homogenize → re-solve** loop. The return type is
**mesh-coupled** (geometry and materials born together) — the space half
of the condense/homogenize asymmetry law
(:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`).

Condensation: per-material representative spectra
-------------------------------------------------

:meth:`Solution.condense <orpheus.sn.solution.Solution.condense>` is the
SN-layer orchestration of energy condensation. It condenses **each
material with its own representative spectrum** — the flux·volume-weighted
flux over the cells where the material appears:

.. math::
   :label: energy-condensation-representative-spectrum

   \varphi^{(m)}_g \;=\;
   \sum_{i:\,\mathrm{mat}(i)=m} V_i\,\phi_{i,g},

.. (vv-status rationale) Representational identity: the per-material
   representative spectrum used as the condense test weight — the
   flux·volume-weighted flux over the material's cells (mirrors how
   ``homogenize`` derives its flux weight). A definition consumed by
   :meth:`Mixture.condense`; the end-to-end rate preservation it feeds is
   the L1 gate, not a separate claim.
.. vv-status: energy-condensation-representative-spectrum documented

used as the test weight in :meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>` — the data-layer
collapse verb, whose spectrum-weighted-collapse theory is
:ref:`sn-energy-condensation` (:doc:`/theory/foundations/frame`) — mirroring how
:meth:`Solution.homogenize` derives its flux weight from the same solved
flux. The result is a **portable** ``dict[int, Mixture]`` keyed by
material id — few-group cross sections carrying the coarse ``eg``, not
bound to any mesh (the **mesh-decoupled** half of the asymmetry law,
:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`). A material with
no flux in a fine group contributes zero weight there; the condense
frame's Moore–Penrose Gram handles any empty coarse group.


.. _sn-gotchas:

Gotchas
=======

Each gotcha is a **consequence → how it manifests → which test / level
catches it** — a trap that hides a solver bug behind a green test.  They
should *shrink* over time as the code hardens.

Degeneracy traps — a passing test that proves nothing
-----------------------------------------------------

.. admonition:: Gotcha — 1-group eigenvalue tests are degenerate
   :class: warning

   :math:`k = \nu\Sigma_f / \Sigma_a` is **flux-shape independent**: a
   1-group eigenvalue is a material-property ratio computable *without*
   solving the transport equation, so it cannot detect any error in the
   spatial, angular, or scattering operators.  **Any verification claim
   needs** :math:`\geq 2` **groups** (``vv-principles`` anti-pattern #3).  A
   1-group eigenvalue is still fine for a *rate* or *convergence-order*
   claim — declare the claim layer.

.. _sn-homogeneous-degeneracy-gotcha:

.. admonition:: Gotcha — homogeneous / uniform-rescale invariance hides coefficient bugs
   :class: warning

   Any eigenvalue problem whose target is the flux :math:`\phi` is invariant
   under a uniform rescale :math:`\phi \to C\phi` (the factor :math:`C`
   cancels in the Rayleigh quotient
   :math:`k = \nu\Sigma_f\,\phi / \Sigma_a\,\phi`).  Homogeneous and
   same-material multi-region problems have a **spatially-uniform** rescale,
   so they are blind to factor-of-two coefficient errors that preserve the
   flux shape — and, in curvilinear geometry, blind to redistribution bugs
   (flat flux → the :math:`\alpha` terms vanish identically).  Only a
   genuine **material interface** makes the rescale factor :math:`C(x)`
   position-dependent and breaks the cancellation.

   This is exactly how **ERR-025** hid.  A missing :math:`1/W` normalisation
   in the 1-D diamond-difference recurrence halved the per-ordinate flux, but
   for Gauss–Legendre :math:`W = \sum_n w_n = 2` the missing factor rescaled
   it back — so every homogeneous test passed at machine precision while the
   heterogeneous eigenvalue was :math:`\sim 1.5\,\%` wrong and did **not**
   converge away under mesh or :math:`S_N`-order refinement (the gap
   plateaued in angle).

   **Catcher:** at least one *absolute*-:math:`\phi` test — the fixed-source
   flat-flux diagnostic (:math:`Q/\Sigma_t`), or an absolute eigenvalue
   comparison against a structurally-independent heterogeneous reference.
   The live pins are the L0 symbolic-recurrence check
   :func:`tests.sn.sweep.slab.test_dd_recurrence.test_dd_per_cell_recurrence_matches_symbolic_derivation`
   and the L1 heterogeneous absolute-:math:`k` regression
   :func:`tests.sn.eigenvalue.test_keff_slab.test_heterogeneous_absolute_keff`
   (a 2-region A+B reflective slab pinned against the Case
   singular-eigenfunction reference; the pre-fix :math:`1.48\times10^{-2}`
   error fails it by two orders of magnitude).

.. admonition:: Gotcha — conservation holds even with wrong per-ordinate balance
   :class: warning

   Global particle balance **telescopes** by construction, so a *scalar*
   balance sum can hold to machine precision while the *per-ordinate*
   flat-flux residual is wrong (``vv-principles`` anti-pattern #8; the
   identity :math:`\sum_n w_n(\alpha_{n+1/2} - \alpha_{n-1/2}) = 0`
   annihilates per-ordinate redistribution errors that cancel in the sum).
   The load-bearing invariant is the **per-ordinate** flat-flux residual
   (= 0), not the telescoped scalar balance.

Curvilinear redistribution
--------------------------

.. admonition:: Gotcha — the α redistribution dome must stay non-negative
   :class: warning

   The curvilinear :math:`\alpha` dome must be non-negative; a negative entry
   drives NaN / overflow through the angular sweep.  The fixed-source
   flat-flux diagnostic (:math:`Q/\Sigma_t`) is the single most powerful
   curvilinear bug detector — a spike at :math:`r = 0` localises a missing
   :math:`\Delta A / w` geometry factor.

Solver-coordination traps
-------------------------

* **Renormalise-then-report ordering.**
  :func:`~orpheus.numerics.eigenvalue.power_iteration` renormalises
  :math:`\phi` to unit production **between**
  :meth:`~orpheus.sn.solver.SNSolver.solve_fixed_source` and
  :meth:`~orpheus.sn.solver.SNSolver.compute_keff`.  So
  ``compute_keff`` sees the *renormalised* :math:`\phi`, while the stored
  ``_psi_typed.boundary`` is the *un-renormalised* trace — the scale
  bridge above is what makes the leakage term consistent across that
  boundary.  Reordering the two (report before renormalise) would break
  the bridge's ``1.0`` shortcut.
* **The outer iterate must stay a bare** ``np.ndarray``.  The Mode-11
  live-arm sentinel in
  ``tests/sn/operators/test_fission_kernel_crosscheck.py`` proves that
  ``power_iteration`` feeds a **bare** :class:`numpy.ndarray` flux to
  :meth:`~orpheus.sn.solver.SNSolver.compute_fission_source`, so the
  bare-``np.ndarray`` dispatch arm of
  :meth:`FissionOperator.apply
  <orpheus.transport.operators.fission.FissionOperator>` is the *live
  production arm* (the sentinel wraps that registered leaf in-process and
  asserts the counter advances).  The estimator's
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  evaluations read the same bare array.  Routing the outer iterate
  through a typed carrier would dark the arm (redding the sentinel) and
  break the estimator's evaluate path — the bare-array outer iterate is
  a load-bearing contract, not an implementation accident.

.. seealso::

   **Sweep-machinery gotchas** — the Krylov ``restart`` sizing bug
   (ERR-053 family), the product-cylinder edge-extrapolation data-flow
   invariant, and the Mode-12 / ERR-067 :math:`G`-reciprocity metric catch
   — are documented alongside the sweep at :ref:`sn-282-gotchas`.


References
==========

.. _bib-bailey-morel-chang-2010:

.. [BaileyMorelChang2010] T.S. Bailey, J.E. Morel, and J.H. Chang,
   "The asymptotic diffusion-limit accuracy of S\ :sub:`N` angular
   differencing schemes," *Nuclear Science and Engineering*,
   165(2):149--169, 2010 (LLNL preprint LLNL-JRNL-420356; OA at
   https://www.osti.gov/servlets/purl/1020346).  **Eq. 43** gives the
   Morel--Montry weight :math:`\tau_n = (\mu_n - \mu_{n-1/2}) /
   (\mu_{n+1/2} - \mu_{n-1/2})` as the unique weight exact for a flux
   linear in :math:`\mu`, with admissible range :math:`\tau \in [0,1]`
   — the W1 source for dropping the spherical clamp
   (:ref:`sn-tau-clamp-vindication`).  The paper's diffusion-limit
   analysis keeps :math:`r` continuous (it deliberately removes spatial
   differencing to isolate the angular error), so it cannot speak to the
   spatial pole-cell order of :ref:`sn-pole-cell-spatial-closure`.

.. [WuXieFischer1999] G. Wu, Z. Wu, and B. Fischer,
   "A discrete ordinates nodal method for one-dimensional neutron
   transport calculation in curvilinear geometries,"
   *Nuclear Science and Engineering*, 133, 1999
   (DOI 10.13182/NSE99-A2095).  Green's-function nodal S\ :sub:`N` with
   Legendre spatial expansion — "very high precision on coarse spatial
   meshes relative to the standard fine-mesh S\ :sub:`N` method with the
   spatial diamond-differencing scheme."  The documented route to
   better-than-diamond central-cell *spatial* order; cited at
   :ref:`sn-pole-cell-spatial-closure` as the class of fix that would
   lift the inherent pole-cell :math:`\mathcal{O}(h)`.

.. [MorelMontry1984] J.E. Morel and G.R. Montry,
   "Analysis and elimination of the discrete-ordinates flux dip,"
   *Transport Theory and Statistical Physics*, 13:5, 1984.

.. [LewisMiller1984] E.E. Lewis and W.F. Miller, Jr.,
   *Computational Methods of Neutron Transport*,
   John Wiley & Sons, 1984.

.. [CaseZweifel1967] K.M. Case and P.F. Zweifel,
   *Linear Transport Theory*,
   Addison-Wesley, 1967.

.. [Lebedev1999] V.I. Lebedev and D.N. Laikov,
   "A quadrature formula for the sphere of the 131st algebraic order
   of accuracy," *Doklady Mathematics*, 59(3):477--481, 1999.

.. [CarlsonLathrop1965] B.G. Carlson and K.D. Lathrop,
   "Transport theory -- the method of discrete ordinates,"
   in *Computing Methods in Reactor Physics*,
   Gordon and Breach, 1968.

.. [AdamsLarsen2002] M.L. Adams and E.W. Larsen,
   "Fast iterative methods for discrete-ordinates particle transport
   calculations," *Progress in Nuclear Energy*, 40(1):3--159, 2002.
   Reviews the SAILOR / Larsen-Adams preconditioned-Krylov framework
   that the Wave E Round 2 inner solver implements.  §II gives the
   source-iteration spectral radius :math:`\rho = c`; §VI reviews the
   KBA / wavefront parallel-sweep ordering whose fan-in discipline the
   octant-group Gauss-Seidel schedule inherits.

.. [TrefethenBau1997] L.N. Trefethen and D. Bau, III,
   *Numerical Linear Algebra*, SIAM, 1997.  §27 (power-iteration
   analysis; the dominance-ratio convergence bound); §3.2 (the
   matrix-free Krylov view).

.. [Polizzi2009] E. Polizzi,
   "Density-matrix-based algorithm for solving eigenvalue
   problems," *Physical Review B*, 79, 115112, 2009.  The FEAST
   contour-integral algorithm — the ``eigenvalue_method``
   forward-hook target.

.. [Pautz2002] S.D. Pautz,
   "An algorithm for parallel S\ :sub:`n` sweeps on unstructured
   meshes," *Nuclear Science and Engineering*, 140(2):111--136, 2002.
   The KBA-style wavefront octant-scheduling reference.  Cited at
   the diagonal-cubature shared-face rule (ERR-056), under the
   boundary Gauss-Seidel schedule above: a boundary face outflowed
   by several octants must be reduced (reflected) only after the
   LAST contributing octant group has swept it.

.. [Pomraning1989] G.C. Pomraning,
   "The transport equation in general geometry,"
   *Nuclear Science and Engineering*, 101:330--340, 1989.
   Page 339 frames the curvilinear pole singularity as **structural**:
   :math:`r = 0` is intrinsically singular in any curvilinear streaming
   operator because the angular-derivative coefficients contain
   :math:`1/r`.  Cited at
   :ref:`sn-phase-d-pomraning-structural-singularity` as the deeper
   reason :math:`\mu = \pm 1` is the only admissible Carlson starting
   direction.

.. [Rahnema2008] F. Rahnema, S. Douglass, and B. Forget,
   "Generalized Energy Condensation Theory,"
   *Nuclear Science and Engineering*, 160:41--58, 2008.
   DOI `10.13182/NSE160-41 <https://doi.org/10.13182/NSE160-41>`_.
   Expands the within-coarse-group flux in a set of orthogonal functions;
   the **zeroth moment** (the piecewise-constant basis function on the
   coarse group) recovers the standard flux-weighted multigroup average
   exactly, and the higher moments add the within-group spectral detail.
   Cited at :ref:`sn-condensation-petrov-galerkin-frame` as the rank-0
   precedent for the spectrum-weighted collapse (rank-:math:`>0` faithful
   reconstruction is deferred, `Issue #275
   <https://github.com/deOliveira-R/ORPHEUS/issues/275>`_).

.. [WIMSD] International Atomic Energy Agency,
   *WIMS-D Library Update*, IAEA-TECDOC / STI/Pub/1264, IAEA, Vienna,
   2007. Tables 11.1 (69-group) and 11.2 (172-group) energy-group
   structures, and Table 11.3 (the 172→69 correspondence). The coarse
   group structures ORPHEUS condenses onto; Table 11.3 is the
   derivation-validation oracle for the containing-interval partition.
   Cited at :ref:`sn-condensation-fractional-overlap`.

.. _sn-chapters:

Chapters in this sub-book
=========================

This page is the S\ :sub:`N` sub-book's index. It currently carries the
bulk of the method's theory inline; the decomposition into chapters is
tracked as issue `#231 <https://github.com/deOliveira-R/ORPHEUS/issues/231>`_.
Chapters split out so far:

.. toctree::
   :maxdepth: 2

   slab_one_group
   slab_multigroup
   cartesian_multid
   curvilinear_one_group
   angular_quadrature
   loss_representation
   boundary_conditions
   verification
   history
