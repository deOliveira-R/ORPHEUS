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
      key_types: [AngularFlux, SNMesh, HarmonicMomentFlux, SweepDependencyGraph]
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
:term:`scalar flux` (contrast the collision-probability integral form).  It resolves
streaming, anisotropic scattering, and interface angular current directly.
ORPHEUS supports three coordinate systems under one balance framework:
**Cartesian** (slab / 2-D, no inter-ordinate coupling), **spherical** (1-D
radial, a single :math:`\alpha`-redistribution dome coupling all ordinates in
:math:`\mu`), and **cylindrical** (1-D radial, an independent :math:`\alpha`
dome per :math:`\mu`-level).  All three share a geometry factor
:math:`\Delta A / w` that guarantees per-ordinate flat-flux consistency; the
curvilinear formulation follows :cite:`MorelMontry1984` in the
:cite:`BaileyMorelChang2010` Eqs. (42)/(43) form (the Morel–Montry
angular-closure weight — unique exact-on-linear-in-:math:`\mu`),
the general framework :cite:`LewisMiller1984`, and the angular discretisation
:cite:`CaseZweifel1967` / :cite:`Hebert2009` (§3.9.4).

The solver is posed as an **operator algebra** over five operators: streaming
:math:`L` (bulk :math:`\hat{\Omega}\cdot\nabla`), collision / removal
:math:`C`, the scattering gain :math:`S`, the boundary law :math:`B` — a
first-class **sibling** operator, *not* folded into :math:`L` — and the rank-1
fission dyad :math:`F`.  They compose the within-group loss operator
:math:`A = L + C - S - B`, so the eigenvalue problem is
:math:`A\,\psi = \tfrac{1}{k}\,F\,\psi` (fixed source: :math:`A\,\psi = q`).
The sub-composite :math:`(L+C)` is lower-triangular under the upwind cell
ordering, which is exactly why :math:`(L+C)^{-1}` **is** the transport :term:`sweep`
(:doc:`/theory/methods/sn/loss_representation`).  :class:`SNSolver` satisfies
the :class:`~numerics.eigenvalue.EigenvalueSolver` protocol and
:func:`solve_sn` returns a :class:`~orpheus.sn.solution.Solution`.  Because the protocol places the
scattering source *inside* ``solve_fixed_source``, the inner source iteration
(in-scatter + anisotropic convergence) stays encapsulated in the SN sweep,
while the outer :func:`~numerics.eigenvalue.power_iteration` loop is the one
shared by CP, MoC, diffusion, and the homogeneous solver (see
:doc:`/api/numerics` for the protocol contract).

The spatial closure is **diamond difference** with the Morel–Montry weight
(:math:`\beta = 0`); the :term:`per-ordinate <ordinate>` discrete transport is the **sweep**,
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
:ref:`sn-direct-seed-solve`.  The "#229 floor" resolution
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
   with an angular :term:`quadrature`, precomputing the coordinate-specific
   streaming stencil.  Its **primary representation is the per-axis
   tuple** :attr:`SNMesh.axes <orpheus.sn.mesh.augmented_mesh.SNMesh.axes>` (the SN phase space factors as a tensor
   product of per-axis 1-D meshes): a legacy ``Mesh1D`` / ``Mesh2D`` is
   converted to axes **once** at the inbound boundary, and
   :meth:`SNMesh.from_axes` stores the caller's tuple verbatim. After
   C5 (:ref:`sn-axis-primary-c5`) the ``mesh`` attribute is *inbound
   provenance only* — ``None`` for an axis-native :math:`d \ge 3` mesh,
   which carries no legacy mesh at all.  (A literal, not an ``:attr:``
   role: the base ``MaterialMesh`` sets it on the instance, so there is
   no autodoc target to link.)  It also **resolves boundary
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
   - **Spherical**: ``face_areas`` (:math:`4\pi r^2`) and ``delta_A``.
   - **Cylindrical**: ``face_areas`` (:math:`2\pi r`) and ``delta_A``.

   The **angular** factor is not a streaming-factory output either
   (2026-08-26): the :math:`\alpha`-dome and the starting direction
   :math:`\mu_{\rm start}` are produced once, per :math:`\mu`-level, by
   :func:`~orpheus.geometry.reduced_operator.angular_redistribution` and
   carried on
   :class:`~orpheus.geometry.reduced_operator.AngularRedistribution`.
   ⛔ ``redist_dAw`` / ``redist_dAw_per_level`` retired with them: that
   array was the *fused product* :math:`\Delta A_i \otimes 1/w_n` of a
   geometric with a quadrature factor, and each of its two consumers
   wanted a different one of the two — so both now form it from
   ``delta_A`` and the measure's weights.

   The Morel--Montry angular weight :math:`\tau` is **not** a factory
   output: it is owned by the angular closure
   (:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`),
   since the geometry-side producer was retired in Issue #236 Phase 2
   Step C (see :ref:`sn-tau-c-on-cellvisit-live`).  The geometry
   factories carry **geometry only** — face areas, the
   :math:`\alpha`-dome, the redistribution factor :math:`\Delta A/w`,
   and the level :term:`starting-direction <starting direction>` edge :math:`\mu_{1/2}`.

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
   solve_sn() --> Solution

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
term requires :term:`diamond difference` in **both space and angle**.

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

For cylindrical 1-D radial sweeps, ordinates with axial direction cosine
:math:`|\mu_z| \to 1` have radial direction cosine
:math:`|\eta| = \sqrt{1 - \mu_z^2} \to 0`.  This is a property of the
polar level, not of any one rule family: it holds on the admitted
:math:`\sigma_y`-folded product family exactly as it did on the
full-circle product and level-symmetric rules those replaced (refused at
cylindrical ``SNMesh`` admission since Q5.6.3 —
:ref:`sn-direct-seed-r12a`).  In this limit the cell
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


.. implements:: dd-slab-scalar
   :by: orpheus.geometry.reduced_operator.slab_streaming

   **Implemented by** 9 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.cell_balance.cell_balance_terms

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.diamond.DiamondDifference.update

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.diamond._DD_W

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average

.. implements:: dd-slab-scalar
   :by: orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission

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


.. implements:: dd-curvilinear-scalar
   :by: orpheus.transport.spatial.cell_balance.cell_balance_for_streaming

   **Implemented by** 6 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: dd-curvilinear-scalar
   :by: orpheus.transport.spatial.cell_balance.cell_balance_terms

.. implements:: dd-curvilinear-scalar
   :by: orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients

.. implements:: dd-curvilinear-scalar
   :by: orpheus.transport.spatial.diamond.DiamondDifference.residual

.. implements:: dd-curvilinear-scalar
   :by: orpheus.transport.spatial.diamond.DiamondDifference.update

.. implements:: dd-curvilinear-scalar
   :by: orpheus.derivations.discrete.sn.balance.derive_wdd_solve

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
   :label: dd-cylindrical-degenerate

   \mathrm{denom} \;=\; (\Delta A / w)\,c_{\rm out} + \Sigma_t\,V_i,
   \qquad
   \mathrm{numer} \;=\; Q_i\,V_i / W
                       + (\Delta A / w)\,c_{\rm in}\,
                          \psi_{n-\tfrac12,\,i},


.. implements:: dd-cylindrical-degenerate
   :by: orpheus.transport.spatial.cell_balance.cell_balance_for_streaming

   **Implemented by** 3 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: dd-cylindrical-degenerate
   :by: orpheus.transport.spatial.cell_balance.cell_balance_terms

.. implements:: dd-cylindrical-degenerate
   :by: orpheus.transport.spatial.diamond.DiamondDifference.update

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

Of the three planned alternatives, one has landed:
:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`
(:math:`\mathcal{O}(\Delta x^2)`, better robustness in optically-thick
cells) ships today under the registry key ``"linear_discontinuous"``,
with its own MMS spatial-convergence gates
(``tests/sn/verification/mms/test_mms_ld_slab.py`` and
``test_mms_ld_2d.py``); see :ref:`ld-ubld-multidim` for the derivation
and the multi-dimensional wiring.  The other two are still **reserved,
not yet implemented**, and are therefore written as literals rather than
``:class:`` roles — a live role would assert a class that does not
exist: ``Step`` (positivity-preserving,
:math:`\mathcal{O}(\Delta x)`) and ``ExponentialCharacteristic``
(positivity-preserving by construction).  Each lands with its own MMS
spatial-convergence verification.

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
LD stress MMS in :doc:`/theory/verification/sn`; the curvilinear machinery — the
sequential ordinate sweep, the pole angular closure (#168 Phase B),
the sweep-frame apply matvec (Phase C), and the direct :math:`\psi_{1/2}`
starting-direction solve (#282 route (a)) — in
:doc:`curvilinear_one_group` (the group axis rides that machinery as
data — :doc:`curvilinear_multigroup`). The #168 Phases A/D/F, ERR-058,
and #196
campaign record is preserved in :doc:`curvilinear_numerics`; the
section below preserves the dispatch-consolidation record (Wave D
Round 2, superseded by the ``LossRepresentation`` polymorphism).

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
``None``.  Wave D replaced it with a geometry-layer boolean,
``ReducedStreamingOperator.requires_upstream_angular_state`` —
``False`` for slab + 2-D Cartesian (no angular redistribution
between successive half-angles), ``True`` for spherical +
cylindrical.  Two-D Cartesian set ``sn_mesh.reduced is None``
(no curvilinear math needed), and the dispatch fell through
to the Cartesian path.

.. note::

   **Both** of those spellings are now retired, so this section is
   two steps of history rather than one.  The boolean was retired on
   2026-08-26: it was exactly ``coord is not CoordSystem.CARTESIAN``
   and had no production reader, the concept having been respelled by
   ``upstream_state.angular_upstream is None`` (what the DD and LD
   cell bodies branch on) and by
   :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.is_cartesian`.
   Strategy selection today is
   :func:`~orpheus.sn.loss_representation.default_for`, which picks the
   first :data:`~orpheus.sn.loss_representation.LOSS_REPRESENTATIONS`
   entry whose ``supports`` admits the mesh, keyed on ``is_1d`` **and**
   ``is_cartesian`` — neither alone is a sufficient discriminator.

Why this mattered:

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

The cell-update strategy lives on the ``scheme`` attribute that
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` realizes in its
constructor (introduced in this round as a constructor argument with
default
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`).  The
default reproduces the inlined sweep math bit-identically — every
regression snapshot at ``tests/sn/regression/snapshots/`` was
generated with DD and continues to match bit-for-bit when the
unified sweep dispatches via ``scheme.update(...)``.  See
:ref:`cell-update-strategies` for the strategy contract and
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` for the
DD scalar form.

:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`
has since landed on exactly this dispatch: a user selects it today by
passing ``scheme=LinearDiscontinuous()`` at
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` construction
(``augmented_mesh.py`` records the same recipe at the ``self.scheme``
default).  ``Step`` and ``ExponentialCharacteristic`` remain
**reserved, not yet implemented** — literals rather than ``:class:``
roles for that reason — and the unified dispatch infrastructure is in
place to receive them.

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


.. _sn-gotchas:

Gotchas
=======

Each gotcha is a **consequence → how it manifests → which test / level
catches it** — a trap that hides a solver bug behind a green test.  They
should *shrink* over time as the code hardens.

.. _sn-symptom-table:

Where to look first — symptom → chapter
----------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 38 32 30

   * - Symptom
     - First suspect
     - Go to
   * - :math:`k` wrong on a vacuum-bounded problem (overshoot
       :math:`\approx L/A`)
     - the reported-:math:`k` functional omitting leakage (the #291
       class)
     - :ref:`sn-keff-estimator`
   * - :math:`k` wrong only when :math:`(n,2n)` is present
     - the emission mis-posed as production (the R7 fork)
     - :ref:`sn-keff-estimator`
   * - :math:`k` right on 1-group / reflective, wrong on multigroup
       heterogeneous
     - scattering-matrix orientation drift — a 1-group green proves
       nothing (degeneracy trap below)
     - :ref:`scattering-matrix-convention`, :doc:`slab_multigroup`
   * - flux spike at :math:`r = 0` on a curvilinear fixed source
     - a missing :math:`\Delta A/w` geometry factor
     - :ref:`balance-curvilinear`
   * - NaN / overflow marching through the angular sweep
     - a negative :math:`\alpha`-dome entry (warning below)
     - :doc:`curvilinear_one_group`
   * - negative or oscillating flux on coarse Cartesian cells
     - the DD closure's unboundedness — refine or change closure
     - :doc:`/theory/foundations/discretization`,
       :ref:`ld-ubld-multidim`
   * - SI iteration count blows up as :math:`c \to 1`
     - :math:`\rho = c` physics, not a bug (acceleration is the fix)
     - :ref:`si-within-group-splitting`
   * - sweep and matvec disagree
     - forbidden post-#206 (one walk) — a representation-seam
       regression
     - :ref:`loss-rep-one-walk-one-instance`
   * - Krylov stalls or diverges after a composite-sizing refactor
     - ``restart`` / ``n_dof`` sizing (the ERR-053 family)
     - :ref:`sn-direct-seed-gotchas`
   * - MMS recovers a lower order than theory
     - the ansatz nulls the term (Mode 7) or the regime is degenerate
     - :doc:`/theory/verification/sn`
   * - the adjoint reciprocity gate reds
     - the three-transposes landmine (Euclidean / Hilbert / walk)
     - :ref:`sn-adjoint`

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
  bridge (:ref:`sn-keff-estimator`) is what makes the leakage term
  consistent across that boundary.  Reordering the two (report before
  renormalise) would break the bridge's ``1.0`` shortcut.
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
   — are documented alongside the sweep at :ref:`sn-direct-seed-gotchas`.


.. _sn-chapters:

Chapters in this sub-book
=========================

This page is the S\ :sub:`N` sub-book's index: the synopsis, the
architecture map, the transport equation, and the shared cell-update
and dispatch contracts; the chapter decomposition is tracked as issue
`#231 <https://github.com/deOliveira-R/ORPHEUS/issues/231>`_.

Several orders through the book serve different jobs (tracks, not one
sequence):

* **Newcomer** — :doc:`placement` first (*why* discrete ordinates —
  the trade-space against CP/MoC/P\ :sub:`N`/diffusion/MC), then the
  broadening progression in toctree order:
  :doc:`slab_one_group` (the whole machine at its simplest) →
  :doc:`slab_multigroup` (energy and the eigenvalue) →
  :doc:`cartesian_multid` (space) → :doc:`curvilinear_one_group` →
  :doc:`curvilinear_multigroup`, with :doc:`angular_quadrature` and
  :doc:`boundary_conditions` as on-ramp references.
* **Modifying the sweep** — :doc:`loss_representation` (the
  representation catalog and the one-walk theorems,
  :ref:`loss-rep-one-walk-one-instance`) → the strategy contract on
  this page (:ref:`cell-update-strategies`) → the multi-D schedule
  (:ref:`sweep-octant-dependency-graph`) → the curvilinear sequential
  sweep (:doc:`curvilinear_one_group`) → the sweep-machinery gotchas
  (:ref:`sn-direct-seed-gotchas`).
* **Porting an equation from the literature** — the machine header at
  the top of this page (sign / normalization / layout / ordering
  conventions) → :doc:`angular_quadrature` (weight normalization) →
  the scattering-matrix convention
  (:ref:`scattering-matrix-convention`) and :ref:`pn-scattering` →
  :doc:`/theory/verification/sn` for the gate the new equation ships with.
* **Debugging a wrong answer** — start at the symptom table
  (:ref:`sn-symptom-table`), then the gotcha catalog below it.

The chapters:

.. toctree::
   :maxdepth: 2

   placement
   slab_one_group
   slab_multigroup
   cartesian_multid
   curvilinear_one_group
   curvilinear_multigroup
   curvilinear_numerics
   angular_quadrature
   loss_representation
   boundary_conditions
   solver
   acceleration
   adjoint
   history
