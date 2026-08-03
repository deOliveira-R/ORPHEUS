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
* :class:`StructuredGeometry` is **pure shape**: a geometry tag
  (``"SLB"`` / ``"CYL"`` / ``"SPH"``), an ordered tuple of
  :class:`~orpheus.geometry.structured_geometry.Region` (each a
  ``(mat_id, outer_thickness_cm)`` pair), and a tuple of
  :class:`~orpheus.geometry.mesh.BC` instances at the geometry's
  endpoints. **No cell counts, no critical-dimension scalars, no
  energy-group count, no infinite-medium kind.**
* The geometry → mesh transition is **the single explicit point**
  where discretization information enters the pipeline.
  :meth:`Mesh1D.from_geometry(geom, region_meshes=...) <orpheus.geometry.mesh.Mesh1D.from_geometry>`
  takes a tuple of
  :class:`~orpheus.geometry.mesh.RegionMesh` (one per region;
  ``n_cells`` + ``method`` ∈ {``"equal-volume"``, ``"uniform"``}) and
  emits a discretized :class:`~orpheus.geometry.mesh.Mesh1D`.
* Reference solvers (``Billiard``, ``MomentSpace``, ``Spectrum``,
  ``BasisSpace``) take ``(geometry: StructuredGeometry, materials,
  **method_kwargs)`` directly via ``__init__``. They never see a
  mesh. They never see ``n_cells``.
* Discrete production solvers take ``(materials, mesh, params)``
  where ``mesh`` is built via ``Mesh1D.from_geometry``.
* Slab convention: :attr:`StructuredGeometry.domain_extent_cm` is
  the **full slab width** (sum of all region thicknesses,
  end-to-end). F_N's natural half-thickness ``a = L / 2`` is
  recovered inside :class:`MomentSpace`.
* The Sood case registry adapter is
  :meth:`La13511Case.to_geometry()
  <orpheus.derivations.continuous.sood_registry.la13511.La13511Case.to_geometry>`,
  which materialises a :class:`StructuredGeometry` from the case's
  ``geometry_kind`` tag and ``truth.critical_dimension_mfp`` (cm = mfp / Σ_t).
  Infinite-medium cases raise — for ``k_\infty`` use
  :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
  or :meth:`MomentSpace.solve_kinf`.
* Two non-trivial classmethods earn their keep on
  :class:`StructuredGeometry`:
  :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
  (CYL with the ``r_cell = pitch / √π`` equal-area transformation)
  and
  :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.pwr_slab_half_cell`
  (SLB with reflective fuel-centre symmetry plane).


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
geometry kind and the radius in cm, it does not need a cell count
and it does not need the published critical dimension's name.

Phase F separates the three concerns into three layers:

1. **Geometry layer** —
   :class:`~orpheus.geometry.structured_geometry.StructuredGeometry`,
   :class:`~orpheus.geometry.structured_geometry.Region`. Pure
   shape + BCs. No cell counts. No scalars from a published table.
2. **Mesh layer** —
   :class:`~orpheus.geometry.mesh.Mesh1D`,
   :class:`~orpheus.geometry.mesh.RegionMesh`. Discrete
   representation. Discretization is supplied at this layer's
   construction step, not pinned to the geometry.
3. **Registry layer** —
   :class:`~orpheus.derivations.continuous.sood_registry.la13511.La13511Case`,
   :class:`~orpheus.derivations.continuous.sood_registry.la13511.La13511Truth`.
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


Geometry-tag semantics
======================

The supported tags are uppercase three-letter mnemonics:

==========  ===============  ==================  ==============
Tag         Coordinate       Endpoints (BC)      Centreline
==========  ===============  ==================  ==============
``"SLB"``   Cartesian        2 (left, right)     n/a
``"CYL"``   Cylindrical      1 (outer)           implicit reflective
``"SPH"``   Spherical        1 (outer)           implicit reflective
==========  ===============  ==================  ==============

The endpoint count comes from the orbit-space classification of the
geometry's billiard table: SLB has two flat surfaces (orbit-space
rank 2); CYL and SPH have one outer surface plus an implicit
centreline reflection (orbit-space rank 1). This is the same
classification the trajectory-resolvent ``Billiard`` class uses.

When future geometries land that genuinely have two surfaces
(``HSPH`` for hollow sphere, ``ANN`` for annulus), they extend
this map with two-endpoint BC tuples — the central abstraction
generalises cleanly.


End-state spot checks
=====================

These are the canonical call sites the reset enables.

Production user — no registry, no truth, no critical anything::

    geom = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=3.6,
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=10),
        RegionMesh(n_cells=3),
        RegionMesh(n_cells=7),
    ))
    materials: Materials = {2: uo2_fuel(), 1: zircaloy(), 0: borated_water()}
    result = solve_cp(materials, mesh, CPParams())
    print(result.keff)

Registry consumer — Sood case::

    case = LA13511_CASES["Ua-1-0-SP"]
    geom = case.to_geometry()
    mesh = Mesh1D.from_geometry(geom, region_meshes=(RegionMesh(n_cells=64),))
    result = solve_cp(case.materials, mesh, CPParams())

Reference solver direct — NO mesh, NO RegionMesh::

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
* Locked decision 5 (``Region`` is geometry-only) — no ``n_cells``
  on a region; that lives on :class:`RegionMesh` at the mesh layer.
* Locked decision 6 (``Mesh1D.from_geometry``) — the single explicit
  point where discretization enters.

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
SN, MoC, and CP curvilinear :term:`sweeps <sweep>` all march through the same data:

* **chord lengths**: cell radial widths
  (:attr:`~orpheus.geometry.mesh.Mesh1D.widths`),
* **face areas**: :math:`A_{i+1/2} = 4\pi r_{i+1/2}^2` (sphere) or
  :math:`2\pi r_{i+1/2}` (cylinder),
* **the geometry factor** :math:`\Delta A_i / w_n` that ensures
  :term:`per-ordinate <ordinate>` flat-flux consistency,
* **the Bailey 2009 dome recursion** for :math:`\alpha`, and
* **the Morel--Montry angular closure** :math:`\tau_{mm}`
  (unclamped raw weight on the sphere; :math:`[1/2, 1]`-clamped on the
  cylinder — see :eq:`morel-montry-clamp` below).

Per Cardinal Rule 2 (architecture is critical), the same data **MUST
NOT** be duplicated across solvers.
:class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
in :mod:`orpheus.geometry.reduced_operator` lifts the math into a
geometry-layer primitive that downstream consumers (SN, MoC, CP)
share.

The Bailey 2009 dome recursion (sphere):

.. math::
   :label: bailey-dome-recursion

   \alpha_{n+\tfrac12} = \alpha_{n-\tfrac12} - w_n\,\mu_n,
   \qquad \alpha_{\tfrac12} = 0.

.. vv-status: bailey-dome-recursion documented

For Gauss--Legendre :term:`quadrature` with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome
(:math:`\alpha_{1/2} = 0 \to \text{peak} \to \alpha_{N+1/2} = 0`)
that closes back to zero at the upper boundary by GL antisymmetry.
The cylindrical analog runs **per-:math:`\mu`-level**: each level
:math:`p` carries its own :math:`(M+1)`-tuple of
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine and :math:`M` is the
number of azimuthal ordinates on that level.

The Morel--Montry closure (Bailey-Morel-Chang 2010 Eq. 43) is the raw
*fractional position* of the ordinate :math:`\mu_m` in its half-angle
interval; the production code applies a **geometry-dependent** clamp on
top of it:

.. math::
   :label: morel-montry-clamp

   \tau_m^{\rm raw} = \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}},
   \qquad
   \tau_m = \begin{cases}
     \tau_m^{\rm raw} & \text{sphere (unclamped, since W1)} \\[4pt]
     \mathrm{clip}\!\left(\tau_m^{\rm raw},\;\tfrac12,\;1\right)
       & \text{cylinder}
   \end{cases}

.. vv-status: morel-montry-clamp documented
   Rationale: this is the literature-transcribed definition of the M-M
   angular closure weight (Bailey-Morel-Chang 2010 Eq. 43) plus the
   production clamp policy; it is a representational identity, not a
   solver claim. The verifiable content is the producer-equivalence
   gate ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py``
   (τ is closure-owned, NOT a reduced-operator field — see the
   τ-ownership note below) and the W1 clamp-silence + positivity gates
   named at :ref:`sn-curvilinear-aniso-norm-reconciliation`.

.. _tau-ownership-note:

.. note:: **Where τ lives.**

   :math:`\tau_m` is **not** a
   :class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
   field.  The geometry-layer primitive carries the curvature
   coefficients (``face_areas``, ``delta_A``, ``alpha_half``,
   ``redist_dAw``, ``mu_start`` and their ``*_per_level`` cylindrical
   siblings) and nothing else; the M-M closure weight is produced by
   :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_raw_per_level`
   (raw) and
   :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
   (the clamp policy applied), which :class:`SNMesh` calls against the
   quadrature and ``self.reduced.coord``.  The split is deliberate: τ
   is a property of the *angular closure scheme*, selectable per mesh,
   while the curvature coefficients are a property of the *geometry*.
   Statements elsewhere in the corpus describing ``tau_mm`` /
   ``tau_mm_per_level`` as factory outputs predate that move.

The raw weight :math:`\tau_m^{\rm raw}` is the **unique** weight exact
for a flux linear in :math:`\mu` (Bailey-Morel-Chang 2010 Eq. 43), with
admissible range :math:`\tau \in [0, 1]`.

* **SPHERE** uses :math:`\tau_m^{\rm raw}` directly (since W1,
  2026-06-13).  On Gauss--Legendre quadrature
  :math:`\tau_m^{\rm raw} \in [0.39, 0.61]` — always interior to
  :math:`[0,1]` — so the closure stays positive without a clamp.  The
  former :math:`[1/2, 1]` clamp was an over-conservative,
  mis-cited (to Lewis & Miller §4.5) positivity floor that was 100 %
  spurious on physical fields and re-floored the anisotropic solution.
  The full vindication is at
  :ref:`sn-curvilinear-aniso-norm-reconciliation` in
  :doc:`/theory/verification/sn`.
* **CYLINDER** retains the clamp :math:`\tau_m =
  \mathrm{clip}(\tau_m^{\rm raw}, \tfrac12, 1)`: product / level-
  symmetric quadratures put the most-inward azimuthal ordinate exactly
  on :math:`\eta = -\sin\theta`, giving :math:`\tau_m^{\rm raw} = 0`
  exactly (a structural :math:`\div 0` block, tracked by
  `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_).

For spherical geometry the cell-edge direction cosines
:math:`\mu_{m+1/2}` are the partial weight sums
:math:`\sum_{n \le m} w_n` shifted by :math:`\mu_{1/2} = -1`; for
cylindrical geometry the weights live in :math:`\varphi`-space rather
than :math:`\eta`-space, so the cell edges are taken at the midpoints
of consecutive :math:`\eta` values with endpoints at :math:`\pm\sin\theta`.

API surface
-----------

The geometry-layer primitive is built by three factory functions, one
per coordinate system:

* :func:`~orpheus.geometry.reduced_operator.slab_streaming(mesh, ang)
  <orpheus.geometry.reduced_operator.slab_streaming>` — Cartesian 1-D;
  no curvature math.  ``requires_upstream_angular_state = False``,
  ``angular_marching_axis = None``.  All ``alpha_*`` and
  ``redist_dAw*`` arrays remain ``None``.
* :func:`~orpheus.geometry.reduced_operator.spherical_streaming(mesh, ang)
  <orpheus.geometry.reduced_operator.spherical_streaming>` — 1-D spherical
  with the dome recursion :eq:`bailey-dome-recursion` and Morel--Montry
  closure :eq:`morel-montry-clamp`.  ``requires_upstream_angular_state =
  True``, ``angular_marching_axis = "mu"``.
* :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming(mesh, ang)
  <orpheus.geometry.reduced_operator.cylindrical_streaming>` — 1-D
  cylindrical with **per-:math:`\mu`-level** :math:`\alpha`,
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
:meth:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.streaming_terms`,
which returns a
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` namedtuple
whose populated fields are geometry-dependent (slab is minimal;
sphere/cylinder carry the full curvature-coefficient bundle).

The two trailing fields ``volume`` and ``abs_mu`` carry the per-cell
volume :math:`V_i` and the absolute primary direction cosine
:math:`|\mu|` (sphere) / :math:`|\eta|` (cylinder, radial) /
:math:`|\mu_x|` (slab).  They are populated by all three factories so
a downstream sweep cell update — see :doc:`/theory/methods/sn/index`,
"Cell update strategies (the strategy contract)" — receives a
self-contained per-cell, per-direction packet and need not reach back
into ``SNMesh`` or the ``Quadrature``.  The ``alpha_in is
None`` test discriminates slab from curvilinear inside cell-update
strategies; the cylindrical pure-azimuthal degenerate case
(``abs_mu < 1e-15``) is the single runtime branch a strategy must
handle for cylindrical sweeps.

Geometric labels, not flow-direction labels
-------------------------------------------

The two face-area fields on
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` are
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
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk`, which encodes
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
``SNMesh._setup_spherical`` and ``SNMesh._setup_cylindrical``.  Hash
equality — :func:`numpy.array_equal`, never ``np.allclose`` — was
enforced at test time by ``tests/geometry/test_reduced_operator.py``
(``foundation``-tagged), so the two paths had to share every
floating-point bit.  That is what made the lift safe at the time: the
then-consumers (the SN sweep and the curvilinear Krylov operator) were
unaffected because the two paths computed the same data.

.. warning:: **That gate has since been DEMOTED — read its green
   accordingly.**

   The two legacy setup methods no longer exist (see
   :ref:`snmesh-as-router` below).  ``SNMesh.__init__`` now calls
   :func:`~orpheus.geometry.reduced_operator.spherical_streaming` /
   :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming`
   itself, so the surviving hash-equality legs compare a fresh factory
   call against ``sn_mesh.reduced`` — *the value that same factory
   produced*, routed through the mesh constructor — and the two
   ``SNMesh.face_areas`` / ``SNMesh.delta_A`` legs are deprecated
   read-throughs to that same object.  The gate therefore pins the
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
   (``tests/sn/regression/test_dd_regression.py``) corroborating but
   nowhere the sole evidence:

   * ``delta_A`` — the closed-form L0 term check
     ``TestL0TermVerification::test_delta_A_magnitude`` in
     ``tests/sn/primitives/test_quadrature.py``, against
     :math:`4\pi\,\Delta(r^2)` / :math:`2\pi\,\Delta r`.  **Sole
     catcher**; the snapshots are blind here, and correctly so —
     ``delta_A`` has no production consumer.
   * ``alpha_half`` — the L0 per-ordinate flat-flux identity
     ``test_per_ordinate_flat_flux_consistency[SPHERICAL]``
     (``catches("ERR-006", "ERR-007")``), plus the sphere snapshots.
   * ``alpha_per_level`` —
     ``tests/sn/sweep/curvilinear/test_alpha_closed_form.py``
     (the Dirichlet-kernel closed form; **cylindrical α only** — every
     fixture there is ``CoordSystem.CYLINDRICAL``), plus the
     cylindrical flat-flux arm and the cylinder snapshots.
   * ``redist_dAw`` / ``redist_dAw_per_level`` —
     ``tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py``,
     the L0 closed-form :math:`\varphi = Q/(\Sigma_t(1-c))` identity,
     plus both snapshot families.  The flat-flux identity does **not**
     cover these: it recomputes :math:`\Delta A / w` rather than
     reading the production array.
   * ``face_areas`` — ``tests/geometry/test_geometry.py`` pins the
     producer :func:`~orpheus.geometry.coord.compute_areas_1d`
     against its closed form; the snapshots pin the forwarding.

   ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` is
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

SNMesh as router
----------------

After Round 1.1 of Wave D of the SN reshape campaign, :class:`SNMesh`
**routes** to :class:`ReducedStreamingOperator` rather than computing
the connection coefficients itself.  The :meth:`SNMesh.__init__`
ladder calls :func:`slab_streaming` / :func:`spherical_streaming` /
:func:`cylindrical_streaming` directly; the historical
``SNMesh._setup_spherical`` and ``SNMesh._setup_cylindrical`` methods
no longer exist.  ``self.reduced`` is the new canonical accessor
every downstream consumer should bind to::

    sn_mesh.reduced.streaming_terms(cell_idx, dir_idx, mu_level_idx)

returns the per-(cell, direction) packet a sweep cell update needs —
no more reaching into ``SNMesh`` for a half-dozen separate arrays.

That migration has since **completed**.  Of the eight legacy attribute
names the lift originally re-exposed as :class:`DeprecationWarning`
``@property`` accessors, only **two** survive on :class:`SNMesh` today
— ``face_areas`` and ``delta_A``, still read-throughs to the matching
field on ``self.reduced``.  The other six (``alpha_half``,
``redist_dAw``, ``alpha_per_level``, ``redist_dAw_per_level``,
``tau_mm``, ``tau_mm_per_level``) are gone: consumers bind to
``streaming_terms(...)`` or to ``sn_mesh.reduced.*`` directly, and the
two ``tau_mm`` names have no ``self.reduced`` field left to route to at
all — τ is closure-owned now, not a factory output (see the
:ref:`τ-ownership note <tau-ownership-note>` above).

The Cartesian path is unchanged: ``SNMesh._setup_cartesian`` still
populates the :math:`2|\mu|/\Delta x` and :math:`2|\mu_y|/\Delta y`
streaming stencils used by the DD-denominator precomputation in the
Cartesian sweep (these are SN-specific and not represented in
:class:`ReducedStreamingOperator`).  Slab geometry additionally gets
a slab :class:`ReducedStreamingOperator` for completeness so
``sn_mesh.reduced`` is always populated.

Migration roadmap
-----------------

This primitive is the foundation for several follow-on issues in the
SN reshape campaign (``.claude/plans/sn_reshape.md``):

* **Issue 10 (Wave D Round 1.1) — DONE**: :class:`SNMesh` consumes
  :class:`ReducedStreamingOperator` via the dispatch ladder above.
  The connection-coefficient math no longer lives in :class:`SNMesh`.
* **SN operator algebra (Depth B, 2026-05)** —
  :class:`~orpheus.sn.operators.streaming.StreamingOperator` /
  :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` consume the
  primitive (as ``SNMesh.reduced``) through the loss-representation walk:
  :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply` reads the
  connection coefficients off ``self.mesh.reduced.coord`` inside the walk.
  (Depth B consumed them through the per-geometry
  ``transport_operator_matvec_*`` matvecs; that family and its unified
  successor were deleted in the typed-field (#197) and walk-unification
  (#280 campaigns) refactors — the primitive itself is unchanged.)
* **MoC and CP campaigns (post-Wave-1)** reuse the same primitive
  with their own consumption patterns (track-segment chord march
  for MoC; ray-traced chord-length integrals for CP).
