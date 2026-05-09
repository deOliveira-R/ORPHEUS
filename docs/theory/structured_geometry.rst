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
  :func:`~orpheus.derivations.common.eigenvalue.solve_homogeneous_infinite`
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
SN, MoC, and CP curvilinear sweeps all march through the same data:

* **chord lengths**: cell radial widths
  (:attr:`~orpheus.geometry.mesh.Mesh1D.widths`),
* **face areas**: :math:`A_{i+1/2} = 4\pi r_{i+1/2}^2` (sphere) or
  :math:`2\pi r_{i+1/2}` (cylinder),
* **the geometry factor** :math:`\Delta A_i / w_n` that ensures
  per-ordinate flat-flux consistency,
* **the Bailey 2009 dome recursion** for :math:`\alpha`, and
* **the Morel--Montry angular closure** :math:`\tau_{mm}` clamped to
  :math:`[1/2, 1]`.

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

For Gauss--Legendre quadrature with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome
(:math:`\alpha_{1/2} = 0 \to \text{peak} \to \alpha_{N+1/2} = 0`)
that closes back to zero at the upper boundary by GL antisymmetry.
The cylindrical analog runs **per-:math:`\mu`-level**: each level
:math:`p` carries its own :math:`(M+1)`-tuple of
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine and :math:`M` is the
number of azimuthal ordinates on that level.

The Morel--Montry closure (Lewis & Miller 1984, §4.5; Bailey
et al. 2009 Eq. 74):

.. math::
   :label: morel-montry-clamp

   \tau_m = \mathrm{clip}\!\left(
       \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}},
       \;\tfrac12,\; 1
   \right)

.. vv-status: morel-montry-clamp documented

clamps the M-M weighting to :math:`[1/2, 1]` so the closure remains
positive (Lewis & Miller §4.5).  For spherical geometry the cell-edge
direction cosines :math:`\mu_{m+1/2}` are the partial weight sums
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
  ``angular_marching_axis = None``.  All ``alpha_*``, ``redist_dAw``,
  ``tau_mm`` arrays remain ``None``.
* :func:`~orpheus.geometry.reduced_operator.spherical_streaming(mesh, ang)
  <orpheus.geometry.reduced_operator.spherical_streaming>` — 1-D spherical
  with the dome recursion :eq:`bailey-dome-recursion` and Morel--Montry
  closure :eq:`morel-montry-clamp`.  ``requires_upstream_angular_state =
  True``, ``angular_marching_axis = "mu"``.
* :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming(mesh, ang)
  <orpheus.geometry.reduced_operator.cylindrical_streaming>` — 1-D
  cylindrical with **per-:math:`\mu`-level** :math:`\alpha`,
  :math:`\Delta A/w`, and :math:`\tau_{mm}` lists.  Requires the
  angular measure to expose ``level_indices`` (e.g.,
  :class:`~orpheus.sn.quadrature.LevelSymmetricSN`,
  :class:`~orpheus.sn.quadrature.ProductQuadrature`).

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
a downstream sweep cell update — see :doc:`discrete_ordinates`,
"Cell update strategies (the strategy contract)" — receives a
self-contained per-cell, per-direction packet and need not reach back
into ``SNMesh`` or the ``AngularQuadrature``.  The ``alpha_in is
None`` test discriminates slab from curvilinear inside cell-update
strategies; the cylindrical pure-azimuthal degenerate case
(``abs_mu < 1e-15``) is the single runtime branch a strategy must
handle for cylindrical sweeps.

Bit-identical contract
----------------------

The factories produce arrays bit-identical to the historical inline
implementations on
:class:`~orpheus.sn.geometry.SNMesh._setup_spherical` and
:class:`~orpheus.sn.geometry.SNMesh._setup_cylindrical`.  Hash
equality is enforced at test time
(``tests/geometry/test_reduced_operator.py``,
``foundation``-tagged) using :func:`numpy.array_equal` (not
``np.allclose``) — the two paths must share every floating-point bit.
This is what makes the lift safe: today's consumers (the SN sweep
and the BiCGSTAB curvilinear operator) are unaffected because the
two paths compute the same data; tomorrow's consumers can bind to
the geometry-layer primitive instead of duplicating the math.

SNMesh as router (post-Round-1.1 of Wave D)
--------------------------------------------

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

The legacy attribute names (``alpha_half``, ``redist_dAw``, ``tau_mm``,
``alpha_per_level``, ``redist_dAw_per_level``, ``tau_mm_per_level``,
``face_areas``, ``delta_A``) survive as ``@property`` accessors that
emit :class:`DeprecationWarning` and route to the matching attribute
on ``self.reduced``.  This preserves the 6 production read sites in
``orpheus/sn/sweep.py`` and ``orpheus/sn/solver.py`` for the
duration of Wave D Round 2 (Issue 12) and Wave E, which then migrate
those call sites to ``streaming_terms(...)`` directly and remove the
deprecated properties.

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
* **Issues 11/12 (Wave D Round 2/3)** make ``SNStreamingOperator.apply``
  consume the primitive directly, eliminating the SN-specific
  curvature attributes from :class:`SNMesh` (the deprecated
  properties retire here).
* **MoC and CP campaigns (post-Wave-1)** reuse the same primitive
  with their own consumption patterns (track-segment chord march
  for MoC; ray-traced chord-length integrals for CP).
