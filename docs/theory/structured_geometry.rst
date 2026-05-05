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
