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
A curvilinear S\ :sub:`N` :term:`sweep <sweep>` marches through this data:

* **chord lengths**: cell radial widths
  (:attr:`~orpheus.geometry.mesh.Mesh1D.widths`),
* **face areas**: :math:`A_{i+1/2} = 4\pi r_{i+1/2}^2` (sphere) or
  :math:`2\pi r_{i+1/2}` (cylinder),
* **the geometry factor** :math:`\Delta A_i / w_n` that ensures
  :term:`per-ordinate <ordinate>` flat-flux consistency,
* **the** :math:`\alpha` **dome recursion**, and
* **the Morel--Montry angular closure** :math:`\tau_{mm}`
  (one geometry-free formula :eq:`morel-montry-closure`, read against a
  per-geometry cell partition :eq:`angular-cell-partition` — see below).

Per Cardinal Rule 2 (architecture is critical), the same data **MUST
NOT** be duplicated across the solvers that need it.
:class:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator`
lifts the math into **one** primitive rather than a per-solver one.

⛔ That sentence read *"in* ``orpheus.geometry.reduced_operator``
*lifts the math into a* **geometry-layer** *primitive rather than a
solver-side one"* until 2026-08-28.  The single-sourcing half stands and
is the Cardinal-Rule-2 point; the **layer** half was refuted by the
un-weld arc's P4.4.  `[M]` the primitive holds no geometry — its one
geometric datum, ``face_areas``, was a verbatim copy of
:attr:`~orpheus.geometry.mesh.Mesh1D.areas`, already single-sourced in
:func:`~orpheus.geometry.coord.compute_areas_1d`, while ``delta_A`` has
zero non-S\ :sub:`N` consumers and
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` carries a
:math:`\Delta A` divided by a *quadrature weight*.  It was also an
island in its own package: every genuine geometry primitive had 1-4
intra-``geometry/`` consumers and this one had **0**.  The layer test it
failed — *a datum belongs to the layer that can define it without naming
a method; everything else is posing, and posing belongs to the method
head that poses it* — moved it to
:mod:`orpheus.sn.mesh.reduced_operator`.

.. important:: **Which solvers is "the solvers that need it"?**  Until
   2026-08-27 this page — and the module's own docstring — read *"SN,
   MoC, and CP curvilinear sweeps all march through the same data"* and
   *"downstream consumers (SN, MoC, CP) share"* it.  That is
   **measurably false, and structurally so**: `[M]` no file under
   ``orpheus/moc/``, ``orpheus/cp/`` or ``orpheus/mc/`` names this
   primitive under **any** of the eight spellings the census at
   :ref:`connection-coefficient-census` enumerates — while both of that
   census's positive controls name every one of them.  That is not a
   migration which has yet to happen, but a term those methods never
   form.  The chart data is still correctly geometry-layer; what was
   wrong is the consumer list.  The reason is worked in
   :ref:`who-needs-a-connection-coefficient` below, and the
   re-runnable predicate — deliberately a predicate rather than a table
   of counts — is at :ref:`connection-coefficient-census`.

The :math:`\alpha` dome recursion (sphere) — Hébert (2009) §3.9.4
Eqs. 3.423-3.424, after Lathrop, K., & Carlson, B. (1966),
*J. Comp. Phys.* 1:173, in the ORPHEUS factor-of-2-absorbed
normalization:

.. math::
   :label: alpha-dome-recursion

   \alpha_{n+\tfrac12} = \alpha_{n-\tfrac12} - w_n\,\mu_n,
   \qquad \alpha_{\tfrac12} = 0.

.. vv-status: alpha-dome-recursion documented
   Rationale: this is the literature-transcribed definition of the
   :math:`\alpha` recursion (Hébert §3.9.4 Eqs. 3.423-3.424, after
   Lathrop & Carlson 1966), a representational identity rather than a
   solver claim.  ⛔ The label read ``bailey-dome-recursion`` until
   2026-08-27; that name encoded the wrong-paper attribution retracted
   at Issue #168 Phase B (see
   :ref:`sn-citation-corrections`).  The verifiable content is the
   dome-closure contract — ``tests/geometry/test_reduced_operator.py``
   (``test_every_shipped_gauss_legendre_dome_closes``,
   ``test_every_shipped_folded_product_dome_closes_on_every_level``,
   and the negative control ``test_a_dome_that_does_not_close_is_refused``)
   plus ``tests/sn/sweep/curvilinear/test_alpha_closed_form.py``
   (``test_production_alpha_is_a_non_negative_closing_dome``).
   ⚠ The SAME recursion is also stated on the S\ :sub:`N` methods page as
   :eq:`alpha-recursion`, which is the label the ``verifies`` markers
   target; see the de-duplication note at the end of this section.

For Gauss--Legendre :term:`quadrature` with :math:`\mu` sorted ascending in
:math:`[-1, 1]`, the recursion produces a non-negative dome
(:math:`\alpha_{1/2} = 0 \to \text{peak} \to \alpha_{N+1/2} = 0`)
that closes back to zero at the upper boundary by GL antisymmetry.
The cylindrical analog runs **per-**\ :math:`\mu`\ **-level**: each level
:math:`p` carries its own :math:`(M+1)`-tuple of
:math:`\alpha^{(p)}_{m+1/2} = \alpha^{(p)}_{m-1/2} - w_m\,\eta_m`,
where :math:`\eta` is the radial direction cosine and :math:`M` is the
number of azimuthal ordinates on that level — **and it closes on that arm
too**, :math:`\alpha^{(p)}_{M+1/2} = 0` on every level, for the same reason
and by the same telescoping.

⭐ **Both ends are zero, and only one of them is an axiom.**  The recursion
is strictly one-sided — it is seeded at :math:`\alpha_{1/2} = 0` and never
consults the far end — so telescoping it over the level gives
:math:`\alpha_{M+1/2} = -\sum_m w_m c_m` in the level's marching cosine
:math:`c` (:math:`\mu` sphere, :math:`\eta` cylinder).  The far endpoint
therefore vanishes **iff the measure's first moment in the marching
coordinate does**, which makes it a *property of the quadrature* — a real
admission contract a bad rule can violate — rather than a property of the
recursion.  One body computes the dome
(:func:`~orpheus.sn.angular.redistribution.alpha_dome`, called by both
curvilinear factories, with the derivations-side name delegating to it) and
one guard refuses a non-closing measure
(``_assert_alpha_dome_closes``, per level on the cylinder so the offending
level is named).  Until ``bea6a367`` (2026-08-12) the contract was a bare
``assert`` on the sphere arm and *nothing* on the cylinder — and a bare
``assert`` is stripped by the canonical ``python -O`` runner, so it did not
run at all.  Full account, including why the fix had to start with
de-duplicating three copies of the recursion:
:ref:`sn-alpha-dome-closes`.

The Morel--Montry closure weight is the **barycentric coordinate of the
ordinate between the two edges of its own angular cell** — predicate
**P2**, :cite:`BaileyMorelChang2010` Eq. 43 = Lathrop 2000 Eq. 23 —
equivalently the UNIQUE closure weight exact for an angular flux affine
in the radial cosine:

.. math::
   :label: morel-montry-closure

   \tau_m
     \;=\; \frac{\mu_m - \mu_{m-1/2}}{\mu_{m+1/2} - \mu_{m-1/2}} .

.. vv-status: morel-montry-closure documented
   Rationale: this is the literature-transcribed definition of the M-M
   angular closure weight (Bailey-Morel-Chang 2010 Eq. 43 = Lathrop 2000
   Eq. 23); it is a representational identity, not a solver claim. The
   verifiable content is the producer-equivalence gate
   ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` (both
   arms now compare against HAND-AUTHORED references — the analytic arc
   closed form on the cylinder, an inline cumulative-weight expression on
   the sphere) plus the ν-closure and P3 gates in
   ``tests/sn/sweep/test_tau_arc_wellposedness.py``.  τ is closure-owned,
   NOT a reduced-operator field — see the τ-ownership note below.

⭐ **There is no "raw" and no "clamped" τ** (Q5.6.4, 2026-08-11): the
:math:`[\tfrac12, 1]` absorber RETIRED, ``morel_montry_tau_raw_per_level``
retired with it, and :eq:`morel-montry-closure` carries **no geometry**.
One generic body serves both arms — the geometry lives entirely in the
**cell partition**, which is a separate object with its own producer.

.. _angular-cell-partition-section:

The angular cell partition — where the geometry actually lives
--------------------------------------------------------------

A level's ordinates each own an angular *cell*; the partition is the
:math:`(M+1)` cell edges in the radial direction cosine
(:math:`\mu` sphere, :math:`\eta` cylinder), and it is the object
:eq:`morel-montry-closure` reads:

.. math::
   :label: angular-cell-partition

   \mu_{m+1/2} \;=\;
   \begin{cases}
     \mu_{m-1/2} + w_m, \qquad \mu_{1/2} = -1
       & \text{sphere: cumulative WEIGHT} \\[8pt]
     \sin\theta\,\cos\omega_{m+1/2},\quad
       \omega_{m+1/2} = \tfrac12\bigl(\omega_m + \omega_{m+1}\bigr),
       \;\; \omega_{1/2} = \pi,\;\; \omega_{M+1/2} = 0
       & \text{cylinder: MIDPOINT in }\omega
   \end{cases}

.. vv-status: angular-cell-partition documented
   Rationale: a geometry-of-the-rule construction, not a physics-equation
   claim with an L0..L3 ladder slot — the partition is a property of the
   quadrature, produced solve-free.  The verifiable content is
   ``tests/sn/sweep/test_angular_cell_partition.py`` — the **direct**
   value gate on the producer, both arms, added 2026-08-11: a
   hand-written cumulative-weight reference (sphere) and the analytic
   equispaced-arc closed form :math:`e_k = \sin\theta\cos(\pi - k\Delta
   \omega)` (cylinder), each with a negative control (the uniform
   partition; the retired chord partition) per ``vv-principles`` #19,
   plus the closing identities, the march-orientation sign law, and a
   labelled control recording that :math:`M = 2` — i.e. every
   ``folded_product(·, 4)`` fixture — is structurally BLIND to the
   partition choice.  Then
   ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py``
   (:math:`\tau` = P2 applied to the partition, same two references),
   ``tests/sn/sweep/test_tau_arc_wellposedness.py`` (the P3 theorem and
   the attainable closed endpoint) and
   ``tests/sn/verification/mms/test_mms_ordering_blindness.py``
   ``::test_the_full_circle_double_cover_is_REFUSED_by_the_cell_partition``
   (the non-monotone-arc refusal).  All are ``foundation`` gates —
   software/structural invariants of a discrete construction.

   ⚠ Until 2026-08-11 the partition producer had **no value gate at
   all**: every listed test read :math:`\tau`, which is P2 *applied* to
   the partition, so a wrong partition that kept :math:`\tau` inside
   :math:`[0,1]` was visible only to the two cylinder :math:`\tau` rows.
   The recurrence's worst partial amplification
   :math:`A(M) = \max_m \prod_{k\le m}(1-\tau_k)/\tau_k` — the number
   quoted as "recurrence error-amplification" in the Q5.6.4
   adjudication — is likewise now committed, in
   ``tests/sn/sweep/curvilinear/test_psi_half_positivity.py``.

Both branches are **derived, not conventional**, and they are derived
from *different* facts:

* **Sphere** — a Gauss--Legendre weight *is* the cell's
  :math:`\mu`-measure, so accumulating weights from :math:`\mu_{1/2} = -1`
  partitions :math:`[-1, 1]` exactly.  This is
  :cite:`BaileyMorelChang2010` Eq. (12) **verbatim**, corroborated
  independently by Lathrop 2000 p. 249 (:math:`\sum \Delta\mu_m = 2`).
  Unchanged by Q5.6.4.
* **Cylinder** — the azimuthal march is a march in :math:`\omega`,
  **arc by arc**, so the cell boundary is the midpoint *in the variable
  the march marches in*.  Taking the midpoint of the **stored**
  :math:`\omega` values keeps it correct for any monotone arc; on an
  equispaced-:math:`\omega` rule it is exactly the half-angle boundary
  :math:`\omega_m \pm \Delta\omega/2`.

Feeding the equispaced-:math:`\omega` case of
:eq:`angular-cell-partition` through :eq:`morel-montry-closure` gives a
closed form, which is the cylinder arm's **structurally independent**
reference in the producer-equivalence gate (it shares no code path with
the producer):

.. math::

   \tau_m \;=\; \tfrac12
     \;+\; \tfrac12\,\cot\omega_m\,\tan\!\bigl(\Delta\omega/4\bigr) .

`[M]` (2026-08-11) agreement with the producer on
``folded_product(n_mu=4, n_phi)``, maximised over all four levels:
:math:`1.1\mathrm{e}{-16}` / :math:`2.2\mathrm{e}{-16}` /
:math:`7.8\mathrm{e}{-16}` / :math:`7.4\mathrm{e}{-15}` /
:math:`2.3\mathrm{e}{-14}` at :math:`n_\varphi = 4/8/16/32/64`.  ⚠ It
**degrades with refinement** (:math:`\arctan2`/:math:`\cos` round-off in
both paths), so the gate asserts ``atol=1e-13`` rather than a
machine-epsilon bound; a row beyond :math:`n_\varphi = 64` must widen it
(``vv-principles`` #16 — never assert tighter than the producer
achieves).

.. warning:: **Do NOT "unify"** :math:`\alpha` **and** :math:`\tau`
   **onto one expression.**

   Both reference the same partition — that is the point of
   :eq:`angular-cell-partition`, and deriving it twice is exactly the
   defect Q5.6.4 fixed.  But they impose **two different conditions** on
   it: :math:`\tau` the **zeroth** moment (P2, above), :math:`\alpha` the
   **first**-moment conservation recursion (P4,
   :eq:`alpha-dome-recursion`; Hébert 3.397--3.399, after Alcouffe &
   O'Dell).  Forcing :math:`\alpha` to equal the geometric tangential
   cosine at these edges silently drives Lathrop's defect
   :math:`\delta \to 0`, i.e. :math:`\tau \to \tfrac12` — the angular
   *diamond* scheme (Hébert 3.406/3.431), a **different method** with a
   different diffusion limit.  ⚠ Hébert's own
   :math:`\eta_{p,q\pm1/2}` is a constant-flux conservation recursion,
   **not** a trig evaluation at a bisected :math:`\omega`; the closed
   form above is a theorem about *our* equispaced-:math:`\omega` rule,
   not the literature's definition of the partition.

.. _tau-ownership-note:

.. note:: **Where τ lives.**

   :math:`\tau_m` is **not** a
   :class:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator`
   field.  The geometry-layer primitive carries the SPATIAL curvature
   coefficients (``face_areas``, ``delta_A``) and a reference to the
   ANGULAR factor
   (:class:`~orpheus.sn.angular.redistribution.AngularRedistribution`,
   which owns the :math:`\alpha`-dome and :math:`\mu_{\rm start}` as of
   the 2026-08-26 un-weld); the M-M closure weight is produced by
   :func:`~orpheus.sn.angular.closure.morel_montry_tau_per_level`
   reading the single partition producer
   :func:`~orpheus.sn.angular.closure.angular_cell_edges_per_level`,
   which :class:`SNMesh` calls against the quadrature and its own
   ``self.coord``.  The split is deliberate: τ is a property of
   the *angular closure scheme*, selectable per mesh, while the curvature
   coefficients are a property of the *geometry*.  Statements elsewhere
   in the corpus describing ``tau_mm`` / ``tau_mm_per_level`` as factory
   outputs predate that move, and any naming
   ``morel_montry_tau_raw_per_level`` predate Q5.6.4.

The admissible range is :math:`\tau \in [0, 1]` — predicate **P3**, an
ordinate lies inside its own angular cell — **enforced since Q5.5
(2026-08-07)**: the producer RAISES on :math:`\tau \notin [0, 1]`.  On a
well-posed monotone march membership certifies the march; a value outside
certifies an ILL-POSED march.  Both realized cases were caught by
measurement: (a) mis-ordered members — T22's ω-ordered mis-ordering
measured :math:`\tau = 1.079` at the producer before surfacing as a NaN
400 lines downstream, the pre-Q5.5 absorption silently laundering it into
a finite wrong answer; (b) a quadrature incompatible with the arm's edge
convention — a raw 3-D ``level_symmetric(4)`` rule fed to the 1-D
spherical arm (24 unsorted ``mu_x`` with duplicates, weights summing to
:math:`4\pi`) measured :math:`\tau \in [-20.3,\, 1.13]` with 23 of 24
ordinates outside, consumed SILENTLY by the *unclamped* sphere closure
until the guard landed — seven operator-equivalence tests ran exactly
this configuration and stayed green because both compared spellings share
the :math:`\tau` (the Mode-12 annihilation).  Issue #336 tracks the
refuse-or-reduce design for ``SNMesh`` on a spherical mesh with a
non-μ-line rule.  The closed endpoints are legal march starts — ``0`` is
an edge-node start and ``1`` an η-degenerate tie
(:func:`~orpheus.sn.angular.closure.march_start_structure_per_level`).
The guard does NOT catch the double cover: a full-circle level's
:math:`[0, 1, 0, 1, \ldots]` fingerprint is entirely inside
:math:`[0, 1]`; that detector is the singular set :math:`\Sigma`, and
since Q5.6.4 the non-monotone arc is refused one frame earlier still, by
:eq:`angular-cell-partition`'s own producer (a full-circle level carries
:math:`\omega` of both signs, so "the midpoint in :math:`\omega`" is not
defined for it).

⭐ **On the cylinder, P3 is now a THEOREM.**  With the partition taken as
the ω-midpoint, a strictly ω-monotone level has
:math:`\omega_{m-1/2} > \omega_m > \omega_{m+1/2}` by construction, and
:math:`\cos` is monotone on :math:`(0, \pi)`, so
:math:`\eta_{m-1/2} < \eta_m < \eta_{m+1/2}` — :math:`\tau \in (0, 1)` is
**forced** (`[M]` 4000 random monotone arcs: :math:`\min\tau =
4.739\mathrm{e}{-7}`, :math:`\min(1-\tau) = 7.599\mathrm{e}{-10}`, never
outside).  Its only equality case is a node ON an arc endpoint, i.e. a
node on :math:`\Sigma` — so **cylinder-P3 reduces to the fold criterion**
:math:`\Sigma = \emptyset`.  P3 keeps its teeth on the **sphere**, where
cumulative-weight edges genuinely need not bracket their nodes (case (b)
above).  Stated plainly so no audit reads cylinder-P3 as live coverage it
is not.

.. _sn-tau-absorber-retirement:

The retired :math:`[\tfrac12, 1]` absorber — what it was compensating for
-------------------------------------------------------------------------

**The clamp was TWO objects welded together** (T27, adjudicated
2026-08-02; membership guard landed Q5.5, 2026-08-07; absorber retired
Q5.6.4, 2026-08-11).  The cylinder-only expression
``max(0.5, min(1.0, τ_raw))`` fused the :math:`[0, 1]` *membership* — now
the P3 guard above — to a :math:`[\tfrac12, 1]` *absorption* whose stated
purpose was blocking an edge-node division.  The
:math:`[\tfrac12, 1]` box was **never** the admissible range of
:math:`\tau`: the sphere ran outside it, unclamped and correct, and `[M]`
at :math:`S_8` Gauss--Legendre **four of eight** M-M τ sit below
:math:`\tfrac12`.  :cite:`BaileyMorelChang2010`'s own :math:`S_2` example
gives :math:`\tau_1 = \mu_1 + 1 = 1 - 1/\sqrt3 \approx 0.4226 < \tfrac12`
(their Eq. 47).  **No source prescribes any limiter on** :math:`\tau`.

.. _sn-tau-absorber-provenance:

Where the number :math:`[\tfrac12, 1]` actually came from
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"No source prescribes it" is the weaker statement.  The stronger one,
established 2026-08-11: the interval is real, it is **Grant's**, and it
is on a **different parameter**.

:cite:`ReedLathrop1970` carry two independent weighted-diamond
parameters — a **spatial** weight :math:`a_{i+1/2}` (their Eqs. 5a/5b
and 11) and an **angular** weight :math:`\tau_m` (their Eqs. 6a/6b).
Their footnote 8 (printed p. 239) discusses Grant's choice for the
*spatial* one: Grant lets that weight depend on the sign of its
argument, and Reed & Lathrop note this "is necessary only to keep
[it] between :math:`\tfrac12` and 1" — then add, decisively, that
**Grant does not determine angular weights at all**.

So :math:`[\tfrac12, 1]` is a bound on the SPATIAL weighted-diamond
parameter of Grant, I. P. (1968), *J. Comp. Phys.* 2(4):381-402,
doi:10.1016/0021-9991(68)90044-2.  ⚠ **Grant 1968 is not in the local
library and has not been read**; the attribution rests on Reed &
Lathrop's footnote, read on the rendered page.  Transplanting that
interval onto the angular :math:`\tau` is exactly what the retired
cylinder absorber did — the number was inherited across a parameter
boundary, which is why no search of the angular literature could ever
find its source.

⭐ **What it was actually compensating for was a WRONG PARTITION, and
that is why retiring it alone made things worse.**  Until 2026-08-11 the
cylinder edges were taken at the midpoint of consecutive :math:`\eta`
values — the **chord** midpoint — with the endpoints pinned at
:math:`\mp\sin\theta`.  That partition is :eq:`angular-cell-partition`
*with its end cells stretched*.  The identity

.. math::

   \tfrac12\bigl[\cos\omega_a + \cos\omega_b\bigr]
     \;=\; \cos\!\Bigl(\tfrac{\omega_a+\omega_b}{2}\Bigr)
            \cos\!\Bigl(\tfrac{\Delta\omega}{2}\Bigr)

(the :math:`\kappa` prefactor's sibling) makes every *interior* chord
edge exactly :math:`\cos(\Delta\omega/2) \times` the arc edge (`[M]`
agreement :math:`10^{-16}`), while the two endpoints stay **unscaled** —
so the outermost cells stretch to absorb the shrink.  The :math:`\eta`
error vanishes as :math:`\Delta\omega \to 0`, but the implied
:math:`\omega`-width spread does **not**: it converges to
:math:`\approx 17.45\,\%` (`[M]` 18.71 / 17.59 / 17.48 / 17.46 % at
:math:`n_\varphi = 8/16/32/64`) against a quadrature whose own cells are
bit-exactly equal.  That :math:`O(1)` inconsistency — one object, the
"boundary between azimuthal cell :math:`m` and :math:`m+1`", derived
independently by :math:`\alpha` (at the real half-angle
:math:`\omega_{m-1/2}`) and by :math:`\tau` (at the chord midpoint), in
disagreement — is what the absorber hid.

**The absorber is condemned on its own terms, with no MMS involved.**  The
ν-closure diagnostic marches the level *implied by* :math:`\tau`
(:math:`\nu_{1/2} = -\sin\theta`,
:math:`\nu_{m+1/2} = (\eta_m - (1-\tau_m)\nu_{m-1/2})/\tau_m`) and asks
whether it lands on :math:`+\sin\theta`.  It is solve-free and it
separates a derived τ from a fabricated one — `[M]`
:math:`\nu/\sin\theta` at close:

.. list-table:: ν-closure: does the march implied BY τ close the level?
   :header-rows: 1
   :widths: 16 21 21 21 21

   * - :math:`n_\varphi`
     - arc ω (production)
     - chord (retired)
     - clamped (retired)
     - :math:`\tau \equiv \tfrac12`
   * - 8
     - ``1.000000``
     - ``1.000000``
     - **``1.016389``**
     - **``1.164784``**
   * - 16
     - ``1.000000``
     - ``1.000000``
     - ``1.001930``
     - ``1.039182``
   * - 32
     - ``1.000000``
     - ``1.000000``
     - ``1.000238``
     - ``1.009677``
   * - 64
     - ``1.000000``
     - ``1.000000``
     - ``1.000030``
     - ``1.002412``

⟹ the clamped values and :math:`\tau \equiv \tfrac12` correspond to **no
partition of the level at all** — their implied march overshoots the
level's own endpoint by 1.6 % and 16.5 % respectively.  (:math:`\tau
\equiv \tfrac12` is Hébert's angular *diamond* scheme; it is listed here
because it is the tempting rescue, and because it optimises truncation
order while breaking the diffusion limit :math:`\tau` exists to fix.)

**Why the sphere's convention cannot simply be transplanted.**
Accumulating weights in :math:`\eta` — even correctly renormalised to
:math:`\sum \bar w = 2\sin\theta`, :cite:`BaileyMorelChang2010` Eq. 52 —
violates **P3** on our rule and **worsens with refinement**: `[M]`
ordinates outside their own cell go 0/4 → 4/8 → 12/16 → 28/32 at
:math:`n_\varphi = 8/16/32/64` **per level** (`[M]` 0/16, 16/32, 48/64,
112/128 over the four levels of ``folded_product(n_mu=4, n_phi)``, with
:math:`\tau` ranging out to :math:`[-2.86,\, 3.86]` at
:math:`n_\varphi = 64`), and the solve diverges (NaN) from
:math:`n_\varphi \ge 16`.  The reason is structural: an arc cell's
:math:`\eta`-measure is
:math:`2\sin\theta\,\sin\omega_m\,\sin(\Delta\omega/2) \propto
\sin\omega_m`, **not** constant — `[M]` at :math:`n_\varphi = 16` it
spans :math:`0.30`--:math:`1.53 \times` the uniform width across one
level (and :math:`0.08`--:math:`1.57\times` at
:math:`n_\varphi = 64`) — while a trapezoid weight is.  ⟹ Eq. (52) is
not a law; it is the
statement that *in their* quadrature the weight equals the cell's
:math:`\eta`-measure.  Ours does not, so we satisfy the same **predicate**
by a different partition.  (This diagnosis is not new — it was written
into the original M-M implementation as *"weights are uniform in*
:math:`\varphi` *-space, not* :math:`\eta` *-space"*, and is preserved
verbatim in
:mod:`orpheus.derivations.discrete.sn.angular_differencing`.  The
diagnosis was right; the chord-midpoint *substitute* was never checked
against a partition predicate.)

**On the σ_y-folded arc the absorption's reason is structurally gone.**
With midpoint nodes the smallest-η ordinate sits at
:math:`\omega_0 = \pi - \Delta\omega/2`, so the closed form above gives

.. math::
   :label: morel-montry-folded-arc

   \tau_0
     \;=\; \tfrac12 - \tfrac12\cot\!\bigl(\Delta\omega/2\bigr)
            \tan\!\bigl(\Delta\omega/4\bigr)
   \;\xrightarrow[\;n_\varphi \to \infty\;]{}\; \tfrac14
   \;\;\text{from inside},
   \qquad
   \tau_m \in \Bigl[\tfrac14,\, \tfrac34\Bigr],
   \qquad
   \tau_m + \tau_{M-1-m} \;=\; 1 ,

since :math:`\cot(\Delta\omega/2)\tan(\Delta\omega/4) \to \tfrac12`.
`[M]` (2026-08-11, :math:`n_\mu = 4`, folded staggered), with
:math:`|\Sigma| = 0` throughout — :math:`\tau \in \{0, 1\}` is
**structurally unreachable** on the fold:

.. list-table:: Folded-arc τ box and the reversal identity
   :header-rows: 1
   :widths: 20 45 35

   * - :math:`n_\varphi`
     - τ range
     - reversal residual
   * - 4
     - ``[0.292893, 0.707107]``
     - 0.5 ULP
   * - 8
     - ``[0.259892, 0.740108]``
     - 0.5 ULP
   * - 16
     - ``[0.252425, 0.747575]``
     - 2.0 ULP
   * - 32
     - ``[0.250603, 0.749397]``
     - 7.0 ULP
   * - 64
     - ``[0.250151, 0.749849]``
     - 12.0 ULP

.. note:: **Retraction (2026-08-11, Q5.6.4).**  Until Q5.6.4 this
   equation read :math:`\tau_{{\rm raw},0} \to \tfrac15`,
   :math:`\tau_{{\rm raw},m} \in [\tfrac15, \tfrac45]`, with the
   reversal identity holding **bit-exactly** (residual :math:`0.0` at
   every :math:`n_\varphi`; measured :math:`[0.2195, 0.7805]` at
   :math:`n_\varphi = 8` falling to :math:`[0.200289, 0.799711]` at 64).
   Those numbers are correct — *for the retired chord partition they were
   measured on*.  The box is now :math:`[\tfrac14, \tfrac34]`.

   ⚠ **And the reversal identity is no longer bit-exact: that is a trade,
   not a regression.**  The chord partition's reversal symmetry was exact
   *because* both end cells were stretched symmetrically — the 17.5 %
   ω-width defect cancelled itself under
   :math:`\omega \to \pi - \omega`.  The ω partition has the correct cells
   and pays 0.5--12 ULP of :math:`\arctan2` / :math:`\cos` round-off.  The
   gate asserts 64 ULP, because the residual grows with arc refinement and
   a bit-exact assertion would be a latent false red
   (``vv-principles`` #16 — never assert tighter than the producer
   achieves).

.. (vv-status rationale) morel-montry-folded-arc: Verified by the Q5.5
   mechanism gates in ``tests/sn/sweep/test_tau_arc_wellposedness.py``,
   re-posed at Q5.6.4 —
   ``test_the_fold_mechanism_is_an_empty_singular_set`` asserts the
   MECHANISM (Σ = ∅, computed via ``singular_set``, never declared) and
   ``test_the_folded_tau_is_bounded_with_the_reversal_identity`` the
   CONSEQUENCE (τ ⊂ [1/4, 3/4] per level plus the reversal identity at
   64 ULP), each at n_φ ∈ {8, 16, 32, 64} on the folded staggered product
   (n_μ = 4).  Teeth measured 2026-08-07: reverting Q5.2's offset to
   δ = 0 reds BOTH legs at every n_φ (8/10 red), which attributes the
   pass to the mechanism; the [0, 1] guard's negative companion reds
   alone under a no-opped guard (1/10 red).  A geometry-of-the-rule
   invariant, not a physics-equation claim with an L0..L3 ladder slot.
.. vv-status: morel-montry-folded-arc documented

**The honest accuracy cost, ratified rather than hidden.**  The
principled partition is not uniformly more accurate on the one
manufactured fixture available: `[M]` on the anisotropic cylindrical MMS
at :math:`n_x = 320` it is BETTER at :math:`n_\varphi = 8`
(:math:`3.128\mathrm{e}{-3}` vs :math:`3.511\mathrm{e}{-3}`) and
:math:`\sim 1.8`--:math:`2\times` WORSE at :math:`n_\varphi = 16/32/64`.
Principled :math:`\ne` more accurate: a scheme satisfying P2/P3 wins over
one with a smaller number on a single MMS, and the L2 norm measures
truncation order — exactly what :math:`\tau \equiv \tfrac12` optimises
and exactly the quantity that is blind to the diffusion limit
:math:`\tau` exists to fix.  The full ladder, and why the MMS is the
wrong instrument for this decision, is at
:ref:`sn-cylinder-angular-floor` in :doc:`/theory/verification/sn`.

The predicate ladder (P0--P4) these decisions are made against, the
:math:`\tau`/:math:`\beta` nomenclature (both letters are overloaded, and
both collisions have cost real time), and a written record of **which
diagnostics are blind on which rules** live in
:mod:`orpheus.derivations.discrete.sn.angular_differencing`.

API surface
-----------

The geometry-layer primitive is built by three factory functions, one
per coordinate system:

* :func:`~orpheus.sn.mesh.reduced_operator.slab_streaming(mesh, ang)
  <orpheus.sn.mesh.reduced_operator.slab_streaming>` — Cartesian 1-D;
  no curvature math.  Both its factors are present and **neutral**: the
  ANGULAR one carries a zero dome and the diameter-ray start, and the
  SPATIAL one carries a unit cross-section with **zero area change**
  (``face_areas == ones(nx+1)``, ``delta_A == zeros(nx)``) — because a
  slab having "no curvature" IS its faces not changing area.

  ⛔ Until P4.1b (2026-08-27) this read *"its spatial arrays remain
  ``None``"*, and they were stored fields.  They are **derived** from
  the mesh now (``face_areas`` is ``mesh.areas`` itself, ``delta_A`` its
  difference), so no factory computes them and the per-coordinate
  ``Optional`` union is dead on **both** factors rather than one.  That
  is what let :meth:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.streaming_terms`
  collapse from three chart arms to one shared body: the retired
  CARTESIAN arm was the spherical body with ``1.0`` / ``1.0`` / ``0.0``
  written out by hand.
* :func:`~orpheus.sn.mesh.reduced_operator.spherical_streaming(mesh, ang)
  <orpheus.sn.mesh.reduced_operator.spherical_streaming>` — 1-D spherical
  with the dome recursion :eq:`alpha-dome-recursion` and Morel--Montry
  closure :eq:`morel-montry-closure`.
* :func:`~orpheus.sn.mesh.reduced_operator.cylindrical_streaming(mesh, ang)
  <orpheus.sn.mesh.reduced_operator.cylindrical_streaming>` — 1-D
  cylindrical with **per-**\ :math:`\mu`\ **-level** :math:`\alpha`,
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
:meth:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.streaming_terms`,
which returns a
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` dataclass
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
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` are
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
   :func:`~orpheus.sn.mesh.reduced_operator.spherical_streaming` /
   :func:`~orpheus.sn.mesh.reduced_operator.cylindrical_streaming`
   itself, so the surviving hash-equality legs compare a fresh factory
   call against ``sn_mesh.reduced`` — *the value that same factory
   produced*, routed through the mesh constructor.

   ⛔ This paragraph used to add *"and the two ``SNMesh.face_areas`` /
   ``SNMesh.delta_A`` legs are deprecated read-throughs to that same
   object"*.  Those accessors **retired at P4.1c** (2026-08-27) — `[M]`
   11 readers, **0** of them in ``orpheus/``, and every one of the
   tests that read them existed to verify the shims themselves.  The
   legs now read ``sn_mesh.reduced.*`` directly.  The gate therefore pins the
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
     :math:`4\pi\,\Delta(r^2)` / :math:`2\pi\,\Delta r`.
     ⛔ This entry read "**sole catcher**; the snapshots are blind here,
     and correctly so — ``delta_A`` has no production consumer."  True
     until 2026-08-26, **false now**: retiring the fused ``redist_dAw``
     cache made ``delta_A`` the spatial factor that BOTH
     ``streaming_terms`` and the angular closure read, so every
     curvilinear snapshot rides on it.
   * ``angular.alpha_per_level`` (was ``alpha_half`` /
     ``alpha_per_level``, until the 2026-08-26 un-weld moved the dome to
     the angular factor) — the L0 per-ordinate flat-flux identity
     ``test_per_ordinate_flat_flux_consistency`` on both arms
     (``catches("ERR-006", "ERR-007")``);
     ``tests/sn/sweep/curvilinear/test_alpha_closed_form.py`` (the
     Dirichlet-kernel closed form; **cylindrical α only** — every
     fixture there is ``CoordSystem.CYLINDRICAL``); plus both snapshot
     families.
   * ``redist_dAw`` / ``redist_dAw_per_level`` — **RETIRED 2026-08-26**
     as a fused product neither of its two consumers owned.  Its catcher
     — ``tests/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py``,
     the L0 closed-form :math:`\varphi = Q/(\Sigma_t(1-c))` identity —
     still covers the QUANTITY, now formed at each consumer from
     ``delta_A`` and the weights.  ⭐ And the historical note that the
     flat-flux identity "recomputes :math:`\Delta A / w` rather than
     reading the production array" is exactly why that gate needed no
     migration: it was already forming the product, not reading the
     cache.
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
``tau_mm``, ``tau_mm_per_level``) are gone — and since 2026-08-26 the
first four have no ``self.reduced`` field left to route to either, the
α-dome having moved to the angular factor and ``redist_dAw`` having
retired as a fused product.  Consumers bind to
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

.. _who-needs-a-connection-coefficient:

Who needs a connection coefficient — and who does not
-----------------------------------------------------

The primitive above is **structurally S**\ :sub:`N`\ **-only**, and that
is a statement about the mathematics, not about how far a migration has
got.  The distinction matters because the two readings license opposite
work: *"MoC and CP have not migrated yet"* invites someone to go and
wire them up, while *"MoC and CP never form this term"* says the
primitive is correctly placed and correctly consumed by exactly one
solver family — S\ :sub:`N`.  (`[M]` the whole
:class:`~orpheus.transport.spatial.scheme.StreamingTerms` /
:class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` chain is
referenced only from ``orpheus/sn/`` and from inside
``orpheus/transport/`` itself, and by **no** file under
``orpheus/moc/``, ``orpheus/cp/``, ``orpheus/mc/`` or
``orpheus/diffusion/``.  As above, the count is deliberately not frozen
here — re-run the predicate at
:ref:`connection-coefficient-census` with
``r"transport\.spatial|DiscretizationScheme|StreamingTerms|CellVisit"``
substituted for its pattern set; the two controls stay non-zero and the
four subjects stay at zero.)

A solver family needs an :math:`\alpha` dome **iff all three** of the
following hold.

1. **It carries an angular unknown** — a :math:`\psi` that survives
   discretisation still wearing a direction index.
2. **That index is read in a local, rotating frame** — a basis that
   turns as the spatial point moves, so a particle streaming in a fixed
   *physical* direction changes its *coordinate* label as it travels.
3. **The resulting angular derivative is discretised by collocation**
   on that index — :math:`\partial_\mu` (sphere) or
   :math:`\partial_\varphi` (cylinder) approximated as a difference
   between neighbouring ordinates, rather than by an expansion in a
   basis that differentiates exactly.

Condition 2 is what mints the term at all.  In a curvilinear chart the
direction cosines are measured against the *local* radial and azimuthal
axes, so straight-line streaming is a continuous relabelling of
:math:`(\mu, \varphi)`; the redistribution term
:math:`(1-\mu^2)/r\,\partial_\mu` (sphere) or
:math:`-(1/r)\,\partial_\varphi(\xi\,\cdot)` (cylinder) is the
bookkeeping for that relabelling.  Condition 3 is what turns its weight
into a *recursion*: collocation supplies no exact derivative, so the
weights must instead be built to preserve the one property that has to
survive — a spatially flat angular flux must redistribute to zero, per
ordinate — and :eq:`alpha-dome-recursion` together with its closure
contract (:ref:`sn-alpha-dome-closes`) is precisely the construction
that delivers it.

Adjudicating the shipped families against those three conditions:

.. list-table:: Which families satisfy the three conditions
   :header-rows: 1
   :widths: 15 19 22 20 24

   * - Family
     - (1) angular unknown
     - (2) local rotating frame
     - (3) collocated angular derivative
     - Consequence
   * - S\ :sub:`N`, curvilinear
     - yes — :math:`\psi_{n,i}`
     - yes — :math:`(\eta, \xi, \mu)` on the local radial frame
     - yes — the half-angle recursion
     - **needs the dome**
   * - MoC
     - yes — :math:`\psi` per track
     - **no** — :math:`\Omega` is fixed in the GLOBAL frame
     - n/a
     - term relocates into track segmentation
   * - CP
     - **no** — angle is integrated out before discretisation
     - n/a
     - n/a
     - term never appears
   * - MC
     - **no** — directions are sampled, not indexed
     - n/a
     - n/a
     - term never appears

**MoC fails condition 2.**  The method of characteristics is *defined*
by choosing the global frame in which :math:`\Omega` is constant along a
track, so :math:`\Omega \cdot \nabla \psi = \mathrm{d}\psi/\mathrm{d}s`
is chart-free and there is no angular derivative left to discretise.
Curvature does not disappear — it moves into *segmentation*, the
ray-region intersection that produces the chord lengths.  `[M]` the
shipped inner loop in :mod:`orpheus.moc.core` forms
:math:`\tau = \Sigma_t \, \ell_{\rm seg} / \sin\theta_p` per segment and
applies plain exponential attenuation
:math:`\Delta\psi = (\psi - Q/\Sigma_t)\,(1 - e^{-\tau})`; no ordinate
couples to its neighbour anywhere in the sweep.  This is also why
:mod:`orpheus.sn.loss_representation.sweep_graph` records that MoC will
define a *per-ray traversal* analog rather than reuse the
S\ :sub:`N` sweep graph.

**CP fails condition 1.**  Collision probability integrates the angular
variable analytically *before* anything is discretised — the transport
kernel is already an angle-integrated function of optical path, so no
angular unknown, and therefore no angular index, ever exists.  `[M]`
:class:`orpheus.cp.solver.CPSolver` dispatches on
:class:`~orpheus.geometry.coord.CoordSystem` to three real setups —
slab, cylinder, and **sphere**, i.e. the curvilinear cases the false
claim was about — and each installs a scalar kernel:
:math:`F(\tau) = e^{-\tau}` with a :math:`y`-weighted quadrature on the
sphere, the :math:`\mathrm{Ki}_3` kernel on the cylinder, :math:`E_3`
on the slab.  There is no :math:`\alpha`, no :math:`\Delta A / w`, and
nothing for either to act on.

**MC fails condition 1 as well, for a different reason.**  Monte Carlo
samples directions from a continuous distribution rather than indexing
a fixed set, so there is no neighbouring ordinate to redistribute *to*.
`[M]` :class:`orpheus.mc.solver.MCMesh` admits ``CARTESIAN`` or
``CYLINDRICAL`` — so it, too, solves a curvilinear problem with no
:math:`\alpha` anywhere.

⭐ **The curvilinear counter-examples are the load-bearing evidence.**
"MoC and CP have not migrated yet" predicts that neither has a
curvilinear capability to migrate.  Both do.  `[M]`
:class:`orpheus.moc.geometry.MOCMesh` wraps a **cylindrical**
:class:`~orpheus.geometry.mesh.Mesh1D` and ray-traces **concentric
annuli** (``_ray_circle_intersections``); :class:`orpheus.cp.solver.CPSolver`
ships a real **sphere**; and, as a third witness,
:class:`orpheus.mc.solver.MCMesh` ships a real **cylinder**.  All three
solve curved geometry, and none carries one line of redistribution
machinery.  A capability that exists *and* declines the primitive
refutes the migration reading in a way that an absent capability never
could — which is exactly why the claim survived unchallenged for so
long: nobody had looked at what those packages already do.

.. _connection-coefficient-census:

Reproducing the census — the predicate, not a table of counts
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The census behind the `[M]` claims above is published as the **recipe,
not as a table of counts**.  Eight independent spellings of the concept
— three symbol spellings, three concept spellings, an attribute-access
spelling and the prose paraphrase — counted as regex *occurrences* over
``*.py`` under each package root:

.. code-block:: python

   import pathlib, re

   PATTERNS = {
       "reduced_operator":         r"reduced_operator",
       "ReducedStreamingOperator": r"\bReducedStreamingOperator\b",
       "AngularRedistribution":    r"\bAngularRedistribution\b",
       "alpha":                    r"\balpha\b",
       "delta_A / face_areas":     r"\b(?:delta_A|face_areas)\b",
       "redistribut*":             r"\bredistribut\w*",
       ".reduced":                 r"\.reduced\b",
       "connection coefficient":   r"connection[ -]coefficient",
   }
   SUBJECTS = ["orpheus/moc", "orpheus/cp", "orpheus/mc"]
   CONTROLS = ["orpheus/sn", "orpheus/geometry"]

   def count(pkg, pat):
       return sum(len(re.findall(pat, f.read_text(encoding="utf-8")))
                  for f in pathlib.Path(pkg).rglob("*.py"))

   for name, pat in PATTERNS.items():
       # POSITIVE CONTROL: a zero here means the filter is broken,
       # not that the tree is clean.
       assert all(count(p, pat) > 0 for p in CONTROLS), name
       # THE FINDING:
       assert all(count(p, pat) == 0 for p in SUBJECTS), name

**Run as written on 2026-08-27 it passes: every one of the eight
patterns is non-zero in both controls and exactly zero in all three
subjects.**  The zeros are the finding and they are falsifiable — the
day one of them stops being zero, the claim on this page is refuted and
should be re-argued rather than patched.

.. note:: **Why the control COUNTS are deliberately not printed here.**
   A positive control has to be **non-zero**; its particular value
   carries no part of the argument.  Freezing it would put a number in
   the corpus that moves under any edit to ``orpheus/sn`` or
   ``orpheus/geometry`` — including the ⛔ correction blocks this very
   pass added, which name the module and so raise several of the control
   counts.  ⛔ An earlier revision of this section did print them, and
   did so wrongly in three independent ways at once: the counts were
   taken **before** this pass's own edits; the column set was a
   **different** partition of "six spellings" from the one the prose
   beside it named (so two of the six numbers belonged to spellings the
   prose never mentioned); and the ``redistribut`` column was
   case-insensitive and unanchored, which silently absorbed every
   ``AngularRedistribution``.  All three are the failure mode this page
   exists to document, so the section now points at a re-runnable
   predicate instead (`plan-authoring` §9: never copy a number the tree
   can re-measure; `plan-authoring` §2: a number without its predicate
   is not re-runnable).

The set of files that name the primitive at all is small enough to
enumerate, and an enumeration — unlike a count — can be checked by
reading it.  `[M]` 2026-08-28, re-run after P4.2 + P4.3 completed the
un-weld (an earlier revision of this list, re-run after P4.4 only, counted
**15** and named ``geometry/reduced_operator.py`` a definer): **14** files
under ``orpheus/`` name any of ``ReducedStreamingOperator``, the three
``*_streaming`` factories, ``StreamingTerms``, ``AngularRedistribution``,
``angular_redistribution``, ``alpha_dome`` or ``AngularMeasure``.  Three
are the definers (``sn/mesh/reduced_operator.py`` — the connection
operator and factories; ``transport/spatial/scheme.py`` —
``StreamingTerms``, beside the contract that consumes it;
``sn/angular/redistribution.py`` — the :math:`\alpha` cluster and
``AngularMeasure``).  The rest: in ``orpheus/sn/``,
``angular/__init__.py``, ``angular/closure.py``, ``mesh/augmented_mesh.py``,
``operators/radial_characteristic.py``, ``solver.py`` and
``sweep/cache.py``; in ``orpheus/transport/spatial/``, ``__init__.py``,
``cell_balance.py``, ``diamond.py`` and ``linear_discontinuous.py``; and
``orpheus/derivations/discrete/sn/angular_differencing.py``.  Not one of
them is under ``orpheus/moc/``, ``orpheus/cp/`` or ``orpheus/mc/`` — the
load-bearing half, and it is unchanged.  ⭐ And for the first time not
one is under ``orpheus/geometry/`` either: the un-weld's own done-when,
now a property of this enumeration.

⚠ ``transport/spatial/linear_discontinuous.py`` is **new to this list and
was missing from it before P4.4**, not added by it: `[M]` it named
``StreamingTerms`` at the previous commit too.  An enumeration is only
checkable by re-running its own predicate, which is how the gap was
found — see `plan-authoring` §2 on a universal owing its denominator.

.. note:: **What WOULD change this answer.**  The three conditions are
   the claim, so the honest way to falsify it is to break one.  A
   discrete-ordinates scheme that expanded the angular flux in a basis
   which differentiates exactly — spherical harmonics, or a
   discontinuous-Galerkin / finite-element discretisation *in angle* —
   would satisfy 1 and 2 and fail 3, and would then need a different
   object entirely (a mass/stiffness pair in :math:`\mu`), not this
   recursion.  No such scheme exists in this codebase.  Conversely,
   neither MoC nor CP is likely to acquire condition 2 or condition 1
   respectively without ceasing to be MoC or CP.

.. note:: **Two labels, one recursion.**  :eq:`alpha-dome-recursion`
   here and :eq:`alpha-recursion` on
   :doc:`/theory/methods/sn/curvilinear_one_group` state the *same*
   recurrence.  This page carries it in the **geometry-primitive**
   register — what the chart owes a consumer, seeded at
   :math:`\alpha_{1/2} = 0` — and the methods page in the
   **discretisation** register, where it enters the cell update with
   :eq:`alpha-cylindrical` as the per-level arm.  Only the methods-page
   label is a ``verifies`` target.  The duplication is a genuine
   single-source-of-truth smell; collapsing it is **recommended and
   deliberately not done here**, because it would move a generated
   V&V-matrix row and re-point ``verifies`` markers that a
   documentation pass may not edit.

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
  connection coefficients off ``self.mesh.coord`` inside the walk.
  (Depth B consumed them through the per-geometry
  ``transport_operator_matvec_*`` matvecs; that family and its unified
  successor were deleted in the typed-field (#197) and walk-unification
  (#280 campaigns) refactors — the primitive itself is unchanged.)
* **MoC and CP campaigns (post-Wave-1)** — ⛔ **retracted 2026-08-27.**
  This item read *"reuse the same primitive with their own consumption
  patterns (track-segment chord march for MoC; ray-traced chord-length
  integrals for CP)"*.  It was never a description of the tree, and it
  is not a migration still owed: neither method forms an angular
  redistribution term at all, for two independent structural reasons
  worked in :ref:`who-needs-a-connection-coefficient` above.  The
  roadmap item is closed as **not applicable**, not as pending.  What
  MoC and CP *do* share with S\ :sub:`N` is the L2 transport layer —
  fields, sources, cross-section data, the scattering kernel and the
  eigenvalue driver (:mod:`orpheus.transport`) — which is a different
  and genuinely satisfied sharing claim.
