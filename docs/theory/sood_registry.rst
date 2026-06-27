.. _theory-sood-registry:

==========================================================================
Sood-family benchmark case registry — method-agnostic truth set
==========================================================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this before modifying the Sood case registry.**

- **What this is**: a registry, not a method. It is the single source
  of truth for benchmark case configurations (cross sections,
  geometry, published reference values) from the Sood-family
  literature — currently LA-13511 (Sood/Forster/Parsons 1999) and
  Atalay 1997. The folder name (``sood_registry/``) preserves the
  author name because it identifies a *case collection*, not a
  method — author-named registries are explicitly allowed under the
  project's folder-naming rule.
- **Architectural intent**: decouple **case configurations**
  (cross sections, geometry, published reference values) from
  **solver methods**. Multiple methods consume the same case via
  different adapters — see :ref:`sood-registry-purpose`.
- **Coverage today**:

  * **LA-13511**: 5 first-slice cases (PUa-1-0-IN, PU-2-0-IN,
    Ua-1-0-SL/CY/SP) + 30 wide-slice active cases + 12 wide-slice
    stub cases (cylinder + 2G bare-critical pending solver
    dispatch). Total: **42 LA-13511 cases.**
  * **Atalay 1997**: 6 slab + 1 sphere = 7 reflected /
    linearly-anisotropic cases.
  * **Total: 49 cases shipped.**
- **Production-protocol-aligned**: every case carries
  ``materials: dict[int, Mixture]`` + ``geometry_spec: GeometrySpec``
  + ``truth: La13511Truth`` — the same objects production solvers
  (:func:`orpheus.cp.solver.solve_cp`, :func:`orpheus.sn.solver.solve_sn`)
  consume directly.
- **Consumers**:

  * Semi-analytical reference solvers (F_N method, PS-1982 wrapper,
    transfer-matrix :math:`k_\infty`, WM-72 cylinder, Atalay
    slab/sphere) — extract numpy arrays via
    :func:`mixture_to_fn_arrays`.
  * Production discrete solvers (CP, SN, MOC) — accept the
    ``materials`` dict + ``geometry_spec.build()`` directly.
- **Cross-references**: :ref:`theory-fn-method` consumes
  ``LA13511_CASES`` for slab/sphere F_N pinning;
  :ref:`theory-singular-eigenfunction` consumes the Atalay catalogue;
  :ref:`theory-trajectory-resolvent` cross-checks Variant α against
  ``Ua-1-0-CY`` / ``Ua-1-0-SP``. All consumers share the same
  ``mixture_to_fn_arrays`` extractor below the trusted-library line.


.. _sood-registry-purpose:

Purpose — why a method-agnostic registry
=========================================

Before this package, each method (F_N, Variant α, Galerkin spectral,
WM-72 cylinder, Atalay slab/sphere) carried its own per-case data
alongside its solvers. This coupling was the wrong factoring for two
load-bearing reasons:

1. **The same Sood case (XS, geometry, published reference value)
   must be consumable by every method that wants to benchmark itself
   against it.** Tests for F_N method needed Sood XS in numpy-array
   form; tests for the production CP solver needed the same XS in
   :class:`Mixture` form; tests for Variant α needed both.
   Duplication of the XS data across each method's test-helper module
   meant the same numerical XS values lived in 5+ places, with
   manual sync required on every change.
2. **Adding a new method to an existing case required touching the
   method's per-case data, not just the method itself.** This
   created a coupling-asymmetry that violated the architectural
   principle "case data should outlive any specific method".

The registry decouples the two: the case lives in
:mod:`sood_registry`, the method lives in its own module
(:mod:`fn_method`, :mod:`singular_eigenfunction`, etc.), and the
adapter functions
(:func:`mixture_to_fn_arrays`, :func:`build_materials`,
:func:`build_mesh`, :func:`build_cp_params`) bridge the two. Adding
a sixth method to verify against an existing case requires only the
method's own implementation and a new test file consuming
``LA13511_CASES["case_id"]``; no edit to ``sood_registry`` itself.

Architectural rule of thumb
----------------------------

If a future method is adding a parameter that *every solver* on the
case would care about (e.g. published flux table, published
:math:`k_\infty`), it goes in :class:`La13511Truth`. If it's a
parameter only the new method cares about (e.g. F_N collocation
order, Variant α quadrature points), it stays in the new method's
config — never in the registry.

The full design rationale lives in the closeout memo
``.claude/agent-memory/method-implementer/sood_registry_phase_a_migration.md``.

Schema
======

.. _sood-registry-schema:

The :class:`La13511Case` dataclass carries six load-bearing fields:

.. list-table::
   :header-rows: 1
   :widths: 22 22 56

   * - Field
     - Type
     - Purpose
   * - ``case_id``
     - ``str``
     - Sood naming convention identifier (``<Material>-<Groups>-<Scattering>-<Geometry>``).
   * - ``problem_number``
     - ``int``
     - Sequential problem number within LA-13511.
   * - ``description``
     - ``str``
     - One-line human-readable description.
   * - ``materials``
     - ``dict[int, Mixture]``
     - Macroscopic cross sections keyed by material ID.
       Single-region cases use ``{0: Mixture(...)}``; multi-region
       cases (none in the first slice yet) add more keys.
   * - ``geometry_spec``
     - :class:`GeometrySpec`
     - Geometry kind + critical dimension (mfp + cm) + per-end
       boundary conditions. The ``build(n_cells)`` method produces
       a :class:`Mesh1D` at the requested refinement.
   * - ``scattering_order``
     - ``int``
     - Legendre order of the scattering kernel (0 = isotropic,
       1 = P_1, 2 = P_2). Most LA-13511 cases are isotropic.
   * - ``truth``
     - :class:`La13511Truth`
     - Published reference values: :math:`k_{\rm eff}` /
       :math:`k_\infty`, flux ratios, critical dimensions, etc.
   * - ``sood_table``
     - ``int``
     - LA-13511 table number where this case is tabulated.
   * - ``primary_reference``
     - ``str``
     - The peer-reviewed paper Sood cites as the source.
   * - ``notes``
     - ``str``
     - Free-form remarks (typo flags, conversion subtleties, …).
   * - ``provenance``
     - :class:`Provenance` ``| None``
     - Phase B addition (2026-05-04). Structured citation metadata
       bundling ``paper_id`` / ``paper_table`` / ``primary_reference``
       / ``notes``. Mirrors the legacy flat ``sood_table`` /
       ``primary_reference`` / ``notes`` triple — Phase B populates
       both forms; Phase F drops the flat fields. ``None`` only on
       cases not yet migrated.

The :class:`GeometrySpec` carries:

* ``geometry`` — one of ``"infinite"``, ``"slab"``, ``"sphere"``,
  ``"cylinder"``, ``"ISLC"``.
* ``critical_dimension_mfp``, ``critical_dimension_cm`` — published
  critical dimension in both unit systems. ``None`` for infinite-
  medium cases.
* ``n_groups`` — energy group count.
* ``mat_id`` — material identifier in the constructed mesh.
* ``bc_left``, ``bc_right`` — per-end boundary conditions
  (defaults: vacuum at outer surface; reflective at sphere /
  cylinder centreline; vacuum at slab face).

The :class:`La13511Truth` carries:

* ``k_eff_or_kinf`` — reference :math:`k_{\rm eff}` (finite cases —
  usually 1.0, critical) or :math:`k_\infty` (infinite cases).
* ``flux_ratios`` — for 1G cases with a published flux table:
  ``Mapping[float, float]`` mapping :math:`r/r_c` to
  :math:`\phi(r)/\phi(0)`.
* ``flux_ratio_groupwise`` — for 2G infinite cases:
  ``Mapping[int, float]`` mapping ORPHEUS group index to
  :math:`\phi_g/\phi_0` (ratio relative to ORPHEUS group 0 = fast).
* ``angular_flux_at_surface`` — reserved for future cases that
  publish surface angular flux :math:`\psi(\mu, r=R)`. ``None`` for
  first-slice cases.
* ``critical_dimension_mfp`` — Phase B addition. Published critical
  dimension in mean free paths. For ``"slab"``: the half-thickness
  :math:`a` (F_N convention). For ``"sphere"`` / ``"cylinder"``: the
  radius :math:`R`. For ``"infinite"``: ``None``. This is a
  registry-truth value (the published critical configuration);
  living on ``Truth`` and not on ``geometry_spec`` mirrors the
  architectural fact that the value is a *truth claim* ("at this
  size, the configuration is critical"), not a geometric description.
  Multiply by :math:`1 / \Sigma_t` to convert to cm — that is what
  :meth:`La13511Case.to_geometry` does internally.
* ``extrapolated_endpoint_mfp`` — Phase B addition. Published
  extrapolated endpoint :math:`z_0` in mean free paths, used in
  diffusion-theory boundary conditions. Optional metadata; not all
  cases publish it. Currently ``None`` on every catalogue case.

Phase B — case → StructuredGeometry adapter
-------------------------------------------

The ``geometry_spec`` field carries both the *structural*
description (kind, regions, BCs) and the *truth* values
(``critical_dimension_mfp`` / ``critical_dimension_cm``). The Phase
B audit (``.claude/plans/dazzling-cuddling-boot.md``) identified
that this entanglement should be split: the truth values move onto
``La13511Truth`` (where they conceptually belong), and the geometry
itself is built on demand from the truth value via
:meth:`La13511Case.to_geometry`:

.. code-block:: python

   geom = case.to_geometry()
   # Returns StructuredGeometry with cm = mfp / Σ_t,
   # tag "SLB" / "SPH" / "CYL", and BCs:
   #   slab    → (BC.vacuum, BC.vacuum) — full-width 2a
   #   sphere  → (BC.vacuum,) — outer; centreline implicit reflective
   #   cylinder→ (BC.vacuum,) — outer; centreline implicit reflective

Slab convention: returns the FULL slab width (``2 *
critical_dimension_mfp / Σ_t``) with vacuum-vacuum BCs — the
production-natural convention. F_N's natural half-thickness is
recovered inside :class:`MomentSpace` from
``geom.domain_extent_cm / 2``. This wastes the slab's natural
half-symmetry; a future improvement would encode Sood symmetric
slabs as half-slabs with reflective+vacuum BCs (half the cells in
production solves at the same accuracy).

Infinite-medium cases (``geometry == "infinite"``) raise
``ValueError`` from :meth:`La13511Case.to_geometry` with a directive
message pointing the caller at :meth:`MomentSpace.solve_kinf` /
:func:`solve_homogeneous_infinite` — those don't have a geometry
to build.

Phase B is *additive*: legacy consumers reading
``case.geometry_spec.critical_dimension_mfp`` continue to work
unchanged. Phase F (after Phases C/D migrate consumers) deletes the
duplicate fields off ``GeometrySpec`` and the flat
``sood_table`` / ``primary_reference`` / ``notes`` triple off
``La13511Case``.

Cross-method consumer pattern
------------------------------

The same case is consumed two ways. **Production solvers** receive
the materials dict and a built mesh:

.. code-block:: python

   from orpheus.derivations.continuous.sood_registry import (
       LA13511_CASES, build_materials, build_mesh,
   )
   from orpheus.cp.solver import solve_cp

   case = LA13511_CASES["Ua-1-0-SL"]
   materials = build_materials(case)        # -> dict[int, Mixture]
   mesh = build_mesh(case, n_cells=64)      # -> Mesh1D
   result = solve_cp(materials, mesh, ...)
   # cross-check against case.truth.k_eff_or_kinf

**Semi-analytical reference solvers** extract numpy arrays:

.. code-block:: python

   from orpheus.derivations.continuous.sood_registry import (
       LA13511_CASES, mixture_to_fn_arrays,
   )
   from orpheus.derivations.continuous.fn_method.multi_group import (
       compute_kinf_mg,
   )

   case = LA13511_CASES["PU-2-0-IN"]
   sigma_t, sigma_s, nu_sigma_f, chi = mixture_to_fn_arrays(
       case.materials[0]
   )
   k_inf = compute_kinf_mg(sigma_t, sigma_s, nu_sigma_f, chi)
   # cross-check against case.truth.k_eff_or_kinf = 2.683767

The two paths above share the same ``case`` object — guaranteeing
the F_N truth value and the production CP solver are checked against
identical XS data.

ORPHEUS vs Sood group convention
=================================

Sood numbers groups :math:`g=N` (fast) → :math:`g=1` (slow), the
reverse of typical nuclear-engineering convention. ORPHEUS uses the
:ref:`canonical fast-first convention <canonical-group-convention>`
:math:`g=0` (fast) → :math:`g=N-1` (slow). The
:mod:`sood_registry.la13511` module does the conversion at
construction time so consumers see XS arrays in **ORPHEUS
convention** directly.

The scattering matrix uses ORPHEUS's ``[from, to]`` convention:
``sigma_s[g, h]`` is :math:`\Sigma_{s, g \to h}`. The
:math:`\Sigma_g^{\rm rem}` removal cross sections used by Sood
Eqs 26-27 are therefore ``sigma_t[g] - sigma_s[g, g]``.

The :func:`...fn_method.origins.k_inf_derivations.derive_*` SymPy
modules use Sood's symbols verbatim so equations match the LA-13511
report letter-for-letter; the conversion is purely a relabeling and
the algebra is identical for either side.

Geometry-spec build conventions
================================

.. _sood-registry-mesh-conventions:

The :meth:`GeometrySpec.build` method constructs a :class:`Mesh1D`
matching the published critical configuration. Conventions per
geometry:

.. list-table:: GeometrySpec.build conventions per geometry
   :header-rows: 1
   :widths: 14 30 28 28

   * - Geometry
     - Domain extent (cm)
     - Default BCs
     - Notes
   * - ``infinite``
     - n/a (raises)
     - n/a
     - Consumers ignore mesh; pass materials to transfer-matrix
       :math:`k_\infty` solver.
   * - ``slab``
     - ``2 * critical_dimension_cm``
     - vacuum at both ends
     - Builds the FULL symmetric slab :math:`[0, 2a]` where :math:`a`
       is the published critical half-thickness. F_N method's
       natural domain.
   * - ``sphere``
     - ``critical_dimension_cm``
     - reflective at :math:`r=0`; vacuum at outer
     - Builds :math:`[0, R]`.
   * - ``cylinder``
     - ``critical_dimension_cm``
     - reflective at :math:`r=0`; vacuum at outer
     - Builds :math:`[0, R]`. Same as sphere convention.
   * - ``ISLC``
     - n/a (raises NotImplementedError)
     - n/a
     - Reserved for future Initial-Slab-then-Layered-Cell cases.

Cross-method dimension conversion
----------------------------------

When the published critical dimension is the **half-thickness**
(slab) vs the **radius** (sphere/cylinder), method consumers each
have their own convention. The registry's
``critical_dimension_cm`` is the **as-published** value; the method
consumer is responsible for the correct interpretation:

.. list-table:: Critical-dimension conversions per consumer
   :header-rows: 1
   :widths: 15 28 28 29

   * - Geometry
     - Sood publishes
     - F_N consumes
     - Production solver consumes
   * - slab
     - half-thickness :math:`a`
     - :math:`a` (Eq 4 BC)
     - Mesh extent :math:`2a` (vacuum BC at ±a)
   * - sphere / cylinder
     - radius :math:`R`
     - :math:`R` (Eq 46 BC)
     - Mesh extent :math:`R` (reflective at 0, vacuum at R)

The registry stores the as-published value; consumers convert
via the unambiguous rules above. The
``test_sood_registry_compatibility.py`` foundation tests pin the
build conventions for every geometry kind.

First slice — five LA-13511 cases
==================================

.. _sood-registry-first-slice:

Phase A migrated the five first-slice cases verbatim from the
legacy ``fn_method.benchmarks.la13511`` module:

.. list-table::
   :header-rows: 1
   :widths: 15 15 25 25 20

   * - Case ID
     - Sood #
     - Material
     - Geometry
     - Truth
   * - ``PUa-1-0-IN``
     - 1
     - Pu-239 (a) 1G
     - infinite
     - :math:`k_\infty = 2.612903`
   * - ``PU-2-0-IN``
     - 44
     - Pu 2G no-upscatter
     - infinite
     - :math:`k_\infty = 2.683767`,
       :math:`\phi_2/\phi_1 = 0.675229`
   * - ``Ua-1-0-SL``
     - 12
     - U-235 (a) 1G
     - slab
     - :math:`a_c = 0.93772556` mfp
   * - ``Ua-1-0-CY``
     - 13
     - U-235 (a) 1G
     - cylinder
     - :math:`R_c = 1.72500292` mfp
   * - ``Ua-1-0-SP``
     - 14
     - U-235 (a) 1G
     - sphere
     - :math:`R_c = 2.4248249802` mfp

Per-case provenance + cross-method status:

* **PUa-1-0-IN**: Sood Eq 19 + 20; verified by V_fn1.1 + V_fn1.2 at
  the algebra level; verified by the Branch-2
  :func:`...fn_method.multi_group.compute_kinf_1g` numerically;
  cross-checked against
  :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous` to
  ≥ 12 digits.
* **PU-2-0-IN**: **Sood Eq 28 has a typo** (caught by V_fn2.1 SymPy
  derivation; the printed Eq 28 reduces to 2.862 instead of the
  published 2.684 at no-upscatter); the corrected form lives in
  the SymPy module, the Branch-2 ``compute_kinf_2g_general`` evaluates
  the corrected form. Sood Eq 29 (printed correctly) is verified via
  V_fn2.2; the flux ratio Eq 32 is verified via V_fn2.3.
* **Ua-1-0-SL**: KLL 1974 NSE 54 truth source. Slab F_N solver
  (:func:`...fn_method.slab.solve_fn_slab_bare_critical`,
  Siewert-Benoist 1979 + Grandjean-Siewert 1979) reproduces at
  ≤ 5e-6 absolute on :math:`a_c`. Cross-checked against Variant α
  slab at 5e-5.
* **Ua-1-0-CY**: WM-72 1973 truth source. Cylinder is **outside the
  F_N pillar** (Mitsis-style Wiener-Hopf is documented as
  non-convergent for bare cylinder); cylinder ships via
  :ref:`theory-singular-eigenfunction` (WM-72 Mitsis-WM Fredholm at
  3e-7 relative). Cross-checked against Variant α cylinder at 8.5e-6.
* **Ua-1-0-SP**: KLL 1974 + Siewert-Thomas 1986 truth sources. Sphere
  F_N solver
  (:func:`...fn_method.sphere.solve_fn_sphere_bare_critical`,
  Siewert-Thomas 1986 1G specialisation of the 2G F_N) reproduces at
  3.6e-8 absolute. Cross-checked against Variant α sphere at 4.2e-6.
  KLL Table VII flux ratios populated from KLL 1974 Table VII c=1.30
  row.

Phase B3 — wide enumeration coverage
=====================================

.. _sood-registry-phase-b3-wide:

Phase B3 expanded the registry from the 5-case Phase A first slice
to **42 LA-13511 cases** spanning every class the existing
``fn_method`` machinery (k_inf 1G/2G/MG, slab F_N, sphere F_N) can
solve TODAY, plus stubs for cases pending solver dispatch:

.. list-table:: Phase B3 wide-enumeration coverage matrix
   :header-rows: 1
   :widths: 35 12 18 35

   * - Class
     - Cases
     - Status
     - Solver / cross-check
   * - 1G k_inf infinite-medium
     - 12
     - active
     - :func:`compute_kinf_1g` (Sood Eq 19)
   * - 2G k_inf infinite-medium (no upscatter)
     - 6
     - active
     - :func:`compute_kinf_2g_no_upscatter` (Sood Eq 29)
   * - 2G k_inf infinite-medium (with upscatter)
     - 2
     - active
     - :func:`compute_kinf_2g_general` (Sood Eq 28 corrected)
   * - 3G k_inf infinite-medium
     - 1
     - active
     - :func:`compute_kinf_mg` (Sood Eq 76)
   * - 6G k_inf infinite-medium
     - 1
     - active
     - :func:`compute_kinf_mg` (handles upscatter)
   * - 1G slab bare-critical
     - 4
     - active
     - :func:`solve_fn_slab_bare_critical` (Siewert-Benoist /
       Grandjean-Siewert F_N at N=12)
   * - 1G sphere bare-critical
     - 3
     - active
     - :func:`solve_fn_sphere_bare_critical` (Siewert-Thomas
       F_N at N=10)
   * - 1G P_1 anisotropic bare-critical (slab + sphere)
     - 5
     - active
     - F_N + KLL flux reconstruction (anisotropic extension)
   * - 1G cylinder bare-critical
     - 3
     - **STUB**
     - **Out of pillar.** Mitsis-style Wiener-Hopf does not
       converge for bare cylinder (WM-72); cylinder shipped
       via :ref:`theory-singular-eigenfunction`.
   * - 2G slab bare-critical
     - 5
     - **STUB**
     - Pending Siewert-Thomas 1986 2G F_N (matrix Λ
       dispersion law)
   * - 2G sphere bare-critical
     - 5
     - **STUB**
     - Same — 2G F_N machinery deferred

Total active: **30 cases**; total stubs: **12 cases** (registered
truth values; activation pending solver dispatch).

The cylinder STUB is intentionally permanent — the F_N machinery
does not extend to bare cylinder, and the WM-72 Mitsis-WM Fredholm
solver in :ref:`theory-singular-eigenfunction` already reaches
3e-7 relative on every cylinder configuration. The 2G STUBs are
pending Siewert-Thomas 1986 2G F_N (matrix Λ dispersion law); the
geometry-sign abstraction that lets slab and sphere share machinery
is already in place at the 1G level (see ``fn_method/core/fn_matrix.py``
``assemble_fn_matrix``). The 2G extension makes Λ a 2×2 matrix, the
Case eigenvalues 2×2 matrix roots, and the F_N matrix entries
matrix-valued — the structural pattern is the same.

Verification gates live in:

* :mod:`tests.derivations.test_sood_registry_wide_kinf` — k_inf gate,
  one parametrised case per :data:`WIDE_SLICE_KINF` entry plus the
  cross-implementation gate against
  :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`.
* :mod:`tests.derivations.test_sood_registry_wide_bare_critical` —
  L1 reference-value gate for slab/sphere F_N.

Tolerances achieved:

* Every k_inf case at ≤ 7e-6 absolute (well within 1e-5 target).
* Slab F_N at N=12: ≤ 3e-6 absolute on :math:`a_c`.
* Sphere F_N at N=10: ≤ 5e-8 absolute on :math:`R_c`.

Atalay 1997 reflected/anisotropic registry
============================================

.. _sood-registry-atalay-cases:

The Atalay 1997 case catalogue lives in
:mod:`orpheus.derivations.continuous.sood_registry.atalay1997`. It is
**method-agnostic** in the same sense as the LA-13511 catalogue —
the case dataclass is :class:`La13511Case` (re-used) carrying
``materials`` + ``geometry_spec`` + ``truth`` — but parametrises
over :math:`(c, R, f_1)` triples rather than material composition,
mirroring how Atalay published the cases.

Cases shipped:

.. list-table:: Atalay 1997 cases
   :header-rows: 1
   :widths: 30 10 10 10 40

   * - Case ID
     - :math:`c`
     - :math:`R`
     - :math:`f_1`
     - Truth (mfp)
   * - ``atalay-1997-slab-c1.30-R0.00-f1_0.00``
     - 1.30
     - 0.00
     - 0.00
     - :math:`2d = 1.87766`
   * - ``atalay-1997-slab-c1.30-R0.25-f1_0.00``
     - 1.30
     - 0.25
     - 0.00
     - :math:`2d = 1.40621`
   * - ``atalay-1997-slab-c1.30-R0.50-f1_0.00``
     - 1.30
     - 0.50
     - 0.00
     - :math:`2d = 0.89317`
   * - ``atalay-1997-slab-c1.30-R0.75-f1_0.00``
     - 1.30
     - 0.75
     - 0.00
     - :math:`2d = 0.40758`
   * - ``atalay-1997-slab-c1.30-R0.00-f1_0.10``
     - 1.30
     - 0.00
     - 0.10
     - :math:`2d = 1.94146`
   * - ``atalay-1997-slab-c1.30-R0.50-f1_0.10``
     - 1.30
     - 0.50
     - 0.10
     - :math:`2d = 0.89831`
   * - ``atalay-1997-sphere-c1.30-R0.00-f1_0.00``
     - 1.30
     - 0.00
     - 0.00
     - :math:`R_c = 2.4248249802`

Atalay 1997 is the **primary source** for the reflected +
linearly-anisotropic cross-product cases that lie outside both the
Sood/Forster/Parsons LA-13511 truth set (which focuses on bare
configurations) and the Burkart-Ishiguro-Siewert 1976 F_N reference
(vacuum-only).

The Atalay catalogue is consumed by the Atalay slab/sphere solvers
in :mod:`...singular_eigenfunction.slab` and
:mod:`...singular_eigenfunction.sphere`. Tolerances achieved are
documented in :ref:`theory-se-precision-floor` (the ERR-038 paper
floor at small slab thicknesses, characterised but not a code bug).

Production-protocol bridge tests
=================================

.. _sood-registry-bridge-tests:

The structural-bridge gates live in
:mod:`tests.derivations.test_sood_registry_compatibility` and verify:

* Every registered case has a well-formed
  ``materials: dict[int, Mixture]`` (foundation).
* The :func:`mixture_to_fn_arrays` extractor returns Sood's published
  XS values exactly (foundation, parametrised over 1G and 2G cases).
* Legacy property accessors (``case.sigma_t``, ``case.sigma_s``,
  ``case.nu_sigma_f``, ``case.chi``) match the extractor — load-
  bearing back-compat with existing F_N tests (foundation).
* :class:`GeometrySpec` build conventions hold for every geometry
  kind (foundation).
* :func:`kinf_homogeneous` consumes the registry directly and
  reproduces Sood truth at both 1G and 2G cases to ≤ 1e-5
  (foundation).
* ``solve_cp`` and ``solve_sn`` accept the registry's ``materials``
  + ``mesh`` without raising (L1 smoke; correctness cross-checks
  land in Phase B).

Why these gates are foundation-tagged and not L1
-------------------------------------------------

The compatibility tests verify **software invariants of the
registry / extractor / build interface** — not solver-correctness
claims about the underlying physics. Examples:

* "The registry has a well-formed ``materials`` dict for every
  case" is a software invariant of the schema. If this breaks, no
  consumer downstream works; it's a foundational property of the
  registry's API contract.
* "The XS extractor returns the published Sood values exactly" is
  a software invariant of the extractor function — it must be
  bit-for-bit because consumers depend on the extractor not
  changing the values silently.
* "GeometrySpec.build returns a Mesh1D with the right
  geometry" is a software invariant of the build interface.

These are all **foundation-tagged** because they don't claim
anything about the *physics* of the cases — they claim only that
the *plumbing* of the registry works.

The L1 / L2 numerical claims (e.g., "F_N slab solver reproduces
:math:`a_c` to 5e-6") live in the consumer test files (e.g.,
``test_fn_la13511_slab.py``, ``test_sood_registry_wide_bare_critical.py``).
Each consumer is responsible for verifying its own physics; the
registry's job is solely to deliver consistent input data.

See :doc:`/testing/architecture` for the V&V level taxonomy
(foundation / L0 / L1 / L2 / L3 / L4) and how each level applies.

Solver-output caching
=====================

.. _sood-registry-cache:

The :mod:`~orpheus.derivations.continuous.sood_registry.cache`
submodule provides a **persistent solver-output cache** for tests
that exercise slow solvers on Sood cases. It is distinct from the
*truth-value* "cache" — published reference values are hard-coded
on each :class:`La13511Case` and need no caching infrastructure.
The cache addresses the complementary problem: when a test calls
e.g. ``solve_fn_sphere_bare_critical(c=1.30, n_modes=15)``, store
the result on disk so subsequent test runs reuse it instead of
recomputing.

Key design choices
-------------------

* **Cache key** — ``(solver_qualified_name, params_hash)`` where
  ``params_hash`` is a 16-char SHA-256 prefix over a canonicalized
  JSON representation of the solver kwargs. Numpy scalars/arrays
  are normalised so ``np.float64(1.30)`` and ``1.30`` produce the
  same key.
* **Versioning** — every cache entry carries a version string. Default
  ``"auto"`` resolves to ``git rev-parse --short HEAD``; a mismatch
  on read is treated as a miss (and the entry is overwritten on the
  next write). This auto-invalidates entries on solver-implementation
  changes that bump the SHA. Callers may pin an explicit string
  instead.
* **Storage** — pickle files under ``.cache/sood_registry/`` at the
  project root, one ``.pkl`` per entry. The directory is gitignored.
  Pickle is chosen over JSON because solver result types are
  dataclasses with numpy arrays — pickle handles them natively.
  Inspecting an entry from a debugger is a single ``pickle.load``
  call.
* **Local-machine state** — nothing about the cache directory should
  ever be committed; the cache lives in ``.cache/sood_registry/``
  per machine.

API
----

Two equivalent forms:

**Decorator** (recommended for one-off test wrappers):

.. code-block:: python

   from orpheus.derivations.continuous.sood_registry import sood_cache

   @sood_cache()
   def fn_sphere_at_n(*, c, n_modes):
       return solve_fn_sphere_bare_critical(c=c, n_modes=n_modes)

   res = fn_sphere_at_n(c=1.30, n_modes=15)  # first call: computes
   res = fn_sphere_at_n(c=1.30, n_modes=15)  # second call: cache hit

**Class** (for tests that need programmatic invalidation /
introspection):

.. code-block:: python

   from orpheus.derivations.continuous.sood_registry import (
       SoodResultCache,
   )

   cache = SoodResultCache()
   result = cache.get_or_compute(
       solver_name="solve_fn_sphere",
       params={"c": 1.30, "n_modes": 15},
       compute=lambda: solve_fn_sphere_bare_critical(c=1.30, n_modes=15),
   )

Where the cache pays off
-------------------------

The cache is **opt-in**: no existing test imports it. Adoption is
per-test, on demand, by the test author. The cache pays off most
where:

* **Sphere F_N at large N** (:math:`N \ge 15`): the prominence-
  filtered first-local-minimum bracket scan + Brent refinement is
  ~10× slower at N=15 than at N=10. Caching at N=15 across the 5+
  test cases that need the high-N flux reconstruction saves
  several minutes per test session.
* **Variant α at fine quadrature**: the Phase 3 cylinder /
  hollow-sphere benchmarks at :math:`(n_r, n_\mu) = (24, 128)` take
  ~30 s per solve. Caching at this resolution is the difference
  between a 1-minute regression suite and a 5-minute one.
* **Future cylinder Westfall-Metcalf F_N**: when the Mitsis-WM
  Fredholm solver is dispatched on the cylinder STUB cases, those
  solves will be ~1 s each at the standard
  :math:`n_{\rm grid} = 24`. Caching is overkill at that speed,
  but the cache infrastructure is in place for future hot paths.

**Cache poisoning hazard** — if a solver has a bug and the cache
stores a wrong answer, the wrong answer survives across runs until
either the version bumps (auto-invalidation on git SHA change) or
the user calls :func:`clear_cache`. **Mitigation**: the test that
consumes the cached output is itself the protection — if the same
test asserts against an external reference (``assert |result -
truth| < tol``), a poisoned cache fails the test exactly as
poison-detection would dictate.

The cache infrastructure is pinned by
:mod:`tests.derivations.test_sood_registry_cache` —
load-bearing invariants: round trip, miss-vs-hit, version
invalidation, hash stability, ``clear()``, decorator integration,
and an L1 smoke against the transfer-matrix :math:`k_\infty`
reference solver on Sood ``PUa-1-0-IN``.

Forward-looking — Wave 3 meta-registry
=========================================

.. _sood-registry-wave-3-foreshadowing:

The current registry has two registries side-by-side:
:class:`LA13511_CASES` (Sood/Forster/Parsons 1999) and the Atalay
catalogue (Atalay 1997). Both share the same case schema
(:class:`La13511Case` is re-used for both, despite the misleading
name) and both feed the same consumer adapter functions
(:func:`mixture_to_fn_arrays`, :func:`build_materials`,
:func:`build_mesh`).

The architectural pattern this surfaces — **method-agnostic case
collections, registered separately, consumed identically** — is
likely to extend in future Waves to:

* **Kaper-Lindeman-Leaf 1974 (KLL)** — interior flux reconstruction
  Tables III + VII (slab + sphere) currently consumed by
  :func:`...fn_method.flux_reconstruction.*` directly; could be
  promoted to a sister registry for consistency with the LA-13511
  + Atalay pattern.
* **Westfall-Metcalf 1973 (WM-72)** — Table II cylinder critical
  radii (six configurations); currently inlined in
  :mod:`tests.derivations.test_singular_eigenfunction_cylinder`;
  could be promoted similarly.
* **Burkart-Ishiguro-Siewert 1976 (BIS)** — two-region anisotropic
  F_N reference values; needed for the future reflected-slab P_N
  expansion.
* **Garcia 2021** — multi-region fixed-source benchmarks (Table 5
  ppP_N); already consumed by
  :mod:`tests.derivations.test_peierls_greens_function_garcia2021`.

The Wave 3 architectural plan
(``.claude/plans/wave3/architecture.md``) discusses promoting these
to a meta-registry where the schema is generalised across all
benchmark sources, but the timing is **deliberately deferred** —
the discipline is "build each new geometry/topology standalone first;
only unify after ≥ 2 working instances" (see
``.claude/agent-memory/feedback_unify_after_two_instances.md``).
The current LA-13511 + Atalay split is the necessary first step;
extension to KLL / WM-72 / BIS / Garcia would require ≥ 4 sister
registries before unification makes architectural sense.

For now: each new published-reference case collection should
follow the LA-13511 / Atalay pattern (drop a new module under
``sood_registry/``, share the case schema, register cases as module-
level constants). When the meta-registry pattern is appropriate, it
will be backwards-compatible with the existing case files.

References
==========

* Sood, A., Forster, R. A., & Parsons, D. K. (1999).
  *Analytical Benchmark Test Set for Criticality Code Verification*,
  LA-13511, Los Alamos National Laboratory. PDF locally available
  at ``scratch/literature/Sood Foster Parsons (1999)Analytical
  Benchmark Test Set for Criticality Code Verification.pdf``.
* Sood, A., Forster, R. A., & Parsons, D. K. (2003).
  *Analytical Benchmark Test Set for Criticality Code Verification.*
  *Progress in Nuclear Energy* **42**, 55.
  Journal-published condensation; verified to be a TEST SET not a
  method paper (see
  ``.claude/agent-memory/literature-researcher/sood_2003_vs_1999_extraction.md``).
* Atalay, M.A. (1997).
  *The reflected slab and sphere criticality problem with
  anisotropic scattering in one-speed neutron transport theory.*
  *Progress in Nuclear Energy* **31**\ (3), 229-252.
  DOI: 10.1016/0149-1970(95)00094-1. PDF locally available at
  ``scratch/literature/Atalay(1997)the reflected slab and sphere
  criticality problem with anisotropic scattering in one-speed
  neutron transport theory.pdf``.
* Kaper, H. G., Lindeman, A. J., & Leaf, G. K. (1974).
  *Benchmark Values for the Slab and Sphere Criticality Problem
  in One-Group Neutron Transport Theory.* *Nuclear Science and
  Engineering* **54**, 94. PDF locally available at
  ``scratch/literature/Kaper-Lindeman-Leaf(1974)Benchmark Values
  for the Slab and Sphere Criticality Problem in One-Group Neutron
  Transport Th.pdf``.
* Westfall, R. M., & Metcalf, J. R. (1973).
  *Singular Eigenfunction Solution of the Monoenergetic Neutron
  Transport Equation for Finite Radially Reflected Critical
  Cylinders.* *Nuclear Science and Engineering* **52**, 1.
* Siewert, C.E., Benoist, P. (1979).
  *The F_N Method in Neutron-Transport Theory. Part I.*
  *Nuclear Science and Engineering* **69**, 156.
* Grandjean, P., Siewert, C.E. (1979).
  *The F_N Method in Neutron-Transport Theory. Part II.*
  *Nuclear Science and Engineering* **69**, 161.
* Siewert, C.E., Thomas, J.R. (1986).
  *On Two-Group Critical Problems in Neutron Transport Theory.*
  *Nuclear Science and Engineering* **94**, 264-270.
* Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976).
  *Neutron transport in two dissimilar media with anisotropic
  scattering.* *Nuclear Science and Engineering* **61**, 72-81.

Internal references:

* Phase A migration closeout:
  ``.claude/agent-memory/method-implementer/sood_registry_phase_a_migration.md``.
* Phase B3 wide-enumeration closeout:
  ``.claude/agent-memory/method-implementer/sood_registry_wide_enumeration_phase_b3.md``.
* Cache infrastructure closeout:
  ``.claude/agent-memory/method-implementer/sood_registry_cache_phase_b4.md``.
* Wave 3 architectural plan:
  ``.claude/plans/wave3/architecture.md``.
* :doc:`fn_method` — primary consumer of the LA-13511 catalogue
  (k_inf cases + slab/sphere bare-critical + reflected slab).
* :doc:`singular_eigenfunction` — primary consumer of the Atalay
  catalogue + cylinder LA-13511 truth values.
* :doc:`trajectory_resolvent` — Variant α cross-checks on shared
  Sood truth values.
