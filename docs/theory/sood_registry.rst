.. _theory-sood-registry:

==========================================================================
Sood-family benchmark case registry — method-agnostic truth set
==========================================================================

.. note::

   **Stub-grade theory page.** This page documents the
   :mod:`orpheus.derivations.continuous.sood_registry` package — the
   single source of truth for benchmark case configurations from the
   Sood-family literature (LA-13511 first slice; LA-UR-03-1987
   cylinder benchmarks and KLL Tables to follow in Phase B).

   The archivist agent will expand each TODO marker into the full
   narrative.

Purpose
========

.. _sood-registry-purpose:

The :mod:`sood_registry` package decouples **case configurations**
(cross sections, geometry, published reference values) from
**solver methods**. Before this package, each method (F_N, Variant α,
etc.) carried its own per-case data alongside its solvers — coupling
that made the same Sood case awkward to consume from multiple
methods. Now a single :class:`La13511Case` is consumed by:

* Semi-analytical reference solvers (F_N method, PS-1982 wrapper,
  transfer-matrix :math:`k_\infty`) — extract numpy XS arrays via
  :func:`mixture_to_fn_arrays`.
* Production discrete solvers (CP, SN, MOC) — accept the
  ``materials: dict[int, Mixture]`` + ``mesh_template.build()``
  directly.

.. note:: TODO — Archivist expansion needed.

   The full design rationale for the method-agnostic split lives in
   the closeout memo
   ``.claude/agent-memory/method-implementer/sood_registry_phase_a_migration.md``.
   Topics: why method coupling was wrong factoring; how the
   production protocol (Mixture + Mesh1D) was already in place but
   unused for cases; what changed in fn_method/benchmarks and how
   the deprecation shim works; how Phase B (Westfall-Metcalf cylinder
   F_N + flux reconstruction + wide enumeration + cache) builds on
   this scaffold.

Schema
======

.. _sood-registry-schema:

The :class:`La13511Case` dataclass carries:

* **Identity** — ``case_id``, ``problem_number``, ``description``.
* **Production-protocol XS** — ``materials: dict[int, Mixture]``,
  keyed by integer material ID. Single-region cases use
  ``{0: Mixture(...)}``; multi-region cases (none yet in Phase A)
  add more keys.
* **Mesh template** — :class:`MeshTemplate` carrying geometry kind
  (``"infinite"``, ``"slab"``, ``"sphere"``, ``"cylinder"``,
  ``"ISLC"``), the published critical dimension in mfp + cm, and
  per-end boundary conditions. The ``build(n_cells)`` method
  produces a :class:`Mesh1D` at the requested refinement.
* **Truth** — :class:`La13511Truth` bundling
  :math:`k_{\rm eff}` / :math:`k_\infty`, flux ratios, and angular
  flux references (whatever the published reference tabulates).
* **Provenance** — Sood table number, primary reference, free-form
  notes.

.. note:: TODO — Archivist expansion needed.

   Walk through the schema field-by-field with code snippets showing
   each field's role. Show concretely how the same case is consumed
   by F_N (via :func:`mixture_to_fn_arrays`) and by production solvers
   (via :func:`build_materials` + :func:`build_mesh`). Include the
   ORPHEUS-vs-Sood group convention conversion that lives at
   construction time. Reference the SymPy module
   :mod:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations`
   for the algebraic identities each :math:`k_\infty` reference value
   is derived from. The implementation is in
   :mod:`orpheus.derivations.continuous.sood_registry.la13511`.

Mesh-template build conventions
================================

.. _sood-registry-mesh-conventions:

The :meth:`MeshTemplate.build` method constructs a :class:`Mesh1D`
matching the published critical configuration:

* **Slab** (``geometry == "slab"``): builds the FULL symmetric slab
  ``[0, 2a]`` where :math:`a` is the published critical
  half-thickness. Default BCs: vacuum at both ends. This matches
  the F_N method's natural domain.
* **Sphere / cylinder**: builds ``[0, R]`` where :math:`R` is the
  critical radius. Default BCs: reflective at :math:`r = 0`
  (centreline / axis), vacuum at the outer surface.
* **Infinite**: no mesh is constructed; consumers should ignore
  ``mesh_template`` and pass ``materials`` to a transfer-matrix
  reference solver.

.. note:: TODO — Archivist expansion needed.

   Document each geometry's BC conventions in detail; include the
   F_N-vs-Variant-α-vs-CP/SN conversion table for "what does the
   slab full-width cm mean for each method". Reference
   ``tests/derivations/test_sood_registry_compatibility.py``
   (test_slab_mesh_builds_full_symmetric_domain etc.) as the load-
   bearing invariants.

First slice — five LA-13511 cases
==================================

.. _sood-registry-first-slice:

Phase A migrates the five first-slice cases verbatim from the
legacy ``fn_method.benchmarks.la13511`` module:

.. list-table::
   :header-rows: 1
   :widths: 15 15 25 25 20

   * - Case ID
     - Sood #
     - Material
     - Geometry
     - Truth
   * - PUa-1-0-IN
     - 1
     - Pu-239 (a) 1G
     - infinite
     - :math:`k_\infty = 2.612903`
   * - PU-2-0-IN
     - 44
     - Pu 2G
     - infinite
     - :math:`k_\infty = 2.683767`
   * - Ua-1-0-SL
     - 12
     - U-235 (a) 1G
     - slab
     - :math:`a_c = 0.93772556` mfp
   * - Ua-1-0-CY
     - 13
     - U-235 (a) 1G
     - cylinder
     - :math:`R_c = 1.72500292` mfp
   * - Ua-1-0-SP
     - 14
     - U-235 (a) 1G
     - sphere
     - :math:`R_c = 2.4248249802` mfp

.. note:: TODO — Archivist expansion needed.

   For each case, document the Sood XS verbatim, show the SymPy
   identity that yields the published truth (for k_inf cases), and
   reference the F_N test that pins the case at the published
   precision. Include cross-method status: Ua-1-0-SP is
   structurally cross-checked at 4.2e-6 by F_N vs Variant α;
   Ua-1-0-CY is cross-checked at 8.5e-6 by Variant α (F_N
   cylinder pending Westfall-Metcalf, Phase B).

Production-protocol bridge tests
=================================

.. _sood-registry-bridge-tests:

The structural-bridge gates live in
``tests/derivations/test_sood_registry_compatibility.py`` and verify:

* Every registered case has a well-formed
  ``materials: dict[int, Mixture]`` (foundation).
* The :func:`mixture_to_fn_arrays` extractor returns Sood's published
  XS values exactly (foundation, parametrised over 1G and 2G cases).
* Legacy property accessors (``case.sigma_t``, ``case.sigma_s``,
  ``case.nu_sigma_f``, ``case.chi``) match the extractor — load-
  bearing back-compat with existing F_N tests (foundation).
* :class:`MeshTemplate` build conventions hold for every geometry
  kind (foundation).
* :func:`kinf_homogeneous` consumes the registry directly and
  reproduces Sood truth at both 1G and 2G cases to ≤ 1e-5
  (foundation).
* ``solve_cp`` and ``solve_sn`` accept the registry's
  ``materials`` + ``mesh`` without raising (L1 smoke; correctness
  cross-checks land in Phase B).

.. note:: TODO — Archivist expansion needed.

   Explain why these gates are foundation-tagged (software
   invariants of the registry / extractor / build interface, not
   solver-correctness claims) vs L1-tagged (the production-solver
   bridge smoke is one rung above foundation because it crosses a
   call boundary into a numerical solver). Reference the V&V
   architecture document
   (``docs/testing/architecture.rst``) for the full ladder.

Phase B preview
================

.. _sood-registry-phase-b:

Phase B (parallel after Phase A lands) will:

1. Build the cylinder Westfall-Metcalf F_N solver and activate
   ``Ua-1-0-CY`` against Sood's published truth.
2. Add flux reconstruction support so the slab/sphere F_N solvers
   can verify the published flux ratios at
   ``r/r_c ∈ {0.25, 0.5, 0.75, 1.0}`` (currently only the critical
   radius is verified).
3. Enumerate the wide LA-13511 catalogue (problems 1-67) into the
   registry — currently only the 5-case first slice is populated.
4. Build a verification cache so re-running the LA-13511 sweep
   doesn't recompute reference solver values from scratch.

.. note:: TODO — Archivist expansion needed.

   Once Phase B lands, this section becomes the "wide enumeration"
   chapter. For now it is a forward-looking placeholder.

References
==========

* Sood, A., Forster, R. A., & Parsons, D. K. (1999).
  *Analytical Benchmark Test Set for Criticality Code Verification*,
  LA-13511, Los Alamos National Laboratory.
* Kaper, H. G., Lindeman, A. J., & Leaf, G. K. (1974).
  Benchmark values for the slab and sphere criticality problem in
  one-group neutron transport theory.
  *Nucl. Sci. Eng.* **54**, 94.
* Westfall, R. M., & Metcalf, J. R. (1973).
  Two-region critical benchmark calculations.
  *Nucl. Sci. Eng.* **52**, 1.
