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

Key components:

* **Cache key** — ``(solver_qualified_name, params_hash)`` where
  ``params_hash`` is a 16-char SHA-256 prefix over a canonicalized
  JSON representation of the solver kwargs. Numpy scalars/arrays
  are normalised so ``np.float64(1.30)`` and ``1.30`` produce the
  same key.
* **Version** — every cache entry carries a version string. Default
  ``"auto"`` resolves to ``git rev-parse --short HEAD``; a
  mismatch on read is treated as a miss (and the entry is
  overwritten on the next write). This auto-invalidates entries
  on solver-implementation changes that bump the SHA.
* **Storage** — pickle files under ``.cache/sood_registry/`` at
  the project root, one ``.pkl`` per entry. The directory is
  gitignored. Pickle is chosen over JSON because solver result
  types are dataclasses with numpy arrays — pickle handles them
  natively. Inspecting an entry from a debugger is a single
  ``pickle.load`` call.
* **API** — two equivalent forms:

  - **Decorator** (recommended for one-off test wrappers):

    .. code-block:: python

       from orpheus.derivations.continuous.sood_registry import sood_cache

       @sood_cache()
       def fn_sphere_at_n(*, c, n_modes):
           return solve_fn_sphere_bare_critical(c=c, n_modes=n_modes)

  - **Class** (for tests that need programmatic invalidation /
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

The cache is **opt-in**: no existing test imports it. Adoption is
per-test, on demand, by the test author. The unit-test gate at
``tests/derivations/test_sood_registry_cache.py`` pins the
load-bearing invariants — round trip, miss-vs-hit, version
invalidation, hash stability, ``clear()``, decorator integration,
and an L1 smoke against the transfer-matrix :math:`k_\infty`
reference solver on Sood ``PUa-1-0-IN``.

.. note:: TODO — Archivist expansion needed.

   Document where caching pays off most: sphere F_N at
   :math:`N \geq 15` (~10x speedup expected); Variant α at fine
   quadrature (Phase 3 cylinder / hollow-sphere benchmarks at
   ``n_r=24, n_mu=128``); future cylinder Westfall-Metcalf F_N
   when Phase B1 lands. Discuss the cache-poisoning hazard
   (a buggy solver's output survives across runs until the SHA
   changes or the user clears manually) and the mitigation:
   external-reference assertions in the same test that would
   themselves fail under poisoning. Reference the closeout memo
   ``.claude/agent-memory/method-implementer/sood_registry_cache_phase_b4.md``
   and the Phase B1/B2/B3 closeout memos for context on the
   parallel work.

Phase B3 — wide enumeration coverage
=====================================

.. _sood-registry-phase-b3-wide:

Phase B3 expanded the registry from the 5-case Phase A first slice to
**42 cases** spanning every LA-13511 class the existing ``fn_method``
machinery (k_inf 1G/2G/MG, slab F_N, sphere F_N) can solve TODAY,
plus stubs for the cases pending solver dispatch. Coverage matrix:

.. list-table::
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
     - :func:`solve_fn_sphere_bare_critical` (Siewert-Thomas F_N at N=10)
   * - 1G cylinder bare-critical
     - 3
     - **STUB**
     - Pending B1 dispatch (Westfall-Metcalf 1973 NSE 52, 1)
   * - 2G slab bare-critical
     - 5
     - **STUB**
     - Pending Siewert-Thomas 1986 2G F_N (matrix Λ dispersion law)
   * - 2G sphere bare-critical
     - 5
     - **STUB**
     - Same — 2G F_N machinery deferred

Total active: **30 cases**; total stubs: **12 cases** (registered
truth values; activation pending solver dispatch).

Verification gates live in:

* ``tests/derivations/test_sood_registry_wide_kinf.py`` — k_inf gate,
  one parametrised case per :data:`WIDE_SLICE_KINF` entry plus the
  cross-implementation gate against
  :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`.
* ``tests/derivations/test_sood_registry_wide_bare_critical.py`` —
  L1 reference-value gate for slab/sphere F_N.

Tolerances achieved: every k_inf case at ≤ 7e-6 absolute (well within
1e-5 target). Slab F_N at N=12 reaches ≤ 3e-6 absolute on
:math:`a_c`; sphere F_N at N=10 reaches ≤ 5e-8 absolute on
:math:`R_c`.

.. note:: TODO — Archivist expansion needed.

   For each Phase-B3 class, document Sood's published cross-section
   set, which equation (Sood Eq 19 / 28 / 29 / 32 / 42 / 59 / 76) gives
   the closed-form k_inf, and how the SymPy proof in
   :mod:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations`
   verifies the algebraic identity. Highlight the URRb/URRc upscatter
   distinction (Σ_21s > 0 — must use Eq 28 general form, not Eq 29).
   For URR-3-0-IN, walk through Forster's worked example (Eqs 60-65 —
   cross sections constructed by inverse design so f_23=4, f_13=15
   give k_inf=1.60 exactly).

   For the 2G bare-critical stubs, document the Siewert-Thomas 1986
   matrix Λ dispersion law and the geometry-sign abstraction that lets
   slab and sphere share machinery (already in place at the 1G level
   — see ``fn_method/core/fn_matrix.py::assemble_fn_matrix``). The
   2G extension makes Λ a 2x2 matrix, the Case eigenvalues are 2x2
   matrix roots, and the F_N matrix entries become matrix-valued —
   the structural pattern is the same.

   Reference the Phase B3 closeout memo
   ``.claude/agent-memory/method-implementer/sood_registry_wide_enumeration_phase_b3.md``
   for the full enumeration log + tolerance table.

Phase B preview (forward-looking)
==================================

.. _sood-registry-phase-b:

Remaining Phase B work after the parallel B1/B2/B3/B4 dispatches:

1. **B1 active**: cylinder Westfall-Metcalf F_N solver. Activates the
   3 cylinder STUB cases (``Ua-1-0-CY`` already cross-checked at
   8.5e-6 by Variant α; B1 will provide a structurally-independent
   reference).
2. **B2 active**: flux reconstruction so the slab/sphere F_N solvers
   can verify the published flux ratios at
   :math:`r/r_c \in \{0.25, 0.5, 0.75, 1.0\}` against
   :attr:`La13511Truth.flux_ratios` (currently only the critical
   radius is verified).
3. **Siewert-Thomas 1986 2G F_N** (separate dispatch, future): the
   load-bearing follow-on after B1/B2 land. Activates 10+ 2G
   bare-critical STUBS (PU-2-0-SL/SP, U-2-0-SL/SP, UAL-2-0-SL/SP,
   URRa-2-0-SL/SP, UD2O-2-0-SL/SP). The 2G sphere F_N is the natural
   structurally-independent reference for the Variant α sphere
   2G/MG extension.

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
