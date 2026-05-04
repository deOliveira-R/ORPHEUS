Cross-Method Regression Protocol
=================================

.. contents::
   :local:
   :depth: 2

Motivation
----------

ORPHEUS ships **multiple continuous-reference solvers** that target
the same physical problems via genuinely different mathematical
frames:

* :mod:`orpheus.derivations.continuous.fn_method` — Case
  singular-eigenfunction representation; collocation system on
  μ-moments; Wiener-Hopf factorisation.
* :mod:`orpheus.derivations.continuous.trajectory_resolvent` —
  bouncing-trajectory angle-resolved Green's function; fixed-point
  iteration on surface inflow.
* :mod:`orpheus.derivations.continuous.peierls_nystrom` — Peierls
  integral-equation kernel; Nyström quadrature with singularity
  subtraction.
* :mod:`orpheus.derivations.continuous.singular_eigenfunction` —
  Mitsis-WM Fredholm method (cylinder-specific).
* (planned) ``spectral_resolvent``, ``galerkin_spectral`` per the
  taxonomy in ``.claude/scratch/folder_naming_taxonomy.md``.

When a method is being **architecturally refactored** (e.g. the
trajectory_resolvent overhaul that motivated this protocol), we need
a regression net that lights up if the refactor breaks numerical
agreement with **structurally-independent** references. The
existing per-method test files already pin each method against its
own analytical-truth source (KLL, Sood, Grandjean-Siewert, NM 1980),
but they were silos — adding one cross-method case required
copy-pasting ~50 lines across two files.

The cross-method protocol fixes that. It defines:

1. an abstract :class:`~tests.cross_method.protocol.CrossMethodCase`
   that bundles a registry case (``La13511Case`` etc.) with
   per-solver tolerances, claim layer, and verification pillar;
2. a :class:`~tests.cross_method.protocol.SolverAdapter` Protocol
   that every solver wraps to (``solve(case) -> ScalarResult``);
3. populated case sets in :mod:`tests.cross_method.cases` covering
   bare-critical slab/sphere, reflected slab, and closed-sphere
   k_inf;
4. parametrised test gates in
   :mod:`tests.cross_method.test_eigenvalue` running every
   (case × adapter) pair plus pairwise agreement gates.

Adding a new solver to the regression net is now ≤ 50 lines: write
the adapter, register it in
:data:`~tests.cross_method.adapters.ADAPTERS_BY_NAME`, opt the
relevant cases in via their ``tolerances`` map.

V&V level mapping
-----------------

Per :doc:`/skills/vv-principles` §"V&V level taxonomy", the
conceptual level of "two solvers agreeing" is **L4** —
"informational, parallel to the ladder; produces zero correctness
info on its own". The L0-L3 ladder is reserved for evidence that
proves *correctness*, not *agreement*.

In practice the ORPHEUS test harness exposes L0..L3 + foundation as
markers; L4 is not registered. The convention in shipped
cross-method gates
(:mod:`tests.derivations.test_fn_la13511_slab_xverif`,
:mod:`tests.derivations.test_fn_la13511_sphere_xverif`) is to tag
them as **L1** because:

* each individual method is L1-verified against analytical truth in
  its own per-method test file;
* the cross-method comparison is a *consumer* of those L1
  references — not a new claim;
* when both methods are L1-grade and structurally independent (e.g.
  F_N's Case-eigenfunction representation vs trajectory_resolvent's
  bouncing-characteristic representation), their pairwise agreement
  is **L1-strength evidence** for either, not just L4
  "implementations agree".

The cross-method gates inherit L1 from this convention. The
``test_*_matches_truth`` gates carry L1 as the truth match. The
``test_fn_*_vs_trajectory_resolvent_*`` gates also carry L1 — the
agreement is L1-strength when both methods' L1 backing exists.

Foundation gates — schema invariants on the protocol metadata —
carry ``foundation``.

Reference contamination — the agreement-tolerance discipline
-----------------------------------------------------------------

A cross-method gate that pairs two L1-verified methods MUST set its
agreement tolerance to **the larger of the two truth tolerances**.
Tighter is the canonical reference-contamination anti-pattern (per
:doc:`/skills/vv-principles` §6 AI failure mode #6). Both methods
reach the truth value to within their respective floors; their
pairwise agreement therefore CANNOT be tighter than the looser of
the two without one method being calibrated against the other (the
contamination mechanism).

The :func:`~tests.cross_method.protocol.agreement_tolerance` helper
implements this rule:

.. code-block:: python

   def agreement_tolerance(case, a, b):
       return max(case.tolerance_for(a), case.tolerance_for(b))

If the actual cross-method agreement is much tighter than this
upper bound, that's interesting numerical evidence (and may mean
both adapters' truth tolerances are pessimistic) but the gate
declares its tolerance at the conservative bound.

Pillar tags — what each truth supports
--------------------------------------

Every :class:`~tests.cross_method.protocol.CrossMethodCase` carries
a ``pillar`` field per
:doc:`/skills/vv-principles` §"The three pillars of verification":

* ``"closed-form"`` — KLL/Sood transcendental dispersion roots,
  Grandjean-Siewert Table XI critical thicknesses, NM 1980
  Burkart-1976 'Exact' values. The bare-critical slab/sphere c-
  sweep is closed-form.
* ``"semi-analytical"`` — PS-1982 Peierls integral-equation
  reference (mpmath quadrature on a SymPy-derived single-integral
  reduction). Used elsewhere in the codebase but not in the
  current cross-method case set.
* ``"MMS"`` — manufactured-solution sources where the trial
  solution and source are derived symbolically. Used elsewhere; not
  in the current cross-method set.
* ``"ancillary"`` — RESERVED for L4-only references that don't
  back a pillar (e.g. another ORPHEUS solver). The protocol
  REJECTS ``ancillary`` truth values via a foundation gate
  (:func:`tests.cross_method.test_eigenvalue.test_case_pillar_is_not_ancillary`)
  to prevent reference contamination.

Truth traceability
------------------

Every :class:`CrossMethodCase` carries ``truth_source`` — the
primary literature citation. Failing tests print this verbatim so
the engineer who fixes the regression knows which *paper* the
agreement is breaking against.

Examples currently shipped:

* "Sood LA-13511 Table 4 (1999), citing Kaper-Lindeman-Leaf 1974
  NSE 54, 94 (Ref. 26)" — KLL Table I via Sood transcription.
  KLL itself is paywalled with no preprint; Sood's transcription is
  the primary access (see
  ``.claude/agent-memory/literature-researcher/kaper_lindeman_leaf_1974_fn_method.md``).
* "Grandjean-Siewert 1979 Table XI" — bare-slab c-sweep extension
  (c=1.10, 1.70, 1.90 not in Sood).
* "Neshat-Maiorino 1980 Table 1 Case 4 / Table 2 Burkart 1976
  'Exact'" — reflected-slab F_N benchmark; F_7 matches Burkart
  'Exact' to all published digits per NM Table 2.
* "Closed-form k_inf = νΣ_f/Σ_a; V_α1 algebraic identity" — V_α1
  closed-sphere identity (algebraic-of-record evidence).

Coverage matrix (current)
--------------------------

The :func:`~tests.cross_method.test_eigenvalue.test_coverage_matrix_diagnostic`
foundation gate prints the (case × adapter) matrix on every run.
As of 2026-05-03:

.. list-table::
   :header-rows: 1
   :widths: 32 12 12 16 16 12

   * - Case set
     - n_cases
     - fn_method
     - trajectory_resolvent
     - L1 truth backing
     - Agreement gate
   * - Bare-critical slab
     - 4
     - ✓
     - ✓
     - Sood / KLL
     - ✓ (slow)
   * - Bare-critical sphere
     - 3
     - ✓
     - ✓
     - Sood / KLL
     - ✓
   * - Reflected slab
     - 4
     - ✓
     - —
     - Sood + NM 1980
     - one-sided
   * - Closed sphere k_inf
     - 1
     - ✗ (uses ``compute_kinf_*`` direct)
     - ✓
     - V_α1 identity
     - one-sided
   * - GS Table XI parametric
     - 5
     - ✓
     - — (no registry case)
     - GS 1979
     - one-sided

Total: **17 cross-method cases**, **84 collected tests** (foundation
+ truth + agreement gates parametrised over the case sets).

Adding a new solver
-------------------

A typical workflow for adding e.g. ``spectral_resolvent`` to the
slab regression net:

1. Write a ``SpectralResolventSlabAdapter`` dataclass in
   :mod:`tests.cross_method.adapters`. Required attributes:
   ``name``, ``method``, ``geometry``. Required method:
   ``solve(case) -> ScalarResult``.

2. The adapter:

   * pulls XS via :func:`mixture_to_fn_arrays` from the case's
     ``registry_case.materials`` (or parses inline parameters from
     ``case.notes``);
   * picks internal numerical parameters (panel count, basis size,
     etc.) appropriate to the case's tolerance floor;
   * performs unit conversions (mfp → cm, half-thickness → full
     slab);
   * returns a :class:`ScalarResult` with the right ``tag``.

3. Register the adapter in
   :data:`~tests.cross_method.adapters.ADAPTERS_BY_NAME`.

4. Opt cases in by adding the adapter name + tolerance to each
   case's ``tolerances`` mapping. The
   :func:`tests.cross_method.test_eigenvalue.test_case_tolerance_adapters_exist`
   foundation gate enforces consistency.

5. Add a per-adapter ``test_*_matches_truth`` rule in
   :mod:`tests.cross_method.test_eigenvalue` if you want a per-
   case truth gate.

6. Optionally add pairwise agreement gates against existing
   adapters. Use
   :func:`~tests.cross_method.protocol.agreement_tolerance` for
   the tolerance.

Adding a new case
-----------------

1. Identify the registry case (``La13511Case`` from
   :mod:`orpheus.derivations.continuous.sood_registry`, or
   inline parameters for cases without a registry entry).

2. Construct a :class:`CrossMethodCase` with:

   * ``case_id`` unique across :data:`ALL_CASES`;
   * ``truth_value`` traced to a primary citation in
     ``truth_source``;
   * ``pillar`` ∈ {closed-form, MMS, semi-analytical} —
     ``ancillary`` is rejected by the foundation gate;
   * ``tolerances`` opt-in mapping; every named adapter must
     exist in :data:`ADAPTERS_BY_NAME`.

3. Append to the appropriate case-set list (or define a new list
   if a new geometry / topology lands).

Architectural seam — relationship to wave3 meta-registry
---------------------------------------------------------

The wave3 plan (``.claude/plans/wave3/architecture.md``) sketches a
per-paper registry + meta-registry under
``orpheus/derivations/registry/``. The cross-method test protocol
**consumes** the existing :class:`La13511Case` schema; when wave3
lands and lifts the schema to the paper-agnostic ``PaperCase``, the
``registry_case`` field of :class:`CrossMethodCase` is the migration
seam.

Until wave3 implementation begins, the cross-method protocol stays
in :mod:`tests.cross_method` and uses the existing
:class:`La13511Case` as the case pointer.

Multi-group cross-method coverage gap (acknowledged)
----------------------------------------------------

Bare-critical slab/sphere is **inherently 1G** — neither fn_method
nor trajectory_resolvent natively ships the multi-group critical-
dimension solve. Per :doc:`/skills/vv-principles` §"1-group
degeneracy" 1G eigenvalue tests are degenerate (k = νΣ_f/Σ_a
shape-independent). The cross-method protocol acknowledges this gap
honestly:

* **closed-sphere α=1 cases** (``CLOSED_SPHERE_KINF_CASES``) extend
  to 2G+ trivially via the V_α1 identity (k_eff = k_inf at any
  group structure when α=1). Currently the case set ships only the
  1G fixture; adding 2G/4G fixtures is the natural next extension.
* **k_inf cases via fn_method's ``compute_kinf_*``** (1G/2G/mG
  closed forms) are tested in the per-method file
  :mod:`tests.derivations.test_fn_la13511_kinf` against
  ``kinf_homogeneous`` (the structurally-independent companion
  identity). Those gates are L1; not duplicated in the cross-
  method protocol.

The current cross-method ≥2G coverage is therefore one case
(closed-sphere 1G that algebraically generalises). Lifting this
gap requires:

* (longer-term) extending fn_method to multi-group critical-
  dimension via Sood Eq 76 generalisation, OR
* extending the cross-method protocol to consume
  ``trajectory_resolvent.solve_greens_function_sphere_mg`` against
  the closed-sphere α=1 / kinf identity at 2G+ — which only needs
  a new ``CrossMethodCase`` row + parameter-set extension to the
  closed-sphere adapter.

References
----------

* Sood, Forster, Parsons (1999), *Analytical Benchmark Test Set
  for Criticality Code Verification*, LA-13511.
* Kaper, Lindeman, Leaf (1974), *Nucl. Sci. Eng.* **54**, 94.
* Grandjean, Siewert (1979), *Nucl. Sci. Eng.* **69**, 161 — Table
  XI critical thicknesses.
* Neshat, Maiorino (1980), *Ann. Nucl. Energy* **7**, 79 — Table 1
  reflected-slab cases + Table 2 Burkart 1976 'Exact' comparison.
* Sanchez (1986), *J. Quant. Spec. Rad. Transfer* **35**, 121 —
  Variant α specular leg.
* Siewert, Thomas (1986), *Nucl. Sci. Eng.* **94**, 264 — sphere
  F_N method.
* :doc:`/skills/vv-principles` — claim taxonomy + pillar
  discipline + reference-contamination anti-pattern.
* :doc:`/skills/algebra-of-record` — structural-independence
  ladder.
* ``.claude/scratch/cross_method_test_protocol_assessment.md`` —
  the Phase-1 architectural assessment that produced this
  protocol.
