.. _theory-fn-method:

==========================================================================
F_N method — analytical benchmark family (Sood/Forster/Parsons LA-13511)
==========================================================================

.. note::

   **Scope of this page (first slice).** This page is a stub-grade
   theory document covering the **k_inf cases** of the LA-13511
   benchmark suite — Sood Eqs 18-32 and 72-76, the closed-form
   infinite-medium fission-eigenvalue formulae for 1G, 2G, and
   general multi-group reactors. The ``F_N method`` proper (slab,
   sphere, cylinder critical-radius via Wiener-Hopf factorisation +
   Case eigenfunctions) lives in *cited* journal papers
   (Kaper-Lindeman-Leaf 1974, Westfall-Metcalf 1973, Siewert-Thomas
   1986, Sanchez-Ganapol 1983) — those references are being acquired
   in parallel with this first slice. As they arrive the
   :mod:`orpheus.derivations.continuous.fn_method.{slab,sphere,
   cylinder}` subpackages will be populated and this page will be
   expanded into a full F_N narrative by the archivist agent.

   Until then, this page documents the algebra-of-record for the
   *easy half* of LA-13511 — the cases that need only the report's
   own Appendix A and serve as the V&V ground floor for the rest of
   the F_N project.

Why F_N at all?
================

The Variant α Green's-function family
(:doc:`peierls_greens`) is ORPHEUS's primary continuous-mu reference
for the angle-resolved transport eigenvalue and flux-shape problems
in compact and 2-surface geometries. It already cross-checks against
external benchmarks: Sood/Forster/Parsons ``Ua-1-0-CY`` cylinder
critical radius at 8.5e-6 (see
:func:`tests.derivations.test_peierls_greens_function_xverif_sood2003_cylinder`).

That cross-check, however, leans on a single published value. The
F_N method gives a **second, structurally-independent reference
solver** — sharing only ``numpy`` (above the trusted-library line per
the project's algebra-of-record discipline) — that produces the same critical radii,
:math:`k_\infty`, and flux ratios via a completely different
mathematical route (singular eigenfunction expansions, Wiener-Hopf
factorisation, F_N rank-N approximation). When both Variant α and
F_N agree at 5+ digits across the LA-13511 catalogue, the
verification chain is structurally far stronger than any single
reference can provide.

This is the same role the PS-1982 Nyström reference plays for
Variant α at vacuum BC: an independent published-method anchor that
exposes any structural bug Variant α might be hiding behind its own
algebra (the canonical worked-example pattern of the
algebra-of-record discipline).

The five-case complexity ramp
==============================

The literature memo
``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``
lays out a 5-case ramp covering structural diversity at monotone
difficulty. The first slice ships **Cases 1 + 5** fully; Cases 2-4
are stubs awaiting the F_N method specifications.

.. list-table:: First-slice case map (LA-13511)
   :header-rows: 1
   :widths: 12 24 16 14 14 20

   * - Case
     - Sood ID
     - Geometry
     - Groups
     - Reference
     - This-slice status
   * - 1
     - PUa-1-0-IN
     - infinite
     - 1
     - k_inf = 2.612903
     - **Branch-1 + Branch-2**
   * - 2
     - Ua-1-0-SL
     - slab
     - 1
     - r_c = 0.93772556 mfp
     - stub (XS only)
   * - 3
     - Ua-1-0-CY
     - cylinder
     - 1
     - r_c = 1.72500292 mfp
     - stub (XS only) — already
       cross-checked by Variant α
   * - 4
     - Ua-1-0-SP
     - sphere
     - 1
     - r_c = 2.4248249802 mfp
     - stub (XS only)
   * - 5
     - PU-2-0-IN
     - infinite
     - 2
     - k_inf = 2.683767;
       phi_2/phi_1 = 0.675229
     - **Branch-1 + Branch-2**

ORPHEUS vs Sood group convention
=================================

Sood numbers groups :math:`g=N` (fast) → :math:`g=1` (slow), the
reverse of typical nuclear-engineering convention. ORPHEUS uses the
standard convention :math:`g=0` (fast) → :math:`g=N-1` (slow). The
:mod:`orpheus.derivations.continuous.fn_method.benchmarks.la13511`
catalogue does the conversion at load time so consumers see ORPHEUS
ordering directly. The Branch-1 SymPy module
(:mod:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations`)
uses Sood's symbols verbatim so equations match the report
letter-for-letter; the conversion is purely a relabeling and the
algebra is identical for either side.

Branch-1 / Branch-2 algebra-of-record
======================================

Per the project's algebra-of-record discipline:

* :mod:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations`
  is the **canonical algebra-of-record**: all closed forms are
  derived symbolically from Sood's transport-equation reductions
  (Eqs 18, 21-25, 72-73), and each ``derive_*()`` function returns
  PASS flags when the symbolic identity closes to zero.
* :mod:`orpheus.derivations.continuous.fn_method.multi_group.k_inf`
  is the **production code**: numpy implementations of the same
  closed forms, structurally independent of
  :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`
  (which solves the dominant eigenvalue problem of :math:`A^{-1}F`
  via :func:`numpy.linalg.eig`).
* :func:`tests.derivations.test_fn_la13511_kinf` pins both branches
  with foundation-tagged tests, plus the Branch-1 ↔ Branch-2 +
  ORPHEUS-cross-implementation agreement gates.

V_fn — the seven verifications
===============================

.. _fn-method-V-fn1-1:

V_fn1.1 — 1G k_inf from balance equation
-----------------------------------------

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_19`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn1_1_kinf_1g_eq_19`.

   Brief: starting from Sood Eq 18,
   :math:`\Sigma_t \phi = \Sigma_s \phi + (\nu\Sigma_f / k_\infty)\phi`,
   solve symbolically for :math:`k_\infty` and verify it equals
   Sood Eq 19, :math:`k_\infty = \nu\Sigma_f / (\Sigma_t - \Sigma_s)`.
   The flux :math:`\phi` cancels — the 1G eigenvalue is flux-shape
   independent (the canonical 1G degeneracy in the
   numerical-bug-signatures catalog).

.. _fn-method-V-fn1-2:

V_fn1.2 — Eq 20 simplifies to Eq 19 (c factor cancels)
-------------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_20_simplifies_to_eq_19`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn1_2_kinf_eq_20_simplifies_to_eq_19`.

   Brief: Sood writes :math:`k_\infty` two ways — Eq 19 and Eq 20.
   Eq 20 includes the explicit "mean number of secondaries"
   :math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`. Substituting
   :math:`c` into Eq 20 leaves :math:`(c \cdot \Sigma_t)/(\Sigma_s +
   \nu\Sigma_f) = 1`, so the :math:`c` and :math:`\Sigma_t` factors
   cancel and Eq 20 collapses to Eq 19. SymPy ``simplify`` makes
   this mechanical; by-hand the cancellation requires four
   factor-cancellation steps.

.. _fn-method-V-fn2-1:

V_fn2.1 — 2G general k_inf from det(M) = 0
-------------------------------------------

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_2g_general_from_matrix`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn2_1_kinf_2g_general_from_matrix`.

   Brief: from Sood Eqs 21-22 (2G balance), rearrange (Eqs 23-24)
   into a 2x2 homogeneous linear system :math:`M(k)\,\vec\phi = 0`
   (Eq 25). Critical fission balance requires :math:`\det(M) = 0`,
   which is a quadratic in :math:`k`. One root is :math:`k=0`
   (trivial); the other is the general 2G :math:`k_\infty`.

   **Discovery: Sood Eq 28 has a typo.** As printed, the chi-numerator
   factors are mis-paired with their respective :math:`\Sigma_g^{\rm
   rem}` removal cross sections. The correct form (verified
   symbolically via :math:`\det(M)=0`) swaps the pairing so that the
   no-upscatter limit (:math:`\Sigma_{21s}=0`) reduces to the
   correctly-printed Eq 29. The numerical reference values in
   LA-13511 (k_inf = 2.683767, phi_ratio = 0.675229) agree with the
   *corrected* Eq 28 + the *printed* Eq 29.

.. _fn-method-V-fn2-2:

V_fn2.2 — Eq 29 makes det(M) = 0 at no-upscatter
-------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn2_2_kinf_2g_no_upscatter_makes_det_zero`.

   Brief: independent verification of Sood Eq 29 — substitute the
   printed Eq 29 closed form into :math:`\det(M)` with
   :math:`\Sigma_{21s}=0` and verify the determinant simplifies to
   zero. Independent of V_fn2.1 (which derives the formula from
   scratch); V_fn2.2 verifies the published formula directly.

.. _fn-method-V-fn2-3:

V_fn2.3 — phi_2/phi_1 from chi-sum + balance (Eq 32)
-----------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn2_3_phi_ratio_eq_32`.

   Brief: adding Eqs 23 + 24 with :math:`\chi_1 + \chi_2 = 1`
   eliminates the :math:`\chi_g` from the resulting relation
   (Sood Eq 30). Solving for :math:`\phi_2/\phi_1` at the no-upscatter
   limit gives Sood Eq 32. SymPy verifies both the chi-elimination
   (the :math:`\chi_1` coefficient must vanish) and the resulting
   ratio identity.

.. _fn-method-V-fn-mg-1:

V_fnMG.1 — Eq 76 for G=2 is the trace of a rank-1 matrix
---------------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn_mg_1_eq_76_g2_form`.

   Brief: the general MG balance (Sood Eq 72) reduces to a single
   matrix inversion (Eq 76), :math:`k_\infty = \overline{\nu\Sigma_f}\,
   (\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})^{-1}\,
   \bar\chi`. The matrix :math:`A^{-1}\,\chi\,(\nu\Sigma_f)^T` is
   rank-1 (outer product structure), so its dominant eigenvalue
   equals its trace, which equals the scalar Eq 76 expression. SymPy
   verifies this for G=2 algebraically; for G ≥ 5 the symbolic
   eigenvalue closes via :func:`sympy.Matrix.eigenvals` no longer
   produces a closed form (Abel-Ruffini), but the matrix-vector
   identity Eq 76 still holds for arbitrary G.

.. _fn-method-V-fn-mg-2:

V_fnMG.2 — Eq 76 with G=1 reduces to Eq 19
-------------------------------------------

.. note:: TODO — Archivist expansion needed.

   Test gate:
   :func:`tests.derivations.test_fn_la13511_kinf.test_v_fn_mg_2_reduces_to_1g`.

   Brief: trivial dimensional-reduction check — Eq 76 with all
   matrices and vectors at G=1 collapses to scalar arithmetic and
   produces Eq 19 exactly. The MG infrastructure must reproduce the
   1G result bit-for-bit; this is enforced via
   :func:`tests.derivations.test_fn_la13511_kinf.test_kinf_mg_reduces_to_kinf_1g_at_n_groups_1`
   on the Branch-2 numpy side as well.

Cross-check claim
==================

The first slice's load-bearing cross-check claim:

   :func:`orpheus.derivations.continuous.fn_method.multi_group.compute_kinf_mg`
   and :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`
   agree on every LA-13511 first-slice case to ≥ 12 digits.

These two solvers are structurally independent above the
trusted-library line: F_N evaluates the closed form Sood Eq 76
directly via :func:`numpy.linalg.solve`, while ``kinf_homogeneous``
solves the dominant-eigenvalue problem of :math:`A^{-1}F` via
:func:`numpy.linalg.eig`. Disagreement would point at a real
implementation bug in one or the other. Foundation-test gate:
:func:`tests.derivations.test_fn_la13511_kinf.test_kinf_mg_agrees_with_existing_orpheus_kinf_homogeneous`.

When the F_N slab/sphere/cylinder solvers are added (after
Kaper-Lindeman-Leaf 1974 + Westfall-Metcalf 1973 are acquired), the
cross-check claim extends to:

   F_N reference solver and Variant α reference solver agree on
   the same physics to ≥ 5 digits for every overlapping LA-13511
   case (sphere primary, cylinder bonus, slab if Variant α slab is
   added).

This is the structurally-independent cross-check Plan 2 needs and
is the strategic motivation for the F_N project as a whole.

References
==========

* Sood, A., Forster, R.A., Parsons, D.K. (1999). *Analytical
  Benchmark Test Set for Criticality Code Verification.* Los
  Alamos National Laboratory report LA-13511. PDF at
  ``scratch/literature/Sood Foster Parsons (1999)Analytical
  Benchmark Test Set for Criticality Code Verification.pdf``.
* Sood, A., Forster, R.A., Parsons, D.K. (2003). *Analytical
  Benchmark Test Set for Criticality Code Verification.*
  *Progress in Nuclear Energy* 42, 55. (Journal-published
  condensation.)
* Kaper, H.G., Lindeman, A.J., Leaf, G.K. (1974). *Benchmark
  Values for the Slab and Sphere Criticality Problem in One-Group
  Neutron Transport Theory.* *Nuclear Science and Engineering* 54,
  94. — pending acquisition for Cases 2 + 4.
* Westfall, R.M., Metcalf, D.R. (1973). *Singular Eigenfunction
  Solution of the Monoenergetic Neutron Transport Equation for
  Finite Radially Reflected Critical Cylinders.* *Nuclear Science
  and Engineering* 52, 1. — pending acquisition for Case 3.
* Siewert, C.E., Thomas, J.R. (1986). *On Two-Group Critical
  Problems in Neutron Transport Theory.* *Nuclear Science and
  Engineering* 94, 264. — pending acquisition for 2G sphere
  bonus case.

Internal references:

* Literature memo:
  ``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``.
* Closeout memo:
  ``.claude/agent-memory/method-implementer/fn_method_kinf_first_slice.md``.
* :doc:`peierls_greens` — companion Variant α reference family.
