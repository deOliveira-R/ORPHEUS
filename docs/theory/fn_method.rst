.. _theory-fn-method:

==========================================================================
F_N method — analytical benchmark family (Sood/Forster/Parsons LA-13511)
==========================================================================

.. note::

   **Scope of this page.** Stub-grade theory document covering:

   * **k_inf cases** (Sood LA-13511 Eqs 18-32, 72-76) — first slice.
   * **Slab F_N** (Siewert-Benoist 1979 Part I + Grandjean-Siewert
     1979 Part II) — second slice; case ``Ua-1-0-SL``.
   * **Sphere F_N** (Siewert-Thomas 1986) — third slice; case
     ``Ua-1-0-SP``. The slab and sphere F_N share a unified assembler
     parametrised by ``geometry_sign ∈ {+1, -1}``.

   Cylinder F_N (Westfall-Metcalf 1972) is pending — the cylinder
   ``Ua-1-0-CY`` is currently cross-checked by Variant α directly
   against the Sood truth value at 8.5e-6.

   The archivist agent will expand each TODO marker into the full
   narrative.

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
     - **Branch-1 + Branch-2** (F_N method)
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
     - **Branch-1 + Branch-2** (true F_N
       method, Siewert-Thomas 1986)
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
:mod:`orpheus.derivations.continuous.sood_registry.la13511`
catalogue (formerly ``fn_method.benchmarks.la13511``; now method-
agnostic — see :doc:`sood_registry`) does the conversion at load
time so consumers see ORPHEUS ordering directly. The Branch-1 SymPy
module
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

Second slice — slab + sphere bare-critical solvers
====================================================

The second slice ships **Cases 2 + 4**: bare-critical slab
(``Ua-1-0-SL``) and bare-critical sphere (``Ua-1-0-SP``) at
:math:`c = 1.30`. The required reference papers
(Siewert-Benoist 1979 Part I, Grandjean-Siewert 1979 Part II,
Kaper-Lindeman-Leaf 1974, Pomraning-Siewert 1982) are all locally
available in ``scratch/literature/``.

Slab F_N method
---------------

The slab F_N solver
(:func:`orpheus.derivations.continuous.fn_method.slab.solve_fn_slab_bare_critical`)
implements the Siewert-Benoist Part I / Grandjean-Siewert Part II
F_N collocation system. For the symmetric critical slab,

.. math::

   \sum_{\alpha=0}^{N} a_\alpha
   \big[B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta}\,A_\alpha(\xi_\beta)\big]
   = 0, \qquad \beta = 0, \ldots, N

with closed-form moment recursions

.. math::

   B_\alpha(\xi) &= \xi B_{\alpha-1}(\xi) - 1/(\alpha+1), \\
   A_\alpha(\xi) &= -\xi A_{\alpha-1}(\xi) + 1/(\alpha+1),

seed by

.. math::

   B_0(\xi) &= 2/c - 1 - \xi\log(1 + 1/\xi), \\
   A_0(\xi) &= 1 - \xi\log(1 + 1/\xi),

and collocation points :math:`\xi_0 = \nu_0,\,\xi_1 = 0,\,\xi_2 = 1`
plus the remaining :math:`N - 2` interior points equally spaced in
:math:`(0, 1)`. For multiplying media (:math:`c > 1`),
:math:`\nu_0 = i u_0` is purely imaginary; the system is built in
complex linear algebra and the bare-critical thickness :math:`a` is
the real positive root of the complex determinant
:math:`\det M(a) = 0`.

Sphere F_N — Siewert-Thomas 1986 unification
---------------------------------------------

The **sphere F_N method** is provided by Siewert-Thomas 1986 (*Nucl.
Sci. Eng.* **94**, 264). The paper develops F_N for the two-group
critical sphere; the **1G specialisation** used here drops the matrix
structure to scalars (:math:`C \to c`, :math:`\Theta \to 1`,
:math:`\Lambda(z) \to 1 - cz\,\mathrm{atanh}(1/z)`) per the cleanly-
mechanical reduction documented in
:ref:`V_fn-sphere-fn.4 <fn-method-V-fn-sphere-fn-4>`.

The load-bearing structural insight (Siewert-Thomas p. 268, verbatim):

   "the critical sphere problem requires only that Eq. (4) be
   changed to read [Eq. 46], and that we interpret a as the critical
   radius, we have incorporated the relevant minus sign in our
   developed equations"

— means **slab and sphere F_N differ by exactly one sign on the
boundary attenuation term**. This is captured in the ORPHEUS
implementation by parametrising the unified F_N matrix assembler
(:func:`orpheus.derivations.continuous.fn_method.core.fn_matrix.assemble_fn_matrix`)
on a ``geometry_sign`` argument: ``+1`` for slab (Eq. 4 symmetric BC),
``-1`` for sphere (Eq. 46 anti-symmetric BC).

.. math::

   M_{\beta, \alpha}(R) = B_\alpha(\xi_\beta)
                          + s\,e^{-2R/\xi_\beta}\,A_\alpha(\xi_\beta),
   \qquad s = \pm 1.

All other machinery is shared between geometries:

* The **moment-integral recursion** for :math:`B_\alpha` and
  :math:`A_\alpha` (verified in :ref:`V_fn-slab.1
  <fn-method-V-fn-slab-1>` through :ref:`V_fn-slab.4
  <fn-method-V-fn-slab-4>`).
* The **dispersion-root** :math:`\nu_0 = i u_0` for :math:`c > 1` (a
  medium property; same root for slab and sphere).
* The **Wiener-Hopf X-function**, which depends only on
  :math:`\Lambda(z)` — no geometry parameter enters its definition.
  Verified in :ref:`V_fn-sphere-fn.5 <fn-method-V-fn-sphere-fn-5>`.

The **only** sphere-specific code beyond the geometry sign is the
**collocation grid**: per Siewert-Thomas Eq. 38a, the sphere uses
shifted Chebyshev-of-the-first-kind nodes strictly inside :math:`(0,
1)`,

.. math::

   \xi_\beta = \frac{1}{2}\left[1 + \cos\!\left(\frac{\beta\pi}{N+1}
   \right)\right], \qquad \beta = 1, \ldots, N,

with :math:`\xi_0 = \nu_0 = i u_0`. This differs from the slab
Grandjean-Siewert grid (which includes :math:`\xi = 0` and
:math:`\xi = 1`). The reason: at :math:`\xi = 0` the geometry-sign-
bearing attenuation term vanishes (:math:`e^{-\infty} = 0`),
collapsing the row to a constant independent of :math:`R` and of
geometry sign. This creates a structural rank deficiency that masks
the genuine sphere root condition.

The Branch-2 production solver
(:func:`orpheus.derivations.continuous.fn_method.sphere.solve_fn_sphere_bare_critical`)
locates :math:`R_c` by minimising :math:`|\det M(R)|` over a
prominence-filtered first-local-minimum bracket scan + Brent
refinement. At :math:`N = 10` the achieved accuracy on Sood
``Ua-1-0-SP`` (:math:`R_c = 2.4248249802` mfp at :math:`c = 1.30`) is
:math:`3.6 \times 10^{-8}` absolute — three decades better than the
:math:`10^{-5}` target.

V_fn-slab — five SymPy verifications
-------------------------------------

.. _fn-method-V-fn-slab-1:

V_fn-slab.1 — B_α moment recursion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_B_recursion`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_1_B_recursion`.

   Brief: SymPy verifies the load-bearing algebraic-long-division
   identity :math:`\mu^{\alpha+1}/(\xi-\mu) = \xi\mu^\alpha/(\xi-\mu)
   - \mu^\alpha`, which is the key step in deriving
   :math:`B_\alpha(\xi) = \xi B_{\alpha-1}(\xi) - 1/(\alpha+1)`.

.. _fn-method-V-fn-slab-2:

V_fn-slab.2 — A_α moment recursion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_A_recursion`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_2_A_recursion`.

   Brief: same long-division pattern as V_fn-slab.1 but at the negative
   collocation argument :math:`-\xi`, giving the sign-flipped recursion
   :math:`A_\alpha(\xi) = -\xi A_{\alpha-1}(\xi) + 1/(\alpha+1)`.

.. _fn-method-V-fn-slab-3:

V_fn-slab.3 — B_0 long-division identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_B0_seed`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_3_B0_seed`.

   Brief: split :math:`\mu/(\xi-\mu) = -1 + \xi/(\xi-\mu)` and verify
   the resulting integral is elementary. The published form
   :math:`B_0(\xi) = 2/c - 1 - \xi\log(1+1/\xi)` includes the
   dispersion-relation :math:`\delta`-mass for :math:`\xi \in (0,1)`;
   for :math:`\xi > 1` the integrand is regular and the closed form
   is :math:`-1 + \xi\log(\xi/(\xi-1))`. The two regimes meet
   continuously through the principal-value branch.

.. _fn-method-V-fn-slab-4:

V_fn-slab.4 — A_0 seed integral
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_A0_seed`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_4_A0_seed`.

   Brief: SymPy evaluates :math:`\int_0^1 \mu/(\xi+\mu)\,d\mu = 1
   - \xi\log(1+1/\xi)` directly (no principal value needed since
   the pole at :math:`\mu = -\xi` is outside the integration domain).
   This gives :math:`A_0(\xi) = 1 - \xi\log(1 + 1/\xi)`.

.. _fn-method-V-fn-slab-5:

V_fn-slab.5 — Critical-slab determinant structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_critical_determinant_structure`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_5_critical_determinant`.

   Brief: the F_N collocation system for the symmetric bare-critical
   slab is :math:`M(a)\,\vec a = 0` with :math:`M(a)_{\beta\alpha}
   = B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta} A_\alpha(\xi_\beta)`.
   Non-trivial solutions exist iff :math:`\det M(a) = 0`. SymPy
   verifies the matrix structure for :math:`N = 0` (1x1 entry) and
   :math:`N = 1` (2x2 polynomial determinant) cases; the full closed-
   form determinant is intractable for general :math:`N`, so the
   Branch-2 numpy implementation evaluates :math:`\det M(a)` numerically
   and bisects.

V_fn-sphere-fn — five SymPy verifications (Siewert-Thomas 1986)
----------------------------------------------------------------

.. _fn-method-V-fn-sphere-fn-1:

V_fn-sphere-fn.1 — Slab/sphere BC sign-flip parameterisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_bc_sign_flip`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_1_bc_sign_flip`.

   Brief: the slab BC :math:`\Psi(-a, \mu) = \Psi(a, -\mu)` (Eq. 4,
   symmetric reflection) versus the sphere BC :math:`\Psi(-a, \mu) =
   -\Psi(a, \mu)` (Eq. 46, anti-symmetric from the
   :math:`r\Phi(r) \to \Psi(x, \mu)` substitution) differ by a single
   geometry sign :math:`s \in \{+1, -1\}`. SymPy verifies the
   parametrised BC and the propagation of :math:`s` into the F_N
   matrix-entry attenuation block.

.. _fn-method-V-fn-sphere-fn-2:

V_fn-sphere-fn.2 — Sphere F_N matrix entry from geometry_sign = -1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_fn_matrix_entry`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_2_matrix_entry`.

   Brief: the unified F_N matrix entry :math:`M_{\beta,\alpha}(R; s)
   = B_\alpha(\xi_\beta) + s\,e^{-2R/\xi_\beta}\,A_\alpha(\xi_\beta)`
   reduces to the published sphere form
   :math:`M_{\beta,\alpha}^{\rm sphere} = B_\alpha - e^{-2R/\xi}A_\alpha`
   when :math:`s = -1` and to the slab form
   :math:`M_{\beta,\alpha}^{\rm slab} = B_\alpha + e^{-2R/\xi}A_\alpha`
   when :math:`s = +1`. SymPy verifies the substitution and the
   distinctness of the two forms.

.. _fn-method-V-fn-sphere-fn-3:

V_fn-sphere-fn.3 — Sphere bare-critical = det M(R) = 0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_critical_condition`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_3_critical_condition`.

   Brief: the sphere F_N collocation system :math:`M(R)\,\vec a = 0`
   is homogeneous (no source for the bare-critical problem); a
   non-trivial solution exists iff :math:`\det M(R) = 0`. Same
   determinantal form as slab, only the matrix entries differ via
   ``geometry_sign``. SymPy verifies the 1×1 + 2×2 worked examples
   and the slab/sphere distinctness.

.. _fn-method-V-fn-sphere-fn-4:

V_fn-sphere-fn.4 — Siewert-Thomas 2G→1G reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_2g_to_1g_reduction`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_4_2g_to_1g_reduction`.

   Brief: Siewert-Thomas 1986 develops F_N for the 2G case. The 1G
   specialisation collapses every 2×2 matrix to a scalar:
   :math:`C \to c`, :math:`\Theta(\mu) \to 1`,
   :math:`\Lambda(z) \to 1 - cz\,\mathrm{atanh}(1/z)`, system size
   :math:`(2N+2)^2 \to (N+1)^2`. SymPy verifies the matrix-dimension
   collapse on a worked 2×2 → scalar example.

.. _fn-method-V-fn-sphere-fn-5:

V_fn-sphere-fn.5 — Wiener-Hopf X-function geometry-independence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_x_function_geometry_independence`.
   Test gate:
   :func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_5_x_function_geometry_independence`.

   Brief: the Wiener-Hopf X-function :math:`X(z) = (1-z)^{-1}\,
   \exp[(1/\pi)\int_0^1 \arg \Lambda^{+}(\tau)\,d\tau/(\tau-z)]`
   depends only on the dispersion function :math:`\Lambda(z) = 1 -
   cz\,\mathrm{atanh}(1/z)`, which is a medium-only quantity (depends
   on :math:`c` and :math:`z`, no geometry parameter). SymPy verifies
   that the integrand has no :math:`R` or :math:`a` symbols. This
   justifies reusing the slab X-function machinery for sphere
   verbatim.

Cross-check claim — second slice
==================================

The second slice's load-bearing cross-check claims:

* **Slab Ua-1-0-SL**: F_N slab solver and Variant α slab solver agree
  on :math:`k_{\rm eff} = 1` at the Sood truth :math:`a_c = 0.93772556`
  mfp to ≤ 5e-5. F_N at :math:`N = 10` reaches ~5e-6; Variant α slab
  at :math:`(n_x, n_\mu) = (48, 128)` reaches ~1e-5. Cross-check
  tolerance 5e-5 is the safe envelope. Foundation-test gate:
  :func:`tests.derivations.test_fn_la13511_slab_xverif.test_fn_slab_vs_variant_alpha_at_sood_ua_1_0_sl`.

* **Sphere Ua-1-0-SP**: F_N sphere (Siewert-Thomas 1986) returns
  :math:`R_c = 2.4248249802` mfp to ≤ 1e-5 (achieved 3.6e-8 at
  :math:`N = 10`); Variant α sphere at the F_N predicted radius gives
  :math:`k_{\rm eff} = 1` to ≤ 1e-5 (achieved 4.2e-6). Foundation-test
  gate:
  :func:`tests.derivations.test_fn_la13511_sphere_xverif.test_fn_sphere_vs_variant_alpha_sphere_at_sood_ua_1_0_sp`.

Together these establish the **structural-independence pillar** for
the Variant α slab + sphere prototypes against the published-method
references that produced the Sood/KLL truth values. The two methods
are mathematically disjoint: F_N method uses Case singular
eigenfunctions + Wiener-Hopf factorisation; Variant α uses bouncing-
trajectory operator with rank-2 closure. Their agreement at the same
physics is the strongest L1 evidence currently available in ORPHEUS
for either method.

The sphere upgrade from PS-1982 wrapper (previous slice) to true F_N
(Siewert-Thomas 1986) makes the cross-check **structurally stronger**:
PS-1982 and Variant α both reduce the Peierls integral equation by
different algebraic paths (procedurally independent only), whereas
F_N method works in the Case singular-eigenfunction representation
and never reduces to an integral equation in :math:`r` (genuinely
structurally independent above the trusted-library line).

Flux reconstruction — deferred
================================

The Sood Case 2 (``Ua-1-0-SL``) flux ratios at :math:`r/r_c \in
\{0.25, 0.5, 0.75, 1.0\}` (KLL Table III: 0.96695, 0.86863, 0.70552,
0.44619) are NOT verified in this slice. The discrete-mode cosine
approximation :math:`\phi(x) \approx \cos(x/u_0)` gives ~20 % error
at the slab edge — too coarse for a flux-shape verification. The
full F_N angular-flux reconstruction
:math:`\psi(\pm a, \mp\mu) = \sum_\alpha a_\alpha \mu^\alpha`
followed by self-consistent scalar-flux integration is documented
as a stub
(:func:`orpheus.derivations.continuous.fn_method.slab.fn_slab_flux_at_x_cosine_only`)
and deferred to a follow-up slice. The KLL Wiener-Hopf iterated
Fredholm representation
:math:`\phi(x) = 2a\,[\cos(x/u_0) + \int_0^1 A(\nu)\,e^{-b/\nu}\cosh(x/\nu)\,d\nu]`
is the alternative path to high-precision flux profiles; future work.

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
  94. **Used for Cases 2 + 4** (slab + sphere bare-critical).
  PDF locally available at ``scratch/literature/Kaper-Lindeman-
  Leaf(1974)Benchmark Values for the Slab and Sphere Criticality
  Problem in One-Group Neutron Transport Th.pdf``.
* Siewert, C.E., Benoist, P. (1979). *The F_N Method in Neutron-
  Transport Theory. Part I: Theory and Applications.* *Nuclear
  Science and Engineering* **69**, 156. **Slab F_N method
  specification**, used for Case 2.
* Grandjean, P., Siewert, C.E. (1979). *The F_N Method in Neutron-
  Transport Theory. Part II: Applications and Numerical Results.*
  *Nuclear Science and Engineering* **69**, 161. **Slab F_N
  numerical-results tables**, including the bare-critical thickness
  Table XI used as the second cross-check anchor for Case 2.
* Pomraning, G.C., Siewert, C.E. (1982). *On the integral form of
  the equation of transfer for a homogeneous sphere.* *J. Quant.
  Spec. Rad. Transfer* **28**, 503. **Spherical Peierls integral
  equation** — alternative kernel derivation (used in the previous
  PS-1982-wrapper sphere F_N stub, now superseded by the true F_N
  method per Siewert-Thomas 1986).
* Westfall, R.M., Metcalf, D.R. (1973). *Singular Eigenfunction
  Solution of the Monoenergetic Neutron Transport Equation for
  Finite Radially Reflected Critical Cylinders.* *Nuclear Science
  and Engineering* 52, 1. — pending acquisition for Case 3.
* Siewert, C.E., Thomas, J.R. (1986). *On Two-Group Critical
  Problems in Neutron Transport Theory.* *Nuclear Science and
  Engineering* **94**, 264-270. **Sphere F_N method specification**
  used for Case 4. The 1G specialisation is a clean reduction of the
  2G machinery (matrix-dimension collapse to scalars). PDF locally
  available at ``scratch/literature/Siewert-Thomas(1986)On Two-Group
  Critical Problems in Neutron-Transport Theory.pdf``.

Internal references:

* Literature memo:
  ``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``.
* Closeout memo:
  ``.claude/agent-memory/method-implementer/fn_method_kinf_first_slice.md``.
* :doc:`peierls_greens` — companion Variant α reference family.
