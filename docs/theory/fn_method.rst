.. _theory-fn-method:

==========================================================================
F_N method — analytical benchmark family (Sood/Forster/Parsons LA-13511)
==========================================================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this before modifying the F_N reference solver.**

- **What this is**: Pillar-2 reference family — Siewert-Benoist 1979 /
  Grandjean-Siewert 1979 / Siewert-Thomas 1986 — projecting the Case
  singular-eigenfunction expansion onto a finite Legendre basis via
  interior collocation. Pillar classification: **semi-analytical**
  (mathematics is symbolic; final root-find on
  :math:`\det M(R) = 0` is numerical).
- **Position in the V&V stack**: structurally-independent cross-check
  for :ref:`theory-singular-eigenfunction` and
  :ref:`theory-trajectory-resolvent` on the Sood LA-13511 truth set.
  The F_N method works in the Case ν-spectrum representation and never
  reduces to an integral equation in :math:`r` — genuinely
  structurally-independent of Variant α (which works on bouncing
  characteristics) above the trusted-library line. The
  boundary-collocated sister method (B_N) is reserved at
  :ref:`theory-bn-method`.
- **Coverage**: slab + sphere via unified assembler (parametrised by
  ``geometry_sign ∈ {+1, -1}``); bare-critical k_inf cases (LA-13511
  Eqs 18-32, 72-76); reflected-slab F_N (Neshat-Maiorino 1980); KLL
  1974 interior flux reconstruction; Atkinson product-Nyström for
  the Path A.i log-singular kernel. Cylinder ships via a different
  pillar in :ref:`theory-singular-eigenfunction` (Westfall-Metcalf
  1973).
- **Critical numerical fact**: the Path A.i flux reconstruction was
  hardened from 1–7 % accuracy to ~1.1e-5 (5000× improvement) via
  Atkinson 1972/1997 § 4.2 product-Simpson rule (ERR-036, fixed
  2026-05-03). Plain Gauss-Legendre on the (z, μ) tensor product
  silently truncates the diagonal of the log-singular kernel, hitting
  the textbook fingerprint :math:`\mathrm{err} \cdot n / \log n \approx
  \mathrm{const}`. Atkinson product-Simpson integrates the singular
  log piece **analytically** against the Lagrange basis. See
  :ref:`fn-method-atkinson-product-nystrom`.
- **Cross-references**: :ref:`theory-singular-eigenfunction` is the
  Atalay-anchored linearly-anisotropic family; :ref:`theory-sood-registry`
  is the truth-value catalogue; :ref:`theory-trajectory-resolvent`
  realises the same physics in :math:`(r, \mu)` phase space.


Capabilities at a glance — what references this module ships
=============================================================

**For readers trying to answer "what F_N method references do we have
for problem class X?"** the table below is the canonical index. The
table body is auto-generated at Sphinx build time by
:mod:`tools.verification.generate_capability_matrices` from the
registry function
:func:`orpheus.derivations.continuous.fn_method.cases.capability_rows`
— so the matrix cannot drift from the shipping registry.

.. include:: _fn_method_capability_matrix.inc.rst


Why F_N at all?
================

The Variant α Green's-function family (:doc:`trajectory_resolvent`) is
ORPHEUS's primary continuous-:math:`\mu` reference for the
angle-resolved transport eigenvalue and flux-shape problems in compact
and 2-surface geometries. It already cross-checks against external
benchmarks: Sood/Forster/Parsons ``Ua-1-0-CY`` cylinder critical
radius at 8.5e-6 (see
:func:`tests.derivations.test_trajectory_resolvent_xverif_sood2003_cylinder`).

That cross-check, however, leans on a single published value. The
F_N method gives a **second, structurally-independent reference
solver** — sharing only ``numpy`` and ``scipy`` (above the
trusted-library line per the project's algebra-of-record
discipline) — that produces
the same critical radii, :math:`k_\infty`, and flux ratios via a
completely different mathematical route (singular eigenfunction
expansions, Wiener-Hopf factorisation, F_N rank-N approximation).

The agreement is structurally meaningful because the two methods are
genuinely disjoint:

* **F_N method** works in the Case ν-spectrum representation. The
  half-space exit-distribution equation is imposed exactly at
  :math:`(N+1)` interior collocation points; the resulting
  homogeneous determinantal equation :math:`\det M(R) = 0` gives the
  critical dimension. There is no integral equation in :math:`r`;
  the radial profile is *output*, not *input*.
* **Variant α** integrates the Boltzmann equation along bouncing
  characteristics in :math:`(r, \mu, \phi)` phase space, with
  analytical bounce-period summation closing the operator at the
  surface. The discrete eigenvalue problem is in the surface
  amplitude space; the radial profile is reconstructed *afterwards*
  by quadrature along characteristics.

When both methods agree at 5+ digits across the LA-13511 catalogue,
the verification chain is structurally far stronger than any single
reference can provide. This is the same role the PS-1982 Nyström
reference plays for Variant α at vacuum BC: an independent
published-method anchor that exposes any structural bug Variant α
might be hiding behind its own algebra (see
:doc:`peierls_nystrom`).

The five-case complexity ramp
==============================

The literature memo
``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``
lays out a 5-case ramp covering structural diversity at monotone
difficulty. The first slice ships **Cases 1 + 5** fully; Cases 2, 4
ship as full F_N solvers; Case 3 (cylinder) is intentionally outside
the F_N pillar — Westfall–Metcalf 1973 explicitly notes the
Mitsis-style Wiener-Hopf is **non-convergent for the bare cylinder**
and ships under :ref:`theory-singular-eigenfunction` instead.

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
     - **Out of pillar** — see
       :ref:`theory-singular-eigenfunction`
       (WM-72 Mitsis-WM Fredholm)
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

Per the project's algebra-of-record discipline (see the
``algebra-of-record`` skill, preloaded for the
``method-implementer``, ``test-architect``, and ``archivist`` agents):

* :mod:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations`,
  :mod:`...origins.fn_slab_derivations`,
  :mod:`...origins.fn_sphere_derivations`,
  :mod:`...origins.fn_slab_reflected_derivations`,
  :mod:`...origins.fn_flux_reconstruction_derivations`, and
  :mod:`...origins.fn_projection_flux_derivations` are the
  **canonical algebra-of-record**: every closed form is derived
  symbolically from primary-source equations (Sood Eqs 18-32 / 72-76,
  Siewert-Benoist Eq 4-6, Grandjean-Siewert Eqs 9-12, Siewert-Thomas
  1986 Eqs 38a-46, KLL Eqs 7+15, Neshat-Maiorino Eqs 10-17), and
  each ``derive_*()`` function returns a PASS flag when the symbolic
  identity closes to zero.
* :mod:`orpheus.derivations.continuous.fn_method.multi_group.k_inf`,
  :mod:`...slab.one_group`, :mod:`...sphere.one_group`,
  :mod:`...slab.reflected`, :mod:`...slab.flux_reconstruction`,
  :mod:`...sphere.flux_reconstruction`, and
  :mod:`...peierls_atkinson_nystrom` are the **Branch-2 production
  code**: numpy/scipy implementations of the same closed forms,
  structurally independent of any ORPHEUS in-house primitive above
  the trusted-library line.
* :func:`tests.derivations.test_fn_la13511_kinf`,
  :func:`tests.derivations.test_fn_la13511_slab`,
  :func:`tests.derivations.test_fn_la13511_sphere`,
  :func:`tests.derivations.test_fn_la13511_slab_reflected`,
  :func:`tests.derivations.test_fn_la13511_slab_flux`,
  :func:`tests.derivations.test_fn_la13511_sphere_flux`,
  :func:`tests.derivations.test_fn_la13511_slab_xverif`,
  :func:`tests.derivations.test_fn_la13511_sphere_xverif`,
  :func:`tests.derivations.test_atkinson_product_nystrom`, and
  :func:`tests.derivations.test_path_ai_legacy_plain_gl_signature`
  pin both branches with foundation-tagged tests, plus the Branch-1 ↔
  Branch-2 + ORPHEUS-cross-implementation agreement gates.

V_fn — the seven verifications
===============================

The seven foundation-tagged SymPy verifications below cover the
infinite-medium :math:`k_\infty` arithmetic. They are individually
trivial in SymPy but tedious by hand; they pin the Sood algebra
letter-for-letter and they immediately catch any future XS
relabelling drift between the registry and the consumers.

.. _fn-method-V-fn1-1:

V_fn1.1 — 1G k_inf from balance equation
-----------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_19`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn1_1_kinf_1g_eq_19`.

Starting from Sood Eq 18 (the 1G integrated transport equation for
an infinite, isotropically-scattering, homogeneous medium),

.. math::
   :label: sood-eq18-1g-balance

   \Sigma_t\,\phi = \Sigma_s\,\phi + \frac{\nu\Sigma_f}{k_\infty}\,\phi

.. (vv-status rationale) derivation: Sood Eq 18 transcription — the 1G balance equation; reduces algebraically to Eq 19 (verified by V_fn1.1 in the SymPy origins module).
.. vv-status: sood-eq18-1g-balance documented

we factor :math:`\phi` and solve for :math:`k_\infty`:

.. math::
   :label: sood-eq19-kinf-1g

   k_\infty = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s} .

.. (vv-status rationale) derivation: Sood Eq 19 — the closed-form 1G k_inf result from Eq 18; verified by V_fn1.1 (test_v_fn1_1_kinf_1g_eq_19) at the SymPy level. The Branch-2 numpy implementation also reproduces Eq 19 bit-for-bit at G=1 reduction (test_kinf_mg_reduces_to_kinf_1g_at_n_groups_1).
.. vv-status: sood-eq19-kinf-1g documented

**The flux** :math:`\phi` **cancels.** This is the canonical 1G
degeneracy: the eigenvalue is a material-property ratio, computable
without solving the transport equation. As a verification claim, 1G
:math:`k_\infty` agreement carries **zero structural information
about any spatial, angular, or scattering operator** (`vv-principles`
§ "1-group degeneracy — canonical statement"). The test exists as a
sanity-anchor for the registry; the Branch-2 multi-group machinery
must reproduce this case bit-for-bit at G=1.

.. _fn-method-V-fn1-2:

V_fn1.2 — Eq 20 simplifies to Eq 19 (c factor cancels)
-------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_1g_eq_20_simplifies_to_eq_19`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn1_2_kinf_eq_20_simplifies_to_eq_19`.

Sood states the same 1G result two ways. Eq 19 is the clean form;
Eq 20 includes the explicit "mean number of secondaries"
:math:`c = (\Sigma_s + \nu\Sigma_f)/\Sigma_t`,

.. math::
   :label: sood-eq20-kinf-1g-c-form

   k_\infty = \frac{c\,\nu\Sigma_f\,\Sigma_t}
                   {(\Sigma_t - \Sigma_s)(\Sigma_s + \nu\Sigma_f)} .

.. (vv-status rationale) derivation: Sood Eq 20 — the same 1G k_inf restated with the explicit "secondaries per collision" c-factor; reduces to Eq 19 (verified by V_fn1.2, test_v_fn1_2_kinf_eq_20_simplifies_to_eq_19).
.. vv-status: sood-eq20-kinf-1g-c-form documented

Substituting :math:`c` into :eq:`sood-eq20-kinf-1g-c-form`:

.. math::

   k_\infty = \frac{\Sigma_s + \nu\Sigma_f}{\Sigma_t}
              \cdot \frac{\nu\Sigma_f\,\Sigma_t}
                         {(\Sigma_t - \Sigma_s)(\Sigma_s + \nu\Sigma_f)}
            = \frac{\nu\Sigma_f}{\Sigma_t - \Sigma_s} ,

so the :math:`c` and :math:`\Sigma_t` factors cancel and Eq 20
collapses to Eq 19. SymPy ``simplify`` makes this mechanical. The
identity is the trivial Branch-1 anchor for Case 1 (PUa-1-0-IN); it
also showcases why SymPy is the right tool — a by-hand reader needs
to mentally cancel four factors, and any single drop creates a
plausible-looking but wrong reduction.

.. _fn-method-V-fn2-1:

V_fn2.1 — 2G general k_inf from det(M) = 0
-------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_2g_general_from_matrix`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn2_1_kinf_2g_general_from_matrix`.

From Sood Eqs 21-22 (2G balance), Eqs 23-24 rearrange to a 2x2
homogeneous linear system :math:`M(k_\infty)\,\vec\phi = 0`
(Sood Eq 25). Critical fission balance requires :math:`\det M = 0`,
which is a quadratic in :math:`k`. One root is :math:`k=0` (trivial);
the other is the desired :math:`k_\infty`.

In Sood notation (:math:`g=2` fast, :math:`g=1` slow,
:math:`\Sigma_g^{\rm rem} = \Sigma_g - \Sigma_{ggs}`):

.. math::
   :label: sood-eq25-2g-matrix

   M = \begin{pmatrix}
         -(\Sigma_{21s} + \tfrac{\chi_2}{k}\nu_1\Sigma_{1f}) &
         \Sigma_2^{\rm rem} - \tfrac{\chi_2}{k}\nu_2\Sigma_{2f} \\[4pt]
         \Sigma_1^{\rm rem} - \tfrac{\chi_1}{k}\nu_1\Sigma_{1f} &
         -(\Sigma_{12s} + \tfrac{\chi_1}{k}\nu_2\Sigma_{2f})
       \end{pmatrix} .

.. (vv-status rationale) derivation: Sood Eq 25 transcription — the 2G fission-balance matrix that yields k_inf via det(M)=0; verified by V_fn2.1 SymPy derivation (the "derived general 2G formula"). The published Eq 28 has a typo; corrected form lives in the SymPy module.
.. vv-status: sood-eq25-2g-matrix documented

The SymPy derivation expands :math:`\det M(k) = 0`, multiplies through
by :math:`k^2` to land on a polynomial in :math:`k`, solves for the
two roots, and discards the :math:`k=0` solution. The surviving root
is the **derived general 2G formula** — the corrected Sood Eq 28.

**Discovery — Sood Eq 28 has a typo.** As printed, the :math:`\chi_g`
numerator factors are mis-paired with their respective
:math:`\Sigma_g^{\rm rem}` removal cross sections. With Eq 28 as
printed and :math:`\Sigma_{21s} = 0`, the result is :math:`k_\infty
\approx 2.862` for Case 5 (PU-2-0-IN) — not the published 2.684. The
correct form (verified symbolically via :math:`\det M = 0`) swaps
the pairing, after which:

* The no-upscatter limit (:math:`\Sigma_{21s} = 0`) reduces to the
  correctly-printed Eq 29.
* The numerical reference values in LA-13511 (k_inf = 2.683767,
  :math:`\phi_2/\phi_1 = 0.675229`) match.

The corrected Eq 28 (output of the SymPy ``derive_*``) is what the
Branch-2 ``compute_kinf_2g_general`` evaluates. The
``sood_registry`` notes field on ``PU_2_0_IN`` flags the typo with a
pointer to this V_fn.

.. _fn-method-V-fn2-2:

V_fn2.2 — Eq 29 makes det(M) = 0 at no-upscatter
-------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_2g_no_upscatter`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn2_2_kinf_2g_no_upscatter_makes_det_zero`.

Independent verification of Sood Eq 29: substitute the printed
Eq 29 closed form

.. math::
   :label: sood-eq29-kinf-2g-no-upscatter

   k_\infty = \frac{\chi_1\,\nu_1\,\Sigma_{1f}}{\Sigma_1^{\rm rem}}
            + \chi_2\!\left[
                  \frac{\nu_1\,\Sigma_{1f}\,\Sigma_{12s}}
                       {\Sigma_1^{\rm rem}\,\Sigma_2^{\rm rem}}
                + \frac{\nu_2\,\Sigma_{2f}}{\Sigma_2^{\rm rem}}
              \right]

.. (vv-status rationale) derivation: Sood Eq 29 — printed correctly in LA-13511; verified by V_fn2.2 (test_v_fn2_2_kinf_2g_no_upscatter_makes_det_zero) which substitutes Eq 29 into det(M) and confirms it vanishes.
.. vv-status: sood-eq29-kinf-2g-no-upscatter documented

into :math:`\det M(k)` with :math:`\Sigma_{21s} = 0` and verify the
determinant simplifies to zero. This is structurally independent of
V_fn2.1 (which derives the formula from scratch); V_fn2.2 verifies
the published formula directly.

The two checks together close the proof that:

1. The general-form derivation in V_fn2.1 is correct (Eq 25 → quadratic
   → non-trivial root).
2. The published Eq 29 is the correct no-upscatter limit (V_fn2.2).
3. Therefore, the printed Eq 28 must be a typo (it disagrees with
   the correct general form whose specialisation is the printed
   Eq 29).

.. _fn-method-V-fn2-3:

V_fn2.3 — phi_2/phi_1 from chi-sum + balance (Eq 32)
-----------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_phi_ratio_2g_no_upscatter`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn2_3_phi_ratio_eq_32`.

Adding the two 2G balance equations Eqs 23 + 24 with
:math:`\chi_1 + \chi_2 = 1` eliminates the :math:`\chi_g` from the
resulting relation (Sood Eq 30). Solving for
:math:`\phi_2/\phi_1` at the no-upscatter limit gives Sood Eq 32:

.. math::
   :label: sood-eq32-phi-ratio

   \frac{\phi_2}{\phi_1}
   = \frac{\Sigma_{12s}}
          {\Sigma_2^{\rm rem} - \dfrac{\nu_2\Sigma_{2f}}{k_\infty}} .

.. (vv-status rationale) derivation: Sood Eq 32 — the 2G no-upscatter flux ratio; verified by V_fn2.3 (test_v_fn2_3_phi_ratio_eq_32). Independent of fission spectrum splitting (chi cancels under chi_1+chi_2=1).
.. vv-status: sood-eq32-phi-ratio documented

SymPy verifies both the chi-elimination (the :math:`\chi_1`
coefficient must vanish identically) and the resulting ratio
identity. The flux-ratio is **independent of the fission spectrum
splitting** :math:`(\chi_1, \chi_2)` — only the total fission rate
appears via :math:`k_\infty`. This is a substantive cross-check
because Sood Eq 32's specific form is not obviously chi-independent
from inspection; SymPy proves it.

.. _fn-method-V-fn-mg-1:

V_fnMG.1 — Eq 76 for G=2 is the trace of a rank-1 matrix
---------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_matrix_form`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn_mg_1_eq_76_g2_form`.

The general G-group balance (Sood Eq 72) reduces to a single
matrix-vector identity (Sood Eq 76):

.. math::
   :label: sood-eq76-kinf-mg

   k_\infty = \overline{\nu\Sigma_f}^{\!\top}\,
              \big(\overline{\overline{\Sigma_t}}
                  - \overline{\overline{\Sigma_s}}\big)^{-1}\,
              \bar\chi .

.. (vv-status rationale) derivation: Sood Eq 76 — the multi-group k_inf as a single matrix-vector identity; verified by V_fnMG.1 (test_v_fn_mg_1_eq_76_g2_form) at G=2 where SymPy rank-1 trace closes; verified numerically at higher G against kinf_homogeneous (np.linalg.eig).
.. vv-status: sood-eq76-kinf-mg documented

SymPy verifies for G=2 that the matrix
:math:`A^{-1}\,\bar\chi\,(\overline{\nu\Sigma_f})^{\!\top}` is
**rank-1** (outer-product structure — the :math:`\bar\chi` column
times the :math:`(\overline{\nu\Sigma_f})^{\!\top}` row). For a
rank-1 matrix the dominant eigenvalue equals the trace, which equals
the scalar Eq 76 expression.

For G ≥ 5 the symbolic eigenvalue closes via
:func:`sympy.Matrix.eigenvals` no longer produces a closed form
(Abel-Ruffini theorem — degree ≥ 5 polynomials have no closed-form
root in radicals). But the matrix-vector identity Eq 76 still holds
for arbitrary G, because the rank-1 structure is preserved. The
Branch-2 ``compute_kinf_mg`` therefore evaluates Eq 76 directly via
:func:`numpy.linalg.solve` — no eigenvalue solver is needed at all.
This is one of the cleanest examples in the project of the
"minimal-SymPy + scaling-argument" discipline (see
``algebra-of-record`` § "Minimal-SymPy + scaling argument"): SymPy
verifies the algebraic structure on the smallest non-trivial case
(G=2), the structure scales mechanically with G via array shape,
and the L1 cross-check at higher G is numerical against
:func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`
(which solves the dominant eigenvalue of :math:`A^{-1}F` via
:func:`numpy.linalg.eig`).

.. _fn-method-V-fn-mg-2:

V_fnMG.2 — Eq 76 with G=1 reduces to Eq 19
-------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.k_inf_derivations.derive_kinf_mg_reduces_to_1g`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_kinf.test_v_fn_mg_2_reduces_to_1g`.

Trivial dimensional-reduction check — Eq 76 with all matrices and
vectors at G=1 collapses to scalar arithmetic and produces Eq 19
exactly. The MG infrastructure must reproduce the 1G result
bit-for-bit; this is enforced via
:func:`tests.derivations.test_fn_la13511_kinf.test_kinf_mg_reduces_to_kinf_1g_at_n_groups_1`
on the Branch-2 numpy side as well. Together V_fnMG.2 + the
foundation Branch-2 reduction test pin the "G=1 is a special case
of G ≥ 2 solver" invariant on both branches.

Cross-check claim — first slice
================================

The first-slice load-bearing cross-check claim:

   :func:`orpheus.derivations.continuous.fn_method.multi_group.compute_kinf_mg`
   and :func:`orpheus.derivations.common.eigenvalue.kinf_homogeneous`
   agree on every LA-13511 first-slice case to ≥ 12 digits.

These two solvers are structurally independent above the trusted-
library line:

* **F_N path** evaluates the closed form Sood Eq 76 directly via
  :func:`numpy.linalg.solve` on
  :math:`(\overline{\overline{\Sigma_t}} - \overline{\overline{\Sigma_s}})\,
  \vec u = \bar\chi`, then dots with :math:`\overline{\nu\Sigma_f}`.
* **kinf_homogeneous path** solves the dominant-eigenvalue problem
  of the production-loss matrix
  :math:`A^{-1}F = (\overline{\overline{\Sigma_t}} -
  \overline{\overline{\Sigma_s}})^{-1}\,(\bar\chi
  \otimes \overline{\nu\Sigma_f})`
  via :func:`numpy.linalg.eig`.

Disagreement would point at a real implementation bug in one or the
other. Foundation-test gate:
:func:`tests.derivations.test_fn_la13511_kinf.test_kinf_mg_agrees_with_existing_orpheus_kinf_homogeneous`.

When the F_N slab/sphere/cylinder solvers are added, the cross-check
extends to:

   F_N reference solver and Variant α reference solver agree on the
   same physics to ≥ 5 digits for every overlapping LA-13511 case
   (sphere primary, cylinder via :ref:`theory-singular-eigenfunction`,
   slab via :func:`...slab.solve_fn_slab_bare_critical`).

This is the structurally-independent cross-check the verification
strategy needs and is the strategic motivation for the F_N project as
a whole.

Second slice — slab + sphere bare-critical solvers
====================================================

The second slice ships **Cases 2 + 4**: bare-critical slab
(``Ua-1-0-SL``) and bare-critical sphere (``Ua-1-0-SP``) at
:math:`c = 1.30`. The required reference papers (Siewert-Benoist
1979 Part I, Grandjean-Siewert 1979 Part II, Kaper-Lindeman-Leaf
1974, Pomraning-Siewert 1982, Siewert-Thomas 1986) are all locally
available in ``scratch/literature/``.

Slab F_N method
---------------

The slab F_N solver
(:func:`orpheus.derivations.continuous.fn_method.slab.solve_fn_slab_bare_critical`)
implements the Siewert-Benoist Part I / Grandjean-Siewert Part II
F_N collocation system. For the symmetric critical slab on
:math:`[-a, a]` with vacuum BC at both faces,

.. math::
   :label: fn-slab-collocation

   \sum_{\alpha=0}^{N} a_\alpha
   \big[B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta}\,A_\alpha(\xi_\beta)\big]
   = 0, \qquad \beta = 0, \ldots, N

.. (vv-status rationale) governing: Defines the Grandjean-Siewert F_N collocation system for the bare-critical slab — the entire slab F_N solver is the verification (test_v_fn_slab_5_critical_determinant + the L1 numerical gate test_fn_slab_solver_at_sood_ua_1_0_sl).
.. vv-status: fn-slab-collocation documented

with closed-form moment recursions derived in V_fn-slab.1 — V_fn-slab.4
(see below):

.. math::

   B_\alpha(\xi) &= \xi B_{\alpha-1}(\xi) - 1/(\alpha+1), \\
   A_\alpha(\xi) &= -\xi A_{\alpha-1}(\xi) + 1/(\alpha+1),

and seeds

.. math::

   B_0(\xi) &= 2/c - 1 - \xi\log(1 + 1/\xi), \\
   A_0(\xi) &= 1 - \xi\log(1 + 1/\xi) .

Collocation points: :math:`\xi_0 = \nu_0` (the discrete dispersion
root), :math:`\xi_1 = 0`, :math:`\xi_2 = 1`, plus the remaining
:math:`N - 2` interior points equally spaced in :math:`(0, 1)`
(Grandjean-Siewert 1979 § II.B). For multiplying media
(:math:`c > 1`), :math:`\nu_0 = i u_0` is purely imaginary; the
system is built in complex linear algebra and the bare-critical
half-thickness :math:`a` is the real positive root of the complex
determinant :math:`\det M(a) = 0`.

The Branch-2 production solver locates :math:`a_c` by minimising
:math:`|\det M(a)|` over a prominence-filtered first-local-minimum
bracket scan + Brent refinement. At :math:`N = 12` the achieved
accuracy on Sood ``Ua-1-0-SL`` (:math:`a_c = 0.93772556` mfp at
:math:`c = 1.30`) is :math:`\le 5 \times 10^{-6}` absolute.

Sphere F_N — Siewert-Thomas 1986 unification
---------------------------------------------

The **sphere F_N method** is provided by Siewert-Thomas 1986 (*Nucl.
Sci. Eng.* **94**, 264). The paper develops F_N for the two-group
critical sphere; the **1G specialisation** used here drops the
matrix structure to scalars (:math:`C \to c`,
:math:`\Theta \to 1`,
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
   :label: fn-unified-matrix-entry

   M_{\beta, \alpha}(R; s)
     = B_\alpha(\xi_\beta)
     + s\,e^{-2R/\xi_\beta}\,A_\alpha(\xi_\beta),
     \qquad s = \pm 1 .

.. (vv-status rationale) governing: The unified slab/sphere F_N matrix entry parametrised by geometry_sign s ∈ {+1, -1}; verified by V_fn-sphere-fn.2 (test_v_fn_sphere_fn_2_matrix_entry).
.. vv-status: fn-unified-matrix-entry documented

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
bearing attenuation term :math:`e^{-2R/\xi}` vanishes
(:math:`e^{-\infty} = 0`), collapsing the row to a constant
independent of :math:`R` and of geometry sign. This creates a
structural rank deficiency that masks the genuine sphere root
condition. Excluding :math:`\xi = 0` from the sphere grid eliminates
the rank deficiency without sacrificing accuracy at the relevant
dispersion-root collocation point :math:`\xi_0 = i u_0`.

The Branch-2 production solver
(:func:`orpheus.derivations.continuous.fn_method.sphere.solve_fn_sphere_bare_critical`)
locates :math:`R_c` by minimising :math:`|\det M(R)|` over a
prominence-filtered first-local-minimum bracket scan + Brent
refinement. At :math:`N = 10` the achieved accuracy on Sood
``Ua-1-0-SP`` (:math:`R_c = 2.4248249802` mfp at :math:`c = 1.30`)
is :math:`3.6 \times 10^{-8}` absolute — three decades better than
the :math:`10^{-5}` target.

V_fn-slab — five SymPy verifications
-------------------------------------

.. _fn-method-V-fn-slab-1:

V_fn-slab.1 — B_α moment recursion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_B_recursion`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_1_B_recursion`.

The :math:`B` moments are the integrals
:math:`B_\alpha(\xi) = \int_0^1 \mu^\alpha/(\xi - \mu)\, d\mu` minus
the dispersion residue (the Cauchy P.V. handling for
:math:`\xi \in (0, 1)`). The recursion
:math:`B_\alpha(\xi) = \xi B_{\alpha-1}(\xi) - 1/(\alpha+1)` is the
load-bearing identity for the F_N matrix assembly. SymPy verifies
the algebraic-long-division step

.. math::
   :label: fn-slab-B-long-division

   \frac{\mu^{\alpha+1}}{\xi - \mu}
   = \frac{\xi\mu^\alpha}{\xi - \mu} - \mu^\alpha

.. (vv-status rationale) derivation: Polynomial-long-division identity used in the B_α moment recursion; verified by V_fn-slab.1 (test_v_fn_slab_1_B_recursion).
.. vv-status: fn-slab-B-long-division documented

(elementary polynomial division), which when integrated against
:math:`d\mu` over :math:`(0, 1)` gives the recursion directly: the
first term yields :math:`\xi B_{\alpha-1}(\xi)`, the second yields
:math:`-1/(\alpha+1)`.

.. _fn-method-V-fn-slab-2:

V_fn-slab.2 — A_α moment recursion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_A_recursion`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_2_A_recursion`.

Same long-division pattern as V_fn-slab.1 but at the negative
collocation argument :math:`-\xi`. The substitution
:math:`\xi \to -\xi` in the long-division identity flips the sign of
the ratio term, giving

.. math::

   A_\alpha(\xi) = -\xi A_{\alpha-1}(\xi) + 1/(\alpha+1) .

The two-sided recursion (B from positive ξ-pole, A from negative
ξ-pole) is the F_N matrix's algebraic spine. It is also why the
moment-integral cost is :math:`O(N)` per :math:`\xi_\beta` rather
than :math:`O(N^2)` — each new α reuses the previous α's value.

.. _fn-method-V-fn-slab-3:

V_fn-slab.3 — B_0 long-division identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_B0_seed`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_3_B0_seed`.

Split :math:`\mu/(\xi-\mu) = -1 + \xi/(\xi-\mu)` and verify that the
resulting integral is elementary. The published seed

.. math::

   B_0(\xi) = \frac{2}{c} - 1 - \xi\log(1 + 1/\xi)

includes the dispersion-relation :math:`\delta`-mass (the
:math:`2/c` factor) for :math:`\xi \in (0, 1)`; for :math:`\xi > 1`
the integrand is regular and the closed form is :math:`-1 +
\xi\log(\xi/(\xi-1))`. The two regimes meet continuously through
the principal-value branch — SymPy verifies the algebraic agreement
at :math:`\xi = 1^+` symbolically (using the limit
:math:`\lim_{\xi\to 1^+} \xi\log(\xi/(\xi-1)) = +\infty` properly
handled by the principal-value definition).

.. _fn-method-V-fn-slab-4:

V_fn-slab.4 — A_0 seed integral
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_A0_seed`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_4_A0_seed`.

SymPy evaluates :math:`\int_0^1 \mu/(\xi+\mu)\,d\mu = 1 -
\xi\log(1+1/\xi)` directly. No principal value is needed since the
pole at :math:`\mu = -\xi` is outside the integration domain
:math:`(0, 1)` for any :math:`\xi > 0`. The clean evaluation gives
:math:`A_0(\xi) = 1 - \xi\log(1 + 1/\xi)`.

.. _fn-method-V-fn-slab-5:

V_fn-slab.5 — Critical-slab determinant structure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_derivations.derive_critical_determinant_structure`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab.test_v_fn_slab_5_critical_determinant`.

The F_N collocation system :eq:`fn-slab-collocation` is
:math:`M(a)\,\vec a = 0` with
:math:`M(a)_{\beta\alpha} = B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta}\,
A_\alpha(\xi_\beta)`. Non-trivial solutions exist iff
:math:`\det M(a) = 0`. SymPy verifies the matrix structure for
:math:`N = 0` (1×1 entry) and :math:`N = 1` (2×2 polynomial
determinant) cases; the full closed-form determinant is intractable
for general :math:`N`, so the Branch-2 numpy implementation
evaluates :math:`\det M(a)` numerically and bisects.

This is a deliberate Branch-1 termination: the SymPy module verifies
the algebraic structure (matrix entries are polynomials in the
moments + an exponential factor; the determinant is algebraically
well-defined), and the Branch-2 numpy implementation evaluates
numerically. The "minimal-SymPy + scaling-argument" pattern: SymPy
proves the structure on the smallest non-trivial dimension; the
structure scales mechanically with :math:`N`; the L1 cross-check at
operating :math:`N` is numerical against published Sood truth.

V_fn-sphere-fn — five SymPy verifications (Siewert-Thomas 1986)
----------------------------------------------------------------

.. _fn-method-V-fn-sphere-fn-1:

V_fn-sphere-fn.1 — Slab/sphere BC sign-flip parameterisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_bc_sign_flip`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_1_bc_sign_flip`.

The slab BC (Siewert-Benoist Eq. 4) is :math:`\Psi(-a, \mu) = \Psi(a,
-\mu)` (symmetric reflection in both space and angle). For the
sphere, the substitution :math:`r\Phi(r) \to \Psi(x, \mu)` (Mitsis
1963) maps the radial 1D transport equation to a slab-like form on
:math:`[-a, a]` with the sphere BC :math:`\Psi(-a, \mu) = -\Psi(a,
\mu)` (Siewert-Thomas Eq. 46, anti-symmetric — the
:math:`r \to -r` reflection picks up the :math:`r` factor's sign).

These two BCs differ by a single geometry sign :math:`s \in \{+1, -1\}`.
SymPy verifies the parametrised BC

.. math::

   \Psi(-a, \mu) = s\,\Psi(a, -s\mu),
   \qquad s = +1 \text{ slab},\; s = -1 \text{ sphere}

and the propagation of :math:`s` into the F_N matrix-entry
attenuation block — the :math:`A_\alpha` term picks up the sign
factor, the :math:`B_\alpha` term does not. This is the algebraic
spine of the unified slab/sphere implementation.

.. _fn-method-V-fn-sphere-fn-2:

V_fn-sphere-fn.2 — Sphere F_N matrix entry from geometry_sign = -1
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_fn_matrix_entry`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_2_matrix_entry`.

The unified F_N matrix entry :eq:`fn-unified-matrix-entry` reduces
to the published sphere form (Siewert-Thomas Eq. 46):

.. math::

   M_{\beta,\alpha}^{\rm sphere}(R)
   = B_\alpha(\xi_\beta) - e^{-2R/\xi_\beta}\,A_\alpha(\xi_\beta)

when :math:`s = -1`, and to the published slab form (Siewert-Benoist
Eq. 4):

.. math::

   M_{\beta,\alpha}^{\rm slab}(a)
   = B_\alpha(\xi_\beta) + e^{-2a/\xi_\beta}\,A_\alpha(\xi_\beta)

when :math:`s = +1`. SymPy verifies the substitution and the
distinctness of the two forms (their pairwise difference is
:math:`-2 e^{-2R/\xi_\beta} A_\alpha(\xi_\beta)`, identically
non-zero for :math:`R, A_\alpha > 0`).

.. _fn-method-V-fn-sphere-fn-3:

V_fn-sphere-fn.3 — Sphere bare-critical = det M(R) = 0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_critical_condition`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_3_critical_condition`.

The sphere F_N collocation system :math:`M(R)\,\vec a = 0` is
homogeneous (no source for the bare-critical problem); a non-trivial
solution exists iff :math:`\det M(R) = 0`. The determinantal form is
the same as slab (V_fn-slab.5), only the matrix entries differ via
``geometry_sign``. SymPy verifies the 1×1 + 2×2 worked examples and
confirms slab/sphere distinctness — the same bracketing scan is used
in both Branch-2 production solvers, with only the
``geometry_sign=-1`` argument changing.

.. _fn-method-V-fn-sphere-fn-4:

V_fn-sphere-fn.4 — Siewert-Thomas 2G→1G reduction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_sphere_2g_to_1g_reduction`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_4_2g_to_1g_reduction`.

Siewert-Thomas 1986 develops F_N for the general 2G case. The 1G
specialisation collapses every 2×2 matrix to a scalar:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Object
     - 2G form
     - 1G specialisation
   * - Number of secondaries
     - :math:`C` (2×2 matrix)
     - :math:`c` (scalar)
   * - Anisotropy factor
     - :math:`\Theta(\mu)` (2×2)
     - :math:`1` (scalar; isotropic)
   * - Dispersion function
     - :math:`\Lambda(z)` (2×2)
     - :math:`1 - cz\,\mathrm{atanh}(1/z)`
   * - F_N system size
     - :math:`(2N+2)^2`
     - :math:`(N+1)^2`

SymPy verifies the matrix-dimension collapse on a worked 2×2 → scalar
example: setting the off-diagonal scattering to zero and the diagonal
:math:`C_{ii} = c` recovers the scalar dispersion function and the
:math:`(N+1)`-sized F_N system. This is the formal justification for
using a *2G paper* as the spec for a 1G solver.

.. _fn-method-V-fn-sphere-fn-5:

V_fn-sphere-fn.5 — Wiener-Hopf X-function geometry-independence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_sphere_derivations.derive_x_function_geometry_independence`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_sphere.test_v_fn_sphere_fn_5_x_function_geometry_independence`.

The Wiener-Hopf X-function (Case 1960; Case-Zweifel 1967 § 4)

.. math::
   :label: fn-x-function

   X(z) = (1 - z)^{-1}\,
          \exp\!\left[\frac{1}{\pi}
                       \int_0^1 \arg\Lambda^{+}(\tau)\,
                       \frac{d\tau}{\tau - z}\right]

.. (vv-status rationale) governing: Wiener-Hopf X-function definition (Case 1960; Case-Zweifel 1967 § 4); verified geometry-independent by V_fn-sphere-fn.5 (test_v_fn_sphere_fn_5_x_function_geometry_independence). Numerical accuracy at the Branch-2 evaluator (fn_method.core.x_function) is the responsibility of the underlying scipy.integrate.quad / mpmath.quad backend.
.. vv-status: fn-x-function documented

depends only on the dispersion function

.. math::

   \Lambda(z) = 1 - cz\,\mathrm{atanh}(1/z) ,

which is a **medium-only quantity** (depends on :math:`c` and
:math:`z`, no geometry parameter). SymPy verifies that the integrand
has no :math:`R` or :math:`a` symbols. This justifies reusing the
slab X-function machinery for sphere verbatim — the X-function is
shared between :ref:`theory-fn-method` and
:ref:`theory-singular-eigenfunction` because both methods depend on
the **same** medium-property dispersion function. This is the
classic case of "shared upstream library OK as long as it is
trusted-library-line" (see ``algebra-of-record`` § "Structural
independence applies above the trusted-library line"). The
X-function evaluator is in :mod:`...fn_method.core.x_function`.

Cross-check claim — second slice
==================================

The second slice's load-bearing cross-check claims:

* **Slab Ua-1-0-SL**: F_N slab solver and Variant α slab solver
  agree on :math:`k_{\rm eff} = 1` at the Sood truth :math:`a_c =
  0.93772556` mfp to ≤ 5e-5. F_N at :math:`N = 10` reaches ~5e-6;
  Variant α slab at :math:`(n_x, n_\mu) = (48, 128)` reaches ~1e-5.
  Cross-check tolerance 5e-5 is the safe envelope. Foundation-test
  gate:
  :func:`tests.derivations.test_fn_la13511_slab_xverif.test_fn_slab_vs_variant_alpha_at_sood_ua_1_0_sl`.

* **Sphere Ua-1-0-SP**: F_N sphere (Siewert-Thomas 1986) returns
  :math:`R_c = 2.4248249802` mfp to ≤ 1e-5 (achieved 3.6e-8 at
  :math:`N = 10`); Variant α sphere at the F_N predicted radius
  gives :math:`k_{\rm eff} = 1` to ≤ 1e-5 (achieved 4.2e-6).
  Foundation-test gate:
  :func:`tests.derivations.test_fn_la13511_sphere_xverif.test_fn_sphere_vs_variant_alpha_sphere_at_sood_ua_1_0_sp`.

Together these establish the **structural-independence pillar** for
the Variant α slab + sphere prototypes against published-method
references that produced the Sood/KLL truth values. The two methods
are mathematically disjoint: F_N method uses Case singular
eigenfunctions + Wiener-Hopf factorisation; Variant α uses
bouncing-trajectory operator with rank-2 closure. Their agreement at
the same physics is the strongest L1 evidence currently available in
ORPHEUS for either method.

The sphere upgrade from PS-1982 wrapper (previous slice) to true F_N
(Siewert-Thomas 1986) makes the cross-check **structurally stronger**:
PS-1982 and Variant α both reduce the Peierls integral equation by
different algebraic paths (procedurally independent only), whereas
F_N method works in the Case singular-eigenfunction representation
and never reduces to an integral equation in :math:`r` (genuinely
structurally independent above the trusted-library line).

Flux reconstruction — KLL 1974 interior-flux extension
=======================================================

.. _fn-method-flux-recon-overview:

Overview — Phase B2 rich-machinery extension
----------------------------------------------

The F_N solvers (:mod:`...slab.one_group`, :mod:`...sphere.one_group`)
return only the **boundary** angular-flux representation
:math:`\psi(\pm a, \mp\mu) = \sum_{\alpha=0}^{N} a_\alpha \mu^\alpha`
plus the critical dimension :math:`r_c`. The Phase B2 extension
(:mod:`...slab.flux_reconstruction`, :mod:`...sphere.flux_reconstruction`)
adds **interior** scalar flux :math:`\phi(z)` / :math:`\phi(r)` and
**interior** angular flux :math:`\psi(z, \mu)` / :math:`\psi(r, \mu)`
at any point.

The interior representation comes from Kaper-Lindeman-Leaf 1974
(KLL) — the canonical published interior-flux benchmark for the
bare-critical slab and sphere. KLL's recipe is **structurally
independent of F_N** (Fredholm integral equation in :math:`A(\nu)`
solved by Wiener-Hopf factorisation, vs F_N's half-range
collocation), giving a cross-check above the trusted-library line.

Branch-1 SymPy module:
:mod:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations`.

The discipline used:

1. F_N gives :math:`r_c` and surface :math:`\psi(\pm a, \mp\mu)`.
2. KLL Fredholm iteration (structurally independent of F_N) gives
   interior :math:`\phi(z)`, :math:`\phi(r)`.
3. BTE characteristic integration (closed-form, given converged
   :math:`\phi`) gives interior :math:`\psi(z, \mu)`,
   :math:`\psi(r, \mu)`.
4. Closure :math:`\phi = \int \psi\,d\mu` self-consistency
   validates the chain.

Each step is verified by a foundation-tagged SymPy module + a
numerical L1 gate against KLL Tables III + VII (the canonical
published benchmarks).

.. _fn-method-flux-recon-V_fn-flux-slab-1:

V_fn-flux-slab.1 — KLL Eq. 7 slab scalar-flux structure
---------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_slab_kll_phi_eq7_structure`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_flux_symbolic.test_v_fn_flux_slab_1_kll_eq7_structure`.

KLL Eq. 7 (the :math:`c > 1` critical-slab scalar-flux
reconstruction) has the structure

.. math::
   :label: kll-1974-slab-phi

   \phi(x) = a\!\left[\cos(x/u_0) + \int_0^1 A(\nu)\,
            e^{-b/\nu}\,\cosh(x/\nu)\,d\nu\right]

.. (vv-status rationale) governing: KLL 1974 Eq. 7 — the bare-critical slab interior scalar flux; verified by V_fn-flux-slab.1 (test_v_fn_flux_slab_1_kll_eq7_structure) symbolically; verified numerically at L1 against KLL Table III (kll-1974-slab-flux label).
.. vv-status: kll-1974-slab-phi documented

where :math:`u_0 = |\nu_0|` is the discrete dispersion-root
magnitude, :math:`b = a/\Sigma_t` is the half-thickness in mfp, and
:math:`A(\nu)` is the continuum amplitude function (the Fredholm
unknown). Both the cosine and :math:`\cosh`-integral terms are even
in :math:`x`, so the bare-critical eigenmode symmetry
:math:`\phi(-x) = \phi(x)` follows from the formula structure (not
from the specific :math:`A(\nu)`). At :math:`x = 0`, the formula
reduces to :math:`\phi(0) = a\,[1 + \int_0^1 A(\nu)\,
e^{-b/\nu}\,d\nu]`, which is the denominator of all flux-ratio
computations.

.. _fn-method-flux-recon-V_fn-flux-slab-2:

V_fn-flux-slab.2 — :math:`\phi(z)/\phi(0)` is normalisation-free
------------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_slab_phi_endpoint_normalization`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_flux_symbolic.test_v_fn_flux_slab_2_endpoint_normalization`.

The multiplicative constant :math:`a` cancels in the ratio
:math:`\phi(z)/\phi(0)`, so the published Sood Table 14 / KLL Table
III ratios are directly computable from the converged Fredholm
solution :math:`A(\nu)` without fixing the normalisation. This is
why the Branch-2 production solver does not need a normalisation
choice; the bare-critical eigenmode is unique only up to a global
multiplicative scalar, and every observable of interest is a ratio.

.. _fn-method-flux-recon-V_fn-flux-slab-3:

V_fn-flux-slab.3 — Interior :math:`\psi(z, \mu)` via characteristics
---------------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_slab_psi_from_phi_characteristic`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_flux_symbolic.test_v_fn_flux_slab_3_psi_from_phi_characteristic`.

Once :math:`\phi(z)` is known via KLL, the BTE
:math:`\mu \partial_z \psi + \psi = (c/2)\phi(z)` integrates along
characteristics with vacuum BC :math:`\psi(\mp b, \pm\mu) = 0`
(:math:`\mu > 0`) to give the closed form

.. math::

   \psi(z, \mu > 0) = \frac{c}{2\mu}
                      \int_{-b}^{z} \phi(z')\,
                      e^{-(z-z')/\mu}\,dz'

(and the symmetric version for :math:`\mu < 0`). SymPy verifies
this expression satisfies the BTE + BC by direct substitution.

.. _fn-method-flux-recon-V_fn-flux-sphere-1:

V_fn-flux-sphere.1 — KLL Eq. 15 sphere scalar-flux structure
--------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_sphere_kll_phi_eq15_structure`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_flux_symbolic.test_v_fn_flux_sphere_1_kll_eq15_structure`.

KLL Eq. 15

.. math::
   :label: kll-1974-sphere-phi

   \phi(r) = \frac{2 b}{r}\!\left[u_0 \sin(r/u_0)
            - \int_0^1 B(\nu)\,e^{-R/\nu}\,\sinh(r/\nu)\,d\nu\right]

.. (vv-status rationale) governing: KLL 1974 Eq. 15 — the bare-critical sphere interior scalar flux; verified by V_fn-flux-sphere.1 (test_v_fn_flux_sphere_1_kll_eq15_structure) symbolically; verified numerically at L1 against KLL Table VII (kll-1974-sphere-flux label).
.. vv-status: kll-1974-sphere-phi documented

is the sphere analog of slab Eq. 7. The :math:`r \to 0` limit is
well-defined via :math:`\sin(r/u_0)/r \to 1/u_0` and
:math:`\sinh(r/\nu)/r \to 1/\nu`. The :math:`r/\nu`-multiplied form
(rather than just :math:`\nu`) is required because the sphere has
the geometric weight factor :math:`r` from the
:math:`r\Phi(r) \to \Psi(x)` substitution, and the :math:`\sinh`
(not :math:`\cosh`) reflects the antisymmetric (sphere) BC vs the
symmetric (slab) BC.

.. _fn-method-flux-recon-V_fn-flux-sphere-2:

V_fn-flux-sphere.2 — Sphere chord-length characteristic
----------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_sphere_psi_from_phi_characteristic`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_flux_symbolic.test_v_fn_flux_sphere_2_psi_from_phi_characteristic`.

Chord length from :math:`(r, \mu)` back to the surface is

.. math::

   s_{\rm in}(r, \mu) = r\mu + \sqrt{R^2 - r^2(1 - \mu^2)} .

The characteristic-integration formula is then

.. math::

   \psi(r, \mu) = \frac{c}{2}
                  \int_0^{s_{\rm in}}
                  \phi\!\left(\sqrt{r^2 - 2 r s' \mu + s'^2}\right)\,
                  e^{-s'}\,ds' .

The radial argument
:math:`\sqrt{r^2 - 2 r s' \mu + s'^2}` is the position at distance
:math:`s'` along the inward-pointing characteristic from
:math:`(r, \mu)`. SymPy verifies the radial argument's geometric
construction and the BTE satisfaction.

.. _fn-method-flux-recon-V_fn-flux-shared-1:

V_fn-flux-shared.1 — Universal angular-flux closure
-----------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_flux_reconstruction_derivations.derive_scalar_flux_angular_integral`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_flux_symbolic.test_v_fn_flux_shared_1_scalar_from_angular`.

:math:`\phi(z) = \int_{-1}^{1} \psi(z, \mu)\,d\mu` is the literal
definition of scalar flux from angular flux. SymPy verifies on three
example trial :math:`\psi`: constant, trigonometric, polynomial.
This is the **closure-consistency test** — the L1 numerical gate
verifies that the reconstructed :math:`\psi(z, \mu)` integrates back
to the converged :math:`\phi(z)` to within quadrature accuracy.

L1 verification gates
----------------------

The verification labels for these gates are:

.. math::
   :label: kll-1974-slab-flux

   \phi_{\rm slab}(z)/\phi_{\rm slab}(0) = \frac{\cos(z/u_0) +
   \int_0^1 A(\nu)\,e^{-b/\nu}\,\cosh(z/\nu)\,d\nu}
   {1 + \int_0^1 A(\nu)\,e^{-b/\nu}\,d\nu}

.. math::
   :label: kll-1974-sphere-flux

   \phi_{\rm sphere}(r)/\phi_{\rm sphere}(0) =
   \frac{(2/r)[u_0 \sin(r/u_0) - \int_0^1 B(\nu)\,e^{-R/\nu}\,
   \sinh(r/\nu)\,d\nu]}{2[1 - \int_0^1 (B(\nu)/\nu)\,
   e^{-R/\nu}\,d\nu]}

The Branch-2 implementation
(:mod:`...slab.flux_reconstruction`, :mod:`...sphere.flux_reconstruction`)
agrees with the KLL benchmark tables at the following tolerances:

* **Slab Table III at c=1.30 (Sood Ua-1-0-SL)**: ≤ 1e-5 absolute on
  :math:`\phi(z)/\phi(0)` at :math:`z/b \in \{0.25, 0.5, 0.75,
  1.0\}` (achieved 4e-7 at :math:`n_{\rm nodes} = 64`).
* **Slab Table III at c \in [1.05, 1.60]**: ≤ 5e-5 absolute on the
  full sweep.
* **Sphere Table VII at c=1.30 (Sood Ua-1-0-SP)**: ≤ 1e-5 absolute
  on :math:`\phi(r)/\phi(0)` at :math:`r/R \in \{0.25, 0.5, 0.75,
  1.0\}` (achieved 2e-8 at :math:`n_{\rm nodes} = 64`).
* **Closure check**: :math:`\int_{-1}^{1} \psi\,d\mu = \phi` to
  ≤ 1e-5 absolute (interior points; near-boundary asymptotically
  the characteristic integral becomes more singular and the
  achievable closure precision is ~1e-5).

.. _fn-method-atkinson-product-nystrom:

Atkinson product-Nyström — the Path A.i hardening (ERR-036)
=============================================================

**This section is the load-bearing narrative for the Atkinson
1972/1997 product-Nyström hardening of the Path A.i flux
reconstruction. It documents what was tried, what failed, and why
Atkinson product-Simpson was the textbook fix.** This is the
canonical worked example for the "Signature 6: log-singular kernel
diagonal truncation" entry in the
``numerical-bug-signatures`` skill.

The Path A.i operator
----------------------

The bare-critical 1G isotropic slab Peierls integral equation is

.. math::
   :label: peierls-slab-bare-critical

   \phi(z) = \frac{c}{2} \int_{-a}^{a} E_1(\Sigma_t |z - z'|)\,
   \phi(z')\,dz', \qquad z \in [-a, a],

.. (vv-status rationale) governing: The bare-critical slab Peierls integral equation; the entire Path A.i + KLL Path B + Atkinson hardening verifies this operator (test_atkinson_product_nystrom + test_path_ai_legacy_plain_gl_signature).
.. vv-status: peierls-slab-bare-critical documented

with vacuum boundaries (the :math:`-a, +a` limits make the
right-hand side a homogeneous Fredholm operator of the second kind
on :math:`L^\infty[-a, a]`). The eigenvalue is exactly 1 by the
criticality condition.

The Path A.i algorithm power-iterates the discrete Peierls operator
to extract the dominant eigenmode :math:`\phi(z)` directly. The
Path B (KLL Wiener-Hopf) algorithm reaches the same eigenmode via
:math:`A(\nu)`-Fredholm iteration. The two paths are *procedurally*
distinct — different code paths reaching the same bare-critical
mode — but they share the underlying integral kernel
:math:`(c/2) E_1`.

.. _fn-method-path-ai-legacy-failure:

What the legacy plain-GL Path A.i did, and how it failed
---------------------------------------------------------

The legacy ``slab_scalar_flux_fn_projection`` discretised the kernel
via plain Gauss-Legendre on the :math:`(z, \mu)` tensor product:

.. math::

   K[i, j]^{\rm legacy}
   = \frac{c}{2} \sum_k \frac{w_{\mu_k}}{\mu_k}\,
                          e^{-|z_i - z_j| / \mu_k}

with all-positive :math:`w_{\mu_k}`. After the :math:`\mu`-integral,
this should equal :math:`(c/2)\,E_1(|z_i - z_j|)`, where
:math:`E_1(0+) = +\infty` (logarithmic divergence).

Off-diagonal (:math:`|z_i - z_j| > 0.05` mfp) the :math:`\mu`-
quadrature converges to :math:`E_1` at machine precision. **At the
diagonal** (:math:`z_i = z_j`), the integrand
:math:`(1/\mu) e^0 = 1/\mu` is itself singular at :math:`\mu \to 0^+`;
with finite-N Gauss-Legendre nodes on :math:`(0, 1)` the smallest
node is :math:`\mu_{\rm min} \sim 1/n_\mu^2`, and the sum saturates
at the *finite* value

.. math::

   \sum_k \frac{w_k}{\mu_k} \approx -\log(\mu_{\rm min})
                               \approx 2 \log(n_\mu),

a finite-but-wrong truncation of the true :math:`+\infty`. The
discrete kernel is therefore **qualitatively wrong at the diagonal**:
a finite-by-:math:`n_\mu`-dependent amount that does not vanish with
refinement.

**Impact** on the bare-critical 1G isotropic slab (Wave 2-A cases
Ua-1-0-SL c=1.30, PUa-1-0-SL c=1.50, PUb-1-0-SL c=1.40), the
discrete eigenmode :math:`\phi(z)/\phi(0)` deviates from the KLL
Path B reference (Wiener-Hopf factorisation, converged to 1e-10) by:

============   =================================
``z/a``        rel err at ``n_quad_z=128``
============   =================================
0.0            0 (by normalization)
0.5            2.0 %
0.9            3.4 %
0.99           4.0 %
============   =================================

The error PEAKS at the slab edges (where the eigenmode's local
spatial structure is most sensitive to short-range kernel
contributions) and is roughly UNIFORM across :math:`c \in [1.3,
1.5]`.

.. _fn-method-path-ai-fingerprint:

The convergence-rate fingerprint — first-order with log correction
-------------------------------------------------------------------

The diagnostic that pinned the bug class was the
**convergence-rate fingerprint**. At :math:`z/a = 0.95`, sweeping
:math:`n_{\rm panels} \in \{32, 64, 128, 256\}`:

==========   ====================   ================
``n``        ``err·n/log(n)``       ``err·sqrt(n)``
==========   ====================   ================
32           1.39                   0.85
64           1.17                   0.61
128          1.05                   0.45
256          1.07                   0.37
==========   ====================   ================

* **`err·n/log(n) ≈ const`** across the 4× n range (1.4 → 1.0,
  stable). This is the textbook signature of **log-singular kernel
  diagonal truncation** (Atkinson 1972 § III; Atkinson 1997
  § 4.2.1).
* **`err·sqrt(n)` decays by ~3×** across the same range. This rules
  out **Schneider** :math:`C^{(0,1)}` **endpoint singularity in the
  solution** (which would give a stable :math:`err \sim 1/\sqrt{n}`
  asymptote).

The discriminator is essential because both bug classes look
similar at low-:math:`n` (1–10 % error, slowly improving) — the
rate fingerprint is the only cheap distinction, and it points
unambiguously at the kernel-diagonal-truncation cause.

The literature memo
``scratch/derivations/peierls_log_singular/atkinson_product_nystrom.md``
predicted a half-order Schneider rate (-1/2); the empirical
dominant rate is first-order with log correction (-0.9). The memo's
recommendation (Atkinson product-Nyström) was right but for a
different reason than the memo emphasised — the diagonal-truncation
fix dwarfs the endpoint-regularity issue at practical :math:`n`.
Graded mesh would help further but is not needed to hit the F_N
moment floor.

The fix — Atkinson 1972 / 1997 § 4.2 product-Simpson rule
----------------------------------------------------------

With the standard small-argument expansion (Abramowitz–Stegun
5.1.11)

.. math::
   :label: e1-small-tau-expansion

   E_1(\tau) = -\gamma_E - \log\tau + \tau - \frac{\tau^2}{4}
              + \frac{\tau^3}{18} - \cdots, \qquad \tau \to 0^+,

.. (vv-status rationale) derivation: Abramowitz–Stegun 5.1.11 small-argument expansion for E_1; standard mathematical identity. The smooth-remainder R(τ) = E_1(τ) + log τ used in the Atkinson decomposition is verified at L0 (test_E1_smooth_remainder_taylor in test_atkinson_product_nystrom).
.. vv-status: e1-small-tau-expansion documented

the kernel decomposes as

.. math::
   :label: peierls-kernel-decomposition

   \frac{c}{2}\,E_1(|z - z'|)
   = \underbrace{-\frac{c}{2}\,\log|z - z'|}_{\text{singular}}
   + \underbrace{\frac{c}{2}\,\bigl[-\gamma_E + R(|z - z'|)\bigr]}
              _{\text{smooth},\;C^\infty},

.. (vv-status rationale) derivation: Algebraic split of the Peierls log-singular kernel into singular log part + smooth remainder; verified at L0 by test_log_decomposition_foundation in test_path_ai_legacy_plain_gl_signature.
.. vv-status: peierls-kernel-decomposition documented

where :math:`R(\tau) = E_1(\tau) + \log\tau` is the
:math:`C^\infty` Taylor remainder (analytic, with :math:`R(0) =
-\gamma_E`).

Atkinson's product-Simpson rule applies as follows:

1. **The singular log piece** is integrated **analytically** against
   the piecewise-quadratic (Simpson) basis on each panel,
   eliminating the log singularity from the discrete operator. This
   is the load-bearing step.
2. **The smooth remainder** is integrated with **standard Simpson**
   on the same nodes — there is no singularity to handle.
3. The two contributions are summed at matrix-assembly time. The
   discrete operator never samples :math:`E_1` at :math:`\tau = 0`.

Closed-form weights F_k(s; t)
------------------------------

The product-Simpson rule against the log kernel needs the
antiderivatives :math:`F_k(s; t) = \int s^k \log|t - s|\, ds` for
:math:`k = 0, 1, 2`. Each integral is constructed by integration by
parts:

.. math::
   :label: fn-Fk-integration-by-parts

   \int s^k \log|t - s|\, ds
     = \frac{s^{k+1}}{k+1} \log|t - s|
       - \frac{1}{k+1} \int \frac{s^{k+1}}{s - t}\, ds .

.. (vv-status rationale) derivation: Integration-by-parts step in the F_k closed-form derivation; standard calculus identity.
.. vv-status: fn-Fk-integration-by-parts documented

Polynomial division of :math:`s^{k+1}/(s-t)` reduces the second
integral to elementary terms plus :math:`t^{k+1} \log|s-t|`, which
combines with the first term. After collection:

.. math::
   :label: fn-Fk-closed-forms

   F_0(s; t) &= (s - t)\,\log|s - t| - s, \\
   F_1(s; t) &= \tfrac{1}{2}(s^2 - t^2)\,\log|s - t|
                - \tfrac{1}{4} s^2 - \tfrac{1}{2} t s, \\
   F_2(s; t) &= \tfrac{1}{3}(s^3 - t^3)\,\log|s - t|
                - \tfrac{1}{9} s^3 - \tfrac{1}{6} t s^2
                - \tfrac{1}{3} t^2 s.

.. (vv-status rationale) governing: Closed-form antiderivatives F_k(s;t) = ∫sᵏ log|t-s|ds — the load-bearing primitives for Atkinson product-Simpson; verified against SymPy in test_F_k_log_primitives_match_sympy and against scipy.integrate.quad with explicit singularity points in test_F_k_integrals_match_scipy_quad.
.. vv-status: fn-Fk-closed-forms documented

Each :math:`F_k` is finite at :math:`s = t` because
:math:`\lim_{u \to 0} u\log|u| = 0`. The implementation uses that
limit in the form of a guard: when :math:`|s - t|` is smaller than
the machine-epsilon scale, the :math:`\log|\cdot|`-prefactored terms
are taken as zero by the limit, leaving the polynomial part.

The closed forms are verified against SymPy's
:func:`sympy.integrate(s**k * log(abs(t - s)), s)` in
:func:`tests.derivations.test_atkinson_product_nystrom.test_F_k_log_primitives_match_sympy`.

Product-Simpson weight construction
------------------------------------

The product-Simpson weights for
:math:`\int_a^b \log|t - s|\,P(s)\, ds` against the unique quadratic
:math:`P` interpolating its three node values at
:math:`(a, m, b)` with :math:`m = (a+b)/2` and :math:`h = (b-a)/2`
are linear combinations of the :math:`F_k`:

.. math::
   :label: fn-product-simpson-weights

   w_a &= \frac{1}{2 h^2}\bigl[F_2 - (m + b) F_1 + m b F_0\bigr], \\
   w_m &= -\frac{1}{h^2}\bigl[F_2 - (a + b) F_1 + a b F_0\bigr], \\
   w_b &= \frac{1}{2 h^2}\bigl[F_2 - (a + m) F_1 + a m F_0\bigr],

.. (vv-status rationale) governing: Atkinson product-Simpson weights for ∫log|t-s|·P(s)ds against the unique quadratic interpolating P(a),P(m),P(b); exact for all quadratics; verified at L0 in test_product_simpson_log_weights_exact_for_quadratics.
.. vv-status: fn-product-simpson-weights documented

where :math:`F_k = F_k(b; t) - F_k(a; t)` (definite integrals over
the panel). These weights are exact for every quadratic :math:`P` —
no quadrature error on the singular log piece.

Setup cost is :math:`O(N^2)` total (every test node looks at every
panel) — the same complexity as the plain-GL operator, with no
runtime quadrature. The singular weights are closed-form via
:func:`...peierls_atkinson_nystrom._F_k_log_primitives` and the
smooth weights are three function evaluations of
:func:`scipy.special.exp1` per panel per test node.

Convergence achieved
---------------------

.. list-table:: Atkinson product-Simpson convergence on Sood
   ``Ua-1-0-SL`` (c=1.30 bare-critical slab, :math:`\sup_z |\phi/\phi_0
   - \phi^{\rm KLL}/\phi^{\rm KLL}_0|`)
   :header-rows: 1
   :widths: 22 26 26 26

   * - :math:`n_{\rm panels}`
     - Legacy plain-GL
     - Atkinson product-Simpson
     - Improvement factor
   * - 32
     - 7 × 10⁻²
     - 5 × 10⁻³
     - 14×
   * - 64
     - 5 × 10⁻²
     - 3.5 × 10⁻⁴
     - 140×
   * - 128
     - 4 × 10⁻²
     - 5 × 10⁻⁵
     - 800×
   * - 256
     - 3 × 10⁻²
     - 1.1 × 10⁻⁵
     - 3000×

* :math:`n_{\rm panels} = 64` reaches sup_err ≤ 5e-4 (~100×
  improvement vs legacy 5e-2).
* :math:`n_{\rm panels} = 256` reaches ≤ 1.1e-5 (the F_N moment
  floor of the surface-coefficient solver — Atkinson is now
  rate-limited by F_N, not by the kernel).

The convergence rate fingerprint is also flipped: Atkinson hits a
de Hoog–Weiss superconvergent :math:`O(h^4 \log h)` regime on
operator norm (Atkinson 1997 Eq. 4.2.83), with empirical doubling
ratios 4×–10× per :math:`n`-doubling at moderate :math:`n` —
consistent with the theoretical rate.

Pinning tests for the failure-mode signature
---------------------------------------------

The legacy plain-GL Path A.i is preserved in
:func:`...slab.flux_reconstruction.slab_scalar_flux_fn_projection` at
its honest 5–7 % tolerance — it is the *plain* baseline against which
Atkinson is compared. The failure-mode signature is pinned by
:mod:`tests.derivations.test_path_ai_legacy_plain_gl_signature`:

* ``test_log_decomposition_foundation`` — kernel decomposition holds.
* ``test_diagonal_truncation_scaling`` — pinning the
  :math:`2 \log(n_\mu)` saturation.
* ``test_off_diagonal_machine_precision`` — μ-quadrature converges
  to :math:`E_1` at :math:`\tau > 0`.
* ``test_first_order_with_log_rate`` — the
  :math:`\mathrm{err} \cdot n / \log n \approx \mathrm{const}`
  fingerprint.

All tagged ``@pytest.mark.catches("ERR-036")``.

The Atkinson fix is pinned by
:mod:`tests.derivations.test_atkinson_product_nystrom`:

* ``test_l1_atkinson_vs_kll_5e_minus_4`` (parametrised over the
  three Wave 2-A cases) — asserts ``sup |err| ≤ 5e-4`` at
  ``n_panels=64``.
* ``test_l1_atkinson_high_n_hits_fn_moment_floor`` — asserts
  ``sup |err| ≤ 5e-5`` at ``n_panels=256``.
* ``test_l1_atkinson_convergence_rate_better_than_first_order`` —
  asserts the n-doubling error ratio is ``> 2.5`` (i.e. faster than
  first order).

Why this pattern matters going forward
---------------------------------------

The same log-singular kernel structure appears in:

* **Sphere Peierls kernel** :math:`r r' / |r-r'|`-style structure
  — has a log singularity after the angular reduction; the
  ``peierls_atkinson_nystrom`` module is structured to be
  sphere-extensible.
* **Cylinder Bickley-Naylor** :math:`\mathrm{Ki}_n` integrals —
  :math:`\mathrm{Ki}_1(\tau) \sim -\tau\log\tau` near 0, log-singular
  in :math:`\tau`.
* **Atalay 1997 X-function** — log-divergent integrand near
  :math:`\nu = 1` (different mechanism — see :ref:`theory-singular-eigenfunction`
  ERR-038 and the related ERR-037 X-function divergence note).

Whenever a new ORPHEUS reference solver builds a discrete operator
from a log-singular kernel, the rule is: **decompose via
:eq:`peierls-kernel-decomposition`, integrate the log piece
analytically against the basis (:eq:`fn-product-simpson-weights`),
integrate the smooth remainder by Simpson/GL.** Plain Gauss-Legendre
on the singular piece is silently incorrect; the only diagnostic
is the convergence-rate fingerprint, which most projects do not
collect.

Reflected-slab F_N — Neshat-Maiorino 1980 extension
=====================================================

The reflected-slab F_N solver extends the Grandjean-Siewert F_N
machinery to a symmetric three-region problem (reflector ‖ core ‖
reflector) with vacuum BC on the outer reflector faces. The
extension is **trivial**: the moment recursions are identical to
the bare-slab case parametrised by the per-region :math:`c_i`.

The architectural insight (from Neshat-Maiorino 1980 § II): **the
X-function is medium-local.** Each region uses its own :math:`c_i`
independently; the coupling lives entirely in the boundary
projection equations (NM Eqs. 10-11). There is **no two-region
X-function** — the coupling is geometric, not spectral.

.. _fn-method-V-fn-slab-refl-1:

V_fn-slab-refl.1 — Reflected-slab moment recursions match bare-slab
-------------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_reflected_derivations.derive_reflected_moment_recursions_match_bare`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_reflected.test_v_fn_slab_refl_1_recursions_match_bare`.

V_fn-slab-refl.1 verifies that the :math:`A_\alpha` and
:math:`B_\alpha^{(i)}` moment recursions of the reflected-slab F_N
method (NM 1980) are **identical** to the bare-slab Grandjean-Siewert
recursions, parametrised by :math:`c_i` per region. The proof is
direct: each region's per-region X-function depends only on its own
:math:`c_i`; the recursions descend from the same long-division
identity (V_fn-slab.1 / V_fn-slab.2). This is the architectural win
that makes the reflected-slab F_N method a trivial extension of
the bare-slab F_N machinery.

.. _fn-method-V-fn-slab-refl-2:

V_fn-slab-refl.2 — NM Eq. 10/11 attenuation signs + Eq. 17 limits
-----------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_reflected_derivations.derive_reflector_attenuation_signs`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_reflected.test_v_fn_slab_refl_2_attenuation_signs`.

V_fn-slab-refl.2 verifies that the exponential signs in NM Eqs. 10-11
(:math:`e^{-\Delta/\hat\xi}` and :math:`e^{+\Delta/\hat\xi}`) are
self-consistent under the reflector-traversal interpretation, and
that NM Eq. 17 (the :math:`F_0` initial-guess for :math:`b_0`) has
the correct limits :math:`b_0 \to 0` as :math:`\Delta \to 0` (no
reflector, no return) and :math:`b_0 \to A_0/B_0^{(2)}` as
:math:`\Delta \to \infty` (infinite-reflector limit).

.. _fn-method-V-fn-slab-refl-3:

V_fn-slab-refl.3 — Critical condition NM Eq. 15 reduces to Eq. 16
-----------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_reflected_derivations.derive_critical_condition_eq15_structure`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_reflected.test_v_fn_slab_refl_3_eq15_critical_condition`.

V_fn-slab-refl.3 verifies that the critical condition NM Eq. 15
(collocated at the core Case discrete eigenvalue
:math:`\xi = \nu_0`) reduces to the closed-form NM Eq. 16

.. math::
   :label: nm1980-eq16-tau-zero

   \tau^{(0)} = -\frac{\nu_0}{2}\,\log\!\left[
                     \frac{b_0 A_0(\nu_0) - B_0^{(1)}(\nu_0)}
                          {A_0(\nu_0) - b_0 B_0^{(1)}(\nu_0)}\right]

.. (vv-status rationale) derivation: Neshat-Maiorino 1980 Eq. 16 — the F_0 closed-form critical core half-thickness; verified by V_fn-slab-refl.3 (test_v_fn_slab_refl_3_eq15_critical_condition).
.. vv-status: nm1980-eq16-tau-zero documented

at :math:`N = 0` with :math:`a_0 = 1`. SymPy verifies the
algebraic equivalence symbolically.

.. _fn-method-V-fn-slab-refl-4:

V_fn-slab-refl.4 — F_0 b_0 (NM Eq. 17) reduces from Eqs. 10-11
--------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_slab_reflected_derivations.derive_F0_initial_guess_structure`.
**Test gate:**
:func:`tests.derivations.test_fn_la13511_slab_reflected.test_v_fn_slab_refl_4_F0_initial_guess`.

V_fn-slab-refl.4 verifies that the :math:`F_0` initial-guess
:math:`b_0` formula NM Eq. 17 follows algebraically from the
:math:`N = 0` truncation of NM Eqs. 10-11 with :math:`a_0 = 1`. The
result is a 2-stream reflector-albedo coefficient — physically, the
fraction of reflector-side neutron current that returns through the
core face, in the lowest-order (single-stream) approximation.

The Branch-2 production solver
(:func:`...slab.reflected.solve_fn_slab_reflected_critical`) reproduces
the NM Table 2 Burkart 1976 "Exact" reference values to ≤ 5e-5
absolute at :math:`F_7` across the 8-case sweep over
:math:`(c_1, c_2, \Delta) \in \{1.02, 1.10, 1.30, 1.50\} \times
\{0.10, \ldots, 0.90\} \times \{0.5, 1, 2, 5\}`.

Path A.i / Path B agreement — the discrete-mode form
=====================================================

The four projection-flux SymPy verifications below close the loop on
the Path A.i (BTE phase-space) ↔ Path B (KLL Wiener-Hopf) agreement.
After the Atkinson hardening (ERR-036), Path A.i and Path B agree on
:math:`\phi(z)/\phi(0)` to within the F_N moment floor at all
:math:`z/a` for the bare-critical slab.

.. _fn-method-V-fn-proj-1:

V_fn-proj.1 — phi(z) = ∫ psi dmu universal closure
--------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_projection_flux_derivations.derive_path_ai_phi_from_psi_integral`.
**Test gate:**
:func:`tests.derivations.test_fn_projection_vs_kll_flux.test_v_fn_proj_1_phi_from_psi_closure`.

V_fn-proj.1 verifies the universal scalar-flux ↔ angular-flux closure
:math:`\phi(z) = \int_{-1}^{1} \psi(z, \mu)\,d\mu`. This holds for
any solution of the BTE independent of geometry, BC, or algorithm —
it is the literal definition. SymPy verifies on three trial
:math:`\psi`: constant, trigonometric, polynomial.

.. _fn-method-V-fn-proj-2:

V_fn-proj.2 — Characteristic propagation satisfies BTE
------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_projection_flux_derivations.derive_psi_characteristic_vacuum_bc_slab`.
**Test gate:**
:func:`tests.derivations.test_fn_projection_vs_kll_flux.test_v_fn_proj_2_characteristic_propagation`.

V_fn-proj.2 verifies that the characteristic-propagation formulas
for :math:`\psi(z, \mu)` (with vacuum BC) satisfy the BTE
identically. Used in Path A.i source iteration.

.. _fn-method-V-fn-proj-3:

V_fn-proj.3 — F_N surface-flux constraint requires non-flat phi
---------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_projection_flux_derivations.derive_fn_surface_flux_constraint`.
**Test gate:**
:func:`tests.derivations.test_fn_projection_vs_kll_flux.test_v_fn_proj_3_surface_flux_constraint`.

V_fn-proj.3 verifies that a constant interior :math:`\phi = \phi_0`
produces a surface outgoing flux

.. math::

   \psi(a, +\mu) = \frac{c\phi_0}{2}\,(1 - e^{-2a/\mu})

which is **NOT polynomial in** :math:`\mu`. Hence the F_N polynomial
:math:`\sum a_\alpha \mu^\alpha` requires a non-flat eigenmode
:math:`\phi(z)` — consistent with the symmetric critical mode
:math:`\phi(z) = \cos(z/u_0)`. This is a useful negative
verification: it rules out the trivial flat-flux solution as a
F_N eigenmode candidate.

.. _fn-method-V-fn-proj-4:

V_fn-proj.4 — Path A.i and Path B share the discrete-mode form
--------------------------------------------------------------

**SymPy derivation:**
:func:`orpheus.derivations.continuous.fn_method.origins.fn_projection_flux_derivations.derive_path_ai_path_b_same_eigenmode`.
**Test gate:**
:func:`tests.derivations.test_fn_projection_vs_kll_flux.test_v_fn_proj_4_path_ai_path_b_same_eigenmode`.

V_fn-proj.4 verifies that Path A.i (BTE phase-space iteration) and
Path B (KLL Wiener-Hopf + Fredholm) share the discrete-mode form
:math:`\phi^{(0)}(z) = \cos(z/u_0)` for the bare-critical 1G slab.
The structural-independence claim is at the *algorithm* level:
agreement on this discrete-mode shape + procedural divergence in
the continuum-correction algorithms is the L1 cross-check evidence.

After the Atkinson hardening, the empirical agreement is at the
:math:`O(10^{-5})` level — the F_N moment floor. The numerical
demonstration is the
:func:`tests.derivations.test_fn_projection_vs_kll_flux.test_l1_path_ai_atkinson_vs_path_b`
gate.

References
==========

* Atkinson, K.E. (1972). *The numerical solution of Fredholm
  integral equations of the second kind with singular kernels.*
  *Numerische Mathematik* **19**, 248–259. — Practical hybrid
  product/standard Simpson algorithm. PDF locally available at
  ``scratch/literature/Atkinson(1972)The numerical solution of
  Fredholm integral equations of the second kind with singular
  kernels.pdf``.
* Atkinson, K.E. (1997). *The Numerical Solution of Integral
  Equations of the Second Kind*, Cambridge University Press. § 4.2
  "Product integration methods", pp. 116–135. PDF locally available
  at ``scratch/literature/Atkinson - The numerical Solution of
  Integral Equations of the Second Kind/``.
* de Hoog, F. & Weiss, R. (1973). *Asymptotic expansions for product
  integration.* *Math. Comp.* **27**, 295–306. — Superconvergence
  rate :math:`O(h^4 \log h)` for product Simpson on log kernels.
* Schneider, C. (1979). *Product integration for weakly singular
  integral equations.* *Math. Comp.* **36**, 207–213. — Solution
  regularity + graded mesh.
* Neshat, K., Maiorino, J.R. (1980). *The F_N method for solving
  the critical problem for a slab with a finite reflector.* *Annals
  of Nuclear Energy* **7**, 79–81. **Reflected-slab F_N method
  specification**. PDF locally available at
  ``scratch/literature/Neshat-Maiorino(1980)The FN method for
  solving the critial problem for a slab with a finite reflector.pdf``.
* Burkart, A.R. (1976). *Trans. Am. Nucl. Soc.* **24**, 190.
  "Exact" reference values cited in NM Table 2 (matched by F_7 to
  all printed digits).
* Sood, A., Forster, R.A., Parsons, D.K. (1999). *Analytical
  Benchmark Test Set for Criticality Code Verification.* Los
  Alamos National Laboratory report LA-13511. PDF at
  ``scratch/literature/Sood Foster Parsons (1999)Analytical
  Benchmark Test Set for Criticality Code Verification.pdf``.
* Sood, A., Forster, R.A., Parsons, D.K. (2003). *Analytical
  Benchmark Test Set for Criticality Code Verification.*
  *Progress in Nuclear Energy* **42**, 55. (Journal-published
  condensation; verified to be a TEST SET not a method paper —
  see ``.claude/agent-memory/literature-researcher/sood_2003_vs_1999_extraction.md``.)
* Kaper, H.G., Lindeman, A.J., Leaf, G.K. (1974). *Benchmark
  Values for the Slab and Sphere Criticality Problem in One-Group
  Neutron Transport Theory.* *Nuclear Science and Engineering*
  **54**, 94. **Used for Cases 2 + 4** (slab + sphere
  bare-critical), AND for the interior flux reconstruction
  (KLL Tables III + VII). PDF locally available at
  ``scratch/literature/Kaper-Lindeman-Leaf(1974)Benchmark Values
  for the Slab and Sphere Criticality Problem in One-Group Neutron
  Transport Th.pdf``.
* Siewert, C.E., Benoist, P. (1979). *The F_N Method in Neutron-
  Transport Theory. Part I: Theory and Applications.* *Nuclear
  Science and Engineering* **69**, 156. **Slab F_N method
  specification**, used for Case 2.
* Grandjean, P., Siewert, C.E. (1979). *The F_N Method in Neutron-
  Transport Theory. Part II: Applications and Numerical Results.*
  *Nuclear Science and Engineering* **69**, 161. **Slab F_N
  numerical-results tables**, including the bare-critical
  thickness Table XI used as the second cross-check anchor for
  Case 2.
* Pomraning, G.C., Siewert, C.E. (1982). *On the integral form of
  the equation of transfer for a homogeneous sphere.* *J. Quant.
  Spec. Rad. Transfer* **28**, 503. **Spherical Peierls integral
  equation** — alternative kernel derivation (used in the previous
  PS-1982-wrapper sphere F_N stub, now superseded by the true F_N
  method per Siewert-Thomas 1986).
* Westfall, R.M., Metcalf, D.R. (1973). *Singular Eigenfunction
  Solution of the Monoenergetic Neutron Transport Equation for
  Finite Radially Reflected Critical Cylinders.* *Nuclear Science
  and Engineering* **52**, 1. — Cylinder is in
  :ref:`theory-singular-eigenfunction`, not F_N.
* Siewert, C.E., Thomas, J.R. (1986). *On Two-Group Critical
  Problems in Neutron Transport Theory.* *Nuclear Science and
  Engineering* **94**, 264-270. **Sphere F_N method specification**
  used for Case 4. The 1G specialisation is a clean reduction of
  the 2G machinery (matrix-dimension collapse to scalars). PDF
  locally available at ``scratch/literature/Siewert-Thomas(1986)On
  Two-Group Critical Problems in Neutron-Transport Theory.pdf``.
* Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976). *Neutron
  transport in two dissimilar media with anisotropic scattering.*
  *Nuclear Science and Engineering* **61**, 72-81. — Two-region F_N
  with anisotropic scattering; reserved for a future Wave on
  reflected-slab P_N anisotropy.
* Case, K.M. (1960). *Elementary solutions of the transport equation
  and their applications.* *Annals of Physics* **9**, 1-23. — The
  foundational source for the singular-eigenfunction expansion;
  shared between F_N and singular_eigenfunction methods (see
  :ref:`theory-singular-eigenfunction`).
* Mitsis, G.F. (1963). *Transport Solutions to the Monoenergetic
  Critical Problems.* Argonne National Laboratory report ANL-6787. —
  The :math:`r\Phi(r) \to \Psi(x, \mu)` substitution that lets
  sphere F_N share slab F_N machinery via the geometry-sign flip.

Internal references:

* Literature memo:
  ``.claude/agent-memory/literature-researcher/sood_fn_method_full_extraction.md``.
* Atkinson method memo:
  ``.claude/agent-memory/literature-researcher/atkinson_product_nystrom_log_kernels.md``.
* Atkinson investigation closeout:
  ``.claude/agent-memory/numerics-investigator/atkinson_path_ai_log_kernel.md``.
* Method-implementer closeout:
  ``.claude/agent-memory/method-implementer/fn_method_kinf_first_slice.md``.
* :doc:`trajectory_resolvent` — companion Variant α reference family.
* :doc:`singular_eigenfunction` — Atalay 1997 + WM-72 reflected /
  cylinder family.
* :doc:`peierls_nystrom` — direct Peierls-integral reference solver
  (the third pillar after F_N and Variant α).
