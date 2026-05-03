.. _theory-singular-eigenfunction:

================================================================================
Singular Eigenfunction Expansion (Mitsis-Westfall-Metcalf family)
================================================================================

.. note::

   **Stub-grade page (Phase B1 hardened).** The archivist agent will
   expand each TODO marker into the full narrative + numerical evidence
   + literature backstory.

   Scope of this page:

   * **Bare-critical infinite cylinder, 1G isotropic scattering** via
     Westfall-Metcalf 1973 (Case ``Ua-1-0-CY`` of the Sood LA-13511
     benchmark family).

   The Branch-1 SymPy verification (V_se-cyl.1 through V_se-cyl.8) is
   the **fully verified algebra-of-record** for the WM-72 derivation.
   It includes the discovery of one typo in the printed paper
   (Eq. 17: single :math:`\nu_0` should be :math:`\nu_0^2` in
   numerator) and the resolution of one structural-formula error in
   the original Phase B1 SymPy (q-formula's second-term denominator
   should be :math:`\mu`, not :math:`R`, forced by the Wronskian
   identity :math:`q(\mu, \mu) = 1`).

   The Branch-2 production solver implements the **full Mitsis-WM
   Fredholm method** (WM-72 Eqs. 30-32) with Mitsis-Zweifel singular
   subtraction + Lagrangian-derivative diagonal handling. It reproduces
   Sood ``Ua-1-0-CY`` to ≤ 3e-7 relative at :math:`n_{\rm grid} = 24`
   and matches all six WM-72 Table II configurations to the same
   precision in ≤ 0.1 s per solve.

Why a separate package from the F_N method
================================================================================

.. _theory-se-vs-fn:

.. note:: TODO — Archivist expansion needed.
   The SymPy module lives at
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations`.
   Brief: both the F_N method (slab/sphere) and the
   singular-eigenfunction expansion (cylinder) belong to the broad
   "Case 1960 / Mitsis 1963" family of methods, but the structural
   mathematics is fundamentally different. The F_N method
   (Siewert 1979-1986) builds a finite (N+1)×(N+1) determinantal
   equation from Wiener-Hopf X-functions; this approach **does not
   converge** for the bare cylinder (Westfall-Metcalf 1973 explicitly
   notes: "we found that the solution as formulated by Mitsis is not
   convergent" for the bare cylinder). WM-72's reformulation uses
   Bessel-kernel expansion + iterative Fredholm scheme — different
   mathematical machinery. The two packages share **only** the
   dispersion-root primitive (``case_nu0``), which is a medium
   property independent of geometry; everything above the
   trusted-library line is structurally distinct, satisfying the
   ``algebra-of-record`` independence requirement.

Bare-critical cylinder via WM-72
================================================================================

.. _theory-se-wm72-derivation:

.. note:: TODO — Archivist expansion needed.
   The complete WM-72 derivation chain is captured in 8 SymPy
   functions; expand each into its full narrative here.

   * V_se-cyl.1 (:func:`...origins.cylinder_derivations.derive_dispersion_function`):
     dispersion function :math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu) = 0`
     is identical to slab/sphere; reflects the medium-property nature
     of Case singular eigenfunctions.
   * V_se-cyl.2 (:func:`...origins.cylinder_derivations.derive_discrete_pseudo_eigenfunction`):
     discrete pseudo-eigenfunction :math:`\eta_0(\mu) = c\nu_0^2/(\nu_0^2 - \mu^2)`
     satisfies WM-72 Eq. 15. **Catches a typo in printed Eq. 17**
     (single :math:`\nu_0` in numerator → should be :math:`\nu_0^2`).
   * V_se-cyl.3 (:func:`...origins.cylinder_derivations.derive_bessel_wronskian_identity`):
     Bessel-Wronskian identity :math:`K_1(z)\,I_0(z) + I_1(z)\,K_0(z) = 1/z`
     used in the Eq. 9 integrodifferential reduction.
   * V_se-cyl.4 (:func:`...origins.cylinder_derivations.derive_bare_cylinder_reduction`):
     bare-cylinder reduction :math:`c_1 = c_2 \implies D(\nu) = 0`,
     :math:`d_0 = 0`, :math:`A(\nu) = B(\nu)` (NOT zero — corrected
     vs the original Phase B1 stylized SymPy).
   * V_se-cyl.5 (:func:`...origins.cylinder_derivations.derive_bare_cylinder_criticality_condition`):
     **(corrected)** bare-cylinder Fredholm criticality structure with
     the corrected q-formula
     :math:`q(\nu, \mu) = (R/\nu)\,K_0(R/\mu)\,I_1(R/\nu) +
     (R/\mu)\,K_1(R/\mu)\,I_0(R/\nu)`. The :math:`(R/\mu)` denominator
     in the second term is forced by the Wronskian identity
     :math:`q(\mu, \mu) = 1`.
   * V_se-cyl.6 (:func:`...origins.cylinder_derivations.derive_discrete_eigenfunction_normalization`):
     discrete normalisation matches WM-72 Eq. 21d.
   * V_se-cyl.7 (:func:`...origins.cylinder_derivations.derive_flux_reconstruction_bare_cylinder`):
     bare-cylinder neutron density profile is :math:`\rho(r) = J_0(r/u_0)`
     (the dominant Case eigenfunction with imaginary
     :math:`\nu_0 = i u_0` for :math:`c > 1`).
   * V_se-cyl.8 (:func:`...origins.cylinder_derivations.derive_singular_subtraction_eq31`):
     **NEW** — Mitsis-Zweifel singular-subtraction structural identity
     for Eq. 31. The load-bearing algebra behind the Branch-2
     production solver's diagonal handling: PV residue and
     :math:`\lambda \delta` collapse via the dispersion identity
     :math:`(1-\lambda) + \lambda = 1`.

   Test gates: :mod:`tests.derivations.test_singular_eigenfunction_cylinder`
   (one foundation test per ``derive_*()``, plus production-solver
   sanity tests, an L1 Sood ``Ua-1-0-CY`` reference-value gate at
   1e-5, and a parametrized L1 gate over all six WM-72 Table II
   configurations).

Production solver — full Mitsis-WM Fredholm method
================================================================================

.. _theory-se-wm72-numerics:

.. note:: TODO — Archivist expansion needed.
   The Branch-2 solver lives at
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.cylinder.one_group`.

   Brief: the production solver implements the full Mitsis-WM Fredholm
   method-of-record (WM-72 Eqs. 30-32) at the bare-cylinder reduction
   (:math:`a_0 = b_0 = 1, d_0 = 0, D(\nu) = 0, A(\nu) = B(\nu)` per
   V_se-cyl.4). The coupled equations

   .. math::

      \Phi'(\mu) = -I_0(R/\mu)\,q(\nu_0,\mu)\,\eta_0(\mu)
                - c \int_0^1 \frac{A'(\nu)\,\nu^2\,H(\nu, \mu)}{\nu + \mu}\,d\nu
      \quad\text{(Eq. 30, bare limit)}

   .. math::

      A'(\nu) = \frac{1}{N_2(\nu)} \int_0^1 \mu^2\,\eta_{2\nu}(\mu)\,\Phi'(\mu)\,d\mu
      \quad\text{(Eq. 31)}

   .. math::

      g(R) := c \int_0^1 \frac{\mu^2\,\Phi'(\mu)}{\mu^2 + u_0^2}\,d\mu = 0
      \quad\text{(Eq. 32, criticality)}

   are discretised on the same Gauss-Legendre grid for both
   :math:`\mu` and :math:`\nu` (same nodes per WM-72 numerical
   scheme). Discretisation transforms the inner Fredholm coupling
   into the linear system

   .. math::

      (\mathbb{I} + c\,M_{A\phi}\,M_{\phi A})\,\mathbf{A}'
      = M_{A\phi}\,\boldsymbol{\Phi}'_0,

   solved by ``numpy.linalg.solve`` (replacing the 1973-era Jacobi
   iteration with a faster + more accurate modern alternative).
   Brent root-find on Eq. 32's residual gives :math:`R_c`.

   **Mitsis-Zweifel singular subtraction** (V_se-cyl.8) handles the
   Cauchy P.V. + :math:`\lambda \delta` of :math:`\eta_{2\nu}` in
   Eq. 31:

   .. math::

      \int_0^1 \mu^2\,\eta_{2\nu}(\mu)\,\Phi'(\mu)\,d\mu
      = \int_0^1 \frac{c\,\nu^2\,[\mu^2 \Phi'(\mu) - \nu^2 \Phi'(\nu)]}
                       {\nu^2 - \mu^2}\,d\mu + \nu^2\,\Phi'(\nu),

   collapsing the singular content into a regular GL-quadrable
   integrand plus a single :math:`\nu^2 \Phi'(\nu)` residue. The
   diagonal point :math:`\mu_j = \nu_i` of the regular integrand is
   evaluated by Lagrangian-interpolation differentiation — exactly
   the technique WM-72 cite on p. 7: "evaluating the derivative term
   by Lagrangian interpolation over all points."

   **Scaled Bessel functions** (``i0e``, ``i1e``, ``k0e``, ``k1e``)
   throughout to avoid overflow at large :math:`R/\mu`. Exponential
   factors in :math:`I_0(R/\mu)\,K_n(R/\mu)` cancel pairwise; the
   :math:`I_n(R/\nu)\,K_n(R/\mu)` products with :math:`\nu \neq \mu`
   carry exponential :math:`e^{R/\nu - R/\mu}` factors that we keep
   in scaled form.

   **Accuracy at :math:`n_{\rm grid} = 24`** (matching WM-72's
   original quadrature order):

   ===========   =====================   ===================   ==================
   ``c``         ``R_c (mfp)`` truth     Solver result          Relative error
   ===========   =====================   ===================   ==================
   1.05          5.411288 (WM-72)        5.4112891              4e-7
   1.10          3.577391 (WM-72)        3.5773921              3e-7
   1.20          2.287209 (WM-72)        2.2872099              4e-7
   1.30          1.72500292 (Sood)       1.72500349             3e-7
   1.40          1.396979 (WM-72)        1.39697910             5e-8
   2.00          0.668613 (WM-72)        0.66861305             8e-8
   ===========   =====================   ===================   ==================

   Wall-clock time per solve: ≤ 0.1 s on a typical container CPU.

   Convergence rate is **near-spectral** for smooth integrands
   (the GL quadrature + barycentric differentiation give the
   expected exponential-rate convergence on smooth analytic
   functions). Empirical: error at :math:`n=12` is ≤ 1e-5; at
   :math:`n=24` ≤ 1e-6; at :math:`n=48` ≤ 1e-7. This is **4-6
   orders of magnitude better** than the original Phase B1 prototype
   (which used direct Nyström on Eq. 6a + single-cell product
   integration on the log-singular kernel diagonal, achieving only
   :math:`O(1/n)` algebraic convergence with a 1e-3 floor at
   :math:`n=128`).

Cross-check vs Variant α at Sood ``Ua-1-0-CY``
================================================================================

.. _theory-se-wm72-xverif:

.. note:: TODO — Archivist expansion needed.
   Test gate:
   :mod:`tests.derivations.test_singular_eigenfunction_cylinder_xverif`.

   Brief: the hardened WM-72 solver agrees with the Variant α cylinder
   solver at the Sood ``Ua-1-0-CY`` configuration to ≤ 1e-5 relative.
   Both solvers reproduce the published Sood :math:`r_c = 1.72500292`
   mfp:

   * **Variant α** at 8.5e-6 (already shipped at
     :mod:`tests.derivations.test_peierls_greens_function_cylinder_xverif_sood2003`,
     via bouncing-characteristic integration with analytical
     bounce-period summation).
   * **WM-72** at ≤ 3e-7 (this module, via singular-eigenfunction
     Fredholm coupling with Mitsis-Zweifel subtraction).

   Both methods share **only** the dispersion-root primitive
   (``case_nu0``) — a medium property independent of geometry. Above
   the trusted-library line, the methods are entirely disjoint:

   * **Variant α**: angle-resolved scalar transport in
     :math:`(r, \mu, \phi)` phase space, with bouncing characteristics
     and analytical bounce-period summation.
   * **WM-72**: scalar-density integral transport equation reduced
     via singular-eigenfunction expansion, with Bessel-kernel matrix
     Fredholm coupling between :math:`\Phi'(\mu)` and :math:`A'(\nu)`.

   Per ``algebra-of-record`` § "Structural independence applies above
   the trusted-library line", the cross-check at the same precision
   level (1e-5) is a true structurally-independent L1 anchor for the
   Sood reference value. A third leg via ``peierls_nystrom``
   (Bickley-Naylor :math:`\mathrm{Ki}_3` integrals) is available for
   future expansion.

Errata caught (V&V publication artifacts)
================================================================================

.. _theory-se-errata:

.. note:: TODO — Archivist expansion needed.

   Brief: this implementation slice caught **two errata**, both
   discovered by the algebra-of-record SymPy re-derivation
   discipline:

   * **WM-72 Eq. 17** (printed): :math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}/
     (\nu_{0l}^2 - \mu^2)`. **Direct substitution into Eq. 15 fails**:
     the LHS reduces to :math:`c\nu_0\mu^2`, the RHS reads
     :math:`c\nu_0^2\mu^2` — mismatched power of :math:`\nu_0`. The
     correct form is :math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}^2/
     (\nu_{0l}^2 - \mu^2)`, which:

     - Satisfies Eq. 15 exactly.
     - Reproduces Eq. 21d's :math:`\nu_0^4` factor on
       :math:`N_0 = \int_0^1 \mu^2\eta_0^2\,d\mu`.
     - Closes the half-range normalisation Eq. 14 under dispersion.

     SymPy V_se-cyl.2 caught this typo on first run.

   * **q-formula in Phase B1 SymPy** (V_se-cyl.5, original): the
     Phase B1 SymPy module wrote
     :math:`q(\nu, \mu) = (R/\nu)\,K_0\,I_1 + R\,K_1\,I_0`, but the
     WM-72 paper (Eq. 28 footnote) actually uses
     :math:`(R/\mu)\,K_1\,I_0` — the :math:`R` in the second term
     should be divided by :math:`\mu`. The Wronskian identity
     :math:`q(\mu, \mu) = 1` (which is structurally required by
     Eq. 29b's :math:`\nu \to \mu` limit) FORCES the :math:`(R/\mu)`
     denominator. The original B1 form gave :math:`q(\mu, \mu) \approx
     0.72`, inconsistent with Eq. 29b. The corrected form (now in
     V_se-cyl.5) closes the Wronskian identity exactly.

     This was caught by the post-hardening cross-check between
     V_se-cyl.5's algebraic claim and the Wronskian identity
     V_se-cyl.3.

   The discipline of re-deriving every published equation in SymPy
   and the discipline of cross-checking every derived identity
   against the Wronskian (or any other independent structural
   relation) are now both empirically validated as load-bearing for
   V&V correctness.

References
==========

* **Westfall & Metcalf 1973** — *Nucl. Sci. Eng.* **52**, 1-11.
  "Singular Eigenfunction Solution of the Monoenergetic Neutron
  Transport Equation for Finite Radially Reflected Critical
  Cylinders." DOI 10.13182/NSE73-A23285.
* **Mitsis 1963** — Argonne National Laboratory report ANL-6787,
  "Transport Solutions to the Monoenergetic Critical Problems."
* **Metcalf & Zweifel 1968** — *Nucl. Sci. Eng.* **33**, 318. — the
  singular-subtraction technique used in WM-72 Eqs. 31 and 33.
* **Atkinson 1976** — SIAM, *A Survey of Numerical Methods for the
  Solution of Fredholm Integral Equations of the Second Kind*,
  Chapter 6 (product integration).
* **Berrut & Trefethen 2004** — *SIAM Review* **46**, 501-517,
  "Barycentric Lagrange Interpolation." — the differentiation
  matrix used for the Eq. 31 diagonal handling.
