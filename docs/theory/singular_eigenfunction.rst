.. _theory-singular-eigenfunction:

================================================================================
Singular Eigenfunction Expansion (Case 1960 family)
================================================================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this before modifying any solver in
:mod:`orpheus.derivations.continuous.singular_eigenfunction`.**

- **What this is**: a Pillar-2 reference family realising
  Meaning (γ) in the :ref:`reference-solvers-three-meanings`
  taxonomy — the angular Green's function
  :math:`G(\tau, \tau'; \mu, \mu')` constructed via Case
  ν-spectrum + half-range completeness.
- **Pillar classification**: closed-form (criticality determinant) /
  semi-analytical (interior flux reconstruction via KLL 1974
  Fredholm iteration; the latter is currently in
  :mod:`orpheus.derivations.continuous.fn_method` while interior
  reconstruction is in production for slab + sphere via KLL 1974;
  cylinder reconstruction is the open extension).
- **Geometry × anisotropy coverage**: cylinder isotropic
  (Westfall-Metcalf 1973); slab + sphere linearly anisotropic with
  Atalay 1997 parity flip; foundational primitives
  (X-function, ν₀ dispersion root, half-range projections) shared
  across geometries via :mod:`...singular_eigenfunction.core`.
- **Status**: Stub-grade theory page. The TODO markers below
  are expansion targets for the archivist; the production code +
  V&V-pinned tests already pass against the Sood LA-13511 truth
  set and the Atalay 1997 Tables 2 + 3 reference values.
- **Cross-references**: :ref:`theory-fn-method` is the structurally-
  independent collocation cross-check; :ref:`theory-galerkin-spectral`
  is the matrix-Galerkin cross-check; :ref:`theory-trajectory-resolvent`
  realises Meaning (α) on the same Sood family.


.. note::

   **Stub-grade page (consolidated 2026-05-03).** This page consolidates
   the previous ``case_method`` page (Atalay 1997 reflected slab + sphere
   with linear anisotropy) and the previous ``singular_eigenfunction``
   page (Westfall–Metcalf 1973 bare radially-reflected cylinder, isotropic
   scattering). The archivist agent will expand each TODO marker into the
   full narrative + numerical evidence + literature backstory.

   Scope of this page:

   * **Bare-critical infinite cylinder, 1G isotropic scattering** via
     Westfall-Metcalf 1973 (Case ``Ua-1-0-CY`` of the Sood LA-13511
     benchmark family).
   * **Reflected slab + linearly anisotropic scattering** — Atalay
     1997 [Atalay1997]_ Eq. 46 (even modes).
   * **Reflected sphere + linearly anisotropic scattering** — Atalay
     1997 Eq. 54 (odd modes; parity-flip of slab Eq. 46).
   * Foundation: Case singular-eigenfunction expansion + first-order
     Fredholm iteration + four NEW parallel half-range relations
     (Atalay Eqs. 28-31) that close the deficit blocking previous
     reflected-slab attempts.

.. _theory-singular-eigenfunction-consolidation:

Why this consolidation
================================================================================

This page (and the underlying
:mod:`orpheus.derivations.continuous.singular_eigenfunction` package)
merges what were historically separate ``case_method/`` (Atalay 1997
slab + sphere with linear anisotropy) and ``singular_eigenfunction/``
(Westfall–Metcalf 1973 cylinder, isotropic) folders/pages. The merge is
justified because **all primary authors of the constituent papers
explicitly call this technique "Case singular eigenfunction expansion"**
and cite Case 1960 [Case1960]_ as the foundational source:

* **Case 1960** [Case1960]_ — the foundational source. Introduces the
  *elementary solutions* (a discrete pair :math:`\pm\nu_0` plus a
  continuum on :math:`(-1, 1)`) that subsequent authors named the
  "Case singular eigenfunctions" after.

* **Mitsis 1963 (ANL-6787)** [Mitsis1963]_, abstract:
  *"Transport solutions to the monoenergetic plane, spherical, and
  cylindrical critical problems with isotropic scattering are developed
  by the method of singular expansion modes."*

* **Westfall–Metcalf 1973** [WestfallMetcalf1973]_, introduction:
  *"Since the introduction of the singular eigenfunction expansion
  technique by Case [Ref. 1] in 1960, a wide variety of transport
  problems have been treated by this method."*

* **Atalay 1997** [Atalay1997]_, abstract: *"Case's singular
  eigenfunction method is used to formulate the criticality conditions.
  In addition to available bi-orthogonality relations in the literature,
  some parallel relations are derived to obtain the solution."*

The two ORPHEUS folders represented **two parametric variations on the
same method**, not two methods:

* :mod:`...singular_eigenfunction.cylinder` (Westfall–Metcalf 1973
  family): cylinder, **isotropic scattering**, modified-Bessel-K radial
  kernels via the addition theorem + full-range completeness theorem.
* :mod:`...singular_eigenfunction.slab` and
  :mod:`...singular_eigenfunction.sphere` (Atalay 1997 family): slab +
  sphere, **linear anisotropy**, half-range bi-orthogonality + first-
  order Fredholm iteration.

Different geometry, different scattering anisotropy, different kernel
reduction (exponential :math:`E_n`'s vs Bessel :math:`K_0`'s), but
**identical mathematical machinery** above the trusted-library line:
discrete eigenvalue :math:`\nu_0` (same dispersion function),
continuum modes on :math:`(-1, 1)`, half-range orthogonality, Fredholm
reduction of the boundary condition. The consolidation also cures the
previous ``case_method/`` folder's violation of the project no-author-
folders rule.

The dispersion function
:math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu)` is shared with
:mod:`...fn_method` via :func:`...fn_method.core.dispersion.case_nu0` —
a *medium* property, acceptable cross-package reuse below the trusted-
library line. Above that line, the F_N method and this package are
structurally independent (see "Why a separate package from the F_N
method" below).

.. _theory-se-vs-fn:
.. _theory-case-vs-fn-method:

Why a separate package from the F_N method
================================================================================

.. note:: TODO — Archivist expansion needed.

   The SymPy modules live at
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations`
   (cylinder)
   and
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`
   (slab + sphere, Atalay).

   Brief: the F_N method (Siewert–Benoist 1979 + Grandjean–Siewert
   1979 for slab; Siewert–Thomas 1986 for sphere) and this package's
   singular-eigenfunction methods (Atalay 1997 slab/sphere; WM-72
   cylinder) belong to the same broad "Case 1960 / Mitsis 1963 /
   McCormick–Kušcer 1965" family, but the structural mathematics
   is fundamentally different.

   * **F_N (slab/sphere).** Imposes the exact half-space exit-
     distribution equation at :math:`(N+1)` collocation points and
     solves a determinantal equation. Slab and sphere are unified by
     a geometry sign :math:`s = \pm 1`.

   * **Atalay (slab/sphere).** Uses the **full normal-mode expansion**
     (discrete + continuum) and reduces the boundary-condition closure
     to a Fredholm integral equation iterated to first order. The two
     methods access the same Wiener-Hopf X-function but exercise it
     differently: F_N uses :math:`X(\xi_\beta)` at collocation points
     :math:`\xi_\beta`; Atalay uses :math:`X(\pm \nu_0)` and
     :math:`X(-\nu)` for :math:`\nu \in (0, 1)` integrated against
     half-range moments.

   * **Westfall–Metcalf (cylinder).** Westfall–Metcalf 1973 explicitly
     notes: *"we found that the solution as formulated by Mitsis is
     not convergent"* for the bare cylinder. WM-72's reformulation
     uses Bessel-kernel expansion + iterative Fredholm scheme —
     different mathematical machinery again. The F_N method is not
     applicable to the cylinder geometry without major modification.

   The two packages are **structurally independent** above the trusted-
   library line (sharing only ``numpy``, ``scipy.integrate.quad``,
   ``mpmath``).

Slab — Atalay 1997 reflected, linearly anisotropic
================================================================================

.. _theory-case-half-range:

The four new parallel half-range relations
--------------------------------------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   Brief: Atalay's load-bearing technical contribution is the
   four NEW half-range relations Eqs. 28-31. The McCormick–Kušcer
   1965 bi-orthogonality relations [McCormickKuscer1965]_ (Atalay
   Eqs. 18-21) integrate the weight :math:`[\phi_{0+}(\mu) + B c
   \nu_0/2] \gamma(\mu) (\nu_0 - \mu)` against the four half-range
   basis members. These suffice for the half-space Milne problem but
   **not** for the reflected-slab boundary condition Atalay Eq. 16,
   which requires a second weight :math:`[\phi_\nu(\mu) + c \nu/(2
   \bar\nu)] \gamma(\mu)`. Eqs. 28-31 close the deficit by integrating
   this parallel weight against the same four basis members.

   The structural parallelism is verified in
   :func:`...origins.slab_sphere_derivations.derive_atalay_half_range_eqs28_to_31`
   (V_case.3): Eqs. 28-31 share the Wiener-Hopf X-function factor
   structure (:math:`X(\pm\nu_0)`, :math:`X(\pm\nu')`, :math:`N(\nu)`)
   of the McCormick–Kušcer set.

.. _theory-case-slab-eq46:

Slab criticality (Eq 46)
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq46

   \pm\frac{\pi}{2} - \arctan\!\frac{R \sin[(d-z_0)/|\nu_0|] + \sin[(d+z_0)/|\nu_0|]}
                                       {R \cos[(d-z_0)/|\nu_0|] - \cos[(d+z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(K_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+K_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`
   (function ``derive_atalay_critical_slab_eq46``).
   Test gate:
   :func:`tests.derivations.test_case_method_symbolic.test_v_case_5_critical_slab_eq46`.
   Closeout memo:
   ``.claude/agent-memory/method-implementer/wave_2b_atalay_case_method.md``.

   Brief: Eq. 46 emerges from Eq. 43 (the boundary-condition closure
   in :math:`(R e^{\pm i a_1}, e^{\pm i a_2})`-form) by observing the
   numerator and denominator are complex conjugates — the log of
   their ratio is :math:`2i \arg(z)`, which gives the arctan form.
   The :math:`\pm \pi/2` ambiguity reflects multiple eigenvalue
   modes; the fundamental mode is the smallest positive :math:`d`.

Sphere — Atalay 1997 reflected, linearly anisotropic, parity flip
================================================================================

.. _theory-case-sphere-eq54:

Sphere criticality (Eq 54) via parity flip
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq54

   \arctan\!\frac{\sin[(d+z_0)/|\nu_0|] - R \sin[(d-z_0)/|\nu_0|]}
                {\cos[(d+z_0)/|\nu_0|] + R \cos[(d-z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(L_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+L_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`
   (function ``derive_atalay_critical_sphere_eq54_via_parity_flip``).
   Test gate:
   :func:`tests.derivations.test_case_method_symbolic.test_v_case_6_critical_sphere_eq54_via_parity_flip`
   plus the numerical parity-flip equivalence at vacuum BC
   (:mod:`tests.derivations.test_case_method_slab_sphere_parity_flip`).

   Brief: the sphere problem is treated as the **odd-mode** of the
   slab problem on :math:`[-R, R]` via the antisymmetric BC
   :math:`\psi(x, \mu) = -\psi(-x, -\mu)` (Atalay Eq. 47). The
   structural change reduces to :math:`T \to T_1` (sign flip on
   second exponential in the T-function), :math:`K_j \to L_j` (same
   integrand structure with :math:`T \to T_1`), and a sin↔cos
   shuffle in the LHS arctan argument. Numerically: at vacuum BC
   :math:`R = 0`, both :math:`T(0,\mu) = T_1(0,\mu) = e^{-2d/\mu}`,
   so :math:`K_j = L_j` bit-for-bit (verified in the parity-flip
   test).

Cylinder — Westfall–Metcalf 1973, bare radially-reflected, isotropic
================================================================================

.. _theory-se-wm72-derivation:

Bare-critical cylinder via WM-72
--------------------------------------------------------------------------------

.. note::

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

.. _theory-se-wm72-numerics:

Production solver — full Mitsis-WM Fredholm method
--------------------------------------------------------------------------------

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

.. _theory-se-wm72-xverif:

Cross-check vs Variant α at Sood ``Ua-1-0-CY``
--------------------------------------------------------------------------------

.. note:: TODO — Archivist expansion needed.
   Test gate:
   :mod:`tests.derivations.test_singular_eigenfunction_cylinder_xverif`.

   Brief: the hardened WM-72 solver agrees with the Variant α cylinder
   solver at the Sood ``Ua-1-0-CY`` configuration to ≤ 1e-5 relative.
   Both solvers reproduce the published Sood :math:`r_c = 1.72500292`
   mfp:

   * **Variant α** at 8.5e-6 (already shipped at
     :mod:`tests.derivations.test_trajectory_resolvent_cylinder_xverif_sood2003`,
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

Anisotropy — linearly anisotropic scattering primitives
================================================================================

.. _theory-case-x-function-eq40:

X-function (Atalay Eq 40)
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq40

   X(\mu) = \exp\!\Bigg\{ -\frac{c}{2} \int_0^1 d\nu\, g_1(c,\nu)\,
       \Big[d^2(\nu^2)\Big(1 + \frac{c\nu^2}{1-\nu^2}\Big)
            + 3 f_1 (1-c)^2 \nu^2 d(-\nu^2)\Big]\,\ln(\nu - \mu) \Bigg\}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`
   (function ``derive_atalay_x_function_eq40``).
   Test gate:
   :func:`tests.derivations.test_case_method_symbolic.test_v_case_8_x_function_eq40`.

.. _theory-case-extrapolated-endpoint-eq42:

Extrapolated endpoint :math:`z_0` (Atalay Eq 42)
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq42

   z_0 = -\frac{\nu_0}{2} \ln\!\frac{d(-\nu_0\bar\nu)}{d(\nu_0\bar\nu)}
        + \frac{c\,\nu_0}{4} \int_0^1 d\mu\, g_1(c,\mu)\,
                                \Big[d^2(\mu^2)\Big(1 + \frac{c\mu^2}{1-\mu^2}\Big)
                                     + 3 f_1 (1-c)^2 \mu^2 d(-\mu^2)\Big]\,
                                \ln\!\frac{\nu_0 + \mu}{\nu_0 - \mu}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`
   (function ``derive_atalay_extrapolated_endpoint_eq42``).

   Implementation note: reproduces Atalay Table 1 :math:`z_0(c, f_1=0)`
   to ≤ 1e-12 absolute after the ERR-037 fix (μ = tanh(t) substitution
   resolving the slowly-cancelling endpoint pole). See
   :mod:`tests.derivations.test_case_method_z0` for the regression
   gate.

.. _theory-case-validity-eq5:

Validity bound (Atalay Eq 5)
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq5

   c \le 1 + \frac{1}{3 f_1}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.slab_sphere_derivations`
   (function ``derive_atalay_validity_bound_eq5``).
   Test gate:
   :func:`tests.derivations.test_case_method_symbolic.test_v_case_9_validity_bound_eq5`.

   Brief: Atalay Eq. 5 is the **one-pair-of-discrete-modes** range.
   For :math:`f_1 = 0` the bound is trivial (all c); for
   :math:`f_1 = 0.30` the upper limit is :math:`c \le 19/9 \approx 2.111`.
   Outside this band complex eigenvalues appear (Dahl–Sjöstrand
   1979; Kohut 1993) and Atalay's first-order Fredholm iteration
   breaks. The slab/sphere solvers raise :class:`ValueError` if
   :math:`c` violates the bound.

Origins — Branch-1 SymPy algebra-of-record
================================================================================

The :mod:`...singular_eigenfunction.origins` subpackage hosts both
SymPy modules:

* :mod:`...singular_eigenfunction.origins.cylinder_derivations` —
  Westfall–Metcalf 1973 cylinder derivations: bare-cylinder integral-
  equation reduction, dispersion function (re-derivation matching
  Case 1960 / WM-72 Eq. 18), Bessel-Wronskian identity used in the
  integrodifferential reduction (Eq. 9), pseudo-eigenfunction
  structure (Eq. 17 + 19), bare-cylinder closure (Eq. 32 + 30
  simplification).

* :mod:`...singular_eigenfunction.origins.slab_sphere_derivations` —
  Atalay 1997 slab + sphere derivations: dispersion relation reduction
  Eq. 11 → 12 (linearly anisotropic), parity flip Eqs. 47-49 mapping
  slab-odd → sphere, half-range relations Eqs. 28-31, Fredholm-form
  prefactor Eq. 32, criticality conditions Eqs. 46 (slab) / 54
  (sphere), X-function Eq. 40, extrapolated endpoint Eq. 42, validity
  bound Eq. 5.

.. _theory-case-sood-coverage:

Sood case coverage (Atalay)
--------------------------------------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   The Atalay-anchored case catalogue lives in
   :mod:`orpheus.derivations.continuous.sood_registry.atalay1997`.

   Brief: Atalay 1997 is the **primary source** for the
   reflected+linearly-anisotropic cross-product cases that lie
   outside both the Sood/Forster/Parsons LA-13511 truth set
   (which focuses on bare configurations) and the
   Burkart–Ishiguro–Siewert 1976 F_N reference [BurkartIshiguroSiewert1976]_ (vacuum-only).
   Specifically, Atalay tabulates :math:`(c, R, f_1)` triples for
   :math:`R \in \{0, 0.25, 0.50, 0.75, 0.99\}` and
   :math:`f_1 \in \{0, 0.10, 0.20, 0.30\}` (slab, even modes, Tables 2-5)
   and :math:`f_1 = 0.10` only (sphere, odd modes, Table 10).

.. _theory-se-errata:

Errata caught (V&V publication artifacts)
================================================================================

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

.. [Case1960] Case, K. M. (1960). "Elementary Solutions of the
   Transport Equation and Their Applications." *Annals of Physics*
   **9**, 1-23.

.. [Mitsis1963] Mitsis, G. F. (1963). "Transport Solutions to the
   Monoenergetic Critical Problems." Argonne National Laboratory
   report ANL-6787.

.. [WestfallMetcalf1973] Westfall, R. M. & Metcalf, D. R. (1973).
   "Singular Eigenfunction Solution of the Monoenergetic Neutron
   Transport Equation for Finite Radially Reflected Critical
   Cylinders." *Nuclear Science and Engineering* **52**, 1-11.
   DOI 10.13182/NSE73-A23285.

.. [Atalay1997] Atalay, M.A. (1997).
   "The reflected slab and sphere criticality problem with anisotropic
   scattering in one-speed neutron transport theory."
   *Progress in Nuclear Energy* **31**\ (3), 229-252.
   DOI: 10.1016/0149-1970(95)00094-1.

.. [McCormickKuscer1965] McCormick, N.J., Kušcer, I. (1965).
   "Bi-orthogonality relations for solving half-space transport problems."
   *J. Math. Phys.* **6**, 1939.

.. [BurkartIshiguroSiewert1976] Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976).
   "Neutron transport in two dissimilar media with anisotropic scattering."
   *Nuclear Science and Engineering* **61**, 72-81.

* **Metcalf & Zweifel 1968** — *Nucl. Sci. Eng.* **33**, 318. — the
  singular-subtraction technique used in WM-72 Eqs. 31 and 33.
* **Atkinson 1976** — SIAM, *A Survey of Numerical Methods for the
  Solution of Fredholm Integral Equations of the Second Kind*,
  Chapter 6 (product integration).
* **Berrut & Trefethen 2004** — *SIAM Review* **46**, 501-517,
  "Barycentric Lagrange Interpolation." — the differentiation
  matrix used for the Eq. 31 diagonal handling.
