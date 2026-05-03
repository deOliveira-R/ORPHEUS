.. _theory-singular-eigenfunction:

================================================================================
Singular Eigenfunction Expansion (Mitsis-Westfall-Metcalf family)
================================================================================

.. note::

   **Stub-grade page (Phase B1).** The archivist agent will expand
   each TODO marker into the full narrative + numerical evidence +
   literature backstory.

   Scope of this page:

   * **Bare-critical infinite cylinder, 1G isotropic scattering**
     via Westfall-Metcalf 1973 (Case ``Ua-1-0-CY`` of the Sood
     LA-13511 benchmark family).

   The Branch-1 SymPy verification (V_se-cyl.1 through V_se-cyl.7)
   is the **fully verified algebra-of-record** for the WM-72
   derivation, and includes the discovery that printed Eq. 17 has
   a typo. The Branch-2 production solver is a **prototype**
   limited to ~1% relative accuracy at the bare-cylinder critical
   radius; tightening to the Sood 1e-5 target requires a future
   hardening pass (graded-mesh refinement near the kernel diagonal,
   or full Mitsis-WM Fredholm iteration).

Why a separate package from the F_N method
================================================================================

.. _theory-se-vs-fn:

.. note:: TODO — Archivist expansion needed.
   The SymPy module lives at
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations`.
   Brief: both the F_N method (slab/sphere) and the
   singular-eigenfunction expansion (cylinder) belong to the
   broad "Case 1960 / Mitsis 1963" family of methods, but the
   structural mathematics is fundamentally different. The F_N
   method (Siewert 1979-1986) builds a finite (N+1)×(N+1)
   determinantal equation from Wiener-Hopf X-functions; this
   approach **does not converge** for the bare cylinder
   (Westfall-Metcalf 1973 explicitly notes: "we found that the
   solution as formulated by Mitsis is not convergent" for the
   bare cylinder). WM-72's reformulation uses Bessel-kernel
   expansion + iterative Fredholm scheme — different mathematical
   machinery. The two packages share **only** the dispersion-root
   primitive (``case_nu0``), which is a medium property
   independent of geometry; everything above the trusted-library
   line is structurally distinct, satisfying the
   ``algebra-of-record`` independence requirement.

Bare-critical cylinder via WM-72
================================================================================

.. _theory-se-wm72-derivation:

.. note:: TODO — Archivist expansion needed.
   The complete WM-72 derivation chain is captured in 7 SymPy
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
     bare-cylinder reduction :math:`c_1 = c_2 \implies D(\nu) = 0`
     (and analogous algebraic cancellations of :math:`d_0` and
     :math:`B(\nu)`).
   * V_se-cyl.5 (:func:`...origins.cylinder_derivations.derive_bare_cylinder_criticality_condition`):
     bare-cylinder criticality reduces to a single integral equation
     after dispersion-root substitution.
   * V_se-cyl.6 (:func:`...origins.cylinder_derivations.derive_discrete_eigenfunction_normalization`):
     discrete normalisation matches WM-72 Eq. 21d (provides confidence
     in the corrected :math:`\eta_0` form via cross-check against
     a separate published formula).
   * V_se-cyl.7 (:func:`...origins.cylinder_derivations.derive_flux_reconstruction_bare_cylinder`):
     bare-cylinder neutron density profile is :math:`\rho(r) = J_0(r/u_0)`
     (the dominant Case eigenfunction with imaginary
     :math:`\nu_0 = i u_0` for :math:`c > 1`).

   Test gates: :mod:`tests.derivations.test_singular_eigenfunction_cylinder`
   (one foundation test per ``derive_*()`` plus production-solver
   sanity tests + L1 Sood ``Ua-1-0-CY`` reference-value gate).

Production solver and accuracy ceiling
================================================================================

.. _theory-se-wm72-numerics:

.. note:: TODO — Archivist expansion needed.
   The Branch-2 solver is at
   :mod:`orpheus.derivations.continuous.singular_eigenfunction.cylinder.one_group`.

   Brief: the production solver discretises the Mitsis cylindrical
   integral transport equation (WM-72 Eq. 6a) directly via
   Gauss-Legendre quadrature on :math:`(0, R)`. The cylinder kernel

   .. math::

      K(r, t) = \int_0^1 \frac{K_0(\max(r,t)/\mu)\,I_0(\min(r,t)/\mu)}
                              {\mu^2}\,d\mu

   has a logarithmic singularity at :math:`r = t` (related to the
   2D Green's function of the streaming operator). Single-cell
   product integration handles the diagonal log-singularity in this
   prototype, achieving :math:`\mathcal{O}(1/n)` algebraic
   convergence. At :math:`n = 128` GL nodes, accuracy is :math:`\sim 0.5\%`
   relative on Sood ``Ua-1-0-CY``.

   **Tightening to the 1e-5 Sood target** requires either:

   1. **Graded mesh refinement** (Atkinson 1976 product integration)
      that places extra nodes near the kernel diagonal where the
      log singularity sits.
   2. **Full Mitsis-WM Fredholm iteration** (WM-72 Eqs 28-33) that
      isolates the discrete Case eigenfunction and iterates the
      continuum amplitude :math:`A'(\nu)` separately. This converges
      faster because the discrete eigenfunction is treated
      analytically and the continuum integration is regularised
      by the principal-value-plus-delta structure of
      :math:`\eta_\nu(\mu)`.

   The Variant α cylinder cross-check at 8.5e-6 (case ``Ua-1-O-CY`` —
   see :mod:`tests.derivations.test_peierls_greens_function_cylinder_xverif_sood2003`)
   already provides a structurally-independent external reference at
   the target accuracy; the WM-72 prototype here is the **second**
   structurally-independent anchor (different mathematical pillar
   than Variant α and the Bickley-Naylor :math:`\mathrm{Ki}_n`
   integrals of the cylinder Nyström solver).

Cross-check vs Variant α at Sood ``Ua-1-0-CY``
================================================================================

.. _theory-se-wm72-xverif:

.. note:: TODO — Archivist expansion needed.
   Test gate:
   :mod:`tests.derivations.test_singular_eigenfunction_cylinder_xverif`.

   Brief: the WM-72 prototype agrees with the Variant α cylinder
   solver at the Sood ``Ua-1-0-CY`` configuration to :math:`\le 2\%`
   relative — the WM-72 accuracy floor. Both solvers reproduce the
   published Sood :math:`r_c = 1.72500292` mfp; Variant α at 8.5e-6,
   WM-72 at ~1%. The agreement at the WM-72 floor is the L1 evidence
   for structural independence: a third method, Bickley-Naylor /
   :math:`\mathrm{Ki}_n` integrals (cylinder Nyström / peierls_nystrom)
   shares NO machinery with either Variant α or WM-72, so a future
   triangle test (WM-72 ↔ Variant α ↔ Nyström) will close the V&V
   audit cleanly when the WM-72 floor is tightened.

Errata caught (V&V publication artifacts)
================================================================================

.. _theory-se-errata:

.. note:: TODO — Archivist expansion needed.

   Brief: this implementation slice caught one typo in the published
   literature, mirroring the discipline-of-record practice from the
   F_N first slice (Sood Eq. 28 typo).

   * **WM-72 Eq. 17** (printed): :math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}/
     (\nu_{0l}^2 - \mu^2)`. **Direct substitution into Eq. 15 fails**:
     the LHS reduces to :math:`c\nu_0\mu^2`, the RHS reads
     :math:`c\nu_0^2\mu^2` — mismatched power of :math:`\nu_0`. The
     correct form is :math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}^2/
     (\nu_{0l}^2 - \mu^2)`, which:

     - Satisfies Eq. 15 exactly.
     - Reproduces Eq. 21d's :math:`\nu_0^4` factor on
       :math:`N_0 = \int_0^1 \mu^2\eta_0^2\,d\mu` (the printed
       single-:math:`\nu_0` form would give :math:`\nu_0^2` instead).
     - Closes the half-range normalisation Eq. 14 under dispersion.

   The Branch-2 production code uses the corrected form. SymPy
   re-derivation in V_se-cyl.2 caught the typo on first run.

References
==========

* **Westfall & Metcalf 1973** — *Nucl. Sci. Eng.* **52**, 1-11.
  "Singular Eigenfunction Solution of the Monoenergetic Neutron
  Transport Equation for Finite Radially Reflected Critical
  Cylinders." DOI 10.13182/NSE73-A23285.
* **Mitsis 1963** — Argonne National Laboratory report ANL-6787,
  "Transport Solutions to the Monoenergetic Critical Problems."
* **Atkinson 1976** — SIAM, *A Survey of Numerical Methods for the
  Solution of Fredholm Integral Equations of the Second Kind*,
  Chapter 6 (product integration).
