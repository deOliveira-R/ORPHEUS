.. _theory-case-method:

================================================================================
Case Singular-Eigenfunction Method (Atalay 1997 family)
================================================================================

.. note::

   **Stub-grade page (Wave 2-B initial drop).** The archivist agent
   will expand each TODO marker into the full narrative + numerical
   evidence + literature backstory.

   Scope of this page:

   * **Reflected slab + linearly anisotropic scattering** — Atalay
     1997 [Atalay1997]_ Eq. 46 (even modes).
   * **Reflected sphere + linearly anisotropic scattering** — Atalay
     1997 Eq. 54 (odd modes; parity-flip of slab Eq. 46).
   * Foundation: Case singular-eigenfunction expansion + first-order
     Fredholm iteration + four NEW parallel half-range relations
     (Atalay Eqs. 28-31) that close the deficit blocking previous
     reflected-slab attempts.

   The Branch-1 SymPy verifications (V_case.1 through V_case.9)
   constitute the **fully verified algebra-of-record** for the Atalay
   1997 reduction chain. Branch-2 production solvers reproduce
   Atalay Tables 2-3 and Sood ``Ua-1-0-SP`` at ~1-7% absolute on the
   critical thickness — looser than the F_N method's 1e-5 anchor
   (driven by the Eq 42 z_0 form's ~1.5% gap) but sufficient for
   the **Atalay-unique cross-product cases** (reflected slab +
   linearly anisotropic, sphere with f_1 = 0.10) that no other
   ORPHEUS solver covers.

Why a separate package from the F_N method
================================================================================

.. _theory-case-vs-fn-method:

.. note:: TODO — Archivist expansion needed.

   The SymPy module lives at
   :mod:`orpheus.derivations.continuous.case_method.origins.derivations`.

   Brief: the F_N method (Siewert-Benoist 1979 + Grandjean-Siewert
   1979 for slab; Siewert-Thomas 1986 for sphere) and the Case
   singular-eigenfunction method (Atalay 1997) belong to the same
   broad "Case 1960 / Mitsis 1963 / McCormick-Kušcer 1965" family,
   but the structural mathematics is fundamentally different. The
   F_N method imposes the exact half-space exit-distribution equation
   at :math:`(N+1)` collocation points and solves a determinantal
   equation; Atalay's method uses the **full normal-mode expansion**
   (discrete + continuum) and reduces the boundary-condition closure
   to a Fredholm integral equation iterated to first order.

   The two methods access the same Wiener-Hopf X-function but
   exercise it differently: F_N uses :math:`X(\xi_\beta)` at
   collocation points :math:`\xi_\beta`; Atalay uses
   :math:`X(\pm \nu_0)` and :math:`X(-\nu)` for :math:`\nu \in (0, 1)`
   integrated against half-range moments. The two are
   **structurally independent** above the trusted-library line
   (sharing only ``numpy``, ``scipy.integrate.quad``, ``mpmath``).

The four new parallel half-range relations
================================================================================

.. _theory-case-half-range:

.. note:: TODO — Archivist expansion needed.

   Brief: Atalay's load-bearing technical contribution is the
   four NEW half-range relations Eqs. 28-31. The McCormick-Kušcer
   1965 bi-orthogonality relations (Atalay Eqs. 18-21) integrate the
   weight :math:`[\phi_{0+}(\mu) + B c \nu_0/2] \gamma(\mu)
   (\nu_0 - \mu)` against the four half-range basis members. These
   suffice for the half-space Milne problem but **not** for the
   reflected-slab boundary condition Atalay Eq. 16, which requires
   a second weight :math:`[\phi_\nu(\mu) + c \nu/(2 \bar\nu)]
   \gamma(\mu)`. Eqs. 28-31 close the deficit by integrating this
   parallel weight against the same four basis members.

   The structural parallelism is verified in
   :func:`...origins.derivations.derive_atalay_half_range_eqs28_to_31`
   (V_case.3): Eqs. 28-31 share the Wiener-Hopf X-function
   factor structure (:math:`X(\pm\nu_0)`, :math:`X(\pm\nu')`,
   :math:`N(\nu)`) of the McCormick-Kušcer set.

Slab criticality (Eq 46)
================================================================================

.. _theory-case-slab-eq46:

.. math::
   :label: case-method-eq46

   \pm\frac{\pi}{2} - \arctan\!\frac{R \sin[(d-z_0)/|\nu_0|] + \sin[(d+z_0)/|\nu_0|]}
                                       {R \cos[(d-z_0)/|\nu_0|] - \cos[(d+z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(K_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+K_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.case_method.origins.derivations`
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

Sphere criticality (Eq 54) via parity flip
================================================================================

.. _theory-case-sphere-eq54:

.. math::
   :label: case-method-eq54

   \arctan\!\frac{\sin[(d+z_0)/|\nu_0|] - R \sin[(d-z_0)/|\nu_0|]}
                {\cos[(d+z_0)/|\nu_0|] + R \cos[(d-z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(L_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+L_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.case_method.origins.derivations`
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

X-function and extrapolated endpoint
================================================================================

.. _theory-case-x-function-eq40:

.. math::
   :label: case-method-eq40

   X(\mu) = \exp\!\Bigg\{ -\frac{c}{2} \int_0^1 d\nu\, g_1(c,\nu)\,
       \Big[d^2(\nu^2)\Big(1 + \frac{c\nu^2}{1-\nu^2}\Big)
            + 3 f_1 (1-c)^2 \nu^2 d(-\nu^2)\Big]\,\ln(\nu - \mu) \Bigg\}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.case_method.origins.derivations`
   (function ``derive_atalay_x_function_eq40``).
   Test gate:
   :func:`tests.derivations.test_case_method_symbolic.test_v_case_8_x_function_eq40`.

.. _theory-case-extrapolated-endpoint-eq42:

.. math::
   :label: case-method-eq42

   z_0 = -\frac{\nu_0}{2} \ln\!\frac{d(-\nu_0\bar\nu)}{d(\nu_0\bar\nu)}
        + \frac{c\,\nu_0}{4} \int_0^1 d\mu\, g_1(c,\mu)\,
                                \Big[d^2(\mu^2)\Big(1 + \frac{c\mu^2}{1-\mu^2}\Big)
                                     + 3 f_1 (1-c)^2 \mu^2 d(-\mu^2)\Big]\,
                                \ln\!\frac{\nu_0 + \mu}{\nu_0 - \mu}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.case_method.origins.derivations`
   (function ``derive_atalay_extrapolated_endpoint_eq42``).

   Implementation note: reproduces Atalay Table 1 :math:`z_0(c, f_1=0)`
   to ~1.5% absolute. The systematic gap is believed to be a
   Case-Zweifel completeness-sum normalisation discrepancy with the
   published Eq. 42 form (where :math:`(c/2) \int_0^1 g_1[\text{bracket}]
   d\nu \approx 0.99` rather than 1.00 in standard Case-Zweifel
   normalisations). Investigation in the closeout memo.

Validity bound (Atalay Eq 5)
================================================================================

.. _theory-case-validity-eq5:

.. math::
   :label: case-method-eq5

   c \le 1 + \frac{1}{3 f_1}

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.case_method.origins.derivations`
   (function ``derive_atalay_validity_bound_eq5``).
   Test gate:
   :func:`tests.derivations.test_case_method_symbolic.test_v_case_9_validity_bound_eq5`.

   Brief: Atalay Eq. 5 is the **one-pair-of-discrete-modes** range.
   For :math:`f_1 = 0` the bound is trivial (all c); for
   :math:`f_1 = 0.30` the upper limit is :math:`c \le 19/9 \approx 2.111`.
   Outside this band complex eigenvalues appear (Dahl-Sjöstrand
   1979; Kohut 1993) and Atalay's first-order Fredholm iteration
   breaks. The slab/sphere solvers raise :class:`ValueError` if
   :math:`c` violates the bound.

Sood case coverage
================================================================================

.. _theory-case-sood-coverage:

.. note:: TODO — Archivist expansion needed.

   The Atalay-anchored case catalogue lives in
   :mod:`orpheus.derivations.continuous.sood_registry.atalay1997`.

   Brief: Atalay 1997 is the **primary source** for the
   reflected+linearly-anisotropic cross-product cases that lie
   outside both the Sood/Forster/Parsons LA-13511 truth set
   (which focuses on bare configurations) and the
   Burkart-Ishiguro-Siewert 1976 F_N reference (vacuum-only).
   Specifically, Atalay tabulates :math:`(c, R, f_1)` triples for
   :math:`R \in \{0, 0.25, 0.50, 0.75, 0.99\}` and
   :math:`f_1 \in \{0, 0.10, 0.20, 0.30\}` (slab, even modes, Tables 2-5)
   and :math:`f_1 = 0.10` only (sphere, odd modes, Table 10).

References
================================================================================

.. [Atalay1997] Atalay, M.A. (1997).
   "The reflected slab and sphere criticality problem with anisotropic
   scattering in one-speed neutron transport theory."
   *Progress in Nuclear Energy* **31**(3), 229-252.
   DOI: 10.1016/0149-1970(95)00094-1.

.. [McCormickKuscer1965] McCormick, N.J., Kušcer, I. (1965).
   "Bi-orthogonality relations for solving half-space transport problems."
   *J. Math. Phys.* **6**, 1939.

.. [BurkartIshiguroSiewert1976] Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976).
   "Neutron transport in two dissimilar media with anisotropic scattering."
   *Nuclear Science and Engineering* **61**, 72-81.
