.. _theory-galerkin-spectral:

Galerkin Spectral Method (slab + sphere, linearly anisotropic)
===============================================================

The :mod:`orpheus.derivations.continuous.galerkin_spectral` package
implements a **Galerkin spectral expansion of Carlvik's integral
equation** [Carlvik1968]_ for the criticality of a multiplying
slab/sphere with linearly anisotropic scattering, following Dahl &
Sjöstrand (1979) [DahlSjostrand1979]_. This is the third
structurally-independent Pillar-2 verification method in the
project, alongside :ref:`F_N collocation <theory-fn-method>` and
:ref:`singular eigenfunction Fredholm <theory-singular-eigenfunction>`.

.. _theory-galerkin-spectral-carlvik-primitives:

Carlvik recurrences as named primitives (folder-naming rationale)
------------------------------------------------------------------

The folder is named ``galerkin_spectral/`` (method-canonical) rather
than ``carlvik_galerkin/`` (its previous half-author/half-method name),
but the integral-form recurrences introduced by Carlvik (1968)
[Carlvik1968]_ remain first-class objects inside the package as
:mod:`...galerkin_spectral.core.carlvik_recurrences`. This deserves
explicit explanation:

* **Author-as-mathematical-primitive-name is fine.** "Carlvik
  recurrences" is a *named mathematical object* in the transport
  literature in the same sense that "Bickley–Naylor functions",
  "Wigner 3-j symbol", or "Gegenbauer polynomials" are — established
  shorthand for a specific construction that subsequent papers
  reference verbatim. Carlvik (1968) introduced the integral-form
  recurrences for the matrix elements
  :math:`\int_0^1 P_n(\mu)\, E_k(d/\mu)\,d\mu`; those recurrences are
  what the ``carlvik_recurrences.py`` filename advertises.

* **Author-as-method-name is not fine.** "Carlvik method" or
  "Carlvik–Galerkin method" violates the project's no-author-folders
  rule because the *method* on top is **Galerkin spectral expansion**
  in Legendre polynomials — the discretisation idea that converts the
  integral equation into a finite eigenvalue problem. Carlvik
  contributed *recurrences*, Dahl–Sjöstrand (1979) contributed the
  *Galerkin matrix construction* on top. The folder name reflects the
  structural method that owns the file layout; the recurrence library
  retains its author-shorthand inside ``core/`` because that is how
  the literature names it.

The same discipline applies elsewhere in ORPHEUS: ``fn_method/`` is a
method (the F_N family), ``peierls_nystrom/`` names the Peierls
integral form + Nyström quadrature (both textbook method names, not
authors), but the inline ``bickley_naylor.py`` would be allowed because
"Bickley–Naylor functions" name a primitive. The boundary is
**method ↔ folder, primitive ↔ filename**.

The Sood/Forster/Parsons LA-13511 (1999) [SoodLA13511_1999]_
``*-1-1-SL/SP`` benchmark cases (anisotropic slab + sphere) are
verified by this package against the Dahl-Sjöstrand 1979 Tables I
and II tabulations.

Distinguishing feature: the Eq. (4) block-matrix linearization yields
the **full eigenvalue spectrum** (typically 6-11 eigenvalues per
configuration, including complex-conjugate pairs at high anisotropy),
unique among the project's verification pillars.

.. _carlvik-galerkin-V_cg-1:

V_cg.1 — Galerkin LHS = 2 F_m via Legendre orthogonality
---------------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   The SymPy derivation lives in
   :mod:`orpheus.derivations.continuous.galerkin_spectral.origins.derivations`
   (function ``derive_galerkin_lhs_identity``). Test gate:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_1_galerkin_lhs_identity`.

   Brief: Multiplying Eq. (1) by :math:`P_m(x)` and integrating over
   :math:`x \in [-1, +1]`, Legendre orthogonality
   :math:`\int_{-1}^{+1} P_m P_n\, dx = 2/(2m+1) \delta_{mn}`
   collapses the LHS to :math:`2 F_m`.

.. _carlvik-galerkin-V_cg-2:

V_cg.2 — Eq. (3) matrix eigenvalue structure
---------------------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_eq3_matrix_eigenvalue``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_2_eq3_matrix_eigenvalue_form`.

   Brief: Galerkin projection of Eq. (1) yields :math:`(A - 3\bar\mu(c-1)
   B) F = (1/(cd)) F`. The structural decomposition into a
   :math:`c`-independent part (:math:`A`) and a :math:`c`-dependent
   correction (:math:`-3\bar\mu(c-1) B`) is the load-bearing form
   that enables the linearization.

.. _carlvik-galerkin-V_cg-3:

V_cg.3 — A_{0,0}(a) closed form (slab fundamental)
---------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_low_order_A_mn_slab``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_3_A_00_closed_form_slab`.

   Brief: SymPy directly integrates the defining double integral
   :math:`A_{0,0}(a) = \tfrac{1}{2}\int\!\!\int E_1(a|x-y|)\,dy\,dx`,
   yielding :math:`A_{0,0}(a) = 2 E_1(2a) + 2/a - e^{-2a}/a -
   1/(2a^2) + e^{-2a}/(2a^2)`. Independently verified by hand
   derivation using the universal :math:`\int E_n = -E_{n+1}` recursion.

.. _carlvik-galerkin-V_cg-4:

V_cg.4 — Slab/sphere basis parity
----------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_low_order_A_mn_sphere``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_4_basis_parity_slab_vs_sphere`.

   Brief: Slab fundamental is even in :math:`x` → uses
   :math:`P_0, P_2, P_4, \ldots`. Sphere fundamental is odd in
   :math:`r` (since :math:`r\phi(r)` is the reduced flux variable)
   → uses :math:`P_1, P_3, P_5, \ldots`.

.. _carlvik-galerkin-V_cg-5:

V_cg.5 — B_{m,n} boundary-chord rank-1 structure
-------------------------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_B_mn_boundary_chord_form``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_5_B_mn_boundary_chord_rank_one`.

   Brief: The boundary-chord term in :math:`B_{m,n}` decomposes as
   a product of two single integrals,
   :math:`(\int P_m(x) K_q(x; a) dx) \cdot (\int P_n(y) y^q dy)`.
   The second factor vanishes by Legendre orthogonality except for
   :math:`n = q` — so the rank-1 contribution affects only one
   column of :math:`B`. This structural simplification reduces the
   :math:`O(N^2)` 2-D integration cost.

.. _carlvik-galerkin-V_cg-6:

V_cg.6 — Eq. (4) block-matrix linearization
--------------------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_eq4_block_linearization``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_6_eq4_block_linearization`.

   Brief: With :math:`G = d(A + 3\bar\mu B)`, :math:`H = -3\bar\mu d B`,
   :math:`K = c F`, the block system
   :math:`\binom{G\ H}{I\ 0}\binom{F}{K} = (1/c)\binom{F}{K}`
   is algebraically equivalent to Eq. (3). The block matrix is
   :math:`2N \times 2N` and admits ALL :math:`2N` eigenvalues
   :math:`1/c_j` in a single standard non-symmetric eigenproblem
   call — replacing Carlvik's iterative :math:`c`-search.

.. _carlvik-galerkin-V_cg-7:

V_cg.7 — Carlvik 1968 Eq. (4b) sign correction
-----------------------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_carlvik_eq4b_corrected_form``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_7_carlvik_eq4b_correction_documented`.

   Brief: Dahl-Sjöstrand 1979 explicitly flags Carlvik 1968 NSE 31,
   295 Eq. (4b) as misprinted: "Note that the sign of the last term
   in Carlvik's expression (4b) has been misprinted." The
   Branch-2 production code in
   :mod:`...galerkin_spectral.core.carlvik_recurrences` does NOT
   transcribe Carlvik 1968 Eq. (4b); it computes :math:`B_{m,n}`
   from the defining double integral of Eq. (1), bypassing the
   typo entirely.

.. _carlvik-galerkin-V_cg-8:

V_cg.8 — μ̄ = 0 isotropic limit
--------------------------------

.. note:: TODO — Archivist expansion needed.

   ``derive_isotropic_limit``. Test:
   :func:`tests.derivations.test_galerkin_spectral_symbolic.test_v_cg_8_isotropic_limit`.

   Brief: At :math:`\bar\mu = 0`, Eq. (3) reduces to :math:`A F =
   (1/(cd)) F` — Carlvik 1968's original isotropic eigenvalue
   equation. Cross-verified at the numerical level against
   :ref:`F_N method <theory-fn-method>` to ≤ 1e-5 absolute on
   :math:`c_{\rm crit}` across the full range
   :math:`c \in (1, 1.5]` for both slab and sphere.

References
----------

.. [DahlSjostrand1979] Dahl, E. B. & Sjöstrand, N. G. (1979).
   "Eigenvalue spectrum of multiplying slabs and spheres for
   monoenergetic neutrons with anisotropic scattering."
   *Nuclear Science and Engineering* **69**, 114-125.
   DOI: 10.13182/NSE69-114.

.. [Carlvik1968] Carlvik, I. (1968).
   *Nuclear Science and Engineering* **31**, 295-300.
   (Original derivation of :math:`A_{m,n}` and :math:`B_{m,n}`
   recurrences. Eq. (4b) printed sign typo corrected by
   Dahl-Sjöstrand 1979.)

.. [SoodLA13511_1999] Sood, A., Forster, R.A., Parsons, D.K. (1999).
   "Analytical Benchmark Test Set for Criticality Code Verification."
   Los Alamos National Laboratory report LA-13511.
