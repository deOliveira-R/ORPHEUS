.. _theory-galerkin-spectral:

==========================================================================
Galerkin Spectral Method (slab + sphere, linearly anisotropic)
==========================================================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this before modifying the Galerkin spectral reference solver.**

- **What this implements**: Galerkin spectral expansion of Carlvik's
  integral equation [Carlvik1968]_ for the criticality of a
  multiplying slab/sphere with linearly anisotropic scattering,
  following Dahl & Sjöstrand (1979) [DahlSjostrand1979]_.
- **Position in the V&V stack**: third structurally-independent
  Pillar-2 verification method, alongside
  :ref:`F_N collocation <theory-fn-method>` and
  :ref:`singular eigenfunction Fredholm <theory-singular-eigenfunction>`.
  Pillar classification: closed-form (matrix eigenvalue).
- **Folder-naming rationale**: ``galerkin_spectral/`` is
  method-canonical (Legendre-Galerkin spectral expansion). The
  Carlvik recurrences live inside as named-mathematical-primitive
  filenames (``core/carlvik_recurrences.py``) — that's the literature
  convention. See :ref:`theory-galerkin-spectral-carlvik-primitives`
  below for the full rule.
- **Cross-references**: Sood case registry consumes this method's
  truth values via :ref:`theory-sood-registry`; F_N method serves as
  the canonical structurally-independent cross-check via the
  ``test_carlvik_galerkin_xverif_fn`` test family.

The page below covers the Carlvik-recurrence primitive design
discipline, the Galerkin assembly, and the V&V pillar position.

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

.. _galerkin-spectral-basis-space-and-variational-principle:

The Galerkin spectral basis space and the variational principle
================================================================

This section is the rich-narrative companion to
:class:`~orpheus.derivations.continuous.galerkin_spectral.BasisSpace`
— the math-heart class for the Galerkin spectral pillar. It develops
the Galerkin orthogonality principle, the Carlvik integral form, the
Dahl–Sjöstrand block-matrix linearisation, and the eigenvalue
spectrum theory at graduate-textbook depth. The class docstring
covers the same material in compressed form for readers who land on
the class definition first; this section is the reference for
everyone else.

The angular flux as an element of L²([0, R] × [-1, +1])
--------------------------------------------------------

The starting object is the angular neutron flux
:math:`\psi(r, \mu)`, defined on the spatial domain
:math:`[0, R]` (sphere) or :math:`[-a, +a]` (slab) and the angular
domain :math:`[-1, +1]` (the polar cosine). For one-speed transport
in a homogeneous, possibly anisotropically-scattering medium with
linearly-anisotropic differential scattering kernel
:math:`\Sigma_s(\mu' \to \mu) = \tfrac{1}{4\pi}[\Sigma_{s,0} + 3
\Sigma_{s,1}\,\mu'\mu]` and isotropic fission, the steady-state
Boltzmann transport equation is

.. math::
   :label: galerkin-spectral-bte

   \mu \frac{\partial\psi}{\partial r} + \Sigma_t \psi(r, \mu)
       = \tfrac{c\Sigma_t}{2}\big[\phi(r) + 3\bar\mu\,J(r)\,\mu\big]

where :math:`\phi(r) = \int_{-1}^{+1} \psi(r, \mu)\, d\mu` is the
scalar flux, :math:`J(r) = \int_{-1}^{+1} \mu\,\psi(r, \mu)\, d\mu`
is the current, :math:`c = (\Sigma_s + \nu\Sigma_f) / \Sigma_t` is
the mean number of secondaries per collision, and
:math:`\bar\mu = \Sigma_{s,1} / \Sigma_{s,0}` is the linearly-
anisotropic mean cosine.

The natural function space is
:math:`L^2([0, R] \times [-1, +1])` — Hilbert-space of square-
integrable :math:`(r, \mu)` functions with the inner product

.. math::

   \langle f, g\rangle = \int_0^R\!\!\int_{-1}^{+1}
       f(r, \mu)\, g(r, \mu)\, d\mu\, dr.

Solutions of (:eq:`galerkin-spectral-bte`) live in this Hilbert
space, and Galerkin spectral methods build approximations as
finite-dimensional subspace projections.

Choice of basis: full-range Legendre polynomials
-------------------------------------------------

Galerkin spectral picks the **Legendre polynomial family
:math:`\{P_n(x)\}`** as the trial-and-test basis on the
rescaled spatial variable :math:`x \in [-1, +1]` (slab) or
:math:`x = r/R` after the reduced-flux substitution :math:`u(r) =
r\,\phi(r)` (sphere). The choice is motivated by three properties
that simultaneously make the scheme:

* **Orthogonal**: :math:`\int_{-1}^{+1} P_m P_n\, dx =
  \tfrac{2}{2n+1}\,\delta_{mn}` — the Galerkin matrix's
  diagonal-only mass matrix is the *inverse* of the
  orthonormalising weights, and the Galerkin LHS reduces (via
  V_cg.1, see :ref:`carlvik-galerkin-V_cg-1`) to a clean
  :math:`2 F_m`.
* **Parity-stratified**: :math:`P_n` is even in :math:`x` for
  even :math:`n` and odd for odd :math:`n`. Slab scalar-flux
  parity :math:`\phi(-x) = \phi(x)` selects only **even-Legendre**
  basis functions :math:`P_0, P_2, \ldots, P_{2(N-1)}`. The
  reduced flux :math:`u(r) = r\,\phi(r)` is antisymmetric in
  :math:`r \to -r`, selecting only **odd-Legendre** basis
  functions :math:`P_1, P_3, \ldots, P_{2N-1}` for the sphere.
  The parity restriction halves the basis size and aligns with
  the geometry's natural parity sector.
* **Kernel-friendly**: the angular-pre-integrated transport
  kernel :math:`E_1(\Sigma_t |x - x'|)` decomposes into the
  Legendre-moment matrix elements :math:`A_{m,n}` and
  :math:`B_{m,n}` (Carlvik 1968) as polynomial-in-:math:`x \otimes`
  log-singular-in-:math:`|x-x'|` integrals, which the recurrences
  evaluate stably.

Why not Chebyshev, Hermite, Jacobi? Chebyshev is the natural choice
when the singularity structure pulls toward the endpoints
(:math:`\sqrt{1-x^2}` weight); the bare-critical transport problem's
endpoint behaviour is logarithmic, not algebraic, so the Legendre
weight (uniform) gives the same convergence rate without the
endpoint pull. Hermite needs the unbounded interval; the bare slab
is bounded. Jacobi with :math:`(1, 1)` weight is equivalent to
shifted-Legendre on :math:`[0, 1]` and would serve the half-range
F_N method (see ``MomentSpace``); the full-range version uses
plain Legendre.

The Galerkin orthogonality principle
-------------------------------------

Galerkin's variational principle (Galerkin 1915) is a
*finite-dimensional best-approximation* statement: among all
Legendre coefficient vectors :math:`(F_0, F_2, \ldots, F_{2(N-1)})`
of length :math:`N`, the one that solves the **Galerkin equations**

.. math::
   :label: galerkin-spectral-orthogonality

   \langle R[\phi_N], P_m\rangle = 0
       \quad\text{for all}\quad m \in \{0, 2, \ldots, 2(N-1)\}

(where :math:`R[\phi_N] = \mathcal L \phi_N - f` is the residual
of the truncated trial function under the operator
:math:`\mathcal L`) is the projection of the true solution onto
the trial subspace **in the energy norm induced by**
:math:`\mathcal L`.

The orthogonality condition (:eq:`galerkin-spectral-orthogonality`)
says: *the residual is orthogonal to the trial space*. Equivalently,
the truncated solution :math:`\phi_N` is the closest element of the
trial space to the true solution :math:`\phi`, measured in the
energy norm :math:`\|\cdot\|_{\mathcal L}` (Strang–Fix 1973
§1.6, Céa's lemma).

Three consequences flow immediately:

1. **The truncation error decays in proportion to
   best-approximation error.** Céa's lemma states
   :math:`\|\phi - \phi_N\|_{\mathcal L} \le C\,
   \inf_{v \in V_N}\|\phi - v\|_{\mathcal L}` where :math:`C`
   depends only on the operator's continuity / coercivity
   constants and :math:`V_N` is the trial subspace. For
   Legendre polynomials with the unit weight, the best
   approximation of an analytic function decays
   *super-algebraically* — faster than any power of :math:`1/N`
   — and Dahl–Sjöstrand observe 7-significant-figure agreement
   already at :math:`N = 9`.

2. **The basis CAN be non-orthogonal — the only requirement is
   linear independence.** The orthogonality of Legendre polynomials
   is not what makes the Galerkin scheme work; it makes the
   *implementation* clean (diagonal mass matrix, cheap
   change-of-basis). A Galerkin scheme on a non-orthogonal basis
   works equally well in principle but with a denser stiffness
   matrix.

3. **The orthogonality condition is symmetric in trial and test
   spaces.** If trial and test bases coincide (true Galerkin),
   the resulting matrix inherits a symmetry property when the
   underlying operator is self-adjoint. The Carlvik integral
   form mixes :math:`E_1` (symmetric kernel) and the rank-1
   boundary-chord term (which is *not* :math:`L^2`-self-adjoint
   under the unit weight), so the matrix is non-symmetric. But
   the *projection* still gives best-approximation in the
   Galerkin norm — this is what guarantees the convergence rate.

True Galerkin vs Petrov–Galerkin vs collocation vs least-squares
-----------------------------------------------------------------

The choice trial = test (true Galerkin) is what distinguishes
:class:`BasisSpace` from neighbouring projection methods used in
the same problem family:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - Method
     - Trial space
     - Test space
     - Implemented by
   * - **Galerkin** (BasisSpace)
     - :math:`P_n(x)`
     - :math:`P_n(x)` — *same*
     - :class:`BasisSpace`
   * - **Petrov–Galerkin**
     - :math:`P_n(x)`
     - :math:`Q_n(x)` — *different*
     - F_N (in full-range view)
   * - **Collocation**
     - :math:`P_n(x)`
     - :math:`\delta(\xi - \xi_\beta)`
     - F_N (Grandjean–Siewert, Siewert–Thomas grids)
   * - **Least-squares**
     - :math:`P_n(x)`
     - :math:`\mathcal L\,P_n(x)`
     - Not implemented in this codebase

Collocation imposes the residual to be exactly zero at :math:`N+1`
discrete points :math:`\xi_\beta`. The trial space and test space
*differ* (test = delta), so collocation is technically
Petrov–Galerkin with delta tests. ``MomentSpace`` (F_N) collocates
on the Grandjean–Siewert grid (slab) or the Siewert–Thomas
Chebyshev-interior grid (sphere). The **half-range** version of the
moment basis (used by F_N) is Petrov–Galerkin even setting aside the
collocation: the test functions are :math:`\mu^\alpha` polynomials
on :math:`[0, 1]`, NOT the same as the trial functions.

Least-squares takes the test space as :math:`\mathcal L\,P_n(x)` —
the residual is projected onto the *operator-applied* basis. The
resulting normal equations are :math:`\mathcal L^* \mathcal L
\phi_N = \mathcal L^* f`, which is symmetric positive-definite by
construction but doubles the operator's condition-number cost.

True Galerkin is the **best-approximation** scheme in the
operator's natural norm: the projected solution is the closest
trial-space element to the true solution. Petrov–Galerkin and
collocation give different projections (best-approximation in a
different norm; exact match at collocation points respectively).
Least-squares is best-approximation in the residual-:math:`L^2`
norm, which is a different functional.

Carlvik's integral form: the Volterra structure
------------------------------------------------

The Boltzmann equation (:eq:`galerkin-spectral-bte`) is a
first-order integro-differential equation in :math:`(r, \mu)`. To
apply Galerkin on the Legendre basis in :math:`x`, Carlvik (1968)
first **integrates out the angular variable** analytically against
the linearly-anisotropic scattering kernel, producing a closed
integral equation in :math:`x` alone:

.. math::
   :label: galerkin-spectral-carlvik-integral

   \phi(x) = \frac{c}{2}\int_{-1}^{+1}
       \big[E_1(a|x-y|) - 3\bar\mu(c-1) E_3(a|x-y|)\big]
       \phi(y)\, dy + B_{\rm bdy}[\phi](x).

The kernel has two parts:

* **Volume term**: :math:`E_1(a|x-y|)` couples the scalar flux
  at :math:`x` to the scalar flux at :math:`y` through the Peierls
  angular average of the streaming attenuation. The kernel is
  *log-singular at :math:`x = y`*: :math:`E_1(\tau) = -\gamma_E
  -\log\tau + R(\tau)` where :math:`R(\tau)` is a smooth
  remainder. The :math:`E_3` term comes from the linearly-
  anisotropic correction.

* **Boundary-chord term**: :math:`B_{\rm bdy}[\phi](x)`
  encodes the contribution from rays that traverse the slab
  end-to-end without scattering. It has a **rank-1 structure** —
  separable as :math:`K_q(x; a)\,M_q[\phi]` where :math:`M_q[\phi]
  = \int_{-1}^{+1} y^q \phi(y)\, dy` is the :math:`q`-th moment
  and :math:`K_q(x; a) = \tfrac{1}{2}[E_3(a|1-x|) + (-1)^q
  E_3(a|1+x|)]` is the boundary-chord kernel. For slab :math:`q
  = 0`; for sphere :math:`q = 1`. The rank-1 structure means
  :math:`B_{\rm bdy}` contributes a single column to the Galerkin
  matrix (V_cg.5) and dramatically simplifies the assembly cost.

The integral form (:eq:`galerkin-spectral-carlvik-integral`) is
*Volterra-type* in the geometry: the kernel is causal in the
characteristic direction (the streaming-attenuation part of
:math:`E_1`) and the boundary-chord term is a finite-rank
correction. Galerkin projection of this Volterra operator onto
:math:`\{P_n(x)\}_n` produces the matrix

.. math::

   2 F_m = c\,d\sum_n \big[A_{m,n}(a) - 3\bar\mu(c-1)B_{m,n}(a)\big] F_n

(Dahl–Sjöstrand Eq. 3, with the Galerkin LHS already reduced via
V_cg.1) where :math:`d = 2a` is the slab thickness in mean free
paths and the matrix elements are

.. math::

   A_{m,n}(a) &= \tfrac{2n+1}{2}
       \int_{-1}^{+1}\!\!\int_{-1}^{+1}
       P_m(x) P_n(y) E_1(a|x-y|)\, dy\, dx, \\
   B_{m,n}(a) &= \tfrac{2n+1}{2}\bigg[
       \int_{-1}^{+1}\!\!\int_{-1}^{+1}
       P_m(x) P_n(y) E_3(a|x-y|)\, dy\, dx \\
       &\quad - \big(\textstyle\int_{-1}^{+1} P_m(x)\,
         K_q(x; a)\, dx\big)
       \big(\textstyle\int_{-1}^{+1} P_n(y) y^q\, dy\big)
       \bigg].

The matrix elements have the *log-singular-on-the-diagonal*
structure inherited from :math:`E_1`. Carlvik's recurrences
(see :ref:`carlvik-galerkin-V_cg-3`) evaluate them stably; the
Branch-2 production code in
:mod:`...galerkin_spectral.core.carlvik_recurrences` uses tensor-
product Gauss-Legendre quadrature on a single square at moderate
:math:`a` and bypasses Carlvik 1968 Eq. (4b) (which has a
documented sign typo, see :ref:`carlvik-galerkin-V_cg-7`).

The block-matrix linearisation: Eq. (4)
----------------------------------------

The matrix Eq. (3) is **nonlinear in the eigenvalue** :math:`c`
because the LHS coefficient :math:`(c - 1)` multiplies :math:`B`.
A direct iterative :math:`c`-search (Carlvik 1968) is feasible but
yields only the dominant eigenvalue. Dahl–Sjöstrand (1979) introduce
the **block-matrix linearisation** that turns the nonlinear-in-
:math:`c` problem into a **standard non-symmetric linear
eigenproblem of double dimension**.

Define the auxiliary variable :math:`K = c F`. Then Eq. (3)
expands as

.. math::

   d\,(A + 3\bar\mu B) F + (-3\bar\mu d B) K
       &= \frac{1}{c}\,F, \\
   F &= \frac{1}{c}\,K.

In block-matrix form,

.. math::
   :label: galerkin-spectral-eq4

   \underbrace{\begin{pmatrix} G & H \\ I & 0 \end{pmatrix}}_{2N \times 2N}
       \begin{pmatrix} F \\ K \end{pmatrix}
       = \frac{1}{c}\,\begin{pmatrix} F \\ K \end{pmatrix}

with :math:`G = d (A + 3\bar\mu B)`, :math:`H = -3\bar\mu d B`,
:math:`I = \mathbf 1_N` the identity, and :math:`0 = \mathbf 0_N`
the zero block.

(:eq:`galerkin-spectral-eq4`) is a **standard non-symmetric matrix
eigenproblem** :math:`M v = \lambda v` with :math:`\lambda = 1/c`
and :math:`v = (F, K)^\top`. Solving it with
:func:`scipy.linalg.eig` returns ALL :math:`2N` eigenvalues
:math:`(1/c_j)` and the corresponding eigenvectors in one call.
Inverting :math:`\lambda_j \to c_j = 1/\lambda_j` and sorting by
:math:`\Re c_j` ascending gives the spectrum: the smallest positive
real :math:`c_j` is the **fundamental** :math:`c_{\rm crit}`, and
the remaining :math:`2N - 1` eigenvalues are higher-mode
secondary-ratio eigenvalues — typically real for low :math:`\bar\mu`,
developing complex-conjugate pairs at high :math:`\bar\mu` for
higher modes (Dahl–Sjöstrand Figs. 3–6).

The :math:`2N`-dimensional spectrum is the **distinguishing feature
of Galerkin spectral relative to all other math-heart classes** in
this codebase:

* :class:`MomentSpace` (F_N) returns only the dominant
  :math:`(c, a)` root from :math:`\det M(a) = 0` — by construction.
* :class:`Billiard` (trajectory_resolvent) returns only the
  power-iteration fundamental :math:`k_{\rm eff}` from the
  Birkhoff resolvent — by construction.
* :class:`BasisSpace` returns the full :math:`2N`-eigenvalue
  spectrum.

The complex-conjugate pairs are physically meaningful. They
correspond to oscillatory higher-mode eigenfunctions; their
appearance threshold (in :math:`(\bar\mu, d)` space) marks the
transition from non-oscillatory to oscillatory residual modes
in the Carlvik integral operator's spectrum. Dahl–Sjöstrand
Figs. 3–6 plot the spectrum over the full
:math:`(\bar\mu, d)` parameter sweep; this kind of spectral
diagnostic is reproducible from
:meth:`BasisSpace.solve_full_spectrum`.

Eigenvalue spectrum convergence and truncation error
-----------------------------------------------------

The convergence rate of Galerkin spectral on the bare-critical
transport problem is set by the smoothness class of the true
scalar flux on the rescaled interval :math:`[-1, +1]`. For
bare-critical configurations:

* The scalar flux is *analytic* on the open interval (no internal
  singularities; the operator is elliptic in the sense that
  :math:`E_1` smooths).
* The endpoint behaviour is logarithmic (the well-known Placzek
  endpoint) — :math:`\phi(\pm 1) \sim \log(1 \mp x)` as
  :math:`x \to \pm 1`.

The endpoint logarithm reduces the convergence rate from
super-algebraic (which Legendre would deliver for a fully analytic
function) to a slower algebraic decay. In practice, the matrix
elements :math:`A_{m,n}` and :math:`B_{m,n}` already absorb the
endpoint structure into the kernel evaluation, so the *effective*
truncation error decays much faster than the endpoint analysis
would predict on the trial functions alone. Dahl–Sjöstrand 1979
Tables I, II demonstrate this empirically: at :math:`N = 9` the
fundamental eigenvalue :math:`c_{\rm crit}` reaches 7
significant figures across the full
:math:`(\bar\mu, d)` sweep tabulated; at :math:`N = 16` the matrix
becomes ill-conditioned (the Legendre columns become
near-linearly-dependent in finite precision), so refinement
saturates rather than continuing to gain digits.

The **practical operating point** is :math:`N = 9`,
:math:`n_{\rm quad} = 128`. These are the values
:class:`BasisSpace` defaults to and that
:func:`solve_galerkin_spectral_slab` /
:func:`solve_galerkin_spectral_sphere` use to reproduce
Dahl–Sjöstrand Tables I, II.

Slab vs sphere reductions: the Bickley-Naylor kernel family
------------------------------------------------------------

The slab and sphere differ by exactly two structural choices:

1. **Basis parity** (V_cg.4): even-Legendre for slab (scalar flux
   :math:`\phi(x)` symmetric), odd-Legendre for sphere (reduced flux
   :math:`r\,\phi(r)` antisymmetric).
2. **Boundary-chord parameter** :math:`q`: :math:`q = 0` for slab,
   :math:`q = 1` for sphere. This sets the rank-1 :math:`B_{m,n}`
   structure.

The kernel :math:`E_1(a|x - y|)` is the same in both geometries
because the angular pre-integration produces the same Peierls form
(the angular kernel for the slab Boltzmann equation and the
spherical Boltzmann equation, when reduced to the :math:`r\phi`
variable, share the :math:`E_1` Bickley-Naylor function). The
**Bickley-Naylor family** :math:`\{E_n, Ki_n\}` parametrises
angular-pre-integrated kernels for slab, sphere, and cylinder:

* **Slab**: :math:`E_n(\tau) = \int_0^1 \mu^{n-1} e^{-\tau/\mu}\,
  d\mu`, the canonical exponential-integral family. :math:`E_1`
  is the Peierls kernel; :math:`E_3` is its second antiderivative
  and appears in the boundary-chord term.
* **Sphere** (after :math:`u = r\phi` substitution): the slab
  :math:`E_1` family applies directly, with :math:`q = 1`
  swapping the parity of the basis.
* **Cylinder**: requires the 3-D Bickley-Naylor :math:`Ki_n`
  family — :math:`Ki_n(\tau) = \int_0^{\pi/2} (\cos\theta)^{n-1}
  e^{-\tau/\cos\theta}\, d\theta`. This is structurally different
  from :math:`E_n` and is **why cylinder is out of pillar** for
  Galerkin spectral. Westfall–Metcalf 1972 documents that the
  Mitsis-style Wiener–Hopf reduction the Carlvik form would
  rely on is non-convergent for the bare cylinder; the cylinder
  benchmark ships through the Mitsis–Westfall–Metcalf Fredholm
  pillar in :mod:`...singular_eigenfunction.cylinder`.

The unification of slab and sphere under the same Galerkin
machinery — sharing the same Eq. (4) block matrix, differing only
in basis parity and :math:`q` — is the architectural payoff of the
formulation. The Branch-2 production code in
:mod:`...core.galerkin_matrix` exposes a single
:func:`solve_eq4_eigenproblem` taking ``geometry="slab"`` or
``geometry="sphere"`` and producing the right matrix; the
geometry-specific entry points
(:func:`solve_galerkin_spectral_slab` /
:func:`solve_galerkin_spectral_sphere`) are thin wrappers.

Distinction from MomentSpace: same moment-basis idea, different orthogonality
------------------------------------------------------------------------------

:class:`BasisSpace` and :class:`MomentSpace`
(:mod:`...fn_method.moment_space`) both attack the same physical
problem with **moment-basis projections**, but they differ
structurally on three axes:

1. **Half-range vs full-range basis.** ``MomentSpace`` projects the
   *boundary angular flux* onto the half-range monomial basis
   :math:`\{\mu^\alpha\}_{\alpha=0}^N` on :math:`[0, 1]`.
   ``BasisSpace`` projects the *interior scalar flux* onto the
   full-range Legendre basis :math:`\{P_n(x)\}_{n=0}^{2(N-1)}` on
   :math:`[-1, +1]`. These are different function spaces; the
   projection operators are not commensurable.

2. **Projection vs collocation.** ``MomentSpace`` collocates the
   moment-basis identity at :math:`N+1` discrete :math:`\xi`-points
   (Grandjean–Siewert grid for slab; Siewert–Thomas Chebyshev grid
   for sphere) — Petrov–Galerkin with delta-function tests, NOT
   true Galerkin. ``BasisSpace`` projects onto the same Legendre
   basis on both sides — true Galerkin with the trial = test
   coincidence that gives best-approximation in the energy norm.

3. **Dominant eigenvalue vs full spectrum.** ``MomentSpace``
   root-finds :math:`\det M(a) = 0` (slab) / :math:`\det M(R) = 0`
   (sphere) on the configuration parameter — a *single-eigenvalue*
   determinantal condition that yields only the fundamental.
   ``BasisSpace`` solves the standard non-symmetric eigenproblem
   :math:`M v = \lambda v` on the Eq. (4) block matrix — yielding
   ALL :math:`2N` eigenvalues including complex pairs.

The two methods cross-check each other at :math:`\bar\mu = 0`
(isotropic limit) on the bare-critical Sood ``*-1-0-SL/SP`` truth
values — agreement at :math:`10^{-6}` is strong evidence of
correctness because the methods descend from completely
different mathematical structures (full-range Legendre Galerkin
vs half-range moment collocation).

The Sanchez–Chandrasekhar three-meanings taxonomy
--------------------------------------------------

The Sanchez–Chandrasekhar three-meanings taxonomy (see
:doc:`/theory/reference_solvers` § three-meanings) locates each of
the math-heart classes within the same Green's-function landscape
of one-speed neutron transport:

* :class:`MomentSpace` (F_N) is the **spectral half-range
  realisation** of the resolvent — the angular flux is projected
  onto the Case singular-eigenfunction half-range moments, then
  the resolvent's spectral kernel is collocated at discrete
  :math:`\xi` points. The eigenvalue falls out of a small dense
  determinantal condition.
* :class:`Billiard` (trajectory_resolvent) is the **path-integral
  (time-domain) realisation** — the Birkhoff transfer operator
  :math:`S` advances trajectories one bounce period at a time;
  its resolvent :math:`T = (I - S)^{-1} = \sum_n S^n` is the
  geometric-series sum over bouncing characteristics weighted by
  streaming attenuation.
* :class:`BasisSpace` (Galerkin spectral) is the **full-range
  polynomial-basis realisation** — the resolvent's interior is
  projected onto a Legendre basis on the spatial variable, with
  the angular variable pre-integrated analytically into the
  kernel before the spatial Galerkin projection. The eigenvalue
  spectrum falls out of a standard non-symmetric matrix
  eigenproblem.

All three are different angles of attack on the same Boltzmann
equation; the algebra-of-record discipline (see
:file:`.claude/skills/algebra-of-record/SKILL.md`) ensures that the Branch-2
implementations of each share NO in-house code primitives above
the trusted-library line. Their cross-checks anchor the Sood
``*-1-0-SL/SP`` and ``*-1-1-SL/SP`` truth values from three
structurally-independent directions.

The shared cross-method result type
:class:`~orpheus.derivations.common.solution_types.CriticalSolution`
is the load-bearing piece of the unification. Each reference
solver returns the same ``CriticalSolution`` shape regardless of
pillar; the ``eigenvalue_kind`` field disambiguates ``"k_eff"``
(MomentSpace, Billiard) from ``"c_critical"`` (BasisSpace) so
cross-method comparators read the right field before comparing.

The pre-Phase-D :class:`TransportSolver` Protocol (in
``orpheus.derivations.common.solver_protocol``) was retired in the
architectural reset — it conflated continuous reference generators
with discrete production solvers, which have functionally different
roles. The reference solvers now consume :class:`StructuredGeometry`
directly via their frozen ``__init__``; the production solvers
consume ``(materials, mesh, params)`` via the canonical free
functions ``solve_cp`` / ``solve_sn`` / ``solve_moc``.

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
