.. _verification-reference-solutions:

Reference Solutions — the ORPHEUS Verification Contract
========================================================

This page is the **binding contract** for every verification
reference solution in ORPHEUS. It establishes (i) the vocabulary
discipline around the words *verification*, *reference solution*,
and *benchmark*; (ii) the operator-form taxonomy every reference
commits to; (iii) the :class:`ContinuousReferenceSolution` dataclass
that carries mesh-independent analytical or semi-analytical
solutions into tests; and (iv) the kernel-primitive identities that
underpin the whole Phase-4 (Peierls) and Phase-5 (analytical
transport) infrastructure.

.. contents::
   :local:
   :depth: 2


.. _vv-vocabulary:

Vocabulary Discipline
---------------------

We follow the tight V&V definitions of Oberkampf & Roache. In this
repository the following words mean:

**Verification**
   The mathematical exercise of proving that a code solves its
   governing equations correctly, by comparison against a reference
   that is derived **from the governing equations themselves**.
   Verification stands alone — it does not require any other code.

**Validation**
   Comparison of a code's output against experiment (ICSBEP,
   IRPhE, …). Out of scope for this page; see the V&V-level
   taxonomy in the project ``CLAUDE.md`` — verification is
   levels L0/L1/L2, validation is level L3, benchmarking is
   level L4.

**Reference solution**
   An analytical or semi-analytical function of the independent
   variables — :math:`\phi(x, g)`, :math:`\psi(x, \mu, g)`,
   :math:`k_{\text{eff}}`, whatever the solver claims to compute
   — derived by pure mathematics (SymPy, mpmath, closed-form
   algebra) from the governing equation. *Mesh-independent*:
   evaluable at any point to arbitrary precision without running
   any solver.

**Benchmark** (forbidden in verification contexts)
   A code-to-code comparison at level L4 of the V&V ladder.
   **Never** used as a verification
   artefact in this repository. Legacy collections such as
   Sood, Parsons & Forster 1999 (LA-13511) and Ganapol 2008
   use the word "benchmark" in their titles for historical reasons,
   but their *contents* are analytical reference solutions derived
   from the transport equation; we cite them as such. When the word
   appears in a literature citation, it always refers to the
   collection's title, never to an ORPHEUS verification artefact.

This discipline is load-bearing: future session agents will read
this page first and must never conflate the two categories.


.. _operator-form-taxonomy:

Operator-Form Taxonomy
----------------------

Every :class:`~orpheus.derivations.ContinuousReferenceSolution`
commits to exactly **one** *operator form* — the mathematical form
of the governing equation it solves. Tests that consume a reference
solution assert that the target solver discretises the same form.
A reference solution in the ``"differential-sn"`` form has **no
business** being consumed by a diffusion test, because diffusion
solves a different equation.

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Tag
     - Equation
     - Target solvers
   * - ``"homogeneous"``
     - :math:`k_\infty = \lambda_{\max}(\mathbf A^{-1}\mathbf F)`
       (infinite medium; no spatial variable)
     - All solvers in the homogeneous-material degenerate limit
   * - ``"differential-sn"``
     - :math:`\mu_n\,\partial\psi_n/\partial x + \Sigma_t\,\psi_n
       = \mathrm{RHS}` (discrete ordinates)
     - :mod:`orpheus.sn`
   * - ``"differential-moc"``
     - :math:`d\psi/ds + \Sigma_t\,\psi = Q` along characteristics
     - :mod:`orpheus.moc`
   * - ``"diffusion"``
     - :math:`-\nabla\!\cdot\!(D\nabla\phi) + \Sigma_r\,\phi = S`
     - :mod:`orpheus.diffusion`
   * - ``"integral-peierls"``
     - :math:`\phi(\tau) = \tfrac{c}{2}\!\int E_1(|\tau-\tau'|)\phi(\tau')\,d\tau' + S(\tau)`
     - :mod:`orpheus.cp`
   * - ``"stochastic-transport"``
     - Integro-differential Boltzmann, sampled by random walks
     - :mod:`orpheus.mc`

The taxonomy is **disjoint by design** — a given reference solution
targets one operator. Homogeneous infinite-medium references are
the one exception: they are degenerate in space and can be consumed
by any solver as a sanity check on the multigroup matrix algebra.


.. _continuous-reference-contract:

The ContinuousReferenceSolution contract
-----------------------------------------

:class:`orpheus.derivations.ContinuousReferenceSolution` is a frozen
dataclass that carries a mesh-independent reference solution into
tests. The contract it exposes:

- **Callable fields** ``phi(x, g)`` (always) and ``psi(x, mu, g)``
  (optional for angle-resolved ansätze). These are closures over
  SymPy or mpmath state and can be evaluated at arbitrary points.
- **Convenience**: :meth:`~orpheus.derivations.ContinuousReferenceSolution.phi_on_mesh`
  for cell-centre evaluation and
  :meth:`~orpheus.derivations.ContinuousReferenceSolution.phi_cell_average`
  for Gauss–Legendre cell averages (needed when comparing against
  a finite-volume solver like CP or diffusion FV where the solver
  output is a cell average, not a point value).
- **Problem specification**: a :class:`~orpheus.derivations.ProblemSpec`
  dataclass with materials, geometry, BCs, optional external source,
  and an ``is_eigenvalue`` flag.
- **Provenance**: a :class:`~orpheus.derivations.Provenance` record
  with literature citation, derivation notes, SymPy string, and
  mpmath working precision.
- **Operator form**: one of the tags in
  :ref:`operator-form-taxonomy`, asserted by consumers.

The full dataclass definition lives in
:mod:`orpheus.derivations._reference`.


.. _reference-solution-registry:

Registry and lookup
-------------------

Two parallel registries coexist during the Phase-0 → Phase-6
migration:

- :func:`orpheus.derivations.get` — legacy
  :class:`~orpheus.derivations._types.VerificationCase` by name.
  Carries only the scalar ``k_inf`` and problem spec. Every
  currently-green test consumes this.
- :func:`orpheus.derivations.continuous_get` — new
  :class:`~orpheus.derivations.ContinuousReferenceSolution` by
  name. Populated incrementally as each existing derivation is
  retrofitted to the Phase-0 contract (see the migration plan).

Both registries use the same key namespace. A derivation that has
been retrofitted registers its continuous reference under the same
name as its legacy case, and
:meth:`~orpheus.derivations.ContinuousReferenceSolution.as_verification_case`
bridges back to the legacy type so existing tests do not break.


.. _kernel-primitives:

Kernel primitives — :math:`E_n` and :math:`\mathrm{Ki}_n`
---------------------------------------------------------

The Phase-4 Peierls slab / cylinder / sphere references rest on two
families of special functions. They are the atomic building blocks
of every integral-transport reference solution, so their
identities and special values are verified as **L0
term-verification tests** in
:mod:`tests.derivations.test_kernels` before any higher-level
reference is built on top of them.

Exponential integral :math:`E_n(x)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canonical definition (Abramowitz & Stegun 5.1.4):

.. math::
   :label: en-definition

   E_n(x) \;=\; \int_1^{\infty} \frac{e^{-xt}}{t^n}\,dt,
   \qquad x \ge 0,\; n \ge 1.

.. math::
   :label: en-kernel-special-values

   E_n(0) \;=\; \frac{1}{n - 1} \quad (n > 1),
   \qquad E_1(0) \;=\; +\infty \text{ (log singularity)}.

.. math::
   :label: en-kernel-derivative

   E_n'(x) \;=\; -E_{n-1}(x).

.. math::
   :label: en-kernel-integral

   \int_0^{\infty} E_n(x)\,dx \;=\; \frac{1}{n}.

Evaluated to arbitrary precision by
:func:`orpheus.derivations._kernels.e_n` (which wraps
:func:`mpmath.expint`) and by :func:`scipy.special.expn` at double
precision. Both engines are exercised against the three identities
above in
:func:`tests.derivations.test_kernels.test_en_closed_form_at_zero`,
:func:`~tests.derivations.test_kernels.test_en_derivative_identity`,
and
:func:`~tests.derivations.test_kernels.test_en_full_line_integral`.


Bickley–Naylor function :math:`\mathrm{Ki}_n(x)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canonical definition (Bickley & Naylor 1935; A&S 11.2):

.. math::
   :label: kin-definition

   \mathrm{Ki}_n(x) \;=\; \int_0^{\pi/2}
       \cos^{n-1}\theta\;\exp\!\left(-\frac{x}{\cos\theta}\right)\,d\theta.

The integrand has an essential singularity at :math:`\theta = \pi/2`
when :math:`x > 0`. ORPHEUS's high-precision evaluator
:func:`orpheus.derivations._kernels.ki_n` resolves this by the
substitution :math:`u = \tan\theta`, which produces a smooth
integrand on :math:`[0, \infty)`:

.. math::

   \mathrm{Ki}_n(x) \;=\; \int_0^{\infty}
       (1 + u^2)^{-(n+1)/2}\,
       \exp\!\bigl(-x\sqrt{1 + u^2}\bigr)\,du.

Special values at :math:`x = 0` follow from the Wallis integrals:

.. math::
   :label: kin-kernel-special-values

   \mathrm{Ki}_1(0) \;=\; \tfrac{\pi}{2},\quad
   \mathrm{Ki}_2(0) \;=\; 1,\quad
   \mathrm{Ki}_3(0) \;=\; \tfrac{\pi}{4},\quad
   \mathrm{Ki}_4(0) \;=\; \tfrac{2}{3},\quad
   \mathrm{Ki}_5(0) \;=\; \tfrac{3\pi}{16}.

Derivative identity (A&S 11.2.11, with the convention
:math:`\mathrm{Ki}_0(x) = K_0(x)`, the modified Bessel function):

.. math::
   :label: kin-kernel-derivative

   \mathrm{Ki}_n'(x) \;=\; -\mathrm{Ki}_{n-1}(x).

Both identities are exercised in
:func:`tests.derivations.test_kernels.test_kin_closed_form_at_zero`,
:func:`~tests.derivations.test_kernels.test_kin_derivative_identity`,
and (for :math:`n = 1`, where the derivative hits the modified
Bessel function)
:func:`~tests.derivations.test_kernels.test_kin1_derivative_is_bessel_k0`.


Legacy naming discrepancy in ``BickleyTables`` (historical)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   **Retired in Phase B.4 (commit 6badbe5, Issue #94).** The
   ``BickleyTables`` class and its ``bickley_tables()`` cache
   function have been deleted from
   :mod:`orpheus.derivations._kernels`. The canonical replacement is
   :func:`~orpheus.derivations._kernels.ki_n_mp` (arbitrary precision,
   A&S naming) and, for double-precision fast paths, the Chebyshev
   interpolant :func:`orpheus.derivations.cp_geometry._ki3_mp`.
   This subsection is retained as **historical documentation** —
   both the tabulation error and the naming convention shaped the
   V&V ladder for the life of the project, and the record of *why*
   the retirement was safe matters for future sessions.

Until Phase B.4 the cylindrical CP modules relied on the legacy
tabulation ``BickleyTables``, whose method names were **off-by-one**
from the Abramowitz & Stegun numbering:

- ``BickleyTables.ki3(x)`` integrated
  :math:`\sin(t)\,\exp(-x/\sin t)` which equals
  :math:`\mathrm{Ki}_2^{\text{A\&S}}(x)` under the substitution
  :math:`\theta = \pi/2 - t` — **not**
  :math:`\mathrm{Ki}_3^{\text{A\&S}}`.
- ``BickleyTables.ki4(x)`` computed a cumulative-trapezoid
  approximation to :math:`\mathrm{Ki}_3^{\text{A\&S}}` on a
  20 000-point grid with ~:math:`10^{-3}` absolute accuracy. The
  cumsum trick (``np.cumsum(ki3[::-1])[::-1] * dx``) was a
  pre-computer-era idiom for realising
  :math:`\mathrm{Ki}_{n+1}(x) = \int_x^\infty \mathrm{Ki}_n(t)\,dt`
  with a single pass over a uniform grid; it kept one antiderivation
  cheap on 1960s hardware at the cost of :math:`O(\Delta x^2)`
  trapezoidal error that propagated through every
  :math:`P_{ij}` evaluation.

The retirement was numerically safe because the downstream P-matrix
assembly in :mod:`orpheus.derivations.cp_cylinder` (and its solver
sibling in :mod:`orpheus.cp.solver`) consumed the legacy ``ki4`` value
under the canonical alias ``Ki3_vec`` (introduced in a Phase-4.2
compatibility shim). The chord-form second-difference formula had
always been :math:`\Delta^2[\mathrm{Ki}_3^{\text{A\&S}}]`; the
legacy naming was purely internal. Swapping
``BickleyTables.ki4_vec`` for the Chebyshev interpolant
:func:`~orpheus.derivations.cp_geometry._ki3_mp` therefore replaced a
~:math:`10^{-3}`-accurate approximation of the canonical
:math:`\mathrm{Ki}_3` with a ~:math:`5\times 10^{-6}`-accurate one —
no convention change, just a precision upgrade of the kernel already
in use.

Cylinder :math:`k_\infty` reference values shifted by up to
~:math:`4\times 10^{-4}` for multi-region 1-group cases as the
legacy tabulation's bias was removed (shift toward the exact mpmath
result, not regression). The nine ``cp_cyl1D_*`` L1 eigenvalue tests
now pass at their declared ``< 1e-5`` tolerance with actual error
~:math:`10^{-7}` — ~100× headroom.

The legacy-convention regression guards
``test_legacy_bt_ki3_equals_kin2`` and
``test_legacy_bt_ki4_approximates_kin3`` that pinned the off-by-one
naming during Phase 0 through Phase B.3 were deleted alongside the
class in commit 6badbe5; they are listed here only as a pointer for
readers tracing through older commits in ``git log``.


.. _verification-campaign-migration:

Verification campaign — audit and migration plan
--------------------------------------------------

Each existing module in :mod:`orpheus.derivations` was audited at
the start of the campaign (Phase 0) and classified by the tier of
its reference:

.. list-table::
   :header-rows: 1
   :widths: 25 25 50

   * - Module
     - Current tier
     - Migration target
   * - ``homogeneous.py``
     - T1 analytical (matrix eigenvalue)
     - Phase 1.1 — expose flat eigenvector as ``phi``.
   * - ``sn.py`` homogeneous
     - T1 analytical
     - Phase 1.1 — fold into ``homogeneous.py``.
   * - ``sn.py`` heterogeneous
     - ~~T3~~ **RETIRED** (Phase 2.1a)
     - Replaced by MMS continuous references in ``sn_mms.py``.
       ``solver_cases()`` deleted from ``sn.py``.
   * - ``sn_heterogeneous.py``
     - T2 semi-analytical (transfer matrix, eigenvalue only)
     - Phase 1.2 — expose continuous :math:`\phi_g(x)` via back-substitution.
   * - ``cp_slab.py``
     - T2.5 (E₃-based P-matrix eigenvalue)
     - Phase 4.1 — replace with Peierls Nyström reference using :func:`e_n`.
   * - ``cp_cylinder.py``
     - T2.5 (legacy Bickley table; naming audit pending)
     - Phase 4.2 — replace with Peierls reference using :func:`ki_n`.
   * - ``cp_sphere.py``
     - T2.5 (E₃-based P-matrix eigenvalue)
     - Phase 4.3 — replace with Davison slab-reduction reference.
   * - ``diffusion.py``
     - T1/T2 (bare slab buckling, transfer-matrix 2-region)
     - Phase 1.3 — expose piecewise-smooth :math:`\phi_g(x)`.
   * - ``moc.py`` homogeneous
     - T1 analytical
     - Phase 1.1 — fold into ``homogeneous.py``.
   * - ``moc.py`` heterogeneous
     - ~~T3~~ **RETIRED** (Phase 2.2a)
     - Replaced by MMS continuous reference in ``moc_mms.py``.
   * - ``mc.py`` homogeneous
     - T1 analytical
     - Phase 1.1 — fold into ``homogeneous.py``.
   * - ``mc.py`` heterogeneous
     - T4 **BANNED** (CP derivation used as reference for MC)
     - Phase 2.3 — replace with Case/Placzek or Peierls continuous reference.
   * - ``sn_mms.py``
     - T1.5 analytical MMS (6 continuous references)
     - Phases 2.1a, 3.1–3.5 complete. Covers: slab 2G hetero,
       2D Cartesian 1G/2G, spherical, cylindrical, P1 anisotropic.
       Phases 3.3–3.4 blocked by sweep-operator inconsistency (#98).
   * - ``moc_mms.py``
     - T1.5 analytical MMS (1 continuous reference)
     - Phase 2.2a complete. Pin-cell MMS spatial-operator test.

Tiers:

T1
   Analytical closed form — fully symbolic, no numerical quadrature.
T2
   Semi-analytical to root-finder tolerance — the reference is the
   limit of a finite mpmath computation.
T2.5
   Semi-analytical but uses a discretisation that overlaps with
   the solver under test (P-matrix with the same approximation
   family). Functional for now, replaced in Phase 4.
T3
   **Methodologically banned** — uses the solver under test
   as the reference (Richardson extrapolation). Phase 2 deletes
   these after drop-in replacements are in place.
T4
   **Methodologically banned** — uses a different solver as the
   reference (cross-solver crutch). Phase 2 deletes these too.

See the verification-campaign plan for the full phased execution
order (Phase 0 scaffolding → Phase 1 retrofit → Phase 2 replacement
→ Phase 3 differential MMS → Phase 4 Peierls → Phase 5 analytical
transport → Phase 6 capstone report).


References
----------

- Oberkampf, W. L. and Roy, C. J., *Verification and Validation in
  Scientific Computing*, Cambridge University Press, 2010.
- Roache, P. J., *Verification and Validation in Computational
  Science and Engineering*, Hermosa Publishers, 1998.
- Case, K. M. and Zweifel, P. F., *Linear Transport Theory*,
  Addison-Wesley, 1967.
- Davison, B., *Neutron Transport Theory*, Oxford, 1957.
- Bell, G. I. and Glasstone, S., *Nuclear Reactor Theory*,
  Van Nostrand Reinhold, 1970.
- Abramowitz, M. and Stegun, I. A., eds., *Handbook of Mathematical
  Functions*, NBS, 1964, §5 (exponential integrals), §11 (Bickley).
- Bickley, W. G. and Naylor, J., "A Short Table of the Functions
  :math:`\mathrm{Ki}_n(x)`," *Philosophical Magazine*, ser. 7,
  **20** (1935) 343.
- Sood, A., Forster, R. A. and Parsons, D. K., *Analytical Benchmark
  Test Set for Criticality Code Verification*, LA-13511, 1999
  (later journal version: *Prog. Nucl. Energy* **42** (2003) 55,
  DOI 10.1016/S0149-1970(02)00098-7). The word "benchmark" in the
  title refers to the published collection; the *contents* are
  analytical reference solutions in the sense of this page.
- Ganapol, B. D., *Analytical Radiation Transport Benchmarks for
  the Next Century*, LA-UR-05-4848, 2005. (Same note on vocabulary.)
