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
transport) infrastructure; and (v) the reference kernel's dense linear
algebra, the eigen and source solves every matrix-reducing reference
stands on (:ref:`verification-reference-kernel`).

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
   IRPhE, …). Out of scope for this page; see the
   :ref:`V&V-level ladder <vv-level-ladder>` — verification is
   levels L0/L1/L2, validation is level L3, benchmarking is
   level L4 (informational, never in place of L0–L3).

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
   Sood, Forster & Parsons 2003 and Ganapol 2008
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


.. _verification-greens-three-meanings:

Three meanings of "Green's function" in this verification suite
----------------------------------------------------------------

When a reference is described as "Green's-function-based" in the V&V
matrix, the description is ambiguous: three structurally-independent
mathematical objects all carry that name in the transport literature,
and this suite consumes all three. The taxonomy is documented in
detail at :ref:`reference-solvers-three-meanings`. In V&V terms:

* **Meaning (α) — trajectory resolvent.** Constructs the scalar
  Green's kernel :math:`G(\rho \to \rho')` by tracing characteristic
  rays and closing multi-bounce trajectories with the resolvent
  :math:`T = (I - S)^{-1}`. Realised in
  :mod:`orpheus.derivations.continuous.trajectory_resolvent`. Pillar:
  semi-analytical (ray-traced quadratures + geometric series).

* **Meaning (β) — spectral resolvent.** Constructs the *same* scalar
  kernel via closed-form spectral μ-integration of the within-medium
  angular Green's function (Sanchez 1986 Eq. A6 / PS-1982 Eq. 21).
  **Reserved, not yet implemented** — the folder
  ``orpheus.derivations.continuous.spectral_resolvent`` is a README-only
  placeholder, so it is named as a literal rather than a ``:mod:`` role;
  the direct evaluator is the headline implementation gap. Pillar:
  semi-analytical (closed-form integrand + 1-D quadrature).

* **Meaning (γ) — singular-eigenfunction angular Green's.**
  Constructs the angular Green's function
  :math:`G(\tau, \tau'; \mu, \mu')` directly via Case ν-spectrum +
  X-function half-range completeness. Realised in
  :mod:`orpheus.derivations.continuous.singular_eigenfunction` and
  :mod:`orpheus.derivations.continuous.fn_method`. Pillar:
  closed-form (criticality determinant) / semi-analytical (interior
  flux reconstruction via KLL 1974 Fredholm iteration).

The verification claim of highest confidence in this suite is the
**triple match** (α) ≈ (β) ≈ (γ) — agreement at the same problem
across three structurally-distinct integrands. This is L1-grade
evidence per the structural-independence rule of the
``vv-principles`` skill. The matrix in :doc:`matrix` indicates which
problems currently realise the triple match (today: vacuum-BC sphere
isotropic; spectrum will widen as Meaning β is implemented).


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
:mod:`orpheus.derivations.common.continuous_reference`.


.. _reference-solution-registry:

Registry and lookup
-------------------

Two parallel registries coexist during the Phase-0 → Phase-6
migration:

- :func:`orpheus.derivations.get` — legacy
  :class:`~orpheus.derivations.common.verification_case.VerificationCase` by name.
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
:mod:`tests.gates.derivations.test_kernels` before any higher-level
reference is built on top of them.

Exponential integral :math:`E_n(x)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canonical definition (Abramowitz & Stegun 5.1.4):

.. math::
   :label: en-definition

   E_n(x) \;=\; \int_1^{\infty} \frac{e^{-xt}}{t^n}\,dt,
   \qquad x \ge 0,\; n \ge 1.

.. (vv-status rationale) definition: the canonical A&S 5.1.4 defining integral
.. of the exponential integral. Definitional — the implementing evaluator
.. ``kernels.e_n`` is pinned by the derived identities (special values,
.. derivative, full-line integral) in ``tests.gates.derivations.test_kernels``.
.. vv-status: en-definition documented

.. math::
   :label: en-kernel-special-values

   E_n(0) \;=\; \frac{1}{n - 1} \quad (n > 1),
   \qquad E_1(0) \;=\; +\infty \text{ (log singularity)}.

.. math::
   :label: en-kernel-derivative

   E_n'(x) \;=\; -E_{n-1}(x).


.. implements:: en-kernel-derivative
   :by: orpheus.derivations.common.kernels.e_n

   **Implemented by** 2 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: en-kernel-derivative
   :by: orpheus.derivations.common.kernels.e_n_derivative

.. math::
   :label: en-kernel-integral

   \int_0^{\infty} E_n(x)\,dx \;=\; \frac{1}{n}.

Evaluated to arbitrary precision by
:func:`orpheus.derivations.common.kernels.e_n` (which wraps
:func:`mpmath.expint`) and by :func:`scipy.special.expn` at double
precision. Both engines are exercised against the three identities
above in
:func:`tests.gates.derivations.test_kernels.test_en_closed_form_at_zero`,
:func:`~tests.gates.derivations.test_kernels.test_en_derivative_identity`,
and
:func:`~tests.gates.derivations.test_kernels.test_en_full_line_integral`.


Bickley–Naylor function :math:`\mathrm{Ki}_n(x)`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Canonical definition (Bickley & Naylor 1935; A&S 11.2):

.. math::
   :label: kin-definition

   \mathrm{Ki}_n(x) \;=\; \int_0^{\pi/2}
       \cos^{n-1}\theta\;\exp\!\left(-\frac{x}{\cos\theta}\right)\,d\theta.

.. (vv-status rationale) definition: the canonical Bickley-Naylor / A&S 11.2
.. defining integral. Definitional — the implementing evaluator
.. ``kernels.ki_n`` (via the ``u = tan(theta)`` substitution) is pinned by the
.. derived special-value and derivative identities in
.. ``tests.gates.derivations.test_kernels``.
.. vv-status: kin-definition documented

The integrand has an essential singularity at :math:`\theta = \pi/2`
when :math:`x > 0`. ORPHEUS's high-precision evaluator
:func:`orpheus.derivations.common.kernels.ki_n` resolves this by the
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
:func:`tests.gates.derivations.test_kernels.test_kin_closed_form_at_zero`,
:func:`~tests.gates.derivations.test_kernels.test_kin_derivative_identity`,
and (for :math:`n = 1`, where the derivative hits the modified
Bessel function)
:func:`~tests.gates.derivations.test_kernels.test_kin1_derivative_is_bessel_k0`.


Legacy naming discrepancy in ``BickleyTables`` (historical)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   **Retired in Phase B.4 (commit 6badbe5, Issue #94).** The
   ``BickleyTables`` class and its ``bickley_tables()`` cache
   function have been deleted from
   :mod:`orpheus.derivations.common.kernels`. The canonical replacement is
   :func:`~orpheus.derivations.common.kernels.ki_n_mp` (arbitrary precision,
   A&S naming) and, for double-precision fast paths, the Chebyshev
   interpolant :func:`orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp`.
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
assembly in :mod:`orpheus.derivations.continuous.flat_source_cp.cylinder` (and its solver
sibling in :mod:`orpheus.cp.solver`) consumed the legacy ``ki4`` value
under the canonical alias ``Ki3_vec`` (introduced in a Phase-4.2
compatibility shim). The chord-form second-difference formula had
always been :math:`\Delta^2[\mathrm{Ki}_3^{\text{A\&S}}]`; the
legacy naming was purely internal. Swapping
``BickleyTables.ki4_vec`` for the Chebyshev interpolant
:func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` therefore replaced a
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


.. _verification-reference-kernel:

The reference kernel's linear algebra — the dense pencil and the least solution
-------------------------------------------------------------------------------

A continuous reference that has reduced its transport problem to a few
hundred unknowns asks two questions of the resulting matrices: the
eigenpairs of a pencil (the fundamental mode, the higher modes, the
adjoint modes), and the solution of a source problem. Both are answered by
one small module, :mod:`orpheus.derivations.common.dense_pencil`, on numpy
and scipy only. It is a kernel primitive in the same sense as
:math:`E_n` and :math:`\mathrm{Ki}_n` above: every reference that reaches a
matrix stands on it, so it is verified on its own, against manufactured
matrices with chosen spectra, before any reference uses it.

**Why the references carry their own.** Production has an eigen-extraction
(:func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`) and a pencil
(:class:`~orpheus.numerics.pencil.OperatorPencil`), built for versatility
and changed whenever production grows. A reference is closed, purpose
built and limited in scope, and is useful only while it stays put, so it
computes with machinery of its own that changes only by ruling; the
principle, in the user's words, and the gate that enforces it are on
:ref:`architecture-reference-insulation`. The Perron–Frobenius extraction
therefore exists twice, once on each side of the branch line. The two
copies are deliberately different in strength (below), and neither may be
merged into the other: a reference agreeing with production would then
share its extraction code, and the agreement would carry less information.

The pencil and its weak form
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let :math:`\{u_1, \dots, u_n\}` be the basis the reference assembled on, and
let :math:`a_L(\cdot,\cdot)` and :math:`a_F(\cdot,\cdot)` be the bilinear
forms of the loss operator :math:`\mathcal{L}` (streaming and collision,
less in-group and down-scattering where they are counted as loss) and the
production operator :math:`\mathcal{F}`, the first slot the test function
and the second the trial function:
:math:`a_L(v, w) = \langle v, \mathcal{L} w\rangle`. A
:class:`~orpheus.derivations.common.dense_pencil.DensePencil` holds the two
matrices of these forms,

.. math::

   L_{ij} = a_L(u_i, u_j), \qquad F_{ij} = a_F(u_i, u_j),

which is what a Galerkin assembly produces. Writing the trial function as
:math:`\phi = \sum_j c_j u_j` and testing
:math:`\mathcal{F}\phi = k\,\mathcal{L}\phi` against each :math:`u_i` gives
the generalised eigenproblem the kernel solves,

.. math::

   F\,c \;=\; k\,L\,c ,

by the QZ reduction (``scipy.linalg.eig(F, L)``). Its eigenvalues are those
of :math:`L^{-1}F`, the matrix of the operator
:math:`\mathcal{L}^{-1}\mathcal{F}` restricted to the trial space, whatever
the basis's Gram matrix is: the Gram matrix multiplies both sides and
cancels.

**Why weak form, and why the adjoint is then the transposed pair.** The
adjoint operator is defined by swapping the slots of each form,
:math:`a_{L^\dagger}(v, w) = a_L(w, v)`, and so on for :math:`F`. Its
matrices on the same basis are therefore

.. math::
   :label: reference-kernel-adjoint-pencil

   (L^\dagger)_{ij} = a_L(u_j, u_i) = L_{ji},
   \qquad
   (F^\dagger)_{ij} = F_{ji},
   \qquad\text{so}\qquad
   (L^\dagger, F^\dagger) = (L^{\mathsf T}, F^{\mathsf T}),

with the same spectrum, since
:math:`\det(F^{\mathsf T} - k L^{\mathsf T}) = \det(F - kL)`.
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.adjoint` returns
exactly this pair, and no metric appears. Had the kernel stored the operator
in strong form instead, as the coefficient matrix :math:`P = W^{-1}L` of
:math:`\mathcal{L}` with :math:`W_{ij} = \langle u_i, u_j\rangle` the Gram
matrix (or the quadrature weights of a nodal rule), the Euclidean transpose
:math:`P^{\mathsf T}` would not be the adjoint: the adjoint of :math:`P` in
the inner product :math:`x^{\mathsf T} W y` is :math:`W^{-1} P^{\mathsf T} W`,
since :math:`(Px)^{\mathsf T} W y = x^{\mathsf T} W (W^{-1}P^{\mathsf T}W) y`.
For :math:`P = W^{-1}L` that is :math:`W^{-1}L^{\mathsf T}`, whose weak form
is :math:`L^{\mathsf T}` again. The weak form keeps the metric inside the
forms, so the only adjoint spelling the kernel offers is the correct one.

.. implements:: reference-kernel-adjoint-pencil
   :by: orpheus.derivations.common.dense_pencil.DensePencil.adjoint

Two identities follow from :eq:`reference-kernel-adjoint-pencil`. Let
:math:`F v_j = k_j L v_j` and :math:`F^{\mathsf T} w_i = k_i L^{\mathsf T}
w_i`. Then :math:`w_i^{\mathsf T} F v_j` equals both
:math:`k_j\, w_i^{\mathsf T} L v_j` and :math:`k_i\, w_i^{\mathsf T} L v_j`,
so for distinct eigenvalues the forward and adjoint modes are
biorthogonal in the loss form (and, multiplying through, in the production
form):

.. math::
   :label: reference-kernel-biorthogonality

   (k_i - k_j)\, w_i^{\mathsf T} L\, v_j = 0
   \quad\Longrightarrow\quad
   w_i^{\mathsf T} L\, v_j = 0 \;\;(k_i \ne k_j).

And a source solve is reciprocal through the adjoint: if
:math:`(L - F)\,x = s` and :math:`(L - F)^{\mathsf T} x^\dagger = d`, then

.. math::
   :label: reference-kernel-reciprocity

   \langle d, x\rangle
   = (x^\dagger)^{\mathsf T}(L - F)\,x
   = \langle x^\dagger, s\rangle .

**The nodal-basis requirement.** The sign test below reads the
eigenvector's coefficients. On a nodal (Lagrange) basis a coefficient is
the function's value at its node, so single-signed coefficients are
single-signed node values. On a modal basis (Legendre polynomials, say) the
coefficients' signs say nothing about the function's sign, in either
direction, and the test would refuse or admit for the wrong reason. A
reference that builds a :class:`~orpheus.derivations.common.dense_pencil.DensePencil`
on a modal basis must change basis before asking for the fundamental mode.

The fundamental mode and the basis of each refusal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The fundamental mode of a criticality problem is defined by a theorem, not
by a selection rule. For a compact operator :math:`K` on a Banach space
ordered by a cone (here: non-negative functions), positive in the sense
that it maps the cone into itself, the **Krein–Rutman theorem** states that
the spectral radius :math:`\rho(K)`, if positive, is itself an eigenvalue
with an eigenvector in the cone; in its strong form (:math:`K` maps every
non-zero element of the cone into the cone's interior) the eigenvalue
:math:`\rho(K)` is simple, it is the only eigenvalue with an eigenvector in
the cone, and every other eigenvalue has modulus strictly below it. The
finite-dimensional case is the **Perron–Frobenius theorem** for
non-negative matrices: an irreducible one has a simple positive eigenvalue
:math:`\rho` with a positive eigenvector, unique up to scale among the
non-negative eigenvectors, and a primitive one (irreducible with period 1)
has every other eigenvalue strictly inside :math:`|\lambda| < \rho`. For
neutron transport :math:`K = \mathcal{L}^{-1}\mathcal{F}` is positive
because the transport inverse and the fission production both preserve
non-negative densities; Bell & Glasstone state the multigroup consequence
(a real, positive :math:`k_0` larger in magnitude than any other
eigenvalue, with a unique non-negative eigenfunction and adjoint
eigenfunction) and rest it on that positivity :cite:`BellGlasstone1970`
(§4.4c, pp. 189–190). Neither Krein and Rutman's paper nor a
matrix-analysis text stating the Perron–Frobenius theorem is in the local
literature folder, so the two statements above carry no theorem numbers;
Bell & Glasstone's is the one checked against its source here.

A discretisation inherits none of this automatically: the discrete
:math:`L^{-1}F` is a positive matrix only if the scheme preserves
positivity, and Bell & Glasstone point out that the multigroup
:math:`P_1` difference equations (§4.4g, p. 197) and the diamond-difference
discrete-ordinates equations (§5.2f, p. 225) need not. So
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.fundamental`
does not *select* the fundamental mode; it reads the eigenvalue of largest
modulus and checks that the matrix behaves as the theorem says a positive
operator's must. Each check is one consequence of the theorem, and its
failure is a :class:`~orpheus.derivations.common.dense_pencil.NoFundamentalMode`
naming that consequence:

.. math::
   :label: reference-kernel-fundamental-contract

   k_0 \in \mathbb{R},\qquad
   k_0 > 0,\qquad
   1 - \frac{|k_1|}{|k_0|} > \delta,\qquad
   \min_i c_{0,i} \;\ge\; -\varepsilon\,\max_i |c_{0,i}|,

with :math:`k_0, k_1` the two eigenvalues of largest modulus, :math:`c_0`
the eigenvector of :math:`k_0` oriented so its coefficients sum to a
non-negative number, :math:`\delta` =
:data:`~orpheus.derivations.common.dense_pencil.DOMINANCE_GAP`
(:math:`10^{-10}`) and :math:`\varepsilon` the pencil's own rounding band,
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.sign_band`
(below).

.. implements:: reference-kernel-fundamental-contract
   :by: orpheus.derivations.common.dense_pencil.DensePencil.fundamental

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - Refusal (message fragment)
     - The consequence it checks
     - What a violation means
   * - ``is complex``
     - :math:`\rho` is an eigenvalue, so the eigenvalue of largest modulus
       is real
     - the matrix is not positive. QZ returns a real eigenvalue of a real
       pencil with an exactly zero imaginary part, so the test is
       ``imag != 0.0``. A complex dominant eigenvalue of a real pencil has
       its conjugate at equal modulus, so the dominance check would refuse
       it too; the complex check runs first and owns the message.
   * - ``is not positive``
     - :math:`\rho > 0` is an eigenvalue
     - a negative eigenvalue leads in modulus (in the gate's fixture,
       :math:`-3` ahead of :math:`2`), so :math:`\rho` is not an eigenvalue
       and the matrix is not positive
   * - ``is not strictly dominant``
     - :math:`\rho` is simple and every other eigenvalue lies strictly
       inside :math:`|\lambda| < \rho` (strong form; primitive matrix)
     - two eigenvalues of equal modulus lead: a double one, as two
       decoupled parts at the same :math:`k` give (the fixture
       :math:`(2, 2, \dots)`), or a pair of opposite sign, as a periodic
       (imprimitive) non-negative matrix gives (the fixture
       :math:`(2, -2, \dots)`). No unique fundamental mode exists
   * - ``is not single-signed``
     - the eigenvector of :math:`\rho` lies in the cone
     - the discretisation produced a sign-changing "fundamental", which is
       the signature of a scheme that is not positivity-preserving

**The bands.** The sign band is a rounding band, not a physical allowance.
A coefficient that is zero in exact arithmetic (a group no fission or
scattering reaches, a node on a vacuum wall) comes out of the QZ reduction
at the level of the backward error, with either sign: in the gate's fixture
an exact zero returns at :math:`-1.2\times10^{-16}`. How far it moves grows
with the conditioning of the loss form and shrinks with the dominance gap,
so the band is the pencil's own,
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.sign_band`:
:data:`~orpheus.derivations.common.dense_pencil.SIGN_SAFETY` (10) times the
first-order error bound of the eigenvector of :math:`A = L^{-1}F` under a
relative backward error :math:`\varepsilon_{\rm mach}` in each form,

.. math::

   \varepsilon \;=\; 10\,\frac{\varepsilon_{\rm mach}\left(\|L^{-1}\|\,\|F\|
   + \kappa(L)\,\|A\|\right)}{|k_0| - |k_1|},

with 2-norms (the eigenvector perturbation bound of a simple eigenvalue,
Golub and Van Loan, *Matrix Computations*, 4th ed.; the section number is
not checked against a local copy). `[M]` 2026-10-06: the first spelling, a
fixed band of :math:`10^{-10}`, refused 9 of 20 pencils with
:math:`\kappa(L)` near :math:`10^8`, whose exact zeros return at up to
:math:`5.6\times10^{-9}` (qa); the derived band admits the exact zeros of
20 of 20 pencils at each of :math:`\kappa(L) = 10^2, 10^6, 10^8` and
refuses a sign change of :math:`-10^{-3}` in 20 of 20 at each; at
:math:`\kappa(L) = 10^{10}` the rounding itself reaches :math:`10^{-4}` and
a :math:`-10^{-3}` change is refused in 15 of 20, the honest limit of
resolution (``scratch/characteristic_architecture/p05_probes/sign_band_probe.py``).
A negative coefficient within the band is therefore set to exactly ``0.0``
and the vector renormalised; outside it, the mode is refused. The returned
:class:`~orpheus.derivations.common.dense_pencil.FundamentalMode` checks its
own invariant at construction (:math:`k` finite and positive, the vector
finite, non-negative and of unit Euclidean norm, and read-only), so a value
of that type is a fundamental mode whoever built it; its scale is the
caller's to set. :data:`~orpheus.derivations.common.dense_pencil.DOMINANCE_GAP`
is the smallest relative gap in modulus between :math:`k_0` and
:math:`k_1` that counts as strict dominance.

**Higher modes.**
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.spectrum`
returns every eigenpair as a
:class:`~orpheus.derivations.common.dense_pencil.SpectralMode`, ordered by
decreasing modulus, complex and sign-changing ones included, with no
refusal: the theorem says nothing about them except that they are not
the fundamental mode. It refuses only a singular loss form, which shows as
an infinite or undefined eigenvalue of the QZ reduction.

**The two twins, side by side.** Production's
:func:`~orpheus.numerics.eigenvalue.dominant_eigenpair` takes the
eigenvalue with the largest real part of a materialised resolvent
:math:`A^{-1}F` and refuses only a complex one; it does not check
positivity, dominance or the sign of the vector, although its docstring
says the criticality contract is enforced
(`#580 <https://github.com/deOliveira-R/ORPHEUS/issues/580>`_). For a
well-posed problem the eigenvalue of largest real part and the eigenvalue
of largest modulus coincide, so production answers correctly there; the
reference kernel is the stricter of the two because a reference is the
instrument that must refuse rather than answer when the posing is wrong.

The least solution of a source problem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Read as a source problem, the same pencil is :math:`L x = F x + s`: the
loss form is the form the unknown is measured in, the production form is
the secondary production (scattering and fission), and :math:`s` is the
external source.
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.least_solution`
returns the sum of its Neumann series. With :math:`T = L^{-1}F`, the
equation is :math:`x = L^{-1}s + T x`; substituting it into itself
:math:`N` times gives
:math:`x = \sum_{n<N} T^n L^{-1}s + T^N x`, the expansion in collision
generations. The series

.. math::
   :label: reference-kernel-least-solution

   x \;=\; \sum_{n=0}^{\infty} \bigl(L^{-1}F\bigr)^{n} L^{-1} s
     \;=\; (L - F)^{-1} s
   \qquad\text{when}\qquad \rho\bigl(L^{-1}F\bigr) < 1

converges for every source exactly when :math:`\rho(T) < 1` (the Neumann
series of a matrix converges iff its spectral radius is below 1), and its
sum is :math:`(I - T)^{-1}L^{-1}s = (L(I - T))^{-1}s = (L - F)^{-1}s`. The
radius of :math:`T` is the modulus of the leading eigenvalue of the pencil,
the :math:`k` of the eigenproblem, so "subcritical" means the same thing in
both readings. The kernel does not apply this condition to the whole
system, though, but only to the part of it the source reaches.

.. implements:: reference-kernel-least-solution
   :by: orpheus.derivations.common.dense_pencil.DensePencil.least_solution

**The reach.** Say that unknown :math:`j` enters equation :math:`i` when
:math:`L_{ij} \ne 0` or :math:`F_{ij} \ne 0`. The reach :math:`R(s)` of a
source is the smallest set of unknowns that contains the support of
:math:`s` and every unknown whose equation contains an unknown already in
the set; :meth:`~orpheus.derivations.common.dense_pencil.DensePencil.reach`
returns it as a boolean mask, by a breadth-first closure. In the language
of the Frobenius normal form (the unknowns permuted so that the coupling
is block triangular, with its strongly connected classes as the diagonal
blocks) it is the union of the source's classes and every class
downstream of them. Write :math:`U` for the unknowns outside the reach. By
construction no equation in :math:`U` contains an unknown in :math:`R`,
so :math:`L_{UR} = F_{UR} = 0` and, ordering :math:`R` first,

.. math::
   :label: reference-kernel-reach

   \begin{pmatrix} L_{RR} - F_{RR} & L_{RU} - F_{RU} \\
                   0 & L_{UU} - F_{UU} \end{pmatrix}
   \begin{pmatrix} x_R \\ x_U \end{pmatrix}
   = \begin{pmatrix} s_R \\ 0 \end{pmatrix}.

:math:`L` and :math:`F` are both block upper triangular, so :math:`L^{-1}`
and :math:`T = L^{-1}F` are too, and a vector that vanishes on :math:`U`
is mapped by :math:`T` to one that vanishes on :math:`U`, with
:math:`(Tv)_R = L_{RR}^{-1}F_{RR}\,v_R`. Every term of the Neumann series
therefore vanishes on :math:`U`, and the least solution is

.. math::

   x_U = 0, \qquad
   x_R = \sum_{n=0}^{\infty} \bigl(L_{RR}^{-1}F_{RR}\bigr)^{n}
         L_{RR}^{-1} s_R
       = (L_{RR} - F_{RR})^{-1} s_R
   \quad\text{when}\quad \rho\bigl(L_{RR}^{-1}F_{RR}\bigr) < 1 .

That is what
:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.least_solution`
computes: it measures the radius of the reached block (the modulus of the
leading eigenvalue of that block's own pencil, by the same QZ reduction
as :meth:`~orpheus.derivations.common.dense_pencil.DensePencil.spectrum`),
refuses when it is not below 1, solves the reached block, and leaves the
rest at zero. A class of radius at least 1 that the source does not reach
is not a refusal: it carries zero. A source that reaches a class of radius
at least 1 is refused, whether or not it also reaches subcritical classes.

.. implements:: reference-kernel-reach
   :by: orpheus.derivations.common.dense_pencil.DensePencil.reach

   **Implemented by** 2 sites: the closure itself, and the least solution
   that restricts the solve to it.

.. implements:: reference-kernel-reach
   :by: orpheus.derivations.common.dense_pencil.DensePencil.least_solution

**Why it is the least solution.** When :math:`L^{-1}`, :math:`F` and
:math:`s` are non-negative, every partial sum is non-negative and
increasing, and any non-negative solution :math:`y` satisfies
:math:`y = \sum_{n<N} T^n L^{-1}s + T^N y \ge \sum_{n<N} T^n L^{-1}s` for
every :math:`N`. So the series' sum lies below every non-negative solution:
it is the least one, and it is the physical one (the neutrons of every
collision generation, counted once).

**Why the radius is checked, not trusted to the solve.** When
:math:`\rho(T) > 1` the matrix :math:`L - F` is usually still invertible,
and a direct solve returns a finite vector without complaint; but for an
irreducible :math:`T \ge 0` the matrix :math:`(I - T)^{-1}` is
non-negative only when :math:`\rho(T) < 1`, and the vector has negative
entries (the gate's control leg: :math:`\rho = 1.5`, a finite
solve with a negative entry). Not a single partial sum of the series is
negative, so that vector is not the limit of anything physical. The
kernel refuses with
:class:`~orpheus.derivations.common.dense_pencil.NoLeastSolution`, naming
the measured radius.

**An unreached class, and the zero source.** A zero source reaches
nothing, so its least solution is zero whatever the radius, even where
:math:`L - F` is singular; this is the empty case of the reach, not a
special case in the code. The case that matters arises in transport along
lines: a line that stays inside a void
region between perfectly reflecting walls carries no loss at all, so its
own :math:`T` has radius 1 and :math:`I - T` is singular, and it carries no
source either. It is a class of its own that no source reaches, so its
least solution is zero, the limit as the wall albedo tends to 1 from
below, while every other line is solved. That is why the characteristic
references' closure is ruled to be the least solution (the user's ruling
"Least solution", 2026-10-06), and why the radius is checked on the reach
and not on the whole system (the user's ruling of the same day): a check
on the whole gain would refuse every source in a geometry holding one
such line assembled into the same pencil.

**The margin.** The radius is computed, so a system exactly at criticality
reads 1 only to within the eigensolver's backward error, with either sign.
`[M]` 2026-10-06, 4 × 4 gains with a chosen spectrum
:math:`(1, 0.5, 0.2, -0.3)` built as :math:`V\Lambda V^{-1}`, 20 seeds of
200 draws: the computed radius spans :math:`1 - 7.2\times10^{-13}` to
:math:`1 + 1.8\times10^{-13}`, and a bare comparison ``radius < 1`` admits
473 of the 4000 (between 18 and 31 per seed); on the gate's own draw
stream it admits 23 of 200, of which 21 return a solution with entries up
to :math:`9.8\times10^{16}` and 2 raise an untyped
``numpy.linalg.LinAlgError`` on the singular :math:`I - G`. The kernel
therefore refuses every radius not below
:math:`1 -` :data:`~orpheus.derivations.common.dense_pencil.SUBCRITICAL_MARGIN`
(:math:`10^{-10}`), about 140 times the worst draw's distance from 1. A
system nearer to criticality than the margin would have a least solution
amplified by more than :math:`10^{10}` along its fundamental mode, which no
reference value survives.

.. dropdown:: What the kernel replaced, and why it failed
   :icon: history

   **The largest-real-part extraction.** Before the kernel, four reference
   sites formed the resolvent ``M = np.linalg.solve(A, F)`` and took its
   eigenvalue with the largest real part (``kinf_and_spectrum_homogeneous``
   and ``kinf_from_cp`` in :mod:`orpheus.derivations.common.eigenvalue`,
   ``derive_1rg`` and ``_bare_slab_spectrum`` in
   :mod:`orpheus.derivations.continuous.cases.diffusion`), and the two
   spectrum variants zeroed coefficients below ``1e-14`` in absolute value
   without refusing anything. That rule answers
   when it should refuse: on the gate's fixture with eigenvalues
   :math:`(-3, 2, 1, 0.5)` it returns :math:`k = 2`, the second eigenvalue
   in modulus, as if it were the fundamental mode, and it never checks
   dominance or the sign of the vector. The kernel reads the eigenvalue of
   largest modulus and refuses, with one message per broken condition.

   **The factor-order trap on the** :math:`k_\infty` **adjoint.** The
   infinite-medium adjoint spectrum
   (:func:`~orpheus.derivations.common.eigenvalue.kinf_and_adjoint_spectrum_homogeneous`)
   must be the dominant eigenvector of
   :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}`, not of
   :math:`(A^{-1}F)^{\mathsf T} = F^{\mathsf T}A^{-\mathsf T}`. The two
   matrices are similar (conjugate by :math:`A^{\mathsf T}`), so every
   :math:`k`-level check passes on both; but the eigenvectors differ, and
   for the rank-one :math:`F = \chi \otimes \nu\Sigma_f` the wrong product's
   dominant eigenvector is exactly :math:`\widehat{\nu\Sigma_f}`, carrying
   no information about :math:`A` at all. The first spelling used
   :math:`\operatorname{eig}(M^{\mathsf T})`; the S\ :sub:`N` daggered solve
   disagreed with it on first contact (campaign #276, gate P1.4), and the
   law was corrected to :math:`(A^{\mathsf T})^{-1}F^{\mathsf T}` by hand,
   with production's ``dominant_eigenpair`` doing the extraction. In the
   kernel the adjoint is
   :meth:`DensePencil.adjoint <orpheus.derivations.common.dense_pencil.DensePencil.adjoint>`,
   the transposed pair :eq:`reference-kernel-adjoint-pencil`, and a
   resolvent is never formed, so neither factor order can be written.

The gates
~~~~~~~~~

The kernel's gates are :file:`tests/gates/derivations/test_reference_kernel.py`,
every matrix in it domain-free (a pencil
:math:`F = L V \Lambda V^{-1}` with a non-symmetric, well-conditioned
:math:`L` and chosen eigenpairs, so :math:`F v_i = \lambda_i L v_i` holds by
construction; a gain with a chosen spectrum), as ``instrument-doctrine``
X4 requires of a primitive that references rest on: it is pinned against
closed forms, never against an in-domain solver.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Claim
     - Gate (``test_reference_kernel.py``)
   * - fundamental mode = the chosen dominant eigenpair (:math:`n = 5, 40`,
       band :math:`10^{-13}` relative on :math:`k`); the within-band clip
     - ``test_k1_the_fundamental_is_the_chosen_dominant_eigenpair``,
       ``test_k1_a_within_band_negative_coefficient_is_clipped_to_exactly_zero``
   * - the four refusals of :eq:`reference-kernel-fundamental-contract`,
       one input each, keyed to disjoint message fragments; a singular loss
     - ``test_k1r_*``
   * - :class:`~orpheus.derivations.common.dense_pencil.FundamentalMode`
       refuses an instance that breaks its invariant
     - ``test_a_fundamental_mode_holds_its_invariant_whoever_builds_it``
   * - the full spectrum, ordered by modulus (a complex pair placed
       between two reals, where an order by real part would misplace it)
     - ``test_k5_the_spectrum_holds_every_mode_by_decreasing_modulus_with_no_refusal``
   * - :eq:`reference-kernel-adjoint-pencil`,
       :eq:`reference-kernel-biorthogonality`,
       :eq:`reference-kernel-reciprocity`
     - ``test_k6i_the_adjoint_has_the_same_spectrum``,
       ``test_k6ii_adjoint_and_forward_modes_are_biorthogonal_in_the_loss_form``,
       ``test_k6iii_the_source_solve_is_reciprocal_through_the_adjoint``
   * - :eq:`reference-kernel-least-solution` against its closed form
       (unit mass, :math:`\rho = 0.5, 0.9, 0.999`) and against the
       explicitly summed series (non-identity mass)
     - ``test_k2_the_least_solution_with_unit_mass``,
       ``test_k2_the_least_solution_with_a_non_identity_mass_is_the_neumann_sum``
   * - :eq:`reference-kernel-reach`: the reach is the downstream closure
       of the source (a chain :math:`0 \to 1 \to 2` and an island); an
       unreached class of radius 1 carries zero, both when it is
       decoupled and when it feeds the reached class; a reached class of
       radius 1 is refused, alone or beside a subcritical one
     - ``test_k2_reach_is_the_downstream_closure_of_the_source``,
       ``test_k2_an_unreached_supercritical_class_carries_zero``,
       ``test_k2r_a_reached_supercritical_class_is_refused``
   * - the supercritical refusal, the radius measured against the mass, the
       margin on both sides, the unit radius, the zero source
     - ``test_k2r_*``
   * - the composite Gauss–Legendre rule
       (:func:`~orpheus.derivations.common.quadrature.composite_gauss_legendre`)
       exact on a piecewise degree-7 polynomial that jumps at every
       breakpoint; each panel's weights sum to its width
     - ``test_k3_*``

The biorthogonality gate is the adjoint's tooth: on its fixture
:math:`L` is non-symmetric and the eigenvectors are not
:math:`L`-orthogonal, so an adjoint that returned the untransposed pencil
leaves off-diagonal entries of order :math:`10^{-1}`. The spectrum gate
alone could not catch that mutation, since a pencil and its transpose have
the same eigenvalues.

`[M]` 2026-10-06, a mutation battery over the kernel's first spelling (the
source solve then a free function of a mass, a gain and a source; textual
mutants of the module patched onto the live objects, the production file
unchanged by ``shasum``): picking the largest real part reddens the
negative-dominant refusal; removing the dominance guard, the sign guard or
the clip, or widening the then-fixed sign band to ``1e-2``, reddens its own row;
ordering the spectrum by real part reddens three rows; an untransposed
adjoint reddens the biorthogonality and reciprocity rows; removing the
radius guard, the zero-source shortcut, or measuring the radius of the
gain without the mass each reddens its ``k2r`` row; the positive control
(solving :math:`L + F` for :math:`L - F`) reddens every least-solution
row. The unreached-class rows' first red is the radius check over the
whole gain, which refuses both of them. Two arms were blind by
construction: removing the complex guard is
caught only by the message fragment, because the dominance guard refuses
a complex pair anyway, and removing the eigenvector normalisation changes
nothing, because ``scipy.linalg.eig`` already returns unit-norm columns.


.. _verification-pillar-2-hardening:

Pillar-2 reference hardening — Atkinson product Nyström as canonical
---------------------------------------------------------------------

Reference solvers in this project are continuously **hardened** —
their evaluation precision is improved to push the L1 cross-check
floor downward — even after a reference is "shipped" and consumed
by tests. The canonical hardening pattern, established by the
Atkinson product Nyström + tanh-substitution work in commit
``4c83e09`` (ERR-036 + ERR-037), is the load-bearing example a future
session should follow when extending or improving any Pillar-2
reference.

The pattern has three steps:

1. **Identify the precision-limiting structure.** The Peierls
   integral kernel :math:`E_1(|\rho - \rho'|)` is log-singular at the
   diagonal. A naïve Gauss-Legendre Nyström quadrature converges
   slowly because the singularity is uniformly bounded but
   non-smooth. ERR-036 captures the fingerprint of the resulting
   precision floor in the verification matrix.
2. **Apply a structurally-justified transformation.** Atkinson's
   product Nyström rule uses pre-computed weights for the singular
   factor, leaving a smooth kernel for standard quadrature. The
   tanh substitution :math:`\mu = \tanh(t)` (ERR-037) further
   smooths the BC-edge integrands by mapping the boundary singularity
   to :math:`t \to \infty`, where Gauss-Hermite-like exponential
   convergence applies.
3. **Pin the hardening with a regression test.** Each ERR entry
   carries ``@pytest.mark.catches("ERR-NNN")`` on the test that
   exercises the failure mode. The matrix in :doc:`matrix` flags
   any catalog entry without a catcher as a publication-blocker.

This is distinct from a *code-bug fix* (ERR catalog entries 1–35
mostly fall in that class). A hardening entry like ERR-036/037
records that the math was right but the *quadrature* was a precision
floor; the fix is in the numerical method, not the algebra.

ERR-038 (Atalay 1997 paper-precision floor, characterised in the
same commit) is a third evidence-class still distinct from the first
two: the *published reference numbers themselves* are at a precision
floor that a higher-precision cross-check cannot improve. The
mitigation is to swap in a structurally-independent reference
(Atkinson-hardened Pillar-2 evaluation of the same problem), which
is what ERR-036 enables.

The three classes — code bug, quadrature floor, paper-precision
floor — populate the error catalog with different lessons and
require different fixes. A future session reading the catalog must
attend to which class an ERR belongs to before trying to "fix" it.


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
T1.5
   Analytical MMS — the reference solution is imposed symbolically
   and exact by construction (T1-grade provenance), but the problem
   is source-driven: it supports convergence-order and flux-shape
   claims, never eigenvalue claims (see the pillar capability table
   at :ref:`verification-evidence-classes`).
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
  Test Set for Criticality Code Verification*, *Prog. Nucl. Energy*
  **42** (2003) 55-106, DOI 10.1016/S0149-1970(02)00098-7; the edition
  the corpus cites (the earlier edition is the 1999 Los Alamos report
  LA-13511, :ref:`sood-registry-editions`). The word "benchmark" in the
  title refers to the published collection; the *contents* are
  analytical reference solutions in the sense of this page.
- Ganapol, B. D., *Analytical Radiation Transport Benchmarks for
  the Next Century*, LA-UR-05-4848, 2005. (Same note on vocabulary.)


.. The Architecture section below (with the registry and cross-section
   library subsections) was moved verbatim at task #10 stage V3 from
   the dissolved ``docs/theory/verification.rst``; headings re-leveled
   to this page's hierarchy. V5+ integrates it with the contract
   sections above.

Architecture
------------

Three interlinked systems flow from one source:

.. code-block:: text

   orpheus/derivations/      SymPy derivations (single source of truth — library)
        │
        ├──→ tests/gates/     pytest imports reference values, runs solvers
        │
        └──→ docs/_generated/ RST fragments with LaTeX + results tables

The ``orpheus/derivations/`` package is the **library** of analytical
and semi-analytical reference solutions: SymPy ``derive_*`` functions
paired with ``test_*`` pytest gates, importable by production code
(see e.g. :mod:`orpheus.derivations.discrete.sn.balance` cited from
``_sweep_1d_cumprod`` (the dissolved ``sweep.py``), or
:mod:`orpheus.derivations.continuous.peierls_nystrom.origins.specular` cited from
:func:`orpheus.derivations.continuous.peierls_nystrom.geometry.reflection_specular`).

Separately, ``scratch/derivations/`` (project root, **not** a Python
package) is the **workbench** — in-flight SymPy drafts, diagnostic
scripts (``scratch/derivations/diagnostics/``), and retired code with
"why archived" notes (``scratch/derivations/archive/``). Workbench
content has no production consumer and is not on the import path; it
graduates to the library by being lifted into ``orpheus/derivations/``
with a paired test (see GitHub Issue #95 for the ongoing migration).


.. _reference-values-lazy-registry:

Reference-Values Registry (Eager vs Lazy)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``orpheus/derivations/reference_values.py`` is the central registry
that every test and every Sphinx page reads from. It holds two
registries — the legacy ``VerificationCase`` table (retrieved by
``get("name")``) and the Phase-0 ``ContinuousReferenceSolution`` table
(``continuous_get("name")``) — both populated by walking the
``orpheus.derivations.*`` modules. **Every case is now analytical or
semi-analytical** (NumPy / SciPy / SymPy only, no solver imports), so
the walk is cheap. The one refinement is that a reference whose
*construction* is expensive is deferred until it is actually
requested:

**Eager tier — cheap analytical / semi-analytical cases**
   Built when the registry is first touched, by walking the
   ``continuous.analytical.homogeneous`` case, the
   ``continuous.cases.{sn, moc, mc, diffusion}`` cases, and the
   ``continuous.flat_source_cp.{slab, cylinder, sphere}`` modules
   (the legacy table) plus every module exposing a
   ``continuous_cases()`` producer (the Phase-0 table). These use only
   NumPy, SciPy, and SymPy, so building them all is fast. The 9-case
   CP grids (:ref:`nine-case-cp-grid`), the homogeneous matrix
   eigenvalues, the diffusion buckling and transcendental cases, and
   the MMS gates all land here.

**Lazy tier — expensive continuous references**
   A reference whose *construction* costs
   :math:`\mathcal{O}(\text{minutes})` — an adaptive-mpmath eigenvalue
   solve such as the Peierls Nyström references — opts into the
   ``continuous_case_builders()`` contract, registering a
   ``name -> thunk`` in ``_CONTINUOUS_BUILDERS`` instead of an
   eagerly-built solution. The thunk runs (and memoises) only when
   that exact name is first requested, so fetching *any* other
   reference does not pay the expensive build (Issue #212 — the full
   Peierls build previously fired on every lookup and looked like a
   solver hang).

``continuous_get("name")`` checks the eager ``_CONTINUOUS`` table
first; if the name is absent it calls the thunk from
``_CONTINUOUS_BUILDERS`` once, memoises the result, and returns it.
The practical effect is that **every test pays only for the reference
values it actually uses**.

.. note::

   **Historical — the retired Richardson side-channel.** Before #290
   the lazy tier had a different tenant: Richardson-extrapolated
   references computed *by the production solvers themselves* — the SN
   heterogeneous slab, the MOC pin cell, and the diffusion fuel +
   reflector case — keyed in a ``_LAZY_LOADERS`` dict and persisted in
   a JSON-backed on-disk cache. Those loaders imported the solver
   packages lazily precisely to dodge the circular dependency they
   created (a solver imports ``orpheus.derivations`` to read its own
   reference values). The side-channel was retired in stages: the SN
   arm at Phase 2.1a (superseded by MMS), MOC never regrew one, and
   the last member — the diffusion 2-region case — at #290 P6,
   together with the MATLAB-port island that computed it.
   ``_LAZY_LOADERS`` no longer exists; **no reference is now derived
   from a production solver under test.** The diffusion 2-region
   reference's successor is the transcendental
   :func:`~orpheus.derivations.continuous.cases.diffusion.derive_2rg_continuous`
   (see the retired-Richardson record at
   :ref:`diffusion-2region-richardson`).

.. _synthetic-xs-library:

Cross-Section Library
~~~~~~~~~~~~~~~~~~~~~

All verification cases use **abstract synthetic cross sections** from
``orpheus/derivations/common/xs_library.py``.  Four regions are defined, each with
{1G, 2G, 4G} variants:

- **Region A** (fissile): fuel-like, with fission and moderate scattering
- **Region B** (scatterer): moderator-like, strong scattering, no fission
- **Region C** (absorber): cladding-like, moderate absorption
- **Region D** (gap): thin, very low opacity

All cross sections satisfy the consistency relation
:math:`\Sigma_t = \Sigma_c + \Sigma_f + \sum_{g'} \Sigma_{s,g \to g'}`.

P1 scattering anisotropy
^^^^^^^^^^^^^^^^^^^^^^^^

Every region carries a P1 scattering matrix in addition to the P0
isotropic form. The P1 matrix is built as
:math:`\Sigma_{s,1}(g \to g') = \bar\mu \cdot \Sigma_{s,0}(g \to g')`,
where :math:`\bar\mu` is the mean lab-frame scattering cosine set
per region to reflect the target nuclide's mass:

.. csv-table::
   :header: Region, Role, Nuclide analogue, :math:`\bar\mu`, Rationale
   :widths: 10, 22, 22, 10, 36

   A, fissile, "heavy uranium/plutonium", 0.05, "A ~ 235 gives :math:`\bar\mu_{cm \to lab} \approx 2/(3A) \approx 0.003`; bumped to 0.05 to expose the anisotropic-correction code path without making it dominant."
   B, moderator, "light H / H\ :sub:`2`\ O", 0.60, "A = 1 is the textbook forward-peaked limit, :math:`\bar\mu_{lab} = 2/3` for s-wave elastic on H. Rounded to 0.60."
   C, cladding, "intermediate Zr-90", 0.10, "A = 90 places it between fuel and moderator; 0.10 gives a mild forward peak consistent with Zr's :math:`2/(3A) \approx 0.007` lab-frame, again slightly amplified for test coverage."
   D, gap, "light He / void", 0.30, "Low-density gas with light nuclei. 0.30 picks a middle value that is neither fully isotropic nor dominated by the streaming peak."

These values are physically motivated but NOT tuned to any specific
isotope. The point is to give every transport solver a non-trivial
:math:`\Sigma_{s,1}` matrix so that the P1 correction term is
exercised by the verification suite; the absolute values are not
intended to match real nuclear data. Production runs always use the
real library from ``orpheus/data/micro_xs/``.

The definitive values and per-group arrays live in
``orpheus/derivations/common/xs_library.py`` (see the ``_MU_BAR`` dict at
the top of the file). Changing them there automatically propagates
to every verification case because all derivation modules import
``make_mixture`` from the same module.

Region layouts:

- **1 region**: A only (homogeneous fissile medium)
- **2 regions**: A + B (fuel + moderator)
- **4 regions**: A + D + C + B (fuel + gap + cladding + moderator)
