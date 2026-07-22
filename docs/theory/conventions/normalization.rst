.. _theory-conventions-normalization:

=============
Normalization
=============

:doc:`notation` fixes *which* symbol means what. This page fixes the
second axis of the wrong-by-a-constant failure class: *what size*
everything is. Three normalization families cover every constant-level
trap recorded in this codebase — the quadrature **weight sum**
:math:`W`, the **harmonic prefactor** :math:`(2\ell+1)/4\pi`, and the
curvilinear :math:`\alpha`-**recursion scale**. For each family this
page states the ORPHEUS law once, shows where the canon spreads, and
names the test that bites.

.. _normalization-weight-sum:

The weight sum :math:`W`
========================

**The law.** :math:`W \equiv \sum_m w_m` is whatever the quadrature
measure says it is — :math:`W = 2` for Gauss–Legendre on
:math:`[-1, 1]`, :math:`W = 4\pi` for the sphere rules — and it is
never silently renormalized. Enforced:
``tests/numerics/test_quadrature_directional.py`` pins the
Gauss–Legendre sum to :math:`2` and the ``lebedev`` /
``level_symmetric`` / ``product`` sums to :math:`4\pi`; the octant
partition test pins that restriction preserves the total.

The measure is the *only* normalization carrier, which fixes where
:math:`W` may and may not appear:

- **Moments carry no** :math:`1/W` **and no** :math:`1/4\pi`. A flux
  moment is a plain weighted sum over ordinates — the projection
  :math:`\Pi` applies the weights and nothing else (next section).
- **An isotropic scalar-rate source enters the per-ordinate equation
  as** :math:`Q/W`, so that its weighted sum returns :math:`Q`. The
  division by :math:`W` happens **once**, at the reconstruction site
  of the scattering application — not inside the moments, not inside
  the basis.
- **Restriction needs no re-normalization**: the ``half_range_clean``
  trait (:doc:`/theory/foundations/discrete_measures`) means
  ``measure.restrict`` yields a valid sub-quadrature whose weights
  keep their meaning.

**The canon trap.** Hébert [Hebert2009]_ §3.9.1 normalizes 1-D
Gauss–Legendre to :math:`\sum w = 2`; five pages later,
Eqs. (3.363)–(3.364) normalize to :math:`\sum w = 1` over the positive
octant — the same symbol :math:`w_n`, two normalizations, no note.
Any weight-bearing formula must be read against its *page-local*
:math:`\sum w` before import.

**The recorded cost — ERR-025.** The 1-D diamond-difference fast
path's precomputed recurrence needed

.. math::
   :label: normalization-dd-source-coefficient

   b \;=\; \frac{2\,\Delta x\,(Q/W)}{2\mu + \Delta x\,
   \Sigma_{\mathrm{t}}},

and the buggy derivation dropped the :math:`1/W` *and* mis-signed the
companion coefficient's numerator. The two errors compensate: the
buggy fixed point is :math:`\psi = Q/(2\Sigma_{\mathrm{t}})`, and the
missing :math:`1/W = 1/2` (Gauss–Legendre!) rescales it back — so
every homogeneous eigenvalue test passed to machine precision, because
:math:`k` is invariant under a uniform rescaling of :math:`\phi`. At a
material interface the rescaling becomes side-dependent and the
cancellation breaks: :math:`k_{\mathrm{eff}}` shifted by
:math:`\sim 1.5 \times 10^{-2}` on the two-region slab. The moral is
structural: **weight-sum errors hide behind eigenvalue invariances**.
The catchers are the per-ordinate flat-flux residual and
heterogeneous multi-group cases — never homogeneous :math:`k` alone
(:ref:`vv-principles anti-patterns 3–4 <verification-index>`).

.. _normalization-prefactor:

The :math:`(2\ell+1)/W` prefactor ledger
========================================

Every constant in the harmonic machinery has exactly one home. The
ledger (laws stated and derived in
:doc:`/theory/foundations/spherical_harmonics`):

.. list-table::
   :header-rows: 1
   :widths: 30 46 24

   * - Object
     - Carries
     - Law
   * - The basis :math:`Y_\ell^m`
     - **Nothing** — the no-prefactor convention
     - :math:`\langle Y_\ell^m, Y_{\ell'}^{m'}\rangle =
       \frac{4\pi}{2\ell+1}\delta_{\ell\ell'}\delta_{mm'}`; the
       addition theorem is clean,
       :eq:`real-sh-addition-theorem`.
   * - The reconstruction :math:`R`
     - :math:`(2\ell+1)`
     - The scattering-source kernel is
       :math:`\sum_\ell (2\ell+1) \sum_m Y_\ell^m\,\phi^{\ell m}` —
       the factor sits *outside* the basis, where the transport
       literature writes it.
   * - The scattering application
     - :math:`1/W`, **once**
     - The reconstructed source is divided by the weight sum at the
       single application site
       (:class:`~orpheus.transport.operators.scattering.ScatteringOperator`).
   * - The projection :math:`\Pi`
     - The weights :math:`w_m`
     - A moment is :math:`\sum_m w_m Y_\ell^m(\hat\Omega_m)\,
       \psi_m` — nothing else.
   * - The adjoint :math:`\Pi^*`
     - **Nothing** — the naked reconstruction
     - Under the :math:`W`-weighted inner product
       :math:`\langle\psi,\varphi\rangle_V = \sum_m w_m \psi_m
       \varphi_m`, the adjoint is
       :math:`(\Pi^* c)_m = \sum_{\ell,m'} Y_\ell^{m'}
       (\hat\Omega_m)\,c_{\ell m'}` with **no** :math:`(2\ell+1)` —
       ERR-039's exact content.
   * - The pair :math:`\Pi R`
     - :math:`4\pi = W`
     - :math:`\Pi R = 4\pi\,I` on band-limited inputs — the
       :math:`(2\ell+1)` in :math:`R` times the
       :math:`4\pi/(2\ell+1)` orthogonality. Asserting
       :math:`\Pi R = I` instead was ERR-051.

**The unification the canon misses.** Hébert carries the prefactor as
:math:`4\pi` in (3.30) and as :math:`2` in (3.425) — silently tied to
dimensionality. Read through the measure, these are **one object**:
:math:`(2\ell+1)/W`, evaluated at :math:`W = 4\pi` (the sphere) and
:math:`W = 2` (1-D Gauss–Legendre). The canon splits by
dimensionality what the weight sum unifies; ORPHEUS's dimension-blind
spelling is why the same scattering operator serves slab, sphere, and
cylinder without a per-geometry prefactor branch.

**The catchers.** All in
``tests/numerics/test_spherical_harmonic_space.py``:
``test_H_equals_g_C_times_S0`` pins the Hilbert adjoint to
:math:`g_C \cdot S_0` at :math:`10^{-12}` — the metric
:math:`g_C = 4\pi/(2\ell+1)` times the *naked* synthesis, with no
bare :math:`(2\ell+1)` on the adjoint side (ERR-039);
``test_R_equals_2l_plus_1_times_S0`` pins the complementary side,
that the :math:`(2\ell+1)` lives in :math:`R`;
``test_pi_R_is_4pi_identity_on_band_limited`` pins
:math:`\Pi R = 4\pi I` (ERR-051). The addition theorem itself is
verified at :math:`\ell \le 3` in
``tests/sn/operators/test_solver_components.py``.

.. _normalization-alpha-crosswalk:

The :math:`\alpha`-recursion crosswalk
======================================

The curvilinear angular-redistribution coefficient satisfies **the
same recursion spelled four ways across three texts** — none of which
acknowledges another spelling exists — plus a fifth spelling in the
review literature under a *different letter*. The factor of two and
the sign **move between the recursion and the balance-equation
divisor**; neither half is meaningful alone. The invariants that hold
across every row: :math:`\alpha_{1/2} = 0` seeds the recursion, and
the resulting system is lower-triangular (a direct one-pass sweep).

.. list-table::
   :header-rows: 1
   :widths: 22 34 16 28

   * - Source
     - Recursion
     - Weights
     - Redistribution scale
   * - **ORPHEUS**
       (:ref:`derivation <sn-alpha-normalization>`)
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - \mu_m w_m`
     - :math:`\sum w = 2`
     - :math:`\Delta A_i / w_m`
   * - Stacey (9.213) [Stacey2007]_
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - \mu_m w_m`
     - :math:`\sum w = 2`
     - :math:`\Delta A / w` form (the Carlson lineage)
   * - Bell & Glasstone (5.21) [BellGlasstone1970]_
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - \mu_m w_m (A_{i+1} - A_i)`
     - —
     - **folded**: :math:`\Delta A` lives inside :math:`\alpha`
   * - Hébert sphere (3.424) [Hebert2009]_
     - :math:`\alpha_{n+\frac{1}{2}} = \alpha_{n-\frac{1}{2}}
       - 2\,\mathcal{W}_n \mu_n`
     - GL on :math:`[-1,1]`
     - :math:`\Delta S_i / (2\,\mathcal{W}_n)` in (3.428)
   * - Hébert cylinder (3.399)
     - :math:`\alpha_{p,q+\frac{1}{2}} = \alpha_{p,q-\frac{1}{2}}
       \mathbin{\boldsymbol{+}} \mathcal{W}_{p,q}\,\mu_{p,q}`
     - per level
     - :math:`-(1/\mathcal{W}_{p,q})[\,\cdot\,]` in (3.400) —
       **sign-flipped against its own sphere, four pages apart**
   * - Bailey–Morel–Chang sphere (11) [BaileyMorelChang2010]_
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - 2\mu_m w_m`
     - :math:`\sum w = 2`
     - :math:`(r/w_m)\,[\alpha\psi]_{m\pm\frac{1}{2}}`
   * - Bailey–Morel–Chang R–Z (50)
     - :math:`\alpha_{m+\frac{1}{2},n} = \alpha_{m-\frac{1}{2},n}
       - \mu_{m,n} w_{m,n}` per :math:`\xi`-level
     - :math:`\sum\sum w = 4\pi`
     - :math:`1/(r\,w_{m,n})`
   * - Larsen–Morel review (1.23b) [LarsenMorel2010]_ — symbol
       :math:`\beta`
     - :math:`\beta_{n+\frac{1}{2}} = \beta_{n-\frac{1}{2}}
       - 2\mu_n w_n` — **identical to the BMC sphere**
       :math:`\alpha` with the sign absorbed
     - :math:`\sum w = 2`
     - :math:`(r/w_n)` in the :math:`r^2`-form (1.23a); ⚠ their
       :math:`\alpha` is the **spatial** weighted-diamond weight
       (1.30)

Every row is correct *in its own convention*. The import rule follows:
**never move a recursion between texts without carrying its divisor**
— the invariant object is the redistribution *term* (recursion scale
:math:`\times` divisor), not the recursion alone. ORPHEUS's arrays
carry :math:`\alpha^{O} = \alpha^{H}/2` relative to Hébert's sphere
form, with the divisor :math:`\Delta A_i / w_n` absorbing the
difference — the equivalence is proven where the machinery lives,
in :ref:`sn-alpha-normalization`
(:doc:`/theory/methods/sn/curvilinear_one_group`), and the named
geometry factor :math:`\Delta A / w` has its own treatment there.

The boundary discipline
=======================

The three import questions — arrow, weight sum, prefactor — live at
:ref:`notation-import-boundary`. Each method's machine header
(:doc:`/theory/methods/sn/index` and siblings) restates the local
normalization answers, so a reader entering through any method page
crosses this part's conventions first.
