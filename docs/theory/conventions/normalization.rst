.. _theory-conventions-normalization:

=============
Normalization
=============

:doc:`notation` fixes *which* symbol means what. This page fixes the
second axis of the wrong-by-a-constant failure class: *what size*
everything is. Three normalization families cover every constant-level
trap recorded in this codebase — the :term:`quadrature` **weight sum**
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
  moment is a plain weighted sum over :term:`ordinates <ordinate>` — the projection
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

**The canon trap.** Hébert :cite:`Hebert2009` §3.9.1 normalizes 1-D
Gauss–Legendre to :math:`\sum w = 2`; five pages later,
Eqs. (3.363)–(3.364) normalize to :math:`\sum w = 1` over the positive
octant — the same symbol :math:`w_n`, two normalizations, no note.
Any weight-bearing formula must be read against its *page-local*
:math:`\sum w` before import.

**The recorded cost — ERR-025.** The 1-D :term:`diamond-difference <diamond difference>` fast
path's precomputed recurrence needed

.. math::
   :label: normalization-dd-source-coefficient

   b \;=\; \frac{2\,\Delta x\,(Q/W)}{2\mu + \Delta x\,
   \Sigma_{\mathrm{t}}},


.. implements:: normalization-dd-source-coefficient
   :by: orpheus.transport.operators.scattering.ScatteringOperator._assemble_per_ordinate_source

   **Implemented by** 5 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: normalization-dd-source-coefficient
   :by: orpheus.transport.spatial.diamond.DiamondDifference.affine_scan_coefficients

.. implements:: normalization-dd-source-coefficient
   :by: orpheus.transport.spatial.diamond.DiamondDifference.update

.. implements:: normalization-dd-source-coefficient
   :by: orpheus.transport.spatial.scheme.DiscretizationSchemeBase.source_emission

.. implements:: normalization-dd-source-coefficient
   :by: orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence

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
(:ref:`vv-principles anti-patterns 3–4 <theory-verification>`).

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
     - :math:`(2\ell+1)/W` — the **Parseval metric**
     - An adjoint is metric-RELATIVE: :math:`\Pi^*` is fixed only once
       *both* inner products are named. The frame's coefficient
       codomain carries :math:`G^{-1}`, the inverse discrete Gram
       (:ref:`frame-parseval-metric`,
       :doc:`/theory/foundations/frame`), which on a degree-exact
       sphere rule IS this ledger's own unifying object
       :math:`(2\ell+1)/W`. Hence
       :math:`(\Pi^* c)_m = \sum_\ell \frac{2\ell+1}{W}
       \sum_{m'} Y_\ell^{m'}(\hat\Omega_m)\,c_{\ell m'}
       = (R c)_m / W` —
       :eq:`hilbert-adjoint-equals-metric-times-S0`. ⛔ **Two earlier
       readings of this row are retired** (2026-08-23, step F-0): "the
       naked reconstruction, **nothing**" is the adjoint under a
       *Euclidean* coefficient metric, and ":math:`g_C \cdot S_0`" the
       adjoint under the *continuum* Gram — the frame installs
       neither.
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

**The unification, once more, from the other end.** The ledger's
:math:`(2\ell+1)/W` is not only a bookkeeping convenience: it is the
frame's **Parseval metric** :math:`G^{-1}`, the inverse of the
discrete trial Gram :math:`G_\ell = W/(2\ell+1)`
(:ref:`frame-parseval-metric`). That is why the same object appears
in the adjoint row and in the :math:`\Pi R` row — the two are
:math:`d_\ell G_\ell = W` read in opposite directions.

.. warning::

   The identification :math:`G_\ell = W/(2\ell+1)` is a property of a
   **degree-exact sphere cubature**, not of the basis. `[M]`
   2026-08-23 on the slab ``gauss_legendre(8)`` measure at
   :math:`L = 2` the discrete Gram is **not diagonal at all** (largest
   live off-diagonal :math:`0.93` of :math:`\sqrt{G_{jj}G_{kk}}`), so
   no :math:`(2\ell+1)/W` diagonal is its inverse and the frame
   refuses the Parseval dressing there. Hébert's :math:`W = 2` row
   remains the right reading of the *prefactor*; it is not a claim
   about the slab frame's metric. See
   :ref:`frame-parseval-dense-refusal`.

**The catchers.** In
``tests/numerics/test_spherical_harmonic_space.py``:
``test_H_equals_parseval_metric_times_S0`` pins the Hilbert adjoint to
:math:`S_0(G^{-1}c)` at :math:`10^{-12}` — the Parseval metric
:math:`(2\ell+1)/4\pi` times the *naked* synthesis, with no bare
:math:`(2\ell+1)` on the adjoint side (ERR-039);
``test_R_equals_2l_plus_1_times_S0`` pins the complementary side,
that the :math:`(2\ell+1)` lives in :math:`R`;
``test_pi_R_is_4pi_identity_on_band_limited`` pins
:math:`\Pi R = 4\pi I` (ERR-051). In ``tests/numerics/test_frame.py``
the ``test_parseval_*`` family pins the metric itself — the isometry
:math:`\|\Pi\psi\|_{G^{-1}} = \|\psi\|_W`, the closure
:math:`\Pi^* = R/W` and :math:`R^* = W\,\Pi` over six sphere families,
the slab ``DENSE`` refusal, and a loaded-not-blind negative leg that
re-installs the pre-F-0 continuum metric. The addition theorem itself
is verified at :math:`\ell \le 3` in
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
across every row: :math:`\alpha_{1/2} = 0` seeds the recursion,
:math:`\alpha_{M+1/2} = 0` closes it at the far end, and
the resulting system is lower-triangular (a direct one-pass :term:`sweep`).
The two zeros are not the same kind of statement — the seed is an axiom of
every spelling, the closure is a *theorem about the quadrature* and hence a
contract a bad rule can violate
(:ref:`sn-alpha-dome-closes`).

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
   * - Stacey (9.213) :cite:`Stacey2007`
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - \mu_m w_m`
     - :math:`\sum w = 2`
     - :math:`\Delta A / w` form (the Carlson lineage)
   * - Bell & Glasstone (5.21) :cite:`BellGlasstone1970`
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - \mu_m w_m (A_{i+1} - A_i)`
     - —
     - **folded**: :math:`\Delta A` lives inside :math:`\alpha`
   * - Hébert sphere (3.424) :cite:`Hebert2009`
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
   * - Bailey–Morel–Chang sphere (11) :cite:`BaileyMorelChang2010`
     - :math:`\alpha_{m+\frac{1}{2}} = \alpha_{m-\frac{1}{2}}
       - 2\mu_m w_m`
     - :math:`\sum w = 2`
     - :math:`(r/w_m)\,[\alpha\psi]_{m\pm\frac{1}{2}}`
   * - Bailey–Morel–Chang R–Z (50)
     - :math:`\alpha_{m+\frac{1}{2},n} = \alpha_{m-\frac{1}{2},n}
       - \mu_{m,n} w_{m,n}` per :math:`\xi`-level
     - :math:`\sum\sum w = 4\pi`
     - :math:`1/(r\,w_{m,n})`
   * - Larsen–Morel review (1.23b) :cite:`LarsenMorel2010` — symbol
       :math:`\beta`
     - :math:`\beta_{n+\frac{1}{2}} = \beta_{n-\frac{1}{2}}
       - 2\mu_n w_n` — **identical to the BMC sphere**
       :math:`\alpha` with the sign absorbed
     - :math:`\sum w = 2`
     - :math:`(r/w_n)` in the :math:`r^2`-form (1.23a); ⚠ their
       :math:`\alpha` is the **spatial** :term:`weighted-diamond <weighted diamond difference>` weight
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
