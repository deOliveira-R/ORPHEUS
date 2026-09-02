.. _spherical-harmonics:

===========================================
Real Spherical Harmonics on a Direction Set
===========================================

Real spherical harmonics :math:`Y_\ell^m(\hat\Omega)` are the Galerkin
basis ORPHEUS uses to project an :term:`angular flux` onto its Pℓ moments and
to reconstruct a :term:`per-ordinate <ordinate>` scattering source from those moments.
This page is the canonical home for the **convention**, the
**addition-theorem identity** that the convention is engineered to
make literal, and the cross-method use of the
:meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.evaluate`
evaluator on
:class:`~orpheus.numerics.basis.SphericalHarmonicBasis`.

The basis is **generic infrastructure** — it lives in
:mod:`orpheus.numerics.basis`, not :mod:`orpheus.sn`, because the same
:math:`Y_\ell^m` table is consumed by SN aniso scattering, by the PN
solver (when it lands; Grand Report v3 §10), and by MC adjoint moment
estimators. It is the synthesis (trial) side of the spherical-harmonic
:class:`~orpheus.numerics.frame.GalerkinFrame` (see :ref:`galerkin-projection`);
the analysis :math:`M` and reconstruction :math:`R` operators are the
frame's two faces, NOT standalone operator classes.

.. contents::
   :local:
   :depth: 2


Key Facts
=========

- The convention is the **no-:math:`4\pi/(2\ell+1)`-prefactor** real
  spherical harmonics (Lewis & Miller 1993, §4.7). The addition
  theorem in this convention reads

  .. (vv-status rationale) The addition-theorem identity is the
     load-bearing structural identity that the convention is
     designed to make literal. Verified at :math:`\ell \le 3` by
     ``tests/sn/operators/test_solver_components.py::TestAnisotropicScattering
     ::test_spherical_harmonics_addition_theorem_L3``.
  .. vv-status: real-sh-addition-theorem documented

  .. math::
     :label: real-sh-addition-theorem

     \sum_{m=-\ell}^{\ell}
     Y_\ell^m(\hat\Omega)\,Y_\ell^m(\hat\Omega')
     \;=\; P_\ell(\hat\Omega \cdot \hat\Omega'),

  with no :math:`(2\ell+1)/4\pi` factor on either side. The
  reconstruction kernel for the SN Pℓ scattering source is
  :math:`\sum_\ell (2\ell+1) \sum_m Y_\ell^m Y_\ell^m`, which is
  the canonical form when the addition theorem is unprefactored.

- The polar axis is :math:`\mu_x` (so :math:`\cos\theta = \mu_x`,
  :math:`\sin\theta = \sqrt{1 - \mu_x^2}`). Azimuth lives in the
  :math:`(\mu_y, \mu_z)` plane:
  :math:`\cos\phi = \mu_y/\sin\theta`,
  :math:`\sin\phi = \mu_z/\sin\theta`.

- For :math:`\ell \le 1` the values are hard-coded and
  bit-identical to the legacy MATLAB
  ``discreteOrdinatesPWR.m`` reference:
  :math:`Y_0^0 = 1`, :math:`Y_1^{-1} = \mu_z`, :math:`Y_1^0 = \mu_x`,
  :math:`Y_1^{+1} = \mu_y`.

- For :math:`\ell \ge 2` the formula uses
  :func:`scipy.special.lpmv` with the Condon–Shortley
  :math:`(-1)^m` phase removed and the norm
  :math:`\sqrt{2(\ell-m)!/(\ell+m)!}` for :math:`m \ne 0`.

- The evaluator returns an ``(N, L+1, 2L+1)`` array. Index
  ``Y[n, ℓ, ℓ+m]`` holds :math:`Y_\ell^m(\hat\Omega_n)`. The
  :math:`m`-axis is shifted by :math:`\ell` so the slice
  ``Y[n, ℓ, ℓ-ℓ : ℓ+ℓ+1]`` covers the :math:`2\ell+1`
  in-range entries; out-of-range slots (``|m| > ℓ``) are zero.


Why this convention
===================

ANSI / textbook real spherical harmonics ship two competing
normalisations:

.. list-table::
   :header-rows: 1
   :widths: 28 36 36

   * - Convention
     - Orthogonality (continuous)
     - Addition theorem
   * - **No prefactor** (this project)
     - :math:`\langle Y_\ell^m, Y_{\ell'}^{m'}\rangle =
       \frac{4\pi}{2\ell+1}\,\delta_{\ell\ell'}\delta_{mm'}`
     - :math:`\sum_m Y_\ell^m\,Y_\ell^m
       = P_\ell(\Omega\cdot\Omega')`
   * - Standard ANSI / Wikipedia
     - :math:`\langle Y_\ell^m, Y_{\ell'}^{m'}\rangle =
       \delta_{\ell\ell'}\delta_{mm'}`
     - :math:`\sum_m Y_\ell^m\,Y_\ell^m
       = \frac{2\ell+1}{4\pi}\,P_\ell(\Omega\cdot\Omega')`

The transport literature (Bell & Glasstone 1970, §1.6; Lewis &
Miller 1993, §4.7) **always** writes the Pℓ scattering reconstruction
in the form

.. math::
   :label: sh-pl-scattering-reconstruction

   q(\hat\Omega) \;=\; \sum_{\ell} (2\ell+1)
     \sum_m Y_\ell^m(\hat\Omega)\,\phi^{\ell m},

.. (vv-status rationale) Literature-transcribed definition: the Pℓ scattering
   reconstruction with the (2ℓ+1) factor outside the basis (Bell & Glasstone
   1970 §1.6; Lewis & Miller 1993 §4.7). The reconstruction face R it defines is
   pinned by the R∘Λ∘M kernel crosscheck
   (``tests/sn/operators/test_scattering_kernel_crosscheck.py``) and
   :eq:`sh-addition-theorem-reconstruction`. A convention definition, not a
   solver claim.
.. vv-status: sh-pl-scattering-reconstruction documented

with the :math:`(2\ell+1)` factor **outside** the spherical-harmonic
basis. Adopting the no-prefactor convention puts the
:math:`(2\ell+1)` factor where the equations want it and removes
:math:`4\pi/(2\ell+1)` factors from the addition-theorem identity. The
discrete projection / reconstruction pair then carries no leftover
constants.

.. warning::

   The SciPy function :func:`scipy.special.sph_harm_y` (the successor
   to ``scipy.special.sph_harm``, which SciPy deprecated in 1.15 and
   removed in 1.17 — note the swapped ``(n, m)`` argument order)
   returns the **complex** spherical harmonics in the standard ANSI
   normalisation,
   :math:`Y_n^m = \sqrt{\tfrac{2n+1}{4\pi}\tfrac{(n-m)!}{(n+m)!}}\,
   P_n^m(\cos\theta)\,e^{im\phi}`. ORPHEUS does NOT consume that
   function. The
   :meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.evaluate`
   method builds real :math:`Y_\ell^m` from
   :func:`scipy.special.lpmv` (associated Legendre values
   :math:`P_\ell^m(\cos\theta)`) so the convention can be controlled
   directly. Mixing the two conventions in one codebase is the
   canonical way to introduce convention-drift bugs (failure mode 6
   in the V&V skill); the project's defense is "one evaluator, one
   convention, on
   :class:`~orpheus.numerics.basis.SphericalHarmonicBasis`".


.. _spherical-harmonics-eigenbasis:

Why spherical harmonics — the eigenbasis of the scattering kernel
=================================================================

The convention above answers *which* normalisation; this section
answers the prior question — *why these functions at all*. The
spherical harmonics are not merely a convenient orthogonal basis on
the sphere: they are the **eigenbasis of the anisotropic scattering
operator**, forced by the rotational symmetry of the scattering
kernel. This is a theorem, not an analogy, and it is the reason the
spherical-harmonic frame is a
:class:`~orpheus.numerics.frame.GalerkinFrame` *owned by the
scattering operator*.

The anisotropic scattering source is the integral operator

.. math::
   :label: sh-aniso-scattering-operator

   (S_{\rm aniso}\,\psi)(\hat\Omega)
   \;=\; \int_{4\pi}
         \Sigma_s(\hat\Omega \cdot \hat\Omega')\,\psi(\hat\Omega')\,
         d\hat\Omega',

.. (vv-status rationale) Literature-transcribed transport definition: the
   anisotropic scattering integral operator (a zonal kernel on S²), the same
   operator as frame.rst :eq:`scattering-zonal-kernel`. Its implementing kernel
   R∘Λ∘M is pinned by the 0-ULP crosscheck
   ``tests/sn/operators/test_scattering_kernel_crosscheck.py``. A definition,
   not a solver claim.
.. vv-status: sh-aniso-scattering-operator documented

whose kernel depends on the directions only through the cosine
:math:`\hat\Omega\cdot\hat\Omega'` — a **zonal** kernel on
:math:`S^2`. By the **Funk–Hecke theorem**, the spherical harmonics
are the eigenfunctions of any zonal kernel, with an eigenvalue that
depends on :math:`\ell` only:

.. math::
   :label: sh-funk-hecke-eigenvalue

   S_{\rm aniso}\,Y_\ell^m \;=\; \Sigma_{s,\ell}\,Y_\ell^m,
   \qquad
   \Sigma_{s,\ell} \;=\; 2\pi\!\int_{-1}^{+1}
         \Sigma_s(t)\,P_\ell(t)\,dt,

.. (vv-status rationale) Classical result: the Funk–Hecke eigenvalue of the
   zonal scattering kernel (Müller 1966), the same identity as frame.rst
   :eq:`funk-hecke-eigenvalue`. The eigenvalues realised in code are the per-ℓ
   Legendre moments Σ_{s,ℓ}, the diagonal of
   :class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`.
   A classical transcription, not a solver claim.
.. vv-status: sh-funk-hecke-eigenvalue documented

and those eigenvalues are exactly the **Legendre moments of the
differential scattering cross section** —
:math:`\Sigma_{s,\ell}` — which are the per-:math:`\ell` block of the
diagonal scattering operator :math:`\Lambda` =
:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`. The
:math:`m`-independence of the eigenvalue is forced by **Schur's
lemma**: the scattering operator commutes with every rotation, so on
each :math:`SO(3)`-irreducible block
:math:`V_\ell = \mathrm{span}\{Y_\ell^m\}` (dimension :math:`2\ell+1`)
it must act as a scalar. The block dimension :math:`2\ell+1` is the
origin of the :math:`(2\ell+1)` reconstruction factor
(:eq:`sh-addition-theorem-reconstruction`); the addition theorem
:eq:`real-sh-addition-theorem` is the *spectral resolution* of the
zonal kernel — the rank-:math:`(2\ell+1)` projector onto the
degree-:math:`\ell` eigenspace.

Consequently the Pℓ scattering kernel
:math:`Q^{(\ell\ge1)} = R\,\Lambda\,M` is the **spectral theorem**
:math:`A = U\Sigma U^*` written out: :math:`M` (analysis) is the
change of basis into the eigenbasis :math:`U^*`, :math:`\Lambda` is
the diagonal spectrum :math:`\Sigma`, and :math:`R` (reconstruction)
is the synthesis :math:`U`. The streaming operator
:math:`\hat\Omega\cdot\nabla`, by contrast, carries the
:math:`\ell=1` direction irrep and does **not** commute with rotations
— it couples :math:`Y_\ell^m` to :math:`Y_{\ell\pm1}^m` (the
block-tridiagonal Pℓ recurrence), so it is *not* diagonalised by this
basis. The harmonics are chosen to diagonalise collision; streaming is
merely tolerated. That asymmetry is what assigns the frame's ownership
to the scattering operator.

.. note::

   The full derivation — Funk–Hecke + Schur, the
   :math:`U\Sigma U^*` reading of :math:`R\,\Lambda\,M`, the
   streaming Clebsch–Gordan asymmetry, the literature corroboration,
   and the **unifying principle** *"an operator owns its frame iff
   the frame is its eigenbasis"* that explains why angular scattering
   is Galerkin while energy condensation and spatial homogenisation
   are Petrov-Galerkin — lives in
   :ref:`frame-eigenbasis-ownership` (:doc:`/theory/foundations/frame`).
   The relocation tripwire (when a second consumer with an
   :math:`L` independent of ``scattering_order`` moves the
   constructor ownership off the scattering operator onto the neutral
   :meth:`Quadrature.angular_frame(L)
   <orpheus.numerics.quadrature.Quadrature.angular_frame>` factory)
   is documented at :ref:`frame-eigenbasis-relocation-tripwire`.


Definitions
===========

For an ordinate :math:`\hat\Omega_n = (\mu_{x,n}, \mu_{y,n}, \mu_{z,n})`
on the unit sphere, with polar axis :math:`\mu_x`, the real
spherical harmonics under the no-prefactor convention are:

.. math::
   :label: real-sh-l0

   Y_0^0(\hat\Omega) \;=\; 1.

.. math::
   :label: real-sh-l1

   Y_1^{-1}(\hat\Omega) \;=\; \mu_z, \qquad
   Y_1^{0}(\hat\Omega)  \;=\; \mu_x, \qquad
   Y_1^{+1}(\hat\Omega) \;=\; \mu_y.

For :math:`\ell \ge 2`, with :math:`P_\ell^m(\cdot)` the unnormalised
associated Legendre function (Condon–Shortley phase removed):

.. math::
   :label: real-sh-l2plus

   Y_\ell^0(\hat\Omega) &\;=\; P_\ell(\mu_x), \\
   Y_\ell^{m}(\hat\Omega) &\;=\;
     \sqrt{\tfrac{2(\ell-m)!}{(\ell+m)!}}\,P_\ell^{m}(\mu_x)\,
     \cos(m\phi),\quad m > 0, \\
   Y_\ell^{-m}(\hat\Omega) &\;=\;
     \sqrt{\tfrac{2(\ell-m)!}{(\ell+m)!}}\,P_\ell^{m}(\mu_x)\,
     \sin(m\phi),\quad m > 0.

.. vv-status: real-sh-l0 documented
.. vv-status: real-sh-l1 documented
.. vv-status: real-sh-l2plus documented


Discrete orthogonality on a quadrature
======================================

For a discrete angular cubature
:math:`\mu_{S^2} = \sum_n w_n\,\delta_{\hat\Omega_n}` whose
``degree_of_exactness`` is at least :math:`2L`, the no-prefactor
real :math:`Y_\ell^m` satisfy the **discrete** orthogonality

.. math::
   :label: real-sh-discrete-orthogonality

   \sum_{n=1}^{N} w_n\,
   Y_\ell^m(\hat\Omega_n)\,Y_{\ell'}^{m'}(\hat\Omega_n)
   \;=\; \frac{4\pi}{2\ell+1}\,
         \delta_{\ell\ell'}\,\delta_{mm'},
   \qquad \ell + \ell' \le 2L.



.. implements:: real-sh-discrete-orthogonality
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.evaluate

   **Implemented by** 4 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: real-sh-discrete-orthogonality
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.mass_matrix

.. implements:: real-sh-discrete-orthogonality
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.metric_per_ell

.. implements:: real-sh-discrete-orthogonality
   :by: orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace

This identity is the discretised form of the continuous orthogonality
on :math:`L^2(S^2)`. Combined with the addition theorem
:eq:`real-sh-addition-theorem`, it produces the central numerical
contract the Galerkin projection / reconstruction pair satisfies on a
sufficiently-exact angular cubature:

.. math::
   :label: pi-r-equals-4pi-i

   \Pi \, R \;=\; 4\pi \, I_{\text{coefficient space}},


.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.addition_theorem_factor

   **Implemented by** 7 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.analyze

.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.reconstruct

.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.frame.FrameBase.gram

.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.frame.GalerkinFrame

.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.frame._FrameAnalysis

.. implements:: pi-r-equals-4pi-i
   :by: orpheus.numerics.frame._FrameReconstruction

where :math:`\Pi` is the spherical-harmonic frame's **analysis
face** (``frame.analysis``, :math:`M = Y^*W`), :math:`R` is its
**reconstruction face** (``frame.reconstruction``, with the
:math:`(2\ell+1)` factor), and the :math:`4\pi` factor comes from
the no-prefactor convention summing the :math:`4\pi/(2\ell+1)`
orthogonality with the :math:`(2\ell+1)` reconstruction weight —
i.e. the frame is **4π-tight** (frame operator :math:`S = T^*T =
4\pi I`). The identity is verified at :math:`L=2,\,3,\,4` against
Lebedev quadratures of order :math:`7,\,13,\,17` by the L1 test
``tests/numerics/test_spherical_harmonic_space.py``.


The four operators ERR-039 originally conflated
================================================

Post-P1.4 of the moment-space + layering plan, the four operators
that share the ``(SH coefficient → angular ordinate)`` signature
are SEPARATELY TYPED with mathematically distinct semantics. Each is
the **naked synthesis** :math:`S_0(c)_n = \sum_{\ell, m}
Y_\ell^m(\hat\Omega_n) c_\ell^m` post-multiplied by a diagonal that
lives in exactly ONE place in the codebase.

.. important::

   **An adjoint is METRIC-RELATIVE.** :math:`\Pi^*` is defined by
   :math:`\langle \Pi\psi, c\rangle_{\rm coeff} =
   \langle \psi, \Pi^* c\rangle_W` — a *pair* of inner products, so
   naming "the" Hilbert adjoint without naming the **coefficient-space
   metric** says nothing. All three metrics below have appeared in this
   corpus, all three are internally consistent, and they induce three
   different :math:`\Pi^*`. The one the code exposes is the third.

   .. list-table:: The three coefficient-space metrics and the :math:`\Pi^*` each induces
      :header-rows: 1
      :widths: 26 34 40

      * - Coefficient metric
        - Where it lives
        - The adjoint it induces
      * - **Euclidean** (no metric)
        - The bare-transpose reading; the frame never installs it on
          the SH codomain.
        - :math:`\Pi^* = S_0` — the naked synthesis, no factor.
      * - **Continuum Gram** :math:`g_C = 4\pi/(2\ell+1)`
        - :meth:`SphericalHarmonicSpace.from_L
          <orpheus.numerics.spaces.SphericalHarmonicSpace.from_L>`
          (:eq:`sh-space-metric`) — the ``project``/``gram``
          cross-Gram vocabulary.
        - :math:`\Pi^* = g_C \cdot S_0`. ⛔ **What the frame exposed
          before F-0** — the wrong side for covariant moments.
      * - **Parseval metric** :math:`G^{-1}`, the inverse *discrete*
          Gram — :math:`(2\ell+1)/4\pi` on a degree-exact sphere rule
        - :attr:`FrameBase.basis_space
          <orpheus.numerics.frame.FrameBase.basis_space>` — the frame
          dresses the basis's space with it (F-0; on a rule whose
          discrete Gram measures ``DENSE`` the dressing is the matrix
          pseudo-inverse :math:`G^{+}` instead — P7,
          :ref:`frame-parseval-dense-arm`).
        - :math:`\Pi^* = S_0 \circ G^{-1} = R/W`
          (:eq:`hilbert-adjoint-equals-metric-times-S0`) — **shipped**,
          and the physical adjoint for the carried moments.

   Why the third is the physical one, in one line: the analysis face's
   output is the **covariant** moment vector :math:`\varphi = Gc`
   (:ref:`frame-parseval-metric`, :doc:`/theory/foundations/frame`), and
   the inner product under which a covariant vector has the same length
   as the field it came from is the *inverse* Gram — Parseval.

* :math:`S_0` itself — the bare synthesis (the frame-theory
  synthesis operator :math:`T^*`), exposed by
  :meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.synthesize`.
* :math:`\Pi^\top = w_n \cdot S_0` — the representation transpose of
  the analysis face, exposed by
  :meth:`frame.analysis.apply_transpose
  <orpheus.numerics.basis.SphericalHarmonicBasis.analyze_transpose>`.
  The :math:`w_n` factor is the quadrature weight carried on the
  analysis face's domain (the frame's ``measure_space``)
  ``inner_product_weights``.

.. math::
   :label: moment-projection-transpose-T

   (\Pi^\top c)_n
   \;=\; w_n \, S_0(c)_n
   \;=\; w_n \sum_{\ell, m} Y_\ell^m(\hat\Omega_n) c_\ell^m.

* :math:`\Pi^* = S_0 \circ G^{-1}` — the **Hilbert adjoint under the
  frame's Parseval metric**, exposed by ``frame.analysis.H`` and
  computed generically by the metric-aware ``AdjointOperator``
  wrapper. The :math:`G^{-1}` factor is the inverse of the frame's
  **discrete** trial Gram
  (:attr:`FrameBase.discrete_gram
  <orpheus.numerics.frame.FrameBase.discrete_gram>`), carried on the
  analysis face's codomain (the frame's ``basis_space``) as its
  ``inner_product_weights``. That codomain is a
  :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace` — but the
  *frame-dressed* copy, not the one
  :meth:`~orpheus.numerics.spaces.SphericalHarmonicSpace.from_L`
  builds: ``from_L`` carries the CONTINUUM Gram
  :math:`g_C` (:eq:`sh-space-metric`) and
  :attr:`FrameBase.basis_space
  <orpheus.numerics.frame.FrameBase.basis_space>` REPLACES it with
  :math:`G^{-1}`. The general statement, and the SH collapse of it:

.. math::
   :label: hilbert-adjoint-equals-metric-times-S0

   (\Pi^* c)_n
   \;=\; \bigl(S_0 \, G^{-1} c\bigr)_n
   \;=\; \sum_\ell \frac{2\ell+1}{4\pi} \sum_m
              Y_\ell^m(\hat\Omega_n)\, c_\ell^m
   \;=\; \frac{(R\,c)_n}{W},
   \qquad W \;=\; \sum_n w_n \;=\; 4\pi .

The **first** equality is the general frame law — true for every
frame whose discrete Gram is diagonal, whatever its values. The
**second** substitutes the SH discrete Gram
:math:`G_\ell = 4\pi/(2\ell+1)`, which a sphere cubature of
``degree_of_exactness`` :math:`\ge 2L` realises exactly
(:eq:`real-sh-discrete-orthogonality`). The **third** — the frame
square closing on the single scalar :math:`W` — additionally needs
the per-:math:`\ell` identity :math:`d_\ell G_\ell = W` with
:math:`d_\ell = 2\ell+1` the addition-theorem factor; it is a
property of *this* basis-measure pairing, not of frames in general
(the indicator frame satisfies Parseval and does **not** satisfy
:math:`M^* = R/W`; see :ref:`frame-parseval-metric`).

.. note::

   ⛔ **Corrected 2026-08-23 (step F-0,**
   ``.claude/plans/frame_square_recarve.md`` **).** Until F-0 this
   equation read

   .. math::

      (\Pi^* c)_n \;=\; \sum_\ell \frac{4\pi}{2\ell+1}
        \sum_m Y_\ell^m(\hat\Omega_n)\, c_\ell^m
        \qquad\text{(pre-F-0 — the CONTINUUM metric } g_C \text{)},

   because the frame exposed ``basis.space`` unchanged and that space
   carries :math:`g_C`. The machinery was self-consistent throughout —
   the adjoint identity :math:`\langle\Pi\psi,c\rangle_{g_C} =
   \langle\psi,\Pi^*c\rangle_W` held at the round-off floor
   (`[M]` 2026-08-23, LS\ :sub:`4`: relative residual
   :math:`9.5\times10^{-16}` at :math:`L=1`, **exactly** :math:`0.0`
   at :math:`L=2`) — so
   nothing could fail; the defect was *which metric was stored*, and
   :math:`g_C` is the **wrong side** for the covariant moments the
   analysis face carries. The pre-F-0 :math:`\Pi^*` is off the
   physical one by exactly :math:`(4\pi/(2\ell+1))^2` per
   :math:`\ell` — :math:`157.9` at :math:`\ell=0`, :math:`17.5` at
   :math:`\ell=1`, :math:`6.3` at :math:`\ell=2`. Measured
   consequence: Parseval read a *ratio* of :math:`81.4` (:math:`L=1`)
   and :math:`65.2` (:math:`L=2`) on an LS\ :sub:`4` rule instead of
   :math:`1.000\,000\,000\,0`. (The ratio is a moment-energy-weighted
   average of those per-:math:`\ell` factors, so its value depends on
   the coefficient draw — what is seed-independent is that it lies
   between the extreme factors PRESENT AT THAT :math:`L`
   (:math:`[17.5,\,157.9]` at :math:`L=1`, :math:`[6.3,\,157.9]` at
   :math:`L=2`) and can therefore never be 1.) The equation label
   survives the correction because it is an API: live
   ``@pytest.mark.verifies`` markers point at it, and every
   :eq:`hilbert-adjoint-equals-metric-times-S0` citer would otherwise
   inherit the retired claim.


.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.analyze_transpose

   **Implemented by** 8 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others. **Re-derived at F-0 (2026-08-23):**
   ``SphericalHarmonicBasis.metric_per_ell`` LEFT this set — it
   produces the CONTINUUM Gram :math:`g_C`, which after F-0 is no
   longer the metric this equation reads; it is declared against
   :eq:`sh-space-metric` instead. The two frame properties that now
   produce the Parseval metric — ``FrameBase.discrete_gram`` (which
   computes :math:`G`) and ``FrameBase.basis_space`` (which installs
   :math:`G^{-1}` as the codomain metric) — JOINED it.

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.synthesize

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.frame.FrameBase.discrete_gram

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.frame.FrameBase.basis_space

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.frame._FrameAnalysis

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.operator.AdjointOperator

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.operator.AdjointOperator.apply

.. implements:: hilbert-adjoint-equals-metric-times-S0
   :by: orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace

* :math:`R = (2\ell+1) \cdot S_0 = W \cdot G^{-1} \cdot S_0
  = W \cdot \Pi^*` — the addition-theorem reconstruction, exposed by the frame's
  **reconstruction face** (``frame.reconstruction``,
  :meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.reconstruct`),
  which reads the :math:`(2\ell+1)` factor live from
  :attr:`SphericalHarmonicBasis.addition_theorem_factor
  <orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor>`
  (re-exposed on the coefficient space as
  :attr:`SphericalHarmonicSpace.addition_theorem_factor`).

.. math::
   :label: sh-addition-theorem-reconstruction

   (R \cdot c)_n
   \;=\; \sum_\ell (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n) c_\ell^m.

Reading that middle equality right-to-left is the whole content of
the frame square: :math:`R` and :math:`\Pi^*` are the **same
operator up to the single scalar** :math:`W = \sum_n w_n`, because
the addition-theorem factor :math:`d_\ell = 2\ell+1` and the
discrete Gram :math:`G_\ell = 4\pi/(2\ell+1)` multiply to
:math:`d_\ell G_\ell = 4\pi = W` for **every** :math:`\ell`. The
:math:`1/W` prefactor the scattering operator applies once
(:eq:`scattering-aniso-composite`,
:doc:`/theory/foundations/operator_algebra`) IS that scalar.

The continuum metric :math:`g_C`
---------------------------------

The metric :math:`g_C` is the single source of truth for the SH
**convention**; it is the continuum Gram of the no-prefactor
harmonics, it lives on :class:`SphericalHarmonicSpace` as built by
:meth:`~orpheus.numerics.spaces.SphericalHarmonicSpace.from_L`, and
it equals :math:`\mathrm{diag}(4\pi/(2\ell+1))` per :math:`\ell`:

.. math::
   :label: sh-space-metric

   \langle c, d \rangle_C
   \;=\; \sum_\ell \frac{4\pi}{2\ell+1} \sum_m c_\ell^m d_\ell^m.


.. implements:: sh-space-metric
   :by: orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.metric_per_ell

   **Implemented by** 3 sites — the formula
   (``SphericalHarmonicBasis.metric_per_ell``), its broadcast into the
   padded ``(L+1, 2L+1)`` storage layout (``_padded_metric_tensor``),
   and the constructor that installs it on the space
   (``SphericalHarmonicSpace.from_L``). Declared at F-0 (2026-08-23),
   when ``metric_per_ell`` left
   :eq:`hilbert-adjoint-equals-metric-times-S0`: the continuum Gram is
   this equation's subject, not the adjoint's factor.

.. implements:: sh-space-metric
   :by: orpheus.numerics.spaces.spherical_harmonic_space._padded_metric_tensor

.. implements:: sh-space-metric
   :by: orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace.from_L

.. warning::

   :eq:`sh-space-metric` is the **continuum** Gram, and after F-0 it is
   *not* the metric a frame's coefficient codomain carries. The frame
   REPLACES it: :attr:`FrameBase.basis_space
   <orpheus.numerics.frame.FrameBase.basis_space>` dresses the space
   with the inverse of the **discrete** Gram measured on that frame's
   own measure — a diagonal ``inner_product_weights`` array when that
   Gram measures ``DIAGONAL``, and since campaign 1 P7 (2026-08-30) the
   matrix pseudo-inverse :math:`G^{+}` as a
   :class:`~orpheus.numerics.metric.DenseMetric` when it does not. On a
   degree-exact sphere cubature the two are reciprocal
   (:math:`G = g_C`, so :math:`G^{-1} = g_C^{-1}`) and the distinction
   is invisible in the *values*; on the slab Gauss–Legendre measure the
   discrete Gram is not even diagonal, so **no diagonal metric
   reproduces** :math:`g_C^{-1}`\ **'s role at all** — which is why
   that frame needs the matrix realization rather than a better choice
   of diagonal (`[M]` the best diagonal candidate reads a Parseval
   ratio of :math:`1.806`; see :ref:`frame-parseval-dense-arm`). Use
   :eq:`sh-space-metric` when you mean the convention (``project`` /
   :attr:`FrameBase.gram <orpheus.numerics.frame.FrameBase.gram>`, the
   cross-Gram :math:`MR`); use
   :eq:`hilbert-adjoint-equals-metric-times-S0` when you mean what
   ``.H`` computes.

.. warning::

   ⛔ **The paragraph above is true and its diagnosis is INCOMPLETE, and
   the missing half is the load-bearing one (ERR-080, 2026-08-31).** On
   the slab Gauss–Legendre measure the discrete Gram is not merely
   *non-diagonal* — at :math:`L = 2` it is **rank-deficient**: 5 live
   slots, rank 4. So :math:`G^{+}` is not "the matrix realization that a
   non-diagonal Gram needs"; it is a pseudo-inverse **silently
   discarding a null direction**, which is a structural degeneracy
   wearing a conditioning costume.

   The null direction has a closed form. The offending harmonic is
   :math:`Y_2^{+2} \propto (1-\mu^2)`, which is exactly
   :math:`\det P = 4(1-\mu^2)`, the squared orbit radius of the
   :math:`SO(2)_x` action on :math:`S^2` (that :math:`\det P` is derived
   from the Procesi–Schwarz condition at :eq:`manifold-s2-mod-so2`, and
   what its three appearances do and do not share is
   :ref:`manifold-one-polynomial`); the predicted null vector
   :math:`[-0.447214,\,0,\,+0.447214,\,0,\,+0.774597]` over the live
   slots :math:`\{(0,0),(1,0),(2,0),(2,1),(2,2)\}` is the SVD-measured
   one: ``[M]`` alignment :math:`\bigl|1 - |\cos\theta|\bigr| =`
   2.75e-14 (the *component-wise* difference reads 1.5e-07 and measures
   only the six digits the prediction was written to — quote the
   alignment, not the difference). **It is a theorem about the quotient,
   not roundoff.**

   The cause is upstream of the metric entirely: a 1-D rule carries no
   azimuthal information, but the measure ``Quadrature.angular_frame``
   integrates against pads the polar nodes with two zero columns and
   then declares the result ``support=SPHERE`` over nodes with
   :math:`\lVert\Omega\rVert \ne 1`. ``_evaluate_real_sh`` duly reads
   :math:`\arctan2(0, 0) = 0` and every :math:`m > 0` harmonic becomes a
   non-zero constant across the ordinate set. The frame is therefore
   **ill-posed on that pairing**, and no choice of metric — diagonal,
   dense, or otherwise — repairs a basis that is linearly dependent on
   its own nodes.

   .. note::

      ⛔ **This paragraph read "**\ ``angular_frame`` **builds its
      measure by** ``column_stack``\ **ing three axis-cosine arrays —
      two of which are the zero FALLBACK, not data" until 2026-09-02.**
      True when written and repealed by #429's own phases 0.1a and 0.2
      (both 2026-09-01), which is the campaign correcting the
      *description* while leaving the *defect*: the padding now runs
      only on the 1-D arm of ``Quadrature._harmonic_frame_measure``
      (`[M]` 10 of 12 shipped rules hand the frame their own measure
      instead), the zeros are written at that one site rather than
      obtained from an accessor named *"direction cosine along axis
      i"*, and ``axis_cosines`` now **refuses** a suppressed axis. None
      of that changes the flux: the fabrication still happens on every
      slab solve.

      ⭐ Since 2026-09-02 (tracker 2.3) the tree also has the
      vocabulary to say what the padding *is*. The map
      :math:`\mu \mapsto \mu\,\hat e_a` is the **orbit barycentre**,
      typed as
      :func:`~orpheus.numerics.manifold.barycentre` with codomain
      :class:`~orpheus.numerics.manifold.Ball`\ ``(3)`` — because
      :math:`1 - \lVert\mu\hat e_a\rVert^2 = 1-\mu^2 = \tfrac14\det P`,
      the same squared orbit radius as above. So ERR-080's level-1 half
      is *that map with a forged codomain*, and `[M]` the forged nodes
      are ``np.array_equal`` to the honest map's image
      (:ref:`manifold-barycentre`). The arithmetic was never wrong; a
      **type** is.

   ⚠ Consequently, campaign 1 P7's reading of the same measurement — one
   ``~1e-16`` *"noise mode"*, recorded in
   :mod:`orpheus.numerics.metric` and used to pin ``_DENSE_METRIC_RCOND``
   — was a correct measurement with an inverted interpretation. The
   numbers in the paragraph above (including the :math:`1.806` diagonal
   candidate) all stand; what does not stand is the inference that a
   dense metric is the whole answer.

   The repair is not a special case: a 1-D angular quadrature is a
   quadrature on the orbit space :math:`S^2/SO(2)_x`, and the surviving
   harmonics are that quotient's **trivial isotypic component**
   :math:`\{Y_\ell^0\} \cong \{P_\ell\}`. Tracked by **#429**;
   planned in
   ``.claude/plans/angular_spaces_derived_from_symmetry.md``. Until it
   lands, :math:`P_{\ge 2}` scattering on any 1-D chart returns a wrong
   answer — gated by
   ``tests/sn/solve/test_pl_order_does_not_move_the_infinite_medium_flux.py``.

   ⭐ **The axis in** :math:`SO(2)_x` **is load-bearing, and naming it
   landed 2026-09-01** (tracker 2.4). It is :math:`x` for the same
   reason this defect exists at all: `[M]` ``_evaluate_real_sh`` takes
   ``cos θ = mu_x``, so the real spherical-harmonic pole is :math:`x`
   and :math:`\{P_\ell(\mu_x)\}` are that basis's :math:`m = 0` members
   — while every product rule's polar factor is :math:`\mu_z`. One
   Gauss–Legendre rule serves both roles, which is why the group cannot
   be spelled without its axis
   (:ref:`manifold-so2-axis-is-a-parameter`). `[M]` the slab's rule now
   declares ``support = S^2/SO2_x``; ⛔ that names the space the repair
   must be posed on and **does not perform the repair** — the frame's
   own measure is still the forged one.

These identities are pinned by
``tests/numerics/test_spherical_harmonic_space.py`` and
``tests/numerics/test_frame.py`` (the
``@pytest.mark.catches("ERR-039")`` suites; the Parseval-metric arm is
the ``test_parseval_*`` family, including the loaded-not-blind
negative leg that re-installs the pre-F-0 continuum metric in-process
and measures the ratio it produces).

.. (vv-status rationale) Both labels are face-DISTINCTION identities:
   they say which diagonal each of the four same-signature operators
   carries (the representation transpose carries w_n; the Hilbert
   adjoint carries the codomain metric), not what any solver computes.
   Each is nevertheless pinned by a live L1 gate against an independent
   closed-form einsum — ``test_T_carries_w_n_and_H_carries_the_parseval_metric``
   and ``test_H_equals_parseval_metric_times_S0`` in
   ``tests/numerics/test_spherical_harmonic_space.py``, plus
   ``test_parseval_frame_square_closes`` (6 sphere families) and
   ``test_analysis_hilbert_adjoint_falls_out_of_the_frame_spaces`` in
   ``tests/numerics/test_frame.py``.
.. vv-status: hilbert-adjoint-equals-metric-times-S0 documented
.. vv-status: moment-projection-transpose-T documented

.. note::

   ERR-039's original confusion was the result of two missing
   abstractions: the :math:`(2\ell+1)` literal lived inline on the
   reconstruction operator with no typed home, and the projection's
   ``apply_transpose`` returned the bare :math:`S_0` but its
   docstring labeled it the W-weighted Hilbert adjoint. The endpoint
   fix gives each of the four operators a distinct construction path
   with the metric / weight diagonals carried by typed spaces; the
   Frame/Basis carve then re-homed the projection and reconstruction
   as the :class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` /
   ``reconstruction`` faces and moved the :math:`(2\ell+1)` factor
   onto
   :attr:`SphericalHarmonicBasis.addition_theorem_factor
   <orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor>`
   (one home). The composition :math:`\Pi R = 4\pi I`
   (the addition-theorem composition, :eq:`pi-r-equals-4pi-i`)
   continues to hold on band-limited inputs and is the genuine
   Galerkin-discipline identity for this SH frame — its 4π-tightness.

   **The F-0 chapter (2026-08-23).** ERR-039's endpoint gave each
   operator a typed construction path and left one question
   unasked: *is the metric the typed space carries the RIGHT one?*
   It was not. The space carried the continuum Gram :math:`g_C`,
   which pairs **contravariant** coefficients :math:`c`; the
   analysis face emits **covariant** moments :math:`\varphi = Gc`,
   whose metric is :math:`G^{-1}`. Every gate stayed green because
   every gate checked *consistency* (the sandwich reproduces the
   pairing it was built from) rather than *Parseval* (does analysis
   preserve length?). This is the same family as ERR-039 — a metric
   / transpose / adjoint conflation — one level deeper: right Gram,
   **wrong side**. The catching instrument is the isometry
   :math:`\|M\psi\|_{\rm codomain} = \|\psi\|_W`, which no
   consistency identity implies. Full derivation:
   :ref:`frame-parseval-metric` (:doc:`/theory/foundations/frame`).

The numerical evidence:

.. list-table:: Galerkin idempotency residuals
   :header-rows: 1
   :widths: 18 18 18 24 22

   * - Lebedev order
     - :math:`L`
     - :math:`N` ordinates
     - Residual on
       :math:`\| \Pi R c - 4\pi c \|_\infty`
     - Convergence floor
   * - 7
     - 2
     - 26
     - :math:`\le 10^{-12}`
     - quadrature-exact for :math:`\ell \le 4`
   * - 13
     - 3
     - 74
     - :math:`\le 10^{-12}`
     - quadrature-exact for :math:`\ell \le 6`
   * - 17
     - 4
     - 110
     - :math:`\le 10^{-12}`
     - quadrature-exact for :math:`\ell \le 8`

The :math:`10^{-12}` floor is the test's tolerance, not the actual
floating-point limit; the agreement is at machine precision (≤
``nulp ≈ 16`` on the multiplications).


Implementation map
==================

The pure-Python reference implementation is
:meth:`SphericalHarmonicBasis.evaluate
<orpheus.numerics.basis.SphericalHarmonicBasis.evaluate>` (and its
per-component sibling
:meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.evaluate_from_components`).
Its return value is the ``(N, L+1, 2L+1)`` **table** consumed by the
spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`. The frame
caches the table once (``frame.table``) and the two faces delegate:

* the **analysis face** ``frame.analysis`` — the Galerkin projection
  :math:`\phi^{\ell m} = \sum_n w_n\,Y_\ell^m(\hat\Omega_n)\,\psi_n`
  (:meth:`SphericalHarmonicBasis.analyze
  <orpheus.numerics.basis.SphericalHarmonicBasis.analyze>`).
* the **reconstruction face** ``frame.reconstruction`` — the
  addition-theorem reconstruction
  :math:`q_n = \sum_\ell (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n)\,
  \phi^{\ell m}`
  (:meth:`SphericalHarmonicBasis.reconstruct
  <orpheus.numerics.basis.SphericalHarmonicBasis.reconstruct>`).

The single home of the :math:`S^2` embedding is
:meth:`Quadrature.angular_frame(L)
<orpheus.numerics.quadrature.Quadrature.angular_frame>`, which binds
the SH basis to the quadrature's angular measure so PN, SN, and MC
consume the same frame without importing from :mod:`orpheus.sn`.

The complete data flow:

.. code-block:: text

   Quadrature (Lebedev / level-symmetric / Gauss-Legendre × φ)
        │   .angular_frame(L)
        ▼
   Frame(SphericalHarmonicBasis(L), s2_measure)
        │
        │   .table = basis.evaluate(measure.nodes) → Y[n, ℓ, ℓ+m]
        │   (cached once; both faces delegate)
        │
        ├──────────────► frame.analysis        (M = Y* W)
        │                    Π : ψ_n  →  φ^{ℓm}
        │
        └──────────────► frame.reconstruction  (R = (2ℓ+1) S₀)
                             R : φ^{ℓm}  →  q_n


Cross-method consumers
======================

The same :math:`Y_\ell^m` table is consumed by multiple solvers; this
is the architectural reason the basis lives in
:mod:`orpheus.numerics.basis`, not :mod:`orpheus.sn`.

.. list-table:: Cross-method consumption of SphericalHarmonicBasis
   :header-rows: 1
   :widths: 22 38 40

   * - Solver
     - How it consumes :math:`Y_\ell^m`
     - Status
   * - SN aniso scattering
     - :math:`Q^{\ell\ge 1}_n = R \Lambda M\,\psi` builds the
       per-ordinate Pℓ source via the frame's analysis /
       reconstruction faces. See
       :class:`~orpheus.transport.operators.scattering.ScatteringOperator`.
     - Live (Frame/Basis carve; the SN scattering operator pulls its
       frame from
       :meth:`Quadrature.angular_frame(L)
       <orpheus.numerics.quadrature.Quadrature.angular_frame>`).
   * - PN solver (§10)
     - Native moment-space basis on the streaming-coupling.
     - Pending (PN solver not yet implemented; the basis is
       ready when it lands).
   * - MC adjoint moments
     - Variance reduction with response moments built against
       :math:`Y_\ell^m`.
     - Pending.
   * - Energy-condensation diagnostics (§17)
     - Within-group anisotropy characterisation.
     - Pending.

The user's architectural rule "**unify after two instances**"
(`feedback_unify_after_two_instances.md`) does not apply here:
the spherical-harmonic evaluator is shared upstream infrastructure,
not a method-specific algorithm. The single SN consumer in production
today is sufficient justification because PN is a queued consumer
with a near-identical use, not a hypothetical one.


Code references
===============

* :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` — the
  basis carrying the :math:`Y_\ell^m` evaluator
  (:meth:`~orpheus.numerics.basis.SphericalHarmonicBasis.evaluate`),
  the no-prefactor convention, and the
  :attr:`~orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor`
  :math:`(2\ell+1)`.
* :class:`~orpheus.numerics.frame.GalerkinFrame` — binds the basis to the
  angular measure; its ``analysis`` face is the Galerkin projection
  :math:`M` that consumes the :math:`Y` table, its ``reconstruction``
  face the addition-theorem :math:`R`.
* :meth:`Quadrature.angular_frame(L)
  <orpheus.numerics.quadrature.Quadrature.angular_frame>` — the
  single home of the :math:`S^2` embedding; builds the SH frame on
  a quadrature.

The pedagogical companion at
``student_resources/02_spherical_harmonics.py`` is a single-:math:`Y_\ell^m`
surface-plot visualisation. It shares the same
:func:`scipy.special.lpmv` machinery and norm
:math:`\sqrt{2(\ell-|m|)!/(\ell+|m|)!}`; do not duplicate the
evaluator there.


History — what was tried and discarded
======================================

The SN solver originally carried
``orpheus.sn.quadrature._build_spherical_harmonics`` — a private
module-level function on the SN quadrature adapter that knew about
the SN ordinate axis layout. Two issues drove the lift to
:mod:`orpheus.numerics`:

1. **Convention drift risk.** The SN-private function had no public
   docstring naming the no-prefactor convention. A future PN
   implementation would have either (a) re-implemented the evaluator
   with the standard ANSI convention, breaking the addition-theorem
   identity that SN relies on, or (b) imported the SN-private
   function and inherited a hidden module dependency. Both paths
   trigger failure mode 6 (definition-site / usage-site convention
   drift).

2. **Cross-method reuse.** The PN solver, MC adjoints, and
   energy-condensation diagnostics all need the same table. Keeping
   it in :mod:`orpheus.sn` would force `import orpheus.sn` from
   modules that have no business depending on SN.

The lift was a pure rename + module relocation; the implementation
is bit-identical to the legacy code (regression snapshots at
``tests/sn/regression/snapshots/`` survive unchanged).


References
==========

* Bell, G. I. and Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.6 (real spherical harmonics in
  transport).
* Lewis, E. E. and Miller, W. F. Jr. (1993). *Computational Methods
  of Neutron Transport*. ANS. §4.7 (the Pℓ Galerkin reconstruction
  with the :math:`(2\ell+1)` factor — the canonical form for the
  no-prefactor convention).
* Beckmann, M. and Wieselquist, W. (2017). *Numerical Recipes for
  Real Spherical Harmonics*. Comp. Phys. Comm. 220, 121–133. (Norm
  conventions and recursion stability for high :math:`\ell`.)
* A Wave 0 step C0.1 plan ("lift ``_build_spherical_harmonics`` to
  ``orpheus/numerics/``") captured the architectural rationale
  reproduced here; that plan file is no longer retained.
