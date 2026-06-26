.. _spherical-harmonics:

===========================================
Real Spherical Harmonics on a Direction Set
===========================================

Real spherical harmonics :math:`Y_\ell^m(\hat\Omega)` are the Galerkin
basis ORPHEUS uses to project an angular flux onto its Pℓ moments and
to reconstruct a per-ordinate scattering source from those moments.
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
     ``tests/sn/test_solver_components.py::TestAnisotropicScattering
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

   q(\hat\Omega) \;=\; \sum_{\ell} (2\ell+1)
     \sum_m Y_\ell^m(\hat\Omega)\,\phi^{\ell m},

with the :math:`(2\ell+1)` factor **outside** the spherical-harmonic
basis. Adopting the no-prefactor convention puts the
:math:`(2\ell+1)` factor where the equations want it and removes
:math:`4\pi/(2\ell+1)` factors from the addition-theorem identity. The
discrete projection / reconstruction pair then carries no leftover
constants.

.. warning::

   The SciPy function :func:`scipy.special.sph_harm` returns the
   **complex** spherical harmonics in the standard ANSI
   normalisation. ORPHEUS does NOT consume that function. The
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

   (S_{\rm aniso}\,\psi)(\hat\Omega)
   \;=\; \int_{4\pi}
         \Sigma_s(\hat\Omega \cdot \hat\Omega')\,\psi(\hat\Omega')\,
         d\hat\Omega',

whose kernel depends on the directions only through the cosine
:math:`\hat\Omega\cdot\hat\Omega'` — a **zonal** kernel on
:math:`S^2`. By the **Funk–Hecke theorem**, the spherical harmonics
are the eigenfunctions of any zonal kernel, with an eigenvalue that
depends on :math:`\ell` only:

.. math::

   S_{\rm aniso}\,Y_\ell^m \;=\; \Sigma_{s,\ell}\,Y_\ell^m,
   \qquad
   \Sigma_{s,\ell} \;=\; 2\pi\!\int_{-1}^{+1}
         \Sigma_s(t)\,P_\ell(t)\,dt,

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
   :ref:`frame-eigenbasis-ownership` (:doc:`galerkin_projection`).
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

.. vv-status: real-sh-discrete-orthogonality documented

This identity is the discretised form of the continuous orthogonality
on :math:`L^2(S^2)`. Combined with the addition theorem
:eq:`real-sh-addition-theorem`, it produces the central numerical
contract the Galerkin projection / reconstruction pair satisfies on a
sufficiently-exact angular cubature:

.. math::
   :label: pi-r-equals-4pi-i

   \Pi \, R \;=\; 4\pi \, I_{\text{coefficient space}},

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

.. vv-status: pi-r-equals-4pi-i documented

The four operators ERR-039 originally conflated
================================================

Post-P1.4 of the moment-space + layering plan, the four operators
that share the ``(SH coefficient → angular ordinate)`` signature
are SEPARATELY TYPED with mathematically distinct semantics. Each is
the **naked synthesis** :math:`S_0(c)_n = \sum_{\ell, m}
Y_\ell^m(\hat\Omega_n) c_\ell^m` post-multiplied by a diagonal that
lives in exactly ONE place in the codebase:

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

* :math:`\Pi^* = g_C \cdot S_0` — the W-weighted Hilbert adjoint,
  exposed by ``frame.analysis.H`` and computed generically by
  the metric-aware ``_AdjointOperator`` wrapper. The
  :math:`g_C^{\ell} = 4\pi/(2\ell+1)` factor is the SH Gram-matrix
  diagonal carried on the analysis face's codomain (the frame's
  ``basis_space``) ``inner_product_weights`` (which IS the canonical
  :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace`).

.. math::
   :label: hilbert-adjoint-equals-metric-times-S0

   (\Pi^* c)_n
   \;=\; \sum_\ell \frac{4\pi}{2\ell+1} \sum_m
              Y_\ell^m(\hat\Omega_n) c_\ell^m.

* :math:`R = (2\ell+1) \cdot S_0 = 4\pi \cdot g_C^{-1} \cdot S_0`
  — the addition-theorem reconstruction, exposed by the frame's
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

The metric :math:`g_C` is the single source of truth for the SH
convention; it lives on :class:`SphericalHarmonicSpace` and equals
:math:`\mathrm{diag}(4\pi/(2\ell+1))` per :math:`\ell`:

.. math::
   :label: sh-space-metric

   \langle c, d \rangle_C
   \;=\; \sum_\ell \frac{4\pi}{2\ell+1} \sum_m c_\ell^m d_\ell^m.

These four identities are pinned by
``tests/numerics/test_spherical_harmonic_space.py`` (the
``@pytest.mark.catches("ERR-039")`` test suite).

.. vv-status: sh-space-metric documented
.. vv-status: sh-addition-theorem-reconstruction documented
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
* `transient-giggling-cake.md`_ — Wave 0 step C0.1, "lift
  ``_build_spherical_harmonics`` to ``orpheus/numerics/``", with
  the architectural rationale captured here.

.. _transient-giggling-cake.md: ../../.claude/plans/transient-giggling-cake.md
