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
:func:`~orpheus.numerics.spherical_harmonics.evaluate_real_sh`
evaluator.

The evaluator is **generic infrastructure** — it lives in
:mod:`orpheus.numerics`, not :mod:`orpheus.sn`, because the same
:math:`Y_\ell^m` table is consumed by SN aniso scattering, by the PN
solver (when it lands; Grand Report v3 §10), and by MC adjoint moment
estimators. Lifting it out of the SN package was Wave 0 step C0.1 of
the SN performance plan (`.claude/plans/transient-giggling-cake.md`).

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
   :func:`~orpheus.numerics.spherical_harmonics.evaluate_real_sh`
   evaluator builds real :math:`Y_\ell^m` from
   :func:`scipy.special.lpmv` (associated Legendre values
   :math:`P_\ell^m(\cos\theta)`) so the convention can be controlled
   directly. Mixing the two conventions in one codebase is the
   canonical way to introduce convention-drift bugs (failure mode 6
   in the V&V skill); the project's defense is "one evaluator, one
   convention, in :mod:`orpheus.numerics.spherical_harmonics`".


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

where :math:`\Pi` is :class:`HarmonicMomentProjection`, :math:`R` is
:class:`HarmonicMomentReconstruction` (with the
:math:`(2\ell+1)` factor), and the :math:`4\pi` factor comes from
the no-prefactor convention summing the :math:`4\pi/(2\ell+1)`
orthogonality with the :math:`(2\ell+1)` reconstruction weight. The
identity is verified at :math:`L=2,\,3,\,4` against Lebedev
quadratures of order :math:`7,\,13,\,17` by the L1 test
``tests/numerics/test_projection_operators.py::
TestGalerkinIdempotencyOnLebedev::test_pi_R_is_identity_on_band_limited``.

.. vv-status: pi-r-equals-4pi-i documented

.. note::

   The strict Hilbert adjoint :math:`\Pi^*` (returned by
   :meth:`HarmonicMomentProjection.apply_transpose`) is **not** the
   same operator as :math:`R` (returned by
   :meth:`HarmonicMomentReconstruction.apply`):

   .. math::

      \Pi^* \ne R, \qquad
      (R \cdot c)_n
        \;=\; \sum_\ell (2\ell+1) \sum_m Y_\ell^m c_\ell^m, \qquad
      (\Pi^* c)_n
        \;=\; \sum_{\ell, m} Y_\ell^m c_\ell^m.

   :math:`R` carries the :math:`(2\ell+1)` factor — the
   addition-theorem weight needed by the Pℓ scattering
   reconstruction — while :math:`\Pi^*` is the naked sum (the
   :math:`W` weight of :math:`\Pi = Y^* W` is on the inner-product
   side, not the operator). Composing :math:`\Pi^* \cdot \Pi`
   yields :math:`4\pi/(2\ell+1)\,\delta_{\ell\ell'}\delta_{mm'}` —
   the Galerkin discipline's diagonal-in-:math:`\ell` adjoint
   composition — while :math:`\Pi R = 4\pi I` (the addition-theorem
   composition) is exactly identity-up-to-:math:`4\pi`. See
   :ref:`galerkin-projection` for the full discipline-level
   explanation.

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

The pure-Python reference implementation lives at
:func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh`. Its
return value is the ``(N, L+1, 2L+1)`` table consumed by:

* :class:`~orpheus.numerics.projection.HarmonicMomentProjection` —
  the Galerkin projection
  :math:`\phi^{\ell m} = \sum_n w_n\,Y_\ell^m(\hat\Omega_n)\,\psi_n`.
* :class:`~orpheus.numerics.projection.HarmonicMomentReconstruction`
  — the addition-theorem reconstruction
  :math:`q_n = \sum_\ell (2\ell+1) \sum_m Y_\ell^m(\hat\Omega_n)\,
  \phi^{\ell m}`.
* The four
  :class:`~orpheus.sn.quadrature.AngularQuadrature` adapters'
  :meth:`spherical_harmonics(L)` method — delegates to
  :func:`evaluate_real_sh` so PN can consume the same table without
  importing from :mod:`orpheus.sn`.

The complete data flow:

.. code-block:: text

   DiscreteMeasure (Lebedev / level-symmetric / Gauss-Legendre × φ)
        │
        │   .nodes  →  (N, 3) direction cosines (μ_x, μ_y, μ_z)
        │   .weights → (N,)
        ▼
   evaluate_real_sh(L, μ_x, μ_y, μ_z) → Y[n, ℓ, ℓ+m]
        │
        ├──────────────► HarmonicMomentProjection(weights, Y, L)
        │                    Π : ψ_n  →  φ^{ℓm}
        │
        └──────────────► HarmonicMomentReconstruction(Y, two_l_plus_one)
                             R : φ^{ℓm}  →  q_n


Cross-method consumers
======================

The same :math:`Y_\ell^m` table is consumed by multiple solvers; this
is the architectural reason the evaluator lives in
:mod:`orpheus.numerics`, not :mod:`orpheus.sn`.

.. list-table:: Cross-method consumption of evaluate_real_sh
   :header-rows: 1
   :widths: 22 38 40

   * - Solver
     - How it consumes :math:`Y_\ell^m`
     - Status
   * - SN aniso scattering
     - :math:`Q^{\ell\ge 1}_n = R \Lambda M\,\psi` builds the
       per-ordinate Pℓ source. See
       :class:`~orpheus.sn.scattering.ScatteringOperator`.
     - Live (Wave 0 of SN performance plan; legacy
       ``_build_spherical_harmonics`` in
       :mod:`orpheus.sn.quadrature` is now a thin re-export of
       :func:`evaluate_real_sh`).
   * - PN solver (§10)
     - Native moment-space basis on the streaming-coupling.
     - Pending (PN solver not yet implemented; the evaluator is
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

* :func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh` — the
  evaluator.
* :class:`orpheus.numerics.projection.HarmonicMomentProjection` — the
  Galerkin projection that consumes the :math:`Y` table.
* :class:`orpheus.numerics.projection.HarmonicMomentReconstruction`
  — the paired reconstruction.
* :class:`orpheus.sn.quadrature.AngularQuadrature` — the SN-side
  adapter that delegates :meth:`spherical_harmonics` to
  :func:`evaluate_real_sh`.

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
