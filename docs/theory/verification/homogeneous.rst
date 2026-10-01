.. _theory-verification-homogeneous:

===========================
Homogeneous Infinite Medium
===========================

.. Seeded at task #10 stage V3 from the per-solver sections of the
   dissolved ``docs/theory/verification.rst``. Prose carried verbatim
   plus one added framing pointer to the theory chapter; the
   homogeneous-solver benchmark table still lives in
   ``foundations/infinite_medium.rst`` until the V7 template-B sweep.
   V5–V7 author the part-level framing and the coverage slices.

For an infinite homogeneous medium, there is no spatial variation and no
leakage.  The neutron balance reduces to a matrix eigenvalue problem:

.. math::

   \mathbf{A} \phi = \frac{1}{k} \mathbf{F} \phi

where :math:`\mathbf{A} = \text{diag}(\Sigma_t) - \Sigma_s^T` is the
removal matrix and :math:`\mathbf{F} = \chi \otimes (\nu\Sigma_f)` is
the fission production matrix.

For 1 group: :math:`k = \nu\Sigma_f / \Sigma_a`.

For multi-group: :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`, which
for the rank-one fission dyad is
:math:`\langle\nu\Sigma_f, \mathbf{A}^{-1}\chi\rangle` (the proof:
:ref:`homogeneous-rank-one-route`).

Two references pin the solver, and they answer different questions.  The
registry's analytical values below are float64 eigenvalues of the oracle's
own fused assembly (:func:`numpy.linalg.solve`, then
:func:`numpy.linalg.eig`), compared at :math:`10^{-12}`: they catch a wrong
assembly or a wrong eigenvalue formula, and their last bits belong to the
platform's LAPACK.  The **exact rational reference**
(``orpheus/derivations/common/exact_homogeneous.py``) is the answer to the
float inputs in exact arithmetic, self-certified against
:math:`\mathbf{F}\phi = k\mathbf{A}\phi` with zero residual, and
``tests/gates/homogeneous/test_kinf_exact_reference.py`` holds
:math:`k_\infty`, the gauged flux and both condensed cross sections to it
within a forward-error bound derived per case from the LU error bound
(3.75 to 43.2 ULP in :math:`k_\infty`; `[M]` 2026-10-01 the solver reads at
most 0.86 ULP off).  The derivation of the bound, its per-case table and
what it cannot resolve are :ref:`homogeneous-exact-reference`.

The full derivation of the infinite-medium balance and the operator
treatment live in the theory chapter, :ref:`theory-homogeneous`; this
page carries the verification evidence built on it.  The analytical
eigenvalue doubles as the shared L1 anchor for every transport method:
each solver's homogeneous cases (see the sibling chapters of this
part) must reproduce these values from its own discretisation.

.. include:: ../../_generated/homogeneous_derivation.rst
