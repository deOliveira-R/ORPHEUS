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

For multi-group: :math:`k = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`.

The full derivation of the infinite-medium balance and the operator
treatment live in the theory chapter, :ref:`theory-homogeneous`; this
page carries the verification evidence built on it.  The analytical
eigenvalue doubles as the shared L1 anchor for every transport method:
each solver's homogeneous cases (see the sibling chapters of this
part) must reproduce these values from its own discretisation.

.. include:: ../../_generated/homogeneous_derivation.rst
