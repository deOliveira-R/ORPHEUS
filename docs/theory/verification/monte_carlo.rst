.. _theory-verification-monte-carlo:

===========
Monte Carlo
===========

.. Seeded at task #10 stage V3: the per-solver sections of the
   dissolved ``docs/theory/verification.rst`` plus the Analytical
   Verification and Verification Suite chapters of
   :doc:`/theory/methods/monte_carlo` (case descriptions moved here;
   the method chapter keeps a thin pointer).  Prose carried verbatim;
   V5–V7 author the part-level framing.

Collision-probability eigenvalue
================================

The Monte Carlo method simulates neutron random walks.  For a homogeneous
infinite medium, the expected multiplication factor can be derived from
collision probabilities:

1. Probability of fission per collision: :math:`\Sigma_f / \Sigma_t`
2. Expected secondaries per collision: :math:`\nu \Sigma_f / \Sigma_t`
3. Mean collisions before absorption: :math:`\Sigma_t / \Sigma_a`
4. Therefore: :math:`k = \nu\Sigma_f / \Sigma_a`

The MC estimator :math:`\hat{k}` converges to this with standard error
:math:`\sigma \sim 1/\sqrt{N_{\rm active}}` by the Central Limit Theorem.

.. include:: ../../_generated/mc_derivation.rst


Analytical Verification
=======================

Homogeneous Infinite Medium
---------------------------

For a single material filling the entire cell, the eigenvalue is derived
from the collision kernel.

**1-group:**  From the random walk analysis in
:doc:`/theory/methods/monte_carlo` (:eq:`fission-weight`):

.. math::
   :label: kinf-1g

   k_\infty = \frac{\nu\Sigma_f}{\Sigma_a}

For region A (1G): :math:`k = 0.75/0.5 = 1.5`.

**Multi-group:**  The eigenvalue is the dominant eigenvalue of

.. math::
   :label: kinf-mg

   k_\infty = \lambda_{\max}\!\left(\mathbf{A}^{-1} \mathbf{F}\right)

where :math:`\mathbf{A} = \text{diag}(\Sigma_t) - \Sigma_s^T` (the
transpose converts from the ``[from, to]`` storage convention to the
in-scatter form needed by the balance equation; see scattering convention
note in :ref:`theory-monte-carlo`) and
:math:`\mathbf{F} = \chi \otimes (\nu\Sigma_f)`.  Computed by
:func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous`.

.. list-table:: Homogeneous eigenvalues
   :header-rows: 1
   :widths: 30 25 25

   * - Case
     - :math:`k_{\text{ref}}`
     - Tolerance
   * - 1G (``mc_cyl1D_1eg_1rg``)
     - 1.500000
     - :math:`\sigma = 0` (exact)
   * - 2G (``mc_cyl1D_2eg_1rg``)
     - 1.875000
     - :math:`z < 5`
   * - 4G (``mc_cyl1D_4eg_1rg``)
     - 1.487762
     - :math:`z < 5`


Heterogeneous Pin Cell
----------------------

For multi-region problems, the reference eigenvalue comes from the **CP
cylinder derivation** (:func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp`).

.. list-table:: MC heterogeneous verification cases
   :header-rows: 1
   :widths: 28 8 8 18 18

   * - Case
     - G
     - Reg
     - :math:`k_{\text{ref}}`
     - Tolerance
   * - ``mc_cyl1D_1eg_2rg``
     - 1
     - 2
     - 0.989750
     - :math:`5\sigma + 0.06 k`
   * - ``mc_cyl1D_2eg_2rg``
     - 2
     - 2
     - 0.739887
     - :math:`5\sigma + 0.06 k`
   * - ``mc_cyl1D_1eg_4rg``
     - 1
     - 4
     - 0.806751
     - :math:`5\sigma + 0.06 k`
   * - ``mc_cyl1D_2eg_4rg``
     - 2
     - 4
     - 0.502725
     - :math:`5\sigma + 0.06 k`
   * - ``mc_cyl1D_4eg_2rg``
     - 4
     - 2
     - 0.648637
     - :math:`5\sigma + 0.06 k`
   * - ``mc_cyl1D_4eg_4rg``
     - 4
     - 4
     - 0.446472
     - :math:`5\sigma + 0.06 k`


Verification Suite
==================

The MC solver is verified by **55 tests** across four levels:

.. list-table:: Verification test summary
   :header-rows: 1
   :widths: 10 8 52

   * - Level
     - Count
     - Description
   * - L0
     - 31
     - Term-level isolation of each algorithmic component
   * - L1
     - 18
     - Eigenvalue: {1,2,4}G × {1,2,4}-region
   * - L2
     - 3
     - Convergence: :math:`\sigma \sim 1/\sqrt{N}`, bias, inactive cycles
   * - XV
     - 2
     - Cross-verification: MC vs CP (cylinder, slab)

**Determinism.**  Verified by ``test_mc_gaps.py::test_seed_reproducibility``
(same seed → identical results) and ``test_different_seeds_differ`` (different
seeds → different histories).


Failure Mode Coverage
---------------------

.. list-table::
   :header-rows: 1
   :widths: 25 45

   * - Failure mode
     - Tests
   * - Sign flip
     - L0-MC-004 (delta-tracking), L0-MC-016 (collision fraction)
   * - Variable swap
     - L0-MC-005, L0-MC-014 (SigS^T)
   * - Missing factor
     - L0-MC-006 (fission weight), L0-MC-015 (free-path)
   * - Factor error
     - L0-MC-012 (direction moments)
   * - Index error
     - L0-MC-005, L0-MC-007 (searchsorted)
   * - Convention drift
     - L0-MC-013 (batch stats), XS consistency
   * - 1G degeneracy
     - L1: 2G/4G high-stats, 4G fast guard
   * - Homogeneous degeneracy
     - L1: heterogeneous {2,4}-region
   * - Geometry mismatch
     - ERR-017 fixed, pitch tested
