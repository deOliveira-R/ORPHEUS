.. _theory-verification-collision-probability:

=====================
Collision Probability
=====================

.. Seeded at task #10 stage V3: the per-solver sections of the
   dissolved ``docs/theory/verification.rst`` plus the Verification
   chapter of :doc:`/theory/methods/collision_probability` (case
   descriptions moved here; the method chapter keeps a thin pointer).
   Prose carried verbatim; V5–V7 author the part-level framing.

Slab Collision Probability
==========================

The slab CP method uses the E₃ exponential integral kernel to compute
first-collision probabilities in a 1D half-cell with reflective centre and
:term:`white boundary condition` at the cell edge.

The CP matrix :math:`P_{\infty}(i,j,g)` gives the probability that a neutron
born uniformly in region *j* has its first collision in region *i*, for energy
group *g*.  It is computed from the E₃ second-difference formula and the
white-BC closure.

The eigenvalue problem dimension is (N_regions × N_groups), solved as a
dense matrix eigenvalue.  This is semi-analytical: E₃ is computed to machine
precision via ``scipy.special.expn``.

.. include:: ../../_generated/cp_slab_derivation.rst


Cylindrical Collision Probability
=================================

The cylindrical CP method uses the Ki₃/Ki₄ Bickley-Naylor kernel for a
Wigner-Seitz cell with annular regions and white boundary condition.

The Ki₄ function is the integrated Bickley-Naylor function:

.. math::

   \text{Ki}_3(x) = \int_0^{\pi/2} e^{-x/\sin\theta} \sin\theta \, d\theta, \qquad
   \text{Ki}_4(x) = \int_x^\infty \text{Ki}_3(t) \, dt

The CP matrix is computed via y-quadrature (Gauss-Legendre with breakpoints
at each radial boundary) and the Ki₄ second-difference formula.

.. include:: ../../_generated/cp_cylinder_derivation.rst


.. _nine-case-cp-grid:

The 9-Case Grid (per geometry)
==============================

Both slab and cylinder derivation modules populate a **3 × 3 grid**
of verification cases, indexed by the number of energy groups and
the number of spatial regions:

.. csv-table::
   :header: Case name, N_groups, N_regions, Layout
   :widths: 25, 12, 12, 51

   ``cp_slab_1eg_1rg``, 1, 1, A (homogeneous fissile)
   ``cp_slab_1eg_2rg``, 1, 2, A + B (fuel + moderator)
   ``cp_slab_1eg_4rg``, 1, 4, A + D + C + B (fuel + gap + clad + moderator)
   ``cp_slab_2eg_1rg``, 2, 1, A
   ``cp_slab_2eg_2rg``, 2, 2, A + B
   ``cp_slab_2eg_4rg``, 2, 4, A + D + C + B
   ``cp_slab_4eg_1rg``, 4, 1, A
   ``cp_slab_4eg_2rg``, 4, 2, A + B
   ``cp_slab_4eg_4rg``, 4, 4, A + D + C + B

The cylinder set has the same structure with ``cp_slab`` replaced by
``cp_cyl1D``. This is **9 + 9 = 18** semi-analytical cases per
Sphinx build.

The grid is deliberate, not incidental:

- **N_groups × N_regions sweep** exposes the two independent failure
  modes — matrix-eigenvalue assembly (groups) and collision-probability
  quadrature (regions). A bug in one is uncovered by the other being
  correct.
- **1/2/4 chosen at each axis** spans the meaningful regimes: 1-group
  is the degenerate :math:`k = \nu\Sigma_f/\Sigma_a` limit, 2-group is the smallest genuine multi-group
  problem, and 4-group exercises the full :math:`\chi` spectrum
  distribution across fast/resonance/thermal ranges.
- **4-region layout A + D + C + B** is the minimal non-trivial
  heterogeneous pin cell: it puts fuel (A), gap (D), clad (C), and
  moderator (B) in the physically correct radial order and forces
  every region-coupling term in the CP matrix to be non-zero.
  Reducing to 2 regions (A + B) removes the gap and cladding; to
  1 region removes heterogeneity entirely.

The per-case ``VerificationCase`` records are populated by SymPy
evaluation of the E\ :sub:`3` (slab) or Ki\ :sub:`4` (cylinder)
second-difference formulas at build time — see
``orpheus/derivations/continuous/flat_source_cp/slab.py`` and
``orpheus/derivations/continuous/flat_source_cp/cylinder.py``. Each record carries:
analytical :math:`k`, materials, geometry, matrix-eigenvalue
context, and (since PR-2) the V&V level and equation-label list
used by the test harness (see :doc:`/theory/verification/harness`).


Test suite
==========

The CP test tree (``tests/cp/``) verifies the implementation against
analytical eigenvalues computed independently by the derivation
modules; the auto-generated :doc:`matrix` carries the current
per-module test counts.

Consolidated Eigenvalue Solvers (CP-20260405-007)
-------------------------------------------------

:func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp` and
:func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous` replaced 10 duplicated
eigenvalue computations across ``homogeneous.py``, ``sn.py``, ``moc.py``,
``mc.py``, ``cp_slab.py``, ``cp_cylinder.py``, and ``cp_sphere.py``.
Both accept optional ``sig_2`` / ``sig_2_mats`` for (n,2n) reactions.
When omitted, (n,2n) is zero and all previous eigenvalues are preserved.


Eigenvalue Verification Cases
-----------------------------

27 core cases: {1, 2, 4} groups × {1, 2, 4} regions ×
{slab, cylinder, sphere}.

.. list-table::
   :header-rows: 1
   :widths: 25 15 20 20

   * - Geometry
     - Tolerance
     - Test file
     - Derivation
   * - Slab (:math:`E_3`)
     - :math:`< 10^{-6}`
     - ``tests/cp/test_slab.py``
     - ``orpheus/derivations/continuous/flat_source_cp/slab.py``
   * - Cylinder (:math:`\text{Ki}_4`)
     - :math:`< 10^{-5}`
     - ``tests/cp/test_cylinder.py``
     - ``orpheus/derivations/continuous/flat_source_cp/cylinder.py``
   * - Sphere (:math:`e^{-\tau}`)
     - :math:`< 10^{-5}`
     - ``tests/cp/test_sphere.py``
     - ``orpheus/derivations/continuous/flat_source_cp/sphere.py``

Historically, the cylinder/sphere tolerances were 10× looser
because the legacy 20 000-point :math:`\text{Ki}_4` table
interpolation introduced :math:`O(\Delta x^2) \approx
6\times 10^{-6}` error. Phase B.4 (commit ``6badbe5``, Issue #94)
replaced the table with the Chebyshev interpolant
:func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` of canonical
:math:`\mathrm{Ki}_3`, which reaches ~:math:`5\times 10^{-6}`
absolute accuracy; the declared ``< 1e-5`` tolerance is now
~100× larger than the actual solver/reference error
(~:math:`10^{-7}`, same kernel on both sides).  See
:ref:`ki-table-construction` for the full postmortem.  The old
convergence-with-table-size regression
``tests/cp/test_verification.py::TestKi4Resolution`` was replaced by
``test_ki3_kernel_is_insensitive_to_n_ki_table``: ``n_ki_table``
is now a no-op and ``keff`` is bit-identical across
``{5000, 20000, 40000}``.

Additionally, **algebraic property tests** (``tests/cp/test_properties.py``)
are run for all three coordinate systems:

- **Row sums** :math:`= 1` (neutron conservation)
- **Reciprocity**: :math:`\Sigt{i} V_i P_{ij} = \Sigt{j} V_j P_{ji}`
- **Non-negativity**: :math:`P_{ij} \ge 0`
- **Homogeneous limit**: 1-region :math:`P = 1`


Extended Verification (CP-20260405-005)
---------------------------------------

31 additional tests in ``tests/cp/test_verification.py`` closing 9 QA gaps:

- **L0 P_inf comparison** (G-1): element-by-element solver vs derivation
  (tolerance :math:`< 10^{-10}`)
- **Multi-group properties** (G-6, W-2): row sums and reciprocity at 2G/4G
- **Upscatter** (G-2, W-1): eigenvalue with nonzero thermal-to-fast transfer
- **(n,2n)** (C-2, G-3, W-3): solver :math:`\keff` matches analytical with
  ``Sig2[0,0] = 0.01``
- **Optical limits** (G-4, G-7): thick (:math:`\tau \gg 1`,
  :math:`P_{ii} > 0.99`) and thin (:math:`\tau \ll 1`, high escape)
- **Convergence rate** (G-5): monotonic error decrease, dominance ratio
  estimation
- **8-region mesh** (G-8): mesh refinement convergence
- **GS inner iterations** (C-1, G-9): thermal > fast inner counts;
  no-self-scatter in 1; GS/Jacobi eigenvalue agreement
- **Ki4 table resolution** (W-6): diminishing returns from 5k to 40k points

Plus 36 diagnostic tests in ``tests/cp/test_diagnostics.py``.

::

   pytest tests/cp/ -v
