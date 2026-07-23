.. _theory-verification-summary:

================
Evidence Summary
================

This page is the **results compilation** of the verification part —
the cross-suite view of what every solver family is verified against
and where each body of evidence lives. It is the middle layer of the
part's three-layer architecture: :doc:`principles` owns the
*doctrine* (what counts as evidence, the L0–L3 ladder, the pillars),
the per-method chapters own the *evidence narratives* (case by case,
with derivations and measured numbers), and this page compiles the
cross-method view — plus the run-book for executing the suite.

Two conventions keep this page honest as the suite grows:

- **Structure here, totals there.** Hand-written tables on this page
  carry only *structural* facts — which reference classes exist,
  which case grids are designed, which properties are pinned. Test
  *counts* are never hand-written: the auto-generated :doc:`matrix`
  reports the live totals on every build, and the reference-case
  table below is emitted from the case registry itself.
- **The authorities are upstream.** The case rows come from
  ``orpheus.derivations.reference_values``; the level tallies from
  ``tests._harness.registry.TEST_REGISTRY``; the per-method accounts
  from the chapters. This page points; it does not duplicate.


The evidence map — what each method is verified against
=======================================================

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Solver family
     - Verified against
   * - :doc:`homogeneous`
     - The closed-form infinite-medium matrix eigenvalue (exact to
       :math:`10^{-12}`) — the shared analytical anchor every
       family's homogeneous limit consumes — plus operator-object
       gates that pin the assembled :math:`K = A^{-1}F` resolvent as
       a *matrix*, not merely its spectrum.
   * - :doc:`sn`
     - The homogeneous exact limit; the MMS ladder (1D slab,
       heterogeneous 2-group, 2D Cartesian, P1 anisotropic,
       curvilinear isotropic and anisotropic, linear-discontinuous,
       prescribed-inflow) for convergence-order and flux-shape
       claims; the Case singular-eigenfunction heterogeneous
       eigenvalue; property gates and the numerical-sensitivity
       record.
   * - :doc:`collision_probability`
     - Semi-analytical :math:`E_3`/:math:`\mathrm{Ki}_4` eigenvalues
       over the designed 27-case grid ({1, 2, 4} groups × {1, 2, 4}
       regions × {slab, cylinder, sphere}); the Peierls–Nyström
       heterogeneous reference (30+ digit Nyström collocation);
       CP-matrix property gates.
   * - :doc:`method_of_characteristics`
     - The flat-source analytical eigenvalue; the manufactured
       continuous solution for the pin-cell spatial operator
       (:ref:`moc-mms-verification`); ray-spacing, azimuthal, and
       polar convergence gates; CP cross-verification.
   * - :doc:`monte_carlo`
     - Statistical z-score gates (:math:`z < 5\sigma`) against
       analytical and CP-referenced eigenvalues; 1-group
       zero-variance determinism; :math:`\sigma \propto 1/\sqrt{N}`
       convergence.
   * - :doc:`diffusion`
     - The bare-slab buckling eigenvalue (analytical); 2-region
       transcendental transfer-matrix references with measured
       :math:`O(h^2)` orders (:ref:`diffusion-2rg-verification`).
       The Richardson heterogeneous reference is retired
       (:ref:`richardson-extrapolation`).

The claim-layer, ladder, and tier vocabulary used throughout is
defined at :ref:`verification-evidence-classes`; the cross sections
feeding these cases come from the real nuclear-data pipeline
(:ref:`theory-cross-section-data`).


Reference cases — the registry view
===================================

One row per case in the continuous-reference registry
(``orpheus.derivations.reference_values``, which auto-discovers every
``continuous_cases()`` entry point under ``orpheus.derivations.*``).
The **Tolerance** column is the declared case tolerance — a fixed
floor for the deterministic solvers, a statistical band
(:math:`z < 5\sigma`) for Monte Carlo, and a convergence-order
declaration (:math:`O(h^2)`) for the diffusion continuous-reference
case; the tolerance ladder's rationale lives at :doc:`principles`.
The table is regenerated at every documentation build (and manually
via ``python -m orpheus.derivations.generate_rst``), so the rows
always match the registry.

This is the *registry* view only. The MMS ladders, property gates,
cross-method agreement gates, and sensitivity studies do not appear
here — they live in the per-method chapters mapped above.

.. include:: ../../_generated/verification_table.rst


.. _verification-cross-method-coverage:

Cross-method coverage
=====================

The cross-method regression net pairs the continuous-reference
families that solve the same benchmark problems through
structurally-independent mathematics — ``fn_method`` (Case
singular-eigenfunction collocation) and ``trajectory_resolvent``
(bouncing-trajectory Green's function) — on shared registry cases
with per-solver tolerances. :doc:`cross_method` owns the protocol:
the adapter contract, the agreement-tolerance discipline (never
tighter than the looser truth tolerance), the L4 ruling, and the
recipes for adding solvers and cases.

The as-designed coverage — the rows mirror
``tests/cross_method/cases.py``, and the live case × adapter matrix
is printed by
:func:`~tests.cross_method.test_eigenvalue.test_coverage_matrix_diagnostic`
on every run:

.. list-table::
   :header-rows: 1
   :widths: 28 8 22 24 18

   * - Case set
     - Cases
     - Families
     - L1 truth backing
     - Agreement gate
   * - Bare-critical slab
     - 4
     - F\ :sub:`N` + trajectory
     - Sood / KLL
     - pairwise (slow)
   * - Bare-critical sphere
     - 3
     - F\ :sub:`N` + trajectory
     - Sood / KLL
     - pairwise
   * - Reflected slab
     - 4
     - F\ :sub:`N`
     - Sood + NM 1980
     - one-sided
   * - Closed sphere :math:`k_\infty`
     - 1
     - trajectory
     - :math:`V_{\alpha 1}` identity
     - one-sided
   * - GS Table XI parametric
     - 5
     - F\ :sub:`N`
     - Grandjean–Siewert 1979
     - one-sided

The two one-sided asymmetries are deliberate, not gaps in waiting:
the closed-sphere :math:`k_\infty` identity is covered on the
F\ :sub:`N` side separately (its ``compute_kinf_*`` closed forms are
gated against ``kinf_homogeneous`` in the per-method file), and the
GS parametric set has no ``trajectory_resolvent`` registry case. The
known multi-group coverage gap is acknowledged and scoped at
:doc:`cross_method`.

Beyond the eigenvalue gates, a foundation-level **polymorphism net**
(``tests/cross_method/test_polymorphism.py``) pins direct
construction of each family's math-heart class (``MomentSpace``,
``Billiard``) against its adapter dispatch — the agreement contract
that survived the Phase-D retirement of the ``TransportSolver``
Protocol.


Structural property tests
=========================

Every math-bearing solver object ships tests of its *defining* laws —
conservation, reciprocity, positivity, symmetry — independent of any
eigenvalue claim. The per-method property files are the source of
truth; the inventory below compiles the pinned classes.

CP Matrix Properties
--------------------

For every CP case (``tests/cp/test_properties.py``), the collision
probability matrix :math:`P_{\infty}` must satisfy:

- **Row sums = 1** (neutron conservation): every neutron born in
  region *i* must have its first collision somewhere.
  :math:`\sum_j P_{\infty}(i,j,g) = 1 \; \forall \, i, g`.
  Tolerance: < 10⁻¹⁰.
- **Reciprocity**: :math:`\Sigma_{t,i} V_i P(i,j) = \Sigma_{t,j} V_j
  P(j,i)`. This is detailed balance — a consequence of time-reversal
  symmetry of the transport equation.  Tolerance: < 10⁻¹⁰.
- **Non-negativity**: :math:`P(i,j,g) \geq 0 \; \forall \, i, j, g`.
- **Homogeneous limit**: a 1-region cell must give
  :math:`P(0,0) = 1`.

SN Properties
-------------

Quadrature and 1-D solution properties
(``tests/sn/primitives/test_properties.py``):

- **GL quadrature weights**: must sum to 2 (measure of [-1, 1]).
- **GL symmetry**: :math:`\mu_i = -\mu_{N-1-i}`.
- **Flux flatness**: a homogeneous slab must have exactly flat flux.
- **Particle balance**: with reflective BCs (no leakage),
  production / absorption = :math:`k_{\text{eff}}`.

The geometry-wide property gates — particle balance across all
geometries, flux non-negativity, curvilinear :math:`r = 0`
positivity, non-flat multi-group eigenvector — are catalogued in the
SN chapter's Property Tests section (:doc:`sn`).

Diffusion Properties
--------------------

(``tests/diffusion/test_properties.py``)

- **Vacuum is the Marshak law, not zero flux**: ``BC("vacuum")``
  means zero *incoming current* — :math:`J^- = 0` is asserted at
  10⁻¹² of the boundary-trace scale, while the boundary-cell flux
  must stay strictly positive (vacuum is not the zero-flux Dirichlet
  idealization).
- **Flux positivity**: all flux values positive in the fundamental
  mode.
- **Flux symmetry**: bare-slab flux is mirror-symmetric about the
  center, with matching outward net currents on the two faces.

MOC Properties
--------------

(``tests/moc/test_properties.py``)

- **Particle balance**: production / absorption =
  :math:`k_{\text{eff}}`.
- **Flux positivity**: :term:`scalar flux` > 0 everywhere.
- **Per-material flux**: the MOC per-material scalar flux matches the
  volume-averaged scalar (``test_flux_per_material_matches_scalar``).
- **Heterogeneous flux depression**: thermal flux is higher in the
  moderator than the fuel (``test_heterogeneous_flux_depression``).

MC Properties
-------------

(``tests/mc/test_properties.py``)

- **Geometry protocol**: ``ConcentricPinCell`` and ``SlabPinCell``
  return correct material IDs at known positions.
- **1G deterministic**: homogeneous 1-group MC has :math:`\sigma = 0`
  (all neutrons see identical cross sections).
- **Determinism**: same seed → identical results; different seeds →
  different histories (``tests/mc/test_gaps.py``).


Convergence studies
===================

Beyond point-value verification, the suite pins discretisation-error
*rates* — with the standing L2 caveat of :ref:`vv-level-ladder`: a
correct rate never certifies the converged-to value.

- **SN spatial (diamond difference)** — the observed order
  approaches 2.0 under mesh refinement, confirming the
  :math:`O(h^2)` truncation error
  (``tests/sn/eigenvalue/test_keff_slab.py::test_spatial_convergence``).
  The MMS ladder measures orders per scheme and geometry
  (:doc:`sn`).
- **SN angular (Gauss–Legendre)** — the eigenvalue error decreases
  faster than any polynomial in :math:`1/N`, confirming spectral
  convergence of the GL quadrature
  (``tests/sn/eigenvalue/test_keff_slab.py::test_angular_convergence``).
- **MOC** — ray-spacing, azimuthal, and polar refinement gates
  (:doc:`method_of_characteristics`).
- **Diffusion** — measured :math:`O(h^2)` orders against the
  2-region transcendental references
  (``tests/diffusion/test_continuous_reference.py``,
  :ref:`diffusion-2rg-verification`).
- **MC statistics** — :math:`\sigma \propto 1/\sqrt{N}`, estimator
  bias, and inactive-cycle gates (``tests/mc/test_convergence.py``).
- **Cross-solver limit** — fine-mesh SN approaches the CP reference
  on a shared slab case (``tests/test_convergence.py``).


The auto-generated verification matrix
======================================

:doc:`matrix` is regenerated on every Sphinx build from the pytest
test registry (``tests._harness.registry.TEST_REGISTRY``) — it is the
*live* census this page deliberately never duplicates: the V&V level
distribution and per-module level grid, the tagging-source
breakdown, per-equation coverage (which tests pin which equation
labels), orphan and documented-only equations, scan-exempt files,
phantom ``verifies`` targets, L0 error-catalog coverage, and the
unmarked-test list. How tests *declare* levels and equation links —
the marker and sentinel contract the matrix renders — is
:doc:`harness`.


Running the Tests
=================

The canonical invocation is ``python -O -m pytest`` — the production
path strips ``assert`` statements, so gate runs match what
production executes. The one deliberate exception: bare-``assert``
canary gates must run *without* ``-O`` (failure mode 8 at
:ref:`verification-failure-modes`), which is why the sentinel net
below drops the flag. Counts are deliberately not quoted here — the
suite grows with every merge; the auto-generated :doc:`matrix`
reports the current totals.

.. code-block:: bash

   # Install test dependencies
   pip install -e ".[test]"

   # The pre-merge gate: the full tree, non-slow, single-process
   python -O -m pytest -m "not slow"

   # The slow tiers: Peierls–Nyström references, SN MMS/curvilinear
   # L1 ladders, MC high-statistics, cross-method agreement gates
   python -O -m pytest -m slow

   # The seconds-fast SN canary net (bare asserts — no -O here)
   python -m pytest -m sentinel

   # A specific solver's tests
   python -O -m pytest tests/homogeneous -v
   python -O -m pytest tests/cp/test_properties.py -v

Two invocation regimes, deliberately distinct: the **pre-merge
gate** is the patient full-tree *serial* run above; for **inner-loop
iteration** on the SN tree, prefer per-capability-tier invocations —
the whole SN suite in one process is memory-heavy. ``pytest-xdist``
is version-pinned and reserved for within-tier parallelism; the
``[test]`` extra's inline notes in ``pyproject.toml`` carry the
operational detail.
