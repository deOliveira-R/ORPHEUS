.. _theory-verification-summary:

================
Evidence Summary
================

.. This page is the designated results-compilation home — the
   cross-suite evidence tables, the structural unit-test inventory,
   the convergence studies, and the run-book. It is authored properly
   at task #10 stage V6; the sections below were moved verbatim at
   stage V3 from the dissolved ``docs/theory/verification.rst`` so
   the part restructure loses nothing.


Reference Cases
===============

.. include:: ../../_generated/verification_table.rst


Unit Tests
==========

Beyond eigenvalue verification, the test suite verifies structural properties
of each solver's intermediate quantities.

CP Matrix Properties
--------------------

For every CP case (slab and cylinder), the collision probability matrix
:math:`P_{\infty}` must satisfy:

- **Row sums = 1** (neutron conservation): every neutron born in region *i*
  must have its first collision somewhere.
  :math:`\sum_j P_{\infty}(i,j,g) = 1 \; \forall \, i, g`.
  Tolerance: < 10⁻¹⁰.

- **Reciprocity**: :math:`\Sigma_{t,i} V_i P(i,j) = \Sigma_{t,j} V_j P(j,i)`.
  This is detailed balance — a consequence of time-reversal symmetry of the
  transport equation.  Tolerance: < 10⁻¹⁰.

- **Non-negativity**: :math:`P(i,j,g) \geq 0 \; \forall \, i, j, g`.

- **Homogeneous limit**: a 1-region cell must give :math:`P(0,0) = 1`.

SN Properties
-------------

- **GL quadrature weights**: must sum to 2 (measure of [-1,1]).
- **GL symmetry**: :math:`\mu_i = -\mu_{N-1-i}`.
- **Flux symmetry**: homogeneous slab must have spatially flat flux.
- **Particle balance**: with reflective BCs (no leakage),
  production / absorption = keff.

Diffusion Properties
--------------------

- **Vacuum BC**: flux at boundary cells is small compared to peak.
- **Flux positivity**: all flux values positive in fundamental mode.
- **Flux symmetry**: bare slab flux is symmetric about the center.

MOC Properties
--------------

- **Particle balance**: production / absorption = keff.
- **Flux positivity**: :term:`scalar flux` > 0 everywhere.
- **Per-material flux**: the MOC per-material scalar flux matches the
  volume-averaged scalar (``test_flux_per_material_matches_scalar``).
- **Heterogeneous flux depression**: thermal flux is higher in the
  moderator than the fuel (``test_heterogeneous_flux_depression``).

MC Properties
-------------

- **Geometry protocol**: ``ConcentricPinCell`` and ``SlabPinCell`` return
  correct material IDs at known positions.
- **1G deterministic**: homogeneous 1-group MC has σ = 0 (all neutrons
  see identical cross sections).


Convergence Studies
===================

Beyond point-value verification, the test suite checks that discretisation
errors decrease at the expected rate:

**SN 1D spatial convergence** (diamond-difference):
  The observed convergence order should approach 2.0 as the mesh is refined,
  confirming the O(h²) truncation error of the diamond-difference scheme.

**SN 1D angular convergence** (Gauss-Legendre):
  The eigenvalue error should decrease faster than any polynomial in
  1/N, confirming spectral convergence of the GL quadrature.

**Diffusion spatial convergence**:
  The three-point finite-difference stencil gives O(h²) convergence,
  verified against the analytical buckling eigenvalue.


.. seealso::

   Verified solvers and their theory pages:

   - :ref:`theory-homogeneous` — matrix eigenvalue (analytical, exact)
   - :ref:`theory-collision-probability` — E₃/Ki₄ semi-analytical eigenvalue
   - :ref:`theory-discrete-ordinates` — homogeneous exact + :ref:`Richardson <richardson-extrapolation>` heterogeneous
   - :ref:`theory-method-of-characteristics` — homogeneous exact + Richardson heterogeneous
   - :ref:`theory-monte-carlo` — z-score against analytical + CP reference

   Cross sections: :ref:`theory-cross-section-data` (real nuclear data pipeline).


Running the Tests
=================

The canonical invocation is ``python -O -m pytest`` (the production
path strips ``assert``; see the harness architecture page for the
``-O`` policy and the marker taxonomy).  Counts are deliberately not
quoted here — the suite grows with every merge; the auto-generated
verification matrix reports the current totals.

.. code-block:: bash

   # Install test dependencies
   pip install -e ".[test]"

   # The pre-merge gate: the full tree, non-slow
   python -O -m pytest -m "not slow"

   # The slow tiers (Richardson, MC high-stats, curvilinear L1)
   python -O -m pytest -m slow

   # A specific solver's tests
   python -O -m pytest tests/homogeneous/ -v
   python -O -m pytest tests/cp/test_properties.py -v
