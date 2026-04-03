Verification Suite
==================

ORPHEUS verifies every transport solver against analytical reference
solutions derived from each method's own mathematical equations using
SymPy. Each derivation is **self-contained**: it starts from the
solver's formulation and derives the expected eigenvalue independently.
The same code that produces the LaTeX equations in this chapter also
produces the reference values consumed by the pytest test suite.

Methodology
-----------

Each verification case defines:

1. **Cross sections** — synthetic macroscopic data for abstract regions
   (not specific materials). This isolates the numerical method from
   the cross-section processing pipeline.

2. **Analytical eigenvalue** — derived from the solver's own equations
   using SymPy (symbolic where possible, numerical for special functions).

3. **Solver tolerance** — method-specific:

   - Infinite-medium solvers: machine precision (< 10\ :sup:`-10`)
   - CP solvers vs analytical CP eigenvalue: < 10\ :sup:`-5`
   - SN (homogeneous): < 10\ :sup:`-4`
   - MOC (homogeneous): < 10\ :sup:`-2`
   - Diffusion: spatial O(h²) convergence
   - Monte Carlo: statistical (z-score < 5σ)

Reference Cases
---------------

.. include:: ../_generated/verification_table.rst

Homogeneous Infinite Medium
---------------------------

For an infinite homogeneous medium, scattering does not change the
eigenvalue because leakage is zero. The eigenvalue depends only on
absorption, fission, and the fission spectrum.

.. include:: ../_generated/homogeneous_derivation.rst

Discrete Ordinates (SN)
-----------------------

The SN method discretises the angular variable into a finite set of
directions. For a homogeneous medium with reflective boundary
conditions, the derivation starts from the 1D SN transport equation
and shows that the spatially-flat, isotropic solution is exact.

.. include:: ../_generated/sn_derivation.rst

Slab Collision Probability
--------------------------

The slab CP method uses the E₃ exponential integral kernel to compute
first-collision probabilities in a 1D half-cell with reflective centre
and white boundary condition at the cell edge.

.. include:: ../_generated/cp_slab_derivation.rst

Cylindrical Collision Probability
---------------------------------

The cylindrical CP method uses the Ki₃/Ki₄ Bickley–Naylor kernel for
a Wigner–Seitz cell with annular regions and white boundary condition.

.. include:: ../_generated/cp_cylinder_derivation.rst

Method of Characteristics (MOC)
-------------------------------

The MOC method solves the transport equation along characteristic rays.
For a homogeneous medium, the derivation starts from the characteristic
ODE and shows that the flat-source solution is exact along every ray.

.. include:: ../_generated/moc_derivation.rst

Monte Carlo
-----------

The Monte Carlo method simulates neutron random walks. For a
homogeneous medium, the derivation uses collision probability theory
to compute the expected multiplication factor from first principles.

.. include:: ../_generated/mc_derivation.rst

Diffusion (Buckling Eigenvalue)
-------------------------------

The 1D finite-difference diffusion solver is verified against the
analytical buckling eigenvalue for a bare homogeneous slab.

.. include:: ../_generated/diffusion_derivation.rst

Convergence Studies
-------------------

Beyond point-value verification, the test suite checks that
discretisation errors decrease at the expected rate:

**SN 1D spatial convergence** (diamond-difference):
  The observed order should approach 2.0 as the mesh is refined,
  confirming the O(h²) truncation error of the diamond-difference scheme.

**SN 1D angular convergence** (Gauss-Legendre):
  The eigenvalue error should decrease faster than any polynomial in
  1/N, confirming spectral convergence of the GL quadrature.

**Diffusion spatial convergence**:
  The three-point finite-difference stencil gives O(h²) convergence,
  verified against the analytical buckling eigenvalue.

Running the Tests
-----------------

.. code-block:: bash

   # Install test dependencies
   pip install -e ".[test]"

   # Run all non-slow tests (~80s)
   pytest

   # Run including slow tests (MC high-stats, 2D SN convergence, ~7min)
   pytest -m "not slow"
   pytest -m "slow"

   # Run a specific solver's tests
   pytest tests/test_homogeneous.py -v
