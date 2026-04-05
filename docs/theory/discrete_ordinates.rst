.. _theory-discrete-ordinates:

====================================
Discrete Ordinates Method (SN)
====================================

.. contents:: Contents
   :local:
   :depth: 3


Overview
========

The discrete ordinates (S\ :sub:`N`) method solves the **integro-differential
form** of the neutron transport equation by discretising the angular variable
:math:`\hat{\Omega}` into a finite set of directions (ordinates).  Unlike the
collision probability method (which works with the integral form and the scalar
flux), the S\ :sub:`N` method retains the **angular flux**
:math:`\psi(\mathbf{r}, \hat{\Omega}, E)` and resolves directional effects
such as streaming, anisotropic scattering, and angular current at interfaces.

The treatment follows [LewisMiller1984]_, with the angular discretisation
framework of [CaseZweifel1967]_.  Three coordinate systems are supported:

- **Cartesian** (slab / 2D) — the simplest case, no angular coupling
- **Spherical** (1D radial) — angular redistribution in :math:`\mu`
- **Cylindrical** (1D radial) — azimuthal redistribution per :math:`\mu`-level

The solver is implemented in :class:`SNSolver`, which satisfies the
:class:`~numerics.eigenvalue.EigenvalueSolver` protocol.  The convenience
function :func:`solve_sn` runs the full calculation and returns an
:class:`SNResult`.


Architecture: Base Geometry and Augmented Geometry
===================================================

The SN solver follows the same two-layer pattern as the CP solver:

1. **Base geometry** — :class:`~geometry.mesh.Mesh1D` or
   :class:`~geometry.mesh.Mesh2D` stores cell edges, material IDs, and
   coordinate system.

2. **Augmented geometry** — :class:`SNMesh` wraps the base mesh and an
   angular quadrature, precomputing the coordinate-specific streaming
   stencil:

   - **Cartesian**: ``streaming_x[n,i] = 2|mu_x|/dx[i]`` — the
     diamond-difference denominator terms (avoids per-cell division
     in the sweep hot loop).
   - **Spherical**: ``face_areas``, ``alpha_half`` — cell face areas
     :math:`4\pi r^2` and angular redistribution coefficients
     :math:`\alpha_{n+1/2} = \sum_{m=0}^{n} w_m \mu_m`.
   - **Cylindrical**: ``face_areas``, ``alpha_per_level`` — cell face
     areas :math:`2\pi r` and per-level azimuthal redistribution
     coefficients.

3. **Solver** — :func:`solve_sn` creates an ``SNMesh``, builds the
   ``SNSolver``, and runs power iteration.

.. code-block:: text

   Mesh1D / Mesh2D (base geometry)
       │
       ▼
   SNMesh (stencil + quadrature + α coefficients)
       │
       ▼
   solve_sn() → SNResult


The Neutron Transport Equation
===============================

Cartesian 1D
--------------

The steady-state transport equation in a 1D slab:

.. math::
   :label: sn-transport-1d

   \mu \frac{\partial \psi(x, \mu)}{\partial x}
   + \Sigt{} \, \psi(x, \mu)
   = \frac{Q}{W}

where :math:`\mu = \cos\theta` is the direction cosine, :math:`Q` is the
total isotropic source, and :math:`W = \sum_n w_n`.

Spherical 1D
--------------

In spherical coordinates, the transport equation acquires an **angular
redistribution term** that couples ordinates:

.. math::
   :label: sn-transport-spherical

   \mu \frac{\partial \psi}{\partial r}
   + \frac{1 - \mu^2}{r} \frac{\partial \psi}{\partial \mu}
   + \Sigt{} \psi = \frac{Q}{W}

The curvature term :math:`(1 - \mu^2)/r \cdot \partial\psi/\partial\mu`
arises because a neutron streaming radially at angle :math:`\mu` *rotates*
its direction as it moves to a different radius.  Discretising this term
requires diamond-difference in **both space and angle**.

Cylindrical 1D
----------------

For an infinitely long cylinder with azimuthal symmetry, the transport
equation in the radial variable :math:`r` is:

.. math::
   :label: sn-transport-cylindrical

   \frac{\eta}{r} \frac{\partial(r\psi)}{\partial r}
   - \frac{1}{r} \frac{\partial(\xi\psi)}{\partial\varphi}
   + \Sigt{} \psi = \frac{Q}{W}

where the direction cosines are:

- :math:`\eta = \sin\theta\cos\varphi` — radial projection (streaming)
- :math:`\xi = \sin\theta\sin\varphi` — azimuthal component
- :math:`\mu = \cos\theta` — axial component

The azimuthal redistribution :math:`-\partial(\xi\psi)/\partial\varphi`
couples ordinates on each :math:`\mu`-level.  Note the **minus sign** —
opposite to the spherical case — which has important numerical implications.


Multi-Group Extension
-----------------------

For :math:`G` energy groups, each transport equation becomes a coupled
system with scattering transfer :math:`\Sigs{g' \to g}` between groups:

.. math::
   :label: sn-multigroup

   \text{streaming} + \Sigt{g} \psi_g
   = \frac{1}{W} \sum_{g'} \Sigs{g' \to g} \phi_{g'}
   + \frac{\chi_g}{W k} \sum_{g'} \nSigf{g'} \phi_{g'}

where the streaming operator depends on the coordinate system.


Angular Quadratures
====================

ORPHEUS provides four angular quadrature types.

Gauss-Legendre (1D)
---------------------

For 1D slab geometry: :math:`N` points on :math:`\mu \in [-1, 1]`,
weights sum to 2.  Optimal for polynomial integrands (degree
:math:`2N-1` exact).  Implemented in :class:`GaussLegendre1D`.

Also used for spherical 1D, where the single direction cosine :math:`\mu`
suffices for the angular redistribution.

Lebedev (Sphere)
------------------

For 2D/3D Cartesian geometry: :math:`N` points on the unit sphere with
octahedral symmetry [Lebedev1999]_.  Weights sum to :math:`4\pi`.
Implemented in :class:`LebedevSphere`.

Level-Symmetric S\ :sub:`N`
------------------------------

For curvilinear geometries requiring :math:`\mu`-level structure.
Standard triangular quadrature with :math:`N/2` distinct :math:`\mu`
values per hemisphere.  Ordinates on each level are permutations of the
direction cosine set satisfying :math:`\eta^2 + \xi^2 + \mu^2 = 1`.

Weights sum to :math:`4\pi`.  Implemented in :class:`LevelSymmetricSN`.

.. warning::

   The current level-symmetric implementation groups :math:`\pm\mu_z`
   hemispheres on the same level, which is incompatible with the
   cylindrical sweep's per-level azimuthal redistribution.  Use
   :class:`ProductQuadrature` for cylindrical geometry until this is fixed.

Product Quadrature (GL × equispaced)
---------------------------------------

Tensor product of Gauss-Legendre in :math:`\mu = \cos\theta` (polar)
and equispaced points in :math:`\varphi` (azimuthal).  Each :math:`\mu`
level has the same number of azimuthal points, giving a clean level
structure ideal for the cylindrical sweep.

Weights sum to :math:`4\pi`.  Implemented in :class:`ProductQuadrature`.

Comparison
-----------

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 15 20

   * - Quadrature
     - Geometry
     - :math:`\sum w`
     - Level structure
     - Best for
   * - Gauss-Legendre
     - Slab, Sphere
     - 2
     - No
     - 1D problems
   * - Lebedev
     - Cartesian 2D
     - :math:`4\pi`
     - No
     - 2D/3D Cartesian
   * - Level-Symmetric
     - Sphere, Cylinder*
     - :math:`4\pi`
     - Yes
     - Curvilinear
   * - Product
     - Cylinder
     - :math:`4\pi`
     - Yes
     - Cylindrical 1D


Spatial Discretisation: Diamond Difference
===========================================

Cartesian Diamond-Difference
------------------------------

For ordinate :math:`n` with :math:`\mu_n > 0` traversing cell :math:`i`
of width :math:`\Delta x_i`, the diamond-difference approximation gives:

.. math::
   :label: diamond-difference

   \psi_{\text{avg},i}
   = \frac{Q_i / W + \frac{2|\mu_n|}{\Delta x_i} \psi_{\text{in}}}
          {\Sigt{i} + \frac{2|\mu_n|}{\Delta x_i}}

The outgoing face flux:

.. math::
   :label: dd-outgoing

   \psi_{\text{out}} = 2 \psi_{\text{avg}} - \psi_{\text{in}}


Spherical Diamond-Difference
------------------------------

For spherical 1D, the balance equation at cell :math:`i`, ordinate
:math:`n` is:

.. math::
   :label: dd-spherical

   \frac{\mu_n}{V_i}
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \frac{1}{V_i}
   \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
   + \Sigt{} \psi_{n,i} = \frac{Q}{W}

where :math:`A = 4\pi r^2`, :math:`V = \frac{4}{3}\pi(r_{\rm out}^3 - r_{\rm in}^3)`,
and the angular redistribution coefficients are:

.. math::
   :label: alpha-spherical

   \alpha_{n+\frac12} = \sum_{m=0}^{n} w_m \mu_m

with :math:`\alpha_{\frac12} = \alpha_{N+\frac12} = 0` (by Gauss-Legendre
antisymmetry).

Diamond-difference closures are applied in **both** space and angle:

- Spatial: :math:`\psi_{n,i} = \frac{1}{2}(\psi_{\rm in} + \psi_{\rm out})`
- Angular: :math:`\psi_{n,i} = \frac{1}{2}(\psi_{n-\frac12} + \psi_{n+\frac12})`

Ordinates are processed from most negative :math:`\mu` to most positive.
Negative-:math:`\mu` ordinates sweep inward (outer → centre); positive
sweep outward (centre → outer).  At :math:`r = 0`, the face area
:math:`A = 0` so spatial streaming vanishes — the angular redistribution
provides the inward-to-outward transition.

Implemented in :func:`_sweep_1d_spherical`.


Cylindrical Diamond-Difference
---------------------------------

For cylindrical 1D, the balance equation at cell :math:`i`,
:math:`\mu`-level :math:`p`, azimuthal ordinate :math:`m` is:

.. math::
   :label: dd-cylindrical

   \frac{\eta_{p,m}}{V_i}
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   - \frac{1}{V_i}
   \bigl[\alpha_{p,m+\frac12}\psi_{m+\frac12} - \alpha_{p,m-\frac12}\psi_{m-\frac12}\bigr]
   + \Sigt{} \psi = \frac{Q}{W}

where :math:`A = 2\pi r`, :math:`V = \pi(r_{\rm out}^2 - r_{\rm in}^2)`,
and the azimuthal redistribution coefficients per level are:

.. math::
   :label: alpha-cylindrical

   \alpha_{p,m+\frac12} = \alpha_{p,m-\frac12} + w_{p,m} \xi_{p,m}

Note the **minus sign** before the :math:`\alpha` terms — opposite to
spherical.  This has important numerical implications: the standard
diamond-difference with absolute values :math:`|\alpha|` (which works for
spherical) gives the wrong sign convention for cylindrical, causing
divergence on heterogeneous problems.

.. warning::

   The cylindrical sweep currently produces exact eigenvalues for
   **homogeneous** problems but incorrect, divergent eigenvalues for
   **heterogeneous** problems.  The fix requires the Lewis & Miller
   starting-direction treatment (§4.5.4), which tracks the product
   :math:`\alpha\psi` instead of :math:`\psi` alone.  See
   ``02.Discrete.Ordinates/TODO_cylindrical_dd.md`` for full analysis.

Implemented in :func:`_sweep_1d_cylindrical`.


Geometry Comparison
--------------------

.. list-table::
   :header-rows: 1
   :widths: 18 27 27 28

   * - Aspect
     - Cartesian
     - Spherical
     - Cylindrical
   * - Streaming
     - :math:`\mu \partial\psi/\partial x`
     - :math:`\mu A \partial\psi/\partial r / V`
     - :math:`\eta A \partial\psi/\partial r / V`
   * - Redistribution
     - None
     - :math:`+\alpha \partial\psi/\partial\mu / V`
     - :math:`-\alpha \partial\psi/\partial\varphi / V`
   * - Face area
     - 1 (per unit area)
     - :math:`4\pi r^2`
     - :math:`2\pi r`
   * - Volume
     - :math:`\Delta x`
     - :math:`\frac{4}{3}\pi \Delta(r^3)`
     - :math:`\pi \Delta(r^2)`
   * - :math:`\alpha` scope
     - N/A
     - Global (all ordinates)
     - Per :math:`\mu`-level
   * - :math:`\alpha` sign
     - N/A
     - Positive
     - Negative
   * - Quadrature
     - GL or Lebedev
     - GL
     - Product or Level-Sym


Sweep Algorithm
================

Because each cell's outgoing flux becomes the next cell's incoming flux,
the equations are solved in the direction of neutron travel.  Three
sweep implementations exist:

1. **Cartesian 1D cumprod**: Exploits GL symmetry to solve the DD
   recurrence via cumulative products — O(N) vectorised numpy ops.
   Very fast (~ms).

2. **Cartesian 2D wavefront**: Sweeps cells along anti-diagonals
   :math:`i + j = k`, vectorised within each diagonal.  Four sweep
   directions for the four octants.

3. **Curvilinear 1D**: Cell-by-cell, ordinate-by-ordinate (angular
   coupling prevents vectorisation across ordinates).  For spherical:
   global ordinate loop.  For cylindrical: per-level azimuthal loop.

Dispatch is automatic based on ``SNMesh.curvature``:

.. code-block:: python

   if sn_mesh.curvature == "spherical":
       return _sweep_1d_spherical(...)
   elif sn_mesh.curvature == "cylindrical":
       return _sweep_1d_cylindrical(...)
   elif is_gl_1d:
       return _sweep_1d_cumprod(...)
   else:
       return _sweep_2d_wavefront(...)


Boundary Conditions
====================

ORPHEUS uses **reflective boundary conditions** on all faces, representing
an infinite lattice.  At a reflective boundary:

.. math::
   :label: reflective-bc

   \psi_n^{\text{in}} = \psi_{n'}^{\text{out}}

where :math:`n'` is the reflected partner ordinate.

For curvilinear geometries:

- **Outer** (:math:`r = R`): reflective, as in Cartesian
- **Inner** (:math:`r = 0`): face area :math:`A = 0`, so no spatial flux
  crosses the origin.  The angular redistribution handles the
  inward-to-outward transition.


Scattering Anisotropy: P\ :sub:`N` Expansion
===============================================

P\ :sub:`0` isotropic scattering adds a direction-independent source.
P\ :sub:`N` (:math:`N \geq 1`) adds per-ordinate sources via Legendre
moments of the angular flux:

.. math::
   :label: pn-scatter-source

   Q_{\text{scatter},g}(\hat{\Omega}_n)
   = \sum_{\ell=0}^{L} (2\ell+1)
     \sum_{g'} \Sigs{g' \to g}^{(\ell)}
     \left[ \sum_{m=-\ell}^{\ell}
       f_{\ell,g'}^m \, Y_\ell^m(\hat{\Omega}_n) \right] / W

See :meth:`SNSolver._build_aniso_scattering` for the implementation.


The Eigenvalue Problem
=======================

The eigenvalue :math:`\keff` is determined by **power iteration**: an outer
loop updates :math:`k` from the production/absorption ratio, with an inner
loop that solves the within-group scattering problem.

Two inner solvers are available:

1. **Source iteration** (sweep-based): fixed-point iteration inverting
   :math:`T` via transport sweeps.  Convergence rate governed by
   spectral radius of :math:`T^{-1}S`.  Works for all geometries.

2. **BiCGSTAB** (direct operator): forms :math:`T` explicitly via
   finite differences and solves :math:`T\psi = b` with a Krylov solver.
   Converges in ~100 iterations regardless of scattering ratio.
   Available for Cartesian and spherical 1D (1G only for spherical).

.. warning::

   The two inner solvers use different spatial discretisations (DD sweep
   vs FD gradient) that produce slightly different :math:`\keff` on
   coarse meshes.  They converge to the same answer as :math:`h \to 0`.


Verification
=============

Homogeneous Infinite Medium
----------------------------

For homogeneous geometry with reflective BCs, the flux is spatially flat
and :math:`\keff = \lambda_{\max}(\mathbf{A}^{-1}\mathbf{F})`.  This is
geometry-independent — Cartesian, spherical, and cylindrical all give
the same :math:`\keff`.

.. list-table::
   :header-rows: 1
   :widths: 12 15 20 20 20

   * - Groups
     - :math:`\kinf`
     - Cartesian (GL S8)
     - Spherical (GL S8)
     - Cylindrical (Prod 4×8)
   * - 1
     - 1.5000
     - exact
     - exact
     - exact
   * - 2
     - 1.8750
     - exact
     - < 2%
     - exact
   * - 4
     - 1.4878
     - exact
     - < 2%
     - exact

Cartesian and cylindrical homogeneous are exact to machine precision.
Spherical has ~1% discretization error from the angular DD on a 20-cell
S8 mesh (1G is exact because :math:`k` is flux-shape-independent).

Spatial and Angular Convergence
--------------------------------

The diamond-difference scheme converges at :math:`O(h^2)` with mesh
refinement.  Gauss-Legendre quadrature shows spectral convergence in
angle.  Both are verified in ``test_sn_1d.py``.

Property Tests
---------------

For all geometries:

- **Particle balance**: :math:`\text{production}/\text{absorption} = \keff`
- **Flux non-negativity**: :math:`\phi \geq 0` everywhere
- **Angular flux at** :math:`r = 0` **all positive** (curvilinear only)
- **Multi-group eigenvector not flat**: flux spectrum differs between
  fuel and moderator (catches 1G-degenerate bugs)

Run the full suite::

   pytest tests/test_sn_1d.py tests/test_sn_properties.py \
          tests/test_sn_solver_components.py tests/test_sn_spherical.py \
          tests/test_sn_cylindrical.py tests/test_sn_quadrature.py \
          tests/test_sn_sweep_regression.py -v -m "not slow"


References
==========

.. [LewisMiller1984] E.E. Lewis and W.F. Miller, Jr.,
   *Computational Methods of Neutron Transport*,
   John Wiley & Sons, 1984.

.. [CaseZweifel1967] K.M. Case and P.F. Zweifel,
   *Linear Transport Theory*,
   Addison-Wesley, 1967.

.. [Lebedev1999] V.I. Lebedev and D.N. Laikov,
   "A quadrature formula for the sphere of the 131st algebraic order
   of accuracy," *Doklady Mathematics*, 59(3):477–481, 1999.
