.. _theory-discrete-ordinates:

==========================================
Discrete Ordinates Method (S\ :sub:`N`)
==========================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before modifying the SN solver.**

- Transport equation: :math:`\mu_m \frac{\partial\psi_m}{\partial x} + \Sigma_t \psi_m = Q/2` (1D slab)
- Curvilinear adds angular redistribution: :math:`\alpha` coefficients couple ordinates
- Diamond difference: :math:`\psi^a = (1+\beta)\psi_{\text{out}} - \beta\psi_{\text{in}}`, Morel-Montry sets :math:`\beta = 0`
- Scattering convention: ``SigS[l][g_from, g_to]`` — source uses **transpose**: ``Q = SigS^T @ phi``
- GL weights sum to 2; Lebedev/LS/Product sum to :math:`4\pi`
- **Gotcha**: 1-group tests are degenerate (k = νΣ_f/Σ_a regardless of flux shape)
- **Gotcha**: homogeneous tests are blind to curvilinear redistribution bugs (flat flux → α terms vanish)
- **Gotcha**: conservation holds even with wrong per-ordinate balance (telescoping sum identity)
- The :math:`\alpha` dome must be non-negative; negative → NaN/overflow
- Fixed-source flat-flux diagnostic (Q/Σ_t) is the most powerful curvilinear bug detector
- Key reference: Bailey, Morel & Chang (2009) — Eq. 50 (α recursion), Eq. 74 (M-M weights)
- Verification uses :ref:`synthetic cross sections <synthetic-xs-library>`, not real nuclear data

.. admonition:: Conventions

   - Scattering matrix: :ref:`scattering-matrix-convention` — ``SigS[g_from, g_to]``, source uses transpose
   - Multi-group balance: :eq:`mg-balance` in :ref:`theory-homogeneous`
   - Cross sections: :ref:`theory-cross-section-data`
   - Verification: :ref:`synthetic-xs-library` — regions A/B/C/D
   - Eigenvalue: :ref:`power-iteration-algorithm` shared with all deterministic solvers


Overview
========

The discrete ordinates (S\ :sub:`N`) method solves the
:ref:`multi-group eigenvalue problem <mg-eigenvalue-problem>` in
integro-differential form by discretising the angular variable
:math:`\hat{\Omega}` into a finite set of directions (ordinates).  Unlike the
collision probability method (which works with the integral form and the scalar
flux), S\ :sub:`N` retains the **angular flux**
:math:`\psi(\mathbf{r}, \hat{\Omega}, E)` and resolves directional effects
such as streaming, anisotropic scattering, and angular current at interfaces.

Three coordinate systems are supported:

- **Cartesian** (slab / 2D) --- the simplest case; no angular coupling
  between ordinates.
- **Spherical** (1D radial) --- angular redistribution in :math:`\mu`;
  a single dome of :math:`\alpha` coefficients couples all ordinates.
- **Cylindrical** (1D radial) --- azimuthal redistribution per
  :math:`\mu`-level; independent :math:`\alpha` domes on each level.

All three share a single balance-equation framework with a geometry
factor :math:`\Delta A / w` that guarantees per-ordinate flat-flux
consistency.  The treatment follows [Bailey2009]_ for the curvilinear
formulation, [LewisMiller1984]_ for the general framework, and
[CaseZweifel1967]_ for the angular discretisation.

The solver is implemented in :class:`SNSolver`, which satisfies the
:class:`~numerics.eigenvalue.EigenvalueSolver` protocol.  The convenience
function :func:`solve_sn` runs the full calculation and returns an
:class:`SNResult`.

Because the protocol puts the scattering source *inside*
``solve_fixed_source``, the inner source iteration (which converges the
in-scatter and anisotropic source) stays encapsulated in the SN-specific
sweep — the outer :func:`~numerics.eigenvalue.power_iteration` loop is
identical to the one used by CP, MOC, diffusion, and the homogeneous
solver.  See :doc:`../api/numerics` for the protocol contract.


Architecture
============

Two-Layer Mesh Pattern
----------------------

The S\ :sub:`N` solver follows the same two-layer pattern as the CP
solver.  This pattern (base :class:`~geometry.mesh.Mesh1D` + augmented
mesh) is shared with :ref:`theory-collision-probability` and
:ref:`theory-method-of-characteristics`.

1. **Base geometry** --- :class:`~geometry.mesh.Mesh1D` or
   :class:`~geometry.mesh.Mesh2D` stores cell edges, material IDs,
   coordinate system, and **boundary condition declarations**.
   Each face carries an optional :class:`~geometry.mesh.BC` field
   (``bc_left``/``bc_right`` for 1-D;
   ``bc_xmin``/``bc_xmax``/``bc_ymin``/``bc_ymax`` for 2-D).
   When ``None`` (the default), the solver applies its own default
   --- for the SN solver, that default is reflective.
   See :ref:`boundary-conditions` for details.

2. **Augmented geometry** --- :class:`SNMesh` wraps the base mesh and
   an angular quadrature, precomputing the coordinate-specific streaming
   stencil.  It also **resolves boundary conditions**: each ``BC`` tag
   on the mesh is looked up in :attr:`SNMesh.BOUNDARY_OPERATOR_REGISTRY` and converted
   to a validated kind string (``"vacuum"`` or ``"reflective"``)
   stored as ``sn_mesh.bc_left``, ``sn_mesh.bc_right``, etc.
   The sweep reads these resolved strings directly --- it never
   inspects the raw :class:`~geometry.mesh.BC` objects.  Precomputed
   stencil contents per coordinate system:

   - **Cartesian**: ``streaming_x[n,i] = 2|mu_x|/dx[i]`` and
     ``streaming_y[n,j] = 2|mu_y|/dy[j]`` --- the diamond-difference
     denominator terms, precomputed to avoid per-cell division in the
     sweep hot loop.
   - **Spherical**: ``face_areas`` (:math:`4\pi r^2`), ``delta_A``,
     ``alpha_half`` (angular redistribution dome),
     ``redist_dAw`` (:math:`\Delta A_i / w_n`, shape ``(nx, N)``),
     and ``tau_mm`` (Morel--Montry closure weights).
   - **Cylindrical**: ``face_areas`` (:math:`2\pi r`), ``delta_A``,
     ``alpha_per_level`` (per-level redistribution domes),
     ``redist_dAw_per_level`` (list of ``(nx, M)`` arrays), and
     ``tau_mm_per_level`` (per-level Morel--Montry weights).

3. **Solver** --- :func:`solve_sn` creates an ``SNMesh``, builds the
   ``SNSolver``, and runs power iteration.

.. code-block:: text

   Mesh1D / Mesh2D (base geometry + BC declarations)
       |
       v
   SNMesh (stencil + quadrature + alpha coefficients + resolved BCs)
       |
       v
   solve_sn() --> SNResult

Quadrature Dispatch
-------------------

The sweep dispatcher in :func:`transport_sweep` routes based on the
``SNMesh.curvature`` attribute and the quadrature type.  Boundary
conditions are **not** passed as a parameter to the sweep --- the
sweep reads the resolved BC kind strings directly from
``sn_mesh.bc_left``, ``sn_mesh.bc_right``, etc.:

.. code-block:: python

   if sn_mesh.curvature == "spherical":
       return _sweep_1d_spherical(...)
   elif sn_mesh.curvature == "cylindrical":
       return _sweep_1d_cylindrical(...)
   elif is_gl_1d:            # ny=1, mu_y=0, no aniso source
       return _sweep_1d_cumprod(...)
   else:
       return _sweep_2d_wavefront(...)

For 1D meshes (``ny=1``):

- **Gauss--Legendre** quadrature takes the fast cumprod path (all
  :math:`\mu_y = 0`, so no y-streaming).
- **Lebedev** quadrature falls through to the 2D wavefront sweep.
  Ordinates with :math:`\mu_x \neq 0` stream along *x*; the
  *y*-streaming terms cancel via reflective BCs on the single-cell
  *y*-dimension.  Ordinates with :math:`\mu_x = \mu_y = 0`
  (z-directed) reduce to pure collision:
  :math:`\psi = Q \cdot w_{\text{norm}} / \Sigt{}`.

Both quadratures recover the analytical eigenvalue exactly on
homogeneous problems (verified to machine precision for 1G/2G/4G).


The Transport Equation
======================

Cartesian 1D (Slab)
--------------------

The steady-state transport equation in a 1D slab:

.. math::
   :label: transport-cartesian

   \mu \frac{\partial \psi(x, \mu)}{\partial x}
   + \Sigt{} \, \psi(x, \mu)
   = \frac{Q}{W}

where :math:`\mu = \cos\theta` is the direction cosine, :math:`Q` is the
total isotropic source (fission + scattering), and :math:`W = \sum_n w_n`
is the quadrature weight sum.

Cartesian 2D
--------------

In two Cartesian dimensions the angular flux depends on two direction
cosines :math:`\mu_x` and :math:`\mu_y`:

.. math::
   :label: transport-cartesian-2d

   \mu_x \frac{\partial \psi}{\partial x}
   + \mu_y \frac{\partial \psi}{\partial y}
   + \Sigt{} \, \psi
   = \frac{Q}{W}

There is no angular coupling between ordinates --- each direction is
solved independently.  The two streaming terms are the only difference
from the 1D case.

Spherical 1D
-------------

In spherical coordinates the transport equation acquires an **angular
redistribution term** that couples ordinates:

.. math::
   :label: transport-spherical

   \mu \frac{\partial \psi}{\partial r}
   + \frac{1 - \mu^2}{r} \frac{\partial \psi}{\partial \mu}
   + \Sigt{} \psi = \frac{Q}{W}

The curvature term :math:`(1 - \mu^2)/r \cdot \partial\psi/\partial\mu`
arises because a neutron streaming radially at angle :math:`\mu` *rotates*
its direction cosine as it moves to a different radius.  Discretising this
term requires diamond difference in **both space and angle**.

Cylindrical 1D
---------------

For an infinitely long cylinder with azimuthal symmetry, the transport
equation in the radial variable :math:`r` is:

.. math::
   :label: transport-cylindrical

   \frac{\eta}{r} \frac{\partial(r\psi)}{\partial r}
   - \frac{1}{r} \frac{\partial(\xi\psi)}{\partial\varphi}
   + \Sigt{} \psi = \frac{Q}{W}

where the direction cosines are:

- :math:`\eta = \sin\theta\cos\varphi` --- radial projection (streaming)
- :math:`\xi = \sin\theta\sin\varphi` --- azimuthal component
- :math:`\mu = \cos\theta` --- axial component

The constraint :math:`\eta^2 + \xi^2 + \mu^2 = 1` holds.  The azimuthal
redistribution :math:`-\partial(\xi\psi)/\partial\varphi` couples ordinates
on each :math:`\mu`-level.

Multi-Group Extension
---------------------

For :math:`G` energy groups, each transport equation becomes a coupled
system with scattering transfer :math:`\Sigs{g' \to g}` between groups:

.. math::
   :label: multigroup

   \text{streaming} + \Sigt{g} \psi_g
   = \frac{1}{W} \left[
       \sum_{g'} \Sigs{g' \to g} \phi_{g'}
       + \frac{\chi_g}{k} \sum_{g'} \nSigf{g'} \phi_{g'}
   \right]

where the streaming operator depends on the coordinate system and
:math:`\phi_g = \sum_n w_n \psi_{g,n}` is the scalar flux.


.. _quadrature-types:

Angular Quadratures
===================

ORPHEUS provides four angular quadrature types for different geometries.

Gauss--Legendre (1D)
--------------------

For 1D slab and spherical geometry: :math:`N` points on
:math:`\mu \in [-1, 1]`, weights sum to 2.  Optimal for polynomial
integrands (degree :math:`2N-1` exact).  Also used for spherical 1D,
where the single direction cosine :math:`\mu` suffices for the angular
redistribution.

Implemented in :class:`GaussLegendre1D`.

Lebedev (Sphere)
-----------------

For 2D/3D Cartesian geometry: :math:`N` points on the unit sphere with
octahedral symmetry [Lebedev1999]_.  Weights sum to :math:`4\pi`.  On a
1D mesh, z-directed ordinates (:math:`\mu_x = \mu_y = 0`) are handled
as pure collision; all others stream along *x* with *y*-terms cancelling
via reflective BCs.

Implemented in :class:`LebedevSphere`.

Level-Symmetric S\ :sub:`N`
----------------------------

Standard triangular quadrature with :math:`N/2` distinct :math:`\mu_z`
values per hemisphere.  Ordinates on each level are permutations of the
direction cosine set satisfying :math:`\eta^2 + \xi^2 + \mu^2 = 1`.
Equal spacing in :math:`\mu^2` is used with :math:`\mu_1^2 = 4/(N(N+2))`
[CarlsonLathrop1965]_.

Weights sum to :math:`4\pi`.  Provides the ``level_indices`` structure
needed by the cylindrical sweep.  Unlike :class:`ProductQuadrature`
(which has one level per :math:`\mu_z` value), the Level-Symmetric
quadrature groups both :math:`+\mu_z` and :math:`-\mu_z` hemispheres
on the same level (grouped by :math:`|\mu_z|`).  Within each level,
ordinates are sorted by increasing :math:`\eta` for the azimuthal sweep.

Implemented in :class:`LevelSymmetricSN`.

Product Quadrature (GL x equispaced)
-------------------------------------

Tensor product of Gauss--Legendre in :math:`\mu = \cos\theta` (polar)
and equispaced points in :math:`\varphi` (azimuthal).  Each :math:`\mu`
level has the same number of azimuthal points, giving a clean level
structure ideal for the cylindrical sweep.  Weights:

.. math::

   w_{p,m} = w_{\text{GL}}(\mu_p) \cdot \frac{2\pi}{N_\varphi}

Sum to :math:`4\pi`.  Within each level, ordinates are sorted by
increasing :math:`\eta = \sin\theta\cos\varphi` to match the
:math:`\alpha` recursion convention from [Bailey2009]_ Eq. 50.

Implemented in :class:`ProductQuadrature`.

Reflection Index
-----------------

Each quadrature implements a :meth:`reflection_index` method that
returns an index array mapping each ordinate :math:`n` to its
**mirror image** :math:`n'` obtained by negating the direction cosine
along a specified axis.  For example, ``reflection_index("x")``
finds the ordinate whose direction cosines match :math:`(-\mu_x, \mu_y, \mu_z)`.

The implementation in :func:`_find_reflections` computes the
Euclidean distance between the target direction (with one component
negated) and all ordinate directions, then returns the closest match:

.. math::

   n' = \arg\min_j \bigl[
       (\mu_{x,j} - (-\mu_{x,n}))^2
       + (\mu_{y,j} - \mu_{y,n})^2
       + (\mu_{z,j} - \mu_{z,n})^2
   \bigr]

For Gauss--Legendre (1D), the reflection in *x* is simply
:math:`n' = N - 1 - n` because the GL points are symmetric about
zero.  Reflection in *y* is the identity since :math:`\mu_y = 0`.

For multi-dimensional quadratures (Lebedev, Level-Symmetric, Product),
the reflection indices are precomputed at construction time for all
three axes (*x*, *y*, *z*) and stored as ``_ref_x``, ``_ref_y``,
``_ref_z``.  These indices are used by the sweep to implement
reflective boundary conditions (see :ref:`boundary-conditions`).

Comparison Table
-----------------

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 15 20

   * - Quadrature
     - Geometry
     - :math:`\sum w`
     - Level structure
     - Best for
   * - Gauss--Legendre
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
     - Sphere, Cylinder
     - :math:`4\pi`
     - Yes
     - Curvilinear
   * - Product
     - Cylinder
     - :math:`4\pi`
     - Yes
     - Cylindrical 1D


The Discrete Balance Equation
=============================

This is the core of the S\ :sub:`N` method.  The balance equations are
presented from simplest to most complex: Cartesian geometries have no
angular redistribution; curvilinear geometries add :math:`\alpha` coupling
and a geometry factor :math:`\Delta A/w`.

.. _balance-cartesian-1d:

Cartesian 1D Balance Equation
------------------------------

Integrating :eq:`transport-cartesian` over a spatial cell
:math:`[x_{i-1/2}, x_{i+1/2}]` of width :math:`\Delta x_i` and applying
the divergence theorem to the streaming term:

.. math::

   \mu_n \bigl[\psi_{i+\frac12} - \psi_{i-\frac12}\bigr]
   + \Sigt{} \Delta x_i\, \psi_{n,i} = S_i \Delta x_i

where :math:`S_i = Q_i / W` and face areas are unity in slab geometry.
Applying the diamond-difference closure
:math:`\psi_{n,i} = \frac{1}{2}(\psi_{\rm in} + \psi_{\rm out})` and
:math:`\psi_{\rm out} = 2\psi_{n,i} - \psi_{\rm in}`, we solve for the
cell-average angular flux:

.. math::
   :label: dd-cartesian-1d

   \psi_{n,i}
   = \frac{S_i + \dfrac{2|\mu_n|}{\Delta x_i}\, \psi_{\rm in}}
          {\Sigt{} + \dfrac{2|\mu_n|}{\Delta x_i}}

This is the simplest balance equation: no :math:`\alpha` redistribution
and no :math:`\Delta A` factor, because slab geometry has no curvature.
The streaming coefficient :math:`2|\mu|/\Delta x` is precomputed by
:class:`SNMesh` as ``streaming_x[n, i]``.

.. _balance-cartesian-2d:

Cartesian 2D Balance Equation
-------------------------------

Integrating :eq:`transport-cartesian-2d` over a rectangular cell
:math:`\Delta x_i \times \Delta y_j`:

.. math::

   \mu_{x,n}\bigl[\psi_{i+\frac12,j} - \psi_{i-\frac12,j}\bigr] \Delta y_j
   + \mu_{y,n}\bigl[\psi_{i,j+\frac12} - \psi_{i,j-\frac12}\bigr] \Delta x_i
   + \Sigt{} \Delta x_i \Delta y_j\, \psi_{n,i,j}
   = S_{i,j}\, \Delta x_i \Delta y_j

Dividing through by :math:`\Delta x_i \Delta y_j` and applying
diamond-difference closures in **both** directions simultaneously:

.. math::

   \psi_{n,i} &= \tfrac{1}{2}(\psi^x_{\rm in} + \psi^x_{\rm out})
   \qquad\text{(x-closure)} \\
   \psi_{n,i} &= \tfrac{1}{2}(\psi^y_{\rm in} + \psi^y_{\rm out})
   \qquad\text{(y-closure)}

yields the 2D DD equation:

.. math::
   :label: dd-cartesian-2d

   \psi_{n,i,j}
   = \frac{S_{i,j}
     + s_x\, \psi^x_{\rm in}
     + s_y\, \psi^y_{\rm in}}
     {\Sigt{} + s_x + s_y}

where the streaming coefficients are:

.. math::

   s_x = \frac{2|\mu_{x,n}|}{\Delta x_i}, \qquad
   s_y = \frac{2|\mu_{y,n}|}{\Delta y_j}

Both outgoing face fluxes are then updated from the DD closure:

.. math::

   \psi^x_{\rm out} = 2\psi_{n,i,j} - \psi^x_{\rm in}, \qquad
   \psi^y_{\rm out} = 2\psi_{n,i,j} - \psi^y_{\rm in}

These are precomputed by :class:`SNMesh` as ``streaming_x[n, i]`` and
``streaming_y[n, j]``, so the inner loop in
:func:`_sweep_2d_wavefront` reduces to a single vectorised division per
diagonal.

.. _balance-curvilinear:

Curvilinear Balance Equation (Spherical and Cylindrical)
---------------------------------------------------------

Derivation from the Continuous PDE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Start with the general 1D curvilinear transport equation.  In
conservative form for a coordinate :math:`r` with face area
:math:`A(r)` and volume element :math:`V`:

.. math::
   :label: conservative-form

   \frac{\mu_n}{V_i}
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \frac{1}{V_i}
   \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
   + \Sigt{} \psi_{n,i} = S_i

.. vv-status: conservative-form documented

where the streaming cosine is :math:`\mu_n` for spherical and
:math:`\eta_m` for cylindrical, and :math:`S_i = Q_i / W` is the
isotropic source density divided by the quadrature weight sum.

**Step 1: Integrate the PDE over a spatial cell.**

For spherical geometry, integrating :eq:`transport-spherical` over the
shell :math:`[r_{i-1/2}, r_{i+1/2}]` and using the divergence theorem
on the radial streaming gives:

.. math::

   \mu_n \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \int_{V_i} \frac{1-\mu^2}{r} \frac{\partial\psi}{\partial\mu}\, dV
   + \Sigt{} V_i \psi_{n,i} = S_i V_i

For cylindrical geometry, integrating :eq:`transport-cylindrical` over
the annular shell gives:

.. math::

   \eta_m \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   - \int_{V_i} \frac{1}{r} \frac{\partial(\xi\psi)}{\partial\varphi}\, dV
   + \Sigt{} V_i \psi_{m,i} = S_i V_i

**Step 2: Discretise the angular redistribution.**

The angular integral is discretised as a finite difference in the
ordinate index.  For spherical:

.. math::

   \int_{V_i} \frac{1-\mu^2}{r}\frac{\partial\psi}{\partial\mu}\, dV
   \;\approx\;
   \alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}

For cylindrical (per :math:`\mu`-level):

.. math::

   -\int_{V_i} \frac{1}{r}\frac{\partial(\xi\psi)}{\partial\varphi}\, dV
   \;\approx\;
   \alpha_{m+\frac12}\psi_{m+\frac12} - \alpha_{m-\frac12}\psi_{m-\frac12}

**Step 3: Apply the geometry factor** :math:`\Delta A / w`.

The raw discretisation above does NOT preserve per-ordinate flat-flux
consistency.  The correct form from [Bailey2009]_ includes the geometry
factor :math:`\Delta A_i / w_n`:

.. math::
   :label: balance-general

   \mu_n
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \frac{\Delta A_i}{w_n}
   \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
   + \Sigt{} V_i \psi_{n,i} = S_i V_i

where :math:`\Delta A_i = A_{i+1/2} - A_{i-1/2}`.  This is
[Bailey2009]_ Eq. 7--10 for spherical and Eq. 50--55 for cylindrical.

Note why :eq:`dd-cartesian-1d` has no :math:`\alpha` or :math:`\Delta A`
terms: in Cartesian geometry the face area is unity (:math:`A = 1`), so
:math:`\Delta A = 0`, and there is no curvature to redistribute angular
flux.

The Alpha Redistribution Coefficients
--------------------------------------

The :math:`\alpha` coefficients encode how the angular flux redistributes
between neighbouring ordinates due to the geometry curvature.  They are
defined recursively:

.. math::
   :label: alpha-recursion

   \alpha_{n+\frac12} = \alpha_{n-\frac12} - w_n \mu_n

with the boundary condition :math:`\alpha_{1/2} = 0`.

For **spherical** geometry, all :math:`N` ordinates form a single
sequence sorted by :math:`\mu` (most negative to most positive).
The :math:`\alpha` values form a **non-negative dome**: they rise while
:math:`\mu < 0`, peak near :math:`\mu = 0`, and fall back to zero at
:math:`\mu = 1`.  The endpoint condition
:math:`\alpha_{N+1/2} = 0` is guaranteed by Gauss--Legendre
antisymmetry: :math:`\sum_n w_n \mu_n = 0`.

For **cylindrical** geometry, each :math:`\mu`-level has its own
independent :math:`\alpha` sequence.  On level :math:`p`, the ordinates
are sorted by increasing :math:`\eta` (radial direction cosine), and the
recursion uses :math:`\eta` instead of :math:`\mu`:

.. math::
   :label: alpha-cylindrical

   \alpha_{p,m+\frac12} = \alpha_{p,m-\frac12} - w_m \eta_m

This is [Bailey2009]_ Eq. 50.  Each level's :math:`\alpha` values form
an independent dome from :math:`\eta = -\sin\theta` to
:math:`\eta = +\sin\theta`.

**Dome shape properties:**

- :math:`\alpha_{n+1/2} \geq 0` for all :math:`n` (non-negative dome).
- The peak occurs near the ordinate where :math:`\mu_n` (or
  :math:`\eta_m`) crosses zero.
- The dome height scales with the quadrature weight sum: higher-order
  quadratures have narrower but taller domes.
- Non-negativity ensures the denominator of the DD equation is
  unconditionally positive, guaranteeing numerical stability.

The code stores these in ``SNMesh.alpha_half`` (spherical, shape
``(N+1,)``) and ``SNMesh.alpha_per_level`` (cylindrical, list of
``(M+1,)`` arrays).

The Geometry Factor and Why It Is Needed
-----------------------------------------

The geometry factor :math:`\Delta A_i / w_n` in :eq:`balance-general`
is the key to correct curvilinear transport.  Without it, the balance
equation violates **per-ordinate flat-flux consistency**: for a spatially
uniform, isotropic flux :math:`\psi = \text{const}`, the streaming and
redistribution terms should cancel exactly for EACH ordinate
individually.

**Proof of consistency.**

Set :math:`\psi_{n,i} = \psi_{n+1/2} = \psi_{n-1/2} = \psi_0` (flat
in both space and angle) and :math:`\psi_{i+1/2} = \psi_{i-1/2} = \psi_0`
(flat in space).  The streaming term becomes:

.. math::

   \mu_n \bigl[A_{i+\frac12} - A_{i-\frac12}\bigr] \psi_0
   = \mu_n \,\Delta A_i\, \psi_0

The redistribution term with the :math:`\Delta A/w` factor becomes:

.. math::

   \frac{\Delta A_i}{w_n}
   \bigl[\alpha_{n+\frac12} - \alpha_{n-\frac12}\bigr] \psi_0
   = \frac{\Delta A_i}{w_n} (-w_n \mu_n) \psi_0
   = -\mu_n \,\Delta A_i\, \psi_0

where we used the recursion :eq:`alpha-recursion`:
:math:`\alpha_{n+1/2} - \alpha_{n-1/2} = -w_n \mu_n`.  The two terms
cancel exactly, giving :math:`\Sigt{} \psi_0 = S_0`, which is the
correct homogeneous solution.

**Without** the :math:`\Delta A/w` factor (i.e., using
:math:`[\alpha_{n+1/2}\psi_{n+1/2} - \alpha_{n-1/2}\psi_{n-1/2}]`
directly), the redistribution term for flat flux is
:math:`(-w_n \mu_n)\psi_0`, but the streaming term is
:math:`\mu_n \Delta A_i \psi_0`.  These differ by a factor of
:math:`\Delta A_i`, so consistency only holds in the limit
:math:`\Delta A_i \to 0` (i.e., at the origin or on an infinitely fine
mesh).

**Consequence of the missing factor:**  The solver creates artificial
angular anisotropy that *worsens* with mesh refinement near :math:`r = 0`
(where :math:`\Delta A_i` is smallest but non-zero).  This manifests as
a flux spike at the origin in fixed-source problems and as divergent
eigenvalues in heterogeneous eigenvalue problems.

The code precomputes this factor as ``SNMesh.redist_dAw`` (spherical,
shape ``(nx, N)``) and ``SNMesh.redist_dAw_per_level`` (cylindrical,
list of ``(nx, M)`` arrays).

The Morel--Montry Flux Dip
----------------------------

Even with the correct :math:`\Delta A/w` factor, the standard
diamond-difference closure (equal weight :math:`\tau = 0.5`) introduces
a flux error near :math:`r = 0` known as the **Morel--Montry flux dip**
[MorelMontry1984]_.

The standard DD angular closure is:

.. math::

   \psi_{n,i} = \frac{1}{2}(\psi_{n-\frac12} + \psi_{n+\frac12})

This can be rewritten as:

.. math::

   \psi_{n+\frac12} = 2\psi_{n,i} - \psi_{n-\frac12}

The contamination factor :math:`\beta` ([Bailey2009]_ Eq. 41) quantifies
the coupling between the leading-order scalar flux and the first-order
current in the asymptotic diffusion limit.  For spherical geometry:

.. math::

   \beta = \frac{1}{2} \sum_{n=1}^{N} \mu_n
   \bigl[\alpha_{n+\frac12}\, \mu_{n+\frac12}
        - \alpha_{n-\frac12}\, \mu_{n-\frac12}\bigr]

where :math:`\mu_{n\pm 1/2}` are the angular cell-edge cosines.  For
cylindrical, the equivalent is a per-level sum using :math:`\eta` and
:math:`\eta_{m\pm 1/2}`.  When :math:`\beta \neq 0`, the discrete
S\ :sub:`N` equations satisfy a **contaminated** diffusion equation near
:math:`r = 0`, producing the artificial flux dip (or spike).

The module :mod:`derivations.sn_contamination` computes :math:`\beta`
for any quadrature and geometry.  With the correct :math:`\Delta A/w`
factor AND Morel--Montry weights, :math:`\beta \sim 10^{-16}`
(machine zero) for both spherical and cylindrical.

Weighted Diamond Difference (WDD) and Morel--Montry Weights
-------------------------------------------------------------

The Morel--Montry (M-M) angular closure replaces the equal-weight DD
with position-dependent weights :math:`\tau_n` [Bailey2009]_ Eq. 74:

.. math::
   :label: wdd-closure

   \psi_{n,i} = (1 - \tau_n)\,\psi_{n-\frac12} + \tau_n\,\psi_{n+\frac12}

Solving for the angular face flux:

.. math::
   :label: wdd-face

   \psi_{n+\frac12}
   = \frac{\psi_{n,i} - (1 - \tau_n)\,\psi_{n-\frac12}}{\tau_n}

The M-M weights are defined as:

.. math::
   :label: mm-weights

   \tau_n = \frac{\mu_n - \mu_{n-\frac12}}{\mu_{n+\frac12} - \mu_{n-\frac12}}

where :math:`\mu_{n\pm 1/2}` are the angular cell edges.

**Spherical cell edges:**  :math:`\mu_{1/2} = -1`,
:math:`\mu_{N+1/2} = +1`, and interior edges by weight-sum:
:math:`\mu_{n+1/2} = \mu_{n-1/2} + w_n`.  This is exact for
Gauss--Legendre quadrature because the weights correspond to the
:math:`\mu`-space widths of the angular cells.

**Cylindrical cell edges:**  :math:`\eta_{1/2} = -\sin\theta`,
:math:`\eta_{M+1/2} = +\sin\theta`, and interior edges at
**midpoints** of consecutive :math:`\eta` values:
:math:`\eta_{m+1/2} = (\eta_m + \eta_{m+1})/2`.
The weight-sum approach is NOT used for cylindrical because the
quadrature weights are uniform in :math:`\varphi`-space (not
:math:`\eta`-space): the Product quadrature spaces :math:`\varphi`
equally, but :math:`\eta = \sin\theta\cos\varphi` is
cosine-distributed, so equal :math:`\varphi`-widths map to unequal
:math:`\eta`-widths.  The midpoint approach gives a proper partition
of :math:`[-\sin\theta, +\sin\theta]`.

For the Product quadrature with equally-spaced :math:`\varphi`,
ordinates come in **pairs** with the same :math:`|\eta|` but opposite
:math:`\xi` (e.g., :math:`\varphi = \pi/4` and :math:`\varphi = 7\pi/4`
both give :math:`\eta = \sin\theta/\sqrt{2}`).  The midpoint between
paired ordinates equals their shared :math:`\eta`, creating zero-width
angular cells.  The resulting :math:`\tau` alternates between 0.5
(DD, for the left member of each pair) and 1.0 (step, for the right
member).  This alternating pattern is correct but could be smoothed by
using a Gauss-type azimuthal quadrature with distinct :math:`\eta`
values (see `GitHub Issue #1 <https://github.com/deOliveira-R/ORPHEUS/issues/1>`_).

The M-M weights force the contamination factor :math:`\beta` to **machine
zero** (verified: :math:`\beta \sim 10^{-16}`), completely eliminating
the Morel--Montry flux dip.

**Clipping:** The code clips :math:`\tau_n` to :math:`[0.5, 1.0]` to
prevent negative angular face fluxes.  This preserves the stability of
the standard DD while gaining the accuracy of the M-M correction.

The code stores these as ``SNMesh.tau_mm`` (spherical, shape ``(N,)``)
and ``SNMesh.tau_mm_per_level`` (cylindrical, list of ``(M,)`` arrays).

Substituting the WDD Closure into the Balance Equation
-------------------------------------------------------

Combining the balance equation :eq:`balance-general` with the WDD
angular closure :eq:`wdd-closure` and the standard spatial DD
(:math:`\psi_{n,i} = \frac{1}{2}(\psi_{\rm in}^s + \psi_{\rm out}^s)`,
:math:`\psi_{\rm out}^s = 2\psi_{n,i} - \psi_{\rm in}^s`), define:

.. math::

   c_{\rm out} &= \frac{\alpha_{n+\frac12}}{\tau_n} \\[6pt]
   c_{\rm in}  &= \frac{(1-\tau_n)}{\tau_n}\,\alpha_{n+\frac12}
                 + \alpha_{n-\frac12}

The cell-average angular flux is then:

.. math::
   :label: dd-solve

   \psi_{n,i} = \frac{
       S_i V_i
       + |\mu_n|(A_{\rm in} + A_{\rm out})\,\psi_{\rm in}^s
       + \dfrac{\Delta A_i}{w_n}\, c_{\rm in}\, \psi_{n-\frac12}
   }{
       2|\mu_n|\, A_{\rm out}^s
       + \dfrac{\Delta A_i}{w_n}\, c_{\rm out}
       + \Sigt{} V_i
   }

where the superscript :math:`s` denotes spatial face fluxes, and
:math:`A_{\rm in}`, :math:`A_{\rm out}` are the cell face areas in the
direction of neutron travel (see :ref:`sweep-algorithm` below for their
definition).  This is the equation solved by both
:func:`_sweep_1d_spherical` and :func:`_sweep_1d_cylindrical`.

Geometry Comparison
--------------------

.. list-table::
   :header-rows: 1
   :widths: 15 28 28 29

   * - Aspect
     - Cartesian
     - Spherical
     - Cylindrical
   * - Streaming cosine
     - :math:`\mu`
     - :math:`\mu`
     - :math:`\eta` (radial)
   * - Face area :math:`A`
     - 1 (per unit area)
     - :math:`4\pi r^2`
     - :math:`2\pi r`
   * - Volume :math:`V`
     - :math:`\Delta x`
     - :math:`\tfrac{4}{3}\pi(r_{\rm out}^3 - r_{\rm in}^3)`
     - :math:`\pi(r_{\rm out}^2 - r_{\rm in}^2)`
   * - :math:`\Delta A`
     - 0
     - :math:`4\pi(r_{\rm out}^2 - r_{\rm in}^2)`
     - :math:`2\pi(r_{\rm out} - r_{\rm in})`
   * - Redistribution
     - None
     - :math:`+(\Delta A/w)\,[\alpha\psi]`
     - :math:`+(\Delta A/w)\,[\alpha\psi]`
   * - :math:`\alpha` scope
     - N/A
     - Global (all :math:`N` ordinates)
     - Per :math:`\mu`-level
   * - :math:`\alpha` recursion variable
     - N/A
     - :math:`\mu`
     - :math:`\eta`
   * - Quadrature required
     - GL or Lebedev
     - GL
     - Product or Level-Sym


.. _cell-update-strategies:

Cell update strategies (the strategy contract)
==============================================

The discrete balance equation derived above (slab DD
:eq:`dd-cartesian-1d`, the M-M-closed curvilinear update) yields, for
each cell, a small algebraic system: combine the upstream face flux
(and, for sphere/cylinder, the upstream angular half-flux) with a
local source and the cell's total cross section to produce the
cell-average flux plus the downstream states.  The closure relating
:math:`\overline{\psi}_i` to :math:`\psi_{i-1/2}` and
:math:`\psi_{i+1/2}` is **not unique** — Diamond Difference (DD),
weighted DD, Linear Discontinuous (LD), Step, and Exponential
Characteristic (EC) are all valid choices, each with different
truncation error, positivity, and cost.  Per Cardinal Rule 2
(architecture), the cell-update math is **the same algebra** in slab,
sphere, and cylindrical 1-D — only the populated fields of the
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` packet
change.  Lifting the closure into a strategy contract makes the
sweep driver thin and lets each closure be unit-tested in isolation.

The strategy contract is owned by
:mod:`orpheus.sn.spatial.cell_update`.

The Protocol
------------

The contract is a ``@runtime_checkable`` ``typing.Protocol`` —
satisfied by structural typing, not inheritance — exposing two class-
level traits and a single :meth:`update` method:

* :class:`~orpheus.sn.spatial.cell_update.CellUpdate`

  - ``is_linear: bool`` — whether the closure is linear in its inputs.
    Diamond Difference is linear; Step's positivity-fixup, weighted-DD
    with a flux-dependent weight, and EC with a clipped argument are
    not.
  - ``is_positivity_preserving: bool`` — whether non-negative inputs
    guarantee non-negative outputs.  Diamond Difference is **not**
    positivity preserving (Lewis & Miller §5.3, where DD's tendency
    to produce negative cell-edge fluxes is exhibited and motivates
    the choice of Step or weighted-DD in stiff cells); Step is.
  - ``update(visit, total_xs, source, upstream_state) ->
    CellResult`` — the cell update itself.  ``visit`` is a
    :class:`~orpheus.sn.spatial.cell_update.CellVisit` packet (see
    next subsection) that combines the geometric
    :class:`~orpheus.geometry.reduced_operator.StreamingTerms` with
    sweep-direction-resolved data.

The two helper dataclasses (frozen, slotted) carry the per-cell
state:

* :class:`~orpheus.sn.spatial.cell_update.UpstreamState`

  - ``spatial_upstream: np.ndarray`` — shape ``(ng,)``.  Face flux
    entering the cell from the upstream face (always populated).
  - ``angular_upstream: np.ndarray | None`` — shape ``(ng,)`` for
    sphere/cylinder; ``None`` for slab.  :math:`\psi_{n-1/2,\,i}`,
    the half-flux at the upstream half-angle.

* :class:`~orpheus.sn.spatial.cell_update.CellResult`

  - ``cell_average_flux: np.ndarray`` — shape ``(ng,)``.  The cell-
    average flux :math:`\overline{\psi}_i = \mathrm{numer}/\mathrm{denom}`
    from the closure's algebraic solve.
  - ``outgoing_spatial_flux: np.ndarray | None`` — shape ``(ng,)`` in
    the typical case; ``None`` for the cylindrical pure-azimuthal
    degenerate case where the cell has no radial face flow (see
    below).
  - ``outgoing_angular_state: np.ndarray | None`` — shape ``(ng,)``
    for sphere/cylinder; ``None`` for slab.  :math:`\psi_{n+1/2,\,i}`
    from the Morel--Montry closure.

The SN sweep DAG and ``CellVisit``
-----------------------------------

The SN sweep is a **topological sort of a directed cell graph**.
For a given ordinate :math:`\Omega_n`, every face :math:`f` of the
mesh is oriented by the sign of :math:`\Omega_n \cdot \hat n_f` — an
edge from cell :math:`A` to cell :math:`B` if :math:`\Omega_n` points
from :math:`A` into :math:`B` across that face.  The sweep walks
cells in a topological order over this DAG so that, when each cell is
visited, all its upstream face fluxes are already known.  This is
the SN-specific graph-theoretic concept; MoC uses a different
mathematical structure (fiber bundles + solution sheaves over
characteristic curves), and CP / diffusion / MC have no sweep at
all.  Per Cardinal Rule 2 (architecture), no shared
``SweepGraph`` Protocol is hoisted across solvers — the sweep DAG
lives in :mod:`orpheus.sn`.

The contract's :meth:`update` consumes a
:class:`~orpheus.sn.spatial.cell_update.CellVisit` packet rather
than a raw
:class:`~orpheus.geometry.reduced_operator.StreamingTerms`.
The :class:`CellVisit` composes:

* ``cell_idx: int`` — the cell being visited.
* ``streaming_terms: StreamingTerms`` — the **purely geometric**
  primitive (``face_area_inner`` / ``face_area_outer`` are
  geometric labels — inner = closer to :math:`r=0`, outer =
  farther — independent of sweep direction).
* ``face_area_downstream: float | None`` — **sweep-direction-
  resolved**.  For an outward sphere or cylinder sweep
  (:math:`\mu \ge 0`) it equals ``streaming_terms.face_area_outer``;
  for an inward sweep (:math:`\mu < 0`) it equals
  ``streaming_terms.face_area_inner``.  ``None`` for slab (slab DD
  does not read face areas) and for the cylindrical pure-azimuthal
  degenerate case (no spatial flow).

The :class:`CellVisit` packets are produced by
:meth:`SNMesh.iter_cell_visits(ordinate_idx, mu_level_idx=None)
<orpheus.sn.geometry.SNMesh.iter_cell_visits>` — a generator that
yields cells in DAG-topological order for the given ordinate.  The
method encapsulates the inward / outward branching, the cylindrical
per-:math:`\mu`-level traversal, and the pure-azimuthal degenerate
handling.  The sweep at :mod:`orpheus.sn.sweep` consumes this
generator::

    for visit in sn_mesh.iter_cell_visits(ordinate_idx=n):
        upstream = UpstreamState(
            spatial_upstream=psi_face,
            angular_upstream=psi_angle[visit.cell_idx],
        )
        result = cell_update.update(
            visit, total_xs, source, upstream,
        )
        ...

The cell-update strategy receives only **resolved** data — no
sign-of-:math:`\mu` branching inside the strategy.  This pattern
moves the graph-theoretic concept to where it belongs (the SN
module) and keeps the geometry-layer
:class:`~orpheus.geometry.reduced_operator.StreamingTerms`
geometry-only and reusable by future MoC / CP / diffusion modules
that have different mathematical structures.

Slab vs curvilinear discrimination
-----------------------------------

A strategy distinguishes slab from curvilinear by a single field test
on the visit's streaming terms:

* **Slab** — ``visit.streaming_terms.alpha_in is None`` (and the rest
  of the curvature bundle, ``alpha_out``, ``delta_A_over_w``,
  ``tau_mm``, ``face_area_inner`` / ``face_area_outer``, are all
  ``None``).  ``upstream_state.angular_upstream is None``.  The
  strategy returns ``CellResult(outgoing_angular_state=None, ...)``.
* **Sphere or cylinder** —
  ``visit.streaming_terms.alpha_in is not None``; the full curvature
  bundle is populated.  ``upstream_state.angular_upstream`` carries
  :math:`\psi_{n-1/2,\,i}`.  The strategy returns the M-M-closed
  ``outgoing_angular_state``.

This single-field discrimination convention is locked in by
foundation-tier protocol-conformance tests in
``tests/sn/spatial/test_cell_update_protocol.py``; concrete
strategies and the Wave D sweep rewrite read this same field to
dispatch.

Cylindrical pure-azimuthal degenerate case
-------------------------------------------

For cylindrical 1-D radial sweeps with a product or level-symmetric
quadrature, ordinates with axial direction cosine
:math:`|\mu_z| \to 1` have radial direction cosine
:math:`|\eta| = \sqrt{1 - \mu_z^2} \to 0`.  In this limit the cell
has **no radial face flow** — the streaming term
:math:`\mu_x \cdot \partial_r` vanishes — and the cell-update
algebra collapses to the redistribution-only form

.. math::

   \mathrm{denom} = (\Delta A / w)\,c_{\rm out} + \Sigma_t\,V_i,
   \qquad
   \mathrm{numer} = Q_i\,V_i + (\Delta A/w)\,c_{\rm in}\,\psi_{n-1/2,\,i},

with no spatial-flux contribution.  The strategy contract signals
this case by setting ``CellResult.outgoing_spatial_flux = None``: the
sweep driver, on receiving ``None``, skips the face-flux update for
that cell.  The angular M-M closure remains active — angular
redistribution physics is still present.

The numerical threshold is ``streaming_terms.abs_mu < 1e-15``, with
``abs_mu`` populated from the **global ordinate**
:math:`|\eta|` on the streaming-terms packet (resolved through
``level_indices`` for cylindrical geometry — see
:doc:`structured_geometry`, "Connection coefficients (reduced
streaming operator)").  In this case
:meth:`SNMesh.iter_cell_visits
<orpheus.sn.geometry.SNMesh.iter_cell_visits>` yields visits with
``face_area_downstream = None`` to signal "no spatial flow" to the
strategy.

The DD recurrence
------------------

For non-degenerate cells, the closure relation reduces — for slab
geometry, the cell-update math is the DD recurrence
:eq:`dd-recurrence` (see :ref:`sweep-cumprod` below); for curvilinear
geometry, it is the curvilinear DD form combining the
:math:`\Delta A/w` redistribution with the M-M angular closure.  The
sweep driver inlines this math today; Wave D (Issue #159) will
rewrite the driver to dispatch through a
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` strategy.

The first concrete strategy — :class:`DiamondDifference` — is shipped
in Round 2 of Wave C (Issue #158) as a bit-identical extraction of
the existing inlined sweep math.  Linear Discontinuous (Lewis &
Miller §5.3 — preview), Exponential Characteristic, and Step
strategies are deferred to a Wave C-extension session, each with its
own MMS spatial-convergence verification.

Diamond Difference
------------------

The first concrete strategy is
:class:`~orpheus.sn.spatial.diamond.DiamondDifference`
(:mod:`orpheus.sn.spatial.diamond`).  It implements the **same**
algebra as the existing inlined sweep — Round 2 of Wave C is a
bit-identical extraction, gated by ``np.array_equal`` hand-calc tests
in ``tests/sn/spatial/test_diamond.py`` against the sweep's scalar
formulas at :mod:`orpheus.sn.sweep`.

Per Wave C decision **D5** (one geometry-polymorphic class), the
strategy is a single :class:`DiamondDifference` that handles slab,
sphere, and cylinder by branching on two
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` fields:
``alpha_in is None`` (slab vs curvilinear) and ``abs_mu < 1e-15``
(cylindrical pure-azimuthal degenerate vs not).

**Slab branch** (``streaming_terms.alpha_in is None``).  The flat /
Cartesian DD recurrence reduces to the per-cell scalar form of
:eq:`dd-recurrence`:

.. math::
   :label: dd-slab-scalar

   \overline{\psi}_i \;=\; \tfrac12\bigl(\psi_{i-1/2}
                                         + \psi_{i+1/2}\bigr),
   \qquad
   \psi_{i+1/2} \;=\; \frac{2|\mu_n| - \Delta x_i\,\Sigma_t}
                              {2|\mu_n| + \Delta x_i\,\Sigma_t}\,
                       \psi_{i-1/2}
                   \;+\; \frac{2\,Q_i\,\Delta x_i / W}
                              {2|\mu_n| + \Delta x_i\,\Sigma_t},

with ``W = Σ_n w_n`` the quadrature weight sum, mirroring the
sweep's vectorised cumprod path
(:func:`orpheus.sn.sweep._sweep_1d_cumprod` lines 117–123) and per-
cell solver (:func:`orpheus.sn.sweep._solve_recurrence` lines 208–
222) at the operation level.  Per the strategy contract, ``source``
arrives at the cell update **already weight-normalised** by the
sweep — for slab, ``source = Q · Δx / W`` (and slab cell volume is
``V = Δx``).  For slab, the strategy sets
:attr:`~orpheus.sn.spatial.cell_update.CellResult.outgoing_angular_state`
to ``None`` — slab geometry has no angular redistribution.

**Curvilinear branch** (``streaming_terms.alpha_in is not None`` and
``abs_mu ≥ 1e-15``).  Sphere or cylinder, away from the cylindrical
pure-azimuthal degenerate case.  The strategy couples the M-M
angular closure :eq:`mm-weights` to the WDD spatial closure
:eq:`wdd-closure`, with the redistribution constants

.. math::
   :label: dd-mm-closure-constants

   c_{\rm out} \;=\; \alpha_{n+\tfrac12}/\tau_n,
   \qquad
   c_{\rm in}  \;=\; \frac{1 - \tau_n}{\tau_n}\,\alpha_{n+\tfrac12}
                       + \alpha_{n-\tfrac12},

built from the Bailey 2009 dome :eq:`alpha-recursion` and the
Morel–Montry weight :eq:`mm-weights`.  The cell-update is then

.. math::
   :label: dd-curvilinear-scalar

   \overline{\psi}_{n,i} \;=\;
       \frac{Q_i\,V_i / W
             + |\mu_n|\,(A_{i-1/2} + A_{i+1/2})\,\psi^s_{n,\,{\rm in}}
             + (\Delta A_i / w_n)\,c_{\rm in}\,\psi_{n-\tfrac12,\,i}}
            {2|\mu_n|\,A^s_{\rm out}
             + (\Delta A_i / w_n)\,c_{\rm out}
             + \Sigma_t\,V_i},

mirroring :func:`orpheus.sn.sweep._sweep_1d_spherical` lines
350–355 (and the structurally identical cylindrical branches at
sweep.py:511–531 / sweep.py:548–575) verbatim, with closures

.. math::

   \psi^s_{\rm out} \;=\; 2\overline{\psi}_{n,i}
                        - \psi^s_{n,\,{\rm in}},
   \qquad
   \psi_{n+\tfrac12,\,i} \;=\;
       (\overline{\psi}_{n,i}
         - (1 - \tau_n)\,\psi_{n-\tfrac12,\,i})/\tau_n.

**Cylindrical pure-azimuthal degenerate branch**
(``streaming_terms.alpha_in is not None`` and ``abs_mu < 1e-15``).
For a level whose axial direction cosine :math:`|\mu_z| \to 1`, the
radial direction cosine :math:`|\eta| \to 0` and the cell has no
radial face flow — the :math:`2|\mu| A_{\rm out}` and
:math:`|\mu|(A_{\rm in} + A_{\rm out})\,\psi^s_{\rm in}`
contributions drop out:

.. math::

   \mathrm{denom} \;=\; (\Delta A / w)\,c_{\rm out} + \Sigma_t\,V_i,
   \qquad
   \mathrm{numer} \;=\; Q_i\,V_i / W
                       + (\Delta A / w)\,c_{\rm in}\,
                          \psi_{n-\tfrac12,\,i},

mirroring :func:`orpheus.sn.sweep._sweep_1d_cylindrical` lines
533–543 verbatim.  The strategy returns
:attr:`~orpheus.sn.spatial.cell_update.CellResult.outgoing_spatial_flux`
``= None`` to signal "no face-flux write" to the sweep driver; the
M-M angular closure remains active.

**Traits and forward references.**  Diamond Difference has

* :attr:`~orpheus.sn.spatial.diamond.DiamondDifference.is_linear`
  ``= True`` — the cell average and downstream states are affine
  combinations of ``source`` and ``upstream_state``;
* :attr:`~orpheus.sn.spatial.diamond.DiamondDifference.is_positivity_preserving`
  ``= False`` — Lewis & Miller §5.3 exhibits the canonical thin-
  cell / large-source counter-example where DD's
  :math:`\psi_{\rm out} = 2\overline{\psi} - \psi_{\rm in}` produces
  negative outgoing flux from positive inputs.

Wave C-extension and Wave D will ship
:class:`Step` (positivity-preserving, :math:`\mathcal{O}(\Delta x)`),
:class:`LinearDiscontinuous` (:math:`\mathcal{O}(\Delta x^2)`,
better robustness in optically-thick cells), and
:class:`ExponentialCharacteristic` (positivity-preserving by
construction) as alternatives, each with its own MMS spatial-
convergence verification.

References
----------

* Lewis, E. E., & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*.  §5.3 covers Diamond Difference, weighted-DD,
  Step, and Linear Discontinuous closures and their positivity /
  truncation properties; §4.5 covers the Morel--Montry angular
  closure used for :math:`\psi_{n+1/2,\,i}`.
* Bailey, T. S., Adams, M. L., Yang, B., & Zika, M. R. (2009).
  *A piecewise linear finite element discretization of the
  diffusion equation for arbitrary polyhedral grids*.
  JCP 227, 3738--3757.  Eq. 50 (dome recursion), Eq. 74
  (Morel--Montry).

See also
--------

* :mod:`orpheus.sn.spatial.cell_update` — the contract module.
* :mod:`orpheus.sn.spatial.diamond` — the
  :class:`~orpheus.sn.spatial.diamond.DiamondDifference` concrete
  strategy.
* :doc:`structured_geometry`, "Connection coefficients (reduced
  streaming operator)" — the upstream side of the contract: where the
  per-cell, per-direction streaming-terms packet is built.


.. _sweep-algorithm:

Sweep Algorithm
===============

Because each cell's outgoing flux becomes the next cell's incoming flux,
the equations must be solved in the direction of neutron travel --- this
is called a **transport sweep**.

.. _sweep-cumprod:

Cartesian 1D: Cumprod Recurrence
---------------------------------

For the 1D slab with Gauss--Legendre quadrature, the DD equation
:eq:`dd-cartesian-1d` defines a recurrence for the outgoing face flux:

.. math::
   :label: dd-recurrence

   \psi_{\rm out} = a_i\, \psi_{\rm in} + b_i

where the coefficients for cell :math:`i` are:

.. math::

   a_i = \frac{2|\mu_n|/\Delta x_i - \Sigt{}}
              {2|\mu_n|/\Delta x_i + \Sigt{}},
   \qquad
   b_i = \frac{S_i}
              {2|\mu_n|/\Delta x_i + \Sigt{}}

This arises from substituting the DD closure
:math:`\psi_{\rm out} = 2\psi_{\rm avg} - \psi_{\rm in}` into
:eq:`dd-cartesian-1d`.  The coefficient :math:`a_i` is the
**stream-to-collision ratio**: it controls how much incoming flux
propagates through cell :math:`i`.

Unrolling the recurrence :math:`\psi_{\rm out}^{(i)} = a_i\, \psi_{\rm out}^{(i-1)} + b_i`
gives a linear first-order relation that can be solved analytically
using **cumulative products**.  Define:

.. math::

   C_i = \prod_{k=0}^{i} a_k, \qquad
   R_i = \sum_{k=0}^{i} \frac{b_k}{C_k}

Then the incoming face flux at cell :math:`i+1` is:

.. math::

   \psi_{\rm in}^{(i+1)} = C_i \bigl(\psi_{\rm in}^{(0)} + R_i\bigr)

and the cell-average flux is :math:`\psi_{\rm avg}^{(i)} = \frac{1}{2}(\psi_{\rm in}^{(i)} + \psi_{\rm out}^{(i)})`.

The implementation in :func:`_sweep_1d_cumprod` computes :math:`C` and
:math:`R` via ``np.cumprod`` and ``np.cumsum``, giving an
:math:`O(N \cdot n_x)` **vectorised** sweep --- all spatial cells for a
given ordinate are resolved simultaneously in numpy array operations,
with no Python-level cell loop.  This typically runs in sub-millisecond
time for practical meshes.

Exploiting GL symmetry, only positive-:math:`\mu` ordinates are swept
forward; negative-:math:`\mu` ordinates are obtained by reversing the
cell array and sweeping with the same coefficients.

.. _sweep-wavefront:

Cartesian 2D: Anti-Diagonal Wavefront Sweep
---------------------------------------------

In 2D, the DD equation :eq:`dd-cartesian-2d` creates a data dependency:
cell :math:`(i, j)` requires incoming face fluxes from its upwind
neighbours in both :math:`x` and :math:`y`.  Cells along an
**anti-diagonal** :math:`i + j = k` are mutually independent because
they share no incoming faces, so they can be solved simultaneously.

**The four quadrant sweeps.**  Each ordinate has a sign pair
:math:`(\text{sgn}(\mu_x), \text{sgn}(\mu_y))` that determines the
sweep direction.  The four combinations define four sweep patterns:

.. list-table::
   :header-rows: 1
   :widths: 20 20 30 30

   * - :math:`\mu_x`
     - :math:`\mu_y`
     - *x*-direction
     - *y*-direction
   * - :math:`+`
     - :math:`+`
     - left :math:`\to` right
     - bottom :math:`\to` top
   * - :math:`-`
     - :math:`+`
     - right :math:`\to` left
     - bottom :math:`\to` top
   * - :math:`+`
     - :math:`-`
     - left :math:`\to` right
     - top :math:`\to` bottom
   * - :math:`-`
     - :math:`-`
     - right :math:`\to` left
     - top :math:`\to` bottom

For each direction pair, the sweep visits anti-diagonals
:math:`k = 0, 1, \ldots, n_x + n_y - 2`.  On diagonal :math:`k`, the
cells :math:`(i, j)` satisfying :math:`i + j = k` (in the swept index
space) are gathered into a numpy batch and solved with a single vectorised
evaluation of :eq:`dd-cartesian-2d`.

**Vectorisation within each diagonal.**  Each diagonal contains up to
:math:`\min(n_x, n_y)` cells.  The incoming face fluxes ``psi_in_x``
and ``psi_in_y`` are gathered by advanced indexing; the DD equation is
evaluated as one numpy operation; and the outgoing face fluxes are
scattered back.  There is no Python-level cell loop within a diagonal.

**Reflective BCs in 2D.**  At each boundary face, the incoming flux for
ordinate :math:`n` is set to the outgoing flux of its reflected partner.
For the left/right boundaries (*x*-reflection), the partner is
``ref_x[n]`` (negating :math:`\mu_x`); for the top/bottom boundaries
(*y*-reflection), the partner is ``ref_y[n]`` (negating :math:`\mu_y`).
The reflection indices are precomputed by the quadrature's
:meth:`reflection_index` method.

Implemented in :func:`_sweep_2d_wavefront`.

Curvilinear 1D: Sequential Ordinate Sweep
------------------------------------------

For spherical and cylindrical geometries, the angular redistribution
couples successive ordinates, preventing vectorisation across the
ordinate dimension.  The sweep proceeds cell-by-cell,
ordinate-by-ordinate:

**Spherical:** Ordinates are processed from most negative :math:`\mu` to
most positive (a single global sequence).  Negative-:math:`\mu` ordinates
sweep **inward** (outer boundary to centre); positive-:math:`\mu`
ordinates sweep **outward** (centre to outer boundary).

**Cylindrical:** For each :math:`\mu`-level, azimuthal ordinates are
processed from most-inward (:math:`\eta = -\sin\theta`) to most-outward
(:math:`\eta = +\sin\theta`).  Negative-:math:`\eta` ordinates sweep
inward; :math:`\eta \approx 0` ordinates have no radial streaming
(pure redistribution); positive-:math:`\eta` ordinates sweep outward.

At each cell, the sweep solves :eq:`dd-solve` for :math:`\psi_{n,i}`,
then updates:

1. **Spatial face flux:**
   :math:`\psi_{\rm out}^s = 2\psi_{n,i} - \psi_{\rm in}^s`
2. **Angular face flux:**
   :math:`\psi_{n+1/2} = (\psi_{n,i} - (1-\tau_n)\psi_{n-1/2})/\tau_n`
   (using M-M weights)
3. **Scalar flux accumulation:**
   :math:`\phi_i \mathrel{+}= w_n \psi_{n,i}`

The spatial face flux propagates to the next cell; the angular face flux
propagates to the next ordinate on the same cell.

Implemented in :func:`_sweep_1d_spherical` and
:func:`_sweep_1d_cylindrical`.

.. _unified-sweep-dispatch:

Unified sweep dispatch
-----------------------

Wave D Round 2 of the SN reshape campaign (Issue #161) consolidates
the four pre-existing sweep paths (1-D cumprod / 2-D wavefront /
spherical / cylindrical) under one
:func:`~orpheus.sn.sweep.transport_sweep` entry point that branches
on a single boolean from the
:class:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator`
primitive (Wave B Issue #6 / Wave D Round 1):

.. code-block:: python

   def transport_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso=None):
       reduced = sn_mesh.reduced
       if reduced is not None and reduced.requires_upstream_angular_state:
           return _curvilinear_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)
       return _cartesian_sweep(Q, sig_t, sn_mesh, psi_bc, Q_aniso)

The pre-Wave-D dispatch did string-equality on
``sn_mesh.curvature == "spherical"`` / ``"cylindrical"`` /
``None``.  The new dispatch reads the canonical geometry-layer
property
:attr:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.requires_upstream_angular_state`
— ``False`` for slab + 2-D Cartesian (no angular redistribution
between successive half-angles), ``True`` for spherical +
cylindrical.  Two-D Cartesian sets ``sn_mesh.reduced is None``
(no curvilinear math is needed), and the dispatch falls through
to the Cartesian path.  Why this matters:

* The :class:`ReducedStreamingOperator` is the primitive that
  already encodes "does this geometry need angular
  redistribution to march the sweep?", so the dispatch reads its
  property directly instead of round-tripping through a
  string tag — Cardinal Rule 2 (architecture).  Consumers
  outside the SN sweep (MoC, CP) read the same property when
  they need the same dispatch.
* The dispatch surface shrinks from four string-equality checks
  to one boolean — a structural simplification that makes the
  control flow easier to reason about and to extend with
  additional cell-update strategies (Wave C-extension).

Cell update strategy parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The curvilinear sweep dispatches per-cell to
:meth:`~orpheus.sn.spatial.cell_update.CellUpdate.update`:

.. code-block:: python

   for cell_idx in march_order:
       st = reduced.streaming_terms(cell_idx, dir_idx, mu_level_idx=p)
       upstream = UpstreamState(
           spatial_upstream=psi_face,
           angular_upstream=psi_angle[cell_idx],
       )
       result = sn_mesh.cell_update.update(
           streaming_terms=st,
           total_xs=sig_t[cell_idx],
           source=QV[cell_idx],
           upstream_state=upstream,
       )
       psi_face = result.outgoing_spatial_flux  # may be None for cylindrical degenerate
       psi_angle[cell_idx] = result.outgoing_angular_state

The cell-update strategy lives on
:attr:`~orpheus.sn.geometry.SNMesh.cell_update` (introduced in
this round as a constructor argument with default
:class:`~orpheus.sn.spatial.diamond.DiamondDifference`).  The
default reproduces the inlined sweep math bit-identically — every
regression snapshot at ``tests/sn/regression/snapshots/`` was
generated with DD and continues to match bit-for-bit when the
unified sweep dispatches via ``cell_update.update(...)``.  See
:ref:`cell-update-strategies` for the strategy contract and
:class:`~orpheus.sn.spatial.diamond.DiamondDifference` for the
DD scalar form.

Wave C-extension will ship :class:`Step`, :class:`LinearDiscontinuous`,
and :class:`ExponentialCharacteristic` as positivity-preserving /
higher-order alternatives; the unified dispatch infrastructure is
in place to receive them — users will pass
``cell_update=LinearDiscontinuous()`` etc. at
:class:`~orpheus.sn.geometry.SNMesh` construction.

The 1-D cumprod fast path (DD-only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Cartesian dispatch checks three preconditions before
selecting the 1-D cumprod fast path:

#. ``cell_update is DiamondDifference`` — the cumprod recurrence
   :eq:`dd-recurrence` is a DD-specific algebraic identity
   (Lewis & Miller §5.3); LD / EC / Step closures do not admit
   the same recurrence.
#. Quadrature is GL1D (``ny == 1`` and all ``mu_y`` vanish).
#. Source is isotropic (``Q_aniso is None``).

If any precondition fails, the Cartesian dispatch routes
through the 2-D wavefront sweep (which handles 1-D as a
special case).  Preserving the cumprod fast path inside the
unified algorithm is required to keep the 1-D regression
snapshots bit-identical and to retain the historical
sub-millisecond sweep time for typical 1-D problems.

The 2-D wavefront sweep retains its inlined DD math (rather
than dispatching through ``cell_update.update``) because anti-
diagonal vectorisation operates on ``(n_diag, ng)`` cell slices,
a shape the per-cell ``CellUpdate`` Protocol does not currently
accept.  Wave C-extension's introduction of LD / EC / Step
strategies will require parameterising 2-D wavefront via the
Protocol while preserving anti-diagonal vectorisation — that's
the open design point for the rollout.  For Wave D, the inlined
DD math is bit-identical to a per-cell
:class:`DiamondDifference` call sequence by construction.

ERR-026 closure status (partial through Wave E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The curvilinear sweep's one-directional WDD closure
:math:`\psi_{n+1/2} = (\overline{\psi} - (1 - \tau_{mm})\,
\psi_{n-1/2})/\tau_{mm}` is preserved bit-identically by
:class:`DiamondDifference` (Wave C extracted it from the
inlined sweep verbatim).  ERR-026 (catalogued in
:doc:`/development` and the V&V matrix at
:doc:`/verification/matrix`) lives in this closure: the
solver's source-iteration path converges to a non-flat
fixed point even though the matrix-free ``apply`` path with
the symmetric closure is exact for **constant** sources.

Wave D's gating contract was bit-identity for ``cell_update =
DiamondDifference`` — the bug is preserved by construction so
the regression snapshots stay green.  Wave E (Issues #98 #99 #164)
took two passes at the closure:

* **Wave E Round 2** wired
  :func:`~orpheus.sn.solver.solve_sn_fixed_source` to route
  through Krylov-on-:meth:`SNStreamingOperator.apply` (the
  symmetric closure) with the sweep-as-``solve`` as preconditioner.
  This closes ERR-026 on **constant-source reflective-BC
  problems** — the canonical
  :file:`tests/sn/test_sweep_operator_inconsistency.py` regression
  suite confirms the krylov path gives the analytical flat flux
  to round-off where the sweep does not.
* **Wave E Round 3** (Issue #98 follow-up) closed the BC-faithfulness
  gap that Round 2 identified: the FD operator's
  :func:`~orpheus.sn.operator.solution_to_angular_flux*` and the
  matvec helpers now consume the
  :class:`~orpheus.geometry.boundary.BoundaryOperator` instances on the
  :class:`~orpheus.sn.geometry.SNMesh` (Wave B Issue 7
  tensor-decomposed BC algebra), dispatching boundary fills via
  :meth:`BoundaryOperator.apply`.  Vacuum, reflective,
  white, periodic, albedo, and mixed BCs are now plumbed
  uniformly through the FD operator; bit-identity to the
  pre-Round 3 hard-coded reflective fill is preserved for
  :class:`SpecularBoundaryOperator` (the standard ``BC.reflective`` factory),
  which is the load-bearing condition for the 11 frozen
  regression snapshots to stay green.

What is **still** open after Round 3: empirically the symmetric-
closure FD operator at the curvilinear outer face uses cell-center
as a face-flux approximation (``psi_right = fi[:, n, i, 0]`` at
``i = nx-1`` for outgoing :math:`\mu > 0`).  This is exact for
constant solutions but only first-order accurate on non-constant
solutions like the manufactured ``A(r) = sin(πr/R)`` ansatz used
by the curvilinear MMS test suite.  Switching the
``solve_sn_fixed_source`` curvilinear default from
``"source_iteration"`` to ``"krylov"`` would *regress* the MMS
convergence rate from the WDD sweep's
~:math:`\mathcal{O}(h^{1.3})` (ERR-026-affected, but a benign
volumetric-error mode for these MMS) to
~:math:`\mathcal{O}(h^{1})` (FD operator's boundary truncation).
Round 3 therefore *keeps* ``inner_solver="source_iteration"`` as
the default for all geometries; ``"krylov"`` is opt-in and
correct for constant-source problems but not the right default
for MMS.

The two ``xfail-strict`` tripwires at
``tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py``
remain ``xfail`` through Round 3 with updated reason strings
reflecting the partial closure.  Full ERR-026 closure on MMS
depends on a follow-up that extrapolates the curvilinear
outer-face flux at second order (DD diamond relation at the
boundary, or analogous ghost-cell technique).

Adams & Larsen 2002 §III.B's "preconditioner correctness vs
operator correctness" frame is the right lens: the sweep's WDD
fixed-point bias is the wrong answer for a *primary solve*, but
as a *preconditioner* the same fixed point is just an effective
scaling of the residual — it does not poison the converged
solution determined by the operator.  The operator must be
correct *and* second-order-accurate; Round 3 closed the
correctness piece (BC-faithfulness), the second-order piece is
the open follow-up.

.. _sn-boundary-face-flux-protocol:

Boundary face-flux strategies (Issue #168 Phase A)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Issue #168 empirical investigation (recorded at
``.claude/agent-memory/numerics-investigator/issue_168_three_defects.md``)
found **three independent O(1) boundary truncation defects** in the
historical curvilinear FD operator.  The single-defect framing of
the original issue was empirically refuted: fixing only the outer-face
extrapolation produces a temporary improvement at fixed mesh count
but the order *degrades* with refinement because the storage
conflation and the sphere-pole stencil dominate as ``h → 0``.

Phase A (this session) addresses **Defects 1 + 2** via a new
:class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`
strategy Protocol — analogous to the Wave C
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` Protocol — and
a structural decoupling of cell-centre storage from BC face-value
storage in
:func:`~orpheus.sn.operator.solution_to_angular_flux_spherical`
(and the cylindrical alias).

* **Defect 1 — outer-face cell-centre truncation.** The pre-Phase-A
  operator used :math:`\psi^{\rm face}_{N-1/2} = \psi_{N-1}`
  (cell-centre as face value), which is :math:`\mathcal{O}(h)`
  accurate on non-constant solutions.  Phase A's
  :class:`~orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`
  strategy uses the one-sided second-order DD diamond extrapolation

  .. math::

     \psi^{\rm face}_{N-1/2} \;=\; \tfrac{3}{2}\,\psi_{N-1}
                                  \;-\; \tfrac{1}{2}\,\psi_{N-2}

  which restores :math:`\mathcal{O}(h^2)` truncation at the outer face.

* **Defect 2 — cell-centre / BC-face-value storage conflation.** The
  pre-Phase-A
  :func:`~orpheus.sn.operator.solution_to_angular_flux_spherical`
  wrote the BC face value into ``fi[..., -1, 0]`` for inward μ < 0
  ordinates — the same slot the matvec at i=N-2 read as a
  *cell-centre* in the symmetric arithmetic-average stencil
  ``0.5·(fi[..., N-2, 0] + fi[..., N-1, 0])``.  Vacuum BC exposed
  the conflation; specular BC accidentally hid it.  Phase A returns
  ``(fi, boundary_face_flux)`` from the decoder — ``fi`` is now
  pure cell-centre storage (the inward slots at i=N-1 are filled
  with the reflected-partner cell-centre, a faithful O(1) value
  independent of BC kind), and the BC face flux lives in its own
  ``boundary_face_flux`` array threaded into the matvec for the
  outer-face read on inward ordinates.

* **Defect 3 — sphere-pole redistribution stencil.** Bailey 2009's
  :math:`\Delta A / w` redistribution at i=0 cancels the streaming
  for *flat* :math:`\psi` but overcorrects for :math:`\psi` varying
  linearly near r=0.  Phase A **preserves the historical Bailey
  treatment** (the
  :class:`~orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`
  strategy returns the cell-centre value at ``cell_idx = 0`` as a
  pass-through).  Defect 3 remains the subject of Phase B with
  literature consultation (Lewis & Miller §4.5; Carlson starting-
  direction; Lathrop pole stencil).  Empirically, Phase A alone
  yields :math:`\mathcal{O}(h^{1.5}\text{--}h^{1.7})` — better than
  the pre-Phase-A :math:`\mathcal{O}(h^{1.25})`, but not yet the
  :math:`\mathcal{O}(h^2)` bar that the four ``xfail-strict``
  tripwires require.  ERR-026 stays at PARTIAL CLOSURE pending
  Phase B.

The Phase A
:class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`
Protocol mirrors the
:class:`~orpheus.sn.spatial.cell_update.CellUpdate` shape — a
``@runtime_checkable Protocol`` plus a parallel concrete ABC
(:class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFluxBase`)
layered on
:class:`~orpheus.numerics.registry.RegistryMixin` so concrete
strategies self-register under a ``key=`` class-creation kwarg.
The ablation/back-bisection counterpart
:class:`~orpheus.sn.spatial.boundary_face_flux.CellCenter` ships
alongside :class:`~orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`
to reproduce the pre-Phase-A Defect-1 behaviour on demand.  See
:mod:`orpheus.sn.spatial.boundary_face_flux` for the full module
docstring and references.

The default boundary face-flux strategy is set on the
:class:`~orpheus.sn.geometry.SNMesh` instance (parallel to the
cell-update strategy):
``SNMesh(mesh, quad, boundary_face_flux=DDExtrapolation())``.  The
Cartesian path is unaffected — the upwind FD stencil there has no
symmetric closure to break — and ignores this attribute.

.. _sn-pole-angular-closure-protocol:

PoleAngularClosure (Issue #168 Phase B)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase B addresses **Defect 3** of Issue #168 — the angular-redistribution
truncation gap on angularly-varying :math:`\psi`.  The pre-Phase-B
operator carried inline τ-symmetric interpolation
:math:`\psi_{n+1/2} \approx \tau_n\,\psi_{n+1} + (1-\tau_n)\,\psi_n`,
which is the **flat-flux collapse** of Hébert (2009) §3.9.4
Eqs. 3.428 + 3.437/3.439 — exact when :math:`\psi` is constant in
:math:`\mu`, but only :math:`\mathcal{O}(1)` accurate on
angularly-varying :math:`\psi`.  Phase B lifts this evaluation into
a :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
strategy Protocol — analogous to Phase A's
:class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux` —
and ships **three concrete strategies** trading off bit-identity,
flat-flux invariance, and asymptotic accuracy:

* :class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
  — bit-for-bit reproduction of the pre-Phase-B inlined math
  (the :math:`\tau`-symmetric form).  **Default**, preserving the
  curvilinear regression-snapshot bit-identity contract and the
  per-ordinate flat-flux invariant the ERR-026 evidence test
  (:file:`tests/sn/test_sweep_operator_inconsistency.py`) relies on.
  Carries Defect 3 by design — the truncation gap is *reproducible*
  so future verification probes can cross-check against the
  documented behaviour.

* :class:`~orpheus.sn.spatial.pole_angular_closure.BaileyFlatFluxRedist`
  — the algebraic flat-flux collapse
  :math:`R_{n,i,g} = (\Delta A/w)\,(\alpha_{n+1/2} - \alpha_{n-1/2})\,
  \psi_{n,i,g} / V_i = -\mu_n\,\Delta A_i\,\psi_{n,i,g} / V_i`
  (using :eq:`bailey-dome-recursion`).  Equivalent to the legacy form
  on flat :math:`\psi` (the
  :func:`tests.sn.l1_analytical.test_pole_closure_flat_flux_identity.test_spherical_flat_flux_legacy_matches_bailey_collapse`
  test pins this), and used as a structurally simpler bridge to the
  flat-flux invariant in unit tests.

* :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
  — the canonical Hébert §3.9.4 per-cell M-M weighted DD angular
  recurrence:

  .. math::
     :label: pole-mm-recurrence

     \phi_{1/2,i,g} &= 0, \\
     \phi_{n+1/2,i,g} &= \frac{\phi_{n,i,g}
                              \;-\; (1 - \tau_n)\,\phi_{n-1/2,i,g}}{\tau_n},
     \qquad n = 1, \ldots, N.

  At :math:`\tau_n = 1/2` (the M-M clamp's lower bound) the recurrence
  reduces to pure DD angular :math:`\phi_{n+1/2,i,g} = 2\,\phi_{n,i,g}
  - \phi_{n-1/2,i,g}` (Hébert Eqs. 3.437 / 3.439).  The
  :math:`\tau \in (1/2, 1]` clamp gives weighted-DD with positive M-M
  weighting per Bailey-Morel-Chang 2010.  The same recurrence runs
  inside :class:`~orpheus.sn.spatial.diamond.DiamondDifference` (the
  sweep's cell update); applying this strategy in the apply matvec
  brings the apply and sweep to the same angular closure, but the
  **spatial** closures still differ (apply uses arithmetic averages
  + DD extrapolation; sweep uses WDD).  Full ERR-026 closure on the
  apply matvec requires aligning the spatial closure also (a
  follow-up beyond Phase B's scope — design memo §6.4 / §11).
  :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
  is therefore **opt-in** in Phase B; the default stays
  :class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`.

α-recursion normalisation
^^^^^^^^^^^^^^^^^^^^^^^^^

Hébert Eq. 3.424 reads :math:`\alpha^{H}_{n+1/2} = \alpha^{H}_{n-1/2}
- 2\,\mathcal{W}_n\,\mu_n` with the corresponding redistribution
divisor :math:`\Delta S_i / (2\,\mathcal{W}_n)` in Eq. 3.428.  The
ORPHEUS arrays carry :math:`\alpha^{O} = \alpha^{H}/2`, absorbing the
factor of 2 into the recurrence; the redistribution divisor reads
:math:`\Delta A_i / w_n` correspondingly.  Both forms are
mathematically equivalent.  This normalisation is documented in
:mod:`orpheus.geometry.reduced_operator` and re-stated explicitly in
:mod:`orpheus.sn.spatial.pole_angular_closure` so the Hébert canonical
form's connection to the ORPHEUS arrays is transparent.

Citation correction
^^^^^^^^^^^^^^^^^^^

Pre-Phase-B the codebase cited "Bailey, T. S., Adams, M. L., Yang,
B., & Zika, M. R. (2009). *A piecewise linear finite element
discretization of the diffusion equation for arbitrary polyhedral
grids*. JCP 227, 3738-3757" for the curvilinear S\ :sub:`N`
:math:`\alpha`-recursion.  This is the **wrong Bailey paper** —
Bailey-Adams-Yang-Zika is a piecewise-linear FE diffusion paper
unrelated to S\ :sub:`N`.  The intended reference is **Bailey,
Morel & Chang (2010)**, NSE 165(2):149-169, "Asymptotic
Diffusion-Limit Accuracy of Sn Angular Differencing Schemes"
(LLNL preprint LLNL-JRNL-420356; OA at
https://www.osti.gov/servlets/purl/1020346).  Phase B corrects the
citations in :mod:`orpheus.geometry.reduced_operator`,
:mod:`orpheus.sn.sweep`, :mod:`orpheus.sn.spatial.diamond`, and the
new :mod:`orpheus.sn.spatial.pole_angular_closure` module.  Hébert
(2009) §3.9.4 is the **primary source** for the curvilinear S\ :sub:`N`
discretization in this codebase; Bailey-Morel-Chang 2010 is the
auxiliary justification for the M-M weighted-diamond :math:`\tau`
clamp.

ERR-026 closure status (Phase B partial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Phase B ships the architectural infrastructure for closing Defect 3
(the :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
Protocol + canonical
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`)
without flipping the default — the canonical form in isolation does
not produce an :math:`\mathcal{O}(h^2)` MMS rate because the apply
matvec's spatial closure still differs from the sweep's WDD form.
The four ``xfail-strict`` curvilinear MMS tripwires therefore stay
xfail through Phase B; ERR-026 stays at **PARTIAL CLOSURE** (Phase A
closed Defects 1+2 spatial, Phase B ships Defect 3 architectural
scaffolding).  The full closure requires a Phase C follow-up that
aligns the apply matvec's spatial closure with the sweep's WDD
form, at which point
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
becomes the natural default.

Starting Direction
-------------------

Both curvilinear sweeps initialise the angular face flux
:math:`\psi_{1/2}` to **zero** for each spatial cell.  This is valid
because :math:`\alpha_{1/2} = 0` by construction, so the product
:math:`\alpha_{1/2} \psi_{1/2}` in the balance equation vanishes
regardless of the value of :math:`\psi_{1/2}`.  The starting-direction
treatment of [LewisMiller1984]_ Section 4.5.4 (which tracks
:math:`\alpha\psi` instead of :math:`\psi` alone) is therefore
unnecessary when the :math:`\alpha` recursion is implemented correctly.


Krylov inner solver
===================

Instead of sweep-based source iteration, the within-group transport
problem can be solved directly using a Krylov method.  Wave E Round 2
(Issue #164) replaces the legacy BiCGSTAB FD-operator path with GMRES
on :meth:`SNStreamingOperator.apply` (the symmetric closure carried
by the operator algebra) with the sweep as a left preconditioner —
the SAILOR / Larsen-Adams preconditioned-Krylov framework
([AdamsLarsen2002]_ §III).

Operator equation
-----------------

The streaming-collision operator :math:`L` is formed explicitly via
:class:`~orpheus.sn.operator.SNStreamingOperator`:

.. math::

   L\psi = \mu_n \nabla\psi + \Sigt{}\psi

For Cartesian geometry, this is a banded matrix with upwind cell-
centre finite-difference gradients.  For curvilinear geometries, the
operator includes the :math:`\Delta A/w` geometry factor and the
**symmetric** Morel--Montry angular face-flux interpolation
:math:`\psi_{n\pm 1/2} = \tau\,\psi_{\rm next} + (1-\tau)\,
\psi_{\rm this}`, which is distinct from the WDD asymmetric closure
:math:`\psi_{n+1/2} = (\overline{\psi} - (1-\tau)\,\psi_{n-1/2})/\tau`
used by the sweep.

The system :math:`L\,\psi = b` is solved with
``scipy.sparse.linalg.gmres`` (replacing the legacy BiCGSTAB), with
the sweep wrapped as a scipy ``LinearOperator`` preconditioner ``M``.
On uniform meshes the symmetric and WDD closures converge to the
same physics in the fine-mesh limit; on curvilinear meshes the
sweep's WDD closure has the ERR-026 closure-bias-driven self-
consistent fixed point that is now bypassed by the Krylov-on-
:meth:`apply` path.

Consistency with the Sweep
---------------------------

Both the sweep and the Krylov path must use the **same** spatial
discretisation to produce identical eigenvalues.  In practice:

- The sweep uses diamond-difference (DD): :math:`L^{-1}` is applied
  implicitly via forward substitution along sweep-direction visits.
- The Krylov path forms :math:`L` explicitly via the symmetric-
  closure FD operator and inverts it via GMRES.

On coarse meshes, DD and FD have different truncation error
constants, so the two paths may give slightly different :math:`\keff`
values.  They converge to the same answer as :math:`h \to 0`.

The curvilinear FD operator reads ``alpha_half`` / ``redist_dAw`` /
``tau_mm`` (spherical) and the per-level analogues (cylindrical)
from ``SNMesh.reduced``.  This ensures both paths share exactly the
same connection-coefficient infrastructure.

Round 2 deviation
-----------------

The campaign plan called for an automatic ``"krylov"`` dispatch on
curvilinear meshes in
:func:`~orpheus.sn.solver.solve_sn_fixed_source` that would silently
close ERR-026 on the curvilinear vacuum-BC MMS cases.
Implementation surfaced an unforeseen coupling: the
:func:`~orpheus.sn.operator.build_equation_map_spherical` /
``build_equation_map_cylindrical`` packed-vector layout that
:meth:`SNStreamingOperator.apply` reuses was designed for the
**reflective** outer-boundary BC only — it has no slot for a vacuum-
BC outer-incoming :math:`\psi`.  Wave E Round 3 owns the equation-map
extension that closes the vacuum-BC path; Round 2 ships the krylov
inner solver as an explicit opt-in for the reflective-BC eigenvalue
case where it is bit-identical math to the legacy BiCGSTAB FD path.


.. _boundary-conditions:

Boundary Conditions
===================

Infrastructure
--------------

Boundary conditions are declared on the **geometry mesh** and resolved
by the **solver's augmented mesh** at construction time.  This two-stage
design separates physics intent (what condition to apply) from solver
mechanics (how to enforce it in the sweep).

**Stage 1 --- Geometry declaration.**
:class:`~geometry.mesh.Mesh1D` carries ``bc_left: BC | None`` and
``bc_right: BC | None``; :class:`~geometry.mesh.Mesh2D` carries
``bc_xmin``/``bc_xmax``/``bc_ymin``/``bc_ymax: BC | None``.
:class:`~geometry.mesh.BC` is a frozen dataclass with two fields:

- ``kind: str`` --- an identifier such as ``"vacuum"``, ``"reflective"``,
  or ``"white"``.
- ``params: dict[str, float]`` --- optional numeric parameters
  (e.g. ``{"albedo": 0.7}``).

Convenience instances are available for the common cases:
:attr:`BC.vacuum <geometry.mesh.BC.vacuum>`,
:attr:`BC.reflective <geometry.mesh.BC.reflective>`, and
:attr:`BC.white <geometry.mesh.BC.white>`.
When a face is left as ``None``, the solver applies its own default
(reflective for the SN solver, matching the infinite-lattice /
eigenvalue convention).

**Stage 2 --- Solver resolution.**
:class:`SNMesh` owns a class-level :attr:`~SNMesh.BOUNDARY_OPERATOR_REGISTRY` dictionary
mapping kind strings to factory callables::

    BOUNDARY_OPERATOR_REGISTRY = {
        "vacuum":     _sn_vacuum_boundary_operator,
        "reflective": _sn_reflective_boundary_operator,
    }

During ``SNMesh.__init__``, each face's :class:`~geometry.mesh.BC` is
looked up in the registry.  If the kind is not found, a ``ValueError``
lists the supported kinds.  For curvilinear geometries (spherical,
cylindrical), only ``"reflective"`` and ``"vacuum"`` are currently
supported --- requesting any other kind on a curvilinear mesh raises
``ValueError``.

As of Wave B Round 3 of the SN reshape campaign (Issue 7 of
``.claude/plans/sn_reshape.md``), the registry factories return
concrete :class:`~orpheus.geometry.boundary.BoundaryOperator` instances
rather than string tags. ``sn_mesh.bc_left``, ``sn_mesh.bc_right``,
``sn_mesh.bc_xmin``, etc. carry tensor-decomposed BC objects whose
``apply(angular_flux_outgoing, quadrature)`` method the
sweep calls directly, with no string-kind branching at the call site.
For backward compatibility with existing tests, the concrete BC
classes still compare equal to their legacy string tag
(``sn_mesh.bc_left == "vacuum"`` evaluates True iff
``sn_mesh.bc_left`` is a ``VacuumBoundaryOperator``); see
:ref:`bc-tensor-decompositions` for the full algebra and the list of
implemented primitives.

**Backward compatibility.**
:func:`solve_sn_fixed_source` still accepts a ``boundary_condition: str``
parameter (default ``"vacuum"``).  Internally it calls
``_apply_default_bcs(mesh, boundary_condition)``, which applies the
string to **all faces** that lack explicit :class:`~geometry.mesh.BC`
declarations.  When the mesh already carries explicit BCs, the parameter
is silently ignored --- mesh-level declarations always take precedence.
:func:`solve_sn` (the eigenvalue entry point) does not expose a
``boundary_condition`` parameter; eigenvalue problems use whatever the
mesh declares (defaulting to reflective on all faces).

.. note::

   Before this infrastructure existed, the SN solver hardcoded
   reflective BCs on all faces and :func:`transport_sweep` accepted
   a ``boundary_condition: str`` parameter.  That parameter has been
   removed --- BCs now flow exclusively through the mesh → SNMesh
   resolution path described above.

Supported Types
---------------

**Reflective** (specular reflection).
At the outer boundary :math:`r = R` (or :math:`x = L`), the incoming
flux for ordinate :math:`n` is set to the outgoing flux of its reflected
partner:

.. math::
   :label: reflective-bc

   \psi_n^{\rm in} = \psi_{n'}^{\rm out}

where :math:`n'` is the reflected partner ordinate (negating the
appropriate direction cosine).  Reflective partner indices are precomputed
by each quadrature's :meth:`reflection_index` method.  This is the
default for eigenvalue problems (infinite lattice / infinite medium).
The CP solver uses white (isotropic) BCs instead; see
:ref:`white-bc-quality` for a comparison showing the ~1% gap between
the two approaches.

**Vacuum** (zero incoming flux).
All incoming angular fluxes at the face are set to zero:

.. math::
   :label: vacuum-bc

   \psi_n^{\rm in} = 0

.. (vv-status rationale) definition: Definitional / notation introduction. Operational rule ψ_n^in = 0 for vacuum boundary; semantics exercised by every vacuum-BC test (test_boundary_conditions, MMS suite); no isolated identity to verify.
.. vv-status: vacuum-bc documented


In the 1-D cumprod path, this means the recurrence starts from zero
instead of the reflected outgoing flux.  In the 2-D wavefront sweep,
the reflective-partner copy is skipped, leaving incoming-face angular
fluxes at their zero initialisation.  Vacuum BCs are the natural
choice for fixed-source MMS verification on finite slabs (see
:ref:`sn-mms-verification`).

.. _bc-tensor-decompositions:

Boundary conditions as tensor decompositions
---------------------------------------------

The boundary conditions used by :class:`SNMesh` are concrete instances
of a more general tensor-decomposed framing, defined in
:mod:`orpheus.geometry.boundary`. A boundary condition is a linear
operator :math:`R` mapping the outgoing angular flux at a face to the
incoming angular flux:

.. math::
   :label: bc-tensor-decomposition

   \psi_{\rm in}(\Omega)
   = (R\,\psi_{\rm out})(\Omega)
   = \sum_\alpha \bigl(G_\alpha\,\psi_{\rm out}\bigr)(\Omega) \cdot A_\alpha,

where :math:`G_\alpha` is a **geometric operator** (permutation,
pushforward, angular average, spatial wrap) and :math:`A_\alpha` is a
**scalar amplitude** (typically an albedo :math:`\in [0, 1]`).

This is the same algebra Lewis & Miller (1984) §3.4 use to introduce
boundary conditions in transport: every BC of practical interest is
either rank-1 (one :math:`G \otimes A` term) or a finite linear
combination of rank-1 primitives (rank-N). The implemented primitives
are:

.. list-table:: Implemented :class:`BoundaryOperator` primitives
   :widths: 20 30 25 25
   :header-rows: 1

   * - Class
     - :math:`G_\alpha`
     - :math:`A_\alpha`
     - Rank / wired into ``solve_sn``
   * - :class:`~orpheus.geometry.boundary.VacuumBoundaryOperator`
     - :math:`0` (no operator)
     - 0
     - 0 / yes
   * - :class:`~orpheus.geometry.boundary.SpecularBoundaryOperator`
     - permutation under reflection axis
     - albedo (1 = perfect)
     - 1 / yes
   * - :class:`~orpheus.geometry.boundary.WhiteBoundaryOperator`
     - cosine-weighted hemispheric average
     - albedo
     - 1 / no (Wave C)
   * - :class:`~orpheus.geometry.boundary.PeriodicBoundaryOperator`
     - spatial pushforward (caller-supplied)
     - 1
     - 1 / no (Wave C/D)
   * - :class:`~orpheus.geometry.boundary.AlbedoBoundaryOperator`
     - identity in angle
     - albedo
     - 1 / no (building block)
   * - :class:`~orpheus.geometry.boundary.MixedBoundaryOperator`
     - sum of components
     - per-component
     - N / no

The Protocol :class:`~orpheus.geometry.boundary.BoundaryOperator` exposes one
method, ``apply(angular_flux_outgoing, quadrature)``, which
the SN sweep calls instead of branching on a string kind. The
:attr:`SNMesh.BOUNDARY_OPERATOR_REGISTRY` factories now return concrete
:class:`BoundaryOperator` instances; the registry pattern is unchanged from a
caller's perspective. White, periodic, albedo, and mixed primitives
ship as building blocks --- the sweep plumbing for them lands in a
later wave of the SN reshape campaign (this issue is the algebra; the
plumbing is the consumption).

The tensor framing pays off architecturally because partial-current
Marshak boundaries (:math:`R = c_1 \, G_{\rm refl} + c_2 \, G_{\rm
diff}`, Bell & Glasstone 1970 §1.5) and multi-region interface
couplings are all instances of the same algebra: pick the geometric
operators, pick the amplitudes, sum. Once this shape is in place, new
BCs are one ``BoundaryOperator`` class and one
``BOUNDARY_OPERATOR_REGISTRY`` entry away --- no sweep edits per BC.

Inner Boundary (Curvilinear)
-----------------------------

At :math:`r = 0`:

- The face area :math:`A(0) = 0`, so **no spatial flux crosses the
  origin**.  The spatial incoming flux for outward-sweeping ordinates is
  zero.
- The **angular redistribution provides the inward-to-outward
  transition**: flux entering as an inward-directed ordinate
  (:math:`\mu < 0` or :math:`\eta < 0`) is redistributed to outward
  ordinates (:math:`\mu > 0` or :math:`\eta > 0`) through the
  :math:`\alpha` coupling.

This means the curvilinear sweep does not need an explicit boundary
condition at :math:`r = 0` --- the geometry handles it naturally.
Curvilinear sweeps currently only support reflective BCs on the outer
face; this is enforced by the validation in
:meth:`SNMesh._resolve_one`.


Scattering
==========

Matrix Convention
-----------------

The ``Mixture.SigS[l]`` matrices use the convention
:math:`\text{SigS}[g_{\rm from}, g_{\rm to}]`:

.. math::

   \text{SigS}[0] = \begin{pmatrix}
       \Sigs{0\to0} & \Sigs{0\to1} \\
       \Sigs{1\to0} & \Sigs{1\to1}
   \end{pmatrix}

For the in-scatter source (total scattering into group :math:`g` from
all groups :math:`g'`):

.. math::

   Q_{\rm scatter}[g]
   = \sum_{g'} \Sigs{g'\to g}\, \phi_{g'}
   = (\text{SigS}^T \cdot \phi)[g]

The vectorised form for batched cells is ``phi @ SigS`` (equivalent to
:math:`(\text{SigS}^T \phi^T)^T` for row-vector :math:`\phi`).

The **analytical eigenvalue problem** uses:

.. math::

   A &= \text{diag}(\Sigt{}) - \text{SigS}^T \quad\text{(removal matrix)} \\
   F &= \chi \otimes \nSigf{} \quad\text{(fission matrix)} \\
   \kinf &= \lambda_{\max}(A^{-1} F)

Note the transpose: :math:`\text{SigS}^T[g, g'] = \Sigs{g'\to g}` gives
the in-scatter contribution, so :math:`\text{diag}(\Sigt{}) - \text{SigS}^T`
is the net removal matrix.  See :ref:`scattering-matrix-convention` for the
full derivation of this convention.

P\ :sub:`0` Isotropic Scattering
----------------------------------

The default mode (``scattering_order=0``).  A direction-independent
source is added to all ordinates equally:

.. math::

   Q_{\rm scatter}(\hat{\Omega}_n)
   = \sum_{g'} \Sigs{g'\to g}^{(0)}\, \phi_{g'} / W

Implemented in :meth:`SNSolver._add_scattering_source`, which performs
``phi @ SigS[0]`` per material.

.. _pn-scattering:

P\ :sub:`N` Anisotropic Scattering
------------------------------------

When ``scattering_order >= 1``, per-ordinate anisotropic sources are
computed from the Legendre moments of the angular flux.  The full
anisotropic scattering source for ordinate :math:`n` and group :math:`g`
is:

.. math::
   :label: pn-scatter

   Q_{\rm scatter}(\hat{\Omega}_n, g)
   = \sum_{\ell=0}^{L} (2\ell+1)
     \sum_{m=-\ell}^{\ell}
     \sum_{g'} \Sigs{g'\to g}^{(\ell)}\,
     f_{\ell,g'}^m \; Y_\ell^m(\hat{\Omega}_n)

where :math:`Y_\ell^m` are real spherical harmonics and the angular flux
moments are computed by quadrature:

.. math::
   :label: flux-moments

   f_{\ell,g}^m = \sum_{n=1}^{N} w_n \, \psi_{n,g} \, Y_\ell^m(\hat{\Omega}_n)

The :math:`(2\ell+1)` factor is the addition theorem normalisation for
real spherical harmonics: it ensures that the P\ :sub:`L` expansion
reproduces the angular flux moments exactly when the angular flux is a
polynomial of degree :math:`\leq L`.

**Implementation in** :meth:`SNSolver._build_aniso_scattering`:

1. **Compute spherical harmonics** at construction time:
   :math:`Y[n, \ell, \ell+m]` for all ordinates, stored as ``self._Y``
   with shape ``(N, L+1, 2L+1)``.

   **Convention.** The polar axis is :math:`\mu_x`, so
   :math:`\cos\theta = \mu_x` and
   :math:`\sin\theta = \sqrt{1 - \mu_x^2}`.  Azimuth is measured in the
   :math:`(\mu_y, \mu_z)` plane:
   :math:`\cos\phi = \mu_y / \sin\theta`,
   :math:`\sin\phi = \mu_z / \sin\theta`.  This matches the MATLAB
   ``discreteOrdinatesPWR.m`` reference for :math:`\ell \le 1`:

   .. math::

      Y_0^0 &= 1 = P_0(\cos\theta)\\
      Y_1^{-1} &= \sin\theta\,\sin\phi = \mu_z\\
      Y_1^0    &= \cos\theta              = \mu_x\\
      Y_1^{1}  &= \sin\theta\,\cos\phi   = \mu_y

   For :math:`\ell \ge 2` the formula extends as standard real
   spherical harmonics in this frame:

   .. math::
      :label: real-spherical-harmonics

      Y_\ell^0(\hat{\Omega}) &= P_\ell(\mu_x)\\
      Y_\ell^{m}(\hat{\Omega}) &= \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}}\,
                                 P_\ell^{m}(\mu_x)\,\cos(m\phi),
                                 \quad m > 0\\
      Y_\ell^{-m}(\hat{\Omega}) &= \sqrt{\frac{2(\ell-m)!}{(\ell+m)!}}\,
                                  P_\ell^{m}(\mu_x)\,\sin(m\phi),
                                  \quad m > 0

   where :math:`P_\ell^m` is the unnormalised associated Legendre
   function (the :math:`(-1)^m` Condon–Shortley phase included by
   ``scipy.special.lpmv`` is removed at the call site).  The
   normalisation is the **"no** :math:`4\pi/(2\ell+1)` **prefactor"**
   convention under which the addition theorem reads

   .. math::
      :label: addition-theorem

      \sum_{m=-\ell}^{\ell} Y_\ell^m(\hat{\Omega})\,Y_\ell^m(\hat{\Omega}')
      = P_\ell(\hat{\Omega} \cdot \hat{\Omega}')

   which is the identity used by Eq. :eq:`pn-scatter` to expand the
   :math:`P_\ell` scattering kernel as a finite tensor product over
   :math:`m`.  Equivalently the discrete orthogonality on a quadrature
   exact for polynomials of degree :math:`\ge 2\ell` reads

   .. math::

      \sum_{n=1}^{N} w_n \, Y_\ell^m(\hat{\Omega}_n)\,
                            Y_{\ell'}^{m'}(\hat{\Omega}_n)
      = \frac{4\pi}{2\ell+1}\,\delta_{\ell\ell'}\,\delta_{mm'}.

   Both identities are verified at :math:`\ell \le 3` by
   ``test_spherical_harmonics_addition_theorem_L3`` and
   ``test_spherical_harmonics_orthogonality_L3`` in
   ``tests/sn/test_solver_components.py``.  The :math:`\ell \le 1`
   branch is kept as bit-identical hardcoded values so existing
   :math:`P_0/P_1` test outputs are preserved at any tolerance
   (``test_spherical_harmonics_l1_unchanged_after_extension``).

2. **Compute flux moments** via an ``einsum`` contraction over the
   ordinate index:

   .. code-block:: python

      fiL[:, :, :, l, l+m] = np.einsum(
          'n,nxyg->xyg', w * Y[:, l, l+m], angular_flux,
      )

   This contracts :math:`\sum_n w_n Y_\ell^m(\hat{\Omega}_n) \psi_n(x,y,g)`
   into a spatial-energy field of shape ``(nx, ny, ng)``.

3. **Reconstruct per-ordinate source**: for each Legendre order
   :math:`\ell \geq 1` (the :math:`\ell = 0` term is handled by
   :meth:`SNSolver._add_scattering_source`) and each :math:`m`, the
   scattered moment ``moment @ sig_s_l[l]`` is multiplied by
   :math:`(2\ell+1) Y_\ell^m(\hat{\Omega}_n)` and accumulated into
   ``Q_aniso[n, :, :, :]``.

4. The resulting ``Q_aniso`` array of shape ``(N, nx, ny, ng)`` is
   passed to :func:`transport_sweep`, which adds it to the isotropic
   source on a per-ordinate basis.

**Equivalence of the code to the mathematical form.**
Equation :eq:`pn-scatter` writes the sum as
:math:`\sum_\ell \sum_m \sum_{g'} \Sigs{}^{(\ell)} f_\ell^m Y_\ell^m`.
The code separates the :math:`\ell = 0` term (isotropic, handled by
``_add_scattering_source``) from the :math:`\ell \geq 1` terms
(anisotropic, handled by ``_build_aniso_scattering``).  For :math:`\ell = 0`,
:math:`Y_0^0 = 1` and :math:`(2 \cdot 0 + 1) = 1`, so the sum reduces to
:math:`\sum_{g'} \Sigs{g' \to g}^{(0)} f_{0,g'}^0 = \sum_{g'} \Sigs{g' \to g}^{(0)} \phi_{g'}`,
which is exactly the P\ :sub:`0` source.  The split is therefore exact
with no double-counting.

The 421-group cross-section library provides both P0 and P1 matrices.

.. _n2n-reactions:

(n,2n) Reactions
-----------------

The :math:`(n,2n)` reaction is a threshold reaction in which a neutron
is absorbed by a nucleus, which then emits **two** neutrons.  The net
effect is a gain of one neutron per reaction (the incident neutron is
consumed, two are produced).

The :math:`(n,2n)` cross section is stored as a group-to-group transfer
matrix ``Mixture.Sig2`` with the same ``[g_from, g_to]`` convention as
the scattering matrix.  The source contribution is:

.. math::

   Q_{(n,2n)}(g) = 2 \sum_{g'} \Sigma_{2,g'\to g}\, \phi_{g'}

The factor of 2 accounts for the two neutrons produced per reaction.
The implementation in :meth:`SNSolver._add_n2n_source` performs:

.. code-block:: python

   Q[ix, iy, :] += 2.0 * (phi[ix, iy, :] @ self.sig2[mid])

This is added to the isotropic source before the transport sweep, on the
same footing as the P\ :sub:`0` scattering source.  The :math:`(n,2n)`
contribution also enters the :math:`\keff` production term in
:meth:`SNSolver.compute_keff`, where row sums of ``Sig2`` (total
:math:`(n,2n)` removal rate) are used.

Normalization Chain
--------------------

The normalization chain in the code ensures consistent scaling:

1. **Fission source** (:meth:`SNSolver.compute_fission_source`):
   :math:`Q_f = \chi \cdot (\nSigf{} \cdot \phi) / k` --- raw,
   un-normalised.

2. **Scattering source** (:meth:`SNSolver._add_scattering_source`):
   :math:`Q_s = \text{SigS}^T \cdot \phi` --- also un-normalised.

3. **Sweep** (:func:`transport_sweep`): applies
   :math:`Q_{\rm scaled} = Q \cdot w_{\rm norm}` where
   :math:`w_{\rm norm} = 1/\sum w_n`.  This is the :math:`1/W` division
   in the S\ :sub:`N` equation.

4. **Scalar flux** (inside sweep):
   :math:`\phi = \sum_n w_n \psi_n` --- standard quadrature integration.

5. **keff** (:meth:`SNSolver.compute_keff`):
   :math:`k = (\nSigf{} \cdot \phi \cdot V) / (\Sigma_a \cdot \phi \cdot V)`
   --- volume-weighted ratio.

The :math:`1/W` in step 3 and the :math:`W` implicit in step 4 cancel:
:math:`\phi = \sum w_n \cdot Q/(W \Sigt{}) = Q/\Sigt{}` for uniform
isotropic source.

**Convention rule:** Sources passed to the sweep must NOT include
:math:`1/W` --- the sweep applies it.  The BiCGSTAB path (direct
operator) must divide sources by :math:`W` itself, since it solves
:math:`T\psi = b` without the sweep.


.. _sn-scattering-fission-operators:

Scattering and fission as LinearOperators
==========================================

Wave A Issue 1 of the SN reshape campaign installed the
:class:`~orpheus.numerics.operator.LinearOperator` Protocol --- a
capability-tagged matrix-free operator algebra (see :ref:`operator-algebra`).
The neutron transport eigenvalue
problem written in operator form is

.. math::

    (L - S - F)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (L - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`L` is the streaming + collision operator, :math:`S` is
the scattering source operator, and :math:`F` is the fission source
operator. Wave D Issue 13 lifts :math:`S` and :math:`F` out of
:class:`~orpheus.sn.solver.SNSolver` and into
:class:`~orpheus.sn.scattering.ScatteringOperator` and
:class:`~orpheus.sn.fission.FissionOperator` respectively. The math is
**moved verbatim** --- the regression contract on the 11 frozen
snapshots at ``tests/sn/regression/snapshots/`` gates the extraction.

Why ``apply``-only
------------------

Both operators advertise a single capability: ``{"apply"}`` only.

* **Scattering (:math:`S`)** is rank-:math:`O(N_{\text{cells}}\cdot
  N_{\text{groups}})`. There is no efficient inverse --- the operator
  is *applied*, never *inverted*. The upper-triangular structure that
  would make a sweep-based ``solve`` tractable does not survive the
  Pℓ Galerkin reconstruction. Any algebraic consumer that asks for
  :math:`S^{-1}` will get a :class:`MissingCapability` at composition
  time, never silently wrong results at call time --- this is the
  load-bearing payoff of the capability-set design (see
  :ref:`operator-algebra`).

* **Fission (:math:`F`)** has rank-1-in-energy structure: the action
  factorises as :math:`(F\phi)_g = \chi_g\,\sum_{g'}\nu\Sigma_{f,g'}
  \phi_{g'}`, an outer product of the emission spectrum with a scalar
  per-cell rate. This rank-1 structure forbids a useful inverse on
  the energy axis (the rate has lost direction information). The
  :math:`1/k` eigenvalue division stays at the **solver** level ---
  :meth:`~orpheus.sn.fission.FissionOperator.apply` returns
  :math:`F\,\phi`; the EigenvalueSolver Protocol's
  ``compute_fission_source`` divides by :math:`k`. This separation
  preserves linearity of the operator (Wave A Protocol contract: an
  operator's ``apply`` is independent of solver state).

Pℓ Galerkin projection on :math:`Y_\ell^m`
-------------------------------------------

The :math:`\ell\ge 1` contribution to :math:`S` is the
:eq:`pn-scatter` Galerkin reconstruction in real spherical
harmonics :math:`Y_\ell^m`, expanded with the discrete-orthogonality
identity :eq:`addition-theorem` (the Lebedev-quadrature L0
verification of the addition theorem lives at
``tests/sn/test_solver_components.py::TestAnisotropicScattering
::test_spherical_harmonics_addition_theorem_L3``).

Per-cell flux moments :eq:`flux-moments` are computed by the
discrete projection

.. math::

    \phi^{\ell m}_g(\vec r)
    \;\approx\; \sum_n w_n\,\psi_{n,g}(\vec r)\,Y_\ell^m(\hat\Omega_n).

The reconstruction back to per-ordinate scattering source uses the
addition theorem:

.. math::

    Q^{(\ell\ge 1)}_{n,g}(\vec r)
    = \sum_{\ell=1}^{L}\,(2\ell+1)\,\sum_m Y_\ell^m(\hat\Omega_n)\,
      \sum_{g'}\Sigma_{s,\ell}^m(g'\to g)\,\phi^{\ell m}_{g'}(\vec r).

The :math:`(2\ell+1)` factor is the discrete-orthogonality
normalisation
:math:`\langle Y_\ell^m | Y_{\ell'}^{m'}\rangle =
(4\pi/(2\ell+1))\,\delta_{\ell\ell'}\delta_{mm'}` working out across
both projection and reconstruction. The Galerkin frame is **real**
spherical harmonics (the
:meth:`~orpheus.sn.quadrature.AngularQuadrature.spherical_harmonics`
implementation), not complex --- this is the convention native to
the Lebedev tabulation and avoids carrying complex arithmetic
through the source-iteration inner loop.

(n,2n) doubling
---------------

The (n,2n) reaction emits **two** neutrons per absorption, producing
a secondary source :math:`2\,\Sigma_{2n}^m(g'\to g)\,\phi_{g'}`. ORPHEUS
folds this into the **scattering** side of the algebra (rather than
giving it its own operator) because:

1. The bookkeeping is identical to in-scatter (vectorise-by-material,
   add-into-:math:`Q`).
2. The legacy code placement
   (:meth:`SNSolver._add_n2n_source`) is inside the same source-
   construction block as scattering. Wave D Issue 13's bit-identical
   extraction needed to preserve that placement to keep the regression
   snapshots bit-identical.
3. Architecturally, both are *secondary-emission scalar-flux-driven*
   sources --- they belong to the same algebra slot.

Backward references — Wave E status
------------------------------------

Wave E Round 2 (Issue #164) wired the operator algebra
:math:`(L, S, F)` into :class:`SNSolver` and replaced the legacy
BiCGSTAB inner-solver path with Krylov-on-:meth:`L.apply` (GMRES
with the sweep as preconditioner).  The
``build_transport_linear_operator*`` and ``build_rhs*`` helpers
were retired; the per-method delegators on
:class:`SNSolver` (:meth:`_add_scattering_source`,
:meth:`_build_aniso_scattering`, :meth:`_add_n2n_source`,
:meth:`compute_fission_source`) remain as thin wrappers over the
new operators for the EigenvalueSolver Protocol surface.

Wave E Round 3 (Issue #98 follow-up) extended the FD operator's
boundary handling to consume the
:class:`~orpheus.geometry.boundary.BoundaryOperator` infrastructure
(Wave B Issue 7), so :func:`solution_to_angular_flux*` and the
matvec helpers now dispatch boundary fills via
:meth:`BoundaryOperator.apply` — vacuum, reflective, white,
albedo, periodic, and mixed BCs are honoured uniformly.


.. _sn-streaming-operator:

SNStreamingOperator: the streaming-collision operator algebra
==============================================================

Wave D Round 3 (Issue #160) installs
:class:`~orpheus.sn.operator.SNStreamingOperator` as the unified
:class:`~orpheus.numerics.operator.LinearOperator` for the
streaming-collision operator
:math:`L = \Omega\cdot\nabla + \Sigma_t`.  This is the Wave D
**capstone**: with :math:`L`, :math:`S`, and :math:`F` all carrying
the Wave A operator-algebra contract, the operator-form transport
equation

.. math::

    (L - S - F)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (L - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

is now expressible in ORPHEUS as a single Python expression composed
from three :class:`LinearOperator` objects.  The Wave A composers
(:class:`~orpheus.numerics.operator.OperatorSum`,
:class:`~orpheus.numerics.operator.OperatorProduct`,
:class:`~orpheus.numerics.operator.ScaledOperator`) compute the
composition's capability set from the constituents per the closure
laws (see :ref:`operator-algebra`).

Three capabilities
------------------

:class:`SNStreamingOperator` advertises ``{"apply", "solve",
"apply_transpose"}`` — every member of the Wave A capability set:

* :meth:`~orpheus.sn.operator.SNStreamingOperator.apply` —
  matrix-free forward action :math:`L\,\psi`.  Reuses the
  symmetric closure math from the existing
  :func:`~orpheus.sn.operator.transport_operator_matvec`,
  :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`,
  :func:`~orpheus.sn.operator.transport_operator_matvec_cylindrical`
  functions (the historical BiCGSTAB FD operator).  The math is
  **extracted verbatim**; the new class is a thin Protocol wrapper.

* :meth:`~orpheus.sn.operator.SNStreamingOperator.solve` —
  inverse action :math:`L^{-1}\,q` via the Wave D Round 2 unified
  sweep (:func:`~orpheus.sn.sweep.transport_sweep`).  Bit-identical
  to a direct :func:`transport_sweep` call on the same arguments.

* :meth:`~orpheus.sn.operator.SNStreamingOperator.apply_transpose` —
  adjoint action :math:`L^*\,\psi`.  Implemented via the explicit
  transpose of the dense matrix assembled by probing
  :meth:`apply` with each unit basis vector.  The construction is
  exact by linear algebra and gates the reciprocity invariant
  :math:`\langle L\,\psi,\,\varphi\rangle = \langle\psi,\,
  L^*\,\varphi\rangle` (see "Reciprocity" subsection below).

Apply and solve use **different** closures by design
----------------------------------------------------

This is the load-bearing architectural fact about
:class:`SNStreamingOperator`, and the reason the operator's
:meth:`apply` is **not** bit-identical to its :meth:`solve`.

The historical SN dispatch in ORPHEUS ships **two distinct
discretisations** of the same continuous operator
:math:`L = \Omega\cdot\nabla + \Sigma_t`, built at different times
for different consumers:

* The **finite-difference operator**
  (:func:`transport_operator_matvec_*` in
  :mod:`orpheus.sn.operator`) was built for the BiCGSTAB inner
  solver path (:meth:`SNSolver._solve_bicgstab_*`).  It uses
  upwind cell-center FD on Cartesian and arithmetic face averages
  with **τ-symmetric Morel-Montry angular interpolation** on
  curvilinear (see the "Explicit Transport Operator" subsection of
  the BiCGSTAB Alternative above and the warning at the head of
  :mod:`orpheus.sn.operator`).

* The **sweep operator**
  (:func:`~orpheus.sn.sweep.transport_sweep`, dispatching through
  the Wave D Round 2 :class:`~orpheus.sn.spatial.cell_update.CellUpdate`
  Protocol with :class:`~orpheus.sn.spatial.diamond.DiamondDifference`
  as the default strategy) uses the **WDD asymmetric closure**
  :math:`\psi_{n+1/2} = (\overline\psi - (1-\tau)\,\psi_{n-1/2})/\tau`.
  This is the historical SN sweep's closure: the upper-triangular
  forward-substitution form that lets a sweep run in
  :math:`O(N\cdot N_{\rm cells})` work.

In the fine-mesh limit both discretisations converge to the same
continuous operator.  On coarse meshes they differ:

* For Cartesian the difference is :math:`O(h)` (upwind FD has
  the same first-order consistency as DD on uniform meshes;
  divergence appears on non-uniform meshes — see the warning at
  the head of :mod:`orpheus.sn.operator`).
* For curvilinear the WDD asymmetric closure has a closure-bias-
  driven self-consistent fixed point that is **not** the
  fine-mesh-limit transport solution (ERR-026).

The reconciliation is **Wave E Issue 15**: the SN solver's inner
loop migrates from "sweep with WDD closure" to "Krylov on apply,
sweep as preconditioner".  Krylov on :meth:`apply` uses the
symmetric closure (the one that agrees with analytical references)
as the system to solve; the WDD sweep is invoked only as a
preconditioner that accelerates the Krylov iteration without
poisoning the converged solution with its closure bias.  ERR-026
closes when Wave E lands; the 2 xfail-strict tripwires at
:file:`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
are the gating bug-catchers for that closure.

Why ship :meth:`solve` at all if it carries the WDD closure?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two reasons:

1. **The Wave D Round 2 unified sweep is bit-identically used as a
   preconditioner.**  Wave E's Krylov-on-apply needs an effective
   preconditioner; the sweep's :math:`O(N\cdot N_{\rm cells})`
   forward substitution is the canonical SN preconditioner (Lewis
   & Miller §4.5; Adams & Larsen 2002 review).  Exposing
   :meth:`solve` keeps that path discoverable through the
   capability set.

2. **The composers (Wave A) need a uniform contract.**  When a
   downstream agent composes :math:`(L - S - F)`, the composition's
   capability set is derived from each operand's capabilities.  If
   :class:`SNStreamingOperator` shipped only :meth:`apply`, the
   composition would lose ``solve`` and the Wave E Krylov-on-apply
   path could not request a sweep-preconditioned matvec without
   bypassing the algebra.

Reciprocity invariant
---------------------

The reciprocity identity is the defining property of the operator-
transpose pairing under the discrete L\ :sup:`2` inner product:

.. math::
    :label: sn-streaming-reciprocity

    \langle L\,\psi,\,\varphi\rangle \;=\;
    \langle\psi,\,L^*\,\varphi\rangle

for any pair :math:`(\psi, \varphi)` in the discrete unknown space
(packed solution vectors of length ``n_unknowns``).  Per Lewis &
Miller §10 (adjoint transport), this identity links forward and
adjoint sources / fluxes; it is the foundation on which detector
sensitivity, perturbation theory, and adjoint-weighted kinetics
all rest.

In :class:`SNStreamingOperator`, :meth:`apply_transpose` is built
from the explicit transpose of the dense matrix assembled by
probing :meth:`apply`, so the reciprocity identity holds **by
construction**: every linear operator has a transpose, and the
transpose of the assembled matrix *is* the operator's transpose.
The reciprocity test
:func:`tests.sn.test_snstreamingoperator.test_reciprocity_round_off`
gates this identity to round-off (``rtol=1e-12``, ``atol=1e-13``)
across slab, spherical, and cylindrical geometries on synthetic
:math:`(\psi, \varphi)` pairs.

A reciprocity-test failure would indicate one of two catastrophic
operator-correctness failures:

(a) :meth:`apply` is non-linear (a sign-flip, an incorrectly
    handled bias term, or a state-dependent operation such as
    re-normalisation) — caught by the linearity gate
    :func:`test_apply_is_linear`.
(b) The dense-assembly probe in
    :meth:`SNStreamingOperator._ensure_dense_matrix` does not
    faithfully assemble the matrix of :meth:`apply` — caught by
    the bit-identical extraction tests against
    :func:`transport_operator_matvec_*`.

Both gates are foundation-tagged in
:file:`tests/sn/test_snstreamingoperator.py`.

Vector layouts (``apply`` vs ``solve``)
---------------------------------------

:meth:`apply` and :meth:`apply_transpose` operate on the **packed
1-D solution vector** used by the legacy BiCGSTAB FD operator
path: an :class:`~orpheus.sn.operator.EquationMap` selects which
``(ordinate, cell)`` combinations are unknowns (the rest are
determined by reflective BCs and the z-hemisphere reduction); the
vector is laid out group-major in Fortran order.  This is the
natural input shape for :func:`scipy.sparse.linalg.bicgstab` and
the canonical layout for the Wave E Krylov-on-apply path.

:meth:`solve` operates on **structured arrays** matching
:func:`transport_sweep`'s contract:

* Source ``Q`` shape ``(nx, ny, ng)``.
* Persistent boundary-flux dict ``psi_bc`` carrying state between
  sweeps.
* Optional anisotropic source ``Q_aniso`` shape
  ``(N, nx, ny, ng)`` for P\ :sub:`ℓ` (:math:`\ell\ge 1`)
  scattering.

The shape mismatch reflects the historical layouts of the two
consumers (BiCGSTAB on packed-vectors / sweep on structured
arrays).  Wave E will normalise these via a single Krylov-on-apply
path; until then, calling :meth:`apply` and :meth:`solve` requires
the caller to be aware of the layout difference.

Why not extract :meth:`apply_transpose` analytically?
-----------------------------------------------------

For a continuous operator :math:`L = \Omega\cdot\nabla + \Sigma_t`,
the L\ :sup:`2`-adjoint is :math:`L^* = -\Omega\cdot\nabla +
\Sigma_t` (the streaming term flips sign under integration by
parts; the collision term is self-adjoint).  Discretising the
adjoint analytically would require:

* Reversing the upwind FD direction (:math:`\Omega \to -\Omega`)
  for the Cartesian streaming term.
* Transposing the τ-symmetric M-M angular interpolation closure.
* Re-deriving the adjoint of the per-level azimuthal redistribution
  for cylindrical.
* Handling adjoint boundary conditions (vacuum and reflective BCs
  are self-adjoint, so this term is benign in current ORPHEUS
  configurations — but a future ``albedo`` or ``albedo_white`` BC
  would need its own adjoint derivation).

Each of these analytical steps is a chance to introduce a sign
flip, a missing factor, or a transposed index — exactly the AI
failure modes the V&V framework names (modes #1, #2, #3, #5 in
``vv-principles``).  The dense-transpose path bypasses all of
them: the transpose of the discrete operator's matrix *is* the
operator's transpose, with no per-geometry derivation to verify.

The cost is :math:`O(n_{\rm unknowns}^2)` time and space for the
dense matrix, dominated by the :math:`n_{\rm unknowns}` calls to
:meth:`apply`.  For verification problems
(``n_unknowns ~ 30-1000``) the cost is negligible; for production-
scale problems (``n_unknowns ~ 10^4-10^6``) the dense path is
unsuitable.  Wave E will ship an :math:`O(n)` analytic-adjoint
matvec when production reciprocity becomes performance-critical
(currently it is not — :meth:`apply_transpose` is only consumed
by adjoint-flux post-processing in the Wave A operator algebra).

Wave E and beyond — landed and forward
--------------------------------------

* **Wave E Issue 14 lifted** the iteration primitives (power
  iteration, source iteration) into stand-alone operator-algebra
  consumers in :mod:`orpheus.numerics.iteration`, decoupling them
  from :class:`SNSolver`.
* **Wave E Issue 15 wired** :class:`SNSolver` to
  :class:`SNStreamingOperator`: the inner solve becomes Krylov on
  :meth:`apply` with sweep as preconditioner, which removes the WDD
  asymmetric closure from the converged-solution path for the
  reflective-BC eigenvalue case.  Wave E Round 3 will extend the
  equation-map layout to the vacuum-BC path that closes ERR-026
  for fixed-source curvilinear MMS.
* **Forward**: when production reciprocity becomes performance-
  critical, an :math:`O(n)` analytic-adjoint matvec replaces the
  dense-transpose fallback; the new path is gated by the same
  reciprocity tests in :file:`tests/sn/test_snstreamingoperator.py`.

References
----------

* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport.*  §10 (adjoint transport, reciprocity, and
  perturbation theory).
* Adams, M. L., & Larsen, E. W. (2002).  *Fast iterative methods
  for discrete-ordinates particle transport calculations.*
  Progress in Nuclear Energy 40 (1), 3–159.  Reviews
  Krylov-on-apply with sweep preconditioning.
* Trefethen, L. N., & Bau, D. (1997).  *Numerical Linear Algebra.*
  §3.2 (matrix-free Krylov view).


.. _sn-iteration-primitives:

Iteration primitives (operator algebra)
========================================

Wave E Round 1 (Issue #163) lifts the iteration primitives out of
:class:`~orpheus.sn.solver.SNSolver` into stand-alone operator-algebra
consumers in :mod:`orpheus.numerics.iteration`.  They consume the
Wave A :class:`~orpheus.numerics.operator.LinearOperator` Protocol
triple :math:`(L, S, F)` directly — no transport-solver knowledge
beyond the operator contract.

The :math:`(L - S - F)\,\psi = q_{\rm ext}` framing
----------------------------------------------------

The Boltzmann transport equation in its operator-algebra form
factors into three pieces (Lewis & Miller 1984 §6.4; Trefethen &
Bau 1997 §3.2 frame the matrix-free Krylov view):

.. math::

    (L - S - F)\,\psi = q_{\rm ext}
    \qquad\text{(within-group fixed source)}

.. math::

    (L - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`L = \Omega\cdot\nabla + \Sigma_t` is the streaming-
collision operator (see :ref:`sn-streaming-operator`),
:math:`S` is the scattering source operator
(see :class:`~orpheus.sn.scattering.ScatteringOperator`),
:math:`F = \chi\otimes\nu\Sigma_f` is the rank-1-in-energy fission
emission operator (see :class:`~orpheus.sn.fission.FissionOperator`),
and :math:`q_{\rm ext}` is an external source.

SourceIteration: discrete fixed-point realisation
--------------------------------------------------

:class:`~orpheus.numerics.iteration.SourceIteration` solves the
fixed-source equation by classical fixed-point iteration:

.. math::

    \psi_{n+1} \;=\; L^{-1}\,(S\,\psi_n + F\,\psi_n + q_{\rm ext}).

The convergence rate is bounded by the spectral radius
:math:`\rho(L^{-1}(S+F))`.  For an SN sweep applied to a homogeneous
infinite medium with isotropic scattering, this radius equals the
scattering ratio :math:`c = \Sigma_s/\Sigma_t` — convergence is
geometric at rate :math:`c` (Lewis & Miller §4.4).

The convergence test mirrors the legacy
:meth:`SNSolver._solve_source_iteration` exactly so Round 2's bit-
identical wiring is straightforward:

.. math::

    {\rm res}_n \;=\; \frac{\|\psi_n - \psi_{n-1}\|_2}
                            {\max(\|\psi_n\|_2,\,10^{-30})}

with the iteration breaking when :math:`{\rm res}_n < {\rm tol}`.

KEigenvalue: outer power iteration
-----------------------------------

:class:`~orpheus.numerics.iteration.KEigenvalue` solves the
eigenvalue problem by classical power iteration on the outer
:math:`k`-update, with :class:`SourceIteration` driving the inner
fixed-source solve at each step:

.. math::

    \psi_{n+1} \;=\; (L - S)^{-1}\,F\,\psi_n / k_n

.. math::

    k_{n+1} \;=\; {\rm keff\_estimator}(L, S, F, \psi_{n+1})

The dominance ratio :math:`|k_1/k_0|` governs outer-loop
convergence (Trefethen & Bau §27).  The inner solve uses
:class:`SourceIteration` with operator triple :math:`(L, S, 0)` —
the fission contribution at the inner level is the **external
source** :math:`F\psi_n/k_n`, NOT a within-group fixed-point term.
Every outer iteration warms up its inner :class:`SourceIteration`
from :math:`\psi_n` (the previous outer iterate); this is the same
amortisation pattern :class:`SNSolver` uses today.

The default ``keff_estimator`` is the generic Rayleigh quotient

.. math::

    k \;=\; \frac{\sum (F\,\psi)}{\sum (L\,\psi) - \sum (S\,\psi)}

which holds for any operator triple where the action carries the
volume measure.  SN consumers that need explicit volume weighting
(matching :meth:`SNSolver.compute_keff`) supply a custom
``keff_estimator`` callable; this is the load-bearing way Round 2
preserves bit-identity with the legacy
:func:`~orpheus.numerics.eigenvalue.power_iteration` path.

The ``inverter`` parameter — closing ERR-026
---------------------------------------------

Both primitives accept an ``inverter: Callable[[ndarray], ndarray]``
that supplies :math:`L^{-1}`.  When ``None``, the primitive routes
through :meth:`L.solve`.  When supplied, the caller controls how
:math:`L^{-1}` is realised — and this is the load-bearing design
choice for closing ERR-026.

* ``inverter = None`` (default):  :math:`L^{-1}\,q = L.solve(q)`.
  For an SN sweep this is the WDD asymmetric closure — which has a
  closure-bias-driven self-consistent fixed point on curvilinear
  meshes that is **not** the fine-mesh-limit transport solution
  (ERR-026).
* ``inverter = lambda q: gmres(as_scipy_linop(L), q, M=...)``:
  Krylov-on-:meth:`apply` (the symmetric closure of
  :class:`~orpheus.sn.operator.SNStreamingOperator`), with the
  sweep injected as a preconditioner :math:`M`.  This is the
  Wave E Round 2 reconciliation that closes ERR-026 for
  curvilinear SN: the converged solution comes from the symmetric
  :meth:`apply` closure (the one that agrees with analytical
  references) while the sweep accelerates the iteration as a
  preconditioner only — its closure bias does not poison the
  converged solution.

By making :math:`L^{-1}` a caller-supplied hook, the iteration
primitives do not need to be re-implemented when the inversion
strategy changes.  The same :class:`SourceIteration` runs in the
synthetic L0 case (where ``L`` is a plain dense matrix and
``inverter`` defaults to a direct solve), in the L1 SN case (where
``inverter`` defaults to the WDD sweep), and in the Wave E
Krylov-on-:meth:`apply` SN case (where ``inverter`` is supplied
explicitly by the caller).

Capability requirements
-----------------------

The primitive constructors enforce the following at construction
time, NEVER mid-iteration (the same Wave A philosophy that gates
:class:`~orpheus.numerics.operator.OperatorSum` etc.):

* ``L`` MUST advertise :pydata:`CAP_APPLY`.
* ``L`` MUST advertise :pydata:`CAP_SOLVE` *or* the caller MUST
  supply ``inverter``.  Without one of those, the iteration cannot
  evaluate :math:`L^{-1}`.
* ``S`` MUST advertise :pydata:`CAP_APPLY`.  Pass
  :class:`~orpheus.numerics.operator.ZeroOperator` for the
  scattering-free case.
* ``F`` MUST advertise :pydata:`CAP_APPLY`.  For
  :class:`SourceIteration` only, pass
  :class:`~orpheus.numerics.operator.ZeroOperator` for the
  fission-free case.

Constructor failure raises
:class:`~orpheus.numerics.operator.MissingCapability` with a message
naming the missing capability and the operand that lacks it.

Forward hook: FEAST and beyond
-------------------------------

:class:`KEigenvalue` accepts ``eigenvalue_method``, currently only
``"power"``.  The hook reserves a path for FEAST-style contour-
integral methods (Polizzi 2009) and Krylov-Schur deflation methods
(Stewart 2001) when accuracy on closely-spaced eigenvalues becomes
load-bearing.  Other values raise :class:`NotImplementedError` at
construction time.

Cross-references
----------------

* :class:`~orpheus.sn.operator.SNStreamingOperator` — the
  :math:`L` operand the SN solver ships, with both ``apply``
  (symmetric closure) and ``solve`` (WDD asymmetric closure).
  See :ref:`sn-streaming-operator` for the design rationale.
* :class:`~orpheus.sn.scattering.ScatteringOperator` — the
  :math:`S` operand carrying P\ :sub:`ℓ` scattering plus (n,2n).
* :class:`~orpheus.sn.fission.FissionOperator` — the rank-1-in-
  energy :math:`F` operand.  Returns :math:`F\,\phi` without the
  :math:`1/k` division (the eigenvalue scaling stays at the outer
  level, see the FissionOperator module docstring for the
  rationale).
* "SNSolver as an operator-algebra coordinator" — Round 2 will
  add this section once the SN solver is wired to consume the
  primitives directly.

.. _cross-solver-migration-sequence:

Cross-solver migration sequence
-------------------------------

The legacy :func:`orpheus.numerics.eigenvalue.power_iteration`
function and the :class:`EigenvalueSolver` Protocol are deprecated
for new SN code (the deprecation warning fires once at module
import).  They stay functional through the cross-solver migration
sequence:

* **CP** (collision-probability) — Issue TBD.  CP currently uses
  :func:`power_iteration` with its own
  ``EigenvalueSolver``-Protocol implementation.  Migration lifts
  CP onto an :math:`L`-equivalent collision-probability matrix
  operator + :class:`KEigenvalue`.
* **Diffusion** — Issue TBD.  The diffusion solver's BiCGSTAB
  inner loop is already a Krylov method; migration wraps it as an
  :math:`L^{-1}` ``inverter`` callable.
* **MoC** (method of characteristics) — Issue TBD.  MoC's track-
  based inner sweep maps onto :math:`L^{-1}` via the same
  facade pattern :class:`SNStreamingOperator.solve` uses.
* **Homogeneous** — Issue TBD.  Already a direct linear solve;
  migration is mostly cosmetic.

Each consumer's wave will land separately to keep regressions
isolated.  When the last consumer migrates, :func:`power_iteration`
and :class:`EigenvalueSolver` retire.

References for iteration primitives
------------------------------------

* Lewis, E. E., & Miller, W. F. (1984).  *Computational Methods of
  Neutron Transport.*  §4.4 (source iteration analysis); §6.4
  (operator-algebra view).
* Trefethen, L. N., & Bau, D. (1997).  *Numerical Linear Algebra.*
  §3.2 (matrix-free Krylov view); §27 (power iteration analysis).
* Polizzi, E. (2009).  *Density-matrix-based algorithm for solving
  eigenvalue problems.*  Phys. Rev. B 79, 115112.  The FEAST
  algorithm — the forward-hook target for ``eigenvalue_method``.
* Adams, M. L., & Larsen, E. W. (2002).  *Fast iterative methods
  for discrete-ordinates particle transport calculations.*
  Progress in Nuclear Energy 40 (1), 3–159.  Reviews
  Krylov-on-apply with sweep preconditioning (the ``inverter``
  hook's Round 2 use case).


.. _sn-solver-operator-algebra-coordinator:

SNSolver as an operator-algebra coordinator
============================================

Wave E Round 2 (Issue #164) closes the campaign loop by rewriting
:class:`~orpheus.sn.solver.SNSolver` to consume the operator triple
:math:`(L, S, F)` directly.  At construction time, the solver builds:

* :attr:`SNSolver.L` — :class:`SNStreamingOperator` carrying the
  symmetric-closure streaming-collision operator.
* :attr:`SNSolver.S` — :class:`~orpheus.sn.scattering.ScatteringOperator`
  carrying the P0 + (n,2n) + Pℓ Galerkin reconstruction (Wave D
  Issue 13).
* :attr:`SNSolver.F` — :class:`~orpheus.sn.fission.FissionOperator`
  carrying the rank-1-in-energy fission emission (Wave D Issue 13).

Each of these is a :class:`~orpheus.numerics.operator.LinearOperator`
in the Wave A operator-algebra sense: capability-tagged, composable
under :class:`~orpheus.numerics.operator.OperatorSum` and
:class:`~orpheus.numerics.operator.OperatorProduct`, and protocol-
conforming so the iteration primitives in
:mod:`orpheus.numerics.iteration` consume them without SN-specific
plumbing.

Adapter, not coordinator
------------------------

:class:`SNSolver` does NOT directly wrap
:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KEigenvalue` for the fixed-source
inner solve.  The reason is the Pℓ anisotropic scattering source:
``ScatteringOperator.build_aniso_source`` requires the **angular
flux of the previous iteration**, not the scalar flux that
:class:`SourceIteration` carries.  Threading that angular state
through the primitive's contract is a future cleanup; Round 2
preserves bit-identity by replicating :class:`SourceIteration`'s
loop structure verbatim inside :meth:`SNSolver._solve_source_iteration`
(the "Approach A" in the
``.claude/skills/algebra-of-record`` discipline — bit-identity now,
architectural cleanup follows).

The :meth:`SNSolver._solve_krylov` path likewise replicates the
Krylov-on-:meth:`apply` pattern inline rather than going through
:class:`SourceIteration` with a custom ``inverter`` hook, for the
same reason (Pℓ angular state) plus the fact that the GMRES outer
iteration has its own warm-start machinery (``self._psi_solution``).

The (L − S − F)·ψ = (1/k)·F·ψ framing at the solver level
-----------------------------------------------------------

Even without direct :class:`SourceIteration` consumption, the
:math:`(L, S, F)` framing organises the solver's API surface:

* :meth:`SNSolver.compute_fission_source` returns
  :math:`F\,\phi/k` — a thin delegator to :meth:`F.apply` with the
  :math:`1/k` outer-loop scaling applied at the solver level.
* :meth:`SNSolver.solve_fixed_source` solves
  :math:`(L - S)\,\psi = q_{\rm ext}` (with :math:`q_{\rm ext}` the
  fission source built by ``compute_fission_source``).  Two paths:

  * ``inner_solver="source_iteration"`` — sweep-driven fixed-point
    iteration; :math:`L^{-1}` is the WDD asymmetric sweep.  ERR-026-
    affected for curvilinear vacuum-BC cases.
  * ``inner_solver="krylov"`` — GMRES on :meth:`L.apply` with the
    sweep as preconditioner.  Reflective-BC equation map only
    (Round 3 owns the vacuum-BC extension).

* :meth:`SNSolver.compute_keff` returns
  :math:`\sum F\,\phi\,V / \sum \Sigma_a\,\phi\,V` (the volume-
  weighted production / absorption ratio); the existing math is
  preserved verbatim because
  :meth:`KEigenvalue.solve`'s default Rayleigh-quotient estimator
  is volume-blind.

The solver-level :math:`1/k` scaling and volume-weighting hooks
are exactly the points where SN's specifics live; the rest of the
solver is operator-algebra coordination.

The eigenvalue :math:`\keff` is determined by **power iteration**: an
outer loop updates :math:`k` from the production/absorption ratio, with
an inner loop that solves the within-group scattering problem.

Two Inner Solvers
-----------------

**Source iteration (sweep-based):**

- Operator: :math:`T^{-1}` (diamond-difference sweep)
- Solution variable: scalar flux :math:`\phi(x, y, g)`
- Fixed-point: :math:`\phi^{(k+1)} = T^{-1}(S \cdot \phi^{(k)} + Q_f)`
- Convergence rate: spectral radius of :math:`T^{-1}S` (~0.97 for 421
  groups)
- Cost per iteration: one transport sweep
- Works for all geometries

**Krylov (direct operator):**

- Operator: :math:`L = \mu \nabla + \Sigt{}` (finite-difference
  gradients, symmetric closure carried by
  :meth:`SNStreamingOperator.apply`)
- Solution variable: angular flux :math:`\psi(x, y, n, g)` (much
  larger than scalar flux)
- System: :math:`L\psi = b` where :math:`b` = fission + scattering
- Convergence: GMRES with sweep preconditioner, typically ~100
  Krylov iterations at ``tol=1e-4`` (always converges)
- Available for all geometries (Cartesian, spherical, cylindrical)
  with reflective outer-BC equation maps

Wave E Round 2 (Issue #164) replaced the legacy BiCGSTAB FD-operator
path with this Krylov path.  See `Krylov inner solver`_ above for
the full discussion.

The two architectures use **different spatial closures** (WDD
asymmetric sweep vs symmetric FD operator) that converge to
different :math:`\keff` on coarse meshes.  They agree in the limit
:math:`h \to 0`.


Verification
============

.. _sn-mms-verification:

Method of Manufactured Solutions (1D slab)
-------------------------------------------

Homogeneous and heterogeneous eigenvalue tests verify :math:`\keff`
--- a scalar. They do not tell us whether the **spatial operator**
itself converges at the design order :math:`\mathcal O(h^{2})` of
diamond difference.  The Method of Manufactured Solutions closes
that gap by constructing a fixed-source problem whose exact angular
flux is known in closed form, so the error against the prescribed
flux is pure spatial-discretisation error.

**Ansatz.**  For a vacuum-BC slab of length :math:`L` in one energy
group, pick an isotropic angular flux

.. math::
   :label: sn-mms-psi

   \psi_n(x) = \frac{1}{W}\,A(x),
   \qquad A(x) = \sin\!\left(\frac{\pi x}{L}\right),

where :math:`W = \sum_n w_n = 2` for Gauss--Legendre.  Because
:math:`A(0) = A(L) = 0`, every ordinate vanishes at both faces ---
the vacuum boundary conditions are satisfied automatically, with no
inflow bookkeeping required on the caller side.  Since :math:`\psi_n`
is independent of ordinate, the scalar flux recovered by any
quadrature order is *exactly* :math:`\phi(x) = A(x)` --- the test
isolates spatial error from angular quadrature error.

**Manufactured source.**  Substituting :eq:`sn-mms-psi` into the
discrete ordinates transport equation :eq:`transport-cartesian`
(with the :math:`1/W` convention ORPHEUS uses),

.. math::

   \mu_n\,\frac{\partial\psi_n}{\partial x} + \Sigma_t\,\psi_n
   = \frac{1}{W}\!\left(\Sigma_s\,\phi + Q^{\text{ext}}_n\right),

and solving algebraically for :math:`Q^{\text{ext}}_n` gives

.. math::
   :label: sn-mms-qext

   Q^{\text{ext}}_n(x)
   = \mu_n\,A'(x) + \bigl(\Sigma_t - \Sigma_s\bigr)\,A(x)
   = \mu_n\,\frac{\pi}{L}\cos\!\left(\frac{\pi x}{L}\right)
     + \bigl(\Sigma_t - \Sigma_s\bigr)\sin\!\left(\frac{\pi x}{L}\right).

The :math:`W` factor cancels cleanly because the ansatz was already
divided by :math:`W`, so what we hand the solver is the full residual
without any additional rescaling.  The expression is per-ordinate and
linear in :math:`\mu_n`: a constant isotropic external source *cannot*
drive a non-trivial manufactured flux because the streaming term
:math:`\mu_n\,\psi'_n` is odd in :math:`\mu`.  That is the fundamental
reason MMS for SN requires the :math:`Q_{\rm aniso}` plumbing path ---
no "cheat" with a cell-by-cell isotropic source exists.

**Why :math:`\sin(\pi x/L)`?**  The ansatz is smooth
(:math:`C^{\infty}`) so all derivatives of the exact solution exist
and DD's :math:`\mathcal O(h^{2})` truncation error dominates.  It
vanishes at both boundaries for free.  Its derivatives do not collapse
to a polynomial --- a cubic ansatz, for instance, has a constant
second derivative so the DD truncation term :math:`\psi'''` would be
zero and the error could disappear for a non-physical reason,
hiding bugs.  Trigonometric or exponential ansätze have bounded
but non-zero derivatives of every order and therefore expose the
leading truncation term cleanly.

**Implementation.**  The case is built by
:func:`orpheus.derivations.continuous.mms.sn.build_1d_slab_mms_case` and
consumed by :func:`orpheus.sn.solve_sn_fixed_source`.  The latter
accepts a per-ordinate external source of shape
:math:`(N, n_x, n_y, n_g)` and threads it through the sweep's
:math:`Q_{\rm aniso}` slot --- merging additively with any P1+
scattering contribution the solver itself builds.  Vacuum boundary
conditions are applied via the mesh-level BC infrastructure
described in :ref:`boundary-conditions`:
:func:`solve_sn_fixed_source` defaults its ``boundary_condition``
parameter to ``"vacuum"`` and the internal helper
``_apply_default_bcs`` stamps :attr:`BC.vacuum <geometry.mesh.BC.vacuum>`
onto every face of the mesh that lacks an explicit BC declaration.
:class:`SNMesh` then resolves these to the ``"vacuum"`` kind string,
which the sweep reads directly.  In the 1-D cumprod path, the
recurrence starts from zero; in the 2-D wavefront path, the
reflective-partner copy is skipped, leaving incoming-face angular
fluxes at their zero initialisation (which is correct because no
code path writes the incoming-face slot of any ordinate except the
reflection step itself).

.. note::

   Before the BC infrastructure was introduced,
   :func:`~orpheus.sn.sweep.transport_sweep` accepted a
   ``boundary_condition: str`` parameter directly.  That parameter
   has been removed --- BCs now flow through the mesh → SNMesh
   resolution path.  The description above reflects the current
   implementation.

**Measured convergence.**  With
:math:`\Sigma_t = 1\ \mathrm{cm^{-1}}`,
:math:`\Sigma_s = 0.5\ \mathrm{cm^{-1}}`,
:math:`L = 5\ \mathrm{cm}`, Gauss--Legendre :math:`S_{16}`:

.. list-table::
   :header-rows: 1
   :widths: 10 20 20

   * - :math:`n_{\rm cells}`
     - :math:`\|\phi_h - \phi_{\rm ex}\|_{L^{2}}`
     - measured order
   * - 10
     - :math:`2.17\!\times\!10^{-3}`
     - ---
   * - 20
     - :math:`5.40\!\times\!10^{-4}`
     - 2.01
   * - 40
     - :math:`1.35\!\times\!10^{-4}`
     - 2.00
   * - 80
     - :math:`3.37\!\times\!10^{-5}`
     - 2.00
   * - 160
     - :math:`8.42\!\times\!10^{-6}`
     - 2.00

Successive ratios hit :math:`4.00\pm0.02`, i.e. the measured order
is exactly the design order of diamond difference.  The L1 test
:func:`tests.sn.test_mms.test_sn_1d_slab_mms_converges_second_order`
asserts a slightly loose ``order > 1.9`` bracket to leave room for
round-off at the finest mesh.

**Risk points / things that can go wrong.**

- *Vacuum BC not honoured.*  If the reflective-partner copy is not
  skipped, incoming-face angular flux at the boundary is non-zero
  (the reflected outgoing from the opposite sweep) and the
  manufactured solution no longer satisfies the discrete problem.
  Symptom: :math:`\mathcal O(1)` error at the coarsest mesh; no
  convergence regardless of refinement.
- *Wrong normalisation for* :math:`Q_{\rm ext}`.  The solver's
  :math:`Q_{\rm aniso}` slot is divided by :math:`W` internally;
  the ansatz has a :math:`1/W` prefactor; the two must cancel.
  If the derivation forgets the :math:`W` cancellation, the
  measured flux is a factor of :math:`W` off but still converges at
  order 2 --- sneaky.  Guard: the second test in ``test_mms.py``
  cross-checks the algebraic symmetry of :eq:`sn-mms-qext`.
- *Non-smooth ansatz.*  A discontinuous material or a piecewise
  linear ansatz degrades the observed order to :math:`\mathcal O(h)`.
  The homogeneous sinusoid avoids both.
- *1-group vs multigroup.*  Because the manufactured flux is isotropic
  and there is no fission in the fixed-source problem, 1 group is
  sufficient --- the degeneracy warning about 1-group eigenvalue
  tests does not apply, since no :math:`\keff` enters.  Multigroup
  and heterogeneous MMS extensions are tracked as follow-ups for
  richer operator coverage.

**Follow-ups.**  MMS for :doc:`method_of_characteristics`, diffusion,
and spherical / cylindrical curvilinear SN is tracked in GitHub
Issues (see ``type:feature level:L1``).  The curvilinear sweeps
need their own ansatz because their vacuum BC plumbing is not
yet wired up. **Heterogeneous and multigroup SN MMS is covered
by the next subsection.**


.. _sn-mms-heterogeneous-verification:

Heterogeneous MMS — 2-group continuous-:math:`\Sigma` slab
-----------------------------------------------------------

The homogeneous MMS case above verifies the Cartesian 1D SN
sweep for a *single-material* slab. To verify the multigroup
operator on a **heterogeneous** problem --- where each cell can
have different cross sections and the scatter matrix couples
groups across positions --- the Method of Manufactured Solutions
is extended in Phase 2.1a of the verification campaign with two
deliberate choices:

1. **Continuous (smooth)** :math:`\Sigma_{t,g}(x)` and
   :math:`\Sigma_{s,g\to g'}(x)` instead of piecewise-constant
   material regions. Discontinuous :math:`\Sigma` at interfaces
   that do not coincide with cell faces degrades diamond
   difference from :math:`\mathcal O(h^{2})` to
   :math:`\mathcal O(h)`, which would contaminate the
   spatial-convergence measurement with interface-treatment
   artefacts. With smooth :math:`\Sigma(x)` the diamond-
   difference operator hits its design :math:`\mathcal O(h^{2})`
   order exactly --- the convergence study becomes a clean test
   of the operator itself. This follows Salari & Knupp
   SAND2000-1444 §6, the canonical MMS reference for
   heterogeneous verification.
2. **Per-group amplitudes** :math:`\mathbf c = (c_1, c_2)` in
   the ansatz, so the scalar flux has a non-trivial group
   spectrum and the downscatter source term in the manufactured
   :math:`Q^{\text{ext}}` is non-zero. A bug that transposes
   the scatter matrix or drops a cross-group source term
   produces an incorrect :math:`\phi_2` that the convergence
   test catches immediately.

**Ansatz.**  The homogeneous ansatz carries over, now with a
per-group amplitude:

.. math::
   :label: sn-mms-hetero-psi

   \psi_{n,g}(x) \;=\; \frac{c_g}{W}\,A(x),
   \qquad A(x) \;=\; \sin\!\left(\frac{\pi x}{L}\right),

where :math:`W = \sum_n w_n` is the quadrature weight sum. The
scalar flux in each group is
:math:`\phi_g(x) = c_g\,A(x)`, so the amplitudes
:math:`\mathbf c` literally are the group fluxes at the slab
midpoint (where :math:`A` peaks). With
:math:`\mathbf c = (1.0, 0.3)` the two groups are linearly
independent and the downscatter coupling is visible.

Both groups share the same *spatial* mode :math:`\sin(\pi x/L)`
--- this is the fundamental mode of the bare slab and is exactly
the shape that emerges from separation of variables in the
diffusion limit. The heterogeneous SN problem would in general
have each group living in its own spatial harmonic, but we
*choose* the shared-mode ansatz as the manufactured target and
derive the non-trivial :math:`Q^{\text{ext}}` that makes it
satisfy the transport equation. The test then measures how well
the numerical SN sweep reproduces this prescribed shape.

**Manufactured source.**  Substituting :eq:`sn-mms-hetero-psi`
into the multigroup discrete-ordinates transport equation

.. math::

    \mu_n\,\frac{\partial\psi_{n,g}}{\partial x}
        + \Sigma_{t,g}(x)\,\psi_{n,g}
    \;=\; \frac{1}{W}\!\left(
        \sum_{g'}\Sigma_{s,g'\to g}(x)\,\phi_{g'}(x)
      + Q^{\text{ext}}_{n,g}(x)
    \right)

and solving algebraically for :math:`Q^{\text{ext}}`:

.. math::
   :label: sn-mms-hetero-qext

   Q^{\text{ext}}_{n,g}(x) \;=\;
       \mu_n\,c_g\,A'(x)
     + c_g\,\Sigma_{t,g}(x)\,A(x)
     \;-\; \sum_{g'}\Sigma_{s,g'\to g}(x)\,c_{g'}\,A(x).

The :math:`W` factor cancels between the ansatz's :math:`1/W`
prefactor and the solver's own :math:`1/W` convention on the
isotropic and anisotropic source slots, so :eq:`sn-mms-hetero-qext`
is the residual hand-delivered to the sweep without any
additional rescaling.

**Structure of the source.**  The streaming term
:math:`\mu_n\,c_g\,A'(x)` is odd in :math:`\mu` and carries the
only angular dependence, which is why SN MMS fundamentally
needs the per-ordinate ``Q_aniso`` plumbing path. The removal
term :math:`c_g\,\Sigma_{t,g}(x)\,A(x)` is diagonal in group
index. The **in-scatter** sum
:math:`\sum_{g'}\Sigma_{s,g'\to g}\,c_{g'}\,A(x)` is the only
term that couples groups, and for :math:`g=2` in the default
2-group setup it contributes
:math:`-\Sigma_{s,1\to 2}(x)\,c_1\,A(x)` --- the thermal source
depends on the fast amplitude through the downscatter cross
section, exactly as the multigroup scatter assembly in the
sweep does.

**Canonical cross sections.**  The reference uses smooth
profiles on :math:`[0, L]`:

.. math::

    \Sigma_{t,1}(x) &= 1.0 + 0.2\sin(\pi x/L), \\
    \Sigma_{t,2}(x) &= 2.0 + 0.3\cos(\pi x/L), \\
    \Sigma_{s,1\to 1}(x) &= 0.3 + 0.1\sin(\pi x/L), \\
    \Sigma_{s,1\to 2}(x) &= 0.2 + 0.05\sin(\pi x/L), \\
    \Sigma_{s,2\to 2}(x) &= 1.5 + 0.15\sin(\pi x/L), \\
    \Sigma_{s,2\to 1}(x) &= 0.

These give :math:`\Sigma_{a,1}(x) = 0.5 + 0.05\sin(\pi x/L) > 0`
trivially and
:math:`\Sigma_{a,2}(x) = 0.5 + 0.3\cos(\pi x/L) - 0.15\sin(\pi x/L)`,
bounded below by :math:`0.5 - \sqrt{0.3^{2} + 0.15^{2}} \approx
0.165 > 0`, so the cross sections are physical everywhere. The
scattering ratios :math:`c_g = \Sigma_{s,\text{tot},g}/\Sigma_{t,g}`
stay around :math:`0.5` for both groups, which means source
iteration converges geometrically at rate :math:`\sim 0.5^n`
per sweep.

**Per-cell material construction.**  The solver consumes the
continuous :math:`\Sigma(x)` by creating **one material per cell**
with cross sections evaluated at the cell centre
:math:`x_i = (x_{i-1/2} + x_{i+1/2})/2`. The midpoint rule for
the cell-average cross section is :math:`\mathcal O(h^{2})`-
accurate on smooth :math:`\Sigma`, matching the diamond-
difference design order and not degrading the measured
convergence rate. The number of materials scales with mesh
refinement, so each mesh in the convergence study builds a
fresh materials dictionary via
:meth:`orpheus.derivations.continuous.mms.sn.SNSlab2GHeterogeneousMMSCase.build_materials`.

**Measured convergence.**  With default parameters
(:math:`L = 5\,\text{cm}`, :math:`\mathbf c = (1.0, 0.3)`,
Gauss--Legendre :math:`S_{16}`):

.. list-table::
   :header-rows: 1
   :widths: 10 20 20 20

   * - :math:`n_{\rm cells}`
     - :math:`\|\phi_1 - \phi_{1,\rm ex}\|_{L^{2}}`
     - :math:`\|\phi_2 - \phi_{2,\rm ex}\|_{L^{2}}`
     - measured order
   * - 20
     - :math:`3.71\!\times\!10^{-4}`
     - :math:`3.38\!\times\!10^{-4}`
     - ---
   * - 40
     - :math:`9.25\!\times\!10^{-5}`
     - :math:`8.45\!\times\!10^{-5}`
     - 2.00
   * - 80
     - :math:`2.31\!\times\!10^{-5}`
     - :math:`2.11\!\times\!10^{-5}`
     - 2.00
   * - 160
     - :math:`5.78\!\times\!10^{-6}`
     - :math:`5.28\!\times\!10^{-6}`
     - 2.00

Both groups hit the design order independently, confirming
that the multigroup scatter coupling is correctly exercised.
The L1 test
:func:`tests.sn.test_mms_heterogeneous.test_sn_heterogeneous_mms_converges_second_order`
asserts ``> 1.9`` to leave round-off headroom at the finest
mesh.

**What this replaces.** Before Phase 2.1a, the heterogeneous
SN verification was
:func:`orpheus.derivations.continuous.cases.sn._derive_sn_heterogeneous`, which
computed the reference :math:`k_{\text{eff}}` by running the
SN solver itself at four mesh refinements and Richardson-
extrapolating the eigenvalue sequence. That is a **T3 circular
self-test** in the verification-campaign taxonomy: the solver
verifies against its own extrapolated output, so any consistent
bug in the SN sweep that affects all mesh refinements the same
way is invisible to the test. The heterogeneous MMS reference
above breaks the circularity: the reference comes from the
manufactured-solution algebra, not from the solver.

**Complementary eigenvalue verification.** The MMS test
verifies the **spatial operator** on a heterogeneous problem
but does not exercise the eigenvalue iteration. Phase 2.1b
lands a Case singular-eigenfunction eigenvalue reference --- see
:ref:`sn-case-heterogeneous-verification` --- that restores
eigenvalue-heterogeneous coverage for the SN solver (T2
semi-analytical, from the first-order Boltzmann equation
itself, no diffusion approximation).


.. _sn-mms-2d-verification:

2D Cartesian MMS — separable sinusoidal ansatz
-----------------------------------------------

Phase 3.1 of the verification campaign extends the MMS spatial-operator
verification to **two Cartesian dimensions**.  The 1D slab MMS tests
verify the :math:`\mu\,\partial\psi/\partial x` streaming term in
isolation; this section adds :math:`\mu_y\,\partial\psi/\partial y`
and confirms that the 2D wavefront sweep
(:func:`orpheus.sn.sweep._sweep_2d_wavefront`) with diamond-difference
closure achieves its design :math:`\mathcal O(h^{2})` convergence rate.

**Ansatz.**  On a rectangle :math:`[0, L_x] \times [0, L_y]` with
vacuum boundary conditions:

.. math::
   :label: sn-mms-2d-psi

   \psi_n(x, y) \;=\; \frac{1}{W}\,A(x, y),
   \qquad A(x, y) \;=\; \sin\!\left(\frac{\pi x}{L_x}\right)
                         \sin\!\left(\frac{\pi y}{L_y}\right).

The ansatz is **isotropic in angle** --- every ordinate carries the
same angular flux amplitude --- so the scalar flux recovered by any
quadrature set equals :math:`\phi(x, y) = A(x, y)` exactly.  This
design is deliberate: it isolates **spatial** discretisation error from
angular quadrature error, exactly as in the 1D case
(:eq:`sn-mms-psi`).

The separable sinusoidal ansatz vanishes on all four domain edges
(:math:`x = 0`, :math:`x = L_x`, :math:`y = 0`, :math:`y = L_y`),
so vacuum BCs are satisfied automatically for every ordinate.

**Manufactured source.**  Substituting :eq:`sn-mms-2d-psi` into the
2D Cartesian transport equation :eq:`transport-cartesian-2d` and
solving for the residual:

.. math::
   :label: sn-mms-2d-qext

   Q^{\text{ext}}_n(x, y) \;=\;
       \mu_{x,n}\,\frac{\partial A}{\partial x}
     + \mu_{y,n}\,\frac{\partial A}{\partial y}
     + (\Sigma_t - \Sigma_s)\,A(x, y)

where the partial derivatives of the separable ansatz are:

.. math::

   \frac{\partial A}{\partial x} =
       \frac{\pi}{L_x}\cos\!\left(\frac{\pi x}{L_x}\right)
       \sin\!\left(\frac{\pi y}{L_y}\right), \qquad
   \frac{\partial A}{\partial y} =
       \sin\!\left(\frac{\pi x}{L_x}\right)
       \frac{\pi}{L_y}\cos\!\left(\frac{\pi y}{L_y}\right).

The manufactured source :eq:`sn-mms-2d-qext` is angle-dependent through
:math:`\mu_{x,n}` and :math:`\mu_{y,n}` (streaming terms) and
angle-independent in the removal term :math:`(\Sigma_t - \Sigma_s) A`.
It enters the solver through the ``Q_aniso`` external source slot in
:func:`orpheus.sn.solve_sn_fixed_source`.

**Quadrature.**  2D problems use Lebedev spherical quadrature
(:class:`orpheus.sn.quadrature.LebedevSphere`, order 17 = 110 ordinates).
Because the ansatz is isotropic in angle, the quadrature-level angular
integration is exact for *any* quadrature set --- the spatial
convergence study isolates spatial error exclusively.

**Measured convergence.**  Four mesh refinements on a
:math:`5 \times 5\,\text{cm}` square domain with
:math:`\Sigma_t = 1.0`, :math:`\Sigma_s = 0.5`:

.. list-table::
   :header-rows: 1

   * - :math:`n_x = n_y`
     - L2 error
     - Order
   * - 10
     - :math:`5.50 \times 10^{-3}`
     -
   * - 20
     - :math:`1.37 \times 10^{-3}`
     - 2.01
   * - 40
     - :math:`3.41 \times 10^{-4}`
     - 2.00
   * - 80
     - :math:`8.53 \times 10^{-5}`
     - 2.00

The measured order is indistinguishable from 2.00 across all
refinements, confirming that the 2D wavefront sweep preserves the
diamond-difference design order.

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SN2DCartesianMMSCase` and
  :func:`orpheus.derivations.continuous.mms.sn.build_2d_cartesian_mms_case`.
- Test:
  :func:`tests.sn.test_mms.test_sn_2d_cartesian_mms_converges_second_order`.
- Sweep:
  :func:`orpheus.sn.sweep._sweep_2d_wavefront` (the 2D diamond-difference
  kernel verified by this test).

**Why this test matters.**  The existing 2D SN tests
(:mod:`tests.sn.test_discrete_ordinates_2d`) are L2 self-convergence
tests with real cross sections that verify the solver as a black box.
This MMS test is more incisive: it provides a **closed-form reference
flux** and asserts the **design convergence order** of the spatial
discretisation.  A bug that corrupts the 2D DD cell-average formula
(e.g. swapping :math:`\Delta x` and :math:`\Delta y`, mis-indexing the
wavefront anti-diagonal, or computing face fluxes with the wrong
sign) would break the :math:`\mathcal O(h^{2})` rate while possibly
still converging at some reduced order — the MMS test catches this
immediately, while a self-convergence test might not.

**Gotchas.**

- *Ordinates with* :math:`\mu_x = \mu_y = 0`.  The Lebedev set
  includes purely :math:`z`-directed ordinates.  For these, the
  streaming terms vanish, and the sweep reduces to
  :math:`\psi = Q/\Sigma_t`.  The manufactured source formula
  handles this correctly because both :math:`\mu_{x,n}` and
  :math:`\mu_{y,n}` multiply the gradient terms.
- *Aspect ratio.*  The test uses :math:`L_x = L_y` (square domain).
  A non-square domain would work identically — the separable ansatz
  is parameterised by :math:`L_x` and :math:`L_y` independently.
  Phase 3.2 extends to 2-group with heterogeneous materials (below).


.. _sn-mms-2d-2g-verification:

2D Cartesian 2-group heterogeneous MMS
----------------------------------------

Phase 3.2 combines the 2D geometry from Phase 3.1 with the
smooth-:math:`\Sigma` heterogeneous approach from Phase 2.1a.  The
cross sections are smooth 2D functions :math:`\Sigma(x, y)` so the
diamond-difference design order :math:`\mathcal O(h^{2})` is preserved
(no interface degradation).

**Ansatz.**  Per-group amplitudes :math:`c_g` with the same 2D shape:

.. math::
   :label: sn-mms-2d-2g-psi

   \psi_{n,g}(x, y) = \frac{c_g}{W}\,A(x, y), \qquad
   A(x, y) = \sin(\pi x/L_x)\,\sin(\pi y/L_y),

giving :math:`\phi_g(x, y) = c_g\,A(x, y)` with
:math:`\mathbf c = (1.0, 0.3)`.

**Manufactured source.**  From the 2D multigroup transport equation:

.. math::
   :label: sn-mms-2d-2g-qext

   Q^{\text{ext}}_{n,g}(x, y) =
       \mu_{x,n}\,c_g\,\partial_x A
     + \mu_{y,n}\,c_g\,\partial_y A
     + \Sigma_{t,g}(x, y)\,c_g\,A
     - \sum_{g'}\Sigma_{s,g'\to g}(x, y)\,c_{g'}\,A.

The thermal (:math:`g = 2`) source couples to :math:`c_1` through
the downscatter term :math:`\Sigma_{s,1\to 2}(x, y)\,c_1\,A`, which
exercises the multigroup scatter assembly in the 2D sweep.

**Cross-section profiles.**  The 2D functions extend the 1D
Phase-2.1a profiles (see :ref:`sn-mms-heterogeneous-verification`)
with a mild :math:`y`-dependent modulation:

- :math:`\Sigma_{t,1}(x,y) = 1.0 + 0.2\sin(\pi x/L_x) + 0.1\cos(\pi y/L_y)`
- :math:`\Sigma_{t,2}(x,y) = 2.0 + 0.3\cos(\pi x/L_x) + 0.1\sin(\pi y/L_y)`

Scattering cross sections carry a :math:`0.05\cos(\pi y/L_y)` modulation.
All :math:`\Sigma_a > 0` bounds from the 1D case are preserved because
the :math:`y`-modulation amplitudes (0.1, 0.05) are smaller than the
1D absorption margin (:math:`\sim 0.165`).

**Measured convergence.**  Four refinements on a :math:`5 \times 5` cm
square:

.. list-table::
   :header-rows: 1

   * - :math:`n_x = n_y`
     - L2 error (g=1)
     - Order (g=1)
     - L2 error (g=2)
     - Order (g=2)
   * - 10
     - :math:`3.79 \times 10^{-3}`
     -
     - :math:`2.85 \times 10^{-3}`
     -
   * - 20
     - :math:`9.41 \times 10^{-4}`
     - 2.01
     - :math:`7.09 \times 10^{-4}`
     - 2.01
   * - 40
     - :math:`2.35 \times 10^{-4}`
     - 2.00
     - :math:`1.77 \times 10^{-4}`
     - 2.00
   * - 80
     - :math:`5.87 \times 10^{-5}`
     - 2.00
     - :math:`4.42 \times 10^{-5}`
     - 2.00

Both groups achieve the design :math:`\mathcal O(h^{2})` rate.

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SN2DCartesian2GHeterogeneousMMSCase`
  and :func:`orpheus.derivations.continuous.mms.sn.build_2d_cartesian_heterogeneous_mms_case`.
- Test:
  :func:`tests.sn.test_mms_2d.test_sn_2d_cartesian_2g_heterogeneous_mms_converges_second_order`.


.. _sn-mms-p1-verification:

P1 anisotropic scattering MMS
-------------------------------

Phase 3.5 verifies that the P\ :sub:`N` anisotropic scattering
source assembly (:ref:`pn-scattering`) preserves
:math:`\mathcal O(h^{2})` convergence. All previous MMS tests use
isotropic (P0) scattering; this test exercises the P1 slot
:math:`\Sigma_s^{(1)}` through a weakly angle-dependent ansatz.

**Ansatz.** On a 1D vacuum-BC slab :math:`[0, L]`:

.. math::
   :label: sn-mms-p1-psi

   \psi_n(x) = \frac{1}{W}\bigl(A(x) + \alpha\,\mu_n\,B(x)\bigr)

with :math:`A(x) = B(x) = \sin(\pi x/L)` and small
:math:`\alpha = 0.1`. The scalar flux is :math:`\phi(x) = A(x)`
(the :math:`\mu`-odd term integrates to zero), and the P1 current
is :math:`J(x) = \alpha\,B(x)/3` (using
:math:`\sum w_n\mu_n^2 = 2/3` for Gauss–Legendre on
:math:`[-1, 1]`).

**Manufactured source.** Substituting :eq:`sn-mms-p1-psi` into
the 1D transport equation with P1 scattering and solving for
the residual:

.. math::
   :label: sn-mms-p1-qext

   Q^{\text{ext}}_n(x) =
       \mu_n\,A'(x)
     + (\Sigma_t - \Sigma_s^{(0)})\,A(x)
     + \alpha\,\mu_n\,(\Sigma_t - \Sigma_s^{(1)})\,B(x)
     + \alpha\,\mu_n^2\,B'(x).

The first two terms are the isotropic MMS source from
:eq:`sn-mms-qext`. The third term comes from the P1 scattering
slot :math:`3\,\Sigma_s^{(1)}\,\mu_n\,J(x)` in the transport
equation, and the fourth from the :math:`\mu_n`-weighted
streaming of :math:`B(x)`.

**Measured convergence.** Four refinements with
:math:`\Sigma_t = 1.0`, :math:`\Sigma_s^{(0)} = 0.5`,
:math:`\Sigma_s^{(1)} = 0.2`, :math:`\alpha = 0.1`:

.. list-table::
   :header-rows: 1

   * - :math:`n_{\text{cells}}`
     - L2 error
     - Order
   * - 20
     - :math:`6.15 \times 10^{-4}`
     -
   * - 40
     - :math:`1.53 \times 10^{-4}`
     - 2.00
   * - 80
     - :math:`3.84 \times 10^{-5}`
     - 2.00
   * - 160
     - :math:`9.59 \times 10^{-6}`
     - 2.00

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SNP1AnisoMMSCase` and
  :func:`orpheus.derivations.continuous.mms.sn.build_p1_aniso_mms_case`.
- Test:
  :func:`tests.sn.test_mms_aniso.test_sn_p1_aniso_mms_converges_second_order`.
- P1 assembly:
  :meth:`orpheus.sn.solver.SNSolver._build_aniso_scattering`.


.. _sn-mms-curvilinear-aniso-verification:

Curvilinear anisotropic MMS — angular redistribution probe
-----------------------------------------------------------

Phase 3.6 closes the **angular-redistribution coverage gap** in the
curvilinear MMS verification chain. The existing isotropic
1D-spherical (:class:`SNSphericalMMSCase`) and 1D-cylindrical
(:class:`SNCylindricalMMSCase`) MMS cases use the ansatz
:math:`\psi_n(r) = A(r)/W` (no :math:`\mu`-dependence). For that
ansatz, **the angular-redistribution operator is identically zero**:
:math:`(1-\mu^2)/r \cdot \partial\psi/\partial\mu = 0` for the sphere,
:math:`-(1/r)\,\partial(\xi\psi)/\partial\varphi = 0` for the
cylinder. The hardest math the curvilinear sweep performs — where
ERR-026 (curvilinear sweep WDD wrong fixed point) lives — is
mathematically invisible to the isotropic MMS because it cancels by
ansatz construction.

This is the ``vv-principles`` failure mode #7 ("MMS simplification
bias") — the MMS test cannot catch a bug class because the ansatz
nulls it. The defence is **not**
to replace the isotropic case (it remains the right probe for the
non-redistribution paths) but to **pair** it with a companion case
whose ansatz activates redistribution. The two cases together let a
narrow-down diagnosis route a failing convergence rate to either
the streaming/removal path (only isotropic fails) or the
redistribution path (only anisotropic fails).

**Spherical anisotropic ansatz**

.. math::
   :label: sn-mms-spherical-aniso-psi

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + B(r)\,\mu_n\bigr),
   \qquad
   A(r) = \sin\!\left(\frac{\pi r}{R}\right),
   \qquad
   B(r) = \frac{r}{R}\Bigl(1 - \frac{r}{R}\Bigr)
            \cos\!\left(\frac{\pi r}{R}\right).

Both :math:`A` and :math:`B` vanish at :math:`r \in \{0, R\}`, so
**every** ordinate satisfies the symmetry BC at :math:`r=0` and the
vacuum BC at :math:`r=R`, regardless of the sign of :math:`\mu_n`.
The :math:`B(r)\,\mu_n` coefficient is non-trivial: ordinates with
opposite sign of :math:`\mu_n` differ in sign of the
angular-flux contribution, but both still vanish at the boundaries.

The choice :math:`B(r) = (r/R)(1-r/R)\cos(\pi r/R)` is **not**
algebraically reducible to a multiple of :math:`A(r)` — the
:math:`(r/R)(1-r/R)` envelope and the :math:`\cos(\pi r/R)` factor
produce a derivative :math:`B'(r)` whose extrema do not co-locate
with :math:`A'(r)`'s extrema, so the redistribution term
:math:`(1-\mu_n^2)\,B/r` cannot be absorbed into a renormalisation
of the streaming term.

The discrete scalar flux is :math:`\phi(r) = A(r)` because
:math:`\sum_n w_n \mu_n = 0` for any symmetric Gauss-Legendre
quadrature — the :math:`B \mu` term integrates to zero in the
scalar moment.

**Spherical manufactured source**

Substituting :eq:`sn-mms-spherical-aniso-psi` into
:eq:`transport-spherical` and solving for the residual:

.. math::
   :label: sn-mms-spherical-aniso-qext

   Q^{\text{ext}}_n(r) =
        \mu_n\,A'(r)
      + \mu_n^2\,B'(r)
      + (1 - \mu_n^2)\,\frac{B(r)}{r}
      + (\Sigma_t - \Sigma_s)\,A(r)
      + \Sigma_t\,\mu_n\,B(r).

The first and fourth terms are the isotropic-MMS source from
:eq:`sn-mms-qext` adapted to spherical. The **second term**
(:math:`\mu_n^2 B'(r)` — :math:`\mu`-weighted streaming of the
anisotropic profile) and the **third term**
(:math:`(1-\mu_n^2)\,B/r` — angular redistribution) are
load-bearing: they are precisely what the isotropic case lacks.
The fifth (:math:`\Sigma_t\,\mu_n B`) comes from the removal
operator acting on the :math:`B \mu` part of :math:`\psi_n`.

**Cylindrical anisotropic ansatz**

The radial direction cosine for cylindrical 1D is :math:`\eta_n =
\sin\theta_n \cos\varphi_n`; the azimuthal partner that drives
the redistribution is :math:`\xi_n = \sin\theta_n \sin\varphi_n`.
Use:

.. math::
   :label: sn-mms-cylindrical-aniso-psi

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + B(r)\,\eta_n\bigr),

with the same :math:`A(r),\,B(r)` shapes. Symmetric ProductQuadrature
gives :math:`\sum_n w_n \eta_n = 0`, so :math:`\phi(r) = A(r)`.

**Cylindrical manufactured source**

Substituting :eq:`sn-mms-cylindrical-aniso-psi` into
:eq:`transport-cylindrical` (treating :math:`\eta_n` and
:math:`\xi_n` as the :math:`\varphi`-dependent functions
:math:`\sin\theta\cos\varphi` and :math:`\sin\theta\sin\varphi`)
and solving for the residual:

.. math::
   :label: sn-mms-cylindrical-aniso-qext

   Q^{\text{ext}}_n(r) =
        \eta_n\,A'(r)
      + \eta_n^2\,B'(r)
      + \xi_n^2\,\frac{B(r)}{r}
      + (\Sigma_t - \Sigma_s)\,A(r)
      + \Sigma_t\,\eta_n\,B(r).

The :math:`\xi_n^2\,B/r` redistribution term is the cylindrical analog
of the sphere's :math:`(1-\mu_n^2)\,B/r`. Both come from the
same operator — angular redistribution of the linearly-:math:`\mu`
(or linearly-:math:`\eta`) ansatz — and both vanish for any
isotropic ansatz.

**Verification chain (Branch 1 / Branch 2)**

Per the ``algebra-of-record`` discipline (Branch-1 SymPy reference,
Branch-2 numpy production, structurally-independent L1 cross-check):

- **Branch 1 (SymPy)**:
  :func:`orpheus.derivations.continuous.mms.sn.derive_spherical_anisotropic_mms`
  and
  :func:`orpheus.derivations.continuous.mms.sn.derive_cylindrical_anisotropic_mms`
  substitute the ansatz into the transport operator symbolically and
  prove ``simplify(LHS - RHS) == 0``. Foundation tests:
  :file:`tests/derivations/test_sn_mms_anisotropic_symbolic.py` (one
  ``@pytest.mark.foundation`` test per ``derive_*`` function).
- **Branch 2 (vectorised numpy)**:
  :class:`orpheus.derivations.continuous.mms.sn.SNSphericalAnisotropicMMSCase`
  and
  :class:`orpheus.derivations.continuous.mms.sn.SNCylindricalAnisotropicMMSCase`
  evaluate :eq:`sn-mms-spherical-aniso-qext` and
  :eq:`sn-mms-cylindrical-aniso-qext` per ordinate using vectorised
  numpy.
- **L1 cross-check (the gate)**: the Branch-2 numerical
  :math:`Q^{\text{ext}}_n(r_i)` agrees with Branch-1 SymPy-evaluated
  :math:`Q^{\text{ext}}_n(r_i)` (via :func:`sympy.lambdify`) to
  :math:`\sim 10^{-16}` (max absolute) on a sample mesh in both
  geometries. Tested in
  :func:`tests.derivations.test_sn_mms_anisotropic_symbolic.test_spherical_aniso_numerical_qext_matches_sympy`
  and the cylindrical sibling.

**Code pointers**

- Derivations: :class:`orpheus.derivations.continuous.mms.sn.SNSphericalAnisotropicMMSCase`,
  :class:`orpheus.derivations.continuous.mms.sn.SNCylindricalAnisotropicMMSCase`,
  :func:`orpheus.derivations.continuous.mms.sn.build_spherical_anisotropic_mms_case`,
  :func:`orpheus.derivations.continuous.mms.sn.build_cylindrical_anisotropic_mms_case`.
- Symbolic factory:
  :func:`orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`,
  :func:`orpheus.derivations.continuous.mms.sn._cylindrical_anisotropic_symbolic`.
- Foundation tests:
  :file:`tests/derivations/test_sn_mms_anisotropic_symbolic.py`.
- Consumer L1 convergence test (Phase-0 work, separate branch):
  ``tests/sn/_l1/test_mms_spherical_anisotropic_dd_convergence_O_h2.py``
  (planned).


.. _sn-case-heterogeneous-verification:

Heterogeneous eigenvalue — Case singular-eigenfunction method
--------------------------------------------------------------

Phase 2.1b of the verification campaign closes the last
heterogeneous gap in the SN verification ladder: the
**eigenvalue iteration** on a 1-group two-region reflective
slab, verified against a semi-analytical reference derived
from the discrete-:math:`S_N` slope matrix itself --- no
diffusion approximation, no cross-code comparison, no
Richardson self-test.

The reference is produced by
:func:`orpheus.derivations.continuous.cases.sn.derive_sn_heterogeneous_continuous`
and consumed by
:func:`tests.sn.test_heterogeneous_transport.test_sn_2region_reflective_case_eigenvalue`
(eigenvalue) and
:func:`tests.sn.test_heterogeneous_transport.test_sn_2region_reflective_flux_shape`
(scalar flux shape). The Phase 2.1a smooth-:math:`\Sigma` MMS
test verifies the **spatial operator** at :math:`\mathcal O(h^{2})`
design order; this section's Case method verifies the
**eigenvalue** iteration at the material-interface-degraded
:math:`\mathcal O(h)` rate expected for diamond-difference on
piecewise-constant :math:`\Sigma`.

**Motivation: why a second verification path.** The Phase 2.1a
MMS test deliberately uses smooth :math:`\Sigma(x)` to avoid
interface degradation and hit the :math:`\mathcal O(h^{2})`
design order of diamond difference. That is the right choice
for verifying the spatial operator, but it **cannot** exercise
the heterogeneous-interface regime where material
discontinuities force the operator into its interface-layer
behaviour --- the regime where a significant fraction of
production solver bugs live (including ERR-025; see
:ref:`investigation-err-025`). The Case singular-eigenfunction
method provides the complementary reference: an eigenvalue
solution with genuine material-interface discontinuities, built
from the transport equation without running the solver.

**Operator.** The 1-group 1D slab SN transport equation in a
single region with cross sections
:math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f)` and reflective BCs
is, per ordinate,

.. math::
   :label: sn-case-per-ordinate

   \mu_n\,\frac{d\psi_n}{dx} + \Sigma_t\,\psi_n
     \;=\; \frac{c_\text{eff}(k)}{W}\,\phi,
   \qquad
   \phi = \sum_m w_m\,\psi_m,
   \qquad
   c_\text{eff}(k) = \Sigma_s + \frac{\nu\Sigma_f}{k},

where :math:`W = \sum_m w_m`. Substituting the scalar-flux
definition and stacking the angular flux into
:math:`\mathbf y \in \mathbb R^N` (for Gauss--Legendre order
:math:`N`), the system becomes a first-order constant-coefficient
ODE

.. math::
   :label: sn-case-slope-matrix

   \frac{d\mathbf y}{dx} \;=\; \mathbf S(k)\,\mathbf y,
   \qquad
   \mathbf S(k)[n, m] \;=\; \frac{1}{\mu_n}
       \left(-\Sigma_t\,\delta_{nm}
             + \frac{c_\text{eff}(k)}{W}\,w_m\right).

Note the **row-scaling** :math:`1/\mu_n`: the slope matrix is
generally non-symmetric even for symmetric GL quadrature,
because the angular ODE has different "speeds" for different
ordinates.

**Per-region spatial modes.** For each region (fuel at
:math:`x \in [0, H_A]` and moderator at :math:`x \in [H_A, L]`),
diagonalise :math:`\mathbf S(k)`:

.. math::
   :label: sn-case-spatial-modes

   \mathbf S(k)\,\mathbf v_i \;=\; \lambda_i\,\mathbf v_i,
   \qquad i = 1,\ldots,N,

via :func:`numpy.linalg.eig`. For subcritical regions
(:math:`c_\text{eff}(k) < 1`, typical moderator) the eigenvalues
come in :math:`\pm` real pairs. For supercritical regions
(:math:`c_\text{eff}(k) > 1`, fuel at :math:`k` below
:math:`k_{\infty,\text{fuel}}`) some pairs are
complex-conjugate. Each real eigenvalue gives one exponential
mode :math:`\exp(\lambda\,x)\,\mathbf v`; each complex-conjugate
pair gives two real modes built from the canonical
:math:`\cos/\sin/\Re/\Im` combination.

**Real bounded basis.** The naive unbounded basis
:math:`\exp(\lambda\,x)\,\mathbf v` is catastrophically
ill-conditioned for optically thick slabs --- the Phase 1.2
diffusion investigation history records the ``expm``-based
transfer-matrix composition dying from :math:`\text{cond}
\sim 10^{17}` on an 80-cm slab, finding spurious roots with
:math:`\mathcal O(10^{-3})` null-vector residuals rather than
machine-precision zeros. The fix, ported verbatim to Phase 2.1b,
is to **anchor each mode at the nearer region edge**:

.. math::
   :label: sn-case-real-basis

   m^{\text{real}}_j(x) &\;=\; \exp(\lambda_j\,\xi_j)\,\mathbf v_j,
       \qquad
       \xi_j = \begin{cases}
         x - L_\text{reg} & \lambda_j \ge 0 \;\;\text{(anchor right)} \\
         x                & \lambda_j < 0 \;\;\text{(anchor left)}
       \end{cases} \\[1mm]
   m^{\text{c}}_j(x) &\;=\; e^{\Re\lambda_j\,\xi_j}\,
       \bigl(\cos(\Im\lambda_j\,\xi_j)\,\mathbf v_{R,j}
          - \sin(\Im\lambda_j\,\xi_j)\,\mathbf v_{I,j}\bigr), \\
   m^{\text{s}}_j(x) &\;=\; e^{\Re\lambda_j\,\xi_j}\,
       \bigl(\sin(\Im\lambda_j\,\xi_j)\,\mathbf v_{R,j}
          + \cos(\Im\lambda_j\,\xi_j)\,\mathbf v_{I,j}\bigr),

where :math:`\mathbf v_j = \mathbf v_{R,j} + i\,\mathbf v_{I,j}`
is the complex eigenvector. Every mode is bounded by
:math:`|\mathbf v_j|` on its region, so the assembled matching
matrix has :math:`\mathcal O(1)` entries.

**Matching matrix.** For the 2-region reflective slab the
coefficient vector has dimension :math:`2N` (one real mode per
eigenvalue per region). The linear constraints are:

.. math::
   :label: sn-case-matching-matrix

   &\text{Reflective at } x = 0:\quad
      \psi^A_n(0) - \psi^A_{N-1-n}(0) = 0,
      \qquad n \in [0, N/2) \\[1mm]
   &\text{Interface at } x = H_A:\quad
      \psi^A_n(H_A) - \psi^B_n(H_A) = 0,
      \qquad n \in [0, N) \\[1mm]
   &\text{Reflective at } x = L:\quad
      \psi^B_n(L) - \psi^B_{N-1-n}(L) = 0,
      \qquad n \in [0, N/2)

:math:`N/2 + N + N/2 = 2N` equations in :math:`2N` unknowns.
The partner index :math:`N-1-n` is the Gauss--Legendre
reflection pairing (ordinates sorted by ascending :math:`\mu`).
The eigenvalue condition is
:math:`\det\mathbf C(k) = 0`.

**Root finding.** :func:`scipy.optimize.brentq` on
:math:`\det\mathbf C(k)` over a coarse :math:`k`-scan, with
sign-change bracketing, refines every candidate to
``xtol=1e-14``. But :func:`numpy.linalg.eig`'s eigenvalue
ordering is not a continuous function of :math:`k` --- at
parameter values where two per-region eigenvalues cross, the
eigenvalue labels permute discontinuously, and
:math:`\det\mathbf C(k)` flips sign by permutation rather than
by passing through zero. brentq will "converge" to such
spurious points.

**Physical validation.** Every candidate root is rebuilt via
SVD of :math:`\mathbf C(k)`, and the null vector's reflective-BC
residuals at :math:`x = 0` and :math:`x = L`, and the interface
continuity residual at :math:`x = H_A`, are explicitly
reconstructed and checked against a dimensionless tolerance
relative to the peak angular flux:

.. math::
   :label: sn-case-physical-validation

   \|\psi(0, +\mu_n) - \psi(0, -\mu_n)\| / \|\psi\|_\text{peak}
     &< \text{tol} \\
   \|\psi^A(H_A) - \psi^B(H_A)\| / \|\psi\|_\text{peak}
     &< \text{tol} \\
   \|\psi(L, +\mu_n) - \psi(L, -\mu_n)\| / \|\psi\|_\text{peak}
     &< \text{tol}

Only candidates passing all three are accepted; the fundamental
is the largest validated root. This is the SN analogue of the
Phase 1.2 diffusion physical validation (same pattern, different
operator).

**Back-substitution.** Once :math:`k_\text{fund}` is found,
the null vector at that :math:`k` is the coefficient vector in
the :math:`2N`-dimensional real basis. Evaluation of
:math:`\phi(x) = \sum_n w_n\,\psi_n(x)` at any point reduces to
a linear combination of a handful of bounded exponential or
trigonometric modes:

.. math::
   :label: sn-case-back-substitution

   \psi(x) = \begin{cases}
     \sum_j c^A_j\,m^A_j(x) & x \le H_A \\[1mm]
     \sum_j c^B_j\,m^B_j(x - H_A) & x > H_A
   \end{cases},
   \qquad
   \phi(x) = \sum_n w_n\,\psi_n(x).

All modes are bounded by :math:`\mathcal O(1)`, so
:math:`\phi(x)` is stable to machine precision.

**The Phase 2.1b diagnostic configuration.** The canonical
test problem is the ``A`` + ``B`` 1-group mixture pair from
:mod:`orpheus.derivations.common.xs_library`:

.. list-table::
   :header-rows: 1
   :widths: 15 15 15 15 15

   * - Region
     - :math:`\Sigma_t`
     - :math:`\Sigma_s`
     - :math:`\nu\Sigma_f`
     - :math:`k_\infty`
   * - A (fuel)
     - 1.0
     - 0.5
     - 0.75
     - 1.5
   * - B (moderator)
     - 2.0
     - 1.9
     - 0
     - ---

with :math:`H_A = H_B = 0.5\,\text{cm}`, reflective BCs on both
outer edges, :math:`S_8` Gauss--Legendre quadrature. The
resulting Case reference is

.. math::

   k_\text{eff}^{\text{Case}}(S_8) = 1.2746160417

--- the exact discrete-:math:`S_8` eigenvalue. For
cross-validation, the same configuration run through ORPHEUS's
:func:`~orpheus.cp.solver.solve_cp` (1D slab E\ :sub:`3` kernel,
completely independent numerical path) gives
:math:`k^{\text{CP}} = 1.2744284665` --- agreement to
:math:`\sim 2\times 10^{-4}`, well below the :math:`\mathcal O(1\%)`
difference that typically exists between discrete-SN and
continuous-angle formulations. This cross-check is used only as
a sanity input, not as a verification crutch.

**Measured convergence.** With :math:`S_8`, refining
:math:`n_\text{per}` per region:

.. list-table::
   :header-rows: 1
   :widths: 15 25 15

   * - :math:`n_\text{per}`
     - :math:`k_\text{solve}`
     - :math:`|k_\text{solve} - k_\text{Case}|`
   * - 20
     - 1.2746074093
     - :math:`\sim 8.6\!\times\!10^{-6}`
   * - 40
     - 1.2746138837
     - :math:`\sim 2.2\!\times\!10^{-6}`
   * - 80
     - 1.2746155022
     - :math:`\sim 5.4\!\times\!10^{-7}`
   * - 160
     - 1.2746159068
     - :math:`\sim 1.3\!\times\!10^{-7}`
   * - 320
     - 1.2746160080
     - :math:`\sim 3.4\!\times\!10^{-8}`

Each refinement roughly halves the error, confirming the
:math:`\mathcal O(h)` rate expected at a material interface with
piecewise-constant :math:`\Sigma`. The finest-mesh residual of
:math:`3.4 \times 10^{-8}` is **machine-precision agreement**
between two independent mathematical constructions (the Case
matching-matrix + back-substitution reference and the
diamond-difference sweep-based power iteration); both
implementations solve the same discrete-:math:`S_N` spectral
problem and agree to within the BiCGSTAB-compatible
truncation.

**Contrast with Phase 2.1a.** The Phase 2.1a MMS section hits
:math:`\mathcal O(h^{2})` because it uses smooth
:math:`\Sigma(x)`; the Phase 2.1b Case section hits
:math:`\mathcal O(h)` because it uses piecewise-constant
:math:`\Sigma(x)` with a genuine material interface. Both are
correct for their respective regimes. The degradation from
:math:`h^{2}` to :math:`h` at the interface is the standard
Salari--Knupp result for DD on discontinuous coefficients, and
is the **reason** Phase 2.1a deliberately chose smooth
:math:`\Sigma` to isolate the spatial operator.


Homogeneous Infinite Medium
----------------------------

For homogeneous geometry with reflective BCs, the flux is spatially flat
and :math:`\keff = \lambda_{\max}(A^{-1}F)`.  This is geometry-independent
--- Cartesian, spherical, and cylindrical must all give the same
:math:`\keff`.

.. list-table::
   :header-rows: 1
   :widths: 10 14 19 19 19 19

   * - Groups
     - :math:`\kinf`
     - Cartesian (GL S8)
     - Spherical (GL S8)
     - Cylindrical (Prod 4x8)
     - Cylindrical (LS S4)
   * - 1
     - 1.5000
     - exact
     - exact
     - exact
     - exact
   * - 2
     - 1.8750
     - exact
     - exact
     - exact
     - exact
   * - 4
     - 1.4878
     - exact
     - exact
     - exact
     - exact

All entries are exact to machine precision.  Spherical 2G/4G results
(previously showing ~1% error) are now exact thanks to the M-M angular
closure weights.

Heterogeneous Convergence
--------------------------

For a cylindrical fuel (r < 0.5) + moderator (r < 1.0) geometry with
Product(4x8) quadrature:

.. list-table::
   :header-rows: 1
   :widths: 20 25 25

   * - Cells/region
     - :math:`\keff` (1G)
     - :math:`\Delta k` from previous
   * - 5
     - 0.9769
     -
   * - 10
     - 0.9842
     - +0.0073
   * - 20
     - 0.9874
     - +0.0032

:math:`\keff` converges monotonically toward the CP reference
(0.9955).  The ~1% residual gap is the white-BC (CP) vs reflective-BC
(SN) approximation difference, consistent with the slab geometry
findings.

For the 2G heterogeneous resolution test, Product(4x8) and Product(8x8)
agree to :math:`< 0.01\%` (keff = 0.7227 for both), confirming
angular convergence.

Why 1-Group Verification Is Degenerate
----------------------------------------

For 1 energy group, the eigenvalue is:

.. math::

   k = \frac{\nSigf{}}{\Sigma_a}

This is a scalar ratio independent of the spatial or angular flux
distribution.  Consequences:

- Angular weight errors scale all flux equally --- cancels in :math:`k`.
- Wrong scattering convention --- no inter-group coupling to distort.
- Wrong flux shape --- does not matter; :math:`k` is a material property.

Only multi-group problems have a flux-shape-dependent eigenvalue:
:math:`k = (\nSigf{} \cdot \phi) / (\Sigma_a \cdot \phi)` where the
dot product weights each group differently.  A wrong group ratio (from
angular errors, normalization errors, or convergence failures) directly
shifts :math:`\keff`.

**Rule:** Every transport solver must be verified on at least 2-group
problems.  1-group success gives false confidence.

Spatial and Angular Convergence
--------------------------------

The diamond-difference scheme converges at :math:`O(h^2)` with mesh
refinement.  Gauss--Legendre quadrature shows spectral convergence in
angle.  Both are verified in ``test_sn_1d.py``.

Property Tests
---------------

For all geometries:

- **Particle balance**: production / absorption :math:`= \keff`
- **Flux non-negativity**: :math:`\phi \geq 0` everywhere
- **Angular flux at** :math:`r = 0` **all positive** (curvilinear only)
- **Multi-group eigenvector not flat**: flux spectrum differs between
  fuel and moderator (catches 1G-degenerate bugs)

Run the full suite::

   pytest tests/sn/ -v -m "not slow"


.. _investigation-history:

Investigation History: Curvilinear Bug
=======================================

This section documents the full history of the cylindrical DD bug and its
resolution.  It is preserved to prevent future sessions from repeating
the same dead ends.

Symptoms
--------

1. Homogeneous eigenvalue problems: exact (1G/2G/4G).
2. Heterogeneous eigenvalue problems: divergent :math:`\keff` with mesh
   refinement (5 cells: 1.15, 10 cells: 0.90, 20 cells: 0.52).
3. Fixed-source flux range: ``[0.59, 5.09]`` (should be near-flat).
4. :math:`\keff` depended strongly on angular quadrature order
   (4x8 vs 8x8 gave a 67% gap on heterogeneous problems).

The Root Cause
--------------

**Two bugs**, both breaking per-ordinate flat-flux consistency:

**Bug 1: Wrong** :math:`\alpha` **recursion.**  The code used
:math:`\alpha = \text{cumsum}(+w\xi)` with the azimuthal cosine
:math:`\xi` (``mu_y``).  The correct recursion ([Bailey2009]_ Eq. 50)
is :math:`\alpha = \text{cumsum}(-w\eta)` with the radial cosine
:math:`\eta` (``mu_x``), and ordinates must be sorted by increasing
:math:`\eta` within each level.

**Bug 2: Missing** :math:`\Delta A/w` **geometry factor.**  The
redistribution term in the balance equation must include
:math:`\Delta A_i / w_m`.  Without this factor, the streaming and
redistribution do NOT cancel per-ordinate for a spatially flat flux,
creating artificial angular anisotropy that worsens with mesh refinement
near :math:`r = 0`.

Six Failed Approaches
----------------------

Before the correct fix was found, six approaches were tested.  All
failed because they addressed symptoms, not the root cause:

1. **Reverse sweep:** Reversed the azimuthal sweep direction.
   No effect --- the dome shape is symmetric.

2. **Step closure:** Replaced DD (:math:`\tau = 0.5`) with step
   (:math:`\tau = 1.0`) in angle.  Reduced the divergence slightly
   but did not eliminate it.

3. **Lewis & Miller starting direction:** Tracked :math:`\alpha\psi`
   instead of :math:`\psi` alone (Section 4.5.4 of [LewisMiller1984]_).
   Unnecessary when the :math:`\alpha` recursion is correct, since
   :math:`\alpha_{1/2} = 0`.

4. **Bidirectional sweep:** Swept azimuthal ordinates in both
   directions and averaged.  Masked the asymmetry but did not fix it.

5. **Scaled** :math:`\alpha`:  Empirically scaled the :math:`\alpha`
   coefficients.  Could not find a consistent scaling.

6. **Zero redistribution:** Set :math:`\alpha = 0` to isolate the
   spatial streaming.  Confirmed the bug was in the redistribution,
   not the spatial DD.

Why the Original Sign Convention Hypothesis Was Wrong
------------------------------------------------------

The initial hypothesis was that the minus sign before the redistribution
in :eq:`transport-cylindrical` required special treatment (tracking
:math:`\alpha\psi`, using :math:`|\alpha|`, etc.).  This was wrong.
The sign convention is absorbed into the :math:`\alpha` definition:
:math:`\text{cumsum}(-\eta w)` with a :math:`+` sign in the balance
equation gives the same physics as :math:`\text{cumsum}(+\xi w)` with
a :math:`-` sign, for symmetric quadratures.

The real issue was the missing :math:`\Delta A/w` geometry factor, which
has nothing to do with signs.

The Fix
-------

Applied the [Bailey2009]_ formulation:

1. Corrected :math:`\alpha` recursion:
   :math:`\alpha_{m+1/2} = \alpha_{m-1/2} - w_m \eta_m`
   with ordinates sorted by increasing :math:`\eta`.

2. Added :math:`\Delta A/w` geometry factor to both the DD sweep and
   the BiCGSTAB operator.

3. Added Morel--Montry angular closure weights (M-M) to eliminate the
   flux dip at :math:`r = 0`.

4. Applied the same :math:`\Delta A/w` fix to the spherical sweep
   (which had the same missing factor).  The fixed-source flux spike
   at :math:`r = 0` dropped from 5.1x to 1.1x.

5. Applied the fix to both the spherical and cylindrical **BiCGSTAB
   operators** (:func:`transport_operator_matvec_spherical`,
   :func:`transport_operator_matvec_cylindrical`).  Multi-group
   spherical BiCGSTAB had been unstable (keff → NaN for 2G+); the
   root cause was the same missing :math:`\Delta A/w` factor in the
   explicit FD operator.  After the fix, 2G and 4G spherical BiCGSTAB
   converge to :math:`< 10^{-6}` of the analytical eigenvalue.

Results After Fix
------------------

.. list-table::
   :header-rows: 1
   :widths: 35 25 25

   * - Test
     - Before
     - After
   * - Homogeneous 1G/2G/4G
     - Exact
     - Exact
   * - Heterogeneous 1G, 5/10/20 cells
     - 1.15 / 0.90 / 0.52 (diverges)
     - 0.977 / 0.984 / 0.987 (converges)
   * - Heterogeneous 2G, 4x8 vs 8x8
     - 0.54 vs 0.91 (67% gap)
     - 0.723 vs 0.723 (<0.01%)
   * - Fixed-source flux range (40 cells)
     - [0.59, 5.09] (spike)
     - [0.51, 1.12] (bounded)
   * - Contamination :math:`\beta` (cylindrical)
     - ~2.0
     - ~10\ :sup:`-16` (machine zero)


.. _investigation-err-025:

Investigation History: ERR-025 and the Phase 2.1b Case reference
=================================================================

This section documents the full diagnostic history of ERR-025 and
the Phase 2.1b Case singular-eigenfunction reference, preserved to
give future sessions a vivid example of how silent
derivation-implementation drift can survive every eigenvalue-based
test a module has.

Symptoms
--------

Phase 2.1b originally set out to land a Case singular-eigenfunction
eigenvalue reference for the 1-group two-region reflective slab
(configuration in :ref:`sn-case-heterogeneous-verification`). The
prototype Case code was self-consistent: homogeneous limits exact,
region-swap invariance exact, :math:`H_B \to 0` limit matching
:math:`k_{\infty,\text{fuel}}`, physical-reconstruction residuals at
machine precision. But the mesh-refined ``solve_sn`` sequence on the
same configuration converged cleanly to a **different** value:

.. list-table::
   :header-rows: 1
   :widths: 20 30 25

   * - Implementation
     - :math:`k_\text{eff}`
     - Residual vs Case
   * - Case singular-eigenfunction
     - 1.27461604
     - ---
   * - CP slab :math:`E_3` kernel (converged)
     - 1.27442847
     - :math:`-1.88\times 10^{-4}`
   * - ``solve_sn`` S\ :sub:`8`, :math:`n=1280`
     - 1.25986980
     - :math:`-1.48\times 10^{-2}`

Two independent transport methods (Case + CP) agreed to
:math:`\sim 2\times 10^{-4}`; ``solve_sn`` disagreed by almost two
orders of magnitude more. The gap did not shrink with mesh
refinement or with increasing quadrature order (:math:`S_4` through
:math:`S_{64}` all converged to the wrong asymptote). It was also
**invisible on every homogeneous test**: pure fuel slab gave
:math:`k_\infty = 1.5` exactly, same-material two-region gave
:math:`k_\infty` exactly, and the Phase 2.1a smooth-:math:`\Sigma`
heterogeneous MMS test hit :math:`\mathcal O(h^{2})`. The only
configuration that showed the discrepancy was the one with a
piecewise-constant material interface under the eigenvalue driver.

Initial Hypothesis (wrong)
--------------------------

The first hypothesis was that the Case implementation was computing
the **continuous-angle** limit of the SN eigenvalue while
``solve_sn`` was solving the discrete-:math:`S_8` problem, and the
:math:`\sim 1.5\%` gap was the quadrature-order error expected for
:math:`S_8` Gauss--Legendre (:math:`\sim 1/N^{2} \approx 1.5\%`).

Testing this hypothesis with an :math:`S_N` order sweep killed it
immediately: ``solve_sn`` at :math:`S_4, S_8, S_{16}, S_{32}, S_{64}`
converged cleanly **in angle** to :math:`\sim 1.2609`, not to
:math:`1.2746` --- plateau of :math:`\sim 3\times 10^{-5}` at
:math:`S_{32}`. The two methods were not approaching the same
continuous-angle answer. The disagreement was structural.

Cross-Check That Localised The Bug
-----------------------------------

The next move was to run the same configuration through a
completely independent transport method: ORPHEUS's collision
probability solver :func:`orpheus.cp.solver.solve_cp` on a 1D slab
with the :math:`E_3` kernel. CP uses the Peierls integral equation
--- no diamond difference, no finite-volume discretisation, no
explicit ordinate sweep --- and its numerical path has nothing in
common with ``solve_sn``. CP converged to :math:`k = 1.27442847`,
essentially matching the Case reference to the CP quadrature
precision.

Result: two independent methods agreed on :math:`k \approx 1.27461`,
``solve_sn`` disagreed by :math:`\sim 1.5\%`. **The SN solver was
the outlier.** Since homogeneous and same-material problems worked,
the bug had to live somewhere that activates specifically at a
**material interface** in the ``solve_sn`` code path.

Root Cause: ERR-025
-------------------

A focused audit of the four SN sweep paths and the source-builder
helpers localised the bug to
:func:`orpheus.sn.sweep._sweep_1d_cumprod`, the 1D Cartesian
Gauss--Legendre fast path. Its diamond-difference face-flux
recurrence coefficients were

.. math::

   a_\text{bug} &= \frac{2\mu}{2\mu + \Delta x\,\Sigma_t}
     \qquad (\text{WRONG: missing } -\Sigma_t\text{ in numerator}) \\
   b_\text{bug} &= \frac{0.5\,\Delta x\,Q}{2\mu + \Delta x\,\Sigma_t}
     \qquad (\text{WRONG: missing } 1/W\text{, extra factor } 0.5)

instead of the canonical diamond-difference recurrence derived
symbolically in
:func:`orpheus.derivations.discrete.sn.balance.derive_cumprod_recurrence`:

.. math::

   a &= \frac{2\mu - \Delta x\,\Sigma_t}{2\mu + \Delta x\,\Sigma_t} \\
   b &= \frac{2\,\Delta x\,(Q/W)}{2\mu + \Delta x\,\Sigma_t}

where :math:`W = \sum_n w_n` is the quadrature weight sum, needed
because :func:`orpheus.sn.solver.SNSolver._add_scattering_source`
produces :math:`Q` in **scalar-flux units** while the per-ordinate
transport equation sees :math:`Q/W` as its right-hand side. The 2D
wavefront sweep :func:`~orpheus.sn.sweep._sweep_2d_wavefront`
already applied this normalisation via its ``weight_norm = 1/W``
factor; the 1D fast path had been independently derived without it
and drifted silently.

Why the Two Errors Cancel For Homogeneous Problems
---------------------------------------------------

The fixed point of the buggy recurrence is
:math:`\psi_n = Q/(2\Sigma_t)`, half the correct
:math:`Q/\Sigma_t`. But for Gauss--Legendre on :math:`[-1, 1]`,
:math:`W = \sum_n w_n = 2`, so the missing :math:`1/W = 1/2`
multiplies the buggy fixed point by exactly :math:`2`, turning
:math:`Q/(2\Sigma_t)` back into the correct
:math:`\psi_n = Q/(W\Sigma_t)` per ordinate. The resulting scalar
flux

.. math::

   \phi = \sum_n w_n\,\psi_n = W\cdot\frac{Q}{W\Sigma_t} = \frac{Q}{\Sigma_t}

is identical in magnitude to the correct value. This is why every
homogeneous test passed at machine precision, including
:math:`k_\infty` assertions to :math:`10^{-8}`.

For an eigenvalue problem, even without the :math:`W=2` coincidence
rescaling, the Rayleigh quotient

.. math::

   k \;=\; \frac{\nu\Sigma_f\,\phi}{\Sigma_a\,\phi}

is **invariant** under a uniform rescaling :math:`\phi \to C\phi`,
because :math:`C` cancels between numerator and denominator.
Homogeneous and same-material-multi-region problems have a
uniform-in-:math:`x` rescaling, so the buggy :math:`k_\text{eff}`
is exact. Only at a **material interface** does the rescale factor
:math:`C(x)` become :math:`x`-dependent (through :math:`\Sigma_t(x)`),
and only then does the cancellation break and a real error appear.

Dead End #1: "It must be S\ :sub:`8` vs continuous angle"
----------------------------------------------------------

The :math:`\sim 1.5\%` magnitude is numerically close to the
typical :math:`1/N^{2}` error for Gauss--Legendre :math:`S_8`,
which made this hypothesis seductive. Cost to refute: 30 seconds
with an :math:`S_N` sweep showing the gap was invariant in
quadrature order. Lesson: **always run the cheapest diagnostic
first**. A 30-second experiment would have saved an hour of
grooming the Case code for non-existent quadrature-convention
bugs.

Dead End #2: "It must be the Case code --- its symmetry checks
passed but it's newer code"
--------------------------------------------------------------

Before the CP cross-check, the natural first suspicion was the
Case prototype. It was ephemeral session code with hand-derived
algebra and no reference implementation to compare against; the
``solve_sn`` path was production code shipping for months with a
full test suite. The CP cross-check inverted this: two
**independent mathematical constructions** (Case via
eigendecomposition, CP via :math:`E_3` integral kernel) agreed,
and the production solver was the outlier. Lesson: **trust
agreement over pedigree**. A test suite that never exercised the
failure mode is not evidence of correctness, no matter how big
it is.

Dead End #3: "Maybe the reflective BC is subtly wrong"
-------------------------------------------------------

One hypothesis was that ``_sweep_1d_cumprod``'s reflective BC
persistence (via the ``psi_bc["bc_1d"]`` dict between outer
iterations) was mishandling the ordinate pairing at a material
interface. The bug is upstream of the BC code --- the coefficients
are wrong before the BC ever touches them --- but the hypothesis
was plausible enough to warrant an audit of the BC application
code. It was clean. Lesson: **trace the data flow in order**.
The BC is applied to the output of the recurrence; if the
recurrence itself is wrong, the BC can't fix it and can't
manifest the bug.

Dead End #4: "Maybe it's the ``compute_keff`` volume integration"
-----------------------------------------------------------------

Another hypothesis was that :func:`~orpheus.sn.solver.SNSolver.compute_keff`
was accumulating :math:`\nu\Sigma_f\,\phi\,V` with a wrong
material-id lookup at the interface cell (e.g., using
``mat_ids[c+1]`` for cell ``c``'s fission). The code turned out to
use the correct per-cell material id. Lesson: **read the code
before blaming it**. The audit took five minutes and ruled out
two more hypotheses for free (``_add_fission_source`` and
``_add_scattering_source`` both checked).

The Fix and the Test
---------------------

The fix is a one-formula correction in
:func:`orpheus.sn.sweep._sweep_1d_cumprod` plus a comment block
pointing at ``sn_balance.derive_cumprod_recurrence`` as the
source of truth. After the fix, Case :math:`\leftrightarrow`
``solve_sn`` agreement went from :math:`1.48\times 10^{-2}` to
:math:`3.4\times 10^{-8}` --- a six-order-of-magnitude
improvement at matching quadrature order. All 165 SN tests pass,
and two new tests landed with the fix:

1. **L0 term verification**
   :func:`tests.sn.test_cartesian.test_sweep_1d_cumprod_recurrence_matches_symbolic_derivation`
   --- a white-box test that calls ``_sweep_1d_cumprod`` on a
   1-cell homogeneous slab with a controlled inflow and source,
   and checks the returned cell-average angular flux against a
   numerical substitution of the **symbolic** expressions
   produced by ``derive_cumprod_recurrence()``. This is the
   minimal isolation of the failure mode; it runs in
   milliseconds and does not need any reference solver.
2. **L1 regression**
   :func:`tests.sn.test_cartesian.test_heterogeneous_absolute_keff`
   --- a black-box test that pins the 2-region A+B reflective
   slab against the Case reference to :math:`5\times 10^{-4}`.
   The pre-fix solver (:math:`1.48\times 10^{-2}` error) would
   fail this assertion by almost two orders of magnitude.

Both tests were verified to **fail** on the buggy code before
the fix was landed.

Aftermath: Issue #95 and the Broader Audit
-------------------------------------------

ERR-025 was a silent derivation-implementation drift: a symbolic
derivation existed (``sn_balance.derive_cumprod_recurrence``) and
was **correct**, but the implementation had independently re-derived
the coefficients and gotten them wrong, with no mechanical link
between the two. GitHub issue #95 tracks the follow-up audit work
to systematically check every solver implementation against its
symbolic derivation module. During the ERR-025 audit, two
pre-existing BiCGSTAB inconsistencies were surfaced (Cartesian
BiCGSTAB using upwind cell-centre FD instead of DD, curvilinear
BiCGSTAB using arithmetic face averages instead of the sweep's DD
closure); these are tracked as issues #96 and #97 respectively.

The four sweep paths audited during ERR-025 diagnosis ---
:func:`~orpheus.sn.sweep._sweep_2d_wavefront`,
:func:`~orpheus.sn.sweep._sweep_1d_spherical`,
:func:`~orpheus.sn.sweep._sweep_1d_cylindrical`, and
post-fix ``_sweep_1d_cumprod`` --- were all verified **clean**
against the ``sn_balance`` symbolic derivation. The source-builder
helpers (:meth:`~orpheus.sn.solver.SNSolver._add_scattering_source`,
:meth:`~orpheus.sn.solver.SNSolver._add_fission_source`,
:meth:`~orpheus.sn.solver.SNSolver._build_aniso_scattering`) were
also audited clean.

Meta-Lessons
------------

1. **Uniform-rescale invariance hides coefficient bugs.** Any
   eigenvalue problem where :math:`\phi` is the target quantity
   is invariant under :math:`\phi \to C\phi`, which makes it
   blind to factor-of-two errors that preserve the flux shape.
   Homogeneous and same-material-multi-region problems have
   spatially uniform rescaling; only genuine material
   interfaces break the cancellation. **Always include at least
   one absolute-:math:`\phi` test** (fixed-source, or an absolute
   eigenvalue comparison against an independent reference) to
   expose rescale-invariant bugs.

2. **Symbolic derivations must be load-bearing, not decorative.**
   ``sn_balance.derive_cumprod_recurrence`` existed, was correct,
   and was not referenced from anywhere in the consuming code.
   It became a museum piece. A comment in the implementation
   pointing at the derivation function would have caught this at
   code review. Issue #95 proposes a CI check to flag orphan
   derivations across the whole codebase.

3. **Cross-check before debugging.** When a self-consistent
   mathematical construction disagrees with a production solver
   and you do not know which is right, spend the 30-minute
   budget to run an **independent third** implementation
   before going deep on either. CP vs SN vs Case is the pattern
   that made the ERR-025 diagnosis possible in under an hour.
   This is explicitly **not** a verification crutch --- the
   final Phase 2.1b reference stands on its own mathematical
   merits --- but as an investigation sanity input it is
   invaluable.

4. **Cheap diagnostic first.** The :math:`S_N` order sweep that
   killed Dead End #1 took 30 seconds; the CP cross-check took
   five minutes. Both ran before any serious debugging. Lesson
   from Phase 1.2 investigation (diffusion): the order of
   operations matters as much as the operations themselves.


Numerical Sensitivities
========================

:math:`\keff` Sensitivity Table (421-Group Heterogeneous PWR Slab)
-------------------------------------------------------------------

All cases: 10 cells, :math:`\delta = 0.2` cm, material layout
``[fuel x 5, clad x 1, cool x 4]``, P0 scattering, 421 energy groups.

.. list-table::
   :header-rows: 1
   :widths: 50 15 35

   * - Configuration
     - :math:`\keff`
     - Notes
   * - 1D GL S16, BiCGSTAB (FD operator)
     - 1.03882
     - True 1D, 16 ordinates
   * - 1D Lebedev 110, source iteration (DD sweep)
     - 1.04294
     - 1D mesh, 2D quadrature
   * - 2D (10x2) Lebedev 110, source iter (DD sweep)
     - 1.04294
     - Pseudo-2D, full volumes
   * - 2D (10x2) Lebedev 110, BiCGSTAB (FD)
     - 1.04007
     - Pseudo-2D, full volumes
   * - 2D (10x2) Lebedev 110, BiCGSTAB, half-volumes
     - 1.04192
     - MATLAB convention
   * - **MATLAB reference**
     - **1.04188**
     - 2D Lebedev, FD, half-volumes

Sources of Variation
---------------------

1. **Angular quadrature** (GL vs Lebedev): ~0.004 difference.
   GL S16 integrates 1D angular flux with 16 points on :math:`[-1,1]`.
   Lebedev 110 integrates over the unit sphere --- more angular
   resolution but different effective weights per :math:`\mu_x`
   direction.  On a coarse heterogeneous mesh, these give different
   eigenvalues.

2. **Spatial discretisation** (DD sweep vs FD gradient): ~0.003
   difference.  Source iteration uses the DD wavefront sweep
   (:math:`T^{-1}`).  BiCGSTAB uses the explicit FD transport operator
   (:math:`T`).  Both are :math:`O(h)` on this mesh but with different
   truncation error constants.

3. **Boundary volume weighting**: ~0.002 difference (full vs half).
   The MATLAB code halves boundary cell volumes.  With ``ny=2`` and
   materials uniform in *y*, only the *x*-direction halving (fuel edge,
   coolant edge) affects :math:`\keff`.  This is an artifact of the
   pseudo-2D implementation: a true 1D calculation has no *y*-volumes.

4. **Inner convergence**: source iteration with ``max_inner=200``,
   ``inner_tol=1e-8`` does not fully converge for 421 groups (spectral
   radius ~0.97).  BiCGSTAB fully converges the inner solve in ~100
   Krylov iterations.

Matching the MATLAB Result
---------------------------

The MATLAB code uses: 2D Lebedev 110 on a 10x2 mesh, explicit FD
operator with BiCGSTAB, boundary half-volumes, P0 scattering.

The BiCGSTAB path with half-volumes reproduces 1.04192 vs MATLAB's
1.04188 (:math:`4 \times 10^{-5}` agreement).  The residual difference
is from floating-point details in cross-section processing.

The cleanest reference is the **1D GL BiCGSTAB** result (1.03882): no
pseudo-2D artifacts, well-conditioned angular quadrature, fully
converged inner solve.


References
==========

.. [Bailey2009] T.S. Bailey, J.E. Morel, and J.H. Chang,
   "A piecewise-linear finite element discretization of the diffusion
   equation for arbitrary polyhedral grids,"
   *Nuclear Science and Engineering*, 162:3, 2009.
   (Eq. 50: :math:`\alpha` recursion; Eq. 53--54: WDD;
   Eq. 74: Morel--Montry weights.)

.. [MorelMontry1984] J.E. Morel and G.R. Montry,
   "Analysis and elimination of the discrete-ordinates flux dip,"
   *Transport Theory and Statistical Physics*, 13:5, 1984.

.. [LewisMiller1984] E.E. Lewis and W.F. Miller, Jr.,
   *Computational Methods of Neutron Transport*,
   John Wiley & Sons, 1984.

.. [CaseZweifel1967] K.M. Case and P.F. Zweifel,
   *Linear Transport Theory*,
   Addison-Wesley, 1967.

.. [Lebedev1999] V.I. Lebedev and D.N. Laikov,
   "A quadrature formula for the sphere of the 131st algebraic order
   of accuracy," *Doklady Mathematics*, 59(3):477--481, 1999.

.. [CarlsonLathrop1965] B.G. Carlson and K.D. Lathrop,
   "Transport theory -- the method of discrete ordinates,"
   in *Computing Methods in Reactor Physics*,
   Gordon and Breach, 1968.

.. [AdamsLarsen2002] M.L. Adams and E.W. Larsen,
   "Fast iterative methods for discrete-ordinates particle transport
   calculations," *Progress in Nuclear Energy*, 40(1):3--159, 2002.
   Reviews the SAILOR / Larsen-Adams preconditioned-Krylov framework
   that the Wave E Round 2 inner solver implements.
