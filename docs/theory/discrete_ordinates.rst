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
- Key reference: Bailey, Morel & Chang (**2010**) NSE 165 ([BaileyMorelChang2010]_)
  — Eq. 43 (the M-M weight, unique exact-on-linear-in-μ); Hébert (2009)
  §3.9.4 ([Hebert2009]_) for the dome recursion + Carlson seed.
- **Curvilinear anisotropic SN** — the "#229 floor" was **three**
  distinct errors (sphere pole-cell O(h) spatial #233; sphere angular
  τ-clamp; cylinder angular floor #229), separated by a norm difference
  (volume-weighted L2 vs L∞).  The sphere τ-clamp is removed (W1, raw
  Bailey Eq. 43); the cylinder keeps it (structural τ=0).  Two unrelated
  "anisotropic" paths: geometric α-dome redistribution (P0) vs Legendre
  P1+ scattering (#9).  Full treatment:
  :ref:`sn-curvilinear-aniso-norm-reconciliation`.
- Verification uses :ref:`synthetic cross sections <synthetic-xs-library>`, not real nuclear data
- 2-D wavefront sweep (Wave 2): per-octant batched dispatch via
  :class:`~orpheus.sn.sweep_graph.SweepDependencyGraph`; mesh-time
  precompute of the per-octant DAG; BC apply once per octant per axis
  (the L7-trap fix). **Since S6.4(e)** the graph exposes TWO storage
  walks —
  :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full`
  (full-cochain oracle) and
  :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`
  (rolling-frontier production) — each parameterised by a LEVEL
  OPERATION object (``_CellSolve`` | ``_CellResidual``); the cell math
  is the discretization's storage-free kernel pair
  (:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
  / :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`).
  See :ref:`sweep-octant-dependency-graph`.
- **Both the 1-D and 2-D sweeps are BARE** (Wave O steps O.4a.2 +
  O.4b, Issue #208): :func:`~orpheus.sn.loss_representation.transport_sweep` (1-D)
  and :func:`~orpheus.sn.loss_representation._sweep_jacobi` (2-D) no longer
  re-apply ``bc`` to the outflow — they read the *seeded* inflow
  trace directly. The reflective coupling :math:`\psi.\text{inflow} =
  B\,\psi.\text{outflow}` is delivered as a sibling :math:`-B` source
  term, NOT re-derived inside the sweep, for every geometry. See
  :ref:`bare-sweep-extraction` and the canonical algebra
  :ref:`bc-extraction` in :doc:`operator_algebra`.
- **2-D Cartesian eigenvalue problems solve via BOTH inner solvers**
  (Wave O "2-D SI Phase A", Issue #208, 2026-06-04):
  :func:`~orpheus.sn.solver.solve_sn` runs the source-iteration inner
  (:meth:`~orpheus.sn.solver.SNSolver._solve_source_iteration`, the
  default for *every* geometry) AND the Krylov inner
  (:meth:`~orpheus.sn.solver.SNSolver._solve_krylov`). The SI inner is
  the geometry-agnostic structural twin of Krylov — same composite RHS,
  same loss decomposition (the resolvent :math:`L + C` plus the
  scattering gain :math:`S` and the boundary reflection gain :math:`B`,
  handed to the variadic driver; zero within-group fission), same
  angular reduction — differing only in the iteration driver
  (:class:`~orpheus.numerics.iteration.SourceIteration` vs
  :class:`~orpheus.numerics.iteration.KrylovAcceleration`). The
  reflective coupling rides the bare 2-D sweep via the sibling
  :math:`-B` on the natively four-face
  :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`. The legacy
  "B1'' face block" (never a code symbol) is retired. Verified
  SI ≡ Krylov ≡ closed-form :math:`k_\infty` (1g → 1.5, 2g → 1.875,
  4g → 1.4878); heterogeneous non-flat 2-D flux shape agrees SI-vs-Krylov
  to :math:`\sim 10^{-9}`. See
  :ref:`bc-extraction-2d-si-krylov-twin` in :doc:`operator_algebra`.
- **The 2-D interior cell-face fluxes are a 1-cochain**
  :math:`C^1_{\rm int}` (Wave O step #205 Phase 5, Issue #208,
  2026-06-04): the 2-D wavefront sweep + matvec no longer carry raw
  ephemeral ``psi_x`` / ``psi_y`` numpy arrays — they are the interior
  1-cochain :math:`C^1_{\rm int}`, and the boundary seed/absorb are the
  typed trace operators :math:`\iota_*` / :math:`\iota^*`. A dedicated
  ``WavefrontFlux`` field carried the cochain from #205 Phase 5 through
  S6.4; it was **retired at S6.4(f)** (#222) when the walk re-layering
  moved the seed/absorb verbs into the shared octant frame, and the
  cochain now lives in the rolling front
  (``_MovingFrontier``, ``orpheus.sn.sweep_graph``) and the
  full-cochain oracle history (``_octant_face_cochain``,
  ``orpheus.sn.loss_representation``).
  See :ref:`wavefront-flux-cochain` in :doc:`operator_algebra` for the
  succession.
- **Storage layout**: angular flux ``(N, ng, nx, ny)``; scalar flux
  and per-cell cross sections ``(ng, nx, ny)``; 1-D problems keep
  ``ny = 1`` (singleton, NOT squeezed).  The canonical statement
  with derivation and migration history lives at
  :ref:`theory-sn-index-convention`.

.. admonition:: Conventions

   - Scattering matrix: :ref:`scattering-matrix-convention` — ``SigS[g_from, g_to]``, source uses transpose
   - **Storage layout**: :ref:`theory-sn-index-convention` — ``(N, ng, nx, ny)`` for ψ, ``(ng, nx, ny)`` for φ / σ
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
consistency.  The treatment follows [BaileyMorelChang2010]_ for the
curvilinear formulation, [LewisMiller1984]_ for the general framework, and
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

2. **Augmented geometry** --- :class:`SNMesh` pairs the spatial mesh
   with an angular quadrature, precomputing the coordinate-specific
   streaming stencil.  Its **primary representation is the per-axis
   tuple** :attr:`SNMesh.axes` (the SN phase space factors as a tensor
   product of per-axis 1-D meshes): a legacy ``Mesh1D`` / ``Mesh2D`` is
   converted to axes **once** at the inbound boundary, and
   :meth:`SNMesh.from_axes` stores the caller's tuple verbatim. After
   C5 (:ref:`sn-axis-primary-c5`) :attr:`SNMesh.mesh` is *inbound
   provenance only* — ``None`` for an axis-native :math:`d \ge 3` mesh,
   which carries no legacy mesh at all.  It also **resolves boundary
   conditions**: each ``BC`` tag
   on the mesh is looked up in :attr:`SNMesh.BOUNDARY_OPERATOR_REGISTRY` and converted
   to a validated kind string (``"vacuum"`` or ``"reflective"``)
   stored in the face-name-keyed :attr:`SNMesh.bc` dict
   (``sn_mesh.bc["xmin"]``, ``sn_mesh.bc["xmax"]``, ... — the dict
   keys are the mesh's true boundary faces; see
   :ref:`bc-face-name-carve`).
   The sweep reads these resolved strings directly --- it never
   inspects the raw :class:`~geometry.mesh.BC` objects.  Precomputed
   stencil contents per coordinate system:

   - **Cartesian**: one per-axis array ``streaming(a)[n,i] =
     2|mu_a|/da[i]`` for every axis ``a < ndim`` (built over
     ``range(ndim)`` from ``quad.axis_cosines(a)`` since C3.6 ---
     no hand-listed x/y pair, no phantom axis on a slab) --- the
     diamond-difference denominator terms, precomputed to avoid
     per-cell division in the sweep hot loop.
   - **Spherical**: ``face_areas`` (:math:`4\pi r^2`), ``delta_A``,
     ``alpha_half`` (angular redistribution dome),
     ``redist_dAw`` (:math:`\Delta A_i / w_n`, shape ``(nx, N)``),
     and ``tau_mm`` (Morel--Montry closure weights).
   - **Cylindrical**: ``face_areas`` (:math:`2\pi r`), ``delta_A``,
     ``alpha_per_level`` (per-level redistribution domes),
     ``redist_dAw_per_level`` (list of ``(nx, M)`` arrays), and
     ``tau_mm_per_level`` (per-level Morel--Montry weights).

3. **Solver** --- :func:`solve_sn` creates an ``SNMesh``, builds the
   ``SNSolver``, and runs power iteration. At :math:`d \le 2` the input
   is a ``Mesh1D`` / ``Mesh2D``; at :math:`d = 3` the input is the
   **axes tuple** itself (the only 3-D entry — there is no ``Mesh3D``;
   see :ref:`sn-c5-3d-admission`).

.. code-block:: text

   Mesh1D / Mesh2D (d<=2)   OR   axes tuple (d=3, axis-native)
       |                              |
       |  axes_from_legacy_mesh       |  (stored verbatim)
       +------------+-----------------+
                    v
   SNMesh.axes  (PRIMARY)  -->  SNMesh (stencil + quadrature
                                + alpha coefficients + resolved BCs)
                    |
                    v
   solve_sn() --> SNResult

Quadrature Dispatch
-------------------

The sweep dispatcher in :func:`transport_sweep` routes based on the
``SNMesh.curvature`` attribute and the quadrature type.  Boundary
conditions are **not** passed as a parameter to the sweep --- the
sweep reads the resolved BC kind strings directly from the
face-name-keyed :attr:`SNMesh.bc` dict (``sn_mesh.bc["xmin"]``,
``sn_mesh.bc["xmax"]``, ...; see :ref:`bc-face-name-carve`):

.. code-block:: python

   if sn_mesh.curvature == "spherical":
       return _sweep_1d_spherical(...)
   elif sn_mesh.curvature == "cylindrical":
       return _sweep_1d_cylindrical(...)
   elif is_gl_1d:            # ny=1, mu_y=0, no aniso source
       return _sweep_1d_cumprod(...)
   else:
       return _sweep_jacobi(...)

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
:math:`\alpha` recursion convention from [BaileyMorelChang2010]_ (the
curvilinear :math:`\alpha`-recursion).

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
:class:`SNMesh` as ``streaming(0)[n, i]``.

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

These are precomputed by :class:`SNMesh` as ``streaming(0)[n, i]`` and
``streaming(1)[n, j]``, so the inner loop in
:func:`_sweep_jacobi` reduces to a single vectorised division per
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
consistency.  The correct form from [BaileyMorelChang2010]_ includes the
geometry factor :math:`\Delta A_i / w_n`:

.. math::
   :label: balance-general

   \mu_n
   \bigl[A_{i+\frac12}\psi_{i+\frac12} - A_{i-\frac12}\psi_{i-\frac12}\bigr]
   + \frac{\Delta A_i}{w_n}
   \bigl[\alpha_{n+\frac12}\psi_{n+\frac12} - \alpha_{n-\frac12}\psi_{n-\frac12}\bigr]
   + \Sigt{} V_i \psi_{n,i} = S_i V_i

where :math:`\Delta A_i = A_{i+1/2} - A_{i-1/2}`.  This is the curvilinear
balance form of [BaileyMorelChang2010]_ for both spherical and cylindrical
geometry.

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

This is the [BaileyMorelChang2010]_ curvilinear :math:`\alpha`-recursion.
Each level's :math:`\alpha` values form an independent dome from
:math:`\eta = -\sin\theta` to :math:`\eta = +\sin\theta`.

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

The Streaming-Equilibrium Identity (canonical L0 gate)
-------------------------------------------------------

The flat-flux consistency proof above is the *per-cell* statement of a
*global* exact solution that the verification suite leans on harder
than any other: for a *homogeneous* medium with a *uniform isotropic*
source :math:`Q` per group and boundaries that sustain a flat flux
(reflective faces, or an infinite/periodic medium), the discrete
transport equation has the exact fixed point

.. math::
   :label: streaming-equilibrium

   \phi \;=\; \frac{Q}{\Sigma_t\,\bigl(1 - c\bigr)},
   \qquad
   \psi_n \;=\; \frac{\phi}{W}
   \quad \forall n,
   \qquad
   W \equiv \sum_n w_n,
   \quad
   c \equiv \frac{\Sigma_{s0}}{\Sigma_t},

per group. For a pure-attenuation configuration (no scattering in the
residual, :math:`c = 0`) this reduces to :math:`\phi = Q/\Sigma_t` with
the per-ordinate angular flux :math:`\psi_n = Q/(W\,\Sigma_t)`.

**Why this is the canonical L0 gate.** Substituting the flat
:math:`\psi_n` into the discrete balance :eq:`balance-general` nulls
the streaming and redistribution terms *per ordinate* (the proof of
consistency above), leaving the pure collision balance
:math:`\Sigma_t\,\psi_n = Q/W + \Sigma_{s0}\,\phi/W`. Every term that
a discretisation can get wrong — a missing :math:`\Delta A/w` factor
(failure mode #3, the flux spike at :math:`r=0`), a wrong
:math:`\alpha` recursion (mode #4), a face-index slip (mode #5), a
weight-normalisation drift (:math:`1/W` vs :math:`1/4\pi`) — breaks
the identity at machine precision, with no discretisation error to
hide behind. The assertion is **per-ordinate**, never
particle-balance: telescoping global balance holds by construction
even when per-ordinate balance is wrong (the canonical ERR-006 hide;
vv-principles anti-pattern #8).

The identity holds in every geometry ORPHEUS supports (slab, sphere,
cylinder, 2-D/3-D Cartesian) and at both algebraic access points —
the sweep (``solve``: given :math:`Q`, recover the flat
:math:`\psi`) and the matvec (``apply``: given the flat :math:`\psi`,
recover the residual :math:`Q` with no spurious boundary or pole
contribution). Tests declare it via
``@pytest.mark.verifies("streaming-equilibrium")``.

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

The contamination factor :math:`\beta` ([BaileyMorelChang2010]_) quantifies
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
with position-dependent weights :math:`\tau_n` [BaileyMorelChang2010]_ Eq. 43:

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

**Clipping (geometry-dependent since W1, 2026-06-13).**  The raw
weight :eq:`mm-weights` is the unique weight exact for a flux linear in
:math:`\mu` (Bailey-Morel-Chang 2010 Eq. 43), admissible range
:math:`\tau \in [0, 1]`.

* **Sphere** uses :math:`\tau_n` **unclamped**.  On Gauss--Legendre
  quadrature :math:`\tau_n \in [0.39, 0.61]` (never 0), so the closure
  is positive without a clamp.  The former :math:`[0.5, 1.0]` clamp was
  an over-conservative, mis-cited positivity floor that re-floored the
  anisotropic solution; the full vindication is at
  :ref:`sn-curvilinear-aniso-norm-reconciliation`.
* **Cylinder** retains the clip :math:`\tau_n \to
  \mathrm{clip}(\tau_n, 0.5, 1.0)`: product / level-symmetric
  quadratures give :math:`\tau_n = 0` exactly on the most-inward
  azimuthal ordinate (a structural :math:`\div 0` block; the alternating
  0.5 / 1.0 pattern described above is the clamp acting on these
  zero-width angular cells).

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
:mod:`orpheus.sn.spatial.scheme`.

The Protocol
------------

The contract is a ``@runtime_checkable`` ``typing.Protocol`` —
satisfied by structural typing, not inheritance — exposing two class-
level traits and a single :meth:`update` method:

* :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme`

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
    :class:`~orpheus.sn.spatial.scheme.CellVisit` packet (see
    next subsection) that combines the geometric
    :class:`~orpheus.geometry.reduced_operator.StreamingTerms` with
    sweep-direction-resolved data.

The two helper dataclasses (frozen, slotted) carry the per-cell
state:

* :class:`~orpheus.sn.spatial.scheme.UpstreamState`

  - ``spatial_upstream: np.ndarray`` — shape ``(ng,)``.  Face flux
    entering the cell from the upstream face (always populated).
  - ``angular_upstream: np.ndarray | None`` — shape ``(ng,)`` for
    sphere/cylinder; ``None`` for slab.  :math:`\psi_{n-1/2,\,i}`,
    the half-flux at the upstream half-angle.

* :class:`~orpheus.sn.spatial.scheme.CellResult`

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
:class:`~orpheus.sn.spatial.scheme.CellVisit` packet rather
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
:meth:`SNMesh.dag_walk(*, ordinate_idx=..., direction_sign=...,
mu_level_idx=None)
<orpheus.sn.geometry.SNMesh.dag_walk>` — a generator that
yields cells in DAG-topological order.  The method takes EXACTLY ONE
of ``ordinate_idx`` (single-ordinate visits, used by the sweep
driver) or ``direction_sign`` (direction-keyed visits, used by the
apply matvec) as a keyword-only argument.  Both invocation modes
encapsulate the inward / outward branching, the cylindrical
per-:math:`\mu`-level traversal, and the pure-azimuthal degenerate
handling.  The sweep at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) consumes this
generator::

    for visit in sn_mesh.dag_walk(ordinate_idx=n):
        upstream = UpstreamState(
            spatial_upstream=psi_face,
            angular_upstream=psi_angle[visit.cell_idx],
        )
        result = scheme.update(
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
``tests/sn/sweep/core/test_discretization_scheme_protocol.py``; concrete
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
:meth:`SNMesh.dag_walk
<orpheus.sn.geometry.SNMesh.dag_walk>` yields visits with
``face_area_downstream = 0.0`` to signal "no spatial flow" to the
strategy (Issue #196 Step 2.5 retired the ``None`` sentinel — the
slab carries ``1.0`` and degenerate cylindrical carries ``0.0`` so
the cell-balance helper consumes one geometry-blind number).

The DD recurrence
------------------

For non-degenerate cells, the closure relation reduces — for slab
geometry, the cell-update math is the DD recurrence
:eq:`dd-recurrence` (see :ref:`sweep-cumprod` below); for curvilinear
geometry, it is the curvilinear DD form combining the
:math:`\Delta A/w` redistribution with the M-M angular closure.  The
sweep driver inlines this math today; Wave D (Issue #159) will
rewrite the driver to dispatch through a
:class:`~orpheus.sn.spatial.scheme.DiscretizationScheme` strategy.

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
formulas at ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``).

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
(``_sweep_1d_cumprod`` (the dissolved ``sweep.py``) lines 117–123) and per-
cell solver (``_solve_recurrence`` (the dissolved ``sweep.py``) lines 208–
222) at the operation level.  Per the strategy contract, ``source``
arrives at the cell update **already weight-normalised** by the
sweep — for slab, ``source = Q · Δx / W`` (and slab cell volume is
``V = Δx``).  For slab, the strategy sets
:attr:`~orpheus.sn.spatial.scheme.CellResult.outgoing_angular_state`
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

mirroring ``_sweep_1d_spherical`` (the dissolved ``sweep.py``) lines
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

mirroring ``_sweep_1d_cylindrical`` (the dissolved ``sweep.py``) lines
533–543 verbatim.  The strategy returns
:attr:`~orpheus.sn.spatial.scheme.CellResult.outgoing_spatial_flux`
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

* :mod:`orpheus.sn.spatial.scheme` — the contract module.
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

.. _sn-affine-outgoing-face-reconstruction:

Generic affine outflow reconstruction
--------------------------------------

.. math::
   :label: sn-affine-outgoing-face-reconstruction-eq

   \psi_{\rm out} = \frac{\bar\psi - (1-w)\,\psi_{\rm in}}{w}

The single-source inverse of the convex cell-average blend
:math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`
(:eq:`dd-recurrence` closure :math:`\psi_{\rm out} = 2\bar\psi - \psi_{\rm in}`
is the :math:`w=\tfrac12` case).  Every consistent affine spatial scheme
reconstructs its downstream face from this one parameterized formula: Diamond
Difference at :math:`w=\tfrac12` (diamond mean), Linear Discontinuous at
:math:`w = 1/(1+k)` with :math:`k = (|\mu|/\theta)/D_2`.  At :math:`w=\tfrac12`
the reconstruction is byte-identical to the inlined :math:`2\bar\psi - \psi_{\rm in}`
(division by :math:`\tfrac12` is the exact power-of-two doubling, which commutes
with round-to-nearest); LD's :math:`w=1/(1+k)` is algebraically equal to its
inlined Schur form :math:`\bar\psi + (|\mu|/\theta)(\bar\psi - \psi_{\rm in})/D_2`
but takes a different floating-point reduction tree (a principled
:math:`\sim`\ 1-ULP re-baseline).

The one parameterized formula above is the **single source** of the
downstream-face reconstruction for *every* affine 1-D spatial scheme.  It is
homed (with its forward partner, the cell-average blend
:math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`, and the source
emission) as a ``@staticmethod`` on the
:class:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase`, so the per-scheme
classes (:class:`~orpheus.sn.spatial.diamond.DiamondDifference`,
:class:`~orpheus.sn.spatial.linear_discontinuous.LinearDiscontinuous`, Step)
carry NO inlined reconstruction of their own — they differ only by the value
they pass for the blend weight :math:`w`.  The #240 Phase 2 Step D1 carve
collapsed the previously-duplicated inline forms (Diamond Difference's
:math:`2\bar\psi - \psi_{\rm in}`, Linear Discontinuous's
:math:`(1+k)\bar\psi - k\,\psi_{\rm in}`) onto this one op
(:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`),
the algebraic inverse of
:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_average`; Step D2
made the trio (``source_emission`` / ``cell_average`` /
``outgoing_face_from_average``) generic advection–reaction reconstructions
(diffusion-consumable, retiring the dangling ``affine_closure`` module).  The
unit gate is :mod:`tests.sn.spatial.test_affine_closure`: the exact-inverse
round-trip :math:`\bar\psi(\,\psi_{\rm in}, \psi_{\rm out}(\psi_{\rm in}, \bar\psi)\,) = \bar\psi`,
the DD :math:`w = \tfrac12` byte-identity, and the LD :math:`w = 1/(1+k)`
algebraic equality.

The Step / DD / LD ladder is the transport-sweep face of a correspondence
familiar from finite-volume advection: each scheme is one choice of how the
**downstream** face value is reconstructed from the cell average and the
**upstream** face value, parameterized by the convex blend weight :math:`w`.

.. list-table:: The affine spatial schemes as advection-flux reconstructions
   :header-rows: 1
   :widths: 16 14 30 40

   * - Scheme
     - Blend :math:`w`
     - Face reconstruction
     - Advection-scheme analog
   * - **Step** (upwind)
     - :math:`w \to 1`
     - :math:`\psi_{\rm out} = \bar\psi` (cell value carried to the
       downstream face)
     - first-order **upwind**: the outgoing flux IS the cell-centre value, no
       slope.  Monotone, unconditionally positive, :math:`O(h)`.
   * - **Diamond Difference**
     - :math:`w = \tfrac12`
     - :math:`\psi_{\rm out} = 2\bar\psi - \psi_{\rm in}` (central / box mean)
     - **box / central** difference: the cell average is the arithmetic mean of
       the two faces.  :math:`O(h^2)`, but can go negative (the DD flux dip).
   * - **Linear Discontinuous**
     - :math:`w = 1/(1+k)`,
       :math:`k = (|\mu|/\theta)/D_2`
     - :math:`\psi_{\rm out} = \bar\psi + \hat\psi` (DG-P1 upstream trace)
     - **DG-P1 with upwind face flux**: a discontinuous linear profile per cell,
       the downstream face is the cell's own :math:`P_0 + P_1` trace.
       :math:`O(h^2)`, diffusion-limit-consistent (slope-source-fed).

The weight :math:`w = 1/(1+k)` is a **Péclet-type blend**: writing
:math:`k = (|\mu|/\theta)/D_2` (the stream-to-collision ratio of
:eq:`ld-ubld-scale-free-invariants`, scaled by the Legendre normalization
:math:`\theta = \tfrac13`), the optically thin cell (:math:`\Sigma_t h \to 0`,
streaming-dominated, :math:`k \to \infty`) drives :math:`w \to 0` so the
downstream face approaches the *upstream* face (the streaming limit, no
collisional re-equilibration within the cell); the optically thick cell
(:math:`\Sigma_t h \to \infty`, collision-dominated, :math:`k \to 0`) drives
:math:`w \to 1` so the downstream face approaches the cell average (the
diffusive limit).  This is exactly the role the Péclet number plays in a
:math:`\kappa`-scheme advection blend — :math:`w` interpolates between upwind
(:math:`w=1`) and central (:math:`w=\tfrac12`) as a function of the cell
optical thickness.  Diamond Difference's fixed :math:`w=\tfrac12` is the
:math:`\Sigma_t = 0` (pure streaming) value of the LD weight only in the
degenerate sense that DD ignores the within-cell collision balance the LD
slope resolves; the two coincide on the *average* but DD has no slope row at
all.

.. note::

   **The spatial closure factors out of the angular index, and the multi-D
   extension factors out of the dimension.**  The reconstruction op above is
   stated per ordinate per axis, so it is a *spectator* to the angular moment
   axis (the angular reduction :eq:`two-moment-angular` rides over it, see
   :ref:`two-moment-axes`) — the same op serves a P0 and a P3 calculation
   unchanged.  In the same spirit, the multi-dimensional LD closure
   (:ref:`ld-ubld-multidim`, next) is the **tensor product** of this 1-D
   per-axis reconstruction across :math:`d` axes: the per-cell
   :math:`2^d \times 2^d` operator is assembled as a Kronecker product of the
   verified 1-D factor operators, so the affine 1-D closure documented here is
   the literal :math:`d=1` building block of the :math:`d`-generic UBLD system
   :eq:`ld-ubld-cell-system`.  Spatial scheme :math:`\otimes` angular order
   :math:`\otimes` dimension: three orthogonal axes of choice, each a tensor
   factor, none special-casing the others.

.. _ld-ubld-multidim:

Multi-dimensional LD: the tensor-product bilinear (UBLD) cell system
---------------------------------------------------------------------

The multi-dimensional analog of Linear Discontinuous on a **Cartesian**
cell is **NOT** the simplex-P1 :math:`\{1, x, y\}` object
(:math:`1+d` moments).  Adams (2001) proved simplex-LD *fails* the thick
diffusion limit on quadrilaterals, while the **bilinear / trilinear
DG-P1** (UBLD) — basis :math:`\{1, x, y, xy\}` (:math:`2^d` moments) —
*passes*.  The :math:`xy` cross moment is diffusion-limit-load-bearing.

The :math:`d`-generic per-cell Galerkin system is assembled as Kronecker
products of the verified 1-D LD factor operators (the streaming
:math:`\Omega\cdot\nabla = \sum_a \mu_a \partial_a` is a sum over axes;
the tensor-product basis separates):

.. math::
   :label: ld-ubld-cell-system

   A\,\vec\psi = \vec R, \qquad
   A = G + F_{\rm out} + \Sigma_t M, \qquad
   \vec R = M\,\vec S + F_{\rm in}\,\psi_{\rm in}^{\rm traces},

a :math:`2^d \times 2^d` dense non-symmetric solve, with
:math:`M = M_1 \otimes \cdots \otimes M_d` (mass),
:math:`G = \sum_a \mu_a\,(M_1 \otimes \cdots \otimes G_{1d} \otimes
\cdots \otimes M_d)` (streaming: gradient on the active axis, mass on the
transverse axes), and :math:`F_{\rm out}` likewise from the per-axis
downstream-face trace.

.. math::
   :label: ld-ubld-d1-reduction

   A\big|_{d=1} =
   \begin{bmatrix} \Sigma_t h + |\mu| & |\mu| \\
                   -|\mu| & \Sigma_t \theta h + |\mu| \end{bmatrix}

The :math:`d=1` reduction (Kronecker-with-one-factor identity) recovers
the production slab 2×2 :eq:`dd-cartesian-1d`-sibling exactly; the
:math:`xy` coupling falls out of the algebra for :math:`d \ge 2`.

.. math::
   :label: ld-ubld-exact-on-bilinear

   \psi(x,y) = a + bx + cy + dxy
   \;\Longrightarrow\;
   \vec\psi_{\rm solved} = \vec\psi_{\rm exact-projections}

The UBLD is **exact on any bilinear flux** (the multi-D analog of the
1-D "exact on linear-in-x" oracle), the :math:`xy` cross moment
exercised — the structurally-independent correctness gate for the
:math:`d \ge 2` closure.

The Branch-1 algebra-of-record (the UBLD weak form)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The canonical symbolic reference for the :math:`d`-generic UBLD system is the
SymPy module :mod:`orpheus.derivations.discrete.sn.ld_ubld` (the
algebra-of-record, State 1A closed-form): the Kronecker assembler
``assemble_ubld`` plus five ``derive_*`` verification functions, each proven by
``sympy.simplify(diff) == 0``.  It is the discrete-SN sibling of
``orpheus.derivations.discrete.sn.balance`` — a *symbolic discretization the
production solver must satisfy*, NOT a continuous reference.

The per-cell system descends from the Galerkin weak form of the within-group
streaming–collision operator (Maginot, Ragusa & Morel 2016, "Non-negative
Methods for Bilinear Discontinuous Differencing of the :math:`S_N` Equations on
Quadrilaterals", NSE 185(1):17–42, Eqs. 1–12).  Multiplying the transport
equation :math:`\Omega\cdot\nabla\psi + \Sigma_t\psi = S` by each basis function
:math:`B_i` and integrating over the cell :math:`K`, then integrating the
streaming term by parts (MRM-2016 Eq. 6),

.. math::
   :label: ld-ubld-weak-form

   \underbrace{(\Omega\cdot)\!\oint_{\partial K} \hat n\,B_i\,\psi\,d\ell}_{\text{surface (upwind)}}
   \;-\; \int_K \psi\,(\Omega\cdot\nabla B_i)\,dV
   \;+\; \Sigma_t\!\int_K B_i\,\psi\,dV
   \;=\; \int_K B_i\,S\,dV,

gives three operators per cell — the **mass** :math:`M_{ij} = \int_K B_i B_j`
(the collision term), the **gradient/stiffness** :math:`G_{ij} = \int_K B_i\,(\Omega\cdot\nabla B_j)`
(the volumetric streaming term, coupling all :math:`2^d` moments), and the
**surface** matrix split per face into an OUTFLOW block (:math:`\Omega\cdot\hat n > 0`,
implicit, the cell's own unknowns) and an INFLOW block (:math:`\Omega\cdot\hat n < 0`,
**upwind**: the incoming face value is the upstream neighbour's outflow trace,
moved to the RHS).  Assembling gives exactly the dense system
:eq:`ld-ubld-cell-system`, :math:`A = G + F_{\rm out} + \Sigma_t M`,
:math:`\vec R = M\vec S + F_{\rm in}\,\psi_{\rm in}^{\rm traces}`.

Why bilinear, not simplex-P1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The naïve multi-D analog of 1-D LD — "cell average plus one slope per axis",
the simplex-P1 basis :math:`\{1, x, y\}` of :math:`1+d` moments — is the
**wrong object on a Cartesian cell**.  Adams (2001, "Discontinuous Finite
Element Transport Solutions in Thick Diffusive Problems", NSE 137(3):298–333)
proved that simplex-LD *fails* the thick diffusion limit on quadrilaterals,
while the **bilinear / trilinear DG-P1** (UBLD) — the tensor-product basis
:math:`\{1, x, y, xy\}` of :math:`2^d` moments — *passes* it.  The reason is the
:math:`xy` **cross moment**: it is exactly what the simplex basis lacks, and it
is the term the leading-order asymptotic diffusion balance needs (Börgers,
Larsen & Adams 1992, "The asymptotic diffusion limit of a linear discontinuous
discretization of a two-dimensional linear transport equation", JCP
98(2):285–300, give the 2-D rectangular analysis explicitly).  The simplex-P1
basis *does* preserve the limit on a genuine simplex (triangle/tetra) mesh
(Wareing, McGhee, Morel & Pautz 2001, NSE 138(3):256–268) — but that is a
different cell topology, not a quadrilateral.  ORPHEUS builds Cartesian cells,
so the :math:`2^d` tensor-product object is the diffusion-limit-consistent
choice; the choice is **load-bearing**, not a convenience (see
:ref:`ld-ubld-scattering-moment-lift` for the companion half of the same
asymptotic argument).

The Kronecker single-source build
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The three matrices are NOT hand-transcribed entry-by-entry (that would be a
:math:`4\times4` / :math:`8\times8` transcription waiting for a sign error).
They are assembled as **Kronecker products of the verified 1-D LD factor
operators** in the Legendre moment basis :math:`\{1, P_1\}` on width :math:`h`:

.. math::
   :label: ld-ubld-kronecker-factors

   M_{1d} = \mathrm{diag}(h,\ \theta h),
   \qquad
   G_{1d} = |\mu|\begin{bmatrix} 0 & 0 \\ -2 & 0 \end{bmatrix},
   \qquad
   F_{\rm out}^{1d} = |\mu|\begin{bmatrix} 1 & 1 \\ 1 & 1 \end{bmatrix},

with the streaming :math:`\Omega\cdot\nabla = \sum_a \mu_a\,\partial_a` a sum
over axes (the tensor-product basis separates), so

.. math::
   :label: ld-ubld-kronecker-assembly

   M = M_1 \otimes \cdots \otimes M_d,
   \qquad
   G = \sum_a \mu_a\,(M_1 \otimes \cdots \otimes G_{1d}^{(a)} \otimes \cdots \otimes M_d),

i.e. the gradient acts on the active axis and **the mass on every transverse
axis** (the volume-integral factorization — this is the load-bearing build
choice).  The :math:`F_{\rm out}` surface matrix is assembled likewise; the
inflow is a :math:`B(-1) = [1, -1]` test-weighting on the active axis (mass on
transverse axes) times :math:`|\mu_{\rm axis}|`.  The diagonal mass weights are
then the Kronecker product of the per-axis diagonals — a power of
:math:`\theta = \tfrac13` equal to the **number of active (slope) axes** of each
moment:

.. math::
   :label: ld-ubld-mass-weights

   M_{ii} = \theta^{|i|},
   \qquad
   |i| = \#\{a : o_a = 1\}
   \;\Longrightarrow\;
   \begin{cases}
     1        & \bar\psi \quad (\text{no slope axis}) \\
     \theta   & \hat\psi_x,\ \hat\psi_y \quad (\text{one slope axis}) \\
     \theta^2 & \hat\psi_{xy} \quad (\text{two slope axes})
   \end{cases}

so the 2-D weights are :math:`(1, \theta, \theta, \theta^2)` and the 3-D
:math:`xyz` cross moment carries :math:`\theta^3`.  These weights re-appear in
the matvec mass-normalization (:eq:`ld-ubld-unified-moment-residual`) — they are
the SAME diagonal Legendre mass.  The :math:`d=1` case is a Kronecker product
with a single factor (an identity), so it reduces EXACTLY to the production
slab 2×2 :eq:`ld-ubld-d1-reduction`; the :math:`xy` coupling *emerges* from the
algebra for :math:`d \ge 2` — no entry is hand-written.

The two oracles
^^^^^^^^^^^^^^^

The module proves the construction with two structurally distinct oracles
(both ``sympy.simplify(diff) == 0``):

* **Oracle (i) — the :math:`d=1` reduction.**
  ``derive_d1_reduction_to_production`` shows the assembled :math:`d=1` system
  equals the production
  :mod:`orpheus.sn.spatial.linear_discontinuous` 2×2 entry-for-entry, with the
  Schur complement :math:`S` and the effective slope denominator
  :math:`D_2' = \Sigma_t h\theta + |\mu|` recovered as the production closed
  forms.  Two further reductions
  (``derive_d1_kernel_view_equals`` / ``derive_d1_scan_view_equals``) prove the
  same :math:`d=1` reduces to BOTH the ÷V DAG kernel ``_kernel_terms`` and the
  ×V scan ``affine_scan_coefficients`` views — the "single-source the math"
  proof that Branch 2's three production views are the SAME algebra
  (:eq:`ld-ubld-rule-of-three-collapse`).

* **Oracle (ii) — exact-on-bilinear at :math:`d=2`.**
  ``derive_d2_exact_on_bilinear`` feeds an exactly-bilinear flux
  :math:`\psi = a + bx + cy + dxy` through the DG-exact upstream face moments and
  the projected source moments, and shows the 4 solved moments equal the exact
  projections (:eq:`ld-ubld-exact-on-bilinear`).  The :math:`xy` cross moment is
  genuinely exercised (:math:`d \ne 0` symbolically) — this is the multi-D
  analog of the 1-D "exact on linear-in-x" oracle and the structurally
  independent correctness gate for the :math:`d \ge 2` closure.

The foundation gate is :mod:`tests.sn.spatial.test_ld_ubld_symbolic` (6
``@pytest.mark.foundation`` tests, one per ``derive_*`` plus an anchor to the
live production ``LinearDiscontinuous.update``); the literature contract is
recorded in
``.claude/agent-memory/literature-researcher/multi_d_ld_closure.md`` (MRM-2016
Eqs. 1–12; the Adams-2001 thick-diffusion verdict; BLA-1992); the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s1_ld_ubld_branch1_closeout.md``.

.. admonition:: ERR-060 — the oracle that earned its keep
   :class: tip

   The first draft of ``assemble_inflow_axis`` dropped the :math:`|\mu_{\rm axis}|`
   streaming factor on the inflow RHS (failure Mode 3, a missing factor).  The
   bug was INVISIBLE to all three :math:`d=1` oracles — the :math:`d=1` RHS is
   built inline, never routed through the multi-axis inflow assembler — and was
   caught by Oracle (ii), the :math:`d=2` exact-on-bilinear gate, which is the
   first consumer of ``assemble_inflow_axis``.  Mutation-verified
   :math:`-O`-safe: re-dropping the factor turns the :math:`d=2` test red while
   the :math:`d=1` tests stay green (proving they are blind to the bug class).
   This is the algebra-of-record discipline working as designed — the bug never
   reached production.  (The ERR-060 marker belongs on the *exact-on-bilinear*
   gates, NOT on the cell-matrix ``A == A`` pin, which checks
   ``assemble_ubld``'s A/M/G/:math:`F_{\rm out}` and is structurally blind to the
   dropped inflow factor — see ``error_catalog.md`` ERR-060.)

.. _ld-ubld-branch2-primitive:

The Branch-2 production primitive + the single-sourced d=1 fast path
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The numpy production counterpart of the symbolic algebra-of-record above
is :mod:`orpheus.sn.spatial._ubld`, in two layers that share ONE source
of truth.  Layer 1 is the :math:`d`-generic dense primitive
(``assemble_ubld`` / ``per_cell_solve``): a batched-over-cells Kronecker
build of the :math:`2^d \times 2^d` system :eq:`ld-ubld-cell-system`,
solved with a batched :func:`numpy.linalg.solve`.  It is the CANONICAL
:math:`d`-generic source for both :math:`d=1` (today) and :math:`d \ge 2`
(S2 wires the bilinear cell-batch kernel onto it); in production
:math:`d=1` does **not** route through this dense solve (that would be
the per-cell-solve performance regression).

Layer 2 is the shared :math:`d=1` closed form ``d1_closed_form`` — the
analytic Schur complement of the primitive's :math:`d=1` 2×2, VECTORIZED
over the cell / ordinate / group stack (no dense solve), so the
production :math:`d=1` path stays on the fast path.  The entire closure
rides two SCALE-FREE invariants:

.. math::
   :label: ld-ubld-scale-free-invariants

   k = \frac{g/\theta}{g/\theta + \Sigma_t}, \qquad
   w = \frac{1}{1 + k}, \qquad g = \frac{|\mu|\,A_{\rm down}}{V},

with ``w`` the cell-average blend weight
(:math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`).  Every
production view's coefficients are an algebraic function of
:math:`(g, \Sigma_t, k, w)` times a power of the cell volume :math:`V`
(the ×V vs ÷V choice applied at the call site).

.. math::
   :label: ld-ubld-rule-of-three-collapse

   \texttt{\_schur\_terms}\;(\times V), \quad
   \texttt{\_kernel\_terms}\;(\div V), \quad
   \texttt{affine\_scan\_coefficients}\;(\times V\ \text{scan})
   \;\longleftarrow\; \texttt{d1\_closed\_form}

The three production 1-D views in
:mod:`orpheus.sn.spatial.linear_discontinuous` — the ×V per-cell Schur
(``_schur_terms``), the ÷V DAG kernel (``_kernel_terms``), and the ×V
scan (``affine_scan_coefficients``) — now ALL derive their coefficients
from ``d1_closed_form``, applying their ×V / ÷V / ×V-scan scaling at the
call site.  The LD 2×2 algebra (the Rule-of-Three) collapses to ONE
place, proven ``==`` the dense primitive's :math:`d=1` reduction
(symbolically by the Branch-1 oracles, numerically end-to-end by the
Branch-2 gate).

The numpy production counterpart descends from the SAME specialized SymPy
ancestor as the Branch-1 algebra-of-record above; only the evaluation strategy
differs (Branch 1 closes the algebra symbolically, Branch 2 evaluates it on
arrays).  The discipline is **construct general, select narrow, specialize only
on measured cost**:

* **Construct general — the dense primitive.**  ``assemble_ubld`` /
  ``per_cell_solve`` build and solve the :math:`2^d \times 2^d` system
  :eq:`ld-ubld-cell-system` for every :math:`d`, batched over the cell /
  ordinate / group stack with a single :func:`numpy.linalg.solve`.  This is the
  canonical :math:`d`-generic source — :math:`d=1` (today), :math:`d \ge 2`
  (S2 wires the bilinear cell batch onto it), :math:`d = 3` (trilinear).

* **Select narrow — the :math:`d=1` closed form.**  ``d1_closed_form`` /
  :class:`~orpheus.sn.spatial._ubld.D1ClosedForm` is the analytic Schur
  complement of the primitive's :math:`d=1` 2×2, VECTORIZED over the stack with
  no dense solve.  Both scale-free invariants of :eq:`ld-ubld-scale-free-invariants`
  drive it.

* **Specialize on measured cost.**  In production :math:`d=1` does **not**
  route through the dense solve — that would be the per-cell-solve performance
  regression (the L16 constraint).  The closed form keeps the production
  :math:`d=1` sweep on the vectorized fast path
  (:class:`~orpheus.sn.loss_representation.CumprodScan` rides the ×V scan view's
  :math:`(a, \mathrm{inverse\_denom}, w)`; the DAG kernel rides the ÷V arrays).

The Rule-of-Three collapse
^^^^^^^^^^^^^^^^^^^^^^^^^^

Before the carve, the LD 2×2 algebra was transcribed in three production views
that had drifted into three independent copies.  All three now derive their
coefficients from the single ``d1_closed_form`` source
(:eq:`ld-ubld-rule-of-three-collapse`), applying only their scaling at the call
site — the Cardinal-Rule-2 / `coding-elegance` Pattern-2 single-source collapse:

.. list-table:: The three production 1-D views — one algebra, three scalings
   :header-rows: 1
   :widths: 30 16 54

   * - Production view (in :mod:`orpheus.sn.spatial.linear_discontinuous`)
     - Scaling
     - Consumer
   * - ``_schur_terms``
     - :math:`\times V`
     - the per-cell Schur (the matvec / ``update`` / ``residual`` path)
   * - ``_kernel_terms``
     - :math:`\div V`
     - the scale-free DAG wavefront kernel (the :math:`d \ge 2` arm rides the
       ÷V system, :eq:`ld-ubld-divv-scale-free-kernel`)
   * - ``affine_scan_coefficients``
     - :math:`\times V` scan
     - the Blelloch parallel-prefix scan (the production :math:`d=1` sweep)

The ×V / ÷V / ×V-scan choice is the volume scaling applied to the same
coefficients: dividing the Galerkin balance by the cell volume :math:`V` leaves
a scale-free system in the per-axis streaming :math:`g_a = |\mu_a| A_{\rm down}/V`
and :math:`\Sigma_t` alone (the form the :math:`d \ge 2` kernel consumes — fed
unit widths and :math:`\mathrm{mus} = (g_0, \ldots)`, it reduces EXACTLY to
``d1_closed_form``); multiplying restores the volume-weighted per-cell Schur
(:math:`D_2' = \theta V\,d_2`, :math:`S_{\times V} = V\cdot\mathrm{eff\_denom}`).
Each view is proven ``==`` the dense primitive's :math:`d=1` reduction —
symbolically by the Branch-1 oracles above (``derive_d1_kernel_view_equals`` /
``derive_d1_scan_view_equals``), numerically end-to-end by the Branch-2 gate.

.. note::

   **A principled ~1-ULP re-baseline, not a bit-identity break.**  Routing the
   three LD views through the shared helper changes the floating-point
   *reduction tree* relative to the legacy inline associations: the helper
   computes :math:`(g, \Sigma_t, k, w)` once and forms each coefficient as an
   algebraic function of them, whereas the old inline code interleaved the
   multiplies and adds differently.  In exact arithmetic the values are
   identical; in IEEE-754 they differ at ~1 ULP because addition is not
   associative.  This satisfies all three `vv-principles` criteria for
   accepting a non-bit-exact change: every intermediate
   (:math:`g`, :math:`k`, :math:`w`) is a *named, inspectable* physics
   quantity; the value is verified against the structurally-independent
   Branch-1 symbolic oracle; and the drift is FP-non-associativity bounded by
   the reduction depth.  The LD gates carry ``rtol = 1e-12`` (far above the
   ULP-scale drift); the DD-only strict gate remains the **bit-identical
   negative control** (DD never reaches the LD helper — its :math:`w=\tfrac12`
   reconstruction is the exact power-of-two doubling that commutes with
   round-to-nearest).

The verification is :mod:`tests.sn.spatial.test_ld_ubld_primitive` (10 tests):
the numpy primitive :math:`==` the SymPy oracle at :math:`d=1` (matrices +
moments) and exact-on-bilinear at :math:`d=2`; the shared closed form
:math:`==` the dense :math:`d=1` solve in all three views; and the LIVE
production scheme (``update`` / ``cell_kernel_batch`` /
``affine_scan_coefficients``) :math:`==` the dense primitive (the link proof).
The closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s1_ld_ubld_branch2_closeout.md``.
The :math:`d \ge 2` hand-off (the bilinear cell-batch kernel + the
:math:`2^{d-1}`-moment face cochain wiring onto this primitive) is S2, the next
subsection.

.. _ld-ubld-d2-wavefront-wiring:

Wiring the d≥2 UBLD kernel onto the DAG wavefront (S2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sub-step **D5b-S2** closes the :math:`d = 1`-only kernel raise so
Linear-Discontinuous runs in :math:`d \ge 2` on the DAG wavefront,
consuming the verified dense primitive.  Three contract widenings, all
GATED on a single scheme trait so Diamond Difference / Step stay
byte-identical:

.. math::
   :label: ld-ubld-n-spatial-moments

   \text{per-cell} = (\text{per\_axis})^{d}, \qquad
   \text{per-face} = (\text{per\_axis})^{d-1}, \qquad
   \text{per\_axis} =
   \begin{cases} 1 & \text{DD / Step} \\ 2 & \text{LD} \end{cases}

The class-level trait
:attr:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.spatial_basis_per_axis`
(the 1-D moment-basis size) indexes the whole contract via the
tensor-product structure: the per-cell unknown is
:math:`(\text{per\_axis})^d` (LD-2D: 4) and each downstream face carries
:math:`(\text{per\_axis})^{d-1}` transverse moments (LD-2D: 2).  The
boolean ``per_axis > 1`` gates the multi-moment face-cochain trailing
axis and the moment-reducing emit; DD/Step at ``per_axis == 1`` keep the
rank-:math:`r` scalar face and rank-3 ``psi_avg`` EXACTLY.

.. math::
   :label: ld-ubld-divv-scale-free-kernel

   A_{\div V}\,\vec\psi = M_{\div V}\,\vec S + \sum_a F_{\rm in}^{(a)},
   \qquad
   \psi_{\rm out}^{(a)}[t] = \psi[o_a{=}0,\,t] + \psi[o_a{=}1,\,t]

The :math:`d \ge 2` arm rides the **scale-free ÷V** form of the dense
system: dividing the Galerkin balance by the cell volume leaves a system
depending only on the per-axis ÷V streaming :math:`g_a = |\mu_a|/\Delta_a`
(the ``s_axes`` the kernel already receives) and :math:`\Sigma_t` — so the
dense assembler is fed unit widths and ``mus = (g_0, \ldots)``, reducing
EXACTLY to ``d1_closed_form`` at :math:`d=1`.  The :math:`d` downstream
faces are the trace of the tensor-Legendre solution at the downstream node
(:math:`P_0(+1) = P_1(+1) = 1` sums the :math:`o_a{=}0` and
:math:`o_a{=}1` blocks), in the :math:`2^{d-1}` transverse-Kronecker order
the next cell's upwind inflow consumes (out-of-cell == in-of-next-cell —
the closure consistency the matvec twin verifies).

The scale-free ÷V system fed to the dense primitive
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :math:`d \ge 2` arm rides the **scale-free ÷V** form of the dense system
:eq:`ld-ubld-divv-scale-free-kernel`.  Dividing the Galerkin balance by the
cell volume leaves a system depending only on the per-axis ÷V streaming
:math:`g_a = |\mu_a|/\Delta_a` (the ``s_axes`` the kernel already receives) and
:math:`\Sigma_t`.  In code this means the dense assembler ``_ubld_system`` /
``per_cell_solve`` is fed **unit widths** (``hs = [1, …]``) and
``mus = (g_0, …, g_{d-1})`` — so at :math:`d = 1` it reduces EXACTLY to
``d1_closed_form`` (the ÷V view of the Rule-of-Three above), and at
:math:`d \ge 2` it is the same dense object the Branch-1 oracle proves
exact-on-bilinear.  The kernel dispatch lives in
:mod:`orpheus.sn.spatial.linear_discontinuous` (``cell_kernel_batch`` /
``residual_kernel_batch``): the :math:`d=1` closed-form fast path vs the
:math:`d \ge 2` dense ``_ubld_system`` / ``per_cell_solve``.

The downstream faces are the trace of the tensor-Legendre cell solution at the
downstream node: since :math:`P_0(+1) = P_1(+1) = 1`, the outgoing face on axis
:math:`a` sums the :math:`o_a = 0` and :math:`o_a = 1` blocks of the cell moment
vector, producing a :math:`2^{d-1}`-moment face object (average + transverse
slopes) in the transverse-Kronecker order the next cell's upwind inflow
consumes.  *Out-of-cell == in-of-next-cell* is the closure consistency the
matvec twin verifies.

The moment-ordering crosswalk
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The cell moment vector is the tensor (Kronecker) product of the per-axis 1-D
Legendre basis, ordered **x-outer / y-inner** so the all-:math:`P_0` cell
average is always slot 0 (the same convention the :ref:`spatial-moment-space`
factor surfaces, :eq:`spatial-moment-kronecker-order`).  The Kronecker layout in
2-D is :math:`[\bar\psi,\ \hat\psi_y,\ \hat\psi_x,\ \hat\psi_{xy}]` (indexing
:math:`[o_x, o_y]` with :math:`o_x` outer); each downstream face carries its
:math:`2^{d-1}` transverse moments in the matching per-axis order.  The
crosswalk between the cell-moment order, the per-face transverse order, and the
downstream-node trace reconstruction is the design record's load-bearing detail
(``.claude/plans/issue_240_d5b_s2_crosswalk.md``; recovery anchor
``.claude/plans/issue_240_phase2_step_d5b_ubld.md``).

The DD bit-identity backward-compat invariant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All three contract widenings — the dense kernel arm, the multi-moment
face-cochain trailing axis (:mod:`orpheus.sn.sweep_graph` ``_MovingFrontier``;
the ``_CellSolve`` / ``_CellResidual`` moment-reducing emit), and the window
zero-pad (:mod:`orpheus.sn.loss_representation`
``FullFieldWavefront._octant_face_cochain``, the ``_inflow_to_moments`` pad) —
are GATED on the single scheme trait
:attr:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.spatial_basis_per_axis`
of :eq:`ld-ubld-n-spatial-moments`.  At ``per_axis == 1`` (DD / Step) the tail
is the EMPTY tuple (:eq:`spatial-moment-append-policy`), so the scalar face and
the rank-3 ``psi_avg`` emit are kept byte-identical — DD backward-compatibility
falls out of the ``face_moment_tail`` formula, NOT an ``if scheme is DD`` branch.
This is the negative control: the DD/Step gate stays at its exact pre-S2 golden.

S2 scope boundaries
^^^^^^^^^^^^^^^^^^^

S2 wires the :math:`d \ge 2` UBLD kernel so 2-D LD *runs* on the DAG wavefront,
but it is deliberately scoped to the **average-moment iterate** only.  Three
things remain owed to the later sub-steps, and naming them here is the honest
boundary:

* The full ``loss_action`` / Krylov 2-D LD needs the spatial-moment iterate
  :math:`\hat\phi` to travel between sweeps so the scattering slope source
  :math:`\Sigma_s\hat\phi` couples the slopes globally — **S3**
  (:ref:`ld-ubld-unified-moment-matvec`).  Without it the S2 closure is
  :math:`O(h^2)` but diffusion-limit-inconsistent (the flat-source signature).

* The non-vanishing domain-inflow moment trace (the ``BoundaryFlux`` /
  ``mesh.trace`` widening to a :math:`2^{d-1}` transverse face moment) — **S4**
  (and its honest-scope caveat, :ref:`ld-cartesian-2d`).

* The strengthened vv Mode-7 stress-ansatz MMS and the thick-diffusive
  tripwire — **S4** and **S3** respectively.

The verification is the kernel round-trip + matvec-twin face reconstruction
(:mod:`tests.sn.spatial.test_linear_discontinuous` ``TestLDKernel``), the
end-to-end two-paths FFW :math:`\equiv` MFW, the DD :math:`\ne` LD routing-flip,
and the :math:`O(h^2)` convergence smoke
(:mod:`tests.sn.verification.mms.test_mms_ld_2d`), plus the :math:`d=2`
numpy↔symbolic entry-wise ``A == A`` cell-assembly pin and the
``test_d2_exact_on_bilinear`` ERR-060 catcher
(:mod:`tests.sn.spatial.test_ld_ubld_primitive`).

.. _two-moment-axes:

Two kinds of "moment": angular vs spatial
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   The word **moment** denotes two ORTHOGONAL things in this solver, and
   the collision is the single most common source of confusion when reading
   the multi-dimensional Linear-Discontinuous (LD) code.  An **angular
   moment** reduces the *direction* dependence; a **spatial moment** resolves
   the *within-cell position* dependence.  They are independent tensor
   factors of the flux, NOT two names for the same axis.

The discrete-ordinates flux is a function of three independent kinds of
variable: direction :math:`\Omega`, position :math:`\vec r` (which the mesh
splits into a cell index plus a *within-cell* coordinate), and energy
group :math:`g`.  Each admits its own moment expansion, and the LD scheme
in :math:`d \ge 2` carries two of them simultaneously.  Distinguishing them
is the prerequisite for reading the next two subsections.

**Angular moment** :math:`\phi_\ell^m` — *how the flux varies with
direction.*  Projecting the per-ordinate angular flux
:math:`\psi(\Omega_n)` onto the real spherical harmonics
:math:`\{Y_\ell^m\}` collapses the :math:`N` discrete directions into the
:math:`(\ell, m)` harmonic coefficients,

.. math::
   :label: two-moment-angular

   \phi_\ell^m(\vec r, g)
   \;=\;
   \sum_{n=1}^{N} w_n\, Y_\ell^m(\Omega_n)\, \psi_n(\vec r, g),

the typed home of which is
:class:`~orpheus.transport.fields.harmonic_moment_field.HarmonicMomentField`
on the space
:math:`\mathrm{SphericalHarmonicSpace}(L) \otimes \mathrm{CellGroupSpace}`.
The angular moment is a **replacement representation** of the angular flux:
a calculation holds EITHER the per-ordinate field
:math:`\psi` of shape :math:`(N, ng, *\text{spatial})` OR its harmonic
moments :math:`\phi_\ell^m` of shape
:math:`(L{+}1, 2L{+}1, ng, *\text{spatial})`, bridged by the Galerkin pair
:class:`~orpheus.numerics.projection.MomentProjection` (:math:`M`, the
:math:`\psi \to \phi` reduction :eq:`two-moment-angular`) and
:class:`~orpheus.numerics.projection.HarmonicMomentReconstruction`
(:math:`R`, the :math:`\phi \to \psi` lift).  You never carry both;
windowing the 2-D Cartesian iterate as :math:`\phi_\ell^m` instead of
:math:`\psi` is the harmonic-moment-projection memory win
(the :math:`N \to (L{+}1)(2L{+}1)` collapse, :eq:`harmonic-moment-projection`).
The :math:`\ell = 0` moment IS the scalar flux exactly.

**Spatial moment** :math:`\hat\psi` — *how the flux varies in space inside
one cell.*  A finite-volume / Diamond-Difference closure carries a single
number per cell (the cell average :math:`\bar\psi`).  The
Linear-Discontinuous closure additionally resolves the **sub-cell slope**:
on a Cartesian cell it expands :math:`\psi` in the tensor product of a 1-D
Legendre basis :math:`\{1, P_1\}` per axis,

.. math::
   :label: two-moment-spatial

   \psi(x, y)\big|_{\rm cell}
   \;=\;
   \bar\psi\,
   + \hat\psi_x\, P_1(\xi_x)
   + \hat\psi_y\, P_1(\xi_y)
   + \hat\psi_{xy}\, P_1(\xi_x) P_1(\xi_y),
   \qquad \xi_a \in [-1, 1],

the four within-cell coefficients of the bilinear (UBLD) basis
:math:`\{1, x, y, xy\}` of :eq:`ld-ubld-cell-system`.  Unlike the angular
moment, the spatial moment is an **additional axis** that rides on whatever
angular representation is in play — it does NOT replace anything.  Its typed
home is the :ref:`spatial-moment-space`
(:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`),
a tensor factor of length :math:`(\text{per\_axis})^d` that composes onto a
field's space alongside the cell/group/angular factors.

The two notions are summarised in the contrast table:

.. list-table:: The two "moment" axes — orthogonal tensor factors
   :header-rows: 1
   :widths: 18 26 28 28

   * - Property
     - Angular moment :math:`\phi_\ell^m`
     - Spatial moment :math:`\hat\psi`
     - Shared
   * - Resolves
     - direction :math:`\Omega`
     - within-cell position :math:`x`
     - —
   * - Basis
     - real spherical harmonics :math:`\{Y_\ell^m\}`
     - tensor-Legendre :math:`\{1, P_1\}` per axis
     - both orthogonal polynomial families
   * - Truncation knob
     - :math:`L` (max Legendre order)
     - ``per_axis`` (1 = DD/Step, 2 = LD)
     - —
   * - Count
     - :math:`(L{+}1)(2L{+}1)`
     - :math:`(\text{per\_axis})^d`
     - —
   * - Typed space
     - :class:`SphericalHarmonicSpace`
     - :class:`SpatialMomentSpace`
     - both :class:`FunctionSpace` factors
   * - Role on the flux
     - **replacement** (hold :math:`\psi` OR :math:`\phi_\ell^m`)
     - **additional** (rides on either)
     - both tensor factors of one field
   * - Set by
     - the angular Pℓ order requested
     - the spatial discretization scheme
     - —

Because they are independent factors, a fully-resolved LD-Pℓ angular flux
carries BOTH indices simultaneously — an angular index :math:`(\ell, m)`
and a spatial-moment index :math:`p`:

.. math::
   :label: two-moment-tensor-product

   \phi_\ell^{m, p}(\vec r_{\rm cell}, g),
   \qquad
   (\ell, m) \in \text{angular harmonics}, \quad
   p \in \{\bar{\,}, \hat x, \hat y, \widehat{xy}\}
   \ \text{(spatial moments)}.

The carrier space is the tensor product of the two moment spaces with the
cell/group space,

.. math::
   :label: two-moment-carrier-space

   \mathrm{SphericalHarmonicSpace}(L)
   \;\otimes\;
   \mathrm{CellGroupSpace}(ng, *\text{spatial})
   \;\otimes\;
   \mathrm{SpatialMomentSpace}(\text{per\_axis}, d),

so the stored ndarray gains a trailing :math:`(\text{per\_axis})^d`
spatial-moment axis after the :math:`(\ell, m, g, *\text{spatial})` prefix.
The orthogonality is what makes the architecture clean: the scattering
operator :math:`\Sigma_s` couples energy groups and (for anisotropic
scattering) angular moments, but is a **spectator** to the spatial-moment
axis — it scatters every spatial moment independently
(:eq:`ld-ubld-scattering-moment-lift`, next subsection).  Conversely the
spatial discretization (the sweep / cell solve) acts on the spatial moments
but is a spectator to the angular index.  Two operators, two axes, no
cross-talk except through the physics each is responsible for.

.. note::

   **Why an LD-P3 calculation needs both.**  Anisotropic scattering up to
   :math:`P_3` is an *angular*-resolution choice — it carries
   :math:`\phi_\ell^m` for :math:`\ell \le 3`.  The Linear-Discontinuous
   spatial closure is a *spatial*-resolution choice — it carries
   :math:`\hat\psi` for the within-cell slope.  An LD-P3 transport
   calculation makes both choices at once and so carries the full
   :math:`\phi_\ell^{m, p}` object of :eq:`two-moment-tensor-product`.
   Collapsing either axis to its average (P0 angular, or DD spatial)
   degrades a *different* physical fidelity: the angular collapse loses the
   flux's directional anisotropy; the spatial collapse loses the
   diffusion-limit accuracy that the :math:`xy` cross-moment provides
   (:eq:`ld-ubld-cell-system`, the load-bearing moment).

.. _ld-ubld-scattering-moment-lift:

The Σ_s ⊗ I spatial-moment scattering lift (S3-A, partial)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sub-step **D5b-S3** completes the *physics* of the multi-dimensional UBLD
Linear-Discontinuous scheme.  Now that the two moment axes are clearly
distinguished (:ref:`two-moment-axes`), the completion is statable in one
sentence: the scattering source must scatter EVERY spatial moment, not just
the cell average.  Where S2 ships an O(h²) but diffusion-limit-INCONSISTENT
closure (it scatters only the spatial-AVERAGE moment — the slope rows of the
source are zero), S3 threads the canonical slope source so the converged
operator becomes the diffusion-limit-CONSISTENT one.

The load-bearing bridge is the scattering operator's
:math:`\Sigma_s \otimes I_{\rm spatial}` lift: :math:`\Sigma_s` carries no
spatial-moment index (it is an energy-group :math:`\to` energy-group matrix
per Legendre order), so it is applied to EVERY spatial moment of the scalar
flux INDEPENDENTLY,

.. math::
   :label: ld-ubld-scattering-moment-lift

   \bigl(\Sigma_s \otimes I_{\rm spatial}\bigr)\,
   (\bar\phi,\ \hat\phi_x,\ \hat\phi_y,\ \hat\phi_{xy})
   \;=\;
   (\Sigma_s\,\bar\phi,\ \Sigma_s\,\hat\phi_x,\
    \Sigma_s\,\hat\phi_y,\ \Sigma_s\,\hat\phi_{xy}),

so the spatial-moment axis is a SPECTATOR to the energy-group matmul,
exactly as the cell axis is.  In code this is the per-material group
contraction with a trailing ``...`` spectator
(``"fg,fc...->gc..."``): at the single-moment closures (Diamond
Difference / Step, ``per_axis == 1``) the trailing axis is ABSENT and the
``...`` matches nothing, so the lift is BYTE-IDENTICAL to the pre-S3
scattering (the negative-control bit-identity, verified rank-2-exact).
At ``per_axis == 1`` :math:`S_{\rm full} \equiv S_{\rm flat}`; only an LD
multi-moment closure activates the slope rows.

.. admonition:: Status — S3-A is PARTIAL
   :class: caution

   The :math:`\Sigma_s \otimes I_{\rm spatial}` lift documented here is the
   LANDED half of S3-A.  The :math:`\hat\phi` spatial-moment **iterate
   carrier** that FILLS the slope rows the lift now accepts is OWED (it was
   blocked on the typed-field space widening — the
   :ref:`spatial-moment-space` subsection, the prerequisite that was minted
   next).  The lift therefore scatters a slope source that, in the production
   path, is still zero (no field carries :math:`\hat\phi` yet); the converged
   fixed point does not change UNTIL the iterate carrier lands.  This page
   marks what is wired (the lift) versus what is owed (the iterate, the
   cell-emit accumulation, the source seams) so a future reader knows the
   S3-A wiring is mid-flight, not complete.

Physics-completion, not an iteration-only change
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The distinction matters for verification, so it is worth stating precisely.
Most changes to the iteration machinery — a Gauss-Seidel splitting, a σ\
:sub:`r`-removal, a synthetic accelerator (DSA), a preconditioner — MUST NOT
change the converged fixed point; they change only the *rate* at which the
iteration reaches it.  The correctness gate for such a change is
**FP-invariance**: the accelerated solve and the plain solve converge to the
same flux (`vv-principles` Mode 9).

The :math:`\Sigma_s \otimes I_{\rm spatial}` slope source is **not** that
kind of change.  S2 and S3 solve DIFFERENT operators:

.. math::
   :label: ld-ubld-s2-s3-operators

   \text{S2:}\quad (L + C - S_{\rm flat})\,\psi = Q_{\rm ext},
   \qquad
   S_{\rm flat} = \Sigma_s \otimes e_0 e_0^{\mathsf T},
   \\[4pt]
   \text{S3:}\quad (L + C - S_{\rm full})\,\psi = Q_{\rm ext},
   \qquad
   S_{\rm full} = \Sigma_s \otimes I_{\rm spatial}.

:math:`S_{\rm flat}` (the rank-1 projector :math:`e_0 e_0^{\mathsf T}` onto
the cell-average moment) scatters ONLY the spatial average — the slope rows
:math:`\hat\phi_x, \hat\phi_y, \hat\phi_{xy}` of the scattering source are
identically zero.  :math:`S_{\rm full}` (the identity on the spatial-moment
axis) scatters all of them.  The two operators have DIFFERENT spectra, hence
DIFFERENT fixed points.  The converged flux CHANGES — and that is the POINT:
the thick-diffusion-limit tripwire (the ``test_ld_thick_diffusive_limit``
xfail) flips xfail :math:`\to` pass *because* the limit becomes correct, not
because the iteration was accelerated.  S3 is therefore **NOT** verified
against the S2 fixed point; verifying it that way would be the Mode-9
mis-application (asserting FP-invariance of a change that legitimately moves
the FP).  The genuine Mode-9 invariant for S3 is the within-group analog:
source-iteration with a lagged moment iterate :math:`\equiv` direct / Krylov
solve of the **same** :math:`(L + C - S_{\rm full})` operator.

Why the slope rows are diffusion-limit-load-bearing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the optically-thick, scattering-dominated (:math:`c \to 1`,
:math:`\Sigma_t h \gg 1`) limit, the transport solution must collapse to the
diffusion solution.  Adams (2001) and Larsen, Morel & Miller (1987) showed
that a spatial discretization passes this asymptotic limit only if its
discrete diffusion limit is a valid (consistent and stable) diffusion
discretization.  For the bilinear UBLD cell the diffusion limit couples the
slope moments :math:`\hat\phi` — the leading-order asymptotic balance is a
relation between the cell-average and the slopes, not the average alone
(Border, Lewis & Adams 1992 give the 2-D asymptotic analysis explicitly).
If the *scattering source* feeds only the cell average (S2's
:math:`S_{\rm flat}`), the slope rows of the within-cell balance see no
scattering re-supply, the discrete diffusion limit is the WRONG diffusion
operator, and the thick-limit error does not vanish under refinement — the
xfail tripwire stays red.  Threading :math:`\Sigma_s \hat\phi` into the slope
rows (:math:`S_{\rm full}`) restores the correct discrete diffusion limit.
This is why the completion is *physics* (the converged answer becomes right
in a regime where it was wrong), not iteration bookkeeping.

.. note::

   This is the same asymptotic reasoning that selects the bilinear
   :math:`\{1, x, y, xy\}` basis over the simplex :math:`\{1, x, y\}` in the
   first place (:eq:`ld-ubld-cell-system` and the parent section): Adams
   (2001) proved the simplex-LD discrete diffusion limit is invalid on
   quadrilaterals.  The :math:`xy` cross-moment carries the limit; the
   :math:`\Sigma_s \otimes I_{\rm spatial}` lift makes sure that cross-moment
   (and the axis slopes) actually receive scattering.  The basis choice and
   the scattering lift are two halves of the *same* diffusion-limit argument.

The producer-side spectator-broadcast (Pattern 7)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The lift is implemented as a one-character change to the einsum subscripts of
the three scattering producers in
:class:`~orpheus.sn.material_xs_field.MaterialXSField`:

.. list-table:: The :math:`\Sigma_s \otimes I_{\rm spatial}` lift in code
   :header-rows: 1
   :widths: 34 30 36

   * - Producer
     - Subscript (pre-S3 :math:`\to` S3)
     - What it scatters
   * - :meth:`~orpheus.sn.material_xs_field.MaterialXSField.apply_p0_in_scatter`
     - ``"fg,fc->gc"`` :math:`\to` ``"fg,fc...->gc..."``
     - the P0 in-scatter :math:`\Sigma_{s,0}^{\mathsf T}\phi`
   * - :meth:`~orpheus.sn.material_xs_field.MaterialXSField.apply_n2n`
     - ``"fg,fc->gc"`` :math:`\to` ``"fg,fc...->gc..."``
     - the :math:`(n,2n)` source :math:`2\Sigma_{2n}^{\mathsf T}\phi`
   * - :meth:`~orpheus.sn.material_xs_field.MaterialXSField.apply_legendre_scattering_moments`
     - ``"mfc,fg->mgc"`` :math:`\to` ``"mfc...,fg->mgc..."``
     - the per-:math:`\ell` block-diagonal :math:`\Lambda\phi`

The trailing ``...`` is the **spectator broadcast**: it matches the
spatial-moment axis (if present) and contracts nothing over it — exactly the
:math:`\otimes I_{\rm spatial}` of :eq:`ld-ubld-scattering-moment-lift`.  This
is the `coding-elegance` Pattern-7 producer-side lift: the convention
(scatter each spatial moment independently) is normalised at the producer, so
no consumer special-cases the axis.  Two properties follow by construction:

* **Byte-identical at the single-moment closures.**  When ``phi`` is rank-2
  (``(ng, n_cells)``, the DD/Step shape) the trailing axis is ABSENT, ``...``
  matches nothing, and ``"fg,fc...->gc..."`` is the SAME contraction as
  ``"fg,fc->gc"`` — verified rank-2-exact
  (``np.array_equal`` of the two einsums when no trailing axis is present).
  No re-baseline of the DD/Step path; this is the negative-control
  bit-identity.

* **The projection pair needed no change.**  The Galerkin projection
  primitives :class:`~orpheus.numerics.projection.MomentProjection` and
  :class:`~orpheus.numerics.projection.HarmonicMomentReconstruction` already
  carry ``...`` for their trailing axes, so :math:`M` and :math:`R` are
  spatial-moment-agnostic out of the box.  The angular reduction
  :eq:`two-moment-angular` and its inverse ride a spatial-moment axis as a
  spectator, which is the architectural payoff of the orthogonal-factor
  framing — the two moment axes never need to know about each other.

.. warning::

   The crosswalk and the original brief ASSUMED ``apply_p0_in_scatter``
   already broadcast over a trailing axis.  It did NOT: the bare ``"fg,fc->gc"``
   hard-codes the cell axis as a single index ``c``, so a rank-3
   ``phi (ng, n_cells, 2^d)`` RAISES (``operand has more dimensions than
   subscripts``).  The fix is the explicit ``...`` spectator, not a reshape.
   Documented so a future reader does not re-derive the (false) assumption.

What is still owed (the iterate carrier and the source seams)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The lift accepts a slope source, but nothing in the production path yet
PRODUCES one.  Filling the slope rows requires (all S3-A proper, owed):

* The :math:`\hat\phi` **iterate carrier** — the between-sweep flux must
  carry the spatial-moment axis.  This was the make-or-break design decision:
  the typed-field spaces validate ``shape == (ng, *spatial)`` with no slot for
  a trailing :math:`(\text{per\_axis})^d` axis, so a slope-carrying field was
  an *illegal state* (Pattern 4 firing correctly).  The resolution — minting
  the first-class :ref:`spatial-moment-space` factor — is the subject of the
  next subsection.

* The **cell-emit moment accumulation** — the wavefront cell solve already
  computes a :math:`(\text{per\_axis})^d`-moment ``psi_avg``, but the
  between-sweep emit currently drops to slot 0 (the cell average); it must
  accumulate the full moment vector.

* The **two source seams** — the :math:`d \ge 2` wavefront genuine
  :math:`(2^d, ng)` moment source through the dense ``_ubld_system``, and the
  :math:`d = 1` scan slope source threaded via
  :meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.kernel_rhs` and the Schur
  ``schur_xV`` term.

The verification chain for the completed S3 is the thick-diffusion-limit
VALUE anchor (the continuous diffusion solution at :math:`\varepsilon \to 0`,
structurally independent of the LD kernel — Adams 2001 / Border-Lewis-Adams
1992 / Larsen-Morel-Miller 1987), the convergence-order MMS smoke (the slope
source exercised), and the genuine Mode-9 SI :math:`\equiv` Krylov on
:math:`(L + C - S_{\rm full})`.  The design record is
``.claude/plans/issue_240_d5b_s3_crosswalk.md``; the verification spec is
``.claude/agent-memory/test-architect/d5b_s3_inc_c_moment_iterate_verification.md``;
the lift's landed-half closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_a_inc_c_closeout.md``.

.. _spatial-moment-space:

The SpatialMomentSpace: a first-class within-cell DG moment carrier (S3-A0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The typed-field-space half of S3-A (the half the scattering-lift TODO above
flagged as a hard prerequisite).  The within-cell tensor-Legendre DG moment
axis — how :math:`\psi` varies in space WITHIN a cell — is minted as a
first-class function space,
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`,
the **spatial** sibling of the **angular**
:class:`~orpheus.numerics.spaces.spherical_harmonic_space.SphericalHarmonicSpace`.
The two "moment" notions are ORTHOGONAL axes (angular harmonics over
direction :math:`\Omega` vs spatial Legendre over within-cell position
:math:`x`); naming each as its own typed factor keeps the distinction
type-visible and dispels the collision.

.. math::
   :label: spatial-moment-space-size

   \dim(\text{SpatialMomentSpace}) \;=\; (\text{per\_axis})^{d},
   \qquad
   \text{per\_axis} =
   \begin{cases}
     1 & \text{DD / Step (cell-average } \{1\}\text{)} \\
     2 & \text{LD (linear } \{1, P_1\}\text{)}
   \end{cases}

The factor composes via the tensor product ``*`` into the bulk-field spaces
EXACTLY as the angular factor does, and is recovered by type via
``space.find_factor(SpatialMomentSpace).per_axis`` (#207).  The
field-space factories
(:meth:`~orpheus.transport.fields._bases.AngularField.from_mesh`,
:meth:`~orpheus.transport.fields._bases.ScalarField.from_mesh`,
:meth:`~orpheus.transport.fields.harmonic_moment_field.HarmonicMomentField.from_mesh_and_L`)
gained an OPTIONAL ``spatial_moments`` parameter (default ``1``) that
appends the factor **iff the within-cell count exceeds 1** — the
"append iff > 1" gate single-sourced from
:func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`
(the cell analogue that delegates to
:func:`orpheus.sn.spatial._ubld.face_moment_tail`, so the cell-moment tail
and the per-face cochain tail can never disagree).  At the default the
field space is BYTE-IDENTICAL to its pre-S3 shape for EVERY scheme (DD,
Step, AND LD): this step builds the CAPABILITY only (construct-general /
select-narrow), and no production field selects the axis yet.

Why a first-class typed factor, not a bare int axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The slope-carrying flux could, in principle, be stored as a plain ndarray
with a trailing :math:`(\text{per\_axis})^d` axis and an integer remembered
somewhere for its width.  That was rejected (the user's design choice
"option b") in favour of a first-class typed
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace`
factor, for the `coding-elegance` Pattern-4 reason
(*make illegal states unrepresentable*).

The typed-field layer validates ``values.shape == space.shape`` at
construction (:meth:`Field.__post_init__`).  Before this step the SN field
spaces were rigidly ``(ng, *spatial)`` / ``(N, ng, *spatial)`` /
``(L+1, 2L+1, ng, *spatial)`` — there was **no slot** for a trailing
spatial-moment axis, so a :math:`\hat\phi`-carrying field FAILED the gate.
A trailing slope axis was, literally, an illegal state (the gate was firing
*correctly* — see the scattering-lift status admonition above, which flagged
this as the make-or-break prerequisite).  Two ways to make the slope-carrying
field legal:

.. list-table:: Bare-int axis vs first-class typed factor
   :header-rows: 1
   :widths: 22 39 39

   * - Aspect
     - Bare ``int`` trailing axis
     - Typed :class:`SpatialMomentSpace` factor
   * - Field validity
     - widen the shape gate to accept *any* trailing axis (loses the
       illegal-state guard)
     - the space DECLARES the axis; the gate stays exact — a slope field
       is now a *legal, declared* shape
   * - Querying the width
     - thread an ``int`` parameter through every call site, or re-derive
       it from a raw ``.shape[-1]``
     - ``space.find_factor(SpatialMomentSpace).per_axis`` — query by
       TYPE, position-independent (#207)
   * - Self-description
     - the axis is anonymous; reading code cannot tell a spatial-moment
       axis from any other trailing axis
     - the factor's type IS its documentation; the
       angular/spatial collision is dispelled at the type level
   * - Precedent
     - none — a one-off convention
     - the EXACT mold of the angular
       :class:`SphericalHarmonicSpace` factor (one architecture, two axes)

The typed factor is the same pattern the angular moment already uses: the
harmonic factor is a :class:`SphericalHarmonicSpace` whose ``L`` is recovered
by ``space.find_factor(SphericalHarmonicSpace).L``, NOT a bare integer
threaded through the API.  Minting the spatial sibling keeps the two axes
*symmetric* — the orthogonal-factor framing of :ref:`two-moment-axes` is then
literally how the carrier space is built (:eq:`two-moment-carrier-space`).

.. note::

   Closing #207 as a side effect.  The
   ``space.find_factor(SphericalHarmonicSpace).L`` query was DOCUMENTED in the
   :class:`HarmonicMomentField` docstrings (issue #207) but had never been
   IMPLEMENTED — three docstrings referenced a method that did not exist.  The
   spatial-moment work needed exactly this composition-tree query, so
   :meth:`~orpheus.numerics.space.TensorProductSpace.find_factor` was minted
   now: it returns the first tensor factor that ``isinstance(factor, T)`` and
   raises :exc:`KeyError` if absent (a structural assertion — the caller
   believes the composed space carries the factor — not a silent ``None``,
   Pattern 4).  Both moment factors (angular and spatial) are now queryable
   by type, and the latent broken claim in the docstrings is made true.

The Kronecker moment ordering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The within-cell basis is the tensor (Kronecker) product of the per-axis 1-D
Legendre basis, ordered **x-outer / y-inner** (matching the UBLD assembler
:func:`orpheus.sn.spatial._ubld.assemble_ubld`).  The convention is fixed so
that the all-:math:`P_0` cell average is ALWAYS slot 0:

.. math::
   :label: spatial-moment-kronecker-order

   d{=}1:&\quad [\,\bar\psi,\ \hat\psi_x\,]
   \\[2pt]
   d{=}2:&\quad [\,\bar\psi,\ \hat\psi_y,\ \hat\psi_x,\ \hat\psi_{xy}\,]
   \\[2pt]
   d{=}3:&\quad [\,\bar\psi,\ \hat\psi_z,\ \hat\psi_y,\ \hat\psi_{yz},\
                  \hat\psi_x,\ \hat\psi_{xz},\ \hat\psi_{xy},\ \hat\psi_{xyz}\,]

The slot-0 (cell-average) convention is single-sourced from
:data:`orpheus.sn.spatial._ubld.AVERAGE_MOMENT` (the constant every moment
consumer reduces on) — the :class:`SpatialMomentSpace` surfaces it via
:attr:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace.average_moment_index`
rather than re-spelling the literal ``0`` (so a layout change happens in ONE
place, not at the scattered ``[..., 0]`` call sites).  Slot 0 is the link
between the two scales: the cell-average moment :math:`\bar\psi` is what the
DD/Step closure carries in full, and it is the moment the
:math:`\Sigma_s \otimes I_{\rm spatial}` lift scatters at every closure; the
remaining slots are the LD-only slope rows the lift activates.

The "append iff > 1" byte-identity policy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The backward-compatibility invariant (#240 D5b) is that DD/Step field shapes
must be UNCHANGED.  This is enforced by a single policy, single-sourced so the
cell-moment tail and the per-face cochain tail can never drift apart:

.. math::
   :label: spatial-moment-append-policy

   \texttt{tail}(n) =
   \begin{cases}
     ()    & n = 1 \quad\text{(DD/Step — NO length-1 axis appended)} \\
     (n,)  & n > 1 \quad\text{(LD — a genuine trailing moment axis)}
   \end{cases}

The critical detail is that ``n == 1`` returns the EMPTY tuple, NOT
``(1,)`` — a length-1 axis is NOT appended.  Appending ``(1,)`` would
broadcast-equal the old shape numerically but would change ``ndarray.shape``
and ``ndim``, breaking every byte-identity gate and every consumer that reads
``.ndim``.  The empty-tuple branch keeps the DD/Step field space *literally
identical* to its pre-S3 self.  :func:`orpheus.sn.spatial._ubld.face_moment_tail`
owns the policy; the cell analogue
:func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`
delegates to it (Pattern 7 — normalise the convention at one site), and
:meth:`BulkField._compose_spatial_moments`
(:mod:`orpheus.transport.fields._bases`) returns the space UNCHANGED when the
tail is ``()``.

Construct-general, select-narrow — what this step does and does NOT do
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This step builds the **capability** to carry the spatial-moment axis and
nothing more.  The discipline is deliberate (`coding-elegance` — construct
general, select narrow, specialize only on measured need):

* **The axis exists.**  The :class:`SpatialMomentSpace` is minted, composes
  into every bulk-field space, and is ``find_factor``-queryable.

* **No production field selects it.**  The ``spatial_moments`` factory
  parameter defaults to ``1`` at EVERY call site and is NOT auto-read from
  ``mesh.scheme.spatial_basis_per_axis``.  So DD, Step, AND LD field shapes
  are unchanged this step — even LD does not yet carry the slope axis.

* **Why default-OFF even for LD.**  Auto-reading the scheme would silently
  widen LD field shapes BEFORE the consumers that FILL the axis exist — the
  iterate carrier, the cell-emit accumulation, the source seams (all S3-A
  proper, owed; see the scattering-lift subsection).  A widened axis that no
  producer fills is precisely the illegal state Pattern 4 forbids; turning the
  capability on before its producers exist would re-introduce it.  The gate
  has teeth on exactly this mistake: making
  :meth:`BulkField._compose_spatial_moments` auto-read the scheme turns the LD
  byte-identity foundation tests RED (mutation-verified).

The S3-A iterate / cell-emit / source seams that thread the scheme's
``spatial_basis_per_axis`` here (selecting the axis for LD) are the NEXT
sub-step.  When they land, the only change at the factory call sites is
passing ``spatial_moments=scheme.spatial_basis_per_axis``; the validator
already accepts the widened space, and the scattering lift already scatters
its slopes.

Verification (foundation-level)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The space and the factory widening are software invariants, so they are
verified at the **foundation** level (data-structure / factory-output
invariants, not an L0/L1/L2 solver claim — they carry no eigenvalue or flux
assertion), in two test modules:

* :mod:`tests.numerics.test_spatial_moment_space` — the space layer: the
  :math:`(\text{per\_axis})^d` size law, the
  :meth:`~orpheus.numerics.space.TensorProductSpace.find_factor` round-trip
  (and the :exc:`KeyError`-when-absent assertion), the composition shape, the
  ``per_axis == 1`` no-widening case, and
  :attr:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace.average_moment_index`
  :math:`==` :data:`~orpheus.sn.spatial._ubld.AVERAGE_MOMENT`.

* ``tests/sn/spatial/test_spatial_moment_field_space.py`` — the factory
  widening: the **byte-identity-at-default negative control** for DD AND LD on
  all three carriers (:class:`AngularField`, :class:`ScalarField`,
  :class:`HarmonicMomentField`), the widened :math:`d{=}1` / :math:`d{=}2`
  shapes, the both-moment-factors-coexist case, and the wrong-shape rejection.
  The mutation check — auto-reading the scheme turns the LD byte-identity
  cases red — is what proves the construct-general gate has teeth.

The design record (the angular-vs-spatial distinction, the FP resolution) is
``.claude/plans/issue_240_d5b_s3_crosswalk.md``; the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_a0_spatial_moment_space_closeout.md``.

.. _ld-ubld-unified-moment-matvec:

The unified moment matvec: a forward apply is intrinsically moment-valued (S3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sub-step **D5b-S3** completes the apply direction: applying the per-cell
:math:`2^d \times 2^d` UBLD operator (:eq:`ld-ubld-cell-system`) to a moment
vector is intrinsically moment-valued, so the matvec carries the full
:math:`(\bar\psi, \hat\psi)` moment vector in every dimension.  The
architectural payoff is a branch removal; the *physics* payoff is the recovery
of the thick-cell diffusion limit, which hinges on a single
frame-consistency identity (ERR-061) that the rest of this subsection derives.
The source files are :mod:`orpheus.sn.spatial.linear_discontinuous`
(``cell_kernel_batch`` / ``residual_kernel_batch`` — now ONE :math:`d`-generic
moment path), :mod:`orpheus.sn.sweep_graph` (``_CellSolve`` / ``_CellResidual``
— the ``len(s_axes) > 1`` moment gate retired), and
:mod:`orpheus.sn.loss_representation` (the ``_spatial_moment_tail`` buffer
widening); the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_unified_matvec_closeout.md``.

The earlier increments (A/B) made the :math:`d{=}1` LD **matvec** Schur-reduce
to a scalar residual: the slope :math:`\hat\psi` was eliminated, leaving a
scalar cell-average unknown.  That was a *flat-source artifact* — with
:math:`\hat Q = 0` the slope had no global coupling, so the Krylov unknown could
stay scalar.  Increment C makes the scattering slope source
:math:`\Sigma_s\hat\phi` couple the slope GLOBALLY (the diffusion-limit-
consistent operator :eq:`ld-ubld-scattering-moment-lift`), so the slope becomes a
genuine global degree of freedom in **every** dimension.

.. math::
   :label: ld-ubld-unified-moment-residual

   (L+C)\,\vec\psi
   \;=\;
   M^{-1}\bigl(A\,\vec\psi - F_{\rm in}\bigr),
   \qquad
   A = (L+C-S)\ \Longleftrightarrow\
   (L+C)\,\vec\psi - S\,\vec\psi = \vec q_{\rm ext}

A matvec is a forward APPLY: applying the per-cell
:math:`2^d \times 2^d` UBLD operator to the moment vector is intrinsically
moment-valued, so ``cell_kernel_batch`` and ``residual_kernel_batch`` collapse to
ONE d-generic dense path for every :math:`d` (the former :math:`d{=}1`
Schur-reduced scalar arm — and the :math:`d \ge 2` raise — are both retired).
The :math:`M^{-1}` factor in :eq:`ld-ubld-unified-moment-residual` is the
matvec/sweep moment-source consistency: the UBLD RHS folds the cell source
mass-weighted (:math:`R = M\vec S`, the test-function projection), but the
operator algebra :math:`A = (L+C) - S` subtracts :math:`S\vec\psi` RAW at the
``OperatorSum`` level, so the residual is divided by the diagonal Legendre mass
to put :math:`(L+C)\vec\psi` in raw per-moment units (the slope rows would
otherwise disagree by :math:`M_{ii} = \theta^{|i|}`).

The architectural headline is **branch removal** (Cardinal Rule 2): the
``len(s_axes) > 1`` moment gate at the cell-solve / cell-residual emit is GONE
(replaced by the pure scheme trait ``spatial_basis_per_axis > 1``), the
:math:`d{=}1` scalar kernel twin is retired, and there are ZERO
``isinstance(scheme, ...)`` branches — dispatch stays via the scheme PROTOCOL +
geometry-keyed ``supports()``.

.. _ld-ubld-sweep-global-frame:

The sweep-frame / global-frame involution (ERR-061 — the diffusion-limit root cause)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is the load-bearing result of the whole multi-dimensional LD campaign, and
it is the kind of "what failed and why" Cardinal Rule 3 demands be archived in
full.  Threading the moment matvec made the operator *internally consistent*
(matvec :math:`\equiv` sweep round-trip to :math:`10^{-16}`, source-iteration
:math:`\equiv` Krylov on the SAME operator) — and yet the converged flux was
**wrong**: on a thick scattering slab (:math:`\Sigma_t = 40`, :math:`c = 0.99`,
:math:`\Sigma_t h = 10`/cell, vacuum) at :math:`n_x = 4` the LD scalar flux was
1.47 against the diffusion solution 2.31 (relative error 36 %), and the error
did not grow under refinement — it *shrank* as the cells thinned (the classic
flat-source-LD signature, persisting THROUGH the slope-source thread).

The cause is a frame-consistency error between two individually-correct
components.  The per-cell LD kernel produces and consumes the :math:`2^d` moment
vector in the per-ordinate **sweep frame**: each axis :math:`a` is oriented so
the *downstream* face is at the local coordinate :math:`+1`.  For an ordinate
sweeping in the NEGATIVE global direction on axis :math:`a`
(:math:`\mathrm{octant\_sign}_a = -1`) the sweep coordinate is the *reverse* of
the global coordinate, so the slope (:math:`P_1`) moment on that axis is
sign-FLIPPED relative to the global-:math:`x` slope.  But the iterate
:math:`\hat\phi` and its scattering source :math:`\Sigma_s\hat\phi` live in the
**global frame** — the angular reduction sums slopes across ordinates of BOTH
sweep directions,

.. math::
   :label: ld-ubld-slope-angular-reduction

   \hat\phi(\vec r, g) \;=\; \sum_{n=1}^{N} w_n\,\hat\psi_n(\vec r, g).

The producer (``_CellSolve.cell`` emit) stored the raw sweep-frame slope; the
consumer (``integrate_angular`` / the scattering apply) summed it as if it were
global-frame.  So the backward ordinates' opposite-signed slopes partially
CANCELLED the forward ones: at a cell with a positive global-:math:`x` gradient
the forward ordinate had :math:`\hat\psi_n = +0.048` but the backward had
:math:`-0.028` — opposite signs, the smoking gun.  The summed
:math:`\hat\phi` was :math:`\sim 6\times` too small to satisfy the LM-1989
discrete diffusion continuity (Larsen & Morel 1989, JCP 83(1):212–236, Eq. 4.9b,
:math:`\bar\phi_j + \hat\phi_j = \bar\phi_{j+1} - \hat\phi_{j+1}`), the slope
source was under-driven, and the discrete diffusion limit was the wrong
diffusion operator.

The fix is a single-sourced :math:`2^d` moment-frame **involution**,
:func:`~orpheus.sn.spatial._ubld.octant_moment_frame_signs`,

.. math::
   :label: ld-ubld-octant-moment-frame-signs

   \mathrm{sign}[o_0, \ldots, o_{d-1}]
   \;=\;
   \prod_{a=0}^{d-1} (\mathrm{octant\_sign}_a)^{\,o_a},
   \qquad o_a \in \{0, 1\},

indexed in the tensor-Legendre Kronecker layout (:math:`o_a` = the :math:`P_0` /
:math:`P_1` selector on axis :math:`a`).  The **average** moment (all
:math:`o_a = 0`) is sign-invariant (the empty product is :math:`1`); a per-axis
**slope** flips once if that axis sweeps backward; the 2-D **cross** moment
:math:`\hat\psi_{xy}` flips when an ODD number of its active axes reverse.  The
map is its own inverse, so the SAME sign vector converts global :math:`\to`
sweep on the source/probe INPUT and sweep :math:`\to` global on the
moment/residual OUTPUT.  It is applied through the shared ``_reframe`` helper at
the cell ops; the OUTGOING FACE (``psi_out``) stays sweep-frame — it propagates
along the wavefront and never crosses into the global-frame iterate (so it is
left untouched).  DD/Step (``spatial_basis_per_axis == 1``) get ``None`` (the
sign-invariant average-only moment), so they pass through ``_reframe``
untouched and stay byte-identical (the negative control); a flat scalar source
(matvec zero / flat external — only the average moment) is frame-invariant and
skipped by the ``arr.shape[-1] != frame_signs.shape[0]`` guard, so it is never
broadcast into a spurious moment axis.

After the fix the diffusion limit is recovered on BOTH the matvec (Krylov) and
the sweep (source-iteration) paths:

.. list-table:: Thick-slab LD vs DD relative error, before/after the frame fix
   :header-rows: 1
   :widths: 14 22 22 22

   * - Mesh
     - Cell optical depth
     - Before (sweep-frame slope)
     - After (global-frame slope)
   * - 1-D, :math:`n_x = 4`
     - :math:`\Sigma_t h = 10`
     - 38.9 %
     - 4.1 %
   * - 1-D, :math:`n_x = 16`
     - :math:`\Sigma_t h = 2.5`
     - 7.9 %
     - 0.2 %
   * - 1-D, :math:`n_x = 64`
     - :math:`\Sigma_t h = 0.6`
     - 0.9 %
     - 0.0 %
   * - 2-D, :math:`n = 4/8/16`
     - thick :math:`\to` thin
     - 8.4 %
     - 1.7 % :math:`\to` 0.4 %

.. warning::

   **The matvec-self-consistency gate is necessary but NEVER sufficient for a
   moment-iterate fold.**  Every component here was individually correct (the
   2×2 matched LM-1989 Eq. 4.3, the dense UBLD matched the analytic 2×2, the
   scattering produced :math:`\Sigma_s\hat\phi` at full strength, the matvec
   round-trip vanished to :math:`10^{-16}`, and source-iteration :math:`\equiv`
   Krylov on the SAME operator).  The bug was the frame consistency *between*
   two correct components — a wrong fixed point that the round-trip and the
   SI :math:`\equiv` Krylov gates are structurally BLIND to: they prove the
   operator is internally consistent, NOT that its fixed point is the
   physically-correct one (`vv-principles` §5 — "O(h²) to the wrong limit is
   still O(h²)").  The decisive evidence was a structurally-independent
   from-scratch LD-SN solver (a direct LM-1989 2×2 + source iteration, no
   ORPHEUS kernel) that reproduced ORPHEUS's WRONG value bit-for-bit when it
   summed sweep-frame slopes and RECOVERED the diffusion limit when it stored
   global-frame slopes — pinning the root cause independent of ORPHEUS's code.
   The lesson: gate the converged VALUE against a structurally-independent
   reference (the continuous diffusion solution + the independent from-scratch
   kernel), never the round-trip.  This is failure Mode 1 (sign flip) +
   Mode 6 (convention drift) — see ``error_catalog.md`` ERR-061.

The thick-cell diffusion tripwire is
``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit``
(1G) and ``::test_ld_thick_diffusive_limit_2g`` (2G heterogeneous, a
group-coupled slope source — Mode 6), both ``@pytest.mark.l1
@pytest.mark.catches("ERR-061")`` and both Mode-8-safe
(``np.testing.assert_array_less`` fires under ``-O``).  The slope-frame
fingerprint is pinned by
``derivations/diagnostics/diag_240_d5b_s3_probe_11_root_cause.py`` (forward and
backward ordinate slopes must share sign in the global frame), and the
structurally-independent confirmation by
``diag_240_d5b_s3_probe_08_independent_ld.py``.

.. _ld-ubld-pure-z-collision-twin:

The pure-z collision-only twin — sweep :math:`\equiv` matvec single source
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. math::
   :label: ld-ubld-pure-z-collision

   \text{(solve)}\quad \bar\psi = \frac{Q}{\sigma_t}
   \qquad\Longleftrightarrow\qquad
   \text{(matvec)}\quad (L+C)\,\bar\psi = \sigma_t\,\bar\psi

The pure-z degenerate ordinates (:math:`\mu_x = \mu_y = 0`, the :math:`\pm z`
poles of a Lebedev or product cubature in a 2-D Cartesian sweep) have no
in-plane streaming, so the cell is **collision-only**: the loss couples to
:math:`\sigma_t` alone, and the sweep balance :math:`\bar\psi = Q/\sigma_t` and
its matvec twin :math:`(L+C)\bar\psi = \sigma_t\bar\psi` are two applications of
the SAME operator (:eq:`ld-ubld-pure-z-collision`).  This is the L21 twin-path
relationship — sweep and matvec are the same physics evaluated in opposite
directions — and it is exactly the kind of paired closure that drifts when a new
axis lands on one side and is forgotten on the other.

At a multi-moment closure (LD) the source :math:`Q` and the probe
:math:`\bar\psi` carry the trailing :math:`2^d` spatial-moment axis that
:math:`\sigma_t` of shape :math:`(ng, *\text{spatial})` lacks; each moment
scales by the SAME scalar (:math:`1/\sigma_t` on the solve, :math:`\sigma_t` on
the apply), so :math:`\sigma_t` must gain a length-1 trailing axis to broadcast.
This reshape is single-sourced through
:func:`~orpheus.sn.loss_representation._moment_broadcast_sigma`
(:math:`\sigma \mapsto \sigma[\ldots, \text{None}]` iff the moment-valued
operand out-ranks :math:`\sigma`), called by BOTH the sweep ``pure_z`` arm
(:math:`Q\,/\,\texttt{\_moment\_broadcast\_sigma}(\sigma_t, Q)`) and the matvec
``pure_z`` arm
(:math:`\texttt{\_moment\_broadcast\_sigma}(\sigma, \bar\psi)\cdot\bar\psi`), so
the twin cannot diverge on the moment-axis convention.  DD/Step (no moment axis)
:math:`\to` :math:`\sigma_t` unchanged, byte-identical.

.. admonition:: ERR-062 — the matvec twin forgot the guard the sweep already had
   :class: warning

   Before this fix the sweep arm HAD the moment-broadcast guard but the matvec
   arm wrote the bare ``sigma * probe[oct_idx]``, so a moment-valued probe
   broadcast-FAILED.  The consequence:
   ``solve_sn_fixed_source(scheme=LinearDiscontinuous(), inner_solver="krylov")``
   on ANY 2-D Cartesian LD mesh whose quadrature carries pure-z ordinates raised
   ``ValueError`` at the first Krylov matvec.  The bug hid through the whole
   D5b-S3 development because every committed 2-D LD test used
   ``level_symmetric`` — which has NO pure-z ordinates — while the production
   MMS uses a Lebedev quadrature that does.  This is the canonical L21
   twin-path asymmetry recurring a THIRD time ("the matvec needs a committed
   gate, not a round-trip"): the round-trip and FFW :math:`\equiv` MFW gates
   ran on ``level_symmetric`` and never exercised the pure-z arm at all.  The
   gate is
   ``tests/sn/verification/mms/test_mms_ld_2d.py::test_ld_2d_krylov_equals_si_pure_z_quadrature``
   (``@pytest.mark.foundation @pytest.mark.catches("ERR-062")``), on a Mode-9
   degeneracy-break config: a pure-z-bearing Lebedev order-5 quadrature
   (:math:`N = 14`, genuine :math:`\mu_y` + the 2 :math:`\pm z` poles),
   heterogeneous 2-material map, 2-group asymmetric XS with non-zero self-scatter,
   NON-SQUARE :math:`5\times4`, vacuum edges.  Mutation-verified: re-introducing
   the bare ``sigma * probe[oct_idx]`` makes the gate FAIL with the exact
   ``ValueError``; with the fix Krylov :math:`\equiv` SI to :math:`\sim10^{-11}`
   (the same :math:`(L+C-S_{\rm full})` fixed point).  See ``error_catalog.md``
   ERR-062.

The source is :mod:`orpheus.sn.loss_representation` (``_moment_broadcast_sigma``,
the ``loss_action`` matvec ``pure_z`` arm and its source-iteration sweep twin);
the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_purez_gate_closeout.md``.

.. _ld-ubld-moment-scan:

The :math:`d{=}1` moment SCAN (the production sweep) — D5b-S3 OWED-2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The unified moment *matvec* above is the APPLY direction; the production
:math:`d{=}1` LD SWEEP (source iteration) rides the fast Blelloch parallel-prefix
scan (:class:`~orpheus.sn.loss_representation.CumprodScan`), NOT the dense
per-cell solve (L16).  Sub-step **D5b-S3 OWED-2** threads the spatial-moment
iterate :math:`\hat\phi` through that scan so the SI path recovers the SAME
diffusion-limit-consistent operator the matvec does.

.. math::
   :label: ld-ubld-moment-scan-source

   b \;=\; \underbrace{\bar S\,\frac{\mathrm{inv}}{w}}_{\text{flat (average) emission}}
       \;+\; \underbrace{\frac{\theta\,\hat S}{D_2'}
              \;-\; \frac{\theta\,|\mu| A_{\rm down}\,\hat S}{D_2'}\,
                    \frac{\mathrm{inv}}{w}}_{\text{slope source}\ \Sigma_s\hat\phi}

The scan propagates the scalar downstream FACE
:math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` along the cell chain with the
**slope-augmented** affine source :eq:`ld-ubld-moment-scan-source`: the flat
(cell-average) emission :math:`\bar S\,\mathrm{inv}/w` plus the slope-source
contribution :math:`\theta\hat S/D_2' - (\theta|\mu|A_{\rm down}\hat S/D_2')\,\mathrm{inv}/w`
that carries :math:`\Sigma_s\hat\phi` into the recurrence.  Then it reconstructs
the per-cell :math:`(\bar\psi, \hat\psi)` moments from the chained upstream face.
The slope-row :math:`\hat S` algebra is single-sourced through
:meth:`~orpheus.sn.spatial._ubld.D1ClosedForm._slope_fold`, shared by the
per-cell matvec Schur (:meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.schur_xV`)
AND the scan
(:meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.scan_slope_face_source` for the
face-chain term, :meth:`~orpheus.sn.spatial._ubld.D1ClosedForm.scan_reconstruct`
for the per-cell moments).

Why the face/cell split is necessary (the load-bearing math)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A moment-carrying parallel-prefix scan is **NOT** a drop-in widening of the
scalar scan.  For *flat-source* LD the cell average is the convex blend of the
two faces, :math:`\bar\psi = (1-w)\psi_{\rm in} + w\,\psi_{\rm out}`, so the
scalar scan can reconstruct :math:`\bar\psi` directly from the chained faces.
With a *slope* source, that closure **decouples**: :math:`\bar\psi` and
:math:`\psi_{\rm out}` no longer satisfy the convex blend, because the slope
source enters the cell balance without entering the face propagation the same
way.  The scan therefore splits the work in two:

#. the FACE chain :math:`\psi_{\rm out} = a\,\psi_{\rm in} + b` propagates with
   the slope-augmented :math:`b` (:eq:`ld-ubld-moment-scan-source`), so the next
   cell's :math:`\psi_{\rm in}` is the correct dense :math:`\bar\psi + \hat\psi`;

#. the CELL moments :math:`(\bar\psi, \hat\psi)` are reconstructed per cell from
   the chained :math:`\psi_{\rm in}` via the per-cell Schur
   (``scan_reconstruct``), **not** via ``cell_average``.

Conflating the two (using ``cell_average`` for :math:`\bar\psi` on the moment
scan) gives the WRONG cell average while the face chain still looks right — the
silent trap this split avoids.  The reconstruction was verified against a
from-scratch dense :math:`d=1` chain (face / :math:`\bar\psi` / :math:`\hat\psi`
all match to :math:`10^{-12}`) and against the live DAG (scan :math:`\equiv` DAG
to :math:`10^{-16}` on a 2G-heterogeneous non-flat config).

The same global-frame involution as the matvec
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Like the matvec, the scan applies the SAME
:func:`~orpheus.sn.spatial._ubld.octant_moment_frame_signs` involution
(:eq:`ld-ubld-octant-moment-frame-signs`) through the shared ``_reframe``
helper: the source moments are mapped global :math:`\to` sweep on INPUT and the
reconstructed :math:`(\bar\psi, \hat\psi)` sweep :math:`\to` global on OUTPUT, so
the angular reduction :math:`\hat\phi = \sum_n w_n \hat\psi_n`
(:eq:`ld-ubld-slope-angular-reduction`) is frame-consistent and the diffusion
limit is recovered on the source-iteration path too (the sweep-side analog of
the ERR-061 matvec fix).  The backward sweep flips the slope so forward and
backward ordinates reinforce rather than cancel.  The scalar OUTGOING FACE stays
sweep-frame — it propagates along the chain and never crosses into the global
iterate (the :math:`d=1` face cochain is :math:`2^{d-1} = 1`, scalar, so it is
not reframed).  The scan is a **consumer** of the matvec's machinery, never a
twin: the same ``_slope_fold`` powers the per-cell Schur and the scan, and the
same involution powers the DAG and the scan.  DD/Step (``per_axis == 1``) get no
moment axis, so the scan runs the existing flat slab body verbatim and stays
byte-identical (the negative control).

The SymPy module / live scheme is :mod:`orpheus.sn.spatial._ubld`
(``D1ClosedForm``),
:meth:`orpheus.sn.spatial.linear_discontinuous.LinearDiscontinuous.moment_scan_closure`,
and :meth:`orpheus.sn.loss_representation._OneDimScanWalk._run` (the slab
joint-batch moment branch); the gates are
``tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_two_paths_scan_equals_dag_oracle``
(scan :math:`\equiv` DAG) and ``::test_ld_thick_diffusive_limit`` (the diffusion
limit on the SI path, the same ERR-061 catcher the matvec uses); the closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s3_owed2_scan_closeout.md``.

.. _ld-cartesian-2d:

The 2-D Cartesian LD stress MMS (D5b-S4)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sub-step **D5b-S4** is the L1 flux-shape verification of the multi-dimensional
bilinear (UBLD) Linear-Discontinuous closure: a Method-of-Manufactured-Solutions
reference whose trial flux is :math:`\mu`-bilinear (so the per-axis SPATIAL
slope rows are genuinely activated, the vv Mode-7 override) with a NON-vanishing
boundary trace (so the prescribed-inflow boundary closure is stressed).

.. math::
   :label: ld-cartesian-2d

   \psi_{n,g}(x,y) = \frac{1}{W}\bigl[\,A_g(x,y)
       + \mu_{x,n}\,B_g(x,y) + \mu_{y,n}\,C_g(x,y)\,\bigr],
   \qquad
   \phi_g(x,y) = A_g(x,y),

with the strengthened drivers (the :math:`b_2,\,c_2` cross-harmonics break the
x↔y reflection so a same-sign slope-row sign bug cannot cancel) and
:math:`a_0 > 0` (non-vanishing at all four edges).  The manufactured source is
the continuous-PDE residual, derived symbolically (Branch 1, the
algebra-of-record) and structurally independent of the LD cell-update code
(L11).

Why this ansatz (the Mode-7 stress design)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The angular structure :math:`\psi = (A + \mu_x B + \mu_y C)/W` is chosen so the
**per-ordinate** field carries a genuine spatial slope on each axis: the
LD :math:`x`-slope row discretizes :math:`\partial_x\psi`, which sees
:math:`\mu_x\,\partial_x B/W` — a :math:`\mu_x`-weighted slope the DD
cell-average path *cannot* represent (DD has no slope moment).  The scalar flux,
however, is :math:`\phi = \int\psi\,d\mu = A` (the :math:`\mu_x B + \mu_y C`
terms integrate to zero over a symmetric quadrature), so the manufactured scalar
solution is :math:`A` alone — the slope is a *genuinely angular-resolved*
forcing, not a trivial consequence of the average.  This is the vv Mode-7
override: the simplest trig that satisfies the BCs would leave the slope rows
nulled by construction (the classic isotropic-ansatz bias); this ansatz
*activates* them deliberately.  Two further design choices:

* **The :math:`b_2, c_2` cross-harmonics break the x↔y reflection.**  Were
  :math:`B` and :math:`C` related by an :math:`x\leftrightarrow y` reflection, a
  *same-sign* slope-row sign bug (the most likely transcription error, since
  both slope rows share the cell-update code path) could leave the measured
  symmetric flux unchanged — a false green.  Adding distinct cross-harmonic
  content to :math:`B` and :math:`C` (so no reflection maps one to the other)
  makes the :math:`x`-slope-source and :math:`y`-slope-source genuinely
  independent, so a same-sign slope error breaks the measured flux.

* **:math:`a_0 > 0` (non-vanishing at all four edges)** stresses the
  prescribed-inflow boundary closure.  A solution that vanished at the boundary
  by construction would test nothing about the BC handling.  (The curvilinear
  pole-regularity constraint :math:`B(0) = 0` does NOT apply here — a Cartesian
  cell has no :math:`1/r` redistribution, so the slope drivers are unconstrained
  at the boundary.)

The manufactured source is the continuous-PDE residual,
:math:`\mu_x\partial_x\psi + \mu_y\partial_y\psi + \Sigma_t\psi = (1/W)(\Sigma_s\phi + Q^{\rm ext})`,
derived symbolically (Branch 1, the algebra-of-record) and **structurally
independent** of the LD cell-update code (L11): the SymPy derivation never
touches the discretization.  Branch 1 and Branch 2 share their spatial
amplitudes as single-sourced :math:`(\text{num}, \text{den})` pairs
(``Rational`` for SymPy, exact float for numpy), so the Branch-2
:math:`\equiv` Branch-1 source cross-check pins the two *evaluators* agree to
machine precision (:math:`1.5\times10^{-16}`), not just the symbolic identity.

The Branch-1 SymPy derivation lives in
:mod:`orpheus.derivations.continuous.mms.sn`
(:func:`~orpheus.derivations.continuous.mms.sn.derive_2d_cartesian_ld_stress_mms`
and the symbolic builder ``_2d_cartesian_ld_stress_symbolic``); the Branch-2
numerical factory is
:class:`~orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase`
(built by ``build_2d_cartesian_ld_stress_mms_case``).

What this verifies (and what it cannot)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The MMS closes the **slope-UNKNOWN** half of the LM-1989 slope-row sign trap in
:math:`d \ge 2`: the bilinear closure genuinely solves :math:`\hat\psi_x` and
:math:`\hat\psi_y` from the average plus the scattering source, and the
convergence is :math:`O(h^2)` to the manufactured value.  The slope-UNKNOWN sign
is **mutation-verified load-bearing**: on this tightly-coupled UBLD wavefront
closure (the slope feeds the propagating face cochain), a slope-row sign error
does not merely shift the limit — it *diverges* the iteration:

.. list-table:: Mutation verification of the slope-row sign (S4)
   :header-rows: 1
   :widths: 50 26 24

   * - Mutation of the UBLD slope-row sign
     - Strengthened result
     - Verdict
   * - full per-axis gradient sign flip
     - NaN
     - CAUGHT
   * - finite-trace slope :math:`[1,-1]\to[1,1]`
     - order :math:`-4.62`
     - CAUGHT
   * - surgical slope-row :math:`-2 \to +2` (both axes, the faithful "same-sign"
       transcription error)
     - inf
     - CAUGHT

(Baseline strengthened: order 2.00–2.14, finest residual :math:`3.5\text{–}6.0\times10^{-3}`.)
Because the catch is catastrophic (NaN/inf), there is no false-green regime for
this closure to hide in — a stronger guarantee than the subtle-cancellation
scenario the strengthening was originally designed against.

The foundation gate is :mod:`tests.derivations.test_sn_mms_ld_2d_stress_symbolic`
(the SymPy substitution identity, an INDEPENDENT finite-difference residual check
that does not reuse SymPy's own ``diff`` — L11, the Branch-2 :math:`\equiv`
Branch-1 source cross-check, and the Mode-7 activation / x↔y-asymmetry checks);
the end-to-end L1 gates are :mod:`tests.sn.verification.mms.test_mms_ld_2d`
(``test_ld_2d_stress_converges_second_order`` — the headline :math:`O(h^2)` +
value band, ``@l1`` ``@verifies("ld-cartesian-2d", "transport-cartesian-2d")``;
``test_ld_2d_stress_krylov_equals_si`` — the L14 matvec twin on the stress
habitat; ``test_ld_2d_stress_two_paths_ffw_equals_mfw`` — the two-DAG-schedule
invariant).  The closeout is
``.claude/agent-memory/method-implementer/issue_240_d5b_s4_ld_2d_stress_mms_closeout.md``.

The slope-SOURCE sign convention has TWO halves of its own (external
:math:`\hat Q` vs the boundary transverse-face-slope).  Both are now closed —
the EXTERNAL half (**Leg A, #247**) and the BOUNDARY half (**Leg B, #251**) —
see the honest-scope note below.

.. _ld-cartesian-2d-slope-source:

.. note:: Honest scope — the slope-SOURCE half of the LM-1989 trap (Leg A
   VERIFIED #247; Leg B VERIFIED #251).

   The LM-1989 slope-row sign trap has two halves: the slope-UNKNOWN sign
   (always exercised when the slope is non-trivially solved — VERIFIED by this
   MMS; mutation-verified — a slope-UNKNOWN sign flip diverges / leaves the
   value band) and the slope-SOURCE sign :math:`\hat Q`.  The slope-SOURCE
   half splits further:

   - **Leg A — the EXTERNAL slope-moment source** :math:`\hat Q` — is now
     VERIFIED (#247).  :func:`orpheus.sn.solver.solve_sn_fixed_source` accepts
     a typed union of TWO bulk ranks — flat ``(N, ng, *spatial)`` (slope rows
     zeroed, the honest default) OR moment-resolved
     ``(N, ng, *spatial, per_axis**ndim)`` (the projected slope rows threaded
     through) — and ``_lift_external_source_to_moments`` threads the
     moment-resolved slope rows into the SI rhs alongside the scattering
     source.  The slope-source sign is pinned STRUCTURALLY (the converged flux
     is only sub-floor sensitive — the vv Mode 10 trap): the lift threads the
     projected moment vector through at machine precision, and a CONSUMED
     slope-row sign flip moves the converged flux :math:`O(1)` above the inner
     tolerance, while the FLAT scalar gate stays GREEN (the Mode-10 asymmetry
     that closed the gap).

   - **The SCATTERING slope source** :math:`\Sigma_s \hat\phi` (the
     Increment-C iterate feedback) IS consumed and is now mutation-verified
     NOT sign-blind (#247, mutation control M4): flipping the iso slope rows of
     the per-ordinate scattering source moves the converged flux above the
     inner tolerance.  (The old value-band MMS WAS empirically blind to this —
     a slope-source-row sign flip left both the :math:`O(h^2)` order and the
     scalar-flux value band unchanged, because :math:`\Sigma_s\hat\phi` is an
     :math:`O(h)`-small DG-internal forcing whose error enters above
     :math:`O(h^2)` — the canonical vv Mode 10 instance.  The #247 mutation
     control replaces the value-band with the consumption proof.)

   - **Leg B — the BOUNDARY transverse-face-slope** — is now VERIFIED
     (**#251**).  The boundary trace ``mesh.trace`` carries the
     :math:`2^{d-1}` transverse face-moments per face per ordinate per group
     (a moment-resolved slot ``(N, ng, *face_shape, 2^{d-1})`` minted by
     :attr:`orpheus.sn.geometry.SNMesh.boundary_face_layout`, appending the
     single-source :func:`orpheus.numerics.moment_layout.face_moment_tail`),
     so a moment-resolved prescribed inflow can carry the along-face
     (transverse) Legendre slope, the sweep outflow STORES the
     :math:`2^{d-1}` moments (the capture is no longer collapsed to the
     average), and ``_inflow_to_moments`` rank-discriminates a scalar inflow
     (seed slot 0, slopes zero — the scalar default) from a moment-resolved
     inflow (thread the projected transverse slope through).  Like Leg A the
     converged near-boundary flux is only sub-floor sensitive to the boundary
     slope (a SHARPER vv Mode 10 — "improves-on-flat" is NOT achievable,
     because the localized :math:`O(h)`-small boundary-trace slope sits below
     the bulk :math:`O(h^2)` discretization floor), so the slope is pinned
     STRUCTURALLY: the producer threads the projected transverse face-slope
     into the cochain at machine precision, and a CONSUMED transverse-slope
     sign flip moves the converged near-boundary flux :math:`O(1)` above the
     inner tolerance (:math:`|\Delta\phi|/|\phi| \approx 4.1\times10^{-3}`
     near-boundary at ``nc=16``, ~5.6 orders above the consumption tolerance;
     linear in the slope magnitude — genuine consumption), while the SCALAR
     inflow gate stays byte-identical (the Mode-10 asymmetry).  DD/Step
     (``per_axis == 1`` → ``face_moment_tail(1) == ()``) leaves the trace
     byte-identical (the negative control); a 1-D slab face is a point
     (``face_shape == ()``), so the transverse face-moment is a 2-D-and-higher
     concern by construction.

     The transverse-slope SIGN under REFLECTION across a face is a separate
     follow-up: the Leg-B MMS is vacuum-BC (which nulls the reflective
     coupling), so the reflective ``B`` operator's moment-axis passthrough
     (verified storage-correct — the ``PermutationOperator(axis=0)``
     broadcasts over the new trailing moment axis without a hard-coded
     trailing-axis assumption) is NOT exercised for its sign.  Physics: a
     normal-flip reflection preserves the tangent-plane (transverse)
     coordinate, so the transverse slope should reflect WITHOUT a sign flip —
     but this is UNVERIFIED (a reflective-LD MMS + an ``op.H`` adjoint check
     on the transverse-slope reflection is the follow-up).

   The full Leg A narrative — the tensor-Legendre projection convention, the
   typed-union bulk widening, the Mode-10 structural-teeth design, and the
   M1–M4 mutation table — is the subsection :ref:`ld-cartesian-2d-legA`
   immediately below.  The full Leg B narrative — the
   :attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout` moment-tail
   storage lever, the ``_inflow_to_moments`` rank-discriminated pass-through,
   the four outflow capture-collapse DROP sites, the
   ``prescribed_inflow`` scalar-or-moment producer, the transverse
   face-moment normalization, the SHARPER Mode-10 (no improves-on-flat leg),
   and the reflective-BC sign follow-up (#252) — is the subsection
   :ref:`ld-cartesian-2d-legB`.

   **S9 (#257)** completes the boundary half: the MMS case
   ``prescribed_inflow`` now itself EMITS the moment-resolved slot (it no longer
   drops the slope at the producer), and the **coherent promise** — *LD is
   second-order at the boundary with no asterisk* — is LOCKED by a dedicated
   first-cell-row convergence gate.  The promise is delivered by the AVERAGE
   moment alone; the transverse slope is a sub-floor inflow-representation
   refinement (the fourth vv Mode-10 companion-unavailable instance).  S9 also
   establishes that the transverse boundary moment is a PROPERTY, not a new
   field type (#263).  The full S9 narrative — the coherent-promise evidence,
   the producer-honesty carve, and the property-vs-type seam — is the subsection
   :ref:`ld-cartesian-2d-coherent-promise`.

.. _ld-cartesian-2d-legA:

Leg A — the external slope-moment source :math:`\hat Q` (#247)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This subsection is the rich record of the EXTERNAL half of the slope-SOURCE
trap, closed under Issue #247.  The change is small in code (two named sites in
:func:`orpheus.sn.solver._build_fixed_source_rhs` and
``_lift_external_source_to_moments``) but it is the first time an LD external
source can carry sub-cell slope information through the public solver, so the
*verification* design — not the code — is the load-bearing content here.  Why
the slope-SOURCE sign is genuinely hard to pin (the vv **Mode 10** trap) and how
the gate gets teeth anyway is the canonical resolution of an
activated-but-unconstrained term; it is recorded in full because the lesson
recurs whenever a term is consumed yet enters the measured quantity below the
discretization floor.

The tensor-Legendre projection convention (the CRUX)
''''''''''''''''''''''''''''''''''''''''''''''''''''

To feed a slope-resolved external source into the UBLD closure, a continuous
:math:`Q^{\rm ext}(x,y)` must be projected onto the per-cell tensor-Legendre
moment vector that the cell update consumes.  The single load-bearing decision
is *which normalization* the projected moments carry — and the answer is locked
by the UBLD mass matrix, not chosen for convenience.

The UBLD cell mass on one axis is :math:`M_{\rm 1d} = \mathrm{diag}(h,\,\theta
h)`, :math:`\theta = 1/3` (Eq. :eq:`ld-ubld-mass-weights`), the L2-Gram of the
Legendre moment basis :math:`\{P_0 = 1,\ P_1 = \xi\}` on :math:`\xi \in [-1,1]`:
:math:`\langle P_0, P_0\rangle = h`, :math:`\langle P_1, P_1\rangle = \theta h =
h/3`.  The cell kernel forms its right-hand side as :math:`R_{\rm source} = M\,
S_{\rm moments}` (the d=1 reduction Eq. :eq:`ld-ubld-d1-reduction` confirms
:math:`R_{\rm prod} = [\,\bar Q\,h,\ \theta\,\hat Q\,h\,]` symbolically).  The
mass matrix therefore *already* supplies the per-volume and the
:math:`\theta` weighting.  The projection must NOT duplicate it.  The projected
moment is the **bare per-volume Legendre coefficient**:

.. math::
   :label: ld-cartesian-2d-projection-coeff

   \bar q \;=\; \frac{1}{V}\!\int_{\rm cell} q \;=\; \text{(cell average)}
   \quad\text{(slot 0)},
   \qquad
   \hat q_a \;=\; \frac{\langle q,\,P_1(\xi_a)\rangle}
                       {\langle P_1,\,P_1\rangle}
              \quad\text{(the }P_1\text{ coefficient on axis }a).

For a cell-linear source :math:`q = a + b\,x`, the axis coefficient is
:math:`\hat q = b\,h/2` — **no** :math:`\theta`, **no** :math:`h`, **no**
:math:`V` in the projected number; the kernel's :math:`M` adds them downstream.
Sharing the :math:`\theta`/:math:`h` weighting between the projection and the
mass would double-count it; this is the apples-to-apples constraint that makes
the projected slope rows match what :math:`M^{-1}R` expects.

For a general bilinear :math:`q = a_{00} + a_{10}x + a_{01}y + a_{11}xy` on a
cell :math:`[x_L,x_R]\times[y_L,y_R]` (:math:`h_x = x_R-x_L`, centre
:math:`x_c`; similarly :math:`y`), the four tensor-Legendre coefficients are
hand-derivable in closed form:

.. math::
   :label: ld-cartesian-2d-bilinear-coeffs

   \bar q   &= a_{00} + a_{10}x_c + a_{01}y_c + a_{11}x_c y_c, \\
   \hat q_y &= \tfrac{h_y}{2}\,(a_{01} + a_{11}x_c), \\
   \hat q_x &= \tfrac{h_x}{2}\,(a_{10} + a_{11}y_c), \\
   \hat q_{xy} &= \tfrac{h_x}{2}\,\tfrac{h_y}{2}\,a_{11}.

These four numbers are the structurally-independent reference for the projector
(see the teeth subsection): a bilinear integrand is integrated exactly by a
2-point Gauss rule, so the quadrature projector reproduces them to machine
precision.

**The d=2 Kronecker moment order is** :math:`[\bar\psi,\ \hat\psi_y,\
\hat\psi_x,\ \hat\psi_{xy}]` — axis 0 (:math:`x`) is the OUTER Kronecker
factor, axis 1 (:math:`y`) the INNER, so the slot order
:math:`[(0,0),(0,1),(1,0),(1,1)]` places the **x-slope at slot 2**, the y-slope
at slot 1, the cross-moment at slot 3 (consistent with
Eq. :eq:`spatial-moment-kronecker-order`).  The cell mass diagonal satisfies
:math:`\mathrm{diag}(M) = \mathrm{diag}(h_x,\theta h_x) \otimes
\mathrm{diag}(h_y,\theta h_y)`, the Kronecker product that fixes which slot is
which.

**The projection supplies GLOBAL-frame coefficients.**  The natural Legendre
coefficients of :math:`q(x,y)` live in the global :math:`x`/:math:`y` frame; the
per-octant sweep frame (where a downwind axis runs the other way) is *not* the
projection's concern.  Production reframes the source global→sweep per octant in
:meth:`~orpheus.sn.sweep_graph._CellSolve.cell` via the slope-sign involution
:math:`\mathrm{octant\_moment\_frame\_signs}`
(Eq. :eq:`ld-ubld-octant-moment-frame-signs`,
:func:`orpheus.sn.spatial._ubld.octant_moment_frame_signs`), exactly as it
reframes the scattering slope source.  So the external :math:`\hat Q` rides the
SAME global→sweep machinery the scattering moments already use — the
:ref:`ld-ubld-sweep-global-frame` involution that the S3 unified matvec had to
get right (ERR-061) is reused unchanged, with no new cell branch.

The projector is structurally independent of production (L11): it evaluates
:math:`\int q\,P_k` with :func:`numpy.polynomial.legendre.leggauss` directly and
NEVER calls ``_lift_external_source_to_moments`` nor any
:class:`~orpheus.sn.spatial.linear_discontinuous.LinearDiscontinuous` cell op.
The reference (Eq. :eq:`ld-cartesian-2d-bilinear-coeffs`) is hand-laid
polynomial algebra, so "the projector matches the reference" is not a production
echo.  One subtlety, pinned by its own sub-gate: the projector returns the cell
**average** in slot 0, whereas the legacy flat producer ``case.external_source``
evaluates :math:`Q` at the cell **centre**; the two differ by :math:`O(h^2)`
(slot-0 ratio :math:`\sim 0.93` at ``nc=8``).  The projector slot-0 is therefore
cross-checked against an *independent* fine-quadrature cell average, NOT against
the cell-centre producer (which would falsely fail by :math:`O(h^2)`).

The typed-union bulk widening
'''''''''''''''''''''''''''''

Before #247, :func:`~orpheus.sn.solver.solve_sn_fixed_source` accepted a single
flat bulk shape :math:`(N, n_g, *\text{spatial})` and rejected everything else.
The widening makes the bulk a **typed union of two ndarray ranks**, discriminated
by RANK, not trailing size:

.. list-table:: The widened bulk-source contract (``_build_fixed_source_rhs``)
   :header-rows: 1
   :widths: 26 24 50

   * - Bulk shape
     - Closure
     - Meaning
   * - :math:`(N, n_g, *\text{spatial})` (flat)
     - any
     - the original path — slope moments :math:`\hat Q` zeroed by the lift (the
       honest default, exact for a region-uniform source).  Byte-identical to
       pre-#247.
   * - :math:`(N, n_g, *\text{spatial},\, \text{per\_axis}^{\,d})`
       (moment-resolved)
     - LD only (``per_axis > 1``)
     - the caller projected :math:`Q^{\rm ext}` onto the tensor-Legendre moment
       vector; the lift threads the slope rows through to join the
       moment-carrying scattering source :math:`\Sigma_s\hat\phi` in the SI rhs.
   * - anything else
     - any
     - ``ValueError`` (see the negative pin)

**Why discriminate by RANK, not trailing size.**  A moment-resolved bulk has
exactly one more axis than a flat bulk; a coincidental spatial dimension could
happen to equal :math:`2^d`, so testing the trailing-axis *length* would
misclassify a flat bulk whose last spatial dim is 4.  The rank is unambiguous: a
flat bulk has :math:`2 + |{\rm spatial}|` axes, a moment-resolved bulk has one
more.  ``_lift_external_source_to_moments`` makes the same rank decision
(``bulk_values.ndim == 2 + len(spatial_shape)`` is the flat-rank test) — the
moment-layout primitive :func:`orpheus.numerics.moment_layout.is_moment_valued_by_rank`
is the canonical discriminator for the cell-internal path.

**Why DD/Step rejects a moment bulk.**  At a flat closure
(:math:`\text{per\_axis} = 1`, hence :math:`n_{\rm cell\ moments} = 1`) there is
NO moment axis — the cell carries a single average, with no slope to fill.  A
moment-resolved input there is a category error, so the validation rejects it
outright (only flat is valid).  This is Pattern 4 (illegal states made
unrepresentable): the relaxation admits exactly the two principled ranks and
nothing in between.

**The negative pin.**  The relaxation must not swallow a real shape bug.  A
moment-resolved bulk whose trailing axis :math:`\neq \text{per\_axis}^{\,d}`
(e.g. a 5-wide axis on a 2-D LD mesh where :math:`2^d = 4`) raises a
``ValueError`` that names the expected :math:`2^d` and the full moment-vector
shape.  The gate
``test_moment_resolved_bulk_still_rejects_wrong_trailing_axis`` pins both arms:
the LD 5-wide reject AND the DD 4-wide reject.

The lift then has three arms, single-sourced for the fixed-source and eigenvalue
paths:

.. list-table:: ``_lift_external_source_to_moments`` (the three arms)
   :header-rows: 1
   :widths: 30 70

   * - Input
     - Action
   * - DD/Step (``tail == ()``)
     - input returned UNCHANGED — byte-identical, the backward-compat negative
       control.
   * - flat :math:`(N, n_g, *\text{spatial})`
     - zero the :math:`2^d` buffer, copy the flat values onto slot 0 (average);
       slope rows stay ZERO (:math:`\hat Q = 0`, the honest default).
   * - moment-resolved :math:`(N, n_g, *\text{spatial},\, 2^d)`
     - thread the moment vector through UNCHANGED (validate the trailing axis);
       the slope rows the caller projected reach the SI rhs.

No callable-projection entry is exposed (Pattern 6, defer abstraction): there is
no production consumer that needs the solver to project a continuous source, and
adding one would make the verification a tautology (the gate would compare the
production projector to itself).  The MMS test does its OWN projection and passes
the array — structurally independent of production by construction (L11).

The Mode-10 structural-teeth design
'''''''''''''''''''''''''''''''''''

The slope-SOURCE sign is the textbook vv **Mode 10** trap (an
*activated-but-unconstrained* term): the slope-source code path is genuinely
exercised — the slope rows are populated, threaded, reframed per octant, and
consumed by the cell update — yet a sign flip on those rows does **not** move the
converged scalar flux above the discretization floor.  The reason is that the
slope-source contribution enters the converged flux as an :math:`O(h^2)`-small
forcing that rides *on top of* the :math:`O(h^2)` discretization error.  Probed
live, the average-moment L2 error vs :math:`\phi_{\rm exact}` under an
x-slope-source sign flip:

.. list-table:: Why the converged flux is sub-floor sensitive to the slope-source sign
   :header-rows: 1
   :widths: 30 40 30

   * - ``nc``
     - correct slope (L2 err / order)
     - x-slope FLIPPED (L2 err / order)
   * - 16
     - :math:`8.18\times10^{-3}`
     - :math:`1.17\times10^{-2}`
   * - 32
     - :math:`1.99\times10^{-3}` (2.04)
     - :math:`2.96\times10^{-3}` (1.98)
   * - 64
     - :math:`4.86\times10^{-4}` (2.03)
     - :math:`7.38\times10^{-4}` (2.01)
   * - 128
     - :math:`1.18\times10^{-4}` (2.05)
     - :math:`1.81\times10^{-4}` (2.03)

Both converge at clean :math:`O(h^2)`; the flipped error is only
:math:`\sim 1.4\text{–}1.5\times` larger at every mesh, and the ratio is roughly
CONSTANT under refinement.  Two consequences for the gate:

* **A convergence-ORDER leg is blind to the sign.**  The order stays 2 both ways
  — :math:`O(h^2)` to the wrong limit is still :math:`O(h^2)` (the vv §5
  warning).
* **A fixed-mesh value-band is too fragile.**  A band that separates the correct
  from the flipped converged flux would need a tolerance tighter than the
  :math:`\sim 1.5\times` gap, and the :math:`O(h^2)` discretization error itself
  eats that margin.  The smallest signal — the :math:`xy` cross-slope — moves
  the converged flux only :math:`\sim 6\times10^{-5}` relative; no value-band
  survives that.

So the teeth do **not** come from the converged flux.  They come from two places
where the sign flip is :math:`O(1)`:

1. **The lift threads the projection through at machine precision** (the
   production-change proof).  ``test_ld_2d_external_slope_source_threaded_through_lift``
   feeds the projected moment vector to the production lift and asserts the
   returned moment source equals the projection EXACTLY — every slope slot, via
   ``np.testing.assert_array_equal``.  A regression that re-zeroes the slope rows
   (the EXACT bug #247 closes) breaks this at :math:`O(1)`, where the converged
   flux would never catch it.  A NEGATIVE-CONTROL leg in the same test pins that
   a FLAT bulk still lifts onto slot 0 with the slope rows EXACTLY ZERO (the
   honest default is preserved).
2. **A consumed slope-row sign flip moves the converged flux :math:`\gg` solver
   tolerance** (the consumption proof).  The inner solve converges to
   :math:`10^{-12}`; a flip of a CONSUMED slope row moves the flux by
   :math:`\sim 3\times10^{-3}` (x), :math:`\sim 10^{-2}` (y), or
   :math:`\sim 6\times10^{-5}` (xy) relative.  The acceptance band
   ``_CONSUMPTION_TOL = 1e-8`` sits :math:`\sim 5\times10^7\times` above the
   fixed point yet far below the §0 trap — the smallest probed flip (xy) clears
   it by :math:`\sim 6000\times`.  This is sharp BECAUSE the test contrasts two
   solves of the *same* problem that differ only in a sign, so the
   discretization floor cancels and the slope-source contribution is the signal.

The convergence-ORDER leg
(``test_ld_2d_external_slope_source_converges_second_order``,
``@verifies("ld-cartesian-2d")``) is kept as a NECESSARY check — it proves the
threaded slope rows are CONSISTENT (the slope-unknown plus the source together
produce a 2nd-order moment, probed :math:`8.18\times10^{-3} \to
1.99\times10^{-3}` at ``nc 16 → 32``) — but it is explicitly **not** the sign
teeth.  A POSITIVE leg
(``test_ld_2d_external_slope_source_improves_on_flat``) closes the loop: the
moment-resolved solve lands strictly closer to :math:`\phi = A` than the
flat-in-moment solve (:math:`3.4\times10^{-3} < 5.9\times10^{-3}` at ``nc=24``),
so the threaded slopes carry real sub-cell information, not noise.

.. note:: **A current-invariant lesson worth recording (vv Mode 10).**

   A Mode-10 gap is closed NOT by tightening the converged-flux value band (the
   :math:`O(h^2)`-small forcing is sub-floor) but by two :math:`O(1)` structural
   teeth: (1) assert the production *producer* threads the projection through at
   machine precision (catches a regression to zeroing), and (2) assert a
   *consumed* source-row sign flip moves the converged answer :math:`\gg` solver
   tol (catches sign-blindness), paired with the FLAT no-op leg that pins the
   asymmetry (the old scalar gate is correctly blind, by construction).  The
   convergence-order leg is necessary for slope consistency but is not the sign
   teeth.  This is the canonical resolution whenever a term is genuinely consumed
   yet its error enters the measured quantity below the convergence floor.

The mutation-control table (M1–M4)
''''''''''''''''''''''''''''''''''

The primary sign-catchers are mutation controls (vv anti-pattern #11 — a
``catches`` claim is verified by re-introducing the exact bug and confirming the
gate reddens).  Two distinct mutations stress two distinct source paths: the
EXTERNAL :math:`\hat Q` (the NEW #247 consumption, M1–M3) and the SCATTERING
:math:`\Sigma_s\hat\phi` (the EXISTING S3 consumption, M4).  Each flips a slope
row and asserts the converged flux changes :math:`\gg` ``_CONSUMPTION_TOL`` while
the FLAT scalar gate stays GREEN — that asymmetry IS the Mode-10 gap being
closed.

.. list-table:: The slope-SOURCE sign mutation controls
   :header-rows: 1
   :widths: 6 30 34 30

   * - \#
     - Source row flipped
     - The NEW moment gate must
     - The FLAT scalar gate
       (``..._stress_converges_second_order``)
   * - M1
     - EXTERNAL :math:`\hat Q` x-slope (slot 2)
     - go RED — converged flux moves :math:`\sim 3\times10^{-3}` (the
       consumption proof)
     - stays GREEN — it feeds a flat source, slope row already 0, flipping zero
       is a no-op
   * - M2
     - EXTERNAL :math:`\hat Q` y-slope (slot 1)
     - go RED — flux moves :math:`\sim 10^{-2}`
     - stays GREEN (flat → no-op)
   * - M3
     - EXTERNAL :math:`\hat Q` cross-moment (slot 3)
     - go RED — flux moves :math:`\sim 6\times10^{-5}` (the weakest signal,
       still :math:`\sim 6000\times` over tol)
     - stays GREEN (flat → no-op)
   * - M4
     - SCATTERING :math:`\Sigma_s\hat\phi` iso slope rows (slots 1:)
     - go RED — flux moves :math:`\sim 2.6\times10^{-3}`
     - stays GREEN — the scalar gate's converged flux is only :math:`\sim
       1.4\times` sensitive (sub-floor)

**M1–M3 verify the NEW external consumption.**  Before #247 the lift zeroed the
external slope rows, so a flip of an already-zero row was a no-op and the
"flipped reddens" assertion could not hold — these mutations only become
catchers once the lift threads the slope rows.  The test flips slot
:math:`\{2,1,3\}` of the projected :math:`\hat Q` and re-solves the full public
solve; the same flip on the FLAT source is then asserted to be a no-op directly
(``flat_lift[..., slot]`` is exactly zero), pinning the asymmetry that closed
the gap.

**M4 verifies the EXISTING scattering consumption was never sign-blind.**  The
scattering slope source :math:`\Sigma_s\hat\phi` (the Increment-C iterate
feedback, Eq. :eq:`ld-ubld-scattering-moment-lift`) has been consumed since S3,
but the OLD value-band MMS was empirically blind to its sign — a slope-source-row
sign flip left both the :math:`O(h^2)` order and the scalar-flux value band
unchanged, because :math:`\Sigma_s\hat\phi` is an :math:`O(h)`-small DG-internal
forcing whose error enters above :math:`O(h^2)`.  M4 monkeypatches
``ScatteringOperator._assemble_per_ordinate_source`` (the :math:`\Sigma_s
\otimes I` over every spatial moment) to negate the iso slope rows and confirms
the converged flux moves :math:`\sim 2.6\times10^{-3}` — the consumption proof
replaces the value-band the old gate relied on.

Each mutation is reverted in a ``finally`` block, and all #247 gates are
``-O``-safe (``np.testing.*`` / ``pytest.fail`` / ``pytest.raises`` only, no
bare ``assert`` that ``python -O`` would strip — vv Mode 8).

No ERR entry was minted: Mode 10 here is a proactive-gap close, not a caught
production bug.  The lift correctly zeroed an unverified-but-honest default
(:math:`\hat Q = 0`); the slope-source sign was UNVERIFIED, not WRONG.  Per the
"log every caught bug" directive, an ``@catches`` marker is added only when a
real production bug surfaces; none did.

Sources and gates
'''''''''''''''''

The production change is in :func:`orpheus.sn.solver._build_fixed_source_rhs`
(the typed-union validation) and ``_lift_external_source_to_moments`` (the
slope-thread arm), both confined to ``solver.py``.  The end-to-end gates live in
:mod:`tests.sn.verification.mms.test_mms_ld_2d` (the #247 block):
``test_ld_2d_external_slope_source_threaded_through_lift`` (the foundation
structural teeth), ``..._converges_second_order`` and ``..._improves_on_flat``
(the L1 necessary + positive legs), ``..._sign_mutation_reddens`` (M1–M3),
``test_ld_2d_scattering_slope_source_sign_mutation_reddens`` (M4),
``test_moment_resolved_bulk_still_rejects_wrong_trailing_axis`` (the negative
pin), with the projection-correctness foundation sub-gates
``test_tensor_legendre_projection_matches_hand_polynomial`` and
``test_projection_slot0_is_cell_average_not_centre``.  The bit-identity of the
flat/DD path is guarded by the strict ``DriftWarning`` regression gate (no golden
moved — the typed-union widening leaves the flat path byte-identical).

.. _ld-cartesian-2d-legB:

Leg B — the boundary transverse-face-slope (#251)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This subsection is the rich record of the BOUNDARY half of the slope-SOURCE
trap, closed under Issue #251 (split from #247).  Where Leg A widened the BULK
external source carried through the *interior*, Leg B widens the BOUNDARY TRACE
``mesh.trace`` so a moment-resolved prescribed inflow can carry the along-face
(transverse) Legendre slope and the sweep outflow can STORE the
:math:`2^{d-1}` transverse face-moments instead of collapsing them to the
average.  It is the boundary twin of Leg A, and the structural-teeth template
(:ref:`ld-cartesian-2d-legA`) carries over almost verbatim — but the
verification design is *sharper*: the boundary slope is sub-floor for ANY value
claim, not just its sign (the "improves-on-flat" leg that closed Leg A is
**unachievable** here).  That sharpening — the first vv **Mode 10** instance
where the O(1)-isolating-companion half of the recipe is genuinely unavailable —
is the load-bearing lesson, and the reason this is recorded in full.

The ``boundary_face_layout`` moment-tail storage lever (the CRUX)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Leg A had a free ride: the BULK iterate :math:`\hat\phi` already carried its
per-cell :math:`2^d`-moment axis (the S3 unified moment matvec), so widening the
external source meant relaxing a single lift's input contract.  Leg B has no
such carrier — *the boundary trace is scalar-per-face end-to-end*.  A new place
must be found to STORE the transverse face-moments, and the elegant answer is a
single attribute on the mesh.

The trace's per-face slot shape is owned not by the
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` itself but by the
:class:`~orpheus.numerics.face_layout.FaceLayout` it is built from, and that
layout is minted by :attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout`.
The widening is therefore ONE site: append the scheme's per-face transverse
moment tail to each slot,

.. math::
   :label: ld-cartesian-2d-face-slot-shape

   \text{slot}(\text{face}) \;=\;
   \bigl(N,\ n_g,\ *\,\text{face\_shape}(\text{label}),\
         \underbrace{*\,\text{face\_moment\_tail}
            \bigl(\text{per\_axis}^{\,d-1}\bigr)}
         _{\text{the new }2^{d-1}\text{ tail}}\bigr),

where the per-face moment count is :math:`n_{\rm face} = \text{per\_axis}^{\,d-1}`
(the FACE tail) — note the exponent :math:`d-1`, NOT the cell-tail exponent
:math:`d` of Eq. :eq:`spatial-moment-kronecker-order`.  A face is a codimension-1
object: it carries a Legendre moment per *transverse* axis only (the
:math:`d-1` axes that run along the face), so its moment count is the cell count
divided by the one normal-axis factor.

Three properties make this the clean lever:

* **The** :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` **needs ZERO
  changes — it was "moment-ready by accident".**  The trace's partial-current
  metric (the :math:`|\Omega\!\cdot\!n|\,w_n` weighting) and its
  ``omega_dot_n`` inflow/outflow table are **per-ordinate** (they classify and
  weight by ordinate, axis 0) and broadcast over ALL trailing axes by
  construction.  A moment axis appended to the slot rides the metric and the
  directional selectors for free.  The trace space was already
  moment-polymorphic; only its layout-supplier was scalar.  This is the
  illegal-states-unrepresentable payoff of the Depth-B field-space refactor: the
  boundary FIELDS (``BoundaryFlux`` / ``BoundarySourceSink`` /
  ``BoundaryResidual`` / ``BoundaryDisplacement``) validate ONLY against
  ``layout.total_size``, never a hardcoded :math:`(N, n_g)`, so they accommodate
  any slot shape the layout dictates.

* **DD/Step is byte-identical (the negative control).**  At a flat closure
  :math:`\text{per\_axis} = 1`, so :math:`n_{\rm face} = 1^{\,d-1} = 1` and
  :func:`~orpheus.numerics.moment_layout.face_moment_tail` returns ``()`` (the
  "append iff > 1" policy — NO length-1 axis).  Every DD/Step slot shape is
  untouched, so every buffer, snapshot, and metric is bit-for-bit unchanged.
  This is the SAME single-source tail policy the interior cell cochain keys on
  (:attr:`_LossRepresentation._n_face_moments` and
  :attr:`_LossRepresentation._spatial_moment_tail`), reused so the storage and
  the cochain can never disagree on the shape.

* **Leg B is a 2-D-and-higher concern by construction.**  A 1-D slab face is a
  *point*: ``face_shape == ()`` and there is no transverse axis, so
  :math:`n_{\rm face} = \text{per\_axis}^{\,0} = 1` even for the 2-basis LD
  closure.  A point has no along-face direction, hence no transverse slope; the
  1-D prescribed-inflow MMS is byte-identical not by coincidence but because the
  exponent :math:`d-1` vanishes.

The scheme is reachable: ``SNMesh`` sets ``self.scheme`` before it builds the
trace, so ``boundary_face_layout`` can read ``self.scheme.spatial_basis_per_axis``
to compute the tail.  Verified live: with LD the face slots are
:math:`(24, 2, 6, 2)` / :math:`(24, 2, 8, 2)` (the trailing ``2`` is the
:math:`2^{d-1}` transverse-moment axis); with DD they are :math:`(24, 2, 6)` /
:math:`(24, 2, 8)` (no axis).

The transverse face-moment normalization (apples-to-apples with the cochain)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The same normalization decision that locked Leg A (Eq.
:eq:`ld-cartesian-2d-projection-coeff`) recurs, transposed to the transverse
axis.  The cochain consumes the upstream face's moments through
:func:`orpheus.sn.spatial._ubld.assemble_inflow_axis`, which weights the
:math:`2^{d-1}`-moment face vector by the **transverse mass** — the Kronecker
product of the per-transverse-axis :func:`~orpheus.sn.spatial._ubld.mass_1d`
:math:`= \mathrm{diag}(h_t,\,\theta h_t)`, :math:`\theta = 1/3` — before applying
the active-axis trace :math:`B(-1) = [1, -1]` and :math:`|\mu_{\rm axis}|`.  The
mass therefore ALREADY supplies the transverse :math:`h_t` and :math:`\theta`
weighting; the trace must NOT duplicate it.

So — exactly as Leg A's cell mass :math:`M = \mathrm{diag}(h, \theta h)` forced
the projected source rows to be bare per-volume coefficients — **the boundary
trace must carry the BARE per-transverse Legendre coefficients**:

.. math::
   :label: ld-cartesian-2d-face-projection-coeff

   b_{\rm bar} \;=\; \frac{\langle\psi_{\rm face},\,P_0\rangle}
                          {\langle P_0, P_0\rangle}
            \;=\; \text{(transverse cell average)}
            \quad\text{(slot 0)},
   \qquad
   b_{\rm slope} \;=\; \frac{\langle\psi_{\rm face},\,P_1(\xi)\rangle}
                            {\langle P_1, P_1\rangle}
            \quad\text{(slot 1, the bare transverse }P_1\text{ coeff)},

with :math:`\xi \in [-1,1]` the transverse coordinate.  For a face-linear inflow
:math:`\psi_{\rm face}(t) = c_0 + c_1 t` on a transverse cell
:math:`[t_L, t_R]` (width :math:`h_t = t_R - t_L`, centre :math:`t_c`), the two
coefficients are hand-derivable in closed form,

.. math::
   :label: ld-cartesian-2d-face-bilinear-coeffs

   b_{\rm bar}   &= c_0 + c_1\,t_c \quad\text{(the cell AVERAGE, not the centre
                    eval — see below)}, \\
   b_{\rm slope} &= \tfrac{h_t}{2}\,c_1 \quad\text{(no }\theta\text{, no }h_t
                    \text{ beyond the }\tfrac{1}{2}\text{; the mass adds them)},

the structurally-independent reference the face projector is pinned against (a
linear integrand is integrated exactly by a 2-point Gauss rule, so the
``leggauss`` projector reproduces them to machine precision).  This is the
1-D-transverse *factor* of Leg A's tensor projector
(Eq. :eq:`ld-cartesian-2d-bilinear-coeffs`): a face projection is the per-axis
Legendre coefficient of the tensor projection along the single transverse axis.

.. note:: **The apples-to-apples trap (the same crux that bit Leg A).**

   If the trace stored a :math:`\theta`- or :math:`h_t`-weighted slope instead
   of the bare coefficient, the cochain's transverse mass would double-apply the
   weighting, and the threaded slope would be wrong by a constant
   :math:`\theta h_t` factor.  The structural threading gate (below) compares the
   trace's slot-1 against the bare projected reference precisely so that this
   double-counting is caught at :math:`O(1)`.

**The slot-0 centre-vs-average subtlety (sharper than Leg A).**  Today's scalar
trace carries the cell-**centre** eval of the manufactured inflow
:math:`\psi_{\rm face}(t_c)/W`, whereas the projection's slot 0 is the cell
**average** :math:`b_{\rm bar}`; the two differ by :math:`O(h^2)`.  The
backward-compat decision is to keep slot 0 = the EXISTING scalar trace (centre)
on the scalar-inflow path — so DD/Step and the existing 1-D prescribed-inflow
MMS stay byte-identical — and to carry the (average-bar, slope) pair only when a
moment-resolved inflow is explicitly supplied.  The structural threading gate
therefore compares the SLOPE slot (slot 1) only; slot 0 may legitimately differ
centre-vs-average.  A dedicated foundation sub-gate
(``test_face_projection_slot0_is_transverse_cell_average``) pins that the
projector's slot 0 IS the average — cross-checked against an *independent*
fine-quadrature transverse cell average, NOT against the cell-centre
``case.prescribed_inflow`` (which would falsely fail by :math:`O(h^2)`).

The face projector is structurally independent of production (L11): it evaluates
:math:`\int \psi_{\rm face}\,P_k` with
:func:`numpy.polynomial.legendre.leggauss` directly and NEVER calls
``_inflow_to_moments`` nor ``assemble_inflow_axis`` nor any LD cell op.

The rank-discriminated inflow lift and the four DROP sites
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The boundary twin of Leg A's ``_lift_external_source_to_moments`` is
:meth:`_LossRepresentation._inflow_to_moments`.  It rank-discriminates the
incoming face against the FLAT (moment-free) face rank :math:`d + 1` — a scalar
face :math:`(N_{\rm oct}, n_g, *\text{transverse})` has rank
:math:`2 + (d-1) = d + 1`, since the transverse part carries :math:`d-1` axes —
using the same single-source primitive Leg A introduced,
:func:`orpheus.numerics.moment_layout.is_moment_valued_by_flat_rank` (no third
rank spelling — the layout-shift hazard is an open-coded ``ndim == flat_ndim``
that silently diverges from its sibling consumers).  Three arms:

.. list-table:: ``_inflow_to_moments`` (the three arms)
   :header-rows: 1
   :widths: 30 70

   * - Input
     - Action
   * - DD/Step (:math:`n_{\rm face} = 1`)
     - identity — the trailing axis is absent, every buffer byte-identical (the
       backward-compat negative control).
   * - scalar inflow :math:`(N_{\rm oct}, n_g, *\text{transverse})`
     - widen — zero the :math:`2^{d-1}` buffer, seed the AVERAGE moment (slot 0)
       from the scalar, leave the transverse SLOPES zero (a scalar trace carries
       no along-face variation; **the scalar default — the Leg-B asymmetry**).
   * - moment-resolved inflow :math:`(N_{\rm oct}, n_g, *\text{transverse},\,
       2^{d-1})`
     - PASS THROUGH — the widened trace already carries the projected transverse
       face-slope; thread it unchanged (validate the trailing width
       :math:`= 2^{d-1}`, ``ValueError`` otherwise — Pattern 4).

No callable-projection entry is exposed (Pattern 6, defer abstraction): the MMS
does its OWN projection and passes the moment array; production accepts the
moment-resolved face, it does not compute it — structurally independent by
construction (L11), exactly Leg A.

The full-field oracle's seed, :meth:`FullFieldWavefront._octant_face_cochain`,
mirrors the same rank discriminator at the IN-edge: a scalar inflow seeds slot 0
(transverse slopes zero), a moment-resolved inflow seeds ALL :math:`2^{d-1}`
moments.  So the two storage strategies — the production
:class:`~orpheus.sn.loss_representation.MovingFrontierWindow` and the
:class:`~orpheus.sn.loss_representation.FullFieldWavefront` oracle — agree on the
seed, which the two-paths bit-identity gate continues to verify.

**The four outflow capture-collapse DROP sites.**  Before #251, the sweep
captured the domain-edge outflow as a :math:`2^{d-1}`-moment object (the interior
cochain is moment-valued throughout the walk) and then collapsed it to slot 0
(``capture = tuple(c[..., AVERAGE_MOMENT] for c in capture)``) before writing
the scalar trace.  Four such collapses — guarded ``if n_face_moments > 1`` —
are removed so the capture RETAINS the :math:`2^{d-1}` axis:

.. list-table:: The four outflow capture-collapse DROP sites
   :header-rows: 1
   :widths: 22 22 56

   * - Path
     - Strategy
     - What it stored / now stores
   * - SOLVE
     - ``MovingFrontierWindow`` (prod)
     - was slot-0 average only → now the full :math:`2^{d-1}` outflow moments
   * - MATVEC
     - ``MovingFrontierWindow`` (prod)
     - same
   * - SOLVE
     - ``FullFieldWavefront`` (oracle)
     - same
   * - MATVEC
     - ``FullFieldWavefront`` (oracle)
     - same

The downstream sheds then land the moments automatically into the now
moment-shaped slot: the SOLVE shed writes ``boundary_flux.face_view(face)``; the
MATVEC shed writes into a ``streamed`` buffer allocated ``zeros_like`` of the
(now widened) ``face_view``, so it auto-widens; and the boundary-residual
``B``-block emit (the outflow defect ``streamed − given`` and the inflow
identity) is ordinate-indexed and spans the whole trailing slot.  The widened
layout makes these writes land the transverse moments with no further edit —
again the illegal-states-unrepresentable dividend of typing the slot shape once.

The producer — both ranks through one slot
''''''''''''''''''''''''''''''''''''''''''

The ergonomic generator for the affine-BC inhomogeneous term :math:`q` is
:meth:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink.prescribed_inflow`.
Its job is to write the inflow ordinate rows of each named face from a given
per-face array and leave everything else zero.  After the trace widened, this
producer must accept **two ranks through the same slot** — and getting this
right was where the architecture bit harder than the storage map anticipated.

.. list-table:: ``prescribed_inflow`` slot assignment (the two-arm relaxation)
   :header-rows: 1
   :widths: 38 62

   * - Supplied array shape
     - Write
   * - full slot — ``arr.shape == view.shape``
     - ``view[inflow] = arr[inflow]`` — the inflow ordinate ROWS (axis 0) span
       ALL trailing axes.  DD/Step's scalar slot AND a moment-resolved full slot
       both take this byte-identical arm.
   * - scalar onto a moment slot — ``arr.shape == view.shape[:-1]``
     - ``view[inflow, ..., AVERAGE_MOMENT] = arr[inflow]`` — seed the average
       moment, the transverse slopes stay zero (the scalar default).
   * - anything else
     - ``ValueError`` naming the expected full slot or the scalar reduction.

.. note:: **Audit the EXISTING scalar producers when you widen a slot, not only
   the new one.**

   The instinct is to fix only the producer line for the new moment-resolved
   caller — change the single trailing ``, :`` to span all trailing axes and
   trust the shape-check to validate the wider shape.  That is *necessary but not
   sufficient*.  The EXISTING scalar callers — the 2-D LD MMS and the 1-D
   prescribed-inflow MMS — supply a SCALAR
   :math:`(N, n_g, *\text{transverse})` array.  Once the LD slot grew to
   :math:`(N, n_g, *\text{transverse}, 2)`, an unconditional
   ``arr.shape != view.shape`` reject fires on *every* existing scalar LD caller
   (it reddened eight LD-2-D gates).  The fix is the SECOND arm above — a
   scalar-onto-a-moment-slot path that seeds slot 0.  This is the same class of
   under-scope as Leg A's field-space layer: a rigid scalar contract sitting
   ABOVE a widened slot needs a *typed-union relaxation*, not just an indexing
   fix.  When a carve widens a trace or field slot, the scalar callers feed the
   SAME widened slot — the producer must accept BOTH ranks.

The sharper Mode-10: structural teeth only, no improves-on-flat leg
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

Leg B is the textbook vv **Mode 10** trap (an *activated-but-unconstrained*
term) one notch sharper than Leg A.  The transverse face-slope is genuinely
exercised — projected, threaded, stored, reframed per octant, and consumed by
the cell update — yet its contribution to the converged near-boundary flux is
**sub-floor for ANY value claim**, not merely for its sign.  Probed live, seeding
the TRUE transverse slope makes the converged near-boundary L2 error against the
manufactured :math:`\phi` SLIGHTLY WORSE, and the FLIPPED slope slightly
BETTER:

.. list-table:: Why the converged near-boundary flux is sub-floor for the
   boundary slope
   :header-rows: 1
   :widths: 14 28 28 30

   * - ``nc``
     - flat (avg/centre) A-err
     - real-slope A-err
     - flipped-slope A-err
   * - 16
     - :math:`1.015\times10^{-2}`
     - :math:`1.030\times10^{-2}`
     - :math:`1.009\times10^{-2}` (improves? NO)
   * - 32
     - :math:`1.943\times10^{-3}`
     - :math:`1.966\times10^{-3}`
     - :math:`1.929\times10^{-3}` (improves? NO)

The reason is geometric: the boundary signal is too LOCALIZED and too SMALL — an
:math:`O(h)` intra-cell slope at the domain edge only — to register against the
bulk :math:`O(h^2)` discretization error that dominates the converged flux.
This is sharper than Leg A, where the bulk slope source carried real sub-cell
information across the WHOLE domain and "improves-on-flat" WAS achievable
(probed :math:`3.40\times10^{-3} < 5.99\times10^{-3}`).

**Consequence for the gate: there is NO converged-value-improvement leg.**  The
positive verification of Leg B is STRUCTURAL only, from the two places where the
sign / threading is :math:`O(1)`:

1. **The producer threads the projected slope through at machine precision** (the
   production-change proof).  ``test_ld_2d_boundary_slope_threaded_through_inflow_to_moments``
   builds a moment-resolved face from the ``leggauss`` projector, feeds it to the
   widened ``_inflow_to_moments``, and asserts (a) the producer RECOGNISES the
   moment-resolved input — it does NOT append a spurious second moment axis — and
   (b) ``np.testing.assert_array_equal`` on slot 1 holds EXACTLY.  A regression
   that re-zeroes the slope (the EXACT #251 bug) breaks this at :math:`O(1)`,
   where the converged flux would never catch it.

2. **A consumed transverse-slope sign flip moves the converged near-boundary
   flux** :math:`\gg` **the consumption tolerance** (the consumption proof).
   ``test_ld_2d_boundary_slope_sign_mutation_reddens``
   (``@verifies("ld-cartesian-2d")``) solves the same problem twice through the
   PUBLIC ``solve_sn_fixed_source`` — once with ``prescribed_inflow`` carrying
   :math:`+b_{\rm slope}` in slot 1, once with :math:`-b_{\rm slope}` — and
   asserts the near-boundary (edge-cell-masked) :math:`|\Delta\phi|/|\phi|`
   exceeds ``_CONSUMPTION_TOL = 1e-8``.  Probed:

   .. list-table::
      :header-rows: 1
      :widths: 14 32 32

      * - ``nc``
        - near-boundary :math:`|\Delta\phi|/|\phi|`
        - global :math:`|\Delta\phi|/|\phi|`
      * - 16
        - :math:`4.10\times10^{-3}`
        - :math:`3.27\times10^{-3}`
      * - 32
        - :math:`8.38\times10^{-4}`
        - :math:`8.99\times10^{-4}`

   The flip clears the tolerance by :math:`\sim 5.6` orders of magnitude and
   HALVES under refinement (:math:`O(h)`, boundary-localized), so this is a
   fixed-mesh consumption test, not a convergence leg.  Linearity is confirmed
   (seed slope :math:`= k\!\cdot\!\bar b`, :math:`k \in \{0.05, 0.1, 0.2\}` →
   flip :math:`|\Delta\phi|/|\phi| = \{2.4, 4.8, 9.5\}\times10^{-3}`, exactly
   linear in :math:`k`), proving the cochain GENUINELY consumes the transverse
   slope and the signal is not a numerical artifact.

These two teeth are PAIRED with the SCALAR-inflow no-op leg
(``test_ld_2d_boundary_scalar_inflow_no_op_negative_control`` and the byte-equal
``slope_sign=0`` vs ``slope_sign=None`` solve), which pins the Leg-B asymmetry:
a scalar inflow has slot 1 :math:`\equiv 0`, so flipping zero is a no-op and the
scalar path is correctly blind — and byte-identical to today's solve, the
bit-identity guard.

.. note:: **A current-invariant lesson worth recording (vv Mode 10 — the
   companion-unavailable branch).**

   Leg A established that a Mode-10 gap is closed by two :math:`O(1)` structural
   teeth (producer threads the projection at machine precision; a consumed
   source-row flip moves the converged answer :math:`\gg` solver tol), paired
   with a no-op control, and that the convergence-order leg is necessary for
   consistency but is not the sign teeth.  Leg B sharpens the recipe's *value*
   half: the canonical Mode-10 resolution suggests adding a companion gate "that
   isolates the term so its error is :math:`O(1)` in the measured quantity (e.g.
   a fixed-source problem where the term is the dominant forcing)".  **For a
   boundary-trace slope that companion is UNAVAILABLE** — there is no regime in
   which a localized :math:`O(h)`-small along-face boundary perturbation becomes
   the dominant forcing of the converged flux; it is intrinsically a boundary
   perturbation that rides below the bulk :math:`O(h^2)` floor.  In that case the
   producer-threading-at-machine-precision plus consumed-flip-:math:`\gg`-tol
   structural pair is the COMPLETE resolution — there is no value-improvement leg
   to add.  This is the first instance where the O(1)-isolating companion half of
   the Mode-10 recipe is genuinely unavailable, and it generalizes: whenever a
   term has no regime where it is the dominant forcing, structural teeth alone
   are the canonical close.

The negative pin, the bit-identity guard, and the reflective follow-up
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

**The negative pin.**  The moment relaxation must not swallow a real shape bug.
``test_ld_2d_boundary_trace_rejects_wrong_transverse_width`` feeds
``_inflow_to_moments`` a moment-resolved inflow with a trailing width
:math:`3 \neq 2^{d-1} = 2` on a 2-D LD mesh and asserts a ``ValueError`` naming
the expected :math:`2^{d-1}` (Pattern 4) — the relaxation admits exactly the two
principled ranks and nothing in between.

**The bit-identity guard.**  The DEFAULT scalar-inflow path must stay
bit-identical after the trace widening, verified three ways: the strict
``DriftWarning`` regression gate over ``tests/sn/sweep/core`` and
``tests/sn/solve`` (no golden moved — the moment-tail widening leaves the
DD/Step and scalar-inflow trace untouched, because ``face_moment_tail(1) == ()``
is the negative control); the existing 1-D prescribed-inflow MMS (byte-identical
by construction — a 1-D face has :math:`n_{\rm face} = 1`); and the explicit
scalar no-op leg above.

No ERR entry was minted: like Leg A, Mode 10 here is a proactive-gap close, not
a caught production bug — the trace correctly zeroed an unverified-but-honest
default (the transverse slope), which was UNVERIFIED, not WRONG.  Per the "log
every caught bug" directive, an ``@catches`` marker is added only when a real
production bug surfaces; none did.

.. note:: **Reflective-BC transverse-slope sign — a tracked follow-up (#252).**

   The Leg-B MMS is vacuum-BC (vacuum nulls the reflective coupling), so the
   reflective ``B`` operator's transverse-slope handling is exercised for
   STORAGE but not for its SIGN.  Storage is verified correct: the realized
   reflective law is a
   :class:`~orpheus.numerics.operator.PermutationOperator` on the ordinate axis
   (axis 0), which broadcasts UNCHANGED over the new trailing moment axis — no
   hard-coded trailing-axis assumption, so the widening introduces NO latent
   storage bug (read at
   :meth:`orpheus.sn.boundary_operator.SNBoundaryOperator._reflect_trace`, and
   confirmed empirically on a reflective-xmin LD-2-D mesh with a seeded slot 1).
   The SIGN, however, is UNVERIFIED.  Physics: a normal-flip reflection across a
   face preserves the tangent-plane (transverse) coordinate, so the transverse
   slope should reflect WITHOUT a sign flip (the permutation-on-axis-0
   pass-through is *probably* correct) — but this is a genuine Mode-1 sign trap
   the vacuum gates cannot see.  It is tracked in **#252**: a reflective-LD MMS
   leg plus an ``op.H`` adjoint check on the transverse-slope reflection.

Sources and gates
'''''''''''''''''

The production change spans three files: the storage lever in
:attr:`orpheus.sn.geometry.SNMesh.boundary_face_layout` (appending
:func:`~orpheus.numerics.moment_layout.face_moment_tail`); the inflow lift
:meth:`_LossRepresentation._inflow_to_moments`, the oracle seed
:meth:`FullFieldWavefront._octant_face_cochain`, and the four outflow
capture-collapse DROP sites in :mod:`orpheus.sn.loss_representation`; and the
producer
:meth:`~orpheus.transport.source_sinks.boundary_source_sink.BoundarySourceSink.prescribed_inflow`.
The end-to-end gates live in :mod:`tests.sn.verification.mms.test_mms_ld_2d` (the
#251 block): ``test_ld_2d_boundary_slope_threaded_through_inflow_to_moments`` (the
foundation structural teeth — the stamp / production-change proof),
``test_ld_2d_boundary_slope_sign_mutation_reddens`` (the L1 consumption proof,
``@verifies("ld-cartesian-2d")``),
``test_ld_2d_boundary_scalar_inflow_no_op_negative_control`` (the Leg-B
asymmetry), ``test_ld_2d_boundary_trace_rejects_wrong_transverse_width`` (the
negative pin), with the face-projection foundation sub-gates
``test_face_transverse_legendre_projection_matches_hand_polynomial`` and
``test_face_projection_slot0_is_transverse_cell_average``.  All #251 gates are
``-O``-safe (``np.testing.*`` / ``pytest.fail`` / ``pytest.raises`` only, no bare
``assert`` that ``python -O`` would strip — vv Mode 8).

.. _ld-cartesian-2d-coherent-promise:

S9 — the coherent boundary promise and the property-vs-type seam (#257)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Leg B (:ref:`ld-cartesian-2d-legB`, #251) landed the boundary CONSUMPTION
path: the trace ``mesh.trace`` carries the transverse moment axis end to end,
``_inflow_to_moments`` threads slot 1, and the sweep outflow stores the
:math:`2^{d-1}` transverse face-moments instead of collapsing them.  But the
MMS *producer* — the case's
:meth:`~orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase.prescribed_inflow`
— still built a SCALAR per-face trace, so it hit the producer's scalar branch
and the slope it could have carried was zeroed by the case, not by the closure.
S9 (Issue #257) closes that producer-blindness — the case now EMITS the
moment-resolved slot — and, in doing so, LOCKS the **coherent promise** the
whole LD-on-the-boundary effort is about with a dedicated convergence gate.

This subsection records three things a future session must not re-derive: (1)
WHY the coherent promise is already TRUE and what delivers it (the average
moment, not the slope); (2) the producer-honesty carve and why it is
byte-identical for DD/Step; and (3) the **property-vs-type** seam (#263) — why
the transverse boundary moment is a PROPERTY of the boundary field and NOT a
new first-class field type, the durable design invariant that bounds the S9
scope.

The coherent promise: LD is second-order at the boundary, no asterisk
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The motivating claim — *"LD gives second-order accuracy EVERYWHERE, including
the boundary, with no 'but not at the boundary' asterisk"* — is **TRUE and
already DELIVERED**.  The subtle and load-bearing point is WHAT delivers it.

It is delivered by the **AVERAGE moment alone**, not by the transverse
boundary-slope moment.  The prescribed inflow is exact at the face cells (the
manufactured trace is evaluated there), and the bulk LD closure
(:eq:`ld-cartesian-2d`) carries that boundary data inward at :math:`O(h^2)`.
The cell that directly consumes the inflow integrates it against the cell's
own linear basis, and the cell-AVERAGE transverse moment is exactly what that
integral needs to :math:`O(h^2)`.  So the boundary cell is second-order
*from the average moment*, before any transverse-slope refinement enters.

This reframes the S9 motivation.  There was never an :math:`O(h)` deficiency in
the converged flux at the boundary to "recover" — the order is already
:math:`O(h^2)` there.  What the transverse boundary-slope moment improves is the
*inflow representation*: it lifts the face trace from :math:`O(h)`-accurate
(bar-only) to :math:`O(h^2)`-accurate (bar + slope), a genuine refinement that
the LD closure genuinely consumes (Leg B's structural teeth prove this).  But a
second-order correction to an already-second-order face balance cannot move the
converged flux above the bulk :math:`O(h^2)` floor.  S9 therefore does NOT
remove an asterisk on the convergence order (there was none); it makes the
producer honest about the moment it could always have carried, and it LOCKS the
no-asterisk promise so a future change to the boundary closure cannot silently
break it.

First-cell-row evidence: the boundary cell is already O(h²)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

The decisive scaling argument is the convergence order of the FIRST CELL ROW —
the cells (:math:`i=0`) that directly consume the ``xmin`` inflow.  If the
average-only inflow were :math:`O(h)`-deficient at the boundary cell, the slope
COULD repair it and the coherent promise would carry an asterisk.  It is in
fact :math:`O(h^2)`, for BOTH the average-only (``flat``) and the
moment-resolved (``mom``) inflow:

.. list-table:: First-cell-row (:math:`i=0`) :math:`L^2` order — already
   :math:`O(h^2)` from the average moment (pure-absorber streaming,
   :math:`\Sigma_t L = 0.1`, the regime most favourable to a boundary-confined
   inflow error)
   :header-rows: 1
   :widths: 24 38 38

   * - inflow treatment
     - first-cell-row :math:`L^2` (``nc`` = 16, 32, 64, 128)
     - observed order
   * - ``flat`` (average only)
     - sequence decays at the design rate
     - :math:`1.993,\ 2.004,\ 2.001`
   * - ``mom`` (average + transverse slope)
     - sequence decays at the design rate
     - :math:`1.998,\ 2.005,\ 2.003`

The orders are indistinguishable: the slope is a sub-floor refinement, not a
deficiency repair.  This is the
:func:`~tests.sn.verification.mms.test_ld_2d_boundary_promise.test_first_cell_row_already_second_order`
gate (``@l1``, ``@verifies("ld-cartesian-2d")`` — the label is REUSED, S9
mints none): it fails iff the first-cell-row flat order drops below
:math:`1.85`, which would mean the average inflow is :math:`O(h)`-deficient at
the boundary cell and the verdict must flip.

This holds across the full optical-depth axis, NOT only at the cheap
:math:`\Sigma_t L = 0.1` headline.  Probed across
:math:`\Sigma_t L \in \{0.1, 0.5, 1.0, 2.0\}` (streaming → thick) and
:math:`c \in \{0, 0.5\}`, ``flat``, ``mom``, and ``flip`` all converge in the
second-order band globally, and the magnitudes track within a documented
sub-floor band (:math:`<20\%`, with a :math:`30\%` guard).  Even in the
streaming limit (where the inflow propagates ballistically and an
inflow-representation error is LEAST boundary-confined) the LD cell balance
integrates the inflow against its own linear basis, so the average moment is
:math:`O(h^2)`-adequate.  Amplifying the :math:`\mu`-dependent (slope-carrying)
drivers up to :math:`20\times` — the user's strongest "make the boundary slope
non-trivial" hypothesis — makes the converged flux MONOTONICALLY WORSE
(:math:`+17\%` at :math:`20\times`, still :math:`O(h^2)`), never better.  The
sub-floor wall is FUNDAMENTAL to a boundary-trace moment, not an artefact of
the cheap regime; these are the
:func:`~tests.sn.verification.mms.test_ld_2d_boundary_promise.test_optical_sweep_slope_never_beats_floor`
and
:func:`~tests.sn.verification.mms.test_ld_2d_boundary_promise.test_amplified_boundary_slope_still_subfloor`
verdict pins (``@slow`` ``@l1``).

Producer honesty: the MMS case emits the moment slot
''''''''''''''''''''''''''''''''''''''''''''''''''''

The S9 production change is a single completion of the existing case (no new
ansatz, no sibling case).
:meth:`~orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase.prescribed_inflow`
now gates on
:func:`~orpheus.numerics.moment_layout.face_moment_count` — the SAME
single-source primitive
:attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout` keys the slot width on:

* When ``face_moment_count == 1`` (DD/Step) it builds the SCALAR per-face trace
  ``(N, ng, n_t)`` by cell-CENTRE evaluation of :math:`(A + \mu_x B + \mu_y
  C)/W`, hits the producer's scalar branch, and is **byte-identical** to the
  pre-S9 build (proven: ``np.array_equal`` against the legacy ``face_coords``
  construction holds — DD/Step has no moment axis, the slots are bit-for-bit
  the old build).
* When ``face_moment_count > 1`` (LD) it builds the FULL moment slot
  ``(N, ng, n_t, face_moment_count)`` via a new case-owned
  :meth:`~orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase._project_inflow_to_face_moments`
  and hits the producer's full-slot branch.

The projector descends ONLY from the case's manufactured harmonics
(``_drivers``) and :func:`numpy.polynomial.legendre.leggauss` — NEVER
``_inflow_to_moments``, ``_ubld``, any LD operator, or the test-side
projectors (the L11 structural-independence discipline of the
:doc:`algebra-of-record </theory/reference_solvers>` pillar).  Per transverse
cell :math:`[t_L, t_R]` mapped to :math:`\xi \in [-1, 1]`, it projects onto the
BARE per-cell Legendre coefficients

.. math::

   \text{slot}_0 \;=\; \frac{\langle \psi, P_0\rangle}{\langle P_0, P_0\rangle}
   \quad(\text{transverse cell AVERAGE}),
   \qquad
   \text{slot}_1 \;=\; \frac{\langle \psi, P_1\rangle}{\langle P_1, P_1\rangle}
   \quad(\text{BARE transverse slope}).

This reuses the exact normalization Leg B locked
(Eq. :eq:`ld-cartesian-2d-face-projection-coeff`): NO :math:`\theta`/:math:`h_t`
weighting, because the cochain's transverse mass :math:`\mathrm{diag}(h_t,
\theta h_t)` applies them downstream — a :math:`\theta`- or :math:`h_t`-weighted
slope would double-apply the mass, a TRUE bug.  Because the case projector and
the test-side ``_face_transverse_legendre`` are deliberately INDEPENDENT
implementations of the same projection, their machine-precision agreement is a
single-source check (Cardinal Rule 2), pinned by the new foundation gate
:func:`~tests.sn.verification.mms.test_mms_ld_2d.test_case_projector_agrees_with_test_face_projector`;
a shared import would make that check tautological and let a double-applied mass
slip through.  The threading is then pinned end to end by the GATE-B leg added
to ``test_ld_2d_boundary_slope_threaded_through_inflow_to_moments`` (the
production producer's slot 1 equals the ``leggauss`` reference at machine
precision — closing the Mode-11 producer-blindness the #251 surrogate had, which
pinned only the consumer).

.. note:: **The second-order interaction the producer honesty introduces.**

   A diagnostic helper whose "average-only" baseline routed through the (now
   slope-honest) production ``prescribed_inflow`` would SILENTLY inherit the new
   slope, collapsing the controlled ``flat``/``mom``/``flip`` toggle and breaking
   the byte-identity no-op control.  The fix is structural: the controlled
   average-only baseline is built TEST-SIDE (the toggle belongs to the test;
   the honesty belongs to production).  This is a general lesson — when a
   refactor makes a production path honest about a quantity a test was
   controlling externally, the test's controlled toggle MUST NOT inherit the
   production change it is meant to be orthogonal to.

The fourth Mode-10 instance — and why no value gate is added
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

S9 is the FOURTH consecutive vv **Mode 10** close-out in the
slope-moment family — #240 D5b-S4 (the bulk slope-source) →
:ref:`ld-cartesian-2d-legA` (#247, the external :math:`\hat Q`) →
:ref:`ld-cartesian-2d-legB` (#251, the boundary transverse slope) → S9 — and it
sits squarely in the **companion-unavailable branch** Leg B opened.  The
transverse boundary-slope moment is *activated-but-unconstrained*: its code path
is genuinely exercised (projected, threaded, stored, reframed per octant,
consumed by the cell update — the Mode-11 sentinel
:func:`~tests.sn.verification.mms.test_ld_2d_boundary_promise.test_slope_toggle_reaches_inflow_to_moments`
confirms slot 1 reaches the production consumer and the converged flux differs,
so the toggle is non-vacuous), yet its contribution to the converged flux is
sub-floor for ANY value claim.

The canonical Mode-10 resolution says: *add a companion gate that isolates the
term so its error is :math:`O(1)` in the measured quantity (a problem where the
term is the dominant forcing)*.  **For a boundary-trace slope that companion is
genuinely UNAVAILABLE.**  The boundary is codimension-1 — measure-zero in the
refinement limit — so a localized :math:`O(h)`-small along-face perturbation has
NO regime in which it becomes the dominant forcing of the converged flux; the
optical-depth and amplitude sweeps above are the empirical proof that no such
regime exists.  In that case the structural pair is the COMPLETE resolution, and
manufacturing a value-improvement leg would be dishonest — it would falsely RED a
correctly-consumed term (probed: seeding the TRUE slope makes the near-boundary
error slightly WORSE, the flipped slope slightly BETTER; both sub-floor).  So
S9 keeps Leg B's structural teeth (machine-precision threading + consumed-flip
:math:`\gg` tolerance + the byte-identical scalar no-op control) and adds NO
value or order gate keyed on the slope.

.. warning::

   Do NOT write "S9 recovers second-order accuracy at the boundary" or
   "the boundary-slope moment makes LD second-order at the boundary".  Both are
   false: the boundary cell is :math:`O(h^2)` from the AVERAGE moment alone, and
   the slope is a sub-floor inflow-representation refinement.  The verifiable
   content of S9 is STRUCTURAL (producer threading + consumption + the
   coherent-promise lock), per the vv Mode-10 companion-unavailable branch — not
   a converged-value claim.

The property-vs-type seam (#263): a boundary moment is a PROPERTY
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

S9 raised a natural design question: if the angular moment representation
(:class:`~orpheus.transport.fields.harmonic_moment_field.HarmonicMomentField`)
is a first-class FIELD TYPE, should the transverse boundary moment be one too —
a ``BoundaryMomentField``?  The answer, tracked in **#263**, is **NO today**:
the transverse boundary moment is a PROPERTY of the boundary field (an untyped
trailing moment axis on the flat face buffer), exactly as the bulk carries its
spatial moments as a property (``spatial_moments`` on the field, rather than a
distinct field type).  The criterion and its trigger live in the
:ref:`field-type-vs-property-criterion` section of the operator-algebra page;
the short version, specialised to the boundary moment:

A representation earns a distinct first-class **type** only when a
**non-canonical dual** coexists with it — two bases that are NOT canonically
isomorphic, connected by a change-of-basis morphism that is itself modelled and
applied (carries truncation error, has an adjoint, participates in the operator
algebra).  Angular order PASSES: the ordinate basis (``AngularFlux``) and the
harmonic-modal basis (``HarmonicMomentField``) are non-canonically isomorphic
(the iso depends on the quadrature :math:`Y_\ell^m(\hat\Omega_n)`), bridged by
the applied :math:`M`/:math:`R` projection/reconstruction pair with truncation
content and adjoints.  The transverse SPATIAL moment FAILS: there is ONE
within-cell basis (the tensor-Legendre tower), no non-canonical dual, so the
only change-of-basis would be the identity.  A ``BoundaryMomentField`` leaf
whose ``_check_partner`` adds nothing beyond class identity would be a vacuous
naming leaf — type-theatrics by the project's own "if the type hint does not
prevent a bug by construction it is theatrics" standard.  So the moment rides as
a PROPERTY (the flat face buffer already holds the moment tail via
:attr:`~orpheus.sn.geometry.SNMesh.boundary_face_layout`), and the first-class
:class:`~orpheus.numerics.spaces.spatial_moment_space.SpatialMomentSpace` field
type is DEFERRED to the collocation trigger (nodal-DG / Lagrange-FEM, where a
nodal point-value basis coexists with the modal coefficients and a Vandermonde
morphism bridges them) — the durable design invariant #263 records.

Sources and gates
'''''''''''''''''

The production change is one file:
:meth:`~orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase.prescribed_inflow`
gated scalar-vs-moment build plus the new case-owned ``leggauss`` projector
:meth:`~orpheus.derivations.continuous.mms.sn.SN2DCartesianLDStressMMSCase._project_inflow_to_face_moments`
(the boundary cochain, the trace, the consumer ``_inflow_to_moments``, and
DD/Step are all UNCHANGED — Leg B already landed the consumption path).  The
coherent-promise gate and the verdict pins live in
:mod:`tests.sn.verification.mms.test_ld_2d_boundary_promise`:
``test_first_cell_row_already_second_order`` (``@l1``
``@verifies("ld-cartesian-2d")`` — the coherent-promise lock),
``test_optical_sweep_slope_never_beats_floor`` and
``test_amplified_boundary_slope_still_subfloor`` (``@slow`` ``@l1`` — the
sub-floor verdict pins guarding the no-value-gate conclusion across optical
depth and amplitude), and ``test_slope_toggle_reaches_inflow_to_moments``
(``@foundation`` — the Mode-11 sentinel proving the toggle is non-vacuous).  The
producer-stamp leg lives in
:mod:`tests.sn.verification.mms.test_mms_ld_2d`:
``test_ld_2d_boundary_slope_threaded_through_inflow_to_moments`` (GATE B — the
production producer's slot 1 equals the ``leggauss`` reference) and
``test_case_projector_agrees_with_test_face_projector`` (GATE C — the
single-source projector agreement).  All gates are ``-O``-safe
(``np.testing.*`` / ``pytest.fail`` only, no bare ``assert`` that ``python -O``
would strip — vv Mode 8).  No ERR entry was minted: like Leg A and Leg B, S9 is
a proactive-gap close (a producer-blindness — the slope was UNVERIFIED at the
producer, not WRONG — the #251 consumer already threaded it correctly), not a
caught production bug.

.. _sweep-wavefront:

Cartesian 2D: Anti-Diagonal Wavefront Sweep
---------------------------------------------

In 2D, the DD equation :eq:`dd-cartesian-2d` creates a data dependency:
cell :math:`(i, j)` requires incoming face fluxes from its upwind
neighbours in both :math:`x` and :math:`y`.  Cells along an
**anti-diagonal** :math:`i + j = k` are mutually independent because
they share no incoming faces, so they can be solved simultaneously.

The wavefront sweep is implemented as a **per-octant batched
forward-substitution** over a precomputed causal cell DAG (Wave 2 of
the SN performance plan, closing Issue #4).  This subsection states
the algebraic framing; the primitives that realise it
(:class:`~orpheus.sn.sweep_graph.OctantLabel`,
:class:`~orpheus.sn.sweep_graph.SweepDependencyGraph` and its two
storage walks, the level-operation pair ``_CellSolve`` /
``_CellResidual``, and the discretization's kernel pair
:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
/ :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`)
are documented in detail at
:ref:`sweep-octant-dependency-graph` immediately below.

The §15A.2 sum-of-tensor-products framing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following Grand Report v3 §15A.2 (lines 2137–2171), the streaming
inverse :math:`L^{-1}` on the 2-D Cartesian SN field space decomposes
as a **direct sum over angular octants**:

.. math::
   :label: streaming-inverse-direct-sum

   L^{-1} \;=\; \bigoplus_{\sigma \in \mathcal{O}} L^{-1}_{\sigma},
   \qquad
   \mathcal{O} \;=\; \{\sigma = (\mathrm{sgn}\,\mu_x,\,
                                  \mathrm{sgn}\,\mu_y) :
                       \sigma \neq (0,0)\}
                  \,\cup\, \{(0,0)\},

.. vv-status: streaming-inverse-direct-sum documented

acting on the octant-restricted tensor space :math:`(N_\sigma,\,n_x,\,n_y,\,n_g)`.
The direction-cosine partition (Eq. :eq:`octant-sign-predicate`) is
the predicate the four
:class:`~orpheus.sn.quadrature.AngularQuadrature` adapters expose as
their cached :attr:`~orpheus.sn.quadrature.AngularQuadrature.octants`
property — a tuple of
:class:`~orpheus.numerics.measure.DiscreteMeasurePartition`
entries realised by
:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
(see :ref:`tensorial-framing` and the
:doc:`discrete_measures` consumer table).

For each non-degenerate octant :math:`\sigma`, the action of
:math:`L^{-1}_\sigma` is a **forward substitution along a per-octant
causal cell DAG** — the topological order is structural (anti-diagonal
sweep on the Cartesian grid), the per-level cell update is one
vectorised einsum.  The pure-:math:`z` degenerate octant
:math:`\sigma = (0,0)` (ordinates with :math:`\mu_x = \mu_y = 0`,
which appear in 3-D angular cubatures projected to the in-plane
2-D problem) has no spatial streaming and reduces to a per-cell
balance :math:`\psi = Q / \Sigma_t` — the wavefront sweep handles
it via a short-circuit and skips the dependency graph.

**The four quadrant sweeps.**  Each non-degenerate octant
:math:`\sigma = (\mathrm{sgn}\,\mu_x, \mathrm{sgn}\,\mu_y) \in \{-1,+1\}^2`
determines a sweep direction:

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

For each octant, the sweep visits topological levels
(anti-diagonals) :math:`k = 0, 1, \ldots, n_x + n_y - 2`.  On level
:math:`k`, the cells :math:`(i, j)` satisfying :math:`i + j = k`
(in the per-octant traversal index space) are gathered into a numpy
batch and solved with a single vectorised evaluation of
:eq:`dd-cartesian-2d` — vectorised across the **ordinate axis**
(:math:`N_\sigma` — every ordinate in the octant), the
**anti-diagonal axis** (:math:`n_{\rm diag}` — number of cells on
this level), and the **group axis** (:math:`n_g`) simultaneously.

**Vectorisation within each level.**  Each level contains up to
:math:`\min(n_x, n_y)` cells.  The incoming face fluxes
``psi_in_x`` and ``psi_in_y`` are gathered by advanced indexing;
the DD equation is evaluated as one numpy operation; and the
outgoing face fluxes are scattered back into the persistent face-
flux buffers.  There is **no Python-level cell loop within a level**
and **no Python-level ordinate loop within an octant** — both axes
are internal to the einsum.

**Reflective BCs in 2D.**  At each boundary face, the incoming flux
for ordinate :math:`n` is set to the outgoing flux of its reflected
partner.  For the left/right boundaries (*x*-reflection), the partner
is ``ref_x[n]`` (negating :math:`\mu_x`); for the top/bottom boundaries
(*y*-reflection), the partner is ``ref_y[n]`` (negating :math:`\mu_y`).
The reflection indices are precomputed by the quadrature's
:meth:`reflection_index` method.  Crucially, the BC apply happens
**once per octant per axis** (not once per ordinate per axis) —
see :ref:`sweep-octant-dependency-graph-l7-trap` for the rationale.

Implemented in :func:`~orpheus.sn.loss_representation._sweep_jacobi`, which
is a thin orchestrator over the
:class:`~orpheus.sn.sweep_graph.SweepDependencyGraph` primitives
described next.

.. _sweep-octant-dependency-graph:

Cartesian 2D: Octant Dependency Graph (Wave 2)
-----------------------------------------------

This section documents the **§15A.2 "upwind trace complex / causal
transport DAG / direction sweep ordering" primitive** as it lives in
:mod:`orpheus.sn.sweep_graph` after Wave 2 of the SN performance plan
(branch ``feature/sn-octant-sweep-graph``, closes Issue #4).  The
shipped architecture replaces the legacy per-ordinate ``for n in
range(N)`` loop in :func:`~orpheus.sn.loss_representation._sweep_jacobi` with
a per-octant batched dispatch, lifting the per-call ``_diag_cache``
build to mesh-time work, and isolating the per-cell DD algebra in the
discretization's pure kernel pair
(:meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`
/ :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`)
that LD / EC / Step closures can override later.

.. note::

   **Architecture history — the dispatch surface re-layered twice.**
   Wave 2 (the original closure of Issue #4) routed the sweep through
   a per-level *packet* (the ``SweepCellSlice`` dataclass) consumed by
   four direction×storage methods — ``update_batch`` / ``residual_batch``
   on the strategy (full-field) plus their ``apply_windowed`` /
   ``residual_windowed`` siblings on the graph.  S6.4(e) **collapsed
   that surface**: the four walk methods became TWO storage walks
   (:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full`,
   :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`)
   each parameterised by a level-operation OBJECT (``_CellSolve`` for
   the solve direction, ``_CellResidual`` for the apply direction —
   direction is never a boolean flag); the per-level ``SweepCellSlice``
   packet was retired (it existed only to feed the now-deleted storage
   adapters); and the strategy's ``update_batch`` / ``residual_batch``
   were replaced by the **storage-free kernel pair**
   ``cell_kernel_batch`` / ``residual_kernel_batch`` (pure cell
   algebra — no gather/scatter).  The historical names ``update_batch``
   / ``residual_batch`` / ``SweepCellSlice`` appear below only as
   *history*; the current contract is the kernel pair + the level
   operations.  See :ref:`sweep-dispatch-relayering` for the WHY.

The primitives
~~~~~~~~~~~~~~

The architecture is a small set of frozen, individually unit-tested
primitives plus a mesh-time precompute step.

.. list-table::
   :header-rows: 1
   :widths: 28 16 56

   * - Primitive
     - Lives in
     - Role
   * - :class:`~orpheus.sn.sweep_graph.OctantLabel`
     - :mod:`orpheus.sn.sweep_graph`
     - Frozen + slotted dataclass carrying one direction sign per
       spatial axis (``signs[axis] ∈ {-1, 0, +1}``) — a single type
       labels a 1-D (``(±1,)``), 2-D (``(±1, ±1)``), or 3-D octant.
       Hashable; used as the key in the per-shape graph family
       :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.for_shape`
       (owned by the ``_DAGWavefront`` representation family since
       S6.4(c) — historically a mesh attribute).  An all-zero
       signature denotes the pure-:math:`z` degenerate octant — no
       graph is built for it
       (:attr:`~orpheus.sn.sweep_graph.OctantLabel.streams` is
       ``False``).  The 3-D ``sign_z`` is dropped by the 2-D Cartesian
       orchestration: the in-plane sweep is invariant under the
       out-of-plane axis, so multiple ordinates with the same in-plane
       ``signs`` but different ``sign_z`` share a single graph instance.
   * - :class:`~orpheus.sn.sweep_graph.SweepDependencyGraph` (+ its
       two storage walks)
     - :mod:`orpheus.sn.sweep_graph`
     - Frozen dataclass holding the per-octant topological levels
       (anti-diagonals) and the per-axis face-index offsets.  Built
       once per ``(shape, octant)`` pair in the
       :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.for_shape`
       cache (S6.4(c); historically at mesh construction); reused
       across every source iteration / Krylov matvec / outer
       iteration.  Exposes TWO storage walks (S6.4(e)):
       :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full`
       carries the COMPLETE per-axis interior face cochain (the
       verification-oracle policy);
       :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`
       advances a rolling :math:`(d{-}1)`-frontier window (the
       production policy, ``O(N·n_g·∏ n_a)`` shrunk to
       ``O(N·n_g·∏_{a<d−1} n_a)`` backing).  The walk owns the level
       loop, the storage, and the per-level operand extraction; it
       dispatches the cell algebra to a level operation (next two rows).
   * - The level-operation pair ``_CellSolve`` / ``_CellResidual``
     - :mod:`orpheus.sn.sweep_graph`
     - The **direction fork, as OBJECTS** (S6.4(e); direction is never
       a boolean flag).  Exactly ONE is constructed per octant walk; the
       storage walk calls ``level_op.cell(...)`` per topological level.
       ``_CellSolve`` runs the solve direction — calls the strategy's
       ``cell_kernel_batch`` then performs the Phase-5c angular-XOR-
       moment per-level emit (write the angular flux + accumulate the
       scalar flux, OR accumulate the harmonic-moment tensor, never
       both).  ``_CellResidual`` runs the apply direction — calls
       ``residual_kernel_batch`` then writes the per-level residual.
       The per-level *emit* expressions and their order are
       bit-identity-load-bearing — relocated verbatim from the four
       retired walk methods.
   * - The kernel pair
       :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
       / :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
     - :mod:`orpheus.sn.spatial.scheme`
     - The **storage-free extension point** on the
       :class:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase` ABC
       (S6.4(e); historically the ``SweepCellSlice``-packeted
       ``update_batch`` / ``residual_batch``).  Each takes the per-axis
       incoming face fluxes + streaming coefficients + the level's cross
       section and source and returns ``(psi_avg, psi_out)`` (solve) or
       ``(residual, psi_out)`` (apply) — PURE cell algebra, no
       gather/scatter (that is the walk's job).  Default raises
       :exc:`NotImplementedError` — additive capability, not a contract
       change.  :class:`~orpheus.sn.spatial.diamond.DiamondDifference`
       overrides the pair; LD / EC / Step closures override it later to
       join the batched wavefront walks (their per-cell
       :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.update`
       stays the canonical reference contract).

Per-shape precompute pattern (family-owned since S6.4(c))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dependency graph is a **derived object** — the
``(shape × octant)`` joint property.  It depends only on cell topology
and the octant sign convention; it does **not** depend on fluxes,
sources, BCs, quadrature, cross sections, or iteration state.  So the
graph build is paid once **per spatial shape** in the cached accessor
:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.for_shape`, owned
by the DAG-consuming ``_DAGWavefront`` representation family:

.. code-block:: python

   class _DAGWavefront(_LossRepresentation):
       @property
       def sweep_graphs(self) -> dict[OctantLabel, SweepDependencyGraph]:
           # cached per shape; same-shape meshes share byte-identical graphs
           return SweepDependencyGraph.for_shape(self.mesh.spatial_shape)

**Ownership history** (two relocations, each a pure refactor):

#. *Wave 2 / C2.4* lifted the per-call ``_diag_cache`` build that
   previously lived inside the 2-D wavefront sweep (rebuilt once per
   sweep call) to **mesh-construction** time — a measurable but
   second-order saving on the 421-group benchmark; the structurally
   important effect was making the graphs named, inspectable state.
#. *S6.4(c)* moved ownership **off the mesh onto the representation
   family**: the mesh is pure geometry, and only the two DAG-walking
   representations (the window + the full-field oracle) ever mention
   the substrate.  This retired the curvilinear
   ``mesh.sweep_graphs = None`` slot — an illegal state (a mesh
   carrying a "no DAG here" marker for a structure it never owned) —
   and replaced mesh-lifetime caching with per-SHAPE caching, so
   same-shape meshes share one graph family (the graphs carry no
   mesh-identity information).  DAG-free representations
   (``CumprodScan``, ``ScanMarch``) and curvilinear meshes simply
   never touch the accessor; curvilinear sweeps walk the cell graph
   differently (per-ordinate march; see
   :meth:`~orpheus.sn.geometry.SNMesh.dag_walk`).

The closed-form precompute lives in
:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.from_cartesian`
and never appears in the sweep loop.  This is structural, not
hand-rolled — the "library version" (a generic topological-sort over
an explicit DAG) would be over-engineering for a regular pattern that
collapses to ~5 lines of ``arange`` + anti-hyperplane extraction.  The
builder is dimension-generic (``d = len(shape) ∈ {1, 2, 3}``): a d=1
chain, the d=2 anti-diagonal, a d=3 anti-hyperplane.

The §15A.2 invariant set
~~~~~~~~~~~~~~~~~~~~~~~~~

The Grand Report v3 §15A.2 (lines 2165–2171) prescribes a fixed set
of L0 invariants every ``SweepDependencyGraph`` instance must satisfy.
These are pinned by ``tests/sn/test_sweep_graph.py`` (63 L0 tests):

* **Upwind orientation** — for each octant
  :math:`\sigma = (\mathrm{sgn}\,\mu_x, \mathrm{sgn}\,\mu_y)`, the
  ``face_in_x`` and ``face_out_x`` offsets satisfy
  ``face_in_x + face_out_x == 1`` and
  ``face_in_x = 0`` iff :math:`\mathrm{sgn}\,\mu_x \ge 0` (and
  analogously on :math:`y`).  Asserted by
  ``test_face_pairing_consistent`` and ``test_upwind_orientation``.
* **Topological sort** — every level's cells depend only on cells in
  strictly earlier levels (under the per-octant orientation).  No
  intra-level dependencies; no back-edges.  Asserted by
  ``test_topologically_sorted``.
* **Cell coverage** — every cell :math:`(i, j) \in [0, n_x) \times
  [0, n_y)` appears in **exactly one** level.  Disjoint union over
  the topological levels reconstructs the full grid.  Asserted by
  ``test_cell_coverage``.
* **Face-pairing consistency** — the incoming-face index of cell
  :math:`(i, j)` on level :math:`k` matches the outgoing-face index
  of its upwind neighbour on level :math:`k - 1` (under the per-
  octant orientation).  Asserted by
  ``test_face_pairing_consistent``.

These four invariants are the **load-bearing correctness floor** of
the wavefront sweep.  Any future closure (LD, EC, Step) plugged in
via the kernel pair consumes the same invariants — they describe
the topology, not the algebra, so the strategy contract is orthogonal
to the graph correctness.

.. _sweep-dispatch-relayering:

The dispatch boundary: walk (scheduler) vs cell update (closure)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A central architectural decision is the **separation between the
scheduler and the closure**.  Three layers stack from storage outward
to algebra (S6.4(e)):

#. **The storage walk** —
   :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full` (full
   cochain) or
   :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`
   (rolling frontier).  Owns the topological-level loop and the
   per-axis face gather/scatter (full cochain) or the frontier
   seed/incoming/emit/shed cochain trace algebra (window).  Storage is
   the walk's concern — the SAME two walks serve every closure and
   both directions.

#. **The level operation** — ``_CellSolve`` or ``_CellResidual``,
   constructed once per octant walk and called as ``level_op.cell(...)``
   per level.  Owns the direction fork (solve vs apply) and the
   per-level *emit* (angular/moment write, or residual write).
   Direction is an OBJECT, never a boolean flag passed down the walk.

#. **The kernel pair** —
   :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.cell_kernel_batch`
   /
   :meth:`~orpheus.sn.spatial.diamond.DiamondDifference.residual_kernel_batch`.
   Owns the **pure cell algebra** and nothing else — no gather, no
   scatter, no storage. This is the ONLY direction-aware math left in
   the SN spatial stack.

**Why this layering (the WHY behind S6.4(e)).**  Wave 2 carried the
storage concern *inside* the strategy: the DD ``update_batch`` /
``residual_batch`` methods gathered the cell's face inputs from the
``SweepCellSlice`` packet, ran the algebra, and scattered the outgoing
faces back — a four-method direction×storage product
(``update_batch`` / ``residual_batch`` full-field +
``apply_windowed`` / ``residual_windowed`` windowed).  That entangled
two orthogonal concerns: a NEW closure (LD / EC / Step) would have had
to re-implement the gather/scatter (storage) plumbing four times just
to supply its cell math.  S6.4(e) lifts storage to the walk layer
**once, above every strategy**, so a closure supplies ONLY its
storage-free kernel pair (pure algebra over the per-axis incoming face
fluxes) and inherits both storage policies (full + window) and both
directions (solve + apply) for free.  The ``SweepCellSlice`` packet —
which existed only to feed the retired storage adapters — is gone with
them.  This is the Cardinal-Rule-2 "build primitives, not products"
discipline: the four-method product collapses to a 2 (walks) × 1
(level-op pair, direction-by-object) × 1 (kernel pair) factoring where
each factor varies independently.

This means: **DD is the only shipping closure today**, but Step / LD
/ EC override the kernel pair later without touching the walk driver
or the level operations.  The Wave C-extension rollout (Issues #157 /
#158) ships the per-cell :meth:`update` method first as the canonical
reference contract; the batched kernel pair is the parallel
level-vectorised capability for closures whose per-cell algebra
vectorises across an ``(N_oct, n_diag, ng)`` slice without per-cell
branching.

The DD ``cell_kernel_batch`` reproduces the legacy 2-D wavefront DD
math **bit-identically** (operation order matters; see
:class:`~orpheus.sn.spatial.diamond.DiamondDifference` docstring on
bit-identity).  The math is the **balance form** of WDD on a 2-D
Cartesian cell:

.. math::
   :label: dd-2d-balance-form

   \overline{\psi}_{i,j}
   \;=\; \frac{Q_{i,j}
               + s_{x,i}\,\psi^{\rm in}_{x,i,j}
               + s_{y,j}\,\psi^{\rm in}_{y,i,j}}
              {\Sigma_{t,i,j} + s_{x,i} + s_{y,j}},
   \qquad
   s_{x,i} = \frac{2|\mu_x|}{\Delta x_i},
   \quad s_{y,j} = \frac{2|\mu_y|}{\Delta y_j},

.. vv-status: dd-2d-balance-form documented

with the spatial closure
:math:`\psi^{\rm out}_x = 2\overline{\psi} - \psi^{\rm in}_x`
(and analogously on :math:`y`).  The operation order is fixed at:

.. code-block:: text

   denom    = sig_t + sx + sy
   psi_avg  = (Q + sx * psi_in_x + sy * psi_in_y) / denom
   psi_out  = 2 * psi_avg - psi_in

Algebraically equivalent rearrangements (e.g., reordering
``sig_t + sx + sy`` to ``sx + sy + sig_t``) break the 1-ULP regression
contract even though the math is identical.  This is the canonical
**bit-identity vs principled-equivalence** instance from the
``vv-principles`` skill — the regression contract is bit-identity
gated on the existing snapshots; deviations are admissible only when
backed by a structurally-independent reference (e.g., :math:`k_\infty`
analytical limit on a homogeneous reflective problem).

.. _sweep-octant-dependency-graph-l7-trap:

The L7-trap fix: BC apply once per octant
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wave 2 closes a class of bugs that the test-architect dispatch
identified as the **L7 trap** — the design pattern where a sweep
driver re-applies a boundary operator at each ordinate iteration.
The legacy ``_sweep_jacobi`` had this shape:

.. code-block:: python

   # legacy — the L7 trap
   for n in range(N):
       psi_x = sn_mesh.bc_xmin.apply(psi_x, quad)[n]   # per-ordinate apply
       psi_y = sn_mesh.bc_ymin.apply(psi_y, quad)[n]   # per-ordinate apply
       # ... walk cells, sweep, etc. ...

Each ``bc.apply`` call sees the FULL ``(N, ny, ng)`` face buffer (so
reflective partners can read across rows) and returns an updated full
buffer.  Calling this :math:`N` times per sweep is wrong on two
counts:

1. **Cost** — :math:`N` redundant invocations of the same boundary
   operator on the same buffer.  For a 2-D ``LS-N`` quadrature with
   :math:`N \sim 30`–:math:`80`, this is the dominant per-sweep cost
   on small meshes.
2. **Correctness** — when reflective BCs interact with mid-sweep
   reflective-buffer state, **the order of BC apply vs ordinate
   sweep matters**.  The legacy code's behaviour is sensitive to
   ordinate iteration order; reorderings that algebraically should
   be no-ops (e.g., octant batching, parallel ordinate evaluation)
   silently change the converged solution.

The Wave-2 form applies each boundary operator **once per octant per
axis** — :math:`O(\text{octants}) = 4` calls, not :math:`O(N)`:

.. code-block:: python

   # Wave 2 — L7-trap closed by construction
   for octant in quad.octants:                    # 4 iterations, structural
       sx, sy = octant.label
       ...
       # Apply BC once for this octant on each axis
       if sx_eff >= 0:
           full_face_x = sn_mesh.bc_xmin.apply(psi_x[:, 0, :, :], quad)
           psi_x[oct_idx, 0, :, :] = full_face_x[oct_idx]
       else:
           full_face_x = sn_mesh.bc_xmax.apply(psi_x[:, nx, :, :], quad)
           psi_x[oct_idx, nx, :, :] = full_face_x[oct_idx]
       # ... analogously on y ...
       sweep_graph.walk_windowed(level_op=_CellSolve(...), ...)  # all N_oct batched

The architectural argument: the boundary operator's *semantics* are
"map outgoing partner-octant fluxes to incoming this-octant fluxes".
That mapping is per-octant by construction — applying it once per
ordinate within an octant is redundant; applying it once per octant
is structurally correct.

.. note::

   The ``sn_mesh.bc_xmin.apply(..., quad)`` spellings in the two
   code blocks above are **historical** (the Wave-2 era 2-arg
   ``apply`` on a per-attribute BC surface). Both spellings are
   retired: the 2-arg ``apply`` by Issue #186 (the law is now a pure
   descriptor; the realizer produces a strict 1-arg operator — see
   :ref:`bc-trace-law-descriptor-model`), and the per-attribute
   ``bc_<face>`` surface by C4 / #220 in favour of the
   face-name-keyed :attr:`SNMesh.bc` dict
   (``sn_mesh.bc["xmin"].apply(psi)`` — see
   :ref:`bc-face-name-carve`). The blocks are preserved verbatim
   because they document the *L7-trap structure* the Wave-2 carve
   closed, which is independent of the storage spelling.

The L7-trap detector test
``tests/sn/test_2d_octant_sweep_equivalence.py::case-3`` is the
load-bearing regression gate — a TESTS-FIRST harness (case 3 with
mixed reflective + vacuum BCs, 2G heterogeneous, ``n_sweeps=2``)
designed to fail if any future refactor reintroduces the per-ordinate
BC apply pattern.  The case-3 design uses the post-sweep ``psi_x`` /
``psi_y`` buffer state as bit-identity oracles (rather than the
converged scalar flux), because the L7-trap is invisible in
single-iteration tests: the FIRST iteration's reflective-buffer state
is zero, so per-ordinate vs once-per-octant give the same answer; the
trap surfaces only on the SECOND iteration when the first iteration's
outgoing-face writes feed the second iteration's BC apply.  The case
also explicitly tags ``@pytest.mark.catches("ERR-003")``: ERR-003 is
the catalogued instance where reflective-BC ordering coupled with
ordinate batching produced a converged-but-wrong solution.

Bit-identity to the legacy implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For LS-family quadratures (``LevelSymmetricSN``,
``ProductQuadrature``) whose ordinate ordering is octant-grouped in
lexicographic order, the Wave-2 implementation is **bit-identical**
to the legacy per-ordinate loop on every regression snapshot — the
existing
``tests/sn/regression/snapshots/2d_1g_LS4_dd_15x15.npz``,
``test_apply_2d_cartesian_bit_identical_to_legacy``, and
``test_unified_sweep_dispatch`` snapshots all pass with
``np.array_equal``.  The argument has three parts:

1. **BC apply equivalence.**  The boundary operator for octant
   :math:`\sigma` reads partner-octant rows of the persistent
   ``psi_x`` / ``psi_y`` buffer.  For LS, the lex order of
   ``quad.octants`` matches the legacy n-order at the
   partner-state granularity, so the same iteration's value is
   observed at the same point.  Per-octant BC apply produces the
   same ``psi_x`` / ``psi_y`` octant-row contents as :math:`N_\sigma`
   copies of the legacy per-ordinate apply.
2. **Per-cell sweep equivalence.**  Within an octant, per-ordinate
   cell sweeps are independent — different rows of ``psi_x`` /
   ``psi_y``, different rows of ``angular_flux``.  Batching is
   therefore bit-identical to any per-ordinate sequencing of the
   same set, modulo the per-cell DD operation order which is
   pinned (see :ref:`sweep-octant-dependency-graph` dispatch
   boundary above).
3. **Lebedev (octant ordering not lex).**  For Lebedev quadrature
   (case 6 in the C2.5 harness) the converged scalar flux matches
   the legacy code, but the iter-to-iter values differ on the
   inner iteration (different traversal order ⇒ different
   Gauss-Seidel updates).  Case 6 uses **vacuum BCs**, where the
   partner-state semantics don't matter; this is a deliberate
   choice in the harness — Lebedev with reflective BCs would
   require redesign of the bit-identity gate, but the converged
   answer is still verified at the snapshot level via
   ``test_unified_sweep_dispatch``.

This taxonomy follows the ``vv-principles`` skill's
**bit-identity vs principled-equivalence** discipline: bit-identity
where structurally trivial (LS-family, octant-grouped lex order);
principled equivalence (closed-form L1 anchor + MMS regression suite)
where bit-identity would require more work than the engineering value
returns.

.. _sweep-octant-architecture-cardinal-rule-2:

Architectural framing (Cardinal Rule 2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per the project memory note ``project_moc_structure.md`` and the
:ref:`cell-update-strategies` discussion above,
:class:`SweepDependencyGraph` is **SN-specific by design**.  MoC
will define its own analog (per-ray traversal) — different DAG
shape, different mathematical structure (fiber bundles + solution
sheaves over characteristic curves rather than a topological sort
over a cell graph).  There is **no shared SweepGraph Protocol**
because there is no shared mathematical structure.  Cardinal Rule 2
(architecture) prefers **late unification** ("unify after two
instances" — see ``feedback_unify_after_two_instances`` in agent
memory) to premature abstraction; the sweep DAG lives in
:mod:`orpheus.sn` and stays there until a second mathematically-
similar consumer arrives, which by current understanding is **never**
for MoC and only conjectural for any other deterministic transport
solver.

By contrast, the **angular octant partition** primitive
(:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`) is
genuinely shared infrastructure — see the cross-method consumer
table at :doc:`discrete_measures` (octant partition consumed by SN
2-D, MoC track-bundle direction grouping, MC boundary-current
hemisphere scoring, future SN boundary realiser).  The split is
**measure-level primitives are shared, sweep-level orchestration
is SN-specific**.

Performance
~~~~~~~~~~~

The Wave-2 plan target for Issue #4 closure was 3–10× speedup on
the 421-group benchmark (the canonical ``test_profile_421g``
smoking-gun probe).  The shipped speedups:

.. list-table::
   :header-rows: 1
   :widths: 35 20 15 30

   * - Configuration
     - Speedup
     - Target?
     - Comment
   * - 421-group LS4 ``31×31``
     - 1.7×
     - Below
     - numpy-dispatch-overhead-dominated regime
   * - 2-group LS4 ``31×31``
     - 2.78×
     - At lower bound
     - per-octant batching wins more for small ``ng``
   * - 1-group LS4 ``15×15`` (regression snapshot)
     - bit-identical
     - n/a
     - regression contract preserved

The headline 421-group speedup is below the 3-10× target.  The
honest analysis: the Wave-2 implementation eliminates the
:math:`N`-fold ordinate loop overhead but the per-octant per-level
kernel calls still number :math:`O(\text{levels} \times
\text{octants}) \approx (n_x + n_y - 1) \times 4 \approx 88` per
sweep on a ``31 × 31`` mesh, each carrying its own numpy dispatch
cost.  At 421 groups, the per-call work scales linearly so the
ratio of useful work to dispatch overhead remains modest.  The
**follow-up direction** noted at Wave 2 was to carry full-:math:`N`
buffers plus an ``octant_indices`` field so the kernel calls become
level-only (~ 60 calls / sweep) rather than ``levels × octants``
(~ 240 calls / sweep), eliminating the per-octant copy round-trip.
The subsequent Phase 5 / S6.4 work took a different route to the
same end: the rolling-frontier window
(:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`)
holds the interior cochain on a contiguous :math:`(d{-}1)`-frontier
slab, turning the per-level gather into a basic-slice zero-copy view
(a measured ``~0.77×`` contiguity speedup AND a ``~3×`` peak-memory
win at d=2) — see :ref:`wavefront-flux-cochain`.

Closing the smoking gun by construction is itself a load-bearing
result: the legacy ``for n in range(N)`` is gone, the metric
(angular-flux tensor; see :ref:`theory-sn-index-convention` for the
canonical ``(N, ng, nx, ny)`` storage) now knows its iterative
structure, and any future numpy-dispatch-cost reduction benefits
all closures uniformly through the strategy contract.

Verification
~~~~~~~~~~~~

The Wave-2 verification chain (per the ``algebra-of-record`` skill
discipline):

* **L0 unit tests** — on the primitives:

  - ``tests/sn/test_octants_property.py`` (60 tests across 8
    quadrature factories) — disjoint union, weight conservation,
    sign-signature correctness, pure-axis ordinates labelled
    ``sign=0``.
  - ``tests/sn/test_cell_kernel_batch.py`` (S6.4(e) successor of
    ``test_cell_update_batch.py``) — term-level L0 on the storage-free
    kernel pair (``cell_kernel_batch`` / ``residual_kernel_batch``):
    bit-identity against per-cell :meth:`update` on a
    single-cell-per-batch reduction; standalone tests against
    analytical DD recurrence on a 1×3 strip; 4-octant bit-identity vs
    the per-ordinate Python loop; plus a ``sha256`` source-of-record
    pin on the two kernel bodies (the explicit left-fold order is
    bit-identity-load-bearing).
  - ``tests/sn/test_sweep_graph.py`` — the §15A.2 invariant set above;
    anti-diagonal cell coverage; topo-order acyclicity per octant sign;
    BC face conventions; and the ``walk_full`` / ``walk_windowed`` ×
    level-operation walks (with ``window ≡ full`` bit-identity oracles).
  - ``tests/sn/primitives/test_dag_ownership.py`` (S6.4(c) successor of ``test_snmesh_sweep_graphs.py``) — graph
    contents agree with hand-derived schedule on a 3×3 mesh; dict
    keys equal ``quad.octants`` labels; cache invalidates when mesh
    changes.

* **L1 closed-form anchor + L7-trap detector** — the C2.5 TESTS-
  FIRST harness ``tests/sn/test_2d_octant_sweep_equivalence.py``
  (7/7 pass), tagged ``@pytest.mark.l1`` and
  ``@pytest.mark.catches("ERR-003")``.  Includes:

  - **case 3 (L7 trap)** — mixed BC + 2G heterogeneous +
    ``n_sweeps=2``, the primary L7-trap detector.
  - **case 7 (closed-form)** — 1G homogeneous reflective with
    :math:`k_\infty = \nu\Sigma_f / \Sigma_a`, the structural-
    independence anchor.
  - cases 1–6 covering BC mixes, ordinate batching corners, and
    Lebedev (vacuum-BC variant).

* **L2 regression** — existing ``tests/sn/test_mms_2d.py``,
  ``test_discrete_ordinates_2d.py``, ``test_streaming_operator.py``,
  ``test_streaming_operator_decomposition.py``,
  ``test_unified_sweep_dispatch.py``, ``tests/sn/regression/``: 56/56
  pass, 6 slow-marked skipped.

The verification chain is the canonical
**L0 (primitive invariants) → L1 (closed-form anchor + bug catcher)
→ L2 (integration regression)** ladder from the ``vv-principles``
skill.

References and pointers
~~~~~~~~~~~~~~~~~~~~~~~

* Grand Report v3 §15A.2 (lines 2137–2171) — the "upwind trace
  complex / causal transport DAG / direction sweep ordering"
  primitive description with the ``assert_*`` invariant set this
  module's tests pin.  Plan file at
  ``.claude/plans/neutron_transport_grand_report_v3.md``.
* Wave 2 plan at ``.claude/plans/transient-giggling-cake.md`` (C2.1
  through C2.8) — the architectural primitives plan, sequencing,
  verification-first harness design.
* Wave 0 :meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
  primitive — the measure-level partition that the SN ``octants``
  property delegates to.  See :doc:`discrete_measures`.
* Wave 1 :math:`R \circ \Lambda \circ M` Galerkin scattering
  composition — the parallel "metric knows its iterative structure"
  refactor for the scattering source build.  See
  :ref:`sn-scattering-fission-operators`.
* :class:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase` — the
  strategy ABC carrying the per-cell :meth:`update` reference contract
  and the storage-free batched kernel pair
  :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
  / :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
  (S6.4(e); was the ``SweepCellSlice``-packeted ``update_batch`` /
  ``residual_batch``).
* :class:`~orpheus.sn.spatial.diamond.DiamondDifference` — the only
  shipping closure that overrides the kernel pair; the reference for
  the bit-identity contract (pure cell algebra — the ONLY
  direction-aware math in the SN spatial stack since S6.4(e) lifted
  storage to the walk layer).
* :mod:`orpheus.sn.sweep_graph` — the two storage walks
  (:meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_full`,
  :meth:`~orpheus.sn.sweep_graph.SweepDependencyGraph.walk_windowed`)
  and the ``_CellSolve`` / ``_CellResidual`` level operations.
* C2.5 TESTS-FIRST harness:
  ``tests/sn/test_2d_octant_sweep_equivalence.py``.

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
:func:`~orpheus.sn.loss_representation.transport_sweep` entry point that branches
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
:meth:`~orpheus.sn.spatial.scheme.DiscretizationScheme.update`:

.. code-block:: python

   for cell_idx in march_order:
       st = reduced.streaming_terms(cell_idx, dir_idx, mu_level_idx=p)
       upstream = UpstreamState(
           spatial_upstream=psi_face,
           angular_upstream=psi_angle[cell_idx],
       )
       result = sn_mesh.scheme.update(
           streaming_terms=st,
           total_xs=sig_t[cell_idx],
           source=QV[cell_idx],
           upstream_state=upstream,
       )
       psi_face = result.outgoing_spatial_flux  # may be None for cylindrical degenerate
       psi_angle[cell_idx] = result.outgoing_angular_state

The cell-update strategy lives on
:attr:`~orpheus.sn.geometry.SNMesh.scheme` (introduced in
this round as a constructor argument with default
:class:`~orpheus.sn.spatial.diamond.DiamondDifference`).  The
default reproduces the inlined sweep math bit-identically — every
regression snapshot at ``tests/sn/regression/snapshots/`` was
generated with DD and continues to match bit-for-bit when the
unified sweep dispatches via ``scheme.update(...)``.  See
:ref:`cell-update-strategies` for the strategy contract and
:class:`~orpheus.sn.spatial.diamond.DiamondDifference` for the
DD scalar form.

Wave C-extension will ship :class:`Step`, :class:`LinearDiscontinuous`,
and :class:`ExponentialCharacteristic` as positivity-preserving /
higher-order alternatives; the unified dispatch infrastructure is
in place to receive them — users will pass
``scheme=LinearDiscontinuous()`` etc. at
:class:`~orpheus.sn.geometry.SNMesh` construction.

The 1-D cumprod fast path (DD-only)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Cartesian dispatch checks three preconditions before
selecting the 1-D cumprod fast path:

#. ``scheme is DiamondDifference`` — the cumprod recurrence
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

The 2-D wavefront sweep dispatches its DD per-cell algebra
through the storage-free kernel pair
:meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
/ :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
on the strategy attached to the
:class:`~orpheus.sn.geometry.SNMesh` (Wave 2 of the SN
performance plan; closes Issue #4 — see
:ref:`sweep-octant-dependency-graph` for the full architecture and
:ref:`sweep-dispatch-relayering` for the S6.4(e) re-layering).
The "inlined DD math" formerly carried inside
:func:`~orpheus.sn.loss_representation._sweep_jacobi` was lifted into
:class:`~orpheus.sn.spatial.diamond.DiamondDifference` as a single
bit-identical extraction, vectorised across the
``(N_oct, n_diag, ng)`` slice — the ordinate axis, anti-diagonal
axis, and group axis simultaneously.  Wave C-extension's LD / EC
/ Step closures override the kernel pair and become drop-in
alternatives at SNMesh construction time:
``SNMesh(mesh, quad, scheme=LinearDiscontinuous())``.  The
open design point of "how to parameterise the 2-D wavefront
without breaking anti-diagonal vectorisation" is now closed: the
storage walk is the scheduler, the level operation owns the
direction fork, and the kernel pair is the closure — the contract
is per-level batched evaluation.

ERR-026 closure status (partial through Wave E)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Superseded (2026-06-12, Issue #195).**

   This subsection records the Wave-E-era reading: ERR-026 PARTIAL,
   the curvilinear ``"krylov"`` default "would regress MMS to
   :math:`\mathcal{O}(h)`", the open second-order follow-up.  That
   reading was the best available *then* and is preserved as history.
   It is now **superseded**: ERR-058 (#195) showed the wrong fixed
   point was the *closure-seed* family, not a boundary-truncation
   order; the curvilinear default returned to ``"source_iteration"``
   (SI :math:`\equiv` Krylov bit-identical post-unification); and the
   isotropic MMS is :math:`\mathcal{O}(h^2)`-consistent.  See
   :ref:`sn-err-058-closure-seed-closeout` for the mechanism,
   structural obstruction, and production decision.  The numerical
   values below stay as bug-era evidence; their *interpretation* is
   carried by the close-out.

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

Wave D's gating contract was bit-identity for ``scheme =
DiamondDifference`` — the bug is preserved by construction so
the regression snapshots stay green.  Wave E (Issues #98 #99 #164)
took two passes at the closure:

* **Wave E Round 2** wired
  :func:`~orpheus.sn.solver.solve_sn_fixed_source` to route
  through Krylov-on-:meth:`InvertibleOperator.apply` (the
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
  :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances on
  the :class:`~orpheus.sn.geometry.SNMesh` (Wave B Issue 7
  tensor-decomposed BC algebra), dispatching boundary fills via the
  realiser-routed 1-arg :meth:`apply` on the resolved
  :class:`~orpheus.numerics.operator.LinearOperator`. Vacuum,
  reflective, white, periodic, albedo, and mixed BCs are now
  plumbed uniformly through the FD operator; bit-identity to the
  pre-Round 3 hard-coded reflective fill is preserved for
  :class:`ReflectiveBoundary(axis=…, albedo=1.0)` (the standard
  ``BC.reflective`` case), which is the load-bearing condition for
  the 11 frozen regression snapshots to stay green. (Note:
  post Issue #186 / B3 + β2 the law itself is a pure descriptor;
  ``BoundaryTraceLaw.apply`` no longer exists. The realiser
  produces the 1-arg :class:`LinearOperator` whose :meth:`apply`
  the matvec calls; the Wave-E Round-3 prose describes the
  contract as it existed at the time, but the architectural
  conclusion — uniform BC consumption through a 1-arg
  ``apply`` — is the same.)

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

Boundary face-flux strategies — Phase A (RETIRED in Phase C)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Issue #168 empirical investigation (recorded at
``.claude/agent-memory/numerics-investigator/issue_168_three_defects.md``)
found **three independent O(1) boundary truncation defects** in the
historical curvilinear FD operator.  Phase A (2026-05-10) addressed
**Defects 1 + 2** via a ``BoundaryFaceFlux`` strategy Protocol — a
one-sided second-order DD diamond extrapolation
:math:`\psi^{\text{face}}_{N-1/2} = \tfrac{3}{2}\,\psi_{N-1} -
\tfrac{1}{2}\,\psi_{N-2}` plus a structural decoupling of
cell-centre storage from BC face-value storage in
:func:`~orpheus.sn.operator.solution_to_angular_flux_spherical`
(returning a ``(fi, boundary_face_flux)`` tuple where ``fi`` was
pure cell-centre storage and the BC face flux lived in its own
companion array).

**Phase C retired the Phase A Protocol entirely.** The sweep-frame
apply matvec subsumes the boundary-face closure into the WDD
propagation chain — the BC trace law owns the boundary edge per
the §16A.3 contract, no separate algebraic extrapolation is needed.
The retired symbols are:

* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`` (Protocol)
* ``orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`` (default)
* ``orpheus.sn.spatial.boundary_face_flux.CellCenter`` (ablation)
* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFluxBase`` (ABC)
* The ``boundary_face_flux`` field on
  :class:`~orpheus.sn.geometry.SNMesh`
* The 21 foundation tests at
  :file:`tests/sn/spatial/test_boundary_face_flux.py`

See :ref:`phase-c-sweep-frame-matvec` for the replacement
architecture. The Phase A subsection is preserved as historical
context for the empirical-defects-investigation reasoning chain.

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
``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``), :mod:`orpheus.sn.spatial.diamond`, and the
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

.. _sn-pole-closure-compute-psi-half:

Half-angle grid exposure (Issue #197 PR-TYPED-6b)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo:: Archivist expansion needed.

   The Issue #197 PR-TYPED-6b dispatch added the public method
   :meth:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep.compute_psi_half_per_level`
   exposing the M-M recurrence's half-angle grid
   :math:`\phi_{m\pm 1/2,i,g}` for one level.  The method is the
   intermediate exposure that lets the unified SN matvec
   (:class:`~orpheus.sn.operator.StreamingOperator`) consume
   :math:`\phi_{m\pm 1/2}` as
   :func:`~orpheus.sn.spatial.cell_balance.cell_balance_for_streaming`'s
   ``psi_angular_upstream`` argument — closing the apply-vs-sweep
   twin path on the curvilinear angular branch.

   Pattern 2 (Single source of truth).  Both the public method AND
   the redistribution body inside
   :meth:`MorelMontryAngularSweep.__call__` route through
   :func:`~orpheus.sn.spatial.pole_angular_closure._mm_psi_half_grid_single_level`.
   The :eq:`pole-mm-recurrence` recurrence body lives once.

   Test gate:
   :file:`tests/sn/spatial/test_compute_psi_half_per_level.py` —
   21 foundation + L0 tests pinning method existence,
   shape contract, recurrence formula, Carlson-context seed
   contract, Pattern 2 round-trip to ``__call__``, linearity, and
   refactor-regression on
   :func:`_mm_weighted_angular_recurrence_single_level`.

   Closeout memo:
   ``.claude/agent-memory/method-implementer/issue_197_pr_typed_6b_closeout.md``.

.. _phase-c-sweep-frame-matvec:

Sweep-frame apply matvec (Issue #168 Phase C)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: Key Facts
   :class: important

   * Phase C (commits ``eae6f05`` → ``d445a8f``, 2026-05-12) rewrites
     :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
     and ``_cylindrical`` as **one sweep iteration semantically**.
     The WDD diamond closure
     :math:`\psi^{\text{face}}_{\text{out}} = 2\,\psi^{\text{cell}}
     - \psi^{\text{face}}_{\text{in}}` propagates the face flux
     cell-by-cell along the direction's DAG; the BC trace law owns
     the boundary edge per :ref:`affine-bc-form`.
   * The pole-face initial condition is
     :math:`\psi^{\text{face}}_{\text{in}}(\text{pole}) =
     \psi^{\text{cell}}[0]` (Lewis–Miller §4.5 Carlson seed), **not**
     :math:`0`. The Carlson seed is the unique anchor that preserves
     the per-ordinate flat-flux invariant under the WDD recurrence.
   * Phase A's
     :class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`
     Protocol (415 LOC + 21 foundation tests) **retires entirely**
     — the boundary-face closure is now inside the WDD propagation
     chain, owned by the BC trace law at the boundary edge.
   * Empirical Gate 1.1 (per-ordinate flat-flux residual on
     reflective curvilinear) finds **spherical** MMS with
     ``MorelMontryAngularSweep`` FAILS, **cylindrical** MMS PASSES.
     The default flip to ``MorelMontryAngularSweep`` is therefore
     DEFERRED to Phase D (`Issue #192
     <https://github.com/deOliveira-R/ORPHEUS/issues/192>`_); the
     four ERR-026 ``xfail-strict`` curvilinear MMS tripwires stay
     xfail; ERR-026 remains at PARTIAL CLOSURE.
   * The structural-frame name for the rewrite is the **sweep /
     wavefront frame** (cross-domain-attacker 2026-05-12 analysis):
     the "ghost cell" idiom is realised as a typed
     :class:`~orpheus.numerics.spaces.trace_space.TraceSpace` vector
     defined by the realised BC operator, not extrapolated from
     interior cell centres.

Phase C resumes the spatial-closure alignment that Phase B
identified as the load-bearing precondition for the
curvilinear-default flip (see
:ref:`sn-pole-angular-closure-protocol` for Phase B's full
:class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
Protocol narrative; this subsection picks up at the point Phase B
left for follow-up). Three unblockers landed before Phase C
resumed:

1. The trajectory_resolvent (Peierls Variant α Green's-function)
   campaign shipped cylinder MR Phase 1b (commits ``37e3e29``,
   ``cf662a6``, ``604f380``, ``e10c33c``), explicitly built to close
   the cylinder-2G ERR-026 gap. trajectory_resolvent now covers
   5 of the 6 deleted curvilinear regression snapshots at
   machine-precision-class precision; the 6th (P1 anisotropic)
   routes to a shape-independent :math:`k_\infty` closed form.
2. The cross-domain-attacker analysis on 2026-05-12 identified a
   **structural inconsistency** between the pre-Phase-C apply
   matvec's boundary closure and the :ref:`affine-bc-form` (§16A.3).
   The apply matvec was passing **cell-centre** values to
   :meth:`bc_outer.apply`, whereas §16A.3 requires the boundary
   face TRACE :math:`\gamma_+ \psi`. The cross-domain-attacker's
   "sweep / wavefront frame" naming gave the rewrite its
   architectural shape: the apply matvec is one sweep iteration
   over the cell-visit DAG, with the BC trace law at the boundary
   edge owning the inflow trace.
3. The Phase A
   :class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`
   Protocol — built to second-order-accuratise the curvilinear
   outer face — was re-classified as **a patch on top of the wrong
   architecture**. Phase A's
   :class:`~orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`
   produces a face-flux extrapolant
   :math:`\psi^{\text{face}}_{N-1/2} = \tfrac{3}{2}\,\psi_{N-1} -
   \tfrac{1}{2}\,\psi_{N-2}` that ignores the BC entirely; Phase B's
   :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
   produces angular closure that is inconsistent with the apply
   matvec's spatial closure (arithmetic interior average + DD
   extrapolation outer face). The fix at all three sites is the
   **same** — make every face closure consume the WDD-propagated
   face value, and let the BC trace law own the boundary edge.
   "Two paths to the same discrete operator over different storage
   conventions" (cross-domain-attacker Smell 16) is the trigger for
   the unification.

The pre-Phase-C arithmetic spatial closure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-Phase-C, the matvec interior face values were computed as the
arithmetic average of cell-centre values
:math:`\psi^{\text{face}}_{i+1/2} = \tfrac{1}{2}(\psi^{\text{cell}}_i
+ \psi^{\text{cell}}_{i+1})`. This is **second-order accurate** on
smooth fields when the cell-centre values are themselves the
analytical values, but it does **NOT** match the sweep's WDD
recurrence which evaluates
:math:`\psi^{\text{face}}_{i+1/2} = 2\,\psi^{\text{cell}}_i -
\psi^{\text{face}}_{i-1/2}`. The two are equivalent only when
:math:`\psi` is constant on a cell; for any angular or spatial
variation the two values diverge by an :math:`\mathcal{O}(\Delta r)`
amount, and the operators are **not the same operator**.

This is the empirical content of Phase B's diagnosis (see
:ref:`sn-pole-angular-closure-protocol` "ERR-026 closure status"):
pairing the canonical
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
angular closure with arithmetic-average spatial closure produces a
**worse** operator on flat :math:`\psi` than the legacy
:math:`\tau`-symmetric interpolation. The canonical M–M form
produces half-angle face fluxes oscillating as
:math:`0, 2c, 0, 2c, \ldots` on flat :math:`\psi` (this is the
correct DD-angular behaviour when seeded with the Carlson starting
direction :math:`\psi_{1/2} = 0`), and the arithmetic spatial
average then combines these oscillating angular face fluxes with
interior-averaged spatial face fluxes into garbage. Phase B's
empirical test ``test_apply_spherical_constant_flux_under_morel_montry_canonical_form``
saw :math:`\phi` range across [0.6, 1.004] on a flat-:math:`\psi`
input that analytical balance demands give exactly :math:`1.0`.

The fix is to align the spatial closure with the sweep: use the
WDD recurrence
:math:`\psi^{\text{face}}_{\text{out}} = 2\,\psi^{\text{cell}}
- \psi^{\text{face}}_{\text{in}}` per cell, propagating face flux
in DAG order. Phase C ships this alignment.

The sweep-frame matvec algebra
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The rewritten matvec is **one sweep iteration semantically**.
For each bulk direction :math:`d \in \{+1, -1\}`, the per-cell
WDD recurrence walks the face flux along the DAG:

.. math::
   :label: phase-c-wdd-recurrence

   \psi^{\text{face}}_{\text{out}}(i) \;=\; 2\,\psi^{\text{cell}}(i)
   \;-\; \psi^{\text{face}}_{\text{in}}(i),
   \qquad
   \psi^{\text{face}}_{\text{in}}(i+1) \;=\;
   \psi^{\text{face}}_{\text{out}}(i),

evaluated cell-by-cell across the direction's DAG order yielded by
:meth:`~orpheus.sn.geometry.SNMesh.dag_walk` invoked with
``direction_sign``. The
per-cell streaming term consumes both the inflow and outflow face
values along with the cell volume and face areas (Hébert §3.9.4
balance, ORPHEUS-normalised):

.. math::
   :label: phase-c-streaming-spherical

   S_{n,i,g} \;=\; \frac{\mu_n}{V_i}\,
   \bigl[ A_{i+1/2}\,\psi^{\text{face}}_{n,i+1/2,g}
        - A_{i-1/2}\,\psi^{\text{face}}_{n,i-1/2,g} \bigr],

with the redistribution term provided by the Phase B
:class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
strategy and the collision term :math:`\Sigma_t \psi^{\text{cell}}_{n,i,g}`
unchanged. The full per-cell update is

.. math::
   :label: phase-c-cell-update

   (T\psi)_{n,i,g} \;=\; S_{n,i,g} \;+\; R_{n,i,g} \;+\;
   \Sigma_t(i,g)\,\psi^{\text{cell}}_{n,i,g},

where :math:`R_{n,i,g}` is the
:meth:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure.__call__`
output (Bailey :math:`\Delta A_i / w_n` redistribution factor with
the strategy-specific :math:`\alpha_{n+1/2}\psi_{n+1/2} -
\alpha_{n-1/2}\psi_{n-1/2}` evaluation). See
:eq:`bailey-dome-recursion` for the :math:`\alpha` recurrence and
:eq:`pole-mm-recurrence` for the M–M angular DD form.

Ordinate vectorisation
^^^^^^^^^^^^^^^^^^^^^^^

Per the user's hard architectural directive (and the precedent at
``orpheus/sn/angular_operator.py:183``), the rewritten matvec
carries **no** ``for n in range(quad.N)`` loop. The per-cell update
operates on whole ordinate subsets via boolean masks:

.. code-block:: python

   eps = 1e-15
   outgoing_mask = quad.mu_x > +eps   # μ > 0 ordinates
   incoming_mask = quad.mu_x < -eps   # μ < 0 ordinates
   mu_out = quad.mu_x[outgoing_mask]
   mu_in  = quad.mu_x[incoming_mask]

A single per-cell statement updates the full outward ordinate
subset:

.. code-block:: python

   psi_cell = fi[:, outgoing_mask, i, 0]               # (ng, n_out)
   psi_face_out = 2.0 * psi_cell - psi_face_in         # WDD diamond
   streaming = (
       mu_out[None, :]
       * (A[i + 1] * psi_face_out - A[i] * psi_face_in)
       / V[i]
   )

This is the canonical "vectorise across ordinates" pattern that
the cross-method ordinate-anti-pattern audit
(`Issue #191 <https://github.com/deOliveira-R/ORPHEUS/issues/191>`_)
tracks systemically. Phase C's contribution is to introduce no new
per-ordinate loop inside any new code; the 14 existing sites stay
untouched as separate work.

The cylindrical case adds an outer loop over the :math:`\mu`-levels
(the per-level azimuthal-DAG topology is intrinsic to cylindrical's
structure), but each level still operates on whole within-level
ordinate subsets via the same masking pattern. A third pass handles
pure-azimuthal degenerate ordinates (:math:`|\eta_n| < 10^{-15}`)
whose streaming term is zero by construction but whose
redistribution + collision contributions still must be scattered to
the equation map.

The new APIs
^^^^^^^^^^^^

Two new APIs surface what the existing infrastructure already knew:

* :meth:`~orpheus.sn.geometry.SNMesh.dag_walk(*, ordinate_idx=...,
  direction_sign=..., mu_level_idx=None)` — Issue #196 Phase G
  Step 2.6 (Q3) canonicalised this as the **single iteration
  primitive** for 1-D sweeps, replacing the legacy pair of
  ordinate-keyed / direction-keyed methods.  Exactly one of
  ``ordinate_idx`` or ``direction_sign`` is supplied (XOR):
  ``ordinate_idx=n`` for the sweep driver's per-ordinate march,
  ``direction_sign=±1`` for the apply matvec's whole-subset walk
  keyed by the **bulk sweep direction sign**. The existing
  cell-visit graph's per-quadrant ``_diag_cache`` is already keyed
  by :math:`(\mathrm{sign}(\mu_x), \mathrm{sign}(\mu_y))`; the
  direction-keyed branch surfaces that sign-only view as a
  first-class API. A foundation test
  (:file:`tests/sn/test_dag_walk.py`) pins
  bit-identity between the two invocation modes across sphere /
  slab / cylindrical for every representative ordinate. For
  cylindrical the per-level
  ``mu_level_idx`` is required (the within-level DAG topology
  differs per level).

* :meth:`~orpheus.sn.operator.EquationMap.unknowns_at_cell_for_mask(cell_idx, ordinate_mask)`
  — a precomputed inverse lookup ``(cell, ordinate) → k``. Lazily
  builds an ``(nx, N) int`` table with :math:`-1` sentinels for
  absent ``(ordinate, cell)`` slots; subsequent calls are O(1) per
  ``(cell, mask)`` pair. Replaces the per-equation O(n_eq) linear
  scan the legacy scatter pattern used. The eq_map still iterates
  ``(spatial outer, ordinate inner)`` at construction time; the
  helper just precomputes the inverse for the sweep-frame matvec.

What retires
^^^^^^^^^^^^

Phase A's
:class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`
Protocol — five symbols, the
:class:`~orpheus.sn.geometry.SNMesh` field, and the 21 foundation
tests — retires entirely. The architectural reasoning is "two paths
to the same operator → unify after the second instance" (per the
:doc:`/development` agent memory ``Unify after two instances``
directive). Phase A was the first instance of a face-closure
strategy; Phase B was the second (angular closure). With Phase C
the **third** instance (the BC at the boundary edge) the unification
has cleaner shape: every face value comes from the WDD recurrence
or the BC trace law applied at the boundary edge, never from an
algebraic extrapolation of cell centres. The retired symbols are:

* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`` (Protocol)
* ``orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFluxBase`` (ABC)
* ``orpheus.sn.spatial.boundary_face_flux.DDExtrapolation`` (default strategy)
* ``orpheus.sn.spatial.boundary_face_flux.CellCenter`` (ablation strategy)
* The ``boundary_face_flux`` constructor field +
  :class:`~orpheus.sn.geometry.SNMesh` attribute
* The ``boundary_face_flux_closure`` keyword argument from
  :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
  and ``_cylindrical``
* :file:`tests/sn/spatial/test_boundary_face_flux.py` (232 LOC,
  21 foundation tests)

Three additional simplifications ship with the rewrite:

* :func:`~orpheus.sn.operator.solution_to_angular_flux_spherical`
  (and its cylindrical alias) returns a single ``fi`` array
  ``(ng, N, nx, 1)`` instead of the Phase A
  ``(fi, boundary_face_flux)`` tuple. Inward-at-boundary cell-centre
  slots ``fi[:, n_inward, -1, 0]`` are filled with the
  **reflected-partner cell-centre value** as an analytical
  extension: the equation map excludes these from unknowns (the BC
  determines them), but the WDD recurrence on flat :math:`\psi`
  requires the cell-centre to be consistent so the per-ordinate
  flat-flux invariant holds.
* :class:`~orpheus.sn.geometry.SNMesh` no longer accepts the
  ``boundary_face_flux=`` keyword (a regression test pins the
  field retirement).
* :class:`~orpheus.sn.operator.InvertibleOperator.apply` dispatch
  drops the ``boundary_face_flux_closure`` plumbing.

What stays
^^^^^^^^^^

* The Phase B
  :class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
  Protocol stays. The sphere centre / cylinder axis is **intrinsic
  geometry** (a coordinate-system singularity, not an external BC),
  so the three-strategy Phase B Protocol is the right shape; only
  the **default** is under question, and that is the Phase D
  decision point.
* The :meth:`~orpheus.sn.operator.InvertibleOperator.apply_transpose`
  machinery via dense-probe construction stays. Linearity of the
  rewritten :meth:`~orpheus.sn.operator.InvertibleOperator.apply`
  (Gate 1.4, pinned to ``rtol=1e-13``) guarantees the transpose is
  correctly tracked.
* The
  :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer` +
  :class:`~orpheus.sn.boundary_realizer.SNMethodSpace` +
  :class:`~orpheus.numerics.operator.LinearOperator`-1-arg
  ``apply`` substrate (Issues #186 + #176 + #188, Waves 0–12) — the
  BC trace law's realised 1-arg ``apply(outflow) → inflow``
  contract is exactly what the matvec consumes at the boundary
  edge.

The pole-face Carlson starting direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The single largest architectural deviation between the Phase C plan
and the shipped code is the **pole-face initial condition** for the
outward WDD sweep. The plan's pseudocode wrote

.. code-block:: python

   psi_face_in = np.zeros((ng, n_out))   # plan's pseudocode

with the comment "Pole face: :math:`\psi^{\text{face}} = 0` by
symmetry (also multiplied by :math:`A_0 = 0`)". The first claim
(symmetry) is wrong; the second (annihilation by zero face area) is
correct but does not help propagation. Empirically, this initial
condition combined with the WDD recurrence on flat :math:`\psi`
produces oscillating face fluxes
:math:`0, 2c, 0, 2c, \ldots` that break the per-ordinate flat-flux
invariant on **all three** Phase B pole-closure strategies:

.. math::
   :label: phase-c-wdd-oscillation

   \psi^{\text{face}}_0 = 0, \quad
   \psi^{\text{cell}}_0 = c \;\Longrightarrow\;
   \psi^{\text{face}}_1 = 2c - 0 = 2c,

   \psi^{\text{cell}}_1 = c \;\Longrightarrow\;
   \psi^{\text{face}}_2 = 2c - 2c = 0,

   \psi^{\text{cell}}_2 = c \;\Longrightarrow\;
   \psi^{\text{face}}_3 = 2c - 0 = 2c, \ldots

The correct initial condition is the **Carlson starting-direction
seed** of [LewisMiller1984]_ §4.5 (paraphrased in Hébert §3.9.4
Eqs. 3.432–3.435 for the angular analogue):
:math:`\psi^{\text{face}}_{\text{in}}(\text{pole}) =
\psi^{\text{cell}}(\text{first cell})`. For true flat :math:`\psi`
this yields

.. math::

   \psi^{\text{face}}_0 = c, \quad
   \psi^{\text{face}}_1 = 2c - c = c, \quad
   \psi^{\text{face}}_2 = 2c - c = c, \quad \ldots

at every cell, so the streaming term
:math:`\mu_n (A_{i+1}\psi^{\text{face}}_{i+1} -
A_i\psi^{\text{face}}_i)/V_i` cancels the redistribution term per
ordinate on flat :math:`\psi` and the per-ordinate flat-flux
invariant holds. The pole-face streaming contribution is still
multiplied by :math:`A_0 = 0` (the pole face has zero area in both
spherical and cylindrical 1-D), so the Carlson seed introduces no
spurious source there; it only **anchors the recurrence** for the
cell-by-cell WDD propagation across the interior.

Why Lewis–Miller §4.5 is the canonical reference: at the spherical
centre :math:`r=0` and the cylindrical axis the angular dependence
of :math:`\psi` becomes **structurally singular** in the
transport-theory sense ([Pomraning1989]_ p. 339; see also
:ref:`sn-phase-d-pomraning-structural-singularity` and earlier
treatments in [LewisMiller1984]_ §4.5) — the angular flux is not a
separable function of
:math:`(\mu, r)` in any neighbourhood of the singular point because
the inward and outward ordinate cones meet there. Lewis–Miller's
"starting direction" handles this by introducing a half-step inward
sweep at :math:`\mu = -1` that initialises the
:math:`\alpha`-cascade and propagates to the outward sweep; the
Carlson seed is the natural anchor for that half-step in the
cell-by-cell WDD formulation. The same logic applies to the
cylindrical axis: the per-level azimuthal-DAG topology has a
half-step inward-zero-weight ordinate at :math:`\mu_x = -1` that
anchors each level's :math:`\alpha`-cascade.

The cylindrical analogue uses the identical Carlson seed per level:
:math:`\psi^{\text{face}}_{\text{in,level}}(\text{pole}) =
\psi^{\text{cell}}(\text{first cell at level})`. The cylindrical
per-level :math:`\alpha`-dome telescoping then absorbs the
half-angle face flux discrepancy in a way the spherical case does
not — see the Gate 1.1 finding below.

.. _phase-c-apply-sweep-equivalence:

apply ↔ sweep structural equivalence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-Phase-C, the matvec's :meth:`apply` and the sweep's :meth:`solve`
were structurally distinct paths to the **same** discrete operator
:math:`L`: :meth:`apply` walked cell-centre storage with arithmetic
face averages, :meth:`solve` walked face storage with WDD
asymmetric propagation. The cross-domain-attacker frame analysis
(Smell 16, ``.claude/agent-memory/cross-domain-attacker/issue_168_phase_c_sweep_frame.md``)
flagged this as the elegance-smell trigger that Phase C resolves:

   *Two paths to the same discrete operator over different storage
   conventions, with order-degradation at boundaries on one path.
   FRAME: Sweep / wavefront — recover the boundary as a DAG edge
   consumed via the trace law, not as a cell-centre algebraic
   closure.*

Under the sweep-frame architecture both paths consume the **same
three primitives**:

#. The WDD diamond closure :eq:`phase-c-wdd-recurrence` per cell.
#. The direction-keyed cell-visit DAG via
   :meth:`~orpheus.sn.geometry.SNMesh.dag_walk` invoked with
   ``direction_sign=±1``.
#. The BC trace law applied **once** at the boundary edge per
   :ref:`affine-bc-form`.

The face-flux propagation identity therefore holds **by
construction** post-Phase-C: extracting the implicit face fluxes
from ``apply`` (by inverting the cell-balance equation) recovers
the same WDD recurrence the sweep walks. The structural-frame
identity is the load-bearing acceptance criterion for
preconditioned-Krylov stability — when ``apply`` is the operator
:math:`L` and the sweep is :math:`L^{-1}` (approximately), they
must agree on what :math:`L` **is**. Foundation tests in
:file:`tests/sn/test_phase_c_gates.py` pin this via
``np.array_equal`` on:

* **Gate 1.2** (apply determinism) — repeat-call invariance
  ``apply(ψ) == apply(ψ)`` bit-identical across two invocations.
* **Gate 1.3** (apply ↔ apply_transpose reciprocity)
  :math:`\langle L\psi, \phi \rangle = \langle \psi, L^T\phi \rangle`
  to ``rtol=1e-12, atol=1e-13``. Free if Gate 1.4 (linearity)
  passes.
* **Gate 1.4** (apply linearity) :math:`L(\alpha\psi + \beta\phi)
  = \alpha L\psi + \beta L\phi` to ``rtol=1e-13``.
  **Precondition** for Gates 1.2 + 1.3 + the dense-probe
  ``apply_transpose`` construction.

.. _bc-trace-contract-respected-by-matvec:

BC trace contract respected by matvec
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note::

   **Superseded by Wave O steps O.4a.2 + O.4b (Issue #208).** This
   subsection documents the **Phase C** matvec, where the boundary law
   was applied *inside* the sweep at the boundary edge
   (``inflow_full = bc_outer.apply(outflow_at_boundary.T)``, the
   "keystone"). That intra-sweep BC re-apply has been **deleted** for
   **every geometry**: O.4a.2 made the **1-D** path bare, and O.4b
   made the **2-D Cartesian wavefront** path bare. The boundary law
   :math:`B` is now a first-class sibling :math:`-B` operator and the
   sweeps read the seeded inflow trace directly. The Phase C
   trace-contract *insight* — the BC must consume the WDD-propagated
   outflow face vector, not cell centres — is **preserved and
   strengthened** by the extraction: the outflow trace is now an
   explicit solved unknown :math:`\psi.\text{outflow}` rather than a
   local variable, and :math:`B` reads it as
   :math:`B\,\psi.\text{outflow}`. See :ref:`bare-sweep-extraction`
   below and the canonical algebra at :ref:`bc-extraction` in
   :doc:`operator_algebra`.

The :doc:`boundary_conditions` :ref:`affine-bc-form` (§16A.3) reads

.. math::

   \gamma_- \psi \;=\; R\,G\,\gamma_+\psi \;+\; q,

requiring the BC operator to consume the **boundary face trace**
:math:`\gamma_+\psi`, not an interior cell-centre approximation
of the trace. The pre-Phase-C ``operator.py:533`` site read

.. code-block:: python

   outgoing = fi[:, :, -1, 0].T              # ← cell-centres
   incoming = bc_outer.apply(outgoing)

silently violating the contract. The contamination was invisible
for **specular reflection at** :math:`\alpha=1` because the
reflection permutation commutes with cell-centre fills (the same
permutation pattern applies to whatever value sits at the boundary
slot), but it surfaces for every other BC: vacuum (:math:`\alpha=0`),
albedo (:math:`0 < \alpha < 1`), prescribed inflow (any nonzero
:math:`q`), and white BC. Each of these is a regime where
higher-order spatial accuracy should appear — and where the
cell-centre approximation degrades the operator to first-order
boundary truncation.

The Phase C matvec honours the contract by construction: the BC
operator's input is the **WDD-propagated outflow face vector**, not
cell centres. The boundary-edge sequence is:

.. code-block:: python

   # ── Phase 1: outward sweep (μ > 0), i = 0 → nx-1 ─────────────
   # Carlson seed at pole: ψ_face_in = ψ_cell[0].
   for visit in sn_mesh.dag_walk(direction_sign=+1):
       i = visit.cell_idx
       psi_cell = fi[:, outgoing_mask, i, 0]
       psi_face_out = 2.0 * psi_cell - psi_face_in          # WDD
       # ... streaming + redistribution + collision scatter ...
       psi_face_in = psi_face_out                            # walk
   # The last cell's ψ_face_out is the boundary outflow face.
   outflow_at_boundary[:, outgoing_mask] = psi_face_out

   # ── BC trace law at boundary edge ────────────────────────────
   # bc_outer is the realised BC (BoundaryTraceLaw) from SNMethodSpace
   # via SNBoundaryRealizer.realize() — a 1-arg LinearOperator
   # whose apply maps Γ_+ → Γ_-, per the affine-bc-form contract.
   inflow_full = bc_outer.apply(outflow_at_boundary.T)

   # ── Phase 2: inward sweep (μ < 0), i = nx-1 → 0 ──────────────
   psi_face_in = inflow_full[incoming_mask, :].T            # BC-set
   for visit in sn_mesh.dag_walk(direction_sign=-1):
       i = visit.cell_idx
       psi_cell = fi[:, incoming_mask, i, 0]
       psi_face_out = 2.0 * psi_cell - psi_face_in          # WDD
       # ... walk ...

The :class:`~orpheus.numerics.operator.LinearOperator` returned by
:meth:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer.realize`
internally consumes its
:class:`~orpheus.numerics.spaces.trace_space.OutflowTraceSpace` mask
(the ordinate slots with :math:`\mu \cdot \hat n > 0` at the face)
and writes only the inflow slots; the outflow slots in the output
are unspecified by the §16A.3 contract and the matvec reads back
only ``inflow_full[incoming_mask, :]``. This is the user's "ghost
cell for higher-order boundary closure" idiom realised as a typed
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` vector
defined by the realised BC operator — not extrapolated from
interior cell centres.

**Gate 1.5** (foundation,
:func:`tests.sn.test_phase_c_gates.test_bc_trace_contract_respected_by_matvec`)
pins this contract: for each
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` concrete kind
(``VacuumInflow`` / ``ReflectiveBoundary`` / ``WhiteBoundary`` /
``AlbedoBoundary`` / ``PrescribedInflow``), the apply matvec's BC
integration consumes the WDD-propagated outflow face value (not
cell centres) and produces the inflow face value consistent with
``bc.realize().apply(outflow_at_boundary)``. The assertion is
bit-identical across 5 random :math:`\psi^{\text{cell}}` inputs per
BC kind × geometry × ordinate count. ``apply(0) = 0`` is pinned
separately under vacuum + reflective BCs; under vacuum BC a flat
cell-centre :math:`\psi` produces a **non-zero** residual (the BC
physically removes the inflow), and that asymmetry is itself a
load-bearing acceptance criterion: a vacuum BC that left flat-flux
residual at zero would be the pre-Phase-C cell-centre contamination
returning silently.


.. _bare-sweep-extraction:

The bare sweep (Wave O step O.4a.2)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Wave O step O.4a.2 (Issue #208, commits ``d7e1316`` / ``4c0ff96`` /
``2bdc66d``, 2026-06-03) **removed the boundary law from the 1-D
sweep entirely**. The boundary-edge ``inflow_full =
bc_outer.apply(outflow_at_boundary.T)`` line shown above (the
"keystone") is **deleted**; the 1-D :func:`~orpheus.sn.loss_representation.transport_sweep`
entry now reads the *seeded* inflow trace directly:

.. code-block:: python

   # PRE-O.4a.2 (bc-in-sweep, the deleted keystone):
   inflow_full = bc_outer.apply(outflow_at_boundary.T)  # re-apply bc
   psi_face_in = inflow_full[incoming_mask, :].T        # backward seed

   # POST-O.4a.2 (1-D bare sweep):
   inflow_full = bc_outer    # incoming-ordinate slots = SEEDED inflow
   psi_face_in = inflow_full[incoming_mask, :].T        # backward seed

The seeded inflow trace is delivered by the caller as the sibling
:math:`-B` source term, in one of two ways depending on the iteration
path:

* **Driver paths** (SI / Krylov): the seed rides in
  :math:`\text{rhs.boundary}` (the boundary source
  :math:`q.\text{boundary} + B\,\psi.\text{outflow}`), delivered by
  :math:`B` as a separate coupling gain to the variadic driver (Wave O
  step O.2a; see :ref:`bc-extraction-variadic-driver`). The bare sweep's
  :meth:`InvertibleOperator._solve_timed_full_field <orpheus.sn.operator.InvertibleOperator._solve_timed_full_field>`
  seeds its boundary buffer from :math:`\text{rhs.boundary}`, **not**
  from the iterate ``initial_guess.boundary`` (the retired
  partner-flux carrier).
* **Direct loops** (direct fixed-source SI, final eigenvalue
  reconstruction sweep): the helper
  :func:`~orpheus.sn.solver._reflect_outflow_into_inflow` fills the
  inflow slots with :math:`B\,\psi.\text{outflow}` in place before the
  sweep, via the canonical
  :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`.

Both routes call the **identical**
:class:`~orpheus.sn.boundary_operator.SNBoundaryOperator` — :math:`B`
is single-sourced. For vacuum :math:`B = 0`, so the bare sweep reads a
zero inflow seed and the result is **bit-identical** to the
pre-extraction ``bc.apply`` of a vacuum law.

The full block-matrix derivation, the three design corrections (keep
the outflow defect; project :math:`B` to the inflow row; seed from
:math:`\text{rhs.boundary}`), the two delivery routes, and the O.2
forcing function all live at :ref:`bc-extraction` in
:doc:`operator_algebra` — the canonical home for the SN transport
operator algebra.

Step **O.4b** extended the bare sweep to the **2-D Cartesian
wavefront** path (both :func:`~orpheus.sn.loss_representation._sweep_jacobi`
and the 2-D matvec
:meth:`StreamingOperator._apply_2d_cartesian <orpheus.sn.operator.StreamingOperator>`):
the intra-octant ``bc.apply`` is gone there too, and the
octant-incoming edge is seeded from the given inflow trace. The
``sn_mesh.reduced is not None`` predicate that guards the dispatch now
selects the **fold shape** (1-D parallel-prefix scan vs 2-D wavefront
DAG), **not** a bare-vs-bc-in-sweep distinction — both folds are bare,
so the sweep body and the helper-guard sites cannot drift. The 2-D
interior face fluxes both folds propagate are the interior 1-cochain
:math:`C^1_{\rm int}`, with the boundary seed/absorb the typed trace
operators :math:`\iota_*` / :math:`\iota^*` (#205 Phase 5; the
``WavefrontFlux`` carrier that named the cochain is retired at
S6.4(f), the cochain now living in ``_MovingFrontier`` /
``_octant_face_cochain`` — see :ref:`wavefront-flux-cochain` in
:doc:`operator_algebra`).

.. _sn-mms-spherical-aniso-spatial-convergence:

Spherical anisotropic-ansatz MMS convergence (Phase C)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: **Retraction (2026-06-13, Issues #229 / #195).**

   The Phase-C/D attribution that follows — "the legacy angular
   closure default leaks :math:`\mathcal{O}(h^{1.3})` shape errors"
   and the rate is held back "until the angular default flips to
   ``MorelMontryAngularSweep``" — was **falsified** by the W1–W4
   root-cause program (2026-06-13).  ERR-058 (Issue #195) was the
   terminal closure-seed fix: the curvilinear *isotropic* MMS now
   converges clean :math:`\mathcal{O}(h^2)` under the existing
   default WITHOUT any default flip (the default-flip and the
   "Phase D pole-face refinement" framing below were never the
   lever).  The residual *anisotropic* floor is the
   angular half-angle-thread INTERPOLATION floor (Issue #229), which
   scales with the angular **quadrature** and is independent of the
   spatial closure, the default, and the :math:`\tau`-clamp.  The
   numbers and the structural-asymmetry reasoning preserved below are
   historical evidence; the *interpretation* is superseded by the
   comprehensive treatment at
   :ref:`sn-curvilinear-aniso-norm-reconciliation`.  W3 (Issue #229)
   removed the xfail markers and migrated the
   :eq:`sn-mms-spherical-psi` / :eq:`sn-mms-spherical-qext` labels to
   green tests; this Gate 3.1 spherical spatial-convergence test was
   **retired** (its claim re-homed on the S32 full-ladder gate — see
   the reconciliation section).

Gate 3.1 (plan §5) is the L1 MMS verification of the spatial
convergence rate on the angularly-non-trivial ansatz

.. math::

   \psi_{\text{chosen}}(r, \mu) \;=\; \frac{A(r) + B(r)\,\mu}{W},
   \qquad
   A(r) = \cos(\pi r / (2R)), \quad B(r) = r/R,

with :math:`W` a normalisation constant. The ansatz **activates**
the angular redistribution term (the linear-:math:`\mu` content
:math:`B(r)\,\mu/W` is the curvilinear sweep's hardest math), in
explicit avoidance of Mode 7 of the ``vv-principles`` skill — an
isotropic-by-construction MMS that would null the redistribution
path by ansatz design and silently miss ERR-026. The companion
isotropic ansatz :math:`\psi = A(r)/W` would still test the
spatial closure but is **insufficient** in isolation; Phase C
ships both as separate test cases (per the Phase B closeout's
"every multi-dim MMS must declare which terms it activates AND
which it nulls" rule).

With the curvilinear default
:class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
the rate stays at the pre-Phase-C
:math:`\mathcal{O}(h^{1.3})` profile. This is the diagnostic that
**the spatial-closure alignment is necessary but not sufficient**
for :math:`\mathcal{O}(h^2)` convergence; the underlying ERR-026
flux-shape drift survives Phase C's WDD sweep-frame alignment
when the angular closure default does not flip to the canonical
M–M form. Gate 3.1 is therefore marked
``@pytest.mark.xfail(strict=False)`` pending Phase D's pole-face
spatial-closure refinement (see
:ref:`sn-curvilinear-trajectory-resolvent-crosscheck-section` for the
Phase D scope summary).

The xfail is intentionally **not strict** at this gate (in
contrast to the four pre-existing ERR-026 ``xfail(strict=True)``
tripwires at the same labels). The non-strict marker reflects an
**empirical** test of an architectural prediction: if Phase C's
sweep-frame matvec accidentally moved the convergence rate past
1.9 on the legacy default (which would be the unexpected outcome
that demands a fresh investigation), the marker would flip to
``xpass`` rather than fail strictly. The strict markers stay on
the four canonical ERR-026 tripwires
(:file:`tests/sn/test_mms_curvilinear.py` and the L1 aniso file)
because those tests cover the closure status that Phase D will
actually close.

.. _sn-mms-cylindrical-aniso-spatial-convergence:

Cylindrical anisotropic-ansatz MMS convergence (Phase C)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: **Retraction (2026-06-13, Issue #229).**

   The claim below that "the Phase D fix is expected to produce a
   clean :math:`\mathcal{O}(h^2)` cylindrical MMS rate" was
   **falsified**.  The cylinder has **NO** pre-floor
   :math:`\mathcal{O}(h^2)` window at any practical quadrature: the
   angular half-angle-thread interpolation floor dominates the spatial
   error before second order can establish (measured: even
   :math:`n_\mu = 16` reaches only order 1.80 on the coarsest segment).
   The cylinder floor scales with the **azimuthal** quadrature
   :math:`n_\varphi` (NOT the polar :math:`n_\mu`), is structurally
   blocked by duplicate azimuthal :math:`\eta` that a 1-D
   :math:`\eta`-thread cannot represent, and would need a 2-D
   :math:`(\eta, \varphi)` closure (out of scope).  W3 (Issue #229)
   replaced the Gate 3.1/3.2 cylindrical spatial-rate claim with a
   verified floor-scaling test
   (``test_cyl_aniso_floor_scales_with_quadrature``) and migrated the
   :eq:`sn-mms-cylindrical-psi` / :eq:`sn-mms-cylindrical-qext` labels
   to green.  See :ref:`sn-curvilinear-aniso-norm-reconciliation` for
   the full treatment; the reasoning preserved below is history.

Gate 3.2 is the cylindrical analogue of Gate 3.1 — same ansatz
structure (linear :math:`\mu_x` content + cosine radial profile)
adapted to the cylindrical level-DAG, parametrised across LS-4 and
Product 2×4 quadratures to surface any quadrature-family-dependent
constants (Signature 4 / ERR-004). Same xfail rationale:
spatial-closure alignment is necessary but not sufficient for
:math:`\mathcal{O}(h^2)` convergence until the angular default
flips to ``MorelMontryAngularSweep``. Cylindrical Gate 1.1
**passes** under the canonical M–M angular closure (see the
Empirical Gate 1.1 finding below), so the Phase D fix is expected
to produce a clean :math:`\mathcal{O}(h^2)` cylindrical MMS rate
without requiring the spherical pole-face refinement — but the
default must flip in unison across both geometries for the
``catches("ERR-026")`` story to be coherent. Phase D will ship
both.

Gate 3.3 (angular convergence at fixed ``nx=80``, varying
``n_ordinates``) **passes** under Phase C — the spatial closure
alignment is sufficient to expose the angular discretisation as
the limiting error when the spatial discretisation is held fine.
This is the inverse signature of Gate 3.1: holding the angular
closure fixed and refining spatially saturates at the angular
discretisation floor; holding the spatial closure fixed and
refining angularly saturates at the spatial discretisation floor;
the legacy default does not produce :math:`\mathcal{O}(h^2)`
spatial because the legacy angular closure leaks
:math:`\mathcal{O}(h^{1.3})` shape errors that the spatial
discretisation cannot resolve away.

.. _sn-curvilinear-homogeneous-kinf-recovery-section:

Homogeneous-reflective k\ :sub:`∞` recovery (Phase C)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gate 4.1 verifies the eigenvalue claim using a closed-form
reference: the 2-group homogeneous reflective sphere recovers the
analytical infinite-medium eigenvalue

.. math::
   :label: sn-curvilinear-homogeneous-kinf-recovery

   \kinf
        \;=\; \rho\bigl(\mathbf{\Sigma}_a^{-1}
              \,\boldsymbol{\chi}\,\boldsymbol{\nu\Sigma_f}^{\top}\bigr)
        \;\stackrel{1\text{G}}{=}\;
        \frac{\nu\Sigma_f}{\Sigma_a}\,,

i.e., the dominant eigenvalue of the multi-group production /
removal transfer matrix on the homogeneous infinite medium
(Lewis--Miller §3.2; reduces to :math:`\nSigf/\Sigma_a` in 1-group).
The reference is computed in closed form by
:func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous`
without any spatial or angular discretisation choice; it is the
State 1A closed-form pillar in the ``algebra-of-record`` taxonomy.
The Phase C sweep-frame matvec recovers it to ``rtol ≤ 5e-4`` on the
2-group homogeneous reflective sphere; pinned at
:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_sn_spherical_homogeneous_kinf_recovery_2g`.

The clean :math:`k_\infty` recovery is **not** a contradiction of
ERR-026 staying at PARTIAL CLOSURE. The eigenvalue is shape-
independent: for a homogeneous reflective problem,
:math:`\kinf` is a material-property ratio
(:math:`\nSigf / \Sigma_a` over the volume-weighted average flux),
and the same ratio falls out of any discretisation that preserves
volume-weighted particle balance. Phase C's WDD spatial closure
preserves balance by construction (the streaming term telescopes
to surface area times average flux on a uniform mesh, and the
redistribution term integrates to zero against the volume weights
across an :math:`\mathcal{R}^4` cell). The shape-dependent ERR-026
flux-shape bug therefore drops out of the eigenvalue but persists
in the **flux shape** — exactly what
:ref:`sn-curvilinear-trajectory-resolvent-crosscheck-section` will
measure in Phase D. Gate 4.1 is therefore the **necessary** but
**not sufficient** evidence chain (per ``vv-principles`` 1-group
degeneracy rule); the sufficient chain requires structurally-
independent flux-shape evidence from Phase D.

.. _sn-curvilinear-trajectory-resolvent-crosscheck-section:

Trajectory-resolvent cross-check (Phase C → Phase D)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gate 4.2 is the **flux-shape cross-check** against the
structurally-independent trajectory_resolvent Green's-function
reference (the Peierls Variant α State 1B semi-analytical pillar
in the ``algebra-of-record`` taxonomy).  The contractual claim is

.. math::
   :label: sn-curvilinear-trajectory-resolvent-crosscheck

   \bigl\|\phi^{\,\text{SN}}_h(r)
        \;-\; \phi^{\,\text{traj.res.}}(r)\bigr\|_{\infty}
        \;\le\; 5\times 10^{-4}
        \quad
        \text{on the 5 P0 curvilinear snapshots,}

with :math:`\phi^{\,\text{SN}}_h` the SN flux at the snapshot's
:math:`n_x` and :math:`\phi^{\,\text{traj.res.}}` the
trajectory-resolvent reference flux at the same radii.  The bare
function entry points cover the 5 P0 deleted curvilinear regression
snapshots:

.. list-table:: trajectory_resolvent reference coverage
   :header-rows: 1
   :widths: 35 50 15

   * - Snapshot
     - Bare entry point
     - Precision
   * - ``sphere_2g_homogeneous_dd_n20``
     - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mg`
     - :math:`k_\infty` exact via V_α1 identity
   * - ``sphere_2g_3reg_dd_n40``
     - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr`
     - MR↔MG reduction ``rtol=1e-9``
   * - ``cyl_1g_homogeneous_LS4_dd_n20``
     - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder.solve_greens_function_cylinder`
     - V_α1_cyl exact; Sood Ua-1-O-CY vacuum ``8.5e-6``
   * - ``cyl_1g_homogeneous_product_dd_n20``
     - same as above (different SN quadrature)
     - same
   * - ``cyl_2g_3reg_LS4_dd_n40``
     - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder.solve_greens_function_cylinder_mr`
     - MR↔MG K=3 2G ``rtol=1e-9``

The 6th snapshot (``sphere_2g_p1_aniso_dd_n20``) routes to Gate 4.1
because :math:`P_1` anisotropic eigenvalue is still
shape-independent for a homogeneous reflective problem.

The test placeholder at
:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_sn_cylinder_homogeneous_vs_trajectory_resolvent_1g`
is marked SKIP pending Phase D's pole-face spatial-closure
refinement. The placeholder is **structurally important**: it pins
the names of the bare entry points so the Phase D session knows
exactly where the reference comes from. The structurally-
independent cross-check is the load-bearing flux-shape evidence
for ERR-026 → CLOSED — without it the closure narrative would rest
on :math:`k_\infty` agreement alone, which is degenerate
(``vv-principles`` 1-group degeneracy rule applied to homogeneous
multi-group: any discretisation that preserves balance gets
:math:`k_\infty` right, so :math:`k_\infty` alone is not flux-shape
evidence). The Phase D L1 acceptance criterion is rtol :math:`\le
5 \times 10^{-4}` against the trajectory_resolvent reference on
each of the 5 P0 snapshots — relaxed from rtol :math:`\le 10^{-9}`
because SN nx-discretisation dominates the error budget at the
practical mesh refinement levels.

Phase B's ``pole-mm-recurrence`` label (:eq:`pole-mm-recurrence`)
**gains a tests edge transitively** through the Phase D fix: once
``MorelMontryAngularSweep`` becomes the default and Gates 3.1 / 3.2
xpass, the canonical Hébert §3.9.4 angular recurrence is exercised
by the apply matvec and pinned by an L1 test chain. Through Phase C
the label remains tested only via the Phase B foundation suite
(:file:`tests/sn/spatial/test_pole_angular_closure.py`); the L1
upgrade is Phase D's responsibility.

Empirical Gate 1.1 finding: spherical-vs-cylindrical structural asymmetry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Gate 1.1 (the canonical curvilinear bug-class L0 diagnostic) is
the per-ordinate flat-flux residual probe: for :math:`\psi` constant
in space and per-ordinate on a reflective-BC homogeneous curvilinear
problem, the apply matvec must produce :math:`(T\psi)_{n,i,g} =
\Sigma_t \cdot \psi_{n,i,g}` per ordinate to ``rtol=1e-12``
(:math:`\Sigma_t = 0` reduces this to bit-zero to ``atol=1e-13``).
Parametrisation: spherical + cylindrical × 4 quadrature variants ×
2 group counts × 3 nx values × 2 :math:`\Sigma_t` values × 3
pole-closure strategies — 288 combinations under the strict
specification.

The empirical outcome decides the **default flip** per the user's
explicit constraint 7 (the "do not flip without empirical
evidence" sequencing). The decisive subset is the (geometry,
pole-closure) crosstab on the canonical Hébert M–M angular closure
strategy
(:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`):

.. list-table:: Empirical Gate 1.1 outcome (Phase C, 2026-05-12)
   :header-rows: 1
   :widths: 20 25 25 30

   * - Geometry
     - Pole closure
     - :math:`\Sigma_t = 0`
     - :math:`\Sigma_t = 0.5`
   * - Sphere
     - ``LegacyTauSymmetricInterpolation``
     - PASS
     - PASS
   * - Sphere
     - ``BaileyFlatFluxRedist``
     - PASS
     - PASS
   * - Sphere
     - ``MorelMontryAngularSweep``
     - **FAIL**
     - **FAIL**
   * - Cylinder
     - ``LegacyTauSymmetricInterpolation``
     - PASS
     - PASS
   * - Cylinder
     - ``BaileyFlatFluxRedist``
     - PASS
     - PASS
   * - Cylinder
     - ``MorelMontryAngularSweep``
     - **PASS**
     - **PASS**

The asymmetry is **structural**: spherical-MMS fails;
cylindrical-MMS passes. The mechanism is the interaction between
the pole-face WDD initial condition (the Carlson seed
:math:`\psi^{\text{face}}_{\text{in}}(\text{pole}) =
\psi^{\text{cell}}[0]`) and the canonical Hébert §3.9.4 angular
recurrence's half-angle face flux at the pole:

* **Cylindrical case** has **per-level :math:`\alpha`-dome
  telescoping**. Each :math:`\mu`-level has its own
  :math:`\alpha_{n+1/2}` recurrence with its own pair of
  starting-direction face fluxes (:math:`\mu_x = -1` inward zero
  weight + :math:`\mu_x = +1` outward zero weight). The
  half-angle face flux discrepancy from the M–M recurrence's
  Carlson seed is absorbed by the level's own :math:`\alpha`-dome
  closure across its azimuthal ordinates — the level integrates to
  zero in the angular flux moment that drives the redistribution,
  and the level-to-level coupling at the level boundaries is
  through pole-azimuth-degenerate ordinates that carry no spatial
  flow. The cancellation is automatic.

* **Spherical case** has **no equivalent telescoping**. The
  spherical pole-face is a single point (the centre :math:`r=0`),
  not a level boundary, and the entire :math:`\alpha`-cascade
  meets there. The half-angle face flux discrepancy from the M–M
  recurrence accumulates across the full ordinate set rather than
  cancelling per level. The Carlson seed at the pole-**face**
  resolves the **outer** sweep direction; the M–M angular
  recurrence's starting-direction face flux at :math:`\mu = -1` is
  a separate seed that must be consistent with the spatial
  closure. The two seeds are not jointly consistent under Phase C
  alone — that is the Phase D scope.

The structural asymmetry is one of the **load-bearing
intellectual findings** of Phase C. The plan §1 had predicted
"sweep-frame architecture more likely to make MMS angular closure
viable (because spatial closure is now WDD throughout, matching
what MMS expects)", but the empirical probe revealed that the
spherical pole is **doubly singular** in a sense the cylindrical
case is not: both the angular :math:`\alpha`-cascade and the
spatial WDD recurrence converge to the same singular point, and
both need consistent starting-direction seeds. The Phase D
follow-up (a Carlson-style **coupled** pole sweep where the
outward-ordinate pole-face initial condition is determined by the
inward-ordinate pole-face propagation, not chosen independently)
is the architectural prescription. This is the symmetry condition
at :math:`r=0` written into the SN discretisation.

Per the user's explicit constraint 7, the default flip to
``MorelMontryAngularSweep`` is **DEFERRED to Phase D**
(`Issue #192 <https://github.com/deOliveira-R/ORPHEUS/issues/192>`_).
Cylindrical-MMS Gate 1.1 PASS is the strong positive signal for
the Phase D fix: the cylindrical structure is already shape-
correct under the canonical Hébert closure with Phase C's
sweep-frame architecture; the Phase D additional refinement
targets the spherical pole-face only and inherits cylindrical
behaviour for free.

ERR-026 closure status (Phase C — sweep-frame matvec aligned)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Phase C ships the architectural alignment — sweep-frame matvec
with WDD spatial closure + BC trace law at the boundary edge +
retired Phase A
:class:`~orpheus.sn.spatial.boundary_face_flux.BoundaryFaceFlux`
Protocol — but per the empirical Gate 1.1 finding above the
curvilinear default stays
:class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
and the four ``xfail-strict`` curvilinear MMS tripwires STAY xfail.
ERR-026 remains at **PARTIAL CLOSURE** through Phase C — the
spatial-closure architecture is aligned (the load-bearing Phase B
precondition); the pole-face spatial-closure refinement is the
Phase D scope.

Verification gate summary
'''''''''''''''''''''''''

.. list-table:: Phase C verification gates
   :header-rows: 1
   :widths: 8 50 22 20

   * - Gate
     - Description
     - Status
     - Pinned at
   * - 1.1
     - Per-ordinate flat-flux residual
     - PASS (Legacy + BFF on both geometries); PASS (MMS on cyl); xfail (MMS on sphere)
     - ``test_phase_c_gates.py``
   * - 1.2
     - apply determinism via ``np.array_equal``
     - PASS
     - ``test_phase_c_gates.py``
   * - 1.3
     - apply ↔ apply_transpose reciprocity
     - PASS (``rtol=1e-12``)
     - ``test_phase_c_gates.py``
   * - 1.4
     - apply linearity (precondition)
     - PASS (``rtol=1e-13``)
     - ``test_phase_c_gates.py``
   * - 1.5
     - BC trace contract honoured by matvec
     - PASS
     - ``test_phase_c_gates.py``
   * - 2.1
     - 5 Cartesian regression snapshots bit-identical
     - PASS (``rtol=1e-12``)
     - ``test_dd_regression.py``
   * - 2.2
     - Phase B 28 foundation tests
     - PASS
     - ``test_pole_angular_closure.py``
   * - 2.3
     - Phase B 5 L1 flat-flux-identity tests
     - PASS
     - ``test_pole_closure_flat_flux_identity.py``
   * - 2.4
     - 21 Phase A ``BoundaryFaceFlux`` tests retired
     - DONE
     - (file deleted)
   * - 3.1
     - Spherical anisotropic MMS spatial convergence
     - xfail (ERR-026 PARTIAL)
     - ``test_phase_c_mms.py``
   * - 3.2
     - Cylindrical anisotropic MMS spatial convergence
     - xfail (ERR-026 PARTIAL)
     - ``test_phase_c_mms.py``
   * - 3.3
     - Angular convergence at fixed nx
     - PASS
     - ``test_phase_c_mms.py``
   * - 4.1
     - :math:`k_\infty` recovery on 2G reflective sphere
     - PASS (``rtol < 5e-4``)
     - ``test_phase_c_crosscheck.py``
   * - 4.2
     - trajectory_resolvent flux-shape cross-check
     - SKIP (Phase D)
     - ``test_phase_c_crosscheck.py``

Phase D scope (shipped 2026-05-12 — Carlson coupled-pole seed)
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

.. note:: **Read this Phase D / Phase F narrative as HISTORY (Issue
   #195 CLOSED 2026-06-12).**  The Carlson coupled-pole *seed concept*
   below is correct and survives; but two of Phase D's terminal
   decisions were later reverted by the ERR-058 fix: (i) the
   curvilinear ``inner_solver`` default flip to ``"krylov"`` was undone
   (curvilinear now defaults to ``"source_iteration"``, SI
   :math:`\equiv` Krylov bit-identical post-unification); and (ii) the
   "magnitude scope OPEN / pre-asymptotic transient" framing was
   falsified — the error PLATEAUED, the dominant defect was the
   *angular* closure seed, and ERR-058 made the isotropic MMS a clean
   :math:`\mathcal{O}(h^2)` ladder.  The
   :class:`CarlsonInwardSweep` *half-angle* seed described here is also
   superseded as the default by
   :class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`.
   The production resolution is at
   :ref:`sn-err-058-closure-seed-closeout`; this section's per-step
   claims are tombstoned inline where they bear on those decisions.

Tracked at `Issue #192
<https://github.com/deOliveira-R/ORPHEUS/issues/192>`_; the
Hébert §3.9.4 inward sweep implementation is `Issue #193
<https://github.com/deOliveira-R/ORPHEUS/issues/193>`_; the
``@pytest.mark.verifies(...)`` wiring for the new equation
labels is `Issue #194
<https://github.com/deOliveira-R/ORPHEUS/issues/194>`_; the
remaining pre-asymptotic-magnitude open question is `Issue #195
<https://github.com/deOliveira-R/ORPHEUS/issues/195>`_. The
shipped deliverables flip ERR-026's **identity-and-rate** scope to
CLOSED while keeping the **magnitude** scope open per
:ref:`sn-phase-d-err-026-closure-narrative` below:

1. **M-M half-angle seed refinement.** A canonical Hébert §3.9.4
   Eqs. (3.432)–(3.435) inward :math:`\mu = -1` sweep seeds the
   M-M angular recurrence's ``psi_half_left`` — replacing the
   hardcoded zero that Phase B had baked in. See
   :ref:`sn-phase-d-carlson-coupled-pole-sweep` for the full
   derivation, including the diagnostic finding that **the seed
   lives in the M-M angular recurrence, not in the WDD spatial
   pole-face IC**.
2. **Default flip.** :attr:`SNMesh.pole_angular_closure` default
   :class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
   → :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`;
   :class:`~orpheus.sn.solver.SNSolver` curvilinear default
   ``"source_iteration"`` → ``"krylov"``. See
   :ref:`sn-phase-d-default-flips`.
3. **Gate 1.1 sphere MMS PASS.** All 4 ``MorelMontryAngularSweep``
   parametrised cases (sphere × cylinder × :math:`\Sigma_t \in
   \{0, 0.5\}`) **xpass** on the per-ordinate flat-flux residual
   probe. See :ref:`sn-phase-d-gate-1-1-empirical`.
4. **Snapshot regeneration deferred.** The 11 DD regression
   snapshots remain bit-identical under Phase D — the SI/sweep
   path is the snapshot generator, and the SI path uses the
   sweep (not the apply matvec), so the Phase D default flip
   does NOT disturb them. Regeneration under a Phase D
   :meth:`InvertibleOperator.apply`-driven Krylov path is
   carried at Issue #195.
5. **Gate 1.5 strengthened (capture-and-compare).** The §16A.3
   BC trace contract now has a stricter parametrised test that
   independently reconstructs the WDD-propagated outflow trace
   and asserts the captured BC apply input matches to
   ``rtol=1e-14``. See
   :ref:`sn-phase-d-gate-1-5-capture-and-compare` and the BC
   companion section :ref:`bc-phase-d-two-bc-applies-per-matvec`.
6. **Marker partial removal — deferred.** The 4 ``xfail-strict``
   ERR-026 tripwires stay through Phase D Step 3 — they will
   ``xpass`` under the new default but require the Step 5
   marker-removal commit. ERR-026 stays at PARTIAL CLOSURE
   pending Issue #195 (pre-asymptotic-magnitude convergence).

.. _sn-phase-d-carlson-coupled-pole-sweep:

Phase D Carlson coupled-pole sweep (Issue #168 Phase D)
-------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   * Phase D (commit landed 2026-05-12 on
     ``refactor/sn-operator-algebra``) closes the structural bug
     in :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
     by replacing the hardcoded ``psi_half_left = 0`` seed with
     the canonical Hébert §3.9.4 Eqs. (3.432)–(3.435) inward
     :math:`\mu = -1` sweep output.
   * The seed lives in the **M-M angular recurrence**
     (``orpheus/sn/spatial/pole_angular_closure.py:411`` —
     ``_mm_weighted_angular_recurrence_single_level``), **NOT**
     in the WDD spatial pole-face initial condition the
     :ref:`Phase C plan <sn-curvilinear-trajectory-resolvent-crosscheck-section>`
     proposed.  The diagnostic memo at
     ``.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md``
     empirically falsified intervention ``[A]`` (WDD pole-face
     replacement) and confirmed intervention ``[B]`` (M-M
     half-angle seed replacement).
   * Architectural choice is **Option α (composition)**: the
     seed lives as a
     :class:`~orpheus.sn.spatial.psi_half_angle_seed.PsiHalfAngleSeed`
     strategy field on :class:`MorelMontryAngularSweep`, not as a
     sibling Protocol on :class:`SNMesh`. The Legacy / Bailey
     closures have no ``psi_half_left`` variable to seed; a
     sibling Protocol would force every consumer to handle an
     irrelevant Protocol.
   * The **L = 0 isotropic-only** assumption is load-bearing: the
     apply matvec's :math:`L` operator currently carries only
     :math:`\Sigma_t \psi` (scattering is composed externally
     via a separate operator).  A future refactor that moves
     scattering INTO :math:`L` MUST extend the moment-folded
     source in :eq:`hebert-3-432-source` to include
     :math:`\ell \ge 1` terms.

The Hébert §3.9.4 equations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Hébert §3.9.4 (pp. 141–144 of [Hebert2009]_) opens the sphere
difference relations at Eq. (3.418) (angularly-integrated
divergence form), introduces the :math:`\alpha`-recursion
:eq:`bailey-dome-recursion` and the cell-balance with
redistribution divisor :math:`\Delta S_i / (2\,\mathcal{W}_n)` at
Eq. (3.428), and then specialises to the auxiliary starting
direction :math:`\mu = -1`.  At this direction the
angular-redistribution coefficient :math:`(1 - \mu^2)` is
**identically zero**, so the streaming–collision balance
decouples from the :math:`\alpha`-cascade and reduces to a plain
DD inward recurrence in radius.

The continuous form at :math:`\mu = -1` (Hébert Eq. (3.432)) is

.. math::
   :label: hebert-3-432

   -\frac{\partial}{\partial r}\,\phi_{-1/2}(r)
   \;+\; \Sigma(r)\,\phi_{-1/2}(r)
   \;=\; \sum_{\ell=0}^{L}
         \frac{2\ell + 1}{2}\,Q_\ell(r)\,P_\ell(-1).

The subscript :math:`-1/2` is Hébert's half-integer index for the
auxiliary starting ordinate — it labels the **inward
zero-weight** direction that sits one half-step above
:math:`\mu = -1` in the :math:`\alpha`-cascade, not a physical
ordinate at :math:`\mu = -0.5`.  The right-hand side is the
Legendre expansion of the scattering source :math:`Q` evaluated
at :math:`\mu = -1`, where :math:`P_\ell(-1) = (-1)^\ell`.

For an **isotropic** operator (:math:`L = 0`, the current ORPHEUS
apply matvec) the source collapses to

.. math::
   :label: hebert-3-432-source

   \bar Q_i \;=\; \sum_\ell \frac{2\ell+1}{2}\,Q_\ell(r_i)\,(-1)^\ell
   \;\;\xrightarrow{L=0}\;\;
   \frac{1}{2}\,\Sigma_t(r_i)\,\phi_0(r_i),

where :math:`\phi_0` is the scalar-flux Legendre :math:`\ell = 0`
moment of the input :math:`\psi`.  See the **L = 0 isotropic-only**
WARNING in the
:class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
class docstring for the failure mode if scattering moves into the
operator.

Discretising Eq. :eq:`hebert-3-432` on a sub-mesh of cell width
:math:`\Delta r_i` gives the DD cell-balance Hébert Eq. (3.433):

.. math::

   -\bigl(\bar\phi_{i+1/2} - \bar\phi_{i-1/2}\bigr)
   \;+\; \Delta r_i \cdot \Sigma_i \cdot \bar\phi_i
   \;=\; \Delta r_i \cdot \bar Q_i,

with Hébert's typographic conventions

.. math::

   \bar\phi_i \;\equiv\; \phi_{1/2,\,i}, \qquad
   \bar Q_i \;\equiv\; Q_{1/2,\,i}, \qquad
   \Delta r_i \;=\; r_{i+1/2} - r_{i-1/2}.

The negative sign on the streaming jump comes from
:math:`\mu = -1 < 0` — particles travel **inward**, so the
discrete jump is :math:`-(\phi_{i+1/2} - \phi_{i-1/2})`.
**Critically**, no :math:`\alpha`-redistribution divisor appears
in this balance because :math:`(1 - \mu^2) = 0` at the endpoint.
This is the entire reason Hébert can solve the :math:`\mu = -1`
sweep in closed form with a plain DD recurrence: the coupled
angular cascade is decoupled at the starting direction.

Combining the DD auxiliary relation
:math:`\phi_{n,i} = \frac{1}{2}(\phi_{n,i-1/2} + \phi_{n,i+1/2})`
specialised to the :math:`-1/2` ordinate with the balance and
solving for :math:`\bar\phi_i` in terms of the known
outgoing-face value :math:`\bar\phi_{i+1/2}` (further from the
centre — known because we sweep **inward** from the outer BC)
yields Hébert Eq. (3.434):

.. math::
   :label: hebert-3-434

   \bar\phi_i \;=\; \frac{\Delta r_i \cdot \bar Q_i
                            + 2 \cdot \bar\phi_{i+1/2}}
                          {\Delta r_i \cdot \Sigma_i + 2}.

Stepping inward to the next face uses the textbook DD auxiliary
relation rearranged (Hébert Eq. (3.435)):

.. math::
   :label: hebert-3-435

   \bar\phi_{i-1/2} \;=\; 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}.

The pair :eq:`hebert-3-434`–:eq:`hebert-3-435` IS the spatial
recurrence.  Together they realise a tridiagonal-style inward
sweep on the radial mesh: outer face :math:`\rightarrow` cell
centre :math:`\rightarrow` inner face :math:`\rightarrow` next
cell centre :math:`\rightarrow \ldots \rightarrow` pole face
:math:`\bar\phi_{1/2}` at :math:`r = 0`.

.. note::

   The three labels :eq:`hebert-3-432-source`,
   :eq:`hebert-3-434`, :eq:`hebert-3-435` are also declared in the
   :mod:`~orpheus.sn.spatial.psi_half_angle_seed` module docstring
   (the canonical algebra-of-record).  Each :math:`:label:` is
   unique across the documentation graph; the Sphinx page is the
   **presentation layer** for the equations the code module owns
   as source-of-truth.  The
   ``@pytest.mark.verifies("hebert-3-43X")`` wiring on the L0
   algebraic-identity tests in
   :file:`tests/sn/spatial/test_psi_half_angle_seed.py` is tracked
   at Issue #194; without that wiring the labels appear in the V&V
   audit as "documented but not tested" (orphan labels).

Why :math:`\mu = -1` is the natural starting direction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The M-M angular closure on sphere is a per-cell
:math:`\alpha`-cascade that **couples** the angular flux across
ordinates within one spatial cell: ordinate :math:`n` reads
:math:`\alpha_{n-1/2}` from the previous (more-inward-:math:`\mu`)
ordinate.  To start the cascade at the smallest-:math:`\mu`
ordinate, one needs a value for :math:`\alpha_{1/2}` AND for the
angular-edge flux :math:`\phi_{1/2,i}` at that seed half-integer.

The :math:`\alpha_{1/2} = 0` seed is **free**: it comes from
:math:`1 - \mu^2` evaluated at :math:`\mu = -1`, i.e.,
"The first value :math:`\alpha` is equal to :math:`1 - (-1)^2 = 0`"
(text below Hébert Eq. (3.422)).  That handles the *angular*
half of the problem.

The flux value :math:`\phi_{1/2,i}`, however, is NOT free.  It
is the **spatial flux profile** at the auxiliary starting
direction, and it must be solved for as a function of position
:math:`i` along the radial mesh.  Eqs. :eq:`hebert-3-432` through
:eq:`hebert-3-435` provide exactly that spatial solve.

At :math:`\mu = -1` the sphere streaming operator collapses to
pure radial divergence **without** the angular-redistribution
coupling.  As Hébert writes (p. 143):

   *"We observe that these directions correspond to particles
   entering the external surface and moving toward the central
   axis with* :math:`\mu = -1`. *The angular redistribution term
   vanishes on these points so that Eq. (3.164) simplifies to
   [Eq. (3.432)]."*

This is the **only** direction on the unit interval :math:`[-1, 1]`
where the spatial 1D-sphere problem reduces to a closed-form
linear recurrence in radius alone, without an inner angular
solve.  Picking any intermediate :math:`\mu` would leave the
coupling term active and re-introduce the cascade
chicken-and-egg.  See also
:ref:`sn-phase-d-pomraning-structural-singularity` for the
deeper structural reason :math:`\mu = \pm 1` is the only
admissible starting direction in any curvilinear geometry.

Why "zero-weight"
~~~~~~~~~~~~~~~~~

In an :math:`N`-point Gauss–Legendre quadrature on :math:`[-1, 1]`
the endpoints :math:`\mu = \pm 1` are **not** base points (the
polynomial is approximated by interior nodes only).  They have no
quadrature weight, hence "zero-weight" — the flux value at
:math:`\mu = -1` does NOT contribute to any
:math:`\sum_n \mathcal{W}_n \phi_n` integral that builds the scalar
flux moments.

The :math:`\mu = -1` ordinate is therefore a **purely auxiliary
numerical construct**: its flux values exist for the sole purpose
of seeding the :math:`\alpha`-cascade for the finite-weight
ordinates that follow.  After the cascade is initialised, the
angular-edge values :math:`\bar\phi_{i\pm 1/2}` are discarded;
only the **cell-centred values** :math:`\bar\phi_i \equiv
\phi_{1/2,i}` are kept (Hébert, p. 143, between Eqs. (3.435) and
(3.436)).  Those cell-centred values feed the finite-weight
ordinates' cell-balance Eq. (3.436) via the
:math:`(\alpha_{n-1/2} + \alpha_{n+1/2})\,\phi_{n-1/2,i} /
(2\,\mathcal{W}_n)` redistribution term, with
:math:`\phi_{n-1/2,i} = \phi_{1/2,i}` at the first
finite-weight ordinate :math:`n = 1`.

The flat-:math:`\psi` algebraic verification trace
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Phase D hypothesis is: *for a flat angular flux*
:math:`\psi_{\text{cell}} = C` *across all cells, the
inward sweep returns* :math:`\bar\phi_{1/2} = C`.  The algebra
verifies this in closed form.

Take a homogeneous problem with constant :math:`\Sigma_t = \Sigma`
and source :math:`\bar Q_i` constructed so the consistent fixed
point is :math:`\bar\phi_i = C` everywhere.  Specialising
Eq. :eq:`hebert-3-432-source` to :math:`L = 0` and applying the
flat-:math:`\psi` ansatz gives the consistent source
:math:`\bar Q = \frac{1}{2} \Sigma \cdot 2C = \Sigma \cdot C` (the
:math:`\phi_0` integral over flat unit-:math:`\psi` against GL
weights summing to 2 returns :math:`2C`; lumped into the discrete
:math:`\bar Q_i = \Sigma \cdot C`).

Substituting into Eq. :eq:`hebert-3-434` with inductive hypothesis
:math:`\bar\phi_{i+1/2} = C`:

.. math::

   \bar\phi_i
   \;=\; \frac{\Delta r \cdot \Sigma \cdot C + 2C}
              {\Delta r \cdot \Sigma + 2}
   \;=\; C \cdot \frac{\Delta r \cdot \Sigma + 2}
                     {\Delta r \cdot \Sigma + 2}
   \;=\; C.

Eq. :eq:`hebert-3-435` then gives
:math:`\bar\phi_{i-1/2} = 2C - C = C`.  The recurrence is
self-similar: every face and cell value stays at :math:`C`.  Hence
:math:`\bar\phi_{1/2}(r = 0) = C` for flat :math:`\psi` on the
consistent flat source — **the hypothesis holds**.

This trace establishes the Phase D fix as a **closed-form
analytical reference** in the
:doc:`algebra-of-record </development>` State-1A pillar sense: the
identity :math:`(L \cdot \psi_{\text{flat}})_{n,i,g} = \Sigma_t
\cdot \psi_{n,i,g}` is verifiable by exact algebra on the discrete
operator, no numerical quadrature required.  The L0 foundation test
:func:`tests.sn.spatial.test_psi_half_angle_seed.test_l0_flat_psi_reflective_identity_at_C_eq_1`
pins this identity at machine precision (``rtol=1e-13``).

The corrected injection-point story
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The single largest **architectural correction** of Phase D is
where the canonical inward-sweep output is injected.  The Phase D
plan (and the literature memo's §7 implementation note) routed
the inward-sweep result :math:`\bar\phi_i` into the **WDD
spatial pole-face initial condition** at
:func:`~orpheus.sn.operator.transport_operator_matvec_spherical`'s
``psi_face_in`` initialisation — the very same site the
:ref:`sn-curvilinear-trajectory-resolvent-crosscheck-section` discussion
identified as the Phase C Carlson seed location.

The numerics-investigator diagnostic
(:file:`tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`)
falsified that hypothesis empirically.  Four interventions tested
against the M-M failing configuration on the flat-:math:`\psi`
probe:

.. list-table:: Phase D injection-point intervention sweep (Σ_t = 0.5)
   :header-rows: 1
   :widths: 8 40 30 22

   * - Probe
     - What it changes
     - Site
     - max\|residual\|
   * - ``[A]``
     - Carlson seed for WDD ``psi_face_in``
     - ``operator.py:738``
     - **1.89e+01 FAIL** (unchanged)
   * - ``[B]``
     - Carlson seed for M-M half-angle ``ψ_{1/2,i}``
     - ``pole_angular_closure.py:411``
     - **1.78e-15 PASS**
   * - ``[C]``
     - BOTH ``[A]`` + ``[B]``
     - both
     - **1.78e-15 PASS** (no extra effect)
   * - ``[D]``
     - M-M half-angle ``ψ_{1/2,i}`` = cell-centre value
     - ``pole_angular_closure.py:411``
     - **1.78e-15 PASS** (degenerate)

Reading the table:

* ``[A]`` confirms the WDD spatial pole-face IC is **not** what's
  wrong.  The Phase C
  ``psi_face_in = fi[:, outgoing_mask, 0, 0]`` Lewis–Miller
  cell-centre seed is already structurally equivalent to the
  Carlson inward-sweep output **on flat ψ** — both equal
  :math:`\psi_{\text{cell}}[0]` in that limit.  Replacing the
  WDD seed changes nothing.
* ``[B]`` is the canonical Carlson intervention: feeding
  :math:`\bar\phi_i` into the M-M recurrence's ``psi_half_left``
  closes the residual to machine precision.
* ``[D]`` is the **falsification check**: on the flat-:math:`\psi`
  reflective probe ``[B]`` and ``[D]`` coincide because the
  inward sweep returns :math:`\bar\phi_i \equiv \psi_{\text{cell}}`
  exactly.  The probe **cannot distinguish** the two.

To prove the Carlson seed is canonical (not merely coincidentally
correct), the diagnostic includes a vacuum-BC structural
independence cross-check.  On a vacuum-BC probe the inward sweep
returns a non-trivial spatial profile

.. math::

   \bar\phi_i \;=\; (0.613, 0.572, 0.527, 0.478, 0.423, 0.362,
                     0.295, 0.220, 0.138, 0.048),

distinctly **not** equal to the cell-centred flat
:math:`\psi_{\text{cell}} = \mathbf{1}`.  The two seeds differ by
up to 0.95 in absolute value, and the resulting operator residuals
differ by max-abs 7.31 — the Carlson seed ``[B]`` is mathematically
distinct and quantitatively superior to the degenerate
broadcast-cell-centre seed ``[D]``.  This is the
**structural-independence evidence** that pins the Phase D fix as
canonical, not as a coincidental match on a degenerate probe.

The pinning test for this structural distinction is
:func:`tests.sn.spatial.test_psi_half_angle_seed.test_l1_carlson_vs_zero_seed_differ_on_vacuum_bc_probe`
— a max-abs difference :math:`> 0.05` between Carlson and Zero
seed residuals on a vacuum-BC probe.  Without this test a future
regression that replaced the Carlson sweep with a naive
broadcast-cell-centre would pass every flat-:math:`\psi`
reflective test silently.

The bug Phase B baked in
~~~~~~~~~~~~~~~~~~~~~~~~

The pre-Phase-D production code at
``orpheus/sn/spatial/pole_angular_closure.py:411`` carried the
hardcoded zero seed:

.. code-block:: python

   psi_half_left = np.zeros((ng, nx), dtype=psi_level.dtype)
   for m in range(M):
       tau_m = tau_level[m]
       psi_half_right = (
           psi_level[:, m, :] - (1.0 - tau_m) * psi_half_left
       ) / tau_m
       redist[:, m, :] = (
           dAw_level[:, m].reshape(1, nx)
           * (alpha_level[m + 1] * psi_half_right
              - alpha_level[m] * psi_half_left)
           / volume.reshape(1, nx)
       )
       psi_half_left = psi_half_right

The Phase B docstring justified the zero seed as: *"for the
forward apply matvec we adopt* :math:`\phi_{1/2,i} = 0`, *the
unique choice that makes the recursion's seed consistent with*
:math:`\alpha_{1/2} = 0` *and that the sweep converges to under
fixed-point iteration."*  This reasoning is wrong — the
:math:`\alpha_{1/2} \psi_{1/2}` product vanishes regardless of
:math:`\psi_{1/2}` because :math:`\alpha_{1/2} = 0`, but the seed
ALSO enters the **denominator-propagation chain**: every
subsequent half-angle face flux
:math:`\psi_{m+1/2,i,g}` depends on :math:`\psi_{m-1/2,i,g}`
recursively, and the chain inherits the seed through the M-M
weighting :math:`(1 - \tau_m)`.  Setting the seed to zero when
Hébert's structural form says :math:`\psi_{1/2,i,g} =
\bar\phi_{1/2,i}` (the inward-sweep output) is a **wrong term
initialisation** — Mode 3 in the
``vv-principles`` 6-failure-mode taxonomy (see ``error_catalog.md``
ERR-026 entry).

How the wrong seed survived Phase B
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The zero seed survived Phase B's L1 flat-flux-identity test
(:file:`tests/sn/l1_analytical/test_pole_closure_flat_flux_identity.py`)
because that test compared the three closures (Legacy / BFF /
M-M) **against each other on flat ψ**, NOT against the closed-form
fixed-point identity :math:`L \cdot \psi = \Sigma_t \cdot \psi`.
All three closures collapse to the same wrong-but-internally-
consistent value on flat :math:`\psi`, so cross-comparison passes
while the absolute closed-form check would have caught it
immediately.

The cylindrical case ALSO carries the zero seed in production but
:ref:`Cylindrical Gate 1.1 <sn-phase-d-gate-1-1-empirical>`
**passes** empirically.  The mechanism is per-:math:`\mu`-level
:math:`\alpha`-dome telescoping: each level's
:math:`\alpha`-cascade ends at :math:`\alpha = 0` by antisymmetry
at the level edges, so the wrong ``psi_half_left = 0`` seed
cancels cleanly per level via the level boundary closure.  The
sphere cascade has no equivalent telescoping — a wrong seed
propagates directly to a wrong fixed point.  Phase D's fix
updates the cylindrical path too for **structural alignment with
the canonical form** (architectural correctness), but cylindrical
behaviour is empirically a regression-stability check, not a new
PASS.

.. _sn-phase-d-pomraning-structural-singularity:

Pomraning structural-singularity cross-reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pomraning (1989) [Pomraning1989]_ frames the curvilinear pole
problem as **geometric**: :math:`r = 0` is structurally singular
in any curvilinear streaming operator because the
angular-derivative coefficients in the streaming term (the
:math:`(1 - \mu^2)/r` factor in the sphere streaming operator)
contain :math:`1/r`.  At :math:`r = 0` the coefficient diverges;
the natural discretisation must somehow handle this.  In his words
(p. 339, right column):

   *"It was pointed out that if the bounding surface of the
   system is used as one of the coordinate surfaces and one
   considers a family of nonintersecting surfaces that starts
   with the bounding surface and progresses inward to fill the
   system, then these surfaces will eventually shrink to a
   surface with a zero area, namely a line or a point. ... A
   special case of this elliptical example is a sphere, where
   the innermost surface is simply a point.  Hence, in general
   there will exist points on the innermost surface where the
   coefficients of the angular derivatives in the streaming term
   are infinite, since these coefficients contain the reciprocal
   of the radii of curvature ... Prime examples of such singular
   points are found in the usual spherical and cylindrical
   geometry formulations where* :math:`1/r` *terms are extant and
   the attendant difficulties are well known, particularly in
   numerical treatments."*

The naive engineering response would be **extrapolation**: pick
:math:`\psi_{\text{face}}(r = 0)` by fitting a polynomial in
:math:`r` through nearby interior cells.  This is what an
incautious starting heuristic does; it is also what produces
the M-M wrong fixed point ERR-026 diagnoses.

The Carlson coupled-pole response is **canonical** because it
sidesteps the singularity entirely: at the auxiliary direction
:math:`\mu = -1` the singular :math:`(1 - \mu^2)/r` term is
**identically zero** (the numerator vanishes), so the spatial
sweep at this direction sees **no singularity at all**.  The
equation tells the discretisation what
:math:`\bar\phi_{1/2}(r = 0)` should be — there is no need to
guess.  The cost is that the :math:`\mu = -1` sweep must be
solved first, then its result used as the seed for the cascade
at finite-weight ordinates (where :math:`(1 - \mu^2) > 0` and
the singularity would otherwise be felt).  This is exactly the
price Pomraning warns about: "difficulties must be dealt with".
The Carlson construction deals with it by **exploiting** the
singularity's vanishing at :math:`\mu = \pm 1` rather than
trying to regularise it at intermediate :math:`\mu`.

Option α: composition over sibling Protocol
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The seed is **M-M-specific**: only
:class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
carries a ``psi_half_left`` variable to seed.  The Legacy and
Bailey closures don't have one — their half-angle face flux
evaluation collapses to cell-centre values unconditionally.  Two
architectures were considered:

* **Option α (composition, shipped)** — the seed strategy lives as
  a :class:`~orpheus.sn.spatial.psi_half_angle_seed.PsiHalfAngleSeed`
  field on
  :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`.
  The abstraction stays local to the closure that consumes it.

* **Option B (sibling Protocol on SNMesh, rejected)** — the seed
  would be a separate Protocol attribute on
  :class:`~orpheus.sn.geometry.SNMesh`, applied by the matvec
  before calling the pole closure.  This would force every
  consumer (Legacy / BFF / M-M) to handle a Protocol that is a
  **no-op** for the non-M-M strategies, violating the
  single-responsibility principle and forcing unrelated tests to
  thread the Protocol through call signatures.

The
:class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonSweepContext`
dataclass bundles the four inputs the Carlson sweep needs that
are NOT in the standard
:class:`~orpheus.sn.spatial.pole_angular_closure.PoleAngularClosure`
``__call__`` signature (``sigma_t``, ``dr``, ``mu_quad``,
``weights``, ``bc_outer_value``), keeping the call-signature
expansion to a single new optional keyword — a minimal
blast-radius extension that Legacy and Bailey closures ignore by
documented Protocol contract.

Linear-operator preservation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both seed strategies — :class:`ZeroSeed` and
:class:`CarlsonInwardSweep` — are **linear in the input** ``psi_cells``
(verified by the ``is_linear: ClassVar[bool] = True`` trait, pinned
by foundation tests).  Linearity is the load-bearing property:
the apply matvec must be a linear operator, otherwise the
operator-algebra capabilities of
:class:`~orpheus.sn.operator.InvertibleOperator`
(apply, apply_transpose, dense matrix probing) break.  The
:class:`CarlsonInwardSweep` is linear because:

* The :math:`\phi_0` moment is a linear projection of input
  :math:`\psi` (Legendre integration is linear).
* :math:`\bar Q = \frac{1}{2} \Sigma_t \cdot \phi_0` is linear
  in :math:`\psi` (:math:`\Sigma_t` is constant).
* The recurrence Eqs. :eq:`hebert-3-434`–:eq:`hebert-3-435` is
  an affine function of :math:`(\bar Q, \bar\phi_{i+1/2})` with
  constant coefficients depending only on
  :math:`(\Sigma_t, \Delta r)`.
* The ``bc_outer_value`` is constructed in the matvec by applying
  the realised BC operator to the cell-centred outer-cell
  :math:`\psi`, then extracting the most-inward ordinate's value
  — both operations are linear in the input :math:`\psi`.

The foundation test
:func:`tests.sn.spatial.test_psi_half_angle_seed.test_carlson_inward_sweep_is_linear_in_input`
pins the linearity directly; the operator-level linearity gate
in :file:`tests/sn/test_streaming_operator.py` pins it transitively
at the matvec boundary (``rtol=1e-12`` — relaxed from the
pre-Phase-D ``rtol=1e-13`` to absorb ~10×ULP non-associativity
drift, justified by the three principled-relaxation criteria of
the ``vv-principles`` bit-identity-vs-principled-equivalence
framework).

The L = 0 isotropic-only limitation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current
:class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
evaluates only the :math:`\ell = 0` (isotropic) Legendre moment
when building the moment-folded source in
:eq:`hebert-3-432-source`.  This is **consistent with the apply
matvec's structure**: the
:class:`~orpheus.sn.operator.InvertibleOperator` apply matvec
carries only an isotropic collision term :math:`\Sigma_t \psi`;
anisotropic scattering (P\ :sub:`1`\ +) is composed externally via
a separate scattering operator, not included in :math:`L`.

.. warning::

   **The L = 0 isotropic-only assumption is load-bearing for the
   Phase D fix.**  If a future refactor moves scattering INTO
   :math:`L` (e.g., to enable a "monolithic" SN apply that
   includes within-group scattering), the Carlson seed becomes
   WRONG: the source at :math:`\mu = -1` (Eq. :eq:`hebert-3-432`)
   needs the full Legendre-moment sum

   .. math::

      \bar Q_i \;=\; \sum_\ell \frac{2\ell+1}{2}\,Q_\ell(r)\,(-1)^\ell,

   not just :math:`\Sigma_t \phi_0`.  This is a Mode-6
   convention-drift risk per the ``vv-principles`` skill (the
   definition-site assumption disagreeing with the usage-site
   intention).  A foundation test pinning the isotropic-only
   assumption (e.g., asserting the apply matvec does NOT couple
   to ``self_scattering``) would catch a future drift; in its
   absence, this WARNING block and the module docstring's
   matching admonition are the only safeguards.  Track the
   future-refactor case under a fresh GitHub issue when the
   monolithic apply work is scheduled.

.. _sn-phase-d-default-flips:

Default flips
~~~~~~~~~~~~~

Phase D ships **two default flips** that activate the full
canonical curvilinear closure path:

#. :attr:`SNMesh.pole_angular_closure
   <orpheus.sn.geometry.SNMesh.pole_angular_closure>` default
   flipped from
   :class:`~orpheus.sn.spatial.pole_angular_closure.LegacyTauSymmetricInterpolation`
   to
   :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`.
   :class:`MorelMontryAngularSweep`'s own constructor default for
   :attr:`psi_half_seed` is
   :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`,
   so the single :class:`SNMesh` flip activates the full Phase D
   fix (canonical M-M closure + canonical Carlson seed) without
   requiring downstream call sites to thread the new strategy
   explicitly.

#. :class:`~orpheus.sn.solver.SNSolver`'s ``inner_solver`` default
   flipped from ``"source_iteration"`` to ``"krylov"`` for
   **curvilinear geometries** (spherical, cylindrical); Cartesian
   stays at ``"source_iteration"``.  The rationale: the Phase D
   fix lives in the apply matvec, and the Krylov path is the one
   that uses
   :meth:`~orpheus.sn.operator.InvertibleOperator.apply`.  The
   sweep path (``"source_iteration"``) uses the spatial WDD
   recurrence and is unaffected by the Phase D fix — leaving its
   ERR-026-affected curvilinear behaviour in place would be wrong
   for the production default.

   .. note:: **Reverted (2026-06-12, Issue #195).**

      This curvilinear ``"krylov"`` default was undone by the ERR-058
      fix.  The premise — that the sweep path was ERR-026-affected
      while Krylov-on-apply was not — held only because the sweep and
      matvec were *distinct* discrete systems at Phase D time.  After
      the Depth-B/Wave-T unification they are ONE system, and the
      ERR-058 closure-seed fix makes that system O(h²)-consistent: SI
      :math:`\equiv` Krylov bit-identical, SI :math:`\sim 10^2\times`
      faster.  The curvilinear default returned to
      ``"source_iteration"``.  See
      :ref:`sn-err-058-closure-seed-closeout`.

.. _sn-phase-d-gate-1-1-empirical:

Empirical Gate 1.1 outcome (Phase D — full 12-cell crosstab)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Phase D acceptance gate is Gate 1.1 on **all three** pole
closures across both curvilinear geometries and both :math:`\Sigma_t`
values.  The parametrised test
:func:`tests.sn.test_phase_c_gates.test_apply_curvilinear_per_ordinate_flat_flux_residual`
produces the 12-cell crosstab:

.. list-table:: Gate 1.1 outcome under Phase D Carlson seed (2026-05-12)
   :header-rows: 1
   :widths: 18 30 26 26

   * - Geometry
     - Pole closure
     - :math:`\Sigma_t = 0`
     - :math:`\Sigma_t = 0.5`
   * - Sphere
     - ``LegacyTauSymmetricInterpolation``
     - PASS
     - PASS
   * - Sphere
     - ``BaileyFlatFluxRedist``
     - PASS
     - PASS
   * - Sphere
     - ``MorelMontryAngularSweep``
     - **XPASS**
     - **XPASS**
   * - Cylinder
     - ``LegacyTauSymmetricInterpolation``
     - PASS
     - PASS
   * - Cylinder
     - ``BaileyFlatFluxRedist``
     - PASS
     - PASS
   * - Cylinder
     - ``MorelMontryAngularSweep``
     - **XPASS**
     - **XPASS**

All 12 cells PASS or XPASS.  The 4 XPASS cells under M-M closure
are the ERR-026 markers — they now flip from FAIL to XPASS on
xfail-strict=False, unblocking the marker-removal commit that
Phase D Step 5 will execute (deferred per the closeout memo's
acceptance gate item 6).

This is the load-bearing **empirical evidence** for the
ERR-026 identity-and-rate scope closure.  The asymmetry between
the Phase C (cylinder PASS / sphere FAIL) and Phase D
(both PASS) crosstabs is the diagnostic mark of the Phase D
intervention: the sphere case required the Carlson seed because
its single-cascade structure has no telescoping; the cylinder case
already passed under Phase C because per-:math:`\mu`-level
:math:`\alpha`-dome telescoping absorbed the zero-seed
inconsistency.

.. _sn-phase-d-gate-1-5-capture-and-compare:

Gate 1.5 strengthening — capture-and-compare BC apply input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase C's Gate 1.5
(:ref:`bc-trace-contract-respected-by-matvec`) was a "round-trip"
check: invoke ``bc.realize().apply(...)`` independently and
compare against the matvec's observable output.  Phase D
strengthens this to a **capture-and-compare** check that pins the
exact value the matvec passes into the BC trace law:

#. Patch ``sn_mesh.bc["xmax"].apply`` (the outer radial face —
   a sphere's ``"outer"`` endpoint renders as ``"xmax"`` since
   C4 / #220, see :ref:`bc-face-name-carve`) to capture every input
   array passed to it during one matvec call.
#. Independently reconstruct the WDD-propagated outflow trace via
   a reference implementation (:func:`_outflow_at_boundary_for_sphere`).
#. Assert the captured BC apply input matches the reference to
   ``rtol=1e-14`` — exactly bit-equal up to FP non-associativity.

The strengthening matters because the Phase D matvec now calls
``bc["xmax"].apply`` **twice** per matvec:

#. **Phase D Carlson context call** — applied to cell-centred
   outer-cell :math:`\psi` to build ``bc_outer_value`` for the
   :class:`CarlsonSweepContext`.  See the BC companion section
   :ref:`bc-phase-d-two-bc-applies-per-matvec`.
#. **Phase C BC trace law call** — applied to the WDD-propagated
   outflow face value at the boundary edge, per the
   :ref:`affine-bc-form` contract.

The capture-and-compare test
:func:`tests.sn.test_phase_c_gates.test_bc_trace_contract_capture_and_compare_sphere`
(parametrised over ``vacuum`` and ``reflective``) **locates the
Phase C call by shape and content matching**: of the two captured
inputs, the one whose shape matches ``(N, ng)`` and whose values
match the independent reference is the Phase C trace law call.
Both vacuum and reflective parametrised cases pass; the test is
foundation-tagged because it pins a software invariant (the
matvec's two-application sequence) rather than a math claim.

.. _sn-phase-d-err-026-closure-narrative:

ERR-026 PARTIAL → PARTIAL (narrowed scope)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Retraction (2026-06-12, Issue #195).**

   The sub-claim table below classified the residual MMS magnitude as
   a benign "pre-asymptotic transient" that finer :math:`n_x` would
   clear (status OPEN, tracked at #195).  **That classification is
   wrong.**  The curvilinear isotropic MMS error did not shrink under
   refinement — it PLATEAUED mesh-independently (orders :math:`\to 0`),
   because the dominant defect was the *angular* closure seed (the
   Carlson proxy source), not a spatial-truncation constant.  The
   "rate :math:`[3.33, 2.46]` already correct, only the constant is
   large" reading was an artefact of the
   :math:`\alpha`-dome-telescoping blindness of the scalar residual.
   ERR-058 (#195) replaced both seeds; the isotropic MMS is now a
   clean :math:`\mathcal{O}(h^2)` ladder.  The table is preserved as
   bug-era evidence — its STATUS / Closed-by interpretation is
   superseded by :ref:`sn-err-058-closure-seed-closeout`.

Phase D **narrows** ERR-026's open scope.  The bug ERR-026
originally diagnosed — *"curvilinear sweep WDD angular closure
converges to wrong fixed-source solution"* — had three
sub-claims, each addressed by a different Wave:

.. list-table:: ERR-026 sub-claim closure tracking
   :header-rows: 1
   :widths: 35 35 30

   * - Sub-claim
     - Status
     - Closed by
   * - Operator identity:
       :math:`(L \cdot \psi_{\text{flat}})_{n,i,g} = \Sigma_t \cdot \psi_{n,i,g}`
       on per-ordinate flat-flux probe
     - **CLOSED**
     - Phase D Carlson seed (Gate 1.1 XPASS)
   * - Convergence rate:
       :math:`\mathcal{O}(h^2)` MMS rate at fixed :math:`N`
     - **CLOSED (rate)**
     - Phase D Carlson seed (empirical rate [3.33, 2.46] across
       refinements; both above the L1 acceptance floor of 1.9)
   * - Convergence magnitude: pre-asymptotic absolute MMS error
       below quadrature floor at practical ``nx`` (:math:`\le 160`)
     - **OPEN**
     - Tracked at `Issue #195
       <https://github.com/deOliveira-R/ORPHEUS/issues/195>`_;
       requires either finer ``nx`` or a higher-order spatial
       closure refinement to fully close

The convergence-rate evidence ``[3.33, 2.46]`` is the slope
sequence measured at successive refinement levels; both values are
above the L1 acceptance floor of 1.9 (second-order accuracy
demonstrated robustly), satisfying the rate sub-claim.  However,
the **absolute magnitude** at the largest tested ``nx`` (=160)
remains above the L1 tolerance that the test architect specified
for full closure on the pre-asymptotic regime.  This is **NOT** a
violation of the Phase D fix — the rate is correct, the asymptotic
regime is the right shape, but the **constant-coefficient** in
front of the :math:`\mathcal{O}(h^2)` term is larger than the
test's pre-asymptotic-magnitude budget at practical mesh
resolutions.

The pre-asymptotic regime is the consequence of the Carlson
sweep's L0-truncated source: at coarse :math:`nx` the Legendre
:math:`\phi_0` moment is computed from the cell-centred input
:math:`\psi` against the GL quadrature — an integration whose own
truncation contributes to the constant in
:math:`\bar Q = \frac{1}{2} \Sigma_t \phi_0`.  Refining
:math:`nx` reduces this contribution, but the rate at which it
reduces is set by the WDD spatial closure's own truncation order,
not by the Carlson sweep itself.  See Issue #195 for the candidate
follow-up paths (higher-order pole-face spatial closure, or a
:math:`\phi_0` recomputation that uses the M-M angular recurrence
output rather than the cell-centred input).

The 4 ``xfail-strict`` ERR-026 tripwires
(:file:`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`,
sphere + cylinder × isotropic + anisotropic ansatz) therefore
**stay xfail** through Phase D Step 3.  They will ``xpass`` under
the Phase D defaults (which is what triggers the deferred Step 5
marker-removal commit); the pre-asymptotic-magnitude regression
that prevents `strict=True` flipping is Issue #195's domain.  The
narrative for ``error_catalog.md`` therefore reads:

   ERR-026 status: **PARTIAL CLOSURE** (was PARTIAL through Phase
   C, narrowed scope through Phase D).  The structural bug (M-M
   recurrence hardcoded ψ\ :sub:`1/2,i` = 0 seed) is closed by the
   Phase D Carlson coupled-pole sweep; Gate 1.1 sphere MMS PASS
   confirms the operator identity and the second-order
   convergence rate is recovered.  The pre-asymptotic-magnitude
   open question (Issue #195) is what keeps the status at PARTIAL
   rather than CLOSED.

.. _sn-phase-d-files-touched:

Files touched by Phase D
~~~~~~~~~~~~~~~~~~~~~~~~

The full Phase D footprint (per the closeout memo at
``.claude/agent-memory/method-implementer/issue_168_phase_d_step3_closeout.md``):

**New modules**

* :mod:`orpheus.sn.spatial.psi_half_angle_seed` — Protocol family
  + ABC + 2 strategies (:class:`ZeroSeed` + :class:`CarlsonInwardSweep`)
  + :class:`CarlsonSweepContext` dataclass.
* :file:`tests/sn/spatial/test_psi_half_angle_seed.py` — 25
  foundation + L0 + L1 tests covering Protocol conformance,
  registry/self-registration, immutability, shape contract,
  bit-identity for :class:`ZeroSeed`, L0 algebraic identities
  (flat-:math:`\psi` at varying C, vacuum-BC nx=3 hand
  calculation, multi-region :math:`\Sigma_t` step), linearity, and
  L1 structural-independence (Carlson vs Zero on vacuum-BC probe).

**Modified files**

* :mod:`orpheus.sn.spatial.pole_angular_closure` —
  :class:`MorelMontryAngularSweep` gains a
  :attr:`psi_half_seed: PsiHalfAngleSeed` field;
  :func:`_mm_weighted_angular_recurrence_single_level` accepts an
  optional ``psi_half_seed`` array; Protocol signatures extended
  with an optional ``carlson_context`` kwarg (Legacy + Bailey
  ignore it).
* :mod:`orpheus.sn.spatial` ``__init__`` re-exports the new
  symbols.
* :mod:`orpheus.sn.operator` — spherical + cylindrical matvecs
  build the
  :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonSweepContext`
  before calling ``pole_angular_closure``.
* :mod:`orpheus.sn.geometry` — :class:`SNMesh` default flipped to
  :class:`MorelMontryAngularSweep`.
* :mod:`orpheus.sn.solver` — curvilinear default ``inner_solver``
  flipped to ``"krylov"``.
* :file:`tests/sn/test_phase_c_gates.py` — Gate 1.5 strengthened
  with capture-and-compare.
* :file:`tests/sn/test_streaming_operator.py` (post-D-K successor
  to the retired ``test_snstreamingoperator.py``) — 3 tests updated
  (one test docstring rewritten to pin the Phase D fix; two
  bit-identity tests threaded with ``sn_mesh.pole_angular_closure``;
  one linearity tolerance relaxed ``rtol=1e-13 → 1e-12``).

The agent-memory trail for Phase D session reproducibility:

* Literature memo:
  ``.claude/agent-memory/literature-researcher/phase_d_carlson_coupled_pole.md``
  — Hébert §3.9.4 derivation + flat-:math:`\psi` algebra +
  architecture-shape correction + open questions.
* Diagnostic memo:
  ``.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md``
  — empirical evidence + 4 plan corrections + structural-
  independence cross-check.
* Step 3 closeout:
  ``.claude/agent-memory/method-implementer/issue_168_phase_d_step3_closeout.md``
  — what shipped + 3 deviations + V&V evidence chain.
* Diagnostic script:
  :file:`tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py`
  — self-contained CLI probe reproducing the diagnostic table.


.. _sn-phase-f-carlson-sweep-path-backport:

Phase F Carlson seed sweep-path backport (Issue #168 Phase F)
-------------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   * Phase F (commit chain landed 2026-05-12 on
     ``refactor/sn-operator-algebra``, atop Phase E ``6708a4a``)
     backports the Phase D Carlson coupled-pole seed
     (:class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`,
     Hébert §3.9.4 Eqs. (3.432)–(3.435)) from the apply-matvec path
     (:func:`~orpheus.sn.operator.transport_operator_matvec_spherical`
     / ``_cylindrical``, fixed in Phase D Step 3) into the SI/sweep
     path
     (``_sweep_1d_spherical`` (the dissolved ``sweep.py``) and
     ``_sweep_1d_cylindrical``).
   * The bug is the **structural twin** of the Phase D defect: the
     SI loop in :file:`orpheus/sn/sweep.py` initialised
     ``psi_angle = np.zeros((nx, ng))`` at the spherical sweep
     entry (line 474, pre-Phase-F) and at the cylindrical per-level
     loop entry (line 634, pre-Phase-F) — the same hardcoded
     zero seed Phase D diagnosed as wrong-term-initialization on
     the apply-matvec twin
     (``orpheus/sn/spatial/pole_angular_closure.py:411``).
   * Phase F factors a **NEW free function**
     :func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`
     ``(Q_bar, sigma_t, dr, bc_outer_value) -> (ng, nx)`` that runs
     :eq:`hebert-3-434`–:eq:`hebert-3-435` driven by the SI
     within-group source ``Q_1d`` rather than by an apply-path
     :math:`\psi` Legendre fold.
     :class:`CarlsonInwardSweep.__call__` is refactored to delegate
     to the same helper after folding ``psi_level → Q̄ = 0.5 ·
     Σ_t · φ_0``.  **One helper, two consumers** — Cardinal Rule 2
     (architecture) enforced via reuse without duplication.
   * Empirical result on ``sphere_2g_3reg`` n=40
     (heterogeneous A|B|A reflective 2-group sphere):
     ``sf[0]/sf[1]`` ratio at the pole was **0.522** (DIVERGING
     to **0.473** under refinement to n=320); post-Phase-F it is
     **0.778** and STABLE under refinement (still 0.777 at n=320).
     The outer-cell reflective-face defect ``sf[-1]/sf[-2]`` was
     **0.887** → **0.997** (essentially CLOSED).
     :math:`\psi(r=0)` quasi-isotropy: ``cv(ψ@i=0)``
     **0.520** → **0.404**, ``max/min(ψ@i=0)`` **6.4×** →
     **1.16×** (Pomraning 1989 prediction substantially approached).
   * **What was logged as open** *(now CLOSED, #196)*: the residual
     O(h) per-cell WDD spatial-closure asymmetry between SI and Krylov
     paths was logged as **manifestation #7 of ERR-026**.  It is now
     **CLOSED** — ERR-058 (#195) showed the gap was a shared
     closure-seed defect, not a discretisation asymmetry; #196 verified
     SI :math:`\equiv` Krylov to the iteration floor on the
     heterogeneous eigenvalue path and added the permanent regression
     gate (see :ref:`sn-phase-f-residual-o-h-open` and
     :ref:`sn-issue-196-eigenvalue-equivalence`).  The Phase E
     flux-shape sentinel
     (:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`)
     **no longer xfails** — it runs as a plain L1 test, the
     structurally-independent Variant-α anchor.

The twin-path bug Phase D left open
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase D's fix lived entirely in the **apply-matvec path**.  The
Phase D Carlson seed is invoked by
:meth:`~orpheus.sn.operator.InvertibleOperator.apply` via the
:attr:`MorelMontryAngularSweep.psi_half_seed` composition; that
covers every Krylov-driven call.  But ORPHEUS's curvilinear
production default is **source iteration**, which dispatches
through :func:`~orpheus.sn.loss_representation.transport_sweep` rather than
through the apply matvec, and the two paths run **different code**
to seed the M-M half-angle recurrence:

.. list-table:: Apply vs SI/sweep dispatch divergence (pre-Phase-F)
   :header-rows: 1
   :widths: 24 38 38

   * - Path
     - Carlson seed site
     - Pre-Phase-F state
   * - Apply matvec (Krylov)
     - :func:`~orpheus.sn.spatial.pole_angular_closure._mm_weighted_angular_recurrence_single_level`
       :math:`\to`
       :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
       via :attr:`MorelMontryAngularSweep.psi_half_seed`
     - **CORRECT** — Phase D Carlson seed installed; Gate 1.1
       XPASS on the per-ordinate flat-flux residual probe
       (residual :math:`\le 10^{-15}`).
   * - SI/sweep (Source Iteration)
     - ``_sweep_1d_spherical`` (the dissolved ``sweep.py``) line 474
       (spherical) and ``_sweep_1d_cylindrical`` line 634
       (cylindrical per-:math:`\mu`-level loop)
     - **WRONG** — hardcoded ``psi_angle = np.zeros((nx, ng))``,
       the very same Phase B zero seed Phase D diagnosed and
       replaced on the apply-matvec twin.  The bug survived
       Phase D's regression suite untouched.

The cylindrical site has its own per-:math:`\mu`-level twin —
each level's azimuthal recurrence enters with the same hardcoded
zero.  Cylindrical Gate 1.1 passed empirically pre-Phase-D
because each level's :math:`\alpha`-dome ends at :math:`\alpha = 0`
by GL antisymmetry, **absorbing** the wrong seed through
level-edge cancellation — but Cardinal Rule 2 (architecture)
demands the structural fix on the sister path even when the
empirical signature is invisible there.

Phase F-Step-2 mesh-refinement evidence (sphere)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Step 2 numerics-investigator probe ran SN on
``sphere_2g_3reg`` (A|B|A reflective 2-group sphere, R=2.0 cm,
GL-8) at :math:`n_{\text{total}} \in \{40, 80, 160, 320\}` and
Variant α (composite-GL trajectory-resolvent reference) at
:math:`n_r \in \{24, 36, 48, 72, 96\}` matching effective
refinements.  The full table from
``.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md``:

.. list-table:: SN sphere pre-Phase-F mesh refinement (g=0 ratios)
   :header-rows: 1
   :widths: 10 18 18 18 18 18

   * - :math:`n_{\text{total}}`
     - :math:`k_{\text{eff}}`
     - :math:`\bigl|\text{sf}[0]/\text{sf}[1] - 1\bigr|`
       (pole)
     - :math:`\bigl|\text{sf}[N{-}1]/\text{sf}[N{-}2] - 1\bigr|`
       (outer)
     - log-log slope (pole)
     - log-log slope (outer)
   * - 40
     - 1.3578153066
     - 4.78e-01
     - 1.13e-01
     - —
     - —
   * - 80
     - 1.3576649296
     - 4.94e-01
     - 9.75e-02
     - −0.049 (**DIV**)
     - +0.21
   * - 160
     - 1.3576295736
     - 5.11e-01
     - 6.59e-02
     - −0.049 (**DIV**)
     - +0.57
   * - 320
     - 1.3576226569
     - 5.27e-01
     - 3.88e-02
     - −0.043 (**DIV**)
     - +0.76

A linear-in-:math:`h` extrapolation of the pole ratio gave
``ratio = 0.473 + 1.06·h`` — a **fixed structural asymptote
at 0.473**, not 1.  The outer ratio converged toward 1 at
:math:`\sim \mathcal{O}(h^{3/4})`, slower than the
:math:`\mathcal{O}(h^2)` DD interior, consistent with a
first-order BC-trace truncation that *vanishes* under
refinement.  Variant α at all five refinements gave inner /
outer ratios **monotonically → 1** (1.001949 → 1.000010 inner;
1.027508 → 1.001004 outer), confirming SN as the outlier and
ruling out the BC-interpretation alternative.

Per the
``vv-principles`` Step 2 decision matrix, the pole cell fires
**Branch 3 (DIVERGENT, high urgency)** and the outer cell fires
**Branch 1 (O(h^p), file follow-up)**.  This made the dispatch
to Step 3 (deep diagnostic) mandatory.

Phase F-Step-3 isolation: SI vs Krylov split
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Step 3 diagnostic ran the **same problem** through
:func:`~orpheus.sn.solver.solve_sn` with the only variable
changed being the ``inner_solver`` kwarg:

.. list-table:: SI vs Krylov on ``sphere_2g_3reg`` n=40 (pre-Phase-F)
   :header-rows: 1
   :widths: 26 18 18 19 19

   * - Inner solver
     - :math:`k_{\text{eff}}`
     - :math:`\text{sf}[0]/\text{sf}[1]`
     - :math:`\text{sf}[N{-}1]/\text{sf}[N{-}2]`
     - cv(ψ\@i=0)
   * - ``"source_iteration"`` (sweep)
     - 1.38069560
     - **0.5223**
     - 0.8871
     - **0.520**
   * - ``"krylov"`` (apply matvec)
     - 1.38464040
     - **1.0288**
     - 0.9745
     - 0.445

The Krylov path **eliminates the pole anomaly entirely** at
n=40, and Krylov's pole ratio converges to 1 cleanly under
refinement (1.029 at n=40 → 1.002 at n=80 → 1.0018 at n=160 —
:math:`\sim\mathcal{O}(h^2)` consistent with second-order DD).
Same materials, same quadrature, same mesh, same
:class:`MorelMontryAngularSweep` pole closure with the Phase D
Carlson seed installed on its ``psi_half_seed`` field — the
**only** difference is which inner-solver dispatch is used.
The Krylov path goes through
:meth:`InvertibleOperator.apply` (which consumes the Phase D
Carlson seed correctly); the SI path goes through
:func:`transport_sweep` (which carries the **legacy zero
seed** untouched by Phase D).

This split is the smoking gun that pins the bug to
:file:`orpheus/sn/sweep.py:474` (and :file:`orpheus/sn/sweep.py:634`
for the cylindrical per-level twin).  See
``.claude/agent-memory/numerics-investigator/phase_f_step3_diagnostic.md``
for the full empirical trail.

Source-driven Hébert (3.434)–(3.435) — the math
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The apply-matvec's Phase D Carlson seed consumes a
:math:`\psi`-current array shaped ``(ng, M, nx)`` and folds it
to :math:`\bar Q = \frac{1}{2} \Sigma_t \phi_0` via the
Legendre :math:`\ell = 0` projection (Eq. :eq:`hebert-3-432-source`
with :math:`P_0(-1) = 1`).  The SI/sweep path has **no such
current :math:`\psi` array at sweep start** — the entire point
of one SI iteration is to *produce* the updated angular flux.
What the SI loop **does** carry at sweep start is the
within-group source

.. math::
   :label: phase-f-q-1d-decomposition

   Q_{\text{1d}}(i, g) \;\equiv\;
   Q^{\text{scatt}}_{\text{within}}(i, g)
   \;+\; \frac{1}{k_{\text{eff}}}\,Q^{\text{fiss}}(i, g)
   \;+\; Q^{\text{ext}}(i, g),

i.e. the **isotropic** within-group source from the previous
power-iteration's scalar flux + fission moment + external
source.

On the fixed-point solution of the SI loop, the operator
identity :math:`L \cdot \psi = Q_{\text{1d}}` is satisfied
ordinate-by-ordinate, with :math:`L` carrying only
:math:`\Sigma_t \psi` for the current isotropic ORPHEUS scope.
The scalar-flux Legendre moment satisfies
:math:`\phi_0 = \sum_n \mathcal{W}_n \psi_n`.  Combining
gives, on the fixed point,

.. math::
   :label: phase-f-source-eq-sigt-phi0

   \Sigma_t(r) \cdot \phi_0(r) \;=\; Q_{\text{1d}}(r),

so the cell-averaged source at :math:`\mu = -1` (Eq.
:eq:`hebert-3-432-source` collapsed to :math:`L = 0`,
isotropic) admits two equivalent expressions:

.. math::
   :label: phase-f-q-bar-twin-forms

   \bar Q_i
   \;=\; \tfrac{1}{2}\,\Sigma_t(r_i) \cdot \phi_0(r_i)
   \quad\text{(apply path: builds }\phi_0\text{ from input }\psi\text{)}

   \bar Q_i
   \;=\; \tfrac{1}{2}\,Q_{\text{1d}}(r_i)
   \quad\text{(sweep path: takes }Q_{\text{1d}}\text{ directly).}

**The two are identical on the fixed point** by Eq.
:eq:`phase-f-source-eq-sigt-phi0`.  Off the fixed point they
differ by the SI residual :math:`r_k = Q_{\text{1d}} -
\Sigma_t \phi_0^{(k)}`, which vanishes as the SI loop
converges.  The sweep path's source-driven Carlson seed is
therefore the **canonically equivalent** invocation of the
same Hébert §3.9.4 math, packaged for a code path that has
:math:`Q_{\text{1d}}` available but not the per-ordinate
:math:`\psi`.

The factor :math:`\tfrac{1}{2}` is the Legendre fold weight
:math:`(2\ell + 1)/2` at :math:`\ell = 0` times
:math:`P_0(-1) = 1`.  For an :math:`L \ge 1` anisotropic
operator (not currently in scope for ORPHEUS's apply
matvec, but flagged in the
:class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
class docstring's L=0 WARNING block), additional terms
:math:`(2\ell + 1) Q_\ell \cdot (-1)^\ell / 2` for
:math:`\ell \ge 1` would enter — the source-driven helper
would need a moment vector ``Q_ell[ell, i, g]`` rather than
the present ``Q_bar[i, g]`` to recover the canonical
construction.

With :math:`\bar Q_i` from either formula, the inward DD
recurrence Eqs. :eq:`hebert-3-434`–:eq:`hebert-3-435`
proceed identically to the apply path:

.. math::
   :label: phase-f-carlson-seed-source-driven

   \bar\phi_i \;=\;
   \frac{\Delta r_i \cdot \tfrac{1}{2}\,Q_{\text{1d}}(r_i)
          \;+\; 2 \cdot \bar\phi_{i+1/2}}
        {\Delta r_i \cdot \Sigma_t(r_i) \;+\; 2},
   \qquad
   \bar\phi_{i-1/2} \;=\; 2 \cdot \bar\phi_i - \bar\phi_{i+1/2}

(sequential in cells from :math:`i = nx - 1` inward to
:math:`i = 0`, vectorised across groups).  The
``bc_outer_value`` at :math:`\bar\phi_{nx+1/2}` is the
outer-face angular flux at :math:`\mu = -1`, realised through
the BC operator on the persistent outflow buffer
``bc_outer``.

Equivalence on the converged eigenmode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Foundation test
:func:`tests.sn.spatial.test_sweep_vs_apply_consistency`
pins the source-vs-:math:`\psi` equivalence directly: for any
flat-:math:`\psi` field ``ψ_const`` with ``bc_outer_value =
ψ_const`` (reflective) and ``Q_1d = Σ_t · Σw · ψ_const`` (the
within-group source built by SI from
``φ_0 = Σw · ψ_const``), the two helpers return
**bit-identical seeds** (up to FP non-associativity).
Apply-path:
:class:`CarlsonInwardSweep`
``(psi_level=ψ_const·ones, ctx)`` produces ``Q̄ = 0.5 · Σ_t · Σw
· ψ_const``; sweep-path:
:func:`carlson_inward_sweep_from_source`
``(Q_bar=0.5·Q_1d, ...)`` produces the same ``Q̄`` — the
recurrence is identical, the bit-equal result is the
**single-invariant property** the test pins.

The architectural choice: one helper, two consumers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase F's structural choice is **factor the helper, delegate
from the strategy** — the Cardinal Rule 2 (architecture)
imperative.  The pre-Phase-F implementation had Eqs.
:eq:`hebert-3-434`–:eq:`hebert-3-435` open-coded inside
:meth:`CarlsonInwardSweep.__call__`.  Naive options for the
backport:

* **Option 1 (REJECTED) — duplicate the recurrence loop in
  the sweep path.**  Two copies of the inward DD recurrence,
  one driven by ``Q̄ = 0.5 · Σ_t · φ_0`` (apply path), one by
  ``Q̄ = 0.5 · Q_1d`` (sweep path).  Equivalent at the
  algorithmic level but a Cardinal-Rule-2 architecture
  violation: a future bug fix to one copy would need to be
  audited against the sister — exactly the failure mode that
  produced the Phase F bug in the first place.
* **Option 2 (REJECTED) — invoke**
  :class:`CarlsonInwardSweep`
  **directly from the sweep path with a synthesized**
  ``psi_level``
  **array.**  The strategy's ``__call__(psi_level,
  context)`` Protocol signature takes ``(ng, M, nx)`` — the
  sweep would have to allocate a flat-:math:`\psi` proxy of
  the right shape just to feed it through the Legendre fold
  that would extract ``φ_0`` from the proxy.  Mathematically
  equivalent but wasteful and obscures intent.
* **Option 3 (SHIPPED) — factor**
  :func:`carlson_inward_sweep_from_source`
  **as a free function that takes** ``Q̄``
  **directly, and have the strategy delegate.**

:meth:`CarlsonInwardSweep.__call__` now reads (in essence):

.. code-block:: python

   def __call__(self, psi_level, context):
       # ψ -> φ_0 -> Q̄ fold (apply-path-specific)
       phi_0 = np.einsum("gmi,m->gi", psi_level, context.weights)
       Q_bar = 0.5 * context.sigma_t * phi_0
       # Delegate to the source-driven recurrence
       return carlson_inward_sweep_from_source(
           Q_bar=Q_bar,
           sigma_t=context.sigma_t,
           dr=context.dr,
           bc_outer_value=context.bc_outer_value,
       )

The sweep path consumes the helper directly with ``Q_bar = 0.5
· Q_1d.T``.  **Single source of truth, two structurally
equivalent invocation points.**  A future bug fix to the
recurrence (e.g., an :math:`L \ge 1` anisotropic extension)
lands in **one** place; both consumers inherit it
automatically.

Why the cylindrical site needed the fix too
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cylindrical Gate 1.1 passes empirically pre-Phase-F (the
per-:math:`\mu`-level :math:`\alpha`-domes telescope back to
:math:`\alpha = 0` at the level boundary, absorbing the wrong
zero seed via cancellation — see the
:ref:`sn-phase-d-gate-1-1-empirical` discussion for the
sphere-vs-cylinder asymmetry).  Phase F nonetheless fixes
both sites for two reasons:

#. **Cardinal Rule 2 (architecture)**: structural alignment of
   the canonical math at both sites prevents a future
   refactor from introducing an asymmetric bug that only the
   sphere catches.  The sweep-path helper is the same code
   regardless of geometry; consuming it consistently from
   both geometries is the architecturally clean choice.
#. **Defense in depth against future stress probes**: if a
   future anisotropic cylindrical case were to break the
   level-edge telescoping (e.g., :math:`L \ge 1` MMS with
   non-trivial azimuthal modes), the wrong zero seed would
   reappear as a failure mode in cylinder.  Fixing it now is
   cheap insurance.

The cylindrical fix sits inside the per-:math:`\mu`-level
loop (lines 678–714 of :file:`orpheus/sn/sweep.py`).  The
helper is invoked **once per level** with the level-specific
``bc_outer_value`` extracted from the persistent outflow
buffer at the most-inward ordinate of the level.  The
linearity of the helper in ``Q_bar`` and ``bc_outer_value``
ensures the per-level invocations remain commutative with
the outer-loop level iteration.

Phase F empirical evidence (post-fix)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The post-Phase-F state recovers the canonical SN behaviour on
the smoking-gun case:

.. list-table:: ``sphere_2g_3reg`` n=40 — pre/post Phase F
   :header-rows: 1
   :widths: 36 22 22 20

   * - Diagnostic
     - Pre-Phase-F (SI)
     - Post-Phase-F (SI)
     - Krylov (reference)
   * - :math:`\text{sf}[0]/\text{sf}[1]`
       (pole ratio, target ~1)
     - **0.522**
     - **0.778**
     - 1.029
   * - :math:`\text{sf}[N{-}1]/\text{sf}[N{-}2]`
       (outer ratio, target ~1)
     - 0.887
     - **0.997**
     - 0.974
   * - :math:`\text{cv}(\psi@i=0)`
       (Pomraning isotropy, target ~0)
     - 0.520
     - **0.404**
     - 0.445
   * - :math:`\max/\min(\psi@i=0)`
       (target ~1)
     - **6.4×**
     - **1.16×**
     - 1.18×
   * - :math:`k_{\text{eff}}`
     - 1.38069560
     - 1.38069560
     - 1.38464040

The pole ratio jumps from a structural plateau at 0.473–0.522
(divergent under refinement) up to a stable ~0.778 plateau
that holds at n=320 — the **structural divergence is
gone**.  ``sf[0]/sf[1] = 0.778`` is not yet ``1`` because the
SI fixed point still differs from the Krylov fixed point by
the residual O(h) WDD spatial-closure asymmetry (see
:ref:`sn-phase-f-residual-o-h-open` below), but the
**diverging-vs-refinement** signature that made the Phase E
flux-shape sentinel xfail-strict is closed.

.. note:: **Retraction (2026-06-13, Issue #196).**

   The table below logs a residual SI :math:`\neq` Krylov
   :math:`\mathcal{O}(h)` gap (pole ratio 0.778 vs Krylov 1.029;
   :math:`\Delta k` 0.286 % at n=40 halving per mesh doubling) and
   reads it as a benign discretisation artefact of "two methods now
   solving the same equation".  **That interpretation is wrong.**  The
   methods did NOT yet solve the same equation at Phase F: the
   *shared* closure seeds were still O(1)-wrong on non-flat fields
   (ERR-058).  After ERR-058 (#195) fixed the seeds, the
   :math:`\mathcal{O}(h)` gap **collapsed to the iteration floor** —
   SI :math:`\equiv` Krylov to :math:`|\Delta k|\approx
   1.9\mathrm{e}{-11}` and L∞ flux-shape :math:`\approx
   2.4\mathrm{e}{-10}` on ``sphere_2g_3reg`` n=40 (from a bug-era
   3.9e-3 / ~30 %); the pole ratio reaches 1 to that floor.  The
   measured numbers stay below as bug-era evidence; the
   production-decision record and post-fix evidence are
   :ref:`sn-issue-196-eigenvalue-equivalence`.

Mesh-refinement convergence (SI vs Krylov, post-Phase-F):

.. list-table:: Post-Phase-F SI-vs-Krylov convergence on ``sphere_2g_3reg`` (bug-era — gap closed by ERR-058/#196)
   :header-rows: 1
   :widths: 10 22 18 22 18 14

   * - :math:`n`
     - :math:`k_{\text{eff}}` (SI)
     - :math:`\text{sf}[0]/\text{sf}[1]` (SI)
     - :math:`k_{\text{eff}}` (Kr)
     - :math:`\text{sf}[0]/\text{sf}[1]` (Kr)
     - :math:`\Delta k`
   * - 40
     - 1.38069560
     - 0.7776
     - 1.38464040
     - 1.0288
     - 0.286 %
   * - 80
     - 1.38075258
     - 0.7771
     - 1.38261730
     - 1.0125
     - 0.135 %
   * - 160
     - 1.38078077
     - 0.7771
     - 1.38167934
     - 1.0018
     - 0.065 %

*(Bug-era reading, retracted — see the note above.)* The
:math:`k_{\text{eff}}` gap between SI and Krylov drops
by a factor of 2 per mesh doubling — apparent clean
:math:`\mathcal{O}(h)` convergence to a shared limit.
Pre-Phase-F the SI sat on the wrong structural fixed point
(0.473–0.522 ratio asymptote diverging from 1) while Krylov
converged to ~1 — the two methods **solved different
equations**, and refinement made it worse for SI.  Phase F
removed the *divergent* signature but, as ERR-058 later
showed, left a *shared* O(1)-on-non-flat seed defect: the two
paths still did not solve the same equation, and the residual
:math:`\mathcal{O}(h)` gap above is the slow trace of that
shared defect, not a discretisation artefact.  ERR-058 (#195)
fixed the seeds and the gap collapsed to the iteration floor
(:ref:`sn-issue-196-eigenvalue-equivalence`); there is no
residual O(h) gap in production.

Files touched by Phase F
~~~~~~~~~~~~~~~~~~~~~~~~

**Modified production code**

* :mod:`orpheus.sn.spatial.psi_half_angle_seed` — NEW free
  function :func:`carlson_inward_sweep_from_source` (lines
  358–419 of the module); :meth:`CarlsonInwardSweep.__call__`
  refactored to delegate after folding ``psi_level → Q̄``;
  ``__all__`` extended.
* ``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) —
  :func:`_sweep_1d_spherical` line ≈ 472–530: replaces the
  legacy ``psi_angle = np.zeros((nx, ng))`` with the Phase F
  Carlson seed call (uses ``bc_outer_obj.apply(bc_outer)`` to
  derive ``bc_outer_value`` at the most-inward ordinate, mirror
  of the apply-path's Phase D logic);
  :func:`_sweep_1d_cylindrical` lines ≈ 678–714: per-level
  Carlson seed inside the :math:`\mu`-level loop, replaces
  the inline level-zero init.

**Tests added**

* :func:`tests.sn.test_phase_c_gates.test_sweep_curvilinear_per_ordinate_flat_flux_residual`
  — **Gate 1.6**, the dual of Gate 1.1 for the SI/sweep path.
  Parametrised over geometry (sphere × cylinder) and
  :math:`\Sigma_t \in \{0.5, 1.5\}`.  Pins
  apply-path-vs-sweep-path bit-identity on the helper output
  AND the flat-:math:`\psi` algebraic identity at :math:`\Sigma
  w = 2` (Hébert convention).  Carries
  ``@pytest.mark.verifies("dd-curvilinear-scalar")`` and
  ``@pytest.mark.catches("ERR-026")`` — see
  :ref:`sn-phase-f-test-wiring` for the proposed extension to
  the Phase F equation labels.
* :file:`tests/sn/spatial/test_sweep_vs_apply_consistency.py` —
  NEW file, **57 foundation tests** pinning:

  #. Apply-path vs sweep-path Carlson seed bit-equivalence on
     matching ``Q̄`` (the load-bearing structural invariant).
  #. Linearity of
     :func:`carlson_inward_sweep_from_source` in ``Q_bar`` and
     ``bc_outer_value`` independently (Protocol-shape contract
     preservation).
  #. SI-vs-Krylov :math:`k_{\text{eff}}` agreement on
     homogeneous reflective spheres (the degenerate case
     where the Phase F fix is invariant — same eigenvalue
     pre- and post-fix).

**Updated tests**

* :func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`
  — *(Phase-F action, since superseded.)* Phase F updated the
  ``xfail-strict`` reason string from *"UNRESOLVED structural
  discrepancy with hypothesised pole issue"* to *"Phase F closed
  gross divergence; residual O(h) drift awaits further work"*, on
  the expectation a future tightening would self-enforce removal.
  **The xfail was removed by ERR-058 (#195):** the canary now runs
  as a plain L1 test (the structurally-independent Variant-α anchor;
  see :ref:`sn-issue-196-eigenvalue-equivalence`).

**Snapshot regeneration**

* 6 curvilinear regression snapshots regenerated under the
  Phase F fix:

  * ``tests/sn/regression/snapshots/sphere_2g_homogeneous_dd_n20.npz``
  * ``tests/sn/regression/snapshots/sphere_2g_3reg_dd_n40.npz``
  * ``tests/sn/regression/snapshots/sphere_2g_p1_aniso_dd_n20.npz``
  * ``tests/sn/regression/snapshots/cyl_1g_homogeneous_LS4_dd_n20.npz``
  * ``tests/sn/regression/snapshots/cyl_1g_homogeneous_product_dd_n20.npz``
  * ``tests/sn/regression/snapshots/cyl_2g_3reg_LS4_dd_n40.npz``

  Bit-identity break is principled per the
  ``vv-principles`` *"Bit-identity vs principled-equivalence"*
  framework: the new seed is the canonical Hébert value
  (replaces the diagnosed wrong zero); the
  structurally-independent verification reference is
  Variant α (composite-GL trajectory-resolvent, accessed via
  Gate 4.2); the drift is algorithmic (intended) and
  well-defined.  All 5 Gate 4.2 snapshots still PASS at the
  Phase E tightened tolerances (sphere
  :math:`r_{\text{tol}} = 2 \times 10^{-2}`, cylinder
  :math:`3 \times 10^{-2}`).

.. _sn-phase-f-residual-o-h-open:

ERR-026 manifestation #7 — CLOSED by ERR-058 (#195), verified + pinned by #196
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: Status — manifestation #7 is CLOSED (Issue #196, 2026-06-13)
   :class: important

   **This was the LAST open manifestation of ERR-026; closing it
   formally retires the curvilinear-SN wrong-fixed-point family.**  The
   "residual :math:`\mathcal{O}(h)` SI-vs-Krylov gap" reading below —
   including the SI :math:`\neq` Krylov tables above (pole ratio 0.778
   vs Krylov 1.029; :math:`\Delta k` converging :math:`\mathcal{O}(h)`)
   and Options (a)/(b)/(c) — was the *two-distinct-systems* picture and
   is **bug-era history**.  The gap was NOT a discretisation artefact;
   it was the shared closure-seed defect (ERR-058), manifest
   differently on the two then-distinct paths.

   * **ERR-048** (Phase G Step 2, 2026-05-13) closed only the **L0
     flat-field** twin-agreement: it patched the SI sweep to MATCH the
     apply-matvec conventions on the homogeneous streaming-equilibrium
     gauntlet (pole-face WDD IC mirror + Carlson seed normalisation).
     The **L1 heterogeneous eigenvalue** :math:`\mathcal{O}(h)`
     asymmetry that manifestation #7 names **PERSISTED** — which is
     exactly why #196 stayed OPEN — because the shared closure seeds
     were still *exact-on-flat / O(1)-wrong-on-non-flat* (the ERR-058
     defect).
   * **ERR-058** (Issue #195, 2026-06-12) was the TERMINAL fix: it
     replaced the shared closure seeds with correct ones — the
     coupled-pole spatial seed :math:`\psi(0,+\mu)=\psi(0,-\mu)` and the
     :class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`
     half-angle seed — so BOTH inner solvers operate on the SAME correct
     discrete operator.
   * **#196** (2026-06-13) VERIFIED the eigenvalue-path equivalence and
     added the permanent regression gate.  See
     :ref:`sn-issue-196-eigenvalue-equivalence` for the measured
     evidence (sphere :math:`|\Delta k|=4.68\mathrm{e}{-12}`, cylinder
     :math:`1.91\mathrm{e}{-11}` on the bug-era snapshot cases) and the
     gate description.

   Option (c) (keep SI, accept an O(h) gap) is moot — there is no gap;
   Option (b) (flip to Krylov) is the opposite of what landed (SI is
   restored as the faster default).  The full production-decision
   record is :ref:`sn-issue-196-eigenvalue-equivalence`.

The bug-era reading (preserved as history) ran: Phase F closed the
**structural** pole defect (the divergent ratio at the pole cell on
heterogeneous MR) and the **outer-cell** defect (sf[-1]/sf[-2]
essentially reaches 1); what was thought to remain was a milder
**convergence-rate** gap between SI and Krylov on heterogeneous MR
snapshots (at n=40 per-cell shape differing by ~5 %, apparently
converging :math:`\mathcal{O}(h)` toward zero under refinement),
logged in ``error_catalog.md`` as **ERR-026 manifestation #7**:

   *"SI-vs-Krylov per-cell agreement (residual O(h) WDD
   asymmetry) — OPEN, new follow-up after Phase F."*

That row now reads **CLOSED by ERR-058 (#195), verified + pinned by
#196**.  The Phase E flux-shape sentinel
:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`
**no longer xfails** — it runs as a plain L1 test (the
structurally-independent Variant-α anchor; see
:ref:`sn-issue-196-eigenvalue-equivalence`).  The two viable
closures that were tracked as Phase F-extensions are recorded here
**only as bug-era history** — neither was taken, because both
presupposed the shared fixed point was correct and only the
arithmetic differed, which the terminal diagnosis refuted:

* **Option (a) — Sweep WDD-closure refinement** *(bug-era,
  not taken).*  Investigate the per-cell WDD recurrence
  :math:`\psi_{n+1/2} = (\psi_n - (1-\tau)\psi_{n-1/2})/\tau`
  in ``_sweep_1d_spherical`` to identify the residual
  numerical asymmetry vs the apply matvec's symmetric closure
  :math:`\psi_{n\pm 1/2} = \tau \psi_{\text{next}} +
  (1-\tau) \psi_{\text{this}}`.  This presumed the seed was
  correct and only the spatial closure differed — false: the
  seed itself was the defect.
* **Option (b) — Flip curvilinear ``inner_solver`` default to
  Krylov** *(bug-era, not taken — in fact reverted).*
  :func:`solve_sn` for spherical / cylindrical would route
  through the Krylov inner (which carried the Phase D Carlson
  seed and produced the cleanly-converging fixed point).  This
  was the Phase-D-era flip; ERR-058 made the SI sweep correct, so
  the default was **reverted to ``source_iteration``** (SI is
  :math:`\sim 10^2\times` faster and now equivalent).

Phase F shipped **option (c)** at the time (keep SI default, achieve
structural alignment of the seed math, accept a residual O(h) gap)
on the reasoning that the methods "now solve the same equation".
The terminal diagnosis (ERR-058) showed they did **not** yet solve
the same equation — the *shared* fixed point was itself wrong on
non-flat fields.  Once the seeds were fixed there is no residual gap
to accept: SI and Krylov agree to the iteration floor at the
eigenvalue level (see :ref:`sn-issue-196-eigenvalue-equivalence`),
and bit-identically at the fixed-source level.

The anti-pattern Phase F surfaced
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Twin-path fix incompleteness** is a Mode-3
(missing-factor / wrong-term-initialization) anti-pattern.
The Phase D fix was scoped to the apply-matvec path because
Gate 1.1 runs through the apply-matvec; the SI/sweep path's
zero seed was untouched.  The bug survived Phase D's entire
regression suite untouched because:

* Phase D's Gate 1.1 MMS xfail-strict marker is on the
  apply-matvec path; the SI/sweep path didn't run that probe.
* The 6 curvilinear regression snapshots were SI-generated
  under the wrong seed; the snapshots **encoded the bug
  bit-identically** and "passed" by tautology.
* Homogeneous degenerate cases (1G, flat-flux reflective) gave
  k = νΣ_f / Σ_a independent of the flux shape, masking the
  structural divergence on the eigenvalue side.
* The heterogeneous-MR case was marked ``xfail`` for **flux
  shape** (Phase E), not for **eigenvalue**, so the
  shape-sentinel signal was deliberately not enforced.

The lesson, **proposed for addition to**
``vv-principles/SKILL.md`` *§ Anti-patterns* (per the Phase F
closeout memo §"Lessons (proposed for skill catalogue)"):

   *Whenever a fix is applied to one of two structurally-
   mirrored production paths (apply-matvec vs SI/sweep,
   prepass vs postpass, etc.), MUST audit the OTHER path for
   the same defect.  Mode 3 wrong-term-initialization
   defects often appear in pairs; fixing one path without
   auditing its sister is a Cardinal Rule 2 (architecture)
   violation that ERR-026 instantiated twice.*

.. _sn-phase-f-test-wiring:

Test wiring proposal — Phase F equation labels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase F declares three new equation labels:
:eq:`phase-f-q-1d-decomposition`,
:eq:`phase-f-source-eq-sigt-phi0`,
:eq:`phase-f-q-bar-twin-forms`, and
:eq:`phase-f-carlson-seed-source-driven`.  The
:eq:`phase-f-carlson-seed-source-driven` label is the
canonical Hébert (3.434)–(3.435) recurrence in
source-driven form — semantically the **same recurrence**
as :eq:`hebert-3-434` and :eq:`hebert-3-435` but with the
sweep-path source substitution made explicit.

The new Gate 1.6 test
:func:`tests.sn.test_phase_c_gates.test_sweep_curvilinear_per_ordinate_flat_flux_residual`
already carries
``@pytest.mark.verifies("dd-curvilinear-scalar")`` and
``@pytest.mark.catches("ERR-026")``.  Per the project's
V&V harness wiring, the test SHOULD additionally declare:

* ``@pytest.mark.verifies("phase-f-carlson-seed-source-driven")``
  — pins the source-driven recurrence on the bit-identity
  helper-vs-strategy probe.
* ``@pytest.mark.verifies("phase-f-q-bar-twin-forms")``
  — pins the apply-vs-sweep source-equivalence identity (the
  load-bearing structural invariant of Phase F).

The other two labels
(:eq:`phase-f-q-1d-decomposition` and
:eq:`phase-f-source-eq-sigt-phi0`) document the
*decomposition* of the SI source and the *fixed-point
identity* :math:`\Sigma_t \phi_0 = Q_{\text{1d}}`; both are
verified transitively by the existing SI convergence
infrastructure (the SI inner-tolerance is the gate that
enforces the fixed-point identity to machine precision).
The proposed wiring is tracked as a follow-up to the V&V
audit harness (see Issue #194 for the sister case of
``hebert-3-43X`` labels — same pattern, same fix).

Pointers
~~~~~~~~

* **Phase F plan**:
  ``.claude/plans/issue_168_phase_f_curvilinear_boundary_eigenvector.md``
  — context, hypothesis, three-step structure, sub-agent
  dispatch chain.
* **Step 2 numerics memo**:
  ``.claude/agent-memory/numerics-investigator/phase_f_step2_mesh_refinement.md``
  — mesh-refinement convergence study, SN-vs-Variant-α
  outlier identification, Step 2 branch-3 decision.
* **Step 3 diagnostic memo**:
  ``.claude/agent-memory/numerics-investigator/phase_f_step3_diagnostic.md``
  — fix-site identification (the smoking gun), SI-vs-Krylov
  isolation, Option-A-vs-B implementation analysis.
* **Phase F closeout memo**:
  ``.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md``
  — what shipped, the empirical evidence tables, files
  touched, residual-open items.
* **ERR-026 catalogue narrative**:
  ``.claude/skills/vv-principles/error_catalog.md`` (§ ERR-026,
  *"What Wave H Phase F added"*) — manifestation table update
  #6 CLOSED, #7 (new) OPEN.
* **Sister section on the BC apply call sequence**:
  :ref:`bc-phase-f-three-bc-applies-per-sweep-iteration` in
  :doc:`boundary_conditions` — extends the Phase D
  two-BC-applies-per-matvec narrative to cover the SI sweep's
  Phase F invocation.


.. _sn-err-058-closure-seed-closeout:

ERR-058 — the curvilinear closure-seed fix (Issue #195 CLOSED)
--------------------------------------------------------------

.. admonition:: Status banner
   :class: important

   **Issue #195 CLOSED 2026-06-12.**  ERR-058 closes the curvilinear
   *wrong-fixed-point* family — the open loop the Phase A–F narrative
   above tracked under the name "ERR-026 PARTIAL CLOSURE".  Two
   independent closure SEEDS in the curvilinear within-group operator
   were wrong on every non-flat field; both are now replaced.  In
   production:

   * The **half-angle thread seed** is
     :class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`
     (the new
     :class:`~orpheus.sn.spatial.pole_angular_closure.MorelMontryAngularSweep`
     ``psi_half_seed`` default).  It replaces
     :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`,
     whose proxy source :math:`\bar Q = \Sigma_t\phi_0/\!\sum w` was the
     dominant defect.
   * The **spatial pole-face seed** of the outward (:math:`\mu>0`)
     sweep is the *mirror inward sweep's pole-face outflow* — the
     Carlson coupled-pole continuity :math:`\psi(0,+\mu)=\psi(0,-\mu)`
     — replacing the historical innermost-cell-centre read
     :math:`\psi(\Delta r/2)`.
   * The curvilinear inner default returned from ``"krylov"`` to
     ``"source_iteration"`` (both
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` and the
     eigenvalue :func:`~orpheus.sn.solver.solve_sn`): post-unification
     the sweep and matvec are ONE discrete system, so SI
     :math:`\equiv` Krylov **bit-identical for fixed-source** and to
     the **iteration floor for eigenvalue** (the eigenvalue solve wraps
     the inner in power iteration — see
     :ref:`sn-issue-196-bit-identical-vs-floor`).  SI is
     :math:`\sim 10^2\times` faster than GMRES (no restart, ERR-053).

   :class:`CarlsonInwardSweep` is **retained** (not deleted) as the
   registered host of the canonical Hébert §3.9.4 recurrence
   (:func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`),
   reachable only by explicit opt-in.  Its class docstring carries a
   ``.. warning::`` block recording the proxy-source caveat by design,
   so a future session cannot re-activate it as a default unaware of
   the falsification.

   The **anisotropic** curvilinear MMS gates improved :math:`\sim 50\times`
   and are now limited by a *fixed-quadrature angular floor* of the
   half-angle thread interpolation — a test-design retune tracked at
   `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_,
   **not** a residual instance of this wrong-fixed-point class.

Motivation preserved — what the Phase A–F loop was chasing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Phase A–F sections above are preserved verbatim as the
*investigation history*.  Their reasoning was sound at the time and is
pedagogically load-bearing — a future reader asking "why did anyone
try a Carlson inward sweep, an apply-vs-sweep twin audit, a
Krylov-default flip?" must find the answer there.  This subsection only
flips the tenses on the *terminal* claims those sections reached:

* Phase D **was expected to** close ERR-026 once the apply-matvec
  half-angle seed was made non-zero; it closed the per-ordinate
  *flat-flux identity* but left the assembled operator wrong on
  non-flat profiles (the Carlson proxy source).  The "PARTIAL CLOSURE /
  pre-asymptotic transient" framing it shipped **was** the best
  available reading of the evidence then; it is now superseded.
* Phase F **was expected to** close the SI-vs-Krylov gap by backporting
  the same Carlson seed to the sweep.  It did make the two paths share
  the *seed strategy*, but both still drove it with the wrong proxy
  source, so the residual "manifestation #7 O(h) gap" it logged was a
  symptom of the shared defect, not a discretisation artefact.  After
  the Depth-B/Wave-T matvec unification, the sweep and matvec became
  ONE discrete system; the gap vanished by construction once the seed
  was fixed.

The premise the *issue itself* carried — a benign "pre-asymptotic
transient" that finer meshes would clear — was **empirically refuted**:
on ``main`` the isotropic curvilinear MMS error PLATEAUS
mesh-independently (sphere :math:`\sim 0.0413`, cylinder
:math:`\sim 0.0494`, :math:`n_x` 20 :math:`\to` 640, orders
:math:`\to 0`), with SI :math:`\equiv` Krylov bit-identical.  No
refinement helps a plateau.

The two manifestations — one class, both flat-field-exact
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ERR-058 is a **closure-seed inconsistency**: two discrete seeds inside
the curvilinear within-group operator were constructed so that they are
*exact* on spatially/angularly flat fields and *O(1)-wrong* on every
other field.  Because a discrete closure seed is part of the operator,
each seed had to be verified per-ordinate on a NON-flat field — and was
not.

.. _sn-err-058-manifestation-a:

Manifestation (a) — the spatial pole-face seed
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The outward (:math:`\mu>0`) radial sweep needs an inflow value at the
pole face :math:`r=0`.  The historical matvec
(now the fused
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`) and
the sweep both read the innermost CELL-CENTRE flux :math:`\psi(\Delta r/2)`
as if it were the pole-FACE value — a half-cell offset.  On a flat
radial profile :math:`\psi(\Delta r/2)=\psi(0)`, so the read is exact;
on the manufactured :math:`A(r)=\sin(\pi r/R)` ansatz it is
:math:`\mathcal{O}(h)`-wrong.  The DD face chain propagates that seed
error as an *undamped odd–even alternation*, and the area-weighted
streaming amplifies it by :math:`\sim A/V \sim 1/r` near the pole.

**The fix — Carlson coupled-pole continuity.**  At :math:`r=0` the
outward characteristic is the *continuation* of the inward one: a
neutron travelling inward along :math:`-\mu` that reaches the centre
emerges travelling outward along :math:`+\mu`, so

.. math::
   :label: sn-err-058-coupled-pole-continuity

   \psi(0,\,+\mu) \;=\; \psi(0,\,-\mu).

.. (vv-status rationale) Representational identity: the r=0
   pole-continuity boundary condition coupling the mirror ordinates.
   Not a solver claim (no eigenvalue / flux value). The verifiable
   content is the per-ordinate operator-admission gate
   (test_curvilinear_operator_admits_mms, catches ERR-058) + the
   strategy-owned seed-adjoint bit-identity (test_g_adjoint_reciprocity).
.. vv-status: sn-err-058-coupled-pole-continuity documented

The :math:`-1` (inward) sweep is therefore run FIRST.  Its pole-face
outflow, read at the *mirror* ordinate
``quad.reflection_index("x")``, seeds the :math:`+1` (outward) sweep.
This is **data** — propagated from the outer boundary, lower-triangular
in cell-visit order — not a self-reference.  It is the
"inward-determines-outward" pole condition deferred at Phase C
(`Issue #192 <https://github.com/deOliveira-R/ORPHEUS/issues/192>`_),
now landed.  The seed is exact on flat :math:`\psi` (so every flat-flux
gate is untouched), lower-triangular (so the operator stays
forward-substitutable), and the matvec and sweep capture/consume it
identically (so the pair stays ONE discrete system).

The **adjoint** routes the :math:`+1` seed cotangent into the
:math:`-1` reversal's initial outflow cotangent at the mirrored
ordinates (see
:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.loss_action_transpose`,
pinned by the dense-probe transpose oracle
``derivations/diagnostics/diag_p42_adjoint_oracle.py``).

.. _sn-coupled-pole-mu-level-invariant:

The μ-level-preservation invariant the mirror seed relies on
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The coupled-pole continuity :eq:`sn-err-058-coupled-pole-continuity`
seeds the outward (:math:`+\mu`) pole face from the inward (:math:`-\mu`)
sweep's pole outflow read **at the mirror ordinate** — concretely,
``pole_face_seed = outflow_at_inner.T[quad.reflection_index("x")]`` in
the fused matvec
(:meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk` and
its adjoint partner) and ``psi_in = pole_outflow[mirror[global_n]]`` in
the SI sweep twin (:meth:`~orpheus.sn.loss_representation.transport_sweep`).
That single index — ``reflection_index("x")[n]`` — is what makes the
seed correct, and it carries a load-bearing assumption that is invisible
in the code but essential to the physics.

**The invariant.**  For the mirror seed to realise
:math:`\psi(0,+\mu_r)=\psi(0,-\mu_r)`, the partner
``reflection_index("x")[n]`` MUST be the *intra-level sign-flip partner*
of ordinate :math:`n`: the ordinate in the **same** :math:`\mu`-level
(same axial cosine :math:`\mu_z` — the level index) with the radial
cosine :math:`\mu_x` negated and :math:`\mu_y,\mu_z` held.

.. math::
   :label: sn-coupled-pole-mu-level-invariant-eq

   m \;=\; \mathrm{reflection\_index}(\text{"x"})[n]
   \;\Longrightarrow\;
   \mu_x[m] = -\,\mu_x[n],\quad
   \mu_y[m] = \mu_y[n],\quad
   \mu_z[m] = \mu_z[n].

.. (vv-status rationale) Structural / representational identity: the
   defining property the x-mirror partner MUST satisfy for the
   coupled-pole seed to be physically correct (intra-level μ_x
   sign-flip). Not a solver claim (no eigenvalue / flux value). The
   verifiable content is the foundation gate
   test_x_reflection_is_intra_level_signflip_partner (asserts the three
   equalities over gauss_legendre/level_symmetric/product cubatures) +
   its involution sibling; documented-only.
.. vv-status: sn-coupled-pole-mu-level-invariant-eq documented

**Why the physics demands it.**  The pole continuity is a statement at a
*fixed axial direction*: a neutron travelling inward along :math:`-\mu_x`
at axial cosine :math:`\mu_z` reaches the centre and emerges travelling
outward along :math:`+\mu_x` at the **same** :math:`\mu_z`.  The axial
component does not turn at the pole — only the radial one reverses.  So
the reflected partner must stay in the same :math:`\mu`-level; a
cross-level partner would couple two different axial directions and seed
the outward sweep with a value from the *wrong* characteristic.

**Why it holds by construction today.**  Two facts conspire.  First,
``reflection_index("x")`` resolves to
``reflection_partners[0] = _find_reflections(-mu_x, mu_y, mu_z, ...)``
(:func:`orpheus.numerics.quadrature.directional._compute_sphere_reflection_partners`),
which nearest-neighbour-matches each node against its image with
**only** :math:`\mu_x` negated — :math:`\mu_y,\mu_z` are passed through
unchanged.  Second, the cylinder/sphere level is grouped on the
**axial** cosine: the level factories key ``level_indices`` on
:math:`|\mu_z|` (sphere / level-symmetric — ``rules_sphere.py``) or hold
:math:`\mu_z=\mu_{\rm GL}` fixed per level (product — ``rules_product.py``),
never on :math:`\mu_x`.  Because the x-mirror holds :math:`\mu_z` and the
level is indexed by :math:`\mu_z`, the x-mirror provably maps an ordinate
to a partner in its own level.  The two facts are *independent* code
sites, so the invariant is an emergent property of their agreement — not
something either site enforces alone.

.. warning::

   **This is a silent-corruption surface — a Mode-7 blind spot at the
   operator-internals level.**  If ``reflection_index("x")`` ever
   returned a **cross-level** partner — a future cubature whose
   reflection table is built differently, or a refactor of
   :func:`~orpheus.numerics.quadrature.directional._compute_sphere_reflection_partners`
   that no longer holds :math:`\mu_z` — then
   ``pole_outflow[mirror[n]]`` would read a *different axial direction's*
   pole value, and the break would be **completely silent under the
   existing solver suite**.  The reason is the same blindness that hid
   ERR-058 itself: on a spatially/angularly **flat** :math:`\psi` field
   the mirror partner's pole value equals the ordinate's own value, so
   the seed is exact regardless of *which* ordinate it reads.  Every
   flat-flux gate, every streaming-equilibrium L0, every reflective
   :math:`k_\infty` check would stay green while the operator quietly
   coupled the wrong characteristics on any non-flat field (``vv-principles``
   Mode 7 — the ansatz-simplification blindness — operating on the
   operator's own internals, exactly the ERR-058 class).  A scalar /
   particle-balance residual cannot see it either, because the
   :math:`\alpha`-dome telescoping (above) sums away per-ordinate seed
   errors.

**Why it now has its own foundation gate.**  Because the solver tests are
structurally blind to a cross-level regression, the invariant is pinned
*directly* — not through any flux or eigenvalue, but as a property of the
quadrature's reflection table itself — by the foundation test
:func:`tests.sn.sweep.curvilinear.test_coupled_pole_mu_level_invariant.test_x_reflection_is_intra_level_signflip_partner`.
It asserts all three equalities of
:eq:`sn-coupled-pole-mu-level-invariant-eq` (intra-level membership,
:math:`\mu_x` sign-flip, :math:`\mu_y,\mu_z` held) over the
``gauss_legendre`` / ``level_symmetric`` / ``product`` cubatures the
curvilinear sweep actually uses; the sibling
:func:`~tests.sn.sweep.curvilinear.test_coupled_pole_mu_level_invariant.test_x_reflection_is_an_involution`
pins ``mirror ∘ mirror = id`` (the partner relation is symmetric — a
necessary corollary of the sign-flip).  Both are
``@pytest.mark.foundation``; the first carries
``@pytest.mark.verifies("sn-err-058-coupled-pole-continuity")``, tying
the table-level invariant to the continuity equation it underwrites.
This gate is the regression tripwire that turns the silent-corruption
surface into a loud one: a cross-level reflection table fails it
immediately, *before* any solver ever runs.

.. note::

   **Re-scope of Issue #193.**  This invariant is what
   `Issue #193 <https://github.com/deOliveira-R/ORPHEUS/issues/193>`_
   now pins.  The issue *originally* targeted a different "level-locality"
   concern — that the cylindrical matvec's
   ``bc_outer.apply``-once-then-per-level-extract pattern stayed correct.
   That concern **dissolved**: Wave O O.4a.2 removed the ``bc_outer.apply``
   keystone from the matvec entirely (the reflective coupling :math:`B`
   moved *outside* the bare sweep as a first-class sibling — see the
   boundary-condition extraction record at :ref:`bc-extraction`), and the
   surviving SI-sweep seed reads the **raw** inflow trace with no
   ``apply`` at all, so the
   apply/restrict commutativity the original test would have guarded is
   now vacuous.  The genuinely load-bearing :math:`\mu`-level invariant
   *moved* to the coupled-pole seed mirror documented here, and that is
   what #193 was re-scoped to gate.

.. _sn-err-058-manifestation-b:

Manifestation (b) — the angular half-angle thread seed (the dominant defect)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Morel–Montry angular recurrence (Hébert §3.9.4 Eqs. 3.437/3.439)
threads the half-angle face fluxes
:math:`\psi_{m\pm 1/2,i}` across a :math:`\mu`-level and needs a starting
seed :math:`\psi_{1/2,i}` at the level's most-inward angular edge.  The
Phase D / Phase F :class:`CarlsonInwardSweep` solved the canonical
*sweep-side* starting-direction ODE (Hébert Eqs. 3.432–3.435) for that
seed, but drove it with the **proxy source**

.. math::
   :label: sn-err-058-proxy-source

   \bar Q_i \;=\; \frac{\Sigma_{t,i}\,\phi_{0,i}}{\sum_n w_n},
   \qquad \phi_{0,i} = \sum_n w_n\,\psi_{n,i},

.. (vv-status rationale) Literature-transcribed definition of the
   falsified proxy source (the CarlsonInwardSweep half-angle seed).
   Recorded as the diagnosed defect, not a solver claim. Its falsity
   is what the per-ordinate operator-admission gate (catches ERR-058)
   detects; documented-only.
.. vv-status: sn-err-058-proxy-source documented

which equals the true within-group source ONLY at the flat-flux
equilibrium :math:`\Sigma_t\phi_0 = \bar Q`.  On any non-equilibrium
field — every MMS reference, every vacuum or heterogeneous problem —
the seed solves the *wrong* starting-direction ODE.  The measured
consequence on the isotropic MMS input
(:math:`\psi_n = A(r)/W`, scalar value :math:`0.5`): the seed returns
:math:`\bar\phi = 0.5777` where the correct angle-flat value is
:math:`0.5000`, and the per-ordinate redistribution residual reaches
:math:`\pm 55` at the pole, :math:`\pm 13` in the bulk, against a
continuous streaming of :math:`\pm 0.31`.  **This was the dominant
defect.**

**The fix —**
:class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`.
For the *operator* (matvec) to be consistent, the seed must approximate
the *input field's* own value at the level edge :math:`\mu_{\rm start}`
— a pure angular-extrapolation problem with NO dependence on
:math:`\Sigma_t`, the source, or the boundary trace.  The new strategy
extrapolates linearly in :math:`\mu` through the level's two most-inward
distinct-:math:`\mu` ordinates :math:`(m_0, m_1)`:

.. math::
   :label: sn-err-058-edge-extrapolation

   \psi_{1/2,i} \;=\; (1-t)\,\psi_{m_0,i} + t\,\psi_{m_1,i},
   \qquad
   t \;=\; \frac{\mu_{\rm start} - \mu_{m_0}}{\mu_{m_1} - \mu_{m_0}}.

.. (vv-status rationale) Representational identity: the
   operator-consistent half-angle thread seed (AngularEdgeExtrapolation,
   the new psi_half_seed default) as a fixed linear map. Not a solver
   claim. The verifiable content is the per-ordinate operator-admission
   gate (catches ERR-058), the isotropic MMS L1 ladders, and the
   strategy-owned seed-adjoint bit-identity; documented-only.
.. vv-status: sn-err-058-edge-extrapolation documented

The starting-direction edge :math:`\mu_{\rm start}` (sphere
:math:`-1`; cylinder :math:`-\sqrt{1-\xi_p^2}`, the level's most-inward
azimuthal edge) is single-sourced from the SAME
:math:`\alpha`/:math:`\tau` construction site as the
:math:`\alpha`-dome (``orpheus.geometry.reduced_operator``) and threaded
to the strategy via the REQUIRED
:attr:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonSweepContext.mu_start`
field — no default, so a forgotten cylinder site cannot silently fall
back to the sphere value.

**Exactness ladder.**  The extrapolation is

* **exact on angle-flat fields**, because the barycentric weights sum
  to one: :math:`(1-t)+t=1`.  Every per-ordinate flat-flux identity
  gate is therefore untouched.
* **exact on linear-in-:math:`\mu` fields**: write the level's input as
  :math:`\psi_{m,i}=a_i+b_i\,\mu_m`.  Then

  .. math::

     \psi_{1/2,i}
       &= (1-t)(a_i + b_i\mu_{m_0}) + t(a_i + b_i\mu_{m_1}) \\
       &= a_i + b_i\bigl[(1-t)\mu_{m_0} + t\mu_{m_1}\bigr]
        = a_i + b_i\,\mu_{\rm start},

  the last bracket collapsing to :math:`\mu_{\rm start}` by the
  definition of :math:`t` in :eq:`sn-err-058-edge-extrapolation`.  The
  M-M recurrence is itself a Möbius/affine map in :math:`\mu`; seeded
  with :math:`\psi_{1/2}=a+b\,\mu_{1/2}` it threads the ENTIRE
  half-angle grid exactly as :math:`\psi_{m+1/2}=a+b\,\mu_{m+1/2}` (for
  *unclamped* :math:`\tau`).  Hence the P1-class anisotropic MMS
  references — whose ansatz is exactly :math:`(A(r)+B(r)\mu)/W` — are
  *admitted* by the operator.
* **O(\Delta\mu^2)-consistent** on general smooth angular profiles —
  the same order as the angular discretisation itself.
* **linear in the input**, so the operator-algebra capabilities
  (:meth:`apply`, :meth:`apply_transpose`, dense probing) are
  preserved.  The strategy OWNS its adjoint
  (:meth:`~orpheus.sn.spatial.psi_half_angle_seed.PsiHalfAngleSeedBase.seed_adjoint`),
  a fixed linear scatter of the seed cotangent onto the two stencil
  ordinates, so a strategy swap on
  :class:`MorelMontryAngularSweep` swaps both the forward and reverse
  maps at once.

.. note::

   The :math:`\tau`-clamp
   (:math:`\tau \to \max(0.5,\min(1.0,\tau_{\rm raw}))`,
   Bailey–Morel–Chang) breaks the *exact* linear-in-:math:`\mu`
   threading wherever it is active.  This is part of the residual
   anisotropic angular floor quantified at :ref:`Issue #229
   <sn-err-058-aniso-floor>` below — NOT a wrong-fixed-point defect.

Why every gate stayed green — the blindness analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both seeds hid behind a regime in which they are exact, and the V&V
suite sat *entirely inside* that regime.  This is
``vv-principles`` Mode 7 (MMS / ansatz simplification bias) operating
not on a manufactured solution but on the *operator's own internals*.

.. list-table:: Which fields each closure seed is exact on (the blind regime)
   :header-rows: 1
   :widths: 30 34 36

   * - Closure seed
     - Exact on
     - Gate that sat in the blind regime
   * - Spatial pole-face (a)
     - flat radial :math:`\psi`
       (:math:`\psi(\Delta r/2)=\psi(0)`)
     - streaming-equilibrium L0; reflective-equilibrium k\ :sub:`∞`
   * - Angular thread (b)
     - flat-flux equilibrium
       (:math:`\Sigma_t\phi_0 = \bar Q`)
     - per-ordinate flat-flux identity (Gate 1.1); homogeneous
       reflective

**The :math:`\alpha`-dome telescoping made scalar checks blind to
(b).**  The M-M redistribution coefficients form a dome that telescopes
under the angular weight sum: :math:`\sum_n w_n\,(\alpha_{n+1/2} -
\alpha_{n-1/2}) = 0` REGARDLESS of the half-angle thread values.  Any
weight-summed (scalar-flux / particle-balance) residual therefore cannot
see a wrong half-angle thread — ``vv-principles`` anti-pattern #8
("NEVER accept particle balance as L0 evidence; require per-ordinate
residual") instantiated *inside a diagnostic*.  During the #195
investigation this telescoping made the scalar residual go
:math:`\mathcal{O}(h^2)` after fixing only (a), while the per-ordinate
residual was still :math:`\mathcal{O}(10)` — which mis-supported a
"near-singular operator / two-solutions gauge mode" hypothesis until a
dense SVD showed :math:`\sigma_{\min}\approx 0.9` (never near-singular)
and the *per-ordinate* check named the real defect.

**Historical compensation explains the Phase-D-era O(h²) reading.**  At
Phase D time the apply path measured :math:`\mathcal{O}(h^2)` under
Krylov (the premise this issue inherited), because its closure internals
compensated differently from the sweep.  The Depth-B/Wave-T matvec
rebuild changed the redistribution assembly and surfaced the latent seed
inconsistency; the SWEEP, by contrast, had ALWAYS plateaued (#195's own
SI data :math:`[0.083, 0.095, 0.098]`).  The wrong fixed point was the
sweep's all along — the same class as #98's original 35 %-at-:math:`r=0`
finding.

The three refuted intermediate hypotheses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Recording the dead paths so a future session does not re-run them
(Sphinx is the brain):

#. **"Pre-asymptotic transient"** (the issue's own premise).  Refuted by
   the :math:`n_x` 20 :math:`\to` 640 plateau — orders :math:`\to 0`, no
   refinement helps.
#. **A pure :math:`r=0`-regularity extrapolation spatial seed**
   (:math:`1.5\,\psi_0 - 0.5\,\psi_1`).  Implemented; it drove the
   *scalar* residual to :math:`\mathcal{O}(h^2)` but the solution still
   plateaued, because the dominant defect was the *angular* seed (b),
   invisible to the scalar residual by the telescoping above.  Superseded
   by the coupled-pole seed :eq:`sn-err-058-coupled-pole-continuity` for
   (a) — which is *data* rather than a one-sided stencil — once (b) was
   diagnosed.
#. **A "near-null gauge mode" theory** (apparent two-solutions paradox).
   Falsified by a dense SVD: :math:`\sigma_{\min}\approx 0.9`, the
   operator is well-conditioned.  The paradox was an artefact of the
   scalar-blind diagnostic, not a property of the operator.

Production closure decision — post-fix evidence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Post-fix (measured 2026-06-12), the isotropic curvilinear MMS solution
error collapses into a clean second-order ladder, with SI and Krylov
bit-identical:

.. list-table:: Post-fix isotropic curvilinear MMS L2 ladders (SI ≡ Krylov)
   :header-rows: 1
   :widths: 16 16 16 16 16 16

   * - :math:`n_x`
     - 20
     - 40
     - 80
     - 160
     - 320
   * - sphere :math:`\|\phi_h-A\|_{L^2}`
     - 1.49e-2
     - 3.73e-3
     - 9.28e-4
     - 2.31e-4
     - 5.74e-5
   * - sphere order
     -
     - 2.00
     - 2.01
     - 2.01
     - 2.01
   * - cylinder :math:`\|\phi_h-A\|_{L^2}`
     - 2.16e-3
     - 5.39e-4
     - 1.35e-4
     - 3.37e-5
     -
   * - cylinder order
     -
     - 2.00
     - 2.00
     - 2.00
     -

The magnitude band :math:`10^{-8} < {\rm err} < 10^{-3}` is satisfied
(sphere :math:`n_x \ge 80`, cylinder :math:`n_x \ge 40`).  SI converges
:math:`\sim 10^2\times` faster than GMRES here (sphere :math:`n_x=160`:
:math:`\sim 0.11\,\mathrm{s}` SI vs :math:`\sim 69\,\mathrm{s}` Krylov),
which is why the curvilinear default returned to source iteration.

The decisive *structural* gate is the **per-ordinate, volume-weighted**
operator-admission residual of :math:`\psi_{\rm ref}` (the scalar
residual is blind, per the telescoping above):

.. list-table:: Per-ordinate volume-weighted residual of ψ_ref under (L+C−S−B), post-fix
   :header-rows: 1
   :widths: 25 25 25 25

   * - Geometry
     - :math:`n_c=40`
     - :math:`n_c=80`
     - measured order
   * - sphere
     - 1.97e-3
     - 9.7e-4
     - :math:`\approx 1.5` (pole-adjacent bounded band under
       the :math:`r^2\,dr` weight)
   * - cylinder
     - 5.50e-5
     - 1.37e-5
     - :math:`\approx 2.0` (pointwise :math:`\mathcal{O}(h^2)`
       everywhere)

The sphere's sub-quadratic residual order is benign: the
pole-adjacent cells legitimately carry a bounded non-decaying
*pointwise* residual where the closure truncation meets the
:math:`\Delta A/V \sim 1/h` geometry factor on cells whose volume
vanishes as :math:`r^2\,dr`; the solution-error ladder above proves it
harmless.  **Bug-era** values for this gate were :math:`\mathcal{O}(10^{-1})`-class
(per-ordinate pointwise up to :math:`\pm 55` at the pole) — three-plus
orders of magnitude above the post-fix bounds, which is the margin the
ERR-058 catcher asserts.

The quadrature/truncation floor is the radial DD closure order itself;
the post-fix sphere/cylinder ladders sit at the DD design order
(2.00–2.01), so "have you tried finer quadrature?" is pre-empted — the
solution-error *is* second-order, and the only residual non-quadratic
behaviour is the volume-weighted pole band, which the solution error
does not inherit.

.. _sn-err-058-aniso-floor:

The anisotropic angular floor — deferred to Issue #229
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Resolved (2026-06-13, W1–W5).**

   This deferral is now fully treated at
   :ref:`sn-curvilinear-aniso-norm-reconciliation`.  The W1–W5
   root-cause program found the single "floor" sketched below was
   **three distinct errors** (a sphere pole-cell spatial closure, a
   sphere angular :math:`\tau`-clamp floor, and a cylinder angular
   floor), separated by a norm difference (volume-weighted L2 vs
   pointwise / L∞).  In particular, the "Open research paths" item
   below — "Unclamped-:math:`\tau` threading on a linear-in-:math:`\mu`
   shell" — was **executed** (W1): the sphere clamp was removed (it was
   mis-cited and 100 % spurious; see
   :ref:`sn-tau-clamp-vindication`), which cleaned the coarse rate but,
   surprisingly, *raised* the S16 fine floor (the prior lower floor was
   a fortuitous cancellation, not a gain).  The numbers preserved below
   are correct historical evidence; the comprehensive treatment and the
   per-error production decisions are at the reconciliation section.

The **anisotropic** curvilinear MMS (ansatz :math:`(A(r)+B(r)\mu)/W`)
dropped :math:`\sim 50\times` under the fix and is now limited by a
*fixed-quadrature angular floor*, NOT by a residual wrong-fixed-point
defect.  The mechanism: the aniso MMS imposes :math:`\psi_n` per
ordinate, so there is no angular error *at the imposed ordinates* — but
the M-M redistribution consumes half-angle THREAD values
:math:`\psi_{m\pm 1/2}` that the recurrence *interpolates*.  On an
angle-varying ansatz the thread's interpolation error is an
angular-quadrature-resolution effect: under spatial refinement at fixed
quadrature the solution converges to an angular floor, and the
pure-spatial rate + magnitude assertions cannot both hold once the
spatial error drops below it.  The floor *scales with quadrature order*
in both geometries — the structural signature confirming the
angular-thread attribution:

.. list-table:: Anisotropic angular floor vs quadrature order (post-ERR-058, SI inner)
   :header-rows: 1
   :widths: 22 24 54

   * - Case
     - Quadrature
     - Floor behaviour
   * - sphere aniso
     - S16 (shipped)
     - :math:`n_x` 10→160: ``[5.9e-2, 1.5e-2, 4.0e-3, 1.15e-3, 7.3e-4]``;
       floor :math:`\approx 7\mathrm{e}{-4}`
   * - sphere aniso
     - S32
     - err@80 = 9.5e-4, err@160 = 2.9e-4 (floor drops :math:`\sim 2.5\times`)
   * - cylinder aniso
     - :math:`n_\mu{=}4` (shipped)
     - :math:`n_x` 40→80: ``1.91e-2 → 1.90e-2`` (hard floor 1.9e-2)
   * - cylinder aniso
     - :math:`n_\mu{=}8`
     - :math:`n_x` 40→80: ``7.50e-3 → 7.39e-3``
       (floor drops :math:`\sim 2.6\times` per :math:`n_\mu` doubling)

The :math:`\tau`-clamp (above) contributes a constant to this floor by
breaking the exact linear-in-:math:`\mu` threading where active.  The
quadrature-aware retune (raise the case quadrature, or split the claim
into a pre-floor spatial-O(h²) segment + a separate
angular-convergence assertion) is `Issue #229
<https://github.com/deOliveira-R/ORPHEUS/issues/229>`_.

Infrastructure retained
~~~~~~~~~~~~~~~~~~~~~~~~~

Per the aggressive-retirement *exception* (a correct primitive that
would be needed if the obstruction is ever bypassed is kept as an
oracle), ERR-058 deletes no correct machinery:

.. list-table:: Curvilinear closure-seed primitives — status after ERR-058
   :header-rows: 1
   :widths: 38 18 44

   * - Primitive
     - Status
     - Why kept
   * - :class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`
     - **production**
     - The new ``psi_half_seed`` default; operator-consistent on
       non-flat fields.
   * - :class:`~orpheus.sn.spatial.psi_half_angle_seed.CarlsonInwardSweep`
       + :func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`
     - retained, opt-in
     - Correct *source-driven* Hébert §3.9.4 recurrence; would seed a
       future TRUE-source sweep-side closure.  Proxy-source caveat
       pinned in its docstring ``warning`` block.
   * - :class:`~orpheus.sn.spatial.psi_half_angle_seed.ZeroSeed`
     - retained, ablation
     - Reproduces the Phase B behaviour for A/B regression-safety
       comparison.
   * - Coupled-pole spatial seed in
       :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
       / :meth:`~orpheus.sn.loss_representation.transport_sweep`
     - **production**
     - The :math:`\psi(0,+\mu)=\psi(0,-\mu)` continuity; matvec + sweep
       share it (one discrete system).

Open research paths
~~~~~~~~~~~~~~~~~~~~~

Two paths could lift the anisotropic angular floor without changing the
isotropic O(h²) result:

#. **TRUE-source-driven sweep-side seed.**  Replace the
   :class:`AngularEdgeExtrapolation` *input-field* extrapolation with the
   canonical Hébert recurrence
   :func:`~orpheus.sn.spatial.psi_half_angle_seed.carlson_inward_sweep_from_source`
   driven by the genuine within-group source
   :math:`\bar Q_i = \sum_\ell \tfrac{2\ell+1}{2}Q_\ell(r_i)(-1)^\ell`
   (the full Legendre fold, not the :math:`\Sigma_t\phi_0` proxy).  On
   the sweep path the true source is available; this would make the
   *starting-direction transport* exact rather than the *input-field
   value* exact.  Likely diagnostic probe: the per-ordinate residual at
   the most-inward ordinate of an aniso ansatz, holding spatial mesh
   fixed and sweeping quadrature order.
#. **Unclamped-:math:`\tau` threading on a linear-in-:math:`\mu` shell.**
   The exact linear-:math:`\mu` threading (above) holds only for
   unclamped :math:`\tau`; quantify the clamp's contribution to the
   floor and, where the cell is well-resolved
   (:math:`\tau_{\rm raw}\in[0.5,1.0]`), thread unclamped to recover the
   exact P1 admission.  Likely probe: the
   :ref:`Issue #229 <sn-err-058-aniso-floor>` floor table with the clamp
   disabled on resolved cells.

Session trail (V&V audit trail)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **ERR-058 catalogue narrative**:
  ``.claude/skills/vv-principles/error_catalog.md`` (§ ERR-058) — the
  authoritative two-manifestation mechanism + post-fix evidence.
* **Re-scope record**: `Issue #195
  <https://github.com/deOliveira-R/ORPHEUS/issues/195>`_ comments
  (2026-06-12) — the premise refutation and the decisive probe-3
  residual evidence.
* **Diagnostics**: ``derivations/diagnostics/diag_195_probe{1,2,3}_*.py``
  (the plateau / error-profile / operator-admission probes), promoted to
  the gate ``tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py``.
* **Investigator memo**:
  ``.claude/agent-memory/numerics-investigator/issue_195_root_cause_2_pole_closure.md``.

Verification chain
~~~~~~~~~~~~~~~~~~~

The ERR-058 fix is pinned by, in order of structural decisiveness:

#. :func:`tests.sn.verification.mms.test_curvilinear_operator_admits_mms.test_operator_admits_isotropic_mms_per_ordinate`
   (``@pytest.mark.l1`` + ``catches("ERR-058")``) — the fast
   per-ordinate volume-weighted operator-admission gate (the structurally
   decisive check, immune to the telescoping blindness).
#. :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_spherical_mms_converges_second_order`
   and
   :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_cylindrical_mms_converges_second_order`
   (``catches("ERR-058")``) — the end-to-end L1 ladders whose
   ``xfail`` markers came off with this fix; they ``verifies`` the
   :eq:`sn-mms-spherical-psi` / :eq:`sn-mms-spherical-qext` /
   :eq:`sn-mms-cylindrical-psi` / :eq:`sn-mms-cylindrical-qext` labels.
#. The flat-flux and streaming-equilibrium gates pin the flat-field
   exactness BOTH fixes preserve (so they did not regress).
#. :func:`tests.sn.operators.test_g_adjoint_reciprocity` — pins the
   strategy-owned seed adjoints.

.. note::

   **vv-status (eq-labels added by this section).**  The labels
   :eq:`sn-err-058-coupled-pole-continuity`,
   :eq:`sn-err-058-proxy-source`, and
   :eq:`sn-err-058-edge-extrapolation` are *structural / representational*
   identities (the pole-continuity boundary condition; the falsified
   proxy-source definition; the operator-consistent edge-extrapolation
   map).  They are NOT solver claims (no eigenvalue / flux-value claim).
   Per the vv-status discipline they are ``documented`` — the
   verifiable content is the per-ordinate operator-admission gate
   (``catches("ERR-058")``) plus the strategy-owned adjoint
   bit-identity, named in the verification chain above.

.. _sn-issue-196-eigenvalue-equivalence:

Issue #196 — eigenvalue-path SI≡Krylov verification and the permanent gate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: Status — manifestation #7 verified and pinned (Issue #196 CLOSED, 2026-06-13)
   :class: important

   ERR-058 (#195, above) replaced the wrong shared closure seeds; **#196
   is the verification and regression-gate step** that confirms the
   replacement closed ERR-026 manifestation #7 at the *eigenvalue*
   level and locks the closure against re-introduction.  This was the
   LAST open manifestation of ERR-026 — with #196 closed, the
   curvilinear-SN wrong-fixed-point family is **formally retired**.

The two-layer history (why the L0 close did not suffice)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Manifestation #7 names a specific defect: a residual :math:`\mathcal{O}(h)`
SI-vs-Krylov WDD spatial-closure asymmetry on curvilinear SN.  Pre-fix,
the source-iteration inner (drives the curvilinear sweep) and the Krylov
inner (drives the apply-matvec) produced eigenvalues differing at
:math:`\mathcal{O}(h)`: **0.286 %** on ``sphere_2g_3reg`` at n=40, **~30
%** per-cell on eigenvector shape, the gap halving under mesh refinement.

Two closures touched this defect, and the honest history distinguishes
them:

* **ERR-048** (Phase G Step 2, 2026-05-13) closed only the **L0
  flat-field** twin-agreement.  It patched the SI sweep
  (``_sweep_1d_spherical`` / ``_sweep_1d_cylindrical``) to MATCH the
  apply-matvec conventions on the homogeneous streaming-equilibrium
  gauntlet (pole-face WDD IC mirror + Carlson seed :math:`\bar Q`
  normalisation).  The **L1 heterogeneous eigenvalue**
  :math:`\mathcal{O}(h)` asymmetry that manifestation #7 names
  **PERSISTED** — which is exactly why #196 stayed OPEN — because the
  shared closure seeds were still *exact-on-flat /
  O(1)-wrong-on-non-flat* (the ERR-058 defect).
* **ERR-058** (Issue #195, 2026-06-12) was the TERMINAL fix.  It
  replaced the shared closure seeds with correct ones (the coupled-pole
  spatial seed :math:`\psi(0,+\mu)=\psi(0,-\mu)` and the
  :class:`~orpheus.sn.spatial.psi_half_angle_seed.AngularEdgeExtrapolation`
  half-angle seed), making BOTH inner solvers operate on the SAME
  correct discrete operator.

.. _sn-issue-196-bit-identical-vs-floor:

Bit-identical (fixed-source) vs floor-equivalent (eigenvalue)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The post-fix SI :math:`\equiv` Krylov agreement is **not the same kind
of agreement** on the two solver entry points, and the distinction is
load-bearing:

* **Fixed-source** (:func:`~orpheus.sn.solver.solve_sn_fixed_source` on
  the curvilinear MMS ladders): SI :math:`\equiv` Krylov is
  **BIT-IDENTICAL**.  Post-unification the sweep and the matvec are ONE
  discrete operator on one quadrature; the within-group inner
  (``L.solve`` vs Krylov-on-``apply``) realises the *same* :math:`L^{-1}`
  arithmetic, so the two paths return the same bits.
* **Eigenvalue** (:func:`~orpheus.sn.solver.solve_sn` with
  ``inner_solver="source_iteration"`` vs ``"krylov"``): SI
  :math:`\equiv` Krylov to the **ITERATION FLOOR** (~:math:`1.9\mathrm{e}{-11}`
  in :math:`k_{\text{eff}}`, ~:math:`2.4\mathrm{e}{-10}` in flux shape),
  **NOT bit-identical**.  The eigenvalue solve wraps the inner in power
  iteration; SI and Krylov are *different iteration schemes* that
  converge to the **same correct fixed point** only to ~``inner_tol``.
  Same physics, not the same arithmetic.

Confusing the two would mis-state the verification claim.  The earlier
close-out's "SI :math:`\equiv` Krylov bit-identical on the curvilinear
MMS ladders" is correct **for the fixed-source ladders specifically**;
the eigenvalue verification below is *floor-equivalence*, which is the
right and sufficient claim for an eigenvalue solve.

Measured eigenvalue-path equivalence (Issue #196)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

All values measured 2026-06-12 under tight iteration tolerances for
BOTH inner solvers (``keff_tol=1e-12``, ``flux_tol=1e-10``,
``inner_tol=1e-10``).  The eigenvalue snapshot cases are the exact
acceptance cases of #196:

.. list-table:: SI≡Krylov eigenvalue equivalence — the bug-era heterogeneous snapshot cases
   :header-rows: 1
   :widths: 30 22 22 26

   * - Case
     - :math:`|\Delta k|` (post-fix)
     - max :math:`|\Delta\varphi|_{L\infty}` (post-fix)
     - Bug-era (pre-ERR-058)
   * - ``sphere_2g_3reg_dd_n40``
     - :math:`4.68\mathrm{e}{-12}`
     - :math:`5.88\mathrm{e}{-11}`
     - :math:`|\Delta k|=3.9\mathrm{e}{-3}` (0.286 %), shape ~30 %
   * - ``cyl_2g_3reg_LS4_dd_n40``
     - :math:`1.91\mathrm{e}{-11}`
     - :math:`4.32\mathrm{e}{-11}`
     - same :math:`\mathcal{O}(h)` family, gap halving under refinement

The homogeneous (k_inf-degenerate, flat-flux) curvilinear snapshots
agree at the rounding floor — as expected, since on a flat eigenmode the
redistribution terms null and SI/Krylov differ only by accumulated
round-off:

.. list-table:: SI≡Krylov eigenvalue equivalence — homogeneous (flat-flux) snapshots
   :header-rows: 1
   :widths: 38 26 26

   * - Case
     - :math:`|\Delta k|`
     - relative :math:`\varphi` diff
   * - ``sphere_2g_homogeneous_dd_n20``
     - :math:`6.92\mathrm{e}{-13}`
     - :math:`2.15\mathrm{e}{-10}`
   * - ``cyl_1g_homogeneous_LS4_dd_n20``
     - :math:`2.22\mathrm{e}{-16}`
     - :math:`2.27\mathrm{e}{-11}`
   * - ``cyl_1g_homogeneous_product_dd_n20``
     - :math:`6.66\mathrm{e}{-16}`
     - :math:`1.10\mathrm{e}{-10}`

.. note::

   The homogeneous cases agree to the floor but, on their own, supply
   **no** evidence for the curvilinear closure — a flat eigenmode is
   degenerate (``flat = flat``; 1-group :math:`k=\nu\Sigma_f/\Sigma_a`
   is flux-shape independent, vv-principles anti-patterns #3/#4).  The
   load-bearing evidence is the **heterogeneous 2-group** cases above,
   where the flux is genuinely non-flat and the angular-redistribution
   terms are exercised.

The permanent regression gate (Issue #196)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The manifestation-#7 catcher is
:func:`tests.sn.eigenvalue.test_keff_curvilinear.test_si_krylov_eigenvalue_equivalence_sphere`
and
:func:`tests.sn.eigenvalue.test_keff_curvilinear.test_si_krylov_eigenvalue_equivalence_cylinder`,
each carrying ``@pytest.mark.catches("ERR-026")``.  Configuration:
heterogeneous 2-group fuel|moderator (region A inner, region B outer,
n=10+10), solved twice under ``inner_solver="source_iteration"`` and
``inner_solver="krylov"`` at the tight tolerances above.  The gate
asserts:

.. list-table:: Manifestation-#7 gate thresholds vs measured vs bug-era
   :header-rows: 1
   :widths: 26 24 24 26

   * - Quantity
     - Asserted bound
     - Measured (post-fix)
     - Bug-era (would trip)
   * - sphere :math:`|\Delta k|`
     - :math:`< 1\mathrm{e}{-7}`
     - :math:`1.9\mathrm{e}{-11}`
     - :math:`3.9\mathrm{e}{-3}`
   * - sphere per-group :math:`|\Delta\varphi|_{L\infty}`
     - :math:`< 1\mathrm{e}{-6}`
     - :math:`2.4\mathrm{e}{-10}`
     - ~30 %
   * - cylinder :math:`|\Delta k|`
     - :math:`< 1\mathrm{e}{-7}`
     - :math:`1.1\mathrm{e}{-11}`
     - same family
   * - cylinder per-group :math:`|\Delta\varphi|_{L\infty}`
     - :math:`< 1\mathrm{e}{-6}`
     - :math:`2.6\mathrm{e}{-11}`
     - ~30 %

A **non-flat-flux guard** (group-0 radial ``max/min > 1.2``) precedes
the equivalence assertion so the test cannot pass vacuously on a
degenerate flat mode — sphere radial ``max/min`` = 3.34, cylinder =
1.67, both well above the guard.  The bug-era values (3.9e-3 / ~30 %)
exceed the asserted bounds by **4–5 orders of magnitude**, so the gate
would have tripped loudly on the pre-fix code.

.. note::

   **Runtime-mode discipline (vv-principles anti-pattern #8).**  The
   canonical ORPHEUS invocation is ``python -O``, under which bare
   ``assert`` statements are stripped to no-ops.  These gates assert via
   bare ``assert`` inside the *collected test module*, which pytest
   rewrites at collection time so the asserts fire under ``-O``.  This
   was confirmed empirically: a negative control with a
   :math:`1\mathrm{e}{-15}` tolerance failed as required under ``-O``.

Structural-independence — why SI≡Krylov alone is not the whole proof
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SI :math:`\equiv` Krylov is **twin-path agreement** — necessary but not
sufficient (vv-principles L11: two implementations agreeing is
cross-implementation evidence, not correctness evidence).  Both inner
solvers could in principle converge to the same *wrong* fixed point.
The independent ground that makes the closure a *correctness* claim, not
merely a *consistency* claim, comes from two structurally-independent
legs:

* The **k_inf homogeneous legs** — on a uniform reflective infinite
  medium :math:`k_\infty=\nu\Sigma_f/\Sigma_a` is an analytical
  (closed-form) eigenvalue the SN snapshots must reproduce.
* The **Variant-α Green's-function cross-check**
  (:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_e_trajectory_resolvent_flux_shape_crosscheck`),
  now a plain L1 test (xfail removed), which compares the SN flux-shape
  snapshot against the composite-GL trajectory-resolvent reference
  within 8 % (sphere) / 12 % (cylinder).  This reference is a
  semi-analytical pillar structurally independent of the SN sweep, so
  agreement pins the *converged-to value*, not just twin-path
  consistency.

Production-decision record — curvilinear default reverted to SI
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The curvilinear :attr:`~orpheus.sn.solver.SNSolver.inner_solver` default
is now ``"source_iteration"``, **reverted from the Phase-D Krylov
flip**.  The Phase-D flip existed ONLY because the sweep's fixed point
was wrong; ERR-058 made it correct, so SI is restored as the default —
it is :math:`\sim 10^2\times` faster than GMRES (no restart) and now
equivalent (bit-identical fixed-source / floor-equivalent eigenvalue).

Crucially, **neither of the old Phase-F closures was taken**:

* *Option (a) — make SI bit-identical to Krylov by refining the WDD
  closure* presupposed the seed was correct and only the spatial
  arithmetic differed.
* *Option (b) — flip the default to Krylov* presupposed the Krylov fixed
  point was the correct one and SI's was a discretisation artefact.

Both presupposed the *shared* fixed point was correct and only the
arithmetic differed.  The terminal diagnosis (ERR-058) showed the shared
fixed point **itself** was wrong on non-flat fields; once the seeds were
fixed, both inner solvers reach the same *correct* fixed point and the
"choose between them" framing dissolves — SI is restored on speed alone.


.. _sn-curvilinear-aniso-norm-reconciliation:

The curvilinear anisotropic-MMS "floor", reconciled (W1–W5)
===========================================================

.. admonition:: Key Facts — curvilinear anisotropic SN
   :class: important

   - **The single "#229 floor" was three distinct errors**, separated
     by a *norm difference*: the production gates measure a
     **volume-weighted L2** :math:`\sqrt{\sum_i V_i\,\Delta_i^2}`; the
     root-cause probes measured **pointwise / L∞**.  An error
     concentrated at the :math:`r \to 0` pole cell is loud in L∞ and
     near-silent in L2.
   - **(a) Sphere central-cell** :math:`\mathcal{O}(h)` **spatial
     closure (#233)** — L∞-only; :math:`\sim 75\,\%` an MMS
     midpoint-vs-shell-volume-average comparison artifact +
     :math:`\sim 25\,\%` genuine but **inherent** first order.  WONTFIX
     for diamond difference (Hébert §3.9.4 / Stacey §9.9 use plain
     diamond at the central cell).  See
     :ref:`sn-pole-cell-spatial-closure`.
   - **(b) Sphere angular** :math:`\tau`-**clamp floor** — fixed by W1
     (the clamp was mis-cited and 100 % spurious on physical fields).
     The sphere now uses the raw Bailey-Morel-Chang 2010 Eq. 43 weight.
     See :ref:`sn-tau-clamp-vindication`.
   - **(c) Cylinder angular floor (#229)** — the half-angle-thread
     INTERPOLATION floor; scales with the **azimuthal** quadrature
     :math:`n_\varphi`, structurally blocked (needs a 2-D
     :math:`(\eta, \varphi)` closure).  The sphere has a pre-floor
     :math:`\mathcal{O}(h^2)` window (clean at S32); the cylinder does
     **not** at any practical quadrature.  See
     :ref:`sn-cylinder-angular-floor`.
   - **Two unrelated "anisotropic" paths** (Issue #9): Path-(I) =
     geometric angular redistribution :math:`(1-\mu^2)/r\,\partial_\mu`
     (the M-M :math:`\alpha`-dome, P0-only — what #229 concerns);
     Path-(II) = :math:`P_1{+}` Legendre SCATTERING
     :math:`R\,\Lambda\,M` (geometry-agnostic).  See
     :ref:`sn-p1-scattering-curvilinear`.
   - **Norm gotcha**: a convergence-rate gate on a volume-weighted L2
     norm CANNOT see a pole-cell defect (the pole sits at one cell of
     :math:`V \sim h^3` → :math:`\sqrt{V} \sim h^{1.5}` →
     :math:`\sim h^{2.5}` contribution → subdominant).  A pointwise /
     L∞ probe is required to surface it.

This section closes the curvilinear-anisotropic-SN investigation
program (W1–W5, branch ``fix/curvilinear-aniso-pole-and-clamp``,
2026-06-13).  It is the sequel to the ERR-058 / #195 / #196 curvilinear
*isotropic* closure-seed family above; that family fixed the
wrong-fixed-point class (now formally retired), and what remained was
the *anisotropic* floor — which this program resolved into three
distinct, separately-actionable errors.

The headline — one floor was three errors, settled by a norm difference
-----------------------------------------------------------------------

The ERR-058 close-out (above) deferred a residual "anisotropic angular
floor" to Issue #229, citing a single
:ref:`floor table <sn-err-058-aniso-floor>`.  The W1–W5 root-cause
program found that the apparent single floor was **three structurally
distinct errors**, and the reason they had been conflated is a **norm
difference** in how two independent investigations measured the same
solves:

* The verification gates (test-architect) measure the **volume-weighted
  L2** norm :math:`\|\Delta\|_{2,V} = \sqrt{\sum_i V_i\,\Delta_i^2}` —
  the natural norm for a finite-volume scheme whose unknown is a
  cell-volume average.
* The diagnostic probes (numerics-investigator) measured **pointwise /
  L∞** — :math:`\max_i |\Delta_i|`.

The two norms weight the :math:`r \to 0` pole cell completely
differently.  Under the spherical volume weight :math:`V \sim h^3`, a
fixed pointwise error at the single pole cell contributes
:math:`\sqrt{V} \sim h^{1.5}` to the L2 sum, so an L∞-:math:`\mathcal{O}(h)`
pole error appears as :math:`\sim h^{2.5}` in L2 — **subdominant to the
interior** :math:`\mathcal{O}(h^2)`, hence invisible.  This is exactly
why the production L2 gate stayed green throughout while a pointwise
probe found a first-order pole cell.

.. list-table:: The three errors behind the "#229 floor"
   :header-rows: 1
   :widths: 6 38 22 16 18

   * - #
     - Error
     - Dominant norm
     - Quadrature-scaling?
     - Status
   * - (a)
     - Sphere pole-cell spatial closure :math:`\mathcal{O}(h)` at
       :math:`r \to 0`
     - L∞ / pointwise central flux (diluted in L2 by :math:`V \propto
       r^2`; invisible in :math:`k_{\rm eff}`)
     - no (spatial)
     - **#233 — documented inherent limitation (ERR-059, WONTFIX-for-DD)**
   * - (b)
     - Sphere angular :math:`\tau`-clamp floor (:math:`\sim 7\mathrm{e}{-4}`
       @ S16)
     - volume-weighted L2 at fine mesh
     - yes (angular)
     - **fixed (W1 unclamp)**
   * - (c)
     - Cylinder angular floor
     - both
     - yes (azimuthal :math:`n_\varphi`)
     - **structurally blocked (#229)**

The remainder of this section treats each error in turn, then the two
unrelated anisotropic paths (Issue #9).

.. _sn-tau-clamp-vindication:

(b) The :math:`\tau`-clamp vindication (W1)
-------------------------------------------

The spherical Morel--Montry weighted-diamond weight is

.. math::
   :label: sn-tau-mm-raw

   \tau_n^{\rm raw}
       \;=\; \frac{\mu_n - \mu_{n-1/2}}{\mu_{n+1/2} - \mu_{n-1/2}}
       \;\in\; [0, 1],

the **unique** weight exact for an angular flux linear in :math:`\mu`
([BaileyMorelChang2010]_ Eq. 43; the same object as
:eq:`mm-weights`).  The production code had wrapped it in a
:math:`[\tfrac12, 1]` clamp,
:math:`\tau_n = \mathrm{clip}(\tau_n^{\rm raw}, \tfrac12, 1)`, cited
to Lewis & Miller §4.5.

W1 established, by three independent lines of evidence, that the clamp
is **mis-cited and 100 % spurious on physical fields**:

#. **Literature.** [BaileyMorelChang2010]_ state the admissible range
   is :math:`\tau \in [0, 1]` and recommend *exactly* the unclamped
   :math:`\tau^{\rm raw}` (their Eq. 43) as the unique exact-on-linear
   weight; Hébert §3.9.4 uses pure diamond (:math:`\tau = \tfrac12`),
   no clamp.  Lewis & Miller §4.5 does **not** prescribe the
   :math:`[\tfrac12, 1]` clamp — the citation was wrong.
#. **Positivity is never needed.** On every realistic converged solve
   (smooth MMS, homogeneous eigenvalue :math:`k_{\rm eff} = 1`, thick
   absorber) there are ZERO negative half-angle fluxes, clamped or
   unclamped — every clamp activation is spurious (measured: 160 / 320
   / 80 / 240 activations across stress configs, 0 protective).  The
   half-flux negativity that *does* transiently appear in early SI
   iterates is inherited from a negative *input* :math:`\psi` and the
   clamp barely reduces it.  On Gauss--Legendre quadrature
   :math:`\tau^{\rm raw} \in [0.39, 0.61]` (never 0), so the unclamped
   weight is always interior to :math:`[0, 1]`.
#. **Stability without it.** Unclamped sphere source iteration
   converges with strictly positive, finite scalar flux on every
   stress config (thick absorber, near-vacuum, :math:`c = 0.999`, S64);
   the clamp costs a few SI iterations on low-scattering problems but is
   dispensable for stability.

.. admonition:: The architectural reason the static removal is correct
   :class: note

   A *dynamic* negative-flux fixup (where :math:`\tau` depends on the
   iterate :math:`\psi`) would make the streaming operator
   **nonlinear**, breaking the linear-Krylov matvec and the SI ≡ Krylov
   twin identity (Pattern-2 discipline, Cardinal Rule 2).  Because the
   fixup is *never needed* on physical fields, the principled W1 fix is
   to **drop the clamp** (a config-time, static change) and use the
   linear unclamped :math:`\tau^{\rm raw}`.  ``tau_mm`` is single-sourced
   in :func:`~orpheus.geometry.reduced_operator.spherical_streaming` and
   inherited by every consumer (the SI sweep and the Krylov matvec
   both), so both twins stay linear and stay identical.

**Geometry split.**  W1 removed the clamp for the **sphere only**.  The
**cylinder keeps** it: product / level-symmetric quadratures place the
most-inward azimuthal ordinate exactly on :math:`\eta = -\sin\theta`,
giving :math:`\tau^{\rm raw} = 0` **exactly** (bit-exact, not "near
zero") — an unclamped recurrence divides by zero there.  This is a
genuine *structural* singularity the sphere provably lacks; the
cylinder's real fix is a 2-D :math:`(\eta, \varphi)` closure
(:ref:`sn-cylinder-angular-floor`), not unclamping.  See
:eq:`morel-montry-clamp` in :doc:`structured_geometry` for the
equation-of-record carrying both branches.

**Mixed accuracy signature (the gotcha).**  Unclamping does NOT
uniformly improve the anisotropic solve.  It *cleans the coarse
convergence rate* (sphere S16 coarse orders 1.978 → 1.995) but *raises
the S16 fine-mesh floor* (:math:`7.3\mathrm{e}{-4} \to 1.2\mathrm{e}{-3}`):

.. list-table:: W1 sphere aniso MMS, matched-quadrature S16 (volume-weighted L2)
   :header-rows: 1
   :widths: 20 40 40

   * - :math:`n_x`
     - Clamped (pre-W1)
     - Unclamped (post-W1)
   * - 10 → 40
     - coarse orders 1.979 / 1.978
     - coarse orders 1.995 / 1.999
   * - 80
     - 1.16e-3
     - 1.40e-3
   * - 160 (floor)
     - **7.3e-4**
     - **1.2e-3**

The lower *clamped* floor was a **fortuitous cancellation**, not a
genuine accuracy gain — the clamp's constant bias happened to partly
offset the angular-thread interpolation floor at S16.  Removing it
exposes the true #229 floor (next subsection), which is what the
unclamped weight should converge to.  Iso solves are unchanged in real
arithmetic (the clamp is silent on flat-in-:math:`\mu` fields) but
**not bit-identical** at IEEE-754: the closure
:math:`(\overline\psi - (1-\tau)\psi_{\rm in})/\tau` returns
:math:`\psi` exactly only :math:`\sim 81\,\%` of the time and within
1 ULP otherwise (reduction-order non-associativity), so the converged
homogeneous-reflective sphere drifts :math:`|\Delta k| = 2.3\mathrm{e}{-13}`,
:math:`\max|\Delta\phi| = 4.4\mathrm{e}{-13}` — an FP-tail, anchored to
the closed-form :math:`k_\infty = 1.875`.  One snapshot
(``sphere_2g_3reg_dd_n40``, genuinely non-flat) was regenerated
(:math:`k\;1.380766 \to 1.381001`); the two flat snapshots drift only
in the FP tail and were not regenerated.

**W1 gates.** ``tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py``
(closure-unit :math:`\tau`-independence on flat fields; converged
homogeneous-reflective iso anchored to :math:`k_\infty`; unclamped
positivity) + the W1 ``@slow`` aniso gates appended to
``tests/sn/verification/mms/test_curvilinear_aniso_convergence.py``
(the S32 clean-:math:`\mathcal{O}(h^2)` full-ladder claim; the S16
coarse-rate-cleaner-unclamped discriminator; the floor-scales-with-
quadrature pin).  Landed in commit ``b2d8a6d``.

.. _sn-cylinder-angular-floor:

(c) The cylinder angular floor (#229) — structurally blocked
------------------------------------------------------------

The anisotropic ansatz :math:`\psi_{\rm chosen} = (A(r) + B(r)\,\mu)/W`
imposes :math:`\psi_n` per ordinate, so there is **no angular error at
the imposed ordinates**.  But the M-M redistribution consumes
half-angle THREAD values :math:`\psi_{m\pm 1/2}` that the recurrence
**interpolates** (they are not imposed).  On an angle-varying ansatz the
thread's interpolation error is an angular-quadrature-resolution effect:
under spatial refinement at fixed quadrature the solution converges to
an angular floor, and a pure-spatial-rate assertion cannot hold once
the spatial error drops below it.

**The cylinder floor scales with the AZIMUTHAL quadrature**
:math:`n_\varphi`, **not the polar** :math:`n_\mu`.  This is the
load-bearing physical fact (and a correction to an earlier mislabel):
the radial direction cosine is :math:`\eta = \sin\theta\,\cos\varphi`,
so the M-M thread marches in azimuth :math:`\varphi` *per polar
:math:`\mu`-level*.  Measured at :math:`n_x = 80`:

.. list-table:: Cylinder aniso floor vs azimuthal quadrature (:math:`n_x = 80`, volume-weighted L2)
   :header-rows: 1
   :widths: 25 25 50

   * - :math:`n_\varphi`
     - L2 error
     - Behaviour
   * - 8
     - 1.90e-2
     - hard floor
   * - 16
     - 7.37e-3
     - drops :math:`2.58\times`
   * - 32
     - 3.10e-3
     - drops :math:`2.38\times`

while :math:`n_\mu` (polar) refinement at fixed :math:`n_\varphi`
leaves the floor **flat**.

**Why it is structurally blocked.**  Product and level-symmetric
quadratures carry **duplicate azimuthal** :math:`\eta`: ordinates come
in :math:`\pm\varphi` symmetry pairs with the same :math:`|\eta|` but
opposite :math:`\xi` (e.g. :math:`\varphi = \pi/4` and
:math:`\varphi = 7\pi/4` both give
:math:`\eta = \sin\theta/\sqrt 2`).  The M-M thread marches in
:math:`\eta` alone, so a field whose true variation is in the full
:math:`(\eta, \varphi)` plane is **not threadable exactly** by a 1-D
:math:`\eta`-march — a structural mismatch, not a tuning problem.  No
partition (midpoint / cumulative-weight / ordinate-interior) gives
:math:`\tau^{\rm raw} \in [\tfrac12, 1]` with bounded edges; the
cumulative-weight partition is exact on level-symmetric but needs
:math:`\tau^{\rm raw} \in [-4.5, 5.5]` (edges outside the level).
Closing the cylinder floor requires a genuine 2-D
:math:`(\eta, \varphi)` angular closure — **out of scope**, tracked by
`Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_.

**The sphere–cylinder asymmetry.**  The sphere DOES have a pre-floor
:math:`\mathcal{O}(h^2)` window: at S16 the coarse orders clear 1.99 and
the floor (:math:`\approx 2.9\mathrm{e}{-4}` at S32, :math:`n_x = 160`)
sits below the segment's finest spatial error, so the clean
second-order window extends to :math:`n_x = 80` at S32.  The cylinder
has **no** such window at *any* practical quadrature — even
:math:`n_\mu = 16` (:math:`N = 512`) reaches only order 1.80 on the
coarsest :math:`\{5, 10, 20\}` segment before the angular floor
dominates.  The mathematics, not runtime, is the blocker.

**W3 gate retune (Issue #229).**  Per the vv-principles anti-pattern
"a claim that cannot hold MUST NOT be asserted; pin what IS true
instead", W3 removed all five aniso xfail markers and migrated the six
equation labels to green tests:

* **Sphere** ``test_sn_spherical_aniso_mms_converges_second_order`` →
  coarse-segment ``orders[:2] > 1.9`` + magnitude band
  :math:`1\mathrm{e}{-8} < \mathrm{err}[-1] < 5\mathrm{e}{-3}`
  (loosened from :math:`1\mathrm{e}{-3}` because the W1 unclamp removed
  the fortuitous-cancellation lower floor) + ``catches("ERR-026")``.
* **Cylinder** ``test_sn_cylindrical_aniso_mms_converges_second_order``
  → floor band :math:`1\mathrm{e}{-3} < \mathrm{err}[-1] <
  5\mathrm{e}{-2}`, **no rate claim** (the floor dominates).  The
  cylinder phase-C spatial test was **repurposed** into
  ``test_cyl_aniso_floor_scales_with_quadrature``
  (:math:`\mathrm{err}(n_\varphi{=}16) < \mathrm{err}(n_\varphi{=}8)/2`
  — the verified-floor second claim that pins the angular attribution).
* The sphere prescribed-inflow redistribution test dropped its
  strict-xfail and rate gate (band :math:`1\mathrm{e}{-8} <
  \mathrm{err} < 5\mathrm{e}{-3}` + a kept converged-value
  ``assert_allclose(2e-2)``).

Landed in commit ``679a1e6`` (audit exit 0; the
:eq:`sn-mms-spherical-psi` / ``-qext`` / :eq:`sn-mms-cylindrical-psi` /
``-qext`` labels and the two spatial-convergence labels are now all
green-tested).

.. _sn-pole-cell-spatial-closure:

(a) The sphere/cylinder pole-cell :math:`\mathcal{O}(h)` spatial closure (#233, ERR-059)
----------------------------------------------------------------------------------------

This is the **new** discovery of the program and the one *surviving*
manifestation in the curvilinear-SN family.  The curvilinear scalar
flux is **first-order** :math:`\mathcal{O}(h)` at the :math:`r \to 0`
central cell in the pointwise / L∞ norm — distinct from #168 (outer
face, CLOSED), ERR-058 (the closure seed, CLOSED), and the angular
floors above.  It decomposes into three parts, **none** of which
warrants a code fix.

Decomposition
~~~~~~~~~~~~~

**Part 1 — :math:`\sim 75\,\%` MMS comparison artifact (not a solver
bug at all).**  The production spherical MMS evaluates the source at the
cell MIDPOINT ``mesh.centers`` and compares
:math:`\phi_{\rm solver}` against :math:`\phi_{\rm exact}(\text{midpoint})`.
But the spherical DD discrete unknown **IS** the cell-volume average

.. math::
   :label: sn-pole-cell-shell-average

   \overline{\phi}_{n,i}
       \;=\; \frac{4\pi}{V_i}\int_{r_{i-1/2}}^{r_{i+1/2}} r^2\,\phi_n(r)\,dr

([Hebert2009]_ Eq. 3.430 — the unknown is *defined* as the shell
average, not a point value; the diamond relation Eq. 3.431 relates it to
the face fluxes).  Under :math:`r^2\,dr` weighting the volume-average and
the midpoint point-value differ by :math:`\mathcal{O}(h)` at the pole
cell, because :math:`r_{\rm lo} = 0` maximally skews the weight (the
volume-centroid sits at :math:`\tfrac34 h`, not :math:`\tfrac12 h`).
Using the *shell-averaged* source AND comparing against the
*shell-volume-average* drops the pole error :math:`\sim 4\times`
(:math:`0.0212 \to 0.00497`) — confirming the bulk of the apparent
error is a comparison subtlety, not solver truncation.

**Part 2 — :math:`\sim 25\,\%` genuine but LITERATURE-ACCEPTED INHERENT
first order.**  Even the fully consistent finite-volume MMS (shell-avg
source + shell-avg reference) leaves the pole at clean
:math:`\mathcal{O}(h^{1.00})`.  The root cause: at :math:`r_{\rm lo} = 0`
the inner face area :math:`A(0) = 0`, so the diamond closure
:math:`\overline\psi = \tfrac12(\psi_{\rm in} + \psi_{\rm out})` gives
:math:`\psi_{\rm out} = 2\overline\psi`, **over-predicting the pole
outer face by exactly +50 %** (mesh-independent rel. error 0.5000), while
the true face is :math:`A(h)` and :math:`2\langle A\rangle_{\rm vol} =
2\cdot\tfrac34 A(h) = 1.5\,A(h)`.  Deeper still, the conservative
*balance itself* is inconsistent at the pole: fed the EXACT cell average
and EXACT inflow it solves for an outer face :math:`-46\,\%` wrong, and
the residual-per-volume plateaus mesh-independently — because
:math:`A_{\rm in} = 0` degenerates the streaming surface integral while
:math:`V \sim h^3`.

[Hebert2009]_ §3.9.4 and Stacey §9.9 **both** use exactly this plain
diamond + Carlson-starting-direction + symmetry scheme at the central
cell with **no special** :math:`\mathcal{O}(h^2)` **closure, and
neither flags reduced order there**.  First-order at the single pole
cell is the accepted, unflagged behaviour of the standard scheme.

**Part 3 — NOT cleanly fixable by a local closure.**  W2 tested the
volume-weighted linear reconstruction
:math:`\overline\psi = \beta\,\psi_{\rm out} + (1-\beta)\,\psi_{\rm in}`
with :math:`\beta = \tfrac34` at the pole (the value that makes
:math:`\overline\psi` :math:`\mathcal{O}(h^3)`-consistent against
:math:`\langle A\rangle_{\rm vol}` at exact faces).  Validated
end-to-end with a faithful production-sweep monkeypatch (and a
:math:`\beta = \tfrac12`-identity regression guard verified to
:math:`3\mathrm{e}{-16}`): :math:`\beta = \tfrac34` does **NOT** restore
:math:`\mathcal{O}(h^2)` — the pole stays :math:`\mathcal{O}(h)`,
magnitude slightly *worse* (:math:`0.0050 \to 0.0106`), and a full-mesh
:math:`\beta` degrades the interior.  Closure-consistency at exact faces
:math:`\neq` fixed-point accuracy: the propagated face flux couples back
through the balance.  A genuine fix needs a non-local higher-order
central-cell reconstruction the canon does not provide — a linear-
discontinuous (Issue #6), cell-update (#158), or nodal scheme
([WuXieFischer1999]_ NSE 133).

Why it is invisible to L2 and to :math:`k_{\rm eff}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The production ``test_sn_spherical_mms_converges_second_order`` uses the
volume-weighted L2 norm.  The pole :math:`\mathcal{O}(h)` at one cell of
:math:`V \sim h^3` contributes :math:`\sqrt V \sim h^{1.5}` →
:math:`h^{2.5}` to L2 — subdominant.  Both midpoint and volume-average
L2 references converge clean :math:`\mathcal{O}(h^{2.00})`; only the L∞
(pole) is :math:`\mathcal{O}(h)`.  For :math:`k_{\rm eff}`: a reflective
sphere recovers :math:`k_\infty = 1.875` exactly, mesh-independent; a
vacuum sphere converges monotone to :math:`\sim 1.78590` at
:math:`\mathcal{O}(h^{1.48})` (combined pole + outer-face first order;
increments :math:`2\mathrm{e}{-5}` at :math:`n_x = 160`, far below
engineering tolerance).  **This is why #233 needed an L∞ / per-cell
probe to surface** — no L2 or eigenvalue gate could see it.

The cylinder shares the same defect, masked
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The cylinder pole vs. **midpoint** is :math:`\mathcal{O}(h^2)`
(1.94 / 1.97 / 1.98) but vs. the **volume average** is
:math:`\mathcal{O}(h)` (0.99 / 0.99 / 1.00) — the SAME diamond
inconsistency, masked by the midpoint comparison: the cylinder's
:math:`r\,dr` (linear) weight puts the volume-centroid at
:math:`\tfrac23 h` while diamond's :math:`\tfrac12 A(h)` happens to
:math:`\approx` the midpoint :math:`A(h/2)`, so the midpoint comparison
the gate uses is *accidentally* :math:`\mathcal{O}(h^2)` for the
cylinder.  The cylinder pole is therefore **not** "clean
:math:`\mathcal{O}(h^2)`" — it is the same :math:`\mathcal{O}(h)`
volume-average defect, hidden by the comparison choice.  Cylinder global
L2 is also clean :math:`\mathcal{O}(h^{2.00})`.

The characterization gate (W2, #233)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Per the vv-principles "pin what is TRUE + protect the floor WITHOUT
calcifying the limitation" discipline, W2 ships a **characterization**
gate, not a fix gate
(``tests/sn/verification/mms/test_curvilinear_pole_cell_characterization.py``,
commit ``255eba4``):

* **Guarantee tests** (carry ``verifies("dd-curvilinear-scalar", ...)``,
  the :eq:`dd-curvilinear-scalar` cell-update label): global
  volume-weighted L2 is :math:`\mathcal{O}(h^2)` (``orders > 1.9``).
  The sphere is asserted under **both** references — midpoint AND the
  Hébert-3.430 shell-volume-average :eq:`sn-pole-cell-shell-average`,
  built from ``scipy.integrate.quad`` (a trusted-library integrator,
  structurally independent of the solver).  Agreement on the order
  across two structurally-different references proves the L2 order is
  REAL, not a midpoint artifact.
* **Characterization tests** (NO ``verifies`` — they pin a *limitation*,
  not a correctness claim): the pole L∞ order is **lower-bounded only**
  (:math:`> 0.8` — "at least first order, does not regress"), the pole
  is the L∞-dominant cell (fraction :math:`> 0.99`), and the interior
  is clean :math:`\mathcal{O}(h^2)` (:math:`> 1.8`).  **No upper bound**
  on the pole order, so a future LD / nodal scheme that lifts the pole
  to :math:`\mathcal{O}(h^2)` keeps the gate green
  (:math:`2.0 > 0.8`) — the characterization gate pins what is true and
  the regression floor without blocking a legitimate improvement
  (vv-principles anti-patterns #5 / #17).

Measured (sphere :math:`n_{\rm ord} = 16`, ladder
:math:`[40, 80, 160, 320]`): L2 midpoint orders :math:`2.01\times3`; L2
shell-average :math:`2.00\times3`; L∞ (pole) :math:`0.91 / 0.95 / 0.97`;
interior :math:`1.84 / 1.92 / 1.96`; pole fraction :math:`1.00` every
mesh.  Cylinder pole-vs-midpoint :math:`1.94 / 1.97 / 1.98`;
pole-vs-volavg :math:`0.99 / 0.99 / 1.00`.

.. _sn-p1-scattering-curvilinear:

Issue #9 — :math:`P_1` anisotropic SCATTERING in curvilinear (the two unrelated paths)
--------------------------------------------------------------------------------------

A persistent source of confusion in this cluster is that "anisotropic"
names **two structurally unrelated** things in a curvilinear SN solve.
Issue #9 is about the *second*; everything above (#229, the
:math:`\alpha`-dome, the :math:`\tau`-clamp) is about the *first*.

* **Path-(I) — geometric angular redistribution.**  The
  :math:`(1-\mu^2)/r\,\partial_\mu\psi` term (sphere) /
  :math:`\xi^2 B / r` (cylinder), threaded by the Morel--Montry Carlson
  :math:`\alpha`-dome.  This is **P0-only**; the "anisotropy" lives in
  the *angular-flux ansatz*, not in the scattering kernel.  The
  existing curvilinear aniso MMS cases
  (:math:`\psi = (A + \zeta B)/W`) exercise ONLY this path.  #229 is a
  Path-(I) test-design floor.
* **Path-(II) — Legendre SCATTERING moments.**  The
  :math:`P_1{+}` scattering source :math:`R\,\Lambda\,M`
  (``scattering.build_aniso_source``, ``scattering_order ≥ 1``),
  geometry-**agnostic**, wired identically for all geometries through
  ``_within_group_triple`` (the :math:`S` in the
  :math:`(L+C),\,S,\,B` triple carries :math:`P_1` when
  ``scattering_order = 1``).  No curvilinear test exercised Path-(II)
  before #9 — it is NEW coverage of an existing capability (NO
  ``orpheus/`` change; Path-(II) works as-is).

L0 — the operator-admits trick
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rather than derive a costly symbolic :math:`P_1`-source MMS, W4 feeds a
*known* anisotropic angular flux :math:`\psi_{{\rm ref},n} = (A + \zeta
B)/W` to the within-group :math:`S` operator at
``scattering_order = 1`` and isolates the :math:`P_1` contribution as
:math:`S_1.\mathrm{apply}(\psi) - S_0.\mathrm{apply}(\psi)`, asserted
**per ordinate** (NOT weight-summed — the :math:`\alpha`-dome
telescopes, vv anti-pattern #8) against a structurally-independent
hand-reference:

* **Sphere** (fully SH-table-INDEPENDENT — :math:`P_1(\mu) = \mu`
  directly):

  .. math::
     :label: sn-p1-sphere-hand-ref

     q_n^{P_1} \;=\; \frac{1}{W}\,3\,\mu_n\,\Sigma_{s1}\,\phi_1,
     \qquad
     \phi_1 \;=\; B(r)\,\frac{\sum_n w_n\,\mu_n^2}{W}.

* **Cylinder** (explicit :math:`Y_1^m` moment-sum, independent of the
  production :math:`R\,\Lambda\,M` einsum):

  .. math::
     :label: sn-p1-cylinder-hand-ref

     q_n^{P_1} \;=\; \frac{1}{W}\,3\,\Sigma_{s1}
                  \sum_m Y_1^m(\Omega_n)\,\phi_1^m,
     \qquad
     \phi_1^m \;=\; \sum_n w_n\,Y_1^m(\Omega_n)\,\psi_n.

Both agree at machine precision (rel. :math:`4.7\mathrm{e}{-15}` sphere /
:math:`5.6\mathrm{e}{-15}` cylinder), with a
``max|S₁−S₀| > 1e-6`` negative control (vv anti-pattern #11 — a dropped
:math:`P_1` makes :math:`S_1 - S_0 \equiv 0` and fails the non-zero
hand-ref match).  **1-group is legitimate here**: this is a
flux-shape / OPERATOR claim (the per-ordinate :math:`P_1` source reads
:math:`\phi_1`, flux-shape-dependent by construction), NOT an eigenvalue
claim — the 1-group-degeneracy rule applies only to *eigenvalue*
verification.

L1 — the directional eigenvalue
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Forward-peaked :math:`P_1` scattering (:math:`\bar\mu > 0`) **lowers**
:math:`k_{\rm eff}` versus :math:`P_0`.  The physics: positive
:math:`\bar\mu` preserves the forward direction, so in a finite,
vacuum-bounded sphere forward-preserved scattered neutrons are more
likely to cross the outer boundary → **enhanced leakage** → lower
:math:`k_{\rm eff}`.  This requires a **vacuum** outer BC (a reflective
sphere has no leakage → :math:`P_0 \equiv P_1`).  Validated robust:

* Homogeneous vacuum sphere :math:`R = 4 / 10 / 25`:
  :math:`\Delta = k_{\rm eff}^{P_1} - k_{\rm eff}^{P_0} =
  -3.76\mathrm{e}{-3} / -1.32\mathrm{e}{-3} / -2.88\mathrm{e}{-4}` — sign
  always negative, :math:`|\Delta|` **grows as the sphere shrinks**
  (the leakage-monotone signature, the structural negative control a
  sign-flipped or absorption-mimicking :math:`P_1` would violate).
* Heterogeneous fuel-core(:math:`r < 5`)+moderator-shell vacuum sphere
  :math:`R = 10`: :math:`\Delta = -1.40\mathrm{e}{-2}` (:math:`140\times`
  the :math:`1\mathrm{e}{-3}` detection bar), with materials
  ``get_mixture("A","2g")`` (the only fissile 2-group mixture;
  asymmetric downscatter-only P0 avoids the 1-group degeneracy) and
  ``get_mixture("C","2g")``.

Two L1 rows pin this: a heterogeneous-sphere
:math:`k_{\rm eff}^{P_1} < k_{\rm eff}^{P_0}` AND
:math:`1\mathrm{e}{-3} < (P_0 - P_1) < 5\mathrm{e}{-2}`; and a
leakage-monotone control
:math:`(P_0 - P_1)|_{R=4} > (P_0 - P_1)|_{R=25} > 0` (the mechanism
pin).  These are the **first curvilinear exercise** of the
geometry-agnostic ``pn-scatter`` / ``flux-moments`` labels (prior tests
were 2-D Cartesian only).  L0 lands in
``tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py``,
L1 in ``tests/sn/eigenvalue/test_keff_curvilinear.py::TestSphereP1DirectionalEigenvalue``
(commit ``d5878e9``).  L2 is deferred (subsumed by L0+L1; a
:math:`P_1`-convergence L2 needs the :math:`\sigma_{s1}`-MMS source and
rides the same #229 floor).

Infrastructure retained
-----------------------

Per the aggressive-retirement exception, the program deletes no correct
machinery:

.. list-table:: Curvilinear aniso program — primitives status
   :header-rows: 1
   :widths: 34 18 48

   * - Primitive
     - Status
     - Why kept / what changed
   * - Spherical :math:`\tau_m^{\rm raw}` (unclamped)
     - **production**
     - W1: the unique exact-on-linear weight; single-sourced in
       :func:`~orpheus.geometry.reduced_operator.spherical_streaming`,
       inherited by SI sweep + Krylov matvec.
   * - Cylindrical :math:`\tau_m` clamp
     - **production**
     - Retained — the :math:`\tau^{\rm raw} = 0` structural
       :math:`\div 0` block; removing it needs a 2-D
       :math:`(\eta,\varphi)` closure (#229).
   * - Pole-cell characterization gate
     - **regression net**
     - Pins the inherent :math:`\mathcal{O}(h)` pole limitation
       (lower-bounded, not calcified) + the global :math:`\mathcal{O}(h^2)`
       guarantee under two independent references.
   * - Shell-volume-average reference :eq:`sn-pole-cell-shell-average`
     - **oracle**
     - The Hébert-3.430 finite-volume unknown; the principled MMS
       reference that removes the :math:`\sim 75\,\%` comparison
       artifact.  Built from ``scipy.integrate.quad``.

Open research paths (research-tag, not production-blocking)
-----------------------------------------------------------

#. **Higher-order central-cell spatial scheme** (lifts the #233 pole
   :math:`\mathcal{O}(h)`).  The canon ([Hebert2009]_ §3.9.4, Stacey
   §9.9) provides no drop-in :math:`\mathcal{O}(h^2)` central-cell
   diamond closure; the documented route is a non-local higher-order
   *spatial* scheme — linear-discontinuous (Issue #6), step-
   characteristic, or the Green's-function nodal method of
   [WuXieFischer1999]_ (NSE 133, "very high precision on coarse meshes
   relative to standard fine-mesh DD").  Likely diagnostic probe: the
   pole-cell per-cell rate under the shell-average reference, holding
   quadrature fixed.
#. **2-D** :math:`(\eta, \varphi)` **cylinder angular closure** (lifts
   the #229 cylinder floor).  The 1-D :math:`\eta`-thread cannot
   represent the duplicate-azimuthal-:math:`\eta` variation of product /
   level-symmetric quadratures; a genuine 2-D angular closure (or a
   Gauss-type azimuthal quadrature with distinct :math:`\eta` values,
   GitHub Issue #1) is required.  Likely probe: the floor-scaling table
   above with the azimuthal quadrature replaced by a distinct-:math:`\eta`
   set.

Session trail (V&V audit trail)
-------------------------------

* **Commits** (branch ``fix/curvilinear-aniso-pole-and-clamp``,
  2026-06-13): ``b2d8a6d`` (W1 sphere unclamp), ``255eba4`` (W2 #233
  pole-cell characterization gate), ``d5878e9`` (W4 #9 :math:`P_1`
  scattering coverage), ``679a1e6`` (W3 #229 gate retune).
* **Diagnostics**: the W1–W4 ``diag_01..31`` decomposition probes (the
  decisive ones: the
  :math:`E_{\rm test} = E_{\rm artifact}(\text{midpoint} -
  \text{volavg}) + E_{\rm true}(\text{solver} - \text{volavg})`
  decomposition; the discrete-balance residual fed exact fields; the
  faithful production-sweep monkeypatch with a :math:`\beta = \tfrac12`
  identity guard).
* **Literature**: [Hebert2009]_ §3.9.4, Stacey §9.9 (plain diamond at
  the central cell, no special closure), [BaileyMorelChang2010]_ Eq. 43
  (the exact-on-linear weight), [WuXieFischer1999]_ (the nodal route to
  :math:`\mathcal{O}(h^2)` at the origin).
* **vv catalogue**: ``error_catalog.md`` — ERR-059 (the pole-cell
  inherent limitation) + the :math:`\tau`-clamp mis-citation finding +
  the ERR-026 surviving-manifestation note.
* **Issues**: #229 (cylinder floor + sphere gate retune), #9
  (:math:`P_1` curvilinear scattering), #233 (pole-cell, stays OPEN to
  track the future higher-order scheme).

.. note:: **vv-status (eq-labels added by this section).**  The labels
   :eq:`sn-tau-mm-raw`, :eq:`sn-pole-cell-shell-average`,
   :eq:`sn-p1-sphere-hand-ref`, and :eq:`sn-p1-cylinder-hand-ref` are
   *structural / representational* identities (the literature-
   transcribed M-M weight; the Hébert-3.430 finite-volume unknown; the
   structurally-independent :math:`P_1` hand-references).  They are NOT
   solver claims.  The :math:`\tau`-clamp / pole-cell / :math:`P_1`
   *verifiable* content is the W1 clamp-silence + positivity gates, the
   W2 ``verifies("dd-curvilinear-scalar")`` guarantee tests, and the W4
   ``verifies("pn-scatter","flux-moments")`` per-ordinate operator-
   admission gates named above — so these eq-labels are ``documented``.


Krylov inner solver
===================

Instead of sweep-based source iteration, the within-group transport
problem can be solved directly using a Krylov method.  Wave E Round 2
(Issue #164) replaces the legacy BiCGSTAB FD-operator path with GMRES
on :meth:`InvertibleOperator.apply` (the symmetric closure carried
by the operator algebra) with the sweep as a left preconditioner —
the SAILOR / Larsen-Adams preconditioned-Krylov framework
([AdamsLarsen2002]_ §III).

Operator equation
-----------------

The streaming-collision operator :math:`L` is formed explicitly via
:class:`~orpheus.sn.operator.InvertibleOperator`:

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
:meth:`InvertibleOperator.apply` reuses was designed for the
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

**Stage 2 --- Solver resolution via the BC realizer.**
:class:`SNMesh` owns a class-level
:attr:`~SNMesh.BOUNDARY_OPERATOR_REGISTRY` mapping kind strings to
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` **subclasses**
(post Wave 8 of the trace-law refactor in
``.claude/plans/transient-giggling-cake.md``)::

    BOUNDARY_OPERATOR_REGISTRY = {
        "vacuum":     VacuumInflow,
        "reflective": ReflectiveBoundary,
    }

The registry values are the law classes themselves, not factory
functions. The pre-refactor ``_sn_vacuum_boundary_operator`` /
``_sn_reflective_boundary_operator`` factories were retired; their
job is now done by :meth:`SNMesh._resolve_one`, which dispatches
through :class:`~orpheus.sn.boundary_realizer.SNBoundaryRealizer`
**uniformly** for every supported mesh (1-D Cartesian, 1-D
spherical, 1-D cylindrical, 2-D Cartesian) — see
:ref:`bc-sn-resolution-table` below. Issue #188 (curvilinear support
in the trace space — then named ``InflowTraceSpace``, now the unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace`) and Issue
#176 (drop 2-arg ``apply`` + simplify shim) collapsed the pre-cleanup
Cartesian-vs-curvilinear bypass into a single realizer-routed
path; details at :ref:`bc-curvilinear-realizer-unification`.

During ``SNMesh.__init__``, each face's :class:`~geometry.mesh.BC`
is looked up in the registry.  If the kind is not found, a
``ValueError`` lists the supported kinds.  For curvilinear
geometries (spherical, cylindrical), only ``"reflective"`` and
``"vacuum"`` are currently supported --- requesting any other kind
on a curvilinear mesh raises ``ValueError``.

The two-arg legacy interface ``bc.apply(angular_flux_outgoing,
quadrature)`` is **retired entirely** post Issue #186 / B3 + β2
(2026-05-11). Concrete BC :meth:`apply` methods no longer exist;
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` is a **pure
descriptor** with no callable interface. Rank-N composition uses
the descriptor-tree algebra
(:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled`) with
:func:`~orpheus.sn.boundary_realize.realize_recursively` as the
sole descriptor→operator type transformer. See
:ref:`bc-trace-law-descriptor-model` for the design rationale.

The resolved BCs at ``sn_mesh.bc["xmin"]`` etc. expose the uniform
1-arg contract through the
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim — internal to the package, not in
:attr:`orpheus.geometry.boundary.__all__`. Post Issue #186 the
shim is a **strict 1-arg passthrough** (extra args raise
:class:`TypeError`); it carries the originating ``BC.kind`` string
so the ``sn_mesh.bc["xmin"] == "vacuum"`` diagnostic
comparison continues to evaluate True iff the underlying law is
:class:`VacuumInflow`. (C4 / #220 re-keyed this surface from the
per-attribute ``sn_mesh.bc_left`` to the face-name-keyed
:attr:`SNMesh.bc` dict — see :ref:`bc-face-name-carve`.)
See :ref:`bc-tensor-decompositions` below
for the operator-algebra view and
:ref:`theory-boundary-conditions` for the full trace-law /
realizer architecture.

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

.. list-table:: Implemented :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` primitives (Wave 7 vocabulary)
   :widths: 20 30 25 25
   :header-rows: 1

   * - Class
     - :math:`G_\alpha`
     - :math:`A_\alpha`
     - Rank / wired into ``solve_sn``
   * - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - :math:`0` (no operator)
     - 0
     - 0 / yes
   * - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - permutation under reflection axis
     - albedo (1 = perfect)
     - 1 / yes
   * - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - cosine-weighted hemispheric average
     - albedo
     - 1 / no (Wave C)
   * - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - spatial pushforward (caller-supplied)
     - 1
     - 1 / no (Wave C/D)
   * - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - identity in angle
     - albedo
     - 1 / no (building block)
   * - :class:`~orpheus.geometry.boundary.PrescribedInflow`
     - 0
     - 0
     - 0 with :math:`q \neq 0` / no (rank-0 source-only affine BC)

The pre-Wave-7 names ``VacuumBoundaryOperator`` /
``SpecularBoundaryOperator`` / ``WhiteBoundaryOperator`` /
``PeriodicBoundaryOperator`` / ``AlbedoBoundaryOperator`` were
**retired in Wave O step O.4a.1**; their canonical successors
(:class:`~orpheus.geometry.boundary.VacuumInflow` /
:class:`~orpheus.geometry.boundary.ReflectiveBoundary` /
:class:`~orpheus.geometry.boundary.WhiteBoundary` /
:class:`~orpheus.geometry.boundary.PeriodicBoundary` /
:class:`~orpheus.geometry.boundary.AlbedoBoundary`) are the sole
live names in :mod:`orpheus.geometry.boundary`.
``MixedBoundaryOperator`` was **retired in Wave 11**; rank-N
(Marshak, partial-current) boundaries are now expressed via the
**descriptor-tree algebra**
(:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled`) on the unrealised
laws:

.. code-block:: python

   tree = 0.3 * spec + 0.7 * white            # LawSum of LawScaled
   op = realize_recursively(tree, method_space)  # OperatorSum of ScaledOperator
   psi_in = op.apply(psi_out)

See :ref:`bc-rank-n-algebra` for the closed algebra and the
:ref:`bc-realize-recursively` walker that lowers a descriptor
tree to a Wave-0 operator tree.

The abstract base :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
is a **pure descriptor** post Issue #186 / B3 + β2 (2026-05-11)
— it has **no** :meth:`apply` method. The :class:`LinearOperatorMixin`
inheritance that historically supplied ``apply`` was removed; the
concrete laws likewise carry no ``apply`` / ``apply_transpose``
methods. The §16A.3 three-layer architecture (descriptor /
realizer / operator) is now enforced by the type system: a static
type checker rejects ``law.apply(...)`` at the linter level
without running the program. The full retrospective on the
predecessor Option A and β1 forms (and why each was rejected) is
at :ref:`bc-trace-law-descriptor-model`.

The tensor framing pays off architecturally because partial-current
Marshak boundaries (:math:`R = c_1 \, G_{\rm refl} + c_2 \, G_{\rm
diff}`, Bell & Glasstone 1970 §1.5) and multi-region interface
couplings are all instances of the same algebra: pick the geometric
operators, pick the amplitudes, sum. New BCs are one
:class:`BoundaryTraceLaw` subclass + one
``BOUNDARY_OPERATOR_REGISTRY`` entry away — no sweep edits per BC.

.. _bc-sn-resolution-table:

SN BC resolution table (post-Issue-#188+#176 wiring)
-----------------------------------------------------

The :meth:`SNMesh._resolve_one` dispatch is summarized below.
Each row maps the user-facing :class:`~orpheus.geometry.mesh.BC`
kind string to (a) the resolved :class:`BoundaryTraceLaw`
subclass, (b) the :class:`SNBoundaryRealizer.realize` output
operator, and (c) the ``creates_sweep_cycle`` flag used by §15A.2
sweep-cycle detection (see :ref:`bc-sweep-cycle`). The realizer
dispatch is **uniform** across every supported mesh — 1-D
Cartesian / spherical / cylindrical and 2-D Cartesian — post
Issue #188 + #176.

.. list-table:: BC.kind → law class → realized SN operator
   :header-rows: 1
   :widths: 16 22 36 12 14

   * - ``BC.kind``
     - Law class
     - Realized SN operator
     - α
     - ``creates_sweep_cycle``
   * - ``"vacuum"``
     - :class:`~orpheus.geometry.boundary.VacuumInflow`
     - :class:`~orpheus.numerics.operator.IncomingOrdinateMaskTensor`
       (per-face inflow indices)
     - —
     - ``False``
   * - ``"reflective"``
     - :class:`~orpheus.geometry.boundary.ReflectiveBoundary`
     - :class:`~orpheus.numerics.operator.PermutationOperator`
       (``quadrature.reflection_index(axis)``)
     - 1 (fast path)
     - **``True``**
   * - ``"reflective"``
     -
     - ``α * PermutationOperator``
       (:class:`~orpheus.numerics.operator.ScaledOperator`)
     - α ≠ 1
     - **``True``**
   * - ``"white"``
     - :class:`~orpheus.geometry.boundary.WhiteBoundary`
     - :class:`~orpheus.sn.angular_operator.AngularAverageOperator`
     - 1 (fast path)
     - ``False``
   * - ``"white"``
     -
     - ``α * AngularAverageOperator``
     - α ≠ 1
     - ``False``
   * - ``"periodic"``
     - :class:`~orpheus.geometry.boundary.PeriodicBoundary`
     - :class:`~orpheus.numerics.operator.PeriodicWrapOperator`
     - 1
     - **``True``**
   * - ``"albedo"``
     - :class:`~orpheus.geometry.boundary.AlbedoBoundary`
     - :class:`~orpheus.numerics.operator.ZeroOperator`
     - 0
     - ``False``
   * - ``"albedo"``
     -
     - :class:`~orpheus.numerics.operator.IdentityOperator`
     - 1
     - ``False``
   * - ``"albedo"``
     -
     - ``α * IdentityOperator``
     - α ∉ {0, 1}
     - ``False``
   * - ``"prescribed_inflow"``
     - :class:`~orpheus.geometry.boundary.PrescribedInflow`
     - :class:`~orpheus.sn.angular_operator.IncomingSourceOperator`
       (source.evaluate; ignores outgoing flux)
     - —
     - ``False``

The :meth:`SNMesh._resolve_one` dispatch constructs the resolved
operator via :meth:`SNBoundaryRealizer.realize(law, method_space)`
where the ``method_space`` is built by
:meth:`SNMethodSpace.for_face` carrying the precomputed unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace` (built
once at :class:`SNMesh` construction for every supported mesh).
The reflective branch derives its reflection axis from the face's
own :class:`~orpheus.sn.axis.FaceLabel` —
``AXIS_NAMES[label.axis_index]`` — so the partner is correct at any
dimension by construction (C4 / #220; see
:ref:`bc-face-name-latent-d3-bug`). The
:class:`~orpheus.geometry.boundary._bound_compat._BoundBoundaryOperator`
shim wraps the result with a ``kind`` tag for the
``sn_mesh.bc["xmin"] == "vacuum"`` string-equality surface.

.. note::

   **The 1-D y-face placeholders were retired in C4 / #220.** Pre-C4,
   a slab :class:`SNMesh` carried a pair of realized no-op
   ``ReflectiveBoundary(axis="y")`` operators at ``bc_ymin`` /
   ``bc_ymax`` so cross-dimensional code could read them without
   coord-system gating — but **no production code ever read them**
   (a 1-D mesh's ``trace.layout.faces`` is ``("xmin", "xmax")``).
   C4 makes them unrepresentable: a slab has no y-axis in its
   :attr:`~SNMesh.axes` tuple, so
   :func:`~orpheus.sn.axis.face_labels` emits no y-label and
   :attr:`SNMesh.bc` has no y-entry — ``slab.bc["ymin"]`` is a
   :class:`KeyError`, not a no-op. See
   :ref:`bc-face-name-carve-what-retired` for the full retirement
   record (the pre-C4 "why the placeholders were once safe"
   rationale is preserved there).

**Pre-cleanup history.** Before Issue #188 + #176 (closed
2026-05-11), curvilinear ``Mesh1D`` bypassed the realizer because the
trace factory (then ``InflowTraceSpace.from_mesh_and_quadrature``, since
C5.3 the geometry-blind
:meth:`TraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.trace_space.TraceSpace.from_quadrature_and_layout>`)
raised :class:`NotImplementedError` on those coord systems; the
``_BoundBoundaryOperator`` shim carried a dual mode where the
``quadrature=`` kwarg, when non-``None``, bound an
:class:`AngularQuadrature` and forwarded ``inner.apply(psi,
bound_quad)`` to the legacy 2-arg :class:`BoundaryTraceLaw` body.
The bypass and dual-mode are gone; details and the algebraic
sequence ("Issue #188 unblocks Issue #176") at
:ref:`bc-curvilinear-realizer-unification`.

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
   into a spatial-energy moment field carrying the principled
   ``(ng, nx, ny)`` layout (see :ref:`theory-sn-index-convention`;
   the codepath presents the full moment field as
   ``(L+1, 2L+1, ng, nx, ny)`` with energy leading the spatial axes).

3. **Reconstruct per-ordinate source**: for each Legendre order
   :math:`\ell \geq 1` (the :math:`\ell = 0` term is handled by
   :meth:`SNSolver._add_scattering_source`) and each :math:`m`, the
   scattered moment ``moment @ sig_s_l[l]`` is multiplied by
   :math:`(2\ell+1) Y_\ell^m(\hat{\Omega}_n)` and accumulated into
   ``Q_aniso[n, :, :, :]``.

4. The resulting ``Q_aniso`` array of shape ``(N, ng, nx, ny)`` is
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

.. note::

   Because the scattering source :eq:`pn-scatter` depends on the
   angular flux **only** through its moments :math:`f_{\ell,g}^m`
   :eq:`flux-moments`, the within-group source iteration's *fixed point
   lives in moment space*: the persistent iterate need not carry all
   :math:`N` ordinates. The 2-D Cartesian SI iterate is therefore held
   as the moment tensor (:math:`N \to (L{+}1)(2L{+}1)`, "angular
   windowing"), with the :math:`\ell\ge 1` reconstruction
   :math:`R\,\Lambda` shared between the windowed and full-angular
   paths. The moments are accumulated **in-sweep** per anti-diagonal
   (:math:`\phi_\ell^m \mathrel{+}= \sum_n w_n Y_\ell^m \psi_n`), so the
   full per-ordinate field is never materialized in the windowed iterate
   (a 3.06× peak-memory win). See :ref:`sn-angular-windowing` in
   :doc:`operator_algebra` for the derivation, the geometry restriction,
   and the bit-identity / principled-equivalence story, and
   :ref:`sn-angular-windowing-in-sweep-accumulation` for the in-sweep
   accumulation.

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

.. note::

   The Pℓ Galerkin reconstruction has named primitives in
   :mod:`orpheus.numerics.projection` (Wave 0 of the SN
   performance plan):
   :class:`~orpheus.numerics.projection.HarmonicMomentProjection`
   (the :math:`\Pi = Y^* W` operator on the angular axis) and
   :class:`~orpheus.numerics.projection.HarmonicMomentReconstruction`
   (the addition-theorem reconstruction with the
   :math:`(2\ell+1)` factor). The full-space projector is the
   tensor product :math:`\Pi \otimes I_x \otimes I_y \otimes I_g`,
   built via the ``&`` dunder of
   :class:`~orpheus.numerics.operator.TensorProductOperator`. See
   :ref:`galerkin-projection` for the cross-method consumer table
   (PN solver, energy condensation, MC adjoint moments) and
   :ref:`spherical-harmonics` for the convention and addition
   theorem. The :math:`Y_\ell^m` evaluator
   :func:`~orpheus.numerics.spherical_harmonics.evaluate_real_sh`
   is the canonical generic infrastructure consumed here; the
   SN-side ``orpheus.sn.quadrature._build_spherical_harmonics``
   is now a thin re-export.

   **Wave 1 (commit ff454f2)**:
   :meth:`~orpheus.sn.scattering.ScatteringOperator.build_aniso_source`
   is now the literal §9 line 1230 operator-algebra composition

   .. math::

      Q^{\rm aniso}_n(\vec r) \;=\; R\,\Lambda\,M\,\psi

   where :math:`\Lambda` is :class:`~orpheus.sn.scattering.LegendreMomentScattering`
   --- the per-ℓ block-diagonal scattering on moment space (the §15.2
   sum-of-tensor-products form
   :math:`\Lambda = \sum_\ell P_\ell \otimes \Sigma_{s,\ell}`).
   The previous ``for n in range(N)`` Python loop over ordinates is
   gone *by construction*: each constituent's :meth:`apply` carries
   the ordinate iteration internally via :func:`numpy.einsum`, not
   via a Python loop. Total flop count is unchanged; the iteration
   is structural rather than buried in a triple-nested loop. The
   refactor is gated by the
   ``slab_2g_p1_aniso_dd_n20`` regression snapshot (rtol=1e-12,
   atol=1e-13) and the full
   :file:`tests/sn/test_mms_aniso.py` Pℓ MMS convergence suite.

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
:class:`~orpheus.geometry.boundary.BoundaryTraceLaw` infrastructure
(Wave B Issue 7), so :func:`solution_to_angular_flux*` and the
matvec helpers now dispatch boundary fills via the realiser-routed
1-arg :meth:`apply` on the resolved
:class:`~orpheus.numerics.operator.LinearOperator` — vacuum,
reflective, white, albedo, periodic, and mixed BCs are honoured
uniformly. (Post Issue #186 / B3 + β2, the law itself is a pure
descriptor; the SN realiser produces the callable. See
:ref:`bc-trace-law-descriptor-model`.)


.. _sn-streaming-operator:

InvertibleOperator: the streaming-collision operator algebra
==============================================================

Wave D Round 3 (Issue #160) installs
:class:`~orpheus.sn.operator.InvertibleOperator` as the unified
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

:class:`InvertibleOperator` advertises ``{"apply", "solve",
"apply_transpose"}`` — every member of the Wave A capability set:

* :meth:`~orpheus.sn.operator.InvertibleOperator.apply` —
  matrix-free forward action :math:`L\,\psi`.  Reuses the
  symmetric closure math from the existing
  :func:`~orpheus.sn.operator.transport_operator_matvec`,
  :func:`~orpheus.sn.operator.transport_operator_matvec_spherical`,
  :func:`~orpheus.sn.operator.transport_operator_matvec_cylindrical`
  functions (the historical BiCGSTAB FD operator).  The math is
  **extracted verbatim**; the new class is a thin Protocol wrapper.

* :meth:`~orpheus.sn.operator.InvertibleOperator.solve` —
  inverse action :math:`L^{-1}\,q` via the Wave D Round 2 unified
  sweep (:func:`~orpheus.sn.loss_representation.transport_sweep`).  Bit-identical
  to a direct :func:`transport_sweep` call on the same arguments.

* :meth:`~orpheus.sn.operator.InvertibleOperator.apply_transpose` —
  adjoint action :math:`L^*\,\psi`.  Implemented via the explicit
  transpose of the dense matrix assembled by probing
  :meth:`apply` with each unit basis vector.  The construction is
  exact by linear algebra and gates the reciprocity invariant
  :math:`\langle L\,\psi,\,\varphi\rangle = \langle\psi,\,
  L^*\,\varphi\rangle` (see "Reciprocity" subsection below).

Apply and solve use **different** closures by design
----------------------------------------------------

This is the load-bearing architectural fact about
:class:`InvertibleOperator`, and the reason the operator's
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
  (:func:`~orpheus.sn.loss_representation.transport_sweep`, dispatching through
  the Wave D Round 2 :class:`~orpheus.sn.spatial.scheme.DiscretizationScheme`
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
   :class:`InvertibleOperator` shipped only :meth:`apply`, the
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

Post-Depth-B (2026-05), :meth:`apply_transpose` on
:class:`InvertibleOperator` is inherited from
:class:`~orpheus.numerics.operator.OperatorSum`'s adjoint-propagation
closure law: each leaf's ``.H`` adjoint composes via the sum/
difference algebra, so the composite transpose is built analytically
(no dense-matrix probing).  The pre-Depth-B implementation assembled
a dense matrix by probing :meth:`apply` with unit basis vectors and
returned the explicit transpose; that path retired with the bundled
``SNStreamingOperator`` class in D-K.

Reciprocity gating today: the foundation linearity gate
:func:`tests.sn.test_streaming_operator.test_apply_is_linear` catches
non-linearity in :meth:`apply`, and the Resolution A bit-exact
decomposition gate
:file:`tests/sn/test_streaming_operator_decomposition.py` catches
:math:`(L+C).{\rm apply} \neq M(\psi;\sigma_t)` drift.  The full
analytical-adjoint Gate 1.3 round-off pin lands with **Wave O**
(`Issue #208 <https://github.com/deOliveira-R/ORPHEUS/issues/208>`_,
the BulkOperator / FullOperator / BoundaryOperator
adjoint algebra); see ``test_phase_c_gates.py`` for the current
xfail-strict placeholder.

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
:func:`transport_sweep`'s contract, in the principled
``(ng, nx, ny)`` / ``(N, ng, nx, ny)`` layout
(see :ref:`theory-sn-index-convention`):

* Source ``Q`` shape ``(ng, nx, ny)``.
* Persistent boundary-flux dict ``psi_bc`` carrying state between
  sweeps.
* Optional anisotropic source ``Q_aniso`` shape
  ``(N, ng, nx, ny)`` for P\ :sub:`ℓ` (:math:`\ell\ge 1`)
  scattering.

The shape mismatch reflects the FD-matvec packed-vector internal
convention (the packed vector is laid out group-major in Fortran
order; see :ref:`theory-sn-index-convention`'s "What stayed
deliberately legacy" subsection).  Wave E will normalise these via
a single Krylov-on-apply path; until then, calling :meth:`apply` and
:meth:`solve` requires the caller to be aware of the layout
difference.  The principled-storage flip
(:ref:`theory-sn-index-convention`) does NOT touch the packed-vector
layout --- that work is deferred to PR-INDEX-7.

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
  :class:`InvertibleOperator`: the inner solve becomes Krylov on
  :meth:`apply` with sweep as preconditioner, which removes the WDD
  asymmetric closure from the converged-solution path for the
  reflective-BC eigenvalue case.  Wave E Round 3 will extend the
  equation-map layout to the vacuum-BC path that closes ERR-026
  for fixed-source curvilinear MMS.
* **Forward**: when production reciprocity becomes performance-
  critical, an :math:`O(n)` analytic-adjoint matvec replaces the
  dense-transpose fallback; the new path is gated by the
  reciprocity tests in :file:`tests/sn/test_streaming_operator.py`
  (post-D-K successor) and :file:`tests/sn/test_phase_c_gates.py`
  Gate 1.3 (xfail pending Wave O).

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
geometric at rate :math:`c` (Lewis & Miller §4.4).  The formal
derivation, the iterations-to-tolerance estimate, and the Phase 3
boundary Gauss-Seidel rate recovery are in
:ref:`si-within-group-splitting`; the labelled spectral-rate
identity is :eq:`si-spectral-rate`.

The convergence test is the relative L2 residual on the iterate — the
same metric :meth:`SNSolver._solve_source_iteration` uses, since that
within-group inner now consumes this primitive directly (Phase 1 R1,
via the ``_within_group_triple`` single-source-of-truth helper):

.. math::

    {\rm res}_n \;=\; \frac{\|\psi_n - \psi_{n-1}\|_2}
                            {\max(\|\psi_n\|_2,\,10^{-30})}

with the iteration breaking when :math:`{\rm res}_n < {\rm tol}`.

.. _si-within-group-splitting:

The within-group SI splitting and its spectral rate
----------------------------------------------------

This subsection derives the source-iteration spectral radius
:math:`\rho_J = c` from the within-group operator splitting, gives
the iterations-to-tolerance estimate the rate implies, and then
documents the Phase 3 *boundary Gauss-Seidel* rate recovery — what
it accelerates, what it does **not** (the load-bearing honest
scope), the failed σ\ :sub:`r`-fold that motivated it (GitHub
issue #215), and the diagonal-cubature shared-face correctness rule
(ERR-056).  The polymorphic schedule lives in
:mod:`orpheus.sn.sweep_schedule`; the SI resolvent that consumes it
is :class:`~orpheus.sn.solver._GaussSeidelResolvent`; the public
entry is :func:`~orpheus.sn.solver.solve_sn_fixed_source` via the
``inner_schedule`` keyword.

.. admonition:: Key Facts (SI rate)
   :class: tip

   * The within-group SI iteration matrix is
     :math:`(L+C)^{-1}(S+B)`; its spectral radius is the scattering
     ratio :math:`\rho_J = c = \max_g \Sigma_{s,g}/\Sigma_{t,g}`
     (:eq:`si-spectral-rate`).  Iterations to relative tolerance
     :math:`\varepsilon`: :math:`n_{\rm Jacobi} \approx
     \log\varepsilon / \log c`.
   * **Boundary Gauss-Seidel** (Phase 3, ``inner_schedule=
     "gauss_seidel"``, default) folds **only** the boundary
     reflection :math:`B` into the resolvent
     (:math:`(L+C-B_{\rm lower})^{-1}` forward substitution).  It
     accelerates the *boundary-layer transient* by a measured,
     regime-independent **~0.86–0.92×**, NOT the dominant flat
     scattering :math:`c`-mode.  **This is NOT the textbook
     scattering-G-S** :math:`c^2`-halving.
   * The dominant within-group scattering rate is recovered ONLY by
     **Krylov** (already production — rate-optimal,
     splitting-invariant) or by **consistent DSA** (future, GitHub
     issue #2).  The scattering :math:`c`-mode **cannot** be folded
     into the directional sweep (the σ\ :sub:`r`-fold trap, issue
     #215).
   * The converged fixed point is **invariant** under the schedule
     (any consistent splitting of :math:`(L+C-S-B)\psi=q` shares
     :math:`\psi^\ast`); only the SI spectral rate changes.  Krylov
     ignores the schedule entirely.
   * 1-D meshes always fall back to Jacobi (the 1-D scan is not a
     wavefront; the regime is scattering-dominated — boundary G-S is
     a no-op).

The four-operator within-group equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Within a single energy group, the steady transport equation factors
into four operators acting on the angular flux :math:`\psi`:

.. math::
   :label: si-within-group-operator-eq

   (L + C - S - B)\,\psi \;=\; q,

where

* :math:`L = \hat\Omega\cdot\nabla` is **streaming** (the sweep's
  spatial derivative — see :ref:`sn-streaming-operator`);
* :math:`C = \Sigma_t\,\mathbb{I}` is the **collision** (total
  removal) operator, diagonal in angle;
* :math:`S = \Sigma_{s,0}\,P_{\rm iso}` is the **within-group
  scattering** gain — it couples back through the scalar flux
  :math:`\phi = \int\psi\,d\Omega`, so :math:`P_{\rm iso}\psi =
  \phi/\!\sum_n\mathcal{W}_n` is the isotropic-projection (rank-1
  in angle) operator (the convention used by
  :class:`~orpheus.sn.scattering.ScatteringOperator`; higher
  Legendre orders add the :math:`P_\ell` blocks);
* :math:`B` is the **boundary reflection** gain — trace-only,
  realised by :class:`~orpheus.sn.boundary_operator.SNBoundaryOperator`,
  delivering :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}` on
  specular faces (see :ref:`bc-extraction` in
  :doc:`operator_algebra`);
* :math:`q` is the external/fission within-group source
  (:eq:`phase-f-q-1d-decomposition`).

Source iteration is a **splitting** of :eq:`si-within-group-operator-eq`:
the **streaming-collision** part :math:`(L+C)` is kept on the LHS
(inverted exactly by **one WDD sweep** — a triangular
forward-substitution, since the sweep visits cells in causal order),
while the **scattering** :math:`S` and the **boundary reflection**
:math:`B` are *lagged* on the RHS, evaluated from the previous
iterate :math:`\psi_n`:

.. math::
   :label: si-jacobi-fixed-point

   \psi_{n+1} \;=\; (L+C)^{-1}\bigl(S\,\psi_n
                    \;+\; B\,\psi_n \;+\; q\bigr).

The iteration matrix is therefore :math:`M = (L+C)^{-1}(S+B)`, and
the asymptotic convergence rate is :math:`\rho(M)`.

Spectral radius :math:`\rho_J = c`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an infinite homogeneous medium with isotropic scattering, the
boundary term vanishes (:math:`B=0`) and the streaming derivative
:math:`L` drops in the flat-flux Fourier mode :math:`k\to0`.  The
spatial operator :math:`(L+C)^{-1}` then reduces to multiplication by
:math:`1/\Sigma_t`, and the isotropic-scattering gain :math:`S`
contributes :math:`\Sigma_{s,0}` per collision.  The dominant
eigenvalue of :math:`(L+C)^{-1}S` is thus the **scattering ratio**:

.. math::
   :label: si-spectral-rate

   \rho_J \;=\; \rho\!\bigl((L+C)^{-1}(S+B)\bigr) \;=\;
   c \;\equiv\; \max_g \frac{\Sigma_{s,g}}{\Sigma_{t,g}},
   \qquad
   n_{\rm Jacobi} \;\approx\; \frac{\log\varepsilon}{\log c}

(the Fourier / mode analysis of Lewis & Miller §4.4, Adams & Larsen
2002 §II).  The right-hand identity gives the iterations
:math:`n_{\rm Jacobi}` needed to drive the relative residual to a
tolerance :math:`\varepsilon`: each iteration multiplies the error
by :math:`c`, so :math:`c^{\,n} = \varepsilon` solves to
:math:`n = \log\varepsilon/\log c`.  Because :math:`c\to1` for a
nearly-pure scatterer, source iteration becomes arbitrarily slow as
:math:`c\to1` — the canonical motivation for acceleration (DSA,
Krylov).

.. note::

   The :math:`c` in :eq:`si-spectral-rate` is the **within-group
   scattering ratio** :math:`\Sigma_{s,0}^{g\to g}/\Sigma_{t,g}` that
   governs the *within-group* fixed point.  The
   :meth:`Mixture.scattering_ratio <orpheus.data.macro_xs.mixture.Mixture.scattering_ratio>`
   property exposes the slightly larger **Case–Zweifel** secondaries-
   per-collision parameter :math:`c_g = (\Sigma_{s,g} +
   \nu\Sigma_{f,g})/\Sigma_{t,g}` (it folds in fission emission for a
   multiplying medium).  The L1 rate anchor
   :func:`tests.sn.verification.analytical.test_si_convergence_rate.test_si_jacobi_rate_matches_scattering_ratio`
   pins :math:`n_{\rm Jacobi}` against ``log(tol)/log(c_max)`` using
   the Case–Zweifel form and accepts a 0.6–1.2 band: the measured
   B-2g slab count was **655** against a predicted
   :math:`\log(10^{-8})/\log(0.975) = 728` (ratio **0.90** — the
   gap is the finite-slab leakage that lowers the effective rate
   below the infinite-medium :math:`c`, plus the multigroup
   coupling).  This is the structurally-independent target the
   recovery improves upon, NOT another ORPHEUS solver
   (``vv-principles`` structural-independence; MMS is **not** paired
   here because MMS does not prove rates against an eigenvalue —
   the rate is a closed-form property of the cross sections).

Jacobi vs Gauss-Seidel — the Phase 3 boundary recovery
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Wave O BC extraction (steps O.4a.2 + O.4b, Issue #208) made the
2-D sweep **bare**: it reads ``psi.boundary.inflow`` as *given* for
the whole sweep, and the reflective coupling is applied externally
via the sibling :math:`-B` term (see :ref:`bare-sweep-extraction`).
A side effect of that architectural win was a *rate regression*: the
retired ``bc.apply``-inside-the-sweep read the **live** boundary
buffer mid-sweep (intra-sweep Gauss-Seidel), whereas the bare sweep
with a fully-lagged external :math:`B` is **inter-sweep Jacobi** —
same converged fixed point, slower SI rate.  Phase 3 recovers the
intra-sweep reflective coupling through a polymorphic, mesh-time
:class:`~orpheus.sn.sweep_schedule.SweepSchedule` without
re-entangling the bare sweep with the BC.  Jacobi and Gauss-Seidel
are the **same** uniform sweep-and-reflect loop — there is *no*
``if jacobi/gs`` branch in the iteration; the splitting is selected
*once* by choosing the schedule:

.. list-table:: The two within-group SI schedules
   :header-rows: 1
   :widths: 16 42 42

   * - Schedule
     - Octant grouping
     - Inter-group reflect
   * - **Jacobi**
       (``"jacobi"``)
     - ONE group containing every octant; :math:`B\,\psi_n` seeded
       once and **frozen** for the whole sweep.  Identical to the
       pre-recovery bare all-octants sweep.
     - None.  All octants read the same lagged inflow seed.
   * - **Gauss-Seidel**
       (``"gauss_seidel"``, default)
     - One group per in-plane octant
       (:class:`~orpheus.sn.sweep_graph.OctantLabel`), in quadrature
       sweep order.
     - After each group, its reflective **outgoing** faces are
       re-reflected (the face-restricted :math:`-B`,
       :meth:`SNBoundaryOperator.reflect_into_inflow`), so a *later*
       group reads the **fresh** current-iterate inflow — the
       :math:`(L+C-B_{\rm lower})^{-1}` forward substitution.

In the Gauss-Seidel schedule, octants swept **before** their
specular partner keep the lagged seed (the cyclic :math:`B_{\rm
upper}` back-edges — a both-faces-reflective axis is a 2-cycle, so
one pass is only *partial* G-S); octants swept **after** read the
fresh value (the order-respecting :math:`B_{\rm lower}` edges).  The
schedule is a **mesh-time derived object** — it depends only on the
quadrature's octant partition and the mesh's reflective-face set,
not on fluxes, sources, or iteration state — so it is built once and
reused across every SI iterate (the same lifetime contract as
:class:`~orpheus.sn.sweep_graph.SweepDependencyGraph`).

The selection lives in :func:`~orpheus.sn.solver._select_si_resolvent`:
``"gauss_seidel"`` on a 2-D Cartesian mesh returns
``(_GaussSeidelResolvent(L+C, B, schedule), (S,))`` — :math:`B` moves
*into* the resolvent while :math:`S` stays a lagged gain;
``"jacobi"`` (and any 1-D mesh) returns ``(L+C, (S, B))`` — both
:math:`S` and :math:`B` lagged.  In **both** cases :math:`S` is a
lagged gain: the sweep never re-scatters mid-sweep.  Only the
boundary coupling gets the Gauss-Seidel treatment.

.. warning::

   **Honest scope — boundary G-S is NOT a scattering accelerator.**
   The recovery folds **only** the boundary reflection :math:`B`.
   It therefore accelerates the *boundary-layer transient*, NOT the
   dominant flat *scattering* :math:`c`-mode of
   :eq:`si-spectral-rate`.  The measured gain is a constant,
   regime-independent **~0.86–0.92×** (e.g. :math:`n_{\rm GS}\approx
   641` vs :math:`n_{\rm Jacobi}\approx 697` on a B-2g reflective
   box), **not** the :math:`c^2`-halving (≈0.5×) one might naively
   expect from "Gauss-Seidel".  The :math:`c^2`-halving is the
   *scattering* Gauss-Seidel result, which does **not** apply to
   boundary-only G-S (the scattering :math:`S` is still fully
   lagged).  The dominant within-group scattering rate is recovered
   ONLY by **Krylov** (already production; rate-optimal,
   splitting-invariant — :math:`n\approx302` vs SI's :math:`n\approx
   860` on the same slab at :math:`\varepsilon=10^{-10}`) or by
   **consistent DSA** (a future feature, GitHub issue #2, with
   :math:`\rho\approx0.22` independent of :math:`c`).  A future
   reader must not mistake boundary-G-S for a scattering accelerator.

Why the scattering :math:`c`-mode cannot be folded into the sweep
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is tempting to try to fold the within-group self-scatter
:math:`\Sigma_{s,0}^{g\to g}` *into* the sweep — to accelerate the
:math:`c`-mode by absorbing the self-scatter into the cell-balance
denominator as a **removal cross-section** :math:`\sigma_r =
\Sigma_t - \Sigma_{s,0}^{g\to g}`, then iterating only on the
residual scattering :math:`\psi_{n+1} = A_{\rm wg}^{-1}(S_{\rm
residual}\psi_n + q)` with a :math:`\sigma_r`-sweep as :math:`A_{\rm
wg}^{-1}`.  **This is a latent correctness trap** (GitHub issue
#215, measured 2026-06-05), and documenting *why* it fails prevents
a future session from re-attempting it.

The σ\ :sub:`r`-sweep inverts :math:`(\hat\Omega\cdot\nabla +
\sigma_r\,\mathbb{I})` — a removal that is **diagonal in angle**.
But the within-group self-scatter is :math:`S_{\rm foldable} =
\Sigma_{s,0}\,P_{\rm iso}` — the **isotropic-projection** operator
(rank-1 in angle, :math:`\phi/\!\sum_n\mathcal{W}_n`).  The two
operators **coincide only for isotropic flux**:

.. math::
   :label: si-sigma-r-fold-mismatch

   \underbrace{\sigma_r\,\mathbb{I}}_{\text{diagonal in angle}}
   \;\ne\;
   \underbrace{\Sigma_{s,0}\,P_{\rm iso}}_{\text{isotropic projection}}
   \qquad\text{unless }\psi\text{ is isotropic, i.e. }
   \psi_n = \tfrac{\phi}{\sum_n\mathcal{W}_n}\ \forall n.

The consequence is a verification-mode-2 (variable-swap / operator
mismatch) bug that the *standard* test regime cannot see:

.. list-table:: σ\ :sub:`r`-fold failure across the BC regime (issue #215)
   :header-rows: 1
   :widths: 30 22 48

   * - Variant
     - Result
     - Why
   * - σ\ :sub:`r`-sweep approximation, **fully-reflective uniform**
       box
     - **exact** (round-off)
     - Flux is isotropic ⟹
       :math:`\sigma_r\mathbb{I}\equiv\Sigma_{s,0}P_{\rm iso}`.  The
       isotropic unit tests pass.
   * - σ\ :sub:`r`-sweep approximation, **anisotropic**
       (vacuum / heterogeneous)
     - **46–56 % flux error**
     - Flux is anisotropic ⟹ the diagonal removal is the wrong
       operator; the error is silent (no crash) and corrupts real
       cases.
   * - "exact" variant (keep the :math:`-\Sigma_{s,0}\!\odot\!\psi`
       remnant on the RHS)
     - **DIVERGES**
     - The remnant gain has spectral radius
       :math:`\Sigma_{s,0}/\sigma_r \approx 39` — the splitting is
       unstable.

This is the textbook reason **DSA needs a *consistent* diffusion
operator**: the correct synthetic acceleration of the
isotropic-projection self-scatter is a diffusion solve whose removal
matches the transport operator's low-order limit, not a directional
sweep with a doctored denominator.  The
:meth:`ScatteringOperator.foldable_part <orpheus.sn.scattering.ScatteringOperator.foldable_part>`
/ :meth:`residual_part <orpheus.sn.scattering.ScatteringOperator.residual_part>`
split (the data API landed under Issue #197 PR-TYPED-1) produces
:math:`\Sigma_{s,0}^{g\to g}` precisely as the input a DSA
preconditioner consumes (the diffusion removal coefficient) — it is
the right input **for DSA**, NOT for a sweep fold.  Any future
within-group accelerator MUST be gated on an **anisotropic** config;
the isotropic box cannot see this error (``vv-principles``
anti-pattern #4: homogeneous/isotropic verification is blind to the
angular structure).

The diagonal-cubature shared-face rule (ERR-056)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Gauss-Seidel schedule must assign each reflective face to the
**LAST** octant group (in sweep order) that outflows through it —
NOT the first.  This is a **correctness** requirement, not merely a
rate optimisation, and the distinction is invisible on axis-aligned
quadratures:

* On an **axis-aligned** quadrature
  (:meth:`Quadrature.product <orpheus.numerics.quadrature.Quadrature.product>`
  — each octant outflows a single face), every reflective face is
  outflowed by exactly **one** octant group, so "last" trivially
  equals "the only one".
* On a **diagonal / spherical** cubature
  (:meth:`Quadrature.lebedev <orpheus.numerics.quadrature.Quadrature.lebedev>`,
  :meth:`level_symmetric <orpheus.numerics.quadrature.Quadrature.level_symmetric>`
  — each octant outflows **two** in-plane faces), a face is shared
  by :math:`\ge 2` octant groups (e.g. ``xmax`` is outflowed by
  every ``+x`` octant: :math:`(+1,0)`, :math:`(+1,+1)`,
  :math:`(+1,-1)`).

Reflecting a shared face after only the **first** outflowing group
absorbs the *not-yet-swept* octants' slots — and because the interior
cochain :math:`C^1_{\rm int}` (carried by the rolling
``_MovingFrontier``) is rebuilt and :math:`\iota_*`-seeded each
``.solve``, those slots still hold the **inflow seed**, not real
outflow.  The reflect then
propagates garbage and the iteration converges to the fixed point of
the **wrong** operator (ERR-056).  Crucially, this seed-contamination
does **not** self-correct at convergence — unlike a lagged-but-real
Jacobi coupling, which reads the previous iterate's genuine value.
Deferring the reflect to the **last** outflowing group guarantees
the face's outflow is complete before it is reduced; octants reading
the face that are swept *before* its reflect simply keep the lagged
seed (the valid cyclic back-edge → partial one-pass G-S).  This is
the general principle that *a face shared by multiple work-units
must be reduced only after the last contributing unit completes* —
the same fan-in discipline KBA wavefront scheduling
([Pautz2002]_; Adams & Larsen 2002 §VI on parallel sweeps) and
multigroup Gauss-Seidel over shared down-scatter targets require.
The full post-mortem (symptom, root cause, the
``test_gs_diagonal_quadrature_shared_face_assigned_to_last_group_only``
structural pin) is catalogued as **ERR-056**.

White and vacuum faces are **excluded** from the reflective set:
vacuum has no coupling (:math:`B=0`), and white reflection couples
*all* ordinates on a face, so the octant-order G-S degenerates to
Jacobi anyway — only **specular** reflection admits the
order-respecting forward-substitution acceleration.

Numerical evidence
~~~~~~~~~~~~~~~~~~~~

All measured 2026-06-05 on this branch
(:func:`tests.sn.verification.analytical.test_si_convergence_rate`,
GL/``product`` :math:`N=8`, ``inner_tol`` as noted).  The Jacobi and
Gauss-Seidel counts are compared **in-process** (no hardcoded
baseline — Jacobi is a permanent live control, so the gates cannot
go stale).

.. list-table:: SI iteration counts — boundary G-S recovery vs the splitting-invariant controls
   :header-rows: 1
   :widths: 34 13 13 13 27

   * - Configuration
     - :math:`n_{\rm Jacobi}`
     - :math:`n_{\rm GS}`
     - ratio
     - Notes
   * - B-2g reflective **box** (2-D, ``inner_tol`` 1e-8)
     - 697
     - 641
     - 0.92×
     - Same fixed point: scalar-flux rel-L\ :sub:`∞`
       :math:`4.86\times10^{-8}` (rate-only, ``vv-principles``
       Mode 9).
   * - B-2g reflective **slab** (1-D, ``inner_tol`` 1e-8)
     - 655
     - 655
     - 1.00×
     - **No-op** by design — the 1-D scan is not a wavefront;
       ``_select_si_resolvent`` falls back to Jacobi.
   * - B-2g **vacuum** slab (G-4 negative control)
     - 128
     - 128
     - 1.00×
     - :math:`B=0` ⟹ the schedule is inert; proves the recovery
       touches *only* reflective coupling.
   * - B-2g reflective slab, **Jacobi vs Krylov**
       (``inner_tol`` 1e-10)
     - ≈860
     - — (302)
     - 2.85×
     - The splitting-invariant **Krylov** floor: an SI splitting can
       never beat it (the rate-optimal anchor, not the target).

The analytic Jacobi anchor (:eq:`si-spectral-rate`) predicts
:math:`\log(10^{-8})/\log(0.975) = 728` for the B-2g slab; the
measured 655 gives ratio **0.90** — the finite-slab leakage +
multigroup correction discussed in the note above.  The
**eigenvalue** path surfaces the analogous measurand via
:attr:`IterationHistory.total_inner_iterations
<orpheus.numerics.eigenvalue.IterationHistory.total_inner_iterations>`
(the Phase 3 measurement seam): A-2g reflective :math:`n=10` gives
SI ``total_inner`` 371, Krylov 310 (with :math:`n_{\rm outer}=3` for
both — the **outer** count is splitting-invariant; the inner SI count
is where the recovery shows).

For comparison, a clean 1-D textbook DSA spike (the future issue-#2
target) gives **8–21×** :math:`c`-independent speed-up on a vacuum
slab (:math:`\rho\approx0.22`), but a *naive* finite-difference DSA
**diverges** on a reflective boundary — the concrete confirmation
that DSA needs the consistent diffusion operator the σ\ :sub:`r`-fold
trap above already implied.

KEigenvalue: outer power iteration
-----------------------------------

:class:`~orpheus.numerics.iteration.KEigenvalue` poses the
k-eigenvalue problem from its operator triple and **delegates** the
outer loop to the canonical
:func:`~orpheus.numerics.eigenvalue.power_iteration` (one loop engine;
see :ref:`eigenvalue-posing`).  Each outer step (run by
``power_iteration``) is classical power iteration on the
:math:`k`-update, with :class:`SourceIteration` driving the inner
fixed-source solve:

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
``keff_estimator`` callable.  This is the K-specific volume-weighting
hook: relocating it into the algorithm as a Rayleigh update on the
resolvent :math:`A_{\rm loss}^{-1}M` (so the loop is *literally*
K/α-agnostic) is the α-eigenvalue wave's first step — see the
:ref:`eigenvalue-posing` honest-scope note.

The ``inverter`` parameter (the Wave-E ERR-026 reconciliation hook)
-------------------------------------------------------------------

.. note:: **Re-framed (2026-06-12, Issue #195).**

   This section's motivating premise — that the WDD sweep ``L.solve``
   has a *closure-bias-driven* fixed point distinct from the symmetric
   ``apply`` closure, so routing to Krylov-on-``apply`` is what "closes
   ERR-026" — was the *two-distinct-closures* picture.  ERR-058 (#195)
   showed the curvilinear wrong fixed point was the *closure-seed*
   family; once the seeds are fixed the sweep and the matvec are ONE
   discrete system and SI ``L.solve`` :math:`\equiv` Krylov-on-``apply``
   bit-identical (see :ref:`sn-err-058-closure-seed-closeout`).  The
   ``inverter`` hook **remains a real, retained generality** — it lets
   the iteration primitives swap :math:`L^{-1}` realisation (direct
   solve, sweep, Krylov-on-apply with a preconditioner) without
   re-implementation — but it is no longer *required* to obtain the
   correct curvilinear fixed point.  Read the WDD-bias-vs-Krylov
   contrast below as Wave-E-era history.

Both primitives accept an ``inverter: Callable[[ndarray], ndarray]``
that supplies :math:`L^{-1}`.  When ``None``, the primitive routes
through :meth:`L.solve`.  When supplied, the caller controls how
:math:`L^{-1}` is realised — the load-bearing design choice in the
Wave-E reconciliation of ERR-026.

* ``inverter = None`` (default):  :math:`L^{-1}\,q = L.solve(q)`.
  For an SN sweep this is the WDD asymmetric closure — which has a
  closure-bias-driven self-consistent fixed point on curvilinear
  meshes that is **not** the fine-mesh-limit transport solution
  (ERR-026).
* ``inverter = lambda q: KrylovAcceleration(L, ...).solve(q)[0]``:
  Krylov-on-:meth:`apply` (the symmetric closure of
  :class:`~orpheus.sn.operator.InvertibleOperator`), with the
  sweep injected as a preconditioner :math:`M`.  The ORPHEUS↔scipy
  boundary is internal to
  :class:`~orpheus.numerics.iteration.KrylovAcceleration` (a single
  ravel-aware adapter wraps :meth:`apply` for scipy's GMRES).  This is the
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

* ``L`` MUST advertise :py:data:`CAP_APPLY`.
* ``L`` MUST advertise :py:data:`CAP_SOLVE` *or* the caller MUST
  supply ``inverter``.  Without one of those, the iteration cannot
  evaluate :math:`L^{-1}`.
* ``S`` MUST advertise :py:data:`CAP_APPLY`.  Pass
  :class:`~orpheus.numerics.operator.ZeroOperator` for the
  scattering-free case.
* ``F`` MUST advertise :py:data:`CAP_APPLY`.  For
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

* :class:`~orpheus.sn.operator.InvertibleOperator` — the
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
* :ref:`sn-solver-operator-algebra-coordinator` — how the SN solver
  consumes :class:`SourceIteration` / :class:`KrylovAcceleration`
  directly through the within-group SSoT helpers.
* :ref:`eigenvalue-posing` — the method-agnostic four-layer
  architecture: where :func:`power_iteration` (the canonical Layer-4
  engine), the K-row posing, and the resolvent inner solve sit, and
  the :math:`\alpha`/adjoint/transient/full-spectrum future seams.

.. _cross-solver-eigenvalue-consumers:

Cross-solver consumers of ``power_iteration``
---------------------------------------------

:func:`orpheus.numerics.eigenvalue.power_iteration` is the
**canonical** Layer-4 power-method algorithm — NOT a legacy
primitive.  It iterates over the method-agnostic
:class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol
boundary (the *late-bound* resolvent layer), which is **strictly more
general** than the operator-triple form: it admits both the
:math:`(L, S, F)`-triple resolvent (SN, MoC) *and* the
**monolithic-matrix resolvent** (CP, diffusion, homogeneous) that has
no separable :math:`(L-S)^{-1}` factor.  All five solver families are
therefore **co-consumers of the same canonical boundary**, each
supplying its own ``EigenvalueSolver``-Protocol realization of the
Layer-3 resolvent.  There is no migration to a single ``KEigenvalue``
engine — and no retirement of ``power_iteration`` — because no narrower
layer can express the no-triple families without manufacturing
fictitious :math:`L`, :math:`S` operators for methods that have no
sweep.  See :ref:`eigenvalue-posing` for the full four-layer
architecture and why the Protocol layer is canonical.

* **SN** (discrete ordinates) — drives ``power_iteration`` directly
  via :func:`~orpheus.sn.solver.solve_sn`; its Layer-3 resolvent is
  the within-group :class:`SourceIteration` /
  :class:`KrylovAcceleration` inner solve built from the
  :func:`_within_group_triple` SSoT (see
  :ref:`sn-solver-operator-algebra-coordinator`).
* **CP** (collision-probability) — drives ``power_iteration`` through
  its own ``EigenvalueSolver``-Protocol implementation; its resolvent
  is **one BiCGSTAB on a monolithic collision-probability matrix**,
  which has no :math:`(L, S, F)` split.  This is exactly the family
  the late-bound Protocol layer exists to admit.
* **Diffusion** — drives ``power_iteration`` with a finite-difference
  resolvent; the BiCGSTAB inner loop *is* the
  :math:`A_{\rm loss}^{-1}` action.
* **MoC** (method of characteristics) — drives ``power_iteration``
  with a track-based inner sweep as its resolvent, via the same
  late-bound boundary :class:`InvertibleOperator.solve` exposes.
* **Homogeneous** — drives ``power_iteration`` over a direct linear
  solve; the analytical
  :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
  is the closed-form algebra-of-record this family also realizes.

The :class:`~orpheus.numerics.iteration.KEigenvalue` adapter is **one
Layer-2b implementer** of this boundary, for callers who *have* a
natural :math:`(L, S, F)` triple and want to skip writing a full
solver class; its :meth:`~orpheus.numerics.iteration.KEigenvalue.solve`
delegates its loop to ``power_iteration`` (one engine — see the
same-morphism analysis in :ref:`eigenvalue-posing`).  Making the
Layer-4 loop *literally* K/α-agnostic — relocating the eigenvalue
scaling out of the K-specific
:meth:`~orpheus.numerics.eigenvalue.EigenvalueSolver.compute_keff` and
renaming the K-flavoured Protocol methods to ``eigen_operator`` /
``mu_to_eigenvalue`` — touches all five families' Protocol surface and
is the **first step of the α-wave**, deferred under *unify-after-two*
because only the k-row exists today (the α-generic agnostic relocation
future seam, :ref:`eigenvalue-posing`).

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
  hook's use case — the ERR-026-closing curvilinear Krylov path).


.. _sn-solver-operator-algebra-coordinator:

SNSolver as an operator-algebra coordinator
============================================

:class:`~orpheus.sn.solver.SNSolver` consumes the operator triple
:math:`(L, S, F)` directly.  At construction time, the solver builds:

* :attr:`SNSolver.L` — :class:`InvertibleOperator` carrying the
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
plumbing.  The within-group inner solve is built once from a single
source of truth — the :func:`_within_group_triple` helper assembles
the honest within-group decomposition :math:`(L+C,\ S,\ B)`: the
invertible resolvent plus its two lagged coupling gains (the bulk
scattering :math:`S` and the trace boundary reflection :math:`B`; zero
within-group fission), handed to the **variadic** driver
:math:`\text{Driver}(L_{\rm resolvent},\,*\text{gains})` (Wave O step
O.2a — the transitional :math:`S + B` fold is retired; see
:ref:`bc-extraction-variadic-driver` in :doc:`operator_algebra`).
:func:`_within_group_krylov` wraps the matching
:class:`~orpheus.numerics.iteration.KrylovAcceleration` — and the
decomposition is shared verbatim across the eigenvalue source-iteration
inner (:meth:`SNSolver._solve_source_iteration`), the eigenvalue Krylov
inner (:meth:`SNSolver._solve_krylov`), and both fixed-source paths.

The within-group inner solve consumes the primitives directly
-------------------------------------------------------------

:class:`SNSolver`'s within-group inner solve **is** the
:class:`~orpheus.numerics.iteration.SourceIteration` /
:class:`~orpheus.numerics.iteration.KrylovAcceleration` primitive — not
a verbatim replica of its loop.
:meth:`SNSolver._solve_source_iteration` constructs a
:class:`SourceIteration` from the :func:`_within_group_triple` SSoT and
runs it; :meth:`SNSolver._solve_krylov` constructs a
:class:`KrylovAcceleration` from :func:`_within_group_krylov` and runs
that.  The Layer-3 resolvent of the SN row in the
:ref:`eigenvalue-posing` architecture is exactly these primitive
instances.

The primitive is **type-agnostic and angular-capable**: it operates on
the typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
composite, which carries the full angular flux on its bulk.  Pℓ
anisotropic scattering therefore rides the angular bulk with no special
plumbing — :meth:`ScatteringOperator.apply` on a ``TimedFullField``
reads the angular moments off the composite and builds the anisotropic
source via :meth:`ScatteringOperator.build_aniso_source`, all inside
the primitive's normal RHS path.  There is **no scalar-flux
limitation** and **no pending "Approach A" cleanup**: the earlier
framing — that :class:`SourceIteration` carried only scalar flux and SN
had to replicate the loop verbatim until the angular state could be
threaded through — was a property of an interim scalar-only carrier
that the typed composite retired.  The
``.claude/skills/algebra-of-record`` "Branch 2 implements the same
operator algebra" discipline is satisfied: SN is the discretized
Branch-2 consumer of the shared primitive, not a parallel loop.

The (L − S − F)·ψ = (1/k)·F·ψ framing at the solver level
-----------------------------------------------------------

Beyond driving the within-group inner solve, the :math:`(L, S, F)`
framing organises the solver's outer API surface:

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
  weighted production / absorption ratio); the SN-specific volume
  weighting lives here rather than in the generic
  :meth:`KEigenvalue.solve` default Rayleigh-quotient estimator, which
  is volume-blind.

The solver-level :math:`1/k` scaling (in
:meth:`~SNSolver.compute_fission_source`) and the volume-weighted
eigenvalue estimate (in :meth:`~SNSolver.compute_keff`) are exactly the
points where SN's specifics live; the rest of the solver is
operator-algebra coordination over the canonical
:func:`~orpheus.numerics.eigenvalue.power_iteration` boundary.  These
two K-specific hooks are also precisely *why* the Layer-4 loop is not
yet literally K/α-agnostic — relocating the eigenvalue scaling into the
algorithm is the first step of the α-wave (see the honest-scope caveat
in :ref:`eigenvalue-posing`).

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
  :meth:`InvertibleOperator.apply`)
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
:math:`(N, n_g, n_x, n_y)` (Issue #196 PR-INDEX-5 principled layout,
the ``g`` axis directly after ``N``) and threads it through the sweep's
:math:`Q_{\rm aniso}` slot --- merging additively with any P1+
scattering contribution the solver itself builds.  This bare-array form
is the **bulk-only / vacuum** special case of the composite source
``q = q_bulk ⊕ q_∂`` the solver also accepts (see
:ref:`sn-composite-fixed-source`); this isotropic slab MMS is
vacuum-automatic, so the boundary leaf is identically zero.  Vacuum
boundary conditions are applied via the mesh-level BC infrastructure
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
   :func:`~orpheus.sn.loss_representation.transport_sweep` accepted a
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
(:func:`orpheus.sn.loss_representation._sweep_jacobi`) with diamond-difference
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
  :func:`orpheus.sn.loss_representation._sweep_jacobi` (the 2D diamond-difference
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


.. _sn-mms-curvilinear-isotropic-verification:

Curvilinear isotropic MMS — radial DD-closure probe
----------------------------------------------------

Phase 3.4 of the verification campaign extends the slab MMS
(:eq:`sn-mms-psi` / :eq:`sn-mms-qext`) to 1-D **spherical** and
1-D **cylindrical** geometries with the simplest non-trivial trial
solution that respects the vacuum-at-outer and symmetry-at-origin
boundary conditions: an **isotropic** ansatz
:math:`\psi_n(r) = A(r)/W`.  By construction the angular
redistribution operator vanishes on this ansatz
(:math:`(1-\mu^2)/r \cdot \partial\psi/\partial\mu = 0` for the
sphere; :math:`-(1/r)\,\partial(\xi\psi)/\partial\varphi = 0` for
the cylinder), so the only spatial-discretisation error that drives
the measured convergence rate is the **radial DD closure**.  The
isotropic case is therefore the focused L1 probe for the
streaming + removal path; the angular redistribution path is
covered by the companion anisotropic case
(:ref:`sn-mms-curvilinear-aniso-verification` below — a deliberate
pairing that defeats the ``vv-principles`` Mode 7 "MMS
simplification bias" failure mode).

**Spherical isotropic ansatz.**  For a vacuum-BC sphere of radius
:math:`R` with reflective inner BC at :math:`r=0` in one energy
group, pick

.. math::
   :label: sn-mms-spherical-psi

   \psi_n(r) = \frac{1}{W}\,A(r),
   \qquad A(r) = \sin\!\left(\frac{\pi r}{R}\right),

with :math:`W = \sum_n w_n = 2` for symmetric Gauss--Legendre.
Because :math:`A(0) = A(R) = 0`, every ordinate vanishes at both
the symmetry centre and the vacuum outer face — both BC kinds are
satisfied automatically.  Since :math:`\psi_n` is independent of
ordinate, the scalar flux recovered by any quadrature order is
exactly :math:`\phi(r) = A(r)`.

**Spherical manufactured source.**  Substituting
:eq:`sn-mms-spherical-psi` into :eq:`transport-spherical` and
using that :math:`(1-\mu^2)\,\partial_\mu\psi/r \equiv 0` for an
isotropic flux gives

.. math::
   :label: sn-mms-spherical-qext

   Q^{\text{ext}}_n(r)
        = \mu_n\,A'(r)
        + (\Sigma_t - \Sigma_s)\,A(r)
        = \mu_n\,\frac{\pi}{R}\cos\!\left(\frac{\pi r}{R}\right)
          + (\Sigma_t - \Sigma_s)\sin\!\left(\frac{\pi r}{R}\right).

This is structurally identical to the slab source
:eq:`sn-mms-qext` with :math:`x \to r` — the spherical
:math:`(2/r)\partial_r` curvature term and the angular
redistribution term both vanish on the isotropic ansatz, leaving
the per-ordinate streaming + removal balance as the residual.

**Cylindrical isotropic ansatz.**  The radial direction cosine for
1-D cylindrical is :math:`\eta_n = \sin\theta_n \cos\varphi_n`.
Use

.. math::
   :label: sn-mms-cylindrical-psi

   \psi_n(r) = \frac{1}{W}\,A(r),
   \qquad A(r) = \sin\!\left(\frac{\pi r}{R}\right),

with the same :math:`W = \sum_n w_n` for the cylindrical Product or
LS quadrature.  Symmetric Product quadrature gives
:math:`\sum_n w_n \eta_n = 0`, so :math:`\phi(r) = A(r)` exactly.

**Cylindrical manufactured source.**

.. math::
   :label: sn-mms-cylindrical-qext

   Q^{\text{ext}}_n(r)
        = \eta_n\,A'(r) + (\Sigma_t - \Sigma_s)\,A(r).

The cylindrical curvature term :math:`-(1/r)\,\partial(\xi\psi)/\partial\varphi`
vanishes by isotropy of :math:`A(r)`, the same way the spherical
:math:`(1-\mu^2)/r \cdot \partial_\mu\psi` vanishes; the radial
streaming :math:`\eta_n A'(r)` and the removal
:math:`(\Sigma_t - \Sigma_s)A(r)` carry the residual.

**Risk point — Mode 7 ansatz bias.**  Per ``vv-principles`` failure
Mode 7 ("MMS simplification bias"), the isotropic ansatz is
deliberately structured to NULL the angular redistribution path.
A passing :math:`\mathcal{O}(h^2)` convergence here is necessary
evidence for the radial DD closure but it is *not* sufficient for
the full curvilinear sweep — ERR-026 (the curvilinear sweep WDD
flux-shape bug) is mathematically invisible to this MMS because the
redistribution term that ERR-026 lives on cancels by ansatz
construction.  The companion anisotropic case
(:ref:`sn-mms-curvilinear-aniso-verification`) is the load-bearing
sufficient evidence for the full sweep; both are required.

**Code pointers.**

- Derivation:
  :class:`orpheus.derivations.continuous.mms.sn.SNSphericalMMSCase`,
  :class:`orpheus.derivations.continuous.mms.sn.SNCylindricalMMSCase`,
  :func:`orpheus.derivations.continuous.mms.sn.build_spherical_mms_case`,
  :func:`orpheus.derivations.continuous.mms.sn.build_cylindrical_mms_case`.
- Tests:
  :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_spherical_mms_converges_second_order`
  (sphere) and
  :func:`tests.sn.verification.mms.test_mms_curvilinear.test_sn_cylindrical_mms_converges_second_order`
  (cylinder), both ``catches("ERR-058")``.  **Their ``xfail`` markers
  came off 2026-06-12 with the ERR-058 closure-seed fix** (Issue #195).
  Post-fix the ladders are clean second-order with SI :math:`\equiv`
  Krylov bit-identical — sphere ``[1.49e-2, 3.73e-3, 9.28e-4, 2.31e-4,
  5.74e-5]`` (orders 2.00–2.01), cylinder ``[2.16e-3, 5.39e-4, 1.35e-4,
  3.37e-5]`` (orders 2.00); the magnitude band
  :math:`[10^{-8}, 10^{-3}]` is met (sphere :math:`n_x\ge 80`,
  cylinder :math:`n_x\ge 40`).  Through the bug era (Wave E Round 3
  2026-05 → ERR-058) they were ``xfail(strict=True)`` — the
  now-superseded "pre-asymptotic transient" reading; see
  :ref:`sn-err-058-closure-seed-closeout`.


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

**Spatial-convergence claims.**  Diamond-Difference is design-order
:math:`\mathcal{O}(h^2)` in the cell width (:eq:`dd-cartesian-1d` /
:eq:`dd-curvilinear-scalar`); the curvilinear anisotropic L1 claim
asserts that the **measured** scalar-flux error against the
manufactured solution :eq:`sn-mms-spherical-aniso-psi` /
:eq:`sn-mms-cylindrical-aniso-psi` falls at the same rate.  For the
sphere,

.. math::
   :label: sn-mms-spherical-aniso-spatial-convergence

   \bigl\|\phi_h(r) - A(r)\bigr\|_{L^2(\Omega)}
        \;=\; \mathcal{O}(h^2)
        \qquad \text{as } h = R/n_x \to 0\,,

with the convergence ORDER (slope of :math:`\log\|\phi_h - A\|`
versus :math:`\log h` over the last two mesh halvings) the L1
acceptance criterion ``min(orders[-2:]) > 1.9``.  The cylindrical
analogue,

.. math::
   :label: sn-mms-cylindrical-aniso-spatial-convergence

   \bigl\|\phi_h(r) - A(r)\bigr\|_{L^2(\Omega)}
        \;=\; \mathcal{O}(h^2)
        \qquad \text{as } h = R/n_x \to 0\,,

uses the same acceptance criterion on the cylindrical-aniso
ansatz.  Both labels are consumed by the
:file:`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py`
gate-3 tests, which **stay ``xfail``** — but, post-ERR-058 (Issue
#195), no longer for the wrong-fixed-point reason.  The ERR-058
closure-seed fix recovered :math:`\mathcal{O}(h^2)` *spatial*
convergence (the isotropic ladders are clean second-order; see
:ref:`sn-err-058-closure-seed-closeout`).  These **anisotropic** rows
remain xfail because the angle-varying ansatz hits the
**fixed-quadrature angular floor** of the half-angle thread
interpolation: under spatial refinement at fixed quadrature the error
converges to a floor (sphere S16 :math:`\approx 7\mathrm{e}{-4}`,
cylinder :math:`n_\mu{=}4` :math:`\approx 1.9\mathrm{e}{-2}`) that
drops only under *quadrature* refinement.  Their markers are re-pinned
to the `Issue #229
<https://github.com/deOliveira-R/ORPHEUS/issues/229>`_
quadrature-aware retune (the regression gate that flips them to
unexpected-pass when the retune lands).  See
:ref:`sn-err-058-aniso-floor` for the floor-vs-quadrature evidence.

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


.. _sn-composite-fixed-source:

The composite fixed-source API — :math:`q = q_{\text{bulk}} \oplus q_\partial`
------------------------------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   - **A fixed source is a source EVERYWHERE.** The right-hand side of
     a fixed-source transport problem is not just a volumetric bulk
     source — it is a source on the whole phase space: a bulk
     :math:`q_{\text{bulk}}` *and* a boundary (prescribed-inflow)
     :math:`q_\partial`. ORPHEUS represents this as the composite
     ``q = q_bulk ⊕ q_∂``, the direct sum of the two role-typed leaves.
   - **The carrier is the object we already have.**
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` accepts the
     composite as a
     :class:`~orpheus.transport.timed_full_field.TimedFullField` — the
     **same** typed direct-sum carrier the SI / Krylov inner already
     flows through internally. This is *not* a new type; it is
     ergonomics to **generate** the right object (Cardinal Rule 2 — we
     already have the right concepts).
   - **A source, role-distinguished from a flux by its leaf types.**
     The composite's bulk leaf is an
     :class:`~orpheus.transport.source_sinks.AngularSourceSink` and its
     boundary leaf a
     :class:`~orpheus.transport.source_sinks.BoundarySourceSink` — the
     *source* column of the role grid (see :ref:`bc-extraction-operator-output-typing`).
     The iterate / solution it produces is a *flux* (``AngularFlux`` ⊕
     ``BoundaryFlux``). Same carrier shape, different role; the class
     gate keeps source and flux arithmetic from silently mixing.
   - **The legacy array is the bulk-only / vacuum special case.** Passing
     the historical ``(N, ng, nx, ny)`` ndarray is *exactly* the
     composite with an all-zero (vacuum) boundary. All 37 pre-existing
     callers keep working bit-unchanged.
   - **One construction point.** The private helper
     :func:`~orpheus.sn.solver._build_fixed_source_rhs` is the single
     place the RHS composite is built (Cardinal Rule 2 — it collapsed a
     ``q_ext_composite`` build that previously lived in *both* the SI
     and Krylov inner paths).
   - **The ergonomic boundary generator.**
     :meth:`~orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow`
     writes ONLY the inflow ordinate slots of the named faces (outflow
     slots of a prescribed inflow are physically meaningless →
     unrepresentable by construction, ``coding-elegance`` Pattern 4),
     leaving everything else zero. It is the known-per-face-array route;
     the lazy ``InflowSourceSpec``-recipe route (``from_spec``) is a
     distinct, still-deferred bridge.

The fixed-source right-hand side is a source on the whole phase space
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A fixed-source SN problem solves the affine within-group system
:eq:`si-within-group-operator-eq`

.. math::

   (L + C - S - B)\,\psi = q,

and the right-hand side :math:`q` is **not** a bulk volumetric source
alone. It has two pieces, one per phase-space locus:

* the **bulk** source :math:`q_{\text{bulk}}(\vec r, \hat\Omega, g)` —
  the per-ordinate volumetric external source :math:`Q^{\text{ext}}_n`
  on every cell;
* the **boundary** source :math:`q_\partial` — the prescribed inflow,
  the inhomogeneous term :math:`q` of the affine boundary law
  :eq:`affine-bc-form` :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`,
  living on the inflow ordinate slots of the boundary trace.

A vacuum boundary is simply :math:`q_\partial \equiv 0`; a non-vacuum
prescribed inflow is a non-zero :math:`q_\partial`. The natural object
is therefore the **direct sum** of the two:

.. math::

   q \;=\; q_{\text{bulk}} \,\oplus\, q_\partial,

an object that "represents the source everywhere". ORPHEUS already has
exactly this carrier: the
:class:`~orpheus.transport.timed_full_field.TimedFullField`, the typed
bulk⊕boundary(⊕history) direct sum that the within-group SI and Krylov
inner paths *already* pass around (the matvec
:math:`(L+C)\psi - (S+B)\psi - F\psi` and the SI rhs
:math:`F\psi + (S+B)\psi + q_{\text{ext}}` are CLOSED ``TimedFullField``
sums). The field-role-typing work did **not** introduce a new source
type — it surfaced the carrier we already had and added the ergonomics
to *generate* it (Cardinal Rule 2: we have the right object, we just
need a better way to build it).

Source vs flux — same carrier, different role
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The composite source and the angular-flux solution share the
``TimedFullField`` carrier *shape* but differ in their leaf **types**,
which encode the role:

.. list-table:: The composite carrier's two roles
   :header-rows: 1
   :widths: 22 39 39

   * - Locus
     - Source role (the RHS ``q``)
     - Flux role (the iterate / solution ``ψ``)
   * - bulk
     - :class:`~orpheus.transport.source_sinks.AngularSourceSink`
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
   * - boundary
     - :class:`~orpheus.transport.source_sinks.BoundarySourceSink`
     - :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`

The role-leaf types are the gate. A *source* and a *flux* are never the
same field even when they ride the same carrier — the
:class:`~orpheus.transport.timed_full_field.TimedFullField` class gate
(via :class:`~orpheus.numerics.field.Field`) rejects
``AngularSourceSink ± AngularFlux`` and the boundary analogue, so the
"RHS is a source, the iterate is a flux" distinction is illegal to mix
by construction. The completed boundary role grid mirrors the bulk
exactly (:ref:`bc-extraction-operator-output-typing`): an operator's
``.apply`` output is a *source/sink* (:math:`A\psi`), its ``.solve``
output is a *flux* (the swept solution trace), and a ``from_balance``
defect is a *residual*. ``q_\partial`` is a ``BoundarySourceSink``
because a prescribed inflow IS a source added to :math:`\gamma_-\psi`,
not the swept solution.

The two accepted forms of ``external_source``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:func:`~orpheus.sn.solver.solve_sn_fixed_source` accepts the
``external_source`` argument in either of two forms, normalised by the
single helper :func:`~orpheus.sn.solver._build_fixed_source_rhs`:

#. **A bare ``np.ndarray`` of shape** :math:`(N, n_g, n_x, n_y)` — the
   per-ordinate-density **bulk** source only, with a **vacuum**
   boundary. This is the original form, and it is *exactly* the
   composite with an all-zero boundary leaf
   (``BoundarySourceSink.zeros_on(sn_mesh)``). Every one of the 37
   pre-existing callers passes this form and keeps working bit-for-bit
   unchanged (the vacuum path is verified bit-identical).
#. **A full** :class:`~orpheus.transport.timed_full_field.TimedFullField`
   **composite** ``q = q_bulk ⊕ q_∂`` — the route for a **non-vacuum
   prescribed inflow**. Its leaf values are re-homed onto the solve's
   own ``sn_mesh``: the trace / grid layout is deterministic from
   ``(mesh, quadrature, materials)``, so this is an exact values-copy
   onto the solve's mesh instance, required because the within-group
   operators are built on ``sn_mesh`` and ``TimedFullField`` algebra
   enforces mesh identity.

.. code-block:: python

   from orpheus.sn import solve_sn_fixed_source
   from orpheus.sn.geometry import SNMesh
   from orpheus.transport.source_sinks import (
       AngularSourceSink, BoundarySourceSink,
   )
   from orpheus.transport.timed_full_field import TimedFullField

   sn = SNMesh(mesh, quadrature, materials)

   # Bulk volumetric source, per-ordinate density (N, ng, nx, ny).
   q_bulk = AngularSourceSink.from_mesh(Q_ext, sn)
   # Prescribed inflow: only the named faces' inflow ordinate slots.
   q_bndry = BoundarySourceSink.prescribed_inflow(
       sn, {"xmin": gamma_minus_xmin, "xmax": gamma_minus_xmax},
   )
   q = TimedFullField(bulk=q_bulk, boundary=q_bndry)

   result = solve_sn_fixed_source(materials, mesh, quadrature, q)

The legacy ``solve_sn_fixed_source(materials, mesh, quadrature, Q_ext)``
with a bare ``Q_ext`` array is identical to the above with a vacuum
``q_bndry`` (``BoundarySourceSink.zeros_on(sn)``).

The single construction point — Cardinal Rule 2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Before the field-role-typing work, the SI inner and the Krylov inner
each built their own ``q_ext_composite`` from the bulk array (the same
``AngularSourceSink.from_isotropic`` / ``from_mesh`` projection paired
with a zero boundary). That was a shared concept living in two places —
precisely the smell Cardinal Rule 2 flags.
:func:`~orpheus.sn.solver._build_fixed_source_rhs` collapses both into
one construction point: ``solve_sn_fixed_source`` calls it once, and
**both** inner paths consume what it returns. The helper:

* validates the bulk shape against :math:`(N, n_g, n_x, n_y)` (Issue
  #196 PR-INDEX-5 principled layout — the ``g`` axis directly after
  ``N``);
* for a bare array, pairs the bulk
  :class:`~orpheus.transport.source_sinks.AngularSourceSink` with a
  vacuum ``BoundarySourceSink``;
* for a composite, re-homes the leaf values onto the solve's
  ``sn_mesh`` (with a layout-size guard on the boundary trace), and
  raises a descriptive ``ValueError`` if the composite was built on an
  incompatible mesh / quadrature / materials.

The validation, the projection, and the vacuum-default boundary now
live in exactly one function. The SI and Krylov paths differ only in
the inner solve they run, not in how the RHS is assembled.

The ergonomic boundary generator — ``prescribed_inflow``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generating the boundary leaf :math:`q_\partial` is the part the
ergonomics target. The classmethod
:meth:`BoundarySourceSink.prescribed_inflow(mesh, {face: (N, ng) values}) <orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow>`
builds the prescribed inflow from known per-face arrays:

* for each face named in the mapping, it writes ONLY the **inflow**
  ordinate slots from the given :math:`(N, n_g)` array;
* **every other slot is left zero** — the outflow ordinate slots of a
  named face, and every slot of an unnamed (vacuum) face.

**Why outflow is unrepresentable (Pattern 4).** The outflow ordinate
slots of a *prescribed-inflow source* are physically meaningless: the
sweep determines the outflow trace, the source does not. Writing them
would be an illegal state. Rather than accept-then-ignore them,
:meth:`~orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow`
makes the illegal state **unrepresentable by construction** — it reads
the inflow ordinate index set
(:meth:`~orpheus.numerics.spaces.trace_space.TraceSpace.inflow_indices_for_face`)
and copies *only* those rows, so an accidentally-populated outflow row
in the caller's array simply cannot reach the field. This is
``coding-elegance`` Pattern 4 (illegal states unrepresentable). It
supersedes the ``zeros_on`` + nested
``face_view(face)[inflow] = …`` slot-fill loop that every
prescribed-inflow consumer (the non-vacuum MMS, the splitting-invariance
probe) previously hand-rolled — the single source of truth for
materialising a prescribed inflow onto the trace (Cardinal Rule 2).

**The recipe → snapshot distinction (vs ``from_spec``).** There are two
distinct routes to a boundary source, related as *recipe → snapshot*,
not as duplicates:

.. list-table:: Two routes to ``q_∂`` — known arrays vs lazy recipe
   :header-rows: 1
   :widths: 28 36 36

   * -
     - :meth:`~orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow`
     - ``BoundarySourceSink.from_spec`` (deferred)
   * - Input
     - known per-face ``(N, ng)`` arrays
     - a lazy
       :class:`~orpheus.geometry.boundary._source.InflowSourceSpec`
       recipe (``evaluate(shape) -> ndarray``)
   * - When
     - the inflow values are already computed (the MMS case)
     - the inflow is described by a per-face recipe evaluated on demand
   * - Status
     - **shipped** — the route the 4.6 MMS and the T4 probe use
     - **deferred** (``unify-after-two`` — no recipe-driven consumer
       that drives a typed boundary-source sweep yet)

The 4.6 MMS uses
:meth:`~orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow`
because it has explicit per-face arrays
(:math:`\gamma_-\psi = (A + \mu_n B)/W`); it does not need the lazy
recipe bridge, which waits for its first genuine consumer per the
``unify-after-two`` discipline.

Why this is ergonomics, not new types
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The entire change is the **generation** of objects that already
existed in the codebase, not the introduction of new ones:

* The carrier :class:`~orpheus.transport.timed_full_field.TimedFullField`
  pre-dates this work — it is what the inner solve already flows.
* The leaf types
  :class:`~orpheus.transport.source_sinks.AngularSourceSink` and
  :class:`~orpheus.transport.source_sinks.BoundarySourceSink` pre-date
  this work — they are the *source* column of the role grid.
* The affine boundary law :eq:`affine-bc-form` and the ``q.boundary``
  slot pre-date this work.

What was missing was the *ergonomic generator* for a non-vacuum
boundary leaf and a *public entry point* that accepts the composite.
:meth:`~orpheus.transport.source_sinks.BoundarySourceSink.prescribed_inflow`
and the second accepted form of
:func:`~orpheus.sn.solver.solve_sn_fixed_source` supply exactly those —
no more. This is the operational meaning of "we already have the right
objects/concepts — we just need better ergonomics to generate them":
the abstraction is unchanged; only the surface that builds it is
better. The first consumer is the non-vacuum prescribed-inflow MMS of
:ref:`sn-mms-nonvacuum`.


.. _sn-mms-nonvacuum:

Non-vacuum prescribed-inflow MMS (Phase 4 / O.2b 4.6)
------------------------------------------------------

.. admonition:: Key Facts
   :class: important

   - **What this section adds.** The entire pre-existing MMS catalog
     (:ref:`sn-mms-curvilinear-isotropic-verification`,
     :ref:`sn-mms-curvilinear-aniso-verification`) is
     *vacuum-automatic*: every manufactured ansatz vanishes at both
     boundaries, so the inflow trace :math:`\gamma_-\psi \equiv 0` on
     every ordinate and the prescribed-inflow source slot
     ``q.boundary`` is **identically zero** in all of them. Phase 4 /
     O.2b sub-step 4.6 fills that gap with a manufactured solution that
     is **non-zero at the outer face**, lighting the
     :math:`q.\text{boundary} \neq 0` path for the first time.
   - **The ansatz is the proven P1 element.** :math:`\psi_n = (A + \mu_n
     B)/W` is the same truncated-Legendre :math:`P_0 \oplus P_1` form
     used by the Phase 3.6 anisotropic cases
     (:eq:`sn-mms-spherical-aniso-psi`). 4.6 changes **only the
     boundary trace** — :math:`A,B` are chosen non-vanishing at the
     outer face — and reuses the verified angular structure. Linear in
     :math:`\mu` *fully* (not partially) activates the curvilinear
     redistribution; the question "do we need :math:`\mu^2`?" is
     answered **no** below.
   - **Two manufactured sources, derived from the continuous
     operator.** The slab source :eq:`sn-mms-nonvacuum-qext` has **no**
     redistribution term (the Cartesian operator lacks the
     :math:`\partial_\mu` coupling); the sphere source
     :eq:`sn-mms-nonvacuum-sph-qext` is the **same closed form** as the
     Phase 3.6 vacuum case — only :math:`A,B` differ. The spherical
     residual therefore lives in **one place**
     (:func:`~orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`,
     Cardinal Rule 2).
   - **HAZARD H1 (sphere pole regularity).** :math:`B(0)=0` is a HARD
     constraint — without it the redistribution :math:`(1-\mu^2)B/r \to
     \infty` at :math:`r=0`. The :math:`(r/R)` prefactor on the sphere
     :math:`B` enforces it; the slab has no pole, so a slab-style
     :math:`B(0)\neq 0` is fine there but **wrong** on the sphere.
   - **The affine-BC-to-RHS framing.** Prescribed inflow IS the
     inhomogeneous term :math:`q` of the affine boundary law
     :eq:`affine-bc-form` :math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q`,
     carried in the ``q.boundary`` slot
     (:class:`~orpheus.transport.source_sinks.BoundarySourceSink`) and
     consumed directly by :math:`(L+C)\text{.solve}` as the sweep
     inflow seed. **No** :class:`~orpheus.geometry.boundary.PrescribedInflow`
     mesh-BC bridge is touched, and **no** ``from_spec`` recipe bridge
     is needed — the inflow is supplied as the boundary leaf of the
     **composite source** ``q = q_bulk ⊕ q_∂`` that
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` now accepts
     directly (see :ref:`sn-composite-fixed-source`).
   - **The convergence rows drive the public composite-source API.**
     The 4.6 MMS no longer assembles the within-group operator triple
     by hand: :func:`~orpheus.sn.solver.solve_sn_fixed_source` accepts
     a :class:`~orpheus.transport.timed_full_field.TimedFullField`
     composite source, so each case bundles its manufactured bulk and
     prescribed-inflow boundary into one ``case.fixed_source(sn)``
     call. The migration off the operator-triple bypass *is* the
     retirement (retirement = test migration).
   - **The load-bearing assertion is the converged VALUE, not the
     rate.** Per the ``vv-principles`` skill anti-pattern #5 (rate is
     necessary, not sufficient), a silently-dropped ``q.boundary``
     converges cleanly at :math:`\mathcal{O}(h^2)` to the **wrong**,
     boundary-zero limit. Only the flux-value-vs-:math:`A(x)` check —
     with :math:`A` non-zero at the boundary (:math:`a_0>0`) — sees it.
   - **T3 (sphere) ships ``xfail(strict)``, re-scoped to Issue #229.**
     The slab rows are clean :math:`\mathcal{O}(h^2)` with value match.
     The sphere row's anisotropic :math:`(A+B\mu)/W` ansatz is
     angle-varying, so after ERR-058 (#195) recovered the spatial
     :math:`\mathcal{O}(h^2)` convergence it now hits the
     fixed-quadrature **angular floor** of the half-angle thread
     interpolation (sphere S16 floor :math:`\approx 7\mathrm{e}{-4}`,
     above the band) — NOT the old #195 plateau, and NOT a non-vacuum
     machinery failure (the boundary value *is* honoured).  The marker
     moved from #195 to the
     `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_
     quadrature-aware retune.  Both the sphere and slab rows now run the
     curvilinear/Cartesian **source-iteration** default of
     :func:`~orpheus.sn.solver.solve_sn_fixed_source` (SI :math:`\equiv`
     Krylov bit-identical post-ERR-058); the composite-source API
     delivers the prescribed inflow identically to every splitting (T4).
     The green companion T3g provides live structural coverage of the
     inflow + redistribution paths now.

This section narrates the Branch-1 SymPy algebra-of-record
(:mod:`orpheus.derivations.continuous.mms.sn`), the Branch-2 numpy
factories, and the L1 / foundation gates that verify the prescribed-
inflow discretisation. The verification chain follows the
``algebra-of-record`` discipline (Branch-1 SymPy reference, Branch-2
numpy production, structurally-independent L1 cross-check).

.. list-table:: 4.6 verification gates (measured, ``-O`` mode)
   :header-rows: 1
   :widths: 6 30 10 12 32

   * - Gate
     - Description
     - Level
     - Pillar
     - Status / evidence
   * - V_nonvac-slab
     - Slab substitution identity ``simplify(W·LHS − Σ_s φ − Q) == 0``
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_slab_substitution_identity`
   * - V_nonvac-sph
     - Sphere substitution identity (reuses the 3.6 spherical residual)
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_sph_substitution_identity`
   * - Decision-A pin
     - Parameterised :math:`A=B=` ``None`` reproduces 3.6 vacuum shapes byte-for-byte
     - foundation
     - regression
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_existing_spherical_aniso_still_passes_after_parameterization`
   * - L1 xcheck (slab)
     - Branch-2 numpy :math:`Q^{\text{ext}}` == lambdified SymPy (≤1e-13)
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_slab_nonvacuum_numerical_qext_matches_sympy`
   * - L1 xcheck (sphere)
     - Branch-2 numpy :math:`Q^{\text{ext}}` == lambdified SymPy (≤1e-13)
     - foundation
     - MMS (1C)
     - PASS — :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_sphere_nonvacuum_numerical_qext_matches_sympy`
   * - T1 (slab 1g)
     - DD :math:`\mathcal{O}(h^2)` + converged value + inflow honoured
     - L1
     - MMS (1C)
     - PASS — orders ``[2.04, 2.01]``, finest L2 ~1.2e-3, max\|φ−A\| ~8e-5
   * - T2 (slab 2g asym)
     - As T1, asymmetric downscatter :math:`\Sigma_s` (ERR-002 hazard)
     - L1
     - MMS (1C)
     - PASS — g0 ``[2.04, 2.01]``, g1 ``[2.05, 2.01]``
   * - T3 (sphere)
     - Curvilinear redistribution under non-vacuum inflow
     - L1
     - MMS (1C)
     - **xfail(strict)** on #229 — aniso angular floor ≈7e-4 (spatial
       O(h²) recovered by ERR-058; boundary value honoured)
   * - T3g (sphere)
     - Inflow honoured at :math:`r=R` + redistribution source live (green now)
     - foundation
     - structural
     - PASS — :func:`tests.sn.verification.analytical.test_mms_prescribed_inflow.test_sphere_nonvacuum_inflow_honoured_and_redistribution_live`
   * - T4 (Mode 9)
     - SI-Jacobi ≡ SI-Gauss-Seidel ≡ Krylov honour ``q.boundary``
     - foundation
     - splitting-invariance
     - PASS — pairwise reldiffs 1.3e-13 … 5.6e-13

Why the existing MMS catalog is vacuum-automatic
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every manufactured solution already in the SN verification ladder —
the isotropic :math:`\psi_n = A(r)/W` cases
(:ref:`sn-mms-curvilinear-isotropic-verification`) and the
anisotropic :math:`\psi_n = (A + \mu_n B)/W` cases
(:ref:`sn-mms-curvilinear-aniso-verification`) — was built with
:math:`A` and :math:`B` chosen to **vanish at both boundaries**. For
the canonical 3.6 sphere, :math:`A(r) = \sin(\pi r/R)` gives
:math:`A(0) = A(R) = 0`, and :math:`B(r) = (r/R)(1-r/R)\cos(\pi r/R)`
gives :math:`B(0) = B(R) = 0`. The slab isotropic case likewise uses
:math:`A(x) = \sin(\pi x/L)`.

The consequence is structural and total. The inflow trace of the
manufactured solution on any face is

.. math::

   \gamma_- \psi_n \big|_{\text{face}}
       = \frac{1}{W}\bigl(A(x_{\text{face}})
                         + \mu_n B(x_{\text{face}})\bigr)
       = \frac{1}{W}\bigl(0 + \mu_n \cdot 0\bigr) = 0
       \qquad \text{for every ordinate } n.

So the affine boundary law :eq:`affine-bc-form`
:math:`\gamma_-\psi = R\,G\,\gamma_+\psi + q` collapses to its
homogeneous (vacuum) form for these cases — and the inhomogeneous
inflow term :math:`q` is identically zero. The existing cases verify
the **interior** spatial / angular operator and the **homogeneous**
vacuum BC, but they say *nothing* about the prescribed-inflow path,
where a non-zero :math:`q` is pushed into the right-hand side. That
path is the one O.2b's field-role-typing work makes a first-class
boundary trace: an inhomogeneous inflow injected as the boundary
*source* slot. Until 4.6, no MMS row exercised it.

The fix is the smallest possible structural delta: keep the proven P1
angular form, keep the proven interior operator, and change **only**
:math:`A,B` so that they are non-zero at the outer face. Then
:math:`\gamma_-\psi \neq 0`, and the converged scalar flux
:math:`\phi(x) = A(x)` is non-zero at the boundary — exactly the
property the verification needs (see "The converged-value assertion"
below).

The ansatz — the P1 element and why linear-in-:math:`\mu` is enough
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slab ansatz and its manufactured source:

.. math::
   :label: sn-mms-nonvacuum-psi

   \psi_n(x) = \frac{1}{W}\bigl(A(x) + \mu_n B(x)\bigr),
   \qquad A(x) = a_0 + a_1\sin(kx),\quad B(x) = b_0\cos(kx),
   \quad a_0 > 0.

.. (vv-status rationale) definition: Definitional ansatz — the
.. manufactured angular flux is *imposed*, not solved for. Its
.. correctness as a reference is established by the source identity
.. :eq:`sn-mms-nonvacuum-qext` (SymPy ``simplify == 0``), not by a
.. property of this expression alone.
.. vv-status: sn-mms-nonvacuum-psi documented

The form :math:`\psi_n = (A + \mu_n B)/W` is the truncated Legendre
:math:`P_0 \oplus P_1` element: :math:`P_0(\mu) = 1` carries the
isotropic amplitude :math:`A`, and :math:`P_1(\mu) = \mu` carries the
first-moment amplitude :math:`B`. This is the **native angular basis**
of the SN closure — the Carlson half-angle pole closure folds the
moment source through :math:`P_\ell(-1) = (-1)^\ell`, so a linear-in-
:math:`\mu` input is exactly the lowest non-trivial moment with a
non-zero :math:`\partial_\mu`.

**Why linear-in-:math:`\mu` fully activates the redistribution.** The
curvilinear angular-redistribution operator is
:math:`\tfrac{1-\mu^2}{r}\,\partial_\mu\psi`. With :math:`\psi` linear
in :math:`\mu`, the angular derivative is a non-zero **constant** in
:math:`\mu`,

.. math::

   \frac{\partial \psi_n}{\partial \mu} = \frac{B(r)}{W} \neq 0,

and multiplying by the redistribution dome :math:`(1-\mu^2)` produces
a genuinely :math:`\mu^2`-structured term :math:`(1-\mu^2)B/r`. The
discrete closure that realises this operator — the Morel–Montry
half-angle recurrence with the Carlson :math:`\mu=-1` seed (see
:ref:`sn-pole-angular-closure-protocol`) — is **linear** in
:math:`\psi`. A linear operator is fully probed by any input that is
non-constant in its argument; the linear-in-:math:`\mu` ansatz is
non-constant in :math:`\mu`, so it exercises the entire linear
redistribution map (including the half-angle recurrence and the
second-moment coupling).

A quadratic-in-:math:`\mu` (P2) ansatz term would add **no** new
structural coverage of the redistribution. Because the closure is
linear, a quadratic input only changes *which point* in the operator's
already-fully-probed range you land on — it does not reach any term
the linear input misses. (A P2 term *would* additionally exercise the
:math:`\sum_n w_n \mu_n^2` quadrature-exactness, but that is a
property of the quadrature, not of the redistribution operator, and is
already covered elsewhere.) This settles the "do we need
:math:`\mu^2`?" question definitively: **no.** The verdict is recorded
in the cross-domain-attacker frame analysis (memo
``phase4_o2b_4_6_mms_ansatz_frame.md``, Q1/Q2) and is empirically
consistent with Phase 3.6, which uses exactly this linear-in-:math:`\mu`
ansatz and whose gate tests carry ``catches("ERR-026")`` — the
redistribution-bug catcher.

**The scalar flux is :math:`A`.** Because Gauss–Legendre (and every
symmetric quadrature on :math:`\mu \in [-1,1]`) satisfies
:math:`\sum_n w_n \mu_n = 0`, the first-moment term integrates out of
the scalar moment:

.. math::

   \phi(x) = \frac{1}{W}\sum_n w_n\bigl(A(x) + \mu_n B(x)\bigr)
           = \frac{1}{W}\Bigl(A(x)\sum_n w_n
                            + B(x)\underbrace{\sum_n w_n \mu_n}_{=\,0}\Bigr)
           = A(x).

This discrete identity is verified directly in
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_slab_nonvacuum_phi_equals_A_under_quadrature`
(≤1e-14 on a sample mesh). The reference scalar flux for the
convergence rows is therefore :math:`\phi_{\text{chosen}}(x) = A(x)`,
which — because :math:`a_0>0` — is **non-zero at the boundary**.

The manufactured slab source, derived from the continuous operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The slab SN transport operator (per ordinate, 1-group) is the
first-order streaming-plus-collision form
:math:`\mu\,\partial_x\psi + \Sigma_t\psi
= \tfrac{1}{W}(\Sigma_s\phi + Q^{\text{ext}})`. The Cartesian operator
has **no** angular-derivative term — the slab geometry produces no
angular redistribution (:eq:`transport-cartesian`). Substituting the
ansatz :eq:`sn-mms-nonvacuum-psi` and solving for the residual source:

.. math::

   \mu\,\frac{\partial}{\partial x}
       \frac{A + \mu B}{W}
   + \Sigma_t\,\frac{A + \mu B}{W}
   &= \frac{1}{W}\bigl(\Sigma_s\,A + Q^{\text{ext}}_n\bigr) \\[1mm]
   \frac{1}{W}\bigl(\mu A' + \mu^2 B'\bigr)
   + \frac{\Sigma_t}{W}\bigl(A + \mu B\bigr)
   &= \frac{1}{W}\bigl(\Sigma_s A + Q^{\text{ext}}_n\bigr),

where :math:`\phi = A` was used on the right. Multiplying through by
:math:`W` and isolating :math:`Q^{\text{ext}}_n` gives the closed form

.. math::
   :label: sn-mms-nonvacuum-qext

   Q^{\text{ext}}_n(x) = \mu_n A'(x) + \mu_n^2 B'(x)
                       + (\Sigma_t - \Sigma_s) A(x)
                       + \Sigma_t\,\mu_n B(x).

.. (vv-status rationale) derivation: A closed form obtained by
.. symbolic substitution into the continuous slab operator; verified
.. by SymPy ``simplify(W·LHS − Σ_s φ − Q) == 0`` and cross-checked
.. against the Branch-2 numpy evaluation to ≤1e-13.
.. vv-status: sn-mms-nonvacuum-qext verified

Note that there is **no** :math:`(1-\mu^2)B/r` term — the slab operator
simply does not generate it. The :math:`\mu^2 B'` term *is* present
(streaming the first-moment piece :math:`\mu B` gives
:math:`\mu \cdot \mu B' = \mu^2 B'`), so the slab still exercises the
second-moment streaming closure, just not the angular redistribution.
The :math:`(\Sigma_t - \Sigma_s)A` term is the within-group removal
net of isotropic self-scatter (:math:`c = \Sigma_s/\Sigma_t < 1`), and
:math:`\Sigma_t\,\mu_n B` is the collision of the first-moment piece.

The Branch-1 algebra-of-record is
:func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_slab_mms`
(building on
:func:`~orpheus.derivations.continuous.mms.sn._nonvacuum_slab_symbolic`),
which performs the substitution symbolically and proves
``simplify(W·LHS − Σ_s·φ − Q_closed) == 0``. Because the slab operator
lacks redistribution, it is a *genuinely different* operator from the
sphere and gets its own fresh symbolic pair — it cannot reuse the
spherical residual (which carries the :math:`\partial_\mu` term the
slab does not have).

**Multi-group generalisation (T2).** The slab case is multi-group-
capable. Each group carries a per-group amplitude :math:`c_g` scaling
the shared shape, :math:`A_g(x) = c_g(a_0 + a_1\sin kx)` and
:math:`B_g(x) = c_g\,b_0\cos kx`, and the source picks up the
in-scatter term

.. math::

   Q^{\text{ext}}_{n,g}(x) = \mu_n A_g'(x) + \mu_n^2 B_g'(x)
       + \Sigma_{t,g}\,A_g(x) + \Sigma_{t,g}\,\mu_n B_g(x)
       - \sum_{g'} \Sigma_s[g', g]\,A_{g'}(x).

The in-scatter sum uses the ORPHEUS scattering convention
``SigS[g_from, g_to]``, so the in-scatter source is
:math:`(\Sigma_s^\top\phi)_g = \sum_{g'}\Sigma_s[g', g]\,A_{g'}` — the
**transpose-active** term where the ERR-002 group-swap hazard lives.
T2 uses a 2-group **asymmetric downscatter-only** :math:`\Sigma_s`
(:math:`\Sigma_s[0,1]\neq 0`, :math:`\Sigma_s[1,0]=0`) so a transposed
scattering matrix would produce a detectably wrong group ratio
(the 1-group-degeneracy rule — multi-group with asymmetric :math:`\Sigma_s` is
mandatory, ``vv-principles`` anti-pattern #3 and failure-mode #6). The
1-group T1 path is the degenerate :math:`c_{\text{groups}} = (1.0,)`
reduction of the same dataclass.

The manufactured spherical source — the Cardinal-Rule-2 reuse
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The spherical ansatz (pole-regular, non-vacuum at :math:`r=R`):

.. math::
   :label: sn-mms-nonvacuum-sph-psi

   \psi_n(r) = \frac{1}{W}\bigl(A(r) + \mu_n B(r)\bigr),
   \quad A(r) = a_0 + a_1\sin(kr),\quad
   B(r) = \frac{r}{R}\bigl[b_0 + b_1\cos(kr)\bigr],\quad B(0)=0.

.. (vv-status rationale) definition: Definitional ansatz — imposed,
.. not solved. Correctness as a reference rests on the source
.. identity :eq:`sn-mms-nonvacuum-sph-qext` (SymPy ``simplify == 0``)
.. and HAZARD H1 (:math:`B(0)=0`), verified in the foundation gate.
.. vv-status: sn-mms-nonvacuum-sph-psi documented

The spherical SN operator carries the angular-redistribution term
(:eq:`transport-spherical`):
:math:`\mu\,\partial_r\psi + \tfrac{1-\mu^2}{r}\,\partial_\mu\psi
+ \Sigma_t\psi = \tfrac{1}{W}(\Sigma_s\phi + Q^{\text{ext}})`.
Substituting :eq:`sn-mms-nonvacuum-sph-psi` (with
:math:`\partial_\mu\psi = B/W`, so the redistribution term is
:math:`\tfrac{1-\mu^2}{r}\cdot\tfrac{B}{W}`) and isolating the source
gives the **same closed form** as the Phase 3.6 vacuum case
(:eq:`sn-mms-spherical-aniso-qext`):

.. math::
   :label: sn-mms-nonvacuum-sph-qext

   Q^{\text{ext}}_n(r) = \mu_n A'(r) + \mu_n^2 B'(r)
                       + (1-\mu_n^2)\,\frac{B(r)}{r}
                       + (\Sigma_t-\Sigma_s) A(r)
                       + \Sigma_t\,\mu_n B(r).

.. (vv-status rationale) derivation: A closed form obtained by
.. symbolic substitution into the continuous spherical operator;
.. verified by SymPy ``simplify == 0`` (reusing the 3.6 spherical
.. residual machinery) and cross-checked against the Branch-2 numpy
.. evaluation to ≤1e-13.
.. vv-status: sn-mms-nonvacuum-sph-qext verified

The structural point is that :eq:`sn-mms-nonvacuum-sph-qext` and the
Phase 3.6 :eq:`sn-mms-spherical-aniso-qext` are *byte-identical* closed
forms — only the radial profiles :math:`A,B` plugged into them differ.
The spherical-operator residual is therefore derived in **exactly one
place**:
:func:`~orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`,
which now takes optional ``A=None, B=None`` arguments. With no
arguments it reproduces the Phase 3.6 vacuum shapes byte-for-byte (the
decision-A regression pin verifies this in
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_existing_spherical_aniso_still_passes_after_parameterization`);
with the 4.6 non-vacuum shapes it re-proves the residual for free
(:func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_spherical_mms`).
This is Cardinal Rule 2 in action — one source of truth for the
spherical transport-operator residual, shared between the vacuum and
non-vacuum cases.

HAZARD H1 — sphere pole regularity demands :math:`B(0)=0`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The redistribution term :math:`(1-\mu^2)B(r)/r` has an explicit
:math:`1/r` factor. As :math:`r\to 0`, this diverges as
:math:`(1-\mu^2)\,B(0)/r \to \infty` **unless** :math:`B(0)=0`. So on
the sphere, :math:`B(0)=0` is a **hard regularity constraint** — not a
stylistic preference. A naive slab-style choice :math:`B(x) = b_0\cos
kx` gives :math:`B(0) = b_0 \neq 0`, which is **fine on the slab** (no
pole, no :math:`1/r`) but **wrong on the sphere** (it manufactures a
non-integrable :math:`1/r` singularity at the centre that the
continuous solution does not actually have).

The 4.6 sphere therefore uses :math:`B(r) = (r/R)[b_0 + b_1\cos kr]`.
The :math:`(r/R)` prefactor forces :math:`B(0) = 0` (pole-regular: the
redistribution :math:`(1-\mu^2)B/r = (1-\mu^2)[b_0+b_1\cos kr]/R` is
*finite* at :math:`r=0`), while leaving :math:`B(R) = b_0 + b_1\cos kR
\neq 0` (the non-vacuum first-moment structure at the outer inflow
face). The amplitude :math:`A(r) = a_0 + a_1\sin kr` needs **no** such
prefactor: :math:`A` has no :math:`1/r` companion in the operator, so
:math:`A(0) = a_0` finite is perfectly regular at the pole, and
:math:`a_0>0` makes :math:`A(R)\neq 0` (non-vacuum). HAZARD H1 is
verified in
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_sph_pole_regularity_and_nonvacuum`
(:math:`B(0)=0`, :math:`A(0)=\tfrac12`, :math:`B(R)\neq 0`; concretely
at :math:`kR=\pi/2`: :math:`A(R)=\tfrac34`, :math:`B(R)=\tfrac{3}{10}`).

The non-vacuum lever — :math:`a_0>0` is the entire novelty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every other choice in 4.6 (the P1 angular form, the interior operator,
the redistribution term, the quadrature, the BC machinery) is shared
with Phase 3.6. The *single* new ingredient is :math:`a_0 > 0`, which
makes :math:`A` — and hence the inflow trace
:math:`\gamma_-\psi_n = (A + \mu_n B)/W` — **non-zero at the outer
face**. That non-zero trace is what lights up the prescribed-inflow
``q.boundary`` path. Strip :math:`a_0` back to zero and 4.6 degenerates
to Phase 3.6 (vacuum-automatic). The non-vacuum-ness is pinned by the
foundation test
:func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_v_nonvac_slab_ansatz_nonvanishing_at_faces`
(:math:`A(0)=a_0>0`) so the verification cannot silently drift back to
the vacuum regime.

The affine-BC-to-RHS framing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Prescribed inflow is **not** a special solver mode — it is the
inhomogeneous term of the universal affine boundary law
(:ref:`affine-bc-form`). The general boundary trace law is

.. math::

   \gamma_-\psi = R\,G\,\gamma_+\psi + q,

where :math:`q \in \Gamma_-` is the **prescribed inflow source**. For a
vacuum boundary :math:`R=0` and :math:`q\equiv 0`. For a manufactured
non-vacuum inflow, :math:`q = \gamma_-\psi_{\text{chosen}}` — the
imposed inflow trace — pushed to the right-hand side of the
discretised within-group system.

In ORPHEUS this :math:`q` is carried by the ``q.boundary`` slot, a
:class:`~orpheus.transport.source_sinks.BoundarySourceSink` field whose
inflow-ordinate entries hold :math:`\gamma_-\psi = (A + \mu_n B)/W` per
face per group. The within-group fixed point is the **affine** system

.. math::

   (L + C - S - B)\,\psi = q,
   \qquad q = q_{\text{ext}}
            + (\text{prescribed inflow in } q.\text{boundary}),

and the inflow term is consumed directly by :math:`(L+C)\text{.solve}`
as the sweep inflow seed. This is the cleanest possible realisation:
the inhomogeneous BC term is *just another source* on the RHS.

**No ``from_spec`` / ``PrescribedInflow``-BC bridge is touched.** A
:class:`~orpheus.geometry.boundary.PrescribedInflow` mesh-BC descriptor
*does* exist (the rank-0 affine BC), but it is a *different surface* —
it declares a prescribed inflow at mesh-construction time as a
first-class boundary condition. The 4.6 MMS deliberately does **not**
use it: the manufactured inflow is injected as the ``q.boundary``
source slot, which is exactly the affine-:math:`q`-to-RHS path, and is
the surface O.2b's field-role-typing work targets. The mesh BCs for the
4.6 cases are plain **vacuum** — the inflow lives entirely in
``q.boundary``.

The public composite-source API drives the convergence rows
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 4.6 convergence rows drive the **public** fixed-source entry point
:func:`~orpheus.sn.solver.solve_sn_fixed_source` directly. Earlier in
this work that entry point hardcoded a vacuum ``q.boundary``
(``zeros_on``) — it had no way to carry a prescribed inflow — and the
rows therefore took an operator-triple *bypass*: assembling
:math:`(L+C)`, :math:`S`, :math:`B` by hand via
:func:`~orpheus.sn.solver._within_group_triple` and driving them with
:func:`~orpheus.numerics.iteration.SourceIteration`. That bypass is
**retired**. The field-role-typing work gave
:func:`~orpheus.sn.solver.solve_sn_fixed_source` a second accepted
source form — the full **composite source** ``q = q_bulk ⊕ q_∂``
represented by a
:class:`~orpheus.transport.timed_full_field.TimedFullField` (see
:ref:`sn-composite-fixed-source` for the API in full). Each case now
bundles its manufactured bulk (:meth:`external_source`) and its
prescribed-inflow boundary (:meth:`prescribed_inflow`) into one
``case.fixed_source(sn)`` and passes it straight to the public solver::

    result = solve_sn_fixed_source(
        materials, mesh, case.quadrature, case.fixed_source(sn),
        max_inner=1000, inner_tol=1e-13,
    )

Migrating the rows off the bypass onto the public API *is* the
retirement (retirement = test migration — the new code is what gets
tested). The slab rows take the SI (1-D Jacobi) inner; the sphere row
takes the curvilinear **source-iteration** default (post-ERR-058, #195;
SI :math:`\equiv` Krylov bit-identical — the curvilinear ``"krylov"``
default was reverted, see :ref:`sn-err-058-closure-seed-closeout`); both
honour the prescribed inflow identically (verified by T4 below).

**The flux/source space bridge — now INTERNAL to the solve (B.5.2).**
The composite RHS lives in **source** space (an
:class:`~orpheus.transport.source_sinks.AngularSourceSink` bulk plus a
:class:`~orpheus.transport.source_sinks.BoundarySourceSink` boundary),
while the iterate :math:`\psi` and the returned solution live in
**flux** space (an :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
bulk plus a :class:`~orpheus.transport.fields.boundary_flux.BoundaryFlux`
boundary). The source-iteration / Krylov inner therefore needs a
flux-typed ``initial_guess`` to template the solution space — without
it, ``S.apply`` would hit an ``AngularSourceSink`` that has no
``integrate_angular`` method. That seed
(``TimedFullField.zeros(bulk=AngularFlux, boundary=BoundaryFlux,
mesh=sn)``) is now built **inside**
:func:`~orpheus.sn.solver.solve_sn_fixed_source`, not hand-passed by
the test. The field-role-typing distinction — the iterate is a *flux*,
the RHS is a *source* — survives intact; it has simply moved behind the
public API where it belongs.

The converged-value assertion — rate is necessary, not sufficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the load-bearing verification design choice, and it is a
direct application of ``vv-principles`` anti-pattern #5 ("NEVER read
'convergence rate is correct' as 'result is correct' — verify the
converged-to value; :math:`\mathcal{O}(h^2)` to the wrong limit is
still :math:`\mathcal{O}(h^2)`") and the necessity hierarchy H4.

Consider the failure mode the test must catch: a bug (or a refactor
regression) that silently **drops the prescribed inflow** — a solve
that runs with ``q.boundary = 0`` despite the manufactured non-vacuum
inflow. That degenerate solve is **still a perfectly consistent
fixed-source problem** — it just solves the *vacuum-BC* version of the
same interior source. It converges cleanly at :math:`\mathcal{O}(h^2)`
to a *different*, boundary-zero scalar flux. A rate-only test passes
it. The only assertion that sees the dropped inflow is the one that
checks the **converged value against** :math:`A(x)`, because
:math:`A(x)` is non-zero at the boundary (:math:`a_0>0`) while the
dropped-inflow limit is zero there — a discrepancy of order
:math:`a_0 \approx 0.5` at the faces, dwarfing the pointwise
convergence error (~8e-5).

The slab T1/T2 rows therefore make **three** assertions per group:
(1) the rate ``orders > 1.9`` (DD design order on a smooth ansatz);
(2) the finest-mesh :math:`\phi_{\text{num}}` matches :math:`A(x)` to
``rtol=atol=5e-3`` — with a guard asserting the reference is genuinely
non-vacuum (:math:`|A(0)|, |A(L)| > 0.1`) so the value check is
discriminating; and (3) an inflow-honoured spot-check that the solved
trace slot equals the imposed :math:`\gamma_-\psi = (A + \mu_n B)/W` to
``rtol=1e-9``. Only the combination is a meaningful test of the
prescribed-inflow path.

The Mode-7 activates/nulls map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``vv-principles`` failure-mode #7 ("MMS simplification bias") requires
every multi-dimensional MMS test to **declare** which operator terms
its ansatz activates and which it nulls — and to ship an
angularly-non-trivial companion whenever the nulled set includes a term
covered by an active ERR-NNN. The 4.6 declaration:

.. list-table:: Mode-7 term map — slab vs sphere under the (A+μB)/W ansatz
   :header-rows: 1
   :widths: 40 30 30

   * - Operator term
     - Slab (Cartesian)
     - Sphere (spherical)
   * - streaming :math:`\mu A'` (isotropic)
     - **activates**
     - **activates**
   * - streaming :math:`\mu^2 B'` (second moment)
     - **activates**
     - **activates**
   * - angular redistribution :math:`(1-\mu^2)B/r`
     - **nulls** (no :math:`\partial_\mu` term)
     - **activates** (the ERR-026 path)
   * - within-group scatter :math:`\Sigma_s\phi` (:math:`c<1`)
     - **activates**
     - **activates**
   * - collision :math:`\Sigma_t\,\mu B` (first moment)
     - **activates**
     - **activates**
   * - 2G group transfer :math:`\Sigma_s^\top` (asymmetric)
     - **activates** (T2)
     - n/a (1g)
   * - prescribed non-vacuum inflow :math:`\gamma_-\psi \neq 0`
     - **activates** (both faces, :math:`a_0>0`)
     - **activates** (:math:`r=R` face)
   * - fission
     - **nulls** (non-fissile; MMS proves no eigenvalue)
     - **nulls**

The slab **nulls the angular redistribution** — the Cartesian operator
has no :math:`\partial_\mu` coupling. Redistribution is exactly where
ERR-026 (the curvilinear sweep WDD wrong-fixed-point bug) lives, so a
slab-only 4.6 would be a textbook Mode-7 trap: it would verify the
prescribed-inflow path while being structurally blind to the hardest
math the curvilinear sweep performs. The **sphere companion is
therefore mandatory** (NEVER ship slab-only — ERR-026 territory). The
sphere activates the redistribution term under non-vacuum inflow,
closing the Mode-7 declaration.

T3 (sphere) — why ``xfail(strict)``, now on Issue #229
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The sphere row
(:func:`tests.sn.verification.analytical.test_mms_prescribed_inflow.test_mms_prescribed_inflow_sphere_activates_redistribution`)
ships ``@pytest.mark.xfail(strict=True)`` with ``catches("ERR-026")``.
The reason is **not** that the non-vacuum machinery fails.

.. note:: **Re-scoped (2026-06-12, Issue #195 CLOSED → Issue #229).**

   This row's xfail was originally attributed to the #195
   "pre-asymptotic transient" plateau (the stagnation table below).
   The ERR-058 closure-seed fix **closed** the curvilinear
   wrong-fixed-point family, so the stagnation is gone — the isotropic
   curvilinear DD interior is now :math:`\mathcal{O}(h^2)`-consistent
   (see :ref:`sn-err-058-closure-seed-closeout`).  T3's ansatz, however,
   is the *anisotropic* :math:`(A(r)+B(r)\mu)/W`, which is angle-varying
   and therefore hits the **fixed-quadrature angular floor** of the
   half-angle thread interpolation (sphere S16 floor
   :math:`\approx 7\mathrm{e}{-4}`, above the band).  The marker
   **stays ``xfail(strict)``** but is now pinned to the
   `Issue #229 <https://github.com/deOliveira-R/ORPHEUS/issues/229>`_
   quadrature-aware retune, NOT the #195 plateau; it flips to
   unexpected-pass when #229 lands (the regression gate for the retune).
   The stagnation table below is preserved as **bug-era evidence**; its
   "pre-asymptotic" interpretation is superseded.

**Bug-era stagnation (pre-ERR-058).**  The slab was pole-free and
converged perfectly: orders ``[2.04, 2.01]``, finest L2 ~1.2e-3,
pointwise ``max|φ−A|`` ~8e-5, boundary value matched.  The sphere L2
(volume-weighted), by contrast, **stagnated** mesh-independently — the
plateau that refuted the "pre-asymptotic transient" premise:

.. list-table:: T3 sphere volume-weighted L2 error (bug-era plateau, pre-ERR-058)
   :header-rows: 1
   :widths: 25 25 25 25

   * - :math:`n_c`
     - 20
     - 40
     - 80
   * - :math:`\|\phi_h - A\|_{L^2(V)}`
     - 2.37e-2
     - 2.42e-2
     - 2.43e-2

The observed "orders" were ≈ :math:`-0.02` to :math:`-0.006` — the
error was *not* decreasing under refinement, the plateau signature
ERR-058 diagnosed.  Post-ERR-058 the spatial convergence is recovered;
the residual gap on this *anisotropic* row is the #229 angular floor,
which DOES drop under quadrature refinement (sphere S16
:math:`\to` S32 halves the floor) — the structural test that the
remaining gap is angular, not a wrong fixed point.

Crucially, the **boundary value is honoured** (always was): the
finest-mesh :math:`\phi[-1] \approx 0.7499 \approx A(R) = 0.75`, and
the inflow-trace spot check passes.  The non-vacuum prescribed-inflow
machinery *works*; the remaining xfail is purely the angular-floor
budget of the anisotropic ansatz under fixed quadrature.

**T3g — the green structural companion.** Because T3 is xfail (now on
#229), it provides *no live* convergence coverage of the 4.6 machinery.
The green companion
:func:`tests.sn.verification.analytical.test_mms_prescribed_inflow.test_sphere_nonvacuum_inflow_honoured_and_redistribution_live`
fills that gap with two non-convergence-dependent claims that pass
*now*: (1) the prescribed inflow at :math:`r=R` is honoured per inflow
ordinate (:math:`\gamma_-\psi = (A(R) + \mu_n B(R))/W` with :math:`A(R)
> 0` non-vacuum); and (2) the redistribution source
:math:`(1-\mu^2)B(r)/r` is non-zero on the mesh interior (the ERR-026
term is live under the 4.6 ansatz — :math:`B(r)\neq 0` on the open
interval, with :math:`B(0)=0` pole-regular). T3g is the live structural
guarantee that the Mode-7 sphere companion exercises the redistribution
path even while the convergence row is parked on #229.

T4 (vv Mode 9) — splitting invariance of the prescribed inflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The consistency floor the convergence rows trust is that a non-zero
prescribed inflow is honoured **identically** by the three operator
splittings of the affine within-group system: SI-Jacobi (the resolvent
:math:`L+C` with lagged gains :math:`S, B`), SI-Gauss–Seidel (the
:math:`B`-folding ``_GaussSeidelResolvent``), and Krylov (the matvec
:math:`L+C-S-B`). All three are different reduction trees of the *same*
affine fixed point :math:`(L+C-S-B)\psi = q`, so they MUST reach the
same :math:`\psi` (``vv-principles`` Mode 9 — verify splittings reach
the same fixed point under anisotropic / :math:`B\neq 0` stressing).
This is a **foundation** test, not an L1 claim: no theory-page
:math:`:label:` is being verified — it pins that three reduction trees
of one affine operator agree on one RHS, which is a software invariant
(``foundation`` NEVER carries ``verifies()``).

T4
(:func:`tests.sn.verification.analytical.test_prescribed_inflow_consistency.test_prescribed_inflow_consistency_si_jacobi_gs_krylov`)
runs two configs. The ``slab_1d`` config (SI is always Jacobi in 1-D)
makes **SI ≡ Krylov** the discriminating pair. The
``cart2d_reflective_y`` config adds reflective-:math:`y` faces so
:math:`B \neq 0` — which is what makes **SI-Jacobi vs SI-Gauss–Seidel**
distinct (G-S folds :math:`B` into the resolvent; Jacobi lags it). The
:math:`B\neq 0`-plus-prescribed-inflow combination is the only config
where the :math:`B`-folding path runs *with* a non-zero boundary source
(the ERR-056 neighbourhood). Measured pairwise reldiffs: 1.3e-13 …
5.6e-13 — comfortably under the 1e-11 ceiling, which itself leaves
headroom for the FP-non-associativity of three reduction trees
(bounded by :math:`\text{iter} \times \text{ULP}` per the
``vv-principles`` bit-identity criteria).

The test carries explicit anti-latent-dud preconditions (the
splitting-invariance check is vacuous if all three trivially agree on
:math:`\psi \equiv 0`): the inflow slot must actually be written
(:math:`>0`), the inflow must non-trivially drive the interior
(:math:`\max|\psi| > 10^{-3}`), and the 2-D row must actually select the
:math:`B`-folding ``_GaussSeidelResolvent`` (not silently fall back to
Jacobi) with an explicit reflective-:math:`y` ``Mesh2D`` BC.

Verification chain — Branch 1 / Branch 2 / L1 cross-check
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following the ``algebra-of-record`` discipline:

- **Branch 1 (SymPy, State 1C — MMS).** The manufactured sources
  :eq:`sn-mms-nonvacuum-qext` (slab) and
  :eq:`sn-mms-nonvacuum-sph-qext` (sphere) are derived by substituting
  the imposed ansatz into the continuous operator and solving for the
  residual, symbolically. The slab pair is
  :func:`~orpheus.derivations.continuous.mms.sn._nonvacuum_slab_symbolic`
  /
  :func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_slab_mms`;
  the sphere reuses
  :func:`~orpheus.derivations.continuous.mms.sn._spherical_anisotropic_symbolic`
  with the 4.6 shapes
  (:func:`~orpheus.derivations.continuous.mms.sn._nonvacuum_spherical_AB`)
  via
  :func:`~orpheus.derivations.continuous.mms.sn.derive_nonvacuum_spherical_mms`.
  Each ``derive_*`` proves ``simplify(W·LHS − Σ_s·φ − Q_closed) == 0``.
  Foundation gate:
  :file:`tests/derivations/test_sn_mms_nonvacuum_symbolic.py`.
- **Branch 2 (vectorised numpy).** The factories
  :class:`~orpheus.derivations.continuous.mms.sn.SNSlabNonVacuumMMSCase`
  and
  :class:`~orpheus.derivations.continuous.mms.sn.SNSphericalNonVacuumMMSCase`
  (built by
  :func:`~orpheus.derivations.continuous.mms.sn.build_slab_nonvacuum_mms_case`,
  :func:`~orpheus.derivations.continuous.mms.sn.build_slab_2g_nonvacuum_mms_case`,
  and
  :func:`~orpheus.derivations.continuous.mms.sn.build_sphere_nonvacuum_mms_case`)
  evaluate the closed-form source per ordinate using vectorised numpy.
  Each carries a ``prescribed_inflow(sn)`` method returning the
  ``q.boundary`` :class:`~orpheus.transport.source_sinks.BoundarySourceSink`,
  and a ``fixed_source(sn)`` bundler returning the composite
  ``q = q_bulk ⊕ q_∂``
  :class:`~orpheus.transport.timed_full_field.TimedFullField` the public
  solver consumes (see :ref:`sn-composite-fixed-source`).
- **L1 cross-check (the gate).** The Branch-2 numpy
  :math:`Q^{\text{ext}}_n` is bit-equal (≤1e-13 max absolute) to the
  Branch-1 SymPy closed form evaluated via :func:`sympy.lambdify` on a
  sample mesh, for both geometries. The two branches are *structurally
  independent above the trusted-library line* — Branch 1 is
  ``lambdify``-d SymPy, Branch 2 is hand-written numpy — so agreement
  catches a copy error between the symbolic derivation and the
  numerical implementation. Tested in
  :func:`tests.derivations.test_sn_mms_nonvacuum_symbolic.test_slab_nonvacuum_numerical_qext_matches_sympy`
  and the spherical sibling.

**Structural independence (L11).** The chosen scalar flux
:math:`\phi = A` is *imposed* analytically; the source :math:`Q^{\text{ext}}`
is SymPy-derived (not generated by the solver's own primitives); the
numpy ``external_source`` is then cross-checked bit-equal to the
lambdified SymPy. The reference is structurally independent of the code
under test — the manufactured source does not pass through any of the
solver's discretisation primitives, so the L1 convergence rows are a
genuine test, not a tautology.

**What this section does NOT verify.** Per the three pillars
(``vv-principles``), MMS is a *source-driven* problem: it verifies the
convergence order (a math claim) and the flux shape (a model claim,
because the source is structurally independent), but it **cannot**
verify an eigenvalue. The 4.6 mixtures are non-fissile by construction
and there is no eigenvalue claim anywhere in this section. The
prescribed-inflow verification is a forward-only, fixed-source result.


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
:math:`\xi` (``mu_y``).  The correct recursion ([BaileyMorelChang2010]_,
the curvilinear :math:`\alpha`-recursion)
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

Applied the [BaileyMorelChang2010]_ formulation:

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
``_sweep_1d_cumprod`` (the dissolved ``sweep.py``), the 1D Cartesian
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
wavefront sweep :func:`~orpheus.sn.loss_representation._sweep_jacobi`
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
``_sweep_1d_cumprod`` (the dissolved ``sweep.py``) plus a comment block
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
:func:`~orpheus.sn.loss_representation._sweep_jacobi`,
``_sweep_1d_spherical`` (the dissolved ``sweep.py``),
``_sweep_1d_cylindrical`` (the dissolved ``sweep.py``), and
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


.. _sn-development-history:

Development history
===================

S\ :sub:`N` is ORPHEUS's **architectural prototype**: the typed-field
algebra and the composable operator form :math:`(L + C - S - F/k)\psi = q`
were developed here *first*, and are the standard the other solvers (CP,
MoC, MC, diffusion) inherit as they are built — which is why this page
carries far more development history than any other theory chapter.

This is a reverse-chronological (latest first) changelog of the major
**architectural** milestones. Iteration-rate work, gate counts, and
intermediate replans are deliberately omitted — see the GitHub issues
and the per-phase plan files for that granularity. Each entry names the
architectural change, its issue, and the commit/merge where the work
lives. Entries marked *(in development)* live on an unmerged feature
branch and have no landed hash yet.

.. list-table::
   :header-rows: 1
   :widths: 8 52 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - in dev
     - **Spatial ⊗ angular product** — τ becomes closure-owned
       (``morel_montry_tau_per_level``); the geometry-τ producers and
       the ``StreamingTerms`` τ/α fields are deleted, and a separability
       MMS gate pins Cartesian-separates / curvilinear-couples. The
       space-angle tensor structure made explicit.
     - #236
     - *(in development)*
       ``feature/sn-spatial-angular-product``
   * - 2026-06
     - **Field-typed operator algebra + data XS layer** — cross sections
       become ``CoefficientField`` leaves; operators become *promotions*
       (``C = M[\sigma_t]``); the carrier is the timeless ``FullField``;
       the §5.6 Operator / Kernel / Functional taxonomy and the
       ``@overload`` "Pattern M" apply-fibration. Closes the typed-field
       campaign and pushes the data layer down to validated invariants
       (χ-simplex, production-weighted χ\ :sub:`mix`, XS-balance).
     - #257
     - ``505e1b7`` → ``f62b895``
   * - 2026-06
     - **Linear-Discontinuous on the DAG + principled polymorphism** —
       LD lands as a polymorphic discretization protocol (the
       coefficient model), the ``matvec_via_kernel`` favoritism is
       reverted, σ is single-sourced from the ``InvertibleOperator``'s
       ``C``, and the unified all-d LD moment matvec + diffusion-limit
       closure ships.
     - #240 / #158
     - ``fde76ac`` → ``e74eafb``
   * - 2026-06
     - **Sweep / matvec re-layering — one walk** — the sweep strategy
       is carved into a first-class ``LossRepresentation``; the
       ScanMarch row-march twin makes *matvec ≡ sweep* literally one
       ``_OctantWalk`` traversal rather than two parallel paths.
     - #222
     - ``8913229`` → ``1b4b0c0``
   * - 2026-06
     - **3-D Cartesian admission** — axis-native ``from_axes`` admits
       *d = 3* end-to-end with no ``Mesh3D``; the constructor data-flow
       inverts so axes are primary, and windowing / G-S gates key on
       genuine dimensionality rather than the reduced-proxy.
     - #225
     - ``1da1e2f``
   * - 2026-06
     - **Foundation-cleanup cluster** — break the numerics→SN import
       cycle (moment-layout policy homed in numerics), retire the
       ``M_spatial`` / ``M_angular_redist`` operator-leaf split (the
       fused ``loss_action`` is the only path), and key the moment-frame
       involution on intent rather than coincidental shape.
     - #243 / #245 / #246 / #238
     - ``66b1bbd`` → ``dd4f542``
   * - 2026-06
     - **LD boundary slope source + transverse face-slope trace** — the
       fixed-source solver consumes a moment-resolved external slope
       source (Leg A) and the boundary trace carries the transverse
       face-slope moment (Leg B); both close the LM-1989 trap as
       structural-teeth vv Mode-10 cases.
     - #247 / #251
     - ``d9396a2`` / ``e5f2b1c``
   * - 2026-06
     - **Wave-O affine-typed operator algebra** — ``FluxDisplacement``
       gives the affine difference space :math:`V`; ``flux + flux``
       becomes a ``TypeError`` while ``flux − flux`` and the typed
       residual ``from_balance`` stay legal, so :math:`(L+C-S-F)\psi = q`
       types coherently and the residual is a typed defect.
     - #208 / #201
     - ``8c2f355`` / ``04e2859``
   * - 2026-06
     - **G-adjoint, reciprocity, and operator role-typing** — the
       analytic reverse-sweep ``StreamingOperator.apply_transpose``; the
       Hilbert-adjoint metric is owned by the ``FunctionSpace``
       (``apply_metric``); composite block-roles are *derived* from the
       operands, retiring the ``InvertibleOperator`` FULL stamp.
     - #20
     - ``0efd233`` → ``7ccc14a``
   * - 2026-06
     - **Angular + face windowing** — the held SI iterate becomes
       ``HarmonicMomentField`` moments (angular window), the interior
       face cochain becomes a rolling ``_MovingFrontier`` (face window),
       and moments are accumulated in-sweep — eliminating the
       full-angular per-sweep transient (a measured ~3× peak-memory win).
       Full-field sweep/matvec oracles retained for the equivalence gate.
     - #205 / #218
     - ``b97d4f9`` → ``c7be111``
   * - 2026-06
     - **Honest variadic L+C−S−B driver** — the transitional ``S+B``
       fold is retired and boundary reflection ``B`` becomes a
       first-class coupling gain, so the within-group drivers generalize
       over a variadic operator list instead of a hard-coded fold.
     - #18
     - ``8563f4b`` → ``83a4ae6``
   * - 2026-06
     - **Source-iteration boundary Gauss-Seidel** — a polymorphic
       ``SweepSchedule`` (Jacobi / octant-group G-S) plus the
       ``_GaussSeidelResolvent`` and ``inner_schedule`` selector; a
       modest reflective-SI accelerator (the dominant scattering c-mode
       stays Krylov / DSA territory). Surfaced ERR-056 (diagonal-cubature
       shared-face fan-in).
     - #19
     - ``514ae21`` → ``a39905a``
   * - 2026-06
     - **1-D forward + transpose matvec carve** — the 1-D SN matvec is
       relocated off the operator into ``_OneDimScanWalk``, making
       *matvec ≡ sweep* a single-sourced code fact for the 1-D path.
     - #206
     - ``eaafbe1`` / ``7300f3e``
   * - 2026-06
     - **Curvilinear closure-seed fix (ERR-058) and eigenvalue
       verification** — coupled-pole spatial + half-angle angular
       closure seeds restore :math:`O(h^2)` for curvilinear sweeps, and
       the SI ≡ Krylov eigenvalue equivalence is pinned, closing the
       last ERR-026 manifestation.
     - #195 / #196
     - ``3b088ee`` → ``8609282``
   * - 2026-06
     - **2-D matvec through the DD closure** — the 2-D forward matvec is
       routed through the diamond-difference cell closure and the legacy
       FD stencil is retired, so the operator and the sweep share one
       discretization.
     - #208 (O.4b)
     - ``2288ea4``
   * - 2026-05
     - **Four-operator typed apply + typed-field foundation** — the
       ``F → S → C → L`` apply overload on ``AngularFlux``, founded on
       the typed ``AngularFlux`` / ``ScalarFlux`` / ``BoundaryFlux``
       fields (``psi_bc`` retired). The birth of the typed-field
       architecture this whole history builds on.
     - #197
     - ``d8ddb03`` → ``eeac45f``
   * - 2026-05
     - **Depth-B field factoring** — ``ScalarFlux`` migrates to a pure
       ``Field`` subtype and the bare-``ndarray`` operator arms are
       retired (D-I / D-J), so typed-field dispatch becomes the *only*
       apply path through the operators.
     - #197 (Depth B)
     - ``c97897d`` → ``4a53737``

.. note::

   The four ``Phase N`` milestones (#18, #19, #20, #205) carry internal
   phase labels in their commit subjects rather than GitHub-issue
   trailers; the issue mapping above is from the SN development-sequence
   campaign record. Issue #236 is the only entry not yet on ``main`` — it
   lives on ``feature/sn-spatial-angular-product`` pending merge.


References
==========

.. _bib-bailey-morel-chang-2010:

.. [BaileyMorelChang2010] T.S. Bailey, J.E. Morel, and J.H. Chang,
   "The asymptotic diffusion-limit accuracy of S\ :sub:`N` angular
   differencing schemes," *Nuclear Science and Engineering*,
   165(2):149--169, 2010 (LLNL preprint LLNL-JRNL-420356; OA at
   https://www.osti.gov/servlets/purl/1020346).  **Eq. 43** gives the
   Morel--Montry weight :math:`\tau_n = (\mu_n - \mu_{n-1/2}) /
   (\mu_{n+1/2} - \mu_{n-1/2})` as the unique weight exact for a flux
   linear in :math:`\mu`, with admissible range :math:`\tau \in [0,1]`
   — the W1 source for dropping the spherical clamp
   (:ref:`sn-tau-clamp-vindication`).  The paper's diffusion-limit
   analysis keeps :math:`r` continuous (it deliberately removes spatial
   differencing to isolate the angular error), so it cannot speak to the
   spatial pole-cell order of :ref:`sn-pole-cell-spatial-closure`.

.. [WuXieFischer1999] G. Wu, Z. Wu, and B. Fischer,
   "A discrete ordinates nodal method for one-dimensional neutron
   transport calculation in curvilinear geometries,"
   *Nuclear Science and Engineering*, 133, 1999
   (DOI 10.13182/NSE99-A2095).  Green's-function nodal S\ :sub:`N` with
   Legendre spatial expansion — "very high precision on coarse spatial
   meshes relative to the standard fine-mesh S\ :sub:`N` method with the
   spatial diamond-differencing scheme."  The documented route to
   better-than-diamond central-cell *spatial* order; cited at
   :ref:`sn-pole-cell-spatial-closure` as the class of fix that would
   lift the inherent pole-cell :math:`\mathcal{O}(h)`.

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
   that the Wave E Round 2 inner solver implements.  §II gives the
   source-iteration spectral radius :math:`\rho = c`; §VI reviews the
   KBA / wavefront parallel-sweep ordering whose fan-in discipline the
   octant-group Gauss-Seidel schedule inherits.

.. [Pautz2002] S.D. Pautz,
   "An algorithm for parallel S\ :sub:`n` sweeps on unstructured
   meshes," *Nuclear Science and Engineering*, 140(2):111--136, 2002.
   The KBA-style wavefront octant-scheduling reference.  Cited at
   :ref:`si-within-group-splitting` for the shared-face fan-in rule
   (ERR-056): a boundary face outflowed by several octants must be
   reduced (reflected) only after the LAST contributing octant group
   has swept it.

.. [Pomraning1989] G.C. Pomraning,
   "The transport equation in general geometry,"
   *Nuclear Science and Engineering*, 101:330--340, 1989.
   Page 339 frames the curvilinear pole singularity as **structural**:
   :math:`r = 0` is intrinsically singular in any curvilinear streaming
   operator because the angular-derivative coefficients contain
   :math:`1/r`.  Cited at
   :ref:`sn-phase-d-pomraning-structural-singularity` as the deeper
   reason :math:`\mu = \pm 1` is the only admissible Carlson starting
   direction.
