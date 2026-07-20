.. _theory-discrete-ordinates:

==========================================
Discrete Ordinates Method (S\ :sub:`N`)
==========================================

.. contents:: Contents
   :local:
   :depth: 3


.. Machine header — the ``nexus-meta`` schema for this page.  Ingestion is
.. PENDING nexus#1 Phase 2: the ``nexus-meta`` directive is NOT yet
.. registered, so the schema is rendered here as a collapsed sphinx-design
.. dropdown and machine-consumed later.

.. dropdown:: Machine header — ``nexus-meta`` schema (module · operators · conventions · invariants)
   :color: muted

   .. code-block:: yaml

      module: sn
      method: discrete-ordinates
      aliases: [SN, discrete ordinates, Sₙ, transport sweep, ordinate method]
      governing_equation: "A ψ = (1/k) F ψ  [eigenvalue] ;  A ψ = q  [fixed source]"
      operators:
        L: streaming, BULK only (Ω·∇) — the boundary law is the sibling B, NOT folded into L
        C: collision / removal (Σ_t)
        S: scattering in-scatter gain (Σ_s0ᵀ φ + anisotropic moments)
        B: boundary law as a first-class SIBLING operator (reflective / vacuum / white trace), every geometry
        F: fission production (χ ⊗ νΣ_f, rank-1 dyad)
      composites:
        A: "L + C - S - B — the within-group loss operator; the Krylov driver applies it"
        (L+C): "lower-triangular under the upwind cell ordering; (L+C)⁻¹ IS the transport sweep"
      key_types: [AngularFlux, SNMesh, HarmonicMomentField, SweepDependencyGraph]
      entry_points:                    # qualnames; Nexus links via implements edges
        - orpheus.sn.solver.solve_sn
        - orpheus.sn.solver.SNSolver
      conventions:
        sign: "μ>0 is +x / outward-radial outflow; inward inflow at the left / outer boundary"
        scattering: "Mixture.SigS[l][g_from, g_to]; the in-scatter source uses the TRANSPOSE  Q = SigSᵀ @ φ"
        diamond_difference: "ψᵃ = (1+β)ψ_out − β ψ_in; Morel–Montry sets β = 0 (Bailey–Morel–Chang 2010 Eq. 43, unique exact-on-linear-in-μ)"
        quadrature_norm: "GL weights sum to 2; Lebedev / level-symmetric / product sum to 4π; moments carry NO 4π prefactor (Σw normalisation)"
        layout: "angular flux (N, ng, nx, ny); scalar flux and per-cell XS (ng, nx, ny); 1-D keeps ny=1 (singleton, not squeezed)"
        group_ordering: "fast → thermal; downscatter makes SigS upper-triangular"
        starting_direction: "curvilinear half-angle seed ψ_{1/2} is first-class typed state (System B), marched directly (Issue #282 route (a)); only levels with first-ordinate raw τ ∈ (0,1) carry the block (R12a)"
      invariants:
        - "particle balance PER ORDINATE (flat-flux residual = 0) — the strong check, NOT the telescoped scalar balance"
        - "sweep ≡ matvec (one loss representation, two applications: solve vs residual)"
        - "α redistribution dome ≥ 0 (negative → NaN / overflow)"
      depends_on: [transport_methods, operator_algebra, spherical_harmonics, frame]
      verification: [L0, L1, L2]       # authored claim; cross-checked vs the Verification slice (§ below)


.. _sn-synopsis:

Synopsis
========

The discrete ordinates (S\ :sub:`N`) method solves the
:ref:`multi-group eigenvalue problem <mg-eigenvalue-problem>` in
integro-differential form by discretising the direction variable
:math:`\hat{\Omega}` into a finite ordinate set
:math:`\{(\hat{\Omega}_m, w_m)\}`, **retaining the angular flux**
:math:`\psi(\mathbf{r}, \hat{\Omega}, E)` rather than collapsing to the
scalar flux (contrast the collision-probability integral form).  It resolves
streaming, anisotropic scattering, and interface angular current directly.
ORPHEUS supports three coordinate systems under one balance framework:
**Cartesian** (slab / 2-D, no inter-ordinate coupling), **spherical** (1-D
radial, a single :math:`\alpha`-redistribution dome coupling all ordinates in
:math:`\mu`), and **cylindrical** (1-D radial, an independent :math:`\alpha`
dome per :math:`\mu`-level).  All three share a geometry factor
:math:`\Delta A / w` that guarantees per-ordinate flat-flux consistency; the
curvilinear formulation follows [BaileyMorelChang2010]_ (Eq. 43, the
Morel–Montry angular-closure weight — unique exact-on-linear-in-:math:`\mu`),
the general framework [LewisMiller1984]_, and the angular discretisation
[CaseZweifel1967]_ / [Hebert2009]_ (§3.9.4).

The solver is posed as an **operator algebra** over five operators: streaming
:math:`L` (bulk :math:`\hat{\Omega}\cdot\nabla`), collision / removal
:math:`C`, the scattering gain :math:`S`, the boundary law :math:`B` — a
first-class **sibling** operator, *not* folded into :math:`L` — and the rank-1
fission dyad :math:`F`.  They compose the within-group loss operator
:math:`A = L + C - S - B`, so the eigenvalue problem is
:math:`A\,\psi = \tfrac{1}{k}\,F\,\psi` (fixed source: :math:`A\,\psi = q`).
The sub-composite :math:`(L+C)` is lower-triangular under the upwind cell
ordering, which is exactly why :math:`(L+C)^{-1}` **is** the transport sweep
(:doc:`/theory/methods/sn/loss_representation`).  :class:`SNSolver` satisfies
the :class:`~numerics.eigenvalue.EigenvalueSolver` protocol and
:func:`solve_sn` returns an :class:`SNResult`.  Because the protocol places the
scattering source *inside* ``solve_fixed_source``, the inner source iteration
(in-scatter + anisotropic convergence) stays encapsulated in the SN sweep,
while the outer :func:`~numerics.eigenvalue.power_iteration` loop is the one
shared by CP, MoC, diffusion, and the homogeneous solver (see
:doc:`/api/numerics` for the protocol contract).

The spatial closure is **diamond difference** with the Morel–Montry weight
(:math:`\beta = 0`); the per-ordinate discrete transport is the **sweep**,
which is byte-identical to the loss-operator **matvec** — one loss
representation, two applications (``solve`` vs residual).  Both the 1-D scan
(:meth:`~orpheus.sn.loss_representation.CumprodScan.sweep`) and the 2-D
wavefront sweep
(:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph`,
per-octant batched dispatch over a mesh-time-precomputed DAG) are **bare**:
the reflective coupling :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}`
rides as a sibling :math:`-B` source term rather than a re-applied boundary
condition (:ref:`bare-sweep-extraction`, and the canonical algebra
:ref:`bc-extraction` in :doc:`/theory/foundations/boundary_conditions`).  2-D Cartesian eigenvalue
problems solve through **both** inner drivers —
:class:`~orpheus.numerics.iteration.SourceIteration` (the geometry-agnostic
default) and :class:`~orpheus.numerics.iteration.KrylovAcceleration` —
verified SI ≡ Krylov ≡ closed-form :math:`k_\infty`
(:ref:`bc-extraction-2d-si-krylov-twin`).  Interior cell-face fluxes are typed
as an interior 1-cochain :math:`C^1_{\rm int}` carried in the rolling front
(:ref:`wavefront-flux-cochain`).

Curvilinear redistribution — the geometric :math:`\alpha`-dome, distinct from
Legendre :math:`P_1^+` scattering anisotropy — and its half-angle
**starting-direction seed** are first-class typed state.  The Issue #282 route
(a) design marches :math:`\psi_{1/2}` directly from the true :math:`q_{1/2}`
source through System B's named resolvent
:meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
(a single-pass exact inverse on the *full* Legendre fold), only on levels
whose first-ordinate raw :math:`\tau \in (0,1)` (**R12a**); see
:ref:`sn-282-direct-starting-direction-solve`.  The "#229 floor" resolution
— three distinct curvilinear errors separated by a
volume-weighted-:math:`L_2`-vs-:math:`L_\infty` norm difference — is
:ref:`sn-curvilinear-aniso-norm-reconciliation`.  The storage / index /
scattering / closure **conventions** and the load-bearing **invariants** are
captured structurally in the machine header above and cross-linked below;
verification is L0–L2 plus semi-analytical (Sood, Case singular-eigenfunction)
benchmarks; the traps that hide solver bugs behind green tests are collected in
:ref:`sn-gotchas`.

.. admonition:: Conventions

   - Scattering matrix: :ref:`scattering-matrix-convention` — ``SigS[g_from, g_to]``, source uses transpose
   - **Storage layout**: :ref:`theory-sn-index-convention` — ``(N, ng, nx, ny)`` for ψ, ``(ng, nx, ny)`` for φ / σ
   - Multi-group balance: :eq:`mg-balance` in :ref:`theory-homogeneous`
   - Cross sections: :ref:`theory-cross-section-data`
   - Verification: :ref:`synthetic-xs-library` — regions A/B/C/D
   - Eigenvalue: :ref:`power-iteration-algorithm` shared with all deterministic solvers


Architecture
============

.. note:: **Implementation map — automation pending.**  The auto-generated
   Nexus filtered flow-graph figure (root symbol + traversal depth →
   graphviz) that will head this section is blocked on the nexus#20
   flow-graph directive; until it ships, the architecture below is
   **hand-authored**.  See :doc:`/api/numerics` for the live
   operator-protocol surface and :doc:`/theory/foundations/operator_algebra` for the algebra.

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
     ``alpha_half`` (angular redistribution dome), and
     ``redist_dAw`` (:math:`\Delta A_i / w_n`, shape ``(nx, N)``).
   - **Cylindrical**: ``face_areas`` (:math:`2\pi r`), ``delta_A``,
     ``alpha_per_level`` (per-level redistribution domes), and
     ``redist_dAw_per_level`` (list of ``(nx, M)`` arrays).

   The Morel--Montry angular weight :math:`\tau` is **not** a factory
   output: it is owned by the angular closure
   (:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`),
   since the geometry-side producer was retired in Issue #236 Phase 2
   Step C (see :ref:`sn-tau-c-on-cellvisit-live`).  The geometry
   factories carry **geometry only** — face areas, the
   :math:`\alpha`-dome, the redistribution factor :math:`\Delta A/w`,
   and the level starting-direction edge :math:`\mu_{1/2}`.

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

The geometry-and-quadrature dispatch is a first-class polymorphism:
:func:`~orpheus.sn.loss_representation.default_for` selects the
:class:`~orpheus.sn.loss_representation.LossRepresentation` whose declared
``supports`` predicate matches the mesh — the 1-D chain scan
(:class:`~orpheus.sn.loss_representation.CumprodScan`, any geometry) or the
multi-D anti-hyperplane wavefront — and the operator then calls it
branchlessly.  (This replaced the pre-carve procedural branch on the
``SNMesh.curvature`` string tags and the since-retired operator-free
``transport_sweep`` entry; the full carve is documented in
:doc:`/theory/methods/sn/loss_representation`.)  Boundary conditions are **not** passed as
a parameter to the sweep --- it reads the resolved BC kind strings
directly from the face-name-keyed :attr:`SNMesh.bc` dict
(``sn_mesh.bc["xmin"]``, ``sn_mesh.bc["xmax"]``, ...; see
:ref:`bc-face-name-carve`).

For 1D meshes (``ny=1``):

- **Gauss--Legendre** quadrature takes the fast 1-D chain-scan
  (:class:`~orpheus.sn.loss_representation.CumprodScan`) path (all
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

The clamp lives **inside the angular closure**, in
:func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
(sphere unclamped, cylinder clipped to :math:`[\tfrac12, 1]`).  The
closure exposes the resulting weight per global ordinate through
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
(spherical: a single :math:`(N,)` array; cylindrical: the per-level
weights gathered to the global ordinate order).  Issue #236 Phase 2
retired the parallel geometry-side :math:`\tau` producer that formerly
baked these weights onto :class:`StreamingTerms`; :math:`\tau` is now
produced **solely** by the closure (see :ref:`sn-tau-c-on-cellvisit-live`).

.. _sn-tau-closure-owned:

τ is an angular-scheme property — the closure owns it (Issue #236 Phase 2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: Archivist expansion needed.

   The Morel--Montry weight :math:`\tau` :eq:`mm-weights` is a function of
   the quadrature :math:`(\mu, w, \text{levels})` ALONE — an
   ANGULAR-scheme property, not a geometry one.  Issue #236 Phase 2
   (Step A) relocates :math:`\tau` PRODUCTION onto the angular closure:
   :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
   produces :math:`\tau` from the quadrature the closure already binds,
   and
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   consumes its OWN :math:`\tau` for the matvec contribution (P0) instead
   of reading it back from the streaming-geometry factory
   (:func:`~orpheus.geometry.reduced_operator.spherical_streaming` /
   :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming`).

   Step A is BIT-IDENTICAL: the producer is a 0-ULP line-for-line replica
   of the factory arithmetic (sphere unclamped, cylinder clamped to
   :math:`[\tfrac12, 1]`), so the geometry factory still bakes an
   IDENTICAL :math:`\tau` for the sweep path while the carve de-risks by
   parallel-run-and-compare.  The producer-equivalence gate (Leg 1)
   ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`` pins
   the closure-produced :math:`\tau` to (a) the geometry-factory value
   (0-ULP) AND (b) the structurally-independent reference
   :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
   (a different code path to the same BMC-2010-Eq.43 weight).  The
   Cartesian :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`
   supplies the neutral :math:`\tau = 1` WITHOUT a geometry branch — the
   closure TYPE is the dispatch.

   Step B (a later dispatch) retires the geometry-side
   :math:`\tau` producer and consolidates the four-site
   ``c_in``/``c_out`` duplication onto the closure.

.. _sn-closure-c-constants-owned:

c_in / c_out are angular-closure constants — Step B1 (one site folded)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: Archivist expansion needed.

   The weighted-diamond constants

   .. math::

      c_{\rm out}[m] &= \frac{\alpha_{m+1/2}}{\tau_m}, \\
      c_{\rm in}[m]  &= \frac{1-\tau_m}{\tau_m}\,\alpha_{m+1/2}
                        + \alpha_{m-1/2}

   are an ANGULAR-closure property: a function of the closure's own
   :math:`\alpha`-dome and :math:`\tau` weight :eq:`mm-weights`
   :eq:`dd-mm-closure-constants`.  Issue #236 Phase 2 Step B consolidates
   the FOUR independent inline rebuilds of this pair onto the closure,
   which already computes it once at construction (per :math:`\mu`-level,
   :math:`(M_p,)` arrays in ``_c_in_per_level`` / ``_c_out_per_level``).

   Step B1 (this dispatch) folds the ONE free seam — the
   :class:`~orpheus.sn.sweep.cache.GeometryCoefficients` populator
   (:meth:`~orpheus.sn.sweep.cache.GeometryCoefficients.from_mesh_and_quad`),
   which already holds ``sn_mesh`` and so reads
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_out_per_ordinate`
   /
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_in_per_ordinate`
   with zero plumbing.  The accessor pair is PUBLIC and polymorphic on the
   base
   :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`:
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   returns its precomputed per-level :math:`c` gathered to the
   :math:`(N,)` global-ordinate order; the Cartesian
   :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`
   returns the NEUTRAL zeros (:math:`\alpha=0,\ \tau=1 \Rightarrow c=0`).
   The dispatch is by closure TYPE, not by a ``coord ==`` branch in the
   cache.

   Step B1 is BIT-IDENTICAL: the closure computes :math:`c` from
   closure-:math:`\tau` (0-ULP equal to geometry-:math:`\tau`, pinned by
   the Step-A Leg-1 gate) and the SAME :math:`\alpha` the populator read,
   so the closure's per-level :math:`c` equals the inline :math:`c`
   bit-for-bit; the per-level :math:`\to (N,)` gather is a pure
   permutation (no arithmetic).  The anchor gate
   ``tests/sn/sweep/core/test_cache.py::test_cache_populator_matches_cell_balance_terms``
   pins the cache ``denom`` (which carries :math:`(\Delta A/w)\,c_{\rm out}`)
   to :func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms` at
   ``rtol=1e-14``, and the curvilinear regression snapshots stay unmoved.

   The remaining THREE inline ``c`` rebuild sites (they need CellVisit
   threading) are later dispatches (Step B2 / B3 / C).  See
   :mod:`orpheus.sn.sweep.pole_angular_closure` for the canonical
   accessor and :mod:`orpheus.sn.sweep.cache` for the folded
   consumer.

.. _sn-closure-c-on-cellvisit:

c_in / c_out reach the stateless DD scheme as CellVisit data — Step B2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. todo:: Archivist expansion needed.

   Step B2 folds the SECOND of the four inline ``c_in`` / ``c_out``
   rebuild sites — the matvec-twin residual
   :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual`
   (formerly rebuilding :math:`c_{\rm out} = \alpha_{\rm out}/\tau`,
   :math:`c_{\rm in} = (1-\tau)/\tau\,\alpha_{\rm out} + \alpha_{\rm in}`
   inline from the geometry-owned :class:`StreamingTerms`).

   The architectural crux: :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
   is deliberately STATELESS — it reads only the
   :class:`~orpheus.transport.spatial.scheme.CellVisit` packet + the
   :class:`~orpheus.transport.spatial.scheme.UpstreamState`, never the mesh or
   the angular closure.  So the closure-owned :math:`c` cannot reach
   ``DD.residual`` by coupling DD to the closure object (that would break
   the spatial :math:`\otimes` angular separation — the SPATIAL scheme
   must not see the ANGULAR closure's type).  Instead the constants travel
   as DATA: the :class:`~orpheus.transport.spatial.scheme.CellVisit` gains two
   angular-closure-owned fields
   (:attr:`~orpheus.transport.spatial.scheme.CellVisit.c_in` /
   :attr:`~orpheus.transport.spatial.scheme.CellVisit.c_out`, distinct in
   provenance from the geometry-owned
   :attr:`~orpheus.transport.spatial.scheme.CellVisit.streaming_terms`), and the
   single production site
   :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit` — through which
   ALL four ``dag_walk`` yield paths funnel (Pattern 2, no per-site
   divergence) — stamps them from
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_in_per_ordinate`
   /
   :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_out_per_ordinate`
   indexed by the GLOBAL ordinate (``direction_idx`` for slab / sphere,
   ``level_indices[p][m]`` for cylinder — mirroring
   :meth:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.streaming_terms`).
   ``DD.residual`` then reads ``visit.c_in`` / ``visit.c_out``; the
   :math:`(\Delta A/w)`-scaled assembly that follows is byte-unchanged —
   only the SOURCE of :math:`c` moved.

   Step B2 also completes the matvec's typed-consumer binding (Issue
   #226): the unified SN matvec reads ``sn_mesh.pole_angular_closure``
   typed against the
   :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
   ABC and drives the angular path through
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.precompute_psi_state`,
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.cell_contribution`,
   and
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.angular_adjoint`.
   These were declared ``@abstractmethod`` on the ABC (matching the
   precedent where
   :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` declares
   ``update`` / ``residual`` abstract so ``mesh.scheme`` consumers see the
   full contract) — making the ABC the COMPLETE strategy contract instead
   of declaring only ``__call__``.

   Step B2 is BIT-IDENTICAL: ``visit.c_in`` / ``visit.c_out``
   (closure-sourced) equal the former inline values 0-ULP — the closure
   computes :math:`c` from closure-:math:`\tau` (0-ULP equal to
   geometry-:math:`\tau`, Step-A Leg-1) and the SAME geometry-:math:`\alpha`,
   and the per-level :math:`\to (N,)` gather is a pure permutation.  The
   matvec residual path (fed by ``DD.residual``) stays bit-for-bit on the
   ``tests/sn/sweep/curvilinear/test_unified_matvec_{sphere,cylinder}.py``
   twin and on the DriftWarning-escalating
   ``tests/sn/sweep/core`` + ``tests/sn/solve`` snapshots.  The remaining
   TWO ``c`` rebuild sites
   (:func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms` for the
   ``DD.update`` solve path; the geometry-side :math:`\tau` producer) are
   Step B3 / C.  See :mod:`orpheus.transport.spatial.scheme` for the CellVisit
   fields and :mod:`orpheus.sn.mesh.augmented_mesh` for the production stamp.

   B2 review fixes (finishing pass).  THREE follow-ups landed after the
   carve, all bit-identical (0-ULP):

   * **Per-ordinate gather cached (L16).** The public accessors
     :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_in_per_ordinate`
     /
     :attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.c_out_per_ordinate`
     re-ran the full :math:`(N,)` per-level :math:`\to` global gather on
     EVERY access, so the per-visit stamp made the visit-producing loop
     :math:`O(N^2\,n_x)`.  The gather is a pure permutation of immutable
     per-level data, so it is now computed ONCE in each mesh-bound
     ``__init__`` (shared
     :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase._build_per_ordinate_cache`,
     called by both ``MorelMontryAngularSweep`` and
     ``IdentityAngularClosure``) and the accessors return the read-only
     cache (``setflags(write=False)`` guards the shared :math:`(N,)` view
     B1's ``GeometryCoefficients`` populator holds).  Measured on a
     ``sphere N=32 nx=200`` walk: :math:`\sim 32\,\text{ms} \to \sim
     22\,\text{ms}` per sweep (:math:`\sim 1.46\times`), value-identical.

   * **Committed production-stamp catcher (vv L11 Mode 11).** The original
     carve had NO committed test exercising ``_make_cell_visit``'s c-stamp
     — the matvec twin reads the closure's ``cell_contribution`` directly
     (never ``DD.residual``), and the diamond / cell-balance fixtures stamp
     visits with a SURROGATE.  A wrong global-ordinate map (a
     ``c_in``:math:`\leftrightarrow`\ ``c_out`` swap, a mis-scattered
     cylinder level block) would ship silently.
     ``tests/sn/sweep/core/test_cell_visit_c_stamp.py`` walks a REAL
     production ``dag_walk`` (sphere + multi-level cylinder + slab) and
     asserts every ``visit.c_in`` / ``visit.c_out`` equals the constants
     recomputed INLINE from that visit's OWN
     :class:`~orpheus.geometry.reduced_operator.StreamingTerms` at 0-ULP
     (the hand-transcribed independent reference, not the closure's own
     ``c`` — vv L11).  Mutation-verified: the ``c_in``\ /\ ``c_out`` swap
     reddens the sphere + cylinder cases.

   * **Test-surrogate dedup (Pattern 2).** The byte-identical
     ``_c_from_streaming_terms`` (``test_diamond.py``) and ``_visit_c``
     (``test_cell_balance_for_streaming.py``) hand-recomputes were unified
     into one shared ``tests/sn/sweep/core/_c_surrogate.py`` consumed by
     both files and the new catcher.

.. _sn-tau-c-on-cellvisit-live:

The live sweep + scan consume closure-owned τ / c — Step B3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Step B1 folded the one free seam (the cache populator); Step B2 carried
the redistribution constants :math:`c_{\rm in}` / :math:`c_{\rm out}`
onto the :class:`~orpheus.transport.spatial.scheme.CellVisit` so the
*apply-direction* residual could read them as data.  Step B3 is the
step that makes the **live** paths consume the closure-owned weight:
the per-cell sweep solve, the matvec solve, and the CumprodScan
fast path now all read the angular weight :math:`\tau` :eq:`mm-weights`
(and the derived constants :eq:`dd-mm-closure-constants`) off the
closure rather than off the geometry-owned ``StreamingTerms.tau_mm``.
After B3 there is **no live reader** of ``StreamingTerms.tau_mm``
anywhere in the sweep, scan, or matvec — precisely the precondition
that let Step C delete the geometry-side :math:`\tau` producer (the two
parallel producers could be reduced to one only once nothing live
depended on the soon-to-be-deleted one; that retirement has now landed
— see the close-out at the end of this section).

This is the third of the four c-fold sites (after the cache populator
in B1 and the residual twin in B2) and the fifth :math:`\tau` consumer.
Like its predecessors it is **bit-identical (0-ULP)**: the carve moves
the *source* of an already-correct number, it does not change the
number.  The sections below derive why :math:`\tau` belongs to the
angular closure, why the constants must travel as visit *data* rather
than as a closure *reference*, why the field default is :math:`1.0`
(and not the more obvious :math:`0.0`), how the three live consumers
share one operator, and what makes the fold provably bit-identical and
therefore a safe regression floor for Step C.

τ is an angular-scheme property the closure owns
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Morel--Montry weight

.. math::

   \tau_n = \frac{\mu_n - \mu_{n-\frac12}}{\mu_{n+\frac12} - \mu_{n-\frac12}}

(:eq:`mm-weights`, Bailey--Morel--Chang 2010 Eq. 43) is built **entirely**
from the angular quadrature: the ordinate cosine :math:`\mu_n`, the
neighbouring angular-cell edges :math:`\mu_{n\pm 1/2}`, and — for the
cylinder — the :math:`\mu`-level partition that groups ordinates.  Not
one of those inputs is a property of the *spatial* streaming geometry:
:math:`\tau` does not depend on the cell volume :math:`V_i`, the face
areas :math:`A_{i\pm 1/2}`, the surface-curvature redistribution area
:math:`\Delta A_i`, or the radial mesh at all.  It is a number attached
to an **ordinate**, not to a **cell**.

That :math:`\tau` had historically lived on
:class:`~orpheus.geometry.reduced_operator.StreamingTerms` (as
``tau_mm``) was an accident of where the curvilinear sweep was first
assembled — the streaming-geometry factory happened to be the object
in scope when the weighted-diamond closure was wired in, so it baked
the angular weight in alongside the genuinely geometric
:math:`\alpha`-dome and face areas.  The architectural correction
(Issue #236 Phase 2) is to give the weight back to its owner: the
**pole-angular closure**
(:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`),
which already binds the quadrature and already computes :math:`\tau`
from it in :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`.
The closure exposes it through one public, polymorphic accessor,
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`,
a read-only :math:`(N,)` array indexed by the global ordinate.  The
two concrete strategies answer it differently *by type*, with no
``coord ==`` branch anywhere:

* :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
  returns its own per-level :math:`\tau` gathered to the global order;
* :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`
  (Cartesian) returns the neutral :math:`\tau \equiv 1` — there is no
  angular redistribution in slab geometry, so the M-M weight reduces to
  its identity element (see below).

This is the same both-sites mint B1 applied to the :math:`c`-accessors:
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
is declared on the
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
ABC, so the contract is complete for every consumer.  (Phase B
originally declared these accessors on **both** the
``@runtime_checkable`` ``PoleAngularClosure`` Protocol and the ABC, to
serve structural-typing and nominal-inheritance consumers alike;
Issue #236 Phase 2 B2 retyped every consumer onto the ABC and Issue
#248 deleted the now-orphaned Protocol, so the ABC is the single
declaration site.)  The gather itself is a pure permutation
of the immutable per-level data, hoisted once into each mesh-bound
``__init__`` via the shared
:meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase._build_per_ordinate_cache`
(renamed from ``_build_c_per_ordinate_cache`` now that it gathers three
constants — :math:`c_{\rm in}`, :math:`c_{\rm out}`, and :math:`\tau`
— rather than two); the accessor returns the cached read-only view, so
the per-visit lookup is :math:`O(1)`.

The spatial ⊗ angular separation forbids coupling DD to the closure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If :math:`\tau` is owned by the angular closure but consumed in the
spatial cell update, the obvious move would be to hand the closure
object to the cell-update scheme so it can ask for
:math:`\tau`.  That move is **forbidden by design**, and the reason is
the load-bearing architectural fact of the SN sweep.

:class:`~orpheus.transport.spatial.diamond.DiamondDifference` is a **stateless
spatial discretization scheme**.  It reads only the per-cell
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet and the
sweep-resolved :class:`~orpheus.transport.spatial.scheme.UpstreamState`; it
never sees the mesh, the quadrature, or the angular closure.  The whole
point of the spatial :math:`\otimes` angular product is that the
spatial scheme is interchangeable (diamond difference, linear
discontinuous, ...) without knowing *which* angular treatment sits on
the other axis of the tensor product, and the angular closure is
interchangeable (Morel--Montry, identity, a future Carlson variant)
without knowing the spatial scheme.  Coupling
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` to
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
would collapse that product into a Cartesian-vs-curvilinear conditional
inside the spatial scheme — exactly the geometry dispatch the unified
body was built to delete.

So the constants travel as **data**, not as a **dependency**.  The
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet — which the
orchestrator already populates per cell and per ordinate — carries the
angular-closure-owned numbers as plain ``float`` fields:
:attr:`~orpheus.transport.spatial.scheme.CellVisit.c_in` and
:attr:`~orpheus.transport.spatial.scheme.CellVisit.c_out` (added in B2) and now
:attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` (B3).  They are stamped
at exactly **one** production site,
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit`, through which all
four ``dag_walk`` yield paths (slab, sphere, cylinder, cylindrical
pure-azimuthal degenerate) funnel — Pattern 2, no per-site divergence.
That site reads the closure's per-global-ordinate accessors and stamps:

.. code-block:: python

   closure = self.pole_angular_closure
   return CellVisit(
       cell_idx=cell_idx,
       streaming_terms=st,
       face_area_downstream=face_area_downstream,
       c_in=float(closure.c_in_per_ordinate[global_ordinate]),
       c_out=float(closure.c_out_per_ordinate[global_ordinate]),
       tau=float(closure.tau_per_ordinate[global_ordinate]),
   )

where ``global_ordinate`` is the global ordinate index resolved the
same way :meth:`~orpheus.geometry.reduced_operator.ReducedStreamingOperator.streaming_terms`
resolves it (``direction_idx`` for slab / sphere,
``level_indices[mu_level_idx][m]`` for cylinder).  The spatial scheme
downstream sees only ``visit.tau`` / ``visit.c_in`` / ``visit.c_out``;
it has no idea a closure produced them.  The provenance is recorded in
the field docstrings (the constants are distinct in origin from the
geometry-owned :attr:`~orpheus.transport.spatial.scheme.CellVisit.streaming_terms`),
but the *type system* never lets the spatial axis reach across to the
angular axis.

Why the CellVisit.tau default is 1.0, not 0.0
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` defaults to
:math:`1.0`, not the zero a numeric field usually defaults to.  The
default is the value the slab / Cartesian path leaves in place — the
identity closure supplies it through ``tau_per_ordinate`` — and the
choice is forced by the angular recurrence the value feeds.

:math:`\tau = 1` is the **neutral element** of the Morel--Montry weight.
With :math:`\tau_n = 1` the WDD angular closure
:eq:`wdd-closure` becomes
:math:`\psi_{n,i} = \tau_n\,\psi_{n+\frac12} + (1-\tau_n)\,\psi_{n-\frac12}
= \psi_{n+\frac12}`, i.e. the step (fully-outgoing) closure, and the
outgoing-face recurrence

.. math::

   \psi^a_{\rm out} = \frac{\bar\psi - (1-\tau)\,\psi^a_{\rm in}}{\tau}
   \;\xrightarrow{\;\tau = 1\;}\;
   \frac{\bar\psi - 0}{1} = \bar\psi

reduces to the **identity** in :math:`\bar\psi` — exactly what slab
geometry needs, where there is no angular redistribution and the
"angular-out" state is just the cell average.  Likewise the
denominator constant :math:`c_{\rm out} = \alpha_{n+1/2}/\tau` and the
scan split :math:`1/\tau`, :math:`(1-\tau)/\tau` are all well-defined
and reduce to their slab values (:math:`0`, :math:`1`, :math:`0`
respectively, since :math:`\alpha = 0` on the slab).

A :math:`0.0` default, by contrast, is a **divide-by-zero landmine**.
Every consumer of :math:`\tau` divides by it:
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update`'s angular
thread divides :math:`(\bar\psi - (1-\tau)\psi^a_{\rm in})` by
:math:`\tau`; the cache derives ``tau_inv = 1/tau``.  A visit that
reached the curvilinear branch with an un-stamped :math:`\tau = 0`
would produce a silent ``inf``/``nan`` rather than a loud error.  The
:math:`1.0` default makes the *un-stamped* visit compute the
**identity** transformation — the safe no-op — which is both correct
for the slab (where the stamp is genuinely neutral) and a benign
fallback that surfaces a missing stamp as a *wrong-but-finite* answer a
regression snapshot catches, rather than a ``nan`` that propagates.

The three live consumers and the L21 framing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three live paths read the closure-owned :math:`\tau` (or the constants
derived from it) after B3.  The first two are the **solve** and
**apply** directions of the same per-cell linear system; the third is
the vectorized scan form of the same recurrence — the apply / sweep
duality this page calls the **L21 twin-path** (two applications of the
*same* operator; cf. :ref:`phase-c-sweep-frame-matvec` for the
apply-direction matvec that is the twin of the curvilinear sweep).

**(1) The scalar solve helper.**
:func:`~orpheus.transport.spatial.cell_balance.cell_balance_terms` — the
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` solve
direction — no longer rebuilds :math:`c_{\rm in}` / :math:`c_{\rm out}`
from ``st.alpha_* / st.tau_mm``.  Its signature now takes them as
keyword inputs, and
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` supplies
them straight off the visit::

   terms = cell_balance_terms(
       visit.streaming_terms, visit.face_area_downstream, total_xs,
       upstream_state, c_in=visit.c_in, c_out=visit.c_out,
   )

The helper now reads :math:`\tau` **not at all** — it consumes only the
already-derived :math:`c` constants, which is the right factoring: the
cell-balance denominator :eq:`dd-solve` needs
:math:`(\Delta A_i / w_n)\,c_{\rm out}`, and the upstream numerator
needs :math:`(\Delta A_i / w_n)\,c_{\rm in}\,\psi_{n-\frac12}`, neither
of which references :math:`\tau` once :math:`c` is in hand.

**(2) The angular recurrence.** The other half of
:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.update` — the
Morel--Montry outgoing-angular-face thread — *does* need the raw
:math:`\tau`:

.. math::
   :label: dd-mm-angular-recurrence

   \psi^a_{\rm out}
   = \frac{\bar\psi - (1 - \tau)\,\psi^a_{\rm in}}{\tau}

and reads it from :attr:`~orpheus.transport.spatial.scheme.CellVisit.tau`
(stamped by
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh._make_cell_visit` from
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`)
rather than from ``visit.streaming_terms.tau_mm``.  This is the line
the :math:`1.0` default protects.

**(3) The CumprodScan split.** The vectorized fast path — the cumulative
product that replaces the per-cell Python loop for the curvilinear
sweep — needs the same recurrence in a form amenable to a forward scan.
:class:`~orpheus.sn.sweep.cache.GeometryCoefficients` sources
:math:`\tau` from
:attr:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.tau_per_ordinate`
and precomputes the split

.. math::
   :label: dd-mm-scan-split

   \texttt{tau\_inv} = \frac{1}{\tau},
   \qquad
   \texttt{mm\_a\_in\_coeff} = \frac{1 - \tau}{\tau},

consumed at the loss-representation scan recurrence (in
:mod:`orpheus.sn.loss_representation`) as

.. math::

   \psi^a_{\rm out}
   = \texttt{tau\_inv}\cdot\bar\psi
     - \texttt{mm\_a\_in\_coeff}\cdot\psi^a_{\rm in}
   = \frac{\bar\psi}{\tau} - \frac{1-\tau}{\tau}\,\psi^a_{\rm in},

which is algebraically identical to :eq:`dd-mm-angular-recurrence` — the
same operator, applied in the vectorized frame.  Precomputing
:math:`1/\tau` and :math:`(1-\tau)/\tau` is a legitimate
perform-once-at-construction hoist (L16): the closure exposes only the
**primitive** :math:`\tau` (Pattern 5 — build the primitive, not the
product), and each consumer derives the trivial :math:`1/\tau`,
:math:`(1-\tau)/\tau`, :math:`\alpha_{\rm out}/\tau` algebra it needs at
its own definition site.  The scan derivation lives in the cache; the
recurrence consumes it.  That keeps the closure's surface minimal (one
weight) while letting two structurally-different consumers (a scalar
Python recurrence and a vectorized NumPy scan) each shape it to their
reduction tree.

Why the fold is bit-identical, and the regression floor for Step C
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Step B3 is **bit-identical (0-ULP)**.  The argument has two legs, both
already established by the earlier carve steps:

#. **Closure-:math:`\tau` is 0-ULP equal to geometry-:math:`\tau`.** The
   closure's
   :func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
   is a line-for-line replica of the geometry factory's :math:`\tau`
   arithmetic (Step A), pinned by the **Leg-1 producer-equivalence
   gate** ``tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py``
   to the geometry-factory value (0-ULP) *and* to the
   structurally-independent reference
   :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
   (a different code path to the same BMC-2010-Eq. 43 weight).  Reading
   :math:`\tau` from the closure therefore yields exactly the same
   ``float64`` bits the former ``st.tau_mm`` read carried.

#. **The per-level :math:`\to (N,)` gather is a pure permutation.** No
   arithmetic happens between the closure's per-level :math:`\tau` and
   the global-ordinate :math:`(N,)` view the accessor returns — only a
   reindex.  So every derived quantity (:math:`c_{\rm in}`,
   :math:`c_{\rm out}`, ``tau_inv``, ``mm_a_in_coeff``) is bit-for-bit
   what the inline rebuilds produced.

This was confirmed not only by the regression snapshots but **in
process**: at sphere (single :math:`\mu`-level) and cylinder
(multi-level) configurations,
:math:`\lvert\texttt{visit.tau} - \texttt{st.tau\_mm}\rvert`,
:math:`\lvert\texttt{closure.tau\_per\_ordinate} - \texttt{st.tau\_mm}\rvert`,
and ``np.array_equal`` on the cache's ``tau_inv`` / ``mm_a_in_coeff``
against the geometry-:math:`\tau`-derived split were **all exactly
zero**.  The DriftWarning-escalating ``tests/sn/sweep/core``,
``tests/sn/sweep``, and ``tests/sn/solve`` snapshots stayed unmoved
(588 + 60 green), with zero drift escalation.

The bit-identity guarantee is what makes B3 the **regression floor for
Step C**.  Two committed catchers pin the new live paths so that the
retirement of the parallel geometry-:math:`\tau` producer (now landed,
Step C) cannot silently break them:

* The **Leg-1 producer-equivalence gate** pins
  ``closure.tau_per_ordinate`` to the BMC-2010-Eq. 43 reference, so the
  closure remains the correct sole producer after the geometry one is
  deleted.

* The **production-stamp catcher**
  ``tests/sn/sweep/core/test_cell_visit_c_stamp.py`` walks a real
  production ``dag_walk`` (sphere, multi-level cylinder, slab) and — in
  its dedicated :math:`\tau` arm — asserts every ``visit.tau`` equals
  the **independently recomputed** Morel--Montry weight for that visit's
  ordinate at 0-ULP.  The reference is
  :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
  (a structurally-independent code path to the same BMC-2010-Eq. 43
  weight, not the closure comparing against itself — vv L11): the test
  pins the stamp's *ordinate map*, the complement of what Leg-1 pins
  (the producers' *values*).  Before Step C this catcher used the
  *geometry-produced* ``st.tau_mm`` as its independent reference; when
  Step C deleted that field the oracle was re-pointed onto
  ``morel_montry_weights`` (with the cylinder clamp replicated), keeping
  it geometry-:math:`\tau`-free and still structurally independent of
  the closure under test.  This arm was added specifically because B3
  made the
  :attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` stamp **live** while
  the existing named twins never call the rewired reader (vv L11
  Mode 11): a mutation stamping ``tau = ... * 1.1`` drifts the converged
  cylinder scalar flux by :math:`\sim 0.2\,\%` with **no** other test
  red, so the dedicated arm is the only committed catcher of a
  :math:`\tau`-stamp ordinate-map error.  Mutation-verified RED on the
  :math:`\times 1.1` stamp across sphere + cylinder + slab; GREEN clean.

* The **seam-6 scan catcher**
  ``tests/sn/sweep/core/test_affine_carve_baseline.py`` reddens on a
  corruption of the CumprodScan :math:`\tau` split, pinning the third
  live consumer.

With these in place, **Step C has now deleted** the geometry-side
:math:`\tau` producer.  The :math:`\tau` blocks inside
:func:`~orpheus.geometry.reduced_operator.spherical_streaming`
(the ``mu_edge`` weight-sum loop)
:math:`\,/\,`
:func:`~orpheus.geometry.reduced_operator.cylindrical_streaming` (the
per-level ``eta_edge`` loop) and the slab synthetic were excised, and
the now-orphaned ``StreamingTerms.tau_mm``, ``StreamingTerms.alpha_in``
/ ``alpha_out`` (whose sole readers were the c-rebuild sites B1--B3 just
retired), ``ReducedStreamingOperator.tau_mm``, and
``ReducedStreamingOperator.tau_mm_per_level`` dataclass fields were
dropped — confident that nothing live depended on them.
See :mod:`orpheus.sn.sweep.pole_angular_closure` for the
``tau_per_ordinate`` accessor and the three-constant cache,
:mod:`orpheus.transport.spatial.scheme` for the
:attr:`~orpheus.transport.spatial.scheme.CellVisit.tau` field,
:mod:`orpheus.sn.mesh.augmented_mesh` for the single production stamp, and
:mod:`orpheus.sn.sweep.cache` for the scan split.

.. _sn-tau-step-c-closeout:

Step C close-out — the geometry-side τ producer is retired
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The retirement is a **surgical field excision plus a test-oracle
migration**, not a blind delete: the τ producer was already dead in
*production* (B3 left zero live readers), but it remained load-bearing
as a *test oracle* — the regression-floor catchers used the
geometry-produced ``st.tau_mm`` as their structurally-independent
reference.  Migrate-then-delete preserved the floor:

#. **The oracles were re-pointed first.** The production-stamp catcher
   and the surviving producer-equivalence legs now read
   :func:`~orpheus.derivations.discrete.sn.contamination.morel_montry_weights`
   — a different code path to the same BMC-2010-Eq. 43 weight, unclamped
   on both geometries, with the cylinder clamp :math:`\mathrm{clip}(\cdot,
   \tfrac12, 1)` replicated in the test surrogate.  This kept the floor
   green *while geometry-:math:`\tau` was still present*, proving the
   migration faithful before any deletion.  The two
   ``*_equals_geometry_factory_0ulp`` legs of
   ``test_tau_producer_equivalence.py`` (which compared closure-:math:`\tau`
   against the soon-to-be-deleted factory output) became vacuous and were
   retired; the independent-reference and clamp-difference legs survive.

#. **The producers were excised surgically.** The τ blocks are
   *interleaved* with still-live outputs:
   :func:`~orpheus.geometry.reduced_operator.spherical_streaming` shares
   its ``mu_edge`` array with the live ``mu_start`` (the Hébert §3.9.4
   starting direction :math:`\mu_{1/2} = -1.0`), and
   :func:`~orpheus.geometry.reduced_operator.cylindrical_streaming`
   shares its per-level loop with the live ``mu_start_per_level``.  A
   whole-function deletion would have been wrong; only the τ statements
   were removed.  The :math:`\alpha`-dome (``alpha_half`` /
   ``alpha_per_level``), the redistribution factor ``redist_dAw``, the
   face areas, and the starting-direction edges all **stay on the
   geometry operator** — they are genuinely geometric.

#. **The deletion was proven inert.** The bit-identity regression gates
   (run under an escalated ``DriftWarning``) showed **zero** failures
   across the sweep / scan / matvec suites, and the test-count delta
   reconciled exactly to the four retired ``*_equals_geometry_factory``
   legs and the two retired ``test_reduced_operator.py``
   τ-bit-identical tests — no silent test loss.  After deletion the
   re-pointed catcher was **mutation-verified RED** (a :math:`\times 1.1`
   stamp and a ``c_in`` :math:`\leftrightarrow` ``c_out`` swap both
   reddened it against the independent oracle), confirming the migrated
   catcher is a real catcher reading the independent reference, not a
   tautology against the closure.

.. note::

   The legacy ``__call__``-argument ``tau_mm`` on the unbound
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   path (``MorelMontryAngularSweep(sn_mesh=None)``, where :math:`\tau` was
   passed as a runtime argument because the closure is not mesh-bound) was
   a **separate surface** that **survived Step C** unchanged — it was the
   closure's own runtime parameter, not the geometry-side field the carve
   retired.  It was subsequently retired under
   `Issue #248 <https://github.com/deOliveira-R/ORPHEUS/issues/248>`_
   (landed in this same re-staging, which deleted the strategy
   ``__call__`` bundle entirely; see the contract-evolution note at
   :ref:`sn-pole-angular-closure-protocol`).

.. _sn-space-angle-separability-section:

Space ⊗ angle separability — the (spatial ⊗ angular) product capstone (Issue #236 Phase 3)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section closes the Issue #236 *(spatial* :math:`\otimes` *angular)
product* narrative on the theory page.  The campaign had three phases:

* **Phase 1 — pairing validity.**  The spatial closure (the diamond /
  weighted-diamond cell update of :eq:`dd-curvilinear-scalar`) and the
  angular closure (the Morel--Montry weight :math:`\tau` of
  :eq:`mm-weights`, the redistribution dome
  :eq:`alpha-recursion`) are two distinct, independently-selectable
  axes — a genuine tensor product, with separate injection points on
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (``cell_update=`` for the spatial
  scheme, ``pole_angular_closure=`` for the angular scheme).
* **Phase 2 — :math:`\tau`-ownership carve.**  The angular weight
  :math:`\tau` was moved off the geometry operator and onto the angular
  closure, so the angular axis literally *owns* its own discretisation
  knob (:ref:`sn-tau-closure-owned` through
  :ref:`sn-tau-step-c-closeout`).  That carve made the product
  *structural in the type system*, not merely conceptual.
* **Phase 3 — separability characterisation (this section).**  Having
  established the product and given each axis its own knob, the final
  question is: *how do the two error contributions combine?*  The answer
  is geometry-dependent, and it is the campaign's headline claim.

The decomposition is pinned permanently by the L1 MMS characterisation
gate :mod:`tests.sn.verification.mms.test_space_angle_separability`,
which carries ``@pytest.mark.verifies("sn-space-angle-separability")``
against :eq:`sn-space-angle-separability` below.

The space–angle error decomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Write the total SN discretisation error of the scalar flux as a function
of the two refinement parameters — the spatial mesh size
:math:`h \sim 1/n_{\rm cells}` and the angular (quadrature) order
:math:`N` (the ordinate count, or the azimuthal order :math:`n_\varphi`
in the cylinder).  Let :math:`E_{\rm space}(h)` be the error of the
spatial closure at infinite quadrature order and :math:`E_{\rm angle}(N)`
the error of the angular closure on an exactly-resolved spatial mesh.
The campaign's headline result is the geometry-split law

.. math::
   :label: sn-space-angle-separability

   E(h, N) \;\approx\;
   \begin{cases}
     E_{\rm space}(h) \,+\, E_{\rm angle}(N),
       & \text{Cartesian (slab / 2-D / 3-D): \textbf{separates}},\\[1.2ex]
     \max\!\bigl(E_{\rm space}(h),\, E_{\rm angle}(N)\bigr),
       & \text{curvilinear (sphere / cylinder): \textbf{gates}}.
   \end{cases}

The two regimes are distinguished operationally by the **mixed second
difference** — the discrete :math:`\partial^2 E / \partial h\,\partial N`
evaluated on a two-quadrature error table over a coarse/fine pair of
mesh sizes:

.. math::
   :label: sn-space-angle-cross-term

   M \;=\; E[h_1, N_1] - E[h_1, N_2] - E[h_2, N_1] + E[h_2, N_2],
   \qquad
   \frac{|M|}{\max(\Delta E_h,\, \Delta E_N)} \;
   \begin{cases}
     \ll 1, & \text{separable (additive),}\\
     = \mathcal{O}(1), & \text{gated (coupled).}
   \end{cases}

For an additively separable error, :math:`E[h,N] = f(h) + g(N)` exactly,
so the cross-term telescopes to zero: :math:`M = f(h_1) + g(N_1) -
f(h_1) - g(N_2) - f(h_2) - g(N_1) + f(h_2) + g(N_2) = 0`.  A non-zero
:math:`M` is therefore a *direct, mechanism-anchored* measurement that
the two axes interact — that the second mixed partial of the error
surface does not vanish.  This is the quantity the ST5 gate measures.

.. (vv-status rationale) Characterisation claim, now tested: both the
   law :eq:`sn-space-angle-separability` and its discriminator
   :eq:`sn-space-angle-cross-term` describe the STRUCTURE of the
   discretisation-error surface (the regime discrimination), not a
   solver eigenvalue or flux VALUE.  This is an L1 MMS-convergence-
   structure (math) claim per vv-principles — MMS does not reach the
   eigenvalue layer, so neither label is, or ever becomes, an
   eigenvalue / flux-value claim.  The ST5 characterisation gate
   ``test_space_angle_separability.py`` now carries the
   ``verifies`` markers for both labels, so each is ``documented`` AND
   ``tested``: the verifying ``tests`` edge is the characterisation
   gate (an L1 MMS gate), not a closed-form / semi-analytical value
   reference.  What each gate leg verifies:
     * :eq:`sn-space-angle-separability` (the geometry-split decomposition
       law) is pinned by all six legs as a *positive* signature — the
       Cartesian legs assert separability (mixed-second-difference
       :math:`\to 0`, N-independent spatial rate); the curvilinear legs
       assert gating (N-gated spatial rate).  The marker is FILE-level.
     * :eq:`sn-space-angle-cross-term` (the mixed-second-difference
       discriminator :math:`M`) is the gate's measured instrument: the
       three legs that assert directly on :math:`|M|/\max`
       (``test_cartesian_slab_iso_space_angle_separable``,
       ``test_cartesian_slab_p1_aniso_floor_n_independent``,
       ``test_sphere_cross_term_large_discriminates_from_cartesian``)
       carry the per-test ``verifies`` marker.  The discriminator is a
       quantity the gate *measures against a declared threshold*, not a
       passive derivation step, so the ``tested`` edge is a real
       coverage claim.
   The posture mirrors the pole-cell characterisation gate this gate is
   modelled on (#233).
.. vv-status: sn-space-angle-separability documented
.. vv-status: sn-space-angle-separability tested
.. vv-status: sn-space-angle-cross-term documented
.. vv-status: sn-space-angle-cross-term tested

Why the two axes factorize: LMM-1987 (spatial) × BMC-2010 (angular)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The decomposition is not an empirical accident — it is forced by the
structure of the asymptotic *diffusion-limit-consistency* literature,
which is **literally split into a spatial paper and an angular paper**.
This split is the strongest possible evidence that the consistency
conditions live on two separate axes.

* **The spatial condition** — Larsen, Morel & Miller, "Asymptotic
  solutions of numerical transport problems in optically thick,
  diffusive regimes," *Journal of Computational Physics*
  **69(2):283--324 (1987)**, DOI
  `10.1016/0021-9991(87)90170-7 <https://doi.org/10.1016/0021-9991(87)90170-7>`_
  (and Part II, Larsen & Morel, JCP **83(2):212--236 (1989)**, DOI
  10.1016/0021-9991(89)90229-5).  LMM analyse the *spatial* differencing
  scheme's diffusion limit (cells scaled so they are not optically
  thin): a scheme whose discrete spatial limit is itself a valid
  diffusion discretisation (linear-discontinuous, weighted-diamond with
  the right closure) is "substantially more accurate" than one without
  (bare diamond difference).  This is a condition **on the spatial axis
  alone** — the angular order does not enter.

* **The angular condition** — Bailey, Morel & Chang,
  [BaileyMorelChang2010]_, "The asymptotic diffusion-limit accuracy of
  S\ :sub:`N` angular differencing schemes," *Nuclear Science and
  Engineering* **165(2):149--169 (2010)**, DOI
  `10.13182/NSE08-66 <https://doi.org/10.13182/NSE08-66>`_.  BMC analyse
  the SN equations **discretised only in angle, with space kept
  continuous** (their analysis deliberately removes spatial differencing
  to isolate the angular error).  They prove that the angular axis
  carries its *own* diffusion-limit condition, independent of the
  spatial one.  Their p. 151 statement is the separability fact in the
  authors' own words: the spatial half "has been shown by Larsen, Morel,
  and Miller," while "retaining full first-order consistency can be
  important for **angular** discretisations" — the angular contribution
  they introduce.

The two conditions factorise.  In the leading-order (:math:`\varepsilon^0`)
diffusion limit, *any* weighted-diamond angular weight (step, diamond,
Morel--Montry) preserves consistency — BMC Eqs. (23)--(25).  The
**first-order** (:math:`\varepsilon^1`) limit carries a contamination
term :math:`\beta` (BMC Eq. 40), a **purely angular** functional of the
redistribution coefficients and quadrature,
:math:`\beta = \sum_m \mu_m\bigl[\alpha_{m+1/2}\mu_{m+1/2} -
\alpha_{m-1/2}\mu_{m-1/2}\bigr]`, which vanishes *only* for the
Morel--Montry weights (BMC Eq. 43, the weight of :eq:`mm-weights`).
Because :math:`\beta` depends on no spatial quantity, the angular
condition is provable on its own axis — exactly as the spatial
condition is provable on its own axis.  The diffusion limit needs
**both**:

.. math::

   \text{accurate diffusion limit}
   \;\;\Longleftrightarrow\;\;
   \underbrace{(\text{LMM spatial condition})}_{\text{depends on the spatial scheme only}}
   \;\;\wedge\;\;
   \underbrace{(\text{BMC angular condition},\ \beta = 0)}_{\text{depends on the angular weights only}}.

This conjunction of two single-axis conditions is *why* the Cartesian
error separates additively (each axis contributes its own consistency
defect, and the two defects add) and *why* a bad pairing can still break
the limit (independence of *selection* is not independence of
*consequence* — both conditions must hold simultaneously).

.. note::

   The literature's double-use of the name "linear-discontinuous" (LD)
   is itself evidence of the two-axis structure: LMM and every
   spatial-scheme paper list LD as a *spatial* scheme, while Lathrop
   (2000) lists "linear-discontinuous" among his *angular* differencing
   schemes.  The same trial-space name applies on either axis; the
   ORPHEUS registries disambiguate by axis (a spatial cell-update vs an
   angular closure), never a single ``LD`` enum.  This is the #158
   (spatial scheme) vs #6 (LD *angular* finite elements) distinction.

Cartesian separates, curvilinear gates
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The geometry split of :eq:`sn-space-angle-separability` follows
mechanically from *where the angular redistribution term lives*.

**Cartesian — additive separation.**  In slab / 2-D / 3-D Cartesian
geometry the curvilinear angular-redistribution term
:math:`\frac{1-\mu^2}{r}\,\partial_\mu\psi` is **absent**
(:math:`r \to \infty`; there is no :math:`1/r`).  The angular closure is
the :class:`~orpheus.sn.sweep.pole_angular_closure.IdentityAngularClosure`,
which contributes no redistribution: each ordinate's spatial sweep is
fully independent of every other ordinate.  The Cartesian cell update
(:eq:`dd-2d-balance-form` / the slab balance) consumes only the per-axis
streaming ratios and :math:`\Sigma_t` — **no** :math:`\tau`, **no**
angular state.  The spatial error and the angular (quadrature) error are
generated by disjoint mechanisms, so they add:
:math:`E(h,N) \approx E_{\rm space}(h) + E_{\rm angle}(N)`, and the
mixed partial vanishes.  The operational signature is that the spatial
convergence **rate** is the same at every quadrature order
(N-independent O(h\ :sup:`2`)).

**Curvilinear — multiplicative gating.**  In the sphere / cylinder the
Morel--Montry angular thread

.. math::

   \psi_{n+\frac12} \;=\;
       \frac{\overline{\psi}_n - (1-\tau_n)\,\psi_{n-\frac12}}{\tau_n}

couples the ordinates *sequentially within a* :math:`\mu`-*level*, and
the coupling enters the **shared cell-balance denominator** of
:eq:`dd-curvilinear-scalar`: the redistribution divisor
:math:`(\Delta A_i / w_n)\,c_{\rm out}` (with
:math:`c_{\rm out} = \alpha_{\rm out}/\tau_n`) sits in the *same*
denominator that produces the spatial cell average
:math:`\overline{\psi}_{n,i}` the spatial closure then uses.  The
angular interpolation error of the :math:`\tau`-thread therefore
**caps** the accuracy the spatial closure can deliver: at a coarse
quadrature, refining :math:`h` cannot drive the cell average below the
angular floor, because the angular term contaminates the denominator
the spatial refinement acts through.  Hence the error *gates*:
:math:`E(h,N) \approx \max(E_{\rm space}(h), E_{\rm angle}(N))`.  You
cannot harvest fine-:math:`h` accuracy at coarse :math:`N`; both axes
must advance together.  The mechanism is documented in detail at
:ref:`sn-tau-c-on-cellvisit-live` (why :math:`\tau` is an angular
property that nonetheless flows through the spatial denominator) and the
shared-denominator algebra is :eq:`dd-curvilinear-scalar`.

.. warning::

   The gating is a property of *today's* curvilinear closure (the 1-D
   :math:`\eta`-march Morel--Montry thread), not a law of nature.  A
   future 2-D angular closure (#229) that resolves the
   :math:`(\eta,\varphi)` azimuthal variation the 1-D march cannot
   thread, or a higher-order spatial scheme (#158 / #6), would *lift*
   the gating — at which point the curvilinear error would begin to
   separate.  The ST5 gate is designed so that lifting the floor reddens
   the gating assertions (the coarse-N saturated h-ratio rises toward
   the O(h\ :sup:`2`) value), signalling that the regime changed; that
   redding is the intended signal to *re-tune* the gate to the new,
   better regime, **not** a regression.

Measured cross-term evidence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The decomposition was established empirically by the four
``diag_sep_*`` probes (reproduced bit-for-bit after the Phase-2
:math:`\tau` carve) and is now pinned by the ST5 gate.  The measured
mixed-second-difference :math:`|M|/\max` spans **three orders of
magnitude** between the separable Cartesian and the gated sphere — a
clean discrimination band, not a brittle exact number.

.. list-table:: Measured space–angle error structure (2026-06-18,
   ``nx ∈ {20, 40, 80}``)
   :header-rows: 1
   :widths: 18 14 30 24 14

   * - Geometry
     - Regime
     - Scalar-flux L2 ladder (coarse → fine quadrature)
     - Spatial h-ratios (coarse-N / fine-N)
     - :math:`|M|/\max`
   * - Slab (isotropic)
     - **separates**
     - N=4 ``[5.42e-4, 1.35e-4, 3.38e-5]`` · N=16 ``[5.40e-4, 1.35e-4, 3.37e-5]``
     - ``[4.01, 4.00]`` / ``[4.01, 4.00]`` (N-independent O(h²))
     - **0.0047**
   * - Slab (P1 aniso)
     - **separates**
     - N=4 floor ``6.80e-3`` · N=16 floor ``6.79e-3`` (flat, angular floor)
     - flat at both N (floor N-independent to <0.3 %)
     - **0.0038**
   * - Cylinder
     - **gates**
     - :math:`n_\varphi`\ =8 ``[1.95e-2, 1.91e-2, 1.90e-2]`` · :math:`n_\varphi`\ =16 ``[8.05e-3, 7.47e-3, 7.37e-3]``
     - ``[1.02, 1.00]`` (saturated); azimuthal floor drops 2.58× at :math:`n_\varphi` 8→16
     - **0.019** (small only because :math:`E \approx E_{\rm angle}` swamps)
   * - Sphere
     - **gates**
     - N=8 ``[1.47e-2, 5.40e-3, 4.69e-3]`` · N=32 ``[1.50e-2, 3.71e-3, 9.29e-4]``
     - ``[2.71, 1.15]`` (saturates) / ``[4.04, 4.00]`` (O(h²) recovers)
     - **0.411**

The reading of the table:

* **Slab (both rows): separable.**  The spatial h-ratio is :math:`\approx
  4` (O(h\ :sup:`2`)) at *every* quadrature order — the spatial rate is
  blind to :math:`N`.  The isotropic row has a genuine O(h\ :sup:`2`)
  window; the P1-anisotropic row sits at a flat MMS/angular floor that
  is the *same* at every :math:`N`.  Both have :math:`|M|/\max \le
  0.005` — the cross-term vanishes whether or not the angular axis is
  active.  The P1 row is the load-bearing control: separability survives
  an *active* angular term, so it is not an artefact of the isotropic
  degeneracy.

* **Sphere: gating, the discriminator.**  At coarse N=8 the finest
  spatial h-ratio collapses to **1.15** — refinement saturates at the
  angular floor.  At fine N=32 the *same* spatial ladder recovers
  :math:`\approx 4.00` (O(h\ :sup:`2`)).  The spatial rate *depends on*
  :math:`N` — the defining gating fact — and the cross-term
  :math:`|M|/\max = 0.411` sits three orders above the Cartesian
  ceiling.

* **Cylinder: the extreme of gating.**  There is no pre-floor
  O(h\ :sup:`2`) window at any practical azimuthal order — the
  :math:`(\eta,\varphi)` variation a 1-D :math:`\eta`-march cannot
  thread exactly (#229) — so :math:`E \approx E_{\rm angle}(n_\varphi)`
  and the spatial h-ratio is :math:`\approx 1` at fixed :math:`n_\varphi`.
  The positive signature is the floor's azimuthal scaling: it drops
  :math:`2.58\times` when :math:`n_\varphi` doubles.  (The cylinder's
  small :math:`|M|/\max` is *not* evidence of separability — it is small
  because the angular floor so dominates that the spatial delta
  :math:`\Delta E_h` in the denominator is itself near zero; the gating
  is read from the *saturation* and the *azimuthal scaling*, not the
  cross-term magnitude.)

.. note::

   The scalar (weight-summed) L2 of the table is, by construction, blind
   to a *wrong angular closure* — the Morel--Montry :math:`\alpha`-dome
   telescopes under :math:`\sum_n w_n \psi_n` (vv-principles L27 / the
   per-ordinate-flat-flux discipline).  Because the curvilinear gating is
   *itself* an angular-closure phenomenon, the ST5 gate adds a
   **per-ordinate** leg
   (``test_curvilinear_gating_per_ordinate_not_blind``) that reproduces
   the sphere gating signature from the max-over-ordinates per-ordinate
   L2 (N=8 finest h-ratio 1.16 saturates; N=32 recovers ≈3) — so the
   gate cannot be telescoped blind to a future angular-closure
   regression.  That leg corrects a measured 1/W normalisation trap:
   ``case.psi_exact(r, μ_n)`` returns :math:`A(r) + B(r)\mu_n` *without*
   the :math:`1/W` factor by its own contract, while the solver stores
   the per-ordinate flux *with* it — the reference must be divided by
   :math:`W = \sum_n w_n` before comparison, else a 2× mismatch swamps
   the metric.

The #233 pole-cell × #229 azimuthal-floor interference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The gating law has a concrete consequence for the two open curvilinear
defects.  The spatial pole-cell :math:`\mathcal{O}(h)` order (#233,
documented at :ref:`sn-pole-cell-spatial-closure` and minted as ERR-059
in :ref:`sn-curvilinear-aniso-norm-reconciliation`) and the azimuthal
angular floor (#229) are **not independent contributors** to the
curvilinear error — they *interfere through the gating*.

The mechanism: the angular thread's interpolation error sets the floor
that spatial refinement saturates at.  So the pole-cell spatial defect
(#233) is only *visible* — only the dominant error — once the angular
floor (#229) has been pushed below it by refining :math:`N`.  At a
coarse quadrature the #229 angular floor *masks* the #233 pole-cell
order entirely (the spatial ladder saturates before the
:math:`\mathcal{O}(h)` pole-cell term emerges); only at a fine
quadrature does the spatial ladder run long enough for the pole-cell
order to surface.  This is precisely why the sphere N=8 ladder saturates
at 1.15 while N=32 recovers O(h\ :sup:`2`): the same spatial closure,
read through two different angular floors.

This interference is the reason the two issues must be characterised
*together* rather than as separate spatial and angular bugs, and the
reason a fix to one cannot be validated in isolation: lifting the #229
angular floor (a 2-D angular closure) would *expose* the #233 pole-cell
order that the floor currently masks.  The gating law
:eq:`sn-space-angle-separability` makes this dependency explicit — the
curvilinear error is :math:`\max(E_{\rm space}, E_{\rm angle})`, so
whichever defect is larger *hides* the other.

The permanent pin: the ST5 characterization gate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The decomposition is pinned permanently by
:mod:`tests.sn.verification.mms.test_space_angle_separability` (Issue
#236 Phase 3, sub-task ST5), an **L1** MMS characterisation gate modelled
on the pole-cell characterisation gate
``test_curvilinear_pole_cell_characterization.py`` (#233).  It carries
``@pytest.mark.verifies("sn-space-angle-separability")`` against
:eq:`sn-space-angle-separability`.  Its six legs pin both regimes as
*positive* signatures (never as an xfail-pending-fix):

* **Cartesian, separable** — ``test_cartesian_slab_iso_space_angle_separable``
  (N-independent O(h\ :sup:`2`) spatial rate, :math:`|M|/\max < 0.05`)
  and ``test_cartesian_slab_p1_aniso_floor_n_independent`` (the active-
  angular-axis control: the P1 floor is N-independent and the cross-term
  stays :math:`\approx 0`).
* **Curvilinear, gating** —
  ``test_sphere_spatial_rate_is_quadrature_gated`` (the discriminator:
  coarse-N saturates, fine-N recovers O(h\ :sup:`2`); the proven
  ``@catches("ERR-026")`` catcher via the fine-N O(h\ :sup:`2`)-recovery
  assertion), ``test_sphere_cross_term_large_discriminates_from_cartesian``
  (:math:`|M|/\max > 0.15`),
  ``test_cylinder_spatial_saturates_at_azimuthal_floor`` (spatial
  saturation + azimuthal floor scaling), and
  ``test_curvilinear_gating_per_ordinate_not_blind`` (the L27
  angular-aware per-ordinate leg, also ``@catches("ERR-026")``).

The gate is *characterisation*, not calcification: if a future 2-D
angular closure (#229) or higher-order spatial scheme (#158 / #6) lifts
the curvilinear gating, the gating assertions are designed to redden so
the regime change is *signalled* and the gate is re-tuned to the new
(better) regime — they are not xfails awaiting a fix.  The ``@slow``
mark reflects that the curvilinear solves dominate the ~2 s wall-clock,
not that the gate is optional.

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
:mod:`orpheus.transport.spatial.scheme`.

The Protocol
------------

The contract is a ``@runtime_checkable`` ``typing.Protocol`` —
satisfied by structural typing, not inheritance — exposing two class-
level traits and a single :meth:`update` method:

* :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme`

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
    :class:`~orpheus.transport.spatial.scheme.CellVisit` packet (see
    next subsection) that combines the geometric
    :class:`~orpheus.geometry.reduced_operator.StreamingTerms` with
    sweep-direction-resolved data.

The two helper dataclasses (frozen, slotted) carry the per-cell
state:

* :class:`~orpheus.transport.spatial.scheme.UpstreamState`

  - ``spatial_upstream: np.ndarray`` — shape ``(ng,)``.  Face flux
    entering the cell from the upstream face (always populated).
  - ``angular_upstream: np.ndarray | None`` — shape ``(ng,)`` for
    sphere/cylinder; ``None`` for slab.  :math:`\psi_{n-1/2,\,i}`,
    the half-flux at the upstream half-angle.

* :class:`~orpheus.transport.spatial.scheme.CellResult`

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
:class:`~orpheus.transport.spatial.scheme.CellVisit` packet rather
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
<orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk>` — a generator that
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

.. note:: **Superseded mechanism (Issue #196 Phase G Step 2.5 → Issue
   #236 Step C).** The ``alpha_in is None`` slab/curvilinear
   discrimination described below was **retired**: Issue #196 Phase G
   Step 2.5 gave slab *neutral* curvature (``face_area_inner =
   face_area_outer = 1.0``, ``delta_A_over_w = 0.0``) so the unified
   cell-balance helper consumes the same packet regardless of geometry,
   and Issue #236 Step C removed the Morel--Montry ``alpha_in`` /
   ``alpha_out`` / ``tau_mm`` fields from
   :class:`~orpheus.geometry.reduced_operator.StreamingTerms` entirely
   (it is now **purely geometric**; :math:`\tau` is closure-owned — see
   :ref:`sn-tau-c-on-cellvisit-live`).  Slab is now distinguished at the
   sweep level by ``upstream_state.angular_upstream is None``, not by a
   ``StreamingTerms`` field test.  The prose below records the historical
   pre-Step-2.5 convention and is pending a dedicated rewrite to the
   neutral-curvature mechanism; the authoritative current description is
   the :class:`~orpheus.geometry.reduced_operator.StreamingTerms`
   docstring.

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
:doc:`/theory/foundations/structured_geometry`, "Connection coefficients (reduced
streaming operator)").  In this case
:meth:`SNMesh.dag_walk
<orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk>` yields visits with
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
:class:`~orpheus.transport.spatial.scheme.DiscretizationScheme` strategy.

The first concrete strategy — :class:`DiamondDifference` — is shipped
in Round 2 of Wave C (Issue #158) as a bit-identical extraction of
the existing inlined sweep math.  Linear Discontinuous (Lewis &
Miller §5.3 — preview), Exponential Characteristic, and Step
strategies are deferred to a Wave C-extension session, each with its
own MMS spatial-convergence verification.

Diamond Difference
------------------

The first concrete strategy is
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`
(:mod:`orpheus.transport.spatial.diamond`).  It implements the **same**
algebra as the existing inlined sweep — Round 2 of Wave C is a
bit-identical extraction, gated by ``np.array_equal`` hand-calc tests
in ``tests/sn/sweep/core/test_diamond.py`` against the sweep's scalar
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
:attr:`~orpheus.transport.spatial.scheme.CellResult.outgoing_angular_state`
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
:attr:`~orpheus.transport.spatial.scheme.CellResult.outgoing_spatial_flux`
``= None`` to signal "no face-flux write" to the sweep driver; the
M-M angular closure remains active.

**Traits and forward references.**  Diamond Difference has

* :attr:`~orpheus.transport.spatial.diamond.DiamondDifference.is_linear`
  ``= True`` — the cell average and downstream states are affine
  combinations of ``source`` and ``upstream_state``;
* :attr:`~orpheus.transport.spatial.diamond.DiamondDifference.is_positivity_preserving`
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

* :mod:`orpheus.transport.spatial.scheme` — the contract module.
* :mod:`orpheus.transport.spatial.diamond` — the
  :class:`~orpheus.transport.spatial.diamond.DiamondDifference` concrete
  strategy.
* :doc:`/theory/foundations/structured_geometry`, "Connection coefficients (reduced
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
:class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase`, so the per-scheme
classes (:class:`~orpheus.transport.spatial.diamond.DiamondDifference`,
:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous`, Step)
carry NO inlined reconstruction of their own — they differ only by the value
they pass for the blend weight :math:`w`.  The #240 Phase 2 Step D1 carve
collapsed the previously-duplicated inline forms (Diamond Difference's
:math:`2\bar\psi - \psi_{\rm in}`, Linear Discontinuous's
:math:`(1+k)\bar\psi - k\,\psi_{\rm in}`) onto this one op
(:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.outgoing_face_from_average`),
the algebraic inverse of
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_average`; Step D2
made the trio (``source_emission`` / ``cell_average`` /
``outgoing_face_from_average``) generic advection–reaction reconstructions
(diffusion-consumable, retiring the dangling ``affine_closure`` module).  The
unit gate is :mod:`tests.transport.spatial.test_affine_closure`: the exact-inverse
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
  :mod:`orpheus.transport.spatial.linear_discontinuous` 2×2 entry-for-entry, with the
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

The foundation gate is :mod:`tests.transport.spatial.test_ld_ubld_symbolic` (6
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
is :mod:`orpheus.transport.spatial._ubld`, in two layers that share ONE source
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
:mod:`orpheus.transport.spatial.linear_discontinuous` — the ×V per-cell Schur
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
  :class:`~orpheus.transport.spatial._ubld.D1ClosedForm` is the analytic Schur
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

   * - Production view (in :mod:`orpheus.transport.spatial.linear_discontinuous`)
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

The verification is :mod:`tests.transport.spatial.test_ld_ubld_primitive` (10 tests):
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
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.spatial_basis_per_axis`
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
:mod:`orpheus.transport.spatial.linear_discontinuous` (``cell_kernel_batch`` /
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
face-cochain trailing axis (:mod:`orpheus.sn.loss_representation.sweep_graph` ``_MovingFrontier``;
the ``_CellSolve`` / ``_CellResidual`` moment-reducing emit), and the window
zero-pad (:mod:`orpheus.sn.loss_representation`
``FullFieldWavefront._octant_face_cochain``, the ``_inflow_to_moments`` pad) —
are GATED on the single scheme trait
:attr:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.spatial_basis_per_axis`
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

* The non-vanishing domain-inflow moment trace (the ``AngularBoundaryFlux`` /
  ``mesh.angular_trace`` widening to a :math:`2^{d-1}` transverse face moment) — **S4**
  (and its honest-scope caveat, :ref:`ld-cartesian-2d`).

* The strengthened vv Mode-7 stress-ansatz MMS and the thick-diffusive
  tripwire — **S4** and **S3** respectively.

The verification is the kernel round-trip + matvec-twin face reconstruction
(:mod:`tests.transport.spatial.test_linear_discontinuous` ``TestLDKernel``), the
end-to-end two-paths FFW :math:`\equiv` MFW, the DD :math:`\ne` LD routing-flip,
and the :math:`O(h^2)` convergence smoke
(:mod:`tests.sn.verification.mms.test_mms_ld_2d`), plus the :math:`d=2`
numpy↔symbolic entry-wise ``A == A`` cell-assembly pin and the
``test_d2_exact_on_bilinear`` ERR-060 catcher
(:mod:`tests.transport.spatial.test_ld_ubld_primitive`).

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
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
on the space
:math:`\mathrm{SphericalHarmonicSpace}(L) \otimes \mathrm{CellGroupSpace}`.
The angular moment is a **replacement representation** of the angular flux:
a calculation holds EITHER the per-ordinate field
:math:`\psi` of shape :math:`(N, ng, *\text{spatial})` OR its harmonic
moments :math:`\phi_\ell^m` of shape
:math:`(L{+}1, 2L{+}1, ng, *\text{spatial})`, bridged by the
spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`'s two faces —
its ``analysis`` face (:math:`M`, the
:math:`\psi \to \phi` reduction :eq:`two-moment-angular`) and its
``reconstruction`` face
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
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`:

.. list-table:: The :math:`\Sigma_s \otimes I_{\rm spatial}` lift in code
   :header-rows: 1
   :widths: 34 30 36

   * - Producer
     - Subscript (pre-S3 :math:`\to` S3)
     - What it scatters
   * - :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_p0_in_scatter`
     - ``"fg,fc->gc"`` :math:`\to` ``"fg,fc...->gc..."``
     - the P0 in-scatter :math:`\Sigma_{s,0}^{\mathsf T}\phi`
   * - :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_n2n`
     - ``"fg,fc->gc"`` :math:`\to` ``"fg,fc...->gc..."``
     - the :math:`(n,2n)` source :math:`2\Sigma_{2n}^{\mathsf T}\phi`
   * - :meth:`~orpheus.transport.mesh.material_xs_field.MaterialXSField.apply_legendre_scattering_moments`
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

* **The projection pair needed no change.**  The spherical-harmonic
  :class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` and
  ``reconstruction`` faces already carry ``...`` for their trailing
  axes, so :math:`M` and :math:`R` are spatial-moment-agnostic out of
  the box.  The angular reduction
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
  :meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.kernel_rhs` and the Schur
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
:meth:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux.from_mesh_and_L`)
gained an OPTIONAL ``spatial_moments`` parameter (default ``1``) that
appends the factor **iff the within-cell count exceeds 1** — the
"append iff > 1" gate single-sourced from
:func:`~orpheus.numerics.spaces.spatial_moment_space.spatial_moment_tail`
(the cell analogue that delegates to
:func:`orpheus.numerics.moment_layout.face_moment_tail`, so the cell-moment tail
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
   :class:`HarmonicMomentFlux` docstrings (issue #207) but had never been
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
:func:`orpheus.transport.spatial._ubld.assemble_ubld`).  The convention is fixed so
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
:data:`orpheus.numerics.moment_layout.AVERAGE_MOMENT` (the constant every moment
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
identical* to its pre-S3 self.  :func:`orpheus.numerics.moment_layout.face_moment_tail`
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
  :math:`==` :data:`~orpheus.numerics.moment_layout.AVERAGE_MOMENT`.

* ``tests/numerics/test_spatial_moment_field_space.py`` — the factory
  widening: the **byte-identity-at-default negative control** for DD AND LD on
  all three carriers (:class:`AngularField`, :class:`ScalarField`,
  :class:`HarmonicMomentFlux`), the widened :math:`d{=}1` / :math:`d{=}2`
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
The source files are :mod:`orpheus.transport.spatial.linear_discontinuous`
(``cell_kernel_batch`` / ``residual_kernel_batch`` — now ONE :math:`d`-generic
moment path), :mod:`orpheus.sn.loss_representation.sweep_graph` (``_CellSolve`` / ``_CellResidual``
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
:func:`~orpheus.transport.spatial._ubld.octant_moment_frame_signs`,

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
:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm._slope_fold`, shared by the
per-cell matvec Schur (:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.schur_xV`)
AND the scan
(:meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.scan_slope_face_source` for the
face-chain term, :meth:`~orpheus.transport.spatial._ubld.D1ClosedForm.scan_reconstruct`
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
:func:`~orpheus.transport.spatial._ubld.octant_moment_frame_signs` involution
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

The SymPy module / live scheme is :mod:`orpheus.transport.spatial._ubld`
(``D1ClosedForm``),
:meth:`orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous.moment_scan_closure`,
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
     (**#251**).  The boundary trace ``mesh.angular_trace`` carries the
     :math:`2^{d-1}` transverse face-moments per face per ordinate per group
     (a moment-resolved slot ``(N, ng, *face_shape, 2^{d-1})`` minted by
     :attr:`orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout`, appending the
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
   :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` moment-tail
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
:meth:`~orpheus.sn.loss_representation.sweep_graph._CellSolve.cell` via the slope-sign involution
:math:`\mathrm{octant\_moment\_frame\_signs}`
(Eq. :eq:`ld-ubld-octant-moment-frame-signs`,
:func:`orpheus.transport.spatial._ubld.octant_moment_frame_signs`), exactly as it
reframes the scattering slope source.  So the external :math:`\hat Q` rides the
SAME global→sweep machinery the scattering moments already use — the
:ref:`ld-ubld-sweep-global-frame` involution that the S3 unified matvec had to
get right (ERR-061) is reused unchanged, with no new cell branch.

The projector is structurally independent of production (L11): it evaluates
:math:`\int q\,P_k` with :func:`numpy.polynomial.legendre.leggauss` directly and
NEVER calls ``_lift_external_source_to_moments`` nor any
:class:`~orpheus.transport.spatial.linear_discontinuous.LinearDiscontinuous` cell op.
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
``mesh.angular_trace`` so a moment-resolved prescribed inflow can carry the along-face
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
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` itself but by the
:class:`~orpheus.numerics.face_layout.FaceLayout` it is built from, and that
layout is minted by :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout`.
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

* **The** :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` **needs ZERO
  changes — it was "moment-ready by accident".**  The trace's partial-current
  metric (the :math:`|\Omega\!\cdot\!n|\,w_n` weighting) and its
  ``omega_dot_n`` inflow/outflow table are **per-ordinate** (they classify and
  weight by ordinate, axis 0) and broadcast over ALL trailing axes by
  construction.  A moment axis appended to the slot rides the metric and the
  directional selectors for free.  The trace space was already
  moment-polymorphic; only its layout-supplier was scalar.  This is the
  illegal-states-unrepresentable payoff of the Depth-B field-space refactor: the
  boundary FIELDS (``AngularBoundaryFlux`` / ``AngularBoundarySourceSink`` /
  ``AngularBoundaryResidual`` / ``AngularBoundaryDisplacement``) validate ONLY against
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
:func:`orpheus.transport.spatial._ubld.assemble_inflow_axis`, which weights the
:math:`2^{d-1}`-moment face vector by the **transverse mass** — the Kronecker
product of the per-transverse-axis :func:`~orpheus.transport.spatial._ubld.mass_1d`
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
:meth:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink.prescribed_inflow`.
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
   :meth:`orpheus.sn.operators.boundary.SNBoundaryOperator._reflect_trace`, and
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
:attr:`orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` (appending
:func:`~orpheus.numerics.moment_layout.face_moment_tail`); the inflow lift
:meth:`_LossRepresentation._inflow_to_moments`, the oracle seed
:meth:`FullFieldWavefront._octant_face_cochain`, and the four outflow
capture-collapse DROP sites in :mod:`orpheus.sn.loss_representation`; and the
producer
:meth:`~orpheus.transport.source_sinks.angular_boundary_source_sink.AngularBoundarySourceSink.prescribed_inflow`.
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
path: the trace ``mesh.angular_trace`` carries the transverse moment axis end to end,
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
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout` keys the slot width on:

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
:doc:`algebra-of-record </theory/references/index>` pillar).  Per transverse
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
(:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`)
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
harmonic-modal basis (``HarmonicMomentFlux``) are non-canonically isomorphic
(the iso depends on the quadrature :math:`Y_\ell^m(\hat\Omega_n)`), bridged by
the applied :math:`M`/:math:`R` projection/reconstruction pair with truncation
content and adjoints.  The transverse SPATIAL moment FAILS: there is ONE
within-cell basis (the tensor-Legendre tower), no non-canonical dual, so the
only change-of-basis would be the identity.  A ``BoundaryMomentField`` leaf
whose ``_check_partner`` adds nothing beyond class identity would be a vacuous
naming leaf — type-theatrics by the project's own "if the type hint does not
prevent a bug by construction it is theatrics" standard.  So the moment rides as
a PROPERTY (the flat face buffer already holds the moment tail via
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.boundary_face_layout`), and the first-class
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
(:class:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel`,
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` and its two
storage walks, the level-operation pair ``_CellSolve`` /
``_CellResidual``, and the discretization's kernel pair
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
/ :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`)
are documented in detail at
:ref:`sweep-octant-dependency-graph` immediately below.

The §15A.2 sum-of-tensor-products framing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Following Grand Report v3 §15A.2 (lines 2137–2171), the loss inverse
(the sweep) :math:`A^{-1} = (L+C)^{-1}` on the 2-D Cartesian SN field
space decomposes as a **direct sum over angular octants** — the block
structure is streaming-induced, since each octant sweeps in a fixed
direction:

.. math::
   :label: streaming-inverse-direct-sum

   A^{-1} \;=\; \bigoplus_{\sigma \in \mathcal{O}} A^{-1}_{\sigma},
   \qquad
   \mathcal{O} \;=\; \{\sigma = (\mathrm{sgn}\,\mu_x,\,
                                  \mathrm{sgn}\,\mu_y) :
                       \sigma \neq (0,0)\}
                  \,\cup\, \{(0,0)\},

.. vv-status: streaming-inverse-direct-sum documented

acting on the octant-restricted tensor space :math:`(N_\sigma,\,n_x,\,n_y,\,n_g)`.
The direction-cosine partition (Eq. :eq:`octant-sign-predicate`) is
the predicate the
:class:`~orpheus.numerics.quadrature.Quadrature` class exposes as
its cached :attr:`~orpheus.numerics.quadrature.Quadrature.octants`
property — a tuple of
:class:`~orpheus.numerics.measure.DiscreteMeasurePartition`
entries realised by
:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
(see :ref:`tensorial-framing` and the
:doc:`/theory/foundations/discrete_measures` consumer table).

For each non-degenerate octant :math:`\sigma`, the action of
:math:`A^{-1}_\sigma` is a **forward substitution along a per-octant
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
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` primitives
described next.

.. _sweep-octant-dependency-graph:

Cartesian 2D: Octant Dependency Graph (Wave 2)
-----------------------------------------------

This section documents the **§15A.2 "upwind trace complex / causal
transport DAG / direction sweep ordering" primitive** as it lives in
:mod:`orpheus.sn.loss_representation.sweep_graph` after Wave 2 of the SN performance plan
(branch ``feature/sn-octant-sweep-graph``, closes Issue #4).  The
shipped architecture replaces the legacy per-ordinate ``for n in
range(N)`` loop in :func:`~orpheus.sn.loss_representation._sweep_jacobi` with
a per-octant batched dispatch, lifting the per-call ``_diag_cache``
build to mesh-time work, and isolating the per-cell DD algebra in the
discretization's pure kernel pair
(:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
/ :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`)
that LD / EC / Step closures can override later.

.. note::

   **Architecture history — the dispatch surface re-layered twice.**
   Wave 2 (the original closure of Issue #4) routed the sweep through
   a per-level *packet* (the ``SweepCellSlice`` dataclass) consumed by
   four direction×storage methods — ``update_batch`` / ``residual_batch``
   on the strategy (full-field) plus their ``apply_windowed`` /
   ``residual_windowed`` siblings on the graph.  S6.4(e) **collapsed
   that surface**: the four walk methods became TWO storage walks
   (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`,
   :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`)
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
   * - :class:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel`
     - :mod:`orpheus.sn.loss_representation.sweep_graph`
     - Frozen + slotted dataclass carrying one direction sign per
       spatial axis (``signs[axis] ∈ {-1, 0, +1}``) — a single type
       labels a 1-D (``(±1,)``), 2-D (``(±1, ±1)``), or 3-D octant.
       Hashable; used as the key in the per-shape graph family
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`
       (owned by the ``_DAGWavefront`` representation family since
       S6.4(c) — historically a mesh attribute).  An all-zero
       signature denotes the pure-:math:`z` degenerate octant — no
       graph is built for it
       (:attr:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel.streams` is
       ``False``).  The 3-D ``sign_z`` is dropped by the 2-D Cartesian
       orchestration: the in-plane sweep is invariant under the
       out-of-plane axis, so multiple ordinates with the same in-plane
       ``signs`` but different ``sign_z`` share a single graph instance.
   * - :class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` (+ its
       two storage walks)
     - :mod:`orpheus.sn.loss_representation.sweep_graph`
     - Frozen dataclass holding the per-octant topological levels
       (anti-diagonals) and the per-axis face-index offsets.  Built
       once per ``(shape, octant)`` pair in the
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`
       cache (S6.4(c); historically at mesh construction); reused
       across every source iteration / Krylov matvec / outer
       iteration.  Exposes TWO storage walks (S6.4(e)):
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`
       carries the COMPLETE per-axis interior face cochain (the
       verification-oracle policy);
       :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
       advances a rolling :math:`(d{-}1)`-frontier window (the
       production policy, ``O(N·n_g·∏ n_a)`` shrunk to
       ``O(N·n_g·∏_{a<d−1} n_a)`` backing).  The walk owns the level
       loop, the storage, and the per-level operand extraction; it
       dispatches the cell algebra to a level operation (next two rows).
   * - The level-operation pair ``_CellSolve`` / ``_CellResidual``
     - :mod:`orpheus.sn.loss_representation.sweep_graph`
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
       :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
       / :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
     - :mod:`orpheus.transport.spatial.scheme`
     - The **storage-free extension point** on the
       :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` ABC
       (S6.4(e); historically the ``SweepCellSlice``-packeted
       ``update_batch`` / ``residual_batch``).  Each takes the per-axis
       incoming face fluxes + streaming coefficients + the level's cross
       section and source and returns ``(psi_avg, psi_out)`` (solve) or
       ``(residual, psi_out)`` (apply) — PURE cell algebra, no
       gather/scatter (that is the walk's job).  Default raises
       :exc:`NotImplementedError` — additive capability, not a contract
       change.  :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
       overrides the pair; LD / EC / Step closures override it later to
       join the batched wavefront walks (their per-cell
       :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.update`
       stays the canonical reference contract).

Per-shape precompute pattern (family-owned since S6.4(c))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dependency graph is a **derived object** — the
``(shape × octant)`` joint property.  It depends only on cell topology
and the octant sign convention; it does **not** depend on fluxes,
sources, BCs, quadrature, cross sections, or iteration state.  So the
graph build is paid once **per spatial shape** in the cached accessor
:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.for_shape`, owned
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
   :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk`).

The closed-form precompute lives in
:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.from_cartesian`
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
   :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full` (full
   cochain) or
   :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`
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
   :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`
   /
   :meth:`~orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`.
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
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` docstring on
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
table at :doc:`/theory/foundations/discrete_measures` (octant partition consumed by SN
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
(:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`)
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

* **L2 regression** — existing ``tests/sn/verification/mms/test_mms_2d.py``,
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
  property delegates to.  See :doc:`/theory/foundations/discrete_measures`.
* Wave 1 :math:`R \circ \Lambda \circ M` Galerkin scattering
  composition — the parallel "metric knows its iterative structure"
  refactor for the scattering source build.  See
  :ref:`sn-scattering-fission-operators`.
* :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase` — the
  strategy ABC carrying the per-cell :meth:`update` reference contract
  and the storage-free batched kernel pair
  :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
  / :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
  (S6.4(e); was the ``SweepCellSlice``-packeted ``update_batch`` /
  ``residual_batch``).
* :class:`~orpheus.transport.spatial.diamond.DiamondDifference` — the only
  shipping closure that overrides the kernel pair; the reference for
  the bit-identity contract (pure cell algebra — the ONLY
  direction-aware math in the SN spatial stack since S6.4(e) lifted
  storage to the walk layer).
* :mod:`orpheus.sn.loss_representation.sweep_graph` — the two storage walks
  (:meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_full`,
  :meth:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed`)
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

.. note::

   **Superseded (coupled-block campaign step 6, R-6.1, 2026-07-12).**
   This section records the Wave-D Round-2 consolidation (#161) — the
   *first* unification of the four sweep paths under a single entry, the
   operator-free ``transport_sweep``.  That entry was itself retired at
   step 6; the dispatch is now the first-class ``LossRepresentation``
   polymorphism (:func:`~orpheus.sn.loss_representation.default_for` selects
   :class:`~orpheus.sn.loss_representation.CumprodScan` for the 1-D scan,
   the anti-hyperplane wavefront for multi-D), documented in
   :doc:`/theory/methods/sn/loss_representation`.  The Wave-D narrative below is preserved as
   the origin of that unification: read ``transport_sweep`` and the
   ``ReducedStreamingOperator`` boolean-dispatch code as the then-current
   implementation, not today's.

Wave D Round 2 of the SN reshape campaign (Issue #161) consolidated
the four pre-existing sweep paths (1-D cumprod / 2-D wavefront /
spherical / cylindrical) under one operator-free ``transport_sweep``
entry point that branched
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
:meth:`~orpheus.transport.spatial.scheme.DiscretizationScheme.update`:

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
:attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.scheme` (introduced in
this round as a constructor argument with default
:class:`~orpheus.transport.spatial.diamond.DiamondDifference`).  The
default reproduces the inlined sweep math bit-identically — every
regression snapshot at ``tests/sn/regression/snapshots/`` was
generated with DD and continues to match bit-for-bit when the
unified sweep dispatches via ``scheme.update(...)``.  See
:ref:`cell-update-strategies` for the strategy contract and
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` for the
DD scalar form.

Wave C-extension will ship :class:`Step`, :class:`LinearDiscontinuous`,
and :class:`ExponentialCharacteristic` as positivity-preserving /
higher-order alternatives; the unified dispatch infrastructure is
in place to receive them — users will pass
``scheme=LinearDiscontinuous()`` etc. at
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` construction.

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
:meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
/ :meth:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
on the strategy attached to the
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (Wave 2 of the SN
performance plan; closes Issue #4 — see
:ref:`sweep-octant-dependency-graph` for the full architecture and
:ref:`sweep-dispatch-relayering` for the S6.4(e) re-layering).
The "inlined DD math" formerly carried inside
:func:`~orpheus.sn.loss_representation._sweep_jacobi` was lifted into
:class:`~orpheus.transport.spatial.diamond.DiamondDifference` as a single
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
  then-``solution_to_angular_flux*`` codec and the
  matvec helpers consumed the
  :class:`~orpheus.geometry.boundary.BoundaryTraceLaw` instances on
  the :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` (Wave B Issue 7
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
the then-``solution_to_angular_flux_spherical`` codec
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
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`
* The 21 foundation tests at
  :file:`tests/sn/sweep/test_boundary_face_flux.py`

See :ref:`phase-c-sweep-frame-matvec` for the replacement
architecture. The Phase A subsection is preserved as historical
context for the empirical-defects-investigation reasoning chain.

.. _sn-pole-angular-closure-protocol:

The pole angular closure (Issue #168 Phase B)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note:: **Contract evolution (Issue #236 Phase 2 B2 → Issue #248).**
   This subsection originally introduced the angular-closure contract
   as a ``@runtime_checkable`` ``PoleAngularClosure`` **Protocol**
   (structural typing, Phase B). Issue #236 Phase 2 B2 retyped every
   production consumer (matvec / sweep / geometry / scheme /
   cell-balance) onto the
   :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
   **ABC** and made the three strategy methods
   (:meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.precompute_psi_state`,
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.cell_contribution`,
   ``angular_adjoint``) ``@abstractmethod`` on it. That left the
   Protocol orphaned and divergent — it carried the ``c_*``/``tau``
   accessors and a legacy ``__call__`` bundle but **not** the three
   strategy methods — so **Issue #248 deleted it** along with the dead
   ``__call__`` bundle (and its ``tau_mm`` argument) and the orphaned
   recurrence helpers. The ABC is now the **sole** angular-closure
   contract; production consumes
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.precompute_psi_state`
   + :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.cell_contribution`,
   never ``__call__``. The narrative below preserves Phase B's
   reasoning; read "strategy ABC" wherever the original said "strategy
   Protocol". The section anchor ``sn-pole-angular-closure-protocol``
   is retained (it is cross-referenced from
   :doc:`/theory/foundations/boundary_conditions` and elsewhere); only the human
   label "protocol" is now loose — the contract is an ABC.

Phase B addresses **Defect 3** of Issue #168 — the angular-redistribution
truncation gap on angularly-varying :math:`\psi`.  The pre-Phase-B
operator carried inline τ-symmetric interpolation
:math:`\psi_{n+1/2} \approx \tau_n\,\psi_{n+1} + (1-\tau_n)\,\psi_n`,
which is the **flat-flux collapse** of Hébert (2009) §3.9.4
Eqs. 3.428 + 3.437/3.439 — exact when :math:`\psi` is constant in
:math:`\mu`, but only :math:`\mathcal{O}(1)` accurate on
angularly-varying :math:`\psi`.  Phase B lifts this evaluation into
a :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
strategy ABC — analogous to Phase A's
``BoundaryFaceFlux`` —
and shipped **three concrete strategies** trading off bit-identity,
flat-flux invariance, and asymptotic accuracy:

.. note:: **Superseded strategy set (PR-TYPED-6c Step 7, 2026-05-18; #248,
   2026-06-18).** Of the three Phase-B strategies below, only
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
   survives — it is now the **sole production curvilinear closure and the
   default**. The two ablation strategies
   ``LegacyTauSymmetricInterpolation`` (the pre-Phase-B inlined
   :math:`\tau`-symmetric form) and ``BaileyFlatFluxRedist`` (the
   algebraic flat-flux collapse) were retired once
   ``MorelMontryAngularSweep`` became the default (no production
   consumer remained), and the divergent
   ``test_spherical_flat_flux_legacy_matches_bailey_collapse`` gate went
   with them. The three-strategy exploration is preserved below as the
   design record of *why* the M-M weighted-DD recurrence is the right
   closure.

* ``LegacyTauSymmetricInterpolation``
  — bit-for-bit reproduction of the pre-Phase-B inlined math
  (the :math:`\tau`-symmetric form).  Was the **default** through
  Phase B, preserving the
  curvilinear regression-snapshot bit-identity contract and the
  per-ordinate flat-flux invariant the ERR-026 evidence test relied on.
  Carried Defect 3 by design — the truncation gap was *reproducible*
  so verification probes could cross-check against the
  documented behaviour.

* ``BaileyFlatFluxRedist``
  — the algebraic flat-flux collapse
  :math:`R_{n,i,g} = (\Delta A/w)\,(\alpha_{n+1/2} - \alpha_{n-1/2})\,
  \psi_{n,i,g} / V_i = -\mu_n\,\Delta A_i\,\psi_{n,i,g} / V_i`
  (using :eq:`bailey-dome-recursion`).  Equivalent to the legacy form
  on flat :math:`\psi` (a now-retired flat-flux-identity test pinned
  this), and used as a structurally simpler bridge to the
  flat-flux invariant in unit tests.

* :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
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
  inside :class:`~orpheus.transport.spatial.diamond.DiamondDifference` (the
  sweep's cell update); applying this strategy in the apply matvec
  brings the apply and sweep to the same angular closure, but the
  **spatial** closures still differ (apply uses arithmetic averages
  + DD extrapolation; sweep uses WDD).  Full ERR-026 closure on the
  apply matvec requires aligning the spatial closure also (a
  follow-up beyond Phase B's scope — design memo §6.4 / §11).
  :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
  was therefore **opt-in** in Phase B (the default was then
  ``LegacyTauSymmetricInterpolation``); it has since become the sole
  production strategy and the default (see the supersession note above).

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
:mod:`orpheus.sn.sweep.pole_angular_closure` so the Hébert canonical
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
``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``), :mod:`orpheus.transport.spatial.diamond`, and the
new :mod:`orpheus.sn.sweep.pole_angular_closure` module.  Hébert
(2009) §3.9.4 is the **primary source** for the curvilinear S\ :sub:`N`
discretization in this codebase; Bailey-Morel-Chang 2010 is the
auxiliary justification for the M-M weighted-diamond :math:`\tau`
clamp.

ERR-026 closure status (Phase B partial)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Phase B ships the architectural infrastructure for closing Defect 3
(the :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
Protocol + canonical
:class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`)
without flipping the default — the canonical form in isolation does
not produce an :math:`\mathcal{O}(h^2)` MMS rate because the apply
matvec's spatial closure still differs from the sweep's WDD form.
The four ``xfail-strict`` curvilinear MMS tripwires therefore stay
xfail through Phase B; ERR-026 stays at **PARTIAL CLOSURE** (Phase A
closed Defects 1+2 spatial, Phase B ships Defect 3 architectural
scaffolding).  The full closure requires a Phase C follow-up that
aligns the apply matvec's spatial closure with the sweep's WDD
form, at which point
:class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
becomes the natural default.

.. _sn-pole-closure-compute-psi-half:

Half-angle grid exposure (Issue #197 PR-TYPED-6b)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. todo:: Archivist expansion needed.

   The Issue #197 PR-TYPED-6b dispatch added the public surface
   :func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
   exposing the M-M recurrence's half-angle grid
   :math:`\phi_{m\pm 1/2,i,g}` for one level.  (Originally an instance
   method on ``MorelMontryAngularSweep``, served by the unbound
   ``sn_mesh=None`` legacy mode; the C5 retirement of that mode,
   2026-07-03, moved it to module level — the surface takes all data
   via arguments and the seed strategy as a keyword, so it never
   needed an instance.)  It is the intermediate exposure that lets the
   unified SN matvec
   (:class:`~orpheus.sn.operators.streaming.StreamingOperator`) consume
   :math:`\phi_{m\pm 1/2}` as
   :func:`~orpheus.transport.spatial.cell_balance.cell_balance_for_streaming`'s
   ``psi_angular_upstream`` argument — closing the apply-vs-sweep
   twin path on the curvilinear angular branch.

   Pattern 2 (Single source of truth).  The :eq:`pole-mm-recurrence`
   recurrence body lives once, in the module-level kernel
   ``pole_angular_closure._psi_half_grid_single_level``.  Both the
   public ``compute_psi_half_per_level`` AND the live production path
   route through it: the matvec/sweep call
   :meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.precompute_psi_state`,
   which dispatches per level through ``_psi_half_grid_for_level``,
   which calls the same kernel.  (Phase B / Phase C drove the
   redistribution through a single ``__call__`` bundle that also routed
   through this kernel; Issue #248 deleted the bundle, so the
   single-source-of-truth invariant now binds
   ``compute_psi_half_per_level`` and ``precompute_psi_state``
   directly.)

   Test gate:
   :file:`tests/sn/sweep/curvilinear/test_compute_psi_half_per_level.py`
   — foundation + L0 tests pinning function existence, shape contract,
   the verbatim Hébert recurrence formula
   :math:`\phi_{m+1/2} = (\phi_m - (1-\tau_m)\phi_{m-1/2})/\tau_m`,
   Carlson-context seed contract, the Pattern-2 round-trip
   (``compute_psi_half_per_level`` against the
   ``_psi_half_grid_single_level`` kernel), and linearity. After
   Issue #248 these gates drive the recurrence through the **live**
   surface the matvec consumes, rather than the retired ``__call__``
   bundle.

   Closeout memo:
   ``.claude/agent-memory/method-implementer/issue_197_pr_typed_6b_closeout.md``.

.. _phase-c-sweep-frame-matvec:

Sweep-frame apply matvec (Issue #168 Phase C)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. admonition:: Key Facts
   :class: important

   * Phase C (commits ``eae6f05`` → ``d445a8f``, 2026-05-12) rewrote
     the then-production ``transport_operator_matvec_spherical``
     and ``_cylindrical`` matvecs (the whole per-geometry family
     since deleted — #197 / #280 campaigns) as **one sweep iteration
     semantically**.
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
     ``BoundaryFaceFlux``
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
     :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` vector
     defined by the realised BC operator, not extrapolated from
     interior cell centres.

Phase C resumes the spatial-closure alignment that Phase B
identified as the load-bearing precondition for the
curvilinear-default flip (see
:ref:`sn-pole-angular-closure-protocol` for Phase B's full
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
closure-contract narrative; this subsection picks up at the point
Phase B left for follow-up). Three unblockers landed before Phase C
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
   ``BoundaryFaceFlux``
   Protocol — built to second-order-accuratise the curvilinear
   outer face — was re-classified as **a patch on top of the wrong
   architecture**. Phase A's
   ``DDExtrapolation``
   produces a face-flux extrapolant
   :math:`\psi^{\text{face}}_{N-1/2} = \tfrac{3}{2}\,\psi_{N-1} -
   \tfrac{1}{2}\,\psi_{N-2}` that ignores the BC entirely; Phase B's
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
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
:class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
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
:meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk` invoked with
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
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
strategy and the collision term :math:`\Sigma_t \psi^{\text{cell}}_{n,i,g}`
unchanged. The full per-cell update is

.. math::
   :label: phase-c-cell-update

   (T\psi)_{n,i,g} \;=\; S_{n,i,g} \;+\; R_{n,i,g} \;+\;
   \Sigma_t(i,g)\,\psi^{\text{cell}}_{n,i,g},

where :math:`R_{n,i,g}` is the strategy's redistribution output
(Bailey :math:`\Delta A_i / w_n` redistribution factor with the
strategy-specific :math:`\alpha_{n+1/2}\psi_{n+1/2} -
\alpha_{n-1/2}\psi_{n-1/2}` evaluation). In current production this
is produced by
:meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.cell_contribution`,
consuming the per-level state that
:meth:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase.precompute_psi_state`
stamps once per sweep. (Phase B / Phase C shipped this redistribution
through a single ``__call__`` bundle on the strategy; that legacy
bundle — and its ``tau_mm`` argument — was retired in Issue #248, and
the two strategy methods are now the sole production surface.) See
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

* :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk`
  (``dag_walk(*, ordinate_idx=..., direction_sign=..., mu_level_idx=None)``)
  — Issue #196 Phase G
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
  (:file:`tests/sn/sweep/core/test_dag_walk.py`) pins
  bit-identity between the two invocation modes across sphere /
  slab / cylindrical for every representative ordinate. For
  cylindrical the per-level
  ``mu_level_idx`` is required (the within-level DAG topology
  differs per level).

* ``EquationMap.unknowns_at_cell_for_mask(cell_idx, ordinate_mask)``
  — a precomputed inverse lookup ``(cell, ordinate) → k``. Lazily
  builds an ``(nx, N) int`` table with :math:`-1` sentinels for
  absent ``(ordinate, cell)`` slots; subsequent calls are O(1) per
  ``(cell, mask)`` pair. Replaces the per-equation O(n_eq) linear
  scan the legacy scatter pattern used. The eq_map still iterates
  ``(spatial outer, ordinate inner)`` at construction time; the
  helper just precomputes the inverse for the sweep-frame matvec.
  (The whole packed-vector ``EquationMap`` codec was subsequently
  retired at D-J — 2026-05-30 — when the typed :class:`FullField`
  composite replaced the packed-vector convention; this bullet
  records the Phase G design as it stood.)

What retires
^^^^^^^^^^^^

Phase A's
``BoundaryFaceFlux``
Protocol — five symbols, the
:class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` field, and the 21 foundation
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
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` attribute
* The ``boundary_face_flux_closure`` keyword argument from
  ``transport_operator_matvec_spherical`` and ``_cylindrical`` (the
  matvec family since deleted — #197 / #280 campaigns)
* :file:`tests/sn/sweep/test_boundary_face_flux.py` (232 LOC,
  21 foundation tests)

Three additional simplifications shipped with the rewrite:

* the then-``solution_to_angular_flux_spherical`` codec
  (and its cylindrical alias) returned a single ``fi`` array
  ``(ng, N, nx, 1)`` instead of the Phase A
  ``(fi, boundary_face_flux)`` tuple. Inward-at-boundary cell-centre
  slots ``fi[:, n_inward, -1, 0]`` are filled with the
  **reflected-partner cell-centre value** as an analytical
  extension: the equation map excludes these from unknowns (the BC
  determines them), but the WDD recurrence on flat :math:`\psi`
  requires the cell-centre to be consistent so the per-ordinate
  flat-flux invariant holds.
* :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` no longer accepts the
  ``boundary_face_flux=`` keyword (a regression test pins the
  field retirement).
* :class:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` dispatch
  drops the ``boundary_face_flux_closure`` plumbing.

What stays
^^^^^^^^^^

* The Phase B
  :class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
  closure contract (the ABC, since Issue #248) stays. The sphere
  centre / cylinder axis is **intrinsic geometry** (a coordinate-system
  singularity, not an external BC), so the three-strategy angular
  closure is the right shape; only the **default** is under question,
  and that is the Phase D decision point.
* The :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply_transpose`
  machinery via dense-probe construction stays. Linearity of the
  rewritten :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply`
  (Gate 1.4, pinned to ``rtol=1e-13``) guarantees the transpose is
  correctly tracked.
* The
  :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` +
  :class:`~orpheus.sn.mesh.method_space.SNMethodSpace` +
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
\psi^{\text{cell}}(\text{first cell at level})`. For a
level-symmetric quadrature the cylinder tolerates a wrong seed via
the **dead first-ordinate weight** (:math:`c_{\rm in}[m_0]=0`) in a
way the spherical case does not — see the Gate 1.1 finding below
(and its #280 Phase 2.5b correction: this is level-symmetric-only,
NOT :math:`\alpha`-dome telescoping, and false for a product
quadrature).

.. _phase-c-apply-sweep-equivalence:

apply ↔ sweep structural equivalence
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pre-Phase-C, the matvec's :meth:`apply` and the sweep's :meth:`solve`
were structurally distinct paths to the **same** discrete operator
:math:`A`: :meth:`apply` walked cell-centre storage with arithmetic
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
   :meth:`~orpheus.sn.mesh.augmented_mesh.SNMesh.dag_walk` invoked with
   ``direction_sign=±1``.
#. The BC trace law applied **once** at the boundary edge per
   :ref:`affine-bc-form`.

The face-flux propagation identity therefore holds **by
construction** post-Phase-C: extracting the implicit face fluxes
from ``apply`` (by inverting the cell-balance equation) recovers
the same WDD recurrence the sweep walks. The structural-frame
identity is the load-bearing acceptance criterion for
preconditioned-Krylov stability — when ``apply`` is the operator
:math:`A` and the sweep is :math:`A^{-1}` (approximately), they
must agree on what :math:`A` **is**. Foundation tests in
:file:`tests/sn/test_phase_c_gates.py` pin this via
``np.array_equal`` on:

* **Gate 1.2** (apply determinism) — repeat-call invariance
  ``apply(ψ) == apply(ψ)`` bit-identical across two invocations.
* **Gate 1.3** (apply ↔ apply_transpose reciprocity)
  :math:`\langle A\psi, \phi \rangle = \langle \psi, A^T\phi \rangle`
  to ``rtol=1e-12, atol=1e-13``. Free if Gate 1.4 (linearity)
  passes.
* **Gate 1.4** (apply linearity) :math:`A(\alpha\psi + \beta\phi)
  = \alpha A\psi + \beta A\phi` to ``rtol=1e-13``.
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
   :doc:`/theory/foundations/operator_algebra`.

The :doc:`/theory/foundations/boundary_conditions` :ref:`affine-bc-form` (§16A.3) reads

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
:meth:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer.realize`
internally consumes the outflow selector
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`
on the unified :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
(the ordinate slots with :math:`\mu \cdot \hat n > 0` at the face)
and writes only the inflow slots; the outflow slots in the output
are unspecified by the §16A.3 contract and the matvec reads back
only ``inflow_full[incoming_mask, :]``. This is the user's "ghost
cell for higher-order boundary closure" idiom realised as a typed
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` vector
defined by the realised BC operator — not extrapolated from
interior cell centres.

**Gate 1.5** (foundation,
:func:`tests.sn.sweep.core.test_phase_c_gates.test_bc_trace_contract_respected_by_matvec_vacuum_sphere` /
:func:`tests.sn.sweep.core.test_phase_c_gates.test_bc_trace_contract_respected_by_matvec_reflective_sphere`)
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
"keystone") is **deleted**; the 1-D ``transport_sweep`` entry (then the
production sweep, since retired at step 6) read the *seeded* inflow trace
directly:

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
  :meth:`InvertibleOperator._solve_timed_full_field <orpheus.sn.operators.streaming.InvertibleOperator._solve_timed_full_field>`
  seeds its boundary buffer from :math:`\text{rhs.boundary}`, **not**
  from the iterate ``initial_guess.boundary`` (the retired
  partner-flux carrier).
* **Direct loops** (direct fixed-source SI, final eigenvalue
  reconstruction sweep): the helper
  :func:`~orpheus.sn.solver._reflect_outflow_into_inflow` fills the
  inflow slots with :math:`B\,\psi.\text{outflow}` in place before the
  sweep, via the canonical
  :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`.

Both routes call the **identical**
:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` — :math:`B`
is single-sourced. For vacuum :math:`B = 0`, so the bare sweep reads a
zero inflow seed and the result is **bit-identical** to the
pre-extraction ``bc.apply`` of a vacuum law.

The full block-matrix derivation, the three design corrections (keep
the outflow defect; project :math:`B` to the inflow row; seed from
:math:`\text{rhs.boundary}`), the two delivery routes, and the O.2
forcing function all live at :ref:`bc-extraction` in
:doc:`/theory/foundations/operator_algebra` — the canonical home for the SN transport
operator algebra.

Step **O.4b** extended the bare sweep to the **2-D Cartesian
wavefront** path (both :func:`~orpheus.sn.loss_representation._sweep_jacobi`
and the 2-D matvec
:meth:`StreamingOperator._apply_2d_cartesian <orpheus.sn.operators.streaming.StreamingOperator>`):
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
:doc:`/theory/foundations/wavefront_cochain`).

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
``LegacyTauSymmetricInterpolation``
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
(:file:`tests/sn/verification/mms/test_mms_curvilinear.py` and the L1 aniso file)
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

The cross-check placeholder landed as the Phase D test
:func:`tests.sn.verification.analytical.test_phase_c_crosscheck.test_phase_d_trajectory_resolvent_crosscheck`
(after the pole-face spatial-closure refinement). It is
**structurally important**: it pins
the names of the bare entry points so the reader knows
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
(:file:`tests/sn/sweep/curvilinear/test_pole_angular_closure.py`); the L1
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
(:class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`):

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

  .. warning::

     **Level-symmetric-only (corrected #280 Phase 2.5b).**  This
     "the cylinder absorbs the seed discrepancy" mechanism holds
     ONLY for a level-symmetric quadrature, where the first-swept
     ordinate's seed weight is exactly zero
     (:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`)
     — a **dead** seed annihilated at source, not a cancellation
     across the azimuthal cascade.  For a **product** quadrature
     the starting direction coincides with the first-swept ordinate
     (:math:`t=0`, #229), so :math:`c_{\rm in}[m_0]\ne 0`, the seed
     is a **live self-coupling**, and the cold cylinder
     ``(L+C).solve`` was seed-**lagged** until the #280 2.5b
     direct-seed fold.  See the ERR-026 crosstab correction note.

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
``BoundaryFaceFlux``
Protocol — but per the empirical Gate 1.1 finding above the
curvilinear default stays
``LegacyTauSymmetricInterpolation``
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
   ``CarlsonInwardSweep`` *half-angle* seed described here is also
   superseded as the default by
   ``AngularEdgeExtrapolation``.
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
   ``LegacyTauSymmetricInterpolation``
   → :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`;
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

.. attention:: **Superseded by Issue #282 route (a) (2026-07-04).**

   The swappable ``PsiHalfAngleSeed`` strategy family
   (``ZeroSeed`` / ``CarlsonInwardSweep`` / ``AngularEdgeExtrapolation``)
   whose design this section — and the Phase F and ERR-058 sections that
   follow — build up was **retired** by Issue #282 route (a).  The
   starting-direction half-angle flux :math:`\psi_{1/2}` is now
   **first-class typed state** the sweep marches *directly* from the true
   within-group source, not a functional of the previous iterate.  Any
   "current default / retained strategy / the seed lives as a strategy
   field" claim in the three sections below is **historical** — read them
   for the *why* (what was tried and the diagnoses that narrowed the
   defect), but for the CURRENT design see
   :ref:`sn-282-direct-starting-direction-solve`.  In particular the
   ``AngularEdgeExtrapolation`` "iterate extrapolation" seed those
   sections land on was itself the #282 walk-order back edge that route
   (a) removes.

.. admonition:: Key Facts
   :class: important

   * Phase D (commit landed 2026-05-12 on
     ``refactor/sn-operator-algebra``) closes the structural bug
     in :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
     by replacing the hardcoded ``psi_half_left = 0`` seed with
     the canonical Hébert §3.9.4 Eqs. (3.432)–(3.435) inward
     :math:`\mu = -1` sweep output.
   * The seed lives in the **M-M angular recurrence**
     (:func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`
     in ``orpheus/sn/sweep/pole_angular_closure.py``), **NOT**
     in the WDD spatial pole-face initial condition the
     :ref:`Phase C plan <sn-curvilinear-trajectory-resolvent-crosscheck-section>`
     proposed.  The diagnostic memo at
     ``.claude/agent-memory/numerics-investigator/phase_d_gate_1_1_sphere_mms_diagnosis.md``
     empirically falsified intervention ``[A]`` (WDD pole-face
     replacement) and confirmed intervention ``[B]`` (M-M
     half-angle seed replacement).
   * Architectural choice is **Option α (composition)**: the
     seed lives as a
     ``PsiHalfAngleSeed``
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
moment of the input :math:`\psi`.  The :math:`L = 0` collapse is the
*apply matvec's* reach (isotropic scattering); the **source** side is
NOT collapsed — Issue #282 route (a) folds **all** Legendre moments of
the true within-group source, because streaming manufactures angular
structure an isotropic flux does not have (an :math:`\ell = 0`-only
fold floored the anisotropic curvilinear MMS).  See
:ref:`sn-282-source-fold` for the full fold and the load-bearing
:math:`\ell = 1` term.

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
   :mod:`~orpheus.sn.sweep.psi_half_angle_seed` module docstring
   (the canonical algebra-of-record).  Each :math:`:label:` is
   unique across the documentation graph; the Sphinx page is the
   **presentation layer** for the equations the code module owns
   as source-of-truth.  The
   ``@pytest.mark.verifies("hebert-3-43X")`` wiring on the L0
   algebraic-identity tests in
   :file:`tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py` is tracked
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
:func:`tests.sn.sweep.curvilinear.test_psi_half_angle_seed.TestCarlsonFlatPsiAlgebraicIdentity.test_carlson_flat_psi_identity_reflective`
pins this identity at machine precision (``rtol=1e-13``).

The corrected injection-point story
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The single largest **architectural correction** of Phase D is
where the canonical inward-sweep output is injected.  The Phase D
plan (and the literature memo's §7 implementation note) routed
the inward-sweep result :math:`\bar\phi_i` into the **WDD
spatial pole-face initial condition** at the then-production
``transport_operator_matvec_spherical`` matvec's (since deleted)
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
:func:`tests.sn.sweep.curvilinear.test_psi_half_angle_seed.TestCarlsonFlatPsiAlgebraicIdentity.test_carlson_vacuum_BC_flat_source_nx_3`
— a vacuum-BC hand calculation on the Carlson inward sweep
(``rtol=1e-13``) whose values are distinct from the degenerate
broadcast-cell-centre seed.  Without this test a future
regression that replaced the Carlson sweep with a naive
broadcast-cell-centre would pass every flat-:math:`\psi`
reflective test silently.

The bug Phase B baked in
~~~~~~~~~~~~~~~~~~~~~~~~

The pre-Phase-D production code at
``orpheus/sn/sweep/pole_angular_closure.py:411`` carried the
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
**passes** empirically.  The mechanism is the **dead first-ordinate
seed** of the level-symmetric quadrature exercised here: the
first-swept ordinate's seed weight is zero
(:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`), so
the wrong ``psi_half_left = 0`` seed is annihilated at source per
level.  (This was originally read as per-:math:`\mu`-level
:math:`\alpha`-dome telescoping "cancelling" the seed; #280 Phase
2.5b corrected that — it is a dead weight, level-symmetric-only, and
**false for a product quadrature**, where the seed is a live
self-coupling and the cold solve was seed-lagged until the
direct-seed fold.)  The sphere cascade has no equivalent dead-seed
weight — a wrong seed propagates directly to a wrong fixed point.  Phase D's fix
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
:class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
carries a ``psi_half_left`` variable to seed.  The Legacy and
Bailey closures don't have one — their half-angle face flux
evaluation collapses to cell-centre values unconditionally.  Two
architectures were considered:

* **Option α (composition, shipped)** — the seed strategy lives as
  a ``PsiHalfAngleSeed``
  field on
  :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`.
  The abstraction stays local to the closure that consumes it.

* **Option B (sibling Protocol on SNMesh, rejected)** — the seed
  would be a separate Protocol attribute on
  :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh`, applied by the matvec
  before calling the pole closure.  This would force every
  consumer (Legacy / BFF / M-M) to handle a Protocol that is a
  **no-op** for the non-M-M strategies, violating the
  single-responsibility principle and forcing unrelated tests to
  thread the Protocol through call signatures.

The
``CarlsonSweepContext``
dataclass bundles the four inputs the Carlson sweep needs that
are NOT in the
:class:`~orpheus.sn.sweep.pole_angular_closure.PoleAngularClosureBase`
strategy's ordinary per-cell call signature (``sigma_t``, ``dr``,
``mu_quad``, ``weights``, ``bc_outer_value``), keeping the
call-signature expansion to a single new optional keyword — a minimal
blast-radius extension that Legacy and Bailey closures ignore by
documented closure contract.

Linear-operator preservation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Both seed strategies — ``ZeroSeed`` and
``CarlsonInwardSweep`` — are **linear in the input** ``psi_cells``
(verified by the ``is_linear: ClassVar[bool] = True`` trait, pinned
by foundation tests).  Linearity is the load-bearing property:
the apply matvec must be a linear operator, otherwise the
operator-algebra operations of
:class:`~orpheus.sn.operators.streaming.InvertibleOperator`
(apply, apply_transpose, dense matrix probing) break.  The
``CarlsonInwardSweep`` is linear because:

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
:func:`tests.sn.sweep.curvilinear.test_psi_half_angle_seed.TestSeedLinearity.test_carlson_inward_sweep_is_linear`
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
``CarlsonInwardSweep``
evaluates only the :math:`\ell = 0` (isotropic) Legendre moment
when building the moment-folded source in
:eq:`hebert-3-432-source`.  This is **consistent with the apply
matvec's structure**: the
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` apply matvec
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
   <orpheus.sn.mesh.augmented_mesh.SNMesh.pole_angular_closure>` default
   flipped from
   ``LegacyTauSymmetricInterpolation``
   to
   :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`.
   :class:`MorelMontryAngularSweep`'s own constructor default for
   ``psi_half_seed`` is
   ``CarlsonInwardSweep``,
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
   :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply`.  The
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
:func:`tests.sn.sweep.core.test_phase_c_gates.test_apply_curvilinear_per_ordinate_flat_flux_residual`
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
already passed under Phase C because — for the **level-symmetric**
quadrature exercised here — the first-swept ordinate's seed weight
is exactly zero (:math:`c_{\rm in}[m_0] = (1-\tau)/\tau = 0` at raw
:math:`\tau = 1`), so the zero-seed inconsistency was annihilated at
source (a **dead** seed), not "absorbed" by any telescoping of the
solve.

.. note::

   **Correction (#280 Phase 2.5b, 2026-07-05).**  An earlier reading
   of this crosstab attributed the cylinder's Phase-C pass to
   ":math:`\alpha`-dome telescoping absorbing the zero-seed
   inconsistency" and generalised it to "the cylinder solve is
   seed-insensitive / was already exact."  That is a **level-symmetric-
   only** artefact and is **false for a product quadrature**: there the
   starting direction coincides with the first-swept ordinate
   (:math:`\mu_{\rm start} \equiv \mu_{m_0}`, :math:`t = 0`, #229), so
   :math:`c_{\rm in}[m_0] \ne 0` and the seed is a **live per-ordinate
   self-coupling** that contributes :math:`O(1)` to the :math:`m_0`
   cell diagonal.  The product-cylinder cold ``(L+C).solve`` was in
   fact seed-**lagged** (cold error :math:`\approx 0.57`) until the
   #280 2.5b direct-seed fold folded that self-coupling into the
   :math:`m_0` diagonal (:math:`c_{\rm out} \to c_{\rm out} -
   c_{\rm in}`), making the cold solve a single-pass direct inverse.
   The augmented :math:`(L+C)` is block-lower-triangular because the
   seed contribution lands **on the block diagonal** (forward
   substitution resolves it) — *not* because the seed "telescopes
   away."  Distinct claim, still valid: the :math:`\alpha`-dome
   telescopes under the angular weight sum
   :math:`\sum_n w_n \psi_n`, which is why **scalar / balance** V&V
   gates are blind to a wrong per-ordinate seed (anti-pattern #8) —
   that blindness statement is unaffected by this correction.

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
   ``CarlsonSweepContext``.  See the BC companion section
   :ref:`bc-phase-d-two-bc-applies-per-matvec`.
#. **Phase C BC trace law call** — applied to the WDD-propagated
   outflow face value at the boundary edge, per the
   :ref:`affine-bc-form` contract.

The capture-and-compare test
:func:`tests.sn.sweep.core.test_phase_c_gates.test_bc_trace_contract_capture_and_compare_sphere`
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

* :mod:`orpheus.sn.sweep.psi_half_angle_seed` — Protocol family
  + ABC + 2 strategies (``ZeroSeed`` + ``CarlsonInwardSweep``)
  + ``CarlsonSweepContext`` dataclass.
* :file:`tests/sn/sweep/curvilinear/test_psi_half_angle_seed.py` — 25
  foundation + L0 + L1 tests covering Protocol conformance,
  registry/self-registration, immutability, shape contract,
  bit-identity for ``ZeroSeed``, L0 algebraic identities
  (flat-:math:`\psi` at varying C, vacuum-BC nx=3 hand
  calculation, multi-region :math:`\Sigma_t` step), linearity, and
  L1 structural-independence (Carlson vs Zero on vacuum-BC probe).

**Modified files**

* :mod:`orpheus.sn.sweep.pole_angular_closure` —
  :class:`MorelMontryAngularSweep` gains a
  ``psi_half_seed: PsiHalfAngleSeed`` field; the per-level M-M
  recurrence (then ``_mm_weighted_angular_recurrence_single_level``,
  now :func:`~orpheus.sn.sweep.pole_angular_closure.compute_psi_half_per_level`)
  accepts an
  optional ``psi_half_seed`` array; Protocol signatures extended
  with an optional ``carlson_context`` kwarg (Legacy + Bailey
  ignore it).
* :mod:`orpheus.sn.sweep` ``__init__`` re-exports the new
  symbols.
* :mod:`orpheus.sn.operators.streaming` — spherical + cylindrical matvecs
  build the
  ``CarlsonSweepContext``
  before calling ``pole_angular_closure``.
* :mod:`orpheus.sn.mesh.augmented_mesh` — :class:`SNMesh` default flipped to
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
     (``CarlsonInwardSweep``,
     Hébert §3.9.4 Eqs. (3.432)–(3.435)) from the apply-matvec path
     (the then-production ``transport_operator_matvec_spherical``
     / ``_cylindrical`` matvec — since deleted, #197 / #280 campaigns —
     fixed in Phase D Step 3) into the SI/sweep
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
     (``orpheus/sn/sweep/pole_angular_closure.py:411``).
   * Phase F factors a **NEW free function**
     :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
     ``(Q_bar, sigma_t, dr, bc_outer_value) -> (ng, nx)`` that runs
     :eq:`hebert-3-434`–:eq:`hebert-3-435` driven by the SI
     within-group source ``Q_1d`` rather than by an apply-path
     :math:`\psi` Legendre fold.
     ``CarlsonInwardSweep.__call__`` is refactored to delegate
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
:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` via the
``MorelMontryAngularSweep.psi_half_seed`` composition; that
covers every Krylov-driven call.  But ORPHEUS's curvilinear
production default is **source iteration**, which (pre-Phase-F) dispatched
through the then-production ``transport_sweep`` entry rather than
through the apply matvec, and the two paths ran **different code**
to seed the M-M half-angle recurrence:

.. list-table:: Apply vs SI/sweep dispatch divergence (pre-Phase-F)
   :header-rows: 1
   :widths: 24 38 38

   * - Path
     - Carlson seed site
     - Pre-Phase-F state
   * - Apply matvec (Krylov)
     - ``_mm_weighted_angular_recurrence_single_level``
       :math:`\to`
       ``CarlsonInwardSweep``
       via ``MorelMontryAngularSweep.psi_half_seed``
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
because — for the **level-symmetric** quadrature exercised here — the
first-swept ordinate's seed weight is the **dead** first-ordinate weight
(:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`), which
annihilates the wrong seed at source.  (This is **not**
:math:`\alpha`-dome ':math:`\alpha=0`' level-edge cancellation
*absorbing* the seed — a level-symmetric-only reading, false for a
product quadrature; see the #280 Phase 2.5b correction at
:ref:`sn-phase-d-gate-1-1-empirical`.)  Cardinal Rule 2 (architecture)
nonetheless demands the structural fix on the sister path even when the
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
The Krylov path went through
:meth:`InvertibleOperator.apply` (which consumes the Phase D
Carlson seed correctly); the SI path went through the then-production
``transport_sweep`` entry (which carried the **legacy zero
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
``CarlsonInwardSweep``
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
:func:`tests.sn.sweep.core.test_sweep_vs_apply_consistency`
pins the source-vs-:math:`\psi` equivalence directly: for any
flat-:math:`\psi` field ``ψ_const`` with ``bc_outer_value =
ψ_const`` (reflective) and ``Q_1d = Σ_t · Σw · ψ_const`` (the
within-group source built by SI from
``φ_0 = Σw · ψ_const``), the two helpers return
**bit-identical seeds** (up to FP non-associativity).
Apply-path:
``CarlsonInwardSweep``
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
``CarlsonInwardSweep.__call__``.  Naive options for the
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
  ``CarlsonInwardSweep``
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

``CarlsonInwardSweep.__call__`` now reads (in essence):

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

Cylindrical Gate 1.1 passes empirically pre-Phase-F (for the
**level-symmetric** quadrature exercised here the first-swept
ordinate's seed weight is the **dead** first-ordinate weight
:math:`c_{\rm in}[m_0]=(1-\tau)/\tau=0` at raw :math:`\tau=1`, which
annihilates the wrong zero seed at source — **not** :math:`\alpha`-dome
telescoping "absorbing" it, a level-symmetric-only reading false for a
product quadrature; see the :ref:`sn-phase-d-gate-1-1-empirical`
discussion and its #280 Phase 2.5b correction for the
sphere-vs-cylinder asymmetry).  Phase F nonetheless fixes
both sites for two reasons:

#. **Cardinal Rule 2 (architecture)**: structural alignment of
   the canonical math at both sites prevents a future
   refactor from introducing an asymmetric bug that only the
   sphere catches.  The sweep-path helper is the same code
   regardless of geometry; consuming it consistently from
   both geometries is the architecturally clean choice.
#. **Defense in depth against future stress probes**: on any
   cylinder rule where the first-ordinate seed weight is **live**
   (a **product** quadrature already is — :math:`c_{\rm in}[m_0]\ne 0`,
   #280 Phase 2.5b), the wrong zero seed enters the fixed point.
   Fixing both sites now is cheap insurance.

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

* :mod:`orpheus.sn.sweep.psi_half_angle_seed` — NEW free
  function :func:`carlson_inward_sweep_from_source` (lines
  358–419 of the module); ``CarlsonInwardSweep.__call__``
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

* :func:`tests.sn.sweep.core.test_phase_c_gates.test_sweep_curvilinear_per_ordinate_flat_flux_residual`
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
* :file:`tests/sn/sweep/core/test_sweep_vs_apply_consistency.py` —
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
     ``AngularEdgeExtrapolation``
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
:func:`tests.sn.sweep.core.test_phase_c_gates.test_sweep_curvilinear_per_ordinate_flat_flux_residual`
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
  :doc:`/theory/foundations/boundary_conditions` — extends the Phase D
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
     ``AngularEdgeExtrapolation``
     (the new
     :class:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep`
     ``psi_half_seed`` default).  It replaces
     ``CarlsonInwardSweep``,
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

   ``CarlsonInwardSweep`` is **retained** (not deleted) as the
   registered host of the canonical Hébert §3.9.4 recurrence
   (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`),
   reachable only by explicit opt-in.  Its class docstring carries a
   ``.. warning::`` block recording the proxy-source caveat by design,
   so a future session cannot re-activate it as a default unaware of
   the falsification.

   .. note:: **Retraction (2026-07-04, Issue #282 route (a)).**  The
      ``CarlsonInwardSweep`` *strategy class* was NOT ultimately
      retained — route (a) deleted the whole ``PsiHalfAngleSeed``
      family.  What survives is the pure Hébert recurrence
      :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
      (a free function, now the SOLVE engine driven by the **true** q½
      source rather than the falsified proxy, no opt-in), plus the
      inlined
      :meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`
      for non-carrying cylinder levels.  See
      :ref:`sn-282-seed-strategy-zoo`.

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
the SI sweep twin (:meth:`~orpheus.sn.loss_representation._OneDimScanWalk.sweep`).
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
Phase D / Phase F ``CarlsonInwardSweep`` solved the canonical
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
``AngularEdgeExtrapolation``.
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
``CarlsonSweepContext.mu_start``
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
* **linear in the input**, so the operator-algebra operations
  (:meth:`apply`, :meth:`apply_transpose`, dense probing) are
  preserved.  The strategy OWNS its adjoint
  (``PsiHalfAngleSeedBase.seed_adjoint``),
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
   * - ``AngularEdgeExtrapolation``
     - **production**
     - The new ``psi_half_seed`` default; operator-consistent on
       non-flat fields.
   * - ``CarlsonInwardSweep``
       + :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
     - retained, opt-in
     - Correct *source-driven* Hébert §3.9.4 recurrence; would seed a
       future TRUE-source sweep-side closure.  Proxy-source caveat
       pinned in its docstring ``warning`` block.
   * - ``ZeroSeed``
     - retained, ablation
     - Reproduces the Phase B behaviour for A/B regression-safety
       comparison.
   * - Coupled-pole spatial seed in
       :meth:`~orpheus.sn.loss_representation._OneDimScanWalk._apply_walk`
       / :meth:`~orpheus.sn.loss_representation._OneDimScanWalk.sweep`
     - **production**
     - The :math:`\psi(0,+\mu)=\psi(0,-\mu)` continuity; matvec + sweep
       share it (one discrete system).

.. note:: **Retraction (2026-07-04, Issue #282 route (a)).**  The three
   ``PsiHalfAngleSeed`` strategy rows above (``AngularEdgeExtrapolation``
   / ``CarlsonInwardSweep`` / ``ZeroSeed``) are superseded: route (a)
   **deleted** the whole strategy family — including the
   ``AngularEdgeExtrapolation`` "production default", which was itself
   the #282 walk-order back edge.  What is genuinely retained is the
   free-function Hébert engine
   :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
   (now the SOLVE driver, on the **true** q½ source) and the inlined
   :meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`
   (non-carrying cylinder levels).  The coupled-pole spatial seed row is
   unaffected.  See :ref:`sn-282-seed-strategy-zoo`.

Open research paths
~~~~~~~~~~~~~~~~~~~~~

Two paths could lift the anisotropic angular floor without changing the
isotropic O(h²) result:

#. **TRUE-source-driven sweep-side seed** — **LANDED as Issue #282 route
   (a) (2026-07-04).**  This path predicted the resolution exactly:
   replace the ``AngularEdgeExtrapolation`` *input-field* extrapolation
   with the canonical Hébert recurrence
   :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
   driven by the genuine within-group source
   :math:`\bar Q_i = \sum_\ell \tfrac{2\ell+1}{2}Q_\ell(r_i)(-1)^\ell`
   (the full Legendre fold, not the :math:`\Sigma_t\phi_0` proxy), making
   the *starting-direction transport* exact rather than the *input-field
   value* exact.  Route (a) shipped precisely this — and the "likely
   diagnostic probe" proposed here (holding spatial mesh fixed and
   sweeping quadrature order) is exactly the **angular-order N-sweep**
   that certified the re-pose principled.  The one refinement over the
   prediction: the seed also had to become **first-class typed state**
   (not just a better strategy) to kill the walk-order back edge that
   made the *solve* non-direct.  See
   :ref:`sn-282-direct-starting-direction-solve`.
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
  ``AngularEdgeExtrapolation``
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
  (``A.solve`` vs Krylov-on-``apply``) realises the *same* :math:`A^{-1}`
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


.. _sn-282-direct-starting-direction-solve:

Route (a) — the direct starting-direction ψ½ solve (Issue #282)
---------------------------------------------------------------

.. admonition:: Status banner
   :class: important

   **Issue #282 — CLOSED by route (a)** (#280 Phase 2.5d, committed
   ``a29ab2d`` on branch ``refactor/sn-walk-unification``, 2026-07-04).
   The lagged Morel–Montry half-angle pole seed — a two-point
   extrapolation of the *previous* source-iteration iterate — was a
   **walk-order back edge**: it made the spherical within-group SOLVE a
   *non-direct* inverse.  Route (a) promotes the starting-direction flux
   :math:`\psi_{1/2}` to **first-class typed state** — System B's own
   :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
   composite (a :math:`V_{\rm cell}`-state-metric ``interior ⊕ boundary``
   of ψ½ role leaves on the split
   :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicInteriorSpace`
   / :class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicBoundarySpace`)
   — and marches it **directly** from the true within-group source
   :math:`\bar q_{1/2}` through System B's named resolvent
   :meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
   (the Hébert §3.9.4 :eq:`hebert-3-434`–:eq:`hebert-3-435` recurrence
   :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`).

   **Keystone:** the sphere cold-start residual
   :math:`\lVert A\cdot\mathrm{solve}(b) - b\rVert_\infty / \lVert b\rVert_\infty`
   collapses from :math:`5.18\times10^{5}` to
   :math:`2.5\times10^{-16}`, and the seed-insensitivity
   :math:`\Delta` from :math:`4.57\times10^{-2}` to :math:`0` **bitwise**.
   The cold solve is now a genuine single-pass exact inverse — the
   posture the DSA program (#2) and the curvilinear Krylov
   preconditioner (#200) require, and the deliverable that lets the
   #280 unified walk build a spherical ``sweep_transpose`` against a
   triangular forward operator.

   This section is the **resolution chapter** of the seed-strategy saga
   documented above (:ref:`sn-phase-d-carlson-coupled-pole-sweep`,
   :ref:`sn-phase-f-carlson-sweep-path-backport`,
   :ref:`sn-err-058-closure-seed-closeout`).  Those sections are
   preserved as the record of *what was tried and why it fell short*;
   the ``PsiHalfAngleSeed`` strategy family they built is **retired**
   here (see :ref:`sn-282-seed-strategy-zoo`).

The lagged pole seed was a walk-order back edge (the #282 defect)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The curvilinear Morel–Montry angular recurrence :eq:`pole-mm-recurrence`
marches half-angle face fluxes :math:`\phi_{n+1/2,i}` up the
:math:`\alpha`-cascade from a **seed** at the level's starting direction
:math:`\mu_{\rm start}` (sphere: :math:`\mu = -1`).  The seed
:math:`\phi_{1/2,i} \equiv \psi_{1/2,i}` is the input field's value at
that closed ray.  Through Phase D–F (:ref:`sn-err-058-manifestation-b`)
the seed was produced by extrapolating the field **linearly in**
:math:`\mu` through the level's two most-inward ordinates —
operator-consistent for the forward *apply*, but structurally poisonous
for the *solve*.

The poison is a **directed-graph** fact.  A within-group solve
:math:`A\psi = q` is a single-pass direct inverse **iff** the operator is
triangular in the sweep-order permutation — every unknown depends only on
already-computed unknowns (:ref:`loss-rep-three-modes`; the #284
object-level discharge).  The two-point extrapolation seed reads
ordinate columns that the sweep visits **later** in the walk (the
inward-marching :math:`\mu<0` sweep needs the seed *before* it has swept
the :math:`\mu` neighbours the extrapolation samples).  In the sweep-order
permutation that places one entry **above the diagonal** — a back edge.
Consequences, all measured on a 4-cell homogeneous sphere:

* the spherical solve is **not** single-pass: the cold-start residual
  :math:`\lVert A\cdot\mathrm{solve}(b) - b\rVert_\infty/\lVert b\rVert_\infty`
  sits at :math:`5.18\times10^{5}` (a direct inverse must be
  :math:`\sim 10^{-16}`);
* the solve is **sensitive to the initial guess** — two random guesses
  give solutions differing by :math:`\Delta = 4.57\times10^{-2}` — because
  the seed reads the *previous iterate*, a lag that is harmless **at** the
  scattering fixed point but pollutes the cold, no-outer-iteration solve;
* a coarse pure-absorber (:math:`c = 0`, no scattering outer loop to mask
  the lag) returns ``NaN``; a coarse S\ :sub:`8` 16-cell fixed source
  returns ``NaN`` under source iteration and **negative flux** under
  Krylov.

The **#282 back edge is spherical-only**: the lagged two-point
extrapolation reads *later* ordinate columns only on the sphere's
Gauss–Legendre cascade (:math:`\tau_{{\rm raw},0}\in(0,1)`, the R12a
trichotomy below).  On a cylinder the starting direction carries **no
independent state** — a *dead* first-ordinate weight on a level-symmetric
rule (:math:`\tau_{{\rm raw},0}=1`, :math:`c_{\rm in}[m_0]=0`), a
:math:`\psi_0` rank-duplicate on a product rule
(:math:`\tau_{{\rm raw},0}=0`) — so no seed row lands **above** the
diagonal (the #282 "0.0-bit" row).  Route (a) therefore touches only the
sphere.

.. note::

   **Not** ':math:`\alpha`-dome telescoping' (#280 Phase 2.5b).  The
   cylinder's seed-insensitivity is the *dead first-ordinate weight* of
   the level-symmetric rule, a level-symmetric-only artefact (see
   :ref:`sn-phase-d-gate-1-1-empirical`).  On a **product** rule the seed
   :math:`\psi_0` is a **live self-coupling** on the :math:`m_0` diagonal
   (:math:`c_{\rm in}[m_0]\ne 0`), and the cold product-cylinder
   ``(L+C).solve`` was itself seed-**lagged** (cold error :math:`\approx
   0.57`) until the #280 Phase 2.5b direct-seed fold
   (:math:`c_{\rm out}\to c_{\rm out}-c_{\rm in}`) folded it onto the
   diagonal, making it a single-pass direct inverse — resolved by the
   SAME forward substitution the sphere route (a) certifies, not by any
   telescoping.

The pole is a straight characteristic — the physics beneath the direct solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The direct solve is not a numerical trick — it is a **physical property
of the closed rays**, and framing it that way (not as a storage choice)
is what makes the rest of route (a) inevitable.  In curvilinear geometry
the streaming operator carries an angular-redistribution term whose
strength is the factor :math:`(1-\mu^2)`:

.. math::

   \Omega\cdot\nabla\psi
   \;=\; \mu\,\frac{\partial\psi}{\partial r}
       \;+\; \frac{1-\mu^2}{r}\,\frac{\partial\psi}{\partial\mu}
   \qquad(\text{sphere}).

At the poles :math:`\mu = \pm 1` the coefficient :math:`(1-\mu^2)`
vanishes.  These are the **radial** directions: a particle at
:math:`\mu = \pm 1` streams straight through the origin, never changing
angular cell.  The redistribution term switches off — equivalently the
:math:`\alpha`-dome endpoints are :math:`\alpha_{1/2} = \alpha_{N+1/2} =
0` — and the transport equation at the pole collapses to a **pure 1-D
spatial ODE in radius alone**:

.. math::
   :label: sn-282-pole-straight-characteristic

   \mu\,\frac{d\psi_{1/2}}{dr} \;+\; \sigma_t(r)\,\psi_{1/2}(r)
   \;=\; \bar q_{1/2}(r),
   \qquad \mu = \mp 1,

with **no coupling to any other ordinate** (Hébert §3.9.4, the Carlson
inward march :eq:`hebert-3-434`–:eq:`hebert-3-435`).  The pole flux is
therefore computable *by itself*, before and independently of the angular
cascade it goes on to seed.

.. (vv-status rationale) Literature-transcribed derivation identity: the
.. sphere transport equation restricted to the straight-characteristic
.. pole μ=∓1, where (1-μ²)=0 kills angular redistribution.  Not a solver
.. claim — the verifiable content is the C(i) direct-solve residual
.. collapse and the block-triangular A_sb=0 certificate (§16.C), tabled
.. under the route-(a) evidence below.
.. vv-status: sn-282-pole-straight-characteristic documented

This is the physics that makes route (a) possible and the augmented
operator triangular.  The seed rows of :eq:`sn-282-block-triangular` are
self-contained (:math:`A_{\rm sb} = 0`) **because** the pole ODE reads no
bulk and no trace unknown — a straight characteristic couples to nothing
downstream.  The representation choice — promote :math:`\psi_{1/2}` to
first-class state — is *downstream* of the physics: any storage would
inherit the same decoupling.  What the lagged seed got wrong was never
the physics; it was reading the *iterate* for a quantity the pole ODE can
solve **directly** from the source (:ref:`sn-282-source-fold`).  The
deeper structural reason :math:`\mu = \pm 1` is the *only* admissible
starting direction in any curvilinear geometry is set out under
:ref:`sn-phase-d-pomraning-structural-singularity`.

ψ½ as first-class state — the augmented composite (System A ⊕ System B)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Route (a) kills the back edge by making the seed a **state variable the
solve computes**, not a functional of the iterate it reads.  The
within-group phase space carries a starting-direction subspace beyond the
bulk and trace:

.. math::
   :label: sn-282-augmented-composite

   V \;=\; V_{\rm bulk} \,\oplus\, V_{\rm trace} \,\oplus\, V_{\rm sd},
   \qquad
   V_{\rm sd} \;=\;
   \bigoplus_{p\,\in\,\mathcal{P}_{\rm carry}}
   \bigl(\underbrace{V_{1/2,p}^{-}}_{\mu=-1\ \rm leg}
         \oplus\, V_{1/2,p}^{+}\bigr),

one starting-direction block :math:`V_{1/2,p}^{\pm}` per **carrying**
:math:`\mu`-level :math:`p` (R12a; on production meshes only the sphere
carries — one level).  Each block holds the level's half-angle flux at
every radial cell (the **interior**
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicInteriorSpace`)
plus its two :math:`r = R` corner slots (the **boundary**
:class:`~orpheus.numerics.spaces.radial_characteristic_space.RadialCharacteristicBoundarySpace`)
— two flat backing buffers with typed ``(level, sign)`` views, mirroring the
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
trace layout.  The single unified ψ½ buffer (which then held the
``RadialCharacteristicSpace`` name) split into this interior ⊕ boundary
pair at 4e-e1.  Like every typed phase-space quantity in this codebase
the seed is realised by a **role family** — here a quadruple of roles,
each split into an ``interior ⊕ boundary`` locus pair (4e-e1) and composed
by System B's ``RadialCharacteristicField`` (role-erased slots — role
identity lives on the members):

.. list-table:: The ψ½ role family — four roles × two split loci, composed by ``RadialCharacteristicField`` (#282 route (a); split at 4e-e1)
   :header-rows: 1
   :widths: 34 30 36

   * - Role — the ``interior ⊕ boundary`` leaf pair
     - Realises
     - Forced by
   * - **flux** —
       :class:`~orpheus.transport.fields.radial_characteristic_interior_flux.RadialCharacteristicInteriorFlux`
       ⊕ :class:`~orpheus.transport.fields.radial_characteristic_boundary_flux.RadialCharacteristicBoundaryFlux`
     - the ψ½ state :math:`\psi_{1/2}` itself (the marched cells ⊕ the
       :math:`r = R` corner)
     - the carrier promotion (2.5d d1)
   * - **source/sink** —
       :class:`~orpheus.transport.source_sinks.radial_characteristic_interior_source_sink.RadialCharacteristicInteriorSourceSink`
       ⊕ :class:`~orpheus.transport.source_sinks.radial_characteristic_boundary_source_sink.RadialCharacteristicBoundarySourceSink`
     - the q½ source block :math:`\bar q_{1/2}` (the fold below) and any
       operator ``.apply`` output on the seed rows
     - the augmented source composite
   * - **displacement** —
       :class:`~orpheus.transport.displacements.radial_characteristic_interior_displacement.RadialCharacteristicInteriorDisplacement`
       ⊕ :class:`~orpheus.transport.displacements.radial_characteristic_boundary_displacement.RadialCharacteristicBoundaryDisplacement`
     - the affine displacement between two ψ½ states (minted per block by ⊖)
     - the composite torsor algebra (2.5d d1)
   * - **residual** —
       :class:`~orpheus.transport.residuals.radial_characteristic_interior_residual.RadialCharacteristicInteriorResidual`
       ⊕ :class:`~orpheus.transport.residuals.radial_characteristic_boundary_residual.RadialCharacteristicBoundaryResidual`
     - System B's typed residual members :math:`r_B = (A\psi)_B - q_B`
     - :func:`~orpheus.sn.solver.evaluate_residual`'s coupled arm (B.2d)

The historical **unified** single-buffer leaf — cells ⊕ corner interleaved
on one ``FaceField[(level, sign, part)]`` — carried the
``RadialCharacteristicField`` name until 4e; that name was reminted onto
System B's composite at 4e-e1b (see the "Where ψ½ lives" note below), and
the eight split role leaves above are the only ψ½ representation.

**Where ψ½ lives — System B, not a block on the bulk composite (the B.2d
eviction).**  The role family above is a **separate composite**, never a
third block bolted onto the bulk ⊕ trace field.  System A — the ordinary
bulk ⊕ trace SN state — is the pure **2-block**
:class:`~orpheus.transport.full_field.FullField`
(``Composite[BulkField, BoundaryField]``); the ψ½ ray is **System B**, its
own 2-block
:class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`
(the marched interior cells ⊕ the :math:`r = R` corner, carrying the same
flux / source-sink / displacement / residual role family), and on a
carrying sphere the driver iterate is the **coupled pair**
:class:`CoupledField[ψ_A, ψ_B] <orpheus.numerics.coupled_system.CoupledField>`
threaded through the within-group 2×2 loss grid (blocks :math:`A_{AA}`,
:math:`A_{AB}`, :math:`A_{BA}`, :math:`A_{BB}` — the :math:`A = M - N`
splitting built by
:func:`~orpheus.sn.coupled_system.build_within_group_system`; the grid
algebra is on :doc:`/theory/foundations/operator_algebra`).  The transitional 2.5d interim —
ψ½ as an **optional third block** on ``FullField`` with a mesh-keyed
*mixed-presence law* and runtime presence pins — is **retired**: a
live-ray :math:`\psi_A` is now **unrepresentable** (the type system is the
guard, not a runtime branch), so the B.2c dead-slot double-count hazard
(the welded seed feed **and** an explicit :math:`A_{AB}` block both firing
on one carrier) dissolved **structurally**.  The coupled flat dimension is
the **honest two-system sum** — no dead padding — which is why the ERR-053
``restart`` sizing (:ref:`sn-282-gotchas`) reads the true count off the
coupled ravel.  The converged ψ½ state is returned as
:attr:`Solution.radial_characteristic <orpheus.sn.solution.Solution.radial_characteristic>`
— System B's **own typed member**, ``None`` exactly when the mesh carries
no seed level (presence validated as a **biconditional** at construction)
— while :attr:`Solution.angular_flux <orpheus.sn.solution.Solution.angular_flux>`
stays the honest 2-block System-A composite.

**The state metric.**  The inner-product weight of :math:`V_{\rm sd}`
is the SPD **state metric** :math:`G_{\rm sd} = V_{\rm cell}` — the
radial cell-volume measure, mirroring the bulk metric
:math:`G_{\rm bulk} = V_{\rm cell}\,w_n` (the SAME spatial measure,
restricted to the single :math:`\mu = \pm 1` ray, without the angular
factor :math:`w_n`).  ψ½ is a **first-class radial state field**, not a
face trace: its operator self-block :math:`A_{\rm ss}` is a *banded
radial transport operator* :math:`\mu\,\partial_r + \sigma_t` (Hébert
Eqs. 3.434–3.435), so — like any state — its Hilbert metric is set by its
**operator role**, not by an integration weight.

Three pole-vanishing quantities were historically conflated into one
"ghost" zero; keep them apart — only the operator coefficient is zero:

.. list-table:: The three pole-vanishing quantities at :math:`\mu = \pm 1`
   :header-rows: 1
   :widths: 8 34 58

   * - Tag
     - Quantity
     - Where it lives / what it governs
   * - **M1**
     - moment / output weight :math:`= 0`
     - the *open* Gauss–Legendre rule has **no node** at the pole, so ψ½
       carries zero weight in :math:`\phi = \sum_n w_n\psi_n` — it lives
       in the **moment reducer** and correctly excludes ψ½ from the
       scalar flux.
   * - **M2**
     - angular through-flux coefficient
       :math:`(1-\mu^2)\big|_{\mu=\pm 1} = 0` (the :math:`\alpha`-dome
       endpoints :math:`\alpha_{1/2} = \alpha_{N+1/2} = 0`)
     - an **operator coefficient inside** :math:`A` — the
       angular-redistribution strength that makes the pole a straight
       characteristic (:eq:`sn-282-pole-straight-characteristic`).
       Correctly zero.
   * - **M3**
     - **state metric** :math:`G_{\rm sd} = V_{\rm cell} \neq 0`
     - *this block's* inner product — governs the G-adjoint reciprocity
       :math:`\langle A\psi,\chi\rangle_G =
       \langle\psi, A^{\dagger}\chi\rangle_G`.

The retired **"ghost metric" bug** installed **M2** (an operator
coefficient) as **M3** (the state metric): it read the angular
through-flux weight :math:`(1-\mu^2)|_{\rm pole} = 0` as the Hilbert
metric and set :math:`G_{\rm sd} \equiv 0`.  Because
:meth:`~orpheus.numerics.operator.SupportsAdjoint.apply_transpose` is the
*exact* Euclidean transpose (:math:`\lVert T - A^{\mathsf T}\rVert =
3.6\times10^{-16}`), the relation :math:`A^{\mathsf T}G = G A^{\dagger}`
behind :math:`A^{\dagger} = G^{-1}A^{\mathsf T}G` holds for **every** SPD
:math:`G_{\rm sd}` (the reciprocity is gauge-free among SPD choices), and
:math:`0` is the **one forbidden value** — it puts the seed rows in
:math:`\ker G`, severing the seed :math:`\to` bulk coupling
:math:`A_{\rm bs}` from :math:`A^{\dagger}` (a wrong adjoint the instant
the seed carries data — a production reciprocity defect of
:math:`1.3\times10^{-2}`, green only on a present-but-zero seed).  We
gauge-fix to :math:`V_{\rm cell}` (dropping the angular :math:`w` — a
single :math:`\mu=\pm 1` ray has no canonical quadrature weight) so the
adjoint's seed block is the physical **backward radial march** and all
bulk/trace observables (:math:`\phi^{\dagger}`, adjoint reaction rates)
are bitwise **gauge-invariant** (the block-upper-triangular
:math:`A^{\dagger}` seats the seed at the top, so only
:math:`\phi^{\dagger}_{\rm seed}` moves with the gauge).  Consequently
:meth:`~orpheus.numerics.space.FunctionSpace.apply_metric` **scales** the
block by :math:`V_{\rm cell}`, its inverse **divides** (empty null
space), and the block contributes :math:`\sum V_{\rm cell}\,x\,y` to the
composite inner product.  This closes a sharp V&V gap — see
:ref:`sn-282-gotchas` (Mode 12, ERR-067).  The gauge derivation of record
is
``.claude/agent-memory/numerics-investigator/radial_characteristic_metric_gauge_derivation.md``.

.. (vv-status rationale) Structural / representational identity: the
.. named-field-typed decomposition of the augmented within-group phase
.. space.  Not a solver claim — the verifiable content is the
.. RadialCharacteristicSpace layout / role-quadruple class-identity
.. arithmetic (Field Layer-1 gate) + the §16.A carrier gates.
.. vv-status: sn-282-augmented-composite documented

.. _sn-282-pole-state-metric:

The through-flux coefficient is not the state metric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The M1/M2/M3 table above turns on one structural fact worth spelling out,
because getting it wrong is exactly what produced the retired ghost
metric.  A codim-1 face carries a **through-flux** coefficient (M2) — the
normal component of the streaming flux across the face, which vanishes
when the characteristic runs *tangent* to it.  Both codim-1 traces have
one, and for one of them the through-flux coincides with the state
metric while for the other it does not:

.. list-table:: The through-flux coefficient (M2) is an OPERATOR coefficient on both faces
   :header-rows: 1
   :widths: 26 22 28 24

   * - Codim-1 face
     - Through-flux M2
     - Vanishes when…
     - State metric M3
   * - **spatial** :math:`r`-face — the boundary trace
       (:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`)
     - :math:`\lvert\Omega\cdot n\rvert\,w`
     - :math:`\Omega` is **tangent** to the surface
       (:math:`\lvert\Omega\cdot n\rvert = 0`, grazing incidence)
     - **equals** the through-flux
       :math:`\lvert\Omega\cdot n\rvert\,w`
   * - **angular** :math:`\mu = \pm 1` edge — the ψ½ seed
       (System B's
       :class:`~orpheus.transport.radial_characteristic_field.RadialCharacteristicField`)
     - :math:`(1-\mu^2)\,w = 0`
     - the ray is a **straight characteristic**
       (:math:`(1-\mu^2) = 0`, the pole)
     - :math:`V_{\rm cell}` (radial cell volume) — **NOT** the
       through-flux

For the pole, :math:`(1-\mu^2)\,w \equiv 0` correctly captures the M2
through-flux: the :math:`\mu = \pm 1` angular face is *entirely grazing*,
so no flux streams *across* it.  That is the straight-characteristic
physics of :eq:`sn-282-pole-straight-characteristic`, and it is exactly
why the augmented operator is triangular.  **What is wrong is reading that
through-flux coefficient as the block's Hilbert STATE metric.**

The through-flux coefficient equals the state metric only when the face's
**operator self-block is trivial**.  The spatial trace's self-block
:math:`A_{\rm tt}` is a pure restriction / reflection map (measured
off-diagonal norm :math:`\approx 2` on a 6-cell sphere, diagonal in
:math:`[-1, 1]` — no interior dynamics), so there the through-flux, the
state metric, and the partial current all coincide.  The pole's
self-block :math:`A_{\rm ss}` is a **banded radial transport operator**
:math:`\mu\,\partial_r + \sigma_t` (off-diagonal norm :math:`\approx 71`,
about :math:`35\times` the trace's, diagonal in :math:`[-1, 3.65]`) with
genuine interior radial dynamics — so its through-flux (:math:`0`) and its
state metric (:math:`V_{\rm cell}`) are *different objects*.  Installing
the M2 through-flux as the M3 state metric is a **category error** — an
operator coefficient placed where the inner product belongs — and that is
precisely the retired ghost :math:`G_{\rm sd} \equiv 0`.  The state metric
is the radial cell volume :math:`V_{\rm cell}` (derived, gauged, and
V&V-consequential above).

.. important::

   **Three measures at the pole; do not conflate them.**  M1, M2, M3 (the
   first table) answer three different questions, and only the operator
   coefficient M2 is zero:

   * **M1 — the scalar-flux moment** :math:`\int\psi\,d\mu`: "how much
     does this ray contribute to :math:`\phi`?"  *Rule-dependent* — under
     the sphere's *open* Gauss–Legendre rule
     (:ref:`sn-282-circle-vs-interval`) the pole has **no interior node**,
     so its moment weight is zero and the seed is a pure auxiliary DOF;
     under a pole-*including* rule (Gauss–Lobatto,
     :ref:`sn-282-lobatto-study`) it would carry a genuine nonzero moment
     weight.
   * **M2 — the angular through-flux** :math:`(1-\mu^2)`: "how much
     streams *across* this angular face?"  Zero at the pole **always**
     (independent of the quadrature) — an operator coefficient, never a
     state metric.
   * **M3 — the state metric** :math:`G_{\rm sd} = V_{\rm cell}`: "what
     is the inner product on the ψ½ state?"  Nonzero **always** — set by
     the operator role.

   The retired bug conflated **M2 with M3** (an operator coefficient read
   as the Hilbert metric).  The moment-vs-through-flux distinction (**M1
   vs M2**) is a *separate* trap — both are quadrature/operator facts, and
   neither is the state metric.

The walk triple — solve marches, apply reads, transpose reverses
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The augmented loss operator :math:`A = L + C` acts on the seed block
through three orientation-coherent paths of the **one** 1-D scan/loop
walk (:class:`~orpheus.sn.loss_representation._OneDimScanWalk`; the #280
unification).  Ordering the unknowns **seed⁻ ≺ seed⁺ ≺ ordinate legs**
makes the augmented operator block-lower-triangular:

.. math::
   :label: sn-282-block-triangular

   A \;=\;
   \begin{bmatrix} A_{\rm ss} & 0 \\[2pt] A_{\rm bs} & A_{\rm bb} \end{bmatrix},
   \qquad
   \begin{aligned}
   A_{\rm ss} &: \text{the seed rows (Hébert DD residual + corner rows),} \\
   A_{\rm bs} &: \text{the seed} \to \text{bulk M-M recurrence coupling,} \\
   A_{\rm bb} &: \text{the bulk} (L+C)\ \text{walk (}A_{\rm sb} = 0\text{).}
   \end{aligned}

The zero upper-right block :math:`A_{\rm sb} = 0` **is** the death of the
back edge: the seed rows are *self-contained* in the seed state (they read
no bulk or trace unknown), and the seed → bulk coupling is one-directional
(the M-M recurrence reads :math:`\psi_{1/2}` into the bulk rows'
``angular_numer``, never the reverse).  The three orientations:

* **SOLVE** — routed through System B's named resolvent
  :meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
  (the **4e-e2 un-weave**: the walk constructs ``A_BB`` over its own
  :math:`\sigma_t` and calls it **once, up front** — System B is
  iterate-independent — instead of inlining the two-leg march).  The two
  ψ½ legs are solved **directly** from the true q½ source: march inward
  from the seeded :math:`r = R` inflow corner (vacuum ⇒ 0; reflective ⇒
  the mirror outflow corner) via the Hébert
  :eq:`hebert-3-434`–:eq:`hebert-3-435` DD recurrence
  (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`,
  the single-sourced DD **engine** — the only place the march lives),
  pole-continue :math:`\psi^{+}_{1/2}(0) = \psi^{-}_{1/2}(0)` (the inward
  march's exit face **is** the pole datum), then march the :math:`\mu = +1`
  leg outward to the :math:`r = R` outflow corner on the **reversed** cell
  data (orientation is carried by the data, never a flag — the 2.5a
  discipline).  The walk then reads the marched inward cells off the
  carrier **as** the M-M recurrence seed; the iterate plays no role.
* **APPLY** (the matvec).  Reads the *given* ψ½ carrier and **emits** the
  seed-block rows: per leg the Hébert (3.434) residual
  :math:`m_{1/2,i} = \sigma_i\psi_{1/2,i} + (2/\Delta r_i)(\psi_{1/2,i} - f_{{\rm in},i})`
  reconstructed from the stored cells (the DD face chain replays the
  solve's arithmetic), plus the inflow-corner *identity* row and the
  outflow-corner *streamed − stored* defect row.  Because the apply
  replays the solve's ops in the same order, ``apply ∘ solve`` closes the
  corner defect to :math:`0` bit.
* **TRANSPOSE** (``loss_action_transpose``).  The exact reverse-mode of
  the forward straight-line program: reverse the corner rows, the +
  leg's chain (descending), the pole continuation, the − leg's chain
  (ascending); the reverse M-M recurrence
  (:meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.angular_adjoint`)
  **stops** at the seed cotangent and lands it on the output composite's
  seed block, rather than scattering it back onto the bulk.

The triangularity is not asserted — it is **certified** by a probed
``triu == 0`` check on the assembled augmented block (the transpose
analogue of the #284 object-level discharge): the sweep is LAPACK forward
substitution on the source subspace.

.. note:: **The ray solve is un-woven from the walk (Cardinal Rule 2 — 4e-e2).**

   The two ψ½ **orchestrations** — the forward SOLVE and its reverse-scan
   solve-transpose — no longer live inline in the walk.  Since 4e-e2 both
   route through System B's named resolvent
   :meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
   / :meth:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve_transpose`
   (``A_BB``, constructed inside the walk over its own :math:`\sigma_t`),
   which the coupled 2×2 grid also exposes at its ``(B, B)`` slot.  The
   two-leg **orchestration** (read source views → inward leg →
   pole-continue → reversed outward leg → write flux views) now lives in
   exactly **one** place —
   :class:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator`
   — and the DD **engine** it drives in exactly one other
   (:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
   / :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_transpose`).
   The walk's ``carlson_inward_sweep_*`` references went **8 → 0**; a
   source-scan tripwire (``test_4e_unweave_walk_source_has_no_carlson_reference``)
   pins that the walk holds zero references to the engines, and the
   :ref:`Mode-11 wrap-sentinels <sn-282-numerical-evidence>` re-aim onto
   the operator's ``solve`` / ``solve_transpose``.  The **forward matvec**
   and its transpose were single-sourced earlier (step 4b) onto
   :func:`~orpheus.sn.sweep.psi_half_angle_seed.radial_characteristic_forward_residual`,
   so ``apply`` and the walk's seed rows already share one kernel.

   **H1-narrow.**  Only the ``A_BB`` *solve* legs were extracted; the
   ``A_AB`` seed → bulk coupling stays **fused** in the within-group
   :math:`M` (the DP-splitting ruling).  In the transpose the reversed
   ordinate loop accumulates the Morel–Montry thread cotangent onto a copy
   of the seed cotangent — *that augmentation is* the fused
   :math:`A_{AB}^{\mathsf T}` feed — and only then does one
   ``A_BB.solve_transpose`` call return the seed-source cotangent.

.. (vv-status rationale) Structural / representational identity: the
.. block-triangular normal form of the augmented (L+C) in the augmented
.. walk order.  The verifiable content is the triu==0 triangularity
.. certificate + the apply∘solve corner-defect=0-bit gate (§16.C), not a
.. flux/eigenvalue claim.
.. vv-status: sn-282-block-triangular documented

.. _sn-282-source-fold:

The starting-direction source fold — why ALL Legendre moments (R14)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The direct solve needs the true within-group source **at the starting
direction**, :math:`\bar q_{1/2}(\mu = \pm 1)` — the value the anisotropic
source :math:`q(r,\mu)` takes at the closed ray, reconstructed from **all**
its Legendre moments (Hébert Eq. (3.432); the "R14 full fold",
:func:`~orpheus.numerics.spaces.radial_characteristic_space.fold_moments_to_radial_characteristic`):

.. math::

   \bar q_{1/2}(\mu = \pm 1)
   \;=\; \sum_{\ell} \frac{2\ell+1}{2}\, q_\ell\,(\pm 1)^\ell,
   \qquad
   q_\ell(r) \;=\; \sum_n w_n\,P_\ell(\mu_n)\,q_n(r),

the exact 1-D addition-theorem weight :math:`(2\ell+1)/2` with
:math:`P_\ell(\pm 1) = (\pm 1)^\ell` — this is Eq. :eq:`hebert-3-432-source`
kept at **full order**, not collapsed to :math:`L = 0`.  The single q½
fold factory
:meth:`RadialCharacteristicField.source_from_angular <orpheus.transport.radial_characteristic_field.RadialCharacteristicField.source_from_angular>`
folds every moment the level resolves; the solver cold-start, the
fixed-source right-hand side, and the within-group joint march (the
System-B q½ source fed to the resolvent grid) all route through
it (Pattern 2 — the :math:`P_1(-1)` sign is spelled **once**).

**The full fold is load-bearing, and this is the subtle part.**  It is
tempting to fold only :math:`\ell = 0` (the apply matvec's isotropic
scattering reach is P0).  That is *wrong for the source*, because
**streaming manufactures angular structure the flux does not have**.
Take an isotropic trial flux :math:`\psi = A(r)` (no :math:`\mu`
dependence).  Applying the spherical streaming–collision operator
(:math:`\Omega\cdot\nabla + \sigma_t`, with
:math:`\Omega\cdot\nabla\psi = \mu\,\partial_r\psi + \tfrac{1-\mu^2}{r}\partial_\mu\psi`
and :math:`\partial_\mu\psi = 0`) gives a source that is **linear in**
:math:`\mu`:

.. math::
   :label: sn-282-anisotropic-source

   q(r,\mu) \;=\; \mu\,A'(r) \;+\; \sigma_t(r)\,A(r),
   \qquad
   q(r,\mu = -1) \;=\; \sigma_t A - A'.

.. (vv-status rationale) Literature-grounded derivation identity: the
.. spherical streaming-collision operator applied to an isotropic trial
.. flux A(r).  Not a solver claim — it is the analytical basis for the
.. full-fold requirement, whose verifiable content is the fold's ℓ=0
.. collapse to ½q₀ bit-identity (isotropic paths unchanged) + the
.. anisotropic-MMS O(h²) convergence gate.
.. vv-status: sn-282-anisotropic-source documented

The value at :math:`\mu = -1` carries the :math:`-A'` term — and that
term lives **entirely in the** :math:`\ell = 1` **moment**.  Working the
fold explicitly on a Gauss–Legendre level
(:math:`\sum_n w_n\mu_n = 0`, :math:`\sum_n w_n = 2`,
:math:`\sum_n w_n\mu_n^2 = 2/3`):

.. math::

   q_0 = 2\sigma_t A, \quad q_1 = \tfrac{2}{3}A',
   \qquad
   \underbrace{\tfrac12 q_0 P_0(-1)}_{\sigma_t A}
   + \underbrace{\tfrac32 q_1 P_1(-1)}_{-A'}
   \;=\; \sigma_t A - A' \;\checkmark

An :math:`\ell = 0`-only fold keeps only the :math:`\tfrac12 q_0 = \sigma_t A`
piece and **drops** the :math:`-A'`.  That dropped term is exactly what
**floored the anisotropic curvilinear MMS**: with the ℓ=0-only seed the
manufactured-solution error refused to converge; with the full fold the
sphere MMS converges :math:`\mathcal{O}(h^2)`-to-exact.

**For an isotropic source the full fold is a no-op.**  The higher moments
vanish (:math:`q_\ell = 0` for :math:`\ell \ge 1`), the fold collapses to
:math:`\tfrac12 q_0` **bit-exactly**, and the eigenvalue and
isotropic-fixed-source paths are **unchanged** by route (a).  The full
fold only "activates" when a >linear-in-:math:`\mu` source is present
(the companion anisotropic-MMS gate); in current production runs — all
isotropic sources — :math:`\ell \ge 1` is manufactured-before-needed and
identically zero (the isotropic-snapshot-blindness discipline: the
machinery is correct for the case the snapshots cannot exercise).

.. _sn-282-r12a:

Which levels carry a ψ½ block — the R12a predicate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A :math:`\mu`-level carries an **independent** starting-direction state
block **iff** the M-M half-angle recurrence genuinely *consumes* it — a
structural predicate on the level's first-ordinate **raw** (unclamped)
Morel–Montry weight
(:func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_raw_per_level`,
read by :attr:`~orpheus.sn.mesh.augmented_mesh.SNMesh.radial_characteristic_levels`):

.. math::

   \text{level } p \text{ carries a ψ½ block}
   \quad\Longleftrightarrow\quad
   \tau_{{\rm raw},0}^{(p)} \,\in\, (0, 1)\ \text{exclusive.}

The trichotomy is **bit-exact** on the production quadratures:

.. list-table:: The R12a carrying-level trichotomy
   :header-rows: 1
   :widths: 18 30 52

   * - :math:`\tau_{{\rm raw},0}`
     - Rule
     - Why the seed is (not) independent state
   * - :math:`= 0`
     - cylinder **product** rules
     - the starting direction coincides with the first ordinate
       (:math:`\eta_0 = \eta_{1/2} = -\sin\theta` bit-exactly, the #229
       clamp fact) — the seed is a rank-duplicate of :math:`\psi_0`.
       **No block.**
   * - :math:`= 1`
     - cylinder **level-symmetric** rules
     - duplicate-:math:`\eta` nodes collapse the midpoint edge onto
       :math:`\eta_0`, so the seed's only consumption path — the
       recurrence weight :math:`(1-\tau_0)` — vanishes.  Dead state.
       **No block.** (This is why the measured cylinder-LS seed
       sensitivity is :math:`0.0` bit.)
   * - :math:`\in (0,1)`
     - sphere **Gauss–Legendre**
     - the recurrence consumes the seed with a genuine weight
       (:math:`\tau_{{\rm raw},0} \approx 0.39\text{–}0.42`).  **Carries.**

So on production meshes **only the sphere carries the block**; every
cylinder inlines the 2-point angular-edge extrapolation
(:meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`
— harmless and bit-exact by the trichotomy above).  R12a **refines** the
earlier R12 letter ("μ_start ∉ the level's μ-nodes"), whose claimed
equivalence to :math:`\tau_{\rm raw} \ne 0` is empirically **false** on
level-symmetric cylinder rules (μ_start ∉ nodes there, yet
:math:`\tau_{\rm raw} = 1` — dead).  The clamp
:math:`0 \mapsto \tfrac12` erases exactly the 0-vs-(0,1) distinction the
predicate needs, which is why the **raw** producer is first-class.

.. _sn-282-circle-vs-interval:

Why the sphere pays for the pole and the cylinder does not — circle vs interval
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The R12a trichotomy raises a deeper question: *why* does the sphere carry
an independent seed while the cylinder never does?  The answer is neither
"sphere vs cylinder" nor "curvilinear vs Cartesian" — it is the
**topology of the redistribution axis**, and it is the single most
clarifying fact about the whole ψ½ apparatus.

Two orthogonal questions must be kept apart:

#. *Does the geometry have an angular-redistribution term*
   :math:`\tfrac{1-\mu^2}{r}\,\partial_\mu`? — this is
   **curvilinear vs Cartesian** (sphere **and** cylinder: yes;
   Cartesian: no).
#. *Does the angular sweep need a separate, off-node starting DOF?* — this
   is the τ_raw predicate (:ref:`sn-282-r12a`), which is
   **quadrature-structural, not geometric**.

ψ½ answers question 2, and the deciding fact is what the redistribution
axis *is*.

**Cylinder — the redistribution axis is a circle.**  At a fixed polar
cosine :math:`\mu_z = \cos\theta`, the cylinder redistributes across the
**azimuthal angle** :math:`\varphi`, which lives on a **circle**
:math:`[0, 2\pi)` — a *periodic* domain.  The production rule
(:func:`~orpheus.numerics.quadrature.rules_product.product_mu_phi`) is
Gauss–Legendre in :math:`\mu_z` **×** *equispaced* in :math:`\varphi`
(``np.linspace(0, 2π, n_φ, endpoint=False)``).  Equispaced sampling of a
smooth periodic function is the trapezoidal rule, which on a circle is
**spectrally accurate** (its error decays faster than any power of
:math:`1/n_\varphi`) — there is no accuracy penalty for the choice.  And
crucially, for **even** :math:`n_\varphi` the grid hits
:math:`\varphi = \pi` exactly (partition node :math:`k = n_\varphi/2`),
where

.. math::

   \mu_x = \sin\theta\cos\pi = -\sin\theta, \qquad
   \mu_y = \sin\theta\sin\pi = 0,

i.e. the **most-inward radial direction** of that level.  The
starting-edge ordinate :math:`\eta_0 = -\sin\theta` therefore lands
*exactly on a quadrature node* — :math:`\tau_{{\rm raw},0} = 0` — and the
seed is a **bulk ordinate for free, at no accuracy cost**.  This is the
structural content of #229: the cylinder's edge-inclusion is a property
of the *circle*, **not** the :math:`[\tfrac12, 1]` Morel–Montry clamp
(the clamp is a separate recurrence-weight stabiliser; the R12a predicate
reads the *un*clamped :math:`\tau_{\rm raw}`).  It is contingent on
**even** :math:`n_\varphi` — an odd azimuthal count would miss
:math:`\varphi = \pi`, and the cylinder *would* then carry a seed.

**Sphere — the redistribution axis is an interval.**  The sphere
redistributes across the **polar cosine** :math:`\mu \in [-1, 1]`, an
**interval** whose two endpoints :math:`\mu = \pm 1` *are* the physical
poles.  The optimal rule on an interval with a smooth integrand is
Gauss–Legendre — but Gauss–Legendre is an **open** rule: it places *no*
node at the endpoints.  So the sphere structurally *cannot* put an
ordinate on :math:`\mu = \pm 1`; the starting edge falls strictly between
nodes (:math:`\tau_{{\rm raw},0} \approx 0.39\text{–}0.42`), and a
**separate off-node seed DOF is unavoidable**.

.. list-table:: The redistribution axis decides the seed
   :header-rows: 1
   :widths: 16 26 24 17 17

   * - Geometry
     - Redistribution axis
     - Optimal rule
     - Edge-inclusive?
     - Seed?
   * - **Cylinder**
     - azimuth :math:`\varphi` — a **circle** (periodic)
     - equispaced (trapezoidal, spectral)
     - **yes** (even :math:`n_\varphi` hits :math:`\varphi=\pi`)
     - no — :math:`\tau_{\rm raw}=0`
   * - **Sphere**
     - polar :math:`\mu` — an **interval** :math:`[-1,1]`
     - Gauss–Legendre (open)
     - **no** (open rule, no endpoint node)
     - yes — :math:`\tau_{\rm raw}\in(0,1)`

The principle in one line: **a periodic redistribution axis gives
edge-inclusion for free; an interval axis makes you pay for it with a
separate seed.**  The cylinder is the existence proof that "seed = a bulk
edge-ordinate" *works* — it works there precisely because the axis is a
circle.  The sphere pays because its axis has physical endpoints and the
best interior rule refuses to stand on them.

.. _sn-282-lobatto-study:

Could the sphere put a node at the pole? — the Gauss–Lobatto study
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The circle-vs-interval framing suggests an obvious question: could the
sphere *buy* the cylinder's free edge-inclusion by switching from
Gauss–Legendre to **Gauss–Lobatto** — a rule that *does* place nodes at
the interval endpoints :math:`\mu = \pm 1` (at the cost of exactness
:math:`2n-3` versus GL's :math:`2n-1`)?  A pole node would give
:math:`\tau_{{\rm raw},0} = 0`, making the seed a bulk ordinate exactly as
on the cylinder — dissolving the whole ψ½ block.  A dedicated empirical
study (scratch, uncommitted — see the note below) answered both halves of
the question.

**Affordable.**  At resolved angular order (:math:`N \ge 8`, and
:math:`N > L` for a :math:`P_L` scattering source) Gauss–Lobatto tracks
Gauss–Legendre at a bounded :math:`\sim 1.2\times` error penalty —
:math:`\sim 1.3\text{–}1.4\times` at S\ :sub:`8`, tightening to
:math:`\sim 1.2\times` at S\ :sub:`16` — i.e. **one to two extra
ordinates** to match GL.  The penalty is **not amplified by anisotropy**
(P0 through P5 all sit at :math:`\sim 1.2\text{–}1.3\times`) and is
**insensitive to the scattering ratio** :math:`c`.  The eigenvalue offset
is :math:`\sim 30\text{–}140` pcm at S\ :sub:`16`, and fine-:math:`N` GL
and GLob agree to :math:`< 6` pcm — the two rules converge to the **same**
:math:`N \to \infty` transport limit (the pole weighting is unbiased; the
straight-characteristic pole handling is redistribution-consistent, with a
per-ordinate flat-flux residual :math:`\sim 10^{-15}`).  Only the
under-resolved :math:`N \lesssim L` corner breaks, and there GL is
rank-deficient too.

**But not a drop-in.**  A pole node lands on the level's lower edge, so
the first-ordinate weight :math:`\tau_{{\rm raw},0} = 0` — and the
production Morel–Montry recurrence :eq:`pole-mm-recurrence` **divides its
first step by that weight**, so the recurrence is *singular*; separately,
the R12a presence predicate keys on :math:`\tau_{\rm raw} \in (0,1)`,
which a pole node also fails.  Adopting a pole-node quadrature is
therefore **not** a quadrature swap — it *requires* the same
seed→bulk-ordinate restructure the cylinder already uses (make the pole
node the seed, straight-characteristic-solved, and start the recurrence
*from* it).

**Ruling: affordable but architecturally declined.**  The fold-in was
**not** adopted, and the reason is architectural, not numerical.  The bulk
is **cell-centred**; a pole ordinate would make it a *mixed* field — an
inert, zero-through-flux, straight-characteristic-solved,
redistribution-special passenger that *every* bulk consumer
(homogenization, condensation, moment extraction, every
``for ordinate in bulk`` loop) would have to know about and skip.  That is
the Cardinal-Rule-2 smell of two concepts in one type forcing a demux
downstream; the zero weight prevents *numerical* corruption but not the
*conceptual* pollution.  The value of the study is precisely that it makes
keeping ψ½ a **separate** object a *chosen* architecture — the pole seed
is kept out of the bulk **because the bulk stays clean**, not because a
pole-node scheme is infeasible.  The clean-bulk / ``FaceField``
architecture this decides is set out on the loss-operator page
(:ref:`loss-rep-facefield-codim1`).

.. note::

   The Gauss–Lobatto study is a set of scratch diagnostics
   (``scratch/experimental/glob_sphere_study/`` and
   ``derivations/diagnostics/diag_glob_0{1..5}_*.py`` — 33 green
   diagnostics covering moment integration, per-ordinate consistency,
   end-to-end penalty, the :math:`\tau_0 = 0` recurrence break, and the
   :math:`k_\infty` anchor).  They are **uncommitted** and are promotion
   targets *only if* a pole-node scheme is ever adopted; do not promote
   them otherwise.  The durable synthesis is
   ``.claude/plans/facefield_codim1_design.md`` §3.5.

.. _sn-282-seed-strategy-zoo:

What was tried and failed — the seed-strategy zoo (retired)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Three swappable ``PsiHalfAngleSeed`` strategies preceded route (a).  All
three are **retired** (2026-07-04); the module
:mod:`~orpheus.sn.sweep.psi_half_angle_seed` shrank from 851 lines to
161 — one engine.  The history matters because it prevents a future
session from re-attempting a known-dead seed:

.. list-table:: The retired seed strategies and why each fell short
   :header-rows: 1
   :widths: 26 40 34

   * - Retired strategy (literal — class deleted)
     - What it did
     - Why it failed
   * - ``ZeroSeed``
     - hardcoded :math:`\psi_{1/2} = 0`
     - the pre-ERR-026 term-initialisation bug (vv-principles failure
       Mode 3, a missing term); wrong off flat flux.
   * - ``CarlsonInwardSweep``
     - this module's Hébert march driven by the **proxy** source
       :math:`\bar Q = \sigma_t\phi_0/\!\sum w`
     - exact only at the flat-flux equilibrium; :math:`\mathcal{O}(1)`
       wrong off equilibrium — **floored the curvilinear MMS at ≈ 0.04**
       L2 independent of mesh (ERR-058b, Issue #195).
   * - ``AngularEdgeExtrapolation``
     - the operator-consistent 2-point angular extrapolation of the
       **iterate**
     - fixed the forward *apply* (#195) but left the *solve* seeded from
       the **previous iterate** — the #282 walk-order **back edge**
       (sphere cold residual :math:`5.18\times10^5`, seed sensitivity
       :math:`4.57\times10^{-2}`).

Route (a) retired the whole family.  The **arithmetic survives** in two
places, both correct on their own:

* :func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
  — the Hébert :eq:`hebert-3-434`–:eq:`hebert-3-435` recurrence, now
  driven by the **true** q½ source (:ref:`sn-282-source-fold`) instead of
  the falsified proxy, and used as the SOLVE engine (not a strategy);
* :meth:`~orpheus.sn.sweep.pole_angular_closure.MorelMontryAngularSweep.edge_extrapolated_seed`
  — the 2-point angular-edge extrapolation, inlined **verbatim** for the
  non-carrying cylinder levels (where the R12a trichotomy makes it
  bit-identical to the retired default: :math:`t = 0` exact on product
  rules, dead seed weight on level-symmetric rules).

.. _sn-282-numerical-evidence:

Numerical evidence — the lag death
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The acceptance gates live in
:file:`tests/sn/sweep/curvilinear/test_282_direct_seed_fixed_point.py`
(the §16.C fixed-point classifiers); every gate measures the **full
coupled state** (System A's bulk ⊕ trace *and* System B's ψ½), because a
bulk-only norm would be blind to any seed error (a Mode-12
functional-invariance point, independent of the metric closure in the
gotcha below).

.. list-table:: #282 route-(a) acceptance evidence
   :header-rows: 1
   :widths: 16 44 40

   * - Gate
     - Measurement
     - Before → after (route a)
   * - **C(i)** cold residual
     - :math:`\lVert A\cdot\mathrm{solve}(b)-b\rVert_\infty/\lVert b\rVert_\infty`
       on a cold start (the keystone: the solve is a single-pass exact
       inverse)
     - sphere :math:`5.18\times10^{5} \to 2.5\times10^{-16}`; slab & cyl
       already :math:`< 10^{-11}` (must **stay**)
   * - **C(ii)** seed-insensitivity
     - two random ``initial_guess`` seeds → bitwise-identical sphere
       solve (the lag signature)
     - :math:`\Delta = 4.57\times10^{-2} \to 0` **bitwise**
   * - **C(ii)** Probe-6
     - :math:`\psi_0` arbitrary, :math:`b = A\psi_0`, cold
       :math:`\mathrm{solve}(b)` recovers :math:`\psi_0`
     - pre-fix only the **warm** sphere solve recovered it; now the cold
       solve does (``rtol`` :math:`10^{-11}`)
   * - **C(iii)** coarse physicality
     - S\ :sub:`8` 16-cell sphere fixed source, finite **and**
       non-negative on **both** inner drivers
     - SI → ``NaN`` / Krylov → negative flux :math:`\to` finite + positive
       on both
   * - **C(iv)** pure absorber
     - :math:`c = 0` sphere (no scattering outer loop to mask the lag) —
       the cold solve **is** the answer
     - ``NaN`` :math:`\to` single-pass exact inverse
       (:math:`< 10^{-11}`, finite, positive)

Three teeth pin that the evidence is not vacuous: **Mode 11** —
a **class-level** wrap-sentinel confirms the sphere cold solve *executes*
:meth:`RadialCharacteristicOperator.solve <orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve>`
while the slab does **not** (no carrying levels), the transpose analogue
wraps :meth:`~orpheus.sn.operators.radial_characteristic.RadialCharacteristicOperator.solve_transpose`,
and a source-scan tripwire
(``test_4e_unweave_walk_source_has_no_carlson_reference``) pins that the
walk holds **zero** ``carlson_inward_sweep_*`` references — the 4e-e2
un-weave (the marches live only behind the operator, so the sentinel wraps
the very method the walk routes through); **Mode 10** — zeroing
the q½ source block **moves** the sphere solve (the carrier is not inert);
**Mode 12** — with the :math:`V_{\rm cell}` state metric the G-reciprocity
gate now *catches* a seed-row flip (the closure, ERR-067; gotcha below).

The eigenvalue re-pose — an N-sweep, not h→0
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Route (a) changed the ψ½ **angular closure**, so it re-posed the sphere
:math:`k`-eigenvalue — by :math:`\sim 4.66\times10^{-4}` at :math:`n = 40`
cells.  Judging whether such a re-pose is **principled** requires the
right discriminator, and this is a genuinely teachable subtlety.

**The** :math:`h \to 0` **continuum test is invalid here.**  A seed *is*
an angular closure: it changes the :math:`\mathcal{O}(N)` **angular
truncation** of the discrete operator.  So the mesh-refinement
(:math:`h \to 0`) limits *at fixed angular order* :math:`N` genuinely
differ between the old and new seed — that difference is not an error, it
is two consistent-but-distinct closures converging to two distinct fixed
angular truncations.  Refining :math:`h` cannot tell them apart.

**The correct discriminator is an angular-order** :math:`N`-**sweep at
fixed mesh** (gate
``test_heterogeneous_1g_angular_order_consistency``).  The retired
edge-extrapolation and the new direct Carlson seed **both converge to the
same transport eigenvalue as** :math:`N \to \infty`: they differ by
:math:`\sim 1.7\times10^{-3}` at Gauss–Legendre order 8 but agree to
:math:`\sim 10^{-6}` by GL32.  A seed that did **not** converge in
:math:`N` would be an *inconsistent* closure (a genuine regression); both
do, so the re-pose is principled.

.. warning::

   **Route (a) is NOT "more accurate" — it is justified structurally.**
   At the *low* angular orders the tests use (GL8), the retired seed is
   actually **closer** to the :math:`N`-limit.  Do **not** frame the
   re-baseline as an accuracy improvement.  Route (a) is justified by
   **structure** — the honest single-pass direct inverse (cold residual
   :math:`5.18\times10^5 \to 2.5\times10^{-16}`) required by the DSA (#2),
   curvilinear-Krylov (#200), and unified-walk (#280) programs — not by
   angular accuracy.

   And the **MMS is blind to the seed**.  Every curvilinear manufactured
   solution in this codebase is :math:`\le` linear-in-:math:`\mu`, which
   is exactly the seed's *exact* regime (:eq:`sn-282-anisotropic-source`
   is the boundary of what the seed can get wrong; vv-principles Mode 7).
   So the MMS-:math:`\mathcal{O}(h^2)` convergence does **not** certify
   the seed — only the :math:`N`-sweep gate does.  (Eigenvalue claims
   need a closed-form or semi-analytical reference regardless; MMS is a
   flux-shape / convergence-order pillar, never an eigenvalue one.)

.. _sn-282-gotchas:

Gotchas
~~~~~~~

**The Krylov restart must be sized from the composite, not the bulk
(ERR-053 family).**  A carve that grows the Krylov composite (adds a
block) **must** resize the GMRES ``restart`` / ``n_dof`` from the
composite ``to_flat().size``, not the bulk :math:`N\cdot n_g\cdot n_x`.
Route (a)'s coupled System B (the ψ½ ray) pushed the coupled ``to_flat``
past the bulk count, so a bulk-sized ``restart`` re-truncated GMRES on the
trace **and** seed degrees of freedom — the restarted subspace could not
represent the coupled iterate and the sphere within-group inner **stalled**
(wrong
:math:`k` under an outer cap, :math:`\sim 868` s).  Fixed at **both**
solver Krylov drivers by sizing ``n_dof = initial_guess.to_flat().size``.
Distinct from #200 (the identity preconditioner); this is a pure
sizing bug.

**The product-cylinder solve consumes the iterate through the
edge-extrapolation stencil — preserve that data flow bit-exactly.**  On a
non-carrying level the seed is still the 2-point extrapolation of the
*iterate* (:math:`t = 0` on product rules: the stencil reads the first
ordinate's iterate column).  That is a *formal* lag, harmless at the fixed
point, and the retirement of the strategy zoo had to keep the
non-carrying data flow **byte-identical** to the pre-2.5d path — a
diverging cylinder solve would trip the §16.D cylinder-unmoved baseline.

**G-reciprocity catches a seed-row error (Mode 12 — CLOSED, ERR-067).**
Under the *retired* **ghost** :math:`G_{\rm sd} = 0` the seed block
carried zero metric weight, so its rows lay **inside** the
metric-weighted G-adjoint reciprocity functional's **invariance group**:
a sign flip on the seed rows left
:math:`\langle A\psi,\chi\rangle_G = \langle\psi, A^{\dagger}\chi\rangle_G`
**exactly** unchanged, at every tolerance, in every regime — a false
green (the classic Mode-12 instance).  The state metric
:math:`G_{\rm sd} = V_{\rm cell}` moves the seed rows **out** of that
invariance group.  With a nonzero ψ½ seed, a seed-row
(:math:`A_{\rm ss}`) sign flip on the forward operator — but not on
:math:`A^{\dagger}`'s independently-coded reverse mode — perturbs
:math:`\langle A\psi,\chi\rangle_G` by
:math:`\sum V_{\rm cell}\,(A_{\rm ss}\psi_{\rm seed})\,\chi_{\rm seed}`
that the unflipped adjoint cannot match, so **reciprocity now REDs** (the
Mode-12-closure gate
``test_mode12_g_reciprocity_catches_a_seed_row_flip``).  The gate carries
**both legs**: a control leg — the *unmutated* nonzero-seed reciprocity
holds :math:`< 10^{-12}`, proving the baseline is the honest
:math:`V_{\rm cell}` adjoint, so the mutated RED is attributable to the
flip (a reverted ghost :math:`G_{\rm sd} = 0` would *also* leave a defect
:math:`\approx 0.107`, a broken baseline mimicking "caught") — and the
mutated leg (:math:`> 10^{-6}`).  This metric-level catch **complements**,
and does not replace, the **object-level** pins that fix the seed
*coefficients*: the C(i) cold residual (forward) and the 2.5b Euclidean
:math:`M^{\mathsf T}` oracle (transpose).  The lesson stands, now
positively — a Mode-12 blindness closes either by gating the object OR by
repairing the functional's **metric** so the error class leaves its
invariance group; here the metric *was itself the bug*, so the correctness
fix and the Mode-12 closure are one and the same.


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
   linear unclamped :math:`\tau^{\rm raw}`.  The weight :math:`\tau` is
   single-sourced in the pole-angular closure
   (:func:`~orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`,
   since Issue #236 Step C — see :ref:`sn-tau-c-on-cellvisit-live`) and
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
:eq:`morel-montry-clamp` in :doc:`/theory/foundations/structured_geometry` for the
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
  :func:`~orpheus.sn.coupled_system.build_within_group_system` (the
  :math:`S` gain of the :math:`(L+C),\,S,\,B` decomposition carries
  :math:`P_1` when
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

The invertible loss operator :math:`A = L + C` (streaming
:math:`L = \mu_n\nabla` plus collision :math:`C = \Sigt{}`) is formed
explicitly via
:class:`~orpheus.sn.operators.streaming.InvertibleOperator`:

.. math::

   A\psi = \mu_n \nabla\psi + \Sigt{}\psi

For Cartesian geometry, this is a banded matrix with upwind cell-
centre finite-difference gradients.  For curvilinear geometries, the
operator includes the :math:`\Delta A/w` geometry factor and the
**symmetric** Morel--Montry angular face-flux interpolation
:math:`\psi_{n\pm 1/2} = \tau\,\psi_{\rm next} + (1-\tau)\,
\psi_{\rm this}`, which is distinct from the WDD asymmetric closure
:math:`\psi_{n+1/2} = (\overline{\psi} - (1-\tau)\,\psi_{n-1/2})/\tau`
used by the sweep.

The system :math:`A\,\psi = b` is solved with
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

- The sweep uses diamond-difference (DD): :math:`A^{-1}` is applied
  implicitly via forward substitution along sweep-direction visits.
- The Krylov path forms :math:`A` explicitly via the symmetric-
  closure FD operator and inverts it via GMRES.

On coarse meshes, DD and FD have different truncation error
constants, so the two paths may give slightly different :math:`\keff`
values.  They converge to the same answer as :math:`h \to 0`.

The curvilinear FD operator reads ``alpha_half`` / ``redist_dAw``
(spherical) and the per-level analogues (cylindrical) from
``SNMesh.reduced``, and the Morel--Montry weight :math:`\tau` from the
pole-angular closure (since Issue #236 Step C — see
:ref:`sn-tau-c-on-cellvisit-live`).  This ensures both paths share
exactly the same connection-coefficient infrastructure.

Round 2 deviation
-----------------

The campaign plan called for an automatic ``"krylov"`` dispatch on
curvilinear meshes in
:func:`~orpheus.sn.solver.solve_sn_fixed_source` that would silently
close ERR-026 on the curvilinear vacuum-BC MMS cases.
Implementation surfaced an unforeseen coupling: the then-
``build_equation_map_spherical`` /
``build_equation_map_cylindrical`` packed-vector layout that
:meth:`InvertibleOperator.apply` reused was designed for the
**reflective** outer-boundary BC only — it has no slot for a vacuum-
BC outer-incoming :math:`\psi`.  Wave E Round 3 owns the equation-map
extension that closes the vacuum-BC path; Round 2 ships the krylov
inner solver as an explicit opt-in for the reflective-BC eigenvalue
case where it is bit-identical math to the legacy BiCGSTAB FD path.


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
   added to the isotropic source on a per-ordinate basis and consumed by
   the within-group sweep (the resolvent ``solve``).

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
   :doc:`/theory/foundations/operator_algebra` for the derivation, the geometry restriction,
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

3. **Sweep** (the within-group resolvent ``solve``,
   :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`): applies
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
predicate-typed matrix-free operator algebra (see :ref:`operator-algebra`).
The neutron transport eigenvalue
problem written in operator form is

.. math::

    (A - S - F)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (A - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`A = L + C` is the invertible loss (streaming + collision)
operator, :math:`S` is the scattering source operator, and :math:`F` is
the fission source operator. Wave D Issue 13 lifts :math:`S` and :math:`F` out of
:class:`~orpheus.sn.solver.SNSolver` and into
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` and
:class:`~orpheus.transport.operators.fission.FissionOperator` respectively. The math is
**moved verbatim** --- the regression contract on the 11 frozen
snapshots at ``tests/sn/regression/snapshots/`` gates the extraction.

Why ``apply``-only
------------------

Both operators report ``is_invertible = False`` and carry no ``solve``
— they are *structural* non-invertibles, declaring no ``inverse()`` /
``solve`` at all (:ref:`design-c-structural-value-split`). (They *are*
adjointable — ``is_adjointable = True`` — since each gained a working
``apply_transpose`` for the outer-layer adjoint / DSA posing; the
"``apply``-only" name refers to the **inverse** axis.)

* **Scattering (:math:`S`)** is rank-:math:`O(N_{\text{cells}}\cdot
  N_{\text{groups}})`. There is no efficient inverse --- the operator
  is *applied*, never *inverted*. The upper-triangular structure that
  would make a sweep-based ``solve`` tractable does not survive the
  Pℓ Galerkin reconstruction. An algebraic consumer that asks for
  :math:`S^{-1}` cannot even spell ``S.inverse()`` (the method is not
  declared --- a *static* error), never silently wrong results at call
  time --- this is the load-bearing payoff of the three-layer operator
  surface (see :ref:`operator-algebra`).

* **Fission (:math:`F`)** has rank-1-in-energy structure: the action
  factorises as :math:`(F\phi)_g = \chi_g\,\sum_{g'}\nu\Sigma_{f,g'}
  \phi_{g'}`, an outer product of the emission spectrum with a scalar
  per-cell rate. This rank-1 structure forbids a useful inverse on
  the energy axis (the rate has lost direction information). The
  :math:`1/k` eigenvalue division stays at the **solver** level ---
  :meth:`~orpheus.transport.operators.fission.FissionOperator.apply` returns
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

   The Pℓ Galerkin reconstruction is realised by the
   spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`
   (Frame/Basis carve), built on a quadrature via
   :meth:`Quadrature.angular_frame(L)
   <orpheus.numerics.quadrature.Quadrature.angular_frame>`. Its
   ``analysis`` face is the :math:`\Pi = Y^* W` projection on the
   angular axis; its ``reconstruction`` face is the addition-theorem
   reconstruction with the :math:`(2\ell+1)` factor. The full-space
   projector is the tensor product
   :math:`\Pi \otimes I_x \otimes I_y \otimes I_g`, built via the
   ``&`` dunder of
   :class:`~orpheus.numerics.operator.TensorProductOperator`. See
   :ref:`galerkin-projection` for the discrete-frame narrative and
   the cross-method consumer table
   (PN solver, energy condensation, MC adjoint moments) and
   :ref:`spherical-harmonics` for the convention and addition
   theorem. The :math:`Y_\ell^m` evaluator
   :meth:`SphericalHarmonicBasis.evaluate
   <orpheus.numerics.basis.SphericalHarmonicBasis.evaluate>`
   is the canonical generic infrastructure consumed here.

   **Wave 1 (commit ff454f2)**:
   :meth:`~orpheus.transport.operators.scattering.ScatteringOperator.build_aniso_source`
   is now the literal §9 line 1230 operator-algebra composition

   .. math::

      Q^{\rm aniso}_n(\vec r) \;=\; R\,\Lambda\,M\,\psi

   where :math:`\Lambda` is :class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`
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
   :file:`tests/sn/verification/mms/test_mms_aniso.py` Pℓ MMS convergence suite.

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
:meth:`~orpheus.numerics.quadrature.Quadrature.spherical_harmonics`
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
:math:`(A, S, F)` into :class:`SNSolver` and replaced the legacy
BiCGSTAB inner-solver path with Krylov-on-:meth:`A.apply` (GMRES
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
(Wave B Issue 7), so the then-``solution_to_angular_flux*`` codec and the
matvec helpers dispatched boundary fills via the realiser-routed
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
:class:`~orpheus.sn.operators.streaming.InvertibleOperator` as the unified
:class:`~orpheus.numerics.operator.LinearOperator` for the invertible
loss operator
:math:`A = L + C = \Omega\cdot\nabla + \Sigma_t`.  This is the Wave D
**capstone**: with :math:`A`, :math:`S`, and :math:`F` all carrying
the Wave A operator-algebra contract, the operator-form transport
equation

.. math::

    (A - S - F)\,\psi = q
    \qquad\text{(fixed source)}

.. math::

    (A - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

is now expressible in ORPHEUS as a single Python expression composed
from three :class:`LinearOperator` objects.  The Wave A composers
(:class:`~orpheus.numerics.operator.OperatorSum`,
:class:`~orpheus.numerics.operator.OperatorProduct`,
:class:`~orpheus.numerics.operator.ScaledOperator`) compute their
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
:attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`
recursively from the constituents per the closure laws (see
:ref:`operator-algebra`).

The full operator surface — apply, solve, apply_transpose
---------------------------------------------------------

:class:`InvertibleOperator` realizes all three verbs — ``apply``,
``solve``, ``apply_transpose`` — and reports both ``is_invertible`` and
``is_adjointable`` ``True`` (the retired ``{"apply", "solve",
"apply_transpose"}`` capability frozenset advertised exactly this;
:ref:`capability-set-semantics`):

* :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply` —
  matrix-free forward action :math:`A\,\psi = (L+C)\,\psi` via the
  operator's own
  :attr:`~orpheus.sn.operators.streaming.InvertibleOperator.loss_representation`
  through the shared apply-direction walk
  (:meth:`~orpheus.sn.loss_representation.LossRepresentation.loss_action`,
  the ``(L+C)ψ`` single emission — the apply-direction twin of
  :meth:`~orpheus.sn.loss_representation.LossRepresentation.sweep`,
  L21 "matvec ≡ sweep"; #206 Phase C).  (Historically ``apply``
  reused the symmetric-closure math extracted verbatim from the
  per-geometry ``transport_operator_matvec*`` functions of the
  BiCGSTAB FD operator; that whole family and its unified successor
  were deleted in the typed-field (#197) and walk-unification
  (#280 campaigns) refactors.)

* :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve` —
  inverse action :math:`A^{-1}\,q` via the operator's own
  :attr:`~orpheus.sn.operators.streaming.InvertibleOperator.loss_representation`
  sweep (the WDD forward substitution; the 1-D scan or multi-D wavefront
  selected by :func:`~orpheus.sn.loss_representation.default_for`).

* :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.apply_transpose` —
  adjoint action :math:`A^{\mathsf T}\,\phi` via the loss-representation's
  named
  :meth:`~orpheus.sn.loss_representation.LossRepresentation.loss_action_transpose`
  (the reverse-direction walk carrying the second angular triangular
  factor on curvilinear; the multi-D Cartesian adjoint is an honest
  deferral raise, never a silent wrong answer).  It gates the
  reciprocity invariant :math:`\langle A\,\psi,\,\varphi\rangle =
  \langle\psi,\,A^{\mathsf T}\,\varphi\rangle` (see "Reciprocity"
  subsection below).  (The pre-Depth-B path assembled a dense matrix by
  probing :meth:`apply` with unit basis vectors and returned its explicit
  transpose; that retired with the bundled ``SNStreamingOperator`` class
  in D-K.)

Apply and solve use **different** closures by design (historical)
-----------------------------------------------------------------

.. note:: **Superseded architecture (Wave D / early Wave E).**  This
   subsection describes a design in which ``apply`` (a separate
   finite-difference operator) and ``solve`` (the WDD sweep) were
   **two distinct discretisations** of :math:`A`.  That split was
   dissolved by the #206 Phase C **matvec ≡ sweep** unification: today
   ``apply`` and ``solve`` run the **one** loss-representation walk
   (:meth:`~orpheus.sn.loss_representation.LossRepresentation.loss_action`
   for the matvec, :meth:`~orpheus.sn.loss_representation.LossRepresentation.sweep`
   for the inverse — "one walk", a code fact per L21), so there is no
   longer a separate FD operator and no by-design bit-difference between
   the two directions.  The FD-operator family
   (``transport_operator_matvec*`` and its unified successor) was deleted
   in the typed-field (#197) and walk-unification (#280 campaigns)
   refactors.  The historical two-closure narrative below is retained for
   the ERR-026 closure-bias reasoning it records.

This was the load-bearing architectural fact about
:class:`InvertibleOperator`, and the reason the operator's
:meth:`apply` was **not** bit-identical to its :meth:`solve`.

The Wave-D-era SN dispatch in ORPHEUS shipped **two distinct
discretisations** of the same continuous operator
:math:`A = \Omega\cdot\nabla + \Sigma_t`, built at different times
for different consumers:

* The **finite-difference operator** (the ``transport_operator_matvec_*``
  family in :mod:`orpheus.sn.operators.streaming`, since deleted) was
  built for the BiCGSTAB inner
  solver path (:meth:`SNSolver._solve_bicgstab_*`).  It used
  upwind cell-center FD on Cartesian and arithmetic face averages
  with **τ-symmetric Morel-Montry angular interpolation** on
  curvilinear (see the "Explicit Transport Operator" subsection of
  the BiCGSTAB Alternative above and the warning at the head of
  :mod:`orpheus.sn.operators.streaming`).

* The **sweep operator**
  (:meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`,
  dispatching its selected representation through the
  :class:`~orpheus.transport.spatial.scheme.DiscretizationScheme`
  Protocol with :class:`~orpheus.transport.spatial.diamond.DiamondDifference`
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
  the head of :mod:`orpheus.sn.operators.streaming`).
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
   operator surface (it reports ``is_invertible = True``).

2. **The composers (Wave A) need a uniform contract.**  When a
   downstream agent composes :math:`(A - S - F)`, the composition's
   :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` /
   :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable` are
   derived recursively from each operand.  If
   :class:`InvertibleOperator` shipped only :meth:`apply`, the
   composition would report ``is_invertible = False`` and the Wave E
   Krylov-on-apply path could not request a sweep-preconditioned matvec
   without bypassing the algebra.

Reciprocity invariant
---------------------

The reciprocity identity is the defining property of the operator-
transpose pairing under the discrete L\ :sup:`2` inner product:

.. math::
    :label: sn-streaming-reciprocity

    \langle A\,\psi,\,\varphi\rangle \;=\;
    \langle\psi,\,A^*\,\varphi\rangle

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
:func:`tests.sn.operators.test_streaming_operator.TestLinearity.test_apply_is_linear` catches
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

Today :meth:`apply`, :meth:`apply_transpose`, and :meth:`solve` all
operate on the **same** typed composite carrier
:class:`~orpheus.transport.full_field.FullField` (bulk
:class:`~orpheus.transport.fields.angular_flux.AngularFlux` +
boundary :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`),
in the principled ``(N, ng, nx, ny)`` layout
(see :ref:`theory-sn-index-convention`).  The source and its state
ride on the composite's typed members:

* **Source** — carried as the composite bulk
  (``rhs.interior.values``, per-ordinate shape ``(N, ng, nx, ny)``);
  the P\ :sub:`ℓ` (:math:`\ell\ge 1`) anisotropic contribution is
  folded into this one per-ordinate source (the separate ``Q_aniso``
  kwarg is gone).
* **Boundary state** — carried as the typed
  :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
  face views on ``rhs.boundary`` (keyed by face name), replacing the
  former stringly-typed ``psi_bc`` dict; the sweep seeds its mutable
  boundary buffer from these trace slots and persists reflective/pole
  state through them between outer iterations.

.. note:: **Superseded packed-vector layout (Wave D / early Wave E).**
   The previous design gave :meth:`apply` / :meth:`apply_transpose` a
   **packed 1-D solution vector** (an ``EquationMap`` selecting the
   unknown ``(ordinate, cell)`` combinations, laid out group-major in
   Fortran order for :func:`scipy.sparse.linalg.bicgstab`) that
   differed from :meth:`solve`'s structured arrays, and ``solve``
   consumed a separate ``Q`` / ``psi_bc`` dict / optional ``Q_aniso``
   triple.  The typed-field campaign (#197, then Depth-B D-H/D-I/D-J)
   retired the packed-vector convention, the ``EquationMap`` codec
   family, and the ``psi_bc`` dict in favour of the single
   :class:`~orpheus.transport.full_field.FullField` composite above; the
   ``Q_aniso`` kwarg folded into the one per-ordinate source.  There is
   no longer a layout difference between :meth:`apply` and
   :meth:`solve`.

Why not extract :meth:`apply_transpose` analytically?
-----------------------------------------------------

For a continuous operator :math:`A = \Omega\cdot\nabla + \Sigma_t`,
the L\ :sup:`2`-adjoint is :math:`A^* = -\Omega\cdot\nabla +
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


.. _sn-scattering-adjoint:

The scattering adjoint, free from the harmonic frame
----------------------------------------------------------

The streaming operator's analytic adjoint is hard (the subsection above):
sign-flipping the upwind direction, transposing the M–M closure, re-deriving
the per-level azimuthal redistribution — each an AI-failure-mode trap — so
:math:`A^*` is taken by the dense-transpose fallback.  The **scattering**
operator :math:`S` is the counterexample: campaign **#276 P3** (commit
``15185e5``, closes
`#118 <https://github.com/deOliveira-R/ORPHEUS/issues/118>`_) made
:math:`S^{T}` fall out **for free**, because :math:`S` is already written as
a harmonic-frame conjugation.

The modernised in-scatter source is ONE frame-conjugated operator
(:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.full_scatter_kernel`):

.. math::

   \mathrm{full\_scatter\_kernel}
   \;=\; R \circ (\Lambda_{\ell\ge 0} + N_{2n}) \circ M ,

where :math:`M` / :math:`R` are the angular frame's analysis /
reconstruction faces, :math:`\Lambda_{\ell\ge 0}` is the per-:math:`\ell`
moment-space group transfer
(:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`),
and :math:`N_{2n}` is the distinct :math:`(n,2n)` multiplication channel
(:class:`~orpheus.transport.operators.scattering.N2NMomentOperator`) —
summed with :math:`\Lambda` in moment space and conjugated by the frame
*together* (one analysis, one reconstruction) for the WHOLE
P0 + :math:`\ell\ge1` + :math:`(n,2n)` source.  Its transpose is therefore
the product transpose

.. math::

   \mathrm{full\_scatter\_kernel}^{T}
   \;=\; M^{T} \circ (\Lambda + N_{2n})^{T} \circ R^{T},

which :meth:`OperatorProduct.apply_transpose
<orpheus.numerics.operator.OperatorProduct.apply_transpose>` assembles from
the leaf transposes — the frame's :math:`M^{T}` / :math:`R^{T}` faces (landed
in the Frame/Basis carve), the per-:math:`\ell` group transpose
:math:`\Lambda^{T}`, and :math:`N_{2n}^{T}` — with **no per-geometry
derivation to verify** (the trap the streaming adjoint above could not
avoid).  The per-ordinate adjoint scattering source is then

.. math::

   S^{T}\chi \;=\; \tfrac{1}{W}\,\mathrm{full\_scatter\_kernel}^{T}\,\chi ,

the producer-side :math:`1/W` transposing as the scalar it is
(:math:`(A/W)^{T} = A^{T}/W`).
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` now
reports ``is_adjointable = True`` (it has a working ``apply_transpose``),
and the old "no ``apply_transpose``" class-docstring confession is
retired.

**Forward fast-path, adjoint frame-path — and why the asymmetry is
principled.**  The production FORWARD source keeps the scalar fast-path
(:attr:`~orpheus.transport.operators.scattering.ScatteringOperator.isotropic_kernel`
for P0 + :math:`(n,2n)`, and the per-:math:`\ell` ``build_aniso_source``)
for SI-sweep performance; the **adjoint** — not the hot path — rides the
validated frame form instead.  The two are thus structurally *different*
representations of the same operator, which is exactly what makes the
verification a genuine cross-check rather than a tautology: the per-group
Euclidean reciprocity
:math:`\langle S\psi, \chi\rangle = \langle\psi, S^{T}\chi\rangle`
(``tests/sn/operators/test_scattering_adjoint.py``,
``TestFullScatterKernel::test_S_euclidean_reciprocity``) pins the frame-form
:math:`S^{T}` against the *independent* scalar fast-path :math:`S`, and the
forward equivalence
:math:`(1/W)\,\mathrm{full\_scatter\_kernel}.\mathrm{apply} \equiv
S.\mathrm{apply}` holds to :math:`\sim 10^{-12}`.

.. note::

   This :math:`S^{T}` is the **Euclidean** transpose (the plain
   group-and-angle matvec adjoint), NOT the metric Hilbert adjoint
   :math:`S^{\dagger} = G^{-1}S^{T}G` — that angular-Gram weighting is the
   :attr:`~orpheus.numerics.operator.LinearOperator.H` wrapper's job.  The
   campaign and commit name it "S†" colloquially; the precise object the
   operator computes is the transpose.

This is the discrete scattering adjoint the SN adjoint chain builds on: the
adjoint flux :math:`\psi^{*}` solving :math:`(L+C-S)^{T}\psi^{*} = q^{*}`,
adjoint-weighted homogenisation, perturbation theory, and detector
sensitivity all need :math:`S^{T}`.  Its companion forward step (campaign
**#276 P2**, commit ``dcea43a``) routes the SN forward *isotropic* source
through the same model-shared
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicScattering`
(:math:`\Sigma_{s0}`) and
:class:`~orpheus.transport.operators.isotropic_scattering.IsotropicN2N`
(:math:`2\Sigma_{2n}`) operators (0-ULP bit-identical), so the
:math:`K_\mathrm{iso}` energy operator — which also assembles the
infinite-medium loss matrix (:ref:`direct-eigensolve-assembly`) — is one
cross-model source.  These model-shared operators live in
:mod:`orpheus.transport.operators`.


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


.. _sn-angular-windowing:

Angular windowing — the SI iterate lives in moment space (Phase 5a)
===================================================================

Wave O step #205 **Phase 5a** (commits ``93807aa`` factoring / ``b97d4f9``
eigenvalue inner / ``13ca001`` fixed-source inner, 2026-06-07) is a
**moment-reduction** of the SN within-group source-iteration *iterate*.
It is **orthogonal** to the :ref:`interior face-flux cochain <wavefront-flux-cochain>`
above: where that types the *per-ordinate* interior face flux a single
sweep propagates (and explicitly frames the interior cochain as
per-ordinate, :math:`\psi^{(1)}_\Omega \in C^1`), Phase 5a observes that
the **persistent** iterate the source iteration carries *between* sweeps
does not need all :math:`N` ordinates — it needs only the
spherical-harmonic moments the scattering operator consumes. The
held iterate's angular dimension drops :math:`N \to (L+1)(2L+1)`, and
the iterate becomes :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
instead of :class:`~orpheus.transport.fields.angular_flux.AngularFlux`.

.. admonition:: Key Facts (angular windowing)
   :class: tip

   - **The SI fixed point lives in moment space.** Within-group source
     iteration is :math:`\psi_{k+1} = (L{+}C)^{-1}(S\,\psi_k + B\,\psi_k
     + q)`. The scattering source :math:`S\,\psi` is a pure function of
     the flux moments :math:`\phi_\ell^m = (M\psi)_\ell^m` — the
     per-ordinate iterate :math:`\psi` carries strictly more than the
     iteration consumes.
   - **Hold the moments, not the ordinates.** The persistent iterate is
     the moment tensor :math:`\phi \in (L{+}1, 2L{+}1, n_g, n_x, n_y)`,
     not :math:`\psi \in (N, n_g, n_x, n_y)`. Measured **18.3×**
     persistent-iterate shrink at :math:`N = 110`, :math:`L = 1`
     (:math:`N / (L{+}1)(2L{+}1) = 110/6`).
   - **The per-step source is bit-identical.** :math:`S` consuming the
     moments equals :math:`S` consuming the full angular flux **bit for
     bit** (0 ULP) under the ORPHEUS unnormalized-harmonic convention.
     The only non-bit-identical change is the SI *convergence test*,
     which moves to the moment :math:`L^2` — *more* principled, not a
     regression.
   - **2-D Cartesian only (load-bearing).** Windowing is valid where the
     sweep is a **direct** solve with no per-ordinate-iterate seed.
     Curvilinear (1-D sphere / cylinder) **must** stay full-angular: the
     Morel–Montry Carlson coupled-pole closure seeds from the previous
     iterate's per-ordinate :math:`\psi` at :math:`\mu = -1`, which the
     moment tensor cannot reconstruct. The Krylov path stays
     full-angular too. Gated on the genuine ``sn_mesh.is_cartesian and
     sn_mesh.ndim == 2`` (C5.4 / #225 — the earlier ``reduced is None``
     proxy was also true at 3-D Cartesian).
   - **Interior-bulk only.** The reflective :math:`B` coupling reads the
     full per-ordinate boundary *trace*; windowing reduces only the
     interior bulk. The biproduct :eq:`wavefront-cochain-biproduct`
     keeps the trace :math:`C^1_\partial` a distinct, **un-reduced**
     summand.
   - **A representation + typed-state win, NOT yet a peak-memory win.**
     5a shrinks the *persistent* iterate (18.3×) and makes its type
     honest. The *peak* memory drops only modestly (~1.2× measured)
     because the per-sweep **transient** full-angular machinery still
     dominates — that transient is the target of Phase 5b (interior-face
     cochain) and Phase 5c
     (:ref:`full-angular output <sn-angular-windowing-in-sweep-accumulation>`,
     the 3.06× linear peak win).


.. _sn-angular-windowing-fixed-point:

The within-group fixed point lives in moment space
--------------------------------------------------

The within-group source iteration solves, for each outer step, the
fixed-point problem

.. (vv-status rationale) The within-group source-iteration fixed point.
   This is the governing iteration the windowing reorganizes; the
   verifiable content is the SI ≡ Krylov cross-check (Krylov stays
   full-angular) and the closed-form k_inf eigenvalue, not the rendered
   recurrence. Documented, not orphan-gated.
.. vv-status: si-within-group-fixed-point documented

.. math::
   :label: si-within-group-fixed-point

   \psi_{k+1}
   \;=\; (L + C)^{-1}\!\left( S\,\psi_k + B\,\psi_k + q \right),

where :math:`L + C` is the within-group **invertible resolvent** (the
streaming + collision the sweep inverts directly), :math:`S` is the
within-group scattering gain (:ref:`pn-scattering` in
:doc:`/theory/methods/sn/index`), :math:`B` is the reflective boundary
coupling (:ref:`bc-extraction`), and :math:`q` the fixed external /
fission source. Both :math:`S` and :math:`B` are **lagged gains** — the
sweep never re-scatters mid-sweep (cf. the variadic driver,
:ref:`bc-extraction-variadic-driver`).

The load-bearing observation is the **arity of the scattering gain**.
:math:`S\,\psi` depends on :math:`\psi` *only* through its
spherical-harmonic flux moments. Writing the moment-projection operator
:math:`M` (the :eq:`flux-moments` quadrature contraction, the SH frame's
analysis face ``frame.analysis`` from
:meth:`Quadrature.angular_frame(L)
<orpheus.numerics.quadrature.Quadrature.angular_frame>`)

.. math::
   :label: angular-windowing-moment-projection

   \phi_\ell^m(\vec r)
   \;=\; (M\psi)_\ell^m(\vec r)
   \;=\; \sum_{n=1}^{N} w_n \, Y_\ell^m(\hat\Omega_n)\,
         \psi_n(\vec r),
   \qquad 0 \le \ell \le L,\; |m| \le \ell,

the scattering source factors **through the moment boundary**:

* the **isotropic** :math:`\ell = 0` (P0) and the **(n,2n)** doubling
  terms (:ref:`pn-scattering`) need only the scalar flux
  :math:`\phi_0 \equiv \phi_0^0`;
* the **anisotropic** :math:`P_{\ell\ge 1}` term needs the higher
  moments :math:`\phi_\ell^m` up to the scattering order :math:`L`.

So the per-ordinate iterate :math:`\psi \in \mathbb{R}^N` is mapped, at
the very first thing the sweep's source assembly does, onto the
:math:`(L{+}1)(2L{+}1)`-dimensional moment space — and **nothing
downstream of that projection ever reads the discarded
:math:`N - (L{+}1)(2L{+}1)` angular degrees of freedom**. The iterate
carries strictly more than the iteration consumes. Angular windowing
holds the iterate at the consumed dimension: the persistent state is
the moment tensor

.. math::
   :label: angular-windowing-moment-iterate

   \phi \;\in\; \mathbb{R}^{(L+1)\times(2L+1)\times n_g \times n_x \times n_y}
   \quad(\texttt{HarmonicMomentFlux}),
   \qquad\text{not}\qquad
   \psi \;\in\; \mathbb{R}^{N \times n_g \times n_x \times n_y}
   \quad(\texttt{AngularFlux}).

.. vv-status: angular-windowing-moment-iterate documented

For :math:`N = 110` (Lebedev order 17) and :math:`L = 1`
(:math:`(L{+}1)(2L{+}1) = 6`) the angular dimension drops **18.3×**.


.. _sn-angular-windowing-factoring:

The scattering factoring — :math:`S_{\rm aniso} = \tfrac{1}{W}\,R\,\Lambda\,M`
------------------------------------------------------------------------------

Phase 5a's commit 1 (``93807aa``) makes the factoring *structural* so
that the windowed and full-angular paths share one source of truth. The
anisotropic in-scatter is the §9 operator composition (the §15.2
sum-of-tensor-products form,
:meth:`~orpheus.transport.operators.scattering.ScatteringOperator.build_aniso_source`)

.. (vv-status rationale) The R·Λ·M anisotropic-scattering factoring.
   A structural identity (associativity of the three-operator
   composition); the verifiable content is the bit-identity of the two
   evaluation arms, pinned by the de-risk probe Q2/Q3 (0 ULP) and the
   independent Bell & Glasstone hand reconstruction Q2b (1.5 ULP).
   Documented.
.. vv-status: angular-windowing-aniso-factoring documented

.. math::
   :label: angular-windowing-aniso-factoring

   Q^{\rm aniso}_n(\vec r)
   \;=\; \frac{1}{W}\,\bigl(R\,\Lambda\,M\,\psi\bigr)_n(\vec r),
   \qquad W = \sum_n w_n,

where (reading right to left):

* :math:`M` is the moment **projection** :eq:`angular-windowing-moment-projection`
  (the SH frame's analysis face ``frame.analysis``);
* :math:`\Lambda` is the per-:math:`\ell` block-diagonal scattering on
  moment space :math:`\Lambda = \sum_\ell P_\ell \otimes \Sigma_{s,\ell}`
  (:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`);
* :math:`R` is the addition-theorem **reconstruction** with the
  :math:`(2\ell+1)` factor
  (the SH frame's reconstruction face ``frame.reconstruction``);
* :math:`1/W` is the producer-side normalization applied at the
  :meth:`~orpheus.transport.operators.scattering.ScatteringOperator.apply` boundary.

The associativity :math:`(R\,\Lambda)\,M` is the whole trick. The
**full-angular** path evaluates the composition as written —
:math:`R\cdot\Lambda\cdot M(\psi)`: project, scatter, reconstruct,
consuming the §5.6 :attr:`kernel <orpheus.transport.operators.scattering.ScatteringOperator.kernel>`
``= frame.conjugate(Λ)`` as a single composed ``np.ndarray`` operator. The
**windowed** path's iterate bulk **is** the moments :math:`\phi = M\psi`,
so :math:`M` is *already done* and only the **moment → source** map
:math:`R\,\Lambda` remains: :math:`R\cdot\Lambda(\phi)`. The dispatch is on
the iterate type: the
:class:`~orpheus.transport.fields.angular_flux.AngularFlux` arm of
:meth:`ScatteringOperator.apply <orpheus.transport.operators.scattering.ScatteringOperator.apply>`
does the :math:`M` projection first; the
:class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`
arm skips it.

.. note::

   **P4 — the windowed arm now takes the explicit typed edges.** Earlier
   the two arms shared the single ``np.ndarray`` primitive
   :meth:`~orpheus.transport.operators.scattering.ScatteringOperator._aniso_source_from_moment_values`
   (``frame.reconstruct_after(Λ)``). The Frame campaign's P4 phase
   re-expressed the windowed arm's :math:`R\,\Lambda` as the **explicit
   typed** carrier path
   ``frame.reconstruct(Λ.apply(φ)) / W`` — :math:`\Lambda` materialises a
   typed :class:`~orpheus.transport.source_sinks.harmonic_moment_source_sink.HarmonicMomentSourceSink`
   (the role-changing edge of the carrier grid,
   :ref:`scattering-carrier-grid`), which the frame's :math:`R` face then
   reconstructs to an :class:`AngularSourceSink`. It threads the *same*
   :math:`\Lambda` kernel and the *same* frame :math:`R` face as the fused
   path, so the two agree numerically; the ndarray
   ``reconstruct_after(Λ)`` primitive is **retained as the 0-ULP
   crosscheck oracle** the typed arm is pinned against, not deleted. See
   :ref:`scattering-carrier-grid` for the explicit-vs-fused design choice.

This factoring **retired** the per-:math:`\ell`
``_PerLegendreOrderScattering`` kernel, which recomputed :math:`M\psi`
independently for every Legendre order — an :math:`L`-fold redundant
projection (aggressive-retirement discipline).

.. admonition:: The :math:`Y_0^0 = 1` convention — the scalar flux is read off
   :class: note

   ORPHEUS uses **unnormalized real harmonics** (the
   "no-:math:`4\pi/(2\ell+1)`-prefactor" convention,
   :ref:`spherical-harmonics`), under which :math:`Y_0^0 = 1` *exactly*.
   Therefore the :math:`\ell = 0` moment **is** the scalar flux,

   .. math::

      \phi_0 \;=\; \sum_n w_n \, Y_0^0(\hat\Omega_n)\,\psi_n
                 \;=\; \sum_n w_n \, \psi_n
                 \;=\; (\texttt{integrate\_angular}\;\psi),

   read directly off the moment tensor's :math:`\ell{=}0` block with no
   rescale (``phi_moments.l_block(0)[0]``). This is what lets the
   windowed P0 + (n,2n) fast path consume the moments with **zero**
   conversion arithmetic, and what makes the eigenvalue outer's scalar
   flux bit-identical to the full-angular ``integrate_angular`` (the
   de-risk Q1 below proves :math:`\Phi[0,0] = \texttt{integrate\_angular}`
   bit-for-bit).


.. _sn-angular-windowing-geometry-restriction:

Why 2-D Cartesian only — the curvilinear seed obstruction
---------------------------------------------------------

Windowing is valid **only** where the within-group sweep is a *direct*
solve that does not seed from the previous iterate's per-ordinate
:math:`\psi`. There are three regimes, and only one admits it:

.. list-table:: Where angular windowing applies
   :header-rows: 1
   :widths: 26 16 58

   * - Path
     - Windowed?
     - Why
   * - **2-D Cartesian DD** (SI inner)
     - **yes**
     - The diamond-difference wavefront sweep is a direct forward
       substitution down the upwind DAG; it inverts :math:`L+C` from
       :math:`q + S\phi + B\psi_\partial` with **no** interior-iterate
       seed. The bulk seed (``initial_guess``) threads through harmlessly
       (the 2-D sweep ignores it). The moment iterate is sufficient.
   * - **Curvilinear 1-D** (sphere / cylinder)
     - **no**
     - The Morel–Montry **Carlson coupled-pole** closure seeds the
       :math:`\mu = -1` starting-direction angular flux from the
       *previous iterate's per-ordinate* :math:`\psi` (the curvilinear
       angular-redistribution recursion is initialized from the inward
       radial sweep's last iterate). A moment tensor cannot reconstruct
       that per-ordinate seed — windowing would lose the closure's
       starting data. Curvilinear stays full-angular.
   * - **Krylov** (any geometry)
     - **no**
     - GMRES iterates the full bulk vector :math:`\psi`; the Krylov
       subspace is built from full-angular matvecs. There is no moment
       sub-iterate to hold.

The gate is the genuine predicate ``sn_mesh.is_cartesian and
sn_mesh.ndim == 2`` (the C5.4 / #225 sharpening of the earlier
``reduced is None`` proxy, which was *also* true at 3-D Cartesian and
would have silently moment-windowed a 3-D solve — vv Mode 9): the
curvilinear meshes carry a non-``None`` ``reduced`` moment-reduction
descriptor and are excluded, and so is 3-D Cartesian (the in-sweep
moment emit is a 2-D kernel). The windowed product ``P @ A.inverse()``
(retyped in :ref:`windowing-retyped` below) — which the SI driver now
holds and applies **directly**, the ``_maybe_window`` factory returning
the plain sweep off that gate (#226 taxonomy step 3,
:ref:`inverse-application-driver`) — is therefore **never even
constructed** off the 2-D-Cartesian gate, so there is no illegal state
to mistype.

The restriction is **interior-bulk-only** in a second sense: the
interior-cochain biproduct :eq:`wavefront-cochain-biproduct`
:math:`C^1 = C^1_{\rm int} \oplus C^1_\partial` keeps the **boundary
trace** :math:`C^1_\partial` (the :class:`AngularBoundaryFlux` summand) a
distinct, *un-reduced* per-ordinate object. The reflective :math:`B`
coupling reads the full per-ordinate face trace via the typed
:math:`\iota_*` / :math:`\iota^*` exchange — windowing reduces only the
interior bulk and never touches the trace. The
``test_2d_windowed_si_reflective_trace_is_nonzero`` guard pins this (a
windowing that zeroed the trace would be a dropped reflective coupling,
invisible to the interior-only scalar-flux snapshot).


.. _sn-angular-windowing-bit-identity:

Bit-identity of the source, principled-equivalence of the convergence test
---------------------------------------------------------------------------

The carve has a clean correctness story split in two
(``vv-principles`` § "Bit-identity vs principled-equivalence").

**The per-step source is bit-identical.** A de-risk probe
(``derivations/diagnostics/diag_p5a_moment_consuming_scatter.py``, a
diagnostic script not retained in the repo) proved
the **moment arm** :math:`S(M\psi)` equals the **full-angular arm**
:math:`S(\psi)` **bit for bit** (``np.array_equal``, 0 ULP) before any
production code was written. Because the moment-projection operator
:math:`M` inside the windowed resolvent is built from the **same**
quadrature harmonics :math:`Y` and weights :math:`w` the scattering
operator uses internally, the stored moments equal :math:`S`'s own
internal projection of the same :math:`\psi` term-for-term — so the
per-sweep re-projection is not just *elided*, its result is *reproduced
exactly*. The probe was cross-checked against a **structurally
independent** Bell & Glasstone §1.6 hand reconstruction of the P1 source
(every factor — :math:`w_n`, :math:`Y`, the :math:`(2\ell+1)=3`, the
:math:`1/W` — written out by hand from the material data, *not* via the
project's projection primitives), which agreed at ~1.5 ULP (the expected
floating-point distance for an independent reduction order). This is the
L11 structural-independence guard the bit-exact comparison alone lacks.

.. list-table:: De-risk probe — moment arm vs full-angular arm
   :header-rows: 1
   :widths: 30 14 28 28

   * - Probe (config: 2-D P1 het, LS-S4)
     - Groups
     - Result
     - Verdict
   * - **Q1** :math:`\Phi[0,0] = \texttt{integrate\_angular}`
       (the :math:`Y_0^0 = 1` scalar-flux read)
     - 2G
     - max :math:`|\Delta| = 0`, ``np.array_equal`` True
     - bit-exact (0 ULP)
   * - **Q2** :math:`R\Lambda(\phi)` vs
       :math:`R\Lambda M(\psi)` (the aniso :math:`\ell\ge 1` arm)
     - 2G
     - max ULP :math:`= 0`
     - bit-exact (0 ULP)
   * - **Q3** full :math:`S(M\psi)` vs :math:`S(\psi)`
       (end-to-end source)
     - 2G **and** 4G
     - max ULP :math:`= 0` (both)
     - bit-exact (0 ULP)
   * - **Q2b** vs INDEPENDENT Bell & Glasstone hand
       reconstruction (L11 structural-independence ground)
     - 2G
     - max rel :math:`= 3.4\times10^{-16}`
     - principled (~1.5 ULP)
   * - **Q4** non-degeneracy: aniso :math:`> 0`,
       P1 :math:`\ne` P0, asymmetric :math:`\Sigma_s` (ERR-002 ground)
     - 2G
     - aniso max :math:`= 0.49`, :math:`|\Sigma_{s,0}-\Sigma_{s,0}^\top| = 0.1`–:math:`0.18`
     - non-degenerate (gate can see Mode 6)

**The convergence test moves to moment space (principled-equivalence).**
The *only* non-bit-identical change is the SI stopping criterion. The
full-angular path tested :math:`\lVert\psi_{k+1} - \psi_k\rVert_2`; the
windowed path necessarily tests :math:`\lVert\phi_{k+1} -
\phi_k\rVert_2`. This is the **physically-meaningful** SI criterion —
scattering iterates the *moments*, so the moment :math:`L^2` is the
natural convergence norm — i.e. the change is *more* principled, not a
regression. Because the stopping point shifts by a fraction of a
tolerance, the converged value agrees with the full-angular path within
``SAFETY × conv_tol``: the measured drift is **2-D eigenvalue
:math:`k_{\rm eff}` :math:`2.4\times10^{-14}` relative** and **scalar
flux :math:`6.3\times10^{-12}` relative**, both well inside the
regression framework's ``kind="iterative"`` gate (drift bounded by
:math:`(\text{iteration count})\times(\text{condition number})\times
\text{ULP}`, criterion 3 of the principled-equivalence test). The
eigenvalue **outer** converges on the scalar flux, which is bit-identical
via :math:`\phi_0` (the :math:`Y_0^0 = 1` read) — so *only* the inner SI
stopping point shifts; the outer power iteration is untouched.

The one *correctness* anchor (not merely *equivalence*) is the
closed-form homogeneous eigenvalue: the windowed 2-D eigenvalue solve
reproduces :math:`k_\infty = \nu\Sigma_f / \Sigma_a` (the closed-form
pillar — MMS does **not** prove eigenvalues), and the
``test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree`` gate
cross-checks the windowed SI flux against the **full-angular Krylov**
flux (a genuinely independent reference, since Krylov is never
windowed), confirming the windowing does not silently drop the
:math:`\ell\ge 1` moments — the central trap.


.. _sn-angular-windowing-honest-scope:

Honest scope — a persistent-iterate + typed-state win, NOT yet a peak win
-------------------------------------------------------------------------

.. warning::

   Phase 5a reduces the **persistent** iterate storage and makes the SI
   state honest. It does **NOT** by itself deliver the full peak-memory
   reduction. State this carefully — the 18.3× number describes the
   *held* iterate, not the *peak*.

   * **What 5a wins.** The held + warm-started ``_psi_typed`` carried
     across the entire eigenvalue solve (and the convergence-test copy)
     shrinks by :math:`N / (L{+}1)(2L{+}1)` — measured **18.3×** at
     :math:`N = 110`, :math:`L = 1`. The iterate **type** also becomes
     honest: the SI state *is* the moments
     (:class:`HarmonicMomentFlux`), so the representation no longer
     over-claims an angular resolution the iteration never uses.
   * **Why the peak win is modest.** The **per-sweep transient**
     full-angular machinery still dominates the peak: the resolvent's
     swept output, the ``psi_x`` / ``psi_y`` interior cochain
     :math:`C^1_{\rm int}` (the cochain section above), and the per-octant
     ``.copy()`` buffers — several full-angular-sized arrays that
     storage-A still materializes within a single sweep. Measured
     **peak** reduction is only **~1.2×** for a :math:`12\times8` config
     (the held :math:`2\times` iterate restored to full-angular size
     against the windowed ``tracemalloc`` peak;
     ``derivations/diagnostics/diag_p5a_peak_memory.py``, a diagnostic
     script not retained in the repo).

   The per-sweep transient has two components, eliminated by two later
   phases. **Phase 5b** (interior-face storage-B — the rolling
   moving-frontier window over the wavefront anti-diagonals, which never
   materializes the whole interior cochain at once) cuts the *interior
   face* cochain transient and is the **3-D enabler** (a 3-D wavefront's
   full interior cochain is prohibitively large to materialize; the moving
   frontier is the only tractable representation). **Phase 5c**
   (:ref:`sn-angular-windowing-in-sweep-accumulation`) cuts the
   *full-angular output* transient by accumulating the harmonic moments
   per anti-diagonal inside the sweep — the **linear peak win** (measured
   3.06×), and the realization of the full peak reduction this honest-scope
   warning deferred.

So Phase 5a is precisely: the **persistent-iterate** shrink + the
**typed-state** correctness win + the **foundation** for 5b/5c. Phase 5b
eliminates the **interior-face transient** and is the **3-D** enabler;
Phase 5c eliminates the **full-angular output transient** (the linear peak
win). The three together deliver the asymptotic peak reduction the moment
iterate makes possible; 5a alone delivers the type and the
persistent-storage win.


.. _sn-angular-windowing-implementation:

Implementation map
-------------------

* The solver's :func:`_maybe_window <orpheus.sn.solver._maybe_window>`
  builds the typed windowed composition ``P @ A.inverse()`` — the
  :class:`~orpheus.sn.operators.windowing.WindowedSweep` — over the
  within-group forward :math:`A`'s inverse (the Jacobi
  :math:`(L+C)^{-1}`, or the reified splitting matrix :math:`M^{-1}`,
  :class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`
  — see :ref:`si-gauss-seidel-reification` in :doc:`/theory/methods/sn/index`)
  and hands it to the :class:`~orpheus.numerics.iteration.SourceIteration`
  driver, which **applies** it directly (#226 taxonomy step 3 —
  :ref:`inverse-application-driver`; the transitional
  ``_MomentWindowedResolvent`` ``.solve`` adapter that once wrapped it is
  gone). Its analysis factor :math:`P` is sourced from the scattering
  operator's **own** frame, which is what makes the stored moments equal
  :math:`S`'s internal projection term-for-term. The composite reports
  :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
  ``= False`` — no round-trip promise (its :math:`P`
  factor is a coisometry) — and accepts the driver's ``initial_guess``
  kwarg accepted-and-ignored (:meth:`WindowedSweep.apply
  <orpheus.sn.operators.windowing.WindowedSweep.apply>`; the multi-D walk
  has no bulk-seed consumer). Gated on ``sn_mesh.is_cartesian and
  sn_mesh.ndim == 2``.
* The **eigenvalue** inner
  (:meth:`SNSolver._solve_source_iteration <orpheus.sn.solver.SNSolver._solve_source_iteration>`)
  is reconstruction-free: it returns the scalar flux read off the
  :math:`\ell{=}0` moment block (``psi_typed.bulk.l_block(0)[0]``, the
  :math:`Y_0^0 = 1` identity), since the outer power iteration rebuilds
  from the scalar flux.
* The **fixed-source** inner
  (:func:`~orpheus.sn.solver.solve_sn_fixed_source`'s SI path) adds a
  **one-shot** full-angular reconstruction for ``Solution.angular_flux``
  — it re-evaluates the converged source ``q + Σ gains·ψ`` once through
  the *un-wrapped* base resolvent (``base_resolvent.solve(...)``).
  Fixed-source must return the full per-ordinate field, and because
  :math:`S`/:math:`B` consume the moments *equal* the full angular's
  moments (de-risk proven), the reconstructed source is the same and one
  sweep reproduces the converged iterate by the fixed point — so the
  reconstruction is bit-identical to the un-windowed converged
  :math:`\psi`.

Cross-references: the moment math :eq:`flux-moments` / :eq:`pn-scatter`
and the :math:`Y_\ell^m` convention live in :ref:`pn-scattering`
(:doc:`/theory/methods/sn/index`); the projection :math:`M`
(equation :math:`\phi_\ell^m = \sum_n w_n Y_\ell^m \psi_n`) and the
addition-theorem reconstruction :math:`R` are the spherical-harmonic
:class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` /
``reconstruction`` faces (:ref:`galerkin-projection`); the interior-face
cochain Phase 5a is orthogonal to is :ref:`wavefront-flux-cochain`; the
SI fixed point it reorganizes is the within-group inner of
:ref:`eigenvalue-posing`.


.. _sn-angular-windowing-in-sweep-accumulation:

Per-anti-diagonal moment accumulation — dropping the full-angular transient (Phase 5c)
--------------------------------------------------------------------------------------

Wave O step #205 **Phase 5c** (commit ``c7be111``, 2026-06-07) closes the
gap Phase 5a's :ref:`honest-scope warning
<sn-angular-windowing-honest-scope>` left open: the **per-sweep
full-angular transient**. Phase 5a held the *persistent* iterate as
moments (the 18.3× shrink) but still **materialized the full
per-ordinate angular field** :math:`\psi \in (N, n_g, n_x, n_y)` inside
*every* sweep, then projected it to moments **post-hoc** — a flat reduce
:eq:`angular-windowing-moment-projection` applied once at the end of the
sweep by the SH frame's analysis face ``frame.analysis.apply`` (the then
Phase-5a moment-windowed resolvent's ``solve`` called ``base.solve`` then
``self._projection.apply(full.bulk.values)``).
That transient full-angular array is the dominant peak-memory cost the
5a warning named.

Phase 5c moves the projection **into** the windowed anti-diagonal walk:
each topological level accumulates the harmonic moment tensor *directly*,
octant-by-octant into one shared global buffer, so the full angular
field is **never materialized** in the windowed iterate. The
``base.solve`` → flat-``apply`` post-projection **leaves production**; the
fused moment emit is what survives — reached at Phase 5c through the
resolvent's ``solve``, at #226 step 2 through the base resolvent's
``solve_moments``, and since #226 step 3 through the typed product's
:meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
that the :class:`~orpheus.numerics.iteration.SourceIteration` driver
applies directly (:ref:`windowing-retyped` and
:ref:`inverse-application-driver`, below).

.. admonition:: Key Facts (in-sweep moment accumulation)
   :class: tip

   - **The projection moves into the sweep.** Where 5a swept the full
     :math:`\psi` then flat-reduced :math:`\phi_\ell^m = \sum_n w_n
     Y_\ell^m \psi_n` once post-hoc, 5c accumulates the moments per
     anti-diagonal *during* the walk:
     :math:`\phi_\ell^m[\text{cells}_k] \mathrel{+}= \sum_{n\in
     \text{octant}} w_n Y_\ell^m(\hat\Omega_n)\,\psi_n[\text{cells}_k]`.
     The full :math:`(N, n_g, n_x, n_y)` field is never allocated for
     the windowed iterate.
   - **Measured peak-memory win: 3.06×.** On S8 / 4g / :math:`24\times24`
     the windowed ``solve`` peak drops 2.26 MB → 0.74 MB — the 1.47 MB
     full-angular transient is eliminated (moment tensor 0.111 MB). The
     win **grows with angular order**: the eliminated transient is
     :math:`N\,n_g\,n_x\,n_y`; the moment tensor is fixed at
     :math:`(L{+}1)(2L{+}1)\,n_g\,n_x\,n_y`.
   - **Principled-equivalence, NOT bit-identity.** The cross-octant
     ``+=`` reorders the ordinate sum vs the flat single-reduce; IEEE-754
     addition is non-associative, so the moments drift at ULP level.
     The per-cell :math:`w\,Y\,\psi` fold is **term-for-term identical**
     to the SH frame's analysis face ``frame.analysis.apply`` — only the
     accumulation *order* differs. Measured max-relative drift
     :math:`2.74\times10^{-16}` (:math:`\le 4N\varepsilon`, 4 ULP).
   - **The scalar is subsumed.** The scalar flux **is**
     :math:`\phi_0^0 = \texttt{moments}[0,0]` (:math:`Y_0^0 = 1`), read
     off the moment tensor — there is no separate scalar reduction
     (``coding-elegance`` Pattern 2: one source of truth).
   - **The fuller view is retained as a verification oracle.** The pre-5c
     "full-angular solve + flat ``frame.analysis.apply``" path is kept
     reachable and pins the optimized path
     (``feedback_aggressive_retirement`` — the "verification oracle"
     exception to retirement).
   - **2-D Cartesian only; the rest is untouched.** Gated on the genuine
     ``sn_mesh.is_cartesian and sn_mesh.ndim == 2`` (the C5.4 / #225
     sharpening of 5a's ``reduced is None`` proxy — the proxy was also
     true at 3-D Cartesian); the 1-D scan representation
     (:class:`~orpheus.sn.loss_representation.CumprodScan`) *raises*
     on a moment frame (illegal-states-
     unrepresentable). Curvilinear / Krylov stay full-angular; both
     one-shot ``Solution.angular_flux`` reconstructions stay separate
     full sweeps.


.. _sn-angular-windowing-in-sweep-transformation:

The transformation — post-hoc reduce → in-sweep accumulate
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5a's resolvent produced its moment iterate in **two arithmetic
stages**: (1) the wavefront walk materialized the full per-ordinate field
:math:`\psi \in (N, n_g, n_x, n_y)`, writing each anti-diagonal's
cell-average into the global angular buffer
(``angular_flux_octant[:, :, ii, jj] = psi_avg``); then (2) a flat
post-sweep reduce collapsed *all* :math:`N` ordinates at once,

.. math::

   \phi_\ell^m(\vec r)
   \;=\;\sum_{n=1}^{N} w_n\, Y_\ell^m(\hat\Omega_n)\,\psi_n(\vec r)
   \qquad\bigl(\texttt{frame.analysis.apply}, \;
   \texttt{einsum}\;\texttt{"n,nlm,n...->lm..."}\bigr),

which is :eq:`angular-windowing-moment-projection` again, evaluated once.
Stage (1)'s full-angular array is the transient that dominates the peak.

Phase 5c fuses the two stages. The :class:`_MovingFrontier
<orpheus.sn.loss_representation.sweep_graph._MovingFrontier>` walk
(:meth:`SweepDependencyGraph.walk_windowed
<orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed>` — at 5c the
solve-direction ``apply_windowed``, collapsed into the level-op walk at
S6.4(e)) already
visited every cell exactly once, anti-diagonal by anti-diagonal, with the
per-level cell-average :math:`\psi_n[\text{cells}_k]` in hand. Phase 5a
threw that local quantity into a global angular buffer and re-read it
later; 5c **projects it on the spot** and discards it. For each level
:math:`k` (anti-diagonal ``cells_k = (ii, jj)``) of each in-plane octant,

.. (vv-status rationale) The in-sweep per-anti-diagonal moment
   accumulation. The verifiable claim is the fuller-view-oracle
   equivalence: the in-sweep accumulation equals the flat post-sweep
   frame.analysis.apply of the same swept psi within the
   reduction-order drift bound (criterion 3 of bit-identity-vs-
   principled-equivalence). Pinned by
   test_2d_windowed_product_equals_post_projection (renamed from
   test_2d_windowed_moments_in_sweep_equal_post_projection when the
   solve_moments cross-reach retired into the typed product, #226 step 2),
   anchored to the structurally-independent SI≡Krylov-full scalar
   cross-check.
.. vv-status: harmonic-moment-projection implemented

.. math::
   :label: harmonic-moment-projection

   \phi_\ell^m[\text{cells}_k]
   \;\mathrel{+}=\;
   \sum_{n\in\text{octant}} w_n\, Y_\ell^m(\hat\Omega_n)\,
       \psi_n[\text{cells}_k],
   \qquad
   \texttt{moment\_buf[:, :, :, ii, jj]}
   \mathrel{+}=
   \texttt{einsum("nlm,ngd,n->lmgd",}\,
       Y_{\rm oct},\,\psi_{\rm avg},\,w_{\rm oct}\texttt{)},

accumulated into **one shared global** ``moment_buf`` of shape
:math:`(L{+}1, 2L{+}1, n_g, n_x, n_y)`. The four in-plane octants are
swept sequentially (the octant loop — since S6.4(b) the shared
``_OctantWalk.sweep_group`` frame, then ``sweep_octant_group`` — each octant
completing its full frontier walk), so the global buffer receives a **cross-octant
``+=``**: octant 1's complete contribution, then octant 2's, and so on.
The output branch (since S6.4(e) inside the ``_CellSolve`` level
operation; at 5c inside ``apply_windowed``) is a single
two-way switch at the per-level output site — *angular mode* writes
``angular_flux_octant[:, :, ii, jj] = psi_avg`` and accumulates the scalar
(``scalar_flux_buf``); *moment mode* does the :eq:`harmonic-moment-projection`
``+=`` instead. The frontier walk and the cell kernel are **untouched**,
which is exactly why the ``window ≡ full-field`` angular-mode oracle
(:ref:`Phase 5b <wavefront-flux-cochain>`) stays bit-identical: 5c adds a
branch only at the *consumer* of ``psi_avg``, never in its *production*.

The reduction-order contrast is the whole story:

.. list-table:: Post-hoc reduce (5a) vs in-sweep accumulate (5c)
   :header-rows: 1
   :widths: 22 39 39

   * -
     - **Phase 5a — post-hoc flat reduce**
     - **Phase 5c — in-sweep accumulation**
   * - Full :math:`\psi` materialized?
     - **Yes** — written to the global angular buffer, then re-read
     - **No** — projected per level and discarded
   * - Reduction shape
     - one flat ``einsum`` over all :math:`N` at once
       (``"n,nlm,n...->lm..."``)
     - per-level ``+=`` per octant
       (``"nlm,ngd,n->lmgd"``), cross-octant accumulate
   * - Per-cell arithmetic
     - :math:`\sum_n w_n Y_\ell^m \psi_n`
     - :math:`\sum_n w_n Y_\ell^m \psi_n` (**term-for-term identical**)
   * - FP reduction tree
     - one fixed order
     - reordered by the octant grouping ⟹ ULP drift
   * - Peak working set
     - moment iterate **+** full-angular transient
     - moment iterate **only**


.. _sn-angular-windowing-in-sweep-equivalence:

Why it is principled-equivalence, not bit-identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5a's per-step source was **bit-identical** to the full-angular arm
(0 ULP — :ref:`sn-angular-windowing-bit-identity`): the moments *stored*
equalled :math:`S`'s *own* internal projection of the same :math:`\psi`,
because the flat reduce and :math:`S`'s internal reduce shared the same
single reduction order. Phase 5c **necessarily** breaks that bit-identity
— and does so for a clean, documented reason rooted in the
``vv-principles`` § "Bit-identity vs principled-equivalence" three-criteria
test. The change is **accepted** because all three criteria hold:

#. **Principled at every step — the intermediate is named.** Each
   per-anti-diagonal partial moment :math:`\sum_{n\in\text{octant}} w_n
   Y_\ell^m(\hat\Omega_n)\,\psi_n[\text{cells}_k]` is **the octant's
   contribution to the harmonic moment** :math:`\phi_\ell^m` — a
   reactor-physics quantity (the per-direction-class flux moment), not a
   "whatever the reduction order produced" array. The cross-octant ``+=``
   composes named partial moments. This is exactly the
   unnamed-intermediate → named-intermediate move the criterion blesses
   (cf. the issue-#169 ``compute_keff`` worked example: per-group
   production rate is principled; the flat cell-product array is not).
#. **Verified against a structurally-independent reference — not
   old-vs-new ULP alone.** Old-vs-new distance is *necessary but not
   sufficient* (both arms could share a systematic offset). The
   structurally-independent ground is the **full-angular Krylov** scalar:
   ``test_2d_p1_aniso_moment_path_carries_signal_and_si_krylov_agree``
   cross-checks the windowed-SI moment :math:`\ell{=}0`
   (``scalar_flux()`` = :math:`\phi_0^0`) against the full-angular Krylov
   flux (Krylov is **never** windowed — it iterates the full bulk vector),
   and the closed-form homogeneous eigenvalue :math:`k_\infty =
   \nu\Sigma_f/\Sigma_a` (the closed-form pillar; MMS does **not** prove
   eigenvalues) anchors the eigenvalue reuse.
#. **The drift is FP-non-associativity, dimensionally explainable.** For
   a single-step computation the bound is :math:`(\text{reduction
   depth})\times\varepsilon`. The reduction depth is the ordinate count
   :math:`N`; the cross-octant ``+=`` reorders the **same** :math:`N`
   summands. The de-risk probe (``derivations/diagnostics/diag_p5c_moment_accum.py``,
   git-excluded, promoted into the permanent gate then deleted per the
   diagnostic-promotion policy) measured max-relative drift
   :math:`2.74\times10^{-16}` — about :math:`78\times` under the
   :math:`4N\varepsilon \approx 2.1\times10^{-14}` bound (:math:`N = 24`
   ordinates, 4× headroom for the partial-sum nesting within octants). A
   drift *above* :math:`4N\varepsilon` would signal an algorithmic change
   masquerading as FP noise — a wrong octant-:math:`Y` slice (Mode 2), a
   missing weight (Mode 3), or an :math:`\ell/m` index drift (Mode 5) —
   and the gate would catch it.

The scalar flux is **subsumed**, not separately reduced. Because
:math:`Y_0^0 = 1` (ORPHEUS unnormalized-harmonic convention,
:ref:`spherical-harmonics`), the :math:`\ell{=}0` slot **is** the scalar
flux, :math:`\phi_0^0 = \sum_n w_n \psi_n` (the :math:`Y_0^0 = 1` read of
:ref:`sn-angular-windowing-factoring`). The moment-mode sweep therefore
returns ``(moment_buf, None)`` — no separate ``scalar_flux_buf``
accumulation — and the eigenvalue inner reads
``moment_buf[0, 0]`` for the outer power iteration's scalar flux. The
windowed scalar is identical to the angular-mode scalar up to the same
reduction-order drift; the existing ``scalar_flux_buf`` einsum
``"ngd,n->gd"`` is literally the :math:`\ell{=}0` case of the moment
einsum ``"nlm,ngd,n->lmgd"`` with :math:`Y_0^0 = 1`.


.. _sn-angular-windowing-in-sweep-oracle:

The fuller-view verification oracle
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Phase 5c **relinquishes a fuller view** — the materialized full-angular
field that the post-hoc reduce consumed. By the
``feedback_aggressive_retirement`` "verification oracle" exception, an
optimization that gives up a fuller view of a concept **keeps that fuller
view reachable** as a verification oracle that pins the optimized path.
The pre-5c "full-angular ``base.solve`` then flat
``frame.analysis.apply``" path is exactly
that fuller view: it is *not* deleted, and — since #226 step 2 — it is
the **inherited**
:meth:`OperatorProduct.apply <orpheus.numerics.operator.OperatorProduct.apply>`
body of the windowed composition itself (``P.apply(A⁻¹.apply(rhs))``, the
un-fused factor-by-factor evaluation the
:class:`~orpheus.sn.operators.windowing.WindowedSweep` overrides). The
permanent gate ``test_2d_windowed_product_equals_post_projection``
asserts the fused in-sweep
:meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
result equals that deforested oracle of the same swept
:math:`\psi`, within the de-risk bound, over the **full moment tensor —
including the** :math:`\ell\ge 1` **block**.

The :math:`\ell\ge 1` coverage is the load-bearing reason this oracle
exists. The :math:`\ell{=}0` scalar cross-check (the SI ≡ Krylov-full
test, criterion 2 above) is **blind to** :math:`\ell\ge 1`: a windowing
that silently swapped or dropped a higher moment would converge the right
scalar flux while corrupting the anisotropic source. The oracle pins the
:math:`\ell\ge 1` block that the scalar cross-check cannot see — it is the
**Mode-5** (:math:`\ell/m` index drift) and **Mode-2** (wrong octant-:math:`Y`
slice) catcher.

.. warning::

   The oracle and the system-under-test share the **same** cell kernel
   (:meth:`~orpheus.transport.spatial.diamond.DiamondDifference.cell_kernel_batch`)
   and the **same** :math:`Y` / ``weights``, differing only in reduction
   *order*. They are therefore **procedurally** independent, **not
   structurally** independent (``vv-principles`` § "Structural
   independence"). The oracle gate alone is a same-math-rearranged
   comparison and is **not sufficient** on its own — it MUST be paired
   with the structurally-independent SI ≡ Krylov-full scalar anchor
   (criterion 2) and the closed-form :math:`k_\infty` eigenvalue. The
   oracle's *value-add* is precisely the :math:`\ell\ge 1` block the
   structurally-independent scalar anchor is blind to; its
   anti-degeneracy leg asserts :math:`\max|\phi_{\ell\ge1}| > 10^{-3}\,
   \max|\phi_0|` so the pinned higher-moment block is genuinely non-zero
   on the canonical config.

The metric the gate uses is the **scale-relative** drift
:math:`\max|\phi^{\rm SUT} - \phi^{\rm oracle}| / \max|\phi^{\rm oracle}|
\le 4N\varepsilon`, **not** an element-wise
``assert_array_almost_equal_nulp``. The moment tensor spans ~3 orders of
magnitude (the :math:`\ell{=}0` scalar dominates; the :math:`\ell\ge 1`
blocks are :math:`\sim 10^{-3}\times`), so an element-wise ULP comparison
would inflate a machine-:math:`\varepsilon` *absolute* difference on a
small :math:`\ell\ge 1` element into hundreds of ULP even when the reorder
is pure FP noise. The scale-relative drift is the principled-equivalence
quantity (criterion 3).


.. _sn-angular-windowing-in-sweep-numerical-evidence:

Numerical evidence
~~~~~~~~~~~~~~~~~~~

.. list-table:: Phase 5c gates (commit ``c7be111``)
   :header-rows: 1
   :widths: 34 18 48

   * - Gate
     - Result
     - What it pins
   * - **De-risk drift** (canonical config
       ``2d_2g_p1_aniso_dd_8x4_het_si``: 8×4 het, vacuum-x /
       reflective-y, mixture B :math:`\bar\mu = 0.6` P1, S4 :math:`N = 24`)
     - **2.74e-16** rel
     - in-sweep accumulation ≡ flat post-projection, :math:`\le 4N\varepsilon`
       (4 ULP) — criterion 3 (FP-reorder, not a bug)
   * - **Peak memory** (S8 / 4g / :math:`24\times24`,
       ``diag_p5c_peak_memory`` tracemalloc)
     - **3.06×** (2.26 MB → 0.74 MB)
     - the 1.47 MB full-angular transient eliminated; moment tensor 0.111 MB
   * - **Windowing L1**
       (``test_2d_anisotropic_windowing.py``, incl. the new equivalence gate)
     - **4 ✓**
     - P1 ≠ P0 + SI ≡ Krylov-full; full-recon self-consistency
       (:math:`\Sigma w\psi = \phi`); reflective trace ≠ 0; in-sweep ≡ oracle
   * - **Eigenvalue 2-D** (``test_keff_2d.py``)
     - **19 ✓**
     - closed-form :math:`k_\infty`, P1-changes-:math:`k_{\rm eff}`,
       SI ≡ Krylov, Jacobi ≡ G-S (Mode 9), refinement convergence
   * - **Regression snapshot**
       (``test_dd_regression.py``, ``2d_2g_p1_aniso_dd_8x4_het_si``)
     - within ``SAFETY × conv_tol``
     - scalar flux drift 6920 ULP / :math:`9.81\times10^{-13}` rel
       (~10× headroom under the :math:`1.0\times10^{-11}` gate)
   * - **5b oracles + matvec consumer**
       (``test_2d_full_field_oracle.py``, window≡full graph,
       ``TestAnisoMomentSourcePath``, A2D-1 hash pin)
     - bit-identical / intact
     - the angular-mode oracle and the moment *consumer* (``S.apply``)
       are untouched — 5c changes only moment *production*

The **peak-memory win grows with angular order** because the eliminated
transient scales as :math:`N\,n_g\,n_x\,n_y` while the moment tensor is
fixed at :math:`(L{+}1)(2L{+}1)\,n_g\,n_x\,n_y`: at S4 :math:`N = 24`,
:math:`L = 1` the moment tensor is :math:`6/24 = 1/4` the angular size;
at S8 :math:`N = 80` it is :math:`6/80 \approx 1/13`; a Lebedev-order-17
:math:`N = 110` it is :math:`6/110 \approx 1/18`. The transient is the
**linear** peak win Phase 5a's honest-scope warning deferred to "Phase
5b/5c": 5a delivered the persistent-iterate shrink, 5b (the
:ref:`moving-frontier interior cochain <wavefront-flux-cochain>`) cut the
interior face transient, and 5c cuts the full-angular *output* transient.


.. _sn-angular-windowing-in-sweep-honest-scope:

Honest scope — what 5c does and does NOT do
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning::

   Phase 5c eliminates the **per-sweep full-angular transient** of the
   *windowed SI* path. State the boundaries carefully:

   * **The persistent iterate was already moments (5a).** 5c does not
     shrink the held iterate — that was 5a's 18.3×. 5c shrinks the
     *transient* the held-moment iterate's sweep still materialized.
   * **Both** ``Solution.angular_flux`` **reconstructions stay full
     sweeps.** The eigenvalue final reconstruction and the fixed-source
     windowed one-shot reconstruction each re-run a *separate*
     full-angular sweep — the within-group resolvent ``solve``
     (:func:`~orpheus.sn.coupled_system.build_within_group_system`,
     applying its ``.resolvent``) — to return the user-facing
     :math:`(N, n_g, n_x, n_y)` field. They are **untouched** — the
     user-facing angular flux is bit-identical (the
     :math:`\Sigma w\psi = \phi` self-consistency gate and the step-3
     ``np.array_equal`` reconstruction pin both hold).
   * **Krylov, 1-D, and curvilinear stay full-angular.** Krylov iterates
     the full bulk vector (no moment sub-iterate); the curvilinear
     Morel–Montry Carlson coupled-pole seed reads the per-ordinate
     iterate at :math:`\mu = -1` (lesson L21), which the moment tensor
     cannot carry. Both are gated out by the genuine
     ``sn_mesh.is_cartesian and sn_mesh.ndim == 2`` (C5.4 / #225 — 3-D
     Cartesian is excluded too); the 1-D scan
     (:class:`~orpheus.sn.loss_representation.CumprodScan`) *raises* if a
     moment frame reaches it, so the unwindowable regime is
     unrepresentable, not merely unreached.

So the three-phase arc is complete: **5a** holds the persistent iterate
as moments (typed-state + persistent-storage win); **5b** carries the
interior face cochain on a rolling moving frontier (interior-transient
elimination + 3-D enabler); **5c** projects the swept flux to moments
in-sweep (full-angular *output*-transient elimination — the linear peak
win). Together they realize the asymptotic peak reduction the moment
iterate makes possible.


.. _sn-angular-windowing-in-sweep-implementation:

Implementation map (Phase 5c)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The moment OUTPUT mode threads as an **optional projection**
(``moment_frame=None``), mirroring the ``reflect=None``
dependency-injection idiom of :func:`_sweep_scheduled <orpheus.sn.loss_representation._sweep_scheduled>` — **not** a boolean flag.

.. note::

   **Retyping (#226 step 2).** At Phase 5c the two output modes were
   **named methods** on the resolvent surface — ``solve`` (full angular)
   vs ``solve_moments`` (moments). That public ``solve_moments`` — a
   method whose output-mode argument silently changed the operator's
   *codomain* — was a composition wearing a config, and it was retired.
   The moment emit is now the typed composition ``P @ A.inverse()``
   (:ref:`windowing-retyped`); the ``moment_frame`` kwarg survives only
   on the **private** ``_solve_timed_full_field`` body, the single
   application-context entry that the fused
   :meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
   drives. The map below is updated to that state.

* the solve-direction windowed walk (at 5c ``apply_windowed``; since
  S6.4(e) :meth:`SweepDependencyGraph.walk_windowed
  <orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph.walk_windowed>` × the
  ``_CellSolve`` level operation) — the single two-way branch at the
  per-level output site (angular write vs the
  :eq:`harmonic-moment-projection` ``moment_buf`` ``+=``). The
  frontier walk and cell kernel are untouched.  The apply-direction walk
  (at 5c ``residual_windowed``; now ``walk_windowed`` × ``_CellResidual``
  — the Krylov matvec) is untouched — Krylov stays full-angular.
* the per-group octant frame (``sweep_octant_group`` then, since S6.4(b),
  ``_OctantWalk.sweep_group``) /
  :func:`_sweep_scheduled <orpheus.sn.loss_representation._sweep_scheduled>` /
  :func:`_sweep_jacobi <orpheus.sn.loss_representation._sweep_jacobi>` — thread the
  optional ``moment_frame`` (2-D Cartesian only; the 1-D scan
  ``CumprodScan.sweep`` raises on a moment frame). Moment mode skips the
  per-octant angular allocation and the scheduled 2-D body returns
  ``(moment_buf, None)``.
* the private
  :meth:`InvertibleOperator._solve_timed_full_field
  <orpheus.sn.operators.streaming.InvertibleOperator._solve_timed_full_field>`
  body (duck-shared by
  :meth:`ScheduledInvertibleOperator._solve_timed_full_field <orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator._solve_timed_full_field>`)
  is the **single** application-context entry: given a ``moment_frame`` it
  emits harmonic moments, else the full angular field. ONE
  representation-sweep call per body for both output modes (since S6.5 the
  operator's own ``loss_representation.sweep`` — the same instance the
  matvec consumes); only the bulk **wrap** differs
  (:class:`~orpheus.transport.fields.angular_flux.AngularFlux` vs
  :class:`~orpheus.transport.fields.harmonic_moment_flux.HarmonicMomentFlux`).
  The former public ``solve_moments`` cross-reach — on
  ``InvertibleOperator`` **and** the dissolved ``_GaussSeidelResolvent``,
  plus the solver-side ``_solve_scheduled`` plumbing twin — retired into
  this one private body (#226 step 2, :ref:`windowing-retyped`).
* the fused entry is
  :meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
  (``P @ A.inverse()``), which calls that private body with its
  ``moment_frame``; the analysis factor
  :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator` carries
  the SH frame (whose basis
  :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` ``L`` carries
  the truncation order). Since #226 step 3 the
  :class:`~orpheus.numerics.iteration.SourceIteration` driver **applies**
  that fused ``apply`` directly (:ref:`inverse-application-driver`); the
  transitional ``_MomentWindowedResolvent`` ``.solve`` adapter that once
  forwarded to it is gone.

Cross-references: the post-hoc reduce 5c replaces is
:eq:`angular-windowing-moment-projection`; the :math:`Y_0^0 = 1`
scalar-flux read is :ref:`sn-angular-windowing-factoring`; the
moving-frontier interior cochain 5c rides is :ref:`wavefront-flux-cochain`
(Phase 5b); the persistent-iterate shrink 5c completes is
:ref:`sn-angular-windowing-honest-scope` (Phase 5a). The
principled-equivalence framework is ``vv-principles`` § "Bit-identity vs
principled-equivalence".


.. _windowing-retyped:

Windowing retyped — the moment emit as a composition (#226 step 2)
------------------------------------------------------------------

Phases 5a–5c built the moment-windowed path as a pair of **named
methods** on the resolvent surface: ``solve`` returned the full angular
field, ``solve_moments`` returned the harmonic moments. The taxonomy
step-2 carve (#226) retired ``solve_moments`` — because a *public method
whose output-mode argument silently changes the operator's codomain is a
composition wearing a config*, not a mode of one morphism.

The two output modes do not share a codomain: ``solve`` lands in the
full-angular composite carrier :math:`V = V_{\rm bulk} \oplus V_\partial`,
while ``solve_moments`` landed in the *moment-bulk* composite
:math:`V^{\rm mom} = \Phi_{\rm bulk} \oplus V_\partial`. A change of
codomain is, by the two-layer law of :ref:`operator-algebra` ("two
operators, one substrate — never two views, one operator"), a
**different arrow**; and an arrow whose target differs from what
:math:`A^{-1}` produces is a *composition* of :math:`A^{-1}` with a
second arrow, never a configuration flag on :math:`A^{-1}`. The honest
object is therefore

.. math::

   \text{windowed} \;=\; P \circ A^{-1},
   \qquad
   P \;=\; \underbrace{M_{\rm frame}}_{\text{analysis on the bulk}}
           \;\oplus\;
           \underbrace{\mathrm{Id}}_{\text{on the trace}},

where :math:`A^{-1}` is the within-group forward's inverse
(:class:`~orpheus.sn.operators.sweep_operator.SweepOperator` on
:math:`L+C`, or on the reified splitting matrix :math:`M` —
:ref:`si-gauss-seidel-reification`) and :math:`M_{\rm frame}` is the
scattering frame's :attr:`~orpheus.numerics.frame.FrameBase.analysis`
face, the angular→moment reduction :math:`\phi_\ell^m = \sum_n w_n
Y_\ell^m \psi_n` (:eq:`harmonic-moment-projection`). The boundary trace
passes through **un-reduced**: windowing is interior-bulk-only (the
reflective :math:`B` coupling reads the full per-ordinate face trace —
:ref:`sn-angular-windowing-geometry-restriction`), so :math:`P` is the
identity on the :math:`V_\partial` summand.

The coisometry factoring of the windowed contract
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`P` — the :class:`~orpheus.sn.operators.windowing.BulkAnalysisOperator`
— is a **block coisometry**, not an isomorphism. Its bulk face satisfies

.. math::

   \text{analysis} \circ \text{reconstruction} \;=\; 4\pi\,\mathrm{I}

under the no-prefactor SH convention (:ref:`spherical-harmonics`) — the
addition-theorem tight-frame identity, pinned by
``test_pi_R_is_4pi_identity_through_the_frame``. It is emphatically **not**
:math:`\mathrm{I}`: asserting :math:`\Pi R = \mathrm{I}` was the ERR-051
mistake — the coisometry carries the :math:`4\pi` frame constant, and a
test that hard-coded :math:`= \mathrm{I}` verified the *wrong* invariant.
In the other order :math:`\text{reconstruction} \circ \text{analysis}
\neq \mathrm{I}` (moments discard the ordinate-resolved angular content
above the truncation order :math:`L`), so :math:`P` is **not
invertible**. By the product's invertibility-closure law the composite
:math:`P \circ A^{-1}` therefore honestly reports ``is_invertible =
False`` (its ``P`` factor is structurally non-invertible), and makes *no
round-trip promise* — which is the whole point: the old ``solve_moments``
name suggested a *solve* (an inverse) where the object is a **projection
composed with an inverse**.

Fusion is an evaluation strategy, not a new operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The composition itself is
:class:`~orpheus.sn.operators.windowing.WindowedSweep` (an
:class:`~orpheus.numerics.operator.OperatorProduct` subclass), spelled by
the operator algebra ``P @ A.inverse()`` — the
:meth:`BulkAnalysisOperator.__matmul__ <orpheus.sn.operators.windowing.BulkAnalysisOperator.__matmul__>`
dispatch recognises a :class:`SweepOperator` right factor and returns the
fused product (mirroring the ``L + C`` fusion precedent, #261: the
dispatch is one-directional on the specific leaf).
:meth:`WindowedSweep.apply <orpheus.sn.operators.windowing.WindowedSweep.apply>`
overrides the inherited factor-by-factor body with the substrate's
MOMENT-emit mode (the per-anti-diagonal accumulation of :ref:`Phase 5c
<sn-angular-windowing-in-sweep-accumulation>`), so the
:math:`(N, n_g, n_x, n_y)` full-angular intermediate is **never
materialized** (the ~3× linear peak-memory win). The **inherited**
:meth:`OperatorProduct.apply <orpheus.numerics.operator.OperatorProduct.apply>`
body — ``P.apply(A⁻¹.apply(rhs))``, which *does* materialize the
intermediate — is retained verbatim as the **fuller-view verification
oracle** (the aggressive-retirement oracle exception): the two
evaluations differ only in the ordinate-sum reduction *order*, so they
agree to principled-equivalence within the scale-relative
:math:`4N\varepsilon` bound (measured :math:`1.8\times10^{-16}` on the
product SUT, re-measured from the pre-migration
:math:`2.7\times10^{-16}`). The permanent gate
``test_2d_windowed_product_equals_post_projection`` (renamed/migrated
from ``test_2d_windowed_moments_in_sweep_equal_post_projection`` when the
SUT became the product) is that fused-≡-deforested pin; it is a
*representation*-equivalence pin (procedurally, not structurally,
independent — shared kernel), whose structurally-independent anchors are
the SI ≡ Krylov-full scalar cross-check and the closed-form
:math:`k_\infty` eigenvalue.

.. admonition:: What was tried and rejected — the moment-proxy residual gate
   :class: warning

   An early proposal verified the windowed path with a **moment-proxy
   residual**: compute a residual *in moment space* and read its
   smallness as evidence that ``solve_moments`` inverted its operator. It
   was rejected as **category-confused**. :math:`P` is a coisometry, so
   :math:`P \circ A^{-1}` has no inverse and no round-trip to take a
   residual *of* — a residual gate presumes a solve, but the windowed
   object is a *projection composed with a solve*. Its only honest
   correctness statements are (i) representation-equivalence to the
   deforested oracle and (ii) the structurally-independent scalar anchor.
   Reifying the composition made the confusion **structural**:
   :class:`~orpheus.sn.operators.windowing.WindowedSweep` reports
   ``is_invertible = False`` (its coisometry ``P`` factor is
   non-invertible), so there is no round-trip an inverse-residual gate
   could measure.

At #226 step 2 one transitional wrapper remained: the driver still spoke
``.solve``, so a thin ``_MomentWindowedResolvent`` adapter held the
product and mapped its ``solve`` to the product's ``apply`` (the
``initial_guess`` kwarg accepted for the
:class:`~orpheus.numerics.iteration.SourceIteration` contract and
dropped). **Taxonomy step 3 removed even that** — the SI driver now
consumes inverse-application operators *directly*, so the production SI
holds ``P @ A.inverse()`` with no wrapper. That driver-consumption model
— how the solver builds the inverse and the driver applies it, why an
apply-only step operator is legitimate, and why the adapter could finally
dissolve — is :ref:`inverse-application-driver`.

Cross-references: the reduction :math:`M_{\rm frame}` is
:eq:`harmonic-moment-projection`; the frame's ``analysis`` /
``reconstruction`` faces are the
:class:`~orpheus.numerics.frame.GalerkinFrame`'s
(:ref:`galerkin-projection`); the two-layer "operators, not views" law is
:ref:`operator-algebra`; the reified :math:`M = (L+C-B_{\rm lower})` the
windowed forward may wrap is :ref:`si-gauss-seidel-reification`
(:doc:`/theory/methods/sn/index`); the invertibility-closure and
principled-equivalence framework is ``vv-principles`` § "Bit-identity vs
principled-equivalence".


.. _sn-iteration-primitives:

Iteration primitives (operator algebra)
========================================

Wave E Round 1 (Issue #163) lifts the iteration primitives out of
:class:`~orpheus.sn.solver.SNSolver` into stand-alone operator-algebra
consumers in :mod:`orpheus.numerics.iteration`.  They consume the
Wave A :class:`~orpheus.numerics.operator.LinearOperator` Protocol
triple :math:`(A, S, F)` directly — no transport-solver knowledge
beyond the operator contract.

The :math:`(A - S - F)\,\psi = q_{\rm ext}` framing
----------------------------------------------------

The Boltzmann transport equation in its operator-algebra form
factors into three pieces (Lewis & Miller 1984 §6.4; Trefethen &
Bau 1997 §3.2 frame the matrix-free Krylov view):

.. math::

    (A - S - F)\,\psi = q_{\rm ext}
    \qquad\text{(within-group fixed source)}

.. math::

    (A - S)\,\psi = \tfrac{1}{k}\,F\,\psi
    \qquad\text{(eigenvalue)}

where :math:`A = L + C = \Omega\cdot\nabla + \Sigma_t` is the
invertible loss (streaming-collision) operator (see
:ref:`sn-streaming-operator`),
:math:`S` is the scattering source operator
(see :class:`~orpheus.transport.operators.scattering.ScatteringOperator`),
:math:`F = \chi\otimes\nu\Sigma_f` is the rank-1-in-energy fission
emission operator (see :class:`~orpheus.transport.operators.fission.FissionOperator`),
and :math:`q_{\rm ext}` is an external source.

SourceIteration: discrete fixed-point realisation
--------------------------------------------------

:class:`~orpheus.numerics.iteration.SourceIteration` solves the
fixed-source equation by classical fixed-point iteration:

.. math::

    \psi_{n+1} \;=\; A^{-1}\,(S\,\psi_n + F\,\psi_n + q_{\rm ext}).

The convergence rate is bounded by the spectral radius
:math:`\rho(A^{-1}(S+F))`.  For an SN sweep applied to a homogeneous
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
via the :func:`~orpheus.sn.coupled_system.build_within_group_system`
single-source-of-truth builder):

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
:mod:`orpheus.sn.loss_representation.sweep_schedule`; since #226 step 2 it
feeds :meth:`SNBoundaryOperator.split <orpheus.sn.operators.boundary.SNBoundaryOperator.split>`,
which reifies the within-group forward as the **regular splitting matrix**
:math:`M = (L+C-B_{\rm lower})`
(:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`,
:ref:`si-gauss-seidel-reification` below — this replaces the dissolved
``_GaussSeidelResolvent``); the public entry is
:func:`~orpheus.sn.solver.solve_sn_fixed_source` via the
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
  :class:`~orpheus.transport.operators.scattering.ScatteringOperator`; higher
  Legendre orders add the :math:`P_\ell` blocks);
* :math:`B` is the **boundary reflection** gain — trace-only,
  realised by :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`,
  delivering :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}` on
  specular faces (see :ref:`bc-extraction` in
  :doc:`/theory/foundations/operator_algebra`);
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
:class:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule` without
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
       (:class:`~orpheus.sn.loss_representation.sweep_graph.OctantLabel`), in quadrature
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
:class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph`).

The selection lives in :func:`~orpheus.sn.solver._select_si_resolvent`:
``"gauss_seidel"`` on a multi-D Cartesian mesh returns
``((L+C) - parts.lower, (S, parts.upper))`` — the regular splitting
:math:`(L+C-B) = M - B_{\rm upper}`: the strictly-lower half
:math:`B_{\rm lower}` folds *into* the reified forward :math:`M`
(:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`,
whose ``solve`` is the octant-group forward substitution) while the
complement :math:`B_{\rm upper}` **lags as an ordinary external gain** —
structurally congruent with the Jacobi arm's ``(S, B)``, so the
:class:`~orpheus.numerics.iteration.SourceIteration` driver needs **no**
case split (the dissolved ``_GaussSeidelResolvent`` needed a bespoke
``(…, (S,))`` gain shape; the reification made the two arms uniform).
``"jacobi"`` (and any 1-D mesh) returns ``(L+C, (S, B))`` — both
:math:`S` and :math:`B` lagged (the degenerate :math:`B_{\rm lower}=0`).
In **both** cases :math:`S` is a lagged gain: the sweep never re-scatters
mid-sweep.  Only the boundary coupling gets the Gauss-Seidel treatment.

.. _si-gauss-seidel-reification:

The reified splitting matrix (#226 step 2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The #226 taxonomy step-2 carve replaced the duck-typed
``_GaussSeidelResolvent`` with an honest **reified splitting matrix**.
The dissolution matters because the old resolvent paired an ``apply``
with a ``solve`` of **different operators**: its ``apply`` computed
:math:`(L+C)\psi` while its ``solve`` inverted :math:`(L+C-B_{\rm
lower})`.  An operator whose forward and inverse faces disagree is not an
operator — and the disagreement is measurable: its round-trip defect was
**O(1)** (:math:`\lVert M^{-1}(M\psi)-\psi\rVert = 2.667`, the §17
falsifier-3 finding).

The boundary Gauss-Seidel is exactly a **regular matrix splitting** of the
within-group loss:

.. math::

   (L+C-B) \;=\; \underbrace{(L+C-B_{\rm lower})}_{M}
             \;-\; \underbrace{B_{\rm upper}}_{N},
   \qquad
   \psi_{k+1} \;=\; M^{-1}\bigl(q + S\,\psi_k + B_{\rm upper}\,\psi_k\bigr).

The reified :math:`M`
(:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`)
is an honest :class:`~orpheus.numerics.operator.OperatorSum` over the
leaves :math:`\{(L+C),\,-B_{\rm lower}\}`, so its ``apply`` is the leaf
sum :math:`(L+C)\psi - B_{\rm lower}\psi` and its ``inverse()`` returns a
genuine :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` on
:math:`M` — the two faces of **one** operator, the way :math:`A` and
:math:`A.H` are.  Because the scheduled walk (the same uniform
sweep-and-reflect loop as Jacobi, differing only in the
:class:`~orpheus.sn.loss_representation.sweep_schedule.SweepSchedule`) is
**exact forward substitution** for :math:`M`, ``M.inverse().apply`` now
round-trips ``M.apply`` at machine precision: the W2-round-trip gate
measures :math:`5.2\times10^{-16}` (bulk) / :math:`4.4\times10^{-16}`
(trace) — the O(1) defect gone.

The row-split law
^^^^^^^^^^^^^^^^^

:meth:`SNBoundaryOperator.split <orpheus.sn.operators.boundary.SNBoundaryOperator.split>`
returns a named :class:`~orpheus.sn.operators.boundary.BoundarySplit`
pair ``(lower, upper)`` of
:class:`~orpheus.sn.operators.boundary.SNMaskedBoundaryOperator` — each
the whole-trace :math:`B` masked to a per-face set of inflow ordinate
**rows**.  Which rows belong to which half is pure **schedule-order**
semantics, computed by
:meth:`SweepSchedule.lower_inflow_rows <orpheus.sn.loss_representation.sweep_schedule.SweepSchedule.lower_inflow_rows>`:

   An inflow row :math:`(f, m')` is in :math:`B_{\rm lower}` **iff**
   ordinate :math:`m'`'s octant group is swept *strictly after* face
   :math:`f`'s reflect group.

A reflective face :math:`f` is reflected exactly **once**, after its
**last** outflowing octant group (the ERR-056 fan-in rule above), at
which point every outflow feeding :math:`f` is complete.  A row swept
strictly after that reflect therefore reads the **fresh** current-iterate
reflection — realized in-sweep by the forward substitution
(:math:`B_{\rm lower}`); a row swept at-or-before the reflect keeps the
**lagged** seed (:math:`B_{\rm upper}`, the cyclic back-edges plus every
row of a never-reflected face — vacuum, white, albedo, periodic).  The
partition is **exact**: the specular map flips one direction-cosine sign,
so a row and its source always sit in *different* octants — :math:`B` has
no octant-diagonal, and :math:`B = B_{\rm lower} + B_{\rm upper}` is a
bit-exact per-face split (the W2-split gate).  The Jacobi schedule yields
an empty lower support (:math:`B_{\rm lower}=0`, :math:`B_{\rm upper}=B`)
— the degenerate that recovers the plain lagged-:math:`B` iteration.

The additive row-masked in-sweep reflect
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Solving :math:`Mz = y` on a strictly-lower inflow row is the inhomogeneous
forward-substitution row :math:`z_{\rm in} = y_{\rm row} + (Bz)_{\rm
row}`.  The buffer already holds the seed :math:`y_{\rm row}` (nothing
writes a lower row before its face's reflect), so the in-sweep reflect
**accumulates** the fresh reflection onto it — the ADDITIVE row-masked
verb
:meth:`SNMaskedBoundaryOperator.reflect_rows_inplace <orpheus.sn.operators.boundary.SNMaskedBoundaryOperator.reflect_rows_inplace>`
(``bf[rows] += (B·bf)[rows]``).  This is what makes ``M.inverse()`` exact
on an **arbitrary** rhs (not merely on production's zero-lower-row
subspace), and it leaves the upper (lagged) rows carrying the seed the
splitting :math:`\psi_{k+1}=M^{-1}(q+B_{\rm upper}\psi_k)` says they carry
— the returned trace **is** the splitting's honest iterate.

.. admonition:: What was tried and rejected — the whole-face-overwrite reflect
   :class: warning

   The dissolved resolvent used a whole-face **ASSIGNMENT**
   (:math:`\psi.{\rm inflow} \leftarrow B\cdot\psi.{\rm outflow}`, the
   verb now named
   :meth:`SNBoundaryOperator.reflect_inflow_inplace <orpheus.sn.operators.boundary.SNBoundaryOperator.reflect_inflow_inplace>`).
   As the in-sweep row update it is **wrong**: it *dropped*
   :math:`y_{\rm row}` and stamped fresh values onto rows the splitting
   defines as **lagged**.  It was benign in production only because a
   reflective inflow row's seed is zero there — but O(1)-wrong as a
   general inverse, which is precisely the round-trip defect the old
   pairing masked.  The whole-face assignment verb is **retained** (single
   source of truth via ``_reflect_trace``) for the two callers whose
   inflow is *wholly recomputed* each sweep and is not a solved unknown of
   a linear row: the direct fixed-source SI loop and the eigenvalue
   reconstruction sweep (:func:`~orpheus.sn.solver._reflect_outflow_into_inflow`
   delegates to it).

The source-subspace domain honesty note
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One subtlety is worth recording so a future reader does not mistake it for
a bug.  The sweep substrate re-derives the **outflow-definition** rows
(the walk's ``shed`` overwrites :math:`x.{\rm out}` with the streamed
value), so the scheduled walk realizes :math:`M^{-1}` **exactly on the
source subspace** :math:`\{y : y.{\rm outflow\text{-}rows}=0\}` — not on
the whole space.  This is not a limitation: every production rhs lands in
that subspace (:math:`q + S\psi + B_{\rm upper}\psi` all write bulk /
inflow rows only), and its :math:`M`-preimage is the set of
**trace-consistent** states :math:`x.{\rm out}={\rm streamed}(x.{\rm
bulk})` — i.e. actual transport states, which is exactly what a solve
output is.  It is the **same** property the already-landed
:math:`(L+C)`\ ``.solve`` has; the :math:`B` feedback of :math:`M` merely
makes it visible.  The W2-round-trip gate therefore round-trips a
trace-**consistent** state and asserts machine precision on **both** the
bulk and the trace — a *stronger* claim than the bulk-only falsifier, and
one the confused pairing still fails at O(1) (its ``apply`` lacks the
:math:`B_{\rm lower}` subtraction entirely).

Verification — the mutation redefinition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reification changed which mutations have teeth.  The spec's original
**M-SPLIT** mutation ("make the mask disagree with the in-sweep fold")
became **unrepresentable** by construction: the masked :math:`B_{\rm
lower}` single-sources the row split for **both** ``M.apply`` **and** the
in-sweep reflect, so there is no second site to make disagree.  It was
replaced by two mutations, one per gate:

* **M-SPLIT-DIR** — flip the split direction (upper-as-lower) in
  ``lower_inflow_rows``.  The flipped rows are read *before* their face's
  reflect fires, so the in-sweep fold never reaches the reader while
  ``apply`` subtracts it: the W2-round-trip defect returns to O(1)
  (``test_mutation_split_direction_reddens_round_trip``).  This is the
  Mode-1-family ``>`` vs ``<`` convention catcher.
* **M-SPLIT-PART** — doctor one half's rows after the split so the
  partition is no longer complementary: :math:`B \neq B_{\rm lower} +
  B_{\rm upper}` and the W2-split gate reddens
  (``test_mutation_partition_break_reddens_split``).

The converged fixed point is **splitting-invariant** (``vv-principles``
Mode 9): ``test_w2_fixed_point_equivalence_diagonal_cubature`` runs G-S ≡
Jacobi to solver tolerance on a config that **breaks** the degenerate
coincidences — a diagonal (``level_symmetric``) cubature with shared
faces (the ERR-056 regime; an axis-aligned ``product`` quad makes
octant-G-S accidentally exact) on a heterogeneous vacuum-x / reflective-y
box (anisotropic flux; the fully-reflective isotropic box is the Mode-9
degenerate).  That gate is **necessary but not sufficient** — the
load-bearing correctness gate is the W2-round-trip; the FP-invariance
pins only the *splitting* claim (same :math:`\psi^*`, only the rate
differs).  All gates live in ``tests/sn/solve/test_gauss_seidel_reification.py``
(``@pytest.mark.foundation`` — software invariants of the splitting, no
theory-page ``:label:`` and no ``verifies()``).

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
:meth:`ScatteringOperator.foldable_part <orpheus.transport.operators.scattering.ScatteringOperator.foldable_part>`
/ :meth:`residual_part <orpheus.transport.operators.scattering.ScatteringOperator.residual_part>`
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
<orpheus.sn.solution.IterationHistory.total_inner_iterations>`
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

    \psi_{n+1} \;=\; (A - S)^{-1}\,F\,\psi_n / k_n

.. math::

    k_{n+1} \;=\; \frac{\sum (F\,\psi_{n+1})}
                       {\sum (A\,\psi_{n+1}) - \sum (S\,\psi_{n+1})}

The dominance ratio :math:`|k_1/k_0|` governs outer-loop
convergence (Trefethen & Bau §27).  The inner solve uses
:class:`SourceIteration` with operator triple :math:`(A, S, 0)` —
the fission contribution at the inner level is the **external
source** :math:`F\psi_n/k_n`, NOT a within-group fixed-point term.
Every outer iteration warms up its inner :class:`SourceIteration`
from :math:`\psi_n` (the previous outer iterate); this is the same
amortisation pattern :class:`SNSolver` uses today.

The :math:`k`-update is the **hardwired** operator-form Rayleigh
quotient (:meth:`KEigenvalue.compute_keff
<orpheus.numerics.iteration.KEigenvalue.compute_keff>`)

.. math::

    k \;=\; \frac{\sum (F\,\psi)}{\sum (A\,\psi) - \sum (S\,\psi)}

— the operator-level spelling of the unified :math:`k` discipline,
fission production over net removal (see :ref:`sn-keff-estimator`).
When ``A`` carries streaming + collision,
:math:`\sum(A\psi) - \sum(S\psi)` is absorption + leakage − the
neutron-multiplying :math:`(n,2n)` emission (the in- and out-group
scatter cancel into :math:`\Sigma_a` via :math:`\Sigma_t - \Sigma_s`),
term-for-term the method-layer functional :eq:`sn-keff-update` with
the volume measure absorbed into the operators' action.  Because the
leakage rides **inside** :math:`A`, this spelling never had the #291
omission.

.. note:: **Estimator injection retired (R8, #259 P1, 2026-07-03).**

   Pre-R8, :class:`KEigenvalue` accepted ``keff_estimator`` /
   ``production_estimator`` *callables* (defaulting to module-level
   ``_default_*`` functions) so a caller could substitute its own
   field-to-scalar functional.  Those kwargs, the ``KeffEstimator`` /
   ``ProductionEstimator`` aliases, and the ``_default_*`` module
   functions are **gone**: the estimators are now hardwired methods
   (:meth:`~orpheus.numerics.iteration.KEigenvalue.compute_keff` /
   :meth:`~orpheus.numerics.iteration.KEigenvalue.compute_production_rate`),
   arithmetic bit-identical to the retired defaults.

   The seam was **dead by design, not dead by being unwired.**  The
   five method-layer solver families (SN / CP / diffusion / MoC /
   homogeneous) implement the
   :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol
   *directly* and never routed through ``KEigenvalue`` at all; and by
   the consistency theorem below, injection could only ever have
   introduced an *inconsistent* functional.  Removing the hook makes
   that illegal estimator/problem pairing unrepresentable
   (``coding-elegance`` Pattern 4 — illegal states unrepresentable).

   **Consistency theorem.**  The loop poses
   :math:`(A - S)\,\psi = F\psi/k` and converges to the fixed point
   :math:`\psi^\star` with
   :math:`(A-S)\,\psi^\star = F\psi^\star/k^\star`.  Applying the
   all-ones covector :math:`\mathbf 1^\top` (the ``\sum``) to both
   sides gives
   :math:`\sum(A\psi^\star) - \sum(S\psi^\star)
   = \sum(F\psi^\star)/k^\star`, so the hardwired ratio returns
   **exactly** :math:`k^\star`.  Every functional that agrees with the
   posed balance at the fixed point returns the same number; the
   "freedom" the injection seam advertised was illusory — all
   *consistent* choices collapse to one value, and any *different*
   injected estimator is by construction inconsistent with the problem
   the loop solves.

   **Honest-**\ ``A.apply``\ **contract.**  The theorem needs
   :math:`\sum(A\psi)` to be the true net-removal rate — i.e.
   ``A.apply`` must compute the real streaming + collision action, not
   a stub.  The injection seam had historically existed to paper over
   a *scalar-level test adapter* whose ``A.apply`` was dishonest; R8
   removes the paper and requires the honest action.  An adapter with
   a stubbed ``A.apply`` now yields a visible non-eigenvalue here — the
   correct failure, surfaced by design rather than masked by a
   substituted functional.

Relocating the :math:`k`-update into the algorithm itself — a Rayleigh
step on the resolvent :math:`A_{\rm loss}^{-1}M` so the loop is
*literally* K/α-agnostic — is the α-eigenvalue wave's first step; see
the :ref:`eigenvalue-posing` honest-scope note.

.. _choosing-inverse-realisation:

Choosing the :math:`A^{-1}` realisation (the inverse-operator family)
---------------------------------------------------------------------

.. note:: **Re-framed (2026-06-12, Issue #195).**

   This section's motivating premise — that the WDD sweep ``A.solve``
   has a *closure-bias-driven* fixed point distinct from the symmetric
   ``apply`` closure, so routing to Krylov-on-``apply`` is what "closes
   ERR-026" — was the *two-distinct-closures* picture.  ERR-058 (#195)
   showed the curvilinear wrong fixed point was the *closure-seed*
   family; once the seeds are fixed the sweep and the matvec are ONE
   discrete system and SI ``A.solve`` :math:`\equiv` Krylov-on-``apply``
   bit-identical (see :ref:`sn-err-058-closure-seed-closeout`).  The
   generality that hook named **survives** — the caller still chooses
   how :math:`A^{-1}` is realised (direct sweep, Green splitting, dense
   factorisation, or Krylov-on-apply preconditioned by the sweep)
   without re-implementing the driver — but since the #226 taxonomy
   **step 3** that choice is a **type** (which inverse-operator family
   member ``A.inverse()`` builds), not a caller-supplied ``inverter``
   callable, and it is no longer *required* to obtain the correct
   curvilinear fixed point.  Read the WDD-bias-vs-Krylov contrast below
   as Wave-E-era history.

Since the #226 taxonomy **step 3** neither iteration primitive takes an
``inverter`` callable.  The *solver* layer builds the inverse **once** —
as an operator, ``A.inverse()`` (or an explicit family member) — and the
drivers **apply** it:

* :class:`SourceIteration` receives the pre-inverted step operator
  ``A_inv`` directly and calls
  ``A_inv.apply(rhs, initial_guess=psi_prev)`` each step — the inverse
  family's canonical seeded-apply signature
  (:class:`~orpheus.numerics.iteration.SupportsSeededApply`).  The former
  caller-supplied ``inverter`` callable — and its per-type
  ``inspect.signature`` dispatch on a duck-typed ``solve`` — dissolved
  into that one contract.
* :class:`~orpheus.numerics.iteration.KrylovAcceleration`'s ``inverter``
  parameter was **renamed** ``preconditioner``: a GMRES *left
  preconditioner* :math:`M \approx (A - S - F)^{-1}` approximates the
  FULL within-group system inverse, a different object from the
  iteration's :math:`A^{-1}` step (the old name was a category mistake).
  Its default is ``A.inverse().apply`` (the sweep) when ``A`` is
  invertible.

The load-bearing principle of the Wave-E reconciliation of ERR-026 —
*the caller controls how* :math:`A^{-1}` *is realised* — therefore
**survives**, but as a **type choice, not a callable injection**: the
caller picks *which inverse-operator-family member* realises
:math:`A^{-1}`, and the driver never changes.  This is the inverse axis
of the :ref:`three-layer operator surface <capability-set-semantics>`
(runtime predicate
:attr:`~orpheus.numerics.operator.LinearOperator.is_invertible` →
operator-returning method ``inverse()`` → realisation verb ``solve``):
inverting *with* an operator IS applying its inverse **object**,
``A.inverse().apply(b)``.

The realisations are **sibling KINDS of** :math:`A^{-1}`, keyed by the
operand's mathematical structure (historically selected by the
``inverter`` argument; now by the type the caller constructs, or that
``A.inverse()`` returns):

* :class:`~orpheus.sn.operators.sweep_operator.SweepOperator` — the
  **direct triangular WDD sweep** of the streaming–collision composite
  :math:`(L+C)`.  Its :meth:`apply` delegates to
  :meth:`~orpheus.sn.operators.streaming.InvertibleOperator.solve`
  (bit-identical), so this IS the historical ``inverter = None`` /
  ``A.solve`` default — the WDD asymmetric closure — now reached as
  ``(L+C).inverse()``.
* :class:`~orpheus.numerics.green_operator.GreenOperator` — the
  **preconditioned-splitting** inverse of a general
  :class:`~orpheus.numerics.operator.OperatorSum` (the Neumann /
  multiple-scattering series :math:`\sum_k (A^{-1}B)^{k} A^{-1}`, wrapping
  the SAME :class:`SourceIteration` driver as its application engine).
  See :ref:`green-operator` (:doc:`/theory/foundations/operator_inverse_family`).
* :class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`
  — the **dense direct LU factorisation** (materialise :math:`[A]`,
  factor once, back-solve per :meth:`apply`).  The natural realisation
  for the synthetic L0 case (``A`` a small dense matrix) and any small
  exactly-solvable block.  See :ref:`matrix-inverse-operator`
  (:doc:`/theory/foundations/operator_inverse_family`).
* :class:`~orpheus.numerics.operator.InverseOperator` — the **generic
  solve-backed wrapper** for a value-bearing leaf (a diagonal or
  cross-section multiplier), whose :meth:`apply` delegates to the leaf's
  own ``solve``.
* :meth:`A.inverse() <orpheus.numerics.operator.OperatorSum.inverse>` —
  the **structure-keyed factory** that lets the operand's structure
  choose for you: a general sum returns the
  :class:`~orpheus.numerics.green_operator.GreenOperator`, while the
  sweep-invertible :math:`(L+C)`
  :class:`~orpheus.sn.operators.streaming.InvertibleOperator` subclass
  **shadows by MRO** and returns the direct
  :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`.

Read as inverse KINDS, the old WDD-bias-vs-Krylov contrast is the
choice between two of these siblings on curvilinear SN:

* The **sweep-inverse** (:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`,
  the WDD asymmetric closure) has a closure-bias-driven self-consistent
  fixed point on curvilinear meshes that — *before the ERR-058 seed
  fix* — was **not** the fine-mesh-limit transport solution (ERR-026).
* The **Krylov/Green-inverse** solves the symmetric :meth:`apply`
  closure of
  :class:`~orpheus.sn.operators.streaming.InvertibleOperator` with the
  sweep injected only as the GMRES ``preconditioner`` :math:`M` — so the
  converged solution comes from the closure that agrees with analytical
  references, while the sweep merely accelerates the iteration.  This
  was the Wave-E Round 2 path that closed ERR-026 for curvilinear SN.

Since ERR-058 (#195) the two are no longer a *correctness* fork: once
the closure **seeds** are fixed, the sweep and the matvec are ONE
discrete system and SI :math:`\equiv` Krylov-on-:meth:`apply`
bit-identical, so the sweep-inverse and the Krylov/Green-inverse
converge to the same fixed point.  The remaining distinction is one of
**rate** (the standard transport-Krylov win as :math:`c \to 1`;
Adams & Larsen 2002).

Because every family member conforms to the ONE seeded-apply contract
(:class:`~orpheus.numerics.iteration.SupportsSeededApply`), the driver
is **inversion-strategy-agnostic**: the same :class:`SourceIteration`
runs the synthetic L0 case (``A`` a dense matrix →
:class:`~orpheus.numerics.matrix_inverse_operator.MatrixInverseOperator`),
the L1 SN case (``(L+C).inverse()`` → the WDD
:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`), and the
Wave-E Krylov-on-:meth:`apply` case — with **no re-implementation**,
because the inversion strategy now rides the *type* of ``A_inv``, not a
branch inside the driver.

Operand requirements
--------------------

The primitive constructors enforce their operands' surface at
construction time, NEVER mid-iteration (the same Wave A philosophy that
gates :class:`~orpheus.numerics.operator.OperatorSum` etc.).  Since the
#226 taxonomy carve the requirement is stated on the operator surface
directly (predicate + verb), not a capability tag:

* :class:`SourceIteration` consumes a **pre-inverted** step operator
  ``A_inv``, which MUST provide a callable ``apply`` — the step operator
  arrives already inverted, so an apply-only object is legitimate by
  design (#226 taxonomy step 3).  A missing ``apply`` raises
  ``TypeError`` at construction.
* :class:`KEigenvalue` poses the ``(A, S, F)`` triple: ``A`` MUST report
  :attr:`~orpheus.numerics.operator.LinearOperator.is_invertible`
  ``= True`` — the posing layer builds :math:`A^{-1}` via
  ``_seeded_inverse(A)``, and a non-invertible ``A`` raises
  :class:`~orpheus.numerics.operator.NotInvertible`.  ``S`` / ``F`` MUST
  provide ``apply`` (pass
  :class:`~orpheus.numerics.operator.ZeroOperator` for the
  scattering- / fission-free case).

Constructor failure raises ``TypeError`` (a missing ``apply``) or
:class:`~orpheus.numerics.operator.NotInvertible` (a non-invertible
``A`` where the posing layer needs its inverse), naming the operand at
fault.

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

* :class:`~orpheus.sn.operators.streaming.InvertibleOperator` — the
  :math:`L` operand the SN solver ships, with both ``apply``
  (symmetric closure) and ``solve`` (WDD asymmetric closure).
  See :ref:`sn-streaming-operator` for the design rationale.
* :class:`~orpheus.transport.operators.scattering.ScatteringOperator` — the
  :math:`S` operand carrying P\ :sub:`ℓ` scattering plus (n,2n).
* :class:`~orpheus.transport.operators.fission.FissionOperator` — the rank-1-in-
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
:math:`(A, S, F)`-triple resolvent (SN, MoC) *and* the
**monolithic-matrix resolvent** (CP, diffusion, homogeneous) that has
no separable :math:`(A-S)^{-1}` factor.  All five solver families are
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
  :func:`~orpheus.sn.coupled_system.build_within_group_system` SSoT (see
  :ref:`sn-solver-operator-algebra-coordinator`).
* **CP** (collision-probability) — drives ``power_iteration`` through
  its own ``EigenvalueSolver``-Protocol implementation; its resolvent
  is **one BiCGSTAB on a monolithic collision-probability matrix**,
  which has no :math:`(A, S, F)` split.  This is exactly the family
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
natural :math:`(A, S, F)` triple and want to skip writing a full
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
  Krylov-on-apply with sweep preconditioning (the ``preconditioner``
  hook's use case — the ERR-026-closing curvilinear Krylov path).


.. _sn-solver-operator-algebra-coordinator:

SNSolver as an operator-algebra coordinator
============================================

:class:`~orpheus.sn.solver.SNSolver` consumes the operator triple
:math:`(A, S, F)` directly.  At construction time, the solver builds:

* :attr:`SNSolver.L` — :class:`InvertibleOperator` carrying the
  symmetric-closure streaming-collision operator.
* :attr:`SNSolver.S` — :class:`~orpheus.transport.operators.scattering.ScatteringOperator`
  carrying the P0 + (n,2n) + Pℓ Galerkin reconstruction (Wave D
  Issue 13).
* :attr:`SNSolver.F` — :class:`~orpheus.transport.operators.fission.FissionOperator`
  carrying the rank-1-in-energy fission emission (Wave D Issue 13).

Each of these is a :class:`~orpheus.numerics.operator.LinearOperator`
in the Wave A operator-algebra sense: predicate-typed, composable
under :class:`~orpheus.numerics.operator.OperatorSum` and
:class:`~orpheus.numerics.operator.OperatorProduct`, and protocol-
conforming so the iteration primitives in
:mod:`orpheus.numerics.iteration` consume them without SN-specific
plumbing.  The within-group inner solve is built once from a single
source of truth — the :func:`~orpheus.sn.coupled_system.build_within_group_system`
builder assembles the :class:`~orpheus.sn.coupled_system.WithinGroupSystem`
record, the honest within-group decomposition :math:`(L+C,\ S,\ B)` as a
named splitting :math:`A = M - N`: the invertible resolvent
:math:`M = (L+C)` plus its two lagged coupling gains :math:`N = (S,\ B_a)`
(the bulk scattering :math:`S` and the trace boundary reflection :math:`B`;
zero within-group fission), handed to the **variadic** driver
:math:`\text{Driver}(L_{\rm resolvent},\,*\text{gains})` (Wave O step
O.2a — the transitional :math:`S + B` fold is retired; see
:ref:`bc-extraction-variadic-driver` in :doc:`/theory/foundations/boundary_conditions`).
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
:class:`SourceIteration` from the :func:`~orpheus.sn.coupled_system.build_within_group_system` SSoT and
runs it; :meth:`SNSolver._solve_krylov` constructs a
:class:`KrylovAcceleration` from :func:`_within_group_krylov` and runs
that.  The Layer-3 resolvent of the SN row in the
:ref:`eigenvalue-posing` architecture is exactly these primitive
instances.

The primitive is **type-agnostic and angular-capable**: it operates on
the typed :class:`~orpheus.transport.timed_full_field.TimedFullField`
composite, which carries the full angular flux on its bulk.  Pℓ
anisotropic scattering therefore rides the angular bulk with no special
plumbing — :meth:`ScatteringOperator.apply` on the timeless
:class:`~orpheus.transport.full_field.FullField` operator carrier (the
driver's :class:`~orpheus.transport.timed_full_field.TimedFullField` iterate
reaches it via MRO) reads the angular moments off the composite and builds
the anisotropic source via :meth:`ScatteringOperator.build_aniso_source`,
all inside the primitive's normal RHS path.  There is **no scalar-flux
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

Beyond driving the within-group inner solve, the :math:`(A, S, F)`
framing organises the solver's outer API surface:

* :meth:`SNSolver.compute_fission_source` returns
  :math:`F\,\phi/k` — a thin delegator to :meth:`F.apply` with the
  :math:`1/k` outer-loop scaling applied at the solver level.
* :meth:`SNSolver.solve_fixed_source` solves
  :math:`(A - S)\,\psi = q_{\rm ext}` (with :math:`q_{\rm ext}` the
  fission source built by ``compute_fission_source``).  Two paths:

  * ``inner_solver="source_iteration"`` — sweep-driven fixed-point
    iteration; :math:`A^{-1}` is the WDD asymmetric sweep.  ERR-026-
    affected for curvilinear vacuum-BC cases.
  * ``inner_solver="krylov"`` — GMRES on :meth:`A.apply` with the
    sweep as preconditioner.  Reflective-BC equation map only
    (Round 3 owns the vacuum-BC extension).

* :meth:`SNSolver.compute_keff` returns **fission production over net
  removal**, :eq:`sn-keff-update` — the volume-weighted method-layer
  functional :math:`R_{\nu\Sigma_f}(\phi) / (R_{\Sigma_a}(\phi) + L -
  E_{2n}(\phi))`, derived in :ref:`sn-keff-estimator` below.  The
  SN-specific volume weighting lives here (in the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  fields) — one honest realization of the same discipline the
  operator-form :meth:`KEigenvalue.compute_keff` spells with the
  measure absorbed into the operators' action.  (Pre-#291 this method
  returned the leakage-blind :math:`\sum F\phi V / \sum \Sigma_a\phi V`
  ratio; see :ref:`sn-keff-estimator` for why that was a
  non-eigenvalue on any vacuum-bounded problem.)

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
outer loop updates :math:`k` from the net-removal balance
:eq:`sn-keff-update` (fission production over absorption + leakage −
:math:`(n,2n)` emission), with an inner loop that solves the
within-group scattering problem.

.. _sn-keff-estimator:

The reported eigenvalue: fission production over net removal
------------------------------------------------------------

:meth:`~orpheus.sn.solver.SNSolver.compute_keff` reports the eigenvalue
of the problem the inner solve **actually poses**.  This is the SN
symptom (#291) and the MoC/CP/homogeneous root (#259) of a single
discipline: *the reported* :math:`k` *must be the eigenvalue of the
fixed-source map every method scales only fission by* :math:`1/k`
*through* — scattering and the :math:`(n,2n)` emission are plain gains
assembled **inside** :meth:`~orpheus.sn.solver.SNSolver.solve_fixed_source`,
never scaled by :math:`1/k`.  An estimator that disagrees with its own
posed problem converges cleanly and silently to a **non-eigenvalue
ratio**.

.. math::
   :label: sn-keff-update

   k \;=\; \frac{R_{\nu\Sigma_f}(\phi)}
                {R_{\Sigma_a}(\phi) \;+\; L \;-\; E_{2n}(\phi)}

.. (vv-status rationale) Governing/definitional identity: the reported k
   IS the eigenvalue of the posed fixed-source map, not a solver
   eigenvalue-correctness claim against an external analytical reference
   (that rests on the multi-group heterogeneous L1/L2 references
   elsewhere on this page). The verifiable content is the cross-engine
   consistency gate tests/sn/eigenvalue/test_keff_estimator_gate.py
   (reported k == the converged fixed-point map ratio k* = P(Mφ*)/P(φ*),
   map-ratio ground-truth noise ≤ 2e-11) with in-file mutation teeth.
   Wiring @pytest.mark.verifies("sn-keff-update") onto that gate is a
   test-side follow-up (this docs pass could not touch tests/).
.. vv-status: sn-keff-update documented

The three terms are typed volume-integrated reaction-rate functionals
and one boundary functional:

* **Numerator** :math:`R_{\nu\Sigma_f}(\phi) = \int_V \nu\Sigma_f\,\phi\,dV`
  — the fission production, the typed
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  over :math:`\nu\Sigma_f` (the :math:`\phi^\dagger\!=\!1` degenerate of
  the homogenization Petrov–Galerkin bilinear).  The :math:`(n,2n)`
  emission is **not** production here — contrast
  :meth:`~orpheus.sn.solver.SNSolver.compute_production_rate`, the
  ERR-052 renormalisation scale anchor, which keeps *total* physical
  production (fission **plus** the :math:`(n,2n)` emission).  The role
  split is the load-bearing #259 correction: the same physical
  :math:`(n,2n)` neutrons are a **scale** quantity for the outer
  renormalisation but a **removal-reduction** in the eigenvalue balance.
* **Absorption** :math:`R_{\Sigma_a}(\phi) = \int_V \Sigma_a\,\phi\,dV`,
  with :math:`\Sigma_a = \Sigma_f + \Sigma_c + \Sigma_L +
  \sum_{g'}\Sigma_{2,g\to g'}` — i.e. ``absorption_xs`` counts the
  :math:`(n,2n)` **collision once** (the neutron is removed from its
  incident group by the collision).  See
  :attr:`~orpheus.data.macro_xs.mixture.Mixture.absorption_xs`.
* **Leakage** :math:`L` — the net vacuum-boundary outflow (below).  On a
  reflective (lattice) problem it is a **structural zero**.
* **Emission** :math:`E_{2n}(\phi) = \int_V \sum_{g,g'} 2\,\Sigma_{2,g'\to
  g}\,\phi_{g'}\,dV` — the :math:`(n,2n)` **emission** (two neutrons out
  per collision; the factor 2).  A gain, so it **reduces** net removal.

The net :math:`(n,2n)` effect on removal is therefore
:math:`\underbrace{\sum_{g'}\Sigma_{2,g\to g'}}_{\text{collision, in }\Sigma_a}
- \underbrace{2\Sigma_2}_{E_{2n}} = -\Sigma_2` — **one extra neutron
gained** per collision, exactly the physics of a neutron-doubling
reaction.

**Balance identity (divergence-telescoping).**  The angle- and
group-summed discrete cell balance for cell :math:`i` in the posed
eigenproblem is

.. math::

   \underbrace{\sum_{f\in\partial i}\!\bigl(\textstyle\sum_g J_g\bigr)\,\Delta A_f}
              _{\text{net face flow}}
   \;+\; \Sigma_{t,i}\,\phi_i\,V_i
   \;=\; \frac{1}{k}\,R_{\nu\Sigma_f,i}
        \;+\; \Sigma_{s,i}\,\phi_i\,V_i
        \;+\; E_{2n,i}

(streaming + total collision on the left; the isotropic fission source
scaled by :math:`1/k`, plus the *unscaled* scatter and :math:`(n,2n)`
gains, on the right).  Summing over **all** cells, every interior face
is shared by two cells with opposite outward normals and equal current
(continuity), so the interior face-flow terms **telescope to zero** —
only the domain-boundary faces survive, and their sum is the net
leakage :math:`L`.  With :math:`\Sigma_t - \Sigma_s = \Sigma_a` this
collapses to

.. math::

   \frac{R_{\nu\Sigma_f}(\phi)}{k} \;=\; R_{\Sigma_a}(\phi) \;+\; L
                                        \;-\; E_{2n}(\phi),

which is :eq:`sn-keff-update` rearranged.  This is the same discrete
divergence discipline the diffusion page states as
:math:`\mathbf 1^{\mathsf T}(C-S)=\Sigma_a` with interior-face
telescoping (see :ref:`diffusion-leakage-boundary-leaves`); SN and
diffusion report the *same* balance-law eigenvalue, differing only in
how the streaming operator is discretised.

The leakage functional
~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

   L \;=\; \sum_{f\,\in\,\text{vacuum}} \oint_{f} dA\,
           \sum_g J_g(\mathbf{r}_f)\,,
   \qquad
   J_g \;=\; \sum_m (\Omega_m\cdot\hat n_f)\, w_m\, \psi_{m,g}

is the face-area integral of the boundary trace's **net outward
current**.  The angular-to-scalar reduction :math:`J_g` is
:meth:`AngularBoundaryFlux.net_current
<orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux.net_current>`
— the single source of the :math:`\Omega\cdot\hat n\,w` contraction, the
angular sibling of the scalar trace's :math:`J = J^+ - J^-`
(:meth:`ScalarBoundaryFlux.net_current
<orpheus.transport.fields.scalar_boundary_flux.ScalarBoundaryFlux.net_current>`).
It is spelled through the trace space's own atoms — the signed
projection table
:attr:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.omega_dot_n`
and the :math:`|\Omega\cdot\hat n|\odot w` partial-current metric (using
the identity :math:`\operatorname{sign}(\Omega\cdot\hat n)\cdot
|\Omega\cdot\hat n|\,w = \Omega\cdot\hat n\,w`) — so no consumer
re-derives the cosine weighting.  Tangential ordinates carry zero weight
and drop out.

The face measure :math:`dA` is supplied by
:meth:`SNSolver._face_area_of`, matching the cell
``volume_measure`` exactly so the balance identity closes:

.. list-table:: Boundary-face measure by geometry
   :header-rows: 1
   :widths: 30 30 40

   * - Geometry
     - Face measure :math:`\Delta A`
     - Source
   * - 1-D slab
     - :math:`1` (per unit cross-section)
     - :attr:`MaterialMesh.areas <orpheus.transport.mesh.material_mesh.MaterialMesh.areas>`
   * - 1-D cylinder
     - :math:`2\pi R` (per unit height)
     - ``MaterialMesh.areas``
   * - 1-D sphere
     - :math:`4\pi R^2`
     - ``MaterialMesh.areas``
   * - 2-D Cartesian
     - transverse edge-cell widths (unit depth)
     - ``mesh.axes`` transverse extent
   * - 3-D Cartesian
     - :math:`\Delta A_{\mathbf c} = \prod_{j\ne a}\Delta_j[c_j]`
       (transverse-area outer product)
     - ``mesh.axes`` transverse extents

The :math:`d \ge 2` Cartesian arms are ONE generic body: the outer
product of the *other* axes' edge widths in **ascending axis order** —
the same codimension-1 enumeration as
:func:`~orpheus.transport.mesh.axis.face_shape`, so the measure array
broadcasts cell-for-cell against the ``(ng, *face_spatial)`` net
current, and the 2-D width vector is just the single-transverse-axis
degenerate (bit-identical to the pre-3-D spelling).

The 3-D arm originally shipped as a **typed refusal**
(``NotImplementedError``): guessing the transverse product's cell
ordering would silently mis-weight the leakage sum, and Cardinal Rule 1
forbids returning a wrong-but-plausible number.  The wire landed
2026-07-13 when the first 3-D vacuum eigenvalue consumer arrived (the
d=3 Mode-9 G-S≡Jacobi gate), with the ordering pinned twice in
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: an **object-level
pin** (face measure ≡ the boundary layer's ``volumes / Δ_axis``, the
mesh's own ascending-axis enumeration — vv Mode-12 discipline: pin the
object, not only the k functional) and the **k* map-ratio gate** on a
Mode-2 asymmetric all-vacuum box, whose teeth are proven by permanent
in-process mutants — a reversed transverse enumeration moves the
reported k by a measured **13.9 %** against the estimator-independent
:math:`k^*` (clean agreement :math:`6\times10^{-10}`), and a transposed
enumeration crash-REDs on the broadcast.  A trace carrying a
``#251`` transverse face-moment tail is refused loudly at the
consumption site (the face integral must consume ONLY the
transverse-average moment — higher Legendre face moments integrate to
zero over each face cell — and that slot-0 read has no consumer yet).

Reflective faces are a structural zero
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The reflective law equates a face's inflow to its reflected outflow
**exactly**, so the net current there vanishes *by construction*.
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` therefore
**skips** reflective faces — it never accumulates them, rather than
accumulating a value that ought to be zero but carries
:math:`\pm`-cancelling angular-sum floating-point noise.

This is a deliberate design choice with a bit-level payoff.  On an
all-reflective (lattice) problem :math:`L` is a structural ``0.0``, and
on a :math:`\Sigma_2`-free mixture :math:`E_{2n}` is exactly ``0.0`` (the
per-material :math:`(n,2n)` loop adds nothing), so
:eq:`sn-keff-update` reduces **bit-identically** to the historical
lattice functional ``production / absorption``.  Every pre-existing
reflective eigenvalue anchor is preserved to the last ULP — the
unification adds terms that vanish structurally, not numerically, on the
cases it must not perturb.

The scale bridge: trace of the last inner solve
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The leakage term reads the **typed** boundary trace of the last inner
solve (``self._psi_typed.boundary``), whereas the numerator/denominator
reaction rates consume the bare-array flux :math:`\phi` the estimator is
handed.  These two representations can be at **different scales**:
:func:`~orpheus.numerics.eigenvalue.power_iteration` renormalises
:math:`\phi` to unit production rate *between* the inner solve and the
:math:`k`-update (ERR-052), so the stored trace belongs to the
**un-renormalised** last iterate while the estimator's :math:`\phi` is
its renormalised multiple.

Leakage is degree-1 homogeneous in :math:`\psi`, so the fix is a single
rescale by the fission-production ratio of the two fluxes
(``self._phi_of_trace``, stored alongside the trace at **both**
inner-path returns) — exactly ``1.0`` when the caller passes the
returned flux itself.  The **contract** is therefore: the flux handed to
:meth:`~orpheus.sn.solver.SNSolver.compute_keff` must be a scalar
multiple of the last inner solve's flux (true for ``power_iteration`` and
for every manual solve-then-estimate loop).

If a vacuum face exists but **no** inner solve has stored a trace,
:meth:`~orpheus.sn.solver.SNSolver._boundary_leakage_rate` raises
``RuntimeError`` — the leakage cannot be answered honestly, and silently
returning it as zero would *reproduce the #291 omission*.  Fail loud;
never return a non-eigenvalue.

The R7 :math:`(n,2n)` convention fork
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The historical spelling put the :math:`(n,2n)` emission in the
**numerator** as production,

.. math::

   k_{\text{old}} \;=\; \frac{R_{\nu\Sigma_f} + E_{2n}}{R_{\Sigma_a}},

which is a **non-eigenvalue** of the posed map whenever
:math:`\Sigma_2 \neq 0` *and* :math:`k \neq 1`.  The reason is exactly
the posing asymmetry: the inner solve scales **only** fission by
:math:`1/k`; the :math:`(n,2n)` emission is an *unscaled* gain in the
sweep source.  So the eigenvalue of that map is
:math:`k^\star = R_{\nu\Sigma_f}/(R_{\Sigma_a} - E_{2n})` (reflective,
:math:`L=0`), and putting the *unscaled* emission in the numerator does
not recover it.  Writing :math:`f = R_{\nu\Sigma_f}`,
:math:`a = R_{\Sigma_a}`, :math:`e = E_{2n} = 2s_2` and substituting
:math:`f = k^\star (a - e)`:

.. math::

   k_{\text{old}}
   \;=\; \frac{k^\star (a - e) + e}{a}
   \;=\; k^\star \;+\; \frac{2 s_2\,(1 - k^\star)}{a}.

The two agree only when :math:`s_2 = 0` (no :math:`(n,2n)`) or
:math:`k^\star = 1` (critical).  For a supercritical
:math:`k^\star > 1` the correction is negative
(:math:`k_{\text{old}} < k^\star`); for subcritical, positive.  The MoC
and CP pages carry the same fork
(:eq:`moc-keff-update`, :eq:`cp-keff-update`); CP was the one member
already spelled on net removal.

What was tried and found
~~~~~~~~~~~~~~~~~~~~~~~~~~

The #291 bias was characterised pre-fix (commit ``d1daaac``, Gauss–
Legendre :math:`n=8`, map-ratio ground truth noise :math:`\le 2\times
10^{-11}`) across the five gate configurations:

.. list-table:: Pre-fix reported :math:`k` vs the posed-problem eigenvalue :math:`k^\star`
   :header-rows: 1
   :widths: 40 18 18 24

   * - Configuration
     - Pre-fix reported
     - Posed :math:`k^\star`
     - Bias
   * - homog. 2G vacuum slab (width 8)
     - 1.83767525
     - 0.98163269
     - :math:`+87.2\%` (:math:`L/A = 0.872`)
   * - het. vacuum sphere P\ :sub:`0`
     - 0.86484694
     - 0.70601977
     - :math:`+22.5\%`
   * - het. vacuum sphere P\ :sub:`1`
     - 0.85080423
     - 0.67876772
     - :math:`+25.3\%`
   * - reflective control (:math:`\Sigma_2=0`)
     - 1.87500000
     - 1.87500000
     - :math:`\equiv` (bias :math:`1.2\times 10^{-10}`)
   * - reflective :math:`\Sigma_2\neq 0`
     - 1.92857143
     - 2.61278195
     - :math:`-26.2\%` (R7 defect)

The two failure classes are visible in one table: the vacuum rows are
the **leakage omission** (#291) — the reported :math:`k` overshoots by
the leakage-to-absorption ratio :math:`L/A`; the last row is the **R7
:math:`(n,2n)` convention** — zero leakage, yet a
:math:`-26.2\%` error because the emission was mis-posed as production.
The reflective-control row is exactly the bit-identity guarantee above.
The exact check on the R7 row is
:math:`0.78/(0.5185 - 0.2200) = 2.61278`, and
:math:`(0.78 + 0.2200)/0.5185 = 1.92857` reproduces the old value —
matching :math:`k_{\text{old}} = k^\star + 2s_2(1-k^\star)/a` term for
term.

Post-fix, reported :math:`k` and the map-ratio :math:`k^\star` agree to
:math:`\le 6\times 10^{-10}` on all five.  The P\ :sub:`0`\ –P\ :sub:`1`
sphere gap :math:`\Delta` roughly **doubled** (:math:`1.404\times
10^{-2} \to 2.725\times 10^{-2}`) but stays inside the diagnostic
:math:`(10^{-3}, 5\times 10^{-2})` band — the P\ :sub:`1` anisotropic
correction is now measured against the correct eigenvalue on both
solves.

The V&V decision was a **principled re-baseline** (per ``vv-principles``
bit-identity-vs-principled-equivalence): the old reported :math:`k` was a
*different functional* from the posed problem's eigenvalue, so the new
value is not a regression to be tolerance-matched but a correction to be
verified against a structurally-independent reference (the fixed-point
map ratio).

Verification
~~~~~~~~~~~~

The permanent gate is
``tests/sn/eigenvalue/test_keff_estimator_gate.py``: it asserts the
reported :math:`k` equals the converged fixed-point map ratio
:math:`k^\star = P(M\phi^\star)/P(\phi^\star)` across the four physics
regimes — {vacuum slab, vacuum sphere (pinning the :math:`4\pi R^2`
face-area convention), reflective bitwise-degenerate, reflective
:math:`\Sigma_2\neq 0`} — with **in-file mutation teeth**: a
leakage-drop mutation reds the vacuum legs while staying bitwise-green
on reflective; a leakage sign-flip crash-reds through the scale-bridge
guard; and the old :math:`(n,2n)`-in-numerator convention reds the
:math:`\Sigma_2\neq 0` leg.

This is a **consistency** gate: the map ratio is the structurally-
independent ground truth for "does the estimator return the eigenvalue
of its own posed map", and is blind by construction to *which*
eigenvalue that is.  The SN solver's eigenvalue **correctness** — that
the posed map's eigenvalue is the *physically right* :math:`k` — rests
on the multi-group heterogeneous L1/L2 references elsewhere on this
page, not on this gate.

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

- Operator: :math:`A = \mu \nabla + \Sigt{}` (finite-difference
  gradients, symmetric closure carried by
  :meth:`InvertibleOperator.apply`)
- Solution variable: angular flux :math:`\psi(x, y, n, g)` (much
  larger than scalar flux)
- System: :math:`A\psi = b` where :math:`b` = fission + scattering
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


.. _sn-consuming-the-frame:

Consuming the frame in SN
=========================

Spatial homogenisation and energy condensation are **discrete-frame
projections** — the Petrov-Galerkin coefficient extraction
:math:`G^{-1}M` of a flux- (or spectrum-) weighted frame. All of that
theory — rate preservation, the source-group / sink-sum matrix rules, the
metric-fold-vs-bilinear adjoint argument, fractional-overlap re-binning,
the condense/homogenize asymmetry law, and the verification gates — is
the frame page's headline **Petrov-Galerkin** consumer; see
:ref:`sn-spatial-homogenization` and :ref:`sn-energy-condensation`
(:doc:`/theory/foundations/frame`). This section keeps only the **SN-layer orchestration**:
how the SN :class:`~orpheus.sn.solution.Solution` drives that machinery
from a converged flux.

Homogenisation: the solve → homogenize → re-solve loop
------------------------------------------------------

:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
takes a coarse mesh (:class:`~orpheus.geometry.mesh.Mesh1D` or
:class:`~orpheus.geometry.mesh.Mesh2D`) and returns a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` — the coarse
geometry already carrying one freshly-homogenised effective
:class:`~orpheus.data.macro_xs.mixture.Mixture` per coarse cell. The SN
:class:`~orpheus.sn.solution.Solution` owns the converged flux, so the SN
layer is what builds the flux-weighted **test** basis the frame consumes;
the frame itself, and the rate-preservation theory that *forces* the flux
weighting (rather than a plain volume average), live in
:ref:`sn-spatial-homogenization` (:doc:`/theory/foundations/frame`). The returned
``MaterialMesh`` is re-promoted to a solvable phase space by
:meth:`SNMesh.from_material_mesh
<orpheus.sn.mesh.augmented_mesh.SNMesh.from_material_mesh>`, closing the
**solve → homogenize → re-solve** loop. The return type is
**mesh-coupled** (geometry and materials born together) — the space half
of the condense/homogenize asymmetry law
(:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`).

Condensation: per-material representative spectra
-------------------------------------------------

:meth:`Solution.condense <orpheus.sn.solution.Solution.condense>` is the
SN-layer orchestration of energy condensation. It condenses **each
material with its own representative spectrum** — the flux·volume-weighted
flux over the cells where the material appears:

.. math::
   :label: energy-condensation-representative-spectrum

   \varphi^{(m)}_g \;=\;
   \sum_{i:\,\mathrm{mat}(i)=m} V_i\,\phi_{i,g},

.. (vv-status rationale) Representational identity: the per-material
   representative spectrum used as the condense test weight — the
   flux·volume-weighted flux over the material's cells (mirrors how
   ``homogenize`` derives its flux weight). A definition consumed by
   :meth:`Mixture.condense`; the end-to-end rate preservation it feeds is
   the L1 gate, not a separate claim.
.. vv-status: energy-condensation-representative-spectrum documented

used as the test weight in :meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>` — the data-layer
collapse verb, whose spectrum-weighted-collapse theory is
:ref:`sn-energy-condensation` (:doc:`/theory/foundations/frame`) — mirroring how
:meth:`Solution.homogenize` derives its flux weight from the same solved
flux. The result is a **portable** ``dict[int, Mixture]`` keyed by
material id — few-group cross sections carrying the coarse ``eg``, not
bound to any mesh (the **mesh-decoupled** half of the asymmetry law,
:ref:`sn-condense-homogenize-asymmetry`, :doc:`/theory/foundations/frame`). A material with
no flux in a fine group contributes zero weight there; the condense
frame's Moore–Penrose Gram handles any empty coarse group.


.. _sn-gotchas:

Gotchas
=======

Each gotcha is a **consequence → how it manifests → which test / level
catches it** — a trap that hides a solver bug behind a green test.  They
should *shrink* over time as the code hardens.

Degeneracy traps — a passing test that proves nothing
-----------------------------------------------------

.. admonition:: Gotcha — 1-group eigenvalue tests are degenerate
   :class: warning

   :math:`k = \nu\Sigma_f / \Sigma_a` is **flux-shape independent**: a
   1-group eigenvalue is a material-property ratio computable *without*
   solving the transport equation, so it cannot detect any error in the
   spatial, angular, or scattering operators.  **Any verification claim
   needs** :math:`\geq 2` **groups** (``vv-principles`` anti-pattern #3).  A
   1-group eigenvalue is still fine for a *rate* or *convergence-order*
   claim — declare the claim layer.

.. _sn-homogeneous-degeneracy-gotcha:

.. admonition:: Gotcha — homogeneous / uniform-rescale invariance hides coefficient bugs
   :class: warning

   Any eigenvalue problem whose target is the flux :math:`\phi` is invariant
   under a uniform rescale :math:`\phi \to C\phi` (the factor :math:`C`
   cancels in the Rayleigh quotient
   :math:`k = \nu\Sigma_f\,\phi / \Sigma_a\,\phi`).  Homogeneous and
   same-material multi-region problems have a **spatially-uniform** rescale,
   so they are blind to factor-of-two coefficient errors that preserve the
   flux shape — and, in curvilinear geometry, blind to redistribution bugs
   (flat flux → the :math:`\alpha` terms vanish identically).  Only a
   genuine **material interface** makes the rescale factor :math:`C(x)`
   position-dependent and breaks the cancellation.

   This is exactly how **ERR-025** hid.  A missing :math:`1/W` normalisation
   in the 1-D diamond-difference recurrence halved the per-ordinate flux, but
   for Gauss–Legendre :math:`W = \sum_n w_n = 2` the missing factor rescaled
   it back — so every homogeneous test passed at machine precision while the
   heterogeneous eigenvalue was :math:`\sim 1.5\,\%` wrong and did **not**
   converge away under mesh or :math:`S_N`-order refinement (the gap
   plateaued in angle).

   **Catcher:** at least one *absolute*-:math:`\phi` test — the fixed-source
   flat-flux diagnostic (:math:`Q/\Sigma_t`), or an absolute eigenvalue
   comparison against a structurally-independent heterogeneous reference.
   The live pins are the L0 symbolic-recurrence check
   :func:`tests.sn.sweep.slab.test_dd_recurrence.test_dd_per_cell_recurrence_matches_symbolic_derivation`
   and the L1 heterogeneous absolute-:math:`k` regression
   :func:`tests.sn.eigenvalue.test_keff_slab.test_heterogeneous_absolute_keff`
   (a 2-region A+B reflective slab pinned against the Case
   singular-eigenfunction reference; the pre-fix :math:`1.48\times10^{-2}`
   error fails it by two orders of magnitude).

.. admonition:: Gotcha — conservation holds even with wrong per-ordinate balance
   :class: warning

   Global particle balance **telescopes** by construction, so a *scalar*
   balance sum can hold to machine precision while the *per-ordinate*
   flat-flux residual is wrong (``vv-principles`` anti-pattern #8; the
   identity :math:`\sum_n w_n(\alpha_{n+1/2} - \alpha_{n-1/2}) = 0`
   annihilates per-ordinate redistribution errors that cancel in the sum).
   The load-bearing invariant is the **per-ordinate** flat-flux residual
   (= 0), not the telescoped scalar balance.

Curvilinear redistribution
--------------------------

.. admonition:: Gotcha — the α redistribution dome must stay non-negative
   :class: warning

   The curvilinear :math:`\alpha` dome must be non-negative; a negative entry
   drives NaN / overflow through the angular sweep.  The fixed-source
   flat-flux diagnostic (:math:`Q/\Sigma_t`) is the single most powerful
   curvilinear bug detector — a spike at :math:`r = 0` localises a missing
   :math:`\Delta A / w` geometry factor.

Solver-coordination traps
-------------------------

* **Renormalise-then-report ordering.**
  :func:`~orpheus.numerics.eigenvalue.power_iteration` renormalises
  :math:`\phi` to unit production **between**
  :meth:`~orpheus.sn.solver.SNSolver.solve_fixed_source` and
  :meth:`~orpheus.sn.solver.SNSolver.compute_keff`.  So
  ``compute_keff`` sees the *renormalised* :math:`\phi`, while the stored
  ``_psi_typed.boundary`` is the *un-renormalised* trace — the scale
  bridge above is what makes the leakage term consistent across that
  boundary.  Reordering the two (report before renormalise) would break
  the bridge's ``1.0`` shortcut.
* **The outer iterate must stay a bare** ``np.ndarray``.  The Mode-11
  live-arm sentinel in
  ``tests/sn/operators/test_fission_kernel_crosscheck.py`` proves that
  ``power_iteration`` feeds a **bare** :class:`numpy.ndarray` flux to
  :meth:`~orpheus.sn.solver.SNSolver.compute_fission_source`, so the
  bare-``np.ndarray`` dispatch arm of
  :meth:`FissionOperator.apply
  <orpheus.transport.operators.fission.FissionOperator>` is the *live
  production arm* (the sentinel wraps that registered leaf in-process and
  asserts the counter advances).  The estimator's
  :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
  evaluations read the same bare array.  Routing the outer iterate
  through a typed carrier would dark the arm (redding the sentinel) and
  break the estimator's evaluate path — the bare-array outer iterate is
  a load-bearing contract, not an implementation accident.

.. seealso::

   **Sweep-machinery gotchas** — the Krylov ``restart`` sizing bug
   (ERR-053 family), the product-cylinder edge-extrapolation data-flow
   invariant, and the Mode-12 / ERR-067 :math:`G`-reciprocity metric catch
   — are documented alongside the sweep at :ref:`sn-282-gotchas`.


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

.. [Rahnema2008] F. Rahnema, S. Douglass, and B. Forget,
   "Generalized Energy Condensation Theory,"
   *Nuclear Science and Engineering*, 160:41--58, 2008.
   DOI `10.13182/NSE160-41 <https://doi.org/10.13182/NSE160-41>`_.
   Expands the within-coarse-group flux in a set of orthogonal functions;
   the **zeroth moment** (the piecewise-constant basis function on the
   coarse group) recovers the standard flux-weighted multigroup average
   exactly, and the higher moments add the within-group spectral detail.
   Cited at :ref:`sn-condensation-petrov-galerkin-frame` as the rank-0
   precedent for the spectrum-weighted collapse (rank-:math:`>0` faithful
   reconstruction is deferred, `Issue #275
   <https://github.com/deOliveira-R/ORPHEUS/issues/275>`_).

.. [WIMSD] International Atomic Energy Agency,
   *WIMS-D Library Update*, IAEA-TECDOC / STI/Pub/1264, IAEA, Vienna,
   2007. Tables 11.1 (69-group) and 11.2 (172-group) energy-group
   structures, and Table 11.3 (the 172→69 correspondence). The coarse
   group structures ORPHEUS condenses onto; Table 11.3 is the
   derivation-validation oracle for the containing-interval partition.
   Cited at :ref:`sn-condensation-fractional-overlap`.

.. _sn-chapters:

Chapters in this sub-book
=========================

This page is the S\ :sub:`N` sub-book's index. It currently carries the
bulk of the method's theory inline; the decomposition into chapters is
tracked as issue `#231 <https://github.com/deOliveira-R/ORPHEUS/issues/231>`_.
Chapters split out so far:

.. toctree::
   :maxdepth: 2

   angular_quadrature
   loss_representation
   boundary_conditions
   verification
   history
