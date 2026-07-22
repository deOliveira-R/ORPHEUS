.. _theory-collision-probability:

====================================
Collision Probability Method
====================================

.. contents:: Contents
   :local:
   :depth: 3


.. _key-facts-cp:

Key Facts
=========

**Read this before modifying the CP solver.**

- CP works with **scalar flux** :math:`\phi`, not :term:`angular flux` :math:`\psi`
- :math:`P_{ij}` = probability that neutron born in *j* has first collision in *i*
- Convention: ``P[i,j]`` = birth in *j*, collision in *i*; flux update uses :math:`P^T` (see :ref:`scattering-matrix-convention`)
- Row sums = 1 (conservation), reciprocity: :math:`\Sigma_{t,i} V_i P(i,j) = \Sigma_{t,j} V_j P(j,i)`
- **Slab**: :math:`E_3` exponential-integral kernel (``scipy.special.expn``)
- **Cylinder**: :math:`\text{Ki}_4` Bickley-Naylor kernel (20,000-point lookup table)
- **Sphere**: exponential kernel
- White-BC approximation: isotropic return at cell boundary → ~1% gap vs reflective SN
- Inner iteration IS needed: self-scatter :math:`\Sigma_s(g \to g) \cdot \phi_g` makes source depend on solution
- **Gotcha**: tautological residual if checking ``denom * phi - transported`` (ERR-016)
- Power iteration tolerance ``keff_tol=1e-7`` bounds eigenvalue error to ~1e-6
- Gauss-Seidel update uses latest flux immediately; Jacobi uses previous iteration
- Verification uses :ref:`synthetic cross sections <synthetic-xs-library>`, not real nuclear data

.. admonition:: Conventions

   - Scattering matrix: :ref:`scattering-matrix-convention` — ``SigS[g_from, g_to]``, source uses transpose
   - Multi-group balance: :eq:`mg-balance` in :ref:`theory-homogeneous`
   - Cross sections: :ref:`theory-cross-section-data`
   - Verification: :ref:`synthetic-xs-library` — regions A/B/C/D
   - Eigenvalue: :ref:`power-iteration-algorithm` shared with all deterministic solvers


Overview
========

The collision probability (CP) method solves the
:ref:`multi-group eigenvalue problem <mg-eigenvalue-problem>` in integral
form.  Rather than tracking the angular flux
:math:`\psi(\mathbf{r}, \hat{\Omega}, E)` as in :term:`discrete ordinates <ordinate>` or
method of characteristics, the CP method works directly with the
**scalar flux** :math:`\phi(\mathbf{r}, E)` by integrating out the
angular variable analytically.

The key quantity is the **collision probability** :math:`P_{ij}`: the
probability that a neutron born uniformly and isotropically in region
:math:`i` has its *first* collision in region :math:`j`
:cite:`Stamm1983`.  Once the
:math:`P_{ij}` matrix is known, the transport problem reduces to a
matrix equation in the region-averaged scalar fluxes.

Three geometries are supported:

- **Slab** (1D Cartesian) --- the :math:`E_3` exponential-integral kernel
- **Concentric cylinders** (1D radial) --- the canonical
  :math:`\mathrm{Ki}_3` Bickley--Naylor kernel (A&S 11.2), evaluated via
  :func:`orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` (Chebyshev
  interpolant of :math:`e^{\tau}\,\mathrm{Ki}_3(\tau)` built from
  :func:`~orpheus.derivations.common.kernels.ki_n_mp` at 30 dps)
- **Concentric spheres** (1D radial) --- the exponential kernel with
  :math:`y`-weighted quadrature

All three share a single eigenvalue solver through the :class:`CPMesh`
augmented-geometry pattern.

**Derivation sources.**  The analytical eigenvalues and CP matrices used
for verification are computed independently by the derivation scripts.
These are the **source of truth** for all equations in this chapter:

- ``orpheus/derivations/continuous/flat_source_cp/geometry.py`` --- the unified geometry-dispatching
  core :class:`~orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry` and
  :func:`~orpheus.derivations.continuous.flat_source_cp.geometry.build_cp_matrix` (Phase B.2 refactor,
  commits ``f1b869b`` → ``bf128d3``)
- ``orpheus/derivations/continuous/flat_source_cp/slab.py`` --- slab :math:`E_3` kernel; thin facade
  over :data:`~orpheus.derivations.continuous.flat_source_cp.geometry.SLAB`
- ``orpheus/derivations/continuous/flat_source_cp/cylinder.py`` --- cylindrical canonical
  :math:`\mathrm{Ki}_3` kernel; thin facade over
  :data:`~orpheus.derivations.continuous.flat_source_cp.geometry.CYLINDER_1D`
- ``orpheus/derivations/continuous/flat_source_cp/sphere.py`` --- spherical :math:`e^{-\tau}` kernel;
  thin facade over :data:`~orpheus.derivations.continuous.flat_source_cp.geometry.SPHERE_1D`
- ``orpheus/derivations/common/kernels.py`` --- :math:`E_3` via
  :func:`~orpheus.derivations.common.kernels.e3_vec` (wraps
  :func:`scipy.special.expn`); arbitrary-precision :math:`\mathrm{Ki}_n`
  via :func:`~orpheus.derivations.common.kernels.ki_n_mp` (wraps
  :func:`mpmath.quad`). Double-precision :math:`\mathrm{Ki}_3`
  goes through :func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` — a
  Chebyshev interpolant of :math:`e^{\tau}\,\mathrm{Ki}_3(\tau)`
  built from ``ki_n_mp`` at 30 dps (~:math:`5\times 10^{-6}`
  accuracy). The legacy ``BickleyTables`` tabulation was retired
  in Phase B.4 (commit ``6badbe5``, Issue #94)
- ``orpheus/derivations/common/eigenvalue.py`` --- shared eigenvalue computation via
  :func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp` and
  :func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous`

Every equation in this chapter can be verified against these scripts.
Every numerical value cited was produced by them.


Architecture: Base Geometry and Augmented Geometry
===================================================

The CP solver separates geometry description from solver logic through
two layers:

1. **Base geometry** --- :class:`~geometry.mesh.Mesh1D` stores cell edges,
   material IDs, and the coordinate system.  It computes volumes and
   face areas via :func:`~geometry.coord.compute_volumes_1d` and
   :func:`~geometry.coord.compute_areas_1d`.

2. **Augmented geometry** --- :class:`CPMesh` wraps a ``Mesh1D`` and adds
   the CP-specific kernel, :term:`quadrature`, and
   :meth:`CPMesh.compute_pinf_group`.  The kernel is selected
   automatically from the mesh's coordinate system via a ``match``
   statement in :meth:`CPMesh.__init__`.

3. **Solver** --- :func:`solve_cp` creates a ``CPMesh``, builds the
   :math:`P^{\infty}` matrices for all energy groups, and runs power
   iteration via :class:`CPSolver` (satisfying the
   :class:`~numerics.eigenvalue.EigenvalueSolver` protocol).

.. code-block:: text

   Mesh1D (edges, mat_ids, coord)
       |
       v
   CPMesh (kernel + quadrature + compute_pinf_group)
       |
       v
   solve_cp() -> CPResult

**Design rationale (CP-20260404-001).**  Adding a new geometry requires
only a new ``_setup_*()`` method and kernel function --- the eigenvalue
solver, post-processing, and plotting are geometry-agnostic.  The
alternative (three separate solver classes) would duplicate ~200 lines of
iteration logic.  The derivation modules mirror this: all three
``_*_cp_matrix`` functions have identical white-BC closure code, differing
only in the kernel and quadrature.


The Integral Transport Equation
================================

Starting point: the steady-state, one-speed transport equation in
integral form.  For a neutron born at :math:`\mathbf{r}'` travelling
toward :math:`\mathbf{r}`, the **uncollided flux** is

.. math::
   :label: first-flight-kernel

   \phi^{\text{unc}}(\mathbf{r})
   = \int_{V} \frac{e^{-\tau(\mathbf{r}', \mathbf{r})}}{4\pi |\mathbf{r} - \mathbf{r}'|^2}
     \, Q(\mathbf{r}') \, dV'

where :math:`Q(\mathbf{r}')` is the total source (fission + scattering +
external) and :math:`\tau(\mathbf{r}', \mathbf{r})` is the **optical
path** (number of mean free paths) between the two points:

.. math::
   :label: optical-path

   \tau(\mathbf{r}', \mathbf{r})
   = \int_0^{|\mathbf{r} - \mathbf{r}'|}
     \Sigt{}\bigl(\mathbf{r}' + s\,\hat{\Omega}\bigr) \, ds


.. _flat-source-approximation-cp:

Flat-Source Approximation
-------------------------

The CP method assumes the source is **spatially flat** within each
sub-region :math:`i`.  The MOC solver uses the same approximation; see
:ref:`flat-source-approximation-moc`.

.. math::
   :label: flat-source

   Q(\mathbf{r}) = Q_i \quad \text{for } \mathbf{r} \in V_i

Under this approximation, the collision rate in region :math:`j` due to
sources in region :math:`i` can be written as

.. math::
   :label: collision-rate

   \Sigt{j} \, \phi_j \, V_j
   = \sum_i P_{ji} \, V_i \, Q_i

where :math:`P_{ji}` is the collision probability and :math:`V_i` is the
volume of region :math:`i`.

.. note::

   **Convention**: in this codebase, :math:`P_{ij}` is indexed as
   :math:`P[\text{birth}_i, \text{collision}_j]`.  The flux update
   uses :math:`P^T`: ``phi = P_inf.T @ source``.
   This is implemented in :meth:`CPSolver.solve_fixed_source`.
   See :ref:`why-p-transpose` for the full derivation.


Definition of Collision Probabilities
======================================

Within-Cell Probabilities and Complementarity
-----------------------------------------------

For a cell with :math:`N` sub-regions, the **within-cell collision
probability** :math:`P_{ij}^{\text{cell}}` is the probability of first
collision in :math:`j` without leaving the cell.  **Complementarity**:

.. math::
   :label: complementarity

   \sum_{j=1}^{N} P_{ij}^{\text{cell}} + P_{i,\text{out}} = 1

where :math:`P_{i,\text{out}}` is the escape probability.  In the code:
``P_out = 1 - P_cell.sum(axis=1)`` (:meth:`CPMesh._apply_white_bc`).
Verified by ``tests/cp/test_properties.py::test_row_sums`` for all three
coordinate systems.


Reciprocity
-----------

From detailed balance :cite:`Hebert2009` section 3.2:

.. math::
   :label: reciprocity

   \Sigt{i} \, V_i \, P_{ij}^{\text{cell}}
   = \Sigt{j} \, V_j \, P_{ji}^{\text{cell}}

**Why reciprocity holds.**  Time-reversal invariance: a neutron born in
:math:`i` colliding in :math:`j` traces a path identical (in reverse) to
one born in :math:`j` colliding in :math:`i`.  The :term:`optical thickness`
along any chord is direction-independent.  The factor :math:`\Sigt{i} V_i`
converts from "per neutron born" (probability) to "per unit source
intensity" (rate), accounting for different source strengths in regions of
different sizes and cross sections.

**Practical consequence.**  The CP matrix need only be computed for
:math:`j \ge i`; the lower triangle follows from:

.. math::
   :label: reciprocity-lower-triangle

   P_{ji}^{\text{cell}} = P_{ij}^{\text{cell}}
   \cdot \frac{\Sigt{i} \, V_i}{\Sigt{j} \, V_j}

This halves the computation cost.  In the code,
:meth:`CPMesh._normalize_rcp` divides the reduced collision probability
by :math:`\Sigt{i} V_i` for each row.  Reciprocity is verified by
``tests/cp/test_properties.py::test_reciprocity`` and extended to multi-group
by ``tests/cp/test_verification.py::TestMultiGroupProperties::test_reciprocity_multigroup``.


Escape and Re-entry (White Boundary Condition)
----------------------------------------------

For an infinite lattice, a neutron escaping one cell immediately enters
an identical neighbour.  The **white boundary condition** assumes the
re-entering angular distribution is **isotropic** --- i.e., the neutron
forgets its direction upon re-entry.

The surface-to-region probability is:

.. math::
   :label: surface-to-region

   P_{\text{in},j} = \frac{\Sigt{j} \, V_j \, P_{j,\text{out}}}{S}

where :math:`S` is the cell surface area, accessed uniformly via
``mesh.areas[-1]`` (:func:`~geometry.coord.compute_areas_1d`):

.. list-table::
   :header-rows: 1
   :widths: 30 30

   * - Coordinate system
     - Surface area :math:`S`
   * - Cartesian (slab)
     - :math:`1` (per unit transverse area)
   * - Cylindrical
     - :math:`2\pi R_{\text{cell}}`
   * - Spherical
     - :math:`4\pi R_{\text{cell}}^2`

**Derivation of** :math:`P_{\text{in},j}`.  By reciprocity between the
surface source and the volume source in region :math:`j`:

.. math::
   :label: pin-from-reciprocity

   \frac{S}{4} \cdot P_{\text{in},j}
   = \Sigt{j} \, V_j \cdot P_{j,\text{out}}

The factor :math:`S/4` is the effective surface source strength for an
isotropic inward flux on a convex surface.  It originates from the
**Cauchy--Dirac mean chord length theorem**: :math:`\bar{\ell} = 4V/S`
for a convex body (:cite:`Hebert2009` section 3.3).  Physically: an isotropic
flux crossing a convex surface sees, on average, a path length of
:math:`4V/S` through the interior.  Solving :eq:`pin-from-reciprocity`:

.. math::

   P_{\text{in},j} = \frac{4 \Sigt{j} V_j P_{j,\text{out}}}{S}

In the standard CP formulation :cite:`Stamm1983` section 3.5, the
surface-to-region probability is defined per unit inward current
:math:`J^-`, and the normalisation convention absorbs the factor of 4.
ORPHEUS uses this convention in the white-BC transform::

    # White-BC closure (geometry-agnostic)
    P_in = sig_t * V * P_out / S_cell

The surface-to-surface probability is:

.. math::
   :label: surface-to-surface

   P_{\text{in,out}} = 1 - \sum_j P_{\text{in},j}

The same formula appears in all three derivation scripts (e.g.,
``orpheus/derivations/continuous/flat_source_cp/slab.py``, line ``P_in = sig_t_g * t_arr * P_out``
with the slab convention :math:`S = 1`, :math:`V = t`; and
``orpheus/derivations/continuous/flat_source_cp/cylinder.py``, line ``S_cell = 2.0 * np.pi * r_cell``
with cylindrical :math:`V = \pi(R_k^2 - R_{k-1}^2)`).


Infinite-Lattice CP Matrix
---------------------------

The infinite-lattice CP accounts for neutrons that escape, re-enter,
possibly escape again (geometric series):

.. math::
   :label: p-inf

   P_{ij}^{\infty}
   = P_{ij}^{\text{cell}}
     + \frac{P_{i,\text{out}} \, P_{\text{in},j}}
            {1 - P_{\text{in,out}}}

This formula is **identical for all three geometries** when expressed
in terms of :math:`V_i` and :math:`S`.  It is implemented in the
white-BC transform selected from :attr:`CPMesh.BC_REGISTRY` and
independently in all three derivation scripts (e.g., ``orpheus/derivations/continuous/flat_source_cp/slab.py``:
``P_inf_g[:,:,g] = P_cell + np.outer(P_out, P_in) / (1.0 - P_inout)``).

.. _cp-bc-registry:

Boundary Condition Infrastructure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The CP solver uses the project-wide ``BC_REGISTRY`` pattern for boundary
condition resolution.  The boundary condition is **declared** on the base
geometry via :class:`~geometry.mesh.BC` on :attr:`Mesh1D.bc_right
<geometry.mesh.Mesh1D.bc_right>` (the outer cell surface), and
**resolved** at :class:`CPMesh` construction time against the registry:

.. code-block:: python

   # Declaration (on geometry)
   mesh = Mesh1D(..., bc_right=BC("vacuum"))

   # Resolution (in CPMesh.__init__)
   factory = CPMesh.BC_REGISTRY[bc.kind]
   self._bc_transform = factory(self, bc)

:attr:`CPMesh.BC_REGISTRY` currently supports two boundary conditions:

.. list-table:: CP Boundary Conditions
   :header-rows: 1
   :widths: 15 30 40

   * - Kind
     - Physics
     - Effect on :math:`P^{\infty}`
   * - ``"white"`` (default)
     - Isotropic re-entry (infinite lattice)
     - Geometric series :eq:`p-inf` applied to :math:`P^{\text{cell}}`
   * - ``"vacuum"``
     - No re-entry (isolated cell)
     - :math:`P^{\infty} = P^{\text{cell}}` (rows sum to < 1)

The default is ``BC("white")``, matching the infinite-lattice assumption
used throughout the CP derivation above.  The vacuum BC is useful for
studying isolated fuel pins where neutrons escaping the cell are lost.

:meth:`CPMesh.compute_pinf_group` calls ``self._bc_transform(P_cell,
sig_t_g)`` to apply whichever BC was resolved at construction time.
Factory docstrings serve as descriptions for programmatic query::

    >>> {k: v.__doc__ for k, v in CPMesh.BC_REGISTRY.items()}

.. plot::
   :caption: White-BC geometric series: a neutron born in region :math:`i` either collides within the cell (:math:`P_{ij}^{\text{cell}}`) or escapes, re-enters isotropically, and the chain repeats.  The infinite sum converges to :math:`P_{ij}^{\infty}`.

   import numpy as np
   import matplotlib.pyplot as plt
   import matplotlib.patches as mpatches

   fig, ax = plt.subplots(figsize=(14, 4))

   cell_w = 2.0
   cell_h = 1.5
   gap = 0.3
   n_cells = 4

   for c in range(n_cells):
       x0 = c * (cell_w + gap)
       rect = mpatches.FancyBboxPatch(
           (x0, 0), cell_w, cell_h, boxstyle="round,pad=0.05",
           facecolor='#ecf0f1', edgecolor='black', linewidth=2)
       ax.add_patch(rect)

       fuel = plt.Circle((x0 + cell_w / 2, cell_h / 2), 0.35,
                         color='#e74c3c', alpha=0.6)
       ax.add_patch(fuel)

       if c == 0:
           ax.text(x0 + cell_w / 2, 0.15, 'cell', fontsize=9,
                   ha='center', color='gray')

   arrow_y = cell_h / 2
   ax.annotate(r'born in $i$', xy=(0.6, arrow_y), fontsize=10,
               ha='center', color='red', fontweight='bold',
               xytext=(0.6, arrow_y + 0.5),
               arrowprops=dict(arrowstyle='->', color='red', lw=1.5))

   ax.annotate('', xy=(1.4, arrow_y + 0.15),
               xytext=(0.8, arrow_y + 0.15),
               arrowprops=dict(arrowstyle='->', color='green', lw=2))
   ax.text(1.1, arrow_y + 0.35, r'$P_{ij}^{\mathrm{cell}}$',
           fontsize=10, ha='center', color='green')

   ax.annotate('', xy=(cell_w + gap / 2, arrow_y - 0.15),
               xytext=(1.5, arrow_y - 0.15),
               arrowprops=dict(arrowstyle='->', color='blue', lw=2))
   ax.text(1.85, arrow_y - 0.4, r'$P_{i,\mathrm{out}}$',
           fontsize=10, ha='center', color='blue')

   x1 = cell_w + gap
   ax.annotate('', xy=(x1 + 0.5, arrow_y),
               xytext=(x1 - gap / 2, arrow_y),
               arrowprops=dict(arrowstyle='->', color='blue', lw=2,
                               linestyle='dashed'))
   ax.text(x1 + 0.15, arrow_y + 0.25, 'white BC\n(isotropic)',
           fontsize=8, ha='center', color='blue', style='italic')

   ax.annotate('', xy=(x1 + 1.4, arrow_y + 0.15),
               xytext=(x1 + 0.6, arrow_y + 0.15),
               arrowprops=dict(arrowstyle='->', color='green', lw=2))
   ax.text(x1 + 1.0, arrow_y + 0.35, r'$P_{\mathrm{in},j}$',
           fontsize=10, ha='center', color='green')

   ax.annotate('', xy=(2 * (cell_w + gap) - gap / 2, arrow_y - 0.15),
               xytext=(x1 + 1.5, arrow_y - 0.15),
               arrowprops=dict(arrowstyle='->', color='orange', lw=2))
   ax.text(x1 + 1.85, arrow_y - 0.4, r'$P_{\mathrm{in,out}}$',
           fontsize=10, ha='center', color='orange')

   x2 = 2 * (cell_w + gap)
   ax.annotate('', xy=(x2 + 0.5, arrow_y),
               xytext=(x2 - gap / 2, arrow_y),
               arrowprops=dict(arrowstyle='->', color='orange', lw=2,
                               linestyle='dashed'))

   x3 = 3 * (cell_w + gap)
   ax.text(x3 + cell_w / 2, cell_h / 2, r'$\cdots$',
           fontsize=24, ha='center', va='center')

   ax.text(
       n_cells * (cell_w + gap) / 2, -0.6,
       r'$P_{ij}^{\infty} = P_{ij}^{\mathrm{cell}}'
       r' + \frac{P_{i,\mathrm{out}} \cdot P_{\mathrm{in},j}}'
       r'{1 - P_{\mathrm{in,out}}}$',
       fontsize=14, ha='center',
       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

   ax.set_xlim(-0.5, n_cells * (cell_w + gap))
   ax.set_ylim(-1.2, cell_h + 0.8)
   ax.set_aspect('equal')
   ax.set_xticks([])
   ax.set_yticks([])
   ax.set_title(
       r'White-BC Geometric Series: $P^{\infty}$ from Escape and Re-entry',
       fontsize=13)
   plt.tight_layout()

**Derivation of the geometric series.**  A neutron born in :math:`i`:

- Collides in :math:`j` on the first pass: probability :math:`P_{ij}^{\text{cell}}`.
- Escapes: probability :math:`P_{i,\text{out}}`.
- After escaping, re-enters isotropically and collides in :math:`j`:
  probability :math:`P_{\text{in},j}`.
- Or traverses without collision: probability :math:`P_{\text{in,out}}`.
- Each subsequent escape-and-re-entry multiplies by :math:`P_{\text{in,out}}`.

Summing:

.. math::
   :label: cp-infinite-lattice-sum

   P_{ij}^{\infty} &= P_{ij}^{\text{cell}}
     + P_{i,\text{out}} P_{\text{in},j}
       (1 + P_{\text{in,out}} + P_{\text{in,out}}^2 + \cdots) \\
   &= P_{ij}^{\text{cell}}
     + \frac{P_{i,\text{out}} P_{\text{in},j}}{1 - P_{\text{in,out}}}

The series converges because :math:`P_{\text{in,out}} < 1` (some
fraction of re-entering neutrons must eventually collide).

**Row sum property:** :math:`\sum_j P_{ij}^{\infty} = 1`.  Proof:
substitute complementarity :eq:`complementarity` for both sums:

.. math::

   \sum_j P_{ij}^{\infty}
   &= \sum_j P_{ij}^{\text{cell}}
     + P_{i,\text{out}} \frac{\sum_j P_{\text{in},j}}{1 - P_{\text{in,out}}} \\
   &= (1 - P_{i,\text{out}})
     + P_{i,\text{out}} \frac{1 - P_{\text{in,out}}}{1 - P_{\text{in,out}}} = 1

Verified numerically by ``tests/cp/test_properties.py::test_row_sums``.

.. warning::

   The white BC is an **approximation**.  The true angular distribution
   at a cell boundary is anisotropic --- neutrons preferentially stream
   in the direction of the flux gradient.  The white BC smears this
   anisotropy into an isotropic re-entry, which overestimates the
   flux in optically thin regions and underestimates it in thick ones.
   For the 1G slab benchmark, the CP method
   (white BC) gives :math:`\keff \approx 1.272` while SN (reflective BC)
   gives :math:`\approx 1.261` --- a ~1% discrepancy entirely due to
   the white-BC approximation.  See :ref:`white-bc-quality`.


The Three CP Kernels
=====================

The CP method requires a geometry-specific **kernel function**
:math:`F(\tau)` that encodes the angular averaging.  The kernel arises
from integrating the point-to-point transmission :math:`e^{-\tau}` over
the angular degrees of freedom eliminated by the geometry's symmetry.
The number of remaining angular integrations determines which special
function appears:

.. list-table::
   :header-rows: 1
   :widths: 15 25 15 25

   * - Geometry
     - Kernel :math:`F(\tau)`
     - :math:`F(0)`
     - Residual quadrature
   * - Slab
     - :math:`E_3(\tau) = \int_0^1 \mu \, e^{-\tau/\mu} \, d\mu`
     - :math:`1/2`
     - None (analytical)
   * - Cylinder
     - :math:`\text{Ki}_4(\tau) = \int_\tau^\infty \text{Ki}_3(t)\,dt`
     - :math:`\approx 0.4244`
     - :math:`y`-quadrature over chord heights
   * - Sphere
     - :math:`e^{-\tau}`
     - :math:`1`
     - :math:`y`-quadrature (:math:`y`-weighted)

In all cases, the **reduced collision probability** is built from a
**second-difference** of the kernel :eq:`second-diff-general`.  This
common structure is why all three geometries are handled by a single
:class:`CPMesh` class.

.. plot::
   :caption: The three CP kernel functions: :math:`E_3` (slab), :math:`\text{Ki}_3` / :math:`\text{Ki}_4` (cylinder), and :math:`e^{-\tau}` (sphere).  All decay exponentially; they differ in :math:`F(0)` and rate of decay.

   import numpy as np
   import matplotlib.pyplot as plt
   from scipy.special import expn
   from scipy.integrate import quad

   x = np.linspace(0.001, 6, 500)

   e3 = expn(3, x)

   ki3 = np.array([
       quad(lambda t, xx=xx: np.exp(-xx / np.sin(t)) * np.sin(t),
            0, np.pi / 2)[0]
       for xx in x
   ])

   dx = x[1] - x[0]
   ki4 = np.cumsum(ki3[::-1])[::-1] * dx

   exp_kernel = np.exp(-x)

   fig, ax = plt.subplots(figsize=(10, 6))
   ax.semilogy(x, e3, '-r', linewidth=2, label=r'$E_3(\tau)$ (slab)')
   ax.semilogy(x, ki3, '-g', linewidth=2,
               label=r'$\mathrm{Ki}_3(\tau)$ (cylinder, integrand)')
   ax.semilogy(x, ki4, '--g', linewidth=2,
               label=r'$\mathrm{Ki}_4(\tau)$ (cylinder, kernel)')
   ax.semilogy(x, exp_kernel, '-b', linewidth=2,
               label=r'$e^{-\tau}$ (sphere)')

   ax.axhline(0.5, color='red', linestyle=':', alpha=0.5)
   ax.text(5.5, 0.55, r'$E_3(0)=\frac{1}{2}$', fontsize=10, color='red')
   ax.axhline(1.0, color='blue', linestyle=':', alpha=0.5)
   ax.text(5.5, 1.1, r'$e^0 = 1$', fontsize=10, color='blue')

   ax.set_xlabel(r'Optical thickness $\tau$', fontsize=12)
   ax.set_ylabel(r'$F(\tau)$', fontsize=12)
   ax.set_title('CP Kernel Functions by Geometry')
   ax.legend(fontsize=11, loc='lower left')
   ax.set_ylim(1e-4, 2)
   ax.grid(True, alpha=0.3)
   plt.tight_layout()

**Derivation source:** The kernels are implemented in
``orpheus/derivations/common/kernels.py`` (slab, arbitrary-precision) and
``orpheus/derivations/continuous/flat_source_cp/geometry.py`` (cylinder double-precision fast
path).  Slab :math:`E_3`: :func:`~orpheus.derivations.common.kernels.e3` and
:func:`~orpheus.derivations.common.kernels.e3_vec` (wrappers over
:func:`scipy.special.expn`). Cylinder :math:`\mathrm{Ki}_3`:
:func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` — a Chebyshev interpolant
of :math:`e^{\tau}\,\mathrm{Ki}_3(\tau)` built from
:func:`~orpheus.derivations.common.kernels.ki_n_mp` at 30 dps. The legacy
``BickleyTables`` tabulation that preceded ``_ki3_mp`` was retired
in Phase B.4 (commit ``6badbe5``, Issue #94).


Why the Kernel Differs by Geometry
-----------------------------------

**Slab (1D):** A neutron in a slab travels at angle :math:`\theta` to
the normal.  The optical path scales as :math:`\tau / \mu` where
:math:`\mu = \cos\theta`.  The angular flux contribution includes a
factor :math:`\mu` (projected area :math:`dA = \mu\,dA_\perp`).
Averaging the transmission :math:`e^{-\tau/\mu}` over the forward
half-space :math:`\mu \in [0, 1]`:

.. math::
   :label: e3-def

   E_3(\tau) = \int_0^1 \mu \, e^{-\tau/\mu} \, d\mu

This is one angular integration, leaving a function of :math:`\tau`
only.  :math:`E_3(0) = 1/2` and :math:`E_3(\tau) \to 0` exponentially.
Computed analytically via :func:`scipy.special.expn` (wrapped as
:func:`_e3` in the solver and :func:`~orpheus.derivations.common.kernels.e3` in the
derivations).

**Cylinder (2D):** The neutron travels at polar angle :math:`\theta` to
the cylinder axis and crosses the annular section at impact parameter
:math:`y`.  The chord length depends on :math:`y` (geometry) and the
optical path depends on :math:`\theta` (as :math:`\tau / \sin\theta`).
Integrating over the polar angle yields the **Bickley--Naylor function**:

.. math::
   :label: ki3-def

   \text{Ki}_3(\tau) = \int_0^{\pi/2}
     e^{-\tau / \sin\theta} \sin\theta \, d\theta

The :math:`\sin\theta` weight arises from the solid-angle measure in
cylindrical coordinates: :math:`d\Omega = \sin\theta\,d\theta\,d\varphi`.
The azimuthal angle :math:`\varphi` integrates trivially by symmetry;
the projected ray length onto the cross-sectional plane is proportional
to :math:`\sin\theta`; and the probability of a ray at angle :math:`\theta`
is proportional to :math:`\sin\theta` from the :math:`d\Omega` measure.

This eliminates one angular dimension, but the :math:`y`-integration over
chord heights **remains as numerical quadrature** --- the second angular
dimension that cylindrical symmetry does NOT eliminate.

The antiderivative :math:`\text{Ki}_4(\tau) = \int_\tau^\infty
\text{Ki}_3(t)\,dt` appears in the second-difference formula because the
derivation involves double antiderivatives of the point-to-point kernel
(see :ref:`second-diff-derivation`).  **Why Ki**\ :sub:`4` **not Ki**\ :sub:`3`:
the slab uses :math:`E_3` (double antiderivative of the point kernel
:math:`E_1`); the cylinder uses :math:`\text{Ki}_4` (single antiderivative
of :math:`\text{Ki}_3`, which is already one integration up from the raw
transmission).  The second-difference formula cancels one antiderivative
level, leaving :math:`\text{Ki}_4`.

The Ki\ :sub:`3` function plays the same role for cylindrical geometry
that :math:`E_3` plays for slab geometry: it represents the probability
of a neutron travelling a certain optical distance in the medium
:cite:`Carlvik1966`.

:math:`\text{Ki}_3(0) = 1` (under the sine-weighted convention used
in :eq:`ki3-def`) and :math:`\text{Ki}_3(x) \to 0` exponentially.

.. note::

   The equation label :eq:`ki3-def` uses the legacy ORPHEUS
   convention :math:`\text{Ki}_3(\tau) = \int_0^{\pi/2} e^{-\tau/\sin\theta}
   \sin\theta\,d\theta`. Under the Abramowitz & Stegun substitution
   :math:`\theta \to \pi/2 - t`, this is actually canonical
   :math:`\mathrm{Ki}_2^{\text{A\&S}}(\tau)`; the kernel consumed by
   the second-difference CP formula is the next function up in the
   hierarchy (called :math:`\text{Ki}_4` below) which equals
   canonical :math:`\mathrm{Ki}_3^{\text{A\&S}}`. The legacy naming
   is preserved in this theory page for cross-consistency with
   existing ``verifies("ki3-def")`` decorators on
   ``tests/cp/test_cylinder.py``,
   ``tests/derivations/test_cp_geometry.py``,
   ``tests/cp/test_verification.py``, and
   ``tests/derivations/test_peierls_cylinder_prefactor.py``. See
   :doc:`/verification/reference_solutions` §"Legacy naming
   discrepancy in ``BickleyTables``" for the full postmortem.

Neither :math:`\text{Ki}_3` nor :math:`\text{Ki}_4` has a
closed-form expression.  The derivation path is now
:func:`~orpheus.derivations.common.kernels.ki_n_mp` (arbitrary precision via
:func:`mpmath.quad`) with a double-precision fast path through
:func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` (a Chebyshev interpolant of
the canonical-:math:`\mathrm{Ki}_3` scaled kernel
:math:`e^{\tau}\,\mathrm{Ki}_3(\tau)`); the solver and the
derivation share the *same* code path so there is no kernel-split
bias across the CP V&V stack. See :ref:`ki-table-construction` for
the full account of what this replaced.

**Sphere (3D):** Full spherical symmetry absorbs all angular variables
into the :math:`y`-quadrature weight.  The kernel is simply
:math:`F(\tau) = e^{-\tau}` --- no special functions needed.  The extra
factor of :math:`y` in the quadrature weight comes from the spherical area
element :math:`2\pi y\,dy` (ring of sources) versus the cylindrical
:math:`2\,dy` (line of sources).  In
:meth:`CPMesh._setup_spherical`::

    # Spherical weight: extra factor of y in the quadrature
    self._y_wts = self._y_wts * self._y_pts


.. _second-diff-derivation:

The Second-Difference Formula: Full Derivation
-------------------------------------------------

The reduced collision probability uses the four-term second-difference:

.. math::
   :label: second-diff-general

   \Delta_2[F](\tau_i, \tau_j, g)
   = F(g) - F(g + \tau_i) - F(g + \tau_j) + F(g + \tau_i + \tau_j)

where :math:`g` is the optical gap, :math:`\tau_i` and :math:`\tau_j` are
the optical thicknesses of the source and target regions.

**Setup.** Regions :math:`i` and :math:`j` separated by optical gap
:math:`g`.  A neutron born at optical depth :math:`s` within region
:math:`i` (:math:`0 \le s \le \tau_i`) must traverse
:math:`d = (\tau_i - s) + g + t` to reach depth :math:`t` in region
:math:`j` (:math:`0 \le t \le \tau_j`).

**Step 1: Source-position average.** Integrate the kernel over the birth
position :math:`s` within region :math:`i`.  Substituting
:math:`u = (\tau_i - s) + g + t` (so :math:`du = -ds`):

.. math::

   I(t) &= \int_0^{\tau_i} F\bigl((\tau_i - s) + g + t\bigr) \, ds \\
        &= \int_{g+t}^{\tau_i + g + t} F(u) \, du \\
        &= \hat{F}(\tau_i + g + t) - \hat{F}(g + t)

where :math:`\hat{F}(x) = \int_0^x F(u)\,du` is the antiderivative.

**Step 2: Collision-position average.** Integrate over the collision
position :math:`t` within region :math:`j`:

.. math::

   \text{rcp}_{ij} &= \int_0^{\tau_j} I(t) \, dt \\
   &= \int_0^{\tau_j} \bigl[\hat{F}(\tau_i + g + t) - \hat{F}(g + t)\bigr] dt

Evaluating each integral using the double antiderivative
:math:`\hat{\hat{F}}(x) = \int_0^x \hat{F}(u)\,du`:

.. math::

   \int_0^{\tau_j} \hat{F}(g + t) \, dt
   &= \hat{\hat{F}}(g + \tau_j) - \hat{\hat{F}}(g) \\
   \int_0^{\tau_j} \hat{F}(\tau_i + g + t) \, dt
   &= \hat{\hat{F}}(\tau_i + g + \tau_j) - \hat{\hat{F}}(\tau_i + g)

Combining:

.. math::
   :label: rcp-from-double-antideriv

   \text{rcp}_{ij} =
     \hat{\hat{F}}(g)
   - \hat{\hat{F}}(g + \tau_i)
   - \hat{\hat{F}}(g + \tau_j)
   + \hat{\hat{F}}(g + \tau_i + \tau_j)

**Step 3: Identify** :math:`\hat{\hat{F}}` **per geometry.**

.. list-table::
   :header-rows: 1
   :widths: 15 25 25 25

   * - Geometry
     - Point-to-point kernel
     - 1st antiderivative
     - :math:`\hat{\hat{F}}`
   * - Slab
     - :math:`E_1(\tau)`
     - :math:`-E_2(\tau)`
     - :math:`E_3(\tau)`
   * - Cylinder
     - :math:`\text{Ki}_3(\tau)`
     - :math:`-\text{Ki}_4(\tau)`
     - :math:`\text{Ki}_5(\tau)`
   * - Sphere
     - :math:`e^{-\tau}`
     - :math:`-e^{-\tau}`
     - :math:`e^{-\tau}`

For the slab, the :math:`E_n` functions satisfy
:math:`E_n'(\tau) = -E_{n-1}(\tau)`.  The code evaluates
:eq:`second-diff-general` directly with :math:`F = E_3` (slab),
:math:`F = \mathrm{Ki}_3^{\text{A\&S}}` (cylinder; the
ORPHEUS-local "Ki\ :sub:`4`\ " of the legacy narrative is canonical
A&S :math:`\mathrm{Ki}_3`), or :math:`F = e^{-\tau}` (sphere).

**Derivation source.** After the Phase B.2 unification, all three
flat-source CP modules dispatch through a single
geometry-invariant operator implemented as the module-level free
function :func:`~orpheus.derivations.continuous.flat_source_cp.geometry._second_difference`::

    def _second_difference(kernel, gap, tau_i, tau_j):
        return (kernel(gap)
                - kernel(gap + tau_i)
                - kernel(gap + tau_j)
                + kernel(gap + tau_i + tau_j))

The per-geometry kernel is supplied by the
:class:`~orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry` singleton:

- :data:`~orpheus.derivations.continuous.flat_source_cp.geometry.SLAB` — ``kernel_F3 = e3_vec``
  (:func:`scipy.special.expn`)
- :data:`~orpheus.derivations.continuous.flat_source_cp.geometry.CYLINDER_1D` —
  ``kernel_F3 = _ki3_mp`` (Chebyshev interpolant of
  canonical :math:`\mathrm{Ki}_3^{\text{A\&S}}`)
- :data:`~orpheus.derivations.continuous.flat_source_cp.geometry.SPHERE_1D` —
  ``kernel_F3 = _exp_kernel`` (``np.exp(-tau)``)

so the four-term structure is written exactly once in the whole
derivations package.

**SymPy verification of the four-term structure.** The following
script verifies the derivation in :eq:`rcp-from-double-antideriv`
symbolically. It builds the exact double integral of the transmission
kernel over source position :math:`s` and collision position :math:`t`,
substitutes the antiderivative identity, and checks that the result
matches the four-term second-difference pattern for a generic
:math:`F`. The same :class:`Delta2` structure emerges independent of
whether :math:`F = E_3`, :math:`\text{Ki}_4`, or :math:`e^{-\tau}` —
proving that the three geometries share one algebraic form::

    import sympy as sp

    # Symbolic setup: optical variables, arbitrary antiderivative F
    s, t, tau_i, tau_j, g = sp.symbols('s t tau_i tau_j g', positive=True)
    F = sp.Function('F')                       # generic kernel
    Fhat  = lambda x: sp.integrate(F(s), (s, 0, x))            # F̂
    Fhh   = lambda x: sp.integrate(Fhat(s),  (s, 0, x))        # F̂̂

    # Step 1: average over birth position s in region i, fixed t
    I_t = sp.integrate(F((tau_i - s) + g + t), (s, 0, tau_i))
    I_t = sp.simplify(I_t.rewrite(sp.Integral))

    # Step 2: average over collision position t in region j
    rcp = sp.integrate(I_t, (t, 0, tau_j))

    # Expected four-term pattern from eq:`rcp-from-double-antideriv`
    rcp_expected = ( Fhh(g)
                   - Fhh(g + tau_i)
                   - Fhh(g + tau_j)
                   + Fhh(g + tau_i + tau_j) )

    # The two must be equal for ANY F (not just E3/Ki4/exp)
    assert sp.simplify(rcp - rcp_expected) == 0

    # Instantiate per-geometry kernels and confirm
    for name, kernel in [
        ('slab',     lambda x: sp.expint(3, x)),         # E_3
        ('sphere',   lambda x: sp.exp(-x)),              # e^(-tau)
        # ('cylinder', lambda x: Ki4(x))                 # Bickley; numerical
    ]:
        Fsym = sp.Function(f'F_{name}')
        check = sp.simplify(
            (rcp_expected.rewrite(sp.Integral)
             .subs(F, kernel))
            .doit()
        )
        print(f'{name}: {check}')

This script is **not** run as part of the test suite (SymPy's
treatment of ``expint`` does not always simplify to a closed form
suitable for boolean equality), but it can be pasted into any Python
session as an independent verification of the derivation. The
generic-kernel assertion (``assert sp.simplify(rcp - rcp_expected)
== 0``) proves the four-term structure without reference to any
specific kernel — which is exactly the claim that slab, cylinder,
and sphere share one algebraic form.

.. plot::
   :caption: Geometric meaning of the four second-difference terms :math:`\Delta_2[F](\tau_i, \tau_j, g)`.  Green arrows mark the two positive terms :math:`+F(g)` and :math:`+F(g+\tau_i+\tau_j)`; orange arrows mark the two negative terms.

   import numpy as np
   import matplotlib.pyplot as plt
   import matplotlib.patches as mpatches

   fig, ax = plt.subplots(figsize=(12, 4))

   colors = ['#e74c3c', '#e67e22', '#bdc3c7', '#bdc3c7', '#3498db']
   labels = ['Region $i$', '', 'Gap', '', 'Region $j$']
   widths = [1.0, 0.3, 0.8, 0.3, 1.0]
   positions = np.cumsum([0] + widths)

   for k in range(5):
       rect = mpatches.FancyBboxPatch(
           (positions[k], 0), widths[k], 1.0,
           boxstyle="square,pad=0",
           facecolor=colors[k], edgecolor='black',
           alpha=0.5 if k in [1, 2, 3] else 0.7)
       ax.add_patch(rect)
       if labels[k]:
           ax.text(positions[k] + widths[k] / 2, 0.5, labels[k],
                   ha='center', va='center', fontsize=11,
                   fontweight='bold')

   bnd_labels = [
       (positions[0], r'$x_{i-1}$'),
       (positions[1], r'$x_i$'),
       (positions[3], r'$x_{j-1}$'),
       (positions[4], r'$x_j$'),
   ]
   for xp, label in bnd_labels:
       ax.plot([xp, xp], [-0.1, 1.1], 'k--', linewidth=1, alpha=0.5)
       ax.text(xp, -0.2, label, ha='center', fontsize=11)

   ax.annotate('', xy=(positions[1], 1.2), xytext=(positions[0], 1.2),
               arrowprops=dict(arrowstyle='<->', color='red', lw=2))
   ax.text((positions[0] + positions[1]) / 2, 1.3, r'$\tau_i$',
           ha='center', fontsize=12, color='red')

   ax.annotate('', xy=(positions[4] + widths[4], 1.2),
               xytext=(positions[3] + widths[3], 1.2),
               arrowprops=dict(arrowstyle='<->', color='blue', lw=2))
   ax.text((positions[3] + widths[3] + positions[4] + widths[4]) / 2, 1.3,
           r'$\tau_j$', ha='center', fontsize=12, color='blue')

   ax.annotate('', xy=(positions[3], 1.55), xytext=(positions[1], 1.55),
               arrowprops=dict(arrowstyle='<->', color='purple', lw=2))
   ax.text((positions[1] + positions[3]) / 2, 1.65,
           r'gap $g$', ha='center', fontsize=12, color='purple')

   y_bot = -0.5
   terms = [
       (positions[1], positions[3], r'$+F(g)$', 'green'),
       (positions[0], positions[3], r'$-F(g+\tau_i)$', 'orange'),
       (positions[1], positions[4] + widths[4], r'$-F(g+\tau_j)$', 'orange'),
       (positions[0], positions[4] + widths[4],
        r'$+F(g+\tau_i+\tau_j)$', 'green'),
   ]
   for k, (x1, x2, label, color) in enumerate(terms):
       y = y_bot - k * 0.35
       ax.annotate('', xy=(x2, y), xytext=(x1, y),
                   arrowprops=dict(arrowstyle='<->', color=color,
                                   lw=1.5))
       ax.text(x2 + 0.1, y, label, fontsize=10, va='center',
               color=color)

   ax.set_xlim(-0.3, 5.0)
   ax.set_ylim(-2.0, 2.0)
   ax.set_aspect('equal')
   ax.set_yticks([])
   ax.set_xticks([])
   ax.set_title('Second-Difference: Four Boundary Evaluations',
                fontsize=13)
   plt.tight_layout()


The Discretised Second-Difference
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The continuous formula :eq:`rcp-from-double-antideriv` is discretised by
evaluating the double antiderivative at region boundary positions.  Define:

.. math::
   :label: s-integral

   S(i, j, g) = \int_0^{R_i}
     \bigl[\text{Ki}_3(\tau_p(y)) - \text{Ki}_3(\tau_m(y))\bigr] \, dy

The within-cell collision probability for :math:`j \ge i` is:

.. math::
   :label: pcell-from-smat

   P_{ij}^{\text{cell}} = \frac{2}{\Sigt{i} V_i}
   \Bigl[ S(i,j) - S(i{-}1,j) - S(i,j{-}1) + S(i{-}1,j{-}1) \Bigr]
   + \delta_{ij}

where :math:`\delta_{ij}` accounts for self-collision and terms with
index 0 vanish (no region interior to region 1).

.. list-table:: Correspondence: continuous vs discrete
   :header-rows: 1
   :widths: 40 40

   * - Continuous term
     - Discrete boundary evaluation
   * - :math:`\hat{\hat{F}}(g)`
     - :math:`S(i{-}1, j{-}1)` --- both inner boundaries
   * - :math:`\hat{\hat{F}}(g + \tau_i)`
     - :math:`S(i, j{-}1)` --- outer of source, inner of target
   * - :math:`\hat{\hat{F}}(g + \tau_j)`
     - :math:`S(i{-}1, j)` --- inner of source, outer of target
   * - :math:`\hat{\hat{F}}(g + \tau_i + \tau_j)`
     - :math:`S(i, j)` --- both outer boundaries

**Implementation.** ORPHEUS does **not** pre-compute the :math:`S` array.
Instead, :meth:`CPMesh._compute_radial_rcp` evaluates the
second-difference at each :math:`y`-quadrature point and integrates
numerically.  This avoids storing :math:`S` and allows the same code path
for all three geometries by parameterising the kernel function.  The
derivation scripts use the same approach.


Self-Collision: Full Derivation
---------------------------------

The self-collision term (:math:`i = j`) cannot use the gap formula because
source and collision positions overlap.  The double integral is:

.. math::
   :label: self-double-integral

   r_{ii} = \int_0^{\tau_i} \int_0^{\tau_i} F_1(|s - t|) \, dt \, ds

where :math:`F_1` is the point-to-point kernel.

**Evaluation (slab case, :math:`F_1 = E_1`).**  Split at :math:`s = t`
and use symmetry (:math:`s > t` half = :math:`s < t` half):

.. math::

   r_{ii} &= 2 \int_0^{\tau_i} \int_0^s E_1(s - t) \, dt \, ds

Evaluate the inner integral using :math:`\int E_1(x) dx = -E_2(x)`:

.. math::

   \int_0^s E_1(s - t) \, dt
   &= \bigl[-E_2(s - t)\bigr]_{t=0}^{t=s}
   = -E_2(0) + E_2(s)
   = E_2(s) - E_2(0)

Since :math:`E_2(0) = 1`, the inner integral is :math:`E_2(s) - 1`.
The outer integral:

.. math::

   r_{ii} &= 2 \int_0^{\tau_i} \bigl[E_2(s) - 1\bigr] ds \\
          &= 2 \bigl[-E_3(s)\bigr]_0^{\tau_i} - 2\tau_i \\
          &= 2 \bigl[-E_3(\tau_i) + E_3(0)\bigr] - 2\tau_i \\
          &= 2 E_3(0) - 2 E_3(\tau_i) - 2\tau_i

Wait --- let me redo this carefully.  :math:`\int_0^s E_1(s-t) dt`:
substituting :math:`u = s - t`, :math:`du = -dt`:

.. math::

   \int_0^s E_1(u) du = \bigl[-E_2(u)\bigr]_0^s = -E_2(s) + E_2(0)

So the inner integral is :math:`E_2(0) - E_2(s)`.  Continuing:

.. math::

   r_{ii} &= 2 \int_0^{\tau_i} \bigl[E_2(0) - E_2(s)\bigr] ds \\
          &= 2 \bigl[\tau_i E_2(0) + E_3(\tau_i) - E_3(0)\bigr]

Using :math:`E_2(0) = 1` and :math:`E_3(0) = 1/2`:

.. math::

   r_{ii} = 2\tau_i + 2 E_3(\tau_i) - 1

The normalised self-collision probability is
:math:`P_{ii} = r_{ii} / (\Sigt{i} t_i)`.  Since
:math:`\tau_i = \Sigt{i} t_i`:

.. math::

   P_{ii} = \frac{2\tau_i + 2 E_3(\tau_i) - 1}{\tau_i}
          = 2 + \frac{2 E_3(\tau_i) - 1}{\tau_i}
          = 1 + \frac{2 E_3(\tau_i) - 2 E_3(0)}{\tau_i}

which is equivalently written as:

.. math::
   :label: self-collision-probability-slab

   P_{ii} = 1 - \frac{2(E_3(0) - E_3(\tau_i))}{\tau_i}

For **thick regions** (:math:`\tau_i \to \infty`): :math:`E_3(\tau_i) \to 0`,
so :math:`P_{ii} \to 1 - 1/\tau_i \to 1`.  For **thin regions**
(:math:`\tau_i \to 0`): :math:`E_3(\tau_i) \to E_3(0) = 1/2`, so
:math:`P_{ii} \to 0`.  Tested by
``tests/cp/test_verification.py::TestOpticalLimits``.

In the solver code (:meth:`CPMesh._compute_slab_rcp`)::

    # Self-collision (slab)
    rcp[i, i] += sti * t[i] - (0.5 - _e3(tau_i))

This is :math:`\Sigt{i} t_i - (E_3(0) - E_3(\tau_i))`, which equals
:math:`r_{ii}/2` (the factor of 2 from two half-spaces is applied later).

For **cylindrical and spherical** geometries, the same structure holds
but with canonical :math:`\mathrm{Ki}_3^{\text{A\&S}}` or
:math:`e^{-\tau}` replacing :math:`E_3`, and with
:math:`y`-quadrature.  In :meth:`CPMesh._compute_radial_rcp`::

    self_same = 2.0 * chords[i, :] - (2.0 / sti) * (
        kernel_zero - kernel(tau_i)
    )
    rcp[i, i] += 2.0 * sti * np.dot(y_wts, self_same)

The term ``kernel_zero - kernel(tau_i)`` is :math:`F(0) - F(\tau_i)`,
the escape fraction.  The same pattern is expressed once in the
unified Phase B.2 derivation core
:mod:`orpheus.derivations.continuous.flat_source_cp.geometry` — the self-term pairing consumes
:meth:`~orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.kernel_F3_at_zero`
(returns :math:`1/2`, :math:`\pi/4`, or :math:`1` for slab / cyl
/ sph) alongside :meth:`~orpheus.derivations.continuous.flat_source_cp.geometry.FlatSourceCPGeometry.kernel_F3`
for the evaluated-at-:math:`\tau_i` term.


.. _optical-path-construction:

Optical Path Construction Along a Chord
=========================================

This section describes how :math:`\tau_m` and :math:`\tau_p` are
constructed at chord height :math:`y`.  This is the most geometrically
intricate part of the cylindrical and spherical CP implementations.


Half-Chord Lengths
-------------------

Consider :math:`N` concentric regions with outer radii
:math:`R_1 < \cdots < R_N`.  A chord at height :math:`y` intersects
a subset of these.  The **half-chord length** through region :math:`k` is:

.. math::
   :label: chord-length

   \ell_k(y) = \begin{cases}
     \sqrt{R_k^2 - y^2} - \sqrt{R_{k-1}^2 - y^2}
       & y < R_{k-1} \\
     \sqrt{R_k^2 - y^2}
       & R_{k-1} \le y < R_k \\
     0 & y \ge R_k
   \end{cases}

with :math:`R_0 = 0`.

- **Case 1** (:math:`y < R_{k-1}`): Chord passes entirely through region
  :math:`k`, entering at the inner boundary and exiting at the outer.
- **Case 2** (:math:`R_{k-1} \le y < R_k`): Chord originates inside
  region :math:`k` --- no inner intersection.
- **Case 3** (:math:`y \ge R_k`): Chord misses region :math:`k` entirely.

Computed by :func:`_chord_half_lengths` (solver) and independently by
``orpheus/derivations/continuous/flat_source_cp/cylinder.py::_chord_half_lengths`` (derivation).  Both
return shape ``(N, n_y)``.

**Optical half-thickness:** ``tau = sig_t_g[:, None] * chords``.

**Gotcha: the innermost region.**  For :math:`k = 1`, :math:`R_0 = 0`
gives :math:`\sqrt{-y^2}` (imaginary).  Both implementations handle this
with a ``r_in = 0`` / ``r_in > 0`` branch that skips the subtraction.

**Implicit zero-padding.**  When :math:`y \ge R_k`,
:func:`_chord_half_lengths` returns :math:`\ell_k = 0`, so
:math:`\tau_k = 0`.  Boundary positions collapse and the
second-difference evaluates to zero --- no explicit conditional logic is
needed to skip non-intersecting regions.

.. plot::
   :caption: Half-chord length :math:`\ell_k(y)` through region :math:`k = 2` (red outline) for each of the three cases.

   import numpy as np
   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(1, 3, figsize=(15, 5))

   radii = [0.4, 0.7, 1.0]
   colors = ['#e74c3c', '#f1c40f', '#3498db']

   for ax_idx, (y_val, case_title) in enumerate([
       (0.2, r'Case 1: $y < R_{k-1}$'),
       (0.5, r'Case 2: $R_{k-1} \leq y < R_k$'),
       (0.8, r'Case 3: $y \geq R_k$'),
   ]):
       ax = axes[ax_idx]

       for r, color in reversed(list(zip(radii, colors))):
           ax.add_patch(plt.Circle((0, 0), r, color=color, alpha=0.3))
           ax.add_patch(plt.Circle((0, 0), r, fill=False,
                                   edgecolor='black', linewidth=1))

       theta = np.linspace(0, 2 * np.pi, 100)
       ax.fill_between(
           radii[1] * np.cos(theta), radii[1] * np.sin(theta),
           alpha=0, edgecolor='red', linewidth=3, linestyle='-')

       if y_val < radii[-1]:
           x_max = np.sqrt(radii[-1] ** 2 - y_val ** 2)
           ax.plot([-x_max, x_max], [y_val, y_val], 'k-', linewidth=2)

           k = 1
           r_out, r_in = radii[k], radii[k - 1]
           if y_val < r_in:
               x_out = np.sqrt(r_out ** 2 - y_val ** 2)
               x_in = np.sqrt(r_in ** 2 - y_val ** 2)
               ax.plot([x_in, x_out], [y_val, y_val], 'r-', linewidth=5,
                       alpha=0.7, label=r'$\ell_k(y)$')
           elif y_val < r_out:
               x_out = np.sqrt(r_out ** 2 - y_val ** 2)
               ax.plot([0, x_out], [y_val, y_val], 'r-', linewidth=5,
                       alpha=0.7, label=r'$\ell_k(y)$')

       if y_val < radii[-1]:
           ax.annotate('', xy=(0, y_val), xytext=(0, 0),
                       arrowprops=dict(arrowstyle='<->', color='black',
                                       lw=1.5))
           ax.text(0.04, y_val / 2, '$y$', fontsize=11)

       ax.text(0.41, -0.15, r'$R_1$', fontsize=10, color='#e74c3c')
       ax.text(0.71, -0.15, r'$R_2$', fontsize=10, color='#f1c40f')
       ax.text(1.01, -0.15, r'$R_3$', fontsize=10, color='#3498db')

       ax.set_xlim(-1.2, 1.2)
       ax.set_ylim(-1.2, 1.2)
       ax.set_aspect('equal')
       ax.set_title(case_title, fontsize=11)
       ax.grid(True, alpha=0.2)
       if ax_idx == 0:
           ax.legend(loc='upper left', fontsize=10)

   fig.suptitle(
       r'Half-chord length $\ell_k(y)$ through region $k=2$ '
       r'(red outline)', fontsize=13, y=1.02)
   plt.tight_layout()


The Two Optical Paths
----------------------

For a neutron born in region :math:`i` targeting region :math:`j \ge i`,
there are two distinct paths along the chord:

1. **Same-side path** (:math:`\tau_m`): gap between the outer boundary
   of :math:`i` and the inner boundary of :math:`j`:

   .. math::
      :label: tau-m

      \tau_m(y) = \sum_{k=i+1}^{j} \tau_k(y)

   For adjacent regions :math:`j = i+1`, this is :math:`\tau_{i+1}`.
   For self-collision :math:`j = i`, :math:`\tau_m = 0`.

2. **Through-centre path** (:math:`\tau_p`): crosses all inner regions
   twice (inward to the centre, back out):

   .. math::
      :label: tau-p

      \tau_p(y) = \tau_m(y) + 2 \sum_{k=1}^{i} \tau_k(y)


Boundary Position Arrays
--------------------------

Cumulative optical distance from the chord midpoint to the outer boundary
of region :math:`k`:

.. math::
   :label: cumulative-optical-path-radial

   x_k(y) = \sum_{m=1}^{k} \tau_m(y), \quad x_0 = 0

In code (:meth:`CPMesh._compute_radial_rcp`, and identically in all three
derivation scripts)::

    bnd_pos = np.zeros((N + 1, n_y))
    for k in range(N):
        bnd_pos[k + 1, :] = bnd_pos[k, :] + tau[k, :]

The gaps: ``gap_d = bnd_pos[j] - bnd_pos[i+1]`` (same-side) and
``gap_c = bnd_pos[i] + bnd_pos[j]`` (through-centre).


Slab Geometry: The :math:`E_3` Kernel
=======================================

The 1D slab half-cell extends from the reflective centre (:math:`x = 0`)
to the cell edge (:math:`x = L`).  Geometry built via
:func:`~geometry.factories.pwr_slab_half_cell` with
``coord = CoordSystem.CARTESIAN``.

.. plot::
   :caption: Slab half-cell: fuel, cladding, and coolant with reflective (left) and white (right) boundary conditions.

   import numpy as np
   import matplotlib.pyplot as plt
   import matplotlib.patches as mpatches

   regions = [
       ('Fuel', 0.9, '#e74c3c'),
       ('Clad', 0.2, '#95a5a6'),
       ('Cool', 0.7, '#3498db'),
   ]

   fig, ax = plt.subplots(figsize=(10, 2.5))
   x = 0
   for name, width, color in regions:
       rect = mpatches.FancyBboxPatch(
           (x, 0), width, 1, boxstyle="square,pad=0",
           facecolor=color, edgecolor='black', alpha=0.7)
       ax.add_patch(rect)
       ax.text(x + width / 2, 0.5, name, ha='center', va='center',
               fontsize=12, fontweight='bold')
       x += width

   ax.set_xlim(-0.1, 2.0)
   ax.set_ylim(-0.2, 1.3)
   ax.set_xlabel('x (cm)')
   ax.set_aspect('equal')
   ax.axvline(0, color='red', linestyle='--', linewidth=2,
              label='Reflective BC')
   ax.axvline(1.8, color='blue', linestyle='--', linewidth=2,
              label='White BC')
   ax.legend(loc='upper right')
   ax.set_yticks([])
   ax.set_title('Slab Half-Cell Geometry')
   plt.tight_layout()

Let :math:`\tau_k = \Sigt{k} \, t_k` be the optical thickness of region
:math:`k`, and define the cumulative optical path from the cell centre:

.. math::
   :label: cumulative-optical-path-slab

   x_0 = 0, \quad x_{k} = \sum_{m=0}^{k-1} \tau_m

The **reduced collision probability** (unnormalised) is built from two
path types:

1. **Direct path** (same direction along the slab): for :math:`j > i`,

   .. math::
      :label: dd-slab

      \delta_d = E_3(g) - E_3(g + \tau_i) - E_3(g + \tau_j) + E_3(g + \tau_i + \tau_j)

   where :math:`g = x_j - x_{i+1}` is the optical gap between regions.

2. **Reflected path** (through the reflective centre at :math:`x = 0`):

   .. math::
      :label: dc-slab

      \delta_c = E_3(g_c) - E_3(g_c + \tau_i) - E_3(g_c + \tau_j) + E_3(g_c + \tau_i + \tau_j)

   where :math:`g_c = x_i + x_j` is the optical path via the centre.

3. **Self-collision** (:math:`i = j`):

   .. math::
      :label: self-slab

      r_{ii} = \Sigt{i} t_i - \bigl(E_3(0) - E_3(\tau_i)\bigr)

The total reduced CP for :math:`i \ne j` is:

.. math::
   :label: rcp-slab-total

   r_{ij} = \frac{1}{2}(\delta_d + \delta_c)

and the within-cell CP is :math:`P_{ij}^{\text{cell}} = r_{ij} /
(\Sigt{i} \, V_i)`, where :math:`V_i = t_i` for slab geometry.

Implemented in :meth:`CPMesh._compute_slab_rcp`.  Verified element-by-element
against ``orpheus/derivations/continuous/flat_source_cp/slab.py::_slab_cp_matrix`` by
``tests/cp/test_verification.py::TestDirectPinfComparison::test_slab_pinf_matches_derivation``
(tolerance :math:`< 10^{-10}`).


Concentric Cylindrical Geometry: The Ki\ :sub:`3`/Ki\ :sub:`4` Kernel
=======================================================================

The Wigner--Seitz cell replaces the square unit cell by a circle of equal
area.  From area equivalence :math:`p^2 = \pi R_{\text{cell}}^2`:

.. math::
   :label: wigner-seitz

   R_{\text{cell}} = \frac{p}{\sqrt{\pi}}

Geometry built via :meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.wigner_seitz_pin_cell`
(``coord = CoordSystem.CYLINDRICAL``).

.. plot::
   :caption: Wigner--Seitz cylindrical cell with concentric annular regions and a chord at height :math:`y`.

   import numpy as np
   import matplotlib.pyplot as plt

   fig, ax = plt.subplots(figsize=(6, 6))

   regions = [
       (0.9, '#e74c3c', 'Fuel'),
       (1.1, '#95a5a6', 'Clad'),
       (2.03, '#3498db', 'Cool'),
   ]

   for r, color, label in reversed(regions):
       ax.add_patch(plt.Circle((0, 0), r, color=color, alpha=0.5,
                               label=label))
   for r, _, _ in regions:
       ax.add_patch(plt.Circle((0, 0), r, fill=False, edgecolor='black',
                               linewidth=1))

   y_chord = 0.6
   x_max = np.sqrt(regions[-1][0] ** 2 - y_chord ** 2)
   ax.plot([-x_max, x_max], [y_chord, y_chord], 'k-', linewidth=2)
   ax.annotate('chord at height $y$', xy=(0, y_chord),
               xytext=(0.5, 1.5), fontsize=11,
               arrowprops=dict(arrowstyle='->', color='black'))
   ax.plot([0, 0], [0, y_chord], 'k--', alpha=0.5)
   ax.text(0.08, y_chord / 2, '$y$', fontsize=12)

   ax.set_xlim(-2.5, 2.5)
   ax.set_ylim(-2.5, 2.5)
   ax.set_aspect('equal')
   ax.legend(loc='upper right', fontsize=11)
   ax.set_xlabel('x (cm)')
   ax.set_ylabel('y (cm)')
   ax.set_title('Wigner\u2013Seitz Cell with Chord Integration')
   ax.grid(True, alpha=0.3)
   plt.tight_layout()

The CP matrix integrates :math:`\text{Ki}_4` second-differences over
chord heights using composite Gauss--Legendre quadrature
(:func:`_composite_gauss_legendre`) with breakpoints at each annular
boundary to capture chord-length discontinuities.

Second-Difference Formula (Cylindrical)
---------------------------------------

The CP matrix is computed by integrating the Ki\ :sub:`4` second
differences over all chord heights:

.. math::
   :label: second-diff-cyl

   r_{ij} = 2 \int_0^{R_{\text{cell}}}
     \bigl[\text{Ki}_4(g) - \text{Ki}_4(g + \tau_i)
           - \text{Ki}_4(g + \tau_j) + \text{Ki}_4(g + \tau_i + \tau_j)\bigr] \, dy

where :math:`g` is the optical gap between regions :math:`i` and
:math:`j` along the chord.  The factor 2 accounts for the two halves
of the chord (by symmetry of the annular geometry, contributions from
the left and right halves are equal).

Two path types contribute (same as the slab):

1. **Same-side path**: source and target on the same side of the chord
   midpoint.
2. **Through-centre path**: the neutron crosses the chord midpoint
   (optical gap :math:`g_c = x_i + x_j` from the cumulative boundary
   positions).

The self-collision term is:

.. math::
   :label: self-cyl

   r_{ii} = 2\Sigt{i} \int_0^{R_{\text{cell}}}
     \left[2\ell_i(y) - \frac{2}{\Sigt{i}}
       \left(\text{Ki}_4(0) - \text{Ki}_4(\tau_i(y))\right)\right] dy

The :math:`y`-integration is performed with composite Gauss--Legendre
quadrature, with breakpoints at each annular boundary to capture the
chord-length discontinuities.

Implemented in :meth:`CPMesh._compute_radial_rcp` with
``self._kernel = Ki_4``.  Verified against
``orpheus/derivations/continuous/flat_source_cp/cylinder.py::_cylinder_cp_matrix``.


Concentric Spherical Geometry: The Exponential Kernel
======================================================

The spherical cell consists of :math:`N` concentric spherical shells.
Volumes: :math:`V_i = \frac{4}{3}\pi(R_i^3 - R_{i-1}^3)`.  Surface:
:math:`S = 4\pi R_{\text{cell}}^2`.  The kernel is simply
:math:`F(\tau) = e^{-\tau}` --- full 3-D symmetry absorbs all angular
variables.

Chord geometry is identical to cylindrical (:eq:`chord-length`), but
the quadrature weight includes :math:`y`.  For a cylinder (2-D), each
chord at height :math:`y` represents a **line of sources**, giving a
weight proportional to :math:`dy`.  For a sphere (3-D), each chord
represents a **ring of sources** with circumference :math:`2\pi y`,
giving a weight proportional to :math:`y\,dy`.
In :meth:`CPMesh._setup_spherical`::

    # Spherical weight: extra factor of y in the quadrature
    self._y_wts = self._y_wts * self._y_pts

Geometry built via :meth:`~orpheus.geometry.mesh.Mesh1D.from_geometry`
(``coord = CoordSystem.SPHERICAL``).

Second-Difference Formula (Spherical)
-------------------------------------

The reduced collision probability has the same structure as cylindrical,
but with :math:`F = \exp(-\cdot)` and :math:`y`-weighted quadrature:

.. math::
   :label: second-diff-sph

   r_{ij} = 2 \int_0^{R_{\text{cell}}}
     \left[e^{-g} - e^{-(g + \tau_i)}
           - e^{-(g + \tau_j)} + e^{-(g + \tau_i + \tau_j)}\right] y \, dy

The self-collision term follows the same pattern:

.. math::
   :label: self-sph

   r_{ii} = 2\Sigt{i} \int_0^{R_{\text{cell}}}
     \left[2\ell_i(y) - \frac{2}{\Sigt{i}}
       \left(1 - e^{-\tau_i(y)}\right)\right] y \, dy

The **same code path** :meth:`CPMesh._compute_radial_rcp` handles both
cylindrical and spherical, parameterised by kernel and weights.
Verified against ``orpheus/derivations/continuous/flat_source_cp/sphere.py::_sphere_cp_matrix``.


Geometry Comparison
====================

.. list-table::
   :header-rows: 1
   :widths: 18 27 27 28

   * - Aspect
     - Slab
     - Cylinder
     - Sphere
   * - Kernel :math:`F(\tau)`
     - :math:`E_3` (:func:`_e3`)
     - :math:`\text{Ki}_4` (:func:`_ki4_lookup`)
     - :math:`e^{-\tau}`
   * - :math:`F(0)`
     - :math:`1/2`
     - :math:`\approx 0.4244`
     - :math:`1`
   * - Quadrature weight
     - none (scalar)
     - :math:`w_i`
     - :math:`w_i \cdot y_i`
   * - Volume :math:`V_i`
     - :math:`t_i`
     - :math:`\pi(R_i^2 - R_{i-1}^2)`
     - :math:`\tfrac{4}{3}\pi(R_i^3 - R_{i-1}^3)`
   * - Angular integration
     - analytical (in :math:`E_3`)
     - numerical (:math:`y`-quadrature)
     - numerical (:math:`y`-quadrature)
   * - Surface :math:`S`
     - :math:`1`
     - :math:`2\pi R_{\text{cell}}`
     - :math:`4\pi R_{\text{cell}}^2`
   * - Prefactor
     - 1/2 (half-space)
     - 2 (chord halves)
     - 2 (chord halves)
   * - Code path (solver)
     - :meth:`CPMesh._compute_slab_rcp`
     - :meth:`CPMesh._compute_radial_rcp`
     - :meth:`CPMesh._compute_radial_rcp`
   * - Derivation script
     - ``orpheus/derivations/continuous/flat_source_cp/slab.py``
     - ``orpheus/derivations/continuous/flat_source_cp/cylinder.py``
     - ``orpheus/derivations/continuous/flat_source_cp/sphere.py``


The Eigenvalue Problem
=======================

Multi-Group Formulation
-----------------------

For :math:`G` energy groups, the neutron balance in region :math:`i`,
group :math:`g` is:

.. math::
   :label: neutron-balance

   \Sigt{ig} V_i \phi_{ig}
   = \sum_j P_{ji,g}^{\infty} V_j \left[
       \frac{\chi_{jg}}{\keff} \sum_{g'} \nSigf{jg'} \phi_{jg'}
       + \sum_{g'} \bigl(\Sigs{j,g' \to g}
       + 2\,\Sigma_{2,j,g' \to g}\bigr) \phi_{jg'}
     \right]

The factor of 2 on :math:`\Sigma_2` accounts for both outgoing (n,2n)
neutrons (the original neutron is already removed by :math:`\Sigt{}`).

.. note::

   **Scattering convention.**  ``Mixture.SigS[l]`` uses
   ``SigS[g_from, g_to]`` convention.  The in-scatter source uses the
   transpose: ``Q_scatter = SigS.T @ phi``.  This is critical for
   correct downscatter/upscatter coupling.


Matrix Form
-------------

:eq:`neutron-balance` in matrix form:
:math:`\mathbf{A}\Phi = (1/\keff)\mathbf{B}\Phi` with:

.. math::
   :label: matrix-A-def

   A[(i,g),(j,g')] = \delta_{ij}\delta_{gg'}\Sigt{ig} V_i
     - P_{ji,g}^{\infty} V_j (\Sigs{j,g' \to g} + 2\Sigma_{2,j,g' \to g})

.. math::
   :label: matrix-B-def

   B[(i,g),(j,g')] = P_{ji,g}^{\infty} V_j \chi_{jg}\, \nSigf{jg'}

The analytical verification eigenvalue is
:math:`\lambda_{\max}(\mathbf{A}^{-1}\mathbf{B})`, computed by
:func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp`, which builds the full
:math:`NG \times NG` matrices and uses ``numpy.linalg.eigvals``.
The solver does NOT form these matrices --- see :ref:`why-not-full-matrices`.


Power Iteration
---------------

Implemented by :func:`~numerics.eigenvalue.power_iteration` via the
:class:`CPSolver` protocol:

1. **Fission source** (:meth:`CPSolver.compute_fission_source`):
   :math:`Q^F_{ig} = \chi_{ig} \sum_{g'} \nSigf{ig'} \phi_{ig'} / \keff`

2. **Fixed-source solve** (:meth:`CPSolver.solve_fixed_source`):
   :math:`\phi_{ig}^{\text{new}} = \sum_j P_{ji,g}^{\infty} V_j Q_{jg}^{\text{total}} / (\Sigt{ig} V_i)`

3. **Update** :math:`\keff` (:meth:`CPSolver.compute_keff`):

   .. math::
      :label: cp-keff-update

      \keff = \frac{\sum_{ig} \nSigf{ig} \phi_{ig} V_i}
                   {\sum_{ig} (\Sigt{ig} - \Sigs{ig}^{\text{out}}
                   - 2\Sigma_{2,ig}^{\text{out}}) \phi_{ig} V_i}

   The denominator is the **net removal rate**: total interactions minus
   all neutrons that remain in the system (scattered or produced by (n,2n)).

   .. admonition:: The (n,2n) subtlety (ERR-015)

      **Bug:** ``compute_keff`` used production/absorption, which is wrong
      when :math:`\Sigma_2 \neq 0` because absorption includes (n,2n)
      removal but does not credit the two outgoing neutrons.  With
      ``Sig2[0,0] = 0.01`` on region A (2G), the solver converged to
      :math:`k = 1.793` instead of analytical :math:`k = 2.045` (12% error).

      **First wrong fix:** Adding (n,2n) production to the numerator
      gave :math:`k = 1.808` --- still wrong because the extra neutron
      has no fission spectrum :math:`\chi` weighting.

      **Correct fix:** Net-removal denominator :eq:`cp-keff-update`.
      When :math:`\Sigma_2 = 0`, this reduces to
      :math:`\nSigf{}\phi V / \Siga{}\phi V`.

      **How it hid:** All test materials had ``Sig2 = 0`` (zero sparse
      matrix).  The ``make_mixture`` API didn't even accept a ``sig_2``
      parameter, making it impossible to construct test materials with
      nonzero (n,2n).

      **Test:** ``tests/cp/test_verification.py::TestN2N::test_n2n_solver_keff_matches_analytical``.

      **Lesson:** When adding a new reaction type, trace it through BOTH
      the transport solve AND the eigenvalue estimate.  Test with the
      term nonzero.

   For lattice models with reflective (or white) boundary conditions,
   the leakage term is zero, so this balance is exact.

   .. note:: **Discipline linkage (#259 P1, 2026-07-03).**

      CP was already the **operator-consistent member** when the unified
      k-estimator discipline landed.  :eq:`cp-keff-update` is fission-only
      production over net removal,
      :math:`\Sigma_t - \Sigma_s^{\text{out}} - 2\Sigma_2^{\text{out}}
      \equiv \Sigma_a + L - E_{2n}`, with the leakage :math:`L`
      structurally zero (:math:`P^{\infty}` is the infinite-lattice
      kernel) and the :math:`(n,2n)` emission on the removal side —
      exactly the eigenvalue of the map the inner solve poses (only
      fission scaled by :math:`1/k`).  The same discipline was
      retrofitted onto SN (the #291 boundary-leakage term) and MoC (the
      R7 numerator/denominator flip); see :ref:`sn-keff-estimator` in
      :doc:`/theory/methods/sn/solver`.  Only the **substrate** still differs —
      CP contracts group-last ``CellXS`` arrays where SN uses the typed
      :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`
      fields — and unifying that is the #259 close-out follow-up.

4. **Converge** (:meth:`CPSolver.converged`) when :math:`|\Delta k| <`
   ``keff_tol`` and :math:`\|\Delta\phi\|_\infty <` ``flux_tol``.

The power iteration algorithm is shared with all ORPHEUS eigenvalue
solvers; see :ref:`power-iteration-algorithm` in the homogeneous theory
for the general formulation.

The power iteration converges to the **dominant eigenvalue**
(:math:`\keff`) and the **fundamental mode** --- the unique non-negative
eigenvector (Perron--Frobenius theorem) :cite:`Hebert2009`.

This is implemented in :class:`CPSolver`, which satisfies the
:class:`~numerics.eigenvalue.EigenvalueSolver` protocol.


Solver Modes: Jacobi and Gauss-Seidel
---------------------------------------

:class:`CPSolver` supports two modes via ``CPParams.solver_mode``.

**Jacobi (default,** :meth:`CPSolver._solve_fixed_source_jacobi` **).**
Scattering source from all groups computed simultaneously using the
previous iteration's flux.  No inner iterations.  One matrix-vector
multiply per group per outer.

**Gauss-Seidel (**\ :meth:`CPSolver._solve_fixed_source_gs` **).**
Groups swept from fast (:math:`g=0`) to thermal (:math:`g=G-1`).  For
each group :math:`g`, **inner iterations** converge within-group
self-scatter before proceeding to :math:`g+1`.

The inner iteration solves the fixed-point equation:

.. math::
   :label: cp-within-group-fixed-point

   \phi_g = T_g\!\bigl[Q_g^{\text{ext}} + \Sigs{g \to g} \phi_g\bigr]

where :math:`T_g[\cdot]` is the CP transport operator for group :math:`g`
(the matrix multiply :math:`P^T V \cdot / (\Sigt{} V)`).  The operator is
exact for a given source, but the source depends on :math:`\phi_g` through
within-group self-scatter :math:`\Sigs{g \to g}`.

**Why inner iterations matter.**  The spectral radius of the inner
iteration operator is approximately :math:`\Sigs{g \to g} / \Sigt{g}`.
Thermal groups (ratio ~0.6--0.9) need 3--8 inner iterations; fast groups
(ratio << 1) converge in 1.  Convergence criterion: relative flux change
:math:`\|\phi_g^{(n+1)} - \phi_g^{(n)}\| / \|\phi_g^{(n+1)}\| < \texttt{inner\_tol}`.

.. admonition:: Tautological inner residual (ERR-016)

   **Bug:** The original inner convergence check computed
   :math:`\|\Sigt{} V \phi_g^{\text{new}} - P^T V Q_g\|`, which is
   **identically zero by construction**:
   :math:`\phi_g^{\text{new}} \equiv P^T V Q_g / (\Sigt{} V)`, so
   the check computes :math:`\|x - x\| = 0`.  The inner loop always
   exited after 1 iteration, making GS functionally identical to
   sequential-group Jacobi.

   **How it hid:** All 27 eigenvalue tests passed (the outer iteration
   converged regardless).  The ``n_inner`` array showed all 1s,
   interpreted as "fast convergence" rather than "broken check".  A QA
   review initially concluded that inner iterations are *fundamentally
   unnecessary* for CP --- this is **wrong**: the transport is exact for
   a given source, but the source depends on :math:`\phi_g` through
   self-scatter.

   **Fix:** Changed to relative flux change
   :math:`\|\phi^{(n+1)} - \phi^{(n)}\| / \|\phi^{(n+1)}\|`.  With the
   corrected residual, thermal groups genuinely require multiple inner
   iterations; fast groups converge in 1.

   **Tests:** ``tests/cp/test_verification.py::TestGSInnerIterations`` ---
   ``test_thermal_needs_more_inner_than_fast`` (thermal > fast inner
   counts), ``test_gs_eigenvalue_matches_jacobi`` (same eigenvalue),
   ``test_no_self_scatter_one_inner`` (zero diagonal in :math:`\Sigma_s`
   converges in 1).

   **Lesson:** A convergence check that compares quantities derived from
   each other by construction tests nothing.  Always verify with a
   problem that *should* require multiple iterations.

**Convergence diagnostics** (stored in :class:`CPResult`):

.. list-table::
   :header-rows: 1
   :widths: 25 15 40

   * - Field
     - Shape
     - Description
   * - ``residual_history``
     - ``(n_outer,)``
     - Neutron balance residual :math:`\|\Sigt{} V \phi - P^T V Q\|_2` per outer (both modes)
   * - ``n_inner``
     - ``(n_outer, ng)``
     - Inner iteration count per group per outer (GS only; ``None`` in Jacobi)

Both modes converge to the **same eigenvalue and flux distribution**
(``tests/cp/test_verification.py::TestGSInnerIterations::test_gs_eigenvalue_matches_jacobi``).

.. list-table:: Solver mode comparison
   :header-rows: 1
   :widths: 35 35

   * - Jacobi
     - Gauss-Seidel
   * - All groups from previous flux
     - Sequential; latest flux from earlier groups
   * - No inner iterations
     - Inner iterations for within-group self-scatter
   * - Simpler (1 matrix multiply/group/outer)
     - Faster convergence for strong upscatter


Cross-Section Data Layout
===========================

Cross sections flow through two layers:

1. :class:`~data.macro_xs.mixture.Mixture` --- per-material: ``(ng,)``
   arrays for :math:`\Sigt{}`, :math:`\Siga{}`, :math:`\Sigf{}`,
   :math:`\nSigf{}`, :math:`\chi`; ``(ng, ng)`` sparse matrices for
   scattering :math:`\Sigs{}` (``SigS[0]``) and (n,2n) :math:`\Sigma_2`
   (``Sig2``).

2. :class:`~data.macro_xs.cell_xs.CellXS` --- per-cell on the mesh:
   ``(N_cells, ng)`` arrays mapped from materials by
   :func:`~data.macro_xs.cell_xs.assemble_cell_xs` using ``mat_ids``.

**Indexing:** ``xs.sig_t[i, g]`` is the total macroscopic cross section
in spatial cell :math:`i` and energy group :math:`g`.  Group 0 = fastest
(highest energy); group :math:`G-1` = slowest (thermal) --- the
:ref:`canonical fast-first energy-group convention
<canonical-group-convention>`.  Cell 0 = innermost; cell :math:`N-1` =
outermost.

**Scattering convention:** ``SigS[g_from, g_to]``.  The in-scatter source
uses the transpose: ``Q += SigS.T @ phi`` (applied per-cell)::

    for k in range(N):
        mid = self.mat_ids[k]
        Q[k, :] += self._scat_mats[mid].T @ flux_distribution[k, :]
        Q[k, :] += 2.0 * (self._n2n_mats[mid].T @ flux_distribution[k, :])


Verification
============

**106 tests** across 6 test files verify the CP implementation against
analytical eigenvalues computed independently by the derivation modules.

Consolidated Eigenvalue Solvers (CP-20260405-007)
---------------------------------------------------

:func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp` and
:func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous` replaced 10 duplicated
eigenvalue computations across ``homogeneous.py``, ``sn.py``, ``moc.py``,
``mc.py``, ``cp_slab.py``, ``cp_cylinder.py``, and ``cp_sphere.py``.
Both accept optional ``sig_2`` / ``sig_2_mats`` for (n,2n) reactions.
When omitted, (n,2n) is zero and all previous eigenvalues are preserved.


Eigenvalue Verification Cases
-------------------------------

27 core cases: {1, 2, 4} groups |times| {1, 2, 4} regions |times|
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
----------------------------------------

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


Implementation Details
=======================

.. _why-p-transpose:

Why P\ :sup:`T` Appears in the Flux Update
--------------------------------------------

With :math:`P_{ij} = P[\text{birth}_i, \text{collision}_j]`, the neutron
balance :eq:`collision-rate` sums over birth regions (first index of
:math:`P`), which is a column sum.  In matrix form:

.. math::
   :label: p-transpose-flux-balance

   \Sigt{} V \phi = P^T V Q

Hence: :math:`\phi = P^T V Q / (\Sigt{} V)`.

In code: ``phi[:, g] = P_inf[:, :, g].T @ source``
(:meth:`CPSolver._solve_fixed_source_jacobi`).

.. admonition:: Historical bug (ERR-009)

   ``P @ source`` instead of ``P.T @ source``.  Correct for homogeneous
   problems (:math:`P` symmetric) but **8% wrong** for the 1G 2-region
   slab (k=1.373 vs analytical 1.272).  Caught by the formal verification
   suite on the first heterogeneous test case.  The bug survived 4 weeks
   because all prior tests were homogeneous.


.. _ki-table-construction:

Ki\ :sub:`3` Kernel Construction (historical: Ki\ :sub:`3`/Ki\ :sub:`4` tables)
-------------------------------------------------------------------------------

**Current path (post-Phase B.4, commit 6badbe5).** Both the
derivation module and the runtime solver consume the **same**
kernel via :func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp`: a Chebyshev
polynomial of degree 63 fit to the scaled kernel
:math:`e^{\tau}\,\mathrm{Ki}_3(\tau)` (canonical A&S) on
:math:`[0, 50]` at Chebyshev-Gauss-Lobatto nodes. The tabulation
values are computed once at module load by
:func:`~orpheus.derivations.common.kernels.ki_n_mp` at 30 dps and cached via
:func:`functools.lru_cache`; Chebyshev fitting of the
:math:`e^{\tau}`-scaled kernel converts the exponentially-decaying
tail into a slowly-varying function the polynomial reaches to
~:math:`5\times 10^{-6}` absolute accuracy. Evaluation at runtime
is one :class:`numpy.polynomial.Chebyshev` call plus one
:func:`numpy.exp` — comparable in cost to the legacy
:func:`numpy.interp` on a 20 000-point grid. For
:math:`\tau > 50`, the clamp to zero is already below double
precision (:math:`\mathrm{Ki}_3(50) \approx 3\times 10^{-23}`).

.. note::

   :func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` uses the **canonical
   A&S convention**:
   :math:`\mathrm{Ki}_n^{\text{A\&S}}(x) = \int_0^{\pi/2}
   \cos^{n-1}\theta\,\exp(-x/\cos\theta)\,d\theta`, so
   :math:`\mathrm{Ki}_3^{\text{A\&S}}(0) = \pi/4`. The legacy
   ORPHEUS convention carried by :eq:`ki3-def` used an extra power
   of :math:`\sin\theta`, making the ORPHEUS-local "Ki\ :sub:`4`\ "
   the canonical :math:`\mathrm{Ki}_3^{\text{A\&S}}`. The
   :math:`P_{ij}` matrix assembled by
   :mod:`orpheus.derivations.continuous.flat_source_cp.cylinder` consumes canonical
   :math:`\mathrm{Ki}_3^{\text{A\&S}}` — the Phase-4.2 alias
   ``BickleyTables.Ki3_vec := BickleyTables.ki4_vec`` made the
   convention uniform before the Phase B.4 retirement deleted both
   names.

**Historical path (retired, pre-Phase B.4).** Until commit
``6badbe5``, the kernel tabulation lived in two places that had to
stay in sync:

1. :func:`_build_ki_tables` on :class:`orpheus.cp.solver.CPMesh`
   (solver side) and
2. The ``BickleyTables`` class in :mod:`orpheus.derivations.common.kernels`
   (derivation side).

Both tabulated the functions identically:

1. :math:`\text{Ki}_3(x_k) = \int_0^{\pi/2} e^{-x_k/\sin\theta} \sin\theta\,d\theta`
   (ORPHEUS-local convention, equivalent to canonical
   :math:`\mathrm{Ki}_2^{\text{A\&S}}`) at 20 000 points on
   :math:`[0, 50]` via :func:`scipy.integrate.quad`. Boundary
   value :math:`\text{Ki}_3(0) = 1` set analytically.

2. :math:`\text{Ki}_4` (ORPHEUS-local) via cumulative trapezoid
   from right to left (a pre-computer-era idiom for the recurrence
   :math:`\mathrm{Ki}_{n+1}(x) = \int_x^\infty \mathrm{Ki}_n(t)\,dt`)::

       ki4_vals = np.cumsum(ki3[::-1])[::-1] * dx

   This realised canonical :math:`\mathrm{Ki}_3^{\text{A\&S}}` to
   :math:`O(\Delta x^2)` trapezoidal error on a uniform grid.

3. Runtime lookup via :func:`numpy.interp` (linear interpolation);
   with :math:`\Delta x = 0.0025`, interpolation error was
   :math:`O(\Delta x^2) \approx 6 \times 10^{-6}` and the
   cumulative-trapezoid bias on the antiderivative pushed the
   effective kernel accuracy to ~:math:`10^{-3}` at :math:`x = 0`
   tapering to ~:math:`10^{-4}` elsewhere — the dominant error in
   every cylindrical flat-source CP reference value for the life
   of the project.

The old ``tests/cp/test_verification.py::TestKi4Resolution`` confirmed
that increasing from 5 000 to 40 000 points produced diminishing
returns (validating the legacy 20 000 default). It was replaced in
Phase B.4 by
``test_cylinder_ki3_kernel_is_insensitive_to_n_ki_table`` — a much
stronger statement than "converges with more points": the Chebyshev
interpolant is built lazily from mpmath and is ignorant of any
``CPParams.n_ki_table`` knob, so ``keff`` is bit-identical across
``{5000, 20000, 40000}`` table-size values. The knob itself is
kept only as an unused no-op argument so that callers constructing
:class:`~orpheus.cp.solver.CPParams` with explicit
``n_ki_table=...`` do not break.

Why the retirement was deferrable until Phase B.4 is covered in
:doc:`/theory/references/peierls_nystrom` §16 (postmortem) and
:doc:`/verification/reference_solutions` §"Legacy naming
discrepancy in ``BickleyTables``" (safety argument for the swap).


Equal-Volume Mesh Subdivision
-------------------------------

:func:`~geometry.factories._subdivide_zone` creates equal-volume cells:

**Cartesian:** :math:`x_k = x_0 + k(x_N - x_0)/N`.

**Cylindrical.**  Equal annular volumes:
:math:`V = \pi(R_k^2 - R_{k-1}^2) = \text{const}`.  Summing:
:math:`\pi R_k^2 = \pi R_0^2 + k \cdot \pi(R_N^2 - R_0^2)/N`, giving:

.. math::
   :label: equal-volume-radius-cylindrical

   R_k = \sqrt{R_0^2 + \frac{k}{N}(R_N^2 - R_0^2)}

For :math:`R_0 = 0`: :math:`R_k = R_N\sqrt{k/N}`.

**Spherical.**  Equal shell volumes:
:math:`V = \frac{4}{3}\pi(R_k^3 - R_{k-1}^3) = \text{const}`, giving:

.. math::
   :label: equal-volume-radius-spherical

   R_k = \left(R_0^3 + \frac{k}{N}(R_N^3 - R_0^3)\right)^{1/3}

For :math:`R_0 = 0`: :math:`R_k = R_N(k/N)^{1/3}`.


.. _white-bc-quality:

Numerical Evidence: White-BC Approximation Quality
-----------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 20 20 15

   * - Geometry
     - CP (white)
     - SN (reflective)
     - Difference
   * - Slab 2-region (fuel + mod)
     - 1.272
     - 1.261
     - +0.9%
   * - Cylinder 2-region (10 cells/zone)
     - 0.995
     - 0.997
     - -0.2%

The slab shows ~1% overestimation because the white BC smears the
anisotropic angular flux at the fuel-moderator interface.  The
cylindrical discrepancy is smaller because the Wigner--Seitz cell has
a larger volume-to-surface ratio, making the boundary angular
distribution closer to isotropic.  The SN solver uses reflective
(specular) BCs; see :ref:`boundary-conditions`.  The ~1% gap between
white and reflective BCs is expected physics, not a bug.


.. _why-not-full-matrices:

Design Decisions
==================

**Single CPMesh class (CP-20260404-001).**  All three geometries share
the eigenvalue iteration.  Only the kernel construction differs.  Adding
a geometry requires one ``_setup_*()`` method and kernel function, not a
new solver class.

**Per-group loop (CP-20260404-008).**  The CP matrix
:math:`P^{\infty}_{ij,g}` is computed independently per group because
:math:`\tau_{ig} = \Sigt{ig} \ell_i` varies with group.  Vectorising
over groups would require restructuring the Ki4 table lookup ---
deferred as a performance improvement.

**Why ORPHEUS does NOT form the full A/B matrices.**  The analytical
verification (:func:`~orpheus.derivations.common.eigenvalue.kinf_from_cp`) builds the
full :math:`NG \times NG` matrices and runs ``numpy.linalg.eigvals``.
The solver does not, for two reasons:

1. **Memory.** For 421 groups, 20 cells: :math:`8420 \times 8420 =`
   567 MB dense.  The per-group approach: :math:`20 \times 20 \times 421 =`
   1.3 MB.
2. **Sparsity.** The scattering matrix is sparse (downscatter only for
   fast groups).  Per-material sparse storage is far more efficient.

**1-group verification is degenerate.** :math:`\keff = \nSigf{}/\Siga{}`
regardless of flux shape.  Wrong CP matrices, weight errors, and
convention drifts are invisible.  The suite tests 1, 2, AND 4 groups.

.. note::

   The 1-group degenerate case masked ERR-001 (z-ordinate weight loss),
   ERR-002 (scattering transpose), ERR-004 (BiCGSTAB normalisation),
   ERR-009 (CP transpose), and ERR-015 ((n,2n) keff) --- all caught
   only by multi-group tests.  See ``.claude/skills/vv-principles/error_catalog.md``.


.. |times| unicode:: U+00D7
