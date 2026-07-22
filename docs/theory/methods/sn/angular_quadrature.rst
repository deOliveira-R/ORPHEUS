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

Built by :meth:`Quadrature.gauss_legendre(n)
<orpheus.numerics.quadrature.Quadrature.gauss_legendre>`.

Lebedev (Sphere)
-----------------

For 2D/3D Cartesian geometry: :math:`N` points on the unit sphere with
octahedral symmetry :cite:`Lebedev1999`.  Weights sum to :math:`4\pi`.  On a
1D mesh, z-directed ordinates (:math:`\mu_x = \mu_y = 0`) are handled
as pure collision; all others stream along *x* with *y*-terms cancelling
via reflective BCs.

Built by :meth:`Quadrature.lebedev(order)
<orpheus.numerics.quadrature.Quadrature.lebedev>`.

Level-Symmetric S\ :sub:`N`
----------------------------

Standard triangular quadrature with :math:`N/2` distinct :math:`\mu_z`
values per hemisphere.  Ordinates on each level are permutations of the
direction cosine set satisfying :math:`\eta^2 + \xi^2 + \mu^2 = 1`.
Equal spacing in :math:`\mu^2` is used with :math:`\mu_1^2 = 4/(N(N+2))`
:cite:`CarlsonLathrop1965`.

Weights sum to :math:`4\pi`.  Provides the ``level_indices`` structure
needed by the cylindrical sweep.  Unlike the product quadrature
(which has one level per :math:`\mu_z` value), the Level-Symmetric
quadrature groups both :math:`+\mu_z` and :math:`-\mu_z` hemispheres
on the same level (grouped by :math:`|\mu_z|`).  Within each level,
ordinates are sorted by increasing :math:`\eta` for the azimuthal sweep.

Built by :meth:`Quadrature.level_symmetric(sn_order)
<orpheus.numerics.quadrature.Quadrature.level_symmetric>`.

Product Quadrature (GL x equispaced)
-------------------------------------

Tensor product of Gauss--Legendre in :math:`\mu = \cos\theta` (polar)
and equispaced points in :math:`\varphi` (azimuthal).  Each :math:`\mu`
level has the same number of azimuthal points, giving a clean level
structure ideal for the cylindrical sweep.  Weights:

.. math::
   :label: quadrature-product-weights

   w_{p,m} = w_{\text{GL}}(\mu_p) \cdot \frac{2\pi}{N_\varphi}

Sum to :math:`4\pi`.  Within each level, ordinates are sorted by
increasing :math:`\eta = \sin\theta\cos\varphi` to match the
:math:`\alpha` recursion convention from :cite:`BaileyMorelChang2010` (the
curvilinear :math:`\alpha`-recursion).

Built by :meth:`Quadrature.product(n_mu, n_phi)
<orpheus.numerics.quadrature.Quadrature.product>`.

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
   :label: quadrature-reflection-index

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


