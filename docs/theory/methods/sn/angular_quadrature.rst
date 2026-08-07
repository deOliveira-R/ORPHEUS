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
1D mesh, z-directed :term:`ordinates <ordinate>` (:math:`\mu_x = \mu_y = 0`) are handled
as pure collision; all others stream along *x* with *y*-terms cancelling
via reflective BCs.

Built by :meth:`Quadrature.lebedev(order)
<orpheus.numerics.quadrature.Quadrature.lebedev>`.

Level-Symmetric S\ :sub:`N`
----------------------------

Standard triangular :term:`quadrature` with :math:`N/2` distinct :math:`\mu_z`
values per hemisphere.  Ordinates on each level are permutations of the
direction cosine set satisfying :math:`\eta^2 + \xi^2 + \mu^2 = 1`.
Equal spacing in :math:`\mu^2` is used with :math:`\mu_1^2 = 4/(N(N+2))`
:cite:`CarlsonLathrop1965`.

.. note::

   ⚠ **The** :math:`\mu_1` **seed's attribution is pending verification.**  The
   *recursion* around it — :math:`\mu_i^2 = \mu_1^2 + (i-1)\,C` with
   :math:`C = 2(1 - 3\mu_1^2)/(N-2)` — is the standard level-symmetric
   construction and is not in doubt.  Whether the particular seed
   :math:`\mu_1^2 = 4/(N(N+2))` (giving :math:`\mu_1 = 0.4082` at
   :math:`S_4`, :math:`0.1179` at :math:`S_{16}`) is the one
   :cite:`CarlsonLathrop1965` tabulates has **not** been checked against the
   source.  Until it is, read the citation as covering the construction, not
   the seed.  Issue **#327**.

Weights sum to :math:`4\pi`.  Provides the ``level_indices`` structure
needed by the cylindrical :term:`sweep`.  Unlike the product quadrature
(which has one level per :math:`\mu_z` value), the Level-Symmetric
quadrature groups both :math:`+\mu_z` and :math:`-\mu_z` hemispheres
on the same level (grouped by :math:`|\mu_z|`).  Within each level,
ordinates are sorted by increasing :math:`\eta` for the azimuthal sweep.

Built by :meth:`Quadrature.level_symmetric(sn_order)
<orpheus.numerics.quadrature.Quadrature.level_symmetric>`.

.. _quadrature-ls-weights:

The weights: one free parameter per :math:`O_h` orbit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The node set is :math:`O_h`-invariant — closed under all 48 signed coordinate
permutations.  Invariance **forces the weight to be constant on each orbit**,
so the free parameters are exactly one per orbit, no more and no fewer.  An
orbit is one multiset :math:`\{|x|, |y|, |z|\}`; since every cosine is drawn
from the level array, the orbit label is the sorted triple of *level indices* —
an integer fact of the construction, not a floating-point classification made
afterwards.

Those free weights are fixed by the lowest independent moment conditions
against the closed form

.. math::
   :label: quadrature-sphere-monomial

   \int_{S^2} x^a y^b z^c \,\mathrm{d}\Omega =
   \begin{cases}
     2\,\dfrac{\Gamma\!\left(\frac{a+1}{2}\right)
               \Gamma\!\left(\frac{b+1}{2}\right)
               \Gamma\!\left(\frac{c+1}{2}\right)}
              {\Gamma\!\left(\frac{a+b+c+3}{2}\right)}
       & a, b, c \text{ all even} \\[2mm]
     0 & \text{otherwise.}
   \end{cases}

Only even triples constrain anything (sign-flip closure kills the odd ones),
and permuted triples give the *same* equation, so the independent conditions
are the even triples up to permutation.  ``[M]`` the resulting degrees:

.. list-table::
   :header-rows: 1
   :widths: 8 10 12 14

   * - :math:`N`
     - orbits
     - degree
     - min weight
   * - 2
     - 1
     - 3
     - 1.570796
   * - 4
     - 1
     - 3
     - 0.523599
   * - 6
     - 2
     - 5
     - 0.201682
   * - 8
     - 3
     - 7
     - 0.131132
   * - 10
     - 4
     - 9
     - 0.057100
   * - 12
     - 5
     - 11
     - 0.027825

⭐ **At** :math:`S_2` **and** :math:`S_4` **the solve is provably a no-op.**
Both node sets are a *single* orbit, so :math:`O_h`-invariance plus
:math:`\sum w = 4\pi` determine the one weight uniquely — the equal-weight
value was already correct, and degree 3 is the ceiling any rule on those nodes
can reach.  This is why :math:`S_4`, the workhorse order, is bit-identical
across the #327 fix.

.. _quadrature-ls-positivity:

Why the family stops at :math:`S_{12}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``[M]`` the moment-matched solve yields a **negative** weight from
:math:`S_{14}` upward on these levels (:math:`-0.027` at :math:`S_{14}`,
:math:`-0.018` at :math:`S_{16}`, :math:`-0.142` at :math:`S_{20}`).
:meth:`Quadrature.level_symmetric
<orpheus.numerics.quadrature.Quadrature.level_symmetric>` **refuses** there.

Positivity is not a preference to trade against accuracy.  The scalar flux is
:math:`\phi = \sum_n w_n \psi_n`, so a negative weight lets a *non-negative*
angular flux integrate to a *negative* scalar flux, and the boundary response
kernels assert non-negativity outright (``assert_response_positive_if_declared``;
the Lambertian's sub-Markov bound).  A rule that violates it is not a transport
quadrature, so the family is defined exactly where it is valid rather than
shipping a silent second construction under one name.

The frontier is **computed, not tabulated** — the builder reads the sign of its
own solution — so it tracks the node set instead of going stale beside it.  If
a different :math:`\mu_1` seed pushes it up, the refusal moves with it and
:func:`~orpheus.numerics.quadrature.registry.select_quadrature` follows,
because its inverter **asks the family** (attempts the construction, returns
the documented ``None`` on refusal) rather than consulting a second copy of the
limit that would drift.  Above :math:`S_{12}`, use Lebedev or the product rule.

.. admonition:: Development history — #327, 2026-08-06
   :class: note

   Until 2026-08-06 the weights were **one equal value for every ordinate**
   (:math:`w = 4\pi/8n_{\rm octant}`; ``[M]`` exactly one distinct weight at
   every order), while the rule advertised degree :math:`N-1`.  ``[M]`` its
   measured degree was **3 at every order** — an over-claim of 12 at
   :math:`S_{16}`.  The degree-3 it did reach was free: *any*
   :math:`O_h`-symmetric node set with :math:`\sum w = 4\pi` reaches 3, so the
   weights were earning none of it.

   The consequence was live but unreached: the discrete spherical-harmonic Gram
   was 9.7 % wrong at :math:`\ell = 2` for :math:`S_8`, 18–25 % for
   :math:`S_{16}`, and **worst of all — 37.5 % — at** :math:`S_4`, the most-used
   order.  It got *worse* with :math:`N`, because refining an order-independent
   deficit converges nothing.  :math:`\ell = 0` and :math:`\ell = 1` are exact,
   which is why isotropic and P1 transport looked healthy and no gate saw it.

   ⭐ **Why no test caught it.**  Every exactness test stopped at degree 2 —
   inside the sector :math:`O_h` symmetry makes free — so they passed on any
   orbit set with the right total mass.  The remaining tests pinned the *tag*,
   asserting the advertised number equals :math:`N-1`, which it did.  The
   measured functional's invariance group contained the error class
   (:ref:`verification-anti-patterns`, Mode 12).

   ⭐ **The tag was wrong in BOTH directions**, and that is what identified it.
   :math:`S_2` *under*-claimed (advertised 1, measured 3) while :math:`S_6`
   upward over-claimed.  A number that is wrong both ways is not a
   mis-calibrated measurement of this rule; it is a **formula describing a
   different one** — :math:`N-1` is the moment-matched construction's degree,
   and the implementation was not moment-matched.  The two coincided only at
   :math:`N = 4`.  The gate that found it therefore asserts *two* claims,
   ``measured >= advertised`` (safety) and ``measured == advertised``
   (accuracy), because a bare equality reddens identically and hides the
   direction.

   Gated by ``tests/numerics/test_advertised_degree_is_measured.py``, which
   measures every production family against
   :eq:`quadrature-sphere-monomial` with the other three rules swept by the
   same body as controls.

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

Ordinate Permutations
---------------------

:meth:`Quadrature.ordinate_permutation(motion)
<orpheus.numerics.quadrature.Quadrature.ordinate_permutation>` answers
"which permutation does a rigid motion induce on this rule's weighted
ordinates?" — the single source read by the boundary realizer's deck
kernel, the specular certifications, and the curvilinear :math:`r = 0`
pole seed.  Writing :math:`Q` for the motion's **linear part**
(ordinates are DIRECTIONS — a translation does not act on them, so a
periodic wrap induces the identity), the returned permutation
:math:`\pi` is defined by

.. math::
   :label: quadrature-ordinate-permutation

   \Omega_{\pi(n)} = Q\,\Omega_n
   \qquad\text{and}\qquad
   w_{\pi(n)} = w_n
   \qquad\text{for every } n,

**certified at derivation**: every image must MATCH a node inside the
certification window (no bare nearest-neighbour), the match must be a
**bijection**, and the matched weights must be equal.  A motion that
fails any clause answers ``None`` — the honest "this motion is not a
symmetry of this rule"; each consumer owns its domain-specific refusal
("no specular pairing" at the BC certification tier, "cannot seed the
r = 0 pole" at the curvilinear sweep tier).

For an axis mirror :math:`\sigma_a` the induced permutation is the
classical *reflection partner* map used by
:term:`reflective boundary conditions <reflective boundary condition>`
(see :ref:`boundary-conditions`).  For Gauss--Legendre (1D, nodes
embedded as :math:`(\mu, 0, 0)`) the x-mirror gives
:math:`\pi(n) = N - 1 - n` by GL-node symmetry and the y/z mirrors fix
every ordinate.  An odd-:math:`n_\varphi` product rule has NO x-mirror
permutation: its mirror planes sit at :math:`k\pi/n_\varphi`, and
:math:`\sigma_x` needs the :math:`k = n_\varphi/2` plane — an integer
only for even :math:`n_\varphi` — so the derivation answers ``None``
there.

.. note::

   **History (ERR-074).** Until G6.3 §7d (2026-08) the three axis-mirror
   maps were a precomputed table (``reflection_index`` /
   ``reflection_partners``), and until campaign Q5.0.1 (2026-08-02) that
   table was built by a bare per-node ``argmin`` — no distance window, no
   injectivity check, no weight comparison — which on a
   :math:`\sigma_x`-unclosed rule silently returned the nearest WRONG
   ordinate: at ``product(4, 5/7/9)`` the shipped axis-0 map was off by
   ``0.58 / 0.42 / 0.33`` in the direction cosines *and still an
   involution*, so a self-inverse check could not see it.  The
   certification above is exactly what that defect class demanded; the
   general motion-derived method then retired the table once every
   consumer read it.

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


