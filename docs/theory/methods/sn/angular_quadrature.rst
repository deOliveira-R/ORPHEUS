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
Equal spacing in :math:`\mu^2` is the Carlson–Lathrop recursion
:math:`\mu_i^2 = \mu_1^2 + (i-1)\,\Delta` with
:math:`\Delta = 2(1 - 3\mu_1^2)/(N-2)` — Eq. (3-52), printed p. 32 of
:cite:`CarlsonLathrop1965` — seeded with the **moment-matched**
:math:`\mu_1` (issue #337, 2026-08-08): the smallest root of *the rule
integrates* :math:`\mu_z^N` *exactly*, root-found by the builder at
construction.  This is the construction behind the published LA-3186
tables, and it reproduces their :math:`\mu_1` to every printed digit at
each tabulated order.

.. note::

   **The seed's provenance — verified against the primary sources and
   then adopted** (issues #327 → #337, 2026-08-08; the LA-3251-MS and
   LA-3186 scans, all load-bearing values render-verified).
   Carlson–Lathrop leave :math:`\mu_1` **free** in the chapter
   (*"Selection of* :math:`\mu_1` *determines the spread of direction
   cosines along the axes"*, printed p. 32); their tabulated family
   lives in :cite:`LathropCarlson1964` and picks :math:`\mu_1` to
   *"integrate as many even powers of* :math:`\mu` *as possible"*.
   Until #337 ORPHEUS carried a **project convention**
   :math:`\mu_1^2 = 4/(N(N+2))` (source-noteless since the first
   quadrature commit, untraced to any publication) whose measured cost
   was two-fold: a degree ceiling of :math:`N-1` (at :math:`S_4` the
   single orbit forces equal weights, so :math:`\mu_1` is the only
   freedom, and the convention's :math:`0.4082483` missed the
   :math:`\int\mu^4` moment by 16.7 % where the matched
   :math:`\mu_1 = \sqrt{(5-\sqrt{10})/15} = 0.3500212` — the classic
   :math:`LQ_4` value — closes it in radicals), and a positivity
   frontier of :math:`S_{12}` where the matched family reaches
   :math:`S_{18}`.  (:math:`S_2` has no freedom at all:
   :math:`\mu_1^2 = 1/3` is forced by the diffusion condition, printed
   p. 45 — bit-identical across both seeds, the carve's regression
   control.)

   **LA-3186 Table I corroboration** (:math:`S_4/S_6/S_8/S_{12}/S_{16}/
   S_{20}` tabulated): our 50-digit re-derivation agrees at every order
   and is the *correctly rounded* value at three where the print
   carries a one-ulp last-digit slip (:math:`S_6` prints …55 for …54,
   :math:`S_{12}` …26 for …27, :math:`S_{16}` …68 for …69).
   :math:`S_{10}/S_{14}/S_{18}` are untabulated anywhere in LA-3186;
   their values rest on the re-derivation, independently reproduced
   through the C-L level system.

   **The point weights are ours, not LA-3186's — measured, deliberate.**
   LA-3186 (p. 15) decomposes level weights by the axis-weight ansatz
   :math:`p_{\{i,j,k\}} = a_i + a_j + a_k`, which constrains only the
   *axial* moment content: `[M]` on the same nodes that decomposition
   reaches full 3-D degree **11 at every order** :math:`\geq 14`, while
   ORPHEUS's per-orbit cross-moment solve reaches **15/15/17** at
   :math:`S_{14}/S_{16}/S_{18}`.  The two decompositions' positivity
   frontiers *interleave*: the ansatz dies at :math:`S_{18}` and serves
   :math:`S_{20}` (at degree 11 — strictly dominated by our
   :math:`S_{14}`), ours serves :math:`S_{18}` and dies at
   :math:`S_{20}`.  And :math:`S_{22}` is **intrinsically dead** for
   the whole even-moment family: an LP over the full decomposition
   kernel certifies no nonnegative point weights exist (best possible
   min :math:`p = -0.0027`).  LA-3251-MS's *"for n > 22 negative
   weights occur"* (printed p. 33) is a claim about LEVEL weights —
   true (first negative level weight at :math:`n = 24`) — not about
   realizable point weights; LA-3186 itself never publishes point
   weights above :math:`n = 16`.  Extraction of record with page
   cites: ``scratch/issue_337_la3186_la4058_extraction.md``.

Weights sum to :math:`4\pi`.  Carries a ``level_indices`` structure of
the same shape the cylindrical :term:`sweep` consumes — though *this*
family is refused there; see the warning below.  Unlike the product
quadrature (which has one level per :math:`\mu_z` value), the
Level-Symmetric quadrature groups both :math:`+\mu_z` and :math:`-\mu_z`
hemispheres on the same level (grouped by :math:`|\mu_z|`).

Within each level, ordinates are ordered by the **fiber's own
coordinate**: primarily by increasing :math:`\eta` (the azimuthal-sweep
convention), with ties broken by increasing :math:`\varphi` and then by
increasing :math:`\operatorname{sign}(\mu_z)`.  Naming the tie-break is
not decoration — because a level here spans both hemispheres, each
:math:`\eta` value is **4-fold degenerate** (the :math:`\pm\xi`,
:math:`\pm\mu_z` sign replications).  So :math:`\eta` alone orders a
level only down to *blocks of four*, and which member of each block
comes first was, until 2026-08-13, :func:`numpy.argsort`'s introsort
partition — an implementation detail that is not even consistent across
array sizes (numpy falls back to insertion sort at :math:`\le 16`
elements, so small levels came out stable and larger ones did not).
The full triple is injective on every level of every shipped rule, and
:meth:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure.from_level_membership`
refuses a level where it is not, so the order is a property of the rule
rather than of a sort algorithm.

.. warning::

   Having the level structure is **necessary but not sufficient** for
   cylindrical use.  Since Q5.6.3 a cylindrical ``SNMesh`` admits only
   **carrying** rules (the R12a march-start predicate,
   :ref:`sn-direct-seed-r12a`), and *every* Level-Symmetric level is
   non-carrying: the :math:`|\mu_z|` grouping puts the hemisphere pair
   on one level, so the two smallest :math:`\eta` values tie
   (``degenerate``) and the seed's :math:`(1-\tau_0)` thread weight
   vanishes.  The :math:`\sigma_y` fold does **not** repair this — the
   hemisphere pair survives the ξ-fold, so a σ_y-folded LS rule is
   still non-carrying (`[M]` Q5.6a).  A carrying LS-family rule needs
   the **double-fold** (quotient by the vertical pair as well) —
   filed-not-built as `Issue #339
   <https://github.com/deOliveira-R/ORPHEUS/issues/339>`_.  On
   cylinders today, use :meth:`Quadrature.folded_product
   <orpheus.numerics.quadrature.Quadrature.folded_product>`.

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
are the even triples up to permutation.  `[M]` the resulting degrees
(2026-08-08, on the moment-matched nodes — #337; the earlier #327
convention-seed table read degrees 3/3/5/7/9/11 with the S12 frontier):

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
     - 5
     - 0.523599
   * - 6
     - 2
     - 7
     - 0.246940
   * - 8
     - 3
     - 9
     - 0.142535
   * - 10
     - 4
     - 11
     - 0.070755
   * - 12
     - 5
     - 11
     - 0.040607
   * - 14
     - 7
     - 15
     - 0.012990
   * - 16
     - 8
     - 15
     - 0.016300
   * - 18
     - 10
     - 17
     - 0.000175

⭐ **At** :math:`S_2` **and** :math:`S_4` **the weight solve is provably a
no-op.**  Both node sets are a *single* orbit, so :math:`O_h`-invariance plus
:math:`\sum w = 4\pi` determine the one weight uniquely.  At :math:`S_2` the
NODES are forced too (:math:`\mu^2 = 1/3`), which made it bit-identical across
the #327 fix AND the #337 seed change — the family's standing regression
control.  :math:`S_4` was bit-identical across #327 (same nodes, forced
weight) but its nodes MOVED at #337 (:math:`0.4082 \to 0.3500`), taking the
workhorse order from degree 3 to degree 5 — the largest single win of the
seed adoption.

.. _quadrature-ls-positivity:

Why the family stops at :math:`S_{18}`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

`[M]` the per-orbit solve on the moment-matched nodes yields a **negative**
weight from :math:`S_{20}` upward (:math:`-2.191\times 10^{-4}` at
:math:`S_{20}`, 50-digit-confirmed — the sign has ~7 orders of margin over
float64, so the flip is the family's, not the arithmetic's).
:meth:`Quadrature.level_symmetric
<orpheus.numerics.quadrature.Quadrature.level_symmetric>` **refuses** there.
(Under the retired convention seed the frontier sat at :math:`S_{12}`:
:math:`-0.027` at :math:`S_{14}` on those nodes — the seed, not the
level-symmetric shape, sets it.)

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
limit that would drift.  Above :math:`S_{18}`, use Lebedev or the product rule
— and :math:`S_{18}` is the end of the road for THIS family, not a solver
limitation: :math:`S_{20}` is servable only by the axis-weight decomposition
at full degree 11 (dominated by our :math:`S_{14}`), and :math:`S_{22}` is
LP-certified infeasible for any nonnegative decomposition (the provenance
note above).

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

.. admonition:: Development history — #337, 2026-08-08
   :class: note

   #327's citation verification (LA-3251-MS) revealed the node seed
   :math:`\mu_1^2 = 4/(N(N+2))` was a source-noteless project convention —
   Carlson–Lathrop leave :math:`\mu_1` free and their tables moment-match
   it.  #337 adopted the moment-matched seed: :math:`\mu_1` is the smallest
   root of the :math:`\int\mu_z^N` condition, root-found by the builder
   (the residual has TWO roots at :math:`N = 6/10/14/18`; the larger one's
   weight solve is strongly negative, which is why *smallest* is the rule).
   Degrees rose :math:`N-1 \to N+1` at :math:`S_4`–:math:`S_{10}` and
   :math:`S_{14}`; the positivity frontier moved :math:`S_{12} \to
   S_{18}`; every LS consumer's nodes moved at :math:`S_4` and above
   (:math:`S_2` bit-identical), re-baselined per principled-equivalence.

   The verification instrument is three-cornered
   (``tests/numerics/test_level_symmetric_nodes.py``): the build-measured
   stamp, an independent monomial sweep, and frozen literals re-solved at
   50 digits — because with a build-measured stamp, "measured ==
   advertised" alone is two re-implementations of one gamma identity
   agreeing.  The degree table is structurally BLIND to the seed at
   :math:`S_{12}/S_{16}/S_{18}` (same degree under both seeds); those
   orders are carried by the :math:`\mu_1` value gate and the axial
   residual gate alone.

Product Quadrature (GL x equispaced)
-------------------------------------

Tensor product of Gauss--Legendre in :math:`\mu = \cos\theta` (polar)
and equispaced points in :math:`\varphi` (azimuthal).  Each :math:`\mu`
level has the same number of azimuthal points, giving a clean level
structure — the **parent** of the cylindrical production family (its
:math:`\sigma_y` quotient below; the full-circle rule itself is refused
at cylindrical ``SNMesh`` admission since Q5.6.3, because every level's
march start is an edge node or a mirror tie —
:ref:`sn-direct-seed-r12a`).  Weights:

.. math::
   :label: quadrature-product-weights

   w_{p,m} = w_{\text{GL}}(\mu_p) \cdot \frac{2\pi}{N_\varphi}

Sum to :math:`4\pi`.  Within each level, ordinates are sorted by
increasing :math:`\eta = \sin\theta\cos\varphi` to match the
:math:`\alpha` recursion convention from :cite:`BaileyMorelChang2010` (the
curvilinear :math:`\alpha`-recursion).

Built by :meth:`Quadrature.product(n_mu, n_phi)
<orpheus.numerics.quadrature.Quadrature.product>`.

.. _quadrature-folded-product:

Folded Product (the σ_y quotient — the cylindrical production family)
---------------------------------------------------------------------

The rule a cylindrical ``SNMesh`` actually admits (Q5.6.3):
the :math:`\sigma_y` **quotient** of the staggered full-circle product.
The 1-D cylindrical transport problem is invariant under
:math:`\mu_y \to -\mu_y`, so the physical angular domain is the
half-sphere orbifold; the fold keeps one representative per mirror
orbit (weights doubled off the mirror plane), turning each level's
azimuthal **circle** into an **arc** :math:`\omega \in (0, \pi)`.
Structurally (Q5.2/Q5.3): the parent is the STAGGERED product (midpoint
azimuthal nodes, offset *derived* — campaign ruling T25), the quotient
descends by **bit-copied selection** (never recomputation), and odd
:math:`n_\varphi` is refused (a staggered odd count puts one node per
level ON the mirror — :math:`|\Sigma| = n_\mu \ne 0`).  The arc's midpoint nodes
are Gauss–Chebyshev (first kind) in :math:`x = \cos\omega` — the arc's
own Gauss family, mirroring the sphere's Legendre structure.

What the fold buys is exactly the R12a carrying property on **every**
level (:ref:`sn-direct-seed-r12a`): the singular set is empty, the
march-start weights satisfy :math:`\tau \in [\tfrac14, \tfrac34]`
with the reversal identity
:math:`\tau_m + \tau_{M-1-m} = 1` to 0.5--12 ULP
(:eq:`morel-montry-folded-arc` in
:doc:`/theory/foundations/structured_geometry`; both numbers were
:math:`[\tfrac15, \tfrac45]` and *bit-exact* until the Q5.6.4 partition
fix, whose retraction note explains why the older reading was exact —
its two end cells were stretched symmetrically), and each level's
ψ½ seed is a genuine independent unknown solved exactly by route (a)'s
forward substitution.  Admission reads this **structure**, never a
provenance tag — a hand-built rule with the same arrays admits
identically
(:func:`~orpheus.sn.angular.closure.assert_carrying_quadrature`).

Built by :meth:`Quadrature.folded_product(n_mu, n_phi)
<orpheus.numerics.quadrature.Quadrature.folded_product>` (the general
quotient primitive is :meth:`Quadrature.quotient
<orpheus.numerics.quadrature.Quadrature.quotient>`).

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


