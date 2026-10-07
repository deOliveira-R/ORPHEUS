.. _theory-chart-and-chord:

=====================================================================
Charts, lines and chords — the geometric kernel of a 1-D geometry
=====================================================================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: geometry
      concept: chart, line, chord, measure density, line domain
      role: "the spatial geometry a 1-D problem keeps: the coordinate map c of a coordinate system as the quotient by its symmetry group G_c, the density of its one measure (the area of a level set), oriented lines in Plücker coordinates, the chord of a line through the level sets of c solved once in the orbit space, point location, the invariant measure on lines pushed to each chart, the directions at a point and the orbit space of oriented lines with its density"
      code: [orpheus.geometry.chart, orpheus.geometry.line, orpheus.geometry.chord, orpheus.geometry.coord]
      depends_on: [structured_geometry, manifolds]
      related: [collision_probability, boundary_conditions]


Key facts
=========

- **A 1-D coordinate system is a chart: the quotient of space by a symmetry
  group.** The coordinate a 1-D problem keeps, :math:`c(x)`, is :math:`x_0`
  on a slab, the distance to the axis on a cylinder and the distance to the
  centre on a sphere (:eq:`geometry-radial-coordinate`). Its symmetry group
  :math:`G_c = \{g \in E(3) : c \circ g = c\}` is the set of rigid motions
  the problem cannot tell apart, and :math:`c` is the quotient map onto the
  orbit space :math:`\mathbb{R}^3/G_c`. :class:`~orpheus.geometry.chart.Chart`
  derives every verb from two data, the number :math:`d` of kept frame
  columns and the linear group :math:`L \subset O(3)` of :math:`G_c`:
  the coordinate is signed where :math:`L` fixes the kept space and a norm
  where it acts on it, the singular stratum (the sphere's centre, the
  cylinder's axis) exists exactly where it acts, and a motion is in
  :math:`G_c` iff its orthogonal part is in :math:`L` and its translation
  has no kept component. The measure is :meth:`CoordSystem.measure
  <orpheus.geometry.coord.CoordSystem.measure>` (one definition)
  (:ref:`chart-and-chord-chart`).
- **The measure has one density, and it is the area of a level set.**
  :meth:`Chart.measure_density <orpheus.geometry.chart.Chart.measure_density>`
  is the derivative of the one definition, :math:`c_d\,d\,r^{d-1}`
  (:eq:`geometry-measure-density`): :math:`1`, :math:`2\pi r` and
  :math:`4\pi r^2`. Because the orbit coordinate has unit gradient, the
  coarea formula makes the same number the area of the level set
  :math:`c = r`, per unit transverse area on the slab and per unit height
  on the cylinder. A basis's mass matrix and a white wall's area both read
  it (:ref:`chart-and-chord-measure-density`).
- **A line is held in Plücker coordinates** (direction :math:`\Omega`,
  moment :math:`m = p \times \Omega`), so the base point it was built from
  is not part of it; its own parameter :math:`t` is measured from its
  foot, the point closest to the origin (:ref:`chart-and-chord-lines`).
- **The chord is solved once, in the orbit space.** Along a line the orbit
  coordinate obeys
  :math:`c(t)^2 = b^2 + \bigl(|P\Omega|(t - t^*)\bigr)^2` on the cylinder
  and the sphere (:eq:`geometry-line-crossing-law`), so every crossing is
  :math:`t^* \pm h_k/|P\Omega|` with the half-chord
  :math:`h_k = \sqrt{(r_k - b)(r_k + b)}`. Every 3-D length is an
  orbit-space length times the **obliquity** :math:`1/|P\Omega|`: one
  factor for the three charts, :math:`1`, :math:`1/\sin\theta` and
  :math:`1/|\mu|` (:eq:`geometry-cylinder-axial-factor`). Segment lengths
  are computed without cancellation (:eq:`geometry-chord-segment-lengths`)
  (:ref:`chart-and-chord-chord`).
- **No point on a chord is ever located.** A segment carries its region
  from the crossing order: an inward crossing of :math:`c = r_k` enters
  region :math:`k-1`, an outward one region :math:`k`
  (:eq:`geometry-crossing-order`). A tangency is not a
  crossing; the centre is a stratum, not a surface. A bare orbit
  coordinate is located inner-owns on the closed domain
  :math:`[r_0, r_n]`. The exteriors are the out-of-range codes :math:`n`
  (inner) and :math:`n + 1` (outer), so indexing a per-region table with
  them raises instead of reading a material, and a line lying in a
  surface is typed as on that interface (:ref:`chart-and-chord-location`).
- **The measure on lines is** :math:`\mathrm{d}A_\perp\,\mathrm{d}\Omega`,
  pushed to each chart as a density over the impact parameter
  :math:`b \ge 0` (:eq:`geometry-measure-on-lines`),
  :math:`S_{d-2}\,b^{d-2}\,|P\Omega|`: :math:`2\pi b` on the sphere,
  :math:`2|P\Omega|` per unit height on the cylinder; :math:`|\Omega_x|`
  per unit area on the slab. Cauchy's mean chord :math:`4V/S`
  (:eq:`geometry-cauchy-mean-chord`) is reproduced on all three charts to
  :math:`10^{-15}` `[M]` 2026-10-05, and the quadrature rule over
  :math:`b` stays with the consumer (:ref:`chart-and-chord-measure`).
- **A transit is a maximal run of a line inside the domain**
  :math:`[r_0, r_n]` (:eq:`geometry-transits`): the slots in which every
  traversed slot has an interior region code; an untraversed slot never
  breaks a run, so only a hollow body's traversed cavity separates two,
  and a line makes at most two. Its walls are breakpoint indices
  (:math:`0` or :math:`n`), never positions; an absent transit carries the
  out-of-range codes (slot count :math:`S`, wall :math:`n + 1`); a tangency
  is not a crossing and a parallel line has none
  (:ref:`chart-and-chord-transits`).
- **The directions at a point are** :math:`S^2/\mathrm{Stab}(x)`, a box of
  measure-uniform coordinates with the constant density :math:`4\pi/|B|`
  (:eq:`geometry-directions-at`): the cosine to the radius under
  :math:`O(2)_x` (sphere, slab; Archimedes), the axial cosine times the
  in-plane angle under :math:`D_{1h}` off the cylinder's axis, one
  coordinate or none on a stratum. The impact parameter is the kernel's
  own :meth:`Chart.image <orpheus.geometry.chart.Chart.image>`; the
  tangencies, the break set of a direction rule, are scale-free and
  include the grazing value at the point's own level; within about
  :math:`\sqrt\epsilon` of grazing there no double-precision :math:`b`
  resolves the side (:ref:`chart-and-chord-directions`).
- **The oriented lines modulo the group are a box too**
  (:eq:`geometry-line-domain`): the impact parameter :math:`b` on the
  sphere, :math:`b` and the polar angle :math:`\theta \in [0, \pi/2]` on
  the cylinder, the cosine :math:`\mu \in [-1, 1]` on the slab, with the
  densities :math:`8\pi^2 b`, :math:`8\pi\sin^2\theta` and
  :math:`2\pi|\mu|`, each the beam density times the directions the
  quotient folds. The cylinder's second coordinate is the polar angle, not
  its cosine, because in :math:`\mu_z` every integrand carries a
  square-root end at :math:`\mu_z = 1` (the user's ruling of 2026-10-06).
  Cauchy's formula ties the density to the chord with no shared formula
  (:ref:`chart-and-chord-line-domain`).
- **One consumer, a reference.** `[M]` 2026-10-07, ``git grep`` over
  ``orpheus/``: the only modules outside the kernel that import
  ``orpheus.geometry.chart``, ``.line`` or ``.chord`` are five modules of
  the characteristic reference
  (:ref:`theory-characteristic-reference`): ``walls``, ``basis``,
  ``closure``, ``transport`` and ``assembly`` under
  ``orpheus/derivations/continuous/characteristic/``. No production
  module imports the kernel: every chord, locator and measure on lines
  production computes is still its own spelling, counted in
  :ref:`chart-and-chord-deferred`.
- **Designs that do not work** (a squared crossing law for all three
  charts, the planar measure on lines, "grazing is derived on every
  chart", a directional point locator, re-solving each 3-D line's
  quadratic, a closed :math:`[0, R]` on a hollow body) are listed with the
  structural reason each fails (:ref:`chart-and-chord-refuted`).


.. _chart-and-chord-role:

Where the kernel sits
=====================

The three modules live in :mod:`orpheus.geometry`, beside
:class:`~orpheus.geometry.transformation.RigidMotion` (the Euclidean group
:math:`E(3)` acting on points and on directions), and are numpy-only: no
transport vocabulary, no cross sections, no neutrons. They answer the
geometric questions every characteristic method asks of a 1-D body (where
a straight line enters and leaves each region, how long it stays, which
region a point is in, and how lines are counted), which today are
answered separately by the Variant-α and Peierls references, the
collision-probability chords, the method of characteristics and Monte
Carlo (:ref:`chart-and-chord-deferred` has the census).

- :mod:`orpheus.geometry.chart` holds
  :class:`~orpheus.geometry.chart.Chart` (the coordinate map and its
  group), :class:`~orpheus.geometry.chart.SingularStratum`, and the two
  images of a line in the orbit space,
  :class:`~orpheus.geometry.chart.RadialImage` and
  :class:`~orpheus.geometry.chart.AxialImage`, the directions at a
  point, :class:`~orpheus.geometry.chart.DirectionDomain` with its
  :class:`~orpheus.geometry.chart.DirectionShape`, and the oriented lines
  modulo the group, :class:`~orpheus.geometry.chart.LineDomain` with its
  :class:`~orpheus.geometry.chart.LineShape`. The density of the measure is
  :meth:`CoordSystem.measure_density <orpheus.geometry.coord.CoordSystem.measure_density>`
  in :mod:`orpheus.geometry.coord`, beside the measure it differentiates.
- :mod:`orpheus.geometry.line` holds :class:`~orpheus.geometry.line.Line`.
- :mod:`orpheus.geometry.chord` holds
  :class:`~orpheus.geometry.chord.ConcentricPartition` (the level sets
  :math:`c = r_0 < \dots < r_n` posed in space),
  :class:`~orpheus.geometry.chord.Chord` (the answer for a batch of lines),
  :class:`~orpheus.geometry.chord.Crossings` and
  :class:`~orpheus.geometry.chord.Transits`.

**Relation to the shape of a 1-D problem.** A
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry`
(:ref:`structured-geometry-value`) is a coordinate system, breakpoints,
materials and boundary laws.
:meth:`ConcentricPartition.of <orpheus.geometry.chord.ConcentricPartition.of>`
reads the first two only: the materials and the laws are the consumer's.
A geometry carries no position; the partition is posed in space by a
:class:`~orpheus.geometry.transformation.RigidMotion` from its canonical
frame, the identity by default (the user's ruling of 2026-10-05). The
canonical frame is the one the angular chart already declares
(:ref:`structured-geometry-angular-chart`): the slab's normal is
:math:`\hat e_x`, the cylinder's axis :math:`\hat e_z`, the sphere's
centre the origin.

**Relation to the directional orbit spaces.** The corpus already carries
orbit spaces on the unit sphere of directions,
:math:`S^2/H` for a closed subgroup :math:`H` of :math:`O(3)`
(:ref:`manifold-orbit-space`;
:class:`~orpheus.numerics.symmetry.SubgroupOfO3`,
:class:`~orpheus.numerics.manifold.Quotient`). The chart is their spatial
twin: the same construction (a group, its orbit map, the strata where the
isotropy jumps) applied to positions under a subgroup of :math:`E(3)`. The
two are tied: the isotropy of a generic point under :math:`G_c` is the
angular symmetry the problem spends (:ref:`chart-and-chord-isotropy`).

**Relation to architecture E.** The boundary-law ontology's architecture
E (#551) names a ``Chart``: the pair (coordinate system, kept
coordinates), realizing :math:`G_c` and deriving its orbit map, its measure
and its singular strata. The user ruled on 2026-10-05 that the geometric
kernel mints that ``Chart`` now, for the three 1-D charts only, and a
second ruling of the same day that it is derived from (kept columns,
group), not matched on the coordinate system. In one dimension both data
are implied by the coordinate system, so the class's one field is the
coordinate system, and :attr:`~orpheus.geometry.chart.Chart.kept_columns`
(the exponent :math:`d` of the measure coordinate :math:`T = r^d`) and
:attr:`~orpheus.geometry.chart.Chart.group` (``O2("x")``, ``Dinfh``,
``O3``) are read from it at one site; architecture E adds the
kept-coordinates field with the :math:`(r, z)` and 2-D charts, the deck
group and the boundary laws that read them.

**The import edges.** :mod:`orpheus.geometry.chart` imports
:mod:`orpheus.geometry.coord`, :mod:`orpheus.geometry.transformation` and
``orpheus.numerics.symmetry`` (for :class:`~orpheus.numerics.symmetry.SubgroupOfO3`
and its realization); :mod:`orpheus.geometry.line` imports the
transformation module; :mod:`orpheus.geometry.chord` imports the chart,
the line and the transformation module. ``orpheus.numerics.symmetry``
itself imports ``orpheus.geometry.transformation``, the package cycle
documented in ``tests/gates/test_layer_imports.py``; the chart imports
submodules by their full names, which is the form that file admits.


.. _chart-and-chord-chart:

The chart and its group
=======================

The coordinate map
------------------

Each 1-D coordinate system names one function on space, the position the
problem keeps. In the canonical frame,

.. math::
   :label: geometry-radial-coordinate

   c(x) \;=\;
   \begin{cases}
     x_0 & \text{slab (CARTESIAN), signed,} \\[2pt]
     \sqrt{x_0^2 + x_1^2} & \text{cylinder (CYLINDRICAL), the distance to the axis } \hat e_z, \\[2pt]
     |x| = \sqrt{x_0^2 + x_1^2 + x_2^2} & \text{sphere (SPHERICAL), the distance to the centre.}
   \end{cases}

:meth:`Chart.orbit_coordinate <orpheus.geometry.chart.Chart.orbit_coordinate>`
evaluates it on a batch of points ``(..., 3)``. The slab's coordinate is
signed: a slab may sit across the origin (breakpoints
:math:`r_0 < 0 < r_n` are admitted), and the reflection
:math:`x_0 \mapsto -x_0` is not a symmetry of the problem. The name is
"orbit coordinate", not "radial coordinate", because on the slab it is
not a radius.

The symmetry group
------------------

The rigid motions under which the coordinate does not change form a group,

.. math::

   G_c \;=\; \{\, g \in E(3) : c \circ g = c \,\},

and :math:`c` is the quotient map :math:`\mathbb{R}^3 \to \mathbb{R}^3/G_c`:
two points have the same coordinate exactly when a motion of
:math:`G_c` carries one to the other. Write a motion as
:math:`g(x) = Qx + t` with :math:`Q` orthogonal.

**Two data determine the chart.** The kept columns, the first :math:`d`
of the canonical frame (:math:`d = 1, 2, 3`), and the linear group
:math:`L` of :math:`G_c`, a closed subgroup of :math:`O(3)`. Then
:math:`G_c` is :math:`L` together with the translations along the
discarded columns, and

.. math::

   (Q, t) \in G_c \iff Q \in L \ \text{and}\ t_i = 0
   \ \text{for every kept column } i,

which is what :meth:`Chart.contains <orpheus.geometry.chart.Chart.contains>`
computes: :math:`Q \in L` by the group's realization (the same
:class:`~orpheus.numerics.symmetry.SubgroupOfO3` machinery that decides
directional symmetry), the translation by its kept entries. Per chart:

.. list-table:: The three groups and their orbit spaces
   :header-rows: 1
   :widths: 14 30 26 30

   * - Chart
     - :math:`d`, :math:`L`
     - :math:`g = (Q, t) \in G_c` iff
     - Orbit space, generic orbit
   * - slab
     - 1, :math:`O(2)_x` (``O2("x")``)
     - :math:`Q\hat e_x = \hat e_x` and :math:`t_x = 0`
     - the line :math:`\mathbb{R}`; a plane
   * - cylinder
     - 2, :math:`D_{\infty h}` (``Dinfh``), the stabiliser of the axis as a line
     - :math:`Q\hat e_z = \pm\hat e_z` and :math:`t_x = t_y = 0`
     - the half-line :math:`r \ge 0`; a cylinder surface
   * - sphere
     - 3, :math:`O(3)` (``O3``)
     - :math:`t = 0`
     - the half-line :math:`r \ge 0`; a sphere

**How** :math:`L` **acts on the kept space decides the rest.** :math:`L`
either fixes the kept space pointwise (the slab: :math:`O(2)_x` holds no
reflection of :math:`\hat e_x`) or acts on it as :math:`O(d)` (the
cylinder and the sphere).
:attr:`Chart.acts_on_kept_space <orpheus.geometry.chart.Chart.acts_on_kept_space>`
decides which by asking the realization whether :math:`L` contains the
reflection of :math:`\hat e_x`. Where :math:`L` fixes the kept space the
orbit coordinate is the kept coordinate itself, signed, and the image of a
line is affine; where it acts, the coordinate is the norm of the kept
components, there is a singular stratum :math:`c = 0` with isotropy
:math:`L`, and the image of a line is a straight line at an impact
parameter (:ref:`chart-and-chord-crossing-law`).

Why the membership conditions are the group, chart by chart, as a check
on the table (the kernel derives them; it does not spell them):

Why each condition is necessary and sufficient, chart by chart:

- **Slab.** :math:`c(Qx + t) = (Qx)_0 + t_0` must equal :math:`x_0` for
  every :math:`x`. At :math:`x = 0` this gives :math:`t_0 = 0`; then
  :math:`\hat e_x^{\mathsf T} Q x = \hat e_x^{\mathsf T} x` for every
  :math:`x` gives :math:`Q^{\mathsf T}\hat e_x = \hat e_x`, which for an
  orthogonal :math:`Q` is :math:`Q\hat e_x = \hat e_x`.
- **Cylinder.** Write :math:`P = \mathrm{diag}(1, 1, 0)`, the projection
  onto the plane normal to the axis, so :math:`c(x) = |Px|`. At
  :math:`x = 0`, :math:`|Pt| = 0`, so :math:`t_x = t_y = 0`. Then
  :math:`|PQx| = |Px|` for every :math:`x`; taking :math:`x = \hat e_z`
  gives :math:`PQ\hat e_z = 0`, so :math:`Q\hat e_z = \pm\hat e_z`.
  Conversely such a :math:`Q` maps the plane normal to the axis onto
  itself isometrically, and a translation along :math:`\hat e_z` changes
  nothing.
- **Sphere.** :math:`|Qx + t| = |x|` at :math:`x = 0` gives :math:`t = 0`,
  and every orthogonal :math:`Q` preserves the norm.

The translation's kept entries are compared with zero at
``_MEMBERSHIP_ATOL`` :math:`= 10^{-12}` absolute; the orthogonal part is
the realization's decision. A motion on :math:`\mathbb{R}^d` with
:math:`d \ne 3` is refused. The decision is closed-form, not sampled:
testing :math:`c \circ g = c` on random points (as the architecture-E
attack did, on 400) decides only what it samples (`vv-principles`
anti-pattern 13), while the realization decides the group. `[M]`
2026-10-05 (the elegance review): before the derivation, a five-arm match
on the coordinate system agreed with "realization of :math:`L` and no
kept translation" on 750 of 750 motions per chart, and the slab tested
against the wrong group ``O2("z")`` disagreed on 100 of 100.

`[M]` 2026-10-05 (the archivist's probe over the kernel, seed 20261005):
``contains`` reads ``True`` on 12 of 12 members (sphere: five random
orthogonal matrices from a QR factorisation, both determinant signs;
cylinder: the rotation about :math:`\hat e_z` by :math:`\sqrt 2` rad, the
mirror :math:`y \mapsto -y`, the mirror :math:`z \mapsto -z`, the
translation by 3.7 along :math:`\hat e_z`; slab: the translation
:math:`(0, 2, -1.3)`, the rotation about :math:`\hat e_x` by :math:`\sqrt 2`
rad, the mirror :math:`y \mapsto -y`) and ``False`` on 3 of 3
non-members (the sphere's centre moved by :math:`0.3\,\hat e_x`, the
cylinder tilted by 0.4 rad about :math:`\hat e_x`, the slab rotated by
0.4 rad about :math:`\hat e_z`). A rotation by :math:`\sqrt 2` rad
generates a dense subgroup of the rotations about its axis; four right
angles would generate :math:`C_4` only.

.. _chart-and-chord-isotropy:

The generic isotropy is the angular symmetry the problem spends
----------------------------------------------------------------

The stabiliser of a point under :math:`G_c`, read for its linear part, is
the set of direction maps that leave the problem at that point unchanged.
At a generic point, in the chart frame (the polar axis
:math:`\hat e_x` pointing along the orbit coordinate's gradient):

- **Sphere**, at :math:`(r, 0, 0)` with :math:`r > 0`: :math:`t = 0` and
  :math:`Q\hat e_x = \hat e_x`, so the isotropy is :math:`O(2)_x`.
- **Slab**, at any point: :math:`Q\hat e_x = \hat e_x` with :math:`t`
  chosen to fix the point, so the isotropy is :math:`O(2)_x`.
- **Cylinder**, at :math:`(r, 0, 0)` with :math:`r > 0`:
  :math:`t_x = t_y = 0` forces :math:`Q\hat e_x = \hat e_x`, and
  :math:`Q\hat e_z = \pm\hat e_z` with a translation along the axis
  restoring :math:`z`, so the isotropy is
  :math:`\{e, \sigma_y, \sigma_z, C_2^x\} = D_{1h}` about the radial
  axis.

These are the spent and unspent groups of the angular-symmetry table
``GEOMETRY_ANGULAR_SYMMETRY`` in :mod:`orpheus.numerics.quadrature.registry`
(slab and sphere: ``O2("x")`` spent; cylinder: ``Dnh(1)`` unspent), read
off the group rather than tabulated. `[M]` 2026-10-01 (the
boundary-ontology attack, ``scratch/boundary_ontology/attack_math.md``, Q1):
4 of 4 rows of that table reproduce from :math:`G_c`. The same attack
measured that the exponents :math:`d = 1, 2, 3` of the measure coordinate
:math:`T(r) = r^d` are the coarea densities of the orbits of :math:`G_c`.
The table does not read the chart: it is keyed by a geometry string and
consumed by the quadrature registry (re-keying it on :math:`G_c` is
#551).

.. _chart-and-chord-strata:

The singular strata
-------------------

The generic orbit has dimension 2 (a plane, a cylinder surface, a sphere),
so the orbit space has dimension 1. Where the isotropy is larger than
generic the orbit drops dimension and the orbit coordinate is not a
submersion:

.. list-table::
   :header-rows: 1
   :widths: 14 22 18 46

   * - Chart
     - Stratum
     - Isotropy
     - Why it matters
   * - sphere
     - the centre, :math:`c = 0`
     - :math:`O(3)`
     - the orbit is a point; :math:`\nabla c = x/|x|` is undefined
   * - cylinder
     - the axis, :math:`c = 0`
     - :math:`D_{\infty h}`, the stabiliser of the axis as a line
     - the orbit is a line; :math:`\nabla c` is undefined
   * - slab
     - none
     - --
     - :math:`c = x_0` is a submersion everywhere

:attr:`Chart.singular_strata <orpheus.geometry.chart.Chart.singular_strata>`
returns them as values,
:class:`~orpheus.geometry.chart.SingularStratum` ``(orbit_value,
isotropy)``: one stratum :math:`c = 0` whose isotropy is the chart's whole
linear group :math:`L`, exactly when :math:`L` acts on the kept space, and
none otherwise. The consequence for the chord is the gotcha below: a level set
:math:`c = 0` is a point or a line, never a surface, so a solid body's
breakpoint :math:`r_0 = 0` is never crossed. A line through the centre
passes through the stratum; it does not cross a circle of radius zero.
The directional analogue on :math:`S^2` is the pole of
:math:`S^2/O(2)_a` (:ref:`manifold-singular-stratum`).

The measure
-----------

:meth:`Chart.measure <orpheus.geometry.chart.Chart.measure>` is
:meth:`CoordSystem.measure <orpheus.geometry.coord.CoordSystem.measure>`,
delegated and never re-derived:
:math:`m_j = c_d\,(r_{j+1}^d - r_j^d)` with
:math:`(c_d, d) = (1, 1), (\pi, 2), (\tfrac43\pi, 3)`, a slab's length per
unit transverse area, a cylinder's area per unit height, a sphere's volume
(:ref:`structured-geometry-one-measure`).

.. _chart-and-chord-measure-density:

The measure density
-------------------

The measure is defined once, by the measure of a cell. With the measure
coordinate :math:`T(r) = r^{d}`
(:class:`~orpheus.geometry.coord.MeasureCoordinate`, its exponent the
number :math:`d` of kept columns) and the constant :math:`c_d`
(``measure_constant``), :math:`m([a, b]) = c_d\,(T(b) - T(a))`. Its density
in the orbit coordinate is the derivative of that one definition, and it
is also the area of the level set at :math:`r`:

.. math::
   :label: geometry-measure-density

   \frac{\mathrm{d}m}{\mathrm{d}r} \;=\; c_d\,T'(r) \;=\; c_d\,d\,r^{d-1}
   \;=\; \bigl|\{x : c(x) = r\}\bigr|
   \;=\;
   \begin{cases}
     1 & \text{slab, per unit transverse area,} \\
     2\pi r & \text{cylinder, per unit height,} \\
     4\pi r^{2} & \text{sphere.}
   \end{cases}

.. implements:: geometry-measure-density
   :by: orpheus.geometry.coord.CoordSystem.measure_density

   **Implemented by** ``CoordSystem.measure_density``, which returns
   :math:`c_d\,T'(r)` with :math:`T'` from
   ``MeasureCoordinate.derivative``; ``Chart.measure_density`` delegates
   to it, as ``Chart.measure`` delegates to ``CoordSystem.measure``.

.. implements:: geometry-measure-density
   :by: orpheus.geometry.coord.MeasureCoordinate.derivative

.. implements:: geometry-measure-density
   :by: orpheus.geometry.chart.Chart.measure_density

**The two readings are one number.** The first equality is calculus:
:math:`m([r, r + \mathrm{d}r]) = c_d\,T'(r)\,\mathrm{d}r`. The second is
the coarea formula, :math:`\mathrm{d}V = |\nabla c|^{-1}\,\mathrm{d}A\,\mathrm{d}c`
on each level set, with :math:`|\nabla c| = 1` on all three charts away
from the singular stratum (:math:`\nabla c = \hat e_x` on the slab,
:math:`Px/|Px|` on the cylinder, :math:`x/|x|` on the sphere). So the
volume between two neighbouring level sets is the area of either times
their distance, and the density is that area: the sphere's surface
:math:`4\pi r^2`, the circumference :math:`2\pi r` of a cylinder per unit
height, and 1 for a plane per unit area. `[M]` 2026-10-07, the kernel at
:math:`r = 1.5`: 28.2743338823081 (:math:`9\pi`), 9.42477796076938
(:math:`3\pi`) and 1.

**Why the density is a kernel verb.** Two consumers read it, and before
the verb each spelled it: the characteristic reference's panel basis
derived :math:`c_d\,d\,r^{d-1}` inside itself for its mass matrix, by the
user's ruling of the second rung (2026-10-06), and the reference's white
walls need the area :math:`A_w` of each wall (the factor
:math:`D = A_w/4` of :ref:`characteristic-wall-coupling`). Two spellings of
one number agree only by construction, so the third rung's ruling of the
same day (the plan's ledger, "on P1 step (b)'s third rung", Q4) moved the
density to the kernel, as the derivative of the measure and not as a
second formula, and retired the basis's own ``PanelBasis.volume_density``
onto it. A third spelling stands in production:
``orpheus.geometry.coord.compute_areas_1d``, the face areas of a 1-D
mesh. Retiring it onto this verb moves bits on the sphere in production
S\ :sub:`N`, so it is filed separately, `#584
<https://github.com/deOliveira-R/ORPHEUS/issues/584>`_.

**The gates.** ``tests/gates/geometry/test_measure_density.py``:
``test_the_measure_density_is_each_charts_level_set_area`` compares the
density with :math:`1`, :math:`2\pi r` and :math:`4\pi r^2` typed per
chart (never from ``measure_constant``) to 4 ulp, and checks that the
chart's verb is the coordinate system's bit for bit;
``test_the_density_integrates_to_the_one_measure`` integrates the density
with two-point Gauss–Legendre (exact for its degree, at most 2) over 2000
seeded intervals and compares with the measure in mpmath at 40 digits.
That comparison needs the conditioning factor
:math:`1 + \max(|a|, |b|)/(b - a)`: the measure :math:`c_d(b^d - a^d)`
cancels on a narrow interval far from 0, so a raw ulp count is not a
statement about the density (`[M]` 2026-10-06, the test-architect's
prototype: 2568 ulp raw on the cylinder, 0.9 ulp over the factor). The
route is gated where the consumers live,
``tests/gates/derivations/test_characteristic_assembly.py::test_the_mass_and_the_wall_area_read_the_one_density``:
with ``MeasureCoordinate.derivative`` doubled in process, the basis's
mass doubles and the white wall's response halves bit for bit, so no
second spelling of the density survives on either side. The rows are
``foundation`` until their ``verifies`` markers name this label.

The invariant-theory reading
----------------------------

The W5 cross-domain review (2026-10-05) gave the chart a second reading,
kept here because it explains the slab's signed coordinate. The
polynomial generator of the invariants of :math:`G_c` is
:math:`\pi(x) = X^{\mathsf T} Q_0 X` with :math:`X = (x, 1)`:
:math:`|x|^2` (sphere) and :math:`x_0^2 + x_1^2` (cylinder), of degree 2,
and :math:`x_0` (slab), of degree 1. The orbit coordinate is
:math:`c = \sqrt{\pi}` on the two curved charts and :math:`c = \pi` on
the slab. A degree-2 invariant satisfies :math:`\pi \ge 0`, and its
boundary :math:`\pi = 0` is the singular stratum; a degree-1 invariant
has no such boundary, which is why the slab has no stratum and a signed
coordinate. Along a line, :math:`\pi` is a polynomial of degree at most
2 in the line parameter: that is the crossing law below. The concentric
partition is then the pencil of quadrics :math:`\pi = r_k^2` (or
:math:`\pi = r_k`), and placing it in space is a congruence of
:math:`Q_0`; `[M]` 2026-10-05 in that review, the congruence moved 600
random lines by at most :math:`4 \times 10^{-12}` in :math:`t` and the
wrong-side congruence :math:`H^{\mathsf T} Q H` differed on 522 of 600.
The kernel computes in the orbit space rather than on the pencil, for the
conditioning reasons of :ref:`chart-and-chord-conditioning`.


.. _chart-and-chord-lines:

Lines
=====

Plücker coordinates
-------------------

An oriented line is the set :math:`\{p + t\,\Omega : t \in \mathbb{R}\}`
of a point :math:`p` and a unit direction :math:`\Omega`. It does not
depend on which of its points is named, and its Plücker coordinates make
that independence explicit:

.. math::

   \bigl(\Omega,\; m\bigr), \qquad m = p \times \Omega, \qquad
   \Omega \cdot m = 0 .

Moving the base point along the line, :math:`p \to p + s\Omega`, leaves
the moment unchanged, because :math:`\Omega \times \Omega = 0`. The
**foot**, the point of the line closest to the origin, is
:math:`\Omega \times m` (expand
:math:`\Omega \times (p \times \Omega) = p - (p\cdot\Omega)\Omega` for a
unit :math:`\Omega`), and its distance from the origin is :math:`|m|`.
The line's own parameter :math:`t` is measured from the foot:
:math:`x(t) = \text{foot} + t\,\Omega` and
:math:`t = (x - \text{foot})\cdot\Omega`. The moment is the angular
momentum about the origin of a unit-speed particle on the line, which is
why the impact parameter of a sphere centred at the origin is
:math:`|m|`.

:class:`~orpheus.geometry.line.Line` holds ``direction`` and ``moment``,
each ``(..., 3)``, one line per leading index: every consumer (a
characteristic sweep over positions and directions, an oracle over a
thousand samples) asks about many lines at once, so a per-line object
would force a Python loop. The verbs:

- :meth:`Line.through <orpheus.geometry.line.Line.through>` builds lines
  from points and directions (broadcast);
- :attr:`~orpheus.geometry.line.Line.foot`,
  :meth:`~orpheus.geometry.line.Line.parameter_of` and
  :meth:`~orpheus.geometry.line.Line.at` move between points and the
  parameter;
- :meth:`~orpheus.geometry.line.Line.reversed` is :math:`(-\Omega, -m)`:
  the same set, the opposite orientation, the same foot, every parameter
  negated;
- :meth:`~orpheus.geometry.line.Line.moved_by` applies a rigid motion:
  the direction moves linearly
  (:meth:`RigidMotion.on_directions
  <orpheus.geometry.transformation.RigidMotion.on_directions>`), the foot
  affinely (:meth:`RigidMotion.on_points
  <orpheus.geometry.transformation.RigidMotion.on_points>`), and the image
  line passes through the image of the foot. The parameter is not carried
  over: the image's foot is again the point closest to the origin, not
  the image of the old foot, so ``moved_by(g).at(t)`` is not
  ``g.on_points(at(t))``. The two parameterisations differ by one shift
  per line, ``moved_by(g).parameter_of(g.on_points(foot))``, because a
  rigid motion preserves length along the line.

A ray is a line with a start parameter, not a separate type: a backward
characteristic from :math:`x` along :math:`\Omega` is the line through
:math:`x` with direction :math:`-\Omega` (or the same line read for
:math:`t \le t_x`), and :meth:`Chord.lengths_beyond
<orpheus.geometry.chord.Chord.lengths_beyond>` reads a chord from a start
parameter on. A direction is a point of the unit sphere, the manifold the
angular quadratures' ordinates live on
(:class:`~orpheus.numerics.manifold.Sphere`); the kernel holds it as the
array of its Cartesian components.

The direction is refused, never renormalised
--------------------------------------------

A direction whose norm departs from 1 by more than 8 units of the float
spacing at 1 (``_UNIT_ULPS``) is refused at construction, and so is a
non-finite direction, base point or moment. A direction computed from
angles, or normalised once, is unit to a few roundings; a vector that is
not unit is a caller error, and renormalising it would hide the error the
caller made (parse at the boundary). The unit test is written as "every
departure within the band", not "some departure beyond it", because a
NaN departure fails every comparison and only the first spelling refuses
it. The direction and the moment must share one shape ``(..., 3)``.

.. warning::

   :class:`~orpheus.geometry.line.Line` is declared with ``eq=False``: two
   lines compare by identity, not by value. "Two base points on one line
   give one line" means that their moments are equal up to the rounding
   of the cross product, which is a numerical statement about the
   coordinates; it is not ``==``, and a line is not usable as a
   content-identity cache key as it stands.


.. _chart-and-chord-chord:

The chord, solved once in the orbit space
=========================================

A :class:`~orpheus.geometry.chord.ConcentricPartition` is a chart, its
breakpoints :math:`r_0 < r_1 < \dots < r_n` (at least two, finite,
:math:`r_0 \ge 0` on the cylinder and the sphere; :math:`r_0 = 0` is a
solid body, :math:`r_0 > 0` a hollow one) and a pose. Region :math:`j`
lies between :math:`r_j` and :math:`r_{j+1}`. A line meets the partition
along its **chord**: the ordered crossings of the level sets
:math:`c = r_k`, and the stretches between them, each in one region.
:meth:`ConcentricPartition.chord
<orpheus.geometry.chord.ConcentricPartition.chord>` moves the line into
the canonical frame by the inverse of the pose, solves there, and reports
every parameter (the slot starts, the closest approach, the crossings) on
the caller's line, measured from the caller's foot: a rigid motion
preserves length along a line, so the two parameterisations differ by one
shift per line, the caller's parameter of the image of the canonical
foot. The chord's ``line`` field is the caller's line, and its ``image``
is the canonical image shifted onto the caller's parameters.

.. _chart-and-chord-crossing-law:

The crossing law
----------------

Let :math:`P` be the orthogonal projection onto the kept columns:
:math:`P = I` on the sphere, :math:`P = \mathrm{diag}(1, 1, 0)` on the
cylinder, :math:`P = \mathrm{diag}(1, 0, 0)` on the slab. Where
:math:`L` acts on the kept space (the cylinder, the sphere),
:math:`c(x) = |Px|`, and along the line
:math:`x(t) = f + t\,\Omega` (with :math:`f` the foot),

.. math::

   c(t)^2 = |Pf|^2 + 2t\,(Pf\cdot P\Omega) + t^2\,|P\Omega|^2 .

For :math:`|P\Omega| > 0`, completing the square gives

.. math::
   :label: geometry-line-crossing-law

   c(t)^2 \;=\; b^2 + \bigl(|P\Omega|\,(t - t^*)\bigr)^2,
   \qquad
   t^* = -\frac{Pf \cdot P\Omega}{|P\Omega|^2},
   \qquad
   b = |Pf + t^*\,P\Omega| ,

so the line crosses :math:`c = r_k` if and only if :math:`b < r_k`, at

.. math::

   t = t^* \pm \frac{h_k}{|P\Omega|},
   \qquad
   h_k = \sqrt{(r_k - b)(r_k + b)} .

:math:`b` is the **impact parameter**, the least orbit coordinate along
the line; :math:`t^*` is the parameter of closest approach; :math:`h_k` is
the half-chord of the circle (or sphere) :math:`c = r_k` in the orbit
space. Geometrically: the line's image in the orbit plane (the plane
normal to the axis, or the line's own plane through the centre) is a
straight line at distance :math:`b` from the axis or centre, traversed at
speed :math:`|P\Omega|`. On the sphere :math:`|P\Omega| = 1` and the foot
is already the closest point, so :math:`t^* = 0` and :math:`b = |m|`.

Where :math:`L` fixes the kept space (the slab) the image is affine,
:math:`c(t) = c_{\text{foot}} + \Omega_x t`, and the crossings are
:math:`t_k = (r_k - c_{\text{foot}})/\Omega_x`. There is no impact
parameter: an affine function on a line has no minimum unless it is
constant.

:meth:`Chart.image <orpheus.geometry.chart.Chart.image>` computes the
push-forward once and returns it as a value of one of two types, so the
question "which kind of chart" is asked at one site:
:class:`~orpheus.geometry.chart.RadialImage` ``(impact_parameter,
origin_position, speed, parameter_origin)`` where :math:`L` acts, and
:class:`~orpheus.geometry.chart.AxialImage` ``(foot_coordinate, rate)``
(with ``rate`` :math:`= \Omega_x`, signed) where it does not. The radial
image is held in orbit-space units: with :math:`s` the signed position
along the image from its closest point, :math:`s(t) = s_0 + |P\Omega|(t -
t_0)`, so every parameter the chord needs is :math:`t_0 + (\pm h -
s_0)/|P\Omega|`, a finite numerator over one division, and the closest
approach :math:`t^* = t_0 - s_0/|P\Omega|` is a derived property. That form
never adds two infinities at a subnormal :math:`|P\Omega|` (the
conditioning table). Each image owns ``orbit_coordinate_at``
(:math:`\sqrt{b^2 + s(t)^2}` as a ``hypot``, or :math:`c_{\text{foot}} +
\dot c\, t`), ``parallel``, ``speed``, the coordinate a parallel line keeps,
and ``shifted(s)``, the same image with every parameter increased by
:math:`s`. A :class:`~orpheus.geometry.chord.Chord` holds the
image, so it has no field that is meaningful on one kind of chart and
``None`` on the other.

.. _chart-and-chord-obliquity:

The obliquity — one factor for three charts
-------------------------------------------

The projection :math:`P` of a unit-speed point on the line moves through
the kept space :math:`\mathbb{R}^d` at speed :math:`|P\Omega|`, so a
length :math:`\ell` measured along the image line in the kept space is
the 3-D length

.. math::
   :label: geometry-cylinder-axial-factor

   \ell_{3\text{D}} \;=\; \frac{\ell_{\text{orbit}}}{|P\Omega|},
   \qquad
   |P\Omega| =
   \begin{cases}
     |\Omega_x| = |\mu| & \text{slab,} \\
     \sqrt{\Omega_x^2 + \Omega_y^2} = \sin\theta & \text{cylinder,} \\
     1 & \text{sphere,}
   \end{cases}

where :math:`\theta` is the angle between the direction and the
cylinder's axis and :math:`\mu` the cosine to the slab's normal. The
reciprocal :math:`1/|P\Omega|` is the **obliquity**, the secant of
radiative transfer: the slab's :math:`1/|\mu|`, the cylinder's axial
:math:`1/\sin\theta`, the sphere's :math:`1`.
:meth:`Chart.projected_speed <orpheus.geometry.chart.Chart.projected_speed>`
evaluates :math:`|P\Omega|`, the norm of the direction's kept components,
scaled by their largest magnitude so that no square underflows or
overflows. The stored quantity is the speed, not the obliquity, because
the speed is finite on a line parallel to the orbit space and the
obliquity is not.

.. note::

   :math:`|P\Omega|` is the speed in the kept space, before the quotient
   by :math:`L`; it is not the rate at which the orbit coordinate
   changes. On the sphere :math:`\mathrm{d}|x|/\mathrm{d}t` runs over
   :math:`[-1, 1]` along a line while :math:`|P\Omega| = 1`; the two
   agree only far from the closest approach.

Today the tree writes this factor three ways
in three places (the Variant-α axial lift, MoC's
``seg.length / sin_p``, and inside the Bickley–Naylor functions of the
cylinder's collision probabilities, which integrate it over the polar
angle); the label names the cylinder
because the cylinder is where the factor is not 1 and not the familiar
:math:`1/|\mu|`. On the cylinder the in-plane chord is solved once and
every axial cosine rescales it, which is the fibration the unmerged hoist
``a336bde4`` measured bit-identical to the oracles it replaced.

.. _chart-and-chord-slots:

Slots, and the segment lengths without cancellation
---------------------------------------------------

The chord of a batch of lines has a fixed shape, so that a batch is one
array: a sequence of **slots**, the stretches between consecutive
potential crossings in the order a line meets them, each in the region
the first of the two crossings enters. A slot the line does not traverse
has length 0.

- **Cylinder and sphere**: :math:`2(n + 1)` potential crossings (each
  surface inbound, then each outbound), so :math:`2n + 1` slots. The
  regions inbound :math:`n-1, \dots, 0`; the inner exterior (the cavity
  of a hollow body, code :math:`n`); the regions outbound
  :math:`0, \dots, n-1`. The region holding the closest approach is
  traversed in its inbound and its outbound slot, split at :math:`t^*`.
- **Slab**: :math:`n + 1` potential crossings and :math:`n` slots, the
  regions in the order the line meets them (ascending for
  :math:`\Omega_x > 0`, descending otherwise).

The slots are not stored beside the crossings; they are read from them.
:attr:`Chord.slot_region <orpheus.geometry.chord.Chord.slot_region>` is
the crossings' ``region_entered`` without its last column, and
:attr:`Chord.slot_start <orpheus.geometry.chord.Chord.slot_start>` the
crossings' ``parameter`` without its last column (:math:`-\infty` on a
parallel line). An absent crossing keeps the parameter where it would be
(the closest approach on the curved charts, where its half-chord is 0), so
the slot it opens has length 0 and a start that is still in order. The
one stored length array is ``traversed_length``, the cancellation-free
lengths below; :attr:`Chord.slot_length <orpheus.geometry.chord.Chord.slot_length>`
reads it, and replaces it on a parallel line (:ref:`chart-and-chord-parallel`).

The orbit-space length of region :math:`j` on one side of the closest
approach is a difference of half-chords, :math:`h_{j+1} - h_j`, and the
kernel never computes it as one. Multiplying by
:math:`(h_{j+1} + h_j)/(h_{j+1} + h_j)` and using
:math:`h_k^2 = r_k^2 - b^2`:

.. math::
   :label: geometry-chord-segment-lengths

   \ell_j \;=\; \frac{1}{|P\Omega|}\,
   \begin{cases}
     \dfrac{(r_{j+1} - r_j)(r_{j+1} + r_j)}{h_{j+1} + h_j}
       & b < r_j \quad \text{(the line crosses both surfaces)}, \\[10pt]
     h_{j+1} & r_j \le b < r_{j+1} \quad \text{(region of closest approach)}, \\[4pt]
     0 & b \ge r_{j+1},
   \end{cases}
   \qquad
   \ell_{\text{cavity}} = \frac{2h_0}{|P\Omega|},
   \qquad
   \ell_j^{\text{slab}} = \frac{r_{j+1} - r_j}{|\Omega_x|},

each for one side (inbound or outbound) on the curved charts; the
cavity length is nonzero only on a hollow body with :math:`b < r_0`. The
numerator depends on the breakpoints alone and the denominator is a sum
of positive terms, so no digits cancel however thin the shell or however
close the line passes to a surface. The slot starts, which are the
crossings, are :math:`t^* - h_{j+1}/|P\Omega|` inbound,
:math:`t^* - h_0/|P\Omega|` for the cavity, and
:math:`t^* + h_j/|P\Omega|` outbound (:math:`t^*` for the region of
closest approach). The lengths are stored, never recovered as
differences of starts (:ref:`chart-and-chord-conditioning`). The
half-chords themselves are not stored: a consumer that needs them
computes :math:`\sqrt{(r_k - b)(r_k + b)}` from the image's impact
parameter.

`[M]` 2026-10-05, the kernel on the spec's fixtures (sphere with
breakpoints :math:`(0, 0.3, 1.1, 2.0)`, line through :math:`(0, -0.5, 0)`
along :math:`\hat e_y`, so :math:`b = 0`): the traversed slots are regions
2, 1, 0, 0, 1, 2 with lengths 0.9, 0.8, 0.3, 0.3, 0.8, 0.9 and starts
:math:`-2.0, -1.1, -0.3, 0, 0.3, 1.1` (parameters from the foot, which is
the origin here). The same partition on a cylinder, the line through
:math:`(0, -0.5, 7)` along :math:`(0, 0.6, 0.8)`: lengths 1.5, 4/3, 0.5,
0.5, 4/3, 1.5, the sphere's in-plane lengths divided by
:math:`|P\Omega| = 0.6`.

.. _chart-and-chord-crossing-order:

The region of a segment comes from the crossing order
-----------------------------------------------------

**No point on a chord is ever located** (the user's ruling of
2026-10-05). A segment's interior lies in exactly one region, and the
region is read from the crossing that opened it. For a crossing of
:math:`c = r_k` with sense :math:`s = \pm 1` (the sign of the rate of
change of :math:`c` through the crossing),

.. math::
   :label: geometry-crossing-order

   \text{region entered}(k, s) \;=\;
   \begin{cases}
     k - 1 & s = -1,\ k \ge 1, \\
     n & s = -1,\ k = 0 \quad \text{(the inner exterior)}, \\
     k & s = +1,\ k \le n - 1, \\
     n + 1 & s = +1,\ k = n \quad \text{(the outer exterior)},
   \end{cases}

and the region of a segment is the region entered at the crossing that
opens it, never the region of a point located on it. On the cylinder and
the sphere an inward crossing (:math:`s = -1`) is one before the closest
approach; on the slab, :math:`s` is the sign of :math:`\Omega_x`.
:class:`~orpheus.geometry.chord.Crossings` carries, for each potential
crossing in traversal order, its ``parameter``, the ``breakpoint`` index
:math:`k`, the ``sense`` and whether the line makes it (``present``),
every array of the batch shape ``(..., m)``, and the partition's
``n_regions``. :attr:`Crossings.region_entered
<orpheus.geometry.chord.Crossings.region_entered>` is the equation above
as a function of ``(breakpoint, sense, n_regions)``, and it is the one
place in the tree where the ruled rule is written: the slots read their
regions from it, so the region of each segment and the order of the
segments cannot disagree.

**The exterior codes are out of range** (the user's ruling of
2026-10-05). The regions are :math:`0, \dots, n - 1`; the inner exterior
is :math:`n` and the outer exterior :math:`n + 1`
(:attr:`ConcentricPartition.inner_exterior
<orpheus.geometry.chord.ConcentricPartition.inner_exterior>`,
:attr:`~orpheus.geometry.chord.ConcentricPartition.outer_exterior`). A
consumer that indexes a per-region table with ``slot_region`` therefore
gets an ``IndexError`` on any chord that enters an exterior, and extends
its table on purpose, for example with a void cavity:
``np.append(sigma_t, [0.0, 0.0])[chord.slot_region]``. `[M]` 2026-10-05:
the hollow sphere :math:`(0.5, 1, 2)`, the line through the centre along
:math:`\hat e_z`, :math:`\Sigma_t = (1, 2)`: indexing raises; the
extended table gives the optical depth 5.0, the true value. Negative
codes, which numpy reads as "count from the end", gave 7.0 silently
(:ref:`chart-and-chord-refuted`).

Why not locate the segment's midpoint, as the reference oracles do today?
A computed crossing lies on :math:`c = r_k` only to within rounding, on
either side; a segment of zero or near-zero length has a midpoint that is
a computed crossing; and every ray algorithm whose answer depends on
which side such a point falls is ill-conditioned whatever the
convention. The crossing order has no convention to choose.

`[M]` 2026-10-05, the fixture above: crossings at
:math:`t = -2.0, -1.1, -0.3, 0.3, 1.1, 2.0` of breakpoints
3, 2, 1, 1, 2, 3 with senses :math:`-,-,-,+,+,+`, entering regions
2, 1, 0, 1, 2 and the outer exterior (code 4). No crossing at
:math:`r_0 = 0`.

.. _chart-and-chord-tangency:

Tangency is not a crossing; the centre is a stratum
---------------------------------------------------

The test is :math:`b < r_k`, strict. A line tangent to :math:`c = r_k`
(:math:`b = r_k`) touches it at one point and stays on the outer side on
both sides of the touch, because :math:`c(t) \ge b` along the line with
equality only at :math:`t^*`. It therefore makes no crossing there and
creates no zero-length segment. A solid body's :math:`r_0 = 0` is never
crossed (no :math:`b` is below 0): a line through the centre or the axis
passes through the singular stratum, and the region holding it is split
only at :math:`t^*`, as every region of closest approach is.

`[M]` 2026-10-05 on the sphere :math:`(0, 0.3, 1.1, 2.0)`: the line at
:math:`b = 1.1` exactly is two slots of region 2 (length
:math:`\sqrt{2.79}` each, the halves of :math:`2\sqrt{2^2 - 1.1^2}`) and
no slot of region 1; the line at :math:`b = 2.0` has an empty chord; on
the hollow sphere :math:`(0.4, 1.1, 2.0)` the line at :math:`b = 0.4` has
no cavity slot, and the line at :math:`b = 0.2` has a cavity of length
:math:`2\sqrt{0.16 - 0.04} = 0.6928`.

.. _chart-and-chord-parallel:

The degenerate parallel line, and a line lying in an interface
--------------------------------------------------------------

When :math:`|P\Omega| = 0` (a cylinder line parallel to the axis; a slab
line parallel to the planes) the crossing law has no quadratic term and
the line keeps one orbit coordinate for its whole length: :math:`b` on
the cylinder, :math:`c_{\text{foot}}` on the slab. The sphere has no such line. The
chord then records ``parallel``, makes no crossing, and is one of:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Where the line lies
     - The chord
   * - strictly inside a region (on the axis of a solid cylinder too: the
       axis is the stratum, not a surface)
     - one slot of that region of infinite length, starting at
       :math:`-\infty`; every other slot 0
   * - in the surface :math:`c = r_k`
     - ``interface`` :math:`= k` (on the interface between regions
       :math:`k - 1` and :math:`k`); every slot 0: the line is assigned
       to neither side
   * - outside :math:`[r_0, r_n]`
     - every slot 0 (a line in a hollow cylinder's cavity: the cavity
       slot, infinite)

:attr:`Chord.interface <orpheus.geometry.chord.Chord.interface>` is the
breakpoint index :math:`k`, or :attr:`ConcentricPartition.no_interface
<orpheus.geometry.chord.ConcentricPartition.no_interface>`
:math:`= n + 1`, out of range for the breakpoints, where the line lies in
no surface; the interface code never coincides with a region or exterior
code that could be mistaken for it.

The typed interface is the user's ruling of 2026-10-05: a line lying in
a surface is on that interface, never folded into the inner-owns
convention, because no derivation picks a side there (the second-order
argument that decides a tangency needs :math:`c` to have a strict minimum
along the line, and here :math:`c` is constant). The degeneracy is
detected by exact equality, :math:`|P\Omega| = 0.0`: a nearly parallel
line is an ordinary line with a very long, finite chord. This is a
decision, recorded 2026-10-05: a line :math:`3 \times 10^{-17}` off
parallel (what a rotated pose makes of an axial line) is not parallel,
and its slot lengths of order :math:`10^{16}` are its true lengths. A
band (":math:`|P\Omega| < \epsilon` counts as parallel") would be a
tolerance the kernel had to derive, and no consumer needs one. The same
exactness decides the interface: a line whose constant coordinate rounds
a hair off a breakpoint (an axial line through a point computed on
:math:`r = 1`; `[M]` the qa review, :math:`b - 1 = -1.1 \times 10^{-16}`)
is inside one region with an infinite slot. "Never assigned to a side"
therefore holds for lines exactly in a surface; the misclassified set has
measure zero for every integral over lines.

`[M]` 2026-10-05 on the cylinder :math:`(0, 0.3, 1.1, 2.0)` along
:math:`\hat e_z`: at :math:`c = 0.7` one infinite slot of region 1; at
:math:`c = 1.1` ``interface`` 2 and every slot 0; at :math:`c = 2.5`
every slot 0. On the slab :math:`(-0.7, 0.3, 1.1, 2.0)` along
:math:`\hat e_y`: at :math:`x = 0.3` ``interface`` 1; at :math:`x = 0`
one infinite slot of region 0.

Half-lines and the coordinate along the chord
---------------------------------------------

:meth:`Chord.lengths_beyond <orpheus.geometry.chord.Chord.lengths_beyond>`
returns, per slot, the length beyond a start parameter: a backward
characteristic's first leg, or a forward one. A slot wholly beyond the
start keeps its cancellation-free length; only the slot containing the
start is cut, as the difference of its end and the start, which is the
one subtraction in the kernel and an unavoidable one. A parallel line is
returned unchanged, being unbounded both ways, and its infinite slot
never enters an arithmetic expression, so a batch that mixes parallel and
crossing lines raises no floating-point warning (`[M]` 2026-10-05, under
warnings as errors).
`[M]` 2026-10-05: from :math:`(0, 0.7, 0)` in region 1 of the sphere
:math:`(0, 0.3, 1.1, 2.0)` along :math:`-\hat e_y`, the half-line meets
regions 1, 0, 1, 2 with lengths 0.4, 0.6 (two slots of 0.3), 0.8, 0.9.

:meth:`Chord.orbit_coordinate_at <orpheus.geometry.chord.Chord.orbit_coordinate_at>`
delegates to the chord's image, which is already on the caller's
parameters (``shifted`` by the pose's parameter shift):
:math:`\sqrt{b^2 + (|P\Omega|(t - t^*))^2}` (a ``hypot``, with no
cancellation near the closest approach) where :math:`L` acts,
:math:`c_{\text{foot}} + \Omega_x t` where it does not.

.. _chart-and-chord-conditioning:

The conditioned forms, measured
-------------------------------

Every formula in the kernel is the one with no cancellation. Each row
below compares the kernel with the textbook form at the same float
inputs; the exact value is mpmath at 50 digits on those inputs.
`[M]` 2026-10-05, the archivist's probe over the kernel, seed 20261005;
the verification spec's prototype measurements
(``scratch/characteristic_architecture/seed_verification_spec.md`` §2,
§8) agree in every order of magnitude.

.. list-table::
   :header-rows: 1
   :widths: 30 26 22 22

   * - Quantity, regime
     - Draws
     - Kernel
     - Textbook form
   * - half-chord near tangency: :math:`\sqrt{(r-b)(r+b)}` against
       :math:`\sqrt{r\cdot r - b\cdot b}`
     - 2000; :math:`r \in [0.5, 3]`, :math:`b = r(1-\delta)`,
       :math:`\delta` log-uniform in :math:`[10^{-15}, 10^{-3}]`
     - :math:`\le 1` ulp
     - up to :math:`1.9 \times 10^{14}` ulp
   * - thin-shell segment: the conditioned quotient against
       :math:`h_{k+1} - h_k`
     - 2000; :math:`r_{k+1} \in [0.8, 2]`,
       :math:`r_k = r_{k+1}(1 - \epsilon)`, :math:`\epsilon` log-uniform
       in :math:`[10^{-12}, 10^{-6}]`
     - :math:`\le 4.4 \times 10^{-16}` relative (1.97 machine epsilon)
     - up to :math:`1.4 \times 10^{-4}` relative
   * - cylinder :math:`|P\Omega|`: the scaled norm of
       :math:`(\Omega_x, \Omega_y)` against :math:`\sqrt{1 - \Omega_z^2}`
     - :math:`|P\Omega| = 10^{-4}, 10^{-6}, 10^{-8}`
     - exact
     - :math:`3 \times 10^{-9}`, :math:`4.4 \times 10^{-5}`, 1 (no
       correct digit) relative
   * - impact parameter, base point far along the line: the Plücker foot
       against :math:`\sqrt{|p|^2 - (p\cdot\Omega)^2}`
     - :math:`b = 0.7`, base point :math:`10^6` and :math:`10^9` away
     - exact
     - :math:`7 \times 10^{-6}`, 0.7 absolute
   * - chord length, the same lines: stored lengths against a difference
       of base-point parameters
     - as above
     - exact
     - :math:`1.5 \times 10^{-12}`, :math:`3.2 \times 10^{-9}` relative
   * - cylinder, a line through :math:`(0.5, 0, 0)` with
       :math:`|P\Omega| = 10^{-160}` and :math:`10^{-170}` (the squares
       underflow): the scaled norm against a plain sum of squares
     - two directions
     - :math:`b = 0` exactly; not parallel
     - :math:`b = 5.6 \times 10^{-6}` at :math:`10^{-160}`; classed
       parallel, with an infinite slot, at :math:`10^{-170}` (`[M]` the qa
       review, on the kernel before the scaled norm)
   * - cylinder, the same line with a subnormal :math:`|P\Omega|`
       (:math:`10^{-310}`, :math:`5 \times 10^{-324}`): the image in
       orbit-space units (a finite numerator :math:`\pm h - s_0` over one
       division; the closest point :math:`\text{foot} - (\text{foot}\cdot u)u`
       through the unit kept direction :math:`u`), and zero lengths kept
       zero, against :math:`t^* = -(\text{foot}\cdot P\Omega)/|P\Omega|^2`
       then :math:`\text{foot} + t^* P\Omega`, and lengths times the obliquity
     - two directions (``test_a_subnormal_projected_speed_lifts_to_infinite_lengths_never_nan``)
     - :math:`b = 0`; traversed slots :math:`\infty`, untraversed 0; no NaN
     - :math:`-\infty \cdot 0`: :math:`b` NaN at :math:`10^{-310}` (`[M]`
       the archivist, 2026-10-05); :math:`0 \cdot \infty`: an untraversed slot
       NaN at :math:`5 \times 10^{-324}` (`[M]` the qa re-review); both rows
       red under the old forms

Near tangency the problem itself is ill-conditioned
(:math:`\mathrm{d}h/\mathrm{d}b = -b/h \to \infty`), so no algorithm makes
:math:`h` accurate when :math:`b` carries rounding; the first row uses the
exact :math:`b` the draw constructed, and it measures the algorithm, not
the problem.

.. _chart-and-chord-pose:

The pose, and the invariance it makes testable
----------------------------------------------

Because the partition carries a pose, two invariances are spellable:

1. the chord does not change when the line moves by an element of
   :math:`G_c` (the partition fixed);
2. the chord does not change when the line and the partition move
   together by any rigid motion.

`[M]` 2026-10-05, 500 lines per chart (base points uniform in
:math:`[-1.5, 1.5]^3`, directions isotropic), the members and
non-members listed under :ref:`chart-and-chord-chart`, over slots longer
than :math:`10^{-3}`:

.. list-table::
   :header-rows: 1
   :widths: 16 28 28 28

   * - Chart
     - Members of :math:`G_c`: max relative slot change
     - Non-member: max relative slot change
     - Line and partition moved together (a random motion)
   * - sphere
     - :math:`2.2 \times 10^{-12}`
     - 20
     - :math:`1.5 \times 10^{-12}`
   * - cylinder
     - :math:`1.0 \times 10^{-12}`
     - 63
     - :math:`1.6 \times 10^{-12}`
   * - slab
     - 0
     - 170
     - :math:`3.5 \times 10^{-13}`

The crossing parameters are on the caller's line: `[M]` 2026-10-05, a
partition posed by a random rigid motion, 200 lines per chart, the point
:math:`\text{foot} + t_k\Omega` of every crossing the line makes (764,
908 and 800 crossings on the sphere, the cylinder and the slab) maps back
to the orbit coordinate :math:`r_k` within :math:`1.3 \times 10^{-15}`,
:math:`1.8 \times 10^{-15}` and :math:`3.5 \times 10^{-14}`.

The member columns are rounding of the moved base point propagated
through the problem's own conditioning (the error in :math:`h` is
:math:`b\,\delta b / h`, large where a slot is short); the non-member
column is the loaded control, the same instrument reading a motion
outside the group.


.. _chart-and-chord-transits:

Transits — the runs of a line inside the domain
===============================================

A characteristic that reflects at the boundary of a 1-D body does not see
the regions one by one; it sees the stretches it spends inside the domain
:math:`[r_0, r_n]` between two walls. :attr:`Chord.transits
<orpheus.geometry.chord.Chord.transits>` reads them off the chord as a
:class:`~orpheus.geometry.chord.Transits` value.

The definition
--------------

Number the slots of a chord :math:`i = 0, \dots, S - 1` in traversal order
(:ref:`chart-and-chord-slots`), with region :math:`\rho_i` (the slot's
``slot_region``) and length :math:`\ell_i` (``slot_length``), and let
:math:`k_i` be the breakpoint index of crossing :math:`i`, which opens slot
:math:`i`; crossing :math:`i + 1` closes it. The traversed slots are

.. math::
   :label: geometry-transits

   \mathcal{T} \;=\; \{\, i : \ell_i > 0 \,\}
   \quad (\text{empty on a parallel line}), \qquad
   \text{a transit is a maximal } [a, z] \text{ with } a, z \in \mathcal{T}
   \text{ and } \rho_i < n \text{ for every } i \in \mathcal{T} \cap [a, z],

.. implements:: geometry-transits
   :by: orpheus.geometry.chord.Chord.transits

   **Implemented by** ``Chord.transits``, which reads the runs from the
   chord's slot lengths, slot regions and crossing breakpoints.


with **entry wall** :math:`k_a` and **exit wall** :math:`k_{z+1}`. In
words: a transit is a maximal run of slots in which every *traversed*
slot lies in an interior region (code :math:`< n`); an untraversed slot,
of length 0, never breaks a run, so only a traversed exterior slot
separates two transits; a transit begins and ends on a traversed slot;
and a parallel line (:ref:`chart-and-chord-parallel`) has no transit,
whether it lies inside a region (one infinite slot), in a surface, or
outside the domain. The transits are ordered along the line.

The value holds five arrays, each ``(..., 2)``, one column per possible
transit in the order the line meets them: ``first_slot`` (:math:`a`),
``stop_slot`` (:math:`z + 1`, one past the last traversed slot, so the
transit is the slice ``first_slot:stop_slot``), ``entry_wall``,
``exit_wall`` and ``present``. The slice may contain untraversed slots,
of length 0: on the solid sphere at :math:`b = 0.7` the transit is slots
0 to 6, and slots 2, 3 and 4 (region 0 inbound, the cavity code, region 0
outbound) have length 0.

Why at most two
---------------

A traversed exterior slot is the only separator, and a chord has at most
one exterior slot. On the cylinder and the sphere the slots are the
regions inbound :math:`n - 1, \dots, 0`, the inner exterior (code
:math:`n`, the cavity), and the regions outbound :math:`0, \dots, n - 1`:
the outer exterior :math:`n + 1` is entered by the last crossing, which
opens no slot. On the slab the slots are the :math:`n` regions and there
is no exterior slot at all. The cavity slot is traversed exactly when the
body is hollow (:math:`r_0 > 0`) and :math:`b < r_0`; a concentric
partition has one cavity, so a line makes 0, 1 or 2 transits, and the
kernel's constant ``_MAX_TRANSITS = 2`` is that count.

Why every wall is a boundary point
----------------------------------

A run begins on the first traversed slot either at the start of the
line's stay in the domain or just after the cavity. On a curved chart the
outermost slot is traversed whenever the line meets the body
(:math:`b < r_n` gives it the length of :eq:`geometry-chord-segment-lengths`,
which is positive), and it is opened by the crossing of :math:`r_n`; the
slot after a traversed cavity is opened by the outward crossing of
:math:`r_0`. The exits mirror the entries. On the slab the first slot is
opened by the crossing of the wall the line meets first. So every wall is
:math:`0` or :math:`n`, never an interface; the gate over 1500 lines per
partition asserts it on every line it draws.

Walls are breakpoint indices, never positions
---------------------------------------------

A wall is named by the breakpoint index of its crossing, :math:`0` (the
inner wall: the cavity's surface, or the slab's :math:`r_0`) or :math:`n`
(the outer wall), and never by where it sits in the value. A position
does not name a wall: on a solid body the entry and the exit wall of the
one transit are the same wall :math:`n`; on a hollow body crossed through
its cavity, the first transit's exit and the second transit's entry are
the same wall :math:`0`; on the slab the entry is :math:`0` rising and
:math:`n` falling. The index names the same wall in every case, and it
indexes a per-boundary-point table (an albedo per wall, a boundary law per
breakpoint) directly. Reading the wall from the slot position instead
(``entry_wall = first_slot``) gives the sphere's wall 0 for 3; it is the
battery's arm T2 and the hand-counted row reds on it.

The absent codes
----------------

An absent transit, the second column of a one-transit line or both
columns of a line that misses the body, is the empty slot range at the
end of the chord, ``first_slot == stop_slot ==`` :math:`S` (the slot
count), with both walls :math:`n + 1`
(:attr:`ConcentricPartition.no_interface
<orpheus.geometry.chord.ConcentricPartition.no_interface>`). Both codes
are out of range, for the same reason as the exterior codes
(:ref:`chart-and-chord-crossing-order`): indexing a slot array with
:math:`S`, or a per-breakpoint table of length :math:`n + 1` with
:math:`n + 1`, raises ``IndexError`` instead of reading a real entry, and
the slice ``first_slot:stop_slot`` of an absent transit is empty. The slot
code is :math:`S`, not :math:`n + 1`: on a curved chart :math:`n + 1` is a
real slot (slot 4 of the 7 on the solid sphere below is region 0
outbound), and a consumer would read its length silently.

Worked cases
------------

`[M]` 2026-10-06, the kernel at the branch head, each line through the
point :math:`(b, 0, 0)` along :math:`\hat e_y` (the closest approach, so
the impact parameter is :math:`b`); a transit is written
``(first_slot, stop_slot, entry_wall, exit_wall)``:

.. list-table::
   :header-rows: 1
   :widths: 34 26 40

   * - Partition and line
     - Transits
     - Why
   * - solid sphere :math:`(0, 0.3, 1.1, 2.0)`, :math:`b = 0.7`
       (:math:`n = 3`, :math:`S = 7`)
     - one, :math:`(0, 7, 3, 3)`; the second column
       :math:`(7, 7, 4, 4)`
     - a solid body has no traversed cavity: one transit, outer wall to
       outer wall, whatever regions it misses (slots 2 to 4 are 0)
   * - the same sphere, :math:`b = 0` and :math:`b = 1.1`
     - one, :math:`(0, 7, 3, 3)`
     - the centre is a stratum, never crossed; at the interior tangency
       :math:`b = 1.1` only slots 0 and 6 are traversed (1.6703 each)
   * - the same sphere, :math:`b = 2.0`
     - none
     - a tangency is not a crossing; every slot is 0
   * - hollow sphere :math:`(0.4, 1.1, 2.0)`, :math:`b = 0.2`
       (:math:`n = 2`, :math:`S = 5`)
     - two, :math:`(0, 2, 2, 0)` and :math:`(3, 5, 0, 2)`
     - the cavity slot 2 is traversed (length 0.6928); the first transit
       ends on the inner wall, the second begins on it
   * - the same hollow sphere, :math:`b = 0.4` exactly
     - one, :math:`(0, 5, 2, 2)`
     - :math:`b = r_0` is a tangency: the cavity slot has length 0 and
       does not break the run
   * - the hollow cylinder, :math:`b = 0.2`, direction
       :math:`(0, 0.6, 0.8)`
     - two, as the sphere
     - the in-plane impact parameter decides; the axial tilt only scales
       the lengths by :math:`1/|P\Omega|`
   * - slab :math:`(-0.7, 0.3, 1.1, 2.0)`, :math:`\Omega_x = 0.6` and
       :math:`-0.6`
     - one, :math:`(0, 3, 0, 3)` rising, :math:`(0, 3, 3, 0)` falling
     - the slab's slots are its regions in the order met; the walls swap
       with the orientation
   * - cylinder :math:`(0, 0.3, 1.1, 2.0)` along :math:`\hat e_z` at
       :math:`c = 0.7`, at :math:`c = 1.1`; slab along :math:`\hat e_y` at
       :math:`x = 0.5`
     - none
     - a parallel line meets no wall: inside a region it is one infinite
       slot, in a surface every slot is 0

Why the kernel carries transits
-------------------------------

The P1 design of the characteristic references closes a reflecting
boundary with the boundary resolvent :math:`P = P_0 + E\,(I - T)^{-1} X` on the boundary
trace space (the user's ruling of 2026-10-06, P1 of the plan
``.claude/plans/characteristic_reference_architecture.md``, "P1 API
sketch", items 1 and 3). The backward characteristic from
:math:`(x, \Omega)` runs to a wall, reflects, and continues; on a
specular wall of the chart's group the reflected line has the same impact
parameter, so the unfolded path is periodic, a cycle of one or two
traversals, each a transit read forward or reversed, and the rank of the
closure per line is the length of that period, derived from the walls'
partners after the boundary laws' deck maps identify them
(:ref:`characteristic-period`). Each leg of the period carries the albedo of the wall at
which the *backward* path reflects, and that pairing is made by the wall's
breakpoint index; `[M]` 2026-10-06 (the elegance review's probe, recorded
in the plan), the other pairing moved :math:`\psi` by 0.34. The transits
are pure geometry, no albedo, no cross section, so they live in the
kernel, where the reference and production can both read them (the
ruling "Transits: in the kernel"). A line lying in a surface has no
transit, and the chord's ``interface`` field names the surface; what to
do with such a line is the consumer's decision (for the characteristic
references the plan's ruling of 2026-10-06 is to refuse it).


.. _chart-and-chord-location:

Point location — three questions, three answers
===============================================

"Which region owns :math:`c = r_k`?" is three questions, told apart by
what is being located (the boundary-point question of the plan, ruled by
the user on 2026-10-05).

**A segment of a line.** Answered by the crossing order
(:ref:`chart-and-chord-crossing-order`): no point is located and no
convention is involved.

**A point with a direction** (a particle on a surface after a
distance-to-boundary step). The region it is in is the region its
direction points into. The kernel answers it with the crossing that put
the particle there: the step that hits a surface produces a crossing
(the breakpoint :math:`k` and the sense), and the region entered is read
from it. A function of (point, direction) cannot do this, because a
computed point is never exactly on :math:`c = r_k`, so it cannot tell
"on the surface" from "one ulp off it" (refuted as a spelling:
:ref:`chart-and-chord-refuted`). On a tangency there is no crossing and
the particle stays on the outer side; on a line lying in a surface the
answer is the typed interface.

**A bare point** (evaluating a field, binning a tally, the region of a
quadrature node). For a field continuous across the interface every
owner gives the same value; for a discontinuous one (cross sections, an
emission density) the value at :math:`r_k` is genuinely two-valued. The
preferred answer is that the producer which places a point records its
region (the multi-region trajectory-resolvent references return
``region_at_node`` beside their nodes), so a consumer never re-locates
it. Where a bare
locator is still needed,
:meth:`ConcentricPartition.region_containing
<orpheus.geometry.chord.ConcentricPartition.region_containing>` is
**inner-owns on the closed domain** :math:`[r_0, r_n]`. A non-finite
coordinate (NaN or :math:`\pm\infty`) is refused rather than located, and
so is a negative one on a chart whose coordinate is a distance (the
cylinder, the sphere), which would otherwise land in a cavity that a
solid body does not have:

.. math::

   \text{region } 0 = [r_0, r_1], \qquad
   \text{region } j = (r_j, r_{j+1}] \quad (j \ge 1),

with the two exteriors out of range: the inner exterior :math:`n` below
:math:`r_0` (the cavity of a hollow cylinder or sphere, the half-space
:math:`x < r_0` of a slab) and the outer exterior :math:`n + 1` above
:math:`r_n`, the same codes the crossing order uses. Inner-owns puts the outer surface :math:`c = r_n` in the
last region, where the boundary condition acts, and partitions the
closed domain exactly; outer-owns would leave the outer surface outside
every region. The domain is :math:`[r_0, r_n]`, not :math:`[0, R]`: on a
hollow body the cavity belongs to no material
(:ref:`chart-and-chord-refuted`).
:meth:`ConcentricPartition.region_at
<orpheus.geometry.chord.ConcentricPartition.region_at>` composes it with
the pose's inverse and the orbit coordinate, for points in space.

`[M]` 2026-10-05, the 1-ulp brackets (``nextafter``) pin the comparison
operator as well as its side. Sphere :math:`(0, 0.3, 1.1, 2.0)`: at
:math:`0, 0.3, 0.3^+, 0.3^-, 1.1, 1.1^+, 2.0, 2.0^+` the regions are
0, 0, 1, 0, 1, 2, 2 and 4 (the outer exterior); at :math:`-0.3`
the call is refused. Hollow sphere :math:`(0.4, 1.1, 2.0)`: at
:math:`0.4` region 0, at :math:`0.4^-` and at 0 the inner exterior (code
2). Slab :math:`(-0.7, 0.3, 1.1, 2.0)`: at :math:`-0.7` region 0, at
:math:`-0.7^-` the inner exterior (code 3), at 0.3 region 0, at 2.0
region 2.

The tree spells point location eight times with two boundary
conventions (`[M]` 2026-10-05, the geometry census at ``a336bde4``, probe
F3: at :math:`r = r_k` exactly, the Peierls ``which_annulus`` answers
:math:`k + 1` and five others :math:`k`); none of them calls this
locator.


.. _chart-and-chord-measure:

The measure on lines
====================

The invariant density
---------------------

The space of oriented lines in :math:`\mathbb{R}^3` carries a density
invariant under every rigid motion, unique up to a constant (the
kinematic density of integral geometry): for each direction, the area
element of the plane normal to it,

.. math::
   :label: geometry-measure-on-lines

   \mathrm{d}L = \mathrm{d}A_\perp\,\mathrm{d}\Omega,
   \qquad
   \int \ell_j(L)\, \mathrm{d}A_\perp = m_j \quad \text{for every } \Omega,

where :math:`\ell_j(L)` is the length the line :math:`L` spends in region
:math:`j` and :math:`m_j` the region's measure: slicing a region into
parallel lines of one direction and summing their lengths gives its
volume. It is the disintegration of the phase-space measure
:math:`\mathrm{d}x\,\mathrm{d}\Omega` along the streaming flow
:math:`x \mapsto x + s\Omega` (the cross-domain review's frame B).

Pushed to each chart, the density over the impact parameter
:math:`b \ge 0` is :meth:`Chart.beam_density
<orpheus.geometry.chart.Chart.beam_density>`. Per unit measure of the
discarded columns, the lines of one direction meet the kept space in a
beam of width :math:`|P\Omega|` per unit of the kept space's cross-section
normal to :math:`P\Omega`. Where :math:`L` acts as :math:`O(d)`, that
cross-section at impact parameter :math:`b` is a :math:`(d-2)`-sphere of
radius :math:`b`, of area :math:`S_{d-2}\,b^{d-2}` (:math:`S_0 = 2`, the
two lines at :math:`\pm b`; :math:`S_1 = 2\pi`), so

.. math::

   w(b, \Omega) \;=\; S_{d-2}\, b^{d-2}\, |P\Omega|
   \qquad (b \ge 0),

and where :math:`L` fixes the kept space the density per unit area of a
plane is :math:`|\Omega_x|`:

.. list-table::
   :header-rows: 1
   :widths: 14 26 60

   * - Chart
     - Density for one direction
     - Why
   * - sphere
     - :math:`2\pi b\,\mathrm{d}b`
     - the lines of one direction through the annulus of radii
       :math:`b, b + \mathrm{d}b` about the centre, in the normal plane
   * - cylinder (per unit height)
     - :math:`2|P\Omega|\,\mathrm{d}b`
     - the two in-plane lines at :math:`\pm b`; a beam of in-plane width
       :math:`\mathrm{d}b` is :math:`|P\Omega|\,\mathrm{d}b` wide normal
       to :math:`\Omega`, per unit height along the axis
   * - slab (per unit area of a plane)
     - :math:`|\Omega_x|`, independent of :math:`b`
     - the lines of one direction crossing unit area of a plane
       :math:`x = \text{const}`; the slab's lines of one direction have
       no impact parameter

On the cylinder the density's :math:`|P\Omega|` and the chord's
obliquity :math:`1/|P\Omega|` cancel, so
:math:`\int_0^R \ell_j\,2|P\Omega|\,\mathrm{d}b = \int_0^R 2\ell_j^{\text{orbit}}\,\mathrm{d}b`
is the same annulus area for every direction. `[M]` 2026-10-05: the
per-region areas of the cylinder :math:`(0, 0.3, 1.1, 2.0)` along four
directions with :math:`|P\Omega| = 0.6, 0.28, 1, 10^{-6}` agree with
:meth:`CoordSystem.measure <orpheus.geometry.coord.CoordSystem.measure>`
to the same printed error for all four directions.

A second reading, measured by the elegance review of 2026-10-05: the
density of a chart's lines is the measure of the coordinate system one
rank down, :math:`2\pi b\,\mathrm{d}b` being the cylindrical measure
over :math:`b` and :math:`2\,\mathrm{d}b` the Cartesian one folded onto
:math:`b \ge 0` (a sphere's
volume integrated with the cylindrical measure over :math:`b` came out
within :math:`9.2 \times 10^{-8}` at 20 000 midpoint cells; with the
Cartesian measure, :math:`-0.78`).

`[M]` 2026-10-05, the per-region identity on the kernel, the impact
parameter integrated by
:func:`~orpheus.derivations.common.quadrature_recipes.chord_quadrature`
with 16 nodes per panel: sphere :math:`(0, 0.3, 1.1, 2.0)`, every region
within :math:`3.7 \times 10^{-16}` of :math:`\tfrac43\pi(r_{j+1}^3 - r_j^3)`;
hollow sphere :math:`(0.4, 1.1, 2.0)`, within :math:`1.7 \times 10^{-16}`;
the same sphere with plain Gauss–Legendre panels at the radii (no
endpoint substitution), :math:`2.9 \times 10^{-5}`. The cylinder's
annulus areas with the same recipe: :math:`5.9 \times 10^{-7}` at 16
nodes per panel, :math:`1.4 \times 10^{-12}` at 32,
:math:`2.5 \times 10^{-16}` at 64. The recipe's substitution suits the
sphere's integrand :math:`b\,\ell(b)` better than the cylinder's
:math:`\ell(b)`; the verification spec's own probe
(``seed_spec_probes/p2_measure.py``, re-run 2026-10-05) reaches
:math:`2.7 \times 10^{-16}` at 16 nodes on the in-plane annuli with the
substitution :math:`b = r\sin\theta`, where plain Gauss–Legendre gives
:math:`4.1 \times 10^{-5}`.

.. _chart-and-chord-cauchy:

Cauchy's mean chord
-------------------

Averaging the chord length over every line that meets a convex body,
weighted by :math:`\mathrm{d}L`:

.. math::
   :label: geometry-cauchy-mean-chord

   \bar\ell
   \;=\; \frac{\int \ell\,\mathrm{d}A_\perp\,\mathrm{d}\Omega}
              {\int \mathbf{1}_{\ell > 0}\,\mathrm{d}A_\perp\,\mathrm{d}\Omega}
   \;=\; \frac{4\pi V}{\int A_{\text{proj}}(\Omega)\,\mathrm{d}\Omega}
   \;=\; \frac{4\pi V}{4\pi \cdot S/4}
   \;=\; \frac{4V}{S},

the numerator by :eq:`geometry-measure-on-lines` for each direction, the
denominator by Cauchy's projection formula (the mean projected area of a
convex body is a quarter of its surface area). Per chart, with the
densities above:

- **sphere** of radius :math:`R`:
  :math:`\int_0^R \ell\,2\pi b\,\mathrm{d}b \,/ \int_0^R 2\pi b\,\mathrm{d}b
  = \tfrac43\pi R^3/(\pi R^2) = \tfrac43 R`;
- **infinite cylinder**, per unit height:
  :math:`\int\!\mathrm{d}\Omega\int_0^R 2\ell^{\text{orbit}}\,\mathrm{d}b
  = 4\pi\cdot\pi R^2` and
  :math:`\int\!\mathrm{d}\Omega\int_0^R 2|P\Omega|\,\mathrm{d}b
  = 2R\cdot 2\pi\!\int_0^\pi \sin^2\theta\,\mathrm{d}\theta = 2\pi^2 R`,
  so :math:`\bar\ell = 2R = 4(\pi R^2)/(2\pi R)`;
- **slab** of width :math:`L`, per unit area:
  :math:`\int\!\mathrm{d}\Omega\,|\mu|\,(L/|\mu|) = 4\pi L` over
  :math:`\int |\mu|\,\mathrm{d}\Omega = 2\pi`, so :math:`\bar\ell = 2L`.

`[M]` 2026-10-05 on the kernel, :math:`R = L = 2`, the impact parameter
by ``chord_quadrature`` at 32 nodes, the polar angle by 64-point
Gauss–Legendre (in :math:`\theta \in [0, \pi]` for the cylinder, in
:math:`\mu \in (0, 1]` for the slab):

.. list-table::
   :header-rows: 1
   :widths: 18 22 20 40

   * - Body
     - :math:`\bar\ell / R`, kernel
     - Exact
     - Control (the wrong density)
   * - sphere
     - 1.333333333333334
     - :math:`4/3`
     - 1.570796 with the planar :math:`\mathrm{d}b`: :math:`\pi/2`, the
       disc's value
   * - cylinder
     - 2.000000000000000
     - 2
     - 2.467401 without the beam factor :math:`|P\Omega|`:
       :math:`\pi^2/4`
   * - slab
     - 2.000000000000001
     - 2
     - --

The sphere's and the disc's mean chords are the per-region identity
divided by the projected area, so they are not independent of it (they
add only the count of lines that hit); the cylinder's adds the direction
average, which the per-direction identity cannot see. Of the three
classical identities of integral geometry for a convex body (the integral
of the chord is :math:`4\pi V`, the measure of the lines that hit is
:math:`\pi S`, their ratio), the surface one is not computed separately
here; it is what would check the normalisation of :math:`\mathrm{d}L` on
its own (`[R]` the cross-domain review). The theorem is stated in the
corpus beside the collision-probability surface source
(:doc:`/theory/methods/collision_probability`); the kernel is the first
place it is computed.

Why the quadrature rule is the consumer's
-----------------------------------------

The measure is a fact about lines; the rule that integrates against it is
a fact about the integrand. A chord-length integrand has a
square-root endpoint :math:`\sqrt{r_k - b}` at every radius, an optical
depth multiplies it by an exponential or a Bickley function, and the rule
that absorbs both depends on what the consumer integrates.
:func:`~orpheus.derivations.common.quadrature_recipes.chord_quadrature`
lives in :mod:`orpheus.derivations`, which imports
:mod:`orpheus.geometry`; a rule inside the geometry package would invert
that edge (the layering: :ref:`architecture-layering`). So
:meth:`Chart.beam_density <orpheus.geometry.chart.Chart.beam_density>`
gives the density and stops there.


.. _chart-and-chord-directions:

The directions at a point
=========================

A reading of the angular flux at a point :math:`x` integrates over the
sphere of directions, and the problem's symmetry makes most of that
sphere redundant: the stabiliser of :math:`x` in :math:`G_c` maps the
problem at :math:`x` to itself, so the angular flux at :math:`x` is
invariant under its linear part, and so is every integrand built from the
chord of the line through :math:`x` (a rigid motion of :math:`G_c` fixing
:math:`x` carries that line to a line through :math:`x` with the same
chord). The integral over :math:`S^2` is then an integral over the orbit
space :math:`S^2/\mathrm{Stab}(x)`. :meth:`Chart.directions_at
<orpheus.geometry.chart.Chart.directions_at>` returns that orbit space as
a :class:`~orpheus.geometry.chart.DirectionDomain`, for the representative
point :math:`x = c\,\hat e_x` of the canonical frame (every point of the
orbit :math:`c` is carried there by an element of :math:`G_c`).

The definition
--------------

The domain is a box :math:`B` in at most two coordinates :math:`q`, each
on a closed interval, with a map :math:`\omega(q)` to a representative
unit direction, such that every orbit of :math:`\mathrm{Stab}(x)` on
:math:`S^2` meets :math:`\omega(B)` once (the box's boundary aside, a set
of measure zero) and the push-forward of :math:`\mathrm{d}\Omega` to the
box is a constant multiple of Lebesgue measure:

.. math::
   :label: geometry-directions-at

   \int_{S^2} f(\Omega)\,\mathrm{d}\Omega
   \;=\; \rho \int_{B} f\bigl(\omega(q)\bigr)\,\mathrm{d}q,
   \qquad
   \rho \;=\; \frac{4\pi}{|B|},
   \qquad
   \text{for every } f \text{ with } f \circ Q = f
   \ \text{for all } Q \in \mathrm{Stab}(x),

.. implements:: geometry-directions-at
   :by: orpheus.geometry.chart.Chart.directions_at

   **Implemented by** ``Chart.directions_at``, which returns the
   :class:`~orpheus.geometry.chart.DirectionDomain` of the point; its
   ``direction`` is the representative map :math:`\omega` and its
   ``density`` is :math:`\rho`.


where :math:`|B|` is the box's coordinate measure (the product of its
widths; 1 for the empty box) and :math:`\rho` is
:attr:`DirectionDomain.density
<orpheus.geometry.chart.DirectionDomain.density>`. Coordinates with a
constant density are **measure-uniform**: a consumer puts any product rule
on the box with the constant weight :math:`\rho` and no Jacobian, and
splits it where the integrand has a kink (the tangencies below). The
density is the size of a generic orbit times the local Jacobian, not the
Jacobian alone: on the cylinder :math:`\mathrm{d}\Omega =
\mathrm{d}w\,\mathrm{d}\alpha` locally, and :math:`\rho = 4` because four
directions share each point of the box.

The four shapes
---------------

:class:`~orpheus.geometry.chart.DirectionShape` names them; the shape is
decided from the chart's pair (kept columns, group) and whether the point
is on the singular stratum (:ref:`chart-and-chord-strata`):

.. list-table::
   :header-rows: 1
   :widths: 17 15 22 26 8 12

   * - Chart and point
     - :math:`\mathrm{Stab}(x)`
     - Coordinates, box
     - Representative :math:`\omega(q)`
     - :math:`\rho`
     - :math:`b(q)`
   * - sphere, centre (``WHOLE``)
     - :math:`O(3)`
     - none; the empty box
     - :math:`\hat e_x`
     - :math:`4\pi`
     - 0
   * - sphere, :math:`c > 0`; slab, every point (``COSINE``)
     - :math:`O(2)_x`
     - ``cosine`` :math:`\mu = \Omega_x \in [-1, 1]`
     - :math:`(\mu, \sqrt{(1-\mu)(1+\mu)}, 0)`
     - :math:`2\pi`
     - :math:`c\sqrt{(1-\mu)(1+\mu)}`; none on the slab
   * - cylinder, axis (``AXIAL_COSINE``)
     - :math:`D_{\infty h}`
     - ``axial_cosine`` :math:`w = |\Omega_z| \in [0, 1]`
     - :math:`(\sqrt{(1-w)(1+w)}, 0, w)`
     - :math:`4\pi`
     - 0
   * - cylinder, :math:`c > 0` (``ANGLE_AXIAL``)
     - :math:`D_{1h}`
     - ``angle`` :math:`\alpha \in [0, \pi]`, ``axial_cosine``
       :math:`w \in [0, 1]`
     - :math:`(s\cos\alpha, s\sin\alpha, w)`,
       :math:`s = \sqrt{(1-w)(1+w)}`
     - 4
     - :math:`c\sin\alpha` (:math:`w < 1`); :math:`c` at :math:`w = 1`

`[M]` 2026-10-06, the kernel at the branch head: the four shapes, axes,
bounds and stabilisers above, with densities :math:`4\pi`,
:math:`2\pi`, :math:`4\pi` and 4 printed by ``density``.

**The sphere off its centre, and the slab: Archimedes.** The stabiliser
is :math:`O(2)_x` (:ref:`chart-and-chord-isotropy`); its orbits on the
sphere of directions are the circles :math:`\Omega_x = \mu`, and the
complete invariant is :math:`\mu`. In polar coordinates about
:math:`\hat e_x`, :math:`\mathrm{d}\Omega = \mathrm{d}\mu\,\mathrm{d}\varphi`,
so the zone between :math:`\mu_1` and :math:`\mu_2` has area
:math:`2\pi(\mu_2 - \mu_1)`, independent of where the zone sits:
Archimedes' hat-box theorem (the radial projection of a zone onto the
circumscribed cylinder preserves its area). The push-forward of
:math:`\mathrm{d}\Omega` to :math:`\mu` is therefore :math:`2\pi\,\mathrm{d}\mu`,
the cosine is measure-uniform, and :math:`\rho = 4\pi/2`. The
representative is the direction of the orbit in the half-plane
:math:`\Omega_y \ge 0,\ \Omega_z = 0`. The slab and the sphere off its
centre share this shape because they share the stabiliser; what differs
is the impact parameter, which the slab does not have.

**The sphere's centre.** The stabiliser is all of :math:`O(3)`, which is
transitive on directions: one orbit, no coordinate, the empty box of
measure 1, and :math:`\rho = 4\pi`, so an invariant integrand, a constant,
integrates to :math:`4\pi f(\hat e_x)`.

**The cylinder's axis.** The stabiliser is the whole linear group
:math:`D_{\infty h}` (the axis is the stratum): the rotations about
:math:`\hat e_z`, the mirrors containing it, and :math:`\sigma_z`. Its
orbits are the pairs of circles :math:`\Omega_z = \pm w`, the invariant
is :math:`w = |\Omega_z|`, and in polar coordinates about
:math:`\hat e_z`, :math:`\mathrm{d}\Omega = \mathrm{d}\Omega_z\,\mathrm{d}\varphi`
gives :math:`2 \cdot 2\pi\,\mathrm{d}w`: :math:`\rho = 4\pi` on
:math:`[0, 1]`.

**The cylinder off its axis: D**\ :sub:`1h`. At :math:`x = (c, 0, 0)` with
:math:`c > 0`, a motion :math:`(Q, t)` of :math:`G_c` has
:math:`t_x = t_y = 0` and :math:`Q\hat e_z = \pm\hat e_z`, so :math:`Q`
maps the plane normal to the axis onto itself; fixing :math:`x` then
forces :math:`t_z = 0` and :math:`Q\hat e_x = \hat e_x`, which leaves
:math:`Q\hat e_y = \pm\hat e_y` and :math:`Q\hat e_z = \pm\hat e_z`. The
stabiliser is the four diagonal sign matrices with first entry 1,

.. math::

   D_{1h} \;=\; \{\, e,\ \sigma_y,\ \sigma_z,\ C_2(x) \,\}
   \;=\; \{\mathrm{diag}(1, \pm 1, \pm 1)\},

``SubgroupOfO3.Dnh(1)`` in the subgroup lattice. `[M]` 2026-10-06, the
realization of ``Dnh(1)`` contains :math:`e`, :math:`\sigma_y`,
:math:`\sigma_z` and :math:`C_2(x)` and refuses :math:`\sigma_x`,
:math:`C_2(z)`, :math:`C_2(y)`, the inversion and the rotation about
:math:`\hat e_x` by :math:`\sqrt 2` rad; the same nine matrices tested as
":math:`Q \in G_c` and :math:`Qx = x`" give the same nine answers. Its
orbits on directions are :math:`(\Omega_x, \pm\Omega_y, \pm\Omega_z)`,
with invariants :math:`\Omega_x`, :math:`|\Omega_y|`, :math:`|\Omega_z|`,
and the fundamental domain is the quarter sphere
:math:`\Omega_y \ge 0,\ \Omega_z \ge 0`. In polar coordinates about
:math:`\hat e_z`, with :math:`w = \Omega_z` and the azimuth :math:`\alpha`
the in-plane angle between :math:`P\Omega` and :math:`\hat e_x`,

.. math::

   \Omega = \bigl(s\cos\alpha,\ s\sin\alpha,\ w\bigr),
   \quad s = \sqrt{(1-w)(1+w)},
   \qquad
   \mathrm{d}\Omega = \mathrm{d}w\,\mathrm{d}\alpha,

and the quarter sphere is the box :math:`\alpha \in [0, \pi]`,
:math:`w \in [0, 1]`, of measure :math:`\pi`; the four copies give
:math:`\rho = 4\pi/\pi = 4`. The cylinder's generic stabiliser is smaller
than the sphere's :math:`O(2)_x` (of the rotations about :math:`\hat e_x`
it holds only the half turn), which is why its box has two coordinates
where the sphere's has one.

The impact parameter, the kernel's own
--------------------------------------

:meth:`DirectionDomain.impact_parameter
<orpheus.geometry.chart.DirectionDomain.impact_parameter>` builds the
lines through the point in the representative directions and returns
:meth:`Chart.image <orpheus.geometry.chart.Chart.image>`'s impact
parameter: one definition of :math:`b` for the reading and the chord
(instrument doctrine X4), never a second closed form. The closed forms in
the table follow from :eq:`geometry-line-crossing-law`: on the sphere
:math:`b = |x \times \Omega| = c\sqrt{1 - \mu^2}`; on the cylinder the
image in the plane normal to the axis passes through :math:`c\,\hat e_x`
along :math:`(\cos\alpha, \sin\alpha)`, at distance
:math:`c\sin\alpha` from the axis, whatever :math:`w`. At :math:`w = 1`
the line is parallel to the axis and keeps :math:`c` for its whole length,
so its impact parameter is :math:`c`, not :math:`c\sin\alpha`. A slab line
has no impact parameter (its image is affine), and the call is refused
with that message.

How exactly the reading and a chord agree, `[M]` 2026-10-06 (2000 seeded
directions at each of the sphere at :math:`c = 1.5` and :math:`0.37` and
the cylinder at :math:`c = 1.5` and :math:`2.0`):

- against ``Chart.image`` of the same lines, 8000 of 8000 bit for bit;
- against the chord of an unposed
  :class:`~orpheus.geometry.chord.ConcentricPartition`, 5320 of 8000 bit
  for bit and the rest within :math:`2.7\,\epsilon c`, because
  ``ConcentricPartition.chord`` first moves the line by the inverse of its
  pose (the identity here), which recomputes the moment;
- against partitions posed by a random rotation and a translation
  :math:`t` (:math:`|t|` from 0 to :math:`2.1 \times 10^{6}`, four poses
  per point, 32 000 lines), within
  :math:`9.5\,\epsilon\,(c + |t|)/|P\Omega|`, the rounding of a moved base
  point (the qa review measured the same growth with :math:`|t|`).

The tangencies: where the chord changes
---------------------------------------

An integrand built from the chord of the line through :math:`x` has a
kink wherever the line becomes tangent to a breakpoint's surface: the
half-chord :math:`h_k = \sqrt{(r_k - b)(r_k + b)}` has an infinite
derivative in :math:`b` at :math:`b = r_k`, and a pair of crossings
appears or disappears. A direction rule at :math:`x` is split there.
:meth:`DirectionDomain.tangencies
<orpheus.geometry.chart.DirectionDomain.tangencies>` returns the values of
the first axis where :math:`b = \ell` for a level :math:`\ell`, sorted:

.. math::

   \text{sphere: } \mu = \pm\sqrt{\bigl(1 - \tfrac{\ell}{c}\bigr)\bigl(1 + \tfrac{\ell}{c}\bigr)},
   \qquad
   \text{cylinder: } \alpha = \arcsin\tfrac{\ell}{c},\ \ \pi - \arcsin\tfrac{\ell}{c},
   \qquad 0 < \ell \le c .

Solving :math:`b(q) = \ell` with the closed forms above gives both. For
:math:`0 < \ell < c` there are two values: directions with :math:`b < \ell`
cross :math:`c = \ell` twice, those with :math:`b > \ell` miss it. At
:math:`\ell = c`, the point's own level (the point on a surface), the two
values merge into the single grazing value, :math:`\mu = 0` or
:math:`\alpha = \pi/2`: every other line through the point crosses the
surface at the point itself, and what changes at the grazing value is the
sense of that crossing, outward for :math:`\Omega_x > 0`, so the backward
characteristic switches the side of the surface it leaves into. The user
ruled on 2026-10-06 that this value is in the break set. Elsewhere the set
is empty: for :math:`\ell > c` every line through the point has
:math:`b \le c < \ell` and crosses; for :math:`\ell \le 0` the level is the
stratum or below it, never crossed; on the slab (no impact parameter; its
one break, the parallel direction :math:`\mu = 0`, is the consumer's own);
and at a stratum, where every line has :math:`b = 0`. On the cylinder the
tangencies do not depend on :math:`w`, so they split the box along
:math:`\alpha` only. `[M]` 2026-10-06: at :math:`c = 1.5`, the level
1.1 gives :math:`\mu = \pm 0.67986927` and :math:`\alpha = 0.82321198,
2.31838068`; the level 1.5 gives :math:`\mu = 0` and
:math:`\alpha = \pi/2` alone.

**Scale-free.** The tangencies depend on :math:`\ell/c` only, and are
computed from the ratio, :math:`\sqrt{(1 - \ell/c)(1 + \ell/c)}`. The form
:math:`\sqrt{(c - \ell)(c + \ell)}/c` underflows or overflows its product:
`[M]` the qa review, at :math:`c = 10^{-200}` it gave :math:`[-0.0]` for
the half level, and at :math:`c \ge 10^{160}` it gave :math:`\pm\infty`.
The kernel's own half-chords still use the unscaled product
(:ref:`chart-and-chord-deferred`, #582).

Grazing at the point's own level is below double precision
----------------------------------------------------------

Near the grazing value at :math:`\ell = c`, a direction :math:`\delta`
off grazing has :math:`b = c\cos\delta`, so :math:`c - b \approx
c\,\delta^2/2`, second order in :math:`\delta`. Below
:math:`\delta \approx \sqrt{2\epsilon} \approx 2 \times 10^{-8}` that is
under an ulp of :math:`c`, and no double-precision :math:`b` resolves the
side: the chord reads a tangency, makes no crossing of the surface and
reports a length of 0 where the exact chord inside the surface is
:math:`2c\sin\delta/|P\Omega|`. `[M]` 2026-10-06, the kernel at
:math:`c = 2` on the partition :math:`(0, 0.3, 1.1, 2.0)`, the in-plane
direction on the cylinder: for :math:`|\delta| \le 3 \times 10^{-9}` both
charts read :math:`c - b = 0`, no transit and length 0 (exact
:math:`1.2 \times 10^{-8}` at :math:`3 \times 10^{-9}`); at
:math:`10^{-8}` the sphere resolves (:math:`c - b = 6.7 \times 10^{-16}`)
and the cylinder does not; where resolved, the length carries the
conditioning of :math:`h` near tangency (:ref:`chart-and-chord-conditioning`),
:math:`1.03 \times 10^{-7}` for the exact :math:`4 \times 10^{-8}` at
:math:`10^{-8}`, and within 1 % at :math:`10^{-7}`. The qa review's probe
found the same band: the correctly rounded :math:`b` equals :math:`c`
for :math:`|\delta| \le 10^{-8}` on both charts. This is the
problem's conditioning, not an algorithm's: the cost is a lost chord of
length about :math:`2c\delta/|P\Omega|` on a band of directions of width
about :math:`\sqrt\epsilon`.

The representative and its refusals
-----------------------------------

:meth:`DirectionDomain.direction
<orpheus.geometry.chart.DirectionDomain.direction>` evaluates
:math:`\omega(q)` on a batch ``(..., k)`` of coordinates. The transverse
component is :math:`\sqrt{(1 - c)(1 + c)}` for a cosine :math:`c`, never
:math:`\sqrt{1 - c^2}`, which loses about :math:`10^{-4}` relative at
:math:`c = 1 - 10^{-12}`; a near-pole representative feeds every grazing
chord, so a lost digit there is a wrong :math:`b`. The call refuses
coordinates of the wrong width, non-finite ones and ones outside the box:
a cosine of 2 would otherwise return the non-unit vector :math:`(2, 0, 0)`,
and a NaN a NaN direction whose impact parameter read :math:`c` (`[M]` the
qa review). :meth:`Chart.directions_at
<orpheus.geometry.chart.Chart.directions_at>` refuses a non-finite orbit
coordinate, and a negative one where the group acts on the kept space (the
cylinder, the sphere); the slab accepts any finite coordinate.

The shape table is a scope boundary
-----------------------------------

The four shapes are written by hand in ``DirectionDomain.shape`` and
``DirectionDomain.stabiliser``, under a ``SCOPE-BOUNDARY[guard]`` tag. The
machinery that would derive them is the point-isotropy computation
:math:`L \cap \mathrm{Stab}(x)` and the orbit-space catalogue of
``orpheus.numerics.manifold`` (:ref:`manifold-orbit-space`), which has
:math:`S^2/O(2)_a` but no :math:`S^2/D_{1h}`. The entry was not built: the
catalogue's lift is the orbit barycentre, the Reynolds projection onto the
group's fixed subspace, which is a right inverse of the quotient map only
when the chart is linear in the ambient coordinates; :math:`D_{1h}`'s
invariants include :math:`\Omega_y^2` and :math:`\Omega_z^2`, its fixed
subspace is the :math:`x` axis, and the barycentre :math:`(\Omega_x, 0, 0)`
forgets the orbit (#581, ruled 2026-10-06). When #581 lands the table
retires onto the catalogue rather than becoming its gate partner, since
two hand-written copies agreeing in a gate would agree by construction.
The gates check the table against :math:`D_{1h}` elements built
independently in the test file.


.. _chart-and-chord-line-domain:

The line domain
===============

The directions at a point are the orbit space a reading at a point
integrates over (:ref:`chart-and-chord-directions`). An operator assembled
over every line of a body, such as the Galerkin block of a characteristic
method (:ref:`characteristic-galerkin-assembly-section`), integrates over
the oriented lines of space instead, against the invariant measure
:math:`\mathrm{d}A_\perp\,\mathrm{d}\Omega` of :eq:`geometry-measure-on-lines`.
The problem's symmetry makes most of that space redundant too: a motion of
:math:`G_c` carries a line to a line with the same chord, so every
integrand built from the chord is invariant under :math:`G_c`, and the
integral over lines is an integral over the orbit space of oriented lines
under :math:`G_c`. :meth:`Chart.line_domain
<orpheus.geometry.chart.Chart.line_domain>` returns that orbit space as a
:class:`~orpheus.geometry.chart.LineDomain`: the lines' counterpart of
:meth:`Chart.directions_at <orpheus.geometry.chart.Chart.directions_at>`,
ruled a kernel verb on ``Chart`` on 2026-10-06 (the plan's ledger, "on P1
step (b)'s third rung", Q2).

The definition
--------------

The domain is a box :math:`B` in at most two coordinates :math:`q`, with a
map :math:`\lambda(q)` to a representative oriented line, such that every
orbit of :math:`G_c` on the oriented lines meets :math:`\lambda(B)` once
(the box's boundary aside) and the invariant measure pushes forward to a
density :math:`\varrho(q)` on the box:

.. math::
   :label: geometry-line-domain

   \int f(L)\,\mathrm{d}A_\perp\,\mathrm{d}\Omega
   \;=\; \int_{B} f\bigl(\lambda(q)\bigr)\,\varrho(q)\,\mathrm{d}q,
   \qquad
   \varrho(q) \;=\; w\bigl(b(q), \Omega(q)\bigr)\cdot\Omega_{\rm fold}(q),
   \qquad
   \text{for every } f \text{ with } f \circ g = f
   \ \text{for all } g \in G_c,

per unit measure of the discarded columns (per unit height on the
cylinder, per unit transverse area on the slab). Here :math:`w` is the
beam density of :eq:`geometry-measure-on-lines`
(:meth:`Chart.beam_density <orpheus.geometry.chart.Chart.beam_density>`)
evaluated on the representative line, and :math:`\Omega_{\rm fold}` is the
measure of the directions that the quotient folds into one point of the
box.

.. implements:: geometry-line-domain
   :by: orpheus.geometry.chart.LineDomain.density

   **Implemented by** ``LineDomain.density``, which multiplies
   ``Chart.beam_density`` of the domain's own representative lines by the
   folded direction measure; ``LineDomain.lines`` is the representative
   map :math:`\lambda` and ``Chart.line_domain`` builds the domain.

.. implements:: geometry-line-domain
   :by: orpheus.geometry.chart.LineDomain.lines

.. implements:: geometry-line-domain
   :by: orpheus.geometry.chart.Chart.line_domain

Unlike the directions at a point, the density is not constant. It is the
invariant measure itself, not divided by :math:`4\pi`, and a consumer that
wants the scalar-flux normalisation divides by :math:`4\pi` itself (the characteristic reference's line weight is the
quadrature weight times :math:`\varrho/4\pi`). The box is unbounded in
:math:`b`: which lines meet a body is the body's question, so a consumer
truncates at its outer radius. The quadrature rule over the box stays the
consumer's, for the reason the beam density's does
(:ref:`chart-and-chord-measure`, "Why the quadrature rule is the
consumer's").

The three shapes
----------------

:class:`~orpheus.geometry.chart.LineShape` names them, decided from the
chart's pair alone (no point is involved): ``COSINE`` where the group fixes
the kept space, otherwise ``IMPACT`` with three kept columns and
``IMPACT_POLAR`` with two.

.. list-table::
   :header-rows: 1
   :widths: 13 22 25 24 16

   * - Chart (shape)
     - Coordinates, box
     - Representative :math:`\lambda(q)`, through, along
     - :math:`w \cdot \Omega_{\rm fold}`
     - :math:`\varrho(q)`
   * - sphere (``IMPACT``)
     - ``impact`` :math:`b \in [0, \infty)`
     - :math:`b\,\hat e_y`, :math:`\hat e_x`
     - :math:`2\pi b \cdot 4\pi`
     - :math:`8\pi^2 b`
   * - cylinder (``IMPACT_POLAR``)
     - ``impact`` :math:`b \in [0, \infty)`, ``polar_angle``
       :math:`\theta \in [0, \pi/2]`
     - :math:`b\,\hat e_y`, :math:`(\sin\theta, 0, \cos\theta)`
     - :math:`2\sin\theta \cdot 4\pi\sin\theta`
     - :math:`8\pi\sin^2\theta`
   * - slab (``COSINE``)
     - ``cosine`` :math:`\mu = \Omega_x \in [-1, 1]`
     - the origin, :math:`(\mu, \sqrt{(1-\mu)(1+\mu)}, 0)`
     - :math:`|\mu| \cdot 2\pi`
     - :math:`2\pi|\mu|`

`[M]` 2026-10-07, the kernel: ``density`` at :math:`b = 0.3` on the sphere
reads 23.6870505626145 (:math:`8\pi^2 \cdot 0.3 = 23.687050562614\ldots`), at
:math:`(b, \theta) = (0.3, 0.7)` on the cylinder 10.430500504411
(:math:`8\pi\sin^2 0.7`), and at :math:`\mu = -0.3` on the slab
1.88495559215388 (:math:`2\pi \cdot 0.3`).

**The sphere: O(3).** An oriented line is a direction :math:`\Omega` and a
foot :math:`p` in the plane normal to :math:`\Omega`. :math:`O(3)` is
transitive on directions, and the stabiliser of a direction is the
:math:`O(2)` of rotations and reflections of the normal plane, whose orbits
there are the circles :math:`|p| = b`. So the complete invariant is
:math:`b`, every direction is folded into one point (:math:`4\pi`), and the
beam density of one direction is :math:`2\pi b`. A line and its reverse are
one orbit: the reflection through the plane normal to :math:`\Omega`
containing the centre reverses :math:`\Omega` and keeps :math:`b`.

**The cylinder: D**\ :sub:`∞h` **with the axial translations.** Write
:math:`\theta` for the angle between the line and the axis and
:math:`\varphi` for the azimuth of :math:`P\Omega`. The rotations about the
axis fold :math:`\varphi` (:math:`2\pi`), the mirror normal to the axis folds
:math:`\theta` and :math:`\pi - \theta` (a factor 2, which also identifies a
line with its reverse, since reversing sends :math:`\theta` to
:math:`\pi - \theta`), the axial translations fold the line's height, and
what is left is :math:`b` and :math:`\theta \in [0, \pi/2]`. With
:math:`\mathrm{d}\Omega = \sin\theta\,\mathrm{d}\theta\,\mathrm{d}\varphi`, the
folded direction measure is :math:`2\pi \cdot 2\sin\theta`, and the beam
density per unit height is :math:`2|P\Omega| = 2\sin\theta`
(:ref:`chart-and-chord-measure`): :math:`\varrho = 8\pi\sin^2\theta`.

**The slab: O(2)**\ :sub:`x` **with the transverse translations.** The
group fixes :math:`\hat e_x`, so the cosine :math:`\mu = \Omega_x` is
invariant and its sign too: a slab line and its reverse are two orbits,
:math:`\mu` and :math:`-\mu`, and the box is :math:`[-1, 1]`. The rotations
about :math:`\hat e_x` fold the azimuth (:math:`2\pi`), the transverse
translations fold the foot, and the beam density per unit area of a plane
is :math:`|\mu|`: :math:`\varrho = 2\pi|\mu|`.

**The representative lines carry their coordinates.** Under
:meth:`Chart.image <orpheus.geometry.chart.Chart.image>`, the
representative's impact parameter is :math:`b` to an ulp, not bit for bit:
a :class:`~orpheus.geometry.line.Line` stores its moment
:math:`p \times \Omega` and returns its foot as :math:`\Omega \times m`,
which rounds (`[M]` 2026-10-06, the test-architect: 1 ulp at :math:`b = 0.3`,
:math:`\theta = 0.2`; exact at :math:`\theta \in \{0, \pi/2\}`). The
cylinder representative's axial cosine is :math:`\cos\theta` and the slab
representative's :math:`\Omega_x` is :math:`\mu`, both bit for bit. A
coordinate outside its bound, a non-finite one (an infinite :math:`b`
included) or a batch of the wrong width is refused with ``ValueError``.
``DirectionDomain`` and ``LineDomain`` share one table of axis bounds and
one validator in ``orpheus/geometry/chart.py``, so the two boxes cannot
disagree on an axis they share.

Why the cylinder's coordinate is the polar angle
-------------------------------------------------

The obvious second coordinate on the cylinder is the axial cosine
:math:`\mu_z = \cos\theta`, the one :meth:`Chart.directions_at
<orpheus.geometry.chart.Chart.directions_at>` uses. In it,
:math:`\mathrm{d}\Omega = \mathrm{d}\mu_z\,\mathrm{d}\varphi`, the folded
measure is the constant :math:`4\pi`, and the density is the beam density
alone, :math:`8\pi\sqrt{1 - \mu_z^2}`: a square-root endpoint at
:math:`\mu_z = 1` (the line parallel to the axis, where the beam's
projected speed vanishes) in every integrand over the domain. Gauss–Legendre
in :math:`\mu_z` then converges algebraically. In :math:`\theta` the
substitution :math:`\mu_z = \cos\theta` absorbs the root:
:math:`\sqrt{1 - \mu_z^2}\,\mathrm{d}\mu_z = \sin^2\theta\,\mathrm{d}\theta`,
an entire function of :math:`\theta`, and Gauss–Legendre converges
geometrically.

The measurement that decided it (`[M]` 2026-10-06, the main agent's
``scratch/characteristic_architecture/p1_step_b3/ta/m12_theta.py``): the
escape probability of a homogeneous white-walled cylinder at
:math:`\tau = 0.5`, against Bickley's closed form, missed by
:math:`5.4 \times 10^{-4}` with 8 Gauss points in :math:`\mu_z`, and by
:math:`2.0 \times 10^{-6}`, :math:`2.3 \times 10^{-9}` and
:math:`1.7 \times 10^{-12}` with 8, 16 and 32 points in :math:`\theta`, in 2
to 5 s; grading :math:`\mu_z` toward 1 instead reached
:math:`7 \times 10^{-11}` at 433 s per block. The test-architect measured the
:math:`\mu_z` rule's algebraic rate on the same escape:
:math:`4.3 \times 10^{-3}`, :math:`5.4 \times 10^{-4}` and
:math:`7.0 \times 10^{-5}` at 4, 8 and 16 points, a ratio of 8 per
doubling. The user ruled the polar angle on 2026-10-06 ("the cylinder's
line coordinate", the plan's ledger). ``directions_at`` keeps
:math:`\mu_z`: at a point the direction measure carries no beam density,
:math:`\mathrm{d}\Omega = \mathrm{d}w\,\mathrm{d}\alpha` is uniform in
:math:`w = |\Omega_z|`, and there is no endpoint to absorb.

.. dropdown:: First got wrong: "the cylinder's axial cosine needs no grading"
   :color: muted

   The third rung's premises were measured on a closed cylinder with a
   mirror wall (``scratch/characteristic_architecture/p1_step_b3/cyl.py``):
   plain Gauss at 8 points in :math:`\mu_z` gave the same conservation,
   :math:`2.6 \times 10^{-10}`, as 16 points or 12 graded layers, and the
   premise read "no grading is needed". Under a mirror each line conserves
   on its own, so the conservation identity is satisfied line by line and
   cannot see the rule over lines at all. The test-architect's verification
   spec moved the question to a white wall, where the escape probability
   reads the rule, and the square-root end showed at once
   (:math:`5.4 \times 10^{-4}` at 8 points). The same blindness of
   conservation to every rule over lines is the hiding mechanism of
   ERR-101.

Cauchy's formula, the independent gate
--------------------------------------

Integrating the chord length through a body over the domain gives
:math:`4\pi` times the body's measure: by :eq:`geometry-measure-on-lines`,
for each direction the lengths of the parallel lines through the body sum
to its measure, and the directions contribute :math:`4\pi`. On the box,

.. math::

   \int_{B} \ell\bigl(\lambda(q)\bigr)\,\varrho(q)\,\mathrm{d}q
   \;=\; 4\pi\,m(\text{body}),

per unit height on the cylinder and per unit area on the slab, with
:math:`\ell` the 3-D chord length inside the body (a hollow body's cavity
is outside it). The gate
``tests/gates/geometry/test_line_domain.py::test_cauchys_formula_on_the_line_domain``
integrates the kernel's own chord lengths with a rule written in the test
(per piece of :math:`b`, the substitution :math:`b = r_j\sin\phi`, which
absorbs each chord's square-root end; Gauss–Legendre in :math:`\theta` and
in :math:`\mu` on each sign; never the reference's line rule) and compares
with :math:`4\pi` times :meth:`Chart.measure
<orpheus.geometry.chart.Chart.measure>` on solid and hollow spheres and
cylinders and on a slab. It shares no formula with the density: the chord
is gated against closed forms in ``test_chord.py`` and the measure is the
one definition, so a wrong fold or a wrong beam factor shows as an O(1)
miss (`[M]` the test-architect's arms: the sphere's fold :math:`2\pi`,
:math:`-0.50`; the cylinder's fold :math:`2\pi`, :math:`-0.50`; the
cylinder's density without :math:`|P\Omega|`, :math:`+0.53`; the slab's
without :math:`|\mu|`, :math:`+13`). Its stabiliser is declared: on the
cylinder the chord carries the obliquity :math:`1/|P\Omega|` and the density
:math:`|P\Omega|`, so their product is blind to a common error in both, and
``test_the_density_is_the_beam_density_times_the_folded_directions`` pins
the density pointwise against the closed forms typed in the test, at L0
under :eq:`geometry-measure-on-lines`. The other rows of the file (the shape
table by hand, the representative lines, the refusals) are ``foundation``;
the Cauchy row waits on this label to move to L0.

The shape table is a scope boundary
-----------------------------------

The table of shapes and folds is written by hand in ``LineDomain.shape``
and ``LineDomain.density`` under a ``SCOPE-BOUNDARY[guard]`` tag, as the
directions' table is: the machinery that would derive it is the orbit
computation of oriented lines under a subgroup of :math:`E(3)`, which does
not exist. When it does, the table and the folds retire onto it.


Verification
============

Nine of the eleven labels on this page are what the kernel's gates under
``tests/gates/geometry/`` name in their ``verifies(...)`` markers (`[M]`
2026-10-07, ``git grep`` of each label's ``verifies`` spelling under
``tests/``); the two minted with the third rung of the characteristic
reference, :eq:`geometry-measure-density` and :eq:`geometry-line-domain`,
are named by none yet, and their rows stay ``foundation`` until the
markers move. The design of the first seven, with every fixture and the
mutation each row must redden, is
``scratch/characteristic_architecture/seed_verification_spec.md``; the
third rung's is ``scratch/characteristic_architecture/p1_step_b3/spec.md``.
Which test carries which label is the generated matrix's to say
(:doc:`/theory/verification/matrix`), so it is not copied here. By file:

- ``test_chart.py``: :eq:`geometry-radial-coordinate`,
  :eq:`geometry-cylinder-axial-factor` and the densities of
  :eq:`geometry-measure-on-lines` at L0 (hand values, with every
  coordinate axis activated and a negative slab coordinate), and the
  group's membership and the singular strata as foundation gates;
- ``test_chord.py``: :eq:`geometry-line-crossing-law`,
  :eq:`geometry-crossing-order`, :eq:`geometry-chord-segment-lengths` and
  :eq:`geometry-cylinder-axial-factor` at L0. Region labels are compared
  exactly and lengths and crossing parameters separately, because every
  concentric chord is symmetric about its closest approach (a reversed
  orientation is invisible to a full-line comparison) and a total length
  cannot see a misplaced interior crossing (it telescopes);
- ``test_line_measure.py``: :eq:`geometry-measure-on-lines` and
  :eq:`geometry-cauchy-mean-chord` at L1, integral identities whose
  agreement is a quadrature converging;
- ``test_line.py``: the line's invariances and refusals, foundation gates
  with no label;
- ``test_chord_transits.py``: :eq:`geometry-transits` at L0, by the hand
  table of walls (every chart, solid and hollow, both sides of every
  tangency, the cylinder at three axial tilts, parallel lines), the same
  table as a mixed ``(2, k)`` batch, and the definition evaluated slot by
  slot over 1500 seeded lines per partition, whose population must draw
  lines with 0, 1 and (on a hollow body) 2 transits; the absent codes are
  a foundation gate;
- ``test_chart_directions.py``: :eq:`geometry-directions-at` at L0: the
  shape table by hand, an invariant integrand over the box against a
  full-sphere rule about :math:`\hat e_z` that shares no coordinate with
  the box (the density is derived from the box's widths, so the
  total-measure check alone is a smoke check), the fundamental-domain
  legs against stabiliser elements built in the test, the round trip of
  the coordinates, the stabiliser by membership, the impact parameter
  against the kernel and the closed form, and the tangencies against the
  closed form in mpmath and against the places where the kernel's
  crossing set changes. The refusals, the conditioning of the
  representative's sine and the scale-free tangencies are foundation
  gates;
- ``test_measure_density.py``: :eq:`geometry-measure-density`, the closed
  forms typed per chart and the density integrated back to the one
  measure under the conditioning factor (:ref:`chart-and-chord-measure-density`),
  ``foundation`` until their markers move;
- ``test_line_domain.py``: :eq:`geometry-line-domain`, the shape table by
  hand, the density against the beam density times the fold at L0 under
  :eq:`geometry-measure-on-lines`, Cauchy's formula on the domain with the
  test's own rule, the representative lines and the refusals
  (:ref:`chart-and-chord-line-domain`);
- ``test_kernel_corroboration.py``: the kernel against today's
  independent spellings, code-to-code agreement, L4, with no correctness
  content. Its value is that a migration which changes an answer shows
  it, and each row asserts by AST that its old spelling does not import
  the kernel, so a migrated row turns red and is deleted.

The gates' constants were re-measured on the kernel (seed 20261005) and
are recorded with their probes in
``scratch/characteristic_architecture/seed_gates/README.md`` (its "v2"
section is the kernel as it stands; where it and a test file disagree, the
test file is the value). `[M]` 2026-10-05 (the test-architect, on the
restructured kernel; 113 rows in the five files):

.. list-table::
   :header-rows: 1
   :widths: 46 26 28

   * - Quantity
     - Measured
     - Tolerance the gate uses
   * - invariance under :math:`G_c` and under a pose: a slot's change over
       :math:`\epsilon R(|p| + R)/(h_{\min}|P\Omega|)`
     - 4.30 (cylinder, group); 3.15 (cylinder, pose)
     - 50 (``_INVARIANCE_C``)
   * - a crossing's orbit coordinate against :math:`r_k`, over
       :math:`\epsilon(|\text{foot}| + |t| + 4)/|P\Omega|`
     - 1.25, 0.96, 0.78
     - 16 (``_SURFACE_C``)
   * - the kernel's exit distance against the Variant-α, Peierls and MoC
       spellings (L4), over :math:`\epsilon R(r + R)/h`
     - 1.80
     - 16 (``_JOIN_C`` in ``test_kernel_corroboration.py``)
   * - multi-region segments against the Variant-α oracle (L4)
     - 0 of 1000 region sequences differ; largest length difference
       :math:`5.1 \times 10^{-15}`
     - absolute :math:`10^{-13}`
   * - slot lengths against mpmath
     - within 4 ulp (8 with the obliquity)
     - 4 / 8 ulp
   * - Cauchy's mean chord, 32-point direction rules
     - :math:`3.3 \times 10^{-16}` (cylinder), :math:`1.6 \times 10^{-15}`
       (slab)
     - :math:`2 \times 10^{-14}`

Each gate's teeth were shown by a mutation battery over the five files
(``seed_gates/battery/``, every mutation in-process, in the kernel's own
algebraic class): `[M]` 2026-10-05, on the restructured kernel, 45 of 45
arms red their target row, and the positive control (every shell segment
scaled by :math:`1 + 10^{-7}`) reds 39 rows, more than any arm (the
largest, the cylinder image reading :math:`\Omega_z` for
:math:`\Omega_y`, reds 34; the crossing-order rule shifted inward reds
33, because every slot region now derives from it). The first battery
found four kernel defects before the review (a NaN direction admitted, a
NaN coordinate located as region :math:`n`, an infinite direction warning
before its refusal, and chord parameters measured on the canonical line
under a pose), and the qa review's own battery found two gaps it closed
(a half-line starting in the cavity; a half-line on a parallel line);
each is repaired and gated. The L4 harness also carries a runtime
independence leg: counting spies on the kernel's entry points must read
zero calls while an old spelling runs, because the AST precondition
misses an indirect import (`[M]` the qa review: 4 of 4 indirect shapes
missed).

The two verbs' gates (`[M]` 2026-10-06, the test-architect's battery
``scratch/characteristic_architecture/p1_step_a/battery/``, each arm an
in-process mutation of the kernel, run over both files): every arm but one
reds at least one row. The transit arms: an untraversed exterior slot
breaking a run (50 rows red), walls read by slot position (64), a parallel
line taken as a transit (8), the absent slot code :math:`n + 1` (77), the
exit wall read at the opening of the last slot (66). The direction arms:
the density's 4 dropped (14), the cylinder's stabiliser taken as
:math:`O(2)_x` (4), :math:`|P\Omega|` as :math:`\sqrt{1 - \Omega_z^2}`
(2), the sine as :math:`\sqrt{1 - c^2}` (4), the reflection
:math:`\pi - \arcsin` dropped (8), the grazing value dropped (8),
:math:`\alpha` on :math:`[0, \pi/2]` (7), a negative coordinate accepted
(2), :math:`\arccos` for :math:`\arcsin` (8), the slab's refusal removed
(1). The positive controls (no transit at all; the sine replaced by the
cosine) red 67 and 25 rows. The one arm green in both files, a tangency
counted as a crossing (:math:`b \le r_k`), is designed-green here: the
crossing pair it adds bounds a slot of length 0, which a transit ignores
by definition; ``test_chord.py`` catches it (`[M]` the qa review, 4 rows
red). The qa review added three rows its probes found missing (the
impact parameter at :math:`w = 1`, the box refusal, the scale-free
tangencies). 168 rows: 89 for the transits, 79 for the directions.

.. _chart-and-chord-deferred:

What the kernel does not do, and what still computes the same thing
===================================================================

**One consumer, a reference.** The characteristic reference
(:ref:`theory-characteristic-reference`) is the only module outside the
three that calls the kernel (Key facts); it is built to replace the
trajectory-resolvent family and is consumed by no production code yet.
The reference family's characteristic oracles, the
collision-probability chords, the method of characteristics and Monte
Carlo each compute chords, crossings, regions and line measures their own
way. `[M]` 2026-10-05, the geometry census over the 385 tracked
``orpheus/**/*.py`` files at ``a336bde4``
(``scratch/characteristic_architecture/geometry_census.md``):

- 13 square-root-of-a-difference-of-squares chord sites and 14 geometric
  discriminants outside the SymPy origins, almost all in reference code
  (production holds two, both in MoC). The census's probe F1 found the
  Variant-α, Peierls and MoC distance-to-surface spellings agreeing to
  :math:`7 \times 10^{-14}` relative over 1000 draws with no shared
  import: independent agreement on one quadratic;
- 8 point-location spellings with two boundary conventions
  (:ref:`chart-and-chord-location`);
- 3 realisations of the measure on lines (``chord_quadrature``, the
  collision-probability :math:`y`-quadrature, MoC's track spacing and
  weights).

Moving the reference family onto the kernel is tracked under #405 (the
slow reference tier, whose chord oracles opened this design); the Monte
Carlo distance to boundary is #534. A comparison of a migrated spelling
with the kernel is the kernel compared with itself through a facade
(`retirement-audit` D.14), so it is deleted when its spelling migrates.

**Not supported, each a named next member of the kernel:**

- a Boolean (constructive solid geometry) partition. The crossing order
  gives a region only when the partition is pulled back from a 1-D orbit
  space; a Boolean partition is not, and needs a membership query per
  segment. The root solver would be shared; the region rule would not.
  The user ruled on 2026-09-24 that general solids come from external
  meshers;
- the finite-volume metrics of unstructured meshes (centroids,
  centroid-to-face distance, non-orthogonality): #322, #335, #539;
- the :math:`(r, z)` and 2-D charts, the deck group and the boundary laws
  read through the chart: #551;
- the directions at a point derived from the group: the shape table of
  :class:`~orpheus.geometry.chart.DirectionDomain` is written by hand, a
  declared scope boundary that retires onto the orbit catalogue's
  :math:`S^2/D_{1h}` entry (#581, :ref:`chart-and-chord-directions`);
- the line domain derived from the group: the shape table and the folds
  of :class:`~orpheus.geometry.chart.LineDomain` are written by hand too,
  a declared scope boundary that retires onto an orbit computation of
  oriented lines under a subgroup of :math:`E(3)`
  (:ref:`chart-and-chord-line-domain`);
- a third spelling of the measure density in production,
  ``compute_areas_1d``, whose retirement onto
  :meth:`CoordSystem.measure_density <orpheus.geometry.coord.CoordSystem.measure_density>`
  moves bits on the sphere in S\ :sub:`N`
  (`#584 <https://github.com/deOliveira-R/ORPHEUS/issues/584>`_,
  :ref:`chart-and-chord-measure-density`).

**A limit at extreme radii (#582).** The chord's half-chord
:math:`\sqrt{(r_k - b)(r_k + b)}` is formed unscaled, so its product
underflows or overflows far from unit radii. `[M]` the qa review
(``scratch/characteristic_architecture/p1_step_a/qa/probe2.log``): a solid
sphere of radius :math:`10^{-200}` reports a line through its centre with
total length 0 and no transit, one of radius :math:`10^{-160}` reports
:math:`1.99999 R` for :math:`2R`, and one of radius :math:`10^{160}`
overflows. The tangencies of :ref:`chart-and-chord-directions` are
already scale-free; the chord is not. Reactor radii are far from these
limits.


.. _chart-and-chord-refuted:

Designs that do not work, and why
=================================

Each row is a design that reads as natural and is wrong, with the
structural reason it fails. The first eight were in the first-pass
design and were removed by the design review before any code (the
cross-domain attacker and the elegance enforcer,
``scratch/characteristic_architecture/w5_cross_domain.md`` and
``w5_elegance.md``); the last four were in the first built kernel and
were removed by the review of the code (qa, ``seed_qa.md``; the elegance
enforcer, ``seed_elegance.md``); the last three belong to the transits
and the directions at a point (2026-10-06). They are kept so that no
later design re-derives them.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - The design
     - Why it fails
   * - "One crossing law :math:`\rho^2 = r_k^2` for the three charts"
     - With the slab's signed coordinate, squaring admits the mirror
       roots: on breakpoints :math:`(-1, 0.5, 2)` it finds 6 roots of
       which 3 are spurious (:math:`x = -2, -0.5, 1`) (`[M]` both
       reviewers). The slab's invariant has degree 1, the curved charts'
       degree 2; the orbit-space law (:eq:`geometry-line-crossing-law`
       and the monotone slab image) is one formula per degree, not one
       squared law.
   * - "The impact parameter is the minimum of the coordinate along the
       line", for every chart
     - On the slab the coordinate is affine along the line and has no
       minimum unless the line is parallel; the slab has no impact
       parameter (its image is an :class:`~orpheus.geometry.chart.AxialImage`,
       which has none).
   * - "The invariant measure on lines in the plane, :math:`b` times the
       line's angle"
     - It is the 2-D density. On the sphere, whose lines of one
       direction fill a disc, it drops the factor :math:`b`: Cauchy's
       mean chord reads :math:`\pi R/2` (the disc's) instead of
       :math:`4R/3`; on the cylinder the density must carry
       :math:`|P\Omega|`, without which the mean chord is
       :math:`\pi^2 R/4` instead of :math:`2R` (`[M]` the cross-domain
       review and this page's table).
   * - "Grazing goes to the outer side, derived, not conventional"
     - True on the sphere only. On a slab line parallel to the planes, or
       a cylinder line along the axis, lying in :math:`c = r_k`, every
       term of :math:`c(t) - r_k` vanishes: the line stays in the
       interface for its whole length and no order of the expansion
       picks a side (`[M]` the cross-domain review). Hence the typed
       interface.
   * - "``locate(point, direction)``" as the directional locator
     - A computed point is never exactly on a surface, so a function of
       the point cannot know it is on one; the datum belongs to the step
       that produced the point, the crossing (the elegance review, item
       3). The cross-domain review's second pass refuted its own first
       draft's shortcut, "the directional answer is the ray's first
       segment": a crossing computed at :math:`t = +10^{-17}` makes a
       sliver segment and picks the wrong side.
   * - "Solve each 3-D line's own quadratic" (leading coefficient
       :math:`1 - \mu_{\text{axial}}^2`)
     - It re-solves the in-plane chord once per axial cosine, discarding
       the fibration the hoist measured bit-identical; it computes the
       leading coefficient as :math:`1 - \Omega_z^2` (no correct digit at
       :math:`|P\Omega| = 10^{-8}`); and it divides by zero on a line
       parallel to the axis.
   * - "Inner-owns, closed on :math:`[0, R]`"
     - On a hollow body this puts the cavity in region 0, a region with
       a material; the domain is :math:`[r_0, r_n]` and the cavity is
       the typed inner exterior (the elegance review, item 4).
   * - "Point" and "Direction" as per-value classes
     - Every consumer is batched (quadratures ``(N, 3)``, a thousand to
       a million oracle samples, MoC tracks), so a per-object value
       forces a Python loop; the values are arrays, and a direction is a
       point of the existing :class:`~orpheus.numerics.manifold.Sphere`.
   * - Exterior codes as negative integers (:math:`-1` for the inner
       exterior, :math:`-2` for the outer)
     - Numpy reads a negative index as "count from the end", so the
       natural idiom ``sigma_t[chord.slot_region]`` gives a hollow body's
       cavity the cross section of the outermost region, and nothing
       raises. `[M]` the qa review: hollow sphere :math:`(0.5, 1, 2)`,
       a line through the centre, :math:`\Sigma_t = (1, 2)`, optical
       depth 7.0 where the void cavity gives 5.0; the elegance review,
       :math:`\Sigma_t = (1, 3)`, 10.0 where 7.0 is right. The interface
       code :math:`-1` also equalled the inner-exterior code. Codes out of
       range (:math:`n`, :math:`n + 1`, and :math:`n + 1` for "no
       interface") make the same idiom raise.
   * - A chart matched on the coordinate system (one ``match`` per verb)
     - Five ``match`` blocks plus a cylinder-or-sphere branch in the chord
       spelled, per verb, what two data determine: the kept columns and
       the linear group. A symmetry group is realised, never tabulated:
       the derived chart asks the group's realization, so a new chart
       (the :math:`(r, z)` one of #551) adds a pair, not an arm in every
       verb. `[M]` the elegance review: the matched ``contains`` equalled
       the realization test on 750 of 750 motions per chart.
   * - Slot starts and slot regions built beside the crossings
     - They were the crossings computed a second time by independent
       expressions (`[M]` the elegance review: equal to 0.0 on 4 of 4
       fixtures), and the crossing-order rule was written twice. Any edit
       to one (dropping the never-crossed :math:`r_0 = 0` from the
       crossings only) would break slot :math:`i` matching crossing
       :math:`i` silently. The slots now read the crossings.
   * - "The orbit coordinate moves at the projected speed", and a
       ``None`` field to say "slab"
     - :math:`|P\Omega|` is the speed in the kept space; on the sphere
       :math:`\mathrm{d}|x|/\mathrm{d}t` runs over :math:`[-1, 1]`
       while :math:`|P\Omega| = 1`. And the slab-or-curved question,
       decided once, was re-asked through ``None`` fields
       (``impact_parameter``, ``closest_approach``), so a chord rebuilt
       with one of them set to ``None`` was silently read as a slab. The
       image types answer it once.
   * - A transit split at every exterior slot, traversed or not
     - A solid body's chord has a cavity slot too (code :math:`n`, length
       0), so every line with :math:`b < r_1` would read two transits, and
       a hollow body at :math:`b = r_0` exactly, a tangency, would read
       two. `[M]` the battery's arm T1: 50 rows red. Only a traversed
       exterior slot separates.
   * - A wall named by its position (the entry and exit of a transit, the
       first and second transit)
     - The same position is a different wall in different cases (a solid
       body enters and leaves by wall :math:`n`; a hollow body's first
       exit and second entry are wall :math:`0`), and no branch on
       :math:`b` against :math:`r_0` is needed once the wall is its
       breakpoint index (the W5 elegance review, E5: the transits read off
       the P0 chord with no inner or outer tag). `[M]` arm T2, the slot
       index read as the wall: 64 rows red.
   * - The :math:`S^2/D_{1h}` entry of the orbit catalogue, lifted by the
       orbit barycentre, as the cylinder's direction domain
     - The barycentre is a right inverse of the quotient map only for a
       chart linear in the ambient coordinates; :math:`D_{1h}`'s fixed
       subspace is the :math:`x` axis, so the barycentre
       :math:`(\Omega_x, 0, 0)` forgets :math:`\Omega_z^2` and does not
       identify the orbit (#581). The kernel's box is the section the
       catalogue lacks.

Frames the cross-domain review found and refuted for the kernel (each
refuted for the question "does the kernel need it?", not as
mathematics):
the geodesic and Christoffel reading (lines are straight; it is decisive
for the curvilinear S\ :sub:`N` streaming term, not here); tensor
networks (the closure rank of a line is at most 2); a Boolean algebra of
half-spaces, bounding-volume hierarchies and non-quadric surfaces (no
consumer, and the user's ruling of 2026-09-24 that general solids come
from external meshers); chain complexes and the discrete Hodge star (no
cells; decisive for the finite-volume metrics below); an invariant-theory
engine (three charts, a hand catalogue suffices); the Volterra operator
(the attenuated line integral belongs to the reference family, not to
geometry).


.. _chart-and-chord-gotchas:

Gotchas
=======

- **A per-region table indexed by** ``slot_region`` **raises on an
  exterior.** The exterior codes :math:`n` and :math:`n + 1` are out of
  range for an :math:`n`-region table on purpose; extend the table with
  the exteriors' values (``np.append(sigma_t, [0.0, 0.0])`` for a void
  cavity and outside) rather than masking.
- **There is no** ``impact_parameter`` **on the chord.** It is
  ``chord.image.impact_parameter``, present only when the image is a
  :class:`~orpheus.geometry.chart.RadialImage`; the half-chords are
  computed from it, not stored.
- **The region of closest approach occupies two slots** on the cylinder
  and the sphere (inbound and outbound, split at :math:`t^*`), so a
  consumer that counts segments by slot sees it twice; the verification
  spec's "five segments" through the solid sphere's centre are six slots
  here.
- **Parameters are measured from the foot of the caller's line**, not
  from the base point the line was built from and not in the partition's
  canonical frame: ``slot_start``, the image's ``closest_approach`` and
  the crossing parameters are :math:`t` with :math:`x = \text{foot} + t\,\Omega` on
  the line passed to ``chord``, whatever the pose. Convert a base point
  with :meth:`Line.parameter_of <orpheus.geometry.line.Line.parameter_of>`.
- **A parallel line is detected by exact zero.** :math:`|P\Omega| = 0.0`
  only; a line one ulp off parallel has a finite chord of enormous length
  and no ``parallel`` flag, and a line one ulp off a surface is inside a
  region, not on the interface.
- **The slab and the cylinder treat a parallel line outside the domain
  differently.** A cylinder line along the axis inside a hollow body's
  cavity fills the cavity slot with an infinite length; a slab line
  parallel to the planes below :math:`r_0` has no slot to fill (the slab
  has no inner-exterior slot) and its chord is all zeros (`[M]`
  2026-10-05).
- **A tangent line is not on the surface.** :math:`b = r_k` makes no
  crossing, and the line is in the outer region on both sides of the
  touch; only a line lying in a surface (:math:`|P\Omega| = 0` at a
  breakpoint) reports an ``interface``.
- **A line is not a value-equality key** (``eq=False``;
  :ref:`chart-and-chord-lines`).
- **A transit's slice holds untraversed slots.** ``first_slot:stop_slot``
  runs from the first traversed slot to the last, and every slot between
  them is in it, including regions the line misses and a solid body's
  cavity slot, all of length 0. Sum lengths over the slice; do not count
  its slots as segments.
- **Walls are breakpoint indices, and absent ones are out of range.**
  ``entry_wall`` and ``exit_wall`` are :math:`0` or :math:`n` for a
  present transit and :math:`n + 1` for an absent one; an absent
  ``first_slot`` is the slot count. Index per-wall tables with them and
  mask with ``present``, never with a sentinel test of your own.
- **The density is not the local Jacobian.** On the cylinder off its axis
  :math:`\mathrm{d}\Omega = \mathrm{d}w\,\mathrm{d}\alpha` locally,
  and ``density`` is 4, because the box is a quarter of the sphere and an
  invariant integrand is counted four times. It is correct only for an
  integrand invariant under the point's stabiliser.
- **The cylinder's impact parameter at** :math:`w = 1` **is** :math:`c`,
  not :math:`c\sin\alpha`: the line is parallel to the axis. A rule
  with a node at :math:`w = 1` sees the parallel line.
- **Grazing at a point on a surface is unresolvable below about**
  :math:`10^{-8}` **rad.** Within that band of the grazing direction the
  chord reads a tangency and reports no length inside the surface
  (:ref:`chart-and-chord-directions`); a direction rule should not put
  nodes there.
- **The reading's** :math:`b` **equals** ``Chart.image`` **bit for bit, not
  a posed partition's chord.** A
  :class:`~orpheus.geometry.chord.ConcentricPartition` moves the line by
  its pose's inverse first, so its :math:`b` differs by a few ulp even
  under the identity pose.
- **The line domain's density is not divided by** :math:`4\pi`. It is
  the invariant measure :math:`\mathrm{d}A_\perp\,\mathrm{d}\Omega` itself;
  a consumer computing a scalar flux divides by :math:`4\pi`
  (:ref:`chart-and-chord-line-domain`).
- **The cylinder's line coordinate is the polar angle, its directions'
  coordinate the axial cosine.** ``line_domain`` uses :math:`\theta`,
  where the density :math:`8\pi\sin^2\theta` is analytic;
  ``directions_at`` uses :math:`w = |\Omega_z|`, where the direction
  measure is uniform. A rule built for one is not a rule for the other.
- **A representative line's** :math:`b` **is its coordinate to an ulp,
  not bit for bit.** The kernel's line returns its foot from its moment,
  which rounds; compare impact parameters with a tolerance.
- **The measure density is a density, and its integral over a narrow
  cell far from 0 is ill-conditioned.** Compare a density integrated over
  :math:`[a, b]` with the measure at a tolerance carrying
  :math:`1 + \max(|a|, |b|)/(b - a)`, never at a bare ulp count.


Development history
===================

.. list-table::
   :header-rows: 1
   :widths: 14 58 14 14

   * - Date
     - Decision
     - Commit
     - Issue
   * - 2026-10-05
     - The geometric kernel seed: :mod:`orpheus.geometry.chart`,
       :mod:`orpheus.geometry.line`, :mod:`orpheus.geometry.chord`. The
       user ruled that the Variant-α chord oracles are not patched further
       ("check if there is a better way to architect the general idea
       before patching old code") and that they are "the seed of something
       bigger"; the design was reviewed (W5) before any code and revised on
       the rulings in ``.claude/plans/characteristic_reference_architecture.md``:
       crossings solved once in the orbit space, lines in Plücker
       coordinates, regions from the crossing order, the typed interface,
       the pose, and architecture E's 1-D ``Chart`` pulled forward.
     - ``fe0ca696``
     - #405, #551
   * - 2026-10-05
     - The review of the built kernel (qa, ``seed_qa.md``; the elegance
       enforcer, ``seed_elegance.md``) and the user's rulings on it. The
       exterior codes, negative integers that numpy read as valid
       indices (a hollow sphere's cavity took the outermost region's cross
       section: optical depth 7.0 for 5.0), became the out-of-range codes
       :math:`n` and :math:`n + 1`, with "no interface" :math:`n + 1`
       ("Out-of-range codes n, n+1"). The chart is derived from (kept
       columns, linear group) and returns the line's image as a
       ``RadialImage`` or an ``AxialImage`` ("Derive it now"). The slots
       read the crossings, so the crossing-order rule has one home. The
       cylinder's beam density is over :math:`b \ge 0`. :math:`|P\Omega|`
       is a scaled norm (an underflowing direction had given
       :math:`b = 5.6 \times 10^{-6}` for 0); a mixed batch of parallel and
       crossing lines no longer evaluates :math:`-\infty + \infty`; a
       negative distance coordinate is refused; ``Line.moved_by`` states
       its parameter shift.
     - ``fe0ca696``
     - #405, #551
   * - 2026-10-06
     - The two verbs P1 of the characteristic references reads first,
       each by ruling ("Transits: in the kernel", "Directions at a point:
       chart verb"): :attr:`Chord.transits <orpheus.geometry.chord.Chord.transits>`,
       with walls by breakpoint index and out-of-range absent codes, and
       :meth:`Chart.directions_at <orpheus.geometry.chart.Chart.directions_at>`,
       the box of :math:`S^2/\mathrm{Stab}(x)` with its density, the
       kernel's impact parameter and the scale-free tangencies, the grazing
       value included. The :math:`S^2/D_{1h}` catalogue entry was not built
       (#581); the chord's own scale limit is #582.
     - ``2b2d7703``
     - #405, #581, #582
   * - 2026-10-07
     - Two verbs for the characteristic reference's third rung, each by
       ruling ("Q2: the line domain is a kernel verb on Chart"; "Q4: the
       volume density moves to the kernel"):
       :meth:`CoordSystem.measure_density <orpheus.geometry.coord.CoordSystem.measure_density>`,
       the derivative of the one measure and the area of a level set, onto
       which the reference's ``PanelBasis.volume_density`` retired; and
       :meth:`Chart.line_domain <orpheus.geometry.chart.Chart.line_domain>`,
       the oriented lines modulo :math:`G_c` with their invariant density,
       gated by Cauchy's formula. The cylinder's second coordinate is the
       polar angle (the user's ruling of 2026-10-06, after the axial
       cosine's square-root end cost the escape probability
       :math:`5.4 \times 10^{-4}` at 8 points). The production twin
       ``compute_areas_1d`` is #584.
     - ``71a207fa``
     - #405, #584
