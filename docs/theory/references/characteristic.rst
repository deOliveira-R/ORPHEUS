.. _theory-characteristic-reference:

=============================================================================
The characteristic reference — transport along the lines of a concentric body
=============================================================================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: derivations
      concept: characteristic reference, boundary resolvent, walls, line period, line closure, panel basis, traversal integrals, Volterra block, hp grading
      role: "the closed reference that integrates transport along the lines of a 1-D concentric body (slab, cylinder, sphere, solid or hollow); this page holds its walls (each boundary point with what its law returns, read from the law's factors), the line part of its boundary resolvent (the period of each line's unfolded path and the least solution of its cycle), its panel basis and the transport along one line on that basis (the traversal integrals, the vacuum Volterra block and the angular flux)"
      code: [orpheus.derivations.continuous.characteristic.walls, orpheus.derivations.continuous.characteristic.closure, orpheus.derivations.continuous.characteristic.basis, orpheus.derivations.continuous.characteristic.transport]
      depends_on: [chart_and_chord, boundary_conditions, reference_solutions]
      related: [trajectory_resolvent, layering]


Key facts
=========

- **What this is.** The characteristic reference is the closed reference
  that solves transport on a 1-D concentric body by integrating along
  the body's lines: :mod:`orpheus.derivations.continuous.characteristic`.
  It is built rung by rung beside the trajectory-resolvent family
  (:ref:`theory-trajectory-resolvent`), the family it is built to
  replace. Two rungs exist. The first is the **walls**
  (:mod:`~orpheus.derivations.continuous.characteristic.walls`) and the
  **line part of the boundary closure**
  (:mod:`~orpheus.derivations.continuous.characteristic.closure`); the
  second is the **panel basis**
  (:mod:`~orpheus.derivations.continuous.characteristic.basis`) and the
  **transport along a line** on it
  (:mod:`~orpheus.derivations.continuous.characteristic.transport`).
  Nothing in the package integrates over the measure on lines: there is
  no assembly, no diffuse part of the resolvent, no question and no
  reading at a point, so the package answers no eigenvalue and no flux
  (:ref:`characteristic-what-is-not-built`).
- **The boundary resolvent has two parts.** The closure of a reflecting
  boundary is :math:`P = P_0 + E\,(I - T)^{-1} X` on the boundary trace
  space. On every wall whose return is specular (a mirror, a partial
  mirror, vacuum as amplitude 0) or a periodic wrap, the returned path is
  a line congruent to the one that left, so :math:`T` is diagonal over
  lines: that is :class:`~orpheus.derivations.continuous.characteristic.closure.LinePeriod`.
  A diffuse wall couples every line to every other and is the second
  part, a finite-rank update over the walls, which is not built
  (:ref:`characteristic-resolvent`).
- **A wall is read from its law's two factors**, the deck
  (``geometry_map``) and the response (``response_kernel``), by one table
  (:ref:`characteristic-walls-factors`). It is keyed by its breakpoint
  index (:math:`0` or :math:`n`), never by a position, and carries a
  specular amplitude, a diffuse amplitude and a partner (the breakpoint
  at which the returned path re-enters).
- **The reference parses tags with its own registry**,
  :data:`~orpheus.derivations.continuous.characteristic.walls.TAG_REGISTRY`,
  never with production's parse (the user's ruling of 2026-10-06, the
  insulation principle of :ref:`architecture-reference-insulation`). Each
  tag kind takes exactly its declared parameters: a missing or an
  undeclared one is refused, never dropped
  (:ref:`characteristic-walls-registry`).
- **The period of a line is derived, never tagged.** A traversal is a
  transit (:ref:`chart-and-chord-transits`) read forward or reversed; the
  next traversal is the one entering at the partner of the wall the last
  one left, the forward one first. The period closes after :math:`m`
  steps, and :math:`m` is the rank (:eq:`characteristic-transit-rank`):
  1 on a solid body, on a shell's line that misses the cavity (the
  tangency :math:`b = r_0` included) and on the periodic slab; 2 on a
  shell's line through the cavity and on a slab between two returning
  walls; 0 on a line that meets no wall (:ref:`characteristic-period`).
- **The line closure is the least solution of a cycle**
  (:eq:`characteristic-closure`): the inflow to traversal :math:`k + 1`
  is the specular amplitude of the wall traversal :math:`k` exits at,
  times what traversal :math:`k` carries to its exit. One rolled
  expression serves every rank; :math:`1 - \Pi` is formed by ``expm1``;
  a lossless trapped line (:math:`\Pi = 1`) carries exactly 0 when it has
  no source and is refused when it has one (:ref:`characteristic-closure-section`).
- **The panel basis** is discontinuous: nodal Lagrange polynomials of
  degree :math:`p` through each panel's Gauss–Legendre points, so a
  coefficient is a value. The panels are graded geometrically (panel ends
  at the depths :math:`w\rho^{j}`) toward every wall and interface and
  never toward a singular stratum, and the mass matrix is exact in the
  chart's volume measure, whose density :math:`\kappa\,d\,r^{d-1}` is
  derived from the measure's one definition (:ref:`characteristic-panel-basis`).
- **The transport along a line is one value**,
  :class:`~orpheus.derivations.continuous.characteristic.transport.TraversalRule`,
  built from the lines, the basis and the walls. The kernel's chord through
  the panel partition gives every piece of every line, and one graded
  attenuated integral gives the traversal integrals :math:`\tau_k`,
  :math:`B_k` and :math:`A_k` (:eq:`characteristic-traversal-integrals`),
  the vacuum Volterra block and the angular flux on a line
  (:ref:`characteristic-transport`).
- **The pieces are graded twice.** Exponentially toward both ends of each
  slot, at 1 to 64 mean free paths; and on a cylinder or a sphere by
  halving toward the complex branch points of the orbit coordinate, which
  lie :math:`c/|P\Omega|` from a point of orbit coordinate :math:`c`, so
  that every piece converges at a rate independent of the impact parameter
  (:ref:`characteristic-branch-grading`).
- **Evidence** `[M]` 2026-10-06: for the walls and the closure, 175 gate
  rows in two files, 127 of them claims on this page's first two labels at
  L0, and a 34-arm mutation battery with a positive control reddening 80
  rows; for the basis and the transport, 328 rows in two more files, all
  ``foundation`` until their markers move under
  :eq:`characteristic-traversal-integrals`, and a battery in which every
  arm meant to redden reddens its target row. All 503 rows pass under
  ``python -O -m pytest`` (:ref:`characteristic-evidence`).


.. _characteristic-place:

Where this reference sits
=========================

The trajectory-resolvent family (:ref:`theory-trajectory-resolvent`) is
seven oracle classes and fifteen solver entry points (the plan's
inventory; ``git grep`` counts the fifteen), one per geometry and
boundary shape, each with its own chord oracle and its own closure;
the user ruled on 2026-10-05 that it is re-architected before it is
patched further (the plan
``.claude/plans/characteristic_reference_architecture.md``, "The ruling
that opened this plan"). The characteristic reference is the result: one
reference posed from the specification, with the geometry read through
the geometric kernel (:ref:`theory-chart-and-chord`: charts, lines,
chords and transits) and the linear algebra through the reference
kernel (:ref:`verification-reference-kernel`). It is a closed reference
in the sense of :ref:`architecture-reference-insulation`: it imports
numpy, scipy, the geometric kernel, the boundary-law declarations of
``orpheus.geometry.boundary`` and the reference kernel, and no
production machinery.

The name says what the method integrates along, the characteristic
lines. During the migration the two packages coexist under different
names, ``characteristic/`` and ``trajectory_resolvent/``, so neither is
a homonym of the other.


.. _characteristic-resolvent:

The boundary resolvent and its two parts
========================================

Write the transport problem on the body with a reflecting boundary as a
problem on the boundary trace space. :math:`P_0` is the vacuum transport
operator (each line integrated from the wall it enters at, nothing
entering); :math:`X` takes a volume emission to the outgoing intensity
at the walls; :math:`E` takes an incoming intensity at the walls to the
flux it produces inside; :math:`T` is the boundary laws' response
composed with their deck map (where the returned intensity re-enters)
and with the wall-to-wall transmission. The closure is

.. math::

   P \;=\; P_0 \;+\; E\,(I - T)^{-1} X ,

the resolvent of the boundary map, ruled as the closure of this
reference on 2026-10-06 (the plan's ledger, "the user, on the W5 and
spec questions (first batch)": "Extend now (Krein form)").

The resolvent splits by what a wall returns. On a wall whose return is
**specular** (a mirror, a polished wall returning a fraction specularly,
vacuum as the fraction 0) or a **periodic wrap**, the returned path is a
line again, congruent to the line that left under the chart's symmetry
group: on the sphere and the cylinder a reflection keeps the impact
parameter :math:`b` and the projected speed :math:`|P\Omega|` (the
cylinder also keeps :math:`\Omega_z`), on the slab it reverses
:math:`\Omega_x`, and a wrap translates the line by the slab's width.
So :math:`T` restricted to these walls maps each line to itself: it is
diagonal over lines, and its block on one line is a cycle of at most two
traversals. That block and its least solution are
:class:`~orpheus.derivations.continuous.characteristic.closure.LinePeriod`
(:ref:`characteristic-period`, :ref:`characteristic-closure-section`).

A **diffuse** wall returns its outflow isotropically, so every line
leaving it feeds every line entering it. On the walls :math:`W` with a
diffuse amplitude the resolvent is a finite-rank update of the
line-diagonal block, ruled on 2026-10-06 ("One resolvent, two parts"):
:math:`K = K_{\rm line} + U\,(I - T_w)^{-1} A\,U^{\mathsf T}`, with
:math:`U` the escape functional from emission to the outgoing partial
current at each diffuse wall, :math:`T_w` the wall-to-wall transmission
with every diffuse wall treated as absorbing in the line part, and
:math:`A` the diffuse amplitudes. The walls already carry the diffuse
amplitude (:ref:`characteristic-walls`); the update is not built
(:ref:`characteristic-what-is-not-built`).


.. _characteristic-walls:

The walls
=========

A **wall** is a boundary point of the domain :math:`[r_0, r_n]` of a body
with :math:`n` regions, together with the three numbers the closure
needs from its boundary law. It is the frozen value
:class:`~orpheus.derivations.continuous.characteristic.walls.Wall` with
four fields:

- ``breakpoint``: the breakpoint index of the boundary point, :math:`0`
  (the inner wall: a shell's cavity surface, or the slab's left face) or
  :math:`n` (the outer wall). A solid sphere or cylinder has one wall,
  at :math:`n`; its centre is a stratum of the chart, not a surface
  (:ref:`chart-and-chord-strata`);
- ``specular``: the amplitude, in :math:`[0, 1]`, returned along the
  reflected line (or, across a wrap, along the translated line);
- ``diffuse``: the amplitude, in :math:`[0, 1]`, returned isotropically
  (Lambertian);
- ``partner``: the breakpoint at which the returned path re-enters, the
  wall itself for a mirror and for any re-emission, the opposite wall
  for a wrap.

A wall is keyed by its breakpoint index for the reason the kernel's
transits are (:ref:`chart-and-chord-transits`, "Walls are breakpoint
indices, never positions"): a position in a tuple does not name a wall.
A solid body's one wall is the first entry of any per-wall tuple and
sits at breakpoint :math:`n`, and on a hollow body the first transit's
exit wall and the second transit's entry wall are both breakpoint
:math:`0`. The index is what the kernel's
:class:`~orpheus.geometry.chord.Transits` reports, so a transit's wall
indexes the walls directly.

Two amplitude fields are kept, rather than one amplitude and a tag naming
its kind, because a wall that returns part of its outflow specularly and
part diffusely is physically legitimate and the resolvent serves it
without change: the line part reads ``specular``, the diffuse part reads
``diffuse``. No shipped law declares such a mixture (a geometry admits a
``BC`` tag or a ``BoundaryTraceLaw`` per
boundary point, and the law sums are not laws), so it cannot be declared
today; it is not illegal.

:class:`~orpheus.derivations.continuous.characteristic.walls.Walls` holds
the walls of one body, inner first, with the body's region count and its
:class:`~orpheus.geometry.chart.Chart`. Its three lookups,
:meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.specular_at`,
:meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.diffuse_at`
and :meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.partner_at`,
take arrays of breakpoint indices and raise ``IndexError`` on a
breakpoint that is not a wall: an interior breakpoint, a solid body's
centre, and the kernel's absent-transit code :math:`n + 1`. Each takes
an optional mask; where the mask is false the lookup returns 0 without
consulting the index, which is how the closure reads an absent
traversal.

Both values check their own invariants, so a directly constructed value
cannot hold a state :meth:`Walls.of
<orpheus.derivations.continuous.characteristic.walls.Walls.of>` refuses:
a ``Wall`` refuses an amplitude outside :math:`[0, 1]`; a ``Walls``
refuses walls that are not at distinct breakpoints among
:math:`\{0, n\}` (``ValueError``), a wrap on a radial chart and a wrap
whose partner does not wrap back (the refusals below).


.. _characteristic-walls-factors:

A wall is read from its law's two factors
-----------------------------------------

Every boundary law states two factors
(:ref:`bc-factor-roles`, :ref:`bc-factor-quotients`): the **deck**
:math:`G`, the measure-preserving map of the boundary phase space that
says where returned intensity re-enters (``geometry_map``), and the
**response** :math:`R`, what the surface does to the intensity that
reaches it (``response_kernel``). A symmetry statement (a mirror, a
wrap) puts its content in the deck and returns everything,
:math:`R = I`; a surface (vacuum, a polished wall, a white wall) has the
identity deck and puts its content in the response.

:func:`Walls.of <orpheus.derivations.continuous.characteristic.walls.Walls.of>`
zips the geometry's boundary declarations with the breakpoint indices of
its boundary points (inner first: :math:`(0, n)` on a hollow body or a
slab, :math:`(n,)` on a solid one), parses a tag into its typed law
(:ref:`characteristic-walls-registry`), and reads every law through one
table on its factors, the private function ``_wall_of``:

.. list-table:: The factor table (``_wall_of``); :math:`b` is the wall's breakpoint, :math:`b'` the opposite one
   :header-rows: 1
   :widths: 24 24 20 32

   * - Deck (``geometry_map``)
     - Response (``response_kernel``)
     - Wall (specular, diffuse, partner)
     - Shipped laws that state it
   * - ``PairedDeck`` (a wrap)
     - ``ScalarResponse(1)``
     - :math:`(1, 0, b')`
     - ``PeriodicBoundary``, tag ``periodic``
   * - ``SelfPairedDeck``, not the identity (a mirror)
     - ``ScalarResponse(1)``
     - :math:`(1, 0, b)`
     - ``ReflectiveBoundary``, tag ``reflective``
   * - ``SelfPairedDeck``, the identity
     - ``ScalarResponse(0)``
     - :math:`(0, 0, b)`
     - ``VacuumInflow``,
       ``AlbedoBoundary(0)``, ``PrescribedInflow()`` with no source, tag
       ``vacuum``
   * - ``SelfPairedDeck``, the identity
     - ``SpecularReemission(a)``
     - :math:`(a, 0, b)`
     - ``AlbedoBoundary(a, SpecularReturn)``, tag ``partial``
   * - ``SelfPairedDeck``, the identity
     - ``LambertianReemission(a)``
     - :math:`(0, a, b)`
     - ``AlbedoBoundary(a, IsotropicReturn)``,
       ``WhiteBoundary``, tag ``white``
   * - anything else
     - anything else
     - refused
     - ``AlbedoBoundary(a)`` with no re-emission shape and
       :math:`a > 0`; ``ZeroFluxBoundary``
       (diffusion's law, response :math:`-1`); a mirror or a wrap with a
       response other than ``ScalarResponse(1)``

`[M]` 2026-10-06, each shipped law constructed and passed to ``_wall_of``
at breakpoint 3 of a three-region body: the rows above are what it
returns, and ``AlbedoBoundary(0.37)`` and ``ZeroFluxBoundary()`` raise
``NotImplementedError``.

**Why the walls are read from the factors.** Three reasons, each a
property the class-name alternative lacks.

1. *One reader.* A tag, a law built by hand and a law a future campaign
   adds are read by what they state, not by what their class is called.
   A tag first becomes its typed law and then meets the same table, so a
   tag and the law it names cannot read differently; the gate
   ``test_a_tag_reads_the_same_walls_as_the_law_it_names`` asserts
   ``Walls.of`` of a tagged slab equal to ``Walls.of`` of the slab with
   the typed twin, per admitted kind.
2. *The table is the physics.* A row is a sentence about the boundary: a
   wrap returns everything at the opposite wall, a mirror everything at
   itself, a surface its response at itself. A cell no row covers is a
   law whose return this reference cannot state as a line return and a
   diffuse return, and it is refused by name rather than approximated.
3. *The quotient decks return everything.* A mirror or a wrap is a
   symmetry statement, and a symmetry adds no physics (:math:`R = I`,
   :ref:`bc-factor-quotients`). A partial return is a surface's, spelled
   ``AlbedoBoundary(a, SpecularReturn)`` on the identity deck. So the
   table admits ``ScalarResponse(1)`` only under a quotient deck and a
   re-emission shape only on the identity deck; a mirror carrying
   ``SpecularReemission`` would reflect twice, and a mirror carrying
   ``ScalarResponse(0.5)`` would state a partial wall in the slot that
   says the body is symmetric. The gate
   ``test_a_quotient_deck_returns_everything_and_a_shape_sits_on_the_identity``
   refuses the six cells no shipped law reaches.

A ``PrescribedInflow`` is read through its factors like every law (the
user's ruling "By its factors", 2026-10-06): with
``NoSource`` its factors are vacuum's,
and it is a vacuum wall. Every law reports a ``source``, ``NoSource`` by
default; ``_wall_of`` refuses any other, before the table
(:ref:`characteristic-walls-refusals`).


.. _characteristic-walls-registry:

The reference's own tag registry
--------------------------------

A ``BC`` tag is a kind and a mapping of parameters, ``BC(kind, params)``.
Its meaning is resolved by whoever reads it. The reference resolves it
with
:data:`~orpheus.derivations.continuous.characteristic.walls.TAG_REGISTRY`,
a mapping from each admitted kind, with the wall's context, to the typed
law the kind names. The context is the axis ``"x"`` (the radial or the
slab coordinate) and the wall's outward sign, :math:`-1` at breakpoint
:math:`0` and :math:`+1` at breakpoint :math:`n`:

.. list-table:: The admitted tag kinds
   :header-rows: 1
   :widths: 16 22 62

   * - Kind
     - Parameters (exactly)
     - Typed law
   * - ``vacuum``
     - none
     - ``VacuumInflow()``
   * - ``reflective``
     - none
     - ``ReflectiveBoundary(axis="x")``
   * - ``partial``
     - ``albedo``
     - ``AlbedoBoundary(albedo, SpecularReturn(axis="x"))``, a polished
       wall returning the fraction ``albedo`` specularly
   * - ``white``
     - none
     - ``WhiteBoundary(axis="x", outward_sign=±1)``, a full Lambertian
       return; a partial white wall is spelled as the typed law
       ``WhiteBoundary(albedo=a)``
   * - ``periodic``
     - none
     - ``PeriodicBoundary(axis="x")``

**Strict parameters.** Each kind takes exactly the parameters its row
names: ``BC("partial")`` with no albedo, ``BC("vacuum", {"albedo":
0.5})``, ``BC("white", {"albedo": 0.6})``, ``BC("periodic", {"shift":
0.0})`` and a ``partial`` tag with an extra parameter are each refused
with the fragment *takes exactly the parameters*. A parameter the reader
does not use is a parameter whose author meant something by it; dropping
it silently would answer a different problem from the one declared.
Production's parse drops the white tag's albedo (``BC("white",
{"albedo": a})`` parses as albedo 1), filed as `#583
<https://github.com/deOliveira-R/ORPHEUS/issues/583>`_; the reference
refuses the same declaration. A kind outside the five
(``albedo``, ``marshak``, …) is refused naming the admitted kinds, read
from the registry.

**The context is carried and not read.** No reading of a wall observes
the axis or the outward sign: a wall is amplitudes and a partner. The
registry still hands the white law its wall's outward sign, so the law
it builds is the law a reader declaring that wall would write; the gate
``test_the_registry_hands_the_white_law_the_walls_outward_sign`` pins it
by content equality, the only instrument that can see it (`[M]` battery
arm W13, the sign negated: 2 of 175 rows red, both that gate's).

**Why the registry is the reference's own.** The step's design
(the plan, "P1 step (b), first rung: API sketch", Q1) recommended moving
the body of production's tag parse,
``orpheus.transport.method._law_from_tag``, down to
``orpheus.geometry.boundary`` as one function shared by the
reference and the transport layer, so that a tag's meaning would be
spelled once instead of three times. The user refused it, on
2026-10-06:

   "mixing tags from reference with production would cause reference
   churn during production development. but we do need a coherent way to
   parse tags for reference"

and ruled the reference's own registry. The argument is the insulation
principle (:ref:`architecture-reference-insulation`): a closed reference
is valuable because it stays put, and production's parse is versatile by
design. It consults each method's admission table, it reads the face's
axis and sign from the face label of a multi-dimensional mesh, and it
admits the ``albedo`` kind as an albedo with no re-emission shape. A
change there for production's sake would change what a reference tag
means. The two parses are therefore a **declared duplicate across the
branch line**, recorded in the conceptual view
(:doc:`/architecture/conceptual_view`, the row "Boundary tag to typed
law"), twins that are never merged. They are not two spellings of one
contract: they differ on purpose, the reference refusing what production
admits or drops.


.. _characteristic-walls-refusals:

The refusals are scope boundaries
---------------------------------

Every declaration the reference does not serve is refused with
``NotImplementedError``, a message beginning *the characteristic
reference does not serve a* and carrying one fragment that names why.
The fragments are disjoint over the refused inputs, so one refusal
cannot pass for another; the gate
``test_an_unserved_wall_is_refused_by_its_own_fragment`` asserts each
input's fragment and the absence of every other (18 inputs). The
source and shape arm of ``_wall_of`` carries the tag
``SCOPE-BOUNDARY[guard]``: it is the declared edge of machinery not
built, not a defect.

.. list-table::
   :header-rows: 1
   :widths: 32 24 44

   * - Declaration
     - Fragment
     - Why it is outside the reference
   * - a law with a boundary source
       (``PrescribedInflow(ConstantInflowSource(1.0))``)
     - *with an inflow source*
     - the reference's questions carry no boundary source; the arm is
       served when the reference vocabulary gains a boundary-source
       question (the tag's revisit clause)
   * - a law whose factors no row reads: an albedo with no re-emission
       shape, ``ZeroFluxBoundary``, a quotient deck with a partial or a
       shaped response
     - *unstated re-emission shape*
     - the law says how much returns but not along which lines, or is no
       transport wall
   * - a wrap whose partner does not wrap back (a periodic face opposite
       a vacuum or a mirror)
     - *whose partner does not wrap back*
     - a wrap identifies two walls, so both carry it: the deck is an
       involution
   * - a wrap on a cylinder or a sphere, solid or hollow
     - *periodic wrap on a cylinder or sphere*
     - a translation maps a radial chart's level sets onto none of them;
       ``Walls`` asks the chart
       (:attr:`Chart.acts_on_kept_space <orpheus.geometry.chart.Chart.acts_on_kept_space>`)
   * - a tag kind outside the registry
     - *the reference admits the tag kinds*
     - the message lists the admitted kinds, read from the registry
   * - a tag with a missing or an undeclared parameter
     - *takes exactly the parameters*
     - :ref:`characteristic-walls-registry`, strict parameters
   * - an amplitude outside :math:`[0, 1]`
       (``AlbedoBoundary(1.5, SpecularReturn)``, a white wall of albedo
       :math:`-0.2`)
     - *not a physical wall*
     - a wall returns a fraction of what reaches it

A line lying in an interface is refused by the reference (ruled
2026-10-06); that refusal is not a wall's, and it belongs to the
reading at a point, which the package does not compute
(:ref:`characteristic-what-is-not-built`).


.. _characteristic-period:

The period of the unfolded line
===============================

Follow a line through its walls. Each time it exits the domain at a wall,
the wall returns it (on a specular wall or a wrap) as a congruent line,
which runs through the domain again and exits at a wall. Unfolding the
returns lays these runs end to end into one path. On a concentric body
that path is periodic: after one or two runs it repeats the run it began
with. The closure needs the runs of one period in order, with the wall
each one exits at; :meth:`LinePeriod.of
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.of>`
computes them for a batch of lines, shaped like the chord.

Traversals
----------

The kernel's :attr:`Chord.transits <orpheus.geometry.chord.Chord.transits>`
gives each line at most two transits, the maximal runs of the line inside
the domain, each with its entry and exit wall as breakpoint indices
(:eq:`geometry-transits`). A **traversal** is one transit read forward or
reversed. The four candidates of a line are numbered
:math:`d \in \{0, 1, 2, 3\}`: candidate :math:`d` reads transit
:math:`d \bmod 2`, reversed when :math:`d \ge 2`. A forward candidate
enters at its transit's entry wall and exits at its exit wall; a reversed
one swaps the two. Write :math:`e(d)` and :math:`x(d)` for the entry and
exit wall of candidate :math:`d`, and :math:`\pi(w)` for the partner of
wall :math:`w`.

The successor rule and the rank
-------------------------------

The traversal after candidate :math:`d` is the present candidate that
enters at the partner of the wall :math:`d` exits at, the lowest-numbered
one (a forward candidate before a reversed one). The period starts at
candidate 0, transit 0 read forward, and closes when it returns there:

.. math::
   :label: characteristic-transit-rank

   \sigma(d) \;=\; \min\,\{\, d' : d' \text{ present},\;
   e(d') = \pi\bigl(x(d)\bigr) \,\}, \qquad
   m \;=\; \min\,\{\, j \ge 1 : \sigma^{j}(0) = 0 \,\} \;\in\; \{1, 2\},

with :math:`m = 0` on a line that has no transit. The period is the
sequence :math:`0, \sigma(0), \dots, \sigma^{m-1}(0)`, and :math:`m` is
the line's **rank**, the size of the block of :math:`T` on that line
(:attr:`LinePeriod.rank
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.rank>`).

.. implements:: characteristic-transit-rank
   :by: orpheus.derivations.continuous.characteristic.closure.LinePeriod.of

   **Implemented by** ``LinePeriod.of``, which chains the successor over
   the four candidates of each line (``argmax`` over the matching
   candidates takes the first, a forward one before a reversed one) and
   derives the rank from whether the second traversal is the first.

.. implements:: characteristic-transit-rank
   :by: orpheus.derivations.continuous.characteristic.closure.LinePeriod.rank

.. implements:: characteristic-transit-rank
   :by: orpheus.derivations.continuous.characteristic.walls.Walls.partner_at

The rule names no chart, no boundary kind and no comparison of
:math:`b` with :math:`r_0`; every rank below follows from it, and no case
is tagged.

The ranks, derived
------------------

Apply the rule to each kind of line. A transit is written
:math:`(\text{entry} \to \text{exit})` in breakpoint indices, and
:math:`d\mathrm{f}`, :math:`d\mathrm{r}` are transit :math:`d` forward and
reversed.

- **A solid sphere or cylinder**, any line that meets the body. One
  transit, :math:`(n \to n)`: the centre is not a surface, so the line
  enters and leaves through the outer wall. Its exit :math:`n` has the
  partner :math:`n`; both :math:`0\mathrm{f}` and :math:`0\mathrm{r}`
  enter at :math:`n`, the forward one is taken, and it is the first
  traversal: :math:`m = 1`.
- **A shell, a line that misses the cavity** (:math:`b \ge r_0`). The
  same: one transit :math:`(n \to n)`, :math:`m = 1`. At :math:`b = r_0`
  exactly the line touches the inner surface, and a tangency is not a
  crossing (:ref:`chart-and-chord-tangency`): the cavity slot has length
  0, it does not break the run, and the line has one transit.
- **A shell, a line through the cavity** (:math:`b < r_0`). Two
  transits, :math:`(n \to 0)` and :math:`(0 \to n)`. Transit 0 exits at
  the inner wall, partner 0; the candidates entering at 0 are
  :math:`1\mathrm{f}` and :math:`0\mathrm{r}`, and :math:`1\mathrm{f}` is
  taken. It exits at :math:`n`, partner :math:`n`; the candidates
  entering there are :math:`0\mathrm{f}` and :math:`1\mathrm{r}`, and
  :math:`0\mathrm{f}` is the first traversal: :math:`m = 2`. One ulp
  inside :math:`r_0` (:math:`b` = ``nextafter(r_0, 0)``) the cavity slot is
  traversed and the rank is 2.
- **A slab between two returning walls** (mirrors, polished walls,
  vacuum as amplitude 0). One transit, :math:`(0 \to n)` rising or
  :math:`(n \to 0)` falling. Rising: :math:`0\mathrm{f}` exits at
  :math:`n`, partner :math:`n` (a mirror re-enters where it left); only
  :math:`0\mathrm{r}` enters at :math:`n`; it exits at 0, partner 0, where
  :math:`0\mathrm{f}` enters: :math:`m = 2`, the transit and its
  reverse.
- **A periodic slab.** One transit, :math:`(0 \to n)` rising. Its exit
  :math:`n` has the partner :math:`0` (the wrap re-enters at the opposite
  face); :math:`0\mathrm{f}` enters at 0 and is the first traversal:
  :math:`m = 1`, and the unfolded line continues in the same direction,
  with no reversal. Falling is the mirror image.
- **A line that meets no wall**: the kernel's parallel line (a cylinder
  line along the axis, a slab line with :math:`\Omega_x = 0`) and a line
  that misses or grazes the body (:math:`b \ge r_n`). No transit,
  :math:`m = 0`, every column absent.

A vacuum wall is a wall with amplitude 0, not a missing wall: it has a
partner (itself), the period runs through it, and the closure then
returns exactly nothing across it (:ref:`characteristic-closure-section`).
So a slab with a vacuum face and a mirror face has rank 2, with the
amplitudes 1 and 0 on its two traversals.

Why the forward candidate, and why that is not a choice
-------------------------------------------------------

Two candidates enter at the same wall only on the radial charts: on the
solid body :math:`0\mathrm{f}` and :math:`0\mathrm{r}`; on a shell's
line through the cavity :math:`1\mathrm{f}` and :math:`0\mathrm{r}` (and
:math:`0\mathrm{f}` and :math:`1\mathrm{r}`). On the slab the orbit
coordinate is monotone along a line, so at a mirror only the reversed
transit enters and across a wrap only the forward one.

The two candidates of a radial chart trace the same path in the orbit
space. Along a line the orbit coordinate is
:math:`c(t)^2 = b^2 + \bigl(|P\Omega|\,(t - t^*)\bigr)^2`
(:eq:`geometry-line-crossing-law`), even in :math:`t - t^*`: the
reflection of the line through its closest approach is a motion that
keeps the impact parameter and the projected speed, maps every level
set of the chart to itself, and carries one candidate onto the other.
So the two cross the same regions with the same lengths in the same
order, have the same optical depth, and, for an emission density that is
a function of the orbit coordinate (as every density on a 1-D body is),
the same source integral. Which one the period names
changes the bookkeeping, not the closure.

The forward-first convention makes the period the shortest one. Taking
the reversed candidate on a solid body gives the period
:math:`(0\mathrm{f}, 0\mathrm{r})`, the rank-1 period traversed twice;
its cycle has the same least solution (`[M]` 2026-10-06, a rank-1 period
with :math:`a = 0.6`, :math:`\tau = 0.83`, :math:`B = 1.7` and the
doubled rank-2 period with both traversals equal: the inflows agree
bitwise, 1.381420). With the convention the rank is the minimal period,
the number the closure's block size and the gates count. The battery's
arm C4 prefers the reversed candidate and reddens 92 of the 175 rows,
through the rank it changes.

Why a period closes within two traversals
-----------------------------------------

A line has at most two transits (:ref:`chart-and-chord-transits`, "Why
at most two"). A slab line has one, which the period can read in both
directions; a radial line's successor is always one of its own transits
read forward (the derivation above). Four candidates and that structure
leave no period longer than two, which is the module's constant
``_MAX_PERIOD = 2``. Every array of a ``LinePeriod`` is ``(..., 2)``, one
column per traversal, with ``present`` masking the line's :math:`m`
columns.

``LinePeriod.of`` checks the structure instead of trusting it. If a line
has a transit and no candidate enters at the partner of the wall it left,
it raises ``RuntimeError`` (*enters at the partner*); that is reachable
from valid inputs by reading a radial chord through a periodic slab's
walls, a chord and walls from two different bodies. If a rank-2 period
does not return to its first traversal, it raises ``RuntimeError``
(*did not close*); a valid ``Walls`` makes the partner map an involution,
so that guard is reached only by a deliberately broken partner map, and
the gate ``test_a_chord_whose_walls_disagree_with_its_transits_is_refused``
reaches it that way.

What a period holds
-------------------

``LinePeriod`` stores the chord, ``candidate`` (each traversal as its
candidate number, 0 where absent), ``amplitude`` (the specular amplitude
of the wall the traversal exits at, read with
:meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.specular_at`;
0 where absent) and ``present``. It derives ``transit``
(:math:`d \bmod 2`), ``reversed``, ``entry_wall`` and ``exit_wall`` (the
transits' absent code :math:`n + 1` where absent) and the rank.

Worked cases
------------

`[M]` 2026-10-06, :meth:`LinePeriod.of
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.of>`
at ``72199f9d``. The bodies: the solid sphere with breakpoints
:math:`(0, 0.3, 1.1, 2.0)` and a polished outer wall of albedo 0.6; the
hollow sphere and cylinder :math:`(0.4, 1.1, 2.0)` with polished walls
of albedo 0.3 inside and 0.6 outside; the slab
:math:`(-0.7, 0.3, 1.1, 2.0)` with the same two walls, or periodic. Each
radial line passes through :math:`(b, 0, 0)` along :math:`\hat e_y`,
each slab line through :math:`(0.5, 0, 0)`. A traversal is written
(transit, reversed, entry wall, exit wall, amplitude).

.. list-table::
   :header-rows: 1
   :widths: 34 8 58

   * - Line
     - Rank
     - Period
   * - solid sphere, :math:`b = 0`, :math:`0.7`, :math:`1.1`
     - 1
     - (0, forward, 3, 3, 0.6); :math:`b = 1.1` is an interior tangency,
       which changes no wall
   * - solid sphere, :math:`b = 2.0`
     - 0
     - none: a tangency to the outer surface is not a crossing
   * - hollow sphere, :math:`b = 0.2` and :math:`b = 0.4 - 1\,\mathrm{ulp}`
     - 2
     - (0, forward, 2, 0, 0.3), (1, forward, 0, 2, 0.6)
   * - hollow sphere, :math:`b = 0.4` and :math:`0.7`
     - 1
     - (0, forward, 2, 2, 0.6)
   * - hollow cylinder, :math:`b = 0.2`, direction :math:`(0, 0.6, 0.8)`
     - 2
     - as the hollow sphere: the in-plane impact parameter decides
   * - slab between the two walls, :math:`\Omega_x = 0.6`
     - 2
     - (0, forward, 0, 3, 0.6), (0, reversed, 3, 0, 0.3)
   * - the same, :math:`\Omega_x = -0.6`
     - 2
     - (0, forward, 3, 0, 0.3), (0, reversed, 0, 3, 0.6)
   * - periodic slab, :math:`\Omega_x = 0.6` and :math:`-0.6`
     - 1
     - (0, forward, 0, 3, 1.0) and (0, forward, 3, 0, 1.0)
   * - slab, vacuum left, mirror right, :math:`\Omega_x = 0.6`
     - 2
     - (0, forward, 0, 3, 1.0), (0, reversed, 3, 0, 0.0)
   * - solid cylinder along :math:`\hat e_z` at :math:`c = 0.7`; slab
       along :math:`(0, 0.6, 0.8)`
     - 0
     - none: a parallel line meets no wall

The gate ``test_the_period_matches_the_hand_counted_table`` holds 73 such
lines with the walls counted by hand, on the sphere, the cylinder at
:math:`\Omega_z \in \{0, 0.8, -0.6\}` and the slab;
``test_the_period_is_the_physically_unfolded_path`` builds the unfolded
path without the period at all, by reflecting the 3-D line at each wall
with a Householder map the test constructs from the wall normal (or
translating it across a wrap), and asserts that each run is the
traversal the period names, region by region and length by length
(:ref:`characteristic-evidence`).


.. _characteristic-closure-section:

The line closure
================

The cycle
---------

Take one line and its period of :math:`m` traversals. Along traversal
:math:`k`, of optical depth :math:`\tau_k` and arc length :math:`L_k`,
the transport equation :math:`\mathrm{d}\psi/\mathrm{d}s + \Sigma_t\psi = q`
integrates from the entry (:math:`s = 0`) to the exit
(:math:`s = L_k`) as

.. math::

   \psi_k(L_k) \;=\; e^{-\tau_k}\,\psi^{\rm in}_k \;+\; B_k, \qquad
   B_k \;=\; \int_0^{L_k} q(s)\,e^{-\tau_k(s, L_k)}\,\mathrm{d}s ,

where :math:`\psi^{\rm in}_k` is the intensity entering the traversal,
:math:`\tau_k(s, L_k)` the optical depth from :math:`s` to the exit and
:math:`B_k` the traversal's **outflow**, its source integral attenuated
to its exit. The wall traversal :math:`k` exits at returns the fraction
:math:`a_k` of that intensity along the next traversal of the period, so

.. math::
   :label: characteristic-closure

   \psi^{\rm in}_{k+1} \;=\; a_k\,\bigl(e^{-\tau_k}\,\psi^{\rm in}_k + B_k\bigr)
   \quad (k \bmod m), \qquad
   \psi^{\rm in}_k \;=\;
   \frac{r_{k-1} + \gamma_{k-1}\,r_{k-2}}{1 - \gamma_0\,\gamma_1}
   \quad (k \bmod 2),

with :math:`r_j = a_j B_j` what traversal :math:`j` returns,
:math:`\gamma_j = a_j e^{-\tau_j}` its gain, and an absent traversal
counted as the cycle's unit, :math:`r = 0` and :math:`\gamma = 1`. The
second form is the least non-negative solution of the first; the
derivation follows.

.. implements:: characteristic-closure
   :by: orpheus.derivations.continuous.characteristic.closure.LinePeriod.inflow

   **Implemented by** ``LinePeriod.inflow``, which takes the optical
   depths :math:`\tau_k` and the outflows :math:`B_k` (with any trailing
   axes, a basis for instance) and returns the inflow at each traversal's
   entry, 0 where absent.

.. implements:: characteristic-closure
   :by: orpheus.derivations.continuous.characteristic.closure.LinePeriod.optical_depth

.. implements:: characteristic-closure
   :by: orpheus.derivations.continuous.characteristic.walls.Walls

The amplitude :math:`a_k` is the wall's specular amplitude only. A
diffuse amplitude returns nothing along this line, and its return is the
second part of the resolvent (:ref:`characteristic-resolvent`).

The albedo pairing
------------------

Which wall's amplitude multiplies which inflow is the step most easily
written wrong, because "the albedo of the wall the transit ends on" names
a different wall depending on whether the transit is read forward or
backward. The rule is the cycle's indexing: **the inflow to traversal**
:math:`k + 1` **carries the amplitude of the wall traversal** :math:`k`
**exits at**. Read backward, that is the wall at which the backward path
from traversal :math:`k + 1` reflects. ``LinePeriod`` stores the
amplitude of each traversal's exit wall and ``inflow`` rolls it one
place forward, so the pairing is the code's indexing and not a
convention a caller supplies.

The pairing is invisible wherever the two walls have the same amplitude,
which is why the gates draw distinct amplitudes (0.3 inside, 0.6
outside, and eight amplitude pairs in the closure row). Its size where
they differ was measured before the code existed, on D2's closed form
for :math:`\psi` (the design review's probe
``scratch/characteristic_architecture/p1_probes/probe_d2_sum.py``,
`[M]` 2026-10-06, re-run here): on the fixture first leg 0.3 of optical
depth 0.7, :math:`B = (0.9, 0.1)`, :math:`\tau = (0.5, 0.5)` and amplitudes
:math:`(0.9, 0.2)`, the correct pairing gives :math:`\psi = 0.7366` and the
swapped pairing :math:`0.4015`, a difference of 0.335, 45 % of
:math:`\psi`. The same probe confirms the correct pairing against the
explicitly unfolded path to :math:`1.4 \times 10^{-15}` relative over
2000 seeded draws of ranks 1 and 2.

The closed form, derived
------------------------

Write the cycle as a fixed point on the :math:`m` inflows,
:math:`\psi^{\rm in} = G\,\psi^{\rm in} + s`, with
:math:`(G\psi)_{k+1} = \gamma_k\,\psi_k` and :math:`s_{k+1} = r_k`.
:math:`G` is a weighted cyclic shift with non-negative weights, and
:math:`G^m = \Pi\,I` with the **cycle product**

.. math::

   \Pi \;=\; \prod_{k=0}^{m-1} a_k\,e^{-\tau_k} \;\in\; [0, 1].

The least non-negative solution of the fixed point is the Neumann series
:math:`\sum_{j \ge 0} G^j s`, the sum over the unfolded path of every
return, each attenuated by the walls and the traversals between it and
the inflow it reaches. Grouping the terms by whole periods,

.. math::

   \sum_{j \ge 0} G^j s
   \;=\; \Bigl(\sum_{p \ge 0} \Pi^p\Bigr) \sum_{j=0}^{m-1} G^j s
   \;=\; \frac{1}{1 - \Pi}\sum_{j=0}^{m-1} G^j s
   \qquad (\Pi < 1).

For :math:`m = 1` the inner sum is :math:`s_0 = r_0`, so
:math:`\psi^{\rm in}_0 = a_0 B_0 / (1 - a_0 e^{-\tau_0})`, the rank-1
resolvent :math:`(1 - \alpha e^{-\tau})^{-1}` of the trajectory-resolvent
family, which that family's page identifies with Sanchez 1986 Appendix
Eq. (A4) :cite:`SanchezTTSP1986` and Pomraning–Siewert 1982 Eq. (14)
:cite:`PomraningSiewert1982`. For :math:`m = 2` it is
:math:`s_k + \gamma_{k-1} s_{k-1}`:

.. math::

   \psi^{\rm in}_0 = \frac{a_1 B_1 + \gamma_1\,a_0 B_0}{1 - \gamma_0\gamma_1},
   \qquad
   \psi^{\rm in}_1 = \frac{a_0 B_0 + \gamma_0\,a_1 B_1}{1 - \gamma_0\gamma_1},

the inflow to :math:`k` being what traversal :math:`k - 1` returns plus
what traversal :math:`k - 2` returns carried once through
:math:`k - 1`. Counting an absent traversal as the unit
(:math:`r = 0`, :math:`\gamma = 1`) makes the :math:`m = 1` case the
:math:`m = 2` expression with the second traversal absent: the column of
the absent traversal returns nothing, and its gain of 1 carries the
present traversal's return around to itself. That is the second form of
:eq:`characteristic-closure`, one expression for every rank, which
``inflow`` evaluates with two ``np.roll`` calls along the period axis and
no branch on the rank. Rank 0 has no present column and returns zeros.

Two exact edges follow, and the gates assert them bitwise:

- **Every amplitude 0** (vacuum everywhere): every :math:`r_k` is 0, so
  the inflow is exactly 0 and finite. The line closure adds nothing to
  the vacuum operator (``test_every_amplitude_zero_adds_nothing``).
- **One wall absorbing**: the traversal whose inflow that wall feeds has
  inflow exactly 0, and the cycle is cut, :math:`\Pi = 0`, so the other
  traversal's inflow is exactly its predecessor's return
  (``test_an_absorbing_wall_zeroes_exactly_the_inflow_it_feeds``). This
  is the edge where the pairing shows bitwise: the swapped pairing puts
  the mirror's 1 on the cut.

The least solution, and the trapped line
----------------------------------------

Every amplitude is at most 1 and every optical depth is non-negative,
so :math:`\Pi = 1` exactly when every amplitude of the period is 1 and
every optical depth is 0: a **lossless trapped line**, which never meets
material and never loses anything at a wall. Such lines exist as soon as
voids are admitted (ruled 2026-10-06): a line of a sphere with a void
outer shell and a mirror outside, whose impact parameter puts it in the
shell only; a slab of void regions between two mirrors; a void periodic
slab. There :math:`I - G` is singular and the closed form is
:math:`0/0` or :math:`B/0`.

The Neumann series is still defined, and it is the answer. With no
source on the line (:math:`s = 0`) every term is 0 and the least
solution is exactly 0, the :math:`\alpha \to 1^-` limit of the
non-singular case and the physically right answer: nothing is emitted on
the line and nothing enters it. With a source (:math:`s \ne 0`) the
terms recur with period :math:`m` (:math:`G^m = I`) and the series
diverges: there is no finite answer, and ``inflow`` raises
:class:`~orpheus.derivations.continuous.characteristic.closure.TrappedSource`.
That is the closure's only refusal; it is a ``ValueError``, raised on
the outflow the caller supplied, and it names the trapped line.

The ruling that chose this form (the plan's ledger, "the user, second
batch", 2026-10-06) settled the question of whether a less specialised
closure would avoid the singularity: it would not on a specular wall,
because the general resolvent has the same singular :math:`I - T` on a
line whose period is all void with every amplitude 1; the singularity is
the problem's. A white wall couples every line to the material, so it
traps nothing.

In code: ``inflow`` forms :math:`1 - \Pi`, marks the trapped lines
(:math:`1 - \Pi = 0`), raises if any trapped line has a non-zero
numerator, and divides only where the line is not trapped, writing 0
elsewhere. The refusal is the trapped line's and not a void's: a rank-2
line lossless on one traversal and lossy on the other has
:math:`\Pi < 1` and a finite, positive inflow
(``test_the_least_solution_on_a_lossless_trapped_line``, which batches a
trapped line, a line through material, and that mixed line).

The cancellation-free :math:`1 - \Pi`
-------------------------------------

On a nearly lossless line, :math:`\Pi` is within a few
:math:`10^{-12}` of 1, and :math:`1 - \Pi` computed as ``1 - Pi`` keeps
only the digits of :math:`\Pi` that differ from 1. ``inflow`` takes the
optical depths, not the transmissions :math:`e^{-\tau_k}`, so that it can
form

.. math::

   1 - \Pi \;=\; -\operatorname{expm1}\Bigl(\sum_k \bigl(\log a_k - \tau_k\bigr)\Bigr)

with every term of the sum small and exact to rounding. `[M]` 2026-10-06,
a rank-2 line between two mirrors with :math:`\tau = (10^{-12},
2.5 \times 10^{-12})` against the closed form in mpmath at 50 digits:
:math:`1 - \Pi` by ``expm1`` carries a relative error of
:math:`3.1 \times 10^{-17}`, by ``1 - exp`` :math:`6.3 \times 10^{-6}`,
and the two inflows by ``expm1`` are within :math:`1.1 \times 10^{-16}`
relative. ``test_a_nearly_lossless_line_keeps_its_digits`` asserts 4 ulp
on three such lines; the battery's arm C8 (``1 - exp``) reddens all
three. An amplitude of 0 gives :math:`\log 0 = -\infty`, a gain of
exactly 0 and :math:`1 - \Pi = 1`, so vacuum walls need no branch; the
divide-by-zero warning of the logarithm is suppressed at that line only.

The optical depth of a traversal
--------------------------------

The closure takes :math:`\tau_k` as data, and
:meth:`LinePeriod.optical_depth
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.optical_depth>`
computes it from the total cross section of each region:

.. math::

   \tau_k \;=\; \sum_{i \in \text{slots of } k} \Sigma_t(\rho_i)\,\ell_i ,

over the slots of the traversal's transit, with :math:`\rho_i` and
:math:`\ell_i` the chord's slot region and slot length
(:ref:`chart-and-chord-slots`). A transit and its reverse cross the same
slots, so they have the same depth. Three details make it exact where it
could be wrong:

- **The exteriors are void.** The codes :math:`n` and :math:`n + 1` read
  :math:`\Sigma_t = 0`. A transit's traversed slots are interior by
  definition (:eq:`geometry-transits`); an exterior slot inside a
  transit's slice is untraversed, of length 0 (the cavity slot of a
  shell's line at :math:`b = r_0`). So the value given to an exterior
  cannot reach the sum: the battery's arm C7b gives the exteriors
  :math:`\Sigma_t = 5` and is declared null, green on every row.
- **A void slot adds exactly 0, even when it is infinitely long.** The
  product :math:`\Sigma_t\,\ell` is formed only where :math:`\Sigma_t > 0`,
  so a void slot of infinite length (a slab line with a subnormal
  :math:`\Omega_x`) adds 0 rather than :math:`0 \times \infty`. `[M]`
  2026-10-06, ``test_a_subnormal_cosine_through_a_void_reads_infinite_depth_and_finite_inflow``:
  :math:`\Omega_x = 5 \times 10^{-324}` through the slab
  :math:`\Sigma_t = (0.5, 0, 0.9)` gives :math:`\tau = +\infty` on both
  traversals, never a NaN, and the inflow is the one-bounce return
  :math:`a_{k-1} B_{k-1}`, finite; with every region void the depth is 0
  and the line is a trapped one.
- **The cross sections are one per region.** A ``sigma_t`` of any other
  length is refused (``ValueError``, *one total cross section per
  region*), never broadcast or truncated.

The gates compare the depths with closed-form segment lengths computed
in mpmath (9 radial lines, a void shell among them, and the slab's
widths over the cosine), to 8 ulp relative. The depth is not a labelled
equation of this page and its rows remain software-invariant rows; it
is part of what ``optical_depth`` implements under
:eq:`characteristic-closure`.


.. _characteristic-panel-basis:

The panel basis
===============

The closure takes the outflows :math:`B_k` as data, and an outflow is the
emission density integrated along a traversal. The density, and the flux
the reference solves for, are represented on a basis over the orbit
coordinate :math:`c` of the body (the radius on a sphere or a cylinder,
the depth on a slab; :eq:`geometry-radial-coordinate`). The basis is
:class:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis`,
built by :meth:`PanelBasis.of
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.of>`
from the body's partition and three numbers: the degree :math:`p`, the
number of layers :math:`L` and the ratio :math:`\rho \in (0, 1)`. It reads
nothing else. Two bodies that differ only in their materials or their
boundary laws have one basis, bitwise, and a moved breakpoint changes it
(``test_the_basis_reads_only_the_bodys_partition_and_the_resolution``).

Panels graded toward walls and interfaces
-----------------------------------------

Each region :math:`[r_k, r_{k+1}]` is cut into **panels**. An end of the
region is **graded** when it is a wall or an interface, and it is not
graded when it is a singular stratum of the chart: the centre of a solid
sphere or the axis of a solid cylinder (:ref:`chart-and-chord-strata`).
Write :math:`w` for the half width of the region when both its ends are
graded and for its whole width when one is. The interior panel ends lie
at the depths

.. math::

   w\,\rho^{j}, \qquad j = 1, \dots, L,

from each graded end: at :math:`r_k + w\rho^j` toward :math:`r_k` and at
:math:`r_{k+1} - w\rho^j` toward :math:`r_{k+1}`. With :math:`L = 0` a
region is one panel. The private function ``_graded_ends`` computes them;
the gate ``test_panels_grade_toward_walls_and_interfaces_and_not_toward_a_singular_stratum``
writes the same law by hand in mpmath and compares at 4 ulp of the outer
radius, for :math:`L \in \{0, 1, 2, 4\}` and :math:`\rho \in \{0.3, 0.5\}`.

**The law places depths, not widths.** Counting from a graded end, the
first panel lies between the depths 0 and :math:`w\rho^{L}`, and the panel
between the depths :math:`w\rho^{j}` and :math:`w\rho^{j-1}` has width
:math:`w\rho^{j-1}(1 - \rho)`. The widths shrink by :math:`\rho` from one
panel to the next except at the end itself: with :math:`\rho = 1/2` the
two panels at the end have the same width :math:`w/2^{L}`. When both ends
are graded, one middle panel :math:`[r_k + w\rho, r_{k+1} - w\rho]` of
width :math:`2w(1 - \rho)` joins the two gradings; when one end is graded,
the last panel, between the depths :math:`w\rho` and :math:`w`, reaches the
ungraded end.

`[M]` 2026-10-06, :meth:`PanelBasis.of
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.of>` at
:math:`L = 2`, :math:`\rho = 1/2`:

.. list-table::
   :header-rows: 1
   :widths: 34 66

   * - Body (breakpoints)
     - Panel ends
   * - solid sphere :math:`(0, 0.5, 1.5, 2.0)`
     - 0, 0.25, 0.375, **0.5**, 0.625, 0.75, 1.25, 1.375, **1.5**,
       1.5625, 1.625, 1.875, 1.9375, **2.0**: the centre is not graded,
       so the first region is graded toward 0.5 alone, with :math:`w = 0.5`
   * - hollow sphere :math:`(0.4, 0.5, 1.5, 2.0)`
     - **0.4**, 0.4125, 0.425, 0.475, 0.4875, **0.5**, 0.625, …: the
       cavity surface is a wall and is graded
   * - solid sphere :math:`(0, 0.01, 1)`
     - 0, 0.005, 0.0075, **0.01**, 0.13375, 0.2575, 0.7525, 0.87625,
       **1.0**

**Why grade, and why not at the centre.** The flux has its boundary layers
and its derivative singularities at the walls and the interfaces, where
the cross section or the boundary condition jumps, so a polynomial on a
panel reaching such an end converges slowly; geometric panels toward the
end confine the singular behaviour to panels that shrink with it. At a
singular stratum nothing jumps: the chart's symmetry group fixes the
centre, the flux is an even and smooth function of :math:`c` there, and
grading would add panels and nothing else. The gate's arms B1 (the centre
graded) and B2 (the inner wall not graded) redden 12 and 14 rows.

**The panel ends are a partition.** They are a
:class:`~orpheus.geometry.chord.ConcentricPartition` on the body's chart
with the body's pose, a refinement of the body's: it contains every
breakpoint of the body and has the same two ends. A directly constructed
``PanelBasis`` checks it (``ValueError``, each case by its own fragment:
another chart, other ends, a missing breakpoint), so each panel lies in
one region, read from :attr:`PanelBasis.region_of_panel
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.region_of_panel>`,
and a per-region table, the cross sections for instance, becomes a
per-panel one with :meth:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis.on_panels`,
which refuses a table without one entry per region. Being a partition is
what lets the geometric kernel chord a line through the panels
(:ref:`characteristic-transport`).

The functions
-------------

On a panel :math:`[a, b]` the basis has the :math:`p + 1` nodal Lagrange
functions through the panel's Gauss–Legendre points: with
:math:`\xi_0 < \dots < \xi_p` the Gauss–Legendre points of
:math:`[-1, 1]` and :math:`x = (2c - (a + b))/(b - a)` the panel's local
coordinate,

.. math::

   u_{P,m}(c) \;=\; \prod_{l \ne m} \frac{x - \xi_l}{\xi_m - \xi_l}
   \quad (c \in [a, b]), \qquad u_{P,m}(c) = 0 \text{ elsewhere}.

Function :math:`m` of panel :math:`P` is basis function
:math:`i = P(p + 1) + m`, and its node is the :math:`m`-th Gauss–Legendre
point of the panel; the nodes come from the reference kernel's
:func:`~orpheus.derivations.common.quadrature.composite_gauss_legendre`,
in panel order, :math:`N = P_{\rm total}(p + 1)` of them. Three
properties follow.

- **A coefficient is a value.** :math:`u_i(c_j) = \delta_{ij}` at the
  nodes, so a vector of coefficients is the vector of the represented
  function's values at the nodes, and the dense pencil's single-sign test
  on a fundamental mode reads coefficients as values
  (:ref:`verification-reference-kernel`). The product form is evaluated
  as written: it is exactly one-hot at a node and divides only by the
  node spacings, never by :math:`x - \xi_l`.
- **The basis is discontinuous.** Functions of different panels have
  disjoint supports, and no node lies on a panel end, because the
  Gauss–Legendre points are interior. So
  :meth:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis.values`
  takes a panel as well as an orbit coordinate: at a panel end the panel
  names the side. A discontinuous basis is what an emission density with
  a jump at every interface needs.
- **It reproduces each per-region polynomial of degree at most** :math:`p`.
  The interpolant at the nodes of any function that is a polynomial of
  degree :math:`p` on each region (jumping at every interface) is that
  function; one of degree :math:`p + 1` is not reproduced, the gate's
  loading control (``test_the_interpolant_reproduces_every_per_region_polynomial_of_degree_p``).

The mass matrix and the volume density
--------------------------------------

The basis's metric is its Gram matrix in the chart's volume measure,
:math:`W_{ij} = \int u_i u_j \,\mathrm{d}V`
(:attr:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis.mass`).
The chart defines the measure once, by the measure of a cell,

.. math::

   m\bigl([r_j, r_{j+1}]\bigr) \;=\; \kappa\,\bigl(T(r_{j+1}) - T(r_j)\bigr),
   \qquad T(r) = r^{d},

(:meth:`CoordSystem.measure <orpheus.geometry.coord.CoordSystem.measure>`,
read by :meth:`Chart.measure <orpheus.geometry.chart.Chart.measure>`),
with the constant :math:`\kappa` (``measure_constant``) and the exponent
:math:`d` (``measure_coordinate.exponent``) equal to :math:`(1, 1)` on the
slab, :math:`(\pi, 2)` on the cylinder and :math:`(4\pi/3, 3)` on the
sphere (`[M]` 2026-10-06, read from each chart's ``CoordSystem``). Its
density in the orbit coordinate is the derivative of that one definition,

.. math::

   \frac{\mathrm{d}V}{\mathrm{d}r} \;=\; \kappa\,d\,r^{d-1}
   \;\in\; \{\,1,\; 2\pi r,\; 4\pi r^{2}\,\},

(:meth:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis.volume_density`),
so :math:`W_{ij} = \int u_i(r)\,u_j(r)\,\kappa\,d\,r^{d-1}\,\mathrm{d}r`. On
a panel the integrand is a polynomial of degree :math:`2p + d - 1 \le
2p + 2`, and Gauss–Legendre with :math:`p + 2` points is exact to degree
:math:`2p + 3`, so each entry is integrated exactly up to rounding.
:math:`W` is block-diagonal, one :math:`(p + 1) \times (p + 1)` block per
panel, because the supports are disjoint.

**Why the density is derived inside the basis.** The kernel exposes the
cell measures and no density. The user ruled on 2026-10-06 (the plan's
ledger, "on P1 step (b)'s second rung", Q4) that the basis derives the
density from ``measure_constant`` and ``measure_coordinate`` rather than
adding a verb ``Chart.volume_density`` to the stable kernel. The density
is then the derivative of the measure's one definition, never a second
spelling of it: the gate ``test_each_panels_mass_is_its_chart_measure``
sums each panel's block, which with a partition of unity is the panel's
volume, against :meth:`Chart.measure <orpheus.geometry.chart.Chart.measure>`
and against the closed-form volume in mpmath. The first leg shares the
upstream constant :math:`\kappa` with the code and is not independent; the
second leg, and ``test_the_volume_density_is_the_charts_by_hand`` (the
three densities written by hand), are. Arm B7 (the factor :math:`d`
dropped) reddens 44 rows.

The walls on the panel partition
--------------------------------

The period of a line (:ref:`characteristic-period`) is read from a chord
and the walls, and the transport reads the chord through the panels, so
the walls must be keyed on the panel partition.
:meth:`Walls.on <orpheus.derivations.continuous.characteristic.walls.Walls.on>`
re-keys them: a wall is an end of the domain, so breakpoint 0 stays 0 and
breakpoint :math:`n` of the body becomes the panel partition's last index
:math:`P_{\rm total}`, and each partner is re-keyed with it; the
amplitudes are kept. It refuses a partition on another chart. It cannot
refuse a partition that is not a refinement with the body's ends,
because a ``Walls`` holds breakpoint indices and no positions; that
invariant is the panel partition's own, checked by ``PanelBasis``, and
:meth:`TraversalRule.of
<orpheus.derivations.continuous.characteristic.transport.TraversalRule.of>`
calls ``walls.on`` with the basis's own partition only
(``test_the_walls_rekey_onto_the_panel_partition``; arms W1, the last
index kept at :math:`n`, and W2, the wrap's partner not re-keyed, redden
87 rows each).


.. _characteristic-transport:

The transport along a line
==========================

Along one line the transport equation is an ordinary differential
equation in arc length, and on the panel basis every quantity the line
closure and a Galerkin assembly need is a linear functional of the
basis coefficients.
:class:`~orpheus.derivations.continuous.characteristic.transport.TraversalRule`
computes them for a batch of lines.

Building the rule
-----------------

:meth:`TraversalRule.of
<orpheus.derivations.continuous.characteristic.transport.TraversalRule.of>`
takes the lines, the basis, the body's walls, the total cross section of
each **region** and two point counts, ``points`` (Gauss–Legendre points
per piece) and ``inner_points`` (per interval of the attenuated
integral). It chords the lines through the basis's panel partition
(:meth:`ConcentricPartition.chord <orpheus.geometry.chord.ConcentricPartition.chord>`),
reads their period through the walls re-keyed onto it,
``LinePeriod.of(basis.partition.chord(lines), walls.on(basis.partition))``,
and reads the cross sections onto the panels with ``basis.on_panels``. One
chord then serves the period, the optical depths and every source
integral: each slot of the chord lies in one panel, its length is the
kernel's cancellation-free length (:eq:`geometry-chord-segment-lengths`),
its panel is its region code in the panel partition, and the line's
closest approach is already a slot end. There is no second crossing
computation. The optical depth of a traversal on the panel chord equals
its depth on the body's chord (the gate
``test_the_panel_chord_unfolds_into_the_hand_counted_period``, 12 lines, 16
ulp).

The user ruled that the pieces come from the kernel's chord through the
panel partition and that :math:`B_k` lives on this new value
(2026-10-06, Q2 and Q3 of the second rung's sketch). The sketch proposed
the factory ``of(period, basis, sigma_per_panel, resolution)``; the
elegance review (its finding C2) showed that it left two constructions
the types still spelled, a period chorded through some other partition
and walls re-keyed onto a partition that is not the basis's, and the
signature ``of(lines, basis, walls, sigma_t, ...)`` makes both
unspellable through the factory. A directly constructed
``TraversalRule`` still checks that its period's chord is through the
basis's own partition object, that ``sigma`` holds one finite,
non-negative value per panel and that each rule has at least one point.
It also refuses a traversed slot whose 3-D length overflowed (a line
within an underflow of parallel to the level sets), a guard tagged
``ELEGANCE-DEBT[guard]`` that retires when the kernel's chord refuses or
resolves such a line (`#582
<https://github.com/deOliveira-R/ORPHEUS/issues/582>`_).

.. _characteristic-traversal-integrals-section:

The traversal integrals
-----------------------

Take traversal :math:`k` of a line's period (:ref:`characteristic-period`),
with arc length :math:`s \in [0, L_k]` from its entry, orbit coordinate
:math:`c(s)` and total cross section :math:`\Sigma_t(c)`, constant on each
panel. Write :math:`\tau_k(s, s') = \int_s^{s'} \Sigma_t(c(s''))\,\mathrm{d}s''`
for the optical depth between two points of the traversal. Its three
integrals are

.. math::
   :label: characteristic-traversal-integrals

   \tau_k \;=\; \tau_k(0, L_k), \qquad
   B_k[u_i] \;=\; \int_0^{L_k} u_i\bigl(c(s)\bigr)\,e^{-\tau_k(s, L_k)}\,\mathrm{d}s,
   \qquad
   A_k[u_i] \;=\; \int_0^{L_k} u_i\bigl(c(s)\bigr)\,e^{-\tau_k(0, s)}\,\mathrm{d}s
   \;=\; B_{\bar k}[u_i],

with :math:`\bar k` the reversed traversal, the same transit read the
other way.

- :math:`\tau_k` is the traversal's **optical depth**, the sum over its
  slots of :math:`\Sigma_t\,\ell` that :meth:`LinePeriod.optical_depth
  <orpheus.derivations.continuous.characteristic.closure.LinePeriod.optical_depth>`
  forms (:ref:`characteristic-closure-section`);
  :attr:`TraversalRule.optical_depth
  <orpheus.derivations.continuous.characteristic.transport.TraversalRule.optical_depth>`
  is that method on the per-panel cross sections.
- :math:`B_k` is the **outflow** of :eq:`characteristic-closure` per basis
  function: an emission density :math:`q = \sum_i q_i u_i` has the outflow
  :math:`\sum_i q_i B_k[u_i]`, by linearity, and the vector
  :math:`(B_k[u_i])_i` is what
  :meth:`~orpheus.derivations.continuous.characteristic.closure.LinePeriod.inflow`
  takes on its trailing basis axis
  (:meth:`TraversalRule.outflow
  <orpheus.derivations.continuous.characteristic.transport.TraversalRule.outflow>`,
  ``(..., 2, N)``, one row per traversal of the period, 0 where absent).
- :math:`A_k` is the **entry response**: a unit intensity entering the
  traversal produces the flux :math:`e^{-\tau_k(0, s)}` at :math:`s`, and
  :math:`A_k[u_i]` is that flux paired with :math:`u_i`, the test against
  which a Galerkin assembly pairs an inflow
  (:meth:`TraversalRule.entry_response
  <orpheus.derivations.continuous.characteristic.transport.TraversalRule.entry_response>`).

**Why** :math:`A_k = B_{\bar k}`. Read traversal :math:`k` backward with
:math:`s' = L_k - s`. The reversed traversal crosses the same slots in
the opposite order, so :math:`c_{\bar k}(s') = c_k(L_k - s')` and
:math:`\tau_{\bar k}(s', L_k) = \tau_k(0, L_k - s')`; substituting in
:math:`B_{\bar k}[u_i]` gives :math:`A_k[u_i]`. The code computes each
transit's two integrals once, read forward, and reading a transit
backward swaps its exit and its entry, so ``entry_response`` is the
outflow of the reversed traversal by selection: the identity holds
bitwise by design, and the gate that asserts it
(``test_the_entry_response_of_a_traversal_is_the_outflow_of_its_reverse``)
only pins the design. The evidence about :math:`A_k` is
``test_the_entry_response_is_the_integral_attenuated_from_the_entry``,
against mpmath.

.. implements:: characteristic-traversal-integrals
   :by: orpheus.derivations.continuous.characteristic.transport.TraversalRule.outflow

   **Implemented by** ``TraversalRule.outflow`` (:math:`B_k`),
   ``TraversalRule.entry_response`` (:math:`A_k`, by selection of the
   reversed reading), ``TraversalRule.optical_depth`` and
   ``LinePeriod.optical_depth`` (:math:`\tau_k`).

.. implements:: characteristic-traversal-integrals
   :by: orpheus.derivations.continuous.characteristic.transport.TraversalRule.entry_response

.. implements:: characteristic-traversal-integrals
   :by: orpheus.derivations.continuous.characteristic.transport.TraversalRule.optical_depth

.. implements:: characteristic-traversal-integrals
   :by: orpheus.derivations.continuous.characteristic.closure.LinePeriod.optical_depth

Every integral here is a polynomial in :math:`c(s)` times an exponential
in :math:`s`, and Gauss–Legendre in arc length is exact only for
polynomials in :math:`s`. Two features defeat one Gauss rule on a slot:
the exponential, which varies over one mean free path, and on a cylinder
or a sphere the orbit coordinate itself, which is not a polynomial in
:math:`s`. The rule cuts each slot into **pieces** on which both are
resolved, and evaluates every integral through one graded body.

.. _characteristic-pieces:

The pieces: exponential grading toward both ends of a slot
----------------------------------------------------------

A slot of the panel chord is a **member** when a transit traverses it,
and it belongs to exactly one transit. Each member slot is cut at its
midpoint, and each half is graded toward its own end of the slot: piece
ends at the distances :math:`2^{k}/\Sigma`, :math:`k = 0, \dots, 6`
(1, 2, 4, …, 64 mean free paths), clipped to the half, when the half's
optical width exceeds 2 (the module constant ``_THIN``). A half of optical
width at most 2 is one piece. So a slot of optical width at most 4 is two
pieces, its halves, and the rest of a thick half beyond 64 mean free paths
is one **middle piece**.

`[M]` 2026-10-06, the piece ends in mean free paths from the slot's
start, on a one-panel slab of width 1 crossed at :math:`\mu = 1`:

.. list-table::
   :header-rows: 1
   :widths: 14 86

   * - :math:`\Sigma l`
     - Piece ends
   * - 1, 3
     - 0, :math:`\Sigma l/2`, :math:`\Sigma l`: the halves only
   * - 5
     - 0, 1, 2, 2.5, 3, 4, 5
   * - 30
     - 0, 1, 2, 4, 8, 15, 22, 26, 28, 29, 30
   * - 1000
     - 0, 1, 2, 4, 8, 16, 32, 64, 500, 936, 968, 984, 992, 996, 998, 999,
       1000: two middle pieces of 436 mean free paths

**Why these numbers.** Inside a slot the cross section is constant, and
the flux an emission density produces there is a smooth part plus a
multiple of :math:`e^{-\Sigma s}`, the transient of what entered at the
slot's start; an integral attenuated to the slot's exit carries the
weight :math:`e^{-\Sigma(\ell - s)}` toward the other end. So the integrands
vary on the scale of one mean free path within a few mean free paths of
either end, and slowly elsewhere. Pieces of one mean free path at each
end, doubling inward, keep the fast variation on short pieces and the slow
remainder on wide ones. The doublings stop at 64 mean free paths because
the transient there is :math:`e^{-64} \approx 1.6 \times 10^{-28}` of its
starting value, far below double precision. Both ends are graded because
the transient and the exit weight sit at opposite ends. The Volterra
block's outer rule integrates on these pieces; the integrals attenuated to
a piece's end are the graded body of :ref:`characteristic-attenuated-integral`,
which does not rely on the piece being narrow.

.. _characteristic-branch-grading:

The pieces: hp grading toward the orbit coordinate's branch points
------------------------------------------------------------------

On a cylinder or a sphere the orbit coordinate along a line is
(:eq:`geometry-line-crossing-law`)

.. math::

   c(t)^{2} \;=\; b^{2} + |P\Omega|^{2}\,(t - t^{*})^{2},

with :math:`b` the impact parameter, :math:`|P\Omega|` the projected speed
and :math:`t^{*}` the parameter of the closest approach. Continued to
complex :math:`t`, :math:`c(t)` has two branch points, where
:math:`c^{2} = 0`:

.. math::

   t \;=\; t^{*} \pm i\,\frac{b}{|P\Omega|} .

For real :math:`t`,
:math:`|t - t^{*} \mp i b/|P\Omega||^{2} = (t - t^{*})^{2} + b^{2}/|P\Omega|^{2}
= c(t)^{2}/|P\Omega|^{2}`, so the distance from a point of the line to the
branch points is

.. math::

   \operatorname{dist}\bigl(t,\ \text{branch points}\bigr) \;=\; \frac{c(t)}{|P\Omega|} .

A basis function is a general polynomial in :math:`c`, odd powers
included, so the integrands are analytic in :math:`t` everywhere except at
the branch points (the even powers of :math:`c` are polynomials in
:math:`t`). On a slot :math:`c` is monotone, because the closest approach
is a slot end, so the point of the slot nearest the branch points is the
end with the smaller orbit coordinate, the slot's **near end**, at the
distance

.. math::

   D \;=\; \frac{c_{\rm near}}{|P\Omega|}:

:math:`D = b/|P\Omega|` on a slot ending at the closest approach, and
:math:`D = r_k/|P\Omega|` on a slot starting at the crossing of a radius
:math:`r_k`, a small cavity or a small inner region for instance.

Gauss–Legendre converges at a rate set by the largest Bernstein ellipse
around the interval, with foci at its ends, inside which the integrand is
analytic: with the ellipse parameter :math:`\rho_E` (the sum of its
semi-axes, the interval mapped to :math:`[-1, 1]`) and :math:`|f| \le M`
inside, the rule with :math:`n + 1` points errs by at most
:math:`64M / \bigl(15(1 - \rho_E^{-2})\,\rho_E^{2n+2}\bigr)`
:cite:`Trefethen2008` (Theorem 4.5, eq. (4.14), p. 77). One piece of
length :math:`\ell \gg D` sees the branch points almost at its end,
:math:`\rho_E \to 1`, and converges slowly.

So the pieces **halve toward the near end**: piece ends at the depths
:math:`\ell\,2^{-k}`, :math:`k = 1, \dots, 52`, from the near end, kept while
the depth exceeds :math:`D` (``TraversalRule._branch_edges``; 52 is
``np.finfo(float).nmant``, below which a depth is under one ulp of the
slot). A piece :math:`[\delta/2, \delta]` (depths from the near end) has
the width :math:`\delta/2` and lies at least :math:`\delta/2` from the
branch points, whatever :math:`b`; the innermost piece
:math:`[0, \delta_{\min}]` has :math:`D < \delta_{\min} \le 2D`. The ellipse parameter of every piece
is therefore bounded below independently of :math:`b`, each piece
converges at a fixed geometric rate in its point count, and the number of
pieces grows only like :math:`\log_2(\ell/D)`. This is **hp grading**:
geometric refinement toward a singularity with a fixed rule per piece.

`[M]` 2026-10-06, the worst ellipse parameter over every position of the
branch points the grading admits, mapped onto each piece (a direct
evaluation of :math:`\rho_E = |z \pm \sqrt{z^2 - 1}|`, the larger root):

.. list-table::
   :header-rows: 1
   :widths: 28 22 22 28

   * - Grading
     - Graded piece, :math:`\rho_E \ge`
     - Innermost piece, :math:`\rho_E \ge`
     - :math:`\rho_E^{-14}`, graded / innermost
   * - halving (the code)
     - :math:`3 + 2\sqrt2 \approx 5.83`
     - 2.89
     - :math:`1.9 \times 10^{-11}` / :math:`3.5 \times 10^{-7}`
   * - quartering
     - 3.00
     - 2.08
     - :math:`2.1 \times 10^{-7}` / :math:`3.5 \times 10^{-5}`

The last column is the bound's decay factor at 7 points, before its
constant: halving gains about four orders of magnitude on
every graded piece. The bound's constant :math:`M` is small on the
innermost piece, where :math:`c` itself is of the order of
:math:`c_{\rm near}`, which is why the measured errors below sit far
under the factor.

Through the centre (:math:`b = 0`) the two slots that touch the closest
approach have :math:`c_{\rm near} = 0`, so :math:`D = 0` and they are not
graded: on each of them :math:`c = |P\Omega|\,|t - t^{*}|` is linear. The
other slots of such a line are graded by the same rule toward the crossing
of their smaller radius, which costs pieces and changes no value (on that
line :math:`c` is linear on every slot). On the slab nothing is graded, and
a cylinder line parallel to the axis has no transit.

`[M]` 2026-10-06, :math:`\psi \cdot q` for :math:`q = c` against mpmath
at 7 points per piece (``probe_ratio2.py`` of
``scratch/characteristic_architecture/p1_step_b2/``, re-run here on the
present code), worst relative error over five points on the line, solid
bodies :math:`(0, 0.5, 1.5, 2.0)` at the resolution
:math:`(p, L, \rho) = (3, 2, 1/2)`:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Line
     - :math:`\psi`, halving
     - :math:`\psi`, quartering (arm T3b)
   * - sphere, :math:`b = 10^{-2}`
     - :math:`1.4 \times 10^{-14}`
     - :math:`5.0 \times 10^{-13}`
   * - cylinder, :math:`b = 10^{-2}`, :math:`\Omega_z = 0.8`
     - :math:`2.8 \times 10^{-14}`
     - :math:`8.9 \times 10^{-13}`
   * - sphere, :math:`b = 3 \times 10^{-3}`
     - :math:`3.9 \times 10^{-15}`
     - :math:`2.1 \times 10^{-12}`
   * - cylinder, :math:`b = 3 \times 10^{-3}`, :math:`\Omega_z = 0.8`
     - :math:`6.3 \times 10^{-15}`
     - :math:`3.9 \times 10^{-12}`

The gate's band is :math:`10^{-13}`, so quartering fails it at 7 points.
At 12 points both ratios meet it (``probe_ratio.py``, :math:`b = 10^{-2}`
and :math:`10^{-3}`: at most :math:`3.4 \times 10^{-15}` either way), and
at :math:`b = 10^{-4}` both stay below :math:`2 \times 10^{-15}` even at 7
points. The ratio is visible only at a small point count, which is why
``test_the_turning_grading_halves_toward_the_closest_approach_at_seven_points``
runs at 7, with the two :math:`b = 10^{-4}` lines as its declared
controls. Each line keeps its live pieces in chord order (slot by slot,
then along the slot), padded with dead pieces to the batch's largest
count.

.. dropdown:: First got wrong: a change of variable, and a coarser ratio
   :color: muted

   **The change of variable.** The step's design mapped each slot ending
   at the closest approach by :math:`c = b + (c_{\rm far} - b)\,u^{2}`,
   :math:`u \in [0, 1]`, so that :math:`c` is a polynomial in :math:`u`
   (the prototype's turning map). The integrand's Jacobian keeps a branch
   point: :math:`s \propto u\sqrt{(c_{\rm far} - b)(2b + (c_{\rm far} - b)u^{2})}`
   vanishes under the root at
   :math:`u = \pm i\sqrt{2b/(c_{\rm far} - b)}`, which at :math:`b = 10^{-4}`
   lies about 0.03 from the end of :math:`[0, 1]` on a turning slot reaching
   :math:`c_{\rm far} = 0.25`. One piece at 16 points missed by :math:`2.5 \times 10^{-12}` at :math:`b = 10^{-4}` and
   :math:`6.5 \times 10^{-13}` at :math:`b = 10^{-3}` (`[M]` 2026-10-06, the
   test-architect's ``probe_turn2.py``). The map also covered only the
   slots ending at the closest approach: a hollow sphere of cavity radius
   0.01 crossed through its cavity missed by :math:`8.5 \times 10^{-7}`
   (qa's ``probe_branch.py``), because there the near end is the crossing
   of the cavity surface, not a closest approach. The main agent retired
   the map for the grading above, which locates the singularity from
   :math:`c^{2} = 0` on every slot (recorded as ERR-099).

   **The ratio.** Pieces shrinking by 4 toward the near end instead of 2
   reach the branch point's scale in half as many pieces, and each piece
   :math:`[\delta/4, \delta]` sees the branch points at a third of its own
   width; the ellipse parameter falls from 5.83 to 3 and the error at 7
   points from :math:`3 \times 10^{-14}` to :math:`9 \times 10^{-13}` (the
   table above). At 12 points both meet the band, which is why the 16-point
   gates could not tell the ratios apart.

.. _characteristic-attenuated-integral:

One attenuated integral
-----------------------

Every quantity of the rule is built from one body,
``TraversalRule._attenuated``: the functions of a slot's panel integrated
between two distances :math:`a` and :math:`z` along the slot, attenuated
to the second,

.. math::

   I_m(a \to z) \;=\; \int_{[a, z]} u_m\bigl(c(s)\bigr)\,
   e^{-\Sigma\,|z - s|}\,\mathrm{d}s ,

on a rule graded exponentially toward :math:`z` (the same cut at
:math:`2^{k}` mean free paths, applied to the whole interval, with
``inner_points`` Gauss–Legendre points per interval). It has three uses:

- a piece's integral attenuated to its **end**,
  :math:`I(s_{\rm lo} \to s_{\rm hi})`, which the outflow and the carried
  flux are made of;
- a piece's integral attenuated to its **start**,
  :math:`I(s_{\rm hi} \to s_{\rm lo})`, which the entry response is made of;
- the integral from a piece's start to a point inside it,
  :math:`I(s_{\rm lo} \to s)`, which the Volterra block (at the piece's
  nodes) and the angular flux (at any point) need.

**Why graded inside a piece.** A piece is graded toward the ends of its
slot, not toward its own ends, and a middle piece of a thick slot is
hundreds of mean free paths wide. Its integral attenuated to its own end
lives in its last mean free path, and that integral is carried at full
weight to the flux just past the piece, so an ungraded rule on the
piece's own nodes misses it at order one. Grading each use toward the end
it is attenuated to resolves every case with one body.

.. dropdown:: First got wrong: the piece integrals on the piece's own nodes
   :color: muted

   The rule first integrated each piece's attenuated integrals on the
   piece's own ``points`` Gauss–Legendre nodes, and graded only the
   integral from a piece's start to a point inside it. The outflow read at
   a slot's exit was right, because a middle piece's contribution reaches
   the exit attenuated by :math:`e^{-64}` or more. The flux just past the
   middle piece was not. `[M]` 2026-10-06, qa's ``probe_thick_psi.py`` on
   that code: the flux just past the middle piece of a 1000-mean-free-path
   slot off by 0.54 relative, :math:`5.2 \times 10^{-6}` at 200 mean free
   paths, and :math:`f \cdot V g` off by :math:`2.4 \times 10^{-3}`. On the
   present piece layout the same defect (battery arm T14, re-dropped here)
   gives, at 1000 mean free paths, :math:`\psi` off by 0.25 half a mean free
   path past a piece start and :math:`f \cdot V g` by :math:`1.3 \times 10^{-3}`;
   at 200 mean free paths :math:`9.0 \times 10^{-11}` and
   :math:`1.6 \times 10^{-12}`; at 100 and fewer, nothing above
   :math:`3 \times 10^{-14}`. The main agent made every carried quantity one
   graded body (recorded as ERR-100).

Along a transit: the outflow, the entry response and the carried flux
---------------------------------------------------------------------

The pieces of a transit, in chord order, have optical depths
:math:`\tau_J = \Sigma_J (s_{J,\rm hi} - s_{J,\rm lo})`. Write
:math:`\tau^{\uparrow}_J` for the depth of the transit's pieces before
:math:`J` and :math:`\tau^{\downarrow}_J` for the depth of those after it.
Both are exclusive cumulative sums of non-negative terms, so no
difference of depths is ever formed. Read forward, the transit's two
integrals are

.. math::

   B[u] \;=\; \sum_J e^{-\tau^{\downarrow}_J}\, I(s_{J,\rm lo} \to s_{J,\rm hi}),
   \qquad
   A[u] \;=\; \sum_J e^{-\tau^{\uparrow}_J}\, I(s_{J,\rm hi} \to s_{J,\rm lo}),

each piece's vector scattered into the columns of its panel's functions
(:meth:`PanelBasis.columns
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.columns>`).
A traversal of the period reads its transit forward or reversed, and
reversing swaps the two, which is how ``outflow`` and ``entry_response``
select their rows.

What a transit's earlier pieces carry into the start of piece :math:`J`
is the vector

.. math::

   C_{J+1} \;=\; e^{-\tau_J}\, C_J + I(s_{J,\rm lo} \to s_{J,\rm hi}),
   \qquad C_{\rm first} = 0,

a scan over the pieces in chord order, one per transit. The scan is a
Python loop, and it is essential: the closed form
:math:`e^{-\tau^{\uparrow}_J}\sum_{J' < J} e^{+\tau^{\uparrow}_{J'+1}} I_{J'}`
overflows on a thick transit and cancels, while the recurrence multiplies
only by transmissions at most 1.

.. _characteristic-volterra:

The Volterra block
------------------

On one line read forward, with nothing entering, basis function
:math:`j` produces the flux

.. math::

   \psi_j(s) \;=\; \int_{\rm entry}^{s} u_j\bigl(c(s')\bigr)\,
   e^{-\tau(s', s)}\,\mathrm{d}s'

along each of the line's transits, a Volterra operator in arc length. Its
Galerkin block, accumulated over a batch of lines with the caller's
weights :math:`w_L`, is

.. math::

   V_{ij} \;=\; \sum_{L} w_L \sum_{\text{transits of } L}
   \int u_i\bigl(c(s)\bigr)\,\psi_j(s)\,\mathrm{d}s ,

(:meth:`TraversalRule.volterra
<orpheus.derivations.continuous.characteristic.transport.TraversalRule.volterra>`,
``(N, N)``, never stored per line). At the node :math:`s_q` of piece
:math:`J`, at the distance :math:`s_q - s_{J,\rm lo}` past its start, the
flux is what the earlier pieces carry in, attenuated from the start, plus
the piece's own integral up to the node:

.. math::

   \psi_j(s_q) \;=\; e^{-\Sigma_J (s_q - s_{J,\rm lo})}\,C_{J,j}
   \;+\; I_j(s_{J,\rm lo} \to s_q),

and :math:`V` gains :math:`w_L \sum_q W_q\,u_i(s_q)\,\psi_j(s_q)` with the
piece's ``points``-point Gauss–Legendre weights :math:`W_q`. The outer rule
integrates the product :math:`u_i \psi_j`, smooth on a piece because the
pieces are graded where :math:`\psi` has its layers; the inner integral is
the graded one.

Only the line's transits read forward enter: a line's flux lives on its
own transits in its own direction, and the reversed traversals of its
period carry its cycle and belong to the opposite line, which a batch
that holds it weights on its own. The inflow part of a line's Galerkin
block, :math:`\sum_k A_k \otimes \psi^{\rm in}_k` with :math:`\psi^{\rm in}`
from :eq:`characteristic-closure` applied to the :math:`B_k`, is not in
:math:`V`; composing the two over the measure on lines is the assembly
(:ref:`characteristic-what-is-not-built`).

**Reciprocity, and where it hides a transposition.** The block of the
reversed line is the transpose, :math:`V(-\Omega) = V(\Omega)^{\mathsf T}`
(``test_the_reversed_lines_triangle_is_the_transpose``, on the slab, where
the gate also checks that :math:`V(\Omega)` is not symmetric). On a
cylinder or a sphere the reversed line is the reflection of the line
through its closest approach, which maps its path in the orbit space onto
itself (:ref:`characteristic-period`), so :math:`V(-\Omega) = V(\Omega)`,
and with reciprocity :math:`V` is **symmetric on every radial line**. A
triangle with :math:`i` and :math:`j` transposed is then invisible on a
radial chord; the battery's arm T8 reddens the slab rows of the triangle
gate only, and the radial rows are declared in its stabiliser.

.. _characteristic-angular-flux:

The angular flux on a line
--------------------------

:meth:`TraversalRule.angular_flux
<orpheus.derivations.continuous.characteristic.transport.TraversalRule.angular_flux>`
reads the flux of every basis function at parameters :math:`t` on the
lines, ``(..., q, N)``, given the inflow at each traversal's entry
``(..., 2, N)`` that :meth:`LinePeriod.inflow
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.inflow>`
returns from the optical depths and the outflows. For a point at the
distance :math:`d` from the start of its slot, in piece :math:`J` of a
transit,

.. math::

   \psi(t) \;=\; e^{-(\tau^{\uparrow}_J + \Sigma_J (d - s_{J,\rm lo}))}\,
   \psi^{\rm in}_{k_f}
   \;+\; e^{-\Sigma_J (d - s_{J,\rm lo})}\,C_J
   \;+\; I(s_{J,\rm lo} \to d),

the inflow attenuated from the transit's entry, what the transit's
earlier pieces carry in, and the piece's own integral up to the point.
:math:`k_f` is the traversal that reads the point's transit forward,
:meth:`LinePeriod.forward_traversal
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.forward_traversal>`:
every transit a line makes is read forward by its period, since the first
traversal is transit 0 forward and the successor rule of
:eq:`characteristic-transit-rank` prefers a forward candidate, so a
shell's second transit is reached forward. ``forward_traversal`` raises
``RuntimeError`` if a transit is not read forward, which only a period
and a chord that disagree can produce.

The point's piece is the number of live piece starts at or before
:math:`t`, minus one. A point that is not on a transit of its line (before
its first piece, in a cavity, beyond a wall) is refused with
``ValueError``. The exit test compares :math:`t` with the parameter of the
slot's closing crossing, the kernel's own spelling of the wall, and the
distance is clipped to the slot's length (:ref:`characteristic-gotchas`,
the exit wall). The same function, the same
pieces and the same attenuated integral serve this reading and the
Volterra block, so the reading of a flux and the assembly share one rule.

The gates compare :math:`\psi` with mpmath at five fractions of each
transit with no inflow; at each transit's exit crossing, where with no
inflow it is the outflow; and with the closure, on mirrors and polished
walls of amplitudes 0.3 and 0.6, with the explicit backward march wall by
wall that the closure's own gate uses
(``test_the_closed_angular_flux_is_the_unfolded_backward_path``). With
every amplitude 0 the inflow is exactly 0 and :math:`\psi` equals the
no-inflow flux bitwise (``test_vacuum_walls_add_nothing_on_a_line``).

.. _characteristic-what-is-not-built:

What the package does not compute
=================================

The package answers, for any batch of lines, each line's period and
closure, the traversal integrals of every basis function, the vacuum
Volterra block with the caller's line weights and the angular flux at
points on the lines (:ref:`characteristic-transport`). It does not
integrate over lines. The measure on lines
(:eq:`geometry-measure-on-lines`,
:meth:`Chart.beam_density <orpheus.geometry.chart.Chart.beam_density>`)
is read by nothing in the package, and three things need it:

- the **assembly**, the Galerkin operator of the line part: each line's
  block :math:`V + \sum_k A_k \otimes \psi^{\rm in}_k`
  (:ref:`characteristic-volterra`) integrated over the lines;
- the **diffuse part of the resolvent**,
  :math:`U\,(I - T_w)^{-1} A\,U^{\mathsf T}` (:ref:`characteristic-resolvent`),
  whose escape functional :math:`U` and wall-to-wall transmission
  :math:`T_w` are integrals over the same measure.
  :class:`~orpheus.derivations.continuous.characteristic.walls.Walls`
  reads the diffuse amplitude and nothing consumes it;
- the **questions** on the dense pencil
  (:ref:`verification-reference-kernel`) and the **reading** at a point,
  whose first leg runs from the point to its first wall and is integrated
  over directions.

The user's ruling on the second rung's sketch (2026-10-06, Q1) put the
line rule, the assembly and the diffuse part (``WallCoupling``) in one
rung, the next, because all three are integrals over that measure. The
plan is ``.claude/plans/characteristic_reference_architecture.md`` ("P1
step (b), second rung: API sketch", item 0, and "P1 API sketch", items 3
to 7); the campaign's issue is #405.

The trajectory-resolvent family (:ref:`theory-trajectory-resolvent`)
is the reference every consumer reads today.


.. _characteristic-refuted:

Designs that do not work, and why
=================================

Each row is a design that reads as natural and is wrong, with the
structural reason it fails, so that no later design re-derives it.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - The design
     - Why it fails
   * - "The rank of a line is the number of distinct walls its transit
       touches" (the first-pass design, D1)
     - On the periodic slab a transit touches two distinct walls and the
       rank is 1: the wrap identifies the two faces, and the line
       continues without reflecting. The rank is the number of walls in
       the period after the deck maps identify them, which the successor
       rule of :eq:`characteristic-transit-rank` counts. The
       verification spec's row A6 asserts it; read as distinct walls, the
       periodic slab would be rank 2 with the transit reversed at each
       face, the periodic face a mirror.
   * - "The period is the transit and the reversed transit" (D1)
     - True of the slab between two returning walls only. A shell's line
       through the cavity has the period (transit 0, transit 1), both
       forward; a solid body's line and the periodic slab have a period
       of one traversal. The successor rule of
       :eq:`characteristic-transit-rank` is the form the code computes.
   * - "The albedo :math:`\alpha_k` of the wall transit :math:`k` ends on"
       (the first-pass closed form, D2)
     - Ambiguous on a backward characteristic: the wall a transit ends on
       read forward is the wall the backward path starts from. The
       closed form equals the unfolded path only with the amplitude of
       the wall at which the backward path reflects into the traversal
       (:ref:`characteristic-closure-section`, "The albedo pairing"); the
       other reading moves :math:`\psi` by 0.335 on the review's fixture.
   * - Albedos paired with walls by position, as a tuple inner first
       (``specular_albedos(geometry)`` of the old family)
     - A solid body's only wall would answer at position 0 for
       breakpoint 0, which is no wall, and the inner and outer albedos
       of a shell can be swapped with no error. Walls are keyed by
       breakpoint index (:ref:`characteristic-walls`); the battery's arm
       W9 (lookup by position) reddens
       ``test_a_body_has_a_wall_at_each_boundary_point_and_none_elsewhere``.
   * - Hoisting production's tag parse
       (``orpheus.transport.method._law_from_tag``) into
       ``orpheus.geometry.boundary``, shared by the reference and
       production (recommended in the step's design, Q1)
     - Refused by the user (:ref:`characteristic-walls-registry`): a
       shared parse makes every production change to a tag's meaning a
       change to the reference. The reference keeps its own registry and
       the duplicate is declared.
   * - A mirror deck with a partial ``ScalarResponse(a)`` read as a
       specular amplitude :math:`a` (the step's design table)
     - A mirror is a symmetry statement and returns everything; a
       partial specular wall is a surface, the identity deck with
       ``SpecularReemission(a)``. Reading the partial scalar under a
       mirror would make "a symmetric body" and "a polished wall" one
       declaration with two meanings. The built table refuses the cell
       (:ref:`characteristic-walls-factors`).
   * - Deciding at build whether a wall summing a specular and a diffuse
       law (``LawSum``) is served
     - Dissolved: a geometry admits a ``BC`` tag or a
       ``BoundaryTraceLaw`` per boundary point, and ``LawSum`` and
       ``LawScaled`` are content nodes, not laws, so no such wall can be
       declared. A ``Wall``'s two amplitude fields can hold such a wall.
   * - Preferring the reversed candidate over the forward one
     - Not wrong in value (both candidates trace one orbit-space path),
       but it doubles the period on the radial charts, so the rank is no
       longer the minimal period; arm C4 reddens 92 rows.
   * - ``1 - Pi`` formed as ``1 - exp(...)``
     - Loses the digits of a nearly lossless line: a relative error of
       :math:`6.3 \times 10^{-6}` on :math:`1 - \Pi` at
       :math:`\sum\tau = 3.5 \times 10^{-12}`
       (:ref:`characteristic-closure-section`).
   * - The change of variable :math:`c = b + (c_{\rm far} - b)u^{2}` on a
       slot ending at the closest approach (the second rung's sketch,
       item 3)
     - It makes :math:`c` a polynomial in :math:`u` and leaves the arc
       length's branch point at :math:`u = \pm i\sqrt{2b/(c_{\rm far} - b)}`,
       near the interval at small :math:`b` (:math:`2.5 \times 10^{-12}` at
       :math:`b = 10^{-4}`, 16 points). It also covers only the slots
       ending at the closest approach, while the branch points of
       :math:`c(t)` sit :math:`c_{\rm near}/|P\Omega|` from every radial
       slot's near end (:math:`8.5 \times 10^{-7}` beside a cavity of radius
       0.01). The grading toward the near end replaces it
       (:ref:`characteristic-branch-grading`, ERR-099).
   * - Pieces shrinking by 4 toward the branch points
     - A graded piece's ellipse parameter falls from 5.83 to 3, and the
       error at 7 points rises from :math:`3 \times 10^{-14}` to
       :math:`9 \times 10^{-13}`, outside the gates' band; at 12 points
       the two ratios cannot be told apart
       (:ref:`characteristic-branch-grading`).
   * - Each piece's attenuated integrals on the piece's own Gauss–Legendre
       nodes
     - Misses the integral carried out of a wide middle piece at order one
       (the flux just past it off by 0.54 at 1000 mean free paths), while
       the outflow at the slot's exit, attenuated by :math:`e^{-64}` or
       more from that piece, hides it. One graded body serves every
       carried quantity (:ref:`characteristic-attenuated-integral`,
       ERR-100).
   * - ``TraversalRule.of(period, basis, sigma_per_panel, resolution)``
       (the second rung's sketch)
     - Leaves two constructions spellable: a period chorded through
       another partition, and walls re-keyed onto a partition that is not
       the basis's. ``of(lines, basis, walls, sigma_t, ...)`` chords
       through the basis's own partition (the elegance review's finding
       C2; :ref:`characteristic-transport`).
   * - The flux at a wall located by the slot's start plus its length
     - The sum rounds an ulp apart from the closing crossing's parameter,
       which is how the kernel spells the wall, and refused 64 of 119
       reads at an exit wall (qa, 2026-10-06). The exit test reads the
       crossing parameter (:ref:`characteristic-angular-flux`).


.. _characteristic-evidence:

Numerical evidence
==================

The gates of the walls and the closure
--------------------------------------

`[M]` 2026-10-06, ``.venv/bin/python -O -m pytest -p no:cacheprovider``
over the two files at the branch ``feature/characteristic-walls-closure``:
175 rows, all passing. The level and label of each test function:

.. list-table::
   :header-rows: 1
   :widths: 52 8 14 26

   * - Test (``tests/gates/derivations/``)
     - Rows
     - Level
     - Claim
   * - ``test_characteristic_closure.py::test_the_period_matches_the_hand_counted_table``
     - 73
     - L0
     - :eq:`characteristic-transit-rank`
   * - ``…closure.py::test_the_period_is_the_physically_unfolded_path``
     - 9
     - L0
     - :eq:`characteristic-transit-rank`, and the reflection invariants
   * - ``…closure.py::test_the_inflow_is_the_unfolded_wall_by_wall_sum``
     - 16
     - L0
     - :eq:`characteristic-closure`
   * - ``…closure.py::test_an_absorbing_wall_zeroes_exactly_the_inflow_it_feeds``
     - 3
     - L0
     - :eq:`characteristic-closure`, the pairing zero
   * - ``…closure.py::test_every_amplitude_zero_adds_nothing``
     - 3
     - L0
     - :eq:`characteristic-closure`, the vacuum edge
   * - ``…closure.py::test_the_least_solution_on_a_lossless_trapped_line``
     - 1
     - L0
     - :eq:`characteristic-closure`, the least solution
   * - ``…closure.py::test_a_nearly_lossless_line_keeps_its_digits``
     - 3
     - L0
     - :eq:`characteristic-closure`, the ``expm1`` form
   * - ``test_characteristic_walls.py::test_every_law_reads_as_the_wall_its_physics_names``
     - 10
     - L0
     - :eq:`characteristic-closure`, the amplitudes :math:`a_k`
   * - ``…walls.py::test_a_served_factor_pair_reads_its_wall``
     - 4
     - L0
     - :eq:`characteristic-closure`, the amplitudes :math:`a_k`
   * - ``…walls.py::test_a_tag_reads_the_same_walls_as_the_law_it_names``
     - 5
     - L0
     - :eq:`characteristic-closure`, the tag and its law
   * - the optical depth (three functions), the batch, the refusals and
       guards of both files, the registry's own rows
     - 48
     - foundation
     - software invariants

The L0 rows of :eq:`characteristic-transit-rank` total 82 and those of
:eq:`characteristic-closure` 45. Each closure row compares
``LinePeriod`` with an independent spelling: the hand-counted table, the
3-D line reflected by Householder in the test, closed-form segment
lengths in mpmath, and the unfolded wall-by-wall sum in mpmath that
marches the backward path one wall at a time with no geometric-series
division, at 16 ulp relative on 40 seeded draws per amplitude pair with
exact zeros of :math:`\tau` mixed in. The walls rows compare with walls
written by hand from each law's physics, not from its factors.

The mutation battery of the walls and the closure
-------------------------------------------------

`[M]` 2026-10-06, the test-architect's battery
``scratch/characteristic_architecture/p1_step_b1/battery/`` (driver
``b1_mutants.py``, results ``summary.txt``): 34 arms, each a one-line
replacement in ``walls.py`` (W) or ``closure.py`` (C) installed in the
process pytest runs, the two files run per arm. Reds out of 175:

.. list-table::
   :header-rows: 1
   :widths: 40 10 50

   * - Arm
     - Reds
     - What it breaks
   * - PC (positive control: every amplitude times 0.999)
     - 80
     - every reading that depends on an amplitude; the control that must redden many rows
   * - C2 rank always 2
     - 58
     - the derived rank
   * - C4 prefer the reversed candidate
     - 92
     - the minimal period
   * - C1 pairing read at the entry wall
     - 26
     - the albedo pairing in the period
   * - C15 inflow pairing swapped (no roll)
     - 20
     - the albedo pairing in the cycle
   * - C14 amplitude not applied
     - 21
     - :math:`r_k = a_k B_k`
   * - C16 gain not rolled
     - 18
     - the once-around term's gain
   * - C10 the ERR-035 denominator (:math:`a_0^2 e^{-2\tau_0}`)
     - 17
     - the cycle product
   * - C9 once-around cross term dropped
     - 15
     - the rank-2 numerator
   * - C7 :math:`\Sigma_t` read in reversed region order
     - 12
     - the optical depth
   * - W1 wrap read as a mirror
     - 12
     - the wrap's partner
   * - W5 vacuum with default amplitude 1
     - 11
     - the vacuum row
   * - C13 refusal widened to any lossless traversal
     - 7
     - the least solution
   * - C5 reversed traversals dropped
     - 6
     - the reversal
   * - W2 Lambertian read as specular; W10 an undeclared parameter dropped
     - 5 each
     - the factor table; strict parameters
   * - W6 radial-wrap check removed; C12 plain division
     - 4 each
     - the refusal; the trapped line's NaN
   * - W4, W7, W8, W12, C8
     - 3 each
     - the vacuum guard, the wrap-back check, the range check, an unknown
       tag read as vacuum, ``1 - exp``
   * - C6 :math:`\tau` over every slot; C17 length check removed; W13
       white sign negated
     - 2 each
     - the transit's slots; the cross-section length; the white law's sign
   * - W3, W9, W11, C11, C18, C19, C20
     - 1 each
     - the source check, lookup by position, ``partial`` without a shape,
       the trapped refusal, the void product, the two ``RuntimeError``
       guards
   * - C7b exteriors given :math:`\Sigma_t = 5`
     - 0
     - declared null: no transit traverses an exterior slot

Every arm meant to redden reddens at least one row, and the one arm
declared null (C7b) stays green, for the structural reason stated in
:ref:`characteristic-closure-section`. The summary's line for C4 is
blank; its log reads 92 failed, 83 passed.

The gates of the basis and the transport
----------------------------------------

`[M]` 2026-10-06, ``.venv/bin/python -O -m pytest -p no:cacheprovider``
over the package's four gate files on the branch
``feature/characteristic-basis-transport``: 503 rows, all passing, in
179 s. The two files of the second rung hold 328 of them,
``test_characteristic_basis.py`` 225 and
``test_characteristic_transport.py`` 103. Every row is ``foundation``:
:eq:`characteristic-traversal-integrals` is the rung's label, and the
verification specification
(``scratch/characteristic_architecture/p1_step_b2/spec.md``) assigns the
rows marked † to ``l0`` under it and the row marked ‡ to ``l0`` under
:eq:`characteristic-closure`, once their markers move.

.. list-table::
   :header-rows: 1
   :widths: 56 8 36

   * - Test (``tests/gates/derivations/``)
     - Rows
     - Reference
   * - ``test_characteristic_basis.py::test_every_breakpoint_is_a_panel_end_and_each_panel_lies_in_one_region``
     - 32
     - the body's breakpoints
   * - ``…basis.py::test_panels_grade_toward_walls_and_interfaces_and_not_toward_a_singular_stratum``
     - 32
     - the grading law written by hand in mpmath
   * - ``…basis.py::test_the_nodes_are_each_panels_gauss_legendre_points``
     - 32
     - the roots of :math:`P_{p+1}` in mpmath
   * - ``…basis.py::test_the_functions_are_cardinal_at_the_nodes_and_sum_to_one``
     - 32
     - the identity
   * - ``…basis.py::test_the_interpolant_reproduces_every_per_region_polynomial_of_degree_p``
     - 32
     - mpmath; degree :math:`p + 1` as the loading control
   * - ``…basis.py::test_the_volume_density_is_the_charts_by_hand``
     - 8
     - the three densities written by hand
   * - ``…basis.py::test_the_mass_matrix_is_the_volume_integral_of_each_product``
     - 18
     - Lagrange products integrated in mpmath
   * - ``…basis.py::test_each_panels_mass_is_its_chart_measure``
     - 32
     - ``Chart.measure`` and the closed-form volume
   * - ``…basis.py``, the signature row, the refusals and the walls re-keyed
       (four functions)
     - 7
     - the walls written by hand; the refusal fragments
   * - ``test_characteristic_transport.py::test_the_panel_chord_unfolds_into_the_hand_counted_period``
     - 12
     - a hand-counted period; the body chord's depth
   * - ``…transport.py::test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit`` †
     - 24
     - mpmath line integrals, per-region cubics, two groups
   * - ``…transport.py::test_each_basis_functions_outflow_is_its_line_integral`` †
     - 3
     - mpmath, each basis function
   * - ``…transport.py::test_the_entry_response_is_the_integral_attenuated_from_the_entry`` †
     - 12
     - mpmath
   * - ``…transport.py::test_the_entry_response_of_a_traversal_is_the_outflow_of_its_reverse``
     - 3
     - itself: the design identity, bitwise
   * - ``…transport.py::test_a_thousand_mean_free_path_slot_integrates_to_the_closed_form`` †
     - 2
     - the closed form, 4 ulp; mpmath
   * - ``…transport.py::test_a_void_region_integrates_its_source_unattenuated`` †
     - 1
     - the closed form :math:`q\ell`; mpmath
   * - ``…transport.py::test_a_turning_slot_integrates_the_square_root_at_the_closest_approach`` †
     - 4
     - mpmath (tanh-sinh, split at :math:`t^{*}`)
   * - ``…transport.py::test_a_slot_starting_at_a_small_radius_integrates_its_branch_point``
     - 5
     - mpmath; the rows the defect of ERR-099 reddens (four lines at one
       panel per region, and cavity 0.001 at the gates' resolution)
   * - ``…transport.py::test_a_small_radius_hidden_by_a_panel_end_is_the_err099_control``
     - 3
     - mpmath; the declared controls of ERR-099, where a panel end of the
       basis sits at the small radius
   * - ``…transport.py::test_the_turning_grading_halves_toward_the_closest_approach_at_seven_points``
     - 6
     - mpmath, at 7 points per piece; two controls at :math:`b = 10^{-4}`
   * - ``…transport.py::test_the_volterra_triangle_is_the_double_integral_along_the_line``
     - 6
     - an mpmath double integral, converged to :math:`10^{-18}` between 24
       and 48 points per piece
   * - ``…transport.py::test_the_reversed_lines_triangle_is_the_transpose``
     - 2
     - reciprocity, with a non-symmetry leg
   * - ``…transport.py::test_the_triangle_of_a_batch_is_the_weighted_sum_of_its_lines``
     - 1
     - the per-line blocks
   * - ``…transport.py::test_a_thick_slot_is_read_by_the_triangle_and_by_psi_in_closed_form``
     - 2
     - closed forms in mpmath, no quadrature
   * - ``…transport.py::test_the_vacuum_angular_flux_is_the_source_integral_since_the_entry``
     - 4
     - mpmath at five fractions of each transit
   * - ``…transport.py::test_the_angular_flux_is_read_at_a_transits_exit_crossing``
     - 6
     - mpmath outflow at the exit crossing
   * - ``…transport.py::test_the_closed_angular_flux_is_the_unfolded_backward_path`` ‡
     - 4
     - the explicit backward march, wall by wall, in mpmath
   * - ``…transport.py::test_vacuum_walls_add_nothing_on_a_line``
     - 3
     - itself, bitwise

The references are written in ``tests/gates/derivations/_characteristic_mp.py``
in mpmath, which imports nothing from ``orpheus``; every traversal a row
reads is written by hand per fixture, never taken from the period under
test. Two inputs come from the code and are declared: the panel ends (the
space being tested, gated on their own by the first two basis rows) and
the nodes, at which a per-region polynomial is sampled into coefficients
(the interpolant reproduces it whatever the nodes, and the nodes are
gated against mpmath roots). The tolerances on thin panels carry the
factor :math:`1 + \kappa_P`, :math:`\kappa_P = \max(|a|, |b|)/(b - a)`,
because the local coordinate of a panel 5e-6 wide at :math:`r = 2` loses
five digits whatever the code does.

The mutation battery of the basis and the transport
---------------------------------------------------

`[M]` 2026-10-06, the test-architect's battery
``scratch/characteristic_architecture/p1_step_b2/battery/`` (plugin
``battery_plugin.py``, results ``summary.txt``, round 3, run when the two
files held 314 rows): 26 arms, each an in-process textual mutant of
``basis.py`` (B), ``walls.py`` (W) or ``transport.py`` (T), the pristine
copies matching the tree under ``diff -q`` afterwards. The honest run is
green; every arm reddens its target row.

.. list-table::
   :header-rows: 1
   :widths: 46 8 46

   * - Arm
     - Reds
     - Target rows reddened
   * - B3 depth taken from the wrong end
     - 190
     - the grading law and everything downstream
   * - B11 a panel straddling an interface (guards off)
     - 158
     - the partition, the interpolant, the transport rows
   * - B4 nodes one short
     - 155
     - the nodes, the cardinal property, the interpolant
   * - B5 reference nodes reversed; B6 local map reversed
     - 118 each
     - the cardinal property, the interpolant, each function's outflow
   * - W1 breakpoint :math:`n` kept as the last index; W2 the wrap's partner
       not re-keyed
     - 87 each
     - the re-keyed walls, the period on the panel chord, the transport
       rows
   * - T1 attenuated from the wrong end (the positive control)
     - 62
     - the outflow, each function's outflow, the entry response, the
       angular flux
   * - B9 mass blocks misplaced; B7 density without :math:`d`
     - 49; 44
     - the mass matrix; the density and the mass rows
   * - B2 inner wall ungraded; B1 centre graded
     - 14; 12
     - the grading law
   * - T12 carried term unattenuated
     - 12
     - the vacuum and closed angular flux, the small-point-count rows
   * - B8 mass rule one point short
     - 11
     - the mass matrix and the panel measure
   * - T3a the grading toward the branch points removed
     - 10
     - the turning rows at :math:`b = 10^{-4}`, the small-radius and the
       small-point-count rows
   * - T6 :math:`A` computed as :math:`B`
     - 8
     - the entry response and the design identity
   * - T9 the carried rows of the triangle dropped
     - 6
     - the triangle
   * - T11 inflow of the wrong traversal; T13 inflow over-attenuated
     - 4 each
     - the closed angular flux
   * - T2 no exponential grading; T3b the branch grading's ratio 1/4; T8
       triangle transposed
     - 2 each
     - the 1000-mean-free-path slot; the small-point-count rows at :math:`b = 10^{-2}`;
       the triangle's slab rows (the radial rows are in its stabiliser)
   * - B10 refinement unchecked; W3 chart unchecked; T7 inner rule from the
       slot's start; T10 weights flipped
     - 1 each
     - the refusal; the chart refusal; the triangle's thick slab row; the
       batch row

Arms on rows added after round 3, re-dropped in this pass (`[M]`
2026-10-06, the same plugin, ``test_characteristic_transport.py`` without
the six slow triangle rows):

- **T14**, each piece's attenuated integrals on the piece's own nodes
  (the defect of ERR-100): 2 red of 91, both rows of
  ``test_a_thick_slot_is_read_by_the_triangle_and_by_psi_in_closed_form``.
- **T16**, the exit test by the slot's start plus its length: 3 red of 91,
  of the 6 rows of ``test_the_angular_flux_is_read_at_a_transits_exit_crossing``.
- **T3a** on the same file: 7 red of 91, the two turning rows at
  :math:`b = 10^{-4}`, the four small-point-count rows of the time and one
  small-radius row.
- **The defect of ERR-099**, the branch grading restricted to the slots
  ending at the closest approach. On the file as it first stood, with the
  small-radius lines at the gates' resolution :math:`(3, 2, 1/2)` only, it
  reddened 1 of 91 rows,
  ``[sphere_cavity_1e-3_b_half_r0]`` (:math:`B_0` off by
  :math:`5.3 \times 10^{-12}`): the basis's own grading toward the cavity
  wall or the interface puts a panel end within a few :math:`r_k` of the
  small radius, and the other three lines miss by at most
  :math:`2.3 \times 10^{-14}`. With one panel per region all four miss, by
  :math:`4 \times 10^{-9}` to :math:`1.1 \times 10^{-7}`. The test-architect
  then split the small-radius gate: the five rows the defect reddens
  (battery arm T17, the four one-panel rows and cavity 0.001 at the gates'
  resolution) carry ``catches("ERR-099")``, and the three it leaves green
  are its declared controls,
  ``test_a_small_radius_hidden_by_a_panel_end_is_the_err099_control``.

Rows no arm reddens, each declared in the specification: the design
identity :math:`A_{\rm fwd} = B_{\rm rev}`, the reciprocity row, the void
row's warning guards, the amplitude-0 row and the signature row.


.. _characteristic-gotchas:

Gotchas
=======

- **A wall is not a position.** Index the walls with the breakpoint a
  transit reports; a solid body's one wall is at :math:`n`, and asking
  at :math:`0` raises ``IndexError``.
- **The partner is not the exit.** The next traversal enters at the
  partner of the wall the last one exits at: the same wall under a
  mirror, the opposite face under a wrap.
- **The amplitude belongs to the exit wall and multiplies the next
  inflow.** A rewrite of ``inflow`` that drops the roll pairs each
  traversal with its own exit amplitude, wrong by 45 % on the review's
  fixture and invisible on every body with equal amplitudes.
- **Pass optical depths, not transmissions.** ``inflow`` takes
  :math:`\tau_k` so that it can form :math:`1 - \Pi` without
  cancellation; passing :math:`-\log` of a rounded transmission loses
  the digits the form exists to keep.
- **A trapped line refuses only with a source.** A void region is
  admitted; a source-free trapped line reads 0. A ``TrappedSource`` on a
  posed problem means an external source sits on a line that never meets
  material and never loses anything at a wall.
- **The registry is the reference's.** Do not route a reference tag
  through production's parse or the reverse; the two differ on purpose
  (``albedo`` is a production kind and a refused one here; production
  drops the white albedo, #583).
- **Diffuse amplitudes are read and not consumed.** Until the diffuse
  part of the resolvent exists, a white wall's return is absent from
  every line closure; the line part alone treats it as absorbing.
- **Read the flux at a wall at the kernel's crossing parameter.**
  ``angular_flux`` locates a point by its slot's closing crossing. A
  parameter formed as the slot's start plus its length can round an ulp
  beyond the wall and is then refused as off the transit. The battery's
  arm T16, which restores the test by start plus length, reddens 3 of the
  6 rows of ``test_the_angular_flux_is_read_at_a_transits_exit_crossing``.
- **A transposed triangle is invisible on a radial line.** The Volterra
  block is symmetric on every cylinder and sphere line
  (:ref:`characteristic-volterra`), so a gate of the triangle's
  orientation needs a slab row.
- **The entry response equals the reversed outflow by design.**
  ``entry_response`` selects the reversed reading of the same two transit
  integrals that ``outflow`` reads, so a bitwise comparison of the two
  cannot fail on a wrong :math:`A_k`. The evidence about :math:`A_k` is
  its comparison with mpmath.
- **The cross sections are per region at the factory and per panel on the
  value.** ``TraversalRule.of`` takes one total cross section per region
  and reads it onto the panels; the direct constructor takes one per
  panel.
- ``Walls.on`` **trusts the refinement.** It re-keys breakpoint :math:`n`
  to the panel count and holds no positions to check against; call it
  with a basis's own partition, as ``TraversalRule.of`` does.
- **A point read by** ``angular_flux`` **lies on a transit of its line.** A
  point in a cavity, beyond a wall, or on a line with no transit is
  refused.
- **The thin threshold is per half slot.** Every slot is cut at its
  midpoint, and each half is graded only when its optical width exceeds
  2: a slot of optical width up to 4 is two pieces.


.. _characteristic-history:

History
=======

.. list-table::
   :header-rows: 1
   :widths: 14 58 14 14

   * - Date
     - Decision
     - Commit
     - Issue
   * - 2026-10-06
     - The first rung of the characteristic reference: the walls, read
       from the laws' factors with the reference's own tag registry (the
       user refused hoisting production's parse), and the line part of
       the boundary resolvent, the period chained through the walls'
       partners and the least solution of its cycle. The design's
       distinct-wall rank and positional albedo were replaced by the
       derived period and the breakpoint-keyed walls.
     - ``72199f9d``
     - #405, #583
   * - 2026-10-06
     - The second rung: the panel basis (discontinuous nodal panels
       graded toward walls and interfaces, the volume density derived
       from the measure's one definition), the walls re-keyed onto the
       panel partition, and the transport along a line
       (``TraversalRule``: the traversal integrals, the Volterra block and
       the angular flux, from one graded attenuated integral). In review
       the sketch's change of variable at the closest approach gave way to
       hp grading toward the branch points, and the piece integrals were
       graded (ERR-099, ERR-100). The ladder was re-cut: the line rule,
       the assembly and the diffuse part move to the next rung.
     - ``ff979520``
     - #405
