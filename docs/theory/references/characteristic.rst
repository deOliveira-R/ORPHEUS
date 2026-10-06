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
      concept: characteristic reference, boundary resolvent, walls, line period, line closure
      role: "the closed reference that integrates transport along the lines of a 1-D concentric body (slab, cylinder, sphere, solid or hollow); this page holds its walls (each boundary point with what its law returns, read from the law's factors) and the line part of its boundary resolvent (the period of each line's unfolded path and the least solution of its cycle)"
      code: [orpheus.derivations.continuous.characteristic.walls, orpheus.derivations.continuous.characteristic.closure]
      depends_on: [chart_and_chord, boundary_conditions, reference_solutions]
      related: [trajectory_resolvent, layering]


Key facts
=========

- **What this is.** The characteristic reference is the closed reference
  that solves transport on a 1-D concentric body by integrating along
  the body's lines: :mod:`orpheus.derivations.continuous.characteristic`.
  It is built rung by rung beside the trajectory-resolvent family
  (:ref:`theory-trajectory-resolvent`), the family it is built to
  replace. What exists is its first rung: the **walls**
  (:mod:`~orpheus.derivations.continuous.characteristic.walls`) and the
  **line part of the boundary closure**
  (:mod:`~orpheus.derivations.continuous.characteristic.closure`). There
  is no source integral, no basis, no assembly, no question and no
  reading yet, so the package answers no eigenvalue and no flux
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
- **Evidence** `[M]` 2026-10-06: 175 gate rows in two files, 127 of them
  claims on this page's two labels at L0 and 48 software invariants and
  refusals; a 34-arm mutation battery with a positive control reddening
  80 rows, every arm meant to redden reddening at least one row, and the
  one arm declared null staying green (:ref:`characteristic-evidence`).


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
reading of a line, which this rung does not compute.


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


.. _characteristic-what-is-not-built:

What this rung does not compute
===============================

The closure needs the outflows :math:`B_k`, and this rung does not
compute them: :math:`B_k` integrates the emission density along a
traversal, so it needs the basis the emission is represented on and the
Galerkin assembly over the measure on lines, the next rungs of the plan
(``.claude/plans/characteristic_reference_architecture.md``, "P1 API
sketch", items 4 and 5). ``inflow`` is tested on abstract :math:`B_k`,
including a trailing basis axis, the shape an assembly over a basis
hands it. Nothing in the package reads a flux at a point or answers a
question; the first leg (from a point to the first wall) and the reading
at a point are the plan's item 7.

The diffuse part of the resolvent (the white and the Lambertian walls)
is not built either. :class:`~orpheus.derivations.continuous.characteristic.walls.Walls`
reads the diffuse amplitude and nothing consumes it. The design, the
finite-rank update of :ref:`characteristic-resolvent`, is recorded in
the plan (item 3, ``WallCoupling``) and in the ledger entry "the user, on
P1's API sketch".

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


.. _characteristic-evidence:

Numerical evidence
==================

The gates
---------

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

The mutation battery
--------------------

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
