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
      concept: characteristic reference, boundary resolvent, walls, line period, line closure, panel basis, even basis at a singular stratum, traversal integrals, Volterra block, hp grading, Galerkin assembly over lines, line rule, white-wall coupling, region cross sections, emission support, multigroup Galerkin system, emission space, k pencil, source pencil, adjoint flux, one-group Rayleigh-Ritz bound, the door, resolution, projection onto the panel basis, role arrows, response as the transposed problem, eigen gauge, flux integral by pairing, grading law at every impact-panel top, weighted line set, reading at a point, point measure, iterated Galerkin, diffuse walls' currents, angular flux at a point
      role: "the closed reference that integrates transport along the lines of a 1-D concentric body (slab, cylinder, sphere, solid or hollow); this page holds its walls (each boundary point with what its law returns, read from the law's factors), the line part of its boundary resolvent (the period of each line's unfolded path and the least solution of its cycle, with an arriving flux), its panel basis (even at a singular stratum) and the transport along one line on that basis (the traversal integrals, the vacuum Volterra block and the angular flux), and one group's transport block: the Galerkin assembly over the lines of the chart's line domain, the white walls' coupling, and the line rule graded from the group's optical scale; and the multigroup Galerkin system on the emission density (the cross sections per region, each group's exact emission support, the k pencil and the source pencil, the adjoint flux by transposition, the one-group Rayleigh-Ritz bound); and the door, the reference posed from a specification and a resolution (the projection of a mesh-free function onto the panel basis, the role arrows a source and a detector enter by, the response as the forward problem on the transposed cross sections, the eigen gauge, the flux integral as a pairing); and the reading at a point, the transported emission integrated over the line domain's lines through the point under the point's measure, with the diffuse walls folded in through their currents, and the angular flux at a point"
      code: [orpheus.derivations.continuous.characteristic.walls, orpheus.derivations.continuous.characteristic.closure, orpheus.derivations.continuous.characteristic.basis, orpheus.derivations.continuous.characteristic.transport, orpheus.derivations.continuous.characteristic.assembly, orpheus.derivations.continuous.characteristic.lines, orpheus.derivations.continuous.characteristic.reading, orpheus.derivations.continuous.characteristic.grading, orpheus.derivations.continuous.characteristic.cross_sections, orpheus.derivations.continuous.characteristic.system, orpheus.derivations.continuous.characteristic.reference, orpheus.derivations.common.eigenvalue, orpheus.derivations.common.angular_measure]
      depends_on: [chart_and_chord, boundary_conditions, reference_solutions]
      related: [characteristic_origins, layering]


Key facts
=========

- **What this is.** The characteristic reference is the closed reference
  that solves transport on a 1-D concentric body by integrating along
  the body's lines: :mod:`orpheus.derivations.continuous.characteristic`.
  It was built rung by rung beside the trajectory-resolvent family
  (:ref:`characteristic-origins-history`), and replaced it: P1 step (e2)
  (2026-10-10) deleted that family, after its tests were re-posed here
  (:ref:`characteristic-successors`). Its algebra of record, the SymPy
  identities and the line laws its tests verify, each with its
  verifier, is :ref:`theory-characteristic-origins`. Its five rungs
  exist. The first
  is the **walls**
  (:mod:`~orpheus.derivations.continuous.characteristic.walls`) and the
  **line part of the boundary closure**
  (:mod:`~orpheus.derivations.continuous.characteristic.closure`); the
  second is the **panel basis**
  (:mod:`~orpheus.derivations.continuous.characteristic.basis`) and the
  **transport along a line** on it
  (:mod:`~orpheus.derivations.continuous.characteristic.transport`); the
  third is **one group's transport block**, assembled over the lines
  (:mod:`~orpheus.derivations.continuous.characteristic.assembly`, with
  its gradings in
  :mod:`~orpheus.derivations.continuous.characteristic.grading`); the
  fourth is the **multigroup Galerkin system**
  (:mod:`~orpheus.derivations.continuous.characteristic.cross_sections`
  and :mod:`~orpheus.derivations.continuous.characteristic.system`),
  which answers the eigen and source questions in nodal coefficients;
  the fifth is **the door**
  (:mod:`~orpheus.derivations.continuous.characteristic.reference`), which
  poses that system from a specification and answers the eigenvalue and
  the flux integrals of its question (5a), and **the reading at a point**
  (:mod:`~orpheus.derivations.continuous.characteristic.reading`, on the
  weighted line sets of
  :mod:`~orpheus.derivations.continuous.characteristic.lines`), which
  answers its point values (5b).
- **The boundary resolvent has two parts, and both are built.** The
  closure of a reflecting boundary is :math:`P = P_0 + E\,(I - T)^{-1} X`
  on the boundary trace space. On every wall whose return is specular (a
  mirror, a partial mirror, vacuum as amplitude 0) or a periodic wrap,
  the returned path is a line congruent to the one that left, so
  :math:`T` is diagonal over lines: that is
  :class:`~orpheus.derivations.continuous.characteristic.closure.LinePeriod`.
  A diffuse wall couples every line to every other: the second part is a
  finite-rank update over the diffuse walls,
  :class:`~orpheus.derivations.continuous.characteristic.closure.WallCoupling`
  (:ref:`characteristic-resolvent`).
- **A wall is read from its law's two factors**, the deck
  (``geometry_map``) and the response (``response_kernel``), by one table
  (:ref:`characteristic-walls-factors`). It is keyed by its breakpoint
  index (:math:`0` or :math:`n`), never by a position, and carries a
  specular amplitude, a diffuse amplitude and a partner (the breakpoint
  at which the returned path re-enters). A wall is specular or diffuse,
  never both, and never returns more than it receives
  (:ref:`characteristic-walls`).
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
  times what traversal :math:`k` carries to its exit, plus any flux
  arriving there from outside the line part (a diffuse wall's re-entry).
  One rolled expression serves every rank; :math:`1 - \Pi` is formed by
  ``expm1``; a lossless trapped line (:math:`\Pi = 1`) carries exactly 0
  when nothing enters it and is refused when something does
  (:ref:`characteristic-closure-section`).
- **The panel basis** is discontinuous: nodal Lagrange polynomials of
  degree :math:`p` through each panel's Gauss–Legendre points, so a
  coefficient is a value. The panels are graded geometrically (panel ends
  at the depths :math:`w\rho^{j}`) toward every wall and interface and
  never toward a singular stratum. On the panel touching the centre or
  the axis the functions are polynomials in :math:`c^2`, because a smooth
  invariant flux is a smooth function of :math:`c^2` there (Schwarz's
  theorem) and an odd mode adds a :math:`b^{2m+2}\log b` term to every
  line integral (:ref:`characteristic-even-basis`). The mass matrix is
  exact in the chart's volume measure, whose density is the kernel's
  :meth:`Chart.measure_density <orpheus.geometry.chart.Chart.measure_density>`
  (:ref:`characteristic-panel-basis`).
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
  lie :math:`c/|P\Omega|` from a point of orbit coordinate :math:`c`, until
  each piece is no wider than its distance to them, so that every piece
  converges at a rate independent of the impact parameter
  (:ref:`characteristic-branch-grading`).
- **One group's transport block is a Galerkin assembly over lines**
  (:eq:`characteristic-galerkin-assembly`):
  :math:`K = \sum_L w_L \bigl(V_L + \sum_k A_k \otimes \mathrm{in}_k\bigr)`
  over each line's forward traversals, with the line weight :math:`w_L`
  the quadrature weight times the line domain's density over
  :math:`4\pi` (:eq:`geometry-line-domain`), so that :math:`K` is the
  scalar-flux operator of an isotropic emission and a closed homogeneous
  body satisfies :math:`K\Sigma_t\mathbf 1 = W\mathbf 1`. Its rows cover
  every panel and its columns the group's **emission support**
  (:ref:`characteristic-galerkin-assembly-section`).
- **The white walls couple through**
  :math:`R\,\alpha\,(I - T\alpha)^{-1}U^{\mathsf T}`
  (:eq:`characteristic-boundary-resolvent`), with reciprocity
  :math:`R = U D^{-1}`, :math:`D = \mathrm{diag}(A_w/4)` and :math:`A_w` the
  wall's area read from the kernel's one density. :math:`I - T\alpha` is
  formed from each injected current's **loss**, with one row replaced by
  the balance, never by the subtraction :math:`1 - T_{ww}`: conservation
  then holds to rounding from :math:`\Sigma_t = 1` down to
  :math:`10^{-12}`, where the subtraction missed by :math:`10^{-4}`, and a
  body that loses nothing is refused exactly (ERR-102,
  :ref:`characteristic-wall-coupling`).
- **Every grading of the line rule is derived from the group's optical
  scale** (the user's ruling of 2026-10-06): the impact parameter in the
  chord half-length, hp toward :math:`b = 0`, exponential at the rim, and
  hp toward each impact-panel top to the nearest feature of the transport
  there, the next radius, the turning slot's layer or the closure's pole
  (the grading law G, ERR-104, ERR-105,
  :ref:`characteristic-grading-law`); the grazing direction halved
  to :math:`\tau_{\min}/64`; the normal direction at :math:`2^k` over the
  body's normal optical depth. A fixed resolution passed every closed-body
  gate and missed the closed forms by up to :math:`5 \times 10^{-1}` near
  void and :math:`7.3 \times 10^{-1}` at :math:`\tau = 1000` (ERR-101,
  ERR-103, :ref:`characteristic-line-rule`). The cylinder's rule over
  :math:`(b, \theta)` is iterated: each polar angle carries its own impact
  rule, graded at its own projected speed :math:`\sin\theta`, which on the
  ``ABA`` cylinder holds 102 320 lines per block against the tensor rule's
  308 480 at the slowest speed, for the same k to
  :math:`8.3 \times 10^{-15}` (`[M]` 2026-10-08, #587,
  :ref:`characteristic-iterated-cylinder-rule`). The cylinder's rule held
  about :math:`10^6` pieces at :math:`\tau = 30` when its cost was
  measured, and that cost was the inner rule's padding, not the count:
  packing each line's live
  intervals brought the white block to 27.8 s per group, from about
  215 s (`[M]` 2026-10-07, #586, :ref:`characteristic-cylinder-cost`). Each
  radial line carries its node's exact level, the panel top and the
  half-chord there, so no half-chord is formed from a rounded :math:`b`
  (#590, ERR-106, :ref:`chart-and-chord-level`).
- **The multigroup system's unknown is the emission density** on each
  group's support, not the flux
  (:eq:`characteristic-pencil`): :math:`W_s q = (S + F/k)Kq + W_s q^{\rm ext}`,
  with the flux read as :math:`W_G\phi = Kq`. :math:`S` and :math:`F` act
  node by node, exactly, because the cross sections are constant on a
  region and no panel crosses one. The emission form is chosen for its
  adjoint: the transpose of the flux form :math:`\phi = \mathcal KE\phi`
  is the adjoint collision density :math:`E^\ast\phi^\dagger`, while the
  transpose of :math:`q = E\mathcal Kq` is the adjoint flux
  :math:`\phi^\dagger = \mathcal KE^\ast\phi^\dagger` itself
  (:eq:`characteristic-adjoint`; measured group ratios 2.147 against
  1.073, the user's ruling of 2026-10-07,
  :ref:`characteristic-galerkin-system`).
- **Each group's emission support is exact**: the regions where a
  scattering, an (n,2n) transfer or a fission enters the group, widened by
  the posed source regions. :math:`\Sigma_t > 0` drops a transparent
  region's emission (a balance off by :math:`2 \times 10^{-1}`), and the
  union over groups puts a source on lossless trapped lines behind a
  mirror (refused although well posed).
- **One operator, two splittings.** The k pencil
  :math:`(W_s - SK, FK)` answers the eigen questions and their adjoints;
  the source pencil :math:`(W_s, (S + F)K)` answers the fixed source and
  the detector (:eq:`characteristic-fixed-source`), so that its
  subcriticality check sees every secondary emission and refuses a body
  supercritical by its (n,2n) emission alone, which the k pencil would
  solve into a negative flux. In one group the Galerkin :math:`k` is a
  Rayleigh–Ritz lower bound, increasing on nested spaces
  (:eq:`characteristic-one-group-bound`).
- **The door** (:ref:`characteristic-door`).
  ``CharacteristicDerivation(specification, resolution)`` refuses at
  construction what the reference does not serve and solves on the first
  evaluation. A source enters by the section and a detector by the
  pullback, so a detector table's rate is :math:`RR^\dagger\Sigma_d =
  4\pi\Sigma_d`, the :math:`4\pi` derived from the angular measure
  (:eq:`characteristic-door-lifts`); a ``Response`` is answered as the
  forward flux of the group-transposed cross sections, which is the
  adjoint scalar flux :math:`R\psi^\dagger` (:eq:`characteristic-door-response`),
  and reciprocity reads
  :math:`\langle\Sigma_d, \phi(Q)\rangle = \langle Q, R\psi^\dagger\rangle/4\pi`.
  The fundamental flux is gauged so that the production the question
  declares (``Eigen.gauge``; by default what fission and the (n,2n)
  reaction emit, S\ :sub:`N`'s functional) is 1 over the body
  (:eq:`characteristic-door-gauge`); a ``Nearest`` answer reads its
  eigenvalue only, because a
  higher mode's net production can vanish (exactly, on a closed
  homogeneous body). A flux integral is the pairing
  :math:`(Wc_w)^{\mathsf T}\phi_h = \int w\,\phi_h\,\mathrm dV`, with no
  point (:eq:`characteristic-door-pairing`).
- **A point value is the transported emission**, :math:`\phi(x) =
  (\mathcal Kq)(x)`, the iterated Galerkin value, not the basis value
  :math:`\phi_h(x)`, which carries the projection
  :math:`\phi_h = P\mathcal Kq`. The directions at the point are the line
  domain's lines through it, built by the line rule's own constructors
  with the point as one more panel end and read at
  :math:`\pm\sqrt{c^2 - b^2}` from their closest approach, weighted by
  :math:`\mathrm d\Omega/4\pi` in line coordinates
  (:eq:`characteristic-quadrature`). One weighted set of lines,
  :class:`~orpheus.derivations.continuous.characteristic.lines.Lines`, has
  two roles: ``LineRule`` (the block) and ``PointRule`` (the point), so
  :math:`\int u_i\,\phi\,\mathrm dV = (Kq)_i` (C6); the white walls enter
  both through one fold of their currents. A point on a cylinder costs
  seconds, far above the verification spec's target (#591,
  :ref:`characteristic-reading`).
- **Evidence** `[M]` 2026-10-06 and 2026-10-07: for the walls and the
  closure, 175 gate rows in two files and a 34-arm mutation battery; for
  the basis and the transport, 328 rows in two more files, the traversal
  integrals' rows at L0; for the third rung, 46 new test functions (223
  cases) in three files, 37 earlier rows re-posed, and a 54-arm battery in
  which 52 arms redden their target rows and 2 are declared blind. The
  battery's honest run over the seven files without the ``slow`` rows:
  754 passed; for the fourth rung, 80 rows in 22 functions, 78 outside
  ``slow`` passed; for the door, 95 rows in 36 functions, 94 outside
  ``slow`` passed (`[M]` 2026-10-08), and a 39-arm battery in which every
  arm reddens at least one row; for the reading at a point, 175 rows in 40
  functions (15 ``slow``) and the level's 12 rows in 7, under a battery
  whose arms each redden their target rows
  (:ref:`characteristic-evidence`).
- **The S\ :sub:`N` rows read it** (:ref:`characteristic-sn-rows`).
  Since P1 step (d) the 13 S\ :sub:`N` cross-check rows that read the
  trajectory resolvent read this reference, the sphere and slab bodies at
  polynomial degree :math:`p = 5` and the cylinders at the door default
  :math:`p = 3`. Its error at each row is a ladder ESTIMATE (geometric,
  :math:`s/(1 - r)`), never a bound, and at every row it is at most the old
  family's and at least three orders below the S\ :sub:`N` residual, so
  every tolerance tightened or held (the A|B|A sphere's :math:`k` from
  :math:`4\times10^{-3}` to :math:`4\times10^{-5}`). Both sides read one
  datum from the posed problem, a partial wall's response factor; its pin
  is ``test_characteristic_walls.py``.


.. _characteristic-place:

Where this reference sits
=========================

The trajectory-resolvent family (:ref:`characteristic-origins-history`) was
seven oracle classes and fifteen solver entry points (the plan's
inventory, counted by ``git grep`` before the deletion), one per geometry
and boundary shape, each with its own chord oracle and its own closure;
the user ruled on 2026-10-05 that it be re-architected before it was
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
lines. During the migration the two packages coexisted under different
names, ``characteristic/`` and ``trajectory_resolvent/``, so neither was
a homonym of the other. P1 step (e2) (2026-10-10) deleted
``trajectory_resolvent/``; the SymPy derivations that grounded it moved
first, at step (e1a), to
:mod:`orpheus.derivations.continuous.characteristic.origins`.


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
line-diagonal block, ruled on 2026-10-06 ("One resolvent, two parts").
The ruling spelled it :math:`K_{\rm line} + U\,(I - T_w)^{-1} A\,U^{\mathsf T}`;
that spelling omits the wall areas. The built form is
:eq:`characteristic-boundary-resolvent`,
:math:`K = K_{\rm line} + R\,\alpha\,(I - T\alpha)^{-1}U^{\mathsf T}`, with
:math:`U` the escape functional from emission to the outgoing partial
current at each diffuse wall, :math:`T` the wall-to-wall transmission of
the line part, :math:`\alpha` the diffuse amplitudes and :math:`R` the
flux of a unit current entering at each wall, which reciprocity makes
:math:`U D^{-1}` with :math:`D = \mathrm{diag}(A_w/4)`
(:class:`~orpheus.derivations.continuous.characteristic.closure.WallCoupling`,
:ref:`characteristic-wall-coupling`). Dropping :math:`D` is not a
rounding: it misses closed-body conservation by 6.6, relative, on a
white sphere (`[M]` 2026-10-06, the test-architect's verification spec).


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
part diffusely is physically legitimate: the line part would read
``specular`` and the diffuse part ``diffuse``. The reference refuses it
all the same. ``Wall`` refuses a wall with both a specular and a diffuse
amplitude (``NotImplementedError``, *a wall is specular or diffuse*),
under a ``SCOPE-BOUNDARY[guard]`` tag, by the user's ruling of 2026-10-06
on the third rung's sketch: the S\ :sub:`N` realizer refuses a
``LawSum`` of two laws too, so production cannot pose such a wall and the
reference, which conforms to production, has no consumer for it. The
refusal is the reader's choice, not the formula's limit: the
test-architect's prototype served a wall of specular 0.4 and diffuse 0.6
through the unchanged coupling and it conserved to
:math:`1.2 \times 10^{-12}` (`[M]` 2026-10-06). Serving it later is one
arm of the reader, its gates, and one term in the white walls' balance,
which today assumes a diffuse wall returns nothing specularly
(:ref:`characteristic-wall-coupling`).

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
a ``Wall`` refuses an amplitude outside :math:`[0, 1]`, a specular and a
diffuse amplitude summing above 1 (a wall returning more than it
receives; `[M]` the main agent, 2026-10-06: ``Wall(1, 1.0, 1.0, 1)``
assembled a block whose smallest entry was :math:`-0.0118` before the
guard) and a wall that is both; a ``Walls``
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
       :math:`-0.2`), or a ``Wall`` built directly with a specular and a
       diffuse amplitude summing above 1
     - *not a physical wall*
     - a wall returns a fraction of what reaches it
   * - a ``Wall`` built directly with both a specular and a diffuse
       amplitude
     - *a wall is specular or diffuse*
     - a ``SCOPE-BOUNDARY``: production poses no such wall
       (:ref:`characteristic-walls`)

The two refusals a ``Wall`` makes of its own amplitudes (the sum above 1
and the mixed wall) are gated beside the coupling that reads them,
``tests/gates/derivations/test_characteristic_assembly.py::test_a_wall_returning_both_ways_or_more_than_it_receives_is_refused``,
each by its own fragment.

A line lying in an interface is refused by the reference (ruled
2026-10-06); that refusal is not a wall's, and it belongs to the
reading at a point (:ref:`characteristic-reading`). The point rule never
places such a line (its polar nodes are interior to
:math:`[0, \pi/2]`), and the angular flux at a cylinder's interface point
along the axis, :math:`w = 1`, is refused by the transport, *a point at
which the angular flux is read lies on a transit of its line*, not by a
fragment of its own (`[M]` 2026-10-08, the archivist's probe, the
two-region cylinder at :math:`c = 0.5`; :math:`w = 1 - 10^{-16}` is read).


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
resolvent :math:`(1 - \alpha e^{-\tau})^{-1}`, which the algebra of
record proves on the sphere (:eq:`peierls-greens-surface-fixed-point`,
:ref:`characteristic-origins-sphere`) and identifies with Sanchez 1986
Appendix Eq. (A4) :cite:`SanchezTTSP1986` and Pomraning–Siewert 1982
Eq. (14) :cite:`PomraningSiewert1982`. For :math:`m = 2` it is
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

.. _characteristic-arriving-flux:

The arriving flux in the cycle
------------------------------

A diffuse wall returns its outflow isotropically, so what it returns
along a line is not that line's own outflow times an amplitude: it is a
flux that arrives at the traversal's entry from outside the line part,
computed by the white walls' coupling (:ref:`characteristic-wall-coupling`).
:meth:`LinePeriod.inflow
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.inflow>`
takes it as the optional argument ``arriving``, :math:`s_k` per traversal,
with the shape of the outflow. It enters the cycle where the return of the
previous traversal enters, and it is **not** multiplied by the wall's
amplitude: the specular amplitude :math:`a` is what the wall returns along
the line, and a wall that returns diffusely has :math:`a = 0`.

With it, the flux entering traversal :math:`k` from outside the cycle's own
carry is

.. math::

   e_k \;=\; a_{k-1} B_{k-1} + s_k \;=\; r_{k-1} + s_k ,

and the fixed point of the cycle is :math:`\psi^{\rm in} = G\,\psi^{\rm in} + e`
with the same weighted cyclic shift :math:`G` as before. Its least
non-negative solution is the Neumann series grouped by whole periods,
exactly as in the closed form above with :math:`e` in place of the shifted
returns:

.. math::

   \psi^{\rm in}_k \;=\; \frac{e_k + \gamma_{k-1}\,e_{k-1}}{1 - \Pi}
   \quad (k \bmod 2),

the flux entering :math:`k` plus the flux that entered :math:`k - 1`,
carried once through it. With :math:`s = 0` it is the second form of
:eq:`characteristic-closure` term for term, so the labelled equation is the
special case and keeps its markers. In the code the change is small: the
returns are rolled once into the entering flux, the arriving flux is added
on the present traversals, and the once-around term rolls the entering
flux, not the returns rolled twice. ``arriving=None`` is the old
expression; a zero ``arriving`` is bitwise equal to it
(``test_a_zero_arriving_flux_is_no_arriving_flux_bitwise``), and the rank-1
and rank-2 closed forms with an arriving flux are
``test_an_arriving_flux_enters_the_cycle_at_its_traversal``, in
``tests/gates/derivations/test_characteristic_closure.py``.

On a lossless trapped line the refusal reads the entering flux, not the
outflow: a flux arriving on a line that never meets material and never
loses anything at a wall has no finite answer either, and
:class:`~orpheus.derivations.continuous.characteristic.closure.TrappedSource`
is raised (``test_an_arriving_flux_on_a_lossless_trapped_line_is_refused``).
The assembly passes the emission's outflows and the unit currents
injected at each diffuse wall as one stacked batch, the emission columns
first, so one call of ``inflow`` gives the line part, the escape, the
response and the transmission of every wall (:ref:`characteristic-galerkin-assembly-section`).

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
That is the line closure's only refusal; it is a ``ValueError``, raised
on the flux entering the line (the outflow the caller supplied, plus any
arriving flux), and it names the trapped line. The white walls' coupling
has its own, for a body that loses nothing behind walls that return
everything (:ref:`characteristic-wall-coupling`).

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
widths over the cosine), to 8 ulp relative. The depth is :math:`\tau_k`
of :eq:`characteristic-traversal-integrals`, and since the second rung
(``ff979520``) the three depth functions' 12 rows carry
``verifies("characteristic-traversal-integrals")`` at L0; it is also
what ``optical_depth`` implements under :eq:`characteristic-closure`.


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
region is one panel.
:func:`~orpheus.derivations.continuous.characteristic.grading.graded_ends`
computes them (the third rung moved it from the basis to the package's
module of gradings, which the slab's grazing rule also calls);
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
grading would add panels and nothing else. What the centre needs instead
is a basis that is even in :math:`c` (:ref:`characteristic-even-basis`). The gate's arms B1 (the centre
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

.. _characteristic-even-basis:

The even basis at a singular stratum
------------------------------------

On the panel whose lower end is a singular stratum of the chart, the
centre of a solid sphere or the axis of a solid cylinder, the functions
are polynomials in :math:`c^2`, not in :math:`c`
(:attr:`PanelBasis.even
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.even>`,
true on that panel and on no other; ruled by the user on 2026-10-06, the
plan's ledger, "on P1 step (b)'s third rung", Q1). The nodes are the same
Gauss–Legendre points :math:`c_m` of the panel :math:`[0, h]`, and the
functions are the Lagrange polynomials in :math:`s = (c/h)^2` through the
squared nodes :math:`s_m = (c_m/h)^2`:

.. math::

   u_{0,m}(c) \;=\; \prod_{l \ne m} \frac{s - s_l}{s_m - s_l},
   \qquad s = (c/h)^{2},

so they span :math:`1, c^{2}, \dots, c^{2p}`, every function of
:math:`c^{2}` of degree :math:`p` in :math:`c^2` and no odd power of
:math:`c`. A coefficient is still a value at a node, the product form is
the same one evaluated in the panel's own coordinate, and nothing
downstream of :meth:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis.values`
changed.

**Why the physics allows it.** The flux is invariant under the isotropy
of the stratum, :math:`O(3)` at the sphere's centre and the
:math:`O(2)` of rotations and reflections about the cylinder's axis
acting on the plane normal to it. A smooth function invariant under
:math:`O(d)` is a smooth function of :math:`|x|^2`: Schwarz's theorem
:cite:`Schwarz1975` for a compact group acting orthogonally, every smooth
invariant function is a smooth function of the generators of the
invariant polynomials, here :math:`|x|^2`. Restricted to a line through
the stratum it is an even function of :math:`c`, and a smooth even
function is a smooth function of :math:`c^2` (Whitney's theorem
:cite:`Whitney1943`). The flux near the stratum is therefore a smooth
function of :math:`c^2`, and the odd powers of :math:`c` a basis in
:math:`c` holds are modes the physics never excites.

**Why the numerics needs it.** The odd modes are not harmless: they make
every line integral through the centre panel non-smooth in the impact
parameter. Along a line of impact parameter :math:`b`, with
:math:`c(s)^2 = b^2 + s^2`, the integral of a power :math:`c^{n}` over the
panel is the Abel transform

.. math::

   \int c^{n}\,\mathrm{d}s \;=\; \int_{b}^{h} \frac{c^{n}\,c}{\sqrt{c^{2} - b^{2}}}\,\mathrm{d}c .

For even :math:`n = 2m` the integrand along the line is
:math:`(b^2 + s^2)^m`, a polynomial in :math:`s` and :math:`b^2`, and the
integral is a polynomial in :math:`b^2` times the chord's square root
:math:`\sqrt{h^2 - b^2}`, which the impact rule's substitution absorbs.
For odd :math:`n = 2m + 1` it carries a term in :math:`b^{2m+2}\log b`:
`[M]` 2026-10-07, SymPy's series of the integral at :math:`b \to 0` gives
:math:`-\tfrac12 b^2\log b` for :math:`n = 1` and
:math:`-\tfrac38 b^4\log b` for :math:`n = 3`. A logarithm at the end of
the impact interval holds Gauss–Legendre in :math:`b` to algebraic
convergence. `[M]` 2026-10-07, the archivist's probe on the built code
(the closed mirror sphere :math:`(0, 0.5, 1)`, :math:`p = 3`, one panel per
region, :math:`\Sigma_t = 0.7`: :math:`K\mathbf 1` from the impact piece
:math:`[0, 0.25]` at :math:`n` points against 128 points, the maximum
over the panels):

.. list-table::
   :header-rows: 1
   :widths: 30 17 17 17 19

   * - Centre panel
     - :math:`n = 4`
     - 6
     - 8
     - 16
   * - even (polynomials in :math:`c^2`)
     - :math:`1.9 \times 10^{-6}`
     - :math:`1.4 \times 10^{-10}`
     - :math:`4.6 \times 10^{-14}`
     - :math:`4.4 \times 10^{-16}`
   * - odd (Lagrange in :math:`c`, the second rung's basis)
     - :math:`8.1 \times 10^{-6}`
     - :math:`3.3 \times 10^{-7}`
     - :math:`3.5 \times 10^{-8}`
     - :math:`1.6 \times 10^{-10}`

The odd basis was measured by switching ``even`` off on every panel in
process. At the operator level the same probe gives closed-body
conservation of the white sphere at 24 impact points per piece of
:math:`1.3 \times 10^{-15}` with the even panel and
:math:`2.6 \times 10^{-13}` with the odd one.

**Why not grade the impact rule toward** :math:`b = 0` **instead.** The
premises of the third rung measured that alternative: grading the
:math:`b` rule geometrically toward 0 took the closed sphere from
:math:`5.3 \times 10^{-8}` to :math:`4.4 \times 10^{-11}` at 16 points
(`[M]` 2026-10-06, the main agent's
``scratch/characteristic_architecture/p1_step_b3/ladder3.py``). It treats
the symptom: the logarithm is still in every integrand and is chased
with points. The even basis removes it at its cause, the basis spanning
functions the physics excludes, and leaves the centre piece smooth in
:math:`b^2` (the impact rule then takes plain Gauss–Legendre on the lower
half of that piece, :ref:`characteristic-line-rule`).

**What follows from it.** The even panel's polynomials have degree
:math:`2p` in :math:`c`, so the mass rule takes :math:`2p + 2` points on
every panel (below). On the even panel a line's integrand is a polynomial
in arc length, so the orbit coordinate's branch points do not affect it
(:ref:`characteristic-branch-grading`): the gates of the turning grading
therefore sit on a hollow body with a cavity of radius :math:`10^{-5}`,
where the line turns in an ordinary panel, and the small-radius rows of
ERR-099 read an even source on the even panel. The gates are ``test_the_even_panel_is_the_one_touching_a_singular_stratum``
(the panel by hand from the body, not from ``singular_strata``),
``test_the_even_panel_spans_the_polynomials_in_c_squared`` (an even
polynomial of degree :math:`2p` reproduced to 64 ulp, :math:`c^{2p+1}`
missed by more than :math:`10^{-3}`), and
``test_the_even_panels_mass_is_the_volume_integral_of_its_even_products``,
in ``tests/gates/derivations/test_characteristic_basis.py``; at the
operator level, ``test_the_centre_impact_piece_converges_geometrically``
and ``test_the_closed_sphere_conserves_to_rounding_at_high_resolution`` in
``test_characteristic_assembly.py``.

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

which the kernel computes,
:meth:`Chart.measure_density <orpheus.geometry.chart.Chart.measure_density>`
(:eq:`geometry-measure-density`), so
:math:`W_{ij} = \int u_i(r)\,u_j(r)\,\kappa\,d\,r^{d-1}\,\mathrm{d}r`. On
the even panel the integrand is a polynomial of degree at most
:math:`4p + d - 1` in :math:`r` (on every other panel :math:`2p + d - 1`),
and Gauss–Legendre with :math:`2p + 2` points, one rule on every panel, is
exact to degree :math:`4p + 3`, so each entry is integrated exactly up to
rounding. :math:`W` is block-diagonal, one :math:`(p + 1) \times (p + 1)`
block per panel, because the supports are disjoint. The second rung's rule
of :math:`p + 2` points was exact only for the odd basis; on the even
panel it reddens the even mass gate on the sphere at :math:`p = 1, 3, 5`
and on the cylinder at :math:`p = 3, 5`, and is blind on the cylinder at
:math:`p = 1`, where degree 5 is exact at 3 points (declared in the gate).

**Why the density is the kernel's.** On the second rung the user ruled
(the plan's ledger, "on P1 step (b)'s second rung", Q4) that the basis
derive the density itself from ``measure_constant`` and
``measure_coordinate`` rather than add a verb to the stable kernel. The
third rung brought a second consumer, the white wall's area
(:ref:`characteristic-wall-coupling`), and the user ruled on 2026-10-06
(Q4 of the third rung) that the density move to the kernel as the
derivative of the measure. The basis's own ``PanelBasis.volume_density``
retired onto it, and the route is gated:
``test_characteristic_assembly.py::test_the_mass_and_the_wall_area_read_the_one_density``
doubles ``MeasureCoordinate.derivative`` in process and requires the mass
to double and the wall's response to halve bit for bit, so no second
spelling survives on either side. ``test_each_panels_mass_is_its_chart_measure``
sums each panel's block, which with a partition of unity is the panel's
volume, against :meth:`Chart.measure <orpheus.geometry.chart.Chart.measure>`
and against the closed-form volume in mpmath. The first leg shares the
upstream constant :math:`\kappa` with the code and is not independent; the
second leg, and ``test_the_volume_density_is_the_charts_by_hand`` (the
three densities written by hand, read now through the basis's chart),
are. Arm B7 of the second rung's battery (the factor :math:`d` dropped)
reddened 44 rows.

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
optical width exceeds 2 (the constant ``_THIN``; both are
:func:`~orpheus.derivations.continuous.characteristic.grading.exponential_ends`
since the third rung moved the gradings into one module). A half of optical
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

On every panel but the even one (:ref:`characteristic-even-basis`) a
basis function is a general polynomial in :math:`c`, odd powers
included, so the integrands are analytic in :math:`t` everywhere except at
the branch points (the even powers of :math:`c` are polynomials in
:math:`t`). On the even panel every function is a polynomial in
:math:`c^2`, hence in :math:`t`, and the branch points do not touch it;
the grading below still runs on its slots, which costs pieces and changes
no value. On a slot :math:`c` is monotone, because the closest approach
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
:math:`\ell\,2^{-k}`, :math:`k = 1, \dots, K`, from the near end, with
:math:`K = \lceil \log_2(\ell/D) \rceil` the number of halvings after
which a piece is no wider than its distance to the branch points
(``TraversalRule._branch_edges``, through
:func:`~orpheus.derivations.continuous.characteristic.grading.halvings`,
the hp law every grading of the package shares; :math:`K` is at most 52,
``np.finfo(float).nmant``, below which a depth is under one ulp of the
slot). A piece :math:`[\delta/2, \delta]` (depths from the near end) has
the width :math:`\delta/2` and lies at least :math:`\delta/2` from the
branch points, whatever :math:`b`; the innermost piece
:math:`[0, \delta_{\min}]` has :math:`D/2 < \delta_{\min} \le D`. The ellipse parameter of every piece
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
     - 4.61
     - :math:`1.9 \times 10^{-11}` / :math:`5.1 \times 10^{-10}`
   * - quartering
     - 3.00
     - 4.61
     - :math:`2.1 \times 10^{-7}` / :math:`5.1 \times 10^{-10}`

The last column is the bound's decay factor at 7 points, before its
constant: halving gains about four orders of magnitude on
every graded piece. `[M]` 2026-10-07, the archivist's re-evaluation of the
innermost column under the law as it now stands: the worst position of
the branch points in the closed half-plane beyond the near end, for
:math:`\delta_{\min}/D` over :math:`(1/2, 1]`, gives 4.61 for either ratio.

.. dropdown:: First got wrong: the halving one short of the law
   :color: muted

   Until the third rung's last review the halving kept a depth only while
   it exceeded :math:`D`, so the innermost piece reached :math:`2D` and its
   ellipse parameter fell to 2.89 (2.08 quartering). The elegance re-review
   of 2026-10-07 found the hp law spelled three ways in the package and
   this one a halving short, and every spelling became
   ``grading.halvings``.

The bound's constant :math:`M` is small on the innermost piece, where
:math:`c` itself is of the order of :math:`c_{\rm near}`, which is why
the measured errors below sit far under the factor.

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
controls. The table was measured on the second rung's code, on solid
bodies whose lines turn in the centre panel. On the third rung that panel
became even and the integrand there a polynomial in arc length, so the
gate's rows moved to a hollow body with a cavity of radius
:math:`10^{-5}`, where the line turns in an ordinary panel; which of
them the ratio-1/4 arm reddens there is in the third rung's battery
(``scratch/characteristic_architecture/p1_step_b3/gates/battery/``). Each line keeps its live pieces in chord order (slot by slot,
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
(:ref:`characteristic-galerkin-assembly-section`).

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

.. _characteristic-galerkin-assembly-section:

The Galerkin assembly over lines
================================

The third rung assembles one group's **transport block**: the operator
that takes an isotropic emission density on the basis to the scalar flux
it produces through the body and its walls, tested against the basis.
:meth:`LineRule.transport
<orpheus.derivations.continuous.characteristic.assembly.LineRule.transport>`
returns it as a
:class:`~orpheus.derivations.continuous.characteristic.assembly.GroupTransport`,
whose :attr:`~orpheus.derivations.continuous.characteristic.assembly.GroupTransport.block`
is the line part plus the white walls' update. The design is the plan's
"P1 step (b), third rung: API sketch" (items 5 and 6), ruled 2026-10-06;
the assembly by Galerkin over lines, rather than collocation at points
and directions, was ruled the same day after its measurement
(``scratch/characteristic_architecture/p1_assembly/report.md``: a
:math:`k` error 28 to 175 times smaller at equal :math:`n`, an observed
order of about 5.5 against 3.6).

The bilinear form, derived
--------------------------

Let :math:`q` be an isotropic emission density, so that a volume element
emits :math:`q/4\pi` per steradian, and let :math:`\phi[q]` be the scalar
flux it produces. The transport block is

.. math::

   K_{ij} \;=\; \int u_i(x)\,\phi[u_j](x)\,\mathrm{d}V
   \;=\; \int_{S^2}\!\mathrm{d}\Omega \int u_i(x)\,\psi_j(x, \Omega)\,\mathrm{d}V .

For one direction :math:`\Omega`, write a point as the foot :math:`p` of
the line through it in the plane normal to :math:`\Omega` plus an arc
length along the line, so that :math:`\mathrm{d}V = \mathrm{d}A_\perp\,\mathrm{d}s`.
The double integral over directions and positions is then an integral
over oriented lines against the invariant measure
:math:`\mathrm{d}L = \mathrm{d}A_\perp\,\mathrm{d}\Omega`
(:eq:`geometry-measure-on-lines`) of an integral along each line:

.. math::

   K_{ij} \;=\; \int \mathrm{d}L \int_{L} u_i\bigl(c(s)\bigr)\,\psi_j(s)\,\mathrm{d}s .

Along the line, :math:`\psi_j` is :math:`1/4\pi` times what the traversal
rule computes from the unit source :math:`u_j`, which carries no
:math:`1/4\pi` (:ref:`characteristic-transport`): the vacuum part, the
Volterra block :math:`V_L` of the line's forward transits
(:ref:`characteristic-volterra`), plus, on each forward traversal
:math:`k`, the inflow :math:`\mathrm{in}_k` the closure returns from the
outflows :math:`B` (:eq:`characteristic-closure`), attenuated from the
entry and paired with :math:`u_i`, which is the entry response
:math:`A_k` (:eq:`characteristic-traversal-integrals`). The integrand is
invariant under the chart's group, because every basis function is a
function of the orbit coordinate, so by :eq:`geometry-line-domain` the
integral over lines is an integral over the line domain's box with the
density :math:`\varrho`, and a quadrature with nodes :math:`q_L` and
weights :math:`W_L` on the box gives:

.. math::
   :label: characteristic-galerkin-assembly

   K \;=\; \underbrace{\sum_{L} w_L \Bigl(V_L + \sum_{k\ \text{forward}} A_k \otimes \mathrm{in}_k\Bigr)}_{K_{\rm line}}
   \;+\; K_{\rm wall},
   \qquad
   w_L \;=\; \frac{W_L\,\varrho(q_L)}{4\pi},

with :math:`K_{\rm wall}` the white walls' update of
:eq:`characteristic-boundary-resolvent` (zero with no diffuse wall). The
rows of :math:`K` run over every basis function and its columns over the
emission support (below).

.. implements:: characteristic-galerkin-assembly
   :by: orpheus.derivations.continuous.characteristic.assembly.LineRule.transport

   **Implemented by** ``LineRule.transport``, which accumulates
   :math:`K_{\rm line}` chunk by chunk from each chunk's
   ``TraversalRule``; ``LineRule.of`` places the lines on the line domain
   and forms :math:`w_L`; ``GroupTransport.block`` adds the walls' update.

.. implements:: characteristic-galerkin-assembly
   :by: orpheus.derivations.continuous.characteristic.assembly.LineRule.of

.. implements:: characteristic-galerkin-assembly
   :by: orpheus.derivations.continuous.characteristic.assembly.GroupTransport.block

**Only the forward traversals.** A point of the line domain is an orbit of
oriented lines, and the integrand is evaluated on its representative
read in its own direction: its transits read forward, and the inflows of
the traversals that read them forward. The reversed traversals of the
period carry the line's cycle and belong to the opposite line. On the
slab the opposite line is another point of the box (:math:`-\mu`), which
the rule holds on its own; on the sphere and the cylinder it is in the
same orbit, and the density :math:`\varrho` already counts it. The battery
arm that sums the reversed traversals too (A5) reddens 26 rows.

**The normalisation, and what it makes checkable.** With
:math:`w_L = W_L\varrho/4\pi`, :math:`K` is the scalar-flux operator of an
isotropic emission. On a closed homogeneous body (every wall a mirror, or
white with :math:`\alpha = 1`, or a periodic wrap), :math:`\psi = 1` solves
the transport problem with the emission :math:`\Sigma_t`, so
:math:`\phi = 1` and

.. math::

   K\,\Sigma_t\mathbf 1 \;=\; W\mathbf 1 ,

the block applied to the nodal :math:`\Sigma_t` (piecewise constant, so in
the basis exactly) equals the volume of each basis function, the mass
matrix applied to ones. The right side reads only the volume measure, the
left side the lines, the closure and the coupling, so the identity is the
assembly's **closed-body conservation** row. `[M]` 2026-10-07, the
archivist's probe on the built code, homogeneous bodies of radius or width
1.3 at :math:`\Sigma_t = 0.7` behind mirrors: :math:`4.0 \times 10^{-15}`
on the sphere at 16 points per piece, :math:`9.8 \times 10^{-16}` on the
slab at 8. On the three-region data of the gates
(:math:`\Sigma_t = (0.6, 1.3, 0.45)`, breakpoints :math:`(0, 0.5, 1.5, 2)`):
the sphere behind a mirror :math:`1.6 \times 10^{-15}` and behind a white
wall :math:`1.5 \times 10^{-15}` at 16 points (:math:`7.9 \times 10^{-13}`
at 8); the hollow sphere :math:`(0.4, 0.5, 1.5, 2)` white inside and out
:math:`1.7 \times 10^{-15}`, mirror inside and white out
:math:`1.8 \times 10^{-15}`; the slab between two white walls
:math:`1.9 \times 10^{-15}` at 4 points.

**What conservation cannot see.** A mirror returns each line's outflow
into the same line and a periodic wrap carries it into a congruent one,
so on a closed body the identity holds **line by line**: every line
conserves on its own, whatever weight the rule gives it. Conservation
therefore cannot see the rule over lines at all, only the transport along
each line and the coupling of the white walls. The test-architect measured
it twice: the slab's cosine rule passed 33 of 33 conservation runs to
:math:`2.5 \times 10^{-15}`, plain 4-point rules included; and the
cylinder's axial-cosine rule read the same under a mirror at 8 points as
at 16 (:ref:`chart-and-chord-line-domain`). The rule over lines is gated
by the closed forms of the walls instead (the escape probability and the
wall transmission, :ref:`characteristic-wall-coupling`), which read every
line with its weight. That every line-measure defect of ERR-101 passed
the conservation rows is this blindness.

**Symmetry is a foundation row, declared blind.** By reciprocity the
block is symmetric where its rows and columns coincide, and the third
rung's sketch proposed the symmetry defect as the alarm for an
under-integrated rule (the first spec's C11). It cannot be one: each
line's quadrature is symmetric under reversing the line (the inbound and
outbound halves of a radial chord mirror each other, and the slab's
cosine rule is symmetric in sign), so each line's block is symmetric
whatever the rule's accuracy. `[M]` 2026-10-06, the main agent: the block
symmetric to :math:`2.5 \times 10^{-16}` at every resolution tried,
including the sphere at 8 impact points, where conservation missed by
:math:`5.2 \times 10^{-5}`, and the slab with a plain cosine rule, which
missed by :math:`1.6 \times 10^{-3}`. Its only teeth are a mismatch
between the outer and the inner arc-length rules
(:math:`1.5 \times 10^{-6}` at 4 outer points against 12 inner), which
``test_the_symmetry_row_sees_the_arc_length_rules`` keeps. The user ruled
(Q3 of the third rung) that conservation is the alarm and symmetry a
declared-blind foundation row.

The emission support
--------------------

The columns of :math:`K` are the basis functions of the regions that
emit in the group: the regions of a mask the caller passes
(``transport(..., support=...)``), which the multigroup system sets to the
group's exact emission support (:ref:`characteristic-galerkin-system`),
and by default, for a block assembled on its own, those with
:math:`\Sigma_t > 0`;
``GroupTransport.support`` holds their indices. The rows run over every panel, so the flux is read
everywhere, and the block is rectangular, :math:`N \times M`.

**Why the columns stop at the emission.** A basis function on a panel
where nothing is emitted into the group is a source in a region that never
emits in it. Under a mirror it is
worse than useless: a line that stays in an outer void shell behind a
mirror is a lossless trapped line, and a source on it has no finite flux
(:ref:`characteristic-closure-section`), so the full block of such a body
does not exist, and ``inflow`` refuses it. `[M]` 2026-10-06, the
test-architect: on a sphere with a void layer :math:`(2.0, 2.6)` under a
mirror of amplitude 1 the full block raises
:class:`~orpheus.derivations.continuous.characteristic.closure.TrappedSource`
(``test_the_full_block_of_a_void_layer_under_a_full_mirror_is_refused``),
while the block on the emission support equals the block of the body
without the layer to :math:`5.4 \times 10^{-17}`
(``test_a_void_outer_layer_is_invisible_on_the_emission_support``). The
main agent took this decision while building the third rung (not a user
ruling); the fourth rung's ruling made the support exact per group. The rows stay full because the white walls' response
:math:`R` must be read on every panel; its reciprocal partner :math:`U`
lives on the support only, which is why :math:`R` is computed directly and
the reciprocity :math:`R = U D^{-1}` is gated rather than assumed (the
elegance review withdrew its objection to the two tallies on that
ground).

A region void in the group that does receive emission into it (another
group scatters or fissions into it there) gets a column; behind a mirror
the lines trapped in it then carry a source, and the block is refused,
because that flux is infinite. Why the support is neither each group's
:math:`\Sigma_t > 0` nor one set shared by every group is in
:ref:`characteristic-galerkin-system`.

The tallies, in one pass
------------------------

:meth:`~orpheus.derivations.continuous.characteristic.assembly.LineRule.transport`
walks the rule's chunks of lines once. For each chunk it builds the
chunk's ``TraversalRule`` and calls
:meth:`LinePeriod.inflow
<orpheus.derivations.continuous.characteristic.closure.LinePeriod.inflow>`
once on a stacked batch of sources: the outflows of the :math:`M` emission
functions, then, for each of the :math:`W` diffuse walls, a unit current
entering at that wall as an arriving flux of :math:`1/D_w` on every
traversal that enters there (:ref:`characteristic-arriving-flux`). From
that one inflow it accumulates, with the line weights,

- the line part, the Volterra block on the support columns plus
  :math:`\sum_k A_k \otimes \mathrm{in}_k` over the forward traversals,
  for the emission columns;
- the walls' response :math:`R`, the same sum for the wall columns;
- the escape :math:`U`, what leaves each forward traversal at a diffuse
  wall, :math:`e^{-\tau_k}\mathrm{in}_k + B_k`, for the emission columns;
- for the wall columns, over **every** traversal of the period (where the
  balance closes line by line): the current injected, the current
  absorbed (:math:`\mathrm{in}\,(1 - e^{-\tau})`, formed by ``expm1``),
  the current leaked at a wall that does not return it, and the current
  reaching each diffuse wall.

The transmission :math:`T` and the loss :math:`\ell` are the last two
divided by the first: ratios of one tally, so the slab's double count of
each orbit over every traversal (its cosine rule is symmetric in sign)
cancels, and so does the rule's quadrature error in the injected current.

Chunk invariance is not bitwise: the chunks re-order the sum over lines.
`[M]` 2026-10-06, the test-architect: the slab's block moves by
:math:`2.9 \times 10^{-15}`, 12.9 ulp of its largest entry, between chunk
sizes, so ``test_the_block_does_not_depend_on_the_piece_budget`` holds the
block to 64 ulp across the default budget, a quarter of it, a budget of
one line per chunk and an unbounded one, and checks that every line lands
in exactly one chunk (the first spec's 1e-15 was refuted).


.. _characteristic-wall-coupling:

The white walls' coupling
=========================

A diffuse wall of amplitude :math:`\alpha_w` re-emits isotropically the
fraction :math:`\alpha_w` of the current that reaches it, so what it
returns depends on every line that reaches it and feeds every line that
leaves it. That couples all the lines through a few numbers, the
currents at the diffuse walls, and the second part of the boundary
resolvent is a finite-rank update over them,
:class:`~orpheus.derivations.continuous.characteristic.closure.WallCoupling`.

The four quantities
-------------------

- The **escape** :math:`U`, :math:`M \times W`: the partial current that
  leaves through wall :math:`w` when the emission is the basis function
  :math:`u_i`. The outgoing current through a surface element
  :math:`\mathrm{d}A` with normal :math:`n` is
  :math:`\int_{\Omega\cdot n > 0} (\Omega\cdot n)\,\psi\,\mathrm{d}\Omega\,\mathrm{d}A`,
  and the lines through :math:`\mathrm{d}A` in the direction :math:`\Omega`
  fill, in the plane normal to :math:`\Omega`, the area
  :math:`\mathrm{d}A_\perp = |\Omega\cdot n|\,\mathrm{d}A`. So the current
  leaving a wall is the invariant measure on lines integrated over the
  lines that exit there, with the intensity they carry out, and the line
  weights of the assembly serve unchanged:

  .. math::

     U_{iw} \;=\; \sum_{L} w_L \sum_{k\ \text{forward, exiting at } w}
     \bigl(e^{-\tau_k}\,\mathrm{in}_k + B_k\bigr)_i .

- The **injection** :math:`1/D_w`. A current :math:`J` entering a wall of
  area :math:`A_w` isotropically has the intensity :math:`\psi` constant
  over the inward hemisphere, and
  :math:`J = A_w\,\psi\int_{\Omega\cdot n < 0}|\Omega\cdot n|\,\mathrm{d}\Omega = \pi A_w\,\psi`,
  so :math:`\psi = J/(\pi A_w)`. The lines carry :math:`4\pi` times an
  intensity (their weights hold the :math:`1/4\pi`), so a unit current
  enters each line as :math:`4\pi/(\pi A_w) = 1/D_w` with

  .. math::

     D_w \;=\; \frac{A_w}{4},

  and :math:`A_w` is the area of the level set at the wall, the kernel's
  :meth:`Chart.measure_density <orpheus.geometry.chart.Chart.measure_density>`
  (:eq:`geometry-measure-density`): :math:`4\pi R^2` on a sphere,
  :math:`2\pi R` per unit height on a cylinder, 1 per unit area on a slab.
  It is read at that one place.
- The **response** :math:`R`, :math:`N \times W`: the flux moments
  :math:`\int u_i\,\phi\,\mathrm{d}V` of a unit current entering at
  :math:`w`, the line part's sum :math:`\sum_k A_k \otimes \mathrm{in}_k`
  for the injection.
- The **transmission** :math:`T`, :math:`W \times W`: the fraction of a
  current entering at :math:`w` that leaves at :math:`w'` through the line
  part (specular returns on the way included), and the **loss**
  :math:`\ell_w`, the fraction absorbed or leaked through a wall that
  returns nothing. Every entering current goes one of the three ways:

  .. math::

     \sum_{w'} T_{w'w} + \ell_w \;=\; 1,

  per line in exact arithmetic, and to 8 ulp in the gate
  ``test_each_injected_current_is_transmitted_or_lost``.

The update, derived
-------------------

Let :math:`j_w` be the current wall :math:`w` returns. What reaches wall
:math:`w` is what the emission sends there, :math:`(U^{\mathsf T} q)_w`,
plus what the returned currents send there through the line part,
:math:`(Tj)_w`, and the wall returns the fraction :math:`\alpha_w` of it:

.. math::

   j \;=\; \alpha\,\bigl(U^{\mathsf T} q + T j\bigr)
   \quad\Longrightarrow\quad
   (I - \alpha T)\,j \;=\; \alpha\,U^{\mathsf T} q
   \quad\Longrightarrow\quad
   j \;=\; \alpha\,(I - T\alpha)^{-1}U^{\mathsf T} q ,

the last step by the push-through identity
:math:`(I - \alpha T)^{-1}\alpha = \alpha\,(I - T\alpha)^{-1}`, which
needs no inverse of :math:`\alpha` (a wall of amplitude 0 is admitted).
The returned currents produce the flux moments :math:`Rj`, so the walls
add to the block

.. math::
   :label: characteristic-boundary-resolvent

   K_{\rm wall} \;=\; R\,\alpha\,(I - T\alpha)^{-1}\,U^{\mathsf T},
   \qquad
   R\big|_{\rm support} \;=\; U\,D^{-1},
   \qquad
   D \;=\; \mathrm{diag}\bigl(A_w/4\bigr),

the diffuse part of :math:`P = P_0 + E\,(I - T)^{-1}X`
(:ref:`characteristic-resolvent`), :math:`E` and :math:`X` restricted to
the diffuse walls' currents being :math:`R` and :math:`U^{\mathsf T}`. The
code splits the update at :math:`j`: the currents
:math:`\alpha(I - T\alpha)^{-1}U^{\mathsf T}`, ``(W, M)``, are one value,
and the walls are folded onto the emission by one map,
``on_emission(e, w) = e + w @ currents``, read by the block with
:math:`w = R` and by the reading at a point with :math:`w = r(x)`, the
point's reading of a unit current entering each wall
(:ref:`characteristic-reading`). The block is
``on_emission(line, response)``, bit for bit
``line + response @ currents``
(``test_the_block_folds_its_walls_through_the_currents``).

.. implements:: characteristic-boundary-resolvent
   :by: orpheus.derivations.continuous.characteristic.closure.WallCoupling.currents

   **Implemented by** ``WallCoupling.currents``, the solve for
   :math:`j = \alpha(I - T\alpha)^{-1}U^{\mathsf T}`, with
   ``WallCoupling.returning`` forming :math:`I - T\alpha` from the loss;
   ``WallCoupling.on_emission``, the fold :math:`K_{\rm line} + R\,j` that
   ``GroupTransport.block`` reads; ``LineRule.transport``, which tallies
   :math:`U`, :math:`R`, :math:`T` and :math:`\ell`; and ``DiffuseWalls``,
   which reads :math:`D` from the chart's density.

.. implements:: characteristic-boundary-resolvent
   :by: orpheus.derivations.continuous.characteristic.closure.WallCoupling.on_emission

.. implements:: characteristic-boundary-resolvent
   :by: orpheus.derivations.continuous.characteristic.closure.DiffuseWalls.of

.. implements:: characteristic-boundary-resolvent
   :by: orpheus.derivations.continuous.characteristic.closure.WallCoupling.returning

.. implements:: characteristic-boundary-resolvent
   :by: orpheus.derivations.continuous.characteristic.assembly.LineRule.transport

**Reciprocity.** The flux at :math:`u_i` from a unit isotropic current
entering at :math:`w` and the current leaving at :math:`w` from the source
:math:`u_i` are the forward and the adjoint readings of one set of lines,
so :math:`R_{iw} = U_{iw}/D_w` on the emission support; and the
wall-to-wall transmission satisfies the surface reciprocity
:math:`A_w T_{w'w} = A_{w'} T_{ww'}`, that is :math:`D^{-1}T` symmetric.
`[M]` 2026-10-07, the archivist's probe: :math:`R = UD^{-1}` on the white
sphere's support to :math:`1.1 \times 10^{-16}` relative, and
:math:`D = \pi R^2` for a sphere of radius 2 exactly. The gates
``test_the_walls_response_is_the_escape_over_the_quarter_area`` (64 ulp)
and ``test_the_transmission_is_reciprocal_in_the_wall_areas`` (8 ulp)
hold both, and the void bodies of
``test_a_void_body_transmits_its_geometric_fractions`` pin the scale that
a reciprocity row cannot see: a hollow void sphere of radii 0.4 and 2
transmits :math:`(0.4/2)^2` of the outer wall's current to the cavity and
the rest back to itself, a void slab everything face to face.

**The correction to the ruled spelling.** The ruling of 2026-10-06 ("One
resolvent, two parts") wrote the update
:math:`U\,(I - T_w)^{-1} A\,U^{\mathsf T}`. With :math:`R = UD^{-1}` the
built form is :math:`U D^{-1}\alpha\,(I - T\alpha)^{-1}U^{\mathsf T}`: the
ruled form omits :math:`D`, an error of the factor :math:`4/A_w`, which is
not 1 on any body (the premises' measurement found it before the sketch;
dropping it misses conservation by 6.6, relative, on a white sphere, the
test-architect's arm, and battery arm A1 reddens 58 rows).

**Why** :math:`R` **stays per nominal unit current.** The line rule
integrates the injected current with a quadrature error, so the current
it actually injects misses 1 slightly: `[M]` 2026-10-06, the elegance
review, :math:`2.7 \times 10^{-14}` on the cylinder at 8 points and
:math:`1.1 \times 10^{-7}` at 4. :math:`T` and :math:`\ell` are fractions
of the injected tally and are unaffected. :math:`R` is the response to
the nominal injection :math:`1/D_w`, so the wall's area enters the block
there and nowhere else. Dividing :math:`R` by its own injected tally
instead was tried in the last review round and reverted: the area then
cancels out of the block entirely, and the route gate on the one density
(``test_the_mass_and_the_wall_area_read_the_one_density``) reddened. The
two normalisations differ by the rule's quadrature error, which is what
the stated asymmetry costs.

The loss, not the difference
----------------------------

Near void the walls exchange almost every neutron: :math:`T_{ww}` is
within the absorption of 1, and :math:`I - T\alpha` is nearly singular in
its total-current mode, whose eigenvalue is about what the body loses.
Formed by subtraction, its diagonal :math:`1 - \alpha_w T_{ww}` keeps only
the digits of :math:`T_{ww}` that differ from 1, and the rounding of the
tally is amplified by the inverse of the absorption. ``WallCoupling``
never subtracts:

- **the diagonal is formed from the loss**: by the balance,
  :math:`1 - \alpha_w T_{ww} = (1 - \alpha_w) + \alpha_w\bigl(\ell_w + \sum_{w' \ne w} T_{w'w}\bigr)`,
  a sum of non-negative terms, with :math:`\ell_w` itself a sum of
  absorbed fractions formed by ``expm1`` and leaked fractions
  (:attr:`~orpheus.derivations.continuous.characteristic.closure.WallCoupling.returning`);
- **one row is replaced by the balance**: the column sums of
  :math:`I - T\alpha` are :math:`1 - \alpha_w\sum_{w'}T_{w'w} = (1 - \alpha_w) + \alpha_w\ell_w`,
  known without cancellation, so the solve replaces the last row of the
  system by the sum of all its rows (and the last row of the right side by
  the sum of its rows). That is an exact row operation, so the solution is
  the same in exact arithmetic; in floating point the total current, the
  nearly singular mode, is read from the balance and not from rounded
  entries.

The diagonal alone was not enough: with it but no balance row, the two
bodies with two white walls still missed conservation by
:math:`7.7 \times 10^{-6}` and :math:`2.2 \times 10^{-5}` at
:math:`\Sigma_t = 10^{-12}` through the rounded off-diagonal entries
(`[M]` 2026-10-06, the main agent). `[M]` 2026-10-07, the archivist's
probe on the built code, closed-body conservation at 16 points per piece
with :math:`\Sigma_t` equal in both regions, the built coupling against
the same tallies solved through the subtraction
:math:`I - T\alpha` (its stored diagonal kept):

.. list-table::
   :header-rows: 1
   :widths: 34 16 16 16 18

   * - Body, :math:`\Sigma_t`
     - 1
     - :math:`10^{-6}`
     - :math:`10^{-9}`
     - :math:`10^{-12}`
   * - white sphere :math:`(0, 0.6, 1)`, built
     - :math:`1.0 \times 10^{-15}`
     - :math:`9.9 \times 10^{-16}`
     - :math:`1.0 \times 10^{-15}`
     - :math:`9.9 \times 10^{-16}`
   * - the same, by subtraction
     - :math:`1.0 \times 10^{-15}`
     - :math:`1.1 \times 10^{-11}`
     - :math:`2.5 \times 10^{-7}`
     - :math:`1.3 \times 10^{-4}`
   * - hollow sphere :math:`(0.3, 0.6, 1)`, white and white, built
     - :math:`9.7 \times 10^{-16}`
     - :math:`9.9 \times 10^{-16}`
     - :math:`9.7 \times 10^{-16}`
     - :math:`9.9 \times 10^{-16}`
   * - the same, by subtraction
     - :math:`9.7 \times 10^{-16}`
     - :math:`1.2 \times 10^{-11}`
     - :math:`1.8 \times 10^{-7}`
     - :math:`1.3 \times 10^{-4}`
   * - slab :math:`(0, 0.6, 1)`, white and white, built
     - :math:`1.1 \times 10^{-15}`
     - :math:`1.2 \times 10^{-15}`
     - :math:`1.0 \times 10^{-15}`
     - :math:`1.0 \times 10^{-15}`
   * - the same, by subtraction
     - :math:`1.1 \times 10^{-15}`
     - :math:`3.1 \times 10^{-10}`
     - :math:`1.3 \times 10^{-7}`
     - :math:`1.9 \times 10^{-4}`

The subtraction's error grows as the inverse of the absorption, three
decades of :math:`\Sigma_t` for three decades of error; the built
coupling is flat. The gate is
``test_a_nearly_void_closed_body_conserves_its_emission`` (four bodies,
:math:`\Sigma_t` from 1 to :math:`10^{-12}`, ERR-102).

The lossless body is refused exactly
------------------------------------

When every diffuse wall returns everything (:math:`\alpha = 1`) and every
loss is exactly 0, the body loses nothing: :math:`I - T\alpha` is exactly
singular (its column sums are the balance, all 0), and a source that
reaches the walls has no finite flux. The loss is exactly 0 only when no
traversal attenuates and nothing leaks, because the absorbed fraction is
formed by ``expm1`` and is positive for any :math:`\Sigma_t > 0`.
``WallCoupling.currents`` tests that predicate, not a computed
determinant, and raises
:class:`~orpheus.derivations.continuous.characteristic.closure.TrappedSource`
(*a source in a body that absorbs nothing, behind walls that return
everything, has no finite flux*), the white analogue of the line part's
trapped line. With no source reaching the walls (a zero or an empty
escape) it returns no current, a zero ``(W, M)``, instead: the
test-architect found on
the first run of ``test_a_lossless_body_with_no_source_has_an_empty_block``
that the solve was attempted for zero columns and raised ``LinAlgError``,
and the main agent added the zero branch (2026-10-06).

Formed by subtraction, the same body is not refused: the computed
transmission misses 1 by a rounding, and the solve divides by it. `[M]`
2026-10-07, the archivist's probe: a void sphere :math:`(0, 0.5, 1.5, 2)`
behind a white wall of amplitude 1, every region in the support, at 8
points: loss exactly 0, :math:`1 - T = 1.1 \times 10^{-16}`, and the
subtraction's block reaches :math:`1.1 \times 10^{16}` for a source that
has no finite flux; the built coupling raises. The refusal lives in
``currents``, beside the balance row that makes the singularity it guards
exact, and not as a predicate on the input (the user's ruling of
2026-10-06 to fix the near-void conditioning in the rung). Gates: ``test_a_source_behind_walls_that_return_everything_in_a_lossless_body_is_refused``,
with one leg per condition dropped (:math:`\alpha = 0.99`, one region
absorbing) building finite blocks, and a void body behind a mirror refused
by the line part's own fragment and not this one.

A wall both specular and diffuse would break the balance as tallied
(it counts a diffuse wall's specular return as zero), and ``Wall`` refuses
it (:ref:`characteristic-walls`).

The white walls' gates
----------------------

Against closed forms written in mpmath in
``tests/gates/derivations/_characteristic_mp.py``, which imports nothing
from ``orpheus``:

- **The escape probability two ways and the wall transmission**
  (``test_the_escape_and_transmission_probabilities_are_the_closed_forms``):
  a homogeneous body behind white walls, :math:`P_{\rm esc}` as
  :math:`1 - \Sigma\,\mathbf 1^{\mathsf T}K_{\rm line}\mathbf 1/V` and as
  :math:`\mathbf 1^{\mathsf T}U/V`, against Hébert's sphere, the slab's
  :math:`(1 - 2E_3(\tau))/(2\tau)` and the cylinder's
  :math:`(1 - P_{ss})/(2\tau)` through Bickley's :math:`\mathrm{Ki}_3`; the
  wall transmission against :math:`P_{ss}` (sphere, cylinder) or
  :math:`2E_3(\tau)` face to face (slab). These read every line with its
  weight, so they are where the rule over lines is gated.
- **The white and the specular laws**
  (``test_the_white_and_specular_laws_are_their_closed_forms_and_differ``):
  :math:`\mathbf 1^{\mathsf T}K\mathbf 1` against the balance of one
  re-emission chain,
  :math:`(V/\Sigma)\,[(1 - P_{\rm esc}) + \alpha P_{\rm esc}(1 - T)/(1 - \alpha T)]`,
  and the specular sphere against its per-line integral, to
  :math:`10^{-13}`; at :math:`\alpha = 0.5` the white and the specular
  spheres differ by more than 100 times the conservation band, the
  discriminator a specular-only code cannot pass.
- the reciprocity rows, the balance row, the void rows, the near-void
  sweep and the refusals above.

`[M]` 2026-10-07, the archivist's probe, homogeneous bodies of radius or
width 1.3 behind white walls (relative errors; :math:`P_{\rm esc}` from the
line part, from :math:`U`, then :math:`T_w`):

.. list-table::
   :header-rows: 1
   :widths: 26 10 21 21 22

   * - Body, :math:`\tau`
     - points
     - :math:`P_{\rm esc}` (line)
     - :math:`P_{\rm esc}` (:math:`U`)
     - :math:`T_w`
   * - sphere, 0.5
     - 16
     - :math:`-1.1 \times 10^{-16}`
     - :math:`2.2 \times 10^{-16}`
     - :math:`-2.2 \times 10^{-16}`
   * - sphere, 2
     - 16
     - :math:`-5.6 \times 10^{-16}`
     - :math:`2.2 \times 10^{-16}`
     - :math:`-2.2 \times 10^{-16}`
   * - sphere, 100
     - 16
     - :math:`-8.8 \times 10^{-14}`
     - :math:`4.4 \times 10^{-16}`
     - :math:`-9.2 \times 10^{-13}`
   * - sphere, 1000
     - 16
     - :math:`-1.3 \times 10^{-12}`
     - :math:`1.1 \times 10^{-15}`
     - :math:`6.1 \times 10^{-11}`
   * - slab, 0.01
     - 8
     - :math:`9.8 \times 10^{-15}`
     - :math:`9.8 \times 10^{-15}`
     - :math:`-1.1 \times 10^{-16}`
   * - slab, 0.5
     - 8
     - :math:`-1.1 \times 10^{-14}`
     - :math:`-1.1 \times 10^{-14}`
     - :math:`1.3 \times 10^{-14}`
   * - slab, 8
     - 8
     - :math:`8.9 \times 10^{-16}`
     - :math:`-1.1 \times 10^{-16}`
     - :math:`4.4 \times 10^{-16}`
   * - slab, 30
     - 8
     - :math:`-1.1 \times 10^{-14}`
     - :math:`-2.2 \times 10^{-16}`
     - :math:`2.2 \times 10^{-15}`

The cylinder's rows are ``slow`` (:ref:`characteristic-line-rule`, the
cost); the main agent measured them on the built code: against Bickley at
:math:`\tau = 0.5`, :math:`1.4 \times 10^{-12}` and
:math:`1.9 \times 10^{-15}` at 8 and 16 points, and
:math:`2.9 \times 10^{-14}` at :math:`\tau = 0.01` at 8 points.


.. _characteristic-line-rule:

The line rule, every grading derived from the optical scale
===========================================================

:class:`~orpheus.derivations.continuous.characteristic.assembly.LineRule`
is the quadrature over the line domain, for one group, the role of a
weighted set of lines
(:class:`~orpheus.derivations.continuous.characteristic.lines.Lines`)
whose functional is the block; the directions at a point are the other
role (:ref:`characteristic-reading`), built by the same one-dimensional
rules, in :mod:`~orpheus.derivations.continuous.characteristic.lines`.
``LineRule.of(basis, walls, sigma_t, points, chunk=512, budget=1024)``
places its nodes on the chart's
:class:`~orpheus.geometry.chart.LineDomain`, ``points`` Gauss–Legendre
points per piece of each coordinate, and gives each line the weight
:math:`w_L = W_L\varrho/4\pi`. It is built per group because its gradings
read the group's total cross sections: the user ruled on 2026-10-06, after
qa found the fixed rules blind, that **every grading is derived from the
group's optical scale**, by the law the traversal rule already follows
along a line, a piece no wider than the feature it must resolve; the
alternative, refusing the regimes the fixed rules could not reach, was
declined. ``LineRule.transport(points, inner_points, support)`` then
assembles the block (:ref:`characteristic-galerkin-assembly-section`).

The three gradings
------------------

Every piece end the package places comes from one of three laws, each
spelled once in
:mod:`~orpheus.derivations.continuous.characteristic.grading` and imported
by every caller:

- **geometric**,
  :func:`~orpheus.derivations.continuous.characteristic.grading.graded_ends`:
  ``layers`` ends at the depths :math:`w\rho^j` toward an end. The panels
  toward walls and interfaces use it (:ref:`characteristic-panel-basis`),
  and so do the grazing halvings below, at :math:`\rho = 1/2`.
- **exponential**,
  :func:`~orpheus.derivations.continuous.characteristic.grading.exponential_ends`:
  ends at :math:`2^k` mean free paths, :math:`k = 0, \dots, 6`, toward where
  an attenuation concentrates, until it vanishes:
  ``VANISHING_DEPTH`` :math:`= 2^6 = 64`, past which :math:`e^{-64}` is
  below double precision; an interval of optical width at most 2 is not
  cut. The traversal rule's pieces (:ref:`characteristic-pieces`) and the
  rim and normal-direction gradings below use it.
- **hp**,
  :func:`~orpheus.derivations.continuous.characteristic.grading.halvings`:
  :math:`\lceil\log_2(w/d)\rceil` halvings of a piece of width :math:`w`
  until it is no wider than its distance :math:`d` to a singularity off the
  interval, at most ``np.finfo(float).nmant``; a positive width at a
  distance of zero or less (a singularity on the interval) is refused.
  Gauss–Legendre then converges geometrically on every piece
  (:ref:`characteristic-branch-grading`). The orbit coordinate's branch
  points (the traversal rule), the next radius and :math:`b = 0` (the
  impact rule) and the grazing direction (the direction rules) use it.

The module exists because the elegance re-review of 2026-10-07 found the
hp law spelled three ways, in the traversal rule, the impact rule and the
direction rules, and the traversal rule's spelling one halving short of
the stated law (:ref:`characteristic-branch-grading`).

The impact rule, in the chord half-length
-----------------------------------------

On the sphere and the cylinder the impact parameter :math:`b` runs over
:math:`[0, r_n]`, cut at every panel end of the basis (a hollow body adds
a void impact panel :math:`[0, r_0]` for the lines through the cavity,
which cross no panel there). On an impact panel :math:`[r_k, r_{k+1}]` the
integration variable is the half-length of the chord through the circle
:math:`r_{k+1}`,

.. math::

   y \;=\; \sqrt{r_{k+1}^{2} - b^{2}} \;\in\; [0,\ y_{\max}],
   \qquad y_{\max} = \sqrt{(r_{k+1} - r_k)(r_{k+1} + r_k)},
   \qquad \mathrm{d}b \;=\; \frac{y}{b}\,\mathrm{d}y ,

the visibility-cone substitution, which absorbs the square-root end at
:math:`b = r_{k+1}` that every chord through that circle carries
(``chord_quadrature``'s, applied here per panel and in the variable the
gradings are stated in). The integrand still has singularities off the
interval and concentrations on it, and :math:`y` is graded toward each
(:func:`~orpheus.derivations.continuous.characteristic.lines.impact_rule`):

- **The next radius's branch point**, hp toward :math:`y = 0`. A line at
  :math:`b` just inside :math:`r_{k+1}` also meets the next radius
  :math:`r_{k+2}`, whose chord
  :math:`\sqrt{r_{k+2}^2 - b^2} = \sqrt{(r_{k+2}^2 - r_{k+1}^2) + y^2}`
  has branch points at :math:`y = \pm i\sqrt{r_{k+2}^2 - r_{k+1}^2}`. In
  :math:`b` that singularity sits a real distance
  :math:`r_{k+2} - r_{k+1}` past the panel's end; in :math:`y` it is
  :math:`\sqrt{(r_{k+2} - r_{k+1})(r_{k+2} + r_{k+1})} \approx \sqrt{2r\,\delta}`
  away, much farther, which is why the substitution absorbs most of it and
  the hp grading is needed only where a wide panel sits just inside a
  thin one. It is one of three features near :math:`y = 0`; the turning
  slot's exponential layer and the closure's pole are the other two, and
  the hp grading stops at the nearest of the three (the grading law,
  :ref:`characteristic-grading-law`).
- **The point** :math:`b = 0`, hp toward :math:`y = y_{\max}` (that is,
  :math:`b = r_k`): it lies at :math:`y = r_{k+1}`, a distance
  :math:`r_{k+1} - y_{\max} = r_k^2/(r_{k+1} + y_{\max})` beyond the
  interval, written so that a tiny :math:`r_k` does not cancel. There the
  substitution's Jacobian :math:`y/b` is singular, and so is the Abel
  transform of the turning panel's odd modes, :math:`b^{2m+2}\log b`
  (:ref:`characteristic-even-basis`).
- **The rim**, exponentially toward :math:`y = 0`, at :math:`2^j` mean
  free paths of the thickest panel the line crosses (those from
  :math:`r_k` outward). The length the line travels inside the circle
  :math:`r_{k+1}` is :math:`2y`, so its optical length there grows from 0
  like :math:`\Sigma y`, and on an optically thick panel the integrand
  changes on the scale of one mean free path in :math:`y` near
  :math:`y = 0`, the rim of that circle.

The panel touching :math:`b = 0` is smooth in :math:`b^2`, because the
panel the line turns in is the even one or a cavity, so its lower half is
plain Gauss–Legendre in :math:`b` and only its upper half is graded as
above. :math:`b` is recovered as
:math:`\sqrt{r_k^2 + (y_{\max} - y)(y_{\max} + y)}`, which tends to
:math:`r_k` without the cancellation of :math:`r_{k+1}^2 - y^2`. Each node
keeps its panel's top :math:`r_{k+1}` and its own :math:`y`
(:class:`~orpheus.derivations.continuous.characteristic.lines.ImpactRule`),
and its line is chorded from that pair, its exact level, never from
:math:`b`: a line stores :math:`b` to an ulp, and near the panel top an
ulp of :math:`b` is the whole half-chord (:ref:`chart-and-chord-level`,
ERR-106).

.. dropdown:: First got wrong: the small first radius read as the centre's singularity
   :color: muted

   qa measured closed-body conservation off by :math:`1.1 \times 10^{-8}` at 16 points on a sphere
   whose first region has radius :math:`10^{-3}`, under the impact rule of
   the time, ``chord_quadrature`` in :math:`b`. It read as the
   :math:`b = 0` singularity. The main agent's measurement placed it
   elsewhere: the wide middle impact panels saw the **next** radius's branch
   point just past their upper end, at a real distance of 0.07; Gauss in
   :math:`b` on the panel :math:`[0.52, 0.88]` left :math:`10^{-6}` at 8
   points, with the neighbours at 0.952 and 1.0. The grading toward the next
   radius is the fix of that, and the hp grading toward :math:`b = 0` the fix
   of the other.

**A19, wide then thin.** The row written as the next-radius grading's
witness, ``test_a_small_first_region_conserves_at_eight_points``, is not
one on the final code: the substitution to :math:`y` absorbs most of that
branch point, and `[M]` 2026-10-07 (the test-architect, ``gates/a19_probe.py``)
the row reads :math:`3.5 \times 10^{-12}` at 8 points with the grading and
:math:`8.1 \times 10^{-12}` without, a margin no tolerance can stand on;
the row is declared blind to that arm and stays as the region's value
witness. The grading matters where a wide panel sits just inside a thin
one, which the basis's own interface grading produces: on a hollow sphere
:math:`(0.2, 0.3, 1.0)` graded 10 layers deep at ratio 0.4, inner mirror,
outer white, 8 points, conservation reads :math:`1.0 \times 10^{-12}` with
the grading and :math:`2.0 \times 10^{-9}` without (the test-architect,
``gates/a19_wide.py``; the main agent's own fixture of that shape,
:math:`9.2 \times 10^{-11}` and :math:`1.7 \times 10^{-3}`).
``test_a_wide_panel_before_a_thin_one_conserves_at_eight_points`` is the
witness, at :math:`10^{-11}`.

.. _characteristic-grading-law:

The grading law at every impact-panel top
-----------------------------------------

The hp grading toward :math:`y = 0` on the impact panel
:math:`[r_k, r_{k+1}]` stops at a distance :math:`d_k`, the distance in
:math:`y` from the tangency :math:`b = r_{k+1}` to the nearest feature of a
line's transport there (:meth:`ImpactPanels.distances
<orpheus.derivations.continuous.characteristic.lines.ImpactPanels.distances>`,
the grading law G). A line just below :math:`r_{k+1}`, at
:math:`y = \sqrt{r_{k+1}^2 - b^2} \to 0`, turns in panel :math:`k` across a
slot of in-plane length :math:`2y`, crosses the shells above, and reaches
the outer wall. Along it, at projected speed :math:`s = |P\Omega|`, three
features sit near :math:`y = 0`:

- **the next radius out**, :math:`b = r_{k+2}`, whose half-chord
  :math:`\sqrt{(r_{k+2}^2 - r_{k+1}^2) + y^2}` has its branch points at
  :math:`y = \pm i\sqrt{r_{k+2}^2 - r_{k+1}^2}` (the bullet above);
- **the layer**: the turning slot transmits :math:`e^{-2\Sigma_k y/s}`,
  which changes over :math:`y \sim s/(2\Sigma_k)`;
- **the pole**: the line's period carries the cycle product
  :math:`\Pi = a\,e^{-(\tau_{\rm out} + 2\Sigma_k y)/s}`
  (:eq:`characteristic-closure`), with :math:`a` the outer wall's specular
  amplitude and :math:`\tau_{\rm out}` the in-plane optical depth of the
  shells above :math:`r_{k+1}` along the tangent line, both sides. The
  closure :math:`1/(1 - \Pi)` has a pole where :math:`\Pi = 1`,

  .. math::

     \ln a - \frac{\tau_{\rm out} + 2\Sigma_k y}{s} = 0
     \quad\Longleftrightarrow\quad
     y \;=\; -\,\frac{s\,(-\ln a) + \tau_{\rm out}}{2\Sigma_k},

  a real distance :math:`(s(-\ln a) + \tau_{\rm out})/(2\Sigma_k)` from
  the interval, on its far side. It nears the tangency where :math:`a \to 1`
  and the shells above are thin or void.

The law takes the nearest:

.. math::

   d_k \;=\; \min\!\Bigl(\sqrt{(r_{k+2} - r_{k+1})(r_{k+2} + r_{k+1})},\;
   \frac{\min\bigl(s,\ s\,(-\ln a) + \tau_{\rm out}(r_{k+1})\bigr)}{2\Sigma_k}\Bigr),

with:

- :math:`s` the projected speed of the lines the impact rule is built for:
  1 on the sphere, and on the cylinder each polar node's own
  :math:`\sin\theta_j`, one impact rule per polar node (the iterated rule,
  :ref:`characteristic-iterated-cylinder-rule`). The layer's and the
  pole's distances are affine in :math:`s`, so each line is graded to the
  features it carries and to no slower line's;
- the :math:`-\ln a` term present only for :math:`0 < a < 1`: a wall with
  no specular return has no pole, and at :math:`a = 1` the closure is
  regular (its source integral vanishes with :math:`1 - \Pi`);
- :math:`d_k = \infty` for the transport terms of a void panel
  (:math:`\Sigma_k = 0`), whose transport does not change with :math:`y`,
  and the next radius's term absent on the outermost panel;
- :math:`\tau_{\rm out}` and the next radius's half-chord read from the
  kernel's chord of the tangent lines :math:`b = r_{k+1}` through the
  impact panels, at unit speed (``traversed_length`` and
  ``half_chord_at``), so no chord length is spelled a second time.

The data that do not depend on the speed (the impact panels, their cross
sections, the next radius's half-chord, :math:`\tau_{\rm out}` and
:math:`-\ln a`) are read once per group by :meth:`ImpactPanels.of
<orpheus.derivations.continuous.characteristic.lines.ImpactPanels.of>`;
:meth:`~orpheus.derivations.continuous.characteristic.lines.ImpactPanels.distances`
evaluates :math:`d_k` at one speed, and
:meth:`~orpheus.derivations.continuous.characteristic.lines.ImpactPanels.rule`
builds the impact rule
(:func:`~orpheus.derivations.continuous.characteristic.lines.impact_rule`)
graded by it, over every panel for the line rule or over the panels below
a point's orbit coordinate for the point rule (its ``below``).

The outermost panel (:math:`\tau_{\rm out} = 0`) is the rim of the body,
where :math:`d = \min(s, s(-\ln a))/(2\Sigma_{n-1})` is the rule's first
form, the rim law. The inner mirror of a hollow body needs no term: its
tangency's panel is the void cavity, which carries neither a layer nor a
pole near :math:`y = 0` (qa's second round measured a hollow sphere with
both mirrors at 0.999, read on the inner wall, to
:math:`8.2 \times 10^{-14}` at 8 line points and :math:`7.8 \times 10^{-16}`
from 16). An under-estimate of :math:`\tau_{\rm out}`, or dropping it,
only grades more finely, so no value gate can redden on it: only a cost
gate could (designed green, qa's second round).

**The evidence.** `[M]` 2026-10-08, qa's second round
(``scratch/characteristic_architecture/p1_step_b5b/qa/r2/``) against the
mpmath routes at 8, 16 and 32 line points: a void outer region under
:math:`a = 0.99`, read at the interface :math:`7.1 \times 10^{-15}`,
:math:`4.4 \times 10^{-16}` and 0, on the wall at most
:math:`4.4 \times 10^{-16}`; a thin outer region :math:`4.1 \times 10^{-14}`
at 8; :math:`a = 1 - 10^{-9}` behind a void outer region at most
:math:`6.7 \times 10^{-16}`; an interior void gap at :math:`a = 0.99`
:math:`3.8 \times 10^{-15}` at 8; the block's
:math:`\mathbf 1^{\mathsf T}K\mathbf 1` under a void outer region at
:math:`a = 0.99` within :math:`3.3 \times 10^{-16}` of the 64-point
block; the three-region vacuum cylinder at its interface
:math:`x = 1.5`, :math:`1.4 \times 10^{-15}` and :math:`2.2 \times 10^{-16}`
at 8 and 16. The gates: the void and thin outer regions read and totalled
(``test_a_void_or_thin_outer_shell_under_a_near_one_mirror_is_read_to_the_bar``,
``test_the_blocks_total_under_a_void_or_thin_outer_shell_is_the_closed_form``,
ERR-104); the partial mirror's wall and the block at 0.6, 0.9 and 0.99
(``test_the_reading_on_and_near_a_partial_mirror_is_the_mpmath_route``,
``test_the_blocks_total_under_a_partial_mirror_is_the_closed_form``); the
small cylinder's layer, the fast catcher of the layer term
(``test_a_small_cylinders_layer_is_graded_outside_slow``, ERR-105), and the
three-region cylinder's interface (``slow``). Removing the law on the final
tree (G returning the next radius alone, which is the third rung's rule)
reddens 21 of the 172 rows of the reading's and the level's files outside
``slow``; removing the pole term alone, 18; removing the layer term alone,
both rows of the small cylinder (the archivist's re-drop, 2026-10-08,
``-O``).

**The cylinder's cost.** The law's cost on the cylinder is set by its
slowest lines: a polar node at :math:`\sin\theta_j` grades each panel top
to :math:`\sin\theta_j/(2\Sigma_k)`, so the impact rule of a near-grazing
polar node holds several times the nodes of a fast one (240 to 1928
impact nodes per polar node on the ``ABA`` cylinder's first group). Each
polar node pays only for its own speed
(:ref:`characteristic-iterated-cylinder-rule`, the measured line counts
and times).

.. dropdown:: First got wrong: the rim law alone, and what the floor hid
   :color: muted

   The first build of the rung graded only the outermost panel, by a rim
   law: the pole :math:`-\ln a\,s_{\min}/(2\Sigma_{n-1})` of the outer
   wall, plus the layer, applied where the panel is the last and infinite
   for a void outer panel ("a void outer panel changes nothing"). It
   repaired the defect the reading found on the wall: `[M]` 2026-10-08, a
   sphere's wall reading under :math:`a = 0.99` went from
   :math:`9.3 \times 10^{-5}` to :math:`1.2 \times 10^{-13}` at 12 line
   points, and the block's :math:`\mathbf 1^{\mathsf T}K\mathbf 1` from
   :math:`5.3 \times 10^{-10}` to :math:`3.8 \times 10^{-16}`. qa's review
   found it incomplete twice (``p1_step_b5b/qa.md``):

   - **F1, the pole at an interior tangency.** Below an outer region that
     absorbs less per length than the one inside it, a void or thin shell
     under a near-1 mirror, :math:`1 - \Pi \approx (-\ln a + \tau_{\rm out})
     + 2\Sigma_k y`, and the pole sits at
     :math:`(-\ln a + \tau_{\rm out})/(2\Sigma_k)` from the interior
     tangency, which the next-radius grading covers only when the shells
     above are at least as thick. `[M]` qa's ``p1_hidden_pole``: a void
     outer region at :math:`a = 0.99` missed by :math:`5.1 \times 10^{-4}`,
     :math:`2.6 \times 10^{-5}` and :math:`5.8 \times 10^{-8}` at the
     interface at 8, 16 and 32 line points, and the block's total by
     :math:`2.6 \times 10^{-6}` at 8;
   - **F2, the layer at an interior radius.** The cylinder's layer
     :math:`e^{-2\Sigma y/\sin\theta}` sits at every tangency, and a point
     on an interface sees it at first order: the three-region vacuum
     cylinder at :math:`x = 1.5` missed by :math:`1.0 \times 10^{-6}` and
     :math:`1.8 \times 10^{-8}` at 8 and 16 line points.

   **A correction the elegance review measured.** The first build was
   credited with fixing the cylinder's vacuum wall (:math:`7.9 \times
   10^{-8}` at 8 line points before it). It had not: the rim law's
   distance on the cylinder was scaled by the slowest polar speed, about
   :math:`10^{-7}`, which drove it below a floor the same build had placed
   on every grading distance in :math:`y`, :math:`\sqrt{2r\,\epsilon(r)}`
   times a margin over the first Gauss node's fraction (``_B_MARGIN``),
   so that no impact node rounded onto a panel top. The floor set the
   cylinder's grading whatever the rim law said: the elegance review's
   probe 3 found the rule built with the computed rim identical to the
   rule built with rim 0 on 2 of 2 cylinder fixtures. The rim law acted on
   the sphere (:math:`s = 1`).

   The floor itself stood where the line rule lacked the exact half-chord:
   a node graded below :math:`b`'s resolution rounds :math:`b` onto the
   panel top and makes a tangent line (on a cylinder's wall under
   :math:`a = 0.9` and 0.99, 4 and 52 impact nodes; the block dropped those
   lines silently and the reading refused them). The image's exact level
   (:ref:`chart-and-chord-level`) removes the cause, and the floor and its
   margin were retired with it; after that the grading law, and not the
   floor, grades the cylinder (the interface reading above). The law
   replaced the rim law (the user's ruling of 2026-10-08, "the general
   grading law lands in 5b"), and its catalogue entries are ERR-104 (the
   pole) and ERR-105 (the layer).

.. _characteristic-iterated-cylinder-rule:

Each polar angle at its own speed: the cylinder's iterated rule (#587)
----------------------------------------------------------------------

The cylinder's lines are the box :math:`(b, \theta) \in [0, R] \times
[0, \pi/2]` of the line domain (:eq:`geometry-line-domain`), with the
density :math:`8\pi\sin^2\theta` of
:meth:`LineDomain.density <orpheus.geometry.chart.LineDomain.density>`.
Both roles integrate a functional of the line over it: the line rule's
block and the point rule's row. The rule over the box is **iterated**,
not a tensor product:

.. math::

   \int_0^{\pi/2}\!\mathrm d\theta \int_0^R\!\mathrm db\; f(b, \theta)
   \;\approx\; \sum_{j} w^\theta_j \sum_{i} w^b_{ij}\, f(b_{ij}, \theta_j),

where :math:`(\theta_j, w^\theta_j)` is the polar rule
(:func:`~orpheus.derivations.continuous.characteristic.lines.polar_rule`,
graded toward grazing and toward the normal over the whole body,
:ref:`characteristic-line-rule`) and :math:`(b_{ij}, w^b_{ij})` is the
impact rule built for the lines of polar angle :math:`\theta_j` alone,
graded by the grading law at that line's own projected speed
:math:`s_j = \sin\theta_j`.
:func:`~orpheus.derivations.continuous.characteristic.lines.impact_per_polar`
builds it from a function of the speed (``ImpactPanels.rule`` at the
group's panels) and the polar rule, and returns the plain
:math:`\mathrm db\,\mathrm d\theta` weights :math:`w^\theta_j w^b_{ij}`;
each role multiplies its own density afterwards, the line rule
:math:`8\pi\sin^2\theta/4\pi`
(:class:`~orpheus.derivations.continuous.characteristic.assembly.LineRule`)
and the point rule its direction measure through the point
(:class:`~orpheus.derivations.continuous.characteristic.reading.PointRule`).
A tensor product is the special case :math:`b_{ij} = b_i` for every
:math:`j`.

**Why each line's own speed is the grading speed.** A line at polar angle
:math:`\theta` from the axis advances :math:`\sin\theta` in the plane per
unit of its own length, so an in-plane length :math:`\ell` in a region of
cross section :math:`\Sigma` is an optical depth
:math:`\Sigma\ell/\sin\theta` along it. The two features of the grading
law that depend on the line are therefore those of that line's speed
:math:`s = \sin\theta`:

- the turning slot, of in-plane length :math:`2y`, transmits
  :math:`e^{-2\Sigma_k y/s}`, a layer of width :math:`s/(2\Sigma_k)` in
  :math:`y`;
- the cycle product :math:`\Pi = a\,e^{-(\tau_{\rm out} + 2\Sigma_k y)/s}`
  makes the closure's pole a distance
  :math:`(s(-\ln a) + \tau_{\rm out})/(2\Sigma_k)` from the interval
  (:ref:`characteristic-grading-law`).

Both distances are affine in :math:`s` with a non-negative slope, and the
third feature, the next radius's branch point, does not depend on
:math:`s`. At a fixed :math:`\theta_j` the integrand
:math:`b \mapsto f(b, \theta_j)` carries the features of the speed
:math:`s_j` only, so the impact rule graded at :math:`s_j` resolves it
exactly as the law requires. A rule graded at a slower speed
:math:`s_{\min} < s_j` subdivides the stretch nearer the tangency than
:math:`d_k(s_j)` at scales finer than this line's integrand varies on, and
buys no accuracy; a rule graded at a faster one misses the line's layer and
pole. The hp grading adds a piece per halving of the distance, so the
slower speed costs about :math:`\log_2(s_j/s_{\min})` pieces per graded
panel top (`[R]`; measured: 30 pieces at :math:`s = 1` and 241 at the
slowest polar node, :math:`\sin\theta = 4.85 \times 10^{-6}`, on the
first group of the ``ABA`` cylinder below, 240 and 1928 impact nodes of
8 points per piece).

A tensor product cannot follow the speed: its impact nodes are shared by
every polar node, so they must resolve the narrowest layer and the
nearest pole of any polar node it holds, those of the slowest. The
iterated rule is the quadrature of the iterated integral
:math:`\int\mathrm d\theta\,\int\mathrm db` (Fubini: the integrand is integrable on the box),
the one whose inner rule may depend on the outer node. Its error is the
polar rule's error on :math:`g(\theta) = \int_0^R f(b, \theta)\,\mathrm db`
plus :math:`\sum_j w^\theta_j E_j`, with :math:`E_j` the impact rule's
error at :math:`\theta_j`, and each :math:`E_j` is controlled by the
grading law at :math:`s_j`. The sphere (:math:`s = 1` for every line) and
the slab (no impact parameter) are unchanged.

**The evidence.** `[M]` 2026-10-08, the main agent
(``scratch/characteristic_architecture/p1_step_c/aba_cyl_cost.py``,
``cost_5a.log``, ``cost_head.log``, ``iterated.log``): the ``ABA``
cylinder (``aba_specification(CYLINDRICAL)``) at resolution
``Resolution(3, 2, 0.4, TransportResolution(8, 12, 12), 8)``, the k solve,
one run each (no repeat protocol), lines per block of the first group:

.. list-table::
   :header-rows: 1
   :widths: 40 16 14 30

   * - Line rule
     - Lines per block
     - Wall time
     - k
   * - the fifth rung's first half, no grading law at the panel tops
       (``fe977a90``)
     - 38 400
     - 384 s
     - 1.231729452890078
   * - the tensor rule, every polar angle at the slowest speed
       (``11b263a5``)
     - 308 480
     - 2296 s
     - 1.231729452890264
   * - the iterated rule, each polar angle at its own speed
     - 102 320
     - 1275 s
     - 1.231729452890254

The iterated rule's k is the tensor rule's to :math:`8.3 \times 10^{-15}`
relative, and the line count is a third of it (the second group: 106 784
against 297 920). Its 1275 s was taken by an in-memory probe of the same
rule (``iterated_probe.py``). `[M]` 2026-10-09, the archivist's re-count
through the shipped ``ImpactPanels``
(``p1_step_c/archivist_587_counts.py`` and its log, no solve): 102 320 and 106 784 lines, 160
and 152 polar nodes, 240 impact nodes per polar node at :math:`s = 1` to
1928 at the slowest; the tensor rule's count is :math:`1928 \times 160 =
308\,480` and :math:`1960 \times 152 = 297\,920`. The k of the first row
differs from the other two by :math:`1.5` and :math:`1.4 \times 10^{-13}`
relative: it lacks the grading law, not the iteration.

`[M]` 2026-10-08, qa (``p1_step_c/qa.md``): no accuracy is lost.
The three-region vacuum cylinder read at :math:`x = 1.5`, the small
cylinder :math:`(0, 0.5, 1)` at its interface and its wall, the white
cylinder's :math:`T_w` and :math:`P_{\rm esc}` at :math:`\tau = 0.5`,
0.01 and :math:`10^{-4}`, and the hollow cylinder's conservation agree
with the tensor rule's within :math:`2 \times 10^{-15}` on every row.
The measure is exact: the weights
over the density sum to :math:`R\pi/2`, and each polar node's impact
weights to :math:`R`, within :math:`3.3 \times 10^{-15}` on 16 cylinder
block rules. With every angle forced to the slowest speed, the iterated
assembler reproduces the tensor rule's readings to all 17 printed digits.
Under in-process mutations of the speed, against Bickley's
:math:`\mathrm{Ki}_2` route:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Mutation
     - three-region cylinder, :math:`x = 1.5`
     - small cylinder, interface
   * - the layer term removed (ERR-105's arm)
     - :math:`1.01 \times 10^{-6}`
     - :math:`7.97 \times 10^{-8}`
   * - every angle graded at speed 1
     - :math:`1.27 \times 10^{-7}`
     - :math:`7.97 \times 10^{-8}`
   * - :math:`\cos\theta` in place of :math:`\sin\theta`
     - :math:`1.27 \times 10^{-7}`
     - :math:`7.97 \times 10^{-8}`
   * - every angle graded at the slowest speed
     - unchanged
     - unchanged

Each angle's own speed is load-bearing: grading faster than it misses the
layer, and grading slower costs lines and changes no value. The sphere's
and the slab's rules are bit-identical to the tensor rule's tree, 86 of 86
rules (the test-architect, ``p1_step_c/ta/bitid.log``: 10 sphere and 4 slab
line rules, 56 sphere and 16 slab point rules), and 176 of 176 arrays in
qa's own comparison.

**The gates.**
``test_each_polar_angle_carries_the_impact_rule_graded_at_its_own_speed``
(``tests/gates/derivations/test_characteristic_assembly.py``, vacuum and
:math:`a = 0.9`) and
``test_each_polar_angle_through_a_point_carries_the_impact_rule_graded_at_its_own_speed``
(``tests/gates/derivations/test_characteristic_reading.py``, the wall, the
interface and the wall under :math:`a = 0.9`) assert, ``array_equal``,
that the lines at each polar node are the impact rule graded by the
grading law written by hand in the test at that node's
:math:`\sin\theta`, with their levels. `[M]` 2026-10-08, the
test-architect (``p1_step_c/ta/README.md``): on ``11b263a5`` 127 of 128
polar nodes differ from their own-speed rule on all five rows, and 0 of
128 on the iterated rule. Of the 16 rows of its battery, the arm that
grades every angle at the slowest speed reddens these five and no value
row, so they are its only catchers: a slower grading is invisible to a
value gate by design (the paragraph above).

**Where the remaining lines are.** `[M]` 2026-10-09, the archivist's
re-count: on the first group 51 of 160 polar nodes have
:math:`\sin\theta < 0.01` and carry 57 536 of the 102 320 lines (the
second group: 43 of 152 nodes, 52 872 of 106 784 lines). They are there
because the polar rule grades toward grazing down to
:math:`\tau_{\min}/64` of the whole body (:ref:`characteristic-line-rule`,
the direction rules), at every impact parameter, whereas a line of short
chord near the rim would not need it. Grading :math:`\theta` per impact
panel was the other half of #587; it is not built, and on this body it
would not cut them, because every line crosses the thin outer-wall panel
that sets the grazing grading (`[M]` 2026-10-09,
:ref:`characteristic-what-is-not-built`).

.. dropdown:: First got wrong: every polar angle graded at the slowest
   :color: muted

   When the grading law landed with the fifth rung's second half
   (``b76b9a9d``), the cylinder's rule was the tensor product of one
   impact rule and the polar rule (the retired ``impact_by_polar``), and
   the law graded that one impact rule at the slowest projected speed the
   polar rule sampled, :math:`\sin\theta_{\min}` (the retired
   ``tangency_distances``, through its ``slowest`` argument). That was
   correct, since the slowest line carries the narrowest layer and the
   nearest pole, and wasteful: every other polar node carried the
   slowest one's impact nodes. On the ``ABA`` cylinder the block grew from
   38 400 lines to 308 480, 8.0 times the ungraded rule's, and the k solve
   from 384 s to 2296 s; the point rule's line counts on the cylinder grew
   3.6 and 5.3 times (one and three regions). The user accepted the cost until #587, which split the
   law's speed-free data (``ImpactPanels``) from its evaluation at a
   speed and made the rule iterated (``200b4233``): 102 320 lines and
   1275 s, the same k
   to :math:`8.3 \times 10^{-15}`.

The direction rules: grazing and normal
---------------------------------------

On the slab the direction coordinate is the cosine :math:`\mu` on each
sign, on the cylinder the polar angle :math:`\theta`; both run from the
**grazing** direction, projected speed :math:`v = |P\Omega| = 0`, to the
**normal** one, :math:`v = 1` (``_grazing_ends`` in ``lines.py``, with
:math:`\mu = v` and :math:`\theta = \arcsin v`). A line crossing a panel of
optical width :math:`\tau_P` normally crosses :math:`\tau_P/v` along
itself, and its attenuation :math:`e^{-\tau_P/v}` has two features:

- **Toward grazing**, the attenuation of a panel changes until
  :math:`\tau_P/v` reaches the depth where it vanishes, 64. So the
  coordinate is halved toward :math:`v = 0` down to
  :math:`v = \tau_{\min}/64`, with :math:`\tau_{\min}` the thinnest
  absorbing panel's optical width: ``halvings(1, tau_min / VANISHING_DEPTH)``
  geometric layers at ratio 1/2. The thinnest **panel**, not the body,
  sets it: a line grazing a thin panel is optically thick in it at
  :math:`v \sim \tau_P`. With no absorbing panel the coordinate is not
  graded.
- **Toward the normal**, :math:`e^{-\tau s}` with :math:`s = 1/v \in [1, \infty)`
  is an exponential layer of width :math:`1/\tau` in :math:`s` on a thick
  body. So :math:`s` is graded at :math:`2^k` over the body's normal
  optical depth, the exponential ends on :math:`[1, \infty)`, mapped back
  to :math:`v = 1/s`. The depth is the sum of the panels' optical widths,
  once across the slab and twice (a diameter) on the cylinder.

Gauss–Legendre on each resulting piece of :math:`\mu` or :math:`\theta`,
``points`` per piece.

The regimes that test a rule over lines
---------------------------------------

A rule over lines is right only for the inputs whose scales it resolves, and
closed-body conservation cannot tell (it holds line by line). The regimes
that read the gradings are the walls' closed forms near void, at a thick
rim, beside a small radius and on a thick slab, and single block entries
rather than totals; each is a row of the gates (ERR-101, ERR-103).

.. dropdown:: First got wrong: fixed resolutions over lines, and what fixed each
   :color: muted

   The first build of the rule had fixed resolutions: 12 halvings toward the
   slab's grazing cosine, plain Gauss in the cylinder's polar angle, and
   ``chord_quadrature`` in :math:`b` with no grading toward the rim. Every
   closed-body gate passed, because conservation holds line by line, and the
   closed-form gates covered :math:`\tau` from 0.5 to 8 only, where a fixed
   rule is fine. qa's review moved the inputs out of that box (its findings
   F1 to F3b, in the plan's order: the slab's grazing cosine, the cylinder's
   polar angle, and the impact parameter at a thick rim and at a small first
   radius; probes under
   ``scratch/characteristic_architecture/p1_step_b3/qa/``), the re-review
   after the fix found three more, and each is a missing input region, not a
   weak tolerance. Before and after, `[M]` 2026-10-06 and 2026-10-07 (qa, the
   test-architect and the main agent, on their probes; references in mpmath
   at 60 digits with the subtractions done there):

   .. list-table::
      :header-rows: 1
      :widths: 30 20 25 25

      * - Rule and input
        - Found by
        - Before
        - After
      * - slab grazing, 12 fixed halvings; slab :math:`(0, 0.9, 2.3)` near
          void, :math:`\Sigma` from :math:`10^{-4}` to :math:`10^{-12}`
        - qa (``p12``, ``p13``)
        - collision probability off by :math:`4 \times 10^{-2}` to
          :math:`5 \times 10^{-1}`
        - at most :math:`7.5 \times 10^{-13}` (vacuum walls),
          :math:`5.5 \times 10^{-13}` (white walls), 8 points
      * - slab grazing, 12 fixed halvings; a basis graded 10 layers deep
        - qa (``p14``)
        - a block entry off by 21 %
        - collision probability :math:`1.3 \times 10^{-15}` at
          :math:`\Sigma = 1`, :math:`2.6 \times 10^{-13}` at :math:`10^{-2}`
      * - slab grazing halved to :math:`\tau_{\min}`, not
          :math:`\tau_{\min}/64`; a basis graded 6 layers,
          :math:`\Sigma = 0.01`
        - qa, re-review (N1, ``q3``)
        - a block entry off by :math:`2.0 \times 10^{-6}` at 8 points, every
          total exact to :math:`5 \times 10^{-16}`
        - :math:`2.6 \times 10^{-10}` with 2 more halvings,
          :math:`9.1 \times 10^{-13}` with 6; the rule now reaches
          :math:`\tau_{\min}/64`
      * - cylinder polar angle, plain Gauss; :math:`\tau = 0.01`
        - qa
        - escape probability off by :math:`3.2 \times 10^{-4}`
        - :math:`2.9 \times 10^{-14}` at 8 points
      * - impact parameter, ungraded rim; sphere :math:`\tau = 100, 1000`
        - qa
        - :math:`P_{ss}` off by :math:`4.6 \times 10^{-3}` and
          :math:`7.3 \times 10^{-1}`
        - :math:`2.1 \times 10^{-13}` and :math:`6.4 \times 10^{-11}` at 16
          points (the main agent)
      * - impact parameter, ``chord_quadrature`` in :math:`b`; first region
          of radius :math:`10^{-3}`
        - qa
        - conservation off by :math:`1.1 \times 10^{-8}` at 16 points
        - :math:`8.8 \times 10^{-11}` at 8 points, :math:`2.9 \times 10^{-13}`
          at 16 (a :math:`10^{-3}` cavity)
      * - normal direction ungraded; slab transmission at
          :math:`\tau = 8, 30`
        - the test-architect, re-review (``gates/slab_tw.py``)
        - :math:`2E_3(\tau)` missed by :math:`4.1 \times 10^{-11}` and
          :math:`1.8 \times 10^{-4}` at 8 points
        - :math:`1.2 \times 10^{-15}` and :math:`2.2 \times 10^{-15}`
      * - the distance to :math:`b = 0` formed as :math:`r_{k+1} - y_{\max}`;
          a cavity of radius :math:`10^{-12}`
        - qa, re-review (N2)
        - refused: the distance cancelled to 0, and the hp law refuses a
          singularity at distance 0
        - assembles and conserves; the distance is
          :math:`r_k^2/(r_{k+1} + y_{\max})` and :math:`b` is formed from
          :math:`(y_{\max} - y)(y_{\max} + y)`

   The archivist's own probe of 2026-10-07 (the table of
   :ref:`characteristic-wall-coupling`) reads the same regimes on the final
   code: the sphere's :math:`T_w` at :math:`\tau = 100` and 1000 to
   :math:`9.2 \times 10^{-13}` and :math:`6.1 \times 10^{-11}` at 16 points,
   the slab's at :math:`\tau = 30` to :math:`2.2 \times 10^{-15}` at 8.
   ERR-101 catalogues the fixed resolutions as one defect class and ERR-103
   the ungraded normal direction.

   **Why the margin study missed two of them.** Before the grazing depth
   was settled the main agent swept the margin, 0 to 20 extra halvings past
   :math:`\tau_{\min}`, and the totals moved by nothing, to
   :math:`5 \times 10^{-16}`, from :math:`\tau = 0.5` down to
   :math:`2.3 \times 10^{-12}`
   (``scratch/characteristic_architecture/p1_step_b3/margin/study.py``); no
   margin was adopted. Totals are integrals over every entry, and the entry
   qa found off by :math:`2.0 \times 10^{-6}` (N1) is invisible in them; and
   the sweep's thickest body was :math:`\tau = 2.3`, so the normal-direction
   layer of a thick slab (:math:`1.8 \times 10^{-4}` at :math:`\tau = 30`) was
   outside it. A study of totals over thin bodies certifies totals over thin
   bodies.

**A reference subtracted in floating point reports defects the code does
not have.** The near-void collision probability is :math:`1 - P_{\rm esc}`
with :math:`P_{\rm esc}` near 1. Evaluated in mpmath and subtracted in
double precision, it reported a defect of :math:`7.6 \times 10^{-7}`;
evaluated at 15 digits, a defect of 100 %. The subtraction belongs in
mpmath, and every near-void reference in the gates forms it there.

**Two limits that are not rule defects.** At
:math:`\Sigma_t = 10^{-315}` the block overflows (qa's F4), because the
flux, about :math:`1/\Sigma_t`, exceeds double precision: that is the
answer. A slab placed at
:math:`x \approx 10^{6}` loses digits to absolute positions (conservation
:math:`5.6 \times 10^{-10}`), filed as
`#585 <https://github.com/deOliveira-R/ORPHEUS/issues/585>`_.

The order of the lines, the chunks and the piece budget
-------------------------------------------------------

The traversal rule pads every line of a chunk to the chunk's largest
piece count, and a line's pieces grow as it nears grazing (the
exponential and hp gradings along it multiply). So ``Lines.ordered``
orders the lines by projected speed, and lines of like cost share a chunk.
``Lines.chunks`` takes the lines ``chunk`` at a time (the line rule's
default 512, the point rule's 4096), and a
chunk whose traversal rule holds more than ``budget`` piece slots
(:attr:`TraversalRule.extent
<orpheus.derivations.continuous.characteristic.transport.TraversalRule.extent>`,
the line rule's default 1024, the point rule's 32 768) is halved until it
fits or holds one line. The budget
bounds a chunk's memory where a line count does not, because a line's
memory follows its pieces. Each line's live pieces are packed first, so
the pieces pad only to the chunk's costliest line; since 2026-10-07 the
Volterra block's inner rule is packed the same way, interval by
interval, instead of padding every attenuated integral to the batch's
thickest stretch (:ref:`characteristic-cylinder-cost`).

`[M]` 2026-10-06, the main agent and the test-architect, before the inner
rule was packed:

.. list-table::
   :header-rows: 1
   :widths: 50 16 16 18

   * - Rule
     - Budget
     - Time
     - Peak memory
   * - 224 cylinder lines at :math:`\mu_z = 0.999`, chunked by count
     - none
     - --
     - 8.2 GB
   * - one slab rule of 66 560 slots
     - none
     - --
     - 36 GB
   * - the same, chunks under 5000 slots
     - 4096
     - --
     - 2.25 GB
   * - two-region white cylinder, 8 points, 8512 lines
     - 4096
     - 126 s
     - 7.9 GB
   * - the same
     - 256
     - 32 s
     - 0.6 GB
   * - one-region white cylinder at :math:`\tau = 0.01`, 8 points, 15 360
       lines graded to :math:`\theta = 3 \times 10^{-7}`
     - 256
     - over 25 min
     - --
   * - the same
     - 1024
     - 57 s
     - 1.3 GB
   * - the same
     - 4096
     - 72 s
     - 3.5 GB

Smaller is faster while the arrays stay in cache, until the chunks shrink
to a few lines and the per-rule overhead dominates: the grazing grading's
lines, each with many pieces, split a budget of 256 into single-line
rules. The default of 1024 is the measured compromise, and it stayed the
fastest after the inner rule was packed. `[M]` 2026-10-07, the main
agent, one process per run, one-region white cylinders at 8 points per
piece and 12 inner points per interval, per group:

.. list-table::
   :header-rows: 1
   :widths: 50 16 16 18

   * - Rule
     - Budget
     - Time
     - Peak memory
   * - :math:`\tau = 30`, 18 432 lines
     - 512
     - 29.4 s
     - 0.7 GB
   * - the same
     - 1024
     - 27.8 s
     - 0.8 GB
   * - the same
     - 4096
     - 31.4 s
     - 2.0 GB
   * - the same
     - 16 384
     - 37.2 s
     - 6.0 GB
   * - :math:`\tau = 0.01`, 15 360 lines
     - 512
     - 17.5 s
     - 0.5 GB
   * - the same
     - 1024
     - 15.6 s
     - 0.7 GB

The geometry of a
chunk (its chord and period) is rebuilt per group, as the sketch ruled;
a cache would follow only if the geometry were a measured fraction of the
time. On the sphere the cost is small: `[M]` 2026-10-07, the archivist's
probe, the three-region sphere at 16 points holds 480 lines and its white
block took 0.8 s.

.. _characteristic-cylinder-cost:

The cylinder's cost on a thick body (#586)
------------------------------------------

The cylinder's rule is built from the impact rule and the polar-angle
rule (iterated since #587, :ref:`characteristic-iterated-cylinder-rule`),
and both grow with the optical scale: the impact rule with the rim's
exponential ends, the polar rule with the grazing and normal gradings.
The measurements of this section were taken before the grading law, when
the rule was the two rules' tensor product. `[M]` 2026-10-07, the archivist's probe: the three-region
cylinder :math:`(0, 0.5, 1.5, 2)` at 8 points holds 38 400 lines, 240
impact nodes times 160 polar nodes; a homogeneous cylinder at
:math:`\tau = 30` holds 18 432 lines, 192 impact nodes times 96 polar
nodes, of 14 to 106 pieces each, about :math:`10^6` pieces.

**The premise that was refuted.** Issue
`#586 <https://github.com/deOliveira-R/ORPHEUS/issues/586>`_ measured the
:math:`\tau = 30` white block at 1117 s, on a machine at a load average of
about 6, and the three-region block at 723 s, conserving to
:math:`5.4 \times 10^{-13}`. It attributed the cost to the line count: in
a tensor product the polar grading that only the lines near grazing need
is applied at every impact parameter, and the rim grading at every polar
angle. It proposed a rule over :math:`(b, \theta)` that is not a tensor
product. A profile refuted that premise FOR the question "what makes the
cylinder slow". `[M]` 2026-10-07, the main agent, the white cylinders at
8 points over a sixteenth of their lines: 78 to 82 % of the time was spent
evaluating the panel basis
(:meth:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis.values`),
and 98 % of those calls came from the Volterra triangle's inner rule, the
attenuated integral of :ref:`characteristic-attenuated-integral`. The fact
it establishes is that the cost of each line was inflated, not that there
were too many lines.

The inner rule grades the stretch between an integral's start and its
stop exponentially toward the stop, so an integral over a thin stretch
has one live interval and one over a thick stretch up to eight. The rule
kept an interval if any integral of the batch used it, so every integral
was evaluated on the batch's thickest stretch, the thin ones on intervals
of zero width. The evaluations were 4.2 times (the three-region cylinder)
to 6.1 times (:math:`\tau = 30`) the live ones.

**The two fixes.** The user ruled on 2026-10-07 to fix the evaluation
first, re-measure, and change the line rule only if a block still took
more than a minute. Neither fix changes the rule.

1. **Packing.** Each line's live entries are packed first, in their
   order along the line; a per-entry field is gathered onto the packed
   entries, and the packed values are scattered back, summed onto their
   owners, by the transpose of the gather. The pieces and the inner
   rule's intervals share this one object (``_Packing`` in the transport
   module), so the owner of a packed entry is computed in one place. An
   array is padded to the chunk's costliest line, not to the batch's
   costliest integral, and the lines of a chunk are ordered by projected
   speed, a proxy for their cost, so the remaining padding is small.
2. **Lagrange tables.** The basis evaluation reads its Lagrange tables
   (the nodes and their differences, for the reference and the even
   panel coordinates), built once per basis, instead of rebuilding the node differences at
   every point. Each factor is the same arithmetic, so the values are
   bit-identical.

`[M]` 2026-10-07, the main agent and qa: on five fixtures (the
three-region and hollow spheres, a slab, the three-region and
:math:`\tau = 30` cylinders) the blocks are bit-identical to the code
before the fix except two, which differ by at most
:math:`2.7 \times 10^{-17}` relative. qa reproduced the old inner rule and showed that the difference
is the order of one summation: 192 of 192 calls are bit-identical once
that order is matched. The basis values are bit-identical on 140 140
entries, degrees 0 to 12. A mutation that shifts the packed owner by one
entry reddens the transport gates.

**The cost now.** `[M]` 2026-10-07, the main agent, one process per run,
budget 1024, white blocks per group:

.. list-table::
   :header-rows: 1
   :widths: 46 14 20 20

   * - Body
     - Points
     - Before
     - After
   * - one region, :math:`\tau = 30`
     - 8
     - about 215 s on an idle machine (1117 s under load)
     - 27.8 s
   * - one region, :math:`\tau = 0.01`
     - 8
     - 57 s
     - 15.6 s
   * - three regions :math:`(0, 0.5, 1.5, 2)`
     - 8
     - 723 s
     - 175 s
   * - one region, :math:`\tau = 100`
     - 8
     - --
     - 38 s
   * - one region, :math:`\tau = 30`
     - 16
     - --
     - 102 s
   * - one region, :math:`\tau = 100`
     - 16
     - --
     - 145 s
   * - three regions
     - 16
     - more than 15 min
     - 792 s

The characteristic gates' slow tier as it stood took 15 min 40 s for 14
rows, from about 62 minutes; the 722 rows outside it took 5 minutes.
#586's acceptance, a :math:`\tau = 30` two-group block in under a minute
at the gates' resolution, is met: two groups at 27.8 s each are 55.6 s.

**What a rule graded per impact panel would still buy.** The rule over
:math:`(b, \theta)` is iterated in one direction: each polar angle
carries its own impact rule (:ref:`characteristic-iterated-cylinder-rule`).
The other direction is not built: the polar rule is graded toward grazing
over the whole body's optical scale at every impact parameter. `[R]` A
polar rule graded on each impact panel by that panel's own optical scale
would skip the grazing grading where a short rim chord does not need it;
on the ``ABA`` cylinder the polar nodes below :math:`\sin\theta = 0.01`
carry 57 536 of the 102 320 lines (`[M]` 2026-10-09). Its saving
multiplies the packed cost; it does not repeat the packing's. It is the
lever if a three-region cylinder at 16 points, 792 s a group before the
grading law, enters a routine path
(`#587 <https://github.com/deOliveira-R/ORPHEUS/issues/587>`_). `[M]`
2026-10-09 (``scratch/characteristic_architecture/p1_step_c/polar_half_scales.py``)
refuted that lever on the ``ABA`` cylinder: the thinnest panel is the
outer-wall panel :math:`[1.96, 2.0]`, 0.02 and 0.04 optical in groups 0
and 1, every line crosses it, and so each impact panel's polar rule keeps
the body's 160 or 152 nodes. What stays open in #587 is a polar weight
absorbing :math:`\sin^2\theta`, unmeasured.

**The slow tier, re-ruled.** The third rung had cut the slow cylinder
rows to one fixture per law at 8 points, keeping the escape rows at
:math:`\tau = 2` and 8 at 16 points (where the 8-point transmission misses
by :math:`9.3 \times 10^{-10}` and :math:`2.3 \times 10^{-8}`) and moving
the thick legs from :math:`\tau = 30` and 100 to :math:`\tau = 8`. With the
packed rule the user ruled on 2026-10-07 to restore the cylinder's escape
and transmission legs at :math:`\tau = 30` and 100 at 16 points, and the
three-region cylinder's closed-body rows for both groups at 8 points.

`[M]` 2026-10-07, the test-architect: the characteristic slow tier holds
20 rows and runs in 36 min 27 s, against 14 rows in 15 min 40 s before the
restoration. The two thick legs miss the transmission :math:`T_w` by
:math:`2.6 \times 10^{-11}` at 16 points (band :math:`3 \times 10^{-10}`),
and both redden when the rim grading is removed, so they catch ERR-101.

The migration's re-pointing step needs two-group cylinder references;
`[R]` from the table, a three-region two-group block at 8 points costs
about 6 minutes.


.. _characteristic-galerkin-system:

The multigroup Galerkin system
==============================

The fourth rung turns the per-group transport blocks of the third rung
into the multigroup problem and answers its questions in nodal
coefficients:
:class:`~orpheus.derivations.continuous.characteristic.system.GalerkinSystem`
(:mod:`~orpheus.derivations.continuous.characteristic.system`) on the
cross sections of
:class:`~orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections`
(:mod:`~orpheus.derivations.continuous.characteristic.cross_sections`),
whose emission matrices come from
``group_emission`` (in ``orpheus.derivations.common.eigenvalue``). It gives
the fundamental mode and the higher modes of the k eigenproblem, their
adjoints, the flux of a fixed source and the adjoint flux of a detector.
The design is the plan's "P1 step (b), fourth rung: API sketch", ruled by
the user on 2026-10-07 (four questions, all as recommended), with the
unknown changed the same day from the flux to the emission density after
a measurement ("the flux-form adjoint refuted", below). The sources, the
detectors and the questions are nodal coefficients on the panel basis;
the question values of the interface vocabulary and the projection of a
mesh-free source onto the basis are the fifth rung's door
(:ref:`characteristic-door`).

The cross sections and the emission matrices
--------------------------------------------

Per region :math:`r` and for :math:`G` groups, the system reads three
arrays, each indexed ``[region, to, from]`` where it is a matrix:

.. math::

   \Sigma_t(r) \in \mathbb R^{G}, \qquad
   S(r) = \bigl(\Sigma_s(r) + 2\Sigma_2(r)\bigr)^{\mathsf T}, \qquad
   F(r) = \chi(r) \otimes \nu\Sigma_f(r),

with :math:`\Sigma_s` and :math:`\Sigma_2` the isotropic (Legendre order
0) scattering and (n,2n) transfer tables of the mixture, stored
``[from, to]``, so that the transpose puts them ``[to, from]``; the
factor 2 is the two neutrons an (n,2n) reaction emits. The two emission
matrices have one assembly site,
``group_emission``, which
returns them as a ``GroupEmission(scattering, fission)``. The
infinite-medium pair of :ref:`theory-homogeneous` reads the same
site, :math:`A = \mathrm{diag}(\Sigma_t) - S` and :math:`F`, so the
characteristic reference and the 0-D reference cannot disagree on the
emission; the re-expression of ``_infinite_medium_matrices`` keeps its
three floating-point operations in their order (one add, one transpose,
one subtraction from the diagonal), so the 0-D values are bitwise
unchanged. The (n,2n) multiplicity literal moved with the emission, and
the census of that literal
(``tests/gates/transport/test_n2n_multiplicity_census.py``,
``_REFERENCE_LITERALS``) names ``group_emission`` as the reference site.

**Why the emission is read on its own.** The transport block of a group
already carries :math:`\Sigma_t` (:ref:`characteristic-galerkin-assembly-section`),
so the system needs the emission alone. Reading it back out of the 0-D
loss matrix by subtracting the diagonal would form :math:`\Sigma_t -
(\Sigma_t - \Sigma_{s,gg})`, a cancellation, and would tie the system to
the loss matrix's spelling; the user ruled on 2026-10-06 (Q4 of the third
rung) that the emission is split out instead.

**The door is the mixtures.**
:meth:`RegionCrossSections.of
<orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.of>`
reads one :class:`~orpheus.data.macro_xs.mixture.Mixture` per region
(``SigT``, ``SigS[0]``, ``Sig2[0]``, ``SigP`` for :math:`\nu\Sigma_f`,
``chi``) through ``group_emission``. The mixture is the input boundary
that checks the data's physics (the balance, the fission spectrum); the
direct constructor checks shapes, finiteness and signs only, and stores
read-only arrays. ``of`` refuses a mixture whose ``SigS`` or ``Sig2``
has a non-zero block at a Legendre order :math:`\ge 1`, naming the region
and the order (``NotImplementedError``): anisotropic emission needs an
angular basis on each line, which the reference does not have. The
refusal is a declared scope boundary read at the mixtures, the one door
(the user's ruling of 2026-10-07, Q2 of the fourth rung); a stack
carrying an all-zero higher block is admitted.

The emission support, per group and exact
-----------------------------------------

The emission density of group :math:`g` can be non-zero only where
something is emitted into :math:`g`. The **emission support** of group
:math:`g` is the set of regions

.. math::

   \mathrm{supp}_g \;=\; \bigl\{\, r : S(r)_{g,\cdot} \ne 0
   \ \text{or}\ F(r)_{g,\cdot} \ne 0 \,\bigr\}
   \;\cup\; \{\, r : \text{a source of group } g \text{ is posed in } r \,\},

read by
:meth:`RegionCrossSections.emission_support
<orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.emission_support>`
(a region-major ``(n, G)`` mask: a non-zero row ``scattering[r, g, :]``
or ``fission[r, g, :]``) and widened by the system's ``source_regions``,
a mask of the same orientation. The columns of group :math:`g`'s block are
the basis functions of :math:`\mathrm{supp}_g`, and the column sets differ
by group. A fissile region whose :math:`\chi` is zero in a group emits in
that group only by scattering: the fission row is zero there, and the
support is read from the scattering. This is the user's ruling of
2026-10-07 (Q1 of the fourth rung), over two alternatives that both fail
on a body the gates pose:

- **Each group's** :math:`\Sigma_t > 0` (the third rung's default) drops a
  region transparent in :math:`g` that receives emission into :math:`g`:
  the region has no column, and its emission is lost without a message.
  `[M]` 2026-10-07, the test-architect (``ta/m4_support.py``): an inner
  region transparent in group 1 that scatters from group 0 into it,
  inside an absorbing, scattering shell, with a group-0 source; the
  absorption balances the source to :math:`-2.2 \times 10^{-2}` on the
  mirror sphere and :math:`-2.0 \times 10^{-1}` on the slab between
  mirrors under the :math:`\Sigma_t > 0` support, and to
  :math:`4.9 \times 10^{-15}`, :math:`2.2 \times 10^{-16}` and
  :math:`5.1 \times 10^{-15}` (the mirror sphere, the slab, the white
  sphere) under the exact one.
- **The union over groups** assembles columns whose emission is
  identically zero, and on a region void in :math:`g` behind a mirror such
  a column is a source on lossless trapped lines
  (:ref:`characteristic-closure-section`): the block is refused with
  :class:`~orpheus.derivations.continuous.characteristic.closure.TrappedSource`
  although the problem is well posed. `[M]` 2026-10-07, the archivist's
  probe (``archivist_rung4_numbers.py``, the session scratchpad), the same
  sphere with an outer layer void in group 1 that emits in group 0 only,
  under a mirror: the exact support (``[[1, 1], [1, 1], [1, 0]]``, 52 and
  32 coefficients) is answered, and the union support is refused with
  "a source on a lossless trapped line". Only the mirror sphere
  discriminates: the lines whose impact parameter exceeds the material's
  radius stay in the layer, lossless in group 1; on the slab and on the
  white sphere a column there is not refused (the test-architect's
  ``m4_support``).

The default of :meth:`LineRule.transport
<orpheus.derivations.continuous.characteristic.assembly.LineRule.transport>`,
:math:`\Sigma_t > 0`, serves a block assembled on its own, where no
emission is known; the system always passes each group's exact support.

**A field off the support is refused.** The emission space has no
coefficient there, so a source with a non-zero coefficient outside its
group's support would be dropped:
:meth:`EmissionSpace.restrict
<orpheus.derivations.continuous.characteristic.system.EmissionSpace.restrict>`
refuses it, naming the regions by group, and the caller poses the system
with those regions in ``source_regions``. The mask is region-major,
as the cross sections; a group-major mask is refused by its shape when
:math:`n \ne G`, and with :math:`n = G` no shape can see a transposed
mask, which is why the orientation is one convention for every
``(n, G)`` array of the package (the elegance review's finding C1, which
measured the transposed mask posing the wrong regions with
:math:`n = G = 2`).

The weak form on the emission density, derived
----------------------------------------------

Write :math:`\mathcal K_g` for the scalar-flux operator of an isotropic
emission in group :math:`g`, the operator whose Galerkin matrix is the
transport block (:eq:`characteristic-galerkin-assembly`), and
:math:`E(k) = S + F/k` for the emission per unit flux, acting point by
point with the region's matrices. The multigroup eigenproblem with an
external source :math:`q^{\rm ext}` is

.. math::

   \phi_g = \mathcal K_g q_g, \qquad
   q_g = \sum_{g'} S_{g \leftarrow g'}\,\phi_{g'}
   + \frac{1}{k}\sum_{g'} F_{g \leftarrow g'}\,\phi_{g'} + q^{\rm ext}_g ,

and eliminating the flux leaves an equation for the emission density
alone,

.. math::

   q \;=\; E(k)\,\mathcal K q \;+\; q^{\rm ext}.

**The trial space.** The emission of group :math:`g` vanishes off
:math:`\mathrm{supp}_g`, so it is expanded in the basis functions there,
:math:`q_g = \sum_{j \in \mathrm{supp}_g} q_{g,j}\,u_j` (``M_g`` of them;
:math:`M = \sum_g M_g` in all). The flux is expanded on every basis
function, :math:`\phi_g = \sum_i \phi_{g,i}\,u_i`, and fixed by testing
:math:`\phi_g = \mathcal K_g q_g` against every :math:`u_i`:

.. math::

   \sum_{i'} W_{i i'}\,\phi_{g,i'} \;=\; \sum_{j \in \mathrm{supp}_g} K_g[i, j]\, q_{g,j},
   \qquad\text{that is}\qquad W_G\,\phi = K q ,

with :math:`W` the mass matrix of the basis, :math:`W_G = I_G \otimes W`,
and :math:`K = \mathrm{blockdiag}(K_g)` of shape :math:`GN \times M`
(:meth:`GalerkinSystem.flux
<orpheus.derivations.continuous.characteristic.system.GalerkinSystem.flux>`).

**Testing the emission equation.** Test
:math:`q = E\mathcal K q + q^{\rm ext}` against :math:`u_i` for
:math:`i \in \mathrm{supp}_g`. The left side gives
:math:`\sum_{j} W_{ij}\,q_{g,j}`, the mass matrix restricted to the
support. On the right, :math:`u_i` lives on one panel, no panel crosses a
region, and the cross sections are constant on a region, so
:math:`E_{g g'}` is a constant on the support of :math:`u_i` and leaves the
integral:

.. math::

   \bigl\langle u_i, (E\,\mathcal K q)_g \bigr\rangle
   \;=\; \sum_{g'} E_{g g'}\bigl(r(i)\bigr)\,\bigl\langle u_i, \mathcal K_{g'} q_{g'} \bigr\rangle
   \;=\; \sum_{g'} E_{g g'}\bigl(r(i)\bigr) \sum_{j \in \mathrm{supp}_{g'}} K_{g'}[i, j]\, q_{g', j},

with :math:`r(i)` the region of node :math:`i`. The step is exact: the
emission acts node by node and adds no projection error. The row
:math:`i` of :math:`K_{g'}` is read whether or not :math:`i` lies in
:math:`\mathrm{supp}_{g'}`, which is why the blocks' rows span every
panel. Collecting the rows gives the Galerkin system of the emission:

.. math::
   :label: characteristic-pencil

   W_s\, q \;=\; \Bigl(S + \frac{1}{k}F\Bigr) K q \;+\; W_s\, q^{\rm ext},
   \qquad
   W_G\,\phi \;=\; K q ,
   \qquad
   W_s = R\,W_G\,R^{\mathsf T},
   \quad
   S = R\,E_{\rm node}(\Sigma_s + 2\Sigma_2),
   \quad
   F = R\,E_{\rm node}(\chi \otimes \nu\Sigma_f),

where :math:`R` (:math:`M \times GN`, a 0/1 selection) restricts a field to
the emission space, group-major, and :math:`E_{\rm node}(e)` is the
:math:`GN \times GN` matrix whose block :math:`(g, g')` is
:math:`\mathrm{diag}_i\, e[r(i), g, g']`: the per-region table placed at
every node. :math:`S` and :math:`F` are :math:`M \times GN`. The k
eigenproblem is the pencil :math:`(W_s - SK,\; FK)` of the source-free
equation, :math:`F K q = k\,(W_s - S K)\,q`
(:attr:`GalerkinSystem.pencil
<orpheus.derivations.continuous.characteristic.system.GalerkinSystem.pencil>`,
on :class:`~orpheus.derivations.common.dense_pencil.DensePencil`): its
fundamental mode, its higher modes (``pencil.spectrum()``) and, by
transposition, their adjoints.

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.pencil

   **Implemented by** the system's matrices: ``emission_mass``
   (:math:`W_s`), ``transport`` (:math:`K`), ``scattering`` and
   ``fission`` (:math:`S`, :math:`F`, through ``_per_node``), ``flux``
   (:math:`W_G\phi = Kq`), the restriction
   ``EmissionSpace.restriction`` (:math:`R`), and the per-region tables of
   ``RegionCrossSections.of`` through ``group_emission``.

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.emission_mass

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.transport

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.scattering

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.fission

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.flux

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.system.EmissionSpace.restriction

.. implements:: characteristic-pencil
   :by: orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.of

.. implements:: characteristic-pencil
   :by: orpheus.derivations.common.eigenvalue.group_emission

**The flux form has the same eigenvalues and the same flux.** Writing the
unknown as the flux instead, :math:`W_G\phi = K(S + F/k)\phi + K q^{\rm ext}`,
the pencil :math:`(W_G - KS,\; KF)` of size :math:`GN`, is the same
problem: put :math:`q = (S + F/k)\phi` and use that :math:`E_{\rm node}`
commutes with :math:`W` panel by panel (both act within a panel, and
:math:`E` is constant there), so :math:`W_s (S + F/k)\phi = (S + F/k) W_G
\phi`. `[M]` 2026-10-07, the main agent's ``main/qform.py``, re-run by the
archivist on the built code (a two-group, two-region mirror sphere with
upscatter, :math:`\Sigma_t = (1, 2)`, degree 3, 8 line points, 12 along
each line): the two forms' :math:`k` agree to :math:`4.4 \times 10^{-16}`
and the flux the emission form gives, :math:`W_G^{-1} K q`, equals the
flux form's mode to :math:`3.2 \times 10^{-15}`. The emission form is
smaller (:math:`M \le GN`; equal on a closed body where every group emits
everywhere). Where the two differ is the adjoint.

Why the unknown is the emission: the adjoint
--------------------------------------------

:meth:`DensePencil.adjoint
<orpheus.derivations.common.dense_pencil.DensePencil.adjoint>` returns the
transposed pair, which is the Galerkin matrix of the adjoint of whatever
operator equation the pencil poses (:ref:`verification-reference-kernel`).
So the question is which equation is posed. The transport operator of an
isotropic emission is self-adjoint by reciprocity, :math:`\mathcal K^\ast =
\mathcal K` in :math:`L^2`, and the emission's adjoint is its transpose
per point, :math:`E^\ast = S^{\mathsf T} + F^{\mathsf T}/k`.

- **The flux form** poses :math:`\phi = \mathcal K E\,\phi`. Its adjoint is
  :math:`\psi = E^\ast \mathcal K\,\psi`, and that is not the adjoint flux.
  The adjoint flux :math:`\phi^\dagger` (the importance) solves
  :math:`\phi^\dagger = \mathcal K E^\ast \phi^\dagger`; applying
  :math:`E^\ast` gives :math:`E^\ast\phi^\dagger = E^\ast \mathcal K
  (E^\ast\phi^\dagger)`, so the flux form's adjoint eigenvector is
  :math:`\psi = E^\ast\phi^\dagger`, the adjoint collision density. In an
  infinite medium :math:`A^{\mathsf T}\phi^\dagger = F^{\mathsf
  T}\phi^\dagger/k` with :math:`A = \mathrm{diag}(\Sigma_t) - S`, so
  :math:`E^\ast\phi^\dagger = \Sigma_t\,\phi^\dagger`: the transpose reads
  the importance multiplied group by group by :math:`\Sigma_t`.
- **The emission form** poses :math:`q = E\,\mathcal K q`. Its adjoint is
  :math:`\phi^\dagger = \mathcal K^\ast E^\ast \phi^\dagger = \mathcal K
  E^\ast \phi^\dagger`: the adjoint flux itself.

.. math::
   :label: characteristic-adjoint

   \phi^\dagger = \mathcal K E^\ast \phi^\dagger
   \quad\Longrightarrow\quad
   W_s\,\phi^\dagger \;=\; K^{\mathsf T}\Bigl(S + \frac{1}{k}F\Bigr)^{\mathsf T}\phi^\dagger ,

the transposed pencil :math:`\bigl((W_s - SK)^{\mathsf T}, (FK)^{\mathsf
T}\bigr)` = ``pencil.adjoint()``, its unknown the coefficients of
:math:`\phi^\dagger` on the emission space.

.. implements:: characteristic-adjoint
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.pencil

   **Implemented by** ``GalerkinSystem.pencil`` read through the kernel's
   ``DensePencil.adjoint``, and by ``GalerkinSystem.response``, the
   adjoint of the source pencil.

.. implements:: characteristic-adjoint
   :by: orpheus.derivations.common.dense_pencil.DensePencil.adjoint

.. implements:: characteristic-adjoint
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.response

**Why the transposed matrix is that Galerkin form.** Transpose the right
side of :eq:`characteristic-pencil`:
:math:`\bigl((S + F/k) K\bigr)^{\mathsf T} = K^{\mathsf T} (S + F/k)^{\mathsf T}`.

1. :math:`(S + F/k)^{\mathsf T}` takes coefficients on the emission
   space to a field on every node: at node :math:`i` of group :math:`g'`
   it is :math:`\sum_g E_{g g'}(r(i))\,\phi^\dagger_{g,i}`, the adjoint
   emission :math:`E^\ast\phi^\dagger` node by node (zero where no group
   has a coefficient).
2. :math:`K^{\mathsf T}` takes a field :math:`y` on every node to
   :math:`\sum_i K_{g'}[i, b]\,y_{g',i}` at the support node :math:`b` of
   group :math:`g'`. By reciprocity :math:`K_{g'}[i, b] = \langle u_i,
   \mathcal K u_b\rangle = \langle u_b, \mathcal K u_i\rangle`, so this is
   :math:`\langle u_b, \mathcal K_{g'} y_{g'}\rangle`, the Galerkin test of
   :math:`\mathcal K y` against the support's functions. For a node
   :math:`i` off :math:`\mathrm{supp}_{g'}` the entry
   :math:`\langle u_b, \mathcal K u_i\rangle` is a column the support never
   assembles, and its value is the row entry :math:`K_{g'}[i, b]`, which
   the full rows do assemble.

So the transpose is the Galerkin matrix of :math:`\mathcal K E^\ast` on
the same basis, tested and expanded on the supports, and
:math:`W_s` is symmetric: no metric enters. The adjoint flux is defined
exactly on the supports, which is where a source can be posed (Q1), so it
is the importance of every source the system admits.

`[M]` 2026-10-07, the evidence that ruled the form (the main agent's
``main/smoke.py`` and ``qform.py``, re-run by the archivist on the built
code, the mirror sphere above, :math:`\chi = (1, 0)`): the flux form's
transposed eigenvector is flat with the group ratio 2.146567717996, the
emission form's 1.07328385899813, and the infinite medium's adjoint
:math:`A^{-\mathsf T}\nu\Sigma_f` 1.07328385899814; the factor between
the first two is :math:`\Sigma_{t,2}/\Sigma_{t,1} = 2`. qa's independent
route (``qa/p1_adjoint_independent.py`` and
``qa/p1b_adjoint_disjoint_support.py``, re-run by the archivist): the
adjoint problem of an isotropic multigroup problem is the forward problem
with each region's transfer matrices transposed (same :math:`\Sigma_t`,
same blocks), so ``response(r)`` was compared with the forward fixed
source ``r`` of the group-transposed problem posed on every region, on a
three-region two-group sphere with up- and downscatter, (n,2n) and
:math:`\chi = (0.7, 0.3)`, behind vacuum and behind a white wall of
amplitude 0.6, two support patterns: 12 of 12 trials agree to
:math:`3.4 \times 10^{-15}` to :math:`8.3 \times 10^{-15}`, relative; the
adjoint's fundamental eigenvalue equals the transposed problem's
:math:`k` to the last digits and its mode equals that problem's flux to
:math:`5.4 \times 10^{-15}`.

Two splittings of one operator
------------------------------

The source-free operator of :eq:`characteristic-pencil` is
:math:`W_s - (S + F/k) K`, and the system splits it two ways:

- **the k pencil** :math:`(W_s - SK,\; FK)`, the fission alone scaled by
  :math:`1/k`: the eigen questions, ``pencil.fundamental()``,
  ``pencil.spectrum()`` and ``pencil.adjoint()``;
- **the source pencil** :math:`(W_s,\; (S + F)K)`, at :math:`k = 1` with
  every secondary emission in the gain
  (:attr:`GalerkinSystem.source_pencil
  <orpheus.derivations.continuous.characteristic.system.GalerkinSystem.source_pencil>`):
  the source questions.

.. math::
   :label: characteristic-fixed-source

   \begin{aligned}
   W_s\, q &= (S + F)\,K q + W_s\, q^{\rm ext}, \qquad W_G\,\phi = K q,\\
   W_s\, \phi^\dagger &= K^{\mathsf T} (S + F)^{\mathsf T} \phi^\dagger + K^{\mathsf T} r,
   \qquad
   \langle r, \phi(q^{\rm ext}) \rangle = r^{\mathsf T} W_G\,\phi
   = \phi^{\dagger\mathsf T} W_s\, q^{\rm ext}.
   \end{aligned}

:meth:`GalerkinSystem.fixed_source
<orpheus.derivations.continuous.characteristic.system.GalerkinSystem.fixed_source>`
restricts the source's nodal coefficients ``(G, N)`` to the emission
space (refusing a coefficient off the support), solves the first line by
the source pencil's least solution and returns the flux ``(G, N)``.
:meth:`GalerkinSystem.response
<orpheus.derivations.continuous.characteristic.system.GalerkinSystem.response>`
takes a detector's nodal coefficients :math:`r` and returns the adjoint
flux on the emission space, the least solution of the second line with
the load :math:`K^{\mathsf T} r`, which is :math:`\langle u_b, \mathcal K
r\rangle` by the reciprocity argument above (the load is not
:math:`K^{\mathsf T} W_G r`: the detector's coefficients are paired with
the flux through :math:`W_G` once, in the reading). The reading follows
from the first line: with :math:`B = W_s - (S + F)K`,
:math:`r^{\mathsf T} W_G \phi = r^{\mathsf T} K q = r^{\mathsf T} K B^{-1}
W_s q^{\rm ext} = (B^{-\mathsf T} K^{\mathsf T} r)^{\mathsf T} W_s
q^{\rm ext}`, and :math:`B^{-\mathsf T} K^{\mathsf T} r` is
:math:`\phi^\dagger`.

.. implements:: characteristic-fixed-source
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.source_pencil

   **Implemented by** ``source_pencil``, ``fixed_source``, ``response`` and
   the emission space's ``restrict``.

.. implements:: characteristic-fixed-source
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.fixed_source

.. implements:: characteristic-fixed-source
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.response

.. implements:: characteristic-fixed-source
   :by: orpheus.derivations.continuous.characteristic.system.EmissionSpace.restrict

**Why the source questions use the source pencil.**
:meth:`DensePencil.least_solution
<orpheus.derivations.common.dense_pencil.DensePencil.least_solution>` sums
the Neumann series over collisions and refuses a gain whose spectral
radius on the reached unknowns is not below 1
(:class:`~orpheus.derivations.common.dense_pencil.NoLeastSolution`). That
check guards the physics only if every secondary emission is in the
gain. On the k pencil the gain is the fission alone, so a closed body made
supercritical by its (n,2n) emission with no fission has a zero gain
there, and the k pencil's solve answers it directly, with a finite and
negative flux. On the source pencil the same body is refused. `[M]`
2026-10-07, the archivist's probe on the gates' fixture (no fission;
:math:`\Sigma_t = (1, 1.5)`, an (n,2n) transfer of 0.35 and 0.4 within
each group; the infinite medium's collision gain has spectral radius
1.155), the three-region mirror sphere: ``fixed_source`` raises
``NoLeastSolution`` with the radius 1.15485837704 on the 104 unknowns the
source reaches, while the direct solve of the k pencil's loss on the same
load returns coefficients down to :math:`-30.0`.

**The two pencils read one system.** On a subcritical body, the source
:math:`q^{\rm ext} = (1/k - 1)\,F\phi_k` returns the fundamental flux
:math:`\phi_k`: the source pencil's emission is then :math:`(S + F)\phi_k
+ (1/k - 1)F\phi_k = (S + F/k)\phi_k`, the mode's own (it needs
:math:`k < 1`, for the source to be non-negative and the source pencil
subcritical). `[M]` 2026-10-07, the test-architect (``ta/m9_source``): to
:math:`3.4 \times 10^{-14}` of the flux's maximum on a vacuum sphere
(:math:`k = 0.0753`) and :math:`3.0 \times 10^{-15}` on a slab
(:math:`k = 0.187`). The parent verification spec's "the source
:math:`F\phi_k/k` returns :math:`\phi_k`" held for a source pencil without
the fission in its gain; with it, the source is the fission deficit.

The one-group Rayleigh–Ritz bound
---------------------------------

In one group the Galerkin :math:`k` is a lower bound on the exact
:math:`k` and increases on nested trial spaces. The derivation:

1. In one group :math:`E = M(k) = \sigma_s + \nu\sigma_f/k`, constant on
   each region, and positive exactly on the emission support (it is the
   support's definition). Write :math:`M` also for the diagonal of its
   values at the support nodes. :eq:`characteristic-pencil` reads
   :math:`W_s q = M K_{ss} q` with :math:`K_{ss}` the block's rows on the
   support (the restriction :math:`R` picks them).
2. :math:`W_s` is block diagonal by panel (the basis is discontinuous) and
   :math:`M` is constant on each panel, so they commute, and
   :math:`W_s M^{-1} q = K_{ss} q`: a symmetric generalised
   eigenproblem, :math:`K_{ss}` symmetric by reciprocity and
   :math:`W_s M^{-1}` symmetric positive definite.
3. This is the Rayleigh–Ritz (Galerkin) form, on the trial space of the
   support's basis functions, of the continuous problem
   :math:`\mathcal K q = \lambda\,M^{-1} q` for the self-adjoint positive
   operator :math:`\mathcal K` in the weight :math:`M^{-1}`:

   .. math::
      :label: characteristic-one-group-bound

      \lambda_h(k) \;=\; \max_{q \in V_h} \frac{\langle q, \mathcal K q\rangle}{\langle q, M(k)^{-1} q\rangle}
      \;\le\; \max_{q} \frac{\langle q, \mathcal K q\rangle}{\langle q, M(k)^{-1} q\rangle} \;=\; \lambda(k),
      \qquad
      V_h \subset V_{h'} \;\Rightarrow\; \lambda_h(k) \le \lambda_{h'}(k).

4. The critical condition is :math:`\lambda(k) = 1`. :math:`M(k)` decreases
   as :math:`k` grows wherever :math:`\nu\sigma_f > 0`, so
   :math:`\lambda(k)` decreases. Then :math:`\lambda(k_h) \ge
   \lambda_h(k_h) = 1 = \lambda(k)` gives :math:`k_h \le k`, and nesting
   gives :math:`k_h \le k_{h'}`.
5. The flux off the support is read from the emission
   (:math:`W_G\phi = Kq`) and does not enter the eigenproblem, so it does
   not enter the bound.

The bound is a property of the Galerkin form with exact integrals. Two
errors compete with its margin: the line rule's integration error, which
can raise :math:`k`, and, against a published truth, the resolution of
that truth's printed digits. It holds in one group only: with two groups
the emission couples the groups through a non-symmetric matrix and no
weight makes the form symmetric. The 1-group degeneracy of
``vv-principles`` applies in full: in one group :math:`k` cannot detect an
error in the scattering's group structure, so these rows are theorem rows
on the Galerkin form, and the multigroup claims rest on the two-group
rows (the closed bodies and Sood's critical sizes, below).

**Nested spaces.** On fixed panels, the degrees :math:`p = 1, 2, 3, 4, 6`
nest (the even panel too: polynomials in :math:`c^2` of degree :math:`p`
lie in those of degree :math:`p + 1`). `[M]` 2026-10-07, the
test-architect (``ta/m5_rr_real``, ``ta/m6_floor``), a heterogeneous
one-group body (two fissile regions round a scatterer) as a vacuum
sphere, a mirror–vacuum slab and a vacuum hollow sphere at 2 grading
layers: the smallest increment of :math:`k` is :math:`2.8 \times 10^{-8}`,
and the line rule's own floor, :math:`k` at 8 line points and 12 along
each line against 16 and 16, is at most :math:`1.7 \times 10^{-13}`. An
under-integrated rule (4 line points, 8 along each line) moves :math:`k`
at :math:`p = 4` by :math:`+1.2 \times 10^{-5}` and :math:`+1.4 \times
10^{-5}` on the spheres and breaks the order.

**Against an independent truth.** Sood's one-group bare critical sizes
:cite:`SoodForsterParsons2003` (Ua-1-0-SP, a sphere of 2.4248249802 mean
free paths; Ua-1-0-SL, a half slab of 0.93772556 behind a central
mirror) are bodies whose exact :math:`k` is 1. The truth resolves
:math:`k` to :math:`|\mathrm d k/\mathrm d(\mathrm{mfp})|` times half a
unit of its last printed digit: :math:`1.7 \times 10^{-11}` for the sphere
and :math:`3.6 \times 10^{-9}` for the slab (:math:`\mathrm d k/\mathrm
d(\mathrm{mfp}) = 0.726`). `[M]` 2026-10-07, the archivist's probe, 8 line
points and 12 along each line, :math:`k - 1` at degrees 1 to 4:

.. list-table::
   :header-rows: 1
   :widths: 34 16 16 16 18

   * - Body, panels
     - :math:`p = 1`
     - :math:`p = 2`
     - :math:`p = 3`
     - :math:`p = 4`
   * - Ua-1-0-SP, one panel
     - :math:`-9.889 \times 10^{-4}`
     - :math:`-4.204 \times 10^{-5}`
     - :math:`-6.706 \times 10^{-6}`
     - :math:`-2.385 \times 10^{-6}`
   * - Ua-1-0-SL, one panel
     - :math:`-4.408 \times 10^{-3}`
     - :math:`-2.380 \times 10^{-5}`
     - :math:`-1.095 \times 10^{-5}`
     - :math:`-3.836 \times 10^{-6}`
   * - Ua-1-0-SP, 2 layers, ratio 0.4
     - :math:`-2.520 \times 10^{-5}`
     - :math:`-2.142 \times 10^{-6}`
     - :math:`-4.942 \times 10^{-7}`
     - :math:`-1.407 \times 10^{-7}`
   * - Ua-1-0-SL, 2 layers, ratio 0.4
     - :math:`-3.411 \times 10^{-4}`
     - :math:`-8.791 \times 10^{-8}`
     - :math:`-1.681 \times 10^{-8}`
     - :math:`+4.275 \times 10^{-10}`

On the ungraded panel every margin is at least :math:`2.4 \times 10^{-6}`,
a thousand times the truth's resolution, and the sequences increase. On
the graded panels the slab's degree 4 reads :math:`+4.3 \times 10^{-10}`:
the graded basis has converged below the slab truth's resolution of
:math:`3.6 \times 10^{-9}`, so the reading is inside the truth's
uncertainty, not a violation, and it has no margin to show. The gate is
posed on the ungraded panel for that reason.

The emission space and the posing
---------------------------------

:class:`~orpheus.derivations.continuous.characteristic.system.EmissionSpace`
is the emission space as an object: each group's support (the basis
indices of its nodes, increasing) and each node's region. Its restriction
:math:`R` is the one spelling of the group-major ragged layout, and the
system's matrices are compositions with it: :math:`W_s = R W_G
R^{\mathsf T}`, :math:`S = R\,E_{\rm node}(\cdot)`, :math:`F = R\,E_{\rm
node}(\cdot)`. ``restrict`` is the parse of a field into it (refusing a
coefficient off the support) and ``split`` cuts a vector on it into one
array per group. A pencil's eigenvector and ``response`` live on this
space; read the flux with ``flux`` and a group's coefficients with
``split``.

**The system is posed, never handed its blocks.** Its fields are the
problem and the resolution: the basis, the walls, the cross sections, the
transport resolution
(:class:`~orpheus.derivations.continuous.characteristic.assembly.TransportResolution`:
the line rule's points per piece and the two rules along each line) and
the source regions. The transport blocks are derived from them
(``GalerkinSystem.groups``, one ``LineRule.of`` and one ``transport`` per
group, each group's rule graded by its own :math:`\Sigma_t`, as the third
rung ruled), so a block assembled for another group or another cross
section cannot be put in.

**Cost.** The dense problem has :math:`M` unknowns, a few hundred at the
gates' resolution, and its eigen and least solves take under 0.03 s
(`[M]` 2026-10-07, the test-architect's specification §4). The blocks
cost what the third rung measured, per group
(:ref:`characteristic-cylinder-cost`): a sphere system of two or three
groups builds in 0.07 to 0.9 s, a slab in 0.7 to 4.1 s, the one-region
cylinder in 26 s for two groups at 8 points.


.. _characteristic-door:

The door: the reference posed from a specification
==================================================

The fifth rung poses the Galerkin system of the fourth rung from the
interface vocabulary and answers the observables of its question (its
first half, 5a; the point value is its second half, 5b,
:ref:`characteristic-reading`):
:class:`~orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation`
(:mod:`~orpheus.derivations.continuous.characteristic.reference`), built
from a
:class:`~orpheus.specification.specification.GeometrySpecification` (the
materials, a finite concentric geometry and one question) and a
:class:`~orpheus.derivations.continuous.characteristic.reference.Resolution`.
The factory
:func:`~orpheus.derivations.continuous.characteristic.reference.characteristic_reference`
wraps it in a
:class:`~orpheus.reference.solution.ReferenceSolution` with no
certificate, the consumer surface the trajectory-resolvent family exposed
until its deletion (:ref:`characteristic-origins-history-reading`). Every reading is
:class:`~orpheus.reference.reading.Uncertified`: the family derives no
error bound yet (#566).

The design is the plan's "P1 step (b), fifth rung: API sketch"
(``.claude/plans/characteristic_reference_architecture.md``), ruled by the
user on 2026-10-07 in four questions, all as recommended: the rung splits
into 5a (the resolution, the projection, the door and the two
observables that need no point, ``Eigenvalue`` and ``FluxIntegral``) and
5b (the reading at a point); a ``Response`` is answered as the
group-transposed forward problem; the eigen flux is gauged by its
production; ``Nearest(tau)`` is served with :math:`\tau` read in the
:math:`k` chart. After the reviews (qa, the elegance review and the
test-architect's measurements, ``scratch/characteristic_architecture/p1_step_b5/``)
the user ruled three corrections the same day: a ``Nearest`` answer reads
its eigenvalue only; the gauge is a total production of 1 over the body,
not the production density 100 the sketch had named; and a ``Response``
answers the adjoint scalar flux :math:`R\psi^\dagger`, so the detector's
role arrow carries a :math:`4\pi`. On 2026-10-08 the user ruled which
production: the one the question declares
(``Eigen.gauge``), by default
what fission and the (n,2n) reaction emit, the functional the
S\ :sub:`N` solver scales by. Each is derived below.

The objects and their roles
---------------------------

**The resolution is one value.**
:class:`~orpheus.derivations.continuous.characteristic.reference.Resolution`
holds the panel basis's degree :math:`p`, its grading (``layers`` and
``ratio``, read by :meth:`PanelBasis.of
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.of>`),
the transport blocks' rules (a
:class:`~orpheus.derivations.continuous.characteristic.assembly.TransportResolution`)
and ``source_points``, the Gauss points per piece that project a
symbolic function onto the basis. It refuses a non-integer count and
``source_points`` below :math:`2(p + 1)`, the count at which the
projection integrates a function of the panel space as exactly as the
mass matrix does (the projection, below). The solve and every reading use
this one value, so a reading cannot be taken at a resolution other than
the solve's. The reading at a point adds no field: its rule over the
lines through the point takes the transport resolution's ``line_points``,
and its traversals the same ``points`` and ``inner_points`` as the blocks
(the user's ruling of 2026-10-08). ``Resolution`` and ``TransportResolution`` are
:class:`~orpheus.numerics.content.ContentIdentity` values, as the
derivation is: its content is its two fields, and the cross sections,
the basis and the system are derived from them.

**The door is construction.** ``CharacteristicDerivation.__post_init__``
is the one place the posing is decided, and it refuses, before anything
is assembled:

- an infinite medium, which has no geometry and so no lines
  (``TypeError``);
- a question at a point other than the physical one (``ValueError``);
- an ``Eigen`` question along any parameter but the fission emission,
  ``CellCoefficient.every(Channel.FISSION_EMISSION)`` resolved on the
  materials (``ValueError``);
- a ``Symbolic`` source or detector that depends on the direction
  :math:`\mu` or :math:`\varphi` (``NotImplementedError``, a declared
  scope boundary: an anisotropic source needs its own first-flight
  transport along each line, the machinery an anisotropic emission
  needs too);
- through the objects it builds, a mixture with anisotropic emission
  (:meth:`RegionCrossSections.of
  <orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.of>`),
  a wall law the reference does not read (:meth:`Walls.of
  <orpheus.derivations.continuous.characteristic.walls.Walls.of>`) and a
  grading the basis refuses.

The cross sections are read with one mixture per interval, in the
geometry's interval order (``geometry.mat_ids``), and the basis on
``ConcentricPartition.of(geometry)``. Construction assembles no
transport block: the system's blocks are derived on first use, and the
gate ``test_a_served_question_constructs_and_assembles_nothing`` counts
zero calls of ``LineRule.of`` for every served question.

**The answer is solved once, typed per question.** ``solve`` resolves
the question into one of three values: the fundamental mode (its :math:`k`, its gauged flux
coefficients ``(G, N)`` and its emission on the emission space, scaled by
the same gauge), a mode nearest :math:`\tau` (its eigenvalue alone) or a
source answer (a flux, ``(G, N)``, and the emission it was formed from).
``solve`` is a traced memo keyed on the derivation alone (#592), and the
``answer`` property holds its decoded value for the instance's life, so
``angular_flux``, which is not memoised, reads it without re-reading the
store. ``evaluate`` then matches on the observable and the answer:

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Observable
     - ``Eigen``, ``Fundamental``
     - ``Eigen``, ``Nearest(tau)``
     - ``FixedSource``, ``Response``
   * - ``Eigenvalue()``
     - :math:`k`
     - the eigenvalue nearest :math:`\tau`
     - refused before any solve (``ValueError``)
   * - ``FluxIntegral(w)``
     - the pairing of :math:`w` with the gauged flux
     - refused after the solve, naming the :math:`k` found
       (``NotImplementedError``)
     - the pairing of :math:`w` with the flux (for ``Response``, the
       adjoint scalar flux)
   * - ``PointValue(x, g)``
     - the reading at :math:`x` of group :math:`g` of the gauged emission
       (:ref:`characteristic-reading`)
     - refused after the solve, naming the :math:`k` found
       (``NotImplementedError``)
     - the reading of the emission (for ``Response``, the adjoint scalar
       flux at :math:`x`)

What only the solve decides is refused at the first evaluation: a body
supercritical for its fixed source or its detector
(:class:`~orpheus.derivations.common.dense_pencil.NoLeastSolution`, the
fourth rung's source pencil), a complex nearest eigenvalue, and a pencil
with no fundamental mode.

**A reading is an entry; the solve it reads is a child entry.** Two
methods are traced memos (:mod:`orpheus.numerics.traced_memo`, #405 P3):
``evaluate``, keyed on the derivation and the observable, and ``solve``,
keyed on the derivation alone. The generation of a reading runs in a fresh
process, rebuilds the derivation from its content and reads ``answer``,
which calls ``solve``; that call is a memo lookup inside the generation,
so the reading's manifest records the solve as one
:class:`~orpheus.numerics.traced_memo.ChildPin`, and every observable of
one derivation reads the same solve entry. This is the unit of caching
the user ruled for P3 (ruling 2 of ``.claude/plans/reference_cache.md``:
the solve is a child entry shared by every observable of that solve), and
the one the trajectory-resolvent family realised through its solver
functions (:ref:`verification-reference-cache-clients`). Until #592 the
solve was a plain ``cached_property`` inside the reading's process, so
each observable re-solved the body in its own generation. `[M]`
2026-10-10, the elegance review's probe on the two-region test sphere at
the tiny resolution (``_hetero_sphere``, ``_TINY``): three observables of
one cold fundamental (its ``Eigenvalue``, one ``FluxIntegral`` of an
indicator and one ``PointValue``) wrote 3 ``evaluate`` entries and 1
``solve`` entry, each ``evaluate`` manifest holding exactly one
``ChildPin``; the first reading took 6.2 s (it includes the solve), each
later one 2.7 s. An equal derivation built afresh reads a stored value
without solving: `[M]` 2026-10-07, qa's probe, 3.64 s cold and 0.013 s
warm on one body. A refusal decided before the solve (an ``Eigenvalue``
of a source question) writes no solve entry; a refusal decided by the
solve (a flux of a ``Nearest`` mode) leaves the solve entry behind,
since the solve itself succeeded. The gates are
``test_m4_3c_the_characteristic_reading_manifest_holds_its_construction_and_pins_one_solve``,
``test_m4_11_two_observables_of_one_cold_derivation_share_one_solve``,
``test_m4_12_a_served_solve_is_read_only_and_its_readings_are_the_in_process_ones``
and ``test_m4_13_a_refusal_before_the_solve_writes_no_solve_entry``
(``tests/gates/numerics/test_traced_memo_clients.py``).

**Gotcha: the two routes hand out arrays of different mutability.** A
solve served from the store decodes fresh, read-only arrays; under
``bypass()`` (the in-process route the value gates use) ``answer`` holds
the generator's own writeable arrays. `[M]` 2026-10-10, qa and the
elegance review: no consumer writes into an answer today (7 of 7
``.answer`` reads in tracked ``.py`` only read), but a consumer that did
would pass every bypass row and raise against the store. Never write into
``answer``'s arrays.

**Gotcha: the solve key needs a set's order fixed.** The gauge of an
eigen question with an (n,2n) material is a ``CellCoefficient`` whose
``cells`` frozenset holds a fission and an (n,2n) cell, and a frozenset's
iteration order follows the process's hash seed. While the exact key kept
iteration order, the test process and a generation (which runs under
``PYTHONHASHSEED=0``) keyed one solve differently, so the parent's
``answer`` missed the solve its readings had written and wrote a second.
The exact key now orders a set by its elements' encodings
(:ref:`verification-reference-cache-key`).

How a function enters: the role arrows
--------------------------------------

The system's unknown is an angle-integrated emission rate per unit
volume (:eq:`characteristic-pencil`), and the interface vocabulary hands
the door mesh-free functions
(:mod:`orpheus.numerics.mesh_free_function`): a ``RegionwiseConstant``
table ``(regions, groups)`` on the angle-integrated space, or a
``Symbolic`` expression per group in :math:`(r, \mu, \varphi)`, which is
a density over :math:`\mathrm d\Omega`. A function's ROLE decides how it
enters phase space (the user's ruling of 2026-10-02,
:ref:`spaces-collapse-pair-two-lifts`), and the rate the transport reads
is the retraction of that lift. The arrows are those of Branch 1's
angular measure, :mod:`orpheus.derivations.common.angular_measure`: the
retraction :math:`Rq = \int_{S^2} q\,\mathrm d\Omega`, its section
:math:`EQ = Q/m` and its adjoint, the pullback :math:`R^\dagger\Sigma =
\Sigma`, with :math:`m = R1 = 4\pi` derived there by integrating 1 over
the sphere, never typed. They satisfy

.. math::
   :label: characteristic-door-lifts

   R \circ E = \mathrm{id}, \qquad R \circ R^\dagger = m = 4\pi ,

and the door reads each role through them
(``_coefficients(function, lift)``, with the lift ``_section``,
``_pullback`` or none):

.. list-table::
   :header-rows: 1
   :widths: 24 36 40

   * - Role and form
     - Lift into phase space
     - Rate the transport reads
   * - a source, table :math:`Q`
     - the section, :math:`EQ`
     - :math:`R E Q = Q`
   * - a detector, table :math:`\Sigma_d`
     - the pullback, :math:`R^\dagger\Sigma_d`
     - :math:`R R^\dagger \Sigma_d = 4\pi\Sigma_d`
   * - a source or a detector, ``Symbolic`` :math:`f`
     - none: :math:`f` is already a density over the directions
     - :math:`R f = \int f\,\mathrm d\Omega`, which is :math:`4\pi f` for
       an isotropic :math:`f`
   * - a weight :math:`w` of ``FluxIntegral``
     - none: it pairs with the scalar flux
     - :math:`w` as given

For a table the door computes the scale as the retraction of the lift of
the unit rate, ``retraction(lift(1))``, so the :math:`4\pi` comes from
the angular measure at every call. A symbolic source and a symbolic
detector are both retracted, since both are densities over the
directions; for a table, the source and the detector differ by
:math:`4\pi`.

.. implements:: characteristic-door-lifts
   :by: orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation

   **Implemented by** the door's ``_coefficients``, with the lifts
   ``_section`` and ``_pullback`` and the retraction ``_retraction``, each
   a call into :mod:`orpheus.derivations.common.angular_measure`
   (``section``, ``pullback``, ``retraction``).

.. implements:: characteristic-door-lifts
   :by: orpheus.derivations.common.angular_measure.retraction

.. implements:: characteristic-door-lifts
   :by: orpheus.derivations.common.angular_measure.section

.. implements:: characteristic-door-lifts
   :by: orpheus.derivations.common.angular_measure.pullback

**Reciprocity follows from the arrows.** Write the forward transport with
its boundary laws as :math:`\mathcal L\psi = E(S + F)R\psi + EQ`
(:math:`\mathcal L = \Omega\cdot\nabla + \Sigma_t`, the walls' laws in
its domain), so the angular flux of the source :math:`Q` is
:math:`\psi = \mathcal T^{-1}EQ` with :math:`\mathcal T = \mathcal L -
E(S + F)R`. A detector reads the scalar flux,
:math:`\langle\Sigma_d, \phi(Q)\rangle = \langle\Sigma_d, R\psi\rangle`.
The importance of the detector is
:math:`\psi^\dagger = \mathcal T^{-\dagger}R^\dagger\Sigma_d`, the
vocabulary's ``Response`` (:math:`E(p_0)^{-\dagger}R` in
:class:`~orpheus.numerics.question.Response`), and its scalar flux is
:math:`\phi^\dagger = R\psi^\dagger`. Then

.. math::

   \langle\Sigma_d, R\psi\rangle
   = \langle R^\dagger\Sigma_d, \mathcal T^{-1}EQ\rangle
   = \langle \mathcal T^{-\dagger}R^\dagger\Sigma_d, EQ\rangle
   = \langle E^\dagger\psi^\dagger, Q\rangle
   = \frac{\langle R\psi^\dagger, Q\rangle}{m},

since :math:`E^\dagger = R/m`. For two tables this is

.. math::
   :label: characteristic-door-reciprocity

   \bigl\langle \Sigma_d,\ \phi(Q) \bigr\rangle
   \;=\; \frac{1}{4\pi}\,\bigl\langle Q,\ R\psi^\dagger(\Sigma_d) \bigr\rangle ,

the reading of ``FixedSource(Q)`` by ``FluxIntegral(Σ_d)`` against the
reading of ``Response(Σ_d)`` by ``FluxIntegral(Q)``. The :math:`4\pi` is
:math:`R\circ R^\dagger`, the difference between the two roles of one
table.

.. implements:: characteristic-door-reciprocity
   :by: orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation

The response as the transposed forward problem
----------------------------------------------

The fourth rung already answers the adjoint, by transposing the source
pencil (:eq:`characteristic-fixed-source`,
``GalerkinSystem.response``). The door answers a ``Response`` by a
second route, the user's ruling of 2026-10-07 (Q2): the forward problem
on the adjoint cross sections, with the detector's rate as its source.

**The derivation.** The adjoint of the transport operator reverses the
direction: :math:`\mathcal L^\dagger = \mathcal P\,\mathcal L\,\mathcal P`,
with :math:`\mathcal P\psi(x, \Omega) = \psi(x, -\Omega)` the parity, for
every wall law the reference serves: a specular return of any amplitude
(vacuum is amplitude 0), the periodic wrap and a diffuse return of any
amplitude each map an outgoing direction to an incoming one by a kernel
symmetric under the reversal of both, so each is its own adjoint up to
:math:`\mathcal P` `[R]`. qa measured it on the built code: the
reciprocity :math:`\langle r, \phi(q)\rangle = \langle q,
\phi^\dagger(r)\rangle` holds in 7 of 7 configurations to between 0 and
:math:`1.1 \times 10^{-15}` (four slabs, between a vacuum and a white
wall, periodic, between a partial and a mirror wall, and between a white
wall of amplitude 0.5 and vacuum; a white sphere; a partially reflecting
cylinder; a hollow sphere; two groups with upscatter, asymmetric
transfers and several regions), and with the
transposition neutered the two sides differ by 0.95 to 1.77 relative in
all seven (`[M]` 2026-10-07, qa's ``p3`` and its control ``p3b``). The retraction does not see the direction, so
:math:`R\mathcal P = R` and :math:`\mathcal P R^\dagger = R^\dagger`. The
adjoint problem reads

.. math::

   \mathcal L^\dagger\psi^\dagger
   = R^\dagger (S + F)^{\mathsf T} E^\dagger \psi^\dagger + R^\dagger\Sigma_d
   = E\,(S + F)^{\mathsf T} R\,\psi^\dagger + E\,(R R^\dagger\Sigma_d),

using :math:`R^\dagger = mE` on an isotropic function and
:math:`E^\dagger = R/m`. Put :math:`\psi^\dagger = \mathcal P\tilde\psi`
and apply :math:`\mathcal P`: :math:`\tilde\psi` solves the FORWARD
problem :math:`\mathcal L\tilde\psi = E(S^{\mathsf T} +
F^{\mathsf T})R\tilde\psi + E(RR^\dagger\Sigma_d)`, with the same
:math:`\Sigma_t`, the transfer matrices transposed and the source rate
:math:`RR^\dagger\Sigma_d`. Its scalar flux is the adjoint scalar flux,
:math:`R\tilde\psi = R\psi^\dagger`:

.. math::
   :label: characteristic-door-response

   R\psi^\dagger(\Sigma_d) \;=\; \phi\bigl[\Sigma_t,\ S^{\mathsf T},\ F^{\mathsf T};\ RR^\dagger\Sigma_d\bigr],
   \qquad
   F^{\mathsf T} = (\chi\otimes\nu\Sigma_f)^{\mathsf T} = \nu\Sigma_f\otimes\chi ,

with :math:`\phi[\cdot\,;\,\cdot]` the forward scalar flux of the cross
sections and the source rate named. This is the self-adjointness of the
isotropic transport (the reciprocity that makes the transport block
symmetric, :ref:`characteristic-galerkin-assembly-section`) written as a
recipe: the importance is computed with the forward machinery unchanged.

.. implements:: characteristic-door-response
   :by: orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.transposed

   **Implemented by** :meth:`RegionCrossSections.transposed
   <orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.transposed>`
   (the same ``total``, ``scattering`` with ``[to, from]`` swapped,
   ``spectrum`` and ``production`` exchanged) and the door's ``Response``
   branch, which poses the system on those cross sections with the
   detector's pullback as its source.

.. implements:: characteristic-door-response
   :by: orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation

**The adjoint cross sections exchange the fission's two factors.** The
fission matrix is the outer product :math:`\chi\otimes\nu\Sigma_f`, so
its transpose is :math:`\nu\Sigma_f\otimes\chi`: in the adjoint problem
the spectrum the fission emits into is the forward production, and the
production that drives it is the forward spectrum.
:class:`~orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections`
therefore stores the two factors, ``spectrum`` (:math:`\chi`) and
``production`` (:math:`\nu\Sigma_f`), and derives ``fission`` as their
product, and ``transposed()`` exchanges them, so each field keeps its
meaning on the adjoint set. The adjoint problem's emission support
(:meth:`RegionCrossSections.emission_support
<orpheus.derivations.continuous.characteristic.cross_sections.RegionCrossSections.emission_support>`
on the transposed set) is where something is emitted OUT of a group in
the forward problem, which is where the importance has an emission.

**Why this route and not the forward system's adjoint.** The alternative,
Q2 (b) of the sketch, is ``system.response`` on the forward system. Its
unknown is the adjoint flux on the FORWARD emission supports
(:eq:`characteristic-adjoint`), so a flux integral of the importance
over a region outside them, and a reading at a point there, would
need a further transport of the adjoint emission
:math:`r + (S + F)^{\mathsf T}\phi^\dagger`, whose support is not the
forward one. The transposed problem's flux is read on every node
(:math:`W_G\phi = Kq` has a row for every basis function), and every
observable reads it as it reads a forward flux. The two routes share the
transport blocks :math:`K_g` and nothing above them, so they remain two
independent checks of each other
(``test_the_doors_response_is_the_forward_systems_adjoint_on_the_same_posing``,
below): one solves :math:`(W_s' - (S' + F')K')q' = W_s'(4\pi r)` on the
transposed supports, the other :math:`(W_s - K^{\mathsf T}(S +
F)^{\mathsf T})\phi^\dagger = K^{\mathsf T}r` on the forward ones, and
the door's reading is :math:`4\pi` times the second's paired with
:math:`W_s q`. ``system.response`` is the importance per unit emission
RATE, the door's answer the adjoint SCALAR flux.

**Why the answer is** :math:`R\psi^\dagger` **and not the per-rate
importance.** The first build read the detector table as a rate and
returned :math:`E^\dagger\psi^\dagger = R\psi^\dagger/4\pi`. That value is
self-consistent (it is the importance per unit source rate), but nothing
in the vocabulary names it: ``Response`` is :math:`\psi^\dagger`, and
:class:`~orpheus.numerics.observable.FluxIntegral` pairs a weight with
the scalar flux, which for :math:`\psi^\dagger` is :math:`R\psi^\dagger`.
Production's S\ :sub:`N` importance is that same scalar flux
(:math:`\sum_n w_n\psi^\dagger_n`). qa's review measured the
inconsistency (`[M]` 2026-10-07, ``qa/p2.py``, a two-group vacuum sphere
with symmetric scattering, where reciprocity makes the ratio a pure
convention): ``Response(f)`` over ``FixedSource(f)`` read
:math:`1/4\pi = 0.0795774715` for a symbolic :math:`f` and exactly 1 for
a table, where :math:`R\psi^\dagger` makes them 1 and :math:`4\pi`. The
first consumer to compare the reference with production would have
differed by :math:`4\pi`, the class of ERR-004. The user ruled
:math:`R\psi^\dagger` (2026-10-07), and the detector's table now enters
by the pullback, its rate :math:`4\pi\Sigma_d`.

The eigen gauge and the modes
-----------------------------

An eigenvector is defined up to a scale, and a flux integral of an
eigen answer reads that scale; a ``Ratio`` of two flux integrals does
not (the gate ``test_the_factory_returns_an_uncertified_reference_that_reads_a_ratio_as_a_quotient``
reads an eigen ratio equal to the unscaled flux's to
:math:`10^{-13}`). The scale is therefore part of the question: an
``Eigen`` carries a ``gauge``, the production its flux is scaled by
(:ref:`structured-geometry-question-values-gauge`). The specification
resolves it: a declared gauge is a
:class:`~orpheus.data.cells.CellCoefficient` (a set of emission cells,
the same kind of key as the parameter) and must resolve on the
materials; no gauge means the default
``EIGEN_GAUGE = CellCoefficient.every(FISSION_EMISSION, N2N_EMISSION)``
(``orpheus.specification.specification.EIGEN_GAUGE``), resolved, or
no gauge at all when no material carries either channel (a pure
scatterer has no production to scale by). The door fixes the scale of
the fundamental mode so that the declared production over the body is 1:

.. math::
   :label: characteristic-door-gauge

   \sum_g \int_V p_g\,\phi_g\,\mathrm dV
   \;=\; \sum_g \bigl(W\,p_g\bigr)^{\mathsf T}\phi_{h,g} \;=\; 1 ,
   \qquad
   p_g(r) = \sum_{\text{channels of the gauge in } r} e_g ,

the pairing below with the weight :math:`p`, read onto the nodes exactly.
:math:`e_g` is what one channel emits per unit flux of group :math:`g`:
:math:`\nu\Sigma_{f,g}` for fission,
:math:`2\sum_{g'}\Sigma_{2,g\to g'}` for the (n,2n) reaction and
:math:`\sum_{g'}\Sigma_{s0,g\to g'}` for scattering
(``orpheus.derivations.common.eigenvalue.production_emission``). The
default is therefore :math:`p = \nu\Sigma_f + 2\sum_{g'}\Sigma_{2,g\to
g'}`, the functional ``SNSolver.compute_production_rate`` integrates. The
flux of ``pencil.fundamental()`` is divided by its declared production.

**The declaration is shared, the physics of each channel is not.** Which
channels count is the question's, one value every reader resolves the
same way. What a channel emits is each side's own:
``production_emission`` sits on the reference side and writes the
(n,2n) multiplicity as the references' own literal 2, beside
``group_emission``'s, so that no reference moves with production's
``N2N_MULTIPLICITY`` (the reason the (n,2n) multiplicity census keeps the
two apart, ``tests/gates/transport/test_n2n_multiplicity_census.py``).
The readers of the declaration, as built (2026-10-08): the door scales
its fundamental flux to a declared production of 1 over the body; the
exact infinite medium reads its flux at a declared production density of
100 per unit volume; the trajectory resolvent accepts the fission gauge
only, since its domain has no (n,2n) emission. The S\ :sub:`N` and
homogeneous solvers do not read the declaration yet (#517): S\ :sub:`N`
scales to the default's functional, the homogeneous solver to a fission
production density of 100, so a test that judges the homogeneous solver
declares the fission gauge.

.. implements:: characteristic-door-gauge
   :by: orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation.solve

   **Implemented by** the door's ``solve`` (the fundamental's flux divided
   by its declared production; ``answer`` reads its memo entry) and
   ``orpheus.derivations.common.eigenvalue.production_emission`` (what each
   channel emits per unit flux).

.. implements:: characteristic-door-gauge
   :by: orpheus.derivations.common.eigenvalue.production_emission

**Why a total of 1 and not a density of 100.** The sketch's Q3 named the
gauge :math:`\langle\nu\Sigma_f, \phi\rangle = 100` "production's gauge,
which the exact homogeneous reference already uses", so that a closed
homogeneous body would read that reference's flux. The premise was
wrong. The 100 of
``orpheus.derivations.common.exact_homogeneous`` and of the homogeneous
solver (``ScaleGauge(production_rate.evaluate, 100.0)``) is a production
DENSITY of an infinite medium, per unit volume, since the medium has no
volume to integrate over; a finite body's gauge is a TOTAL. Gauging the
body's total by 100 equates a total with a density, and the agreement
the sketch wanted holds only by that pairing: on a closed homogeneous
sphere the body's group integral equalled the medium's density reading,
444.44, while the mean flux was :math:`444.44/V = 13.263`, and at a fixed
physical state the mean flux of such a gauge scales as :math:`1/V` (qa's
F3, `[M]` 2026-10-07, ``main/smoke.py``, a sphere of volume 33.51). The
user ruled the total production 1 (2026-10-07), the target the
S\ :sub:`N` eigenvalue entries give their ``ScaleGauge``
(``orpheus/sn/solver.py``, ``ScaleGauge(..., 1.0)``).

**Why the production is declared and not fixed in the reference.** The
5a build fixed the gauge at the fission production
:math:`\langle\nu\Sigma_f, \phi\rangle`, and its docstring called that
"the gauge the S\ :sub:`N` solver fixes". It was not:
``SNSolver.compute_production_rate`` adds the (n,2n) emission,
:math:`2\int\sum\Sigma_{2}\phi\,\mathrm dV`, so on a body with (n,2n)
the two fluxes differ by the ratio of the two productions, a convention
neither side names (the archivist's finding, 2026-10-08). A gauge is a
choice of functional, not physics, so the user ruled it part of the
question (2026-10-08): declared once, resolved by the specification, read
by every reference and, after #517, by every solver. `[M]` 2026-10-08,
the archivist's probe (``archivist_gauge2.py``, the session scratchpad),
a mirror sphere of radius 1 of ``_UP2N`` (the gates' mixture with (n,2n)
emission) at the working resolution: the default gauge's production per
unit flux is (0.15, 0.2) against fission's (0.05, 0.2); under the
default the declared production reads 1 and the fission production
0.647; the group-0 integral is 3.5294 under the default and 5.4545 with
the fission gauge declared, a ratio of 1.5455, the inverse of 0.647.

With the total gauge a closed homogeneous body of volume :math:`V` reads
the medium's flux divided by :math:`100V`, so 100 times its group
integral equals the medium's reading. `[M]` 2026-10-08, the archivist's
probe (``archivist_rung5a_numbers.py``, the session scratchpad;
``PYTHONPATH`` set to the tree, run with ``.venv/bin/python -O`` from
outside the repository), a mirror sphere of radius 1 of the gates' PU2
mixture at the gates' working resolution: the group integrals are
1.5130723862 and 2.2408273111, 100 times them are the medium's 151.307 and
224.083 to :math:`2.2 \times 10^{-16}` and :math:`8.9 \times 10^{-16}`,
and the total production reads 1.

**A higher mode has no production gauge.** Let :math:`\phi_n` be a mode
of :math:`(\mathcal L - S)\phi = F\phi/k` with :math:`k_n \ne k_0`, and
:math:`\phi^\dagger_0` the fundamental adjoint. Pairing the forward
equation with :math:`\phi^\dagger_0` and the adjoint equation with
:math:`\phi_n` and subtracting gives the biorthogonality
:math:`(1/k_n - 1/k_0)\,\langle\phi^\dagger_0, F\phi_n\rangle = 0`, so
:math:`\langle\phi^\dagger_0, F\phi_n\rangle = 0`. The fission is rank one
at each point, :math:`F = \chi\otimes\nu\Sigma_f`, so

.. math::

   \langle\phi^\dagger_0, F\phi_n\rangle
   = \int_V \bigl(\chi\cdot\phi^\dagger_0(x)\bigr)\,\bigl(\nu\Sigma_f\cdot\phi_n(x)\bigr)\,\mathrm dV = 0 .

On a closed homogeneous body the fundamental adjoint is flat in space,
the infinite medium's importance :math:`\phi^\dagger_\infty`, and
:math:`\chi\cdot\phi^\dagger_\infty > 0`, so the factor leaves the
integral and the net production :math:`\int\nu\Sigma_f\cdot\phi_n\,\mathrm
dV` of every higher mode is exactly zero. On a heterogeneous body it is
not forced to zero, but nothing keeps it away from zero either: the
production of a higher mode can take either sign and vanish. A gauge
that divides by it divides by rounding. The default gauge adds the
(n,2n) production, which the biorthogonality does not constrain, so on a
body with (n,2n) the gauge of a higher mode is not forced to zero but is
still not kept away from it. `[M]` 2026-10-08, the same probe,
the gates' mirror slab (width 2 cm, PU2, the working resolution), the
mode nearest the closed form :math:`k(\pi/a) = 0.2862762486` (read
0.2862762300): net production :math:`-1.28 \times 10^{-16}` against an
absolute production :math:`\sum_g\int|\nu\Sigma_{f,g}\phi_g|\,\mathrm dV
= 0.1124`, a ratio of :math:`-1.1 \times 10^{-15}`, and the left-half
group-0 integral divided by it reads :math:`-8.5 \times 10^{14}`, the
:math:`-8.4 \times 10^{16}` the test-architect measured under the
withdrawn gauge 100 (``ta/probe_n4.py``). The fundamental's ratio of
net to absolute production is 1. qa and the elegance review found the
same independently: modes 2 to 6 of a closed two-group sphere served
maximum coefficients :math:`3.8 \times 10^{16}` to :math:`3.1 \times
10^{17}` (qa's F1), and a symmetric two-region slab's odd mode served a
flux integral of :math:`-1.4 \times 10^{17}` (the elegance review's V1).

The user ruled (2026-10-07) that a ``Nearest`` answer reads its
eigenvalue only and that a flux integral of it is refused, naming the
:math:`k` found, as a declared scope boundary: a higher mode's flux needs
a scale that cannot vanish on it (a norm of its emission, for instance),
which belongs to the eigen answer's contract (#529). The refusal is keyed
on the question's mode, not on the value found: ``Nearest`` at
:math:`\tau` near the fundamental is refused too
(``test_a_nearest_answer_reads_its_eigenvalue_and_refuses_a_flux_integral``).

**The chart of** :math:`\tau`. ``Nearest(tau)`` picks the eigenvalue of
the :math:`k` pencil nearest :math:`\tau` in :math:`k`, the chart in which
every reference reads ``Eigenvalue()``; the chart of a parameter is owed
to #529. Between the harmonic and the arithmetic mean of two eigenvalues
the :math:`k` and :math:`1/k` charts pick different modes, and the gate
``test_tau_is_read_in_the_k_chart`` sits there.

**The null-fission cluster.** The production matrix :math:`FK` of the
pencil has rank below its size whenever a group receives no fission
emission (:math:`\chi_g = 0`) or a region produces none, and the pencil
then has a cluster of eigenvalues at :math:`k = 0`, each a vector
:math:`FK` annihilates. Rounding scatters them around zero, some complex.
They are not modes of the fission problem, and the first build's
``Nearest(0)`` served one: `[M]` 2026-10-07, qa's F1 (``qa/p6``,
``qa/p7``), the two-group slab, ``Nearest(0)``, ``Nearest(-1)`` and
``Nearest(1e-12)`` read :math:`k = 7.7 \times 10^{-23}` and
:math:`6.8 \times 10^{-18}`, and 38 of the 80 cluster eigenvalues were
complex, so the same :math:`\tau` was refused when its nearest happened
to be complex. A mode is now kept when :math:`\lVert FKv\rVert` exceeds
the rank tolerance of :math:`FK`, its larger dimension times the machine
epsilon times its 2-norm (``_nearest``). `[M]` 2026-10-08, the probe, the
mirror slab of PU2 (:math:`\chi \otimes \nu\Sigma_f` of rank one per
node): of the 40 eigenvalues 20 are filtered (10 of them complex, the
largest :math:`|k| = 1.4 \times 10^{-16}`, the tolerance
:math:`2.2 \times 10^{-15}`), the smallest kept is :math:`k = 4.27
\times 10^{-3}`, and ``Nearest(0)`` reads it, where the unfiltered
nearest is :math:`2.8 \times 10^{-18}`.

The projection onto the panel basis
-----------------------------------

A source, a detector or a weight given as a ``Symbolic`` expression is
projected onto the panel basis in the chart's volume measure,

.. math::
   :label: characteristic-projection

   c = W^{-1}\,\bigl(\langle u_i, f\rangle\bigr)_i,
   \qquad
   \langle u_i, f\rangle = \int u_i\,f\,\mathrm dV ,

panel by panel, since :math:`W` is block diagonal by panel
(:meth:`PanelBasis.project
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.project>`).
The function is evaluated at ``Symbolic.r`` after
``Symbolic.without(mu, phi)`` (the door has already refused one that
depends on the direction). The load is integrated by Gauss–Legendre with
``source_points`` per piece on the panel ends and on the function's own
``steps`` (the points where a ``Piecewise`` expression jumps), in the
measure the mass matrix is built on: ``project`` and the mass matrix's
panel blocks (``panel_mass``) read one rule, ``_panel_rule``, so the
load and the mass share one measure, which is what makes a projection
return a function of the space to itself (the elegance review's C7).

**What is exact.** The mass matrix integrates with :math:`2(p + 1)`
points per panel. With ``source_points`` at least that, the load
integrates the product of a basis function with any function of the
panel space as exactly as the mass matrix does, so such a function is
returned to rounding: every per-region polynomial of degree :math:`p`,
except on an even panel (:ref:`characteristic-even-basis`), whose space
is the polynomials in :math:`c^2`, where the function :math:`c` is not
returned (the gate reads it off by more than :math:`10^{-4}`).
``Resolution`` refuses fewer points. A step inside a panel is integrated
exactly only once named in ``steps``; unnamed, the gate's step at
:math:`r = 1.06` moved the projection by 0.22, named it agrees to
:math:`4.0 \times 10^{-15}` (`[M]` 2026-10-07, the test-architect).

**A table is read, not projected.** A ``RegionwiseConstant`` is read onto
the nodes (:meth:`PanelBasis.on_nodes
<orpheus.derivations.continuous.characteristic.basis.PanelBasis.on_nodes>`,
``per_region[region of each node]``): the basis is nodal Lagrange, so a
constant's coefficients are its value at every node, exactly, on every
panel including an even one. Projecting it instead made exactness
depend on ``source_points``, which ``Resolution`` then admitted down to 1:
`[M]` 2026-10-07, a region-wise flux integral against the exact 444.444
was off by :math:`1.4 \times 10^{-2}` at 1 point, :math:`2.0 \times
10^{-13}` at 2 and about :math:`5 \times 10^{-15}` from 3 (qa's F4), and
a projected constant was off by 5.74 at 1 point (the elegance review's
V2); at any count it was exact to rounding and never bitwise. Read onto
the nodes, the readings at 8 and 10 points are equal bit for bit
(``test_a_regionwise_source_is_read_onto_the_nodes_whatever_the_projection_points``).

**The source's regions.** A posed source widens its group's emission
support (the fourth rung's Q1). The door reads the regions from the
coefficients: a region and a group are posed where any coefficient is
non-zero (``_regions_of``). The projection solves panel by panel, so a
zero load stays exactly zero and a zero row of a table poses no region.

.. implements:: characteristic-projection
   :by: orpheus.derivations.continuous.characteristic.basis.PanelBasis.project

   **Implemented by** ``PanelBasis.project`` on the rule ``_panel_rule``
   shared with ``panel_mass``, and ``PanelBasis.on_nodes`` for a table.

.. implements:: characteristic-projection
   :by: orpheus.derivations.continuous.characteristic.basis.PanelBasis.on_nodes

A flux integral needs no point
------------------------------

``FluxIntegral(w)`` is the pairing of the weight's coefficients with the
flux coefficients through the mass matrix:

.. math::
   :label: characteristic-door-pairing

   \bigl\langle w, \phi \bigr\rangle_h \;=\; \sum_g \bigl(W c_{w,g}\bigr)^{\mathsf T}\phi_{h,g}
   \;=\; \sum_g \int_V w_g\,\phi_{h,g}\,\mathrm dV ,

with :math:`c_w = Pw` the projection (or the nodal reading of a table)
and :math:`\phi_h` the Galerkin flux, :math:`W_G\phi_h = Kq`. The
equality holds for any :math:`w`, not only one in the basis space:
:math:`(Wc_w)^{\mathsf T}\phi_h = \langle Pw, \phi_h\rangle`, and since
:math:`\phi_h` lies in the basis space and :math:`P` is the orthogonal
projection onto it, :math:`\langle Pw, \phi_h\rangle = \langle w,
\phi_h\rangle`. No point loop and no volume rule are needed.

.. implements:: characteristic-door-pairing
   :by: orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation.evaluate

**What the reading carries.** The Galerkin flux is the projection of the
transported emission: :math:`W\phi_h = Kq` with
:math:`K_{ij} = \langle u_i, \mathcal K u_j\rangle` says
:math:`\phi_h = P\,\mathcal Kq`. So the reading reproduces every moment
of the transported flux :math:`\mathcal Kq` against the basis, and

.. math::

   \langle w, \mathcal K q\rangle - \langle w, \phi_h\rangle
   = \langle w, (I - P)\mathcal Kq\rangle
   = \bigl\langle (I - P)w,\ (I - P)\mathcal Kq \bigr\rangle ,

zero for a weight in the basis space (every ``RegionwiseConstant``) and
otherwise the product of the weight's projection error and the
transported flux's, second order. `[M]` 2026-10-07, qa's ``p8``, a step
inside a panel in the weight, against a reading on a body split at the
step: :math:`1.7 \times 10^{-7}`, :math:`3.6 \times 10^{-10}`,
:math:`3.3 \times 10^{-12}` and :math:`7 \times 10^{-14}` at
:math:`p = 2` to 5. The emission :math:`q` is itself the Galerkin
solution, so the reading carries its error too, which no weight removes.

**A correction to the sketch.** The sketch said the reading "carries the
weight's projection error". As the volume integral of :math:`w` against
:math:`\phi_h` it carries none: the test-architect's gate
``test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux``
integrates the reconstructed Galerkin flux against :math:`w` by an
independent rule (numpy Gauss–Legendre, 24 points per panel, split at
the weight's step) and reads the pairing equal to
:math:`2.2 \times 10^{-16}`, :math:`1.1 \times 10^{-16}` and
:math:`1.1 \times 10^{-16}` for a region-wise, a smooth and a stepped
symbolic weight. What the sketch named is the second form above, the
reading against the transported flux, where the weight's projection
error enters only multiplied by the flux's.


.. _characteristic-reading:

The reading at a point
======================

The second half of the fifth rung (5b) reads the scalar flux of one group
at one point, the observable
:class:`~orpheus.numerics.observable.PointValue` ``(position, group)``, and
exposes the angular flux at a point for the gates. Its module is
:mod:`~orpheus.derivations.continuous.characteristic.reading`; the lines it
integrates over are built by the same one-dimensional rules as the line
rule's, in :mod:`~orpheus.derivations.continuous.characteristic.lines`. The
design is the plan's "P1 step (b), rung 5b: API sketch"
(``.claude/plans/characteristic_reference_architecture.md``), ruled by the
user on 2026-10-08 as recommended in its four items: the point value is the
transported emission and not the basis value; the directions at the point
are the line rule's own lines read at the point's parameters, from the same
constructors, with no new ``Resolution`` field; the angular flux is exposed
for gates and mints no observable; the flux integral stays the pairing.
After the reviews the user ruled the same day that the grading law at every
impact-panel top (:ref:`characteristic-grading-law`) and the image's exact
level (:ref:`chart-and-chord-level`, #590) land in this rung,
with the cylinder's extra cost accepted until #587.

What is read: the transported emission, not the basis value
------------------------------------------------------------

The Galerkin flux is the L2 projection of the transported emission. The
block is :math:`K_{ij} = \langle u_i, \mathcal K u_j\rangle`
(:eq:`characteristic-galerkin-assembly`), with :math:`\mathcal K` the
scalar-flux operator of an isotropic emission through the body and its
walls, so the flux equation :math:`W\phi_h = Kq` reads, row by row,

.. math::

   \bigl(W\phi_h\bigr)_i \;=\; \int u_i\,\phi_h\,\mathrm dV
   \;=\; \int u_i\,(\mathcal Kq)\,\mathrm dV
   \quad\text{for every } i
   \qquad\Longleftrightarrow\qquad
   \phi_h \;=\; P\,\mathcal K q ,

:math:`P` the orthogonal projection onto the panel space in the volume
measure. The value of :math:`\phi_h` at a point therefore carries the
basis's projection error :math:`((I - P)\mathcal Kq)(x)`, which is first
order in the basis's resolution of the flux, whatever the accuracy of
:math:`q`. The reading transports the converged emission once more, to the
point:

.. math::

   \phi(x) \;=\; (\mathcal K q)(x) \;=\; \int_{S^2} \psi(x, \Omega)\,\mathrm d\Omega ,

with :math:`\psi(x, \Omega)` the angular flux the emission :math:`q`
produces at :math:`x` along the line through :math:`x` in direction
:math:`\Omega`, its walls' returns included. This is the construction
of the iterated Galerkin solution :cite:`Atkinson1997` (§3.4, Eq.
(3.4.80), after Sloan): the Galerkin solution of a second-kind integral
equation is passed once more through the integral operator, and the
projection of the result is the Galerkin solution again (Eq. (3.4.81);
here :math:`P(\mathcal Kq) = \phi_h`). The exact flux is
:math:`\mathcal K` of the exact emission, so the reading's error is
:math:`\mathcal K(q_{\rm exact} - q)`, the emission's error transported,
with no projection error added. The alternative the
sketch named, :math:`\phi_h(x)` read from the basis, costs nothing and
carries the projection error, so the gates' :math:`10^{-12}` bar against
the mpmath routes could not be met without a basis that resolves the flux
to that bar.

The reading and the block are one transport under two test measures. The
block integrates :math:`u_i\,\psi` over the lines of space; the reading
integrates :math:`\psi` over the directions at one point. Integrating the
reading of :math:`u_j` against :math:`u_i` over the body recovers the
block,

.. math::

   \int u_i(x)\,\bigl(\mathcal K u_j\bigr)(x)\,\mathrm dV \;=\; K_{ij},

which is the verification spec's C6. It holds for the transport the two
share, so it is a check of the reading's measure and of the walls' fold,
not of the transport (C6 below).

The lines through a point are the line domain's lines
-----------------------------------------------------

A line through the point :math:`x` of orbit coordinate :math:`c` has an
impact parameter :math:`b \le c`. Under the chart's group it is congruent
to the line domain's canonical line of the same coordinates
(:meth:`LineDomain.lines <orpheus.geometry.chart.LineDomain.lines>`,
:ref:`chart-and-chord-line-domain`): the group moves any line onto the
canonical one with the same impact parameter (and, on the cylinder, the
same polar angle), and it moves :math:`x` to a point of the canonical line
at the same orbit coordinate. By the crossing law
(:eq:`geometry-line-crossing-law`), the canonical line reaches the orbit
coordinate :math:`c` at the two orbit-space distances

.. math::

   s \;=\; \pm\sqrt{c^2 - b^2}

from its closest approach, one inward of it and one outward. So the
directions at :math:`x` are integrated over the line domain's own
coordinates, and each line is read at the two points where it reaches the
level :math:`c`: both are images of :math:`x` under the group, and the
two directions they stand for are the two branches of the direction map
(the sphere's :math:`\pm\mu`, the cylinder's :math:`\alpha` and
:math:`\pi - \alpha`). The parameters are
``RadialImage.parameters_at(c, side)``, the chord's own expression for its
crossing of the radius :math:`c` (:ref:`chart-and-chord-level`), so a point
on a wall or an interface is read at the chord's crossing to the last bit.
On the slab a line meets each level once, ``AxialImage.parameters_at``.

The coordinates come from the line rule's constructors, on the panel
partition with :math:`c` inserted as one more end (each new panel keeping
the cross section of the panel it splits), and the impact panels are kept
below :math:`c` (a line with :math:`b > c` does not pass through the
point). The gradings the point needs are then the gradings the line rule
already makes:

- a **tangency** :math:`b = r_k` of the line to a radius below the point
  is an impact-panel end, under the visibility substitution;
- the **grazing direction** at the point is the top impact panel's end
  :math:`b = c`, where the substitution's variable
  :math:`y = \sqrt{c^2 - b^2}` is :math:`c|\mu|` exactly (the sphere), so
  the grazing direction is graded as the tangency to the top panel is;
- a point **near a wall or an interface** is graded by the law at every
  impact-panel top (:ref:`characteristic-grading-law`): the next radius
  out, the turning slot's exponential layer and the closure's pole are its
  features near :math:`y = 0`, and the point's own level is a panel top;
- on the slab the point splits its panel, so the thinnest optical width
  that grades the cosine toward grazing includes the point's distance to
  each face.

The alternative the P1 sketch named, a separate rule over the direction
box of :meth:`Chart.directions_at
<orpheus.geometry.chart.Chart.directions_at>` split at its tangencies, is
refuted below (:ref:`characteristic-refuted`): it re-derives each grading
in a second variable.

The point's measure in line coordinates
---------------------------------------

The weight of a line is the point's measure :math:`\mathrm d\Omega/4\pi`,
written in the line coordinates. The directions at :math:`x` are the box
of :eq:`geometry-directions-at` with its constant density :math:`\rho`, so
:math:`\mathrm d\Omega/4\pi = (\rho/4\pi)\,\mathrm dq` in the box
coordinates :math:`q`, and the weight is that density times the Jacobian
from the line coordinates to the box's, times the line coordinates'
quadrature weights. Per chart, with :math:`x = c\,\hat e_x`:

- **The sphere off its centre.** The box is the cosine
  :math:`\mu = \Omega_x \in [-1, 1]`, :math:`\rho = 2\pi`, so
  :math:`\mathrm d\Omega/4\pi = \tfrac12\mathrm d\mu`. A direction of
  cosine :math:`\mu` has :math:`b = c\sqrt{1 - \mu^2}`, so
  :math:`\mu = \pm\sqrt{c^2 - b^2}/c` and
  :math:`2\mu\,\mathrm d\mu = -2b\,\mathrm db/c^2`:

  .. math::

     \tfrac12\,|\mathrm d\mu| \;=\; \frac{b\,\mathrm db}{2c\sqrt{c^2 - b^2}}
     \qquad\text{on each branch } \pm\mu .

- **The cylinder off its axis.** The box is :math:`(\alpha, w)` with
  :math:`\alpha \in [0, \pi]` the in-plane angle and :math:`w = |\Omega_z|`,
  :math:`\rho = 4`, so :math:`\mathrm d\Omega/4\pi = \mathrm d\alpha\,\mathrm dw/\pi`.
  The line domain's coordinates are :math:`b` and the polar angle
  :math:`\theta`. The in-plane image passes at :math:`b = c\sin\alpha`,
  so :math:`\mathrm db = c\cos\alpha\,\mathrm d\alpha` and
  :math:`\mathrm d\alpha = \mathrm db/\sqrt{c^2 - b^2}` on each branch
  :math:`\alpha` and :math:`\pi - \alpha`; and :math:`w = \cos\theta`, so
  :math:`\mathrm dw = \sin\theta\,\mathrm d\theta`:

  .. math::

     \frac{\mathrm d\alpha\,\mathrm dw}{\pi}
     \;=\; \frac{\mathrm db\,\sin\theta\,\mathrm d\theta}{\pi\sqrt{c^2 - b^2}}
     \qquad\text{on each branch.}

- **The slab.** The box is the cosine, :math:`\rho = 2\pi`, and it is the
  line domain's own coordinate: :math:`\tfrac12\mathrm d\mu`, each line
  read once.
- **The strata.** At the sphere's centre every direction is a diameter:
  the box is empty (``WHOLE``, :math:`\rho = 4\pi`), the rule is the one
  line :math:`b = 0` of weight 1, read at its closest approach. On the
  cylinder's axis every direction has :math:`b = 0`: the box is
  :math:`w` (``AXIAL_COSINE``, :math:`\rho = 4\pi`), the rule is the lines
  :math:`b = 0` at the polar rule's angles, of weight
  :math:`\sin\theta\,\mathrm d\theta`, each read at its closest approach.
  The measure is a Dirac in :math:`b` there, so the stratum is its own arm
  of the ``match`` on the pair (line shape, direction shape), as the
  line rule's is on the line shape; the predicate is
  :attr:`DirectionDomain.on_stratum
  <orpheus.geometry.chart.DirectionDomain.on_stratum>`.

The square root :math:`\sqrt{c^2 - b^2}` in the two Jacobians is the point's
half-chord, and it is formed from each node's level, not from :math:`b`
(below). Summed over the lines and their branches the weights are 1: the
gate ``test_the_point_rule_integrates_the_direction_moments`` holds the
total, :math:`\langle(\Omega\cdot\hat r)^2\rangle = 1/3` and
:math:`\langle(\Omega\cdot\hat r)^4\rangle = 1/5` (and
:math:`\langle\Omega_z^2\rangle = 1/3` on the cylinder) at 14 points, the
centre, the axis, walls, interfaces and a cavity's wall among them, with
:math:`\hat r` read from the point the rule places on each line, to
:math:`10^{-13}` (`[M]` 2026-10-08, worst :math:`2.2 \times 10^{-15}`).

**The reading is a quadrature over the lines through the point.** The
lines carry every flux as :math:`4\pi` times a flux per steradian (their
weights hold the :math:`1/4\pi`,
:data:`~orpheus.derivations.continuous.characteristic.closure.FULL_SOLID_ANGLE`),
so with :math:`\Psi_L = 4\pi\psi` along line :math:`L`, :math:`t_\sigma(x)`
its parameter at the point on side :math:`\sigma`, and
:math:`\omega_L^x` the weight above,

.. math::
   :label: characteristic-quadrature

   \phi(x) \;=\; \int_{S^2}\psi(x, \Omega)\,\mathrm d\Omega
   \;=\; \int 4\pi\psi\;\frac{\mathrm d\Omega}{4\pi}
   \;\approx\; \sum_{L} \omega_L^x \sum_{\sigma} \Psi_L\bigl(t_\sigma(x)\bigr),
   \qquad
   \omega_L^x \;=\; \frac{\rho}{4\pi}\,J_L\,W_L ,

with :math:`W_L` the line coordinates' quadrature weight, :math:`J_L` the
Jacobian of the chart above, and :math:`\sigma` over the two sides off a
radial chart's stratum and the one point otherwise. As a functional of the
group's emission coefficients it is a row, :math:`\phi_g(x) = r_g(x)\cdot
q_g` (``PointRule.row``), the discrete measure per evaluation point of the
design.

.. implements:: characteristic-quadrature
   :by: orpheus.derivations.continuous.characteristic.reading.PointRule.of

   **Implemented by** ``PointRule.of``, which places the lines and their
   weights :math:`\omega_L^x`, and ``PointRule.row``, which reads each
   line's transport at the point and sums; ``GalerkinSystem.point_flux``
   applies each group's row to its emission.

.. implements:: characteristic-quadrature
   :by: orpheus.derivations.continuous.characteristic.reading.PointRule.row

.. implements:: characteristic-quadrature
   :by: orpheus.derivations.continuous.characteristic.system.GalerkinSystem.point_flux

One weighted set of lines, two test measures
--------------------------------------------

:class:`~orpheus.derivations.continuous.characteristic.lines.Lines` is a
weighted set of the oriented lines through a body for one group: the
basis, the walls, the group's :math:`\Sigma_t`, the line coordinates, their
weights, their exact levels (required on a radial chart, absent on the
slab, refused otherwise: *a radial chart's lines carry their exact levels
and a slab's carry none*), and the chunk and the piece budget. It holds the
guard, the order by projected speed (``ordered``) and the chunking
(``chunks``, a chunk over the budget halved until it fits). Two roles hold
one, each with its own measure:

- :class:`~orpheus.derivations.continuous.characteristic.assembly.LineRule`,
  the lines of space, weight :math:`W_L\varrho/4\pi`
  (:eq:`geometry-line-domain`), whose functional is the block
  (``LineRule.transport``);
- :class:`~orpheus.derivations.continuous.characteristic.reading.PointRule`,
  the directions at one point, weight :math:`\omega_L^x`, whose functional
  is the point's row (``PointRule.row``), with the point's orbit coordinate
  beside the lines and its sides derived from the direction domain's shape.

Both are built by the one-dimensional rules of
:mod:`~orpheus.derivations.continuous.characteristic.lines`:
:func:`~orpheus.derivations.continuous.characteristic.lines.impact_rule`
graded at a projected speed by
:class:`~orpheus.derivations.continuous.characteristic.lines.ImpactPanels`
(the grading law's speed-free data per impact panel, and its ``rule`` at
one speed),
:func:`~orpheus.derivations.continuous.characteristic.lines.polar_rule`,
:func:`~orpheus.derivations.continuous.characteristic.lines.cosine_rule`
(the mirrored rule over both signs) and
:func:`~orpheus.derivations.continuous.characteristic.lines.impact_per_polar`
(the cylinder's iterated rule, each polar angle's impact rule built at
its own speed), so a grading is written once and both measures read it. Holding the roles apart keeps a block from being
computed on a point's measure, which was spellable when the point's lines
were a ``LineRule`` (the elegance review's second round: a "block" with
:math:`\mathbf 1^{\mathsf T}K\mathbf 1 = 1.49` and no meaning).

**C6, the identity of the two measures.** One constructor is checked
structurally: on a wall, where no end is inserted, the point rule's line
coordinates are the line rule's, the same set bit for bit
(``test_on_the_outer_wall_the_point_rule_has_the_line_rules_lines``, the
sphere's and the cylinder's walls and both slab faces). One transport is
checked by value: the reading of each :math:`u_j`, integrated against
:math:`u_i` over the body by a graded volume rule written in the test (8
geometric layers of ratio 1/4 toward each panel end, 12 points per piece),
equals :math:`K_{ij}` to :math:`10^{-11}` of :math:`\max|K|`
(``test_the_volume_integral_of_the_reading_is_the_galerkin_block``,
vacuum, a partial mirror 0.6 and white walls; `[M]` 2026-10-08, the
test-architect, :math:`1.2 \times 10^{-12}` to :math:`1.3 \times 10^{-12}`;
at 4 layers and 8 points the volume rule's own floor was
:math:`2 \times 10^{-9}` to :math:`3 \times 10^{-9}`, so the band is the
test's rule, not the reading). The regional transfers read through points
alone are symmetric to :math:`10^{-11}`
(``test_the_region_transfers_read_through_points_are_reciprocal``, the P1
spec's C8).

The diffuse walls through their currents
----------------------------------------

A white wall returns, per emission function, the current
:math:`j = \alpha(I - T\alpha)^{-1}U^{\mathsf T}q`
(:ref:`characteristic-wall-coupling`). That vector is
:attr:`WallCoupling.currents
<orpheus.derivations.continuous.characteristic.closure.WallCoupling.currents>`,
``(W, M)``, the balance-row solve. Every functional of the transport is
first evaluated on the stacked sources
(:class:`~orpheus.derivations.continuous.characteristic.lines.StackedSources`:
the :math:`M` emission functions, then a unit current entering each
diffuse wall, injected as :math:`1/D_w` on the traversals entering it), and
then folded onto the emission by
:meth:`WallCoupling.on_emission
<orpheus.derivations.continuous.characteristic.closure.WallCoupling.on_emission>`,

.. math::

   f(q) \;=\; f_{\rm emission} \;+\; f_{\rm walls}\,j ,

``on_emission(e, w) = e + w @ currents``, the one fold of the walls onto
the emission. The block is
``on_emission(line, response)``, :math:`K = K_{\rm line} + R\,j`, which is
:eq:`characteristic-boundary-resolvent`; the point's row is
``on_emission(row_line, r(x))``, with :math:`r_w(x)` the reading at the
point of a unit current entering wall :math:`w`. So the block and the
reading share one definition of the returned currents and one fold
(``test_the_block_folds_its_walls_through_the_currents``: the block is
``line + response @ currents`` bit for bit, and with no diffuse wall the
currents are ``(0, M)``). Behind white walls a homogeneous body reads the
first flight plus each wall's re-entering current from Hébert's escape and
transmission probabilities, attenuated to the point, to :math:`10^{-12}`,
and a closed white body reads :math:`q/\Sigma` flat
(``test_the_reading_behind_white_walls_is_the_escape_closed_form``).

The angular flux at a point
---------------------------

:meth:`GalerkinSystem.angular_flux
<orpheus.derivations.continuous.characteristic.system.GalerkinSystem.angular_flux>`
returns :math:`\psi(x, \Omega)` per steradian, ``(..., G)``, at coordinates
of the point's direction box (:meth:`Chart.directions_at
<orpheus.geometry.chart.Chart.directions_at>`), and
:meth:`CharacteristicDerivation.angular_flux
<orpheus.derivations.continuous.characteristic.reference.CharacteristicDerivation.angular_flux>`
returns it on the answer's scale. Each direction is read on the line
through the point in that direction (:meth:`Line.through
<orpheus.geometry.line.Line.through>`), by one traversal rule, at the
point's parameter. The line carries its exact level: in the orbit space
the point :math:`x = c\,\hat e_x` sits at the signed distance

.. math::

   s \;=\; \frac{Px\cdot P\Omega}{|P\Omega|} \;=\; \frac{c\,\Omega_x}{|P\Omega|}

from the line's closest approach, so the half-chord at the point's own
radius is :math:`c|\Omega_x|/|P\Omega|` and the level is
:math:`(c,\ c|\Omega_x|/|P\Omega|)`, exact from the direction; the point
lies past the closest approach when :math:`\Omega_x > 0` and before it
otherwise. A direction grazing a wall at the point is then read: at
:math:`\mu = \pm 2^{-40}` on a sphere's wall, both branches equal the mpmath
backward path to :math:`4.4 \times 10^{-16}` (partial mirror 0.6) and
:math:`7.8 \times 10^{-16}` (mirror)
(``test_a_direction_grazing_the_wall_is_read_through_its_level``; `[M]`
2026-10-08, the test-architect; the mpmath route itself needed 40 extra
digits there). It is the one consumer of the direction box in the package,
and it mints no observable: an angular observable is a change to the
interface vocabulary, owed to the posing sequence (#529).

Its gates: :math:`\psi` against the mpmath backward path
:math:`\psi_{\rm disk}/4\pi` at the box's cosines (normal, oblique,
grazing at :math:`10^{-3}`, beside each tangency), at four points, to
:math:`10^{-12}`
(``test_the_angular_flux_is_the_backward_paths_integral``); the box
integral of :math:`\psi`, by a rule written in the test and split at every
panel end's tangency, equal to the reading to :math:`10^{-12}`
(``test_the_angular_flux_integrated_over_the_box_is_the_reading``:
`[M]` 2026-10-08, :math:`4.4 \times 10^{-16}` at :math:`x = 1.1`,
:math:`1.2 \times 10^{-14}` on the wall, graded to :math:`2^{-40}`); at the
centre :math:`4\pi\psi` equal to the reading; and the grazing limit on the
wall.

The door reads the point value
------------------------------

The answers carry the emission beside the flux, scaled together: the
fundamental mode's gauge scales its emission vector, and a source answer
holds the emission its flux was formed from
(:meth:`GalerkinSystem.source_emission
<orpheus.derivations.continuous.characteristic.system.GalerkinSystem.source_emission>`).
``evaluate(PointValue(x, g))`` returns
``Uncertified(system.point_flux(x, emission)[g])``; a ``Nearest`` answer
refuses it after the solve, as it refuses a flux integral, naming the
:math:`k` found; ``admit_observable`` refuses a position outside
:math:`[r_0, r_n]` and a group outside the problem's before any solve; and
a ``Ratio`` of two point values is their quotient
(``test_a_ratio_of_point_values_is_the_quotient_and_gauge_free``). The gate
``test_a_point_value_is_read_from_the_solved_emission_and_refused_on_a_mode_after_the_solve``
in the door's file holds the wiring bitwise. Through the door the
readings meet closed forms: a closed homogeneous body's fundamental reads
the infinite medium's flux in the declared gauge at the centre, inside and
on the wall; a closed body with a uniform source reads
:math:`(\mathrm{diag}\,\Sigma_t - S - F)^{-1}Q`; a pure absorber's point
values equal the mpmath route; and Garcia's Case 1 is read per point
(below).

The refusals
------------

- **A point outside the body**, beyond the wall or in a hollow body's
  cavity: both readings refuse it with one message, *a point is read inside
  the body*, through one door
  (:func:`~orpheus.derivations.continuous.characteristic.reading.refuse_outside`),
  naming the interval :math:`[r_0, r_n]`
  (``test_a_point_outside_the_body_is_refused_with_one_message``).
- **A point nearer the centre than** :math:`\sqrt{\text{tiny}} \approx
  1.5 \times 10^{-154}`, tiny the smallest normal float, other than the
  centre itself: its half-chords and weights square below the normal
  range, so it is refused (*half-chords that underflow*), and the centre is
  the reading to ask for. Points :math:`10^{-12}` to :math:`10^{-150}` from
  the centre read the centre's closed form to
  :math:`5.0 \times 10^{-14}` (``test_a_point_near_the_centre_reads_the_centre``);
  :math:`10^{-160}` to the smallest subnormal are refused
  (``test_a_point_below_the_smallest_normal_float_is_refused``). The guard
  carries ``ELEGANCE-DEBT[guard] #582``: it retires when the kernel's
  half-chords are scaled so that their squares neither underflow nor
  overflow.
- **A slab direction within** :math:`\sqrt\epsilon \approx 1.5 \times
  10^{-8}` **of grazing** in the angular flux: a slab line is parametrised
  from its foot, so a grazing direction's crossings sit at parameters of
  order the width over :math:`|\Omega_x|`, whose spacing outgrows a mean
  free path; it is refused (*too large to locate the point on its line*)
  where it was read as 0 for :math:`q/(4\pi\Sigma)` at
  :math:`10^{-100}` and off by :math:`1.6 \times 10^{-3}` at
  :math:`10^{-15}` (qa's F5;
  ``test_the_slabs_grazing_direction_below_its_resolution_is_refused``).
  The mechanism is #585's: a slab line works in absolute positions along
  itself, so its crossings are differences of large numbers. The guard
  carries ``ELEGANCE-DEBT[guard] #585`` and retires when a slab line's
  crossings are formed relative to a point inside the body. The radial
  charts have no such refusal: their grazing directions are read through
  the level (above).

Performance, as measured
------------------------

The reading builds no Volterra block, the assembly's dominant cost, but
transports every basis function of the support to the point along every
line, :math:`L \times J \times N` per chunk. `[M]` 2026-10-08, the main
agent (``PointRule.of``'s docstring), before the grading law: a point at
:math:`c = 0.37` of white cylinders at 8 line points, 12 along each line:
one region, 5440 lines,
2.37 s at the line rule's chunk 512 and budget 1024, 1.49 s at chunk 4096
and budget 32 768 (1.25 GB), 1.76 s at 8192 and 131 072 (2.8 GB); three
regions, 8512 lines, 15.0 s, 8.1 s and 8.3 s. The point rule takes 4096 and
32 768, a factor 1.6 to 1.8 over the line rule's chunking, because its
pieces hold no Volterra arrays and the per-chunk overhead dominates
sooner.

`[M]` 2026-10-09, the archivist's re-measure on the iterated rule
(``scratch/characteristic_architecture/p1_step_c/archivist_587_point_cost.py``
and ``archivist_587_point_cost2.py``, with their logs; ``PointRule.of`` and
``row`` at :math:`c = 0.37`, 8 line points, 12 along each line, chunk 4096
and budget 32 768, the minimum of 3 calls in one process, the host at a
load average of about 6 to 10 on 10 cores), beside the same probe on the
tensor rule at the slowest polar speed (`[M]` 2026-10-08, the archivist,
``p1_step_b5b/archivist_point_cost.py``, a load average of about 3.5):

.. list-table::
   :header-rows: 1
   :widths: 34 13 13 13 13 14

   * - Body
     - Lines, tensor rule
     - Lines
     - Time, tensor rule
     - Time
     - Peak memory
   * - three-region sphere :math:`(0, 0.5, 1.5, 2)`, white
     - 56
     - 56
     - 0.02 s
     - 0.02 s
     - 0.2 GB
   * - one-region cylinder of radius 1, white
     - 19 584
     - 8544
     - 4.94 s
     - 3.08 s
     - 1.2 GB
   * - the same, vacuum
     - 19 584
     - 8544
     - 5.13 s
     - 3.08 s
     - 1.2 GB
   * - three-region cylinder :math:`(0, 0.5, 1.5, 2)`, vacuum
     - 44 992
     - 15 704
     - 42.8 s
     - 23.4 s
     - 1.4 GB

The peak memory is the process's, measured on the iterated rule (the
white three-region cylinder's process, which also assembled its blocks,
peaked at 2.3 GB and is not in the table). Before the grading law the
same points held 5440 and 8512 lines (above). The cylinder's point cost
is the line count the grading law produces at each polar node's own speed
(:ref:`characteristic-iterated-cylinder-rule`); the white walls add a few
stacked sources and no measurable time. Contracting the emission before
the scan (#591) divides the cost per line; grading the polar angle per
impact panel (#587's other half, not built) was measured not to cut the
``ABA`` cylinder's near-grazing lines (:ref:`characteristic-what-is-not-built`).

**This is not yet measured against the target.** The verification spec's
§8 target is the hoisted trajectory-resolvent's one transport at a point,
0.073 s per call (``scratch/characteristic_architecture/p1_verification_spec.md``
§8), a figure that was not taken by a repeat protocol; §8's protocol (the
same question at matched accuracy, alternating subprocesses, the minimum
of at least 15) has not been run. A cylinder point at 3.1 to 23 s is 42
to 320 times that figure, so the target is not met on any reading of it. The
lever named is to contract the emission before the scan, so that one
function is transported along each line instead of :math:`N` (#591).
The grading law's cost on a reading at an interface was measured on the
tensor rule only: `[M]` 2026-10-08, qa's second round, the three-region
vacuum cylinder read at its interface took 136 s and 560 s at 8 and 16
line points, against 59 s and 216 s before the law (the test with its
route; 2.3 and 2.6 times); it has not been re-timed on the iterated rule
(:ref:`characteristic-grading-law`).


.. _characteristic-sn-rows:

The S\ :sub:`N` rows that read it, and their error estimates
============================================================

Since P1 step (d) of the campaign (``d9425977``, 2026-10-10) the
characteristic reference is the reference of every S\ :sub:`N` row that
compared against the trajectory-resolvent family: 13 rows, 17 test ids.
This section is what they read, at which resolution, how the reference's
error at each row is estimated, and what the step changed. The design is
the plan's "Step (d): the context and three rulings"
(``.claude/plans/characteristic_reference_architecture.md``, the user's
rulings of 2026-10-09); the measurements are the test-architect's
(``scratch/characteristic_architecture/p1_step_d/ta_d/``) and qa's
(``…/p1_step_d/qa_d/``).

The rows and their working points
---------------------------------

.. list-table:: The S\ :sub:`N` rows on the characteristic reference
   :header-rows: 1
   :widths: 30 34 18 18

   * - Rows (``tests/gates/sn/…``)
     - Problem
     - Route to the reference
     - Working point
   * - ``verification/analytical/test_phase_c_crosscheck.py``: the
       A|B|A sphere's k and shape
     - fuel A | moderator B | fuel A at 0.5, 1.5, 2.0 cm, two groups,
       reflective, scattering order 0
     - ``_aba_reference.aba_reference``
     - :math:`p = 5`
   * - the same file: the A|B|A cylinder's k and shape (strict xfails),
       its RECORD; ``test_l1_standoff_slab_cylinder.py`` and
       ``sweep/curvilinear/test_unified_matvec_cylinder.py`` (strict
       xfails and RECORDs)
     - the same body as a cylinder
     - ``aba_reference``, ``verify_cylinder_k``
     - :math:`p = 3`, the door default
   * - ``test_phase_c_crosscheck.py``: the three homogeneous edge rows
     - mixture A, :math:`R = 2` cm, reflective; the sphere at two groups,
       the cylinder at one
     - the door, ``edge_specification``
     - :math:`p = 3`
   * - ``verification/analytical/test_partial_reflector_resolvent.py``
       (ERR-094)
     - mixture A at P0, two groups: a slab of 4 cm between albedos 0.3 and
       0.7, a sphere of radius 4 cm under albedo 0.7, both specular
     - the door, ``partial_reflector_specification``
     - :math:`p = 5`

A working point is a rung of one joint ladder,
``tests/gates/derivations/_characteristic_ladders.py``'s ``rung(p)``,
which moves every axis that enters an eigenvalue together:

.. list-table:: The joint ladder, ``rung(p)`` (grading ratio 0.4 on every rung)
   :header-rows: 1

   * - degree :math:`p`
     - grading layers
     - transport rule (line, traversal, inner points)
     - source points
   * - 2
     - 1
     - (4, 8, 8)
     - 6
   * - 3 (the door default, the door gates' resolution)
     - 2
     - (8, 12, 12)
     - 8
   * - 4
     - 3
     - (12, 16, 16)
     - 10
   * - 5
     - 4
     - (16, 20, 20)
     - 12
   * - 6, 7, 8
     - 5, 6, 7
     - (20, 24, 24), (24, 28, 28), (28, 32, 32)
     - 14, 16, 18

The source points project a symbolic source, detector or weight onto the
panel basis and do not enter :math:`k`. The sphere and slab bodies run at
:math:`p = 5` because a solve there takes seconds; the A|B|A cylinder
stays at the door default because one solve there takes about 21 minutes
(`[M]` 1259 s, 2026-10-09). Step (d) renamed no row: the rows kept the
old reference's names, because the ``#516`` census of
``test_crosscheck_harness.py`` keys on them and the user ruled the names
unchanged for that step. Step (e2), which deleted the old family,
renamed ``trajectory_resolvent`` to ``characteristic_reference`` in the
seven ids (``…_against_characteristic_reference``,
``test_phase_d_characteristic_reference_crosscheck``,
``test_cylinder_l1_sweep_vs_characteristic_reference``,
``test_unified_cylinder_l1_mr_2g_characteristic_reference`` and its
``_record`` row) and re-keyed the census with them.

How the reference's error is estimated
--------------------------------------

The module ``tests.gates.derivations._characteristic_ladders`` holds, for
each S\ :sub:`N` fixture, the reference's :math:`k` (and, on the A|B|A
bodies, its 80-ratio shape) at the working point and at rungs around it,
as `[M]` tables, with a command that re-measures each,

.. code-block:: text

   python -O -m tests.gates.derivations._characteristic_ladders <ladder>

with ``<ladder>`` one of ``aba-sphere``, ``aba-cylinder [degree ...]``,
``partial-reflector``, ``edge``, ``sphere-sn``, ``cylinder-sn``,
``records``, ``old-sphere-reference``, ``old-cylinder-reference``, run from
the repository root, serially. ``reference_error(fixture)`` and the two
shape functions turn the tables into the reference's relative error
estimate at the working point. A ladder **estimates**: it can refute a
reference and never certifies one (the user's step-5 ruling of
2026-10-03, #405 P2), so every reading stays ``Uncertified`` (#566) and the
estimate is the documented provenance of a tolerance, never a bound.

**Geometric convergence.** The hp ladder converges exponentially, so the
steps between consecutive rungs shrink geometrically: with
:math:`s_n = \lvert k_{n+1} - k_n\rvert / \lvert k_{n+1}\rvert` and
:math:`s_{n+j} = s_n\, r^j`, the working point's error is the sum of every
step above it,

.. math::

   e_n = \sum_{j \ge 0} s_n\, r^j = \frac{s_n}{1 - r},
   \qquad r = \frac{s_n}{s_{n-1}},

with :math:`r` read off the two measured steps around the working point
(``_ladder_rules.geometric_error``, which refuses :math:`r \ge 1`). The
last step alone would understate the error by the factor
:math:`1/(1 - r)`. `[M]` 2026-10-09 (``ladders_p5.log``,
``probe_joint_high.log``), the steps up from :math:`p = 3` to 7:

.. list-table:: Relative steps in :math:`k` up the joint ladder
   :header-rows: 1

   * - body
     - 3→4
     - 4→5
     - 5→6
     - 6→7
     - 7→8
     - :math:`r` at :math:`p = 5`
   * - A|B|A sphere
     - 1.48e-8
     - 9.30e-11
     - 2.75e-11
     - 1.59e-13
     - 3.17e-13
     - 0.296
   * - partial-reflector slab
     - 2.48e-7
     - 4.31e-9
     - 8.81e-10
     - 3.26e-12
     - 4.06e-12
     - 0.204
   * - partial-reflector sphere
     - 1.55e-5
     - 1.02e-6
     - 5.72e-8
     - 2.61e-9
     - 1.28e-10
     - 0.056

The last two A|B|A sphere and slab steps sit at the solve's rounding. On
these three bodies the ladder runs to :math:`p = 8`, so the model can be
checked against the ladder's own tail: each :math:`p = 5` estimate lies
between its distance to the :math:`p = 8` rung and twice that distance
(A|B|A sphere :math:`3.90\times10^{-11}` against :math:`2.80\times10^{-11}`;
slab :math:`1.11\times10^{-9}` against :math:`8.88\times10^{-10}`;
sphere :math:`6.05\times10^{-8}` against :math:`5.99\times10^{-8}`), and
``test_crosscheck_harness.py::test_each_reference_estimate_covers_the_finest_measured_rung``
asserts it from the tables. Its declared limit: :math:`p = 8` is not the
limit, so "covers" holds to the :math:`p = 8` rung's own error, at least
1.3 orders below each estimate.

**The A|B|A sphere's shape converges in pairs of degrees.** Its 80-ratio
shape steps up from :math:`p = 3` are :math:`2.31\times10^{-4}`,
:math:`9.30\times10^{-6}`, :math:`9.07\times10^{-6}`,
:math:`4.41\times10^{-7}` and :math:`4.72\times10^{-7}`: a step of one rung
does not contract from :math:`p = 4` to 5 (ratio 0.975), and its geometric
estimate would read :math:`3.6\times10^{-4}`. Steps of two rungs do
(3→5 :math:`2.21\times10^{-4}`, 5→7 :math:`9.00\times10^{-6}`, ratio
0.041), so the shape estimate is :math:`9.4\times10^{-6}`. `[M]`
(``probe_shape_axes.log``, ``probe_shape_sp.log``) the pairing is on the
degree axis alone: at :math:`p = 5` raising the grading layers moves the
shape :math:`4.9\times10^{-11}`, the line and traversal points less than
:math:`2\times10^{-14}`, the degree :math:`9.07\times10^{-6}`, and the
reading's source points 12 → 48 move it :math:`1.3\times10^{-15}`.

**The A|B|A cylinder climbs a panel ladder.** The joint ladder's rung above
the door default would cost about 2.7 hours: on the homogeneous edge
cylinder the joint step 3→4 multiplies the cost by 7.6 (29 s to 220 s).
So the cylinder's ladder moves the panel axes only, ``cylinder_panel_rung(p)``
= ``Resolution(p, p - 1, 0.4, TransportResolution(8, 12, 12), 2(p + 1))``,
with the transport rule held at the door default's, and the transport
axis enters through its own measured step. `[M]` 2026-10-09/10, one solve
each, serially (``ladder_cylinder.log``): :math:`p = 2` 1.2317299451673347
in 375 s, :math:`p = 3` 1.2317294528902538 in 1259 s, :math:`p = 4`
1.231729422811464 in 3399 s; steps :math:`4.00\times10^{-7}` and
:math:`2.44\times10^{-8}`. The cylinder's own ratio, 0.061, is not used:
its step 2→3 moves the grading layers from 1 to 2 and is dominated by
them, so it says little about the steps above (qa's finding F3). The
estimate divides the step 3→4 by :math:`1 - r` with :math:`r` = 0.2957,
the largest contraction ratio measured on any :math:`p = 5` ladder (the
A|B|A sphere's :math:`k`), and adds the transport rule's step one rung
below the default, :math:`(6, 9, 9) \to (8, 12, 12)`,
:math:`6.9\times10^{-11}` (an over-estimate of the default's own
transport error): :math:`3.5\times10^{-8}` in :math:`k`. The shape is
estimated the same way from its panel steps (:math:`8.55\times10^{-4}`,
:math:`1.79\times10^{-4}`): :math:`2.5\times10^{-4}`, with the
:math:`k` transport step standing in for a shape step nobody measured
(`[R]`; it is six orders below the panel term).

**The homogeneous edge bodies have a closed form.** Their :math:`k` is the
medium's :math:`k_\infty` (``kinf_homogeneous`` on the library's arrays),
so the reference's error is its distance from it; the estimate is the
larger of that distance and the step up, :math:`7.2\times10^{-15}` on the
two-group sphere and :math:`3.1\times10^{-15}` on the one-group cylinder.
The S\ :sub:`N` side's residual on each edge row is its frozen snapshot's
own distance from :math:`k_\infty` (:math:`4.46\times10^{-12}`,
:math:`1.48\times10^{-16}`, :math:`2.96\times10^{-16}`).

The tolerance rules (``tolerance_for``, ``summed_tolerance``,
``geometric_error``, ``richardson_error``) moved at step (d) into
``tests/gates/derivations/_ladder_rules.py``, which reads no reference, so
that both families' ladders could use them: the corroboration rows'
independence leg refused an old-side module whose imports reached the new
reference, and the old family's ladders imported the rules. Step (e2)
deleted the old ladders and the corroboration rows; the rules stay, read
by ``_characteristic_ladders.py``.

Statement (ii): the new estimate is at most the old one
-------------------------------------------------------

The verification specification's statement (ii)
(``scratch/characteristic_architecture/p1_verification_spec.md`` §7) is
that, at every row's fixture, the characteristic reference's error
estimate is at most the trajectory resolvent's, so that no tolerance
recomputed on the new reference loosens. It holds on every row, and every
tolerance tightened or held:

.. list-table:: Statement (ii), and the tolerances before and after step (d)
   :header-rows: 1
   :widths: 28 16 22 16 18

   * - Row
     - new estimate
     - old estimate
     - S\ :sub:`N` residual
     - tolerance, old → new
   * - A|B|A sphere, :math:`k`
     - 3.9e-11
     - 3.2e-4
     - 1.5e-5
     - 4e-3 → 4e-5
   * - A|B|A sphere, shape
     - 9.4e-6
     - 1.4e-3
     - 7.3e-3
     - 2e-2 → 2e-2
   * - A|B|A cylinder, :math:`k` (strict xfail)
     - 3.5e-8
     - 1.9e-3
     - 3.0e-5
     - 8e-5 → 6e-5
   * - A|B|A cylinder, shape (strict xfail)
     - 2.5e-4
     - none (#516)
     - 5.0e-3
     - 2e-2 → 2e-2
   * - A|B|A cylinder at folded 4×8 (standoff, unified; strict xfails)
     - 3.5e-8
     - none (#516)
     - 5.9e-4
     - 2e-3 → 2e-3
   * - partial-reflector slab
     - 1.1e-9
     - 4.6e-4
     - 3.0e-4
     - 1e-3 → 4e-4
   * - partial-reflector sphere
     - 6.1e-8
     - 0.9e-4
     - 6.0e-4
     - 1e-3 → 7e-4
   * - edge, sphere (two groups)
     - 7.2e-15
     - 1e-10, declared
     - 4.5e-12
     - 1e-9 → 9e-12
   * - edge, cylinder (one group, 2 ids)
     - 3.1e-15
     - 1e-10, declared
     - 1.5e-16, 3.0e-16
     - 1e-9 → 4e-14

Every tolerance is ``tolerance_for(e, b)``, the smallest
one-significant-figure :math:`T` with :math:`T \ge 10\,b` and
:math:`T \ge 2(e + b)` for the S\ :sub:`N` residual :math:`e` and the
reference's estimate :math:`b` (with no estimate, :math:`b` is assumed at
the floor and :math:`T \ge 2.5\,e`: the old cylinder rows), except the
partial-reflector rows, which keep their own sum rule,
``summed_tolerance(e, b)`` = :math:`e + b` rounded up (the user's ruling 2
of 2026-10-09), each error there being a measured distance to a limit.
On every row the reference's estimate is now at least three orders below
the S\ :sub:`N` residual, so the S\ :sub:`N` side governs every tolerance;
the A|B|A sphere's :math:`k` row lost the old family's :math:`10\,b` floor
and tightened a hundredfold. `[M]` 2026-10-09 it reads
:math:`1.8\times10^{-5}`, a margin of 2.2, and the partial-reflector
sphere reads :math:`6.54\times10^{-4}` against :math:`7\times10^{-4}`, a
margin of 1.07, which the row's docstring records.

Statement (i), the step before (step (c), ``17708ac8``), is the mirror
condition, checked by code-to-code rows deleted with the old family at
step (e2) (the ``test_characteristic_reference_corroboration`` module,
then under ``tests/gates/derivations/``, L4, no level marker): the new
reference sat inside the old one's own error estimate at each row's
problem. `[M]` 2026-10-09 (the plan's "Step
(c) landed"), new against old, against the old estimate: A|B|A sphere
:math:`k` 8.3e-5 against 3.2e-4, its shape 4.0e-5 against 1.4e-3, the
A|B|A cylinder :math:`k` 5.6e-4 against 1.9e-3, the partial-reflector
slab 3.28e-5 against 4e-5 and sphere 4.1e-5 against 1e-4. Those rows also
asserted, statically and at run time, that the old side executed nothing
of this package (``tests/gates/_corroboration.py``, which P0's kernel rows
still use).

The cylinder RECORD's re-baseline
---------------------------------

The cylinder rows cannot verify (the reference has no certificate,
#566), so their strict xfails expect the verification verbs' refusal,
``ReferenceNotValid``, and a RECORD pins what both sides read today:
``_aba_reference.CYLINDER_3REG_RECORD``, held within
:math:`2\times10^{-5}` (relative for the eigenvalues, absolute for the
gap). Step (d) re-baselined its two reference keys through the
``records`` re-measure, and its three S\ :sub:`N` keys did not move:

.. list-table:: ``CYLINDER_3REG_RECORD`` across step (d) (`[M]` 2026-10-10, ``records.log``)
   :header-rows: 1

   * - key
     - before (the trajectory resolvent at (24, 16, 32))
     - after (the characteristic reference at the door default)
   * - ``k_ref``
     - 1.231036749830859
     - 1.2317294528902538
   * - ``phase_c_k_gap``
     - 5.97e-4
     - 3.4606e-5
   * - ``phase_c_k_sn``
     - 1.2317720792844793
     - unchanged (re-read within 1 ulp)
   * - ``unified_k``
     - 1.2310184196907399
     - unchanged (re-read within 1.5e-10)
   * - ``standoff_sweep_k_nx40``
     - 1.23101841974857
     - unchanged (re-read identical)

The new ``k_ref`` agrees with step (c)'s reading of the same resolution
(1.231729452890254) to one ulp. The gap between the folded 16×32 S\ :sub:`N`
solve and the reference fell from :math:`5.97\times10^{-4}` to
:math:`3.46\times10^{-5}`, inside the row's :math:`6\times10^{-5}`. The 7
cylinder xfails stay strict, with ``raises=ReferenceNotValid``; their
reason names #566 and keeps the ``#516`` token the harness checks.

What both sides read from one object
------------------------------------

The reference shares no project primitive with the S\ :sub:`N` sweep
above the trusted-library line: transport along lines, Galerkin over
them and a dense pencil against discrete ordinates, a spatial scheme and
a sweep. What both sides read from one object is the posed problem: the
cross sections, the geometry and the boundary law. For the partial
reflectors that includes the law's response factor, ``SpecularReturn.kernel``,
which S\ :sub:`N` reads through its leakage predicate and its
curvilinear corner, and this reference through ``law.response_kernel``
in its walls (:ref:`characteristic-walls-factors`). A defect there moves
both sides together (X4: one input, so their agreement cannot see it).
`[M]` 2026-10-10 (qa, ``p1_step_d/qa_d/half_alpha_plugin.py``):
:math:`\alpha \to \alpha/2` in that factor left the partial-reflector slab
row green. Its pin is external to the comparison,
``tests/gates/derivations/test_characteristic_walls.py``, red three times
under the same mutation, which compares each wall with a wall written by
hand from the law's physics. The partial-reflector and phase C docstrings
stopped claiming full independence at step (d) and now say this.

The re-pointed supports, and what step (e) owed
-----------------------------------------------

The 8 ``rests_on`` ids that named old-family tests now name this
reference's gates that carry the same dependency: the line integral of a
per-region source that jumps at both interfaces, on the sphere and on the
cylinder (``test_characteristic_transport.py``, in place of the old
``test_mr_oracle_first_leg_matches_the_line_integral``); the Rayleigh–Ritz
nesting of the one-group :math:`k` and the flux integral of a step weight
(in place of the old sphere's radial-convergence row and its fine-rule
reading); the fundamental mode satisfying its pencil (in place of the old
reading's fixed-point identity, R7b2.2.2, under the sphere's shape row); an interface between equal materials being invisible and the
cylinder's escape and transmission probabilities against their closed
forms (in place of the old cylinder's MR↔MG reduction and its WM-72
vacuum row); and the unfolded wall-by-wall sum of a partial mirror's cycle
(``test_characteristic_closure.py``, in place of the old slab's method of
images). `[M]` 2026-10-09/10 (the commit's record): a counting spy saw 17
of 17 ids reach this reference and 0 calls of the old family; scaling the
S\ :sub:`N` reading by :math:`1 + 2\,\mathrm{tol}` reddens 10 of the 10
ids that are not xfails.

The old family's spelling on these problems moved at step (d) to the old
side of the corroboration rows (``_trajectory_resolvent_aba.py``), and its
ladder tables to ``_trajectory_resolvent_ladders.py``, both then under
``tests/gates/derivations/``. Step (e2) deleted the family with them, and
re-posed P0's corroboration rows that imported its chord oracle
(:ref:`characteristic-successors`).


.. _characteristic-successors:

The old family's tests, re-posed on this reference
==================================================

P1 step (e) retired the trajectory-resolvent family in two commits, by
the user's ruling of 2026-10-10 ("migrate, then delete"; the plan's
"Step (e): the audit and four rulings"). Step (e1b) (``44303919``) built
the successors on this reference while the old family still existed;
step (e2) deleted the family's 12 numeric modules and 22 test files in
one commit, in the order of the dependency audit
(``.claude/plans/characteristic_reference_architecture_dependency_audit.md``).
The record of (e1b), with a row per old test, its battery and its runs,
is ``scratch/characteristic_architecture/p1_step_e/ta_e1b/README.md``.

**What moved.** Of the 22 deleted files' rows, 123 were to be kept or
re-posed. Every one has a successor `[M]` (the record's section 1): 65
are new rows, 10 existing gates widened by a parameter, 29 existing
gates that already asserted the contract, 12 existing gates that
asserted it and were also widened, 2 an existing gate plus a new row,
and 5 re-homed unchanged because they never read the old family. The new
files are ``test_characteristic_albedos.py`` (the method of images, the
ordering of :math:`k` in each wall's albedo, the vacuum mode, interface
continuity, the closing cavity), ``test_characteristic_convergence.py``
(self-convergence on every resolution axis, five bodies),
``test_characteristic_independent_references.py`` (PS-1982 and the
one-group cylinder of Sood and of Westfall–Metcalf),
``test_characteristic_nystrom_withdrawn.py`` (the three Peierls–Nyström
rows withdrawn under #506), all under ``tests/gates/derivations/``, and
``tests/gates/numerics/test_symbolic_steps.py``.

**The catalogued defects.** ERR-034, ERR-035, ERR-090 and ERR-091 lived
in the old family's code; each now has a catcher on this reference, each
measured by re-dropping the defect into this reference
(:doc:`/theory/verification/error_catalog`, each entry's "caught by").
ERR-034 and ERR-090 are no longer quiet here: the old collocation hid
ERR-034 under a flat emission, while the Galerkin assembly transports
each basis function, so the closed slabs of the :math:`k_\infty` floor
red under it; and ERR-090's class, one emission piece across an
interface, is refused by the panel basis.

**The labels that moved.** Eleven ``peierls-greens-*`` labels, now on
:ref:`theory-characteristic-origins`, lost their only verifier with the
deleted files. Each was read on the family's page and placed on a row
that asserts its equation:

.. list-table::
   :header-rows: 1
   :widths: 34 36 30

   * - Label (``peierls-greens-`` plus)
     - The equation, in short
     - Its verifier now
   * - ``annulus-impact-parameter-partition``,
       ``hollow-sph-impact-parameter-partition``
     - the line through a shell meets the cavity iff its impact
       parameter is below the inner radius
     - ``test_characteristic_closure.py::test_the_period_matches_the_hand_counted_table``
       (its hollow rows, at, below and just below the inner radius)
   * - ``cylinder-T``
     - the cylinder's rank-1 specular closure
     - ``test_characteristic_transport.py::test_the_closed_angular_flux_is_the_unfolded_backward_path``
       (case ``cylinder_solid_partial``)
   * - ``cylinder-architecture``
     - :math:`\psi = F + e^{-\Sigma L}\psi_{\rm surf}` on the cylinder,
       :math:`\phi` its angular reduction
     - ``test_characteristic_reading.py::test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route``
       (slow)
   * - ``cylinder-mr-interface-continuity``
     - :math:`\phi` continuous across an interface
     - ``test_characteristic_albedos.py::test_the_modes_scalar_flux_is_continuous_across_a_material_interface``
   * - ``cylinder-mr-kinf``
     - :math:`k_\infty = \rho(A^{-1}F)`, :math:`A = \mathrm{diag}(\Sigma_t) - \Sigma_s^{\rm T}`,
       :math:`F = \chi\otimes\nu\Sigma_f`
     - ``test_characteristic_system.py::test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio``
       (its cylinder rows, slow)
   * - ``cylinder-mr-quadrature-convergence``
     - each resolution axis contracts under refinement
     - ``test_characteristic_convergence.py::test_k_contracts_on_every_resolution_axis_and_the_working_point_is_within_the_old_floor``
       (its cylinder row, slow)
   * - ``cylinder-mr-trajectory-segments``
     - the in-plane conic :math:`r(s)^2` of a cylinder line
     - ``tests/gates/geometry/test_chord.py::test_the_multi_region_segments_are_the_hand_written_table``
   * - ``cylinder-mr-wm72-vacuum``
     - :math:`k = 1` at Westfall–Metcalf's critical radius, to
       :math:`10^{-5}`
     - ``test_characteristic_independent_references.py::test_k_is_one_at_the_one_group_cylinders_published_critical_radius``
       (slow)
   * - ``mr-regionwise-source``
     - each segment reads its own region's piece of the emission
     - ``test_characteristic_transport.py::test_the_outflow_of_a_per_region_polynomial_is_its_line_integral_attenuated_to_the_exit``
   * - ``slab-asym-method-of-images``
     - :math:`k` of the slab :math:`[0, L]` with a mirror at 0 equals
       :math:`k` of the vacuum slab :math:`[0, 2L]`
     - ``test_characteristic_albedos.py::test_a_mirror_is_the_symmetry_plane_of_the_doubled_vacuum_slab``
       and ``::test_the_mirrored_slabs_flux_is_twice_the_doubled_slabs_right_half``

Paths are under ``tests/gates/derivations/`` unless given. Five more
labels were verified by 26 SymPy rows that called five ``derive_*``
functions duplicating the kept ones; step (e2) deleted the rows and the
functions, and the four labels with no other verifier moved:
``annulus-through-rank2`` and ``hollow-sph-through-rank2`` to the
angular-flux row above (cases ``cylinder_hollow_cavity`` and
``sphere_hollow_cavity``), ``cylinder-trajectory`` to the segment-table
row, and ``slab-trajectory`` to
``tests/gates/geometry/test_chord.py::test_the_slab_chord_in_both_orientations``.
``annulus-architecture`` keeps its kept SymPy verifier.

Step (f) moved the three V_α2 rows of the Peierls Nyström primitives and
the slab's substitution row off the closure labels they named while
asserting the :math:`T_{00}` identity: they now verify
``V-alpha-2``, ``slab-V-alpha-2`` and the new ``cylinder-V-alpha-2``
(:ref:`characteristic-origins-census`). ``cylinder-T`` keeps the
angular-flux row above. ``slab-T`` holds only for a source symmetric
about the slab's mid-plane, so it has a row of its own,
``test_characteristic_transport.py::test_a_symmetric_slabs_surface_inflow_is_the_one_transit_closed_form``,
and the angular-flux row's ``slab_symmetric_half`` case, whose source is
not symmetric, verifies the general cycle only.

**What the migration measured.** Four findings, each `[M]` 2026-10-10
in the record's section 5:

- *The old family's slab method-of-images value was 5.3e-4 off.* Its
  one-group mirror-vacuum slab :math:`[0, 1]` read k = 0.0584228 at its
  row's resolution; this reference reads 0.0584537 at rung 4, and
  0.0584536527 at rungs 5 and 6, where its two posings agree to
  1.6e-14. The Peierls–Nyström slab :math:`[0, 2]`, an E\ :sub:`1`
  kernel with no lines, converges onto this reference's value (0.05845381
  at 8 panels of order 6, 0.05845368 at 12 of order 8). The old row's
  1e-7 equality of the two posings was an agreement of two readings that
  share the reference's error (``vv-principles`` anti-pattern #7); the
  successor pair compares two posings too, so it carries the same
  blindness: it catches ERR-034 through what the defect does to the
  identity between the two posings, not through a value against an
  independent reference.
- *A vacuum cavity's effect is second order in its radius on a sphere,
  first order on a cylinder.* A cavity absorbs the lines through it, a
  fraction :math:`(R_{\rm in}/R)^2` on the sphere: :math:`k`'s distance
  to the solid sphere's is 1.9e-2, 2.2e-4 and 1.8e-6 at
  :math:`R_{\rm in}/R = 10^{-1}, 10^{-2}, 10^{-3}`, and the cylinder's
  9.4e-2, 1.1e-2 and 1.2e-3. The old row's 1e-7 band at
  :math:`R_{\rm in} = 10^{-3}R` was below the physics; the successor
  asserts the order, two-sided, and :math:`k_{\rm hollow} < k_{\rm solid}`.
- *PS-1982 does not converge at an optical radius of 25* (200 power
  iterations at 30 and 40 nodes), so the old thick-sphere row's 2e-3
  agreement measured PS-1982's iteration. The successor bands the thick
  sphere against :math:`k_\infty` only, and compares with PS-1982 at
  optical radii 2.5 and 5 (:math:`k` to 1e-5, the shape to 1e-4).
- *A monotone law is a weak instrument.* The ordering of :math:`k` as a
  slab thickens stays green under a vacuum wall returning half, the
  scattering matrix transposed and :math:`\Sigma_t` larger by 1e-9, each
  of which keeps the approach monotone; only a vacuum wall read as a
  mirror reddens it. Self-convergence's contraction leg alone is blind
  to the tangency substitution dropped (an algebraic convergence still
  contracts); the rate leg (each angular step at most 1e-2 of the last:
  1.4e-3 to 2.8e-4 honest, 0.26 to 0.57 under the mutation) catches it.

**What the default tier kept, and what it did not.** The cylinder's
eigenvalue stays in the tier that runs without ``slow``: the closed
cylinder and the closed annulus read :math:`k_\infty` at rung 2 to 1e-6,
and Ua-1-0-CY reads :math:`k = 1` within 3e-5 at rung 2 (its measured
error 8.4e-6). Both red by value under a line speed scaled by 1.0001.
Eight old rows that ran in the default tier have successors only in the
slow tier, each for a measured cost: the annulus's closing cavity (four
annulus solves, about 75 s), four rows of equal-material interfaces on
the cylinder (the smallest cylinder leg that keeps the contract costs
140 s), the cylinder's :math:`P_{ss}` against Bickley, the asymmetric
slab's convergence at intermediate albedos (about 200 s) and the
cylinder's pencil residual; the sphere and slab legs of each are in the
default tier.


.. _characteristic-what-is-not-built:

What the package does not compute
=================================

The package answers, for any batch of lines, each line's period and
closure, the traversal integrals of every basis function, the vacuum
Volterra block with the caller's line weights and the angular flux at
points on the lines (:ref:`characteristic-transport`); for one group, the
transport block over the emission support, its line part and its white
walls' coupling, on a line rule graded from the group's optical scale
(:ref:`characteristic-galerkin-assembly-section`); for the multigroup
problem, the emission matrices per region, each group's exact emission
support, the k pencil and the source pencil on the emission space, the
fundamental and higher modes and their adjoints, the flux of a fixed
source and the adjoint flux of a detector, all in nodal coefficients
(:ref:`characteristic-galerkin-system`); and, through the door, the
eigenvalue, the flux integrals and the point values of an ``Eigen``, a
``FixedSource`` or a ``Response`` question posed from a specification,
with a mesh-free source, detector or weight projected onto the basis
(:ref:`characteristic-door`, :ref:`characteristic-reading`), and the
angular flux at a point for the gates. What is missing (the campaign's
issue is #405):

- **a point reading at the cost of the verification spec's target**: a
  point on a cylinder costs seconds against the target's 0.073 s, and
  the spec's §8 protocol has not been run (#591,
  :ref:`characteristic-reading`, "Performance");
- **the flux of a higher mode**: a ``Nearest`` answer reads its
  eigenvalue only, and refuses a flux integral and a point value, because
  the production gauge can vanish on it; a
  scale that cannot vanish (a norm of the mode's emission) belongs to the
  eigen answer's contract, #529 (:ref:`characteristic-door`);
- **an anisotropic source or detector**: a ``Symbolic`` function that
  depends on the direction is refused, as an anisotropic emission is;
  both need an angular basis on each line;
- **an error bound**: every reading is ``Uncertified`` (#566);
- **a cylinder polar rule graded per impact panel**: the rule over
  :math:`(b, \theta)` is iterated in :math:`b`, each polar angle carrying
  an impact rule graded at its own speed
  (:ref:`characteristic-iterated-cylinder-rule`), but the polar rule is
  the same at every impact parameter, graded toward grazing over the whole
  body, and on the ``ABA`` cylinder the polar nodes below
  :math:`\sin\theta = 0.01` carry 57 536 of its 102 320 lines per block.
  Grading :math:`\theta` per impact panel is not built, and it would not
  cut them on that body: `[M]` 2026-10-09
  (``scratch/characteristic_architecture/p1_step_c/polar_half_scales.py``)
  every line crosses the thinnest panel, the outer-wall panel
  :math:`[1.96, 2.0]` (0.02 and 0.04 optical in groups 0 and 1), which
  sets the polar grazing grading, so each impact panel's polar rule keeps
  the body's 160 or 152 nodes. A polar weight absorbing
  :math:`\sin^2\theta` is unmeasured and stays in #587
  (:ref:`characteristic-cylinder-cost`);
- **a point extremely near the sphere's centre or the cylinder's axis**:
  a point closer than :math:`1.5 \times 10^{-154}` is refused, its squares
  underflowing (#582, :ref:`characteristic-reading`);
- **a wall both specular and diffuse**, refused as a scope boundary
  because production poses none (:ref:`characteristic-walls`);
- **a body far from the origin**: a slab at :math:`x \approx 10^{6}` loses
  digits to absolute positions, and for the same reason a slab direction
  within :math:`\sqrt\epsilon` of grazing is refused by the angular flux
  at a point (#585).

Since P1 step (d) (``d9425977``) this reference is the one the
S\ :sub:`N` cross-check rows read (:ref:`characteristic-sn-rows`); the
trajectory-resolvent family (:ref:`characteristic-origins-history`) was read,
outside its own gates, only by the step (c) corroboration rows that
compared the two, until step (e2) deleted it (2026-10-10).


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
     - Dissolved on the first rung: a geometry admits a ``BC`` tag or a
       ``BoundaryTraceLaw`` per boundary point, and ``LawSum`` and
       ``LawScaled`` are content nodes, not laws, so no such wall can be
       declared. On the third rung the user ruled the refusal at the type
       (2026-10-06): ``Wall`` refuses a wall that is both, as the
       S\ :sub:`N` realizer refuses a ``LawSum``
       (:ref:`characteristic-walls`).
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
   * - The symmetry of the block as the alarm for an under-integrated
       rule over lines (the first verification spec's C11)
     - Designed green: each line's quadrature is symmetric under reversing
       the line, so each line's block is symmetric whatever the rule's
       accuracy (:math:`2.5 \times 10^{-16}` where conservation missed by
       :math:`1.6 \times 10^{-3}`). Closed-body conservation is the alarm
       for the transport along a line and the coupling, and the closed
       forms of the walls are the alarm for the rule over lines
       (:ref:`characteristic-galerkin-assembly-section`).
   * - Gauss–Legendre in the cylinder's axial cosine :math:`\mu_z`, with
       "no grading needed" (the third rung's premises)
     - Measured under a mirror, where every line conserves on its own and
       no rule over lines can show. Under a white wall the square-root end
       of :math:`|P\Omega| = \sqrt{1 - \mu_z^2}` at :math:`\mu_z = 1` costs the
       escape :math:`5.4 \times 10^{-4}` at 8 points; the polar angle is
       analytic (:ref:`chart-and-chord-line-domain`).
   * - Grading the impact rule toward :math:`b = 0` to absorb the centre's
       :math:`b^{2m+2}\log b`
     - Chases a term the physics excludes: an odd mode at a singular
       stratum is not in the flux (Schwarz). The even basis removes it at
       its cause (:ref:`characteristic-even-basis`).
   * - A fixed resolution for the rule over lines (12 halvings toward the
       slab's grazing cosine; plain Gauss in the polar angle;
       ``chord_quadrature`` in :math:`b`)
     - Every feature of the integrand over lines sits at a scale the
       optical widths set: grazing at :math:`\tau_{\min}`, the rim at a
       mean free path, the normal direction at :math:`1/\tau`. A fixed
       rule is right only in the band it was measured in, and conservation
       cannot see it leave that band (ERR-101,
       :ref:`characteristic-line-rule`).
   * - The grazing halving stopped at :math:`\tau_{\min}`
     - A panel's attenuation :math:`e^{-\tau_P s}` changes until
       :math:`\tau_P s` reaches 64; stopping at :math:`\tau_{\min}` left
       one block entry off by :math:`2.0 \times 10^{-6}` while every total
       was exact (:ref:`characteristic-line-rule`).
   * - Chunks of lines bounded by a line count
     - A line's pieces grow without bound near grazing: 224 cylinder lines
       at :math:`\mu_z = 0.999` took 8.2 GB. The chunk is bounded by its
       piece slots (:ref:`characteristic-line-rule`).
   * - :math:`I - T\alpha` formed by subtraction (the ruled sketch's
       :math:`I - T_w`)
     - Rounding amplified by the inverse of the absorption
       (:math:`1.3 \times 10^{-4}` at :math:`\Sigma_t = 10^{-12}`), and a
       lossless body returned :math:`10^{16}` instead of a refusal
       (ERR-102, :ref:`characteristic-wall-coupling`).
   * - The update :math:`U\,(I - T_w)^{-1}A\,U^{\mathsf T}` (the ruling's
       spelling)
     - Omits :math:`D = \mathrm{diag}(A_w/4)`, the injection of a unit
       isotropic current; conservation misses by 6.6 on a white sphere
       (:ref:`characteristic-wall-coupling`).
   * - :math:`R` divided by its own injected tally
     - Cancels the wall's area out of the block; the route gate on the one
       density reddened. :math:`R` is per nominal unit current
       (:ref:`characteristic-wall-coupling`).
   * - The flux at a wall located by the slot's start plus its length
     - The sum rounds an ulp apart from the closing crossing's parameter,
       which is how the kernel spells the wall, and refused 64 of 119
       reads at an exit wall (qa, 2026-10-06). The exit test reads the
       crossing parameter (:ref:`characteristic-angular-flux`).

   * - The flux as the pencil's unknown, :math:`(W_G - KS, KF)` (the P1
       sketch's item 6, ruled 2026-10-06)
     - Right for :math:`k` and the flux, wrong for the adjoint: its
       transpose poses :math:`\psi = E^\ast\mathcal K\psi`, whose
       eigenvector is the adjoint collision density
       :math:`E^\ast\phi^\dagger`, not the importance; the group ratio read
       2.147 against 1.073. The emission form's transpose is the adjoint
       flux (:ref:`characteristic-galerkin-system`).
   * - Each group's emission support read from :math:`\Sigma_t > 0` (the
       third rung's default)
     - A region transparent in a group that receives emission into it
       has no column, and the emission is dropped without a message: the
       absorption balance missed the source by
       :math:`2.2 \times 10^{-2}` (sphere) and
       :math:`2.0 \times 10^{-1}` (slab).
   * - One emission support shared by every group, the union
     - Columns whose emission is identically zero; behind a mirror a
       column in a region void in the group is a source on lossless
       trapped lines, refused with ``TrappedSource`` although the problem
       is well posed.
   * - The transport blocks handed to the system as a stored field, with a
       check of their supports
     - A block carries neither its group nor its :math:`\Sigma_t`, so two
       illegal states were admitted and answered silently wrong (the
       elegance review's V1, `[M]` 2026-10-07, a two-region two-group
       vacuum sphere at degree 2): blocks assembled on another mixture's
       :math:`\Sigma_t` gave :math:`k = 0.108657` against 0.107926, and
       the blocks in reversed group order 0.168. The system derives its
       blocks from its posing (:ref:`characteristic-galerkin-system`).
   * - The source questions answered on the k pencil,
       :math:`(W_s - SK)^{-1}` applied to the source
     - Its gain is the fission alone, so the least solution's
       subcriticality check cannot see a body made supercritical by its
       (n,2n) emission: the solve returns a finite, negative flux
       (:math:`-30` on the gates' mirror sphere). The source pencil puts
       every secondary emission in the gain and refuses it.
   * - The source-region mask group-major, ``(G, n)``
     - Every other ``(n, G)`` array of the package is region-major; with
       :math:`n = G` a transposed mask has the right shape and poses the
       wrong regions silently. The mask is region-major and a
       group-major one is refused by its shape when :math:`n \ne G`.
   * - The eigen flux of every eigen answer gauged by its fission
       production, ``Nearest`` included (the fifth rung's sketch, item 3)
     - A higher mode's net production is forced to zero on a closed
       homogeneous body (biorthogonality to the flat fundamental adjoint,
       with the fission rank one per point) and is not kept away from zero
       anywhere, so the gauge divides by rounding: the mirror slab's first
       spatial mode has a net production of
       :math:`-1.3 \times 10^{-16}`, and its left-half integral read
       :math:`-8.4 \times 10^{16}` under the gauge 100. A ``Nearest``
       answer reads its eigenvalue only (:ref:`characteristic-door`).
   * - The gauge :math:`\langle\nu\Sigma_f, \phi\rangle = 100` over the
       body (the sketch's Q3, "production's gauge")
     - The 100 is an infinite medium's production density, per unit
       volume; a finite body's gauge is a total. Gauging the total by 100
       made a closed body's integral equal the medium's density reading,
       and the mean flux of a body at a fixed physical state scale as
       :math:`1/V`. The gauge is the total production 1.
   * - The eigen gauge fixed inside the reference as the fission
       production (the 5a build)
     - It differed from S\ :sub:`N`'s production rate, which adds the
       (n,2n) emission, by a ratio nobody named; on the gates' ``_UP2N``
       sphere the two fluxes differ by 1.5455. A gauge is a choice of
       functional, so it is declared on the question
       (``Eigen.gauge``) and every reader resolves the same declaration.
   * - The fission's production recovered from the fission matrix
       :math:`\chi\otimes\nu\Sigma_f` by summing over the emitted group
     - Correct on the forward set only when :math:`\sum\chi = 1`, and
       wrong on the transposed set, whose fission matrix is
       :math:`\nu\Sigma_f\otimes\chi`: the sum returns
       :math:`\chi\sum\nu\Sigma_f` under the name :math:`\nu\Sigma_f`.
       `[M]` 2026-10-08, qa's adjoint fuel: forward production
       (0.125, 0.24) and spectrum (0.7, 0.3); the transposed set's
       production is the spectrum (0.7, 0.3), and the sum read
       (0.2555, 0.1095). The cross sections store the two factors and
       ``transposed()`` exchanges them.
   * - The detector table read as a rate, the response's answer the
       per-rate importance :math:`E^\dagger\psi^\dagger`
     - Self-consistent but named by nothing: ``Response`` is
       :math:`\psi^\dagger` and a flux integral reads the scalar flux,
       :math:`R\psi^\dagger = 4\pi E^\dagger\psi^\dagger`. A symbolic
       detector and a table then disagreed by :math:`4\pi`, and the first
       comparison with production's importance would have too (the class
       of ERR-004).
   * - A region-wise table projected by the quadrature rather than read
       at the nodes
     - Exact only to rounding and only from enough points: a projected
       constant was off by 5.74 at 1 point, which the resolution admitted.
       A nodal basis reads a constant's coefficients exactly
       (``PanelBasis.on_nodes``).
   * - ``Nearest(tau)`` searched over the whole spectrum of the
       :math:`k` pencil
     - The pencil carries a cluster of eigenvalues at :math:`k = 0` that
       the fission does not see whenever its production matrix is rank
       deficient; ``Nearest(0)`` served
       :math:`k = 7.7 \times 10^{-23}` with a flux of order
       :math:`10^{18}`, and the same :math:`\tau` was refused when the
       nearest residue happened to be complex. A mode is kept when its
       production exceeds the production matrix's rank tolerance.
   * - The role of a function as a boolean ``rate`` flag
     - The flag selected an arrow, the table branch ignored it, and the
       detector's arrow pointed the wrong way. The role is the arrow
       itself, the section, the pullback or none (the elegance review's
       C5).
   * - The point value read from the basis, :math:`\phi_h(x)` (the
       alternative in the fifth rung's sketch)
     - :math:`W\phi_h = Kq` makes :math:`\phi_h = P\mathcal Kq`, the L2
       projection of the transported emission, so its value at a point
       carries the basis's projection error, first order in the basis's
       resolution of the flux, whatever the accuracy of :math:`q`. The
       reading transports the emission once more,
       :math:`(\mathcal Kq)(x)`, the iterated Galerkin value
       (:ref:`characteristic-reading`).
   * - A separate direction rule at the point, over the box of
       ``Chart.directions_at``, split at ``DirectionDomain.tangencies`` and
       substituted there, with its own point count (the P1 sketch's
       alternative)
     - It re-derives the tangency, grazing and near-wall gradings in a
       second variable, so a grading added to one rule is missing from the
       other: the two-angular-rules hazard of #516 moved one level down.
       The point's directions are the line domain's lines through it,
       built by the line rule's constructors, and on a wall the two rules'
       lines are the same set bit for bit
       (:ref:`characteristic-reading`).
   * - The point's lines held as a ``LineRule`` (the first build's
       ``PointRule``, a copy of the line rule's fields)
     - Two copies of the guard, the sort and the chunking, already
       diverged at review (a chunk of 0 refused by one and admitted by the
       other, failing later), and a "block" computable on a point's
       measure (:math:`\mathbf 1^{\mathsf T}K\mathbf 1 = 1.49`, with no
       meaning). One weighted set, ``Lines``, with two roles
       (:ref:`characteristic-reading`).
   * - The rim law alone: the closure's pole graded at the outermost panel
       only, and the layer at the wall only (the rung's first build)
     - The pole sits at every interior tangency whose outer shells are
       thin or void, a distance :math:`(s(-\ln a) + \tau_{\rm out})/(2\Sigma_k)`
       from it (a void outer region at :math:`a = 0.99` read off by
       :math:`5.1 \times 10^{-4}` at the interface), and the cylinder's
       layer at every radius (an interface point off by
       :math:`1.0 \times 10^{-6}`). The law grades every impact-panel top
       (:ref:`characteristic-grading-law`, ERR-104, ERR-105).
   * - A floor on every grading distance in :math:`y`,
       :math:`\sqrt{2r\,\epsilon(r)}` times a margin over the first Gauss
       node's fraction (the first build's ``_B_MARGIN``), and a clamp
       :math:`\max(c - b, 0)` in the point's parameters
     - Patches on the extreme of one defect, a line chorded from its
       rounded :math:`b`: the floor kept nodes from rounding onto a panel
       top but not the half-chord's drift near it (a wall reading off by
       about :math:`3 \times 10^{-15}/(1 - a)`), it coarsened as the points
       grew, and on the cylinder it, not the rim law, set the grading. The
       line carries its exact level (:ref:`chart-and-chord-level`,
       ERR-106).
   * - The diffuse walls' update solved twice, once in the block and once
       in the reading
     - Two spellings of :math:`\alpha(I - T\alpha)^{-1}U^{\mathsf T}` that a
       change to one balance row would split. The currents are one value
       and the walls are folded by one map, read by both
       (:ref:`characteristic-wall-coupling`).

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
:eq:`characteristic-closure`, once their markers move. They moved in the
second rung's own commit (``ff979520``): the six † functions carry
``verifies("characteristic-traversal-integrals")`` and the ‡ function
``verifies("characteristic-closure")``, at ``l0``, as do the closure
file's three optical-depth functions (`[M]` 2026-10-07, ``git grep`` of
the markers in the two files).

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


The gates of the line rule, the assembly and the walls' coupling
-----------------------------------------------------------------

`[M]` 2026-10-07, ``pytest --collect-only`` on branch
``feature/characteristic-rung3`` at ``71a207fa``: the third rung's three
files hold 210 rows in 41 functions,
``tests/gates/derivations/test_characteristic_assembly.py`` 174 (34
functions, 14 rows ``slow``), ``tests/gates/geometry/test_line_domain.py``
30 (5) and ``tests/gates/geometry/test_measure_density.py`` 6 (2); the
seven files of the package and its kernel verbs hold 769 rows, 755 of them
outside ``slow``. The commit counts its additions as 46 new test functions
(223 cases), the third rung's own and the even basis's and the arriving
flux's rows in the earlier files, and 37 earlier rows re-posed for the
even basis and the retired density. Every new row is ``foundation`` except
``test_the_density_is_the_beam_density_times_the_folded_directions``, at
``l0`` under :eq:`geometry-measure-on-lines`; each row's planned level is
in its docstring, waiting on this page's labels: closed-body conservation,
the escape and the walls' closed forms, the reciprocity rows and the
operator rows ``l1`` under :eq:`characteristic-galerkin-assembly` or
:eq:`characteristic-boundary-resolvent`; the self-convergence ladders
(the centre piece, the impact and arc-length rules, the slab's cosine,
the cylinder's polar angle) ``l2``; the measure density's two rows and
Cauchy's formula ``l0`` under :eq:`geometry-measure-density` and
:eq:`geometry-line-domain`. The markers are the test-architect's to move.

The rows by what they gate, with the reference each compares against:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Rows
     - Reference
   * - closed-body conservation on twelve sphere and slab bodies at three
       cross-section columns and three cylinder bodies at one (``slow``,
       8 points), the
       white sphere at 24 points to :math:`10^{-14}`, each group through
       its own cross section, and the near-void sweep from
       :math:`\Sigma_t = 1` to :math:`10^{-12}`
     - the mass matrix applied to ones, :math:`W\mathbf 1`
   * - the escape probability two ways and the wall transmission, sphere
       and slab at :math:`\tau` from 0.01 to 1000, the cylinder ``slow``
     - Hébert, :math:`(1 - 2E_3)/(2\tau)`, Bickley's :math:`P_{ss}`, in
       mpmath
   * - the white and the specular laws, the void transmissions, the
       balance of each injected current, reciprocity two ways
     - one re-emission chain's balance; the per-line specular integral;
       the solid angle :math:`(r_0/R)^2`; the identity
       :math:`\sum T + \ell = 1`; :math:`R = UD^{-1}` and
       :math:`A_w T_{w'w} = A_{w'}T_{ww'}`
   * - an interface between equal materials, a transparent cavity, a void
       outer layer on the emission support, the full block of a void layer
       refused
     - the same body without the interface, with an inner mirror, without
       the layer; the refusal
   * - the slab's derived grazing against 12 fixed halvings, the near-void
       slab with open walls, thin graded panels and one block entry (both
       ``slow``), a small first region, a wide panel before a thin one, a
       tiny cavity
     - closed forms in mpmath; a block with 30 more halvings at 16 points;
       conservation
   * - the self-convergence ladders below the working point
     - the rule at 128 points, or the ladder's own monotone ratio
   * - the route of the one density, the cylinder's polar rule by
       structure, the budget as a partition, no diffuse wall adding
       nothing, symmetry declared blind with its arc-length teeth, the
       refusals
     - bitwise identities and refusal fragments

The battery: `[M]` 2026-10-07, the test-architect's
``scratch/characteristic_architecture/p1_step_b3/gates/battery/`` (plugin
``battery_plugin.py``, table ``battery_table.md``, each arm an in-process
textual mutant run over the files that can redden, ``-O``, ``slow``
deselected): the honest run 754 passed; 54 arms, 52 redden their target
rows, and 2 are declared blind. The arms, by what they break (reds in
their scope):

.. list-table::
   :header-rows: 1
   :widths: 58 12 30

   * - Arm
     - Reds
     - Where
   * - A2 the first group's coupling reused for every group; K1 the
       measure's derivative one power too many
     - 158; 157
     - every coupled block; every density reader
   * - A6 a chunk over the budget dropping a line; A16 the rim grading
       removed; K2 the density's constant dropped
     - 112; 106; 100
     - the assembly; every thick or rim-reading row; the mass and the walls
   * - I2 the arriving flux scaled; P2 the white wall read as specular;
       L5 the representative's foot along :math:`\hat e_x`
     - 87; 80; 69
     - the closure and the coupling; the lines
   * - A10 plain Gauss in :math:`b`; A12 the impact rule quartered; A14 the
       carried attenuation restarted at a region change
     - 59 each
     - conservation and the closed forms
   * - A1 :math:`D` dropped; L2 the sphere's fold :math:`2\pi`
     - 58; 57
     - every white wall; the sphere's density
   * - I1 the arriving flux on the wrong traversal; E2 the even panel off;
       L4 the slab's density without :math:`|\mu|`
     - 50; 37; 37
     - the closure; the even rows; the slab
   * - A3 a mirror read as diffuse; A5 the reversed traversals summed; E1
       the even panel keyed on any zero; P1 the pairing at the entry wall
     - 26; 26; 23; 21 (and 1 collection error)
     - the walls' rows; the even rows; the closure
   * - A11 the slab's grazing grading removed; L6 the box guard deleted;
       T3a the branch grading removed; A8 the quarter area of the other
       wall; N12 the subtraction's diagonal and no balance row
     - 15; 13; 13; 12; 12
     - the slab's closed forms; the refusals; the turning rows; the hollow
       white bodies; the near-void sweep
   * - E3 the mass rule at :math:`p + 2`; A20 the old ``chord_quadrature``
       rule in :math:`b`; A17 the hp toward :math:`b = 0` removed; A13 a
       void attenuating; N2 the balance row removed
     - 9; 9; 7; 6; 6
     - the even mass; the small radii and the thick sphere; the near-void
       sweep
   * - A18 12 fixed grazing halvings; I3, I4, L7, T17 (ERR-099's defect),
       L1, L3 the cylinder's polar bound :math:`\pi` or one
       :math:`\sin\theta`, K3, K4, T3b, A7
     - 5 to 2
     - their target rows
   * - A4 the lossless refusal deleted; A15 the zero update falling
       through; A21 plain Gauss in :math:`\theta`; A22 the normal ends
       dropped; A24, A25 the two cancellations; W1, W2 the wall guards
     - 1 each
     - their target rows
   * - A19 the hp toward the next radius removed; A23 the grazing halving
       stopped at :math:`\tau_{\min}`
     - 0 in scope; 1 each on their witnesses
     - A19 on the wide-then-thin row (:math:`1.0 \times 10^{-12}` against
       :math:`2.0 \times 10^{-9}`), A23 on the ``slow`` block-entry row
   * - N1 the diagonal formed as :math:`1 - \alpha T_{ww}` alone; T2 the
       inner exponential grading removed
     - 0, declared blind
     - N1 is masked by the balance row (with it removed too, N12 reds 12);
       T2's catcher is the second rung's thousand-mean-free-path row,
       outside the arm's scope



The gates of the Galerkin system
--------------------------------

`[M]` 2026-10-07, ``.venv/bin/python -O -m pytest -p no:cacheprovider``
on branch ``feature/characteristic-rung4`` (uncommitted):
``tests/gates/derivations/test_characteristic_system.py`` holds 80 rows in
22 test functions, 2 of them ``slow`` (the one-region cylinder, mirror and
white); the 78 others passed in 78 s. Each row rests on the third rung's
block rows (``rests_on``: closed-body conservation and each group's own
cross section). The closed forms are written in the test from the
``Mixture`` arrays, never through ``group_emission`` or the 0-D helpers
(``instrument-doctrine`` X4). The bands are the test-architect's
measurements (``scratch/characteristic_architecture/p1_step_b4/ta/``):

.. list-table::
   :header-rows: 1
   :widths: 34 40 26

   * - Rows
     - Claim and reference
     - Measured
   * - ``test_the_cross_sections_are_the_mixtures_tables``,
       ``test_the_emission_acts_node_by_node``
     - the per-region tables and their per-node layout, bitwise, against
       tables built in the test (an (n,2n) transfer, asymmetric
       :math:`\chi\otimes\nu\Sigma_f`, supports differing by group)
     - ``array_equal``
   * - ``test_the_emission_support_is_where_something_is_emitted_into_the_group``
     - the support against a hand table, one region per edge (transparent
       but receiving; emitting in one group only; scattering out of a
       group with nothing in; :math:`\chi = (1, 0)`)
     - ``array_equal``
   * - ``test_an_anisotropic_emission_is_refused_naming_the_region_and_the_order``,
       ``test_a_field_off_the_support_is_refused_and_answered_once_the_support_is_widened``,
       ``test_a_source_region_mask_in_group_major_orientation_is_refused``
     - the refusals, each with a positive leg
     - message fragments
   * - ``test_the_fundamental_mode_satisfies_its_pencil``
     - the pencil residual of the fundamental mode
     - :math:`\le 1.4 \times 10^{-14}`, band :math:`10^{-12}`
   * - ``test_a_closed_body_reads_k_inf_with_a_flat_flux_in_the_infinite_mediums_group_ratio``
     - seven closed bodies (mirror, white, periodic, hollow) times three
       mixtures (downscatter with :math:`\chi` in both groups; upscatter
       with :math:`\chi = (1, 0)`; three groups), and the cylinder
       (``slow``): :math:`k = \rho(A^{-1}F)`, a flat flux with the infinite
       medium's group ratio
     - :math:`k` within :math:`4.4 \times 10^{-14}`, flat within
       :math:`8.8 \times 10^{-12}`, ratio within :math:`1.9 \times
       10^{-14}`; bands :math:`10^{-12}`, :math:`10^{-10}`,
       :math:`10^{-12}`
   * - ``test_a_closed_sphere_reads_soods_k_inf_with_upscatter``
     - URRb-2-0-IN and URRc-2-0-IN, Sood's printed :math:`k_\infty`
       1.365821 and 1.633380 :cite:`SoodForsterParsons2003`
     - :math:`4.1 \times 10^{-7}`, :math:`6.6 \times 10^{-8}`; band
       :math:`5 \times 10^{-7}`
   * - ``test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux``
     - :math:`\phi = (\mathrm{diag}\,\Sigma_t - S - F)^{-1} q` with
       upscatter, (n,2n) and subcritical fission
     - :math:`3.5 \times 10^{-12}`, :math:`1.7 \times 10^{-14}`,
       :math:`8.4 \times 10^{-12}`; band :math:`10^{-10}`
   * - ``test_a_region_transparent_in_a_group_keeps_its_emission_into_it``,
       ``test_a_layer_void_in_a_group_under_a_mirror_is_answered``
     - the support's two edges, by the absorption balance
     - :math:`\le 5.1 \times 10^{-15}`; band :math:`10^{-12}`
   * - ``test_one_group_k_increases_on_nested_spaces``,
       ``test_one_group_k_is_below_one_at_soods_critical_size``
     - :eq:`characteristic-one-group-bound`, on nested degrees and against
       Ua-1-0-SP and Ua-1-0-SL
     - smallest step :math:`2.8 \times 10^{-8}`; smallest margin
       :math:`2.4 \times 10^{-6}`
   * - ``test_k_is_one_at_soods_two_group_critical_sizes``
     - PU-2-0, U-2-0 and URRa-2-0, slab and sphere: :math:`k = 1` within
       :math:`|\mathrm dk/\mathrm d(\mathrm{mfp})|` times half a unit of
       the printed critical size, the slope measured in the test
     - URRa-2-0-SL is the tightest, :math:`6.4 \times 10^{-7}` of
       :math:`7.2 \times 10^{-7}`
   * - ``test_the_source_question_returns_the_mode_from_its_own_fission_deficit``,
       ``test_a_source_in_a_body_supercritical_by_its_n2n_emission_is_refused``
     - the two pencils on one system; the (n,2n) refusal with the k
       pencil's negative direct solve as its control
     - :math:`3.4 \times 10^{-14}`; refusal
   * - ``test_the_closed_body_adjoint_is_the_infinite_medium_importance``,
       ``test_the_one_group_adjoint_is_the_forward_flux``,
       ``test_the_forward_and_adjoint_modes_are_biorthogonal``,
       ``test_a_detectors_reading_of_a_source_is_its_response_paired_with_the_source``
     - :eq:`characteristic-adjoint`: the closed-body importance by hand;
       the one-group adjoint equal to the forward flux (self-adjoint);
       biorthogonality of four modes in the production form; the reading
       :math:`\langle r, \phi(q)\rangle = \phi^{\dagger\mathsf T} W_s q`
       on 20 seeded pairs
     - bands :math:`10^{-12}`, :math:`10^{-10}`, :math:`10^{-10}`,
       :math:`10^{-12}`

**Sood's UAL-2-0 sizes are excluded.** At UAL-2-0-SL and UAL-2-0-SP the
reference converges to :math:`k - 1 = +7.364 \times 10^{-6}` and
:math:`+9.402 \times 10^{-6}`, stable to :math:`10^{-9}` over five
resolutions up to degree 6, six layers and 24 points: 2.9 and 7.5
half-units of the printed digit. Production S\ :sub:`N` (diamond
difference, Gauss–Legendre, 400 to 1600 cells and S\ :sub:`32` to
S\ :sub:`128`) converges at second order on the slab to the Richardson
limit :math:`+7.4 \times 10^{-6}`, beside the reference, so the printed
slab size is the likely error; the sphere has no second method yet. The
rows wait on #588 (`[M]` 2026-10-07, the test-architect, ``ta/m2b_ladder``
and ``ta/m8_sn_ual``). No tolerance was widened.

**The two-sided reading row has one-sided partners.** The reciprocity of
the source and the detector questions is blind to an error both share
(``vv-principles`` anti-pattern 37); the closed-body source row (a closed
form) and the closed-body adjoint row are its one-sided partners.

The gates of the door
---------------------

`[M]` 2026-10-08, the archivist's run, ``.venv/bin/python -O -m pytest
-p no:cacheprovider -m "not slow"`` on branch
``feature/characteristic-rung5a`` (uncommitted):
``tests/gates/derivations/test_characteristic_reference.py`` holds 95
rows in 36 test functions, 1 of them ``slow`` (the closed cylinder's
fixed source); the 94 others passed in 63 s. The rows are the
test-architect's specification
(``scratch/characteristic_architecture/p1_step_b5/spec.md``), each with
its ``rests_on`` into the fourth rung's rows and its own lower rows, and
the bands are the test-architect's measurements on the tree after the
review fixes (``ta/measure_5.log``, ``ta/measure_7.log``). Every closed
form is written in the test from the ``Mixture`` arrays or the
geometry's numbers, never through the door's objects; the reading row
says that it reads the door's own flux coefficients, since its claim is
the reading.

**After the declared gauge** (2026-10-08, the same day):
``test_a_heterogeneous_eigen_flux_produces_one_neutron`` was re-posed in
``fe977a90`` onto the declared production. On a sphere one of whose
materials carries (n,2n) emission, the default production (fission plus
(n,2n)) reads 1 to :math:`10^{-13}`, and the fission production alone
reads below 1 by more than :math:`10^{-4}` (0.991), as
:eq:`characteristic-door-gauge` requires.

.. list-table::
   :header-rows: 1
   :widths: 34 40 26

   * - Rows
     - Claim and reference
     - Measured
   * - ``test_the_cross_sections_carry_chi_and_production_and_transpose_by_swapping_them``
     - ``production`` is ``SigP`` bitwise; ``transposed()`` keeps the
       total, transposes the scattering, exchanges the spectrum and the
       production and transposes the fission
     - ``array_equal``
   * - ``test_a_per_region_constant_projects_to_its_value_at_every_node``,
       ``test_a_function_of_the_panel_space_projects_to_its_coefficients``,
       ``test_the_projection_of_a_smooth_function_converges_with_the_degree``,
       ``test_a_step_inside_a_panel_is_integrated_once_it_is_named``
     - :eq:`characteristic-projection` on a sphere, a slab, a hollow
       sphere and a cylinder: a constant returns its value; a random
       element of the panel space returns its coefficients (the function
       :math:`c` is not returned on the even panel); a smooth function's
       error falls by more than 3 per degree to below :math:`10^{-6}` at
       :math:`p = 6`; a step inside a panel is integrated once named
     - :math:`3.3 \times 10^{-15}`; :math:`1.9 \times 10^{-15}` (the even
       panel :math:`2.8 \times 10^{-2}`); from
       :math:`2.3 \times 10^{-2}` to :math:`2.9 \times 10^{-7}`, steps
       :math:`\ge 5.0`; :math:`4.0 \times 10^{-15}` named, 0.22 not
   * - ``test_the_door_refuses_at_construction_naming_what_it_refused``,
       ``test_the_refusal_fragments_are_disjoint``,
       ``test_a_served_question_constructs_and_assembles_nothing``,
       ``test_a_malformed_resolution_is_refused``,
       ``test_a_point_value_is_read_from_the_solved_emission_and_refused_on_a_mode_after_the_solve``,
       ``test_an_eigenvalue_of_a_source_answer_is_refused_before_any_solve``
     - the eleven refusals at construction, each by its type and its
       shortest message fragment, the fragments pairwise disjoint; every
       served question constructs with zero ``LineRule.of`` calls; the
       resolution's refusals in six rows (an untyped transport, fewer
       projection points than the mass rule, none, a real degree, layer
       count or point count); the source eigenvalue refused before any block is
       assembled; a point value read from the solved emission bitwise
       (``system.point_flux``) and refused on a ``Nearest`` answer after
       its solve (re-posed in rung 5b from the refusal it pinned before)
     - message fragments and call counts
   * - ``test_a_flux_integral_is_the_volume_integral_of_the_galerkin_flux``
     - :eq:`characteristic-door-pairing` against the reconstructed
       Galerkin flux integrated by numpy Gauss–Legendre, 24 points per
       panel, for a region-wise, a smooth and a stepped symbolic weight
     - :math:`2.2 \times 10^{-16}`, :math:`1.1 \times 10^{-16}`,
       :math:`1.1 \times 10^{-16}`; band :math:`10^{-13}`
   * - ``test_the_doors_k_is_the_fundamental_of_the_system_written_by_hand``
     - the door's :math:`k` equals the :math:`k` of a ``GalerkinSystem``
       posed by hand (mixtures in interval order, ids ``(1, 0, 2)``,
       walls as ``Wall`` tuples) on a slab and a sphere
     - bitwise
   * - ``test_a_closed_homogeneous_body_reads_the_exact_infinite_mediums_k_and_flux``,
       ``test_a_heterogeneous_eigen_flux_produces_one_neutron``
     - :eq:`characteristic-door-gauge`: four closed bodies times three
       mixtures against the exact infinite medium (rational arithmetic),
       :math:`k` and 100 times each group's integral; a heterogeneous
       vacuum sphere's production reads 1
     - :math:`k` within :math:`4.6 \times 10^{-14}`, flux within
       :math:`3.4 \times 10^{-14}`; the declared production 0 off; band
       :math:`10^{-12}`, :math:`10^{-13}`
   * - ``test_the_door_reads_k_one_at_soods_critical_sphere``
     - Sood's PU-2-0-SP at its printed critical radius
       :cite:`SoodForsterParsons2003`, posed through a ``BC.vacuum`` tag
     - :math:`|k - 1| = 7.0 \times 10^{-7}`; band
       :math:`3.6 \times 10^{-6}`
   * - ``test_the_mode_nearest_tau_is_the_mirror_slabs_first_spatial_mode``,
       ``test_tau_is_read_in_the_k_chart``,
       ``test_tau_near_the_fundamental_reads_the_fundamentals_k``,
       ``test_tau_zero_reads_a_genuine_mode_not_the_null_fission_cluster``,
       ``test_a_nearest_answer_reads_its_eigenvalue_and_refuses_a_flux_integral``
     - the mirror slab's first spatial mode against the closed form
       :math:`k(B)`, the dominant eigenvalue of
       :math:`(\mathrm{diag}(\Sigma_t/L(B)) - S)^{-1}F` with
       :math:`L(B) = \arctan(B/\Sigma_t)/(B/\Sigma_t)` (the transport
       kernel's Fourier transform) at :math:`B = \pi/a`; the :math:`k`
       chart; the fundamental; the null-fission cluster excluded; the
       flux of a ``Nearest`` answer refused
     - :math:`6.5 \times 10^{-8}` (PU2), :math:`5.0 \times 10^{-7}`
       (URRb), band :math:`10^{-6}`; ``Nearest(0)`` reads
       :math:`4.3 \times 10^{-3}`; refusal
   * - ``test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_through_the_door``,
       ``test_a_source_where_its_group_emits_nothing_widens_the_support_and_balances``,
       ``test_a_source_question_on_a_body_supercritical_by_n2n_is_refused_through_the_door``
     - the fourth rung's closed-body source row through the door
       (region integral over the volume against
       :math:`(\mathrm{diag}\,\Sigma_t - S - F)^{-1}q`, the cylinder
       ``slow``); a source widening the support balances its absorption;
       the (n,2n)-supercritical body refused for a source and a detector
     - :math:`\le 5.8 \times 10^{-15}` (cylinder
       :math:`2.2 \times 10^{-14}`, 36 s); :math:`1.1 \times 10^{-15}`;
       ``NoLeastSolution``
   * - ``test_a_symbolic_source_is_a_density_over_directions_and_a_symbolic_weight_or_detector_is_not``,
       ``test_a_regionwise_source_is_read_onto_the_nodes_whatever_the_projection_points``
     - :eq:`characteristic-door-lifts`: a constant symbolic source reads
       as the table :math:`4\pi q`, a symbolic weight as the table
       :math:`w`, a symbolic detector as the table :math:`d`; a table is
       read onto the nodes, bitwise at 8 and 10 projection points
     - :math:`2.2 \times 10^{-16}`, 0, :math:`2.2 \times 10^{-16}`;
       bitwise
   * - ``test_a_detectors_reading_of_a_source_is_the_sources_reading_of_the_detectors_importance``,
       ``test_a_closed_bodys_importance_is_the_infinite_mediums_adjoint_solve``,
       ``test_the_doors_response_is_the_forward_systems_adjoint_on_the_same_posing``
     - :eq:`characteristic-door-reciprocity` on three bodies, two seeded
       pairs each; :eq:`characteristic-door-response` on closed bodies
       against :math:`4\pi(\mathrm{diag}\,\Sigma_t - S - F)^{-\mathsf T}r`
       (the untransposed solve 1.5 off); the door's route against the
       fourth rung's adjoint pencil on the forward system
     - :math:`\le 7.8 \times 10^{-16}`; :math:`\le 6.2 \times 10^{-15}`;
       0; band :math:`10^{-13}`
   * - ``test_the_factory_returns_an_uncertified_reference_that_reads_a_ratio_as_a_quotient``,
       ``test_c1_*`` to ``test_c4_*``
     - the factory; a ``Ratio`` the quotient of its two readings, bitwise;
       content identity of the two resolutions and the derivation; the
       traced memo (an equal derivation reads with no spawn, another
       resolution is absent; a refusal crosses the memo with its type)
     - bitwise; counts

`[M]` 2026-10-08, the archivist's probe re-read the headline numbers on
the built code (``archivist_rung5a_numbers.py``): reciprocity
``Response(r)`` over ``FixedSource(q)`` is
:math:`4\pi = 12.566370614359174` to :math:`2.2 \times 10^{-16}` on the
first two seeded pairs of the vacuum sphere; the gauge and the
null-fission cluster as quoted in :ref:`characteristic-door`.

**The battery.** `[M]` 2026-10-08, the test-architect
(``ta/battery/run.sh``, ``ta/battery/battery_table.md``): 39 arms, each
installed in the pytest process by a plugin that rebinds one function to
a transformed source, over the file under ``-O -m "not slow"``; the
baseline is clean and 39 of 39 arms redden at least one row. The control
(the pairing scaled by 1.1) reddens 10 rows; the widest arms (the
pairing's groups reversed; the projection without :math:`W^{-1}`)
redden 26. Production files were ``diff -q`` identical to their pristine
copies after the run. Declared blind: ``test_tau_near_the_fundamental_reads_the_fundamentals_k``
cannot see the chart (``test_tau_is_read_in_the_k_chart`` does); the
reciprocity rows share the transport blocks, so an error in :math:`K`
symmetric under :math:`i \leftrightarrow j` is invisible to them (the
third rung's rows carry :math:`K`); a uniform scale of the pairing is
inside a ratio's stabiliser; a projection row on the slab cannot see the
measure density, which is 1 there. The memo rows redden under every arm
that rewrites ``evaluate``, because the traced memo pins the function's
source and a rewritten function has none: they are not catchers of
those arms (``vv-principles`` anti-pattern 17 (h)).

The gates of the reading at a point
-----------------------------------

`[M]` 2026-10-08, ``pytest --collect-only`` on branch
``feature/characteristic-point-reading``:
``tests/gates/derivations/test_characteristic_reading.py`` holds 175 rows
in 40 test functions, 15 rows ``slow`` (the cylinder rows, the closed
white bodies of C6 and C8, the slab's box integral);
``tests/gates/geometry/test_chord_level.py`` holds 12 rows in 7 functions,
none ``slow``; the door's file holds the re-posed wiring row. `[M]`
2026-10-08, the archivist's run, ``.venv/bin/python -O -m pytest -p
no:cacheprovider -m "not slow"`` over the two files and the door's row:
173 passed, 15 deselected, in 155 s. The
rows are the test-architect's specification
(``scratch/characteristic_architecture/p1_step_b5b/spec.md``, §4 and §8),
each with its ``rests_on`` (the point measure on the landed rungs, the
readings of a given emission on the measure, the block identities and the
angular flux on the readings, the door on the readings). Every reference
is written in mpmath or by hand in
``tests/gates/derivations/_characteristic_point_mp.py``, which shares no
line, no rule and no closure spelling with the code: the backward path in
the plane of the point and the direction (a disk billiard), its first leg
to the first wall and its period closed by the geometric sum; the slab's
:math:`E_2` image series; the cylinder's polar integral of the in-plane
diameter's closure and Bickley's :math:`\mathrm{Ki}_2` in its Struve form;
the white walls' escape and transmission probabilities. Seven functions
carry ``verifies("characteristic-quadrature")``: the singular points, the
cylinder's axis, the general points, the cylinder's general points, the
partial mirror's wall and its neighbourhood, the void or thin outer
shells, and the cylinder's interface.

.. list-table::
   :header-rows: 1
   :widths: 34 40 26

   * - Rows
     - Claim and reference
     - Measured
   * - ``test_the_point_rule_integrates_the_direction_moments``,
       ``test_on_the_outer_wall_the_point_rule_has_the_line_rules_lines``,
       ``test_a_line_set_carries_levels_iff_its_chart_is_radial``,
       ``test_a_point_outside_the_body_is_refused``
     - the point measure (total 1, the second and fourth moments of
       :math:`\Omega\cdot\hat r`, :math:`\Omega_z^2` on the cylinder, the
       point on each line) at 14 points; one constructor on a wall; the
       levels required iff radial; the point outside refused
     - worst :math:`2.2 \times 10^{-15}`, band :math:`10^{-13}`; bitwise;
       message fragments
   * - ``test_the_reading_at_a_singular_point_is_the_closed_form``,
       ``test_the_reading_on_a_cylinders_axis_is_the_polar_integral``
       (``slow``)
     - the sphere's centre under :math:`a = 0, 0.6, 1` and the slab's faces
       under two pairs of albedos, two groups with :math:`q` and
       :math:`\Sigma_t` jumping at both interfaces; the cylinder's axis
     - centre :math:`\le 1.4 \times 10^{-16}`, faces
       :math:`\le 4.4 \times 10^{-16}`; band :math:`10^{-12}`
   * - ``test_the_reading_at_a_general_point_is_the_mpmath_route``,
       ``test_the_reading_at_a_general_point_of_a_cylinder_is_the_mpmath_route``
       (``slow``)
     - the three-region sphere at five points under three laws, the hollow
       sphere at four points under four, the slab at three; the cylinder
       under vacuum and a partial mirror
     - :math:`6.7 \times 10^{-14}` (sphere, :math:`x = 1.5`),
       :math:`5.4 \times 10^{-14}` (hollow); band :math:`10^{-12}`
   * - ``test_the_reading_on_and_near_a_partial_mirror_is_the_mpmath_route``,
       ``test_the_blocks_total_under_a_partial_mirror_is_the_closed_form``,
       ``test_on_a_cylinders_partial_mirror_the_block_and_the_wall_reading_are_the_closed_forms``
       (``slow``)
     - on and a millionth inside a partial mirror at 0.6, 0.9 and 0.99,
       the hollow sphere's outer wall, the slab :math:`10^{-6}` and
       :math:`10^{-3}` from a face; the homogeneous sphere's and cylinder's
       :math:`\mathbf 1^{\mathsf T}K\mathbf 1` and the cylinder's wall
     - :math:`\le 2.8 \times 10^{-14}` at 8 line points; block
       :math:`\le 5.1 \times 10^{-16}`; cylinder :math:`1.6 \times 10^{-15}`,
       wall :math:`4.3 \times 10^{-14}` and :math:`3.4 \times 10^{-13}`
   * - ``test_a_void_or_thin_outer_shell_under_a_near_one_mirror_is_read_to_the_bar``,
       ``test_the_blocks_total_under_a_void_or_thin_outer_shell_is_the_closed_form``,
       ``test_a_cylinders_interface_point_is_read_to_the_bar`` (``slow``),
       ``test_a_small_cylinders_layer_is_graded_outside_slow``
     - qa's F1 and F2 under the grading law: the void and thin outer
       shells under :math:`a = 0.99`, read and totalled; the three-region
       cylinder's interface; a two-region cylinder's interface and wall
       (the wall ``slow``)
     - :math:`7.1 \times 10^{-15}`, :math:`2.2 \times 10^{-16}`,
       :math:`4.1 \times 10^{-14}`, :math:`6.7 \times 10^{-16}`; block
       :math:`1.1 \times 10^{-16}` and :math:`2.2 \times 10^{-16}`; small
       cylinder :math:`4.5 \times 10^{-14}` (interface),
       :math:`4.8 \times 10^{-15}` (wall)
   * - ``test_the_wall_reading_as_the_mirror_tends_to_one_holds_the_bar``,
       ``test_a_direction_grazing_the_wall_is_read_through_its_level``
     - qa's F3 under the level: the wall at :math:`a = 0.999` and 0.99999;
       :math:`\mu = \pm 2^{-40}` on the wall
     - :math:`1.1 \times 10^{-15}`, :math:`4.4 \times 10^{-16}`;
       :math:`4.4 \times 10^{-16}`, :math:`7.8 \times 10^{-16}`
   * - ``test_a_point_outside_the_body_is_refused_with_one_message``,
       ``test_the_slabs_grazing_direction_below_its_resolution_is_refused``,
       ``test_a_point_near_the_centre_reads_the_centre``,
       ``test_a_point_below_the_smallest_normal_float_is_refused``
     - qa's F5: one message for both readings; the slab's grazing refused;
       :math:`c` from :math:`10^{-12}` to :math:`10^{-150}` reads the
       centre; below :math:`\sqrt{\text{tiny}}` refused
     - fragments; :math:`5.0 \times 10^{-14}`
   * - ``test_the_reading_of_an_emission_outside_the_basis_converges_in_the_degree``
     - an emission outside the basis (:math:`1.5 + \cos 3c` per region,
       L2-projected) read at three points converges monotonically in the
       degree, 1 to 5 against 6; ``l2``, no value leg (the general-point
       rows carry the value)
     - a ratio of :math:`10^{-3}` overall
   * - ``test_the_volume_integral_of_the_reading_is_the_galerkin_block``,
       ``test_the_region_transfers_read_through_points_are_reciprocal``,
       ``test_an_interface_between_equal_materials_is_invisible_to_the_reading``
     - C6 under vacuum, a partial mirror and white walls (the white bodies
       ``slow``); C8; the middle region split in two and three equal
       regions, readings at five points unchanged
     - :math:`1.2 \times 10^{-12}` to :math:`1.3 \times 10^{-12}` of
       :math:`\max|K|`, band :math:`10^{-11}`; band :math:`10^{-11}`;
       band :math:`10^{-13}`
   * - ``test_the_block_folds_its_walls_through_the_currents``,
       ``test_the_reading_behind_white_walls_is_the_escape_closed_form``
     - the one fold, bitwise, on four white bodies and none; white walls
       against the escape closed forms at seven points
     - bitwise; band :math:`10^{-12}`
   * - ``test_the_angular_flux_is_the_backward_paths_integral``,
       ``test_the_angular_flux_integrated_over_the_box_is_the_reading``,
       ``test_the_grazing_limit_on_the_wall``,
       ``test_at_the_centre_the_angular_flux_is_isotropic_and_four_pi_of_it_is_the_reading``
     - :math:`\psi` against :math:`\psi_{\rm disk}/4\pi` at four points; the
       box integral (the slab ``slow``); the grazing limit at
       :math:`|\mu| = 10^{-6}`; the centre
     - band :math:`10^{-12}`; :math:`4.4 \times 10^{-16}`,
       :math:`1.2 \times 10^{-14}` (wall), :math:`8.9 \times 10^{-16}`
       (slab); :math:`\le 4.4 \times 10^{-16}` under a mirror and
       :math:`1.6 \times 10^{-4}` under 0.6, within the band
       :math:`10\epsilon/\mu^2` the row was written with; band
       :math:`10^{-14}`
   * - ``test_a_pure_absorbers_point_values_are_the_mpmath_route``,
       ``test_garcias_case_1_per_point``
     - E4: through the door, the vacuum and the partial-0.6 absorber
       sphere at three points; E1: Garcia's Case 1 at 15 radii
       :cite:`Garcia2021`, half the published flux, within the published
       rounding plus :math:`6 \times 10^{-5}` (ten times the ladder's
       step), ``catches("ERR-090")``
     - band :math:`10^{-12}`; at the finest rung every point within the
       rounding alone (worst :math:`r = 3.0`: :math:`3.0 \times 10^{-5}` of
       :math:`3.5 \times 10^{-5}`); the trajectory-resolvent family's
       worst was :math:`3 \times 10^{-3}` and :math:`4 \times 10^{-2}`
   * - ``test_a_closed_bodys_fundamental_reads_the_infinite_medium_flux_in_the_declared_gauge``,
       ``test_the_fundamentals_readings_carry_the_declared_production``,
       ``test_a_closed_body_with_a_uniform_source_reads_the_infinite_medium_flux_at_every_point``,
       ``test_a_detectors_importance_read_at_points_is_reciprocal_to_the_forward_flux``,
       ``test_in_one_group_the_importance_is_four_pi_times_the_forward_flux``,
       ``test_a_point_value_of_a_mode_nearest_tau_is_refused_after_the_solve``,
       ``test_a_ratio_of_point_values_is_the_quotient_and_gauge_free``
     - D1 to D7: the door's point values against closed forms (the
       gauged infinite medium at the centre, inside and on the wall; the
       declared production; a uniform source; reciprocity through
       readings; the one-group :math:`4\pi`), the ``Nearest`` refusal, a
       ``Ratio``
     - bands :math:`10^{-12}`, :math:`10^{-9}` (measured
       :math:`1.9 \times 10^{-11}`), :math:`10^{-12}`, :math:`10^{-9}`
       (measured :math:`1.2 \times 10^{-10}`), :math:`10^{-13}`
   * - the seven rows of ``test_chord_level.py``
     - the image's exact level (:ref:`chart-and-chord-level`)
     - bitwise; :math:`10^{-15}` of :math:`y`; fragments

**The battery.** `[M]` 2026-10-08, the test-architect
(``ta/battery/run.sh``, plugin ``battery_plugin.py``, verdicts
``ta/battery/verdicts.md``): each arm an in-process textual mutant rebound
in every module binding, over the reading's and the level's files and the
door's re-posed row outside ``slow`` (173 rows), production files
``diff -q`` identical to pristine copies after the run. On the final tree
the baseline is clean and every arm reddens its target rows; the 15
``slow`` rows run unmutated are green. The positive control (each group's
row built on the other group's :math:`\Sigma_t`) reddens 85 rows. The arms
by what they break, with their reds: the sphere's Jacobian without its
:math:`b` (83), one side dropped (82), the cylinder's polar weight without
:math:`\sin\theta` (4); the level dropped by ``chord`` (9: the level's
rows, the wall at :math:`a \to 1`, the grazing directions, the wall's box
integral), the default level moved off :math:`(b, 0)` (4), a tangent line
counted as crossing (2), the agreement check deleted (2), the levels
unchecked on ``Lines`` (3); the grading law at the outermost panel top
only (7), its layer term removed (1); the outside refusal split (2), the
slab's grazing unrefused (2), the subnormal point unrefused (4); the
currents without :math:`\alpha` (5), the diffuse part dropped from the
reading (7), the block's fold bypassed by one ulp (4); the angular flux's
side from :math:`\Omega_y` (7) and without its :math:`1/4\pi` (11); the
panel tangencies lost (5), a second constructor (2); the door's emission
unscaled (4), the source emission alone (19), the response untransposed
(1), the detector lifted by the section (4), the door's point scaled by
:math:`1 + 10^{-7}` (11). Declared blind: the panel tangencies' arm is
caught only by the one-constructor and the equal-materials rows; the
response untransposed is invisible in one group by construction; an
under-estimate of :math:`\tau_{\rm out}` in the grading law only grades
more finely (qa). Not armed: the axial image's parameters, the reading's
parameter re-spelt (each pinned bitwise by a ``test_chord_level.py``
row), the equal-materials transport, the ``Nearest`` point refusal.

`[M]` 2026-10-08, the archivist's re-drops on the final tree (the same
plugin, ``-O``, the reading's and the level's files outside ``slow``,
172 rows): the grading law reduced to the next radius alone, the third
rung's rule, reddens 21 rows; its pole term removed reddens 18 (the
partial mirror's walls and blocks at 0.9 and 0.99, the void and thin
outer shells, the walls at :math:`a \to 1`, the wall's box integral, the
partial absorber through the door); its layer term removed reddens both
rows of ``test_a_small_cylinders_layer_is_graded_outside_slow``, the
``slow`` wall row included (run alone, 37 s).

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
- **A white wall's return is not in the line closure.** The line part
  treats a diffuse wall as absorbing (its specular amplitude is 0), and
  what it returns enters through the white walls' coupling, as an
  arriving flux on the lines that leave it
  (:ref:`characteristic-wall-coupling`). Reading
  ``GroupTransport.line`` as the block of a white-walled body drops the
  walls; read ``block``.
- **The block is rectangular.** Its columns are the emission support
  (``GroupTransport.support``), its rows every panel; apply it to the
  coefficients on the support, ``block @ x[support]``.
- **The rule is per group.** ``LineRule.of`` takes one group's
  :math:`\Sigma_t`, because every grading reads it; a rule built for one
  group is under-resolved for a thicker or a near-void one.
- **The lines are ordered by projected speed, not by coordinate.** A
  caller that wants the :math:`(b, \theta)` grid back sorts
  ``LineRule.lines.coordinates`` itself (the polar-rule gate does); the
  levels are permuted with them (``Lines.ordered``), so sort the three
  arrays by one key.
- **Conservation cannot see the rule over lines.** It holds line by line
  on a closed body; gate a change to a line rule against the walls'
  closed forms, not against conservation (ERR-101).
- **Totals cannot see one entry.** A grading that is exact on every
  total can leave one block entry off by :math:`10^{-6}` (qa's N1); a gate
  of a rule over lines compares entries.
- :math:`R` **is per nominal unit current.** It is the response to the
  injection :math:`1/D`, not to the current the rule actually injects; do
  not renormalise it by the injected tally, which cancels the wall's area
  out of the block.
- **A wall is specular or diffuse.** ``Wall`` refuses both at once, and
  one returning more than it receives.
- **Subtract in mpmath.** A near-void reference such as
  :math:`1 - P_{\rm esc}` evaluated in mpmath and subtracted in double
  precision reports a defect the code does not have.
- **Size a cylinder fixture before adding it.** A white cylinder block
  costs, per group at 8 points, 15.6 s at :math:`\tau = 0.01`, 27.8 s at
  :math:`\tau = 30` and 175 s on the gates' three-region body; at 16 points,
  102 s at :math:`\tau = 30` and 792 s on the three-region body (`[M]`
  2026-10-07). Most of it is the inner rule's basis evaluations. A change
  that pads an integral to the batch's thickest stretch again multiplies
  that part 4 to 6 times and reddens no correctness gate, because the
  values do not change (:ref:`characteristic-cylinder-cost`).
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

- **The pencil's vector is the emission, not the flux.** An eigenvector of
  ``GalerkinSystem.pencil`` lives on the emission space, one coefficient
  per support node of each group; read the flux with ``flux(vector)`` and
  one group's coefficients with ``emission.split(vector)``. The adjoint's
  vector and ``response`` are the adjoint flux on the same space.
- **Source questions go to the source pencil.** A fixed source solved on
  the k pencil's loss is not refused when (n,2n) alone makes the body
  supercritical; it returns a negative flux. Use ``fixed_source``.
- **A source must sit on its group's support.** ``fixed_source`` refuses a
  coefficient outside it; pose the system with those regions in
  ``source_regions``. The mask is region-major, ``(n, G)``, as the cross
  sections; with as many regions as groups a transposed mask cannot be
  detected by its shape.
- **The detector is paired with the flux once.** ``response`` takes the
  detector's nodal coefficients :math:`r` with the load
  :math:`K^{\mathsf T}r`, and the reading is
  :math:`r^{\mathsf T}W_G\phi = \phi^{\dagger\mathsf T}W_s q`; a load
  :math:`K^{\mathsf T}W_G r` applies the metric twice.
- **The one-group bound has no margin on a graded basis at a printed
  critical size.** The graded basis converges below the truth's digits
  (the half slab reads :math:`+4.3 \times 10^{-10}` against a resolution
  of :math:`3.6 \times 10^{-9}`); test the bound on ungraded panels.
- **A table and a symbolic function of one role enter differently.** A
  ``RegionwiseConstant`` is a value on the angle-integrated space and is
  lifted by its role's arrow; a ``Symbolic`` source or detector is a
  density over the directions and is retracted, so a constant symbolic
  source :math:`q` is the table :math:`4\pi q`, while a symbolic weight
  is the table :math:`w` as given.
- **A response reads the adjoint scalar flux.** ``Response(r)`` read by
  ``FluxIntegral(q)`` is :math:`4\pi` times ``FixedSource(q)`` read by
  ``FluxIntegral(r)`` for two tables, and :math:`4\pi` times the fourth
  rung's ``system.response(r)`` paired with :math:`W_s q`, which is the
  importance per unit emission rate.
- **The door's cross sections are the adjoint ones for a response.**
  ``CharacteristicDerivation.cross_sections`` is
  ``RegionCrossSections.of(...).transposed()`` when the question is a
  ``Response``: its ``production`` is the forward :math:`\chi`.
- **The gauge is the question's.** The eigen flux is scaled so that the
  declared production is 1; with no gauge declared that production
  counts the (n,2n) emission as well as fission. To compare with a reader
  that scales by fission alone (today the homogeneous solver, until
  #517), declare ``Eigen(..., gauge=CellCoefficient.every(Channel.FISSION_EMISSION))``;
  on a body with (n,2n) the two fluxes differ by the ratio of the two
  productions.
- **A** ``Nearest`` **answer has no flux.** Read its ``Eigenvalue``; a
  flux integral of it is refused even when :math:`\tau` lands on the
  fundamental.
- **A jump inside a panel is split where** ``Symbolic.steps`` **finds
  it.** ``project`` cuts the load's rule at those points (the sign
  changes of a ``Piecewise`` condition, a ``Heaviside``, an ``Abs``, a
  ``sign``, a ``Max`` or a ``Min``), and ``steps`` refuses any other
  non-smooth head by name; a jump at a transcendental root is outside its
  scope.
- **A point value is not the basis value.** ``PointValue`` reads
  :math:`(\mathcal Kq)(x)`; the flux coefficients reconstructed at
  :math:`x` are :math:`\phi_h(x) = (P\mathcal Kq)(x)` and differ from it by
  the basis's projection error. Compare a point reading with a point
  reference, and a flux integral with the pairing.
- **The two roles of a line set are not interchangeable.** A
  ``PointRule``'s weights are :math:`\mathrm d\Omega/4\pi` at one point,
  not the line measure: a block assembled on its lines has no meaning, and
  ``PointRule`` has no ``transport``. Build a block with ``LineRule.of``
  and a point's lines with ``PointRule.of``; both read the same
  one-dimensional rules.
- **A radial line near a tangency must carry its level.** A line built
  from :math:`b` alone loses its half-chord below
  :math:`\sqrt{2r\,\epsilon(r)}`; ``Lines`` refuses a radial set without
  levels, and a line built outside the rules is chorded with
  ``chord(line, level=(radius, half_chord))``
  (:ref:`chart-and-chord-level`).
- **The angular flux is per steradian; the lines carry** :math:`4\pi`
  **times it.** ``angular_flux`` divides by
  :data:`~orpheus.derivations.continuous.characteristic.closure.FULL_SOLID_ANGLE`;
  a functional of the lines' transport is :math:`4\pi` times a flux per
  steradian, and their weights hold the :math:`1/4\pi`.
- **The white walls enter a reading only through their currents.** A
  functional of the stacked sources is folded onto the emission by
  ``WallCoupling.on_emission``; reading its emission part alone is the
  reading of a body whose diffuse walls absorb.
- **A point within** :math:`1.5 \times 10^{-154}` **of the centre is
  refused**, other than the centre itself: read the centre, which is
  the limit to :math:`5 \times 10^{-14}` from :math:`10^{-12}` down.
- **Size a cylinder point before adding it to a gate.** A point reading
  on a cylinder costs seconds (:ref:`characteristic-reading`,
  "Performance"; #591), and the grading law's cost grows with the number
  of panel tops; the sphere's costs hundredths of a second.
- **The grading law cannot be gated against an under-estimate of**
  :math:`\tau_{\rm out}`. Dropping it only grades more finely, so no
  value gate reddens; a change that over-estimates it (or the speed) is
  what the void and thin outer shells and the small cylinder catch.

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
   * - 2026-10-07
     - The third rung: one group's transport block. The panel basis became
       even at a singular stratum (Schwarz; the user's ruling over grading
       the impact rule toward the centre), its density moved to the
       kernel's ``measure_density``, and the line closure took an arriving
       flux. ``LineRule`` assembles the block by Galerkin over the lines of
       the kernel's new line domain (the cylinder's coordinate the polar
       angle, by ruling), on the emission support, and ``WallCoupling``
       adds the white walls through
       :math:`R\,\alpha\,(I - T\alpha)^{-1}U^{\mathsf T}`, with the area
       factor :math:`D` the ruled spelling had omitted, and with
       :math:`I - T\alpha` formed from the loss and a balance row after
       the elegance review measured the subtraction's
       :math:`10^{-4}` near void (ERR-102). Conservation was re-posed as
       the alarm and symmetry declared blind. qa found the line rule's
       fixed resolutions blind on regimes the gates did not reach, and
       the user ruled every grading derived from the optical scale
       (ERR-101; the normal direction, ERR-103); the gradings moved into
       one module. The thick cylinder's cost was ruled the next step
       (#586).
     - ``71a207fa``
     - #405, #584, #585, #586
   * - 2026-10-07
     - The cylinder's cost (#586) was the Volterra triangle's inner rule,
       not the line count of the tensor rule: every attenuated integral
       was padded to the batch's thickest stretch, 4.2 to 6.1 times the
       live evaluations. Each line's live entries are now packed first
       (one packing shared by the pieces and the inner rule), and the
       basis reads Lagrange tables built once. The rule is unchanged; the
       :math:`\tau = 30` white block went from about 215 s to 27.8 s per
       group. The non-tensor rule the issue proposed was not needed, and
       the user restored the thick cylinder legs and the three-region
       closed-body rows to the slow tier.
     - ``95d1a511``
     - closes #586
   * - 2026-10-07
     - The fourth rung: the multigroup Galerkin system. The emission
       matrices were split out of the 0-D pair into one assembly site,
       ``group_emission``; ``RegionCrossSections`` reads them from one
       mixture per region and refuses anisotropy there; each group's
       emission support is exact (the user's ruling over each group's
       :math:`\Sigma_t > 0` and over one shared support).
       ``GalerkinSystem`` derives its blocks from its posing (the elegance
       review measured blocks handed in answering silently wrong) and
       carries the k pencil and the source pencil on the emission space.
       The unknown was first the flux; the main agent measured its
       transpose to be the adjoint collision density, and the user ruled
       the emission density the unknown the same day. Sood's UAL-2-0
       critical sizes were found off beyond their digits (#588).
     - ``42222c47``
     - #405, #588
   * - 2026-10-07
     - The fifth rung's first half (5a): the door. ``Resolution`` holds
       the solve's and every reading's resolution; ``PanelBasis.project``
       projects a symbolic function on the mass matrix's own rule and
       ``on_nodes`` reads a table exactly; ``CharacteristicDerivation``
       refuses at construction and answers ``Eigenvalue`` and
       ``FluxIntegral`` by pairing, with ``PointValue`` refused for 5b (the
       user's ruling on the split). A ``Response`` is the forward problem
       on the transposed cross sections, which now store :math:`\chi` and
       :math:`\nu\Sigma_f` as two factors. The reviews found that a
       production gauge divides by rounding on a higher mode; the user
       ruled a ``Nearest`` answer eigenvalue-only, corrected the sketch's
       gauge (a density of 100) to the total production 1, and ruled the
       response's answer the adjoint scalar flux :math:`R\psi^\dagger`,
       its :math:`4\pi` the detector's role arrow. ``Nearest`` excludes the
       null-fission cluster.
     - ``fe977a90``
     - #405, #529, #566
   * - 2026-10-08
     - The eigen gauge became part of the question. The 5a door had fixed
       it at the fission production and called that S\ :sub:`N`'s, whose
       production rate adds the (n,2n) emission. The user ruled the gauge
       declared: ``Eigen.gauge``, resolved by the specification, with the
       default ``EIGEN_GAUGE`` (fission and (n,2n)) and no gauge where
       nothing produces; ``production_emission`` gives what each channel
       emits on the reference side, in the references' own (n,2n) literal.
       The door, the exact infinite medium and the trajectory resolvent
       read it; the production solvers follow in #517.
     - ``fe977a90``
     - #405, #517
   * - 2026-10-08
     - The fifth rung's second half (5b): the reading at a point.
       ``PointValue`` reads the transported emission
       :math:`(\mathcal Kq)(x)`, the iterated Galerkin value, not the
       basis value (the user's ruling on the sketch); the directions at
       the point are the line domain's lines through it, from the line
       rule's constructors with the point as one more end, weighted by
       :math:`\mathrm d\Omega/4\pi`; the answers carry their emission and
       ``_refuse_point_reading`` retired; the angular flux at a point is
       exposed with no observable. The reading found two grading defects
       in the third rung's line rule, the partial mirror's closure pole
       and the cylinder's polar-speed layer, graded first by a rim law at
       the outermost panel; qa's review found both at interior tangencies
       and the half-chord lost by storing :math:`b`, and the user ruled the
       grading law at every impact-panel top and the image's exact level
       into the rung (ERR-104, ERR-105, ERR-106; the rim law, its floor
       and margin, and the tangent clamp retired). The elegance review
       made the line set one value with two roles (``Lines``,
       ``LineRule``, ``PointRule``) in a module of its own,
       :mod:`~orpheus.derivations.continuous.characteristic.lines`, and the
       white walls one fold of their currents (``WallCoupling.update``
       retired onto ``currents`` and ``on_emission``). The slab's grazing
       refusal is #585's; the cylinder's extra cost is accepted until #587;
       the point's cost against the spec's target is #591.
     - ``b76b9a9d``
     - #405, #585, #587, #590, #591
   * - 2026-10-08
     - The cylinder's rule over :math:`(b, \theta)` became iterated
       (#587, its first half): the grading law's speed-free data per
       impact panel became one value, ``ImpactPanels`` (its ``of``,
       ``distances`` and ``rule``), and ``impact_per_polar`` builds each
       polar angle's impact rule at its own projected speed, so the
       ``slowest`` argument, ``tangency_distances`` and ``impact_by_polar``
       retired. On the ``ABA`` cylinder the block went from 308 480 lines
       and 2296 s to 102 320 lines and 1275 s, the same k to
       :math:`8.3 \times 10^{-15}`; the sphere and the slab are
       bit-identical. Grading the polar angle per impact panel, #587's
       other half, was left open.
     - ``200b4233``
     - #405, #587
   * - 2026-10-09
     - P1 step (c), the corroboration: code-to-code rows (L4, no level
       marker, deleted with the old family at step (e2)) showing the new
       reference inside the trajectory resolvent's own error estimate on
       every S\ :sub:`N` row's problem and on Garcia 2021, with static and
       run-time legs asserting that the old side executes nothing of this
       package (``tests/gates/_corroboration.py``, shared with P0's kernel
       rows). Just before it, grading the cylinder's polar angle per impact
       panel (#587's second half) was refuted for cutting the ``ABA``
       cylinder's lines.
     - ``17708ac8``
     - #405
   * - 2026-10-10
     - P1 step (d): the 13 S\ :sub:`N` rows (17 ids) re-pointed onto this
       reference (:ref:`characteristic-sn-rows`). Its error estimates live in
       ``tests/gates/derivations/_characteristic_ladders.py`` (geometric over
       the joint ladder; the A|B|A sphere's shape over pairs of degrees; the
       cylinder over a panel ladder plus the transport step, with the
       largest ratio any ladder showed, qa's F3), the tolerance rules in
       ``_ladder_rules.py``. Statement (ii) held on every row and every
       tolerance tightened or held; the cylinder RECORD's reference keys
       were re-baselined and its S\ :sub:`N` keys did not move. The rows'
       docstrings stopped claiming full independence: both sides read the
       partial wall's ``SpecularReturn.kernel`` (qa's F1), pinned outside
       the comparison. ERR-094's corner was measured invisible to its band at
       the sphere fixture (qa's F2). The labelled equation
       :eq:`sn-curvilinear-characteristic-reference-crosscheck` was rewritten for
       this reference (it kept its old name, naming the retired family,
       until step (e2) renamed it).
     - ``d9425977``
     - #405, #566
   * - 2026-10-10
     - #592: the solve became a traced memo of its own, ``solve``, keyed on
       the derivation alone, and ``answer`` the per-instance hold of its
       value, so every reading records the solve as one child entry and
       the observables of one derivation share one solve
       (:ref:`characteristic-door`). The exact key began to order a set by
       its elements, after qa measured one eigen gauge's solve keyed two
       ways across hash seeds (:ref:`verification-reference-cache-key`).
       The ``implements`` edge of :eq:`characteristic-door-gauge` moved
       from ``answer`` to ``solve``.
     - ``45405f86``, ``23e7711e``
     - #405, #592
   * - 2026-10-10
     - P1 step (e): the trajectory-resolvent family retired. Step (e1a)
       moved its SymPy derivations to
       :mod:`orpheus.derivations.continuous.characteristic.origins`; step
       (e1b) re-posed its 123 kept rows on this reference, with catchers
       for ERR-034, ERR-035 and ERR-091 and the eleven labels' verifiers;
       step (e2) deleted its 12 numeric modules and 22 test files in one
       commit (:ref:`characteristic-successors`).
     - ``f66c1c45``, ``44303919``; (e2) *(in development)*
     - #405
