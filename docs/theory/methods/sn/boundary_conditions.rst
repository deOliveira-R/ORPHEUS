.. _boundary-conditions:

Boundary Conditions
===================

Infrastructure
--------------

Boundary conditions are declared in two places and joined in a third. The
**geometry** says which faces of the body exist and what each one *is*: a
face the domain continues across (a **deck**: a mirror, or the translation
wrap of a periodic pair) or a face where the domain ends (an integer
**tag**). The **problem's caller** says what each tag *means*, by a
``BoundaryConditions`` mapping
from tag to condition (``orpheus.transport.boundary_conditions``). The **SN Problem**
(:class:`~orpheus.sn.problem.SNProblem`) joins the two at its construction and
realizes the result as one operator per face. The split exists so that
one geometry serves every method and every exterior: the geometry never
states what a tag means (a tag is a declaration, not an observation, the
exact twin of a material id), and a method never reads a tag's content, only the
condition it is bound to. The method-agnostic theory (the affine trace map,
the cone of walls, the deck family, the colouring stabiliser, the admission
lattice) is :ref:`bc-architecture-of-record` on the foundations page, and
each retired name of the trace-law architecture maps to its successor in
:ref:`bc-wave2-crosswalk`; this page is the SN half.

**Stage 1 --- Geometry declares face roles.**
The roles are declared on the
:class:`~orpheus.geometry.structured_geometry.StructuredGeometry` (its
``face_roles``), one per boundary point, and the
:class:`~orpheus.mesh.mesher.Mesher` carries them onto the mesh:
:class:`~orpheus.mesh.structured.Mesh1D` and
:class:`~orpheus.mesh.structured.Mesh2D` both carry ``face_roles``, a
:class:`~orpheus.mesh.face_roles.FaceRoles` mapping from face name to role
over exactly their boundary faces (``xmin`` and ``xmax`` on a slab or a hollow
body, ``xmax`` alone on a solid one, whose centre :math:`r = 0` is an interior
point and carries no role, and ``ymin``/``ymax`` added in 2-D;
:ref:`structured-geometry-face-laws`). A role is a
:class:`~orpheus.geometry.deck.MirrorDeck`, a
:class:`~orpheus.geometry.deck.TranslationDeck` or an ``int`` tag, and
``None`` is no role, refused on every face. A deck has no field: the face fixes
its mirror plane, its partner and its translation direction, so a deck
declared with the wrong axis is unspellable. The decks a face admits are those
of its chart (:func:`~orpheus.geometry.deck.admitted_decks`): a Cartesian axis
admits the mirror and the translation, and a **curved axis admits no deck**
(the curved-face rule below). The inner face of a hollow radial body is a tag
like any other, with one difference that decides its admission: the radial
axis has one trace slot, the outer surface, so the inner tag has no trace and
no operator to land in (:ref:`structured-geometry-hollow-inner-law`).

**Stage 2 --- the Problem binds and realizes.**
:class:`SNProblem` takes ``boundary_conditions``, a
``BoundaryConditions`` or a plain
mapping from tag to ``WallCondition``
(the one shared exterior type, an affine pair of a response network and an
inflow). Construction runs two steps:

1. ``bind(problem, boundary_conditions)`` joins the declared conditions to the
   problem's tagged faces and keeps the result as :attr:`SNProblem.binding`.
   The join is an exact cover: a tag no condition names raises
   ``UnboundTagError``, a condition no face carries raises
   ``UnusedConditionError``, and a condition that is not invariant under the
   problem's colouring symmetry raises ``NonInvariantConditionError``. The
   carrier is the problem itself, so the symmetry group the binding is checked
   against is the problem's own and never another carrier's.
2. :func:`~orpheus.transport.method.realize_boundary`, the one method-generic
   face loop, hands each face's declaration (its deck, or the condition its tag
   is bound to) to :meth:`SNProblem.realize_boundary_role
   <orpheus.sn.problem.SNProblem.realize_boundary_role>`. That method builds an
   ``SNMethodSpace`` (the quadrature and the trace ladder, two fields and no
   face) and calls ``SN_REALIZERS.realize``, which admits first and realizes
   second. A declaration the table does not admit on that face raises
   ``BoundaryRefused`` at construction, carrying the refusal's axis and the
   issue that tracks it; an admitted one returns a typed
   ``RealizedBoundary``. The results are stored face by face in
   :attr:`SNProblem.bc`, a ``dict[str, RealizedBoundary]`` over the faces that
   have a trace (the hollow inner tag has none, so it meets the table only
   through its admission).

``SN_REALIZERS`` (``orpheus/sn/boundary/realizer.py``) is a sealed
``RealizerTable`` and the ONE door for SN: the question "is this declaration
admitted on this face" and the act "realize it" are answered by the same
table, so the query and the realization cannot disagree (the string registry
that preceded it admitted two tags while its realizer realized six laws). The
table is keyed by the units of the response network, never by
``WallCondition`` itself, so a method never writes a branch over the shared
family's cone point; the realization is a functor on the wall's network
(:ref:`bc-sn-resolution-table`).

**The declaration is the only source of a boundary condition.**
No S\ :sub:`N` entry point carries a default: the eigenvalue entries
(:func:`~orpheus.sn.solver.solve_sn`, ``solve_sn_adjoint``) and the source
entries (:func:`~orpheus.sn.solver.solve_sn_fixed_source`,
``solve_sn_adjoint_fixed_source``, ``solve_sn_multiplying_source``) each take
a keyword ``boundary_conditions`` (default
``BoundaryConditions.empty()``, the declaration of a carrier whose faces are all
decks) and hand it to :class:`SNProblem`. A fixed-source problem on a
vacuum-bounded domain binds ``WallCondition.vacuum()`` to its tag; an
infinite-lattice eigenvalue problem declares mirror decks on its faces (or
binds ``WallCondition.specular(1.0)``, the same operator;
:ref:`sn-loss-kernel-gauge`). A tagged face with no condition is refused
loudly, so a face cannot be left undeclared
(:ref:`structured-geometry-no-default-law`).

.. code-block:: python

   from orpheus.derivations.common.xs_library import get_mixture
   from orpheus.geometry import CoordSystem, StructuredGeometry
   from orpheus.geometry.deck import MirrorDeck
   from orpheus.mesh import CellsByCount, Mesher
   from orpheus.numerics.quadrature import Quadrature
   from orpheus.sn.problem import SNProblem
   from orpheus.transport.boundary_conditions import WallCondition

   def slab(roles):
       geometry = StructuredGeometry(coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 2.0), mat_ids=(0,), face_roles=roles)
       return Mesher(geometry).partition(CellsByCount.uniform_volume(8)).mesh

   quadrature = Quadrature.gauss_legendre(n_ordinates=8)
   materials = {0: get_mixture("B", "2g")}

   # a mirror at x = 0 (a deck) and a vacuum tag (0) at x = 2
   problem = SNProblem(slab((MirrorDeck(), 0)), quadrature, materials, boundary_conditions={0: WallCondition.vacuum()})
   for face, op in problem.bc.items():
       print(face, type(op).__name__, op.return_fraction, op.is_vacuum, op.closes_axis)

   # a partial wall: 30 % specular plus 50 % diffuse return on both faces (one tag)
   wall = WallCondition.specular(0.3) + WallCondition.white(0.5)
   problem = SNProblem(slab((0, 0)), quadrature, materials, boundary_conditions={0: wall})
   for face, op in problem.bc.items():
       print(face, type(op).__name__, op.return_fraction, [type(t).__name__ for t in op.terms])

Output (``[M]`` 2026-10-24, ``/tmp`` script run with ``PYTHONPATH`` set to the
worktree, ``orpheus.__file__`` printed first):

.. code-block:: text

   xmin ScaledOrdinatePermutation 1.0 False True
   xmax ZeroMap 0.0 True False
   xmin RealizedSum 0.8 ['ScaledOrdinatePermutation', 'RankOneReemission']
   xmax RealizedSum 0.8 ['ScaledOrdinatePermutation', 'RankOneReemission']

The mirror deck realizes to a unit scaled permutation that closes its axis, the
vacuum tag to the zero map, and the partial wall to the sum of its two chains'
images with the return fractions adding (:math:`0.3 + 0.5 = 0.8`).

**The curved-face rule.** A reflecting curved wall is
``WallCondition.specular(1.0)`` bound to a tag, never a deck. The admissible
decks of a chart are :math:`N(G_c)/G_c`, the normaliser of the chart's own
continuous symmetry group modulo that group, and for the radial axis of a
cylinder or a sphere that quotient is trivial: the rotation group already acts
on the face and the mirror it would add is not a new symmetry of the domain.
The geometry refuses the declaration at construction (``[M]`` 2026-10-24, the
script of :ref:`bc-sn-resolution-table`):

.. code-block:: text

   DeckError - StructuredGeometry.face_roles['xmax']: a MirrorDeck on axis 0 of a cylindrical chart: W = N(G_c)/G_c is trivial on the

(truncated at 118 characters of the message; the message goes on to name the
response that replaces the deck). The specular wall the curved face takes is
admitted at the outer face, and is the only specular wall admitted at a hollow
body's inner face (below).

.. note::

   Before the wave-2 carve (2026-10-24) the declaration was made on the mesh
   as ``face_laws`` and the laws were string-tagged ``BC`` values or typed
   law objects, resolved through a class-level registry by
   ``realize_boundary_law``. That architecture is recorded, in the past tense
   and with each retired name, in :ref:`bc-sn-retired-architecture` at the end
   of this page.

Supported Types
---------------

**Reflective** (specular reflection).
At the outer boundary :math:`r = R` (or :math:`x = L`), the incoming
flux for :term:`ordinate` :math:`n` is set to the outgoing flux of its reflected
partner:

.. math::
   :label: reflective-bc

   \psi_n^{\rm in} = \psi_{n'}^{\rm out}

where :math:`n'` is the reflected partner ordinate (negating the
appropriate direction cosine).  The pairing is derived from the mirror
motion by each :term:`quadrature`'s
:meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`
(:eq:`quadrature-ordinate-permutation`) and **certified** by
``ordinate_pairing(quadrature, mirror_motion(axis, dimension))``
(``orpheus/numerics/quadrature/ordinate_pairing.py``), which folds onto that
permutation and asserts three independent claims on it: the trace measure
:math:`w\,|\mu_a|` is preserved to ``MEASURE_RTOL`` (ERR-042), the pairing is
an involution (ERR-044), and every non-tangential ordinate pairs with the
opposite sign class on the axis (ERR-045). Each passes on a table the other
two reject, so none is dropped for another, and the function returns the
permutation it certified, so the realizer has one computation and not a
certification beside a second derivation. A rule that is not closed under the
mirror (an odd-``n_phi`` product rule has no x-mirror closure) is refused on the
angular-resolution axis (:ref:`bc-refusal-axes`).

There are two spellings of the same wall, and one operator. On a Cartesian face
the mirror is a deck, ``MirrorDeck()`` (a symmetry plane of the domain: no
amplitude, nothing to scale); on any face, and the only way on a curved one, it
is ``WallCondition.specular(1.0)`` (a wall that returns everything). Both realize
to a ``ScaledOrdinatePermutation`` at weight 1.0, whose ``closes_axis`` is true,
and the two operators are bit-identical (``[M]`` 2026-10-24, ``np.array_equal`` of
the two realized ``inner`` operators applied to one random vector, an 8-ordinate
Gauss-Legendre slab: 2 of 2 faces). This is the law an infinite-lattice or infinite-medium eigenvalue problem declares on its
faces.
The CP solver uses white (isotropic) BCs instead; see
:ref:`white-bc-quality` for a comparison showing the ~1% gap between
the two approaches.

.. warning::

   ⚠ **Two reflective axis pairs plus diamond differencing make the
   within-group loss operator EXACTLY SINGULAR** (#344).  A specular
   face is a *closed* boundary for the diamond closure's undamped face
   sawtooth, so with :math:`\ge 2` axes reflective at both ends
   :math:`A = L+C-S-N_{2n}-B` (:eq:`sn-within-group-with-n2n`) has a
   non-trivial kernel — the kernel is a *pure-trace* object, so the two
   bulk gains :math:`S` and :math:`N_{2n}` annihilate it and neither
   enters the count — ``[M]``
   :math:`\dim\ker A = 12` on a 2-D ``level_symmetric`` :math:`S_4`
   2-group box, :math:`138` at :math:`d=3` :math:`(3,4,5)` — and the
   returned boundary trace is one member of a solution manifold.  This
   is the STANDARD infinite-lattice :math:`k_\infty` declaration (every
   face reflective) for :func:`~orpheus.sn.solver.solve_sn` in two or
   more dimensions.  Nothing a user normally checks reveals it — every
   mirror-even functional is blind by theorem — and the solver projects
   the trace onto the canonical member and says so.  The gauge reads the
   REALIZED operator (``closes_axis``), so a mirror deck and a unit
   specular wall are gauged alike.  Full treatment:
   :ref:`sn-loss-kernel-gauge` in :doc:`cartesian_multid`.

**Vacuum** (zero incoming flux).
All incoming angular fluxes at the face are set to zero:

.. math::
   :label: vacuum-bc

   \psi_n^{\rm in} = 0

.. (vv-status rationale) definition: Definitional / notation introduction. Operational rule ψ_n^in = 0 for vacuum boundary; semantics exercised by every vacuum-BC test (test_boundary_conditions, MMS suite); no isolated identity to verify.
.. vv-status: vacuum-bc documented


Vacuum is the zero network, ``WallCondition.vacuum()``, and realizes to the
``ZeroMap`` :math:`\Gamma_+ \to \Gamma_-`; the question "is this face a vacuum"
is answered by the *realized* map (``op.is_vacuum``), which is also the image of
a pure inflow (its source travels the boundary-source channel) and of any wall of
zero weight (``WallCondition.specular(0.0)`` is the vacuum: one cone point, one
value). The zero map reads its half-space off the trace ladder, so it fills only
the incoming trace. Vacuum BCs are the natural choice for fixed-source MMS
verification on finite slabs (see :ref:`sn-mms-verification`).

**White and partial walls** (diffuse return).
``WallCondition.white(d)`` is the Lambertian wall: it returns the fraction
:math:`d` of the outgoing partial current, re-emitted uniformly (isotropically in
the face's half space) on :math:`\Gamma_-`. It realizes to a
``RankOneReemission``, the factored rank-one kernel :math:`d\,(E \circ R_t)` with
:math:`R_t` the contraction of :math:`\Gamma_+` to its partial current and
:math:`E J = J/M_-` the uniform section, where :math:`M_- = \sum_{\Omega\cdot n<0}
w\,|\Omega\cdot n|` is the INFLOW mass. The divisor is the inflow mass because that
is what makes :math:`R_t^{\rm in}\circ E = \mathrm{id}`, hence a return fraction
that is exactly :math:`d` and not :math:`d\,M_-/M_+`; on every shipped rule
:math:`M_+ = M_-` bitwise (``[M]``
``tests/gates/sn/architecture/test_realized_boundary.py::test_the_mass_census_inflow_equals_outflow_on_every_registry_rule``,
70 of 70 (rule, face) pairs), so the choice is invisible on the shipped rules and
decisive on any rule that breaks the symmetry. ``WallCondition.specular(s)`` for
:math:`s < 1` is the partially specular wall, a ``ScaledOrdinatePermutation``
with :math:`\alpha = s`; a sum of the two, as in the code block above, is a
``RealizedSum``.

**Prescribed inflow.**
An inflow is a nonzero source :math:`q \ge 0` over a network (``WallCondition.uniform(v)``,
``WallCondition.framed(profile)``; both over the zero network, so the linear part is
the zero map). The linear factor realizes to the ``ZeroMap``; :math:`q` is not an operator
and travels the boundary-source channel (:ref:`bc-affine-source-channel`): the
``SNSource`` realizer samples it on :math:`\Gamma_-` of its face
(:meth:`SNProblem.declared_boundary_source
<orpheus.sn.problem.SNProblem.declared_boundary_source>` reads the problem's
binding through ``AngularBoundarySourceSink.from_method`` and packs each face's
sample). ``SNSource.sample`` is the one place SN checks, at the call, that the trace
it samples on is an inflow half-trace (ERR-047): a source written anywhere else is
discarded by the sweep. The two eigenvalue entries (``solve_sn``, ``solve_sn_adjoint``) refuse a
declared source (``_refuse_a_boundary_source``: a k-eigenvalue problem has no
external inflow), and ``solve_sn_adjoint_fixed_source`` builds a zero boundary
source and never reads the declared :math:`q`, which is correct: :math:`q`
enters the reciprocity identity as :math:`\langle q, |\Omega\cdot n|\gamma_-\psi^+\rangle`,
not through the operator.

**Periodic** (the translation wrap).
A ``TranslationDeck`` on two opposite faces is the identity relabelling between two
*distinct* index sets, :math:`\Gamma_+(f') \to \Gamma_-(f)` with :math:`f'` the
opposite face, realized (``SNDecks``) as the unit ``ScaledOrdinatePermutation`` whose
domain face is the partner. The crossing is earned, not assumed: a translation moves
no direction and the outward normals are opposite, so :math:`n_f = -n_{f'}` sends
outgoing to incoming, and the type's constructor refuses a pairing that does not
cross or that moves measure. A pair whose crossing sets differ is refused on the
topological axis (it is not a translation quotient and needs an explicit gluing
map, #178).

.. _bc-tensor-decompositions:

Boundary conditions as tensor decompositions
---------------------------------------------

The boundary conditions used by :class:`SNProblem` are concrete instances
of a more general tensor-decomposed framing, defined in
``orpheus.transport.boundary_conditions`` (the response network and its chains;
:ref:`bc-architecture-of-record`). A boundary condition is a linear
operator :math:`B` mapping the outgoing angular flux at a face to the
incoming angular flux:

.. math::
   :label: bc-tensor-decomposition

   \psi_{\rm in}(\Omega)
   = (B\,\psi_{\rm out})(\Omega)
   = \sum_\alpha \bigl(G_\alpha\,\psi_{\rm out}\bigr)(\Omega) \cdot A_\alpha,

.. (vv-status rationale) Definitional/literature-transcribed framing (a BC
   as R = sum_alpha G_alpha A_alpha, Lewis & Miller 1984 §3.4). The concrete
   rank-N primitives (vacuum/reflective/white) carry their own verification.
.. vv-status: bc-tensor-decomposition documented

where :math:`G_\alpha` is a **per-term operator** (permutation,
pushforward, angular average, spatial wrap) and :math:`A_\alpha` is a
**scalar amplitude** (typically an :term:`albedo` :math:`\in [0, 1]`).

This is the same algebra Lewis & Miller (1984) §3.4 use to introduce
boundary conditions in transport: every BC of practical interest is
either rank-1 (one :math:`G_\alpha \otimes A_\alpha` term) or a finite
linear combination of rank-1 primitives (rank-N). In the shared family
the terms are the two **chains** of the response network on the rungs of the
face's trace ladder (outflow half-trace :math:`\Gamma_+`, partial current
:math:`S`, inflow half-trace :math:`\Gamma_-`): the specular chain :math:`P`
(a mirror-plane ordinate permutation) and the diffuse chain :math:`E\circ R_t`,
so that

.. math::
   R \;=\; s\,P \;+\; d\,(E\circ R_t),
   \qquad s,\,d \ge 0,\; s + d \le 1,

a point of the convex cone of walls, and the wall's return fraction is
:math:`\rho = s + d` (exact by the split pair :math:`R_t^{\rm in}\circ E = \mathrm{id}`
and the measure preservation of :math:`P`; the derivation, the cone, and why a
third chain is an edit of a closed sum are
:ref:`bc-architecture-of-record`). The :math:`G_\alpha` of this section is the
whole per-term map, which is the composite response-after-geometry, never the
response alone; the letters :math:`G`, :math:`R` of the affine trace form
:eq:`affine-bc-form` name the two *factors* of one term
(:ref:`bc-factor-roles`).

The tensor framing pays off because partial-current Marshak boundaries
(:math:`B = c_1 \, B_{\rm refl} + c_2 \, B_{\rm diff}` — a specular term plus a
Lambertian term, Bell & Glasstone 1970 §1.5) and multi-region interface
couplings are all instances of the same algebra: pick the per-term maps, pick the
amplitudes, sum. Here the Marshak mix is the sum of two walls,
``WallCondition.specular(c1) + WallCondition.white(c2)``, whose arithmetic is the
cone's own: a sum or a scalar multiple is a ``WallCondition`` exactly when its
coordinates stay in the cone, and is otherwise a displacement that is never
bindable. A mirror *deck* is not an operand of the arithmetic: it has no
amplitude, and ``MirrorDeck() * 0.5`` and ``MirrorDeck() + MirrorDeck()`` are
``TypeError`` (``[M]`` 2026-10-24), the language's own refusal in place of a
guard. A new condition of the shared family is a new chain in the closed sum, not a
registry entry; the SN consequence is that every ``match`` over the chains
(``SN_REALIZERS``' registrations) reddens until it handles the new member.

.. _bc-albedo-reemission-closure:

Why a partial wall is two chains, and why a bare albedo has no spelling
------------------------------------------------------------------------

A surface's re-emission law is a product of two **independent** degrees
of freedom: *how much* of the arriving flux returns (:math:`\alpha`)
and *in what angular shape*. Every shape admits every amplitude, so the
shared type carries them as the two coordinates :math:`(s, d)` of the cone
rather than as a menu of their product:

.. code-block:: python

   from orpheus.transport.boundary_conditions import WallCondition

   polished = WallCondition.specular(0.7)          # 0.7 of the outflow, mirrored
   matte = WallCondition.white(0.7)                # 0.7 of the outflow, Lambertian
   both = WallCondition.from_amplitudes(0.3, 0.4)  # the cone point (0.3, 0.4)
   print(polished.cone_point, matte.cone_point, both.cone_point, both.return_fraction)

which prints (``[M]`` 2026-10-24):

.. code-block:: text

   (0.7, 0.0) (0.0, 0.7) (0.3, 0.4) 0.7

There is **no spelling of an amplitude without a shape**: the facades name the
shape, and ``from_amplitudes`` takes both. This is the SN-side consequence of a
fact about the trace. A **scalar** trace has one boundary degree of freedom,
:math:`J^- = \rho J^+`; there is no angular distribution to fix, so the return
fraction is the whole law, and a scalar method (diffusion) realizes a chain of
weight :math:`w` as the albedo :math:`w` whichever the chain is. An **angular** trace resolves the full hemisphere: there
:math:`R = \alpha\,I` is a map :math:`\Gamma_+ \to \Gamma_+` and supplies no
crossing to :math:`\Gamma_-`, so nothing says which outgoing direction feeds which
incoming one, and composing it anyway pairs incoming ordinate :math:`j` with
outgoing ordinate :math:`j`, an artefact of **array position** and not of geometry
(:ref:`the angular-resolution axis <bc-method-realizability>`). The shared type
makes that composition unspellable, so SN has no refusal to raise and the
under-determined law has no value.

.. dropdown:: Why the positional pairing is worse than meaningless (a first attempt, measured)

   A bare albedo law existed until wave 2, and SN refused it at the realizer
   because composing it positionally is a **configuration-dependent accident**.
   Measured `[M]` on the ``xmax`` face while that law existed, comparing the
   positional pairing's index map against the specular one:

   .. list-table::
      :header-rows: 1
      :widths: 40 30

      * - quadrature
        - positional :math:`=` specular?
      * - ``product(2, 4)``
        - **True**
      * - ``level_symmetric(6)``
        - **True**
      * - ``gauss_legendre(4)`` / ``(8)``
        - False (the slab mirror reverses order)
      * - ``lebedev(17)``
        - False

   So before the refusal, a bare albedo behaved *exactly like a mirror* on two of
   the tree's quadratures and like nothing in particular on the others,
   silently and with no seam to notice. A user who validated on
   ``level_symmetric`` and ran production on ``lebedev`` got different physics
   from the same law object. The old answer was not a defensible default that a
   closure would merely make explicit; it was a coincidence of index order that
   *looked* defensible on the half of the fixture set where the coincidence held.
   That is the argument for the design that replaced it, in which the shape is a
   chain and the pairing is *certified*, never positional.

**Where the pairing's invariants live.** The specular pairing :math:`\pi` carries
three independent invariants (measure preservation, ERR-042; involution, ERR-044;
inflow to outflow, ERR-045). They are certified in ONE place, ``ordinate_pairing``
(above), and the realized type checks the two it needs when it is constructed:
``ScaledOrdinatePermutation`` refuses at construction a pairing that does not cross
(ERR-045) and one that does not carry the outflow trace measure onto the inflow
trace measure, which is what makes ``return_fraction == alpha``. Before
that check a mirror composed with a shift constructed, reported a return fraction of
1.0 and measured 1.198 on a non-uniform outflow (the record in the
``realized_operators`` module docstring); the involution is not checked by the
type because the gather uses the inverse and the facts do not depend on it. Both
carriers of a specular return (the deck and the constitutive wall) therefore fire
the same certification, because both are built by the same realizer.

.. _bc-sn-resolution-table:

SN BC resolution table
----------------------

The realizer table of SN, ``SN_REALIZERS``, answers admission and realization
together. It is keyed by the *units* of the declared network (a chain of the closed
sum ``Chain``, the marker ``ZeroNetwork``, the marker ``InflowTerms``) and by the
decks, and it is filled by five realizer classes that register by class keyword:

.. list-table:: Declared unit → SN realizer → realized SN operator
   :header-rows: 1
   :widths: 24 16 36 24

   * - declared unit
     - realizer
     - realized operator :math:`\Gamma_+ \to \Gamma_-`
     - return fraction
   * - ``ZeroNetwork`` (``WallCondition.vacuum()``, a zero weight, the linear
       part of an inflow)
     - ``SNZero``
     - ``ZeroMap``: the zero morphism between the two half-traces, whose
       space hooks emit :math:`|\Gamma_-|` rows forward and :math:`|\Gamma_+|`
       rows on the transpose
     - 0
   * - ``s * SpecularChain`` (``WallCondition.specular(s)``)
     - ``SNSpecularChain``
     - ``ScaledOrdinatePermutation`` :math:`s\,P_\pi`, with :math:`\pi` the
       certified mirror pairing of the face's own plane
     - :math:`s`
   * - ``d * DiffuseChain`` (``WallCondition.white(d)``)
     - ``SNDiffuseChain``
     - ``RankOneReemission`` :math:`d\,(E\circ R_t)`, the factored Lambertian:
       ``PartialCurrentOperator`` contracts :math:`\Gamma_+` to the partial
       current and ``IsotropicEmissionOperator`` re-emits it on :math:`\Gamma_-`
     - :math:`d`
   * - a wall of both chains
     - the meet of the two admissions
     - ``RealizedSum`` of the two images, in the canonical chain order
     - :math:`s + d`
   * - ``MirrorDeck()``
     - ``SNDecks``
     - ``ScaledOrdinatePermutation`` at :math:`\alpha = 1`, the same matrix as
       ``specular(1.0)`` (a symmetry has no amplitude)
     - 1
   * - ``TranslationDeck()``
     - ``SNDecks``
     - ``ScaledOrdinatePermutation`` at :math:`\alpha = 1` with the identity
       permutation, bound :math:`\Gamma_+(f') \to \Gamma_-(f)` (the PARTNER
       face's outflow)
     - 1
   * - ``InflowTerms`` (the source :math:`q`)
     - ``SNSource``
     - no operator: ``sample`` returns an array on :math:`\Gamma_-` of the face
     - n/a

The realizer is a functor on the response network, so the image of a sum is the sum
of the images. The facts the consumers read of a realized face (the return fraction,
whether it relabels ordinates, whether it closes its axis into a mirror loop, whose
outflow it consumes, whether it is the zero map) are properties of the *type*,
derived from its parameters and never stored beside them: the leakage functional
:eq:`sn-leakage-functional` skips a face exactly when ``op.return_fraction == 1.0``,
the gauge reads ``closes_axis``, and diffusion-synthetic acceleration reads
``is_vacuum or closes_axis``. A lying ``alpha`` is unspellable because the operator
is built from spaces read off the trace ladder.

There is deliberately **no sweep-cycle column**. Whether a
configuration sweeps in one pass is a property of the whole face
configuration, not of the law on any single face — a reflecting
face opposite a vacuum is acyclic, two reflecting faces are not —
so it cannot be tabulated per row. The per-law
``creates_sweep_cycle`` ``ClassVar`` that once occupied such a column
was retired 2026-07-30 for exactly that reason; the configuration-
level criterion, and the record of why the flag could not work, are
at :ref:`bc-sweep-cycle`.

**Admission, by site.** A face reads as one of three site classes: ``CART`` (any
Cartesian face), ``RAD_OUTER`` (the outer face of a solid or hollow radial body)
and ``RAD_INNER`` (a hollow radial body's inner face, a tag with no trace). The SN
columns, typed by hand in
``tests/gates/transport/test_admission_matrix.py`` (``_SN_COLUMNS``) and pinned
against ``SN_REALIZERS`` exactly, are:

.. list-table:: SN admission (R = realized; refusals name axis and issue)
   :header-rows: 1
   :widths: 34 16 24 26

   * - declared value
     - ``CART``
     - ``RAD_OUTER``
     - ``RAD_INNER``
   * - vacuum
     - R
     - R
     - PLUMBING #511
   * - specular 1.0
     - R
     - R
     - R
   * - specular 0.5
     - R
     - R
     - PLUMBING #511
   * - white 0.5
     - R
     - PLUMBING #653
     - PLUMBING #511
   * - both chains (0.3, 0.5)
     - R
     - PLUMBING #653
     - PLUMBING #511
   * - inflow 1.0
     - R
     - PLUMBING #653
     - PLUMBING #511
   * - specular 0.5 over inflow 1.0
     - R
     - PLUMBING #653
     - PLUMBING #511

Why each refusal is there:

* **#653, the radial corner action.** The corner action of the radial
  characteristic is ruled only for a permuting return or a zero return, so the
  diffuse chain and every source are refused at ``RAD_OUTER`` *at admission*
  (they used to construct and refuse at the first sweep: 12 of 12 declarations in
  plan section 4.3). A mixed wall at a radial outer site is refused through its
  diffuse chain, because a wall of two chains is the meet of the two chain
  admissions.
* **#511, the hollow inner face.** The radial axis has one trace slot, the outer
  surface, so a condition declared at the inner face has no operator to land in,
  and a method built on the axes computes a reflective cavity silently. SN
  therefore admits the inner site only at a unit specular return, the cavity its
  angular closure computes, and refuses the rest. ``[M]`` before the refusal an
  inner vacuum and an inner reflective law gave bit-identical SN fluxes on a
  cylinder and a sphere with :math:`r_0 = 0.5` (2026-09-25); the gates are
  ``tests/gates/mesh/test_hollow_inner_law.py`` (``test_sn_refuses_a_declared_inner_law_it_would_drop``,
  ``test_a_direct_mesh_is_refused_too``, and the last gate, which verifies that the
  admitted ``specular(1.0)`` is the void-cavity answer). Diffusion's table refuses
  the same declarations with the same shared refusal.
* **The angular-resolution axis.** With a method space, the specular chain and the
  mirror deck ask ``ordinate_pairing``; a rule with no closure under the mirror is
  ``Refused(ANGULAR_RESOLUTION)``. ``[M]`` 2026-10-24 (the admission of
  ``WallCondition.specular(0.5)`` on ``xmin`` of a 2-D face layout, product rules
  with ``n_mu = 2``): ``n_phi = 3`` is ``Refused`` on ``ANGULAR_RESOLUTION``,
  ``n_phi = 4`` is ``Realized``, and the diffuse chain, which asks no pairing, is
  ``Realized`` on both. The gate is
  ``tests/gates/sn/architecture/test_sn_realizer_table.py::test_a_certification_failure_names_its_axis_at_the_table``.
* **The unregistered-unit door.** SN registers every shared unit, so its table
  meets only another method's own condition there: diffusion's zero flux is
  ``Refused(STATE_CONE)`` because the table's cone is ``NONNEGATIVE_FLUX`` and that
  unit's return fraction is negative. SN reaches this verdict without importing the
  diffusion type (an L3 to L3 edge the layer gate forbids).

The admission is the table's; whether the problem constructs is therefore the
same question as whether the table admits, and it is read at the binding of the
carrier's faces, so a mesh built without a geometry (the bare constructor) is
covered too.

**The method space.** :class:`SNProblem` constructs the table's argument as
``SNMethodSpace(quadrature=..., trace=...)`` from the quadrature and the unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` that the
construction body builds once for every supported mesh (1-D Cartesian, spherical,
cylindrical, 2-D Cartesian). The space names no face: the face a realization acts
on is the ``FaceSite`` the realizer receives beside it, so there is exactly one face
value per realization (two independently supplied encodings of one orientation are
ERR-041's defect class, and the invariant itself, that the inflow half-space holds
exactly the ordinates with :math:`\Omega\cdot n_f < 0`, is witnessed on the
trace ladder in ``tests/gates/numerics/test_angular_trace_space.py``). The reflective
arm reads its reflection axis from the face's own name, ``AXIS_NAMES[axis_index]``,
so the partner is correct at any dimension by construction
(:ref:`bc-face-name-latent-d3-bug`).

.. note::

   **A realized SN face is typed** :math:`\Gamma_+ \to \Gamma_-`: it consumes the
   outflow half-trace and produces the inflow half-trace, and the consumer composes
   :math:`\iota_-\circ\text{law}\circ\gamma_+`. A shape assertion cannot tell the
   two typings apart (:math:`|\Gamma_+| = |\Gamma_-|` on every quadrature × face in
   the tree), so read the *declared spaces*, never the output shape. Full
   derivation at :ref:`bc-domain-narrowing`.

A slab :class:`SNProblem` has no y-axis in its
:attr:`~orpheus.sn.problem.SNProblem.axes` tuple, so
:func:`~orpheus.mesh.axis.face_labels` emits no y-label and :attr:`SNProblem.bc`
has no y-entry: ``slab.bc["ymin"]`` is a :class:`KeyError`, not a no-op
(:ref:`bc-face-name-carve-what-retired`).

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
condition at :math:`r = 0` --- the geometry handles it naturally, and a solid
radial body carries no face role there (:ref:`structured-geometry-face-laws`).

The **outer face** of a curvilinear body admits the vacuum and any specular wall,
``WallCondition.specular(alpha)``, and refuses the diffuse chain and every inflow
(#653; the table above). A reflecting curved wall is ``specular(1.0)``, never a
mirror deck (the curved-face rule). The **inner face of a hollow body**
(:math:`r_0 > 0`) is a tag with no trace and admits only ``specular(1.0)`` (#511).
Both are enforced by the table at construction. ``[M]`` 2026-10-24, a script of the
three declarations over the cylinder (a solid body ``(0, 1)`` and a hollow one
``(0.5, 1)``, ``folded_product(4, 8)``, 8 cells, 2-group mixture B):

.. code-block:: python

   from orpheus.derivations.common.xs_library import get_mixture
   from orpheus.geometry import CoordSystem, StructuredGeometry
   from orpheus.geometry.deck import MirrorDeck
   from orpheus.mesh import CellsByCount, Mesher
   from orpheus.numerics.quadrature import Quadrature
   from orpheus.sn.problem import SNProblem
   from orpheus.transport.boundary_conditions import BoundaryRefused, WallCondition

   def cylinder(breakpoints, roles):
       geometry = StructuredGeometry(coord=CoordSystem.CYLINDRICAL, breakpoints=breakpoints, mat_ids=(0,), face_roles=roles)
       return Mesher(geometry).partition(CellsByCount.uniform_volume(8)).mesh

   quadrature = Quadrature.folded_product(n_mu=4, n_phi=8)
   materials = {0: get_mixture("B", "2g")}

   # a curved face has no mirror deck ...
   try:
       cylinder((0.0, 1.0), (MirrorDeck(),))
   except ValueError as error:
       print(type(error).__name__, "-", str(error)[:118])

   # ... its reflecting wall is a response, admitted at the outer face with the vacuum;
   # the diffuse return is refused there (#653)
   for name, wall in (("vacuum", WallCondition.vacuum()), ("specular(1.0)", WallCondition.specular(1.0)), ("white(0.5)", WallCondition.white(0.5))):
       try:
           problem = SNProblem(cylinder((0.0, 1.0), (0,)), quadrature, materials, boundary_conditions={0: wall})
           print("outer", name, "->", type(problem.bc["xmax"]).__name__)
       except BoundaryRefused as refusal:
           print("outer", name, "-> refused:", refusal.admission.axis.name, refusal.admission.issue)

   # the hollow body's inner face is a tag with no trace: only the unit specular cavity is admitted (#511)
   for name, wall in (("vacuum", WallCondition.vacuum()), ("specular(1.0)", WallCondition.specular(1.0))):
       try:
           problem = SNProblem(cylinder((0.5, 1.0), (1, 0)), quadrature, materials, boundary_conditions={1: wall, 0: WallCondition.vacuum()})
           print("inner", name, "-> admitted; faces with a trace:", sorted(problem.bc))
       except BoundaryRefused as refusal:
           print("inner", name, "-> refused:", refusal.admission.axis.name, refusal.admission.issue, "on", refusal.face)

.. code-block:: text

   DeckError - StructuredGeometry.face_roles['xmax']: a MirrorDeck on axis 0 of a cylindrical chart: W = N(G_c)/G_c is trivial on the
   outer vacuum -> ZeroMap
   outer specular(1.0) -> ScaledOrdinatePermutation
   outer white(0.5) -> refused: PLUMBING 653
   inner vacuum -> refused: PLUMBING 511 on xmin
   inner specular(1.0) -> admitted; faces with a trace: ['xmax']

.. _bc-sn-retired-architecture:

The retired trace-law architecture (history)
--------------------------------------------

This section records the trace-law architecture that wave 2 retired (2026-10-24);
the present architecture is :ref:`bc-architecture-of-record`, and each retired name
maps to its successor in :ref:`bc-wave2-crosswalk`. Everything below is past tense
and is kept because the reasons it records are why the present design is as it is.

**The two-stage declaration, as it was.** The laws were declared on the mesh as
``face_laws`` (a ``FaceLaws`` value keyed by face name), each law either a ``BC`` tag
(a frozen dataclass of a ``kind`` string and a ``FrozenMapping`` of numeric
parameters, with the convenience instances ``BC.vacuum``, ``BC.reflective`` and
``BC.white``) or an already-typed ``BoundaryTraceLaw`` from the package
``orpheus.geometry.boundary``. ``SNProblem`` owned a class-level
``BOUNDARY_OPERATOR_REGISTRY`` mapping kind strings to law classes
(``"vacuum"`` to ``VacuumInflow``, ``"reflective"`` to ``ReflectiveBoundary``), and
``SNProblem.realize_boundary_law`` dispatched through ``SNBoundaryRealizer``
uniformly for every supported mesh, after Issues #188 and #176 (2026-05-11)
removed the curvilinear bypass: before them a curvilinear ``Mesh1D`` skipped the
realizer because the trace factory raised ``NotImplementedError`` on those
coordinate systems, and a ``_BoundBoundaryOperator`` shim carried a dual mode that
forwarded a legacy 2-arg ``apply(psi, quadrature)`` to the law's own body
(:ref:`bc-curvilinear-realizer-unification`). The registry was met in
``SNProblem.__init__`` by looking each face's tag up and raising ``ValueError``
listing the supported kinds, and a curvilinear body accepted only ``"reflective"``
and ``"vacuum"``. The shared ``resolve_boundary_conditions`` carried the face loop,
the tag-to-law parse and, until 2026-09-30, a reflective default for an undeclared
face (:ref:`bc-dual-registry`).

**The shim and the descriptor model.** The 2-arg ``bc.apply(...)`` interface was
retired entirely after Issue #186 (2026-05-11); ``BoundaryTraceLaw`` became a pure
descriptor with no callable interface, composed by the descriptor-tree algebra
``LawSum`` / ``LawScaled`` and lowered by ``realize_recursively``
(:ref:`bc-trace-law-descriptor-model`). The resolved boundaries at
``problem.bc["xmin"]`` were wrapped in ``_BoundBoundaryOperator``, a strict 1-arg
passthrough that carried the originating law (``bc[face].law.geometry_map``,
``.response_kernel``) and a ``kind`` tag reading the law's registry key, so that
``problem.bc["xmin"] == "vacuum"`` was true iff the law was ``VacuumInflow`` (the
string arm of ``__eq__``, #553, which the successor does not carry). The consumers
asked the law what the operator did (``law_permutes_ordinates``,
``response_kernel.is_zero``, ``.amplitude``); the typed realized operators answer
the same questions from their own type.

**The law primitives and the resolution table.** The table of that page mapped each
kind to a law class and a realized operator: ``"vacuum"`` to ``VacuumInflow`` and a
``ZeroOperator`` narrowed to :math:`\Gamma_+ \to \Gamma_-`; ``"reflective"`` to
``ReflectiveBoundary`` and ``PermutationOperator(local_perm) & IdentityOperator()``;
``"white"`` to ``WhiteBoundary`` and the factored Lambertian
(``IsotropicEmissionOperator @ PartialCurrentOperator``, with the law's declared
``axis`` and ``outward_sign`` cross-checked against the installation face before
construction); ``"periodic"`` to ``PeriodicBoundary`` and the partner-face
permutation (``PairedDeck.domain_face`` named the partner); ``"albedo"`` to
``AlbedoBoundary(alpha, closure)`` with the closure ``SpecularReturn`` or
``IsotropicReturn`` selecting the realized body, and the closure-free spelling
*refused* on an angular trace; and ``"prescribed_inflow"`` to ``PrescribedInflow``,
whose linear factor collapsed onto the zero morphism while its source travelled the
boundary-source channel (assembled by ``AngularBoundarySourceSink.from_problem_laws``).
The narrowing to :math:`\Gamma_+ \to \Gamma_-` landed in campaign phases B3.2
(``vacuum``, ``reflective``), B3.4a (``white``, ``prescribed_inflow``), B3.4b
(``albedo``, with the closure) and B3.4c (``periodic``'s partner-face channel; its
arrow derived at G6.3 step 7), and the strict xfails that pinned the deferred rows
retired with them (:ref:`bc-domain-narrowing`). The reasoning of the closure and of the
refusal survives in :ref:`bc-albedo-reemission-closure`.

**The albedo closure and the pairing invariants, as they were.** ``AlbedoBoundary``
carried the amplitude and the shape as two parameters (``AlbedoBoundary(0.7,
SpecularReturn(axis="x"))`` polished, ``IsotropicReturn(axis="x", outward_sign=+1)``
matte, ``AlbedoBoundary(0.7)`` unstated), the closure being amplitude-free because
the response-kernel protocol's ``amplitude`` member already carried :math:`\alpha`
for the diffusion realizer. The three pairing invariants (ERR-042, 044, 045) were
methods on ``ReflectiveBoundary`` while that was the only law standing on the
pairing, and moved to ``orpheus.geometry.boundary._specular`` once a specular
closure existed on albedo, so that both carriers fired the same certification;
wave 2 moved them again, to ``ordinate_pairing``. The pre-Wave-7 names
``VacuumBoundaryOperator``, ``SpecularBoundaryOperator``, ``WhiteBoundaryOperator``,
``PeriodicBoundaryOperator`` and ``AlbedoBoundaryOperator`` were retired in Wave O
step O.4a.1, and ``MixedBoundaryOperator`` in Wave 11 (rank-N boundaries became the
descriptor-tree algebra). The Marshak mix was
``0.3 * AlbedoBoundary(1.0, SpecularReturn("x")) + 0.7 * WhiteBoundary(axis="x",
outward_sign=+1, albedo=1.0)``, a ``LawSum`` of ``LawScaled``, lowered by
``realize_recursively(tree, method_space, SNBoundaryRealizer())``; the mirror
``ReflectiveBoundary`` was refused as an operand of the composition nodes
(``bc-deck-law-not-an-operand``), which the language's ``TypeError`` on a deck now
does. B3.0 had corrected a misfiling in the page's table: the cosine-weighted
hemispheric average (``LambertianReemission``) was filed as a geometric operator, and an average is
not multiplicative and not a bijection, so it is a response kernel and not a deck.

**The 1-D y-face placeholders.** Until C4 (#220) a slab ``SNProblem`` carried a pair
of realized no-op ``ReflectiveBoundary(axis="y")`` operators at ``bc_ymin`` and
``bc_ymax`` so that cross-dimensional code could read them without coordinate-system
gating; no production code ever read them, and C4 made them unrepresentable
(:ref:`bc-face-name-carve-what-retired`).

**Before this infrastructure.** The SN solver hardcoded reflective boundaries on
all faces and the then-production ``transport_sweep`` accepted a
``boundary_condition: str`` parameter, which was removed. The source entries later
carried a ``boundary_condition`` parameter of their own (default ``"vacuum"``),
which filled every face of an all-undeclared ``Mesh2D`` or axis tuple through a
helper, ``_apply_default_bcs``, while the resolver read any face left undeclared as
reflective; P1 step 3c of #405 (2026-09-30) retired the parameter, the helper and
the resolver's default together with ``None`` as a declaration. Until #412
(2026-09-18) the object that resolved the laws was called the *solver's augmented
mesh* (the class was ``SNMesh``); it is the Problem.

.. list-table:: Boundary-layer milestones of the S\ :sub:`N` page
   :header-rows: 1
   :widths: 16 62 22

   * - date
     - decision
     - issue / commit
   * - 2026-05-11
     - The 2-arg ``apply`` and the curvilinear bypass retired: every mesh routed
       through the realizer.
     - #186, #188, #176
   * - 2026-09-18
     - The solver's augmented mesh renamed the Problem (``SNMesh`` to ``SNProblem``).
     - #412
   * - 2026-09-30
     - The entries took the declared laws verbatim; the reflective default and
       ``None`` as a declaration retired.
     - #405, ``4318adc5``
   * - 2026-10-24
     - Wave-2 carve: the trace-law architecture replaced by face roles, one
       ``WallCondition`` and ``SN_REALIZERS``; WF2.4 retirement tip.
     - ``471c4737``
