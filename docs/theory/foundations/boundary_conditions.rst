.. _theory-boundary-conditions:

===============================================================
Boundary Conditions — Face Roles, Wall Conditions and Admission
===============================================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before touching** the face roles of a geometry or a mesh (:mod:`orpheus.geometry.deck`), the conditions package
(:mod:`orpheus.transport.boundary_conditions`) or any method's realizer table. The sections from
:ref:`bc-architecture-of-record` to :ref:`bc-wave2-crosswalk` state the architecture as it stands. The sections after the
crosswalk record the trace-law architecture that wave 2 of the layer DAG replaced (2026-10-24): their mathematics still holds
where they re-spell it, and a name in double backticks there is retired.

- **A boundary face closes the transport equation in exactly one of two ways, and each face carries one role saying which.**
  Either the domain *continues* across the face (a deck: the inflow trace is the outflow trace of a partner face carried by an
  isometry) or the domain *ends* there (the inflow trace is an affine function of the face's own outflow trace, supplied by an
  exterior model bound to an integer tag). The role is the sum ``FaceRole = DeckTransformation | int``, parsed by
  :func:`~orpheus.geometry.deck.parse_face_role` (a ``bool`` is not a tag) and stored once per face by
  :class:`~orpheus.mesh.face_roles.FaceRoles`. Gates: ``tests/gates/geometry/test_deck.py``,
  ``tests/gates/mesh/test_face_roles.py``. Derivation: :ref:`bc-aor-trace-closures`.
- **A deck has no amplitude, no source and nothing to scale.** The decks are a closed sum of two types with no field,
  :class:`~orpheus.geometry.deck.MirrorDeck` and :class:`~orpheus.geometry.deck.TranslationDeck`; each is the Koopman operator
  of an isometry, invertible and measure preserving, so its return fraction is exactly one. The admissible decks of a chart are
  :math:`W = N(G_c)/G_c`, which is trivial on the radial axis of a cylinder or a sphere: no 1-D curved face carries a deck, and a
  reflecting curved wall is a condition, ``WallCondition.specular(1.0)``. Gates: ``test_a_curved_face_admits_no_deck`` and
  ``test_every_deck_match_ends_in_assert_never`` in ``tests/gates/geometry/test_deck.py``. Derivation: :ref:`bc-aor-deck-family`.
- **A wall is one affine trace map on the sub-Markov simplex.** :class:`~orpheus.transport.boundary_conditions.WallCondition` is the pair
  :math:`(R, q)` with :math:`\gamma_-\psi|_f = R\,\gamma_+\psi|_f + q`, the linear part
  :math:`R = s\,P + d\,(E\circ R_t)` a point :math:`(s, d)` of the sub-Markov simplex :math:`K = \{s, d \ge 0,\ s + d \le 1\}` (compact, so not a cone: a positive multiple of a point of :math:`K` can leave it) and :math:`q \ge 0` (the source cone).
  Vacuum, the specular wall, the white wall and a prescribed inflow are facade constructors over this one type, so
  ``WallCondition.specular(0.0) == WallCondition.white(0.0) == WallCondition.vacuum()``, one value with one digest. Gates:
  ``tests/gates/transport/test_wall_condition.py``, ``tests/gates/transport/test_response_network.py``. Derivation:
  :ref:`bc-aor-wall-condition`.
- **Whether an operator on the face trace is a deck or a response is decided by two tests** (:ref:`bc-factor-roles`).
  *Multiplicativity* is necessary: the Koopman operator of a measure-preserving bijection satisfies
  :math:`G(\psi\varphi) = (G\psi)(G\varphi)` and an average never does, so the Lambertian average is a response by this test
  alone. A specular *wall* passes it (a permutation is multiplicative), and the second test, *is it a quotient of the physical
  domain*, decides: a wall is not a quotient, so its mirror is constitutive. It is the chain
  :class:`~orpheus.transport.boundary_conditions.SpecularChain` of a wall condition and never a deck. Gates:
  ``test_a_condition_is_not_a_deck_slot`` and ``test_a_deck_is_not_an_operand`` in
  ``tests/gates/transport/test_wall_condition_algebra.py``.
- **The response carries the crossing** :math:`\Gamma_+ \to \Gamma_-`. For a deck the crossing is the isometry's. For a wall no
  ambient isometry exists to provide it, so the response performs it: the permutation :math:`P` of the specular chain, or the
  retraction :math:`R_t` followed by the section :math:`E` of the diffuse chain. Every realized SN response is typed
  :math:`\Gamma_+ \to \Gamma_-` of its own face (:class:`~orpheus.sn.boundary.realized_operators.ZeroMap`,
  :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`,
  :class:`~orpheus.sn.boundary.realized_operators.RankOneReemission`). Gate:
  ``TestVacuumRealizesOnTheTrace.test_the_zero_map_runs_outflow_to_inflow_of_its_own_face`` in
  ``tests/gates/sn/operators/test_sn_boundary_realizer.py`` (``ERR-041``).
- **The return fraction of a wall is** :math:`\rho = s + d` **exactly, by the split pair, and the trace measure is a left
  eigenvector of the realized map**: :math:`R^{\mathsf T}(w|\mu|)_- = \rho\,(w|\mu|)_+`. The proof needs :math:`E` to divide by
  the *inflow* mass, so that :math:`R_t^{\rm in}\circ E = \mathrm{id}`. Gates: ``test_the_left_eigenvector_rider_through_the_realizer``
  and ``test_return_fraction_is_the_declared_and_the_measured`` in ``tests/gates/sn/architecture/test_sn_realizer_table.py`` pin the
  identity on the registry rules, where the two masses agree bitwise, so they are blind to which mass ``E`` divides by. The divisor
  has one witness, ``test_the_rider_on_a_trace_whose_masses_differ`` in ``tests/gates/sn/architecture/test_realized_boundary.py``
  (`[M]` mutation, ``RankOneReemission._derive_inner`` dividing by the outflow mass: 1 failed, 17 passed, and the one red is that
  test; the mass census of the same file is green under it). Derivation: :ref:`bc-aor-wall-condition`.
- **Membership of :math:`K \times \{q \ge 0\}` is a value predicate, so a wall outside it does not exist.** Direct construction outside :math:`K`
  (``WallCondition.specular(1.2)``, ``WallCondition.white(-0.1)``) raises
  :class:`~orpheus.transport.boundary_conditions.SubmarkovViolationError` or
  :class:`~orpheus.transport.boundary_conditions.BoundaryResponseNotPositiveError` (``ERR-046``, ``ERR-043``). The arithmetic
  ``+``, ``-``, unary ``-`` and scalar ``*`` returns a ``WallCondition`` if and only if the result lies in :math:`K`, and
  otherwise a :class:`~orpheus.transport.boundary_conditions.WallDisplacement`: a sensitivity direction that is never
  bindable. Gate: ``test_the_algebra_and_the_cone`` in ``tests/gates/transport/test_wall_condition_algebra.py``.
- **A binding joins the declared conditions to a carrier's tags, and is lawful only if it is an exact cover that the colouring
  group leaves invariant.** :func:`~orpheus.transport.boundary_conditions.bind` takes the carrier (which owns the faces and
  their colouring symmetry) and the conditions
  (:class:`~orpheus.transport.boundary_conditions.BoundaryConditions`, the exact twin of ``Materials``), refuses an unbound tag
  and an unused condition, and refuses a condition that a non-identity element of the problem's colouring group moves. A 2-D
  mesh has no group machinery yet, and refuses a non-face-local condition on a tag of several faces (#658). Gate:
  ``tests/gates/transport/test_binding.py``. Derivation: :ref:`bc-aor-roles-binding`, :ref:`bc-aor-colouring`.
- **Realization is a partial map and admission is its domain, and both are read from ONE table per method.**
  :class:`~orpheus.transport.boundary_conditions.RealizerTable` answers :meth:`~orpheus.transport.boundary_conditions.RealizerTable.admit`
  (``Realized(faithful)``, ``Approximated`` or ``Refused(axis, reason, issue)``) and
  :meth:`~orpheus.transport.boundary_conditions.RealizerTable.realize` through one door, so ``admit`` is ``Refused`` if and only if
  ``realize`` raises :class:`~orpheus.transport.boundary_conditions.BoundaryRefused`. Five tables exist (``SN_REALIZERS``,
  ``DIFFUSION_REALIZERS``, ``CP_REALIZERS``, ``MOC_REALIZERS``, ``MC_REALIZERS``) and
  ``tests/gates/transport/test_admission_matrix.py`` pins every shared value class at every site (111 rows, 90 of them over the five real tables, read
  at the carve tip). Gates: ``test_one_door`` in ``tests/gates/transport/test_admission_mechanism.py`` and the pin. Derivation and the
  tables: :ref:`bc-aor-admission`.
- **A refusal is a statement about the declaration, with its axis and its issue as data.** A refusal lies on one of three
  axes (angular resolution, topology, state cone) or is *plumbing*, which has no theory behind it and therefore owes the issue
  that tracks the missing machinery (``Refused(PLUMBING, ...)`` without an issue does not construct). Today's limits stand, and
  a wrong answer stays a visible failure: CP's vacuum-bounded eigenvalue exceeds its white-bounded one (#546), carried by three
  strict ``xfail`` rows in ``tests/gates/cp/test_cp_vacuum_not_above_white.py``, and MC on a cylinder is
  ``Refused(TOPOLOGY, #547)``, carried by strict ``xfail`` rows in ``tests/gates/mc``.
- **A forward solve and an adjoint solve share one binding.** Every shared kernel is reciprocal in the trace metric
  :math:`M = \mathrm{diag}(w|\mu|)`, the adjoint of a deck is the deck of its partner face, and the inflow :math:`q` is not
  adjointed: ``solve_sn_adjoint_fixed_source`` never reads it, while ``solve_sn`` and ``solve_sn_adjoint`` refuse a declared
  source. Gates: ``test_pairing_of_and_its_adjoint`` in ``tests/gates/geometry/test_deck.py``, ``test_reciprocity`` in
  ``tests/gates/sn/architecture/test_realized_boundary.py``, ``test_the_adjoint_ignores_q`` in
  ``tests/gates/sn/architecture/test_sn_boundary_binding.py``. Statement: :ref:`bc-aor-adjoint`.
- **The inflow** :math:`q` **travels the boundary-source channel, never an operator slot.** A source is an array on
  :math:`\Gamma_-`, sampled by the method's table through the one verb
  :meth:`~orpheus.transport.boundary_conditions.RealizerTable.sample_inflow` and carried as the boundary leaf of the composite
  source :class:`~orpheus.transport.fields.AngularBoundarySourceSink`. Vacuum and a prescribed inflow realize to the same zero
  map and differ only in :math:`q`. Gates: ``test_the_functor_law`` and ``test_the_source_verb`` in ``tests/gates/transport/test_admission_mechanism.py``,
  ``tests/gates/transport/test_boundary_source_from_specs.py``. Derivation: :ref:`bc-aor-inflow`, :ref:`bc-affine-source-channel`.
- **The face ordinate partition is THREE-way**, not two: :math:`\{1..N\} = I_f \sqcup O_f \sqcup T_f` with
  :math:`T_f = \{|\Omega\cdot\hat n| \le \texttt{TANGENTIAL\_EPS}\}` (:eq:`ordinate-partition-inflow-outflow`). **"Not inflow"
  is NOT "outflow".** Measured: a cylinder under ``product(n_mu=2, n_phi=4)`` carries **4 of 8** ordinates tangential at
  ``xmax``; ``gauss_legendre(5)`` carries 1; every ``lebedev`` carries 4 to 8; only ``gauss_legendre(4)`` is the clean two-way
  case, so **the slab is the unrepresentative fixture.** Measured too: :math:`|\Gamma_+| = |\Gamma_-|` on every quadrature and
  face in the tree, so **a shape assertion cannot distinguish** :math:`\Gamma_+ \to \Gamma_+` **from**
  :math:`\Gamma_+ \to \Gamma_-` (the measured functional has the error class in its stabiliser, ``vv-principles`` mode 12).
- **The realized boundary is a first-class sibling operator** :math:`B` **in the SN algebra.** For every SN geometry the
  realized per-face images are assembled into the whole-trace
  :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`, the :math:`A_{ss}` boundary block of the canonical loss operator
  :math:`(L_{\rm full} + C - S - F - B)`. The reflection :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}` is not re-applied
  inside the streaming sweep for any geometry; it is delivered as the off-diagonal :math:`-B` source term and the outer Krylov
  or source-iteration loop drives the boundary consistency residual to zero. The block-matrix derivation is at
  :ref:`bc-extraction`.
- **Whether a configuration sweeps in one pass is a property of the whole face configuration, never of a single condition.**
  No condition and no deck carries a sweep-cycle flag: a boolean on the boundary *kind* cannot express a
  configuration-dependent property, since ``reflective|vacuum`` is acyclic while ``reflective|reflective`` is not. The criterion
  is a strongly-connected-component decomposition of the :math:`(\text{face}, \text{ordinate})` trace digraph, computed by
  :mod:`orpheus.derivations.discrete.sn.sweep_acyclicity` and gated by ``tests/gates/sn/sweep/test_sweep_acyclicity.py``. See
  :ref:`bc-sweep-cycle`.

.. admonition:: V&V status

   This page states structure, not a solver claim: its equations are definitional or are theorems about the types, and no
   equation on it asks for a closed-form or MMS reference. The verification load is carried by three kinds of gate.

   - **Foundation gates** pin the structural laws (the closed deck family, the closed chain sum, the cone as a value predicate,
     one door for admission and realization, the exact cover of a binding): ``tests/gates/geometry/test_deck.py``,
     ``tests/gates/geometry/test_colouring_symmetry.py``, ``tests/gates/transport/test_wall_condition.py``,
     ``tests/gates/transport/test_wall_condition_algebra.py``, ``tests/gates/transport/test_response_network.py``,
     ``tests/gates/transport/test_binding.py`` and ``tests/gates/transport/test_admission_mechanism.py``.
   - **Reference rows** pin values typed by hand from the plan of record: the admission pin
     ``tests/gates/transport/test_admission_matrix.py``, which covers every shared value class at every site of every method's
     table (seven value classes; the "white 1.0" cells of SN and diffusion and the whole deck and ``ZeroFluxCondition`` table are
     not in it, and their printed values were confirmed by a probe, :ref:`bc-aor-admission-probe`).
   - **A regression capture** is not a reference: ``tests/gates/sn/architecture/test_wave2_carry.py`` holds the eigenvalues,
     fluxes and dense face matrices (SN on four bodies, an S16 slab) captured from ``main`` at ``5d04c296`` by a script before the
     carve and reproduced bit-equal after it. It detects drift and says nothing about correctness. Only its seven
     ``bc_equivalence`` snapshot images are independent derivations from the law's equation: they are the ones re-anchored from
     recording to reference on 2026-08-01 (:ref:`bc-snapshot-reanchoring`).
   - **Honest reds**: a defect that the tree has and no fix addresses is carried by a strict ``xfail`` that fails when the
     defect is repaired (#546, #547), never by a masked admission. The two strict #654 xfails in ``tests/gates/mc/test_properties.py`` (a slab half-cell's mirror plane and a Wigner-Seitz cell's white outer face, both refused by MC's table; #180 tracks the capability) are honest reds too.

   Last run at the carve tip ``471c4737``, serial, ``python -O -m pytest`` over the gate files named above together with
   ``tests/gates/sn/architecture/test_sn_realizer_table.py``, ``tests/gates/sn/architecture/test_realized_boundary.py``,
   ``tests/gates/sn/architecture/test_sn_boundary_binding.py``, ``tests/gates/cp/test_cp_vacuum_not_above_white.py``,
   ``tests/gates/mc/test_mc_boundary_table.py``, ``tests/gates/mesh/test_face_roles.py``,
   ``tests/gates/mesh/test_hollow_inner_law.py`` and ``tests/gates/numerics/test_ordinate_pairing.py`` (17 files in all): 494 passed and
   3 xfailed (the three CP rows of #546). `[M]` command: ``python -O -m pytest -p no:cacheprovider -q`` over those 17 files, serial,
   ``494 passed, 3 xfailed, 6 warnings in 41.96s``. The law gates of the laws section that are not among the 17 files
   (``tests/gates/test_layer_imports.py``, ``tests/gates/transport/test_wave2_s4_tables.py``,
   ``tests/gates/transport/test_content_identity_transport.py``, ``tests/gates/geometry/test_content_identity_geometry.py``,
   ``tests/gates/sn/operators/test_sn_boundary_realizer.py``, ``tests/gates/transport/test_boundary_source_from_specs.py``,
   ``tests/gates/cp/test_cp_realizer_table.py``, ``tests/gates/sn/sweep/test_sweep_acyclicity.py``) are a second run of the same
   invocation: ``785 passed, 1 warning in 19.52s``.

.. _bc-architecture-of-record:

The architecture of record
==========================

This section is the method-agnostic theory of how a transport problem closes its boundary, as the code of the carve tip states it.
Everything in it is true of the tree at ``471c4737`` and is read from the code named beside each part. The older sections of the
page, after the crosswalk, carry the mathematics that was found while the architecture was a family of string-tagged trace laws; the
parts of that mathematics that still hold are cross-linked where this text leans on them.

The architecture has five objects, each a type, and one rule that places them:

- the **face role** (:class:`~orpheus.geometry.deck.MirrorDeck`, :class:`~orpheus.geometry.deck.TranslationDeck` or an integer
  tag): what a boundary face *is*, a property of the geometry and of nothing else;
- the **wall condition** (:class:`~orpheus.transport.boundary_conditions.WallCondition`): what the exterior behind a tag *does*,
  one affine trace map on the sub-Markov simplex, supplied where the method is known;
- the **binding** (:func:`~orpheus.transport.boundary_conditions.bind`): the join of tags to conditions, lawful only when it is an
  exact cover that the problem's colouring group leaves invariant;
- the **realizer table** (:class:`~orpheus.transport.boundary_conditions.RealizerTable`): per method, the one partial map from
  declarations to operators, whose domain *is* the admission;
- the **inflow terms** (:class:`~orpheus.transport.boundary_conditions.InflowTerms`): the affine part :math:`q`, which realizes
  into an array and never into an operator.

The rule is the ontological one of the project: each concept sits in the lowest layer whose vocabulary defines it. A deck is
geometry (it names isometries and faces), so it lives in ``orpheus/geometry/deck.py`` and imports no transport vocabulary. A
condition is transport (it names half-traces, kernels and sources), so it lives in ``orpheus/transport/boundary_conditions/``. A
realized operator is a method's, so it lives with its method. ``tests/gates/test_layer_imports.py`` enforces the layers.

.. _bc-aor-trace-closures:

The face trace and its two closures
-----------------------------------

Let :math:`F` be the set of boundary faces of a domain :math:`D`; a face is named by the vocabulary of
``orpheus/geometry/faces.py`` (``xmin``, ``xmax``, ``ymin``, ...). For a face :math:`f` with outward unit normal :math:`\hat n_f`
the **outflow half-trace** is :math:`\Gamma_+(f) = f \times \{\Omega : \Omega\cdot\hat n_f > 0\}` and the **inflow half-trace**
:math:`\Gamma_-(f)` is its complement :math:`\{\Omega\cdot\hat n_f < 0\}`, both carrying the trace measure
:math:`|\Omega\cdot\hat n_f|\,d\Omega\,dA`. Ordinates with :math:`|\Omega\cdot\hat n_f|` at most ``TANGENTIAL_EPS`` belong to neither
half (:eq:`ordinate-partition-inflow-outflow`). A transport equation on :math:`D` is closed when :math:`\gamma_-\psi` is given on
every face, and :math:`\gamma_-\psi|_f` is determined in exactly two ways, which are different objects.

**The domain continues across the face (a deck).** There is an isometry :math:`g` of the covering space that carries the partner
face :math:`f'` onto :math:`f` (:math:`f' = f` for a mirror), and

.. math::

   \gamma_-\psi|_f \;=\; G_g\,\gamma_+\psi|_{f'}, \qquad (G_g\varphi)(x) = \varphi(g^{-1}x),

with :math:`G_g` the Koopman operator of :math:`g`. It is invertible, preserves the trace measure, and has no parameter beyond
:math:`g`: no amplitude, no source, nothing to scale. Its return fraction is exactly one.

**The domain ends at the face (a boundary condition).** The exterior is modelled by an affine map on the face's own half-traces,

.. math::

   \gamma_-\psi|_f \;=\; R\,\gamma_+\psi|_f + q, \qquad R : L^2(\Gamma_+) \to L^2(\Gamma_-),\quad q \ge 0 \text{ on } \Gamma_-,

with :math:`R` a kernel and :math:`q` a function on the inflow half-trace. Every kernel the methods share is written with
:math:`\hat n_f` alone, so it carries no axis and no outward sign: the face supplies both, and one condition object is correct on
every face of its tag.

.. list-table:: The two closures side by side
   :header-rows: 1
   :widths: 24 38 38

   * - property
     - deck (the domain continues)
     - condition (the domain ends)
   * - what it is
     - the Koopman operator of an isometry; a relabelling of the trace
     - an affine map :math:`(R, q)`; a kernel plus a source
   * - parameters
     - none (the face fixes the plane, the partner and the direction)
     - the cone coordinates :math:`(s, d)` and the inflow weights
   * - return fraction
     - exactly 1
     - :math:`\rho = s + d \le 1` for the shared family
   * - source
     - none
     - :math:`q \ge 0` on :math:`\Gamma_-`
   * - may be scaled or added
     - no: ``0.5 * MirrorDeck()`` is Python's own ``TypeError``
     - yes, inside the cone (:ref:`bc-aor-wall-condition`)
   * - lives in
     - ``orpheus/geometry/deck.py``
     - ``orpheus/transport/boundary_conditions/``
   * - declared by
     - the geometry, as the role of a face
     - the problem, bound to a tag where the method is known

**The criterion that separates them is the two-test criterion of** :ref:`bc-factor-roles`. An operator that fails multiplicativity
:math:`G(\psi\varphi) = (G\psi)(G\varphi)` is a kernel and belongs to a condition; one that passes it is still a condition unless
it is the deck transformation of an actual quotient of the physical domain. A wall is not a quotient (the domain does not continue
on its far side; the table of which boundary is an orbifold, a torus or neither is :ref:`bc-factor-quotients`), so its specular return,
although a permutation, is constitutive. That is why a polished wall and a symmetry
plane have the same matrix and different types, and why neither type can be written into the other's slot
(``test_a_condition_is_not_a_deck_slot``, ``test_a_deck_is_not_an_operand`` in
``tests/gates/transport/test_wall_condition_algebra.py``).

.. _bc-aor-deck-family:

The deck family is derived, not chosen
--------------------------------------

**Which isometries can glue a face.** The argument, `[R]` reasoned (the code states its result in ``admitted_decks`` and
``test_a_curved_face_admits_no_deck`` pins that result as a table, which is not a proof): a deck is an isometry :math:`g` of the
covering space. The chart of the problem is the
quotient of that space by the chart's continuous symmetry group :math:`G_c` (translations tangent to a slab, rotations of a
cylinder or a sphere), so :math:`g` must descend to the chart, that is, normalise :math:`G_c`; and :math:`G_c` itself acts
trivially on the chart. The admissible decks are therefore

.. math::

   W \;=\; N(G_c)\,/\,G_c .

On a Cartesian axis :math:`W` is generated by the reflection in the face's own plane and the translation to the opposite face; these
are the two minted decks. On the radial axis of a cylinder or a sphere :math:`G_c` absorbs the rotational symmetry and :math:`W` is
trivial. The consequence is a design rule: **no 1-D curved face carries a deck**. A void cavity gives, in the 1-D radial reduction, the same direction map as a
specular wall (a ray leaving inward re-enters with the opposite cosine), but the cavity's surface is a true boundary, because what
lies beyond it is not known to the geometry: it could be specular, a vacuum (an absorber in the centre) or a prescribed inflow
(a beam through the centre). A reflecting curved wall, and the hollow cavity's inner face (ruling D5), are therefore *tags* bound to
:meth:`~orpheus.transport.boundary_conditions.WallCondition.specular` at 1.0, never decks.

.. list-table:: The classes of isometry and the decks that realize them (Coxeter)
   :header-rows: 1
   :widths: 20 12 12 12 44

   * - class
     - :math:`E^1`
     - :math:`E^2`
     - :math:`E^3`
     - the deck
   * - translation
     - yes
     - yes
     - yes
     - :class:`~orpheus.geometry.deck.TranslationDeck` (minted)
   * - reflection
     - yes
     - yes
     - yes
     - :class:`~orpheus.geometry.deck.MirrorDeck` (minted)
   * - rotation
     - no
     - yes
     - yes
     - ``RotationDeck``, deferred: a rotation's partner face is not derivable from the face (a quarter turn about a corner pairs
       ``xmin`` with ``ymin`` or ``ymax``)
   * - glide reflection
     - no
     - yes
     - yes
     - ``GlideDeck``, deferred: it carries a translation part
   * - screw motion
     - no
     - no
     - yes
     - ``ScrewDeck``, deferred: it carries an axis and a translation part
   * - rotoreflection
     - no
     - no
     - yes
     - ``RotoreflectionDeck``, deferred: it carries an order

The non-identity isometries of :math:`E^n` with :math:`n \le 3` fall into exactly six types (the classes up to conjugacy and parameter: a translation length, a rotation angle and so on; Coxeter, *Introduction to Geometry*, chapters 3 and 7), two in dimension
one, four in dimension two, six in dimension three. :class:`~orpheus.geometry.deck.IsometryClass` enumerates them and
``ISOMETRY_CLASSES_BY_DIMENSION`` records which exist where. The classification is complete because the four classes without a
minted type are rows of ``DEFERRED_DECKS``, keyed by the enumeration, each naming the type it will have and the machinery that mints
it (``minted | deferred == set(IsometryClass)``, ``test_isometry_classification_is_closed`` in
``tests/gates/geometry/test_deck.py``). Minting the four now would have minted four types with guessed fields and four unreachable
arms in every method.

**A deck has no field.** :class:`~orpheus.geometry.deck.MirrorDeck` and :class:`~orpheus.geometry.deck.TranslationDeck` store
nothing and compare equal to themselves by type. The face fixes the mirror plane, the partner and the direction of the translation,
so a deck that disagrees with its face (the class of defect in which a mirror carried an axis different from its face's) cannot be
spelled. The position of a mirror plane is not stored either (a 1-D mirror acts on directions only, and the offset is invisible to
every discrete quantity), and the length of a translation is not content (a deck does not know the domain's extent).

**What a realizer reads of a deck** is its pairing, derived and never stored:
:func:`~orpheus.geometry.deck.pairing_of` takes a deck, a face and a dimension and returns a
:class:`~orpheus.geometry.deck.FacePairing` ``(face, partner, motion)`` through one exhaustive ``match`` that ends in
``assert_never``, so the day a third deck is minted the type checker fails that match until it decides
(``test_every_deck_match_ends_in_assert_never``). The motion is :func:`~orpheus.geometry.deck.mirror_motion` for a mirror (the
reflection whose plane is normal to the face's axis; the partner is the face itself) and ``translation_motion`` for a translation
(signed by the face, so ``xmin`` carries :math:`-\hat e_x` and ``xmax`` carries :math:`+\hat e_x`, and the partner is
``face_opposite``). :class:`~orpheus.geometry.deck.FacePairing` states the whole law, so a hand-built value cannot spell what
``pairing_of`` cannot: a face paired with itself carries exactly the reflection in its own plane; a pair of distinct faces joined by
a translation is the periodic gluing, with the partner opposite and the translation along the face's axis with its outward sign; a
distinct-face pair joined by anything else belongs to a deferred class and is refused. The pairing seen from the partner is the
adjoint's pairing, :math:`G_g^\sharp = G_{g^{-1}}`:

.. code-block:: python

   from orpheus.geometry.coord import CoordSystem
   from orpheus.geometry.deck import MirrorDeck, TranslationDeck, admitted_decks, pairing_of

   for deck in (MirrorDeck(), TranslationDeck()):
       for face in ("xmin", "xmax"):
           pairing = pairing_of(deck, face, dimension=1)
           partner_pairing = pairing_of(deck, pairing.partner, dimension=1)
           print(type(deck).__name__, face, pairing.partner, pairing.motion.translation.tolist(), pairing.inverse == partner_pairing)
   print(sorted(d.__name__ for d in admitted_decks(CoordSystem.CARTESIAN, 0)))
   print(sorted(d.__name__ for d in admitted_decks(CoordSystem.CYLINDRICAL, 0)))

.. code-block:: text

   MirrorDeck xmin xmin [0.0] True
   MirrorDeck xmax xmax [0.0] True
   TranslationDeck xmin xmax [-1.0] True
   TranslationDeck xmax xmin [1.0] True
   ['MirrorDeck', 'TranslationDeck']
   []

The last two lines are :func:`~orpheus.geometry.deck.admitted_decks` on the first axis of a slab and of a cylinder: the quotient
:math:`W` as a table. ``mirror_motion`` is the one axis-letter to mirror spelling; ``orpheus/geometry/symmetry.py`` delegates to it.

**A face role is parsed once, at the chart.** :func:`~orpheus.geometry.deck.parse_face_role_on` is the one place a deck on a face
whose axis admits none is refused, with a :class:`~orpheus.geometry.deck.DeckError` that names the quotient; the geometry, both
meshes and both axis classes read their roles through it, and :func:`~orpheus.geometry.deck.check_pairings_close` requires that a
deck gluing a face to a different face have the same deck type there.

.. _bc-aor-wall-condition:

The wall condition: an affine trace map on the sub-Markov simplex
-----------------------------------------------------------------

**The response is a tensor network of typed factors.** The linear part :math:`R` of a wall is a network of factor tensors on the
rungs of one face's trace ladder, the ladder of :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`: the
outflow half-trace :math:`\Gamma_+` (tier ``"outflow"``), the partial current :math:`S` (tier ``"current"``) and the inflow
half-trace :math:`\Gamma_-` (tier ``"inflow"``); the tier marker is numerics' own ``TraceRole``, and "embedding" is not an
operator name in this project (the reading of an operator as a tensor network is :doc:`/theory/foundations/operator_tensor_network`).
The factors are field-less frozen dataclasses (a factor is a type, its realization is the method's):

.. list-table:: The factors of the response network
   :header-rows: 1
   :widths: 24 12 12 52

   * - factor
     - from
     - to
     - action
   * - :class:`~orpheus.transport.boundary_conditions.MirrorPermutation` :math:`P`
     - :math:`\Gamma_+`
     - :math:`\Gamma_-`
     - :math:`(P\varphi)(x,\Omega) = \varphi(x, \Omega - 2(\Omega\cdot\hat n_f)\hat n_f)`: the ordinate permutation of the face's
       mirror plane
   * - :class:`~orpheus.transport.boundary_conditions.CurrentRetraction` :math:`R_t`
     - :math:`\Gamma_+`
     - :math:`S`
     - :math:`R_t\varphi = \sum_{\Omega\cdot\hat n_f > 0} |\Omega\cdot\hat n_f|\,w\,\varphi`: integration against the trace measure
   * - :class:`~orpheus.transport.boundary_conditions.IsotropicSection` :math:`E`
     - :math:`S`
     - :math:`\Gamma_-`
     - :math:`E J = J / M_-`, with :math:`M_- = \sum_{\Omega\cdot\hat n_f < 0} |\Omega\cdot\hat n_f|\,w` the INFLOW mass: uniform in
       :math:`\Omega`

:math:`R_t` and :math:`E` are the project's canonical split pair, a retraction and its section with
:math:`R_t^{\rm in}\circ E = \mathrm{id}` on :math:`\Gamma_-` (in production as ``AxisRetractionOperator`` and
``AxisSectionOperator`` of ``orpheus/numerics/axis_marginal.py``), where :math:`R_t^{\rm in}` is the retraction over the inflow
rung. Two chains compose from :math:`\Gamma_+` to :math:`\Gamma_-`, the **specular chain** :math:`P` and the **diffuse chain**
:math:`E\circ R_t`
(:class:`~orpheus.transport.boundary_conditions.SpecularChain`, :class:`~orpheus.transport.boundary_conditions.DiffuseChain`; the
factors are listed in the order they are applied). The weights of the two chains are the scalars :math:`s` and :math:`d`. The linear
part of a wall is the weighted sum

.. math::

   R \;=\; s\,P \;+\; d\,(E\circ R_t), \qquad (s, d) \in K = \{\,s, d \ge 0,\ s + d \le 1\,\},

a point of the sub-Markov simplex :math:`K` (the positive quadrant cut by the sub-Markov constraint :math:`s + d \le 1`; a compact convex set, and not a cone, since :math:`2\,\mathrm{specular}(0.6)` leaves it), with :math:`q \ge 0` a function on :math:`\Gamma_-`, which is a cone. The code names the admissible set "the cone" (``_cone_refusal``, ``cone_point``, the state-cone axis); the name is inherited and the set is the simplex. The network is
:class:`~orpheus.transport.boundary_conditions.ResponseNetwork`: a canonical weight map keyed by chain instances, whose weights
*are* the point :math:`(s, d)` (``cone_point`` reads it by one exhaustive match, so the chart is total and not lossy), with no
collapse of a one-term sum to a leaf and so no stored intermediate rounding.

**Why the chains are a closed sum.** ``Chain = SpecularChain | DiffuseChain`` is a ``typing.Union``, which has no creation hook, so
an open family could not be enforced. The family is closed the way the decks are: every reader of a chain does an exhaustive
``match`` ending in ``assert_never``, so a third chain is an edit of the alias and every ``match`` reddens until it handles the
new member; a chain enters with its declared ``unit_return`` (the return fraction at weight one) and the left-eigenvector gate
below. :func:`~orpheus.transport.boundary_conditions.chain_composes` decides whether a tuple of factors composes from ``"outflow"``
to ``"inflow"``, and the module applies it to every member of the closed sum when it is imported, so a chain such as
``(IsotropicSection, CurrentRetraction)`` is unspellable (law L-M, ``test_chains_compose_and_the_family_is_closed`` in
``tests/gates/transport/test_response_network.py``). A method-only kernel (the neutron-guide wall, an anisotropic SN kernel, zero
flux) is another :class:`~orpheus.transport.boundary_conditions.BoundaryCondition` type and never a chain.

**The return fraction, and why it is exact.** The partial current of the outflow trace is :math:`J^+(\varphi) = R_t\varphi` and that
of an inflow trace is :math:`J^-(\chi) = R_t^{\rm in}\chi`. For :math:`\varphi \ge 0` on :math:`\Gamma_+`,

.. math::

   J^-(P\varphi) = J^+(\varphi), \qquad
   J^-\big(E\,R_t\varphi\big) = R_t^{\rm in}E\,(R_t\varphi) = R_t\varphi = J^+(\varphi),

the first because :math:`P` is a measure-preserving permutation (``ordinate_pairing`` certifies it at realization, below), the second
because :math:`R_t^{\rm in}\circ E = \mathrm{id}`. Hence

.. math::

   J^-(R\varphi) \;=\; s\,J^+(\varphi) + d\,J^+(\varphi) \;=\; (s + d)\,J^+(\varphi) \;=\; \rho\,J^+(\varphi),
   \qquad \rho = s + d .

Written against the trace measure :math:`m_\pm = (w|\mu|)` on the two rungs, this says that :math:`m_-` is a *left eigenvector* of
the realized map, :math:`R^{\mathsf T} m_- = \rho\,m_+`. For the shared family it holds by the construction, with one precondition
the shipped code once failed to meet: :math:`E` must divide by the INFLOW mass :math:`M_-`. A section that divided by the outflow
mass :math:`M_+` would give :math:`\rho_{\rm diffuse} = d\,M_-/M_+`, silently. The carry is bit-identical because the two masses
agree bitwise on every quadrature the census reads (``test_the_mass_census_inflow_equals_outflow_on_every_registry_rule`` in
``tests/gates/sn/architecture/test_realized_boundary.py``, which reads 70 of 70 (rule, face) pairs over 35 rules of four registry
families, not every rule the families can build), so the repair moved no number. Because the masses agree on the registry, no registry
fixture can see the divisor: the witness of the inflow-mass divisor at the three sites is
``test_the_rider_on_a_trace_whose_masses_differ`` in the same file, the only gate that reddens when ``E`` divides by :math:`M_+`
(`[M]` an in-process mutation of ``RankOneReemission._derive_inner`` to ``float(np.sum(out_weights))``: 1 failed, 17 passed over
``test_sn_realizer_table.py`` and ``test_realized_boundary.py``; the same divisor scaled by :math:`1+10^{-7}` reddens 17 of 18,
every rider row, and leaves the census green). The consumers of :math:`\rho` (the curvilinear corner action of ``SNBoundaryOperator`` and the leakage term of the
eigenvalue, whose omission for a partial reflector was ``ERR-094``) read it from the realized type
(:attr:`~orpheus.sn.boundary.realized_operators.RealizedBoundary.return_fraction`), never from the declaration. A ratio gate (":math:`J^-/J^+` of one positive
outflow equals ``return_fraction``") cannot discriminate a wrong :math:`\rho`, since a type that reports the one measured ratio
passes; the law is therefore gated as the array identity
``inner.apply_transpose(w|mu| on Gamma_-) == rho * (w|mu| on Gamma_+)`` (``array_equal`` for the specular image; for the diffuse
image up to the split pair's own residual, at most :math:`N/3` ulp for :math:`N` slots)
(``test_the_left_eigenvector_rider_through_the_realizer`` in ``tests/gates/sn/architecture/test_sn_realizer_table.py``; it pins the identity on the registry rules and is blind to :math:`M_+` versus :math:`M_-`). A
method-only kernel that is not built from this pair owes the same gate, and where it fails it has no scalar return fraction and its
consumers must compute :math:`J^+ - J^-(R\psi)` on the actual outflow.

**One type, with facades (rulings D1 and D2).** A general condition with a full parameter range is the type, and a conventional name
is a thin named constructor over it, never a second type (ruling D1: white folds into the albedo family). Combinations of
conditions are one normalised type with non-negative coefficients and total return at most one, checked at construction (ruling D2);
signed combinations stay values but cannot be bound. The consequence is that the two spellings of one physical wall cannot differ:
:math:`R = 0` is one point of :math:`K` however it is reached, so the vacuum has one value and one digest. The constructors of
:class:`~orpheus.transport.boundary_conditions.WallCondition` are:

.. list-table:: The facades of ``WallCondition``
   :header-rows: 1
   :widths: 34 26 40

   * - constructor
     - network and source
     - note
   * - ``vacuum()``
     - empty network, no source
     - the zero network :math:`(0, 0)`
   * - ``specular(albedo=1.0)``
     - :math:`s = \text{albedo}`
     - ``specular(0.0)`` is the vacuum
   * - ``white(albedo=1.0)``
     - :math:`d = \text{albedo}`
     - the first positional argument is the albedo; ``white(0.0)`` is the vacuum
   * - ``uniform(value)``
     - empty network, uniform ``q = value``
     - a prescribed inflow; the facade is named after the coordinate ``InflowTerms.uniform``
   * - ``framed(profile, weight=1.0)``
     - empty network, sampled ``q = weight * profile``
     - a prescribed inflow written as a profile (:ref:`bc-aor-inflow`)
   * - ``from_amplitudes(specular, diffuse)``
     - :math:`(s, d)`
     - the plain-data edge: the section of ``cone_point``, taking two primitives and never the reference's ``WallAmplitudes``

**Membership is a value predicate.** One function, ``_cone_refusal``, defines "inside" (the simplex :math:`K` and the source cone :math:`q \ge 0`): every network weight is
non-negative, the return fraction (:math:`\sum_c w_c\,\mathrm{unit\_return}(c)` by ``math.fsum``, correctly rounded) is at most one,
and every inflow weight is non-negative. ``WallCondition.__post_init__`` raises what it returns, so a wall outside the cone does not
exist, however it is spelled. The arithmetic ``+``, ``-``, unary ``-`` and scalar ``*`` (one definition, shared with
:class:`~orpheus.transport.boundary_conditions.WallDisplacement`) computes the pair of the result and exits through one function,
which returns a ``WallCondition`` if and only if the pair lies in the cone and otherwise a ``WallDisplacement``: a point of the affine
hull of the cone times the source space, a sensitivity direction (the parameter coordinates of a perturbation study), never
bindable, refused by :class:`~orpheus.transport.boundary_conditions.BoundaryConditions` with a ``TypeError`` naming the type. The
sign of a scalar and the presence of ``-`` decide nothing: ``x + h * delta`` for a finite difference re-enters the cone without
re-typing. There is no division (no consumer, and it would leave the dyadic coordinates), and a deck or a method-only condition
defines no operator at all (``WallCondition.vacuum() + TranslationDeck()`` is Python's own ``TypeError``, law L-F). The code below
was run at the carve tip, and its output is pasted:

.. code-block:: python

   from orpheus.geometry.coord import CoordSystem
   from orpheus.sn.boundary.realizer import SN_REALIZERS
   from orpheus.transport.boundary_conditions import BoundaryError, BoundaryRefused, FaceSite, WallCondition as W

   print(W.specular(0.0) == W.white(0.0) == W.vacuum())
   mixed = W.specular(0.5) + W.white(0.25)
   print(type(mixed).__name__, mixed.cone_point, mixed.return_fraction)
   print(type(2.0 * W.specular(0.6)).__name__)
   try:
       W.specular(1.2)
   except BoundaryError as error:
       print(type(error).__name__)
   outer = FaceSite("xmax", CoordSystem.CYLINDRICAL, 1)
   print(SN_REALIZERS.admit(W.specular(1.0), outer))
   print(SN_REALIZERS.admit(W.white(0.5), outer).axis, SN_REALIZERS.admit(W.white(0.5), outer).issue)
   try:
       SN_REALIZERS.require(W.white(0.5), outer)
   except BoundaryRefused as error:
       print(error.method, error.face, error.admission.axis.name, error.admission.issue)

.. code-block:: text

   True
   WallCondition (0.5, 0.25) 0.75
   WallDisplacement
   SubmarkovViolationError
   Realized(faithful=True)
   RefusalAxis.PLUMBING 653
   sn xmax PLUMBING 653

The first line is the ruling D1 ("one value, one digest"); the second shows that a mixed wall is an ordinary condition (the cone
admits it, ruling D2); the third that ``2.0 * specular(0.6)`` leaves the cone and becomes a displacement; the last three are the
admission of :ref:`bc-aor-admission` read through the door.

**Method-only responses stay their own types** and are never terms of the shared family.
:class:`~orpheus.diffusion.zero_flux.ZeroFluxCondition` is the Dirichlet idealisation :math:`\phi_\Gamma = 0`, which in the
partial-current basis is :math:`J^- = -J^+`, :math:`\mathcal{A} = -1`, so ``return_fraction == -1.0``, outside :math:`K`: no
non-negative angular flux produces it, which is why it is not a wall and why a table whose state is a non-negative flux refuses
it on the state-cone axis. It is face-local, has no fields (so every instance is equal to every other), and defines no arithmetic.
A type under the old name could not be minted without two nameable classes sharing a qualname, which the content-identity
relocation gate refuses and which takes no exemption, so the new type carries a fresh name.

**Face-local, with no orientation field.** ``WallCondition`` and the factors hold no ``axis`` and no ``outward_sign``: the realizer
reads the face's normal from the face name. Face-local means the content is stated in the face's own frame, so every isometry
carries it from face to face (:ref:`bc-aor-colouring`).

**What a reference may state (ruling D3).** The reference specification states a wall as
``WallAmplitudes(specular, diffuse)`` (``orpheus/specification/exterior.py``), exactly the coordinate pair :math:`(s, d)`, as plain
data; the one translation into production is ``WallCondition.from_amplitudes``, which imports neither the specification nor the
derivations. The reference vocabulary is closed: vacuum :math:`(0, 0)`, specular :math:`(s, 0)` or diffuse :math:`(0, d)`. A combined
wall (both coordinates positive), a prescribed inflow and zero flux are not stated by a reference; a combined pair raises
``NotImplementedError`` naming #650, which widens the vocabulary. The production type is wider than the reference's: a combined wall
is an ordinary ``WallCondition``.

.. _bc-aor-roles-binding:

Face roles, tags and the binding
--------------------------------

**A tag is a colour of the true boundary, the exact twin of a material id.** A geometry says which faces share an exterior model
and never what the model is: a face carries an integer tag, as a cell carries an integer material id
(``orpheus/data/materials.py``: an ``int``, not a ``bool``), and the model behind a tag is declared by the problem, where the method
is known. The role of a face is the sum ``FaceRole = DeckTransformation | int``, parsed by
:func:`~orpheus.geometry.deck.parse_face_role`: a deck, or an integer (numpy integers are coerced to ``int``, no sign restriction),
and anything else is a ``TypeError`` naming where. A ``bool`` is not a tag, for the reason it is not a material id.

The geometry declares its roles at its constructors, which take a role per end face:
:meth:`~orpheus.geometry.structured_geometry.StructuredGeometry.slab` takes ``left`` and ``right``, ``cylinder`` and ``sphere`` take
``outer`` and, for a hollow body, ``inner``. The inner face of a hollow body is a tag (ruling D5: a deck folds the domain onto
itself and the cavity's face does not, so it is a true boundary, and no void region stands in for it; in 1-D radial geometry no face
is a deck at all). A mesh stores the same roles, one per face over the inventory the coordinate system and the extent determine
(:class:`~orpheus.mesh.face_roles.FaceRoles`, built by ``FaceRoles.over``; ``face_inventory`` is the geometry's, written once in
``orpheus/geometry/faces.py``): along the first axis the faces are the boundary points of its edges (both ends on a slab or a hollow
body, the outer surface alone on a solid one, whose centre is an interior point and carries no role), and along every further
axis both ends. A role per face, and not per side, is the seed for a side whose faces carry different roles.

**The declaration of the exterior is a twin of** ``Materials``. :class:`~orpheus.transport.boundary_conditions.BoundaryConditions`
maps a tag to a :class:`~orpheus.transport.boundary_conditions.BoundaryCondition`: keys are ``int`` (a ``bool`` is a ``TypeError``),
every value must be a ``BoundaryCondition`` (a ``WallDisplacement``, a deck and a ``str`` are refused with a ``TypeError`` naming the
type: the one narrowing site for displacements), and the mapping is stored in declared order as a frozen mapping. One difference from
``Materials``: *empty is admitted*, because a carrier of decks only (a half-cell slab whose two faces are mirrors) has no tag.

**Binding joins them.** A binding is a map :math:`\beta : T \to \text{BoundaryCondition}` from the set :math:`T` of tags the faces
carry. :func:`~orpheus.transport.boundary_conditions.bind` takes a carrier, an object that owns the faces **and** their colouring
symmetry (:class:`~orpheus.transport.boundary_conditions.BoundaryCarrier`: ``face_roles`` and ``symmetry``), and the conditions.
Passing the carrier and not its two derived values makes a binding checked against another carrier's group unspellable (a G-set is
the group together with the faces it acts on). It refuses, in this order:

#. a tag a face carries with no condition: :class:`~orpheus.transport.boundary_conditions.UnboundTagError`, naming the face and the
   tag;
#. a condition whose tag no face carries: :class:`~orpheus.transport.boundary_conditions.UnusedConditionError` (the exact cover
   catches a typo that a reachable-subset rule would hide);
#. a condition that a non-identity element of the colouring group moves:
   :class:`~orpheus.transport.boundary_conditions.NonInvariantConditionError` (:ref:`bc-aor-colouring`).

The result, :class:`~orpheus.transport.boundary_conditions.Binding`, holds ``decks`` and ``tags`` (two typed maps with disjoint keys),
the ``declared`` conditions and the ``symmetry``; the join of a face to its condition is *derived* from ``(tags, declared)`` and never
stored, so a binding that joins a face to a condition of another tag cannot be built. The parse of the value re-checks the exact cover,
the invariance and the face domain, so a direct ``Binding(...)`` cannot build a state that is not invariant under the group it was
handed; that the group *is* the carrier's is ``bind``'s guarantee, not the type's. ``Binding.sources`` reads
``BoundaryCondition.source`` of every condition type, with no ``isinstance`` on a leaf. Every refusal is a
:class:`~orpheus.transport.boundary_conditions.BoundaryError`, a ``ValueError`` carrying ``.declaration``: the object, or its type,
that was refused, so a test compares the object and not a spelling. The code below was run at the carve tip:

.. code-block:: python

   from orpheus.geometry import StructuredGeometry
   from orpheus.geometry.deck import MirrorDeck
   from orpheus.mesh import CellsByCount, Mesher
   from orpheus.transport.boundary_conditions import WallCondition as W, bind, BoundaryError

   def mesh(left, right):
       geometry = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=left, right=right)
       return Mesher(geometry).partition((CellsByCount.uniform_volume(2), CellsByCount.uniform_volume(2))).mesh

   symmetric = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 0), left=0, right=0)
   print(symmetric.face_roles, len(symmetric.symmetry.elements))
   m = mesh(0, 1)
   print(dict(m.face_roles), len(m.symmetry.elements))
   b = bind(m, {0: W.vacuum(), 1: W.specular(1.0)})
   print(dict(b.tags), type(b.condition_of("xmax")).__name__)
   try:
       bind(m, {0: W.vacuum()})
   except BoundaryError as error:
       print(type(error).__name__, error.tag, error.face)
   print(dict(mesh(MirrorDeck(), 1).face_roles))

.. code-block:: text

   (0, 0) 2
   {'xmin': 0, 'xmax': 1} 1
   {'xmin': 0, 'xmax': 1} WallCondition
   UnboundTagError 1 xmax
   {'xmin': MirrorDeck(), 'xmax': 1}

The first line is a symmetric slab with one tag on both faces, whose colouring group has two elements (the identity and the mirror
:math:`x \mapsto 2 - x`); the second, the same slab with different tags on its two faces, whose group is the identity alone (the
roles differ). The third and fourth are the binding and the refusal of an unbound tag, and the last is a slab whose left face is a
mirror deck and whose right face a tag.

.. _bc-aor-colouring:

The colouring stabiliser, and when a binding is lawful
------------------------------------------------------

**The group.** The symmetry group of a problem is the subgroup of isometries that respect both the material colouring and the roles of
the faces,

.. math::

   H \;=\; \{\,h \in \mathrm{Isom}(D) \;:\; \mathrm{mat}(h\,c) = \mathrm{mat}(c)\ \text{for every region } c,\quad
   \mathrm{role}(h\,f) = h\cdot\mathrm{role}(f)\cdot h^{-1}\ \text{for every face } f\,\},

with :math:`h \cdot t = t` for a tag :math:`t` and conjugation for a deck. It is computed from the GEOMETRY's material runs and end
roles and never from the mesh's cell edges, so one ``(geometry, binding)`` pair is admitted or refused for every mesh, and refining a
mesh cannot change admission (the raw partition of a geometry, a mesh or a material mesh is handed to
:func:`~orpheus.geometry.colouring.colouring_symmetry`, which merges equal adjacent material ids once, so the group is a function of
the colouring and not of the cell partition).

**A binding is lawful** when, for every non-identity :math:`h \in H` that moves a face of tag :math:`t`,
:math:`h_*\beta(t) = \beta(t)`. Every shared kernel is face-local, so the equation holds for it by a theorem; it can fail only for a
condition whose content is stated in a frame the face does not supply, which is a prescribed inflow whose profile is written in
laboratory directions (:class:`~orpheus.transport.boundary_conditions.LabFrameInflow`, whose ``face_local`` is False, so that
``invariant_under`` holds for the identity alone). The test is conservative and deterministic: ``function o h == function`` is
undecidable for an open callable, and sampling the two faces at a tolerance would be a check keyed on nothing numerical. The authored
witness is the symmetric slab with tag 0 on both faces and ``WallCondition.framed(LabFrameInflow(f))``, refused at ``bind``; the
controls (two tags, a one-face tag, an asymmetric slab, a uniform source, a face-frame source) bind
(``test_bind_invariance`` in ``tests/gates/transport/test_binding.py``).

**"Face-local implies invariant" is a one-dimensional theorem.** On a Cartesian slab the group acts on faces by the mirror
:math:`x \mapsto a + b - x`, which has no tangential coordinate to reverse. In higher dimension a mirror fixing a face reverses or
permutes its tangent frame, so a ``FaceFrameInflow`` reading a tangential coordinate would not be invariant.

**What the value is.** :class:`~orpheus.geometry.colouring.ColouringSymmetry` holds ``elements`` with the identity first, one element
or two for a slab: the *stabiliser of the fundamental domain* (a slab whose two faces carry a ``TranslationDeck`` has an orbifold that
is a circle, which has more reflections than :math:`x \mapsto a + b - x`; the value stays the stabiliser until a full deck group
lands). The reflection is in :math:`H` when all three of these hold, for the breakpoints :math:`b_0 < \dots < b_n` of the material
*runs*:

#. every mirrored breakpoint :math:`(b_0 + b_n) - b_{n-i}` agrees with :math:`b_i` under the geometry's one coordinate-agreement
   convention, ``coordinates_agree``, at 64 ulp of the end magnitude :math:`\max(|b_0|, |b_n|)`;
#. the material runs read the same reversed, compared exactly (integer ids);
#. the two end roles are the same tag or the same deck type.

`[M]` The float rule was measured at WF2.1 (the memo ``measure_float_rule.md``, untracked scratch, not re-run for this page): the
exact comparison accepts 2 of 13 symmetric slabs, a width-relative tolerance refuses an offset slab, and 64 ulp at the end
magnitude accepts 13 of 13 and refuses 6 of 6 asymmetric ones (the genuine ceiling being 3.0 ulp). The gates that pin the rule are
the three ``test_colouring_symmetry_float_rule_*`` rows (``_accepts_the_symmetric``, ``_refuses_the_asymmetric``, ``_constant_is_64_ulp_at_the_end_magnitude``) in ``tests/gates/geometry/test_colouring_symmetry.py`` (they pass at the carve tip, in the
494-passed run of the V&V admonition above). The activation is an interior breakpoint, so every row
carries the palindromic id triple ``(0, 1, 0)`` (with one id the run merge deletes the interior breakpoints and the float rule has
nothing to read). A radial chart returns the identity only (:math:`G_c` absorbs the rotation and the quotient is trivial).

**The scope boundary, #658.** A multi-axis carrier (``Mesh2D``) has no group machinery yet and answers
:class:`~orpheus.geometry.colouring.UncomputedSymmetries`, the typed pessimistic answer, never an optimistic :math:`\{e\}`. ``bind``
then refuses, for any tag with two or more faces, a condition that is not face-local, with a message naming the 2-D colouring group as
the missing machinery. This is a declared scope boundary (``SCOPE-BOUNDARY[guard]``, revisit when a consumer of the problem's group
exists: a quotient or a cache): reducing the group instead of refusing is not built, because nothing consumes the problem's group
today, so a reduced group would have no reader. The issue that tracks the 2-D group is #658.

.. _bc-aor-admission:

Realization is a partial map, and admission is its domain
---------------------------------------------------------

**The definition.** A transport method :math:`M` is a discretization of the trace, a projection :math:`\Pi_M : \Gamma \to \Gamma_h`
onto whatever :math:`M` can represent (:ref:`bc-method-realizability`). Realization is then a partial map

.. math::

   \mathrm{real}_M : (\text{declared},\ \text{site},\ \text{method space}) \;\longmapsto\; B_f : \Gamma_{+,h}(f') \to \Gamma_{-,h}(f),

with :math:`\mathrm{real}_M(g) = G_{g,h}` for a deck and :math:`\mathrm{real}_M(\beta) = R_{f,h}` for a condition, extended linearly
over combinations. The face action that a consumer composes from a realized image is
:math:`B_{\rm face} = \iota_-\circ\mathrm{real}\circ\gamma_+`, the narrowed trace maps of :ref:`bc-domain-narrowing` (nothing is
computed on the outflow rows and then discarded). :math:`q` realizes separately into the boundary-source channel (:ref:`bc-aor-inflow`). **Admission answers the
question "is** :math:`\mathrm{real}_M` **defined at this value on this face, and how good is it".** It uses the vocabulary of
:ref:`bc-method-realizability` and :ref:`bc-refusal-axes`:

.. list-table:: The admission lattice
   :header-rows: 1
   :widths: 30 16 54

   * - member
     - tier
     - meaning
   * - ``Realized(faithful=True)``
     - 1
     - exact and faithful: :math:`R = \alpha I` for any :math:`\alpha` (a scalar commutes with every projection), and :math:`\alpha`
       is recoverable downstairs
   * - ``Realized(faithful=False)``
     - 2
     - exact but not faithful: the projection identifies laws that differ upstairs (at :math:`P_1`, specular and Lambertian both
       give :math:`J^- = \alpha J^+`)
   * - ``Approximated(name, issue)``
     - 3
     - not exact: the law depends on structure the method does not represent; the member holds a set of ``(name, issue)`` pairs
   * - ``Refused(axis, reason, issue)``
     - none
     - refused on exactly one axis of angular resolution, topology (spatial or topological) or state cone, or on ``PLUMBING``

:class:`~orpheus.transport.boundary_conditions.Admission` is a sealed class with exactly these three direct subclasses
(:class:`~orpheus.transport.boundary_conditions.Realized`, :class:`~orpheus.transport.boundary_conditions.Approximated`,
:class:`~orpheus.transport.boundary_conditions.Refused`), so ``isinstance``, ``match`` and annotations read as for a union while
``Admission.meet`` lives on the type. The lattice is ``Realized(True) > Realized(False) > Approximated > Refused``, and ``meet`` is
its one spelling: the lesser member; the meet of two ``Approximated`` is the union of their ``(name, issue)`` pairs (a set, because
a meet on a total order would drop one, and a pair, because which issue tracks which approximation is the information two separate
sets would lose); the meet of two ``Refused`` is the first argument's; the meet of no members is ``Realized(True)``. **A wall's
admission is the meet over its nonzero chains and its source**, so "both chains" is never an independent entry
(``test_both_chains_is_the_meet`` in ``tests/gates/transport/test_admission_matrix.py``).

``PLUMBING`` is the fourth thing, not an axis: *not built yet* has no theory behind it, so ``Refused(PLUMBING, reason)`` without an
issue number does not construct (``test_plumbing_needs_an_issue``). The axes are about the METHOD'S CARRIER: topology asks "can
its operator couple two faces" (diffusion's per-face scalar has no slot for a face pair, and refuses a translation deck however
well it resolves angle), state cone asks whether the value is in the method's state cone (an angular flux is non-negative, so SN
refuses :math:`\mathcal{A} = -1`; diffusion's partial current is signed and realizes it), angular resolution asks whether the
method can represent the kernel's angular structure. The present instance of the third is the quadrature itself, below.

**Two levels.** The *type level* is the front end's drop-down and the upper bound of the instance level:
:meth:`~orpheus.transport.boundary_conditions.RealizerTable.catalogue` returns, per key of the table, its ``status``. The *instance
level* is per face and per value: CP admits a mirror deck at a slab's ``xmin`` only; SN admits at a hollow cylinder's inner face only
the cavity its closure computes; and the quadrature's equivariance is certified through ``ordinate_pairing``. Both levels read ONE
table per method, the table the realizer dispatches through, so the query and the realization cannot disagree (law L-D).

**The table and its one door.** A :class:`~orpheus.transport.boundary_conditions.RealizerTable` belongs to one method and is keyed
by TYPE: the chain types ``SpecularChain`` and ``DiffuseChain`` (exactly the members of the closed sum), the marker ``ZeroNetwork`` (the
key of the zero image), the marker ``InflowTerms`` (the source channel, whose entry is a
:class:`~orpheus.transport.boundary_conditions.SourceRealizer`, the others
:class:`~orpheus.transport.boundary_conditions.ResponseRealizer`) and the method-only condition types. ``WallCondition`` itself is not
a key: it is realized generically by the functor, so a method never writes a branch over the shared family's cone point. Decks have
exactly one :class:`~orpheus.transport.boundary_conditions.DeckRealizer` per table (an exhaustive match over ``DeckTransformation``).

- **Realization is a functor on the network.** The image of the empty network is the method's zero map, the image of one chain is
  that chain's realizer's, and the image of a sum is the sum of the images, summed in ``CHAIN_ORDER`` so the realization does not
  depend on how the wall was built (law L-L, ``test_the_functor_law`` and, on the SN table, ``test_the_realizer_is_a_homomorphism``).
- **Registration is a class keyword plus** ``seal()``: ``class SNSpecularChain(ResponseRealizer, table=SN_REALIZERS,
  realizes=SpecularChain)`` registers at class creation and refuses a duplicate, a sealed table and a key that is not a legal unit.
  The table and every realizer that registers into it live in one module per method, which calls ``table.seal()`` after its last registration
  (the module's last line for SN, diffusion and MoC; mid-module in ``cp/solver.py`` and ``mc/solver.py``):
  anyone holding the table has executed every registration. An unsealed table refuses ``admit``, ``realize``, ``require`` and
  ``catalogue`` with :class:`~orpheus.transport.boundary_conditions.TableNotSealedError`.
- **Lookup is by exact type, never the MRO.** A subclass is a new declaration and owes its own entry; a method-only condition of
  another method is never realized silently as its parent.
- **One door.** :meth:`~orpheus.transport.boundary_conditions.RealizerTable.realize` calls
  :meth:`~orpheus.transport.boundary_conditions.RealizerTable.admit` (through ``require``) and raises
  :class:`~orpheus.transport.boundary_conditions.BoundaryRefused` if and only if the admission is ``Refused``; no realizer body
  raises on its own for an admission reason (law L-D, ``test_one_door``). The refusal carries the method, the face and the
  admission, so the axis and the issue are data. The one declared second door is POSING: the eigenvalue entries refuse a boundary
  source before any realization, a refusal of the question and not of the declaration (``ERR-095``).
- **The unregistered door.** A unit with no entry is ``Refused(STATE_CONE, ...)`` when its ``return_fraction`` is negative and the
  table's :class:`~orpheus.transport.boundary_conditions.StateCone` is ``NONNEGATIVE_FLUX`` (so SN refuses diffusion's zero flux
  without importing it: SN may not import diffusion), and ``Refused(PLUMBING, ..., table.unregistered_issue)`` otherwise, naming the
  admitted keys. Which issue tracks a method's missing entries is the method's to say: CP, MoC and MC cite #654 (their tables hold
  only what the retired string registries held), SN and diffusion cite #659 (their tables register every shared unit and meet
  this door only with another method's own condition).

The shared body that realizes a whole binding is :func:`~orpheus.transport.method.realize_boundary`: for each face of the binding in
inventory order, the declaration is the face's deck if it is glued and otherwise the condition its tag is bound to; a face with a
trace slot is realized through the method's :meth:`~orpheus.transport.method.TransportMethod.realize_boundary_role`, which calls
the method's own table; a face with no trace slot (the inner face of a hollow radial body: a tag with no inner-surface trace) meets
the table through ``require`` at the site level only, which raises the same refusal. The face inventory of the result is the trace
inventory by construction.

**The quadrature's equivariance is certified, and the certificate is the permutation.** The specular chain and the mirror deck
realize through the ordinate permutation of the face's mirror plane. ``ordinate_pairing(quadrature, motion)``
(``orpheus/numerics/quadrature/ordinate_pairing.py``) answers one question, keyed on that pair and nothing else: *is the ordinate
permutation this motion induces a certified pairing?* It folds onto ``Quadrature.ordinate_permutation`` (no bijective,
weight-preserving match becomes ``QuadratureNotClosedUnderMotionError``) and asserts three independent claims on the returned
permutation: the trace measure :math:`m_n = w_n|\mu_{a,n}|` satisfies :math:`m_{\pi(n)} = m_n` (``ERR-042``), :math:`\pi\circ\pi =
\mathrm{id}` (``ERR-044``), and every non-tangential ordinate pairs with the opposite sign class on the mirror's axis (``ERR-045``);
each passes on a table the other two reject, which is why none is dropped for another. It returns the permutation, so the realizer
uses the certified array and there is no certification beside a second derivation. A rule that is not closed under the mirror (an
odd-``n_phi`` product rule has no x-mirror closure) is refused by SN's table on the angular-resolution axis when a method space
is given, and admitted at site level without one (``test_a_certification_failure_names_its_axis_at_the_table`` in
``tests/gates/sn/architecture/test_sn_realizer_table.py``). The oblique reflection is a declared scope boundary (#657).

**What the five tables say today.** The notation of the cells is ``R(T)`` for ``Realized(faithful=True)``, ``R(F)`` for
``Realized(faithful=False)``, ``A(#n)`` for ``Approximated`` citing issue ``n`` and ``F(axis, #n)`` for ``Refused``; the cells below
were read from the live tables by one probe at the carve tip (:ref:`bc-aor-admission-probe`) and agree with the pin
``tests/gates/transport/test_admission_matrix.py`` where the pin has the cell (it has no "white 1.0" cell of SN or diffusion and no deck row, which this probe alone reads) (111 rows and 13 further test items, 124 passed). The value classes are those of the
pin: ``vacuum`` (the empty network, ``specular(0.0)`` and ``white(0.0)`` included), ``specular 1.0`` and ``specular 0.5``,
``white 0.5`` and ``white 1.0``, ``both 0.3 0.5`` (``from_amplitudes(0.3, 0.5)``), ``inflow 1.0`` (``uniform(1.0)``) and
``specular 0.5 over inflow 1.0`` (``specular(0.5) + uniform(1.0)``). The sites are the three geometric classes of
:class:`~orpheus.transport.boundary_conditions.FaceSite` (``CART``, ``RAD_OUTER``, ``RAD_INNER``, derived from the face name and
the coordinate system and never stored).

.. list-table:: SN and diffusion (Cartesian, radial outer and hollow inner faces)
   :header-rows: 1
   :widths: 30 14 14 14 14 14

   * - value class
     - SN ``CART``
     - SN ``RAD_OUTER``
     - SN ``RAD_INNER``
     - diffusion ``CART`` and ``RAD_OUTER``
     - diffusion ``RAD_INNER``
   * - vacuum
     - R(T)
     - R(T)
     - F(PLUMBING, #511)
     - R(T)
     - F(PLUMBING, #511)
   * - specular 1.0
     - R(T)
     - R(T)
     - R(T)
     - R(F)
     - R(F)
   * - specular 0.5
     - R(T)
     - R(T)
     - F(PLUMBING, #511)
     - R(F)
     - F(PLUMBING, #511)
   * - white 0.5, white 1.0
     - R(T)
     - F(PLUMBING, #653)
     - F(PLUMBING, #511)
     - R(F)
     - F(PLUMBING, #511)
   * - both 0.3 0.5
     - R(T)
     - F(PLUMBING, #653)
     - F(PLUMBING, #511)
     - R(F)
     - F(PLUMBING, #511)
   * - inflow 1.0, specular 0.5 over inflow 1.0
     - R(T)
     - F(PLUMBING, #653)
     - F(PLUMBING, #511)
     - F(PLUMBING, #656)
     - F(PLUMBING, #511)
   * - ``ZeroFluxCondition``
     - F(STATE_CONE)
     - F(STATE_CONE)
     - F(STATE_CONE)
     - R(T)
     - F(PLUMBING, #511)

.. list-table:: CP, MoC and MC
   :header-rows: 1
   :widths: 25 15 15 15 15 15

   * - value class
     - CP ``CART`` ``xmin``
     - CP ``CART`` ``xmax``, and ``RAD_OUTER``
     - MoC ``RAD_OUTER``
     - MC ``CART``
     - MC ``RAD_OUTER``
   * - vacuum
     - F(PLUMBING, #513)
     - R(T)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
   * - specular 1.0
     - F(PLUMBING, #513)
     - F(PLUMBING, #654)
     - A(#547)
     - F(PLUMBING, #654)
     - F(TOPOLOGY, #547)
   * - specular 0.5
     - F(PLUMBING, #513)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(TOPOLOGY, #547)
   * - white 0.5
     - F(PLUMBING, #513)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
   * - white 1.0
     - F(PLUMBING, #513)
     - R(T)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
   * - both 0.3 0.5
     - F(PLUMBING, #513)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(TOPOLOGY, #547)
   * - inflow 1.0
     - F(PLUMBING, #513)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
   * - specular 0.5 over inflow 1.0
     - F(PLUMBING, #513)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(PLUMBING, #654)
     - F(TOPOLOGY, #547)

.. list-table:: The decks, and the method-only unit, per table
   :header-rows: 1
   :widths: 18 82

   * - declaration
     - admission
   * - ``MirrorDeck``
     - SN: R(T) at a Cartesian face, certified by ``ordinate_pairing`` (else F(ANGULAR_RESOLUTION)). Diffusion: R(F) (a mirror plane
       is :math:`\mathcal{A} = 1`, tier 2). CP: R(T) at a slab's ``xmin`` only, where its kernel's image term is a forced mirror, and
       F(PLUMBING, #654) elsewhere. MC: F(PLUMBING, #654). MoC: no deck realizer (the unregistered door, #654)
   * - ``TranslationDeck``
     - SN: R(T), earned by the crossing check (a translation moves no direction and the outward normals of the two faces are
       opposite), else F(TOPOLOGY) (a pair whose crossing sets differ is not a translation quotient; the realizer's docstring names #178 for
       the explicit gluing map it would need). Diffusion: F(TOPOLOGY) (a per-face scalar has no slot for cross-face coupling). CP:
       F(PLUMBING, #654), and F(PLUMBING, #513) at a face CP drops. MC: R(T) (the pitch is not compared with the declared extent,
       #655). MoC: the unregistered door, #654
   * - ``ZeroFluxCondition``
     - diffusion: R(T), :math:`\mathcal{A} = -1`. Every other table: F(STATE_CONE) through the unregistered door

This page states which declarations each table admits and why; what the realized operators do is the subject of the method pages:
SN's :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`, ``ScaledOrdinatePermutation`` and ``RankOneReemission`` as terms of
the boundary operator :math:`B` are in :doc:`/theory/methods/sn/boundary_conditions`, and diffusion's albedo operator on the
scalar trace is in :doc:`/theory/methods/diffusion_1d`.

The tables are the data of the following rulings and measurements, each of which is also a reason a cell is what it is.

- **SN refuses on the radial outer face what its corner action cannot carry** (#653). The radial characteristic's corner action is
  ruled only for a permuting return or a zero return, and a source is a free parameter that the corner block (a linear operator on
  the ray alone) cannot carry. Before the carve the solve constructed and refused at the first sweep (12 of 12 solves, `[M]` the
  WF2.1 measurement of plan section 4.3, the memo ``measure_admission_table.md``, untracked scratch and not re-run); the refusal is
  now at admission. The zero network and the specular chain are realized there
  (a deck on a radial face is not spellable).
- **The hollow inner face has a refusal of its own** (#511). The radial axis carries one trace slot, the outer surface, so a method
  assembled on the axes' trace inventory has nowhere to put an inner condition: it computes the reflective cavity (a unit specular
  return) whatever is declared. SN and diffusion therefore refuse every inner declaration but ``specular(1.0)`` through the one
  :func:`~orpheus.transport.boundary_conditions.hollow_inner_refusal`, with the method's reason. `[M]` 2026-10-23, the S3 build's
  probe, recorded in the docstring of ``tests/gates/diffusion/test_diffusion_hollow_inner_law.py`` and not re-run for this page: before this admission an inner ``vacuum()``, ``specular(1.0)`` and ``white(0.5)`` gave a bit-identical eigenvalue on a hollow
  sphere and on a hollow cylinder in diffusion (3 of 3 conditions on 2 of 2 coordinate systems); SN had dropped it until the P1 guard
  of 2026-09-25. The admission is per table because the reason is per method: CP refuses a hollow body's inner face on its own guard
  (#513: it reads only the outer face), MoC and MC refuse a hollow body outright (#514, #513), and a method that gains an
  inner-surface trace lifts its own refusal alone. The gates are ``tests/gates/sn/architecture/test_sn_hollow_inner_law.py``, ``tests/gates/diffusion/test_diffusion_hollow_inner_law.py`` and, for CP,
  ``tests/gates/cp/test_cp_realizer_table.py``.
- **Today's limits stand (ruling D3).** The focus of the carve was reorganisation, not capability. A limit is a table cell citing
  its issue; widening it is that issue's work (#650 for the reference's vocabulary, #653, #654, #656 for the tables).
- **No masking (ruling D4).** CP realizes the vacuum as what it computes. Its eigenvalue omits the leakage term of the balance, so a
  vacuum-bounded eigenvalue exceeds a white-bounded one on the same cell (#546); the admission is ``Realized(True)`` and the defect
  is kept visible by three strict ``xfail`` rows asserting the physics (``k_vacuum <= k_white`` on the two-group ``A|B`` body, on
  3 of 3 charts) with ``raises=AssertionError`` (a bare strict ``xfail`` would also swallow an ``ImportError``), which turn red
  (XPASS) the day #546 is repaired. `[M]` the figures recorded in the docstring of
  ``tests/gates/cp/test_cp_vacuum_not_above_white.py`` (``_probe_546_xfail.py`` at ``5d04c296``, ``CPParams(keff_tol=1e-11)``, not
  re-run for this page): slab 1.4003162 against 1.2308420, cylinder 0.9354010 against 0.7393875, sphere
  0.5854991 against 0.4140867. The CP table also refuses an albedo below one on a diffuse chain (F(PLUMBING, #654)): the retired
  string registry keyed its closure by the word ``white`` whatever the albedo, a silent drop that is now a refusal.
- **MC computes a torus, not the declared cylinder** (#547). MC samples a periodic equal-area square lattice cell, a translation
  whose opposite faces are coupled by an image lattice, and the declared cylinder has no such coupling: its specular wall on a cylinder
  is ``Refused(TOPOLOGY, #547)`` for every weight, and the two tests that declare it are strict ``xfail`` rows
  (``tests/gates/mc/test_gaps.py``, ``tests/gates/mc/test_properties.py``) that turn red when MC takes a 2-D geometry value that
  describes what it computes. The axis is topology and not plumbing because theory stands behind the refusal: the carrier MC answers
  with and the declared body differ in exactly the question the topological axis asks. MoC realizes the same declaration as a named
  domain approximation, ``A(#547)`` (a specular return on the equal-area square cell), for weight one only.
- **The mixed cells follow the meet.** MC at ``RAD_OUTER``, "both chains" and "specular over a source", read
  ``F(TOPOLOGY, #547)``, because the meet keeps the first refusal in ``CHAIN_ORDER`` (specular before diffuse, chains before the
  source).

.. _bc-aor-admission-probe:

The probe behind the tables
~~~~~~~~~~~~~~~~~~~~~~~~~~~

One probe produced the cells above. It read each of the five sealed tables (``SN_REALIZERS``, ``DIFFUSION_REALIZERS``,
``CP_REALIZERS``, ``MOC_REALIZERS``, ``MC_REALIZERS``) through ``admit`` at the sites above, with the value classes above, and
printed one line per value class. `[M]` ``python -O`` run from outside the repository with ``PYTHONPATH`` set to the worktree and
``orpheus.__file__`` printed first; the ``catalogue`` of each table printed beside it read: SN, all four keys R(T);
diffusion: ``ZeroNetwork`` R(T), ``SpecularChain`` and ``DiffuseChain`` R(F), ``ZeroFluxCondition`` R(T), ``InflowTerms``
F(PLUMBING, #656); CP: ``ZeroNetwork`` R(T), ``SpecularChain`` F(PLUMBING, #654), ``DiffuseChain`` R(T); MoC: ``ZeroNetwork``
F(PLUMBING, #654), ``SpecularChain`` A(#547), ``DiffuseChain`` F(PLUMBING, #654); MC: ``ZeroNetwork`` and ``DiffuseChain``
F(PLUMBING, #654), ``SpecularChain`` F(TOPOLOGY, #547). The catalogue is the upper bound of the instance level: no sampled cell above exceeds its key's entry
(law L-D, ``test_catalogue_is_an_upper_bound`` samples the same property over a test-local table).

.. _bc-aor-adjoint:

The adjoint, stated once
------------------------

The adjoint entries of the SN solvers dagger the same objects the forward solve uses, so a forward and an adjoint solve share one
binding. Four statements make that true, and the fourth is the one a consumer reads.

#. **Every shared kernel is reciprocal in the trace metric** :math:`M = \mathrm{diag}(w|\mu|)`. The specular return is an
   involutive measure-preserving permutation :math:`P`, so :math:`M^{-1}P^{\mathsf T}M = P`; the Lambertian is
   :math:`u \otimes (w|\mu|)/M_-` with :math:`M_- = \sum_{\Omega\cdot n<0} w|\mu|`, self-adjoint under :math:`M`. The adjoint of a shared condition's realization is
   therefore its realization on the reversed ordinates (``test_reciprocity`` in
   ``tests/gates/sn/architecture/test_realized_boundary.py``).
#. **A deck's adjoint is its partner's deck**: :math:`G_g^\sharp = G_{g^{-1}}`, which is the deck on the partner face,
   ``pairing_of(d, f).inverse == pairing_of(d, partner(f))``; it holds for a translation because ``translation_motion`` is signed by
   the face (``test_pairing_of_and_its_adjoint`` in ``tests/gates/geometry/test_deck.py``; the code block of
   :ref:`bc-aor-deck-family` prints the equality).
#. **The inflow** :math:`q` **is not adjointed.** ``solve_sn_adjoint_fixed_source`` builds a zero ``AngularBoundarySourceSink`` and
   never reads the declared :math:`q`, and the two eigenvalue entries, ``solve_sn`` and ``solve_sn_adjoint``, refuse a declared boundary source
   (``_refuse_a_boundary_source``, ``ERR-095``). This is correct: :math:`q` enters the reciprocity identity as the inner product
   :math:`\langle q, |\Omega\cdot n|\,\gamma_-\psi^\dagger\rangle` (`[R]` the identity is the plan's statement, not a gate's: the
   gate pins the behaviour, not the identity), a term of the adjoint problem's functional and not of its operator. The statement
   had never been written down; ``test_the_adjoint_ignores_q`` in
   ``tests/gates/sn/architecture/test_sn_boundary_binding.py`` pins it (the same adjoint flux for a vacuum face, ``uniform(1.0)`` and
   ``uniform(50.0)``, while ``solve_sn`` and ``solve_sn_adjoint`` refuse a bound source).
#. **The left-eigenvector rider.** Reciprocity in :math:`M` is a statement about :math:`B^{\mathsf T}`. The scalar :math:`\rho` a
   consumer reads exists exactly when the trace measure is a left eigenvector of the realized map: :math:`J^-(R\varphi) = \rho
   J^+(\varphi)` for every :math:`\varphi \ge 0` if and only if :math:`R^{\mathsf T}(w|\mu|)_- = \rho\,(w|\mu|)_+`. It holds by
   construction for the shared family (:ref:`bc-aor-wall-condition`); a method-only kernel owes it as a gate.

The eigenvalue entries (:func:`~orpheus.sn.solver.solve_sn` and ``solve_sn_adjoint``) refuse a declared boundary source by posing:
an eigenproblem is homogeneous, so a source has no term to land in, and answering the source-free problem would drop the declaration
in silence.

.. _bc-aor-inflow:

The inflow channel, in brief
----------------------------

The affine part :math:`q \ge 0` is a function on the inflow half-trace :math:`\Gamma_-` and the one factor the realizer tier does not
carry. :class:`~orpheus.transport.boundary_conditions.InflowTerms` states it as a weighted sum of profiles: ``uniform`` (every uniform
contribution folded into one number, so ``uniform(1.0) + uniform(2.0) == uniform(3.0)`` with no profile object for a constant) plus
``sampled``, a tuple of ``(weight, profile)`` in first-appearance order, whose parse drops a zero weight, merges equal profiles and
refuses anything that is not an inflow profile. **The frame is the type and the content is open**: the profiles are Python callables
over a space, not closed expressions, so what is decidable, and what the binding's invariance refusal needs, is the frame the profile
is written in. :class:`~orpheus.transport.boundary_conditions.FaceFrameInflow` is written in the face's own frame and is face-local;
:class:`~orpheus.transport.boundary_conditions.LabFrameInflow` is written in laboratory directions and is not. A profile
implements the open protocol :class:`~orpheus.transport.boundary_conditions.InflowSampler` (``sample(space)`` returns an array with
exactly ``space.shape``, in the units of :math:`\psi`).

**Sampling is the method's, not the declaration's.** The single verb that evaluates :math:`q` on a method's trace is
:meth:`~orpheus.transport.boundary_conditions.RealizerTable.sample_inflow`: ``require`` the source (a source the method cannot carry
is refused at the door, so the refusal and the sampling cannot disagree), then the table's registered
:class:`~orpheus.transport.boundary_conditions.SourceRealizer`, whose default body is
:func:`~orpheus.transport.boundary_conditions.sample_terms`: ``q = uniform * ones + sum_i w_i sampler_i.sample(space)``, with a
``ValueError`` when a sampler returns anything but ``space.shape`` (a wrong-shape array would otherwise be broadcast silently). SN
overrides ``sample`` to run the placement check of ``ERR-047``: the space sampled on must be the inflow half-trace of the face, since
a source written anywhere else is discarded by the sweep. The sampled arrays are packed into the boundary leaf of the composite
source :math:`q = q_{\rm bulk} \oplus q_\partial`,
:class:`~orpheus.transport.fields.AngularBoundarySourceSink` (:meth:`~orpheus.transport.fields.AngularBoundarySourceSink.from_method`
reads the faces whose bound condition carries a nonzero source), and
:func:`~orpheus.sn.solver._build_fixed_source_rhs` is the single construction point of the right-hand side for both inner solvers.

Two consequences follow, and both are derived and measured in the older sections. A prescribed inflow realizes to the **zero map**:
the table consults the ``ZeroNetwork`` realizer for an empty network whatever the source, so vacuum and a prescribed inflow are the
same operator and differ only in :math:`q` (:ref:`bc-affine-source-channel`). And :math:`q` never rides an operator slot, because an
operator that ignores its input and returns :math:`q` is an affine map in a linear slot, which doubles the delivery under source
iteration and raises a convergence-claim error under Krylov (:ref:`bc-affine-channel-the-defect`). The sites where the source is
refused are table cells: diffusion has no channel for it (#656, plumbing), and SN refuses it on a radial outer face (#653) and on a
hollow inner face (#511).

.. _bc-aor-laws:

The laws and the gates that pin them
------------------------------------

Each law below is a statement the types obey, with the gate that fails first when it is false. A law is a *claim* and the gate is its
*instrument*; the gates were re-run at the carve tip in the two runs of the V&V admonition (494 passed and 3 xfailed over its 17 files; 785 passed over the eight further law-gate files it names).

.. list-table:: The laws of the architecture
   :header-rows: 1
   :widths: 14 46 40

   * - law
     - statement
     - gate
   * - L-A the adjoint
     - reciprocity in the trace metric; a deck's adjoint is its partner's; :math:`q` is not adjointed; the left-eigenvector rider
     - ``test_reciprocity``, ``test_pairing_of_and_its_adjoint``, ``test_the_adjoint_ignores_q``,
       ``test_the_left_eigenvector_rider_through_the_realizer`` (the identity on registry rules), and
       ``test_the_rider_on_a_trace_whose_masses_differ`` (the inflow-mass divisor)
   * - L-B the cone
     - construction outside :math:`K` raises; arithmetic returns a condition if and only if the pair lies in :math:`K`; the albedo
       twin is gone
     - ``test_wall_condition_parse``, ``test_the_algebra_and_the_cone``
   * - L-C invariance
     - a condition a non-identity element of :math:`H` moves is refused at ``bind``
     - ``test_bind_invariance``
   * - L-D admission equals realization
     - ``admit`` is ``Refused`` if and only if ``realize`` raises ``BoundaryRefused``; ``catalogue`` is an upper bound
     - ``test_one_door``, ``test_catalogue_is_an_upper_bound``
   * - L-E naturality of the mirror
     - the slab mirror :math:`\sigma` carries the realization at ``xmin`` onto that at ``xmax``,
       :math:`P_\sigma\,\mathrm{real}(\beta, \texttt{xmin})\,P_\sigma^{-1} = \mathrm{real}(\beta, \texttt{xmax})`, for every shared
       unit: two realizations compared through a third object, so it does not read ``face_normal``
     - ``test_mirror_naturality``
   * - L-F a deck is not an operand, a condition is not a deck
     - ``0.5 * MirrorDeck()`` and ``vacuum() + TranslationDeck()`` raise Python's own ``TypeError``; no value puts a response in the
       deck slot or a deck in the condition slot
     - ``test_a_deck_is_not_an_operand``, ``test_a_condition_is_not_a_deck_slot``
   * - L-G one definition per quantity
     - the coordinate-agreement convention has one definition; ``white(a)`` and ``from_amplitudes(0.0, a)`` have one digest; the
       wall-bound message has one owner in transport and one in the references (the insulation ruling: references import no
       production)
     - ``test_white_is_a_facade``, ``test_one_agreement_convention``
   * - L-H layering
     - ``geometry`` imports nothing from numerics, mesh, specification, reference, transport or the L3 packages
       (``FORBIDDEN_EDGES["geometry"]``)
     - ``tests/gates/test_layer_imports.py``
   * - L-I no string dispatch in the core
     - no ``dict[str, ...]`` registry on a method class, no ``.kind`` read on a declaration, ``BC`` does not exist
     - ``test_the_string_registries_and_the_typed_law_reads_are_gone`` in ``tests/gates/transport/test_wave2_s4_tables.py``
   * - L-J identity
     - equal declarations digest equal; every deck and condition type has its row in the content-identity rosters
     - ``tests/gates/transport/test_content_identity_transport.py``, ``tests/gates/geometry/test_content_identity_geometry.py``
   * - L-K the carried numbers
     - every eigenvalue and flux the suite pins is unchanged by the carve except the declared exceptions
     - ``tests/gates/sn/architecture/test_wave2_carry.py``
   * - L-L the functor law
     - realization is a homomorphism of the algebra: exact on parameters, within a few ulp on dense matrices
     - ``test_the_functor_law``, ``test_the_realizer_is_a_homomorphism``
   * - L-M chains compose, the family is closed
     - every chain's factor tiers compose from outflow to inflow; every ``match`` over ``Chain`` ends in ``assert_never``
     - ``test_chains_compose_and_the_family_is_closed``

.. _bc-wave2-crosswalk:

From the trace-law architecture to the architecture of record
=============================================================

Until 2026-10-24 a boundary condition was a *trace law*: a ``BoundaryTraceLaw`` subclass chosen by a string key, split into a
geometry factor :math:`G` and a response factor :math:`R`, composed by a descriptor algebra (``LawSum``, ``LawScaled``), realized by a
``BoundaryRealizer`` that a string registry reached, and declared on the geometry through the tag ``BC``. Wave 2 of the layer DAG
replaced it, in five slices of one carve and a retirement, by the architecture of :ref:`bc-architecture-of-record`: face roles that
are a deck or an integer tag, one affine wall condition on the sub-Markov simplex, a binding with the colouring stabiliser, and one realizer table
per method that answers admission and realization together. The record of what landed is the development-history row at the end of
the page.

**How to read the rest of the page.** The sections below this table record the architecture that wave 2 replaced. Their mathematics
still holds where they re-spell it (the two-test criterion of :ref:`bc-factor-roles`, the quotient picture of
:ref:`bc-factor-quotients`, the three tiers of :ref:`bc-method-realizability`, the trace maps of :ref:`bc-domain-narrowing`, the
affine source channel of :ref:`bc-affine-source-channel`), and a name in double backticks there is retired. A section wholly about
retired machinery says so in its first sentence and points back here. The table maps every retired symbol to its successor as a live
role, or to ``none`` with the reason (a successor that is a property of a type, such as the trace ladder, has no symbol of its own).
The steps ``R1`` to ``R20`` are the rows of the retirement's dependency table, in the order of deletion (leaves first); the row
``R-ERR041`` is the one the retirement added.

.. list-table:: Every retired symbol and its successor
   :header-rows: 1
   :widths: 8 34 58

   * - step
     - retired (``literal``)
     - successor (a live role), or ``none`` and the reason
   * - R1
     - ``assert_specular_pairing_valid``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (the composition of the three checks)
   * - R1
     - ``assert_specular_pairing_involutive``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (its involution check,
       ``PairingNotInvolutiveError``, ``ERR-044``)
   * - R1
     - ``assert_specular_pairing_measure_preserving``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (its measure check,
       ``PairingNotMeasurePreservingError``, ``ERR-042``)
   * - R1
     - ``assert_specular_pairing_maps_inflow_to_outflow``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (its crossing check,
       ``PairingDoesNotCrossError``, ``ERR-045``)
   * - R1
     - ``_specular_pairing``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (it returns the permutation the realizer
       uses)
   * - R1
     - ``_axis_cosines``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` reads ``Quadrature.mean_axis_cosine``:
       no successor symbol
   * - R2
     - ``IncomingOutgoingTraceClassificationError``
     - none: the classification is carried by the trace ladder,
       :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`, whose half-spaces exclude the
       tangential ordinates from both halves (``ERR-040``)
   * - R2
     - ``VacuumAppliedToOutgoingTraceError``
     - none: the zero map reads both half-traces of its own face off the trace ladder, so the mismatch is unspellable
       (``ERR-041``, :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`)
   * - R2
     - ``BoundaryGeometryMapNotMeasurePreservingError``
     - :class:`~orpheus.numerics.quadrature.ordinate_pairing.PairingNotMeasurePreservingError` (``ERR-042``)
   * - R2
     - ``ReflectionNotInvolutiveError``
     - :class:`~orpheus.numerics.quadrature.ordinate_pairing.PairingNotInvolutiveError` (``ERR-044``)
   * - R2
     - ``ReflectionDidNotMapInflowToOutflowError``
     - :class:`~orpheus.numerics.quadrature.ordinate_pairing.PairingDoesNotCrossError` (``ERR-045``)
   * - R2
     - ``BoundarySourceNotOnIncomingTraceError (the old home)``
     - :class:`~orpheus.sn.boundary.errors.BoundarySourceNotOnIncomingTraceError` (SN's own, rooted at the transport
       :class:`~orpheus.transport.boundary_conditions.BoundaryError`; ``ERR-047``)
   * - R3
     - ``TAG_REGISTRY``
     - :meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.of` (the exterior of a reference is plain
       data, ``Walls.of(geometry, exterior)``)
   * - R3
     - ``_law_of``
     - :meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.of`
   * - R3
     - ``_wall_of``
     - :meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.of`
   * - R3
     - ``_parameters``
     - :meth:`~orpheus.derivations.continuous.characteristic.walls.Walls.of`
   * - R3
     - ``specular_albedo (the law form)``
     - :func:`~orpheus.derivations.common.reference_body.specular_albedos` (``specular_albedos(geometry, exterior)``)
   * - R3
     - ``BC.to_alpha``
     - :func:`~orpheus.derivations.common.reference_body.specular_albedos` (the amplitudes are read from
       ``WallAmplitudes``)
   * - R4
     - ``BC_REGISTRY (CP)``
     - ``CP_REALIZERS`` (``orpheus/cp/solver.py``)
   * - R4
     - ``bc_kind (CP)``
     - ``CP_REALIZERS`` (``orpheus/cp/solver.py``)
   * - R4
     - ``BC_REGISTRY (MoC)``
     - ``MOC_REALIZERS`` (``orpheus/moc/geometry.py``)
   * - R4
     - ``bc_kind (MoC)``
     - ``MOC_REALIZERS`` (``orpheus/moc/geometry.py``)
   * - R4
     - ``BC_REGISTRY (MC)``
     - ``MC_REALIZERS`` (``orpheus/mc/solver.py``)
   * - R5
     - ``BOUNDARY_OPERATOR_REGISTRY (SN)``
     - ``SN_REALIZERS`` (``orpheus/sn/boundary/realizer.py``)
   * - R5
     - ``BOUNDARY_OPERATOR_REGISTRY (diffusion)``
     - ``DIFFUSION_REALIZERS`` (``orpheus/diffusion/boundary_realizer.py``)
   * - R5
     - ``_law_from_tag``
     - :func:`~orpheus.transport.boundary_conditions.bind`
   * - R5
     - ``resolve_boundary_conditions``
     - :func:`~orpheus.transport.method.realize_boundary` after :func:`~orpheus.transport.boundary_conditions.bind`
   * - R5
     - ``realize_boundary_law``
     - :meth:`~orpheus.transport.method.TransportMethod.realize_boundary_role`
   * - R6
     - ``assert_realizable``
     - :meth:`~orpheus.transport.boundary_conditions.RealizerTable.admit` and
       :meth:`~orpheus.transport.boundary_conditions.RealizerTable.require` (admission is the table's)
   * - R6
     - ``assert_inflow_outflow_classification``
     - the trace ladder, :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` (``ERR-040``)
   * - R6
     - ``assert_outgoing_leakage_unconstrained``
     - none: a realized response is typed :math:`\Gamma_+ \to \Gamma_-`, so no outflow row is in its codomain
   * - R6
     - ``assert_geometry_map_measure_preserving``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (``ERR-042``)
   * - R6
     - ``assert_is_involutive``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (``ERR-044``)
   * - R6
     - ``assert_reflection_maps_inflow_to_outflow``
     - :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` (``ERR-045``) and the crossing check of
       :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`
   * - R6
     - ``assert_response_positive_if_declared``
     - :class:`~orpheus.transport.boundary_conditions.BoundaryResponseNotPositiveError`, raised by the cone parse of
       :class:`~orpheus.transport.boundary_conditions.WallCondition` (``ERR-043``)
   * - R6
     - ``assert_submarkov``
     - :class:`~orpheus.transport.boundary_conditions.SubmarkovViolationError`, raised by the cone parse of
       :class:`~orpheus.transport.boundary_conditions.WallCondition` (``ERR-046``)
   * - R6
     - ``assert_source_is_placeable``
     - :meth:`~orpheus.sn.boundary.realizer.SNSource.sample` (the placement check, ``ERR-047``)
   * - R6
     - ``law_permutes_ordinates``
     - :attr:`~orpheus.sn.boundary.realized_operators.RealizedBoundary.permutes_ordinates`
   * - R6
     - ``geometry_map``
     - :class:`~orpheus.geometry.deck.FacePairing`, read through :func:`~orpheus.geometry.deck.pairing_of`
   * - R6
     - ``response_kernel``
     - :class:`~orpheus.transport.boundary_conditions.ResponseNetwork` and
       :attr:`~orpheus.transport.boundary_conditions.WallCondition.return_fraction`
   * - R7
     - ``_BoundBoundaryOperator``
     - :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`
   * - R8
     - ``realize_recursively``
     - :meth:`~orpheus.transport.boundary_conditions.RealizerTable.realize` (realization is a functor on the response
       network)
   * - R8
     - ``BoundaryRealizer (the Protocol)``
     - :class:`~orpheus.transport.boundary_conditions.RealizerTable`
   * - R8
     - ``SNBoundaryRealizer``
     - ``SN_REALIZERS``
   * - R8
     - ``DiffusionBoundaryRealizer``
     - ``DIFFUSION_REALIZERS``
   * - R8
     - ``BoundaryRealizerRegistry (dissolved earlier, #290)``
     - none: every method holds its own sealed table
   * - R9
     - ``LawSum``
     - :meth:`~orpheus.transport.boundary_conditions.WallCondition.__add__` (a sum of conditions);
       :class:`~orpheus.transport.boundary_conditions.WallDisplacement` outside the cone
   * - R9
     - ``LawScaled``
     - :meth:`~orpheus.transport.boundary_conditions.WallCondition.__mul__` (a scalar multiple);
       :class:`~orpheus.transport.boundary_conditions.WallDisplacement` outside the cone
   * - R9
     - ``LawNode``
     - none: a condition or a displacement is the node
   * - R9
     - ``_admit_as_operand``
     - none: a deck defines no operator, so Python's own ``TypeError`` refuses ``0.5 * MirrorDeck()`` (law L-F)
   * - R10
     - ``ScalarResponse``
     - a weight of :class:`~orpheus.transport.boundary_conditions.ResponseNetwork`
   * - R10
     - ``LambertianReemission``
     - :class:`~orpheus.transport.boundary_conditions.DiffuseChain`
   * - R10
     - ``SpecularReemission``
     - :class:`~orpheus.transport.boundary_conditions.SpecularChain`
   * - R10
     - ``ReemissionClosure``
     - the pair of chains, :class:`~orpheus.transport.boundary_conditions.SpecularChain` and
       :class:`~orpheus.transport.boundary_conditions.DiffuseChain`
   * - R10
     - ``BoundaryResponseKernel``
     - :class:`~orpheus.transport.boundary_conditions.ResponseNetwork` and
       :attr:`~orpheus.transport.boundary_conditions.WallCondition.return_fraction`
   * - R11
     - ``NoSource``
     - :meth:`~orpheus.transport.boundary_conditions.InflowTerms.empty`
   * - R11
     - ``ConstantInflowSource``
     - :class:`~orpheus.transport.boundary_conditions.InflowTerms` (the ``uniform`` coordinate), spelled
       ``WallCondition.uniform``
   * - R11
     - ``InflowSourceSpec``
     - :class:`~orpheus.transport.boundary_conditions.InflowSampler` (``sample``), evaluated by
       :meth:`~orpheus.transport.boundary_conditions.RealizerTable.sample_inflow`
   * - R12
     - ``SelfPairedDeck``
     - :class:`~orpheus.geometry.deck.MirrorDeck`
   * - R12
     - ``PairedDeck``
     - :class:`~orpheus.geometry.deck.TranslationDeck`
   * - R12
     - ``BoundaryGeometryMap``
     - :class:`~orpheus.geometry.deck.FacePairing`, read through :func:`~orpheus.geometry.deck.pairing_of`
   * - R12
     - ``_mirror_motion``
     - :func:`~orpheus.geometry.deck.mirror_motion`
   * - R13
     - ``WhiteBoundary``
     - :meth:`~orpheus.transport.boundary_conditions.WallCondition.white`
   * - R13
     - ``ReflectiveBoundary``
     - :class:`~orpheus.geometry.deck.MirrorDeck` on a Cartesian face; on a curved face a tag bound to
       :meth:`~orpheus.transport.boundary_conditions.WallCondition.specular` at 1.0
   * - R13
     - ``PeriodicBoundary``
     - :class:`~orpheus.geometry.deck.TranslationDeck`
   * - R13
     - ``VacuumInflow``
     - :meth:`~orpheus.transport.boundary_conditions.WallCondition.vacuum`
   * - R13
     - ``PrescribedInflow``
     - :meth:`~orpheus.transport.boundary_conditions.WallCondition.uniform` and
       :meth:`~orpheus.transport.boundary_conditions.WallCondition.framed`
   * - R13
     - ``AlbedoBoundary (the old type)``
     - :meth:`~orpheus.transport.boundary_conditions.WallCondition.specular` or
       :meth:`~orpheus.transport.boundary_conditions.WallCondition.white` at the albedo; ``albedo = 0`` is the vacuum
   * - R13
     - ``SpecularReturn (the old type)``
     - :class:`~orpheus.transport.boundary_conditions.SpecularChain`
   * - R13
     - ``IsotropicReturn``
     - :class:`~orpheus.transport.boundary_conditions.DiffuseChain`
   * - R13
     - ``ZeroFluxBoundary (the old type)``
     - :class:`~orpheus.diffusion.zero_flux.ZeroFluxCondition`
   * - R13
     - ``BoundaryTraceLaw``
     - :class:`~orpheus.transport.boundary_conditions.BoundaryCondition` for a response and the ``FaceRole`` sum for a
       deck
   * - R13
     - RegistryMixin ``key=``, ``kind`` and ``.registry`` on the seven laws
     - none: dispatch is by type, through the realizer table
   * - R14
     - ``parse_boundary_law``
     - :func:`~orpheus.geometry.deck.parse_face_role`
   * - R14
     - ``parse_boundary_laws``
     - :func:`~orpheus.geometry.deck.parse_face_roles`
   * - R14
     - ``FaceLaws``
     - :class:`~orpheus.mesh.face_roles.FaceRoles`
   * - R14
     - ``face_laws``
     - ``face_roles`` (of a mesh, :class:`~orpheus.mesh.face_roles.FaceRoles`; of a geometry,
       ``StructuredGeometry.face_roles``)
   * - R14
     - ``boundaries (of a geometry)``
     - :attr:`~orpheus.geometry.structured_geometry.StructuredGeometry.face_roles`
   * - R14
     - ``bc_low``
     - :attr:`~orpheus.mesh.axis.AxisMesh.role_low`
   * - R14
     - ``bc_high``
     - :attr:`~orpheus.mesh.axis.AxisMesh.role_high`
   * - R14
     - ``bc_outer (of a radial axis)``
     - :attr:`~orpheus.mesh.axis.RadialAxisMesh.role_outer`
   * - R14
     - ``outer_law (of a mesh)``
     - :attr:`~orpheus.mesh.structured.Mesh1D.outer_role`
   * - R15
     - ``BoundaryError (the old home)``
     - :class:`~orpheus.transport.boundary_conditions.BoundaryError` (carries ``.declaration``, the object, where the
       old one carried a string ``.law``)
   * - R15
     - ``BoundaryResponseNotPositiveError (the old home)``
     - :class:`~orpheus.transport.boundary_conditions.BoundaryResponseNotPositiveError`
   * - R15
     - ``SubmarkovViolationError (the old home)``
     - :class:`~orpheus.transport.boundary_conditions.SubmarkovViolationError`
   * - R16
     - ``BC``
     - integer tags (``FaceRole``) and :class:`~orpheus.transport.boundary_conditions.BoundaryConditions`
   * - R17
     - the package ``orpheus/geometry/boundary/`` (15 modules)
     - ``orpheus/geometry/deck.py`` and ``orpheus/transport/boundary_conditions/``
   * - R17
     - ``orpheus/transport/boundary_realization.py``
     - ``orpheus/transport/boundary_conditions/realization.py``
       (:func:`~orpheus.transport.boundary_conditions.stamp_boundary_role`)
   * - R17
     - ``orpheus/sn/boundary/bound_compat.py``
     - ``orpheus/sn/boundary/realized_operators.py``
   * - R18
     - the ``Literal`` alias ``BoundaryCondition`` of ``continuous_reference``
     - ``ReferenceBoundaryKind`` (the reference's own string vocabulary)
   * - R18
     - ``ProblemSpec.boundary_conditions``
     - :attr:`~orpheus.derivations.common.continuous_reference.ProblemSpec.boundary_kinds` (written by 19 cases and
       read by none: #660)
   * - R19
     - the string arm of ``__eq__`` (``VacuumInflow() == "vacuum"``, #553)
     - none: equality is content identity, and the string spelling retired with ``kind``
   * - R20
     - ``_narrowed_zero_operator``
     - :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`
   * - R20
     - ``_deck_kernel``
     - :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`
   * - R20
     - ``_specular_crossed_diagnosis``
     - none: the crossing check of ``ScaledOrdinatePermutation`` at construction (``ERR-045``)
   * - R20
     - ``_checked_angular_average``
     - :class:`~orpheus.sn.boundary.realized_operators.RankOneReemission`
   * - R20
     - ``_attenuated_kernel_operator (the realizer's)``
     - moved to ``realized_operators.py`` as the one private helper
   * - R20
     - ``_outflow_restriction``
     - the trace ladder's ``trace.outflow_space(face)``
   * - R20
     - ``_partner_outflow_restriction``
     - the trace ladder's ``trace.outflow_space(partner)``
   * - R-ERR041
     - ``SNMethodSpace.inflow_indices and SNMethodSpace.outflow_indices``
     - the trace ladder: ``trace.inflow_space(face)``, ``trace.outflow_space(face)``
   * - R-ERR041
     - ``SNMethodSpace.face, for_face, minimal and inflow_indices_for_face``
     - :class:`~orpheus.sn.mesh.method_space.SNMethodSpace` is now ``(quadrature, trace)``; the face travels in the
       :class:`~orpheus.transport.boundary_conditions.FaceSite`


.. _bc-overview-three-layers:

The §16A three-layer decomposition
==================================

A boundary condition in transport-theory codes is, in the
discrete-form-typical mathematical sense, a **single linear operator**
that takes the outgoing :term:`angular flux` at a face and returns the
incoming angular flux. In ORPHEUS we explicitly factor that single
operator into three layers because each layer carries different
mathematical, physical, and architectural responsibilities. The split
is taken verbatim from Grand Report v3 §16A.3 and the source-of-record
plan ``.claude/plans/transient-giggling-cake.md``.

The mathematics of this section holds; its type carriers were those of the
trace-law architecture that wave 2 retired (2026-10-24). The present
architecture is :ref:`bc-architecture-of-record`, and each retired name maps to
its successor in :ref:`bc-wave2-crosswalk`. In the retired trace-law
architecture one object, the ``BoundaryTraceLaw``, held the three factors below as three properties of
one value; wave 2 split them by the dichotomy this section proves
(:ref:`bc-factor-roles`): :math:`G` is the **face role** of a glued face, a
deck (:class:`~orpheus.geometry.deck.MirrorDeck` or
:class:`~orpheus.geometry.deck.TranslationDeck`); :math:`R` and :math:`q` are
the two halves of one :class:`~orpheus.transport.boundary_conditions.WallCondition`,
the condition a tag is bound to; and the realization of each is the method's
realizer table (:ref:`bc-realizer-layer`). A face has exactly one role, the sum
``FaceRole = DeckTransformation | int``, so a face carries a :math:`G` and no
:math:`(R, q)`, or a tag and no :math:`G`.

.. _affine-bc-form:

Layer 2 — the affine law on the boundary trace
----------------------------------------------

The full mathematical form of every boundary law in this codebase is
the **affine map**

.. math::
   :label: affine-bc-form

   \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q,

.. (vv-status rationale) Structural/definitional framing: the master affine
   form of every boundary law (Grand Report §16A.3). Per this page's own note,
   its equations are definitional / structural-architecture statements; the
   concrete rank-0 / rank-n realisations (:eq:`bc-rank-n-tensor-decomposition`,
   ``WallCondition.uniform``) are the tested forms. Not a solver claim.
.. vv-status: affine-bc-form documented

where:

* :math:`\gamma_\pm` are the **trace operators** that restrict the
  angular flux :math:`\psi(\mathbf{r}, \Omega)` from the volumetric
  function space to the inflow / outflow boundary trace spaces
  :math:`\Gamma_\pm` (see :ref:`bc-trace-structure` below for the
  formal definition).
* :math:`G` is the **deck transformation**: the composition (Koopman)
  operator of an isometry :math:`g` of the covering space, a
  measure-preserving permutation (the mirror) or spatial wrap-around (the
  translation). It carries pure geometry (it changes nothing about the
  physical interaction at the boundary; it just relabels the angular fluxes
  that meet there). When it is the non-trivial factor it is a map
  :math:`\Gamma_+(f') \to \Gamma_-(f)` (:math:`f' = f` for a mirror) and it
  carries the **crossing**, because the mirror that exchanges the two
  hemispheres is an ambient isometry — see :ref:`bc-factor-roles` below. It is
  realized by two types that have no field and no parameter, because the face
  fixes the plane, the partner and the direction of the motion
  (:class:`~orpheus.geometry.deck.MirrorDeck`,
  :class:`~orpheus.geometry.deck.TranslationDeck`, read through
  :func:`~orpheus.geometry.deck.pairing_of`).

* :math:`R` is the **response kernel** — the constitutive law, a
  :class:`~orpheus.transport.boundary_conditions.ResponseNetwork`
  :math:`R = s\,P + d\,(E \circ Rt)`: a weighted sum of the two chains
  that compose from :math:`\Gamma_+` to :math:`\Gamma_-`, the specular chain
  :math:`P` (a scalar amplitude on a mirror permutation, for a polished
  wall) and the diffuse chain :math:`E \circ Rt` (the rank-one angular kernel
  of Lambertian re-emission). The weights are a point :math:`(s, d)` of the
  sub-Markov simplex :math:`K = \{s, d \ge 0,\ s + d \le 1\}`, so the response is
  sub-Markov by construction. When it is the non-trivial factor it is the one
  that crosses, and every realized response in the tree is typed
  :math:`\Gamma_+ \to \Gamma_-` accordingly. A kernel outside the family (a
  general weak-form angular kernel, an angle-dependent reflectivity) is not a
  chain: it would be another
  :class:`~orpheus.transport.boundary_conditions.BoundaryCondition` type with
  its own realizer. The one shipped today is diffusion's
  :class:`~orpheus.diffusion.zero_flux.ZeroFluxCondition`, with return
  fraction :math:`-1`, outside :math:`K`.

  .. note::

     This bullet typed :math:`R` as an endomorphism
     :math:`\Gamma_- \to \Gamma_-` until 2026-08-04, and the
     :math:`G` bullet above claimed the crossing for :math:`G`
     unconditionally. As a **classification** the pair
     :math:`G : \Gamma_+ \to \Gamma_-`,
     :math:`R : \Gamma_- \to \Gamma_-` is coherent and it answers the
     question the taxonomy exists to answer; it is *not* the realized
     typing of either factor, and reading it as such is what the new
     :ref:`bc-taxonomy-vs-factorization` separates out.

* :math:`q \in \Gamma_-` is the **prescribed inflow source** — a
  vector-valued quantity on :math:`\Gamma_-` only, stored as the
  :class:`~orpheus.transport.boundary_conditions.InflowTerms` of a
  :class:`~orpheus.transport.boundary_conditions.WallCondition` and read for
  every condition type through ``BoundaryCondition.source``. The empty case
  :math:`q \equiv 0` is the homogeneous BC; the inhomogeneous case
  :math:`q \neq 0` is the rank-0 affine BC, a nonzero source over the empty
  network (``WallCondition.uniform(value)`` and ``WallCondition.framed(profile)``).
  It is **not** realized as an operator: :math:`q` is the one factor the
  realizer tier does not carry, travelling instead as the boundary leaf of the
  composite source (the table's ``SourceRealizer.sample`` returns an array) —
  see :ref:`bc-affine-source-channel`, which is also where the consequence of
  that split (vacuum and prescribed inflow realize to the *same* operator) is
  derived and measured.

On a face that ends the domain the factor :math:`G` is absent (the face's role is
a tag); on a glued face :math:`R = I` and :math:`q = 0` (its role is a deck). The
product form classifies the law's content; it is not the recipe by which the law
is evaluated (:ref:`bc-taxonomy-vs-factorization`).

Three remarks make this form load-bearing:

1. **Method-agnostic.** Nothing in :eq:`affine-bc-form` is SN-specific.
   The same affine map describes how MoC track-bundles, MC particle
   histories, CP boundary-to-region coupling matrices, and diffusion
   bilinear-form weak BCs all interact with the geometry. Each
   method's *realization* of the operators :math:`G`, :math:`R`,
   :math:`q` differs (see :ref:`bc-realizer-layer`); the algebraic
   shape of the law itself is universal.
2. **Affine, not linear.** The :math:`q` term is what makes the map
   affine. Most published transport-theory references treat the
   homogeneous case (:math:`q \equiv 0`) and never give the affine
   form an explicit name; ORPHEUS does because the two rank-0 cases, the
   vacuum (:math:`R = q = 0`, ``WallCondition.vacuum()``) and the prescribed
   inflow (:math:`R = 0` but :math:`q \neq 0`, ``WallCondition.uniform(value)``
   or ``.framed(profile)``), need a single uniform contract, and
   ``WallCondition`` is that contract: one type over the pair
   :math:`(\text{network}, \text{inflow})`. Note that it is :math:`R` alone
   that vanishes in both. In the trace-law architecture :math:`G` was the
   **identity** deck element there, not zero (the zero map is not a bijection
   and so cannot be a geometry map at all); today a tagged face has no deck
   element at all, and the same fact reads as "the face role is a tag".
   Writing ":math:`R = G = 0`" once spelled the same vanishing twice, once in
   the wrong tier; campaign phase B3 corrected it.
3. **The three factors are three values, and no object holds all three.**
   :math:`G` is the face's role (a deck), :math:`R` is
   ``WallCondition.network`` and :math:`q` is ``WallCondition.source``. The
   retired ``BoundaryTraceLaw`` ABC exposed ``geometry_map``,
   ``response_kernel`` and ``source`` as properties on every concrete law,
   defaulting to ``None / None / NoSource()``, so a cross-method realizer could
   introspect the geometric and response components separately; the SN realizer
   nevertheless dispatched on the law's class. The successor makes the
   separation the type structure: a realizer table is keyed by UNIT, the chain
   types of the response network, ``ZeroNetwork`` for the empty network,
   ``InflowTerms`` for the source, a method-only condition type where a method has
   one, and the one deck realizer for the decks, so a method never branches on
   the cone point of a wall. `[M]` 2026-10-24: ``RealizerTable.catalogue()`` of
   the five tables (``SN_REALIZERS``, ``DIFFUSION_REALIZERS``, ``CP_REALIZERS``,
   ``MC_REALIZERS``, ``MOC_REALIZERS``) lists only ``ZeroNetwork``, ``SpecularChain``
   and ``DiffuseChain`` in all five, ``InflowTerms`` in SN and diffusion and
   ``ZeroFluxCondition`` in diffusion; ``WallCondition in catalogue()`` is false for 5
   of 5.
   The geometric component, the response component and the source are therefore
   dispatched independently, which is what the retired properties were meant to
   allow and the dispatch on the law's class did not.

.. _bc-factor-roles:

Which factor owns what — the decidable criterion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two factors are not "the first thing that happens" and "the second
thing that happens". They are **different kinds of mathematical
object**, and membership is decidable — but it takes **two** tests, and
the first published form of this section shipped only one.

**The necessary test — multiplicativity.**

    :math:`G` is the **composition (Koopman) operator of a
    measure-preserving bijection of the boundary phase space** —
    :math:`(G\psi)(x) = \psi(g^{-1}x)` for some :math:`g` acting on
    :math:`\partial\Omega \times S^d`.

Such operators are invertible, preserve the trace measure
:math:`|\Omega\cdot\hat n|\,d\Omega\,dA`, form a **group**, and are
**multiplicative**: :math:`G(\psi\varphi) = (G\psi)(G\varphi)`. A
relabeling satisfies that identity; **an averaging operator never
does**. Anything *failing* multiplicativity is a kernel, and kernels
are :math:`R` — the argument that moved the Lambertian out of the
geometry slot at campaign phase B3.

**Why that test is not sufficient.** A *specular kernel* is a
permutation, hence multiplicative too. So multiplicativity alone cannot
separate

* a polished wall that returns :math:`\alpha` of the flux specularly,
  from
* a symmetry plane, across which the domain genuinely continues.

Both relabel; only one is geometry. Reading multiplicativity as *the*
criterion puts a surface's re-emission law in :math:`G`, and the two
objects have nothing in common but their matrix.

**The sufficient test — is it a quotient?**

    :math:`G` is the deck transformation **of an actual quotient of the
    physical domain** (:ref:`bc-factor-quotients`).

A physical surface is not a quotient: the domain does *not* continue on
the other side of it, so nothing is identified with anything and there
is no deck group to be an element of. A surface's specular pairing is
therefore **constitutive** — it is :math:`R`.

.. admonition:: The law this yields
   :class: important

   **Exactly one of** :math:`G`, :math:`R` **is non-trivial.**

   It is the contrapositive of *"*:math:`R = I` *exactly when the BC is
   a pure symmetry statement adding no physics"* — a law that asserts
   any physics at all has :math:`G = \mathrm{id}`. See
   :ref:`bc-factor-quotients` for the table, which every shipped value
   satisfies, and :ref:`bc-deck-law-not-an-operand` for why the algebra of
   conditions cannot rebuild a violation by scaling a deck.

In the present architecture the law is the shape of the role sum, not a
property each value is checked for. A face's role is a deck (:math:`G`, with
:math:`R = I` and no :math:`q`) or an integer tag whose bound condition is
:math:`(R, q)` (with :math:`G` absent), and ``FaceRole = DeckTransformation | int``
has no member that holds both. Law L-F gates the disjointness of the two slots
(``test_a_condition_is_not_a_deck_slot`` in
``tests/gates/transport/test_wall_condition_algebra.py``). `[M]` 2026-10-24, a
probe outside the suite that calls ``parse_face_role`` on seven condition values
(``WallCondition.vacuum()``, ``.specular(1.0)``, ``.specular(0.4)``,
``.white(0.4)``, ``.uniform(1.0)``, ``.from_amplitudes(0.3, 0.5)`` and
diffusion's ``ZeroFluxCondition()``) and ``BoundaryConditions({1: deck})`` on both
decks: 7 of 7 and 2 of 2 raise ``TypeError``; the positive controls
``parse_face_role(MirrorDeck(), ...)`` and
``BoundaryConditions({1: WallCondition.vacuum()})`` are accepted.

Physically: a change of direction caused by the **geometry** is
:math:`G`; a change of direction caused by the **constitutive
assumption of the BC** is :math:`R`. Absorption, accommodation and
diffusivity are :math:`R`. Mirrors, translations and rotations are
:math:`G` — but only when they are mirrors *of the domain*, not of a
wall standing in it. (A rotation is a deck class that is named and not yet
minted: ``DEFERRED_DECKS`` lists it, with the geometry value whose arrival mints it.)

**Whichever factor is non-trivial carries the crossing.** For a
**quotient** law the crossing is geometric, and the argument is exact:
the specular mirror
:math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n` is the unique
ambient isometry fixing the face; it exchanges the hemispheres and
preserves :math:`|\Omega\cdot\hat n|`. So across a symmetry plane the
passage from :math:`\Gamma_+` to :math:`\Gamma_-` is not something the
*physics* does — it is something the *geometry* provides, and
:math:`G`, not :math:`R`, carries it.

For a **constitutive** law there is no such isometry to lean on. A wall
is not a quotient — that is the sufficient test above — so nothing
identifies an outgoing direction with an incoming one geometrically,
and the crossing is performed by the *physics*: the response integrates
the outgoing flux and re-emits an incoming one. The realization is the
witness. The diffuse chain :math:`E \circ Rt` of the response network (the
white / diffuse-albedo wall) realizes as a **two-link chain** whose links state
the crossing out loud:
:class:`~orpheus.sn.boundary.angular.PartialCurrentOperator` collapses
the *outflow* angle axis, and
:class:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator`
broadcasts the result over the *inflow* one, so the composite types
itself :math:`\Gamma_+ \to \Gamma_-` (:ref:`bc-response-adjoint`); the
specular chain :math:`P` is likewise a
:math:`\Gamma_+ \to \Gamma_-` permutation. `[M]` 2026-10-24: **no realized
response is an endomorphism of** :math:`\Gamma_-`, because the declaration itself is
typed by the rungs of the trace ladder, and ``chain_composes`` (the predicate
"the factors' tiers compose from ``"outflow"`` to ``"inflow"``", evaluated at the
import of ``orpheus.transport.boundary_conditions.network``, which raises
``TypeError`` otherwise) is true of ``SpecularChain`` and ``DiffuseChain``, 2 of 2
members of the closed sum ``Chain = SpecularChain | DiffuseChain``; the control
``chain_composes((IsotropicSection, CurrentRetraction))`` is false. The claim about the
*realized* operators is pinned by ``test_realized_law_has_gamma_plus_domain`` in
``tests/gates/sn/operators/test_b3_domain_narrowing.py`` (leg A: the operator consumes
:math:`\Gamma_+`; leg B: fed a full-face probe it does not emit ``N`` rows, which is what
separates :math:`\Gamma_+ \to \Gamma_-` from :math:`\Gamma_+ \to \Gamma_+` when
:math:`|\Gamma_+| = |\Gamma_-|`); `[M]` 2026-10-24, ``-k gamma_plus_domain`` gives 10 passed.

Since exactly one of :math:`G`, :math:`R` is non-trivial, "the
non-trivial factor crosses" is well defined. (The law read "at most one" from
2026-10-01 until the carve: the closure-free ``AlbedoBoundary(1.0)`` of the retired
trace-law family had :math:`G = I` and :math:`R = I`, the second case below at
:math:`\alpha = 1`, and so it was the one value the law had to be weakened for. That value
has no successor, so "exactly one" holds again: a deck face has :math:`G \neq \mathrm{id}`
and :math:`R = I`, and a tag face has :math:`G` absent and a response, which is never the
identity, because the empty network is :math:`0` and every chain permutes or contracts.)
Two boundary cases sharpen the law rather than break it:

* **Rank-0 values** (the vacuum and the prescribed inflow, ``WallCondition.vacuum()``
  and the inflow facades). :math:`R = 0`, the empty network, so the composite is
  the zero map on :math:`\Gamma_+` landing in :math:`\Gamma_-`
  (``ZeroMap`` in SN): it *does* cross, vacuously, because no pairing
  information is needed to emit nothing.
* **A bare scalar response on an ANGULAR trace.** :math:`R = \alpha I`
  is non-trivial in *magnitude* and trivial in *angular structure* — it
  commutes with everything and therefore cannot pair directions — while
  :math:`G = \mathrm{id}` supplies no crossing either. So neither
  factor crosses, and the law is genuinely under-determined: this was
  the closure-free ``AlbedoBoundary(α)`` spelling that the SN realizer
  refused, and the same law is *complete* on a scalar trace, where
  :math:`J^- = \alpha J^+` exhausts the single degree of freedom. The spelling
  does not exist today: the linear part of a ``WallCondition`` is a combination
  of the two chains, and each names its crossing (:math:`P` or
  :math:`E \circ Rt`), so a response that crosses nothing cannot be written
  on the shared type. See :ref:`bc-method-realizability` for the scalar-vs-angular
  axis, which is the same distinction one level up.

.. note::

   This paragraph read *"The crossing is geometric … which is why*
   :math:`G` *and not* :math:`R` *carries it"* until 2026-08-04. The
   mirror argument it gives is correct and is retained above; what was
   wrong is that it was stated for **every** law, and a law with no
   isometry has no geometric crossing to inherit.

**Why a misassignment can hide — a theorem.** If :math:`R` is rank-one,
:math:`R = u \otimes v`, then

.. math::

   R \circ G \;=\; u \otimes \bigl(G^{\mathsf T} v\bigr),

a :math:`\Gamma_+ \to \Gamma_-` operator — the same type as the boundary
law's composite action, with :math:`u \in \Gamma_-` and
:math:`G^{\mathsf T} v \in \Gamma_+^{*}`. Now :math:`G` is the
composition operator of a measure-preserving bijection :math:`g`, and
that is precisely the hypothesis under which the three candidate
"transposes" coincide: :math:`G^{-1} = G_{g^{-1}}` because :math:`g` is
a bijection, and preserving the trace measure makes that inverse both
the Euclidean transpose and the metric adjoint (the same theorem the
deck-transformation row of
:ref:`the evaluation table <bc-taxonomy-vs-factorization>` cites). So
:math:`G^{\mathsf T} v = v \circ g`, and the Lambertian's
:math:`v = |\Omega\cdot\hat n|` is **preserved** by both the mirror and
the periodic translation — so :math:`v \circ g` is *the same function*
:math:`|\Omega\cdot\hat n|`, merely read on :math:`\Gamma_+` rather than
on :math:`\Gamma_-`. Hence :math:`R \circ G` comes out the SAME operator
whichever admissible :math:`G` was used:

    :math:`G` is **unobservable exactly when** :math:`R` is rank-one.

Read "unobservable" precisely: it means the composite does not depend on
*which* :math:`G` it was given. The hypothesis doing the work is that
:math:`v` is :math:`G`-invariant — true of every response this codebase
ships, because the only rank-one response is the Lambertian and its
:math:`v` *is* the trace measure's own weight. A hypothetical rank-one
response with a direction-selective :math:`v` (say
:math:`v = \delta_{\Omega_0}`) would make :math:`G` fully observable, so
do not read the slogan as "rank-one alone suffices".

The white BC is precisely the invariant case. Its :math:`G` slot
therefore had no observable consequence, and the physics drifted into it
— the cosine-weighted Lambertian average shipped as a geometry map until
campaign phase B3 moved it to the response tier where it belongs. The
same theorem is why that correction leaves the composite :math:`R\,G`
unchanged by construction.

In the present architecture the hiding place is gone, not merely policed: a tag
face has no :math:`G` slot, and a deck face has no :math:`R`, so a response
cannot be filed under geometry and the unobservability of :math:`G` is no longer a
way for a misassignment to survive. The theorem stands as the explanation of why
that slot was one.


.. note::

   **This step concluded** ":math:`R \circ G = R`" **until 2026-08-04.**
   That spelling type-checked only while the two half-traces were the
   same space in practice: with :math:`G : \Gamma_+ \to \Gamma_-` and
   the classifying :math:`R : \Gamma_- \to \Gamma_-`, the left side is
   :math:`\Gamma_+ \to \Gamma_-` and the right side is
   :math:`\Gamma_- \to \Gamma_-`, so they cannot be equal as operators.
   The step silently identified :math:`G^{\mathsf T} v` with :math:`v`
   by treating :math:`v` as the *function* :math:`|\Omega\cdot\hat n|`
   without tracking which half-trace it is restricted to — harmless
   before campaign phase **B3.2** (:ref:`bc-domain-narrowing`) narrowed
   the SN law onto :math:`\gamma_\pm` and made the two halves genuinely
   distinct spaces, a type abuse afterwards. **The theorem itself is
   untouched**: its content is the :math:`G`-independence of the
   composite, which is what makes :math:`G` unobservable and what makes
   the B3 correction safe. Nothing numerical changes.

.. _bc-factor-quotients:

The quotient picture — which BC is the orbifold
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`G` is the **deck transformation of the quotient by which the
physical domain is represented**:

.. list-table::
   :header-rows: 1
   :widths: 18 30 20 32

   * - BC
     - quotient
     - fixed points
     - what it is
   * - periodic
     - :math:`\mathbb{R}^d/\Lambda` by a translation
     - none — **free**
     - a torus; a genuine **covering space**, a manifold
   * - reflective
     - by a reflection
     - the mirror plane
     - an **orbifold** (Thurston *reflector* boundary)
   * - rotational (⅛-core)
     - by a finite rotation
     - the axis
     - an **orbifold** (cone points); a deck class named in ``DEFERRED_DECKS``
       and not minted

The orbifold label therefore attaches to **reflective**, not to
periodic — periodic is the free / covering-space case. And
:math:`R = I` **exactly when** the BC is a pure symmetry statement
adding no physics. Vacuum, white and albedo are not symmetry statements
at all, which is why none of them is a deck: each is a wall bound to a tag,
and all its content sits in :math:`R`.

Read as a *test* rather than as a remark, that sentence is the
sufficient criterion of :ref:`bc-factor-roles`, and it partitions every
implemented value:

.. list-table:: Which tier owns each value's content
   :header-rows: 1
   :widths: 28 22 26 24

   * - value
     - face role
     - :math:`R` (and :math:`q`)
     - what it asserts
   * - ``MirrorDeck()``
     - the deck itself, :math:`G`
     - :math:`I`
     - a symmetry plane — a quotient, **zero physics**
   * - ``TranslationDeck()``
     - the deck itself, :math:`G`
     - :math:`I`
     - a torus — a quotient
   * - ``WallCondition.specular(α)``
     - an integer tag, :math:`G` absent
     - :math:`\alpha P`
     - a **surface** returning :math:`\alpha` specularly
   * - ``WallCondition.white(α)``
     - an integer tag
     - :math:`\alpha\,(E \circ Rt)`
     - a surface returning :math:`\alpha` diffusely (Lambertian)
   * - ``WallCondition.from_amplitudes(s, d)``
     - an integer tag
     - :math:`s\,P + d\,(E \circ Rt)`, :math:`s + d \le 1`
     - a surface returning a specular fraction and a diffuse fraction (the Marshak mix)
   * - ``WallCondition.vacuum()``
     - an integer tag
     - :math:`0`
     - a surface returning nothing
   * - ``WallCondition.uniform(v)``, ``.framed(profile)``
     - an integer tag
     - :math:`0`, with :math:`q \neq 0`
     - a surface that emits inward and returns nothing
   * - ``ZeroFluxCondition()`` (diffusion's)
     - an integer tag
     - a method-only response, return fraction :math:`-1`
     - the Dirichlet idealization :math:`\phi = 0` at the face

**One pair realizes to one matrix and is still two values.**
``WallCondition.specular(α)`` realizes in SN to :math:`\alpha` times the matrix
of the mirror deck ``MirrorDeck()``. At :math:`\alpha = 1` the polished wall and
the mirror realize to one permutation matrix and are still two values: the wall
states how a surface returns the neutrons that reach it, the mirror states that
the solution beyond the plane is the mirror image of the solution inside it.
Keeping the *types* distinct is what makes "put a wall's response in the
geometry slot" unspellable, the exact error this section's earlier form
permitted. `[M]` 2026-10-24, a probe outside the suite (the realization of both
values through ``SN_REALIZERS.realize`` on the Cartesian ``xmin`` and ``xmax``
faces, the dense matrices read as ``tests/gates/sn/architecture/_realized_boundary_helpers.py``
reads them, compared with ``np.array_equal``): the realized matrix of
``WallCondition.specular(α)`` equals :math:`\alpha` times that of ``MirrorDeck()``
in 40 of 40 cases (five quadratures, ``gauss_legendre(8)``, ``gauss_legendre(16)``,
``product(2,4)``, ``level_symmetric(6)`` and ``lebedev(17)``, times two faces,
times :math:`\alpha \in \{0, 0.3, 0.7, 1\}`); the control, ``WallCondition.white(α)``
against the same product, differs in 30 of 30 cases with :math:`\alpha > 0`. Two
gates pin the equality against a reference that shares no code with the realizer:
``TestRealizeSpecular`` in ``tests/gates/sn/operators/test_sn_boundary_realizer.py`` asserts, with
``np.testing.assert_array_equal``, that ``WallCondition.specular(1.0)`` and ``MirrorDeck()`` each realize
to the hand-computed narrowed gather of the mirror partner on ``lebedev(17)``
(``test_specular_unit_weight_lebedev_matches_hand_computed``) and that the partial weight
:math:`0.7` realizes to :math:`0.7` times it (``test_specular_partial_weight_matches_hand_computed``);
`[M]` 2026-10-24, ``python -O -m pytest tests/gates/sn/operators/test_sn_boundary_realizer.py -k
"TestRealizeSpecular or uncertified"`` gives 7 passed.

The retired architecture had three more spellings of two of these values, and the
fold is one of the carve's decisions. ``AlbedoBoundary(α, SpecularReturn(a))`` was
the specular wall, ``AlbedoBoundary(α, IsotropicReturn(a, s))`` and
``WhiteBoundary(a, s, α)`` were the diffuse wall under two names, and each pair
shared one realization body so that the equivalence held by construction; the
gate ``tests/gates/geometry/test_reemission_closure.py::TestEquivalenceTheorems::test_specular_closure_equals_reflective``
asserted the specular one with ``np.array_equal`` (`[M]` 2026-10-01, a probe outside
the suite: 20 of 20 realized actions bitwise equal, five quadratures from
``gauss_legendre(4)`` to ``lebedev(17)`` on an x or y face, times
:math:`\alpha \in \{0, 0.3, 0.7, 1\}`). Today each pair is one value:
``WallCondition.specular(α)`` and ``WallCondition.white(α)`` are the only
spellings, and ``WallCondition.specular(0.0)``, ``.white(0.0)`` and ``.vacuum()``
are one value with one digest (`[M]` 2026-10-24: ``specular(0.0) == white(0.0) == vacuum()``
is true), because :math:`R = 0` is one point of the cone :math:`K` however it is
reached.

Every shipped value satisfies the criterion, by the shape of the role sum rather
than by a check per class: a face's role is a deck or a tag, so no value has both
factors different from the identity (the disjointness gate and its measurement are
in :ref:`bc-factor-roles`). The two decks have :math:`R = I`, and every condition
type, ``WallCondition`` and diffusion's ``ZeroFluxCondition``, has :math:`G`
absent. The one value on which both factors were the identity, the closure-free
``AlbedoBoundary(1.0)``, stated no crossing at all and the SN realizer refused it
(:ref:`bc-method-realizability`); it has no successor, because a ``WallCondition``'s
network is a combination of the two chains, each of which names its crossing. In
particular ``MirrorDeck`` has no field (``test_a_deck_is_its_type_and_has_no_field``
in ``tests/gates/geometry/test_deck.py``), so a mirror with an amplitude is not
constructible: it would have both factors non-trivial, the partial wall
``WallCondition.specular(α)`` spelled a second time in the geometry slot, and a
consumer that recognises a mirror by its type would read that second spelling as a
perfect mirror. That is the defect class ERR-094 records (three such readers in
S\ :sub:`N`); the history of the retired parameter is in :ref:`bc-history`.

.. _bc-deck-law-not-an-operand:

A deck law is not an operand
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A deck is a symmetry statement, and a symmetry has no amplitude. The first form of
this page let the descriptor algebra of :ref:`bc-rank-n-algebra` rebuild an attenuated
mirror even after the mirror's own albedo parameter was retired:
``0.7 * ReflectiveBoundary("x")`` was a ``LawScaled`` whose realization was
:math:`0.7` times the mirror's permutation, the same matrix as the partial wall and
the same type confusion. The trace-law family closed that route with a guard: the two
composition nodes, ``LawScaled`` and ``LawSum``, refused a **deck law** (a law whose
geometry factor :math:`G` was not the identity deck element: ``ReflectiveBoundary`` and
``PeriodicBoundary`` among the shipped laws) with a ``TypeError`` naming the law and
the spelling to use instead. Wave 2 (2026-10-24) retired the family, and the refusal
with it, because the successor makes the construction unspellable instead of
refused.

**The present structure.** A deck is a role, not a condition. ``MirrorDeck`` and
``TranslationDeck`` define no operator, and the arithmetic of conditions
(``+``, ``-``, unary ``-`` and scalar ``*``, one definition shared by
``WallCondition`` and ``WallDisplacement``) accepts only a pair of that arithmetic
and answers ``NotImplemented`` to anything else, so the scaled or summed deck is
Python's own ``TypeError``, raised from the operator protocol with no frame of
``orpheus`` in its traceback (law L-F; the gate is ``test_a_deck_is_not_an_operand``
in ``tests/gates/transport/test_wall_condition_algebra.py``, which also holds the
positive control ``0.5 * WallCondition.vacuum()``). `[M]` 2026-10-24, run from the
worktree root of this commit:

.. code-block:: text

   0.7 * MirrorDeck()                    -> TypeError: unsupported operand type(s) for *: 'float' and 'MirrorDeck'
   MirrorDeck() + MirrorDeck()           -> TypeError: unsupported operand type(s) for +: 'MirrorDeck' and 'MirrorDeck'
   WallCondition.vacuum() + MirrorDeck() -> TypeError: unsupported operand type(s) for +: 'WallCondition' and 'MirrorDeck'
   0.5 * WallCondition.specular(1.0)     -> WallCondition   (== WallCondition.specular(0.5))

The induction that the retired nodes needed (below) is not needed here: the nodes could
nest, so "no admitted tree holds a deck law at any depth" had to be proved on the
height of the tree, whereas the operand set of the arithmetic is one closed set,
``WallCondition``, ``WallDisplacement`` and the real scalars that are not ``bool``,
and every result is a member of it. There is no tree to descend.

**Why a symmetry cannot be scaled.** A deck law asserts that the domain
is a fundamental domain of a symmetry group: the solution on the far
side of the face is the image of the solution on this side under the
group element. That assertion has no amplitude. Neutrons are not
returned by anything; the face is not a surface at all. Scaling it by
:math:`c` has two readings, and both are something else. Read as
physics, :math:`c\,G` with :math:`c < 1` is a surface that loses a
fraction :math:`1 - c` of its outflow and returns the rest along the
mirrored direction, which is ``WallCondition.specular(c)``, a response
whose role is a tag and whose :math:`G` is absent. Read as a symmetry, it is
nothing: the image of a solution under a rigid motion is the solution itself,
not a fraction of it.

**Why a symmetry cannot be mixed.** A sum containing a deck is one
of two things. If the other terms are also decks on the same face,
the sum is a scaled deck in disguise (``0.5 * R + 0.5 * R`` with
``R`` the mirror realizes to ``R``'s matrix, the costume again). If any other
term is a response, the sum is not sub-Markov in general: a mirror plus a white
wall of albedo 1 returns the outflow twice, once specularly and once
diffusely, so the incoming partial current is twice the outgoing one. The same
arithmetic is what the cone predicate catches among responses: the sum
``WallCondition.specular(1.0) + WallCondition.white(1.0)`` has return fraction
:math:`2`, outside :math:`K`, so it is a ``WallDisplacement`` and never bindable,
while a legitimate mix of specular and diffuse return (the Marshak boundary) is
a mix of two responses with amplitudes inside the cone,
``WallCondition.specular(0.5) + WallCondition.white(0.5)``, which equals
``WallCondition.from_amplitudes(0.5, 0.5)`` and is a ``WallCondition``
(`[M]` 2026-10-24, the probe of the code block above).

**The predicate is the type, never the realized matrix.**
``WallCondition.specular(1.0)`` permutes ordinates exactly as the mirror does and
composes freely (``0.5 * WallCondition.specular(1.0)`` is ``specular(0.5)``), because
it is a condition; the question "is this role a deck" is
``isinstance(role, DeckTransformation)``, which ``parse_face_role`` asks. What a
deck *does* is read off its pairing (``FacePairing.permutes_ordinates`` is derived
from the rigid motion, ``not motion.is_translation``), and that property is for
realizers to read, never for a consumer to decide what a deck is. The retired
predicate was ``law.geometry_map.is_identity``, a property of ``SelfPairedDeck`` and
``PairedDeck`` read off the rigid motion (``RigidMotion.is_identity``, the
conjunction of a linear part equal to :math:`I` and a zero translation), with
``PairedDeck.is_identity`` ``False`` for every constructible value, since a paired
deck that fixed every face pointwise would be the self-paired case its constructor
refused.

**The retired guard.** The refusal in the composition nodes was tagged
``ELEGANCE-DEBT[guard] #551`` in ``orpheus.geometry.boundary._composition``: the deck
laws were to leave the law family for the geometry (the boundary-law ontology's
architecture E), after which the composition nodes would take responses only and
scaling a symmetry would be unspellable rather than refused. Wave 2 is that
retirement (issue #551 is the issue that tracked it). The retired check on direct
children was sufficient by induction on the height of the tree: a leaf law is a tree
of height 0 and is not itself an operand, and a node of height :math:`h + 1` was
admitted only if each direct child was a law that passed the check or a node of height
at most :math:`h`, which by the induction hypothesis held no deck law; the check had
to live in ``__post_init__`` because every construction route reached it (the
arithmetic dunders, a direct ``LawSum(a, b)`` or ``LawScaled(c, x)``,
``dataclasses.replace``, ``sum`` and ``functools.reduce``) and a check in the dunders
alone was bypassed by a direct construction. Two limits, both measured by the review
(2026-10-01): the nodes did not type their children, so a child that was not a law at
all (``sum()``'s starting ``0``) was admitted there and refused later, by
``realize_recursively``; and a route that bypasses construction (``object.__new__``, or
loading a pickle written before the check existed) bypassed the check too. The first
limit does not survive: ``sum([...])`` from the integer ``0`` is a ``TypeError`` of the
operator protocol (``WallCondition`` defines no ``__radd__``, `[M]` 2026-10-24:
``sum([specular(1.0), white(1.0)])`` raises ``unsupported operand type(s) for +: 'int'
and 'WallCondition'``). The second does: ``object.__new__`` skips ``__post_init__``, which
is where the cone check of ``WallCondition`` lives, so a condition built that way can lie
outside :math:`K` (`[M]` 2026-10-24: ``object.__new__(WallCondition)`` with a specular weight
of ``1.5`` reports ``return_fraction`` 1.5, where ``WallCondition(network)`` with the same
network raises ``SubmarkovViolationError``); the standing defence is that every constructor and every arithmetic
operator goes through ``__post_init__`` or ``_combine``.

**Gates.** The successor gates are the L-F gates of
``tests/gates/transport/test_wall_condition_algebra.py``
(``test_a_condition_is_not_a_deck_slot`` and ``test_a_deck_is_not_an_operand``) and
``tests/gates/geometry/test_deck.py`` (``test_a_deck_is_its_type_and_has_no_field``);
`[M]` 2026-10-24: ``python -O -m pytest tests/gates/transport/test_wall_condition_algebra.py
tests/gates/geometry/test_deck.py tests/gates/geometry/test_face_roles.py`` gives
52 passed. The retired gates were ``tests/gates/geometry/test_reflective_is_a_mirror.py``
(36 rows: an albedo refused as a positional argument, as a keyword and through
``dataclasses.replace``; ``axis`` the only field; ``kind`` the registry key on every
axis; the two amplitude checks retired while the certification still ran; :math:`G`
the mirror and :math:`R` the identity; value semantics, including that a mirror and a
perfect specular wall are different values) and
``tests/gates/geometry/test_deck_laws_do_not_compose.py`` (217 rows: three deck laws,
``ReflectiveBoundary``, ``PeriodicBoundary`` and a synthetic mirror, against 21
composition routes covering the arithmetic dunders of the laws and of both node types,
the direct constructors and the nested spellings; ``dataclasses.replace``; the nested
refusal located at the inner node; and the positive legs, seven responses including the
polished wall composing on every route and a response mix realizing through
S\ :sub:`N` and through diffusion as before). `[M]` 2026-10-01: before that carve 12 of
36 and 67 of 217 failed; after it, 253 of 253 passed.

.. note::

   **SN apply matvec honours the affine BC contract (Issue #168
   Phase C, 2026-05-12).** The then-production per-geometry matvecs
   (``transport_operator_matvec_spherical`` and ``_cylindrical`` —
   since deleted in the typed-field campaign (#197), their successor
   ``_transport_operator_matvec_unified`` in turn retired at the walk
   unification (#280 campaigns)) were rewritten as one sweep
   iteration semantically: the realized boundary operator is applied **at least
   once** per matvec at the boundary edge on the WDD-propagated
   outflow face values (:math:`\gamma_+ \psi`), not on cell-centre
   approximations.  The live forward action is now
   :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`
   through the loss-representation walk. The pre-Phase-C cell-centre-as-face-value
   contamination — and the Phase A ``BoundaryFaceFlux`` Protocol
   that patched it — both retire in Phase C. See
   :ref:`bc-trace-contract-respected-by-matvec` for the
   verification gate that pins this contract, and
   :ref:`bc-two-bc-applies-per-matvec` for the Phase D
   strengthening that audits **two** BC apply calls per matvec
   (Phase D Carlson context + Phase C realized boundary operator).

.. _bc-taxonomy-vs-factorization:

The taxonomy and the factorization are different questions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Everything the two sections above say about :math:`R \circ G` is a
**taxonomy**. Two successive designs in this campaign were built on the
premise that it is also a *recipe for evaluation*, and both were refuted
by the tree's own realization — which is why the distinction gets its
own section rather than a parenthesis.

**What** :math:`R \circ G` **is.** A classification of the law's
*content*: the multiplicativity test plus the quotient test decide
whether a boundary law asserts geometry or physics
(:ref:`bc-factor-roles`, :ref:`bc-factor-quotients`). In that register
the typing :math:`G : \Gamma_+ \to \Gamma_-`,
:math:`R : \Gamma_- \to \Gamma_-` is perfectly coherent and it answers
the only question the taxonomy exists to answer. What makes the
classification load-bearing is that it is **decidable** and that it
makes "put a wall's response in the geometry slot" unspellable.

**What it is NOT.** How a face's condition is *evaluated* — and,
decisively, how its adjoint is obtained — follows its **kind** (a deck, or a
response), not the two-factor product:

.. list-table:: Evaluation follows the KIND, not the classifying product
   :header-rows: 1
   :widths: 16 40 44

   * - kind
     - structure
     - adjoint
   * - **deck transformation**
     - **Atomic.** A measure-preserving bijection does not factor into
       two meaningful pieces. Pure geometry — a law imposed by theorem
       and transport-method-agnostic apart from needing to know *which
       space* it acts on. It is :math:`\Gamma_+(f) \to \Gamma_-(f)` for a
       self-paired face and :math:`\Gamma_+(f') \to \Gamma_-(f)` through
       a genuine face pair.
     - A **theorem**: the composition operator of a bijection :math:`g`
       is invertible with :math:`G^{-1} = G_{g^{-1}}`, and
       measure-preservation makes that inverse the transpose. Nothing to
       verify per law.
   * - **response**
     - :math:`N` **composed** operations
       :math:`\Gamma_+(f) \to \cdots \to \Gamma_-(f)`. The Lambertian is
       an outflow **angle contraction** :math:`C` followed by an
       **isotropic broadcast** :math:`B`, so :math:`R = B \circ C` (:math:`B` is the section :math:`E` and :math:`C` the current retraction :math:`R_t` of the response network of :ref:`bc-architecture-of-record`; the section writes :math:`E \circ R_t` and :math:`B \circ C` for the one chain), with
       an intermediate state in the angle-integrated per-face
       **scalar-current space** :math:`S(f)`. Constitutive — an
       assumption about a real surface.
     - **Conditional** on :math:`S(f)` carrying a non-degenerate metric
       — see :ref:`bc-response-adjoint`.

:math:`S(f)` is not bookkeeping: it is the **outgoing partial current**
:math:`J^+_g = \sum_{\Gamma_+} w\,\lvert\Omega\cdot\hat n\rvert\,\psi`,
a real physical quantity with units, and it is exactly what the
Lambertian's realization used to compute and immediately throw away as
the anonymous local ``psi_avg`` inside the welded
``AngularAverageOperator``'s apply. Both halves shipped at **G6.3**
(2026-08-04): the per-face accessor is
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.current_space`
(step 2) and the factored :math:`B \circ C` spelling is
:class:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator`
:math:`\circ`
:class:`~orpheus.sn.boundary.angular.PartialCurrentOperator` (step 3).

.. note::

   **The host is a new tier of the** :math:`\Gamma` **ladder, not**
   :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`.
   This paragraph predicted the latter — "the type that will host it
   exists" — and the prediction did **not** hold. ``S(f)`` is not a
   subset of ordinates but the *collapse of the ordinate axis*, one
   value per group per boundary cell, carrying a deliberately **unit**
   metric because the angular measure
   :math:`\lvert\Omega\cdot\hat n\rvert \odot w` has already been
   consumed by the contraction that lands there.
   :class:`~orpheus.numerics.spaces.scalar_trace_space.ScalarTraceSpace`
   carries the :math:`(J^+, J^-)` **pair** for the whole boundary under
   the face-**area** metric — diffusion's P1 Cauchy data. Same physical
   name, different object, different metric, different scope; hosting
   one in the other would have double-counted the area weight.

.. important::

   **Why the conflation was invisible, and what it cost.** Reading the
   classifying product as the computational one is what made
   ":math:`R \circ G = R`" look sound (it is a type error the moment
   :math:`\Gamma_\pm` are distinct spaces —
   :ref:`see the theorem <bc-factor-roles>`) and what made "the crossing
   is geometric" look general (it is proven for the mirror only). Both
   claims shipped — on this page and in the
   ``orpheus.geometry.boundary`` Protocols — and both were corrected
   on 2026-08-04.
   The transferable lesson is the campaign's own: **a declaration tier's
   typing is a declaration, not a measurement — check it against the
   realization before designing on it.** In this case one read of the
   then-welded ``AngularAverageOperator``'s first line, which typed
   itself :math:`\Gamma_+ \to \Gamma_-`, refuted both — and the typing
   survived the operator: the chain that replaced it composes
   :math:`\Gamma_+ \to S(f) \to \Gamma_-`, so the refutation now reads
   off the two links' bound spaces instead of one docstring.

.. _bc-response-adjoint:

When a response has an adjoint — the intermediate metric cancels
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A deck transformation's adjoint is free (the theorem above). A
**response**'s is not, because factoring it introduces an intermediate
space :math:`S(f)` whose metric appears in each factor's Hilbert
adjoint. Write :math:`G_+`, :math:`G_-`, :math:`G_S` for the metrics on
:math:`\Gamma_+`, :math:`\Gamma_-`, :math:`S(f)`, and recall the Hilbert
adjoint of a map between metric spaces,
:math:`A^{*} = G_{\mathrm{dom}}^{-1} A^{\mathsf T} G_{\mathrm{cod}}`.
For :math:`R = B \circ C` with :math:`C : \Gamma_+ \to S` and
:math:`B : S \to \Gamma_-`,

.. math::
   :label: bc-response-factored-adjoint

   R^{*} \;=\; C^{*} B^{*}
   \;=\; \bigl(G_+^{-1} C^{\mathsf T} G_S\bigr)
         \bigl(G_S^{-1} B^{\mathsf T} G_-\bigr)
   \;=\; G_+^{-1} C^{\mathsf T} B^{\mathsf T} G_-
   \;=\; G_+^{-1} R^{\mathsf T} G_- .

.. (vv-status history) This equation carried
   ``.. vv-status: bc-response-factored-adjoint documented`` until G6.3 step 3b
   (2026-08-04). The sentinel's rationale was: "The FACTORED spelling does not
   exist in code yet — the Lambertian ships as one operator with
   ``is_adjointable=False``, so there is no function for a ``verifies`` marker
   to point at." Both halves shipped — ``S(f)`` is
   ``AngularTraceSpace.current_space``, the factorization is
   ``IsotropicEmissionOperator @ PartialCurrentOperator`` — so the precondition
   expired and the directive is removed rather than left as a standing claim
   that the equation is documented-only. It is now GATED on the shipped chain:
   ``tests/gates/sn/operators/test_lambertian_chain.py::
   TestReciprocityAgainstTheMirrorFace::test_H_is_pointwise_the_mirror_face_kernel``
   carries ``verifies("bc-response-factored-adjoint")``, alongside six
   abstract-matrix gates in ``tests/gates/numerics/test_factored_adjoint_identity.py``.
   Its PRECONDITION is separately gated:
   ``tests/gates/numerics/test_angular_face_trace_space.py::
   test_the_half_trace_metric_is_strictly_positive`` pins the non-degeneracy
   this equation requires of the intermediate, and
   ``test_the_metric_is_not_euclidean`` pins that the metric is load-bearing.

**The intermediate metric cancels.** :math:`G_S` appears once as a
factor and once as its inverse, so it drops out of the composite
entirely — which means the requirement on :math:`S(f)` is exactly **one
binary condition: the metric must be non-degenerate.** Not "the
physically correct metric". Measured with :math:`|\Gamma_+| = 7 \neq
|\Gamma_-| = 5` deliberately (so no accidental index pairing can mask
an error), sweeping :math:`G_S` over eleven orders of magnitude:

.. list-table:: :eq:`bc-response-factored-adjoint`, measured 2026-08-04
   :header-rows: 1
   :widths: 34 30 36

   * - :math:`G_S`
     - :math:`\max\lvert R^{*}_{\text{factored}} - R^{*}_{\text{direct}}\rvert`
     - reading
   * - identity (Euclidean)
     - ``0.000e+00``
     - exact
   * - ``1e-6``
     - ``1.110e-16``
     - one ULP — cancels
   * - ``3.7e5``
     - ``1.110e-16``
     - one ULP — cancels
   * - **``0`` (degenerate)**
     - **``7.628e-01``**
     - ⚠ **BROKEN** — no adjoint exists

and the weighted adjoint law
:math:`\langle Rx, y\rangle_{G_-} = \langle x, R^{*}y\rangle_{G_+}`
holding at exactly ``0.0``.

The probe is a G6.3 design-session artefact and is reproducible in a
dozen lines, so it is described rather than pinned to a path: build
:math:`C = (\cos\!w/\operatorname{norm})` as a :math:`(1, 7)` row and
:math:`B = \mathbf{1}` as a :math:`(5, 1)` column (so
:math:`R = BC` is the rank-one Lambertian shape; the *shipped* split
puts the :math:`1/\operatorname{norm}` in :math:`B` instead, since
:math:`C` must produce a current and :math:`B` an intensity — the
composite, and therefore this measurement, is unaffected), take
:math:`G_+`,
:math:`G_-` diagonal with entries drawn from :math:`[0.5, 2]`, and
compare :math:`C^{*}B^{*}` against :math:`G_+^{-1}R^{\mathsf T}G_-`
for each :math:`G_S` in the table. The asymmetric sizes are the
load-bearing part of the design: at :math:`|\Gamma_+| = |\Gamma_-|` an
index-pairing error would be invisible, and equal sizes are `[M]` what
**every** production quadrature happens to give — measured
:math:`|\Gamma_+| = |\Gamma_-|` on every face of
``gauss_legendre(4)`` (2/2), ``gauss_legendre(8)`` (4/4),
``product(2,4)`` (2/2), ``product(4,4)`` (4/4),
``level_symmetric(6)`` (24/24) and ``lebedev(17)`` (49/49). That
equality is an **accident, not a contract**, which is why the retired
``_narrowed_zero_operator`` refused to lean on it. Its successors take three positions.
The zero image (``ZeroMap``) reads both half-spaces off the trace ladder
(``trace.outflow_space(face)``, ``trace.inflow_space(face)``) and needs no index pairing. The
permutation image (``ScaledOrdinatePermutation``) pairs by index and therefore ASSERTS the
bijection at construction instead of assuming it from the accidental equality: equal sizes,
inflow-to-outflow membership and measure carried (ERR-045 is the crossing check; a size
mismatch raises ``ValueError``). A factored response, routed through :math:`S(f)`, never pairs
by index at all.

Two consequences to carry:

1. ⚠ **Only the COMPOSITE is metric-free.** :math:`C^{*}` and
   :math:`B^{*}` *individually* depend on :math:`G_S`, and factoring
   exists precisely so that the factors become usable — so
   :math:`S(f)`'s metric must still be chosen deliberately. What the
   cancellation guarantees is that no admissible choice can make the
   composite wrong.
2. **The full per-face tier can never be the intermediate.**
   :math:`\Gamma(f)`'s metric :math:`|\Omega\cdot\hat n| \odot w_n`
   **vanishes** on the tangential ordinates
   (:math:`|\Omega\cdot\hat n| \le \varepsilon_{\text{tan}}`), so it is
   only semi-definite; the two half-traces are strictly positive
   precisely because they exclude those rows
   (:ref:`bc-trace-structure`). That retroactively gives the G6.1 gate
   ``test_the_half_trace_metric_is_strictly_positive`` its real reason:
   it reads as hygiene and is in fact the precondition for
   :eq:`bc-response-factored-adjoint` to hold.

**The payoff — the deferred transpose fell out, and shipped.** The
welded ``AngularAverageOperator`` reported ``is_adjointable = False``
and deferred its transpose to campaign phase **B5**. Factored at
**G6.3 step 3b** there was nothing left to defer — each link has ONE
honest transpose:
:math:`C^{\mathsf T}(s) = \cos\!w \otimes s` (the outer product, on
:meth:`~orpheus.sn.boundary.angular.PartialCurrentOperator.apply_transpose`)
and :math:`B^{\mathsf T}(\varphi) = \bigl(\sum_{\Gamma_-} \varphi\bigr)
/ \operatorname{norm}` (the sum over inflow, on
:meth:`~orpheus.sn.boundary.angular.IsotropicEmissionOperator.apply_transpose`),
so :math:`R^{\mathsf T} = C^{\mathsf T} B^{\mathsf T}` is

.. math::

   R^{\mathsf T}(\varphi) \;=\;
   \frac{\cos\!w}{\operatorname{norm}} \sum_{\Gamma_-} \varphi ,

verified `[M]` bit-exactly (``max|Rᵀφ − (cos_w/norm)·Σφ| = 0.0`` on
``product(2,4)``, ``xmax``) against the dense transpose of the chain's
own ``apply``, and :math:`\le` 1 ULP on ``gauss_legendre(8)``,
``level_symmetric(6)`` and ``lebedev(17)``. Note where the
normalisation sits: in :math:`B`, not :math:`C`, because :math:`C`
produces a *current* and :math:`B` must produce an *intensity* — which
is what leaves :math:`S(f)` carrying an honest :math:`J^+` and lets an
albedo enter as the pure scalar law :math:`J^- = \alpha J^+`.
``LambertianReemission.is_adjointable`` (an attribute of the retired
trace-law type) flipped to ``True`` in the same step.

.. warning::

   That one-liner is the **Euclidean transpose** :math:`R^{\mathsf T}`,
   not the Hilbert adjoint. The two differ by exactly the re-weighting
   in :eq:`bc-response-factored-adjoint`:
   :math:`R^{*} = G_+^{-1} R^{\mathsf T} G_-`, and on a trace the metric
   is **never** Euclidean (it is :math:`|\Omega\cdot\hat n| \odot w_n`
   — the ERR-067 family is what happens when a half-trace pairing drops
   it). Advertising :math:`R^{\mathsf T}` under a name that reads as
   "the adjoint" is exactly the two-``.T``-semantics ambiguity the
   welded operator declined to resolve — and note how the factoring
   *dissolved* rather than settled it: each link's ``apply_transpose``
   is unambiguous on its own, so no link ever had to choose a reading,
   and the composite's Hilbert adjoint follows from its bound spaces.
   The honest channel is therefore still ``.H``, carrying the metrics
   of the
   :class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`
   the chain is bound to — measured `[M]` equal to
   :math:`G_+^{-1} R^{\mathsf T} G_-` to :math:`\le` 1 ULP, with
   :math:`\langle Rx, y\rangle_{G_-} - \langle x, R^{*}y\rangle_{G_+}`
   at :math:`\le` 1 ULP, on ``gauss_legendre(8)``, ``product(2,4)``,
   ``level_symmetric(6)`` and ``lebedev(17)`` at ``xmax``.

**The adjoint of the shared family, stated once.**
The factored adjoint above is the identity every response obeys; the family of wall
conditions and decks adds four statements, each gated.

* *Reciprocity.* Every shared kernel is reciprocal through the reversal of the
  ordinates. Write :math:`B` for the realized map :math:`\Gamma_+ \to \Gamma_-` of one
  face, :math:`M_\pm = \operatorname{diag}(w\,|\mu|)` for the half-trace metrics read off
  the ladder (:math:`M_+` on :math:`\Gamma_+`, :math:`M_-` on :math:`\Gamma_-`), and
  :math:`P_{\rm rev}` for the :math:`\mu \to -\mu` reversal, built by matching direction
  cosines and never from the realizer. The gated identity is
  :math:`M_+^{-1} B^{\mathsf T} M_- = P_{\rm rev}\,B\,P_{\rm rev}` on each realized shared
  image. For the specular return :math:`B = \alpha P_\pi`, and the identity is the
  involutivity and the measure preservation of :math:`\pi`; for the Lambertian
  :math:`B = d\,\mathbf 1_- \otimes (w|\mu|)_+ / M_-`, with :math:`M_- = \sum_{\Gamma_-} w|\mu|`
  the INFLOW mass. (On the full space, where the certified mirror pairing :math:`\pi` maps
  the trace to itself, the statement reads :math:`M^{-1} P^{\mathsf T} M = P`. That is a
  statement about :math:`\pi`; the realized map goes between two different spaces and has no
  self-adjointness of its own, which is the point of the retirement of the flag described
  in :ref:`bc-narrowed-involution`.) The adjoint of a shared condition's
  realization is therefore its realization on the reversed ordinates, and a forward solve
  and an adjoint solve share one binding (law L-A (i), ``test_reciprocity`` in
  ``tests/gates/sn/architecture/test_realized_boundary.py``: ``array_equal`` for the
  specular image, ``1e-15`` for the Lambertian, with a non-reciprocal control,
  :math:`R(\mu' \to \mu) \sim \mu'^2`, that reads a defect orders above the bound).
* *The left-eigenvector rider.* The scalar :math:`\rho` that consumers read (the
  return fraction, read by DSA, the corner swap, the leakage term and ERR-094) exists
  exactly when the trace measure is a left eigenvector of the realized map,
  :math:`R^{\mathsf T}(w|\mu|)_- = \rho\,(w|\mu|)_+`, equivalently
  :math:`J^-(R\varphi) = \rho\,J^+(\varphi)` for every :math:`\varphi \ge 0`. For
  the shared family it holds by construction, because the section :math:`E` divides by the INFLOW
  mass :math:`M_-` (so :math:`Rt_{\rm in} \circ E = \mathrm{id}`) and :math:`P` preserves the
  partial current. It is gated as the array identity
  ``op.apply_transpose(w|mu| on Gamma_-) == rho * (w|mu| on Gamma_+)``, which a wrong
  :math:`\rho` cannot pass and a one-vector ratio can. `[M]` 2026-10-24, a probe outside the
  suite that runs the gate's body (``test_the_left_eigenvector_rider_through_the_realizer``
  in ``tests/gates/sn/architecture/test_sn_realizer_table.py``) through ``SN_REALIZERS.realize``:
  the identity holds in 16 of 16 cases (four quadratures, ``gauss_legendre(8)``,
  ``product(2,4)``, ``level_symmetric(6)`` and ``lebedev(17)``, times faces ``xmin`` and
  ``xmax``, times ``specular(0.7)`` with ``array_equal`` and ``white(0.3)`` within
  ``ceil(N/3)`` ulp; the maximum measured deviation of the white wall is 0 ulp); the
  control, the same transposed image compared with :math:`0.5\,(w|\mu|)_+` for the
  ``specular(0.7)`` wall, is not equal.
* *The deck's adjoint is its partner's deck.* :math:`G_g^{\#} = G_{g^{-1}}`, which is the deck on
  the partner face: ``pairing_of(d, f).inverse == pairing_of(d, partner(f))`` (law L-A (ii),
  ``test_pairing_of_and_its_adjoint`` in ``tests/gates/geometry/test_deck.py``).
* :math:`q` *is not adjointed.* ``solve_sn_adjoint_fixed_source`` builds a zero boundary
  source (``AngularBoundarySourceSink.zeros``) and never reads the declared :math:`q`. This
  is correct, because :math:`q` enters the reciprocity identity as
  :math:`\langle q, |\Omega\cdot n|\,\gamma_-\psi^{+}\rangle` against the adjoint flux; the
  gate is ``test_the_adjoint_ignores_q`` in
  ``tests/gates/sn/architecture/test_sn_boundary_binding.py`` (law L-A (iii)). The refusal
  is a separate matter: the two eigenvalue entries, ``solve_sn`` and ``solve_sn_adjoint``,
  refuse a declared boundary source because the pencil :math:`(A, F)` has no source term
  (``_refuse_a_boundary_source``, ERR-095).

That is the campaign's thesis in one line: **the adjoint falls out of
well-posedness instead of being hand-rolled.**

.. _bc-deck-length-one-chain:

The deck transformation is the same chain, one link long
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two-link chain above invites a taxonomy — *composed* responses
against *atomic* deck transformations — and that taxonomy is a trap. It
would make "atomic" a **kind**, and a kind needs its own code arm. It is
not a kind; it is a **degenerate length**. A boundary law is a sequence
whose first link has domain :math:`\Gamma_+(f)` and whose last has
codomain :math:`\Gamma_-(f)`, with the interior determined by whatever
the physics needs; a measure-preserving bijection has nothing to factor,
so its sequence has one link. Specular reflection is therefore built by
the same body as the Lambertian, bound to the same two spaces, and
composed by the same ``@`` — there is no second path (G6.3 step 5, issue
**#330**).

The realizer tables state it in one type. The SN image of ``alpha * SpecularChain``,
of ``MirrorDeck`` and of ``TranslationDeck`` is the one class
``ScaledOrdinatePermutation`` (a mirror deck at ``alpha = 1.0``, the wrap at
``alpha = 1.0`` with the identity permutation across the two faces, a partial specular
wall with its weight), built through one crossing function bound to the same two trace
spaces. A symmetry has no amplitude, so the deck is that type at weight one and
nothing else distinguishes it from the unit specular wall at the level of the
realized operator; what distinguishes them is the role they were declared in
(:ref:`bc-factor-quotients`).

.. _bc-narrowed-involution:

The involution that isn't — a flag retired in favour of the algebra
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A mirror is an involution — :math:`\sigma^2 = \mathrm{id}` is the whole
content of a symmetry plane, and ERR-044 guards it. It is therefore
natural to expect the *realized* specular operator to report itself an
involution, and until G6.3 step 5
:class:`~orpheus.numerics.elementary_operators.PermutationOperator` carried an
``is_involution`` attribute — set at construction from
``perm[perm] == np.arange(N)`` — to say so. **It could not, and the reason
is structural rather than numerical.**

The realized operator does not carry the full reflection table; it carries the
table **narrowed to** :math:`\Gamma_+`-local indices (today the crossing function
of ``ScaledOrdinatePermutation`` reads it off the certified pairing), because the
boundary map was narrowed to :math:`\Gamma_+ \to \Gamma_-` in B3.4a. Whether that narrowed
array satisfies ``perm[perm] == arange`` depends on how ``to_local``'s
``searchsorted`` happens to order the locals — which is a property of the
quadrature, not of the mirror. `[M]` for **one** value, a mirror about
``x`` on ``xmin`` (in today's spelling, ``MirrorDeck()`` on ``xmin``):

.. list-table::
   :header-rows: 1
   :widths: 34 12 24 30

   * - quadrature
     - :math:`|\Gamma_+|`
     - full-space table
     - **narrowed** ``perm[perm] == arange``
   * - ``gauss_legendre(4)``
     - 2
     - involution
     - ``True``
   * - ``gauss_legendre(8)``
     - 4
     - involution
     - ``True``
   * - ``product(4, 4)``
     - 4
     - involution
     - ``True``
   * - ``level_symmetric(6)``
     - 24
     - involution
     - ``True``
   * - ``lebedev(17)``
     - 49
     - involution
     - ``False``

The physics does not vary with the quadrature, so the raw index answer
carries no invariant meaning — and Lebedev is only the row that makes
that *visible*. The deeper point holds at all five: the flag's documented
purpose is **self-adjointness in the unweighted inner product**, and
self-adjointness is undefined for a map between two different spaces. The
four ``True`` rows were never right; they were unfalsifiable.

Binding the operator is what converts this from a wrong value into a
**refused question**. Once the ends are :math:`\Gamma_+(f)` and
:math:`\Gamma_-(f)`, the composition :math:`P \circ P` is not an
expression at all, and ``@`` says so at construction:

.. code-block:: text

   P @ P            -> IncompatibleOperatorComposition:
                       OperatorProduct requires each factor's domain to equal the codomain of the
                       factors to its right; ScaledOrdinatePermutation declares domain
                       AngularFaceTraceSpace('angular_trace[xmin:outflow]', shape=(8,)) against
                       AngularFaceTraceSpace('angular_trace[xmin:inflow]', shape=(8,)).
   P_xmin @ P_xmax  -> IncompatibleOperatorComposition   (cross-face: ... against 'angular_trace[xmax:inflow]')
   P @ P.inverse()  -> composes: OperatorProduct (Γ₋ → Γ₊ → Γ₋)

   `[M]` 2026-10-24, ``P`` the specular ``ScaledOrdinatePermutation`` of ``gauss_legendre(16)`` on a slab
   face (built with the helpers of ``tests/gates/sn/architecture/_realized_boundary_helpers.py``),
   the three expressions evaluated, the messages printed.

⭐ So the flag was **retired rather than corrected**, and that is the
load-bearing choice. A refined flag would have answered ``False`` for a
bound cross-space permutation — the right value, still stored, still
obliged to answer for every future caller in every future binding state.
Asking ``P @ P`` instead replaces a value that *can* be wrong with a
composition that *cannot be formed*: the same claim, delivered by the
algebra, with no second clause to keep in step. The involution that IS
real lives one tier up, on the full-space mirror pairing
(the permutation :math:`\pi` that
:func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` returns,
certified before it is returned), where domain and codomain coincide — asserted there
by the ERR-044 involution check (``PairingNotInvolutiveError``). **Two different claims had been sharing one
predicate name across two tiers**, and the fix is to make the tier
explicit rather than to tune the value.

.. note::

   The retirement had zero production consumers and two test assertions,
   both on unbound abstract permutations. Those were *rewired, not
   deleted* — from "the flag matches the index test" to "the square is the
   identity", asked of the algebra
   (``tests/gates/numerics/test_permutation_operator.py``). The behavioural
   content survived the attribute; only the caching of it did not.

One consequence worth stating, because dropping it would have been the
quieter bug. :math:`P : \Gamma_+ \to \Gamma_-` has
:math:`P^{-1} : \Gamma_- \to \Gamma_+`, so
:meth:`~orpheus.numerics.elementary_operators.PermutationOperator.inverse` **inverts
the binding** rather than carrying or discarding it. Discarding it would
have left the return leg composable with anything while the forward leg
was fully typed — an asymmetry no gate on the forward leg could see.

.. _bc-method-realizability:

When a law realizes in a method — the three tiers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A transport **method** is a *discretization of the trace* — a
projection

.. math::

   \Pi \;:\; \Gamma \longrightarrow \Gamma_h

onto whatever that method can represent on the boundary. SN keeps
per-ordinate values on the quadrature; a P1 / diffusion closure keeps
the two half-range moments :math:`(J^+, J^-)`; MoC keeps track angles;
MC keeps nothing at all (it resolves the trace by sampling).

Whether a boundary law **realizes** in a method is therefore not a
question about the law's complexity. It is the question of whether the
naturality square commutes:

.. math::
   :label: bc-realizability-square

   \Pi \circ (R\,G) \;=\; (R_h\,G_h) \circ \Pi ,

.. (vv-status rationale) Structural/definitional: the commuting square that
   DEFINES "this law realizes in this method". It is the framing statement the
   realizer tables' dispatch and refusals instantiate; the verifiable
   content is per-realizer (the diffusion 𝒜-table pinned unit-by-unit in
   ``tests/gates/diffusion/test_diffusion_realizer_table.py`` and the admission rows
   in ``tests/gates/transport/test_admission_matrix.py``, the SN narrowing pinned
   bit-identical in ``tests/gates/sn/operators/test_b3_domain_narrowing.py``), not a
   solver claim of this equation's own.
.. vv-status: bc-realizability-square documented

i.e. *does discretizing the action of the face's factor agree with acting
with the discretized factor?* (On a face exactly one of :math:`G`, :math:`R` is
non-trivial, so :math:`R\,G` is that factor and the square is the square of the
face's own factor.) Three tiers follow, and the middle one is the one
that surprises readers.

**Tier 1 — exact and faithful.** A scalar law :math:`R = \alpha I` for **any**
:math:`\alpha`, not merely :math:`0` or :math:`1`, **when it is complete upstairs**. A
scalar commutes with every linear projection, so :eq:`bc-realizability-square` holds
identically, and :math:`\alpha` is *recoverable downstairs* — the projected law still
carries the number that distinguishes it. On diffusion's scalar trace the members of
tier 1 are the empty network (the vacuum, :math:`\mathcal{A} = 0`,
``Realized(True)``) and the method-only
:class:`~orpheus.diffusion.zero_flux.ZeroFluxCondition` (:math:`\mathcal{A} = -1`,
``Realized(True)``): each is a complete scalar law, with nothing angular for the
projection to lose. A wall of the shared family that is not empty is *not* in tier 1
there, although its image is also :math:`\mathcal{A} I`: the response is angular
upstairs (a permutation or a cosine-weighted average), so the projection loses
something, and diffusion's table states the lower tier (below). The diffusion
realization is the composition of two maps, the unit to its scalar :math:`\mathcal{A}`
and the scalar to the operator :math:`\mathcal{A}\,I`, which is why one scalar read
per unit is the whole realizer (``DIFFUSION_REALIZERS``).

**Tier 2 — exact but NOT faithful.** The square commutes — the
realization is *correct*, nothing is approximated — and the projection
nonetheless **identifies laws that differ upstairs**. At P1, specular
and Lambertian both give :math:`J^- = \alpha J^+`; diffusion simply
cannot tell them apart. The diffusion table's module docstring states it for the response network of a wall:
*"Diffusion annihilates the shape. It sees* :math:`\Pi_{P_1}\colon (s, d) \mapsto s + d`
*on the response network of a wall: the angular redistribution that tells a specular
return from a Lambertian one (a permutation against a cosine-weighted average in SN) is
integrated out of the half-range* :math:`\ell = 0` *moments. Both chains therefore realize
as the scalar WEIGHT times the identity, and a sum of chains is the scalar sum."*
Both chains are ``Realized(faithful=False)``. The distinction is real in transport (SN
realizes them as a permutation against a cosine-weighted average) and vanishes in any
P1-closed method by construction.

This is the same fact as :ref:`the rank-one theorem <bc-factor-roles>`
read from the other side. Where the response destroys directional
information, the composite :math:`R \circ G` comes out the **same
operator** for any admissible measure-preserving :math:`G` — upstairs,
in the continuum. (That is the theorem's honest conclusion; the equality
":math:`R \circ G = R`" this sentence carried until 2026-08-04 does not
type-check once the two half-traces are distinct spaces.) Tier 2 is the *method's*
version of the same collapse: the projection, not the response, is what
destroys the distinction. A reader who conflates "the realization is exact"
with "the realization is a faithful record of the law" will mis-read the
diffusion table's identical-looking rows (the mirror deck, the specular chain and
the diffuse chain, each :math:`\mathcal{A}` times the identity) as a coincidence.
They are not: they are the image of distinct values under a non-injective
:math:`\Pi`.

**Tier 3 — not exact.** The law's action depends on structure the
method does not represent: an anisotropic response kernel below P1 is
the canonical case. Here :eq:`bc-realizability-square` genuinely fails
and no amount of care in the realizer recovers it. The table states it as
``Approximated(names, issues)``, each approximation carried with the issue that tracks it. The
one shipped row is MoC's unit specular wall on a cylinder, which MoC realizes on a different
domain, an equal-area square cell with reflecting flat faces
(``Approximated.of("an equal-area square cell with reflecting flat faces, not the declared
cylinder", 547)``; `[M]` 2026-10-24, ``MOC_REALIZERS.admit(WallCondition.specular(1.0), ...)`` at
the outer face of a cylinder returns it, while the same table refuses weight :math:`0.5`, the white
wall and the vacuum as ``Refused(PLUMBING, ...)``).

So the dividing line is **scalar vs angular**, NOT trivial vs
non-trivial. An :math:`\alpha = 0.37` **scalar** law (a complete scalar
response, as zero flux's :math:`\mathcal{A} = -1` is) is "non-trivial" and
lands in tier 1; a Lambertian wall is "simple" and lands in tier 2, at any
weight.

.. admonition:: Tier 2, measured — the shape of a wall is INVISIBLE to diffusion
   :class: tip

   The retired re-emission closure turned tier 2 from an argument into a
   measurement: ``AlbedoBoundary(0.4)``, ``AlbedoBoundary(0.4, SpecularReturn("x"))``
   and ``AlbedoBoundary(0.4, IsotropicReturn(...))`` all realized to
   ``0.4 · IdentityOperator()`` on the scalar trace (campaign phase B3.4b). The
   successor states it on the shared type, `[M]` 2026-10-24, the diffusion table's
   realization of each value on a Cartesian face, applied to one positive outflow
   probe (the setup of ``tests/gates/diffusion/test_boundary_realizer.py``,
   ``DiffusionMethodSpace.minimal()``):

   .. code-block:: text

      declaration                  repr of the realized operator                                                  out / J⁺
      WallCondition.specular(0.4)  (0.4·<IdentityOperator domain='?' codomain='?' invertible adjointable>)         [0.4]
      WallCondition.white(0.4)     (0.4·<IdentityOperator domain='?' codomain='?' invertible adjointable>)         [0.4]
      MirrorDeck()                 <IdentityOperator domain='?' codomain='?' invertible adjointable>               [1.]
      WallCondition.vacuum()       <ZeroOperator domain='?' codomain='?' adjointable>                              [0.]
      ZeroFluxCondition()          (-1.0·<IdentityOperator domain='?' codomain='?' invertible adjointable>)        [-1.]

   (The last column is ``np.unique(realized.apply(J⁺) / J⁺)``: one value means the
   operator acts as that scalar on the probe.)

   ``specular(0.4)`` and ``white(0.4)`` give equal images (``np.array_equal`` of their
   actions on the probe is true): two distinct values, one image. The shape is not
   *ignored* by the diffusion realizer — it is **annihilated by** :math:`\Pi`, exactly as
   :eq:`bc-realizability-square` permits. The same two values are pairwise
   distinguishable under SN, where one selects a permutation and the other a
   cosine-weighted average. The third retired spelling, the closure-free
   ``AlbedoBoundary(0.4)``, is the one the SN realizer refused, and it has no successor.

   Read that refusal in this frame too: the SN realizer rejected the closure-free
   spelling not because :math:`\alpha\,I` is *wrong* but because it was
   **under-determined** at SN's resolution. Tier 1's guarantee — "a scalar commutes with
   every projection" — is a statement about laws that are *complete* upstairs. A bare
   :math:`\alpha` is not a complete angular law; it is a complete *scalar* one, and the
   tier argument applies only where its own hypothesis holds. The response network
   closes the gap on the shared type by never offering the bare scalar on an angular
   trace: its weights always multiply a chain that names its crossing.

.. _bc-equivariance:

Equivariance — when the deck transformation has a realization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:math:`G` is **method-independent as a geometric object** — the mirror
:math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n` and the
translation :math:`x \mapsto x + L` are ambient facts about the domain,
not about the solver. Its **realization** :math:`G_h` is not. A
discrete :math:`G_h` exists exactly when the discretization is
**equivariant** under :math:`g`, i.e. when :math:`\Pi` intertwines the
two:

.. math::

   \Pi \circ G \;=\; G_h \circ \Pi .

And that condition splits cleanly by **which coordinate** :math:`g`
touches:

* **Specular acts on the ANGULAR coordinate.** So :math:`G_h` exists
  only if the *quadrature* is symmetric under the reflection: there
  must be an index map :math:`\pi` with
  :math:`\Omega_{\pi(n)} = \Omega_n - 2(\Omega_n\cdot\hat n)\hat n` and
  matching weights. That is precisely what the three certification checks of
  :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing` verify, each with its own error
  class: **ERR-042** (``PairingNotMeasurePreservingError``: the pushforward of the
  discrete face measure :math:`m_n = w_n |\mu_{a,n}|` under :math:`\pi` is that measure,
  :math:`m_{\pi(n)} = m_n`), **ERR-044** (``PairingNotInvolutiveError``:
  :math:`\pi \circ \pi = \mathrm{id}`) and **ERR-045** (``PairingDoesNotCrossError``:
  every non-tangential ordinate's partner is non-tangential with the opposite sign on
  the face's axis). A fourth refusal, ``QuadratureNotClosedUnderMotionError``, is the
  case where no weight-preserving bijection exists at all (ERR-073). Read the three
  correctly: they are **discretization-admits-the-symmetry** checks, not physics checks.
  The physics — that specular reflection is an isometry — is not in question; what is
  in question is whether *this quadrature* can express it. All three are failure mode
  **#5 (index error)** in the V&V taxonomy, and their independence is **measured, not
  assumed**: the GL-8 neighbour-pair table passes involution while redding the measure,
  and the identity table passes involution *and* measure while redding only the
  inflow → outflow check (``test_the_three_errors_keep_their_teeth`` in
  ``tests/gates/numerics/test_ordinate_pairing.py``: each broken table raises exactly its own error class). The function takes the face's
  mirror motion and no axis argument; the axis is the one the motion names
  (``test_ordinate_pairing_takes_no_axis_argument``). The certification fires at the table's
  door: the SN realizer table asks ``ordinate_pairing`` when it admits a ``SpecularChain`` or
  a ``MirrorDeck`` and answers ``Refused(ANGULAR_RESOLUTION)`` when the quadrature is not
  closed under the mirror, so every :class:`~orpheus.sn.problem.SNProblem` construction
  certifies them (the problem realizes its boundary through the shared
  ``realize_boundary`` body, which admits first).
* **A spatial wrap acts only on the SPATIAL coordinate.** Ordinate
  :math:`n` at face :math:`f'` feeds ordinate :math:`n` at face
  :math:`f`, untouched — which is why the wrap's pairing
  (``pairing_of(TranslationDeck(), face, dimension=d)``) answers
  ``permutes_ordinates = False``, DERIVED from its motion's identity
  linear part (``not motion.is_translation``) rather than declared. **Every angular
  discretization is therefore trivially equivariant under it.** Periodic is the more
  method-agnostic of the two deck transformations, and for a sharper
  reason than "it is a trace connection": there is no angular symmetry
  requirement to fail.

.. _bc-refusal-axes:

Three independent axes of method-dependence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every refusal of a realizer table is a ``Refused(axis, reason, issue)`` whose axis is
one of three **independent** axes, or the fourth thing that is not an axis, plumbing
(``RefusalAxis``: ``ANGULAR_RESOLUTION``, ``TOPOLOGY``, ``STATE_CONE`` and ``PLUMBING``).
A method can sit anywhere on each of the three, and reading them as one axis ("this method
is too coarse") is what kept the taxonomy invisible. The table that realizes a declaration
is the table that refuses it, so the refusal and its axis are read from the one door
(``RealizerTable.admit``) and cannot disagree with the realization.

.. list-table:: The three axes of method-dependence
   :header-rows: 1
   :widths: 20 38 42

   * - axis
     - question
     - the refusal that shows it
   * - **angular resolution**
     - can the method represent :math:`R`'s angular structure, or the symmetry
       :math:`G` acts by?
     - ``SN_REALIZERS`` refuses a ``SpecularChain`` or a ``MirrorDeck`` when
       ``ordinate_pairing`` cannot certify the quadrature as closed under the face's
       mirror, ``Refused(ANGULAR_RESOLUTION)``
       (``test_a_certification_failure_names_its_axis_at_the_table``); tier 3
       above is the other member, an anisotropic :math:`R` below P1, which no shipped
       wall is.
   * - **spatial / topological**
     - can the method's operator express **cross-face** coupling?
     - ``DIFFUSION_REALIZERS`` **refuses a** ``TranslationDeck`` — "the wrap couples
       'xmin' to its opposite face 'xmax': a per-face scalar has no slot for cross-face
       coupling", ``Refused(TOPOLOGY)``, the one geometry P1 cannot integrate away into a
       per-face albedo :math:`J^- = \mathcal{A} J^+`. Its codomain is a per-face scalar,
       a block-diagonal object with no slot for a face *pair*. **Nothing to do with
       angle**: as above, a translation is trivially equivariant for every angular
       discretization. A method could resolve angle perfectly and still refuse this. MC
       refuses a specular wall on a cylinder on this axis for the converse reason (it
       realizes a periodic square lattice cell, whose faces ARE coupled by an image
       lattice, not the declared cylinder, #547).
   * - **state-cone / sign**
     - is the value representable in the method's **state cone**?
     - ``SN_REALIZERS`` **refuses diffusion's** ``ZeroFluxCondition``:
       :math:`\mathcal{A} = -1` needs a *signed* current, and the SN state is an angular
       flux with :math:`\psi \ge 0`, which admits no negative angular inflow
       (:math:`\psi \ge 0 \Rightarrow J^\pm \ge 0`). The table's own ``StateCone`` is
       ``NONNEGATIVE_FLUX`` and the unit's return fraction is negative, so the unregistered-unit door
       answers ``Refused(STATE_CONE)`` without SN importing the diffusion type. Diffusion realizes the
       very same value without difficulty, as ``-IdentityOperator()`` (the one-term
       combination :math:`(-1)\,I`, ``DiffusionZeroFlux``).

Note the shape of that table: on the topological axis diffusion refuses
what SN accepts, and on the state-cone axis SN refuses what diffusion
accepts. Neither method is "coarser". They are incomparable. `[M]` 2026-10-24, the
admissions on a Cartesian ``xmin`` face read through ``SN_REALIZERS.admit`` and
``DIFFUSION_REALIZERS.admit`` (SN, diffusion): ``TranslationDeck()``
``Realized(True)`` and ``Refused(TOPOLOGY)``; ``ZeroFluxCondition()``
``Refused(STATE_CONE)`` and ``Realized(True)``; ``MirrorDeck()`` ``Realized(True)`` and
``Realized(False)``; ``WallCondition.specular(0.4)`` ``Realized(True)`` and
``Realized(False)``. The two refusals of the pair are each read beside a ``Realized``
answer of the other table for the same value, which is the control that the refusal is
the method's and not the value's.

:math:`q` **is a fourth thing, and not on any of these axes.** It is a
**vector in** :math:`\Gamma_-`, not an operator, so the only question
it asks of a method is whether :math:`\Gamma_-` is represented at the
fidelity the source demands. Diffusion's refusal of any nonzero source
(``WallCondition.uniform`` or ``.framed``, ``DiffusionSource``) is therefore a
**plumbing** refusal, not a representability one — :math:`\Gamma_-` on
a scalar trace is one number per face, which the trace carries
perfectly well. The table says so in place: ``Refused(PLUMBING, #656)``
(`[M]` 2026-10-24: ``DIFFUSION_REALIZERS.admit(WallCondition.uniform(1.0), ...)`` on a
Cartesian face returns it, while ``SN_REALIZERS.admit`` of the same value is
``Realized(True)``), and the refusal disappears the day the diffusion fixed-source arm is
wired, **with no theory changing**. Stated explicitly so a future reader does
not mis-file it beside the three structural refusals above. Plumbing is the axis of
every refusal that has no theory behind it, and ``Refused`` itself refuses to be built as
``PLUMBING`` without the issue that tracks the missing machinery, so "not built yet"
cannot be an unexplained refusal: the hollow radial body's inner face (#511), the
radial characteristic's corner action for a diffuse return or a source (#653) and the
missing registry entries of the CP, MoC and MC tables (#654) are the others, read from the
same door (`[M]` 2026-10-24: ``SN_REALIZERS.admit(WallCondition.white(1.0), ...)`` on a cylinder's outer
face is ``Refused(PLUMBING, #653)``, and on its inner face ``white(0.5)`` is
``Refused(PLUMBING, #511)``).


Layer 1 — trace structure
-------------------------

The trace operators :math:`\gamma_\pm` carry their domain
information on **one** typed
:class:`~orpheus.numerics.space.FunctionSpace` subclass,
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`, which stores
the whole boundary :math:`\Gamma = \partial\Omega \times S^d` once and
exposes inflow / outflow as two directional **selectors** over it:

* :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`
  selects :math:`\Gamma_- = \{(\mathbf{r}, \Omega) \in \partial\Omega
  \times S^2 : \Omega \cdot \hat n(\mathbf{r}) < 0\}` — the per-face
  directional half of the boundary on which the incoming angular flux
  is constrained by the boundary condition.
* :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`
  selects :math:`\Gamma_+` symmetrically — the boundary half on which
  the outgoing flux is *not* constrained by the BC but is *consumed* by
  it (as :math:`\gamma_+ \psi`).

.. note:: **One space, two selectors (Issues #205 / #201).** The
   pre-#188 design carried two separate typed spaces,
   ``InflowTraceSpace`` and ``OutflowTraceSpace``, one per direction.
   The View-G field-vocabulary refactor (#205 / #201) collapsed them
   into the single :class:`AngularTraceSpace
   <orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace>` on the observation
   that **inflow and outflow are operations on one space, not two
   spaces**: whether an :term:`ordinate` is incoming or outgoing at a face is a
   *predicate* — :math:`\mathrm{sign}(\Omega \cdot \hat n_f)` —
   evaluated against the same signed-projection data, not a property of
   the space's identity. :class:`AngularTraceSpace` stores the signed
   projection :math:`\Omega \cdot \hat n_f` once per face; the two
   ``*_indices_for_face`` methods are selectors over it (see
   :ref:`bc-trace-structure`).

The signed-projection table is read **once**.
:func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`
renders :math:`\Omega\cdot\hat n_f` from the quadrature and the face
names, and
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout`
derives from it every half-space of every face, so no consumer holds a
second classifier of "incoming" and "outgoing". A realized face
operator of :ref:`bc-realizer-layer` names its two half-traces by
reading them off that ladder: its domain :math:`\Gamma_+` is
``trace.outflow_space(face)`` and its codomain :math:`\Gamma_-` is
``trace.inflow_space(face)``, and each type holds its face, the trace and, where it has one, its
amplitude and partner (``RankOneReemission`` holds ``(alpha, face,
trace)``, ``ZeroMap`` ``(face, trace)``, ``ScaledOrdinatePermutation``
``(alpha, permutation, face, partner, trace)``), and no axis or outward
sign of its own
(``orpheus/sn/boundary/realized_operators.py``). Of the two
certifications that sat in the realizer, ERR-041 (vacuum fills the
incoming trace only) moved to where the data lives and is witnessed on
the ladder, in
``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``,
and ERR-047 (a source lies on :math:`\Gamma_-`) remains the one SN-side
check, :meth:`~orpheus.sn.boundary.realizer.SNSource.sample`, which
refuses a space that is not an inflow half-trace; the face's own
:math:`\Gamma_-` is what ``from_method`` hands it.


.. _bc-domain-narrowing:

The trace maps as operators — the narrowed domain
-------------------------------------------------

Until campaign phase **B3.2** the trace operators :math:`\gamma_\pm`
existed on this page and in the affine form :eq:`affine-bc-form` but
had **no type in the code**. They were spelled three different ways and
typed as none of them, and the SN boundary law was consequently handed
the *whole face slot* — all ``quad.n_ordinates`` ordinate rows — with the outflow
rows thrown away afterwards by a slice-write at the consumer. B3.1
minted the missing primitive and B3.2 narrowed the law onto it.

:class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator` is the
gather / scatter pair. Given a **sorted, unique** index set
:math:`S \subset \{0,\dots,N-1\}` of size :math:`m < N`:

.. math::
   :label: bc-trace-restriction-pair

   \gamma_S : \mathbb{R}^N \to \mathbb{R}^m,\quad
   (\gamma_S x)_j = x_{S(j)}
   \qquad\text{and}\qquad
   \iota_S = \gamma_S^{\mathsf T} : \mathbb{R}^m \to \mathbb{R}^N,\quad
   (\iota_S y)_i =
   \begin{cases} y_j & i = S(j) \\ 0 & i \notin S \end{cases}

.. (vv-status rationale) Definitional: the gather/scatter pair that types the
   affine form's γ±. Its DEFINING LAWS (γι = I, ιγ idempotent and symmetric,
   ι materialised against the dense γᵀ, γ₋∘ι₊ = 0, and the three-way partition
   resolving I) are pinned by the nine foundation tests opening
   ``tests/gates/numerics/test_trace_restriction_operator.py``; this equation states
   the definition those tests verify, not a solver claim.
.. vv-status: bc-trace-restriction-pair documented

with :math:`\gamma_S \iota_S = I` on the restricted space and
:math:`\iota_S \gamma_S = P_S` the orthogonal projector onto it. It is
a **sibling of** :class:`~orpheus.numerics.elementary_operators.PermutationOperator`
and deliberately **not a subclass**: same ``np.take`` mechanism with a
non-square index array, but different algebra in kind — a permutation
is a bijection (invertible, with an algebra-closed
:meth:`~orpheus.numerics.elementary_operators.PermutationOperator.inverse`) while a
restriction is rank-deficient by construction and has a *scatter*
transpose rather than an inverse. Inheriting would promise what the
type cannot honour.

The two named instances live on the trace space, cached per face
(:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_restriction`
:math:`= \gamma_+` and
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_restriction`
:math:`= \gamma_-`), and the SN face action is their composition with
the face's realized response :math:`R_f` (the operator :math:`\Gamma_+ \to \Gamma_-`
of :ref:`bc-realizer-layer`):

.. math::
   :label: bc-face-action-narrowed

   B_{\rm face} \;=\; \iota_- \circ R_f \circ \gamma_+ ,
   \qquad
   B_{\rm face}^{\mathsf T}
   \;=\; \iota_+ \circ R_f^{\mathsf T} \circ \gamma_- .

.. (vv-status rationale) Structural: the composition ``_reflect_trace`` spells
   for every face. Its verifiable content is the B3.2 bit-identity gate —
   the composition reproduces the retired full-face-then-slice expression
   exactly (``np.array_equal``, against a numpy reference materialised off the
   DECLARED kind, over slab-asym / slab-sym / sphere / cyl ``product(2,4)`` /
   2-D Cartesian LS4), in ``tests/gates/sn/operators/test_b3_domain_narrowing.py``.
   Not an independent solver claim.
.. vv-status: bc-face-action-narrowed documented

**Nothing is computed and then discarded.** The rows the pre-B3.2
slice-write dropped are simply not in the operator's domain, so a
non-zero outflow emission — which would corrupt the outflow-definition
residual :math:`r_{\rm outflow}`, a quantity that carries no :math:`B`
term at all — is now **unrepresentable** rather than projected away.
That is the ``coding-elegance`` Pattern-4 form of what
:ref:`bc-extraction-design-corrections` correction 2 previously
achieved by projection.

.. _bc-narrowing-what-it-removed:

What the narrowing removed — three spellings, one pair
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every restriction-shaped expression the boundary review catalogued is a
composition of :eq:`bc-trace-restriction-pair`, not a primitive of its
own:

.. list-table::
   :header-rows: 1
   :widths: 46 24 30

   * - spelling found in the subsystem
     - is
     - status
   * - the slice-write ``out[sel] = full[sel]`` in ``_reflect_trace``
     - :math:`\iota_- \circ \gamma_-` (:math:`= P_{\rm in}`)
     - **dissolved** at B3.2 — the realized response's codomain *is*
       :math:`\Gamma_-`, so :math:`\iota_-` is the honest scatter, not
       a projection of a wider image
   * - ``IncomingSourceOperator``'s dense inflow-mask multiply
     - :math:`\iota_- \circ \gamma_-` (measured bit-identical to the
       slice-write)
     - **dissolved** at B3.4a — the operator was thereafter asked to
       fill :math:`|\Gamma_-|` rows directly, so nothing off
       :math:`\Gamma_-` was left to erase. That operator itself retired
       at **P3**, and the declared inflow is now sampled at the same
       :math:`|\Gamma_-|` shape one tier out, by
       :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.from_method`
       through ``SNSource.sample`` — see :ref:`bc-affine-source-channel`
   * - ``IncomingOrdinateMaskTensor`` (**retired**)
     - :math:`I - \iota_- \circ \gamma_-`
     - taken **off the vacuum path** at B3.2 and **deleted at B3.3**.
       Of its thirteen tests, twelve asserted either the masking
       semantics the narrowing removed or a law
       :class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator`'s
       own battery already gates; only the non-aliasing claim was
       genuinely uncovered, and it migrated

and ``P_in ∘ P_out = 0`` stops being a curiosity: it is
:math:`\gamma_- \circ \iota_+ = 0`, true because two disjoint index
sets have nothing to hand each other.

.. warning::

   **Two traps the narrowing exposes, both measured.**

   1. **The index remap is** ``searchsorted``, **not** ``arange``.
      Mapping a subset of *global* rows into positions inside a
      restricted space is
      :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace.to_local`
      (owned by the half-trace SPACE since G6.5 — the embedding data,
      ``ordinate_indices``, is the space's), which needs a sorted index
      set — which is why sortedness is a **construction guard** on the
      space, not tidiness. The naive
      ``arange(sel.size)`` is right only when the subset is a *prefix*
      of the index set, and the two sites where the remap appears are
      discriminated by **different, complementary fixtures**: the
      reflective narrowing (``perm[inflow]`` into :math:`\Gamma_+`) is
      wrong on the **slab**, where the mirror reverses order —
      ``gauss_legendre(4)`` at ``xmax`` gives ``perm[inflow] = [3, 2]``
      → local ``[1, 0]`` where ``arange`` says ``[0, 1]`` — and right
      on the cylinder; the schedule split (a rows-subset into
      :math:`\Gamma_-`) is wrong in **2-D**, where the lower-half rows
      are not a prefix, and right in 1-D. **1-D coverage is not
      sufficient and neither is 2-D.**
   2. **A shape assertion cannot detect a mis-typed domain.** Measured:
      :math:`|\Gamma_+| = |\Gamma_-|` on each of the 46 (rule, face) pairs
      measured (9 rules; the census of the signed-projection gotcha
      below), so an un-narrowed endomorphism :math:`\Gamma_+ \to \Gamma_+`
      has *exactly the right output shape*. This is a textbook vv
      **Mode 12** blindness — the measured functional (shape) has the
      error class in its invariance group — and it is what let three
      un-narrowed realizer arms survive B3.2's first pass. Only the
      anti-Mode-12 leg, which interrogates the emitted operator's
      declared spaces rather than its output shape, found them.

This passage records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

⛔ **This passage, down to the partial-constructor paragraph, records
the state between B3.4a and B3.4b, and both of its deferrals have since
landed as it predicted.** B3.4b (``943b37c1``, 2026-08-01) gave albedo
the explicit re-emission closure carried in :math:`R` that the passage
calls for, and made the closure-free spelling a refusal; B3.4c
(``79b5affe``, the same day) built periodic's partner-face channel, and
G6.3 step 7 (``c3bb7341``, 2026-08-07) derived its arrow
:math:`\Gamma_+(f') \to \Gamma_-(f)`. Every operator the realizer table
builds is narrowed (domain :math:`\Gamma_+`, codomain :math:`\Gamma_-`),
and the strict xfails that pinned the deferred rows are retired
(``tests/gates/sn/operators/test_b3_domain_narrowing.py``, ``[M]``
2026-10-24: 27 of 27 pass under ``python -O``, ``python -O -m pytest
tests/gates/sn/operators/test_b3_domain_narrowing.py -q``, 471c4737). The
first consequence below is moot for the same reason (a closure-free
albedo cannot be built, so ``0.3·specular + 0.7·albedo`` cannot be
spelled), and the second is re-measured in place. The design analysis is
kept because it is why the two repairs took the shape they did.

**Four of the seven laws were narrowed at B3.4a; two remained.** B3.2
narrowed the two laws SN reaches from a Problem — ``vacuum`` and
``reflective`` — and measured the remainder then as *six* realizer rows
across four law kinds. **B3.4a** took two of those kinds:

* ``white``, whose Lambertian kernel now contracts over
  :math:`\Gamma_+` and re-emits on :math:`\Gamma_-`
  (:ref:`bc-narrowing-b34a`); and
* ``prescribed_inflow``, whose rank-0 source was thereafter asked for
  :math:`|\Gamma_-|` rows directly. **P3** then collapsed that arm
  further — onto the zero morphism vacuum already returned, with the
  source itself moved to the boundary-source channel
  (:ref:`bc-affine-source-channel`), so the narrowing survives while the
  operator that carried it does not.

What remained was **four rows across two law kinds** — ``albedo`` at
*three* rows (:math:`\alpha = 0` and :math:`\alpha = 1` take fast paths
returning a bare :class:`~orpheus.numerics.operator.ZeroOperator` /
:class:`~orpheus.numerics.operator.IdentityOperator`, which are
**endomorphisms**, and :math:`\alpha \notin \{0,1\}` scales an
identity), plus ``periodic``. Both are blocked on a **design ruling**,
not on plumbing: albedo is under-determined on an angular trace
(:math:`R = \alpha\,I` is a :math:`\Gamma_+ \to \Gamma_+` endomorphism
and its :math:`G = \mathrm{id}` supplies no crossing), so
**B3.4b** must give it an explicit re-emission closure carried in
:math:`R`; periodic's :math:`G` reads the PARTNER face's
:math:`\Gamma_+`, which **B3.4c** was to build (#183, #189). Measured on
both then: they silently accepted a :math:`\Gamma_+` input and echoed it
back — i.e. :math:`\Gamma_+ \to \Gamma_+`, the wrong codomain, invisible
to a shape check (vv Mode 12 again). All four were unreachable in
production — the SN registry admitted only ``{vacuum, reflective}`` — so
the tree was green; they were pinned by strict xfails, each
``--runxfail``-verified to red for *its own* documented reason.

Two consequences to carry meanwhile. First, a narrowed law **cannot
honestly compose with an un-narrowed one**: the sum of a
:math:`\Gamma_+ \to \Gamma_-` leaf and a :math:`\Gamma_+ \to \Gamma_+`
one has no well-typed codomain. Note that this **no longer announces
itself as a raise** — since B3.4a, ``0.3·specular + 0.7·white`` is a
sum of two narrowed leaves and is simply correct, while
``0.3·specular + 0.7·albedo`` *runs* and returns
:math:`|\Gamma_-|`-shaped output that is silently wrong, because
:math:`|\Gamma_+| = |\Gamma_-|` swallows the mismatch. The Mode-12
lesson applies to the *algebra* as well as to the leaves.

Second, ``SNMethodSpace.minimal`` became a **partial constructor**: a
quadrature alone cannot name a *face's* :math:`\Gamma_+`, so it stopped
sufficing for each law as that law was narrowed. At B3.4a only
``albedo`` and ``periodic`` still realized from it, precisely because
they had not been narrowed. By 2026-09-22 it realized **none**: ``[M]``
2026-09-22 (``scratch/_definitions/round2/minimal_probe.py``), 0 of 10
law spellings realized on ``SNMethodSpace.minimal(gauss_legendre(8))``,
the eight narrowed ones refusing for want of the face's inflow or
outflow indices (periodic for want of a face), and so did
``zero_flux`` and the closure-free ``albedo``, which a face method
space refused as well, while the same ten on a face method space
realized 8. No production code called it (0 of 352 modules under
``orpheus/``); 6 test modules did. The constructor retired in WF2.4
(``c2a14005``): the SN method space is the pair ``(quadrature, trace)``,
the face is an argument of the realization, and face orientation is a
structural demand of realization, not an implementation detail.


.. _bc-narrowing-b34a:

B3.4a — narrowing white and prescribed inflow, and what dissolved
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

The two laws B3.4a narrowed are the campaign's cleanest illustration of
its thesis: **a too-wide codomain does not merely cost rows, it buys
bugs, and narrowing does not fix them — it makes them unspellable.**
Each law arrived carrying a correction for a defect that the narrowing
then removed the possibility of, so in both cases the shipped diff
*deletes* the correction rather than repairing it.

**White — the outflow classifier that had a twin.** Pre-B3.4a the
then-welded ``AngularAverageOperator`` was a
full-face endomorphism whose ``cos_w`` carried :math:`N` entries zeroed
off the outgoing hemisphere by its OWN test,
``(outward_sign * mu_n) > 0.0``. That is a **second outflow
classifier**, and it disagreed with the trace space's wherever a
quadrature carries tangential ordinates, because the trace space
classifies against ``TANGENTIAL_EPS`` while a strict ``> 0.0`` compare
does not. Measured at ``xmax``:

.. list-table:: Tangential ordinates, and where the two classifiers diverge
   :header-rows: 1
   :widths: 30 14 22 34

   * - Quadrature
     - :math:`N`
     - tangential at a face
     - rows the ``> 0.0`` test claims that :math:`\Gamma_+` does not
   * - ``gauss_legendre(8)``
     - 8
     - 0
     - 0 — which is why the twin never surfaced in slab tests
   * - ``product(2, 4)``
     - 8
     - **4**
     - **2** — the two whose :math:`\Omega\cdot\hat n` round-off is
       *positive*; the composite diverged by ``6.1e-05`` on a cylinder
   * - ``lebedev(9)``, ``lebedev(17)``
     - 38, 110
     - 12 on *every* face
     - **0 on every face** — the band's round-off falls on the
       non-positive side, so the two classifiers happen to agree
   * - ``level_symmetric(4)``, ``level_symmetric(6)``
     - 24, 48
     - 0
     - 0

Read that last column carefully: a quadrature can carry many tangential
ordinates and still expose no disagreement, because what matters is not
how many sit in the band but **which side of zero their round-off falls
on**. Measured over the whole production inventory
(``gauss_legendre`` 4/8, ``product`` 2×4 / 3×4 / 4×8, ``lebedev`` 9/17,
``level_symmetric`` 4/6) across all six face names, the disagreement
occurs **only for the** ``product`` **family, and there only on**
``xmax`` / ``xmin`` / ``ymax`` — ``ymin`` carries the same tangential
count with zero mis-admissions, because the sign flip moves the
round-off across zero. ``lebedev`` has twelve tangential ordinates per
face and mis-admits none anywhere; ``level_symmetric`` has none at all.

Two consequences worth carrying. A tangential-count audit is **not** a
sufficient screen for this bug class — only using one classifier is.
And the exposure is *face-asymmetric within a single quadrature*, so a
fixture that exercises one face of a ``product`` rule can be green
while its opposite face is wrong.

B3.4a's fix was to route the classification through the codebase's
**single** face-name :math:`\to` signed-projection primitive:
``(axis, outward_sign)`` is rendered as the face NAME and handed to
:func:`~orpheus.numerics.spaces.angular_trace_space.build_omega_dot_n`,
which classifies against ``TANGENTIAL_EPS`` exactly as the trace space
does. It landed on the welded operator's own ``from_quadrature``
constructor, and after **G6.3 step 3b** retired that operator the same
lines lived one layer out, in the realizer's ``_checked_angular_average``
factory, until wave 2 replaced that factory by the typed
:class:`~orpheus.sn.boundary.realized_operators.RankOneReemission`,
whose hemisphere is ``trace.outflow_space(face)`` read off the ladder.
With the domain narrowed the kernel classifies **nothing**: it
is *handed* :math:`\Gamma_+`. The twin is not fixed, it is unspellable,
and every ``cos_w`` entry is strictly positive by construction rather
than mostly zero — which is why the guard tightened from ``>= 0`` to
``> 0``, and why it is
:class:`~orpheus.sn.boundary.angular.PartialCurrentOperator`'s
constructor that today refuses a full-face weight vector by name.

**Prescribed inflow — the mask that dissolved.** An
``InflowSourceSpec`` filled whatever
block *shape* it was handed; it carried no trace knowledge. So pre-B3.4a
``IncomingSourceOperator`` emitted a
full-face block and then zeroed every non-inflow ordinate — outflow AND
tangential — to make :math:`q \in \Gamma_-` hold. With the codomain
narrowed it simply asked the spec to fill
``(|Γ₋|,) + psi_out.shape[1:]``:
the rows the mask used to zero were not in the codomain to be emitted
on, so **ERR-047 was closed by the TYPE rather than by an erasure**.
That is precisely the dissolution vacuum's projector underwent at B3.2
— the mask was never load-bearing physics, only the cost of a codomain
that was too big. Its companion, an *unmasked* fallback branch legal
only for :math:`q \equiv 0`, retired with it: post-B3.2 every
realization needs :math:`\gamma_+` too, so a method space with no face
data could not reach the operator at all and the branch was already
unreachable on the realize path.

.. note::

   **That operator no longer exists** (campaign phase **P3**,
   2026-08-05), so this paragraph is history — but the *typing* half of
   it survives verbatim one tier out. The inflow is still sampled at
   exactly ``(|Γ₋|,) + trailing``, now by
   :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.from_method`
   (through ``SNSource.sample`` and ``sample_terms``) on the way into the
   boundary-source channel, and ERR-047 is still closed by the type rather
   than by an erasure. What B3.4a could not
   see is that the operator was ALSO affine — it returned :math:`q`
   whatever it was handed — and that this put an affine term inside the
   linear :math:`B` block, which the missing ``BlockRole.BOUNDARY``
   stamp did **not** prevent. See :ref:`bc-affine-source-channel` for
   the collapse, the measurement, and why vacuum and prescribed inflow
   realize to the same object.

**Equivalence — and the reason two DIFFERENT effects both measure at
1 ULP.** White's operator was measured against a reconstruction of the
pre-B3.4a body over an :math:`\mathcal{O}(1)` random probe, six seeds:

.. list-table:: Old-vs-new, and which effect produces the difference
   :header-rows: 1
   :widths: 26 10 16 22 26

   * - Quadrature
     - face
     - rows the ``> 0.0`` test mis-admits
     - old vs new on an :math:`\mathcal{O}(1)` probe
     - cause
   * - ``gauss_legendre(8)``
     - ``xmax``
     - 0
     - **bit-identical**, every seed
     - — (padding with exact zeros is bit-neutral here)
   * - ``level_symmetric(6)``
     - ``xmax``
     - 0
     - **bit-identical**, every seed
     - —
   * - ``level_symmetric(6)``
     - ``ymax``
     - 0
     - ``1.11e-16`` / ``5.55e-17`` (:math:`\le` 1 ULP)
     - **reduction order only**
   * - ``lebedev(17)``
     - ``xmax``, ``ymax``
     - 0
     - ``1.11e-16`` / ``5.55e-17`` (:math:`\le` 1 ULP)
     - **reduction order only**
   * - ``product(2, 4)``
     - ``xmax``, ``ymax``
     - **2**
     - :math:`\le` 1 ULP — *and this is the trap*
     - **the classifier twin**, masquerading as noise

Two mechanically different things are being measured, and on a
well-scaled probe they are **indistinguishable**.

Where the classifiers agree (every quadrature but ``product``), the
only change is that the sum now runs over the restricted
:math:`|\Gamma_+|`-entry array instead of a zero-padded
full-:math:`N` one, so the floating-point reduction ORDER changed and
addition is not associative. That is *not* an error and *not* a
regression — it is the ``vv-principles`` **principled-equivalence**
case: a named intermediate (the cosine-weighted current over
:math:`\Gamma_+`) replacing an unnamed zero-padded one, drift bounded
by reduction depth :math:`\times` ULP.

Where the classifiers **disagree**, the change is a genuine VALUE
correction that merely *looks* like noise. At ``product(2, 4)``,
``xmax``, the four tangential ordinates carry
:math:`\Omega\cdot\hat n` of :math:`+5.0\times10^{-17}` (twice) and
:math:`-1.5\times10^{-16}` (twice) — round-off, not exact zero — so a
strict ``> 0.0`` test admits the two positive ones into
:math:`\Gamma_+` while the trace space calls all four tangential. The
resulting spurious weights are :math:`7.85\times10^{-17}` against a
normalization of :math:`2.5651`, so the **denominator is unchanged to
the last bit** (measured :math:`\Delta\text{norm} = 0.0` exactly) and
the whole discrepancy lives in the NUMERATOR — where it is
:math:`\psi`-weighted. It therefore scales with the flux carried on
the mis-admitted rows and is **not bounded by floating point**:
measured **6.1e-05** relative when those rows carry :math:`10^{12}`
times the genuine outflow flux, which is the regime the cylinder
composite hit.

.. warning::

   **Do not read the 1-ULP row for** ``product(2, 4)`` **as evidence
   of equivalence.** An :math:`\mathcal{O}(1)`-scaled probe cannot
   separate a reduction-order artefact from a mis-classified-ordinate
   bug, because the offending weight is itself :math:`\mathcal{O}(\epsilon)`
   — the error only becomes visible once the *flux ratio* between the
   tangential and outflow rows is large. This is why the twin survived
   in-tree: every equivalence measurement was taken on a well-scaled
   probe, and ``gauss_legendre`` (the slab default) carries no
   tangential ordinates at all, so it could not see the disagreement
   under any probe. The narrowing is justified structurally — one
   classifier, not two — not by the ULP table.

**And one guard was added, not removed, and wave 2 then made it
unspellable.** ``WhiteBoundary`` declared its own ``axis`` /
``outward_sign`` while the method space independently named the face —
two encodings of ONE orientation, and until B3.4a nothing compared them.
A white law declared for ``x`` and installed on ``ymin`` averaged over
the wrong hemisphere and reported nothing. B3.4a's helper
``_checked_angular_average`` cross-checked them and raised a
``BoundaryError`` with ``law="white"`` on a mismatch. It compared index
**SETS**, not sizes, because :math:`|\Gamma_+| = |\Gamma_-|` on every
quadrature in the tree would have made a size comparison Mode-12 blind.
Wave 2 removed the second encoding rather than keep the guard: a
``DiffuseChain`` carries no axis and no outward sign, the realized
:class:`~orpheus.sn.boundary.realized_operators.RankOneReemission` holds
``(alpha, face, trace)``, and its hemisphere is
``trace.outflow_space(face)``, read from the one face name, so a
declaration and a method space can no longer disagree about the
orientation. The cross-check retired with the index sets (WF2.4, the
realizer module's "Retired with the isinstance ladder").

**Two things B3.4a deliberately did NOT do.** It did not type the
Lambertian kernel as the rank-one :math:`u \otimes v` it now visibly is
(with :math:`u = \mathbf{1}_{\Gamma_-}` and :math:`v =
\cos w / \mathrm{norm}` on :math:`\Gamma_+`) — deferred to phase
**B5** as what would make its adjoint structurally available, carrying
the Euclidean-vs-cosine-metric transpose question with it. In the event
that deferral was discharged at **G6.3 step 3b**, and by a different
mechanism: the kernel was *factored* into
:math:`B \circ C` rather than typed as :math:`u \otimes v`, which
dissolved the transpose question instead of answering it
(:ref:`bc-response-adjoint`). And it did
not give the narrowed laws a domain **validator**: white happened to
refuse a full-face input (because
:meth:`PartialCurrentOperator.apply
<orpheus.sn.boundary.angular.PartialCurrentOperator.apply>`, the
chain's first link, checks
``psi.shape[0]``), while prescribed inflow ignored its input entirely
and so had nothing to validate against. Refusal and non-endomorphism
are separate properties; see the warning at :ref:`bc-worked-example`.


.. _bc-affine-source-channel:

The affine source channel — which tier carries :math:`q`
--------------------------------------------------------

:eq:`affine-bc-form` is an **affine** map, and an affine map has two
halves. The realizer of :ref:`bc-realizer-layer` realizes the
**linear** half and only the linear half; the inhomogeneous half
travels a separate, typed channel and is added by the solve. Campaign
phase **P3** (2026-08-05) is the carve that made that sentence true of
the code as well as of the design; before it, one law realized both
halves into a single object that was declared linear and was not. The
wave-2 carve (2026-10-24) then made the split the data model itself: the
exterior of a face is the pair ``WallCondition(network, inflow)``, which
**is** the affine pair :math:`(R, q)`, so the linear half is the realized
``network``, the source half is the ``inflow``, and no object holds both.

Write the affine trace condition with the linear factor named:

.. math::
   :label: bc-affine-linear-factor

   \gamma_-\psi\big|_f
   \;=\; L_f\;\gamma_+\psi\big|_{f'}
   \;+\; q_f ,
   \qquad
   \operatorname{realize}(\text{declaration}_f) \;=\; L_f ,
   \qquad
   q_f \in \Gamma_-(f).

.. (vv-status rationale) Structural/tier-assignment: the same affine form as
   :eq:`affine-bc-form`, with the realizer's output NAMED as its linear factor
   L. It is a statement about which tier carries which half, not a solver
   claim; its verifiable content is the linearity of every realized leaf
   (``B(0) = 0``, ``B(2x) = 2B(x)``) plus the source-tier claim
   :eq:`bc-single-delivery`. Per this page's V&V-status note, its equations are
   definitional / structural-architecture statements.
.. vv-status: bc-affine-linear-factor documented

Here :math:`L_f = R_f` and :math:`f' = f` for a wall (a
:class:`~orpheus.transport.boundary_conditions.WallCondition`), while for
a deck :math:`L_f` is the Koopman operator of the deck's isometry (unit
weight, no amplitude) with :math:`f'` the partner face and
:math:`q_f = 0`. Then every declaration the SN realizer table admits
splits cleanly, and the split is the whole content of the carve.
``ZeroFluxCondition``, diffusion's method-only condition, is REFUSED:
its return fraction is :math:`-1`, outside the non-negative cone of an
angular method, so it meets the table's unregistered-unit door on the
state-cone axis (``[M]`` ``scratch/wf25_bcB_measurements.py`` section S8:
``BoundaryRefused``, ``RefusalAxis.STATE_CONE``; a negative angular
inflow is unrepresentable, so an angular method says ``vacuum``
instead; zero flux is realizable only on a scalar trace).

.. list-table:: Every SN-realizable declaration as a linear factor plus a source
   :header-rows: 1
   :widths: 28 30 14 28

   * - declaration
     - :math:`L = R_f` (wall) or deck operator
     - :math:`q`
     - realized leaf
   * - :meth:`WallCondition.vacuum()
       <orpheus.transport.boundary_conditions.WallCondition.vacuum>`
     - :math:`0` — the zero morphism :math:`\Gamma_+ \to \Gamma_-`
     - :math:`0`
     - :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`, both
       spaces bound (``ZeroMorphism`` over the face's two half-traces)
   * - :meth:`WallCondition.uniform(v)
       <orpheus.transport.boundary_conditions.WallCondition.uniform>` with
       :math:`v \neq 0`, or :meth:`~orpheus.transport.boundary_conditions.WallCondition.framed`
     - :math:`0` — **the same zero morphism**
     - :math:`\neq 0`
     - :class:`~orpheus.sn.boundary.realized_operators.ZeroMap`, the same
       type and the same two bound spaces
   * - :class:`~orpheus.geometry.deck.MirrorDeck`
     - :math:`P_\pi` — the ordinate permutation of the face's own mirror,
       unscaled (a symmetry has no amplitude)
     - :math:`0`
     - :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`,
       :math:`\alpha = 1`
   * - :meth:`WallCondition.specular(a)
       <orpheus.transport.boundary_conditions.WallCondition.specular>`
     - :math:`a\,P_\pi` — the same permutation at the wall's own amplitude
     - :math:`0`
     - :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`,
       :math:`\alpha = a`
   * - :meth:`WallCondition.white(d)
       <orpheus.transport.boundary_conditions.WallCondition.white>`
     - :math:`d\,(E \circ R_t)` — the factored Lambertian, emit
       :math:`\circ` contract
     - :math:`0`
     - :class:`~orpheus.sn.boundary.realized_operators.RankOneReemission`
   * - a wall with both chains
       (:meth:`~orpheus.transport.boundary_conditions.WallCondition.from_amplitudes`)
     - the sum of the two images, in the canonical chain order
     - :math:`0`
     - :class:`~orpheus.sn.boundary.realized_operators.RealizedSum`
   * - :class:`~orpheus.geometry.deck.TranslationDeck`
     - the identity relabelling between two distinct index sets,
       :math:`\Gamma_+(f') \to \Gamma_-(f)`, **earned** by the crossing
       check
     - :math:`0`
     - :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`,
       :math:`\alpha = 1`, partner :math:`f'`

**Vacuum and prescribed inflow differ ONLY in** :math:`q`. That is not
a tidiness observation, it is what makes a separate affine operator
**unnecessary** rather than merely untidy: there is no linear content
to distinguish them by, so a type that existed to carry prescribed
inflow's realization was carrying an *empty* linear factor and a
source that belongs to another tier. Measured on one
:class:`~orpheus.sn.problem.SNProblem`, one face, both conditions
realized through :meth:`SNProblem.realize_boundary_role
<orpheus.sn.problem.SNProblem.realize_boundary_role>` (the problem's
``bc[face]`` entry is that realization):

.. math::
   :label: bc-prescribed-zero-linear-factor

   \operatorname{realize}\bigl(\texttt{WallCondition.uniform}(q)\bigr)
   \;=\;
   \operatorname{realize}\bigl(\texttt{WallCondition.vacuum}()\bigr)
   \;=\;
   0_{\,\Gamma_+(f) \to \Gamma_-(f)}
   \qquad\text{for every } q .

.. (vv-status rationale) Structural/typing identity, not a solver claim: the
   two conditions realize to the SAME expression, so it is verified by object-
   and space-identity assertions, not by a value. Gated at the production tier
   by ``tests/gates/sn/operators/test_operator_block_role.py``
   (``test_prescribed_inflow_realizes_the_same_object_vacuum_does``:
   ``type(...) is type(...)`` plus ``domain is`` / ``codomain is`` against the
   TRACE's own ``outflow_space`` / ``inflow_space``) and by
   ``tests/gates/sn/operators/test_capability_survival.py``'s re-posed capability
   rows; the ``ZeroMorphism`` two-space contract itself is pinned by the
   foundation suite ``tests/gates/numerics/test_zero_operator_spaces.py``. An
   identity row against a SIBLING operator would prove only that the two AGREE,
   not that either names the right END; this row compares against the trace, so
   the hole is closed (gotcha 5 below, with the swap mutation).
.. vv-status: bc-prescribed-zero-linear-factor documented

.. admonition:: The fixture every ``[M]`` on this page's section uses
   :class: note

   A **2-group, heterogeneous, scattering-active slab**, chosen so that
   no result below rides on a degeneracy: two regions of mixtures
   ``A`` | ``D`` from the cross-section library at ``2g``
   (:math:`c \approx 0.90`–:math:`0.96`, with a non-symmetric
   :math:`\Sigma_{s0}` so the group transfer is transpose-active), 1.0
   and 2.0 cm thick, 6 + 6 cells, each face a deck or an integer tag
   bound to ``WallCondition.vacuum()`` unless a row says it declares an
   inflow, ``Quadrature.gauss_legendre(n_ordinates=8)`` — hence
   :math:`|\Gamma_+| = |\Gamma_-| = 4` per face and a 32-value trace
   (2 faces :math:`\times` 8 ordinates :math:`\times` 2 groups). The
   prescribed inflow is declared where a user declares it, as
   ``SNProblem(mesh, quad, materials, boundary_conditions={tag:
   WallCondition.uniform(2.5)})`` with the face's role that tag
   (``tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py``,
   ``_het_slab``, builds exactly this). The bulk source
   is a spatially uniform isotropic
   :math:`Q = 1` through
   :meth:`AngularSourceSink.from_isotropic
   <orpheus.transport.fields.AngularSourceSink.from_isotropic>`,
   identical on every leg of every comparison (gotcha 1 below explains
   why that matters), and ``inner_tol = 1e-13`` unless a row says
   otherwise. **Not** the campaign's ``placeholder_materials`` slab:
   :math:`\Sigma_s \equiv 0` there, so a delivered inflow never feeds a
   scattering source and the SI solve is a single sweep — the
   delivery-count claim survives that degeneracy but a flux-level claim
   on it would be configuration-blind.

``[M]`` on that fixture, ``xmin`` declaring
``WallCondition.uniform(2.5)`` (``scratch/wf25_bcB_measurements.py``,
section S1; the untracked script is re-run from outside the repository,
``PYTHONPATH=<worktree>``): the realized leaf is a ``ZeroMap`` and so is
the leaf of the vacuum face, ``type(prescribed) is type(vacuum)`` is
``True``, and each leaf is bound to the TRACE's own half-spaces:
``prescribed.domain is trace.outflow_space("xmin")`` and
``prescribed.codomain is trace.inflow_space("xmin")`` are both ``True``,
with ``domain = AngularFaceTraceSpace('angular_trace[xmin:outflow]',
shape=(4, 2))`` and ``codomain = …[xmin:inflow]…`` — i.e. exactly
:math:`\Gamma_+(f)` and :math:`\Gamma_-(f)`. There is no "default-source
spelling" to compare: ``WallCondition.uniform(0.0) ==
WallCondition.vacuum()`` is ``True`` (the drop rule: a zero weight is not
stored), the pair realizes to the same ``ZeroMap`` type, and its
``q.linf`` is ``0.0``, which is the honest reading of ":math:`L` does not
depend on :math:`q`".

One consequence is worth stating before the mechanics, because it is
the reason the collapse is *safe*: since the realized operator no
longer distinguishes the two conditions, **anything that must still
distinguish them has to read the DECLARATION** (the binding's sources),
not the operator, or be protected at admission. Two consumers are in
that position, and :ref:`bc-affine-channel-law-tier-consumers` says
how each is.


.. _bc-affine-channel-where-q-travels:

Where :math:`q` travels instead — the composite source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The typed :math:`q` is
:class:`~orpheus.transport.fields.AngularBoundarySourceSink` —
the *eager, whole-boundary, hub-bound* snapshot of the inflow,
packed into the unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
flat layout (one flat vector of ``layout.total_size`` values, with
per-face views). It carries angular-flux units, because :math:`q` is
added to :math:`\gamma_-\psi`, not to a volumetric rate.

It is not a second, parallel source object.
:class:`~orpheus.transport.full_field.FullField` **is**
``Composite[BulkField, BoundaryField]``, so the bulk and boundary
source are **one** object

.. math::
   :label: bc-composite-source

   q \;=\; q_{\rm bulk} \,\oplus\, q_\partial ,
   \qquad
   q_{\rm bulk} \in \text{BulkField},\quad
   q_\partial \in \text{BoundaryField},

.. (vv-status rationale) Representational: names the composite the solve's RHS
   already is (the ``TimedFullField`` pair). Its verifiable content is the
   class gate on every operator-output boundary plus the RHS construction
   rows in ``tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py``
   (``TestTheDeclarationIsNotInert``,
   ``TestTheSourceCannotBeSpecifiedTwice``). Not a solver claim.
.. vv-status: bc-composite-source documented

and the single construction point for it is
:func:`~orpheus.sn.solver._build_fixed_source_rhs`. **Both** inner
paths consume what that helper returns — the source-iteration driver
and the Krylov driver — so there is exactly one place where a boundary
source can be attached and exactly one place where a bug in the
attachment could live. A caller may hand
:func:`~orpheus.sn.solver.solve_sn_fixed_source` a bare
``(N, ng, *spatial)`` bulk array or the full composite
:class:`~orpheus.transport.timed_full_field.TimedFullField`; either
way the helper normalises to :eq:`bc-composite-source`.

The path from a declared condition to a packed snapshot is a **ladder**, each
rung delegating to the next so the packing rule is stated exactly once:

.. list-table:: The declaration → snapshot ladder on :class:`~orpheus.transport.fields.AngularBoundarySourceSink`
   :header-rows: 1
   :widths: 30 44 26

   * - rung
     - takes the source in the form …
     - delegates to
   * - :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.from_method`
     - **the form the problem holds** — a boundary condition someone
       declared. Reads ``binding.sources`` (the tagged faces whose
       condition has a non-zero source; vacuum, specular, white and every
       deck contribute nothing because their whole content is
       :math:`L`) and samples each face's
       :class:`~orpheus.transport.boundary_conditions.InflowTerms` on
       that face's own :math:`\Gamma_-(f)`
     - ``sample_inflow``, then ``_pack_inflow``
   * - :meth:`~orpheus.transport.boundary_conditions.RealizerTable.sample_inflow`
     - **the one source verb** — ``require`` the method table's admission
       of the source on the face's site, then the registered
       ``SourceRealizer.sample`` (SN: ``SNSource.sample``, which checks
       the ERR-047 placement and calls ``sample_terms``,
       ``q = uniform·1 + Σᵢ wᵢ samplerᵢ.sample(space)``), evaluated at
       :math:`\Gamma_-(f)` itself with exactly ``space.shape``
     - :func:`~orpheus.transport.boundary_conditions.realization.sample_terms`
   * - :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.prescribed_inflow`
     - **known per-face arrays** — the inflow ordinate slots of the
       named faces are written from them
     - ``_pack_inflow``
   * - ``_pack_inflow``
     - **the one packing rule** — writes the **inflow** ordinate slots of
       the named faces and leaves every other slot zero (a prescribed
       source's outflow rows are physically meaningless, so they are
       unrepresentable rather than overwritten)
     - ``zeros``
   * - :meth:`AngularBoundarySourceSink.zeros(mesh.angular_trace) <orpheus.numerics.field.Field.zeros>`
     - **nothing but the trace SPACE** — the all-zero source every
       sourceless condition needs, and the base the rungs above allocate
       before filling. The shared
       :meth:`Field.zeros <orpheus.numerics.field.Field.zeros>` body
       allocates ``np.zeros(space.shape)`` and nothing else
     - — (the bottom)
   * - *(sibling, not a rung)*
       :meth:`~orpheus.transport.fields._bases.AngularBoundaryField.from_face_arrays`
     - **every** face's **full** slot, outflow rows included — the
       general inherited constructor for non-inflow uses; nothing in the
       ladder delegates to it, and a prescribed inflow must not use it
     - —

Only the top rung starts from a **declaration**, and that is why it is
the top rung. A driver or a test that reaches for ``prescribed_inflow``
directly is supplying by hand what the declaration should have
produced, and is therefore exercising a path no user travels — see
:ref:`verification-user-path`, whose worked example this boundary
source is.

``[M]`` the ladder is not decorative:
``from_method`` (``SNProblem.declared_boundary_source()``) and
``prescribed_inflow(problem, {face: slot})`` return **bit-identical**
arrays for a declared ``WallCondition.uniform(2.5)``
(``np.array_equal`` on the flat values, section S8 of
``scratch/wf25_bcB_measurements.py``), and ``from_method`` returns the
zero trace field (``q.linf == 0.0``) for ``WallCondition.vacuum()``, a
``MirrorDeck``, ``WallCondition.white(1.0)``,
``WallCondition.specular(0.5)`` and ``WallCondition.uniform(0.0)``
(section S1). Only ``WallCondition.uniform(2.5)`` gives
``q.linf == 2.5``. On the documented fixture the packed :math:`q` has 8
non-zeros in 32 values: :math:`|\Gamma_-(\texttt{xmin})| = 4` inflow
ordinates :math:`\times` 2 groups, and nothing anywhere else on the
trace.

The claim the channel exists to support is the affine law **evaluated
on the answer**. For prescribed inflow :math:`L = 0`, so at
convergence

.. math::
   :label: bc-single-delivery

   \gamma_-\psi\big|_f \;=\; q_f
   \qquad\text{exactly, and once.}


.. implements:: bc-single-delivery
   :by: orpheus.transport.boundary_conditions.realization.sample_terms

   **Implemented by** 8 sites. Every symbol that executes this
   equation's arithmetic is declared, not only the canonical one: a
   test is adjudicated against the transcription it actually ran, so
   declaring a single site would refute the tests that exercise the
   others.

.. implements:: bc-single-delivery
   :by: orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face

.. implements:: bc-single-delivery
   :by: orpheus.sn.solver._build_fixed_source_rhs

.. implements:: bc-single-delivery
   :by: orpheus.sn.boundary.realizer.SNSource.sample

.. implements:: bc-single-delivery
   :by: orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink

.. implements:: bc-single-delivery
   :by: orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink.from_method

.. implements:: bc-single-delivery
   :by: orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink._pack_inflow

.. implements:: bc-single-delivery
   :by: orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink.prescribed_inflow

Three failure directions are distinguishable in one reading of that
equation, with no reference solver and no discretization dependence:
a doubled delivery reads :math:`2q`, a lost channel reads :math:`0`,
and a wrong linear factor reads :math:`q + L\gamma_+\psi \neq q`.

.. warning::

   ⚠ **That reading is bit-exact on source iteration; on Krylov it is
   not guaranteed exact, and the gate bounds it.** The difference is a
   property of the inner solver, not of the delivery. It is the one
   measurement on this page most likely to be gated wrongly, in either
   direction.

   The sweep *writes* the inflow seed into :math:`\gamma_-\psi`, so on
   SI the converged trace is :math:`q` **as a copy**: ``[M]``
   ``0x1.4000000000000p+1``, :math:`\lVert \gamma_-\psi - q
   \rVert_\infty = 0.0` exactly, on all five fixtures of section S2 below (SI reads
   ``0.0`` on every one, at ``inner_tol`` ``1e-13`` and ``1e-10``).
   GMRES instead *solves* for the trace rows along with everything else,
   so the converged trace can carry the iteration residual:

   .. list-table:: :math:`\lVert \gamma_-\psi - q \rVert_\infty` at :math:`q = 2.5`, Krylov, at ``471c4737``
      :header-rows: 1
      :widths: 34 12 18 10 12 14

      * - fixture
        - ``inner_tol``
        - :math:`\lVert \gamma_-\psi - q \rVert_\infty`
        - ULP
        - exact?
        - :math:`/\,(q\cdot\texttt{tol})`
      * - GL-8, ``vac|vac``
        - ``1e-13``
        - ``0.0``
        - 0
        - **yes**
        - ``0.0``
      * - GL-16, ``vac|vac``
        - ``1e-13``
        - ``4.441e-16``
        - 1
        - no
        - ``1.78e-3``
      * - GL-8, ``vac|refl``
        - ``1e-13``
        - ``8.882e-16``
        - 2
        - no
        - ``3.55e-3``
      * - GL-8, 12+12 cells
        - ``1e-13``
        - ``8.882e-16``
        - 2
        - no
        - ``3.55e-3``
      * - GL-8, ``vac|vac``
        - ``1e-12``
        - ``2.220e-15``
        - 5
        - no
        - ``8.88e-4``
      * - GL-8, ``vac|vac``
        - ``1e-11``
        - ``4.619e-14``
        - 104
        - no
        - ``1.85e-3``
      * - GL-8, ``vac|vac``
        - ``1e-10``
        - ``1.732e-14``
        - 39
        - no
        - ``6.93e-5``
      * - GL-8, ``vac|vac``
        - ``1e-8``
        - ``8.744e-11``
        - 196 894
        - no
        - ``3.50e-3``

   ``[M]`` ``scratch/wf25_bcB_measurements.py`` sections S2 (the five
   fixtures, both solvers) and the tolerance sweep printed under it: the
   fixture is the 2-group heterogeneous slab above, ``xmin`` declaring
   ``WallCondition.uniform(2.5)``, run from outside the repository with
   ``PYTHONPATH`` the worktree. Read the column of ULP counts down the
   tolerance sweep: it is ``0, 5, 104, 39, 196 894`` as ``inner_tol``
   goes from ``1e-13`` to ``1e-8``, not monotone, so the ULP count is
   not a stable quantity of this fixture, while the ratio to
   :math:`q\cdot\texttt{inner\_tol}` stayed below ``3.6e-3`` in all eight
   Krylov rows. Seven of the eight Krylov rows are inexact. The first
   measurement of the same five fixtures (P3, 2026-08) read ``18``, ``7``,
   ``1``, ``23`` and ``27 580`` ULP; the tree at ``471c4737`` reads ``0``,
   ``1``, ``2``, ``2`` and ``39``. The change that moved the readings
   is in part the sweep-left preconditioner of #200 (2026-10-04,
   ``orpheus/sn/solver.py``, ``Krylov(preconditioner=(L+C)^{-1})``):
   ``[M]`` ``scratch/wf25_bcB_precond_probe.py``, run from outside the
   repository under ``python -O`` with ``PYTHONPATH`` the worktree on the
   keystone fixture (GL-8, ``xmin`` declaring ``WallCondition.uniform(2.5)``,
   ``inner_tol`` ``1e-13``), reads ``0`` ULP with the sweep preconditioner
   and ``11`` ULP (``dev = 4.88e-15``) with ``Krylov.__init__`` forced to
   ``IdentityOperator()`` (2026-10-24). The remaining gap, 11 against the
   ``18`` ULP of P3, is unattributed. This is one substitution on one
   fixture, not a bisect. It does not matter to the gate below, which is
   built not to depend on the reading being exact.

   **What the gate ships, and why it is TWO legs.** A single
   ``assert_array_equal`` parameterized over both inner solvers asserts an
   exactness GMRES does not promise: it is green on the keystone's own
   fixture today (GL-8, ``vac|vac``, ``0`` ULP) and would be red on seven
   of the eight Krylov rows above, for a reason that has nothing to do
   with the delivery count. The keystone
   ``test_the_declared_boundary_law_holds_on_the_answer`` therefore
   splits the claim:

   * **the delivery COUNT**, ``assert_allclose(rtol=0, atol=1e-9)`` on
     both legs — eleven orders of margin over readings of ``0.0`` /
     ``2.5`` / ``5.0``, so the three cases stay exactly
     distinguishable while the failure message stays readable; and
   * **the exactness**, per path: ``assert_array_equal`` on SI, and
     ``assert_array_almost_equal_nulp(nulp=64)`` on Krylov, a budget
     that admits every ``1e-13`` reading above (at most 2 ULP).

   Splitting is what keeps the two claims from eroding each other. Do
   **not** relax the SI leg to a tolerance to "match" Krylov: that
   discards the one exact reading in the whole claim. And note what the
   table above implies about the ULP budget — a **fixed** ``nulp`` is a
   floating-point claim about a quantity that is really an *iteration
   residual*, so it is stable only at the fixture's ``inner_tol``.
   ``nulp=64`` has headroom at ``1e-13`` (at most 2 ULP) and is exceeded
   at ``1e-11`` (104 ULP) and ``1e-8``. That is a deliberate trade: a ULP
   budget cannot be quietly loosened the way an ``rtol`` can, so if
   someone changes the fixture's tolerance the row may go red and asks
   to be re-derived, which is the intended conversation.

The pre-P3 double delivery is likewise exact on SI —
``5.000000000000``, i.e. :math:`2.5 + 2.5` in IEEE — and on Krylov it
did not produce a reading at all, because the solve raised (see
:ref:`bc-affine-channel-the-defect`). See
:ref:`bc-affine-channel-two-channels` for the *flux*-level comparison
and which assertion belongs there.

.. note::

   :eq:`bc-single-delivery` carries **no** ``vv-status`` sentinel
   because it needs none: it is a genuine L1 equation claim with a
   committed gate. ``tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py``'s
   ``test_the_declared_boundary_law_holds_on_the_answer`` carries
   ``@pytest.mark.verifies("bc-single-delivery")`` on both inner-solver
   parameters, plus ``@pytest.mark.catches("ERR-075")`` — the error
   catalog's entry for the double delivery described in
   :ref:`bc-affine-channel-the-defect`. Its negative control
   (``test_an_undeclared_face_has_a_zero_inflow_trace``) pins the
   complement: an all-vacuum mesh under the identical bulk source must
   read :math:`\gamma_-\psi = 0` on **both** faces, without which the
   keystone would pass against an implementation that wrote :math:`q`
   onto every inflow trace unconditionally.


.. _bc-affine-channel-the-defect:

The defect the collapse closed — an affine map in a linear slot
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Until P3 the SN realizer answered ``PrescribedInflow`` with an
``IncomingSourceOperator`` whose ``apply`` **ignored its input** and
returned :math:`q`. That is an affine map declared as a
:class:`~orpheus.numerics.operator.LinearOperator`. It carried no
:attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY` stamp, and the
stamp's absence was *believed* to fence it out of the :math:`B` block.

**It did not.** The property that assembled :math:`B` was
``SNBoundaryOperator._face_laws``
(:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`), and it
read

.. code-block:: text

   return {
       face: self.problem.bc[face]
       for face in self.problem.angular_trace.layout.faces
   }

— every face in the trace layout, with **no** ``block_role`` filter.
:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` carried the
stamp itself, applied every law it collected, and never asked a leaf
what its role was. So the affine leaf reached the block, and :math:`B`
stopped being linear. Its successor, ``SNBoundaryOperator._face_operators``,
reads ``self.problem.bc[face]`` over the same trace layout and still
applies no ``block_role`` filter (gotcha 3 below).

``[M]`` in a worktree pinned at the pre-carve commit (``ef4c3537``,
which already carried the P2′ source wiring), same 2-group
heterogeneous slab, ``xmin`` declaring
``PrescribedInflow(ConstantInflowSource(2.5))`` (today
``WallCondition.uniform(2.5)``):

.. list-table:: The operator tier, before and after the collapse
   :header-rows: 1
   :widths: 42 29 29

   * - measured on :math:`B`
     - pre-P3
     - post-P3
   * - faces collected (``_face_laws``, today ``_face_operators``)
     - ``['xmax', 'xmin']``
     - ``['xmax', 'xmin']`` — **unchanged; there was never a filter**
   * - ``xmin`` leaf's ``block_role``
     - ``None``
     - ``BlockRole.BOUNDARY``
   * - realized ``xmin`` leaf
     - ``IncomingSourceOperator``
     - ``ZeroOperator`` at P3; ``ZeroMap`` at wave 2
   * - :math:`\lVert B(0) \rVert_\infty`, bulk block
     - ``0.0``
     - ``0.0``
   * - :math:`\lVert B(0) \rVert_\infty`, **boundary block**
     - ``2.5`` — a linear operator has :math:`B(0) = 0`
     - ``0.0``
   * - :math:`\lVert B(2x) - 2B(x) \rVert_\infty`
     - ``2.5``
     - ``0.0``
   * - ``B.is_adjointable``
     - ``False`` — the affine leaf declined a transpose, and the
       per-face intersection rule made the **whole** block
       non-adjointable
     - ``True``

``[M]`` the post column re-read at ``471c4737``
(``scratch/wf25_bcB_measurements.py`` section S4, the same fixture with the
declaration bound to a tag): the leaf is a ``ZeroMap``, ``block_role`` is
``BlockRole.BOUNDARY``, :math:`\lVert B(0) \rVert_\infty = 0.0` and
:math:`\lVert B(2x) - 2B(x) \rVert_\infty = 0.0` on both faces, and
``B.is_adjointable`` is ``True``.

Two end-to-end consequences followed, and they are different on the
two inner solvers.

**On source iteration the inflow was delivered TWICE.** :math:`q_\partial`
enters the SI right-hand side and :math:`B\psi` enters as a coupling
gain, so once P2′ wired the declared law into the source channel while
the affine operator was still in the block, both carried it. ``[M]``
on the documented fixture, comparing a declared inflow (``D``) against
the same inflow hand-supplied on an all-vacuum mesh (``C``) against a
no-inflow control (``V``):

.. list-table:: Converged inflow trace :math:`\gamma_-\psi` at ``xmin``, and the delivery ratio
   :header-rows: 1
   :widths: 34 22 22 22

   * - configuration
     - pre-P3, SI
     - pre-P3, Krylov
     - post-P3, both
   * - ``D`` — declared ``PrescribedInflow(2.5)`` (today
       ``WallCondition.uniform(2.5)``)
     - ``5.0`` (bit-exact: ``0x1.4000000000000p+2``)
     - raises (below)
     - ``2.5``
   * - ``D0`` — declared, source channel disabled (the pre-P2′
       behaviour, reached at the time by monkeypatching
       ``from_problem_laws`` back to the then-current mesh-keyed
       ``zeros_on``; today's spelling of the same probe is
       ``zeros(mesh.angular_trace)``)
     - ``2.5``
     - raises (below)
     - —
   * - ``C`` — all-vacuum mesh, :math:`q_\partial` supplied directly
     - ``2.5``
     - ``2.5``
     - ``2.5``
   * - ``V`` — all-vacuum, no inflow (control)
     - ``0.0``
     - ``0.0``
     - ``0.0``
   * - :math:`\lVert \varphi_D - \varphi_V \rVert_\infty /
       \lVert \varphi_C - \varphi_V \rVert_\infty`
     - ``2.0000000000``
     - —
     - ``1.0000000000``

Every entry in that table is a trace value read at the SI convergence
point, where it is bit-exact; the ``post-P3, both`` column reads
``2.5`` on SI bit-exactly and on Krylov to the iteration residual (the
warning at :eq:`bc-single-delivery` carries that scan). The ratio is
**exactly** two pre-P3 and **exactly** one post-P3, on
every fixture measured. Note what row ``D0`` says about the history:
pre-P2′ a declared inflow was *not* inert — it was delivered, through
the affine operator inside :math:`B`. P2′ then added the channel it
should have travelled all along, and the two deliveries summed. P3
removes the wrong one.

**On the Krylov path it was worse, and had been for longer.** GMRES's
residual bookkeeping rests on the Arnoldi relation
:math:`A V_k = V_{k+1} H_k`, which requires :math:`A` to be linear. An
affine :math:`A(x) = A_{\rm lin}(x) + c` breaks it, so the residual
SciPy tracks internally becomes meaningless: it reports convergence
while the iterate does not solve the equation.
:func:`~orpheus.sn.solver._check_convergence_claim` catches exactly
that — it recomputes the honest equation residual
:math:`\lVert A\psi - q \rVert / \lVert q \rVert` after any convergence
claim and raises
:class:`~orpheus.sn.solver.ConvergenceClaimError`:

.. code-block:: text

   solve_sn_fixed_source[krylov]: the within-group solve claimed
   convergence (running residual 0.000e+00 < tol 1.0e-13) but the
   honest equation residual is ‖Aψ − q‖/‖q‖ = 2.313e+00 — the
   iteration's fixed point does not solve the equation (the #282
   lag-death class: a stale/lagged block inside M; the free-identity
   stop is structurally blind to it).

``[M]`` pre-P3 **both** rows raise: ``D`` at a defect of ``2.313e+00``
and ``D0`` at ``2.504e+00`` on the documented fixture (P3's
verification plan measured ``1.718e+00`` on its own probe — the
magnitude is fixture-dependent, the :math:`\mathcal{O}(1)` scale is
not). So **a declared prescribed inflow combined with the Krylov inner
solver was simply UNUSABLE**, at the pre-carve commit *and* before the
source channel existed. P3 is a bug **fix**, not prophylaxis; ``[M]``
post-P3 the same configuration converges on both inner solvers and the
convergence-claim check is silent.

.. warning::

   **Only one fence was ever real, and it is the reason a green suite
   never saw any of this.** ``prescribed_inflow`` was not a registered
   ``BC`` kind: the ``BC`` tag's ``params`` was ``dict[str, float]``, so a
   tag could carry ``{"albedo": 0.7}`` but never an inflow distribution,
   and the SN tag registry admitted only ``{vacuum, reflective}``. A
   non-trivial prescribed inflow was reachable **only** by constructing
   the law object and installing it, which no production driver did.
   Issue **#189** was opened for registering the kind in the constant
   case. Wave 2 removed the fence by removing the registry: a face's
   role is a deck or an integer tag, the tag is bound to a
   ``WallCondition`` at posing, and ``boundary_conditions={tag:
   WallCondition.uniform(2.5)}`` **is** the declaration, so the path a
   user travels now reaches the source channel and the solve (the gates
   of ``tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py``).

   The two defects protected each other. The bypass — every existing
   non-vacuum test supplied :math:`q_\partial` by hand — was why no test
   could observe the law tier at all; and the law tier's inertness (as
   far as any driver was concerned) was why the bypass looked harmless.
   That is the :ref:`verification-user-path` doctrine's canonical
   instance: the gates were operator-tier or right-hand-side-tier, and
   never a *solve* started from a *declaration*.


.. _bc-affine-channel-two-channels:

The two channels, and what "equivalent" means
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two user-visible ways to put flux into :math:`\Gamma_-`:

1. **Declare the condition** — bind a tag to
   ``WallCondition.uniform(v)`` (or ``WallCondition.framed(profile)``) in the
   problem's ``boundary_conditions`` and let
   :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.from_method`
   find it (``SNProblem.declared_boundary_source``). This is the user path.
2. **Supply** :math:`q_\partial` **directly** — build the source with
   :meth:`~orpheus.transport.fields.AngularBoundarySourceSink.prescribed_inflow`
   and pass the composite. This is what the existing non-vacuum
   manufactured-solution coverage does, deliberately bypassing the
   declaration tier.

Post-P3 the two reach the same fixed point, and the *strength* of that
agreement is the sharpest single statement about what the carve did:

.. list-table:: Declared (``D``) versus hand-supplied (``C``), post-P3
   :header-rows: 1
   :widths: 34 17 15 17 17

   * - fixture
     - inner
     - :math:`\lVert \varphi_D - \varphi_C \rVert_\infty`
     - ``array_equal``
     - delivery ratio
   * - slab 2G het, ``xmin`` only
     - SI
     - ``0.000e+00``
     - ``True``
     - ``1.0000000``
   * - slab 2G het, ``xmin`` only
     - Krylov
     - ``0.000e+00``
     - ``True``
     - ``1.0000000``
   * - slab 2G het, **both** faces declaring
     - SI / Krylov
     - ``0.000e+00``
     - ``True``
     - ``1.0000000``
   * - slab 2G het, ``xmax`` reflective
     - SI / Krylov
     - ``0.000e+00``
     - ``True``
     - ``1.0000000``
   * - slab 2G het, ``gauss_legendre(16)``
     - SI / Krylov
     - ``0.000e+00``
     - ``True``
     - ``1.0000000``

``[M]`` re-measured at ``471c4737`` (section S3 of
``scratch/wf25_bcB_measurements.py``, the declared problem against a
vacuum problem handed the same :math:`q_\partial` by ``prescribed_inflow``):
``array_equal`` is ``True`` and the maximum difference ``0.000e+00`` in
all eight (fixture, inner solver) cells above, at :math:`\lVert \varphi
\rVert_\infty = 6.838442` for the first fixture, and the delivery ratio
:math:`\lVert \varphi_D - \varphi_V \rVert_\infty /
\lVert \varphi_C - \varphi_V \rVert_\infty` is ``1.0000000000`` in each.

**Bit-identical, on every fixture and both inner solvers — and that is
structural, not luck.** Post-P3 the two configurations are the *same
floating-point program*: the realized :math:`B` is a ``ZeroMap``
in both (by :eq:`bc-prescribed-zero-linear-factor`), every other
operator is built from the same ``(mesh, quadrature, materials)``, and
the two right-hand sides carry bit-equal :math:`q_\partial` arrays. The
only precondition is that last one — the declared spec must *evaluate*
to the same floats the hand-built slot was filled with, which holds by
construction for a constant source compared against a constant slot,
and is exactly the property
:func:`~orpheus.transport.boundary_conditions.realization.sample_terms`
guarantees (it evaluates the declared terms at :math:`\Gamma_-(f)`
itself and refuses any sample that is not exactly ``space.shape``, so
the declared and the hand-built slot are the same array shape).

.. warning::

   ⚠ **The pre-P3 number for this same comparison is 1.998e-13 — and
   a tolerance-based two-channel flux gate is therefore BLIND to the
   defect P3 removed.** ``[M]`` at ``ef4c3537``, rows ``D0`` versus
   ``C``: :math:`\lVert \varphi_{D0} - \varphi_C \rVert_\infty =
   1.9984\times10^{-13}` with ``array_equal = False``, against
   :math:`\lVert \varphi \rVert_\infty = 6.838442` — a relative
   :math:`2.92\times10^{-14}`, i.e. :math:`0.3 \times
   \texttt{inner\_tol}` at ``inner_tol = 1e-13``.

   That non-zero difference was never an error. It was the
   **signature** of the two channels being two different computations
   reaching the same fixed point: the delivery travelled through
   :math:`B` as an off-diagonal gain rather than through
   :math:`q_{\rm ext}`, so the reduction tree differed. Post-P3 it
   collapses to ``0.0`` because there is one computation.

   The consequence for gate design is concrete. A two-channel flux
   comparison at ``rtol = 10 × inner_tol`` — a defensible tolerance
   derived from the iteration floor — **passes** on a
   :math:`2.9\times10^{-14}` relative difference. It can see the
   :math:`2\times` double delivery and it cannot see a re-introduced
   *structural* difference between the channels. Assert
   ``np.array_equal`` on the flux for the constant-source case, and
   state the precondition; if a future non-constant spec makes the two
   :math:`q_\partial` arrays differ in their last bit, pin **that**
   (a source-tier ``array_equal`` on the two :math:`q_\partial`
   arrays) and let the flux comparison relax — do not silently swap
   an exact gate for a tolerance and call it the same claim.

**Which gate asserts what.** The distinction matters in both
directions and is easy to get backwards:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - quantity
     - assertion
     - why
   * - converged inflow trace
       :math:`\gamma_-\psi\big|_f` against the declared :math:`q_f`
       (:eq:`bc-single-delivery`) — the **delivery count**
     - ``assert_allclose(rtol=0, atol=1e-9)``, both inner solvers
     - eleven orders of margin over ``0.0`` / ``2.5`` / ``5.0``, so the
       three readings stay exactly distinguishable. This is the leg
       that IS the claim; the exactness legs below are about *how* the
       trace is reached, not *how much* arrived
   * - the same trace — the **exactness**, SI leg
     - ``assert_array_equal`` — bit-exact
     - the sweep *writes* the seed, so :math:`q` arrives as a copy;
       measured ``0.0`` deviation on every fixture and every tolerance.
       If a future moment-resolved trace makes it inexact, descend to
       ``assert_array_almost_equal_nulp`` and say why — **not** to
       ``rtol``
   * - the same trace — the **exactness**, Krylov leg
     - ``assert_array_almost_equal_nulp(nulp=64)``
     - GMRES *solves* for the trace rows, so the reading carries the
       iteration residual (18 ULP measured, and tolerance-dependent —
       see the warning above). A fixed ULP budget is chosen over an
       ``rtol`` precisely because it cannot be quietly loosened
   * - converged **flux**, declared versus hand-supplied
     - ``np.array_equal`` post-P3 (measured on both inner solvers; the
       gate asserts the SI leg), with the bit-equal
       :math:`q_\partial` precondition stated
     - the two solves are the same float program; a tolerance here is
       blind to the ``1.998e-13`` structural-difference signature
       above
   * - the operator tier: :math:`B(0) = 0`,
       :math:`B(2x) = 2B(x)`, :math:`B(x{+}y) = B(x){+}B(y)` over the
       whole shipped declaration set
     - exact, and reachable **without a solve**
     - it catches an affine term anywhere in the block at
       *construction* time — the day the next affine operator is
       written, before it reaches a solver. It is blind to the
       delivery count (a source doubled in :math:`q_\partial` leaves
       :math:`B` perfectly linear) and blind to a wrong-but-linear
       :math:`L`

.. note::

   A **superposition** gate — "the flux is affine in the inflow
   amplitude with slope :math:`s`" — was considered for this claim and
   **rejected**. The double delivery *is* :math:`q \mapsto 2q`, which
   is still exactly affine in :math:`q`: :math:`\varphi(a) =
   \varphi(0) + a\,s` holds for every :math:`s`, including
   :math:`s = 2 s_{\rm true}`. The measured functional's invariance
   group contains the scale factor, so the gate is designed-green on
   the entire defect class at every tolerance and every refinement —
   a textbook ``vv-principles`` **Mode 12** non-catcher. Only the
   *coefficient*, pinned against an independent reference, has teeth,
   and that is what :eq:`bc-single-delivery` does directly.

   A frozen snapshot of the converged flux was rejected for a
   different reason: it asserts that production equals a recording of
   production, cannot distinguish "delivered twice" from "the spec
   changed", and would have to be re-baselined by the very carve it is
   meant to gate.


.. _bc-affine-channel-gotchas:

Five gotchas, all measured
~~~~~~~~~~~~~~~~~~~~~~~~~~

**1. The two bulk-source arms take different types, and mixing them
reads as a boundary difference.**
:func:`~orpheus.sn.solver._build_fixed_source_rhs` accepts either a
per-ordinate ``(N, ng, *spatial)`` array *or* a composite whose bulk
leaf is an
:class:`~orpheus.transport.fields.AngularSourceSink`. The two are
**not** interchangeable spellings of one source. The array form is
already a per-ordinate density :math:`Q_n`, while
:meth:`AngularSourceSink.from_isotropic
<orpheus.transport.fields.AngularSourceSink.from_isotropic>`
takes a *scalar* :math:`Q(\vec r, g)` and applies the
:math:`1/\sum_n w_n` projection at the producer boundary (Pattern 7).
So the same-looking pair of spellings differs by exactly the
quadrature's total weight. ``[M]`` on the fixture above,
``np.ones((N, ng, …))`` against
``AngularSourceSink.from_isotropic(np.ones((ng, …)))``:
:math:`\varphi[0,0] = 4.524718` versus :math:`2.262359`, ratio
**2.000000** exactly — which is :math:`\sum_n w_n = 2` for a
Gauss–Legendre rule on :math:`\mu \in [-1, 1]` (it would be
:math:`4\pi` for a full-sphere rule).

A factor of :math:`\sum_n w_n` in the **bulk** presents as a
disagreement of the *boundary* channel, because that is the only thing
the comparison was varying. Any two-channel comparison must therefore
supply the *identical* bulk spelling on both legs; building both legs'
source through one helper makes the trap unspellable rather than
avoided by care.

**2. Double specification is REFUSED, not resolved.** If the problem
declares a non-zero source (a tag bound to ``WallCondition.uniform(v)``
with ``v != 0``) **and** the composite
``external_source`` carries a non-zero boundary leaf, there are two
answers to one question.
:func:`~orpheus.sn.solver._build_fixed_source_rhs` raises. Neither
alternative is acceptable: adding double-counts (the exact defect this
section documents), and overriding makes one of the two inputs a silent
no-op. **Both wrong answers are quiet**, which is why the refusal is
loud. A composite whose boundary leaf is all-zero alongside a
declaration is fine — that is the normal declared-condition path.
``[M]`` section S8 of ``scratch/wf25_bcB_measurements.py``: the double
specification raises ``ValueError`` ("specified TWICE"), the row
``TestTheSourceCannotBeSpecifiedTwice`` of
``tests/gates/sn/solve/test_declared_inflow_reaches_the_rhs.py``.

**3. The** ``BlockRole.BOUNDARY`` **stamp is honest metadata; it is
NOT a fence.** ``[M]`` and stated plainly because a reader will
otherwise assume it is load-bearing: ``_face_operators`` (``_face_laws``
before it) never filtered on ``block_role``, so the stamp's presence or
absence changes **nothing** about :math:`B`'s behaviour. Making
:func:`~orpheus.transport.boundary_conditions.realization.stamp_boundary_role`
the identity in ``orpheus.sn.boundary.realized_operators`` reddens 12 of
the 51 rows of
``test_operator_block_role.py``, ``test_declared_inflow_reaches_the_rhs.py``
and ``test_declared_law_is_linear.py`` (12 in the first, 0 in each of the
other two): every red is a direct stamp assertion and none is numerical
(``MUT=stamp PYTHONPATH=scratch:<worktree> python -O -m pytest -p
wf25_bcB_mutplug <the three files>``; the control, unmutated, passes 51
of 51). The stamp claim needs its own assertion, and a value gate
structurally cannot supply one. What the stamp *does* do is describe a
leaf's role for a reader and for the block-role marker classes; what it
never did is keep anything out of the block. Post-P3 the stamp is also
no longer an exception: every realized face type carries it, the zero
map of a prescribed inflow included, because the leaf is genuinely
linear.

**4. A declared prescribed inflow is REFUSED on a radial outer face
(the carrying 1-D curvilinear mesh), at posing — and the hand-supplied
channel is not.** This is the one place where the two channels of
:ref:`bc-affine-channel-two-channels` are *not* interchangeable, and it
is a deliberate loud deferral (issue #653) rather than a defect. ``[M]``
section S4 of ``scratch/wf25_bcB_measurements.py``, a 2-group
heterogeneous sphere with the outer tag bound to
``WallCondition.uniform(2.5)``: ``SNProblem(...)`` raises
``BoundaryRefused`` carrying ``issue == 653``, before any sweep:

.. code-block:: text

   BoundaryRefused: sn realizes ['DiffuseChain', 'InflowTerms',
       'SpecularChain', 'ZeroNetwork']; not WallCondition (PLUMBING: the
       radial characteristic's corner action is ruled only for a
       permuting return or a zero return)

while the hand-supplied channel delivers: section S5, the same sphere
with a vacuum outer tag and ``q_∂`` handed to ``prescribed_inflow``,
reads :math:`\gamma_-\psi(\texttt{xmax}) = [2.5, 2.5, 2.5, 2.5]` on both
inner solvers.

A seed-carrying mesh closes its ray boundary at the off-quadrature
:math:`\mu = \pm 1` corner, and
:class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`'s
corner block is a **linear** operator on the ray alone. A prescribed
inflow's corner value is a free parameter, not a function of the
outflow ray, so a linear corner block structurally cannot carry it,
whereas the *source* side has a ruled three-arm inflow-corner law and
delivers it fine. The refusal is the realizer table's, keyed on the
site class: ``SNSource.admit`` (and ``SNDiffuseChain.admit``) returns
the corner refusal at ``SiteClass.RAD_OUTER`` (until wave 2 such a
declaration constructed, and refused with ``NotImplementedError`` at the
first sweep). What reaches the corner block is therefore never a
sourced face, and its predicate, ``_has_ruled_corner_action``, asks only
the realized TYPE (``permutes_ordinates or is_vacuum``). The family test
the retired code made is not needed because a zero source is not a
member of a family waiting for a value: ``WallCondition.uniform(0.0)``
**is** the vacuum (``[M]`` section S4: ``uniform(0.0)`` and ``vacuum``
both pose on the sphere and realize a ``ZeroMap``; ``uniform(2.5)``
and ``white`` are refused with issue 653).

**5. An identity row against a sibling proves the two conditions AGREE
— not that either names the right END; the row that does is the one
against the trace.** :eq:`bc-prescribed-zero-linear-factor` is gated
by ``test_prescribed_inflow_realizes_the_same_object_vacuum_does``. A
row that only asserted ``prescribed.domain is vacuum.domain`` would stay
**green** if **both** leaves were bound swapped
(:math:`\text{domain} := \Gamma_-`, :math:`\text{codomain} := \Gamma_+`),
because it compares the two operators against each other and never
against the trace; and a size check cannot catch the swap either:
:math:`|\Gamma_+| = |\Gamma_-|` on each of the 46 (rule, face) pairs
measured (9 rules from five families; ``[M]`` section S7 of ``scratch/wf25_bcB_measurements.py``:
46 of 46 (quadrature, face) pairs that build a signed-projection table,
9 rules :math:`\times` 6 face names = 54 attempted and the 8 refused
being the two one-dimensional rules on the four ``y``/``z`` faces;
control: a hand-built lopsided table reads unequal). This is a textbook
``vv-principles`` **Mode 12** blindness: the measured functional (mutual
agreement) has the error class in its invariance group.

The production row is not that row. It asserts ``op.domain is
sn.angular_trace.outflow_space(face)`` and ``op.codomain is
sn.angular_trace.inflow_space(face)`` for the vacuum face and the
prescribed face, against the **trace's** own spaces and not against a
sibling, and the ``ZeroMap`` derives its two ends from the trace ladder
at one site (``_derive_inner``). ``[M]`` the swap mutation reddens it:
``MUT=swap PYTHONPATH=scratch:<worktree> python -O -m pytest -p
wf25_bcB_mutplug tests/gates/sn/operators/test_operator_block_role.py::TestBoundaryLeaves::test_prescribed_inflow_realizes_the_same_object_vacuum_does``
passes unmutated (1 passed) and fails with
``ZeroMap._derive_inner`` rebound to a ``ZeroMorphism`` of swapped ends
(1 failed). The same discrimination for the hand-built type
(``tests/gates/numerics/test_zero_operator_spaces.py``,
``test_the_ends_are_not_swapped``, on an intentionally *unequal*
:math:`3 \times 7` fixture, which is what makes a swap observable there)
and for the specular chain
(``tests/gates/sn/operators/test_specular_deck_chain.py``) stands beside
it.


.. _bc-affine-channel-law-tier-consumers:

The two consumers that must still tell a source apart
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Until wave 2, collapsing the operator raised an obvious question: what
breaks in the consumers that needed to tell prescribed inflow apart?
Nothing did, because both read the LAW family (``isinstance(law,
PrescribedInflow)``), not the realized operator, and both said so in
their own comments before P3 was planned. The law objects are gone
(:ref:`bc-wave2-crosswalk` maps each retired name to its successor, and
:ref:`bc-architecture-of-record` states the present architecture), and
so is that discriminator, but the question stays: the realized operator
of a sourced face is the zero map of a vacuum face, so a consumer of the
realized operator cannot see :math:`q`. Today each of the two is
protected differently, and ``[M]`` both hold
(``scratch/wf25_bcB_measurements.py``, sections S4 and S6):

.. list-table::
   :header-rows: 1
   :widths: 30 40 30

   * - consumer
     - what it asks
     - now
   * - :class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`'s
       corner block, via ``_has_ruled_corner_action``
     - of the realized TYPE: does it permute ordinates, or is it the
       zero map? The corner block is linear on the off-quadrature ray,
       so a free-parameter inflow has no expressible action there
     - **never reached by a source**: SN's table refuses every source
       and the diffuse chain at a radial outer face at admission
       (``SNSource.admit``, issue #653), at posing; ``uniform(2.5)``
       and ``white(1.0)`` are refused there, ``uniform(0.0)`` and
       ``vacuum()`` are admitted (gotcha 4 above)
   * - :class:`~orpheus.sn.acceleration.dsa.DSALowOrderSystem`
     - is the face's low-order edge row proven? A zero response (a
       ``ZeroMap``, ``is_vacuum``) gives the Marshak row and an
       ordinate-permuting face that returns everything
       (``closes_axis``) gives the zero-net-current row (a partial
       specular reflector permutes the ordinates too, but its net current
       is not zero, so it is refused: ERR-094); and a face whose
       condition declares a SOURCE has **no row** in the low-order edge
       system
     - **the declaration is read**: the face is excluded when its name is
       in ``problem.binding.sources``, because the realized operator is
       the linear part only. ``[M]`` DSA admits ``vacuum()``,
       ``MirrorDeck`` and ``uniform(0.0)`` and refuses ``uniform(2.5)``
       (source in the binding), ``specular(0.5)`` and ``white(1.0)``

What changed with the carve is the discriminator, not the discipline.
The family test existed because a ``PrescribedInflow`` at its default
zero source was a different object from the vacuum, and a value test on
:math:`q` would have admitted it and then silently dropped the source the
day one was set. A condition is now a frozen value: ``WallCondition.uniform(0.0)
== WallCondition.vacuum()`` (the drop rule: a zero weight is not stored;
``[M]`` section S1: equal, both realize a ``ZeroMap``, ``q.linf`` is
``0.0``), so there is no zero-source member of the family waiting to
receive a source, and a different source is a different condition. The
discipline that remains is the one the DSA row shows: **a consumer of
the realized operator, which carries the linear part only, must read the
declaration for the source and never infer it from the operator.** A
response-factor test alone would admit a sourced face, because its
:math:`R` *is* zero.

.. _bc-realizer-layer:

Layer 3 — the method realizer
-----------------------------

A declared wall condition
(:class:`~orpheus.transport.boundary_conditions.WallCondition`, an affine
pair of a response network and an inflow) or a deck describes the physics at
a face but is **not** by itself ready for consumption by a transport sweep.
The conversion from the method-agnostic declaration to a method-specific
object is the job of the method's **realizer table**, a
:class:`~orpheus.transport.boundary_conditions.RealizerTable`, and each
transport method owns exactly one: ``SN_REALIZERS``,
``DIFFUSION_REALIZERS``, ``CP_REALIZERS``, ``MOC_REALIZERS`` and
``MC_REALIZERS`` (:data:`~orpheus.sn.boundary.realizer.SN_REALIZERS` is the
one this page follows; the others are built the same way and live in
``orpheus/diffusion/boundary_realizer.py``, ``orpheus/cp/solver.py``,
``orpheus/moc/geometry.py`` and ``orpheus/mc/solver.py``).

**What a table is.** A table is keyed by the *unit* of a declaration and
filled by realizer classes that register by a class keyword
(``class SNSpecularChain(ResponseRealizer, table=SN_REALIZERS,
realizes=SpecularChain)``); each table is sealed after its last registration (on the
module's last line for SN, diffusion and MoC; at the end of the realizer
block for CP and MC), so anyone holding the table has executed every
registration. The units are
the two chains of the response network
(:class:`~orpheus.transport.boundary_conditions.SpecularChain`, the
permuting return, and
:class:`~orpheus.transport.boundary_conditions.DiffuseChain`, the Lambertian
return), the marker
:class:`~orpheus.transport.boundary_conditions.ZeroNetwork` (the empty
network, vacuum), the marker
:class:`~orpheus.transport.boundary_conditions.InflowTerms` (the source
channel :math:`q`, whose entry realizes to an array and never to an
operator), and any condition type that belongs to one method alone (the
diffusion zero flux). ``WallCondition`` itself is not a key: a wall is
realized by a **functor** on its response network, the image of the empty
network being the method's zero map, the image of ``s * SpecularChain`` and
``d * DiffuseChain`` being those chains' realizers at weights :math:`s` and
:math:`d`, and the image of a sum being the sum of the images, so a method
never writes a branch over the shared family. A deck (a mirror plane or a
translation wrap) meets the table's
:class:`~orpheus.transport.boundary_conditions.DeckRealizer`, at most one
per table (MoC registers none, so a deck there meets the unregistered door,
#654).

**One door.** :meth:`RealizerTable.realize
<orpheus.transport.boundary_conditions.realization.RealizerTable.realize>`
takes the declaration, the face's
:class:`~orpheus.transport.boundary_conditions.FaceSite` and the **method
space** (for SN an :class:`~orpheus.sn.mesh.method_space.SNMethodSpace`, the
container holding the :term:`quadrature` and the trace ladder the realizer
reads), calls
:meth:`~orpheus.transport.boundary_conditions.realization.RealizerTable.admit`
first and raises
:class:`~orpheus.transport.boundary_conditions.BoundaryRefused` iff the
admission is a refusal, then returns the method's operator type (for SN a
:class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`, a
:math:`\Gamma_+ \to \Gamma_-` operator carrying
:attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY` whose facts are
properties of its type). Admission and realization are one table read, so
the question "is this declaration admitted on this face" and the realization
cannot disagree (the admission tables and their axes are
:ref:`bc-method-realizability` and the refusal axes there). The per-method
hook is the :class:`~orpheus.transport.method.TransportMethod` Protocol: a
method-mesh declares its table as the class variable ``boundary_realizers``
and implements ``realize_boundary_role``, and the one shared
:func:`~orpheus.transport.method.realize_boundary` body owns the face loop
and the choice of the declaration (a deck, or the condition the face's tag is
bound to).

Why this third layer? Because the same affine wall is realized by
**structurally different** operators in each transport method, and the
table is where that difference is written down, once per method and per
unit. The image of the empty network (vacuum) and the realizer of each chain,
read off the five tables at this commit:

.. list-table:: What each method's table makes of a unit
   :header-rows: 1
   :widths: 12 22 22 22 22

   * - Method
     - ``ZeroNetwork`` (vacuum)
     - ``SpecularChain``
     - ``DiffuseChain``
     - Decks
   * - SN
     - ``ZeroMap``: the zero map :math:`\Gamma_+ \to \Gamma_-` on that
       face's two half-traces (the forward emits the zero of
       :math:`\Gamma_-`, the transpose the zero of :math:`\Gamma_+`; see
       :ref:`bc-domain-narrowing`)
     - ``ScaledOrdinatePermutation``, :math:`s\,P_\pi`; refused on the
       angular-resolution axis when the quadrature does not certify the
       mirror
     - ``RankOneReemission``, :math:`d\,(E\circ R_t)`; refused on a radial
       outer face (#653)
     - ``ScaledOrdinatePermutation`` at weight 1.0 (a mirror by
       ``ordinate_pairing``, a wrap by its crossing check)
   * - Diffusion
     - the albedo :math:`\mathcal{A} = 0`, the Marshak
       zero-incoming-current condition :math:`J^- = 0` on the scalar
       partial-current trace: the *same* declaration, a structurally
       different operator
     - :math:`w\cdot I` (tier 2: the permutation is annihilated by the
       :math:`P_1` projection)
     - :math:`w\cdot I` (tier 2)
     - mirror :math:`I`; translation refused (topology)
   * - CP
     - the vacuum closure, :math:`P_\infty = P_{\rm cell}` (escaped
       neutrons are lost), at the outer face only (the slab's left face and
       a hollow inner face are refused, #513)
     - refused (plumbing, #654)
     - the white closure, at weight 1.0 only (below one the albedo would be
       dropped, so it is refused, #654)
     - a mirror only at a slab's left face; no wrap
   * - MoC
     - refused (#654): MoC links specular tracks only
     - the track linking at the square cell's flat faces, weight 1 only;
       a named domain approximation (#547)
     - refused (#654)
     - no deck realizer (unregistered door, #654)
   * - MC
     - refused (#654): MC wraps its faces periodically
     - refused (topology on a solid cylinder, #547; plumbing at a Cartesian
       face, #654)
     - refused (#654)
     - the translation wrap; a mirror is refused (#654)

``[M]`` read at 471c4737 from the ``status`` and ``admit`` bodies of the
realizer classes in those five modules
(``grep -n 'status = \|realizes=' orpheus/sn/boundary/realizer.py
orpheus/diffusion/boundary_realizer.py orpheus/cp/solver.py
orpheus/moc/geometry.py orpheus/mc/solver.py`` lists them), not from the
carve plan. The "refused (#654)" cells are the capability the table
records as not built: the old string registries held a closed list per method
(CP the white wall and the vacuum, MoC the specular tracks, MC the periodic
wrap), and the tables keep that coverage and state the rest as refusals that
name the issue.

Splitting the realization out of the declaration makes each piece
independently testable and gives every method a single bolt-in point for its
boundary treatment: a method adopts the architecture by building and sealing
its own table (see :ref:`bc-cross-method-stubs`).

This layer used to be a ``BoundaryRealizer`` Protocol, one realizer class per
method that took a trace law and the method space, and a string-keyed
registry that named which laws each method admitted, with the mesh's
``realize_boundary_law`` arm as the per-method hook (SN first; diffusion
adopted it at #290 P3; the hook and the dissolution of the Wave-5 stubs were
#290 P7b). Wave 2 (2026-10-24)
retired that architecture with the trace laws it consumed; the present one is
:ref:`bc-architecture-of-record`, and each retired name maps to its
successor in :ref:`bc-wave2-crosswalk`.


.. _bc-extraction:

Boundary-condition extraction — :math:`B` as a sibling operator
===============================================================

The composite metric-correct G-adjoint that closes this extraction
narrative is documented at :ref:`g-adjoint` in
:doc:`/theory/foundations/operator_adjoint`.

Wave O step O.4a.2 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, three commits
``d7e1316`` / ``4c0ff96`` / ``2bdc66d``, 2026-06-03) made the realized
boundary condition :math:`B` a **first-class sibling** of the streaming +
collision operator :math:`L + C`. The canonical SN transport algebra
became

.. (vv-status rationale) Structural framing of the post-extraction SN
   loss operator. The verifiable claim — that the matvec/SI driver
   path with the realized ``B`` folded in agrees with the analytical
   infinite-medium balance and the homogeneous closed-form :math:`k_\infty`
   — is verified by the reflective convergence-equivalence gates
   catalogued below, not by this label directly.
.. vv-status: bc-extraction-loss-operator documented

.. math::
   :label: bc-extraction-loss-operator

   (L_{\rm full} + C - S - F - B)\,\psi \;=\; q,

acting on the **direct-sum transport state**

.. math::
   :label: bc-extraction-direct-sum-state

   V \;=\; V_{\rm bulk} \;\oplus\; V_{\rm inflow} \;\oplus\;
           V_{\rm outflow},

where :math:`V_{\rm bulk}` is the cell-centre angular flux
(:class:`~orpheus.transport.fields.angular_flux.AngularFlux`) and
:math:`V_{\rm inflow} \oplus V_{\rm outflow}` is the boundary trace
(:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`),
partitioned per face by the sign of :math:`\Omega\cdot\hat n` into the
inflow (:math:`\Omega\cdot\hat n < 0`) and outflow
(:math:`\Omega\cdot\hat n > 0`) ordinate slots (the
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` selectors
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`
/ :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`,
single source of truth — see :ref:`trace-spaces-doc`).

.. vv-status: bc-extraction-direct-sum-state documented

This is the realisation, for the boundary block, of the Wave T
prediction (:ref:`tensor-network-decomposition`): "Wave O typing must accept
non-SOTP summands." :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`
is the :math:`A_{ss}` leaf — a bespoke
:class:`~orpheus.numerics.operator.LinearOperator` carrying
:attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY`, NOT a
:class:`~orpheus.numerics.operator.TensorProductOperator`.


The block matrix
----------------

On :math:`V = V_{\rm bulk} \oplus V_{\rm boundary}` the operator families
occupy disjoint blocks of the :math:`2\times 2` block matrix, grouped by
the three :class:`~orpheus.numerics.operator.BlockRole` values
(:math:`b` = bulk row/column, :math:`s` = surface/trace row/column):

.. math::
   :label: bc-extraction-block-matrix

   \underbrace{
   \begin{bmatrix} A_{bb} & 0 \\ 0 & 0 \end{bmatrix}
   }_{C,\,S,\,N_{2n},\,F\;(\text{BULK})}
   \;+\;
   \underbrace{
   \begin{bmatrix} A_{bb} & A_{bs} \\ A_{sb} & A_{ss} \end{bmatrix}
   }_{L_{\rm full}\;(\text{FULL})}
   \;-\;
   \underbrace{
   \begin{bmatrix} 0 & 0 \\ 0 & A_{ss} \end{bmatrix}
   }_{B\;(\text{BOUNDARY})}

.. vv-status: bc-extraction-block-matrix documented

The position :math:`A_{ss}` is occupied by **both** :math:`L_{\rm full}`
and :math:`B` — with complementary triangle structure on the trace
splitting :math:`V_{\rm boundary} = V_{\rm inflow} \oplus V_{\rm outflow}`:

.. math::
   :label: bc-extraction-trace-blocks

   \underbrace{
   \begin{bmatrix} I & 0 \\ -T & I \end{bmatrix}
   }_{A_{ss}\ \text{of}\ L_{\rm full}}
   \qquad\text{vs.}\qquad
   \underbrace{
   \begin{bmatrix} 0 & R \\ 0 & 0 \end{bmatrix}
   }_{A_{ss}\ \text{of}\ B}
   \qquad\text{on}\quad
   \begin{bmatrix} \psi.\text{inflow} \\ \psi.\text{outflow} \end{bmatrix}.

.. vv-status: bc-extraction-trace-blocks documented

:math:`L_{\rm full}`'s trace–trace block is **unit-lower-triangular**,
with the identity on BOTH diagonal sub-blocks: the inflow row is the
carried identity :math:`I\cdot\psi.\text{inflow}` (it reads nothing
else), and the outflow row is the self-consistency defect's
stored-unknown identity :math:`I\cdot\psi.\text{outflow}` (design
correction 1 below; the per-row sign is free because
:math:`q.\text{outflow} \equiv 0` — the 2-D path spells the outflow
row with the opposite sign, same diagonal-bearing structure).
:math:`T` is the closure's **direct inflow→outflow transmission**: the
strictly-sub-diagonal coefficient the sweep chain hands the outflow
face in terms of the inflow face that seeded it — nonzero whenever a
direction's chain runs face-to-face (e.g. diamond difference on a
slab), zero when the chain terminates on the :math:`r = 0`
pole-regularity seed instead of a face, or under a pure-upwind
closure whose faces read cells only.

That identity diagonal is **load-bearing twice over** (Issue #298):

* It is what makes every trace row of :math:`L_{\rm full} + C`
  **diagonal-bearing**, so the augmented within-group operator is
  block lower-triangular in the ordering inflow → bulk → outflow —
  and fully triangular once the bulk is ordered by the sweep
  schedule. Its direct solve IS the sweep, literally **forward
  substitution**: read the given inflow, sweep the bulk, define the
  outflow.
* It is the diagonal the sibling :math:`-B` leans on: :math:`B`'s
  :math:`A_{ss}` (:math:`R` the realized per-face operator of the wall's response
  network, :ref:`bc-realizer-layer`) sits **strictly upper** on the same ordering —
  the inflow row reading the outflow column — the ONE up-edge that
  closes the boundary cycle. Vacuum kills it (:math:`B = 0`, pure
  forward substitution); SI/Krylov iterates it
  (:eq:`bc-extraction-two-residuals`); and :math:`B`'s low rank is
  what lets a Woodbury closure solve it in closed form instead
  (Issue #300).

(The block matrix above wrote :math:`L_{\rm full}`'s :math:`(s,s)`
entry as :math:`0` until Issue #298 — silently contradicting the role
table below, whose inflow-identity and outflow-defect rows are both
:math:`(s,s)` content. The :math:`0` was never load-bearing — no
downstream derivation consumed it — but the root-page triangularity
argument was about to.)

The three operator roles (Wave O block-role typing, Issue #208) read
off the block structure directly:

.. list-table:: Block occupancy by operator role
   :header-rows: 1
   :widths: 16 16 30 38

   * - Operator(s)
     - Role
     - Reads / writes
     - Block content
   * - :math:`L_{\rm full}`
       (:class:`~orpheus.sn.operators.streaming.StreamingOperator`,
       via :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` ``L+C``)
     - ``FULL``
     - Reads :math:`\psi.\text{bulk}` **and** the *given*
       :math:`\psi.\text{inflow}` trace; writes :math:`\psi.\text{bulk}`
       and the :math:`\psi.\text{outflow}` trace.
     - :math:`A_{bb}` (streaming) + :math:`A_{bs}` (inflow seeds the
       sweep) + :math:`A_{sb}` (sweep produces outflow) +
       :math:`A_{ss}` (the trace rows' unit-triangular structure,
       :eq:`bc-extraction-trace-blocks`). The
       **outflow row keeps the self-consistency defect**
       :math:`\psi.\text{outflow} - \text{streamed}`; the **inflow
       row carries the identity** :math:`I\cdot\psi.\text{inflow}`.
       **No BC logic.**
   * - :math:`C, S, F`
       (the collision multiplier :math:`C = M[\sigma_t]`
       — :class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`,
       :class:`~orpheus.transport.operators.scattering.ScatteringOperator`,
       :class:`~orpheus.transport.operators.fission.FissionOperator`)
     - ``BULK``
     - Bulk → bulk only.
     - :math:`A_{bb}` only; the boundary block is zero.
   * - :math:`B`
       (:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`)
     - ``BOUNDARY``
     - Maps :math:`V_{\rm outflow} \to V_{\rm inflow}` via the
       realized per-face operator :math:`\psi.\text{inflow} =
       B\,\psi.\text{outflow}`.
     - :math:`A_{ss}` only — strictly upper on the trace splitting
       (:eq:`bc-extraction-trace-blocks`); emits on the **inflow row
       ONLY** (see design correction 2 below).

The outer Krylov / SI loop drives **two residuals to zero**
simultaneously:

.. math::
   :label: bc-extraction-two-residuals

   \begin{aligned}
   r_{\rm inflow} &\;=\; \psi.\text{inflow}
                          \;-\; B\,\psi.\text{outflow}
                          \;-\; q.\text{inflow}
                    \;\longrightarrow\; 0
                    \quad (\text{boundary consistency}), \\
   r_{\rm outflow} &\;=\; \psi.\text{outflow}
                           \;-\; \text{streamed}(\psi.\text{bulk},
                           \psi.\text{inflow})
                    \;\longrightarrow\; 0
                    \quad (\text{outflow definition}).
   \end{aligned}

.. vv-status: bc-extraction-two-residuals documented

For vacuum :math:`B = 0` (no inflow); for any other wall :math:`B` is the
realized operator of that wall's response network, the reflection that the
pre-extraction sweep applied at the boundary edge: a scaled ordinate
permutation :math:`s\,P_\pi` for a specular return, the factored Lambertian
re-emission :math:`d\,(E\circ R_t)` for a diffuse one, and the sum of the two
images for a mixed wall, and for a deck face the unit permutation of its
mirror or its wrap (the wrap being :math:`B`'s only off-diagonal block)
(:ref:`bc-realizer-layer`).


The deleted keystone
--------------------

Before O.4a.2, the streaming sweep re-applied the boundary law
**inside one matvec call**: the backward (inward) sweep seeded its
boundary-face inflow from the forward sweep's own reflected outflow,

.. code-block:: python

   # PRE-O.4a.2 — the "keystone" (operator.py _compute_LpC, now DELETED):
   outflow_at_boundary = _sweep_direction(+1, pole_face_seed)  # forward
   inflow_full = bc_outer.apply(outflow_at_boundary.T)         # ← KEYSTONE
   outflow_at_inner = _sweep_direction(-1, inflow_full[...])   # backward

This single line **coupled bulk ↔ boundary within the matvec**: the
operator :math:`L` secretly contained the boundary reflection, so
:math:`L` was not a pure streaming operator and :math:`B` had no
independent existence. O.4a.2 **deletes the keystone**: the backward
sweep now reads the *given* outer inflow trace
(:math:`\text{face\_outer}`'s :math:`\mu < 0` ordinates) directly, so
one matvec call computes a pure-streaming residual with no BC
reflection. The reflective coupling moves entirely to the sibling
:math:`-B`:

.. code-block:: python

   # POST-O.4a.2 — bare, keystone-free (operator.py _compute_LpC):
   outflow_at_boundary = _sweep_direction(+1, pole_face_seed)  # forward
   outflow_at_inner   = _sweep_direction(-1, face_outer)       # GIVEN inflow

The curvilinear pole seed (``psi_view[:, :, 0, 0]``) survives the
deletion because it is the **r = 0 regularity condition**, NOT a
boundary condition — it reads the innermost cell flux, a geometric
case-split on ``curvature != "cartesian"``, never a ``bc.apply``.


.. _bc-extraction-design-corrections:

Three design corrections (what was tried and corrected)
-------------------------------------------------------

The extraction surfaced three subtle traps. All three are preserved
here per Cardinal Rule 3 so a future session re-deriving the block
matrix does not re-make them.

**1. Keep the outflow defect — NOT the raw outflow.**
The in-flight plan prose said :math:`L_{\rm full}` should emit the
*raw* streamed outflow on the outflow row. This is **wrong**.
:math:`\psi.\text{outflow}` is a *stored unknown* that the sibling
:math:`-B` reads as its input (:math:`B\,\psi.\text{outflow}`).
Emitting the raw outflow would make the outflow row
:math:`-\,\text{streamed}` (an off-diagonal-only row with no diagonal
on :math:`\psi.\text{outflow}`), which **singularises** that row: the
:math:`A_{ss}` outflow-column diagonal disappears and :math:`-B` is no
longer a well-posed sibling. The fix is to keep the row as the
self-consistency defect :math:`\psi.\text{outflow} - \text{streamed}`
— the identity :math:`I\cdot\psi.\text{outflow}` diagonal stays on the
outflow row, and the outflow-definition residual
:math:`r_{\rm outflow}` of :eq:`bc-extraction-two-residuals` is the
quantity the outer loop drives to zero. Keeping the row as
``computed − stored`` also makes the vacuum path **bit-identical** to
the pre-extraction matvec (the per-row sign is free because
:math:`q.\text{outflow} \equiv 0` — the outflow trace is a pure
definition with no source).

**2.** :math:`B` **must not emit on the outflow row** — first solved by
projection, since B3.2 solved by **typing**.
At the time of the extraction the realized per-face law was a
**full-face operator**: a specular
:class:`~orpheus.numerics.elementary_operators.PermutationOperator` for reflective,
the then-welded ``AngularAverageOperator`` for
white. Its permutation mapped the input's *inflow* slots onto the
*output's outflow* slots (a spurious :math:`R\cdot\psi.\text{inflow}`),
because the permutation was defined on the whole face, not just the
:math:`A_{ss}` :math:`V_{\rm outflow} \to V_{\rm inflow}` map. In the
legacy sweep this was harmless — the sweep only ever read the
inflow slots of ``bc.apply(face)``, discarding the outflow output.
But as a sibling :math:`-B` on the direct-sum state, a non-zero
outflow emission corrupts the outflow-definition residual
:math:`r_{\rm outflow}` (which must carry **no** :math:`B` term).
*Empirically confirmed before the fix*: the outflow slots carried
nonzero :math:`R\cdot\psi.\text{inflow}`.

The extraction's fix was an **output projection** —
``_reflect_trace`` wrote only the ``inflow_indices_for_face`` slots of
the law's full-face image. Campaign phase **B3.2** replaced that with
the honest domain: the law is typed :math:`\Gamma_+ \to \Gamma_-`, so
the composition is :math:`\iota_-\circ\text{law}\circ\gamma_+` and the
spurious emission is **unrepresentable** rather than projected away
(:ref:`bc-domain-narrowing`, :eq:`bc-face-action-narrowed`). The
diagnosis in this correction is what B3.2 acted on; only the remedy
changed, from a mask at the consumer to a domain at the producer.

.. warning::

   **The transpose scatters over** :math:`\Gamma_+` **— never over**
   :math:`\Gamma_-`. This trap **survives the narrowing** and is the
   one piece of correction 2 that is still live discipline rather than
   history. Because :math:`(\iota_-\circ\text{law}\circ\gamma_+)^{\mathsf T}
   = \iota_+\circ\text{law}^{\mathsf T}\circ\gamma_-`, the transpose's
   *input* is masked to the forward's codomain :math:`\Gamma_-` and its
   *output* lands on :math:`\Gamma_+`. Scattering that output over
   :math:`\Gamma_-` instead — "projecting :math:`\text{law}^{\mathsf T}`
   onto the law's own codomain", which reads like the natural mirror of
   the forward — extracts the operator's **diagonal block**. For vacuum
   that spells a spurious :math:`+1` where the forward is the zero map.
   It was caught **only** by the A2a grid-reciprocity arm on the
   heterogeneous-**vacuum** sphere: off-diagonal permutation laws are
   bit-identical under either spelling, so every reflective fixture
   stayed green over the wrong one. (This Euclidean ``apply_transpose``
   is the un-weighted shadow of the metric-correct Hilbert adjoint
   ``B.H`` under :math:`|\Omega\cdot\hat n|\,w`; the two are separate —
   see :ref:`g-adjoint`.)

**3. The bare sweep seeds inflow from** :math:`\text{rhs.boundary}`,
**not** :math:`\text{initial\_guess.boundary}`.
Under the extraction the WDD sweep ``(L+C).solve`` is **bare** (see
:ref:`bare-sweep-extraction` in
:doc:`/theory/methods/sn/curvilinear_one_group`): it reads
the seeded inflow trace directly instead of re-applying ``bc``.
:meth:`StreamingCollisionOperator._solve_timed_full_field <orpheus.sn.operators.streaming.StreamingCollisionOperator._solve_timed_full_field>`
must therefore seed the sweep's boundary buffer from
:math:`\text{rhs.boundary}` — the *boundary source*
:math:`q.\text{boundary} + B\,\psi.\text{outflow}` — **not** from the
iterate ``initial_guess.boundary`` (the retired partner-flux carrier).
The ``initial_guess`` still threads the bulk Carlson warm start;
only the boundary seed moved.


.. _bc-extraction-variadic-driver:

The honest :math:`L+C-S-N_{2n}-B` driver via variadic couplings
---------------------------------------------------------------

The within-group inner solve no longer hands the drivers a fixed
:math:`(A, S, F)` operator *triple*. Wave O step O.2a generalised both
:class:`~orpheus.numerics.iteration.SourceIteration` and
:class:`~orpheus.numerics.iteration.KrylovAcceleration` to the
**variadic** shape :math:`\text{Driver}(A_{\rm resolvent},\,*\text{gains})`:
one invertible resolvent :math:`A` plus a homogeneous bag of lagged
coupling operators :math:`g_i`. The two consume the gains identically —

.. math::
   :label: bc-extraction-variadic-matvec

   \text{matvec} \;=\; A.\text{apply} \;-\; \sum_i g_i.\text{apply}
   \,,\qquad
   \text{rhs} \;=\; q_{\rm ext} \;+\; \sum_i g_i.\text{apply}\,.

.. vv-status: bc-extraction-variadic-matvec documented

The driver is now **problem-type-agnostic**: it sees only the resolvent
operator and a bag of operators it must lag. (Since #226 taxonomy step 3,
:class:`~orpheus.numerics.iteration.SourceIteration` receives that
resolvent **already inverted** — the pre-built inverse operator it
*applies* — while :class:`~orpheus.numerics.iteration.KrylovAcceleration`
keeps the *forward* resolvent for its GMRES matvec; the
:eq:`bc-extraction-variadic-matvec` form above is the Krylov matvec, and
the rhs term is the shared source assembly. See
:ref:`inverse-application-driver`.) *Which* leaves are gains is a
**posing-layer** decision, not an iteration-layer one (see
:ref:`eigenvalue-posing`) — the gains are exactly the posing's coupling
terms.

For the SN **k-eigenvalue** within-group inner the posing's couplings
are the bulk scattering :math:`S` and the boundary reflection
:math:`B`; the fission :math:`F` is zero within-group (it enters as the
external source :math:`q_{\rm ext}` per the eigenvalue
outer / within-group split, Lewis & Miller §6.4). So the within-group
loss decomposition is the honest

.. math::
   :label: bc-extraction-within-group-decomposition

   (L+C,\; S,\; N_{2n},\; B)
   \quad\Longrightarrow\quad
   \underbrace{(L+C).\text{apply} - S.\text{apply} - N_{2n}.\text{apply}
               - B.\text{apply}}_{\equiv\,(L+C-S-N_{2n}-B)\,\psi}
   \,,\qquad
   \text{rhs} = q_{\rm ext} + S.\text{apply}(\psi)
              + N_{2n}.\text{apply}(\psi) + B.\text{apply}(\psi)

.. vv-status: bc-extraction-within-group-decomposition documented

assembled once by the single-source-of-truth builder
:func:`~orpheus.sn.coupled_system.build_within_group_system`, which returns
the frozen :class:`~orpheus.sn.coupled_system.WithinGroupSystem` record —
the loss grid together with the bound leaves it is the signed sum of
(:class:`~orpheus.sn.coupled_system.SNLossFactors`).  The **splitting**
:math:`A = M - N` (Hackbusch 2016 §11 — block partitionings; a
*splitting*, **not** a *regular* splitting in Varga's sense, so the
comparison theorem does **not** bound the boundary Gauss-Seidel rate
against Jacobi's — :ref:`sn-boundary-gs-not-regular`) is not on that
record: since 2026-09-13 it is a Strategy VALUE labelled from those
leaves, :class:`~orpheus.sn.splitting.Splitting`
(:ref:`sn-splitting-is-a-strategy-value`). On a seedless
(slab / cylinder /
Cartesian) mesh the Jacobi value degrades to exactly this triple: its
:attr:`~orpheus.sn.splitting.Splitting.implicit`
is :math:`M = (L+C)` — the invertible resolvent
(:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator`, ``.solve`` = the WDD
sweep) — and its :attr:`~orpheus.sn.splitting.Splitting.explicit`
pieces are :math:`N = (S,\ N_{2n},\ B_a)`,
the three lagged
couplings the driver applies: the bulk scattering gain
(:class:`~orpheus.transport.operators.scattering.ScatteringOperator`,
:attr:`block_role <orpheus.numerics.operator.BlockRole>` ``BULK``), the
bulk :math:`(n,2n)` emission gain
(:class:`~orpheus.transport.operators.n2n.N2NOperator`, ``block_role``
``BULK`` — a third member since CS4c step 3;
:ref:`the two collision gains <operator-algebra-two-gains>`) and the
boundary reflection gain
(:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`,
``block_role`` ``BOUNDARY``). On a *carrying* sphere the record poses the
ψ½ System B as the second row of a coupled :math:`M - N` block grid — the
full 2×2 coupled block operator documented in
:ref:`coupled-block-operator` (its starting-direction physics in
:ref:`sn-direct-seed-solve` in
:doc:`/theory/methods/sn/curvilinear_one_group`);
the former ``_within_group_triple`` / ``_lagged_gains`` construction pair
retired into this one builder at B.2d. The :math:`B\,\psi.\text{outflow}`
term lands on :math:`\text{rhs.boundary}`, which the bare ``(L+C).solve``
sweep reads as the inflow seed (:ref:`bare-sweep-extraction` in
:doc:`/theory/methods/sn/curvilinear_one_group`).

**This retires the transitional** :math:`S + B` **fold.** The
predecessor packed the boundary reflection into the *middle slot* of
the fixed triple by returning a summed operator
:math:`S + B` — the now-deleted ``SNSolver._scattering_with_boundary_op``
property. The honest composition keeps :math:`S` and :math:`B` as two
**separate first-class gains**.

**Why variadic — the fixed triple encoded a false posing distinction.**
:math:`S`, :math:`F` and :math:`B` are *homogeneous* in the driver:
each is subtracted in the matvec and summed in the rhs, exactly as
:eq:`bc-extraction-variadic-matvec` shows. The fixed :math:`(A, S, F)`
triple gave :math:`S` and :math:`F` named slots the *resolvent layer*
never uses — it was encoding a posing-layer role assignment (which
operator is loss, which is the eigen-operator) at the iteration layer,
where it does not belong. Collapsing the triple to a homogeneous
:math:`*\text{gains}` bag moves the role distinction back to the
posing layer (its proper home) and lets a fourth gain (a future
:math:`B`-trace term, an :math:`\alpha`-time term) slot in as a data
addition rather than a new named slot. Existing positional
:math:`(A, S, F)` callers stay source-compatible — ``gains = (S, F)``.

**Why** :math:`B` **is a SEPARATE gain, not folded into** :math:`S`.
Two structural reasons forbid the old fold:

#. **The adjoint metric lives on the trace.** :math:`B` acts only on the
   boundary-trace block of the composite carrier (its
   :attr:`domain <orpheus.sn.operators.boundary.SNBoundaryOperator.domain>`
   and ``codomain`` are ``problem.full_field_space``, the space
   :math:`L`, :math:`C`, :math:`S` and :math:`F` report, and its bulk output
   is zero), and the cosine-weighted
   :math:`|\Omega\cdot\hat n|\,w` adjoint metric (Wave O step O.2 — the
   codomain inner product of :math:`L`'s boundary-trace block) lives on
   that **trace** block, not the bulk. Folding :math:`B` into the
   bulk :math:`S` would erase the trace typing the adjoint
   ``.H`` reads.
#. :math:`B` **cannot join the** :math:`L+C` **preconditioner.** A
   generic :class:`~orpheus.numerics.operator.OperatorSum` carries **no**
   direct-sweep ``solve`` — its inverse action is the *iterative*
   :class:`~orpheus.numerics.green_operator.GreenOperator` splitting
   (:ref:`green-operator`), not the ``(L+C)``
   :class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` subclass's
   :math:`O(N\cdot N_{\rm cells})` forward-substitution sweep. So folding
   :math:`B` into an :math:`L + C - B` sum would **demote** the cheap
   direct sweep the SI step and the Krylov preconditioner depend on to an
   iterated solve. :math:`B` must stay a *gain* (lagged, applied) — never
   a summand of the resolvent.

The old fold type-checked at the time **only because**
:attr:`ScatteringOperator.domain` *was* ``None`` (the pre-W-D bulk
operators inherited the
:class:`~orpheus.numerics.operator.LinearOperator` default). The
:class:`~orpheus.numerics.operator.OperatorSum` domain-compatibility
check fires only when both operands declare non-``None`` domains that
differ; with :math:`S` untagged the check skipped, so the
trace-typed :math:`B` summed silently with the bulk :math:`S`. The
structural reason the fold stays gone is the **variadic-driver
redesign** below — :math:`B` is a lagged *gain*, never a resolvent
summand — not a typing tripwire: P4.5 W-D gave :math:`S` the
**composite** :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
(the *same* instance :math:`L`/:math:`C`/:math:`B` carry, so the
within-group ``(L+C) - S`` guard *validates* the ``- S`` arm), NOT a
bulk-only :math:`V_{\rm bulk}` domain, so an ``OperatorSum`` of the
composite-typed :math:`S` and the composite-typed :math:`B` would
still compose — the once-envisioned "bulk-S ≠ trace-B" rejection seam
(see :ref:`bc-extraction-scope-future`) was **not** the shape W-D
landed.

.. note::

   The drivers' :class:`~orpheus.numerics.iteration.KrylovAcceleration`
   matvec :eq:`bc-extraction-variadic-matvec` and the
   :class:`~orpheus.numerics.iteration.SourceIteration` rhs are now the
   *honest* :math:`(L+C-S-N_{2n}-B)\,\psi` and
   :math:`q_{\rm ext}+S\psi+N_{2n}\psi+B\psi`
   — the reassociation :math:`A-(S+B)\to(A-S)-B` is documented as a
   **principled-equivalence** change in
   :ref:`bc-extraction-numerical-evidence` (criterion 3 of the
   ``vv-principles`` bit-identity-vs-principled-equivalence gate), not
   a bug.


.. _bc-extraction-two-routes:

The two :math:`-B` delivery routes
----------------------------------

The same :math:`-B` coupling reaches the sweep two ways, both calling
the **identical**
:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` (single
source of truth, Cardinal Rule 2):

.. list-table:: The two delivery routes for :math:`-B`
   :header-rows: 1
   :widths: 22 40 38

   * - Route
     - Mechanism
     - Used by
   * - **Variadic gain**
       (:func:`~orpheus.sn.coupled_system.build_within_group_system`
       returns :math:`B` as a gain)
     - :math:`B` is one of the ``*gains`` the variadic driver lags:
       the matvec subtracts :math:`B.\text{apply}`, the SI rhs adds it
       (:eq:`bc-extraction-variadic-matvec`).
       :math:`B\,\psi.\text{outflow}` lands on
       :math:`\text{rhs.boundary}`, which the bare ``(L+C).solve``
       sweep reads as the inflow seed.
     - The eigenvalue SI inner driver
       (:meth:`SNSolver._solve_source_iteration <orpheus.sn.solver.SNSolver._solve_source_iteration>`),
       the eigenvalue Krylov inner
       (:meth:`SNSolver._solve_krylov <orpheus.sn.solver.SNSolver._solve_krylov>`),
       and both fixed-source paths
       (:func:`_solve_fixed_source_si <orpheus.sn.solver._solve_fixed_source_si>` /
       :func:`_solve_fixed_source_krylov <orpheus.sn.solver._solve_fixed_source_krylov>`)
       — every solve that routes through a driver.
   * - **Masked additive reflect**
       (:meth:`SNMaskedBoundaryOperator.reflect_rows_inplace
       <orpheus.sn.operators.boundary.SNMaskedBoundaryOperator.reflect_rows_inplace>`,
       on the mask's rows)
     - **Accumulates** :math:`B\,\psi.\text{outflow}` onto the mask's
       inflow rows, in place on the boundary buffer
       (``bf[f][rows] += (B·bf)[f][rows]``), through the same
       :class:`SNBoundaryOperator` the mask wraps.  Additive, **not** a
       whole-face assignment — the row it completes is the
       inhomogeneous forward-substitution row :math:`z_{\rm in} =
       y_{\rm row} + (Bz)_{\rm row}`, whose seed :math:`y_{\rm row}`
       an overwrite would drop (:ref:`the rejected whole-face
       overwrite <gs-whole-face-overwrite-rejected>`).
     - The octant-group Gauss–Seidel forward substitution: the reified
       :math:`M = (L+C-B_{\rm lower})` passes
       ``reflect=self.lower.reflect_rows_inplace`` into its scheduled
       walk (:class:`~orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator`),
       so each just-swept octant group's reflective faces feed the next
       one — the inter-group reflect of :ref:`the reified splitting
       matrix <si-gauss-seidel-reification>`.
       Through the **Jacobi** split's ``upper`` half — which IS the
       full-inflow mask, because a Jacobi schedule reflects no face
       in-sweep (``[M]`` 4/4 geometries: ``lower`` is empty and
       ``upper`` carries every inflow row of every face) — the same
       verb is also the sweep-tier gates' *inter-sweep* reflect
       (``tests/gates/sn/_test_helpers.py::reflect_outflow_into_inflow``).

       ⛔ **Until 2026-09-07 this row named a whole-trace ASSIGNMENT
       verb instead** — ``SNBoundaryOperator.reflect_inflow_inplace``,
       which filled each face's inflow slots with
       :math:`B\,\psi.\text{outflow}` outright before a bare sweep.  It
       had lost its last production caller at #448 (2026-09-06): the
       direct fixed-source SI loop moved onto the variadic driver at
       Wave O O.2a, and the eigenvalue reconstruction sweep in
       :func:`solve_sn <orpheus.sn.solver.solve_sn>` became one step of
       that same driven map (:ref:`sn-finalize-one-step`).  CS4c step 6
       item 6.5 retired it, with its trace-only leaf
       ``reflect_into_inflow``, rather than keeping a second verb for a
       test-tree consumer: *zero the inflow rows, then reflect
       additively on the full-inflow mask* reproduces the assignment
       **bit-for-bit** — ``[M]`` ``np.array_equal`` on 40/40 seeds ×
       4/4 geometries against a NON-zero-inflow buffer, and dropping
       the zeroing moves the answer by :math:`O(1)` (the positive
       control that the reading is not a zero-inflow artefact).  So the
       semantics survived without the surface; the gate is
       ``tests/gates/sn/operators/test_reflect_helper_reexpression.py``.

The masked reflect is **not** a fold of :math:`B` into :math:`S`: it is
the trace-only :math:`A_{ss}` action of the *same* :math:`B`, restricted
to a subset of that operator's block ROWS and expressed on the boundary
trace alone. Both routes therefore deliver the identical :math:`-B`
coupling, and cannot drift, because both descend from
:meth:`SNBoundaryOperator._reflect_trace <orpheus.sn.operators.boundary.SNBoundaryOperator>`
(:ref:`bc-extraction-reflect-trace`).  Under the Gauss-Seidel schedule
they are not even alternatives but the two halves of ONE splitting
:math:`(L+C-B) = M - B_{\rm upper}`: the masked reflect carries
:math:`B_{\rm lower}` *inside* the reified forward :math:`M`, the
variadic gain carries :math:`B_{\rm upper}` lagged
(:ref:`si-gauss-seidel-reification`).


.. _bc-extraction-reflect-trace:

The trace-only :math:`A_{ss}` core — ``_reflect_trace``
-------------------------------------------------------

:math:`B` is the :math:`A_{ss}` block :math:`V_{\rm outflow} \to
V_{\rm inflow}`: it maps the *outflow* trace to the *inflow* trace.
Both delivery routes ultimately need the same per-face action —
restrict the face slot to :math:`\Gamma_+`, apply that face's realized
operator (the scaled ordinate permutation
:class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation` for a
specular wall, the factored Lambertian
:class:`~orpheus.sn.boundary.realized_operators.RankOneReemission` for a diffuse
one, the zero map :class:`~orpheus.sn.boundary.realized_operators.ZeroMap` for
vacuum), and
scatter the image back over :math:`\Gamma_-`. To guarantee
they cannot drift, that action is the single
:meth:`SNBoundaryOperator._reflect_trace <orpheus.sn.operators.boundary.SNBoundaryOperator>`
core, and **every** public route descends from it — there are exactly
two families.  The ``_apply_faces`` lift onto a zero-bulk
:class:`~orpheus.transport.full_field.FullField` carrier serves
:meth:`B.apply <orpheus.sn.operators.boundary.SNBoundaryOperator.apply>`,
its Euclidean transpose, and the masked operator's row-restricted
``apply``; the in-place additive
:meth:`SNMaskedBoundaryOperator.reflect_rows_inplace
<orpheus.sn.operators.boundary.SNMaskedBoundaryOperator.reflect_rows_inplace>`
calls the core directly, because it writes onto a caller's trace buffer
rather than returning a field.

⛔ **Until 2026-09-07 the sentence above named a different second
consumer** — the trace-only *leaf* ``B.reflect_into_inflow``, added at
Wave O step O.2a (commit ``8563f4b``) with its mutating façade
``reflect_inflow_inplace``, and retired together with it at CS4c step 6
item 6.5.  The no-drift argument is unchanged in kind and stronger in
reach: the second consumer is now the verb production actually binds
(the Gauss-Seidel resolvent's row update), not a leaf whose only
remaining caller was a test-tree helper.

**Why the leaf existed, and why that reason expired.**
:meth:`B.apply <orpheus.sn.operators.boundary.SNBoundaryOperator.apply>` operates on a :class:`~orpheus.transport.full_field.FullField`
(zero bulk, trace populated) — the timeless, history-blind operator carrier
(:meth:`SNBoundaryOperator.apply <orpheus.sn.operators.boundary.SNBoundaryOperator.apply>`
is the base arrow ``FullField -> FullField``; the comonad lives on the
driver), the bulk only a carrier to reach the :math:`A_{ss}` boundary block. The pre-extraction direct helper
fabricated a throwaway zero-bulk field purely to call ``B.apply`` and
then discarded the (zero) bulk output.  ``reflect_into_inflow`` removed
that probe: it took a bare
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
trace and returned the boundary-only
:class:`~orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink`
directly.  The saving was real for a *whole-trace* caller — and by 2026-09-07
there was none: the only trace-shaped consumer left was the sweep-tier gates'
inter-sweep reflect, which the masked verb serves through the Jacobi split's
``upper`` half.  A verb whose whole justification is an ergonomic saving for
callers that no longer exist is a second spelling of one action, so item 6.5
retired it rather than keeping the two in step.

**The full-field route's input contract** (CS4c step 6 item 6.3,
2026-09-07).  Because ``B.apply`` is reached through a composite carrier,
``_apply_faces`` parses that carrier before touching it, and it does so
through the SAME body its :math:`L` and :math:`L+C` siblings use —
:meth:`FullField.require_member
<orpheus.transport.full_field.FullField.require_member>`: a ``TypeError``
naming the refusing surface and the carrier it wanted for a foreign object
(until this item ``_apply_faces`` read ``psi.interior`` unguarded and leaked
a raw ``AttributeError``), a ``ValueError`` carrying the greppable
*space-content invariant* vocabulary for a content mismatch.  The parse is
keyed on the operator's **carrier mesh**, not on its bound end, and
deliberately so: :math:`B_a` is fed the windowed *moment* iterate as well as
the angular composite its ends name (``[M]`` 59 / 58 / 47
``HarmonicMomentFlux``-interior composites per 2-D windowed solve), and a
bound-end comparison would refuse every one of them.  ⚠ The guard is
tagged ``ELEGANCE-DEBT[guard]`` (`#457
<https://github.com/deOliveira-R/ORPHEUS/issues/457>`_) — it is a
protection, not the target state.  It retires when every leaf is bound on
the end it acts on, after which an alien carrier cannot be typed at all and
the admission is the ordinary composability guard.

Keeping :math:`B`'s emission off the outflow row is load-bearing: as
the sibling :math:`-B` reading the *whole* boundary block, a non-zero
outflow emission corrupts the outflow-definition residual
:math:`r_{\rm outflow}` (which must carry **no** :math:`B` term —
:ref:`bc-extraction-design-corrections`). The extraction achieved that
by *projecting* a full-face law's image onto the inflow rows; since
campaign phase **B3.2** it is achieved by *typing* — the realized face operator
:math:`R`'s domain is :math:`\Gamma_+` and its codomain is :math:`\Gamma_-`, so
:meth:`_reflect_trace <orpheus.sn.operators.boundary.SNBoundaryOperator>`
composes :math:`\iota_-\circ R\circ\gamma_+` and there is no
wider image to project (:ref:`bc-domain-narrowing`). The transpose
composes :math:`\iota_+\circ R^{\mathsf T}\circ\gamma_-`, so it
masks its INPUT to :math:`\Gamma_-` and scatters its OUTPUT over
:math:`\Gamma_+` — the asymmetry that
:ref:`correction 2's warning <bc-extraction-design-corrections>`
explains and that no reflective fixture can detect.


.. _bc-extraction-scope:

Scope — both 1-D and 2-D are now bare (O.4b complete)
-----------------------------------------------------

O.4a.2 made the **1-D** sweep bare (slab / sphere / cylinder). Step
**O.4b** then made the **2-D Cartesian wavefront sweep bare as well**
(both :func:`~orpheus.sn.loss_representation._sweep_jacobi` and the 2-D matvec
:meth:`StreamingOperator._apply_2d_cartesian <orpheus.sn.operators.streaming.StreamingOperator>`):
the intra-sweep ``bc.apply`` is **gone** for every geometry. The
octant-incoming face edge is seeded from the *given* inflow trace and
the reflective coupling :math:`\psi.\text{inflow} = B\,\psi.\text{outflow}`
is delivered by the sibling :math:`-B`
(:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`) for the
2-D trace exactly as for the 1-D trace. The 2-D matvec emits the
boundary block as an active-trace residual (outflow slots carry the
self-consistency defect ``streamed − ψ.outflow``; inflow slots carry
the identity ``ψ.inflow``), wired into the composed Krylov matvec as
the boundary gain :math:`B` of
:func:`~orpheus.sn.coupled_system.build_within_group_system`.
The interior face fluxes the bare 2-D sweep + matvec propagate are the
interior 1-cochain :math:`C^1_{\rm int}` — carried since S6.4(f) by the
rolling ``_MovingFrontier`` front (the ``WavefrontFlux`` type that named
it through #205 Phase 5 is retired; see :ref:`wavefront-flux-cochain`).

The dispatch is still guarded by a **single predicate** so the two
geometry paths cannot drift: ``problem.reduced is not None`` is the
**same** predicate the representation dispatch
(:func:`~orpheus.sn.loss_representation.default_for`, via each
representation's ``supports``) reads to select the 1-D scan
(:class:`~orpheus.sn.loss_representation.CumprodScan`) vs the 2-D
wavefront body, and the **same**
predicate the two entries' reconstructions read.  ⛔ This clause named
*"the direct-helper guards … before calling
``_reflect_outflow_into_inflow``"* until 2026-09-06: neither entry calls a
direct helper any more (the fixed-source SI loop stopped at Wave O O.2a,
the eigenvalue finalize at #448 — :ref:`sn-finalize-one-step`), so the
predicate now selects a *representation*, which was always what it was
about.  Both branches are now bare-sweep + sibling :math:`-B`; the predicate
selects the *fold shape* (1-D parallel-prefix scan vs 2-D wavefront
DAG), **not** a bare-vs-bc-in-sweep distinction.


.. _bc-extraction-scope-future:

Closed typing-completion seam — :attr:`ScatteringOperator.domain` (P4.5 W-D)
----------------------------------------------------------------------------

This Wave-O typing-completion **landed in P4.5 W-D** (commit
``0610b39``); it is recorded here because it was a documented seam at
the time of Wave O's close-out. The Wave-O framing envisioned giving
:class:`~orpheus.transport.operators.scattering.ScatteringOperator`
(and the other bulk leaves) a **bulk** :math:`V_{\rm bulk}` domain so
that :class:`~orpheus.numerics.operator.OperatorSum` would **reject** a
re-introduced :math:`S + B` fold (the domain-compatibility check
throwing ``IncompatibleOperatorComposition`` on a bulk :math:`S`
summed with a trace :math:`B`).

W-D closed the seam, but with a **different and stronger** choice: the
bulk leaves :math:`C`/:math:`S`/:math:`F` carry the **composite**
:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` — the
*same* instance :math:`L`/:math:`B` advertise — not a bulk-only
:math:`V_{\rm bulk}` space. The motivation shifted from the negative
"reject :math:`S + B`" tripwire to the positive **within-group
composition guard**: every operand of the within-group loss
:math:`(L + C) - S` now reports the same composite domain, so the
:class:`~orpheus.numerics.operator.OperatorSum` guard **validates** the
build (equal domains AND codomains) on every solve instead of silently
skipping a ``None``-spaced summand (W-D — see
:ref:`g-adjoint`). A consequence is that the original
"bulk-:math:`S` ≠ trace-:math:`B`" rejection no longer applies: a
composite-typed :math:`S` and a composite-typed :math:`B` speak the
same space, so an ``OperatorSum`` of the two would compose. The reason
the :math:`S + B` fold stays gone is the **variadic-driver redesign**
(:ref:`bc-extraction-variadic-driver`) — :math:`B` is a lagged gain,
never a resolvent summand — not a typing rejection.

The space-anonymous ``domain = None`` survives only for the **bare /
test constructor** (a :class:`ScatteringOperator` built without
``from_solver_data``'s ``full_field_space=`` thread): then the guard
skips that operand, preserving the legacy backward-compatible contract
for direct callers.


.. _bc-extraction-2d-si-krylov-twin:

The 2-D Cartesian eigenvalue SI inner is the geometry-agnostic twin of Krylov
-----------------------------------------------------------------------------

Because the variadic :math:`-S - N_{2n} - B` gains ride the **bare** sweep
for
every geometry (above), the two within-group eigenvalue inner solvers
are **structural twins** — they share every operator and every
reduction, differing only in the iteration driver. This holds for 2-D
Cartesian exactly as it does for slab / sphere / cylinder, so a 2-D
Cartesian eigenvalue problem solves through **both** inner solvers:

- :meth:`SNSolver._solve_source_iteration <orpheus.sn.solver.SNSolver._solve_source_iteration>`
  — the source-iteration inner, the :func:`~orpheus.sn.solver.solve_sn`
  **default** ``inner_solver="source_iteration"`` for *every* geometry,
  driven by :class:`~orpheus.numerics.iteration.SourceIteration`;
- :meth:`SNSolver._solve_krylov <orpheus.sn.solver.SNSolver._solve_krylov>`
  — the Krylov inner, opt-in ``inner_solver="krylov"``, driven by
  :class:`~orpheus.numerics.iteration.KrylovAcceleration`.

The two inners are identical except for that driver. Both build the
same composite right-hand side
(:meth:`AngularSourceSink.from_isotropic <orpheus.transport.fields.angular_source_sink.AngularSourceSink.from_isotropic>`
bulk + ``AngularBoundarySourceSink.zeros(problem.angular_trace)``
boundary inside a
:class:`~orpheus.transport.timed_full_field.TimedFullField`), the same
loss decomposition (the resolvent :math:`L + C` from
:class:`~orpheus.sn.operators.streaming.StreamingOperator` + the collision
multiplier :math:`C = M[\sigma_t]`
(:class:`~orpheus.transport.operators.multiplication_operator.MultiplicationOperator`),
plus the scattering
gain :math:`S` and the boundary reflection gain :math:`B` from
:func:`~orpheus.sn.coupled_system.build_within_group_system`;
zero within-group fission), and the
same :meth:`integrate_angular <orpheus.transport.fields.angular_flux.AngularFlux.integrate_angular>`
angular reduction. Neither driver carries any geometry dependence.

The reflective coupling reaches both drivers on the **bare** 2-D
wavefront sweep through the sibling :math:`-B` (the **variadic-gain**
route of :ref:`bc-extraction-two-routes`), never through an in-sweep
``bc.apply``. The :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`
is natively **four-face** (``xmin`` / ``xmax`` / ``ymin`` / ``ymax``)
and is the *same* operator the 2-D Krylov path uses — there is no
separate per-geometry boundary closure.

.. admonition:: The "B1'' face block" is retired legacy
   :class: note

   The 2-D Cartesian eigenvalue path was historically described as
   needing a distinct "B1'' face block" that was "1-D-only", which is
   why the source-iteration inner was once guarded against 2-D meshes.
   That guard is **removed**. "B1''" was never a code symbol — it was a
   1-D boundary-closure *name* in docstrings and comments, fully
   superseded by the L2
   :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` +
   :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`
   bare-boundary architecture (O.4a.2 / O.4b above), which realises the
   boundary handling for *all* geometries. The 2-D path never required
   a separate 1-D-only block. Because
   :func:`~orpheus.sn.solver.solve_sn` defaults to the source-iteration
   inner for every geometry, the now-removed guard had meant the
   **default** 2-D Cartesian eigenvalue entry point raised; the carve
   restores it.

**Numerical evidence (SI ≡ Krylov ≡ closed-form** :math:`k_\infty`
**).** The twin is pinned at the production
:func:`~orpheus.sn.solver.solve_sn` entry (not a hand-rolled power
loop) by ``tests/gates/sn/eigenvalue/test_keff_2d.py::TestSIKrylov2DEquivalence``:

.. list-table:: 2-D Cartesian SI/Krylov verification (Wave O step #208)
   :header-rows: 1
   :widths: 38 34 28

   * - Leg
     - Reference (pillar)
     - Result
   * - Default-entry homogeneous (1G / 2G / 4G)
       (:func:`test_default_entry_hits_kinf <tests.gates.sn.eigenvalue.test_keff_2d.TestSIKrylov2DEquivalence.test_default_entry_hits_kinf>`)
     - Closed-form :math:`k_\infty = \lambda_{\max}(A^{-1}F)`
       (1g → 1.5, 2g → 1.875, 4g → 1.4878)
     - SI hits :math:`k_\infty` to :math:`< 10^{-8}`
   * - Heterogeneous 2G fuel\|moderator, **non-flat** flux
       (:func:`test_si_krylov_heterogeneous_2g_nonflat_flux <tests.gates.sn.eigenvalue.test_keff_2d.TestSIKrylov2DEquivalence.test_si_krylov_heterogeneous_2g_nonflat_flux>`)
     - SI vs Krylov flux **shape** + eigenvalue
     - flux shape agrees to :math:`\sim 10^{-9}`
   * - 2-D SI :math:`k_{\rm eff}` Cauchy under refinement
       (:func:`test_si_2d_keff_converges_under_refinement <tests.gates.sn.eigenvalue.test_keff_2d.TestSIKrylov2DEquivalence.test_si_2d_keff_converges_under_refinement>`)
     - Self-convergence (consistency regression catcher)
     - monotone, single fixed point

The structural-independence discipline (``vv-principles`` L11 /
anti-pattern #1) applies: SI ≡ Krylov *alone* is twin-path agreement —
necessary but **not** sufficient, since both could share a defect. It
becomes correctness evidence only because the homogeneous leg
independently anchors the same production path to the **closed-form**
:math:`k_\infty` eigenvalue (the closed-form pillar; per ``vv-principles``,
twin-implementation agreement is L4-class on its own and MMS does not
prove eigenvalues). The heterogeneous leg carries a genuinely non-flat
(≥2G, fuel\|moderator) flux so the angular / wavefront redistribution
terms are active rather than nulled (``vv-principles`` anti-patterns
#3 / #4), and the un-xfailed L2 mesh-convergence pin
(``tests/gates/sn/sweep/cartesian_2d/test_discrete_ordinates_2d.py::test_do_mesh_convergence``,
the ERR-003 catcher) plus the ``2d_2g_LS4_dd_8x4_het_si`` regression
snapshot round out the catch surface.


.. _bc-extraction-numerical-evidence:

Numerical evidence
------------------

The extraction is verified by three independent grounds (per the
``vv-principles`` skill's three pillars and the bit-identity
vs principled-equivalence gate).

**1. Vacuum bit-identity.** With :math:`B = 0` the boundary gain
contributes nothing (:math:`B.\text{apply} \equiv 0`), so the variadic
matvec :math:`(L+C).\text{apply} - S.\text{apply} - B.\text{apply}`
reduces exactly to :math:`(L+C).\text{apply} - S.\text{apply}` and the
vacuum path is **bit-identical** to the pre-extraction matvec.  (The
two-gain spelling here is Wave O's, which is what the cited captures
were taken against; the :math:`(n,2n)` gain :math:`N_{2n}` joined the
tuple at CS4c step 3 and rides this argument unchanged, contributing
exactly zero on the :math:`\Sigma_{2n} \equiv 0` fixtures the captures
use.) Verified
by:

- the matvec 18-baseline snapshot
  (:func:`numpy.array_equal` against the pre-O.4a.2 captures across
  slab / sphere / cylinder × 1G / 2G / asymmetric :math:`\Sigma_s` ×
  vacuum / white / specular), and
- the end-to-end regression snapshots.

This is the bit-identity-by-inheritance gate: vacuum keeps the
verified pre-extraction value for free (``vv-principles``
§"Bit-identity vs principled-equivalence", criterion: implementation
unchanged on the vacuum path because the bare sweep reads a zero
inflow seed).

**2. Reflective convergence-equivalence (closed-form pillar).** The
reflective path relocates the reflection from inside the sweep to the
sibling :math:`-B`, so it is **not** bit-identical but
*convergence-equivalent* to a structurally-independent analytical
reference:

.. list-table:: Reflective convergence-equivalence gates
   :header-rows: 1
   :widths: 40 30 30

   * - Test
     - Reference (pillar)
     - Both solvers?
   * - Curvilinear streaming-equilibrium
       (``tests/gates/sn/sweep/curvilinear/test_streaming_equilibrium_curvilinear.py``)
     - Analytical infinite-medium balance
       :math:`\phi = q/\Sigma_t` (closed-form)
     - ``source_iteration`` AND ``krylov``
   * - Reflective :math:`k_\infty` homogeneous
       (``tests/gates/sn/verification/analytical/test_kinf_homogeneous.py``)
     - :math:`k_\infty = \nu\Sigma_f / \Sigma_a` (closed-form
       eigenvalue — MMS does NOT prove eigenvalues)
     - both
   * - ``test_si_carve_recovers_analytical_kinf``
       (``tests/gates/sn/operators/test_invertible_operator.py``)
     - Analytical :math:`k_\infty` via the SI path with :math:`B`
       folded (closed-form)
     - SI path
   * - Invertible-operator :math:`Q/\Sigma_t` recovery
       (``tests/gates/sn/operators/test_invertible_operator.py``)
     - Flat-flux fixed-source balance (closed-form)
     - direct ``−B`` drive

**3. Reflective eigenvalue regression (principled-equivalence ULP).**
The reflective cylinder eigenvalue regression snapshot **drifts within
tolerance**: :math:`4\times 10^{-13}` on :math:`k_{\rm eff}` and
:math:`7\times 10^{-12}` relative on the scalar flux. This is **not** a
bug — it is FP-non-associativity from relocating the reflection
(``vv-principles`` § criterion 3: the reduction-tree changes because
the reflection now happens in :math:`-B` rather than fused into the
sweep, so additions occur in a different IEEE-754 order). The drift is
bounded by ``iteration_count × condition_number × ULP``, well under the
existing ``rtol`` regression tolerance. The new value is
convergence-equivalent to the analytical references above (criterion
2), so the regression contract is satisfied without relaxation beyond
the snapshot tolerance.

**The O.2a variadic reassociation is a second principled-equivalence
instance.** Splitting the matvec from :math:`(L+C) - (S+B)` (the
retired fold) to :math:`(L+C) - S - B` (the two separate gains
:ref:`bc-extraction-variadic-driver` had at the time — :math:`N_{2n}`
became a third at CS4c step 3), and the rhs symmetrically,
re-associates the same additions into a different IEEE-754 order. The
regression snapshots drift at FP-noise level — reflective cylinder
:math:`4.2\times 10^{-13}` on :math:`k_{\rm eff}` and :math:`6.8\times
10^{-12}` relative on the scalar flux, anisotropic 3–5 ULP — all within
the existing tolerances (:math:`10^{-11}` / :math:`10^{-9}` /
:math:`10^{-12}`). Per ``vv-principles`` criterion 2 the new value is
verified against **structurally-independent** references, not merely
shown close to the old value: the NEW-1 closed-form :math:`Q/\Sigma_t`
flat-flux balance, the SI ≡ Krylov twin (:ref:`bc-extraction-2d-si-krylov-twin`),
and the ``keff_2d`` closed-form :math:`k_\infty`. The reassociation
satisfies all three criteria (named intermediates — each gain's output
is a principled source/sink; structurally-independent reference;
dimensionally-explainable drift), so no contract relaxation is needed.


.. _bc-extraction-operator-output-typing:

Operator-output role typing — :math:`A\psi` is a source/sink
------------------------------------------------------------

Wave O step B.5.2 (`Issue #208
<https://github.com/deOliveira-R/ORPHEUS/issues/208>`_, commit
``6ef5063``, 2026-06-03) retyped every SN operator's ``.apply`` output
``.boundary`` from
:class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux` to
:class:`~orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink`
— the *source/sink* role leaf. This completes the **boundary** half of
the B.5 operator-output "dimensional-sin" carve; the **bulk** half
(``.apply.bulk`` →
:class:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink`)
landed earlier in commit ``f400743``. The two halves make the boundary
role grid a clean parallel of the bulk.

.. list-table:: The completed role grid (bulk ‖ boundary)
   :header-rows: 1
   :widths: 18 28 28 26

   * - Block
     - ``.apply`` (operator output :math:`A\psi`)
     - ``.solve`` (swept solution trace)
     - ``from_balance`` (the defect)
   * - **bulk** (:math:`V_{\rm bulk}`)
     - :class:`~orpheus.transport.fields.angular_source_sink.AngularSourceSink`
       (``f400743``)
     - :class:`~orpheus.transport.fields.angular_flux.AngularFlux`
     - :class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
       (consumed by :func:`~orpheus.sn.solver.evaluate_residual`,
       O.2 — :ref:`affine-typed-residual`)
   * - **boundary** (:math:`V_{\rm inflow} \oplus V_{\rm outflow}`)
     - :class:`~orpheus.transport.fields.angular_boundary_source_sink.AngularBoundarySourceSink`
       (``6ef5063``, B.5.2)
     - :class:`~orpheus.transport.fields.angular_boundary_flux.AngularBoundaryFlux`
     - :class:`~orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual`
       (consumed by :func:`~orpheus.sn.solver.evaluate_residual`,
       O.2 — :ref:`affine-typed-residual`)

The governing principle (the load-bearing rationale)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   *A residual only arises after we compare an operator output against
   something else and get a defect (a balance). The output of an
   operator is NOT a residual straightaway.*

Each operator's ``.apply`` emits :math:`A\psi` — a **source/sink**
(a signed reaction-rate / flux density: a *source* when produced, a
*sink* when it is an operator-loss output such as :math:`L\psi`; the
single role leaf holds both, hence ``SourceSink``). The residual is
**only** the named composition
:meth:`AngularBoundaryResidual.from_balance(lhs, rhs) <orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual.from_balance>`
of the affine boundary balance
:math:`r_\Gamma = \gamma_-\psi - (R\,G\,\gamma_+\psi + q)`. The GMRES
*flat* residual :math:`b - A\psi` is formed internally on the **raveled
vector** (via :meth:`TimedFullField.to_flat <orpheus.transport.timed_full_field.TimedFullField.to_flat>`)
and is **never typed as a field** — so at B.5.2
:class:`AngularBoundaryResidual` had no operator-output consumer; its first
consumer is the honest :math:`L+C-S-N_{2n}-F-B` driver of Wave O step
**O.2** (see :ref:`bc-extraction-operator-output-o2`). That consumer
has since landed:
:func:`~orpheus.sn.solver.evaluate_residual` types the balance defect
:math:`(L+C-S-N_{2n}-B)\psi - q` via ``from_balance`` (Wave O step O.2
close-out, :ref:`affine-typed-residual`), so
:class:`AngularBoundaryResidual` and its bulk sibling
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
are now consumed, not merely minted.


The two-hat tension and why ``AngularBoundarySourceSink`` dissolves it
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The earlier in-flight plan
(``.claude/plans/b52_boundary_residual_retype.md``) proposed typing the
matvec output boundary as :class:`AngularBoundaryResidual`. That choice was
**rejected** for two reasons:

1. **It breaks consistency with the already-landed bulk.** The bulk
   ``.apply.bulk`` uses the source/sink leaf
   (:class:`AngularSourceSink`), **not** a residual, for operator
   outputs. Typing the boundary output as a residual would make the two
   halves of the same carve disagree on what an ``.apply`` output *is*.

2. **It creates a "two-hat" tension that the class gate cannot
   satisfy.** The realized boundary operator
   :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` (:math:`B`)
   emits :math:`B\,\psi.\text{outflow}`, and that **same** emission is
   consumed two ways:

   .. list-table:: :math:`B`'s two consumers — the "two hats"
      :header-rows: 1
      :widths: 24 38 38

      * - Consumer
        - Composition
        - The hat :math:`B\,\psi.\text{outflow}` would wear
      * - Krylov matvec
        - :math:`(L+C).\text{apply} - S.\text{apply} -
          N_{2n}.\text{apply} - B.\text{apply}`
        - a **residual** term (subtracted from the diagonal)
      * - SI rhs
        - :math:`q_{\rm ext} + S.\text{apply} + N_{2n}.\text{apply} +
          B.\text{apply}`
        - a **source** term (the inflow seed the bare sweep reads)

   One operator cannot emit :class:`AngularBoundaryResidual` for the matvec
   **and** :class:`AngularBoundarySourceSink` for the SI rhs — the
   :class:`TimedFullField <orpheus.transport.timed_full_field.TimedFullField>`
   class gate (strict class identity:
   ``type(self.boundary) is not type(other.boundary)`` ⟹ ``TypeError``)
   throws on ``AngularBoundaryResidual + AngularBoundarySourceSink`` the moment the SI
   rhs tries to add :math:`B.\text{apply}` (a residual, under OPT-BR)
   to :math:`S.\text{apply}` and :math:`q_{\rm ext}` (sources). The
   variadic driver (:ref:`bc-extraction-variadic-driver`) makes this
   sharper than the retired fold: each gain's output is summed
   *individually*, so :math:`B`'s lone hat must be a source/sink for the
   rhs sum :eq:`bc-extraction-variadic-matvec` to close.

Choosing :class:`AngularBoundarySourceSink` for **all** operator outputs
dissolves the two-hat: :math:`B` wears **one** hat (it always emits a
source/sink), and **both** sums close as homogeneous
:class:`AngularBoundarySourceSink` sums —

.. math::
   :label: bc-extraction-two-hat-closed-sums

   \underbrace{(L+C).\text{apply} - S.\text{apply}
               - B.\text{apply}}_{\text{Krylov matvec}}
   \quad\text{and}\quad
   \underbrace{q_{\rm ext} + S.\text{apply}
               + B.\text{apply}}_{\text{SI rhs}}

both stay within the single :class:`AngularBoundarySourceSink` class. This
needs **no SI-driver restructure** and **no partial-O.2**:
:class:`AngularBoundaryResidual` stays reserved for the named
``from_balance`` composition exactly as
:class:`AngularResidual` waits on the bulk.

.. vv-status: bc-extraction-two-hat-closed-sums documented

A throwaway **decision instrument**
(``derivations/diagnostics/diag_b52_boundary_typing_decision.py``, the
B0 de-risk) proved on a 1-D reflective slab **and** a 2-D reflective
box that the OPT-BSS choice (``AngularBoundarySourceSink`` for the matvec
output) closes both sums, while the OPT-BR choice
(:class:`AngularBoundaryResidual` for the matvec output) throws the two-hat
``TypeError`` on the SI rhs.


Why the Krylov path is safe with a ``AngularBoundarySourceSink`` matvec output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The matvec output never escapes scipy as a :class:`AngularBoundarySourceSink`,
so the *solution* side stays :class:`AngularBoundaryFlux`. The mechanism is the
flat round-trip:

* :meth:`TimedFullField.to_flat <orpheus.transport.timed_full_field.TimedFullField.to_flat>`
  ravels the composite to ``[bulk.values.ravel(), boundary.values]`` —
  a **type-agnostic** 1-D vector (the class of ``.boundary`` is erased).
* scipy's GMRES iterate is reconstructed via
  :meth:`TimedFullField.from_flat <orpheus.transport.timed_full_field.TimedFullField.from_flat>`,
  which rebuilds the boundary with
  ``replace(template.boundary, values=...)`` off the flux
  ``solution_template``. Because the template's boundary is a
  :class:`AngularBoundaryFlux`, the reconstructed iterate's boundary is a
  :class:`AngularBoundaryFlux`.

So the matvec's *internal* :class:`AngularBoundarySourceSink` boundary class
lives only inside one ``op.apply`` call; the moment the result is
raveled and handed to scipy, the class is gone, and the iterate scipy
hands back is reconstructed as the flux type. The solve/iterate/trace
sites are therefore correctly **kept** :class:`AngularBoundaryFlux`:
:meth:`MultiplicationOperator.solve <orpheus.transport.operators.multiplication_operator.MultiplicationOperator.solve>`
(the collision multiplier :math:`C = M[\sigma_t]`),
the boundary buffer of
:meth:`StreamingCollisionOperator._solve_timed_full_field <orpheus.sn.operators.streaming.StreamingCollisionOperator._solve_timed_full_field>`,
the cold-start ``initial_guess`` iterates
(``TimedFullField.zeros(..., boundary=AngularBoundaryFlux, ...)``), the
converged traces, and the sweep's persistent boundary buffer.


The 13 retyped sites
~~~~~~~~~~~~~~~~~~~~

Thirteen sites (operator outputs + ``q_ext`` sources) flipped from
:class:`AngularBoundaryFlux` to :class:`AngularBoundarySourceSink`:

.. list-table:: B.5.2 retyped sites
   :header-rows: 1
   :widths: 34 38 28

   * - Module / symbol
     - Site
     - Emission
   * - :mod:`orpheus.sn.loss_representation`
     - ``_OneDimScanWalk._apply_walk`` (``m_boundary``) — the fused
       1-D matvec body; #206 relocated it here, #238 folded the former
       ``_compute_LpC`` / ``_compute_decomposition`` /
       ``_SpatialSweepDirection`` sites into this single walk
     - :math:`L+C` boundary block
   * - :mod:`orpheus.sn.operators.streaming`
     - :meth:`StreamingOperator._apply_2d_cartesian <orpheus.sn.operators.streaming.StreamingOperator>`
     - 2-D boundary block
   * - :mod:`orpheus.transport.operators.multiplication_operator`
     - :meth:`MultiplicationOperator.apply <orpheus.transport.operators.multiplication_operator.MultiplicationOperator.apply>`
       (the collision multiplier :math:`C = M[\sigma_t]`)
     - bulk → bulk; boundary zero
   * - :mod:`orpheus.transport.operators.scattering`
     - :meth:`ScatteringOperator.apply <orpheus.transport.operators.scattering.ScatteringOperator>`
     - boundary zero
   * - :mod:`orpheus.transport.operators.n2n`
     - :meth:`N2NOperator.apply <orpheus.transport.operators.n2n.N2NOperator>`
     - boundary zero
   * - :mod:`orpheus.transport.operators.fission`
     - :meth:`FissionOperator.apply <orpheus.transport.operators.fission.FissionOperator>`
     - boundary zero
   * - :mod:`orpheus.sn.operators.boundary`
     - :meth:`SNBoundaryOperator._apply_faces <orpheus.sn.operators.boundary.SNBoundaryOperator>`
       (``apply`` **and** ``apply_transpose``)
     - :math:`B\,\psi.\text{outflow}` on the inflow slots
   * - :mod:`orpheus.sn.solver`
     - ``q_ext.boundary`` at the **3 source builds** (eigenvalue SI,
       eigenvalue Krylov, fixed-source SI / reconstruction)
     - prescribed inflow (zero for vacuum / reflective)

(The B.5.2 ``_zero_within_group_fission`` slot — a ``ZeroOperator``
``codomain_zero`` emitting the boundary zero for an explicit ``F = 0``
within-group operator — was designed here but NEVER wired: the Wave-O
within-group decomposition (:func:`~orpheus.sn.coupled_system.build_within_group_system`)
routes within-group fission through ``q_ext`` instead, so no zero-fission
operator is ever constructed. The dead helper retired 2026-07-03 (C4).)

The change is **type-only**: the zero-trace allocation (then spelled
``AngularBoundarySourceSink.zeros_on(mesh)``, since CS4b S5
``.zeros(mesh.angular_trace)``)
and the per-face-view writes produce **bit-identical** ``.values`` —
only the wrapping role-type differs. The dead :class:`AngularBoundaryFlux`
runtime imports were retired from the retyped sites.


.. _bc-extraction-operator-output-o2:

The extraction close-out — where the remaining items landed
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Wave O step **O.2a has landed the honest** :math:`L+C-S-B` **driver**
via the variadic couplings of :ref:`bc-extraction-variadic-driver`:
the transitional :math:`S + B` fold is **retired** and :math:`B` is now
a first-class coupling gain. Of the items B.5.2 left for the rest of
O.2, **the adjoint metric and its gate landed in step O.2b
R5** (:ref:`g-adjoint`), and **the residual column landed
in the O.2 close-out** (:ref:`affine-typed-residual`):

* :meth:`AngularResidual.from_balance <orpheus.transport.residuals.angular_residual.AngularResidual.from_balance>`
  and
  :meth:`AngularBoundaryResidual.from_balance <orpheus.transport.residuals.angular_boundary_residual.AngularBoundaryResidual.from_balance>`
  now have their first production-reachable consumer:
  :func:`~orpheus.sn.solver.evaluate_residual` types the within-group
  balance defect :math:`r = (L+C-S-N_{2n}-B)\psi - q` as the composite
  ``FullField(bulk=AngularResidual, boundary=AngularBoundaryResidual)`` — the
  **timeless** carrier, because a residual is a one-shot balance defect
  carrying no iteration history (the ``history_depth = 0`` degenerate;
  P4.5 W-C confines the timed type to the driver iterate)
  (see :ref:`affine-typed-residual`). The honest variadic driver still
  emits each gain's output as a :class:`AngularBoundarySourceSink` and the
  GMRES defect is still the *flat* :math:`b - A\psi` on the raveled
  vector — the typed residual is an **additive diagnostic + DSA
  substrate**, never in the convergence path (so the converged flux
  stays bit-identical).
* the :math:`|\Omega\cdot\hat n|\,w` :math:`G`-metric adjoint ``.H``
  (the boundary-weighted inner product for the transpose) **landed in
  R5** — ``op.H`` is now the metric-correct G-adjoint :math:`G^{-1}
  A^{\mathsf T} G` over the composite :math:`V_{\rm bulk}\oplus
  V_{\rm trace}` (:ref:`g-adjoint`). This is exactly why
  :math:`B` stays trace-typed as a separate gain
  (:ref:`bc-extraction-variadic-driver`): a bulk-folded :math:`B` would
  erase the trace metric the adjoint needs.
* **Gate-1.3** (the O.2 adjoint verification gate) **landed in R5** —
  the dense-probe oracle + the L11 wrong-metric control
  (:ref:`g-adjoint`).

The direct-helper survived O.2a (the driver stopped routing through it,
but the final eigenvalue reconstruction sweep still did) and has **no
production caller at all since #448**, when that sweep became one step of
the driven map (:ref:`sn-finalize-one-step`); it now lives with the
sweep-tier gates (:ref:`bc-extraction-two-routes`).  The
:attr:`ScatteringOperator.domain` typing completion that was once a
documented seam **landed in P4.5 W-D** — :math:`S` now carries the
composite full-field space (:ref:`bc-extraction-scope-future`).

The residual column is now wired: :class:`AngularBoundaryResidual` and
:class:`~orpheus.transport.residuals.angular_residual.AngularResidual`
are consumed by :func:`~orpheus.sn.solver.evaluate_residual` (Wave O
step O.2 close-out — :ref:`affine-typed-residual`). ⛔ That close-out
also completed the *affine* flux algebra — the iterate increment was a
typed ``AngularDisplacement`` and ``flux + flux`` was a
:class:`TypeError` — and **that half was overturned on 2026-08-19**:
flux lives in the positive cone of a vector space, so ``flux + flux`` is
legal and the increment is the flux type carrying a signed value
(:ref:`cone-typed-field-algebra`). The residual half above is unaffected
— the residual role was always a plain vector role. The canonical home
for the field algebra is now :ref:`cone-typed-field-algebra`.


.. vv-status note: the operator-output role typing of B.5.2 is a
   type-only refactor (bit-identical ``.values``); its correctness is
   verified by the same gates that verify the extraction
   (:ref:`bc-extraction-numerical-evidence`) plus the type-residual
   gates catalogued below. B.5.2's verification ground:

* the **B0 decision instrument** (``diag_b52_boundary_typing_decision.py``,
  a throwaway probe that was not preserved in git), which proved OPT-BSS
  closes both sums while OPT-BR throws the two-hat;
* the core operator / boundary / 2-D suite (324 passed);
* SI eigenvalue slab / sphere / cylinder × 1 / 2 / 4-group — the
  **two-hat exerciser** (the SI rhs sum that OPT-BR would throw on);
* Krylov :math:`k_\infty` (14 cases);
* the type-residual gates (``test_native_matvec`` boundary-output type
  assert migrated; positive type asserts added for the 2-D matvec
  ``test_bc_extraction_2d`` and for :math:`B` in
  ``test_sn_boundary_operator``);
* the dimensional-check sentinel suite (36 / 36, run without ``-O``);
* MMS L1 1-D + 2-D + curvilinear (8 passed, 6 xfail — flux-shape /
  convergence-order pillar; MMS does **not** prove the eigenvalue).

The change was reviewed by the ``elegance-enforcer`` (PASS, no
conditions).


.. _bc-trace-structure:

Trace structure (Γ\_-, Γ\_+)
============================

The transport equation lives on a phase space :math:`\Omega \times
S^d` where :math:`\Omega \subset \mathbb{R}^d` is the spatial domain
and :math:`S^d` is the unit sphere of directions. The boundary
:math:`\partial\Omega` carries an outward unit normal :math:`\hat
n(\mathbf{r})` at every regular point. For an angular flux
:math:`\psi(\mathbf{r}, \Omega)` defined on the full phase space, the
**boundary trace** splits naturally into two pieces by the sign of
:math:`\Omega \cdot \hat n`, :math:`\Gamma_-` (inward) and
:math:`\Gamma_+` (outward). The equation is stated once, as
:eq:`trace-half-decomposition` on the operator-algebra page, where the
two selectors that implement it are bound; this page carries its
discrete realisation, :eq:`inflow-mask-discrete`.

Points with :math:`\Omega \cdot \hat n = 0` are **tangential** —
they belong to neither half. For axis-aligned ordinates on
axis-aligned faces these arise exactly (no round-off) for face
normals perpendicular to the ordinate's direction cosine; for
general curvilinear faces or generic ordinates they arise only
at a measure-zero subset that the discrete representation
identifies via a small tolerance
(``TANGENTIAL_EPS = 4 * np.finfo(np.float64).eps``, i.e.
:math:`\approx 8.9\times 10^{-16}`, defined in ``orpheus.numerics.tangential``
and re-exported by :mod:`orpheus.numerics.spaces.angular_trace_space`). That constant is
the codebase's ONE outflow/inflow classifier: campaign phase **B3.4a**
retired the white operator's private ``> 0.0`` twin precisely because a
strict compare disagrees with it wherever a quadrature carries
tangential ordinates (:ref:`bc-narrowing-b34a`).

In the discrete setting, the spatial boundary is a union of finite
faces :math:`\{f_1, \ldots, f_F\}` and the angular variable is a
finite ordinate set :math:`\{\Omega_n : n = 1, \ldots, N\}`. The
sign predicate :eq:`trace-half-decomposition` then collapses to a
**per-face boolean mask** of shape :math:`(F, N)`:

.. math::
   :label: inflow-mask-discrete

   \mathrm{inflow\_mask}[f, n]
   \;=\; \bigl(\Omega_n \cdot \hat n_f < -\epsilon\bigr),
   \qquad
   \mathrm{outflow\_mask}[f, n]
   \;=\; \bigl(\Omega_n \cdot \hat n_f > +\epsilon\bigr).

.. _bc-half-trace:

**A half-trace is the restriction of the whole-boundary trace space to
one face's inflow or outflow selector**: :math:`\Gamma_-(f)`, the rows of
face :math:`f`'s ordinate slot with :math:`\Omega_n\cdot\hat n_f <
-\epsilon`, or :math:`\Gamma_+(f)`, the rows with :math:`\Omega_n\cdot\hat
n_f > +\epsilon`, each carrying the trace metric
:math:`|\Omega_n\cdot\hat n_f|\,w_n` restricted to its own rows. It is a
space in its own right,
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`,
minted and cached by the whole-boundary space
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
(:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_space`,
:meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_space`),
and this paragraph is the corpus's one definition of the term. A
realized boundary operator is an arrow from an outflow half-trace to
the inflow half-trace of the face it fills,
:math:`\Gamma_+(f') \to \Gamma_-(f)`, with :math:`f' = f` for every wall
condition and for a mirror deck and :math:`f'` the opposite face for a
translation wrap (the realized type's ``domain_face``, :ref:`bc-worked-example`
Step 5). Every unit the S\ :sub:`N` realizer table realizes is realised
narrowed in this sense. On a Cartesian face with a mirror-closed quadrature the
one declaration it refuses is a method-only condition of another method,
``ZeroFluxCondition``, refused on the state-cone axis (Step 4 of
:ref:`bc-worked-example` runs the six spellings and prints the refusal); the
radial sites refuse more (the diffuse chain and every source at a radial outer face,
#653; everything but the unit specular cavity at a hollow body's inner face,
#511), and a mirror on a quadrature that ``ordinate_pairing`` cannot certify is
refused on the angular-resolution axis. Two half-traces do not make up the face: the tangential rows belong
to neither, and the count of them depends on the quadrature:

.. code-block:: python

   import numpy as np
   from orpheus.numerics.quadrature import Quadrature
   from orpheus.numerics.tangential import TANGENTIAL_EPS

   for name, rule in [("gauss_legendre(4)", Quadrature.gauss_legendre(4)),
                      ("product(4,4)", Quadrature.product(4, 4)),
                      ("lebedev(17)", Quadrature.lebedev(17))]:
       tangential = int(np.sum(np.abs(rule.axis_cosines(0)) <= TANGENTIAL_EPS))
       print(f"{name:18s} {tangential:3d} of {rule.n_ordinates:3d} ordinates are tangential on an x face")

.. code-block:: text

   gauss_legendre(4)    0 of   4 ordinates are tangential on an x face
   product(4,4)         8 of  16 ordinates are tangential on an x face
   lebedev(17)         12 of 110 ordinates are tangential on an x face

The run above is the command behind the counts, `[M]` 2026-10-24: 0 of 4, 8 of 16 and 12 of
110 ordinates on an :math:`x` face of ``gauss_legendre(4)``, ``product(4,4)`` and ``lebedev(17)``. Those
rows are, on every shipped quadrature, exactly the kernel of the face's metric (the
measurement is recorded on the class), so the whole face tier is the direct
sum of its two half-traces only in the quotient by that kernel, and a
half-trace, whose metric is strictly positive, can be the intermediate
space of a factored boundary response where the whole face tier cannot.
The face and the direction are part of a half-trace's identity (its name
is ``angular_trace[<face>:<role>]``), because the two faces of a slab
carry half-traces of equal size over different ordinates.

The mask :eq:`inflow-mask-discrete` is the discrete realization of
:math:`\Gamma_\pm`. The trace space stores the signed projection
:math:`\Omega_n\cdot\hat n_f` and derives the two masks from it on demand
(Step 3 of :ref:`bc-worked-example` prints it). It is the load-bearing primitive that downstream
consumers need:

* **The SN realized operators read BOTH half-traces off the trace ladder**:
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_space`
  and
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_space`
  give a face's :math:`\Gamma_+` and :math:`\Gamma_-` (the index sets
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`
  and
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_indices_for_face`
  are what those spaces are built from), and the sweep composes a face's
  action through the restrictions
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_restriction`
  and
  :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.inflow_restriction`
  (:ref:`bc-worked-example`, Steps 3 and 6). The method space the realizer
  receives, :class:`~orpheus.sn.mesh.method_space.SNMethodSpace`, holds
  exactly two fields, the quadrature and the trace. It held the index sets
  as ``inflow_indices`` / ``outflow_indices`` **fields** until WF2.4 (so that a
  hand-built space could name its own half-traces, :ref:`bc-domain-narrowing`),
  together with an ERR-041 cross-check of the claimed ``inflow_indices`` against the
  signed projection the face name implies. Both retired: the typed operators
  read their half-spaces off the ladder, so ``outflow_indices`` had no production
  reader and ``inflow_indices`` was read only by the check, which compared an
  index set the zero map never consumed. ERR-041's invariant, that the zero map fills
  only the incoming trace, is carried by the ladder itself:
  ``tests/gates/numerics/test_angular_trace_space.py::TestTheInflowHalfSpaceIsTheIncomingOrdinates``
  witnesses that ``trace.inflow_space(face)`` holds exactly the ordinates with
  :math:`\Omega\cdot\hat n_f < 0`.
* **A source is placed on** :math:`\Gamma_-` **by typing and checked at the
  one place it is sampled.** ``SNSource.sample``
  (in ``orpheus.sn.boundary.realizer``) raises
  :class:`~orpheus.sn.boundary.errors.BoundarySourceNotOnIncomingTraceError`
  (ERR-047) when the space it samples on is not an inflow half-trace. The
  realized linear part never carries a source (a pure inflow's linear image
  is the zero map; :math:`q` is an array on :math:`\Gamma_-`), so
  :math:`q \in \Gamma_-` holds by construction and the check guards the
  one caller that could pass another space. The retired form was the universal invariant
  ``BoundaryTraceLaw.assert_source_is_placeable``, which read the inflow set as a
  **structural** check (can the realization NAME :math:`\Gamma_-`?). ⭐ Until
  campaign phase **P6** it probed the source's *support* on the per-ordinate
  shape first, and a source that answered that probe with zeros **skipped the
  certification entirely** while still delivering. Since a spec received
  :math:`\Gamma_-(f)` itself from P6 on, :math:`q \in \Gamma_-` held by
  construction, the probe was retired, and the discriminator was whether the
  law carried a source at all (``NoSource`` or not) rather than what that
  source currently evaluated to; since **B3.4a** the realizer sized the
  source's block from :math:`|\Gamma_-|` (before B3.4a it held because the
  realized operator masked a full-face evaluation, :ref:`bc-narrowing-b34a`).
* The SN curvilinear sweep (1-D spherical / cylindrical) consumes
  the same realizer-routed trace as the Cartesian path — Issue #188
  (C188.1+C188.2 in :mod:`orpheus.numerics.spaces.angular_trace_space`, C188.3 in
  :mod:`orpheus.sn.problem`) lifted the curvilinear deferral and
  Issue #176 then dropped the legacy 2-arg shim that existed only
  to bridge that deferral.

The class :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
carries the signed projection :math:`\Omega\cdot\hat n_f` as the ``omega_dot_n``
array of shape ``(n_faces, N)``, beside the ordinate ``directions`` and the
``layout``; the three are required keyword-only fields excluded from
:meth:`__eq__` and :meth:`__hash__` (``compare=False``) — preserving the
:class:`~orpheus.numerics.space.FunctionSpace` identity convention
``(name, shape)``, which for this axes-less class IS content identity
(its factory folds the quadrature and layout into an
``angular_trace#<digest>`` name) and which the 2026-09-07 identity flip
therefore leaves unchanged: the flip made identity structural only for
spaces that declare their ``axes``. Construction goes through the classmethod factory
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`;
since the three fields are required, a trace space with an unpopulated
projection is not constructible.

**Geometry-blind, layout-driven (post Issue #225 / C5.3).** The factory
is **geometry-blind**: it takes only the angular quadrature and a
:class:`~orpheus.numerics.face_layout.FaceLayout` (canonically
:attr:`SNProblem.boundary_face_layout <orpheus.sn.problem.SNProblem>`), and
reads every datum from those two — the layout's ``"{axis}{min|max}"``
face names imply axis-aligned outward normals, so the
:math:`\Omega\cdot\hat n` row for an axis-:math:`a` face is
:math:`\pm\mu_a` (sign from ``min`` / ``max``). It works for every
constructible :class:`~orpheus.mesh.structured.Mesh1D` coord system
(``CARTESIAN`` / ``SPHERICAL`` / ``CYLINDRICAL`` — all share the
``("xmin", "xmax")`` radial-axis face structure, with the
:math:`\mu_x` of a :meth:`Quadrature.gauss_legendre
<orpheus.numerics.quadrature.Quadrature.gauss_legendre>` quadrature as
the direction cosine along that axis), for 2-D Cartesian, and — since
C5.5 — for 3-axis Cartesian. The factory has **no** geometry refusal:
the former ``mesh`` parameter (on the retired
``from_mesh_and_quadrature``) was gate-only and is gone (see
:ref:`sn-c5-geometry-blind-trace`). The 2-D cylindrical
(axisymmetric :math:`(r, z)`) case never reaches the factory because
such a :class:`~orpheus.mesh.structured.Mesh2D` cannot become an
:class:`SNProblem` (no 2-D cylindrical SN sweep exists); the refusal lives
at the :class:`SNProblem` construction surface, not in the trace factory.


.. _bc-law-layer:

Boundary law (``BoundaryTraceLaw`` ABC + concretes)
===================================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Until WF2.4 a boundary was a typed LAW. The base class
``BoundaryTraceLaw``
was an ``abc.ABC`` whose MRO was exactly ``[BoundaryTraceLaw,
RegistryMixin, ABC, object]``: it mixed in
:class:`~orpheus.numerics.registry.RegistryMixin` (so each concrete
subclass self-registered under its ``key=`` class-creation kwarg) and
nothing else. It did **not** inherit
:class:`~orpheus.numerics.operator.LinearOperator` — that inheritance
was dropped at Issue #186 / B3 + β2 (see
:ref:`bc-trace-law-descriptor-model`), and the class carried no ``@`` /
``__matmul__``. The dunders it did carry (``+``, ``-``,
``*``, ``/``, unary ``-``) were the **descriptor-tree** algebra: they
returned ``LawSum`` / ``LawScaled`` nodes, never operators.
The ABC shipped:

1. Three properties named for the :eq:`affine-bc-form` factors:
   ``geometry_map``, ``response_kernel``, ``source``. The ABC's
   defaults were ``None``, ``None``,
   ``NoSource``, but **all seven
   concrete laws populated all three** since campaign phase B1 — the
   two factor tiers were typed specifications
   (``SelfPairedDeck`` /
   ``PairedDeck`` for :math:`G`;
   ``ScalarResponse`` /
   ``LambertianReemission`` for
   :math:`R`), never realized matrices. As measured at the last tree
   that carried them:

   .. list-table::
      :header-rows: 1
      :widths: 30 34 36

      * - law
        - ``geometry_map``
        - ``response_kernel``
      * - ``VacuumInflow``
        - ``SelfPairedDeck.identity()``
        - ``ScalarResponse(alpha=0.0)``
      * - ``ReflectiveBoundary``
        - ``SelfPairedDeck.mirror(axis)``
        - ``ScalarResponse(alpha=1.0)``
      * - ``WhiteBoundary``
        - ``SelfPairedDeck.identity()``
        - ``LambertianReemission(alpha, axis, outward_sign)``
      * - ``AlbedoBoundary``
        - ``SelfPairedDeck.identity()``
        - by closure (**B3.4b**): ``SpecularReemission(alpha, axis)`` /
          ``LambertianReemission(alpha, axis, outward_sign)`` /
          ``ScalarResponse(alpha=albedo)`` with none
      * - ``PeriodicBoundary``
        - ``PairedDeck.wrap(axis)``
        - ``ScalarResponse(alpha=1.0)``
      * - ``ZeroFluxBoundary``
        - ``SelfPairedDeck.identity()``
        - ``ScalarResponse(alpha=-1.0)``
      * - ``PrescribedInflow``
        - ``SelfPairedDeck.identity()``
        - ``ScalarResponse(alpha=0.0)`` (plus ``source``)

   They were **read by production**: phase B2.2 repointed five sites
   from string comparison to
   ``law_permutes_ordinates(law)`` (a function over BOTH tiers since
   **B3.4b**, when a specular pairing became legal in :math:`R`) /
   ``law.response_kernel.is_zero`` /
   ``law.response_kernel.amplitude``, and the diffusion realizer's
   whole law → :math:`\mathcal{A}` table collapsed to the single read
   ``law.response_kernel.amplitude``. The ABC kept the ``None``
   defaults, so the diffusion realizer guarded against a subclass that
   never populated them rather than dying on an ``AttributeError``.
2. Five **universal** ``assert_*`` invariants and three **specific**
   invariants on the BCs that needed them — but four of the five
   universals were empty and two of those were overridden by nobody.
   :ref:`bc-universal-invariants` holds the inventory.
3. A ``realize`` hook that raised with guidance (route through a
   method realizer) — see :ref:`bc-realizer-layer-detail`.
4. **No ``apply`` method at all** (Issue #186 / B3 + β2,
   2026-05-11). The descriptor model that survived the C176.3
   Option A interim was retired in favour of a pure-descriptor
   contract: ``BoundaryTraceLaw`` was no longer a
   :class:`~orpheus.numerics.operator.LinearOperator` subclass,
   no concrete law carried ``apply`` / ``apply_transpose``
   methods, and none reported the operator predicates
   ``is_invertible`` / ``is_adjointable``.
   The §16A.3 three-layer split (descriptor / realizer / operator)
   was enforced by the **type system**, not by convention —
   ``law.apply(psi)`` on a raw law was an ``AttributeError`` at
   runtime and a static-type error at the linter level. The
   ``SNBoundaryRealizer`` was the
   sole bridge from descriptor to callable; see
   :ref:`bc-trace-law-descriptor-model` for the design rationale
   and the predecessor approaches that were tried and rejected.

What took the two factor tiers' place is one split of a different kind. The
:math:`G` tier, an isometry that carries a partner face onto the face, is the
deck of ``orpheus/geometry/deck.py`` (a face role, geometry's business); the
:math:`R` tier and the source :math:`q` are the response network and the inflow
of one shared exterior type in ``orpheus/transport/boundary_conditions/``
(``WallCondition``, an affine pair whose linear part lies in the sub-Markov simplex). Both are realized through
one table per method, ``SN_REALIZERS`` for SN, into typed operators
:math:`\Gamma_+ \to \Gamma_-` whose facts are properties of their types
(:ref:`bc-worked-example`, Steps 4 and 5).

Seven concrete laws shipped under ``orpheus.geometry.boundary``,
one per submodule. The Grand Report v3 vocabulary was used verbatim
for the class names; the pre-refactor names were deprecated aliases
only until Wave O step O.4a.1 retired them (see :ref:`bc-naming-audit`
for the historical index).

The table below covers the **six affine trace laws** — those that
decomposed into the :eq:`affine-bc-form` factors
:math:`(G_\alpha, R_\alpha, q)`. The seventh,
``ZeroFluxBoundary``, was deliberately
absent: :math:`\phi_\Gamma = 0` is a *relation* between the two
traces, :math:`A_-\gamma_- + A_+\gamma_+ = 0`, not a map from one to
the other, so it had no :math:`(G, R)` pair to tabulate. It rode the
same ABC as a pragmatic placement — the diffusion realizer collapsed
it to the albedo-family :math:`\mathcal{A} = -1` (outside the
sub-Markov range by construction) and the SN realizer refused it
outright, with a hand-written ``isinstance`` guard rather than a type
distinction (issue **#177** was to give the relation tier its own type).
The successor gives it that type: ``ZeroFluxCondition`` in
``orpheus/diffusion/zero_flux.py`` is a method-only ``BoundaryCondition``, and the
SN table refuses it on the state-cone axis because its return fraction is
:math:`-1` (:ref:`bc-worked-example`, Step 4).

.. list-table:: Concrete ``BoundaryTraceLaw`` subclasses (retired)
   :header-rows: 1
   :widths: 17 14 18 11 11 29

   * - Class
     - Registry key
     - :math:`G_\alpha`
     - :math:`R_\alpha`
     - :math:`q`
     - Trace-edge family
   * - ``VacuumInflow``
     - ``"vacuum"``
     - identity (fixes no geometry)
     - 0
     - 0
     - **none** — the inflow is data
   * - ``ReflectiveBoundary``
     - ``"reflective"``
     - axis-reflection permutation
     - albedo
     - 0
     - same-face back-edge, mirror-partner map
   * - ``WhiteBoundary``
     - ``"white"``
     - identity (fixes no geometry)
     - albedo × cosine-weighted hemispheric average (Lambertian) —
       **rank-one**
     - 0
     - same-face back-edge, all-to-all on the face
   * - ``PeriodicBoundary``
     - ``"periodic"``
     - spatial wrap along ``axis``; the realizer derived the partner
       face from the installation face
     - 1
     - 0
     - opposite-face **pair**, mutually feeding
   * - ``AlbedoBoundary``
     - ``"albedo"``
     - identity (fixes no geometry)
     - albedo — **magnitude only**; the angular re-emission closure was
       the gap (see the census note below)
     - 0
     - same-face back-edge for :math:`\alpha \neq 0` (degenerated
       to none at :math:`\alpha = 0`); not modelled in
       ``sweep_acyclicity``
   * - ``PrescribedInflow``
     - ``"prescribed_inflow"``
     - identity (fixes no geometry)
     - 0
     - :math:`q \in \Gamma_-`
     - **none** — the inflow is data

The last column is the **structure** each law contributed to the
:math:`(\text{face}, \text{ordinate})` trace digraph, which *is*
intrinsic to the law. Whether the resulting digraph carries a
**cycle** is not: that depends on the other faces, and only the
strongly-connected-component decomposition described under
:ref:`bc-sweep-cycle` can answer it. A same-face back-edge is a
*forward* edge when it is the only one; it closes a loop only when
a second law feeds back the other way. The opposite-face pair is
the exception — periodic closed a loop from a single law. The classification
outlives the laws: it is a statement about the realized maps, and
:ref:`bc-sweep-cycle` states it for the successors.

Per-law rank census: vacuum and prescribed inflow carried no response
kernel (rank-0) — the identity in their :math:`G_\alpha` column was the
*deck* identity, i.e. "this law fixes no geometry", and the vanishing
was :math:`R`'s. White was **rank-1 in angle** — one isotropic re-entry
mode, fed by the cosine-weighted (Lambertian) outflow average, written
to every inflow ordinate; the whole of it lived in :math:`R_\alpha`.
Reflective was an **angular permutation** (the axis reflection, scaled
by the albedo) — rank :math:`N/2` per slab face: structured,
trace-sized, NOT rank-1. The ranks are properties of the maps and carry
over to the realized types: ``ZeroMap`` (rank 0),
``RankOneReemission`` (rank 1 in angle) and ``ScaledOrdinatePermutation``
(an ordinate permutation scaled by a weight in :math:`(0, 1]`).

.. note:: **Albedo's angular closure was the one genuine gap** — the
   *closure* landed at campaign phase **B3.4b**; wiring the law into
   the SN law-admission registry was still open at that point (the registry
   admitted ``{vacuum, reflective}``, and it retired in WF2.4). Its
   :math:`R_\alpha` stated the *magnitude* with which flux returned but
   not the *distribution*, and on an angular trace those are
   independent — the scalar is complete only on a scalar
   (partial-current) trace, where the distribution has no degrees of
   freedom. The closure made
   :math:`\text{albedo}(\alpha, \text{isotropic}) \equiv
   \text{white}(\alpha)` and
   :math:`\text{albedo}(\alpha, \text{specular}) \equiv
   \text{reflective}(\alpha)` **theorems** rather than coincidences —
   in the code both routes executed *one* realization body rather than
   two transcriptions agreeing. The successor states the same facts as
   identities of values: ``WallCondition.specular(a)`` and
   ``WallCondition.white(a)`` are the two chains at weight ``a``, and
   ``a = 0`` of either is the one value ``WallCondition.vacuum()``.

   **Retraction (2026-08-04).** This note used to forecast that the
   specular closure would move its content across into
   :math:`G_\alpha`, "which is exactly what the
   :ref:`membership criterion <bc-factor-roles>` predicts". That was
   **backwards**: it read the criterion's *necessary* half
   (multiplicativity — which a permutation passes) as the whole test.
   The *sufficient* half is the quotient test, and a polished wall is
   not a quotient of the domain, so its specular return is
   **constitutive**. The 2026-08-01 ruling and the shipped code both
   put it in :math:`R_\alpha`
   (``SpecularReemission``), with
   ``AlbedoBoundary.geometry_map`` staying
   ``SelfPairedDeck.identity``
   **unconditionally** — whatever the closure. The ruling stands in the
   successor: a specular wall is the constitutive ``SpecularChain`` of a
   ``WallCondition``, and only a mirror plane the domain continues across
   is a deck.

   What blocked albedo's **domain narrowing** was unchanged and still
   read correctly: without a closure :math:`R = \alpha\,I` was a
   :math:`\Gamma_+ \to \Gamma_+` endomorphism and albedo's
   :math:`G = \mathrm{id}` supplied no crossing, so there was nothing
   to narrow *onto* — which is why the SN realizer **refused** the
   closure-free spelling on an angular trace.

Periodic was a **spatial pushforward** pairing opposite faces. Marshak /
partial-current boundaries were rank-N via the
**descriptor-tree algebra** on the unrealised laws (``LawSum``
/ ``LawScaled`` over ``BoundaryTraceLaw`` leaves) —
realised once per face by
``realize_recursively``; see
:ref:`bc-rank-n-algebra`. Every one of these kernels was **trace-sized**
— tiny against the bulk dimension — which is
exactly the low-rank/structured shape a Woodbury boundary closure
exploits (:ref:`the scoped SMW statement <smw-low-rank-exception>`,
Issue #300); the realized successors are trace-sized in the same way.


.. _bc-naming-audit:

Naming audit: pre-refactor vs Grand Report v3 vocabulary
--------------------------------------------------------

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Wave 7 of the refactor renamed every concrete BC to match the
Grand Report v3 vocabulary verbatim. During the deprecation window
the pre-refactor names were re-exported as deprecated aliases from
``orpheus.geometry.boundary.__init__`` so existing import sites
kept working unchanged. Those aliases were **retired in Wave O step
O.4a.1** once every code and test consumer had migrated; the
canonical names in the middle column were then the sole importable
symbols, until wave 2 retired the whole package (WF2.4). The table is
retained as the historical naming index for readers tracing pre-Wave-O
commits:

.. list-table:: Wave 7 BC renames (pre-refactor name → Wave-7 name; both retired)
   :header-rows: 1
   :widths: 35 35 30

   * - Pre-refactor name (retired Wave O O.4a.1)
     - Wave-7 name
     - Why renamed
   * - ``VacuumBoundaryOperator``
     - ``VacuumInflow``
     - Emphasizes "inflow set to zero", not "operator that vacuums";
       distinguishes from the rank-N case
       ``PrescribedInflow`` which also writes only the inflow
       trace.
   * - ``SpecularBoundaryOperator``
     - ``ReflectiveBoundary``
     - "Specular" is one specific axis-aligned reflection;
       "Reflective" is the family name that the Grand Report uses.
       A future ``SymmetryBoundary`` (deferred) will share the
       reflective-family base but apply on non-physical octant
       boundaries.
   * - ``WhiteBoundaryOperator``
     - ``WhiteBoundary``
     - Drops the redundant "Operator" suffix that pre-dated the
       law / realizer split. The law is no longer "the operator";
       it's the abstract physical statement that gets *realized*
       to an operator.
   * - ``PeriodicBoundaryOperator``
     - ``PeriodicBoundary``
     - Same rationale: pre-refactor "Operator" suffix is
       structurally misleading.
   * - ``AlbedoBoundaryOperator``
     - ``AlbedoBoundary``
     - Same rationale.
   * - ``MixedBoundaryOperator``
     - **retired Wave 11** — see :ref:`bc-rank-n-algebra`.
     - Replaced by the descriptor-tree algebra
       (``LawSum`` / ``LawScaled`` over
       ``BoundaryTraceLaw`` leaves; Issue #186 / B3 + β2);
       the dedicated class added no value over the inherited
       algebra dunders.


.. _bc-sweep-cycle:

Sweep cycles: a configuration property, not a per-law flag
-----------------------------------------------------------

The SN sweep visits cells in a **topological order** of the directed
cell-visit graph, whose edges are oriented by
:math:`\mathrm{sign}(\Omega_n \cdot \hat n_f)`. (The production code
does not *run* a topological sort to find that order — it is closed-form
index arithmetic, ``level_of = local.sum(axis=0)`` in
:mod:`orpheus.sn.loss_representation.sweep_graph` — which is precisely
why nothing in the solver ever builds a digraph and nothing in the
solver can detect a cycle.) For most BCs the boundary is the *root*
of that order — inflow values come from the BC, get propagated
through the cells, and exit as outflow values that the BC consumes
but doesn't feed back. For the wall conditions whose response network
is not empty this is no longer true (a specular wall, a mirror deck and a
diffuse wall at nonzero weight are the same-face case; a translation-wrap
pair is the two-face case):

* **A specular wall, a mirror deck or a diffuse wall at nonzero weight.**
  The outflow flux at a face is mapped to an inflow flux at the **same**
  face (under the reflection permutation for the first two, under the
  rank-one re-emission for the diffuse wall), so the condition adds a
  *same-face trace back-edge*: an outflow slot feeds an inflow slot.
  Whether that edge closes a **cycle** depends on the rest of the
  configuration — see the warning below.
* **A translation-wrap pair.** The edge spans two different faces (outflow
  at face A maps to inflow at face B, **and vice versa**). Because the pair
  of edges is mutually feeding, a translation wrap closes a cycle from a
  *single pair of decks*. The retired SN registry refused it outright (it admitted
  only ``vacuum`` and ``reflective``), so no periodic problem was
  constructible from a tag; the realizer table that replaced the registry
  admits the translation deck, and both faces of a periodic slab realize as
  a ``ScaledOrdinatePermutation`` whose ``domain_face`` is the opposite face
  (`[M]` 2026-10-24, the ``periodic`` listing under Step 4 of
  :ref:`bc-worked-example`). Whether a given sweep schedule can march a
  configuration that carries such a pair is the question of the digraph
  below, not of the realizer.

.. warning::

   **A reflecting face does not, by itself, create a sweep cycle.**
   The proposition "reflective ⟹ the sweep DAG acquires a cycle" is
   *false as stated*, and it is worth stating precisely because the
   opposite reading would forbid something ORPHEUS actually does: the
   curvilinear :math:`r = 0` pole mirror is a specular reflection kept
   **inside** the walk, certified lower-triangular, because the sweep
   visits :math:`\mu < 0` before :math:`\mu > 0` and the mirror feeds
   strictly downstream.

   A cycle requires a closed **loop** — e.g. *both* faces of a slab
   reflecting, so the left mirror feeds the :math:`\mu>0` sweep and
   the right mirror feeds it back. The honest criterion is therefore a
   **strongly-connected-component decomposition** of the
   :math:`(\text{face}, \text{ordinate})` trace digraph, not a boolean
   on the boundary *kind*.

   This is executable, not editorial:
   :mod:`orpheus.derivations.discrete.sn.sweep_acyclicity` is the
   algebra of record, and
   ``tests/gates/sn/sweep/test_sweep_acyclicity.py`` gates it. Measured on
   an S\ :sub:`4` slab:

   .. csv-table::
      :header: left, right, acyclic?, why
      :widths: 16, 16, 12, 56

      vacuum, vacuum, **yes**, no boundary edge at all
      reflective, vacuum, **yes**, one mirror — a forward edge
      vacuum, reflective, **yes**, the mirror image of the above
      reflective, reflective, **no**, "2 SCCs, one per mirror pair"

   A second distinction the same module makes, easy to conflate with
   the first: acyclicity says *some* one-pass order exists, not that a
   **given** one does. A left-reflecting slab is one-pass in the
   :math:`\mu<0`-first order; a right-reflecting slab is equally
   acyclic but needs :math:`\mu>0` first. Triangularity is a property
   of an (operator, **order**) pair.

What the trace digraph replaced: the retired per-law flag
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

   **Retired 2026-07-30.** Until then
   ``BoundaryTraceLaw`` carried a
   ``creates_sweep_cycle`` ``ClassVar``, defaulting to ``False`` and
   declared explicitly on three laws: ``True`` on
   ``ReflectiveBoundary`` and ``PeriodicBoundary``,
   ``False`` (redundantly) on
   ``PrescribedInflow``; the other
   four laws simply inherited the default. It has been removed from
   the ABC, from those three laws, and from the tests that asserted
   its values. The archaeology is kept here because the *shape* of
   the mistake recurs, not because the flag might come back.

Three findings retired it, in increasing order of importance.

1. **It had zero production readers.** The §15A.2 sweep-cycle
   detector it was designed to feed was never built, and the
   invariant named for that job, ``assert_cycles_are_declared``, was
   likewise never implemented. The flag was purely declarative for
   its whole life: a value asserted by tests and read by nothing.

2. **The claim attached to it was false.** Its docstring asserted
   "reflective ⟹ the sweep DAG acquires a cycle" — the proposition
   the warning above refutes. A reflecting face adds a trace
   *back-edge*; a **cycle** needs a closed *loop*.

3. **It could not have worked even in principle** — and this is the
   part worth carrying forward. Whether a face's back-edge closes a
   cycle depends on the **whole face configuration**, not on the
   kind of any one law. The measured truth table above makes the
   point in one line: ``reflective|vacuum`` is acyclic,
   ``reflective|reflective`` is not, and the flag reads ``True`` on
   *reflective* in both. Meanwhile ``periodic|vacuum`` **is** cyclic
   from a single law. So the one value ``True`` was carrying two
   structurally different facts — "this law can take part in a loop
   that *other* faces may close" (reflective) and "this law closes a
   loop by itself" (periodic) — which is the tell that the property
   does not live on the law at all. A ``ClassVar`` on the law
   *class* is evaluated once, at class-creation time, with no
   knowledge of the mesh or of the opposite face; there is no value
   it could have held that would be correct.

The general design rule, which is the reason to keep this record:

.. admonition:: A face declaration may carry only what is intrinsic to it
   :class: tip

   A boundary law was a **descriptor of one face**. It could honestly
   declare its own algebraic content — its :math:`G_\alpha`,
   :math:`R_\alpha`, :math:`q`, and the *family* of trace edge it
   contributed (see the
   :ref:`concrete-subclass table <bc-law-layer>`). It could not declare
   a property of the *configuration* it was placed in. The successor obeys
   the rule by what it stores: a ``WallCondition`` holds a response network
   and an inflow and nothing about orientation or placement (no axis and no
   outward sign: the realizer reads the face's normal from the face name),
   so one condition object is correct on every face of its tag.
   Sweepability, triangularity, and cycle-freedom are all
   configuration properties — the same relinquishment recorded for
   operator names in
   :mod:`orpheus.derivations.discrete.sn.sweep_acyclicity`: an
   operator's name may state what the object **is**, never that it
   sweeps, because a new mesh or a new opposite face can falsify the
   latter without touching the object.

The honest replacement is not a flag but a computation, and it
already exists: build the trace digraph for the *configuration* and
decompose it. Vacuum and a pure inflow are the only closures that
are unconditionally cycle-free, and for a structural reason — they
supply the inflow as **data**, contributing no trace edge at all, so
they cannot participate in any loop (the empty network is the zero
map whatever the source: ``WallCondition.vacuum()`` and
``WallCondition.uniform(2.5)`` both realize as a ``ZeroMap``, `[M]`
2026-10-24, Step 4 of :ref:`bc-worked-example`). A specular or a diffuse
wall at nonzero weight (and the mirror deck, the same operator at weight
one) each add a same-face back-edge: harmless
alone, loop-closing when a second such closure faces it (a pair of white
walls is cyclic for exactly the reason a pair of mirrors is —
the coupling is all-to-all rather than mirror-partnered, but it
still feeds an inflow slot from an outflow slot on the same face).
Periodic pairs opposite faces mutually and so is the one closure that is
cyclic on its own.

The gates on that computation are ``@pytest.mark.foundation``: the
claim is a software/structural invariant of a discrete construction,
not an equation claim, so the tests carry no ``verifies(...)`` marker
(the ``verifies`` ⊥ level doctrine). What they pin is the SCC
decomposition itself, the acyclic ⟹ lower-triangular certificate, and
— as mutation teeth — that dropping the boundary edge *falsely*
certifies a cyclic configuration as acyclic.


.. _bc-realizer-layer-detail:

Realizer (``BoundaryRealizer`` Protocol + ``SNBoundaryRealizer``)
=================================================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

The ``BoundaryRealizer`` Protocol was
``@runtime_checkable`` and lived at
``orpheus.transport.boundary_realization`` (until layer-DAG unit
WP-E1 it was ``orpheus.geometry.boundary._realizer``; realization
builds operators on a method space, which is transport vocabulary, and
geometry only declares). Its contract was one
attribute and one method:

.. code-block:: text

   @runtime_checkable
   class BoundaryRealizer(Protocol):
       method_name: str
       def realize(
           self,
           law: BoundaryTraceLaw,
           method_space: Any,
       ) -> LinearOperator: ...

The Protocol intentionally did *not* prescribe how
``realize`` dispatched over law types — different methods were to
have different optimal dispatch strategies. The SN realizer used
``isinstance`` because the law set was small and stable; a future
realizer that needed runtime extension might have used the
:class:`~orpheus.numerics.registry.RegistryMixin` machinery instead.
The successor drops the Protocol and the ladder together: a method owns one
``RealizerTable``, keyed by the TYPE of the unit it realizes (a chain of the
response network, the zero-network marker, the inflow marker, or a method-only
condition), filled by realizer classes that register by class keyword, and
sealed by the module that defines them, after its last registration. Its one door,
``realize``, calls ``admit`` first, so the answer to "is it admitted" and the
realization cannot disagree (``orpheus/transport/boundary_conditions/realization.py``;
for SN, ``SN_REALIZERS`` in ``orpheus/sn/boundary/realizer.py``).

The Wave 5 SN dispatch table was the documented standard — the §15.2
:math:`G_\alpha` geometric primitives, one per boundary law. The table is
kept as the record of what the ``isinstance`` ladder produced; each cell
describes the ``SNBoundaryRealizer`` as it stood on the last tree that carried
it (2026-10-23), in the past tense:

.. _bc-tensor-primitives:

.. list-table:: SN realization map (law → Wave-0 / Wave-1 primitive), retired
   :header-rows: 1
   :widths: 24 38 38

   * - Law
     - Realized representation (α = 1 fast path)
     - Realized representation (α ∉ {0, 1})
   * - ``VacuumInflow`` — **narrowed**
     - the **zero map** :math:`\Gamma_+ \to \Gamma_-`: a
       :class:`~orpheus.numerics.operator.ZeroMorphism` bound to the two
       half-traces, which derived the forward's :math:`|\Gamma_-|` rows
       and the transpose's :math:`|\Gamma_+|` rows from its **declared
       pair** (the 2026-08-22 S4-amendment; before it, two per-site
       ``codomain_zero`` / ``transpose_zero`` closures carried those
       shapes).  Both spaces were load-bearing: relying on the
       endomorphic ``0.0 * x`` echo would have been right only by accident
       (:math:`|\Gamma_+| = |\Gamma_-|` on every reachable fixture — a
       coincidence, not a contract).
     - n/a (vacuum had no α parameter)
   * - ``ReflectiveBoundary(axis)`` — **narrowed** (B3.4a),
       **bound** (G6.3 step 5)
     - ``PermutationOperator(local_perm, axis=0,
       domain=Γ₊(f), codomain=Γ₋(f)) & IdentityOperator()``
       — a 2-factor
       :class:`~orpheus.numerics.operator.TensorProductOperator` on the
       **reduced** ordinate axis, with
       ``local_perm = Γ₊(f).to_local(π⁻¹[inflow])`` where :math:`\pi` was
       the mirror's derived ordinate permutation
       (:meth:`~orpheus.numerics.quadrature.Quadrature.ordinate_permutation`;
       :math:`\pi^{-1} = \pi` for a mirror) and the remap was the
       half-trace SPACE's own (G6.5 — the space owns its row order).
       Row :math:`j` read the mirror of the :math:`j`-th inflow
       ordinate *at that ordinate's position inside* :math:`\Gamma_+`;
       ``to_local`` was mandatory because a slab mirror **reverses**
       order. The ERR-045 inflow → outflow invariant was *consumed*
       there rather than re-assumed — ``to_local`` raised if the mirror
       sent an inflow ordinate anywhere but :math:`\Gamma_+`. The
       binding made this the **one-link** case of the white arm's chain
       (:ref:`bc-deck-length-one-chain`) and cost zero arithmetic
       — `[M]` ``apply`` and ``apply_transpose`` bit-identical to the
       unbound build on all five shipped quadratures. ⚠ **That reads as
       the two verbs it names, not as "inert".** ``.H`` is a THIRD verb
       and the binding did change it — from the Euclidean transpose to
       the metric-weighted Hilbert adjoint
       (:ref:`tensor-product-spaces`). It was bit-identical
       *for the mirror*, whose metric cancels, and **87 %** away for
       the Lambertian row below.
     - n/a (a symmetry has no amplitude; the bare TP was the
       representation at every axis)
   * - ``WhiteBoundary(axis, outward_sign, α)`` — **narrowed**
       (B3.4a), **factored** (G6.3 step 3b)
     - ``(IsotropicEmissionOperator(...) @ PartialCurrentOperator(...))
       & IdentityOperator()`` — the Lambertian kernel as a two-link
       chain :math:`\Gamma_+ \to S(f) \to \Gamma_-`: the first link
       contracted over :math:`\Gamma_+` to the outgoing partial current,
       the second broadcast it over
       :math:`\Gamma_-`. Both half-traces came from the single
       face-name :math:`\to` signed-projection primitive classified
       against ``TANGENTIAL_EPS``, so the pre-B3.4a private
       ``> 0.0`` outflow test — a second classifier that disagreed with
       the trace space on any quadrature carrying tangential ordinates
       — was gone. The law's declared ``axis`` / ``outward_sign`` was
       cross-checked against the installation face's :math:`\Gamma_+`
       (index SETS, not sizes) before construction.
     - ``α · <that TP>`` (the one-term combination)
   * - ``AlbedoBoundary(α, SpecularReturn(a))`` — **narrowed**
       (B3.4b)
     - the reflective row's body (``_deck_kernel``, shared with
       reflective and periodic), so the law realized to α times the
       mirror ``ReflectiveBoundary(a)``'s matrix (bitwise; the gate was
       ``TestEquivalenceTheorems::test_specular_closure_equals_reflective``);
       α = 1 gave the bare TP, α = 0 the narrowed zero map
     - ``α · <that TP>`` (the one-term combination)
   * - ``AlbedoBoundary(α, IsotropicReturn(a, s))`` — **narrowed**
       (B3.4b)
     - the white row's body (``_checked_angular_average``), so the law
       was ``WhiteBoundary(a, s, α)`` as a matrix; α = 0 gave the
       narrowed zero map
     - ``α · <that TP>`` (the one-term combination)
   * - ``AlbedoBoundary(α)``, closure-free — **refused**
       (``BoundaryError``, B3.4b)
     - — (:math:`R = \alpha I` was an endomorphism of :math:`\Gamma_+`
       and :math:`G` supplied no crossing, so on an angular trace nothing
       paired an outgoing ordinate with an incoming one; until B3.4b the
       arm returned the endomorphisms ``ZeroOperator`` /
       ``IdentityOperator`` / ``α·(I & I)``. A scalar method needed no
       closure, so the diffusion realizer took the same object)
     - —
   * - ``PeriodicBoundary`` — **narrowed** (B3.4c); arrow derived
       (G6.3 step 7)
     - ``PermutationOperator(arange) & IdentityOperator()``, bound
       :math:`\Gamma_+(f') \to \Gamma_-(f)` and consuming the
       PARTNER face's :math:`\Gamma_+`. This was the one law whose domain
       was a different face, and that IS what made it a quotient rather
       than a wall: ``PairedDeck.domain_face`` named the partner and
       ``SNBoundaryOperator._face_domains`` supplied it, so the
       pushforward the spec named lived in the CHANNEL and the action on
       the trace was the identity relabelling between two DISTINCT index
       sets — the ordinate permutation the wrap MOTION induced, through
       the same kernel as the mirror's. Earned, not assumed — the
       kernel's membership and size checks certified
       :math:`\Gamma_+(f') \equiv \Gamma_-(f)` (opposite outward
       normals) by construction. Issue **#183** was exactly this gap and
       closed here; the ``PeriodicWrapOperator`` type retired with it,
       and step 7 retired the unbound-identity stand-in (an endomorphism
       :math:`V \to V` can never be an isomorphism between two different
       spaces — the one link the composability check could not police)
     - n/a (periodic had no α parameter)
   * - ``PrescribedInflow(source)`` — the rank-0 **affine** law;
       **narrowed** (B3.4a), **collapsed** (P3)
     - the **zero map** :math:`\Gamma_+ \to \Gamma_-` — literally the
       same :class:`~orpheus.numerics.operator.ZeroOperator` expression
       the vacuum row above built, from the same
       ``_narrowed_zero_operator`` body, stamped
       :attr:`~orpheus.numerics.operator.BlockRole.BOUNDARY` like every
       other law. The law stayed affine; this tier realized its LINEAR
       factor, and for prescribed inflow that factor was zero, so
       **vacuum and prescribed inflow differed only in** :math:`q` — which
       travelled the boundary-source channel instead
       (:ref:`bc-affine-source-channel`). Until **P3** this arm returned
       an ``IncomingSourceOperator`` whose ``apply`` ignored the
       outgoing flux and asked the source spec to fill
       ``(|Γ₋|,) + psi_out.shape[1:]`` — an AFFINE map in a linear slot,
       measured :math:`\lVert B(0) \rVert_\infty = 2.5` at
       :math:`q = 2.5`; the spec was then evaluated at that same shape by
       ``AngularBoundarySourceSink.from_specs``
       one tier out. The dense inflow **mask** had already dissolved
       with the codomain at B3.4a, so :math:`q \in \Gamma_-` held by
       TYPING rather than by an erasure (ERR-047), exactly as vacuum's
       projector dissolved at B3.2. The law was unreachable from a ``BC`` tag
       in production (the SN registry admitted only
       ``{vacuum, reflective}``, #189), declarable only by
       constructing it directly.
     - n/a

The successors of these rows are the three realized types of
``orpheus/sn/boundary/realized_operators.py``. ``ZeroMap`` is the zero map
:math:`\Gamma_+ \to \Gamma_-` of the first and last rows;
``ScaledOrdinatePermutation`` is the reflective, specular-closure and periodic rows
(a mirror deck, a unit specular wall and a translation wrap are the same type, at
weight one); ``RankOneReemission`` is the white and isotropic-closure rows;
``RealizedSum`` is the image of a sum of chains. They build a
``PermutationOperator`` (scaled and role-stamped by the realizer, not a tensor
product with ``IdentityOperator``) and the factored Lambertian
``emission @ contraction`` from spaces read off the trace ladder, so the
``& IdentityOperator()`` tensor-product spelling above is retired with the
ladder's ``isinstance`` arms (see ``ScaledOrdinatePermutation._derive_inner`` and
``RankOneReemission._derive_inner`` in ``realized_operators.py``).

.. note::

   **Every law the S**\ :sub:`N` **realizer admitted was narrowed.** Since
   campaign phase **B3.2** a realized SN law was typed
   :math:`\Gamma_+ \to \Gamma_-` (:ref:`bc-domain-narrowing`). B3.2
   landed ``vacuum`` and ``reflective``; **B3.4a** landed ``white`` and
   ``prescribed_inflow`` (:ref:`bc-narrowing-b34a`); **B3.4b** completed
   ``albedo`` with a re-emission closure and made the closure-free
   spelling a refusal; **B3.4c** built periodic's partner-face channel
   and G6.3 step 7 derived its arrow :math:`\Gamma_+(f') \to
   \Gamma_-(f)`. The strict xfails that pinned the deferred rows were
   retired (``tests/gates/sn/operators/test_b3_domain_narrowing.py``). `[M]`
   2026-09-22: on a face method space the two closure rows realized
   ``array_equal`` to the reflective and white rows at α ∈ {0, 0.5, 1}
   (``scratch/_definitions/round2/minimal_probe.py``; the reflective row
   then carried an albedo, and since its retirement on 2026-10-01 the
   specular row equalled α times the mirror, `[M]` 20 of 20 bitwise, see
   :ref:`bc-factor-quotients`). A shape assertion
   could not tell the two typings apart, since
   :math:`|\Gamma_+| = |\Gamma_-|` on every quadrature × face in the
   tree, so the *declared spaces* were the thing to read, not the output shape.

The mirror, which had no amplitude, and the α = 1.0 white and albedo
laws returned the **bare** primitive: :math:`1\cdot X` IS :math:`X` by the
combination's one-term collapse (:eq:`combination-zero-and-unit`), so no
coefficient was multiplied in. This was load-bearing for bit-identity:
without it, the mirror
``ReflectiveBoundary`` (pre-refactor
``SpecularBoundaryOperator(axis="x", albedo=1.0)``) would have shifted by
one ULP under the realizer relative to its pre-refactor
``np.take(psi_out, reflection_index, axis=0)`` body — see
:ref:`bc-numerical-evidence` for the bit-equivalence pin (whose
``specular_x_lebedev17`` row is still ``assert_array_equal``, now
against the mirror isometry rather than a recording). The narrowing
preserved that bit-identity and was **gated against the retired
expression, not against the new code called twice**: the reference was
materialised in numpy off the law *descriptor*
(``α·np.take(ψ, mirror_partner, 0)[inflow]``, the partner map from the
independent geometric reference in ``tests/_harness/references.py``,
transpose via ``argsort``), compared with
``np.array_equal`` over slab-asymmetric, slab-symmetric, sphere,
cylinder ``product(2,4)`` and 2-D Cartesian ``level_symmetric(4)``.
The gate was falsified independently — forcing the naive ``arange``
remap reddened 6 rows across slab and sphere in both directions, with a
positive control confirming 26 interceptions.

The :class:`~orpheus.sn.mesh.method_space.SNMethodSpace` dataclass was the
realizer's second argument. At the last tree that carried the ladder it held:

* ``quadrature`` — the angular quadrature (mandatory).
* ``face`` — the face-name label (``"xmin"`` … ``"zmax"``) so the vacuum
  branch could look up the right inflow indices. (The pre-C4 ``"left"`` /
  ``"right"`` spellings were aliases of ``"xmin"`` / ``"xmax"``; since the C4
  face-name carve — :ref:`bc-face-name-carve` — every face was keyed by
  its canonical ``"{axis}{min|max}"`` name.)
* ``inflow_indices`` — the per-face inflow ordinate indices: a realized law's
  **codomain** :math:`\Gamma_-` (derived from the held trace at ``for_face``
  time).
* ``outflow_indices`` — its **domain** :math:`\Gamma_+`, added at campaign
  phase **B3.2** as the sibling of ``inflow_indices``. The codomain had been a
  *field* precisely so a hand-built space could name its own trace without a
  whole trace space; a law's domain deserved the same. The realizer raised a
  ``BoundaryError`` naming ``for_face`` when a narrowed law was realized
  without it.
* ``mesh`` and ``trace`` — the (optional) spatial mesh and the single unified
  :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` for any
  realizer branch that needed more than the per-face index list. Since
  the #205 / #201 unification this was **one** ``trace`` attribute, not a
  separate ``inflow_trace`` / ``outflow_trace`` pair — inflow and
  outflow are selectors over the same trace space (see the
  one-space-two-selectors note above). The ``mesh`` slot was optional
  metadata (C5.3, #225): nothing in the realizer chain read it, and an
  axis-native ``SNProblem`` with no legacy mesh adapter passed ``None``.

The present class holds only ``quadrature`` and ``trace``
(:class:`~orpheus.sn.mesh.method_space.SNMethodSpace`; the face is the
``FaceSite`` the table receives beside it, :ref:`bc-worked-example` Step 4).

``SNMethodSpace.for_face`` was the standard construction site inside
``SNProblem.realize_boundary_law`` (per face, driven by the shared
``resolve_boundary_conditions`` body, #290 P7b), and ``SNMethodSpace.minimal``
returned a quadrature-only method space for unit tests that did not need mesh +
face metadata. **Since B3.2** ``minimal`` was a *partial* constructor:
a quadrature alone cannot name a particular face's :math:`\Gamma_+`,
so it did not suffice for a narrowed law. B3.4a widened that from two
laws to four, leaving ``albedo`` and ``periodic`` as the two that still
realized from a ``minimal`` space because they awaited narrowing, and
that paragraph predicted that after B3.4b / B3.4c it would realize
nothing at all. It realized nothing: `[M]` 2026-09-22, 0 of 10 law
spellings on ``SNMethodSpace.minimal(gauss_legendre(8))``, against 8 of
10 on a face method space
(``scratch/_definitions/round2/minimal_probe.py``). What it still served
was the refusal side: 6 test modules constructed it to assert that a law
could not be realized without face data, and no production code called it.
Both factories retired in WF2.4 with the fields above.


.. _bc-dual-registry:

The law registry (and the realizer registry that was dissolved)
===============================================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

ONE registry connected the tag layer to the law layer:

**Law registry** — keyed by ``BC.kind`` string
(``"vacuum"``, ``"reflective"``, ``"white"``, ``"periodic"``,
``"albedo"``, ``"prescribed_inflow"``, ``"zero_flux"``). The
registry WAS ``BoundaryTraceLaw.registry`` (a class-level dict
maintained by :class:`~orpheus.numerics.registry.RegistryMixin`).
Concrete laws self-registered at module import time via the
``key=`` class-creation kwarg:

.. code-block:: text

   class VacuumInflow(BoundaryTraceLaw, key="vacuum"):
       ...

Lookup was ``BoundaryTraceLaw.create("vacuum")`` or direct
dictionary access ``BoundaryTraceLaw.registry["vacuum"]``. Each
method-mesh additionally carried its own **admission table**
(``BOUNDARY_OPERATOR_REGISTRY: dict[str, type[BoundaryTraceLaw]]``
— the subset of laws its realizer could honestly realize, e.g.
``zero_flux`` was diffusion-only), and the shared
``resolve_boundary_conditions`` body
used THAT table to recover the law class from a mesh-declared
``BC`` — an unsupported tag refused at
phase-space construction with the method's supported list.

**The realizer registry was dissolved at #290 P7b.** The Grand
Report §16A.11 design (lines 3252–3257) paired the law registry
with a second, method-name-keyed ``BoundaryRealizerRegistry``
(``"SN"``, ``"MoC"``, …; realizers self-registered via a decorator
at import time). It shipped in Wave 5 and was retired when the
second functional realizer (diffusion, #290 P3) made the real
consumption pattern visible: **no consumer ever resolved a realizer
by method name**. Production held a method-mesh, and the mesh's
``realize_boundary_law`` arm — the per-method hook of the
:class:`~orpheus.transport.method.TransportMethod` Protocol —
instantiated its own realizer directly; the rank-N walker took the
realizer as an explicit argument. The string indirection carried a
real hazard class for zero payoff: a registry populated by import
side-effects is EMPTY in a fresh process until the right module
happens to be imported, a timing miss invisible to in-suite tests
(process-global state masks it). Dissolving the registry deleted
the hazard class instead of gating it.

Wave 2 retired the other registry too. The string-keyed law registry, the
per-mesh ``BOUNDARY_OPERATOR_REGISTRY`` tables (SN and diffusion) and the
tag-to-law parse are gone, and what answers "is it admitted" and "realize it"
is one realizer table per method (``SN_REALIZERS``, ``DIFFUSION_REALIZERS``,
``CP_REALIZERS``, ``MOC_REALIZERS``, ``MC_REALIZERS``), keyed by unit type and
sealed. The registration-timing hazard named above is answered there by a
rule stated in the table's own docstring: the table and every realizer that
registers into it live in ONE module, which calls ``seal()`` after its last
registration, so
anyone holding the table has executed every registration, and an unsealed table
refuses ``admit``, ``realize``, ``require`` and ``catalogue``.

The two extension axes the dual-registry design named survived
without it:

* Adding a new BC type meant adding one
  ``BoundaryTraceLaw`` subclass with a ``key=`` registration
  (+ the admission-table entries of the methods that supported it).
  Each existing realizer added a dispatch branch; no existing law
  changed.
* Adding a new transport method meant minting one method-mesh
  (structurally conforming to ``TransportMethod``), one
  ``BoundaryRealizer``, and one admission table. **M** realizer
  branches needed to be implemented (one per admitted BC), but no
  existing law changed — and no central registration step existed.


.. _bc-cross-method-stubs:

Adopting the architecture in a new method (and the retired stubs)
-----------------------------------------------------------------

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

A method adopted the unified BC architecture by shipping three
pieces (the diffusion adoption at #290 P3–P7b was the worked
example):

1. a **method-mesh** (``DiffusionMesh(MaterialMesh)``) carrying the
   method's trace + an admission table + the
   ``realize_boundary_law`` arm — it conformed structurally to
   :class:`~orpheus.transport.method.TransportMethod` and got the
   whole tag → law → realized-``bc``-dict pipeline from the shared
   ``resolve_boundary_conditions``
   body;
2. a **realizer** (``DiffusionBoundaryRealizer``) mapping each
   admitted law onto the method's operators;
3. a **method space** (``DiffusionMethodSpace``) carrying whatever
   discretization metadata the realizer read.

The same TransportMethod Protocol survives, and its hook is now
``realize_boundary_role`` with the shared body ``realize_boundary``
(:func:`~orpheus.transport.method.realize_boundary`). The three pieces
above are replaced by one table per method; the crosswalk
(:ref:`bc-wave2-crosswalk`) states each replacement.

.. note::

   **Historical: the Wave-5 stub scaffolding (retired #290 P7b).**
   Between Wave 5 and #290, ``MoCBoundaryRealizer`` /
   ``MCBoundaryRealizer`` / ``CPBoundaryRealizer`` (and, until #290
   P3, ``DiffusionBoundaryRealizer``) existed as
   ``NotImplementedError`` stubs auto-registered by each method's
   ``__init__.py``, "holding the dispatch architecture in place."
   With the registry dissolved there was no dispatch table to hold a
   place in — a stub realizer that cannot realize anything served no
   consumer — so the three stub modules, their auto-import lines,
   and their stub-invariant tests were deleted. MoC / MC / CP kept
   their legacy solver-side BC validation until each adopted the
   architecture per the recipe above.

The SN realizer was **not** auto-imported by ``orpheus.sn.__init__``
(it was a heavy module that every SN consumer paid for); the
:class:`~orpheus.sn.problem.SNProblem` construction imported
it explicitly when it needed it. ``orpheus.sn.boundary`` now imports its table,
``SN_REALIZERS``, on package import.


.. _bc-worked-example:

Worked example — end to end
===========================

The example walks one declaration from the geometry to the sweep: a face role
declared on the geometry (an integer tag or a deck), the tag bound to a
``WallCondition`` in ``BoundaryConditions``, the SN problem realizing the binding
through ``SN_REALIZERS``, the realized operator, and its consumption by the
sweep's face action. The blocks are one script run in order in a single
namespace (each block uses the names of the earlier ones) and every printed line
is the real output of that run, `[M]` 2026-10-24, on the tree of this page's
commit. The example is a 1-D Cartesian slab; the same chain runs on a 2-D
Cartesian mesh with the faces ``xmin`` / ``xmax`` / ``ymin`` / ``ymax``. The
tree before wave 2 ran a different chain for the same slab (``BC("vacuum")``,
a law class looked up in a registry, ``SNBoundaryRealizer.realize``, a shim
pairing the operator with its law); its steps are recorded as history in the
notes under the steps they replaced.

Step 1 — declaration on the geometry
------------------------------------

Every boundary face of the geometry carries exactly one **role**: a deck (the
domain continues across the face: a mirror plane or a translation wrap) or an
``int`` tag (the domain ends at the face, and the problem says later what the
exterior is). The geometry says which faces share an exterior model and never
what the model is: a tag is a declaration, not an observation. The slab below
puts tag ``1`` on its left face and a mirror plane on its right face, and the
mesher carries both onto the mesh's first and last face, keyed by the face names
``xmin`` and ``xmax``:

.. code-block:: python

   from orpheus.geometry import StructuredGeometry
   from orpheus.geometry.deck import MirrorDeck
   from orpheus.mesh import CellsByCount, Mesher

   geometry = StructuredGeometry.slab((0.0, 1.0), (0,), left=1, right=MirrorDeck())
   mesh = Mesher(geometry).partition(CellsByCount.uniform_width(10)).mesh
   print(geometry.face_roles, sorted(geometry.tags))
   print(dict(mesh.face_roles))

.. code-block:: text

   (1, MirrorDeck()) [1]
   {'xmin': 1, 'xmax': MirrorDeck()}


``FaceRole = DeckTransformation | int`` is defined in ``orpheus/geometry/deck.py``
(the parse of a declaration into roles is ``parse_face_roles``; the set of tags a
carrier holds is ``tags_of``), and the mesh's roles are the frozen mapping
``orpheus/mesh/face_roles.py`` builds over the face inventory. The mesh is
method-agnostic: it holds no condition.

.. note::

   **History: the tag was a** ``BC`` **and the geometry carried the law's
   name.** Until wave 2 the same step read
   ``StructuredGeometry.slab((0.0, 1.0), (0,), left=BC("vacuum"),
   right=BC("reflective"))`` and ``mesh.face_laws`` held
   ``{"xmin": BC("vacuum"), "xmax": BC("reflective")}``. ``BC`` was a thin
   wrapper ``BC(kind: str, params: Mapping[str, float])`` whose kind named a
   boundary law, so the geometry already said what the exterior was and a
   misspelt kind was a runtime lookup failure of a later step. Its equality and
   hash were its content (:ref:`structured-geometry-content-identity`).

Step 2 — binding the tags to conditions
---------------------------------------

The method (or the user) supplies ``BoundaryConditions``, the exact twin of the
materials: it maps a **tag** (an ``int``) to a condition, as ``Materials`` maps a
material id to a material. ``bind`` joins it to a carrier's faces, where the
carrier is the object that owns the faces AND their colouring symmetry (the
reflection :math:`x \mapsto a + b - x` when the layers and the two end roles are
mirror-symmetric; ``orpheus/geometry/colouring.py``). The binding is lawful when
it is an exact cover of the tags (every tag a face carries has a condition and
every condition's tag is carried by a face) and when each condition is invariant
under the colouring group; the listing shows the binding, the one-value
identity of the empty network, and the two exact-cover refusals:

.. code-block:: python

   from orpheus.transport.boundary_conditions import BoundaryConditions, WallCondition, bind

   conditions = BoundaryConditions({1: WallCondition.vacuum()})
   binding = bind(mesh, conditions)
   print(dict(binding.decks), dict(binding.tags), binding.faces)
   print(binding.condition_of("xmin") == WallCondition.specular(0.0) == WallCondition.white(0.0))
   for bad in (BoundaryConditions.empty(), BoundaryConditions({1: WallCondition.vacuum(), 2: WallCondition.white(1.0)})):
       try:
           bind(mesh, bad)
       except Exception as error:
           print(type(error).__name__ + ":", error)

.. code-block:: text

   {'xmax': MirrorDeck()} {'xmin': 1} ('xmin', 'xmax')
   True
   UnboundTagError: face 'xmin' carries tag 1 and no condition is declared for it
   UnusedConditionError: a condition is declared for tag 2 and no face carries it (carried tags: [1])


The vacuum condition printed ``True`` against ``specular(0.0)`` and ``white(0.0)``:
the albedo is a coordinate of the response network's cone point, :math:`R = 0` is
one point of the cone however it is reached, and the three spellings are one value
with one digest. A deck needs no condition (a carrier of decks only binds with
``BoundaryConditions.empty()``). The deck of the right face stays in ``binding.decks``,
the tagged face in ``binding.tags``; the face-to-condition join is derived from
``(tags, declared)`` and never stored, so a binding that joins a face to a
condition of another tag cannot be built.

.. note::

   **History: law resolution in** ``SNProblem.__init__``. The same step called the
   shared ``resolve_boundary_conditions`` body, which looked the kind up in
   ``SNProblem.BOUNDARY_OPERATOR_REGISTRY["vacuum"]`` (the class ``VacuumInflow``) and
   instantiated it. The registry carried only ``"vacuum"`` and ``"reflective"``
   because those were the only kinds the production sweep had been wired for
   (white, periodic and albedo were realizable but waited on sweep-side
   plumbing, tracked in separate issues); adding a kind was one dict-entry edit.

Step 3 — the method space and the trace ladder
----------------------------------------------

``SNProblem`` binds the conditions to its own faces at construction (the problem
is the carrier) and builds ONE unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace` for the
whole problem, through the geometry-blind factory
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`,
which takes the angular quadrature and the problem's
:attr:`~orpheus.sn.problem.SNProblem.boundary_face_layout` (the single source of
truth for which faces exist) and nothing else. The trace stores the **signed
projection** :math:`\Omega\cdot\hat n_f` once per face, as an ``(n_faces, N)``
array, not two boolean masks: inflow and outflow are **selectors** over this one
table, derived on demand by the sign of :math:`\Omega\cdot\hat n_f` against
``TANGENTIAL_EPS`` (:ref:`bc-trace-structure`; the one-space-two-selectors design
of Issues #205 / #201). Whether an ordinate is incoming or outgoing at a face is a
predicate evaluated against the same boundary data, not a property of two distinct
domains; one space removes the class of bugs where the inflow and outflow
descriptions of one boundary drift apart.

.. code-block:: python

   import numpy as np
   from scipy.sparse import csr_matrix
   from orpheus.data.macro_xs.mixture import Mixture
   from orpheus.numerics.quadrature import Quadrature
   from orpheus.sn.problem import SNProblem

   mixture = Mixture(SigT=np.array([1.0]), SigC=np.array([0.2]), SigL=np.zeros(1), SigF=np.zeros(1), SigP=np.zeros(1),
                     SigS=[csr_matrix(np.array([[0.8]]))], Sig2=[csr_matrix(np.zeros((1, 1)))], chi=np.zeros(1))
   problem = SNProblem(mesh, Quadrature.gauss_legendre(4), {0: mixture}, boundary_conditions=conditions)
   trace = problem.angular_trace
   print(type(trace).__name__, trace.face_names)
   print(trace.omega_dot_n)
   for face in trace.face_names:
       print(face, "inflow", trace.inflow_indices_for_face(face), "outflow", trace.outflow_indices_for_face(face))
   print(trace.outflow_space("xmin").name, trace.outflow_space("xmin").shape, "->", trace.inflow_space("xmin").name)

.. code-block:: text

   AngularTraceSpace ('xmin', 'xmax')
   [[ 0.86113631  0.33998104 -0.33998104 -0.86113631]
    [-0.86113631 -0.33998104  0.33998104  0.86113631]]
   xmin inflow [2 3] outflow [0 1]
   xmax inflow [0 1] outflow [2 3]
   angular_trace[xmin:outflow] (2, 1) -> angular_trace[xmin:inflow]


The printed projection has the ``xmin`` row :math:`-\mu_x` (outward normal
:math:`-\hat x`) and the ``xmax`` row :math:`+\mu_x`; on ``xmin`` the inflow
ordinates (:math:`\Omega\cdot\hat n < -\epsilon`, that is :math:`\mu_x > \epsilon`:
rightward-pointing) are the indices ``[2 3]`` and the outflow ordinates ``[0 1]``,
and each half-trace is a space named by its face and direction. Neither index set is
the complement of the other: the tangential band belongs to neither
(:ref:`the half-trace <bc-half-trace>`; on ``gauss_legendre(4)`` it is empty). The method space the
realizer receives is
:class:`~orpheus.sn.mesh.method_space.SNMethodSpace`: the quadrature and this trace,
two fields, and no face (the face is the ``FaceSite`` that travels beside it).

.. note::

   **History: the split trace and the indices on the method space.** Before Issue
   #188 and the #205 / #201 unification this step built a per-face
   ``InflowTraceSpace`` / ``OutflowTraceSpace`` *pair* through
   ``InflowTraceSpace.from_mesh_and_quadrature(mesh, quad, faces=("left", "right"))``
   and ``SNMethodSpace.for_face`` took two arguments, ``inflow_trace=`` and
   ``outflow_trace=``. After the unification it took one ``trace=`` argument and
   derived BOTH per-face index sets from it (``inflow_indices`` = the codomain
   :math:`\Gamma_-`, ``outflow_indices`` = the domain :math:`\Gamma_+`), and the
   ``mesh=`` argument was optional metadata that nothing in the realizer chain read
   (C5.3, #225). Wave 2 retired ``for_face``, ``minimal`` and the three index fields.

Step 4 — realization through the table
--------------------------------------

The problem holds its method's realizer table,
:data:`~orpheus.sn.boundary.realizer.SN_REALIZERS`, and the table is the ONE door:
``realize`` calls ``admit`` first and raises ``BoundaryRefused`` iff the admission is
a refusal, so the question "is it admitted" and the realization read one table and
cannot disagree. It is called per face by the shared
:func:`~orpheus.transport.method.realize_boundary` body, through the SN arm
``SNProblem.realize_boundary_role``. The listing calls the table by hand on the left
face, then realizes five spellings of an exterior on it and prints the refusal of a sixth,
a condition that belongs to another method:

.. code-block:: python

   from orpheus.sn.boundary.realizer import SN_REALIZERS
   from orpheus.sn.mesh.method_space import SNMethodSpace
   from orpheus.transport.boundary_conditions import BoundaryRefused
   from orpheus.diffusion.zero_flux import ZeroFluxCondition

   assert problem.boundary_realizers is SN_REALIZERS
   method_space = SNMethodSpace(quadrature=problem.quad, trace=trace)
   site = problem.face_site("xmin")
   print(SN_REALIZERS.admit(WallCondition.vacuum(), site, method_space))
   spellings = {"vacuum": WallCondition.vacuum(), "specular(0.5)": WallCondition.specular(0.5), "white(0.5)": WallCondition.white(0.5),
                "0.3 specular + 0.4 white": WallCondition.specular(0.3) + WallCondition.white(0.4),
                "uniform(2.5)": WallCondition.uniform(2.5), "zero flux": ZeroFluxCondition()}
   for label, declared in spellings.items():
       try:
           print(f"{label:26s}", type(SN_REALIZERS.realize(declared, site, method_space)).__name__)
       except BoundaryRefused as refusal:
           print(f"{label:26s} refused:", str(refusal)[:150].rstrip())
   print({face: type(op).__name__ for face, op in problem.bc.items()})

.. code-block:: text

   Realized(faithful=True)
   vacuum                     ZeroMap
   specular(0.5)              ScaledOrdinatePermutation
   white(0.5)                 RankOneReemission
   0.3 specular + 0.4 white   RealizedSum
   uniform(2.5)               ZeroMap
   zero flux                  refused: sn realizes ['DiffuseChain', 'InflowTerms', 'SpecularChain', 'ZeroNetwork']; not ZeroFluxCondition (STATE_CONE: a return fraction -1.0 < 0 is outside
   {'xmin': 'ZeroMap', 'xmax': 'ScaledOrdinatePermutation'}


The wall is realized as a FUNCTOR on its response network: the empty network goes to
``ZeroMap`` (and so do vacuum, ``specular(0.0)``, ``white(0.0)`` and a pure inflow,
whose linear part is zero: the affine :math:`q` travels the boundary-source channel,
:ref:`bc-affine-source-channel`), one specular chain at weight :math:`s` to
``ScaledOrdinatePermutation``, one diffuse chain at weight :math:`d` to
``RankOneReemission``, and a sum of chains to ``RealizedSum``. A zero-flux condition
is refused: SN's table is keyed by the shared units, so a condition of another method
meets the table's unregistered-unit door, and because its return fraction is
:math:`-1 < 0` and the table's state cone is non-negative flux, the refusal is on the
state-cone axis without SN importing the diffusion type. The right face of the slab is
a deck, realized by the table's one deck realizer: a mirror is certified by
``ordinate_pairing`` (a quadrature that is not closed under the face's mirror is
refused on the angular-resolution axis) and a translation wrap by the equality of its
crossing sets. A periodic slab, whose two faces are translation decks, realizes the
pair of wraps:

.. code-block:: python

   from orpheus.geometry.deck import TranslationDeck

   periodic = StructuredGeometry.slab((0.0, 1.0), (0,), left=TranslationDeck(), right=TranslationDeck())
   periodic_mesh = Mesher(periodic).partition(CellsByCount.uniform_width(10)).mesh
   wrapped = SNProblem(periodic_mesh, Quadrature.gauss_legendre(4), {0: mixture})
   print({face: (type(op).__name__, op.domain_face) for face, op in wrapped.bc.items()})

.. code-block:: text

   {'xmin': ('ScaledOrdinatePermutation', 'xmax'), 'xmax': ('ScaledOrdinatePermutation', 'xmin')}


.. note::

   **History: the realizer was a class with an** ``isinstance`` **ladder.** The step
   instantiated ``SNBoundaryRealizer()`` (stateless; #290 P7b had removed the registry
   lookup that used to sit there) and called ``.realize(law, method_space)``. Its
   vacuum arm fired after the ERR-041 guard cross-checked the claimed
   ``inflow_indices`` against the signed projection the face name alone implies, built
   the zero map with ``_narrowed_zero_operator`` and stamped it with
   ``stamp_boundary_role``. Before B3.2 that arm returned
   ``IncomingOrdinateMaskTensor(inflow_indices=..., n_ordinates=quad.n_ordinates) &
   IdentityOperator()``, a **full-face** operator that zeroed the inflow rows and
   *preserved* the outflow rows, which the consumer then discarded with a slice-write;
   two campaign phases had documented that survival as having "no consumer today". The
   narrowing did not answer that question, it removed it: those rows left the operator's
   domain. The two spaces of the zero map are load-bearing, not ceremony: the endomorphic
   ``ZeroOperator`` returns ``0.0 * x``, an echo of its input's shape, right here only
   by accident because :math:`|\Gamma_+| = |\Gamma_-|` on every reachable fixture. The
   mechanism that carries the shapes was re-spelled on 2026-08-22 (the S4 amendment):
   until then two per-site closures on ``ZeroOperator``, ``codomain_zero`` and
   ``transpose_zero``, passed by each caller, duplicated what the bound spaces already
   knew, and the amendment split the zero into the stateless ``ZeroOperator`` and the
   born-bound ``ZeroMorphism``, which derives both zeros from its declared pair
   (`[M]` ``grep -rn "codomain_zero" orpheus/`` returned four hits at that time, every
   one past-tense prose). ``ZeroMap`` is built on the ``ZeroMorphism``.

.. _bc-step5-pair-with-law:

Step 5 — the facts are properties of the realized type
------------------------------------------------------

Every ``problem.bc[<face>]`` entry has a uniform one-argument ``apply(psi)`` contract
and carries what its consumers ask of it as **properties of its type**, derived from
the type's parameters and never stored beside them: ``return_fraction`` (the partial
current :math:`J^-/J^+` of the realized map), ``permutes_ordinates`` (a scaled ordinate
permutation that is not the identity), ``domain_face`` (the face whose outflow it
consumes: the face itself, or the opposite face for a translation wrap), ``is_vacuum``
(the type is the zero map), ``closes_axis`` (a mirror loop on its own face that returns
everything), and ``domain`` / ``codomain``, the two half-trace spaces. The consumers (the
sweep schedule's reflective set, the solver's leakage face list, the diffusion-synthetic
low-order admission check and the question of which axes close into a mirror loop: `[M]`
2026-10-24, ``grep -rn -e .permutes_ordinates -e .is_vacuum -e .closes_axis -e .return_fraction orpheus/sn``
returns, besides the definition site ``realized_operators.py``, readers in ``sweep_schedule.py``, ``dsa.py``, ``solver.py``, ``operators/boundary.py``
and ``problem.py``) ask the operator what the face DOES and not which declaration it came from:

.. code-block:: python

   import dataclasses
   for face, op in problem.bc.items():
       print(face, type(op).__name__, op.domain.name, "->", op.codomain.name)
       print("   return_fraction", op.return_fraction, "permutes_ordinates", op.permutes_ordinates,
             "domain_face", op.domain_face, "is_vacuum", op.is_vacuum, "closes_axis", op.closes_axis)
   unit_specular = SN_REALIZERS.realize(WallCondition.specular(1.0), problem.face_site("xmax"), method_space)
   print("mirror deck == unit specular wall:", problem.bc["xmax"] == unit_specular)
   half = dataclasses.replace(unit_specular, alpha=0.5)
   print("replace(alpha=0.5) rebuilds the operator:", half.return_fraction, half.apply(np.array([[3.0], [4.0]])).ravel())

.. code-block:: text

   xmin ZeroMap angular_trace[xmin:outflow] -> angular_trace[xmin:inflow]
      return_fraction 0.0 permutes_ordinates False domain_face xmin is_vacuum True closes_axis False
   xmax ScaledOrdinatePermutation angular_trace[xmax:outflow] -> angular_trace[xmax:inflow]
      return_fraction 1.0 permutes_ordinates True domain_face xmax is_vacuum False closes_axis True
   mirror deck == unit specular wall: True
   replace(alpha=0.5) rebuilds the operator: 0.5 [2.  1.5]


Two consequences are visible in the listing. A mirror deck and a unit specular wall
realize the SAME operator (``True``): a symmetry has no amplitude, and the operator is
a function of the realized map, not of how the face was declared, so the kernel of the
loss operator, the proven mirror row of the diffusion-synthetic low-order system and
the reflective face list all read one fact. And ``dataclasses.replace`` of the
weight rebuilds the operator (``alpha = 0.5`` gives return fraction ``0.5`` and half of
the gathered values) instead of editing a fact beside an operator that disagrees: the
type holds only its parameters, and the operator is built from the half-trace spaces
read off the trace ladder, so a lying ``alpha`` is unspellable. The types are frozen
and compare and hash by content (the type, the face, the partner, ``alpha``, the
permutation's index bytes and the two half-spaces' ordinates and metric weights).

.. note::

   **History: step 5 paired the operator with its law.** The step wrapped the realized
   operator in the internal shim ``_BoundBoundaryOperator(realized, law)``, a strict
   one-argument passthrough (Wave 9 migrated 13 production sites from two-argument to
   one-argument; C4 / #220 re-keyed ``bc_<face>`` attributes to the face-name-keyed
   ``bc`` dict). Campaign phase B2.0 made it hold the **law** the operator was realized
   from: until then it kept a copy of ``law.key`` (a string) and nothing else, so
   ``problem.bc[face]`` could say what it was declared as and not what its law DOES,
   and five production sites (``sweep_schedule``'s reflective set, the two ruled-corner
   gates, ``solver``'s leakage face list, the DSA low-order admission check) re-derived
   the answer by comparing the string against literals. Phase B2.2 repointed them to
   ``law_permutes_ordinates(bc[face].law)`` (spelled
   ``bc[face].law.geometry_map.permutes_ordinates`` inline until B3.4b collapsed the four
   sites into one function that asks both factor tiers, because a polished wall's
   specular return lives in :math:`R`), ``bc[face].law.response_kernel.is_zero`` and
   ``bc[face].law.response_kernel.amplitude`` (named ``.scalar`` until B3.0 minted the
   kernel tier: with a rank-one ``LambertianReemission`` in play "scalar" was actively
   wrong, while "amplitude" was true of both). The tag frozensets those sites had keyed
   on (``_RULED_CORNER_KINDS``, ``_SUPPORTED_BC``, both ``{"vacuum", "reflective"}``)
   retired with them. The shim's ``kind`` was a read-through of the law's registry key
   (7 of 7 registered laws agreed with their class's own ``kind``, `[M]` 2026-10-01,
   by walking ``BoundaryTraceLaw.__subclasses__``), load-bearing for the
   ``problem.bc["xmin"] == "vacuum"`` string-equality surface until phase B2.2 retired
   it; the one override that made the two differ, and why it was dangerous, is in
   :ref:`bc-history`. ``law`` was required, because a realized boundary law that cannot
   say which law it realizes is precisely the state B2 existed to delete. The shim was
   internal to its package (not in ``__all__``; a test pinned that). Before Issue #176
   / #186 its Wave-8/9 form carried an optional ``quadrature=`` kwarg that bound an
   ``AngularQuadrature`` and forwarded ``inner.apply(psi, bound_quad)`` to a legacy
   two-argument law body, a dual mode that existed only because the trace factory (then
   ``InflowTraceSpace.from_mesh_and_quadrature``) raised ``NotImplementedError`` for a
   curvilinear ``Mesh1D`` and forced the per-face resolution (then
   ``SNProblem._resolve_one``) to bypass the realizer for spherical and cylindrical
   meshes. Issue #188 lifted that deferral, #176 (C176.1) dropped the bound-quadrature
   mode, and #186 (C-B3.4) dropped the residual ``*_extra, **_kw`` argument swallow.

Step 6 — consumption by the sweep
---------------------------------

The sweep does not apply a realized operator to a whole face slot. The face action
is the **composition** of the realized operator with the restriction of the DOMAIN
face's outflow, ``face_action = law @ gamma_out``, and the result is scattered into
the face's inflow rows by the transpose of the inflow restriction. This is what
``SNBoundaryOperator._reflect_trace`` does for every face of every call (the sum with
the sibling :math:`-B` is :ref:`bc-extraction`); building the product runs the
composability check, so an operator handed the wrong face's half-trace raises at
construction instead of computing a plausible wrong answer (the B3.4c defect, invisible
to a shape check because :math:`|\Gamma_+| = |\Gamma_-|` on every reachable face).
The listing runs that composition on the problem of Step 3, one group, one value
per ordinate:

.. code-block:: python

   psi = np.arange(1.0, 5.0)[:, None]                    # one group, four ordinates of the face slot
   for face in ("xmin", "xmax"):
       law = problem.bc[face]
       gamma_out = trace.outflow_restriction(law.domain_face)   # gamma_+ of the domain face
       gamma_in = trace.inflow_restriction(face)                # gamma_- of the face
       face_action = law @ gamma_out                            # the composed face action; the composability check runs here
       image = face_action.apply(psi)                           # (|Gamma_+|, ng) -> (|Gamma_-|, ng)
       print(face, "image on Gamma_-:", image.ravel(), " scattered into the face slot:", gamma_in.apply_transpose(image).ravel())

.. code-block:: text

   xmin image on Gamma_-: [0. 0.]  scattered into the face slot: [0. 0. 0. 0.]
   xmax image on Gamma_-: [4. 3.]  scattered into the face slot: [4. 3. 0. 0.]


The vacuum face returns zero on :math:`\Gamma_-`; the mirror face returns the outflow
values at the mirrored ordinates (the ``xmax`` outflow values ``[3. 4.]`` come back
reversed as ``[4. 3.]``, which is why a slab mirror needs the half-trace space's own
row remap, ``to_local``, and not an ``arange``), and both are scattered into the two
inflow rows of a four-row slot with the outflow rows left at zero, so that a non-zero
outflow emission, which would corrupt the outflow-definition residual
:math:`\psi_{\text{out}} - \text{streamed}`, is unrepresentable. The transpose leg is
the same composed operator's ``apply_transpose``: :math:`(\text{law}\circ\gamma_+)^T =
\gamma_+^T\circ\text{law}^T` falls out of the product, so the scatter is over
:math:`\Gamma_+` because :math:`\gamma_+` is the operator that was composed.

.. note::

   **History: what the sweep consumed changed twice.** The Wave-8 sweep called
   ``problem.bc["xmin"].apply(psi_out_full)`` on the WHOLE face slot, an ``(N, ng)``
   array. Wave O steps O.4a.2 / O.4b removed ``bc.apply`` from the sweep for every
   geometry (the reflective coupling is delivered by the sibling :math:`-B`,
   :ref:`bc-extraction`), and campaign phase B3.2 narrowed the law's domain to
   :math:`\Gamma_+` (:ref:`bc-domain-narrowing`). G6.3 step 8 then replaced the
   sequence of ``.apply`` calls by the composition above. The public ``bc[face].apply``
   survives, so a caller can still reach a realized operator directly; it must hand it a
   :math:`\Gamma_+`-shaped argument.

.. warning::

   **Being narrowed and validating one's own domain are separate properties.**
   Measured, and still true with the new names: fed a full-face input instead of
   :math:`\Gamma_+`, ``ZeroMap`` and ``ScaledOrdinatePermutation`` return
   :math:`|\Gamma_-|` rows with **no raise** (the zero map's rows are zeros whatever it is fed, so
   its answer is right by accident; the mirror's are **wrong values**), while ``RankOneReemission``
   refuses, because its first link, ``PartialCurrentOperator.apply``, checks
   ``psi.shape[0]`` against :math:`|\Gamma_+|`:

   .. code-block:: python

      wrong = problem.bc["xmax"].apply(psi)                  # the whole face slot, not Gamma_+: no raise
      print("mirror fed the whole slot:", wrong.ravel(), "(the correct image above is [4. 3.])")
      white = SNProblem(mesh, Quadrature.gauss_legendre(4), {0: mixture}, boundary_conditions={1: WallCondition.white(0.5)})
      try:
          white.bc["xmin"].apply(psi)
      except ValueError as error:
          print("white fed the whole slot:", str(error).split(". ")[0])

   .. code-block:: text

      mirror fed the whole slot: [2. 1.] (the correct image above is [4. 3.])
      white fed the whole slot: PartialCurrentOperator.apply: psi has 4 rows along the leading axis, expected |Γ₊| = 2


   The mirror of the listing returned ``[2. 1.]`` where the correct image is ``[4. 3.]``;
   the white face raised. The construction guard lives on the trace restriction
   operator and does not travel to the operator the realizer *emits*. The wrong-shape
   call is unreachable through ``_reflect_trace``, which always feeds a guarded
   ``gamma_out.apply(...)``, and **reachable through** ``problem.bc[face].apply``.
   Both are narrowed; only one validates. Never read "narrowed" as "guarded", and never
   credit the weaker property as the stronger one. (Before the carve the same split
   held across vacuum's ``ZeroMorphism`` (a ``ZeroOperator`` before 2026-08-22), reflective's ``TensorProductOperator`` and the
   white chain; prescribed inflow ignored its input entirely, so it had nothing to
   validate against.)


.. _bc-universal-invariants:

Universal ``assert_*`` invariants
=================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Before wave 2 the abstract class ``BoundaryTraceLaw`` declared **five universal** assertion methods plus **three specific** assertions on the laws that needed them, and the aggregate ``assert_realizable`` fired them when a law was realized. The suite ``tests/gates/geometry/test_bc_universal_invariants.py`` exercised them; it retired with them in WF2.4 (commit ``ece6a281``). The split between "universal" and "specific" followed the Grand Report v3 §16A.12 + §27.6 catalog: universal invariants are properties the affine law :eq:`affine-bc-form` should satisfy for **any** physically meaningful boundary condition, and specific invariants pinned properties that only a subset of laws claimed (involution is meaningful for a mirror and not for a white wall).

.. warning::

   **The declared surface was wider than the implemented one** (measured on the tree before the carve, as this page recorded; the methods no longer exist to re-measure): four of the five universal methods had an **empty** base body, and two of those four (``assert_inflow_outflow_classification``, ``assert_outgoing_leakage_unconstrained``) were overridden by nobody, so they were permanent no-ops that asserted nothing about any law. Four of the seven concrete laws (``VacuumInflow``, ``PeriodicBoundary``, ``PrescribedInflow``, ``ZeroFluxBoundary``) overrode no universal invariant at all, and only ``assert_source_is_placeable`` fired for every law. A second gap, orthogonal to the empty bodies: ``assert_realizable`` had exactly **one** production caller, the SN realizer. The diffusion realizer never called it, so the claim that "every law arrives at its primitive construction already certified" held on the SN arm only.

The carve removed both gaps by construction and not by assertion: each transport method (SN, diffusion, CP, MoC, MC) holds one sealed ``RealizerTable`` whose one door answers admission and realization together (:ref:`bc-realize-recursively`), and each invariant that survived is now **parsed** (a value that violates it is unconstructible), **certified** (a function returns the structure that carries it) or **admitted** (the table refuses with a named axis), never a method a caller may forget to call. The table lists, for each retired method, where its invariant lives at the tip.

.. list-table:: The universal invariants of the retired surface, and where each lives now
   :header-rows: 1
   :widths: 22 14 64

   * - Retired method
     - Pinned error
     - Where the invariant lives now, and what it asserts there
   * - ``assert_inflow_outflow_classification``
     - ERR-040
     - **The trace ladder**, :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`: the partition is a property of two SPACES and not an assertion. The inflow half-trace of a face holds exactly the ordinates with :math:`\Omega\cdot\hat n_f < -\epsilon`, the outflow half-trace those with :math:`\Omega\cdot\hat n_f > +\epsilon`, and a tangential ordinate (:math:`|\Omega\cdot\hat n_f| \le \epsilon`, with :math:`\epsilon` the one constant ``TANGENTIAL_EPS``) is in neither (:eq:`ordinate-partition-inflow-outflow`). The error class ``IncomingOutgoingTraceClassificationError`` had no production raise site and was not carried over: no class is minted for it. Catcher of ERR-040: ``tests/gates/numerics/test_angular_trace_space.py::test_axis_aligned_ordinates_excluded_from_both_selectors``.
   * - ``assert_outgoing_leakage_unconstrained``
     - none
     - **Nowhere: the method retired, and its claim became unspellable.** It asserted that the BC does not constrain the outgoing trace flux, and asserted nothing. A realized face operator is a map :math:`\Gamma_+ \to \Gamma_-` that reads the outgoing trace as its input and writes only the inflow trace, so a constraint on the outgoing trace has no slot in its codomain. A future Dirichlet-outflow law would be a different codomain, hence a new type.
   * - ``assert_geometry_map_measure_preserving`` (empty base body; overridden by ``ReflectiveBoundary`` only)
     - ERR-042
     - **The certification** :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing`, keyed on ``(quadrature, motion)`` alone: the trace measure :math:`m_n = w_n\,|\mu_{a,n}|` of the pairing :math:`\pi` satisfies :math:`m_{\pi(n)} = m_n` to relative tolerance ``MEASURE_RTOL`` = 1e-12, else ``PairingNotMeasurePreservingError``. The check is independent of involution, as it was: an involutive table that pairs ordinates of different weight classes passes involution and breaks the measure, the hole the pre-#52 delegation left open. The SN type ``ScaledOrdinatePermutation`` re-checks, on the trace it is given, that the pairing crosses and carries the outflow trace measure onto the inflow trace measure (``tests/gates/sn/architecture/test_realized_boundary.py::test_a_pairing_that_does_not_cross_or_moves_measure_is_refused_at_construction``).
   * - ``assert_response_positive_if_declared`` (empty base body; overridden by ``WhiteBoundary`` and ``AlbedoBoundary`` only)
     - ERR-043
     - **The parse of the wall's cone**, :class:`~orpheus.transport.boundary_conditions.WallCondition`: a negative response weight or a negative inflow weight raises ``BoundaryResponseNotPositiveError`` at construction (``WallCondition.white(-0.1)``), so a wall outside the cone does not exist. The error is no longer a check on a realized kernel's output.
   * - ``assert_source_is_placeable`` (the only universal with a real base body; no law overrode it)
     - ERR-047
     - **The SN source realizer**, ``SNSource.sample``: a space that is not an inflow half-trace is refused with ``BoundarySourceNotOnIncomingTraceError`` (``orpheus.sn.boundary.errors``). The old method asked that a source-carrying law be realized with a nameable :math:`\Gamma_-`, which needed an inflow-index set; at its P6 re-pose it had stopped probing the values of :math:`q` for support on :math:`\Gamma_-` and required that index set instead, which is why the ``NoSource`` default passed; there is no index set now, and the sampler is handed :math:`\Gamma_-(f)` itself and must return exactly its shape (``sample_terms`` raises ``ValueError`` otherwise), so the support of :math:`q` holds by construction and the one residual fact is that the space IS an inflow half-trace. Catcher: ``tests/gates/transport/test_boundary_source_from_specs.py::TestTheRecipeLandsOnTheInflowSlotsOnly``.

.. list-table:: The specific invariants of the retired surface, and where each lives now
   :header-rows: 1
   :widths: 30 14 56

   * - Retired method
     - Pinned error
     - Where the invariant lives now
   * - ``ReflectiveBoundary.assert_is_involutive``
     - ERR-044
     - The second check of ``ordinate_pairing``: :math:`\pi\circ\pi = \mathrm{id}` (``perm[perm] == arange``), else ``PairingNotInvolutiveError``.
   * - ``ReflectiveBoundary.assert_reflection_maps_inflow_to_outflow``
     - ERR-045
     - The third check of ``ordinate_pairing``: every non-tangential ordinate pairs with the opposite sign class on the mirror's axis, else ``PairingDoesNotCrossError``; and, on the SN trace, the gather of ``ScaledOrdinatePermutation``, which raises ``ValueError`` (ERR-045) when a paired ordinate is not in :math:`\Gamma_+`.
   * - ``WhiteBoundary.assert_submarkov`` and ``AlbedoBoundary.assert_submarkov``
     - ERR-046
     - The same cone parse as ERR-043: the return fraction of the wall (the correctly rounded ``math.fsum`` of its weights) is at most one, else ``SubmarkovViolationError`` at construction (``WallCondition.specular(1.2)``). The realized sum of chains refuses a total above one a second time (``RealizedSum``), and a realized leaf takes its amplitude in :math:`(0, 1]`.

The three checks of ``ordinate_pairing`` stay three, because each passes on a table the other two reject (the catalogue's lesson of ERR-045: the identity-table mutant passes involution and measure and reddens only the crossing check). The catchers are ``tests/gates/numerics/test_ordinate_pairing.py::test_the_three_errors_keep_their_teeth`` and ``::test_a_poisoned_quadrature_pairing_reddens_at_use_time``.

``ordinate_pairing`` reads the mirror plane from the motion alone (the axis whose diagonal entry of the motion's linear part is :math:`-1`), so a pair of a reflection about x and ``axis="y"`` is unspellable. An oblique reflection is refused with the base ``OrdinatePairingError``: the per-axis mean cosine its crossing check reads has no oblique analogue yet (#657).

The nature of the claim changed with the mechanism. An ``assert_*`` method was a promise that a caller ran it, and a method that nothing overrode asserted nothing; a parsed value, a returned permutation and a table's refusal each assert themselves, which is the repair that ``instrument-doctrine`` X3 asks for ("return the structure so that it asserts itself").


.. _bc-named-error-catalog:

Named-error catalog (ERR-040..ERR-047)
======================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Per Grand Report v3 §26A.4 and the `vv-principles` skill's "Log every caught bug" directive, every typed error carries a matching ``@pytest.mark.catches("ERR-NNN")`` decorator on the test that proves it fires under the right fault injection. The eight errors were shipped under ``orpheus.geometry.boundary._errors``, a module that retired in WF2.4. The entries themselves, which are the record of each defect and say where its invariant lives at the tip, are in :doc:`/theory/verification/error_catalog`; the table below re-points this page's account of each to what the entry names now.

.. list-table::
   :header-rows: 1
   :widths: 10 25 8 22 35

   * - Entry
     - Trigger condition (unchanged)
     - Mode
     - Retired class
     - What the entry names now
   * - ERR-040
     - Tangential ordinate at a face that requires strict partition.
     - #5 (index)
     - ``IncomingOutgoingTraceClassificationError``
     - No class: it was never raised. The trace ladder's half-spaces exclude the tangential band from both halves (:ref:`bc-universal-invariants`).
   * - ERR-041
     - Vacuum law applied to an outgoing trace.
     - #6 (convention)
     - ``VacuumAppliedToOutgoingTraceError``
     - No class: the realizer guard and the two index fields it cross-checked retired together. The zero map reads its half-spaces off the trace ladder, :math:`\Gamma_+(f) \to \Gamma_-(f)` of its own face, and the invariant is witnessed on the ladder and on the realized operator (two catchers).
   * - ERR-042
     - Geometric map :math:`G` does not preserve :math:`w(\Omega) |\Omega \cdot \hat n|`.
     - #5 + #6
     - ``BoundaryGeometryMapNotMeasurePreservingError``
     - ``orpheus.numerics.quadrature.ordinate_pairing.PairingNotMeasurePreservingError``.
   * - ERR-043
     - Response kernel produces negative output (a sign-flipped construction, e.g. :math:`-\alpha`).
     - #1 (sign)
     - ``BoundaryResponseNotPositiveError`` (old home)
     - The same name at ``orpheus.transport.boundary_conditions.errors``, raised when a ``WallCondition`` is constructed outside the cone.
   * - ERR-044
     - Reflection permutation is not its own inverse: :math:`\pi \circ \pi \neq \mathrm{id}`.
     - #5 (index)
     - ``ReflectionNotInvolutiveError``
     - ``orpheus.numerics.quadrature.ordinate_pairing.PairingNotInvolutiveError``.
   * - ERR-045
     - Reflection maps an inflow ordinate to itself.
     - #5 (index)
     - ``ReflectionDidNotMapInflowToOutflowError``
     - ``orpheus.numerics.quadrature.ordinate_pairing.PairingDoesNotCrossError``. The old mechanism note ("a non-axis-aligned reflection mislabelled as ``ReflectiveBoundary``; the right family would be a future ``SymmetryBoundary``") has no successor type to name: ``SymmetryBoundary`` was never minted (`[M]` 2026-10-24: ``git grep -c SymmetryBoundary -- orpheus tests`` finds no file, and this page is the control that the word is greppable), and an oblique reflection is refused outright by ``ordinate_pairing`` with the base ``OrdinatePairingError`` (#657).
   * - ERR-046
     - Sub-Markov BC with :math:`\alpha > 1`: an albedo or white return fraction above one would imply a source on the boundary surface.
     - #4 (factor)
     - ``SubmarkovViolationError`` (old home)
     - The same name at ``orpheus.transport.boundary_conditions.errors``, raised when a ``WallCondition`` is constructed with a return fraction above one.
   * - ERR-047
     - Boundary source :math:`q` has nonzero outflow entries.
     - #6 (convention)
     - ``BoundarySourceNotOnIncomingTraceError`` (old home)
     - The same name at ``orpheus.sn.boundary.errors``, raised by ``SNSource.sample`` for a space that is not the inflow half-trace of its face.

Before wave 2 all eight extended :class:`ValueError` through one base, ``BoundaryError``, and a catch site could recover the offending law's NAME from ``BoundaryError.law``. The roots now differ by layer, and every one still extends ``ValueError``, so an ``except ValueError`` written against the pre-refactor code keeps working. Transport's ``BoundaryError`` carries ``.declaration``, the refused OBJECT or its type, so a test compares objects (``err.declaration == MirrorDeck()``) and not a spelling; it roots ERR-043, ERR-046 and ERR-047 and the door refusals (``BoundaryRefused``, which carries the method, the face and the admission with its axis and issue; ``UnboundTagError``; ``UnusedConditionError``; ``NonInvariantConditionError``). The three quadrature certifications root in numerics' ``OrdinatePairingError``, and geometry's refusals in ``DeckError``, because neither layer holds a module named for boundary conditions.


.. _bc-rank-n-algebra:

.. _bc-descriptor-tree-vs-operator-tree:

Descriptor-tree algebra for rank-N boundaries
=============================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Rank-N (Marshak, partial-current) boundary conditions are **not** a special class: they are points of one sub-Markov simplex, combined by the arithmetic of one type. Until wave 2 the combination was a closed **descriptor-tree algebra** over ``BoundaryTraceLaw | LawSum | LawScaled`` nodes, pure declarative structure with **no** ``apply`` method on any node. The dunders of the law (``+``, ``-``, ``*``, ``/``, unary ``-``) returned ``LawSum`` and ``LawScaled`` instances, never operators, and ``realize_recursively`` was the sole path from the descriptor tree to the operator tree (:ref:`bc-realize-recursively`). Wave 2 replaced the tree by the arithmetic of the declaration itself, and the replacement is a smaller object with a stronger guarantee:

* **A wall is a point, not a tree.** :class:`~orpheus.transport.boundary_conditions.WallCondition` is the affine pair (response network, inflow). The response :math:`R = s\,P + d\,(E\circ R_t)` has the weights :math:`(s, d)` of its two chains (the specular chain :math:`P`, a permutation of ordinates, and the diffuse chain :math:`E\circ R_t`, a partial-current average followed by a uniform re-emission), and the source :math:`q \ge 0` is the inflow.
* **Admissibility is a value predicate** (decision 50 of the carve plan). The weights live in :math:`K = \{s, d \ge 0,\ s + d \le 1\}` and the inflow weights are non-negative; one function, ``_cone_refusal``, is the only definition of "inside", and ``__post_init__`` raises what it returns. A wall outside :math:`K` (``specular(1.2)``, ``white(-0.1)``, a direct two-chain ``(0.6, 0.6)``) does not exist.
* **The arithmetic exits through one function.** ``+``, ``-``, unary ``-`` and scalar ``*`` act on the pair and return, through ``_combine``, a ``WallCondition`` when the result lies in the cone and a :class:`~orpheus.transport.boundary_conditions.WallDisplacement` otherwise: a point of the affine hull of the cone times the source space, a direction in which a parameter may move (a sensitivity), never a wall and never bindable. The sign of a scalar and the presence of ``-`` decide nothing: ``x + h * delta`` for a finite difference re-enters the cone without re-typing.
* **No tree shape is kept, because none exists.** There are no ``LawSum`` and ``LawScaled`` nodes to flatten; the network is a mapping from chain to weight, held in the canonical chain order, so ``a + b`` and ``b + a`` are one value.

The §15.2 sum-of-tensor-products form

.. math::
   :label: bc-rank-n-tensor-decomposition

   B \;=\; \sum_{\alpha} c_{\alpha}\, G_{\alpha},
   \qquad c_{\alpha} \in \mathbb{R},
   \quad G_{\alpha} \in
   \{\text{permutation, average, mask, wrap, identity, source}\},

.. vv-status: bc-rank-n-tensor-decomposition documented

is the menu of building blocks, and the live type realizes the following members of it. The coefficients :math:`c_\alpha` of the two response chains are the weights of the wall (:math:`s` on the permutation :math:`P`, :math:`d` on the average-then-re-emission chain), and they are not free in :math:`\mathbb{R}`: the cone restricts them to :math:`K`, which is the sub-Markov condition of ERR-046. The source is the ``InflowTerms`` of the wall. The wrap is not a weight of a wall: a periodic pair is a ``TranslationDeck``, realized as a permutation at unit weight. The mask retired with the inflow-only vacuum (:ref:`bc-vacuum-semantic-correction`): vacuum is the empty network, whose image is the zero map. A method-only condition (zero flux, an anisotropic SN kernel) is a separate type of the same base and sits outside this sum.

The standard Marshak boundary (Bell & Glasstone 1970 §1.5), a mix of specular reflection (weight :math:`c_1`) and diffuse white reflection (weight :math:`c_2`), is the wall :math:`c_1\,\mathrm{specular}(1) + c_2\,\mathrm{white}(1)`, a point of the cone when :math:`c_1 + c_2 \le 1`. Its specular part is a wall's constitutive specular return and never the mirror: a mirror is a symmetry of the domain, a deck, and is refused as an operand of the arithmetic (:ref:`bc-deck-law-not-an-operand`).

.. code-block:: python

   import numpy as np

   from orpheus.geometry import CoordSystem
   from orpheus.numerics.face_layout import FaceLayout
   from orpheus.numerics.quadrature import Quadrature
   from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
   from orpheus.sn.boundary.realizer import SN_REALIZERS
   from orpheus.sn.mesh.method_space import SNMethodSpace
   from orpheus.transport.boundary_conditions import FaceSite, WallCondition

   # Declare: the Marshak wall is a point of the cone; nothing is realized yet.
   marshak = 0.3 * WallCondition.specular(1.0) + 0.7 * WallCondition.white(1.0)
   print(type(marshak).__name__, marshak.return_fraction)
   print(marshak == WallCondition.specular(0.3) + WallCondition.white(0.7))
   print(hasattr(marshak, "apply"))       # a declaration is not callable

   # Realize it on the face xmax of a slab, through the SN table.
   quad = Quadrature.lebedev(17)
   layout = FaceLayout.from_named_shapes([(f, (int(quad.n_ordinates),)) for f in ("xmin", "xmax")])
   space = SNMethodSpace(quadrature=quad, trace=AngularTraceSpace.from_quadrature_and_layout(quad, layout))
   site = FaceSite("xmax", CoordSystem.CARTESIAN, 1)
   op = SN_REALIZERS.realize(marshak, site, space)
   print(type(op).__name__, [(type(t).__name__, t.return_fraction) for t in op.terms])
   print(op.domain.shape, op.codomain.shape)
   psi_out = np.random.default_rng(0).random(op.domain.shape)   # Gamma_+ shaped
   psi_in = op.apply(psi_out)                                   # Gamma_- shaped
   print(psi_in.shape)

   # The partition of the 110 ordinates at this face: inflow, outflow, tangential.
   n_in = len(space.trace.inflow_indices_for_face("xmax"))
   n_out = len(space.trace.outflow_indices_for_face("xmax"))
   print(n_in, n_out, quad.n_ordinates - n_in - n_out)

   # Nesting folds at construction; a wall written two ways is one value.
   nested = 0.5 * (0.3 * WallCondition.specular(1.0) + 0.7 * WallCondition.white(1.0))
   print(nested == WallCondition.specular(0.15) + WallCondition.white(0.35))
   print([(type(t).__name__, t.return_fraction) for t in SN_REALIZERS.realize(nested, site, space).terms])
   print(0.0 * WallCondition.specular(1.0) == WallCondition.vacuum())
   a, b = WallCondition.specular(0.3), WallCondition.white(0.5)
   same = np.array_equal(SN_REALIZERS.realize(a + b, site, space).apply(psi_out),
                         SN_REALIZERS.realize(b + a, site, space).apply(psi_out))
   print((a + b) == (b + a), same)

`[M]` 2026-10-24 at the tip ``471c4737``: the block above run as ``PYTHONPATH=<worktree root> python`` from outside the repository (``orpheus.__file__`` resolved into the worktree), its output pasted:

.. code-block:: text

   WallCondition 1.0
   True
   False
   RealizedSum [('ScaledOrdinatePermutation', 0.3), ('RankOneReemission', 0.7)]
   (49,) (49,)
   (49,)
   49 49 12
   True
   [('ScaledOrdinatePermutation', 0.15), ('RankOneReemission', 0.35)]
   True
   True True

The realization is the formal sum of the realized leaves, one per chain, in the canonical chain order, consumable by the SN sweep and Krylov path through the uniform one-argument ``apply``. A declaration and a realized operator are **separate type families** and never inter-compose, as the descriptor tree and the operator tree never did: ``WallCondition.white(0.5) + realized`` is a ``TypeError`` (the usage block of :ref:`bc-realize-recursively`), and the user realizes first.

**Why the tree was replaced.** The tree admitted what the cone forbids. The docstring of ``tests/gates/transport/test_wall_condition_algebra.py`` records two inputs measured at the base tip ``4ec95bc6``: ``LawSum(LawScaled(0.6, spec), LawScaled(0.6, white))`` constructed (a super-Markov wall), and ``(a + b) + c == a + (b + c)`` was ``False`` for three laws, because the tree kept its association. Neither is spellable now: the first is refused at construction, and a wall holds one weight per chain, so a third term merges into an existing chain. The same docstring states the rounding that remains: the two associations of a sum are exact on dyadic coordinates and differ by one ulp for the ``0.1, 0.2, 0.3`` triple, so a wall written two ways can differ in digest, a cache miss and never a wrong answer. The source travels with the response through every combination (the defect of the first draft, an algebra on responses that dropped it).

Closed-algebra guarantees
-------------------------

Each retired guarantee, and what it became:

* **Constant folding on scalars.** ``LawScaled(α, LawScaled(β, x))`` collapsed to ``LawScaled(α * β, x)`` at construction, which kept the tree shallow. A scalar now multiplies the weights of the network, so there is no nesting: the nested wall of the block above IS ``0.15 * specular + 0.35 * white`` (printed ``True``), and a zero weight is not stored (``0.0 * specular(1.0)`` is ``vacuum()``, printed ``True``).
* **No associativity flattening on law sums.** ``(a + b) + c`` was ``LawSum(LawSum(a, b), c)``, distinct from ``LawSum(a, LawSum(b, c))``, and only the OPERATOR flattened them (`[M]` 2026-10-20, a three-leaf tree on Lebedev-17: terms ``[0.3, 0.5, 0.1]`` from both shapes, ``array_equal``). The declaration is flat now. A wall has at most one weight per chain, so a realized sum has at most two terms (``RealizedSum`` refuses a repeated chain), and the realization is summed in the canonical chain order, so it does not depend on how the wall was built (printed ``True True`` for ``a + b`` against ``b + a``: equal as declarations, bit-equal as applied).
* **Subtraction rewrote via** ``LawScaled(-1, ...)``. ``-`` and unary ``-`` are the pair arithmetic: ``WallCondition.white(0.3) - WallCondition.specular(0.5)`` is a ``WallDisplacement`` (the usage block of :ref:`bc-realize-recursively`), and a difference that stays in the cone is a ``WallCondition``. There was no ``LawDifference`` and there is none.
* **Division rewrote via** ``LawScaled(1/α, ...)``. There is no ``__truediv__`` on a wall: it has no consumer, and division would leave the dyadic coordinates in which the algebra is exact.
* **Operands.** A mirror deck and a method-only condition define no operator at all, so ``0.5 * MirrorDeck()`` and ``ZeroFluxCondition() + x`` are Python's own ``TypeError`` (the usage block of :ref:`bc-realize-recursively`), where the law family needed a guard (``_admit_as_operand``, an ``ELEGANCE-DEBT`` tag owed to #551).

The pre-refactor implementations
--------------------------------

Two prior approaches converged on the descriptor-tree design
that wave 2 later retired, through empirical falsification:

**Wave 11 (~2026-03)** — ``MixedBoundaryOperator(components:
list[tuple[float, BoundaryOperator]])`` class whose
``apply`` body looped over ``components`` and summed
``coeff * primitive.apply(psi, quad)``. The SN realizer
dispatched on it via an ``isinstance(law, MixedBoundaryOperator)``
branch that ran the same loop with
``coeff * realize(primitive, ms)`` summed via
:class:`~orpheus.numerics.operator.OperatorSum`. Wave 11 deleted this code because the
delayed-realization-by-container pattern broke down once vacuum
needed per-face inflow indices that the bare-law container had
no access to.

**β1 interim landing (Issue #186 / B3, ~2026-04)** — every
``BoundaryTraceLaw`` inherited the Wave-0 operator-algebra
dunders from :class:`~orpheus.numerics.operator.LinearOperator`,
so writing ``0.3 * spec + 0.7 * white`` directly produced an
:class:`~orpheus.numerics.operator.OperatorSum` of ``ScaledOperator``-wrapped raw
``BoundaryTraceLaw`` leaves (NOT realized). The
``realize_recursively`` walker then traversed the Wave-0
composer tree, realized each leaf, and emitted a parallel
operator tree. This achieved the right algebraic shape but
**violated the type system**: the resulting tree was an
:class:`~orpheus.numerics.operator.OperatorSum` instance (a :class:`~orpheus.numerics.operator.LinearOperator`!) whose
``apply`` could not actually be called before realization —
calling it raised ``BoundaryError`` at apply-time because the
leaves were laws, not operators. The convention "you must realize
this OperatorSum before calling apply" was a runtime contract that
the type system did nothing to enforce. β1 retained the
ergonomic of "the same ``+`` operator before and after
realization" but at the cost of conflating two type families.

**β2 (Issue #186 / B3 + β2, the architecture until wave 2)** — separated the two
type families explicitly: ``LawSum`` / ``LawScaled`` for
the descriptor tree (no ``apply``); the operator algebra's
composites for the operator (with ``apply``; since package S2,
2026-10-20, the formal linear combination
:class:`~orpheus.numerics.operator.OperatorSum`, where at β2 it was the
binary ``OperatorSum`` / ``ScaledOperator`` tree). The static type system enforced "you cannot call
this until it's realized" — ``LawSum`` had no ``apply``
method on the class, so the linter flagged ``tree.apply(...)``
without running the program. The ergonomic of "the same ``+``"
survived because both the law-tree and the operator-tree use the
same Python ``+`` syntax; the runtime dispatch on type tells the
reader (and the type checker) which algebra is in effect.

The reference-image harness verified at the β2 transition that the Marshak case
``0.3 * spec + 0.7 * white`` survived: the realized
combination ``0.3·X + 0.7·Y`` (until S2 the binary
``OperatorSum(ScaledOperator, ScaledOperator)`` tree) reproduced the derived convex combination
:math:`0.3\,R_{\text{spec}} + 0.7\,R_{\text{diff}}` inside the
reduction-order bound, because the operator-tree shape after
``realize_recursively`` was algebraically identical. At the tip the same wall is the
harness's mixed case, ``mixed_30spec_70white_LS4`` of
``tests/gates/sn/architecture/test_wave2_carry.py``, realized through ``SN_REALIZERS``. (Until the
2026-08-01 re-anchoring the claim was the narrower one that the
realized output matched the β1-era *recorded* output; see
:ref:`bc-snapshot-reanchoring` for why the reference moved.)

**Wave 2 (2026-10-24)** — retired the descriptor tree itself. The two type families survive as a declaration (``WallCondition``, whose arithmetic is closed in the cone, and ``WallDisplacement`` outside it) and a realized operator (``RealizedBoundary``), and the table's functor replaces the walker (:ref:`bc-realize-recursively`). Nothing in the walker's role was lost: the tree's shape, which only the walker consumed, is the wall's normal form, and the invariant it could not state, that a wall stay in the cone, is a parse.


.. _bc-realize-recursively:

The ``realize_recursively`` walker — a descriptor → operator type transformer
=============================================================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

``realize_recursively`` (``orpheus.transport.boundary_realization``, deleted in WF2.4) was the **type transformer** from the descriptor-tree algebra (``BoundaryTraceLaw | LawSum | LawScaled``) to the operator algebra (``LinearOperator``): the only path from a non-callable descriptor to a callable operator. Its successor is one method of one object per transport method, :meth:`~orpheus.transport.boundary_conditions.RealizerTable.realize`, and the walker's two jobs are split between the type and the table. The structural recursion needs no walker, because a wall is flat (:ref:`bc-rank-n-algebra`); and the lookup of the leaf realizer, which was an explicit argument of the walker, is the table's keyed dispatch.

**The realization is a functor on the response network.** For the shared family, realization is a partial map from declarations to operators on a method's trace space. On a wall it is determined by the images of the units the wall decomposes into: the image of the empty network is :math:`Z_f`, the zero map :math:`\Gamma_+(f)\to\Gamma_-(f)` of the face; the image of a chain :math:`c` with weight :math:`w` is :math:`w` times the chain realizer's operator; and the image of a sum is the sum of the images:

.. math::

   \mathcal R(0) = Z_f,
   \qquad
   \mathcal R(w\,c) = w\,\mathrm{realize}(c),
   \qquad
   \mathcal R(x + y) = \mathcal R(x) + \mathcal R(y),

summed in the canonical chain order, so the realization does not depend on how the wall was built. This is the homomorphism the walker was, with the descriptor tree's shape replaced by the network's normal form. The affine part :math:`q` realizes separately, as an array, through ``sample_inflow``, and never as an operator.

The output of the walker was already normal: the operator ``+`` and scalar ``α·`` build the formal linear combination (:ref:`composition-algebra`), so the law tree's shape was not kept. The nested law ``0.5 · (0.3·X + 0.7·Y)`` realized as the two-term combination ``0.15·X + 0.35·Y`` (`[M]` 2026-10-20, Lebedev-17 at the ``xmax`` face: terms ``(0.15, …), (0.35, …)``; the dyadic 0.5 makes the products exact). The declaration now carries that normal form, so the same wall realizes to a ``RealizedSum`` whose leaves have return fractions 0.15 and 0.35 (the first block of :ref:`bc-rank-n-algebra`, the line after the one printing ``True``).

**The amplitude.** Before wave 2 the SN realizer's ``_attenuated_kernel_operator`` returned ``α · (K_ω ⊗ I)`` at every amplitude, the bare tensor product at :math:`\alpha = 1` (the one-term collapse) and the zero map :math:`\Gamma_+ \to \Gamma_-` at :math:`\alpha = 0`, with no branch on :math:`\alpha` (`[M]` 2026-10-20 by the orchestrator, a probe over the geometry and S\ :sub:`N` operator gates: in 1619 of 1619 calls, ``0 · (K_ω ⊗ I)`` was a ``ZeroMorphism`` with exactly the spaces the dedicated vacuum builder bound). The helper moved into ``orpheus/sn/boundary/realized_operators.py`` and its amplitude is parsed into :math:`(0, 1]` (``_parse_alpha``): :math:`\alpha = 1` is still the bare tensor product, and a zero weight is no longer a leaf at all, because the wall does not store it, so the zero map is the image of the EMPTY network (``ZeroMap``, whose inner operator is the ``ZeroMorphism`` bound to the face's two half-traces).

**The door.** ``RealizerTable.realize(declared, site, space)`` is, in the source (``orpheus/transport/boundary_conditions/realization.py``, quoted and abridged by two ``assert`` lines; not a standalone script):

.. code-block:: none

   def realize(self, declared, site, space):
       self.require(declared, site, space)
       if isinstance(declared, DeckTransformation):
           return self._deck.realize(declared, pairing_of(declared, site.face, dimension=site.dimension), site, space)
       return reduce(operator.add, (self._responses[key].realize(unit, weight, site, space) for key, unit, weight in _units(declared)))

``require`` is ``admit`` followed by the raise of ``BoundaryRefused`` iff the admission is ``Refused``, so the question "is this admitted here" and the construction cannot disagree: no realizer body raises on its own for an admission reason. ``_units`` is the one decomposition that ``admit`` (the meet over the units) and ``realize`` (the sum over the images) share: one ``(key, unit, weight)`` triple per chain of a wall, or the ``ZeroNetwork`` marker when the network is empty, or the one method-only condition under its own exact type. A table is keyed by TYPE, filled by realizers that register through a class keyword (``class SNSpecularChain(ResponseRealizer, table=SN_REALIZERS, realizes=SpecularChain)``) and closed by ``seal()`` at module scope after the last registration of the method's defining module (the last line of the module for SN, diffusion and MoC; in ``orpheus/cp/solver.py`` and ``orpheus/mc/solver.py`` the ``seal()`` call precedes later code, but every registration of those two tables comes before it), so anyone holding the table has executed every registration. Lookup is by exact type and never through the method resolution order: a subclass is a new declaration and owes its own entry. The five tables are ``SN_REALIZERS``, ``DIFFUSION_REALIZERS``, ``CP_REALIZERS``, ``MOC_REALIZERS`` and ``MC_REALIZERS``.

A usage block exercising the door at its refusals and at the two ``TypeError`` boundaries (``python`` run from the worktree root with ``PYTHONPATH`` set to it; the output is pasted):

.. code-block:: python

   from orpheus.diffusion.zero_flux import ZeroFluxCondition
   from orpheus.geometry import CoordSystem
   from orpheus.geometry.deck import MirrorDeck
   from orpheus.numerics.face_layout import FaceLayout
   from orpheus.numerics.quadrature import Quadrature
   from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
   from orpheus.sn.boundary.realizer import SN_REALIZERS
   from orpheus.sn.mesh.method_space import SNMethodSpace
   from orpheus.transport.boundary_conditions import BoundaryRefused, FaceSite, WallCondition

   quad = Quadrature.lebedev(17)
   layout = FaceLayout.from_named_shapes([(f, (int(quad.n_ordinates),)) for f in ("xmin", "xmax")])
   space = SNMethodSpace(quadrature=quad, trace=AngularTraceSpace.from_quadrature_and_layout(quad, layout))
   site = FaceSite("xmax", CoordSystem.CARTESIAN, 1)

   # A right-kind declaration the method does not realize: a refusal that names its axis.
   try:
       SN_REALIZERS.realize(ZeroFluxCondition(), site, space)
   except BoundaryRefused as err:
       print(type(err).__name__, err.admission.axis.name)

   # The wrong kind of object: Python's own TypeError.
   displacement = WallCondition.white(0.3) - WallCondition.specular(0.5)
   print(type(displacement).__name__)
   for attempt in (lambda: SN_REALIZERS.realize(displacement, site, space), lambda: 0.5 * MirrorDeck()):
       try:
           attempt()
       except TypeError as err:
           print("TypeError:", str(err).split(", got")[0])

   # A declaration and a realized operator do not inter-compose.
   realized = SN_REALIZERS.realize(WallCondition.white(0.5), site, space)
   try:
       WallCondition.white(0.5) + realized
   except TypeError as err:
       print("TypeError:", err)

`[M]` 2026-10-24 at the tip ``471c4737``: the block above run as ``PYTHONPATH=<worktree root> python`` from outside the repository (``orpheus.__file__`` resolved into the worktree), its output pasted:

.. code-block:: text

   BoundaryRefused STATE_CONE
   WallDisplacement
   TypeError: sn: admit takes a DeckTransformation, a BoundaryCondition or an InflowTerms
   TypeError: unsupported operand type(s) for *: 'float' and 'MirrorDeck'
   TypeError: OperatorSum requires apply on every operand; WallCondition lacks 'apply'.

Type-system contract
--------------------

The walker's input was intentionally narrow: ``BoundaryTraceLaw`` instances and the two composer dataclasses ``LawSum`` and ``LawScaled`` were the only valid node shapes; operator-algebra composites (``OperatorProduct``, ``TensorProductOperator``, ``OperatorSum``) were not recognised, because they belong to the operator tree; and an unknown node raised ``TypeError`` and not ``BoundaryError``, because that is a **typing** failure (the caller passed the wrong kind of object) and not a BC-domain failure. The contract the table keeps has the same two-way split.

* ``admit``, ``require`` and ``realize`` take a ``DeckTransformation`` or a ``BoundaryCondition``, and ``admit`` and ``require`` also an ``InflowTerms`` (the source channel's own unit); a ``WallDisplacement``, a realized operator and any other non-declaration handed to ``admit`` or ``require`` is a ``TypeError`` whose message names the accepted types (the first ``TypeError`` of the block above; the door's message begins ``sn: admit takes a DeckTransformation``). The exclusion of an ``InflowTerms`` from ``realize`` is made by the type annotation only: ``require`` admits an ``InflowTerms``, and ``realize`` then reaches a bare ``assert isinstance(declared, BoundaryCondition)`` in ``orpheus/transport/boundary_conditions/realization.py``, so the same call raises ``AssertionError`` under plain ``python`` and, with the assertion stripped by ``python -O``, ``KeyError`` on ``self._responses[InflowTerms]`` (``[M]`` 2026-10-24, ``SN_REALIZERS.realize(WallCondition.uniform(1.0).source, site, space)`` under both invocations at the worktree tip). That path is a defect owed a ``TypeError`` raise and a gate on ``realize`` (an open item, not the contract).
* A declaration of the right kind that the method does not realize is not a typing failure: it is a ``Refused`` admission, and ``realize`` raises ``BoundaryRefused`` carrying the method, the face and the admission (the first line of the block above: ``STATE_CONE``, because the SN table's cone is non-negative flux and diffusion's zero flux reports a negative return fraction). An unregistered unit is refused on the ``PLUMBING`` axis with the issue the table names for it, or on ``STATE_CONE``, which the table decides without importing the other method's type.
* The type-level and value-level questions are one door: ``admit(declared, site, space=None)`` is the type- and site-level answer, the front end's drop-down before a problem exists, and ``catalogue()`` is its upper bound over the table's keys.

The β1 to β2 transition (:ref:`bc-rank-n-algebra`) had sharpened the walker's type signature from "any Wave-0 composer tree with ``BoundaryTraceLaw`` leaves" to "any descriptor-tree node": the dispatch table shrank from five Wave-0 composers plus the leaf to three descriptor types plus the leaf (the leaf counted as the same ``BoundaryTraceLaw`` branch). The eliminated branches (``OperatorProduct``, ``TensorProductOperator`` and ``SumOfTensorProductsOperator``, the last retired at package S2) handled operator composition patterns with no descriptor-tree analog, so they were dead dispatch paths once ``LawSum`` and ``LawScaled`` replaced the in-tree Wave-0 algebra. Removing them clarified the walker's role, which was exactly the type transformer between the two algebras. Wave 2 removed the transformer's reason to exist, the descriptor tree, and kept its one theorem as a test: ``tests/gates/sn/architecture/test_sn_realizer_table.py::test_the_realizer_is_a_homomorphism`` (law L-L) asserts, on the parameters, that the image of ``s·S + d·D`` is the ``RealizedSum`` of the two scaled leaves, that ``real(c·x)`` scales the parameters of ``real(x)``, and that the image of the empty network is the ``ZeroMap``; and, on dense matrices, bit equality for the specular image and for dyadic coefficients, and fewer than 5 ulp per entry otherwise. Its docstring records why the Marshak fixture ``(0.3, 0.7)`` is not the law's only point: the mutation that rescales the diffuse image by the wall's total return fraction is inert there and reddens on ``(0.3, 0.5)``.

Placement — the deferral that fired at the second method
--------------------------------------------------------

The walker lived in ``orpheus.transport.boundary_realization`` (the realization seam, next to the ``BoundaryRealizer`` Protocol it dispatched through), and it was **method-blind**: the leaf realizer was a REQUIRED argument, and ``method_space`` was threaded verbatim to that realizer without inspection. It was *not* on the production single-BC path: production realized one BC directly (each method-mesh's ``realize_boundary_law`` arm → ``<Method>BoundaryRealizer().realize``). The walker was the **rank-N composition entry point**: the only thing that realized a *descriptor tree* (the Marshak ``0.3 * Reflective + 0.7 * White`` partial-current BC of :eq:`affine-bc-form`) rather than a single leaf law. Production did not wire a rank-N BC then, so the walker ran only from the rank-N tests.

**How it got here is a worked example of the
defer-until-two-instances rule.** Until #290 the walker lived in
:mod:`orpheus.sn.boundary.realizer`, honestly SN-specific: it
threaded an :class:`~orpheus.sn.mesh.method_space.SNMethodSpace` and
hardcoded ``SNBoundaryRealizer`` at the leaf. (Before that it sat in
a separate ``boundary_realize`` module — the near-twin filename next
to ``boundary_realizer`` was a standing legibility hazard; merging
the two retired it.) The cross-method generalization was explicitly
**deferred until the second functional realizer ships**, with the
recorded insight that the deferral trigger was not local to boundary
realization: the same event — a second transport method arriving —
would also mint the ``TransportMethod`` Protocol flagged in
:mod:`orpheus.transport.mesh.material_mesh`, because the
**boundary-realizer seam** and the **homogenization method-layer**
were two independent witnesses to one missing type, to be typed
*together* at method #2 rather than as string-keyed half-steps.

That trigger fired when diffusion adopted the architecture (#290
P3), and the landing (#290 P7b) matched the prophecy on every point
but one: the deferral-era sketch had the walker "resolving its leaf
realizer through ``BoundaryRealizerRegistry``" — the actual carve
**dissolved the registry instead**. With a real second method in
hand, the consumption pattern was visible: production holds a
method-mesh and therefore holds its realizer; nobody resolves
realizers by name. The
:class:`~orpheus.transport.method.TransportMethod` Protocol landed
on the **method-mesh layer** (``SNProblem`` / ``DiffusionMesh`` — the
user ruling: the method-mesh IS the method's behavior carrier, not a
stateless singleton), the twin per-mesh ``_resolve_bcs`` loops
collapsed into the ONE shared
``resolve_boundary_conditions`` body,
and the walker moved to geometry's boundary package with an
explicit-realizer signature. Only
the **leaf** dispatch was ever method-bound; the **composer**
dispatch (``LawSum`` / ``LawScaled`` → operator composers)
survived the move byte-identically, exactly as recorded.

The second move was a layering correction, not a deferral firing.
The realization builds operators (:mod:`orpheus.numerics.operator`)
on a method space, and geometry is declarative: a geometry module that
imports the operator algebra at runtime is one of the edges of the
``{geometry, numerics}`` package cycle (numerics reads geometry's
rigid motions and groups). Layer-DAG unit WP-E1
(``.claude/plans/layer_dag.md`` ruling 7) therefore moved the seam
(the Protocol, the role stamp and the walker) to
``orpheus.transport.boundary_realization``, the lowest layer whose
vocabulary defines it, and the SN shim ``_BoundBoundaryOperator`` to
``orpheus.sn.boundary.bound_compat``, beside its one production
consumer ``SNProblem``. No shim or re-export was left in geometry, and
``tests/gates/test_layer_imports.py`` refuses a runtime import of
:mod:`orpheus.numerics.operator` by any geometry module.

**Wave 2 (2026-10-24) closed the thread.** The seam retired in WF2.4: the Protocol and the walker with the descriptor tree they served, the role stamp as ``stamp_boundary_role`` in ``orpheus/transport/boundary_conditions/realization.py`` (its one definition), the SN shim ``_BoundBoundaryOperator`` as the typed ``RealizedBoundary``, and ``resolve_boundary_conditions`` as ``realize_boundary`` (``orpheus/transport/method.py``), the one shared body that realizes each traced face of a method-mesh's binding (and requires the untraced ones, such as the inner face of a hollow body, at the site level only) through the method's ``realize_boundary_role`` and so through ``RealizerTable.realize``. The "rank-N composition entry point" is no longer a separate path: a mixed wall is an ordinary ``WallCondition``, and ``SNProblem.realize_boundary_role`` hands it to ``SN_REALIZERS.realize`` like any other (read from ``orpheus/sn/problem.py``; the executed form is the first block of :ref:`bc-rank-n-algebra`). The ``TransportMethod`` Protocol survives at the method-mesh layer, where #290 P7b put it; its table member is ``boundary_realizers``.


.. _bc-vacuum-semantic-correction:

The vacuum semantic correction (§16A.5) — and its dissolution at B3.2
=====================================================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`. Its subject, the vacuum law, is the empty response network at the tip: its image is the ``ZeroMap`` :math:`\Gamma_+(f) \to \Gamma_-(f)`.

.. important::

   **Status (campaign phase B3.2, 2026-07-31).** This whole section is
   now **history**. Vacuum realizes to the **zero map**
   :math:`\Gamma_+ \to \Gamma_-`
   (:ref:`bc-domain-narrowing`); the inflow-only mask, the outflow rows
   it preserved, and the question of what preserving them was *for* are
   all off the vacuum path. The section is kept because the reasoning
   below is the reasoning a future session would otherwise re-derive —
   and because the way it was resolved (by narrowing the domain rather
   than by answering the question) is the more general lesson.

The most subtle design decision of the trace-law refactor concerned
vacuum. Pre-Wave-7 the legacy ``VacuumBoundaryOperator.apply`` body
was:

.. code-block:: python

   def apply(self, psi_out: np.ndarray, quadrature) -> np.ndarray:
       return np.zeros_like(psi_out)

This returns **all zeros**, including the outflow ordinates that
the BC has no physical interpretation for (vacuum says nothing
about :math:`\gamma_+ \psi`; it only sets :math:`\gamma_- \psi = 0`).

The post-Wave-8 SN realizer's vacuum branch therefore returned an
``IncomingOrdinateMaskTensor`` — a
sparse mask that zeroed **only the inflow ordinates** and preserved the
outflow trace. That was read as the §16A.10 *trace-correct*
representation: a projector onto the outflow ordinate subspace, which
looked like the right algebraic object for the affine law
:eq:`affine-bc-form` to read.

Why it seemed to matter — and why each reason dissolved
--------------------------------------------------------

.. note::

   **Retraction (2026-07-31, campaign phase B3.2).** The three
   consequences below were the published justification for preserving
   the outflow rows. **All three are now moot, and the first was
   measurably false.** They are preserved verbatim in substance,
   past-tensed, with each disposition stated — deleting them would
   destroy the record of *why* a wrong contract looked right for two
   campaign phases.

Three downstream consequences were argued to make the inflow-only mask
the right contract:

1. **Sensitivity adjoints.** *The argument was:* a future adjoint
   sensitivity path needs the outflow trace preserved to compute the
   response of an outgoing-current functional to the inflow BC; the
   zeros-all body would silently lose the gradient at the boundary.

   **Disposition: the consumer was measured not to exist**, and B3.2
   removed the preservation entirely. The rows are outside the
   operator's domain now, so there is nothing to preserve. This was the
   load-bearing justification and it was a *promise about a future
   consumer* — precisely the "declared capability, no consumer" pattern
   the boundary review catalogued five times in this subsystem. The
   adjoint work that did land (#276, :ref:`g-adjoint`) reaches the
   boundary through :math:`B^{\mathsf T}` and the metric-correct
   ``B.H``, neither of which needs a vacuum law to echo its input.
2. **Compositional clarity.** *The argument was:* the realized vacuum
   operator is self-adjoint and idempotent — a projector — whereas the
   zeros-all operator is the ``ZeroOperator`` projector, the wrong type
   tag for "inflow-mask only".

   **Disposition: inverted.** With the domain narrowed to
   :math:`\Gamma_+`, the vacuum law is not an endomorphism at all, so
   "idempotent" is not even a well-typed thing to ask of it — the
   composite :math:`\text{law} \circ \text{law}` does not exist.
   :class:`~orpheus.numerics.operator.ZeroOperator` was then exactly the
   right type tag, carrying two space hooks rather than an echo. At the tip the
   image is ``ZeroMap``, whose inner operator is a
   :class:`~orpheus.numerics.operator.ZeroMorphism` bound to the face's two
   half-traces, so it emits the zero of the space it lands in and not an echo
   of the one it came from.
3. **Algebraic uniformity.** *The argument was:* every other rank-1 law
   acts on the inflow ordinates only and leaves the outflow rows
   untouched, so the legacy vacuum was the asymmetric special case.

   **Disposition: the uniformity is real and B3.2 made it structural
   rather than conventional.** Every law now *maps*
   :math:`\Gamma_+ \to \Gamma_-`; "leaves the outflow rows untouched"
   is no longer a behaviour a law could get wrong, because those rows
   are not in its codomain. The observation was right; the mechanism
   was one layer too shallow.

The algebra: where the two legacy semantics agreed and where they diverged
--------------------------------------------------------------------------

Decompose the angular ordinate set at a given face :math:`f` as

.. math::
   :label: ordinate-partition-inflow-outflow

   \{1, \ldots, N\} \;=\; I_f \,\sqcup\, O_f \,\sqcup\, T_f,

with :math:`I_f = \{n : \Omega_n \cdot \hat n_f < -\epsilon\}` the
inflow set, :math:`O_f = \{n : \Omega_n \cdot \hat n_f > +\epsilon\}`
the outflow set, and :math:`T_f` the (measure-zero in the continuum,
:math:`\epsilon`-band in the discretisation) tangential set. For
:math:`\psi_{\text{out}} \in \mathbb{R}^N` representing the trace
ordinate values at face :math:`f`, the two vacuum representations
are:

.. math::
   :label: vacuum-legacy-vs-trace-correct

   \mathrm{zeros\_all}(\psi_{\text{out}})[n] &= 0
       \qquad \forall\, n \in \{1, \ldots, N\}, \\[2pt]
   \mathrm{inflow\_mask}(\psi_{\text{out}})[n] &=
       \begin{cases}
           0 & n \in I_f, \\[2pt]
           \psi_{\text{out}}[n] & n \in O_f \cup T_f.
       \end{cases}

.. (vv-status rationale) HISTORICAL structural/explanatory identity: contrasts
   the two pre-B3.2 vacuum spellings (zeros-all vs the inflow-only mask) —
   they agree on the inflow set and diverge on the outflow set. NEITHER is the
   live realization: since campaign phase B3.2 vacuum realizes to the zero map
   Γ₊ → Γ₋, and the live gate is the re-posed vacuum snapshot case in
   ``tests/gates/sn/architecture/test_wave2_carry.py`` (the retired
   ``test_bc_equivalence_snapshot.py`` asserted Γ₋ shape + all zero). Kept as the record of a contract that looked right for two campaign
   phases. An explanatory comparison, never a solver claim.
.. vv-status: vacuum-legacy-vs-trace-correct documented

The two functions agree **on** :math:`I_f` (both give 0) and
**diverge** on :math:`O_f` (legacy gives 0, trace-correct gives
:math:`\psi_{\text{out}}[n]`). They diverge on :math:`T_f` too.

.. warning::

   **The published mitigation for** :math:`T_f` **was an
   initialisation claim, not a structural guarantee — and the two
   spellings of "outflow" genuinely differ.** This page previously
   argued that ORPHEUS's quadrature adapters carry every tangential
   ordinate at :math:`\mu = 0`, so :math:`\psi_{\text{out}}[n] = 0` on
   :math:`T_f` "for a properly-initialised flux", making the divergence
   physically restricted to :math:`O_f`. Two corrections, both
   **measured** during the B3 crosswalk:

   1. :math:`\mu_{\rm axis} \approx 0` is what *defines* a tangential
      ordinate; it does not make :math:`\psi` vanish there. What is
      true is weaker and conditional — the boundary operators write
      **only** the inflow and outflow slots, so a zero-initialised
      trace carrier keeps its tangential rows at zero. That is a
      property of the *carrier's initialisation*, not of the algebra.
      (The tangential rows do carry zero *weight* in the
      :math:`|\Omega\cdot\hat n|\,w` trace metric, which is a separate
      and genuinely structural fact — see :ref:`g-adjoint`.)
   2. The mask preserved :math:`O_f \cup T_f` while
      :meth:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.outflow_indices_for_face`
      selects **strict** :math:`O_f`. So
      :math:`I - P_{\rm in} \neq P_{\rm out}`, and the mask was never a
      projector onto :math:`\Gamma_+` — it is
      :math:`I - \iota_-\gamma_-`, a different map whenever the
      quadrature carries tangential ordinates, which is **every**
      production quadrature except ``gauss_legendre(4)``. Measured on a
      cylinder under ``product(n_mu=2, n_phi=4)``: **4 of 8** ordinates
      at ``xmax`` are tangential.

   Both corrections are moot for vacuum since B3.2 — the map is the
   zero map on :math:`\Gamma_+` — but the second is *live discipline*
   everywhere else: **never spell "not inflow" as "outflow"**.

The §16A.5 production-relevant subset was **the inflow rows**.
Every SN sweep call site read :math:`\psi_{\text{in}}[n]` only for
:math:`n \in I_f` — outflow rows were never consumed downstream.
The Wave 8 close-out audited all 13 ``bc.apply(...)`` sites in
``orpheus.sn.loss_representation`` (the dissolved ``sweep.py``) and :mod:`orpheus.sn.operators.streaming`
(the file:line references below are frozen at that audit and do not
resolve against the live tree — Wave O later removed ``bc.apply`` from
the sweep entirely, see :ref:`bc-extraction`):

* ``sweep.py:334,351`` (1-D slab) read
  ``psi_face_left_in[n_half + n]`` for positive-μ ordinates
  only.
* ``sweep.py:508,654`` (spherical / cylindrical) gate the read
  by ``mu_n < 0`` / ``eta_n < 0`` (inflow at the outer face).
* ``sweep.py:843-854`` (2-D wavefront) reads
  ``full_face_x[oct_idx]`` only, where ``oct_idx ↔ inflow
  ordinates`` by construction.
* ``operator.py:230-256`` (2-D FD matvec) has explicit
  ``mu_x[n] > 1e-15`` / ``mu_y[n] < -1e-15`` gates per face.
* ``operator.py:530`` (spherical FD matvec) gates on
  ``quad.mu_x[n] < -1e-15``.

No call site read :math:`\psi_{\text{in}}[n]` for
:math:`n \in O_f`, so the two semantics produced **bit-identical
observable output** for every production consumer of that era.
That is what made the §16A.5 change a safe semantic upgrade — and,
read a second time, it is exactly the evidence that the preserved
outflow rows had **no consumer at all**. B3.2 acted on the second
reading.

Post Issue #188 + #176 (2026-05-11) the **realizer path is
uniform** across every supported mesh — 1-D Cartesian, 1-D
spherical, 1-D cylindrical, and 2-D Cartesian; that uniformity
survives B3.2 unchanged (every geometry's vacuum face realizes to the
same zero map, every geometry's reflective face to a permutation on
its own reduced axis). **Empirical confirmation at the time**:
spherical 26/26 + cylindrical 25/25 + MMS curvilinear 2/2 xfail
(pre-existing ERR-026) green on C188.3 — the curvilinear sweeps had
been consuming the zeros-all legacy body via the bound-quadrature shim
and moved onto the realizer output, bit-identical on inflow rows (the
only rows the sweep read), confirming the Wave 8 call-site audit
empirically.

The Wave 6 snapshot harness gated the vacuum case explicitly, and its
assertion was **re-posed, not weakened, at B3.2**: the
``vacuum_lebedev17`` case cross-checked the live ``inflow_indices``
against the frozen snapshot's index set, fed the realized law
``psi_out[outflow_indices]``, and asserted *both* that the emission had
:math:`\Gamma_-` shape and that it was identically zero — "the narrowed
vacuum law is the ZERO map :math:`\Gamma_+ \to \Gamma_-`". The old
inflow-rows-only comparison, which documented the intentional semantic
divergence on the outflow rows, had nothing left to compare. The two index fields
retired in WF2.4, and the case is one of the seven of
``tests/gates/sn/architecture/test_wave2_carry.py::test_the_equivalence_snapshots_carry``:
``WallCondition.vacuum()`` realized through ``SN_REALIZERS`` on ``xmin`` of a
Lebedev-17 trace, its image ``array_equal`` to the image the HEAD realizer gave and to
the frozen reference ``psi_in``, so the claim is unchanged.

"Option a" (Wave 7) — historical context, retired Issue #186
-------------------------------------------------------------

The Wave 7 brief considered three migration strategies for the
2-arg legacy path. **Option (a)** ("vacuum-stays-legacy") landed:
``VacuumInflow`` carried a standalone
``apply(psi_out, quad)`` whose body was
``np.zeros_like(psi_out)`` (the pre-§16A.5 zeros-all form), and
the realizer path produced the inflow-only-mask form via
``IncomingOrdinateMaskTensor``.
The two paths agreed on inflow rows (the production-relevant
subset) and diverged on outflow rows.

Options (b) ("face-aware BC" — add ``face`` to the constructor)
and (c) ("combined ABC merge") were rejected because they would
have distorted the law's interface to better serve a transitional
path the refactor was retiring anyway.

**Status after Issue #186 (2026-05-11):** Option (a)'s standalone
``apply`` body was **deleted**. ``VacuumInflow`` (like
every other concrete BC) was a pure descriptor with no ``apply``
method; the only path to vacuum action was
``realize_recursively`` or ``SNBoundaryRealizer``. The
**two-paths-divergence was therefore eliminated by design** — there
was no longer a "second path" that could disagree with the realizer
path. What that single path *produces* has since changed again: the
§16A.5 inflow-only mask was the unique vacuum semantics from #186
until campaign phase **B3.2**, which replaced it with the zero map
:math:`\Gamma_+ \to \Gamma_-` (:ref:`bc-vacuum-semantic-correction`).
The uniqueness claim is what survived; the operator behind it did not. It holds at the tip for the same structural reason: the only path to vacuum action is ``RealizerTable.realize``, and its image is the ``ZeroMap``.

This is the load-bearing architectural payoff of B3 + β2: the
documentation no longer needs to caveat which path you're on,
because there is only one path. See
:ref:`bc-trace-law-descriptor-model` for the design rationale.


.. _bc-trace-law-descriptor-model:

The trace-law descriptor model
==============================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is :ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Issue #186 / Scope B3 + β2 (landed 2026-05-11 on branch ``feature/bc-curvilinear-realizer-cleanup``) was the architectural **closure** of the BC trace-law refactor, and wave 2 retired the architecture it closed. B3 + β2 collapsed the remaining 2-arg ``apply`` affordance of the Wave-8/9 era into a **pure-descriptor** contract:

* ``BoundaryTraceLaw`` stopped inheriting :class:`~orpheus.numerics.operator.LinearOperator`.
* The abstract ``apply`` the mixin used to provide was removed, and so was ``apply_transpose``, and so was any operator-surface advertisement: the two-axis ``is_invertible`` / ``is_adjointable`` predicates (before the #226 carve P4, the retired ``capabilities: ClassVar[frozenset[str]]`` frozenset).
* Every concrete BC (vacuum, reflective, white, albedo, periodic, prescribed-inflow, zero-flux) became a **frozen dataclass** carrying only its parameters (axis, albedo, source, ...) and the relevant ``assert_*`` invariants, with **no** ``apply`` method on any concrete BC. All three affine-factor properties were overridden on every concrete law since campaign phase B1 (the measured table is under :ref:`bc-law-layer`); at the time of Issue #186 only ``source`` was (by prescribed-inflow). ``kind`` was derived once on the ABC from the registry key, and no concrete law overrode it (7 of 7 registered laws, `[M]` 2026-10-01).
* The base class carried a **minimal algebra** that returned ``LawSum`` and ``LawScaled`` nodes, the descriptor-tree composition algebra of :ref:`bc-rank-n-algebra`: the dunders ``+``, ``-``, ``*``, ``/``, unary ``-`` and their reflected variants, each returning a new descriptor-tree node and none an operator.
* ``SNBoundaryRealizer`` was the **sole** bridge from descriptor to callable. Calling ``law.apply(psi)`` raised ``AttributeError`` at runtime, and a static type checker flagged the call without running it.

The §16A.3 three-layer architecture (descriptor / realizer / operator) was thereby **enforced by the type system**, not by convention.

**What replaced it.** The three layers stand and each is a different object. The descriptor is a frozen ``ContentIdentity`` value, ``WallCondition`` (or a deck), with no ``apply``: ``hasattr(marshak, "apply")`` printed ``False`` in the first block of :ref:`bc-rank-n-algebra`, `[M]` by that block, and no gate asserts it of a wall (``git grep -n 'hasattr([^)]*apply' -- tests orpheus`` finds such assertions only on the functionals and integral kernels of ``tests/gates/transport/test_functional_category.py`` and ``test_integral_kernel_category.py`` and a private-method check in ``tests/gates/mesh/test_mesh1d.py``, none on a wall). The realizer is the method's sealed ``RealizerTable``, keyed by unit type, and the operator is a typed ``RealizedBoundary`` (``ZeroMap``, ``ScaledOrdinatePermutation``, ``RankOneReemission``, ``RealizedSum``) whose facts (its return fraction, whether it permutes ordinates, whose outflow it consumes) are properties of its type and not of a stored flag. What the enforcement mechanism lost is the registry: a law was named by a key string, and a declaration is now an object whose class is the key.

What was tried and rejected before B3 + β2 landed
-------------------------------------------------

This subsection preserves the design history because the rejected
paths are the load-bearing intellectual content of the close-out
— future sessions asking "why does the realizer exist?" need to
see why every alternative failed.

**Option A** (Issue #176 / C176.3, ~2026-04). Concrete BC
``apply`` methods kept a keyword-optional
``quadrature: AngularQuadrature | None = None`` parameter with
defensive errors. Every ``apply`` named in this subsection is written
as a literal, not a role: the descriptor classes carry no ``apply``
today — that is the whole point of the bullet above — so a live
``:meth:`` role here would advertise a link to a method the type
system deliberately removed.

* ``ReflectiveBoundary.apply`` / ``WhiteBoundary.apply``
  raised ``BoundaryError`` when ``quadrature is None``
  because their geometric / response operators needed the
  quadrature to construct themselves.
* ``VacuumInflow`` / ``AlbedoBoundary`` /
  ``PeriodicBoundary`` / ``PrescribedInflow`` accepted
  and ignored the ``quadrature`` parameter.

Option A was the **interim** landing — it preserved the
direct-call convenience pattern ("sketching code can write
``ReflectiveBoundary(axis='x').apply(psi, quad)``") while routing
production through the realizer. The C176.3 audit identified
three architectural costs that made Option A unsustainable:

1. **Asymmetric semantics on ``quadrature=None``.** Two BCs
   raised; four accepted-and-ignored. The behaviour was
   inconsistent across the law family and required per-BC
   documentation of "when is this method usable".
2. **Vacuum two-paths-divergence.** Direct ``VacuumInflow.apply(psi)``
   returned ``np.zeros_like(psi)`` (the pre-§16A.5 zeros-all body).
   The realizer-routed path returned
   ``IncomingOrdinateMaskTensor``
   output (the §16A.5 inflow-only mask). The two paths agreed on
   inflow rows (the production-relevant subset; see
   :ref:`bc-vacuum-semantic-correction`) but diverged on outflow
   rows. The divergence was harmless at every existing production
   call site; the cost argued at the time was a documentation-burden
   landmine for *future adjoint-sensitivity consumers that read
   outflow rows*. **That consumer was later measured not to exist**,
   and B3.2 removed the outflow rows from the law's domain
   altogether — so the cost was real (two paths, one contract) but the
   reason given for it was not. The retirement of Option A stands on
   costs 1 and 3.
3. **Liskov violation.** The abstract
   ``BoundaryTraceLaw.apply(self, psi_out)`` was strict 1-arg
   (post Issue #176 / C176.4). The concrete
   ``apply(self, psi_out, quadrature=None)`` was technically
   Liskov-substitutable (the optional parameter has a default), but
   calling ``bc.apply(psi)`` polymorphically on a
   ``BoundaryTraceLaw``-typed parameter could fail at runtime
   for Reflective/White — the type signature said "this works" but
   the runtime behaviour said "this raises". The static type
   system could not catch the failure because the contract was
   carried only in the docstring.

The B3 audit cataloged every remaining 2-arg ``bc.apply(psi,
quad)`` call site in production AND tests; none was load-bearing
for correctness. The Wave-6 snapshot regression carried legacy
halves (regenerated through the realizer path), the
``PrescribedInflow`` invariant tests ignored the
``quadrature`` argument anyway, and the realizer-vs-legacy
equivalence assertions could be replaced by hand-computed
expressions strictly stronger than legacy-agreement. The B3 sweep
rewrote every such site in one commit cycle (C-B3.7 through
C-B3.12) and deleted the ``apply`` methods.

**β1 interim** (sub-option within Issue #186 B3, considered but
not landed as a final state). Keep
:class:`~orpheus.numerics.operator.LinearOperator` inheritance on
``BoundaryTraceLaw`` and drop only the abstract
``apply``. Rank-N composition would then build an
:class:`~orpheus.numerics.operator.OperatorSum` of ``ScaledOperator``-wrapped *unrealized*
laws, and ``realize_recursively`` would walk the resulting
Wave-0 composer tree. This β1 form was **algebraically equivalent**
to β2 but conflated two type families: the
:class:`~orpheus.numerics.operator.OperatorSum` instance representing a not-yet-realized
descriptor tree was structurally identical to the
:class:`~orpheus.numerics.operator.OperatorSum` representing an actual operator composition,
and only runtime inspection of the leaves could tell which was
which. β2 was preferred because the type system could then enforce
"this is a law, that is an operator" by static inspection — a
type checker rejects ``law_tree.apply(...)`` at the linter level
without ever running the code. β2 is more verbose (two new
dataclasses) but is the architecturally-checkable form.

**The vacuum two-paths-divergence was eliminated by design.** Under
B3 + β2 there was no "direct path" any more — every consumer must
realize the law before applying it, so there is exactly one vacuum
semantics at any moment (the inflow-only mask until campaign phase
B3.2, the zero map :math:`\Gamma_+ \to \Gamma_-` since). The Wave-6
snapshot harness's ``vacuum_lebedev17`` case tracked whichever it was —
at the tip it is one of the seven cases of
``tests/gates/sn/architecture/test_wave2_carry.py``, carrying the ``ZeroMap`` image — and
no caller can accidentally invoke the pre-§16A.5 zeros-all path
because that path no longer exists in the code.

Empirical justification
-----------------------

The 18-test ``tests/gates/geometry/test_law_composition.py`` suite pinned the descriptor-tree contract (foundation and L1 tests):

* algebra closure on every dunder for every node-type pairing (laws × ``LawSum`` × ``LawScaled``);
* ``LawScaled`` constant folding (``2 * (3 * spec) == LawScaled(6, spec)``);
* the walker's exhaustive dispatch on the three node types and its ``TypeError`` on unknown nodes;
* walker value-correctness against hand-composed expectation (``realize_recursively(law_tree, ms, realizer).apply(psi) == 0.3 * realize(spec).apply(psi) + 0.7 * realize(white).apply(psi)``);
* the absence of ``apply`` on any descriptor-tree node (``not hasattr(tree, "apply")``).

The suite retired with the tree in WF2.4 (commit ``ece6a281``). Each claim has a successor, or became unspellable:

* closure of the arithmetic and the cone: ``tests/gates/transport/test_wall_condition_algebra.py::test_the_algebra_and_the_cone`` (gate 4), with the source summed through every combination (gate 5);
* constant folding: no nesting exists to fold, and the identity is a value equality of two walls (gate 4);
* exhaustive dispatch and the typing ``TypeError``, three separate claims. The binding refusal (a displacement is a value and not a condition, so ``BoundaryConditions({1: displacement})`` is refused) is ``test_displacements_are_values_not_conditions`` of the same file (gate 6). The operator-protocol refusal (a condition is not a deck slot, a deck is not an operand), which must come from Python's operator protocol with no production frame in the traceback, is ``test_a_condition_is_not_a_deck_slot`` and ``test_a_deck_is_not_an_operand`` (gate 7). The door's own ``TypeError`` (``sn: admit takes a DeckTransformation ...``), raised from a production frame (``realization.py``) and so outside gate 7's predicate, is witnessed on ``admit`` only, by ``tests/gates/transport/test_admission_mechanism.py::test_the_source_door_is_met_for_every_condition``; no test asserts that ``realize`` of a displacement raises ``TypeError`` (unwitnessed);
* walker value-correctness: ``test_the_realizer_is_a_homomorphism`` (:ref:`bc-realize-recursively`), whose reference is the parameters of the leaves and a dense matrix of each image, a stronger statement than agreement with a hand-composed sum of two images that share the realizer;
* absence of ``apply``: no successor gate (see above).

The reference-image harness verifies that the realized output matches the independently derived image on every case, bit-identically for vacuum, the two specular rows and the periodic row, and inside the reduction-order bound for the Lambertian and Marshak rows (the operator tree is structurally the same; only the route to it changed). Until the 2026-08-01 re-anchoring the comparison was against the pre-B3 realizer-path *output*; :ref:`bc-snapshot-reanchoring` records why that was a drift lock rather than a reference, and what replaced it.

Call-site contract
------------------

There is **one** way to realize a declaration, and it is the table of the method (the usage block of :ref:`bc-realize-recursively` is the executed form):

.. code-block:: python

   # ``site``, ``space`` and ``psi_out`` as in the blocks above.
   op = SN_REALIZERS.realize(WallCondition.specular(1.0), site, space)   # a RealizedBoundary on this face
   psi_in = op.apply(psi_out)                                            # Gamma_+ shaped argument, Gamma_- shaped result

A face needs a **site** and a **method space**: since B3.2 a realized face's domain is that face's :math:`\Gamma_+`, which a quadrature alone cannot name, and the site names the face (the face is never a field of the space). The Cartesian mirror is not a wall at all but a deck, ``MirrorDeck()``, realized by the same table through its deck realizer: ``SN_REALIZERS.realize(MirrorDeck(), site, space)`` gives the same image as a unit specular wall, because a symmetry has no amplitude. For a mixed wall the call is the same, and the tree is the wall itself:

.. code-block:: python

   tree = 0.3 * WallCondition.specular(1.0) + 0.7 * WallCondition.white(1.0)
   op_tree = SN_REALIZERS.realize(tree, site, space)
   psi_in = op_tree.apply(psi_out)

No call site routes through a putative ``law.apply(psi)``: that method does not exist on a wall either. The retired spelling, for the reader who meets it in the history above, was ``SNBoundaryRealizer().realize(law, SNMethodSpace.for_face(...))`` followed by ``realize_recursively(tree, ms, SNBoundaryRealizer())`` for a tree; ``SNMethodSpace.for_face`` and the ``minimal`` constructor retired with the index sets, and ``SNMethodSpace`` is built from a quadrature and a trace.


.. _bc-numerical-evidence:

Numerical evidence
==================

The reference-image harness
(:mod:`tests.gates.sn.architecture.test_wave2_carry`, which took the seven frozen images over from the retired
``tests/gates/geometry/test_bc_equivalence_snapshot.py`` at wave 2 slice S2b) is the widest
mutation net in the boundary subsystem. Seven cases compare the
realized operator against a frozen ``.npz`` image at a per-case
tolerance.

.. _bc-snapshot-reanchoring:

The 2026-08-01 re-anchoring — from recording to reference
---------------------------------------------------------

Until campaign phase **B3.4b** these artefacts were **recordings of
production output**: the generator called
``SNBoundaryRealizer().realize(...)`` and froze ``op.apply(psi_out)``.
Two consequences followed mechanically. The assertion was
``production == a recording of production`` — a regression LOCK, worth
exactly what the recorded code was right, and unable to say the value
is *correct* (``vv-principles`` §bit-identity criterion 2). And it
broke on every signature change, which it did twice: at **B3.2**, when
the SN law narrowed to :math:`\Gamma_+ \to \Gamma_-`, and at
**B3.4a**.

The interim repair kept the frozen full-face artefacts and *restricted*
them, arguing that an artefact frozen before the narrowing is an
independent statement. That is only **procedural** independence: it
certifies *"the new path agrees with the old path"*, and the premise of
the whole narrowing campaign is that the old path read the wrong
half-trace.

So every case was **re-anchored against an expression derived from the
mathematics** — :eq:`affine-bc-form`, the mirror isometry
:math:`\Omega \mapsto \Omega - 2(\Omega\cdot\hat n)\hat n`, and the
Lambertian partial-current balance :math:`J^- = \alpha J^+`. The
generator imported nothing from
:mod:`orpheus.sn.boundary.realizer`, which made it refactor-proof by
construction, and the artefact **IS** the reference rather than a
recording (WF2.4 deleted the generator; see the file paragraph below).

The file stays frozen for one reason: a reference the harness
recomputed would let the generator's expression and production drift
*together*. The committed artefact is the barrier — an expression
change moves the gate only through a regeneration, which is a
reviewable diff; since wave 2 the carry fixture also records each file's
sha256, and ``test_the_snapshots_were_not_rebaselined`` in
:mod:`tests.gates.sn.architecture.test_wave2_carry` reddens on a file that
no longer matches it.

**What the migration measured.** Against the retired pre-B3.2
recordings restricted to :math:`\Gamma_-`: ``vacuum``, ``specular_x``
and ``specular_y`` are **bit-identical**; ``white_xmax``,
``white_xmin`` and ``mixed`` differ by **1–2 ULP** (reduction order
only — the reference contracts with :func:`numpy.tensordot` where
production runs a broadcast-multiply-then-``sum``); ``periodic``
differs by **98 %**, by design. Six of seven claims did not move; what
changed is where they come from.

.. list-table:: BC reference-image cases
   :header-rows: 1
   :widths: 24 26 16 34

   * - Case
     - Declaration today (the one the recording was made under, in parentheses)
     - Quadrature / face
     - Reference, and its tolerance
   * - ``vacuum_lebedev17``
     - ``WallCondition.vacuum()`` (``VacuumInflow()``)
     - Lebedev 17, ``xmin``
     - :math:`R = q = 0`, so the **zero map**: an image with
       :math:`\Gamma_-` shape, identically zero.
       ``assert_array_equal`` — no arithmetic is performed.
   * - ``specular_x_lebedev17``
     - ``WallCondition.specular(1.0)`` (``ReflectiveBoundary(axis="x")``); the same image is also asserted for a ``MirrorDeck()`` on the Cartesian face
     - Lebedev 17, ``xmax``
     - the mirror gather :math:`\psi^-(\Omega) = \psi^+(\Omega')`.
       ``assert_array_equal`` — reduction depth 0.
   * - ``specular_y_partial_07_LS6``
     - ``WallCondition.specular(0.7)`` (``AlbedoBoundary(0.7, SpecularReturn(axis="y"))``)
     - LevelSymmetricSN(6), ``ymax``
     - the same gather scaled by α — the α-fold row, and the only
       one on a non-``x`` axis. ``assert_array_equal``. The partially
       specular wall is a response, so this row exercises the wall's
       specular chain: ``SNSpecularChain`` realizes it as α times the
       ordinate permutation of the mirror, a
       :class:`~orpheus.sn.boundary.realized_operators.ScaledOrdinatePermutation`
       (the recording shared the deck kernel of the then-law); its
       frozen image is the one recorded when the row was spelled as
       the then-attenuated mirror, bit-identical across both
       respellings (:ref:`bc-history`).
   * - ``white_xmax_LS4``
     - ``WallCondition.white(1.0)`` (``WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)``)
     - LevelSymmetricSN(4), ``xmax``
     - isotropic re-emission with :math:`J^- = J^+` in the
       :math:`w_n|\Omega_n\cdot\hat n|` measure.
       ``rtol = |Γ₊|·ε``.
   * - ``white_xmin_partial_03_GL``
     - ``WallCondition.white(0.3)`` (``WhiteBoundary(axis="x", outward_sign=-1, albedo=0.3)``)
     - GaussLegendre1D(8), ``xmin``
     - the same law at α = 0.3 on a quadrature whose :math:`\sum w`
       is 2, not :math:`4\pi` — the canary against a hard-coded
       normalisation. ``rtol = |Γ₊|·ε``.
   * - ``mixed_30spec_70white_LS4``
     - ``WallCondition.from_amplitudes(0.3, 0.7)`` (``0.3 * AlbedoBoundary(1.0, SpecularReturn("x")) + 0.7 * white``,
       the Wave-0 algebra over the realized leaves; the sum is now the wall's own, realized to a ``RealizedSum``)
     - LevelSymmetricSN(4), ``xmax``
     - the pointwise convex combination of the two images above.
       ``rtol = (|Γ₊| + 2)·ε``.
   * - ``periodic_lebedev17``
     - ``TranslationDeck()`` (``PeriodicBoundary()``)
     - Lebedev 17, ``xmin`` ← ``xmax``
     - the PARTNER face's outflow at the same ordinates (the harness reads the stored
       ``psi_out_partner`` when the realized face's ``domain_face`` is not the installation face).
       ``assert_array_equal``; the row ran inside ``xfail(strict=True)`` until the partner-face
       wrap landed, and passes plainly today.

Every tolerance above is **derived, not measured-and-rounded-up**. The
old ``nulp = 4`` / ``nulp = 64`` constants encoded no claim; the bound
now states one. Recursive summation of :math:`n` terms carries a
relative error :math:`\le (n-1)\,u\,\kappa` with
:math:`u = \varepsilon/2`; every summand in the Lambertian is positive
(the probe is :math:`U(0,2)`, and :math:`w_n|\Omega_n\cdot\hat n| > 0`
on :math:`\Gamma_+`), so :math:`\kappa = 1` and :math:`|\Gamma_+|\cdot
\varepsilon` is that bound with room for the trailing division. The
mixed row adds one rounding for its scaling and one for the sum.
Measured against those bounds: 1.9e-16 (2.7e-15), 1.2e-16 (8.9e-16),
3.2e-16 (3.1e-15).

.. note::

   **The periodic row was a deliberate strict xfail, not a passing
   gate, until the partner-face wrap landed.** A translation identifies the two faces without touching
   direction, so a particle leaving through ``xmax`` in direction
   :math:`\Omega` re-enters at ``xmin`` in the SAME :math:`\Omega`:
   :math:`\Gamma_-(\text{xmin})` and :math:`\Gamma_+(\text{xmax})` are
   the same ordinate SET, and the law is a **two-face coupling**.
   Production then realized it as a per-face angular identity, so the
   narrowed composition
   :math:`\iota_- \circ \text{law} \circ \gamma_+` hands it this
   face's own :math:`\Gamma_+ = \{\mu_x < 0\}` and it feeds
   outgoing-left flux back in as incoming-right flux — measured 98 %
   relative error, an :math:`\mathcal{O}(1)` wrong answer.

   The artefact therefore carries TWO probes, drawn from independent
   seeds: with one shared draw a per-face endomorphism would look
   correct on the rows that coincide. **B3.4c** built the partner-face
   wrap; the marker's deletion was forced by the XPASS(strict) failure,
   which is this campaign's standing technique — the xfail set IS the
   todo list. The body itself is already correct: a companion test
   shows the identity reproduces the reference exactly when fed the
   PARTNER's :math:`\Gamma_+`, so the defect was entirely in which
   half-trace the composition supplied. Today a ``TranslationDeck`` realizes to a face whose
   :attr:`~orpheus.sn.boundary.realized_operators.RealizedBoundary.domain_face` is the PARTNER
   (``pairing_of`` gives the pairing :math:`(f, f')` with :math:`f'` the opposite face), and the row
   reads that partner's outflow.

.. note::

   **A narrowed row is re-posed onto** :math:`\Gamma_+`\ **, never
   weakened.** The B3.2 / B3.4a recipe was to feed the realized law
   ``snapshot["psi_out"][space.outflow_indices]`` and compare against
   the frozen pre-narrowing image restricted to :math:`\Gamma_-`. The
   2026-08-01 re-anchoring made that restriction unnecessary — the
   artefacts are now stored on the narrowed spaces directly, with the
   face and BOTH index sets alongside — and it kept the probe: each
   stored ``psi_out`` is the same deterministic full-face draw the
   retired artefacts used, restricted to :math:`\Gamma_+`. The
   migration changed the schema and the reference; it never changed
   what the operators are probed with.

   **The Marshak row stopped being an honest red at B3.4a.** While
   white was un-narrowed, ``0.3 * spec + 0.7 * white`` mixed a narrowed
   leaf with a full-\ :math:`N` one, so the sum's factors disagreed and
   its ``apply`` raised on a :math:`\Gamma_+` probe —
   ``AngularAverageOperator.apply: psi.shape[0] = |Γ₊|, expected N``.
   The row was pinned by ``xfail(strict=True)`` precisely so that the
   narrowing would force the marker's deletion. Both leaves are now
   narrowed and the composition is measured correct: it returns
   :math:`|\Gamma_-|` rows equal to
   ``0.3 * spec.apply(ψ₊) + 0.7 * white.apply(ψ₊)``. So the rank-N
   composition claim, which was **unchanged and un-weakened** but
   *unstateable*, is stateable again for this pair.

   .. warning::

      **A mix with an un-narrowed leaf no longer announces itself.**
      Do NOT generalize the old raise into a safety property. Measured
      post-B3.4a: ``0.3 * spec + 0.7 * albedo(0.5)`` **runs** and
      returns :math:`|\Gamma_-|`-shaped output that is silently
      wrong — the albedo leaf is a :math:`\Gamma_+ \to \Gamma_+` echo
      and :math:`|\Gamma_+| = |\Gamma_-|` swallows the mismatch. The
      raise was an accident of white's shape check, not a guarantee of
      the algebra; until B3.4b / B3.4c land, a mixed tree containing
      albedo or periodic is unsound and shape-invisible (vv Mode 12).
      ⛔ Both landed (B3.4b and B3.4c on 2026-08-01, periodic's arrow at
      G6.3 step 7 on 2026-08-07), so the last clause is history: every
      leaf the realizer builds is narrowed, and the closure-free albedo
      that made this mix spellable is refused. The lesson above, that
      the raise was never a safety property, stands. What replaces the
      accident is two checks. First, ``RealizedSum`` refuses at
      construction terms that fill different faces, consume different
      outflows or live on different trace spaces (``ValueError``).
      Second, every realized leaf is a :math:`\Gamma_+ \to \Gamma_-`
      image by type: the realizer's leaves are the four realized
      boundary types, each built from the outflow and inflow spaces of
      the trace ladder, and the permutation leaf raises (ERR-045) on an
      inflow ordinate whose partner is not an outflow ordinate.

The specular rows are ``assert_array_equal`` and the choice is
structural, not optimistic: a gather introduces no re-association and
an α-fold is one multiplication of the same two floats, so the
predicted drift is EXACTLY zero — a tolerance there would admit the
failure mode the row exists to catch. (Before the re-anchoring these
rows carried ``nulp ≤ 4``, on a reasoning about the then
``ScaledOperator`` wrapper's extra rounding step. That reasoning is right about the extra
rounding and wrong about its effect: the reference performs the same
multiplication, so both sides round identically.)

The ``.npz`` files live at
``tests/gates/sn/operators/snapshots/bc_equivalence_*.npz`` and are committed to
the repository — the artefacts ARE the verification reference.
``tests/gates/geometry/_generate_bc_equivalence_snapshots.py`` regenerated
them until WF2.4 deleted it (recover it with ``git show
197549f6:tests/gates/geometry/_generate_bc_equivalence_snapshots.py``); no generator is in the
tree now, and the carry fixture's recorded sha256 of each file is the gate on a
regeneration. Since the re-anchoring, regeneration was legitimate only when the
reference EXPRESSION or the probe changed, never to make a red go
away — that would re-anchor the gate on the very code it gates, which
is the shape the inversion removed.

.. note::

   **The** ``albedo_05_lebedev17`` **case was retired at B3.4b.** It
   pinned ``psi_in == 0.5 * psi_out`` on the WHOLE FACE, an artefact of
   a law SN no longer realizes: a bare ``AlbedoBoundary(α)`` was refused
   because its :math:`R = \alpha I` is a :math:`\Gamma_+ \to \Gamma_+`
   endomorphism and its :math:`G` supplies no crossing (the account
   written then is :ref:`bc-method-realizability`). Today the closure-free
   return is unspellable, for a structural reason: the linear part of a
   :class:`~orpheus.transport.boundary_conditions.WallCondition` is a
   ``ResponseNetwork`` whose keys are members of the closed sum
   ``Chain = SpecularChain | DiffuseChain``, and an import-time check
   verifies that each chain composes from outflow to inflow, so no
   chain spells a :math:`\Gamma_+ \to \Gamma_+` return. The
   successor of the case is the second method
   on the ``specular_x_lebedev17`` case, which pinned the same α-fold on
   the same quadrature through the theorem that
   ``AlbedoBoundary(α, SpecularReturn(a))`` realized to :math:`\alpha`
   times the matrix of the mirror ``ReflectiveBoundary(a)`` (it asserted
   ``array_equal`` against :math:`0.5` times the frozen mirror image) —
   and, since the re-anchoring, against the mirror isometry rather
   than against a sibling implementation's recorded output. The live
   harness keeps the α-fold on the ``specular_y_partial_07_LS6`` row
   (``WallCondition.specular(0.7)``) and the mirror-equals-unit-wall
   claim on ``specular_x_lebedev17`` (a ``MirrorDeck`` and
   ``WallCondition.specular(1.0)`` give the same image,
   ``array_equal``).


.. _bc-two-bc-applies-per-matvec:

Two BC apply calls per curvilinear matvec
=========================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

.. important::

   **Historical — this describes the Phase-D-era curvilinear matvec,
   which no longer exists.** Wave O step O.4a.2 deleted the keystone
   re-apply and made the sweep read ``ψ.boundary.inflow`` as *given*, so
   the matvec calls ``bc.apply`` **zero** times; the reflective coupling
   moved to the sibling :math:`-B` (:ref:`bc-extraction`). The Gate-1.5
   test named below was re-posed accordingly — it now pins the
   0-call extraction *and* the emitted WDD outflow, and asserting both
   together is what prevents a silent regression that re-absorbs the BC
   into the matvec. The vacuum realization named in the two call
   descriptions is likewise the pre-B3.2 one; see
   :ref:`bc-vacuum-semantic-correction`. The section is kept because
   the two-call decomposition it analyses is the reason the extraction
   needed a capture-and-compare gate rather than a round-trip check.

Phase C (:ref:`bc-trace-contract-respected-by-matvec`) established
that the SN curvilinear matvec applies the BC trace law **once per
matvec** at the boundary edge, consuming the WDD-propagated outflow
trace :math:`\gamma_+\psi`.  Phase D (Issue #168 Phase D, 2026-05-12,
landed on ``refactor/sn-operator-algebra``) extends the matvec with
a **second** BC apply call, used to build the Carlson coupled-pole
seed's ``bc_outer_value`` (see
:ref:`sn-phase-d-carlson-coupled-pole-sweep` in
:doc:`/theory/methods/sn/curvilinear_numerics` for the math).  The §16A.3 affine trace
law contract is therefore exercised **twice per matvec** in the
post-Phase-D code path:

.. list-table:: BC apply call sequence inside one curvilinear matvec
   :header-rows: 1
   :widths: 12 28 30 30

   * - Order
     - Caller / purpose
     - Input shape & meaning
     - Output use
   * - **#1**
     - Phase D Carlson context build
       (the then-production ``transport_operator_matvec_spherical``
       / ``_cylindrical`` matvec — the whole per-geometry family since
       deleted in the typed-field campaign (#197), successor retired at
       the walk unification (#280 campaigns) — early in the call)
     - ``(N, ng)`` — outer-cell cell-centred :math:`\psi` (NOT the
       face trace; a first-order proxy used only to construct the
       linear-in-:math:`\psi` ``bc_outer_value`` scalar at
       :math:`\mu = -1`)
     - Extract the most-inward-ordinate row; scalar feeds into
       ``CarlsonSweepContext.bc_outer_value`` (this Phase-D context
       object was later retired by Issue #282 route (a) — the
       starting-direction inflow corner is now a typed carrier slot; see
       :ref:`sn-direct-seed-solve`)
   * - **#2**
     - Phase C BC trace law application (at the boundary edge after
       the WDD sweep completes)
     - ``(N, ng)`` — WDD-propagated outflow face values
       :math:`\gamma_+\psi` (the §16A.3 contract input)
     - Fill the inflow rows used as
       :math:`\psi^{\text{face}}_{\text{in}}` for the inward sweep
       phase

The two calls are **structurally distinct**:

* **Call #1** is a Phase D-specific use of the BC operator as a
  *linear-in-ψ* construction of the inward-zero-weight ordinate's
  outer-face flux.  For vacuum BC the then-realized
  ``IncomingOrdinateMaskTensor`` zeroed the inflow ordinate,
  giving ``bc_outer_value = 0``; for reflective BC the
  :class:`PermutationOperator` mirrors outgoing :math:`\leftrightarrow`
  incoming, giving ``bc_outer_value = ψ_cell[N-1]`` (i.e. the
  cell-centred outer-cell value).  Both behaviours preserve operator
  linearity in the input :math:`\psi`.

  The input shape ``(N, ng)`` for cell-centres is a structural
  proxy: the BC trace law expects a trace ``(N, ng)`` shape, and
  feeding it the outer cell-centre array IS the right shape for
  Call #1's linear extraction.  The §16A.3 contract is not
  literally honoured here in the *interpretation* sense (the input
  is not a face trace), but the resulting scalar is linear in the
  matvec's input :math:`\psi` and gives the correct value at
  :math:`\mu = -1` on the only configurations the apply matvec
  cares about (reflective + flat :math:`\psi` :math:`\rightarrow` C,
  vacuum :math:`\rightarrow` 0).  This is the **principled
  shortcut**: a linearly-compatible scalar extraction whose values
  match the canonical inward-zero-weight ordinate's flux on the
  load-bearing test configurations.

* **Call #2** is the canonical Phase C use — the BC trace law
  consumes the **actual** WDD-propagated outflow face trace and
  produces the inflow trace per the §16A.3 affine-bc-form
  contract.  This is the call the Phase C
  :ref:`Gate 1.5 <bc-trace-contract-respected-by-matvec>` test
  pins.

Capture-and-compare Gate 1.5 strengthening
------------------------------------------

The pre-Phase-D Gate 1.5 test was a "round-trip" check: invoke
``bc.realize().apply(...)`` independently and compare against the
matvec's observable output.  Phase D **strengthens** the gate to a
capture-and-compare check that audits the *exact* value the matvec
passes into the BC trace law — necessary because the matvec now
calls ``bc.apply`` twice and the test must locate Call #2 (the
§16A.3 call) unambiguously.

The Phase D test
:func:`tests.gates.sn.sweep.core.test_phase_c_gates.test_bc_trace_contract_capture_and_compare_sphere`
(parametrised over ``vacuum`` and ``reflective``):

#. Monkey-patches ``problem.bc["xmax"].apply`` (the outer radial
   face — a sphere's ``"outer"`` endpoint renders as ``"xmax"``)
   with a recorder wrapper that appends every input array passed to
   it during one matvec call.
#. Independently reconstructs the WDD-propagated outflow trace via
   a reference implementation
   (``_outflow_at_boundary_for_sphere(problem, sig_t, psi_input)``).
#. **Locates Call #2** by matching shape ``(N, ng)`` AND content
   (the captured input that bit-matches the independent reference
   IS the Phase C call).
#. Asserts the located input matches the reference to
   ``rtol=1e-14`` — bit-equal up to FP non-associativity.

The test is **foundation-tagged** (``@pytest.mark.foundation``)
because it pins a software invariant about the matvec's
two-application sequence, not a math claim about the BC operator.
The matching strategy by *both* shape and content protects against a
future regression that adds a third BC apply with the right shape
but wrong content — the test would still locate the canonical Phase
C call provided its content matches the WDD reference.

Both ``vacuum`` and ``reflective`` parametrised cases pass.  The
``vacuum`` case was the load-bearing check because Call #1 produced
non-trivial output under vacuum (the then-realized
``IncomingOrdinateMaskTensor`` zeroed inflow ordinates, so the
extracted ``bc_outer_value`` was zero — but the **input** to Call #1
was the outer cell-centre value, which is **not** zero on a
non-trivial :math:`\psi`).  Locating Call #2 unambiguously required
content matching, not just shape matching.


.. _bc-three-bc-applies-per-sweep-iteration:

BC applies in the SI sweep path
===============================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

Phase D (Issue #168 Phase D, :ref:`bc-two-bc-applies-per-matvec`
above) instituted the *two BC apply calls per curvilinear matvec*
contract on the apply-matvec path (the within-group operator,
:class:`~orpheus.sn.operators.streaming.StreamingCollisionOperator` post-Depth-B; the
matvec then lived at ``_transport_operator_matvec_unified`` — since
deleted at the walk unification (#280 campaigns); the live forward action
is now :meth:`~orpheus.sn.operators.streaming.StreamingOperator.apply`
through the loss-representation walk).  Phase F
(Issue #168 Phase F, 2026-05-12, also landed on
``refactor/sn-operator-algebra``) propagates the same pattern to the
**SI/sweep path** (the then-production ``transport_sweep`` entry — since
retired at the coupled-block campaign step 6 (R-6.1, 2026-07-12) —
dispatching to ``_sweep_1d_spherical`` (the dissolved ``sweep.py``) /
``_sweep_1d_cylindrical``).  See
:ref:`sn-phase-f-carlson-sweep-path-backport` in
:doc:`/theory/methods/sn/curvilinear_numerics` for the math and the
twin-path-fix-incompleteness anti-pattern that motivated the
backport.

The post-Phase-F SI sweep iteration applies the BC operator in
**three** distinct invocations per sweep call — one for the Phase F
Carlson seed, plus :math:`N_{\text{inward ordinates}}` legacy inflow
applications inside the per-ordinate loop (which are not
fundamentally new — they predated Phase D and Phase F).  The
load-bearing addition is the **Phase F seed call**:

.. list-table:: BC apply call sequence inside one SI sphere sweep
   :header-rows: 1
   :widths: 12 28 30 30

   * - Order
     - Caller / purpose
     - Input shape & meaning
     - Output use
   * - **#1**
     - Phase F Carlson seed
       (``_sweep_1d_spherical`` (the dissolved ``sweep.py``) early in
       the call, before the per-ordinate loop)
     - ``(N, ng)`` — persistent outer-face outflow buffer
       ``bc_outer`` carrying the previous outward sweep's
       outgoing flux per ordinate (zero on the first SI
       iteration)
     - Extract the most-inward-ordinate row; scalar feeds into
       :func:`carlson_inward_sweep_from_source` as
       ``bc_outer_value`` (the seed for Hébert (3.434)–(3.435)
       at the outer face)
   * - **#2 … #N₋**
     - Per-inward-ordinate inflow read inside the
       per-ordinate sweep loop
     - ``(N, ng)`` — same ``bc_outer`` buffer (re-read each
       iteration, since intervening outward ordinates may have
       updated it)
     - Read the inflow row ``psi_in_full[n]`` for the current
       inward ordinate :math:`\mu_n < 0` as the spatial
       sweep's incoming flux

**Comparison with the apply-matvec twin** (per
:ref:`bc-two-bc-applies-per-matvec`):

* The apply matvec consolidates its inflow logic into the **single
  Phase C trace law call** at the boundary edge — the BC operator
  is invoked once on the WDD-propagated outflow trace, producing
  the inflow trace per the §16A.3 affine-bc-form contract for ALL
  ordinates simultaneously.
* The SI sweep, by contrast, processes ordinates **sequentially**;
  each inward ordinate independently reads its inflow row from
  the persistent ``bc_outer`` buffer.  The per-ordinate apply
  calls are **not** a §16A.3 trace-law application of the same
  semantic kind — they are *consumer reads* of an already-updated
  buffer.  The Phase F seed call (Call #1) is the only call that
  semantically mirrors Phase D's apply-matvec Call #1 (a
  Phase-specific use of the BC operator as a linear-in-:math:`\psi`
  construction of the inward-zero-weight ordinate's outer-face
  flux).

The Phase F seed call's role is exactly analogous to the
apply-matvec's Phase D Call #1 (per the
:ref:`bc-two-bc-applies-per-matvec` table): a
linear-in-:math:`\psi` extraction of the inward-zero-weight
ordinate's outer-face flux.  Vacuum BC zeros it; reflective BC
yields the most-inward ordinate's mirrored outflow.  In both
contexts the BC operator's **linearity** is the load-bearing
contract — the Phase F helper
:func:`~orpheus.sn.sweep.psi_half_angle_seed.carlson_inward_sweep_from_source`
must remain linear in the input to preserve the SI loop's
fixed-point convergence properties.

The cylindrical path has the analogous structure with a
**per-:math:`\mu`-level** seed call: one BC apply per level,
each invocation extracting the level-specific most-inward
ordinate's row from the same persistent ``bc_outer_cyl``
buffer.  The
``_sweep_1d_cylindrical`` (the dissolved ``sweep.py``) body invokes
the BC operator ``len(quad.level_indices) + N_inward`` times
total per sweep — once per :math:`\mu`-level for the Phase F
seed, plus once per inward ordinate inside each level's
azimuthal loop.

Phase F leaves the §16A.3 BC trace contract semantics
**unchanged**: the SI sweep path does not call the §16A.3 trace
law application that the apply matvec uses at the boundary
edge — the SI sweep updates ``bc_outer`` directly from each
outward ordinate's last visit (line ~593 of
:file:`orpheus/sn/sweep.py`), then the next inward ordinate
reads it via ``apply``.  The semantic contract is the same
(BC operator maps outflow trace :math:`\gamma_+\psi` to inflow
trace :math:`\gamma_-\psi`), but the *invocation pattern* is
per-ordinate sequential rather than once-per-matvec
collective.

No new Gate 1.5 test variant is needed for the Phase F seed
call.  The Phase F bit-identity test module
:mod:`tests.gates.sn.sweep.core.test_sweep_vs_apply_consistency`
(57 foundation tests) pins that the sweep-path's
``bc_outer_value`` extraction matches the apply-path's Phase D
Call #1 result on every test configuration — the structural
invariant that the two paths' Carlson seeds agree on matching
inputs subsumes the BC-apply-input pinning.


.. _bc-curvilinear-realizer-unification:

Curvilinear realizer unification
================================

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

The pre-cleanup architecture carried a **Cartesian / curvilinear
split** at the hub-side resolver then named ``SNProblem._resolve_one``
(retired at #290 P7b — the shipped hook was
``SNProblem.realize_boundary_law`` until wave 2 replaced it by
:meth:`~orpheus.sn.problem.SNProblem.realize_boundary_role`):
the slab and 2-D Cartesian
paths constructed a trace space (then named ``InflowTraceSpace``,
unified into :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
post-#188 and made geometry-blind in C5.3) and routed the BC through
``SNBoundaryRealizer``, while spherical and cylindrical ``Mesh1D``
bypassed the realizer entirely because that factory
(``from_mesh_and_quadrature``, since C5.3
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`)
raised :class:`NotImplementedError` on those coord systems. The bypass
wrapped the bare law instance in
``_BoundBoundaryOperator(law, quadrature=self.quad)``, a
dual-mode shim whose ``apply(psi)`` forwarded
``law.apply(psi, bound_quad)`` to the legacy 2-arg body.

**Why the split existed.** The Wave 2 ``InflowTraceSpace``
factory deferred curvilinear support because no curvilinear-Krylov
consumer at the time needed the per-face mask — the curvilinear
sweep paths computed the inflow / outflow predicate on the fly
inside the inner loop. The deferral was load-bearing for the
shim's dual-mode design: ``quadrature=None`` (Cartesian, wraps a
realized 1-arg op) and ``quadrature=<set>`` (curvilinear, wraps a
legacy 2-arg BC).

**Why the split has been removed.** Issue #188 / C188.1+C188.2
discovered that all three :class:`Mesh1D` coord systems
(``CARTESIAN`` / ``SPHERICAL`` / ``CYLINDRICAL``) share the same
``("left", "right")`` face structure with the radial axis as the
outward normal and the ``GaussLegendre1D``
adapter's :math:`\mu_x` as the direction cosine along that axis —
identically to the slab case. The mask predicate
:math:`\Omega \cdot \hat n < -\epsilon` therefore applies
unchanged. The factory's curvilinear guard was lifted; only 2-D
cylindrical :class:`Mesh2D` (which has no SN sweep in ORPHEUS)
continued to raise :class:`NotImplementedError` (C5.3 later made the
factory geometry-blind, :ref:`sn-c5-geometry-blind-trace`).

Issue #188 / C188.3 then collapsed ``SNProblem._resolve_one``
to a single path: every supported mesh (1-D Cartesian / spherical
/ cylindrical + 2-D Cartesian) built an
``SNMethodSpace.for_face`` and routed through
``SNBoundaryRealizer``. Issue #176 / C176.1 then dropped the
``quadrature=`` kwarg from the shim because no
production-issued shim carried ``_quadrature is not None`` after
C188.3. Issue #176 / C176.3+C176.4 then trimmed the concrete BC
``apply`` signatures to the Option-A interim (keyword-optional
``quadrature=None`` with defensive errors). Issue #186 / B3 + β2
then retired Option A entirely — every concrete BC ``apply``
method was deleted; see
:ref:`bc-trace-law-descriptor-model` for the retrospective.

The architectural sequence is therefore:

* **Issue #188 unblocks Issue #176.** The shim's dual mode existed
  ONLY because curvilinear ``InflowTraceSpace`` could not be
  constructed. Once #188 lifted that, #176's "drop the 2-arg form"
  cleanup became possible without breaking curvilinear sweeps.
* **Issue #176 unblocks Issue #186.** The Option-A interim was
  necessary because dropping ``apply`` outright before the
  curvilinear sweeps consumed realizer output (#188) and the test
  fleet migrated to the realizer-routed contract (#176 / C176.5
  cleanup commits) would have broken curvilinear regression.
  Once those landed, the descriptor cleanup (#186 / B3 + β2)
  became the next step on the architectural ladder.

.. _bc-face-name-carve:

The face-name carve — one crosswalk, one face-keyed BC dict
===========================================================

Wave 8 through Issue #186 settled *how a single boundary law is
realized* (the law / realizer / shim split of
:ref:`bc-overview-three-layers`). What they did **not** settle is
*how the set of resolved laws is keyed and stored on the* ``SNProblem``.
Pre-C4 that storage was a hand-listed per-geometry construction with
named attributes — ``bc_xmin`` / ``bc_xmax`` / ``bc_ymin`` /
``bc_ymax`` (2-D), ``bc_left`` / ``bc_right`` aliases (1-D), plus a
pair of degenerate ``bc_ymin`` / ``bc_ymax`` placeholders on a slab
that no production code ever read. Three separate hand-lists carried
the same ``(axis, endpoint) → "{axis}{min|max}"`` knowledge, and a
fourth hand-list mapped a face name back to a reflection axis. C4
(part of the N-D layout campaign, Issue #220) collapses all four to
**one crosswalk function and one dict-comprehension loop**, keyed by
the same :class:`~orpheus.mesh.axis.FaceLabel` inventory the trace
layout already derives from.

This is the storage-layer counterpart of the realizer unification
in :ref:`bc-curvilinear-realizer-unification`: that section made the
*realization path* uniform across geometries; this one makes the
*storage and keying* uniform across **dimensions**. After C4 a
3-axis mesh yields six face slots and six BC entries with **no edit**
to either producer — the pre-C4 3-branch body would have been
silently wrong the day it was reached. C5 (:ref:`sn-axis-primary-c5`)
admits exactly that 3-axis mesh — not via a hypothetical ``Mesh3D``
dataclass but via the axes tuple directly — and the face-name keying
built here carries it through unchanged.

.. note::

   **What wave 2 (2026-10-24) changed under this carve, and what it kept.**
   The crosswalk, the face inventory as the key set of both the trace layout
   and ``bc``, the structural absence of the pole and the key identity below
   hold as written. What changed is the declaration side and the value side:
   each face's declaration is its role (a deck, or an integer tag bound to a
   :class:`~orpheus.transport.boundary_conditions.WallCondition`), joined by
   :func:`~orpheus.transport.boundary_conditions.bind` into a
   :class:`~orpheus.transport.boundary_conditions.Binding` and realized by
   :func:`~orpheus.transport.method.realize_boundary` through the table
   ``SN_REALIZERS``, and the values of :attr:`SNProblem.bc` are typed
   :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary` operators.
   The passages below that describe the C4-era resolver, its laws and its
   shim are history and say so. The present architecture is
   :ref:`bc-architecture-of-record`; each retired name maps to its successor
   in :ref:`bc-wave2-crosswalk`.

The three string-named layers and the single crosswalk
-------------------------------------------------------

Three SN-side structures key on the same boundary-face string names
``"xmin"`` / ``"xmax"`` / ``"ymin"`` / ``"ymax"`` (and, in C5,
``"zmin"`` / ``"zmax"``):

.. list-table:: The three face-name-keyed layers and their shared crosswalk
   :header-rows: 1
   :widths: 26 30 44

   * - Layer
     - Structure
     - Role
   * - **Face layout**
     - :attr:`SNProblem.boundary_face_layout`
       (:class:`~orpheus.numerics.face_layout.FaceLayout`)
     - The flat-buffer descriptor: which faces exist, each face's
       per-face shape ``(N, ng, *codim-1 cells)``, and the offsets
       that pack them into one backing array.
   * - **Trace space**
     - :attr:`SNProblem._trace`
       (:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`)
     - The inner-product geometry on the boundary: per-face
       inflow / outflow ordinate masks over the signed
       :math:`\Omega\cdot\hat n` it carries (``trace.layout.faces``
       reproduces the same names).
   * - **BC dict**
     - :attr:`SNProblem.bc`
       (``dict[str, RealizedBoundary]``,
       :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`)
     - The realized boundary operator per face with a trace — the typed
       :math:`\Gamma_+ \to \Gamma_-` map of the face's declaration (a deck,
       or the condition its tag is bound to). At C4 the value was the
       ``_BoundBoundaryOperator`` shim around the 1-arg realized law.

Pre-C4 each of these grew its key set from its own per-geometry
hand-list. The crosswalk knowledge — "axis 0 is ``x``; the
``min`` / ``max`` endpoints suffix the axis name; a solid radial
axis's single ``outer`` endpoint renders as the ``max`` face"
— was duplicated at the layout builder, the BC resolver, and the
reflective-axis dispatch. C4 lifts that knowledge into **one**
single-sourced rendering on the structural face key:

.. code-block:: python

   # orpheus/mesh/axis.py
   from orpheus.geometry.faces import AXIS_NAMES        # ("x", "y", "z")
   _ENDPOINT_SUFFIX = {"min": "min", "max": "max", "outer": "max"}

   @dataclass(frozen=True, slots=True)
   class FaceLabel:
       axis_index: int
       endpoint: str   # "min" / "max" / "outer"

       @property
       def face_name(self) -> str:
           suffix = _ENDPOINT_SUFFIX.get(self.endpoint)
           if suffix is None:                       # fail loud
               raise ValueError(...)
           return f"{AXIS_NAMES[self.axis_index]}{suffix}"

:attr:`FaceLabel.face_name <orpheus.mesh.axis.FaceLabel.face_name>`
is THE rendering of the structural identity ``(axis_index,
endpoint)`` into the ``"{axis}{min|max}"`` string world. Both
producers — :attr:`SNProblem.boundary_face_layout` and the shared
:func:`~orpheus.transport.method.realize_boundary` body (its ``traced`` set
is ``{label.face_name for label in face_labels(method.axes)}``) — read it, so
a key drift between the face layout and the BC dict is **unrepresentable by
construction**: they cannot disagree because they read the same function over
the same :func:`~orpheus.mesh.axis.face_labels` inventory. (Until wave 2 the
second producer was the C4-era ``SNProblem.realize_boundary_law`` hook, called
per face by ``resolve_boundary_conditions``.)

.. note::

   ``AXIS_NAMES`` moved **down** from
   :mod:`orpheus.sn.loss_representation.sweep_graph` to the axis module
   (then ``orpheus.transport.mesh.axis``, now :mod:`orpheus.mesh.axis`) in
   C4 — to the bottom of the SN dependency graph, next to the axis
   primitives it names; a later step moved it down again, to
   :mod:`orpheus.numerics.face_layout` (below), and since layer DAG WP-B
   (2026-10-22) it lives in :mod:`orpheus.geometry.faces`, from which the
   axis module re-exports it. ``sweep_graph`` re-exported it only outward;
   ``sweep_schedule`` and ``loss_representation`` now import it
   downward. At C4 this put the single source of the axis↔name crosswalk
   in the same module as :class:`~orpheus.mesh.axis.FaceLabel`, the
   walk's in/outflow-face derivation, and the schedule's
   outgoing-face derivation. Today the crosswalk is single-sourced in
   :mod:`orpheus.geometry.faces`; ``FaceLabel.face_name`` and the schedule
   read it from there, and no consumer hand-lists
   ``("x", ...), ("y", ...)`` pairs.

The ``"outer" → "max"`` convention and fail-loud
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A solid sphere or cylinder is a :class:`~orpheus.mesh.axis.RadialAxisMesh`
with a single ``"outer"`` endpoint (the pole at :math:`r=0` is
**not** an endpoint — see :ref:`bc-pole-structural-absence`). The
crosswalk renders ``"outer"`` as the ``max`` face of its axis, so a
sphere's outer radius is keyed ``"xmax"``. This is the **historical
curvilinear convention** — every curvilinear boundary operator,
trace face name, and sweep schedule already keys the outer radius on
``"xmax"`` — preserved verbatim by ``_ENDPOINT_SUFFIX["outer"] =
"max"`` rather than re-derived.

Any endpoint label that is **not** one of the three canonical
strings (``"min"`` / ``"max"`` / ``"outer"``) raises
:class:`ValueError`. An :class:`~orpheus.mesh.axis.AxisMesh` exposes
user-overridable ``label_low`` / ``label_high`` fields (a slab user
may rename them ``"left"`` / ``"right"`` for convention); such a
renamed endpoint has **no face name** and must fail loud rather than
silently desynchronize from the ``"{axis}{min|max}"`` world that the
three layers key on. The failure surfaces at L0 — at the crosswalk
itself — not three layers up as a mis-keyed boundary operator or a
``KeyError`` deep inside a sweep.

The derivation chain — face_labels → {layout, bc}
-------------------------------------------------

The SN phase space factors as a tensor product of per-axis 1-D
meshes (grand report §15.1). The axis tuple
:attr:`SNProblem.axes <orpheus.sn.problem.SNProblem.axes>` is therefore the root of every boundary-keyed
structure:

.. code-block:: text

   SNProblem.axes  (tuple[Axis1D, ...])
        │
        │  face_labels(axes) — one FaceLabel per (axis, endpoint),
        │  iterated axis-ascending then endpoint-in-axis-order
        ▼
   SNProblem.face_labels  (tuple[FaceLabel, ...])
        │
        ├─── boundary_face_layout : one slot per label,
        │       named  label.face_name,  shaped (N, ng, *face_shape(label))
        │
        └─── bc : one entry per label,
                keyed  label.face_name,  realized from that face's role
                axes[label.axis_index].roles[label.endpoint]
                (a deck, or a tag bound to a condition)

The **declaration** source for each face is its role: each axis carries one
role per endpoint (``AxisMesh.roles``, a deck or an ``int`` tag, on
:class:`~orpheus.mesh.axis.AxisMesh`), keyed by the *same* axes that
``face_labels`` derives the labels from, and the material mesh's ``face_roles``
joins them under the face names. The face inventory **is** the role inventory:
a face that exists has exactly one role; a face that does not (the pole) has no
label and no role. A tag is bound to its condition by the problem's declared
conditions: :func:`~orpheus.transport.boundary_conditions.bind` joins the two
into a :class:`~orpheus.transport.boundary_conditions.Binding`, which is an
exact cover of the tags (a tag no condition names, or a condition no face's tag
names, is refused) and, where the colouring group is computed, invariant under
it: a condition that is not face-local on a carrier whose group is not computed
(``UncomputedSymmetries``, every ``Mesh2D``) is refused rather than checked
(SCOPE-BOUNDARY, #658; ``binding.py``). The
resolution is one loop over the binding's faces, the body of
:func:`~orpheus.transport.method.realize_boundary` in
``orpheus/transport/method.py`` (quoted with its annotations and docstring
elided):

.. code-block:: python

   def realize_boundary(method):
       table, binding = method.boundary_realizers, method.binding
       traced = {label.face_name for label in face_labels(method.axes)}
       realized = {}
       for face in binding.faces:
           declared = binding.decks[face] if face in binding.decks else binding.condition_of(face)
           if face in traced:
               realized[face] = method.realize_boundary_role(declared, face)
           else:
               table.require(declared, method.face_site(face))
       return realized

A glued face is realized from its deck, a tagged face from the condition its tag
is bound to. A face with a trace slot goes to the method's hook,
:meth:`~orpheus.sn.problem.SNProblem.realize_boundary_role`; a face without one
(the inner face of a hollow radial body: a tag with no inner-surface trace) only
meets the table through ``require``, which raises the same refusal. Every axis
carries a role on each endpoint (``None`` is not a role, and the axis
constructors refuse it), so the declaration is read verbatim and there is no
default to fall back on (:ref:`structured-geometry-no-default-law`).

The C4-era body of this loop was ``resolve_boundary_conditions``: for each label
it read ``axes[label.axis_index].bc[label.endpoint]``, turned the tag into a law
with ``_law_from_tag`` and called ``method.realize_boundary_law(law,
label.face_name)``. Wave 2 retired it with the string tables it read (see
:ref:`bc-wave2-crosswalk`). Its first form was the mesh-private
``SNProblem._resolve_one``, which owned the per-face loop itself (#290 P7b moved
the loop up and renamed the hook), and it read ``… or default`` with ``default =
BC("reflective")``, filling an undeclared axis law as reflective, until P1 step
3c of #405 (2026-09-30) retired that default together with ``None`` as a
declaration. The realizer plumbing it drove (registry →
``SNMethodSpace.for_face`` → ``SNBoundaryRealizer`` → ``_BoundBoundaryOperator``)
was **unchanged** from the pre-C4 path: C4 changed only the *keying and
storage*, not the *realization*, so the resolved operators were bit-identical
objects to the pre-C4 ones (see :ref:`bc-face-name-carve-verification`). Today
the hook builds ``SNMethodSpace(quadrature, trace)`` and calls
``SN_REALIZERS.realize(role, face_site(face), method_space)``, which admits the
declaration on the face first and then returns the face's
:class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`.

The :attr:`boundary_face_layout` producer is the dual comprehension:

.. code-block:: python

   # orpheus/sn/problem.py — boundary_face_layout
   return FaceLayout.from_named_shapes([
       (label.face_name, (N, self.ng, *self.face_shape(label), *moment_tail))
       for label in self.face_labels
   ])

(``moment_tail`` is the trailing transverse-moment axis of a multi-moment
scheme, empty for diamond difference and step; it did not exist at C4.)

The pre-C4 body of each producer was a 3-branch ``isinstance`` /
coord-system split (1-D slab / 1-D curvilinear / 2-D Cartesian) that
hand-listed the face names. The dict-comprehension slot order
reproduces the historical hand-listed order **byte-for-byte** —
``face_labels`` iterates axis-ascending then endpoint-in-axis-order,
which is exactly the order the hand-lists used. The affine
``sha256`` goldens stayed byte-identical across the carve.

.. _bc-pole-structural-absence:

The pole is structurally absent, not null (Pattern 4 sharpened)
---------------------------------------------------------------

A :class:`~orpheus.mesh.axis.RadialAxisMesh` has
``endpoints = ("outer",)`` — exactly one BC-bearing endpoint. A
solid sphere or cylinder therefore has a ``bc`` dict with **exactly
one entry** (``"xmax"``) and **no pole entry**. The geometric pole
at :math:`r=0` is the angular closure's regularity condition (the
:math:`1-\mu^2` redistribution coefficient vanishes at
:math:`\mu=\pm 1`; the inward sweep seeds from a moment-folded
source at :math:`\mu=-1` — see
:ref:`sn-pole-angular-closure-protocol` in
:doc:`/theory/methods/sn/curvilinear_one_group`), not a boundary condition.

C4 **sharpens** the pre-existing Pattern-4 treatment of the pole.
Pre-C4 the pole-as-non-BC was spelled by an explicit null:
``bc_left = bc_xmin = None`` — a named attribute that *exists* and
*holds* ``None``. Post-C4 the pole-BC is **structurally absent**: a
dict key that does not exist. Asking for it is a :class:`KeyError`,
not a ``None`` that a consumer might forget to guard:

.. list-table:: Pole-BC representation before and after C4
   :header-rows: 1
   :widths: 30 35 35

   * - Aspect
     - Pre-C4 (``None`` placeholder)
     - Post-C4 (structural absence)
   * - Sphere ``bc`` surface
     - ``bc_xmin = None``, ``bc_xmax = <op>``
     - ``bc = {"xmax": <op>}``
   * - Asking for the pole
     - ``problem.bc_xmin`` → ``None``
     - ``problem.bc["xmin"]`` → :class:`KeyError`
   * - Failure mode of a buggy consumer
     - silent ``None.apply`` → ``AttributeError`` deep
       in a sweep, or a guard that *should* exist but
       doesn't
     - immediate ``KeyError`` at the access site — the
       illegal access is unrepresentable rather than null

This is the illegal-states-unrepresentable principle (Pattern 4)
applied to the dict: a pole-BC is not "a BC that is ``None``", it is
"not a face at all". The :func:`~orpheus.mesh.axis.face_labels`
inventory simply does not emit a label for the pole, so neither
producer writes a slot or entry for it.

.. _bc-face-name-latent-d3-bug:

The latent d=3 axis-dispatch bug, closed by construction
--------------------------------------------------------

Before C4, ``SNProblem._resolve_one`` derived a reflective law's
reflection axis from a hand-listed membership test::

    axis = "y" if face in ("ymin", "ymax") else "x"

This is correct at :math:`d \le 2` by string coincidence — every
non-``y`` face is on the ``x`` axis. But a ``"zmin"`` / ``"zmax"``
face (a 3-axis mesh) would have fallen into the ``else`` branch and
silently built the **X-axis** reflection permutation — the *wrong
reflection partner*. A reflective law that reflects across the wrong
axis is a ``vv-principles`` **Mode-9 class** error (a law that is
wrong on a configuration the degenerate lower-dimensional test never
reaches): it would produce a plausible-but-wrong converged flux on a
3-D reflective problem, invisible to any :math:`d \le 2` test.

C4 derived the axis from the label's own
``AXIS_NAMES[label.axis_index]``, so the reflection partner was
correct at **any** dimension by construction::

    law = ReflectiveBoundary(axis=AXIS_NAMES[label.axis_index])

This is the boundary-resolution sibling of the C3.6 finding that a
z-face never sheds in the in-plane projection — both are latent
3-D correctness traps closed *before* ``Mesh3D`` exists, so that the
N-D layout campaign reaches C5 with the boundary keying already
3-D-correct.

Wave 2 made the class unspellable instead of correct by derivation: no
declaration carries an axis name. A reflecting Cartesian face is a
``MirrorDeck()`` and a reflecting curved wall a ``WallCondition.specular(...)``;
:func:`~orpheus.geometry.deck.pairing_of` reads the mirror's plane from the
face's own name (``face_normal(face)``, then ``mirror_motion(AXIS_NAMES[axis],
dimension)``), and the ordinate permutation is ``ordinate_pairing(quadrature,
motion)``. The d=2 pin ``test_2d_reflective_y_face_builds_y_axis_permutation``
(below) still asserts that the ``ymin`` / ``ymax`` faces get the y mirror.

.. note::

   **No ERR entry was filed for the latent d=3 axis bug.** When C4
   landed, no 3-axis mesh was constructible (the axis-native d=3
   admission did not arrive until C5 — :ref:`sn-axis-primary-c5`), so
   no 3-axis mesh had ever reached ``_resolve_one`` and **no production
   bug ever shipped**: the hand-listed dispatch was retired *before*
   any caller could exercise its wrong branch. The d=2 observable proxy
   (``test_2d_reflective_y_face_builds_y_axis_permutation``, with a
   non-vacuity guard that the x- and y-reflection maps differ under
   Lebedev) pins the *structural* correctness of the per-label axis
   derivation, which is what makes the d=3 extension correct by
   construction. The error catalog records *shipped* L0-caught bugs;
   a defect closed by construction before its triggering type exists
   is documented here, not in ``error_catalog.rst``.

Why string-keyed, not FaceLabel-keyed
-------------------------------------

Issue #220 allowed either a ``dict[str, ...]`` keyed by face name or
a ``dict[FaceLabel, ...]`` keyed by the structural label directly.
C4 chose **string-keyed**, isomorphic to
:attr:`FaceLayout.faces <orpheus.numerics.face_layout.FaceLayout>`,
for one reason: **every consumer iterates ``trace.layout.faces``
(strings)**. The within-group operator
:class:`~orpheus.sn.operators.boundary.SNBoundaryOperator` (it reads
``problem.bc[face]``; its C4-era helper was ``_face_laws``) and the schedule's
``reflective_faces`` both walk the trace layout's face-name strings
and index the BC by that string. A ``FaceLabel``-keyed dict would
force a reverse ``name → label`` lookup at *every* consumer, re-deriving
the very crosswalk C4 single-sources.

:class:`~orpheus.mesh.axis.FaceLabel` remains the **structural source**
— it is the load-bearing key for the dim-agnostic face inventory,
the outflow-ordinate mask cache, and the sweep DAG's face-trace
state. ``face_name`` is its *single rendering* into the string world
the consumers already speak. One crosswalk function, called by both
producers, makes the keys of the layout and the BC dict identical by
construction:

.. math::
   :label: bc-face-name-key-identity

   \operatorname{set}(\texttt{sn\_mesh.bc})
   = \operatorname{set}(\texttt{boundary\_face\_layout.faces})
   = \{\, \ell.\texttt{face\_name} : \ell \in \texttt{face\_labels} \,\}

.. (vv-status rationale) Structural by-construction identity: the BC-dict keys,
   the boundary-face-layout faces, and the FaceLabel.face_name renderings are
   the same set because both producers call one crosswalk. The crosswalk is
   pinned by the foundation gate ``tests/gates/sn/primitives/test_face_name_crosswalk.py``
   (the exhaustive (axis,endpoint)→face_name table + fail-loud negatives). A
   single-source-of-truth structural identity, not a solver claim.
.. vv-status: bc-face-name-key-identity documented

— a set equality that *cannot drift* because both sides are the same
comprehension over the same inventory.

.. _bc-face-name-carve-what-retired:

What C4 retired
---------------

C4 retires every pre-C4 named-attribute spelling of the resolved BC
surface, with no deprecation shim (aggressive retirement — a
read-through ``@property`` outliving its merge cycle would be the
very desync the carve removes):

* **The named instance attributes** ``bc_xmin`` / ``bc_xmax`` /
  ``bc_ymin`` / ``bc_ymax`` (2-D) and the ``bc_left`` / ``bc_right``
  aliases (1-D). Consumers now key into :attr:`SNProblem.bc` by face
  name. Accessing a retired attribute is an :class:`AttributeError`.
* **The degenerate 1-D y-face placeholders.** Pre-C4, a slab
  :class:`SNProblem` carried a pair of realized no-op
  ``ReflectiveBoundary(axis="y")`` operators at ``bc_ymin`` /
  ``bc_ymax``, routed through ``SNMethodSpace.minimal`` so
  cross-dimensional code could read them without coord-system
  gating. **No production code ever read them**: a 1-D Problem's
  ``trace.layout.faces`` is ``("xmin", "xmax")``, so the generic
  consumers (which iterate the trace layout) never asked for a
  y-face. The placeholders were a uniformity affordance with no
  consumer — exactly the kind of dead realized state the
  face-labels-derived dict makes unrepresentable: a slab has no
  y-axis in its :attr:`~orpheus.sn.problem.SNProblem.axes` tuple, so
  :func:`~orpheus.mesh.axis.face_labels` emits no y-label, so
  :attr:`bc` has no y-entry, so ``slab.bc["ymin"]`` is a
  :class:`KeyError`. (Pre-C4 design rationale for *why the
  placeholders were once safe* is preserved in the
  :ref:`bc-curvilinear-realizer-unification` history above; C4
  removes the need for the rationale by removing the placeholders.)
* **The hand-listed reflective-axis dispatch**
  (``"y" if face in (...) else "x"``) — see
  :ref:`bc-face-name-latent-d3-bug`.

Consumers migrated in C4: ``SNBoundaryOperator._face_laws``
(within-group operator) and ``sweep_schedule.reflective_faces``
(the schedule) both changed from ``getattr(mesh, f"bc_{face}")`` to
``mesh.bc[face]``, iterating over ``trace.layout.faces`` exactly as
before.

.. _bc-face-name-carve-verification:

Verification — bit-identity by inheritance + new L0 pins
--------------------------------------------------------

C4 is a **structural** carve: it changes keying and storage, not
the realized operators or any numerical value. The verification
strategy reflects that:

* **Bit-identity by inheritance.** A BC realization is object
  construction; the realizer plumbing (registry →
  ``SNMethodSpace.for_face`` → ``SNBoundaryRealizer`` →
  ``_BoundBoundaryOperator``) was unchanged — it was reached through
  the C4-era ``_resolve_one`` (today through
  :meth:`~orpheus.sn.problem.SNProblem.realize_boundary_role` and the
  table ``SN_REALIZERS``). The resolved
  operators were the same objects as before, so every solver test
  that exercised them inherited its prior verification. The affine
  ``sha256`` goldens stayed byte-identical; the broad sweep /
  operators / primitives / solve suite was green.
* **L0 crosswalk pins**
  (:mod:`tests.gates.sn.primitives.test_face_name_crosswalk`,
  foundation-tagged). An exhaustive **hand-transcribed**
  ``(axis, endpoint) → face-name`` table for :math:`d \in \{1,2,3\}`
  (mirror-not-import, so the test is not tautological against the
  production derivation it verifies), the d=3 ``z``-face admission
  (the crosswalk is a pure function, so the 3-axis rendering is
  verifiable **now** with no ``Mesh3D``), and both fail-loud
  negatives (a non-canonical endpoint → :class:`ValueError`; an
  axis beyond the named inventory → :class:`IndexError`).
* **L0 bc-dict / face-layout inventory pins**
  (:mod:`tests.gates.sn.operators.test_snmesh_realizer_wiring`,
  foundation-tagged):

  - ``test_bc_inventory_equals_face_layout_across_geometries`` —
    ``set(sn.bc) == set(boundary_face_layout.faces)`` across slab
    (2 faces), 2-D Cartesian (4), sphere (1), cylinder (1), the
    Issue #220 acceptance set.
  - ``test_2d_reflective_y_face_builds_y_axis_permutation`` — the
    d=2 observable proxy for the latent d=3 axis-dispatch bug, with
    a non-vacuity guard asserting the x- and y-reflection maps
    differ under Lebedev (else the pin would be vacuous).
  - ``test_bc_dict_misses_and_retired_attributes_fail_loud`` — a
    face that does not exist is a :class:`KeyError` (plain dict, no
    masking default); every retired named attribute is an
    :class:`AttributeError` (no silent ``None``-shim survives).

The test-architect verification design memo is
``.claude/agent-memory/test-architect/c4_snmesh_bc_dict_verification.md``.

C4 closure
----------

The face-name carve lands in two parts under Issue #220 (the SN N-D
layout campaign): C4.1 (the :attr:`FaceLabel.face_name` crosswalk +
the :attr:`boundary_face_layout` loop collapse) and the bc-dict
resolution loop + named-attribute retirement. The carve is byte-identical
on every numerical output (affine ``sha256`` goldens unchanged) and
leaves the SN boundary keying **3-D-correct by construction** ahead of
``Mesh3D`` (C5). The realizer plumbing it sits on top of was unified by
the predecessor close-out below.

Predecessor closure — curvilinear realizer unification (Issue #188 + #176 + #186)
---------------------------------------------------------------------------------

This section records the trace-law architecture that wave 2 retired (2026-10-24); the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its successor in :ref:`bc-wave2-crosswalk`.

The realizer path that C4 keys into was made uniform across geometries
by the curvilinear-realizer-unification arc
(:ref:`bc-curvilinear-realizer-unification`), closed by branch
``feature/bc-curvilinear-realizer-cleanup``
(2026-05-11). Three GitHub issues converged on that branch:

* **Issue #188** — curvilinear ``InflowTraceSpace`` support
  (commits ``9cf2b0a`` + ``17067d5``). Lifted the
  :class:`NotImplementedError` guard on spherical / cylindrical
  Mesh1D so every supported mesh can build a per-face inflow mask.
* **Issue #176** — drop 2-arg ``apply`` + simplify shim (commits
  ``cf29ce4`` + ``a4a43c2`` + ``913e501`` + ``188bf9a``). Collapsed
  the dual-mode shim into a strict 1-arg passthrough; landed the
  Option-A interim with keyword-optional ``quadrature=None`` on
  the concrete laws.
* **Issue #186 (B3 + β2)** — pure-descriptor cleanup (commits
  ``f71a32c`` + ``da414eb`` + ``89d09a4`` + ``633cc69`` +
  ``bb674da`` + the test-migration trail). Retired the Option-A
  ``apply`` methods, dropped :class:`LinearOperator`
  inheritance from ``BoundaryTraceLaw``, and formalised the
  descriptor-tree composition algebra via the new
  ``LawSum`` / ``LawScaled`` types. The architectural
  sequence is therefore Issue #188 → #176 → #186: each step
  unblocked the next.

The ``_BoundBoundaryOperator`` shim survived that closure because a resolved BC
had to carry **both faces of its realization** — what it does
(``inner``) and what it means (``law``, since phase B2.0). Its
``kind``-string tag was the older, weaker reason: load-bearing for the
BC-resolution diagnostic and several ``problem.bc["xmin"] ==
"vacuum"``-style test sites (the dict-keyed spelling since C4 / #220;
this was ``problem.bc_left == "vacuum"`` pre-C4), and then a read-through
of ``law`` rather than a stored copy — phase B2.2 retired it, and the
shim survived that retirement on the descriptor alone. The dual-mode
bound-quadrature
backing was gone (#176), and the ``*_extra, **_kw`` swallow was
gone (#186 / C-B3.4). Every supported mesh consumed a strict
1-arg :class:`LinearOperator` produced by
``SNBoundaryRealizer`` for single BCs, or by
``realize_recursively`` for rank-N descriptor trees. Wave 2 retired the shim (see
:ref:`bc-wave2-crosswalk`): its successor
:class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary` is itself the
operator, carries ``inner``, and answers the face's facts
(``return_fraction``, ``permutes_ordinates``, ``domain_face``, ``is_vacuum``) as
properties of its type, so no consumer asks for a law.

Plan documents:

* ``.claude/plans/transient-giggling-cake.md`` — the foundational
  12-wave BC trace-law refactor plan (Waves 0–12 close-out
  documented at :ref:`theory-boundary-conditions`).
* ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` —
  the #188 + #176 cleanup plan (Option-A landing).
* ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` — the
  Issue #186 B3 + β2 cleanup plan (descriptor-model landing).

At this closure Grand Report v3 §16A.3 (the three-layer architecture) was
**enforced by the type system** — descriptors had no ``apply``, operators did —
and §16A.5 (the trace-correct vacuum representation) was uniform across coord
systems, the legacy zeros-all path gone. Both hold under the wave-2 types: a
:class:`~orpheus.transport.boundary_conditions.WallCondition` is data and has no
``apply`` (the realized operator has), and a vacuum wall realizes to the zero
map :math:`\Gamma_+ \to \Gamma_-` of the face it is installed on
(``ZeroMap``, the image of the empty response network).


.. _sn-axis-primary-c5:

The axis-primary inversion and 3-D admission
============================================

C4 (:ref:`bc-face-name-carve`) made the *boundary keying*
dimension-agnostic. C5 makes the **whole mesh** dimension-agnostic and
then admits the first 3-axis Cartesian :class:`SNProblem` — *without* a
``Mesh3D`` dataclass. The design fork (resolved by the user,
2026-06-11) is **axis-native**: a 3-D problem enters ORPHEUS only
through :meth:`SNProblem.from_axes` with a 3-tuple of
:class:`~orpheus.mesh.axis.AxisMesh`. :class:`~orpheus.mesh.structured.Mesh1D`
and :class:`~orpheus.mesh.structured.Mesh2D` stay the :math:`d \le 2`
user-facing surface, bit-identical to before
(``sha256`` affine goldens unchanged, no regeneration). A ``Mesh3D``
would have had **exactly one consumer** (SN — ``cp`` / ``mc`` / ``moc`` /
``diffusion`` consume zero ``Mesh2D``); the "Unify after two instances"
discipline forbids minting a base type for a single consumer.

The campaign's keystone insight, surfaced by the C5 elegance audit, was
that the d=3 admission could not be a clean *extension* until a
pre-existing **data-flow inversion** in the constructor was repaired. C5
is therefore sequenced *clean before extend*: C5.1–C5.4 invert and
de-phantom the hub, and only then C5.5 admits d=3 as a
one-line gate removal.

.. _sn-c5-lossy-roundtrip:

Pre-C5: the lossy axes → mesh → axes round-trip
-----------------------------------------------

The SN phase space factors as a tensor product of per-axis 1-D meshes
(grand report §15.1); the natural primary representation of an
:class:`SNProblem` is therefore its **axes tuple**
:attr:`SNProblem.axes <orpheus.sn.problem.SNProblem.axes>`. Pre-C5.1 the constructor did not treat it that
way. :meth:`SNProblem.from_axes` *synthesized a legacy*
:class:`~orpheus.mesh.structured.Mesh1D` / :class:`~orpheus.mesh.structured.Mesh2D`
from the caller's axes (via ``legacy_mesh_from_axes``), handed that
mesh to ``__init__``, and ``__init__`` then **discarded the caller's
tuple and re-derived the axes from the synthesized mesh**:

.. code-block:: text

   from_axes(axes)                     __init__(mesh, ...)
        │                                    │
        │  legacy_mesh_from_axes(axes)       │  axes = axes_from(mesh)   ← re-derived
        ▼                                    ▼
   Mesh1D / Mesh2D  ───────────────────►  self.axes   (NOT the caller's tuple)

This ``axes → mesh → axes`` round-trip is **lossy** in two ways, and
its existence was the structural reason d=3 appeared to need a "third
construction arm":

1. **Custom endpoint labels were silently reset.** An
   :class:`~orpheus.mesh.axis.AxisMesh` carries user-overridable
   ``label_low`` / ``label_high`` fields (a slab user may name them
   ``"left"`` / ``"right"``). The legacy mesh has no slot for those
   labels, so the round-trip dropped them and the re-derived axes came
   back with default labels — a silent desync of exactly the kind C4's
   :attr:`FaceLabel.face_name <orpheus.mesh.axis.FaceLabel.face_name>`
   crosswalk relies on never happening.
2. **d=3 had nowhere to round-trip *through*.** A 3-axis tuple cannot
   synthesize a ``Mesh1D`` or ``Mesh2D``, so the inverted flow
   *mandated* a legacy mesh at every dimension — which is exactly the
   ``d \ge 3`` blocker the user directive named ("clean before
   extending").

C5.1 inverts the flow: the axes tuple becomes primary, stored verbatim.

.. _sn-c5-axis-primary-construction:

The axis-primary construction — one body, verbatim axes
-------------------------------------------------------

Post-C5.1, **both** entry surfaces funnel into one private body,
``_init_core``, which stores the axes tuple as-is:

.. code-block:: text

   SNProblem(mesh, ...)   ──►  axes = axes_from_legacy_mesh(mesh)   (convert ONCE,
                                                                  at the inbound
                                                                  boundary)
   from_axes(axes, ...) ─►  axes  (stored verbatim — the caller's tuple)
                                            │
                                            ▼
                                       _init_core(axes, ...)

The legacy ``SNProblem(mesh, ...)`` surface converts via
``axes_from_legacy_mesh`` **once**, at the inbound boundary
(parse-don't-validate); :meth:`SNProblem.from_axes` stores the caller's
tuple directly. There is no longer an ``axes → mesh → axes``
round-trip — the conversion is one-directional, ``mesh → axes``, and
only on the legacy surface.

The pre-C5.1 constructor branched on ``isinstance(mesh, Mesh1D)`` vs
``isinstance(mesh, Mesh2D)`` to compute per-dimension metadata (cell
widths, spatial shape). That **isinstance metadata branch dissolves**
into axis-derived properties. The single load-bearing identity is that
per-axis cell widths come from the axis edges:

.. math::
   :label: sn-axis-widths

   \texttt{axis\_widths}[a] = \operatorname{np.diff}(\texttt{axes}[a].\texttt{edges})

.. (vv-status rationale) Representational bit-identity: per-axis cell widths are
   np.diff of the axis edges, byte-identical to the retired Mesh1D.widths /
   Mesh2D.dx / Mesh2D.dy spellings. Pinned by the bit-identity gates
   ``tests/gates/sn/primitives/test_axis_native_construction.py``
   (``test_d2_metadata_byte_identical_axis_vs_legacy`` /
   ``test_1d_slab_metadata_byte_identical_axis_vs_legacy``). A carve
   representational identity, not a solver claim.
.. vv-status: sn-axis-widths documented

This is **bitwise identical** to the legacy per-dataclass spellings it
replaces — :attr:`Mesh1D.widths <orpheus.mesh.structured.Mesh1D>`,
:attr:`Mesh2D.dx <orpheus.mesh.structured.Mesh2D>`, and
:attr:`Mesh2D.dy <orpheus.mesh.structured.Mesh2D>` are each
``np.diff(edges)`` over the same edge arrays (``mesh.py:287`` /
``:567`` / ``:572``), so the carve produces the same floating-point
bytes. The whole-mesh coordinate system is likewise derived from the
per-axis coordinates by a new pure primitive
:func:`~orpheus.mesh.axis.coord_system` (a multi-axis mesh must be
all-Cartesian); the constructor's reduced-operator dispatch and the
angular-closure default now read the **axis-derived** :attr:`SNProblem.coord <orpheus.sn.problem.SNProblem.coord>`,
not ``mesh.coord``.

After C5.1, the ``mesh`` attribute of
:class:`~orpheus.sn.problem.SNProblem` is **inbound provenance
only** — it records *which legacy mesh the caller passed, if any*. It is
``None`` when the Problem was built from axes at :math:`d \ge 3` (no legacy
mesh exists to record; ``problem.py`` spells the branch
``legacy_mesh_from_axes(axes, mat_map=mat_map) if len(axes) <= 2 else
None``). It is written here as a literal rather than an ``:attr:`` role
because the base ``MaterialMesh`` sets it on the *instance* with no
class-level annotation, so autodoc mints no target for it and a role
would render as unlinked text. A handful of :math:`d \le 2` consumers (the 1-D
reduced streaming constructors, the trace build, realizer metadata)
still read ``self.mesh`` at C5.1; each dissolves across C5.2–C5.5 as
its datum is repointed to an axis-native source. ``legacy_mesh_from_axes``
narrows from a round-trip *source* to a :math:`d \le 2` **adapter**
synthesis for those remaining consumers.

Custom endpoint labels now fail loud (C4 doctrine)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

With the axes stored verbatim, a custom endpoint label survives
construction — and therefore reaches the
:attr:`FaceLabel.face_name <orpheus.mesh.axis.FaceLabel.face_name>`
crosswalk. A label that is **not** one of the canonical strings
(``"min"`` / ``"max"`` / ``"outer"``) now raises :class:`ValueError`
**at construction** (the crosswalk's fail-loud — see
:ref:`bc-face-name-carve`), rather than being silently normalized away
by the round-trip. This is the C4 doctrine made operative by C5.1:
**overridable labels cannot silently desync the crosswalk**, so a label
the three face-keyed layers cannot key on must fail at L0, at the
construction site, not three layers up as a mis-keyed boundary
operator or a ``KeyError`` deep inside a sweep. (Pre-C5.1 the
round-trip *masked* this class of user error by discarding the label;
that masking is gone.)

.. _sn-c5-phantom-retirement:

The phantom shims retire (ny, dy, dx)
-------------------------------------

With per-axis widths and the rank-generic :attr:`SNProblem.spatial_shape <orpheus.sn.problem.SNProblem.spatial_shape>`
now native, the legacy phantom-bearing metadata retires. Every spatial
read becomes rank-honest:

.. list-table:: The retired spellings and their rank-honest replacements
   :header-rows: 1
   :widths: 24 30 46

   * - Retired spelling
     - The phantom it carried
     - Replacement
   * - ``SNProblem.ny`` / ``SNProblem.dy``
     - At :math:`d = 1` these **lied** — ``ny`` returned a phantom
       ``1`` and ``dy`` a phantom ``[1.0]`` (the Issue #214 phantom
       class), and at :math:`d \ge 3` they underspecify the mesh.
     - :attr:`SNProblem.spatial_shape <orpheus.sn.problem.SNProblem.spatial_shape>` (the per-axis cell counts) and
       :attr:`SNProblem.axis_widths <orpheus.sn.problem.SNProblem.axis_widths>` (per-axis widths). ``AttributeError``
       on the retired names.
   * - ``SNProblem.dx``
     - A duplicate spelling of the per-axis widths.
     - :attr:`SNProblem.axis_widths <orpheus.sn.problem.SNProblem.axis_widths>` — promoted from the private
       ``_axis_widths`` to **the** single public spelling of per-axis
       cell widths.
   * - ``SNProblem.nx``
     - (kept) Documented :attr:`spatial_shape[0] <orpheus.sn.problem.SNProblem.spatial_shape>`
       sugar — honest at any :math:`d`, with a broad legitimate 1-D
       consumer base.
     - unchanged.

The phantom ``ny`` / ``dy`` at :math:`d = 1` is the **same defect
class** the N-D layout campaign closed at C2 / Issue #214: a trailing
singleton that masquerades as a real axis. The C5.2 retirement removes
the masquerade at the metadata source.

The two production ``dr`` consumers (the
:mod:`~orpheus.sn.loss_representation` 1-D bare sweep and the
:mod:`~orpheus.sn.angular.closure` Carlson preamble) repoint
from ``.dx`` to :attr:`SNProblem.axis_widths <orpheus.sn.problem.SNProblem.axis_widths>`. The
field / cross-section / scattering read-through chains collapse to the
rank-generic :attr:`spatial_shape <orpheus.sn.problem.SNProblem.spatial_shape>`:

* :class:`~orpheus.transport.fields.angular_flux.AngularFlux` (and the
  ``BulkField`` base) **retire** their ``nx`` / ``ny`` read-throughs.
  This is a live :math:`d = 3` landmine, not cosmetic: an
  ``(nx, ny)``-keyed field read **silently truncates** a 3-D tensor to
  its first two axes (a ``vv-principles`` Mode-2 / Mode-5 class
  index error that the degenerate :math:`d \le 2` test never reaches).
* :class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField` and
  :class:`~orpheus.transport.operators.scattering.ScatteringOperator` collapse their
  ``nx`` / ``ny`` reads to **one** rank-generic
  :attr:`spatial_shape <orpheus.sn.problem.SNProblem.spatial_shape>` read-through each.

Finally, a new :attr:`SNProblem.volume_measure <orpheus.sn.problem.SNProblem.volume_measure>` property gives the
SN-side ``keff`` rate consumers (the production / absorption rates in
:mod:`~orpheus.sn.solver`) a native source: they read it instead of
reaching through ``problem.mesh.volume_measure``. While the
:math:`d \le 2` adapter is present it delegates to the dataclass
measure (bit-identical, including the stored equal-volume cell
measures, which the ``precomputed_volumes`` hatch carried until P1 step
3b of #405 and the required ``Mesh1D.volumes`` field carries since);
the axis-native arm lands with C5.5.

.. _sn-c5-geometry-blind-trace:

The geometry-blind trace space (z faces admitted)
-------------------------------------------------

The trace layer (:ref:`bc-trace-structure`) carried two C5-blockers:
a hand-listed face-normal table that **silently lacked the z faces**,
and a gate-only ``mesh`` parameter on the trace factory.

**The face-normal source collapses onto** ``AXIS_NAMES``. Pre-C5.3,
``trace_space._FACE_NORMALS`` was a hand-listed **four-entry**
transcription (``xmin`` / ``xmax`` / ``ymin`` / ``ymax``) — it had no
``zmin`` / ``zmax`` rows, the silent :math:`d = 3` blocker. C5.3
derives the table from :data:`~orpheus.geometry.faces.AXIS_NAMES`
so every axis-aligned face (all **six** at :math:`d = 3`) is present by
construction: face ``"{axis}min"`` has outward normal
:math:`-\hat e_{\text{axis}}`, face ``"{axis}max"`` has :math:`+\hat
e_{\text{axis}}`. The :math:`\Omega\cdot\hat n` row for ``zmax`` is
then exactly :math:`+\mu_z` and for ``zmin`` exactly :math:`-\mu_z`.

To share the ``"{axis}{min|max}"`` crosswalk without an ``sn``-ward
import from the ``numerics`` layer, ``AXIS_NAMES`` **moved down** to
:mod:`orpheus.numerics.face_layout` — the home of
:class:`~orpheus.numerics.face_layout.FaceLayout`, keeper of the face
string-name world. :mod:`orpheus.mesh.axis` re-exports it, so SN consumers
are unchanged; the trace space (a ``numerics`` leaf) then read it
without depending on ``sn``. Since layer DAG WP-B (2026-10-22) the naming half
of that crosswalk (``AXIS_NAMES``, ``face_name``, ``face_normal``,
``face_opposite``) lives in :mod:`orpheus.geometry.faces`: a face's name and
normal are statements about the domain and need no measure, so geometry is
their home, and :mod:`~orpheus.numerics.face_layout` keeps only the
flat-buffer descriptor and the streaming-normal measure.

**The trace factory is geometry-blind.** The mesh parameter on the old
``from_mesh_and_quadrature`` factory was **gate-only** — its single use
was an ``isinstance`` check that refused a curvilinear ``Mesh2D``. That
refusal is **unreachable**: a curvilinear ``Mesh2D`` cannot become an
:class:`SNProblem` in the first place (2-D cylindrical SN has no sweep), so
no such mesh ever reached the factory. The ``isinstance`` check carried
no data the factory used. C5.3 therefore renames the factory to
:meth:`AngularTraceSpace.from_quadrature_and_layout
<orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace.from_quadrature_and_layout>`
and **retires the dead mesh parameter** (aggressive retirement; callers
and the bare-constructor error message migrated). Every datum now comes
from the quadrature (the :math:`\mu_x` / :math:`\mu_y` / :math:`\mu_z`
cosines) and the layout's face names (the axis-aligned normals implied
by the ``"{axis}{min|max}"`` convention).

With the gate gone, the trace is built **unconditionally** (the pre-C5.3
``isinstance`` gate excluded only the unconstructible 2-D cylindrical mesh): at
C5.3 inside ``SNProblem.realize_boundary_law``, today in the :class:`SNProblem`
constructor (``self._trace = AngularTraceSpace.from_quadrature_and_layout(...)``),
ahead of the binding's realization. :attr:`SNProblem.angular_trace` is typed
and documented as **always non-None**. The ``mesh`` parameter of
``SNMethodSpace.for_face`` became **optional metadata** — nothing in the
realizer chain read it, and an axis-native :class:`SNProblem` passed
``None``; its successor
:class:`~orpheus.sn.mesh.method_space.SNMethodSpace` is the plain pair
``(quadrature, trace)`` and holds no mesh.

.. note::

   **The fail-loud rank-mismatch raise (the C5.3 elegance-review
   carry, landed in C5.5).** A trace built for an axis-:math:`k` face
   on a quadrature that lacks ``mu_k`` previously **zero-padded** the
   :math:`\Omega\cdot\hat n` row to all-tangential — silently
   producing a face on which *no* ordinate is inflow or outflow. C5.5
   makes ``_build_omega_dot_n`` raise loudly on that rank mismatch:
   asking for a ``"zmax"`` face from a 2-D quadrature is a
   construction error, not a silent all-tangential face.

.. _sn-c5-windowing-gs-gate:

Windowing and Gauss–Seidel gate on genuine dimensionality (vv Mode 9)
---------------------------------------------------------------------

C5.4 is the **highest-risk edit of the campaign** — a textbook
``vv-principles`` **Mode 9** case (a splitting / optimization verified
only in a regime where the wrong gate is *accidentally* satisfied). Two
gates inside the SN source-iteration driver keyed on
``problem.reduced is None``:

* the **moment-windowing** gate (then ``_maybe_window``; since package
  S4 of the numerics untangling the solver's ``_angular_representation``,
  which returns the
  :class:`~orpheus.sn.operators.windowing.AngularMomentWindow`
  representation on the gate), which decides whether the SI iterate is held
  as compact harmonic moments rather than the full angular flux; and
* the **Gauss–Seidel splitting** selector — then the solver-private
  ``_select_si_splitting``, since 2026-09-13 the schedule resolution
  :func:`~orpheus.sn.splitting.resolve_schedule` that the
  :class:`~orpheus.sn.splitting.Splitting` value is labelled by
  (:ref:`sn-splitting-is-a-strategy-value`) —
  which decides whether the boundary-G-S accelerator is used.

``reduced is None`` is a **coincidence proxy**: it is ``None`` for
*every* multi-D Cartesian mesh, including a :math:`d = 3` one. The
moment-windowing path's in-sweep moment-emission kernel is **2-D
only** (it indexes a ``(N_oct, ng, nx, ny)`` block; see Issue #227).
So at :math:`d = 3` the old proxy would have **silently
moment-windowed the SI iterate** — the
:class:`~orpheus.sn.loss_representation.FullFieldWavefront` spine
refuses moment mode on the Jacobi path and ``None``-subscripts on the
G-S path: a **corrupted iterate, not a principled refusal**. This is
precisely the Mode-9 failure: a :math:`d \le 2` test cannot observe it
because at :math:`d \le 2` ``reduced is None`` *and* the 2-D kernel is
correct, so the proxy and the truth coincide.

C5.4 retargets both gates to the **genuine** dimensionality predicate:

.. list-table:: The Mode-9 gate retarget
   :header-rows: 1
   :widths: 34 30 36

   * - Gate
     - Pre-C5.4 (coincidence proxy)
     - Post-C5.4 (genuine condition)
   * - Moment windowing (then ``_maybe_window``; now
       ``_angular_representation``)
     - ``reduced is None``
     - ``is_cartesian and ndim == 2`` — the genuine
       windowing-eligibility condition (the 2-D moment kernel's exact
       domain).
   * - Boundary-G-S (then ``_select_si_splitting``; now
       :func:`~orpheus.sn.splitting.resolve_schedule`)
     - ``reduced is None``
     - ``is_cartesian and not is_1d`` — multi-D Cartesian.  The
       retargeted predicate moved verbatim with the gate when the
       splitting became a value.

The G-S resolvent's old ``"2-D Cartesian ONLY"`` docstring was **stale
Phase-3 narration**: :attr:`SweepSchedule.gauss_seidel
<orpheus.sn.loss_representation.sweep_schedule>` and the scheduled sweep
(``_sweep_scheduled``) have been **d-generic since C3**, so the
resolvent is constructible at :math:`d = 3`. The narration is corrected
in C5.4; the actual :math:`d = 3` boundary-G-S *fixed-point invariance*
(that G-S and Jacobi converge to the **same bulk** flux) is
**value-gated**
by the C5.5 Mode-9 mixed-BC box (:ref:`sn-c5-value-gates`) before any
:math:`d = 3` G-S solve is trusted — the Mode-9 discipline made
operative: never gate a splitting's FP-invariance on a degenerate box.

.. note::

   ⚠ **The word "bulk" above is load-bearing, added 2026-08-15 (#344),
   and this gate does not need re-scoping — but its measurands do need
   naming.**  The C5.5 box is x-reflective / y-vacuum / z-reflective,
   i.e. it closes **two** reflective axis pairs, and ``[M]`` that makes
   :math:`A = L+C-S-N_{2n}-B` **exactly singular** there:
   :math:`\dim\ker A =
   36` at cells :math:`(5,3,4)`, ``level_symmetric`` :math:`S_4`,
   :math:`n_g = 2` (:math:`= n_g (N/4)\, n_y`).  So G-S and Jacobi do
   **not** return the same boundary trace on it.  The gate is sound
   because what it asserts — :math:`k_{\rm eff}` and the normalized
   flux *shape* — is mirror-even, and every mirror-even functional
   annihilates :math:`\ker A` by theorem
   (:eq:`sn-kernel-mirror-blindness`).  ⟹ a future strengthening of
   this gate must not reach for the raw trace without gauging it
   first.  Derivation: :ref:`sn-loss-kernel-gauge`.

The one ``reduced is not None`` branch that **stays** is the 1-D
sweep-cache: there the predicate keys on the *availability of the data
it reads* (the reduced-operator cache), not on dimensionality. That is
not a proxy — it is the genuine guard.

.. _sn-c5-3d-admission:

3-D Cartesian admission: the axes tuple is the only 3-D entry
-------------------------------------------------------------

After the C5.1–C5.4 cleanup, the :math:`d = 3` admission is an
**extension, not a new arm**. A 3-axis Cartesian :class:`SNProblem` now
constructs and **solves** through the same generic body as
:math:`d \le 2`, **mesh-adapter-free from birth** (``self.mesh is
None``) on the d-generic
:class:`~orpheus.sn.loss_representation.FullFieldWavefront` spine.

* **The gate retires.** :meth:`SNProblem.from_axes` drops the
  ``d \ge 3`` admission guard; :math:`d \le 2` still synthesizes the
  legacy adapter for its remaining consumers.
* **Axis-native arms.** The cell-volume array is the iterated outer
  product of the per-axis widths,
  :math:`V[i,j,k] = \mathrm{d}x_i\,\mathrm{d}y_j\,\mathrm{d}z_k`; the
  :attr:`volume_measure <orpheus.sn.problem.SNProblem.volume_measure>` is the rank-:math:`d`
  meshgrid-of-centers
  :class:`~orpheus.numerics.measure.DiscreteMeasure` (the natural
  rank-:math:`d` generalization of the ``Mesh2D`` analogue).
* **Entry surface.** :func:`~orpheus.sn.solver.solve_sn` and
  :func:`~orpheus.sn.solver.solve_sn_fixed_source` accept the **axes
  tuple** — the *only* 3-D entry — through one inbound seam
  (``_as_problem``). A new ``mat_map`` keyword is the axes-entry
  material channel (it raises if combined with a legacy mesh, which
  carries its own material map). The entry takes every face's role from
  the declaration, on both surfaces: each axis declares a role (a deck or a
  tag) on each endpoint, the exterior conditions come from the
  ``boundary_conditions`` keyword, and no entry fills a face
  (:ref:`structured-geometry-no-default-law`).

.. note::

   **Every face is declared.** The value gates below declare their
   faces: every face a mirror deck for the headline
   :math:`k_\infty` identity and the closed-form flux identities, and a
   mixed declaration (mirror decks and a vacuum tag) for the Mode-9 box.
   ``test_axes_declared_laws_are_taken_verbatim`` checks that the
   entry realises a mixed 3-D declaration face by face as declared.
   Until P1 step 3c of #405 (2026-09-30) an undeclared face had two
   readings, vacuum under a fixed-source entry's fill and reflective
   under the resolver's default; both retired.

.. _sn-c5-value-gates:

Numerical evidence — the d=3 value gates
----------------------------------------

C5.5's admission is gated by four value tests
(:mod:`tests.gates.sn.solve.test_d3_admission`), **all driven through the
production entry points** (``np.testing.assert_*`` only — Mode-8 safe
under ``python -O``, where bare ``assert`` is stripped). Each probes a
distinct failure class:

.. list-table:: The d=3 admission value gates
   :header-rows: 1
   :widths: 30 16 54

   * - Gate
     - V&V level
     - Evidence
   * - **k_inf 3-D ≡ 2-D ≡ 1-D**
     - L1 (closed-form eigenvalue)
     - Homogeneous all-reflective boxes at :math:`d = 1, 2, 3`. The
       reference is the closed-form matrix eigenvalue
       :math:`k_\infty = \lambda_{\max}(A^{-1}F)`,
       :math:`A = \operatorname{diag}(\Sigma_t) - \Sigma_{s0}^{\mathsf T}`
       — **never** the sweep. Each dimension matches ``case.k_inf`` to
       ``atol=1e-8``; the d=3 box solved
       :math:`1.8750000050` against the closed form :math:`1.875`
       (2g). Run at **2 groups and 4 groups, never 1 group** — a 1-G
       eigenvalue is the flux-shape-independent ratio
       :math:`\nu\Sigma_f/\Sigma_a` and is degenerate.
   * - **Per-ordinate ψ = Q/(W·Σₜ)**
     - L1 (closed-form flux)
     - Pure absorber (:math:`c = 0`), all-reflective. DD is flat-flux
       *exact* and :math:`c = 0` needs no iteration, so **every**
       ordinate must carry the closed-form value
       :math:`\psi_{n,g} = Q_g/(W\,\Sigma_{t,g})` to ``rtol=1e-10``.
       Per-group distinct :math:`Q` and :math:`\Sigma_t` make a group
       swap (Mode-2) observable; the per-ordinate residual is the
       sharpest Mode-1 / Mode-3 / Mode-4 probe.
   * - **Scattering multigroup balance**
     - L1 (closed-form flux)
     - Scattering medium, all-reflective:
       :math:`\phi = (\operatorname{diag}(\Sigma_t) -
       \Sigma_{s0}^{\mathsf T})^{-1} Q`. The group-coupling companion —
       a **Mode-6 convention-drift catcher** because mixture C's
       scattering matrix is **asymmetric**, so :math:`\Sigma_s` vs
       :math:`\Sigma_s^{\mathsf T}` (the ``SigS`` / ``SigS^T``
       convention, see :ref:`theory-discrete-ordinates`) is observable.
       Measured max relative error :math:`2.6\times 10^{-9}` — this is
       **SI-convergence-limited, not a discretization error** (DD is
       flat-flux exact on a homogeneous box).
   * - **Mode-9 G-S ≡ Jacobi FP-invariance** (of the **bulk** — the
       box is kernel-BEARING, see the note above)
     - L2 (integration)
     - Boundary-Gauss–Seidel and Jacobi converge to the **same** d=3
       bulk on a box that **breaks every degenerate
       coincidence**: mixed BCs (x-reflective / y-vacuum / z-reflective
       — axis-asymmetric, so a wrong reflection partner shifts the
       answer), ``nx ≠ ny ≠ nz`` (5, 3, 4), a heterogeneous 2-G split
       across x (a non-flat-flux guard), and a **diagonal**
       level-symmetric cubature (ERR-056 shared-face discipline —
       diagonal cubatures share faces between octants, the regime where
       a wrong G-S shared-face reflect is observable). :math:`k_{\rm
       eff}` agrees to ``atol=1e-8`` *and* the normalized flux shape to
       ``rtol=1e-6``.

The four gates together cover the verification ladder for the new
capability: a closed-form eigenvalue (L1, the only pillar that can
verify :math:`k`), two closed-form flux identities (L1, isolating the
streaming / collision / scattering operators per-term), and a
splitting FP-invariance on the degenerate-breaking box (L2 / Mode-9).
The eigenvalue gate's reference is **structurally independent** of the
sweep (a matrix eigensolve, not a transport solve), satisfying the
``vv-principles`` requirement that an eigenvalue claim rest on a
closed-form or semi-analytical reference rather than MMS.

What runs the d=3 path, and what is deferred
--------------------------------------------

The :math:`d = 3` admission runs on the **d-generic
FullFieldWavefront ORACLE spine** — the never-stuck full-field
representation that is correct from day one (the four value gates), but
**not** the optimized sweep kernels. Two kernel widenings are deferred
to Issue #227, gated on *measurement* against the spine (the C3.6
principle: "construct general, select narrow, specialize only on
measured cost"):

* **ScanMarch** :math:`d \ge 3` — the row-march kernels currently
  unpack 2-D pairs; the
  :math:`\text{scan}(x)\circ\text{march}(y, z)` generalization widens
  the predicate **only with** the kernel and a profile showing it beats
  the spine.
* **MovingFrontierWindow** :math:`d \ge 3` — the rolling-frontier
  window is built ``frontier_dim = d-1`` and its ``supports`` is
  conservatively ``is_cartesian and ndim == 2``; the :math:`d = 3`
  windowed *walk* is graph-layer-pinned but the window kernels need
  their own profile (the 2-D window was a ~0.71–0.80× **speedup** plus
  a peak-memory win — the :math:`d = 3` economics need separate
  numbers).

Separately, the **multi-D adjoint** (``loss_action_transpose`` raised
:class:`NotImplementedError` at any multi-D) was a **pre-existing
deferral** orthogonal to C5 (G-adjoint campaign territory), noted here
only so the :math:`d = 3` capability map was complete. Its two arms
differ and are discriminated. The transpose **matvec**
(``loss_action_transpose``) is lifted: since #310 C4/C5 every registered
loss representation has a reverse walk (the ``has_transpose_walk``
docstring: every registered leaf overrides it to true), and the 2-D adjoint
row is no longer excluded (:doc:`/theory/methods/sn/loss_representation`
carries the account). At :math:`d = 3` it is gated by
``tests/gates/sn/sweep/core/test_multi_d_reverse_walk.py::test_d3_dense_mt_and_pairing_on_the_spine``.
``[M]`` 2026-10-24, ``python -O -m pytest`` on that node id from this
worktree: ``1 passed, 1 warning in 5.47s``. The transpose **solve**
(``has_transpose_sweep``, the sweep :math:`S^{H}`) of the multi-D wavefront
is still deferred: the base class returns ``False`` and the multi-D
representations (``FullFieldWavefront``, ``MovingFrontierWindow``,
``ScanMarch``) do not override it; only the 1-D scan realises it. The four
gates above are forward solves.

C5 closure
----------

C5 lands in six substeps under Issue #225 (the SN N-D layout campaign):
C5.1 (axis-primary inversion), C5.2 (phantom-shim retirement +
native ``coord`` / ``volume_measure``), C5.3 (geometry-blind trace),
C5.4 (Mode-9 windowing / G-S gate retarget), C5.5 (3-D admission), and
C5.6 (structure-pin flips to the now-constructible mesh). The
:math:`d \le 2` path is **byte-identical** on every numerical output
(the affine ``sha256`` goldens are unchanged across the whole carve);
the :math:`d = 3` path is correct from birth on the FullFieldWavefront
spine, value-gated by the four tests above. The campaign reaches its
3-D admission **without** a ``Mesh3D`` dataclass — the axes tuple,
made primary by C5.1, *is* the N-D entry.


Anti-pattern catalog
====================

A short list of patterns the refactor's authors considered and
rejected, with the reasoning preserved so future sessions don't
re-attempt them. Each entry says what, if anything, the wave-2 types
(2026-10-24) made of it; a retired name is a literal.

1. **Single ``BoundaryOperator`` ABC carrying both law and
   realizer responsibilities.** This is the pre-refactor shape.
   Rejected because the law / realizer split is the architectural
   point of the refactor: laws are method-agnostic, realizers are
   method-specific. Keeping them in one class would force every
   law to know about the SN sweep's half-trace plumbing (an inflow
   mask at the time; the :math:`\gamma_\pm` restrictions since B3.2),
   which is precisely what Waves 0 to 2 of the twelve-wave
   refactor were supposed to abstract away. The split survives wave 2 in
   other objects: a boundary condition (a
   :class:`~orpheus.transport.boundary_conditions.WallCondition`, or a deck)
   is method-agnostic data, and each method's realizer table
   (``SN_REALIZERS``, ``DIFFUSION_REALIZERS``, ``CP_REALIZERS``,
   ``MOC_REALIZERS``, ``MC_REALIZERS``) realizes it, so a condition cannot
   know a sweep's plumbing.
2. **Dedicated ``MixedBoundaryOperator`` class.** Pre-Wave-11
   shape. Rejected (deleted Wave 11) — the Wave-0 algebra dunders
   on every ``BoundaryTraceLaw`` already produced
   ``OperatorSum`` shapes; the dedicated class added a
   second realization path with no semantic difference. See
   :ref:`bc-rank-n-algebra`. Under wave 2 the sum of walls is
   :class:`~orpheus.transport.boundary_conditions.WallCondition`'s own ``+`` and
   ``*`` (``WallCondition.from_amplitudes(0.3, 0.7)`` spells the common pair),
   and its realization is a ``RealizedSum``: one mixed-wall type at each
   layer, no dedicated class.
3. **Shared ``BoundaryRealizerBase`` ABC for cross-method
   realizers.** Considered Wave 5, when only one functional
   realizer existed. Rejected per the "Unify after two instances"
   architectural discipline: building the abstraction on a single
   instance would force a particular shape on every future method
   based on SN's dispatch idiom (``isinstance``). Vindicated at
   method #2: the diffusion realizer (#290 P3) chose a DIFFERENT
   shape (law → albedo scalar → structure-keyed collapse, not an
   isinstance ladder over per-law primitives), and the structural
   ``BoundaryRealizer`` Protocol stayed the only shared contract through
   wave 1. Wave 2 then minted the shared object, after the diffusion shape had
   shown what SN's dispatch idiom lacked: the
   :class:`~orpheus.transport.boundary_conditions.RealizerTable`, one per
   method (five tables exist), which answers admission and realization
   together and, after one split by declaration kind (deck, source, condition),
   keys its realizers by the unit realized, not by an
   ``isinstance`` ladder over the laws. The lesson stands in the
   corrected form: the abstraction was extracted when the instances had shown
   what they share (a table keyed by the unit, a refusal that names the face),
   and it does not have SN's shape.
4. **Adding ``face`` to ``VacuumInflow``'s constructor for
   semantic correctness on the standalone-apply path.** Option (b)
   in :ref:`bc-vacuum-semantic-correction`. Rejected because it
   would have inflated every test signature for one wave to fix a
   path that the refactor was retiring anyway. The transitional
   legacy-vacuum body was the right cost. Both the path and the constructor are
   gone: vacuum is ``WallCondition.vacuum()``, a condition with no face
   parameter, and a realized vacuum is the ``ZeroMap`` of the face it is
   realized on, whose domain and codomain are read off that face's trace
   ladder, so a vacuum built against the wrong half-trace is unspellable
   (ERR-041 is the catalogued case).
5. **Auto-importing every cross-method realizer at
   ``orpheus.geometry.boundary`` import time.** Considered in
   the registry era to make ``BoundaryRealizerRegistry.get("MoC")``
   work without the caller having to ``import orpheus.moc`` first.
   Rejected because :mod:`orpheus.sn` is a heavy module that every
   consumer of the boundary package would then pay for. #290 P7b
   made the whole question moot by dissolving the registry: each
   method-mesh imports its own realizer explicitly, so there is no
   name-lookup to keep populated and no import-side-effect timing
   to defend.
6. **Cartesian-vs-curvilinear bypass in
   ``SNProblem._resolve_one`` + dual-mode shim.** Pre Issue #188
   shape: curvilinear ``Mesh1D`` bypassed the realizer and wrapped
   the bare 2-arg law in
   ``_BoundBoundaryOperator(law, quadrature=self.quad)``, while
   Cartesian routed through the realizer with
   ``_BoundBoundaryOperator(realized)``. **Retired Issue #188 +
   #176**; documented here because the seductive trap is to
   "preserve flexibility" by keeping the dual mode after the
   curvilinear deferral lifts. The right move is to delete the
   bypass and consolidate on one path — see
   :ref:`bc-curvilinear-realizer-unification`. Wave 2 retired the shim too; the one
   path is ``SNProblem.realize_boundary_role`` into ``SN_REALIZERS``.
7. **Option A (keyword-optional ``quadrature=None`` on the
   concrete laws).** Landed Issue #176 / C176.3 as the interim
   form; **retired Issue #186 / B3 + β2** in favour of the
   pure-descriptor model (no ``apply`` on any law). The
   architectural costs of Option A (asymmetric semantics on
   ``quadrature=None``, vacuum two-paths-divergence, Liskov
   violation under polymorphic typing) made it unsustainable as
   the long-term contract; the interim was kept only long enough
   to land curvilinear realizer unification (Issue #188 first)
   before the descriptor cleanup could ship. See
   :ref:`bc-trace-law-descriptor-model` for the full retrospective.
8. **Calling ``apply`` on a raw BC declaration.** Under the
   pre-#186 contract this either worked (with surprising
   semantics — see Option A entry above) or raised
   ``BoundaryError``. Under post-#186 it was a **static type
   error** — ``BoundaryTraceLaw`` had no ``apply``
   method on the class, and neither did ``LawSum`` /
   ``LawScaled``. Wave 2 keeps the split in other types: a
   :class:`~orpheus.transport.boundary_conditions.WallCondition` and the decks
   are data and have no ``apply``, and the realized operator has. The correct
   contract is the problem's own door, ``problem.bc[face].apply(psi)``
   (the face's :class:`~orpheus.sn.boundary.realized_operators.RealizedBoundary`,
   on a :math:`\Gamma_+`-shaped argument), or, for a face outside a problem,
   ``SN_REALIZERS.realize(declared, site, method_space).apply(psi)`` (the call
   ``tests/gates/sn/architecture/test_wave2_carry.py`` makes as
   ``H.sn_table().realize(declared, where, space)``). The table is the **sole**
   bridge; the §16A.3 three-layer split is enforced by the type system.
9. **In-tree Wave-0 operator algebra over unrealized
   ``BoundaryTraceLaw`` instances (β1 form).** Considered as
   the rank-N composition mechanism during Issue #186 B3
   exploration. Rejected in favour of the separate-type-family
   approach (β2 / ``LawSum`` / ``LawScaled``). β1
   produced ``OperatorSum`` trees whose leaves were laws,
   not operators — the type checker could not distinguish a
   not-yet-realized "operator" from a real operator, and calling
   ``apply`` on the β1 tree raised at the leaf realization
   step. β2 made the law-vs-operator distinction inspectable
   statically: ``LawSum`` had no ``apply`` method, full
   stop. See :ref:`bc-rank-n-algebra` for the detailed
   comparison. Wave 2 replaced the β2 family by one declared type: a sum of
   walls is a :class:`~orpheus.transport.boundary_conditions.WallCondition`
   (its ``+`` and ``*``; a combination that leaves the cone of walls is a
   ``WallDisplacement``, which is never bindable), data with no ``apply``,
   and its realization is a ``RealizedSum``. The law-vs-operator distinction
   is still inspectable statically, and a sum is one type, not a type family.
10. **A string registry beside the realizer, answering admission separately.**
    Until wave 2 each method kept a registry of admitted tags
    (``BOUNDARY_OPERATOR_REGISTRY``, ``BC_REGISTRY``) next to the realizer that
    built the operator: the SN registry admitted two kinds while the realizer
    realized six laws, so the query and the realization could disagree.
    Rejected: the realizer table
    (:class:`~orpheus.transport.boundary_conditions.RealizerTable`) answers
    "is this declaration admitted on this face" and "realize it" from one set
    of registrations, and a refusal is raised at construction naming the face
    (``BoundaryRefused``). The registries are retired
    (see :ref:`bc-wave2-crosswalk`).


References
==========

* Grand Report v3 §16, §16A.1–5 (affine boundary form + trace
  structure), §16A.10 (sparse trace primitives), §16A.11 (dual
  registry), §16A.12 + §27.6 (universal invariants), §26A.4
  (named-error catalog). Source: ``.claude/plans/neutron_transport_grand_report_v3.md``.
* The 12-wave refactor plan:
  ``.claude/plans/transient-giggling-cake.md``.
* The post-Wave-12 cleanup plan (Issue #188 + #176):
  ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md``.
* Lewis, E. E. & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*. American Nuclear Society. §3.4 (boundary
  conditions in transport).
* Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.5 (albedo, white, and Marshak
  boundary conditions).
* The tensor decomposition equation :eq:`bc-tensor-decomposition`
  at :ref:`bc-tensor-decompositions` (in
  :doc:`/theory/methods/sn/boundary_conditions`) shows the algebra
  :math:`B = \sum_\alpha G_\alpha \otimes A_\alpha` that the older sections of this page refine into the trace-law affine form
  :eq:`affine-bc-form`; the form of record is the pair :math:`(R, q)` of :ref:`bc-aor-trace-closures`.
* :ref:`operator-algebra` for the Wave-0 primitives the realized
  BCs decompose into.
* The V&V error catalog in the ``vv-principles`` skill
  (``docs/theory/verification/error_catalog.rst``) carries
  the ERR-040..ERR-047 entries in canonical form.

The sources of :ref:`bc-architecture-of-record`:

* Wave 2 of the layer DAG, the specification of record: ``.claude/plans/wave2_boundary_carve.md`` (section 1, the mathematics;
  section 2, the types; section 3, the laws and their gates; section 4, the admission tables; section 7.2, the retirement's
  dependency table; section 11, the records of the slices), and the rulings D1 to D5 and the thirteenth exchange of
  ``.claude/plans/boundary_law_ontology.md``.
* Coxeter, H. S. M. (1969). *Introduction to Geometry*, 2nd ed., Wiley, chapters 3 and 7: the classification of the non-identity
  isometries of :math:`E^n`, :math:`n \le 3`, into six types (translation, reflection, rotation, glide reflection, screw motion,
  rotoreflection), as the code states and uses it (``orpheus/geometry/deck.py``, ``IsometryClass``).
* Marshak, R. E. (1947). "Note on the spherical harmonic method as applied to the Milne problem for a sphere". *Phys. Rev.* 71,
  443; and Bell, G. I. & Glasstone, S. (1970), §3.4 (diffusion-theory boundary conditions): the :math:`P_1` closure that the
  diffusion table's tier statements refer to, as cited in ``orpheus/diffusion/boundary_realizer.py``.
* ``.claude/plans/diffusion_crosswalk.md``: the :math:`(J^+, J^-)` convention contract of the scalar trace.
* The operator-as-tensor-network reading of the response :math:`R`: :doc:`/theory/foundations/operator_tensor_network`.
* The issues the admission cells cite, each read as ``OPEN`` by ``gh issue list --state all --limit 700 --json number,state`` on 2026-10-24: #511 (the hollow shell's inner boundary law
  is dropped), #513 (CP, MoC and MC ignore a slab's left boundary law), #514 (MoC ignores the mesh's coordinate system),
  #546 (CP's vacuum eigenvalue omits leakage), #547 (MoC and MC realize their law on a square cell), #650 (the reference
  vocabulary), #653 (the SN radial corner action), #654 (the CP, MoC and MC tables), #655 (MC's pitch), #656 (diffusion's inflow
  channel), #657 (``RigidMotion`` reads its own reflection plane), #658 (the 2-D colouring group), #659 (the SN and diffusion
  unregistered units), #660 (``ProblemSpec.boundary_kinds``) and #661 (``DiffusionMethodSpace.face``): 15 of 15 OPEN.


.. _bc-history:

Development history
===================

Reverse-chronological changelog of the boundary architecture's
milestones recorded since this table was opened (the earlier campaign
phases, B1 to B3.4c and the Wave refactors, are recorded at the sections
they concern). Entries marked *(in development)* live on an unmerged
branch; trust ``git`` over this table for merge status.

.. list-table::
   :header-rows: 1
   :widths: 10 58 8 24

   * - When
     - Milestone
     - Issue
     - Where
   * - 2026-10-22 to 2026-10-24
     - **A boundary condition is a face role, one wall condition on the sub-Markov simplex, a binding and a realizer table (wave 2 of the layer
       DAG, workflow WF2.5 documents it).** The trace-law architecture, ``BoundaryTraceLaw`` and its seven concretes chosen by
       string, the tag ``BC``, the descriptor algebra ``LawSum`` and ``LawScaled``, the ``BoundaryRealizer`` Protocol and the
       string registries, was replaced by the architecture of :ref:`bc-architecture-of-record`; every retired symbol and its
       successor is in :ref:`bc-wave2-crosswalk`. What landed, from the plan of record (``.claude/plans/wave2_boundary_carve.md``,
       section 11):

       * *The leaves (WF2.2, additive, integrated 2026-10-23).* L1, the deck family, the colouring symmetry and the one
         coordinate-agreement convention (``487d2d19``; amendments ``8a3e6c16`` and ``f6690aeb``). L4, ``ordinate_pairing`` and
         the ``TraceRole`` ladder (``21a5e0e7``; ``afc8ed6d``). L5, ``WallAmplitudes``, the reference's plain-data wall
         (``40aa15b5``).
       * *The carve (WF2.3, sequential on one branch, no dual spelling).* S1, face roles in the geometry and the mesh
         (``6b090ff9``, ``5b7f1dee``, ``4ec95bc6``). S2a, the conditions package (``f488034c``, ``7dcb9384``, ``cefc8339``). S2b,
         the SN consumers and the typed realized operators (``a90d10c2``, ``be353f8c``, ``0757a8d9``). S3, diffusion
         (``9c0e8fde``, ``508c6478``, ``d4da4328``). S4, CP, MoC and MC (``0f5e1c5d``, ``2a7ea40d``, ``de354496``). S5, the
         references and the specification state their exterior as plain data (``25bae6a3``, ``96b5c675``).
       * *The retirement (WF2.4).* ``ece6a281`` deleted ``orpheus/geometry/boundary/`` (15 modules),
         ``orpheus/transport/boundary_realization.py`` and ``orpheus/sn/boundary/bound_compat.py``, carries ``Closes #551`` and
         ``Closes #189``; `[M]` ``git show --shortstat ece6a281``: 96 files changed, 421 insertions and 12957 deletions. The amendments
         are ``c2a14005`` and ``0707a332``, and ``8f2cefa4`` made Sphinx ``-W`` green at that tip.
       * *The rulings that bind.* D1: a general condition is the type and a conventional name a facade over it. D2: one
         normalised combination type, non-negative, total return at most one, signed combinations unbindable. D3: today's
         limits stand and the reference states ``WallAmplitudes`` only (#650). D4: no masking, so CP's vacuum is realized and
         #546 is a strict ``xfail``. D5: the hollow inner face is a tag. The standing rulings of the plan's section 0: the
         thirteenth exchange (one shared exterior type whose response is a tensor network of typed factors; MC on a cylinder is
         a strict ``xfail`` naming #547), plain data at the references edge with one translation that imports neither the
         specification nor the derivations, and a sequential carve with no transit exemption on the qualname gate.
       * *Departures from the draft, each ruled at integration.* The uniform-inflow facade is ``WallCondition.uniform``, not
         ``inflow``; ``Admission`` is a sealed class and not a ``TypeAlias``; no ``OrdinatePairing`` carrier (the type checks its
         own facts); the gauge reads the realized operator and not ``binding.decks``; diffusion and SN refuse the hollow inner face
         except a unit specular return (#511), and CP refuses it (#513); the MC meets follow ``CHAIN_ORDER``.
       * *State at the retirement tip.* `[M]` the plan's record (``python -O -m pytest -m "not slow"``, serial, per directory):
         0 failed and 0 collection errors in 18 of 18 directories of ``tests/gates``; Sphinx ``-W`` exit 0 at ``8f2cefa4``; the
         {geometry, numerics} package cycle gone from the layer census.
     - #551, #189 (``Closes`` trailers of ``ece6a281``); #546, #547, #650, #653, #654, #655, #656, #658, #659, #660, #661 cited
     - ``refactor/wave2-boundary-carve``; ``ece6a281`` (the retirement)
   * - 2026-10-20
     - **The walk is a homomorphism into the operator algebra, and the
       S\ :sub:`N` amplitude law has no branch.** With package S2 of the
       numerics untangling, ``+`` and ``α·`` build the normalised formal
       linear combination, so ``realize_recursively``
       maps ``LawSum`` to ``+`` and ``LawScaled`` to ``α·`` and no longer
       keeps the law tree's shape (``0.5·(0.3X + 0.7Y)`` realizes as
       ``0.15X + 0.35Y``). The S\ :sub:`N` realizer's three α-regimes
       became one expression, ``α · (K_ω ⊗ I)``: :math:`0\cdot X` is the
       zero map and :math:`1\cdot X` is :math:`X` (`[M]` 1619 of 1619
       calls). ``ScaledOperator``, whose zero-scalar refusal forced the
       :math:`\alpha = 0` branch, retired. See
       :ref:`bc-realize-recursively`.
     - —
     - ``[LANDED 4e5ee36b]``
   * - 2026-10-01
     - **The mirror carries no amplitude, and a deck law is not an
       operand.** ``ReflectiveBoundary`` lost its ``albedo`` parameter,
       its ``kind`` override and its two amplitude checks; it is now
       the parameter-free deck mirror, :math:`G` the mirror about
       ``axis`` and :math:`R` = ``ScalarResponse(1.0)``. Before this,
       ``ReflectiveBoundary(axis, α)`` with :math:`\alpha \ne 1` had
       both factors non-trivial (the one shipped row that violated the
       criterion of :ref:`bc-factor-quotients`): it was the partial
       wall ``AlbedoBoundary(α, SpecularReturn(axis))`` spelled a second
       time in the geometry slot, the fossil of the 2026-05
       tensor-decomposition framing, whose amplitude slot it filled.
       Its ``kind`` reported ``"partial"`` at :math:`\alpha \ne 1` while
       its registry key stayed ``"reflective"``, which is why the
       realized shim read the key and not ``kind``: reading ``kind``
       would have dropped such faces from the sweep schedule's
       reflective set. S\ :sub:`N` read the attenuated mirror as a
       perfect one in three places (ERR-094, whose review on
       2026-09-30 added the two amplitude checks ``assert_submarkov``
       and ``assert_response_positive_if_declared`` to the class; they
       retired with the field). The tag parser always built it at
       :math:`\alpha = 1`, so only a typed declaration reached it.
       ``LawSum`` and
       ``LawScaled`` now refuse a deck
       law as a direct child (:ref:`bc-deck-law-not-an-operand`), on
       the new ``SelfPairedDeck.is_identity`` / ``PairedDeck.is_identity``
       read off ``RigidMotion.is_identity``; the refusal is tagged
       ``ELEGANCE-DEBT[guard] #551``. Every documented composition
       example and the rank-N Marshak example now take
       ``AlbedoBoundary(1.0, SpecularReturn("x"))`` as the specular leaf.
       Gates: ``tests/gates/geometry/test_reflective_is_a_mirror.py``
       (36 rows) and ``tests/gates/geometry/test_deck_laws_do_not_compose.py``
       (217 rows). `[M]` 2026-10-01, ``.venv/bin/python -O -m pytest``
       on the two files: before the carve 12 of 36 and 67 of 217
       failed, after it 253 of 253 pass. A carve-time carry gate (kept
       in ``scratch/``, not in the suite) showed the frozen arrays of
       the partial wall (75) and of the mirror (22) bit-identical
       across the respelling, because ``AlbedoBoundary(α,
       SpecularReturn)`` and the old ``ReflectiveBoundary(axis, α)``
       already shared one realization body. Ruling and plan of record:
       ``.claude/plans/boundary_law_ontology.md``, the fourth exchange
       (why the albedo existed) and the tenth (this cleanup).
     - #551
     - ``65f1dc7a`` (code) and ``99b9847d`` (documentation) on ``main``
