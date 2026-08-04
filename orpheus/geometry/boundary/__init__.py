r"""Boundary conditions for transport solvers — trace-law architecture.

This package ships the post-12-wave **trace-law / realizer / tensor-
algebra** architecture for transport-equation boundary conditions
described by Grand Report v3 §16, §16A, and §16A.10. The narrative
theory page is :doc:`/theory/foundations/boundary_conditions`; the design plan
is ``.claude/plans/transient-giggling-cake.md`` (the 12-wave
"transient-giggling-cake" effort).


Architecture: three layers, two registries
==========================================

A boundary condition on a transport interface is, in the
method-agnostic sense, a single affine map on the boundary trace

.. math::

    \gamma_- \psi \;=\; R\,G\,\gamma_+ \psi \;+\; q,
    :label: bc-affine-form-init

where :math:`\gamma_\pm` are the **inflow / outflow trace operators**
that restrict :math:`\psi(\mathbf{r}, \Omega)` to the directional
halves of :math:`\partial\Omega`; :math:`G : \Gamma_+ \to \Gamma_-`
is the **deck transformation** (a mirror, a spatial wrap, a rotation —
the composition operator of a measure-preserving bijection of the
boundary phase space, which is why it carries the CROSSING: the mirror
that exchanges the two hemispheres is an ambient isometry);
:math:`R : \Gamma_- \to \Gamma_-` is the **constitutive response** (a
scalar amplitude in :math:`[0, 1]` for sub-Markov BCs, or a rank-one
angular kernel for diffuse re-emission); and :math:`q \in \Gamma_-` is
the optional **prescribed inflow source** (zero for the homogeneous
case).

Membership in :math:`G` takes **two** tests. Multiplicativity —
:math:`G(\psi\varphi) = (G\psi)(G\varphi)` — is **necessary**: it holds
for a relabeling and never for an average, so anything failing it is a
kernel and belongs to :math:`R`. That is what campaign phase **B3.0**
used to move the Lambertian out of the geometry slot. It is **not
sufficient**: a specular *kernel* is a permutation, hence multiplicative,
and is still constitutive. The sufficient test is the quotient one —
:math:`G` is the deck transformation of an **actual quotient of the
physical domain**, and a wall standing in the domain is not one (**B3.4b**,
user ruling 2026-08-01). Whence **exactly one of** :math:`G`, :math:`R`
**is non-trivial**. See :ref:`bc-factor-roles` on the theory page and the
``_factors`` module docstring.

Which methods can realize which laws — three tiers and three axes
=================================================================

.. note::

   The referenceable anchor for this taxonomy is ``bc-method-realizability`` on
   the boundary theory page, which owns it. This module is not autodoc'd, so a
   label defined here would never register — and would collide the day the
   package is added to ``docs/api/geometry.rst``. What follows is the code-side
   statement of the same content, sited here because it is what a reader of the
   realizers actually has open.

A transport method is a **discretization of the trace**, i.e. a
projection :math:`\Pi : \Gamma \to \Gamma_h` (SN keeps ordinates, P1
the half-range moments, MoC track angles, MC nothing at all). Whether a
law realizes in a method is whether the naturality square commutes:
:math:`\Pi \circ (R\,G) = (R_h\,G_h) \circ \Pi`. Three tiers:

1. **Exact and faithful** — :math:`R = \alpha I` for **any**
   :math:`\alpha`, not merely 0 or 1. Scalars commute with every linear
   projection and :math:`\alpha` is recoverable downstairs. This is why
   :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
   is one line.
2. **Exact but NOT faithful** — the square commutes, the realization is
   correct, and the projection nonetheless **identifies laws that differ
   upstairs**. At P1 a specular mirror and a Lambertian average both give
   :math:`J^- = \alpha J^+`: neither is approximated, P1 simply cannot
   tell them apart. This is the same fact as "on a scalar trace
   :math:`G` is forced to the identity by dimension", read from the
   other side — and it is why a rank-one :math:`R` makes :math:`G`
   unobservable.
3. **Not exact** — the law's action depends on structure the method does
   not represent (an anisotropic kernel below P1).

So the dividing line is **scalar vs angular**, NOT trivial vs
non-trivial.

:math:`G` is method-independent *as a geometric object* but its
realization is conditional on **equivariance**, and that condition
splits by which coordinate the map touches. A specular mirror acts on
the ANGULAR coordinate, so :math:`G_h` exists only if the quadrature is
symmetric under the reflection — which is exactly what ERR-042
(measure-preserving), ERR-044 (involutive) and ERR-045
(inflow :math:`\to` outflow) verify: they are **discretization-admits-
the-symmetry** checks, not physics checks. A spatial wrap acts only on
the SPATIAL coordinate, so every angular discretization is trivially
equivariant under it — periodic is the more method-agnostic of the two,
for a sharper reason than "it is a trace connection".

**Three INDEPENDENT axes** decide a refusal. Each realizer's guards name
their axis explicitly, because reading them as one axis is what made
this taxonomy invisible for so long:

.. list-table::
   :header-rows: 1
   :widths: 22 40 38

   * - axis
     - question
     - the refusal that shows it
   * - angular resolution
     - can the method represent :math:`R`'s angular structure?
     - tier 3 — an anisotropic :math:`R` below P1
   * - spatial / topological
     - can the method's operator express cross-face coupling?
     - **diffusion refuses periodic** — its codomain is a per-face
       scalar with no slot for a face pair. Nothing to do with angle.
   * - state-cone / sign
     - is the value representable in the method's state cone?
     - **SN refuses zero-flux** — :math:`\mathcal{A} = -1` needs a
       signed current, and :math:`\psi \ge 0` admits no negative
       angular inflow.

:math:`q` is a fourth thing entirely: a **vector in** :math:`\Gamma_-`,
not an operator, so it needs :math:`\Gamma_-` represented at whatever
fidelity the source demands. The diffusion arm refuses it for a
**plumbing** reason (#290 P5 does not exist yet), not a representability
one — that refusal will vanish with no theory changing.

The §16A.3 decomposition splits this map into three concrete layers:

1. **Trace structure** — the unified typed
   :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`
   (every supported mesh — 1-D Cartesian / spherical / cylindrical
   + 2-D Cartesian — post Issue #188; 2-D cylindrical
   :class:`Mesh2D` is the only mesh that still raises, deferred
   until a 2-D cylindrical SN sweep ships). It carries the signed
   :math:`\Omega_n \cdot \hat n_f` per face; inflow and outflow are
   *selectors* over the sign predicate
   :math:`\mathrm{sign}(\Omega_n \cdot \hat n_f)`.
2. **Boundary law** — :class:`BoundaryTraceLaw` ABC + seven concrete
   subclasses. The ABC declares three properties for the
   :eq:`bc-affine-form-init` factors
   (:attr:`~BoundaryTraceLaw.geometry_map`,
   :attr:`~BoundaryTraceLaw.response_kernel`,
   :attr:`~BoundaryTraceLaw.source`), and **all three are real and
   read**: campaign phase **B1** minted the typed ``G`` / ``R``
   specification objects and populated them on all seven concretes,
   and **B2** repointed five production sites off ``law.kind`` string
   comparisons onto them (the sweep schedule's reflective set, the
   ruled-corner predicate, the solver's leakage list, DSA admission
   and row selection, and the diffusion albedo). The ABC's ``None``
   default survives only so a stub law is detectable as unpopulated.
   Method-agnostic.
3. **Method realisation** —
   :class:`~orpheus.geometry.boundary.BoundaryRealizer` Protocol +
   two functional implementations
   (:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` and, since
   #290 P3 / issue #182,
   :class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`),
   each OWNED by its method-mesh: the mesh's ``realize_boundary_law``
   arm — the per-method hook of the
   :class:`~orpheus.transport.method.TransportMethod` Protocol —
   instantiates its own realizer directly.

One registry connects the tag layer to the law layer:

* **Law registry** — :attr:`BoundaryTraceLaw.registry` (a class-level
  dict maintained by
  :class:`~orpheus.numerics.registry.RegistryMixin`). Keyed by the
  ``BC.kind`` string (``"vacuum"``, ``"reflective"``, ``"white"``,
  ``"periodic"``, ``"albedo"``, ``"prescribed_inflow"``,
  ``"zero_flux"``). Concrete laws self-register at import time via
  the ``key=`` class-creation kwarg.

(The Wave-5 *realizer registry* — string-keyed ``method_name →
realizer`` with ``NotImplementedError`` stubs for MoC/MC/CP — was
dissolved at #290 P7b: no consumer ever resolved a realizer by name.
You hold the method-mesh, therefore you hold its realizer. The two
extension axes survive without it: adding a new BC type is one law
class + a ``key=`` registration (every realizer adds a dispatch
branch); adding a new transport method is one method-mesh + one
realizer + one per-method admission table
(``BOUNDARY_OPERATOR_REGISTRY``) — no central registration step.)


Concrete laws (Wave 7 vocabulary, post-#186 descriptor model)
=============================================================

Each of the seven concrete :class:`BoundaryTraceLaw` subclasses is a
**frozen dataclass descriptor** — instantiate it to declare the
law's parameters; realise it (via
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`) to
obtain a 1-arg callable :class:`LinearOperator`. None of the
concrete laws has an :meth:`apply` method; the realiser is the
sole bridge. The canonical SN-realised representation per law.

.. note::

   :math:`R` below is the **response factor alone** — the crossing
   :math:`\Gamma_+ \to \Gamma_-`, a scalar for every law in this
   family — matching :eq:`bc-affine-form-init` and the theory page's
   typing. :math:`G` is the geometry, an endomorphism of
   :math:`\Gamma_+`. The *composite* the realizer produces is
   :math:`R \circ G`, never called :math:`R`.

* :class:`VacuumInflow` (registry key ``"vacuum"``) — the rank-0 case
  :math:`R = 0`, :math:`q = 0`. :math:`G` is the **identity deck element**,
  not zero: the zero map is not a bijection, so it cannot be a geometry map
  at all, and ":math:`R = G = 0`" wrote one vanishing twice, once in the
  wrong tier (corrected at B3.0; **this entry outlived that correction until
  B3.4b**, the same two-phase lag the white entry had). SN realises to the
  narrowed **zero map** :math:`\Gamma_+ \to \Gamma_-` — a
  :class:`~orpheus.numerics.operator.ZeroOperator` carrying both space hooks,
  so the forward emits the zero of :math:`\Gamma_-` and the transpose the
  zero of :math:`\Gamma_+`. It was an ``IncomingOrdinateMaskTensor`` (a
  full-face projector whose preserved rows the consumer discarded) until
  **B3.2** narrowed the domain and left nothing for a projector to do.
* :class:`ReflectiveBoundary(axis, albedo)` (registry key
  ``"reflective"``) — :math:`G = G_{\text{refl}}`, the ordinate
  permutation ``quadrature.reflection_index(axis)``; :math:`R =
  \alpha`. SN realises the composite to
  :class:`~orpheus.numerics.operator.PermutationOperator` (α=1
  fast path) or
  ``ScaledOperator(α, PermutationOperator)`` (α ≠ 1).
* :class:`WhiteBoundary(axis, outward_sign, albedo)` (registry key
  ``"white"``) — :math:`G = I` (a white face fixes NO geometry) and
  :math:`R = ` :class:`LambertianReemission(\alpha)`, the cosine-weighted
  diffuse re-emission. **This entry said** :math:`G = G_{\text{diff}}`
  **until B3.4a** — the exact misassignment campaign phase B3.0 corrected
  and this docstring outlived by two phases. An average is not
  multiplicative and not a bijection, so it cannot be a deck
  transformation; the physics settled into the empty :math:`G` slot
  because a rank-one response annihilates :math:`G` entirely, making the
  error unobservable. SN realises the composite to
  :class:`~orpheus.sn.boundary.angular.AngularAverageOperator` (α=1 fast
  path) or scaled — since **B3.4a** typed :math:`\Gamma_+ \to \Gamma_-`.
* :class:`PeriodicBoundary(axis)` (registry key ``"periodic"``) —
  :math:`G = ` :class:`SpatialWrap(axis)`, the translation carrying the
  partner face onto this one; :math:`R = 1` (a periodic face is loss-free
  by construction, and a pure symmetry statement adds no physics). SN
  realises (campaign phase **B3.4c**) to a bare
  :class:`~orpheus.numerics.operator.IdentityOperator` fed the PARTNER
  face's :math:`\Gamma_+`: :meth:`SpatialWrap.domain_face` names the partner
  and the composition supplies it, so the pushforward lives in the CHANNEL and
  the action on the trace is the identity — earned by the opposite-normals
  identification :math:`\Gamma_+(f') \equiv \Gamma_-(f)`, which the realizer
  certifies. (Until B3.4c it realised to a ``PeriodicWrapOperator`` that never
  read the partner face; issue #183 was that gap, and the type retired with it
  — its body was the identity and its content is now on the factor.) The law
  gained its ``axis`` field at B1; the partner face
  is deliberately NOT a field, since which face is opposite depends on
  where the law is installed (configuration) while "wrap along x" is
  intrinsic.
* :class:`AlbedoBoundary(albedo, reemission)` (registry key ``"albedo"``) —
  an absorbing **surface**. :math:`G = ` :meth:`SelfPairedDeck.identity`
  **always**: a wall is not a quotient of the domain, so it fixes no geometry. All the
  content is in :math:`R = \alpha\,C`, where the *re-emission closure*
  :math:`C` names the angular shape of the return —
  :class:`SpecularReturn(axis)` :math:`\Rightarrow`
  :class:`SpecularReemission`, :class:`IsotropicReturn(axis, sign)`
  :math:`\Rightarrow` :class:`LambertianReemission`, and no closure
  :math:`\Rightarrow` :class:`ScalarResponse`. The closure is
  amplitude-FREE so :math:`\alpha` has one home, on the law.

  SN realises the two completions through the **same bodies** as the
  geometry-tier laws (the pairing, and the Lambertian average), so
  ``AlbedoBoundary(α, SpecularReturn(a)) ≡ ReflectiveBoundary(a, α)`` and
  ``AlbedoBoundary(α, IsotropicReturn(a, s)) ≡ WhiteBoundary(a, s, α)`` as
  matrices while asserting different physics. SN **REFUSES** the
  closure-free spelling: :math:`\alpha\,I` is an endomorphism of
  :math:`\Gamma_+` and :math:`G` supplies no crossing, so on an angular
  trace nothing says which outgoing direction feeds which incoming one.
  Diffusion realises that same object unchanged — on a scalar trace
  :math:`J^- = \alpha J^+` is the complete law (**B3.4b**; the
  angular-resolution axis of :ref:`bc-method-realizability`). Before B3.4b
  the SN arm answered with full-face endomorphisms (``ZeroOperator`` /
  ``IdentityOperator`` / ``α·(I & I)``) that the composite then read
  positionally.
* :class:`PrescribedInflow(source)` (registry key
  ``"prescribed_inflow"``, Wave 7 addition) — the rank-0 affine BC
  :math:`R = 0`, :math:`q \neq 0`. :math:`G` is the **identity deck
  element**, not zero: the zero map is not a bijection, so it cannot be a
  geometry map at all, and the older spelling ":math:`R = G = 0`" wrote
  one vanishing twice, once in the wrong tier (corrected at B3.0; this
  entry outlived that correction until B3.4a). SN realises to
  :class:`~orpheus.sn.boundary.angular.IncomingSourceOperator`, whose
  :meth:`apply` ignores the outgoing flux and asks the source to fill
  ``(|Γ₋|,) + psi_out.shape[1:]``. Since **B3.4a** the delivered
  :math:`q` lives on :math:`\Gamma_-` **by typing** — there are no other
  rows in the codomain to write — where it previously lived there because
  an inflow MASK erased them (#52 / ERR-047). Both the mask and its
  unmasked companion branch are retired. Consumes an
  :class:`InflowSourceSpec` (typically :class:`NoSource` for the
  homogeneous case or :class:`ConstantInflowSource(value)` for a uniform
  inflow).
* :class:`ZeroFluxBoundary` (registry key ``"zero_flux"``, #290 P3 /
  user ruling 3) — the homogeneous Dirichlet idealization
  :math:`\phi_\Gamma = 0`: albedo-family member :math:`\mathcal{A} =
  -1` in the partial-current basis (:math:`J^- = -J^+`), deliberately
  outside the sub-Markov :math:`[0, 1]` range. The diffusion realizer
  maps it to ``ScaledOperator(-1, IdentityOperator)``; the SN realizer
  REFUSES it (a negative angular inflow is unrepresentable — use
  :class:`VacuumInflow` for the physical zero-incoming law).


Rank-N boundaries via the descriptor-tree algebra
=================================================

Rank-N (Marshak, partial-current) boundary conditions are
**not** a dedicated class. They are built directly through the
**descriptor-tree algebra** that
:class:`BoundaryTraceLaw` exposes via its ``+`` / ``-`` / ``*`` /
``/`` / unary ``-`` dunders. These return
:class:`~orpheus.geometry.boundary._composition.LawSum` /
:class:`~orpheus.geometry.boundary._composition.LawScaled` nodes —
pure descriptor structures with **no** ``apply`` method. The
:func:`realize_recursively` walker (in :mod:`._realizer` since
#290 P7b; method-blind — the leaf realizer is the caller's) is
the **sole** type transformer from descriptor tree to operator
tree (Issue #186 / B3 + β2, 2026-05-11):

.. code-block:: python

    from orpheus.geometry.boundary import (
        ReflectiveBoundary, WhiteBoundary, realize_recursively,
    )
    from orpheus.sn.boundary.realizer import SNBoundaryRealizer
    from orpheus.sn.mesh.method_space import SNMethodSpace

    # Build the descriptor tree (no realisation yet).
    spec = ReflectiveBoundary(axis="x", albedo=1.0)
    white = WhiteBoundary(axis="x", outward_sign=+1, albedo=1.0)
    marshak_law = 0.3 * spec + 0.7 * white
    # marshak_law is LawSum(LawScaled(0.3, ReflectiveBoundary(...)),
    #                       LawScaled(0.7, WhiteBoundary(...))).
    # Not callable: marshak_law.apply does NOT exist.

    # Realise the tree at one face, with the method's own realizer.
    ms = SNMethodSpace.for_face(...)
    marshak_op = realize_recursively(marshak_law, ms, SNBoundaryRealizer())
    # marshak_op is:
    #   OperatorSum(ScaledOperator(0.3, PermutationOperator(...)),
    #               ScaledOperator(0.7, AngularAverageOperator(...)))
    psi_in = marshak_op.apply(psi_out)   # 1-arg LinearOperator

The output is a Wave-0
:class:`~orpheus.numerics.operator.OperatorSum` of
:class:`~orpheus.numerics.operator.ScaledOperator` wrappers around
realised Wave-0 primitives, consumable by the SN sweep / Krylov
path via the uniform 1-arg :meth:`apply`.

The descriptor-tree and operator-tree algebras are **two separate
type families**: the descriptor tree
(:class:`~orpheus.geometry.boundary.LawSum` /
:class:`~orpheus.geometry.boundary.LawScaled` over
:class:`BoundaryTraceLaw` leaves) carries no ``apply``; the
operator tree (:class:`OperatorSum` / :class:`ScaledOperator` over
:class:`LinearOperator` leaves) carries :meth:`apply`. The two
algebras never inter-compose — mixing a :class:`LawNode` with an
already-realised :class:`LinearOperator` via ``+`` is a static-type
error; the caller MUST :func:`realize_recursively` the descriptor
tree first. The static type checker enforces this distinction —
no runtime contract is needed.

Two retired predecessors converged on this design:

* **Wave 11**'s ``MixedBoundaryOperator(components: list[tuple[
  float, BoundaryOperator]])`` class delayed realisation until
  apply-time; deleted because the delayed-realisation pattern
  broke down once vacuum needed per-face inflow indices that the
  bare-law container had no access to.
* **β1** (the Issue #186 / B3 interim that kept
  :class:`LinearOperator` on :class:`BoundaryTraceLaw`):
  ``0.3 * spec + 0.7 * white`` produced an :class:`OperatorSum`
  with raw-law leaves whose realisation was deferred to
  apply-time. β1 was algebraically equivalent to β2 but conflated
  the two type families — the resulting :class:`OperatorSum`
  instance looked like an operator tree to the type checker
  while behaving as a descriptor tree at runtime.

See Bell & Glasstone 1970 §1.5 for the physics; see
:doc:`/theory/foundations/boundary_conditions` § "Descriptor-tree algebra for
rank-N boundaries" and § "The trace-law descriptor model
(Issue #186 / B3 + β2)" for the architectural retrospective.


Descriptor-tree composition module
==================================

The :mod:`_composition` submodule houses the descriptor-tree
composition dataclasses introduced in Issue #186:

* :class:`LawSum(a, b)` — frozen dataclass representing
  :math:`a + b` over two :data:`LawNode` operands. Closed under
  ``+`` / ``-`` / ``*`` / ``/`` / unary ``-``.
* :class:`LawScaled(scalar, inner)` — frozen dataclass
  representing :math:`\alpha \cdot \mathrm{inner}`. Closed under
  the same dunders; ``LawScaled.__mul__`` performs **constant
  folding** (``α * (β * x) → (α β) * x`` at construction).
* :data:`LawNode` — type alias
  ``Union[BoundaryTraceLaw, LawSum, LawScaled]``, the input type
  of :func:`realize_recursively`.

All three are re-exported below for downstream consumers. None
has an :meth:`apply` method — the descriptor tree becomes
callable only after :func:`realize_recursively`. See
:mod:`._composition` for the algebra closure details.


Universal invariants + named-error catalog
==========================================

The ABC declares five universal ``assert_*`` invariants (per Grand
Report v3 §16A.12 + §27.6). What is actually implemented, measured:

* ``assert_source_lives_on_incoming_trace`` — the only one with a
  base body (the ERR-047 gate). No law overrides it; every law is
  certified by that one body.
* ``assert_geometry_map_measure_preserving`` — empty base;
  overridden by :class:`ReflectiveBoundary` alone.
* ``assert_response_positive_if_declared`` — empty base; overridden
  by :class:`WhiteBoundary` and :class:`AlbedoBoundary`.
* ``assert_inflow_outflow_classification`` and
  ``assert_outgoing_leakage_unconstrained`` — empty base,
  **overridden by nobody**. They are permanent no-ops that assert
  nothing about any law today.

So four of the seven concrete laws (:class:`VacuumInflow`,
:class:`PeriodicBoundary`, :class:`PrescribedInflow`,
:class:`ZeroFluxBoundary`) override no universal invariant at all.
Eight typed errors in :mod:`._errors` replace the pre-refactor
generic :class:`ValueError` raises:

* :class:`IncomingOutgoingTraceClassificationError` — ERR-040, fires
  on tangential ordinates where strict partition was required.
* :class:`VacuumAppliedToOutgoingTraceError` — ERR-041, fires when
  vacuum is mistakenly applied on :math:`\Gamma_+`.
* :class:`BoundaryGeometryMapNotMeasurePreservingError` — ERR-042,
  fires when the geometric map :math:`G` does not preserve
  :math:`w(\Omega)|\Omega \cdot \hat n|`.
* :class:`BoundaryResponseNotPositiveError` — ERR-043, fires on
  negative response kernel output.
* :class:`ReflectionNotInvolutiveError` — ERR-044, fires when
  ``perm[perm] != arange``.
* :class:`ReflectionDidNotMapInflowToOutflowError` — ERR-045, fires
  when a reflection maps an inflow ordinate to itself.
* :class:`SubmarkovViolationError` — ERR-046, fires when
  :math:`\alpha > 1` on a sub-Markov BC (albedo, white).
* :class:`BoundarySourceNotOnIncomingTraceError` — ERR-047, fires on
  a :class:`InflowSourceSpec` with nonzero entries on :math:`\Gamma_+`.

All eight extend :class:`ValueError` (via the
:class:`BoundaryError` base) for backward compatibility with
existing ``except ValueError`` consumers. The
``@pytest.mark.catches("ERR-NNN")`` decorators on the relevant
fault-injection tests in
:mod:`tests.geometry.test_bc_universal_invariants` pin the error
firing under the right conditions.


Cross-method realizers: two functional, and no stubs
====================================================

The :class:`BoundaryRealizer` Protocol has exactly two
implementations —
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer` (angular:
per-ordinate masks / permutations / averages) and
:class:`~orpheus.diffusion.boundary_realizer.DiffusionBoundaryRealizer`
(#290 P3 / issue #182 — scalar: every law collapses to the
albedo-family response :math:`J^- = \mathcal{A} J^+` on the
partial-current trace). MoC, MC and CP have **none**: their
``NotImplementedError`` stub realizers existed only to populate the
Wave-5 string-keyed registry, and went with it at #290 P7b (see the
parenthetical under "One registry connects the tag layer to the law
layer" above). A method adopting the unified BC model ships its
realizer next to its method-mesh — there is no central registration
step to hold open, so an empty slot costs nothing and a stub bought
nothing.

The "unify after two instances" trigger fired at the second
functional realizer (#290 P3): the shared surface turned out to be
the :class:`BoundaryRealizer` Protocol itself plus the
:func:`stamp_boundary_role` helper — no ``BoundaryRealizerBase`` ABC
is warranted (the two realize bodies share no algorithm, only the
contract). The genuinely shared *machinery* — the descriptor-tree
walker generalised over realizers and a typed ``MethodSpace``
Protocol — is the ``TransportMethod`` mint tracked as #290 P7.


Package layout (Wave 4 source-layout split, post-#186 descriptor cleanup)
=========================================================================

* :mod:`_base` -- :class:`BoundaryTraceLaw` ABC (Wave 3 / Wave 7
  ABC merge; Issue #186 / B3 + β2 dropped
  :class:`~orpheus.numerics.operator.LinearOperator`
  inheritance and the abstract ``apply``).
* :mod:`_composition` -- :class:`LawSum` /
  :class:`LawScaled` / :data:`LawNode` descriptor-tree composition
  dataclasses (Issue #186 / B3 + β2).
* :mod:`_errors` -- :class:`BoundaryError` + 8 typed subclasses
  (Wave 3 / ERR-040..ERR-047).
* :mod:`_source` -- :class:`InflowSourceSpec` Protocol +
  :class:`NoSource` + :class:`ConstantInflowSource` (Wave 3).
* :mod:`_realizer` -- the realization seam: :class:`BoundaryRealizer`
  Protocol + :func:`stamp_boundary_role` + the
  :func:`realize_recursively` rank-N walker (Wave 5; the walker moved
  in from ``orpheus.sn.boundary`` and the Wave-5
  ``BoundaryRealizerRegistry`` dissolved at #290 P7b).
* :mod:`_bound_compat` -- ``_BoundBoundaryOperator`` shim wrapping
  the realised operator with a ``kind`` string tag (Wave 8
  introduced; Wave 9 added curvilinear bound-quadrature path,
  retired Issue #176 / C176.1; Issue #186 / C-B3.4 dropped the
  ``*_extra, **_kw`` swallow). Now a **strict 1-arg** passthrough.
  **Internal** — not in this module's :attr:`__all__`.
* :mod:`vacuum` -- :class:`VacuumInflow`.
* :mod:`reflective` -- :class:`ReflectiveBoundary`.
* :mod:`white` -- :class:`WhiteBoundary`.
* :mod:`periodic` -- :class:`PeriodicBoundary`.
* :mod:`albedo` -- :class:`AlbedoBoundary`.
* :mod:`prescribed_inflow` -- :class:`PrescribedInflow` (Wave 7).
* :mod:`zero_flux` -- :class:`ZeroFluxBoundary` (#290 P3).

The pre-Wave-7 ``mixed.py`` submodule (carrying the now-retired
``MixedBoundaryOperator``) was **deleted in Wave 11**; the registry
no longer contains a ``"mixed"`` key and a test pins that absence
(:func:`tests.geometry.test_boundary.test_registry_contains_all_primitives`).

The Wave-7 deprecated aliases (``VacuumBoundaryOperator``,
``SpecularBoundaryOperator``, ``WhiteBoundaryOperator``,
``PeriodicBoundaryOperator``, ``AlbedoBoundaryOperator``) were
**retired in Wave O step O.4a.1** once every consumer had migrated.
None is importable; the canonical names above are the sole exported
symbols. The pre-refactor → canonical naming index is kept on the
theory page (:doc:`/theory/foundations/boundary_conditions` §
"Naming audit") for readers tracing pre-Wave-O commits.


Cross-references
================

* :doc:`/theory/foundations/boundary_conditions` — the master theory page
  carrying the full §16A.3 three-layer decomposition, the affine
  form derivation, the dual-registry design, the universal
  invariants, the named-error catalog, the vacuum semantic
  correction (§16A.5), the Wave-0 rank-N algebra, the
  worked end-to-end example, and the Cartesian / curvilinear
  split.
* :doc:`/theory/methods/sn/index` § "Boundary Conditions" —
  the SN-side consumption: the shared
  :func:`~orpheus.transport.method.resolve_boundary_conditions` body
  walks the mesh axes' declarations and dispatches each parsed law
  through ``SNMesh.realize_boundary_law`` →
  :class:`SNBoundaryRealizer` to produce the resolved 1-arg
  :class:`~orpheus.numerics.operator.LinearOperator`.
* :doc:`/theory/foundations/operator_algebra` § "Boundary conditions as
  Wave-0 / Wave-1 primitives" — the operator-algebra view of the
  realisation: every realised BC IS a Wave-0
  :class:`LinearOperator` composable with the rest of the algebra.
* :mod:`orpheus.numerics.operator` — Wave-0 primitives
  (:class:`PermutationOperator`, :class:`IncomingOrdinateMaskTensor`,
  :class:`IdentityOperator`, :class:`ZeroOperator`, plus the
  :class:`OperatorSum` / :class:`ScaledOperator` composers that
  appear in the *operator* tree).
* :mod:`._composition` — :class:`LawSum` / :class:`LawScaled` /
  :data:`LawNode` composers that appear in the *descriptor* tree
  (Issue #186 / B3 + β2).
* :mod:`._base` — :class:`BoundaryTraceLaw` ABC: pure descriptor
  with no ``apply``; carries the minimal algebra dunders that
  build descriptor-tree nodes.
* :mod:`orpheus.sn.boundary.realizer` —
  :class:`SNBoundaryRealizer` (functional realizer dispatching by
  ``isinstance``; the leaf descriptor → operator transformer).
* :mod:`._realizer` —
  :func:`realize_recursively` walker, the **type transformer**
  from a descriptor tree
  (``BoundaryTraceLaw | LawSum | LawScaled``) to an operator tree
  (Wave-0 :class:`OperatorSum` / :class:`ScaledOperator` around
  realised 1-arg :class:`LinearOperator` leaves).


References
----------

* Grand Report v3 §15.2, §16, §16A.1–5, §16A.10–12, §26A.4, §27.6,
  §28. Source: ``.claude/plans/neutron_transport_grand_report_v3.md``.
* The 12-wave refactor plan:
  ``.claude/plans/transient-giggling-cake.md`` (Waves 0-12).
* The post-Wave-12 cleanup plans:

  - ``.claude/plans/curvilinear-realizer-and-2arg-cleanup.md`` —
    Issue #188 + #176 (curvilinear realiser unification + Option-A
    interim).
  - ``.claude/plans/bc-trace-law-descriptor-cleanup.md`` —
    Issue #186 / B3 + β2 (descriptor-model cleanup).
* Lewis, E. E. & Miller, W. F. (1984). *Computational Methods of
  Neutron Transport*. American Nuclear Society. §3.4 (boundary
  conditions in transport).
* Bell, G. I. & Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.5 (albedo, white, and Marshak boundary
  conditions).
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Abstract base — Wave 7 merged the legacy ``BoundaryOperator`` ABC into
# :class:`BoundaryTraceLaw`. Wave O step O.4a.1 retired the deprecated
# ``BoundaryOperator`` alias (all consumers moved to ``BoundaryTraceLaw``;
# the name is now free for the ``BlockRole.BOUNDARY`` marker in
# :mod:`orpheus.numerics.operator`).
# ---------------------------------------------------------------------------

from ._base import BoundaryTraceLaw, law_permutes_ordinates

# ---------------------------------------------------------------------------
# Issue #186 (B3 + β2) -- descriptor-tree composition. LawSum / LawScaled
# form a closed algebra over BoundaryTraceLaw | LawSum | LawScaled, used
# for rank-N boundary composition (e.g. ``0.3 * spec + 0.7 * white``).
# Realised to a Wave-0 operator tree via :func:`realize_recursively`
# (in ._realizer since #290 P7b).
# ---------------------------------------------------------------------------

from ._composition import LawNode, LawScaled, LawSum
from ._factors import (
    BoundaryGeometryMap,
    BoundaryResponseKernel,
    IsotropicReturn,
    LambertianReemission,
    ReemissionClosure,
    ScalarResponse,
    SelfPairedDeck,
    SpatialWrap,
    SpecularReemission,
    SpecularReturn,
)

# ---------------------------------------------------------------------------
# Wave 3 -- typed error catalog (ERR-040..ERR-047).
# ---------------------------------------------------------------------------

from ._errors import (
    BoundaryError,
    BoundaryGeometryMapNotMeasurePreservingError,
    BoundaryResponseNotPositiveError,
    BoundarySourceNotOnIncomingTraceError,
    IncomingOutgoingTraceClassificationError,
    ReflectionDidNotMapInflowToOutflowError,
    ReflectionNotInvolutiveError,
    SubmarkovViolationError,
    VacuumAppliedToOutgoingTraceError,
)

# ---------------------------------------------------------------------------
# Wave 3 -- prescribed-inflow source.
# ---------------------------------------------------------------------------

from ._source import InflowSourceSpec, ConstantInflowSource, NoSource

# ---------------------------------------------------------------------------
# The realization seam (Wave 5 Protocol + role stamp; the rank-N walker
# moved in from ``orpheus.sn.boundary`` and the Wave-5 registry
# dissolved at #290 P7b).
# ---------------------------------------------------------------------------

from ._realizer import (
    BoundaryRealizer,
    realize_recursively,
    stamp_boundary_role,
)

# ---------------------------------------------------------------------------
# Concrete BCs -- split into per-BC submodules in Wave 4, renamed to
# Grand Report v3 vocabulary in Wave 7. Wave O step O.4a.1 retired the
# deprecated ``*BoundaryOperator`` aliases (all consumers moved to the
# canonical names below).
# ---------------------------------------------------------------------------

from .albedo import AlbedoBoundary
from .periodic import PeriodicBoundary
from .prescribed_inflow import PrescribedInflow
from .reflective import ReflectiveBoundary
from .vacuum import VacuumInflow
from .white import WhiteBoundary
from .zero_flux import ZeroFluxBoundary


__all__ = [
    # Abstract base
    "BoundaryTraceLaw",
    # The composite question BOTH tiers can answer (B3.4b): does the
    # realized R∘G permute the angular index? Four production sites ask
    # it; it is a function, not a Protocol member, so the two factor
    # tiers stay structurally disjoint.
    "law_permutes_ordinates",
    # The affine form's two operator factors, as typed SPECIFICATIONS
    # (Grand Report v3 §16A.2; campaign phases B1 + B3). ``G : Γ₊ → Γ₋`` is the
    # DECK TRANSFORMATION — it carries the crossing, because the mirror that
    # exchanges the hemispheres is geometry — and ``R : Γ₋ → Γ₋`` is the
    # CONSTITUTIVE response; the realized composite is ``B``. Membership takes
    # TWO tests (``_factors`` docstring): multiplicativity is NECESSARY — a
    # relabeling satisfies ``G(ψφ) = (Gψ)(Gφ)``, an average never does, which
    # is what moved the Lambertian across the line in B3
    # (``HemisphericalAverage`` -> ``LambertianReemission``) and retired
    # ``NullMap``, whose ``G = 0`` is not a bijection and merely respelled
    # ``ScalarResponse(0.0)`` a tier too high. It is NOT sufficient: a specular
    # kernel is a permutation, hence multiplicative, and still constitutive.
    # The sufficient test is the quotient one — ``G`` is the deck
    # transformation of an ACTUAL quotient of the domain — whence EXACTLY ONE
    # of ``G``, ``R`` is non-trivial (B3.4b, user ruling 2026-08-01).
    "BoundaryGeometryMap",
    "BoundaryResponseKernel",
    "LambertianReemission",
    "ScalarResponse",
    # The self-paired half of the deck tier (G5): one type carrying the
    # involutive rigid motion of a face paired with ITSELF. Its two
    # inhabitants are the trivial pairing and the mirror.
    "SelfPairedDeck",
    "SpatialWrap",
    "SpecularReemission",
    # The re-emission CLOSURE tier: amplitude-free angular shapes a surface
    # law instantiates at its own α. ``AlbedoBoundary`` is the one consumer —
    # its ``reemission`` field — and the shapes exist because α and the
    # angular distribution are independent degrees of freedom (B3.4b).
    "ReemissionClosure",
    "IsotropicReturn",
    "SpecularReturn",
    # Descriptor-tree composition (Issue #186)
    "LawNode",
    "LawScaled",
    "LawSum",
    # Errors
    "BoundaryError",
    "BoundaryGeometryMapNotMeasurePreservingError",
    "BoundaryResponseNotPositiveError",
    "BoundarySourceNotOnIncomingTraceError",
    "IncomingOutgoingTraceClassificationError",
    "ReflectionDidNotMapInflowToOutflowError",
    "ReflectionNotInvolutiveError",
    "SubmarkovViolationError",
    "VacuumAppliedToOutgoingTraceError",
    # Source
    "InflowSourceSpec",
    "ConstantInflowSource",
    "NoSource",
    # The realization seam: Protocol + rank-N walker + role stamp
    "BoundaryRealizer",
    "realize_recursively",
    "stamp_boundary_role",
    # Concrete BCs (Wave 7 canonical names + the #290 P3 addition)
    "AlbedoBoundary",
    "PeriodicBoundary",
    "PrescribedInflow",
    "ReflectiveBoundary",
    "VacuumInflow",
    "WhiteBoundary",
    "ZeroFluxBoundary",
]
